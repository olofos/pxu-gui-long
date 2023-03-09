//XXX::new(xp(C::from(p0), 1.0, consts)).goto_im(xm(C::from(p0), 1.0, consts)).goto_p(p_start+1.5*PI).goto_m(0.0).generate_contour()
use crate::{
    kinematics::{dxm_dp, dxp_dp, xm, xp, CouplingConstants},
    nr,
};

use std::collections::VecDeque;
use std::f64::consts::TAU;

use itertools::Itertools;
use num::complex::Complex64;

#[derive(Debug, Copy, Clone)]
pub enum InterpolationPoint {
    Xp(f64, f64),
    Xm(f64, f64),
    C(Complex64),
}

impl InterpolationPoint {
    pub fn evaluate(&self, consts: CouplingConstants) -> Complex64 {
        match self {
            Self::Xp(p, m) => xp(Complex64::from(p), *m, consts),
            Self::Xm(p, m) => xm(Complex64::from(p), *m, consts),
            Self::C(c) => *c,
        }
    }
}

enum InterpolationStrategy {
    XpConstM(f64, f64, f64),
    XpConstP(f64, f64, f64),
    XmConstM(f64, f64, f64),
    XmConstP(f64, f64, f64),
    Lerp(InterpolationPoint, InterpolationPoint),
}

impl InterpolationStrategy {
    fn new(pt1: InterpolationPoint, pt2: InterpolationPoint) -> Self {
        use InterpolationPoint::{Xm, Xp};
        use InterpolationStrategy::*;

        let pts = (pt1, pt2);

        if let (Xp(p1, m1), Xp(p2, m2)) = pts {
            if p1 == p2 {
                return XpConstP(p1, m1, m2);
            } else if m1 == m2 {
                return XpConstM(m1, p1, p2);
            }
        } else if let (Xm(p1, m1), Xm(p2, m2)) = pts {
            if p1 == p2 {
                return XmConstP(p1, m1, m2);
            } else if m1 == m2 {
                return XmConstM(m1, p1, p2);
            }
        }
        Lerp(pts.0, pts.1)
    }

    fn max_step(&self) -> f64 {
        match self {
            Self::XpConstP(_, _, _) | Self::XmConstP(_, _, _) => 1.0 / 8.0,
            Self::XpConstM(_, p1, p2) | Self::XmConstM(_, p1, p2) => {
                0.5f64.min(1.0 / 128.0 / (p2 - p1).abs())
            }
            Self::Lerp(_, _) => 1.0 / 16.0,
        }
    }

    fn evaluate(&self, t: f64, consts: CouplingConstants) -> InterpolationPoint {
        match self {
            Self::XpConstP(p, m1, m2) => {
                let m = m1 + t * (m2 - m1);
                InterpolationPoint::Xp(*p, m)
            }
            Self::XmConstP(p, m1, m2) => {
                let m = m1 + t * (m2 - m1);
                InterpolationPoint::Xm(*p, m)
            }
            Self::XpConstM(m, p1, p2) => {
                let p = p1 + t * (p2 - p1);
                InterpolationPoint::Xp(p, *m)
            }
            Self::XmConstM(m, p1, p2) => {
                let p = p1 + t * (p2 - p1);
                let pt = InterpolationPoint::Xm(p, *m);
                pt
            }
            Self::Lerp(c1, c2) => {
                let z1 = c1.evaluate(consts);
                let z2 = c2.evaluate(consts);
                InterpolationPoint::C((1.0 - t) * z1 + t * z2)
            }
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum InterpolationComponent {
    Xp,
    Xm,
}

pub struct XInterpolator {}

#[derive(Debug)]
enum Asymptote {
    Num(Complex64),
    Infinity,
}

impl Asymptote {
    fn is_infinite(&self) -> bool {
        match self {
            Self::Infinity => true,
            _ => false,
        }
    }

    fn num(n: impl Into<Complex64>) -> Self {
        Self::Num(n.into())
    }
}

impl XInterpolator {
    pub fn generate_xp_full(p_range: i32, m: f64, consts: CouplingConstants) -> Vec<Complex64> {
        Self::generate_xp(p_range as f64, (p_range + 1) as f64, m, consts)
    }

    pub fn generate_xm_full(p_range: i32, m: f64, consts: CouplingConstants) -> Vec<Complex64> {
        Self::generate_xm(p_range as f64, (p_range + 1) as f64, m, consts)
    }

    pub fn generate_xp(
        p_start: f64,
        p_end: f64,
        m: f64,
        consts: CouplingConstants,
    ) -> Vec<Complex64> {
        if p_start > p_end {
            return Self::generate_xp(p_end, p_start, m, consts);
        }
        if p_start.floor() != p_end.ceil() - 1.0 {
            log::error!("Trying to generate xp for more than one momentum region ({p_end:2} > {p_start:2}+1)");
        }

        let x_start: Asymptote;
        let x_end: Asymptote;

        if p_start.floor() == p_start {
            let p_range = p_start as i32;
            // if (p_range == 0 && m == 0.0) || (p_range == -1 && m == consts.k() as f64) {
            if (p_start * consts.k() as f64) + m == 0.0 {
                x_start = Asymptote::num(consts.s());
            } else if p_range < 0 {
                x_start = Asymptote::num(0.0);
            } else {
                x_start = Asymptote::Infinity;
            }
        } else {
            x_start = Asymptote::num(xp(p_start, m, consts));
        }

        if p_end.ceil() == p_end {
            let p_range = p_end as i32 - 1;
            if (p_range == -1 && m == 0.0) || (p_range == -2 && m == consts.k() as f64) {
                x_end = Asymptote::num(-1.0 / consts.s());
            } else if p_range < -1 {
                x_end = Asymptote::num(0.0);
            } else {
                x_end = Asymptote::Infinity;
            }
        } else {
            x_end = Asymptote::num(xp(p_end, m, consts));
        }

        const INFINITY: f64 = 256.0;
        let mut points = VecDeque::new();

        if let Asymptote::Num(n) = x_start {
            points.push_back((p_start + 0.0, n));
        }

        let step_size = 1.0 / 8.0;
        if p_end - p_start > 2.0 * step_size {
            let i1 = ((p_start + step_size) / step_size).floor() as i32;
            let i2 = ((p_end - step_size) / step_size).floor() as i32;

            for i in i1..=i2 {
                let p = i as f64 * step_size;
                points.push_back((p, xp(p, m, consts)));
            }
        } else {
            let p = (p_start + p_end) / 2.0;
            points.push_back((p, xp(p, m, consts)));
        }

        if let Asymptote::Num(n) = x_end {
            points.push_back((p_end, n));
        }

        if x_start.is_infinite() {
            loop {
                let (p, x) = points.front().unwrap();
                let dp = p - p_start;
                if x.norm_sqr() > INFINITY * INFINITY || dp < 1.0 / 16384.0 {
                    break;
                }
                let p = p_start + dp / 4.0;
                points.push_front((p, xp(p, m, consts)));
            }
        }

        if x_end.is_infinite() {
            loop {
                let (p, x) = points.back().unwrap();
                let dp = p_end - p;
                if x.norm_sqr() > INFINITY * INFINITY || dp < 1.0 / 16384.0 {
                    break;
                }
                let p = p_end - dp / 4.0;
                points.push_back((p, xp(p, m, consts)));
            }
        }

        let mut points = Vec::from(points);

        let min_cos = (2.0 * TAU / 360.0).cos();

        for i in 0.. {
            let mut refinements = vec![];

            let mut prev = false;
            for ((p1, x1), (p2, x2), (p3, x3)) in points.iter().tuple_windows::<(_, _, _)>() {
                let z1 = x2 - x1;
                let z2 = x3 - x2;
                let cos = (z1.re * z2.re + z1.im * z2.im) / (z1.norm() * z2.norm());
                if cos < min_cos {
                    if !prev {
                        refinements.push((p1 + p2) / 2.0);
                    }
                    refinements.push((p2 + p3) / 2.0);
                } else {
                    prev = false;
                }
            }

            if refinements.is_empty() {
                break;
            }
            points.extend(refinements.into_iter().map(|p| (p, xp(p, m, consts))));
            points.sort_by(|(p1, _), (p2, _)| p1.partial_cmp(p2).unwrap());

            if i >= 5 {
                break;
            }
        }
        points.into_iter().map(|(_, x)| x).collect()
    }

    pub fn generate_xm(
        p_start: f64,
        p_end: f64,
        m: f64,
        consts: CouplingConstants,
    ) -> Vec<Complex64> {
        Self::generate_xp(p_start, p_end, m, consts)
            .into_iter()
            .map(|x| x.conj())
            .collect()
    }
}

#[derive(Debug, Clone)]
pub struct PInterpolator {
    valid: bool,
    pub component: InterpolationComponent,
    pub p: Complex64,
    pub pt: InterpolationPoint,
    pub consts: CouplingConstants,

    pub p_path: Vec<Complex64>,
    pub x_path: Vec<Complex64>,
}

/// let contour = MomentumInterpolator::xp(0.25, consts).goto(InterpolationPoint::Xm(0.25, 1.0)).goto(InterpolationPoint::Xm(0.25, 0.0))).full_contour();

impl PInterpolator {
    pub fn xp(p: f64, consts: CouplingConstants) -> Self {
        let pt = InterpolationPoint::Xp(p, 1.0);
        let p = Complex64::from(p);
        Self {
            valid: true,
            component: InterpolationComponent::Xp,
            p,
            pt,
            consts,
            p_path: vec![p],
            x_path: vec![pt.evaluate(consts)],
        }
    }

    pub fn xm(p: f64, consts: CouplingConstants) -> Self {
        let pt = InterpolationPoint::Xm(p, 1.0);
        let p = Complex64::from(p);
        Self {
            valid: true,

            component: InterpolationComponent::Xm,
            p,
            pt,
            consts,
            p_path: vec![p],
            x_path: vec![pt.evaluate(consts)],
        }
    }

    pub fn clear_path(self) -> Self {
        let p_path = vec![*self.p_path.last().unwrap()];
        let x_path = vec![*self.x_path.last().unwrap()];
        Self {
            p_path,
            x_path,
            ..self
        }
    }

    pub fn goto_re(self, re: f64) -> Self {
        let im = self.pt.evaluate(self.consts).im;
        let x = Complex64::new(re, im);
        let pt = InterpolationPoint::C(x);
        self.goto(pt)
    }

    pub fn goto_im(self, im: f64) -> Self {
        let re = self.pt.evaluate(self.consts).re;
        let x = Complex64::new(re, im);
        let pt = InterpolationPoint::C(x);
        self.goto(pt)
    }

    pub fn goto_xp(self, p: f64, m: f64) -> Self {
        let pt = InterpolationPoint::Xp(p, m);
        let mut result = self.goto(pt);
        result.pt = pt;
        result
    }

    pub fn goto_xm(self, p: f64, m: f64) -> Self {
        let pt = InterpolationPoint::Xm(p, m);
        let mut result = self.goto(pt);
        result.pt = pt;
        result
    }

    pub fn go_towards_xp(self, p: f64, m: f64) -> Self {
        self.go_towards(InterpolationPoint::Xp(p, m))
    }

    pub fn go_towards_xm(self, p: f64, m: f64) -> Self {
        self.go_towards(InterpolationPoint::Xm(p, m))
    }

    pub fn go_towards(self, pt: InterpolationPoint) -> Self {
        if !self.valid {
            return self;
        }
        let (_, int_pt) = self.do_goto(pt);
        int_pt
    }

    pub fn goto(self, pt: InterpolationPoint) -> Self {
        if !self.valid {
            return self;
        }
        let (valid, int_pt) = self.do_goto(pt);
        Self { valid, ..int_pt }
    }

    fn do_goto(mut self, pt: InterpolationPoint) -> (bool, Self) {
        if !self.valid {
            return (false, self);
        }

        let strategy = InterpolationStrategy::new(self.pt, pt);
        let mut pt = self.pt;
        let mut p = self.p;

        let mut p_path = std::mem::take(&mut self.p_path);
        let mut x_path = std::mem::take(&mut self.x_path);

        let mut t = 0.0;

        'outer: while t < 1.0 {
            let mut step = strategy.max_step().min(1.0 - t);

            for i in 0.. {
                pt = strategy.evaluate(t + step, self.consts);
                let w = pt.evaluate(self.consts);
                let next_p = nr::find_root(|z| self.f(z) - w, |z| self.df(z), p, 1.0e-3, 50);
                if let Some(next_p) = next_p {
                    if (next_p.re - p.re).abs() < 1.0 / 8.0 && (next_p.im - p.im).abs() < 1.0 / 4.0
                    {
                        t += step;
                        p = next_p;
                        p_path.push(p);
                        x_path.push(w);
                        break;
                    }
                }
                if i > 5 {
                    break 'outer;
                }
                step /= 2.0;
            }
        }

        (
            t == 1.0,
            Self {
                pt,
                p,
                p_path,
                x_path,
                ..self
            },
        )
    }

    pub fn contour(&self) -> Vec<Complex64> {
        if !self.valid {
            return vec![];
        }

        if let InterpolationPoint::C(_) = self.pt {
            log::info!("Can only generate contour from Xp or Xm");
            return vec![];
        }

        let mut result = vec![];

        result
    }

    fn f(&self, z: Complex64) -> Complex64 {
        match self.component {
            InterpolationComponent::Xp => xp(z, 1.0, self.consts),
            InterpolationComponent::Xm => xm(z, 1.0, self.consts),
        }
    }

    fn df(&self, z: Complex64) -> Complex64 {
        match self.component {
            InterpolationComponent::Xp => dxp_dp(z, 1.0, self.consts),
            InterpolationComponent::Xm => dxm_dp(z, 1.0, self.consts),
        }
    }
}
