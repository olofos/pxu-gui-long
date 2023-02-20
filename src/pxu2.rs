//XXX::new(xp(C::from(p0), 1.0, consts)).goto_im(xm(C::from(p0), 1.0, consts)).goto_p(p_start+1.5*PI).goto_m(0.0).generate_contour()
use crate::{
    kinematics::{dxm_dp, dxp_dp, xm, xp, CouplingConstants},
    nr,
};

use std::{collections::VecDeque, f64::consts::TAU};

use itertools::Itertools;
use num::{complex::Complex64, Zero};

#[derive(Debug, Copy, Clone)]
pub enum InterpolationPoint {
    Xp(f64, f64),
    Xm(f64, f64),
    C(Complex64),
}

impl InterpolationPoint {
    fn evaluate(&self, consts: CouplingConstants) -> Complex64 {
        match self {
            Self::Xp(p, m) => xp(Complex64::from(p) * TAU, *m, consts),
            Self::Xm(p, m) => xm(Complex64::from(p) * TAU, *m, consts),
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
        use InterpolationPoint::{Xm, Xp, C};
        use InterpolationStrategy::*;

        let pts = (pt1, pt2);

        if let (Xp(p1, m1), Xp(p2, m2)) = pts {
            if p1 == p2 {
                log::debug!("const p xp {p1} {m1} {m2}");
                return XpConstP(p1, m1, m2);
            } else if m1 == m2 {
                log::debug!("const m xp {m1} {p1} {p2}");
                return XpConstM(m1, p1, p2);
            }
        } else if let (Xm(p1, m1), Xm(p2, m2)) = pts {
            if p1 == p2 {
                log::debug!("const p xm {p1} {m1} {m2}");
                return XmConstP(p1, m1, m2);
            } else if m1 == m2 {
                log::debug!("const m xm {m1} {p1} {p2}");
                return XmConstM(m1, p1, p2);
            }
        }
        log::debug!("lerp {:?} {:?}", pts.0, pts.1);
        Lerp(pts.0, pts.1)
    }

    fn max_step(&self) -> f64 {
        match self {
            Self::XpConstP(_, _, _) | Self::XmConstP(_, _, _) => 1.0 / 32.0,
            Self::XpConstM(_, _, _) | Self::XmConstM(_, _, _) => 1.0 / 32.0,
            Self::Lerp(_, _) => 1.0 / 32.0,
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
enum XAsymptote {
    Num(Complex64),
    Infinity,
}

impl XAsymptote {
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
    pub fn generate_xp(p_range: i32, m: f64, consts: CouplingConstants) -> Vec<Complex64> {
        Self::generate_x(p_range, m, consts, |p| {
            xp(Complex64::from(TAU * p), m, consts)
        })
    }

    fn generate_x(
        p_range: i32,
        m: f64,
        consts: CouplingConstants,
        f: impl Fn(f64) -> Complex64,
    ) -> Vec<Complex64> {
        let p_start = p_range as f64;
        let mut points = VecDeque::new();
        const INFINITY: f64 = 256.0;

        let start: XAsymptote;
        let end: XAsymptote;

        if (p_range == 0 && m == 0.0) || (p_range == -1 && m == consts.k() as f64) {
            start = XAsymptote::num(consts.s());
        } else if p_range < 0 {
            start = XAsymptote::num(0.0);
        } else {
            start = XAsymptote::Infinity;
        }

        if (p_range == -1 && m == 0.0) || (p_range == -2 && m == consts.k() as f64) {
            end = XAsymptote::num(-1.0 / consts.s());
        } else if p_range < -1 {
            end = XAsymptote::num(0.0);
        } else {
            end = XAsymptote::Infinity;
        }

        if let XAsymptote::Num(n) = start {
            points.push_back((p_start + 0.0, n));
        }

        let steps = 16;
        for i in 1..steps {
            let p = p_start + i as f64 / steps as f64;
            points.push_back((p, f(p)));
        }

        if let XAsymptote::Num(n) = end {
            points.push_back((p_start + 1.0, n));
        }

        if start.is_infinite() {
            loop {
                let (p, x) = points.front().unwrap();
                if x.norm_sqr() > INFINITY * INFINITY || (p - p_start) < 1.0 / 16384.0 {
                    break;
                }
                let p = p_start + (p - p_start) / 4.0;
                points.push_front((p, f(p)));
            }
        }

        if end.is_infinite() {
            loop {
                let (p, x) = points.back().unwrap();
                let dp = 1.0 - (p - p_start);
                if x.norm_sqr() > INFINITY * INFINITY || dp < 1.0 / 16384.0 {
                    break;
                }
                let p = p_start + 1.0 - dp / 4.0;
                points.push_back((p, f(p)));
            }
        }

        let mut points = Vec::from(points);

        let min_cos = (2.0 * TAU / 360.0).cos();

        for i in 0.. {
            points.sort_by(|(p1, _), (p2, _)| p1.partial_cmp(p2).unwrap());

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

            if refinements.is_empty() || i >= 6 {
                break;
            }
            points.extend(refinements.into_iter().map(|p| (p, f(p))));
        }

        points.into_iter().map(|(_, x)| x).collect()
    }
}

#[derive(Debug, Clone)]
pub struct PInterpolator {
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
            component: InterpolationComponent::Xp,
            p,
            pt,
            consts,
            p_path: vec![TAU * p],
            x_path: vec![pt.evaluate(consts)],
        }
    }

    pub fn xm(p: f64, consts: CouplingConstants) -> Self {
        let pt = InterpolationPoint::Xm(p, 1.0);
        let p = Complex64::from(p);
        Self {
            component: InterpolationComponent::Xm,
            p,
            pt,
            consts,
            p_path: vec![TAU * p],
            x_path: vec![pt.evaluate(consts)],
        }
    }

    pub fn get_x(&self) -> Complex64 {
        self.pt.evaluate(self.consts)
    }

    pub fn get_p(&self) -> Complex64 {
        self.p
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
        log::debug!("{:?} {:?}", pt, result.pt);
        result.pt = pt;
        result
    }

    pub fn goto(mut self, pt: InterpolationPoint) -> Self {
        let strategy = InterpolationStrategy::new(self.pt, pt);
        let mut pt = self.pt;
        let mut p = self.p;

        let mut p_path = std::mem::take(&mut self.p_path);
        let mut x_path = std::mem::take(&mut self.x_path);

        let mut t = 0.0;
        while t < 1.0 {
            let mut step = strategy.max_step().min(1.0 - t);

            loop {
                pt = strategy.evaluate(t + step, self.consts);
                let w = pt.evaluate(self.consts);
                let next_p = nr::find_root(|z| self.f(z) - w, |z| self.df(z), p, 1.0e-3, 50);
                if let Some(next_p) = next_p {
                    t += step;
                    p = next_p;
                    p_path.push(p);
                    x_path.push(w);
                    break;
                }
                step /= 2.0;
            }
        }

        log::debug!("{pt:?}");

        Self {
            pt,
            component: self.component,
            p,
            consts: self.consts,
            p_path,
            x_path,
        }
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

    fn next(&self, pt: InterpolationPoint) -> Option<Complex64> {
        let n = nr::find_root(
            |z| self.f(z) - pt.evaluate(self.consts),
            |z| self.df(z),
            self.p,
            1.0e-3,
            50,
        );

        unimplemented!()
    }
}
