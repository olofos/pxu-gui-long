use crate::{
    kinematics::{den2_dp, dxp_dp, en2, xm, xp, CouplingConstants},
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
    Re(f64),
}

impl InterpolationPoint {
    pub fn evaluate(&self, consts: CouplingConstants) -> Complex64 {
        match self {
            Self::Xp(p, m) => xp(Complex64::from(p), *m, consts),
            Self::Xm(p, m) => xm(Complex64::from(p), *m, consts),
            Self::C(c) => *c,
            Self::Re(re) => Complex64::from(re),
        }
    }
}

enum InterpolationStrategy {
    XpConstM(f64, f64, f64),
    XpConstP(f64, f64, f64),
    XmConstM(f64, f64, f64),
    XmConstP(f64, f64, f64),
    ReLerp(f64, f64),
    Lerp(InterpolationPoint, InterpolationPoint),
}

impl InterpolationStrategy {
    fn new(pt1: InterpolationPoint, pt2: InterpolationPoint) -> Self {
        use InterpolationPoint::{Re, Xm, Xp};
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
        } else if let (Re(r1), Re(r2)) = pts {
            return ReLerp(r1, r2);
        }
        Lerp(pts.0, pts.1)
    }

    fn max_step(&self) -> f64 {
        match self {
            Self::XpConstP(_, _, _) | Self::XmConstP(_, _, _) => 1.0 / 8.0,
            Self::XpConstM(_, p1, p2) | Self::XmConstM(_, p1, p2) => {
                0.5f64.min(1.0 / 8.0 / (p2 - p1).abs())
            }
            Self::ReLerp(_, _) => 1.0 / 16.0,
            Self::Lerp(_, _) => 1.0 / 16.0,
        }
    }

    fn argument(&self, t: f64) -> f64 {
        match self {
            Self::XpConstP(_, m1, m2) => m1 + t * (m2 - m1),
            Self::XmConstP(_, m1, m2) => m1 + t * (m2 - m1),
            Self::XpConstM(_, p1, p2) => p1 + t * (p2 - p1),
            Self::XmConstM(_, p1, p2) => p1 + t * (p2 - p1),
            Self::ReLerp(r1, r2) => r1 + t * (r2 - r1),
            Self::Lerp(_, _) => t,
        }
    }

    fn evaluate_from_argument(&self, arg: f64, consts: CouplingConstants) -> InterpolationPoint {
        match self {
            Self::XpConstP(p, _, _) => {
                let m = arg;
                InterpolationPoint::Xp(*p, m)
            }
            Self::XmConstP(p, _, _) => {
                let m = arg;
                InterpolationPoint::Xm(*p, m)
            }
            Self::XpConstM(m, _, _) => {
                let p = arg;
                InterpolationPoint::Xp(p, *m)
            }
            Self::XmConstM(m, _, _) => {
                let p = arg;
                InterpolationPoint::Xm(p, *m)
            }
            Self::ReLerp(_, _) => InterpolationPoint::Re(arg),
            Self::Lerp(c1, c2) => {
                let z1 = c1.evaluate(consts);
                let z2 = c2.evaluate(consts);
                let x = (1.0 - arg) * z1 + arg * z2;
                if x.im.abs() < 1.0e-5 {
                    InterpolationPoint::Re(x.re)
                } else {
                    InterpolationPoint::C(x)
                }
            }
        }
    }

    fn evaluate(&self, t: f64, consts: CouplingConstants) -> InterpolationPoint {
        self.evaluate_from_argument(self.argument(t), consts)
    }
}

pub struct XInterpolator {}

#[derive(Debug)]
enum Asymptote {
    Num(Complex64),
    Infinity,
}

impl Asymptote {
    fn is_infinite(&self) -> bool {
        matches!(self, Self::Infinity)
    }

    fn num(n: impl Into<Complex64>) -> Self {
        Self::Num(n.into())
    }

    fn new_xp_start(p: f64, m: f64, consts: CouplingConstants) -> Self {
        if p == p.floor() {
            let m_eff = m + p * consts.k() as f64;
            if m_eff == 0.0 {
                Self::num(consts.s())
            } else if m_eff > 0.0 {
                Self::Infinity
            } else {
                Self::num(0.0)
            }
        } else {
            Self::num(xp(p, m, consts))
        }
    }

    fn new_xp_end(p: f64, m: f64, consts: CouplingConstants) -> Self {
        if p == p.ceil() {
            let m_eff = m + p * consts.k() as f64;
            if m_eff == 0.0 {
                Self::num(-1.0 / consts.s())
            } else if m_eff > 0.0 {
                Self::Infinity
            } else {
                Self::num(0.0)
            }
        } else {
            Self::num(xp(p, m, consts))
        }
    }
}

trait Refiner<T>
where
    T: Clone,
{
    fn mid_point(t1: &T, t2: &T) -> T;
    fn cmp(t1: &T, t2: &T) -> Option<std::cmp::Ordering>;

    fn refine_raw(
        points: impl Into<Vec<(T, Complex64)>>,
        eval: impl Fn(T, Complex64) -> Option<Complex64>,
    ) -> Vec<(T, Complex64)> {
        let mut points: Vec<(T, Complex64)> = points.into();

        let min_cos = (2.0 * TAU / 360.0).cos();

        for i in 0.. {
            let mut refinements: Vec<(T, Complex64)> = vec![];

            let mut prev = false;
            for ((p1, x1), (p2, x2), (p3, x3)) in points.iter().tuple_windows::<(_, _, _)>() {
                let z1 = x2 - x1;
                let z2 = x3 - x2;
                let cos = (z1.re * z2.re + z1.im * z2.im) / (z1.norm() * z2.norm());
                if cos < min_cos {
                    if !prev {
                        refinements.push((Self::mid_point(p1, p2), (x1 + x2) / 2.0));
                    }
                    refinements.push((Self::mid_point(p2, p3), (x2 + x3) / 2.0));
                } else {
                    prev = false;
                }
            }

            if refinements.is_empty() {
                break;
            }

            for (point, z) in refinements {
                let Some(x) = eval(point.clone(), z) else {continue;};
                points.push((point, x));
            }

            points.sort_by(|(p1, _), (p2, _)| Self::cmp(p1, p2).unwrap());

            if i >= 5 {
                break;
            }
        }

        points
    }

    fn refine(
        points: impl Into<Vec<(T, Complex64)>>,
        eval: impl Fn(T, Complex64) -> Option<Complex64>,
    ) -> Vec<Complex64> {
        Self::refine_raw(points, eval)
            .into_iter()
            .map(|(_, x)| x)
            .collect()
    }
}

impl XInterpolator {
    pub fn generate_xp_full(p_range: i32, m: f64, consts: CouplingConstants) -> Vec<Complex64> {
        Self::generate_xp(p_range as f64, (p_range + 1) as f64, m, consts)
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

        let x_start = Asymptote::new_xp_start(p_start, m, consts);
        let x_end = Asymptote::new_xp_end(p_end, m, consts);

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

        Self::refine(points, |p, _| Some(xp(p, m, consts)))
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

impl Refiner<f64> for XInterpolator {
    fn mid_point(t1: &f64, t2: &f64) -> f64 {
        (t1 + t2) / 2.0
    }

    fn cmp(t1: &f64, t2: &f64) -> Option<std::cmp::Ordering> {
        t1.partial_cmp(t2)
    }
}

const MAX_RE_P_JUMP: f64 = 1.0 / 8.0;
const MAX_IM_P_JUMP: f64 = 1.0 / 4.0;

#[derive(Debug, Clone)]
pub struct PInterpolatorMut {
    valid: bool,
    p: Complex64,
    pt: InterpolationPoint,
    consts: CouplingConstants,
}

impl PInterpolatorMut {
    pub fn xp(p: f64, consts: CouplingConstants) -> Self {
        let pt = InterpolationPoint::Xp(p, 1.0);
        let p = Complex64::from(p);
        Self {
            valid: true,
            p,
            pt,
            consts,
        }
    }

    pub fn goto_re(&mut self, re: f64) {
        let im = self.pt.evaluate(self.consts).im;

        let pt = if im.abs() < 1.0e-5 {
            InterpolationPoint::Re(re)
        } else {
            let x = Complex64::new(re, im);
            InterpolationPoint::C(x)
        };

        self.goto(pt);
    }

    pub fn goto_im(&mut self, im: f64) {
        let re = self.pt.evaluate(self.consts).re;
        let pt = if im.abs() < 1.0e-5 {
            InterpolationPoint::Re(re)
        } else {
            let x = Complex64::new(re, im);
            InterpolationPoint::C(x)
        };
        self.goto(pt);
    }

    pub fn goto_xp(&mut self, p: f64, m: f64) {
        let pt = InterpolationPoint::Xp(p, m);
        if self.goto(pt) {
            self.pt = pt;
        }
    }

    pub fn goto_xm(&mut self, p: f64, m: f64) {
        let pt = InterpolationPoint::Xm(p, m);
        if self.goto(pt) {
            self.pt = pt;
        }
    }

    pub fn goto_p(&mut self, p: f64) {
        if !self.valid {
            return;
        }

        match self.pt {
            InterpolationPoint::Xp(_, m) => self.goto_xp(p, m),
            InterpolationPoint::Xm(_, m) => self.goto_xm(p, m),
            _ => {
                self.valid = false;
            }
        }
    }

    pub fn goto_m(&mut self, m: f64) {
        if !self.valid {
            return;
        }

        match self.pt {
            InterpolationPoint::Xp(p, _) => self.goto_xp(p, m),
            InterpolationPoint::Xm(p, _) => self.goto_xm(p, m),
            _ => {
                self.valid = false;
            }
        }
    }

    fn find_point(&self, w: Complex64, guess: Complex64) -> Option<Complex64> {
        nr::find_root(|z| self.f(z) - w, |z| self.df(z), guess, 1.0e-5, 50)
    }

    fn goto(&mut self, pt: InterpolationPoint) -> bool {
        if !self.valid {
            return false;
        }

        let strategy = InterpolationStrategy::new(self.pt, pt);
        let mut pt = self.pt;
        let mut p = self.p;

        let mut t = 0.0;

        'outer: while t < 1.0 {
            let mut step = strategy.max_step().min(1.0 - t);

            for i in 0.. {
                pt = strategy.evaluate(t + step, self.consts);
                let w = pt.evaluate(self.consts);
                if let Some(next_p) = self.find_point(w, p) {
                    if (next_p.re - p.re).abs() < MAX_RE_P_JUMP
                        && (next_p.im - p.im).abs() < MAX_IM_P_JUMP
                    {
                        t += step;
                        p = next_p;
                        break;
                    }
                }
                if i > 5 {
                    break 'outer;
                }
                step /= 2.0;
            }
        }

        self.p = p;
        self.pt = pt;
        t == 1.0
    }

    fn generate_path(&self, pt: InterpolationPoint) -> Vec<(f64, Complex64)> {
        let strategy = InterpolationStrategy::new(self.pt, pt);
        let mut p = self.p;

        let mut path: Vec<(f64, Complex64)> = vec![(strategy.argument(0.0), p)];

        let mut t = 0.0;

        'outer: while t < 1.0 {
            let mut step = strategy.max_step().min(1.0 - t);

            for i in 0.. {
                let pt = strategy.evaluate(t + step, self.consts);
                let w = pt.evaluate(self.consts);
                if let Some(next_p) = self.find_point(w, p) {
                    if (next_p.re - p.re).abs() < MAX_RE_P_JUMP
                        && (next_p.im - p.im).abs() < MAX_IM_P_JUMP
                    {
                        t += step;
                        p = next_p;
                        path.push((strategy.argument(t), p));
                        break;
                    }
                }
                if i > 8 {
                    break 'outer;
                }
                step /= 2.0;
            }
        }

        path
    }

    pub fn contour_re(&self, x: f64) -> Vec<Complex64> {
        let (pt1, pt2) = if x > 0.0 {
            (
                InterpolationPoint::Re(1.0 / 8192.0),
                InterpolationPoint::Re(256.0),
            )
        } else {
            (
                InterpolationPoint::Re(-1.0 / 8192.0),
                InterpolationPoint::Re(-256.0),
            )
        };

        let mut path = VecDeque::new();

        path.extend(self.generate_path(pt1).into_iter().rev());
        path.pop_back();
        path.extend(self.generate_path(pt2));

        Self::refine(path, |t, p| {
            let pt = InterpolationPoint::Re(t);
            let w = pt.evaluate(self.consts);
            self.find_point(w, p)
        })
    }

    pub fn contour(&self) -> Vec<Complex64> {
        if !self.valid {
            return vec![];
        }

        let (p1, p2) = match self.pt {
            InterpolationPoint::C(z) => {
                if z.im.abs() < 1.0 / 128.0 {
                    return self.contour_re(z.re);
                } else {
                    log::info!("Can't draw contour for non real C (found {:?})", self.pt);
                    return vec![];
                }
            }
            InterpolationPoint::Re(re) => {
                return self.contour_re(re);
            }
            InterpolationPoint::Xp(p, _) | InterpolationPoint::Xm(p, _) => {
                (p.floor() + 1.0 / 256.0, p.ceil() - 1.0 / 256.0)
            }
        };

        let pt_at = |p| match self.pt {
            InterpolationPoint::C(_) | InterpolationPoint::Re(_) => {
                unreachable!();
            }
            InterpolationPoint::Xp(_, m) => InterpolationPoint::Xp(p, m),
            InterpolationPoint::Xm(_, m) => InterpolationPoint::Xm(p, m),
        };

        let pt1 = pt_at(p1);
        let pt2 = pt_at(p2);

        let mut path = VecDeque::new();

        path.extend(self.generate_path(pt1).into_iter().rev());
        path.pop_back();
        path.extend(self.generate_path(pt2));

        fn zero_asymptote_dist((_, p): (f64, Complex64)) -> f64 {
            if (p.re - p.re.floor()) < 0.5 {
                (p - p.re.floor()).norm_sqr()
            } else {
                (p - p.re.ceil()).norm_sqr()
            }
        }

        fn zero_asymptote_value((_, p): (f64, Complex64)) -> Complex64 {
            if (p.re - p.re.floor()) < 0.5 {
                Complex64::from(p.re.floor())
            } else {
                Complex64::from(p.re.ceil())
            }
        }

        if path.len() < 2 {
            return vec![];
        }

        if (path.front().unwrap().0 - p1).abs() < 1.0 / 128.0
            && zero_asymptote_dist(*path.front().unwrap()) < 1.0 / 64.0
        {
            path.push_front((p1.floor(), zero_asymptote_value(*path.front().unwrap())));
        }

        if (path.back().unwrap().0 - p2).abs() < 1.0 / 128.0
            && zero_asymptote_dist(*path.back().unwrap()) < 1.0 / 64.0
        {
            path.push_back((p2.ceil(), zero_asymptote_value(*path.back().unwrap())));
        }

        Self::refine(path, |t, p| {
            let pt = pt_at(t);
            let w = pt.evaluate(self.consts);
            self.find_point(w, p)
        })
    }

    fn f(&self, z: Complex64) -> Complex64 {
        xp(z, 1.0, self.consts)
    }

    fn df(&self, z: Complex64) -> Complex64 {
        dxp_dp(z, 1.0, self.consts)
    }
}

impl Refiner<f64> for PInterpolatorMut {
    fn mid_point(t1: &f64, t2: &f64) -> f64 {
        (t1 + t2) / 2.0
    }

    fn cmp(t1: &f64, t2: &f64) -> Option<std::cmp::Ordering> {
        t1.partial_cmp(t2)
    }
}

#[derive(Debug)]
pub struct EPInterpolator {
    branch_point_p: Option<Complex64>,
    starting_path_p: Option<Vec<(f64, Complex64)>>,
    p_start: f64,
    consts: CouplingConstants,
}

impl EPInterpolator {
    pub fn new(p_range: i32, consts: CouplingConstants) -> Self {
        let p_start = p_range as f64;
        Self {
            branch_point_p: None,
            starting_path_p: None,
            p_start,
            consts,
        }
    }

    fn cut_x(p: Complex64, im: f64, consts: CouplingConstants) -> Complex64 {
        let sin = (std::f64::consts::PI * p).sin();
        let m_eff = 1.0 + consts.k() as f64 * p;

        let numerator = m_eff + Complex64::i() * im;
        let denominator = 2.0 * consts.h * sin;

        numerator / denominator
    }

    fn cut_xp(p: Complex64, im: f64, consts: CouplingConstants) -> Complex64 {
        Self::cut_x(p, im, consts) * (Complex64::i() * std::f64::consts::PI * p).exp()
    }

    fn cut_xm(p: Complex64, im: f64, consts: CouplingConstants) -> Complex64 {
        Self::cut_x(p, im, consts) * (-Complex64::i() * std::f64::consts::PI * p).exp()
    }

    fn cut_u(p: Complex64, im: f64, consts: CouplingConstants, p_branch: f64) -> Complex64 {
        let xp = Self::cut_xp(p, im, consts);

        let up = xp + 1.0 / xp - 2.0 * consts.kslash() / consts.h * xp.ln();
        let branch_shift = p_branch * consts.k() as f64 * Complex64::i() / consts.h;

        up - Complex64::i() / consts.h - branch_shift
    }

    pub fn get_cut_p(&mut self) -> (Option<Complex64>, Option<Vec<Complex64>>) {
        self.get_cut_f(|p, _, _| p)
    }

    pub fn get_cut_xp(&mut self) -> (Option<Complex64>, Option<Vec<Complex64>>) {
        self.get_cut_f(Self::cut_xp)
    }

    pub fn get_cut_xm(&mut self) -> (Option<Complex64>, Option<Vec<Complex64>>) {
        self.get_cut_f(Self::cut_xm)
    }

    pub fn get_cut_u(&mut self) -> (Option<Complex64>, Option<Vec<Complex64>>) {
        let p_start = self.p_start;
        self.get_cut_f(|p, im, consts| Self::cut_u(p, im, consts, p_start))
    }

    fn get_cut_f(
        &mut self,
        cut_f: impl Fn(Complex64, f64, CouplingConstants) -> Complex64,
    ) -> (Option<Complex64>, Option<Vec<Complex64>>) {
        let consts = self.consts;

        let branch_point_p = self.compute_branch_point_p();
        let Some(starting_path) = self.compute_starting_path_p() else {return (None, None)};

        let branch_point = branch_point_p.map(|p| cut_f(p, 0.0, consts));
        let mut path = VecDeque::new();

        path.extend(
            starting_path
                .iter()
                .rev()
                .map(|(im, p)| ((-im, *p), cut_f(*p, -im, consts))),
        );

        path.pop_back();

        path.extend(
            starting_path
                .iter()
                .map(|(im, p)| ((*im, *p), cut_f(*p, *im, consts))),
        );

        let eval = |(im, p_guess), _| self.find_p_at_im(im, p_guess).map(|p| cut_f(p, im, consts));

        let path = Self::refine(path, eval);

        (branch_point, Some(path))
    }

    fn find_p_at_im(&self, im: f64, guess: Complex64) -> Option<Complex64> {
        nr::find_root(
            |p| en2(p, 1.0, self.consts) + im * im,
            |p| den2_dp(p, 1.0, self.consts),
            guess,
            1.0e-3,
            50,
        )
    }

    fn compute_branch_point_p(&mut self) -> Option<Complex64> {
        if self.branch_point_p.is_none() {
            let im_guess = if self.p_start >= 0.0 { 2.5 } else { -2.5 };
            self.branch_point_p = self.find_p_at_im(0.0, Complex64::new(self.p_start, im_guess));
        }
        self.branch_point_p
    }

    fn compute_starting_path_p(&mut self) -> &Option<Vec<(f64, Complex64)>> {
        if self.starting_path_p.is_some() {
            return &self.starting_path_p;
        }

        let Some(p0) = self.compute_branch_point_p() else { self.starting_path_p = None; return &self.starting_path_p; };

        let mut path = vec![];

        path.push((0.0, p0));
        let mut p_prev = p0;

        for i in 1.. {
            let im = i as f64 * i as f64 * i as f64 / 8.0;

            let Some(p) = self.find_p_at_im(im, p_prev) else {break;};

            path.push((im, p));
            p_prev = p;

            if p.im.abs() > 2.0 {
                break;
            }
        }
        self.starting_path_p = Some(path);
        &self.starting_path_p
    }
}

impl Refiner<(f64, Complex64)> for EPInterpolator {
    fn mid_point(t1: &(f64, Complex64), t2: &(f64, Complex64)) -> (f64, Complex64) {
        ((t1.0 + t2.0) / 2.0, (t1.1 + t2.1) / 2.0)
    }

    fn cmp(t1: &(f64, Complex64), t2: &(f64, Complex64)) -> Option<std::cmp::Ordering> {
        t1.0.partial_cmp(&t2.0)
    }
}
