use crate::kinematics::{den2_dp, du_dp, dxm_dp, dxp_dp, en2, u, xm, xp, CouplingConstants};
use crate::nr::{self, OneParameterFunction};
use num::complex::Complex;
use num::Zero;
use std::f64::consts::PI;

type C = Complex<f64>;

struct XpFixedM {
    m: f64,
    consts: CouplingConstants,
}

impl XpFixedM {
    fn new(m: f64, consts: CouplingConstants) -> Self {
        Self { m, consts }
    }
}

impl nr::OneParameterFunction for XpFixedM {
    fn evaluate(&self, p: f64) -> C {
        xp(C::from(p), self.m, self.consts)
    }
}

struct XmFixedM {
    m: f64,
    consts: CouplingConstants,
}

impl XmFixedM {
    fn new(m: f64, consts: CouplingConstants) -> Self {
        Self { m, consts }
    }
}

impl nr::OneParameterFunction for XmFixedM {
    fn evaluate(&self, p: f64) -> C {
        xm(C::from(p), self.m, self.consts)
    }
}

struct XpFixedP {
    p: f64,
    consts: CouplingConstants,
}

impl XpFixedP {
    fn new(p: f64, consts: CouplingConstants) -> Self {
        Self { p, consts }
    }
}

impl nr::OneParameterFunction for XpFixedP {
    fn evaluate(&self, m: f64) -> C {
        xp(C::from(self.p), m, self.consts)
    }
}

struct XmFixedP {
    p: f64,
    consts: CouplingConstants,
}

impl XmFixedP {
    fn new(p: f64, consts: CouplingConstants) -> Self {
        Self { p, consts }
    }
}

impl nr::OneParameterFunction for XmFixedP {
    fn evaluate(&self, m: f64) -> C {
        xm(C::from(self.p), m, self.consts)
    }
}

struct FixedRe {
    re: f64,
}

impl FixedRe {
    fn new(re: f64) -> Self {
        Self { re }
    }
}

impl nr::OneParameterFunction for FixedRe {
    fn evaluate(&self, im: f64) -> C {
        C::new(self.re, im)
    }
}

struct FixedIm {
    im: f64,
}

impl FixedIm {
    fn new(im: f64) -> Self {
        Self { im }
    }
}

impl nr::OneParameterFunction for FixedIm {
    fn evaluate(&self, re: f64) -> C {
        C::new(re, self.im)
    }
}

struct Lerp {
    z0: C,
    z1: C,
}

impl Lerp {
    fn new(z0: C, z1: C) -> Self {
        Self { z0, z1 }
    }
}

impl nr::OneParameterFunction for Lerp {
    fn evaluate(&self, t: f64) -> C {
        self.z0 + (self.z1 - self.z0) * t
    }
}

struct XpFunc {
    consts: CouplingConstants,
}

impl nr::Func for XpFunc {
    fn f(&self, z: C) -> C {
        xp(z, 1.0, self.consts)
    }
    fn df(&self, z: C) -> C {
        dxp_dp(z, 1.0, self.consts)
    }
}

struct XmFunc {
    consts: CouplingConstants,
}

impl nr::Func for XmFunc {
    fn f(&self, z: C) -> C {
        xm(z, 1.0, self.consts)
    }
    fn df(&self, z: C) -> C {
        dxm_dp(z, 1.0, self.consts)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Component {
    P,
    Xp,
    Xm,
    U,
}

impl PxuPoint {
    fn get(&self, component: Component) -> C {
        match component {
            Component::P => self.p,
            Component::U => self.u,
            Component::Xp => self.xp,
            Component::Xm => self.xm,
        }
    }
}

#[derive(Debug)]
pub struct PxuGrid {
    pub points: Vec<Vec<PxuPoint>>,
}

#[derive(Debug)]
pub struct Grid {
    pub p: Vec<Vec<C>>,
    pub x: Vec<Vec<C>>,
    pub u: Vec<Vec<C>>,
}

impl Grid {
    pub fn new(p_range: i32, consts: CouplingConstants) -> Self {
        let p = Self::fill_p(p_range, consts);
        let x = Self::fill_x(p_range, consts);
        let u = Self::fill_u(p_range, consts);

        Self { p, x, u }
    }

    fn fill_x(p_range: i32, consts: CouplingConstants) -> Vec<Vec<C>> {
        let mut lines = vec![];

        let p_start = p_range as f64 * 2.0 * PI;

        for m in 1..consts.k() {
            let mut xp_points = vec![];
            let mut xm_points = vec![];

            if p_range < 0 {
                xp_points.push(C::zero());
                xm_points.push(C::zero());
            }

            let steps = 64;

            for i in 1..=(steps - 1) {
                let p = p_start + 2.0 * PI * i as f64 / (steps as f64);
                xp_points.push(xp(C::from(p), m as f64, consts));
                xm_points.push(xm(C::from(p), m as f64, consts));
            }

            if p_range < -1 {
                xp_points.push(C::zero());
                xm_points.push(C::zero());
            }

            lines.push(xp_points);
            lines.push(xm_points);
        }

        if p_range == 0 {
            lines.push(vec![C::from(consts.s()), C::from(1000.0)]);
        }

        if p_range == -1 {
            lines.push(vec![C::from(-1000.0), C::from(-1.0 / consts.s())]);
        }
        lines
    }

    fn fill_u(p_range: i32, consts: CouplingConstants) -> Vec<Vec<C>> {
        let mut lines = vec![];

        let k = consts.k() as i32;
        for y in (-k + 1)..=(k - 1) {
            lines.push(vec![C::new(-1000.0, y as f64), C::new(1000.0, y as f64)]);
        }

        lines
    }

    fn fill_p(p_range: i32, consts: CouplingConstants) -> Vec<Vec<C>> {
        let mut lines = vec![];

        lines
    }
}

struct SheetData {
    p_branch: i32,
    xp_branch: i32,
    xm_branch: i32,
}

struct GridLine {}

struct Cut {}

impl PxuGrid {
    pub fn new_pm(consts: CouplingConstants) -> Self {
        let p_range = 0;
        let p_start = 2. * PI * p_range as f64;

        let p0 = C::from(p_start + PI / 6.0);

        let mut points = vec![];

        let xp_func = XpFunc { consts };
        let xm_func = XmFunc { consts };
        {
            let p0 = p_start + PI / 4.0;
            let xp_fixed_p = XpFixedP::new(p0, consts);

            let starts =
                nr::shoot_two_sided(&xp_func, &xp_fixed_p, 0.5, 1.0, 1.5, C::from(p0), 0.5);

            for (m, p1) in starts {
                let fixed_m = XpFixedM::new(m, consts);

                let pts = nr::shoot_two_sided(
                    &xp_func,
                    &fixed_m,
                    p_start + 1.0 * PI / 64.0,
                    p0,
                    p_start + 63.0 * PI / 32.0,
                    p1,
                    PI / 128.0,
                )
                .into_iter()
                .map(|(_, p)| PxuPoint::new(p, consts))
                .collect::<Vec<_>>();
                points.push(pts);
            }
        }

        {
            let p0 = p_start + PI / 4.0;
            let xm_fixed_p = XmFixedP::new(p0, consts);

            let starts = nr::shoot(&xm_func, &xm_fixed_p, 1.0, 0.0, C::from(p0), 0.5);

            let (m, p1) = *starts.last().unwrap();

            {
                let fixed_m = XmFixedM::new(m, consts);

                let pts = nr::shoot_two_sided(
                    &xm_func,
                    &fixed_m,
                    p_start + 1.0 * PI / 64.0,
                    p0,
                    p_start + 63.0 * PI / 32.0,
                    p1,
                    PI / 128.0,
                )
                .into_iter()
                .map(|(_, p)| PxuPoint::new(p, consts))
                .collect::<Vec<_>>();
                points.push(pts);
            }
        }

        if p_range == 0 {
            let x0 = xp(p0, 1.0, consts);
            let fixed_re = FixedRe::new(x0.re);
            let start = nr::shoot(&xp_func, &fixed_re, x0.im, 0.0, p0, x0.im / 1.0)
                .last()
                .unwrap()
                .1;
            let fixed_im = FixedIm::new(0.0);
            let t0 = consts.s();
            let t1 = x0.re;
            let t2 = 8.0 * consts.s();
            let pts = nr::shoot_two_sided(
                &xp_func,
                &fixed_im,
                t0,
                t1,
                t2,
                start,
                (t1 - t0).abs() / 4.0,
            )
            .into_iter()
            .map(|(_, p)| PxuPoint::new(p, consts))
            .collect::<Vec<_>>();

            points.push(pts);
        }

        if p_range == 0 {
            let p0 = C::from(p_start + 0.1 * PI);
            let x0 = xp(p0, 1.0, consts);
            let fixed_re = FixedRe::new(x0.re);

            let pts = nr::shoot(&xp_func, &fixed_re, x0.im, -x0.im, p0, x0.im / 8.0);
            let p1 = pts.last().unwrap().1;

            points.push(
                pts.into_iter()
                    .map(|(_, p)| PxuPoint::new(p, consts))
                    .collect::<Vec<_>>(),
            );

            let xm_fixed_p = XmFixedP::new(p0.re, consts);
            let pts = nr::shoot(&xp_func, &xm_fixed_p, 1.0, 0.0, p1, 0.5);
            let p3 = pts.last().unwrap().1;

            points.push(
                pts.into_iter()
                    .map(|(_, p)| PxuPoint::new(p, consts))
                    .collect::<Vec<_>>(),
            );

            let xm_fixed_m = XmFixedM::new(0.0, consts);
            let pts = nr::shoot_two_sided(
                &xp_func,
                &xm_fixed_m,
                p_start + PI / 32.0,
                p0.re,
                p_start + 2.0 * PI - PI / 16.0,
                p3,
                PI / 32.0,
            )
            .into_iter()
            .map(|(_, p)| PxuPoint::new(p, consts))
            .collect::<Vec<_>>();

            points.push(pts);
        }

        {
            // let p0 = C::from(p_start + 5.5);
            let p0 = C::from(p_start + 1.0);
            let x0 = xp(p0, 1.0, consts);
            let fixed_re = FixedRe::new(x0.re);
            // let pts = nr::shoot(&xp_func, &fixed_re, x0.im, -x0.im, p0, x0.im / 8.0);
            let pts = nr::shoot(&xp_func, &fixed_re, x0.im, -x0.im, p0, x0.im / 8.0);
            let p1 = pts.last().unwrap().1;
            let pts = pts
                .into_iter()
                .map(|(_, p)| PxuPoint::new(p, consts))
                .collect::<Vec<_>>();

            // points.push(pts);

            let xm_fixed_p = XmFixedP::new(p0.re, consts);

            let starts = nr::shoot(&xp_func, &xm_fixed_p, 1.0, 4.0, p1, 1.0);
            for (m, p) in starts {
                let xm_fixed_m = XmFixedM::new(m, consts);
                let pts = nr::shoot_two_sided(
                    &xp_func,
                    &xm_fixed_m,
                    p_start + PI / 32.0,
                    p0.re,
                    p_start + 2.0 * PI - PI / 16.0,
                    p,
                    PI / 32.0,
                )
                .into_iter()
                .map(|(_, p)| PxuPoint::new(p, consts))
                .collect::<Vec<_>>();

                points.push(pts);
            }
        }

        {
            let p0 = nr::find_root(
                |p| en2(p, 1.0, consts) - C::new(0.0, 0.0),
                |p| den2_dp(p, 1.0, consts),
                C::new(0.0, 0.5),
                1.0e-3,
                50,
            );

            let mut cut_p = vec![p0.unwrap()];
            for i in 1..64 {
                let im = i as f64 * i as f64 / 64.0;

                let p = nr::find_root(
                    |p| en2(p, 1.0, consts) - C::new(-im, 0.001),
                    |p| den2_dp(p, 1.0, consts),
                    *cut_p.last().unwrap(),
                    1.0e-3,
                    50,
                );

                cut_p.push(p.unwrap());
            }

            cut_p.reverse();
            for i in 1..64 {
                let im = i as f64 * i as f64 / 64.0;

                let p = nr::find_root(
                    |p| en2(p, 1.0, consts) - C::new(-im, -0.001),
                    |p| den2_dp(p, 1.0, consts),
                    *cut_p.last().unwrap(),
                    1.0e-3,
                    50,
                );

                cut_p.push(p.unwrap());
            }

            let cut_p = cut_p
                .into_iter()
                .map(|p| PxuPoint::new(p, consts))
                .collect::<Vec<_>>();

            points.push(cut_p);
        }

        Self { points }
    }
}

#[derive(Debug, Clone)]
pub struct PxuPoint {
    pub p: C,
    pub xp: C,
    pub xm: C,
    pub u: C,
}

impl PxuPoint {
    pub fn new(p: C, consts: CouplingConstants) -> Self {
        let xp = xp(p, 1.0, consts);
        let xm = xm(p, 1.0, consts);
        let u = u(p, consts);
        Self { p, xp, xm, u }
    }

    fn limit_p(&self, p: Option<C>, consts: CouplingConstants) -> Option<Self> {
        if let Some(p) = p {
            if (self.p - p).norm_sqr() > 4.0 {
                None
            } else {
                Some(Self::new(p, consts))
            }
        } else {
            None
        }
    }

    pub fn shift_xp(&self, new_xp: C, consts: CouplingConstants) -> Option<Self> {
        let p = nr::find_root(
            |p| xp(p, 1.0, consts) - new_xp,
            |p| dxp_dp(p, 1.0, consts),
            self.p,
            1.0e-6,
            50,
        );
        self.limit_p(p, consts)
    }

    pub fn shift_xm(&self, new_xm: C, consts: CouplingConstants) -> Option<Self> {
        let p = nr::find_root(
            |p| xm(p, 1.0, consts) - new_xm,
            |p| dxm_dp(p, 1.0, consts),
            self.p,
            1.0e-6,
            50,
        );
        self.limit_p(p, consts)
    }

    pub fn shift_u(&self, new_u: C, consts: CouplingConstants) -> Option<Self> {
        let p = nr::find_root(
            |p| u(p, consts) - new_u,
            |p| du_dp(p, consts),
            self.p,
            1.0e-6,
            50,
        );
        self.limit_p(p, consts)
    }
}
