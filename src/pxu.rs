use crate::kinematics::{den2_dp, du_dp, dxm_dp, dxp_dp, en2, u, xm, xp, CouplingConstants};
use crate::nr::{self};
use crate::pxu2::{PInterpolator, XInterpolator};
use num::complex::Complex;
use num::Zero;

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

// struct Lerp {
//     z0: C,
//     z1: C,
// }

// impl Lerp {
//     fn new(z0: C, z1: C) -> Self {
//         Self { z0, z1 }
//     }
// }

// impl nr::OneParameterFunction for Lerp {
//     fn evaluate(&self, t: f64) -> C {
//         self.z0 + (self.z1 - self.z0) * t
//     }
// }

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

// #[derive(Debug, Clone, Copy, PartialEq)]
// pub enum Component {
//     P,
//     Xp,
//     Xm,
//     U,
// }

// impl PxuPoint {
//     fn get(&self, component: Component) -> C {
//         match component {
//             Component::P => self.p,
//             Component::U => self.u,
//             Component::Xp => self.xp,
//             Component::Xm => self.xm,
//         }
//     }
// }

#[derive(Debug)]
pub struct Grid {
    pub p: Vec<Vec<C>>,
    pub x: Vec<Vec<C>>,
    pub u: Vec<Vec<C>>,
}

impl Grid {
    pub fn new(p_range: i32, consts: CouplingConstants) -> Self {
        // let p = Self::fill_p(p_range, consts);
        let x = Self::fill_x(p_range, consts);
        let u = Self::fill_u(p_range, consts);

        let mut p = vec![];
        for i in 0..=0 {
            p.extend(Self::fill_p(p_range + i, consts));
        }
        Self { p, x, u }
    }

    fn fill_x(p_range: i32, consts: CouplingConstants) -> Vec<Vec<C>> {
        let mut lines = vec![];

        // let p_start = p_range as f64;

        for m in 0..=consts.k() {
            // for m in 1..=1 {
            // let mut xp_points = vec![];

            // if m == 0 && p_range == 0 {
            //     xp_points.push(C::from(consts.s()));
            // } else if m == consts.k() && p_range == -1 {
            //     xp_points.push(C::from(consts.s()));
            // } else if p_range < 0 {
            //     xp_points.push(C::zero());
            // }

            // let steps = 64;

            // for i in 1..=(steps - 1) {
            //     let p = p_start + i as f64 / (steps as f64);
            //     xp_points.push(xp(C::from(p), m as f64, consts));
            // }

            // if m == 0 && p_range == -1 {
            //     xp_points.push(C::from(-1.0 / consts.s()));
            // } else if m == consts.k() && p_range == -2 {
            //     xp_points.push(C::from(-1.0 / consts.s()));
            // } else if p_range < -1 {
            //     xp_points.push(C::zero());
            // }
            let xp_points = XInterpolator::generate_xp(p_range, m as f64, consts);

            lines.push(xp_points.iter().map(|z| z.conj()).collect::<Vec<_>>());
            lines.push(xp_points);
        }

        if p_range == 0 {
            lines.push(vec![C::from(consts.s()), C::from(1000.0)]);
        }

        if p_range == -1 {
            lines.push(vec![C::from(-1000.0), C::from(-1.0 / consts.s())]);
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
            for i in 1..256 {
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
            for i in 1..256 {
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

            let mut cut_xp = vec![C::zero()];
            cut_xp.extend(cut_p.iter().map(|p| xp(*p, 1.0, consts)));
            cut_xp.push(C::zero());

            let cut_xm = cut_p
                .iter()
                .map(|p| xm(*p, 1.0, consts))
                .collect::<Vec<_>>();

            lines.push(cut_xp.iter().map(|z| z.conj()).collect::<Vec<_>>());
            lines.push(cut_xm.iter().map(|z| z.conj()).collect::<Vec<_>>());

            lines.push(cut_xp);
            lines.push(cut_xm);
        }

        lines
    }

    fn fill_u(_p_range: i32, consts: CouplingConstants) -> Vec<Vec<C>> {
        let mut lines = vec![];

        let k = consts.k() as i32;
        for y in (-k)..=(k) {
            let y = y as f64 / consts.h;
            lines.push(vec![C::new(-1000.0, y), C::new(1000.0, y)]);
        }

        lines
    }

    fn fill_p(p_range: i32, consts: CouplingConstants) -> Vec<Vec<C>> {
        let mut lines = vec![];

        let p_start = p_range as f64;
        lines.push(vec![C::from(p_start), C::from(p_start + 1.0)]);

        let xp_func = XpFunc { consts };

        let p_s = {
            let p0 = 1.0 / 8.0;
            let x0 = xp(C::from(p0), 1.0, consts);
            let fixed_re = FixedRe::new(x0.re);
            let start = nr::shoot(&xp_func, &fixed_re, x0.im, 0.0, C::from(p0), x0.im / 1.0)
                .last()
                .unwrap()
                .1;
            let fixed_im = FixedIm::new(0.0);
            let t0 = x0.re;
            let t1 = consts.s();
            nr::shoot(&xp_func, &fixed_im, t0, t1, start, (t1 - t0).abs() / 4.0)
                .last()
                .unwrap()
                .1
        };

        let p_min_one_over_s = {
            let p0 = -1.0 / 64.0;
            let x0 = xp(C::from(p0), 1.0, consts);
            let fixed_re = FixedRe::new(x0.re);
            let pts = nr::shoot(&xp_func, &fixed_re, x0.im, 0.0, C::from(p0), x0.im / 1.0);

            let start = pts.last().unwrap().1;
            let fixed_im = FixedIm::new(0.0);
            let t0 = x0.re;
            let t1 = -1.0 / consts.s();
            let pts = nr::shoot(&xp_func, &fixed_im, t0, t1, start, (t1 - t0).abs() / 4.0);

            pts.last().unwrap().1
        };

        {
            // xp m = 0
            let p0 = p_start + 1.0 / 8.0;
            let xp_fixed_p = XpFixedP::new(p0, consts);

            let (m, p1) = *nr::shoot(&xp_func, &xp_fixed_p, 1.0, 0.0, C::from(p0), 0.5)
                .last()
                .unwrap();

            let fixed_m = XpFixedM::new(m, consts);

            let mut pts = vec![];
            if p_range != 0 {
                pts.push(C::from(p_start));
            } else {
                pts.push(p_s);
            }
            pts.extend(
                nr::shoot_two_sided(
                    &xp_func,
                    &fixed_m,
                    p_start + 1.0 * 1.0 / 64.0,
                    p0,
                    p_start + 63.0 * 1.0 / 64.0,
                    p1,
                    1.0 / 32.0,
                )
                .into_iter()
                .map(|(_, p)| p),
            );
            if p_range != -1 {
                pts.push(C::from(p_start + 1.0));
            } else {
                pts.push(p_min_one_over_s);
            }

            lines.push(pts.iter().map(|z| z.conj()).collect::<Vec<_>>());
            lines.push(pts);
        }

        if p_range == 0 {
            let fixed_im = FixedIm::new(0.0);
            {
                let t0 = consts.s();
                let t1 = 32.0 * consts.s();
                let mut pts = nr::shoot(&xp_func, &fixed_im, t0, t1, p_s, 1.0 / 4.0)
                    .into_iter()
                    .map(|(_, p)| p)
                    .collect::<Vec<_>>();
                pts.push(C::zero());

                lines.push(pts.iter().map(|z| z.conj()).collect::<Vec<_>>());
                lines.push(pts);
            }
        }

        if p_range == -1 {
            let fixed_im = FixedIm::new(0.0);
            {
                let t0 = -1.0 / consts.s();
                let t1 = -32.0 * consts.s();
                let mut pts = nr::shoot(&xp_func, &fixed_im, t0, t1, p_min_one_over_s, 1.0 / 8.0)
                    .into_iter()
                    .map(|(_, p)| p)
                    .collect::<Vec<_>>();
                pts.push(C::zero());
                lines.push(pts.iter().map(|z| z.conj()).collect::<Vec<_>>());
                lines.push(pts);
            }
        }

        if p_range != -1 {
            let p0 = C::from(p_start + 1.0 / 4.0);
            let x0 = xp(p0, 1.0, consts);
            let fixed_re = FixedRe::new(x0.re);
            let pts = nr::shoot(&xp_func, &fixed_re, x0.im, -x0.im, p0, x0.im / 8.0);
            let p1 = pts.last().unwrap().1;
            let xm_fixed_p = XmFixedP::new(p0.re, consts);

            let starts = nr::shoot_two_sided(
                &xp_func,
                &xm_fixed_p,
                0.0,
                1.0,
                1.0f64.max(consts.k() as f64 - 2.0),
                p1,
                1.0,
            );
            for (m, p) in starts.into_iter() {
                let xm_fixed_m = XmFixedM::new(m, consts);
                let mut pts = vec![];
                if m == 0.0 && p_range == 0 {
                    pts.push(p_s);
                } else {
                    pts.push(C::from(p_start));
                }
                pts.extend(
                    nr::shoot_two_sided(
                        &xp_func,
                        &xm_fixed_m,
                        p_start + 1.0 / 128.0,
                        p0.re,
                        p_start + 1.0 - 1.0 / 128.0,
                        p,
                        1.0 / 128.0,
                    )
                    .into_iter()
                    .map(|(_, p)| p),
                );
                pts.push(C::from(p_start));

                lines.push(pts.iter().map(|z| z.conj()).collect::<Vec<_>>());
                lines.push(pts);
            }
        }

        if p_range == -1 {
            let p0 = C::from(-1.0 / 64.0);
            let x0 = xp(p0, 1.0, consts);
            let fixed_re = FixedRe::new(x0.re);
            let pts = nr::shoot(&xp_func, &fixed_re, x0.im, -x0.im, p0, x0.im / 8.0);
            let p1 = pts.last().unwrap().1;

            let xm_fixed_p = XmFixedP::new(p0.re, consts);

            // let starts = nr::shoot(&xp_func, &xm_fixed_p, 1.0, 0.0, p1, 1.0);
            let starts = nr::shoot(&xp_func, &xm_fixed_p, 1.0, (consts.k() - 1) as f64, p1, 1.0);
            let starts = [*starts.last().unwrap()];
            for (m, p) in starts.into_iter() {
                let xm_fixed_m = XmFixedM::new(m, consts);
                let mut pts = vec![];
                if m == 0.0 && p_range == 0 {
                    // pts.push(p_s);
                } else {
                    // pts.push(C::from(p_start));
                }
                pts.extend(
                    nr::shoot_two_sided(
                        &xp_func,
                        &xm_fixed_m,
                        p_start + 1.0 / 32.0,
                        p0.re,
                        p_start + 1.0 - 1.0 / 128.0,
                        p,
                        1.0 / 128.0,
                    )
                    .into_iter()
                    .map(|(_, p)| p),
                );

                // pts.push(C::from(p_start));

                lines.push(pts.iter().map(|z| z.conj()).collect::<Vec<_>>());
                lines.push(pts);
            }
        }

        {
            let p0 = nr::find_root(
                |p| en2(p, 1.0, consts) - C::new(0.0, 0.0),
                |p| den2_dp(p, 1.0, consts),
                C::new(0.0, 2.5),
                // C::new(0.0, 0.5),
                1.0e-3,
                50,
            );
            // log::info!("{:?}", p0);

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

            let cut_p = cut_p.into_iter().map(|p| p).collect::<Vec<_>>();

            lines.push(cut_p);
        }

        lines
    }
}

// struct SheetData {
//     p_branch: i32,
//     xp_branch: i32,
//     xm_branch: i32,
// }

// struct GridLine {}

pub struct Cut {
    pub p: Vec<Vec<C>>,
    pub x: Vec<Vec<C>>,
    pub u: Vec<Vec<C>>,

    pub branch_points_p: Vec<C>,
    pub branch_points_x: Vec<C>,
    pub branch_points_u: Vec<C>,
}

impl Cut {
    pub fn get(p_range: i32, consts: CouplingConstants) -> Vec<Self> {
        vec![
            Self::x(p_range, consts),
            Self::e(p_range, consts),
            Self::x_log(p_range, consts),
        ]
    }

    fn x(p_range: i32, consts: CouplingConstants) -> Self {
        let p_start = p_range as f64;
        let xp_func = XpFunc { consts };

        let p_s = {
            let p0 = 1.0 / 8.0;
            let x0 = xp(C::from(p0), 1.0, consts);
            let fixed_re = FixedRe::new(x0.re);
            let start = nr::shoot(&xp_func, &fixed_re, x0.im, 0.0, C::from(p0), x0.im / 1.0)
                .last()
                .unwrap()
                .1;
            let fixed_im = FixedIm::new(0.0);
            let t0 = x0.re;
            let t1 = consts.s();
            nr::shoot(&xp_func, &fixed_im, t0, t1, start, (t1 - t0).abs() / 4.0)
                .last()
                .unwrap()
                .1
        };

        let p_min_one_over_s = {
            let p0 = -1.0 / 64.0;
            let x0 = xp(C::from(p0), 1.0, consts);
            let fixed_re = FixedRe::new(x0.re);
            let pts = nr::shoot(&xp_func, &fixed_re, x0.im, 0.0, C::from(p0), x0.im / 1.0);

            let start = pts.last().unwrap().1;
            let fixed_im = FixedIm::new(0.0);
            let t0 = x0.re;
            let t1 = -1.0 / consts.s();
            let pts = nr::shoot(&xp_func, &fixed_im, t0, t1, start, (t1 - t0).abs() / 4.0);

            pts.last().unwrap().1
        };

        let mut p_points = vec![];
        {
            if p_range != -1 {
                let p0 = C::from(p_start + 1.0 / 4.0);
                let x0 = xp(p0, 1.0, consts);
                let fixed_re = FixedRe::new(x0.re);
                let p1 = nr::shoot(&xp_func, &fixed_re, x0.im, -x0.im, p0, x0.im / 8.0)
                    .last()
                    .unwrap()
                    .1;
                let xm_fixed_p = XmFixedP::new(p0.re, consts);

                let p2 = nr::shoot(&xp_func, &xm_fixed_p, 1.0, 0.0, p1, 1.0)
                    .last()
                    .unwrap()
                    .1;

                let xm_fixed_m = XmFixedM::new(0.0, consts);
                p_points.push(C::from(p_start));
                p_points.extend(
                    nr::shoot_two_sided(
                        &xp_func,
                        &xm_fixed_m,
                        p_start + 1.0 / 128.0,
                        p0.re,
                        p_start + 1.0 - 1.0 / 128.0,
                        p2,
                        1.0 / 128.0,
                    )
                    .into_iter()
                    .rev()
                    .map(|(_, p)| p),
                );
                if p_range == 0 {
                    p_points.push(p_s);
                } else {
                    p_points.push(C::from(p_start));
                }
            }

            {
                // xp m = 0
                let p0 = p_start + 1.0 / 8.0;
                let xp_fixed_p = XpFixedP::new(p0, consts);

                let (m, p1) = *nr::shoot(&xp_func, &xp_fixed_p, 1.0, 0.0, C::from(p0), 0.5)
                    .last()
                    .unwrap();

                let fixed_m = XpFixedM::new(m, consts);

                p_points.extend(
                    nr::shoot_two_sided(
                        &xp_func,
                        &fixed_m,
                        p_start + 1.0 * 1.0 / 64.0,
                        p0,
                        p_start + 63.0 * 1.0 / 64.0,
                        p1,
                        1.0 / 512.0,
                    )
                    .into_iter()
                    .map(|(_, p)| p),
                );
                if p_range != -1 {
                    p_points.push(C::from(p_start + 1.0));
                } else {
                    p_points.push(p_min_one_over_s);
                }
            }
        }

        let mut x_points = vec![];

        if p_range == 0 {
            x_points.push(C::from(consts.s()));
        } else if p_range < 0 {
            x_points.push(C::zero());
        }

        let steps = 64;

        for i in 1..=(steps - 1) {
            let p = p_start + i as f64 / (steps as f64);
            x_points.push(xp(C::from(p), 0.0, consts));
        }

        if p_range == -1 {
            x_points.push(C::from(-1.0 / consts.s()));
        } else if p_range < -1 {
            x_points.push(C::zero());
        }

        let mut u_points = vec![];

        let u_s = {
            let s = consts.s();
            s + 1.0 / s - 2.0 * consts.kslash() / consts.h * s.ln()
        };
        if p_range == 0 {
            u_points.push(C::new(-100.0, -1.0 / consts.h));
            u_points.push(C::new(u_s, -1.0 / consts.h));
        } else if p_range == -1 {
            u_points.push(C::new(-u_s, -1.0 / consts.h));
            u_points.push(C::new(100.0, -1.0 / consts.h));
        } else {
            u_points.push(C::new(-100.0, -1.0 / consts.h));
            u_points.push(C::new(100.0, -1.0 / consts.h));
        }

        let x = vec![
            x_points.iter().rev().map(|z| *z).collect::<Vec<_>>(),
            x_points.iter().map(|z| z.conj()).collect::<Vec<_>>(),
        ];

        let p = vec![
            p_points.iter().map(|z| z.conj()).collect::<Vec<_>>(),
            p_points,
        ];
        let u = vec![
            u_points.iter().map(|z| z.conj()).collect::<Vec<_>>(),
            u_points,
        ];

        let branch_points_p = match p_range {
            0 => vec![p_s, p_s.conj()],
            -1 => vec![p_min_one_over_s, p_min_one_over_s.conj()],
            _ => vec![],
        };
        let branch_points_x = match p_range {
            0 => vec![C::from(consts.s())],
            -1 => vec![C::from(-1.0 / consts.s())],
            _ => vec![],
        };
        let branch_points_u = match p_range {
            0 => vec![C::new(u_s, 1.0 / consts.h), C::new(u_s, -1.0 / consts.h)],
            -1 => vec![C::new(-u_s, 1.0 / consts.h), C::new(-u_s, -1.0 / consts.h)],
            _ => vec![],
        };

        Self {
            p,
            x,
            u,
            branch_points_p,
            branch_points_x,
            branch_points_u,
        }
    }

    fn e(p_range: i32, consts: CouplingConstants) -> Self {
        let p0 = nr::find_root(
            |p| en2(p, 1.0, consts) - C::new(0.0, 0.0),
            |p| den2_dp(p, 1.0, consts),
            C::new(0.0, 2.5),
            // C::new(0.0, 0.5),
            1.0e-3,
            50,
        )
        .unwrap();

        let mut cut_p = vec![p0];
        for i in 0..=256 {
            let im = i as f64 * i as f64 / 64.0;

            let p = nr::find_root(
                |p| en2(p, 1.0, consts) - C::new(-im, 0.0),
                |p| den2_dp(p, 1.0, consts),
                *cut_p.last().unwrap(),
                1.0e-3,
                50,
            );

            cut_p.push(p.unwrap());
        }

        let mut xp_points = vec![];

        xp_points.push(C::zero());
        xp_points.extend(
            cut_p
                .iter()
                .rev()
                .map(|p| xp(*p + C::from(1.0e-5), 1.0, consts)),
        );
        xp_points.extend(cut_p.iter().map(|p| xp(*p + C::from(-1.0e-5), 1.0, consts)));
        xp_points.push(C::zero());

        let mut xm_points = vec![];

        xm_points.extend(
            cut_p
                .iter()
                .rev()
                .map(|p| xm(*p + C::from(1.0e-5), 1.0, consts)),
        );
        xm_points.extend(cut_p.iter().map(|p| xm(*p + C::from(-1.0e-5), 1.0, consts)));

        let mut u_points = vec![];

        u_points.extend(
            cut_p
                .iter()
                .rev()
                .map(|p| u(*p + C::from(1.0e-5), p_range, consts)),
        );
        u_points.extend(
            cut_p
                .iter()
                .map(|p| u(*p + C::from(-1.0e-5), p_range, consts)),
        );

        let p = vec![cut_p.iter().map(|z| z.conj()).collect::<Vec<_>>(), cut_p];
        let x = vec![
            xp_points.iter().map(|z| z.conj()).collect::<Vec<_>>(),
            xp_points,
            xm_points.iter().map(|z| z.conj()).collect::<Vec<_>>(),
            xm_points,
        ];
        let u_cuts = vec![
            u_points.iter().map(|z| z.conj()).collect::<Vec<_>>(),
            u_points,
        ];

        let branch_points_p = vec![p0, p0.conj()];
        let branch_points_x = vec![
            xp(p0, 1.0, consts),
            xm(p0, 1.0, consts),
            xp(p0, 1.0, consts).conj(),
            xm(p0, 1.0, consts).conj(),
        ];
        let branch_points_u = vec![u(p0, p_range, consts), u(p0, p_range, consts).conj()];

        Self {
            p,
            x,
            u: u_cuts,
            branch_points_p,
            branch_points_x,
            branch_points_u,
        }
    }

    fn x_log(p_range: i32, consts: CouplingConstants) -> Self {
        // let x_points;
        let branch_points_x;
        if p_range == 0 {
            //     x_points = vec![C::from(-100.0), C::zero()];
            branch_points_x = vec![C::zero()];
        } else if p_range == -1 {
            //     x_points = vec![C::zero(), C::from(100.0)];
            branch_points_x = vec![C::zero()];
        } else {
            //     x_points = vec![C::from(-100.0), C::from(100.0)];
            branch_points_x = vec![]
        }

        let mut p = vec![];

        let mut x = vec![];

        if p_range == 0 {
            let p0 = p_range as f64 + 1.0 / 32.0;
            let p1 = p_range as f64 + 24.0 / 32.0;

            let mut p_int = PInterpolator::xp(p0, consts)
                .goto_xm(p0, 1.0)
                .goto_xm(p1, 1.0)
                .goto_im(0.0);

            p_int.p_path = vec![*p_int.p_path.last().unwrap()];
            p_int.x_path = vec![*p_int.x_path.last().unwrap()];

            let p_int2 = p_int.clone().goto_re(-1000.0);
            x.push(p_int2.x_path);
            p.push(p_int2.p_path.iter().map(|p| p.conj()).collect());
            p.push(p_int2.p_path);

            let p_int2 = p_int.goto_re(0.01);
            x.push(p_int2.x_path);
            p.push(p_int2.p_path.iter().map(|p| p.conj()).collect());
            p.push(p_int2.p_path);

            let mut p_int = PInterpolator::xp(p1, consts).goto_im(0.0);

            p_int.p_path = vec![*p_int.p_path.last().unwrap()];
            p_int.x_path = vec![*p_int.x_path.last().unwrap()];

            let p_int2 = p_int.clone().goto_re(-1000.0);
            x.push(p_int2.x_path);
            p.push(p_int2.p_path.iter().map(|p| p.conj()).collect());
            p.push(p_int2.p_path);

            let p_int2 = p_int.goto_re(0.01);
            x.push(p_int2.x_path);
            p.push(p_int2.p_path.iter().map(|p| p.conj()).collect());
            p.push(p_int2.p_path);
        }

        if p_range > 0 {
            let p0 = p_range as f64 + 8.0 / 32.0;
            let p1 = p_range as f64 + 24.0 / 32.0;

            let mut p_int = PInterpolator::xp(p0, consts).goto_xp(p1, 1.0);

            // p_int.p_path = vec![*p_int.p_path.last().unwrap()];
            // p_int.x_path = vec![*p_int.x_path.last().unwrap()];

            let p_int2 = p_int.clone();
            x.push(p_int2.x_path);
            p.push(p_int2.p_path.iter().map(|p| p.conj()).collect());
            p.push(p_int2.p_path);

            // let p_int2 = p_int.goto_re(0.01);
            // x.push(p_int2.x_path);
            // p.push(p_int2.p_path.iter().map(|p| p.conj()).collect());
            // p.push(p_int2.p_path);
        }

        let u = vec![];

        let branch_points_p = vec![];
        let branch_points_u = vec![];

        Self {
            p,
            x,
            u,
            branch_points_p,
            branch_points_x,
            branch_points_u,
        }
    }
}

#[derive(Debug, Clone)]
pub struct PxuPoint {
    pub p: C,
    pub xp: C,
    pub xm: C,
    pub u: C,
    pub consts: CouplingConstants,
}

impl PxuPoint {
    pub fn new(p: C, consts: CouplingConstants) -> Self {
        let xp = xp(p, 1.0, consts);
        let xm = xm(p, 1.0, consts);
        let u = u(p, p.re.floor() as i32, consts);
        Self {
            p,
            xp,
            xm,
            u,
            consts,
        }
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

    pub fn shift_xp(&self, new_xp: C) -> Option<Self> {
        let p = nr::find_root(
            |p| xp(p, 1.0, self.consts) - new_xp,
            |p| dxp_dp(p, 1.0, self.consts),
            self.p,
            1.0e-6,
            50,
        );
        self.limit_p(p, self.consts)
    }

    pub fn shift_xm(&self, new_xm: C) -> Option<Self> {
        let p = nr::find_root(
            |p| xm(p, 1.0, self.consts) - new_xm,
            |p| dxm_dp(p, 1.0, self.consts),
            self.p,
            1.0e-6,
            50,
        );
        self.limit_p(p, self.consts)
    }

    pub fn shift_u(&self, new_u: C) -> Option<Self> {
        let p = nr::find_root(
            |p| u(p, p.re.floor() as i32, self.consts) - new_u,
            |p| du_dp(p, self.consts),
            self.p,
            1.0e-6,
            50,
        );
        self.limit_p(p, self.consts)
    }
}
