use crate::kinematics::{dxm_dp, dxp_dp, u, xm, xp, CouplingConstants};
use crate::newton_raphson::find_root;
use num::complex::Complex;
use std::f64::consts::PI;

type C = Complex<f64>;

#[derive(Debug)]
pub struct PxuGrid {
    pub points: Vec<Vec<PxuPoint>>,
}

impl PxuGrid {
    pub fn new_pm(consts: CouplingConstants) -> Self {
        let p0 = C::new(PI / 6.0, 0.0);
        let pxu_real = PxuPoint::new(p0, consts);

        let mut points = vec![];

        for i in 0..=3 {
            let m = i as f64 / 2.0;

            let dp = PI / 128.0;

            let pxu1 = pxu_real.shift_xp(xp(p0, m, consts), consts).unwrap();
            let mut left_points = vec![pxu1.clone()];

            for j in 1..(((2.0 * PI - dp) - p0.re) / dp).floor() as i32 {
                println!("{m} {i} L");
                let p = p0 + j as f64 * dp;
                let pxu = left_points
                    .last()
                    .unwrap()
                    .shift_xp(xp(p, m, consts), consts);
                if pxu.is_none() {
                    break;
                }
                left_points.push(pxu.unwrap());
            }

            let mut right_points = vec![pxu1.clone()];

            for j in 1..((p0.re - (0.0 + dp)) / dp).floor() as i32 {
                println!("{m} {i} R");

                let p = p0 - j as f64 * dp;
                let pxu = right_points
                    .last()
                    .unwrap()
                    .shift_xp(xp(p, m, consts), consts);
                if pxu.is_none() {
                    break;
                }
                right_points.push(pxu.unwrap());
            }

            left_points.reverse();
            left_points.pop();

            left_points.extend(right_points.into_iter());

            points.push(left_points);
        }

        for i in 0..=0 {
            let m = i as f64 / 2.0;

            let dp = PI / 128.0;

            let pxu1 = pxu_real.shift_xm(xm(p0, m, consts), consts).unwrap();
            let mut left_points = vec![pxu1.clone()];

            for j in 1..(((2.0 * PI - dp) - p0.re) / dp).floor() as i32 {
                println!("{m} {i} L");
                let p = p0 + j as f64 * dp;
                let pxu = left_points
                    .last()
                    .unwrap()
                    .shift_xm(xm(p, m, consts), consts);
                if pxu.is_none() {
                    break;
                }
                left_points.push(pxu.unwrap());
            }

            let mut right_points = vec![pxu1.clone()];

            for j in 1..((p0.re - (0.0 + dp)) / dp).floor() as i32 {
                println!("{m} {i} R");

                let p = p0 - j as f64 * dp;
                let pxu = right_points
                    .last()
                    .unwrap()
                    .shift_xm(xm(p, m, consts), consts);
                if pxu.is_none() {
                    break;
                }
                right_points.push(pxu.unwrap());
            }

            left_points.reverse();
            left_points.pop();

            left_points.extend(right_points.into_iter());

            points.push(left_points);
        }

        for i in (0..=14).step_by(2) {
            let m = i as f64 / 2.0;

            let dp = PI / 128.0;

            let mut pxu0 = pxu_real.shift_xp(xp(p0, 0.0, consts), consts).unwrap();
            let mut n = 0.0;

            while n < m {
                n += 0.25;
                pxu0 = pxu0.shift_xp(xm(p0, n, consts), consts).unwrap();
            }

            let pxu1 = pxu0.shift_xp(xm(p0, m, consts), consts).unwrap();
            let mut left_points = vec![pxu1.clone()];

            for j in 1..(((2.0 * PI - dp) - p0.re) / dp).floor() as i32 {
                println!("{m} {i} L");
                let p = p0 + j as f64 * dp;
                let pxu = left_points
                    .last()
                    .unwrap()
                    .shift_xp(xm(p, m, consts), consts);
                if pxu.is_none() {
                    break;
                }
                left_points.push(pxu.unwrap());
            }

            let mut right_points = vec![pxu1.clone()];

            for j in 1..((p0.re - (0.0 + dp)) / dp).floor() as i32 {
                println!("{m} {i} R");

                let p = p0 - j as f64 * dp;
                let pxu = right_points
                    .last()
                    .unwrap()
                    .shift_xp(xm(p, m, consts), consts);
                if pxu.is_none() {
                    break;
                }
                right_points.push(pxu.unwrap());
            }

            left_points.reverse();
            left_points.pop();

            left_points.extend(right_points.into_iter());

            points.push(left_points);
        }

        for i in (0..=4).step_by(2) {
            let m = i as f64 / 2.0;

            let dp = PI / 128.0;

            let mut pxu0 = pxu_real.shift_xm(xm(p0, 0.0, consts), consts).unwrap();
            let mut n = 0.0;

            while n < m {
                n += 0.25;
                pxu0 = pxu0.shift_xm(xp(p0, n, consts), consts).unwrap();
            }

            let pxu1 = pxu0.shift_xm(xp(p0, m, consts), consts).unwrap();
            let mut left_points = vec![pxu1.clone()];

            for j in 1..(((2.0 * PI - dp) - p0.re) / dp).floor() as i32 {
                println!("{m} {i} L");
                let p = p0 + j as f64 * dp;
                let pxu = left_points
                    .last()
                    .unwrap()
                    .shift_xm(xp(p, m, consts), consts);
                if pxu.is_none() {
                    break;
                }
                left_points.push(pxu.unwrap());
            }

            let mut right_points = vec![pxu1.clone()];

            for j in 1..((p0.re - (0.0 + dp)) / dp).floor() as i32 {
                println!("{m} {i} R");

                let p = p0 - j as f64 * dp;
                let pxu = right_points
                    .last()
                    .unwrap()
                    .shift_xm(xp(p, m, consts), consts);
                if pxu.is_none() {
                    break;
                }
                right_points.push(pxu.unwrap());
            }

            left_points.reverse();
            left_points.pop();

            left_points.extend(right_points.into_iter());

            points.push(left_points);
        }

        {
            let s = consts.s();
            let x0 = 2.0 * s;
            let pxu0 = pxu_real
                .shift_xp(C::new(pxu_real.xp.re, 0.0), consts)
                .unwrap();
            let pxu1 = pxu0.shift_xp(C::new(x0, 0.0), consts).unwrap();
            let mut right_points = vec![pxu1.clone()];

            for j in 0..64 {
                let x = C::new(x0 + j as f64 / 16.0, 0.0);
                let pxu = right_points.last().unwrap().shift_xp(x, consts);
                if pxu.is_none() {
                    break;
                }
                right_points.push(pxu.unwrap());
            }

            let mut left_points = vec![pxu1.clone()];

            for j in 0..=32 {
                let x = C::new(x0 + (s - x0) * j as f64 / 32.0, 0.0);
                let pxu = left_points.last().unwrap().shift_xp(x, consts);
                if pxu.is_none() {
                    break;
                }
                left_points.push(pxu.unwrap());
            }

            left_points.reverse();
            left_points.pop();

            left_points.extend(right_points.into_iter());

            points.push(left_points);
        }

        {
            let s = consts.s();
            let x0 = 2.0 * s;
            let pxu0 = pxu_real
                .shift_xm(C::new(pxu_real.xm.re, 0.0), consts)
                .unwrap();
            let pxu1 = pxu0.shift_xm(C::new(x0, 0.0), consts).unwrap();
            let mut right_points = vec![pxu1.clone()];

            for j in 0..64 {
                let x = C::new(x0 + j as f64 / 16.0, 0.0);
                let pxu = right_points.last().unwrap().shift_xm(x, consts);
                if pxu.is_none() {
                    break;
                }
                right_points.push(pxu.unwrap());
            }

            let mut left_points = vec![pxu1.clone()];

            for j in 0..=32 {
                let x = C::new(x0 + (s - x0) * j as f64 / 32.0, 0.0);
                let pxu = left_points.last().unwrap().shift_xm(x, consts);
                if pxu.is_none() {
                    break;
                }
                left_points.push(pxu.unwrap());
            }

            left_points.reverse();
            left_points.pop();

            left_points.extend(right_points.into_iter());

            points.push(left_points);
        }

        for i in 1..=0 {
            let p = C::new(i as f64 * PI / 32.0, 0.0);
            let mut up_points = vec![PxuPoint::new(p, consts)];

            for j in 1..=64 {
                let m = 1.0 + j as f64 / 32.0;
                let pxu = up_points.last().unwrap().shift_xp(xp(p, m, consts), consts);
                if pxu.is_none() {
                    break;
                }
                up_points.push(pxu.unwrap());
            }

            let mut down_points = vec![PxuPoint::new(p, consts)];

            for j in 1..=16 {
                let m = 1.0 - j as f64 / 16.0;
                let pxu = down_points
                    .last()
                    .unwrap()
                    .shift_xp(xp(p, m, consts), consts);
                if pxu.is_none() {
                    break;
                }
                down_points.push(pxu.unwrap());
            }

            down_points.reverse();
            down_points.pop();
            down_points.extend(up_points.into_iter());

            points.push(down_points);
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
        let u = u(xp, consts) - C::i() / consts.h;
        Self { p, xp, xm, u }
    }

    pub fn shift_xp(&self, new_xp: C, consts: CouplingConstants) -> Option<Self> {
        let p = find_root(
            |p| xp(p, 1.0, consts) - new_xp,
            |p| dxp_dp(p, 1.0, consts),
            self.p,
            1.0e-6,
            50,
        );
        if p.is_none() {
            return None;
        }
        Some(Self::new(p.unwrap(), consts))
    }

    pub fn shift_xm(&self, new_xm: C, consts: CouplingConstants) -> Option<Self> {
        let p = find_root(
            |p| xm(p, 1.0, consts) - new_xm,
            |p| dxm_dp(p, 1.0, consts),
            self.p,
            1.0e-6,
            50,
        );
        if p.is_none() {
            return None;
        }
        Some(Self::new(p.unwrap(), consts))
    }
}
