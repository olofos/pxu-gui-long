use crate::contours::{Component, Contours};
use crate::interpolation::PInterpolatorMut;
use crate::kinematics::{xm, xm_crossed, xp, xp_crossed, CouplingConstants};
use crate::point::Point;
use num::complex::Complex64;

#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct State {
    pub points: Vec<Point>,
}

impl State {
    pub fn new(m: usize, consts: CouplingConstants) -> Self {
        let mut points = vec![];

        let mut p_int = PInterpolatorMut::xp(0.025, consts);
        p_int
            .goto_m(m as f64)
            .goto_p(0.025 + 0.022 * (m - 1) as f64);
        let mut pt = Point::new(p_int.p(), consts);

        let s = consts.s();
        let us = s + 1.0 / s - (s - 1.0 / s) * s.ln();

        let u0 = us + 3.0;
        let step_size = 1.0 / 4.0;
        let max_steps = 2 * ((u0 - pt.u.re).abs() / step_size) as usize;
        for _ in 0..max_steps {
            let du = u0 - pt.u.re;
            let u = pt.u.re + du.abs().min(step_size).copysign(du);
            pt.update(Component::U, Complex64::new(u, pt.u.im), &[], consts);
            if (u0 - pt.u.re).abs() < 0.01 {
                break;
            }
        }
        if (u0 - pt.u.re).abs() >= 0.01 {
            log::warn!(
                "Could not find u (h={} k={} du={})",
                consts.h,
                consts.k(),
                u0 - pt.u.re
            );
        }
        points.push(pt);

        for i in 1..m {
            let mut pt = points[i - 1].clone();
            let xm = pt.xm;
            let steps = 4;
            for _ in 1..=steps {
                pt.update(Component::Xp, xm, &[], consts);
            }
            points.push(pt);
        }

        Self { points }
    }

    pub fn update_points(
        points: &mut Vec<Point>,
        active_point: usize,
        component: Component,
        new_value: Complex64,
        contours: &Contours,
        consts: CouplingConstants,
    ) {
        let crossed_cuts = contours
            .get_crossed_cuts(&points[active_point], component, new_value, consts)
            .collect::<Vec<_>>();

        points[active_point].update(component, new_value, &crossed_cuts, consts);

        for i in (active_point + 1)..points.len() {
            let new_value = if points[i - 1].sheet_data.e_branch > 0 {
                xm(points[i - 1].p, 1.0, consts)
            } else {
                xm_crossed(points[i - 1].p, 1.0, consts)
            };
            let crossed_cuts = contours
                .get_crossed_cuts(&points[i], Component::Xp, new_value, consts)
                .collect::<Vec<_>>();
            points[i].update(Component::Xp, new_value, &crossed_cuts, consts);
        }

        for i in (0..active_point).rev() {
            let new_value = if points[i + 1].sheet_data.e_branch > 0 {
                xp(points[i + 1].p, 1.0, consts)
            } else {
                xp_crossed(points[i + 1].p, 1.0, consts)
            };
            let crossed_cuts = contours
                .get_crossed_cuts(&points[i], Component::Xm, new_value, consts)
                .collect::<Vec<_>>();
            points[i].update(Component::Xm, new_value, &crossed_cuts, consts);
        }
    }

    pub fn update(
        &mut self,
        active_point: usize,
        component: Component,
        new_value: Complex64,
        contours: &Contours,
        consts: CouplingConstants,
    ) {
        Self::update_points(
            &mut self.points,
            active_point,
            component,
            new_value,
            contours,
            consts,
        );
    }

    pub fn p(&self) -> Complex64 {
        self.points.iter().map(|pxu| pxu.p).sum::<Complex64>()
    }

    pub fn en(&self, consts: CouplingConstants) -> Complex64 {
        self.points
            .iter()
            .map(|pt| pt.en(consts))
            .sum::<Complex64>()
    }
}
