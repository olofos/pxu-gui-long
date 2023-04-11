use crate::contours::{Component, Contours, UCutType};
use crate::interpolation::PInterpolatorMut;
use crate::kinematics::{xm, xm_crossed, xp, xp_crossed, CouplingConstants};
use crate::point::Point;
use num::complex::Complex64;

#[derive(Debug, Clone)]
pub struct State {
    pub points: Vec<Point>,
    pub consts: CouplingConstants,
    pub active_point: usize,
}

impl State {
    pub fn new(m: usize, consts: CouplingConstants) -> Self {
        let mut points = vec![];

        let mut p_int = PInterpolatorMut::xp(0.05, consts);
        p_int.goto_m(m as f64).goto_p(0.1);
        let mut pt = Point::new(p_int.p(), consts);

        let s = consts.s();
        let us = s + 1.0 / s - (s - 1.0 / s) * s.ln();

        let u0 = us + 3.0;
        let step_size = 0.25;
        let max_steps = 2 * ((u0 - pt.u.re).abs() / 0.25) as usize;
        for _ in 0..max_steps {
            let du = u0 - pt.u.re;
            let u = pt.u.re + du.abs().min(step_size).copysign(du);
            pt.update(Component::U, Complex64::new(u, pt.u.im), &[]);
            if (u0 - pt.u.re).abs() < 0.01 {
                break;
            }
        }
        log::info!("{} {max_steps}", u0 - pt.u.re);
        points.push(pt);

        for i in 1..m {
            let mut pt = points[i - 1].clone();
            let u = pt.u;
            let steps = 4;
            for n in 1..=steps {
                pt.update(
                    Component::U,
                    u - 2.0 * Complex64::i() / consts.h * (n as f64 / steps as f64),
                    &[],
                );
            }
            points.push(pt);
        }

        let active_point = points.len() / 2;

        Self {
            points,
            consts,
            active_point,
        }
    }

    pub fn active_point(&self) -> &Point {
        &self.points[self.active_point]
    }

    pub fn update(
        &mut self,
        active_point: usize,
        component: Component,
        new_value: Complex64,
        contour_generator: &Contours,
        u_cut_type: UCutType,
    ) {
        let crossed_cuts = contour_generator
            .get_crossed_cuts(&self.points[active_point], component, new_value, u_cut_type)
            .collect::<Vec<_>>();

        self.active_point = active_point;
        self.points[self.active_point].update(component, new_value, &crossed_cuts);

        for i in (self.active_point + 1)..self.points.len() {
            let new_value = if self.points[i - 1].sheet_data.e_branch > 0 {
                xm(self.points[i - 1].p, 1.0, self.points[i - 1].consts)
            } else {
                xm_crossed(self.points[i - 1].p, 1.0, self.points[i - 1].consts)
            };
            let crossed_cuts = contour_generator
                .get_crossed_cuts(&self.points[i], Component::Xp, new_value, u_cut_type)
                .collect::<Vec<_>>();
            self.points[i].update(Component::Xp, new_value, &crossed_cuts);
        }

        for i in (0..self.active_point).rev() {
            let new_value = if self.points[i + 1].sheet_data.e_branch > 0 {
                xp(self.points[i + 1].p, 1.0, self.points[i + 1].consts)
            } else {
                xp_crossed(self.points[i + 1].p, 1.0, self.points[i + 1].consts)
            };
            let crossed_cuts = contour_generator
                .get_crossed_cuts(&self.points[i], Component::Xm, new_value, u_cut_type)
                .collect::<Vec<_>>();
            self.points[i].update(Component::Xm, new_value, &crossed_cuts);
        }
    }
}
