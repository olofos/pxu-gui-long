use itertools::Itertools;
use num::complex::Complex64;

use crate::Component;
use crate::Contours;
use crate::CouplingConstants;
use crate::State;

// pub struct PathSegment {
//     pub p: Vec<Vec<Complex64>>,
//     pub xp: Vec<Vec<Complex64>>,
//     pub xm: Vec<Vec<Complex64>>,
//     pub u: Vec<Vec<Complex64>>,
//     pub sheet_data: SheetData,
// }

#[derive(Default, serde::Deserialize, serde::Serialize)]
#[serde(default)]
pub struct Path {
    pub p: Vec<Vec<Complex64>>,
    pub xp: Vec<Vec<Complex64>>,
    pub xm: Vec<Vec<Complex64>>,
    pub u: Vec<Vec<Complex64>>,
    pub base_path: Option<BasePath>,
}

impl Path {
    pub fn from_base_path(
        base_path: BasePath,
        contours: &Contours,
        consts: CouplingConstants,
    ) -> Self {
        let m = base_path.start.points.len();
        let mut state = base_path.start.clone();

        let mut p = vec![vec![]; m];
        let mut xp = vec![vec![]; m];
        let mut xm = vec![vec![]; m];
        let mut u = vec![vec![]; m];

        for (i, point) in state.points.iter().enumerate() {
            p[i].push(point.get(Component::P));
            xp[i].push(point.get(Component::Xp));
            xm[i].push(point.get(Component::Xm));
            u[i].push(point.get(Component::U));
        }

        let max_step = match base_path.component {
            Component::P => 0.05,
            Component::Xp | Component::Xm => 0.1,
            Component::U => 0.5 / consts.h,
        };

        for (segment_start, segment_end) in base_path.path.iter().tuple_windows() {
            let steps = ((segment_end - segment_start).norm() / max_step).ceil() as usize;

            let mut segment = vec![];
            for step in 0..=steps {
                let t = step as f64 / steps as f64;
                let z = segment_start * (1.0 - t) + segment_end * t;
                state.update(
                    base_path.excitation,
                    base_path.component,
                    z,
                    contours,
                    consts,
                );

                segment.push((t, state.clone()));
            }

            let min_cos = (2.0 * std::f64::consts::TAU / 360.0).cos();

            for _ in 0..5 {
                let mut refinements: Vec<(f64, State)> = vec![];

                let mut prev = false;
                for ((t1, s1), (t2, s2), (t3, s3)) in segment.iter().tuple_windows::<(_, _, _)>() {
                    let mut refine = false;
                    'comp: for comp in [Component::P, Component::Xp, Component::Xm, Component::U] {
                        for i in 0..m {
                            let x1 = s1.points[i].get(comp);
                            let x2 = s2.points[i].get(comp);
                            let x3 = s3.points[i].get(comp);

                            let z1 = x2 - x1;
                            let z2 = x3 - x2;
                            let cos = (z1.re * z2.re + z1.im * z2.im) / (z1.norm() * z2.norm());
                            if cos < min_cos {
                                refine = true;
                                break 'comp;
                            }
                        }
                    }

                    if refine {
                        if !prev {
                            refinements.push(((t1 + t2) / 2.0, s1.clone()));
                        }
                        refinements.push(((t2 + t3) / 2.0, s2.clone()));
                        prev = true;
                    } else {
                        prev = false;
                    }
                }

                if refinements.is_empty() {
                    break;
                }

                for (t, state) in refinements.into_iter() {
                    let z = segment_start * (1.0 - t) + segment_end * t;
                    let mut state = state;
                    state.update(
                        base_path.excitation,
                        base_path.component,
                        z,
                        contours,
                        consts,
                    );

                    segment.push((t, state));
                }

                segment.sort_unstable_by(|(t1, _), (t2, _)| {
                    t1.partial_cmp(t2).unwrap_or(std::cmp::Ordering::Greater)
                });
            }

            let mut extra_states = vec![];
            for ((t1, s1), (t2, s2)) in segment.iter().tuple_windows() {
                for i in 0..s1.points.len() {
                    if s1.points[i].sheet_data != s2.points[i].sheet_data {
                        let crossed_cuts = contours.get_crossed_cuts(
                            &s1.points[i],
                            base_path.component,
                            s2.points[i].get(base_path.component),
                            consts,
                        );
                        let z1 = s1.points[i].get(base_path.component);
                        let z2 = s2.points[i].get(base_path.component);
                        for cut in crossed_cuts {
                            if matches!(cut.typ, crate::CutType::ULongNegative(_)) {
                                continue;
                            }
                            let Some((_, _, u)) = cut.intersection(
                                z1,z2,
                                consts,
                            ) else {
                                continue;
                            };

                            let mut state = s2.clone();
                            state.update(
                                base_path.excitation,
                                base_path.component,
                                z1 * (1.0 - u) + z2 * u,
                                contours,
                                consts,
                            );
                            let t = (1.0 - u) * t1 + u * t2;
                            extra_states.push((t, state));
                        }
                    }
                }
            }

            segment.extend(extra_states);
            segment.sort_unstable_by(|(t1, _), (t2, _)| {
                t1.partial_cmp(t2).unwrap_or(std::cmp::Ordering::Greater)
            });

            for (_, state) in segment.into_iter().skip(1) {
                for (i, point) in state.points.iter().enumerate() {
                    p[i].push(point.get(Component::P));
                    xp[i].push(point.get(Component::Xp));
                    xm[i].push(point.get(Component::Xm));
                    u[i].push(point.get(Component::U));
                }
            }
        }
        log::info!("len = {}", p[0].len());

        Self {
            p,
            xp,
            xm,
            u,
            base_path: Some(base_path),
        }
    }
}

#[derive(serde::Deserialize, serde::Serialize)]
pub struct BasePath {
    pub start: State,
    pub end: State,
    pub path: Vec<Complex64>,
    pub component: Component,
    pub excitation: usize,
}

impl BasePath {
    pub fn from_editable_path(
        editable_path: &EditablePath,
        component: Component,
        excitation: usize,
    ) -> Self {
        let start = editable_path.states.first().unwrap().clone();
        let end = editable_path.states.last().unwrap().clone();
        let path = editable_path
            .states
            .iter()
            .map(|state| state.points[excitation].get(component))
            .collect();

        Self {
            start,
            end,
            path,
            component,
            excitation,
        }
    }
}

#[derive(Default)]
pub struct EditablePath {
    pub states: Vec<State>,
}

impl EditablePath {
    pub fn clear(&mut self) {
        self.states = vec![];
    }

    pub fn get(&self, component: Component) -> Vec<Vec<Complex64>> {
        if self.states.is_empty() {
            return vec![];
        }

        let mut result = vec![vec![]; self.states[0].points.len()];

        for state in self.states.iter() {
            for (i, point) in state.points.iter().enumerate() {
                result[i].push(point.get(component));
            }
        }

        result
    }

    pub fn push(&mut self, state: &State) {
        self.states.push(state.clone());
    }
}
