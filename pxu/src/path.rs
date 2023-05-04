use itertools::Itertools;
use num::complex::Complex64;

use crate::kinematics::SheetData;
use crate::Component;
use crate::Contours;
use crate::CouplingConstants;
use crate::State;

#[derive(Clone, serde::Deserialize, serde::Serialize)]
pub struct Path {
    pub segments: Vec<Vec<Segment>>,
    pub base_path: BasePath,
}

#[derive(Clone, serde::Deserialize, serde::Serialize)]
pub struct Segment {
    pub p: Vec<Complex64>,
    pub xp: Vec<Complex64>,
    pub xm: Vec<Complex64>,
    pub u: Vec<Complex64>,
    pub sheet_data: SheetData,
}

struct ConstructedSegment {
    start: Complex64,
    end: Complex64,
    path: Vec<(f64, State)>,
}

#[derive(Clone, serde::Deserialize, serde::Serialize)]
pub struct BasePath {
    pub start: State,
    pub path: Vec<Complex64>,
    pub component: Component,
    pub excitation: usize,
    pub name: String,
}

#[derive(Clone, serde::Deserialize, serde::Serialize)]
pub struct SavedPath {
    pub start: State,
    pub deltas: Vec<[i32; 2]>,
    pub component: Component,
    pub excitation: usize,
    pub consts: crate::CouplingConstants,
    pub name: String,
}

const SCALE_FACTOR: f64 = 100_000.0;

impl From<SavedPath> for BasePath {
    fn from(saved_path: SavedPath) -> Self {
        let SavedPath {
            start,
            deltas,
            component,
            excitation,
            consts: _consts,
            name,
        } = saved_path;

        let mut z = start.points[excitation].get(component);
        let mut path = vec![z];
        for dz in deltas {
            z += Complex64::new(dz[0] as f64 / SCALE_FACTOR, dz[1] as f64 / SCALE_FACTOR);
            path.push(z);
        }

        BasePath {
            start,
            path,
            component,
            excitation,
            name,
        }
    }
}

impl From<(BasePath, CouplingConstants)> for SavedPath {
    fn from((base_path, consts): (BasePath, CouplingConstants)) -> Self {
        let BasePath {
            path,
            start,
            component,
            excitation,
            name,
        } = base_path;
        let deltas = path
            .iter()
            .tuple_windows()
            .map(|(a, b)| (b - a))
            .map(|z| {
                [
                    (z.re * SCALE_FACTOR).round() as i32,
                    (z.im * SCALE_FACTOR).round() as i32,
                ]
            })
            .collect();

        SavedPath {
            start,
            deltas,
            component,
            excitation,
            consts,
            name,
        }
    }
}

impl SavedPath {
    pub fn new(
        name: impl Into<String>,
        path: Vec<Complex64>,
        start: State,
        component: Component,
        excitation: usize,
        consts: CouplingConstants,
    ) -> Self {
        let deltas = path
            .into_iter()
            .tuple_windows()
            .map(|(a, b)| b - a)
            .map(|z| {
                [
                    (z.re * 100000.0).round() as i32,
                    (z.im * 100000.0).round() as i32,
                ]
            })
            .collect::<Vec<_>>();

        let name = name.into();

        SavedPath {
            start,
            deltas,
            component,
            excitation,
            consts,
            name,
        }
    }
    pub fn encode(&self) -> Option<String> {
        ron::to_string(&self).ok()
    }

    pub fn encode_compressed(&self) -> Option<String> {
        use base64::Engine;
        use std::io::Write;

        let str = self.encode()?;
        let mut enc = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::best());
        enc.write_all(str.as_bytes()).ok()?;
        let data = enc.finish().ok()?;
        Some(base64::engine::general_purpose::URL_SAFE.encode(data))
    }

    pub fn decode(input: &str) -> Option<Self> {
        use base64::Engine;
        use std::io::Write;

        let input = input.trim();

        if let Ok(path) = ron::from_str(input) {
            return Some(path);
        }
        log::info!("Could not decode RON, trying JSON");
        if let Ok(path) = serde_json::from_str(input) {
            return Some(path);
        }
        log::info!("Could not decode JSON, trying base64");

        let Ok(data) = base64::engine::general_purpose::URL_SAFE.decode(input) else {
            log::warn!("Could not decode base64");
            return None;
        };

        let mut dec = flate2::write::DeflateDecoder::new(Vec::new());
        let Ok(()) = dec.write_all(&data[..]) else {
            log::warn!("Could not deflate");
            return None;
        };
        let Ok(data) = dec.finish() else {
            log::warn!("Could not deflate");
            return None;
        };
        let Ok(input) = String::from_utf8(data) else {
            log::warn!("Resulting data is not a string");
            return None;
        };
        if let Ok(saved_path) = ron::from_str::<SavedPath>(&input) {
            return Some(saved_path);
        }
        log::warn!("Could not decode RON");
        None
    }

    pub fn save(paths: &Vec<Self>) -> Option<String> {
        ron::to_string(paths).ok()
    }

    pub fn save_compressed(paths: &Vec<Self>) -> Option<String> {
        use base64::Engine;
        use std::io::Write;

        let str = Self::save(paths)?;
        let mut enc = flate2::write::DeflateEncoder::new(Vec::new(), flate2::Compression::best());
        enc.write_all(str.as_bytes()).ok()?;
        let data = enc.finish().ok()?;
        Some(base64::engine::general_purpose::URL_SAFE.encode(data))
    }

    pub fn load(input: &str) -> Option<Vec<Self>> {
        use base64::Engine;
        use std::io::Write;

        let input = input.trim();

        if let Ok(saved_paths) = ron::from_str(input) {
            return Some(saved_paths);
        }
        log::info!("Could not decode RON, trying JSON");
        if let Ok(saved_paths) = serde_json::from_str(input) {
            return Some(saved_paths);
        }
        log::info!("Could not decode JSON, trying base64");

        let Ok(data) = base64::engine::general_purpose::URL_SAFE.decode(input) else {
            log::warn!("Could not decode base64");
            return None;
        };

        let mut dec = flate2::write::DeflateDecoder::new(Vec::new());
        let Ok(()) = dec.write_all(&data[..]) else {
            log::warn!("Could not deflate");
            return None;
        };
        let Ok(data) = dec.finish() else {
            log::warn!("Could not deflate");
            return None;
        };
        let Ok(input) = String::from_utf8(data) else {
            log::warn!("Resulting data is not a string");
            return None;
        };
        if let Ok(saved_paths) = ron::from_str(&input) {
            return Some(saved_paths);
        }
        log::warn!("Could not decode RON");
        None
    }
}

impl ConstructedSegment {
    fn segments(self) -> Vec<Segment> {
        let mut segments = vec![];
        for i in 0..self.path[0].1.points.len() {
            let mut p = vec![];
            let mut xp = vec![];
            let mut xm = vec![];
            let mut u = vec![];

            let sheet_data = self.path[0].1.points[i].sheet_data.clone();

            for (_, state) in self.path.iter() {
                if state.points[i].sheet_data != sheet_data {
                    log::warn!(
                        "Inconsistent sheet data ({i}): {:?} vs {:?}",
                        sheet_data,
                        state.points[i].sheet_data
                    );
                }
                p.push(state.points[i].p);
                xp.push(state.points[i].xp);
                xm.push(state.points[i].xm);
                u.push(state.points[i].u);
            }

            segments.push(Segment {
                p,
                xp,
                xm,
                u,
                sheet_data,
            });
        }

        segments
    }

    fn normalize(&mut self) {
        let Some(&(t0,_)) = self.path.first() else {return};
        let Some(&(t1,_)) = self.path.last() else {return};

        for entry in self.path.iter_mut() {
            entry.0 = (entry.0 - t0) / (t1 - t0);
        }
    }
    fn refine(&mut self, base_path: &BasePath, contours: &Contours, consts: CouplingConstants) {
        let m = self.path[0].1.points.len();

        let min_cos = (2.0 * std::f64::consts::TAU / 360.0).cos();

        for _ in 0..10 {
            let mut refinements: Vec<(f64, State)> = vec![];

            let mut prev = false;
            for ((t1, s1), (t2, s2), (t3, s3)) in self.path.iter().tuple_windows::<(_, _, _)>() {
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
                let z = self.start * (1.0 - t) + self.end * t;
                let mut state = state;
                state.update(
                    base_path.excitation,
                    base_path.component,
                    z,
                    contours,
                    consts,
                );

                self.path.push((t, state));
            }

            self.path.sort_unstable_by(|(t1, _), (t2, _)| {
                t1.partial_cmp(t2).unwrap_or(std::cmp::Ordering::Greater)
            });
        }
    }

    fn refine_all(
        mut segments: Vec<Self>,
        base_path: &BasePath,
        contours: &Contours,
        consts: CouplingConstants,
    ) -> Vec<Self> {
        segments
            .iter_mut()
            .for_each(|segment| segment.refine(base_path, contours, consts));
        segments
    }

    fn split(
        self,
        base_path: &BasePath,
        contours: &Contours,
        consts: CouplingConstants,
    ) -> Vec<Self> {
        let mut segments: Vec<ConstructedSegment> = vec![ConstructedSegment {
            start: self.start,
            end: self.end,
            path: vec![self.path[0].clone()],
        }];

        for ((t1, s1), (t2, s2)) in self.path.into_iter().tuple_windows::<(_, _)>() {
            let mut state = s1.clone();

            let z1 = s1.points[base_path.excitation].get(base_path.component);
            let z2 = s2.points[base_path.excitation].get(base_path.component);

            loop {
                let s = std::iter::zip(state.points.iter(), s2.points.iter())
                    .filter_map(|(pt1, pt2)| {
                        contours
                            .get_crossed_cuts(
                                pt1,
                                base_path.component,
                                pt2.get(base_path.component),
                                consts,
                            )
                            .first()
                            .map(|&(t, _)| t)
                    })
                    .fold(1.0_f64, |a, b| a.min(b));

                if s == 1.0 {
                    segments.last_mut().unwrap().path.push((t2, s2));
                    break;
                }

                let zs = state.points[base_path.excitation].get(base_path.component);

                let z = zs + (s - 0.01).clamp(0.0, 1.0) * (z2 - zs);

                if !state.update(
                    base_path.excitation,
                    base_path.component,
                    z,
                    contours,
                    consts,
                ) {
                    log::warn!(
                        "Couldn't update #1 ({} {:?})",
                        base_path.excitation,
                        base_path.component
                    );
                    break;
                }

                let t = t1 + ((z - z1) / (z2 - z1)).re * (t2 - t1);

                segments.last_mut().unwrap().path.push((t, state.clone()));
                segments.last_mut().unwrap().end =
                    state.points[base_path.excitation].get(base_path.component);

                let z = zs + (s + 0.01).clamp(0.0, 1.0) * (z2 - zs);

                if !state.update(
                    base_path.excitation,
                    base_path.component,
                    z,
                    contours,
                    consts,
                ) {
                    log::warn!(
                        "Couldn't update #2 ({} {:?})",
                        base_path.excitation,
                        base_path.component,
                    );
                    break;
                }

                let t = t1 + ((z - z1) / (z2 - z1)).re * (t2 - t1);

                segments.push(ConstructedSegment {
                    start: state.points[base_path.excitation].get(base_path.component),
                    end: self.end,
                    path: vec![(t, state.clone())],
                });
            }
        }

        segments.iter_mut().for_each(|segment| segment.normalize());

        segments
    }

    fn split_all(
        segments: Vec<Self>,
        base_path: &BasePath,
        contours: &Contours,
        consts: CouplingConstants,
    ) -> Vec<Self> {
        let mut result = vec![];
        for segment in segments {
            result.extend(segment.split(base_path, contours, consts));
        }
        result
    }
}

impl Path {
    pub fn from_base_path(
        base_path: BasePath,
        contours: &Contours,
        consts: CouplingConstants,
    ) -> Self {
        let mut state = base_path.start.clone();

        let mut segments = vec![];

        let max_step = match base_path.component {
            Component::P => 0.05,
            Component::Xp | Component::Xm => 0.1,
            Component::U => 0.5 / consts.h,
        };

        for (start, end) in base_path.path.iter().tuple_windows() {
            let mut path = vec![];
            let steps = ((end - start).norm() / max_step).ceil() as usize;

            for step in 0..=steps {
                let t = step as f64 / steps as f64;
                let z = start * (1.0 - t) + end * t;
                state.update(
                    base_path.excitation,
                    base_path.component,
                    z,
                    contours,
                    consts,
                );

                path.push((t, state.clone()));
            }

            segments.push(ConstructedSegment {
                start: *start,
                end: *end,
                path,
            })
        }

        segments = ConstructedSegment::split_all(segments, &base_path, contours, consts);
        segments = ConstructedSegment::refine_all(segments, &base_path, contours, consts);

        let segments = segments
            .into_iter()
            .map(|segment| segment.segments())
            .collect::<Vec<_>>();

        let rows = segments.len();
        let cols = segments[0].len();

        let mut segments: Vec<Vec<_>> = (0..cols)
            .map(|col| (0..rows).map(|row| segments[row][col].clone()).collect())
            .collect();

        segments
            .iter_mut()
            .for_each(|excitation| excitation.iter_mut().for_each(|segment| segment.simplify()));

        log::info!(
            "{},{},{},{} points",
            segments[0].iter().map(|s| s.p.len()).sum::<usize>(),
            segments[0].iter().map(|s| s.xp.len()).sum::<usize>(),
            segments[0].iter().map(|s| s.xm.len()).sum::<usize>(),
            segments[0].iter().map(|s| s.u.len()).sum::<usize>(),
        );

        Self {
            base_path,
            segments,
        }
    }
}

impl Segment {
    fn simplify_line(points: &mut Vec<Complex64>) {
        if points.len() < 2 {
            return;
        }
        let mut new_points = vec![points[0]];
        for (z1, z2, z3) in points.iter().tuple_windows::<(_, _, _)>() {
            let w1 = z2 - z1;
            let w2 = z3 - z2;
            let cross = w1.re * w2.im - w1.im * w2.re;

            if cross * cross > w1.norm_sqr() * w2.norm_sqr() * 1.0e-6 {
                new_points.push(*z2);
            }
        }
        new_points.push(*points.last().unwrap());
        *points = new_points;
    }

    fn simplify(&mut self) {
        Self::simplify_line(&mut self.p);
        Self::simplify_line(&mut self.xp);
        Self::simplify_line(&mut self.xm);
        Self::simplify_line(&mut self.u);
    }

    pub fn get(&self, component: Component) -> &Vec<Complex64> {
        match component {
            Component::P => &self.p,
            Component::Xp => &self.xp,
            Component::Xm => &self.xm,
            Component::U => &self.u,
        }
    }
}

impl BasePath {
    pub fn from_editable_path(
        editable_path: &EditablePath,
        component: Component,
        excitation: usize,
    ) -> Self {
        let start = editable_path.states.first().unwrap().clone();
        let path = editable_path
            .states
            .iter()
            .map(|state| state.points[excitation].get(component))
            .collect();

        Self {
            start,
            path,
            component,
            excitation,
            name: "(Unnamed)".to_owned(),
        }
    }
}

#[derive(Default, serde::Deserialize, serde::Serialize)]
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
