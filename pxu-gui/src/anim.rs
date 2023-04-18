#![allow(dead_code)]

use itertools::Itertools;
use num::complex::Complex64;

#[derive(Debug)]
struct PathSegment {
    t_start: f64,
    t_end: f64,
    z_start: Complex64,
    z_end: Complex64,
}

#[derive(Debug, Default, PartialEq, Eq)]
pub enum Status {
    #[default]
    Stopped,
    Playing,
    Paused,
}

pub struct Anim {
    pub status: Status,
    segment: usize,
    prev_frame: chrono::DateTime<chrono::Utc>,
    path: Vec<PathSegment>,
    z: Complex64,
    pub active_point: usize,
    pub speed: f64,
    pub t: f64,
    pub component: pxu::Component,
    pub total_len: f64,
    pub start: Vec<pxu::Point>,
    pub end: Vec<pxu::Point>,
}

impl Default for Anim {
    fn default() -> Self {
        Self {
            status: Default::default(),
            segment: Default::default(),
            prev_frame: Default::default(),
            path: Default::default(),
            z: Default::default(),
            active_point: Default::default(),
            speed: 50.0,
            t: Default::default(),
            component: pxu::Component::Xp,
            total_len: Default::default(),
            start: Default::default(),
            end: Default::default(),
        }
    }
}

impl Anim {
    pub fn update(&mut self) -> Option<Complex64> {
        if self.status == Status::Stopped {
            return None;
        }

        let now = chrono::Utc::now();
        if self.status == Status::Playing {
            let speed = self.speed
                * match self.component {
                    pxu::Component::P => 0.01,
                    pxu::Component::Xp | pxu::Component::Xm | pxu::Component::U => 0.1,
                };
            self.t +=
                (now - self.prev_frame).num_nanoseconds().unwrap() as f64 / 1_000_000_000.0 * speed;
        }
        if self.t >= self.total_len {
            self.t = self.total_len;
            self.status = Status::Paused;
        }
        self.prev_frame = now;

        // let mut segment = 0;
        // while self.path[segment].t_end < self.t {
        //     segment += 1;
        //     if segment >= self.path.len() {
        //         self.status = Status::Paused;
        //         return Some(self.path.last().unwrap().z_end);
        //     }
        // }

        if self.path[self.segment].t_end < self.t {
            self.segment += 1;
        }
        if self.path[self.segment].t_start > self.t {
            self.segment -= 1;
        }
        if self.segment >= self.path.len() {
            self.status = Status::Paused;
            return Some(self.path.last().unwrap().z_end);
        }

        let seg = &self.path[self.segment];

        let dt = self.t - seg.t_start;
        let dp = (seg.z_end - seg.z_start) / (seg.z_end - seg.z_start).norm();
        Some(seg.z_start + dt * dp)
    }

    pub fn start(&mut self, paths: &[Vec<pxu::Point>]) -> Vec<pxu::Point> {
        self.active_point = paths
            .iter()
            .map(|path| {
                path.iter()
                    .map(|pt| pt.get(self.component))
                    .tuple_windows::<(_, _)>()
                    .map(|(p1, p2)| (p2 - p1).norm())
                    .sum::<f64>()
            })
            .enumerate()
            .max_by(|(_, x), (_, y)| x.partial_cmp(y).unwrap())
            .map(|(n, _)| n)
            .unwrap();

        let p = paths[self.active_point][0].get(self.component);
        self.z = p;
        self.path = vec![];

        let ps = paths[self.active_point]
            .iter()
            .map(|pt| pt.get(self.component))
            .collect::<Vec<_>>();
        let mut t = 0.0;
        for i in 0..ps.len() - 1 {
            let p_start = ps[i];
            let p_end = ps[i + 1];
            let dp = p_end - p_start;
            let dt = dp.norm();

            let t_start = t;
            let t_end = t + dt;

            self.path.push(PathSegment {
                z_start: p_start,
                z_end: p_end,
                t_start,
                t_end,
            });

            t += dt;
        }

        self.status = Status::Playing;
        self.prev_frame = chrono::Utc::now();
        self.t = 0.0;
        self.segment = 0;
        self.total_len = t;

        self.start = paths
            .iter()
            .map(|path| path.first().unwrap().clone())
            .collect();
        self.end = paths
            .iter()
            .map(|path| path.last().unwrap().clone())
            .collect();

        self.start.clone()
    }

    pub fn stop(&mut self) {
        self.status = Status::Stopped;
    }

    pub fn pause(&mut self) {
        self.status = Status::Paused;
    }

    pub fn unpause(&mut self) {
        self.status = Status::Playing;
    }

    pub fn _is_playing(&self) -> bool {
        self.status == Status::Playing
    }

    pub fn is_stopped(&self) -> bool {
        self.status == Status::Stopped
    }

    pub fn is_paused(&self) -> bool {
        self.status == Status::Paused
    }

    pub fn goto_start(&mut self) -> Vec<pxu::Point> {
        self.t = 0.0;
        self.segment = 0;
        self.start.clone()
    }

    pub fn goto_end(&mut self) -> Vec<pxu::Point> {
        self.t = self.total_len;
        self.segment = self.path.len() - 1;
        self.end.clone()
    }
}
