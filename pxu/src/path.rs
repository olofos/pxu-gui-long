use num::complex::Complex64;

use crate::kinematics::SheetData;
use crate::Component;
use crate::Point;

pub struct PathSegment {
    pub p: Vec<Vec<Complex64>>,
    pub xp: Vec<Vec<Complex64>>,
    pub xm: Vec<Vec<Complex64>>,
    pub u: Vec<Vec<Complex64>>,
    pub sheet_data: SheetData,
}

#[derive(Default)]
pub struct Path {
    pub segments: Vec<PathSegment>,
}

#[derive(Default)]
pub struct EditablePath {
    pub paths: Vec<Vec<Point>>,
}

impl EditablePath {
    pub fn clear(&mut self) {
        self.paths = vec![];
    }

    pub fn get(&self, component: Component) -> Vec<Vec<Complex64>> {
        if self.paths.is_empty() || self.paths[0].is_empty() {
            return vec![];
        }

        let mut result = vec![vec![]; self.paths[0].len()];

        for points in self.paths.iter() {
            for (i, point) in points.iter().enumerate() {
                result[i].push(point.get(component));
            }
        }

        // self.points
        //     .iter()
        //     .map(|pts| pts.iter().map(|pt| pt.get(component)).collect())
        //     .collect()

        result
    }

    pub fn push(&mut self, points: &Vec<Point>) {
        if !self.paths.is_empty() && points.len() != self.paths[0].len() {
            self.clear();
        }

        self.paths.push(points.clone());
    }
}
