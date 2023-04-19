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

        for z in base_path.path.iter() {
            state.update(
                base_path.excitation,
                base_path.component,
                *z,
                contours,
                consts,
            );

            for (i, point) in state.points.iter().enumerate() {
                p[i].push(point.get(Component::P));
                xp[i].push(point.get(Component::Xp));
                xm[i].push(point.get(Component::Xm));
                u[i].push(point.get(Component::U));
            }
        }

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
