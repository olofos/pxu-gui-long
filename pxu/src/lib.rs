#![warn(clippy::all, rust_2018_idioms)]

mod contours;
mod cut;
pub mod interpolation;
pub mod kinematics;
mod nr;
pub mod path;
mod point;
mod state;

pub use contours::{
    compute_branch_point, BranchPointType, Component, Contours, GridLine, GridLineComponent,
    UCutType,
};
pub use cut::{Cut, CutType};
pub use kinematics::CouplingConstants;
pub use point::Point;
pub use state::State;

#[derive(serde::Deserialize, serde::Serialize)]
pub struct Pxu {
    pub consts: CouplingConstants,
    #[serde(skip)]
    pub contours: Contours,
    pub state: State,
    // pub active_point: usize,
}

impl Pxu {
    pub fn new(consts: CouplingConstants) -> Self {
        Self {
            consts,
            contours: Default::default(),
            state: Default::default(),
            // active_point: Default::default(),
        }
    }

    // pub fn active_point(&self) -> &Point {
    //     &self.state.points[self.active_point]
    // }

    // pub fn set_active_point(&mut self, active_point: usize) {
    //     if active_point <= self.state.points.len() {
    //         self.active_point = active_point;
    //     }
    // }
}
