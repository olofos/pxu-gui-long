#![warn(clippy::all, rust_2018_idioms)]

mod contours;
mod cut;
pub mod interpolation;
pub mod kinematics;
mod nr;
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
}

impl Pxu {
    pub fn new(consts: CouplingConstants) -> Self {
        Self {
            consts,
            contours: Contours::default(),
            state: State::default(),
        }
    }
}
