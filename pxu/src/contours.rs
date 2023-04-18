use std::collections::VecDeque;

use crate::cut::{Cut, CutType, CutVisibilityCondition};
use crate::interpolation::{EPInterpolator, InterpolationPoint, PInterpolatorMut, XInterpolator};
use crate::kinematics::{xp, CouplingConstants, UBranch};
use crate::{nr, Point, State};
use itertools::Itertools;
use num::complex::Complex64;

const P_RANGE_MIN: i32 = -3;
const P_RANGE_MAX: i32 = 3;

const INFINITY: f64 = 100.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Component {
    P,
    Xp,
    Xm,
    U,
}

impl Component {
    pub fn conj(&self) -> Self {
        match self {
            Self::P => Self::P,
            Self::Xp => Self::Xm,
            Self::Xm => Self::Xp,
            Self::U => Self::U,
        }
    }
}

impl UBranch {
    pub fn cross_scallion(self) -> Self {
        match self {
            Self::Outside => Self::Between,
            Self::Between => Self::Outside,
            Self::Inside => {
                log::error!("Can't cross scallion from inside kidney");
                Self::Outside
            }
        }
    }

    pub fn cross_kidney(self) -> Self {
        match self {
            Self::Inside => Self::Between,
            Self::Between => Self::Inside,
            Self::Outside => {
                log::error!("Can't cross kidney from outside scallion");
                Self::Inside
            }
        }
    }

    pub fn cross(self, cut_typ: &CutType) -> Self {
        match cut_typ {
            CutType::UShortScallion(_) => self.cross_scallion(),
            CutType::UShortKidney(_) => self.cross_kidney(),
            _ => self,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
pub enum UCutType {
    Long,
    SemiShort,
    #[default]
    Short,
}

impl std::fmt::Display for UCutType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let name = match self {
            Self::Long => "Long",
            Self::SemiShort => "Semi short",
            Self::Short => "Short",
        };
        write!(f, "{name}")
    }
}

impl UCutType {
    pub fn all() -> impl Iterator<Item = Self> {
        [Self::Long, Self::SemiShort, Self::Short].into_iter()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[allow(clippy::enum_variant_names)]
pub enum BranchPointType {
    XpPositiveAxisImXmNegative,
    XpPositiveAxisImXmPositive,
    XpNegativeAxisFromAboveWithImXmNegative,
    XpNegativeAxisFromBelowWithImXmNegative,
    XpNegativeAxisFromAboveWithImXmPositive,
    XpNegativeAxisFromBelowWithImXmPositive,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BranchPointData {
    pub p: f64,
    pub m: f64,
    pub typ: BranchPointType,
}

#[derive(Debug)]
enum CutDirection {
    Positive,
    Negative,
}

#[derive(Debug)]
enum XCut {
    Scallion,
    Kidney,
}

#[derive(Debug, PartialEq, Eq)]
enum SplitCutBranchPoint {
    Old,
    New,
}

#[derive(Debug, PartialEq, Eq)]
enum SplitCutOrder {
    OldFirst,
    NewFirst,
}

#[derive(Debug)]
enum GeneratorCommand {
    AddGridLineU {
        y: f64,
    },
    AddGridLineX {
        m: f64,
    },
    AddGridLineP,

    ComputeBranchPoint {
        p_range: i32,
        branch_point_type: BranchPointType,
    },
    ClearCut,
    ComputeCutX(CutDirection),
    ComputeCutXFull(XCut),
    ComputeCutP {
        reverse: bool,
    },
    ComputeCutEP,
    ComputeCutEXp,
    ComputeCutEXm,
    ComputeCutEU,
    SetCutPath {
        path: Vec<Complex64>,
        branch_point: Option<Complex64>,
    },
    PushCut {
        p_range: i32,
        component: Component,
        cut_type: CutType,
        periodic: bool,
        visibility: Vec<CutVisibilityCondition>,
        pre_shift: Complex64,
    },
    PopCut,
    SwapCuts,
    SplitCut(Component, CutType, SplitCutBranchPoint, SplitCutOrder),

    EStart {
        p_range: i32,
    },

    PStartXp {
        p: f64,
    },
    PGotoXp(f64, f64),
    PGotoXm(f64, f64),
    PGotoP(f64),
    PGotoM(f64),
    PGotoIm(f64),
    PGotoRe(f64),
}

#[derive(Default)]
struct RuntimeCutData {
    branch_point: Option<Complex64>,
    path: Option<Vec<Complex64>>,
}

#[derive(Default)]
struct ContourGeneratorRuntimeContext {
    p_int: Option<PInterpolatorMut>,
    e_int: Option<EPInterpolator>,
    branch_point_data: Option<BranchPointData>,
    cut_data: RuntimeCutData,
}

struct ContourCommandGenerator {
    component: Option<Component>,
    cut_type: Option<CutType>,
    periodic: bool,
    visibility: Vec<CutVisibilityCondition>,
    group_visibility: Vec<CutVisibilityCondition>,
    pre_shift: Complex64,
    commands: VecDeque<GeneratorCommand>,
}

pub enum GridLineComponent {
    Real,
    Xp(f64),
    Xm(f64),
}

pub struct GridLine {
    pub path: Vec<Complex64>,
    pub component: GridLineComponent,
    #[cfg(feature = "egui")]
    pub bounding_box: egui::Rect,
}

impl GridLine {
    fn new(path: Vec<Complex64>, component: GridLineComponent) -> Self {
        #[cfg(feature = "egui")]
        {
            let mut x0 = path[0].re as f32;
            let mut y0 = -path[0].im as f32;
            let mut x1 = path[0].re as f32;
            let mut y1 = -path[0].im as f32;

            for p in path.iter() {
                x0 = x0.min(p.re as f32);
                y0 = y0.min(-p.im as f32);

                x1 = x1.max(p.re as f32);
                y1 = y1.max(-p.im as f32);
            }

            let bounding_box = egui::Rect {
                min: egui::pos2(x0, y0),
                max: egui::pos2(x1, y1),
            };

            Self {
                path,
                component,
                bounding_box,
            }
        }
        #[cfg(not(feature = "egui"))]
        {
            Self { path, component }
        }
    }
}

#[derive(Default)]
pub struct Contours {
    cuts: Vec<Cut>,
    commands: VecDeque<GeneratorCommand>,
    grid_p: Vec<GridLine>,
    grid_x: Vec<GridLine>,
    grid_u: Vec<GridLine>,

    rctx: ContourGeneratorRuntimeContext,

    num_commands: usize,
    loaded: bool,
}

fn branch_point_mass(p_start: f64, k: f64, branch_point_type: BranchPointType) -> f64 {
    match branch_point_type {
        BranchPointType::XpPositiveAxisImXmNegative => 2.0 * p_start * k + 2.0,
        BranchPointType::XpPositiveAxisImXmPositive => -(2.0 * p_start * k + 2.0),

        BranchPointType::XpNegativeAxisFromAboveWithImXmNegative => (2.0 * p_start + 1.0) * k + 2.0,
        BranchPointType::XpNegativeAxisFromBelowWithImXmNegative => (2.0 * p_start - 1.0) * k + 2.0,

        BranchPointType::XpNegativeAxisFromAboveWithImXmPositive => {
            -((2.0 * p_start + 1.0) * k + 2.0)
        }

        BranchPointType::XpNegativeAxisFromBelowWithImXmPositive => {
            -((2.0 * p_start - 1.0) * k + 2.0)
        }
    }
}

pub fn compute_branch_point(
    p_range: i32,
    branch_point_type: BranchPointType,
    consts: CouplingConstants,
) -> Option<BranchPointData> {
    let p_start = p_range as f64;
    let k = consts.k() as f64;
    let s = consts.s();
    let u_of_x = |x: Complex64| -> Complex64 { x + 1.0 / x - (s - 1.0 / s) * x.ln() };
    let du_dx = |x: Complex64| -> Complex64 { (x - s) * (x + 1.0 / s) / (x * x) };

    let u_of_s = u_of_x(Complex64::from(s))
        * match branch_point_type {
            BranchPointType::XpPositiveAxisImXmNegative
            | BranchPointType::XpPositiveAxisImXmPositive => 1.0,
            BranchPointType::XpNegativeAxisFromAboveWithImXmNegative
            | BranchPointType::XpNegativeAxisFromAboveWithImXmPositive
            | BranchPointType::XpNegativeAxisFromBelowWithImXmNegative
            | BranchPointType::XpNegativeAxisFromBelowWithImXmPositive => -1.0,
        };

    let m = branch_point_mass(p_start, k, branch_point_type);
    let guess = xp(0.5, m, consts);

    let x_branch_point = nr::find_root(
        |x| u_of_x(x) - u_of_s - m * Complex64::i() / consts.h,
        du_dx,
        guess,
        1.0e-3,
        10,
    );

    if let Some(x_branch_point) = x_branch_point {
        let p = x_branch_point.arg().abs() / std::f64::consts::PI;
        Some(BranchPointData {
            p,
            m,
            typ: branch_point_type,
        })
    } else {
        log::info!("Could not find branch point");
        None
    }
}

impl Contours {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn update(&mut self, p_range: i32, consts: CouplingConstants) -> bool {
        if self.num_commands == 0 {
            self.clear();
            self.commands = ContourCommandGenerator::generate_commands(p_range, consts);
            self.num_commands = self.commands.len();
            log::debug!("Generated {} commands", self.num_commands,)
        }

        if !self.loaded {
            if let Some(command) = self.commands.pop_front() {
                self.execute(command, consts);
            } else {
                self.cuts.sort_unstable_by_key(|cut| match cut.typ {
                    CutType::Log(_) => 2,
                    CutType::ULongNegative(_) => 3,
                    CutType::ULongPositive(_) => 4,
                    CutType::UShortScallion(_) => 5,
                    CutType::UShortKidney(_) => 6,
                    CutType::E => {
                        if cut.component == Component::P {
                            7
                        } else {
                            1
                        }
                    }
                    CutType::DebugPath => 8,
                });
                self.loaded = true;
            }
        }
        self.loaded
    }

    pub fn clear(&mut self) {
        log::debug!("Clearing grid and cuts");
        self.commands.clear();
        self.num_commands = 0;
        self.grid_x.clear();
        self.grid_u.clear();
        self.cuts.clear();
        self.loaded = false;

        self.grid_p = vec![GridLine::new(
            vec![
                Complex64::from(P_RANGE_MIN as f64),
                Complex64::from(P_RANGE_MAX as f64 + 1.0),
            ],
            GridLineComponent::Real,
        )];
    }

    pub fn progress(&self) -> (usize, usize) {
        if self.num_commands > 0 {
            (self.num_commands - self.commands.len(), self.num_commands)
        } else {
            (0, 1)
        }
    }

    pub fn get_grid(&self, component: Component) -> &Vec<GridLine> {
        match component {
            Component::P => &self.grid_p,
            Component::Xp | Component::Xm => &self.grid_x,
            Component::U => &self.grid_u,
        }
    }
    pub fn get_visible_cuts(
        &self,
        state: &State,
        component: Component,
        u_cut_type: UCutType,
        consts: CouplingConstants,
    ) -> impl Iterator<Item = &Cut> {
        let mut pt = state.active_point().clone();
        pt.u += 2.0 * (pt.sheet_data.log_branch_p * consts.k()) as f64 * Complex64::i() / consts.h;

        self.cuts
            .iter()
            .filter(move |c| c.component == component && c.is_visible(&pt, u_cut_type))
    }

    pub fn get_crossed_cuts(
        &self,
        pt: &Point,
        component: Component,
        new_value: Complex64,
        consts: CouplingConstants,
        u_cut_type: UCutType,
    ) -> impl Iterator<Item = &Cut> {
        let mut pt = pt.clone();
        pt.u += 2.0 * (pt.sheet_data.log_branch_p * consts.k()) as f64 * Complex64::i() / consts.h;

        let new_value = if component == Component::U {
            new_value
                + 2.0 * (pt.sheet_data.log_branch_p * consts.k()) as f64 * Complex64::i() / consts.h
        } else {
            new_value
        };

        self.cuts.iter().filter(move |c| {
            c.component == component
                && c.is_visible(&pt, u_cut_type)
                && c.intersection(pt.get(component), new_value, consts)
                    .is_some()
        })
    }

    fn execute(&mut self, command: GeneratorCommand, consts: CouplingConstants) {
        use GeneratorCommand::*;

        match command {
            AddGridLineU { y } => {
                self.grid_u.push(GridLine::new(
                    vec![Complex64::new(-INFINITY, y), Complex64::new(INFINITY, y)],
                    GridLineComponent::Real,
                ));
            }

            AddGridLineX { m } => {
                let path = XInterpolator::generate_xp_full(0, m, consts);
                if path.len() > 1 {
                    self.grid_x.push(GridLine::new(
                        path.iter().map(|x| x.conj()).collect(),
                        GridLineComponent::Xm(m),
                    ));
                    self.grid_x
                        .push(GridLine::new(path, GridLineComponent::Xp(m)));
                }
            }

            EStart { p_range } => {
                self.rctx.e_int = Some(EPInterpolator::new(p_range, consts));
            }

            PStartXp { p } => {
                self.rctx.p_int = Some(PInterpolatorMut::xp(p, consts));
            }

            PGotoXp(p, m) => {
                let Some(ref mut p_int) = self.rctx.p_int else { return };
                p_int.goto_xp(p, m);
            }

            PGotoXm(p, m) => {
                let Some(ref mut p_int) = self.rctx.p_int else { return };
                p_int.goto_xm(p, m);
            }

            PGotoRe(x) => {
                let Some(ref mut p_int) = self.rctx.p_int else { return };
                p_int.goto_re(x);
            }

            PGotoIm(x) => {
                let Some(ref mut p_int) = self.rctx.p_int else { return };
                p_int.goto_im(x);
            }

            PGotoP(p) => {
                let Some(ref mut p_int) = self.rctx.p_int else { return };
                p_int.goto_p(p);
            }

            PGotoM(m) => {
                let Some(ref mut p_int) = self.rctx.p_int else { return };
                p_int.goto_m(m);
            }

            AddGridLineP => {
                let Some(ref mut p_int) = self.rctx.p_int else { return };
                let path = p_int.contour();

                let (component, conj_component) = match p_int.pt() {
                    InterpolationPoint::Xp(_, m) => {
                        (GridLineComponent::Xp(m), GridLineComponent::Xm(m))
                    }
                    InterpolationPoint::Xm(_, m) => {
                        (GridLineComponent::Xm(m), GridLineComponent::Xp(m))
                    }
                    InterpolationPoint::Re(_) => (GridLineComponent::Real, GridLineComponent::Real),
                    _ => {
                        log::warn!("Cannot draw grid line for C");
                        return;
                    }
                };

                if path.len() > 1 {
                    self.grid_p.push(GridLine::new(
                        path.iter().map(|p| p.conj()).collect(),
                        conj_component,
                    ));
                    self.grid_p.push(GridLine::new(path, component));
                }
            }

            ClearCut => {
                self.rctx.cut_data.path = None;
                self.rctx.cut_data.branch_point = None;
            }

            ComputeCutP { reverse } => {
                let Some(ref mut p_int) = self.rctx.p_int else { return };
                let new_path = if reverse {
                    p_int.contour().into_iter().rev().collect()
                } else {
                    p_int.contour()
                };

                if let Some(ref mut path) = self.rctx.cut_data.path {
                    path.extend(new_path);
                } else {
                    self.rctx.cut_data.path = Some(new_path);
                }
            }

            ComputeBranchPoint {
                p_range,
                branch_point_type,
            } => {
                self.rctx.branch_point_data =
                    compute_branch_point(p_range, branch_point_type, consts);
            }

            ComputeCutX(cut_direction) => {
                self.rctx.cut_data.path = None;
                self.rctx.cut_data.branch_point = None;

                let Some(BranchPointData { p: p_branch_point, m, typ: branch_point_type }) = self.rctx.branch_point_data else {
                    log::warn!("No branch point set");
                    return;
                };

                let (p_start, p_end) = match cut_direction {
                    CutDirection::Positive => (0.0, p_branch_point),
                    CutDirection::Negative => (p_branch_point, 1.0),
                };

                let path = match branch_point_type {
                    BranchPointType::XpPositiveAxisImXmNegative
                    | BranchPointType::XpNegativeAxisFromAboveWithImXmNegative
                    | BranchPointType::XpNegativeAxisFromBelowWithImXmNegative => {
                        XInterpolator::generate_xm(p_start, p_end, m, consts)
                    }

                    BranchPointType::XpPositiveAxisImXmPositive
                    | BranchPointType::XpNegativeAxisFromAboveWithImXmPositive
                    | BranchPointType::XpNegativeAxisFromBelowWithImXmPositive => {
                        XInterpolator::generate_xp(p_start, p_end, m, consts)
                    }
                };

                let branch_point = *match cut_direction {
                    CutDirection::Positive => path.last().unwrap(),
                    CutDirection::Negative => path.first().unwrap(),
                };

                self.rctx.cut_data.path = Some(path);
                self.rctx.cut_data.branch_point = Some(branch_point);
            }

            ComputeCutXFull(xcut) => {
                self.rctx.cut_data.path = None;
                self.rctx.cut_data.branch_point = None;

                let m = match xcut {
                    XCut::Scallion => 0.0,
                    XCut::Kidney => -consts.k() as f64,
                };

                let half_path = XInterpolator::generate_xp_full(0, m, consts);
                let mut path = half_path.clone();
                path.extend(half_path.iter().map(|x| x.conj()).rev());

                self.rctx.cut_data.path = Some(path);
                self.rctx.cut_data.branch_point = Some(match xcut {
                    XCut::Scallion => Complex64::from(consts.s()),
                    XCut::Kidney => Complex64::from(-1.0 / consts.s()),
                });
            }

            ComputeCutEP => {
                let Some(ref mut e_int) = self.rctx.e_int else {return};
                let (branch_point, path) = e_int.get_cut_p();
                self.rctx.cut_data.path = path;
                self.rctx.cut_data.branch_point = branch_point;
            }

            ComputeCutEXp => {
                let Some(ref mut e_int) = self.rctx.e_int else {return};
                let (branch_point, path) = e_int.get_cut_xp();
                self.rctx.cut_data.path = path;
                self.rctx.cut_data.branch_point = branch_point;
            }

            ComputeCutEXm => {
                let Some(ref mut e_int) = self.rctx.e_int else {return};
                let (branch_point, path) = e_int.get_cut_xm();
                self.rctx.cut_data.path = path;
                self.rctx.cut_data.branch_point = branch_point;
            }

            ComputeCutEU => {
                let Some(ref mut e_int) = self.rctx.e_int else {return};
                let (branch_point, path) = e_int.get_cut_u();
                self.rctx.cut_data.path = path;
                self.rctx.cut_data.branch_point = branch_point;
            }

            SetCutPath { path, branch_point } => {
                self.rctx.cut_data.path = Some(path);
                self.rctx.cut_data.branch_point = branch_point;
            }

            PushCut {
                p_range,
                component,
                cut_type,
                periodic,
                visibility,
                pre_shift,
            } => {
                let Some(ref path) = self.rctx.cut_data.path else {
                    log::warn!("No path for cut");
                    return;
                };

                if self.rctx.cut_data.path.is_none() {
                    log::warn!("No path for cut");
                    return;
                };

                let shift = match component {
                    Component::U => Complex64::new(0.0, (p_range * consts.k()) as f64 / consts.h),
                    _ => Complex64::from(0.0),
                };

                let paths = path.iter().map(|z| z + pre_shift).collect();

                let cut = Cut::new(
                    component,
                    paths,
                    self.rctx.cut_data.branch_point.map(|z| z + pre_shift),
                    cut_type,
                    p_range,
                    periodic,
                    visibility,
                );

                self.cuts.push(cut.conj().shift(shift));
                self.cuts.push(cut.shift(shift));
            }

            SplitCut(component, cut_typ, branch_point, order) => {
                let Some(mut cut) = self.cuts.pop() else {return};
                let Some(_) = self.cuts.pop() else {return};
                let Some(ref path) = self.rctx.cut_data.path else { return };

                let shift = match cut.component {
                    Component::U => {
                        Complex64::new(0.0, (cut.p_range * consts.k()) as f64 / consts.h)
                    }
                    _ => Complex64::from(0.0),
                };

                for (p1, p2) in path
                    .iter()
                    .map(|p| {
                        if component == Component::Xp {
                            p + shift
                        } else {
                            p.conj() + shift
                        }
                    })
                    .tuple_windows::<(_, _)>()
                {
                    if let Some((j, x)) = cut.intersection(p1, p2, consts) {
                        let mut new_path = vec![x];
                        new_path.extend(cut.path.split_off(j + 1));
                        cut.path.push(x);

                        let mut new_cut = Cut {
                            path: new_path,
                            branch_point: None,
                            typ: cut.typ.clone(),
                            p_range: cut.p_range,
                            component: cut.component,
                            periodic: false,
                            visibility: vec![],
                        };
                        if branch_point == SplitCutBranchPoint::New && cut.branch_point.is_some() {
                            new_cut.branch_point = cut.branch_point;
                            cut.branch_point = None;
                        }
                        for vis in cut.visibility.iter() {
                            let vis = match vis {
                                CutVisibilityCondition::UpBranch(b) => {
                                    if component == Component::Xp {
                                        let b = b.clone().cross(&cut_typ);
                                        CutVisibilityCondition::UpBranch(b)
                                    } else {
                                        CutVisibilityCondition::UpBranch(b.clone())
                                    }
                                }
                                CutVisibilityCondition::UmBranch(b) => {
                                    if component == Component::Xm {
                                        let b = b.clone().cross(&cut_typ);
                                        CutVisibilityCondition::UmBranch(b)
                                    } else {
                                        CutVisibilityCondition::UmBranch(b.clone())
                                    }
                                }
                                _ => vis.clone(),
                            };
                            new_cut.visibility.push(vis);
                        }
                        log::debug!("Intersection found");

                        if order == SplitCutOrder::NewFirst {
                            self.cuts.push(new_cut.shift_conj(shift));
                            self.cuts.push(new_cut);
                            self.cuts.push(cut.shift_conj(shift));
                            self.cuts.push(cut);
                        } else {
                            self.cuts.push(cut.shift_conj(shift));
                            self.cuts.push(cut);
                            self.cuts.push(new_cut.shift_conj(shift));
                            self.cuts.push(new_cut);
                        }

                        return;
                    }
                }

                log::warn!("No intersection found");

                self.cuts.push(cut.conj());
                self.cuts.push(cut);
            }

            PopCut => {
                let Some(cut) = self.cuts.pop() else {return};
                let Some(_) = self.cuts.pop() else {return};

                let shift = match cut.component {
                    Component::U => {
                        Complex64::new(0.0, (cut.p_range * consts.k()) as f64 / consts.h)
                    }
                    _ => Complex64::from(0.0),
                };

                self.rctx.cut_data.path = Some(cut.path.iter().map(|z| z - shift).collect());
                self.rctx.cut_data.branch_point = cut.branch_point.map(|z| z - shift);
            }

            SwapCuts => {
                if self.cuts.len() < 4 {
                    return;
                }
                let len = self.cuts.len();
                self.cuts.swap(len - 4, len - 2);
                self.cuts.swap(len - 3, len - 1);
            }
        }
    }
}

impl ContourCommandGenerator {
    fn generate_commands(p_range: i32, consts: CouplingConstants) -> VecDeque<GeneratorCommand> {
        let bctx = Self::new();
        bctx.do_generate_commands(p_range, consts)
    }

    fn new() -> Self {
        Self {
            component: None,
            cut_type: None,
            periodic: false,
            visibility: vec![],
            group_visibility: vec![],
            pre_shift: Complex64::from(0.0),
            commands: VecDeque::new(),
        }
    }

    fn reset(&mut self) {
        self.component = None;
        self.cut_type = None;
        self.periodic = false;
        self.visibility = self.group_visibility.clone();
        self.pre_shift = Complex64::from(0.0);
    }

    fn add(&mut self, command: GeneratorCommand) -> &mut Self {
        self.commands.push_back(command);
        self
    }

    fn e_start(&mut self, p_range: i32) -> &mut Self {
        self.add(GeneratorCommand::EStart { p_range })
    }

    fn compute_cut_e_p(&mut self) -> &mut Self {
        self.add(GeneratorCommand::ComputeCutEP)
    }

    fn compute_cut_e_xp(&mut self) -> &mut Self {
        self.add(GeneratorCommand::ComputeCutEXp)
    }

    fn compute_cut_e_xm(&mut self) -> &mut Self {
        self.add(GeneratorCommand::ComputeCutEXm)
    }

    fn compute_cut_e_u(&mut self) -> &mut Self {
        self.add(GeneratorCommand::ComputeCutEU)
    }

    fn p_start_xp(&mut self, p: f64) -> &mut Self {
        self.add(GeneratorCommand::PStartXp { p })
    }

    fn goto_xp(&mut self, p: f64, m: f64) -> &mut Self {
        self.add(GeneratorCommand::PGotoXp(p, m))
    }

    fn goto_xm(&mut self, p: f64, m: f64) -> &mut Self {
        self.add(GeneratorCommand::PGotoXm(p, m))
    }

    fn goto_re(&mut self, re: f64) -> &mut Self {
        self.add(GeneratorCommand::PGotoRe(re))
    }

    fn goto_im(&mut self, im: f64) -> &mut Self {
        self.add(GeneratorCommand::PGotoIm(im))
    }

    fn goto_p(&mut self, p: f64) -> &mut Self {
        self.add(GeneratorCommand::PGotoP(p))
    }

    fn goto_m(&mut self, m: f64) -> &mut Self {
        self.add(GeneratorCommand::PGotoM(m))
    }

    fn p_grid_line(&mut self) -> &mut Self {
        self.add(GeneratorCommand::AddGridLineP)
    }

    fn do_generate_commands(
        mut self,
        p_range: i32,
        consts: CouplingConstants,
    ) -> VecDeque<GeneratorCommand> {
        self.generate_u_grid(consts);

        let max = P_RANGE_MAX - P_RANGE_MIN;

        self.generate_cuts(p_range, consts);

        for i in 1..max {
            if p_range - i >= P_RANGE_MIN {
                self.generate_cuts(p_range - i, consts);
            }

            if p_range + i <= P_RANGE_MAX {
                self.generate_cuts(p_range + i, consts);
            }
        }

        self.generate_x_grid(p_range, consts);
        for i in 1..max {
            if p_range - i >= P_RANGE_MIN {
                self.generate_x_grid(p_range - i, consts);
            }

            if p_range + i <= P_RANGE_MAX {
                self.generate_x_grid(p_range + i, consts);
            }
        }
        self.generate_p_grid(p_range, consts);

        for i in 1..max {
            if p_range - i >= P_RANGE_MIN {
                self.generate_p_grid(p_range - i, consts);
            }

            if p_range + i <= P_RANGE_MAX {
                self.generate_p_grid(p_range + i, consts);
            }
        }

        self.commands
    }

    fn generate_u_grid(&mut self, consts: CouplingConstants) {
        self.add(GeneratorCommand::AddGridLineU { y: 0.0 });

        for y in 1..=100 {
            self.add(GeneratorCommand::AddGridLineU {
                y: y as f64 / consts.h,
            });
            self.add(GeneratorCommand::AddGridLineU {
                y: -y as f64 / consts.h,
            });
        }
    }

    fn generate_x_grid(&mut self, p_range: i32, _consts: CouplingConstants) {
        if p_range == 0 {
            for m in -50..50 {
                let m = m as f64;
                self.add(GeneratorCommand::AddGridLineX { m });
            }
        }
    }

    fn generate_p_grid(&mut self, p_range: i32, consts: CouplingConstants) {
        let p_start = p_range as f64;
        let k = consts.k() as f64;
        const M_MAX: i32 = 60;
        const M_MIN: i32 = 20;
        {
            let p0 = p_start + 1.0 / 16.0;
            let p2 = p_start + 15.0 / 16.0;

            self.p_start_xp(p0).goto_xp(p0, -p_start * k).p_grid_line();

            self.p_start_xp(p0)
                .goto_xm(p0, 1.0)
                .goto_xm(p0, -p_start * k)
                .p_grid_line();

            self.p_start_xp(p2)
                .goto_xp(p2, -(p_start + 1.0) * k)
                .p_grid_line();

            self.p_start_xp(p0)
                .goto_xm(p0, 1.0)
                .goto_xm(p2, 1.0)
                .goto_xm(p2, -(p_start + 1.0) * k)
                .p_grid_line();
        }

        if p_range == 0 {
            let p0 = p_start + 1.0 / 16.0;
            let p2 = p_start + 15.0 / 16.0;

            self.p_start_xp(p0);

            for m in 3..=M_MIN {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            self.p_start_xp(p2).goto_xp(p2, 3.0);

            for m in 3..=(consts.k() + 1) {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            for m in (consts.k() + 3)..=M_MAX {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            self.p_start_xp(p0).goto_xp(p0, 3.0).goto_xp(p2, 3.0);

            for m in ((3 - consts.k())..=1).rev() {
                self.goto_xp(p2, m as f64).p_grid_line();
            }
        }

        if p_range > 0 {
            let p0 = p_start + 1.0 / 16.0;
            let p2 = p_start + 15.0 / 16.0;

            self.p_start_xp(p0);

            for m in 2..=(p_range * consts.k() + 1) {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            for m in (p_range * consts.k() + 3)..=(2 + (2 * p_range + 2) * consts.k()) {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            if k > 0.0 {
                self.p_start_xp(p2).goto_xp(p2, p_start * k + 3.0);

                for m in (p_range * consts.k() + 3)..=((p_range + 1) * consts.k() + 1) {
                    self.goto_xp(p0, m as f64).p_grid_line();
                }

                for m in ((p_range + 1) * consts.k() + 3)..=M_MAX {
                    self.goto_xp(p0, m as f64).p_grid_line();
                }
            } else {
                self.p_start_xp((p0 + p2) / 2.0).goto_m(3.0);

                for m in 3..=M_MAX {
                    self.goto_m(m as f64).p_grid_line();
                }
            }

            self.p_start_xp(p0)
                .goto_xp(p0, p_start * k + 3.0)
                .goto_xp(p2, p_start * k + 3.0)
                .goto_xp(p2, p_start * k + 1.0);

            for m in (((p_range - 1) * consts.k() + 3)..=(p_range * consts.k() + 1)).rev() {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            for m in (1..=((p_range - 1) * consts.k() + 1)).rev() {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            self.goto_xp(p2, 1.0);

            for m in ((-(consts.k()) + 2)..=0).rev() {
                self.goto_xp(p2, m as f64).p_grid_line();
            }
        }

        if p_range == -1 {
            let p0 = p_start + 1.0 / 16.0;

            self.p_start_xp(p0);

            for m in 3..=(consts.k() - 1) {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            for m in (consts.k() + 1)..=M_MAX {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            self.p_start_xp(p0).goto_xm(p0, 1.0);

            for m in 1..=(consts.k() - 1) {
                self.goto_xm(p0, m as f64).p_grid_line();
            }

            for m in (consts.k() + 1)..=(2 * consts.k() - 2) {
                self.goto_xm(p0, m as f64).p_grid_line();
            }

            self.p_start_xp(p0).goto_xm(p0, 1.0);
            for m in ((-2 * consts.k())..=-1).rev() {
                self.goto_xm(p0, m as f64).p_grid_line();
            }
        }

        if p_range < -1 {
            let p0 = p_start + 1.0 / 16.0;

            self.p_start_xp(p0);

            for m in 2..=(-(p_range + 1) * consts.k() - 1) {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            for m in (-(p_range + 1) * consts.k() + 1)..=(-p_range * consts.k() - 1) {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            for m in (-p_range * consts.k() + 1)..=M_MAX {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            self.p_start_xp(p0);

            for m in 1..=(-(p_range + 1) * consts.k() - 1) {
                self.goto_xm(p0, m as f64).p_grid_line();
            }

            for m in (-(p_range + 1) * consts.k() + 1)..=(-p_range * consts.k() - 1) {
                self.goto_xm(p0, m as f64).p_grid_line();
            }

            for m in (-p_range * consts.k() + 1)..=M_MAX {
                self.goto_xm(p0, m as f64).p_grid_line();
            }

            self.p_start_xp(p0).goto_xm(p0, 1.0);

            for m in ((-2 * consts.k())..=0).rev() {
                self.goto_xm(p0, m as f64).p_grid_line();
            }
        }

        {
            // Real positive line
            let p0 = p_start + 1.0 / 8.0;
            self.p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .p_grid_line();
        }
        {
            // Real negative line
            let p0 = p_start + 7.0 / 8.0;
            self.p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .p_grid_line();
        }
    }

    fn begin_group(&mut self, group_visibility: &[CutVisibilityCondition]) -> &mut Self {
        self.group_visibility = Vec::from(group_visibility);
        self.reset();
        self
    }

    fn clear_cut(&mut self) -> &mut Self {
        self.add(GeneratorCommand::ClearCut)
    }

    fn compute_branch_point(
        &mut self,
        p_range: i32,
        branch_point_type: BranchPointType,
    ) -> &mut Self {
        self.add(GeneratorCommand::ComputeBranchPoint {
            p_range,
            branch_point_type,
        })
    }

    fn compute_cut_path_x(&mut self, direction: CutDirection) -> &mut Self {
        self.add(GeneratorCommand::ComputeCutX(direction))
    }

    fn compute_cut_path_x_full(&mut self, xcut: XCut) -> &mut Self {
        self.add(GeneratorCommand::ComputeCutXFull(xcut))
    }

    fn compute_cut_path_p(&mut self) -> &mut Self {
        self.add(GeneratorCommand::ComputeCutP { reverse: false })
    }

    fn compute_cut_path_p_rev(&mut self) -> &mut Self {
        self.add(GeneratorCommand::ComputeCutP { reverse: true })
    }

    fn set_cut_path(&mut self, path: Vec<Complex64>, branch_point: Option<Complex64>) -> &mut Self {
        self.add(GeneratorCommand::SetCutPath { path, branch_point })
    }

    fn push_cut(&mut self, p_range: i32) -> &mut Self {
        let Some(component) = std::mem::replace(&mut self.component, None) else {
            log::warn!("Can't push cut without component");
            self.reset();
            return self;
        };
        let Some(cut_type) = std::mem::replace(&mut self.cut_type, None) else {
            log::warn!("Can't push cut without type");
            self.reset();
            return self;
        };
        self.add(GeneratorCommand::PushCut {
            p_range,
            component,
            cut_type,
            periodic: self.periodic,
            visibility: self.visibility.clone(),
            pre_shift: self.pre_shift,
        });
        self.reset();
        self
    }

    fn split_cut(
        &mut self,
        component: Component,
        cut_type: CutType,
        branch_point: SplitCutBranchPoint,
        order: SplitCutOrder,
    ) -> &mut Self {
        self.add(GeneratorCommand::SplitCut(
            component,
            cut_type,
            branch_point,
            order,
        ))
    }

    fn pop_cut(&mut self) -> &mut Self {
        self.add(GeneratorCommand::PopCut)
    }

    fn swap_cuts(&mut self) -> &mut Self {
        self.add(GeneratorCommand::SwapCuts)
    }

    fn create_cut(&mut self, component: Component, cut_type: CutType) -> &mut Self {
        if self.component.is_some() || self.cut_type.is_some() {
            log::warn!("New cut created before previous cut was pushed");
        }
        self.reset();
        self.component = Some(component);
        self.cut_type = Some(cut_type);
        self
    }

    fn log_branch(&mut self, p_range: i32) -> &mut Self {
        self.visibility
            .push(CutVisibilityCondition::LogBranch(p_range));
        self
    }

    fn im_xm_negative(&mut self) -> &mut Self {
        self.visibility.push(CutVisibilityCondition::ImXm(-1));
        self
    }

    fn im_xm_positive(&mut self) -> &mut Self {
        self.visibility.push(CutVisibilityCondition::ImXm(1));
        self
    }

    fn im_xp_positive(&mut self) -> &mut Self {
        self.visibility.push(CutVisibilityCondition::ImXp(1));
        self
    }

    fn im_xp_negative(&mut self) -> &mut Self {
        self.visibility.push(CutVisibilityCondition::ImXp(-1));
        self
    }

    fn xp_outside(&mut self) -> &mut Self {
        self.visibility
            .push(CutVisibilityCondition::UpBranch(UBranch::Outside));
        self
    }

    fn xp_between(&mut self) -> &mut Self {
        self.visibility
            .push(CutVisibilityCondition::UpBranch(UBranch::Between));
        self
    }

    fn xp_inside(&mut self) -> &mut Self {
        self.visibility
            .push(CutVisibilityCondition::UpBranch(UBranch::Inside));
        self
    }

    fn xm_outside(&mut self) -> &mut Self {
        self.visibility
            .push(CutVisibilityCondition::UmBranch(UBranch::Outside));
        self
    }

    fn xm_between(&mut self) -> &mut Self {
        self.visibility
            .push(CutVisibilityCondition::UmBranch(UBranch::Between));
        self
    }

    fn xm_inside(&mut self) -> &mut Self {
        self.visibility
            .push(CutVisibilityCondition::UmBranch(UBranch::Inside));
        self
    }

    fn e_branch(&mut self, branch: i32) -> &mut Self {
        self.visibility
            .push(CutVisibilityCondition::EBranch(branch));
        self
    }

    fn periodic(&mut self) -> &mut Self {
        self.periodic = true;
        self
    }

    fn pre_shift(&mut self, z: Complex64) -> &mut Self {
        self.pre_shift = z;
        self
    }

    fn generate_cuts(&mut self, p_range: i32, consts: CouplingConstants) {
        self.generate_long_cuts(p_range, consts);
        self.generate_semishort_cuts(p_range, consts);
        self.generate_short_cuts(p_range, consts);
    }

    fn generate_long_cuts(&mut self, p_range: i32, consts: CouplingConstants) {
        let p_start = p_range as f64;
        let k = consts.k() as f64;
        let s = consts.s();

        let us = s + 1.0 / s - (s - 1.0 / s) * s.ln();

        self.begin_group(&[CutVisibilityCondition::LongCuts]);

        {
            // Log
            self.clear_cut()
                .set_cut_path(
                    vec![Complex64::from(-INFINITY), Complex64::from(0.0)],
                    Some(Complex64::from(0.0)),
                )
                .create_cut(Component::Xp, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(-INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range + 1) as f64 * k) / consts.h,
                    )),
                )
                .create_cut(Component::U, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(-INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range - 1) as f64 * k) / consts.h,
                    )),
                )
                .create_cut(Component::U, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            // Real negative line
            let p0 = p_start + 7.0 / 8.0;
            self.clear_cut()
                .p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .goto_re(-consts.s() * 4.0)
                .compute_cut_path_p()
                .create_cut(Component::P, CutType::Log(Component::Xp))
                .e_branch(1)
                .push_cut(p_range);
        }

        {
            // U long positive
            self.clear_cut()
                .set_cut_path(
                    vec![Complex64::from(0.0), Complex64::from(INFINITY)],
                    Some(Complex64::from(s)),
                )
                .create_cut(Component::Xp, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmNegative)
                .compute_cut_path_x(CutDirection::Positive)
                .create_cut(Component::Xm, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmPositive)
                .compute_cut_path_x(CutDirection::Positive)
                .create_cut(Component::Xm, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                        Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
                )
                .create_cut(Component::U, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            // Real positive line
            let p0 = p_start + 1.0 / 8.0;
            self.clear_cut()
                .p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .goto_re(consts.s() * 4.0)
                .compute_cut_path_p()
                .create_cut(Component::P, CutType::ULongPositive(Component::Xp))
                .e_branch(1)
                .push_cut(p_range);

            self.create_cut(Component::P, CutType::Log(Component::Xm))
                .e_branch(-1)
                .push_cut(p_range);

            self.create_cut(Component::P, CutType::ULongNegative(Component::Xm))
                .e_branch(-1)
                .push_cut(p_range);
        }

        {
            // U long negative
            self.clear_cut()
                .set_cut_path(
                    vec![Complex64::from(-INFINITY), Complex64::from(0.0)],
                    Some(-1.0 / Complex64::from(s)),
                )
                .create_cut(Component::Xp, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(-INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range + 1) as f64 * k) / consts.h,
                    )),
                )
                .create_cut(Component::U, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(-INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range - 1) as f64 * k) / consts.h,
                    )),
                )
                .create_cut(Component::U, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            // Real negative line
            let p0 = p_start + 7.0 / 8.0;
            self.clear_cut()
                .p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .goto_re(-consts.s() * 4.0)
                .compute_cut_path_p()
                .create_cut(Component::P, CutType::ULongNegative(Component::Xp))
                .e_branch(1)
                .push_cut(p_range);

            self.create_cut(Component::P, CutType::ULongPositive(Component::Xm))
                .e_branch(-1)
                .push_cut(p_range);
        }

        {
            // Scallion
            self.clear_cut()
                .compute_cut_path_x_full(XCut::Scallion)
                .create_cut(Component::Xp, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            if p_range != 0 {
                // For p = 0 we add this cut after the E cut
                self.clear_cut()
                    .compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmPositive)
                    .compute_cut_path_x(CutDirection::Negative)
                    .create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                    .log_branch(p_range)
                    .push_cut(p_range);
            }

            if p_range != -1 {
                // For p = -1 we add this cut after the E cut
                self.clear_cut()
                    .compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmNegative)
                    .compute_cut_path_x(CutDirection::Negative)
                    .create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                    .log_branch(p_range)
                    .push_cut(p_range);
            }

            if p_range != 0 && p_range != -1 {
                // For p = 0, -1 we add this cut after the E cut
                self.clear_cut()
                    .set_cut_path(
                        vec![
                            Complex64::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                            Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                        ],
                        Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
                    )
                    .create_cut(Component::U, CutType::UShortScallion(Component::Xp))
                    .log_branch(p_range)
                    .push_cut(p_range);
            }

            let p0 = p_start + 1.0 / 8.0;
            self.clear_cut();

            self.p_start_xp(p0)
                .goto_m(-(p_range * consts.k()) as f64)
                .compute_cut_path_p_rev();

            self.p_start_xp(p0)
                .goto_xm(p0, 1.0)
                .goto_m(-(p_range * consts.k()) as f64)
                .compute_cut_path_p();

            self.create_cut(Component::P, CutType::UShortScallion(Component::Xp))
                .e_branch(1)
                .push_cut(p_range);

            self.create_cut(Component::P, CutType::UShortKidney(Component::Xm))
                .e_branch(-1)
                .push_cut(p_range);

            if p_range == 0 {
                let p0 = p_start + 1.0 / 8.0;
                let p1 = p_start + 7.0 / 8.0;

                self.clear_cut()
                    .p_start_xp(p0)
                    .goto_m(3.0)
                    .goto_p(p1)
                    .goto_m(0.0)
                    .compute_cut_path_p();

                self.create_cut(Component::P, CutType::UShortScallion(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);

                self.create_cut(Component::P, CutType::UShortKidney(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);
            }
        }

        {
            // Kidney

            self.clear_cut()
                .compute_cut_path_x_full(XCut::Kidney)
                .create_cut(Component::Xp, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Positive)
                .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Positive)
                .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);

            if p_range != 0 {
                // For p = 0 we add this after E
                self.clear_cut()
                    .compute_branch_point(
                        p_range,
                        BranchPointType::XpNegativeAxisFromBelowWithImXmNegative,
                    )
                    .compute_cut_path_x(CutDirection::Positive)
                    .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                    .log_branch(p_range)
                    .xp_between()
                    .push_cut(p_range);
            }

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Positive)
                .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);

            if p_range != -1 {
                self.clear_cut()
                    .set_cut_path(
                        vec![
                            Complex64::new(INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                            Complex64::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                        ],
                        Some(Complex64::new(
                            -us,
                            -(1.0 + (p_range + 1) as f64 * k) / consts.h,
                        )),
                    )
                    .create_cut(Component::U, CutType::UShortKidney(Component::Xp))
                    .log_branch(p_range)
                    .xp_between()
                    .push_cut(p_range);
            }

            if p_range != 0 {
                self.clear_cut()
                    .set_cut_path(
                        vec![
                            Complex64::new(INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                            Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        ],
                        Some(Complex64::new(
                            -us,
                            -(1.0 + (p_range - 1) as f64 * k) / consts.h,
                        )),
                    )
                    .create_cut(Component::U, CutType::UShortKidney(Component::Xp))
                    .log_branch(p_range)
                    .xp_between()
                    .push_cut(p_range);
            }

            let p0 = p_start + 1.0 / 8.0;
            let p1 = p_start + 7.0 / 8.0;

            if p_range != 0 {
                self.clear_cut();
                self.p_start_xp(p0)
                    .goto_m(-((p_range + 1) * consts.k()) as f64)
                    .compute_cut_path_p_rev();

                self.p_start_xp(p0)
                    .goto_xm(p0, 1.0)
                    .goto_p(p1)
                    .goto_m(-((p_range + 1) * consts.k()) as f64)
                    .compute_cut_path_p();

                self.create_cut(Component::P, CutType::UShortKidney(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);
                self.create_cut(Component::P, CutType::UShortScallion(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);
            } else {
                self.clear_cut();
                self.p_start_xp(p0)
                    .goto_m(-((p_range + 1) * consts.k()) as f64)
                    .compute_cut_path_p_rev();

                self.create_cut(Component::P, CutType::UShortKidney(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);
                self.create_cut(Component::P, CutType::UShortScallion(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);

                self.clear_cut();

                self.p_start_xp(p0)
                    .goto_xm(p0, 1.0)
                    .goto_p(p1)
                    .goto_m(-((p_range + 1) * consts.k()) as f64)
                    .compute_cut_path_p();

                self.create_cut(Component::P, CutType::UShortKidney(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);
                self.create_cut(Component::P, CutType::UShortScallion(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);
            }

            if p_range == -1 {
                let p0 = p_start + 1.0 / 8.0;
                let p1 = p_start + 6.0 / 8.0;

                self.clear_cut()
                    .p_start_xp(p0)
                    .goto_m(-consts.k() as f64 + 1.0)
                    .goto_p(p1)
                    .goto_m(-consts.k() as f64 + 3.0)
                    .goto_p(p0)
                    .goto_m(0.0)
                    .compute_cut_path_p();
                self.create_cut(Component::P, CutType::UShortKidney(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);
                self.create_cut(Component::P, CutType::UShortScallion(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);
            }
        }

        self.e_start(p_range);

        self.compute_cut_e_p()
            .create_cut(Component::P, CutType::E)
            .push_cut(p_range);

        self.compute_cut_e_xp();

        self.create_cut(Component::Xp, CutType::E)
            .log_branch(p_range)
            .im_xm_negative()
            .push_cut(p_range);

        #[allow(clippy::comparison_chain)]
        if p_range == 0 {
            self.compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmPositive)
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);
        }

        self.compute_cut_e_xm();
        self.create_cut(Component::Xm, CutType::E)
            .log_branch(p_range)
            .im_xp_negative()
            .push_cut(p_range);

        if p_range == 0 {
            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Positive)
                .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);
        } else if p_range == -1 {
            self.compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmNegative)
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);
        }

        self.compute_cut_e_u()
            .create_cut(Component::U, CutType::E)
            .log_branch(p_range)
            .im_xp_negative()
            .im_xm_negative()
            .push_cut(p_range);

        if p_range == 0 {
            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                        Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
                )
                .create_cut(Component::U, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range - 1) as f64 * k) / consts.h,
                    )),
                )
                .create_cut(Component::U, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);
        } else if p_range == -1 {
            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range + 1) as f64 * k) / consts.h,
                    )),
                )
                .create_cut(Component::U, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                        Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
                )
                .create_cut(Component::U, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);
        }
    }

    fn generate_semishort_cuts(&mut self, p_range: i32, consts: CouplingConstants) {
        let p_start = p_range as f64;
        let k = consts.k() as f64;
        let s = consts.s();

        let us = s + 1.0 / s - (s - 1.0 / s) * s.ln();

        self.begin_group(&[CutVisibilityCondition::SemiShortCuts]);

        {
            // Log
            self.clear_cut()
                .set_cut_path(
                    vec![Complex64::from(-INFINITY), Complex64::from(0.0)],
                    Some(Complex64::from(0.0)),
                )
                .create_cut(Component::Xp, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);

            // Real negative line
            let p0 = p_start + 7.0 / 8.0;
            self.clear_cut()
                .p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .goto_re(-consts.s() * 4.0)
                .compute_cut_path_p()
                .create_cut(Component::P, CutType::Log(Component::Xp))
                .e_branch(1)
                .push_cut(p_range);
        }

        {
            // U long positive
            self.clear_cut()
                .set_cut_path(
                    vec![Complex64::from(0.0), Complex64::from(INFINITY)],
                    Some(Complex64::from(s)),
                )
                .create_cut(Component::Xp, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmNegative)
                .compute_cut_path_x(CutDirection::Positive)
                .create_cut(Component::Xm, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmPositive)
                .compute_cut_path_x(CutDirection::Positive)
                .create_cut(Component::Xm, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                        Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
                )
                .create_cut(Component::U, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            // Real positive line
            let p0 = p_start + 1.0 / 8.0;
            self.clear_cut()
                .p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .goto_re(consts.s() * 4.0)
                .compute_cut_path_p()
                .create_cut(Component::P, CutType::ULongPositive(Component::Xp))
                .e_branch(1)
                .push_cut(p_range);

            self.create_cut(Component::P, CutType::Log(Component::Xm))
                .e_branch(-1)
                .push_cut(p_range);

            self.create_cut(Component::P, CutType::ULongNegative(Component::Xm))
                .e_branch(-1)
                .push_cut(p_range);
        }

        {
            // U long negative
            self.clear_cut()
                .set_cut_path(
                    vec![Complex64::from(-INFINITY), Complex64::from(0.0)],
                    Some(-1.0 / Complex64::from(s)),
                )
                .create_cut(Component::Xp, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Negative)
                .create_cut(Component::Xm, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(-INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range + 1) as f64 * k) / consts.h,
                    )),
                )
                .create_cut(Component::U, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(-INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range - 1) as f64 * k) / consts.h,
                    )),
                )
                .create_cut(Component::U, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            // Real negative line
            let p0 = p_start + 7.0 / 8.0;
            self.clear_cut()
                .p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .goto_re(-consts.s() * 4.0)
                .compute_cut_path_p()
                .create_cut(Component::P, CutType::ULongNegative(Component::Xp))
                .e_branch(1)
                .push_cut(p_range);

            self.create_cut(Component::P, CutType::ULongPositive(Component::Xm))
                .e_branch(-1)
                .push_cut(p_range);
        }

        {
            // Scallion
            self.clear_cut()
                .compute_cut_path_x_full(XCut::Scallion)
                .create_cut(Component::Xp, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            if p_range != 0 {
                // For p = 0 we add this cut after the E cut
                self.clear_cut()
                    .compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmPositive)
                    .compute_cut_path_x(CutDirection::Negative)
                    .create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                    .log_branch(p_range)
                    .push_cut(p_range);
            }

            if p_range != -1 {
                // For p = -1 we add this cut after the E cut
                self.clear_cut()
                    .compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmNegative)
                    .compute_cut_path_x(CutDirection::Negative)
                    .create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                    .log_branch(p_range)
                    .push_cut(p_range);
            }

            if p_range != 0 && p_range != -1 {
                // For p = 0, -1 we add this cut after the E cut
                self.clear_cut()
                    .set_cut_path(
                        vec![
                            Complex64::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                            Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                        ],
                        Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
                    )
                    .create_cut(Component::U, CutType::UShortScallion(Component::Xp))
                    .log_branch(p_range)
                    .push_cut(p_range);
            }

            let p0 = p_start + 1.0 / 8.0;
            self.clear_cut();

            self.p_start_xp(p0)
                .goto_m(-(p_range * consts.k()) as f64)
                .compute_cut_path_p_rev();

            self.p_start_xp(p0)
                .goto_xm(p0, 1.0)
                .goto_m(-(p_range * consts.k()) as f64)
                .compute_cut_path_p();

            self.create_cut(Component::P, CutType::UShortScallion(Component::Xp))
                .e_branch(1)
                .push_cut(p_range);

            self.create_cut(Component::P, CutType::UShortKidney(Component::Xm))
                .e_branch(-1)
                .push_cut(p_range);

            if p_range == 0 {
                let p0 = p_start + 1.0 / 8.0;
                let p1 = p_start + 7.0 / 8.0;

                self.clear_cut()
                    .p_start_xp(p0)
                    .goto_m(3.0)
                    .goto_p(p1)
                    .goto_m(0.0)
                    .compute_cut_path_p();

                self.create_cut(Component::P, CutType::UShortScallion(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);

                self.create_cut(Component::P, CutType::UShortKidney(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);
            }
        }

        {
            // Kidney

            self.clear_cut()
                .compute_cut_path_x_full(XCut::Kidney)
                .create_cut(Component::Xp, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Positive)
                .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Positive)
                .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);

            if p_range != 0 {
                // For p = 0 we add this after E
                self.clear_cut()
                    .compute_branch_point(
                        p_range,
                        BranchPointType::XpNegativeAxisFromBelowWithImXmNegative,
                    )
                    .compute_cut_path_x(CutDirection::Positive)
                    .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                    .log_branch(p_range)
                    .xp_between()
                    .push_cut(p_range);
            }

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Positive)
                .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);

            if p_range != -1 {
                self.clear_cut()
                    .set_cut_path(
                        vec![
                            Complex64::new(INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                            Complex64::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                        ],
                        Some(Complex64::new(
                            -us,
                            -(1.0 + (p_range + 1) as f64 * k) / consts.h,
                        )),
                    )
                    .create_cut(Component::U, CutType::UShortKidney(Component::Xp))
                    .log_branch(p_range)
                    .xp_between()
                    .push_cut(p_range);
            }

            if p_range != 0 {
                self.clear_cut()
                    .set_cut_path(
                        vec![
                            Complex64::new(INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                            Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        ],
                        Some(Complex64::new(
                            -us,
                            -(1.0 + (p_range - 1) as f64 * k) / consts.h,
                        )),
                    )
                    .create_cut(Component::U, CutType::UShortKidney(Component::Xp))
                    .log_branch(p_range)
                    .xp_between()
                    .push_cut(p_range);
            }

            let p0 = p_start + 1.0 / 8.0;
            let p1 = p_start + 7.0 / 8.0;

            if p_range != 0 {
                self.clear_cut();
                self.p_start_xp(p0)
                    .goto_m(-((p_range + 1) * consts.k()) as f64)
                    .compute_cut_path_p_rev();

                self.p_start_xp(p0)
                    .goto_xm(p0, 1.0)
                    .goto_p(p1)
                    .goto_m(-((p_range + 1) * consts.k()) as f64)
                    .compute_cut_path_p();

                self.create_cut(Component::P, CutType::UShortKidney(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);
                self.create_cut(Component::P, CutType::UShortScallion(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);
            } else {
                self.clear_cut();
                self.p_start_xp(p0)
                    .goto_m(-((p_range + 1) * consts.k()) as f64)
                    .compute_cut_path_p_rev();

                self.create_cut(Component::P, CutType::UShortKidney(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);
                self.create_cut(Component::P, CutType::UShortScallion(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);

                self.clear_cut();

                self.p_start_xp(p0)
                    .goto_xm(p0, 1.0)
                    .goto_p(p1)
                    .goto_m(-((p_range + 1) * consts.k()) as f64)
                    .compute_cut_path_p();

                self.create_cut(Component::P, CutType::UShortKidney(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);
                self.create_cut(Component::P, CutType::UShortScallion(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);
            }

            if p_range == -1 {
                let p0 = p_start + 1.0 / 8.0;
                let p1 = p_start + 6.0 / 8.0;

                self.clear_cut()
                    .p_start_xp(p0)
                    .goto_m(-consts.k() as f64 + 1.0)
                    .goto_p(p1)
                    .goto_m(-consts.k() as f64 + 3.0)
                    .goto_p(p0)
                    .goto_m(0.0)
                    .compute_cut_path_p();
                self.create_cut(Component::P, CutType::UShortKidney(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);
                self.create_cut(Component::P, CutType::UShortScallion(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);
            }
        }

        self.e_start(p_range);

        self.compute_cut_e_p()
            .create_cut(Component::P, CutType::E)
            .push_cut(p_range);

        self.compute_cut_e_xp();

        #[allow(clippy::comparison_chain)]
        if p_range == 0 {
            self.create_cut(Component::Xp, CutType::E)
                .log_branch(p_range)
                .xm_between()
                .push_cut(p_range);

            self.compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmPositive)
                .compute_cut_path_x(CutDirection::Negative)
                .split_cut(
                    Component::Xm,
                    CutType::UShortScallion(Component::Xp),
                    SplitCutBranchPoint::Old,
                    SplitCutOrder::NewFirst,
                );

            self.create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);
        } else if p_range < 0 {
            self.create_cut(Component::Xp, CutType::E)
                .log_branch(p_range)
                .xm_between()
                .push_cut(p_range);
        } else {
            self.create_cut(Component::Xp, CutType::E)
                .log_branch(p_range)
                .xm_outside()
                .push_cut(p_range);
        }

        self.compute_cut_e_xm();

        if p_range == 0 {
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range)
                .xp_inside()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Positive);

            self.split_cut(
                Component::Xp,
                CutType::UShortKidney(Component::Xp),
                SplitCutBranchPoint::New,
                SplitCutOrder::NewFirst,
            );

            self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);
        } else if p_range == -1 {
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range)
                .xp_outside()
                .push_cut(p_range);

            self.compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmNegative)
                .compute_cut_path_x(CutDirection::Negative)
                .split_cut(
                    Component::Xp,
                    CutType::UShortScallion(Component::Xp),
                    SplitCutBranchPoint::Old,
                    SplitCutOrder::NewFirst,
                );

            self.create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);
        } else if p_range < 0 {
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range)
                .xp_outside()
                .push_cut(p_range);
        } else {
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range)
                .xp_between()
                .xm_outside()
                .push_cut(p_range);
        }

        if p_range == 0 {
            self.compute_cut_e_u()
                .create_cut(Component::U, CutType::E)
                .log_branch(p_range)
                .xp_between()
                .xm_between()
                .push_cut(p_range);

            self.set_cut_path(
                vec![
                    Complex64::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                    Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                ],
                Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
            )
            .split_cut(
                Component::Xm,
                CutType::UShortScallion(Component::Xp),
                SplitCutBranchPoint::Old,
                SplitCutOrder::NewFirst,
            );

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range - 1) as f64 * k) / consts.h,
                    )),
                )
                .split_cut(
                    Component::Xp,
                    CutType::UShortKidney(Component::Xm),
                    SplitCutBranchPoint::New,
                    SplitCutOrder::NewFirst,
                );

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                        Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
                )
                .create_cut(Component::U, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range - 1) as f64 * k) / consts.h,
                    )),
                )
                .create_cut(Component::U, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);
        } else if p_range == -1 {
            self.compute_cut_e_u()
                .create_cut(Component::U, CutType::E)
                .log_branch(p_range)
                .xp_outside()
                .xm_between()
                .push_cut(p_range);

            self.set_cut_path(
                vec![
                    Complex64::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                    Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                ],
                Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
            )
            .split_cut(
                Component::Xp,
                CutType::UShortScallion(Component::Xp),
                SplitCutBranchPoint::Old,
                SplitCutOrder::NewFirst,
            );

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range + 1) as f64 * k) / consts.h,
                    )),
                )
                .split_cut(
                    Component::Xm,
                    CutType::UShortKidney(Component::Xp),
                    SplitCutBranchPoint::New,
                    SplitCutOrder::NewFirst,
                );

            self.create_cut(Component::U, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_between()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                        Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
                )
                .create_cut(Component::U, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);
        } else if p_range > 0 {
            self.compute_cut_e_u()
                .create_cut(Component::U, CutType::E)
                .log_branch(p_range)
                .xp_between()
                .xm_outside()
                .push_cut(p_range);
        } else if p_range < -1 {
            self.compute_cut_e_u()
                .create_cut(Component::U, CutType::E)
                .log_branch(p_range)
                .xp_outside()
                .xm_between()
                .push_cut(p_range);
        }
    }

    fn generate_short_cuts(&mut self, p_range: i32, consts: CouplingConstants) {
        let p_start = p_range as f64;
        let k = consts.k() as f64;
        let s = consts.s();

        let us = s + 1.0 / s - (s - 1.0 / s) * s.ln();

        self.begin_group(&[CutVisibilityCondition::ShortCuts]);

        {
            // Log

            self.clear_cut()
                .set_cut_path(
                    vec![Complex64::from(-INFINITY), Complex64::from(0.0)],
                    Some(Complex64::from(0.0)),
                )
                .create_cut(Component::Xp, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Negative);
            self.create_cut(Component::Xm, CutType::Log(Component::Xp))
                .xp_between()
                .push_cut(p_range);
            self.create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .xp_inside()
                .im_xp_positive()
                .im_xm_negative()
                .push_cut(p_range);
            self.create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range + 1)
                .xp_inside()
                .im_xp_negative()
                .im_xm_negative()
                .push_cut(p_range);
            self.create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range + 1)
                .xp_inside()
                .im_xm_positive()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Negative);
            self.create_cut(Component::Xm, CutType::Log(Component::Xp))
                .xp_between()
                .push_cut(p_range);
            self.create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .xp_inside()
                .im_xp_positive()
                .push_cut(p_range);
            self.create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .xp_inside()
                .im_xp_negative()
                .im_xm_negative()
                .push_cut(p_range);
            self.create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range + 1)
                .xp_inside()
                .im_xp_negative()
                .im_xm_positive()
                .push_cut(p_range);

            self.clear_cut().set_cut_path(
                vec![
                    Complex64::new(-INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                    Complex64::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                ],
                None,
            );

            self.create_cut(Component::U, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .xp_inside()
                .im_xp_positive()
                .push_cut(p_range);

            if p_range == 0 {
                self.create_cut(Component::U, CutType::Log(Component::Xp))
                    .xp_between()
                    .periodic()
                    .push_cut(p_range);
            }

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(-INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                    ],
                    None,
                )
                .create_cut(Component::U, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .xp_inside()
                .im_xp_negative()
                .push_cut(p_range);

            // Real negative line
            let p0 = p_start + 7.0 / 8.0;
            self.clear_cut()
                .p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .goto_re(-consts.s() * 4.0)
                .compute_cut_path_p()
                .create_cut(Component::P, CutType::Log(Component::Xp))
                .e_branch(1)
                .push_cut(p_range);
        }

        {
            // U long positive
            self.clear_cut()
                .set_cut_path(
                    vec![Complex64::from(0.0), Complex64::from(INFINITY)],
                    Some(Complex64::from(s)),
                )
                .create_cut(Component::Xp, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmNegative)
                .compute_cut_path_x(CutDirection::Positive);
            self.create_cut(Component::Xm, CutType::ULongPositive(Component::Xp))
                .xp_between()
                .push_cut(p_range);
            self.create_cut(Component::Xm, CutType::ULongPositive(Component::Xp))
                .xp_outside()
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmPositive)
                .compute_cut_path_x(CutDirection::Positive);
            self.create_cut(Component::Xm, CutType::ULongPositive(Component::Xp))
                .xp_between()
                .push_cut(p_range);

            if p_range <= 0 {
                self.create_cut(Component::Xm, CutType::ULongPositive(Component::Xp))
                    .log_branch(p_range)
                    .xp_outside()
                    .push_cut(p_range);
            } else {
                self.create_cut(Component::Xm, CutType::ULongPositive(Component::Xp))
                    .log_branch(p_range)
                    .xp_outside()
                    .im_xm_positive()
                    .push_cut(p_range);

                self.create_cut(Component::Xm, CutType::ULongPositive(Component::Xp))
                    .log_branch(p_range - 1)
                    .xp_outside()
                    .im_xm_negative()
                    .push_cut(p_range);
            }

            self.clear_cut().set_cut_path(
                vec![
                    Complex64::new(INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                    Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                ],
                None,
            );

            self.create_cut(Component::U, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .xp_outside()
                .push_cut(p_range);

            if p_range == 0 {
                self.create_cut(Component::U, CutType::ULongPositive(Component::Xp))
                    .xp_between()
                    .periodic()
                    .push_cut(p_range);
            }

            // Real positive line
            let p0 = p_start + 1.0 / 8.0;
            self.clear_cut()
                .p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .goto_re(consts.s() * 4.0)
                .compute_cut_path_p()
                .create_cut(Component::P, CutType::ULongPositive(Component::Xp))
                .e_branch(1)
                .push_cut(p_range);

            self.create_cut(Component::P, CutType::Log(Component::Xm))
                .e_branch(-1)
                .push_cut(p_range);

            // Real negative line
            let p0 = p_start + 7.0 / 8.0;
            self.clear_cut()
                .p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .goto_re(-consts.s() * 4.0)
                .compute_cut_path_p();

            self.create_cut(Component::P, CutType::ULongPositive(Component::Xm))
                .e_branch(-1)
                .push_cut(p_range);
        }

        {
            // Scallion
            self.clear_cut()
                .compute_cut_path_x_full(XCut::Scallion)
                .create_cut(Component::Xp, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmPositive)
                .compute_cut_path_x(CutDirection::Negative);
            self.create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                .xp_between()
                .push_cut(p_range);
            if p_range <= 0 {
                self.create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                    .log_branch(p_range)
                    .xp_outside()
                    .push_cut(p_range);
            } else {
                self.create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                    .log_branch(p_range)
                    .xp_outside()
                    .im_xm_positive()
                    .push_cut(p_range);
                self.create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                    .log_branch(p_range - 1)
                    .xp_outside()
                    .im_xm_negative()
                    .push_cut(p_range);
            }

            self.clear_cut()
                .compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmNegative)
                .compute_cut_path_x(CutDirection::Negative);
            self.create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                .xp_between()
                .push_cut(p_range);
            if p_range >= 0 {
                self.create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                    .xp_outside()
                    .log_branch(p_range)
                    .push_cut(p_range);
            } else {
                self.create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                    .log_branch(p_range)
                    .xp_outside()
                    .im_xm_negative()
                    .push_cut(p_range);

                self.create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                    .log_branch(p_range + 1)
                    .xp_outside()
                    .im_xm_positive()
                    .push_cut(p_range);
            }

            self.clear_cut().set_cut_path(
                vec![
                    Complex64::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                    Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                ],
                Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
            );

            self.create_cut(Component::U, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .xp_outside()
                .push_cut(p_range);

            if p_range == 0 {
                self.create_cut(Component::U, CutType::UShortScallion(Component::Xp))
                    .xp_between()
                    .periodic()
                    .push_cut(p_range);
            }

            let p0 = p_start + 1.0 / 8.0;
            self.clear_cut();

            self.p_start_xp(p0)
                .goto_m(-(p_range * consts.k()) as f64)
                .compute_cut_path_p_rev();

            self.p_start_xp(p0)
                .goto_xm(p0, 1.0)
                .goto_m(-(p_range * consts.k()) as f64)
                .compute_cut_path_p();

            self.create_cut(Component::P, CutType::UShortScallion(Component::Xp))
                .e_branch(1)
                .push_cut(p_range);

            self.create_cut(Component::P, CutType::UShortKidney(Component::Xm))
                .e_branch(-1)
                .push_cut(p_range);

            if p_range == 0 {
                let p0 = p_start + 1.0 / 8.0;
                let p1 = p_start + 7.0 / 8.0;

                self.clear_cut()
                    .p_start_xp(p0)
                    .goto_m(3.0)
                    .goto_p(p1)
                    .goto_m(0.0)
                    .compute_cut_path_p();

                self.create_cut(Component::P, CutType::UShortScallion(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);

                self.create_cut(Component::P, CutType::UShortKidney(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);
            }
        }

        {
            // Kidney

            self.clear_cut()
                .compute_cut_path_x_full(XCut::Kidney)
                .create_cut(Component::Xp, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Positive);

            self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .xp_between()
                .push_cut(p_range);

            if p_range >= -1 {
                self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                    .log_branch(p_range)
                    .xp_inside()
                    .im_xp_positive()
                    .push_cut(p_range);

                self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                    .log_branch(p_range + 1)
                    .xp_inside()
                    .im_xp_negative()
                    .push_cut(p_range);
            } else {
                self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                    .log_branch(p_range)
                    .xp_inside()
                    .im_xp_positive()
                    .im_xm_negative()
                    .push_cut(p_range);

                self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                    .log_branch(p_range + 1)
                    .xp_inside()
                    .im_xm_positive()
                    .push_cut(p_range);

                self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                    .log_branch(p_range + 1)
                    .xp_inside()
                    .im_xp_negative()
                    .im_xm_negative()
                    .push_cut(p_range);
            }

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Positive);
            self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .xp_between()
                .push_cut(p_range);

            if true {
                if p_range == -1 {
                    self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                        .log_branch(p_range)
                        .xp_inside()
                        .im_xp_positive()
                        .push_cut(p_range);

                    self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                        .log_branch(p_range + 1)
                        .xp_inside()
                        .im_xp_negative()
                        .push_cut(p_range);
                } else if p_range >= 0 {
                    self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                        .log_branch(p_range)
                        .xp_inside()
                        .im_xp_positive()
                        .im_xm_positive()
                        .push_cut(p_range);

                    self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                        .log_branch(p_range)
                        .xp_inside()
                        .im_xp_negative()
                        .im_xm_negative()
                        .push_cut(p_range);

                    self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                        .log_branch(p_range - 1)
                        .xp_inside()
                        .im_xp_positive()
                        .im_xm_negative()
                        .push_cut(p_range);

                    self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                        .log_branch(p_range + 1)
                        .xp_inside()
                        .im_xp_negative()
                        .im_xm_positive()
                        .push_cut(p_range);
                } else {
                    self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                        .log_branch(p_range)
                        .xp_inside()
                        .im_xp_positive()
                        .push_cut(p_range);

                    self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                        .log_branch(p_range)
                        .xp_inside()
                        .im_xp_negative()
                        .im_xm_negative()
                        .push_cut(p_range);

                    self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                        .log_branch(p_range + 1)
                        .xp_inside()
                        .im_xp_negative()
                        .im_xm_positive()
                        .push_cut(p_range);
                }
            }

            self.clear_cut().set_cut_path(
                vec![
                    Complex64::new(INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                    Complex64::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                ],
                Some(Complex64::new(
                    -us,
                    -(1.0 + (p_range + 1) as f64 * k) / consts.h,
                )),
            );

            self.create_cut(Component::U, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_inside()
                .im_xp_positive()
                .push_cut(p_range);

            if p_range == 0 {
                self.create_cut(Component::U, CutType::UShortKidney(Component::Xp))
                    .xp_between()
                    .periodic()
                    .push_cut(p_range);
            }

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range - 1) as f64 * k) / consts.h,
                    )),
                )
                .create_cut(Component::U, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_inside()
                .im_xp_negative()
                .push_cut(p_range);

            let p0 = p_start + 1.0 / 8.0;
            let p1 = p_start + 7.0 / 8.0;

            if p_range != 0 {
                self.clear_cut();
                self.p_start_xp(p0)
                    .goto_m(-((p_range + 1) * consts.k()) as f64)
                    .compute_cut_path_p_rev();

                self.p_start_xp(p0)
                    .goto_xm(p0, 1.0)
                    .goto_p(p1)
                    .goto_m(-((p_range + 1) * consts.k()) as f64)
                    .compute_cut_path_p();

                self.create_cut(Component::P, CutType::UShortKidney(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);
                self.create_cut(Component::P, CutType::UShortScallion(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);
            } else {
                self.clear_cut();
                self.p_start_xp(p0)
                    .goto_m(-((p_range + 1) * consts.k()) as f64)
                    .compute_cut_path_p_rev();

                self.create_cut(Component::P, CutType::UShortKidney(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);
                self.create_cut(Component::P, CutType::UShortScallion(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);

                self.clear_cut();

                self.p_start_xp(p0)
                    .goto_xm(p0, 1.0)
                    .goto_p(p1)
                    .goto_m(-((p_range + 1) * consts.k()) as f64)
                    .compute_cut_path_p();

                self.create_cut(Component::P, CutType::UShortKidney(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);
                self.create_cut(Component::P, CutType::UShortScallion(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);
            }

            if p_range == -1 {
                let p0 = p_start + 1.0 / 8.0;
                let p1 = p_start + 6.0 / 8.0;

                self.clear_cut()
                    .p_start_xp(p0)
                    .goto_m(-consts.k() as f64 + 1.0)
                    .goto_p(p1)
                    .goto_m(-consts.k() as f64 + 3.0)
                    .goto_p(p0)
                    .goto_m(0.0)
                    .compute_cut_path_p();
                self.create_cut(Component::P, CutType::UShortKidney(Component::Xp))
                    .e_branch(1)
                    .push_cut(p_range);
                self.create_cut(Component::P, CutType::UShortScallion(Component::Xm))
                    .e_branch(-1)
                    .push_cut(p_range);
            }
        }

        self.e_start(p_range);

        self.compute_cut_e_p()
            .create_cut(Component::P, CutType::E)
            .push_cut(p_range);

        self.compute_cut_e_xp();

        #[allow(clippy::comparison_chain)]
        if p_range == 0 {
            self.create_cut(Component::Xp, CutType::E)
                .log_branch(p_range)
                .xm_between()
                .push_cut(p_range);

            self.compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmPositive)
                .compute_cut_path_x(CutDirection::Negative)
                .split_cut(
                    Component::Xm,
                    CutType::UShortScallion(Component::Xp),
                    SplitCutBranchPoint::Old,
                    SplitCutOrder::NewFirst,
                );

            self.pop_cut();
            self.create_cut(Component::Xp, CutType::E)
                .xm_between()
                .push_cut(p_range);
        }
        if p_range == -1 {
            self.create_cut(Component::Xp, CutType::E)
                .xm_between()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Positive)
                .split_cut(
                    Component::Xm,
                    CutType::UShortKidney(Component::Xm),
                    SplitCutBranchPoint::New,
                    SplitCutOrder::OldFirst,
                );
            self.pop_cut();
            self.create_cut(Component::Xp, CutType::E)
                .log_branch(p_range)
                .xm_inside()
                .im_xm_negative()
                .push_cut(p_range);
            self.create_cut(Component::Xp, CutType::E)
                .log_branch(p_range + 1)
                .xm_inside()
                .im_xm_positive()
                .push_cut(p_range);
        }
        if p_range < -1 {
            self.create_cut(Component::Xp, CutType::E)
                .log_branch(p_range)
                .xm_inside()
                .im_xm_negative()
                .push_cut(p_range);

            self.create_cut(Component::Xp, CutType::E)
                .log_branch(p_range + 1)
                .xm_inside()
                .im_xm_positive()
                .push_cut(p_range);
        }
        if p_range > 0 {
            self.create_cut(Component::Xp, CutType::E)
                .log_branch(p_range)
                .xm_outside()
                .im_xp_negative()
                .push_cut(p_range);

            self.create_cut(Component::Xp, CutType::E)
                .log_branch(p_range - 1)
                .xm_outside()
                .im_xp_positive()
                .push_cut(p_range);
        }

        self.compute_cut_e_xm();

        if p_range == 0 {
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range)
                .xp_inside()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromBelowWithImXmNegative,
                )
                .compute_cut_path_x(CutDirection::Positive);

            self.split_cut(
                Component::Xp,
                CutType::UShortKidney(Component::Xp),
                SplitCutBranchPoint::New,
                SplitCutOrder::OldFirst,
            );

            self.pop_cut();
            self.create_cut(Component::Xm, CutType::E)
                .xp_between()
                .push_cut(p_range);

            self.swap_cuts();
            self.pop_cut();
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range)
                .xp_inside()
                .im_xp_negative()
                .push_cut(p_range);
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range - 1)
                .xp_inside()
                .im_xp_positive()
                .push_cut(p_range);
        }
        if p_range == -1 {
            self.create_cut(Component::Xm, CutType::DebugPath)
                .push_cut(p_range);

            self.compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmNegative)
                .compute_cut_path_x(CutDirection::Negative)
                .split_cut(
                    Component::Xp,
                    CutType::UShortScallion(Component::Xp),
                    SplitCutBranchPoint::Old,
                    SplitCutOrder::OldFirst,
                );

            self.pop_cut();
            self.create_cut(Component::Xm, CutType::E)
                .xp_between()
                .push_cut(p_range);
            self.swap_cuts();
            self.pop_cut();
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range)
                .xp_outside()
                .im_xm_negative()
                .push_cut(p_range);
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range + 1)
                .xp_outside()
                .im_xm_positive()
                .push_cut(p_range);
        }
        if p_range < -1 {
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range)
                .xp_outside()
                .im_xm_negative()
                .push_cut(p_range);

            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range + 1)
                .xp_outside()
                .im_xm_positive()
                .push_cut(p_range);
        }
        if p_range > 0 {
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range)
                .xp_inside()
                .im_xp_negative()
                .push_cut(p_range);

            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range - 1)
                .xp_inside()
                .im_xp_positive()
                .push_cut(p_range);
        }

        if p_range == 0 {
            self.compute_cut_e_u()
                .create_cut(Component::U, CutType::DebugPath)
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range - 1) as f64 * k) / consts.h,
                    )),
                )
                .split_cut(
                    Component::Xp,
                    CutType::UShortKidney(Component::Xp),
                    SplitCutBranchPoint::New,
                    SplitCutOrder::NewFirst,
                );

            self.pop_cut();

            for p_range in 0..=P_RANGE_MAX {
                self.create_cut(Component::U, CutType::E)
                    .log_branch(p_range)
                    .xp_inside()
                    .xm_between()
                    .im_xp_negative()
                    .pre_shift(-p_range as f64 * k * Complex64::i() / consts.h)
                    .push_cut(p_range);
                self.swap_cuts();
            }

            for p_range in P_RANGE_MIN..=0 {
                self.create_cut(Component::U, CutType::E)
                    .log_branch(p_range)
                    .xp_inside()
                    .xm_between()
                    .im_xp_positive()
                    .pre_shift(-(2 + p_range) as f64 * k * Complex64::i() / consts.h)
                    .push_cut(p_range);
                self.swap_cuts();
            }

            self.set_cut_path(
                vec![
                    Complex64::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                    Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                ],
                Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
            )
            .split_cut(
                Component::Xm,
                CutType::UShortScallion(Component::Xm),
                SplitCutBranchPoint::Old,
                SplitCutOrder::NewFirst,
            );

            self.pop_cut();

            self.create_cut(Component::U, CutType::E)
                .xp_between()
                .xm_between()
                .periodic()
                .push_cut(p_range);
            self.swap_cuts();

            self.pop_cut();

            for p_range in P_RANGE_MIN..=P_RANGE_MAX {
                self.create_cut(Component::U, CutType::E)
                    .log_branch(p_range)
                    .xp_between()
                    .xm_outside()
                    .pre_shift(p_range as f64 * k * Complex64::i() / consts.h)
                    .push_cut(p_range);
            }
        } else if p_range == -1 {
            self.compute_cut_e_u()
                .create_cut(Component::U, CutType::E)
                .xp_outside()
                .xm_between()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        Complex64::new(INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                        Complex64::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                    ],
                    Some(Complex64::new(
                        -us,
                        -(1.0 + (p_range + 1) as f64 * k) / consts.h,
                    )),
                )
                .split_cut(
                    Component::Xm,
                    CutType::UShortKidney(Component::Xp),
                    SplitCutBranchPoint::New,
                    SplitCutOrder::NewFirst,
                );

            self.pop_cut();

            for p_range in P_RANGE_MIN..=P_RANGE_MAX {
                self.create_cut(Component::U, CutType::E)
                    .log_branch(p_range)
                    .xp_outside()
                    .xm_between()
                    .pre_shift(-(p_range + 1) as f64 * k * Complex64::i() / consts.h)
                    .push_cut(p_range);
                self.swap_cuts();
            }

            self.set_cut_path(
                vec![
                    Complex64::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                    Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                ],
                Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
            )
            .split_cut(
                Component::Xp,
                CutType::UShortScallion(Component::Xp),
                SplitCutBranchPoint::Old,
                SplitCutOrder::OldFirst,
            );
            self.pop_cut();

            for p_range in P_RANGE_MIN..=-1 {
                self.create_cut(Component::U, CutType::E)
                    .log_branch(p_range)
                    .xp_between()
                    .xm_inside()
                    .pre_shift((p_range + 1) as f64 * k * Complex64::i() / consts.h)
                    .push_cut(p_range);
                self.swap_cuts();
            }

            for p_range in 0..=P_RANGE_MAX {
                self.create_cut(Component::U, CutType::E)
                    .log_branch(p_range)
                    .xp_between()
                    .xm_inside()
                    .pre_shift((p_range - 1) as f64 * k * Complex64::i() / consts.h)
                    .push_cut(p_range);
                self.swap_cuts();
            }

            self.pop_cut();

            self.create_cut(Component::U, CutType::E)
                .log_branch(p_range)
                .xp_outside()
                .xm_inside()
                .im_xm_negative()
                .push_cut(p_range);

            self.create_cut(Component::U, CutType::E)
                .log_branch(p_range + 1)
                .xp_outside()
                .xm_inside()
                .im_xm_positive()
                .pre_shift(-1.0 * k * Complex64::i() / consts.h)
                .push_cut(p_range + 1);
        } else if p_range > 0 {
            self.compute_cut_e_u();
            self.create_cut(Component::U, CutType::E)
                .log_branch(p_range)
                .xp_inside()
                .xm_outside()
                .im_xp_negative()
                .push_cut(p_range);
            self.create_cut(Component::U, CutType::E)
                .log_branch(p_range - 1)
                .xp_inside()
                .xm_outside()
                .im_xp_positive()
                .pre_shift(-1.0 * k * Complex64::i() / consts.h)
                .push_cut(p_range - 1);
        } else if p_range < -1 {
            self.compute_cut_e_u();
            self.create_cut(Component::U, CutType::E)
                .log_branch(p_range)
                .xp_outside()
                .xm_inside()
                .im_xm_negative()
                .push_cut(p_range);
            self.compute_cut_e_u();
            self.create_cut(Component::U, CutType::E)
                .log_branch(p_range + 1)
                .xp_outside()
                .xm_inside()
                .im_xm_positive()
                .pre_shift(-1.0 * k * Complex64::i() / consts.h)
                .push_cut(p_range + 1);
        }
    }
}
