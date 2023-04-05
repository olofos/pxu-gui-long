use std::collections::VecDeque;

use crate::interpolation::{EPInterpolator, InterpolationPoint, PInterpolatorMut, XInterpolator};
use crate::kinematics::{
    du_crossed_dp, du_dp, dxm_crossed_dp, dxm_dp, dxp_crossed_dp, dxp_dp, u, u_crossed, xm,
    xm_crossed, xp, xp_crossed, CouplingConstants, SheetData, UBranch,
};
use crate::nr::{self};
use itertools::Itertools;
use num::complex::Complex64;

const P_RANGE_MIN: i32 = -3;
const P_RANGE_MAX: i32 = 3;

const INFINITY: f64 = 100.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Component {
    P,
    Xp,
    Xm,
    U,
}

impl Component {
    fn conj(&self) -> Self {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum UCutType {
    Long,
    SemiShort,
    #[default]
    Short,
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
enum GeneratorCommands {
    AddGridLineU(f64),
    AddGridLineX(f64),
    AddGridLineP,

    ComputeBranchPoint(i32, BranchPointType),
    ClearCut,
    ComputeCutX(CutDirection),
    ComputeCutXFull(XCut),
    ComputeCutP(bool),
    ComputeCutEP,
    ComputeCutEXp,
    ComputeCutEXm,
    ComputeCutEU,
    SetCutPath(Vec<Complex64>, Option<Complex64>),
    PushCut(
        i32,
        Component,
        CutType,
        Option<Complex64>,
        Vec<CutVisibilityCondition>,
        Complex64,
    ),
    SplitCut(i32, Component, CutType, SplitCutBranchPoint, SplitCutOrder),

    EStart(i32),

    PStartXp(f64),
    PGotoXp(f64, f64),
    PGotoXm(f64, f64),
    PGotoP(f64),
    PGotoM(f64),
    PGotoIm(f64),
    PGotoRe(f64),
}

struct RuntimeCutData {
    branch_point: Option<Complex64>,
    path: Option<Vec<Complex64>>,
}

struct BuildTimeCutData {
    component: Option<Component>,
    cut_type: Option<CutType>,
    period: Option<Complex64>,
    visibility: Vec<CutVisibilityCondition>,
    group_visibility: Vec<CutVisibilityCondition>,
    pre_shift: Complex64,
}

struct ContourGeneratorRuntimeContext {
    p_int: Option<PInterpolatorMut>,
    e_int: Option<EPInterpolator>,
    branch_point_data: Option<BranchPointData>,
    cut_data: RuntimeCutData,
}

impl ContourGeneratorRuntimeContext {
    fn new() -> Self {
        Self {
            p_int: None,
            e_int: None,
            branch_point_data: None,
            cut_data: RuntimeCutData {
                branch_point: None,
                path: None,
            },
        }
    }
}

struct ContourGeneratorBuildTimeContext {
    cut_data: BuildTimeCutData,
}

impl ContourGeneratorBuildTimeContext {
    fn new() -> Self {
        Self {
            cut_data: BuildTimeCutData {
                component: None,
                cut_type: None,
                period: None,
                visibility: vec![],
                group_visibility: vec![],
                pre_shift: Complex64::from(0.0),
            },
        }
    }

    fn reset(&mut self) {
        self.cut_data.component = None;
        self.cut_data.cut_type = None;
        self.cut_data.period = None;
        self.cut_data.visibility = self.cut_data.group_visibility.clone();
        self.cut_data.pre_shift = Complex64::from(0.0);
    }
}

pub enum GridLineComponent {
    Real,
    Xp(f64),
    Xm(f64),
}

pub struct GridLine {
    pub path: Vec<Complex64>,
    pub component: GridLineComponent,
}

pub struct ContourGenerator {
    cuts: Vec<Cut>,
    commands: VecDeque<GeneratorCommands>,
    pub consts: Option<CouplingConstants>,

    grid_p: Vec<GridLine>,
    grid_x: Vec<GridLine>,
    grid_u: Vec<GridLine>,

    rctx: ContourGeneratorRuntimeContext,
    bctx: ContourGeneratorBuildTimeContext,

    num_commands: usize,
}

impl Default for ContourGenerator {
    fn default() -> Self {
        Self {
            cuts: vec![],
            commands: VecDeque::new(),
            consts: None,
            grid_p: vec![],
            grid_x: vec![],
            grid_u: vec![],
            rctx: ContourGeneratorRuntimeContext::new(),
            bctx: ContourGeneratorBuildTimeContext::new(),
            num_commands: 0,
        }
    }
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

impl ContourGenerator {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn generate_all(consts: CouplingConstants) -> Self {
        let pt = Point::new(0.5, consts);
        let mut generator = Self::new();
        while !generator.update(&pt) {}
        generator
    }

    pub fn update(&mut self, pt: &Point) -> bool {
        if let Some(consts) = self.consts {
            if consts != pt.consts {
                self.consts = None;
            }
        }

        if self.consts.is_none() {
            self.clear();
            self.consts = Some(pt.consts);
            self.generate_commands(pt);
            self.num_commands = self.commands.len();
            log::debug!("Generated {} commands", self.num_commands,)
        }

        if let Some(command) = self.commands.pop_front() {
            self.execute(command);
        }

        self.commands.is_empty()
    }

    fn clear(&mut self) {
        log::debug!("Clearing grid and cuts");
        self.commands.clear();
        self.grid_x.clear();
        self.grid_u.clear();
        self.cuts.clear();

        self.grid_p = vec![GridLine {
            path: vec![
                Complex64::from(P_RANGE_MIN as f64),
                Complex64::from(P_RANGE_MAX as f64 + 1.0),
            ],
            component: GridLineComponent::Real,
        }];
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
        pt: &Point,
        component: Component,
        u_cut_type: UCutType,
    ) -> impl Iterator<Item = &Cut> {
        let mut pt = pt.clone();
        pt.u += 2.0 * (pt.sheet_data.log_branch_p * pt.consts.k()) as f64 * Complex64::i()
            / pt.consts.h;

        self.cuts
            .iter()
            .filter(move |c| c.component == component && c.is_visible(&pt, u_cut_type))
    }

    pub fn get_crossed_cuts(
        &self,
        pt: &Point,
        component: Component,
        new_value: Complex64,
        u_cut_type: UCutType,
    ) -> impl Iterator<Item = &Cut> {
        let mut pt = pt.clone();
        pt.u += 2.0 * (pt.sheet_data.log_branch_p * pt.consts.k()) as f64 * Complex64::i()
            / pt.consts.h;

        let new_value = if component == Component::U {
            new_value
                + 2.0 * (pt.sheet_data.log_branch_p * pt.consts.k()) as f64 * Complex64::i()
                    / pt.consts.h
        } else {
            new_value
        };

        self.cuts.iter().filter(move |c| {
            c.component == component
                && c.is_visible(&pt, u_cut_type)
                && c.intersection(pt.get(component), new_value).is_some()
        })
    }

    fn execute(&mut self, command: GeneratorCommands) {
        use GeneratorCommands::*;

        let Some(consts) = self.consts else {
            log::warn!("Executing commands when consts is not set!");
            return;
        };

        match command {
            AddGridLineU(y) => {
                self.grid_u.push(GridLine {
                    path: vec![Complex64::new(-INFINITY, y), Complex64::new(INFINITY, y)],
                    component: GridLineComponent::Real,
                });
            }

            AddGridLineX(m) => {
                let path = XInterpolator::generate_xp_full(0, m, consts);
                self.grid_x.push(GridLine {
                    path: path.iter().map(|x| x.conj()).collect(),
                    component: GridLineComponent::Xm(m),
                });
                self.grid_x.push(GridLine {
                    path,
                    component: GridLineComponent::Xp(m),
                });
            }

            EStart(p_range) => {
                self.rctx.e_int = Some(EPInterpolator::new(p_range, consts));
            }

            PStartXp(p) => {
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

                self.grid_p.push(GridLine {
                    path: path.iter().map(|p| p.conj()).collect(),
                    component: conj_component,
                });
                self.grid_p.push(GridLine { path, component });
            }

            ClearCut => {
                self.rctx.cut_data.path = None;
                self.rctx.cut_data.branch_point = None;
            }

            ComputeCutP(reverse) => {
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

            ComputeBranchPoint(p_range, branch_point_type) => {
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

            SetCutPath(path, branchpoint) => {
                self.rctx.cut_data.path = Some(path);
                self.rctx.cut_data.branch_point = branchpoint;
            }

            PushCut(p_range, component, cut_type, period, visibility, pre_shift) => {
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

                let paths = vec![path.iter().map(|z| z + pre_shift).collect()];

                let cut = Cut::new(
                    component,
                    paths,
                    self.rctx.cut_data.branch_point.map(|z| z + pre_shift),
                    cut_type,
                    p_range,
                    period,
                    visibility,
                );

                self.cuts.push(cut.conj().shift(shift));
                self.cuts.push(cut.shift(shift));
            }

            SplitCut(p_range, component, cut_typ, branch_point, order) => {
                let Some(mut cut) = self.cuts.pop() else {return};
                let Some(_) = self.cuts.pop() else {return};
                let Some(ref path) = self.rctx.cut_data.path else { return };

                let shift = match cut.component {
                    Component::U => Complex64::new(0.0, (p_range * consts.k()) as f64 / consts.h),
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
                    if let Some((i, j, x)) = cut.intersection(p1, p2) {
                        let mut new_path = vec![x];
                        new_path.extend(cut.paths[i].split_off(j + 1));
                        cut.paths[i].push(x);

                        let mut new_cut = Cut {
                            paths: vec![new_path],
                            branch_point: None,
                            typ: cut.typ.clone(),
                            p_range,
                            component: cut.component,
                            period: None,
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

                log::debug!("No intersection found");

                self.cuts.push(cut.conj());
                self.cuts.push(cut);
            }
        }
    }

    fn add(&mut self, command: GeneratorCommands) -> &mut Self {
        self.commands.push_back(command);
        self
    }

    fn e_start(&mut self, p_range: i32) -> &mut Self {
        self.add(GeneratorCommands::EStart(p_range))
    }

    fn compute_cut_e_p(&mut self) -> &mut Self {
        self.add(GeneratorCommands::ComputeCutEP)
    }

    fn compute_cut_e_xp(&mut self) -> &mut Self {
        self.add(GeneratorCommands::ComputeCutEXp)
    }

    fn compute_cut_e_xm(&mut self) -> &mut Self {
        self.add(GeneratorCommands::ComputeCutEXm)
    }

    fn compute_cut_e_u(&mut self) -> &mut Self {
        self.add(GeneratorCommands::ComputeCutEU)
    }

    fn p_start_xp(&mut self, p: f64) -> &mut Self {
        self.add(GeneratorCommands::PStartXp(p))
    }

    fn goto_xp(&mut self, p: f64, m: f64) -> &mut Self {
        self.add(GeneratorCommands::PGotoXp(p, m))
    }

    fn goto_xm(&mut self, p: f64, m: f64) -> &mut Self {
        self.add(GeneratorCommands::PGotoXm(p, m))
    }

    fn goto_re(&mut self, re: f64) -> &mut Self {
        self.add(GeneratorCommands::PGotoRe(re))
    }

    fn goto_im(&mut self, im: f64) -> &mut Self {
        self.add(GeneratorCommands::PGotoIm(im))
    }

    fn goto_p(&mut self, p: f64) -> &mut Self {
        self.add(GeneratorCommands::PGotoP(p))
    }

    fn goto_m(&mut self, m: f64) -> &mut Self {
        self.add(GeneratorCommands::PGotoM(m))
    }

    fn p_grid_line(&mut self) -> &mut Self {
        self.add(GeneratorCommands::AddGridLineP)
    }

    fn generate_commands(&mut self, pt: &Point) {
        let consts = pt.consts;
        self.generate_u_grid(consts);

        let p_range = pt.p.re.floor() as i32;

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
    }

    fn generate_u_grid(&mut self, consts: CouplingConstants) {
        self.add(GeneratorCommands::AddGridLineU(0.0));

        for y in 1..=100 {
            self.add(GeneratorCommands::AddGridLineU(y as f64 / consts.h));
            self.add(GeneratorCommands::AddGridLineU(-y as f64 / consts.h));
        }
    }

    fn generate_x_grid(&mut self, p_range: i32, _consts: CouplingConstants) {
        if p_range == 0 {
            for m in -50..50 {
                self.add(GeneratorCommands::AddGridLineX(m as f64));
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
        self.bctx.cut_data.group_visibility = Vec::from(group_visibility);
        self.bctx.reset();
        self
    }

    fn clear_cut(&mut self) -> &mut Self {
        self.add(GeneratorCommands::ClearCut)
    }

    fn compute_branch_point(
        &mut self,
        p_range: i32,
        branch_point_type: BranchPointType,
    ) -> &mut Self {
        self.add(GeneratorCommands::ComputeBranchPoint(
            p_range,
            branch_point_type,
        ))
    }

    fn compute_cut_path_x(&mut self, direction: CutDirection) -> &mut Self {
        self.add(GeneratorCommands::ComputeCutX(direction))
    }

    fn compute_cut_path_x_full(&mut self, xcut: XCut) -> &mut Self {
        self.add(GeneratorCommands::ComputeCutXFull(xcut))
    }

    fn compute_cut_path_p(&mut self) -> &mut Self {
        self.add(GeneratorCommands::ComputeCutP(false))
    }

    fn compute_cut_path_p_rev(&mut self) -> &mut Self {
        self.add(GeneratorCommands::ComputeCutP(true))
    }

    fn set_cut_path(&mut self, path: Vec<Complex64>, branchpoint: Option<Complex64>) -> &mut Self {
        self.add(GeneratorCommands::SetCutPath(path, branchpoint))
    }

    fn push_cut(&mut self, p_range: i32) -> &mut Self {
        let Some(component) = std::mem::replace(&mut self.bctx.cut_data.component, None) else {
            log::warn!("Can't push cut without component");
            self.bctx.reset();
            return self;
        };
        let Some(cut_type) = std::mem::replace(&mut self.bctx.cut_data.cut_type, None) else {
            log::warn!("Can't push cut without type");
            self.bctx.reset();
            return self;
        };
        self.add(GeneratorCommands::PushCut(
            p_range,
            component,
            cut_type,
            self.bctx.cut_data.period,
            self.bctx.cut_data.visibility.clone(),
            self.bctx.cut_data.pre_shift,
        ));
        self.bctx.reset();
        self
    }

    fn split_cut(
        &mut self,
        p_range: i32,
        component: Component,
        cut_type: CutType,
        branch_point: SplitCutBranchPoint,
        order: SplitCutOrder,
    ) -> &mut Self {
        self.add(GeneratorCommands::SplitCut(
            p_range,
            component,
            cut_type,
            branch_point,
            order,
        ))
    }

    fn create_cut(&mut self, component: Component, cut_type: CutType) -> &mut Self {
        if self.bctx.cut_data.component.is_some() || self.bctx.cut_data.cut_type.is_some() {
            log::warn!("New cut created before previous cut was pushed");
        }
        self.bctx.cut_data.component = Some(component);
        self.bctx.cut_data.cut_type = Some(cut_type);
        self.bctx.cut_data.visibility = self.bctx.cut_data.group_visibility.clone();
        self
    }

    fn log_branch(&mut self, p_range: i32) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::LogBranch(p_range));
        self
    }

    fn im_xm_negative(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::ImXm(-1));
        self
    }

    fn im_xm_positive(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::ImXm(1));
        self
    }

    fn im_xp_positive(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::ImXp(1));
        self
    }

    fn im_xp_negative(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::ImXp(-1));
        self
    }

    fn xp_outside(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::UpBranch(UBranch::Outside));
        self
    }

    fn xp_between(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::UpBranch(UBranch::Between));
        self
    }

    fn xp_inside(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::UpBranch(UBranch::Inside));
        self
    }

    fn xm_outside(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::UmBranch(UBranch::Outside));
        self
    }

    fn xm_between(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::UmBranch(UBranch::Between));
        self
    }

    fn xm_inside(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::UmBranch(UBranch::Inside));
        self
    }

    fn e_branch(&mut self, branch: i32) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::EBranch(branch));
        self
    }

    fn period(&mut self, z: Complex64) -> &mut Self {
        self.bctx.cut_data.period = Some(z);
        self
    }

    fn pre_shift(&mut self, z: Complex64) -> &mut Self {
        self.bctx.cut_data.pre_shift = z;
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

            // self.clear_cut()
            //     .set_cut_path(
            //         vec![
            //             Complex64::new(-INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
            //             Complex64::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
            //         ],
            //         Some(Complex64::new(
            //             -us,
            //             -(1.0 + (p_range + 1) as f64 * k) / consts.h,
            //         )),
            //     )
            //     .create_cut(Component::U, CutType::Log(Component::Xp))
            //     .log_branch(p_range)
            //     .xp_between()
            //     .push_cut(p_range);

            // self.clear_cut()
            //     .set_cut_path(
            //         vec![
            //             Complex64::new(-INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
            //             Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
            //         ],
            //         Some(Complex64::new(
            //             -us,
            //             -(1.0 + (p_range - 1) as f64 * k) / consts.h,
            //         )),
            //     )
            //     .create_cut(Component::U, CutType::Log(Component::Xp))
            //     .log_branch(p_range)
            //     .xp_between()
            //     .push_cut(p_range);

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
                    p_range,
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
                p_range,
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
                    p_range,
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
                p_range,
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
                    p_range,
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
                p_range,
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
                    p_range,
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

            {
                // self.clear_cut().set_cut_path(
                //     vec![
                //         Complex64::new(-INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                //         Complex64::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                //     ],
                //     Some(Complex64::new(
                //         -us,
                //         -(1.0 + (p_range - 1) as f64 * k) / consts.h,
                //     )),
                // );
                // self.create_cut(Component::U, CutType::Log(Component::Xp))
                //     .xp_between()
                //     .push_cut(p_range);
            }

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
                    .period(2.0 * Complex64::i() * k / consts.h)
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
                    .period(2.0 * Complex64::i() * k / consts.h)
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
                    // self.create_cut(Component::Xm, CutType::DebugPath)
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
            if p_range >= -1 {
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
                    .period(2.0 * Complex64::i() * k / consts.h)
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
            self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .xp_inside()
                .im_xp_positive()
                .im_xm_negative()
                .push_cut(p_range);
            self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range + 1)
                .xp_inside()
                .im_xp_negative()
                .im_xm_negative()
                .push_cut(p_range);
            self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range + 1)
                .xp_inside()
                .im_xm_positive()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPointType::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(CutDirection::Positive);
            self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .xp_between()
                .push_cut(p_range);
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
                    .period(2.0 * Complex64::i() * k / consts.h)
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
                    p_range,
                    Component::Xm,
                    CutType::UShortScallion(Component::Xp),
                    SplitCutBranchPoint::Old,
                    SplitCutOrder::NewFirst,
                );
        }
        if p_range == -1 {
            // self.create_cut(Component::Xp, CutType::E)
            //     .log_branch(p_range)
            //     // .xm_between()
            //     .push_cut(p_range);
        }
        if p_range < -1 {
            self.create_cut(Component::Xp, CutType::E)
                .log_branch(p_range)
                .xm_inside()
                .im_xm_negative()
                .push_cut(p_range);
        }
        if p_range <= -1 {
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
                p_range,
                Component::Xp,
                CutType::UShortKidney(Component::Xp),
                SplitCutBranchPoint::New,
                SplitCutOrder::NewFirst,
            );

            // self.create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
            //     .log_branch(p_range)
            //     .xp_between()
            //     .push_cut(p_range);
        }
        if p_range == -1 {
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range)
                .xp_outside()
                .push_cut(p_range);

            self.compute_branch_point(p_range, BranchPointType::XpPositiveAxisImXmNegative)
                .compute_cut_path_x(CutDirection::Negative)
                .split_cut(
                    p_range,
                    Component::Xp,
                    CutType::UShortScallion(Component::Xp),
                    SplitCutBranchPoint::Old,
                    SplitCutOrder::NewFirst,
                );

            // self.create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
            //     .log_branch(p_range)
            //     .push_cut(p_range);
        }
        if p_range < -1 {
            self.create_cut(Component::Xm, CutType::E)
                .log_branch(p_range)
                .xp_outside()
                .im_xm_negative()
                .push_cut(p_range);
        }
        if p_range < -1 {
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
                .create_cut(Component::U, CutType::E)
                // .log_branch(p_range)
                .xp_inside()
                .xm_between()
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
                    p_range,
                    Component::Xp,
                    CutType::UShortKidney(Component::Xm),
                    SplitCutBranchPoint::New,
                    SplitCutOrder::OldFirst,
                );

            self.set_cut_path(
                vec![
                    Complex64::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                    Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                ],
                Some(Complex64::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
            )
            .split_cut(
                p_range,
                Component::Xm,
                CutType::UShortScallion(Component::Xp),
                SplitCutBranchPoint::Old,
                SplitCutOrder::NewFirst,
            );
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
                p_range,
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
                    p_range,
                    Component::Xm,
                    CutType::UShortKidney(Component::Xp),
                    SplitCutBranchPoint::New,
                    SplitCutOrder::NewFirst,
                );
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
                .push_cut(p_range - 2);
        } else if p_range < -1 {
            self.compute_cut_e_u()
                .create_cut(Component::U, CutType::E)
                .log_branch(p_range)
                .xp_outside()
                .xm_inside()
                .push_cut(p_range);
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CutType {
    E,
    DebugPath,

    Log(Component),
    ULongPositive(Component),
    ULongNegative(Component),
    UShortScallion(Component),
    UShortKidney(Component),
}

impl CutType {
    fn conj(&self) -> Self {
        match self {
            Self::E => Self::E,
            Self::DebugPath => Self::DebugPath,

            Self::ULongPositive(component) => Self::ULongPositive(component.conj()),
            Self::ULongNegative(component) => Self::ULongNegative(component.conj()),
            Self::UShortScallion(component) => Self::UShortScallion(component.conj()),
            Self::UShortKidney(component) => Self::UShortKidney(component.conj()),
            Self::Log(component) => Self::Log(component.conj()),
        }
    }
}

#[derive(Debug, Clone)]
enum CutVisibilityCondition {
    ImXp(i8),
    ImXm(i8),
    LogBranch(i32),
    EBranch(i32),
    UpBranch(UBranch),
    UmBranch(UBranch),

    ImXpOrUpBranch(i8, UBranch),
    ImXmOrUmBranch(i8, UBranch),

    LongCuts,
    SemiShortCuts,
    ShortCuts,
    ShortOrSemiShortCuts,
}

impl CutVisibilityCondition {
    fn check(&self, pt: &Point, u_cut_type: UCutType) -> bool {
        let long_cuts = u_cut_type == UCutType::Long;
        match self {
            Self::ImXp(sign) => pt.xp.im.signum() as i8 == sign.signum(),
            Self::ImXm(sign) => pt.xm.im.signum() as i8 == sign.signum(),
            Self::LogBranch(b) => *b == (pt.sheet_data.log_branch_p + pt.sheet_data.log_branch_m),
            Self::EBranch(b) => pt.sheet_data.e_branch == *b,
            Self::UpBranch(b) => {
                if u_cut_type != UCutType::Short {
                    if *b == UBranch::Between || *b == UBranch::Inside {
                        pt.sheet_data.u_branch.0 != UBranch::Outside
                    } else {
                        pt.sheet_data.u_branch.0 == UBranch::Outside
                    }
                } else {
                    pt.sheet_data.u_branch.0 == *b
                }
            }
            Self::UmBranch(b) => {
                if u_cut_type != UCutType::Short {
                    if *b == UBranch::Between || *b == UBranch::Inside {
                        pt.sheet_data.u_branch.1 != UBranch::Outside
                    } else {
                        pt.sheet_data.u_branch.1 == UBranch::Outside
                    }
                } else {
                    pt.sheet_data.u_branch.1 == *b
                }
            }

            Self::ImXpOrUpBranch(b1, b2) => match u_cut_type {
                UCutType::Long => Self::ImXp(*b1).check(pt, u_cut_type),
                UCutType::SemiShort | UCutType::Short => {
                    Self::UpBranch(b2.clone()).check(pt, u_cut_type)
                }
            },

            Self::ImXmOrUmBranch(b1, b2) => match u_cut_type {
                UCutType::Long => Self::ImXm(*b1).check(pt, u_cut_type),
                UCutType::SemiShort | UCutType::Short => {
                    Self::UmBranch(b2.clone()).check(pt, u_cut_type)
                }
            },

            Self::LongCuts => long_cuts,
            Self::SemiShortCuts => u_cut_type == UCutType::SemiShort,
            Self::ShortCuts => u_cut_type == UCutType::Short,
            Self::ShortOrSemiShortCuts => !long_cuts,
        }
    }

    fn conj(&self) -> Self {
        match self {
            Self::ImXp(sign) => Self::ImXm(-sign),
            Self::ImXm(sign) => Self::ImXp(-sign),
            Self::LogBranch(b) => Self::LogBranch(*b),
            Self::EBranch(b) => Self::EBranch(*b),
            Self::UpBranch(b) => Self::UmBranch(b.clone()),
            Self::UmBranch(b) => Self::UpBranch(b.clone()),
            Self::ImXpOrUpBranch(b1, b2) => Self::ImXmOrUmBranch(-*b1, b2.clone()),
            Self::ImXmOrUmBranch(b1, b2) => Self::ImXpOrUpBranch(-*b1, b2.clone()),

            Self::LongCuts => Self::LongCuts,
            Self::SemiShortCuts => Self::SemiShortCuts,
            Self::ShortCuts => Self::ShortCuts,
            Self::ShortOrSemiShortCuts => Self::ShortOrSemiShortCuts,
        }
    }
}

#[derive(Debug, Clone)]
pub struct Cut {
    pub component: Component,
    pub paths: Vec<Vec<Complex64>>,
    pub branch_point: Option<Complex64>,
    pub typ: CutType,
    pub p_range: i32,
    pub period: Option<Complex64>,
    visibility: Vec<CutVisibilityCondition>,
}

impl Cut {
    fn new(
        component: Component,
        paths: Vec<Vec<Complex64>>,
        branch_point: Option<Complex64>,
        typ: CutType,
        p_range: i32,
        period: Option<Complex64>,
        visibility: Vec<CutVisibilityCondition>,
    ) -> Self {
        Self {
            component,
            paths,
            branch_point,
            typ,
            p_range,
            period,
            visibility,
        }
    }

    fn conj(&self) -> Self {
        let paths = self
            .paths
            .iter()
            .map(|path| path.iter().rev().map(|z| z.conj()).collect())
            .collect();
        let branch_point = self.branch_point.map(|z| z.conj());
        let visibility = self.visibility.iter().map(|v| v.conj()).collect();

        Cut {
            component: self.component.conj(),
            paths,
            branch_point,
            typ: self.typ.conj(),
            visibility,
            period: self.period,
            p_range: self.p_range,
        }
    }

    fn shift_conj(&self, dz: Complex64) -> Self {
        let paths = self
            .paths
            .iter()
            .map(|path| path.iter().map(|z| (z - dz).conj() + dz).collect())
            .collect();
        let branch_point = self.branch_point.map(|z| (z - dz).conj() + dz);
        let visibility = self.visibility.iter().map(|v| v.conj()).collect();
        Cut {
            component: self.component.conj(),
            paths,
            branch_point,
            typ: self.typ.conj(),
            visibility,
            period: self.period,
            p_range: self.p_range,
        }
    }

    fn shift(mut self, dz: Complex64) -> Self {
        for path in self.paths.iter_mut() {
            for z in path.iter_mut() {
                *z += dz;
            }
        }

        if let Some(ref mut z) = self.branch_point {
            *z += dz;
        }
        self
    }

    pub fn intersection(&self, p1: Complex64, p2: Complex64) -> Option<(usize, usize, Complex64)> {
        if let Some(period) = self.period {
            (-5..=5).find_map(|n| {
                let shift = n as f64 * period;
                self.find_intersection(p1 + shift, p2 + shift)
            })
        } else {
            self.find_intersection(p1, p2)
        }
    }

    fn find_intersection(&self, p1: Complex64, p2: Complex64) -> Option<(usize, usize, Complex64)> {
        fn cross(v: Complex64, w: Complex64) -> f64 {
            v.re * w.im - v.im * w.re
        }

        let p = p1;
        let r = p2 - p1;

        for (i, path) in self.paths.iter().enumerate() {
            for (j, (q1, q2)) in path.iter().tuple_windows::<(_, _)>().enumerate() {
                let q = q1;
                let s = q2 - q1;

                if cross(r, s) != 0.0 {
                    let t = cross(q - p, s) / cross(r, s);
                    let u = cross(q - p, r) / cross(r, s);

                    if (0.0..=1.0).contains(&t) && (0.0..=1.0).contains(&u) {
                        return Some((i, j, p + t * r));
                    }
                }
            }
        }
        None
    }

    pub fn is_visible(&self, pt: &Point, u_cut_type: UCutType) -> bool {
        self.visibility
            .iter()
            .all(|cond| cond.check(pt, u_cut_type))
    }
}

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
        contour_generator: &ContourGenerator,
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

#[derive(Debug, Clone)]
pub struct Point {
    pub p: Complex64,
    pub xp: Complex64,
    pub xm: Complex64,
    pub u: Complex64,
    pub consts: CouplingConstants,
    pub sheet_data: SheetData,
}

impl Point {
    pub fn new(p: impl Into<Complex64>, consts: CouplingConstants) -> Self {
        let p: Complex64 = p.into();
        let log_branch_p = 0;
        let log_branch_m = p.re.floor() as i32;
        let u_branch = if log_branch_m >= 0 {
            (UBranch::Outside, UBranch::Outside)
        } else if log_branch_m == -1 {
            (UBranch::Between, UBranch::Between)
        } else {
            (UBranch::Inside, UBranch::Inside)
        };

        let sheet_data = SheetData {
            log_branch_p,
            log_branch_m,
            e_branch: 1,
            u_branch,
        };

        let xp = xp(p, 1.0, consts);
        let xm = xm(p, 1.0, consts);
        let u = u(p, consts, &sheet_data);
        Self {
            p,
            xp,
            xm,
            u,
            consts,
            sheet_data,
        }
    }

    pub fn set_coupling_constants(&mut self, consts: CouplingConstants) {
        self.consts = consts;
        self.set(self.p);
    }

    fn set(&mut self, p: Complex64) {
        self.p = p;
        if self.sheet_data.e_branch > 0 {
            self.xp = xp(p, 1.0, self.consts);
            self.xm = xm(p, 1.0, self.consts);

            self.u = u(p, self.consts, &self.sheet_data);
        } else {
            self.xp = xp_crossed(p, 1.0, self.consts);
            self.xm = xm_crossed(p, 1.0, self.consts);
            self.u = u_crossed(p, self.consts, &self.sheet_data);
        };
    }

    fn try_set(&mut self, p: Option<Complex64>, sheet_data: &SheetData) -> bool {
        let Some(p) = p else {return false};
        let new_xp: Complex64;
        let new_xm: Complex64;
        let new_u: Complex64;

        if sheet_data.e_branch > 0 {
            new_xp = xp(p, 1.0, self.consts);
            new_xm = xm(p, 1.0, self.consts);
            new_u = u(p, self.consts, sheet_data);
        } else {
            new_xp = xp_crossed(p, 1.0, self.consts);
            new_xm = xm_crossed(p, 1.0, self.consts);
            new_u = u_crossed(p, self.consts, sheet_data);
        }

        if (self.p - p).re.abs() > 0.125 || (self.p - p).im.abs() > 0.25 {
            log::debug!(
                "p jump too large {} {}",
                (self.p - p).norm_sqr(),
                (self.p - p).re.abs()
            );
            return false;
        }

        if (self.xp - new_xp).norm_sqr() > 16.0 / (self.consts.h * self.consts.h) {
            log::debug!(
                "xp jump too large: {} ({}) {} ({})",
                (self.xp - new_xp).norm_sqr(),
                (self.xp - new_xp).norm_sqr() * (self.consts.h * self.consts.h),
                self.xp.norm_sqr(),
                self.xp.norm_sqr() * (self.consts.h * self.consts.h)
            );
            // return false;
        }

        if (self.xm - new_xm).norm_sqr() > 16.0 / (self.consts.h * self.consts.h) {
            log::debug!(
                "xm jump too large: {} ({}) {} ({})",
                (self.xm - new_xm).norm_sqr(),
                (self.xm - new_xm).norm_sqr() * (self.consts.h * self.consts.h),
                self.xm.norm_sqr(),
                self.xm.norm_sqr() * (self.consts.h * self.consts.h)
            );

            // return false;
        }

        if (self.u - new_u).norm_sqr() > 16.0 / (self.consts.h * self.consts.h) {
            log::debug!("u jump too large");
            // return false;
        }

        self.sheet_data = sheet_data.clone();
        self.p = p;
        self.xp = new_xp;
        self.xm = new_xm;
        self.u = new_u;

        true
    }

    fn shift_xp(
        &self,
        new_xp: Complex64,
        sheet_data: &SheetData,
        guess: Complex64,
    ) -> Option<Complex64> {
        if sheet_data.e_branch > 0 {
            nr::find_root(
                |p| xp(p, 1.0, self.consts) - new_xp,
                |p| dxp_dp(p, 1.0, self.consts),
                guess,
                1.0e-6,
                50,
            )
        } else {
            nr::find_root(
                |p| xp_crossed(p, 1.0, self.consts) - new_xp,
                |p| dxp_crossed_dp(p, 1.0, self.consts),
                guess,
                1.0e-6,
                50,
            )
        }
    }

    fn shift_xm(
        &self,
        new_xm: Complex64,
        sheet_data: &SheetData,
        guess: Complex64,
    ) -> Option<Complex64> {
        if sheet_data.e_branch > 0 {
            nr::find_root(
                |p| xm(p, 1.0, self.consts) - new_xm,
                |p| dxm_dp(p, 1.0, self.consts),
                guess,
                1.0e-6,
                50,
            )
        } else {
            nr::find_root(
                |p| xm_crossed(p, 1.0, self.consts) - new_xm,
                |p| dxm_crossed_dp(p, 1.0, self.consts),
                guess,
                1.0e-6,
                50,
            )
        }
    }

    fn shift_u(
        &self,
        new_u: Complex64,
        sheet_data: &SheetData,
        guess: Complex64,
    ) -> Option<Complex64> {
        if sheet_data.e_branch > 0 {
            nr::find_root(
                |p| u(p, self.consts, sheet_data) - new_u,
                |p| du_dp(p, self.consts, sheet_data),
                guess,
                1.0e-6,
                50,
            )
        } else {
            nr::find_root(
                |p| u_crossed(p, self.consts, sheet_data) - new_u,
                |p| du_crossed_dp(p, self.consts, sheet_data),
                guess,
                1.0e-6,
                50,
            )
        }
    }

    pub fn get(&self, component: Component) -> Complex64 {
        match component {
            Component::P => self.p,
            Component::U => self.u,
            Component::Xp => self.xp,
            Component::Xm => self.xm,
        }
    }

    pub fn update(&mut self, component: Component, new_value: Complex64, crossed_cuts: &[&Cut]) {
        let mut new_sheet_data = self.sheet_data.clone();
        for cut in crossed_cuts {
            match cut.typ {
                CutType::E => {
                    new_sheet_data.e_branch = -new_sheet_data.e_branch;
                }
                CutType::UShortScallion(Component::Xp) => {
                    new_sheet_data.u_branch = (
                        new_sheet_data.u_branch.0.cross_scallion(),
                        new_sheet_data.u_branch.1,
                    );
                }
                CutType::UShortScallion(Component::Xm) => {
                    new_sheet_data.u_branch = (
                        new_sheet_data.u_branch.0,
                        new_sheet_data.u_branch.1.cross_scallion(),
                    );
                }
                CutType::UShortKidney(Component::Xp) => {
                    new_sheet_data.u_branch = (
                        new_sheet_data.u_branch.0.cross_kidney(),
                        new_sheet_data.u_branch.1,
                    );
                }
                CutType::UShortKidney(Component::Xm) => {
                    new_sheet_data.u_branch = (
                        new_sheet_data.u_branch.0,
                        new_sheet_data.u_branch.1.cross_kidney(),
                    );
                }
                CutType::Log(Component::Xp) => {
                    if self.xp.im >= 0.0 {
                        new_sheet_data.log_branch_p += 1;
                    } else {
                        new_sheet_data.log_branch_p -= 1;
                    }
                }
                CutType::Log(Component::Xm) => {
                    if self.xm.im <= 0.0 {
                        new_sheet_data.log_branch_m += 1;
                    } else {
                        new_sheet_data.log_branch_m -= 1;
                    }
                }
                _ => {}
            }
            log::info!("Intersection with {:?}: {:?}", cut.typ, new_sheet_data);
        }

        for guess in vec![
            self.p,
            self.p - 0.01,
            self.p + 0.01,
            self.p - 0.05,
            self.p + 0.05,
        ]
        .into_iter()
        {
            let p = match component {
                Component::P => Some(new_value),
                Component::Xp => self.shift_xp(new_value, &new_sheet_data, guess),
                Component::Xm => self.shift_xm(new_value, &new_sheet_data, guess),
                Component::U => self.shift_u(new_value, &new_sheet_data, guess),
            };

            if self.try_set(p, &new_sheet_data) {
                return;
            }
        }
    }

    pub fn same_sheet(&self, other: &Point, component: Component, u_cut_type: UCutType) -> bool {
        let sd1 = &self.sheet_data;
        let sd2 = &other.sheet_data;

        match component {
            Component::P => sd1.e_branch == sd2.e_branch,
            Component::U => {
                if (sd1.log_branch_p + sd1.log_branch_m) != (sd2.log_branch_p + sd2.log_branch_m)
                    || self.consts.k() > 0
                        && (sd1.log_branch_p - sd1.log_branch_m)
                            != (sd2.log_branch_p - sd2.log_branch_m)
                {
                    false
                } else if u_cut_type == UCutType::Long {
                    self.xp.im.signum() == other.xp.im.signum()
                        && self.xm.im.signum() == other.xm.im.signum()
                } else {
                    sd1.u_branch == sd2.u_branch
                }
            }
            Component::Xp => {
                if (sd1.log_branch_p + sd1.log_branch_m) != (sd2.log_branch_p + sd2.log_branch_m) {
                    false
                } else if u_cut_type == UCutType::Long {
                    self.xm.im.signum() == other.xm.im.signum()
                } else {
                    sd1.u_branch.1 == sd2.u_branch.1
                }
            }
            Component::Xm => {
                if (sd1.log_branch_p + sd1.log_branch_m) != (sd2.log_branch_p + sd2.log_branch_m) {
                    false
                } else if u_cut_type == UCutType::Long {
                    self.xp.im.signum() == other.xp.im.signum()
                } else {
                    sd1.u_branch.0 == sd2.u_branch.0
                }
            }
        }
    }
}
