use std::collections::VecDeque;

use crate::kinematics::{
    den2_dp, du_crossed_dp, du_dp, dxm_crossed_dp, dxm_dp, dxp_crossed_dp, dxp_dp, en2, u,
    u_crossed, xm, xm_crossed, xp, xp_crossed, CouplingConstants, SheetData,
};
use crate::nr::{self};
use crate::pxu2::{PInterpolatorMut, XInterpolator};
use itertools::Itertools;
use num::complex::Complex;

type C = Complex<f64>;

const P_RANGE_MIN: i32 = -3;
const P_RANGE_MAX: i32 = 3;

const INFINITY: f64 = 1000.0;

#[derive(Debug, Clone, Copy, PartialEq)]
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

#[derive(Debug, Clone, Copy)]
enum BranchPoint {
    XpPositiveAxisImXmNegative,
    XpPositiveAxisImXmPositive,
    XpNegativeAxisFromAboveWithImXmNegative,
    XpNegativeAxisFromBelowWithImXmNegative,
    XpNegativeAxisFromAboveWithImXmPositive,
    XpNegativeAxisFromBelowWithImXmPositive,
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

#[derive(Debug)]
enum GeneratorCommands {
    AddGridLineU(f64),
    AddGridLineXReal(f64),
    AddGridLineX(f64),
    AddGridLineP,

    ComputeBranchPoint(i32, BranchPoint),
    ClearCut,
    ComputeCutX(i32, CutDirection),
    ComputeCutXFull(XCut),
    ComputeCutP(bool),
    ComputeCutEP(i32),
    SetCutPath(Vec<C>, Option<C>),
    PushCut(i32, Component, CutType, Vec<CutVisibilityCondition>),
    PushCutFromP(i32, Component, CutType, Vec<CutVisibilityCondition>),

    PStartXp(f64),
    PGotoXp(f64, f64),
    PGotoXm(f64, f64),
    PGotoP(f64),
    PGotoM(f64),
    PGotoIm(f64),
    PGotoRe(f64),
}

struct RuntimeCutData {
    branch_point: Option<C>,
    path: Option<Vec<C>>,
}

struct BuildTimeCutData {
    component: Option<Component>,
    cut_type: Option<CutType>,
    visibility: Vec<CutVisibilityCondition>,
}

struct ContourGeneratorRuntimeContext {
    p_int: Option<PInterpolatorMut>,
    branch_point_data: Option<(f64, f64, BranchPoint)>, // (p, m, branch_point_type)
    cut_data: RuntimeCutData,
}

impl ContourGeneratorRuntimeContext {
    fn new() -> Self {
        Self {
            p_int: None,
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
                visibility: vec![],
            },
        }
    }

    fn clear(&mut self) {
        self.cut_data.component = None;
        self.cut_data.cut_type = None;
        self.cut_data.visibility.clear();
    }
}

pub struct ContourGenerator {
    cuts: Vec<Cut>,
    commands: VecDeque<GeneratorCommands>,
    consts: Option<CouplingConstants>,

    grid_p: Vec<Vec<C>>,
    grid_x: Vec<Vec<C>>,
    grid_u: Vec<Vec<C>>,

    rctx: ContourGeneratorRuntimeContext,
    bctx: ContourGeneratorBuildTimeContext,

    num_commands: usize,
}

impl ContourGenerator {
    pub fn new() -> Self {
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

    pub fn update(&mut self, pt: &PxuPoint) -> bool {
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
            log::info!("Generated {} commands", self.num_commands,)
        }

        let start = chrono::Utc::now();
        while !self.commands.is_empty()
            && (chrono::Utc::now() - start).num_milliseconds() < (1000.0 / 20.0f64).floor() as i64
        {
            if let Some(command) = self.commands.pop_front() {
                self.execute(command);
            }
        }

        self.commands.is_empty()
    }

    fn clear(&mut self) {
        log::info!("Clearing grid and cuts");
        self.commands.clear();
        self.grid_x.clear();
        self.grid_u.clear();
        self.cuts.clear();

        self.grid_p = vec![vec![
            C::from(P_RANGE_MIN as f64),
            C::from(P_RANGE_MAX as f64),
        ]];
    }

    pub fn progress(&self) -> f32 {
        if self.num_commands > 0 {
            1.0 - self.commands.len() as f32 / self.num_commands as f32
        } else {
            0.0
        }
    }

    pub fn get_grid(&self, component: Component) -> &Vec<Vec<C>> {
        match component {
            Component::P => &self.grid_p,
            Component::Xp | Component::Xm => &self.grid_x,
            Component::U => &self.grid_u,
        }
    }
    pub fn get_visible_cuts(
        &mut self,
        pt: &PxuPoint,
        component: Component,
    ) -> impl Iterator<Item = &Cut> {
        let mut pt = pt.clone();
        pt.u += 2.0 * (pt.sheet_data.log_branch_p * pt.consts.k()) as f64 * C::i() / pt.consts.h;

        self.cuts
            .iter()
            .filter(move |c| c.component == component && c.is_visible(&pt))
    }

    pub fn get_crossed_cuts(
        &mut self,
        pt: &PxuPoint,
        component: Component,
        new_value: C,
    ) -> impl Iterator<Item = &Cut> {
        let mut pt = pt.clone();
        pt.u += 2.0 * (pt.sheet_data.log_branch_p * pt.consts.k()) as f64 * C::i() / pt.consts.h;

        let new_value = if component == Component::U {
            new_value + (pt.sheet_data.log_branch_p * pt.consts.k()) as f64 * C::i() / pt.consts.h
        } else {
            new_value
        };

        self.cuts.iter().filter(move |c| {
            c.component == component
                && c.is_visible(&pt)
                && c.intersection(pt.get(component), new_value).is_some()
        })
    }

    fn execute(&mut self, command: GeneratorCommands) {
        use GeneratorCommands::*;

        let Some(consts) = self.consts else {
            log::warn!("Executing commands when consts is not set!");
            return;
        };

        fn branch_point_mass(p_start: f64, k: f64, branch_point_type: BranchPoint) -> f64 {
            match branch_point_type {
                BranchPoint::XpPositiveAxisImXmNegative => 2.0 * p_start * k + 2.0,
                BranchPoint::XpPositiveAxisImXmPositive => -(2.0 * p_start * k + 2.0),

                BranchPoint::XpNegativeAxisFromAboveWithImXmNegative => {
                    (2.0 * p_start + 1.0) * k + 2.0
                }
                BranchPoint::XpNegativeAxisFromBelowWithImXmNegative => {
                    (2.0 * p_start - 1.0) * k + 2.0
                }

                BranchPoint::XpNegativeAxisFromAboveWithImXmPositive => {
                    -((2.0 * p_start + 1.0) * k + 2.0)
                }

                BranchPoint::XpNegativeAxisFromBelowWithImXmPositive => {
                    -((2.0 * p_start - 1.0) * k + 2.0)
                }
            }
        }

        match command {
            AddGridLineU(y) => {
                self.grid_u
                    .push(vec![C::new(-INFINITY, y), C::new(INFINITY, y)]);
            }

            AddGridLineX(m) => {
                let path = XInterpolator::generate_xp_full(0, m, consts);
                self.grid_x.push(path.iter().map(|x| x.conj()).collect());
                self.grid_x.push(path);
            }

            AddGridLineXReal(x) => {
                if x > 0.0 {
                    self.grid_x.push(vec![C::from(x), C::from(INFINITY)]);
                } else {
                    self.grid_x.push(vec![C::from(x), C::from(-INFINITY)]);
                }
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
                self.grid_p.push(path.iter().map(|p| p.conj()).collect());
                self.grid_p.push(path);
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
                let p_start = p_range as f64;
                let k = consts.k() as f64;
                let s = consts.s();
                let u_of_x = |x: C| -> C { x + 1.0 / x - (s - 1.0 / s) * x.ln() };
                let du_dx = |x: C| -> C { (x - s) * (x + 1.0 / s) / (x * x) };

                let u_of_s = u_of_x(C::from(s))
                    * match branch_point_type {
                        BranchPoint::XpPositiveAxisImXmNegative
                        | BranchPoint::XpPositiveAxisImXmPositive => 1.0,
                        BranchPoint::XpNegativeAxisFromAboveWithImXmNegative
                        | BranchPoint::XpNegativeAxisFromAboveWithImXmPositive
                        | BranchPoint::XpNegativeAxisFromBelowWithImXmNegative
                        | BranchPoint::XpNegativeAxisFromBelowWithImXmPositive => -1.0,
                    };

                let m = branch_point_mass(p_start, k, branch_point_type);
                let guess = xp(0.5, m, consts);

                let x_branch_point = nr::find_root(
                    |x| u_of_x(x) - u_of_s - m * C::i() / consts.h,
                    du_dx,
                    guess,
                    1.0e-3,
                    10,
                );

                if let Some(x_branch_point) = x_branch_point {
                    let p = x_branch_point.arg().abs() / std::f64::consts::PI;
                    self.rctx.branch_point_data = Some((p, m, branch_point_type));
                } else {
                    log::info!("Could not find branch point");
                    self.rctx.branch_point_data = None;
                };
            }

            ComputeCutX(p_range, cut_direction) => {
                self.rctx.cut_data.path = None;
                self.rctx.cut_data.branch_point = None;

                let Some((p_branch_point, m, branch_point_type)) = self.rctx.branch_point_data else {
                    log::info!("No branch point set");
                    return;
                };

                let (p_start, p_end) = match cut_direction {
                    CutDirection::Positive => (0.0, p_branch_point),
                    CutDirection::Negative => (p_branch_point, 1.0),
                };

                let path = match branch_point_type {
                    BranchPoint::XpPositiveAxisImXmNegative
                    | BranchPoint::XpNegativeAxisFromAboveWithImXmNegative
                    | BranchPoint::XpNegativeAxisFromBelowWithImXmNegative => {
                        XInterpolator::generate_xm(p_start, p_end, m, consts)
                    }

                    BranchPoint::XpPositiveAxisImXmPositive
                    | BranchPoint::XpNegativeAxisFromAboveWithImXmPositive
                    | BranchPoint::XpNegativeAxisFromBelowWithImXmPositive => {
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
                    XCut::Scallion => C::from(consts.s()),
                    XCut::Kidney => C::from(-1.0 / consts.s()),
                });
            }

            ComputeCutEP(p_range) => {
                self.rctx.cut_data.path = None;
                self.rctx.cut_data.branch_point = None;

                let p_start = p_range as f64;

                let p0 = nr::find_root(
                    |p| en2(p, 1.0, consts),
                    |p| den2_dp(p, 1.0, consts),
                    C::new(p_start, 2.5),
                    1.0e-3,
                    50,
                );

                self.rctx.cut_data.branch_point = p0;
                let Some(p0) = p0 else { return };

                let mut path = vec![];

                path.push((0.0, p0));
                let mut p_prev = p0;

                const STEP: i32 = 16;

                for i in 1.. {
                    let im = i as f64 * i as f64 / (STEP as f64);

                    let p = nr::find_root(
                        |p| en2(p, 1.0, consts) - C::new(-im, 0.0),
                        |p| den2_dp(p, 1.0, consts),
                        p_prev,
                        1.0e-3,
                        50,
                    );

                    let Some(p) = p else {break;};

                    path.push((im, p));
                    p_prev = p;

                    if p.im.abs() > 2.0 {
                        break;
                    }
                }

                self.rctx.cut_data.path = Some(path.into_iter().map(|(_, p)| p).collect());
            }

            SetCutPath(path, branchpoint) => {
                self.rctx.cut_data.path = Some(path);
                self.rctx.cut_data.branch_point = branchpoint;
            }

            PushCut(p_range, component, cut_type, visibility) => {
                let Some(ref path) = self.rctx.cut_data.path else {
                    log::info!("No path for cut");
                    return;
                };

                if self.rctx.cut_data.path.is_none() {
                    log::info!("No path for cut");
                    return;
                };

                let shift = match component {
                    Component::U => C::new(0.0, (p_range * consts.k()) as f64 / consts.h),
                    _ => C::from(0.0),
                };

                let mut cut = Cut::new(
                    component,
                    vec![path.clone()],
                    self.rctx.cut_data.branch_point.clone(),
                    cut_type,
                );
                cut.visibility = visibility;

                self.cuts.push(cut.conj().shift(shift));
                self.cuts.push(cut.shift(shift));
            }

            PushCutFromP(p_range, component, cut_type, visibility) => {
                let Some(ref p_path) = self.rctx.cut_data.path else {
                    log::info!("No path for cut");
                    return;
                };

                if self.rctx.cut_data.path.is_none() {
                    log::info!("No path for cut");
                    return;
                };

                let sheet_data = SheetData {
                    log_branch_p: 0,
                    log_branch_m: p_range,
                    e_branch: 1,
                    u_branch: (1, 1),
                };

                let mut path = p_path
                    .iter()
                    .rev()
                    .map(|&p| match component {
                        Component::P => unimplemented!(),
                        Component::Xp => xp(p + 0.00001, 1.0, consts),
                        Component::Xm => xm(p + 0.00001, 1.0, consts),
                        Component::U => u(p + 0.00001, consts, &sheet_data),
                    })
                    .collect::<Vec<_>>();

                path.extend(p_path.iter().map(|&p| match component {
                    Component::P => unimplemented!(),
                    Component::Xp => xp(p - 0.00001, 1.0, consts),
                    Component::Xm => xm(p - 0.00001, 1.0, consts),
                    Component::U => u(p - 0.00001, consts, &sheet_data),
                }));

                let branch_point = if let Some(p) = self.rctx.cut_data.branch_point {
                    Some(match component {
                        Component::P => unimplemented!(),
                        Component::Xp => xp(p, 1.0, consts),
                        Component::Xm => xm(p, 1.0, consts),
                        Component::U => u(p, consts, &sheet_data),
                    })
                } else {
                    None
                };

                let mut cut = Cut::new(component, vec![path], branch_point, cut_type);
                cut.visibility = visibility;

                self.cuts.push(cut.conj());
                self.cuts.push(cut);
            }

            _ => {
                log::info!("Command {:?} not implemented", command);
            }
        }
    }

    fn add(&mut self, command: GeneratorCommands) -> &mut Self {
        self.commands.push_back(command);
        self
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

    fn generate_commands(&mut self, pt: &PxuPoint) {
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

        for y in 1..=(32 * consts.k()) {
            self.add(GeneratorCommands::AddGridLineU(y as f64 / consts.h));
            self.add(GeneratorCommands::AddGridLineU(-y as f64 / consts.h));
        }
    }

    fn generate_x_grid(&mut self, p_range: i32, consts: CouplingConstants) {
        for m in (p_range * consts.k())..=((p_range + 1) * consts.k()) {
            self.add(GeneratorCommands::AddGridLineX(m as f64));
        }

        if p_range == 0 {
            self.add(GeneratorCommands::AddGridLineXReal(consts.s()));
        }

        if p_range == -1 {
            self.add(GeneratorCommands::AddGridLineXReal(-1.0 / consts.s()));
        }
    }

    fn generate_p_grid(&mut self, p_range: i32, consts: CouplingConstants) {
        let p_start = p_range as f64;
        let k = consts.k() as f64;
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

            for m in 3..=(2 * consts.k()) {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            self.p_start_xp(p2).goto_xp(p2, 3.0);

            for m in 3..=(consts.k() + 1) {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            for m in (consts.k() + 3)..=(6 * consts.k()) {
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

            self.p_start_xp(p2).goto_xp(p2, p_start * k + 3.0);

            for m in (p_range * consts.k() + 3)..=((p_range + 1) * consts.k() + 1) {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            for m in ((p_range + 1) * consts.k() + 3)..=(6 * consts.k()) {
                self.goto_xp(p0, m as f64).p_grid_line();
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

            for m in (consts.k() + 1)..=(6 * consts.k()) {
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

            for m in (-p_range * consts.k() + 1)..=(6 * consts.k()) {
                self.goto_xp(p0, m as f64).p_grid_line();
            }

            self.p_start_xp(p0);

            for m in 1..=(-(p_range + 1) * consts.k() - 1) {
                self.goto_xm(p0, m as f64).p_grid_line();
            }

            for m in (-(p_range + 1) * consts.k() + 1)..=(-p_range * consts.k() - 1) {
                self.goto_xm(p0, m as f64).p_grid_line();
            }

            for m in (-p_range * consts.k() + 1)..=(6 * consts.k()) {
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

    fn clear_cut(&mut self) -> &mut Self {
        self.add(GeneratorCommands::ClearCut)
    }

    fn compute_branch_point(&mut self, p_range: i32, branch_point_type: BranchPoint) -> &mut Self {
        self.add(GeneratorCommands::ComputeBranchPoint(
            p_range,
            branch_point_type,
        ))
    }

    fn compute_cut_path_x(&mut self, p_range: i32, direction: CutDirection) -> &mut Self {
        self.add(GeneratorCommands::ComputeCutX(p_range, direction))
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

    fn set_cut_path(&mut self, path: Vec<C>, branchpoint: Option<C>) -> &mut Self {
        self.add(GeneratorCommands::SetCutPath(path, branchpoint))
    }

    fn push_cut(&mut self, p_range: i32) -> &mut Self {
        let Some(component) = std::mem::replace(&mut self.bctx.cut_data.component, None) else {
            log::info!("Can't push cut without component");
            self.bctx.clear();
            return self;
        };
        let Some(cut_type) = std::mem::replace(&mut self.bctx.cut_data.cut_type, None) else {
            log::info!("Can't push cut without type");
            self.bctx.clear();
            return self;
        };
        self.add(GeneratorCommands::PushCut(
            p_range,
            component,
            cut_type,
            self.bctx.cut_data.visibility.clone(),
        ));
        self.bctx.clear();
        self
    }

    fn create_cut(&mut self, component: Component, cut_type: CutType) -> &mut Self {
        if self.bctx.cut_data.component.is_some() || self.bctx.cut_data.cut_type.is_some() {
            log::info!("New cut created before previous cut was pushed");
        }
        self.bctx.cut_data.component = Some(component);
        self.bctx.cut_data.cut_type = Some(cut_type);
        self
    }

    fn log_branch(&mut self, p_range: i32) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::LogBranch(p_range));
        self
    }

    fn im_xm_positive(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::ImXm(1));
        self
    }

    fn im_xm_negative(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::ImXm(-1));
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
            .push(CutVisibilityCondition::UpBranch(1));
        self
    }

    fn xp_inside(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::UpBranch(-1));
        self
    }

    fn xm_outside(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::UmBranch(1));
        self
    }

    fn xm_inside(&mut self) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::UmBranch(-1));
        self
    }

    fn e_branch(&mut self, branch: i32) -> &mut Self {
        self.bctx
            .cut_data
            .visibility
            .push(CutVisibilityCondition::EBranch(branch));
        self
    }

    fn generate_cuts(&mut self, p_range: i32, consts: CouplingConstants) {
        log::info!("{p_range}");

        let p_start = p_range as f64;
        let k = consts.k() as f64;
        let s = consts.s();

        let us = s + 1.0 / s - (s - 1.0 / s) * s.ln();

        {
            // Log
            self.clear_cut()
                .set_cut_path(vec![C::from(-INFINITY), C::from(0.0)], Some(C::from(0.0)))
                .create_cut(Component::Xp, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPoint::XpNegativeAxisFromAboveWithImXmNegative,
                )
                .compute_cut_path_x(p_range, CutDirection::Negative)
                .create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPoint::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(p_range, CutDirection::Negative)
                .create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPoint::XpNegativeAxisFromBelowWithImXmNegative,
                )
                .compute_cut_path_x(p_range, CutDirection::Negative)
                .create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPoint::XpNegativeAxisFromBelowWithImXmPositive,
                )
                .compute_cut_path_x(p_range, CutDirection::Negative)
                .create_cut(Component::Xm, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        C::new(-INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                        C::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                    ],
                    Some(C::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h)),
                )
                .create_cut(Component::U, CutType::Log(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        C::new(-INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        C::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                    ],
                    Some(C::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h)),
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
                .set_cut_path(vec![C::from(0.0), C::from(INFINITY)], Some(C::from(s)))
                .create_cut(Component::Xp, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(p_range, BranchPoint::XpPositiveAxisImXmNegative)
                .compute_cut_path_x(p_range, CutDirection::Positive)
                .create_cut(Component::Xm, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(p_range, BranchPoint::XpPositiveAxisImXmPositive)
                .compute_cut_path_x(p_range, CutDirection::Positive)
                .create_cut(Component::Xm, CutType::ULongPositive(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        C::new(INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                        C::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                    ],
                    Some(C::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
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
                    vec![C::from(-INFINITY), C::from(0.0)],
                    Some(-1.0 / C::from(s)),
                )
                .create_cut(Component::Xp, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPoint::XpNegativeAxisFromAboveWithImXmNegative,
                )
                .compute_cut_path_x(p_range, CutDirection::Negative)
                .create_cut(Component::Xm, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPoint::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(p_range, CutDirection::Negative)
                .create_cut(Component::Xm, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPoint::XpNegativeAxisFromBelowWithImXmNegative,
                )
                .compute_cut_path_x(p_range, CutDirection::Negative)
                .create_cut(Component::Xm, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPoint::XpNegativeAxisFromBelowWithImXmPositive,
                )
                .compute_cut_path_x(p_range, CutDirection::Negative)
                .create_cut(Component::Xm, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_negative()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        C::new(-INFINITY, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                        C::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h),
                    ],
                    Some(C::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h)),
                )
                .create_cut(Component::U, CutType::ULongNegative(Component::Xp))
                .log_branch(p_range)
                .im_xp_positive()
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        C::new(-INFINITY, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                        C::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h),
                    ],
                    Some(C::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h)),
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

            self.clear_cut()
                .compute_branch_point(p_range, BranchPoint::XpPositiveAxisImXmPositive)
                .compute_cut_path_x(p_range, CutDirection::Negative)
                .create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(p_range, BranchPoint::XpPositiveAxisImXmNegative)
                .compute_cut_path_x(p_range, CutDirection::Negative)
                .create_cut(Component::Xm, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .set_cut_path(
                    vec![
                        C::new(-INFINITY, -(1.0 + p_range as f64 * k) / consts.h),
                        C::new(us, -(1.0 + p_range as f64 * k) / consts.h),
                    ],
                    Some(C::new(us, -(1.0 + p_range as f64 * k) / consts.h)),
                )
                .create_cut(Component::U, CutType::UShortScallion(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

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
                    BranchPoint::XpNegativeAxisFromAboveWithImXmNegative,
                )
                .compute_cut_path_x(p_range, CutDirection::Positive)
                .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPoint::XpNegativeAxisFromAboveWithImXmPositive,
                )
                .compute_cut_path_x(p_range, CutDirection::Positive)
                .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPoint::XpNegativeAxisFromBelowWithImXmNegative,
                )
                .compute_cut_path_x(p_range, CutDirection::Positive)
                .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
                .push_cut(p_range);

            self.clear_cut()
                .compute_branch_point(
                    p_range,
                    BranchPoint::XpNegativeAxisFromBelowWithImXmPositive,
                )
                .compute_cut_path_x(p_range, CutDirection::Positive)
                .create_cut(Component::Xm, CutType::UShortKidney(Component::Xp))
                .log_branch(p_range)
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

        self.add(GeneratorCommands::ComputeCutEP(p_range));
        self.create_cut(Component::P, CutType::E).push_cut(p_range);

        self.create_cut(Component::Xp, CutType::E)
            .add(GeneratorCommands::PushCutFromP(
                p_range,
                Component::Xp,
                CutType::E,
                vec![
                    CutVisibilityCondition::LogBranch(p_range),
                    CutVisibilityCondition::ImXm(-1),
                ],
            ));

        self.create_cut(Component::Xm, CutType::E)
            .add(GeneratorCommands::PushCutFromP(
                p_range,
                Component::Xm,
                CutType::E,
                vec![
                    CutVisibilityCondition::LogBranch(p_range),
                    CutVisibilityCondition::ImXp(-1),
                ],
            ));

        self.create_cut(Component::U, CutType::E)
            .add(GeneratorCommands::PushCutFromP(
                p_range,
                Component::U,
                CutType::E,
                vec![
                    CutVisibilityCondition::LogBranch(p_range),
                    CutVisibilityCondition::ImXp(-1),
                    CutVisibilityCondition::ImXm(-1),
                ],
            ));
    }
}

#[derive(Debug)]
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
    LogBranchGT(i32),
    LogBranchLE(i32),
    LogBranchSum(i32),
    EBranch(i32),
    UpBranch(i32),
    UmBranch(i32),
}

impl CutVisibilityCondition {
    fn check(&self, pt: &PxuPoint) -> bool {
        match self {
            Self::ImXp(sign) => pt.xp.im.signum() as i8 == sign.signum(),
            Self::ImXm(sign) => pt.xm.im.signum() as i8 == sign.signum(),
            Self::LogBranch(b) => *b == (pt.sheet_data.log_branch_p + pt.sheet_data.log_branch_m),
            Self::LogBranchGT(b) => (pt.sheet_data.log_branch_p + pt.sheet_data.log_branch_m) > *b,
            Self::LogBranchLE(b) => (pt.sheet_data.log_branch_p + pt.sheet_data.log_branch_m) <= *b,
            Self::LogBranchSum(b) => *b == pt.sheet_data.log_branch_p - pt.sheet_data.log_branch_m,
            Self::EBranch(b) => pt.sheet_data.e_branch == *b,
            Self::UpBranch(b) => pt.sheet_data.u_branch.0 == *b,
            Self::UmBranch(b) => pt.sheet_data.u_branch.1 == *b,
        }
    }

    fn conj(&self) -> Self {
        match self {
            Self::ImXp(sign) => Self::ImXm(-sign),
            Self::ImXm(sign) => Self::ImXp(-sign),
            Self::LogBranch(b) => Self::LogBranch(*b),
            Self::LogBranchGT(b) => Self::LogBranchGT(*b),
            Self::LogBranchLE(b) => Self::LogBranchLE(*b),
            Self::LogBranchSum(b) => Self::LogBranchSum(*b),
            Self::EBranch(b) => Self::EBranch(*b),
            Self::UpBranch(b) => Self::UmBranch(*b),
            Self::UmBranch(b) => Self::UpBranch(*b),
        }
    }
}

#[derive(Debug)]
pub struct Cut {
    pub component: Component,
    pub paths: Vec<Vec<C>>,
    pub branch_point: Option<C>,
    pub typ: CutType,
    visibility: Vec<CutVisibilityCondition>,
}

impl Cut {
    fn new(
        component: Component,
        paths: Vec<Vec<C>>,
        branch_point: Option<C>,
        typ: CutType,
    ) -> Self {
        Self {
            component,
            paths,
            branch_point,
            typ,
            visibility: vec![],
        }
    }

    fn conj(&self) -> Self {
        let paths = self
            .paths
            .iter()
            .map(|path| path.iter().map(|z| z.conj()).collect())
            .collect();
        let branch_point = self.branch_point.map(|z| z.conj());
        let visibility = self.visibility.iter().map(|v| v.conj()).collect();
        Cut {
            component: self.component.conj(),
            paths,
            branch_point,
            typ: self.typ.conj(),
            visibility,
        }
    }

    fn shift(mut self, dz: C) -> Self {
        let dz: C = dz.into();

        for path in self.paths.iter_mut() {
            for z in path.iter_mut() {
                *z += dz;
            }
        }

        for z in self.branch_point.iter_mut() {
            *z += dz;
        }
        self
    }

    pub fn intersection(&self, p1: C, p2: C) -> Option<C> {
        fn cross(v: C, w: C) -> f64 {
            v.re * w.im - v.im * w.re
        }

        let p = p1;
        let r = p2 - p1;

        for path in self.paths.iter() {
            for (q1, q2) in path.iter().tuple_windows::<(_, _)>() {
                let q = q1;
                let s = q2 - q1;

                if cross(r, s) != 0.0 {
                    let t = cross(q - p, s) / cross(r, s);
                    let u = cross(q - p, r) / cross(r, s);

                    if 0.0 <= t && t <= 1.0 && 0.0 <= u && u <= 1.0 {
                        return Some(p + t * r);
                    }
                }
            }
        }
        None
    }

    pub fn is_visible(&self, pt: &PxuPoint) -> bool {
        self.visibility.iter().all(|cond| cond.check(pt))
    }
}

#[derive(Debug, Clone)]
pub struct PxuPoint {
    pub p: C,
    pub xp: C,
    pub xm: C,
    pub u: C,
    pub consts: CouplingConstants,
    pub sheet_data: SheetData,
}

impl PxuPoint {
    pub fn new(p: impl Into<C>, consts: CouplingConstants) -> Self {
        let p: C = p.into();
        let log_branch_p = 0;
        let log_branch_m = p.re.floor() as i32;

        let sheet_data = SheetData {
            log_branch_p,
            log_branch_m,
            e_branch: 1,
            u_branch: (1, 1),
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

    fn set(&mut self, p: C) {
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

    fn try_set(&mut self, p: Option<C>, sheet_data: &SheetData) -> bool {
        let Some(p) = p else {return false};
        let new_xp: C;
        let new_xm: C;
        let new_u: C;

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
            log::info!(
                "p jump too large {} {}",
                (self.p - p).norm_sqr(),
                (self.p - p).re.abs()
            );
            return false;
        }

        if (self.xp - new_xp).norm_sqr() > 16.0 / (self.consts.h * self.consts.h) {
            log::info!(
                "xp jump too large: {} ({}) {} ({})",
                (self.xp - new_xp).norm_sqr(),
                (self.xp - new_xp).norm_sqr() * (self.consts.h * self.consts.h),
                self.xp.norm_sqr(),
                self.xp.norm_sqr() * (self.consts.h * self.consts.h)
            );
            // return false;
        }

        if (self.xm - new_xm).norm_sqr() > 16.0 / (self.consts.h * self.consts.h) {
            log::info!(
                "xm jump too large: {} ({}) {} ({})",
                (self.xm - new_xm).norm_sqr(),
                (self.xm - new_xm).norm_sqr() * (self.consts.h * self.consts.h),
                self.xm.norm_sqr(),
                self.xm.norm_sqr() * (self.consts.h * self.consts.h)
            );

            // return false;
        }

        if (self.u - new_u).norm_sqr() > 16.0 / (self.consts.h * self.consts.h) {
            log::info!("u jump too large");
            // return false;
        }

        self.sheet_data = sheet_data.clone();
        self.p = p;
        self.xp = new_xp;
        self.xm = new_xm;
        self.u = new_u;

        true
    }

    fn shift_xp(&self, new_xp: C, sheet_data: &SheetData, guess: C) -> Option<C> {
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

    fn shift_xm(&self, new_xm: C, sheet_data: &SheetData, guess: C) -> Option<C> {
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

    fn shift_u(&self, new_u: C, sheet_data: &SheetData, guess: C) -> Option<C> {
        if sheet_data.e_branch > 0 {
            nr::find_root(
                |p| u(p, self.consts, &sheet_data) - new_u,
                |p| du_dp(p, self.consts, &sheet_data),
                guess,
                1.0e-6,
                50,
            )
        } else {
            nr::find_root(
                |p| u_crossed(p, self.consts, &sheet_data) - new_u,
                |p| du_crossed_dp(p, self.consts, &sheet_data),
                guess,
                1.0e-6,
                50,
            )
        }
    }

    pub fn get(&self, component: Component) -> C {
        match component {
            Component::P => self.p,
            Component::U => self.u,
            Component::Xp => self.xp,
            Component::Xm => self.xm,
        }
    }

    pub fn update(&mut self, component: Component, new_value: C, crossed_cuts: &[&Cut]) {
        let mut new_sheet_data = self.sheet_data.clone();
        for cut in crossed_cuts {
            match cut.typ {
                CutType::E => {
                    new_sheet_data.e_branch = -new_sheet_data.e_branch;
                }
                CutType::UShortScallion(Component::Xp) => {
                    new_sheet_data.u_branch =
                        (-self.sheet_data.u_branch.0, self.sheet_data.u_branch.1);
                }
                CutType::UShortScallion(Component::Xm) => {
                    new_sheet_data.u_branch =
                        (self.sheet_data.u_branch.0, -self.sheet_data.u_branch.1);
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
                break;
            }
        }
    }
}
