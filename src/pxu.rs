use std::collections::VecDeque;

use crate::kinematics::{
    den2_dp, du_crossed_dp, du_dp, dxm_crossed_dp, dxm_dp, dxp_crossed_dp, dxp_dp, en2, u,
    u_crossed, xm, xm_crossed, xp, xp_crossed, CouplingConstants, SheetData,
};
use crate::nr::{self};
use crate::pxu2::{PInterpolator, PInterpolatorMut, XInterpolator};
use itertools::Itertools;
use num::complex::Complex;
use num::Zero;

type C = Complex<f64>;

const P_RANGE_MIN: i32 = -3;
const P_RANGE_MAX: i32 = 3;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Component {
    P,
    Xp,
    Xm,
    U,
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
enum GeneratorCommands {
    AddGridLineU(f64),
    AddGridLineXReal(f64),
    AddGridLineX(f64),
    AddGridLineP,

    ComputeBranchPoint(i32, BranchPoint),
    AddCutX(i32, CutDirection),
    AddCutPFull,
    AddLogCutXReal,

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
    old_cuts: Cuts,
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
            old_cuts: Cuts::new(),
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

    pub fn get_visible_old_cuts(
        &mut self,
        pt: &PxuPoint,
        component: Component,
    ) -> impl Iterator<Item = &Cut> {
        self.old_cuts.visible(pt, component, true)
    }

    pub fn get_visible_cuts(
        &mut self,
        pt: &PxuPoint,
        component: Component,
    ) -> impl Iterator<Item = &Cut> {
        let mut pt = pt.clone();
        pt.u += (pt.sheet_data.log_branch + pt.sheet_data.log_branch_sum) as f64
            * pt.consts.k() as f64
            * C::i()
            / pt.consts.h;

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
        self.old_cuts.crossed(pt, component, new_value, true)
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
                    .push(vec![C::new(-1000.0, y), C::new(1000.0, y)]);
            }

            AddGridLineX(m) => {
                let path = XInterpolator::generate_xp_full(0, m, consts);
                self.grid_x.push(path.iter().map(|x| x.conj()).collect());
                self.grid_x.push(path);
            }

            AddGridLineXReal(x) => {
                if x > 0.0 {
                    self.grid_x.push(vec![C::from(x), C::from(1000.0)]);
                } else {
                    self.grid_x.push(vec![C::from(x), C::from(-1000.0)]);
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

            AddCutPFull(reverse) => {
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

            AddLogCutXReal => {
                self.cuts.push(
                    Cut::new(
                        Component::Xp,
                        vec![vec![C::from(-1000.0), C::from(0.0)]],
                        vec![C::from(0.0)],
                        CutType::LogX(Component::Xp, 1),
                    )
                    .im_xp_positive(),
                );

                self.cuts.push(
                    Cut::new(
                        Component::Xp,
                        vec![vec![C::from(-1000.0), C::from(0.0)]],
                        vec![C::from(0.0)],
                        CutType::LogX(Component::Xp, -1),
                    )
                    .im_xp_negative(),
                );

                self.cuts.push(
                    Cut::new(
                        Component::Xm,
                        vec![vec![C::from(-1000.0), C::from(0.0)]],
                        vec![C::from(0.0)],
                        CutType::LogX(Component::Xm, 1),
                    )
                    .im_xm_negative(),
                );

                self.cuts.push(
                    Cut::new(
                        Component::Xm,
                        vec![vec![C::from(-1000.0), C::from(0.0)]],
                        vec![C::from(0.0)],
                        CutType::LogX(Component::Xm, -1),
                    )
                    .im_xm_positive(),
                );
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
                    self.ctx.branch_point = Some((p, m, branch_point_type));
                } else {
                    log::info!("Could not find branch point");
                    self.ctx.branch_point = None;
                };
            }

            AddCutX(p_range, cut_direction) => {
                let Some((p_branch_point, m, branch_point_type)) = self.ctx.branch_point else {
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

                    _ => {
                        unreachable!();
                    }
                };

                let branch_points = vec![*match cut_direction {
                    CutDirection::Positive => path.last().unwrap(),
                    CutDirection::Negative => path.first().unwrap(),
                }];

                self.cuts.push(
                    Cut::new(
                        Component::Xm,
                        vec![path],
                        branch_points,
                        CutType::LogX(Component::Xp, 1),
                    )
                    .log_branch(p_range),
                );
            }

            _ => {
                log::info!("Command {:?} not implemented", command);
            }
        }
    }

    fn add(&mut self, command: GeneratorCommands) {
        self.commands.push_back(command);
    }

    fn p_start_xp(&mut self, p: f64) -> &mut Self {
        self.add(GeneratorCommands::PStartXp(p));
        self
    }

    fn goto_xp(&mut self, p: f64, m: f64) -> &mut Self {
        self.add(GeneratorCommands::PGotoXp(p, m));
        self
    }

    fn goto_xm(&mut self, p: f64, m: f64) -> &mut Self {
        self.add(GeneratorCommands::PGotoXm(p, m));
        self
    }

    fn goto_re(&mut self, re: f64) -> &mut Self {
        self.add(GeneratorCommands::PGotoRe(re));
        self
    }

    fn goto_im(&mut self, im: f64) -> &mut Self {
        self.add(GeneratorCommands::PGotoIm(im));
        self
    }

    fn goto_p(&mut self, p: f64) -> &mut Self {
        self.add(GeneratorCommands::PGotoP(p));
        self
    }

    fn goto_m(&mut self, m: f64) -> &mut Self {
        self.add(GeneratorCommands::PGotoM(m));
        self
    }

    fn p_grid_line(&mut self) -> &mut Self {
        self.add(GeneratorCommands::AddGridLineP);
        self
    }

    fn generate_commands(&mut self, pt: &PxuPoint) {
        let consts = pt.consts;
        self.generate_u_grid(consts);

        let p_range = pt.p.re.floor() as i32;

        let max = P_RANGE_MAX - P_RANGE_MIN;

        self.generate_log_cuts(p_range, consts);

        for i in 1..max {
            if p_range - i >= P_RANGE_MIN {
                self.generate_log_cuts(p_range - i, consts);
            }

            if p_range + i <= P_RANGE_MAX {
                self.generate_log_cuts(p_range + i, consts);
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
        self.generate_real_log_cuts(consts);
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

    fn compute_branch_point(&mut self, p_range: i32, branch_point_type: BranchPoint) -> &mut Self {
        self.add(GeneratorCommands::ComputeBranchPoint(
            p_range,
            branch_point_type,
        ));
        self
    }

    fn generate_real_log_cuts(&mut self, consts: CouplingConstants) {
        self.add(GeneratorCommands::AddLogCutXReal);
    }

    fn create_cut(&mut self) -> &mut Self {
        self
    }

    fn push_cut(&mut self) -> &mut Self {
        self
    }

    fn generate_p_cut(
        &mut self,
        p_range: i32,
        branch_point_type: BranchPoint,
        consts: CouplingConstants,
    ) {
    }

    fn generate_log_cuts(&mut self, p_range: i32, consts: CouplingConstants) {
        log::info!("{p_range}");

        let p_start = p_range as f64;
        let k = consts.k() as f64;

        self.compute_branch_point(p_range, BranchPoint::XpPositiveAxisImXmNegative)
            .add(GeneratorCommands::AddCutX(p_range, CutDirection::Positive));

        {
            // Scallion
            let p0 = p_start + 1.0 / 8.0;
            self.p_start_xp(p0)
                .goto_m(-(p_range * consts.k()) as f64)
                .add(GeneratorCommands::AddCutPFull);
            self.p_start_xp(p0)
                .goto_xm(p0, 1.0)
                .goto_m(-(p_range * consts.k()) as f64)
                .add(GeneratorCommands::AddCutPFull);
        }

        {
            // Real positive line
            let p0 = p_start + 1.0 / 8.0;
            self.p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .goto_re(consts.s() * 4.0)
                .add(GeneratorCommands::AddCutPFull);
        }

        self.compute_branch_point(p_range, BranchPoint::XpPositiveAxisImXmPositive)
            .add(GeneratorCommands::AddCutX(p_range, CutDirection::Positive));

        if p_range == 0 {
            let p0 = p_start + 1.0 / 8.0;
            let p1 = p_start + 7.0 / 8.0;

            self.p_start_xp(p0)
                .goto_m(3.0)
                .goto_p(p1)
                .goto_m(0.0)
                .add(GeneratorCommands::AddCutPFull);
        } else if p_range == -1 {
            let p0 = p_start + 1.0 / 8.0;
            let p1 = p_start + 6.0 / 8.0;

            self.p_start_xp(p0)
                .goto_m(-consts.k() as f64 + 1.0)
                .goto_p(p1)
                .goto_m(-consts.k() as f64 + 3.0)
                .goto_p(p0)
                .goto_m(0.0)
                .add(GeneratorCommands::AddCutPFull);
        }

        self.compute_branch_point(
            p_range,
            BranchPoint::XpNegativeAxisFromAboveWithImXmNegative,
        )
        .add(GeneratorCommands::AddCutX(p_range, CutDirection::Positive));

        {
            // Kidney
            let p0 = p_start + 1.0 / 8.0;
            let p1 = p_start + 7.0 / 8.0;
            self.p_start_xp(p0)
                .goto_m(-((p_range + 1) * consts.k()) as f64)
                .add(GeneratorCommands::AddCutPFull);
            self.p_start_xp(p0)
                .goto_xm(p0, 1.0)
                .goto_p(p1)
                .goto_m(-((p_range + 1) * consts.k()) as f64)
                .add(GeneratorCommands::AddCutPFull);
        }

        {
            // Real negative line
            let p0 = p_start + 7.0 / 8.0;
            self.p_start_xp(p0)
                .goto_m(-p_start * k + 1.0)
                .goto_im(0.0)
                .goto_re(-consts.s() * 4.0)
                .add(GeneratorCommands::AddCutPFull);
        }

        self.compute_branch_point(
            p_range,
            BranchPoint::XpNegativeAxisFromBelowWithImXmNegative,
        )
        .add(GeneratorCommands::AddCutX(p_range, CutDirection::Negative));

        self.compute_branch_point(
            p_range,
            BranchPoint::XpNegativeAxisFromAboveWithImXmPositive,
        )
        .add(GeneratorCommands::AddCutX(p_range, CutDirection::Negative));

        self.compute_branch_point(
            p_range,
            BranchPoint::XpNegativeAxisFromBelowWithImXmPositive,
        )
        .add(GeneratorCommands::AddCutX(p_range, CutDirection::Negative));
    }
}

fn get_branch_point(m: f64, consts: CouplingConstants, branch: f64) -> f64 {
    {
        let s = consts.s();
        let u_of_x = |x: C| -> C { x + 1.0 / x - (s - 1.0 / s) * x.ln() };
        let du_dx = |x: C| -> C { (x - s) * (x + 1.0 / s) / (x * x) };
        let sign = if branch == 0.0 { 1.0 } else { -1.0 };
        let guess = if m > 0.0 {
            C::new(0.0, 1.0)
        } else {
            C::new(0.0, 0.01)
        };
        let x = nr::find_root(
            |x| {
                u_of_x(x) - sign * u_of_x(C::from(consts.s()))
                    + branch * C::i() * consts.k() as f64 / consts.h
                    - 2.0 * m * C::i() / consts.h
            },
            du_dx,
            guess,
            1.0e-3,
            10,
        );
        let x = x.unwrap();
        x.arg() / std::f64::consts::PI
    }
}

fn get_branch_point_x(m: f64, consts: CouplingConstants, branch: f64) -> C {
    {
        let s = consts.s();
        let u_of_x = |x: C| -> C { x + 1.0 / x - (s - 1.0 / s) * x.ln() };
        let du_dx = |x: C| -> C { (x - s) * (x + 1.0 / s) / (x * x) };
        let sign = if branch == 0.0 { 1.0 } else { -1.0 };
        let guess = if m > 0.0 {
            C::new(0.0, 1.0)
        } else {
            C::new(0.0, 0.01)
        };
        let x = nr::find_root(
            |x| {
                u_of_x(x) - sign * u_of_x(C::from(consts.s()))
                    + branch * C::i() * consts.k() as f64 / consts.h
                    - 2.0 * m * C::i() / consts.h
            },
            du_dx,
            guess,
            1.0e-3,
            10,
        );
        x.unwrap()
    }
}

#[derive(Debug)]
pub enum CutType {
    U(Component),
    LogX(Component, i32),
    E,
    DebugPath,
}

#[derive(Debug)]
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
            Self::LogBranch(b) => *b == pt.sheet_data.log_branch,
            Self::LogBranchGT(b) => pt.sheet_data.log_branch > *b,
            Self::LogBranchLE(b) => pt.sheet_data.log_branch <= *b,
            Self::LogBranchSum(b) => *b == pt.sheet_data.log_branch_sum,
            Self::EBranch(b) => pt.sheet_data.e_branch == *b,
            Self::UpBranch(b) => pt.sheet_data.u_branch.0 == *b,
            Self::UmBranch(b) => pt.sheet_data.u_branch.1 == *b,
        }
    }
}

#[derive(Debug)]
pub struct Cuts {
    cuts: Vec<Cut>,
    consts: Option<CouplingConstants>,
}

impl Cuts {
    pub fn new() -> Self {
        // let data = HashMap::new();
        let cuts = vec![];
        let consts = None;
        Self { cuts, consts }
    }

    fn populate(&mut self, pt: &PxuPoint) {
        if let Some(consts) = self.consts {
            if consts != pt.consts {
                log::info!("Clearing cuts");
                self.cuts.clear();
            }
        }
        self.consts = Some(pt.consts);

        if self.cuts.is_empty() {
            let start = chrono::Utc::now();
            for p_range in P_RANGE_MIN..=P_RANGE_MAX {
                self.cuts.extend(Cut::get(p_range, pt.consts));
            }
            let end = chrono::Utc::now();
            let elapsed = end - start;

            log::info!(
                "Created {} cuts in {} ms",
                self.cuts.len(),
                elapsed.num_milliseconds()
            );
        }
    }

    pub fn visible(
        &mut self,
        pt: &PxuPoint,
        component: Component,
        populate: bool,
    ) -> impl Iterator<Item = &Cut> {
        if populate {
            self.populate(pt);
        }

        let mut pt = pt.clone();
        pt.u += (pt.sheet_data.log_branch + pt.sheet_data.log_branch_sum) as f64
            * pt.consts.k() as f64
            * C::i()
            / pt.consts.h;

        self.cuts
            .iter()
            .filter(move |c| c.component == component && c.is_visible(&pt))
    }

    pub fn crossed(
        &mut self,
        pt: &PxuPoint,
        component: Component,
        new_value: C,
        populate: bool,
    ) -> impl Iterator<Item = &Cut> {
        if populate {
            self.populate(pt);
        }

        let mut pt = pt.clone();
        pt.u += (pt.sheet_data.log_branch + pt.sheet_data.log_branch_sum) as f64
            * pt.consts.k() as f64
            * C::i()
            / pt.consts.h;

        let new_value = if component == Component::U {
            new_value
                + (pt.sheet_data.log_branch + pt.sheet_data.log_branch_sum) as f64
                    * pt.consts.k() as f64
                    * C::i()
                    / pt.consts.h
        } else {
            new_value
        };

        self.cuts.iter().filter(move |c| {
            c.component == component
                && c.is_visible(&pt)
                && c.intersection(pt.get(component), new_value).is_some()
        })
    }
}

#[derive(Debug)]
pub struct Cut {
    pub component: Component,
    pub paths: Vec<Vec<C>>,
    pub branch_points: Vec<C>,
    pub typ: CutType,
    visibility: Vec<CutVisibilityCondition>,
}

impl Cut {
    fn new(component: Component, paths: Vec<Vec<C>>, branch_points: Vec<C>, typ: CutType) -> Self {
        Self {
            component,
            paths,
            branch_points,
            typ,
            visibility: vec![],
        }
    }

    fn im_xp_positive(mut self) -> Self {
        self.visibility.push(CutVisibilityCondition::ImXp(1));
        self
    }

    fn im_xp_negative(mut self) -> Self {
        self.visibility.push(CutVisibilityCondition::ImXp(-1));
        self
    }

    fn im_xm_positive(mut self) -> Self {
        self.visibility.push(CutVisibilityCondition::ImXm(1));
        self
    }

    fn im_xm_negative(mut self) -> Self {
        self.visibility.push(CutVisibilityCondition::ImXm(-1));
        self
    }

    fn log_branch(mut self, branch: i32) -> Self {
        self.visibility
            .push(CutVisibilityCondition::LogBranch(branch));
        self
    }

    fn log_branch_sum(mut self, branch: i32) -> Self {
        self.visibility
            .push(CutVisibilityCondition::LogBranchSum(branch));
        self
    }

    fn log_branch_gt(mut self, branch: i32) -> Self {
        self.visibility
            .push(CutVisibilityCondition::LogBranchGT(branch));
        self
    }

    fn log_branch_le(mut self, branch: i32) -> Self {
        self.visibility
            .push(CutVisibilityCondition::LogBranchLE(branch));
        self
    }

    fn e_branch(mut self, branch: i32) -> Self {
        self.visibility
            .push(CutVisibilityCondition::EBranch(branch));
        self
    }

    fn up_branch(mut self, branch: i32) -> Self {
        self.visibility
            .push(CutVisibilityCondition::UpBranch(branch));
        self
    }

    fn um_branch(mut self, branch: i32) -> Self {
        self.visibility
            .push(CutVisibilityCondition::UmBranch(branch));
        self
    }

    pub fn get(p_range: i32, consts: CouplingConstants) -> Vec<Cut> {
        let mut cuts = vec![];
        cuts.extend(Self::x_cuts_x(p_range, consts));
        cuts.extend(Self::x_cuts_p(p_range, consts));
        cuts.extend(Self::x_p_cuts_u(p_range, consts));
        cuts.extend(Self::p_cuts_x(p_range, consts));
        cuts.extend(Self::p_cuts_p(p_range, consts));
        cuts.extend(Self::e_cuts(p_range, consts));

        cuts
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

    fn x_cuts_p(p_range: i32, consts: CouplingConstants) -> Vec<Cut> {
        let p_start = p_range as f64;

        let p_s = {
            let p0 = 1.0 / 8.0;
            let p_int = PInterpolator::xp(p0, consts)
                .goto_im(0.0)
                .goto_re(consts.s())
                .clear_path();
            *p_int.p_path.last().unwrap()
        };

        let p_min_one_over_s = {
            let p0 = -1.0 / 64.0;
            let p_int = PInterpolator::xp(p0, consts)
                .goto_im(0.0)
                .goto_re(-1.0 / consts.s())
                .clear_path();
            *p_int.p_path.last().unwrap()
        };

        let mut cuts = vec![];

        let mut p_points = vec![];
        p_points.push(C::from(p_start));
        if p_range != -1 {
            let p0 = p_start + 0.25;
            let p_int = PInterpolator::xp(p0, consts)
                .goto_xm(p0, 1.0)
                .goto_xm(p0, 0.0)
                .clear_path();

            let p_int2 = p_int.clone().goto_xm(p_start + 127.0 / 128.0, 0.0);
            p_points.extend(p_int2.p_path.into_iter().rev());

            let p_int2 = p_int.clone().goto_xm(p_start + 1.0 / 128.0, 0.0);
            p_points.extend(p_int2.p_path);

            if p_range == 0 {
                p_points.push(p_s);
            } else {
                p_points.push(C::from(p_start));
            }
        }

        {
            let p0 = p_start + 1.0 / 8.0;
            let p_int = PInterpolator::xp(p0, consts).goto_xp(p0, 0.0);
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p_start + 1.0 / 16.0, 0.0);
            p_points.extend(p_int2.p_path.into_iter().rev());

            let p_int2 = p_int.goto_xp(p_start + 1.0 - 1.0 / 16.0, 0.0);
            p_points.extend(p_int2.p_path);

            if p_range != -1 {
                p_points.push(C::from(p_start + 1.0));
            } else {
                p_points.push(p_min_one_over_s);
            }
        }

        cuts.push(
            Cut::new(
                Component::P,
                vec![p_points.clone()],
                vec![],
                CutType::U(Component::Xp),
            )
            .e_branch(1),
        );

        cuts.push(
            Cut::new(
                Component::P,
                vec![p_points.iter().map(|x| x.conj()).collect()],
                vec![],
                CutType::U(Component::Xm),
            )
            .e_branch(1),
        );

        if p_range == -1 {
            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.clone()],
                    vec![],
                    CutType::U(Component::Xm),
                )
                .e_branch(-1),
            );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.iter().map(|x| x.conj()).collect()],
                    vec![],
                    CutType::U(Component::Xp),
                )
                .e_branch(-1),
            );
        }

        if p_range == 0 {
            let mut p_points: Vec<C> = vec![];
            // let mut x_cuts = vec![];

            let p0 = p_range as f64 + 1.0 / 8.0;
            let p2 = p_range as f64 + 7.0 / 8.0;

            let p_int = PInterpolator::xp(p0, consts)
                .goto_xm(p0, 1.0)
                .goto_xm(p2, 1.0)
                .goto_xm(p2, -2.0);

            // x_cuts.push(p_int.x_path.clone());

            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xm(p_start + 0.01, -2.0);
            // x_cuts.push(p_int2.x_path.clone());
            p_points.extend(p_int2.p_path.into_iter().rev());

            let p_int2 = p_int.clone().goto_xm(p_start + 0.99, -2.0);
            // x_cuts.push(p_int2.x_path.clone());
            p_points.extend(p_int2.p_path);

            // cuts.push(
            //     Cut::new(Component::Xp, x_cuts, vec![], CutType::DebugPath).log_branch(p_range).e_branch(1),
            // );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.iter().map(|p| p.conj()).collect()],
                    vec![],
                    CutType::U(Component::Xp),
                )
                .e_branch(1),
            );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.clone()],
                    vec![],
                    CutType::U(Component::Xm),
                )
                .e_branch(1),
            );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.iter().map(|p| p.conj()).collect()],
                    vec![],
                    CutType::U(Component::Xm),
                )
                .e_branch(-1),
            );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.clone()],
                    vec![],
                    CutType::U(Component::Xp),
                )
                .e_branch(-1),
            );

            let mut p_points: Vec<C> = vec![];
            let mut x_cuts = vec![];

            let p_int = PInterpolator::xp(p0, consts)
                .goto_xm(p0, 1.0)
                .goto_xm(p2, 1.0)
                .goto_xm(p2, -(consts.k() as f64));

            x_cuts.push(p_int.x_path.clone());

            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xm(p_start + 0.01, -(consts.k() as f64));
            x_cuts.push(p_int2.x_path.clone());
            p_points.extend(p_int2.p_path.into_iter().rev());

            let p_int2 = p_int.clone().goto_xm(p_start + 0.99, -(consts.k() as f64));
            x_cuts.push(p_int2.x_path.clone());
            p_points.extend(p_int2.p_path);

            p_points.push(p_min_one_over_s);

            // cuts.push(
            //     Cut::new(Component::Xp, x_cuts, vec![], CutType::DebugPath)
            //         .log_branch(p_range)
            //         .e_branch(1),
            // );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.iter().map(|p| p.conj()).collect()],
                    vec![],
                    CutType::U(Component::Xm),
                )
                .e_branch(1),
            );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.clone()],
                    vec![],
                    CutType::U(Component::Xp),
                )
                .e_branch(1),
            );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.iter().map(|p| p.conj()).collect()],
                    vec![],
                    CutType::U(Component::Xp),
                )
                .e_branch(-1),
            );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.clone()],
                    vec![],
                    CutType::U(Component::Xm),
                )
                .e_branch(-1),
            );
        }

        if p_range == -1 {
            let mut p_points: Vec<C> = vec![];
            // let mut x_cuts = vec![];

            let p0 = p_range as f64 + 1.0 / 8.0;

            // let x = get_branch_point_x(m - 1.0 - consts.k() as f64, consts, -1.0);

            let p_int = PInterpolator::xp(p0, consts).goto_xp(p0, consts.k() as f64);
            // .goto_xm(p0, consts.k() as f64)
            // .goto_xm(p0, 0.0);

            // x_cuts.push(p_int.x_path.clone());

            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p_start + 0.9999, consts.k() as f64);
            // x_cuts.push(p_int2.x_path.clone());
            p_points.extend(p_int2.p_path.into_iter().rev());

            let p_int2 = p_int.clone().goto_xp(p_start + 0.01, consts.k() as f64);
            // x_cuts.push(p_int2.x_path.clone());
            p_points.extend(p_int2.p_path);

            let p_int = PInterpolator::xp(p0, consts)
                .goto_xp(p0, consts.k() as f64)
                .goto_xm(p0, consts.k() as f64);

            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xm(p_start + 0.01, consts.k() as f64);
            // x_cuts.push(p_int2.x_path.clone());
            p_points.extend(p_int2.p_path.into_iter().rev());

            let p_int2 = p_int.clone().goto_xm(p_start + 0.9999, consts.k() as f64);
            // x_cuts.push(p_int2.x_path.clone());
            p_points.extend(p_int2.p_path);

            // cuts.push(
            //     Cut::new(Component::Xp, x_cuts, vec![x], CutType::DebugPath)
            //         .log_branch(p_range).e_branch(1),
            // );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.iter().map(|p| p.conj()).collect()],
                    vec![],
                    CutType::U(Component::Xp),
                )
                .e_branch(-1),
            );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.clone()],
                    vec![],
                    CutType::U(Component::Xm),
                )
                .e_branch(-1),
            );
        }

        if p_range > 0 {
            // let mut x_cuts = vec![];

            let p0 = p_range as f64 + 1.0 / 8.0;
            let p2 = p_range as f64 + 7.0 / 8.0;

            let mut p_points: Vec<C> = vec![];
            let m = -2.0 - 2.0 * p_range as f64 * consts.k() as f64;

            // let x = get_branch_point_x(m - 1.0 - consts.k() as f64, consts, -1.0);

            let p_int = PInterpolator::xp(p0, consts)
                .goto_xm(p0, 1.0)
                .goto_xm(p2, 1.0)
                .goto_xm(p2, m);

            // x_cuts.push(p_int.x_path.clone());

            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xm(p_start + 0.01, m);
            // x_cuts.push(p_int2.x_path.clone());
            p_points.extend(p_int2.p_path.into_iter().rev());

            let p_int2 = p_int.clone().goto_xm(p_start + 0.9999, m);
            // x_cuts.push(p_int2.x_path.clone());
            p_points.extend(p_int2.p_path);

            // cuts.push(
            //     Cut::new(Component::Xp, x_cuts, vec![], CutType::DebugPath).log_branch(p_range).e_branch(1),
            // );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.iter().map(|p| p.conj()).collect()],
                    vec![],
                    CutType::U(Component::Xp),
                )
                .e_branch(1),
            );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.clone()],
                    vec![],
                    CutType::U(Component::Xm),
                )
                .e_branch(1),
            );
        }

        if p_range <= -1 {
            let mut p_points: Vec<C> = vec![];
            // let mut x_cuts = vec![];

            let p0 = p_range as f64 + 1.0 / 32.0;
            let p2 = p_range as f64 + 1.0 / 8.0;

            let m = -2.0 - 2.0 * p_range as f64 * consts.k() as f64;
            let m1 = 2.0 - 1.0 * p_range as f64 * consts.k() as f64;

            // let x = get_branch_point_x(m - 1.0 - consts.k() as f64, consts, -1.0);

            let p_int = PInterpolator::xp(p2, consts)
                .goto_xm(p2, 1.0)
                .goto_xm(p2, m1)
                .goto_xm(p0, m1)
                .goto_xm(p0, m);

            // x_cuts.push(p_int.x_path.clone());

            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xm(p_start + 0.01, m);
            // x_cuts.push(p_int2.x_path.clone());
            p_points.extend(p_int2.p_path.into_iter().rev());

            let p_int2 = p_int.clone().goto_xm(p_start + 0.99, m);
            // x_cuts.push(p_int2.x_path.clone());
            p_points.extend(p_int2.p_path);

            // cuts.push(
            //     Cut::new(Component::Xp, x_cuts, vec![], CutType::DebugPath).log_branch(p_range).e_branch(1),
            // );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.iter().map(|p| p.conj()).collect()],
                    vec![],
                    CutType::U(Component::Xp),
                )
                .e_branch(1),
            );

            cuts.push(
                Cut::new(
                    Component::P,
                    vec![p_points.clone()],
                    vec![],
                    CutType::U(Component::Xm),
                )
                .e_branch(1),
            );
        }

        cuts
    }

    fn x_cuts_x(p_range: i32, consts: CouplingConstants) -> Vec<Cut> {
        let mut cuts = vec![];

        let branch_points = if p_range == 0 {
            vec![C::from(consts.s())]
        } else if p_range == -1 {
            vec![C::from(-1.0 / consts.s())]
        } else {
            vec![]
        };

        cuts.push(
            Cut::new(
                Component::Xp,
                vec![XInterpolator::generate_xp_full(p_range, 0.0, consts)],
                branch_points.clone(),
                CutType::U(Component::Xp),
            )
            .log_branch(p_range), // .e_branch(1),
        );

        cuts.push(
            Cut::new(
                Component::Xm,
                vec![XInterpolator::generate_xm_full(p_range, 0.0, consts)],
                branch_points.clone(),
                CutType::U(Component::Xm),
            )
            .log_branch(p_range), // .e_branch(1),
        );

        cuts.push(
            Cut::new(
                Component::Xp,
                vec![XInterpolator::generate_xm_full(p_range, 0.0, consts)],
                branch_points.clone(),
                CutType::U(Component::Xp),
            )
            .log_branch(p_range), // .e_branch(1),
        );

        cuts.push(
            Cut::new(
                Component::Xm,
                vec![XInterpolator::generate_xp_full(p_range, 0.0, consts)],
                branch_points.clone(),
                CutType::U(Component::Xm),
            )
            .log_branch(p_range), // .e_branch(1),
        );

        if p_range == 0 {
            let ps = get_branch_point(1.0, consts, 0.0);

            let paths = vec![XInterpolator::generate_xp(ps, 1.0, 2.0, consts)];
            let branch_points = vec![*paths[0].first().unwrap()];

            cuts.push(
                Cut::new(
                    Component::Xp,
                    paths,
                    branch_points,
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range), // .e_branch(1),
            );

            let paths = vec![XInterpolator::generate_xm(ps, 1.0, 2.0, consts)];
            let branch_points = vec![*paths[0].first().unwrap()];

            cuts.push(
                Cut::new(
                    Component::Xm,
                    paths,
                    branch_points,
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range), // .e_branch(1),
            );

            let ps = get_branch_point(-1.0, consts, 0.0);

            let paths = vec![XInterpolator::generate_xm(ps, 1.0, -2.0, consts)];
            let branch_points = vec![*paths[0].first().unwrap()];

            cuts.push(
                Cut::new(
                    Component::Xp,
                    paths,
                    branch_points,
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range), // .e_branch(1),
            );

            let paths = vec![XInterpolator::generate_xp(ps, 1.0, -2.0, consts)];
            let branch_points = vec![*paths[0].first().unwrap()];

            cuts.push(
                Cut::new(
                    Component::Xm,
                    paths,
                    branch_points,
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range), // .e_branch(1),
            );

            let paths = vec![XInterpolator::generate_xm(
                0.0,
                1.0,
                -(consts.k() as f64),
                consts,
            )];
            let branch_points = vec![];

            cuts.push(
                Cut::new(
                    Component::Xp,
                    paths,
                    branch_points,
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range), // .e_branch(1),
            );

            let paths = vec![XInterpolator::generate_xp(
                0.0,
                1.0,
                -(consts.k() as f64),
                consts,
            )];
            let branch_points = vec![];

            cuts.push(
                Cut::new(
                    Component::Xm,
                    paths,
                    branch_points,
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range), // .e_branch(1),
            );

            let m = 1.0 - p_range as f64 * consts.k() as f64;
            let p_minus_one_over_s = get_branch_point(m, consts, 1.0);
            // let x = get_branch_point_x(m, consts, 1.0);
            let m = 2.0 - (p_range + 1) as f64 * consts.k() as f64;

            let paths = vec![XInterpolator::generate_xp(
                p_minus_one_over_s.floor(),
                p_minus_one_over_s,
                m,
                consts,
            )];
            let branch_points = vec![*paths[0].last().unwrap()];

            cuts.push(
                Cut::new(
                    Component::Xp,
                    paths.clone(),
                    branch_points.clone(),
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range)
                .im_xm_positive(), // .e_branch(1),
            );

            let paths = vec![XInterpolator::generate_xm(
                p_minus_one_over_s.floor(),
                p_minus_one_over_s,
                m,
                consts,
            )];
            let branch_points = vec![*paths[0].last().unwrap()];

            cuts.push(
                Cut::new(
                    Component::Xm,
                    paths.clone(),
                    branch_points.clone(),
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range)
                .im_xp_negative(), // .e_branch(1),
            );

            let m = -1.0 + (p_range + 1) as f64 * consts.k() as f64;
            let p_minus_one_over_s = get_branch_point(m, consts, 1.0);
            let x = get_branch_point_x(m, consts, 1.0);
            let m = -2.0 + (p_range + 1) as f64 * consts.k() as f64;

            let paths = vec![XInterpolator::generate_xp(
                p_minus_one_over_s.floor(),
                p_minus_one_over_s,
                m,
                consts,
            )];
            let branch_points = vec![*paths[0].last().unwrap()];

            cuts.push(
                Cut::new(
                    Component::Xm,
                    paths.clone(),
                    branch_points.clone(),
                    // vec![x],
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range)
                .im_xp_negative(), // .e_branch(1),
            );

            let paths = vec![XInterpolator::generate_xm(
                p_minus_one_over_s.floor(),
                p_minus_one_over_s,
                m,
                consts,
            )];
            let branch_points = vec![*paths[0].last().unwrap()];

            cuts.push(
                Cut::new(
                    Component::Xp,
                    paths.clone(),
                    branch_points.clone(),
                    // vec![x],
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range)
                .im_xm_positive(), // .e_branch(1),
            );
        } else if p_range == -1 {
            let ps = get_branch_point(1.0, consts, 1.0);

            let paths = vec![XInterpolator::generate_xp(
                ps,
                0.0,
                2.0 - consts.k() as f64,
                consts,
            )];
            let branch_points = vec![*paths[0].last().unwrap()];

            cuts.push(
                Cut::new(
                    Component::Xp,
                    paths,
                    branch_points,
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range)
                //.e_branch(1)
                .im_xm_negative(),
            );

            let paths = vec![XInterpolator::generate_xm(
                ps,
                0.0,
                2.0 - consts.k() as f64,
                consts,
            )];
            let branch_points = vec![*paths[0].last().unwrap()];

            cuts.push(
                Cut::new(
                    Component::Xm,
                    paths,
                    branch_points,
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range)
                //.e_branch(1)
                .im_xp_positive(),
            );

            let m =
                1.0 + if p_range < 0 { p_range - 1 } else { p_range } as f64 * consts.k() as f64;

            let p_minus_one_over_s = get_branch_point(m, consts, m.signum());
            let m = 2.0 + (2 * p_range - 1) as f64 * consts.k() as f64;

            let paths = vec![XInterpolator::generate_xp(
                p_minus_one_over_s.floor(),
                p_minus_one_over_s,
                m,
                consts,
            )];
            let branch_points = vec![paths[0][0]];

            cuts.push(
                Cut::new(
                    Component::Xp,
                    paths,
                    branch_points,
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range)
                //.e_branch(1)
                .im_xm_positive(),
            );

            let paths = vec![XInterpolator::generate_xm(
                p_minus_one_over_s.floor(),
                p_minus_one_over_s,
                m,
                consts,
            )];
            let branch_points = vec![paths[0][0]];

            cuts.push(
                Cut::new(
                    Component::Xm,
                    paths,
                    branch_points,
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range)
                //.e_branch(1)
                .im_xp_negative(),
            );
        } else {
            let paths = vec![XInterpolator::generate_xp(
                0.0,
                1.0,
                2.0 + p_range as f64 * consts.k() as f64,
                consts,
            )];
            let branch_points = vec![];

            cuts.push(
                Cut::new(
                    Component::Xp,
                    paths,
                    branch_points,
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range)
                //.e_branch(1)
                .im_xm_negative(),
            );

            let paths = vec![XInterpolator::generate_xm(
                0.0,
                1.0,
                2.0 + p_range as f64 * consts.k() as f64,
                consts,
            )];
            let branch_points = vec![];

            cuts.push(
                Cut::new(
                    Component::Xm,
                    paths,
                    branch_points,
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range)
                //.e_branch(1)
                .im_xp_positive(),
            );
        }

        if p_range > 0 {
            let paths = vec![XInterpolator::generate_xp(
                0.0,
                1.0,
                2.0 + (3 * p_range) as f64 * consts.k() as f64,
                consts,
            )];
            let branch_points = vec![];

            cuts.push(
                Cut::new(
                    Component::Xp,
                    paths,
                    branch_points,
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range)
                //.e_branch(1)
                .im_xm_positive(),
            );

            let paths = vec![XInterpolator::generate_xm(
                0.0,
                1.0,
                -2.0 - p_range as f64 * consts.k() as f64,
                consts,
            )];
            let branch_points = vec![];

            cuts.push(
                Cut::new(
                    Component::Xp,
                    paths,
                    branch_points,
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range)
                .im_xm_negative(), //.e_branch(1),
            );

            let paths = vec![XInterpolator::generate_xm(
                0.0,
                1.0,
                2.0 + (3 * p_range) as f64 * consts.k() as f64,
                consts,
            )];
            let branch_points = vec![];

            cuts.push(
                Cut::new(
                    Component::Xm,
                    paths,
                    branch_points,
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range)
                //.e_branch(1)
                .im_xp_negative(),
            );

            let paths = vec![XInterpolator::generate_xp(
                0.0,
                1.0,
                -2.0 - p_range as f64 * consts.k() as f64,
                consts,
            )];
            let branch_points = vec![];

            cuts.push(
                Cut::new(
                    Component::Xm,
                    paths,
                    branch_points,
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range)
                .im_xp_positive(), //.e_branch(1),
            );
        }

        if p_range < -1 {
            let paths = vec![XInterpolator::generate_xp(
                0.0,
                1.0,
                2.0 + 3.0 * p_range as f64 * consts.k() as f64,
                consts,
            )];
            let branch_points = vec![];

            cuts.push(
                Cut::new(
                    Component::Xp,
                    paths,
                    branch_points,
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range)
                //.e_branch(1)
                .im_xm_positive(),
            );

            let paths = vec![XInterpolator::generate_xm(
                0.0,
                1.0,
                2.0 + 3.0 * p_range as f64 * consts.k() as f64,
                consts,
            )];
            let branch_points = vec![];

            cuts.push(
                Cut::new(
                    Component::Xm,
                    paths,
                    branch_points,
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range)
                //.e_branch(1)
                .im_xp_negative(),
            );
        }

        if p_range <= -1 {
            let paths = vec![XInterpolator::generate_xm(
                0.0,
                1.0,
                -2.0 - p_range as f64 * consts.k() as f64,
                consts,
            )];
            let branch_points = vec![];

            cuts.push(
                Cut::new(
                    Component::Xp,
                    paths,
                    branch_points,
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range), //.e_branch(1),
            );

            let paths = vec![XInterpolator::generate_xp(
                0.0,
                1.0,
                -2.0 - p_range as f64 * consts.k() as f64,
                consts,
            )];
            let branch_points = vec![];

            cuts.push(
                Cut::new(
                    Component::Xm,
                    paths,
                    branch_points,
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range), //.e_branch(1),
            );
        }

        cuts
    }

    fn x_p_cuts_u(p_range: i32, consts: CouplingConstants) -> Vec<Cut> {
        let s = consts.s();
        const INFINITY: f64 = 1000.0;

        let us = s + 1.0 / s - (s - 1.0 / s) * s.ln();

        let k = consts.k() as f64;

        let shift = p_range as f64 * k / consts.h;

        let mut cuts = vec![];

        cuts.push(
            Cut::new(
                Component::U,
                vec![vec![
                    C::new(INFINITY, -(1.0 + p_range as f64 * k) / consts.h + shift),
                    C::new(us, -(1.0 + p_range as f64 * k) / consts.h + shift),
                ]],
                vec![C::new(us, -(1.0 + p_range as f64 * k) / consts.h + shift)],
                CutType::LogX(Component::Xp, 0),
            )
            .log_branch(p_range),
        );

        cuts.push(
            Cut::new(
                Component::U,
                vec![vec![
                    C::new(INFINITY, (1.0 + p_range as f64 * k) / consts.h + shift),
                    C::new(us, (1.0 + p_range as f64 * k) / consts.h + shift),
                ]],
                vec![C::new(us, (1.0 + p_range as f64 * k) / consts.h + shift)],
                CutType::LogX(Component::Xm, 0),
            )
            .log_branch(p_range),
        );

        if p_range == 0 {
            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(-INFINITY, 1.0 / consts.h + shift),
                        C::new(us, 1.0 / consts.h + shift),
                    ]],
                    vec![C::new(us, 1.0 / consts.h + shift)],
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range),
            );

            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(-INFINITY, -1.0 / consts.h + shift),
                        C::new(us, -1.0 / consts.h + shift),
                    ]],
                    vec![C::new(us, -1.0 / consts.h + shift)],
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range),
            );

            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(INFINITY, (1.0 - k) / consts.h + shift),
                        C::new(-us, (1.0 - k) / consts.h + shift),
                    ]],
                    vec![C::new(-us, (1.0 - k) / consts.h + shift)],
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range)
                .im_xm_positive(),
            );

            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(INFINITY, (-1.0 + k) / consts.h + shift),
                        C::new(-us, (-1.0 + k) / consts.h + shift),
                    ]],
                    vec![C::new(-us, (-1.0 + k) / consts.h + shift)],
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range)
                .im_xp_negative(),
            );
        } else if p_range == -1 {
            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(-us, 1.0 / consts.h + shift),
                        C::new(INFINITY, 1.0 / consts.h + shift),
                    ]],
                    vec![C::new(-us, 1.0 / consts.h + shift)],
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range),
            );

            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(-us, -1.0 / consts.h + shift),
                        C::new(INFINITY, -1.0 / consts.h + shift),
                    ]],
                    vec![C::new(-us, -1.0 / consts.h + shift)],
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range),
            );
        } else {
            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(-INFINITY, 1.0 / consts.h + shift),
                        C::new(INFINITY, 1.0 / consts.h + shift),
                    ]],
                    vec![],
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range)
                .im_xm_negative(),
            );

            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(-INFINITY, -1.0 / consts.h + shift),
                        C::new(INFINITY, -1.0 / consts.h + shift),
                    ]],
                    vec![],
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range)
                .im_xp_positive(),
            );

            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(
                            -INFINITY,
                            (1.0 + (2 * p_range) as f64 * k) / consts.h + shift,
                        ),
                        C::new(
                            INFINITY,
                            (1.0 + (2 * p_range) as f64 * k) / consts.h + shift,
                        ),
                    ]],
                    vec![],
                    CutType::U(Component::Xm),
                )
                .log_branch(p_range)
                .im_xm_positive(),
            );

            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(
                            -INFINITY,
                            -(1.0 + (2 * p_range) as f64 * k) / consts.h + shift,
                        ),
                        C::new(
                            INFINITY,
                            -(1.0 + (2 * p_range) as f64 * k) / consts.h + shift,
                        ),
                    ]],
                    vec![],
                    CutType::U(Component::Xp),
                )
                .log_branch(p_range)
                .im_xp_negative(),
            );
        }

        {
            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(
                            -INFINITY,
                            -(1.0 + (p_range + 1) as f64 * k) / consts.h + shift,
                        ),
                        C::new(-us, -(1.0 + (p_range + 1) as f64 * k) / consts.h + shift),
                    ]],
                    vec![C::new(
                        -us,
                        -(1.0 + (p_range + 1) as f64 * k) / consts.h + shift,
                    )],
                    CutType::LogX(Component::Xp, 1),
                )
                .log_branch(p_range)
                .im_xp_positive(),
            );

            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(
                            -INFINITY,
                            (1.0 + (p_range + 1) as f64 * k) / consts.h + shift,
                        ),
                        C::new(-us, (1.0 + (p_range + 1) as f64 * k) / consts.h + shift),
                    ]],
                    vec![C::new(
                        -us,
                        (1.0 + (p_range + 1) as f64 * k) / consts.h + shift,
                    )],
                    CutType::LogX(Component::Xm, 1),
                )
                .log_branch(p_range)
                .im_xm_negative(),
            );

            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(
                            -INFINITY,
                            -(1.0 + (p_range - 1) as f64 * k) / consts.h + shift,
                        ),
                        C::new(-us, -(1.0 + (p_range - 1) as f64 * k) / consts.h + shift),
                    ]],
                    vec![C::new(
                        -us,
                        -(1.0 + (p_range - 1) as f64 * k) / consts.h + shift,
                    )],
                    CutType::LogX(Component::Xp, -1),
                )
                .log_branch(p_range)
                .im_xp_negative(),
            );

            cuts.push(
                Cut::new(
                    Component::U,
                    vec![vec![
                        C::new(
                            -INFINITY,
                            (1.0 + (p_range - 1) as f64 * k) / consts.h + shift,
                        ),
                        C::new(-us, (1.0 + (p_range - 1) as f64 * k) / consts.h + shift),
                    ]],
                    vec![C::new(
                        -us,
                        (1.0 + (p_range - 1) as f64 * k) / consts.h + shift,
                    )],
                    CutType::LogX(Component::Xm, -1),
                )
                .log_branch(p_range)
                .im_xm_positive(),
            );
        }

        cuts
    }

    fn p_cuts_x(p_range: i32, consts: CouplingConstants) -> Vec<Cut> {
        let mut cuts = vec![];

        // xp negative axis from above
        let paths = vec![vec![C::from(-100.0), C::zero()]];
        let branch_points = vec![C::from(-1.0 / consts.s()), C::zero()];

        cuts.push(
            Cut::new(
                Component::Xp,
                paths,
                branch_points,
                CutType::LogX(Component::Xp, 1),
            )
            .log_branch(p_range)
            .im_xp_positive(),
        );

        // xp negative axis from below
        let paths = vec![vec![C::from(-100.0), C::zero()]];
        let branch_points = vec![C::from(-1.0 / consts.s()), C::zero()];

        cuts.push(
            Cut::new(
                Component::Xp,
                paths,
                branch_points,
                CutType::LogX(Component::Xp, -1),
            )
            .log_branch(p_range)
            .im_xp_negative(),
        );

        // xm negative axis from below
        let paths = vec![vec![C::from(-100.0), C::zero()]];
        let branch_points = vec![C::from(-1.0 / consts.s()), C::zero()];

        cuts.push(
            Cut::new(
                Component::Xm,
                paths,
                branch_points,
                CutType::LogX(Component::Xm, 1),
            )
            .log_branch(p_range)
            .im_xm_negative(),
        );

        // xm negative axis from above
        let paths = vec![vec![C::from(-100.0), C::zero()]];
        let branch_points = vec![C::from(-1.0 / consts.s()), C::zero()];

        cuts.push(
            Cut::new(
                Component::Xm,
                paths,
                branch_points,
                CutType::LogX(Component::Xm, -1),
            )
            .log_branch(p_range)
            .im_xm_positive(),
        );

        // xp image of xm negative real axis from below

        let m = 1.0 + if p_range < 0 { p_range } else { p_range + 1 } as f64 * consts.k() as f64;

        let p_minus_one_over_s = get_branch_point(m, consts, m.signum());
        let m = 2.0 + (2 * p_range + 1) as f64 * consts.k() as f64;

        let paths = vec![XInterpolator::generate_xp(
            p_minus_one_over_s,
            p_minus_one_over_s.ceil(),
            m,
            consts,
        )];
        let branch_points = vec![paths[0][0]];

        cuts.push(
            Cut::new(
                Component::Xp,
                paths,
                branch_points,
                CutType::LogX(Component::Xm, 1),
            )
            .log_branch(p_range)
            // .e_branch(1)
            .im_xm_negative(),
        );

        // xm image of xp negative real axis from above

        let paths = vec![XInterpolator::generate_xm(
            p_minus_one_over_s,
            p_minus_one_over_s.ceil(),
            m,
            consts,
        )];
        let branch_points = vec![paths[0][0]];

        cuts.push(
            Cut::new(
                Component::Xm,
                paths,
                branch_points,
                CutType::LogX(Component::Xp, 1),
            )
            .log_branch(p_range)
            // .e_branch(1)
            .im_xp_positive(),
        );

        // xp image of xm negative real axis from above

        let m = 1.0 + if p_range < 0 { p_range - 1 } else { p_range } as f64 * consts.k() as f64;

        let p_minus_one_over_s = get_branch_point(m, consts, m.signum());
        let m = 2.0 + (2 * p_range - 1) as f64 * consts.k() as f64;

        let paths = vec![XInterpolator::generate_xp(
            p_minus_one_over_s,
            p_minus_one_over_s.ceil(),
            m,
            consts,
        )];
        let branch_points = vec![paths[0][0]];

        cuts.push(
            Cut::new(
                Component::Xp,
                paths,
                branch_points,
                CutType::LogX(Component::Xm, -1),
            )
            .log_branch(p_range)
            // .e_branch(1)
            .im_xm_positive(),
        );

        // xm image of xp negative real axis from below

        let paths = vec![XInterpolator::generate_xm(
            p_minus_one_over_s,
            p_minus_one_over_s.ceil(),
            m,
            consts,
        )];
        let branch_points = vec![paths[0][0]];

        cuts.push(
            Cut::new(
                Component::Xm,
                paths,
                branch_points,
                CutType::LogX(Component::Xp, -1),
            )
            .log_branch(p_range)
            // .e_branch(1)
            .im_xp_negative(),
        );

        // xp real positive axis

        let paths = vec![vec![C::zero(), C::from(100.0)]];
        let branch_points = vec![C::from(consts.s()), C::zero()];

        cuts.push(
            Cut::new(
                Component::Xp,
                paths,
                branch_points,
                CutType::LogX(Component::Xp, 0),
            )
            .log_branch(p_range),
        );

        // xm real positive axis

        let paths = vec![vec![C::zero(), C::from(100.0)]];
        let branch_points = vec![C::from(consts.s()), C::zero()];

        cuts.push(
            Cut::new(
                Component::Xm,
                paths,
                branch_points,
                CutType::LogX(Component::Xm, 0),
            )
            .log_branch(p_range),
        );

        // xp image of xm real positive axis

        let m = 2.0 + (2 * p_range) as f64 * consts.k() as f64;
        let p_s = get_branch_point(m / 2.0, consts, 0.0);

        let paths = vec![XInterpolator::generate_xp(p_s.floor(), p_s, m, consts)];
        let branch_points = vec![*paths[0].last().unwrap()];

        cuts.push(
            Cut::new(
                Component::Xp,
                paths,
                branch_points,
                CutType::LogX(Component::Xm, 0),
            )
            .log_branch(p_range), // .e_branch(1),
        );

        // xm image of xp real positive axis

        let paths = vec![XInterpolator::generate_xm(p_s.floor(), p_s, m, consts)];
        let branch_points = vec![*paths[0].last().unwrap()];

        cuts.push(
            Cut::new(
                Component::Xm,
                paths.clone(),
                branch_points.clone(),
                CutType::LogX(Component::Xp, 0),
            )
            .log_branch(p_range), // .e_branch(1),
        );

        // xp image of xm real positive axis

        let m = -2.0 - (2 * (p_range)) as f64 * consts.k() as f64;
        let p_s = get_branch_point(m / 2.0, consts, 0.0);

        let paths = vec![XInterpolator::generate_xm(p_s.floor(), p_s, m, consts)];
        let branch_points = vec![*paths[0].last().unwrap()];

        cuts.push(
            Cut::new(
                Component::Xp,
                paths.clone(),
                branch_points.clone(),
                CutType::LogX(Component::Xm, 0),
            )
            .log_branch(p_range), // .e_branch(1),
        );

        // xm image of xp real positive axis

        let paths = vec![XInterpolator::generate_xp(p_s.floor(), p_s, m, consts)];
        let branch_points = vec![*paths[0].last().unwrap()];

        cuts.push(
            Cut::new(
                Component::Xm,
                paths.clone(),
                branch_points.clone(),
                CutType::LogX(Component::Xp, 0),
            )
            .log_branch(p_range), // .e_branch(1),
        );

        // xp image of xm real negative axis

        let m = -1.0 - (p_range) as f64 * consts.k() as f64;
        let p_minus_one_over_s = get_branch_point(m, consts, 1.0);
        let x = get_branch_point_x(m, consts, 1.0);
        let m = -2.0 - (2 * p_range + 1) as f64 * consts.k() as f64;

        let paths = vec![XInterpolator::generate_xm(
            p_minus_one_over_s,
            p_minus_one_over_s.ceil(),
            m,
            consts,
        )];
        let branch_points = vec![*paths[0].first().unwrap()];

        cuts.push(
            Cut::new(
                Component::Xp,
                paths.clone(),
                branch_points.clone(),
                CutType::LogX(Component::Xm, 1),
            )
            .log_branch(p_range)
            .im_xm_negative(), // .e_branch(1),
        );

        // xm image of xp real negative axis

        let paths = vec![XInterpolator::generate_xp(
            p_minus_one_over_s,
            p_minus_one_over_s.ceil(),
            m,
            consts,
        )];
        let branch_points = vec![*paths[0].first().unwrap()];

        cuts.push(
            Cut::new(
                Component::Xm,
                paths.clone(),
                branch_points.clone(),
                CutType::LogX(Component::Xp, 1),
            )
            .log_branch(p_range)
            .im_xp_positive(), // .e_branch(1),
        );

        // xp image of xm real negative axis

        let m = -1.0 - (p_range - 1) as f64 * consts.k() as f64;
        let p_minus_one_over_s = get_branch_point(m, consts, 1.0);
        let m = -2.0 - (2 * p_range - 1) as f64 * consts.k() as f64;

        let paths = vec![XInterpolator::generate_xm(
            p_minus_one_over_s,
            p_minus_one_over_s.ceil(),
            m,
            consts,
        )];
        let branch_points = vec![*paths[0].first().unwrap()];

        cuts.push(
            Cut::new(
                Component::Xp,
                paths.clone(),
                branch_points.clone(),
                CutType::LogX(Component::Xm, -1),
            )
            .log_branch(p_range)
            .im_xm_positive(), // .e_branch(1),
        );

        // xm image of xp real negative axis

        let paths = vec![XInterpolator::generate_xp(
            p_minus_one_over_s,
            p_minus_one_over_s.ceil(),
            m,
            consts,
        )];
        let branch_points = vec![*paths[0].first().unwrap()];

        cuts.push(
            Cut::new(
                Component::Xm,
                paths.clone(),
                branch_points.clone(),
                CutType::LogX(Component::Xp, -1),
            )
            .log_branch(p_range)
            .im_xp_negative(), // .e_branch(1),
        );

        cuts
    }

    fn p_cuts_p(p_range: i32, consts: CouplingConstants) -> Vec<Cut> {
        let mut cuts = vec![];
        let p_start = p_range as f64;

        if p_range == 0 {
            // Real negative axis

            let m = consts.k() as f64 + 2.0;
            let p1 = get_branch_point(m - 1.0, consts, 1.0);

            let mut p_points = vec![];
            // let mut x_cuts = vec![];

            let p0 = 1.0 / 8.0;
            let p2 = 7.0 / 8.0;

            let p_int = PInterpolator::xp(p2, consts)
                .goto_xp(p2, 3.0)
                .goto_xp(p0, 3.0)
                .goto_xp(p0, m + 1.0)
                .goto_xp(p2, m + 1.0)
                .goto_xp(p2, m);

            // x_cuts.push(p_int.x_path.clone());

            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(0.99, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p1, m);
            let branch_point = *p_int2.p_path.last().unwrap();
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            let p_int = PInterpolator::xp(p2, consts).goto_xp(p2, m);
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p1, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(0.99, m);
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            // cuts.push(Cut::new(
            //     Component::Xp,
            //     x_cuts,
            //      vec![],
            //      CutType::DebugPath,
            //
            // ));

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xp, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xm, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xp, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xm, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(1),
                );
            }

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xm, 0),
                    )
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xp, 0),
                    )
                    .e_branch(-1),
                );
            }

            // Real positive axis

            let m = 2.0;
            let p1 = get_branch_point(m - 1.0, consts, 0.0);

            let mut p_points = vec![];
            // let mut x_cuts = vec![];

            let p0 = 1.0 / 8.0;
            let p2 = 7.0 / 8.0;

            let p_int = PInterpolator::xp(p2, consts)
                .goto_xp(p2, m + 1.0)
                .goto_xp(p0, m + 1.0)
                .goto_xp(p0, m);

            // x_cuts.push(p_int.x_path.clone());

            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(0.01, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p1, m);
            let branch_point = *p_int2.p_path.last().unwrap();
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            let p_int = PInterpolator::xp(p0, consts).goto_xp(p0, m);
            // x_cuts.push(p_int.x_path.clone());

            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p1, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(0.01, m);
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            // cuts.push(Cut::new(
            //     Component::Xp,
            //     x_cuts,
            //      vec![],
            //      CutType::DebugPath,
            //
            // ));

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xp, 0),
                    )
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xm, 0),
                    )
                    .e_branch(1),
                );
            }

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xm, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xp, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xm, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xp, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(-1),
                );
            }
        } else if p_range == -1 {
            // Real negative axis

            let mut p_points = vec![];
            // let mut x_cuts = vec![];

            let p0 = p_range as f64 + 1.0 / 8.0;
            let p2 = p_range as f64 + 1.0 - 1.0 / 8.0;

            let m = 2.0;

            let p1 =
                get_branch_point(m - 1.0 + p_start * consts.k() as f64, consts, -1.0) + p_start;

            let p_int = PInterpolator::xp(p2, consts).goto_xp(p2, m);

            // x_cuts.push(p_int.x_path.clone());
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p_start + 0.99, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p1, m);
            let branch_point = *p_int2.p_path.last().unwrap();
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            let p_int = PInterpolator::xp(p0, consts)
                .goto_xp(p0, m + 1.0)
                .goto_xp(p2, m + 1.0)
                .goto_xp(p2, m);

            // x_cuts.push(p_int.x_path.clone());
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p1, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p_start + 0.99, m);
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            // cuts.push(Cut::new(
            //     Component::Xp,
            //     x_cuts,
            //      vec![],
            //      CutType::DebugPath,
            //
            // ));

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xp, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xm, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xp, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xm, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(1),
                );
            }

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xm, 0),
                    )
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xp, 0),
                    )
                    .e_branch(-1),
                );
            }

            // Real positive axis

            let mut p_points = vec![];
            // let mut x_cuts = vec![];

            let p0 = p_range as f64 + 1.0 / 8.0;
            let p2 = p_range as f64 + 1.0 - 1.0 / 8.0;

            let m = 2.0 + p_start * consts.k() as f64;

            let p1 = get_branch_point(m - 1.0, consts, 0.0) + p_start;

            let p_int = PInterpolator::xp(p2, consts)
                .goto_xp(p2, m - 1.0)
                .goto_xp(p0, m - 1.0)
                .goto_xp(p0, m);

            // x_cuts.push(p_int.x_path.clone());
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p_start + 0.01, m);
            p_points.extend(p_int2.p_path.into_iter().rev());

            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p1, m);
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            let p_int = PInterpolator::xp(p2, consts)
                .goto_xp(p2, m + 1.0)
                .goto_xp(p0, m + 1.0)
                .goto_xp(p0, m);

            // x_cuts.push(p_int.x_path.clone());
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p1, m);
            let branch_point = *p_int2.p_path.last().unwrap();
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p_start + 0.01, m);
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            // cuts.push(Cut::new(
            //     Component::Xp,
            //     x_cuts,
            //      vec![],
            //      CutType::DebugPath,
            //
            // ));

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xp, 0),
                    )
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xm, 0),
                    )
                    .e_branch(1),
                );
            }

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xm, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xp, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xm, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xp, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(-1),
                );
            }
        } else if p_range > 0 {
            // Real negative axis

            let mut p_points = vec![];
            // let mut x_cuts = vec![];

            let p0 = p_range as f64 + 1.0 / 8.0;
            let p2 = p_range as f64 + 1.0 - 1.0 / 8.0;

            let m = (p_start + 1.0) * consts.k() as f64 + 2.0;

            let p1 = get_branch_point(m - 1.0 - consts.k() as f64, consts, -1.0) + p_start;

            // let x = get_branch_point_x(m - 1.0 - consts.k() as f64, consts, -1.0);

            let p_int = PInterpolator::xp(p2, consts).goto_xp(p2, m);

            // x_cuts.push(p_int.x_path.clone());
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p_start + 0.99, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p1, m);
            let branch_point = *p_int2.p_path.last().unwrap();
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            let p_int = PInterpolator::xp(p2, consts)
                .goto_xp(p2, m - 1.5)
                .goto_xp(p0, m - 1.5)
                .goto_xp(p0, m + 1.5)
                .goto_xp(p2, m + 1.5)
                .goto_xp(p2, m);

            // x_cuts.push(p_int.x_path.clone());
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p1, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p_start + 0.99, m);
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            // cuts.push(Cut::new(
            //     Component::Xp,
            //     x_cuts,
            //      vec![x],
            //      CutType::DebugPath,
            //     log_branch(p_range)
            // );

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xp, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xm, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xp, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xm, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(1),
                );
            }

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xm, 0),
                    )
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xp, 0),
                    )
                    .e_branch(-1),
                );
            }

            // Real positive axis

            let mut p_points = vec![];
            // let mut x_cuts = vec![];

            let p0 = p_range as f64 + 1.0 / 8.0;
            let p2 = p_range as f64 + 1.0 - 1.0 / 8.0;

            let m = p_start * consts.k() as f64 + 2.0;

            let p1 = get_branch_point(m - 1.0, consts, 0.0) + p_start;

            // let x = get_branch_point_x(m - 1.0, consts, 0.0);

            let p_int = PInterpolator::xp(p0, consts).goto_xp(p0, m);

            // x_cuts.push(p_int.x_path.clone());
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p_start + 0.01, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p1, m);
            let branch_point = *p_int2.p_path.last().unwrap();
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            let p_int = PInterpolator::xp(p2, consts)
                .goto_xp(p2, m + 1.5)
                .goto_xp(p0, m + 1.5)
                .goto_xp(p0, m);

            // x_cuts.push(p_int.x_path.clone());
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p1, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p_start + 0.01, m);
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            // cuts.push(Cut::new(
            //     Component::Xp,
            //     x_cuts,
            //      vec![x],
            //      CutType::DebugPath,
            //     log_branch(p_range)
            // );

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xp, 0),
                    )
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xm, 0),
                    )
                    .e_branch(1),
                );
            }

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xm, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xp, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xm, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xp, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(-1),
                );
            }
        } else if p_range < -1 {
            // Real negative axis

            let mut p_points = vec![];
            // let mut x_cuts = vec![];

            let p0 = p_range as f64 + 1.0 / 8.0;
            let p2 = p_range as f64 + 1.0 - 1.0 / 8.0;

            let m = (p_start + 1.0) * consts.k() as f64 + 2.0;

            let p1 = get_branch_point(m - 1.0 - consts.k() as f64, consts, -1.0) + p_start;

            // let x = get_branch_point_x(m - 1.0 - consts.k() as f64, consts, -1.0);

            let p_int = PInterpolator::xp(p2, consts).goto_xp(p2, m);

            // x_cuts.push(p_int.x_path.clone());
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p_start + 0.99, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p1, m);
            let branch_point = *p_int2.p_path.last().unwrap();
            p_points.extend(p_int2.p_path);
            // // x_cuts.push(p_int2.x_path);

            let p_int = PInterpolator::xp(p2, consts)
                .goto_xp(p2, m - 1.5)
                .goto_xp(p0, m - 1.5)
                .goto_xp(p0, m + 1.5)
                .goto_xp(p2, m + 1.5)
                .goto_xp(p2, m);

            // x_cuts.push(p_int.x_path.clone());
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p1, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p_start + 0.99, m);
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            // cuts.push(Cut::new(
            //     Component::Xp,
            //     x_cuts,
            //      vec![x],
            //      CutType::DebugPath,
            //     log_branch(p_range)
            // );

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xp, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xm, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xp, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xm, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(1),
                );
            }

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xm, 0),
                    )
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xp, 0),
                    )
                    .e_branch(-1),
                );
            }

            // Positive real axis

            let mut p_points = vec![];
            // let mut x_cuts = vec![];

            let p0 = p_range as f64 + 1.0 / 8.0;
            let p2 = p_range as f64 + 1.0 - 1.0 / 8.0;

            let m = (p_start) * consts.k() as f64 + 2.0;

            let p1 = get_branch_point(m - 1.0, consts, 0.0) + p_start;

            let x = get_branch_point_x(m - 1.0, consts, 0.0);

            let p_int = PInterpolator::xp(p0, consts).goto_xp(p0, m);

            // x_cuts.push(p_int.x_path.clone());
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p_start + 0.01, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p1, m);
            let branch_point = *p_int2.p_path.last().unwrap();
            p_points.extend(p_int2.p_path);
            // x_cuts.push(p_int2.x_path);

            let p_int = PInterpolator::xp(p0, consts)
                .goto_xp(p0, m + 1.5)
                .goto_xp(p2, m + 1.5)
                .goto_xp(p2, m - 1.5)
                .goto_xp(p0, m - 1.5)
                .goto_xp(p0, m);

            // x_cuts.push(p_int.x_path.clone());
            let p_int = p_int.clear_path();

            let p_int2 = p_int.clone().goto_xp(p1, m);
            p_points.extend(p_int2.p_path.into_iter().rev());
            // // x_cuts.push(p_int2.x_path);

            let p_int2 = p_int.clone().goto_xp(p_start + 0.01, m);
            p_points.extend(p_int2.p_path);
            // // x_cuts.push(p_int2.x_path);

            // cuts.push(
            //     Cut::new(Component::Xp, x_cuts, vec![x], CutType::DebugPath).log_branch(p_range).e_branch(1),
            // );

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xp, 0),
                    )
                    .e_branch(1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xm, 0),
                    )
                    .e_branch(1),
                );
            }

            {
                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xm, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xp, 1),
                    )
                    .log_branch_le(p_range)
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.iter().map(|p| p.conj()).collect()],
                        vec![branch_point.conj()],
                        CutType::LogX(Component::Xm, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(-1),
                );

                cuts.push(
                    Cut::new(
                        Component::P,
                        vec![p_points.clone()],
                        vec![branch_point],
                        CutType::LogX(Component::Xp, -1),
                    )
                    .log_branch_gt(p_range)
                    .e_branch(-1),
                );
            }
        }

        cuts
    }

    fn e_cuts(p_range: i32, consts: CouplingConstants) -> Vec<Cut> {
        let p_start = p_range as f64;

        let p0 = nr::find_root(
            |p| en2(p, 1.0, consts),
            |p| den2_dp(p, 1.0, consts),
            C::new(p_start, 2.5),
            1.0e-3,
            50,
        );

        let Some(p0) = p0 else { return vec![] };

        let mut cut = vec![];

        cut.push((0.0, p0));
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

            cut.push((im, p));
            p_prev = p;

            if p.im.abs() > 2.0 {
                break;
            }
        }

        let mut cuts = vec![];

        let sheet_data = SheetData {
            log_branch: p_range,
            log_branch_sum: -p_range,
            e_branch: 1,
            u_branch: (1, 1),
        };

        let cut = cut
            .iter()
            .map(|&(im, p)| {
                (
                    im,
                    p,
                    xp(p + 0.00001, 1.0, consts),
                    xp(p - 0.00001, 1.0, consts),
                    xm(p + 0.00001, 1.0, consts),
                    xm(p - 0.00001, 1.0, consts),
                    u(p + 0.00001, consts, &sheet_data),
                    u(p - 0.00001, consts, &sheet_data),
                    u(p.conj() + 0.00001, consts, &sheet_data),
                    u(p.conj() - 0.00001, consts, &sheet_data),
                )
            })
            .collect::<Vec<_>>();

        // p

        let paths = vec![cut
            .iter()
            .map(|(_, p, _, _, _, _, _, _, _, _)| *p)
            .collect()];
        let branch_points = vec![p0];

        cuts.push(Cut::new(Component::P, paths, branch_points, CutType::E));

        let paths = vec![cut
            .iter()
            .map(|(_, p, _, _, _, _, _, _, _, _)| p.conj())
            .collect()];
        let branch_points = vec![p0.conj()];

        cuts.push(Cut::new(Component::P, paths, branch_points, CutType::E));

        // xp

        let mut paths = vec![
            cut.iter()
                .map(|(_, _, xp, _, _, _, _, _, _, _)| *xp)
                .collect::<Vec<_>>(),
            cut.iter()
                .map(|(_, _, _, xp, _, _, _, _, _, _)| *xp)
                .collect::<Vec<_>>(),
        ];
        paths.push(vec![paths[0][0], paths[1][0]]);
        let branch_points = vec![(paths[0][0] + paths[1][0]) / 2.0];
        cuts.push(if branch_points[0].im > 0.0 {
            Cut::new(Component::Xp, paths, branch_points, CutType::E)
                .log_branch(p_range)
                .im_xm_positive()
        } else {
            Cut::new(Component::Xp, paths, branch_points, CutType::E)
                .log_branch(p_range)
                .im_xm_negative()
        });

        let mut paths = vec![
            cut.iter()
                .map(|(_, _, _, _, xm, _, _, _, _, _)| xm.conj())
                .collect::<Vec<_>>(),
            cut.iter()
                .map(|(_, _, _, _, _, xm, _, _, _, _)| xm.conj())
                .collect::<Vec<_>>(),
        ];
        paths.push(vec![paths[0][0], paths[1][0]]);
        let branch_points = vec![(paths[0][0] + paths[1][0]) / 2.0];
        cuts.push(if branch_points[0].im > 0.0 {
            Cut::new(Component::Xp, paths, branch_points, CutType::E)
                .log_branch(p_range)
                .im_xm_positive()
        } else {
            Cut::new(Component::Xp, paths, branch_points, CutType::E)
                .log_branch(p_range)
                .im_xm_negative()
        });

        // xm

        let mut paths = vec![
            cut.iter()
                .map(|(_, _, xp, _, _, _, _, _, _, _)| xp.conj())
                .collect::<Vec<_>>(),
            cut.iter()
                .map(|(_, _, _, xp, _, _, _, _, _, _)| xp.conj())
                .collect::<Vec<_>>(),
        ];
        paths.push(vec![paths[0][0], paths[1][0]]);
        let branch_points = vec![(paths[0][0] + paths[1][0]) / 2.0];

        let xm_lower = branch_points[0].im < 0.0;

        cuts.push(if branch_points[0].im < 0.0 {
            Cut::new(Component::Xm, paths, branch_points, CutType::E)
                .log_branch(p_range)
                .im_xp_negative()
        } else {
            Cut::new(Component::Xm, paths, branch_points, CutType::E)
                .log_branch(p_range)
                .im_xp_positive()
        });

        let mut paths = vec![
            cut.iter()
                .map(|(_, _, _, _, xm, _, _, _, _, _)| *xm)
                .collect::<Vec<_>>(),
            cut.iter()
                .map(|(_, _, _, _, _, xm, _, _, _, _)| *xm)
                .collect::<Vec<_>>(),
        ];
        paths.push(vec![paths[0][0], paths[1][0]]);
        let branch_points = vec![(paths[0][0] + paths[1][0]) / 2.0];

        let xm_lower = branch_points[0].im < 0.0;

        cuts.push(if branch_points[0].im < 0.0 {
            Cut::new(Component::Xm, paths, branch_points, CutType::E)
                .log_branch(p_range)
                .im_xp_negative()
        } else {
            Cut::new(Component::Xm, paths, branch_points, CutType::E)
                .log_branch(p_range)
                .im_xp_positive()
        });

        let mut paths = vec![
            cut.iter()
                .map(|(_, _, _, _, _, _, _, _, u, _)| *u)
                .collect::<Vec<_>>(),
            cut.iter()
                .map(|(_, _, _, _, _, _, _, _, _, u)| *u)
                .collect::<Vec<_>>(),
        ];
        paths.push(vec![paths[0][0], paths[1][0]]);
        let branch_points = vec![(paths[0][0] + paths[1][0]) / 2.0];

        if p_range >= 0 {
            cuts.push(
                Cut::new(Component::U, paths, branch_points, CutType::E)
                    .log_branch(p_range)
                    .im_xm_positive()
                    .im_xp_positive(),
            );
        } else {
            cuts.push(
                Cut::new(Component::U, paths, branch_points, CutType::E)
                    .log_branch(p_range)
                    .im_xp_negative()
                    .im_xm_negative(),
            );
        }

        let mut paths = vec![
            cut.iter()
                .map(|(_, _, _, _, _, _, u, _, _, _)| *u)
                .collect::<Vec<_>>(),
            cut.iter()
                .map(|(_, _, _, _, _, _, _, u, _, _)| *u)
                .collect::<Vec<_>>(),
        ];
        paths.push(vec![paths[0][0], paths[1][0]]);
        let branch_points = vec![(paths[0][0] + paths[1][0]) / 2.0];

        if p_range < 0 {
            cuts.push(
                Cut::new(Component::U, paths, branch_points, CutType::E)
                    .log_branch(p_range)
                    .im_xm_positive()
                    .im_xp_positive(),
            );
        } else {
            cuts.push(
                Cut::new(Component::U, paths, branch_points, CutType::E)
                    .log_branch(p_range)
                    .im_xp_negative()
                    .im_xm_negative(),
            );
        }

        cuts
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
        let log_branch = p.re.floor() as i32;
        let log_branch_sum = -log_branch;

        let sheet_data = SheetData {
            log_branch,
            log_branch_sum,
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
                CutType::LogX(Component::Xp, branch) => {
                    new_sheet_data.log_branch += branch;
                    new_sheet_data.log_branch_sum += branch;
                }
                CutType::LogX(Component::Xm, branch) => {
                    new_sheet_data.log_branch += branch;
                    new_sheet_data.log_branch_sum -= branch;
                }
                CutType::E => {
                    new_sheet_data.e_branch = -new_sheet_data.e_branch;
                }
                CutType::U(Component::Xp) => {
                    new_sheet_data.u_branch =
                        (-new_sheet_data.u_branch.0, new_sheet_data.u_branch.1);
                }
                CutType::U(Component::Xm) => {
                    new_sheet_data.u_branch =
                        (new_sheet_data.u_branch.0, -new_sheet_data.u_branch.1);
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
