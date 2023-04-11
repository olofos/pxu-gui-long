use crate::contours::{Component, UCutType};
use crate::cut::{Cut, CutType};
use crate::kinematics::{
    du_crossed_dp, du_dp, dxm_crossed_dp, dxm_dp, dxp_crossed_dp, dxp_dp, u, u_crossed, xm,
    xm_crossed, xp, xp_crossed, CouplingConstants, SheetData, UBranch,
};
use crate::nr;
use num::complex::Complex64;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
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
