use crate::kinematics::{CouplingConstants, UBranch};
pub use crate::point::Point;

use crate::pxu::{Component, UCutType};
use itertools::Itertools;

use num::complex::Complex64;

#[derive(Debug, Clone)]
pub struct Cut {
    pub component: Component,
    pub paths: Vec<Vec<Complex64>>,
    pub branch_point: Option<Complex64>,
    pub typ: CutType,
    pub p_range: i32,
    pub periodic: bool,
    pub(crate) visibility: Vec<CutVisibilityCondition>,
}

impl Cut {
    pub fn new(
        component: Component,
        paths: Vec<Vec<Complex64>>,
        branch_point: Option<Complex64>,
        typ: CutType,
        p_range: i32,
        periodic: bool,
        visibility: Vec<CutVisibilityCondition>,
    ) -> Self {
        Self {
            component,
            paths,
            branch_point,
            typ,
            p_range,
            periodic,
            visibility,
        }
    }

    pub fn conj(&self) -> Self {
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
            periodic: self.periodic,
            p_range: self.p_range,
        }
    }

    pub fn shift_conj(&self, dz: Complex64) -> Self {
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
            periodic: self.periodic,
            p_range: self.p_range,
        }
    }

    pub fn shift(mut self, dz: Complex64) -> Self {
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

    pub fn intersection(
        &self,
        p1: Complex64,
        p2: Complex64,
        consts: CouplingConstants,
    ) -> Option<(usize, usize, Complex64)> {
        if self.periodic {
            let period = 2.0 * Complex64::i() * consts.k() as f64 / consts.h;
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

#[derive(Debug, Clone, PartialEq)]
pub enum CutVisibilityCondition {
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
