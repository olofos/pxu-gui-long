use num::complex::Complex64;
use std::f64::consts::{PI, TAU};

#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct CouplingConstants {
    pub h: f64,
    k: f64,
}

impl CouplingConstants {
    pub fn new(h: f64, k: i32) -> Self {
        Self { h, k: k as f64 }
    }

    pub fn k(&self) -> i32 {
        self.k.round() as i32
    }

    pub fn kslash(&self) -> f64 {
        self.k / TAU
    }

    pub fn s(&self) -> f64 {
        ((self.kslash() * self.kslash() + self.h * self.h).sqrt() + self.kslash()) / self.h
    }

    pub fn get_set_k(&mut self, k: Option<f64>) -> f64 {
        if let Some(k) = k {
            self.k = k;
        }
        self.k() as f64
    }
}

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum UBranch {
    Outside,
    Between,
    Inside,
}

impl std::fmt::Display for UBranch {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Self::Outside => "outside",
            Self::Between => "between",
            Self::Inside => "inside",
        };
        write!(f, "{}", s)
    }
}

#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub struct SheetData {
    pub log_branch_p: i32,
    pub log_branch_m: i32,
    pub log_branch_x: i32,
    pub e_branch: i32,
    pub u_branch: (UBranch, UBranch),
    pub im_x_sign: (i8, i8),
}

pub fn en(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let sin = (PI * p).sin();
    let m_eff = m + consts.k() as f64 * p;

    (m_eff * m_eff + 4.0 * consts.h * consts.h * sin * sin).sqrt()
}

pub fn den_dp(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let sin = (PI * p).sin();
    let cos = (PI * p).cos();
    let m_eff = m + consts.k() as f64 * p;

    TAU * (consts.kslash() * m_eff + 2.0 * consts.h * consts.h * sin * cos) / en(p, m, consts)
}

pub fn den_dm(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let m_eff = m + consts.k() as f64 * p;
    m_eff / en(p, m, consts)
}

pub fn en2(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let sin = (PI * p).sin();
    let m_eff = m + consts.k() as f64 * p;

    m_eff * m_eff + 4.0 * consts.h * consts.h * sin * sin
}

pub fn den2_dp(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let sin = (PI * p).sin();
    let cos = (PI * p).cos();
    let m_eff = m + consts.k() as f64 * p;

    TAU * (2.0 * consts.kslash() * m_eff + 4.0 * consts.h * consts.h * sin * cos)
}

const SIGN: f64 = 1.0;

fn xpm_common(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let sin = (PI * p).sin();
    let m_eff = m + consts.k() as f64 * p;

    let numerator = m_eff + SIGN * en(p, m, consts);
    let denominator = 2.0 * consts.h * sin;

    numerator / denominator
}

fn dxpm_common_dp(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let sin = (PI * p).sin();
    let cos = (PI * p).cos();

    let term1 = -xpm_common(p, m, consts) * (cos / sin) / 2.0;
    let term2 = consts.kslash() / (2.0 * consts.h * sin);
    let term3 = (consts.kslash() * (m + consts.k() as f64 * p)
        + 2.0 * consts.h * consts.h * sin * cos)
        / (en(p, m, consts) * 2.0 * consts.h * sin);

    TAU * (term1 + term2 + SIGN * term3)
}

pub fn xp(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    xpm_common(p, m, consts) * (Complex64::i() * PI * p).exp()
}

pub fn dxp_dp(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let exp = (Complex64::i() * PI * p).exp();
    dxpm_common_dp(p, m, consts) * exp + (Complex64::i() * PI) * xpm_common(p, m, consts) * exp
}

pub fn xm(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    xpm_common(p, m, consts) * (-Complex64::i() * PI * p).exp()
}

pub fn dxm_dp(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let exp = (-Complex64::i() * PI * p).exp();
    dxpm_common_dp(p, m, consts) * exp - (Complex64::i() * PI) * xpm_common(p, m, consts) * exp
}

pub fn u(p: impl Into<Complex64>, consts: CouplingConstants, sheet_data: &SheetData) -> Complex64 {
    let p = p.into();
    let xp = xp_on_sheet(p, 1.0, consts, sheet_data);

    let up = xp + 1.0 / xp - 2.0 * consts.kslash() / consts.h * xp.ln();
    let branch_shift =
        2.0 * (sheet_data.log_branch_p * consts.k()) as f64 * Complex64::i() / consts.h;

    up - Complex64::i() / consts.h - branch_shift
}

pub fn du_dp(
    p: impl Into<Complex64>,
    consts: CouplingConstants,
    sheet_data: &SheetData,
) -> Complex64 {
    let p = p.into();
    let cot = 1.0 / (PI * p).tan();
    let sin = (PI * p).sin();

    let term1;
    let term2;
    let term3;

    if sheet_data.e_branch > 0 {
        term1 = den_dp(p, 1.0, consts) * cot;
        term2 = -TAU * en(p, 1.0, consts) / (2.0 * sin * sin);
        term3 =
            -2.0 * consts.kslash() * dxpm_common_dp(p, 1.0, consts) / xpm_common(p, 1.0, consts);
    } else {
        term1 = -den_dp(p, 1.0, consts) * cot;
        term2 = TAU * en(p, 1.0, consts) / (2.0 * sin * sin);
        term3 = -2.0 * consts.kslash() * dxpm_common_crossed_dp(p, 1.0, consts)
            / xpm_common_crossed(p, 1.0, consts);
    }

    (term1 + term2 + term3) * consts.h
}

fn xpm_common_crossed(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let sin = (PI * p).sin();
    let m_eff = m + consts.k() as f64 * p;

    let numerator = m_eff - SIGN * en(p, m, consts);
    let denominator = 2.0 * consts.h * sin;

    numerator / denominator
}

fn dxpm_common_crossed_dp(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let sin = (PI * p).sin();
    let cos = (PI * p).cos();

    let term1 = -xpm_common_crossed(p, m, consts) * (cos / sin) / 2.0;
    let term2 = consts.kslash() / (2.0 * consts.h * sin);
    let term3 = (consts.kslash() * (m + consts.k() as f64 * p)
        + 2.0 * consts.h * consts.h * sin * cos)
        / (en(p, m, consts) * 2.0 * consts.h * sin);

    TAU * (term1 + term2 - SIGN * term3)
}

pub fn xp_crossed(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    xpm_common_crossed(p, m, consts) * (Complex64::i() * PI * p).exp()
}

pub fn dxp_crossed_dp(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let exp = (Complex64::i() * PI * p).exp();
    dxpm_common_crossed_dp(p, m, consts) * exp
        + (Complex64::i() * PI) * xpm_common_crossed(p, m, consts) * exp
}

pub fn xm_crossed(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    xpm_common_crossed(p, m, consts) * (-Complex64::i() * PI * p).exp()
}

pub fn dxm_crossed_dp(p: impl Into<Complex64>, m: f64, consts: CouplingConstants) -> Complex64 {
    let p = p.into();
    let exp = (-Complex64::i() * PI * p).exp();
    dxpm_common_crossed_dp(p, m, consts) * exp
        - (Complex64::i() * PI) * xpm_common_crossed(p, m, consts) * exp
}

pub fn xp_on_sheet(
    p: impl Into<Complex64>,
    m: f64,
    consts: CouplingConstants,
    sheet_data: &SheetData,
) -> Complex64 {
    if sheet_data.e_branch > 0 {
        xp(p, m, consts)
    } else {
        xp_crossed(p, m, consts)
    }
}

pub fn xm_on_sheet(
    p: impl Into<Complex64>,
    m: f64,
    consts: CouplingConstants,
    sheet_data: &SheetData,
) -> Complex64 {
    if sheet_data.e_branch > 0 {
        xm(p, m, consts)
    } else {
        xm_crossed(p, m, consts)
    }
}

pub fn dxp_dp_on_sheet(
    p: impl Into<Complex64>,
    m: f64,
    consts: CouplingConstants,
    sheet_data: &SheetData,
) -> Complex64 {
    if sheet_data.e_branch > 0 {
        dxp_dp(p, m, consts)
    } else {
        dxp_crossed_dp(p, m, consts)
    }
}

pub fn dxm_dp_on_sheet(
    p: impl Into<Complex64>,
    m: f64,
    consts: CouplingConstants,
    sheet_data: &SheetData,
) -> Complex64 {
    if sheet_data.e_branch > 0 {
        dxm_dp(p, m, consts)
    } else {
        dxm_crossed_dp(p, m, consts)
    }
}

pub fn u_of_x(x: impl Into<Complex64>, consts: CouplingConstants) -> Complex64 {
    let s = consts.s();

    let x: Complex64 = x.into();

    x + 1.0 / x - (s - 1.0 / s) * x.ln()
}

pub fn du_dx(x: impl Into<Complex64>, consts: CouplingConstants) -> Complex64 {
    let s = consts.s();

    let x: Complex64 = x.into();

    (x - s) * (x + 1.0 / s) / (x * x)
}
