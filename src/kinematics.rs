use num::complex::Complex;
use std::f64::consts::{PI, TAU};

type C = Complex<f64>;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CouplingConstants {
    pub h: f64,
    kslash_: f64,
}

struct _Context {
    pub h: f64,
    pub kslash: f64,

    pub p_branch: i32,
    pub log_branch: i32,
}

impl CouplingConstants {
    pub fn new(h: f64, k: u32) -> Self {
        Self {
            h,
            kslash_: k as f64 / TAU,
        }
    }

    pub fn k(&self) -> u32 {
        (TAU * self.kslash() + 0.1).floor() as u32
    }

    pub fn kslash(&self) -> f64 {
        self.kslash_
    }

    pub fn s(&self) -> f64 {
        ((self.kslash() * self.kslash() + self.h * self.h).sqrt() + self.kslash()) / self.h
    }

    pub fn get_set_k(&mut self, k: Option<f64>) -> f64 {
        if let Some(k) = k {
            self.kslash_ = k / std::f64::consts::TAU;
        }
        self.k() as f64
    }

    pub fn get_set_s(&mut self, s: Option<f64>) -> f64 {
        if let Some(s) = s {
            self.h = 2.0 * self.kslash() * s / (s * s - 1.0);
        }
        self.s()
    }
}

const DP_DP: f64 = TAU;
const PP: f64 = 1.0;

pub fn en(p: C, m: f64, consts: CouplingConstants) -> C {
    let p = p / PP;
    let sin = (PI * p).sin();
    let m_eff = m + consts.k() as f64 * p;

    (m_eff * m_eff + 4.0 * consts.h * consts.h * sin * sin).sqrt()
}

pub fn den_dp(p: C, m: f64, consts: CouplingConstants) -> C {
    let p = p / PP;

    let sin = (PI * p).sin();
    let cos = (PI * p).cos();
    let m_eff = m + consts.k() as f64 * p;

    DP_DP * (consts.kslash() * m_eff + 2.0 * consts.h * consts.h * sin * cos)
        / en(PP * p, m, consts)
}

pub fn en2(p: C, m: f64, consts: CouplingConstants) -> C {
    let p = p / PP;

    let sin = (PI * p).sin();
    let m_eff = m + consts.k() as f64 * p;

    m_eff * m_eff + 4.0 * consts.h * consts.h * sin * sin
}

pub fn den2_dp(p: C, m: f64, consts: CouplingConstants) -> C {
    let p = p / PP;

    let sin = (PI * p).sin();
    let cos = (PI * p).cos();
    let m_eff = m + consts.k() as f64 * p;

    DP_DP * (2.0 * consts.kslash() * m_eff + 4.0 * consts.h * consts.h * sin * cos)
}

fn x(p: C, m: f64, consts: CouplingConstants) -> C {
    let p = p / PP;

    let sin = (PI * p).sin();
    let m_eff = m + consts.k() as f64 * p;

    let numerator = m_eff + en(PP * p, m, consts);
    let denominator = 2.0 * consts.h * sin;

    numerator / denominator
}

fn dx_dp(p: C, m: f64, consts: CouplingConstants) -> C {
    let p = p / PP;

    let sin = (PI * p).sin();
    let cos = (PI * p).cos();

    let term1 = -x(PP * p, m, consts) * (cos / sin) / 2.0;
    let term2 = consts.kslash() / (2.0 * consts.h * sin);
    let term3 = (consts.kslash() * (m + consts.k() as f64 * p)
        + 2.0 * consts.h * consts.h * sin * cos)
        / (en(PP * p, m, consts) * 2.0 * consts.h * sin);

    DP_DP * (term1 + term2 + term3)
}

pub fn xp(p: C, m: f64, consts: CouplingConstants) -> C {
    let p = p / PP;
    x(PP * p, m, consts) * (C::i() * PI * p).exp()
}

pub fn dxp_dp(p: C, m: f64, consts: CouplingConstants) -> C {
    let p = p / PP;

    let exp = (C::i() * PI * p).exp();
    DP_DP * (dx_dp(PP * p, m, consts) / DP_DP * exp + (C::i() / 2.0) * x(PP * p, m, consts) * exp)
}

pub fn xm(p: C, m: f64, consts: CouplingConstants) -> C {
    let p = p / PP;

    x(PP * p, m, consts) * (-C::i() * PI * p).exp()
}

pub fn dxm_dp(p: C, m: f64, consts: CouplingConstants) -> C {
    let p = p / PP;

    let exp = (-C::i() * PI * p).exp();
    DP_DP * (dx_dp(PP * p, m, consts) / DP_DP * exp - (C::i() / 2.0) * x(PP * p, m, consts) * exp)
}

pub fn u(p: C, p_range: i32, consts: CouplingConstants) -> C {
    let p = p / PP;

    let xp = xp(PP * p, 1.0, consts);

    xp + 1.0 / xp
        - 2.0 * consts.kslash() / consts.h * xp.ln()
        - C::i() * (1.0 + consts.k() as f64 * p_range as f64) / consts.h
}

pub fn du_dp(p: C, consts: CouplingConstants) -> C {
    let p = p / PP;

    let cot = 1.0 / (PI * p).tan();
    let sin = (PI * p).sin();

    let term1 = den_dp(PP * p, 1.0, consts) / DP_DP * cot;
    let term2 = -en(PP * p, 1.0, consts) / (2.0 * sin * sin);
    let term3 =
        -2.0 * consts.kslash() * dx_dp(PP * p, 1.0, consts) / DP_DP / x(PP * p, 1.0, consts);

    DP_DP * (term1 + term2 + term3) * consts.h
}
