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

pub fn en(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    let sin = (PI * p).sin();
    let m_eff = m + consts.k() as f64 * p;

    (m_eff * m_eff + 4.0 * consts.h * consts.h * sin * sin).sqrt()
}

pub fn den_dp(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    let sin = (PI * p).sin();
    let cos = (PI * p).cos();
    let m_eff = m + consts.k() as f64 * p;

    TAU * (consts.kslash() * m_eff + 2.0 * consts.h * consts.h * sin * cos) / en(p, m, consts)
}

pub fn den_dm(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    let m_eff = m + consts.k() as f64 * p;
    m_eff / en(p, m, consts)
}

pub fn en2(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    let sin = (PI * p).sin();
    let m_eff = m + consts.k() as f64 * p;

    m_eff * m_eff + 4.0 * consts.h * consts.h * sin * sin
}

pub fn den2_dp(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    let sin = (PI * p).sin();
    let cos = (PI * p).cos();
    let m_eff = m + consts.k() as f64 * p;

    TAU * (2.0 * consts.kslash() * m_eff + 4.0 * consts.h * consts.h * sin * cos)
}

const SIGN: f64 = 1.0;

fn x(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    let sin = (PI * p).sin();
    let m_eff = m + consts.k() as f64 * p;

    let numerator = m_eff + SIGN * en(p, m, consts);
    let denominator = 2.0 * consts.h * sin;

    numerator / denominator
}

fn dx_dp(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    let sin = (PI * p).sin();
    let cos = (PI * p).cos();

    let term1 = -x(p, m, consts) * (cos / sin) / 2.0;
    let term2 = consts.kslash() / (2.0 * consts.h * sin);
    let term3 = (consts.kslash() * (m + consts.k() as f64 * p)
        + 2.0 * consts.h * consts.h * sin * cos)
        / (en(p, m, consts) * 2.0 * consts.h * sin);

    TAU * (term1 + term2 + SIGN * term3)
}

fn dx_dm(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    let sin = (PI * p).sin();
    (1.0 + den_dm(p, m, consts)) * (2.0 * consts.h * sin)
}

pub fn xp(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    x(p, m, consts) * (C::i() * PI * p).exp()
}

pub fn dxp_dp(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    let exp = (C::i() * PI * p).exp();
    dx_dp(p, m, consts) * exp + (C::i() * PI) * x(p, m, consts) * exp
}

pub fn dxp_dm(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    let exp = (C::i() * PI * p).exp();
    dx_dm(p, m, consts) * exp
}

pub fn dxm_dm(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    let exp = (-C::i() * PI * p).exp();
    dx_dm(p, m, consts) * exp
}

pub fn xm(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    x(p, m, consts) * (-C::i() * PI * p).exp()
}

pub fn dxm_dp(p: impl Into<C>, m: f64, consts: CouplingConstants) -> C {
    let p = p.into();
    let exp = (-C::i() * PI * p).exp();
    dx_dp(p, m, consts) * exp - (C::i() * PI) * x(p, m, consts) * exp
}

pub fn u(p: impl Into<C>, consts: CouplingConstants, log_branch: i32) -> C {
    let p = p.into();
    let xp = xp(p, 1.0, consts);
    let xm = xm(p, 1.0, consts);

    let up = xp + 1.0 / xp - 2.0 * consts.kslash() / consts.h * xp.ln();

    let um = xm + 1.0 / xm - 2.0 * consts.kslash() / consts.h * xm.ln();

    let branch_shift = log_branch as f64 * consts.k() as f64 * C::i() / consts.h;

    let t = 0.0;
    let s = 0.5;
    let u0p = up - C::i() / consts.h - (1.0 + t) * branch_shift;
    let u0m = um + C::i() / consts.h + (1.0 - t) * branch_shift;

    if (u0p - u0m).norm_sqr() > 0.01 {
        log::info!("{:.2}", u0p - u0m);
    }

    s * u0p + (1.0 - s) * u0m

    // ((xp + 1.0 / xp - 2.0 * consts.kslash() / consts.h * xp.ln() - C::i() * (1.0) / consts.h)
    //     + (xm + 1.0 / xm - 2.0 * consts.kslash() / consts.h * xm.ln() + C::i() * (1.0) / consts.h))
    //     / 2.0

    // - 0.0 * log_branch as f64 * consts.k() as f64 * C::i() / consts.h)
}

pub fn du_dp(p: impl Into<C>, consts: CouplingConstants) -> C {
    let p = p.into();
    let cot = 1.0 / (PI * p).tan();
    let sin = (PI * p).sin();

    let term1 = den_dp(p, 1.0, consts) * cot;
    let term2 = -TAU * en(p, 1.0, consts) / (2.0 * sin * sin);
    let term3 = -2.0 * consts.kslash() * dx_dp(p, 1.0, consts) / x(p, 1.0, consts);

    (term1 + term2 + term3) * consts.h
}
