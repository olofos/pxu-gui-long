use num::complex::Complex;
use std::f64::consts::TAU;

type C = Complex<f64>;

#[derive(Debug, Clone, Copy)]
pub struct CouplingConstants {
    pub h: f64,
    pub kslash: f64,
}

impl CouplingConstants {
    pub fn new(h: f64, k: u32) -> Self {
        Self {
            h,
            kslash: k as f64 / TAU,
        }
    }

    pub fn k(&self) -> u32 {
        (TAU * self.kslash + 0.1).floor() as u32
    }

    pub fn s(&self) -> f64 {
        ((self.kslash * self.kslash + self.h * self.h).sqrt() + self.kslash) / self.h
    }

    pub fn get_set_k(&mut self, k: Option<f64>) -> f64 {
        if let Some(k) = k {
            self.kslash = k / std::f64::consts::TAU;
        }
        self.k() as f64
    }

    pub fn get_set_s(&mut self, s: Option<f64>) -> f64 {
        if let Some(s) = s {
            self.h = 2.0 * self.kslash * s / (s * s - 1.0);
        }
        self.s()
    }
}

fn dispersion_relation(p: C, m: f64, consts: CouplingConstants) -> C {
    let sin = (p / 2.0).sin();
    let m_eff = m + consts.kslash * p;

    (m_eff * m_eff + 4.0 * consts.h * consts.h * sin * sin).sqrt()
}

fn x(p: C, m: f64, consts: CouplingConstants) -> C {
    let sin = (p / 2.0).sin();
    let m_eff = m + consts.kslash * p;

    let numerator = m_eff + dispersion_relation(p, m, consts);
    let denominator = 2.0 * consts.h * sin;

    numerator / denominator
}

fn dx_dp(p: C, m: f64, consts: CouplingConstants) -> C {
    let sin = (p / 2.0).sin();
    let cos = (p / 2.0).cos();

    let term1 = -x(p, m, consts) * (cos / sin) / 2.0;
    let term2 = consts.kslash / (2.0 * consts.h * sin);
    let term3 = (consts.kslash * (m + consts.kslash * p) + 2.0 * consts.h * consts.h * sin * cos)
        / (dispersion_relation(p, m, consts) * 2.0 * consts.h * sin);

    term1 + term2 + term3
}

pub fn xp(p: C, m: f64, consts: CouplingConstants) -> C {
    x(p, m, consts) * (C::i() * p / 2.0).exp()
}

pub fn dxp_dp(p: C, m: f64, consts: CouplingConstants) -> C {
    let i_half = C::new(0.0, 0.5);
    let exp = (i_half * p).exp();
    dx_dp(p, m, consts) * exp + (i_half) * x(p, m, consts) * exp
}

pub fn xm(p: C, m: f64, consts: CouplingConstants) -> C {
    x(p, m, consts) * (-C::i() * p / 2.0).exp()
}

pub fn dxm_dp(p: C, m: f64, consts: CouplingConstants) -> C {
    let i_half = C::new(0.0, 0.5);
    let exp = (-i_half * p).exp();
    dx_dp(p, m, consts) * exp - (i_half) * x(p, m, consts) * exp
}

pub fn u(x: C, consts: CouplingConstants) -> C {
    x + 1.0 / x - 2.0 * consts.kslash / consts.h * x.ln()
}
