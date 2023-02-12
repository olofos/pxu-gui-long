use num::complex::{Complex, ComplexFloat};

type C = Complex<f64>;

pub trait Func {
    fn f(&self, z: C) -> C;
    fn df(&self, z: C) -> C;
}

pub trait OneParameterFunction {
    fn evaluate(&self, t: f64) -> C;
}

pub fn find_root(
    f: impl Fn(C) -> C,
    df: impl Fn(C) -> C,
    guess: C,
    precision_goal: f64,
    max_iterations: usize,
) -> Option<C> {
    let mut result = guess;
    for _ in 0..max_iterations {
        result = result - f(result) / df(result);
        if f(result).abs() < precision_goal {
            return Some(result);
        }
    }
    None
}

const PRECISION_GOAL: f64 = 1.0e-6;
const MAX_ITERATIONS: usize = 50;

pub fn shoot(
    func: &impl Func,
    one_param_func: &impl OneParameterFunction,
    t0: f64,
    t1: f64,
    z0: C,
    max_step: f64,
) -> Vec<(f64, C)> {
    assert!(
        (one_param_func.evaluate(t0) - func.f(z0)).abs() < PRECISION_GOAL,
        "{}",
        (one_param_func.evaluate(t0) - func.f(z0)).abs()
    );

    let mut result = vec![(t0, z0)];
    let dt = t1 - t0;
    let mut s = 0.0;

    let max_step = max_step / dt.abs();

    while s < 1.0 {
        let last = *result.last().unwrap();

        let mut step = max_step.min(1.0 - s);
        s += step;

        let mut i = 0;

        loop {
            let t = t0 + s * dt;
            let next = find_root(
                |z| func.f(z) - one_param_func.evaluate(t),
                |z| func.df(z),
                last.1,
                PRECISION_GOAL,
                MAX_ITERATIONS,
            );
            if let Some(next) = next {
                result.push((t, next));
                break;
            } else {
                step /= 2.0;
                s -= step;
            }
            i += 1;
            if i > 6 {
                log::info!("> ({s})");
                return result;
            }
        }
    }
    result
}

pub fn shoot_two_sided(
    func: &impl Func,
    one_param_func: &impl OneParameterFunction,
    t0: f64,
    t1: f64,
    t2: f64,
    z1: C,
    max_step: f64,
) -> Vec<(f64, C)> {
    let mut result = shoot(func, one_param_func, t1, t0, z1, max_step);
    let result2 = shoot(func, one_param_func, t1, t2, z1, max_step);
    result.reverse();
    result.pop();
    result.extend(result2);
    result
}
