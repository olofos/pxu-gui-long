use num::complex::{Complex, ComplexFloat};

type C = Complex<f64>;

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
