use std::f64::consts::TAU;

use indicatif::{ProgressBar, ProgressStyle};
use num::complex::Complex64;
use pxu::kinematics::CouplingConstants;

fn main() -> std::io::Result<()> {
    let consts = CouplingConstants::new(2.0, 5);

    let mut contours = pxu::Contours::new();

    let spinner_style = ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
    )
    .unwrap()
    .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");

    println!("[1/?] Generating pxu.contours");
    let pb = ProgressBar::new(1);
    pb.set_style(spinner_style);
    loop {
        pb.set_length(contours.progress().1 as u64);
        pb.set_position(contours.progress().0 as u64);
        if contours.update(0, consts) {
            pb.finish_and_clear();
            break;
        }
    }

    if false {
        let mut state = pxu::State::new(1, consts);
        state.update(
            0,
            pxu::Component::P,
            Complex64::new(0.03, 0.03),
            &contours,
            consts,
        );
        state.update(
            0,
            pxu::Component::P,
            Complex64::new(-0.03, 0.03),
            &contours,
            consts,
        );
        state.update(
            0,
            pxu::Component::P,
            Complex64::new(-0.06, 0.0),
            &contours,
            consts,
        );

        let center = Complex64::new(-0.3, 0.5);
        let radius = 1.2;
        let steps = 128;

        let mut path = vec![];

        state.update(0, pxu::Component::Xp, center - radius, &contours, consts);
        for i in 0..=(4 * steps) {
            let theta = 6.0 * (i as f64 / steps as f64 - 0.5);
            let xp = center + Complex64::from_polar(radius, theta);
            path.push(xp);
        }

        let saved_path = pxu::path::SavedPath {
            base_path: pxu::path::BasePath {
                path,
                start: state,
                component: pxu::Component::Xp,
                excitation: 0,
            },
            consts,
        };

        let s = serde_json::to_string(&saved_path)?;
        println!("{s}");
    }

    if false {
        let center = Complex64::new(0.0, 0.0);
        let radius = 0.05;
        let steps = 128;

        let mut state = pxu::State::new(1, consts);
        state.update(0, pxu::Component::P, center + radius, &contours, consts);

        let mut path = vec![];

        for i in 0..=(steps) {
            let theta = TAU * (i as f64 / steps as f64);
            let z = center + Complex64::from_polar(radius, theta);
            path.push(z);
        }

        let saved_path = pxu::path::SavedPath {
            base_path: pxu::path::BasePath {
                path,
                start: state,
                component: pxu::Component::P,
                excitation: 0,
            },
            consts,
        };

        let s = serde_json::to_string(&saved_path)?;
        println!("{s}");
    }

    if true {
        let center = Complex64::new(0.0, 0.0);
        let radius = 0.10;
        let steps = 128;

        let mut state = pxu::State::new(1, consts);
        state.update(0, pxu::Component::P, center + radius, &contours, consts);

        let mut path = vec![];

        for i in 0..=(steps) {
            let theta = TAU * (i as f64 / steps as f64);
            let z = center + Complex64::from_polar(radius, theta);
            path.push(z);
        }

        let saved_path = pxu::path::SavedPath {
            base_path: pxu::path::BasePath {
                path,
                start: state,
                component: pxu::Component::P,
                excitation: 0,
            },
            consts,
        };

        let s = serde_json::to_string(&saved_path)?;
        println!("{s}");
    }

    Ok(())
}
