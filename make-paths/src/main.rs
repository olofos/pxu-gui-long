use std::f64::consts::{PI, TAU};

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

    if false {
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

    if true {
        fn goto(
            state: &mut pxu::State,
            component: pxu::Component,
            new_value: impl Into<Complex64>,
            contours: &pxu::Contours,
            consts: CouplingConstants,
            steps: usize,
        ) {
            let z0 = state.points[0].get(component);
            let z1 = new_value.into();

            for i in 0..=steps {
                let z = z0 + (i as f64 / steps as f64) * (z1 - z0);
                state.update(0, component, z, contours, consts);
            }

            if (state.points[0].get(component) - z1).norm() > 1.0e-6 {
                println!(
                    "Could not goto ({})",
                    (state.points[0].get(component) - z1).norm()
                );
            }
        }

        let mut state = pxu::State::new(1, consts);

        let x0 = 2.8;
        let y0 = -1.5;
        let r = 0.25;

        goto(&mut state, pxu::Component::U, 0.0, &contours, consts, 10);
        goto(
            &mut state,
            pxu::Component::U,
            Complex64::new(0.0, y0),
            &contours,
            consts,
            3,
        );

        goto(
            &mut state,
            pxu::Component::U,
            -x0 + r,
            &contours,
            consts,
            16,
        );

        let mut path = vec![Complex64::new(-x0 + r, y0)];

        let steps = (0..=4)
            .map(|n| PI / 2.0 * n as f64 / 4.0)
            .collect::<Vec<_>>();
        for y in 0..=2 {
            let c = Complex64::new(x0 - r, y0 + r + 5.0 * y as f64);
            for theta in steps.iter() {
                path.push(c + Complex64::from_polar(r, *theta - PI / 2.0));
            }

            let c = Complex64::new(x0 - r, y0 - r + 5.0 * (y as f64 + 0.5));
            for theta in steps.iter() {
                path.push(c + Complex64::from_polar(r, *theta));
            }

            let c = Complex64::new(-x0 + r, y0 + r + 5.0 * (y as f64 + 0.5));
            for theta in steps.iter() {
                path.push(c + Complex64::from_polar(r, -PI / 2.0 - *theta));
            }

            let c = Complex64::new(-x0 + r, y0 - r + 5.0 * (y as f64 + 1.0));
            for theta in steps.iter() {
                path.push(c + Complex64::from_polar(r, PI - *theta));
            }
        }

        let saved_path = pxu::path::SavedPath {
            base_path: pxu::path::BasePath {
                path,
                start: state,
                component: pxu::Component::U,
                excitation: 0,
            },
            consts,
        };

        let s = serde_json::to_string(&saved_path)?;
        println!("{s}");
    }

    Ok(())
}
