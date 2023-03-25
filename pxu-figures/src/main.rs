use ::pxu::{interpolation::PInterpolatorMut, kinematics::CouplingConstants, pxu};
use itertools::Itertools;
use num::complex::Complex64;
use std::fs::File;
use std::io::{prelude::*, BufWriter};
use std::process::{Command, Stdio};


fn crop(
    contour: &Vec<Complex64>,
    x_min: f64,
    x_max: f64,
    y_min: f64,
    y_max: f64,
) -> Vec<Complex64> {
    let x_range = x_min..=x_max;
    let y_range = y_min..=y_max;
    let inside = |z: &Complex64| x_range.contains(&z.re) && y_range.contains(&z.im);

    let crosses = |z1: &Complex64, z2: &Complex64| {
        (z1.re < x_min) && (z2.re > x_max)
            || (z2.re < x_min) && (z1.re > x_max)
            || (z1.im < y_min) && (z2.im > y_max)
            || (z2.im < y_min) && (z1.im > y_max)
    };

    if contour.len() < 2 {
        return vec![];
    }

    let mut coordinates: Vec<Complex64> = vec![];

    if inside(&contour[1]) || crosses(&contour[0], &contour[1]) {
        coordinates.push(contour[0]);
    }

    for (z1, z2, z3) in contour.iter().tuple_windows::<(_, _, _)>() {
        if inside(z1) || inside(z2) || inside(z3) || crosses(z1, z2) || crosses(z2, z3) {
            coordinates.push(*z2);
        }
    }

    if inside(&contour[contour.len() - 2])
        || crosses(&contour[contour.len() - 2], &contour[contour.len() - 1])
    {
        coordinates.push(contour[contour.len() - 1]);
    }

    coordinates
}

fn format_coordinate(p: Complex64) -> String {
    format!("({:.3},{:.3})", p.re, p.im)
}

fn format_contour(contour: Vec<Complex64>) -> Vec<String> {
    let mut coordinates = contour
        .into_iter()
        .map(format_coordinate)
        .collect::<Vec<_>>();
    coordinates.dedup();
    coordinates
}

fn main() -> std::io::Result<()> {
    tracing_subscriber::fmt::fmt()
        .with_writer(std::io::stderr)
        .init();

    let mut contour_generator = pxu::ContourGenerator::new();

    let consts = CouplingConstants::new(2.0, 5);
    let pt = pxu::PxuPoint::new(1.5, consts);

    while !contour_generator.update(&pt) {}

    {
        let the_file = File::create("test.tex")?;
        let mut file = BufWriter::new(the_file);

        let begin = r#"
\nonstopmode
\documentclass[10pt,a4paper]{article}
\usepackage{pgfplots}
\pgfplotsset{compat=1.18}
\usepackage[active,tightpage]{preview}
\PreviewEnvironment{tikzpicture}
\setlength\PreviewBorder{0pt}
\pdfvariable suppressoptionalinfo \numexpr 1023 \relax
\begin{document}
\pagestyle{empty}
\begin{tikzpicture}
"#;

        let end = r#"
\end{tikzpicture}
\end{document}
"#;

        file.write_all(begin.as_bytes())?;

        let xmin = -2.6;
        let xmax = 2.6;
        let ymin = -0.7;
        let ymax = 0.7;

        writeln!(file, "\\begin{{axis}}[hide axis,width=10cm,height=6cm,ticks=none,xmin={xmin},xmax={xmax},ymin={ymin},ymax={ymax},clip]")?;
        // writeln!(file, "\\begin{{axis}}[hide axis,width=10cm,height=6cm,ticks=none,xmin={xmin},xmax={xmax},ymin={ymin},ymax={ymax},clip mode=individual]")?;

        for contour in contour_generator.get_grid(pxu::Component::P)
        // .iter()
        // .take(200)
        {
            let coordinates = format_contour(crop(contour, xmin, xmax, ymin, ymax));

            if !coordinates.is_empty() {
                writeln!(
                    file,
                    "\\addplot [very thin,gray] coordinates {{ {} }};",
                    coordinates.join(" ")
                )?;
            }
        }

        for cut in contour_generator
            .get_visible_cuts(&pt, pxu::Component::P, true)
            .filter(|cut| matches!(cut.typ, pxu::CutType::E))
        {
            let color = match cut.typ {
                pxu::CutType::E => "black",
                pxu::CutType::Log(pxu::Component::Xp) => "red",
                pxu::CutType::Log(pxu::Component::Xm) => "green",
                _ => {
                    continue;
                }
            };

            let coordinates = format_contour(crop(&cut.paths[0], xmin, xmax, ymin, ymax));

            if !coordinates.is_empty() {
                writeln!(
                    file,
                    "\\addplot [{color}] coordinates {{ {} }};",
                    coordinates.join(" ")
                )?;

                if let Some(branch_point) = cut.branch_point {
                    writeln!(
                        file,
                        "\\addplot [{color}, only marks, mark size=0.02cm] coordinates {{ {} }};",
                        format_coordinate(branch_point)
                    )?;
                }
            }
        }

        let mut write_m = |p_int: &PInterpolatorMut, p_range: i32, p0: f64, m, anchor| {
            let mut p_int2 = p_int.clone();
            p_int2.goto_p(p0 - 0.001);
            let p1 = p_int2.p();
            p_int2.goto_p(p0 + 0.001);
            let p2 = p_int2.p();
            let dp = p2 - p1;

            let dx = dp.re * 10.0 / (xmax - xmin);
            let dy = dp.im * 6.0 / (ymax - ymin);

            let angle = dy.atan2(dx) * 180.0 / std::f64::consts::PI;

            p_int2.goto_p(p0);
            p_int2.goto_m(m as f64);

            let dir = |angle: f64| {
                // if (-45.0..=135.0).contains(&angle) {
                //     angle
                // } else if angle > 0.0 {
                //     angle - 180.0
                // } else {
                //     angle + 180.0
                // }
                angle
            };

            {
                let rotation = dir(angle);
                writeln!(
                file,
                "\\node[scale=0.35,anchor={anchor},inner sep=1.0pt,rotate={rotation:.1}] at ({:.3}, {:.3}) {{$\\scriptscriptstyle {}$}};",
                    p_int2.p().re,
                    p_int2.p().im,
                    m + p_range * consts.k(),
                ).expect("error");
            }

            if p_int2.p().im.abs() > 0.001 {
                let rotation = dir(-angle);
                writeln!(
                file,
                "\\node[scale=0.35,anchor={anchor},inner sep=1.0pt,rotate={rotation:.1}] at ({:.3}, {:.3}) {{$\\scriptscriptstyle {}$}};",
                    p_int2.p().re,
                    -p_int2.p().im,
                    p_range * consts.k() + 2 - m,
                ).expect("error");
            }
        };

        // let p_range = 1;
        for p_range in -3..=2 {
            let p_start = p_range as f64;

            let bp = pxu::compute_branch_point(
                p_range,
                pxu::BranchPointType::XpPositiveAxisImXmNegative,
                consts,
            )
            .unwrap();

            log::info!(
                "bp({p_range}) = ({},{}) -> ({},{})",
                bp.p,
                bp.m,
                bp.p + p_start,
                bp.m - p_start * consts.k() as f64
            );

            let p0 = p_start + 0.45;
            let mut p_int = PInterpolatorMut::xp(p0, consts);

            for m in 1..=15 {
                p_int.goto_m(m as f64);
                write_m(&p_int, p_range, p0, m, "south");
            }

            let p0 = p_start + 0.5;
            let p1 = p_start + bp.p - 0.003;

            let m = p_range * consts.k() + 2;

            let mut p_int = PInterpolatorMut::xp(p0, consts);
            p_int.goto_m(m as f64 + 1.0 * (p_start + 0.5).signum());
            p_int.goto_p(p1);

            p_int.goto_m(m as f64);
            write_m(&p_int, p_range, p1, m, "north");

            let p2 = p1 + 0.1;
            p_int.goto_m(m as f64 - 1.0).goto_p(p2);

            if p_range > 0 {
                for dm in (1..=4).rev() {
                    p_int.goto_m((m - dm) as f64);
                    write_m(&p_int, p_range, p2, m - dm, "north");
                }
            }

            if p_range > 0 {
                let p0 = p_start + 0.4;

                let mut p_int = PInterpolatorMut::xp(p0, consts);
                p_int.goto_m(m as f64 - 1.0);
                p_int.goto_p(p1);
                p_int.goto_m(m as f64 + 1.0);

                let p2 = p1 + 0.1;
                p_int.goto_p(p2);

                for m in ((m + 1)..(m + 4)).rev() {
                    p_int.goto_m(m as f64);
                    write_m(&p_int, p_range, p2, m, "north");
                }
            }
        }

        writeln!(file, "\\end{{axis}}")?;

        file.write_all(end.as_bytes())?;
    }

    {
        let mut cmd = Command::new("lualatex.exe");
        cmd.args(["--interaction=nonstopmode", "--output-format=pdf", "test"])
            .stderr(Stdio::null())
            .stdout(Stdio::null());

        log::info!("Running Lualatex",);
        let mut child = cmd.spawn()?;

        if child.wait()?.success() {
            log::info!("done.");
        } else {
            log::error!("Lualatex failed.");
        }
    }

    // contour_generator.get_visible_cuts(pt, pxu::Component::P, true)

    Ok(())
}
