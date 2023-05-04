use crate::cache;
use crate::fig_compiler::FigureCompiler;
use crate::fig_writer::{FigureWriter, Node};
use crate::utils::{error, Settings, Size};

use num::complex::Complex64;
use pxu::GridLineComponent;
use pxu::{interpolation::PInterpolatorMut, kinematics::UBranch, Pxu};
use std::io::Result;
use std::sync::Arc;

fn fig_p_xpl_preimage(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "p-xpL-preimage",
        -2.6..2.6,
        0.0,
        Size {
            width: 15.0,
            height: 6.0,
        },
        pxu::Component::P,
        pxu::UCutType::Long,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::P, pxu::UCutType::Long, 0)
        .filter(|cut| matches!(cut.typ, pxu::CutType::E))
    {
        figure.add_cut(cut, &[], pxu.consts)?;
    }

    let k = pxu.consts.k() as f64;

    for p_range in -3..=2 {
        let p_start = p_range as f64;

        let bp1 = pxu::compute_branch_point(
            p_range,
            pxu::BranchPointType::XpPositiveAxisImXmNegative,
            pxu.consts,
        )
        .unwrap();

        let bp2 = pxu::compute_branch_point(
            p_range,
            pxu::BranchPointType::XpNegativeAxisFromAboveWithImXmNegative,
            pxu.consts,
        )
        .unwrap();

        let p0 = p_start + (bp1.p + bp2.p) / 2.0;
        let mut p_int = PInterpolatorMut::xp(p0, pxu.consts);

        for m in 1..=15 {
            p_int.goto_m(m as f64);
            p_int.write_m_node(&mut figure, "south", 1, pxu.consts)?;
        }

        let mut p_int = PInterpolatorMut::xp(p0, pxu.consts);

        for m in (-15..=0).rev() {
            p_int.goto_m(m as f64);
            p_int.write_m_node(&mut figure, "south", 1, pxu.consts)?;
        }

        // let p0 = p_start + 0.5;
        let p1 = p_start + bp1.p - 0.003;

        let m = p_range * pxu.consts.k() + 2;

        let mut p_int = PInterpolatorMut::xp(p0, pxu.consts);
        p_int.goto_m(m as f64 + 1.0 * (p_start + 0.5).signum());
        p_int.goto_p(p1);

        p_int.goto_m(m as f64);
        if p_range >= 0 {
            p_int.write_m_node(&mut figure, "north", 1, pxu.consts)?;
        } else {
            p_int.write_m_node(&mut figure, "north", -1, pxu.consts)?;
        }

        let p2 = p_start + bp1.p + 0.01;
        if p_range > 0 {
            p_int.goto_m(m as f64 - 1.0).goto_p(p2);

            for dm in (1..=4).rev() {
                p_int.goto_m((m - dm) as f64);
                p_int.write_m_node(&mut figure, "south", -1, pxu.consts)?;
            }
        }

        // if p_range > 0 {
        //     let p0 = p_start + 0.4;

        //     let mut p_int = PInterpolatorMut::xp(p0, consts);
        //     p_int.goto_m(m as f64 - 1.0);
        //     p_int.goto_p(p1);
        //     p_int.goto_m(m as f64 + 1.0);

        //     let p2 = p1;
        //     p_int.goto_p(p2);

        //     for m in ((m + 1)..(m + 4)).rev() {
        //         p_int.goto_m(m as f64);
        //         p_int.write_m_node(&mut figure, "south", consts)?;
        //     }
        // }

        let p2 = p_start + 0.2;
        if p_range != 0 {
            let mut p_int = PInterpolatorMut::xp(p0, pxu.consts);
            p_int
                .goto_m(-p_start * k + 1.0)
                .goto_p(p_start + 0.1)
                .goto_conj()
                .goto_p(p2);
            for dm in 0..=pxu.consts.k() {
                p_int.goto_m(-p_start * k - dm as f64);
                p_int.write_m_node(&mut figure, "south", -1, pxu.consts)?;
            }
        }
    }

    figure.finish(cache, settings)
}

fn fig_xpl_cover(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xpL-cover",
        -5.0..5.0,
        1.9,
        Size {
            width: 6.0,
            height: 3.0,
        },
        pxu::Component::Xp,
        pxu::UCutType::Long,
        settings,
    )?;

    figure.add_axis()?;
    for contour in pxu.contours.get_grid(pxu::Component::Xp).iter().filter(
        |line| matches!(line.component, GridLineComponent::Xp(m) if (-8.0..=6.0).contains(&m)),
    ) {
        figure.add_grid_line(contour, &["thin", "black"])?;
    }
    figure.finish(cache, settings)
}

fn fig_p_plane_long_cuts_regions(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "p-plane-long-cuts-regions",
        -2.6..2.6,
        0.0,
        Size {
            width: 15.0,
            height: 6.0,
        },
        pxu::Component::P,
        pxu::UCutType::Long,
        settings,
    )?;

    let color_physical = "blue!10";
    let color_mirror_p = "red!10";
    let color_mirror_m = "green!10";

    {
        let x1 = figure.bounds.x_range.start - 0.25;
        let x2 = figure.bounds.x_range.end + 0.25;
        let y1 = figure.bounds.y_range.start - 0.25;
        let y2 = figure.bounds.y_range.end + 0.25;

        figure.add_plot_all(
            &[format!("fill={color_physical}").as_str()],
            vec![
                Complex64::new(x1, y1),
                Complex64::new(x1, y2),
                Complex64::new(x2, y2),
                Complex64::new(x2, y1),
            ],
        )?;
    }

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::P, pxu::UCutType::Long, 0)
    {
        let color_mirror = match cut.typ {
            pxu::CutType::ULongPositive(pxu::Component::Xp)
            | pxu::CutType::ULongNegative(pxu::Component::Xp) => color_mirror_p,
            pxu::CutType::ULongPositive(pxu::Component::Xm)
            | pxu::CutType::ULongNegative(pxu::Component::Xm) => color_mirror_m,
            _ => {
                continue;
            }
        };

        let mut cropped_path = figure.crop(&cut.path);
        if cropped_path.len() >= 2 {
            let len = cropped_path.len();
            let start = cropped_path[0];
            let mid = cropped_path[len / 2];
            let end = cropped_path[len - 1];

            cropped_path.push(Complex64::new(mid.re.round(), end.im));
            cropped_path.push(Complex64::new(mid.re.round(), start.im));

            figure.add_plot_all(
                &["draw=none", format!("fill={color_mirror}").as_str()],
                cropped_path,
            )?;
        }
    }

    figure.add_grid_lines(&pxu, &[])?;
    figure.add_cuts(&pxu, &[])?;

    figure.finish(cache, settings)
}

fn fig_p_plane_short_cuts(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "p-plane-short-cuts",
        -2.6..2.6,
        0.0,
        Size {
            width: 25.0,
            height: 10.0,
        },
        pxu::Component::P,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;
    figure.add_cuts(&pxu, &[])?;

    figure.finish(cache, settings)
}

fn fig_xp_cuts_1(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xp-cuts-1",
        -4.0..4.0,
        0.0,
        Size {
            width: 12.0,
            height: 12.0,
        },
        pxu::Component::Xp,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_axis()?;
    for contour in pxu.contours
        .get_grid(pxu::Component::Xp)
        .iter()
        .filter(|line| matches!(line.component, GridLineComponent::Xp(m) | GridLineComponent::Xm(m) if (-10.0..).contains(&m)))
    {
        figure.add_grid_line(contour, &[])?;
    }

    figure.add_cuts(&pxu, &[])?;

    figure.finish(cache, settings)
}

fn fig_u_period_between_between(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "u-period-between-between",
        -6.0..4.0,
        0.25,
        Size {
            width: 5.0,
            height: 12.5,
        },
        pxu::Component::U,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;

    let mut pxu = (*pxu).clone();
    pxu.state.points[0].sheet_data.u_branch = (
        ::pxu::kinematics::UBranch::Between,
        ::pxu::kinematics::UBranch::Between,
    );

    let path = pxu
        .get_path_by_name("U period between/between")
        .ok_or_else(|| error("Path not found"))?;

    figure.add_cuts(&pxu, &[])?;
    figure.add_path(&pxu, path, &[])?;

    figure.finish(cache, settings)
}

fn fig_u_band_between_outside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "u-band-between-outside",
        -6.0..4.0,
        0.25,
        Size {
            width: 5.0,
            height: 12.5,
        },
        pxu::Component::U,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;

    let mut pxu = (*pxu).clone();
    pxu.state.points[0].sheet_data.u_branch = (
        ::pxu::kinematics::UBranch::Between,
        ::pxu::kinematics::UBranch::Outside,
    );

    let path = pxu
        .get_path_by_name("U band between/outside")
        .ok_or_else(|| error("Path not found"))?;

    figure.add_cuts(&pxu, &[])?;
    figure.add_path(&pxu, path, &[])?;

    figure.finish(cache, settings)
}

fn fig_u_band_between_inside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "u-band-between-inside",
        -6.0..4.0,
        0.25,
        Size {
            width: 5.0,
            height: 12.5,
        },
        pxu::Component::U,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;

    let mut pxu = (*pxu).clone();
    pxu.state.points[0].sheet_data.u_branch = (
        ::pxu::kinematics::UBranch::Between,
        ::pxu::kinematics::UBranch::Inside,
    );

    let path = pxu
        .get_path_by_name("U band between/inside")
        .ok_or_else(|| error("Path not found"))?;

    figure.add_cuts(&pxu, &[])?;
    figure.add_path(&pxu, path, &[])?;

    figure.finish(cache, settings)
}

fn fig_p_band_between_outside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "p-band-between-outside",
        -2.6..2.6,
        0.0,
        Size {
            width: 15.0,
            height: 6.0,
        },
        pxu::Component::P,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;

    let path = pxu
        .get_path_by_name("U band between/outside")
        .ok_or_else(|| error("Path not found"))?;

    let mut pxu = (*pxu).clone();
    pxu.state = path.base_path.start.clone();

    figure.add_cuts(&pxu, &[])?;
    figure.add_path(&pxu, path, &[])?;

    figure.finish(cache, settings)
}

fn fig_p_band_between_inside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "p-band-between-inside",
        -2.6..2.6,
        0.0,
        Size {
            width: 15.0,
            height: 6.0,
        },
        pxu::Component::P,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;
    let path = pxu
        .get_path_by_name("U band between/inside")
        .ok_or_else(|| error("Path not found"))?;

    figure.add_cuts(&pxu, &[])?;
    figure.add_path(&pxu, path, &[])?;

    figure.finish(cache, settings)
}

fn fig_xp_band_between_inside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xp-band-between-inside",
        -2.6..2.6,
        0.0,
        Size {
            width: 8.0,
            height: 8.0,
        },
        pxu::Component::Xp,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;
    let path = pxu
        .get_path_by_name("U band between/inside (single)")
        .ok_or_else(|| error("Path not found"))?;

    let mut pxu = (*pxu).clone();
    pxu.state.points[0].sheet_data.u_branch = (UBranch::Between, UBranch::Inside);
    pxu.state.points[0].sheet_data.log_branch_p = 0;
    pxu.state.points[0].sheet_data.log_branch_m = -1;
    pxu.state.points[0].sheet_data.im_x_sign = (1, -1);

    figure.add_cuts(&pxu, &[])?;
    figure.add_path(&pxu, path, &["solid"])?;

    figure.finish(cache, settings)
}

fn fig_xp_band_between_outside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xp-band-between-outside",
        -3.1..2.1,
        0.0,
        Size {
            width: 8.0,
            height: 8.0,
        },
        pxu::Component::Xp,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;

    let path = pxu
        .get_path_by_name("U band between/outside (single)")
        .ok_or_else(|| error("Path not found"))?;

    let mut pxu = (*pxu).clone();
    pxu.state.points[0].sheet_data.u_branch = (UBranch::Between, UBranch::Outside);
    pxu.state.points[0].sheet_data.log_branch_p = 0;
    pxu.state.points[0].sheet_data.log_branch_m = -1;
    pxu.state.points[0].sheet_data.im_x_sign = (1, -1);

    figure.add_cuts(&pxu, &[])?;
    figure.add_path(&pxu, path, &["solid"])?;

    figure.finish(cache, settings)
}

fn fig_xm_band_between_inside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xm-band-between-inside",
        -1.3..1.3,
        0.0,
        Size {
            width: 8.0,
            height: 8.0,
        },
        pxu::Component::Xm,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;

    let path = pxu
        .get_path_by_name("U band between/inside")
        .ok_or_else(|| error("Path not found"))?;

    let mut pxu = (*pxu).clone();
    pxu.state.points[0].sheet_data.u_branch = (UBranch::Between, UBranch::Inside);
    pxu.state.points[0].sheet_data.log_branch_p = 0;
    pxu.state.points[0].sheet_data.log_branch_m = -1;
    pxu.state.points[0].sheet_data.im_x_sign = (1, -1);

    figure.add_cuts(&pxu, &[])?;
    figure.add_path(&pxu, path, &[])?;

    figure.finish(cache, settings)
}

fn fig_xp_period_between_between(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xp-period-between-between",
        -3.1..2.1,
        0.0,
        Size {
            width: 8.0,
            height: 8.0,
        },
        pxu::Component::Xp,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;
    let path = pxu
        .get_path_by_name("U period between/between (single)")
        .ok_or_else(|| error("Path not found"))?;

    let mut pxu = (*pxu).clone();
    pxu.state = path.base_path.start.clone();

    figure.add_cuts(&pxu, &[])?;
    figure.add_path(&pxu, path, &[])?;

    figure.finish(cache, settings)
}

fn fig_xm_period_between_between(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xm-period-between-between",
        -3.1..2.1,
        0.0,
        Size {
            width: 8.0,
            height: 8.0,
        },
        pxu::Component::Xm,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;
    let path = pxu
        .get_path_by_name("U period between/between (single)")
        .ok_or_else(|| error("Path not found"))?;

    let mut pxu = (*pxu).clone();
    pxu.state = path.base_path.start.clone();

    figure.add_cuts(&pxu, &[])?;
    figure.add_path(&pxu, path, &[])?;

    figure.finish(cache, settings)
}

fn fig_p_period_between_between(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "p-period-between-between",
        -0.15..0.15,
        0.0,
        Size {
            width: 8.0,
            height: 8.0,
        },
        pxu::Component::P,
        pxu::UCutType::Short,
        settings,
    )?;

    figure.add_grid_lines(&pxu, &[])?;
    let path = pxu
        .get_path_by_name("U period between/between (single)")
        .ok_or_else(|| error("Path not found"))?;

    let mut pxu = (*pxu).clone();
    pxu.state = path.base_path.start.clone();

    figure.add_cuts(&pxu, &[])?;
    figure.add_path(&pxu, path, &[])?;

    figure.finish(cache, settings)
}

type FigureFunction =
    fn(pxu: Arc<Pxu>, cache: Arc<cache::Cache>, settings: &Settings) -> Result<FigureCompiler>;

pub const ALL_FIGURES: &[FigureFunction] = &[
    fig_p_xpl_preimage,
    fig_xpl_cover,
    fig_p_plane_long_cuts_regions,
    fig_p_plane_short_cuts,
    fig_xp_cuts_1,
    fig_u_period_between_between,
    fig_u_band_between_outside,
    fig_u_band_between_inside,
    fig_p_band_between_outside,
    fig_p_band_between_inside,
    fig_xp_band_between_inside,
    fig_xp_band_between_outside,
    fig_xm_band_between_inside,
    fig_xp_period_between_between,
    fig_xm_period_between_between,
    fig_p_period_between_between,
];
