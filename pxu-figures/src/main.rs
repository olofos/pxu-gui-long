use ::pxu::pxu::{ContourGenerator, GridLine, GridLineComponent};
use ::pxu::{
    interpolation::{InterpolationPoint, PInterpolatorMut},
    kinematics::CouplingConstants,
    pxu,
};
use itertools::Itertools;
use num::complex::Complex64;
use std::fs::File;
use std::io::{prelude::*, BufWriter, Result};
use std::ops::Range;
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};
use std::sync::Arc;
use std::thread;

mod cache;

fn error(message: &str) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::Other, message)
}

#[derive(Debug, Clone, PartialEq)]
struct Bounds {
    x_range: Range<f64>,
    y_range: Range<f64>,
}

impl Bounds {
    fn new(x_range: Range<f64>, y_range: Range<f64>) -> Self {
        Self { x_range, y_range }
    }

    fn width(&self) -> f64 {
        self.x_range.end - self.x_range.start
    }

    fn height(&self) -> f64 {
        self.y_range.end - self.y_range.start
    }

    fn inside(&self, z: &Complex64) -> bool {
        self.x_range.contains(&z.re) && self.y_range.contains(&z.im)
    }

    fn crosses(&self, z1: &Complex64, z2: &Complex64) -> bool {
        (z1.re < self.x_range.start) && (z2.re > self.x_range.end)
            || (z2.re < self.x_range.start) && (z1.re > self.x_range.end)
            || (z1.im < self.y_range.start) && (z2.im > self.y_range.end)
            || (z2.im < self.y_range.start) && (z1.im > self.y_range.end)
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Size {
    width: f64,
    height: f64,
}

#[derive(Debug)]
struct Figure {
    name: String,
    bounds: Bounds,
    size: Size,
    writer: BufWriter<File>,
}

const LUALATEX: &str = "lualatex.exe";
const FIGURE_PATH: &str = "./figures";
const TEX_EXT: &str = "tex";
const PDF_EXT: &str = "pdf";

impl Figure {
    const FILE_START: &str = r#"
\nonstopmode
\documentclass[10pt,a4paper]{article}
\usepackage{pgfplots}
\pgfplotsset{compat=1.18}
\usepgfplotslibrary{fillbetween}
\usepackage[active,tightpage]{preview}
\PreviewEnvironment{tikzpicture}
\setlength\PreviewBorder{0pt}
\pdfvariable suppressoptionalinfo \numexpr 1023 \relax
\begin{document}
\pagestyle{empty}
\begin{tikzpicture}
"#;

    const FILE_END: &str = r#"
\end{tikzpicture}
\end{document}
"#;

    fn new(name: &str, bounds: Bounds, size: Size) -> std::io::Result<Self> {
        let mut path = PathBuf::from(FIGURE_PATH).join(name);
        path.set_extension(TEX_EXT);

        log::info!("[{name}]: Creating file {}", path.to_string_lossy());

        let file = File::create(&path)?;
        let mut writer = BufWriter::new(file);

        writer.write_all(Self::FILE_START.as_bytes())?;

        let x_min = bounds.x_range.start;
        let x_max = bounds.x_range.end;

        let y_min = bounds.y_range.start;
        let y_max = bounds.y_range.end;

        let width = size.width;
        let height = size.height;

        writeln!(writer, "\\begin{{axis}}[hide axis,scale only axis,ticks=none,xmin={x_min},xmax={x_max},ymin={y_min},ymax={y_max},clip,width={width}cm,height={height}cm]")?;
        // writeln!(writer, "\\begin{{axis}}[hide axis,scale only axis,ticks=none,width={width}cm,height={height}cm]")?;

        Ok(Self {
            name: name.to_owned(),
            writer,
            bounds,
            size,
        })
    }

    fn format_coordinate(p: Complex64) -> String {
        format!("({:.3},{:.3})", p.re, p.im)
    }

    fn format_contour(&self, contour: Vec<Complex64>) -> Vec<String> {
        let mut coordinates = contour
            .into_iter()
            .map(Self::format_coordinate)
            .collect::<Vec<_>>();
        coordinates.dedup();
        coordinates
    }

    fn crop(&self, contour: &Vec<Complex64>) -> Vec<Complex64> {
        // return contour.clone();
        if contour.len() < 2 {
            return vec![];
        }

        let mut coordinates: Vec<Complex64> = vec![];

        let include = |z1, z2| {
            self.bounds.inside(z1) || self.bounds.inside(z2) || self.bounds.crosses(z1, z2)
        };

        if let [z1, z2] = &contour[0..=1] {
            if include(z1, z2) {
                coordinates.push(*z1);
            }
        }

        for (z1, z2, z3) in contour.iter().tuple_windows::<(_, _, _)>() {
            if include(z1, z2) || include(z2, z3) {
                coordinates.push(*z2);
            }
        }

        if let [z1, z2] = &contour[(contour.len() - 2)..=(contour.len() - 1)] {
            if include(z1, z2) {
                coordinates.push(*z2);
            }
        }

        coordinates
    }

    fn add_plot(&mut self, options: &[&str], contour: &Vec<Complex64>) -> Result<()> {
        self.add_plot_all(options, self.crop(contour))
    }

    fn add_plot_all(&mut self, options: &[&str], contour: Vec<Complex64>) -> Result<()> {
        let coordinates = self.format_contour(contour);

        if !coordinates.is_empty() {
            writeln!(
                self.writer,
                "\\addplot [{}] coordinates {{ {} }};",
                options.join(","),
                coordinates.join(" ")
            )?;
        }
        Ok(())
    }

    fn add_grid_line(&mut self, grid_line: &GridLine, options: &[&str]) -> Result<()> {
        self.add_plot(&[&["very thin", "gray"], options].concat(), &grid_line.path)?;

        Ok(())
    }

    fn add_cut(&mut self, cut: &pxu::Cut, options: &[&str]) -> Result<()> {
        let (color, style) = match cut.typ {
            pxu::CutType::E => ("black", ""),
            pxu::CutType::Log(pxu::Component::Xp) => ("red", ""),
            pxu::CutType::Log(pxu::Component::Xm) => ("green", ""),
            pxu::CutType::ULongPositive(pxu::Component::Xp) => ("red", ""),
            pxu::CutType::ULongNegative(pxu::Component::Xp) => ("red", "dashed"),
            pxu::CutType::ULongPositive(pxu::Component::Xm) => ("green", ""),
            pxu::CutType::ULongNegative(pxu::Component::Xm) => ("green", "dashed"),
            _ => {
                return Ok(());
            }
        };

        for path in cut.paths.iter() {
            self.add_plot(&[&[color, style], options].concat(), path)?
        }

        if let Some(branch_point) = cut.branch_point {
            self.add_plot_all(
                &[&[color, "only marks", "mark size=0.02cm"], options].concat(),
                vec![branch_point],
            )?;
        }

        Ok(())
    }

    fn add_axis(&mut self) -> Result<()> {
        let options = ["very thin", "black"];
        self.add_plot(
            &options,
            &vec![
                Complex64::new(self.bounds.x_range.start - 1.0, 0.0),
                Complex64::new(self.bounds.x_range.end + 1.0, 0.0),
            ],
        )?;
        self.add_plot(
            &options,
            &vec![
                Complex64::new(0.0, self.bounds.y_range.start - 1.0),
                Complex64::new(0.0, self.bounds.y_range.end + 1.0),
            ],
        )
    }

    fn finnish(mut self, cache: Arc<cache::Cache>) -> std::io::Result<FigureCompiler> {
        writeln!(self.writer, "\\end{{axis}}\n")?;
        self.writer.write_all(Self::FILE_END.as_bytes())?;
        self.writer.flush()?;

        FigureCompiler::new(self, cache)
    }

    fn transform_vec(&self, v: Complex64) -> Complex64 {
        Complex64::new(
            v.re * self.size.width / self.bounds.width(),
            v.im * self.size.height / self.bounds.height(),
        )
    }
}

struct FigureCompiler {
    name: String,
    child: Child,
}

impl FigureCompiler {
    fn new(figure: Figure, cache: Arc<cache::Cache>) -> Result<Self> {
        let name = figure.name;
        if cache.check(&name)? {
            log::info!("[{name}] Matches cached entry");
            let child = Command::new("/bin/true").spawn()?;
            Ok(Self { name, child })
        } else {
            let mut path = PathBuf::from(FIGURE_PATH).join(name.clone());
            path.set_extension(TEX_EXT);

            let mut cmd = Command::new(LUALATEX);
            cmd.arg(format!("--output-directory={}", FIGURE_PATH))
                .args(["--interaction=nonstopmode", "--output-format=pdf"])
                .arg(path.as_os_str())
                .stderr(Stdio::null())
                .stdout(Stdio::null());

            log::info!("[{name}]: Running Lualatex");
            let child = cmd.spawn()?;

            Ok(Self { name, child })
        }
    }

    fn wait(mut self) -> Result<String> {
        if self.child.wait()?.success() {
            log::info!("[{}]: Lualatex done.", self.name);
        } else {
            log::error!("[{}]: Lualatex failed.", self.name);
        }
        Ok(self.name)
    }
}

trait Node {
    fn write_m_node(
        &mut self,
        figure: &mut Figure,
        anchor: &str,
        rot_sign: i32,
        consts: CouplingConstants,
    ) -> Result<()>;
}

impl Node for PInterpolatorMut {
    fn write_m_node(
        &mut self,
        figure: &mut Figure,
        anchor: &str,
        rot_sign: i32,
        consts: CouplingConstants,
    ) -> Result<()> {
        let p = match self.pt() {
            InterpolationPoint::Xp(p, _) | InterpolationPoint::Xm(p, _) => p,
            _ => {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::Other,
                    format!("Expected xp or xm, found {:?}", self.pt()),
                ));
            }
        };

        let mut p_int2 = self.clone();
        p_int2.goto_p(p - 0.001);
        let p1 = p_int2.p();
        p_int2.goto_p(p + 0.001);
        let p2 = p_int2.p();
        let dp = figure.transform_vec(p2 - p1);
        let rotation = dp.im.atan2(dp.re) * 180.0 / std::f64::consts::PI
            + if rot_sign >= 0 { 0.0 } else { 180.0 };

        let (color, m) = match self.pt().normalized(consts) {
            InterpolationPoint::Xp(_, m) => ("black", m),
            InterpolationPoint::Xm(_, m) => ("blue", m),
            _ => unreachable!(),
        };

        writeln!(figure.writer,"\\node[scale=0.3,anchor={anchor},inner sep=0.4pt,rotate={rotation:.1},{color}] at ({:.3}, {:.3}) {{$\\scriptscriptstyle {}$}};",
                 self.p().re,
                 self.p().im,
                 m)
    }
}

fn fig_xpl_preimage(
    contour_generator: Arc<ContourGenerator>,
    cache: Arc<cache::Cache>,
) -> Result<FigureCompiler> {
    let mut figure = Figure::new(
        "xpL-preimage",
        Bounds::new(-2.6..2.6, -0.7..0.7),
        Size {
            width: 10.0,
            height: 4.0,
        },
    )?;

    let Some(consts) = contour_generator.consts else {
        return Err(error("No consts set"));
    };

    let pt = &pxu::PxuPoint::new(0.5, consts);

    for contour in contour_generator.get_grid(pxu::Component::P).iter()
    // .take(100)
    {
        figure.add_grid_line(contour, &[])?;
    }

    for cut in contour_generator
        .get_visible_cuts(pt, pxu::Component::P, true)
        .filter(|cut| matches!(cut.typ, pxu::CutType::E))
    {
        figure.add_cut(cut, &[])?;
    }

    let consts = pt.consts;
    let k = consts.k() as f64;

    for p_range in -3..=2 {
        let p_start = p_range as f64;

        let bp1 = pxu::compute_branch_point(
            p_range,
            pxu::BranchPointType::XpPositiveAxisImXmNegative,
            consts,
        )
        .unwrap();

        let bp2 = pxu::compute_branch_point(
            p_range,
            pxu::BranchPointType::XpNegativeAxisFromAboveWithImXmNegative,
            consts,
        )
        .unwrap();

        let p0 = p_start + (bp1.p + bp2.p) / 2.0;
        let mut p_int = PInterpolatorMut::xp(p0, consts);

        for m in 1..=15 {
            p_int.goto_m(m as f64);
            p_int.write_m_node(&mut figure, "south", 1, consts)?;
        }

        let mut p_int = PInterpolatorMut::xp(p0, consts);

        for m in (-15..=0).rev() {
            p_int.goto_m(m as f64);
            p_int.write_m_node(&mut figure, "south", 1, consts)?;
        }

        // let p0 = p_start + 0.5;
        let p1 = p_start + bp1.p - 0.003;

        let m = p_range * consts.k() + 2;

        let mut p_int = PInterpolatorMut::xp(p0, consts);
        p_int.goto_m(m as f64 + 1.0 * (p_start + 0.5).signum());
        p_int.goto_p(p1);

        p_int.goto_m(m as f64);
        if p_range >= 0 {
            p_int.write_m_node(&mut figure, "north", 1, consts)?;
        } else {
            p_int.write_m_node(&mut figure, "north", -1, consts)?;
        }

        let p2 = p_start + bp1.p + 0.01;
        if p_range > 0 {
            p_int.goto_m(m as f64 - 1.0).goto_p(p2);

            for dm in (1..=4).rev() {
                p_int.goto_m((m - dm) as f64);
                p_int.write_m_node(&mut figure, "south", -1, consts)?;
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
            let mut p_int = PInterpolatorMut::xp(p0, consts);
            p_int
                .goto_m(-p_start * k + 1.0)
                .goto_p(p_start + 0.1)
                .goto_conj()
                .goto_p(p2);
            for dm in 0..=consts.k() {
                p_int.goto_m(-p_start * k - dm as f64);
                p_int.write_m_node(&mut figure, "south", -1, consts)?;
            }
        }
    }

    figure.finnish(cache)
}

fn fig_xpl_cover(
    contour_generator: Arc<ContourGenerator>,
    cache: Arc<cache::Cache>,
) -> Result<FigureCompiler> {
    let mut figure = Figure::new(
        "xpL-cover",
        Bounds::new(-6.0..6.0, -0.2..4.0),
        Size {
            width: 6.0,
            height: 4.0,
        },
    )?;

    figure.add_axis()?;
    for contour in contour_generator
        .get_grid(pxu::Component::Xp)
        .iter()
        .filter(
            |line| matches!(line.component, GridLineComponent::Xp(m) if (-8.0..=5.0).contains(&m)),
        )
    {
        figure.add_grid_line(contour, &["thin", "black"])?;
    }
    figure.finnish(cache)
}

fn fig_p_plane_long_cuts_regions(
    contour_generator: Arc<ContourGenerator>,
    cache: Arc<cache::Cache>,
) -> Result<FigureCompiler> {
    let mut figure = Figure::new(
        "p_plane_long_cuts_regions",
        Bounds::new(-2.6..2.6, -0.7..0.7),
        Size {
            width: 10.0,
            height: 4.0,
        },
    )?;

    let Some(consts) = contour_generator.consts else {
        return Err(error("No consts set"));
    };

    let pt = &pxu::PxuPoint::new(0.5, consts);

    let color_physical = "yellow!10";
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

    for cut in contour_generator.get_visible_cuts(pt, pxu::Component::P, true) {
        let color_mirror = match cut.typ {
            pxu::CutType::ULongPositive(pxu::Component::Xp)
            | pxu::CutType::ULongNegative(pxu::Component::Xp) => color_mirror_p,
            pxu::CutType::ULongPositive(pxu::Component::Xm)
            | pxu::CutType::ULongNegative(pxu::Component::Xm) => color_mirror_m,
            _ => {
                continue;
            }
        };

        let mut cropped_path = figure.crop(&cut.paths[0]);
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

    for contour in contour_generator.get_grid(pxu::Component::P).iter() {
        figure.add_grid_line(contour, &[])?;
    }

    for cut in contour_generator
        .get_visible_cuts(pt, pxu::Component::P, true)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::ULongNegative(_) | pxu::CutType::ULongPositive(_)
            )
        })
    {
        figure.add_cut(cut, &[])?;
    }

    figure.finnish(cache)
}

fn main() -> std::io::Result<()> {
    tracing_subscriber::fmt::fmt()
        .with_writer(std::io::stderr)
        .init();

    let cache = cache::Cache::load(FIGURE_PATH)?;

    let consts = CouplingConstants::new(2.0, 5);
    let contour_generator = pxu::ContourGenerator::generate_all(consts);

    // let figures = vec![
    //     fig_xpl_preimage(contour_generator_ref.clone())?,
    //     fig_xpl_cover(contour_generator_ref.clone())?,
    //     fig_p_plane_long_cuts_regions(contour_generator_ref.clone())?,
    // ];

    // for c in figures {
    //     c.wait()?;
    // }

    let fig_functions = [
        fig_xpl_preimage,
        fig_xpl_cover,
        fig_p_plane_long_cuts_regions,
    ];

    let contour_generator_ref = Arc::new(contour_generator);
    let cache_ref = Arc::new(cache);

    let handles = fig_functions
        .into_iter()
        .map(|f| {
            let contour_generator_ref = contour_generator_ref.clone();
            let cache_ref = cache_ref.clone();
            thread::spawn(move || {
                let figure = f(contour_generator_ref, cache_ref)?;
                let result = figure.wait()?;
                Result::Ok(result)
            })
        })
        .collect::<Vec<_>>();

    let mut new_cache = cache::Cache::new(FIGURE_PATH);

    for handle in handles {
        let name = handle
            .join()
            .map_err(|err| error(format!("Join error: {err:?}").as_str()))??;
        new_cache.update(&name)?;
    }

    log::info!("Saving cache");
    new_cache.save();

    Ok(())
}
