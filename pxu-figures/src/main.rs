use itertools::Itertools;
use num::complex::Complex64;
use pxu::{
    interpolation::{InterpolationPoint, PInterpolatorMut},
    kinematics::CouplingConstants,
    State,
};
use pxu::{Contours, GridLine, GridLineComponent};
use std::fs::File;
use std::io::{prelude::*, BufWriter, Result};
use std::ops::Range;
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};
use std::sync::Arc;
use std::thread;

use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

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
struct FigureWriter {
    name: String,
    bounds: Bounds,
    size: Size,
    writer: BufWriter<File>,
    plot_count: u64,
    component: pxu::Component,
}

const LUALATEX: &str = "lualatex";
const FIGURE_PATH: &str = "./figures";
const TEX_EXT: &str = "tex";
const SUMMARY_NAME: &str = "all-figures";

impl FigureWriter {
    const FILE_START_1: &str = r#"
\nonstopmode
\documentclass[10pt,a4paper]{article}
\usepackage{luatextra}
\begin{luacode}
progress_file=io.open(""#;
    const FILE_START_2: &str = r#"","w")
\end{luacode}
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
\directlua{progress_file:write("!")}
\directlua{io.close(progress_file)}
\end{document}
"#;

    fn new(
        name: &str,
        bounds: Bounds,
        size: Size,
        component: pxu::Component,
    ) -> std::io::Result<Self> {
        let mut path = PathBuf::from(FIGURE_PATH).join(name);
        path.set_extension(TEX_EXT);

        log::info!("[{name}]: Creating file {}", path.to_string_lossy());

        let file = File::create(&path)?;
        let mut writer = BufWriter::new(file);

        let mut progress_path = path.clone();
        progress_path.set_extension("prg");
        writer.write_all(Self::FILE_START_1.as_bytes())?;
        write!(writer, "{}", progress_path.to_string_lossy())?;
        writer.write_all(Self::FILE_START_2.as_bytes())?;

        let _ = std::fs::remove_file(progress_path);

        let x_min = bounds.x_range.start;
        let x_max = bounds.x_range.end;

        let y_min = bounds.y_range.start;
        let y_max = bounds.y_range.end;

        let width = size.width;
        let height = size.height;

        writeln!(writer, "\\begin{{axis}}[hide axis,scale only axis,ticks=none,xmin={x_min},xmax={x_max},ymin={y_min},ymax={y_max},clip,width={width}cm,height={height}cm]")?;

        Ok(Self {
            name: name.to_owned(),
            writer,
            bounds,
            size,
            plot_count: 0,
            component,
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
            writeln!(self.writer, r#"\directlua{{progress_file:write(".")}}"#)?;
            writeln!(self.writer, r#"\directlua{{progress_file:flush()}}"#)?;
            self.plot_count += 1;
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
            pxu::CutType::UShortScallion(pxu::Component::Xp) => ("red", ""),
            pxu::CutType::UShortKidney(pxu::Component::Xp) => ("red", "dashed"),
            pxu::CutType::UShortScallion(pxu::Component::Xm) => ("green", ""),
            pxu::CutType::UShortKidney(pxu::Component::Xm) => ("green", "dashed"),
            _ => {
                return Ok(());
            }
        };

        self.add_plot(&[&[color, style, "thick"], options].concat(), &cut.path)?;

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

    fn finish(mut self, cache: Arc<cache::Cache>) -> std::io::Result<FigureCompiler> {
        writeln!(self.writer, "\\end{{axis}}\n")?;
        let component_name = match self.component {
            pxu::Component::P => "p",
            pxu::Component::Xp => "x^+",
            pxu::Component::Xm => "x^-",
            pxu::Component::U => "u",
        };
        writeln!(
            self.writer,
            "\\node at (current bounding box.north east) [anchor=north east,fill=white,outer sep=0.1cm,draw] {{$\\scriptstyle {component_name}$}};"
        )?;
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
    plot_count: u64,
}

impl FigureCompiler {
    fn new(figure: FigureWriter, cache: Arc<cache::Cache>) -> Result<Self> {
        let name = figure.name;
        if cache.check(&name)? {
            log::info!("[{name}]: Matches cached entry");
            let child = Command::new("/bin/true").spawn()?;
            Ok(Self {
                name,
                child,
                plot_count: 0,
            })
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

            Ok(Self {
                name,
                child,
                plot_count: figure.plot_count,
            })
        }
    }

    fn wait(mut self, pb: &ProgressBar) -> Result<String> {
        pb.set_length(self.plot_count);
        let mut path = PathBuf::from(FIGURE_PATH).join(&self.name);
        path.set_extension("prg");
        loop {
            if let Ok(meta) = path.metadata() {
                pb.set_position(meta.len());
            }

            if let Some(result) = self.child.try_wait()? {
                if result.success() {
                    log::info!("[{}]: Lualatex done.", self.name);
                } else {
                    log::error!("[{}]: Lualatex failed.", self.name);
                }
                break;
            }
            thread::sleep(std::time::Duration::from_millis(250));
        }
        let _ = std::fs::remove_file(path);
        Ok(self.name)
    }
}

#[derive(Debug, Default)]
struct Summary {
    names: Vec<String>,
}

impl Summary {
    const START: &str = r#"\nonstopmode
    \documentclass[12pt,a4paper]{article}
    \usepackage{graphicx}
    \usepackage{cprotect}
    \usepackage{caption}
    \captionsetup{labelformat=empty}
    \begin{document}
    \pagestyle{empty}
    "#;

    const END: &str = r#"\end{document}"#;

    fn add(&mut self, name: String) {
        self.names.push(name);
    }

    fn finish(self) -> Result<Child> {
        let mut path = PathBuf::from(FIGURE_PATH).join(SUMMARY_NAME);
        path.set_extension(TEX_EXT);

        let mut writer = BufWriter::new(File::create(path.clone())?);

        writer.write_all(Self::START.as_bytes())?;

        for name in self.names {
            writeln!(writer,"\\begin{{figure}}\\centering\\includegraphics{{{FIGURE_PATH}/{name}}}\\cprotect\\caption{{\\verb|\\includegraphics{{{FIGURE_PATH}/{name}}}|}}\\end{{figure}}")?;
        }

        writer.write_all(Self::END.as_bytes())?;

        writer.flush()?;

        let mut cmd = Command::new(LUALATEX);
        cmd.arg(format!("--output-directory={}", FIGURE_PATH))
            .args(["--interaction=nonstopmode", "--output-format=pdf"])
            .arg(path.as_os_str())
            .stderr(Stdio::null())
            .stdout(Stdio::null());

        log::info!("[{SUMMARY_NAME}]: Running Lualatex");
        cmd.spawn()
    }
}

trait Node {
    fn write_m_node(
        &mut self,
        figure: &mut FigureWriter,
        anchor: &str,
        rot_sign: i32,
        consts: CouplingConstants,
    ) -> Result<()>;
}

impl Node for PInterpolatorMut {
    fn write_m_node(
        &mut self,
        figure: &mut FigureWriter,
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

        writeln!(figure.writer,"\\node[scale=0.5,anchor={anchor},inner sep=0.4pt,rotate={rotation:.1},{color}] at ({:.3}, {:.3}) {{$\\scriptscriptstyle {}$}};",
                 self.p().re,
                 self.p().im,
                 m)
    }
}

fn fig_xpl_preimage(
    contours: Arc<Contours>,
    cache: Arc<cache::Cache>,
    consts: CouplingConstants,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xpL-preimage",
        Bounds::new(-2.6..2.6, -0.7..0.7),
        Size {
            width: 15.0,
            height: 6.0,
        },
        pxu::Component::P,
    )?;

    let state = State::new(1, consts);

    for contour in contours.get_grid(pxu::Component::P).iter() {
        figure.add_grid_line(contour, &[])?;
    }

    for cut in contours
        .get_visible_cuts(&state, pxu::Component::P, pxu::UCutType::Long)
        .filter(|cut| matches!(cut.typ, pxu::CutType::E))
    {
        figure.add_cut(cut, &[])?;
    }

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

    figure.finish(cache)
}

fn fig_xpl_cover(
    contours: Arc<Contours>,
    cache: Arc<cache::Cache>,
    _consts: CouplingConstants,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xpL-cover",
        Bounds::new(-6.0..6.0, -0.2..4.0),
        Size {
            width: 6.0,
            height: 4.0,
        },
        pxu::Component::Xp,
    )?;

    figure.add_axis()?;
    for contour in contours.get_grid(pxu::Component::Xp).iter().filter(
        |line| matches!(line.component, GridLineComponent::Xp(m) if (-8.0..=6.0).contains(&m)),
    ) {
        figure.add_grid_line(contour, &["thin", "black"])?;
    }
    figure.finish(cache)
}

fn fig_p_plane_long_cuts_regions(
    contours: Arc<Contours>,
    cache: Arc<cache::Cache>,
    consts: CouplingConstants,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "p_plane_long_cuts_regions",
        Bounds::new(-2.6..2.6, -0.7..0.7),
        Size {
            width: 15.0,
            height: 6.0,
        },
        pxu::Component::P,
    )?;

    let state = pxu::State::new(1, consts);

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

    for cut in contours.get_visible_cuts(&state, pxu::Component::P, pxu::UCutType::Long) {
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

    for contour in contours.get_grid(pxu::Component::P).iter() {
        figure.add_grid_line(contour, &[])?;
    }

    for cut in contours
        .get_visible_cuts(&state, pxu::Component::P, pxu::UCutType::Long)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::ULongNegative(_) | pxu::CutType::ULongPositive(_)
            )
        })
    {
        figure.add_cut(cut, &[])?;
    }

    figure.finish(cache)
}

fn fig_p_plane_short_cuts(
    contours: Arc<Contours>,
    cache: Arc<cache::Cache>,
    consts: CouplingConstants,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "p-plane-short-cuts",
        Bounds::new(-2.6..2.6, -0.7..0.7),
        Size {
            width: 25.0,
            height: 6.0 * 25.0 / 15.0,
        },
        pxu::Component::P,
    )?;

    let state = pxu::State::new(1, consts);

    for contour in contours.get_grid(pxu::Component::P).iter()
    // .take(100)
    {
        figure.add_grid_line(contour, &[])?;
    }

    for cut in contours
        .get_visible_cuts(&state, pxu::Component::P, pxu::UCutType::SemiShort)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::UShortScallion(_) | pxu::CutType::UShortKidney(_)
            )
        })
    {
        figure.add_cut(cut, &[])?;
    }

    figure.finish(cache)
}

fn fig_xp_cuts_1(
    contours: Arc<Contours>,
    cache: Arc<cache::Cache>,
    consts: CouplingConstants,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xp-cuts-1",
        Bounds::new(-4.0..4.0, -6.0..6.0),
        Size {
            width: 12.0,
            height: 18.0,
        },
        pxu::Component::Xp,
    )?;

    figure.add_axis()?;
    for contour in contours
        .get_grid(pxu::Component::Xp)
        .iter()
        .filter(|line| matches!(line.component, GridLineComponent::Xp(m) | GridLineComponent::Xm(m) if (-10.0..).contains(&m)))
    {
        figure.add_grid_line(contour, &["thin", "black"])?;
    }

    let mut state = pxu::State::new(1, consts);
    state.points[0].sheet_data.u_branch.1 = ::pxu::kinematics::UBranch::Between;

    for cut in contours
        .get_visible_cuts(&state, pxu::Component::Xp, pxu::UCutType::SemiShort)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::UShortScallion(_) | pxu::CutType::UShortKidney(_)
            )
        })
    {
        figure.add_cut(cut, &["very thick"])?;
    }

    figure.finish(cache)
}

fn main() -> std::io::Result<()> {
    // tracing_subscriber::fmt::fmt()
    //     .with_writer(std::io::stderr)
    //     .init();

    let spinner_style = ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
    )
    .unwrap()
    .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");

    let cache = cache::Cache::load(FIGURE_PATH)?;

    let consts = CouplingConstants::new(2.0, 5);

    let mut contours = pxu::Contours::new();

    println!("[1/3] Generating contours");
    let pb = ProgressBar::new(1);
    pb.set_style(spinner_style.clone());
    loop {
        pb.set_length(contours.progress().1 as u64);
        pb.set_position(contours.progress().0 as u64);
        if contours.update(0, consts) {
            pb.finish_and_clear();
            break;
        }
    }

    let fig_functions = [
        fig_xpl_preimage,
        fig_xpl_cover,
        fig_p_plane_long_cuts_regions,
        fig_p_plane_short_cuts,
        fig_xp_cuts_1,
    ];

    let contours_ref = Arc::new(contours);
    let cache_ref = Arc::new(cache);

    let mut handles = vec![];

    println!("[2/3] Builing figures");
    let mb = MultiProgress::new();

    for f in fig_functions {
        let contours_ref = contours_ref.clone();
        let cache_ref = cache_ref.clone();
        let pb = mb.add(ProgressBar::new_spinner());
        pb.set_style(spinner_style.clone());
        handles.push(thread::spawn(move || {
            pb.set_message("Generating tex file");
            let figure = f(contours_ref, cache_ref, consts)?;
            pb.set_message(format!("Compiling {}.tex", figure.name));
            let result = figure.wait(&pb);
            pb.finish_and_clear();
            result
        }));
    }

    let names = handles
        .into_iter()
        .map(|handle| {
            handle
                .join()
                .map_err(|err| error(format!("Join error: {err:?}").as_str()))?
        })
        .collect::<Result<Vec<_>>>()?;

    let mut new_cache = cache::Cache::new(FIGURE_PATH);
    let mut summary = Summary::default();

    for name in names {
        new_cache.update(&name)?;
        summary.add(name);
    }

    println!("[3/4] Saving cache");
    new_cache.save()?;

    println!("[4/4] Building summary");

    if summary.finish()?.wait()?.success() {
        log::info!("[{SUMMARY_NAME}] Done.");
    } else {
        log::error!("[{SUMMARY_NAME}] Error.");
    }

    Ok(())
}
