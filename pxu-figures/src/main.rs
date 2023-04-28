use itertools::Itertools;
use num::complex::Complex64;
use pxu::{
    interpolation::{InterpolationPoint, PInterpolatorMut},
    kinematics::{CouplingConstants, UBranch},
    Pxu,
};
use pxu::{GridLine, GridLineComponent};
use std::fs::File;
use std::io::{prelude::*, BufWriter, Result};
use std::ops::Range;
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};
use std::sync::Arc;
use std::thread;

use clap::Parser;
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

#[derive(Parser, Clone)]
#[command(author, version, about, long_about = None)]
struct Settings {
    #[arg(short, long, default_value = "lualatex")]
    lualatex: String,
    #[arg(short, long, default_value = "./figures")]
    output_dir: String,
    #[arg(short, long)]
    rebuild: bool,
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,
    #[arg(short, long)]
    jobs: Option<usize>,
}

#[derive(Debug)]
struct FigureWriter {
    name: String,
    bounds: Bounds,
    size: Size,
    writer: BufWriter<File>,
    plot_count: u64,
    component: pxu::Component,
    y_shift: Option<f64>,
}

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
\pgfplotsset{compat=1.17}
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
        settings: &Settings,
    ) -> std::io::Result<Self> {
        let mut path = PathBuf::from(&settings.output_dir).join(name);
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
            y_shift: None,
        })
    }

    fn format_coordinate(&self, p: Complex64) -> String {
        format!(
            "({:.3},{:.3})",
            p.re,
            p.im + self.y_shift.unwrap_or_default()
        )
    }

    fn format_contour(&self, contour: Vec<Complex64>) -> Vec<String> {
        let mut coordinates = contour
            .into_iter()
            .map(|z| self.format_coordinate(z))
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

    fn add_cut(
        &mut self,
        cut: &pxu::Cut,
        options: &[&str],
        consts: CouplingConstants,
    ) -> Result<()> {
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

        let shifts = if cut.component == pxu::Component::U && cut.periodic {
            let period = 2.0 * consts.k() as f64 / consts.h;
            (-5..=5).map(|n| Some(period * n as f64)).collect()
        } else {
            vec![None]
        };

        for shift in shifts {
            self.y_shift = shift;

            self.add_plot(&[&[color, style, "thick"], options].concat(), &cut.path)?;

            if let Some(branch_point) = cut.branch_point {
                self.add_plot_all(
                    &[&[color, "only marks", "mark size=0.02cm"], options].concat(),
                    vec![branch_point],
                )?;
            }
        }

        self.y_shift = None;

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

    fn finish(
        mut self,
        cache: Arc<cache::Cache>,
        settings: &Settings,
    ) -> std::io::Result<FigureCompiler> {
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

        FigureCompiler::new(self, cache, settings)
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
    size: Size,
}

#[derive(Debug)]
struct FinishedFigure {
    name: String,
    size: Size,
}

impl FigureCompiler {
    fn new(figure: FigureWriter, cache: Arc<cache::Cache>, settings: &Settings) -> Result<Self> {
        let FigureWriter {
            name,
            size,
            plot_count,
            ..
        } = figure;
        if !settings.rebuild && cache.check(&name)? {
            log::info!("[{name}]: Matches cached entry");
            let child = Command::new("/bin/true").spawn()?;
            Ok(Self {
                name,
                child,
                plot_count: 0,
                size,
            })
        } else {
            let mut path = PathBuf::from(&settings.output_dir).join(name.clone());
            path.set_extension(TEX_EXT);

            let mut cmd = Command::new(&settings.lualatex);
            cmd.arg(format!("--output-directory={}", settings.output_dir))
                .args(["--interaction=nonstopmode", "--output-format=pdf"])
                .arg(path.as_os_str())
                .stderr(Stdio::null())
                .stdout(Stdio::null());

            log::info!("[{name}]: Running Lualatex");
            let child = cmd.spawn()?;

            Ok(Self {
                name,
                child,
                plot_count,
                size,
            })
        }
    }

    fn wait(mut self, pb: &ProgressBar, settings: &Settings) -> Result<FinishedFigure> {
        pb.set_length(self.plot_count + 1);
        let mut path = PathBuf::from(&settings.output_dir).join(&self.name);
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
        Ok(FinishedFigure {
            name: self.name,
            size: self.size,
        })
    }
}

#[derive(Debug, Default)]
struct Summary {
    finished_figures: Vec<FinishedFigure>,
}

impl Summary {
    const START: &str = r#"\nonstopmode
    \documentclass[12pt,a4paper]{article}
    \usepackage{graphicx}
    \usepackage{cprotect}
    \usepackage{caption}
    \captionsetup{labelformat=empty}
    \usepackage{pdflscape}
    \begin{document}
    \pagestyle{empty}
    "#;

    const END: &str = r#"\end{document}"#;

    fn add(&mut self, finished_figure: FinishedFigure) {
        self.finished_figures.push(finished_figure);
    }

    fn finish(self, settings: &Settings) -> Result<Child> {
        let mut path = PathBuf::from(&settings.output_dir).join(SUMMARY_NAME);
        path.set_extension(TEX_EXT);

        let mut writer = BufWriter::new(File::create(path.clone())?);

        writer.write_all(Self::START.as_bytes())?;

        let output_dir = &settings.output_dir;

        for finished_figure in self.finished_figures {
            let name = &finished_figure.name;
            let Size { width, height } = finished_figure.size;

            let landscape = width > 20.0;

            if landscape {
                write!(writer, "\\begin{{landscape}}")?;
            }

            let includegraphics = format!(
                "\\includegraphics[width={width}cm,height={height}cm]{{{output_dir}/{name}}}"
            );
            write!(writer, "\\begin{{figure}}\\centering")?;
            write!(writer, "{includegraphics}")?;
            write!(writer, "\\cprotect\\caption{{\\verb|")?;
            write!(writer, "{includegraphics}")?;
            write!(writer, "|}}\\end{{figure}}")?;

            if landscape {
                write!(writer, "\\end{{landscape}}")?;
            }

            writeln!(writer)?;
        }

        writer.write_all(Self::END.as_bytes())?;

        writer.flush()?;

        let mut cmd = Command::new(&settings.lualatex);
        cmd.arg(format!("--output-directory={}", settings.output_dir))
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
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xpL-preimage",
        Bounds::new(-2.6..2.6, -0.7..0.7),
        Size {
            width: 15.0,
            height: 6.0,
        },
        pxu::Component::P,
        settings,
    )?;

    for contour in pxu.contours.get_grid(pxu::Component::P).iter() {
        figure.add_grid_line(contour, &[])?;
    }

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
        Bounds::new(-6.0..6.0, -0.2..4.0),
        Size {
            width: 6.0,
            height: 4.0,
        },
        pxu::Component::Xp,
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
        Bounds::new(-2.6..2.6, -0.7..0.7),
        Size {
            width: 15.0,
            height: 6.0,
        },
        pxu::Component::P,
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

    for contour in pxu.contours.get_grid(pxu::Component::P).iter() {
        figure.add_grid_line(contour, &[])?;
    }

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::P, pxu::UCutType::Long, 0)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::ULongNegative(_) | pxu::CutType::ULongPositive(_)
            )
        })
    {
        figure.add_cut(cut, &[], pxu.consts)?;
    }

    figure.finish(cache, settings)
}

fn fig_p_plane_short_cuts(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "p-plane-short-cuts",
        Bounds::new(-2.6..2.6, -0.7..0.7),
        Size {
            width: 25.0,
            height: 6.0 * 25.0 / 15.0,
        },
        pxu::Component::P,
        settings,
    )?;

    for contour in pxu.contours.get_grid(pxu::Component::P).iter()
    // .take(100)
    {
        figure.add_grid_line(contour, &[])?;
    }

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::P, pxu::UCutType::SemiShort, 0)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::UShortScallion(_) | pxu::CutType::UShortKidney(_)
            )
        })
    {
        figure.add_cut(cut, &[], pxu.consts)?;
    }

    figure.finish(cache, settings)
}

fn fig_xp_cuts_1(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xp-cuts-1",
        Bounds::new(-4.0..4.0, -6.0..6.0),
        Size {
            width: 12.0,
            height: 18.0,
        },
        pxu::Component::Xp,
        settings,
    )?;

    figure.add_axis()?;
    for contour in pxu.contours
        .get_grid(pxu::Component::Xp)
        .iter()
        .filter(|line| matches!(line.component, GridLineComponent::Xp(m) | GridLineComponent::Xm(m) if (-10.0..).contains(&m)))
    {
        figure.add_grid_line(contour, &["thin", "black"])?;
    }

    let mut state = pxu::State::new(1, pxu.consts);
    state.points[0].sheet_data.u_branch.1 = ::pxu::kinematics::UBranch::Between;

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::Xp, pxu::UCutType::SemiShort, 0)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::UShortScallion(_) | pxu::CutType::UShortKidney(_)
            )
        })
    {
        figure.add_cut(cut, &["very thick"], pxu.consts)?;
    }

    figure.finish(cache, settings)
}

fn fig_u_period_between_between(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "u-period-between-between",
        Bounds::new(-6.0..4.0, -12.25..12.75),
        Size {
            width: 5.0,
            height: 12.5,
        },
        pxu::Component::U,
        settings,
    )?;

    for contour in pxu.contours.get_grid(pxu::Component::U).iter()
    // .filter(|line| matches!(line.component, GridLineComponent::Xp(m) | GridLineComponent::Xm(m) if (-10.0..).contains(&m)))
    {
        figure.add_grid_line(contour, &["very thin", "gray"])?;
    }

    let mut pxu = (*pxu).clone();
    pxu.state.points[0].sheet_data.u_branch = (
        ::pxu::kinematics::UBranch::Between,
        ::pxu::kinematics::UBranch::Between,
    );

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::U, pxu::UCutType::Short, 0)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::UShortScallion(_) | pxu::CutType::UShortKidney(_)
            )
        })
    {
        figure.add_cut(cut, &["very thick"], pxu.consts)?;
    }

    let path = pxu
        .get_path_by_name("U period between/between")
        .ok_or(error("Path not found"))?;

    for segment in &path.segments[0] {
        figure.add_plot(&["very thick", "blue"], &segment.u)?;
    }

    figure.finish(cache, settings)
}

fn fig_u_band_between_outside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "u-band-between-outside",
        Bounds::new(-6.0..4.0, -12.25..12.75),
        Size {
            width: 5.0,
            height: 12.5,
        },
        pxu::Component::U,
        settings,
    )?;

    for contour in pxu.contours.get_grid(pxu::Component::U).iter()
    // .filter(|line| matches!(line.component, GridLineComponent::Xp(m) | GridLineComponent::Xm(m) if (-10.0..).contains(&m)))
    {
        figure.add_grid_line(contour, &["very thin", "gray"])?;
    }

    let mut pxu = (*pxu).clone();
    pxu.state.points[0].sheet_data.u_branch = (
        ::pxu::kinematics::UBranch::Between,
        ::pxu::kinematics::UBranch::Outside,
    );

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::U, pxu::UCutType::Short, 0)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::UShortScallion(_) | pxu::CutType::UShortKidney(_)
            )
        })
    {
        figure.add_cut(cut, &["very thick"], pxu.consts)?;
    }

    let path = pxu
        .get_path_by_name("U band between/outside")
        .ok_or(error("Path not found"))?;

    for segment in &path.segments[0] {
        figure.add_plot(&["very thick", "blue"], &segment.u)?;
    }

    figure.finish(cache, settings)
}

fn fig_u_band_between_inside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "u-band-between-inside",
        Bounds::new(-6.0..4.0, -12.25..12.75),
        Size {
            width: 5.0,
            height: 12.5,
        },
        pxu::Component::U,
        settings,
    )?;

    for contour in pxu.contours.get_grid(pxu::Component::U).iter()
    // .filter(|line| matches!(line.component, GridLineComponent::Xp(m) | GridLineComponent::Xm(m) if (-10.0..).contains(&m)))
    {
        figure.add_grid_line(contour, &["very thin", "gray"])?;
    }

    let mut pxu = (*pxu).clone();
    pxu.state.points[0].sheet_data.u_branch = (
        ::pxu::kinematics::UBranch::Between,
        ::pxu::kinematics::UBranch::Inside,
    );

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::U, pxu::UCutType::Short, 0)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::UShortScallion(_) | pxu::CutType::UShortKidney(_)
            )
        })
    {
        figure.add_cut(cut, &["very thick"], pxu.consts)?;
    }

    let path = pxu
        .get_path_by_name("U band between/inside")
        .ok_or(error("Path not found"))?;

    for segment in &path.segments[0] {
        figure.add_plot(&["very thick", "blue"], &segment.u)?;
    }

    figure.finish(cache, settings)
}

fn fig_p_band_between_outside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "p-band-between-outside",
        Bounds::new(-2.6..2.6, -0.7..0.7),
        Size {
            width: 15.0,
            height: 6.0,
        },
        pxu::Component::P,
        settings,
    )?;

    for contour in pxu.contours.get_grid(pxu::Component::P).iter()
    // .filter(|line| matches!(line.component, GridLineComponent::Xp(m) | GridLineComponent::Xm(m) if (-10.0..).contains(&m)))
    {
        figure.add_grid_line(contour, &["very thin", "gray"])?;
    }

    let path = pxu
        .get_path_by_name("U band between/outside")
        .ok_or(error("Path not found"))?;

    let mut pxu = (*pxu).clone();
    pxu.state = path.base_path.start.clone();

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::P, pxu::UCutType::Short, 0)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::UShortScallion(_) | pxu::CutType::UShortKidney(_)
            )
        })
    {
        figure.add_cut(cut, &["very thick"], pxu.consts)?;
    }

    let mut straight_segments = vec![];
    let mut dotted_segments = vec![];

    let mut same_branch = false;
    let mut points = vec![];

    for segment in &path.segments[0] {
        let segment_same_branch = segment.sheet_data.is_same(
            &pxu.state.points[0].sheet_data,
            pxu::Component::P,
            pxu::UCutType::Short,
        );

        if segment_same_branch != same_branch && !points.is_empty() {
            if same_branch {
                straight_segments.push(points);
            } else {
                dotted_segments.push(points);
            }
            points = vec![];
        }

        points.extend(&segment.p);
        same_branch = segment_same_branch;
    }

    if same_branch {
        straight_segments.push(points);
    } else {
        dotted_segments.push(points);
    }

    for points in dotted_segments {
        figure.add_plot(&["very thick", "blue", "dotted"], &points)?;
    }

    for points in straight_segments {
        figure.add_plot(&["very thick", "blue"], &points)?;
    }

    figure.finish(cache, settings)
}

fn fig_p_band_between_inside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "p-band-between-inside",
        Bounds::new(-2.6..2.6, -0.7..0.7),
        Size {
            width: 15.0,
            height: 6.0,
        },
        pxu::Component::P,
        settings,
    )?;

    for contour in pxu.contours.get_grid(pxu::Component::P).iter()
    // .filter(|line| matches!(line.component, GridLineComponent::Xp(m) | GridLineComponent::Xm(m) if (-10.0..).contains(&m)))
    {
        figure.add_grid_line(contour, &["very thin", "gray"])?;
    }

    let path = pxu
        .get_path_by_name("U band between/inside")
        .ok_or(error("Path not found"))?;

    let mut pxu = (*pxu).clone();
    pxu.state = path.base_path.start.clone();

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::P, pxu::UCutType::Short, 0)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::UShortScallion(_) | pxu::CutType::UShortKidney(_)
            )
        })
    {
        figure.add_cut(cut, &["very thick"], pxu.consts)?;
    }

    let mut straight_segments = vec![];
    let mut dotted_segments = vec![];

    let mut same_branch = false;
    let mut points = vec![];

    for segment in &path.segments[0] {
        let segment_same_branch = segment.sheet_data.is_same(
            &pxu.state.points[0].sheet_data,
            pxu::Component::P,
            pxu::UCutType::Short,
        );

        if segment_same_branch != same_branch && !points.is_empty() {
            if same_branch {
                straight_segments.push(points);
            } else {
                dotted_segments.push(points);
            }
            points = vec![];
        }

        points.extend(&segment.p);
        same_branch = segment_same_branch;
    }

    if same_branch {
        straight_segments.push(points);
    } else {
        dotted_segments.push(points);
    }

    for points in dotted_segments {
        figure.add_plot(&["very thick", "blue", "dotted"], &points)?;
    }

    for points in straight_segments {
        figure.add_plot(&["very thick", "blue"], &points)?;
    }

    figure.finish(cache, settings)
}

fn fig_xp_band_between_inside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xp-band-between-inside",
        Bounds::new(-2.6..2.6, -2.6..2.6),
        Size {
            width: 8.0,
            height: 8.0,
        },
        pxu::Component::Xp,
        settings,
    )?;

    for contour in pxu.contours.get_grid(pxu::Component::Xp).iter()
    // .filter(|line| matches!(line.component, GridLineComponent::Xp(m) | GridLineComponent::Xm(m) if (-10.0..).contains(&m)))
    {
        figure.add_grid_line(contour, &["very thin", "gray"])?;
    }

    let path = pxu
        .get_path_by_name("U band between/inside")
        .ok_or(error("Path not found"))?;

    let mut pxu = (*pxu).clone();
    pxu.state.points[0].sheet_data.u_branch = (UBranch::Between, UBranch::Inside);
    pxu.state.points[0].sheet_data.log_branch_p = 0;
    pxu.state.points[0].sheet_data.log_branch_m = -1;
    pxu.state.points[0].sheet_data.im_x_sign = (1, -1);

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::Xp, pxu::UCutType::Short, 0)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E
                    | pxu::CutType::UShortScallion(_)
                    | pxu::CutType::UShortKidney(_)
                    | pxu::CutType::Log(pxu::Component::Xp)
                    | pxu::CutType::ULongPositive(pxu::Component::Xp)
            )
        })
    {
        figure.add_cut(cut, &["very thick"], pxu.consts)?;
    }

    let mut straight_segments = vec![];
    let mut dotted_segments = vec![];

    let mut same_branch = false;
    let mut points = vec![];

    for segment in &path.segments[0] {
        let segment_same_branch = segment.sheet_data.is_same(
            &pxu.state.points[0].sheet_data,
            pxu::Component::Xp,
            pxu::UCutType::Short,
        );

        if segment_same_branch != same_branch && !points.is_empty() {
            if same_branch {
                straight_segments.push(points);
            } else {
                dotted_segments.push(points);
            }
            points = vec![];
        }

        points.extend(&segment.xp);
        same_branch = segment_same_branch;
    }

    if same_branch {
        straight_segments.push(points);
    } else {
        dotted_segments.push(points);
    }

    for points in dotted_segments {
        figure.add_plot(&["very thick", "blue", "dotted"], &points)?;
    }

    for points in straight_segments {
        figure.add_plot(&["very thick", "blue"], &points)?;
    }

    figure.finish(cache, settings)
}

fn fig_xp_band_between_outside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xp-band-between-outside",
        Bounds::new(-3.1..2.1, -2.6..2.6),
        Size {
            width: 8.0,
            height: 8.0,
        },
        pxu::Component::Xp,
        settings,
    )?;

    for contour in pxu.contours.get_grid(pxu::Component::Xp).iter() {
        figure.add_grid_line(contour, &["very thin", "gray"])?;
    }

    let path = pxu
        .get_path_by_name("U band between/outside")
        .ok_or(error("Path not found"))?;

    let mut pxu = (*pxu).clone();
    pxu.state.points[0].sheet_data.u_branch = (UBranch::Between, UBranch::Outside);
    pxu.state.points[0].sheet_data.log_branch_p = 0;
    pxu.state.points[0].sheet_data.log_branch_m = -1;
    pxu.state.points[0].sheet_data.im_x_sign = (1, -1);

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::Xp, pxu::UCutType::Short, 0)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::UShortScallion(_) | pxu::CutType::UShortKidney(_)
            )
        })
    {
        figure.add_cut(cut, &["very thick"], pxu.consts)?;
    }

    let mut straight_segments = vec![];
    let mut dotted_segments = vec![];

    let mut same_branch = false;
    let mut points = vec![];

    for segment in &path.segments[0] {
        let segment_same_branch = segment.sheet_data.is_same(
            &pxu.state.points[0].sheet_data,
            pxu::Component::Xp,
            pxu::UCutType::Short,
        );

        if segment_same_branch != same_branch && !points.is_empty() {
            if same_branch {
                straight_segments.push(points);
            } else {
                dotted_segments.push(points);
            }
            points = vec![];
        }

        points.extend(&segment.xp);
        same_branch = segment_same_branch;
    }

    if same_branch {
        straight_segments.push(points);
    } else {
        dotted_segments.push(points);
    }

    for points in dotted_segments {
        figure.add_plot(&["very thick", "blue", "dotted"], &points)?;
    }

    for points in straight_segments {
        figure.add_plot(&["very thick", "blue"], &points)?;
    }

    figure.finish(cache, settings)
}

fn fig_xm_band_between_inside(
    pxu: Arc<Pxu>,
    cache: Arc<cache::Cache>,
    settings: &Settings,
) -> Result<FigureCompiler> {
    let mut figure = FigureWriter::new(
        "xm-band-between-inside",
        Bounds::new(-1.3..1.3, -1.3..1.3),
        Size {
            width: 8.0,
            height: 8.0,
        },
        pxu::Component::Xm,
        settings,
    )?;

    for contour in pxu.contours.get_grid(pxu::Component::Xm).iter() {
        figure.add_grid_line(contour, &["very thin", "gray"])?;
    }

    let path = pxu
        .get_path_by_name("U band between/inside")
        .ok_or(error("Path not found"))?;

    let mut pxu = (*pxu).clone();
    pxu.state.points[0].sheet_data.u_branch = (UBranch::Between, UBranch::Inside);
    pxu.state.points[0].sheet_data.log_branch_p = 0;
    pxu.state.points[0].sheet_data.log_branch_m = -1;
    pxu.state.points[0].sheet_data.im_x_sign = (1, -1);

    for cut in pxu
        .contours
        .get_visible_cuts(&pxu, pxu::Component::Xm, pxu::UCutType::Short, 0)
        .filter(|cut| {
            matches!(
                cut.typ,
                pxu::CutType::E | pxu::CutType::UShortScallion(_) | pxu::CutType::UShortKidney(_)
            )
        })
    {
        figure.add_cut(cut, &["very thick"], pxu.consts)?;
    }

    let mut straight_segments = vec![];
    let mut dotted_segments = vec![];

    let mut same_branch = false;
    let mut points = vec![];

    for segment in &path.segments[0] {
        let segment_same_branch = segment.sheet_data.is_same(
            &pxu.state.points[0].sheet_data,
            pxu::Component::Xm,
            pxu::UCutType::Short,
        );

        if segment_same_branch != same_branch && !points.is_empty() {
            if same_branch {
                straight_segments.push(points);
            } else {
                dotted_segments.push(points);
            }
            points = vec![];
        }

        points.extend(&segment.xm);
        same_branch = segment_same_branch;
    }

    if same_branch {
        straight_segments.push(points);
    } else {
        dotted_segments.push(points);
    }

    for points in dotted_segments {
        figure.add_plot(&["very thick", "blue", "dotted"], &points)?;
    }

    for points in straight_segments {
        figure.add_plot(&["very thick", "blue"], &points)?;
    }

    figure.finish(cache, settings)
}

fn main() -> std::io::Result<()> {
    let settings = Settings::parse();

    if settings.verbose > 0 {
        tracing_subscriber::fmt::fmt()
            .with_writer(std::io::stderr)
            .init();
        log::set_max_level(log::LevelFilter::Debug);
    }

    if settings.rebuild {
        println!(" ---  Rebuilding all figures");
    }

    let spinner_style = ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
    )
    .unwrap()
    .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");

    let cache = cache::Cache::load(&settings.output_dir)?;

    let consts = CouplingConstants::new(2.0, 5);

    let mut contours = pxu::Contours::new();

    let pb = if settings.verbose == 0 {
        println!("[1/5] Generating contours");
        ProgressBar::new(1)
    } else {
        ProgressBar::hidden()
    };
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
        fig_u_period_between_between,
        fig_u_band_between_outside,
        fig_u_band_between_inside,
        fig_p_band_between_outside,
        fig_p_band_between_inside,
        fig_xp_band_between_inside,
        fig_xp_band_between_outside,
        fig_xm_band_between_inside,
    ];

    let mut pxu = Pxu::new(consts);
    pxu.contours = contours;
    pxu.state = pxu::State::new(1, consts);

    let pb = if settings.verbose == 0 {
        println!("[2/5] Loading paths");
        ProgressBar::new(1)
    } else {
        ProgressBar::hidden()
    };

    let saved_paths = pxu::path::SavedPath::load("7VxLayTJEf4rYk8SVLTz_dDR4INP3suCYRkGjaaZER61hNTC8_P9RT4qI7WSGQZjWEhdOrMrKzIeX3xRWZmt3y-fzzdP5-vLx4e70_n5-vfLx-tLUgeVrbZK6RhdtNFZs_GXynqtrTM6ZG21dslfbd_5Bn3wiv9ccNpEm6Pb1MHl8ueUMdFpHy0G35fBLluTrQ1BZa2cSiwcw0IZG1SKyrLoFww2B8j0iYeZoLXezMH7gF420EMHa9PV9vz1eDx__Hxzvrm-_Pbw5eOnp5vT7dePj9ekN9G_v1bbsXWu9fbSm5d_PZ7_fTyetvZ5td3df_z-8fnuy4nV3fTV1dWHq-3z8dv5hl2knd_Ip4SZL5213I4O7Wg0t11AWysVuaMsd0xW6MCX3PGRh4VYOsnBeu8zDzPa8DCXPXdsYgHFEZcmlGHwAHey4jm14StWe0hzOfKk1mbuJF2uBOPQCY5F2-QgzbmgWGcVeJgJfA9_oKOKaB6wEQ_mDm5FB2K4A6HoYALuYDp0MDU6rMhGrBR3DItmdbkD5dGBIdyBWejARO7A4I3YeO7AFejALdyJrCg7jDuJLWVXlo7DpPBxaUe2zZYpEQgMcr63exND-gjc2W-ExC6Qp-ozsQ5dBVau68Zad6XZnG4N29nNZAc0-4tnmmOKy7rH2Jfdlezk7mP2fnc-h6VFpcSrhasEsseRI9wDzKFvkS-YaJAoYOlYYRR1EDG8GroK7hrsCiA7HhmpHagM4YZghnZDNiO-AZ4ToeUBybYYQuJOkhJJTkVSB5LKkdSapDkk7STpAJKeIekykr4k6WSS3icZFpLxIhlIkhEmGXqSmCAJFpIoIgkvkrgjCUiSSCUJYRLQJoF4EolAMkFIZg7JlKIp2WjKQ5pSlKbspSmxacp5muiAJqagiURo4heaqIcmVqKJsGjiMppojiYGpIkcaeJNmiiVJraliYhp4mia6JsmZqeJ9GmqByQrBckSQrK2rDqz6syqM6vOrDqz6syqM6vOrDqz6syqM6vOrDqz6syqM6vOrDqz6syqM6_rzIft9uH-8eF0PJ2v__m4Hb_f3p1vzncPJ-z-3D6cnrHLdfn12hzU9q9rf4DBp5v74_Uv3x8vbu-ebr8dLz7VXaC_tM9fIP6NTTJ10NjjUnUDzByCy9iLiiljgwy6BFxLvL2VjIJWEXtauW5_vTGU3hr7wkN19MFrnyyGWWN9nfHd7S71Q7td_3g5P999Pm7tc9rt-uNmF-IBStAVF-CDzC3gylWkJqROqhHAPp6L5TLsQrPFlnmjshhyAVt5pZkzCKSMNQips-Vb7OiBSWKFkt1cRZxJaKqqgDJbyxCL22w0FaBps76i0zOvFGG4uNmaoQxE3FphD0ZJNY-QOJimpkbcMHkFVdqMafyVNyhamh6sEBuRBfBLza8IQrF1bMRYVUkcLkmN5uLWSDJlkF_NXrjY9Ba1Ji5Tu85J029CFjZJEE9dfmQSadNySndtPLNVU9Lx-KY7LKJuEgylbikTR3eAYrpJtnuLmrvYidS96AvRpO5y6j7nnOyhUDy-RojjRi1wHE7q8USUqYeZ723RZ0xQAwVDhTpWfCGj0IFFHVm8b9wAx55qMGTKa9i0zHcVsJrnLCgerXF53DMECfFiVqGM0FGoLiwShgr7hVuEt4QThW-Fy0UkRIBE3EQ4RZRF8AUmBFQEggSwBN4GCgc0B14HiHdg71jf4b9nxMiSkTojn0aSjcwb6ThydCTuyOaR4iPvBxkMhhi0MbhkEMxgnUFFg58GaQ0mG_Q2OG8Q4U6OO1_uFNpYVRarX3-4Vv3aS9XN08PL6fPFw9Pdl7vTxfkrul--Xvzt4hbU_m7dUn4ULocSk70yzvnonTGJr4Wgc0C99xH1SDtdC9cbQ-mtsShc9oBHgey0Sco7b1L29fjI9Ke1VsYgDC4Y1Ljg_9-FLQBP7Wmkt6xCq7Kp6RUBETa2Vj1-AC3Mgvib-iiYGJ_lKtBh6nMXDMPDaeMdIFS7zl2mPiTCT2jqym7QrD0GASU6VZ4LSIBY6S8iK0IVlpAqtTZofuStT4AQtBXHczP2goRJkY6tCWFaVb51W25lN2910QFDtlo1YNwW2nohb7Wq1ESoLYcUqy0PRXsr9AaNVr_qmKgE6XfZ1IRbZhTftaCmBlSjphs_iDaFNVNSMwPGUbeOnw-70UyHzRe6PGtXF7HjqHsuMvc1hwaW2fzMRaC7n59xe1QcU2wLli3t2CNLPbQIOLWIcyloOEiFe3MDDDXEhPKc35Yx3MwNb70wWxZWUahZQoVm6K29MS6OO4aYIXtMOLQYqg19hRXCOGGzcIXwkHCc8Kdws_C-CIqIlQihiKwI-IDBwMYAzEDRgNbA2wDhQOaO1o7fHdE7yHfc76mwZ8eeMHsO7Wk1Um3k30jKkakjfUdOj0Qf2T8oYfDEII_BKINmBvfsfLRT1M5aO5Ht3LbT3c6AOylWmvy5QrWvqVqFOj2c5yr1TpGyB5VUsKi2AYY7lcsiyasUrNZJoboan91--BB1xClvssVCMigPf-P-ZLCg1Nbw6olXq-XsIaqQxuI1ahWxitTwsrFcDT0ERos7rIdkHv5SJMdyhFErbwPKlGEcH2I966hxT-TS_27Vsu9WLXrr9OHfT29ULXpdtmzgybeCGGO8wSOWqQiwQFxS5YJOHP2cCkpQnPH4oUvkAQ88nvlGIPAZ7reKAZaTB6xtgVPmQfAMQ5A_Nv6CAYTLGw9kWOC2rUgor_U4gXQdxVNtZdbygi0BuKxPwRQ0L5qWZyKzFQvqU7voTMNokkCTcJrmpUklksqStIKkeSTtpskjNDmLJj_S5GKavE8yLiTiNV94dc8s79Vckx6TjpP6k2GTyZMzZj_NPpz9O_t-josMmWxPoyYBk-xp2kkjqas0QlonzZ4cMgN8hr5MiildRFhWeq30Wum10mul10qvlV4rvVZ6rfRa6bXSa6XXSq-VXiu9_rTpJV_g_vbDL3B_u_h0gw3GfiTmrryifPelLfY_seeHl7F4OYsNPT4Eg93CrLAHqGJW2eG7-tIWr1wddhEj3rhqaGg1tNQHrU02OuBdK9-STD8xYxV8oqO13jlly1v8g1E2JIN3tnBO2cd9KVLLX_bWRDggYwMSyujy2rb9HJ1_kY5Xv-6n3tu-9dr2rd3GN340buKeLAb7yfxWXdX9Bd7tjlrXk6x4J11wVPcGMt7Lo1Mjzp9b-a7Fn9_5YzzjCrdvRVIRr_htvk1tP072poE0C6F5Bpqnp1k3mvSmySRswQximK-8umsW-WrC1-rMys6WTFZOnWncLOOVR6fZZ9XmIEmLhKkrrCusK6wrrCusK6wrrCusK6xvhfV_sg55qM_c7yxE-JSHCXj2z1ZrLDCU8fUAiQrJZ2-88irH6PpSRIfyH6d4qAt8flQfLO7yEGCCwprFiuMjNlmdIFYhYMalwAcmUw4pp2SSzVjT8Cmdlzq4rkawrlNYlIUQgzN8TFf7Q1uKaJUVH_j_8dUIhZ_6H1b0X9YjOsMfm86mZlD0WGWGelg6Giw5s3HtzJfnQ02unO4tn1v5jhetGLKV0Zy3gc8Cs5wi3PDJJcxQT-XJ3jSQJhk0y6d5cpoUo0lnmsyZEwvxdJtJqSZW0harbl8lAiO8hAdOa65iaVnOXrXpubGVL-tBwMQ_TeJbyvsIyNiKuCo3IZd4Kka-7EzjZhnzBPPss2qT2pNFwtQV0hXSFdIV0hXSFdIV0hXSFdJtPPH-7BLk8fh09_D5j78P_vAf").unwrap();
    pb.set_style(spinner_style.clone());
    pb.set_length(saved_paths.len() as u64);

    pxu.paths = saved_paths
        .into_iter()
        .map(|saved_path| {
            pb.set_message(saved_path.name.clone());
            let path =
                pxu::path::Path::from_base_path(saved_path.into(), &pxu.contours, pxu.consts);
            pb.inc(1);
            path
        })
        .collect();
    pb.finish_and_clear();

    let pxu_ref = Arc::new(pxu);
    let cache_ref = Arc::new(cache);

    if settings.verbose == 0 {
        if settings.rebuild {
            println!("[3/5] Builing figures (ignoring cache)");
        } else {
            println!("[3/5] Builing figures");
        }
    }
    let mb = Arc::new(MultiProgress::new());

    let num_threads = if let Some(jobs) = settings.jobs {
        jobs
    } else {
        num_cpus::get()
    }
    .min(fig_functions.len());

    let pool = threadpool::ThreadPool::new(num_threads);
    let (tx, rx) = std::sync::mpsc::channel();

    for (i, f) in fig_functions.into_iter().enumerate() {
        let pxu_ref = pxu_ref.clone();
        let cache_ref = cache_ref.clone();
        let spinner_style = spinner_style.clone();
        let settings = settings.clone();
        let mb = mb.clone();
        let tx = tx.clone();
        pool.execute(move || {
            let pb = if settings.verbose == 0 {
                mb.add(ProgressBar::new_spinner())
            } else {
                ProgressBar::hidden()
            };
            pb.set_style(spinner_style);
            pb.set_message("Generating tex file");
            match f(pxu_ref, cache_ref, &settings) {
                Ok(figure) => {
                    pb.set_message(format!("Compiling {}.tex", figure.name));
                    let result = figure.wait(&pb, &settings);
                    pb.finish_and_clear();
                    tx.send(result.map(|r| (i, r))).unwrap();
                }
                Err(e) => {
                    tx.send(Err(e)).unwrap();
                }
            }
        });
    }

    pool.join();

    let mut finished_figures = rx
        .into_iter()
        .take(fig_functions.len())
        .collect::<Result<Vec<_>>>()?;
    finished_figures.sort_by_key(|&(n, _)| n);
    let finished_figures = finished_figures.into_iter().map(|(_, r)| r);

    let mut new_cache = cache::Cache::new(&settings.output_dir);
    let mut summary = Summary::default();

    for finished_figure in finished_figures {
        new_cache.update(&finished_figure.name)?;
        summary.add(finished_figure);
    }

    if settings.verbose == 0 {
        println!("[4/5] Saving cache");
    }
    new_cache.save()?;

    if settings.verbose == 0 {
        println!("[5/5] Building summary");
    }

    if summary.finish(&settings)?.wait()?.success() {
        log::info!("[{SUMMARY_NAME}] Done.");
    } else {
        log::error!("[{SUMMARY_NAME}] Error.");
    }

    Ok(())
}
