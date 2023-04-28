use itertools::Itertools;
use num::complex::Complex64;
use pxu::GridLine;
use pxu::{
    interpolation::{InterpolationPoint, PInterpolatorMut},
    kinematics::CouplingConstants,
};
use std::fs::File;
use std::io::{prelude::*, BufWriter, Result};
use std::ops::Range;
use std::path::PathBuf;
use std::sync::Arc;

use crate::cache;
use crate::fig_compiler::FigureCompiler;
use crate::utils::{Settings, Size, TEX_EXT};

#[derive(Debug, Clone, PartialEq)]
pub struct Bounds {
    pub x_range: Range<f64>,
    pub y_range: Range<f64>,
}

impl Bounds {
    pub fn new(x_range: Range<f64>, y_range: Range<f64>) -> Self {
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

#[derive(Debug)]
pub struct FigureWriter {
    pub name: String,
    pub bounds: Bounds,
    pub size: Size,
    writer: BufWriter<File>,
    pub plot_count: u64,
    component: pxu::Component,
    y_shift: Option<f64>,
    u_cut_type: pxu::UCutType,
}

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

    pub fn new(
        name: &str,
        bounds: Bounds,
        size: Size,
        component: pxu::Component,
        u_cut_type: pxu::UCutType,
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
            u_cut_type,
        })
    }

    fn format_coordinate(&self, p: Complex64) -> String {
        format!(
            "({:.5},{:.5})",
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

    pub fn crop(&self, contour: &Vec<Complex64>) -> Vec<Complex64> {
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

    pub fn add_plot(&mut self, options: &[&str], contour: &Vec<Complex64>) -> Result<()> {
        self.add_plot_all(options, self.crop(contour))
    }

    pub fn add_plot_all(&mut self, options: &[&str], contour: Vec<Complex64>) -> Result<()> {
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

    pub fn add_grid_line(&mut self, grid_line: &GridLine, options: &[&str]) -> Result<()> {
        self.add_plot(&[&["very thin", "gray"], options].concat(), &grid_line.path)?;

        Ok(())
    }

    pub fn add_grid_lines(&mut self, pxu: &pxu::Pxu, options: &[&str]) -> Result<()> {
        for contour in pxu.contours.get_grid(self.component).iter() {
            self.add_grid_line(contour, options)?;
        }
        Ok(())
    }

    pub fn add_cut(
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

    pub fn add_cuts(&mut self, pxu: &pxu::Pxu, options: &[&str]) -> Result<()> {
        use pxu::{kinematics::UBranch, CutType::*, UCutType::*};

        let u_cut_type = self.u_cut_type;
        for cut in pxu
            .contours
            .get_visible_cuts(pxu, self.component, self.u_cut_type, 0)
            .filter(|cut| match cut.typ {
                Log(comp) | ULongPositive(comp) => {
                    u_cut_type != Short
                        || (comp == pxu::Component::Xp
                            && cut.component == pxu::Component::Xp
                            && pxu.state.points[0].sheet_data.u_branch.1 != UBranch::Between)
                        || (comp == pxu::Component::Xm
                            && cut.component == pxu::Component::Xm
                            && pxu.state.points[0].sheet_data.u_branch.0 != UBranch::Between)
                }
                ULongNegative(_) => u_cut_type == Long,
                UShortScallion(_) | UShortKidney(_) => u_cut_type != Long,
                E => true,
                DebugPath => false,
            })
        {
            self.add_cut(cut, options, pxu.consts)?;
        }
        Ok(())
    }

    pub fn add_axis(&mut self) -> Result<()> {
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

    pub fn finish(
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

pub trait Node {
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
