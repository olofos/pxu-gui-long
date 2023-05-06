use std::fs::File;
use std::io::{prelude::*, BufWriter, Result};
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};

use clap::Parser;
use indicatif::ProgressBar;

use crate::fig_compiler::FinishedFigure;

pub const TEX_EXT: &str = "tex";
pub const SUMMARY_NAME: &str = "all-figures";

pub fn error(message: &str) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::Other, message)
}

#[derive(Debug, Clone, PartialEq)]
pub struct Size {
    pub width: f64,
    pub height: f64,
}

#[derive(Parser, Clone)]
#[command(author, version, about, long_about = None)]
pub struct Settings {
    #[arg(short, long, default_value = "lualatex")]
    pub lualatex: String,
    #[arg(short, long, default_value = "./figures")]
    pub output_dir: String,
    #[arg(short, long)]
    pub rebuild: bool,
    #[arg(short, long, action = clap::ArgAction::Count)]
    pub verbose: u8,
    #[arg(short, long)]
    pub jobs: Option<usize>,
}

#[derive(Debug, Default)]
pub struct Summary {
    finished_figures: Vec<FinishedFigure>,
}

impl Summary {
    const START: &str = r#"\nonstopmode
    \documentclass[12pt,a4paper]{article}
    \usepackage[bottom=0.5cm, right=0.5cm, left=0.5cm, top=0.5cm]{geometry}
    \usepackage{graphicx}
    \usepackage{cprotect}
    \usepackage[justification=centering]{caption}
    \captionsetup{labelformat=empty}
    \usepackage{pdflscape}
    \pagestyle{empty}
    \begin{document}
    "#;

    const END: &str = r#"\end{document}"#;

    pub fn add(&mut self, finished_figure: FinishedFigure) {
        self.finished_figures.push(finished_figure);
    }

    pub fn finish(self, settings: &Settings, pb: &ProgressBar) -> Result<Child> {
        pb.set_message(format!("Creating {}.{}", SUMMARY_NAME, TEX_EXT));
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
            write!(writer, "\\cprotect\\caption{{")?;
            write!(writer, "\\verb|{includegraphics}|")?;
            if !finished_figure.caption.is_empty() {
                write!(writer, "\\\\{}", finished_figure.caption)?;
            }
            write!(writer, "}}\\end{{figure}}")?;

            if landscape {
                write!(writer, "\\end{{landscape}}")?;
            }

            writeln!(writer)?;
        }

        writer.write_all(Self::END.as_bytes())?;

        writer.flush()?;

        pb.set_message(format!("Compiling {}.{}", SUMMARY_NAME, TEX_EXT));

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
