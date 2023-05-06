use std::io::Result;
use std::path::PathBuf;
use std::process::{Child, Command, Stdio};
use std::sync::Arc;
use std::thread;

use indicatif::ProgressBar;

use crate::cache;
use crate::fig_writer::FigureWriter;
use crate::utils::{Settings, Size, TEX_EXT};

pub struct FigureCompiler {
    pub name: String,
    pub caption: String,
    child: Child,
    plot_count: u64,
    size: Size,
}

#[derive(Debug)]
pub struct FinishedFigure {
    pub name: String,
    pub caption: String,
    pub size: Size,
}

impl FigureCompiler {
    pub fn new(
        figure: FigureWriter,
        cache: Arc<cache::Cache>,
        settings: &Settings,
    ) -> Result<Self> {
        let FigureWriter {
            name,
            caption,
            size,
            plot_count,
            ..
        } = figure;
        if !settings.rebuild && cache.check(&name)? {
            log::info!("[{name}]: Matches cached entry");
            let child = Command::new("/bin/true").spawn()?;
            Ok(Self {
                name,
                caption,
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
                caption,
                child,
                plot_count,
                size,
            })
        }
    }

    pub fn wait(mut self, pb: &ProgressBar, settings: &Settings) -> Result<FinishedFigure> {
        pb.set_length(self.plot_count + 1);
        let mut path = PathBuf::from(&settings.output_dir).join(&self.name);
        path.set_extension("prg");
        loop {
            pb.tick();
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
            caption: self.caption,
            size: self.size,
        })
    }
}
