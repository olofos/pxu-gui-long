use pxu::{kinematics::CouplingConstants, Pxu};
use std::io::Result;
use std::sync::Arc;

use clap::Parser;
use indicatif::{MultiProgress, ProgressBar, ProgressStyle};

mod cache;
mod fig_compiler;
mod fig_writer;
mod figures;
mod utils;

use crate::figures::ALL_FIGURES;
use crate::utils::{error, Settings, Summary, SUMMARY_NAME};

fn main() -> std::io::Result<()> {
    let settings = Settings::parse();

    if settings.verbose > 0 {
        tracing_subscriber::fmt()
            .with_max_level(tracing::Level::INFO)
            .with_file(true)
            .with_line_number(true)
            .with_writer(std::io::stderr)
            .without_time()
            .init();
        log::set_max_level(log::LevelFilter::Debug);
    }

    let num_threads = if let Some(jobs) = settings.jobs {
        jobs
    } else {
        num_cpus::get()
    }
    .min(ALL_FIGURES.len());

    if settings.rebuild {
        println!(" ---  Rebuilding all figures");
    }

    let spinner_style = ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
    )
    .unwrap();
    let spinner_style_no_progress =
        ProgressStyle::with_template("[{elapsed_precise}] {spinner} {msg}")
            .unwrap()
            .tick_chars("⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏");

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

    let mut pxu = Pxu::new(consts);
    pxu.contours = contours;
    pxu.state = pxu::State::new(1, consts);

    let mb = Arc::new(MultiProgress::new());
    let pb = if settings.verbose == 0 {
        println!("[2/5] Loading paths");
        mb.add(ProgressBar::new(1))
    } else {
        ProgressBar::hidden()
    };

    let saved_paths = make_paths::get_plot_paths(&pxu.contours, pxu.consts);
    pb.set_style(spinner_style.clone());
    pb.set_length(saved_paths.len() as u64);

    let pool = threadpool::ThreadPool::new(num_threads);
    let (tx, rx) = std::sync::mpsc::channel();

    let saved_paths_len = saved_paths.len();

    for saved_path in saved_paths {
        let tx = tx.clone();
        let contours = pxu.contours.clone();
        let consts = pxu.consts;
        let spinner_style = spinner_style_no_progress.clone();
        let settings = settings.clone();
        let mb = mb.clone();
        pool.execute(move || {
            let pb = if settings.verbose == 0 {
                mb.add(ProgressBar::new_spinner())
            } else {
                ProgressBar::hidden()
            };
            pb.set_style(spinner_style);
            pb.enable_steady_tick(std::time::Duration::from_millis(100));

            pb.set_message(saved_path.name.clone());
            let path = pxu::path::Path::from_base_path(saved_path.into(), &contours, consts);
            tx.send(path).unwrap();
            pb.finish_and_clear();
        });
    }

    pxu.paths = rx
        .into_iter()
        .take(saved_paths_len)
        .map(|r| {
            pb.inc(1);
            r
        })
        .collect::<Vec<_>>();

    pool.join();
    pb.finish_and_clear();

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

    let pool = threadpool::ThreadPool::new(num_threads);
    let (tx, rx) = std::sync::mpsc::channel();

    let pb = if settings.verbose == 0 {
        mb.add(ProgressBar::new_spinner())
    } else {
        ProgressBar::hidden()
    };

    pb.set_style(spinner_style.clone());
    pb.set_message("Building figures");
    pb.set_length(ALL_FIGURES.len() as u64);
    pb.enable_steady_tick(std::time::Duration::from_millis(250));

    for (i, f) in ALL_FIGURES.iter().enumerate() {
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

    let mut finished_figures = rx
        .into_iter()
        .take(ALL_FIGURES.len())
        .map(|r| {
            pb.inc(1);
            r
        })
        .collect::<Result<Vec<_>>>()?;
    pool.join();
    pb.finish_and_clear();

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

    let pb = if settings.verbose == 0 {
        ProgressBar::new_spinner()
    } else {
        ProgressBar::hidden()
    };

    pb.set_style(spinner_style_no_progress);
    pb.enable_steady_tick(std::time::Duration::from_millis(100));

    if summary.finish(&settings, &pb)?.wait()?.success() {
        log::info!("[{SUMMARY_NAME}] Done.");
    } else {
        log::error!("[{SUMMARY_NAME}] Error.");
        return Err(error("Error compiling summary"));
    }

    pb.finish_and_clear();

    Ok(())
}
