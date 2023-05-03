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

const PATHS: &str = r"7VxLayTJEf4rYk8SVLTz_dDR4INP3suCYRkGjaaZER61hNTC8_P9RT4qI7WSGQZjWEhdOrMrKzIeX3xRWZmt3y-fzzdP5-vLx4e70_n5-vfLx-tLUgeVrbZK6RhdtNFZs_GXynqtrTM6ZG21dslfbd_5Bn3wiv9ccNpEm6Pb1MHl8ueUMdFpHy0G35fBLluTrQ1BZa2cSiwcw0IZG1SKyrLoFww2B8j0iYeZoLXezMH7gF420EMHa9PV9vz1eDx__Hxzvrm-_Pbw5eOnp5vT7dePj9ekN9G_v1bbsXWu9fbSm5d_PZ7_fTyetvZ5td3df_z-8fnuy4nV3fTV1dWHq-3z8dv5hl2knd_Ip4SZL5213I4O7Wg0t11AWysVuaMsd0xW6MCX3PGRh4VYOsnBeu8zDzPa8DCXPXdsYgHFEZcmlGHwAHey4jm14StWe0hzOfKk1mbuJF2uBOPQCY5F2-QgzbmgWGcVeJgJfA9_oKOKaB6wEQ_mDm5FB2K4A6HoYALuYDp0MDU6rMhGrBR3DItmdbkD5dGBIdyBWejARO7A4I3YeO7AFejALdyJrCg7jDuJLWVXlo7DpPBxaUe2zZYpEQgMcr63exND-gjc2W-ExC6Qp-ozsQ5dBVau68Zad6XZnG4N29nNZAc0-4tnmmOKy7rH2Jfdlezk7mP2fnc-h6VFpcSrhasEsseRI9wDzKFvkS-YaJAoYOlYYRR1EDG8GroK7hrsCiA7HhmpHagM4YZghnZDNiO-AZ4ToeUBybYYQuJOkhJJTkVSB5LKkdSapDkk7STpAJKeIekykr4k6WSS3icZFpLxIhlIkhEmGXqSmCAJFpIoIgkvkrgjCUiSSCUJYRLQJoF4EolAMkFIZg7JlKIp2WjKQ5pSlKbspSmxacp5muiAJqagiURo4heaqIcmVqKJsGjiMppojiYGpIkcaeJNmiiVJraliYhp4mia6JsmZqeJ9GmqByQrBckSQrK2rDqz6syqM6vOrDqz6syqM6vOrDqz6syqM6vOrDqz6syqM6vOrDqz6syqM6_rzIft9uH-8eF0PJ2v__m4Hb_f3p1vzncPJ-z-3D6cnrHLdfn12hzU9q9rf4DBp5v74_Uv3x8vbu-ebr8dLz7VXaC_tM9fIP6NTTJ10NjjUnUDzByCy9iLiiljgwy6BFxLvL2VjIJWEXtauW5_vTGU3hr7wkN19MFrnyyGWWN9nfHd7S71Q7td_3g5P999Pm7tc9rt-uNmF-IBStAVF-CDzC3gylWkJqROqhHAPp6L5TLsQrPFlnmjshhyAVt5pZkzCKSMNQips-Vb7OiBSWKFkt1cRZxJaKqqgDJbyxCL22w0FaBps76i0zOvFGG4uNmaoQxE3FphD0ZJNY-QOJimpkbcMHkFVdqMafyVNyhamh6sEBuRBfBLza8IQrF1bMRYVUkcLkmN5uLWSDJlkF_NXrjY9Ba1Ji5Tu85J029CFjZJEE9dfmQSadNySndtPLNVU9Lx-KY7LKJuEgylbikTR3eAYrpJtnuLmrvYidS96AvRpO5y6j7nnOyhUDy-RojjRi1wHE7q8USUqYeZ723RZ0xQAwVDhTpWfCGj0IFFHVm8b9wAx55qMGTKa9i0zHcVsJrnLCgerXF53DMECfFiVqGM0FGoLiwShgr7hVuEt4QThW-Fy0UkRIBE3EQ4RZRF8AUmBFQEggSwBN4GCgc0B14HiHdg71jf4b9nxMiSkTojn0aSjcwb6ThydCTuyOaR4iPvBxkMhhi0MbhkEMxgnUFFg58GaQ0mG_Q2OG8Q4U6OO1_uFNpYVRarX3-4Vv3aS9XN08PL6fPFw9Pdl7vTxfkrul--Xvzt4hbU_m7dUn4ULocSk70yzvnonTGJr4Wgc0C99xH1SDtdC9cbQ-mtsShc9oBHgey0Sco7b1L29fjI9Ke1VsYgDC4Y1Ljg_9-FLQBP7Wmkt6xCq7Kp6RUBETa2Vj1-AC3Mgvib-iiYGJ_lKtBh6nMXDMPDaeMdIFS7zl2mPiTCT2jqym7QrD0GASU6VZ4LSIBY6S8iK0IVlpAqtTZofuStT4AQtBXHczP2goRJkY6tCWFaVb51W25lN2910QFDtlo1YNwW2nohb7Wq1ESoLYcUqy0PRXsr9AaNVr_qmKgE6XfZ1IRbZhTftaCmBlSjphs_iDaFNVNSMwPGUbeOnw-70UyHzRe6PGtXF7HjqHsuMvc1hwaW2fzMRaC7n59xe1QcU2wLli3t2CNLPbQIOLWIcyloOEiFe3MDDDXEhPKc35Yx3MwNb70wWxZWUahZQoVm6K29MS6OO4aYIXtMOLQYqg19hRXCOGGzcIXwkHCc8Kdws_C-CIqIlQihiKwI-IDBwMYAzEDRgNbA2wDhQOaO1o7fHdE7yHfc76mwZ8eeMHsO7Wk1Um3k30jKkakjfUdOj0Qf2T8oYfDEII_BKINmBvfsfLRT1M5aO5Ht3LbT3c6AOylWmvy5QrWvqVqFOj2c5yr1TpGyB5VUsKi2AYY7lcsiyasUrNZJoboan91--BB1xClvssVCMigPf-P-ZLCg1Nbw6olXq-XsIaqQxuI1ahWxitTwsrFcDT0ERos7rIdkHv5SJMdyhFErbwPKlGEcH2I966hxT-TS_27Vsu9WLXrr9OHfT29ULXpdtmzgybeCGGO8wSOWqQiwQFxS5YJOHP2cCkpQnPH4oUvkAQ88nvlGIPAZ7reKAZaTB6xtgVPmQfAMQ5A_Nv6CAYTLGw9kWOC2rUgor_U4gXQdxVNtZdbygi0BuKxPwRQ0L5qWZyKzFQvqU7voTMNokkCTcJrmpUklksqStIKkeSTtpskjNDmLJj_S5GKavE8yLiTiNV94dc8s79Vckx6TjpP6k2GTyZMzZj_NPpz9O_t-josMmWxPoyYBk-xp2kkjqas0QlonzZ4cMgN8hr5MiildRFhWeq30Wum10mul10qvlV4rvVZ6rfRa6bXSa6XXSq-VXiu9_rTpJV_g_vbDL3B_u_h0gw3GfiTmrryifPelLfY_seeHl7F4OYsNPT4Eg93CrLAHqGJW2eG7-tIWr1wddhEj3rhqaGg1tNQHrU02OuBdK9-STD8xYxV8oqO13jlly1v8g1E2JIN3tnBO2cd9KVLLX_bWRDggYwMSyujy2rb9HJ1_kY5Xv-6n3tu-9dr2rd3GN340buKeLAb7yfxWXdX9Bd7tjlrXk6x4J11wVPcGMt7Lo1Mjzp9b-a7Fn9_5YzzjCrdvRVIRr_htvk1tP072poE0C6F5Bpqnp1k3mvSmySRswQximK-8umsW-WrC1-rMys6WTFZOnWncLOOVR6fZZ9XmIEmLhKkrrCusK6wrrCusK6wrrCusK6xvhfV_sg55qM_c7yxE-JSHCXj2z1ZrLDCU8fUAiQrJZ2-88irH6PpSRIfyH6d4qAt8flQfLO7yEGCCwprFiuMjNlmdIFYhYMalwAcmUw4pp2SSzVjT8Cmdlzq4rkawrlNYlIUQgzN8TFf7Q1uKaJUVH_j_8dUIhZ_6H1b0X9YjOsMfm86mZlD0WGWGelg6Giw5s3HtzJfnQ02unO4tn1v5jhetGLKV0Zy3gc8Cs5wi3PDJJcxQT-XJ3jSQJhk0y6d5cpoUo0lnmsyZEwvxdJtJqSZW0harbl8lAiO8hAdOa65iaVnOXrXpubGVL-tBwMQ_TeJbyvsIyNiKuCo3IZd4Kka-7EzjZhnzBPPss2qT2pNFwtQV0hXSFdIV0hXSFdIV0hXSFdJtPPH-7BLk8fh09_D5j78P_vAf";

fn main() -> std::io::Result<()> {
    let settings = Settings::parse();

    if settings.verbose > 0 {
        tracing_subscriber::fmt::fmt()
            .with_writer(std::io::stderr)
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

    let saved_paths = pxu::path::SavedPath::load(PATHS).unwrap();
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
