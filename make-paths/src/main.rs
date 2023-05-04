use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use pxu::kinematics::CouplingConstants;

#[derive(Parser, Clone)]
#[command(author, version, about, long_about = None)]
struct Settings {
    #[arg(short, long)]
    compressed: bool,
    #[arg(short, long, action = clap::ArgAction::Count)]
    verbose: u8,
    path_number: Option<usize>,
}

fn main() -> std::io::Result<()> {
    let settings = Settings::parse();

    let consts = CouplingConstants::new(2.0, 5);
    let mut contours = pxu::Contours::new();

    let spinner_style = ProgressStyle::with_template(
        "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
    )
    .unwrap()
    .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ");

    eprintln!("[1/3] Generating contours");
    let pb = ProgressBar::new(1);
    pb.set_style(spinner_style);
    loop {
        pb.set_length(contours.progress().1 as u64);
        pb.set_position(contours.progress().0 as u64);
        if contours.update(0, consts) {
            pb.finish_and_clear();
            break;
        }
    }

    eprintln!("[2/3] Generating paths");
    let saved_paths = make_paths::get_interactive_paths(&contours, consts);

    eprintln!("[3/3] Saving paths");

    let result = if settings.compressed {
        pxu::path::SavedPath::save_compressed(&saved_paths)
    } else {
        pxu::path::SavedPath::save(&saved_paths)
    }
    .unwrap();
    println!("{result}");

    Ok(())
}
