#[derive(Default, Debug, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct Arguments {
    pub show_fps: bool,
    pub show_dev: bool,
    pub continuous_mode: bool,
    pub paths: Option<String>,
}

#[cfg(target_arch = "wasm32")]
impl From<url::Url> for Arguments {
    fn from(url: url::Url) -> Self {
        let Some(query) = url.query() else { return Default::default(); };
        let Ok(settings) = serde_urlencoded::from_str(query) else { return Default::default(); };
        settings
    }
}

#[cfg(target_arch = "wasm32")]
impl From<Option<url::Url>> for Arguments {
    fn from(url: Option<url::Url>) -> Self {
        let Some(url) = url else { return Default::default();};
        Self::from(url)
    }
}

#[cfg(not(target_arch = "wasm32"))]
impl Arguments {
    pub fn parse() -> Self {
        let matches = clap::command!() // requires `cargo` feature
            .arg(
                clap::Arg::new("fps")
                    .short('f')
                    .long("show-fps")
                    .help("Show fps")
                    .action(clap::ArgAction::SetTrue)
                    .required(false),
            )
            .arg(
                clap::Arg::new("dev")
                    .short('d')
                    .long("show-dev")
                    .help("Show dev gui")
                    .action(clap::ArgAction::SetTrue)
                    .required(false),
            )
            .arg(
                clap::Arg::new("continuous_mode")
                    .short('c')
                    .long("continuous")
                    .help("Draw in continuous mode")
                    .action(clap::ArgAction::SetTrue)
                    .required(false),
            )
            .arg(
                clap::Arg::new("paths")
                    .long("paths")
                    .help("Load paths")
                    .required(false),
            )
            .get_matches();

        Self {
            show_fps: matches.get_flag("fps"),
            show_dev: matches.get_flag("dev"),
            continuous_mode: matches.get_flag("continuous_mode"),
            paths: matches.get_one::<String>("paths").map(|s| s.clone()),
        }
    }
}
