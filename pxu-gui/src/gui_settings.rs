#[cfg(not(target_arch = "wasm32"))]
use clap::Parser;

#[derive(Default, Debug, serde::Serialize, serde::Deserialize)]
#[serde(default)]
pub struct GuiSettings {
    pub show_fps: bool,
    pub show_dev: bool,
}

#[cfg(target_arch = "wasm32")]
impl From<url::Url> for GuiSettings {
    fn from(url: url::Url) -> Self {
        let Some(query) = url.query() else { return Default::default(); };
        let Ok(settings) = serde_urlencoded::from_str(query) else { return Default::default(); };
        settings
    }
}

#[cfg(target_arch = "wasm32")]
impl From<Option<url::Url>> for GuiSettings {
    fn from(url: Option<url::Url>) -> Self {
        let Some(url) = url else { return Default::default();};
        Self::from(url)
    }
}

#[cfg(not(target_arch = "wasm32"))]
#[derive(Default, Debug, serde::Serialize, serde::Deserialize, Parser)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    #[arg(short = 'f', long = "fps")]
    show_fps: bool,
    #[arg(short = 'd', long = "dev")]
    show_dev: bool,
}

#[cfg(not(target_arch = "wasm32"))]
impl From<Cli> for GuiSettings {
    fn from(cli: Cli) -> Self {
        Self {
            show_fps: cli.show_fps,
            show_dev: cli.show_dev,
        }
    }
}
