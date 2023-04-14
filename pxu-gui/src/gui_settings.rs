#[derive(Default, Debug, serde::Serialize, serde::Deserialize)]
pub struct GuiSettings {
    pub show_fps: bool,
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
