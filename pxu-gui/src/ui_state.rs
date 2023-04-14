use pxu::UCutType;

#[derive(serde::Deserialize, serde::Serialize)]
#[serde(default)]
pub struct UiState {
    #[serde(skip)]
    pub fullscreen_component: Option<pxu::Component>,
    pub u_cut_type: UCutType,
    pub show_cuts: bool,
    #[serde(skip)]
    pub show_side_panel: bool,
}

impl Default for UiState {
    fn default() -> Self {
        Self {
            fullscreen_component: Default::default(),
            u_cut_type: Default::default(),
            show_cuts: true,
            show_side_panel: true,
        }
    }
}

impl UiState {
    pub fn toggle_fullscreen(&mut self, component: pxu::Component) {
        if self.fullscreen_component.is_some() {
            if self.fullscreen_component != Some(component) {
                log::warn!(
                    "Toggling the wrong fullscreen component ({:?} vs {component:?})",
                    self.fullscreen_component
                );
            }
            self.fullscreen_component = None;
        } else {
            self.fullscreen_component = Some(component);
        }
    }

    pub fn close_fullscreen(&mut self) {
        self.fullscreen_component = None;
    }
}
