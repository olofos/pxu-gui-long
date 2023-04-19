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
    pub active_point: usize,
    #[serde(skip)]
    pub edit_path: bool,
    #[serde(skip)]
    pub edit_path_component: pxu::Component,
}

impl Default for UiState {
    fn default() -> Self {
        Self {
            fullscreen_component: Default::default(),
            u_cut_type: Default::default(),
            show_cuts: true,
            show_side_panel: true,
            active_point: 0,
            edit_path: false,
            edit_path_component: pxu::Component::P,
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

    pub fn menu(&mut self, ui: &mut egui::Ui, component: Option<pxu::Component>) {
        ui.menu_button("Cut type", |ui| {
            if UCutType::all()
                .map(|typ| ui.radio_value(&mut self.u_cut_type, typ, typ.to_string()))
                .any(|r| r.clicked())
            {
                ui.close_menu();
            }
        });

        if ui
            .button(if self.show_side_panel {
                "Hide side panel"
            } else {
                "Show side panel"
            })
            .clicked()
        {
            self.show_side_panel = !self.show_side_panel;
            ui.close_menu();
        }

        if let Some(component) = component {
            if ui
                .button(if self.fullscreen_component.is_none() {
                    "Fullscreen"
                } else {
                    "Close fullscreen"
                })
                .clicked()
            {
                self.toggle_fullscreen(component);
                ui.close_menu();
            }
        } else if ui
            .add_enabled(
                self.fullscreen_component.is_some(),
                egui::Button::new("Close fullscreen"),
            )
            .clicked()
        {
            self.close_fullscreen();
            ui.close_menu();
        }
    }
}
