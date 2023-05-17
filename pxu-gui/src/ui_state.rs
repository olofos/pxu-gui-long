use pxu::UCutType;

use crate::arguments::Arguments;

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
    pub show_fps: bool,
    #[serde(skip)]
    pub show_dev: bool,
    #[serde(skip)]
    pub continuous_mode: bool,
    #[serde(skip)]
    pub path_index: Option<usize>,
    #[serde(skip)]
    pub saved_paths_to_load: Option<Vec<pxu::path::SavedPath>>,
    #[serde(skip)]
    pub path_load_progress: Option<(usize, usize)>,
    #[serde(skip)]
    pub show_x: bool,
}

impl Default for UiState {
    fn default() -> Self {
        Self {
            fullscreen_component: Default::default(),
            u_cut_type: Default::default(),
            show_cuts: true,
            show_side_panel: true,
            active_point: 0,
            show_fps: false,
            show_dev: false,
            continuous_mode: false,
            path_index: None,
            saved_paths_to_load: None,
            path_load_progress: None,
            show_x: false,
        }
    }
}

impl UiState {
    pub fn set(&mut self, arguments: Arguments) {
        self.show_fps = arguments.show_fps;
        self.show_dev = arguments.show_dev;
        self.continuous_mode = arguments.continuous_mode;
        self.show_x = arguments.show_x;

        if let Some(ref paths) = arguments.paths {
            let mut saved_paths_to_load = pxu::path::SavedPath::load(paths);
            if let Some(ref mut paths) = saved_paths_to_load {
                self.path_load_progress = Some((0, paths.len()));
                paths.reverse();
            }
            self.saved_paths_to_load = saved_paths_to_load
        }
    }

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
