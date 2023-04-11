use egui::{vec2, Pos2};
use num::complex::Complex64;
use pxu::kinematics::CouplingConstants;
use pxu::UCutType;

use crate::plot::Plot;

/// We derive Deserialize/Serialize so we can persist app state on shutdown.
#[derive(serde::Deserialize, serde::Serialize)]
#[serde(default)] // if we add new fields, give them default values when deserializing old state

pub struct PxuGuiApp {
    pxu: pxu::State,
    show_cuts: bool,
    u_cut_type: UCutType,
    #[serde(skip)]
    contours: pxu::Contours,
    p_plot: Plot,
    xp_plot: Plot,
    xm_plot: Plot,
    u_plot: Plot,
    #[serde(skip)]
    #[cfg(debug_assertions)]
    frame_history: crate::frame_history::FrameHistory,
}

impl Default for PxuGuiApp {
    fn default() -> Self {
        let bound_state_number = 1;

        let consts = CouplingConstants::new(2.0, 5);
        let state = pxu::State::new(bound_state_number, consts);

        Self {
            pxu: state,
            contours: pxu::Contours::new(),
            p_plot: Plot {
                component: pxu::Component::P,
                height: 0.75,
                width_factor: 1.5,
                origin: Pos2::new(1.5, 0.0),
            },
            xp_plot: Plot {
                component: pxu::Component::Xp,
                height: (8.0 * consts.s()) as f32,
                width_factor: 1.0,
                origin: Pos2::ZERO,
            },
            xm_plot: Plot {
                component: pxu::Component::Xm,
                height: (8.0 * consts.s()) as f32,
                width_factor: 1.0,
                origin: Pos2::ZERO,
            },
            u_plot: Plot {
                component: pxu::Component::U,
                height: ((4 * consts.k() + 1) as f64 / consts.h) as f32,
                width_factor: 1.0,
                origin: Pos2::ZERO,
            },
            show_cuts: true,
            u_cut_type: Default::default(),
            #[cfg(debug_assertions)]
            frame_history: Default::default(),
        }
    }
}

impl PxuGuiApp {
    /// Called once before the first frame.
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        // This is also where you can customize the look and feel of egui using
        // `cc.egui_ctx.set_visuals` and `cc.egui_ctx.set_fonts`.

        // Load previous app state (if any).
        // Note that you must enable the `persistence` feature for this to work.
        if let Some(storage) = cc.storage {
            return eframe::get_value(storage, eframe::APP_KEY).unwrap_or_default();
        }

        Default::default()
    }
}

impl eframe::App for PxuGuiApp {
    /// Called by the frame work to save state before shutdown.
    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        eframe::set_value(storage, eframe::APP_KEY, self);
    }

    /// Called each time the UI needs repainting, which may be many times per second.
    /// Put your widgets into a `SidePanel`, `TopPanel`, `CentralPanel`, `Window` or `Area`.
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        #[cfg(debug_assertions)]
        {
            self.frame_history
                .on_new_frame(ctx.input().time, _frame.info().cpu_usage);
        }

        #[cfg(not(target_arch = "wasm32"))] // no File->Quit on web pages!
        egui::TopBottomPanel::top("top_panel").show(ctx, |ui| {
            // The top panel is often a good place for a menu bar:
            egui::menu::bar(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Quit").clicked() {
                        _frame.close();
                    }
                });
            });
        });

        let old_consts = self.pxu.consts;
        let mut new_consts = self.pxu.consts;

        egui::SidePanel::right("side_panel").show(ctx, |ui| {
            ui.heading("Side Panel");

            ui.add(
                egui::Slider::new(&mut new_consts.h, 0.1..=10.0)
                    .text("h")
                    .logarithmic(true),
            );

            ui.add(
                egui::Slider::from_get_set(0.0..=10.0, |v| new_consts.get_set_k(v))
                    .integer()
                    .text("k"),
            );
            ui.add(
                egui::Slider::from_get_set(1.0..=8.0, |n| {
                    if let Some(n) = n {
                        self.pxu = pxu::State::new(n as usize, self.pxu.consts);
                    }
                    self.pxu.points.len() as f64
                })
                .integer()
                .text("M"),
            );

            #[cfg(debug_assertions)]
            ui.add(egui::Checkbox::new(&mut self.show_cuts, "Show cuts"));

            ui.horizontal(|ui| {
                ui.label("U cuts: ");
                ui.radio_value(&mut self.u_cut_type, UCutType::Long, "Long");
                ui.radio_value(&mut self.u_cut_type, UCutType::SemiShort, "Semi short");
                ui.radio_value(&mut self.u_cut_type, UCutType::Short, "Short");
            });

            if ui.add(egui::Button::new("Reset")).clicked() {
                *self = Self::default();
            }

            if ui.add(egui::Button::new("Test")).clicked() {
                let s = serde_json::to_string(&self).unwrap();
                log::info!("json: {s}");
            }

            ui.separator();

            {
                let p = self.pxu.points.iter().map(|pxu| pxu.p).sum::<Complex64>();
                ui.label(format!("Momentum: {:.2}", p));

                let en = self
                    .pxu
                    .points
                    .iter()
                    .map(|pt| {
                        let xp = pt.xp;
                        let xm = pt.xm;
                        -Complex64::i() * self.pxu.consts.h / 2.0 * (xp - 1.0 / xp - xm + 1.0 / xm)
                    })
                    .sum::<Complex64>();
                ui.label(format!(
                    "Energy: {:.2} / {:.2}",
                    en,
                    ::pxu::kinematics::en(p, self.pxu.points.len() as f64, self.pxu.consts)
                ));
            }

            ui.separator();

            {
                ui.label("Active excitation:");

                ui.label(format!("Momentum: {:.2}", self.pxu.active_point().p));

                let xp = self.pxu.active_point().xp;
                let xm = self.pxu.active_point().xm;

                ui.label(format!(
                    "Energy: {:.2}",
                    -Complex64::i() * self.pxu.consts.h / 2.0 * (xp - 1.0 / xp - xm + 1.0 / xm)
                ));

                ui.add_space(10.0);
                ui.label(format!("x+: {:.3}", self.pxu.active_point().xp));
                ui.label(format!("x-: {:.3}", self.pxu.active_point().xm));
                ui.label(format!("u: {:.3}", self.pxu.active_point().u));

                ui.add_space(10.0);
                ui.label("Branch info:");

                ui.label(format!(
                    "Log branches: {:+} {:+}",
                    self.pxu.active_point().sheet_data.log_branch_p,
                    self.pxu.active_point().sheet_data.log_branch_m
                ));

                ui.label(format!(
                    "E branch: {:+} ",
                    self.pxu.active_point().sheet_data.e_branch
                ));
                ui.label(format!(
                    "U branch: ({:+},{:+}) ",
                    self.pxu.active_point().sheet_data.u_branch.0,
                    self.pxu.active_point().sheet_data.u_branch.1
                ));
            }

            #[cfg(debug_assertions)]
            {
                ui.separator();
                {
                    ui.label(format!("FPS: {}", self.frame_history.fps()));

                    self.frame_history.ui(ui);

                    let (current, total) = self.contours.progress();
                    ui.add(egui::ProgressBar::new(current as f32 / total as f32).show_percentage());
                }
            }

            ui.with_layout(egui::Layout::bottom_up(egui::Align::LEFT), |ui| {
                ui.horizontal(|ui| {
                    ui.spacing_mut().item_spacing.x = 0.0;
                    ui.label("powered by ");
                    ui.hyperlink_to("egui", "https://github.com/emilk/egui");
                    ui.label(" and ");
                    ui.hyperlink_to(
                        "eframe",
                        "https://github.com/emilk/egui/tree/master/crates/eframe",
                    );
                    ui.label(".");
                });
            });
        });

        if old_consts != new_consts {
            self.pxu = pxu::State::new(self.pxu.points.len(), new_consts);
            self.contours.clear();
        }

        {
            let start = chrono::Utc::now();
            while (chrono::Utc::now() - start).num_milliseconds()
                < (1000.0 / 20.0f64).floor() as i64
            {
                if self
                    .contours
                    .update(self.pxu.active_point().p.re.floor() as i32, self.pxu.consts)
                {
                    break;
                }
                ctx.request_repaint();
            }
        }

        egui::CentralPanel::default().show(ctx, |ui| {
            let available_size = ui.available_size();
            ui.horizontal(|ui| {
                self.p_plot.draw(
                    ui,
                    available_size * vec2(0.49, 0.49),
                    &mut self.contours,
                    self.show_cuts,
                    self.u_cut_type,
                    &mut self.pxu,
                );

                self.u_plot.draw(
                    ui,
                    available_size * vec2(0.49, 0.49),
                    &mut self.contours,
                    self.show_cuts,
                    self.u_cut_type,
                    &mut self.pxu,
                );
            });
            ui.horizontal(|ui| {
                self.xp_plot.draw(
                    ui,
                    available_size * vec2(0.49, 0.49),
                    &mut self.contours,
                    self.show_cuts,
                    self.u_cut_type,
                    &mut self.pxu,
                );

                self.xm_plot.draw(
                    ui,
                    available_size * vec2(0.49, 0.49),
                    &mut self.contours,
                    self.show_cuts,
                    self.u_cut_type,
                    &mut self.pxu,
                );
            });
        });
    }
}
