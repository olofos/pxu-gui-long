use egui::{vec2, Pos2};
use itertools::Itertools;
use pxu::kinematics::CouplingConstants;
use pxu::UCutType;

use crate::plot::Plot;

#[derive(Debug)]
struct AnimPathSegment {
    t_start: f64,
    t_end: f64,
    p_start: num::complex::Complex64,
    p_end: num::complex::Complex64,
}

struct AnimData {
    playing: bool,
    segment: usize,
    prev_frame: chrono::DateTime<chrono::Utc>,
    path: Vec<AnimPathSegment>,
    p: num::complex::Complex64,
    active_point: usize,
    speed: f64,
    t: f64,
}

impl Default for AnimData {
    fn default() -> Self {
        Self {
            playing: Default::default(),
            segment: Default::default(),
            prev_frame: Default::default(),
            path: Default::default(),
            p: Default::default(),
            active_point: Default::default(),
            speed: 10.0,
            t: Default::default(),
        }
    }
}

impl AnimData {
    fn next(&mut self) -> Option<(usize, num::complex::Complex64)> {
        if !self.playing {
            return None;
        }
        let t1 = chrono::Utc::now();
        self.t += (t1 - self.prev_frame).num_nanoseconds().unwrap() as f64 / 1_000_000_000.0
            * self.speed
            / 100.0;
        self.prev_frame = t1;

        while self.path[self.segment].t_end < self.t {
            self.segment += 1;
            if self.segment >= self.path.len() {
                self.playing = false;
                return Some((self.active_point, self.path.last().unwrap().p_end));
            }
        }

        let seg = &self.path[self.segment];

        let dt = self.t - seg.t_start;
        let dp = (seg.p_end - seg.p_start) / (seg.p_end - seg.p_start).norm();
        Some((self.active_point, seg.p_start + dt * dp))
    }

    pub fn start(&mut self, paths: &[Vec<pxu::Point>]) -> Vec<pxu::Point> {
        self.active_point = paths
            .iter()
            .map(|path| {
                path.iter()
                    .map(|pt| pt.p)
                    .tuple_windows::<(_, _)>()
                    .map(|(p1, p2)| (p2 - p1).norm())
                    .sum::<f64>()
            })
            .enumerate()
            .max_by(|(_, x), (_, y)| x.partial_cmp(y).unwrap())
            .map(|(n, _)| n)
            .unwrap();

        let p = paths[self.active_point][0].p;
        self.p = p;
        self.path = vec![];

        let ps = paths[self.active_point]
            .iter()
            .map(|pt| pt.p)
            .collect::<Vec<_>>();
        let mut t = 0.0;
        for i in 0..ps.len() - 1 {
            let p_start = ps[i];
            let p_end = ps[i + 1];
            let dp = p_end - p_start;
            let dt = dp.norm();

            let t_start = t;
            let t_end = t + dt;

            self.path.push(AnimPathSegment {
                p_start,
                p_end,
                t_start,
                t_end,
            });

            t += dt;
        }

        self.segment = 0;
        self.playing = true;
        self.prev_frame = chrono::Utc::now();
        self.t = 0.0;

        if self.speed == 0.0 {
            self.speed = 10.0;
        }

        paths.iter().map(|path| path[0].clone()).collect()
    }
}

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
    #[serde(skip)]
    anim_data: AnimData,
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
                origin: Pos2::new(0.5, 0.0),
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
            anim_data: AnimData::default(),
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
                        self.anim_data.playing = false;
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

            {
                let enabled = !self.anim_data.playing
                    && !self.pxu.paths.is_empty()
                    && self.pxu.paths[0].len() > 1;
                if ui.add_enabled(enabled, egui::Button::new("Play")).clicked() {
                    self.pxu.points = self.anim_data.start(&self.pxu.paths);
                }

                if ui
                    .add_enabled(self.anim_data.playing, egui::Button::new("Stop"))
                    .clicked()
                {
                    self.anim_data.playing = false;
                }

                ui.add(egui::Slider::new(&mut self.anim_data.speed, 1.0..=100.0).text("Speed"));
            }

            ui.separator();

            {
                ui.label(format!("Momentum: {:.3}", self.pxu.p()));
                ui.label(format!("Energy: {:.3}", self.pxu.en()));
            }

            ui.separator();

            {
                ui.label("Active excitation:");

                ui.label(format!("Momentum: {:.3}", self.pxu.active_point().p));

                ui.label(format!(
                    "Energy: {:.3}",
                    self.pxu.active_point().en(self.pxu.consts)
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

                {
                    let xp = self.pxu.active_point().xp;
                    let xm = xp.conj();
                    let h = self.pxu.consts.h;
                    let k = self.pxu.consts.k() as f64;
                    let p = xp.arg() / std::f64::consts::PI;
                    let m = h / 2.0
                        * (xp + 1.0 / xp
                            - xm
                            - 1.0 / xm
                            - 2.0 * num::complex::Complex64::i() * (k * p) / h)
                            .im;
                    ui.label(format!("p = {p:.3} m = {m:.3}"));
                }
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
            self.anim_data.playing = false;
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

        if self.anim_data.playing {
            if let Some((active_point, p)) = self.anim_data.next() {
                self.pxu.update(
                    active_point,
                    pxu::Component::P,
                    p,
                    &self.contours,
                    self.u_cut_type,
                );
                ctx.request_repaint();
            } else {
                self.anim_data.playing = false;
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

        if ctx.input().key_pressed(egui::Key::Space) {
            for (i, pt) in self.pxu.points.iter().enumerate() {
                self.pxu.paths[i].push(pt.clone());
            }
        }
    }
}
