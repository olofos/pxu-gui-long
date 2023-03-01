// use std::f64::consts::PI;

const PI: f64 = 0.5;

use egui::style::Margin;
use egui::{vec2, Color32, Pos2, Rect, Stroke, Ui, Vec2};
use num::complex::Complex64;

use crate::kinematics::CouplingConstants;
use crate::pxu;
use crate::pxu::PxuPoint;

/// We derive Deserialize/Serialize so we can persist app state on shutdown.
#[derive(serde::Deserialize, serde::Serialize)]
#[serde(default)] // if we add new fields, give them default values when deserializing old state

pub struct TemplateApp {
    #[serde(skip)]
    consts: CouplingConstants,
    #[serde(skip)]
    pxu: PxuPoint,
    #[serde(skip)]
    z: num::complex::Complex64,
    #[serde(skip)]
    branch: i32,
    #[serde(skip)]
    grid: pxu::Grid,
    #[serde(skip)]
    cuts: Vec<pxu::Cut>,
    #[serde(skip)]
    p_plot: Plot,
    #[serde(skip)]
    xp_plot: Plot,
    #[serde(skip)]
    xm_plot: Plot,
    #[serde(skip)]
    u_plot: Plot,
    #[serde(skip)]
    p_range: i32,
    #[serde(skip)]
    show_dots: bool,
}

struct Plot {
    component: pxu::Component,
    height: f32,
    width_factor: f32,
    origin: Pos2,
}

impl Plot {
    fn draw(
        &mut self,
        ui: &mut Ui,
        desired_size: Vec2,
        grid: &mut pxu::Grid,
        cuts: &Vec<pxu::Cut>,
        show_dots: bool,
        pxu: &mut PxuPoint,
    ) {
        let contours = match self.component {
            pxu::Component::P => &grid.p,
            pxu::Component::Xp | pxu::Component::Xm => &grid.x,
            pxu::Component::U => &grid.u,
        };

        egui::Frame::canvas(ui.style())
            .outer_margin(Margin::same(0.0))
            .inner_margin(Margin::same(0.0))
            .show(ui, |ui| {
                let (response, painter) = ui.allocate_painter(desired_size, egui::Sense::drag());

                if response.hovered() {
                    let zoom = ui.input().zoom_delta();
                    self.zoom(zoom);

                    if response.dragged_by(egui::PointerButton::Secondary) {
                        let delta = response.drag_delta();
                        self.origin -= Vec2::new(
                            delta.x * (self.height / desired_size.y) * (self.width_factor),
                            delta.y * (self.height / desired_size.y),
                        );
                    }
                }

                let rect = response.rect;

                let visible_rect = Rect::from_center_size(
                    self.origin,
                    vec2(
                        self.height * self.width_factor * desired_size.x / desired_size.y,
                        self.height,
                    ),
                );

                let to_screen = eframe::emath::RectTransform::from_to(visible_rect, rect);
                ui.set_clip_rect(rect);

                let origin = to_screen * egui::pos2(0.0, 0.0);

                let mut shapes = vec![
                    egui::epaint::Shape::line(
                        vec![
                            egui::pos2(rect.left(), origin.y),
                            egui::pos2(rect.right(), origin.y),
                        ],
                        Stroke::new(0.75, Color32::GRAY),
                    ),
                    egui::epaint::Shape::line(
                        vec![
                            egui::pos2(origin.x, rect.bottom()),
                            egui::pos2(origin.x, rect.top()),
                        ],
                        Stroke::new(0.75, Color32::GRAY),
                    ),
                ];

                for points in contours {
                    let points = points
                        .iter()
                        .map(|z| to_screen * egui::pos2(z.re as f32, -z.im as f32))
                        .collect::<Vec<_>>();

                    shapes.push(egui::epaint::Shape::line(
                        points.clone(),
                        Stroke::new(0.5, Color32::GRAY),
                    ));

                    if show_dots {
                        for center in points {
                            shapes.push(egui::epaint::Shape::Circle(egui::epaint::CircleShape {
                                center,
                                radius: 1.5,
                                fill: Color32::RED,
                                stroke: Stroke::NONE,
                            }));
                        }
                    }
                }

                let z = match self.component {
                    pxu::Component::P => pxu.p,
                    pxu::Component::U => pxu.u,
                    pxu::Component::Xp => pxu.xp,
                    pxu::Component::Xm => pxu.xm,
                };

                let size = egui::epaint::Vec2::splat(8.0);
                let center = to_screen * egui::pos2(z.re as f32, -z.im as f32);
                let point_rect = egui::Rect::from_center_size(center, size);
                let point_id = response.id.with(0);
                let point_response = ui.interact(point_rect, point_id, egui::Sense::drag());

                let stroke = if point_response.hovered() || point_response.dragged() {
                    egui::epaint::Stroke::new(2.0, Color32::LIGHT_BLUE)
                } else {
                    egui::epaint::Stroke::NONE
                };

                let radius = if point_response.hovered() || point_response.dragged() {
                    6.0
                } else {
                    4.0
                };

                if point_response.dragged() {
                    let new_value = to_screen.inverse() * (center + point_response.drag_delta());
                    let new_value = Complex64::new(new_value.x as f64, -new_value.y as f64);

                    let mut new_log_branch = pxu.log_branch;

                    for cut in cuts {
                        if cut.component == self.component
                            && cut.is_visible(&pxu, pxu.p.re.floor() as i32)
                            && cut.intersection(z, new_value).is_some()
                        {
                            match cut.typ {
                                pxu::CutType::LogX(_, branch) => {
                                    new_log_branch = branch;
                                }
                                _ => {}
                            }
                            log::info!("Intersection with {:?}", cut.typ);
                        }
                    }

                    log::info!("new_log_branch = {new_log_branch}");

                    match self.component {
                        pxu::Component::P => {
                            *pxu = PxuPoint::new(new_value, new_log_branch, pxu.consts);
                        }
                        pxu::Component::Xp => {
                            if let Some(new_pxu) = pxu.shift_xp(new_value, new_log_branch) {
                                *pxu = new_pxu;
                            }
                        }
                        pxu::Component::Xm => {
                            if let Some(new_pxu) = pxu.shift_xm(new_value, new_log_branch) {
                                *pxu = new_pxu;
                            }
                        }
                        pxu::Component::U => {
                            if let Some(new_pxu) = pxu.shift_u(new_value, new_log_branch) {
                                *pxu = new_pxu;
                            }
                        }
                    };
                }

                let center = to_screen * egui::pos2(z.re as f32, -z.im as f32);

                shapes.push(egui::epaint::Shape::Circle(egui::epaint::CircleShape {
                    center,
                    radius,
                    fill: Color32::BLUE,
                    stroke,
                }));

                for cut in cuts.iter().filter(|c| {
                    c.component == self.component && c.is_visible(&pxu, pxu.p.re.floor() as i32)
                }) {
                    for points in cut.paths.iter() {
                        let points = points
                            .iter()
                            .map(|z| to_screen * egui::pos2(z.re as f32, -z.im as f32))
                            .collect::<Vec<_>>();

                        let color = match cut.typ {
                            pxu::CutType::U(comp) => {
                                if comp == pxu::Component::Xp {
                                    Color32::from_rgb(255, 0, 0)
                                } else {
                                    Color32::from_rgb(0, 192, 0)
                                }
                            }
                            pxu::CutType::LogX(comp, _) => {
                                if comp == pxu::Component::Xp {
                                    Color32::from_rgb(255, 128, 128)
                                } else {
                                    Color32::from_rgb(128, 255, 128)
                                }
                            }
                            pxu::CutType::E => Color32::BLACK,
                        };

                        shapes.push(egui::epaint::Shape::line(
                            points.clone(),
                            Stroke::new(3.0, color),
                        ));

                        if show_dots {
                            for center in points {
                                shapes.push(egui::epaint::Shape::Circle(
                                    egui::epaint::CircleShape {
                                        center,
                                        radius: 2.5,
                                        fill: Color32::RED,
                                        stroke: Stroke::NONE,
                                    },
                                ));
                            }
                        }

                        let branch_points = cut
                            .branch_points
                            .iter()
                            .map(|z| to_screen * egui::pos2(z.re as f32, -z.im as f32))
                            .collect::<Vec<_>>();

                        for center in branch_points {
                            shapes.push(egui::epaint::Shape::Circle(egui::epaint::CircleShape {
                                center,
                                radius: 4.0,
                                fill: color,
                                stroke: Stroke::NONE,
                            }));
                        }
                    }
                }

                {
                    let f = ui.fonts();

                    let text = match self.component {
                        pxu::Component::P => "p",
                        pxu::Component::U => "u",
                        pxu::Component::Xp => "X+",
                        pxu::Component::Xm => "X-",
                    };

                    let text_shape = egui::epaint::Shape::text(
                        &f,
                        rect.right_top() + vec2(-10.0, 10.0),
                        egui::Align2::RIGHT_TOP,
                        text,
                        egui::TextStyle::Monospace.resolve(ui.style()),
                        Color32::BLACK,
                    );
                    shapes.push(egui::epaint::Shape::rect_filled(
                        text_shape.visual_bounding_rect().expand(6.0),
                        egui::Rounding::none(),
                        Color32::WHITE,
                    ));
                    shapes.push(egui::epaint::Shape::rect_stroke(
                        text_shape.visual_bounding_rect().expand(4.0),
                        egui::Rounding::none(),
                        egui::Stroke::new(0.5, Color32::BLACK),
                    ));
                    shapes.push(text_shape);
                }

                painter.extend(shapes);
            });
    }

    fn zoom(&mut self, zoom: f32) {
        self.height /= zoom;
    }
}

impl Default for TemplateApp {
    fn default() -> Self {
        let consts = CouplingConstants::new(2.0, 5);
        let p_range = 1;
        Self {
            consts,
            pxu: PxuPoint::new(
                num::complex::Complex::from(p_range as f64 + 0.25),
                p_range,
                consts,
            ),
            z: num::complex::Complex::new(0.0, 0.5),
            branch: 1,
            grid: pxu::Grid::new(p_range, consts),
            cuts: pxu::Cut::get(p_range, consts),
            p_plot: Plot {
                component: pxu::Component::P,
                height: 0.75,
                width_factor: 1.5,
                origin: Pos2::new(((2 * p_range + 1) as f32) * PI as f32, 0.0),
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
                height: ((2 * consts.k() + 1) as f64 / consts.h) as f32,
                width_factor: 1.0,
                origin: Pos2::ZERO,
            },
            p_range,
            show_dots: false,
        }
    }
}

impl TemplateApp {
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

impl eframe::App for TemplateApp {
    /// Called by the frame work to save state before shutdown.
    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        eframe::set_value(storage, eframe::APP_KEY, self);
    }

    /// Called each time the UI needs repainting, which may be many times per second.
    /// Put your widgets into a `SidePanel`, `TopPanel`, `CentralPanel`, `Window` or `Area`.
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
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

        let old_consts = self.consts;
        let old_p_range = self.p_range;

        egui::SidePanel::right("side_panel").show(ctx, |ui| {
            ui.heading("Side Panel");

            ui.add(
                egui::Slider::new(&mut self.consts.h, 0.1..=10.0)
                    .text("h")
                    .logarithmic(true),
            );

            ui.add(
                egui::Slider::from_get_set(0.0..=10.0, |v| self.consts.get_set_k(v))
                    .integer()
                    .text("k"),
            );

            ui.add(
                egui::Slider::from_get_set(1.00001..=5.0, |v| self.consts.get_set_s(v)).text("s"),
            );

            ui.add(egui::Slider::new(&mut self.p_range, -10..=5).text("Range"));

            ui.add(egui::Checkbox::new(&mut self.show_dots, "Show dots"));

            if ui.add(egui::Button::new("Reset")).clicked() {
                *self = Self::default();
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

        if old_consts != self.consts || old_p_range != self.p_range {
            // self.pxu = PxuPoint::new(
            //     self.pxu.p + 2.0 * PI * (self.p_range - old_p_range) as f64,
            //     self.consts,
            // );

            self.grid = pxu::Grid::new(self.p_range, self.consts);
            self.cuts = pxu::Cut::get(self.p_range, self.consts);
        }

        if old_consts != self.consts {
            self.pxu = PxuPoint::new(self.pxu.p, self.pxu.log_branch, self.consts);

            // self.xp_plot.height *= (self.consts.s() / old_consts.s()) as f32;

            // self.u_plot.height /= (self.consts.h / old_consts.h) as f32;
            // if self.consts.k() > 1 && old_consts.k() > 1 {
            //     self.xp_plot.height *= (2 * self.consts.k()) as f32 / (2 * old_consts.k()) as f32;
            //     self.u_plot.height *= (2 * self.consts.k()) as f32 / (2 * old_consts.k()) as f32;
            // }
        }

        if old_p_range != self.p_range {
            self.p_plot.origin.x += (self.p_range - old_p_range) as f32 * (2.0 * PI) as f32;
        }

        self.grid = pxu::Grid::new(self.pxu.log_branch, self.consts);
        self.cuts = pxu::Cut::get(self.pxu.log_branch, self.consts);

        egui::CentralPanel::default().show(ctx, |ui| {
            let available_size = ui.available_size();
            ui.horizontal(|ui| {
                self.p_plot.draw(
                    ui,
                    available_size * vec2(0.49, 0.49),
                    &mut self.grid,
                    &self.cuts,
                    self.show_dots,
                    &mut self.pxu,
                );

                self.u_plot.draw(
                    ui,
                    available_size * vec2(0.49, 0.49),
                    &mut self.grid,
                    &self.cuts,
                    self.show_dots,
                    &mut self.pxu,
                );
            });
            ui.horizontal(|ui| {
                self.xp_plot.draw(
                    ui,
                    available_size * vec2(0.49, 0.49),
                    &mut self.grid,
                    &self.cuts,
                    self.show_dots,
                    &mut self.pxu,
                );

                self.xm_plot.draw(
                    ui,
                    available_size * vec2(0.49, 0.49),
                    &mut self.grid,
                    &self.cuts,
                    self.show_dots,
                    &mut self.pxu,
                );
            });
        });
    }
}
