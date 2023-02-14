use std::f64::consts::PI;

use egui::style::Margin;
use egui::{vec2, Color32, Pos2, Rect, Stroke, Ui, Vec2};
use num::complex::Complex64;

use crate::kinematics::CouplingConstants;
use crate::pxu;
use crate::pxu::{PxuGrid, PxuPoint};

/// We derive Deserialize/Serialize so we can persist app state on shutdown.
#[derive(serde::Deserialize, serde::Serialize)]
#[serde(default)] // if we add new fields, give them default values when deserializing old state

pub struct TemplateApp {
    #[serde(skip)]
    consts: CouplingConstants,
    #[serde(skip)]
    grid: PxuGrid,
    #[serde(skip)]
    pxu: PxuPoint,
    #[serde(skip)]
    z: num::complex::Complex64,
    #[serde(skip)]
    branch: i32,
    #[serde(skip)]
    new_grid: pxu::Grid,
    #[serde(skip)]
    p_plot: Plot,
    #[serde(skip)]
    x_plot: Plot,
    #[serde(skip)]
    u_plot: Plot,
    #[serde(skip)]
    p_range: i32,
    #[serde(skip)]
    show_dots: bool,
}

enum PlotComponent {
    P,
    X,
    U,
}

struct Plot {
    component: PlotComponent,
    height: f32,
    width_factor: f32,
    origin: Pos2,
}

impl Plot {
    fn draw(
        &mut self,
        ui: &mut Ui,
        desired_size: Vec2,
        grid: &pxu::Grid,
        show_dots: bool,
        pxu: &mut PxuPoint,
    ) {
        let contours = match self.component {
            PlotComponent::P => &grid.p,
            PlotComponent::X => &grid.x,
            PlotComponent::U => &grid.u,
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
                        Stroke::new(1.0, Color32::BLUE),
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

                /*
                let size = egui::epaint::Vec2::splat(8.0);

                    let x = extract(&self.pxu);
                    let center = to_screen * egui::pos2(x.re as f32, -x.im as f32);

                    let point_rect = egui::Rect::from_center_size(center, size);
                    let point_id = response.id.with(i);
                    let point_response = ui.interact(point_rect, point_id, egui::Sense::drag());

                    let mut stroke = egui::epaint::Stroke::NONE;

                    if point_response.hovered() || point_response.dragged() {
                        stroke = egui::epaint::Stroke::new(1.0, Color32::GREEN);
                    }
                    if point_response.dragged() {
                        let new_value =
                            to_screen.inverse() * (center + point_response.drag_delta());

                        self.pxu = set(
                            &self.pxu,
                            num::complex::Complex::new(new_value.x as f64, -new_value.y as f64),
                            self.consts,
                        );
                    }
                    */

                match self.component {
                    PlotComponent::P => {
                        let z = pxu.p;

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
                            let new_value =
                                to_screen.inverse() * (center + point_response.drag_delta());

                            *pxu = PxuPoint::new(
                                Complex64::new(new_value.x as f64, -new_value.y as f64),
                                pxu.consts,
                            );

                            log::info!("E = {}", crate::kinematics::en(pxu.p, 1.0, pxu.consts));
                        }

                        let center = to_screen * egui::pos2(z.re as f32, -z.im as f32);

                        shapes.push(egui::epaint::Shape::Circle(egui::epaint::CircleShape {
                            center,
                            radius,
                            fill: Color32::BLUE,
                            stroke,
                        }));
                    }
                    PlotComponent::X => {
                        let z = pxu.xp;

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

                        if point_response.dragged() {
                            let new_value =
                                to_screen.inverse() * (center + point_response.drag_delta());

                            if let Some(new_pxu) = pxu
                                .shift_xp(Complex64::new(new_value.x as f64, -new_value.y as f64))
                            {
                                *pxu = new_pxu;

                                let tst = pxu.xp.im - pxu.xm.im;

                                // let tst = pxu.xp + 1.0 / pxu.xp
                                //     - pxu.xm
                                //     - 1.0 / pxu.xm
                                //     - 2.0 * Complex64::i() * (1.0 + pxu.consts.kslash * pxu.p)
                                //         / pxu.consts.h;
                                log::info!("{:?}", tst);
                            }
                        }

                        shapes.push(egui::epaint::Shape::Circle(egui::epaint::CircleShape {
                            center: to_screen * egui::pos2(z.re as f32, -z.im as f32),
                            radius: 4.0,
                            fill: Color32::BLUE,
                            stroke,
                        }));

                        let z = pxu.xm;

                        let size = egui::epaint::Vec2::splat(8.0);
                        let center = to_screen * egui::pos2(z.re as f32, -z.im as f32);
                        let point_rect = egui::Rect::from_center_size(center, size);
                        let point_id = response.id.with(1);
                        let point_response = ui.interact(point_rect, point_id, egui::Sense::drag());

                        let stroke = if point_response.hovered() || point_response.dragged() {
                            egui::epaint::Stroke::new(2.0, Color32::LIGHT_BLUE)
                        } else {
                            egui::epaint::Stroke::NONE
                        };

                        if point_response.dragged() {
                            let new_value =
                                to_screen.inverse() * (center + point_response.drag_delta());

                            if let Some(new_pxu) = pxu
                                .shift_xm(Complex64::new(new_value.x as f64, -new_value.y as f64))
                            {
                                *pxu = new_pxu;
                            }
                        }
                        shapes.push(egui::epaint::Shape::Circle(egui::epaint::CircleShape {
                            center: to_screen * egui::pos2(z.re as f32, -z.im as f32),
                            radius: 4.0,
                            fill: Color32::BLUE,
                            stroke,
                        }));
                    }
                    PlotComponent::U => {
                        let z = pxu.u;

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

                        if point_response.dragged() {
                            let new_value =
                                to_screen.inverse() * (center + point_response.drag_delta());

                            if let Some(new_pxu) =
                                pxu.shift_u(Complex64::new(new_value.x as f64, -new_value.y as f64))
                            {
                                *pxu = new_pxu;
                            }
                        }
                        shapes.push(egui::epaint::Shape::Circle(egui::epaint::CircleShape {
                            center: to_screen * egui::pos2(z.re as f32, -z.im as f32),
                            radius: 4.0,
                            fill: Color32::BLUE,
                            stroke,
                        }));
                    }
                };

                painter.extend(shapes);
            });
    }

    fn zoom(&mut self, zoom: f32) {
        self.height /= zoom;
        // let Rect { min, max } = self.visible_rect;
        // self.visible_rect = Rect::from_min_max(
        //     Pos2 {
        //         x: min.x / zoom,
        //         y: min.y / zoom,
        //     },
        //     Pos2 {
        //         x: max.x / zoom,
        //         y: max.y / zoom,
        //     },
        // );
    }
}

impl Default for TemplateApp {
    fn default() -> Self {
        let consts = CouplingConstants::new(2.0, 5);
        let p_range = 0;
        Self {
            consts,
            grid: PxuGrid::new_pm(consts),
            pxu: PxuPoint::new(num::complex::Complex::from(2.0), consts),
            z: num::complex::Complex::new(0.0, 0.5),
            branch: 1,
            new_grid: pxu::Grid::new(p_range, consts),
            p_plot: Plot {
                component: PlotComponent::P,
                height: 2.0,
                width_factor: 4.0,
                origin: Pos2::new(((2 * p_range + 1) as f32) * PI as f32, 0.0),
            },
            x_plot: Plot {
                component: PlotComponent::X,
                height: (8.0 * consts.s()) as f32,
                width_factor: 1.0,
                origin: Pos2::ZERO,
            },
            u_plot: Plot {
                component: PlotComponent::U,
                height: ((2 * consts.k() + 1) as f64 / consts.h) as f32,
                width_factor: 1.0,
                origin: Pos2::ZERO,
            },
            p_range,
            show_dots: true,
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

fn _extract_grid_component(
    the_grid: &PxuGrid,
    extract: &fn(&PxuPoint) -> num::complex::Complex<f64>,
) -> Vec<Vec<(f64, f64)>> {
    the_grid
        .points
        .iter()
        .map(|l| {
            l.iter()
                .map(extract)
                .map(|z| (z.re, z.im))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>()
}

impl TemplateApp {
    fn _plot(
        &mut self,
        ui: &mut Ui,
        desired_size: Vec2,
        extracts: Vec<fn(&PxuPoint) -> num::complex::Complex64>,
        setters: Vec<fn(&PxuPoint, num::complex::Complex64, CouplingConstants) -> PxuPoint>,
        visible_rect: Rect,
    ) -> () {
        egui::Frame::canvas(ui.style())
            .outer_margin(Margin::same(0.0))
            .inner_margin(Margin::same(0.0))
            .show(ui, |ui| {
                let (response, painter) = ui.allocate_painter(desired_size, egui::Sense::hover());

                let zoom;
                {
                    if response.hovered() {
                        let inp = ui.input();
                        zoom = inp.zoom_delta();
                    } else {
                        zoom = 1.0;
                    }
                }

                let rect = response.rect;

                let to_screen =
                    eframe::emath::RectTransform::from_to(visible_rect.expand(zoom), rect);
                ui.set_clip_rect(rect);

                let origin = to_screen * egui::pos2(0.0, 0.0);

                let mut grid = vec![];
                for extract in extracts.iter() {
                    grid.extend(_extract_grid_component(&self.grid, extract).into_iter());
                }

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

                for points in grid {
                    let points = points
                        .iter()
                        .map(|(x, y)| to_screen * egui::pos2(*x as f32, -*y as f32))
                        .collect::<Vec<_>>();

                    shapes.push(egui::epaint::Shape::line(
                        points,
                        Stroke::new(1.0, Color32::BLUE),
                    ));
                }

                for (i, (extract, set)) in extracts.iter().zip(setters.iter()).enumerate() {
                    let size = egui::epaint::Vec2::splat(8.0);

                    let x = extract(&self.pxu);
                    let center = to_screen * egui::pos2(x.re as f32, -x.im as f32);

                    let point_rect = egui::Rect::from_center_size(center, size);
                    let point_id = response.id.with(i);
                    let point_response = ui.interact(point_rect, point_id, egui::Sense::drag());

                    let mut stroke = egui::epaint::Stroke::NONE;

                    if point_response.hovered() || point_response.dragged() {
                        stroke = egui::epaint::Stroke::new(1.0, Color32::GREEN);
                    }
                    if point_response.dragged() {
                        let new_value =
                            to_screen.inverse() * (center + point_response.drag_delta());

                        self.pxu = set(
                            &self.pxu,
                            num::complex::Complex::new(new_value.x as f64, -new_value.y as f64),
                            self.consts,
                        );
                    }

                    // let stroke = ui.style().interact(&point_response).fg_stroke;
                    let x = extract(&self.pxu);
                    let center = to_screen * egui::pos2(x.re as f32, -x.im as f32);

                    shapes.push(egui::epaint::Shape::Circle(egui::epaint::CircleShape {
                        center,
                        radius: 4.0,
                        fill: Color32::RED,
                        stroke,
                    }));
                }

                painter.extend(shapes);
            });
        // });
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

            ui.add(egui::Slider::new(&mut self.consts.h, 0.5..=5.5).text("h"));

            ui.add(
                egui::Slider::from_get_set(0.0..=10.0, |v| self.consts.get_set_k(v))
                    .integer()
                    .text("k"),
            );

            ui.add(
                egui::Slider::from_get_set(1.00001..=5.0, |v| self.consts.get_set_s(v)).text("s"),
            );

            ui.add(egui::Slider::new(&mut self.p_range, -5..=5).text("Range"));

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
            self.grid = PxuGrid::new_pm(self.consts);
            self.pxu = PxuPoint::new(self.pxu.p, self.consts);

            self.new_grid = pxu::Grid::new(self.p_range, self.consts);
        }

        if old_consts != self.consts {
            self.x_plot.height *= (self.consts.s() / old_consts.s()) as f32;

            self.u_plot.height /= (self.consts.h / old_consts.h) as f32;
            if self.consts.k() > 1 && old_consts.k() > 1 {
                self.x_plot.height *= (2 * self.consts.k()) as f32 / (2 * old_consts.k()) as f32;
                self.u_plot.height *= (2 * self.consts.k()) as f32 / (2 * old_consts.k()) as f32;
            }
        }

        if old_p_range != self.p_range {
            self.x_plot.height *= 2.0f32.powi(self.p_range - old_p_range);
            self.p_plot.origin.x += (self.p_range - old_p_range) as f32 * (2.0 * PI) as f32;
        }

        // let the_grid = PxuGrid::new_pm(self.consts);

        egui::CentralPanel::default()
            // .frame(Frame::default().inner_margin(Margin::same(5.0)))
            .show(ctx, |ui| {
                let available_size = ui.available_size();
                ui.horizontal(|ui| {
                    // self.plot(
                    //     // ctx,
                    //     ui,
                    //     available_size * vec2(0.49, 0.49),
                    //     // "xp/xm",
                    //     vec![|pxu: &PxuPoint| pxu.xp, |pxu: &PxuPoint| pxu.xm],
                    //     vec![
                    //         |pxu: &PxuPoint,
                    //          x: num::complex::Complex64,
                    //          consts: CouplingConstants| {
                    //             pxu.shift_xp(x, consts).unwrap_or_else(|| {
                    //                 log::info!("!");
                    //                 pxu.clone()
                    //             })
                    //         },
                    //         |pxu: &PxuPoint,
                    //          x: num::complex::Complex64,
                    //          consts: CouplingConstants| {
                    //             pxu.shift_xm(x, consts).unwrap_or(pxu.clone())
                    //         },
                    //     ],
                    //     eframe::emath::Rect::from_x_y_ranges(
                    //         -2.0 * s..=2.0 * s,
                    //         -1.5 * s..=1.5 * s,
                    //     ),
                    // );

                    self.p_plot.draw(
                        ui,
                        available_size * vec2(0.49, 0.49),
                        &self.new_grid,
                        self.show_dots,
                        &mut self.pxu,
                    );

                    self.x_plot.draw(
                        ui,
                        available_size * vec2(0.49, 0.49),
                        &self.new_grid,
                        self.show_dots,
                        &mut self.pxu,
                    );

                    // ui.with_layout(layout, add_contents)

                    // self.plot(
                    //     //     ctx,
                    //     ui,
                    //     available_size * vec2(0.49, 0.49),
                    //     // "u",
                    //     vec![|pxu: &PxuPoint| pxu.u],
                    //     vec![|pxu: &PxuPoint,
                    //           u: num::complex::Complex64,
                    //           consts: CouplingConstants| {
                    //         pxu.shift_u(u, consts).unwrap_or(pxu.clone())
                    //     }],
                    //     eframe::emath::Rect::from_x_y_ranges(-2.0..=4.0, -3.0..=3.0),
                    // );
                });
                ui.horizontal(|ui| {
                    // self.plot(
                    //     // ctx,
                    //     ui,
                    //     available_size * vec2(0.49, 0.49),
                    //     // "p",
                    //     vec![|pxu: &PxuPoint| pxu.p],
                    //     vec![|_: &PxuPoint,
                    //           p: num::complex::Complex64,
                    //           consts: CouplingConstants| {
                    //         PxuPoint::new(p, consts)
                    //     }],
                    //     eframe::emath::Rect::from_x_y_ranges(-7.0..=7.0, -1.0..=1.0),
                    // );

                    self.u_plot.draw(
                        ui,
                        available_size * vec2(0.49, 0.49),
                        &self.new_grid,
                        self.show_dots,
                        &mut self.pxu,
                    );

                    egui::Frame::canvas(ui.style())
                        .outer_margin(Margin::same(0.0))
                        .inner_margin(Margin::same(0.0))
                        .show(ui, |ui| {
                            let (response, painter) = ui.allocate_painter(
                                available_size * vec2(0.49, 0.49),
                                egui::Sense::hover(),
                            );

                            let rect = response.rect;
                            let to_screen = eframe::emath::RectTransform::from_to(
                                eframe::emath::Rect::from_x_y_ranges(-1.0..=1.0, -1.0..=1.0),
                                rect,
                            );
                            ui.set_clip_rect(rect);

                            let size = egui::epaint::Vec2::splat(8.0);

                            let x = self.z;

                            let cut_points = vec![
                                num::complex::Complex64::new(-0.5, 0.0),
                                num::complex::Complex64::new(0.5, 0.0),
                            ];

                            let center = to_screen * egui::pos2(x.re as f32, -x.im as f32);

                            let point_rect = egui::Rect::from_center_size(center, size);
                            let point_id = response.id.with(1);
                            let point_response =
                                ui.interact(point_rect, point_id, egui::Sense::drag());

                            let mut stroke = egui::epaint::Stroke::NONE;

                            if point_response.hovered() || point_response.dragged() {
                                stroke = egui::epaint::Stroke::new(1.0, Color32::LIGHT_BLUE);
                            }

                            if point_response.dragged() {
                                fn cross(a: Complex64, b: Complex64) -> f64 {
                                    a.re * b.im - a.im * b.re
                                }

                                let new_value =
                                    to_screen.inverse() * (center + point_response.drag_delta());

                                let w = num::complex::Complex64::new(
                                    new_value.x as f64,
                                    -new_value.y as f64,
                                );

                                let p = cut_points[0];
                                let r = cut_points[1] - cut_points[0];

                                let q = self.z;
                                let s = w - self.z;

                                if cross(r, s) != 0.0 {
                                    let t = cross(q - p, s) / cross(r, s);
                                    let u = cross(q - p, r) / cross(r, s);

                                    if 0.0 < t && t < 1.0 && 0.0 < u && u < 1.0 {
                                        self.branch *= -1;
                                    }
                                }

                                self.z = w;
                            }

                            let mut shapes = vec![];

                            for i in 0..cut_points.len() - 1 {
                                let z1 = cut_points[i];
                                let z2 = cut_points[i + 1];

                                let p1 = to_screen * egui::pos2(z1.re as f32, -z1.im as f32);
                                let p2 = to_screen * egui::pos2(z2.re as f32, -z2.im as f32);

                                shapes.push(egui::epaint::Shape::LineSegment {
                                    points: [p1, p2],
                                    stroke: egui::epaint::Stroke::new(
                                        2.0,
                                        if self.branch == 1 {
                                            Color32::GREEN
                                        } else {
                                            Color32::RED
                                        },
                                    ),
                                })
                            }

                            // let stroke = ui.style().interact(&point_response).fg_stroke;
                            let x = self.z;
                            let center = to_screen * egui::pos2(x.re as f32, -x.im as f32);

                            shapes.push(egui::epaint::Shape::Circle(egui::epaint::CircleShape {
                                center,
                                radius: 4.0,
                                fill: Color32::BLUE,
                                stroke,
                            }));

                            painter.extend(shapes);
                        });
                });
            });
    }
}
