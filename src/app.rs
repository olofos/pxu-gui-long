use std::f64::consts::PI;

use egui::style::Margin;
use egui::{pos2, vec2, Color32, Frame, Pos2, Rect, Stroke, Ui, Vec2};
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
}

enum PlotComponent {
    P,
    X,
    U,
}

struct Plot {
    visible_rect: Rect,
    component: PlotComponent,
}

impl Plot {
    fn draw(&mut self, ui: &mut Ui, desired_size: Vec2, grid: &pxu::Grid) {
        let contours = match self.component {
            PlotComponent::P => &grid.p,
            PlotComponent::X => &grid.x,
            PlotComponent::U => &grid.u,
        };

        egui::Frame::canvas(ui.style())
            .outer_margin(Margin::same(0.0))
            .inner_margin(Margin::same(0.0))
            .show(ui, |ui| {
                let (response, painter) = ui.allocate_painter(desired_size, egui::Sense::hover());

                if response.hovered() {
                    let zoom = ui.input().zoom_delta();
                    self.zoom(zoom);
                }

                let rect = response.rect;

                let to_screen = eframe::emath::RectTransform::from_to(self.visible_rect, rect);
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

                    for center in points {
                        shapes.push(egui::epaint::Shape::Circle(egui::epaint::CircleShape {
                            center,
                            radius: 2.0,
                            fill: Color32::RED,
                            stroke: Stroke::NONE,
                        }));
                    }
                }

                painter.extend(shapes);
            });
    }

    fn zoom(&mut self, zoom: f32) {
        let Rect { min, max } = self.visible_rect;
        self.visible_rect = Rect::from_min_max(
            Pos2 {
                x: min.x / zoom,
                y: min.y / zoom,
            },
            Pos2 {
                x: max.x / zoom,
                y: max.y / zoom,
            },
        );
    }
}

impl Default for TemplateApp {
    fn default() -> Self {
        let consts = CouplingConstants::new(2.0, 5);
        let p_range = 0;
        Self {
            consts,
            grid: PxuGrid::new_pm(consts),
            pxu: PxuPoint::new(num::complex::Complex::new(2.0 - 2.0 * PI, 0.0), consts),
            z: num::complex::Complex::new(0.0, 0.5),
            branch: 1,
            new_grid: pxu::Grid::new(p_range, consts),
            p_plot: Plot {
                component: PlotComponent::P,
                visible_rect: Rect::from_x_y_ranges(-7.0..=7.0, -2.0..=2.0),
            },
            x_plot: Plot {
                component: PlotComponent::X,
                visible_rect: Rect::from_x_y_ranges(-4.0..=4.0, -4.0..=4.0),
            },
            u_plot: Plot {
                component: PlotComponent::U,
                visible_rect: Rect::from_x_y_ranges(-4.0..=4.0, -4.0..=4.0),
            },
            p_range,
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

fn extract_grid_component(
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
    fn plot(
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
                    grid.extend(extract_grid_component(&self.grid, extract).into_iter());
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

        // let the_grid = PxuGrid::new_pm(self.consts);

        let s = self.consts.s() as f32;

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

                    self.p_plot
                        .draw(ui, available_size * vec2(0.49, 0.49), &self.new_grid);

                    self.x_plot
                        .draw(ui, available_size * vec2(0.49, 0.49), &self.new_grid);

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

                    self.u_plot
                        .draw(ui, available_size * vec2(0.49, 0.49), &self.new_grid);

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
