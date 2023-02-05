use egui::{vec2, Color32, Stroke};

use crate::kinematics::CouplingConstants;
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
}

impl Default for TemplateApp {
    fn default() -> Self {
        let consts = CouplingConstants::new(2.0, 5);
        Self {
            consts,
            grid: PxuGrid::new_pm(consts),
            pxu: PxuPoint::new(num::complex::Complex::new(2.0, 0.0), consts),
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
    fn plot_window(
        &mut self,
        ctx: &egui::Context,
        title: impl Into<egui::WidgetText>,
        extracts: Vec<fn(&PxuPoint) -> num::complex::Complex64>,
        setters: Vec<fn(&PxuPoint, num::complex::Complex64, CouplingConstants) -> PxuPoint>,
        the_grid: &PxuGrid,
        visible_rect: eframe::emath::Rect,
    ) -> () {
        egui::Window::new(title)
            .default_size(vec2(512.0, 512.0))
            .show(ctx, |ui| {
                egui::Frame::canvas(ui.style()).show(ui, |ui| {
                    let desired_size = ui.available_width() * vec2(0.5, 0.5);
                    // let (_id, rect) = ui.allocate_space(desired_size);

                    // let painter = ui.painter();

                    let (response, painter) =
                        ui.allocate_painter(desired_size, egui::Sense::hover());

                    let rect = response.rect;

                    let to_screen = eframe::emath::RectTransform::from_to(visible_rect, rect);
                    ui.set_clip_rect(rect);

                    let origin = to_screen * egui::pos2(0.0, 0.0);

                    let mut grid = vec![];
                    for extract in extracts.iter() {
                        grid.extend(extract_grid_component(the_grid, extract).into_iter());
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
                        if point_response.dragged() {
                            let new_value =
                                to_screen.inverse() * (center + point_response.drag_delta());

                            self.pxu = set(
                                &self.pxu,
                                num::complex::Complex::new(new_value.x as f64, -new_value.y as f64),
                                self.consts,
                            );
                        }

                        let stroke = ui.style().interact(&point_response).fg_stroke;

                        let x = extract(&self.pxu);
                        let center = to_screen * egui::pos2(x.re as f32, -x.im as f32);

                        shapes.push(egui::epaint::Shape::Circle(egui::epaint::CircleShape {
                            center,
                            radius: 2.0,
                            fill: Color32::RED,
                            stroke,
                        }));
                    }

                    painter.extend(shapes);
                });
            });
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

        let the_grid = PxuGrid::new_pm(self.consts);
        let s = self.consts.s() as f32;

        egui::CentralPanel::default().show(ctx, |_ui| {});

        self.plot_window(
            ctx,
            "p",
            vec![|pxu: &PxuPoint| pxu.p],
            vec![
                |_: &PxuPoint, p: num::complex::Complex64, consts: CouplingConstants| {
                    PxuPoint::new(p, consts)
                },
            ],
            &the_grid,
            eframe::emath::Rect::from_x_y_ranges(-1.0..=7.0, -1.0..=1.0),
        );

        self.plot_window(
            ctx,
            "xp/xm",
            vec![|pxu: &PxuPoint| pxu.xp, |pxu: &PxuPoint| pxu.xm],
            vec![
                |pxu: &PxuPoint, x: num::complex::Complex64, consts: CouplingConstants| {
                    pxu.shift_xp(x, consts).unwrap_or_else(|| {
                        log::info!("!");
                        pxu.clone()
                    })
                },
                |pxu: &PxuPoint, x: num::complex::Complex64, consts: CouplingConstants| {
                    pxu.shift_xm(x, consts).unwrap_or(pxu.clone())
                },
            ],
            &the_grid,
            eframe::emath::Rect::from_x_y_ranges(-2.0 * s..=4.0 * s, -3.0 * s..=3.0 * s),
        );

        self.plot_window(
            ctx,
            "u",
            vec![|pxu: &PxuPoint| pxu.u],
            vec![
                |pxu: &PxuPoint, u: num::complex::Complex64, consts: CouplingConstants| {
                    pxu.shift_u(u, consts).unwrap_or(pxu.clone())
                },
            ],
            &the_grid,
            eframe::emath::Rect::from_x_y_ranges(-2.0..=4.0, -3.0..=3.0),
        );
    }
}
