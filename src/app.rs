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
}

impl Default for TemplateApp {
    fn default() -> Self {
        let consts = CouplingConstants::new(2.0, 5);
        Self {
            consts,
            grid: PxuGrid::new_pm(consts),
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
    extract: fn(&PxuPoint) -> num::complex::Complex<f64>,
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
        &self,
        ctx: &egui::Context,
        title: impl Into<egui::WidgetText>,
        extracts: Vec<fn(&PxuPoint) -> num::complex::Complex64>,
        the_grid: &PxuGrid,
        visible_rect: eframe::emath::Rect,
    ) -> () {
        egui::Window::new(title)
            .default_size(vec2(512.0, 512.0))
            .show(ctx, |ui| {
                egui::Frame::canvas(ui.style()).show(ui, |ui| {
                    let desired_size = ui.available_width() * vec2(0.5, 0.5);
                    let (_id, rect) = ui.allocate_space(desired_size);
                    let to_screen = eframe::emath::RectTransform::from_to(visible_rect, rect);

                    ui.set_clip_rect(rect);

                    let origin = to_screen * egui::pos2(0.0, 0.0);

                    let mut grid = vec![];
                    for extract in extracts {
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
                    ui.painter().extend(shapes);
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
                egui::Slider::from_get_set(0.0..=10.0, |v: Option<f64>| {
                    if let Some(v) = v {
                        self.consts.kslash = v / std::f64::consts::TAU;
                    }
                    std::f64::consts::TAU * self.consts.kslash
                })
                .integer()
                .text("k"),
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
        let s = ((self.consts.kslash * self.consts.kslash + self.consts.h * self.consts.h).sqrt()
            + self.consts.kslash)
            / self.consts.h;
        let s = s as f32;

        egui::CentralPanel::default().show(ctx, |_ui| {});

        self.plot_window(
            ctx,
            "p",
            vec![|pxu: &PxuPoint| pxu.p],
            &the_grid,
            eframe::emath::Rect::from_x_y_ranges(-1.0..=7.0, -1.0..=1.0),
        );

        self.plot_window(
            ctx,
            "xp",
            vec![|pxu: &PxuPoint| pxu.xp, |pxu: &PxuPoint| pxu.xm],
            &the_grid,
            eframe::emath::Rect::from_x_y_ranges(-2.0 * s..=4.0 * s, -3.0 * s..=3.0 * s),
        );

        egui::Window::new("Window 2")
            .default_size(vec2(512.0, 512.0))
            .show(ctx, |ui| {
                egui::Frame::canvas(ui.style()).show(ui, |ui| {
                    // ui.ctx().request_repaint();

                    let desired_size = ui.available_width() * vec2(0.5, 0.5);
                    let (_id, rect) = ui.allocate_space(desired_size);
                    let to_screen = eframe::emath::RectTransform::from_to(
                        eframe::emath::Rect::from_x_y_ranges(
                            -2.0 * s..=4.0 * s,
                            -3.0 * s..=3.0 * s,
                        ),
                        rect,
                    );

                    ui.set_clip_rect(rect);

                    let mut shapes = vec![
                        egui::epaint::Shape::line(
                            vec![
                                to_screen * egui::pos2(-f32::INFINITY, 0.0),
                                to_screen * egui::pos2(f32::INFINITY, 0.0),
                            ],
                            Stroke::new(0.75, Color32::GRAY),
                        ),
                        egui::epaint::Shape::line(
                            vec![
                                to_screen * egui::pos2(0.0, -f32::INFINITY),
                                to_screen * egui::pos2(0.0, f32::INFINITY),
                            ],
                            Stroke::new(0.75, Color32::GRAY),
                        ),
                    ];

                    let grid = extract_grid_component(&the_grid, |pxu| pxu.xp);

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

                    let grid = extract_grid_component(&the_grid, |pxu| pxu.xm);

                    for points in grid {
                        let points = points
                            .iter()
                            .map(|(x, y)| to_screen * egui::pos2(*x as f32, -*y as f32))
                            .collect::<Vec<_>>();

                        shapes.push(egui::epaint::Shape::line(
                            points,
                            Stroke::new(1.0, Color32::LIGHT_RED),
                        ));
                    }

                    ui.painter().extend(shapes);
                });
            });

        self.plot_window(
            ctx,
            "u",
            vec![|pxu: &PxuPoint| pxu.u],
            &the_grid,
            eframe::emath::Rect::from_x_y_ranges(-2.0..=4.0, -3.0..=3.0),
        );
    }
}
