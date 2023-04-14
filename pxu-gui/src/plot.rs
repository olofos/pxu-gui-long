use egui::style::Margin;
use egui::{vec2, Color32, Pos2, Rect, Stroke, Ui, Vec2};
use num::complex::Complex64;

use pxu::kinematics::UBranch;
use pxu::UCutType;

#[derive(serde::Deserialize, serde::Serialize)]
pub struct Plot {
    pub component: pxu::Component,
    pub height: f32,
    pub width_factor: f32,
    pub origin: Pos2,
}

#[allow(clippy::too_many_arguments)]
impl Plot {
    pub fn draw(
        &mut self,
        ui: &mut Ui,
        desired_size: Vec2,
        contours: &mut pxu::Contours,
        show_cuts: bool,
        u_cut_type: UCutType,
        pxu: &mut pxu::State,
    ) {
        egui::Frame::canvas(ui.style())
            .outer_margin(Margin::same(0.0))
            .inner_margin(Margin::same(0.0))
            .show(ui, |ui| {
                let (response, painter) = ui.allocate_painter(desired_size, egui::Sense::drag());

                if response.hovered() {
                    let zoom = ui.input().zoom_delta();
                    self.zoom(zoom);

                    if response.dragged() {
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
                        self.height * self.width_factor * rect.aspect_ratio(),
                        self.height,
                    ),
                );

                let to_screen = eframe::emath::RectTransform::from_to(visible_rect, rect);

                ui.set_clip_rect(rect);

                let origin = to_screen * egui::pos2(0.0, 0.0);

                let mut shapes = if self.component != pxu::Component::P {
                    vec![
                        egui::epaint::Shape::line(
                            vec![
                                egui::pos2(rect.left(), origin.y),
                                egui::pos2(rect.right(), origin.y),
                            ],
                            Stroke::new(0.75, Color32::DARK_GRAY),
                        ),
                        egui::epaint::Shape::line(
                            vec![
                                egui::pos2(origin.x, rect.bottom()),
                                egui::pos2(origin.x, rect.top()),
                            ],
                            Stroke::new(0.75, Color32::DARK_GRAY),
                        ),
                    ]
                } else {
                    vec![]
                };

                let mut hovered_point = None;
                let mut dragged_point = None;

                for j in 0..pxu.points.len() {
                    let z = pxu.points[j].get(self.component);

                    let size = egui::epaint::Vec2::splat(8.0);
                    let center = to_screen * egui::pos2(z.re as f32, -z.im as f32);
                    let point_rect = egui::Rect::from_center_size(center, size);
                    let point_id = response.id.with(j);
                    let point_response = ui.interact(point_rect, point_id, egui::Sense::drag());

                    if point_response.hovered() {
                        hovered_point = Some(j);
                    }

                    if point_response.dragged() {
                        dragged_point = Some(j);
                    }

                    if point_response.dragged() {
                        let new_value =
                            to_screen.inverse() * (center + point_response.drag_delta());
                        let new_value = Complex64::new(new_value.x as f64, -new_value.y as f64);

                        pxu.set_active_point(j);
                        pxu.update(j, self.component, new_value, contours, u_cut_type);
                    }
                }

                let grid_contours = contours.get_grid(self.component);

                for grid_line in grid_contours {
                    let points = grid_line
                        .path
                        .iter()
                        .map(|z| to_screen * egui::pos2(z.re as f32, -z.im as f32))
                        .collect::<Vec<_>>();

                    shapes.push(egui::epaint::Shape::line(
                        points.clone(),
                        Stroke::new(0.75, Color32::GRAY),
                    ));
                }

                let mut branch_point_shapes = vec![];

                if show_cuts {
                    let shift = if self.component == pxu::Component::U {
                        2.0 * (pxu.active_point().sheet_data.log_branch_p * pxu.consts.k()) as f32
                            / pxu.consts.h as f32
                    } else {
                        0.0
                    };

                    let visible_cuts = contours
                        .get_visible_cuts(pxu, self.component, u_cut_type)
                        .collect::<Vec<_>>();

                    let long_cuts = u_cut_type == UCutType::Long;

                    for cut in visible_cuts {
                        let color = match cut.typ {
                            pxu::CutType::E => Color32::BLACK,

                            pxu::CutType::Log(comp) => {
                                if u_cut_type == UCutType::Short
                                    && (comp != cut.component
                                        || pxu.active_point().sheet_data.u_branch
                                            == (UBranch::Between, UBranch::Between))
                                {
                                    continue;
                                } else if comp == pxu::Component::Xp {
                                    Color32::from_rgb(255, 128, 128)
                                } else {
                                    Color32::from_rgb(128, 255, 128)
                                }
                            }

                            pxu::CutType::ULongNegative(comp) => {
                                if !long_cuts {
                                    continue;
                                } else if comp == pxu::Component::Xp {
                                    Color32::from_rgb(255, 0, 0)
                                } else {
                                    Color32::from_rgb(0, 192, 0)
                                }
                            }

                            pxu::CutType::ULongPositive(comp) => {
                                if u_cut_type == UCutType::SemiShort
                                    || u_cut_type == UCutType::Short
                                        && (comp != cut.component
                                            || pxu.active_point().sheet_data.u_branch
                                                == (UBranch::Between, UBranch::Between))
                                {
                                    continue;
                                } else if comp == pxu::Component::Xp {
                                    Color32::from_rgb(255, 0, 0)
                                } else {
                                    Color32::from_rgb(0, 192, 0)
                                }
                            }

                            pxu::CutType::UShortScallion(comp) => {
                                if long_cuts {
                                    continue;
                                } else if comp == pxu::Component::Xp {
                                    Color32::from_rgb(255, 0, 0)
                                } else {
                                    Color32::from_rgb(0, 192, 0)
                                }
                            }

                            pxu::CutType::UShortKidney(comp) => {
                                if long_cuts {
                                    continue;
                                } else if comp == pxu::Component::Xp {
                                    Color32::from_rgb(255, 0, 0)
                                } else {
                                    Color32::from_rgb(0, 192, 0)
                                }
                            }
                            _ => Color32::from_rgb(255, 128, 0),
                        };

                        let period_shifts = if cut.periodic {
                            let period = 2.0 * pxu.consts.k() as f64 / pxu.consts.h;
                            (-5..=5).map(|n| period as f32 * n as f32).collect()
                        } else {
                            vec![0.0]
                        };

                        for period_shift in period_shifts.iter() {
                            for points in cut.paths.iter() {
                                let points = points
                                    .iter()
                                    .map(|z| {
                                        to_screen
                                            * egui::pos2(
                                                z.re as f32,
                                                -(z.im as f32 - shift + period_shift),
                                            )
                                    })
                                    .collect::<Vec<_>>();

                                match cut.typ {
                                    pxu::CutType::UShortKidney(_)
                                    | pxu::CutType::ULongNegative(_) => {
                                        egui::epaint::Shape::dashed_line_many(
                                            &points.clone(),
                                            Stroke::new(3.0, color),
                                            4.0,
                                            4.0,
                                            &mut shapes,
                                        );
                                    }
                                    _ => {
                                        shapes.push(egui::epaint::Shape::line(
                                            points.clone(),
                                            Stroke::new(3.0, color),
                                        ));
                                    }
                                }
                            }

                            if let Some(ref z) = cut.branch_point {
                                let center = to_screen
                                    * egui::pos2(
                                        z.re as f32,
                                        -(z.im as f32 - shift + period_shift),
                                    );
                                branch_point_shapes.push(egui::epaint::Shape::Circle(
                                    egui::epaint::CircleShape {
                                        center,
                                        radius: 4.0,
                                        fill: color,
                                        stroke: Stroke::NONE,
                                    },
                                ));
                            }
                        }
                    }
                }

                for (i, pt) in pxu.points.iter().enumerate() {
                    let is_hovered = matches!(hovered_point, Some(n) if n == i);
                    let is_dragged = matches!(dragged_point, Some(n) if n == i);
                    let is_active = pxu.active_point == i;

                    let z = pt.get(self.component);
                    let center = to_screen * egui::pos2(z.re as f32, -z.im as f32);

                    let radius = if is_hovered || is_dragged {
                        6.0
                    } else if is_active {
                        5.0
                    } else {
                        4.0
                    };

                    let stroke = if is_active {
                        egui::epaint::Stroke::new(2.0, Color32::LIGHT_BLUE)
                    } else {
                        egui::epaint::Stroke::NONE
                    };

                    let fill = if is_active {
                        Color32::BLUE
                    } else if pxu.points[i].same_sheet(
                        pxu.active_point(),
                        self.component,
                        u_cut_type,
                    ) {
                        Color32::BLACK
                    } else {
                        Color32::GRAY
                    };

                    shapes.push(egui::epaint::Shape::Circle(egui::epaint::CircleShape {
                        center,
                        radius,
                        fill,
                        stroke,
                    }));
                }

                for path in pxu.paths.iter() {
                    let points = path
                        .iter()
                        .map(|pt| pt.get(self.component))
                        .map(|z| to_screen * egui::pos2(z.re as f32, -(z.im as f32)))
                        .collect::<Vec<_>>();

                    shapes.push(egui::epaint::Shape::line(
                        points,
                        Stroke::new(3.0, Color32::GRAY),
                    ));
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
                painter.extend(branch_point_shapes);
            });
    }

    fn zoom(&mut self, zoom: f32) {
        self.height /= zoom;
    }
}
