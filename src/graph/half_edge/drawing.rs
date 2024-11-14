use core::f64;
use std::{fmt::Display, mem};

use cgmath::{Angle, Basis2, InnerSpace, Matrix3, Rad, Rotation, SquareMatrix, Vector2, Zero};

use super::{layout::FancySettings, Flow};

#[derive(Clone, Debug)]
pub enum CetzArc {
    CenterRadius {
        center: Vector2<f64>,
        radius: f64,
        start_angle: Rad<f64>,
        arc_angle_relative_to_start: Rad<f64>,
    },
    Line {
        start: Vector2<f64>,
        end: Vector2<f64>,
    },
}

impl CetzArc {
    pub fn angle(&self) -> String {
        match self {
            CetzArc::CenterRadius {
                center,
                radius,
                start_angle,
                arc_angle_relative_to_start,
            } => {
                let start =
                    center + Vector2::new(radius * start_angle.cos(), radius * start_angle.sin());

                let end_angle = start_angle + arc_angle_relative_to_start;

                let end = center + Vector2::new(radius * end_angle.cos(), radius * end_angle.sin());

                format!(
                    "cetz.angle.angle({},{},{}, label: $ alpha $, radius: {})",
                    center.to_cetz(),
                    start.to_cetz(),
                    end.to_cetz(),
                    radius
                )
            }
            CetzArc::Line { start, end } => {
                format!(
                    "hobby({},{},stroke:stroke,mark: (end: \">\"))",
                    start.to_cetz(),
                    end.to_cetz()
                )
            }
        }
    }

    pub fn hobby_arrow(&self) -> String {
        match self {
            CetzArc::CenterRadius {
                center,
                radius,
                start_angle,
                arc_angle_relative_to_start,
            } => {
                let start =
                    center + Vector2::new(radius * start_angle.cos(), radius * start_angle.sin());
                let end_angle = start_angle + arc_angle_relative_to_start;
                let end = center + Vector2::new(radius * end_angle.cos(), radius * end_angle.sin());

                // println!("start: {:?}, end: {:?}", start_angle, end_angle);
                // println!("middle: {:?}", ((start_angle + end_angle) / 2.));

                let middle_angle = start_angle + arc_angle_relative_to_start / 2.;

                let middle = center
                    + Vector2::new(radius * (middle_angle).cos(), radius * (middle_angle).sin());

                format!(
                    "cetz.draw.hobby({},{},{},stroke:stroke,mark: (end: \">\"))",
                    start.to_cetz(),
                    middle.to_cetz(),
                    end.to_cetz()
                )
            }
            CetzArc::Line { start, end } => {
                format!(
                    "cetz.draw.hobby({},{},stroke:stroke,mark: (end: \">\"))",
                    start.to_cetz(),
                    end.to_cetz()
                )
            }
        }
    }
}

#[derive(Clone, Debug)]
pub enum EdgeGeometry {
    Fancy {
        pos: Vector2<f64>,
        label_pos: Vector2<f64>,
        label_angle: Rad<f64>,
    },
    Simple {
        pos: Vector2<f64>,
        angle: Rad<f64>,
    },
    FancyArrow {
        pos: Vector2<f64>,
        arrow_arc: CetzArc,
        label_pos: Vector2<f64>,
        label_angle: Rad<f64>,
    },
}

pub trait CetzString {
    fn to_cetz(&self) -> String;
}

impl<F: Display> CetzString for Vector2<F> {
    fn to_cetz(&self) -> String {
        format!("({:.2}, {:.2})", self.x, self.y)
    }
}
impl<F: Display> CetzString for Rad<F> {
    fn to_cetz(&self) -> String {
        format!("{:.2}rad", self.0)
    }
}
#[test]
fn fancyness() {
    let sink = Vector2::new(10., 0.);
    println!("node({})", sink.to_cetz());
    let source = Vector2::new(0., 0.);
    println!("node({})", source.to_cetz());
    let middle = Vector2::new(3., 2.);
    println!("node({})", middle.to_cetz());

    let (center, radius) = EdgeGeometry::circumcircle([&source, &sink, &middle]);

    println!("circle({},radius:{})", center.to_cetz(), radius);

    let simple = EdgeGeometry::Simple {
        pos: middle,
        angle: Rad(0.),
    };

    let fancy = simple.clone().to_fancy(
        source,
        None,
        Flow::Sink,
        &FancySettings {
            label_shift: 0.0,
            arrow_angle_percentage: None,
            arrow_shift: 0.1,
        },
    );

    if let EdgeGeometry::Fancy {
        pos,
        label_pos,
        label_angle,
    } = fancy
    {
        println!("hobby({},{})", source.to_cetz(), pos.to_cetz(),);
        println!(
            "content({},angle:{}rad,[sr])",
            label_pos.to_cetz(),
            label_angle.0
        );
    }

    let fancy = simple.clone().to_fancy(
        source,
        Some(sink),
        Flow::Sink,
        &FancySettings {
            label_shift: 0.1,
            arrow_angle_percentage: None,
            arrow_shift: 0.1,
        },
    );
    if let EdgeGeometry::Fancy {
        pos,
        label_pos,
        label_angle,
    } = fancy
    {
        println!(
            "hobby({},{},{})",
            source.to_cetz(),
            pos.to_cetz(),
            sink.to_cetz()
        );
        println!(
            "content({},angle:{}rad,[sr])",
            label_pos.to_cetz(),
            label_angle.0
        );
    }

    let fancy = simple.to_fancy(
        source,
        Some(sink),
        Flow::Source,
        &FancySettings {
            label_shift: 0.1,
            arrow_angle_percentage: Some(0.8),
            arrow_shift: 0.1,
        },
    );

    if let EdgeGeometry::FancyArrow {
        pos,
        label_pos,
        label_angle,
        arrow_arc,
    } = fancy
    {
        println!(
            "hobby({},{},{})",
            source.to_cetz(),
            pos.to_cetz(),
            sink.to_cetz()
        );
        println!(
            "content({},angle:{}rad,[sr])",
            label_pos.to_cetz(),
            label_angle.0
        );
        println!("{}", arrow_arc.hobby_arrow());
        println!("{}", arrow_arc.angle());
    }
}

impl EdgeGeometry {
    pub fn circumcircle(points: [&Vector2<f64>; 3]) -> (Vector2<f64>, f64) {
        let [a, b, c] = points;

        let [a_norm, b_norm, c_norm]: [f64; 3] = points.map(|a| a.magnitude2());
        let cx =
            Matrix3::new(a_norm, b_norm, c_norm, a.y, b.y, c.y, 1., 1., 1.).determinant() / 2.0;

        let cy = Matrix3::new(a.x, b.x, c.x, a_norm, b_norm, c_norm, 1., 1., 1.).determinant() / 2.;

        let q = Matrix3::new(a.x, b.x, c.x, a.y, b.y, c.y, 1., 1., 1.).determinant();

        let s = Matrix3::new(a.x, b.x, c.x, a.y, b.y, c.y, a_norm, b_norm, c_norm).determinant();

        let center = Vector2 {
            x: cx / q,
            y: cy / q,
        };

        let radius = (s / q + (cx.powi(2) + cy.powi(2)) / (q.powi(2))).sqrt();

        (center, radius)
    }

    pub fn to_fancy(
        self,
        source: Vector2<f64>,
        sink: Option<Vector2<f64>>,
        flow: Flow,
        settings: &FancySettings, //arrow shift and arc percentage
    ) -> Self {
        match self {
            EdgeGeometry::Simple { pos, angle } => {
                if let Some(sink) = sink {
                    if sink != source {
                        let (center, radius) = Self::circumcircle([&source, &pos, &sink]);

                        let start_angle = Vector2::unit_x().angle(source - center);

                        let arrow_arc = settings.arrow_angle_percentage.map(|length| {
                            let angle_between = (source - center).angle(sink - center).normalize();

                            let source_to_middle =
                                (source - center).angle(pos - center).normalize();

                            let (start, arc_angle) = if angle_between > source_to_middle {
                                //trig_dir
                                let angle_shift =
                                    ((angle_between - angle_between * length) / 2.).normalize();

                                let start = start_angle + angle_shift;
                                (start, angle_between - angle_shift * 2.)
                            } else {
                                //clockwise
                                let true_angle_between = Rad::full_turn() - angle_between;
                                let angle_shift =
                                    ((true_angle_between - true_angle_between * length) / 2.)
                                        .normalize();

                                let arc_angle = -true_angle_between + angle_shift * 2.;

                                let start = start_angle - angle_shift;
                                (start, arc_angle)
                            };

                            CetzArc::CenterRadius {
                                center,
                                radius: radius + settings.arrow_shift,
                                start_angle: start,
                                arc_angle_relative_to_start: arc_angle,
                            }
                        });

                        let center_to_center = pos - center;
                        let center_to_center_norm = center_to_center.normalize();
                        let label_pos = center
                            + center_to_center
                            + center_to_center_norm * settings.label_shift();
                        if let Some(arrow) = arrow_arc {
                            EdgeGeometry::FancyArrow {
                                pos,
                                arrow_arc: arrow,
                                label_pos,
                                label_angle: angle,
                            }
                        } else {
                            EdgeGeometry::Fancy {
                                pos,
                                label_pos,
                                label_angle: angle,
                            }
                        }
                    } else {
                        let center = (source + pos) / 2.;

                        let radius = (source - center).magnitude();
                        let angle = Vector2::unit_y().angle(pos - center);
                        let arrow_arc = settings.arrow_angle_percentage.map(|length| {
                            let start = Rad::zero() + Rad::turn_div_2() * length + angle;
                            let end = Rad::turn_div_2() * length;

                            CetzArc::CenterRadius {
                                center,
                                radius: radius + settings.arrow_shift,
                                start_angle: start.normalize(),
                                arc_angle_relative_to_start: end,
                            }
                        });

                        let center_to_center = pos - center;
                        let center_to_center_norm = center_to_center.normalize();

                        let label_pos = center
                            + center_to_center
                            + center_to_center_norm * settings.label_shift();
                        if let Some(arrow) = arrow_arc {
                            EdgeGeometry::FancyArrow {
                                pos,
                                arrow_arc: arrow,
                                label_pos,
                                label_angle: angle,
                            }
                        } else {
                            EdgeGeometry::Fancy {
                                pos,
                                label_pos,
                                label_angle: angle,
                            }
                        }
                    }
                } else {
                    let edge = pos - source;

                    let right_angle: Basis2<f64> =
                        cgmath::Rotation2::from_angle(Rad(0.5f64 * f64::consts::PI));
                    let shift_norm = right_angle.rotate_vector(edge).normalize();

                    // let label_angle = edge.angle(Vector2::unit_x());

                    if let Some(length_proportion) = settings.arrow_angle_percentage {
                        let label_pos = source + (settings.label_shift()) * shift_norm + 0.5 * edge;
                        let edge_norm = edge.normalize();

                        let mut start = source
                            + edge_norm * (1. - length_proportion) / 2.
                            + shift_norm * settings.arrow_shift;

                        let mut end = pos - edge_norm * (1. - length_proportion) / 2.
                            + shift_norm * settings.arrow_shift;

                        if let Flow::Sink = flow {
                            mem::swap(&mut start, &mut end)
                        }

                        let arrow_arc = CetzArc::Line { start, end };

                        EdgeGeometry::FancyArrow {
                            pos,
                            arrow_arc,
                            label_pos,
                            label_angle: angle,
                        }
                    } else {
                        let label_pos = source + settings.label_shift() * shift_norm + 0.5 * edge;
                        EdgeGeometry::Fancy {
                            pos,
                            label_pos,
                            label_angle: angle,
                        }
                    }
                }
            }
            _ => self,
        }
    }
}
#[derive(Debug, Clone, Copy)]
pub enum Decoration {
    Coil,
    Wave,
    Zigzag,
    Dashed,
    Arrow,
    None,
}

impl Decoration {
    pub fn to_cetz(&self) -> String {
        match self {
            Decoration::Coil => "\"coil\"",
            Decoration::Wave => "\"wave\"",
            Decoration::Zigzag => "\"zigzag\"",
            Decoration::Dashed => "\"dashed\"",
            Decoration::Arrow => "\"arrow\"",
            Decoration::None => "\"none\"",
        }
        .to_string()
    }
}

pub trait CetzEdge {
    fn label(&self) -> String;
    fn decoration(&self) -> Decoration;
}
