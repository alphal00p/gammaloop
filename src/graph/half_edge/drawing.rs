use core::f64;
use std::fmt::Display;

use cgmath::{
    Angle, Basis2, InnerSpace, Matrix3, Rad, RelativeEq, Rotation, SquareMatrix, Vector2, Zero,
};

#[derive(Clone, Debug)]
pub enum CetzArc {
    CenterRadius {
        center: Vector2<f64>,
        radius: f64,
        start_angle: Rad<f64>,
        end_angle: Rad<f64>,
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
                end_angle,
            } => {
                let start =
                    center + Vector2::new(radius * start_angle.cos(), radius * start_angle.sin());

                let end = center + Vector2::new(radius * end_angle.cos(), radius * end_angle.sin());

                let middle = center
                    + Vector2::new(
                        radius * ((start_angle + end_angle) / 2.).cos(),
                        radius * ((start_angle + end_angle) / 2.).sin(),
                    );

                format!(
                    "cetz.angle.angle(({},{}),({},{}), ({},{}), label: $ alpha $, radius: {})",
                    center.x, center.y, start.x, start.y, end.x, end.y, radius
                )
            }
            CetzArc::Line { start, end } => {
                format!(
                    "hobby(({:.4},{:.4}),({:.4},{:.4}),stroke:stroke,mark: (end: \">\"))",
                    start.x, start.y, end.x, end.y
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
                end_angle,
            } => {
                let start =
                    center + Vector2::new(radius * start_angle.cos(), radius * start_angle.sin());

                let end = center + Vector2::new(radius * end_angle.cos(), radius * end_angle.sin());

                println!("start: {:?}, end: {:?}", start_angle, end_angle);
                println!("middle: {:?}", ((start_angle + end_angle) / 2.));

                let middle = center
                    + Vector2::new(
                        radius * ((start_angle + end_angle) / 2.).cos(),
                        radius * ((start_angle + end_angle) / 2.).sin(),
                    );

                format!(
                    "hobby(({:.4},{:.4}),({:.4},{:.4}),({:.4},{:.4}),stroke:stroke,mark: (end: \">\"))",
                    start.x, start.y, middle.x, middle.y, end.x, end.y
                )
            }
            CetzArc::Line { start, end } => {
                format!(
                    "hobby(({:.4},{:.4}),({:.4},{:.4}),stroke:stroke,mark: (end: \">\"))",
                    start.x, start.y, end.x, end.y
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
    },
    FancyArrow {
        pos: Vector2<f64>,
        arrow_arc: CetzArc,
        label_pos: Vector2<f64>,
        label_angle: Rad<f64>,
    },
}

#[test]
fn circum() {
    let a = Vector2::new(1., 0.);
    let b = Vector2::new(0., 1.);
    let c = Vector2::new(0., -1.);

    let (center, radius) = EdgeGeometry::circumcircle([&b, &c, &a]);

    assert!(center.relative_eq(&Vector2::zero(), 0.1, 1.));
}

pub trait CetzString {
    fn to_cetz(&self) -> String;
}

impl<F: Display> CetzString for Vector2<F> {
    fn to_cetz(&self) -> String {
        format!("({}, {})", self.x, self.y)
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

    let simple = EdgeGeometry::Simple { pos: middle };

    let fancy = simple.clone().to_fancy(source, None, 0., None);

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

    let fancy = simple.clone().to_fancy(source, Some(sink), 0.1, None);
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

    let fancy = simple.to_fancy(source, Some(sink), 0.1, Some((0.1, 0.8)));

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
        mut label: f64,            //label shift
        arrow: Option<(f64, f64)>, //arrow shift and arc percentage
    ) -> Self {
        match self {
            EdgeGeometry::Simple { pos } => {
                if let Some(sink) = sink {
                    if sink != source {
                        let (center, radius) = Self::circumcircle([&source, &pos, &sink]);

                        // println!("{:?}{:?}{:?}", source, pos, sink);

                        // println!("{:?} {:?}", center, radius);

                        let mut angle_between =
                            (source - center).angle((sink - center)).normalize();

                        if angle_between > Rad::full_turn() / 2. {
                            angle_between = Rad::full_turn() - angle_between;
                        }

                        let end_angle = Vector2::unit_x().angle(sink - center);

                        let start_angle = Vector2::unit_x().angle(source - center);

                        let arrow_arc = arrow.map(|(shift, length)| {
                            let angle_shift =
                                ((angle_between - angle_between * length) / 2.).normalize();

                            println!("shift:{:?}", angle_shift);
                            // println!("angle between {:?}", angle_between);

                            let (start, end) = if start_angle > end_angle {
                                (start_angle - angle_shift, end_angle + angle_shift)
                            } else {
                                (start_angle + angle_shift, end_angle - angle_shift)
                            };

                            label += shift;

                            let arc = CetzArc::CenterRadius {
                                center,
                                radius: radius + shift,
                                start_angle: start,
                                end_angle: end,
                            };
                            arc
                        });

                        let center_to_center = pos - center;
                        let center_to_center_norm = center_to_center.normalize();
                        let label_angle = Vector2::unit_x().angle(center_to_center_norm)
                            - <Rad<f64> as Angle>::turn_div_4();
                        let label_pos = center + center_to_center + center_to_center_norm * label;
                        if let Some(arrow) = arrow_arc {
                            EdgeGeometry::FancyArrow {
                                pos,
                                arrow_arc: arrow,
                                label_pos,
                                label_angle,
                            }
                        } else {
                            EdgeGeometry::Fancy {
                                pos,
                                label_pos,
                                label_angle,
                            }
                        }
                    } else {
                        let center = (source + pos) / 2.;

                        let radius = (source - center).magnitude();
                        let angle = Vector2::unit_y().angle(pos - center);
                        let arrow_arc = arrow.map(|(shift, length)| {
                            let start = Rad::zero() + Rad::turn_div_2() * length + angle;
                            let end = Rad::turn_div_2() - Rad::turn_div_2() * length + angle;
                            label += shift;

                            CetzArc::CenterRadius {
                                center,
                                radius: radius + shift,
                                start_angle: start.normalize(),
                                end_angle: end.normalize(),
                            }
                        });

                        let center_to_center = pos - center;
                        let center_to_center_norm = center_to_center.normalize();
                        let label_angle = Vector2::unit_x().angle(center_to_center_norm)
                            - <Rad<f64> as Angle>::turn_div_4();
                        let label_pos = center + center_to_center + center_to_center_norm * label;
                        if let Some(arrow) = arrow_arc {
                            EdgeGeometry::FancyArrow {
                                pos,
                                arrow_arc: arrow,
                                label_pos,
                                label_angle,
                            }
                        } else {
                            EdgeGeometry::Fancy {
                                pos,
                                label_pos,
                                label_angle,
                            }
                        }
                    }
                } else {
                    let edge = pos - source;

                    let right_angle: Basis2<f64> =
                        cgmath::Rotation2::from_angle(Rad(0.5f64 * f64::consts::PI));
                    let shift_norm = right_angle.rotate_vector(edge).normalize();

                    let label_angle = edge.angle(Vector2::unit_x());

                    if let Some((arrow_shift, length_proportion)) = arrow {
                        let label_pos = source + (label + arrow_shift) * shift_norm + 0.5 * edge;
                        let edge_norm = edge.normalize();

                        let arrow_arc = CetzArc::Line {
                            start: source
                                + edge_norm * (1. - length_proportion) / 2.
                                + shift_norm * arrow_shift,
                            end: pos - edge_norm * (1. - length_proportion) / 2.
                                + shift_norm * arrow_shift,
                        };

                        EdgeGeometry::FancyArrow {
                            pos,
                            arrow_arc,
                            label_pos,
                            label_angle,
                        }
                    } else {
                        let label_pos = source + label * shift_norm + 0.5 * edge;
                        EdgeGeometry::Fancy {
                            pos,
                            label_pos,
                            label_angle,
                        }
                    }
                }
            }
            _ => self,
        }
    }
}

pub enum Decoration {
    Coil,
    Wave,
    Zigzag,
    Dashed,
    None,
}

impl Decoration {
    pub fn to_cetz(&self) -> String {
        match self {
            Decoration::Coil => "coil",
            Decoration::Wave => "wave",
            Decoration::Zigzag => "zigzag",
            Decoration::Dashed => "dashed",
            Decoration::None => "none",
        }
        .to_string()
    }
}

pub trait CetzEdge {
    fn label(&self) -> String;
    fn decoration(&self) -> Decoration;
}
