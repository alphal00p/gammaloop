use std::collections::{BTreeMap, HashMap};

use cgmath::Point2;
use linnet::half_edge::layout::spring::{Constraint, PointConstraint, ShiftDirection};
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub enum PinConstraint {
    Fixed(f64, f64),
    FixX(f64),
    FixY(f64),
    LinkX(String),
    LinkY(String),
    LinkBoth(String),
    Combined(Box<PinConstraint>, Box<PinConstraint>),
}

impl PinConstraint {
    pub fn point_constraint(
        &self,
        index: usize,
        map: &mut HashMap<String, usize>,
    ) -> (Point2<f64>, PointConstraint) {
        match self {
            PinConstraint::Fixed(x, y) => (
                Point2::new(*x, *y),
                PointConstraint {
                    x: Constraint::Fixed,
                    y: Constraint::Fixed,
                },
            ),
            PinConstraint::FixX(x) => (
                Point2::new(*x, 0.0),
                PointConstraint {
                    x: Constraint::Fixed,
                    y: Constraint::Free,
                },
            ),
            PinConstraint::FixY(y) => (
                Point2::new(0.0, *y),
                PointConstraint {
                    x: Constraint::Free,
                    y: Constraint::Fixed,
                },
            ),
            PinConstraint::LinkX(group) => {
                let (group_name, direction) = Self::parse_direction(group);
                let reference = *map
                    .entry(format!("link_x_{}", group_name))
                    .or_insert_with(|| index);
                (
                    Point2::new(0.0, 0.0),
                    PointConstraint {
                        x: Constraint::Grouped(reference, direction),
                        y: Constraint::Free,
                    },
                )
            }
            PinConstraint::LinkY(group) => {
                let (group_name, direction) = Self::parse_direction(group);
                let reference = *map
                    .entry(format!("link_y_{}", group_name))
                    .or_insert_with(|| index);
                (
                    Point2::new(0.0, 0.0),
                    PointConstraint {
                        y: Constraint::Grouped(reference, direction),
                        x: Constraint::Free,
                    },
                )
            }
            PinConstraint::LinkBoth(group) => {
                let (group_name, direction) = Self::parse_direction(group);
                let reference = *map
                    .entry(format!("link_{}", group_name))
                    .or_insert_with(|| index);
                (
                    Point2::new(0.0, 0.0),
                    PointConstraint {
                        x: Constraint::Grouped(reference, direction),
                        y: Constraint::Grouped(reference, direction),
                    },
                )
            }
            PinConstraint::Combined(x_constraint, y_constraint) => {
                let (mut pos, constraintx) = x_constraint.point_constraint(index, map);
                let (pos_y, constrainty) = y_constraint.point_constraint(index, map);

                pos.y = pos_y.y;
                (
                    pos,
                    PointConstraint {
                        x: constraintx.x,
                        y: constrainty.y,
                    },
                )
            }
        }
    }

    fn parse_direction(group: &str) -> (&str, ShiftDirection) {
        if let Some(stripped) = group.strip_prefix('+') {
            (stripped, ShiftDirection::PositiveOnly)
        } else if let Some(stripped) = group.strip_prefix('-') {
            (stripped, ShiftDirection::NegativeOnly)
        } else {
            (group, ShiftDirection::Any)
        }
    }

    pub fn parse(input: &str) -> Option<Self> {
        let input = input
            .trim()
            .trim_matches('"')
            .trim_matches(|c| c == '(' || c == ')');

        if let Some(group) = Self::parse_link_group(input) {
            return Some(PinConstraint::LinkBoth(group));
        }

        if input.contains(',') || input.split_whitespace().count() == 2 {
            let parts: Vec<&str> = if input.contains(',') {
                input.split(',').map(|s| s.trim()).collect()
            } else {
                input.split_whitespace().collect()
            };

            if parts.len() == 2 {
                let x_constraint = Self::parse_single_constraint(parts[0]);
                let y_constraint = Self::parse_single_constraint(parts[1]);

                match (x_constraint, y_constraint) {
                    (Some(x), Some(y)) => {
                        return Some(PinConstraint::Combined(Box::new(x), Box::new(y)));
                    }
                    _ => {
                        if let (Ok(x), Ok(y)) = (parts[0].parse::<f64>(), parts[1].parse::<f64>()) {
                            return Some(PinConstraint::Fixed(x, y));
                        }
                    }
                }
            }
        }

        Self::parse_single_constraint(input)
    }

    fn parse_single_constraint(input: &str) -> Option<Self> {
        let input = input.trim();

        if let Some(value) = input.strip_prefix("x:") {
            if let Some(group) = Self::parse_link_group(value) {
                Some(PinConstraint::LinkX(group))
            } else if let Ok(x) = value.parse::<f64>() {
                Some(PinConstraint::FixX(x))
            } else {
                None
            }
        } else if let Some(value) = input.strip_prefix("y:") {
            if let Some(group) = Self::parse_link_group(value) {
                Some(PinConstraint::LinkY(group))
            } else if let Ok(y) = value.parse::<f64>() {
                Some(PinConstraint::FixY(y))
            } else {
                None
            }
        } else {
            Self::parse_link_group(input).map(PinConstraint::LinkBoth)
        }
    }

    fn parse_link_group(input: &str) -> Option<String> {
        if let Some(group) = input.strip_prefix('@') {
            Some(group.to_string())
        } else if let Some(group) = input.strip_prefix("+@") {
            Some(format!("+{group}"))
        } else {
            input.strip_prefix("-@").map(|group| format!("-{group}"))
        }
    }
}

impl std::fmt::Display for PinConstraint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PinConstraint::Fixed(x, y) => write!(f, "{},{}", x, y),
            PinConstraint::FixX(x) => write!(f, "x:{}", x),
            PinConstraint::FixY(y) => write!(f, "y:{}", y),
            PinConstraint::LinkX(group) => write!(f, "x:@{}", group),
            PinConstraint::LinkY(group) => write!(f, "y:@{}", group),
            PinConstraint::LinkBoth(group) => write!(f, "@{}", group),
            PinConstraint::Combined(x, y) => write!(f, "{},{}", x, y),
        }
    }
}

pub fn expand_template(template: &str, statements: &BTreeMap<String, String>) -> String {
    let chars = template.chars().collect::<Vec<_>>();
    let mut result = String::new();
    let mut index = 0;

    while index < chars.len() {
        match chars[index] {
            '{' if chars.get(index + 1) == Some(&'{') => {
                result.push('{');
                index += 2;
            }
            '{' => {
                let mut end = index + 1;
                while end < chars.len() && chars[end] != '}' {
                    end += 1;
                }

                if end < chars.len() {
                    let key = chars[index + 1..end].iter().collect::<String>();
                    if let Some(value) = statements.get(&key) {
                        result.push_str(value.trim().trim_matches('"'));
                    } else {
                        result.push('{');
                        result.push_str(&key);
                        result.push('}');
                    }
                    index = end + 1;
                } else {
                    result.push('{');
                    index += 1;
                }
            }
            '}' if chars.get(index + 1) == Some(&'}') => {
                result.push('}');
                index += 2;
            }
            ch => {
                result.push(ch);
                index += 1;
            }
        }
    }

    result
}
