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

        if let Some(stripped) = input.strip_prefix('@') {
            return Some(PinConstraint::LinkBoth(stripped.to_string()));
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
            if let Some(stripped) = value.strip_prefix('@') {
                Some(PinConstraint::LinkX(stripped.to_string()))
            } else if let Ok(x) = value.parse::<f64>() {
                Some(PinConstraint::FixX(x))
            } else {
                None
            }
        } else if let Some(value) = input.strip_prefix("y:") {
            if let Some(stripped) = value.strip_prefix('@') {
                Some(PinConstraint::LinkY(stripped.to_string()))
            } else if let Ok(y) = value.parse::<f64>() {
                Some(PinConstraint::FixY(y))
            } else {
                None
            }
        } else {
            input
                .strip_prefix('@')
                .map(|a| PinConstraint::LinkBoth(a.to_string()))
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
    let mut result = template.to_string();
    let mut chars = template.chars().peekable();
    let mut placeholders = Vec::new();
    let mut current_pos = 0;

    while let Some(ch) = chars.next() {
        if ch == '{' {
            let start = current_pos;
            let mut key = String::new();
            let mut found_closing = false;

            for inner_ch in chars.by_ref() {
                current_pos += inner_ch.len_utf8();
                if inner_ch == '}' {
                    found_closing = true;
                    break;
                } else if inner_ch == '{' {
                    break;
                } else {
                    key.push(inner_ch);
                }
            }

            if found_closing && !key.is_empty() {
                placeholders.push((start, current_pos + 1, key));
            }
        }
        current_pos += ch.len_utf8();
    }

    for (start, end, key) in placeholders.iter().rev() {
        if let Some(value) = statements.get(key) {
            let clean_value = value.trim().trim_matches('"');
            result.replace_range(*start..*end, clean_value);
        }
    }

    result
}
