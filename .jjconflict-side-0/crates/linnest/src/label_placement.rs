//! Numeric geometry for the Typst label-placement search.
//!
//! Keep the operation order in sync with CeTZ's matrix/vector arithmetic. Text
//! measurement and candidate enumeration remain in Typst; this batches the
//! per-candidate transforms, collision costs and search without changing positions,
//! costs, or tie-breaking. Typst certifies every exponential used by the search.
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Deserialize)]
struct Frame {
    point: [f64; 2],
    tangent: [f64; 2],
}

#[derive(Deserialize)]
#[serde(transparent)]
struct Transform([[f64; 4]; 4]);

impl Transform {
    fn point(&self, [x, y]: [f64; 2]) -> [f64; 3] {
        if self.0
            == [
                [1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]
        {
            return [x, y, 0.0];
        }
        std::array::from_fn(|i| {
            let [a, b, c, d] = self.0[i];
            a * x + b * y + c * 0.0 + d * 1.0
        })
    }
}

#[derive(Clone, Copy, Deserialize)]
struct Point(
    #[serde(deserialize_with = "crate::deserialize_f64")] f64,
    #[serde(deserialize_with = "crate::deserialize_f64")] f64,
);

impl Point {
    fn lerp(self, other: Self, at: f64) -> Self {
        Self(
            self.0 + (other.0 - self.0) * at,
            self.1 + (other.1 - self.1) * at,
        )
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct CubicSegment {
    start: Point,
    control_start: Point,
    control_end: Point,
    end: Point,
}

impl CubicSegment {
    fn point(&self, at: f64) -> [f64; 2] {
        // Match Kurvst's Typst cubic-point de Casteljau evaluation exactly.
        let ab = self.start.lerp(self.control_start, at);
        let bc = self.control_start.lerp(self.control_end, at);
        let cd = self.control_end.lerp(self.end, at);
        let abc = ab.lerp(bc, at);
        let bcd = bc.lerp(cd, at);
        let Point(x, y) = abc.lerp(bcd, at);
        [x, y]
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct ObstacleSpec {
    segments: Vec<CubicSegment>,
    transform: Transform,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    pad_x: f64,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    pad_y: f64,
}

impl ObstacleSpec {
    fn boxes(&self) -> Vec<Bounds> {
        let mut boxes = Vec::with_capacity(self.segments.len() * 12);
        for segment in &self.segments {
            let mut start = self.transform.point(segment.point(0.0));
            for step in 1..=12 {
                let end = self.transform.point(segment.point(f64::from(step) / 12.0));
                boxes.push(Bounds {
                    left: start[0].min(end[0]) - self.pad_x,
                    right: start[0].max(end[0]) + self.pad_x,
                    bottom: start[1].min(end[1]) - self.pad_y,
                    top: start[1].max(end[1]) + self.pad_y,
                });
                start = end;
            }
        }
        boxes
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct CandidateSpec {
    frames: Vec<Frame>,
    positions: Vec<f64>,
    transform: Transform,
    origin: [f64; 3],
    corners: [[f64; 3]; 4],
    total: f64,
    preferred: f64,
    accuracy: f64,
    side: f64,
    preferred_side: f64,
    path_index: usize,
    path_shift: f64,
    gap: f64,
    clear_box: bool,
    fixed: bool,
}

#[derive(Clone, Copy, Debug, Deserialize, PartialEq, Serialize)]
struct Bounds {
    #[serde(deserialize_with = "crate::deserialize_f64")]
    left: f64,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    right: f64,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    bottom: f64,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    top: f64,
}

impl Bounds {
    fn expanded(self, padding: f64) -> Self {
        Self {
            left: self.left - padding,
            right: self.right + padding,
            bottom: self.bottom - padding,
            top: self.top + padding,
        }
    }

    fn overlap(self, other: Self, padding: f64) -> f64 {
        let width = (self.right.min(other.right) - self.left.max(other.left) + padding).max(0.0);
        if width == 0.0 {
            return 0.0;
        }
        width * (self.top.min(other.top) - self.bottom.max(other.bottom) + padding).max(0.0)
    }

    fn union(self, other: Self) -> Self {
        Self {
            left: self.left.min(other.left),
            right: self.right.max(other.right),
            bottom: self.bottom.min(other.bottom),
            top: self.top.max(other.top),
        }
    }
}

#[derive(Deserialize)]
struct CollisionCandidate {
    bounds: Bounds,
    cost: ciborium::Value,
}

impl CollisionCandidate {
    fn cost(&self) -> Result<f64, String> {
        Cost::from_value(&self.cost).map(Cost::float)
    }
}

#[derive(Deserialize)]
struct CollisionLabel {
    candidates: Vec<CollisionCandidate>,
    edge: Option<ciborium::Value>,
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct Obstacle {
    #[serde(flatten)]
    bounds: Bounds,
    edge: Option<ciborium::Value>,
    #[serde(default)]
    self_loop: bool,
}

impl Obstacle {
    fn belongs_to(&self, label: &CollisionLabel) -> bool {
        use ciborium::Value::{Float, Integer};
        match (&self.edge, &label.edge) {
            (_, None) => false,
            // Typst compares integer and float identifiers numerically.
            (Some(Integer(left)), Some(Float(right)))
            | (Some(Float(right)), Some(Integer(left))) => i128::from(*left) as f64 == *right,
            (left, right) => left == right,
        }
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct CollisionSpec {
    placements: Vec<CollisionLabel>,
    obstacles: Vec<Obstacle>,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    label_padding: f64,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    obstacle_padding: f64,
}

#[derive(Deserialize, Serialize)]
struct CollisionGeometry {
    boxes: Vec<Vec<Bounds>>,
    areas: Vec<Vec<f64>>,
    extents: Vec<Bounds>,
    costs: Vec<Vec<ciborium::Value>>,
}

#[derive(Clone, Copy)]
enum Cost {
    Integer(i64),
    Float(f64),
}

impl Cost {
    fn from_value(value: &ciborium::Value) -> Result<Self, String> {
        match value {
            ciborium::Value::Integer(value) => i64::try_from(*value)
                .map(Self::Integer)
                .map_err(|_| "label costs exceed the Typst integer range".to_owned()),
            ciborium::Value::Float(value) => Ok(Self::Float(*value)),
            _ => Err("label candidate costs must be numbers".to_owned()),
        }
    }

    fn float(self) -> f64 {
        match self {
            Self::Integer(value) => value as f64,
            Self::Float(value) => value,
        }
    }

    fn add(self, other: Self) -> Result<Self, String> {
        match (self, other) {
            (Self::Integer(left), Self::Integer(right)) => left
                .checked_add(right)
                .map(Self::Integer)
                .ok_or_else(|| "label search integer overflow".to_owned()),
            _ => Ok(Self::Float(self.float() + other.float())),
        }
    }

    fn subtract(self, other: Self) -> Result<Self, String> {
        match (self, other) {
            (Self::Integer(left), Self::Integer(right)) => left
                .checked_sub(right)
                .map(Self::Integer)
                .ok_or_else(|| "label search integer overflow".to_owned()),
            _ => Ok(Self::Float(self.float() - other.float())),
        }
    }
}

#[derive(Deserialize)]
#[serde(rename_all = "kebab-case")]
struct SearchMath {
    temperatures: Vec<f64>,
    exponentials: Vec<(f64, f64)>,
    #[serde(deserialize_with = "crate::deserialize_f64")]
    pair_padding: f64,
}

#[derive(Serialize)]
struct SearchResult {
    choices: Vec<usize>,
    costs: Vec<ciborium::Value>,
    missing: Vec<f64>,
}

impl CollisionGeometry {
    fn search(&self, math: &SearchMath) -> Result<SearchResult, String> {
        let count = self.boxes.len();
        if math.temperatures.len() != 84
            || self.areas.len() != count
            || self.extents.len() != count
            || self.costs.len() != count
            || (0..count).any(|i| {
                self.boxes[i].is_empty()
                    || self.areas[i].len() != self.boxes[i].len()
                    || self.costs[i].len() != self.boxes[i].len()
            })
        {
            return Err("invalid label search geometry or temperature schedule".to_owned());
        }
        let costs: Vec<Vec<_>> = self
            .costs
            .iter()
            .map(|costs| costs.iter().map(Cost::from_value).collect::<Result<_, _>>())
            .collect::<Result<_, _>>()?;
        let exponentials: BTreeMap<_, _> = math
            .exponentials
            .iter()
            .map(|&(argument, value)| (argument.to_bits(), value))
            .collect();
        let mut missing = Vec::new();
        let mut pending = BTreeSet::new();
        let mut choices = vec![0; count];
        let mut best_choices = choices.clone();
        let mut energy = Cost::Integer(0);
        let mut best_energy = energy;
        let mut random = 42_u32;
        for (sweep, &temperature) in math.temperatures.iter().enumerate() {
            if sweep == 80 {
                choices.clone_from(&best_choices);
                energy = best_energy;
            }
            for i in 0..count {
                let candidates = self.boxes[i].len();
                if candidates <= 1 {
                    continue;
                }
                random = random.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                let first = if sweep < 80 {
                    (random >> 16) as usize % candidates
                } else {
                    0
                };
                let last = if sweep < 80 { first + 1 } else { candidates };
                let others: Vec<_> = choices
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| *j != i)
                    .map(|(j, &choice)| (self.boxes[j][choice], self.areas[j][choice]))
                    .filter(|(bounds, _)| self.extents[i].overlap(*bounds, math.pair_padding) > 0.0)
                    .collect();
                let score = |choice: usize| {
                    let mut cost = costs[i][choice];
                    for &(bounds, area) in &others {
                        cost = Cost::Float(
                            cost.float()
                                + 4.0 * self.boxes[i][choice].overlap(bounds, math.pair_padding)
                                    / self.areas[i][choice].min(area).max(1e-9),
                        );
                    }
                    cost
                };
                let mut current_score = score(choices[i]);
                for proposal in first..last {
                    if proposal == choices[i] {
                        continue;
                    }
                    let proposal_score = score(proposal);
                    let delta = proposal_score.subtract(current_score)?;
                    random = random.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
                    let accept = if delta.float() < -1e-12 {
                        true
                    } else if sweep < 80 {
                        // Typst's max(0, delta) keeps the integer zero on ties.
                        let argument = if delta.float() <= 0.0 {
                            0.0
                        } else {
                            -delta.float() / temperature
                        };
                        let probability = match exponentials.get(&argument.to_bits()) {
                            Some(&value) => value,
                            None => {
                                if pending.insert(argument.to_bits()) {
                                    missing.push(argument);
                                }
                                // Prediction only. The caller must supply Typst's
                                // results and replay until `missing` is empty.
                                // Every replay extends the verified execution
                                // prefix even if this platform's exp differs.
                                argument.exp()
                            }
                        };
                        f64::from(random) / 4_294_967_296.0 < probability
                    } else {
                        false
                    };
                    if accept {
                        choices[i] = proposal;
                        current_score = proposal_score;
                        energy = energy.add(delta)?;
                        if energy.float() < best_energy.float() - 1e-12 {
                            best_energy = energy;
                            best_choices.clone_from(&choices);
                        }
                    }
                }
            }
        }
        Ok(SearchResult {
            costs: best_choices
                .iter()
                .enumerate()
                .map(|(i, &choice)| self.costs[i][choice].clone())
                .collect(),
            choices: best_choices,
            missing,
        })
    }
}

impl CollisionSpec {
    fn geometry(&self) -> Result<CollisionGeometry, String> {
        let mut geometry = CollisionGeometry {
            boxes: Vec::with_capacity(self.placements.len()),
            areas: Vec::with_capacity(self.placements.len()),
            extents: Vec::with_capacity(self.placements.len()),
            costs: Vec::with_capacity(self.placements.len()),
        };
        for label in &self.placements {
            let boxes: Vec<_> = label
                .candidates
                .iter()
                .map(|candidate| candidate.bounds.expanded(self.label_padding))
                .collect();
            let extent = boxes
                .iter()
                .copied()
                .reduce(Bounds::union)
                .ok_or("labels must have at least one candidate")?;
            let relevant: Vec<_> = self
                .obstacles
                .iter()
                .filter(|obstacle| !obstacle.belongs_to(label) || obstacle.self_loop)
                .collect();
            let mut areas = Vec::with_capacity(boxes.len());
            let mut costs = Vec::with_capacity(boxes.len());
            for (candidate, bounds) in label.candidates.iter().zip(&boxes) {
                let raw = candidate.bounds;
                let area = (raw.right - raw.left) * (raw.top - raw.bottom);
                areas.push(area);
                let mut cost = candidate.cost()?;
                // Keep zero terms too: adding them converts integer costs to
                // floats. Preserve both obstacle order and numeric types.
                for obstacle in &relevant {
                    cost += bounds.overlap(obstacle.bounds, self.obstacle_padding) / area.max(1e-9);
                }
                costs.push(if relevant.is_empty() {
                    candidate.cost.clone()
                } else {
                    ciborium::Value::Float(cost)
                });
            }
            geometry.boxes.push(boxes);
            geometry.areas.push(areas);
            geometry.extents.push(extent);
            geometry.costs.push(costs);
        }
        Ok(geometry)
    }
}

#[derive(Serialize)]
#[serde(rename_all = "kebab-case")]
struct Candidate {
    position: [f64; 2],
    at: f64,
    side: f64,
    path_index: usize,
    path_shift: f64,
    bounds: Bounds,
    cost: f64,
}

impl CandidateSpec {
    fn candidates(&self) -> Result<Vec<Candidate>, String> {
        if self.frames.len() != self.positions.len() {
            return Err("label frames and positions must have equal lengths".to_owned());
        }
        Ok(self
            .positions
            .iter()
            .zip(&self.frames)
            .map(|(&at, frame)| {
                let [tx, ty] = frame.tangent;
                let length = (tx * tx + ty * ty).sqrt();
                let normal = if length <= 1e-9 {
                    [0.0, self.side]
                } else {
                    let scale = self.side / length;
                    [-ty * scale, tx * scale]
                };
                let transformed = self.transform.point(normal);
                let outward: [f64; 3] = std::array::from_fn(|i| transformed[i] - self.origin[i]);
                let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
                let normal_squared = dot(outward, outward);
                let nearest = if self.clear_box && normal_squared > 1e-18 {
                    self.corners
                        .iter()
                        .map(|&corner| dot(corner, outward))
                        .fold(f64::INFINITY, f64::min)
                        / normal_squared
                } else {
                    0.0
                };
                let position = if self.fixed {
                    frame.point
                } else {
                    std::array::from_fn(|i| frame.point[i] + normal[i] * (self.gap - nearest))
                };
                let center = self.transform.point(position);
                let bounds = self
                    .corners
                    .map(|corner| [center[0] + corner[0], center[1] + corner[1]]);
                let relative = (at - self.preferred) / self.total;
                Candidate {
                    position,
                    at,
                    side: self.side,
                    path_index: self.path_index,
                    path_shift: at - self.total / 2.0 + self.path_shift,
                    bounds: Bounds {
                        left: bounds.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min),
                        right: bounds
                            .iter()
                            .map(|p| p[0])
                            .fold(f64::NEG_INFINITY, f64::max),
                        bottom: bounds.iter().map(|p| p[1]).fold(f64::INFINITY, f64::min),
                        top: bounds
                            .iter()
                            .map(|p| p[1])
                            .fold(f64::NEG_INFINITY, f64::max),
                    },
                    cost: (if self.total <= self.accuracy {
                        0.0
                    } else {
                        0.002 * relative.powi(2)
                    }) + if self.side == self.preferred_side {
                        0.0
                    } else {
                        0.00002
                    },
                }
            })
            .collect())
    }
}

pub fn candidates_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: CandidateSpec = ciborium::de::from_reader(arg)
        .map_err(|error| format!("Invalid label candidate geometry: {error}"))?;
    crate::graph_api::encode_cbor(&spec.candidates()?)
}

pub fn collision_geometry_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: CollisionSpec = ciborium::de::from_reader(arg)
        .map_err(|error| format!("Invalid label collision geometry: {error}"))?;
    crate::graph_api::encode_cbor(&spec.geometry()?)
}

pub fn search_bytes(geometry: &[u8], math: &[u8]) -> Result<Vec<u8>, String> {
    let geometry: CollisionGeometry = ciborium::de::from_reader(geometry)
        .map_err(|error| format!("Invalid label collision geometry: {error}"))?;
    let math: SearchMath = ciborium::de::from_reader(math)
        .map_err(|error| format!("Invalid label search math: {error}"))?;
    crate::graph_api::encode_cbor(&geometry.search(&math)?)
}

pub fn obstacle_boxes_bytes(arg: &[u8]) -> Result<Vec<u8>, String> {
    let spec: ObstacleSpec = ciborium::de::from_reader(arg)
        .map_err(|error| format!("Invalid label obstacle geometry: {error}"))?;
    crate::graph_api::encode_cbor(&spec.boxes())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn host_exponentials_override_speculation_and_finish_a_verified_search() {
        let candidate = |left, cost| CollisionCandidate {
            bounds: Bounds {
                left,
                right: left + 1.0,
                bottom: 0.0,
                top: 1.0,
            },
            cost: ciborium::Value::Float(cost),
        };
        // Either single move raises the energy, but taking both lowers it.
        let geometry = CollisionSpec {
            placements: vec![
                CollisionLabel {
                    candidates: vec![candidate(0.0, 0.001), candidate(10.9995, 0.0)],
                    edge: None,
                },
                CollisionLabel {
                    candidates: vec![candidate(10.0, 0.001), candidate(0.9995, 0.0)],
                    edge: None,
                },
            ],
            obstacles: vec![],
            label_padding: 0.0,
            obstacle_padding: 0.0,
        }
        .geometry()
        .unwrap();
        let mut math = SearchMath {
            temperatures: (0..84).map(|sweep| 0.15 * 0.9_f64.powi(sweep)).collect(),
            exponentials: vec![],
            pair_padding: 0.0,
        };
        let predicted = geometry.search(&math).unwrap();
        assert_eq!(predicted.choices, [1, 1]);
        assert!(!predicted.missing.is_empty());
        let mut result = predicted;
        // Deliberately supply host answers that reject every uphill proposal.
        // This forces a different path, exercising correction after a prediction
        // differs, rather than only the usual two identical executions.
        for _ in 0..=160 {
            if result.missing.is_empty() {
                break;
            }
            math.exponentials
                .extend(result.missing.iter().map(|&argument| (argument, 0.0)));
            result = geometry.search(&math).unwrap();
        }
        assert!(result.missing.is_empty());
        assert_eq!(result.choices, [0, 0]);
        assert_eq!(
            result.costs,
            [ciborium::Value::Float(0.001), ciborium::Value::Float(0.001)]
        );
        math.temperatures.pop();
        assert!(geometry.search(&math).is_err());
        assert_eq!(
            Cost::Integer(9_007_199_254_740_993)
                .subtract(Cost::Integer(9_007_199_254_740_992))
                .unwrap()
                .float(),
            1.0
        );
    }

    #[test]
    fn obstacle_boxes_keep_degenerate_segments_and_handle_empty_paths() {
        let mut spec = ObstacleSpec {
            segments: vec![CubicSegment {
                start: Point(2.0, 3.0),
                control_start: Point(2.0, 3.0),
                control_end: Point(2.0, 3.0),
                end: Point(2.0, 3.0),
            }],
            transform: Transform([
                [2.0, 0.5, 0.0, 1.0],
                [-0.25, 3.0, 0.0, -4.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]),
            pad_x: 0.25,
            pad_y: 0.5,
        };
        assert_eq!(
            spec.boxes(),
            vec![
                Bounds {
                    left: 6.25,
                    right: 6.75,
                    bottom: 4.0,
                    top: 5.0,
                };
                12
            ]
        );
        spec.segments.clear();
        assert!(spec.boxes().is_empty());
    }

    #[test]
    fn collision_geometry_preserves_edge_ownership_and_self_loops() {
        let bounds = Bounds {
            left: 0.0,
            right: 1.0,
            bottom: 0.0,
            top: 1.0,
        };
        let edge = Some(ciborium::Value::Integer(7.into()));
        let mut spec = CollisionSpec {
            placements: vec![CollisionLabel {
                candidates: vec![CollisionCandidate {
                    bounds,
                    cost: ciborium::Value::Integer(0.into()),
                }],
                edge: edge.clone(),
            }],
            obstacles: vec![Obstacle {
                bounds,
                edge,
                self_loop: false,
            }],
            label_padding: 0.25,
            obstacle_padding: 0.08,
        };
        let own_edge = spec.geometry().unwrap();
        assert_eq!(own_edge.costs[0][0], ciborium::Value::Integer(0.into()));
        assert_eq!(own_edge.boxes[0][0], bounds.expanded(0.25));
        assert_eq!(own_edge.areas[0], [1.0]);
        spec.obstacles[0].edge = Some(ciborium::Value::Float(7.0));
        assert_eq!(spec.geometry().unwrap().costs, own_edge.costs);
        spec.obstacles[0].self_loop = true;
        let expected = ciborium::Value::Float(1.08 * 1.08);
        assert_eq!(spec.geometry().unwrap().costs[0][0], expected);
        spec.obstacles[0].self_loop = false;
        spec.obstacles[0].edge = None;
        assert_eq!(spec.geometry().unwrap().costs[0][0], expected);
        spec.obstacles[0].bounds.left = 10.0;
        spec.obstacles[0].bounds.right = 11.0;
        assert_eq!(
            spec.geometry().unwrap().costs[0][0],
            ciborium::Value::Float(0.0)
        );
        spec.obstacles.clear();
        assert_eq!(
            spec.geometry().unwrap().costs[0][0],
            ciborium::Value::Integer(0.into())
        );
    }

    #[test]
    fn collision_cost_additions_retain_the_original_order() {
        let bounds = Bounds {
            left: 0.0,
            right: 1.0,
            bottom: 0.0,
            top: 1.0,
        };
        let cost = 9_007_199_254_740_992.0;
        let spec = CollisionSpec {
            placements: vec![CollisionLabel {
                candidates: vec![CollisionCandidate {
                    bounds,
                    cost: ciborium::Value::Float(cost),
                }],
                edge: None,
            }],
            obstacles: (0..3)
                .map(|_| Obstacle {
                    bounds,
                    edge: None,
                    self_loop: false,
                })
                .collect(),
            label_padding: 0.0,
            obstacle_padding: 0.0,
        };
        // Each individual addition rounds back to the initial value. Summing
        // obstacle penalties before adding the initial cost would change it.
        assert_eq!(
            spec.geometry().unwrap().costs[0][0],
            ciborium::Value::Float(cost)
        );
        assert_ne!(cost + 3.0, cost);
    }

    #[test]
    fn transformed_label_clears_its_box_and_preserves_fixed_position() {
        let mut spec = CandidateSpec {
            frames: vec![Frame {
                point: [1.0, 2.0],
                tangent: [3.0, 0.0],
            }],
            positions: vec![5.0],
            transform: Transform([
                [2.0, 0.0, 0.0, 10.0],
                [0.0, 3.0, 0.0, -4.0],
                [0.0, 0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0, 1.0],
            ]),
            origin: [10.0, -4.0, 0.0],
            corners: [
                [-2.0, -3.0, 0.0],
                [2.0, -3.0, 0.0],
                [2.0, 3.0, 0.0],
                [-2.0, 3.0, 0.0],
            ],
            total: 10.0,
            preferred: 5.0,
            accuracy: 1e-6,
            side: 1.0,
            preferred_side: 1.0,
            path_index: 2,
            path_shift: 0.5,
            gap: 0.25,
            clear_box: true,
            fixed: false,
        };
        let candidate = spec.candidates().unwrap().remove(0);
        assert_eq!(candidate.position, [1.0, 3.25]);
        assert_eq!(
            (candidate.bounds.left, candidate.bounds.right),
            (10.0, 14.0)
        );
        assert_eq!(
            (candidate.bounds.bottom, candidate.bounds.top),
            (2.75, 8.75)
        );
        assert_eq!(candidate.path_shift, 0.5);
        assert_eq!(candidate.cost, 0.0);
        spec.fixed = true;
        assert_eq!(spec.candidates().unwrap()[0].position, [1.0, 2.0]);
        spec.positions.clear();
        assert!(matches!(spec.candidates(), Err(error) if error.contains("equal lengths")));
    }
}
