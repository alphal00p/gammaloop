//! Topology-preserving ImPrEd relaxation of an explicitly embedded carrier.
//!
//! Physical nodes, edge anchors and route samples use one coordinate space. The
//! movement solver preserves native affine constraints jointly with old-state
//! separators. Refinement changes drawing samples only, never graph incidence.

pub mod movement;

use cgmath::{InnerSpace, Vector2, Zero};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

use movement::{AxisConstraint, Halfplane, MovementDofs, PointRecord};

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct ImpredConfig {
    /// Nominal duration of the cooling schedule, in reference steps.
    pub steps: usize,
    /// Maximum integration stride after the initial fine-step settling phase.
    pub step_scale: usize,
    pub target: f64,
    pub repulsion: f64,
    pub attraction: f64,
    pub parallel_attraction_balance: f64,
    pub pull: f64,
    pub pull_balance: f64,
    pub external_max_points: usize,
    pub split_length_ratio: f64,
    pub contract_chord_ratio: f64,
    pub edge_clearance: f64,
    pub node_edge_strength: f64,
    /// Lay out measured labels as tethered points beside their carriers. Off
    /// by default; label-aware rendering opts in.
    #[serde(default = "ImpredConfig::default_labels")]
    pub labels: bool,
    /// Cooling progress to resume from, in `[0, 1)`: a warm pass refines an
    /// already relaxed layout without reheating it.
    #[serde(default)]
    pub warm_start: f64,
    /// Rotate the drawing to its external-pull optimum at every refinement
    /// checkpoint. Overdamped forces turn a whole drawing only very slowly.
    pub level: bool,
}

impl Default for ImpredConfig {
    fn default() -> Self {
        Self {
            steps: 500,
            step_scale: 2,
            target: 2.4,
            repulsion: 2.5,
            attraction: 2.5,
            parallel_attraction_balance: 1.0,
            pull: 0.45,
            pull_balance: 1.0,
            external_max_points: 2,
            split_length_ratio: 1.5,
            contract_chord_ratio: 1.25,
            edge_clearance: 0.4,
            node_edge_strength: 4.0,
            labels: Self::default_labels(),
            warm_start: 0.0,
            level: true,
        }
    }
}

impl ImpredConfig {
    fn default_labels() -> bool {
        false
    }

    fn validate(self) -> Result<(), String> {
        if self.step_scale == 0 {
            return Err("ImPrEd step scale must be positive".into());
        }
        for value in [
            self.repulsion,
            self.attraction,
            self.parallel_attraction_balance,
            self.pull,
            self.pull_balance,
            self.node_edge_strength,
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err("ImPrEd force multipliers must be finite and nonnegative".into());
            }
        }
        if !self.target.is_finite()
            || self.target <= 0.0
            || !self.edge_clearance.is_finite()
            || self.edge_clearance <= 0.0
        {
            return Err("ImPrEd spacing and edge clearance must be finite and positive".into());
        }
        if !self.split_length_ratio.is_finite()
            || !self.contract_chord_ratio.is_finite()
            || self.contract_chord_ratio <= 0.0
            || self.contract_chord_ratio >= self.split_length_ratio
        {
            return Err("ImPrEd refinement requires 0 < contraction < subdivision".into());
        }
        if !(0.0..1.0).contains(&self.warm_start) {
            return Err("ImPrEd warm start must lie in [0, 1)".into());
        }
        if self.external_max_points > 3 {
            return Err("External routes support at most three intermediate points".into());
        }
        Ok(())
    }

    fn integration_stride(self, epoch: usize) -> usize {
        // Large early forces repeatedly hit the movement constraints. Settle
        // with the reference step first, then advance the same cooling schedule
        // faster. Land on every original refinement checkpoint and the final
        // temperature instead of skipping route changes or truncating cooling.
        // Checkpoints run after their update, so that update also stays fine.
        let scale = if epoch < self.steps.div_ceil(5) || epoch.is_multiple_of(25) {
            1
        } else {
            self.step_scale
        };
        scale
            .min(25 - epoch % 25)
            .min(self.steps.saturating_sub(epoch).saturating_sub(1).max(1))
    }
}

/// External momentum flow. Mixed flow uses horizontal pull; a single flow uses
/// radial pull from the centroid of the physical interaction nodes.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub enum ExternalFlow {
    Incoming,
    Outgoing,
}

impl ExternalFlow {
    /// Signed mixed-flow pull: incoming tips go left, outgoing tips right.
    fn horizontal(self, pull: f64) -> f64 {
        match self {
            Self::Incoming => -pull,
            Self::Outgoing => pull,
        }
    }
}

/// A measured label box beside an internal carrier. ImPrEd moves it as a free
/// point tethered to its edge's anchor. It never joins routes or the crossing
/// certificate, so it may cross its own carrier to change sides.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq)]
pub struct EdgeLabel {
    /// Half-width and half-height of the label box.
    pub extents: [f64; 2],
    /// Box center; seeded beside the carrier when absent.
    #[serde(default)]
    pub center: Option<[f64; 2]>,
    /// Hold a drawn label where it was placed: at its arc fraction and signed
    /// normal offset on its own carrier, without sliding or changing sides.
    #[serde(default)]
    pub pinned: bool,
    /// The pinned placement, `[arc fraction, normal offset]`, taken from
    /// `center` when solving starts.
    #[serde(default)]
    pub pin: Option<[f64; 2]>,
}

/// Fractions of the spacing: the clear gap between a carrier and its label's
/// tether rest, and the contact distance a label keeps from other geometry.
const LABEL_GAP: f64 = 0.25;
const LABEL_CONTACT: f64 = 0.1;
/// Label tethers and contacts relative to the edge attraction and repulsion.
const LABEL_STIFFNESS: f64 = 3.0;

/// A label's position relative to its own carrier polyline.
struct LabelFrame {
    foot: Vector2<f64>,
    tangent: Vector2<f64>,
    normal: Vector2<f64>,
    arc: f64,
    length: f64,
}

/// Extent of a box with half-extents `[x, y]` along a unit direction.
fn box_support([x, y]: [f64; 2], direction: Vector2<f64>) -> f64 {
    x * direction.x.abs() + y * direction.y.abs()
}

/// A dense, portable view of existing native layout coordinates.
///
/// Routes follow physical source-to-target orientation. Internal anchors occur
/// once inside a route; external anchors are its dangling endpoint. `samples`
/// distinguishes drawing samples from physical vertices and preserves their
/// fixed charge budget, including an anchor selected from a seed route knot.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ImpredLayout {
    pub positions: Vec<[f64; 2]>,
    pub routes: Vec<Vec<usize>>,
    pub node_points: Vec<usize>,
    pub anchor_points: Vec<usize>,
    pub external: Vec<Option<ExternalFlow>>,
    pub constraints: Vec<PointRecord>,
    pub pull_scales: Vec<f64>,
    pub samples: Vec<bool>,
    pub nodes_fixed: bool,
    /// Per-edge labels, empty without labels. External edges carry none.
    #[serde(default)]
    pub labels: Vec<Option<EdgeLabel>>,
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ImpredReport {
    /// Actual force/projection iterations; cooling may advance multiple steps.
    pub iterations: usize,
    pub geometry_checks: usize,
    pub inserted_points: usize,
    pub removed_points: usize,
    pub initial_crossings: usize,
    pub max_primal_residual: f64,
    pub max_stationarity_residual: f64,
    pub max_complementarity_residual: f64,
    pub max_separator_residual: f64,
    pub max_cap_residual: f64,
}

#[derive(Clone, Copy, Debug)]
struct Segment {
    edge: usize,
    index: usize,
    a: usize,
    b: usize,
}

#[derive(Clone, Copy)]
struct PointSegment {
    point: usize,
    segment: usize,
    t: f64,
    distance: f64,
    difference: Vector2<f64>,
}

impl PointSegment {
    fn normal(&self) -> Vector2<f64> {
        self.difference / self.distance.max(1e-15)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Geometry {
    crossings: Vec<((usize, usize), (usize, usize))>,
    contacts: usize,
    overlaps: usize,
    degenerate: usize,
}

impl Geometry {
    fn certify(&self, initial: &Self) -> Result<(), String> {
        if self.contacts > 0
            || self.overlaps > 0
            || self.degenerate > 0
            || self.crossings != initial.crossings
        {
            return Err(format!("ImPrEd carrier topology changed: {self:?}"));
        }
        Ok(())
    }
}

struct Model {
    segments: Vec<Segment>,
    pair_weights: Vec<Vec<f64>>,
    attraction_weights: Vec<f64>,
    /// Points on each edge's route; a label never contacts its own carrier.
    route_points: Vec<BTreeSet<usize>>,
}

impl ImpredLayout {
    fn point(&self, point: usize) -> Vector2<f64> {
        self.positions[point].into()
    }

    fn validate(&self) -> Result<(), String> {
        let count = self.positions.len();
        if self.constraints.len() != count
            || self.samples.len() != count
            || self.routes.len() != self.anchor_points.len()
            || self.routes.len() != self.external.len()
            || self.routes.len() != self.pull_scales.len()
        {
            return Err("ImPrEd coordinate and physical ownership arrays disagree".into());
        }
        if !self.labels.is_empty() && self.labels.len() != self.routes.len() {
            return Err("ImPrEd labels must cover every edge".into());
        }
        for (edge, label) in self.labels.iter().enumerate() {
            let Some(label) = label else { continue };
            if self.external[edge].is_some()
                || label.extents.iter().any(|x| !x.is_finite() || *x < 0.0)
                || label.center.is_some_and(|c| c.iter().any(|x| !x.is_finite()))
            {
                return Err("ImPrEd labels need finite boxes on internal edges".into());
            }
        }
        if self.positions.iter().flatten().any(|x| !x.is_finite()) {
            return Err("ImPrEd positions must be finite".into());
        }
        if self
            .node_points
            .iter()
            .chain(&self.anchor_points)
            .any(|p| *p >= count)
        {
            return Err("Invalid physical point identity".into());
        }
        let nodes: BTreeSet<_> = self.node_points.iter().copied().collect();
        if nodes.len() != self.node_points.len() {
            return Err("Physical vertices need unique drawing identities".into());
        }
        let mut owned = vec![false; count];
        for &node in &self.node_points {
            owned[node] = true;
        }
        for (edge, route) in self.routes.iter().enumerate() {
            if route.len() < 2
                || route.iter().any(|p| *p >= count)
                || route
                    .iter()
                    .filter(|&&p| p == self.anchor_points[edge])
                    .count()
                    != 1
            {
                return Err("Invalid ImPrEd route or anchor identity".into());
            }
            let anchor = self.anchor_points[edge];
            if self.external[edge].is_some() {
                if anchor != route[0] && anchor != *route.last().unwrap() {
                    return Err("External anchor must be a route endpoint".into());
                }
                if !self.pull_scales[edge].is_finite() || self.pull_scales[edge] <= 0.0 {
                    return Err("External pull scales must be positive finite values".into());
                }
            } else if anchor == route[0] || anchor == *route.last().unwrap() {
                return Err("Internal anchor must be inside its route".into());
            }
            for &point in route {
                if !nodes.contains(&point) {
                    if owned[point] {
                        return Err("A route point has multiple physical owners".into());
                    }
                    owned[point] = true;
                }
            }
        }
        if owned.iter().any(|x| !x) {
            return Err("ImPrEd contains unowned drawing points".into());
        }
        Ok(())
    }

    fn segments(&self) -> Vec<Segment> {
        self.routes
            .iter()
            .enumerate()
            .flat_map(|(edge, route)| {
                route.windows(2).enumerate().map(move |(index, p)| Segment {
                    edge,
                    index,
                    a: p[0],
                    b: p[1],
                })
            })
            .collect()
    }

    fn geometry(&self) -> Geometry {
        let segments = self.segments();
        let mut out = Geometry {
            crossings: Vec::new(),
            contacts: 0,
            overlaps: 0,
            degenerate: 0,
        };
        let orient = |a: Vector2<f64>, b: Vector2<f64>, c: Vector2<f64>| (b - a).perp_dot(c - a);
        let on = |a: Vector2<f64>, b: Vector2<f64>, p: Vector2<f64>| {
            orient(a, b, p).abs() <= 1e-8
                && (0..2).all(|k| p[k] >= a[k].min(b[k]) - 1e-8 && p[k] <= a[k].max(b[k]) + 1e-8)
        };
        for (i, first) in segments.iter().enumerate() {
            let (a, b) = (self.point(first.a), self.point(first.b));
            if (a - b).magnitude() <= 1e-8 {
                out.degenerate += 1;
            }
            for second in &segments[i + 1..] {
                let (c, d) = (self.point(second.a), self.point(second.b));
                let values = [
                    orient(a, b, c),
                    orient(a, b, d),
                    orient(c, d, a),
                    orient(c, d, b),
                ];
                if values.iter().all(|x| x.abs() <= 1e-8) {
                    let axis = usize::from((a.x - b.x).abs() < (a.y - b.y).abs());
                    let extent = a[axis].max(b[axis]).min(c[axis].max(d[axis]))
                        - a[axis].min(b[axis]).max(c[axis].min(d[axis]));
                    if extent > 1e-8 {
                        out.overlaps += 1;
                        continue;
                    }
                }
                let signs = values.map(|x| if x.abs() <= 1e-8 { 0.0 } else { x.signum() });
                if signs[0] * signs[1] < 0.0 && signs[2] * signs[3] < 0.0 {
                    out.crossings
                        .push(((first.edge, first.index), (second.edge, second.index)));
                } else if first.a != second.a
                    && first.a != second.b
                    && first.b != second.a
                    && first.b != second.b
                    && (on(a, b, c) || on(a, b, d) || on(c, d, a) || on(c, d, b))
                {
                    out.contacts += 1;
                }
            }
        }
        out
    }

    fn model(&self, flexible: &[bool], balance: f64) -> Model {
        let count = self.positions.len();
        let segments = self.segments();
        let mut charges: Vec<_> = self
            .samples
            .iter()
            .map(|sample| if *sample { 0.35 } else { 1.0 })
            .collect();
        let mut pairs = BTreeMap::new();
        for (edge, route) in self.routes.iter().enumerate() {
            let samples: Vec<_> = route
                .iter()
                .copied()
                .filter(|&point| {
                    self.samples[point]
                        || (self.external[edge].is_some() && point == self.anchor_points[edge])
                })
                .collect();
            if !samples.is_empty() {
                let weight = if self.external[edge].is_some() {
                    1.0
                } else {
                    0.35
                } / samples.len() as f64;
                for point in samples {
                    charges[point] = weight;
                }
            }
            if self.external[edge].is_none() && route[0] != *route.last().unwrap() {
                let pair = (
                    route[0].min(*route.last().unwrap()),
                    route[0].max(*route.last().unwrap()),
                );
                *pairs.entry(pair).or_insert(0usize) += 1;
            }
        }
        let mut pair_weights = vec![vec![0.0; count]; count];
        for i in 0..count {
            for j in 0..count {
                if i != j {
                    pair_weights[i][j] = charges[i] * charges[j];
                }
            }
        }
        for (edge, route) in self
            .routes
            .iter()
            .enumerate()
            .filter(|(edge, _)| flexible[*edge])
        {
            if self.external[edge].is_some() {
                for (i, &first) in route.iter().enumerate() {
                    for &second in &route[i + 1..] {
                        pair_weights[first][second] = 0.0;
                        pair_weights[second][first] = 0.0;
                    }
                }
            } else {
                for pair in route.windows(2) {
                    pair_weights[pair[0]][pair[1]] = 0.0;
                    pair_weights[pair[1]][pair[0]] = 0.0;
                }
            }
        }
        let attraction_weights = segments
            .iter()
            .map(|segment| {
                let route = &self.routes[segment.edge];
                if self.external[segment.edge].is_none() && route[0] != *route.last().unwrap() {
                    let pair = (
                        route[0].min(*route.last().unwrap()),
                        route[0].max(*route.last().unwrap()),
                    );
                    (pairs[&pair] as f64).powf(-balance)
                } else {
                    1.0
                }
            })
            .collect();
        Model {
            segments,
            pair_weights,
            attraction_weights,
            route_points: self
                .routes
                .iter()
                .map(|route| route.iter().copied().collect())
                .collect(),
        }
    }

    fn point_segments(&self, model: &Model) -> Vec<PointSegment> {
        let mut result =
            Vec::with_capacity(self.positions.len().saturating_sub(2) * model.segments.len());
        let geometry: Vec<_> = model
            .segments
            .iter()
            .map(|segment| {
                let a = self.point(segment.a);
                let ab = self.point(segment.b) - a;
                (a, ab, ab.magnitude2().max(1e-20))
            })
            .collect();
        for point in 0..self.positions.len() {
            for (index, segment) in model.segments.iter().enumerate() {
                if point == segment.a || point == segment.b {
                    continue;
                }
                let (a, ab, squared) = geometry[index];
                let t = ((self.point(point) - a).dot(ab) / squared).clamp(0.0, 1.0);
                let difference = self.point(point) - (a + t * ab);
                let distance = difference.magnitude();
                result.push(PointSegment {
                    point,
                    segment: index,
                    t,
                    distance,
                    difference,
                });
            }
        }
        result
    }

    fn segment_forces(
        &self,
        model: &Model,
        config: ImpredConfig,
        exponent: f64,
    ) -> Vec<Vector2<f64>> {
        let mut force = vec![Vector2::zero(); self.positions.len()];
        for (segment, &weight) in model.segments.iter().zip(&model.attraction_weights) {
            let vector = self.point(segment.b) - self.point(segment.a);
            let attraction = config.attraction
                * weight
                * (vector.magnitude() / config.target).powf(exponent)
                * vector;
            force[segment.a] += attraction;
            force[segment.b] -= attraction;
        }
        force
    }

    fn forces(
        &self,
        model: &Model,
        pairs: &[PointSegment],
        config: ImpredConfig,
        progress: f64,
        pull: &[f64],
    ) -> Vec<Vector2<f64>> {
        let rep_exponent = 2.0 + 2.0 * progress;
        let mut force = vec![Vector2::zero(); self.positions.len()];
        // Visit each symmetric pair once, while each point still accumulates
        // neighbors in ascending order (including its original diagonal term).
        // Every finite self-pair has the same zero displacement and weight.
        // Keep its addition (and exceptional coefficient values) in order, but
        // evaluate the shared power only once per epoch.
        let diagonal =
            config.repulsion * (0.0 * (config.target / (config.target * 0.02)).powf(rep_exponent));
        for i in 0..self.positions.len() {
            force[i] += diagonal * Vector2::zero();
            for j in i + 1..self.positions.len() {
                let difference = self.point(i) - self.point(j);
                let weight = model.pair_weights[i][j];
                // With a positive distance floor the ratio is bounded (at most
                // 100 even when the floor is subnormal), and the exponent is
                // between 2 and 4. A zero weight therefore has a zero coefficient.
                // Retain the vector additions, including their signed zeros and
                // overflow behavior; keep the original arithmetic if the floor
                // itself underflows to zero.
                let coefficient = if weight == 0.0 && config.target * 0.02 > 0.0 {
                    weight
                } else {
                    weight
                        * (config.target / difference.magnitude().max(config.target * 0.02))
                            .powf(rep_exponent)
                };
                force[i] += config.repulsion * coefficient * difference;
                // Recompute the reversed difference rather than negating:
                // equal coordinates can carry different signed zeros.
                force[j] += config.repulsion * coefficient * (self.point(j) - self.point(i));
            }
        }
        for (total, attraction) in
            force
                .iter_mut()
                .zip(self.segment_forces(model, config, 1.0 - 0.6 * progress))
        {
            *total += attraction;
        }
        let clearance = config.target * config.edge_clearance;
        for pair in pairs {
            if pair.t > 0.0 && pair.t < 1.0 && pair.distance < clearance {
                let repulsion = config.node_edge_strength
                    * config.target
                    * (1.0 - pair.distance / clearance)
                        .max(0.0)
                        .powf(rep_exponent)
                    * pair.normal();
                let segment = model.segments[pair.segment];
                force[pair.point] += repulsion;
                force[segment.a] -= (1.0 - pair.t) * repulsion;
                force[segment.b] -= pair.t * repulsion;
            }
        }
        let mixed = self.mixed_flow();
        let center = self.node_center();
        let mut reaction = Vector2::zero();
        for (edge, flow) in self.external.iter().enumerate() {
            if let Some(flow) = flow {
                let anchor = self.anchor_points[edge];
                let outward = if mixed {
                    Vector2::new(flow.horizontal(pull[edge]), 0.0)
                } else {
                    let vector = self.point(anchor) - center;
                    let length = vector.magnitude();
                    if length <= 1e-9 {
                        continue;
                    }
                    pull[edge] * vector / length
                };
                force[anchor] += outward;
                reaction += outward;
            }
        }
        if !self.node_points.is_empty() {
            for &point in &self.node_points {
                force[point] -= reaction / self.node_points.len() as f64;
            }
        }
        force
    }

    // The movement certificate also bounds each component by 2*cap. If both
    // normal components have magnitude at most one, their rounded dot product
    // cannot exceed 4*cap. Such a distant separator has nonpositive residual
    // and needs neither matrix storage nor the later residual scan. Infinite
    // cap retains all separators for geometry inspection and reference checks.
    fn separators(&self, model: &Model, pairs: &[PointSegment], cap: f64) -> Vec<Halfplane> {
        let mut lengths = Vec::with_capacity(model.segments.len());
        let perpendiculars: Vec<_> = model
            .segments
            .iter()
            .map(|segment| {
                let vector = self.point(segment.b) - self.point(segment.a);
                let length = vector.magnitude();
                lengths.push(length);
                Vector2::new(-vector.y, vector.x) / length.max(1e-30)
            })
            .collect();
        lengths.sort_by(f64::total_cmp);
        let median = if lengths.is_empty() {
            1.0
        } else if lengths.len() % 2 == 1 {
            lengths[lengths.len() / 2]
        } else {
            (lengths[lengths.len() / 2 - 1] + lengths[lengths.len() / 2]) * 0.5
        };
        let clearance = 1e-6 * median.max(1.0);
        let count = self.positions.len();
        let mut rows = Vec::with_capacity(3 * pairs.len() + count * count.saturating_sub(1));
        for pair in pairs {
            let segment = model.segments[pair.segment];
            let perpendicular = perpendiculars[pair.segment];
            let interior = pair.t > 0.0 && pair.t < 1.0;
            let signed = if interior {
                (self.point(pair.point) - self.point(segment.a)).dot(perpendicular)
            } else {
                0.0
            };
            let distance = if interior {
                signed.abs()
            } else {
                pair.distance
            };
            let room = 0.45 * (distance - clearance).max(0.0);
            // Both endpoint bounds are at least this shared room, so the
            // same certificate omits all three rows before endpoint arithmetic.
            if cap.is_finite()
                && room >= 4.0 * cap
                && if interior {
                    perpendicular.x.abs() <= 1.0 && perpendicular.y.abs() <= 1.0
                } else {
                    // A finite positive divisor at least as large as either
                    // component guarantees the rounded normal is bounded by
                    // one, without evaluating either division.
                    let divisor = pair.distance.max(1e-15);
                    divisor.is_finite()
                        && pair.difference.x.abs() <= divisor
                        && pair.difference.y.abs() <= divisor
                }
            {
                continue;
            }
            let normal = if interior {
                perpendicular * if signed >= 0.0 { 1.0 } else { -1.0 }
            } else {
                pair.normal()
            };
            rows.push(Halfplane {
                point: pair.point,
                normal: (-normal).into(),
                bound: room,
            });
            for endpoint in [segment.a, segment.b] {
                let offset = normal.dot(self.point(pair.point) - self.point(endpoint));
                let available = if interior {
                    room
                } else {
                    (offset - distance).max(0.0) + room
                };
                rows.push(Halfplane {
                    point: endpoint,
                    normal: normal.into(),
                    bound: available,
                });
            }
        }
        for i in 0..self.positions.len() {
            for j in 0..self.positions.len() {
                if i != j {
                    let difference = self.point(i) - self.point(j);
                    let distance = difference.magnitude();
                    let bound = 0.45 * (distance - clearance).max(0.0);
                    let divisor = distance.max(1e-30);
                    if cap.is_finite()
                        && bound >= 4.0 * cap
                        && divisor.is_finite()
                        && difference.x.abs() <= divisor
                        && difference.y.abs() <= divisor
                    {
                        continue;
                    }
                    let normal = -difference / divisor;
                    rows.push(Halfplane {
                        point: i,
                        normal: normal.into(),
                        bound,
                    });
                }
            }
        }
        rows
    }

    /// Unit normal of an internal carrier at its anchor, left of the route.
    fn carrier_normal(&self, edge: usize) -> Option<Vector2<f64>> {
        let route = &self.routes[edge];
        let slot = route.iter().position(|&p| p == self.anchor_points[edge])?;
        let tangent = self.point(route[(slot + 1).min(route.len() - 1)])
            - self.point(route[slot.saturating_sub(1)]);
        let length = tangent.magnitude();
        (length > 0.0).then(|| Vector2::new(-tangent.y, tangent.x) / length)
    }

    /// The tether rest of a label: beside its anchor, on `side` of the carrier.
    fn label_rest(
        &self,
        edge: usize,
        label: &EdgeLabel,
        normal: Vector2<f64>,
        side: f64,
        config: ImpredConfig,
    ) -> Vector2<f64> {
        let offset = box_support(label.extents, normal) + LABEL_GAP * config.target;
        self.point(self.anchor_points[edge]) + side * offset * normal
    }

    /// Seat every label at its tether rest, keeping its side of the carrier.
    /// New labels face away from the vertex centroid.
    fn seat_labels(&mut self, config: ImpredConfig) {
        let center = self.node_center();
        for edge in 0..self.labels.len() {
            let (Some(label), Some(normal)) = (self.labels[edge], self.carrier_normal(edge)) else {
                continue;
            };
            if label.pinned && label.center.is_some() {
                continue;
            }
            let anchor = self.point(self.anchor_points[edge]);
            let away = label.center.map_or(anchor - center, |c| Vector2::from(c) - anchor);
            let side = if away.dot(normal) < 0.0 { -1.0 } else { 1.0 };
            let rest = self.label_rest(edge, &label, normal, side, config);
            self.labels[edge] = Some(EdgeLabel {
                center: Some(rest.into()),
                ..label
            });
        }
    }

    /// The nearest point on an edge's own carrier polyline to `center`, with
    /// the route's unit tangent, left normal, arc position and total length.
    fn label_frame(&self, edge: usize, center: Vector2<f64>) -> Option<LabelFrame> {
        let mut arc = 0.0;
        let mut best: Option<(f64, LabelFrame)> = None;
        for pair in self.routes[edge].windows(2) {
            let a = self.point(pair[0]);
            let ab = self.point(pair[1]) - a;
            let length = ab.magnitude();
            if length > 0.0 {
                let t = ((center - a).dot(ab) / (length * length)).clamp(0.0, 1.0);
                let foot = a + t * ab;
                let distance = (center - foot).magnitude();
                if best.as_ref().is_none_or(|(nearest, _)| distance < *nearest) {
                    let tangent = ab / length;
                    best = Some((
                        distance,
                        LabelFrame {
                            foot,
                            tangent,
                            normal: Vector2::new(-tangent.y, tangent.x),
                            arc: arc + t * length,
                            length: 0.0,
                        },
                    ));
                }
            }
            arc += length;
        }
        best.map(|(_, frame)| LabelFrame {
            length: arc,
            ..frame
        })
    }

    /// The point, unit tangent and left normal at an arc distance along an
    /// edge's route polyline.
    fn route_point_at(&self, edge: usize, arc: f64) -> Option<(Vector2<f64>, Vector2<f64>, Vector2<f64>)> {
        let mut remaining = arc.max(0.0);
        let mut last = None;
        for pair in self.routes[edge].windows(2) {
            let a = self.point(pair[0]);
            let ab = self.point(pair[1]) - a;
            let length = ab.magnitude();
            if length <= 0.0 {
                continue;
            }
            let tangent = ab / length;
            let normal = Vector2::new(-tangent.y, tangent.x);
            if remaining <= length {
                return Some((a + remaining * tangent, tangent, normal));
            }
            remaining -= length;
            last = Some((a + ab, tangent, normal));
        }
        last
    }

    /// Record each pinned label's placement relative to its carrier.
    fn pin_labels(&mut self) {
        for edge in 0..self.labels.len() {
            let Some(label) = self.labels[edge] else { continue };
            let (true, Some(center)) = (label.pinned, label.center.map(Vector2::from)) else {
                continue;
            };
            let Some(frame) = self.label_frame(edge, center) else { continue };
            let pin = [frame.arc / frame.length.max(1e-300), (center - frame.foot).dot(frame.normal)];
            self.labels[edge] = Some(EdgeLabel { pin: Some(pin), ..label });
        }
    }

    /// Total overlap depth of a label box at `center` with other geometry:
    /// other edges' points and segments, and other labels.
    fn label_overlap(
        &self,
        model: &Model,
        edge: usize,
        label: &EdgeLabel,
        center: Vector2<f64>,
        config: ImpredConfig,
    ) -> f64 {
        let contact = LABEL_CONTACT * config.target;
        let route = &model.route_points[edge];
        let depth = |difference: Vector2<f64>, extra: f64| {
            let distance = difference.magnitude();
            let unit = if distance > 0.0 {
                difference / distance
            } else {
                Vector2::new(1.0, 0.0)
            };
            (box_support(label.extents, unit) + extra + contact - distance).max(0.0)
        };
        let mut total = 0.0;
        for point in (0..self.positions.len()).filter(|point| !route.contains(point)) {
            total += depth(center - self.point(point), 0.0);
        }
        for segment in model.segments.iter().filter(|segment| segment.edge != edge) {
            let a = self.point(segment.a);
            let ab = self.point(segment.b) - a;
            let t = ((center - a).dot(ab) / ab.magnitude2().max(1e-20)).clamp(0.0, 1.0);
            total += depth(center - (a + t * ab), 0.0);
        }
        for (other, label) in self.labels.iter().enumerate() {
            if let Some(EdgeLabel {
                extents,
                center: Some(c),
                ..
            }) = label.filter(|_| other != edge)
            {
                let difference = center - Vector2::from(c);
                let unit = difference / difference.magnitude().max(1e-300);
                total += depth(difference, box_support(extents, unit));
            }
        }
        total
    }

    /// At refinement checkpoints, move each label to its carrier's other side
    /// when its rest there overlaps less than half as much, as the label
    /// search would choose.
    fn flip_labels(&mut self, model: &Model, config: ImpredConfig) {
        for edge in 0..self.labels.len() {
            let Some(label) = self.labels[edge].filter(|label| !label.pinned) else { continue };
            let Some(center) = label.center.map(Vector2::from) else { continue };
            let Some(frame) = self.label_frame(edge, center) else { continue };
            let side = if (center - frame.foot).dot(frame.normal) < 0.0 {
                -1.0
            } else {
                1.0
            };
            let offset = box_support(label.extents, frame.normal) + LABEL_GAP * config.target;
            let here = frame.foot + side * offset * frame.normal;
            let there = frame.foot - side * offset * frame.normal;
            if self.label_overlap(model, edge, &label, there, config)
                < 0.5 * self.label_overlap(model, edge, &label, here, config)
            {
                self.labels[edge] = Some(EdgeLabel {
                    center: Some(there.into()),
                    ..label
                });
            }
        }
    }

    /// Forces on each label, with their reactions added to `force`. A label
    /// rides a rail beside its own carrier: held at its gap across the carrier
    /// but free to slide along it within the label search's central range,
    /// so only contact that sliding cannot relieve reaches the carrier.
    /// Contact with other edges' points and segments and with other labels
    /// follows the point repulsion law in excess of its value at the reach.
    /// The tether's reaction spreads over the anchor and endpoints, so the
    /// carrier moves aside rather than bowing.
    fn label_forces(
        &self,
        model: &Model,
        config: ImpredConfig,
        rep_exponent: f64,
        force: &mut [Vector2<f64>],
    ) -> Vec<Vector2<f64>> {
        let mut result = vec![Vector2::zero(); self.labels.len()];
        let labelled: Vec<_> = self
            .labels
            .iter()
            .enumerate()
            .filter_map(|(edge, label)| {
                let label = (*label)?;
                Some((edge, label, Vector2::from(label.center?)))
            })
            .collect();
        if labelled.is_empty() {
            return result;
        }
        let contact = LABEL_CONTACT * config.target;
        let floor = config.target * 0.02;
        let excess = |reach: f64, distance: f64| {
            LABEL_STIFFNESS
                * config.repulsion
                * ((reach / distance.max(floor)).powf(rep_exponent) - 1.0).max(0.0)
                * distance.max(floor)
        };
        let direction = |difference: Vector2<f64>| {
            let distance = difference.magnitude();
            let unit = if distance > 0.0 {
                difference / distance
            } else {
                Vector2::new(1.0, 0.0)
            };
            (distance, unit)
        };
        for &(edge, label, center) in &labelled {
            let route = &self.routes[edge];
            let own = &model.route_points[edge];
            let pinned = label.pin.and_then(|[fraction, offset]| {
                let length = self.label_frame(edge, center)?.length;
                let (point, _, normal) = self.route_point_at(edge, fraction * length)?;
                Some(point + offset * normal)
            });
            if let Some(rest) = pinned {
                let tether = LABEL_STIFFNESS * 2.0 * config.attraction * (rest - center);
                result[edge] += tether;
                force[self.anchor_points[edge]] -= 0.5 * tether;
                force[route[0]] -= 0.25 * tether;
                force[*route.last().unwrap()] -= 0.25 * tether;
            } else if let Some(frame) = self.label_frame(edge, center) {
                let side = if (center - frame.foot).dot(frame.normal) < 0.0 {
                    -1.0
                } else {
                    1.0
                };
                let offset =
                    box_support(label.extents, frame.normal) + LABEL_GAP * config.target;
                let rest = frame.foot + side * offset * frame.normal;
                let across = (rest - center).dot(frame.normal) * frame.normal;
                let range = (frame.length * 3.0 / 32.0, frame.length * 29.0 / 32.0);
                let along = (frame.arc.clamp(range.0, range.1) - frame.arc) * frame.tangent;
                let tether = LABEL_STIFFNESS * 2.0 * config.attraction * (across + along);
                result[edge] += tether;
                force[self.anchor_points[edge]] -= 0.5 * tether;
                force[route[0]] -= 0.25 * tether;
                force[*route.last().unwrap()] -= 0.25 * tether;
            }
            for point in (0..self.positions.len()).filter(|point| !own.contains(point)) {
                let (distance, unit) = direction(center - self.point(point));
                let reach = box_support(label.extents, unit) + contact;
                if distance < reach {
                    let push = excess(reach, distance) * unit;
                    result[edge] += push;
                    force[point] -= push;
                }
            }
            for segment in model.segments.iter().filter(|segment| segment.edge != edge) {
                let a = self.point(segment.a);
                let ab = self.point(segment.b) - a;
                let t = ((center - a).dot(ab) / ab.magnitude2().max(1e-20)).clamp(0.0, 1.0);
                let (distance, unit) = direction(center - (a + t * ab));
                let reach = box_support(label.extents, unit) + contact;
                if distance < reach {
                    let push = excess(reach, distance) * unit;
                    result[edge] += push;
                    force[segment.a] -= (1.0 - t) * push;
                    force[segment.b] -= t * push;
                }
            }
        }
        for (i, &(first, a, ca)) in labelled.iter().enumerate() {
            for &(second, b, cb) in &labelled[i + 1..] {
                let (distance, unit) = direction(ca - cb);
                let reach = box_support(a.extents, unit) + box_support(b.extents, unit) + contact;
                if distance < reach {
                    let push = excess(reach, distance) * unit;
                    result[first] += push;
                    result[second] -= push;
                }
            }
        }
        result
    }

    fn mixed_flow(&self) -> bool {
        self.external.contains(&Some(ExternalFlow::Incoming))
            && self.external.contains(&Some(ExternalFlow::Outgoing))
    }

    fn node_center(&self) -> Vector2<f64> {
        self.node_points
            .iter()
            .map(|&point| self.point(point))
            .fold(Vector2::zero(), |a, b| a + b)
            / (self.node_points.len().max(1) as f64)
    }

    /// A rigid rotation respects only free axes and directed groups that
    /// reference their own point; pins and alignments fix the drawing frame.
    fn levelable(&self) -> bool {
        self.mixed_flow()
            && !self.nodes_fixed
            && !self.node_points.is_empty()
            && self.constraints.iter().enumerate().all(|(index, record)| {
                !record.fixed_axes.iter().any(|fixed| *fixed)
                    && record.constraints.iter().all(|axis| match axis {
                        AxisConstraint::Free => true,
                        AxisConstraint::Fixed => false,
                        AxisConstraint::Grouped { point, .. } => *point == index,
                    })
            })
    }

    /// Rotation about the physical-vertex centroid c that minimizes the
    /// horizontal pull potential -Σ F·(r - c) = -(A cos θ - B sin θ).
    fn pull_rotation(&self, pull: &[f64]) -> f64 {
        let center = self.node_center();
        let (mut along, mut across) = (0.0, 0.0);
        for ((flow, &anchor), &pull) in self.external.iter().zip(&self.anchor_points).zip(pull) {
            if let Some(flow) = flow {
                let force = flow.horizontal(pull);
                let offset = self.point(anchor) - center;
                along += force * offset.x;
                across += force * offset.y;
            }
        }
        -across.atan2(along)
    }

    /// Every force other than the pull depends only on distances, so the pull
    /// rotation is the exact optimum of the rotation mode. Halve the angle
    /// until directed coordinates keep their movement margins and the crossing
    /// witnesses still hold.
    fn level(&mut self, pull: &[f64], dofs: &MovementDofs, initial: &Geometry) {
        let angle = self.pull_rotation(pull);
        if !angle.is_finite() || angle == 0.0 {
            return;
        }
        let center = self.node_center();
        let original = self.positions.clone();
        for halving in 0..12 {
            let (sin, cos) = (angle * 0.5_f64.powi(halving)).sin_cos();
            let rotate = |point: [f64; 2]| -> [f64; 2] {
                let offset = Vector2::from(point) - center;
                (center
                    + Vector2::new(
                        cos * offset.x - sin * offset.y,
                        sin * offset.x + cos * offset.y,
                    ))
                .into()
            };
            for (point, &start) in self.positions.iter_mut().zip(&original) {
                *point = rotate(start);
            }
            if dofs.validate(&self.positions).is_ok()
                && dofs
                    .sign_rows(&self.positions)
                    .iter()
                    .all(|row| row.bound >= 0.0)
                && self.geometry().certify(initial).is_ok()
            {
                // Labels turn with the drawing they annotate.
                for label in self.labels.iter_mut().flatten() {
                    label.center = label.center.map(rotate);
                }
                return;
            }
        }
        self.positions = original;
    }

    fn protected(&self, point: usize) -> bool {
        self.node_points.contains(&point)
            || self.anchor_points.contains(&point)
            || self.constraints[point].fixed_axes.iter().any(|x| *x)
            || self.constraints[point]
                .constraints
                .iter()
                .any(|axis| !matches!(axis, AxisConstraint::Free))
            || self.constraints.iter().any(|record| {
                record.constraints.iter().any(|axis| {
                    matches!(axis, AxisConstraint::Grouped { point: reference, .. } if *reference == point)
                })
            })
    }

    fn remove_point(&mut self, point: usize) {
        self.positions.remove(point);
        self.constraints.remove(point);
        self.samples.remove(point);
        let remap = |index: &mut usize| {
            if *index > point {
                *index -= 1;
            }
        };
        for index in self.node_points.iter_mut().chain(&mut self.anchor_points) {
            remap(index);
        }
        for route in &mut self.routes {
            for index in route {
                remap(index);
            }
        }
        for record in &mut self.constraints {
            for axis in &mut record.constraints {
                if let AxisConstraint::Grouped {
                    point: reference, ..
                } = axis
                {
                    remap(reference);
                }
            }
        }
    }

    fn empty_corner(&self, corner: [usize; 3], target: f64) -> bool {
        let [a, b, c] = corner.map(|point| self.point(point));
        let low = Vector2::new(a.x.min(b.x).min(c.x), a.y.min(b.y).min(c.y));
        let high = Vector2::new(a.x.max(b.x).max(c.x), a.y.max(b.y).max(c.y));
        let tolerance = 1e-10 * target.max(1.0);
        for point in 0..self.positions.len() {
            if corner.contains(&point) {
                continue;
            }
            let p = self.point(point);
            if (0..2)
                .any(|axis| p[axis] < low[axis] - tolerance || p[axis] > high[axis] + tolerance)
            {
                continue;
            }
            let orientations = [
                (b - a).perp_dot(p - a),
                (c - b).perp_dot(p - b),
                (a - c).perp_dot(p - c),
            ];
            if orientations.iter().all(|value| *value >= -tolerance)
                || orientations.iter().all(|value| *value <= tolerance)
            {
                return false;
            }
        }
        true
    }

    fn contract_external(
        &mut self,
        config: ImpredConfig,
        eligible: &[bool],
        initial: &Geometry,
    ) -> usize {
        let mut removed = 0;
        for (edge, &eligible) in eligible.iter().enumerate() {
            if !eligible {
                continue;
            }
            let mut slot = 1;
            while slot + 1 < self.routes[edge].len() {
                let corner = [
                    self.routes[edge][slot - 1],
                    self.routes[edge][slot],
                    self.routes[edge][slot + 1],
                ];
                let chord = (self.point(corner[2]) - self.point(corner[0])).magnitude();
                if !self.protected(corner[1])
                    && chord > 1e-7
                    && chord < config.contract_chord_ratio * config.target
                    && self.empty_corner(corner, config.target)
                {
                    self.routes[edge].remove(slot);
                    if self.geometry().certify(initial).is_ok() {
                        self.remove_point(corner[1]);
                        removed += 1;
                        continue;
                    }
                    self.routes[edge].insert(slot, corner[1]);
                }
                slot += 1;
            }
        }
        removed
    }

    fn insert_external(&mut self, config: ImpredConfig, eligible: &[bool]) -> usize {
        let mut inserted = 0;
        for (edge, &eligible) in eligible.iter().enumerate() {
            if !eligible {
                continue;
            }
            let old = self.routes[edge].clone();
            let mut room = config.external_max_points.saturating_sub(old.len() - 2);
            let mut route = vec![old[0]];
            for pair in old.windows(2) {
                if room > 0
                    && (self.point(pair[1]) - self.point(pair[0])).magnitude()
                        > config.target * config.split_length_ratio
                {
                    let midpoint: [f64; 2] =
                        ((self.point(pair[0]) + self.point(pair[1])) * 0.5).into();
                    let anchor = &self.constraints[self.anchor_points[edge]];
                    let fixed_axes = std::array::from_fn(|axis| {
                        anchor.fixed_axes[axis]
                            || matches!(anchor.constraints[axis], AxisConstraint::Fixed)
                    });
                    route.push(self.positions.len());
                    self.positions.push(midpoint);
                    self.samples.push(true);
                    self.constraints.push(PointRecord {
                        constraints: fixed_axes.map(|fixed| {
                            if fixed {
                                AxisConstraint::Fixed
                            } else {
                                AxisConstraint::Free
                            }
                        }),
                        shift: [0.0, 0.0],
                        reference: midpoint,
                        fixed_axes,
                        is_node: false,
                    });
                    inserted += 1;
                    room -= 1;
                }
                route.push(pair[1]);
            }
            self.routes[edge] = route;
        }
        inserted
    }

    /// Relax the supplied embedding without changing its crossing witnesses.
    /// Native affine degrees of freedom are eliminated before every protected
    /// movement solve. Existing crossing routes keep their fixed sample counts.
    pub fn solve(&mut self, config: ImpredConfig) -> Result<ImpredReport, String> {
        self.solve_with_progress(config, |_, _| Ok(()))
    }

    /// The callback observes the same hundred-step cooling snapshots as the
    /// reference solver, reports nominal progress, and can cancel with an error.
    pub fn solve_with_progress(
        &mut self,
        config: ImpredConfig,
        mut progress: impl FnMut(usize, &Self) -> Result<(), String>,
    ) -> Result<ImpredReport, String> {
        config.validate()?;
        self.validate()?;
        let initial = self.geometry();
        initial.certify(&initial)?;
        let mut report = ImpredReport {
            initial_crossings: initial.crossings.len(),
            ..Default::default()
        };
        let crossed: BTreeSet<_> = initial
            .crossings
            .iter()
            .flat_map(|&(a, b)| [a.0, b.0])
            .collect();
        let eligible: Vec<_> = self
            .external
            .iter()
            .enumerate()
            .map(|(edge, flow)| flow.is_some() && !crossed.contains(&edge))
            .collect();
        let flexible: Vec<_> = self
            .routes
            .iter()
            .enumerate()
            .map(|(edge, route)| {
                !crossed.contains(&edge)
                    && ((self.external[edge].is_none() && route.len() > 2)
                        || (self.external[edge].is_some() && config.external_max_points > 0))
            })
            .collect();
        let pull: Vec<_> = self
            .external
            .iter()
            .enumerate()
            .map(|(edge, flow)| {
                if flow.is_none() || config.pull == 0.0 {
                    0.0
                } else {
                    config.pull * config.target * self.pull_scales[edge].powf(config.pull_balance)
                }
            })
            .collect();
        if pull.iter().any(|value| !value.is_finite()) {
            return Err("Topology balancing produces non-finite external force".into());
        }
        let mut dofs = MovementDofs::new(&self.constraints, self.nodes_fixed)?;
        dofs.validate(&self.positions)?;
        if !config.labels {
            self.labels.clear();
        }
        self.seat_labels(config);
        self.pin_labels();
        let mut model = self.model(&flexible, config.parallel_attraction_balance);
        if model.segments.is_empty() {
            return Ok(report);
        }
        let cap_normals: [[f64; 2]; 16] = std::array::from_fn(|i| {
            let angle = (i as f64 + 0.5) * (2.0 * std::f64::consts::PI / 16.0);
            [angle.cos(), angle.sin()]
        });
        let cap_apothem = (std::f64::consts::PI / 16.0).cos();
        let mut epoch = 0;
        while epoch < config.steps {
            let stride = config.integration_stride(epoch);
            let u = config.warm_start
                + (1.0 - config.warm_start) * epoch as f64 / config.steps.saturating_sub(1).max(1) as f64;
            let pairs = self.point_segments(&model);
            let mut forces = self.forces(&model, &pairs, config, u, &pull);
            let label_forces = self.label_forces(&model, config, 2.0 + 2.0 * u, &mut forces);
            let proposal: Vec<[f64; 2]> = forces
                .into_iter()
                .map(|force| (force * (0.045 * stride as f64)).into())
                .collect();
            let cap = config.target * (0.10 * (1.0 - u) + 0.002) * stride as f64;
            let separators = self.separators(&model, &pairs, cap);
            let mut rows: Vec<_> = separators
                .iter()
                .filter(|row| row.bound < cap * (1.0 + 1e-12))
                .cloned()
                .collect();
            rows.reserve(cap_normals.len() * self.positions.len());
            for point in 0..self.positions.len() {
                for &normal in &cap_normals {
                    rows.push(Halfplane {
                        point,
                        normal,
                        bound: cap * cap_apothem,
                    });
                }
            }
            rows.extend(dofs.sign_rows(&self.positions));
            let (step, certificate) = dofs.project(&proposal, &rows, cap)?;
            let separator_error = separators
                .iter()
                .map(|row| Vector2::from(row.normal).dot(step[row.point].into()) - row.bound)
                .fold(0.0, f64::max);
            let cap_error = step
                .iter()
                .map(|&point| Vector2::from(point).magnitude() - cap)
                .fold(0.0, f64::max);
            // The Euclidean check is tighter for ordinary magnitudes. This
            // component bound also certifies omitted far separators when squared
            // lengths underflow, without relying on unit-normal roundoff bounds.
            let component_error = step.iter().flatten().any(|value| value.abs() > 2.0 * cap);
            if separator_error > 2e-9 * cap || cap_error > 2e-9 * cap || component_error {
                return Err(format!(
                    "ImPrEd movement geometry certificate failed: separators={separator_error},cap={cap_error},component_limit_exceeded={component_error}"
                ));
            }
            report.max_primal_residual =
                report.max_primal_residual.max(certificate.primal_residual);
            report.max_stationarity_residual = report
                .max_stationarity_residual
                .max(certificate.stationarity_residual);
            report.max_complementarity_residual = report
                .max_complementarity_residual
                .max(certificate.complementarity_residual);
            report.max_separator_residual = report.max_separator_residual.max(separator_error);
            report.max_cap_residual = report.max_cap_residual.max(cap_error);
            for (point, delta) in self.positions.iter_mut().zip(step) {
                point[0] += delta[0];
                point[1] += delta[1];
            }
            // Labels move freely under the same cap; they certify nothing.
            for (label, force) in self.labels.iter_mut().zip(label_forces) {
                if let Some(center) = label.as_mut().and_then(|label| label.center.as_mut()) {
                    let delta = force * (0.045 * stride as f64);
                    let delta = delta * (cap / delta.magnitude().max(cap));
                    center[0] += delta.x;
                    center[1] += delta.y;
                }
            }
            dofs.validate(&self.positions)?;
            report.iterations += 1;
            if epoch % 25 == 0 || epoch + 1 == config.steps {
                // Refinement can pin new samples, so eligibility is rechecked.
                if config.level && self.levelable() {
                    self.level(&pull, &dofs, &initial);
                }
                self.flip_labels(&model, config);
                self.geometry().certify(&initial)?;
                report.geometry_checks += 1;
                if epoch % 100 == 0 || epoch + 1 == config.steps {
                    progress(epoch + 1, self)?;
                }
                if config.external_max_points > 0 && epoch > 0 {
                    let removed = self.contract_external(config, &eligible, &initial);
                    let inserted = self.insert_external(config, &eligible);
                    if removed > 0 || inserted > 0 {
                        self.geometry().certify(&initial)?;
                        report.removed_points += removed;
                        report.inserted_points += inserted;
                        dofs = MovementDofs::new(&self.constraints, self.nodes_fixed)?;
                        dofs.validate(&self.positions)?;
                        model = self.model(&flexible, config.parallel_attraction_balance);
                    }
                }
            }
            epoch += stride;
        }
        Ok(report)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FIXTURES: [(&str, &str); 4] = [
        (
            "higgs-sunset",
            include_str!("impred/fixtures/higgs-sunset.json"),
        ),
        (
            "xs-scattering",
            include_str!("impred/fixtures/xs-extra-higgs-scattering.json"),
        ),
        (
            "eight-gluon-ladder",
            include_str!("impred/fixtures/multi-leg-eight-gluon-ladder.json"),
        ),
        (
            "nonplanar-k5",
            include_str!("impred/fixtures/nonplanar-gluon-k5.json"),
        ),
    ];

    #[derive(Deserialize)]
    struct ForceSample {
        progress: f64,
        forces: Vec<[f64; 2]>,
    }
    #[derive(Deserialize)]
    struct Fixture {
        layout: ImpredLayout,
        config: ImpredConfig,
        forces: Vec<ForceSample>,
        separators: Vec<Halfplane>,
        final_state: FinalState,
    }

    #[derive(Deserialize)]
    struct FinalState {
        routes: Vec<Vec<[f64; 2]>>,
        insertions: usize,
        removals: usize,
    }

    #[test]
    fn native_forces_and_separators_match_frozen_python_reference() {
        for (name, source) in FIXTURES {
            let fixture: Fixture = serde_json::from_str(source).unwrap();
            let layout = fixture.layout;
            layout.validate().unwrap();
            let crossed: BTreeSet<_> = layout
                .geometry()
                .crossings
                .iter()
                .flat_map(|&(a, b)| [a.0, b.0])
                .collect();
            let flexible: Vec<_> = layout
                .routes
                .iter()
                .enumerate()
                .map(|(edge, route)| {
                    !crossed.contains(&edge) && (layout.external[edge].is_some() || route.len() > 2)
                })
                .collect();
            let model = layout.model(&flexible, fixture.config.parallel_attraction_balance);
            let pairs = layout.point_segments(&model);
            let pull: Vec<_> = layout
                .pull_scales
                .iter()
                .map(|value| {
                    fixture.config.pull
                        * fixture.config.target
                        * value.powf(fixture.config.pull_balance)
                })
                .collect();
            for sample in fixture.forces {
                let actual = layout.forces(&model, &pairs, fixture.config, sample.progress, &pull);
                for (point, (actual, expected)) in actual.iter().zip(&sample.forces).enumerate() {
                    let expected = Vector2::from(*expected);
                    let tolerance = 1e-11 * expected.magnitude().max(1.0);
                    assert!(
                        (*actual - expected).magnitude() <= tolerance,
                        "{name} point {point} at {}: {actual:?} != {expected:?}",
                        sample.progress
                    );
                }
            }
            let actual = layout.separators(&model, &pairs, f64::INFINITY);
            assert_eq!(actual.len(), fixture.separators.len(), "{name}");
            let mut used = vec![false; fixture.separators.len()];
            for row in actual {
                let found = fixture
                    .separators
                    .iter()
                    .enumerate()
                    .position(|(i, expected)| {
                        !used[i]
                            && row.point == expected.point
                            && (Vector2::from(row.normal) - Vector2::from(expected.normal))
                                .magnitude()
                                < 1e-12
                            && (row.bound - expected.bound).abs() < 1e-12
                    });
                let Some(index) = found else {
                    panic!("{name} separator missing from reference: {row:?}");
                };
                used[index] = true;
            }
        }
    }

    #[test]
    fn infinite_diagnostic_cap_retains_even_infinite_bound_rows() {
        let layout = dangling(&[[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [1e200, 0.0]]);
        let model = layout.model(&[true], 1.0);
        let pairs = layout.point_segments(&model);
        let rows = layout.separators(&model, &pairs, f64::INFINITY);
        assert!(rows.iter().any(|row| row.bound.is_infinite()));
        assert_eq!(
            rows.len(),
            3 * pairs.len() + layout.positions.len() * (layout.positions.len() - 1)
        );
    }

    #[test]
    fn omitted_separators_are_satisfied_at_every_component_cap_corner() {
        for (name, source) in FIXTURES {
            let fixture: Fixture = serde_json::from_str(source).unwrap();
            let layout = fixture.layout;
            let model = layout.model(&vec![false; layout.routes.len()], 1.0);
            let pairs = layout.point_segments(&model);
            let full = layout.separators(&model, &pairs, f64::INFINITY);
            let key =
                |row: &Halfplane| (row.point, row.normal.map(f64::to_bits), row.bound.to_bits());
            let mut omitted = 0;
            for cap in [f64::from_bits(1), f64::MIN_POSITIVE, 0.001, 0.1, 1.0] {
                let retained: BTreeSet<_> = layout
                    .separators(&model, &pairs, cap)
                    .iter()
                    .map(key)
                    .collect();
                for row in full.iter().filter(|row| !retained.contains(&key(row))) {
                    omitted += 1;
                    // A linear separator attains its box maximum at a corner.
                    // Test the actual floating-point evaluation, including caps
                    // so small that the Euclidean norm would underflow.
                    for x in [-2.0 * cap, 2.0 * cap] {
                        for y in [-2.0 * cap, 2.0 * cap] {
                            assert!(
                                row.normal[0] * x + row.normal[1] * y <= row.bound,
                                "{name}: omitted {row:?} rejects the certified corner {x},{y}"
                            );
                        }
                    }
                }
            }
            assert!(omitted > 0, "{name}: fixture must exercise far-row pruning");
        }
    }

    #[test]
    fn complete_relaxation_matches_the_accepted_python_layouts() {
        for (name, source) in FIXTURES {
            let mut fixture: Fixture = serde_json::from_str(source).unwrap();
            let report = fixture.layout.solve(fixture.config).unwrap();
            assert_eq!(report.iterations, fixture.config.steps, "{name}");
            assert_eq!(
                report.inserted_points, fixture.final_state.insertions,
                "{name}"
            );
            assert_eq!(
                report.removed_points, fixture.final_state.removals,
                "{name}"
            );
            assert_eq!(
                fixture.layout.routes.len(),
                fixture.final_state.routes.len(),
                "{name}"
            );
            // Near active separator intersections, roundoff can amplify along
            // the nonplanar trajectory while all movement certificates hold.
            let tolerance = if name == "nonplanar-k5" { 1e-6 } else { 1e-9 };
            for (edge, (route, expected)) in fixture
                .layout
                .routes
                .iter()
                .zip(fixture.final_state.routes)
                .enumerate()
            {
                assert_eq!(route.len(), expected.len(), "{name} edge {edge}");
                for (&point, expected) in route.iter().zip(expected) {
                    let distance =
                        (fixture.layout.point(point) - Vector2::from(expected)).magnitude();
                    assert!(
                        distance <= tolerance,
                        "{name} edge {edge}: coordinate difference {distance}"
                    );
                }
            }
        }
    }

    #[test]
    fn larger_steps_preserve_routes_and_cooling_checkpoints() {
        for (name, source) in FIXTURES {
            let fixture: Fixture = serde_json::from_str(source).unwrap();
            let mut reference = fixture.layout.clone();
            let mut accelerated = fixture.layout;
            let mut reference_progress = Vec::new();
            let reference_report = reference
                .solve_with_progress(fixture.config, |epoch, _| {
                    reference_progress.push(epoch);
                    Ok(())
                })
                .unwrap();
            let mut accelerated_progress = Vec::new();
            let report = accelerated
                .solve_with_progress(
                    ImpredConfig {
                        step_scale: 2,
                        ..fixture.config
                    },
                    |epoch, _| {
                        accelerated_progress.push(epoch);
                        Ok(())
                    },
                )
                .unwrap();
            assert!(
                report.iterations < 3 * reference_report.iterations / 4,
                "{name}"
            );
            assert_eq!(reference_progress, accelerated_progress, "{name}");
            assert_eq!(
                reference_report.geometry_checks, report.geometry_checks,
                "{name}"
            );
            assert_eq!(reference.routes, accelerated.routes, "{name}");
            assert_eq!(reference.geometry(), accelerated.geometry(), "{name}");

            // Translation does not change the rendered shape. Retain a tight
            // physical-distance bound on every vertex and route point, without
            // permitting rotation, scaling, or a change of route correspondence.
            let translation = reference
                .node_points
                .iter()
                .fold(Vector2::zero(), |sum, &p| {
                    sum + accelerated.point(p) - reference.point(p)
                })
                / reference.node_points.len() as f64;
            for point in 0..reference.positions.len() {
                let error =
                    (accelerated.point(point) - translation - reference.point(point)).magnitude();
                assert!(
                    error < 0.01 * fixture.config.target,
                    "{name}: point {point} moved {error}"
                );
            }
        }
    }

    #[test]
    fn levelling_reaches_the_external_pull_optimum_within_constraints() {
        let fixture: Fixture = serde_json::from_str(FIXTURES[3].1).unwrap();
        let config = ImpredConfig {
            steps: 200,
            ..fixture.config
        };
        let pull: Vec<_> = fixture
            .layout
            .pull_scales
            .iter()
            .map(|scale| config.pull * config.target * scale.powf(config.pull_balance))
            .collect();
        // Tilt the drawing together with its directed references. Overdamped
        // forces alone do not turn it back within a short schedule.
        let (sin, cos) = 25_f64.to_radians().sin_cos();
        let turn = |[x, y]: [f64; 2]| [cos * x - sin * y, sin * x + cos * y];
        let mut tilted = fixture.layout;
        for point in &mut tilted.positions {
            *point = turn(*point);
        }
        for record in &mut tilted.constraints {
            record.reference = turn(record.reference);
        }
        assert!(tilted.levelable());

        let mut drifting = tilted.clone();
        drifting
            .solve(ImpredConfig {
                level: false,
                ..config
            })
            .unwrap();
        let residual = drifting.pull_rotation(&pull);
        assert!(residual.abs() > 1_f64.to_radians(), "residual {residual}");

        let mut levelled = tilted.clone();
        levelled
            .solve(ImpredConfig {
                level: true,
                ..config
            })
            .unwrap();
        let residual = levelled.pull_rotation(&pull);
        assert!(residual.abs() < 1e-9, "residual {residual}");

        // Shared-coordinate external groups define the drawing frame.
        let aligned: Fixture = serde_json::from_str(FIXTURES[2].1).unwrap();
        assert!(!aligned.layout.levelable());
    }

    #[test]
    fn integration_stride_handles_small_budgets_and_rejects_zero_scale() {
        let fixture: Fixture = serde_json::from_str(FIXTURES[0].1).unwrap();
        let accelerated = ImpredConfig {
            step_scale: 2,
            ..fixture.config
        };
        assert_eq!(accelerated.integration_stride(750), 1);
        assert_eq!(accelerated.integration_stride(751), 2);
        assert_eq!(accelerated.integration_stride(774), 1);
        assert_eq!(accelerated.integration_stride(775), 1);
        for steps in [0, 1, 2, 24, 25, 26, 100, 101] {
            let mut layout = fixture.layout.clone();
            let mut progress = Vec::new();
            let report = layout
                .solve_with_progress(
                    ImpredConfig {
                        steps,
                        step_scale: usize::MAX,
                        ..fixture.config
                    },
                    |epoch, _| {
                        progress.push(epoch);
                        Ok(())
                    },
                )
                .unwrap();
            assert!(report.iterations <= steps);
            if steps == 0 {
                assert!(progress.is_empty());
            } else {
                assert_eq!(progress.last(), Some(&steps));
            }
        }
        assert!(
            ImpredConfig {
                step_scale: 0,
                ..fixture.config
            }
            .validate()
            .is_err()
        );
    }

    fn dangling(points: &[[f64; 2]]) -> ImpredLayout {
        ImpredLayout {
            positions: points.to_vec(),
            routes: vec![(0..points.len()).collect()],
            node_points: vec![0],
            anchor_points: vec![points.len() - 1],
            external: vec![Some(ExternalFlow::Outgoing)],
            constraints: points
                .iter()
                .map(|&reference| PointRecord {
                    constraints: [AxisConstraint::Free, AxisConstraint::Free],
                    shift: [0.0, 0.0],
                    reference,
                    fixed_axes: [false, false],
                    is_node: false,
                })
                .collect(),
            pull_scales: vec![1.0],
            samples: (0..points.len())
                .map(|i| i > 0 && i + 1 < points.len())
                .collect(),
            nodes_fixed: false,
            labels: Vec::new(),
        }
    }

    #[test]
    fn subdivision_retains_charge_budget_and_serial_compliance() {
        let config = ImpredConfig::default();
        let mut layout = dangling(&[[0.0, 0.0], [8.0, 0.0]]);
        let before = layout.segment_forces(&layout.model(&[true], 1.0), config, 1.0);
        assert_eq!(layout.insert_external(config, &[true]), 1);
        assert_eq!(layout.routes[0], vec![0, 2, 1]);
        assert_eq!(layout.positions[2], [4.0, 0.0]);
        let model = layout.model(&[true], 1.0);
        assert!(
            model
                .pair_weights
                .iter()
                .flatten()
                .all(|&value| value == 0.0)
        );
        let after = layout.segment_forces(&model, config, 1.0);
        assert_eq!(after[0], before[0] * 0.25);
        assert_eq!(after[2], Vector2::zero());
    }

    #[test]
    fn contraction_cannot_sweep_across_an_isolated_physical_vertex() {
        let config = ImpredConfig::default();
        let mut layout = dangling(&[[0.0, 0.0], [1.0, 1.0], [2.0, 0.0]]);
        let isolated = layout.positions.len();
        layout.positions.push([1.0, 0.25]);
        layout.samples.push(false);
        layout.node_points.push(isolated);
        layout.constraints.push(PointRecord {
            constraints: [AxisConstraint::Free, AxisConstraint::Free],
            shift: [0.0, 0.0],
            reference: [1.0, 0.25],
            fixed_axes: [false, false],
            is_node: true,
        });
        let geometry = layout.geometry();
        assert_eq!(layout.contract_external(config, &[true], &geometry), 0);
        assert_eq!(layout.routes[0], vec![0, 1, 2]);
        layout.positions[isolated] = [1.0, 2.0];
        assert_eq!(layout.contract_external(config, &[true], &geometry), 1);
        assert_eq!(layout.routes[0], vec![0, 1]);
        assert_eq!(layout.node_points, vec![0, 2]);
        assert_eq!(layout.anchor_points, vec![1]);
    }
}
