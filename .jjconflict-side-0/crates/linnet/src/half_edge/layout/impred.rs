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
    pub steps: usize,
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
}

impl Default for ImpredConfig {
    fn default() -> Self {
        Self {
            steps: 1500,
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
        }
    }
}

impl ImpredConfig {
    fn validate(self) -> Result<(), String> {
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
        if self.external_max_points > 3 {
            return Err("External routes support at most three intermediate points".into());
        }
        Ok(())
    }
}

/// External momentum flow. Mixed flow uses horizontal pull; a single flow uses
/// radial pull from the centroid of the physical interaction nodes.
#[derive(Clone, Copy, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub enum ExternalFlow {
    Incoming,
    Outgoing,
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
}

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct ImpredReport {
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
    normal: Vector2<f64>,
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
        }
    }

    fn point_segments(&self, model: &Model) -> Vec<PointSegment> {
        let mut result =
            Vec::with_capacity(self.positions.len().saturating_sub(2) * model.segments.len());
        for point in 0..self.positions.len() {
            for (index, segment) in model.segments.iter().enumerate() {
                if point == segment.a || point == segment.b {
                    continue;
                }
                let a = self.point(segment.a);
                let ab = self.point(segment.b) - a;
                let t =
                    ((self.point(point) - a).dot(ab) / ab.magnitude2().max(1e-20)).clamp(0.0, 1.0);
                let difference = self.point(point) - (a + t * ab);
                let distance = difference.magnitude();
                result.push(PointSegment {
                    point,
                    segment: index,
                    t,
                    distance,
                    normal: difference / distance.max(1e-15),
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
        for i in 0..self.positions.len() {
            for j in i..self.positions.len() {
                let difference = self.point(i) - self.point(j);
                let coefficient = model.pair_weights[i][j]
                    * (config.target / difference.magnitude().max(config.target * 0.02))
                        .powf(rep_exponent);
                force[i] += config.repulsion * coefficient * difference;
                if i != j {
                    // Recompute the reversed difference rather than negating:
                    // equal coordinates can carry different signed zeros.
                    force[j] += config.repulsion * coefficient * (self.point(j) - self.point(i));
                }
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
                    * pair.normal;
                let segment = model.segments[pair.segment];
                force[pair.point] += repulsion;
                force[segment.a] -= (1.0 - pair.t) * repulsion;
                force[segment.b] -= pair.t * repulsion;
            }
        }
        let mixed = self.external.contains(&Some(ExternalFlow::Incoming))
            && self.external.contains(&Some(ExternalFlow::Outgoing));
        let center = self
            .node_points
            .iter()
            .map(|&point| self.point(point))
            .fold(Vector2::zero(), |a, b| a + b)
            / (self.node_points.len().max(1) as f64);
        let mut reaction = Vector2::zero();
        for (edge, flow) in self.external.iter().enumerate() {
            if let Some(flow) = flow {
                let anchor = self.anchor_points[edge];
                let outward = if mixed {
                    Vector2::new(
                        if *flow == ExternalFlow::Incoming {
                            -pull[edge]
                        } else {
                            pull[edge]
                        },
                        0.0,
                    )
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

    fn separators(&self, model: &Model, pairs: &[PointSegment]) -> Vec<Halfplane> {
        let mut lengths: Vec<_> = model
            .segments
            .iter()
            .map(|s| (self.point(s.b) - self.point(s.a)).magnitude())
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
            let vector = self.point(segment.b) - self.point(segment.a);
            let perpendicular = Vector2::new(-vector.y, vector.x) / vector.magnitude().max(1e-30);
            let signed = (self.point(pair.point) - self.point(segment.a)).dot(perpendicular);
            let interior = pair.t > 0.0 && pair.t < 1.0;
            let normal = if interior {
                perpendicular * if signed >= 0.0 { 1.0 } else { -1.0 }
            } else {
                pair.normal
            };
            let distance = if interior {
                signed.abs()
            } else {
                pair.distance
            };
            let room = 0.45 * (distance - clearance).max(0.0);
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
                    rows.push(Halfplane {
                        point: i,
                        normal: (-difference / distance.max(1e-30)).into(),
                        bound: 0.45 * (distance - clearance).max(0.0),
                    });
                }
            }
        }
        rows
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

    /// The callback observes the same hundred-iteration snapshots as the
    /// reference solver and can cancel by returning an error.
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
        let mut model = self.model(&flexible, config.parallel_attraction_balance);
        if model.segments.is_empty() {
            return Ok(report);
        }
        let cap_normals: [[f64; 2]; 16] = std::array::from_fn(|i| {
            let angle = (i as f64 + 0.5) * (2.0 * std::f64::consts::PI / 16.0);
            [angle.cos(), angle.sin()]
        });
        let cap_apothem = (std::f64::consts::PI / 16.0).cos();
        for epoch in 0..config.steps {
            let u = epoch as f64 / config.steps.saturating_sub(1).max(1) as f64;
            let pairs = self.point_segments(&model);
            let proposal: Vec<[f64; 2]> = self
                .forces(&model, &pairs, config, u, &pull)
                .into_iter()
                .map(|force| (force * 0.045).into())
                .collect();
            let cap = config.target * (0.10 * (1.0 - u) + 0.002);
            let separators = self.separators(&model, &pairs);
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
            if separator_error > 2e-9 * cap || cap_error > 2e-9 * cap {
                return Err(format!(
                    "ImPrEd movement geometry certificate failed: separators={separator_error},cap={cap_error}"
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
            dofs.validate(&self.positions)?;
            report.iterations = epoch + 1;
            if epoch % 25 == 0 || epoch + 1 == config.steps {
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
            let actual = layout.separators(&model, &pairs);
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
