use std::collections::HashMap;

use cgmath::{EuclideanSpace, MetricSpace, Point2, Vector2, Zero};

use rand::{distributions::Uniform, prelude::Distribution, Rng};
use serde::{Deserialize, Serialize};

use crate::{
    half_edge::{
        involution::{EdgeIndex, EdgeVec, Flow, Hedge, HedgePair, HedgeVec},
        layout::simulatedanneale::{Energy, Neighbor},
        nodestore::NodeStorageOps,
        subgraph::{subset::SubSet, Inclusion, ModifySubSet, SuBitGraph, SubSetLike},
        swap::Swap,
        HedgeGraph, NodeIndex, NodeVec,
    },
    parser::GlobalData,
};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
#[cfg_attr(feature = "rkyv", archive(check_bytes))]
pub struct PointConstraint {
    pub x: Constraint,
    pub y: Constraint,
}

impl Default for PointConstraint {
    fn default() -> Self {
        PointConstraint {
            x: Constraint::Free,
            y: Constraint::Free,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
#[cfg_attr(feature = "rkyv", archive(check_bytes))]
pub enum ShiftDirection {
    Any,
    PositiveOnly,
    NegativeOnly,
}

/// A node or edge-control-point index in the shared layout coordinate space.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Hash)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
#[cfg_attr(feature = "rkyv", archive(check_bytes))]
pub enum LayoutPointIndex {
    Node(NodeIndex),
    Edge(EdgeIndex),
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
#[cfg_attr(feature = "rkyv", archive(check_bytes))]
pub enum Constraint {
    Fixed,
    #[default]
    Free,
    Grouped(LayoutPointIndex, ShiftDirection),
}

impl Constraint {
    pub(crate) fn force_target(self, index: LayoutPointIndex) -> Option<LayoutPointIndex> {
        match self {
            Constraint::Fixed => None,
            Constraint::Free => Some(index),
            Constraint::Grouped(reference, _) => Some(reference),
        }
    }
}

pub trait HasPointConstraint {
    fn point_constraint(&self) -> &PointConstraint;
}

impl HasPointConstraint for PointConstraint {
    fn point_constraint(&self) -> &PointConstraint {
        self
    }
}

pub(crate) fn directional_force_shift(
    constraints: &PointConstraint,
    index: LayoutPointIndex,
    point: Point2<f64>,
    magnitude: f64,
) -> Vector2<f64> {
    if magnitude == 0.0 {
        return Vector2::zero();
    }

    let component = |constraint, coordinate| match constraint {
        Constraint::Grouped(reference, ShiftDirection::PositiveOnly)
            if reference == index && coordinate <= 0.0 =>
        {
            magnitude
        }
        Constraint::Grouped(reference, ShiftDirection::NegativeOnly)
            if reference == index && coordinate >= 0.0 =>
        {
            -magnitude
        }
        _ => 0.0,
    };

    Vector2::new(
        component(constraints.x, point.x),
        component(constraints.y, point.y),
    )
}

/// A physical edge's normalized charge at its anchor or a half-edge route point.
pub(super) type EdgeChargeSample = (Option<(Hedge, usize)>, Point2<f64>, f64);

/// Geometry-independent load and extra demand from distributed shared-X attachments.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ExternalPullTopology {
    pub load_share: f64,
    pub attachment_deficit: f64,
}

impl Default for ExternalPullTopology {
    fn default() -> Self {
        Self {
            load_share: 1.0,
            attachment_deficit: 0.0,
        }
    }
}

pub struct LayoutState<'a, E, V, H, N: NodeStorageOps<NodeData = V>> {
    pub graph: &'a HedgeGraph<E, V, H, N>,
    pub ext: SuBitGraph,
    active_subgraph: Option<SuBitGraph>,
    /// External load sharing and shared-axis attachment demand, independent of geometry.
    pub(super) external_pull_topology: EdgeVec<ExternalPullTopology>,
    pub vertex_points: NodeVec<Point2<f64>>,
    pub edge_points: EdgeVec<Point2<f64>>,
    /// Interior drawing points, ordered from each half-edge's vertex to its edge anchor.
    /// These carry no interaction-node identity and do not change the graph topology.
    pub route_points: HedgeVec<Vec<Point2<f64>>>,
    /// Raw auxiliary depths: optional seeds before force layout, final raw values afterwards.
    /// Effective depth collapses to zero without changing hard-pinned raw values.
    pub vertex_depths: NodeVec<Option<f64>>,
    pub edge_depths: EdgeVec<Option<f64>>,
    /// Independent of XY constraints; callers must also pin depths outside an active subgraph.
    pub vertex_depth_pins: NodeVec<bool>,
    pub edge_depth_pins: EdgeVec<bool>,
    /// Dimensionless rest-length multipliers, defaulting to one; set before starting the layout.
    pub edge_spring_length_scales: EdgeVec<f64>,
    pub delta: f64,
    pub directional_force: f64,
    // Tracks which node/edge entries were mutated during proposal generation so
    // the energy function can update only the affected terms.
    pub changed_nodes: SubSet<NodeIndex>,
    pub changed_edges: SubSet<EdgeIndex>,
    pub incremental: bool,
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    pub fn new_layout_state(
        &self,
        vertex_points: NodeVec<Point2<f64>>,
        edge_points: EdgeVec<Point2<f64>>,
        delta: f64,
        directional_force: f64,
        incremental: bool,
    ) -> LayoutState<'_, E, V, H, N>
    where
        E: HasPointConstraint,
        V: HasPointConstraint,
    {
        let ext = self.external_filter();
        let len_v = vertex_points.len().0;
        let len_e = edge_points.len().0;
        let mut state = LayoutState {
            graph: self,
            ext,
            active_subgraph: None,
            external_pull_topology: vec![ExternalPullTopology::default(); len_e].into(),
            vertex_points,
            edge_points,
            route_points: vec![Vec::new(); self.n_hedges()].into(),
            vertex_depths: vec![None; len_v].into(),
            edge_depths: vec![None; len_e].into(),
            vertex_depth_pins: vec![false; len_v].into(),
            edge_depth_pins: vec![false; len_e].into(),
            edge_spring_length_scales: vec![1.0; len_e].into(),
            delta,
            directional_force,
            changed_nodes: SubSet::empty(len_v),
            changed_edges: SubSet::empty(len_e),
            incremental,
        };
        state.rebuild_external_pull_topology();
        state
    }
}

impl<'a, E, V, H, N: NodeStorageOps<NodeData = V>> Clone for LayoutState<'a, E, V, H, N> {
    fn clone(&self) -> Self {
        LayoutState {
            graph: self.graph,
            ext: self.ext.clone(),
            active_subgraph: self.active_subgraph.clone(),
            external_pull_topology: self.external_pull_topology.clone(),
            vertex_points: self.vertex_points.clone(),
            edge_points: self.edge_points.clone(),
            route_points: self.route_points.clone(),
            vertex_depths: self.vertex_depths.clone(),
            edge_depths: self.edge_depths.clone(),
            vertex_depth_pins: self.vertex_depth_pins.clone(),
            edge_depth_pins: self.edge_depth_pins.clone(),
            edge_spring_length_scales: self.edge_spring_length_scales.clone(),
            delta: self.delta,
            directional_force: self.directional_force,
            changed_nodes: self.changed_nodes.clone(),
            changed_edges: self.changed_edges.clone(),
            incremental: self.incremental,
        }
    }
}

impl<'a, E, V, H, N: NodeStorageOps<NodeData = V>> LayoutState<'a, E, V, H, N> {
    pub fn with_active_subgraph(mut self, active_subgraph: SuBitGraph) -> Self
    where
        E: HasPointConstraint,
        V: HasPointConstraint,
    {
        assert_eq!(active_subgraph.size(), self.graph.n_hedges());
        let mut ext: SuBitGraph = self.graph.empty_subgraph();
        for (pair, _, _) in self.graph.iter_edges_of(&active_subgraph) {
            let hedge = match pair {
                HedgePair::Paired { .. } => None,
                HedgePair::Split {
                    source,
                    sink,
                    split,
                } => Some(match split {
                    Flow::Source => source,
                    Flow::Sink => sink,
                }),
                HedgePair::Unpaired { hedge, .. } => Some(hedge),
            };
            if let Some(hedge) = hedge {
                ext.add(hedge);
            }
        }
        self.ext = ext;
        self.active_subgraph = Some(active_subgraph);
        self.rebuild_external_pull_topology();
        self
    }

    /// Load sharing and shared-X attachment demand, recomputed when the active
    /// graph changes and independent of route geometry or pull parameters.
    pub fn external_pull_topology(&self, edge: EdgeIndex) -> ExternalPullTopology {
        self.external_pull_topology[edge]
    }

    fn rebuild_external_pull_topology(&mut self)
    where
        E: HasPointConstraint,
        V: HasPointConstraint,
    {
        self.external_pull_topology =
            vec![ExternalPullTopology::default(); self.edge_points.len().0].into();
        if self.ext.included_iter().next().is_none() {
            return;
        }
        let n = self.vertex_points.len().0;
        let mut adjacency = vec![Vec::new(); n];
        for edge in self.active_edges() {
            let (_, pair) = &self.graph[&edge];
            if let HedgePair::Paired { source, sink } = *pair {
                if self.hedge_is_active(source) && self.hedge_is_active(sink) {
                    let a = self.graph.node_id(source).0;
                    let b = self.graph.node_id(sink).0;
                    if a != b {
                        adjacency[a].push(b);
                        adjacency[b].push(a);
                    }
                }
            }
        }
        let mut external = vec![Vec::new(); n];
        for hedge in self.ext.included_iter() {
            external[self.graph.node_id(hedge).0]
                .push((self.graph[&hedge], self.graph.flow(hedge)));
        }
        let mut visited = vec![false; n];
        let mut local_index = vec![0; n];
        let roots = self.active_nodes().collect::<Vec<_>>();
        for root in roots {
            if visited[root.0] {
                continue;
            }
            let mut nodes = vec![root.0];
            visited[root.0] = true;
            let mut cursor = 0;
            while cursor < nodes.len() {
                for &neighbor in &adjacency[nodes[cursor]] {
                    if !visited[neighbor] {
                        visited[neighbor] = true;
                        nodes.push(neighbor);
                    }
                }
                cursor += 1;
            }
            let mut incoming = 0;
            let mut outgoing = 0;
            for (i, &node) in nodes.iter().enumerate() {
                local_index[node] = i;
                for &(_, flow) in &external[node] {
                    match flow {
                        Flow::Source => outgoing += 1,
                        Flow::Sink => incoming += 1,
                    }
                }
            }
            if incoming == 0 || outgoing == 0 {
                continue;
            }

            // Unit total load on each side. Local incoming/outgoing attachments
            // cancel at their vertex before the internal spring network carries it.
            let demand = nodes
                .iter()
                .map(|&node| {
                    let local_outgoing = external[node]
                        .iter()
                        .filter(|(_, flow)| *flow == Flow::Source)
                        .count();
                    let local_incoming = external[node].len() - local_outgoing;
                    // Cancel equal rational loads exactly, independently of
                    // half-edge order, before converting to floating point.
                    let outward = local_outgoing as u128 * incoming as u128;
                    let inward = local_incoming as u128 * outgoing as u128;
                    let numerator = if outward >= inward {
                        (outward - inward) as f64
                    } else {
                        -((inward - outward) as f64)
                    };
                    numerator / (incoming as f64 * outgoing as f64)
                })
                .collect::<Vec<f64>>();
            let local_adjacency = nodes
                .iter()
                .map(|&node| {
                    adjacency[node]
                        .iter()
                        .map(|&neighbor| local_index[neighbor])
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>();
            let (width, potential) = Self::component_pull_distribution(
                &local_adjacency,
                &demand,
                (1.0 / incoming as f64).max(1.0 / outgoing as f64),
            );
            let mut groups: HashMap<_, Vec<(EdgeIndex, f64)>> = HashMap::new();
            for node in nodes {
                for &(edge, flow) in &external[node] {
                    let side_count = match flow {
                        Flow::Source => outgoing,
                        Flow::Sink => incoming,
                    };
                    self.external_pull_topology[edge].load_share =
                        (width / side_count as f64).min(1.0);
                    let Some(reference) = self.graph[edge]
                        .point_constraint()
                        .x
                        .force_target(LayoutPointIndex::Edge(edge))
                    else {
                        continue;
                    };
                    if !self.point_is_active(reference) {
                        continue;
                    }
                    if matches!(self.constraints(reference).x, Constraint::Fixed) {
                        continue;
                    }
                    let outgoing = flow == Flow::Source;
                    let sign = if outgoing { 1.0 } else { -1.0 };
                    groups
                        .entry((reference, outgoing))
                        .or_default()
                        .push((edge, sign * potential[local_index[node]]));
                }
            }
            // A common external X coordinate must clear the outermost owner.
            // In the unit-current model, each external stem has voltage drop
            // 1 / side_count. Moving the shared line from the mean owner to the
            // outermost owner adds this mean potential deficit to every stem's
            // load. Only actual shared coordinates receive the correction;
            // unshared endpoints form singleton groups and contribute exactly zero.
            for members in groups.values() {
                let outermost = members
                    .iter()
                    .map(|(_, value)| *value)
                    .fold(f64::NEG_INFINITY, f64::max);
                let deficit = members
                    .iter()
                    .map(|(_, value)| outermost - value)
                    .sum::<f64>()
                    / members.len() as f64;
                for &(edge, _) in members {
                    self.external_pull_topology[edge].attachment_deficit = width * deficit;
                }
            }
        }
    }

    /// Normalize by the greatest spring load, including the external stems.
    /// Unit-conductance edges share load in parallel; series edges retain the
    /// same load, so longer chains can still stretch naturally. Rest lengths
    /// do not enter this topology calculation, and self-loops conduct no load.
    /// Return the capacity and unit-current potentials used to compare owners.
    fn component_pull_distribution(
        adjacency: &[Vec<usize>],
        demand: &[f64],
        boundary_load: f64,
    ) -> (f64, Vec<f64>) {
        let n = adjacency.len();
        let ground = (0..n).max_by_key(|&i| adjacency[i].len()).unwrap_or(0);
        let apply_laplacian = |values: &[f64], result: &mut [f64]| {
            for (i, neighbors) in adjacency.iter().enumerate() {
                result[i] = if i == ground {
                    0.0
                } else {
                    neighbors.iter().map(|&j| values[i] - values[j]).sum()
                };
            }
        };
        let mut potential = vec![0.0; n];
        let mut residual = demand.to_vec();
        residual[ground] = 0.0;
        let rhs_norm_sq = residual.iter().map(|r| r * r).sum::<f64>();
        if rhs_norm_sq == 0.0 {
            return (1.0 / boundary_load, potential);
        }
        let tolerance_sq = rhs_norm_sq * 1e-22;
        let mut direction = residual
            .iter()
            .enumerate()
            .map(|(i, &r)| {
                if i == ground {
                    0.0
                } else {
                    r / adjacency[i].len() as f64
                }
            })
            .collect::<Vec<_>>();
        let mut residual_dot_preconditioned = residual
            .iter()
            .zip(&direction)
            .map(|(r, z)| r * z)
            .sum::<f64>();
        let mut product = vec![0.0; n];
        let mut converged = false;
        // Grounded Jacobi-preconditioned conjugate gradients use O(V + E)
        // storage. A bounded solve avoids introducing a dense matrix dependency.
        for _ in 0..n.saturating_mul(4).saturating_add(32) {
            apply_laplacian(&direction, &mut product);
            let curvature = direction
                .iter()
                .zip(&product)
                .map(|(p, q)| p * q)
                .sum::<f64>();
            if !curvature.is_finite() || curvature <= 0.0 {
                break;
            }
            let alpha = residual_dot_preconditioned / curvature;
            for i in 0..n {
                potential[i] += alpha * direction[i];
                residual[i] -= alpha * product[i];
            }
            let norm_sq = residual.iter().map(|r| r * r).sum::<f64>();
            if !norm_sq.is_finite() {
                break;
            }
            if norm_sq <= tolerance_sq {
                // Verify the actual residual, including the grounded equation;
                // recurrence roundoff must not silently overestimate capacity.
                let actual_sq = adjacency
                    .iter()
                    .enumerate()
                    .map(|(i, neighbors)| {
                        let current: f64 =
                            neighbors.iter().map(|&j| potential[i] - potential[j]).sum();
                        (demand[i] - current).powi(2)
                    })
                    .sum::<f64>();
                converged = actual_sq <= tolerance_sq * 4.0;
                break;
            }
            let next_dot = residual
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != ground)
                .map(|(i, &r)| r * r / adjacency[i].len() as f64)
                .sum::<f64>();
            let beta = next_dot / residual_dot_preconditioned;
            for i in 0..n {
                if i != ground {
                    direction[i] = residual[i] / adjacency[i].len() as f64 + beta * direction[i];
                }
            }
            residual_dot_preconditioned = next_dot;
        }
        if !converged {
            // One channel and no attachment correction are conservative if
            // conditioning prevents convergence. Discard unverified potentials.
            potential.fill(0.0);
            return (1.0, potential);
        }
        let max_load = adjacency
            .iter()
            .enumerate()
            .flat_map(|(i, neighbors)| neighbors.iter().map(move |&j| (i, j)))
            .map(|(i, j)| (potential[i] - potential[j]).abs())
            .fold(boundary_load, f64::max);
        (1.0 / max_load, potential)
    }

    pub(super) fn has_active_subgraph(&self) -> bool {
        self.active_subgraph.is_some()
    }

    pub(super) fn hedge_is_active(&self, hedge: Hedge) -> bool {
        self.active_subgraph
            .as_ref()
            .is_none_or(|active| active.includes(&hedge))
    }

    pub(super) fn node_is_active(&self, node: NodeIndex) -> bool {
        self.active_subgraph.as_ref().is_none_or(|active| {
            self.graph
                .iter_crown(node)
                .any(|hedge| active.includes(&hedge))
        })
    }

    pub(super) fn edge_is_active(&self, edge: EdgeIndex) -> bool {
        self.active_subgraph.as_ref().is_none_or(|active| {
            let (_, pair) = &self.graph[&edge];
            active.intersects(pair)
        })
    }

    /// Incoming and outgoing flows select horizontal pulling; a single flow is radial.
    pub(super) fn external_flows_are_mixed(&self) -> bool {
        let mut flows = self.ext.included_iter().map(|hedge| self.graph.flow(hedge));
        flows
            .next()
            .is_some_and(|first| flows.any(|flow| flow != first))
    }

    fn point_is_active(&self, point: LayoutPointIndex) -> bool {
        match point {
            LayoutPointIndex::Node(node) => self.node_is_active(node),
            LayoutPointIndex::Edge(edge) => self.edge_is_active(edge),
        }
    }

    fn constraints(&self, index: LayoutPointIndex) -> PointConstraint
    where
        E: HasPointConstraint,
        V: HasPointConstraint,
    {
        match index {
            LayoutPointIndex::Node(index) => *self.graph[index].point_constraint(),
            LayoutPointIndex::Edge(index) => *self.graph[index].point_constraint(),
        }
    }

    pub(super) fn active_nodes(&self) -> impl Iterator<Item = NodeIndex> + '_ {
        (0..self.vertex_points.len().0)
            .map(NodeIndex)
            .filter(|&node| self.node_is_active(node))
    }

    pub(super) fn active_edges(&self) -> impl Iterator<Item = EdgeIndex> + '_ {
        (0..self.edge_points.len().0)
            .map(EdgeIndex)
            .filter(|&edge| self.edge_is_active(edge))
    }

    /// Complete vertex-to-anchor path for one physical incidence.
    pub(super) fn half_route_points(&self, hedge: Hedge) -> Vec<Point2<f64>> {
        std::iter::once(self.vertex_points[self.graph.node_id(hedge)])
            .chain(self.route_points[hedge].iter().copied())
            .chain(std::iter::once(self.edge_points[self.graph[&hedge]]))
            .collect()
    }

    pub(super) fn half_route_length(&self, hedge: Hedge) -> f64 {
        self.half_route_points(hedge)
            .windows(2)
            .map(|points| points[0].distance(points[1]))
            .sum()
    }

    pub(super) fn edge_route_hedges(&self, edge: EdgeIndex) -> Vec<Hedge> {
        let (_, pair) = &self.graph[&edge];
        let hedges = match *pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                vec![source, sink]
            }
            HedgePair::Unpaired { hedge, .. } => vec![hedge],
        };
        hedges
            .into_iter()
            .filter(|&hedge| self.hedge_is_active(hedge))
            .collect()
    }

    /// One unit of edge charge shared by its anchor and active interior points.
    /// Position-independent weights keep force gradients consistent with energy.
    pub(super) fn edge_charge_samples(&self, edge: EdgeIndex) -> Vec<EdgeChargeSample> {
        if !self.edge_is_active(edge) {
            return Vec::new();
        }
        let mut samples = vec![(None, self.edge_points[edge], 1.0)];
        for hedge in self.edge_route_hedges(edge) {
            samples.extend(
                self.route_points[hedge]
                    .iter()
                    .enumerate()
                    .map(|(index, &point)| (Some((hedge, index)), point, 1.0)),
            );
        }
        let weight = 1.0 / samples.len() as f64;
        for sample in &mut samples {
            sample.2 = weight;
        }
        samples
    }

    fn active_route_points(&self) -> Vec<(Hedge, usize)> {
        self.route_points
            .iter()
            .filter(|(hedge, _)| self.hedge_is_active(*hedge))
            .flat_map(|(hedge, points)| (0..points.len()).map(move |index| (hedge, index)))
            .collect()
    }

    fn shift_route_point(&mut self, hedge: Hedge, index: usize, shift: Vector2<f64>) -> bool {
        if !self.hedge_is_active(hedge) || shift == Vector2::zero() {
            return false;
        }
        self.route_points[hedge][index] += shift;
        self.mark_edge_changed(self.graph[&hedge]);
        true
    }

    fn mark_node_changed(&mut self, index: NodeIndex) {
        self.changed_nodes.add(index);
    }

    fn mark_edge_changed(&mut self, index: EdgeIndex) {
        self.changed_edges.add(index);
    }

    pub fn clear_changes(&mut self) {
        self.changed_nodes.clear();
        self.changed_edges.clear();
    }
}

pub struct LayoutNeighbor;

impl<'a, E, V, H, N: NodeStorageOps<NodeData = V> + Clone> Neighbor<LayoutState<'a, E, V, H, N>>
    for LayoutNeighbor
{
    fn propose(
        &self,
        s: &LayoutState<'a, E, V, H, N>,
        rng: &mut impl Rng,
        step: f64,
        _temp: f64,
    ) -> LayoutState<'a, E, V, H, N> {
        let mut st = s.clone();
        let active_nodes = st.active_nodes().collect::<Vec<_>>();
        let active_edges = st.active_edges().collect::<Vec<_>>();
        let active_routes = st.active_route_points();
        let step_range: Uniform<f64> = Uniform::from(-step..step);
        if active_nodes.is_empty() && active_edges.is_empty() {
            return st;
        }

        let mut didnothing = true;
        let mut attempts = 0;
        while didnothing {
            attempts += 1;
            if st.has_active_subgraph() && attempts > 1024 {
                return st;
            }
            if !active_routes.is_empty() && rng.gen_bool(0.25) {
                let (hedge, index) = active_routes[rng.gen_range(0..active_routes.len())];
                let shift = LayoutNeighbor::axis_shift(&step_range, rng);
                didnothing = !st.shift_route_point(hedge, index, shift);
                continue;
            }
            match rng.gen_range(0..100) {
                0..=69 => {
                    // single-DOF
                    if rng.gen_bool(0.6) {
                        if active_nodes.is_empty() {
                            continue;
                        }
                        let v = active_nodes[rng.gen_range(0..active_nodes.len())];

                        let shift = LayoutNeighbor::axis_shift(&step_range, rng);
                        let changed = apply_vertex_shift(&mut st, v, shift);
                        didnothing = !changed;
                    } else {
                        if active_edges.is_empty() {
                            continue;
                        }
                        let e = active_edges[rng.gen_range(0..active_edges.len())];
                        let shift = LayoutNeighbor::axis_shift(&step_range, rng);
                        let changed = apply_edge_shift(&mut st, e, shift);
                        didnothing = !changed;
                    }
                }
                _ => {
                    // vertex block
                    if active_nodes.is_empty() {
                        continue;
                    }
                    let v = active_nodes[rng.gen_range(0..active_nodes.len())];

                    let shift = LayoutNeighbor::diagonal_shift(&step_range, rng, 0.6);

                    // Cache whether a vertex move succeeded so we can mark it.
                    let mut changed_any = apply_vertex_shift(&mut st, v, shift);
                    let incident_edges = st
                        .graph
                        .iter_crown(v)
                        .filter(|&hedge| st.hedge_is_active(hedge))
                        .map(|hedge| st.graph[&hedge])
                        .collect::<Vec<_>>();
                    for index in incident_edges {
                        // Propagate to incident edge control points; any change gets recorded.
                        changed_any |= apply_edge_shift(&mut st, index, shift);
                    }

                    didnothing = !changed_any;
                } // _ => {
                  //     // everything
                  //     let e = EdgeIndex(rng.gen_range(0..n_e.0));
                  //     st.edge_points[e].x += step_range.sample(rng) * 0.5;
                  //     st.edge_points[e].y += step_range.sample(rng) * 0.5;
                  // }
            }
        }
        st
    }
}

impl LayoutNeighbor {
    fn axis_shift(step_range: &Uniform<f64>, rng: &mut impl Rng) -> Vector2<f64> {
        if rng.gen_bool(0.5) {
            Vector2::from((step_range.sample(rng), 0.0))
        } else {
            Vector2::from((0.0, step_range.sample(rng)))
        }
    }

    fn diagonal_shift(step_range: &Uniform<f64>, rng: &mut impl Rng, scale: f64) -> Vector2<f64> {
        let mut sample = || step_range.sample(rng) * scale;
        Vector2::from((sample(), sample()))
    }
}

pub struct PinnedLayoutNeighbor;

impl<'a, E, V, H, N> Neighbor<LayoutState<'a, E, V, H, N>> for PinnedLayoutNeighbor
where
    E: HasPointConstraint,
    V: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    fn prepare(&self, state: &mut LayoutState<'a, E, V, H, N>) {
        state.synchronize_grouped_coordinates();
        state.clear_changes();
    }

    fn propose(
        &self,
        s: &LayoutState<'a, E, V, H, N>,
        rng: &mut impl Rng,
        step: f64,
        _temp: f64,
    ) -> LayoutState<'a, E, V, H, N> {
        let mut st = s.clone();
        st.synchronize_grouped_coordinates();
        let active_nodes = st.active_nodes().collect::<Vec<_>>();
        let active_edges = st.active_edges().collect::<Vec<_>>();
        let active_routes = st.active_route_points();
        let step_range: Uniform<f64> = Uniform::from(-step..step);
        if active_nodes.is_empty() && active_edges.is_empty() {
            return st;
        }

        let mut didnothing = true;
        let mut attempts = 0;
        while didnothing {
            attempts += 1;
            if st.has_active_subgraph() && attempts > 1024 {
                return st;
            }
            if !active_routes.is_empty() && rng.gen_bool(0.25) {
                let (hedge, index) = active_routes[rng.gen_range(0..active_routes.len())];
                let shift = LayoutNeighbor::axis_shift(&step_range, rng);
                didnothing = !apply_route_point_shift(&mut st, hedge, index, shift);
                continue;
            }
            match rng.gen_range(0..100) {
                0..=6 => {
                    // Reorder constrained lines without crossing their repulsive barrier.
                    didnothing = !st.swap_on_constrained_line(rng);
                }
                7..=69 => {
                    // single-DOF
                    if rng.gen_bool(0.6) {
                        if active_nodes.is_empty() {
                            continue;
                        }
                        let v = active_nodes[rng.gen_range(0..active_nodes.len())];

                        let mut shift = LayoutNeighbor::axis_shift(&step_range, rng);
                        let bias = directional_force_shift(
                            st.graph[v].point_constraint(),
                            LayoutPointIndex::Node(v),
                            st.vertex_points[v],
                            st.directional_force * step,
                        );
                        shift += bias;
                        let changed = apply_vertex_shift_with_groups(&mut st, v, shift);
                        didnothing = !changed;
                    } else {
                        if active_edges.is_empty() {
                            continue;
                        }
                        let e = active_edges[rng.gen_range(0..active_edges.len())];
                        let mut shift = LayoutNeighbor::axis_shift(&step_range, rng);
                        let bias = directional_force_shift(
                            st.graph[e].point_constraint(),
                            LayoutPointIndex::Edge(e),
                            st.edge_points[e],
                            st.directional_force * step,
                        );
                        shift += bias;
                        let changed = apply_edge_shift_with_groups(&mut st, e, shift);
                        didnothing = !changed;
                    }
                }
                _ => {
                    // vertex block
                    if active_nodes.is_empty() {
                        continue;
                    }
                    let v = active_nodes[rng.gen_range(0..active_nodes.len())];

                    let shift = LayoutNeighbor::diagonal_shift(&step_range, rng, 0.6);
                    let vertex_bias = directional_force_shift(
                        st.graph[v].point_constraint(),
                        LayoutPointIndex::Node(v),
                        st.vertex_points[v],
                        st.directional_force * step,
                    );

                    // Cache whether a vertex move succeeded so we can mark it.
                    let mut changed_any =
                        apply_vertex_shift_with_groups(&mut st, v, shift + vertex_bias);

                    let incident_edges = st
                        .graph
                        .iter_crown(v)
                        .filter(|&hedge| st.hedge_is_active(hedge))
                        .map(|hedge| st.graph[&hedge])
                        .collect::<Vec<_>>();
                    for index in incident_edges {
                        // Propagate to incident edge control points; any change gets recorded.
                        let edge_bias = directional_force_shift(
                            st.graph[index].point_constraint(),
                            LayoutPointIndex::Edge(index),
                            st.edge_points[index],
                            st.directional_force * step,
                        );
                        let edge_shift = shift + edge_bias;
                        changed_any |= apply_edge_shift_with_groups(&mut st, index, edge_shift);
                    }

                    didnothing = !changed_any;
                }
            }
        }
        st
    }
}

/// Move an interior point without changing its physical edge or anchor constraints.
/// Fixed owner axes apply to its bends; grouped coordinates refer only to the anchor.
pub(crate) fn apply_route_point_shift<'a, E, V, H, N>(
    state: &mut LayoutState<'a, E, V, H, N>,
    hedge: Hedge,
    index: usize,
    mut shift: Vector2<f64>,
) -> bool
where
    E: HasPointConstraint,
    N: NodeStorageOps<NodeData = V>,
{
    let constraints = state.graph[state.graph[&hedge]].point_constraint();
    if matches!(constraints.x, Constraint::Fixed) {
        shift.x = 0.0;
    }
    if matches!(constraints.y, Constraint::Fixed) {
        shift.y = 0.0;
    }
    state.shift_route_point(hedge, index, shift)
}

pub(crate) fn apply_vertex_shift<'a, E, V, H, N: NodeStorageOps<NodeData = V> + Clone>(
    state: &mut LayoutState<'a, E, V, H, N>,
    idx: NodeIndex,
    shift: Vector2<f64>,
) -> bool {
    if shift == Vector2::zero() {
        return false;
    }
    state.vertex_points[idx] += shift;
    state.mark_node_changed(idx);
    true
}

pub(crate) fn apply_edge_shift<'a, E, V, H, N: NodeStorageOps<NodeData = V> + Clone>(
    state: &mut LayoutState<'a, E, V, H, N>,
    idx: EdgeIndex,
    shift: Vector2<f64>,
) -> bool {
    if shift == Vector2::zero() {
        return false;
    }
    state.edge_points[idx] += shift;
    state.mark_edge_changed(idx);
    true
}

#[derive(Clone, Copy)]
enum LayoutAxis {
    X,
    Y,
}

impl<'a, E, V, H, N> LayoutState<'a, E, V, H, N>
where
    E: HasPointConstraint,
    V: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    fn group_members(
        &self,
        reference: LayoutPointIndex,
        axis: LayoutAxis,
    ) -> Vec<LayoutPointIndex> {
        let matches_reference = |constraints: PointConstraint| {
            matches!(
                match axis {
                    LayoutAxis::X => constraints.x,
                    LayoutAxis::Y => constraints.y,
                },
                Constraint::Grouped(group_reference, _) if group_reference == reference
            )
        };
        let mut members = Vec::new();
        for i in 0..self.vertex_points.len().0 {
            let index = NodeIndex(i);
            if self.node_is_active(index)
                && matches_reference(*self.graph[index].point_constraint())
            {
                members.push(LayoutPointIndex::Node(index));
            }
        }
        for i in 0..self.edge_points.len().0 {
            let index = EdgeIndex(i);
            if self.edge_is_active(index)
                && matches_reference(*self.graph[index].point_constraint())
            {
                members.push(LayoutPointIndex::Edge(index));
            }
        }
        members
    }

    fn shift_coordinate(&mut self, index: LayoutPointIndex, axis: LayoutAxis, shift: f64) {
        match index {
            LayoutPointIndex::Node(index) => {
                match axis {
                    LayoutAxis::X => self.vertex_points[index].x += shift,
                    LayoutAxis::Y => self.vertex_points[index].y += shift,
                }
                self.mark_node_changed(index);
            }
            LayoutPointIndex::Edge(index) => {
                match axis {
                    LayoutAxis::X => self.edge_points[index].x += shift,
                    LayoutAxis::Y => self.edge_points[index].y += shift,
                }
                self.mark_edge_changed(index);
            }
        }
    }

    fn coordinate(&self, index: LayoutPointIndex, axis: LayoutAxis) -> f64 {
        match (index, axis) {
            (LayoutPointIndex::Node(index), LayoutAxis::X) => self.vertex_points[index].x,
            (LayoutPointIndex::Node(index), LayoutAxis::Y) => self.vertex_points[index].y,
            (LayoutPointIndex::Edge(index), LayoutAxis::X) => self.edge_points[index].x,
            (LayoutPointIndex::Edge(index), LayoutAxis::Y) => self.edge_points[index].y,
        }
    }

    fn set_coordinate(&mut self, index: LayoutPointIndex, axis: LayoutAxis, value: f64) -> bool {
        if self.coordinate(index, axis) == value {
            return false;
        }
        let coordinate = match index {
            LayoutPointIndex::Node(index) => {
                self.mark_node_changed(index);
                &mut self.vertex_points[index]
            }
            LayoutPointIndex::Edge(index) => {
                self.mark_edge_changed(index);
                &mut self.edge_points[index]
            }
        };
        match axis {
            LayoutAxis::X => coordinate.x = value,
            LayoutAxis::Y => coordinate.y = value,
        }
        true
    }

    pub(crate) fn synchronize_grouped_coordinates(&mut self) {
        for i in 0..self.vertex_points.len().0 {
            let index = NodeIndex(i);
            let point = LayoutPointIndex::Node(index);
            if !self.point_is_active(point) {
                continue;
            }
            let constraints = self.constraints(point);
            if let Constraint::Grouped(reference, _) = constraints.x {
                if self.point_is_active(reference) {
                    let value = self.coordinate(reference, LayoutAxis::X);
                    self.set_coordinate(point, LayoutAxis::X, value);
                }
            }
            if let Constraint::Grouped(reference, _) = constraints.y {
                if self.point_is_active(reference) {
                    let value = self.coordinate(reference, LayoutAxis::Y);
                    self.set_coordinate(point, LayoutAxis::Y, value);
                }
            }
        }
        for i in 0..self.edge_points.len().0 {
            let index = EdgeIndex(i);
            let point = LayoutPointIndex::Edge(index);
            if !self.point_is_active(point) {
                continue;
            }
            let constraints = self.constraints(point);
            if let Constraint::Grouped(reference, _) = constraints.x {
                if self.point_is_active(reference) {
                    let value = self.coordinate(reference, LayoutAxis::X);
                    self.set_coordinate(point, LayoutAxis::X, value);
                }
            }
            if let Constraint::Grouped(reference, _) = constraints.y {
                if self.point_is_active(reference) {
                    let value = self.coordinate(reference, LayoutAxis::Y);
                    self.set_coordinate(point, LayoutAxis::Y, value);
                }
            }
        }
    }

    /// Movable coordinates on a common constrained line. Nodes and edge control
    /// points share the same coordinate space and can participate in the same move.
    fn constrained_lines(&self) -> Vec<(LayoutAxis, Vec<LayoutPointIndex>)> {
        let mut lines = Vec::new();
        for (axis, fixed_axis) in [
            (LayoutAxis::X, LayoutAxis::Y),
            (LayoutAxis::Y, LayoutAxis::X),
        ] {
            let mut points = self
                .active_nodes()
                .map(LayoutPointIndex::Node)
                .chain(self.active_edges().map(LayoutPointIndex::Edge))
                .filter(|&point| {
                    let constraint = self.constraints(point);
                    let (free, fixed) = match axis {
                        LayoutAxis::X => (constraint.x, constraint.y),
                        LayoutAxis::Y => (constraint.y, constraint.x),
                    };
                    matches!(free, Constraint::Free)
                        && matches!(fixed, Constraint::Fixed | Constraint::Grouped(_, _))
                        && self.group_members(point, axis).is_empty()
                        && self.coordinate(point, axis).is_finite()
                        && self.coordinate(point, fixed_axis).is_finite()
                })
                .collect::<Vec<_>>();
            points.sort_by(|&a, &b| {
                self.coordinate(a, fixed_axis)
                    .total_cmp(&self.coordinate(b, fixed_axis))
            });
            let mut start = 0;
            while start < points.len() {
                let fixed = self.coordinate(points[start], fixed_axis);
                let mut end = start + 1;
                while end < points.len()
                    && (self.coordinate(points[end], fixed_axis) - fixed).abs() <= 1e-9
                {
                    end += 1;
                }
                if end - start >= 2 {
                    let mut line = points[start..end].to_vec();
                    line.sort_by(|&a, &b| {
                        self.coordinate(a, axis)
                            .total_cmp(&self.coordinate(b, axis))
                    });
                    lines.push((axis, line));
                }
                start = end;
            }
        }
        lines
    }

    fn swap_line_coordinates(
        &mut self,
        a: LayoutPointIndex,
        b: LayoutPointIndex,
        axis: LayoutAxis,
    ) -> bool {
        let pa = self.coordinate(a, axis);
        let pb = self.coordinate(b, axis);
        if pa == pb {
            return false;
        }
        self.set_coordinate(a, axis, pb);
        self.set_coordinate(b, axis, pa);
        true
    }

    // Swaps along a pinned/grouped axis allow reordering on constrained lines
    // without forcing points through the short-range repulsion between them.
    fn swap_on_constrained_line(&mut self, rng: &mut impl Rng) -> bool {
        let lines = self.constrained_lines();
        if lines.is_empty() {
            return false;
        }
        let (axis, points) = &lines[rng.gen_range(0..lines.len())];
        let a = rng.gen_range(0..points.len());
        let mut b = rng.gen_range(0..points.len() - 1);
        if b >= a {
            b += 1;
        }
        self.swap_line_coordinates(points[a], points[b], *axis)
    }

    /// Try one deterministic forward/backward sweep of adjacent line exchanges.
    /// Only a strict decrease of the current planar energy is accepted. Bounding
    /// the sweep limits full energy evaluations; later epochs can improve further.
    pub(super) fn improve_constrained_lines(&mut self, energy: &SpringChargeEnergy) -> f64 {
        let mut lines = self.constrained_lines();
        if lines.is_empty() {
            return 0.0;
        }
        let mut current = energy.total_energy(self);
        if !current.is_finite() {
            return 0.0;
        }
        let original_nodes = self.vertex_points.clone();
        let original_edges = self.edge_points.clone();
        for (axis, points) in &mut lines {
            let count = points.len() - 1;
            for index in (0..count).chain((0..count).rev()) {
                let mut proposal = self.clone();
                if !proposal.swap_line_coordinates(points[index], points[index + 1], *axis) {
                    continue;
                }
                // Evaluate the full energy: node swaps also change incident
                // segment crossings, which a control-point-only delta can miss.
                let next = energy.total_energy(&proposal);
                let tolerance = 1e-12 * (1.0 + current.abs().max(next.abs()));
                if next.is_finite() && next < current - tolerance {
                    *self = proposal;
                    current = next;
                    points.swap(index, index + 1);
                }
            }
        }
        self.active_nodes()
            .map(|node| self.vertex_points[node].distance(original_nodes[node]))
            .chain(
                self.active_edges()
                    .map(|edge| self.edge_points[edge].distance(original_edges[edge])),
            )
            .fold(0.0, f64::max)
    }

    fn shift_axis(
        &mut self,
        index: LayoutPointIndex,
        axis: LayoutAxis,
        constraint: Constraint,
        shift: f64,
    ) -> bool {
        if shift == 0.0 {
            return false;
        }
        match constraint {
            Constraint::Fixed => false,
            Constraint::Free => {
                self.shift_coordinate(index, axis, shift);
                true
            }
            Constraint::Grouped(reference, _) if reference == index => {
                self.shift_coordinate(index, axis, shift);
                let value = self.coordinate(index, axis);
                for member in self.group_members(reference, axis) {
                    self.set_coordinate(member, axis, value);
                }
                true
            }
            Constraint::Grouped(_, _) => false,
        }
    }

    fn shift_constrained(&mut self, index: LayoutPointIndex, shift: Vector2<f64>) -> bool {
        if !self.point_is_active(index) {
            return false;
        }
        let constraints = self.constraints(index);
        let changed_x = self.shift_axis(index, LayoutAxis::X, constraints.x, shift.x);
        let changed_y = self.shift_axis(index, LayoutAxis::Y, constraints.y, shift.y);
        changed_x || changed_y
    }
}

pub(crate) fn apply_vertex_shift_with_groups<'a, E, V, H, N>(
    state: &mut LayoutState<'a, E, V, H, N>,
    idx: NodeIndex,
    shift: Vector2<f64>,
) -> bool
where
    E: HasPointConstraint,
    V: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    state.shift_constrained(LayoutPointIndex::Node(idx), shift)
}

pub(crate) fn apply_edge_shift_with_groups<'a, E, V, H, N>(
    state: &mut LayoutState<'a, E, V, H, N>,
    idx: EdgeIndex,
    shift: Vector2<f64>,
) -> bool
where
    E: HasPointConstraint,
    V: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    state.shift_constrained(LayoutPointIndex::Edge(idx), shift)
}

#[derive(Clone, Copy)]
pub struct SpringChargeEnergy {
    pub spring_length: f64,            // L
    pub k_spring: f64,                 // 1.0
    pub c_vv: f64,                     // vertex-vertex charge (≈ 0.14*L^3)
    pub dangling_charge: f64,          // dangling edge charge (≈ 0.14*L^3)
    pub dangling_centroid_charge: f64, // dangling edge vs node-centroid charge
    pub external_pull: f64,            // constant outward force on external endpoints
    /// Exponent of each topology weight: 0 gives uniform per-leg pull, 1 gives
    /// normalized pull, and larger values strengthen the correction.
    pub external_pull_balance: f64,
    /// Multiplier for distributed shared-X attachment demand; zero disables it.
    pub external_pull_attachment: f64,
    pub c_ev: f64,             // edge-vertex (≈ 0.028*L^3)
    pub c_ee_local: f64,       // edge-edge local (≈ 0.014*L^3)
    pub c_center: f64,         // central pull (dimensionless relative strength)
    pub crossing_penalty: f64, // crossing energy penalty (≈ penalty*L^2)
    pub eps: f64,              // softened distance (≈ eps*L)
}

impl<'a, E, V, H, N: NodeStorageOps<NodeData = V> + Clone> Energy<LayoutState<'a, E, V, H, N>>
    for SpringChargeEnergy
{
    fn energy(
        &self,
        prev: Option<(&LayoutState<'a, E, V, H, N>, f64)>,
        next: &LayoutState<'a, E, V, H, N>,
    ) -> f64 {
        if !next.incremental
            || next
                .route_points
                .iter()
                .any(|(_, points)| !points.is_empty())
            || prev.is_some_and(|(state, _)| {
                state
                    .route_points
                    .iter()
                    .any(|(_, points)| !points.is_empty())
            })
        {
            return self.total_energy(next);
        }

        if let Some((prev_state, prev_energy)) = prev {
            if next.changed_nodes.is_empty() && next.changed_edges.is_empty() {
                // No mutations since the last evaluation, so the cached value is exact.
                return prev_energy;
            }

            // Compute the delta in a single pass over affected terms.
            prev_energy
                + self.delta_energy(prev_state, next, &next.changed_nodes, &next.changed_edges)
        } else {
            // First evaluation falls back to the full O(n²) pass.
            self.total_energy(next)
        }
    }

    fn on_accept(&self, state: &mut LayoutState<'a, E, V, H, N>) {
        state.clear_changes();
    }
}

#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
pub struct ParamTuning {
    pub length_scale: f64,            // scales L: default 1.0
    pub k_spring: f64,                // spring stiffness: default 1.0
    pub beta: f64,                    // vertex–vertex strength
    pub gamma_dangling: f64,          // dangling edge vs vertex–vertex
    pub gamma_dangling_centroid: f64, // dangling edge vs node centroid
    pub external_pull: f64,           // constant external pull relative to k_spring * L
    /// Finite nonnegative topology-weight exponent, defaulting to 1.
    pub external_pull_balance: f64,
    /// Multiplier for distributed shared-X attachment demand; zero disables it.
    pub external_pull_attachment: f64,
    pub gamma_ev: f64,         // edge–vertex vs vertex–vertex
    pub gamma_ee: f64,         // local edge–edge vs vertex–vertex
    pub g_center: f64,         // central vs vertex–vertex
    pub crossing_penalty: f64, // fixed penalty per crossing
    pub eps: f64,              // softening epsilon
}

impl ParamTuning {
    pub fn add_to_global(&self, global_data: &mut GlobalData) {
        global_data
            .statements
            .insert("length_scale".to_string(), self.length_scale.to_string());
        global_data
            .statements
            .insert("k_spring".to_string(), self.k_spring.to_string());
        global_data
            .statements
            .insert("beta".to_string(), self.beta.to_string());
        global_data
            .statements
            .insert("gamma_ev".to_string(), self.gamma_ev.to_string());

        global_data.statements.insert(
            "gamma_dangling".to_string(),
            self.gamma_dangling.to_string(),
        );
        global_data.statements.insert(
            "gamma_dangling_centroid".to_string(),
            self.gamma_dangling_centroid.to_string(),
        );
        global_data
            .statements
            .insert("external_pull".to_string(), self.external_pull.to_string());
        global_data.statements.insert(
            "external_pull_balance".to_string(),
            self.external_pull_balance.to_string(),
        );
        global_data.statements.insert(
            "external_pull_attachment".to_string(),
            self.external_pull_attachment.to_string(),
        );
        global_data
            .statements
            .insert("gamma_ee".to_string(), self.gamma_ee.to_string());
        global_data
            .statements
            .insert("g_center".to_string(), self.g_center.to_string());
        global_data.statements.insert(
            "crossing_penalty".to_string(),
            self.crossing_penalty.to_string(),
        );
        global_data
            .statements
            .insert("eps".to_string(), self.eps.to_string());
    }
    pub fn parse(global_data: &GlobalData) -> Self {
        let mut tune = ParamTuning::default();

        for (key, value) in &global_data.statements {
            match key.as_str() {
                "length_scale" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.length_scale = v;
                    }
                }
                "gamma_dangling" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.gamma_dangling = v;
                    }
                }
                "gamma_dangling_centroid" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.gamma_dangling_centroid = v;
                    }
                }
                "external_pull" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.external_pull = v;
                    }
                }
                "external_pull_balance" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.external_pull_balance = v;
                    }
                }
                "external_pull_attachment" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.external_pull_attachment = v;
                    }
                }
                "k_spring" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.k_spring = v;
                    }
                }
                "beta" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.beta = v;
                    }
                }
                "gamma_ev" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.gamma_ev = v;
                    }
                }
                "gamma_ee" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.gamma_ee = v;
                    }
                }
                "g_center" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.g_center = v;
                    }
                }
                "crossing_penalty" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.crossing_penalty = v;
                    }
                }
                "eps" => {
                    if let Ok(v) = value.parse::<f64>() {
                        tune.eps = v;
                    }
                }
                _ => {}
            }
        }

        tune
    }
}
impl Default for ParamTuning {
    fn default() -> Self {
        Self {
            length_scale: 1.0,
            k_spring: 1.0,
            beta: 0.14,
            gamma_dangling: 0.14,
            gamma_dangling_centroid: 0.0,
            external_pull: 0.0,
            external_pull_balance: 1.0,
            external_pull_attachment: 1.0,
            gamma_ev: 0.20,
            gamma_ee: 0.10,
            g_center: 0.05,
            crossing_penalty: 0.0,
            eps: 1e-4,
        }
    }
}

#[cfg(feature = "energy_trace")]
mod energy_trace {
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::time::Duration;

    #[derive(Debug, Clone, Copy, Default)]
    pub struct EnergyTiming {
        pub total_energy_ns: u64,
        pub partial_energy_ns: u64,
        pub vv_ns: u64,
        pub ev_ns: u64,
        pub spring_ns: u64,
        pub ee_local_ns: u64,
        pub dangling_ns: u64,
        pub center_ns: u64,
        pub crossing_ns: u64,
    }

    static TOTAL_ENERGY_NS: AtomicU64 = AtomicU64::new(0);
    static PARTIAL_ENERGY_NS: AtomicU64 = AtomicU64::new(0);
    static VV_NS: AtomicU64 = AtomicU64::new(0);
    static EV_NS: AtomicU64 = AtomicU64::new(0);
    static SPRING_NS: AtomicU64 = AtomicU64::new(0);
    static EE_LOCAL_NS: AtomicU64 = AtomicU64::new(0);
    static DANGLING_NS: AtomicU64 = AtomicU64::new(0);
    static CENTER_NS: AtomicU64 = AtomicU64::new(0);
    static CROSSING_NS: AtomicU64 = AtomicU64::new(0);

    #[inline]
    pub fn record_total(d: Duration) {
        TOTAL_ENERGY_NS.fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
    }

    #[inline]
    pub fn record_partial(d: Duration) {
        PARTIAL_ENERGY_NS.fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
    }

    #[inline]
    pub fn record_vv(d: Duration) {
        VV_NS.fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
    }

    #[inline]
    pub fn record_ev(d: Duration) {
        EV_NS.fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
    }

    #[inline]
    pub fn record_spring(d: Duration) {
        SPRING_NS.fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
    }

    #[inline]
    pub fn record_ee_local(d: Duration) {
        EE_LOCAL_NS.fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
    }

    #[inline]
    pub fn record_dangling(d: Duration) {
        DANGLING_NS.fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
    }

    #[inline]
    pub fn record_center(d: Duration) {
        CENTER_NS.fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
    }

    #[inline]
    pub fn record_crossing(d: Duration) {
        CROSSING_NS.fetch_add(d.as_nanos() as u64, Ordering::Relaxed);
    }

    pub fn snapshot() -> EnergyTiming {
        EnergyTiming {
            total_energy_ns: TOTAL_ENERGY_NS.load(Ordering::Relaxed),
            partial_energy_ns: PARTIAL_ENERGY_NS.load(Ordering::Relaxed),
            vv_ns: VV_NS.load(Ordering::Relaxed),
            ev_ns: EV_NS.load(Ordering::Relaxed),
            spring_ns: SPRING_NS.load(Ordering::Relaxed),
            ee_local_ns: EE_LOCAL_NS.load(Ordering::Relaxed),
            dangling_ns: DANGLING_NS.load(Ordering::Relaxed),
            center_ns: CENTER_NS.load(Ordering::Relaxed),
            crossing_ns: CROSSING_NS.load(Ordering::Relaxed),
        }
    }

    pub fn reset() {
        TOTAL_ENERGY_NS.store(0, Ordering::Relaxed);
        PARTIAL_ENERGY_NS.store(0, Ordering::Relaxed);
        VV_NS.store(0, Ordering::Relaxed);
        EV_NS.store(0, Ordering::Relaxed);
        SPRING_NS.store(0, Ordering::Relaxed);
        EE_LOCAL_NS.store(0, Ordering::Relaxed);
        DANGLING_NS.store(0, Ordering::Relaxed);
        CENTER_NS.store(0, Ordering::Relaxed);
        CROSSING_NS.store(0, Ordering::Relaxed);
    }
}

#[cfg(feature = "energy_trace")]
pub use energy_trace::{
    reset as energy_timing_reset, snapshot as energy_timing_snapshot, EnergyTiming,
};

impl SpringChargeEnergy {
    /// Topology multiplier before the balancing exponent, shared with external
    /// layout solvers that apply their own force strength and exponent.
    pub fn external_pull_scale(&self, topology: ExternalPullTopology) -> f64 {
        let scale =
            topology.load_share + self.external_pull_attachment * topology.attachment_deficit;
        assert!(
            scale.is_finite(),
            "external pull attachment produces an unrepresentable topology scale"
        );
        scale
    }

    pub(super) fn external_pull_strength(&self, topology: ExternalPullTopology) -> f64 {
        if self.external_pull == 0.0 {
            return 0.0;
        }
        let scale = if self.external_pull_balance == 0.0 {
            1.0
        } else {
            let topology_scale = self.external_pull_scale(topology);
            match self.external_pull_balance {
                1.0 => topology_scale,
                balance => topology_scale.powf(balance),
            }
        };
        let strength = self.external_pull * scale;
        assert!(
            strength.is_finite(),
            "external pull and balancing produce an unrepresentable force"
        );
        strength
    }

    #[cfg_attr(feature = "energy_trace", inline(never))]
    fn vv_term(&self, dist: f64) -> f64 {
        0.5 * self.c_vv / (dist + self.eps)
    }

    #[cfg_attr(feature = "energy_trace", inline(never))]
    fn ev_term(&self, dist: f64) -> f64 {
        self.c_ev / (dist + self.eps)
    }

    #[cfg_attr(feature = "energy_trace", inline(never))]
    fn spring_term_with_length(&self, dist: f64, length: f64) -> f64 {
        let t = length - dist;
        0.5 * self.k_spring * (t * t)
    }

    #[cfg_attr(feature = "energy_trace", inline(never))]
    fn ee_local_term(&self, dist: f64) -> f64 {
        0.5 * self.c_ee_local / (dist + self.eps)
    }

    #[cfg_attr(feature = "energy_trace", inline(never))]
    fn dangling_term(&self, dist: f64) -> f64 {
        0.5 * self.dangling_charge / (dist + self.eps)
    }

    #[cfg_attr(feature = "energy_trace", inline(never))]
    fn dangling_centroid_term(&self, dist: f64) -> f64 {
        self.dangling_centroid_charge / (dist + self.eps)
    }

    fn dangling_centroid_energy<'a, E, V, H, N>(&self, state: &LayoutState<'a, E, V, H, N>) -> f64
    where
        N: NodeStorageOps<NodeData = V> + Clone,
    {
        let active_nodes = state.active_nodes().collect::<Vec<_>>();
        if (self.dangling_centroid_charge == 0.0 && self.external_pull == 0.0)
            || active_nodes.is_empty()
        {
            return 0.0;
        }

        let centroid = Point2::from_vec(
            active_nodes.iter().fold(Vector2::zero(), |sum, &node| {
                sum + state.vertex_points[node].to_vec()
            }) / active_nodes.len() as f64,
        );
        let horizontal = state.external_flows_are_mixed();
        state
            .ext
            .included_iter()
            .map(|hedge| {
                let edge = state.graph[&hedge];
                let point = state.edge_points[edge];
                let repulsion = if self.dangling_centroid_charge == 0.0 {
                    0.0
                } else {
                    self.dangling_centroid_term(point.distance(centroid))
                };
                let extension = if horizontal {
                    let side = match state.graph.flow(hedge) {
                        Flow::Source => 1.0,
                        Flow::Sink => -1.0,
                    };
                    side * (point.x - centroid.x)
                } else {
                    point.distance(centroid)
                };
                repulsion
                    - self.external_pull_strength(state.external_pull_topology[edge]) * extension
            })
            .sum()
    }

    #[cfg_attr(feature = "energy_trace", inline(never))]
    fn center_term(&self, r: f64) -> f64 {
        0.5 * self.c_center * r.powi(2)
    }

    #[cfg_attr(feature = "energy_trace", inline(never))]
    fn crossing_term(&self) -> f64 {
        self.crossing_penalty
    }

    fn edge_segments<'a, E, V, H, N>(
        s: &LayoutState<'a, E, V, H, N>,
        edge: EdgeIndex,
    ) -> Vec<(Point2<f64>, Point2<f64>, Option<NodeIndex>)>
    where
        N: NodeStorageOps<NodeData = V> + Clone,
    {
        let hedges = s.edge_route_hedges(edge);
        // An unrouted self-loop's out-and-back span does not describe its drawn loop.
        if hedges.len() == 2
            && s.graph.node_id(hedges[0]) == s.graph.node_id(hedges[1])
            && hedges.iter().all(|&hedge| s.route_points[hedge].is_empty())
        {
            return Vec::new();
        }
        hedges
            .into_iter()
            .flat_map(|hedge| {
                let node = s.graph.node_id(hedge);
                s.half_route_points(hedge)
                    .windows(2)
                    .enumerate()
                    .map(|(index, points)| (points[0], points[1], (index == 0).then_some(node)))
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    pub(super) fn edge_spring_length<'a, E, V, H, N>(
        s: &LayoutState<'a, E, V, H, N>,
        edge: EdgeIndex,
        base: f64,
    ) -> f64
    where
        N: NodeStorageOps<NodeData = V> + Clone,
    {
        let base = base * s.edge_spring_length_scales[edge];
        let (_, pair) = &s.graph[&edge];
        if let Some(active) = &s.active_subgraph {
            match pair.with_subgraph(active) {
                Some(HedgePair::Unpaired { .. } | HedgePair::Split { .. }) => base * 2.0,
                Some(HedgePair::Paired { .. }) | None => base,
            }
        } else if matches!(pair, HedgePair::Unpaired { .. }) {
            base * 2.0
        } else {
            base
        }
    }

    fn segments_cross(a: Point2<f64>, b: Point2<f64>, c: Point2<f64>, d: Point2<f64>) -> bool {
        fn orient(a: Point2<f64>, b: Point2<f64>, c: Point2<f64>) -> f64 {
            let ab = b - a;
            let ac = c - a;
            ab.x * ac.y - ab.y * ac.x
        }

        let o1 = orient(a, b, c);
        let o2 = orient(a, b, d);
        let o3 = orient(c, d, a);
        let o4 = orient(c, d, b);
        const EPS: f64 = 1e-12;

        if o1.abs() <= EPS || o2.abs() <= EPS || o3.abs() <= EPS || o4.abs() <= EPS {
            return false;
        }

        (o1 > 0.0) != (o2 > 0.0) && (o3 > 0.0) != (o4 > 0.0)
    }

    fn edge_crosses<'a, E, V, H, N>(
        s: &LayoutState<'a, E, V, H, N>,
        ei: EdgeIndex,
        ej: EdgeIndex,
    ) -> bool
    where
        N: NodeStorageOps<NodeData = V> + Clone,
    {
        let seg_i = Self::edge_segments(s, ei);
        if seg_i.is_empty() {
            return false;
        }
        let seg_j = Self::edge_segments(s, ej);
        if seg_j.is_empty() {
            return false;
        }
        for (a, b, ai) in &seg_i {
            for (c, d, aj) in &seg_j {
                if ai.is_some() && ai == aj {
                    continue;
                }
                if Self::segments_cross(*a, *b, *c, *d) {
                    return true;
                }
            }
        }
        false
    }

    fn total_energy<'a, E, V, H, N>(&self, s: &LayoutState<'a, E, V, H, N>) -> f64
    where
        N: NodeStorageOps<NodeData = V> + Clone,
    {
        #[cfg(feature = "energy_trace")]
        let total_start = std::time::Instant::now();

        let n = s.vertex_points.len().0;
        let m = s.edge_points.len().0;

        let mut energy = 0.0;

        #[cfg(feature = "energy_trace")]
        let vv_start = std::time::Instant::now();
        for i in 0..n {
            let ni = NodeIndex(i);
            if !s.node_is_active(ni) {
                continue;
            }
            let np = s.vertex_points[ni];
            for j in (i + 1)..n {
                let nj = NodeIndex(j);
                if !s.node_is_active(nj) {
                    continue;
                }
                let vj = s.vertex_points[nj];
                energy += self.vv_term(np.distance(vj));
            }
        }
        #[cfg(feature = "energy_trace")]
        energy_trace::record_vv(vv_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let ev_start = std::time::Instant::now();
        for i in 0..n {
            let ni = NodeIndex(i);
            if !s.node_is_active(ni) {
                continue;
            }
            let np = s.vertex_points[ni];
            if self.c_ev != 0.0 {
                for e in 0..m {
                    let ei = EdgeIndex(e);
                    if !s.edge_is_active(ei) {
                        continue;
                    }
                    for (_, point, weight) in s.edge_charge_samples(ei) {
                        energy += weight * self.ev_term(np.distance(point));
                    }
                }
            }
        }
        #[cfg(feature = "energy_trace")]
        energy_trace::record_ev(ev_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let spring_start = std::time::Instant::now();
        #[cfg(feature = "energy_trace")]
        let ee_local_start = std::time::Instant::now();
        for i in 0..n {
            let ni = NodeIndex(i);
            if !s.node_is_active(ni) {
                continue;
            }
            for e in s
                .graph
                .iter_crown(ni)
                .filter(|&hedge| s.hedge_is_active(hedge))
            {
                let ei = s.graph[&e];
                let length = Self::edge_spring_length(s, ei, self.spring_length);
                energy += self.spring_term_with_length(s.half_route_length(e), length);
                for e in s
                    .graph
                    .iter_crown(ni)
                    .filter(|&hedge| s.hedge_is_active(hedge))
                {
                    let ej = s.graph[&e];
                    if ei == ej {
                        continue;
                    }
                    for (_, first, first_weight) in s.edge_charge_samples(ei) {
                        for (_, second, second_weight) in s.edge_charge_samples(ej) {
                            energy += first_weight
                                * second_weight
                                * self.ee_local_term(first.distance(second));
                        }
                    }
                }
            }
        }
        #[cfg(feature = "energy_trace")]
        {
            energy_trace::record_spring(spring_start.elapsed());
            energy_trace::record_ee_local(ee_local_start.elapsed());
        }

        #[cfg(feature = "energy_trace")]
        let center_start = std::time::Instant::now();
        for i in 0..n {
            let ni = NodeIndex(i);
            if !s.node_is_active(ni) {
                continue;
            }
            let np = s.vertex_points[ni];
            if self.c_center != 0.0 {
                energy += self.center_term(np.distance(EuclideanSpace::origin()));
            }
        }
        #[cfg(feature = "energy_trace")]
        energy_trace::record_center(center_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let dangling_start = std::time::Instant::now();
        for hi in s.ext.included_iter() {
            for hj in s.ext.included_iter() {
                if hi >= hj {
                    continue;
                }
                let hi_idx = s.graph[&hi];
                let hj_idx = s.graph[&hj];
                let pi = s.edge_points[hi_idx];
                let pj = s.edge_points[hj_idx];
                energy += self.dangling_term(pi.distance(pj));
            }
        }
        energy += self.dangling_centroid_energy(s);
        #[cfg(feature = "energy_trace")]
        energy_trace::record_dangling(dangling_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let crossing_start = std::time::Instant::now();
        if self.crossing_penalty != 0.0 {
            let m = s.edge_points.len().0;
            for i in 0..m {
                let ei = EdgeIndex(i);
                if !s.edge_is_active(ei) {
                    continue;
                }
                for j in (i + 1)..m {
                    let ej = EdgeIndex(j);
                    if !s.edge_is_active(ej) {
                        continue;
                    }
                    if Self::edge_crosses(s, ei, ej) {
                        energy += self.crossing_term();
                    }
                }
            }
        }
        #[cfg(feature = "energy_trace")]
        energy_trace::record_crossing(crossing_start.elapsed());

        #[cfg(feature = "energy_trace")]
        energy_trace::record_total(total_start.elapsed());
        energy
    }

    /// Compute the energy delta for mutated nodes/edges in a single pass.
    fn delta_energy<'a, E, V, H, N>(
        &self,
        prev: &LayoutState<'a, E, V, H, N>,
        next: &LayoutState<'a, E, V, H, N>,
        node_changes: &SubSet<NodeIndex>,
        edge_changes: &SubSet<EdgeIndex>,
    ) -> f64
    where
        N: NodeStorageOps<NodeData = V> + Clone,
    {
        #[cfg(feature = "energy_trace")]
        let total_start = std::time::Instant::now();

        let n = next.vertex_points.len().0;
        let m = next.edge_points.len().0;
        let mut delta = 0.0;

        #[cfg(feature = "energy_trace")]
        let vv_start = std::time::Instant::now();
        for i in 0..n {
            let ni = NodeIndex(i);
            if !next.node_is_active(ni) {
                continue;
            }
            let node_changed = node_changes.includes(&ni);
            let prev_np = prev.vertex_points[ni];
            let next_np = next.vertex_points[ni];
            for j in (i + 1)..n {
                let nj = NodeIndex(j);
                if !next.node_is_active(nj) {
                    continue;
                }
                if !(node_changed || node_changes.includes(&nj)) {
                    continue;
                }
                let prev_vj = prev.vertex_points[nj];
                let next_vj = next.vertex_points[nj];
                delta += self.vv_term(next_np.distance(next_vj))
                    - self.vv_term(prev_np.distance(prev_vj));
            }
        }
        #[cfg(feature = "energy_trace")]
        energy_trace::record_vv(vv_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let ev_start = std::time::Instant::now();
        if self.c_ev != 0.0 {
            for i in 0..n {
                let ni = NodeIndex(i);
                if !next.node_is_active(ni) {
                    continue;
                }
                let node_changed = node_changes.includes(&ni);
                let prev_np = prev.vertex_points[ni];
                let next_np = next.vertex_points[ni];
                for e in 0..m {
                    let ei = EdgeIndex(e);
                    if !next.edge_is_active(ei) {
                        continue;
                    }
                    if !(node_changed || edge_changes.includes(&ei)) {
                        continue;
                    }
                    let prev_ep = prev.edge_points[ei];
                    let next_ep = next.edge_points[ei];
                    delta += self.ev_term(next_np.distance(next_ep))
                        - self.ev_term(prev_np.distance(prev_ep));
                }
            }
        }
        #[cfg(feature = "energy_trace")]
        energy_trace::record_ev(ev_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let spring_start = std::time::Instant::now();
        #[cfg(feature = "energy_trace")]
        let ee_local_start = std::time::Instant::now();
        for i in 0..n {
            let ni = NodeIndex(i);
            if !next.node_is_active(ni) {
                continue;
            }
            let node_changed = node_changes.includes(&ni);
            let prev_np = prev.vertex_points[ni];
            let next_np = next.vertex_points[ni];

            let include_node = node_changed;
            for hedge in next
                .graph
                .iter_crown(ni)
                .filter(|&hedge| next.hedge_is_active(hedge))
            {
                let ei = next.graph[&hedge];
                let edge_changed = edge_changes.includes(&ei);
                if !(include_node || edge_changed) {
                    continue;
                }
                let prev_ep = prev.edge_points[ei];
                let next_ep = next.edge_points[ei];
                let length = Self::edge_spring_length(next, ei, self.spring_length);
                delta += self.spring_term_with_length(next_np.distance(next_ep), length)
                    - self.spring_term_with_length(prev_np.distance(prev_ep), length);

                for other in next
                    .graph
                    .iter_crown(ni)
                    .filter(|&hedge| next.hedge_is_active(hedge))
                {
                    let ej = next.graph[&other];
                    if ei == ej {
                        continue;
                    }
                    let other_edge_changed = edge_changes.includes(&ej);
                    if !(include_node || edge_changed || other_edge_changed) {
                        continue;
                    }
                    let prev_ejp = prev.edge_points[ej];
                    let next_ejp = next.edge_points[ej];
                    delta += self.ee_local_term(next_ejp.distance(next_ep))
                        - self.ee_local_term(prev_ejp.distance(prev_ep));
                }
            }
        }
        #[cfg(feature = "energy_trace")]
        {
            energy_trace::record_spring(spring_start.elapsed());
            energy_trace::record_ee_local(ee_local_start.elapsed());
        }

        #[cfg(feature = "energy_trace")]
        let center_start = std::time::Instant::now();
        if self.c_center != 0.0 {
            for i in 0..n {
                let ni = NodeIndex(i);
                if !next.node_is_active(ni) || !node_changes.includes(&ni) {
                    continue;
                }
                let prev_np = prev.vertex_points[ni];
                let next_np = next.vertex_points[ni];
                let prev_r = prev_np.distance(EuclideanSpace::origin());
                let next_r = next_np.distance(EuclideanSpace::origin());
                delta -= self.center_term(prev_r);
                delta += self.center_term(next_r);
            }
        }
        #[cfg(feature = "energy_trace")]
        energy_trace::record_center(center_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let dangling_start = std::time::Instant::now();
        for hi in next.ext.included_iter() {
            let edge_i = next.graph[&hi];
            let edge_i_changed = edge_changes.includes(&edge_i);
            for hj in next.ext.included_iter() {
                if hi >= hj {
                    continue;
                }
                let edge_j = next.graph[&hj];
                let edge_j_changed = edge_changes.includes(&edge_j);
                if !(edge_i_changed || edge_j_changed) {
                    continue;
                }
                let prev_pi = prev.edge_points[edge_i];
                let prev_pj = prev.edge_points[edge_j];
                let next_pi = next.edge_points[edge_i];
                let next_pj = next.edge_points[edge_j];
                delta += self.dangling_term(next_pi.distance(next_pj))
                    - self.dangling_term(prev_pi.distance(prev_pj));
            }
        }
        if (self.dangling_centroid_charge != 0.0 || self.external_pull != 0.0)
            && (next.active_nodes().any(|node| node_changes.includes(&node))
                || next
                    .ext
                    .included_iter()
                    .any(|hedge| edge_changes.includes(&next.graph[&hedge])))
        {
            delta += self.dangling_centroid_energy(next) - self.dangling_centroid_energy(prev);
        }
        #[cfg(feature = "energy_trace")]
        energy_trace::record_dangling(dangling_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let crossing_start = std::time::Instant::now();
        if self.crossing_penalty != 0.0 {
            let m = next.edge_points.len().0;
            for i in 0..m {
                let ei = EdgeIndex(i);
                if !next.edge_is_active(ei) {
                    continue;
                }
                for j in (i + 1)..m {
                    let ej = EdgeIndex(j);
                    if !next.edge_is_active(ej)
                        || !(edge_changes.includes(&ei) || edge_changes.includes(&ej))
                    {
                        continue;
                    }
                    let prev_cross = Self::edge_crosses(prev, ei, ej);
                    let next_cross = Self::edge_crosses(next, ei, ej);
                    if prev_cross != next_cross {
                        if next_cross {
                            delta += self.crossing_term();
                        } else {
                            delta -= self.crossing_term();
                        }
                    }
                }
            }
        }
        #[cfg(feature = "energy_trace")]
        energy_trace::record_crossing(crossing_start.elapsed());

        #[cfg(feature = "energy_trace")]
        energy_trace::record_partial(total_start.elapsed());
        delta
    }

    pub fn from_graph(n_nodes: usize, viewport_w: f64, viewport_h: f64, tune: ParamTuning) -> Self {
        assert!(tune.external_pull_balance.is_finite() && tune.external_pull_balance >= 0.0);
        assert!(tune.external_pull_attachment.is_finite() && tune.external_pull_attachment >= 0.0);
        let area = (viewport_w * viewport_h).max(1e-9);
        let spring_length = tune.length_scale * (area / (n_nodes.max(1) as f64)).sqrt();
        let spring_length_sq = spring_length.powi(2);
        let repulsion_scale = spring_length.powi(3);

        SpringChargeEnergy {
            spring_length,
            k_spring: tune.k_spring,
            c_vv: tune.beta * repulsion_scale,
            c_ev: tune.beta * tune.gamma_ev * repulsion_scale,
            c_ee_local: tune.beta * tune.gamma_ee * repulsion_scale,
            c_center: tune.beta * tune.g_center,
            dangling_charge: tune.gamma_dangling * tune.beta * repulsion_scale,
            dangling_centroid_charge: tune.gamma_dangling_centroid * tune.beta * repulsion_scale,
            external_pull: tune.external_pull * tune.k_spring * spring_length,
            external_pull_balance: tune.external_pull_balance,
            external_pull_attachment: tune.external_pull_attachment,
            crossing_penalty: tune.crossing_penalty * spring_length_sq,
            eps: tune.eps * spring_length,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::half_edge::{
        builder::HedgeGraphBuilder, involution::Flow, nodestore::DefaultNodeStore, NoData,
    };
    use rand::{rngs::SmallRng, SeedableRng};

    fn normalized_pull_scales(
        n: usize,
        internal: &[(usize, usize)],
        external: &[(usize, Flow)],
    ) -> Vec<f64> {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let nodes = (0..n)
            .map(|_| builder.add_node(PointConstraint::default()))
            .collect::<Vec<_>>();
        for &(a, b) in internal {
            builder.add_edge(nodes[a], nodes[b], PointConstraint::default(), false);
        }
        for &(node, flow) in external {
            builder.add_external_edge(nodes[node], PointConstraint::default(), false, flow);
        }
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let state = graph.new_layout_state(
            vec![Point2::origin(); n].into(),
            vec![Point2::origin(); graph.n_edges()].into(),
            1.0,
            0.0,
            false,
        );
        external
            .iter()
            .enumerate()
            .map(|(i, _)| {
                test_energy(0.0).external_pull_scale(
                    state.external_pull_topology(EdgeIndex(internal.len() + i)),
                )
            })
            .collect()
    }

    #[test]
    fn external_pull_balance_defaults_roundtrips_and_preserves_endpoints() {
        assert_eq!(ParamTuning::default().external_pull_balance, 1.0);
        let mut global = GlobalData {
            name: String::new(),
            payload: None,
            statements: Default::default(),
            edge_statements: Default::default(),
            node_statements: Default::default(),
        };
        assert_eq!(ParamTuning::parse(&global).external_pull_balance, 1.0);
        for (balance, half_scale) in [
            (0.0, 1.0),
            (0.5, 0.5_f64.sqrt()),
            (1.0, 0.5),
            (2.0, 0.25),
            (2.5, 0.25 * 0.5_f64.sqrt()),
        ] {
            let tune = ParamTuning {
                external_pull: 2.0,
                external_pull_balance: balance,
                ..ParamTuning::default()
            };
            tune.add_to_global(&mut global);
            let parsed = ParamTuning::parse(&global);
            assert_eq!(parsed.external_pull_balance, balance);
            let energy = SpringChargeEnergy::from_graph(2, 3.0, 4.0, parsed);
            assert_eq!(energy.external_pull_balance, balance);
            assert!(
                (energy.external_pull_strength(ExternalPullTopology {
                    load_share: 0.5,
                    attachment_deficit: 0.0
                }) - energy.external_pull * half_scale)
                    .abs()
                    < 1e-12
            );
            if balance == 0.0 || balance == 1.0 {
                for scale in [0.0, 0.25, 1.0 / 3.0, 1.0] {
                    let expected = if balance == 0.0 {
                        energy.external_pull
                    } else {
                        energy.external_pull * scale
                    };
                    assert_eq!(
                        energy
                            .external_pull_strength(ExternalPullTopology {
                                load_share: scale,
                                attachment_deficit: 0.0
                            })
                            .to_bits(),
                        expected.to_bits()
                    );
                }
            }
        }
        for invalid in [-0.1, f64::NEG_INFINITY, f64::NAN, f64::INFINITY] {
            assert!(std::panic::catch_unwind(|| SpringChargeEnergy::from_graph(
                2,
                3.0,
                4.0,
                ParamTuning {
                    external_pull_balance: invalid,
                    ..ParamTuning::default()
                },
            ))
            .is_err());
        }
    }

    #[test]
    fn external_pull_attachment_roundtrips_and_short_circuits_unused_demand() {
        let mut global = GlobalData {
            name: String::new(),
            payload: None,
            statements: Default::default(),
            edge_statements: Default::default(),
            node_statements: Default::default(),
        };
        assert_eq!(ParamTuning::parse(&global).external_pull_attachment, 1.0);
        for attachment in [0.0, 1.0, 16.0] {
            let tune = ParamTuning {
                external_pull_attachment: attachment,
                ..ParamTuning::default()
            };
            tune.add_to_global(&mut global);
            let parsed = ParamTuning::parse(&global);
            assert_eq!(parsed.external_pull_attachment, attachment);
            assert_eq!(
                SpringChargeEnergy::from_graph(2, 3.0, 4.0, parsed).external_pull_attachment,
                attachment
            );
        }
        for invalid in [-0.1, f64::NEG_INFINITY, f64::NAN, f64::INFINITY] {
            assert!(std::panic::catch_unwind(|| SpringChargeEnergy::from_graph(
                2,
                3.0,
                4.0,
                ParamTuning {
                    external_pull_attachment: invalid,
                    ..ParamTuning::default()
                }
            ))
            .is_err());
        }
        let topology = ExternalPullTopology {
            load_share: 0.25,
            attachment_deficit: 2.0,
        };
        let mut energy = SpringChargeEnergy {
            external_pull: 2.0,
            external_pull_attachment: f64::MAX,
            ..test_energy(0.0)
        };
        assert!(std::panic::catch_unwind(|| energy.external_pull_strength(topology)).is_err());
        energy.external_pull_balance = 0.0;
        assert_eq!(energy.external_pull_strength(topology), 2.0);
        energy.external_pull_balance = 1.0;
        energy.external_pull = 0.0;
        assert_eq!(energy.external_pull_strength(topology), 0.0);
    }

    #[test]
    fn external_pull_balance_strengthens_correction_without_reversing_force() {
        let mut energy = SpringChargeEnergy {
            external_pull: 3.0,
            ..test_energy(0.0)
        };
        for topology_scale in [0.0, 1e-4, 0.25, 0.5, 1.0] {
            let mut previous = energy.external_pull;
            for balance in [0.0, 0.25, 0.5, 1.0, 2.0, 2.5, 100.0, f64::MAX] {
                energy.external_pull_balance = balance;
                let force = energy.external_pull_strength(ExternalPullTopology {
                    load_share: topology_scale,
                    attachment_deficit: 0.0,
                });
                assert!(force.is_finite() && force >= 0.0 && force <= previous);
                previous = force;
            }
        }
    }

    #[test]
    fn normalized_pull_shares_parallel_load_and_preserves_long_series_chains() {
        for length in [1, 4, 257] {
            for multiplicity in [1, 3] {
                let internal = (0..length)
                    .flat_map(|i| std::iter::repeat_n((i, i + 1), multiplicity))
                    .collect::<Vec<_>>();
                for legs in [1, 2, 3, 5] {
                    let external = (0..legs)
                        .map(|_| (0, Flow::Sink))
                        .chain((0..legs).map(|_| (length, Flow::Source)))
                        .collect::<Vec<_>>();
                    let scales = normalized_pull_scales(length + 1, &internal, &external);
                    let expected = (multiplicity as f64 / legs as f64).min(1.0);
                    for scale in scales {
                        assert!((scale - expected).abs() < 1e-8,
                            "length={length}, multiplicity={multiplicity}, legs={legs}: {scale} != {expected}");
                        assert!(scale <= 1.0);
                    }
                }
            }
        }
    }

    #[test]
    fn normalized_pull_uses_edge_currents_and_ignores_unloaded_cycles() {
        let external = [
            (0, Flow::Sink),
            (0, Flow::Sink),
            (0, Flow::Sink),
            (0, Flow::Sink),
            (1, Flow::Source),
            (1, Flow::Source),
            (1, Flow::Source),
            (1, Flow::Source),
        ];
        // One direct path and a path with three springs split load 3:1.
        // A mincut would count two channels equally; conductance alone would
        // incorrectly reduce force when both paths gain the same series length.
        let internal = [(0, 1), (0, 2), (2, 3), (3, 1)];
        let scales = normalized_pull_scales(4, &internal, &external);
        for scale in scales {
            assert!((scale - 1.0 / 3.0).abs() < 1e-12);
        }
        // Relabeling vertices, reversing edge orientation/order, and adding
        // self-loops or a pendant cycle cannot change load through the bridge.
        let bare = [(0, 1)];
        let with_cycles = [(0, 1), (0, 0), (1, 1), (0, 2), (2, 3), (3, 0)];
        for internal in [&bare[..], &with_cycles[..]] {
            for permutation in [[0, 1, 2, 3], [3, 1, 0, 2], [1, 2, 3, 0]] {
                let renamed = internal
                    .iter()
                    .rev()
                    .map(|&(a, b)| (permutation[b], permutation[a]))
                    .collect::<Vec<_>>();
                let external = external
                    .iter()
                    .map(|&(v, f)| (permutation[v], f))
                    .collect::<Vec<_>>();
                for scale in normalized_pull_scales(4, &renamed, &external) {
                    assert!((scale - 0.25).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn normalized_pull_handles_contacts_distributed_attachments_and_components() {
        for (incoming, outgoing) in [(1, 1), (2, 2), (2, 3), (51, 50)] {
            let external = (0..incoming)
                .map(|_| (0, Flow::Sink))
                .chain((0..outgoing).map(|_| (0, Flow::Source)))
                .collect::<Vec<_>>();
            let scales = normalized_pull_scales(1, &[], &external);
            for (i, scale) in scales.into_iter().enumerate() {
                let count = if i < incoming { incoming } else { outgoing };
                let expected = incoming.min(outgoing) as f64 / count as f64;
                assert!((scale - expected).abs() < 1e-12);
            }
        }
        // Mixed attachments at b cancel locally; the a-c bridge carries 2/3
        // of the unit side load. No singular source/sink contraction is needed.
        let external = [
            (0, Flow::Sink),
            (1, Flow::Sink),
            (1, Flow::Source),
            (2, Flow::Source),
            (2, Flow::Source),
        ];
        let expected = [0.75, 0.75, 0.5, 0.5, 0.5];
        for (scale, expected) in normalized_pull_scales(3, &[(0, 1), (0, 2)], &external)
            .into_iter()
            .zip(expected)
        {
            assert!((scale - expected).abs() < 1e-12);
        }
        // Rational cancellation is exact even when endpoint insertion order
        // would accumulate a nonzero rounding residual from repeated 1/6 terms.
        for legs in [3, 50] {
            let mut repeated_contacts = (0..2)
                .flat_map(|node| {
                    std::iter::repeat_n((node, Flow::Source), legs)
                        .chain(std::iter::repeat_n((node, Flow::Sink), legs))
                })
                .collect::<Vec<_>>();
            for _ in 0..2 {
                assert_eq!(
                    normalized_pull_scales(2, &[(0, 1)], &repeated_contacts),
                    vec![1.0; 4 * legs]
                );
                repeated_contacts.reverse();
            }
        }
        let contact_each_end = [
            (0, Flow::Sink),
            (0, Flow::Source),
            (1, Flow::Sink),
            (1, Flow::Source),
        ];
        assert_eq!(
            normalized_pull_scales(2, &[(0, 1)], &contact_each_end),
            vec![1.0; 4]
        );
        // A separate component, a homogeneous-flow component and an isolated
        // vertex must not change the load sharing of the first bridge.
        let internal = [(0, 1), (2, 3), (2, 3), (2, 3)];
        let external = [
            (0, Flow::Sink),
            (0, Flow::Sink),
            (1, Flow::Source),
            (1, Flow::Source),
            (2, Flow::Sink),
            (3, Flow::Source),
            (4, Flow::Source),
            (4, Flow::Source),
        ];
        assert_eq!(
            normalized_pull_scales(6, &internal, &external),
            vec![0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 1.0, 1.0]
        );
        assert!(normalized_pull_scales(3, &[], &[]).is_empty());
    }

    #[test]
    fn normalized_pull_accounts_for_actual_shared_x_attachments() {
        // Unit side current on 0--1--2 gives potentials (0, 1, 1.5).
        // Each outgoing stem carries 1/4, while a common X line must cover
        // an additional mean owner deficit of 1/4. The incoming owners coincide.
        for reverse in [false, true] {
            for mode in [
                "free",
                "fixed",
                "y",
                "shared",
                "separate",
                "free-reference",
                "fixed-reference",
                "node-reference",
                "fixed-node-reference",
            ] {
                let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
                let nodes = (0..3)
                    .map(|i| {
                        let mut constraint = PointConstraint::default();
                        if mode == "fixed-node-reference" && i == 0 {
                            constraint.x = Constraint::Fixed;
                        }
                        builder.add_node(constraint)
                    })
                    .collect::<Vec<_>>();
                for i in 0..2 {
                    builder.add_edge(nodes[i], nodes[i + 1], PointConstraint::default(), false);
                }
                let incoming = if reverse { Flow::Source } else { Flow::Sink };
                let outgoing = if reverse { Flow::Sink } else { Flow::Source };
                for _ in 0..2 {
                    builder.add_external_edge(
                        nodes[0],
                        PointConstraint::default(),
                        false,
                        incoming,
                    );
                }
                for (i, owner) in [1, 1, 2, 2].into_iter().enumerate() {
                    let reference = if matches!(mode, "node-reference" | "fixed-node-reference") {
                        LayoutPointIndex::Node(nodes[0])
                    } else if mode == "separate" {
                        LayoutPointIndex::Edge(EdgeIndex(4 + 2 * (owner - 1)))
                    } else {
                        LayoutPointIndex::Edge(EdgeIndex(4))
                    };
                    let group = Constraint::Grouped(reference, ShiftDirection::Any);
                    let mut constraint = PointConstraint::default();
                    match mode {
                        "fixed" => constraint.x = Constraint::Fixed,
                        "y" => constraint.y = group,
                        "shared" | "separate" | "node-reference" | "fixed-node-reference" => {
                            constraint.x = group;
                        }
                        "free-reference" if i != 0 => constraint.x = group,
                        "fixed-reference" => {
                            constraint.x = if i == 0 { Constraint::Fixed } else { group };
                        }
                        _ => {}
                    }
                    builder.add_external_edge(nodes[owner], constraint, false, outgoing);
                }
                let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
                let state = graph.new_layout_state(
                    vec![Point2::origin(); 3].into(),
                    vec![Point2::origin(); 8].into(),
                    1.0,
                    0.0,
                    false,
                );
                let corrected = matches!(mode, "shared" | "free-reference" | "node-reference");
                for attachment in [0.0, 1.0, 16.0] {
                    let energy = SpringChargeEnergy {
                        external_pull_attachment: attachment,
                        ..test_energy(0.0)
                    };
                    for edge in 2..8 {
                        let topology = state.external_pull_topology(EdgeIndex(edge));
                        let load_share = if edge < 4 { 0.5 } else { 0.25 };
                        let deficit = if corrected && edge >= 4 { 0.25 } else { 0.0 };
                        assert!((topology.load_share - load_share).abs() < 1e-12);
                        assert!((topology.attachment_deficit - deficit).abs() < 1e-12);
                        assert!(
                            (energy.external_pull_scale(topology)
                                - (load_share + attachment * deficit))
                                .abs()
                                < 1e-12,
                            "mode={mode}, reverse={reverse}, edge={edge}, attachment={attachment}"
                        );
                    }
                }
                assert_eq!(
                    state.external_pull_topology,
                    state.clone().external_pull_topology
                );
                if mode == "shared" {
                    // Removing one attachment at each owner changes the base
                    // load from 1/4 to 1/2, while retaining the 1/4 rail deficit.
                    let mut selected = graph.empty_subgraph::<SuBitGraph>();
                    for edge in [0, 1, 2, 3, 4, 6] {
                        selected.add(graph[&EdgeIndex(edge)].1);
                    }
                    let isolated = state.clone().with_active_subgraph(selected);
                    for edge in [4, 6] {
                        assert!(
                            (test_energy(0.0).external_pull_scale(
                                isolated.external_pull_topology(EdgeIndex(edge))
                            ) - 0.75)
                                .abs()
                                < 1e-12
                        );
                    }
                    for edge in [5, 7] {
                        assert_eq!(
                            test_energy(0.0).external_pull_scale(
                                isolated.external_pull_topology(EdgeIndex(edge))
                            ),
                            1.0
                        );
                    }
                    // An excluded group reference cannot move; its remaining
                    // members retain only the component's load-sharing baseline.
                    let mut selected = graph.empty_subgraph::<SuBitGraph>();
                    for edge in [0, 1, 2, 3, 5, 6, 7] {
                        selected.add(graph[&EdgeIndex(edge)].1);
                    }
                    let isolated = state.with_active_subgraph(selected);
                    for edge in [5, 6, 7] {
                        assert!(
                            (test_energy(0.0).external_pull_scale(
                                isolated.external_pull_topology(EdgeIndex(edge))
                            ) - 1.0 / 3.0)
                                .abs()
                                < 1e-12
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn normalized_pull_can_exceed_one_for_long_distributed_attachments() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let nodes = (0..=10)
            .map(|_| builder.add_node(PointConstraint::default()))
            .collect::<Vec<_>>();
        for i in 0..10 {
            builder.add_edge(nodes[i], nodes[i + 1], PointConstraint::default(), false);
        }
        builder.add_external_edge(nodes[0], PointConstraint::default(), false, Flow::Sink);
        let shared = PointConstraint {
            x: Constraint::Grouped(LayoutPointIndex::Edge(EdgeIndex(11)), ShiftDirection::Any),
            y: Constraint::Free,
        };
        for &node in &nodes[1..] {
            builder.add_external_edge(node, shared, false, Flow::Source);
        }
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let state = graph.new_layout_state(
            vec![Point2::origin(); 11].into(),
            vec![Point2::origin(); 21].into(),
            1.0,
            0.0,
            false,
        );
        for edge in 11..21 {
            assert!(
                (test_energy(0.0)
                    .external_pull_scale(state.external_pull_topology(EdgeIndex(edge)))
                    - 1.75)
                    .abs()
                    < 1e-10
            );
        }
        let mut energy = SpringChargeEnergy {
            external_pull: 2.0,
            ..test_energy(0.0)
        };
        for balance in [0.0, 0.25, 1.0, 2.0, 4.0] {
            energy.external_pull_balance = balance;
            assert_eq!(
                energy.external_pull_strength(ExternalPullTopology {
                    load_share: 1.75,
                    attachment_deficit: 0.0
                }),
                2.0 * 1.75_f64.powf(balance)
            );
        }
        energy.external_pull_balance = f64::MAX;
        assert!(
            std::panic::catch_unwind(|| energy.external_pull_strength(ExternalPullTopology {
                load_share: 1.75,
                attachment_deficit: 0.0
            }))
            .is_err()
        );
        energy.external_pull = 0.0;
        assert_eq!(
            energy.external_pull_strength(ExternalPullTopology {
                load_share: 1.75,
                attachment_deficit: 0.0
            }),
            0.0
        );
    }

    #[test]
    fn normalized_pull_rebuilds_for_active_edges_and_exposed_half_edges() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        builder.add_edge(a, b, PointConstraint::default(), false);
        for _ in 0..2 {
            builder.add_external_edge(a, PointConstraint::default(), false, Flow::Sink);
        }
        for _ in 0..2 {
            builder.add_external_edge(b, PointConstraint::default(), false, Flow::Source);
        }
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let state = graph.new_layout_state(
            vec![Point2::origin(); 2].into(),
            vec![Point2::origin(); 6].into(),
            1.0,
            0.0,
            false,
        );
        assert_eq!(
            state.external_pull_topology,
            vec![ExternalPullTopology::default(); 6].into()
        );
        let mut selected: SuBitGraph = graph.external_filter();
        let pair = graph[&EdgeIndex(0)].1;
        selected.add(pair);
        let isolated = state.clone().with_active_subgraph(selected);
        for i in 2..6 {
            assert_eq!(
                isolated.external_pull_topology[EdgeIndex(i)].load_share,
                0.5
            );
        }
        assert_eq!(
            state.external_pull_topology,
            vec![ExternalPullTopology::default(); 6].into()
        );
        let HedgePair::Paired { source, .. } = pair else {
            panic!("paired edge");
        };
        let mut half = graph.empty_subgraph::<SuBitGraph>();
        half.add(source);
        for hedge in graph.external_filter::<SuBitGraph>().included_iter() {
            if graph.node_id(hedge) == a {
                half.add(hedge);
            }
        }
        let half = state.with_active_subgraph(half);
        assert_eq!(half.external_pull_topology[EdgeIndex(0)].load_share, 1.0);
        for i in 2..4 {
            assert_eq!(half.external_pull_topology[EdgeIndex(i)].load_share, 0.5);
        }
        for i in 4..6 {
            assert_eq!(half.external_pull_topology[EdgeIndex(i)].load_share, 1.0);
        }
        assert_eq!(
            half.external_pull_topology,
            half.clone().external_pull_topology
        );
    }

    #[test]
    fn route_spring_length_is_invariant_under_collinear_subdivision() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let HedgePair::Paired { source, sink } = graph[&EdgeIndex(0)].1 else {
            panic!("expected paired edge");
        };
        let mut state = graph.new_layout_state(
            vec![Point2::new(-3.0, 0.0), Point2::new(3.0, 0.0)].into(),
            vec![Point2::origin()].into(),
            1.0,
            0.0,
            true,
        );
        let energy = test_energy(0.0);
        let original = energy.energy(None, &state);
        state.route_points[source] = vec![Point2::new(-2.0, 0.0), Point2::new(-0.3, 0.0)];
        state.route_points[sink] = vec![Point2::new(2.0, 0.0), Point2::new(1.0, 0.0)];
        assert_eq!(state.half_route_length(source), 3.0);
        assert_eq!(energy.energy(None, &state), original);
        assert_eq!(state.clone().route_points, state.route_points);
        state.route_points[source][0].y = 2.0;
        let length = 5.0_f64.sqrt() + 6.89_f64.sqrt() + 0.3;
        let expected = 0.5 * (length - 1.0).powi(2) + 2.0;
        assert!((energy.energy(None, &state) - expected).abs() < 1e-12);
    }

    #[test]
    fn route_charge_samples_normalize_active_physical_edge_charge() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let HedgePair::Paired { source, sink } = graph[&EdgeIndex(0)].1 else {
            panic!("expected paired edge");
        };
        let mut state = graph.new_layout_state(
            vec![Point2::new(-3.0, 0.0), Point2::new(3.0, 0.0)].into(),
            vec![Point2::origin()].into(),
            1.0,
            0.0,
            true,
        );
        assert_eq!(
            state.edge_charge_samples(EdgeIndex(0)),
            vec![(None, Point2::origin(), 1.0)]
        );
        state.route_points[source] = vec![Point2::new(-2.0, 0.0), Point2::new(-1.0, 0.0)];
        state.route_points[sink] = vec![Point2::new(1.0, 0.0)];
        let samples = state.edge_charge_samples(EdgeIndex(0));
        assert_eq!(samples.len(), 4);
        assert_eq!(samples.iter().map(|sample| sample.2).sum::<f64>(), 1.0);
        let energy = SpringChargeEnergy {
            k_spring: 0.0,
            c_ev: 3.0,
            ..test_energy(0.0)
        };
        let expected = state
            .vertex_points
            .iter()
            .flat_map(|(_, vertex)| {
                samples
                    .iter()
                    .map(|(_, point, weight)| weight * energy.ev_term(vertex.distance(*point)))
            })
            .sum::<f64>();
        assert_eq!(energy.energy(None, &state), expected);
        let mut active = graph.empty_subgraph::<SuBitGraph>();
        active.add(source);
        state = state.with_active_subgraph(active);
        let samples = state.edge_charge_samples(EdgeIndex(0));
        assert_eq!(samples.len(), 3);
        assert!(samples
            .iter()
            .all(|(address, _, _)| address.is_none_or(|(hedge, _)| hedge == source)));
        assert!((samples.iter().map(|sample| sample.2).sum::<f64>() - 1.0).abs() < 1e-12);
        let before = energy.energy(None, &state);
        state.route_points[sink][0] = Point2::new(1000.0, 1000.0);
        assert_eq!(energy.energy(None, &state), before);
    }

    #[test]
    fn route_point_shifts_preserve_fixed_axes_and_inactive_half_edges() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(
            a,
            b,
            PointConstraint {
                x: Constraint::Fixed,
                y: Constraint::Free,
            },
            false,
        );
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let HedgePair::Paired { source, sink } = graph[&EdgeIndex(0)].1 else {
            panic!("expected paired edge");
        };
        let mut state = graph.new_layout_state(
            vec![Point2::new(-3.0, 0.0), Point2::new(3.0, 0.0)].into(),
            vec![Point2::origin()].into(),
            1.0,
            0.0,
            true,
        );
        state.route_points[source] = vec![Point2::new(-1.0, 2.0)];
        assert!(apply_route_point_shift(
            &mut state,
            source,
            0,
            Vector2::new(4.0, 5.0)
        ));
        assert_eq!(state.route_points[source][0], Point2::new(-1.0, 7.0));
        assert!(state.changed_edges.includes(&EdgeIndex(0)));
        state.clear_changes();
        let mut active = graph.empty_subgraph::<SuBitGraph>();
        active.add(sink);
        state = state.with_active_subgraph(active);
        assert!(!apply_route_point_shift(
            &mut state,
            source,
            0,
            Vector2::new(4.0, 5.0)
        ));
        assert_eq!(state.route_points[source][0], Point2::new(-1.0, 7.0));
        assert!(state.changed_edges.is_empty());
    }

    #[test]
    fn route_annealing_moves_bends_and_evaluates_their_complete_energy() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let HedgePair::Paired { source, sink } = graph[&EdgeIndex(0)].1 else {
            panic!("expected paired edge");
        };
        let mut state = graph.new_layout_state(
            vec![Point2::new(-3.0, 0.0), Point2::new(3.0, 0.0)].into(),
            vec![Point2::origin()].into(),
            1.0,
            0.0,
            true,
        );
        state.route_points[source] = vec![Point2::new(-1.0, 1.0)];
        state.route_points[sink] = vec![Point2::new(1.0, -1.0)];
        let energy = SpringChargeEnergy {
            c_ev: 0.4,
            c_ee_local: 0.2,
            ..test_energy(0.1)
        };
        let mut rng = SmallRng::seed_from_u64(42);
        let mut moves = 0;
        for _ in 0..100 {
            let before = energy.energy(None, &state);
            let mut next = PinnedLayoutNeighbor.propose(&state, &mut rng, 0.1, 1.0);
            moves += usize::from(next.route_points != state.route_points);
            assert_eq!(
                energy.energy(Some((&state, before)), &next),
                energy.total_energy(&next)
            );
            energy.on_accept(&mut next);
            state = next;
        }
        assert!(moves > 0);
    }

    #[test]
    fn route_crossing_penalty_includes_routed_self_loops_only() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        let c = builder.add_node(PointConstraint::default());
        builder.add_edge(a, a, PointConstraint::default(), false);
        builder.add_edge(b, c, PointConstraint::default(), false);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut state = graph.new_layout_state(
            vec![
                Point2::new(-3.0, 0.0),
                Point2::new(-1.0, 2.0),
                Point2::new(1.0, 2.0),
            ]
            .into(),
            vec![Point2::origin(), Point2::new(0.5, 2.0)].into(),
            1.0,
            0.0,
            true,
        );
        assert!(SpringChargeEnergy::edge_segments(&state, EdgeIndex(0)).is_empty());
        let HedgePair::Paired { source, sink } = graph[&EdgeIndex(0)].1 else {
            panic!("expected paired loop");
        };
        state.route_points[source] = vec![Point2::new(-3.0, 3.0), Point2::new(0.0, 3.0)];
        state.route_points[sink] = vec![Point2::new(-3.0, -3.0), Point2::new(0.0, -3.0)];
        assert!(SpringChargeEnergy::edge_crosses(
            &state,
            EdgeIndex(0),
            EdgeIndex(1)
        ));
    }

    #[test]
    fn route_crossing_penalty_checks_bends_beyond_a_shared_vertex() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        let c = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        builder.add_edge(a, c, PointConstraint::default(), false);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut state = graph.new_layout_state(
            vec![
                Point2::origin(),
                Point2::new(4.0, 0.0),
                Point2::new(4.0, 4.0),
            ]
            .into(),
            vec![Point2::new(2.0, 0.0), Point2::new(2.0, 4.0)].into(),
            1.0,
            0.0,
            true,
        );
        assert!(!SpringChargeEnergy::edge_crosses(
            &state,
            EdgeIndex(0),
            EdgeIndex(1)
        ));
        let first = graph
            .iter_crown(a)
            .find(|&hedge| graph[&hedge] == EdgeIndex(0))
            .unwrap();
        let second = graph
            .iter_crown(a)
            .find(|&hedge| graph[&hedge] == EdgeIndex(1))
            .unwrap();
        state.route_points[first] = vec![Point2::new(0.0, 3.0), Point2::new(3.0, 3.0)];
        state.route_points[second] = vec![Point2::new(0.0, 1.0), Point2::new(3.0, 1.0)];
        assert!(SpringChargeEnergy::edge_crosses(
            &state,
            EdgeIndex(0),
            EdgeIndex(1)
        ));
        let energy = SpringChargeEnergy {
            k_spring: 0.0,
            crossing_penalty: 7.0,
            ..test_energy(0.0)
        };
        assert_eq!(energy.energy(None, &state), 7.0);
    }

    #[test]
    fn directional_force_only_restores_wrong_side() {
        let reference = LayoutPointIndex::Node(NodeIndex(0));
        let constraints = PointConstraint {
            x: Constraint::Grouped(reference, ShiftDirection::PositiveOnly),
            y: Constraint::Grouped(reference, ShiftDirection::NegativeOnly),
        };

        assert_eq!(
            directional_force_shift(&constraints, reference, Point2::new(-1.0, 1.0), 2.0),
            Vector2::new(2.0, -2.0)
        );
        assert_eq!(
            directional_force_shift(&constraints, reference, Point2::new(0.0, 0.0), 2.0),
            Vector2::new(2.0, -2.0)
        );
        assert_eq!(
            directional_force_shift(&constraints, reference, Point2::new(1.0, -1.0), 2.0),
            Vector2::zero()
        );
        assert_eq!(
            directional_force_shift(
                &constraints,
                LayoutPointIndex::Node(NodeIndex(1)),
                Point2::new(-1.0, 1.0),
                2.0,
            ),
            Vector2::zero()
        );
    }

    fn test_energy(c_center: f64) -> SpringChargeEnergy {
        SpringChargeEnergy {
            spring_length: 1.0,
            k_spring: 1.0,
            c_vv: 0.0,
            dangling_charge: 0.0,
            dangling_centroid_charge: 0.0,
            external_pull: 0.0,
            external_pull_balance: 1.0,
            external_pull_attachment: 1.0,
            c_ev: 0.0,
            c_ee_local: 0.0,
            c_center,
            crossing_penalty: 0.0,
            eps: 1e-4,
        }
    }

    #[test]
    fn constrained_line_swaps_lower_energy_on_fixed_and_grouped_rotated_lines() {
        let fixed = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Fixed,
        };
        for (rotated, reflection) in [(false, 1.0), (true, 1.0), (false, -1.0), (true, -1.0)] {
            for grouped in [false, true] {
                let rail = if grouped {
                    Constraint::Grouped(LayoutPointIndex::Edge(EdgeIndex(0)), ShiftDirection::Any)
                } else {
                    Constraint::Fixed
                };
                let constraint = if rotated {
                    PointConstraint {
                        x: Constraint::Free,
                        y: rail,
                    }
                } else {
                    PointConstraint {
                        x: rail,
                        y: Constraint::Free,
                    }
                };
                let point = |x, y| {
                    if rotated {
                        Point2::new(reflection * y, x)
                    } else {
                        Point2::new(x, reflection * y)
                    }
                };
                let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
                let a = builder.add_node(fixed);
                let b = builder.add_node(fixed);
                builder.add_external_edge(a, constraint, false, Flow::Source);
                builder.add_external_edge(b, constraint, false, Flow::Source);
                let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
                let mut state = graph.new_layout_state(
                    vec![point(-1.0, -2.0), point(-1.0, 2.0)].into(),
                    vec![point(1.0, 2.0), point(1.0, -2.0)].into(),
                    1.0,
                    0.0,
                    true,
                );
                state.vertex_depths = vec![Some(3.0), Some(-4.0)].into();
                state.edge_depths = vec![Some(5.0), Some(-6.0)].into();
                state.edge_depth_pins = vec![true; 2].into();
                state.edge_spring_length_scales = vec![1.5; 2].into();
                let before = state.clone();
                let energy = SpringChargeEnergy {
                    crossing_penalty: 10.0,
                    ..test_energy(0.0)
                };
                let initial_energy = energy.total_energy(&state);
                assert_eq!(state.improve_constrained_lines(&energy), 4.0);
                assert!(energy.total_energy(&state) < initial_energy);
                assert_eq!(
                    state.edge_points,
                    vec![point(1.0, -2.0), point(1.0, 2.0)].into()
                );
                assert_eq!(state.vertex_points, before.vertex_points);
                assert_eq!(state.vertex_depths, before.vertex_depths);
                assert_eq!(state.edge_depths, before.edge_depths);
                assert_eq!(state.edge_depth_pins, before.edge_depth_pins);
                assert_eq!(
                    state.edge_spring_length_scales,
                    before.edge_spring_length_scales
                );
                assert_eq!(state.improve_constrained_lines(&energy), 0.0);
                // Annealing uses the same coordinate-space move and preserves
                // the same constraints, with its own acceptance policy.
                let mut proposal = before;
                assert!(proposal.swap_on_constrained_line(&mut SmallRng::seed_from_u64(4)));
                assert_eq!(proposal.vertex_points, state.vertex_points);
                assert_eq!(proposal.edge_points, state.edge_points);
            }
        }
    }

    #[test]
    fn constrained_line_swaps_can_exchange_node_and_edge_but_protect_groups_and_complement() {
        let fixed = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Fixed,
        };
        let rail = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Free,
        };
        let y_group = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Grouped(LayoutPointIndex::Edge(EdgeIndex(2)), ShiftDirection::Any),
        };
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(rail);
        let b = builder.add_node(fixed);
        let c = builder.add_node(fixed);
        let inactive = builder.add_node(rail);
        builder.add_edge(a, b, fixed, false);
        builder.add_external_edge(c, rail, false, Flow::Source);
        builder.add_external_edge(c, y_group, false, Flow::Source);
        builder.add_external_edge(c, y_group, false, Flow::Source);
        builder.add_external_edge(inactive, rail, false, Flow::Source);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut selected = graph.empty_subgraph::<SuBitGraph>();
        for (pair, edge, _) in graph.iter_edges() {
            if edge != EdgeIndex(4) {
                selected.add(pair);
            }
        }
        let mut state = graph
            .new_layout_state(
                vec![
                    Point2::new(1.0, 2.0),
                    Point2::new(-1.0, -2.0),
                    Point2::new(-1.0, 2.0),
                    Point2::new(1.0, -20.0),
                ]
                .into(),
                vec![
                    Point2::new(0.0, -2.0),
                    Point2::new(1.0, -2.0),
                    Point2::new(1.0, 5.0),
                    Point2::new(1.0, 5.0),
                    Point2::new(1.0, 20.0),
                ]
                .into(),
                1.0,
                0.0,
                true,
            )
            .with_active_subgraph(selected);
        state.vertex_depths = vec![Some(7.0); 4].into();
        state.edge_depths = vec![Some(-8.0); 5].into();
        let before = state.clone();
        let energy = test_energy(0.0);
        let initial_energy = energy.total_energy(&state);
        assert_eq!(state.improve_constrained_lines(&energy), 4.0);
        assert!(energy.total_energy(&state) < initial_energy);
        assert_eq!(state.vertex_points[a], Point2::new(1.0, -2.0));
        assert_eq!(state.edge_points[EdgeIndex(1)], Point2::new(1.0, 2.0));
        for node in [b, c, inactive] {
            assert_eq!(state.vertex_points[node], before.vertex_points[node]);
        }
        for edge in [EdgeIndex(0), EdgeIndex(2), EdgeIndex(3), EdgeIndex(4)] {
            assert_eq!(state.edge_points[edge], before.edge_points[edge]);
        }
        assert_eq!(state.vertex_depths, before.vertex_depths);
        assert_eq!(state.edge_depths, before.edge_depths);
        assert_eq!(
            state.changed_nodes.included_iter().collect::<Vec<_>>(),
            vec![a]
        );
        assert_eq!(
            state.changed_edges.included_iter().collect::<Vec<_>>(),
            vec![EdgeIndex(1)]
        );
        assert_eq!(state.improve_constrained_lines(&energy), 0.0);
    }

    #[test]
    fn constrained_line_swaps_leave_free_layouts_and_equal_energy_orders_unchanged() {
        for constraint in [
            PointConstraint::default(),
            PointConstraint {
                x: Constraint::Fixed,
                y: Constraint::Free,
            },
        ] {
            let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
            let a = builder.add_node(PointConstraint {
                x: Constraint::Fixed,
                y: Constraint::Fixed,
            });
            builder.add_external_edge(a, constraint, false, Flow::Source);
            builder.add_external_edge(a, constraint, false, Flow::Source);
            let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
            let mut state = graph.new_layout_state(
                vec![Point2::origin()].into(),
                vec![Point2::new(1.0, -2.0), Point2::new(1.0, 2.0)].into(),
                1.0,
                0.0,
                true,
            );
            let before = state.clone();
            let energy = test_energy(0.0);
            assert_eq!(state.improve_constrained_lines(&energy), 0.0);
            assert_eq!(state.vertex_points, before.vertex_points);
            assert_eq!(state.edge_points, before.edge_points);
            assert!(state.changed_nodes.is_empty());
            assert!(state.changed_edges.is_empty());
        }
    }

    #[test]
    fn edge_spring_length_scales_are_local_and_preserve_dangling_factor() {
        for split in [None, Some(Flow::Source), Some(Flow::Sink)] {
            let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
            let a = builder.add_node(PointConstraint::default());
            let b = builder.add_node(PointConstraint::default());
            builder.add_edge(a, b, PointConstraint::default(), false);
            builder.add_external_edge(a, PointConstraint::default(), false, Flow::Source);
            builder.add_external_edge(b, PointConstraint::default(), false, Flow::Sink);
            let mut graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
            if let Some(split) = split {
                let HedgePair::Paired { source, sink } = graph[&EdgeIndex(0)].1 else {
                    panic!("expected paired edge");
                };
                graph[&EdgeIndex(0)].1 = HedgePair::Split {
                    source,
                    sink,
                    split,
                };
            }
            let mut state = graph.new_layout_state(
                vec![Point2::origin(), Point2::new(6.0, 0.0)].into(),
                vec![
                    Point2::new(2.0, 0.0),
                    Point2::new(0.0, 4.0),
                    Point2::new(6.0, 4.0),
                ]
                .into(),
                1.0,
                0.0,
                false,
            );
            let energy = SpringChargeEnergy {
                spring_length: 3.0,
                ..test_energy(0.0)
            };
            assert_eq!(state.edge_spring_length_scales, vec![1.0; 3].into());
            for (scales, lengths, expected_energy) in [
                ([1.0, 1.0, 1.0], [3.0, 6.0, 6.0], 5.0),
                ([0.5, 1.0, 1.0], [1.5, 6.0, 6.0], 7.25),
                ([0.5, 2.0, 1.0], [1.5, 12.0, 6.0], 37.25),
            ] {
                state.edge_spring_length_scales = scales.to_vec().into();
                for (edge, length) in lengths.into_iter().enumerate() {
                    assert_eq!(
                        SpringChargeEnergy::edge_spring_length(
                            &state,
                            EdgeIndex(edge),
                            energy.spring_length,
                        ),
                        length,
                    );
                }
                assert_eq!(energy.energy(None, &state), expected_energy);
                let mut cloned = state.clone();
                assert_eq!(
                    cloned.edge_spring_length_scales,
                    state.edge_spring_length_scales
                );
                assert_eq!(energy.energy(None, &cloned), expected_energy);
                cloned.edge_spring_length_scales[EdgeIndex(0)] = 4.0;
                assert_eq!(state.edge_spring_length_scales[EdgeIndex(0)], scales[0]);
            }
        }
    }

    #[test]
    fn heterogeneous_edge_spring_length_scales_incremental_energy_matches_full() {
        for second_flow in [Flow::Sink, Flow::Source] {
            let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
            let a = builder.add_node(PointConstraint::default());
            let b = builder.add_node(PointConstraint::default());
            let c = builder.add_node(PointConstraint::default());
            builder.add_edge(a, b, PointConstraint::default(), false);
            builder.add_edge(b, c, PointConstraint::default(), false);
            builder.add_external_edge(a, PointConstraint::default(), false, Flow::Source);
            builder.add_external_edge(c, PointConstraint::default(), false, second_flow);
            let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
            let mut state = graph.new_layout_state(
                vec![
                    Point2::new(-2.0, 0.0),
                    Point2::new(1.0, 2.0),
                    Point2::new(4.0, -1.0),
                ]
                .into(),
                vec![
                    Point2::new(-1.0, 1.0),
                    Point2::new(2.5, 0.5),
                    Point2::new(-4.0, -1.0),
                    Point2::new(6.0, 1.0),
                ]
                .into(),
                1.0,
                0.0,
                true,
            );
            state.edge_spring_length_scales = vec![0.5, 2.0, 1.5, 0.75].into();
            let energy = SpringChargeEnergy {
                external_pull: 2.0,
                ..test_energy(0.2)
            };
            let mut rng = SmallRng::seed_from_u64(17);
            let mut cached = energy.energy(None, &state);
            for iteration in 0..100 {
                let next = LayoutNeighbor.propose(&state, &mut rng, 0.2, 0.3);
                let incremental = energy.energy(Some((&state, cached)), &next);
                let exact = energy.total_energy(&next);
                assert!(
                    (incremental - exact).abs() <= 1e-9 * (1.0 + exact.abs()),
                    "iteration {iteration}: incremental={incremental}, exact={exact}",
                );
                state = next;
                energy.on_accept(&mut state);
                cached = incremental;
                assert_eq!(energy.energy(Some((&state, cached)), &state), cached);
            }
        }
    }

    #[test]
    fn isolated_energy_ignores_complement_geometry() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        let c = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        builder.add_edge(b, c, PointConstraint::default(), false);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let selected_pair = graph.iter_edges().next().unwrap().0;
        let mut selected = graph.empty_subgraph::<SuBitGraph>();
        selected.add(selected_pair);
        let state = graph
            .new_layout_state(
                vec![
                    Point2::new(-1.0, 0.0),
                    Point2::new(1.0, 0.0),
                    Point2::new(4.0, 3.0),
                ]
                .into(),
                vec![Point2::origin(), Point2::new(3.0, 2.0)].into(),
                1.0,
                0.0,
                false,
            )
            .with_active_subgraph(selected);
        let energy = SpringChargeEnergy {
            spring_length: 1.5,
            k_spring: 2.0,
            c_vv: 3.0,
            dangling_charge: 4.0,
            dangling_centroid_charge: 5.0,
            external_pull: 2.0,
            external_pull_balance: 1.0,
            external_pull_attachment: 1.0,
            c_ev: 6.0,
            c_ee_local: 7.0,
            c_center: 8.0,
            crossing_penalty: 9.0,
            eps: 1e-4,
        };
        let baseline = energy.energy(None, &state);
        let mut moved = state.clone();
        moved.vertex_points[NodeIndex(2)] = Point2::new(10_000.0, -20_000.0);
        moved.edge_points[EdgeIndex(1)] = Point2::new(-30_000.0, 40_000.0);
        assert_eq!(energy.energy(None, &moved), baseline);

        let mut next = state.clone();
        apply_vertex_shift(&mut next, NodeIndex(0), Vector2::new(0.25, -0.5));
        let incremental = energy.energy(Some((&state, baseline)), &next);
        let exact = energy.total_energy(&next);
        assert!((incremental - exact).abs() <= 1e-9 * (1.0 + exact.abs()));
    }

    #[test]
    fn one_sided_isolated_edge_has_one_dangling_spring() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let HedgePair::Paired { source, .. } = graph.iter_edges().next().unwrap().0 else {
            panic!("expected paired edge");
        };
        let mut selected = graph.empty_subgraph::<SuBitGraph>();
        selected.add(source);
        let state = graph
            .new_layout_state(
                vec![Point2::origin(), Point2::new(100.0, 0.0)].into(),
                vec![Point2::new(3.0, 0.0)].into(),
                1.0,
                0.0,
                false,
            )
            .with_active_subgraph(selected);
        let energy = test_energy(0.0);

        assert_eq!(
            SpringChargeEnergy::edge_spring_length(&state, EdgeIndex(0), 1.0),
            2.0
        );
        assert_eq!(energy.energy(None, &state), 0.5);
    }

    #[test]
    fn center_term_penalizes_distance_from_origin() {
        let energy = test_energy(2.0);

        assert_eq!(energy.center_term(0.0), 0.0);
        assert!(energy.center_term(2.0) > energy.center_term(1.0));
    }

    #[test]
    fn length_scale_rescales_dimensional_coefficients() {
        let tune = ParamTuning {
            length_scale: 0.5,
            k_spring: 11.0,
            beta: 3.0,
            gamma_dangling: 0.7,
            gamma_dangling_centroid: 1.3,
            external_pull: 2.0,
            external_pull_balance: 1.0,
            external_pull_attachment: 1.0,
            gamma_ev: 0.2,
            gamma_ee: 0.4,
            g_center: 0.05,
            crossing_penalty: 13.0,
            eps: 1e-4,
        };
        let small = SpringChargeEnergy::from_graph(4, 4.0, 4.0, tune);
        assert_eq!(small.external_pull, 22.0);
        let large = SpringChargeEnergy::from_graph(
            4,
            4.0,
            4.0,
            ParamTuning {
                length_scale: 1.0,
                ..tune
            },
        );

        assert_eq!(small.spring_length, 1.0);
        assert_eq!(large.spring_length, 2.0);
        assert_eq!(large.external_pull / small.external_pull, 2.0);
        assert_eq!(large.c_vv / small.c_vv, 8.0);
        assert_eq!(large.c_ev / small.c_ev, 8.0);
        assert_eq!(large.c_ee_local / small.c_ee_local, 8.0);
        assert_eq!(large.dangling_charge / small.dangling_charge, 8.0);
        assert_eq!(
            large.dangling_centroid_charge / small.dangling_centroid_charge,
            8.0
        );
        assert_eq!(large.crossing_penalty / small.crossing_penalty, 4.0);
        assert_eq!(large.eps / small.eps, 2.0);
        assert_eq!(large.c_center, small.c_center);
    }
}
