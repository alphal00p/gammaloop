use cgmath::{EuclideanSpace, InnerSpace, Point2, Point3, Vector2, Vector3, Zero};

use crate::half_edge::{
    involution::{EdgeIndex, EdgeVec, Flow, Hedge, HedgeVec},
    layout::spring::{
        apply_edge_shift_with_groups, apply_route_point_shift, apply_vertex_shift_with_groups,
        directional_force_shift, Constraint, HasPointConstraint, LayoutPointIndex, LayoutState,
        PointConstraint, SpringChargeEnergy,
    },
    nodestore::NodeStorageOps,
    subgraph::SubSetLike,
    swap::Swap,
    NodeIndex, NodeVec,
};
use rand::{rngs::SmallRng, Rng, SeedableRng};
#[derive(Debug, Clone, Copy)]
pub struct ForceLayoutConfig {
    pub steps: usize,
    pub epochs: usize,
    pub step: f64,
    pub cool: f64,
    pub max_delta: f64,
    pub early_tol: f64,
    pub seed: u64,
    /// Initial multiplier of raw auxiliary depths (normally 1).
    pub depth_scale: f64,
    /// Fraction of the iteration budget used to collapse depth, in 0..=1 (normally 0.5).
    /// Zero starts planar; one still evaluates the final iteration on the exact plane.
    pub flattening_end: f64,
    /// Initial fraction of the final repulsion coefficients, in 0..=1 (normally 1).
    pub initial_repulsion: f64,
    /// Fraction of the iteration budget used to grow repulsion, in 0..=1 (normally 0.7).
    /// Zero starts with final repulsion; one still reaches it on the final iteration.
    pub repulsion_growth: f64,
}

impl ForceLayoutConfig {
    fn repulsion_scale(self, iteration: usize, total_iterations: usize) -> f64 {
        let finish = self.repulsion_growth * total_iterations.saturating_sub(1) as f64;
        if iteration as f64 >= finish || self.initial_repulsion == 1.0 {
            1.0
        } else {
            let u = iteration as f64 / finish;
            self.initial_repulsion + (1.0 - self.initial_repulsion) * u * u * (3.0 - 2.0 * u)
        }
    }
}

pub fn force_directed_layout<'a, E, V, H, N>(
    state: &mut LayoutState<'a, E, V, H, N>,
    energy: &SpringChargeEnergy,
    cfg: ForceLayoutConfig,
) where
    E: HasPointConstraint,
    V: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    let mut session = ForceLayoutSession::new(state.clone(), *energy, cfg);
    session.run_to_end();
    *state = session.into_state();
}

/// A force solver that can be advanced without reinitializing its random state.
///
/// The session owns the layout state while borrowing the graph that state lays out.
/// Higher-level wrappers can use it for synchronous frame iterators.
pub struct ForceLayoutSession<'a, E, V, H, N>
where
    E: HasPointConstraint,
    V: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    state: LayoutState<'a, E, V, H, N>,
    energy: SpringChargeEnergy,
    cfg: ForceLayoutConfig,
    workset: ForceWorkSet,
    forces: ForceBuffers,
    node_z: NodeVec<f64>,
    edge_z: EdgeVec<f64>,
    // Interior depths start on the plane and persist across streamed steps.
    route_z: HedgeVec<Vec<f64>>,
    z_bound: f64,
    iteration: usize,
    step_size: f64,
    flattening_finish: f64,
    total_iterations: usize,
    done: bool,
    last_max_move: f64,
}

#[derive(Debug, Clone)]
pub struct ForceLayoutSnapshot {
    pub vertex_points: NodeVec<Point2<f64>>,
    pub edge_points: EdgeVec<Point2<f64>>,
    pub route_points: HedgeVec<Vec<Point2<f64>>>,
    pub iteration: usize,
    pub done: bool,
    pub max_move: f64,
}

impl<'a, E, V, H, N> ForceLayoutSession<'a, E, V, H, N>
where
    E: HasPointConstraint,
    V: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    pub fn new(
        mut state: LayoutState<'a, E, V, H, N>,
        energy: SpringChargeEnergy,
        cfg: ForceLayoutConfig,
    ) -> Self {
        assert!(cfg.depth_scale.is_finite() && cfg.depth_scale >= 0.0);
        assert!((0.0..=1.0).contains(&cfg.flattening_end));
        assert!((0.0..=1.0).contains(&cfg.initial_repulsion));
        assert!((0.0..=1.0).contains(&cfg.repulsion_growth));
        state.synchronize_grouped_coordinates();
        let mut rng = SmallRng::seed_from_u64(cfg.seed);
        let workset = ForceWorkSet::new(&state);
        let total_iterations = cfg.steps.saturating_mul(cfg.epochs);
        let step_size = cfg.step.max(0.0);

        // Break perfect symmetry (e.g., all x=0) so radial forces can separate axes.
        // Directional forces need no jitter and must leave unconstrained axes unchanged.
        let radial_forces = [
            energy.k_spring,
            energy.c_vv,
            energy.dangling_charge,
            energy.dangling_centroid_charge,
            energy.c_ev,
            energy.c_ee_local,
            energy.c_center,
        ]
        .into_iter()
        .any(|strength| strength != 0.0)
            || energy.external_pull != 0.0 && !state.external_flows_are_mixed();
        let jitter = 1e-3 * energy.spring_length * step_size;
        if total_iterations > 0 && jitter > 0.0 && radial_forces {
            apply_initial_jitter(&mut state, &mut rng, jitter, &workset);
        }
        let z_spread = 10.0 * energy.spring_length.abs();
        // Unspecified entries without a directly movable planar axis stay on the
        // layout plane instead of being stranded at random z. Explicit depths are
        // independent of XY constraints, and hard pins are never randomized.
        let node_z = init_node_z(
            &mut rng,
            &state.vertex_depths,
            z_spread,
            &workset.movable_node_depth,
        );
        let edge_z = init_edge_z(
            &mut rng,
            &state.edge_depths,
            z_spread,
            &workset.movable_edge_depth,
        );
        let route_z = state
            .graph
            .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]);
        let z_bound = workset
            .active_nodes
            .iter()
            .map(|&node| &node_z[node])
            .chain(workset.active_edges.iter().map(|&edge| &edge_z[edge]))
            .fold(z_spread, |bound, z| bound.max(z.abs()));
        let no_active_points = workset.active_nodes.is_empty() && workset.active_edges.is_empty();

        let forces = ForceBuffers::new(&state);
        Self {
            state,
            energy,
            cfg,
            workset,
            forces,
            node_z,
            edge_z,
            route_z,
            z_bound,
            iteration: 0,
            step_size,
            flattening_finish: cfg.flattening_end * total_iterations.saturating_sub(1) as f64,
            total_iterations,
            done: total_iterations == 0 || no_active_points,
            last_max_move: 0.0,
        }
    }

    pub fn step(&mut self, iterations: usize) -> ForceLayoutSnapshot {
        for _ in 0..iterations {
            if self.done {
                break;
            }
            self.advance();
        }
        self.snapshot()
    }

    pub fn run_to_end(&mut self) {
        while !self.done {
            self.advance();
        }
    }

    pub fn snapshot(&self) -> ForceLayoutSnapshot {
        ForceLayoutSnapshot {
            vertex_points: self.state.vertex_points.clone(),
            edge_points: self.state.edge_points.clone(),
            route_points: self.state.route_points.clone(),
            iteration: self.iteration,
            done: self.done,
            max_move: self.last_max_move,
        }
    }

    pub fn into_state(mut self) -> LayoutState<'a, E, V, H, N> {
        for &idx in &self.workset.active_nodes {
            self.state.vertex_depths[idx] = Some(self.node_z[idx]);
        }
        for &idx in &self.workset.active_edges {
            self.state.edge_depths[idx] = Some(self.edge_z[idx]);
        }
        self.state
    }

    fn advance(&mut self) {
        // Virtual depth breaks symmetry, but its separation disappears in the drawing.
        // Collapse smoothly within the original cooling schedule, leaving the remaining
        // iterations to relax projected overlaps on the exact plane without a restart.
        let iteration = self.iteration;
        let scale = if iteration as f64 >= self.flattening_finish {
            0.0
        } else {
            let u = iteration as f64 / self.flattening_finish;
            self.cfg.depth_scale * (1.0 - u).powi(2) * (1.0 + 2.0 * u)
        };
        // Continue from a compact layout to the final repulsion, without changing
        // coordinates, spring lengths, endpoint pull, or the cooling schedule.
        let repulsion_scale = self.cfg.repulsion_scale(iteration, self.total_iterations);
        let mut energy = self.energy;
        energy.c_vv *= repulsion_scale;
        energy.c_ev *= repulsion_scale;
        energy.c_ee_local *= repulsion_scale;
        energy.dangling_charge *= repulsion_scale;
        energy.dangling_centroid_charge *= repulsion_scale;
        self.forces.compute(
            &self.state,
            &energy,
            (&self.node_z, &self.edge_z, &self.route_z),
            scale,
            &self.workset,
        );
        let forces_v = &mut self.forces.vertices;
        let forces_e = &mut self.forces.edges;
        let forces_r = &self.forces.routes;
        if self.state.directional_force != 0.0 {
            for &idx in &self.workset.movable_nodes {
                let bias = directional_force_shift(
                    self.state.graph[idx].point_constraint(),
                    LayoutPointIndex::Node(idx),
                    self.state.vertex_points[idx],
                    self.state.directional_force,
                );
                forces_v[idx] += Vector3::new(bias.x, bias.y, 0.0);
            }
            for &idx in &self.workset.movable_edges {
                let bias = directional_force_shift(
                    self.state.graph[idx].point_constraint(),
                    LayoutPointIndex::Edge(idx),
                    self.state.edge_points[idx],
                    self.state.directional_force,
                );
                forces_e[idx] += Vector3::new(bias.x, bias.y, 0.0);
            }
        }

        let mut max_move = 0.0_f64;
        for &idx in &self.workset.movable_nodes {
            let mut shift3 = clamp_shift3(forces_v[idx] * self.step_size, self.cfg.max_delta);
            if scale != 0.0 && self.workset.movable_node_depth[idx] {
                let z = (self.node_z[idx] + shift3.z).clamp(-self.z_bound, self.z_bound);
                shift3.z = z - self.node_z[idx];
                self.node_z[idx] = z;
            }
            apply_vertex_shift_with_groups(&mut self.state, idx, Vector2::new(shift3.x, shift3.y));
            max_move = max_move.max(shift3.magnitude());
        }
        for &idx in &self.workset.movable_edges {
            let mut shift3 = clamp_shift3(forces_e[idx] * self.step_size, self.cfg.max_delta);
            if scale != 0.0 && self.workset.movable_edge_depth[idx] {
                let z = (self.edge_z[idx] + shift3.z).clamp(-self.z_bound, self.z_bound);
                shift3.z = z - self.edge_z[idx];
                self.edge_z[idx] = z;
            }
            apply_edge_shift_with_groups(&mut self.state, idx, Vector2::new(shift3.x, shift3.y));
            max_move = max_move.max(shift3.magnitude());
        }

        for &(hedge, point) in &self.workset.movable_routes {
            let shift = clamp_shift3(forces_r[hedge][point] * self.step_size, self.cfg.max_delta);
            let mut actual = Vector3::zero();
            if scale != 0.0 && !self.state.edge_depth_pins[self.state.graph[&hedge]] {
                let z = (self.route_z[hedge][point] + shift.z).clamp(-self.z_bound, self.z_bound);
                actual.z = z - self.route_z[hedge][point];
                self.route_z[hedge][point] = z;
            }
            let before = self.state.route_points[hedge][point];
            apply_route_point_shift(
                &mut self.state,
                hedge,
                point,
                Vector2::new(shift.x, shift.y),
            );
            let delta = self.state.route_points[hedge][point] - before;
            actual.x = delta.x;
            actual.y = delta.y;
            max_move = max_move.max(actual.magnitude());
        }

        self.iteration += 1;
        let epoch_end = self.cfg.steps != 0 && self.iteration.is_multiple_of(self.cfg.steps);
        let budget_end = self.iteration >= self.total_iterations;
        let stationary = max_move < self.cfg.early_tol || self.step_size <= 0.0;
        // Continuous forces cannot easily exchange points sharing a constrained
        // line. Try the same discrete coordinate move as annealing once planar,
        // including before an apparent convergence or the final budgeted step.
        let swap_move = if scale == 0.0 && (epoch_end || budget_end || stationary) {
            self.state.improve_constrained_lines(&energy)
        } else {
            0.0
        };
        max_move = max_move.max(swap_move);
        self.last_max_move = max_move;
        // Even a zero/underflowed step must reach the final planar force settings.
        // An accepted exchange gets another chance to settle or improve ordering.
        if scale == 0.0 && repulsion_scale == 1.0 && stationary && swap_move == 0.0 || budget_end {
            self.done = true;
        }
        if epoch_end {
            self.step_size = (self.step_size * self.cfg.cool).max(0.0);
        }
    }
}

struct ForceWorkSet {
    active_nodes: Vec<NodeIndex>,
    active_edges: Vec<EdgeIndex>,
    movable_nodes: Vec<NodeIndex>,
    movable_edges: Vec<EdgeIndex>,
    movable_routes: Vec<(Hedge, usize)>,
    movable_node_depth: NodeVec<bool>,
    movable_edge_depth: EdgeVec<bool>,
    force_nodes: Vec<NodeIndex>,
    force_node: NodeVec<bool>,
    force_edge: EdgeVec<bool>,
    node_targets: NodeVec<[Option<LayoutPointIndex>; 2]>,
    edge_targets: EdgeVec<[Option<LayoutPointIndex>; 2]>,
    incident_hedges: NodeVec<Vec<Hedge>>,
    dangling_edges: Vec<EdgeIndex>,
}

impl ForceWorkSet {
    fn new<'a, E, V, H, N>(state: &LayoutState<'a, E, V, H, N>) -> Self
    where
        E: HasPointConstraint,
        V: HasPointConstraint,
        N: NodeStorageOps<NodeData = V> + Clone,
    {
        let n = state.vertex_points.len().0;
        let m = state.edge_points.len().0;

        let active_nodes = state.active_nodes().collect::<Vec<_>>();
        let active_edges = state.active_edges().collect::<Vec<_>>();
        let mut movable_nodes = Vec::new();
        let mut movable_node_depth = NodeVec::with_capacity(n);
        let mut force_nodes = Vec::new();
        let mut force_node = NodeVec::with_capacity(n);
        let mut node_targets = NodeVec::with_capacity(n);
        for i in 0..n {
            let idx = NodeIndex(i);
            let active = state.node_is_active(idx);
            let constraints = state.graph[idx].point_constraint();
            let point = LayoutPointIndex::Node(idx);
            node_targets.push([
                constraints.x.force_target(point),
                constraints.y.force_target(point),
            ]);
            let movable = active && can_shift_directly(constraints, LayoutPointIndex::Node(idx));
            let movable_depth = active
                && !state.vertex_depth_pins[idx]
                && (movable || state.vertex_depths[idx].is_some());
            movable_node_depth.push(movable_depth);
            if movable || movable_depth {
                movable_nodes.push(idx);
            }
            let receives_force = active
                && (can_receive_force(constraints, LayoutPointIndex::Node(idx)) || movable_depth);
            if receives_force {
                force_nodes.push(idx);
            }
            force_node.push(receives_force);
        }

        let mut movable_edges = Vec::new();
        let mut movable_edge_depth = EdgeVec::with_capacity(m);
        let mut force_edge = EdgeVec::with_capacity(m);
        let mut edge_targets = EdgeVec::with_capacity(m);
        for i in 0..m {
            let idx = EdgeIndex(i);
            let active = state.edge_is_active(idx);
            let constraints = state.graph[idx].point_constraint();
            let point = LayoutPointIndex::Edge(idx);
            edge_targets.push([
                constraints.x.force_target(point),
                constraints.y.force_target(point),
            ]);
            let movable = active && can_shift_directly(constraints, LayoutPointIndex::Edge(idx));
            let movable_depth = active
                && !state.edge_depth_pins[idx]
                && (movable || state.edge_depths[idx].is_some());
            movable_edge_depth.push(movable_depth);
            if movable || movable_depth {
                movable_edges.push(idx);
            }
            let receives_force = active
                && (can_receive_force(constraints, LayoutPointIndex::Edge(idx)) || movable_depth);
            force_edge.push(receives_force);
        }

        let mut incident_hedges = NodeVec::with_capacity(n);
        for i in 0..n {
            let idx = NodeIndex(i);
            incident_hedges.push(
                state
                    .graph
                    .iter_crown(idx)
                    .filter(|&h| state.hedge_is_active(h))
                    .collect(),
            );
        }
        let mut movable_routes = Vec::new();
        for (hedge, points) in &state.route_points {
            let constraints = state.graph[state.graph[&hedge]].point_constraint();
            if state.hedge_is_active(hedge)
                && (!matches!(constraints.x, Constraint::Fixed)
                    || !matches!(constraints.y, Constraint::Fixed))
            {
                movable_routes.extend((0..points.len()).map(|i| (hedge, i)));
            }
        }

        let dangling_edges = state.ext.included_iter().map(|h| state.graph[&h]).collect();

        ForceWorkSet {
            active_nodes,
            active_edges,
            movable_nodes,
            movable_edges,
            movable_routes,
            movable_node_depth,
            movable_edge_depth,
            force_nodes,
            force_node,
            force_edge,
            node_targets,
            edge_targets,
            incident_hedges,
            dangling_edges,
        }
    }
}

fn can_shift_directly(constraints: &PointConstraint, index: LayoutPointIndex) -> bool {
    can_shift_axis(constraints.x, index) || can_shift_axis(constraints.y, index)
}

fn can_shift_axis(constraint: Constraint, index: LayoutPointIndex) -> bool {
    constraint.force_target(index) == Some(index)
}

fn can_receive_force(constraints: &PointConstraint, index: LayoutPointIndex) -> bool {
    constraints.x.force_target(index).is_some() || constraints.y.force_target(index).is_some()
}

fn clamp_shift3(shift: Vector3<f64>, max_delta: f64) -> Vector3<f64> {
    if max_delta <= 0.0 {
        return shift;
    }

    let mag = shift.magnitude();
    if mag > max_delta {
        shift * (max_delta / mag)
    } else {
        shift
    }
}

fn apply_initial_jitter<'a, E, V, H, N>(
    state: &mut LayoutState<'a, E, V, H, N>,
    rng: &mut impl Rng,
    jitter: f64,
    workset: &ForceWorkSet,
) where
    E: HasPointConstraint,
    V: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    for &idx in &workset.movable_nodes {
        let shift = Vector2::new(
            rng.gen_range(-jitter..=jitter),
            rng.gen_range(-jitter..=jitter),
        );
        if shift != Vector2::zero() {
            apply_vertex_shift_with_groups(state, idx, shift);
        }
    }

    for &idx in &workset.movable_edges {
        let shift = Vector2::new(
            rng.gen_range(-jitter..=jitter),
            rng.gen_range(-jitter..=jitter),
        );
        if shift != Vector2::zero() {
            apply_edge_shift_with_groups(state, idx, shift);
        }
    }
    for &(hedge, point) in &workset.movable_routes {
        let shift = Vector2::new(
            rng.gen_range(-jitter..=jitter),
            rng.gen_range(-jitter..=jitter),
        );
        apply_route_point_shift(state, hedge, point, shift);
    }
}

type ProjectedEdgeChargeSample = (Option<(Hedge, usize)>, Point3<f64>, f64);

struct ForceBuffers {
    vertices: NodeVec<Vector3<f64>>,
    edges: EdgeVec<Vector3<f64>>,
    routes: HedgeVec<Vec<Vector3<f64>>>,
    samples: EdgeVec<Vec<ProjectedEdgeChargeSample>>,
    half_routes: HedgeVec<Vec<Point3<f64>>>,
    projected_vertices: NodeVec<Vector3<f64>>,
    projected_edges: EdgeVec<Vector3<f64>>,
    node_points: NodeVec<Point3<f64>>,
    edge_points: EdgeVec<Point3<f64>>,
}

impl ForceBuffers {
    fn new<E, V, H, N>(state: &LayoutState<'_, E, V, H, N>) -> Self
    where
        E: HasPointConstraint,
        V: HasPointConstraint,
        N: NodeStorageOps<NodeData = V> + Clone,
    {
        let nodes = state.vertex_points.len().0;
        let edges = state.edge_points.len().0;
        Self {
            vertices: vec![Vector3::zero(); nodes].into(),
            edges: vec![Vector3::zero(); edges].into(),
            routes: state
                .graph
                .new_hedgevec(|h, _| vec![Vector3::zero(); state.route_points[h].len()]),
            samples: state.graph.new_edgevec(|_, edge, _| {
                state
                    .edge_charge_samples(edge)
                    .into_iter()
                    .map(|(address, _, weight)| (address, Point3::origin(), weight))
                    .collect()
            }),
            half_routes: state
                .graph
                .new_hedgevec(|h, _| vec![Point3::origin(); state.route_points[h].len() + 2]),
            projected_vertices: vec![Vector3::zero(); nodes].into(),
            projected_edges: vec![Vector3::zero(); edges].into(),
            node_points: vec![Point3::origin(); nodes].into(),
            edge_points: vec![Point3::origin(); edges].into(),
        }
    }

    fn compute<'a, E, V, H, N>(
        &mut self,
        state: &LayoutState<'a, E, V, H, N>,
        energy: &SpringChargeEnergy,
        depths: (&NodeVec<f64>, &EdgeVec<f64>, &HedgeVec<Vec<f64>>),
        scale: f64,
        workset: &ForceWorkSet,
    ) where
        E: HasPointConstraint,
        V: HasPointConstraint,
        N: NodeStorageOps<NodeData = V> + Clone,
    {
        let (node_z, edge_z, route_z) = depths;
        for &node in &workset.active_nodes {
            self.node_points[node] =
                point3_from_point(state.vertex_points[node], node_z[node], scale);
        }
        for &edge in &workset.active_edges {
            self.edge_points[edge] =
                point3_from_point(state.edge_points[edge], edge_z[edge], scale);
        }
        // The session keeps topology and route counts fixed. Cache sample weights
        // once, and reuse projected path and force storage at every iteration.
        for &edge in &workset.active_edges {
            for (address, point, _) in &mut self.samples[edge] {
                *point = address.map_or(self.edge_points[edge], |(h, i)| {
                    point3_from_point(state.route_points[h][i], route_z[h][i], scale)
                });
            }
        }
        for &node in &workset.active_nodes {
            for &hedge in &workset.incident_hedges[node] {
                let points = &mut self.half_routes[hedge];
                points[0] = self.node_points[node];
                for (i, &point) in state.route_points[hedge].iter().enumerate() {
                    points[i + 1] = point3_from_point(point, route_z[hedge][i], scale);
                }
                let last = points.len() - 1;
                points[last] = self.edge_points[state.graph[&hedge]];
            }
        }
        let samples = &self.samples;
        let forces_v = &mut self.vertices;
        let forces_e = &mut self.edges;
        let forces_r = &mut self.routes;
        for (_, forces) in forces_r.iter_mut() {
            forces.fill(Vector3::zero());
        }
        for (_, force) in forces_v.iter_mut() {
            *force = Vector3::zero();
        }
        for (_, force) in forces_e.iter_mut() {
            *force = Vector3::zero();
        }

        // Vertex-vertex repulsion. Fixed vertices still act as static sources, but
        // we only accumulate forces that can reach a planar or raw-depth degree of freedom.
        for (i, &ni) in workset.active_nodes.iter().enumerate() {
            let pi = self.node_points[ni];
            for &nj in &workset.active_nodes[i + 1..] {
                if !workset.force_node[ni] && !workset.force_node[nj] {
                    continue;
                }
                let pj = self.node_points[nj];
                let d = pi - pj;
                let dist = d.magnitude();
                if dist <= 1e-9 {
                    continue;
                }
                let dir = d / dist;
                let fmag = 0.5 * energy.c_vv / (dist + energy.eps).powi(2);
                let f = dir * fmag;
                // Lower-index reactions arrive before this node's own pairs,
                // retaining the original neighbor order in every sum.
                if workset.force_node[ni] {
                    forces_v[ni] += f;
                }
                if workset.force_node[nj] {
                    forces_v[nj] -= f;
                }
            }
        }

        // Each physical edge keeps one unit of charge, distributed over its anchor
        // and interior samples. Refining a route therefore does not multiply charge.
        if energy.c_ev != 0.0 {
            for &ni in &workset.active_nodes {
                let pi = self.node_points[ni];
                for &ei in &workset.active_edges {
                    for &(address, pe, weight) in &samples[ei] {
                        let d = pi - pe;
                        let dist = d.magnitude();
                        if dist <= 1e-9 {
                            continue;
                        }
                        let f = d / dist * (weight * energy.c_ev / (dist + energy.eps).powi(2));
                        if workset.force_node[ni] {
                            forces_v[ni] += f;
                        }
                        match address {
                            None if workset.force_edge[ei] => forces_e[ei] -= f,
                            Some((h, i)) => forces_r[h][i] -= f,
                            _ => {}
                        }
                    }
                }
            }
        }

        // One spring per incidence uses the whole half-route length. Its chain
        // gradient acts on every segment without adding springs when points are added.
        for &ni in &workset.active_nodes {
            for &hedge in &workset.incident_hedges[ni] {
                let ei = state.graph[&hedge];
                let points = &self.half_routes[hedge];
                let length: f64 = points.windows(2).map(|p| (p[0] - p[1]).magnitude()).sum();
                let rest = SpringChargeEnergy::edge_spring_length(state, ei, energy.spring_length);
                let tension = energy.k_spring * (rest - length);
                for (i, pair) in points.windows(2).enumerate() {
                    let d = pair[0] - pair[1];
                    let dist = d.magnitude();
                    if dist <= 1e-9 {
                        continue;
                    }
                    let f = d / dist * tension;
                    if i == 0 {
                        if workset.force_node[ni] {
                            forces_v[ni] += f;
                        }
                    } else {
                        forces_r[hedge][i - 1] += f;
                    }
                    if i + 2 == points.len() {
                        if workset.force_edge[ei] {
                            forces_e[ei] -= f;
                        }
                    } else {
                        forces_r[hedge][i] -= f;
                    }
                }
            }
            let hedges = &workset.incident_hedges[ni];
            for a in 0..hedges.len() {
                for b in (a + 1)..hedges.len() {
                    let ea = state.graph[&hedges[a]];
                    let eb = state.graph[&hedges[b]];
                    if ea == eb {
                        continue;
                    }
                    for &(address_a, pa, weight_a) in &samples[ea] {
                        for &(address_b, pb, weight_b) in &samples[eb] {
                            let d = pa - pb;
                            let dist = d.magnitude();
                            if dist <= 1e-9 {
                                continue;
                            }
                            let f = d / dist
                                * (weight_a * weight_b * energy.c_ee_local
                                    / (dist + energy.eps).powi(2));
                            match address_a {
                                None if workset.force_edge[ea] => forces_e[ea] += f,
                                Some((h, i)) => forces_r[h][i] += f,
                                _ => {}
                            }
                            match address_b {
                                None if workset.force_edge[eb] => forces_e[eb] -= f,
                                Some((h, i)) => forces_r[h][i] -= f,
                                _ => {}
                            }
                        }
                    }
                }
            }
        }

        // Dangling edge repulsion.
        if energy.dangling_charge != 0.0 {
            let ext_edges = &workset.dangling_edges;

            for i in 0..ext_edges.len() {
                for j in (i + 1)..ext_edges.len() {
                    let ei = ext_edges[i];
                    let ej = ext_edges[j];
                    if !workset.force_edge[ei] && !workset.force_edge[ej] {
                        continue;
                    }
                    let pi = self.edge_points[ei];
                    let pj = self.edge_points[ej];
                    let d = pi - pj;
                    let dist = d.magnitude();
                    if dist <= 1e-9 {
                        continue;
                    }
                    let dir = d / dist;
                    let fmag = 0.5 * energy.dangling_charge / (dist + energy.eps).powi(2);
                    let f = dir * fmag;
                    if workset.force_edge[ei] {
                        forces_e[ei] += f;
                    }
                    if workset.force_edge[ej] {
                        forces_e[ej] -= f;
                    }
                }
            }
        }

        // Pull mixed external flows horizontally with constant, topology-normalized
        // force. A bottleneck-load baseline prevents narrow connections from being
        // overstretched by their external leg count; shared X groups add the spring
        // demand needed to clear distributed attachments.
        // A single external flow instead pulls radially outward. Neither mode imposes
        // a target extent, so springs along a long chain can each stretch. Sharing the
        // opposite reaction over active nodes prevents translational drift. This term
        // is planar, so auxiliary depth cannot reduce pressure on visible endpoints.
        if (energy.dangling_centroid_charge != 0.0 || energy.external_pull != 0.0)
            && !workset.active_nodes.is_empty()
        {
            let centroid = workset
                .active_nodes
                .iter()
                .fold(Vector2::zero(), |sum, &node| {
                    sum + state.vertex_points[node].to_vec()
                })
                / workset.active_nodes.len() as f64;
            let horizontal = state.external_flows_are_mixed();
            let mut reaction = Vector2::zero();
            for hedge in state.ext.included_iter() {
                let ei = state.graph[&hedge];
                let d = state.edge_points[ei].to_vec() - centroid;
                let dist = d.magnitude();
                let mut force = if dist > 1e-9 && energy.dangling_centroid_charge != 0.0 {
                    d / dist * (energy.dangling_centroid_charge / (dist + energy.eps).powi(2))
                } else {
                    Vector2::zero()
                };
                if horizontal {
                    force.x += energy.external_pull_strength(state.external_pull_topology[ei])
                        * match state.graph.flow(hedge) {
                            Flow::Source => 1.0,
                            Flow::Sink => -1.0,
                        };
                } else if dist > 1e-9 {
                    force +=
                        d / dist * energy.external_pull_strength(state.external_pull_topology[ei]);
                }
                if workset.force_edge[ei] {
                    forces_e[ei] += Vector3::new(force.x, force.y, 0.0);
                }
                reaction += force;
            }
            reaction /= workset.active_nodes.len() as f64;
            for &ni in &workset.force_nodes {
                forces_v[ni] -= Vector3::new(reaction.x, reaction.y, 0.0);
            }
        }

        // Center gravity (if enabled).
        if energy.c_center != 0.0 {
            for &ni in &workset.force_nodes {
                forces_v[ni] += center_gravity_force(self.node_points[ni], energy.c_center);
            }
        }

        // Chain rule for effective z = scale * raw z. Mask pins before the movement
        // clamp so a forbidden depth force cannot consume the planar movement budget.
        for (ni, force) in forces_v.iter_mut() {
            force.z = if scale != 0.0 && workset.movable_node_depth[ni] {
                force.z * scale
            } else {
                0.0
            };
        }
        for (ei, force) in forces_e.iter_mut() {
            force.z = if scale != 0.0 && workset.movable_edge_depth[ei] {
                force.z * scale
            } else {
                0.0
            };
        }

        for (hedge, forces) in forces_r.iter_mut() {
            let edge = state.graph[&hedge];
            let constraints = state.graph[edge].point_constraint();
            let active = state.hedge_is_active(hedge);
            for force in forces {
                if !active || matches!(constraints.x, Constraint::Fixed) {
                    force.x = 0.0;
                }
                if !active || matches!(constraints.y, Constraint::Fixed) {
                    force.y = 0.0;
                }
                force.z = if active
                    && scale != 0.0
                    && !state.edge_depth_pins[edge]
                    && (!matches!(constraints.x, Constraint::Fixed)
                        || !matches!(constraints.y, Constraint::Fixed))
                {
                    force.z * scale
                } else {
                    0.0
                };
            }
        }

        // A grouped axis is one shared degree of freedom. Its generalized force is
        // the sum of every dependent point's force along that axis.
        project_grouped_forces(
            workset,
            forces_v,
            forces_e,
            &mut self.projected_vertices,
            &mut self.projected_edges,
        );
    }
}

fn project_grouped_forces(
    workset: &ForceWorkSet,
    forces_v: &mut NodeVec<Vector3<f64>>,
    forces_e: &mut EdgeVec<Vector3<f64>>,
    projected_v: &mut NodeVec<Vector3<f64>>,
    projected_e: &mut EdgeVec<Vector3<f64>>,
) {
    for (node, force) in projected_v.iter_mut() {
        *force = Vector3::new(0.0, 0.0, forces_v[node].z);
    }
    for (edge, force) in projected_e.iter_mut() {
        *force = Vector3::new(0.0, 0.0, forces_e[edge].z);
    }

    let mut add = |target: LayoutPointIndex, x: Option<f64>, y: Option<f64>| {
        let force = match target {
            LayoutPointIndex::Node(index) => &mut projected_v[index],
            LayoutPointIndex::Edge(index) => &mut projected_e[index],
        };
        force.x += x.unwrap_or(0.0);
        force.y += y.unwrap_or(0.0);
    };

    for (index, targets) in &workset.node_targets {
        if let Some(target) = targets[0] {
            add(target, Some(forces_v[index].x), None);
        }
        if let Some(target) = targets[1] {
            add(target, None, Some(forces_v[index].y));
        }
    }
    for (index, targets) in &workset.edge_targets {
        if let Some(target) = targets[0] {
            add(target, Some(forces_e[index].x), None);
        }
        if let Some(target) = targets[1] {
            add(target, None, Some(forces_e[index].y));
        }
    }

    std::mem::swap(forces_v, projected_v);
    std::mem::swap(forces_e, projected_e);
}

fn point3_from_point(p: Point2<f64>, z: f64, scale: f64) -> Point3<f64> {
    Point3::new(p.x, p.y, if scale == 0.0 { 0.0 } else { scale * z })
}

fn center_gravity_force(point: Point3<f64>, c_center: f64) -> Vector3<f64> {
    (Point3::origin() - point) * c_center
}

fn init_node_z(
    rng: &mut impl Rng,
    depths: &NodeVec<Option<f64>>,
    spread: f64,
    movable: &NodeVec<bool>,
) -> NodeVec<f64> {
    let mut out = NodeVec::with_capacity(depths.len().0);
    for (idx, z) in depths.iter() {
        out.push(z.unwrap_or_else(|| {
            if movable[idx] {
                rng.gen_range(-spread..=spread)
            } else {
                0.0
            }
        }));
    }
    out
}

fn init_edge_z(
    rng: &mut impl Rng,
    depths: &EdgeVec<Option<f64>>,
    spread: f64,
    movable: &EdgeVec<bool>,
) -> EdgeVec<f64> {
    let mut out = EdgeVec::with_capacity(depths.len().0);
    for (idx, z) in depths.iter() {
        out.push(z.unwrap_or_else(|| {
            if movable[idx] {
                rng.gen_range(-spread..=spread)
            } else {
                0.0
            }
        }));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::half_edge::{
        builder::HedgeGraphBuilder,
        involution::{Flow, HedgePair},
        layout::{
            simulatedanneale::Energy,
            spring::{ParamTuning, ShiftDirection},
        },
        nodestore::DefaultNodeStore,
        subgraph::{ModifySubSet, SuBitGraph},
        HedgeGraph, NoData,
    };

    type LayoutForces = (
        NodeVec<Vector3<f64>>,
        EdgeVec<Vector3<f64>>,
        HedgeVec<Vec<Vector3<f64>>>,
    );

    fn compute_forces<'a, E, V, H, N>(
        state: &LayoutState<'a, E, V, H, N>,
        energy: &SpringChargeEnergy,
        node_z: &NodeVec<f64>,
        edge_z: &EdgeVec<f64>,
        route_z: &HedgeVec<Vec<f64>>,
        scale: f64,
        workset: &ForceWorkSet,
    ) -> LayoutForces
    where
        E: HasPointConstraint,
        V: HasPointConstraint,
        N: NodeStorageOps<NodeData = V> + Clone,
    {
        let mut buffers = ForceBuffers::new(state);
        buffers.compute(state, energy, (node_z, edge_z, route_z), scale, workset);
        (buffers.vertices, buffers.edges, buffers.routes)
    }

    #[test]
    fn routed_spring_depth_force_uses_effective_depth_chain_rule() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let node = builder.add_node(PointConstraint::default());
        builder.add_external_edge(node, PointConstraint::default(), false, Flow::Source);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut state = graph.new_layout_state(
            vec![Point2::new(-1.0, 0.0)].into(),
            vec![Point2::new(2.0, 0.0)].into(),
            1.0,
            0.0,
            false,
        );
        state.route_points[Hedge(0)] = vec![Point2::new(0.5, 1.0)];
        let energy = SpringChargeEnergy {
            k_spring: 1.3,
            ..SpringChargeEnergy::from_graph(
                1,
                2.0,
                2.0,
                ParamTuning {
                    beta: 0.0,
                    ..ParamTuning::default()
                },
            )
        };
        let scale = 0.4;
        let raw_depth = 0.3;
        let rest =
            SpringChargeEnergy::edge_spring_length(&state, EdgeIndex(0), energy.spring_length);
        let potential = |depth| {
            let a = point3_from_point(state.vertex_points[node], 2.0, scale);
            let p = point3_from_point(state.route_points[Hedge(0)][0], depth, scale);
            let b = point3_from_point(state.edge_points[EdgeIndex(0)], -1.0, scale);
            0.5 * energy.k_spring * ((a - p).magnitude() + (p - b).magnitude() - rest).powi(2)
        };
        let gradient = (potential(raw_depth + 1e-5) - potential(raw_depth - 1e-5)) / 2e-5;
        let forces = compute_forces(
            &state,
            &energy,
            &vec![2.0].into(),
            &vec![-1.0].into(),
            &vec![vec![raw_depth]].into(),
            scale,
            &ForceWorkSet::new(&state),
        );
        assert!(gradient.abs() > 1e-3);
        assert!((forces.2[Hedge(0)][0].z + gradient).abs() < 1e-8);
        state.edge_depth_pins[EdgeIndex(0)] = true;
        let pinned = compute_forces(
            &state,
            &energy,
            &vec![2.0].into(),
            &vec![-1.0].into(),
            &vec![vec![raw_depth]].into(),
            scale,
            &ForceWorkSet::new(&state),
        );
        assert_eq!(pinned.2[Hedge(0)][0].z, 0.0);
        assert_eq!(pinned.2[Hedge(0)][0].x, forces.2[Hedge(0)][0].x);
        assert_eq!(pinned.2[Hedge(0)][0].y, forces.2[Hedge(0)][0].y);
    }

    #[test]
    fn routed_forces_match_every_planar_energy_gradient() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        builder.add_edge(a, b, PointConstraint::default(), false);
        builder.add_edge(a, a, PointConstraint::default(), false);
        builder.add_external_edge(a, PointConstraint::default(), false, Flow::Sink);
        builder.add_external_edge(b, PointConstraint::default(), false, Flow::Source);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut state = graph.new_layout_state(
            vec![Point2::new(-2.0, 0.5), Point2::new(3.0, -0.4)].into(),
            vec![
                Point2::new(0.2, 2.0),
                Point2::new(0.7, -2.0),
                Point2::new(-4.0, 3.0),
                Point2::new(-5.0, -1.0),
                Point2::new(5.0, 1.0),
            ]
            .into(),
            1.0,
            0.0,
            false,
        );
        for (hedge, points) in state.route_points.iter_mut() {
            let node = graph.node_id(hedge);
            let anchor = state.edge_points[graph[&hedge]];
            let vertex = state.vertex_points[node];
            // Distinct routed self-loop incidences, with unequal sample counts.
            for i in 0..1 + hedge.0 % 2 {
                let t = (i + 1) as f64 / (2 + hedge.0 % 2) as f64;
                points.push(
                    vertex + (anchor - vertex) * t + Vector2::new(0.13 * hedge.0 as f64, 0.27),
                );
            }
        }
        let zero = SpringChargeEnergy {
            spring_length: 1.3,
            k_spring: 0.0,
            c_vv: 0.0,
            dangling_charge: 0.0,
            dangling_centroid_charge: 0.0,
            external_pull: 0.0,
            external_pull_balance: 1.0,
            external_pull_attachment: 1.0,
            c_ev: 0.0,
            c_ee_local: 0.0,
            c_center: 0.0,
            crossing_penalty: 0.0,
            eps: 0.01,
        };
        for (term, energy) in [
            SpringChargeEnergy {
                k_spring: 0.9,
                ..zero
            },
            SpringChargeEnergy { c_vv: 0.8, ..zero },
            SpringChargeEnergy { c_ev: 0.7, ..zero },
            SpringChargeEnergy {
                c_ee_local: 0.6,
                ..zero
            },
            SpringChargeEnergy {
                dangling_charge: 0.5,
                dangling_centroid_charge: 0.4,
                external_pull: 0.3,
                ..zero
            },
        ]
        .into_iter()
        .enumerate()
        {
            let (fv, fe, fr) = compute_forces(
                &state,
                &energy,
                &vec![0.0; 2].into(),
                &vec![0.0; 5].into(),
                &graph.new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                0.0,
                &ForceWorkSet::new(&state),
            );
            let total = fv
                .iter()
                .map(|(_, f)| *f)
                .chain(fe.iter().map(|(_, f)| *f))
                .chain(fr.iter().flat_map(|(_, fs)| fs.iter().copied()))
                .fold(Vector3::zero(), |a, b| a + b);
            assert!(
                total.magnitude() < 1e-10,
                "internal forces must balance: term {term}"
            );
            for axis in 0..2 {
                for (node, force) in &fv {
                    let mut plus = state.clone();
                    let mut minus = state.clone();
                    plus.vertex_points[node][axis] += 1e-5;
                    minus.vertex_points[node][axis] -= 1e-5;
                    let gradient =
                        (energy.energy(None, &plus) - energy.energy(None, &minus)) / 2e-5;
                    assert!(
                        (gradient + force[axis]).abs() < 1e-6,
                        "node {node:?}, term {term}: {gradient} / {force:?}"
                    );
                }
                for (edge, force) in &fe {
                    let mut plus = state.clone();
                    let mut minus = state.clone();
                    plus.edge_points[edge][axis] += 1e-5;
                    minus.edge_points[edge][axis] -= 1e-5;
                    let gradient =
                        (energy.energy(None, &plus) - energy.energy(None, &minus)) / 2e-5;
                    assert!(
                        (gradient + force[axis]).abs() < 1e-6,
                        "edge {edge:?}, term {term}: {gradient} / {force:?}"
                    );
                }
                for (hedge, forces) in &fr {
                    for (i, force) in forces.iter().enumerate() {
                        let mut plus = state.clone();
                        let mut minus = state.clone();
                        plus.route_points[hedge][i][axis] += 1e-5;
                        minus.route_points[hedge][i][axis] -= 1e-5;
                        let gradient =
                            (energy.energy(None, &plus) - energy.energy(None, &minus)) / 2e-5;
                        assert!(
                            (gradient + force[axis]).abs() < 1e-6,
                            "route {hedge:?}/{i}, term {term}: {gradient} / {force:?}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn routed_spring_refinement_keeps_endpoint_forces_and_stiffness() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let node = builder.add_node(PointConstraint::default());
        builder.add_external_edge(node, PointConstraint::default(), false, Flow::Source);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut state = graph.new_layout_state(
            vec![Point2::origin()].into(),
            vec![Point2::new(6.0, 0.0)].into(),
            1.0,
            0.0,
            false,
        );
        let energy = SpringChargeEnergy {
            k_spring: 2.0,
            ..SpringChargeEnergy::from_graph(
                1,
                2.0,
                2.0,
                ParamTuning {
                    beta: 0.0,
                    ..ParamTuning::default()
                },
            )
        };
        let baseline_energy = energy.energy(None, &state);
        let baseline = compute_forces(
            &state,
            &energy,
            &vec![0.0].into(),
            &vec![0.0].into(),
            &vec![vec![]].into(),
            0.0,
            &ForceWorkSet::new(&state),
        );
        for count in [1, 3, 7] {
            state.route_points[Hedge(0)] = (1..=count)
                .map(|i| Point2::new(6.0 * i as f64 / (count + 1) as f64, 0.0))
                .collect();
            let forces = compute_forces(
                &state,
                &energy,
                &vec![0.0].into(),
                &vec![0.0].into(),
                &vec![vec![0.0; count]].into(),
                0.0,
                &ForceWorkSet::new(&state),
            );
            assert_eq!(energy.energy(None, &state), baseline_energy);
            assert_eq!(forces.0, baseline.0);
            assert_eq!(forces.1, baseline.1);
            assert!(forces.2[Hedge(0)].iter().all(|f| *f == Vector3::zero()));
        }
    }

    #[test]
    fn routed_sessions_preserve_fixed_axes_inactive_halves_and_streaming() {
        let rail = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Free,
        };
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, rail, false);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut state = graph.new_layout_state(
            vec![Point2::new(-2.0, 0.0), Point2::new(2.0, 0.0)].into(),
            vec![Point2::new(0.0, 0.3)].into(),
            1.0,
            0.0,
            false,
        );
        state.route_points[Hedge(0)] = vec![Point2::new(-1.0, 1.0), Point2::new(-0.4, 1.5)];
        state.route_points[Hedge(1)] = vec![Point2::new(1.0, -1.0)];
        let mut selection = graph.empty_subgraph::<SuBitGraph>();
        selection.add(Hedge(0));
        let state = state.with_active_subgraph(selection);
        let energy = SpringChargeEnergy::from_graph(2, 2.0, 2.0, ParamTuning::default());
        let cfg = ForceLayoutConfig {
            steps: 13,
            epochs: 2,
            step: 0.02,
            cool: 0.8,
            max_delta: 0.1,
            early_tol: 0.0,
            seed: 41,
            depth_scale: 0.4,
            flattening_end: 0.5,
            initial_repulsion: 0.2,
            repulsion_growth: 0.7,
        };
        let mut streamed = ForceLayoutSession::new(state.clone(), energy, cfg);
        let initial = streamed.step(0);
        assert_eq!(initial.route_points[Hedge(1)], state.route_points[Hedge(1)]);
        assert!(streamed.route_z[Hedge(0)].iter().all(|z| *z == 0.0));
        for _ in 0..5 {
            streamed.step(2);
        }
        assert!(streamed.route_z[Hedge(0)].iter().any(|z| *z != 0.0));
        streamed.run_to_end();
        let mut batch = ForceLayoutSession::new(state.clone(), energy, cfg);
        batch.run_to_end();
        let actual = streamed.snapshot();
        let expected = batch.snapshot();
        assert_eq!(actual.vertex_points, expected.vertex_points);
        assert_eq!(actual.edge_points, expected.edge_points);
        assert_eq!(actual.route_points, expected.route_points);
        assert_eq!(actual.route_points[Hedge(1)], state.route_points[Hedge(1)]);
        assert_eq!(actual.vertex_points[b], state.vertex_points[b]);
        assert_ne!(
            actual.route_points[Hedge(0)],
            initial.route_points[Hedge(0)]
        );
        for (before, after) in state.route_points[Hedge(0)]
            .iter()
            .zip(&actual.route_points[Hedge(0)])
        {
            assert_eq!(before.x, after.x);
        }
        assert_eq!(
            actual.edge_points[EdgeIndex(0)].x,
            state.edge_points[EdgeIndex(0)].x
        );
    }

    #[test]
    fn constrained_line_swaps_precede_stopping_and_preserve_streamed_execution() {
        let fixed = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Fixed,
        };
        let rail = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Free,
        };
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(fixed);
        let b = builder.add_node(fixed);
        builder.add_external_edge(a, rail, false, Flow::Source);
        builder.add_external_edge(b, rail, false, Flow::Source);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut state = graph.new_layout_state(
            vec![Point2::new(-1.0, -2.0), Point2::new(-1.0, 2.0)].into(),
            vec![Point2::new(1.0, 2.0), Point2::new(1.0, -2.0)].into(),
            1.0,
            0.0,
            true,
        );
        state.vertex_depths = vec![Some(3.0), Some(-4.0)].into();
        state.edge_depths = vec![Some(5.0), Some(-6.0)].into();
        let energy = SpringChargeEnergy {
            k_spring: 0.0,
            crossing_penalty: 2.0,
            ..SpringChargeEnergy::from_graph(
                2,
                2.0,
                2.0,
                ParamTuning {
                    beta: 0.0,
                    ..ParamTuning::default()
                },
            )
        };
        let cfg = ForceLayoutConfig {
            steps: 100,
            epochs: 2,
            step: 0.0,
            cool: 1.0,
            max_delta: 0.2,
            early_tol: 1e-6,
            seed: 7,
            depth_scale: 0.0,
            flattening_end: 0.0,
            initial_repulsion: 1.0,
            repulsion_growth: 0.7,
        };
        let initial_energy = energy.energy(None, &state);
        let mut streamed = ForceLayoutSession::new(state.clone(), energy, cfg);
        let unchanged = streamed.step(0);
        assert_eq!(unchanged.edge_points, state.edge_points);
        let reordered = streamed.step(1);
        assert_eq!(reordered.iteration, 1);
        assert_eq!(reordered.max_move, 4.0);
        assert!(!reordered.done);
        assert_eq!(
            reordered.edge_points,
            vec![Point2::new(1.0, -2.0), Point2::new(1.0, 2.0)].into()
        );
        let settled = streamed.step(1);
        assert!(settled.done);
        assert_eq!(settled.iteration, 2);
        assert_eq!(settled.max_move, 0.0);
        let mut batch = ForceLayoutSession::new(state.clone(), energy, cfg);
        batch.run_to_end();
        assert_eq!(batch.snapshot().vertex_points, settled.vertex_points);
        assert_eq!(batch.snapshot().edge_points, settled.edge_points);
        let out = streamed.into_state();
        let mut full = out.clone();
        full.incremental = false;
        assert!(energy.energy(None, &full) < initial_energy);
        assert_eq!(out.vertex_depths, state.vertex_depths);
        assert_eq!(out.edge_depths, state.edge_depths);
        // Final-budget steps still attempt an exchange, while a zero budget
        // remains an identity operation and never initializes a discrete move.
        let mut last = ForceLayoutSession::new(
            state.clone(),
            energy,
            ForceLayoutConfig {
                steps: 1,
                epochs: 1,
                ..cfg
            },
        );
        last.run_to_end();
        assert_eq!(last.snapshot().edge_points, out.edge_points);
        assert_eq!(last.snapshot().max_move, 4.0);
        let mut zero =
            ForceLayoutSession::new(state.clone(), energy, ForceLayoutConfig { steps: 0, ..cfg });
        zero.run_to_end();
        assert_eq!(zero.snapshot().edge_points, state.edge_points);
        assert_eq!(zero.snapshot().iteration, 0);
    }

    #[test]
    fn external_pull_balance_matches_energy_and_is_translation_invariant() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        for _ in 0..2 {
            builder.add_external_edge(a, PointConstraint::default(), false, Flow::Sink);
        }
        for _ in 0..3 {
            builder.add_external_edge(b, PointConstraint::default(), false, Flow::Source);
        }
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let original = graph.new_layout_state(
            vec![Point2::new(-1.0, -2.0), Point2::new(3.0, 2.0)].into(),
            vec![
                Point2::new(1.0, 1.0),
                Point2::new(-4.0, -1.0),
                Point2::new(-3.0, 2.0),
                Point2::new(6.0, -3.0),
                Point2::new(5.0, 1.0),
                Point2::new(4.0, 4.0),
            ]
            .into(),
            1.0,
            0.0,
            true,
        );
        let node_z = vec![0.0; 2].into();
        let edge_z = vec![0.0; 6].into();
        for (balance, incoming, outgoing, reaction) in [
            (0.0, -6.0, 6.0, -3.0),
            (
                0.5,
                -6.0 * 0.5_f64.sqrt(),
                6.0 * (1.0_f64 / 3.0).sqrt(),
                6.0 * 0.5_f64.sqrt() - 9.0 * (1.0_f64 / 3.0).sqrt(),
            ),
            (1.0, -3.0, 2.0, 0.0),
            (2.0, -1.5, 2.0 / 3.0, 0.5),
            (
                2.5,
                -1.5 * 0.5_f64.sqrt(),
                (2.0 / 3.0) * (1.0_f64 / 3.0).sqrt(),
                1.5 * 0.5_f64.sqrt() - (1.0_f64 / 3.0).sqrt(),
            ),
        ] {
            for centroid_charge in [0.0, 1.7] {
                let energy = SpringChargeEnergy {
                    k_spring: 0.0,
                    external_pull: 6.0,
                    external_pull_balance: balance,
                    dangling_centroid_charge: centroid_charge,
                    ..SpringChargeEnergy::from_graph(
                        2,
                        2.0,
                        2.0,
                        ParamTuning {
                            beta: 0.0,
                            ..ParamTuning::default()
                        },
                    )
                };
                let original_energy = energy.energy(None, &original);
                for translation in [Vector2::zero(), Vector2::new(-53.0, 18.0)] {
                    let mut state = original.clone();
                    for (_, point) in state.vertex_points.iter_mut() {
                        *point += translation;
                    }
                    for (_, point) in state.edge_points.iter_mut() {
                        *point += translation;
                    }
                    assert!((energy.energy(None, &state) - original_energy).abs() < 1e-10);
                    let (forces_v, forces_e, _) = compute_forces(
                        &state,
                        &energy,
                        &node_z,
                        &edge_z,
                        &state
                            .graph
                            .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                        0.0,
                        &ForceWorkSet::new(&state),
                    );
                    let total = forces_v
                        .iter()
                        .map(|(_, f)| f)
                        .chain(forces_e.iter().map(|(_, f)| f))
                        .fold(Vector3::zero(), |sum, f| sum + f);
                    assert!(total.magnitude() < 1e-12);
                    if centroid_charge == 0.0 {
                        let expected_force = |force: Vector3<f64>, x| {
                            let expected = Vector3::new(x, 0.0, 0.0);
                            if balance == 0.0 || balance == 1.0 {
                                assert_eq!(force, expected);
                            } else {
                                assert!((force - expected).magnitude() < 1e-12);
                            }
                        };
                        for (_, force) in &forces_v {
                            expected_force(*force, reaction);
                        }
                        for i in 1..3 {
                            expected_force(forces_e[EdgeIndex(i)], incoming);
                        }
                        for i in 3..6 {
                            expected_force(forces_e[EdgeIndex(i)], outgoing);
                        }
                    }
                    for (index, force) in forces_v
                        .iter()
                        .map(|(i, f)| (LayoutPointIndex::Node(i), f))
                        .chain(forces_e.iter().map(|(i, f)| (LayoutPointIndex::Edge(i), f)))
                    {
                        for axis in 0..2 {
                            let mut plus = state.clone();
                            let mut minus = state.clone();
                            for (sample, offset) in [(&mut plus, 1e-4), (&mut minus, -1e-4)] {
                                match index {
                                    LayoutPointIndex::Node(i) => {
                                        sample.vertex_points[i][axis] += offset;
                                        sample.changed_nodes.add(i);
                                    }
                                    LayoutPointIndex::Edge(i) => {
                                        sample.edge_points[i][axis] += offset;
                                        sample.changed_edges.add(i);
                                    }
                                }
                            }
                            let plus_energy = energy.energy(None, &plus);
                            let gradient = (plus_energy - energy.energy(None, &minus)) / 2e-4;
                            assert!((force[axis] + gradient).abs() < 1e-7);
                            let incremental = energy.energy(Some((&state, original_energy)), &plus);
                            assert!((incremental - plus_energy).abs() < 1e-10);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn distributed_attachment_gain_matches_shared_coordinate_energy_gradient() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let nodes = (0..3)
            .map(|_| builder.add_node(PointConstraint::default()))
            .collect::<Vec<_>>();
        for i in 0..2 {
            builder.add_edge(nodes[i], nodes[i + 1], PointConstraint::default(), false);
        }
        for _ in 0..2 {
            builder.add_external_edge(nodes[0], PointConstraint::default(), false, Flow::Sink);
        }
        for owner in [1, 1, 2, 2] {
            builder.add_external_edge(
                nodes[owner],
                PointConstraint {
                    x: Constraint::Grouped(
                        LayoutPointIndex::Edge(EdgeIndex(4)),
                        ShiftDirection::Any,
                    ),
                    y: Constraint::Free,
                },
                false,
                Flow::Source,
            );
        }
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let state = graph.new_layout_state(
            vec![
                Point2::new(0.0, 0.0),
                Point2::new(1.0, 0.0),
                Point2::new(2.0, 0.0),
            ]
            .into(),
            vec![
                Point2::new(0.5, 0.0),
                Point2::new(1.5, 0.0),
                Point2::new(-1.0, -1.0),
                Point2::new(-1.0, 1.0),
                Point2::new(3.0, -2.0),
                Point2::new(3.0, -1.0),
                Point2::new(3.0, 1.0),
                Point2::new(3.0, 2.0),
            ]
            .into(),
            1.0,
            0.0,
            true,
        );
        for attachment in [0.0, 1.0, 16.0] {
            for balance in [0.0, 0.25, 1.0, 2.0] {
                let energy = SpringChargeEnergy {
                    k_spring: 0.0,
                    external_pull: 2.0,
                    external_pull_balance: balance,
                    external_pull_attachment: attachment,
                    ..SpringChargeEnergy::from_graph(
                        3,
                        4.0,
                        4.0,
                        ParamTuning {
                            beta: 0.0,
                            ..ParamTuning::default()
                        },
                    )
                };
                let (forces_v, forces_e, _) = compute_forces(
                    &state,
                    &energy,
                    &vec![0.0; 3].into(),
                    &vec![0.0; 8].into(),
                    &graph.new_hedgevec(|_, _| Vec::new()),
                    0.0,
                    &ForceWorkSet::new(&state),
                );
                let expected = 8.0 * (0.25 + attachment * 0.25).powf(balance);
                assert!((forces_e[EdgeIndex(4)].x - expected).abs() < 1e-12);
                for edge in 5..8 {
                    assert_eq!(forces_e[EdgeIndex(edge)].x, 0.0);
                }
                let total = forces_v
                    .iter()
                    .map(|(_, f)| f.x)
                    .chain(forces_e.iter().map(|(_, f)| f.x))
                    .sum::<f64>();
                assert!(total.abs() < 1e-12);
                let mut plus = state.clone();
                let mut minus = state.clone();
                for edge in 4..8 {
                    plus.edge_points[EdgeIndex(edge)].x += 1e-4;
                    minus.edge_points[EdgeIndex(edge)].x -= 1e-4;
                    plus.changed_edges.add(EdgeIndex(edge));
                }
                let reference = energy.energy(None, &state);
                let plus_energy = energy.energy(None, &plus);
                let gradient = (plus_energy - energy.energy(None, &minus)) / 2e-4;
                assert!((forces_e[EdgeIndex(4)].x + gradient).abs() < 1e-7);
                assert!(
                    (energy.energy(Some((&state, reference)), &plus) - plus_energy).abs() < 1e-10
                );
            }
        }
    }

    #[test]
    fn external_pull_balance_leaves_radial_and_single_leg_series_forces_unchanged() {
        for length in [1, 4] {
            for last_flow in [Flow::Sink, Flow::Source] {
                let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
                let nodes = (0..=length)
                    .map(|_| builder.add_node(PointConstraint::default()))
                    .collect::<Vec<_>>();
                for i in 0..length {
                    builder.add_edge(nodes[i], nodes[i + 1], PointConstraint::default(), false);
                }
                builder.add_external_edge(nodes[0], PointConstraint::default(), false, Flow::Sink);
                builder.add_external_edge(
                    nodes[length],
                    PointConstraint::default(),
                    false,
                    last_flow,
                );
                let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
                let state = graph.new_layout_state(
                    (0..=length)
                        .map(|i| Point2::new(i as f64, 0.0))
                        .collect::<Vec<_>>()
                        .into(),
                    (0..length + 2)
                        .map(|i| Point2::new(i as f64 - 1.0, 2.0))
                        .collect::<Vec<_>>()
                        .into(),
                    1.0,
                    0.0,
                    false,
                );
                let mut energy = SpringChargeEnergy {
                    k_spring: 0.0,
                    external_pull: 5.0,
                    external_pull_balance: 0.0,
                    ..SpringChargeEnergy::from_graph(
                        length + 1,
                        4.0,
                        4.0,
                        ParamTuning {
                            beta: 0.0,
                            ..ParamTuning::default()
                        },
                    )
                };
                let node_z = vec![0.0; length + 1].into();
                let edge_z = vec![0.0; length + 2].into();
                let workset = ForceWorkSet::new(&state);
                let baseline = compute_forces(
                    &state,
                    &energy,
                    &node_z,
                    &edge_z,
                    &state
                        .graph
                        .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                    0.0,
                    &workset,
                );
                let baseline_energy = energy.energy(None, &state);
                for balance in [0.5, 1.0, 2.0, 2.5] {
                    energy.external_pull_balance = balance;
                    for attachment in [0.0, 1.0, 16.0, f64::MAX] {
                        energy.external_pull_attachment = attachment;
                        let forces = compute_forces(
                            &state,
                            &energy,
                            &node_z,
                            &edge_z,
                            &state
                                .graph
                                .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                            0.0,
                            &workset,
                        );
                        assert!((energy.energy(None, &state) - baseline_energy).abs() < 1e-10);
                        for (i, force) in &forces.0 {
                            assert!((*force - baseline.0[i]).magnitude() < 1e-12);
                        }
                        for (i, force) in &forces.1 {
                            assert!((*force - baseline.1[i]).magnitude() < 1e-12);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn normalized_pull_preserves_pins_and_streamed_execution() {
        let fixed = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Fixed,
        };
        let free = PointConstraint::default();
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(fixed);
        let b = builder.add_node(free);
        builder.add_edge(a, b, free, false);
        builder.add_external_edge(a, fixed, false, Flow::Sink);
        builder.add_external_edge(a, free, false, Flow::Sink);
        for _ in 0..3 {
            builder.add_external_edge(b, free, false, Flow::Source);
        }
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let state = graph.new_layout_state(
            vec![Point2::new(-1.0, -2.0), Point2::new(3.0, 2.0)].into(),
            vec![
                Point2::new(1.0, 1.0),
                Point2::new(-4.0, -1.0),
                Point2::new(-3.0, 2.0),
                Point2::new(6.0, -3.0),
                Point2::new(5.0, 1.0),
                Point2::new(4.0, 4.0),
            ]
            .into(),
            1.0,
            0.0,
            false,
        );
        let energy = SpringChargeEnergy::from_graph(
            2,
            4.0,
            4.0,
            ParamTuning {
                external_pull: 1.0,
                ..ParamTuning::default()
            },
        );
        let cfg = ForceLayoutConfig {
            steps: 13,
            epochs: 3,
            step: 0.02,
            cool: 0.85,
            max_delta: 0.1,
            early_tol: 0.0,
            seed: 17,
            depth_scale: 1.0,
            flattening_end: 0.5,
            initial_repulsion: 0.15,
            repulsion_growth: 0.7,
        };
        let mut batch = ForceLayoutSession::new(state.clone(), energy, cfg);
        batch.run_to_end();
        let mut streamed = ForceLayoutSession::new(state.clone(), energy, cfg);
        for amount in [1, 7, 0, 3, 11] {
            streamed.step(amount);
        }
        streamed.run_to_end();
        let batch = batch.into_state();
        let streamed = streamed.into_state();
        assert_eq!(streamed.vertex_points, batch.vertex_points);
        assert_eq!(streamed.edge_points, batch.edge_points);
        assert_eq!(streamed.vertex_depths, batch.vertex_depths);
        assert_eq!(streamed.edge_depths, batch.edge_depths);
        assert_eq!(
            streamed.external_pull_topology,
            state.external_pull_topology
        );
        assert_eq!(streamed.vertex_points[a], state.vertex_points[a]);
        assert_eq!(
            streamed.edge_points[EdgeIndex(1)],
            state.edge_points[EdgeIndex(1)]
        );
        assert_eq!(state.external_pull_topology[EdgeIndex(1)].load_share, 0.5);
        assert_eq!(state.external_pull_topology[EdgeIndex(2)].load_share, 0.5);
        for i in 3..6 {
            assert_eq!(
                state.external_pull_topology[EdgeIndex(i)].load_share,
                1.0 / 3.0
            );
        }
    }

    #[test]
    fn same_direction_external_pull_is_radial_balanced_and_matches_energy() {
        for flow in [Flow::Source, Flow::Sink] {
            let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
            let a = builder.add_node(PointConstraint::default());
            let b = builder.add_node(PointConstraint::default());
            builder.add_external_edge(a, PointConstraint::default(), false, flow);
            builder.add_external_edge(b, PointConstraint::default(), false, flow);
            let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
            let energy = SpringChargeEnergy {
                k_spring: 0.0,
                external_pull: 5.0,
                ..SpringChargeEnergy::from_graph(
                    2,
                    2.0,
                    2.0,
                    ParamTuning {
                        beta: 0.0,
                        ..ParamTuning::default()
                    },
                )
            };
            let mut state = graph.new_layout_state(
                vec![Point2::new(-1.0, -2.0), Point2::new(1.0, 2.0)].into(),
                vec![Point2::new(3.0, 4.0), Point2::new(-12.0, 5.0)].into(),
                1.0,
                0.0,
                false,
            );
            let workset = ForceWorkSet::new(&state);
            let (forces_v, forces_e, _) = compute_forces(
                &state,
                &energy,
                &vec![0.0; 2].into(),
                &vec![0.0; 2].into(),
                &state
                    .graph
                    .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                0.0,
                &workset,
            );
            assert!((forces_e[EdgeIndex(0)] - Vector3::new(3.0, 4.0, 0.0)).magnitude() < 1e-12);
            assert!((forces_e[EdgeIndex(1)].magnitude() - 5.0).abs() < 1e-12);
            let total = forces_v
                .iter()
                .map(|(_, force)| force)
                .chain(forces_e.iter().map(|(_, force)| force))
                .fold(Vector3::zero(), |sum, force| sum + force);
            assert!(total.magnitude() < 1e-12);
            assert_eq!(energy.energy(None, &state), -90.0);
            for (index, force) in forces_v
                .iter()
                .map(|(i, f)| (LayoutPointIndex::Node(i), f))
                .chain(forces_e.iter().map(|(i, f)| (LayoutPointIndex::Edge(i), f)))
            {
                for axis in 0..2 {
                    let mut plus = state.clone();
                    let mut minus = state.clone();
                    for (sample, shift) in [(&mut plus, 1e-5), (&mut minus, -1e-5)] {
                        let point = match index {
                            LayoutPointIndex::Node(i) => &mut sample.vertex_points[i],
                            LayoutPointIndex::Edge(i) => &mut sample.edge_points[i],
                        };
                        point[axis] += shift;
                    }
                    let gradient =
                        (energy.energy(None, &plus) - energy.energy(None, &minus)) / 2e-5;
                    assert!((force[axis] + gradient).abs() < 1e-8);
                }
            }
            // Increasing distance and changing only rest lengths cannot change the pull.
            state.edge_points[EdgeIndex(0)] = Point2::new(30.0, 40.0);
            state.edge_spring_length_scales = vec![10.0, 0.1].into();
            let (_, stretched, _) = compute_forces(
                &state,
                &energy,
                &vec![0.0; 2].into(),
                &vec![0.0; 2].into(),
                &state
                    .graph
                    .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                0.0,
                &workset,
            );
            assert!((stretched[EdgeIndex(0)] - forces_e[EdgeIndex(0)]).magnitude() < 1e-12);
        }
    }

    #[test]
    fn repulsion_continuation_is_smooth_monotone_and_reaches_final_settings() {
        let cfg = ForceLayoutConfig {
            steps: 1001,
            epochs: 1,
            step: 0.1,
            cool: 1.0,
            max_delta: 0.0,
            early_tol: 0.0,
            seed: 7,
            depth_scale: 0.0,
            flattening_end: 0.0,
            initial_repulsion: 0.15,
            repulsion_growth: 0.7,
        };
        assert_eq!(cfg.repulsion_scale(0, 1001), 0.15);
        assert_eq!(cfg.repulsion_scale(700, 1001), 1.0);
        assert_eq!(cfg.repulsion_scale(1000, 1001), 1.0);
        assert!((cfg.repulsion_scale(350, 1001) - 0.575).abs() < 1e-12);
        let mut previous = 0.15;
        for iteration in 0..=1000 {
            let scale = cfg.repulsion_scale(iteration, 1001);
            assert!(scale >= previous && scale <= 1.0);
            previous = scale;
        }
        // Smoothstep has zero slope at both joins to constant settings.
        assert!(cfg.repulsion_scale(1, 1001) - 0.15 < 1e-5);
        assert!(1.0 - cfg.repulsion_scale(699, 1001) < 1e-5);
        for total in [0_usize, 1, 11] {
            assert_eq!(cfg.repulsion_scale(total.saturating_sub(1), total), 1.0);
            for iteration in 0..total {
                assert_eq!(
                    ForceLayoutConfig {
                        initial_repulsion: 1.0,
                        ..cfg
                    }
                    .repulsion_scale(iteration, total),
                    1.0
                );
                assert_eq!(
                    ForceLayoutConfig {
                        repulsion_growth: 0.0,
                        ..cfg
                    }
                    .repulsion_scale(iteration, total),
                    1.0
                );
            }
        }
    }

    #[test]
    fn repulsion_continuation_preserves_session_coordinates_and_default_trajectory() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        builder.add_external_edge(a, PointConstraint::default(), false, Flow::Sink);
        builder.add_external_edge(b, PointConstraint::default(), false, Flow::Source);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let state = graph.new_layout_state(
            vec![Point2::new(-1.0, -0.5), Point2::new(1.0, 0.5)].into(),
            vec![
                Point2::new(0.0, 0.5),
                Point2::new(-2.0, -1.0),
                Point2::new(2.0, 1.0),
            ]
            .into(),
            1.0,
            0.0,
            false,
        );
        let energy = SpringChargeEnergy::from_graph(
            2,
            2.0,
            2.0,
            ParamTuning {
                gamma_dangling_centroid: 1.0,
                external_pull: 0.5,
                ..ParamTuning::default()
            },
        );
        let cfg = ForceLayoutConfig {
            steps: 11,
            epochs: 1,
            step: 0.02,
            cool: 1.0,
            max_delta: 0.1,
            early_tol: 0.0,
            seed: 7,
            depth_scale: 0.0,
            flattening_end: 0.0,
            initial_repulsion: 0.15,
            repulsion_growth: 0.7,
        };
        let mut session = ForceLayoutSession::new(state.clone(), energy, cfg);
        let initial = session.snapshot();
        let unchanged = session.step(0);
        assert_eq!(unchanged.vertex_points, initial.vertex_points);
        assert_eq!(unchanged.edge_points, initial.edge_points);
        // The first continuation step equals the same ordinary solver with all
        // repulsive coefficients reduced, while springs and pulling stay fixed.
        let mut compact = energy;
        compact.c_vv *= cfg.initial_repulsion;
        compact.c_ev *= cfg.initial_repulsion;
        compact.c_ee_local *= cfg.initial_repulsion;
        compact.dangling_charge *= cfg.initial_repulsion;
        compact.dangling_centroid_charge *= cfg.initial_repulsion;
        let mut fixed = ForceLayoutSession::new(
            state.clone(),
            compact,
            ForceLayoutConfig {
                initial_repulsion: 1.0,
                ..cfg
            },
        );
        let first = session.step(1);
        let expected = fixed.step(1);
        assert_eq!(first.vertex_points, expected.vertex_points);
        assert_eq!(first.edge_points, expected.edge_points);
        session.step(3);
        session.run_to_end();
        let mut continuous = ForceLayoutSession::new(state.clone(), energy, cfg);
        continuous.run_to_end();
        assert_eq!(
            session.snapshot().vertex_points,
            continuous.snapshot().vertex_points
        );
        assert_eq!(
            session.snapshot().edge_points,
            continuous.snapshot().edge_points
        );
        for initial_repulsion in [0.0, 0.15, 1.0] {
            let base_cfg = ForceLayoutConfig {
                initial_repulsion,
                repulsion_growth: 0.0,
                ..cfg
            };
            let mut instant = ForceLayoutSession::new(state.clone(), energy, base_cfg);
            let mut default = ForceLayoutSession::new(
                state.clone(),
                energy,
                ForceLayoutConfig {
                    initial_repulsion: 1.0,
                    ..cfg
                },
            );
            instant.run_to_end();
            default.run_to_end();
            assert_eq!(
                instant.snapshot().vertex_points,
                default.snapshot().vertex_points
            );
            assert_eq!(
                instant.snapshot().edge_points,
                default.snapshot().edge_points
            );
        }
        // A frozen or apparently converged compact state is not a final solution.
        for (depth_scale, flattening_end, expected_iterations) in [(0.0, 0.0, 8), (1.0, 0.9, 10)] {
            let mut frozen = ForceLayoutSession::new(
                state.clone(),
                energy,
                ForceLayoutConfig {
                    step: 0.0,
                    early_tol: f64::MAX,
                    depth_scale,
                    flattening_end,
                    ..cfg
                },
            );
            frozen.run_to_end();
            assert_eq!(frozen.snapshot().iteration, expected_iterations);
            assert_eq!(frozen.snapshot().vertex_points, state.vertex_points);
            assert_eq!(frozen.snapshot().edge_points, state.edge_points);
        }
    }

    #[test]
    fn external_pull_is_constant_horizontal_balanced_and_matches_energy() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_external_edge(a, PointConstraint::default(), false, Flow::Source);
        builder.add_external_edge(b, PointConstraint::default(), false, Flow::Sink);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let energy = SpringChargeEnergy {
            k_spring: 0.0,
            ..SpringChargeEnergy::from_graph(
                2,
                2.0,
                4.0,
                ParamTuning {
                    beta: 0.0,
                    external_pull: 3.0,
                    ..ParamTuning::default()
                },
            )
        };
        for translation in [-137.0, 0.0, 213.0] {
            let state = graph.new_layout_state(
                vec![
                    Point2::new(99.0 + translation, -2.0),
                    Point2::new(101.0 + translation, 2.0),
                ]
                .into(),
                vec![
                    Point2::new(100.0 + translation, 0.0),
                    Point2::new(104.0 + translation, 7.0),
                ]
                .into(),
                1.0,
                0.0,
                false,
            );
            assert!((energy.energy(None, &state) - 24.0).abs() < 1e-10);
            let workset = ForceWorkSet::new(&state);
            for scale in [0.0, 1.0] {
                let (forces_v, forces_e, _) = compute_forces(
                    &state,
                    &energy,
                    &vec![100.0, -20.0].into(),
                    &vec![-50.0, 70.0].into(),
                    &state
                        .graph
                        .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                    scale,
                    &workset,
                );
                // One endpoint starts at the centroid; the incoming one starts on the wrong side.
                assert_eq!(forces_e[EdgeIndex(0)], Vector3::new(6.0, 0.0, 0.0));
                assert_eq!(forces_e[EdgeIndex(1)], Vector3::new(-6.0, 0.0, 0.0));
                for (_, force) in &forces_v {
                    assert_eq!(*force, Vector3::zero());
                }
                for i in 0..4 {
                    let mut plus = state.clone();
                    let mut minus = state.clone();
                    let h = 1e-4;
                    let force = if i < 2 {
                        plus.vertex_points[NodeIndex(i)].x += h;
                        minus.vertex_points[NodeIndex(i)].x -= h;
                        forces_v[NodeIndex(i)].x
                    } else {
                        plus.edge_points[EdgeIndex(i - 2)].x += h;
                        minus.edge_points[EdgeIndex(i - 2)].x -= h;
                        forces_e[EdgeIndex(i - 2)].x
                    };
                    let gradient =
                        (energy.energy(None, &plus) - energy.energy(None, &minus)) / (2.0 * h);
                    assert!((force + gradient).abs() < 1e-7);
                }
            }
        }
    }

    #[test]
    fn isolated_half_edge_radial_pull_matches_energy_and_ignores_complement() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let HedgePair::Paired { source, sink } = graph.iter_edges().next().unwrap().0 else {
            panic!("expected paired edge");
        };
        let energy = SpringChargeEnergy {
            k_spring: 0.0,
            ..SpringChargeEnergy::from_graph(
                1,
                1.0,
                1.0,
                ParamTuning {
                    beta: 0.0,
                    external_pull: 3.0,
                    ..ParamTuning::default()
                },
            )
        };
        for hedge in [source, sink] {
            let selected_node = graph.node_id(hedge);
            let mut selected = graph.empty_subgraph::<SuBitGraph>();
            selected.add(hedge);
            let mut positions: NodeVec<_> = vec![Point2::new(1000.0, -2000.0); 2].into();
            positions[selected_node] = Point2::new(10.0, 2.0);
            let state = graph
                .new_layout_state(
                    positions,
                    vec![Point2::new(14.0, 3.0)].into(),
                    1.0,
                    0.0,
                    false,
                )
                .with_active_subgraph(selected);
            let workset = ForceWorkSet::new(&state);
            let (forces_v, forces_e, _) = compute_forces(
                &state,
                &energy,
                &vec![0.0; 2].into(),
                &vec![0.0].into(),
                &state
                    .graph
                    .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                0.0,
                &workset,
            );
            // A single exposed half-edge has no incoming/outgoing partition.
            assert!(
                (forces_e[EdgeIndex(0)]
                    - Vector3::new(12.0 / 17.0_f64.sqrt(), 3.0 / 17.0_f64.sqrt(), 0.0))
                .magnitude()
                    < 1e-12
            );
            assert_eq!(forces_v[selected_node], -forces_e[EdgeIndex(0)]);
            for i in 0..3 {
                let mut plus = state.clone();
                let mut minus = state.clone();
                let h = 1e-4;
                let force = if i < 2 {
                    plus.vertex_points[NodeIndex(i)].x += h;
                    minus.vertex_points[NodeIndex(i)].x -= h;
                    forces_v[NodeIndex(i)].x
                } else {
                    plus.edge_points[EdgeIndex(0)].x += h;
                    minus.edge_points[EdgeIndex(0)].x -= h;
                    forces_e[EdgeIndex(0)].x
                };
                let gradient =
                    (energy.energy(None, &plus) - energy.energy(None, &minus)) / (2.0 * h);
                assert!((force + gradient).abs() < 1e-7);
            }
        }
    }

    #[test]
    fn external_pull_respects_pins_and_preserves_free_y() {
        let fixed = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Fixed,
        };
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let node = builder.add_node(fixed);
        builder.add_external_edge(node, fixed, false, Flow::Sink);
        builder.add_external_edge(node, PointConstraint::default(), false, Flow::Source);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut state = graph.new_layout_state(
            vec![Point2::new(20.0, 30.0)].into(),
            vec![Point2::new(25.0, 35.0), Point2::new(15.0, 32.0)].into(),
            1.0,
            0.0,
            false,
        );
        let energy = SpringChargeEnergy {
            k_spring: 0.0,
            ..SpringChargeEnergy::from_graph(
                1,
                1.0,
                1.0,
                ParamTuning {
                    beta: 0.0,
                    external_pull: 3.0,
                    ..ParamTuning::default()
                },
            )
        };
        force_directed_layout(
            &mut state,
            &energy,
            ForceLayoutConfig {
                steps: 100,
                epochs: 1,
                step: 0.1,
                cool: 1.0,
                max_delta: 0.5,
                early_tol: 0.0,
                seed: 42,
                depth_scale: 0.0,
                flattening_end: 0.0,
                initial_repulsion: 1.0,
                repulsion_growth: 0.7,
            },
        );
        assert_eq!(state.vertex_points[NodeIndex(0)], Point2::new(20.0, 30.0));
        assert_eq!(state.edge_points[EdgeIndex(0)], Point2::new(25.0, 35.0));
        assert!((state.edge_points[EdgeIndex(1)].x - 45.0).abs() < 1e-6);
        assert_eq!(state.edge_points[EdgeIndex(1)].y, 32.0);
    }

    #[test]
    fn edge_spring_length_scales_change_forces_consistently_with_energy() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        builder.add_external_edge(a, PointConstraint::default(), false, Flow::Source);
        builder.add_external_edge(b, PointConstraint::default(), false, Flow::Sink);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut state = graph.new_layout_state(
            vec![Point2::origin(), Point2::new(3.0, 4.0)].into(),
            vec![
                Point2::new(0.0, 4.0),
                Point2::new(4.0, 0.0),
                Point2::new(3.0, 8.0),
            ]
            .into(),
            1.0,
            0.0,
            false,
        );
        let energy = SpringChargeEnergy {
            spring_length: 3.0,
            k_spring: 2.0,
            c_vv: 0.0,
            dangling_charge: 0.0,
            dangling_centroid_charge: 0.0,
            external_pull: 0.0,
            external_pull_balance: 1.0,
            external_pull_attachment: 1.0,
            c_ev: 0.0,
            c_ee_local: 0.0,
            c_center: 0.0,
            crossing_penalty: 0.0,
            eps: 1e-4,
        };
        let workset = ForceWorkSet::new(&state);
        let node_z = vec![0.0; 2].into();
        let edge_z = vec![0.0; 3].into();
        let mut buffers = ForceBuffers::new(&state);
        for (scales, expected_energy, node_forces, edge_forces) in [
            (
                [1.0, 1.0, 1.0],
                9.0,
                [Vector3::new(-4.0, 2.0, 0.0), Vector3::new(0.0, -4.0, 0.0)],
                [
                    Vector3::new(0.0, -2.0, 0.0),
                    Vector3::new(4.0, 0.0, 0.0),
                    Vector3::new(0.0, 4.0, 0.0),
                ],
            ),
            (
                [0.5, 1.0, 1.0],
                16.5,
                [Vector3::new(-4.0, 5.0, 0.0), Vector3::new(-3.0, -4.0, 0.0)],
                [
                    Vector3::new(3.0, -5.0, 0.0),
                    Vector3::new(4.0, 0.0, 0.0),
                    Vector3::new(0.0, 4.0, 0.0),
                ],
            ),
            (
                [0.5, 2.0, 1.0],
                76.5,
                [Vector3::new(-16.0, 5.0, 0.0), Vector3::new(-3.0, -4.0, 0.0)],
                [
                    Vector3::new(3.0, -5.0, 0.0),
                    Vector3::new(16.0, 0.0, 0.0),
                    Vector3::new(0.0, 4.0, 0.0),
                ],
            ),
        ] {
            state.edge_spring_length_scales = scales.to_vec().into();
            assert_eq!(energy.energy(None, &state), expected_energy);
            buffers.compute(
                &state,
                &energy,
                (
                    &node_z,
                    &edge_z,
                    &state
                        .graph
                        .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                ),
                0.0,
                &workset,
            );
            let forces_v = &buffers.vertices;
            let forces_e = &buffers.edges;
            assert_eq!(*forces_v, node_forces.to_vec().into());
            assert_eq!(*forces_e, edge_forces.to_vec().into());

            for (index, force) in forces_v
                .iter()
                .map(|(i, f)| (LayoutPointIndex::Node(i), f))
                .chain(forces_e.iter().map(|(i, f)| (LayoutPointIndex::Edge(i), f)))
            {
                for axis in 0..2 {
                    let h = 1e-5;
                    let mut plus = state.clone();
                    let mut minus = state.clone();
                    for (sample, offset) in [(&mut plus, h), (&mut minus, -h)] {
                        let point = match index {
                            LayoutPointIndex::Node(i) => &mut sample.vertex_points[i],
                            LayoutPointIndex::Edge(i) => &mut sample.edge_points[i],
                        };
                        point[axis] += offset;
                    }
                    let gradient =
                        (energy.energy(None, &plus) - energy.energy(None, &minus)) / (2.0 * h);
                    assert!((force[axis] + gradient).abs() < 1e-8);
                }
            }
        }
    }

    #[test]
    fn edge_spring_length_scales_do_not_change_other_energy_or_forces() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(PointConstraint::default());
        let b = builder.add_node(PointConstraint::default());
        builder.add_edge(a, b, PointConstraint::default(), false);
        builder.add_external_edge(a, PointConstraint::default(), false, Flow::Source);
        builder.add_external_edge(b, PointConstraint::default(), false, Flow::Sink);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut state = graph.new_layout_state(
            vec![Point2::origin(), Point2::new(3.0, 4.0)].into(),
            vec![
                Point2::new(0.0, 4.0),
                Point2::new(4.0, 0.0),
                Point2::new(3.0, 8.0),
            ]
            .into(),
            1.0,
            0.0,
            false,
        );
        let mut energy = SpringChargeEnergy::from_graph(
            2,
            6.0,
            8.0,
            ParamTuning {
                k_spring: 0.0,
                gamma_dangling_centroid: 1.0,
                crossing_penalty: 2.0,
                ..ParamTuning::default()
            },
        );
        let workset = ForceWorkSet::new(&state);
        energy.external_pull = 2.0;
        let node_z = vec![1.0, -2.0].into();
        let edge_z = vec![3.0, 4.0, -5.0].into();
        let baseline_energy = energy.energy(None, &state);
        let baseline_forces = compute_forces(
            &state,
            &energy,
            &node_z,
            &edge_z,
            &state
                .graph
                .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
            0.25,
            &workset,
        );
        assert!(baseline_energy > 0.0);
        state.edge_spring_length_scales = vec![0.5, 2.0, 3.0].into();
        assert_eq!(energy.energy(None, &state), baseline_energy);
        assert_eq!(
            compute_forces(
                &state,
                &energy,
                &node_z,
                &edge_z,
                &state
                    .graph
                    .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                0.25,
                &workset
            ),
            baseline_forces
        );
    }

    #[test]
    fn auxiliary_depth_seeds_and_pins_are_independent_of_xy() {
        let fixed = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Fixed,
        };
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        for _ in 0..3 {
            let node = builder.add_node(fixed);
            builder.add_external_edge(node, fixed, false, Flow::Source);
        }
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut state = graph.new_layout_state(
            vec![Point2::origin(); 3].into(),
            vec![Point2::origin(); 3].into(),
            1.0,
            0.0,
            false,
        );
        assert_eq!(state.vertex_depths, vec![None; 3].into());
        assert_eq!(state.edge_depths, vec![None; 3].into());
        assert_eq!(state.vertex_depth_pins, vec![false; 3].into());
        assert_eq!(state.edge_depth_pins, vec![false; 3].into());
        state.vertex_depths = vec![Some(1000.0), Some(2.0), None].into();
        state.edge_depths = vec![Some(-1000.0), Some(-2.0), None].into();
        state.vertex_depth_pins[NodeIndex(0)] = true;
        state.edge_depth_pins[EdgeIndex(0)] = true;
        let original = state.clone();
        assert_eq!(original.vertex_depths, state.vertex_depths);
        assert_eq!(original.edge_depths, state.edge_depths);
        assert_eq!(original.vertex_depth_pins, state.vertex_depth_pins);
        assert_eq!(original.edge_depth_pins, state.edge_depth_pins);
        let workset = ForceWorkSet::new(&state);
        assert_eq!(workset.movable_nodes, vec![NodeIndex(1)]);
        assert_eq!(workset.movable_edges, vec![EdgeIndex(1)]);
        assert_eq!(workset.force_nodes, vec![NodeIndex(1)]);
        assert!(!workset.force_edge[EdgeIndex(0)]);
        assert!(workset.force_edge[EdgeIndex(1)]);
        assert!(!workset.force_edge[EdgeIndex(2)]);
        let energy = SpringChargeEnergy {
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
            c_center: 0.0,
            crossing_penalty: 0.0,
            eps: 1e-4,
        };
        force_directed_layout(
            &mut state,
            &energy,
            ForceLayoutConfig {
                steps: 4,
                epochs: 2,
                step: 0.1,
                cool: 0.9,
                max_delta: 0.2,
                early_tol: 1e-6,
                seed: 2,
                depth_scale: 1.0,
                flattening_end: 0.5,
                initial_repulsion: 1.0,
                repulsion_growth: 0.7,
            },
        );
        assert_eq!(state.vertex_points, original.vertex_points);
        assert_eq!(state.edge_points, original.edge_points);
        assert_eq!(state.vertex_depths[NodeIndex(0)], Some(1000.0));
        assert_eq!(state.edge_depths[EdgeIndex(0)], Some(-1000.0));
        assert_ne!(state.vertex_depths[NodeIndex(1)], Some(2.0));
        assert_ne!(state.edge_depths[EdgeIndex(1)], Some(-2.0));
        assert_eq!(state.vertex_depths[NodeIndex(2)], Some(0.0));
        assert_eq!(state.edge_depths[EdgeIndex(2)], Some(0.0));
    }

    #[test]
    fn free_auxiliary_depths_are_bounded_by_spread_and_supplied_magnitudes() {
        let fixed = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Fixed,
        };
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        for _ in 0..2 {
            let node = builder.add_node(fixed);
            builder.add_external_edge(node, fixed, false, Flow::Source);
        }
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let energy = SpringChargeEnergy {
            spring_length: 1.0,
            k_spring: 0.0,
            c_vv: 0.0,
            dangling_charge: 0.0,
            dangling_centroid_charge: 0.0,
            external_pull: 0.0,
            external_pull_balance: 1.0,
            external_pull_attachment: 1.0,
            c_ev: 1e8,
            c_ee_local: 0.0,
            c_center: 0.0,
            crossing_penalty: 0.0,
            eps: 1e-4,
        };
        for (seed, bound) in [(1.0, 10.0), (50.0, 50.0)] {
            let mut state = graph.new_layout_state(
                vec![Point2::origin(); 2].into(),
                vec![Point2::origin(); 2].into(),
                1.0,
                0.0,
                false,
            );
            state.vertex_depths = vec![Some(seed), Some(-seed)].into();
            state.edge_depths = vec![Some(seed), Some(-seed)].into();
            force_directed_layout(
                &mut state,
                &energy,
                ForceLayoutConfig {
                    steps: 10,
                    epochs: 1,
                    step: 100.0,
                    cool: 1.0,
                    max_delta: 0.0,
                    early_tol: 0.0,
                    seed: 2,
                    depth_scale: 1.0,
                    flattening_end: 0.5,
                    initial_repulsion: 1.0,
                    repulsion_growth: 0.7,
                },
            );
            assert_eq!(state.vertex_depths, vec![Some(bound), Some(-bound)].into());
            assert_eq!(state.edge_depths, vec![Some(bound), Some(-bound)].into());
            assert_eq!(state.vertex_points, vec![Point2::origin(); 2].into());
            assert_eq!(state.edge_points, vec![Point2::origin(); 2].into());
        }
    }

    #[test]
    fn effective_depth_forces_obey_chain_rule_and_mask_pins_before_clamping() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let node = builder.add_node(PointConstraint::default());
        builder.add_external_edge(node, PointConstraint::default(), false, Flow::Source);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let mut state = graph.new_layout_state(
            vec![Point2::origin()].into(),
            vec![Point2::new(3.0, 0.0)].into(),
            1.0,
            0.0,
            false,
        );
        let energy = SpringChargeEnergy {
            spring_length: 1.0,
            k_spring: 2.0,
            c_vv: 0.0,
            dangling_charge: 0.0,
            dangling_centroid_charge: 0.0,
            external_pull: 0.0,
            external_pull_balance: 1.0,
            external_pull_attachment: 1.0,
            c_ev: 0.0,
            c_ee_local: 0.0,
            c_center: 0.5,
            crossing_penalty: 0.0,
            eps: 1e-4,
        };
        let node_z = vec![2.0].into();
        let edge_z = vec![-6.0].into();
        let edge = EdgeIndex(0);
        for pinned in [false, true] {
            state.vertex_depth_pins[node] = pinned;
            state.edge_depth_pins[edge] = pinned;
            let workset = ForceWorkSet::new(&state);
            let (forces_v, forces_e, _) = compute_forces(
                &state,
                &energy,
                &node_z,
                &edge_z,
                &state
                    .graph
                    .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                0.5,
                &workset,
            );
            let expected_v = Vector3::new(3.6, 0.0, if pinned { 0.0 } else { -2.65 });
            let expected_e = Vector3::new(-3.6, 0.0, if pinned { 0.0 } else { 2.4 });
            assert!((forces_v[node] - expected_v).magnitude() < 1e-12);
            assert!((forces_e[edge] - expected_e).magnitude() < 1e-12);
            if pinned {
                assert!((clamp_shift3(forces_v[node], 0.2).x - 0.2).abs() < 1e-12);
                assert!((clamp_shift3(forces_e[edge], 0.2).x + 0.2).abs() < 1e-12);
            }
            let (planar_v, planar_e, _) = compute_forces(
                &state,
                &energy,
                &node_z,
                &edge_z,
                &state
                    .graph
                    .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                0.0,
                &workset,
            );
            assert_eq!(planar_v[node], Vector3::new(2.0, 0.0, 0.0));
            assert_eq!(planar_e[edge], Vector3::new(-2.0, 0.0, 0.0));
        }
        assert_eq!(point3_from_point(Point2::origin(), 1000.0, 0.0).z, 0.0);

        state.vertex_depth_pins[node] = false;
        state.edge_depth_pins[edge] = false;
        let energy = SpringChargeEnergy {
            k_spring: 0.0,
            c_center: 0.0,
            dangling_centroid_charge: 125.0,
            external_pull: 0.0,
            eps: 0.0,
            ..energy
        };
        let workset = ForceWorkSet::new(&state);
        for scale in [0.0, 0.5, 1.0, 2.0] {
            let (forces_v, forces_e, _) = compute_forces(
                &state,
                &energy,
                &node_z,
                &edge_z,
                &state
                    .graph
                    .new_hedgevec(|h, _| vec![0.0; state.route_points[h].len()]),
                scale,
                &workset,
            );
            assert_eq!(forces_v[node], Vector3::new(-125.0 / 9.0, 0.0, 0.0));
            assert_eq!(forces_e[edge], Vector3::new(125.0 / 9.0, 0.0, 0.0));
        }
    }

    #[test]
    fn depth_schedule_uses_total_budget_without_early_stop_or_cooling_restart() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let node = builder.add_node(PointConstraint::default());
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let energy = SpringChargeEnergy {
            spring_length: 1.0,
            k_spring: 0.0,
            c_vv: 0.0,
            dangling_charge: 0.0,
            dangling_centroid_charge: 0.0,
            external_pull: 0.0,
            external_pull_balance: 1.0,
            external_pull_attachment: 1.0,
            c_ev: 0.0,
            c_ee_local: 0.0,
            c_center: 1.0,
            crossing_penalty: 0.0,
            eps: 1e-4,
        };
        for (steps, epochs, flattening_end, scales) in [
            (0, 5, 0.5, &[][..]),
            (5, 0, 0.5, &[][..]),
            (1, 1, 0.5, &[0.0][..]),
            (1, 5, 0.5, &[1.0, 0.5, 0.0, 0.0, 0.0][..]),
            (5, 1, 0.5, &[1.0, 0.5, 0.0, 0.0, 0.0][..]),
            (1, 5, 1.0, &[1.0, 0.84375, 0.5, 0.15625, 0.0][..]),
            (1, 5, 0.0, &[0.0; 5][..]),
        ] {
            for depth_scale in [0.0, 1.0, 2.0] {
                for cool in [0.0, 0.5, 1.0, -1.0] {
                    for early_tol in [0.0, f64::MAX] {
                        let mut state = graph.new_layout_state(
                            vec![Point2::new(1.0, 2.0)].into(),
                            vec![].into(),
                            1.0,
                            0.0,
                            false,
                        );
                        state.vertex_depths[node] = Some(2.0);
                        let cfg = ForceLayoutConfig {
                            steps,
                            epochs,
                            step: 0.1,
                            cool,
                            max_delta: 0.0,
                            early_tol,
                            seed: 2,
                            depth_scale,
                            flattening_end,
                            initial_repulsion: 1.0,
                            repulsion_growth: 0.7,
                        };
                        let mut expected = state.clone();
                        if !scales.is_empty() {
                            let workset = ForceWorkSet::new(&expected);
                            apply_initial_jitter(
                                &mut expected,
                                &mut SmallRng::seed_from_u64(cfg.seed),
                                1e-3 * energy.spring_length * cfg.step,
                                &workset,
                            );
                        }
                        let mut raw_z = 2.0;
                        let mut step = cfg.step;
                        for (iteration, scale) in scales.iter().enumerate() {
                            let scale = scale * depth_scale;
                            let point = expected.vertex_points[node];
                            expected.vertex_points[node] -= point.to_vec() * step;
                            raw_z -= raw_z * scale * scale * step;
                            if scale == 0.0 && (early_tol > 0.0 || step <= 0.0) {
                                break;
                            }
                            if (iteration + 1) % steps == 0 {
                                step = (step * cool).max(0.0);
                            }
                        }
                        force_directed_layout(&mut state, &energy, cfg);
                        assert!(
                            (state.vertex_points[node] - expected.vertex_points[node]).magnitude()
                                < 1e-12,
                            "{cfg:?}"
                        );
                        assert!(
                            (state.vertex_depths[node].unwrap() - raw_z).abs() < 1e-12,
                            "{cfg:?}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn center_gravity_force_points_toward_origin() {
        let force = center_gravity_force(Point3::new(2.0, -3.0, 4.0), 0.5);

        assert_eq!(force, Vector3::new(-1.0, 1.5, -2.0));
    }

    #[test]
    fn depth_collapse_restores_visible_spring_length() {
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let node = builder.add_node(PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Fixed,
        });
        builder.add_external_edge(
            node,
            PointConstraint {
                x: Constraint::Fixed,
                y: Constraint::Free,
            },
            false,
            Flow::Source,
        );
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();
        let energy = SpringChargeEnergy::from_graph(
            1,
            1.0,
            1.0,
            ParamTuning {
                beta: 0.0,
                ..ParamTuning::default()
            },
        );
        for depth_scale in [0.0, 0.05, 1.0, 2.0] {
            let mut state = graph.new_layout_state(
                vec![Point2::origin()].into(),
                vec![Point2::new(0.0, 0.1)].into(),
                0.2,
                0.0,
                false,
            );
            force_directed_layout(
                &mut state,
                &energy,
                ForceLayoutConfig {
                    steps: 100,
                    epochs: 30,
                    step: 0.1,
                    cool: 0.95,
                    max_delta: 0.2,
                    early_tol: 1e-6,
                    seed: 2,
                    depth_scale,
                    flattening_end: 0.5,
                    initial_repulsion: 1.0,
                    repulsion_growth: 0.7,
                },
            );
            let endpoint = state.edge_points[EdgeIndex(0)];
            assert_eq!(state.vertex_points[node], Point2::origin());
            assert_eq!(endpoint.x, 0.0);
            assert!((endpoint.y.abs() - 2.0 * energy.spring_length).abs() < 1e-4);
        }
    }

    #[test]
    fn constrained_points_start_on_virtual_layout_plane() {
        let mut rng = SmallRng::seed_from_u64(7);
        // Unspecified entries without a directly movable planar axis stay on the
        // layout plane instead of being stranded at random z. Explicit depths are
        // independent of XY constraints, and hard pins are never randomized.
        let node_z = init_node_z(
            &mut rng,
            &vec![None, None, None, Some(30.0), Some(-40.0)].into(),
            10.0,
            &vec![false, true, false, false, true].into(),
        );
        let edge_z = init_edge_z(
            &mut rng,
            &vec![None, None, None, Some(-30.0), Some(40.0)].into(),
            10.0,
            &vec![false, false, true, true, false].into(),
        );

        assert_eq!(node_z[NodeIndex(0)], 0.0);
        assert_ne!(node_z[NodeIndex(1)], 0.0);
        assert_eq!(node_z[NodeIndex(2)], 0.0);
        assert_eq!(edge_z[EdgeIndex(0)], 0.0);
        assert_eq!(edge_z[EdgeIndex(1)], 0.0);
        assert_ne!(edge_z[EdgeIndex(2)], 0.0);
        assert_eq!(node_z[NodeIndex(3)], 30.0);
        assert_eq!(node_z[NodeIndex(4)], -40.0);
        assert_eq!(edge_z[EdgeIndex(3)], -30.0);
        assert_eq!(edge_z[EdgeIndex(4)], 40.0);
    }

    #[test]
    fn cross_kind_group_forces_are_summed_at_reference() {
        let reference = LayoutPointIndex::Node(NodeIndex(0));
        let grouped_node = PointConstraint {
            x: Constraint::Grouped(reference, ShiftDirection::Any),
            y: Constraint::Free,
        };
        let grouped_edge = PointConstraint {
            x: Constraint::Grouped(reference, ShiftDirection::Any),
            y: Constraint::Fixed,
        };
        let fixed = PointConstraint {
            x: Constraint::Fixed,
            y: Constraint::Fixed,
        };
        let mut builder = HedgeGraphBuilder::<PointConstraint, PointConstraint>::new();
        let a = builder.add_node(grouped_node);
        let b = builder.add_node(fixed);
        builder.add_edge(a, b, grouped_edge, false);
        let graph: HedgeGraph<_, _, NoData, DefaultNodeStore<_>> = builder.build();

        assert!(can_shift_directly(&grouped_node, reference));
        assert!(!can_shift_directly(
            &grouped_edge,
            LayoutPointIndex::Edge(EdgeIndex(0)),
        ));
        assert!(can_receive_force(
            &grouped_edge,
            LayoutPointIndex::Edge(EdgeIndex(0)),
        ));

        let mut node_points = NodeVec::new();
        node_points.push(Point2::origin());
        node_points.push(Point2::origin());
        let mut edge_points = EdgeVec::new();
        edge_points.push(Point2::origin());
        let mut state = graph.new_layout_state(node_points, edge_points, 1.0, 0.0, false);

        let mut forces_v = NodeVec::new();
        forces_v.push(Vector3::new(1.0, 10.0, 100.0));
        forces_v.push(Vector3::new(2.0, 20.0, 200.0));
        let mut forces_e = EdgeVec::new();
        forces_e.push(Vector3::new(3.0, 30.0, 300.0));

        let mut buffers = ForceBuffers::new(&state);
        project_grouped_forces(
            &ForceWorkSet::new(&state),
            &mut forces_v,
            &mut forces_e,
            &mut buffers.projected_vertices,
            &mut buffers.projected_edges,
        );

        assert_eq!(forces_v[NodeIndex(0)], Vector3::new(4.0, 10.0, 100.0));
        assert_eq!(forces_v[NodeIndex(1)], Vector3::new(0.0, 0.0, 200.0));
        assert_eq!(forces_e[EdgeIndex(0)], Vector3::new(0.0, 0.0, 300.0));

        state.vertex_points[b] = Point2::new(4.0, 0.0);
        state.edge_points[EdgeIndex(0)] = Point2::new(0.0, 2.0);
        state.vertex_depths = vec![Some(3.0), Some(-3.0)].into();
        state.vertex_depth_pins = vec![true; 2].into();
        state.edge_depths[EdgeIndex(0)] = Some(5.0);
        let workset = ForceWorkSet::new(&state);
        assert_eq!(workset.movable_edges, vec![EdgeIndex(0)]);
        assert!(workset.movable_edge_depth[EdgeIndex(0)]);
        let energy = SpringChargeEnergy::from_graph(
            2,
            1.0,
            1.0,
            ParamTuning {
                beta: 0.0,
                ..ParamTuning::default()
            },
        );
        force_directed_layout(
            &mut state,
            &energy,
            ForceLayoutConfig {
                steps: 2,
                epochs: 1,
                step: 0.1,
                cool: 1.0,
                max_delta: 0.2,
                early_tol: 0.0,
                seed: 2,
                depth_scale: 1.0,
                flattening_end: 0.5,
                initial_repulsion: 1.0,
                repulsion_growth: 0.7,
            },
        );
        assert_eq!(state.vertex_points[a].x, state.edge_points[EdgeIndex(0)].x);
        assert_ne!(state.vertex_points[a].x, 0.0);
        assert_eq!(state.vertex_points[b], Point2::new(4.0, 0.0));
        assert_eq!(state.edge_points[EdgeIndex(0)].y, 2.0);
        assert_eq!(state.vertex_depths, vec![Some(3.0), Some(-3.0)].into());
        assert_ne!(state.edge_depths[EdgeIndex(0)], Some(5.0));
    }
}
