use cgmath::{EuclideanSpace, InnerSpace, Point2, Point3, Vector2, Vector3, Zero};

use crate::half_edge::{
    involution::{EdgeIndex, EdgeVec},
    layout::spring::{
        apply_edge_shift_with_groups, apply_vertex_shift_with_groups, directional_force_shift,
        Constraint, HasPointConstraint, LayoutPointIndex, LayoutState, PointConstraint,
        SpringChargeEnergy,
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
    assert!(cfg.depth_scale.is_finite() && cfg.depth_scale >= 0.0);
    assert!((0.0..=1.0).contains(&cfg.flattening_end));
    state.synchronize_grouped_coordinates();
    let mut rng = SmallRng::seed_from_u64(cfg.seed);
    let workset = ForceWorkSet::new(state);
    let iterations = cfg.steps.saturating_mul(cfg.epochs);
    let mut step = cfg.step.max(0.0);

    // Break perfect symmetry (e.g., all x=0) so forces can separate axes.
    let jitter = 1e-3 * energy.spring_length * step;
    if iterations > 0 && jitter > 0.0 {
        apply_initial_jitter(state, &mut rng, jitter, &workset);
    }
    let z_spread = 10.0 * energy.spring_length.abs();
    // Unspecified entries without a directly movable planar axis stay on the
    // layout plane instead of being stranded at random z. Explicit depths are
    // independent of XY constraints, and hard pins are never randomized.
    let mut node_z = init_node_z(
        &mut rng,
        &state.vertex_depths,
        z_spread,
        &workset.movable_node_depth,
    );
    let mut edge_z = init_edge_z(
        &mut rng,
        &state.edge_depths,
        z_spread,
        &workset.movable_edge_depth,
    );
    let z_bound = node_z
        .iter()
        .map(|(_, z)| z)
        .chain(edge_z.iter().map(|(_, z)| z))
        .fold(z_spread, |bound, z| bound.max(z.abs()));

    // Virtual depth breaks symmetry, but its separation disappears in the drawing.
    // Collapse smoothly within the original cooling schedule, leaving the remaining
    // iterations to relax projected overlaps on the exact plane without a restart.
    let flattening_finish = cfg.flattening_end * iterations.saturating_sub(1) as f64;
    for iteration in 0..iterations {
        let scale = if iteration as f64 >= flattening_finish {
            0.0
        } else {
            let u = iteration as f64 / flattening_finish;
            cfg.depth_scale * (1.0 - u).powi(2) * (1.0 + 2.0 * u)
        };
        let (mut forces_v, mut forces_e) =
            compute_forces(state, energy, &node_z, &edge_z, scale, &workset);
        if state.directional_force != 0.0 {
            for &idx in &workset.movable_nodes {
                let bias = directional_force_shift(
                    state.graph[idx].point_constraint(),
                    LayoutPointIndex::Node(idx),
                    state.vertex_points[idx],
                    state.directional_force,
                );
                forces_v[idx] += Vector3::new(bias.x, bias.y, 0.0);
            }
            for &idx in &workset.movable_edges {
                let bias = directional_force_shift(
                    state.graph[idx].point_constraint(),
                    LayoutPointIndex::Edge(idx),
                    state.edge_points[idx],
                    state.directional_force,
                );
                forces_e[idx] += Vector3::new(bias.x, bias.y, 0.0);
            }
        }

        let mut max_move: f64 = 0.0;
        for &idx in &workset.movable_nodes {
            let mut shift3 = clamp_shift3(forces_v[idx] * step, cfg.max_delta);
            if scale != 0.0 && workset.movable_node_depth[idx] {
                let z = (node_z[idx] + shift3.z).clamp(-z_bound, z_bound);
                shift3.z = z - node_z[idx];
                node_z[idx] = z;
            }
            apply_vertex_shift_with_groups(state, idx, Vector2::new(shift3.x, shift3.y));
            max_move = max_move.max(shift3.magnitude());
        }
        for &idx in &workset.movable_edges {
            let mut shift3 = clamp_shift3(forces_e[idx] * step, cfg.max_delta);
            if scale != 0.0 && workset.movable_edge_depth[idx] {
                let z = (edge_z[idx] + shift3.z).clamp(-z_bound, z_bound);
                shift3.z = z - edge_z[idx];
                edge_z[idx] = z;
            }
            apply_edge_shift_with_groups(state, idx, Vector2::new(shift3.x, shift3.y));
            max_move = max_move.max(shift3.magnitude());
        }

        // Even a zero/underflowed step must not terminate before depth collapses.
        if scale == 0.0 && (max_move < cfg.early_tol || step <= 0.0) {
            break;
        }
        if (iteration + 1) % cfg.steps == 0 {
            step = (step * cfg.cool).max(0.0);
        }
    }
    for (idx, z) in node_z.iter() {
        state.vertex_depths[idx] = Some(*z);
    }
    for (idx, z) in edge_z.iter() {
        state.edge_depths[idx] = Some(*z);
    }
}

struct ForceWorkSet {
    movable_nodes: Vec<NodeIndex>,
    movable_edges: Vec<EdgeIndex>,
    movable_node_depth: NodeVec<bool>,
    movable_edge_depth: EdgeVec<bool>,
    force_nodes: Vec<NodeIndex>,
    force_edges: Vec<EdgeIndex>,
    force_node: NodeVec<bool>,
    force_edge: EdgeVec<bool>,
    incident_edges: NodeVec<Vec<EdgeIndex>>,
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

        let mut movable_nodes = Vec::new();
        let mut movable_node_depth = NodeVec::with_capacity(n);
        let mut force_nodes = Vec::new();
        let mut force_node = NodeVec::with_capacity(n);
        for i in 0..n {
            let idx = NodeIndex(i);
            let constraints = state.graph[idx].point_constraint();
            let movable = can_shift_directly(constraints, LayoutPointIndex::Node(idx));
            let movable_depth =
                !state.vertex_depth_pins[idx] && (movable || state.vertex_depths[idx].is_some());
            movable_node_depth.push(movable_depth);
            if movable || movable_depth {
                movable_nodes.push(idx);
            }
            let receives_force =
                can_receive_force(constraints, LayoutPointIndex::Node(idx)) || movable_depth;
            if receives_force {
                force_nodes.push(idx);
            }
            force_node.push(receives_force);
        }

        let mut movable_edges = Vec::new();
        let mut movable_edge_depth = EdgeVec::with_capacity(m);
        let mut force_edges = Vec::new();
        let mut force_edge = EdgeVec::with_capacity(m);
        for i in 0..m {
            let idx = EdgeIndex(i);
            let constraints = state.graph[idx].point_constraint();
            let movable = can_shift_directly(constraints, LayoutPointIndex::Edge(idx));
            let movable_depth =
                !state.edge_depth_pins[idx] && (movable || state.edge_depths[idx].is_some());
            movable_edge_depth.push(movable_depth);
            if movable || movable_depth {
                movable_edges.push(idx);
            }
            let receives_force =
                can_receive_force(constraints, LayoutPointIndex::Edge(idx)) || movable_depth;
            if receives_force {
                force_edges.push(idx);
            }
            force_edge.push(receives_force);
        }

        let mut incident_edges = NodeVec::with_capacity(n);
        for i in 0..n {
            let idx = NodeIndex(i);
            incident_edges.push(
                state
                    .graph
                    .iter_crown(idx)
                    .map(|h| state.graph[&h])
                    .collect(),
            );
        }

        let dangling_edges = state.ext.included_iter().map(|h| state.graph[&h]).collect();

        ForceWorkSet {
            movable_nodes,
            movable_edges,
            movable_node_depth,
            movable_edge_depth,
            force_nodes,
            force_edges,
            force_node,
            force_edge,
            incident_edges,
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
}

fn compute_forces<'a, E, V, H, N>(
    state: &LayoutState<'a, E, V, H, N>,
    energy: &SpringChargeEnergy,
    node_z: &NodeVec<f64>,
    edge_z: &EdgeVec<f64>,
    scale: f64,
    workset: &ForceWorkSet,
) -> (NodeVec<Vector3<f64>>, EdgeVec<Vector3<f64>>)
where
    E: HasPointConstraint,
    V: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    let n = state.vertex_points.len().0;
    let m = state.edge_points.len().0;

    let mut forces_v = NodeVec::with_capacity(n);
    for _ in 0..n {
        forces_v.push(Vector3::zero());
    }

    let mut forces_e = EdgeVec::with_capacity(m);
    for _ in 0..m {
        forces_e.push(Vector3::zero());
    }

    // Vertex-vertex repulsion. Fixed vertices still act as static sources, but
    // we only accumulate forces that can reach a planar or raw-depth degree of freedom.
    for &ni in &workset.force_nodes {
        let pi = point3_from_point(state.vertex_points[ni], node_z[ni], scale);
        for j in 0..n {
            if ni.0 == j {
                continue;
            }
            let nj = NodeIndex(j);
            let pj = point3_from_point(state.vertex_points[nj], node_z[nj], scale);
            let d = pi - pj;
            let dist = d.magnitude();
            if dist <= 1e-9 {
                continue;
            }
            let dir = d / dist;
            let fmag = 0.5 * energy.c_vv / (dist + energy.eps).powi(2);
            let f = dir * fmag;
            forces_v[ni] += f;
        }
    }

    // Edge-vertex repulsion.
    if energy.c_ev != 0.0 {
        for &ni in &workset.force_nodes {
            let pi = point3_from_point(state.vertex_points[ni], node_z[ni], scale);
            for e in 0..m {
                let ei = EdgeIndex(e);
                let pe = point3_from_point(state.edge_points[ei], edge_z[ei], scale);
                let d = pi - pe;
                let dist = d.magnitude();
                if dist <= 1e-9 {
                    continue;
                }
                let dir = d / dist;
                let fmag = energy.c_ev / (dist + energy.eps).powi(2);
                let f = dir * fmag;
                forces_v[ni] += f;
            }
        }

        for i in 0..n {
            let ni = NodeIndex(i);
            let pi = point3_from_point(state.vertex_points[ni], node_z[ni], scale);
            for &ei in &workset.force_edges {
                let pe = point3_from_point(state.edge_points[ei], edge_z[ei], scale);
                let d = pi - pe;
                let dist = d.magnitude();
                if dist <= 1e-9 {
                    continue;
                }
                let dir = d / dist;
                let fmag = energy.c_ev / (dist + energy.eps).powi(2);
                let f = dir * fmag;
                forces_e[ei] -= f;
            }
        }
    }

    // Springs and local edge-edge repulsion around nodes.
    for i in 0..n {
        let ni = NodeIndex(i);
        let pi = point3_from_point(state.vertex_points[ni], node_z[ni], scale);
        let edges = &workset.incident_edges[ni];

        for &ei in edges {
            if !workset.force_node[ni] && !workset.force_edge[ei] {
                continue;
            }
            let pe = point3_from_point(state.edge_points[ei], edge_z[ei], scale);
            let d = pi - pe;
            let dist = d.magnitude();
            if dist <= 1e-9 {
                continue;
            }
            let dir = d / dist;
            let length = SpringChargeEnergy::edge_spring_length(state, ei, energy.spring_length);
            let fmag = energy.k_spring * (length - dist);
            let f = dir * fmag;
            if workset.force_node[ni] {
                forces_v[ni] += f;
            }
            if workset.force_edge[ei] {
                forces_e[ei] -= f;
            }
        }

        for a in 0..edges.len() {
            for b in (a + 1)..edges.len() {
                let ea = edges[a];
                let eb = edges[b];
                if !workset.force_edge[ea] && !workset.force_edge[eb] {
                    continue;
                }
                let pa = point3_from_point(state.edge_points[ea], edge_z[ea], scale);
                let pb = point3_from_point(state.edge_points[eb], edge_z[eb], scale);
                let d = pa - pb;
                let dist = d.magnitude();
                if dist <= 1e-9 {
                    continue;
                }
                let dir = d / dist;
                let fmag = energy.c_ee_local / (dist + energy.eps).powi(2);
                let f = dir * fmag;
                if workset.force_edge[ea] {
                    forces_e[ea] += f;
                }
                if workset.force_edge[eb] {
                    forces_e[eb] -= f;
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
                let pi = point3_from_point(state.edge_points[ei], edge_z[ei], scale);
                let pj = point3_from_point(state.edge_points[ej], edge_z[ej], scale);
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

    // Repel each dangling endpoint from the current node centroid. Sharing the
    // opposite reaction over all nodes prevents this internal force from
    // introducing translational drift. This term is intentionally planar, so
    // auxiliary depth cannot reduce its pressure on the visible endpoints.
    if energy.dangling_centroid_charge != 0.0 && n > 0 {
        let centroid = state
            .vertex_points
            .iter()
            .fold(Vector2::zero(), |sum, (_, point)| sum + point.to_vec())
            / n as f64;
        let mut reaction = Vector2::zero();
        for &ei in &workset.dangling_edges {
            let d = state.edge_points[ei].to_vec() - centroid;
            let dist = d.magnitude();
            if dist <= 1e-9 {
                continue;
            }
            let force = d / dist * (energy.dangling_centroid_charge / (dist + energy.eps).powi(2));
            if workset.force_edge[ei] {
                forces_e[ei] += Vector3::new(force.x, force.y, 0.0);
            }
            reaction += force;
        }
        reaction /= n as f64;
        for &ni in &workset.force_nodes {
            forces_v[ni] -= Vector3::new(reaction.x, reaction.y, 0.0);
        }
    }

    // Center gravity (if enabled).
    if energy.c_center != 0.0 {
        for &ni in &workset.force_nodes {
            forces_v[ni] += center_gravity_force(
                point3_from_point(state.vertex_points[ni], node_z[ni], scale),
                energy.c_center,
            );
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

    // A grouped axis is one shared degree of freedom. Its generalized force is
    // the sum of every dependent point's force along that axis.
    project_grouped_forces(state, &mut forces_v, &mut forces_e);

    (forces_v, forces_e)
}

fn project_grouped_forces<'a, E, V, H, N>(
    state: &LayoutState<'a, E, V, H, N>,
    forces_v: &mut NodeVec<Vector3<f64>>,
    forces_e: &mut EdgeVec<Vector3<f64>>,
) where
    E: HasPointConstraint,
    V: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    let mut projected_v = state
        .graph
        .new_nodevec(|node, _, _| Vector3::new(0.0, 0.0, forces_v[node].z));
    let mut projected_e = state
        .graph
        .new_edgevec(|_, edge, _| Vector3::new(0.0, 0.0, forces_e[edge].z));

    let mut add = |target: LayoutPointIndex, x: Option<f64>, y: Option<f64>| {
        let force = match target {
            LayoutPointIndex::Node(index) => &mut projected_v[index],
            LayoutPointIndex::Edge(index) => &mut projected_e[index],
        };
        force.x += x.unwrap_or(0.0);
        force.y += y.unwrap_or(0.0);
    };

    for i in 0..forces_v.len().0 {
        let index = NodeIndex(i);
        let point = LayoutPointIndex::Node(index);
        let constraints = state.graph[index].point_constraint();
        if let Some(target) = constraints.x.force_target(point) {
            add(target, Some(forces_v[index].x), None);
        }
        if let Some(target) = constraints.y.force_target(point) {
            add(target, None, Some(forces_v[index].y));
        }
    }
    for i in 0..forces_e.len().0 {
        let index = EdgeIndex(i);
        let point = LayoutPointIndex::Edge(index);
        let constraints = state.graph[index].point_constraint();
        if let Some(target) = constraints.x.force_target(point) {
            add(target, Some(forces_e[index].x), None);
        }
        if let Some(target) = constraints.y.force_target(point) {
            add(target, None, Some(forces_e[index].y));
        }
    }

    *forces_v = projected_v;
    *forces_e = projected_e;
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
        involution::Flow,
        layout::{
            simulatedanneale::Energy,
            spring::{ParamTuning, ShiftDirection},
        },
        nodestore::DefaultNodeStore,
        HedgeGraph, NoData,
    };

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
            c_ev: 0.0,
            c_ee_local: 0.0,
            c_center: 0.0,
            crossing_penalty: 0.0,
            eps: 1e-4,
        };
        let workset = ForceWorkSet::new(&state);
        let node_z = vec![0.0; 2].into();
        let edge_z = vec![0.0; 3].into();
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
            let (forces_v, forces_e) =
                compute_forces(&state, &energy, &node_z, &edge_z, 0.0, &workset);
            assert_eq!(forces_v, node_forces.to_vec().into());
            assert_eq!(forces_e, edge_forces.to_vec().into());

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
        let energy = SpringChargeEnergy::from_graph(
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
        let node_z = vec![1.0, -2.0].into();
        let edge_z = vec![3.0, 4.0, -5.0].into();
        let baseline_energy = energy.energy(None, &state);
        let baseline_forces = compute_forces(&state, &energy, &node_z, &edge_z, 0.25, &workset);
        assert!(baseline_energy > 0.0);
        state.edge_spring_length_scales = vec![0.5, 2.0, 3.0].into();
        assert_eq!(energy.energy(None, &state), baseline_energy);
        assert_eq!(
            compute_forces(&state, &energy, &node_z, &edge_z, 0.25, &workset),
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
        assert_eq!(workset.force_edges, vec![EdgeIndex(1)]);
        let energy = SpringChargeEnergy {
            spring_length: 1.0,
            k_spring: 1.0,
            c_vv: 0.0,
            dangling_charge: 0.0,
            dangling_centroid_charge: 0.0,
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
            let (forces_v, forces_e) =
                compute_forces(&state, &energy, &node_z, &edge_z, 0.5, &workset);
            let expected_v = Vector3::new(3.6, 0.0, if pinned { 0.0 } else { -2.65 });
            let expected_e = Vector3::new(-3.6, 0.0, if pinned { 0.0 } else { 2.4 });
            assert!((forces_v[node] - expected_v).magnitude() < 1e-12);
            assert!((forces_e[edge] - expected_e).magnitude() < 1e-12);
            if pinned {
                assert!((clamp_shift3(forces_v[node], 0.2).x - 0.2).abs() < 1e-12);
                assert!((clamp_shift3(forces_e[edge], 0.2).x + 0.2).abs() < 1e-12);
            }
            let (planar_v, planar_e) =
                compute_forces(&state, &energy, &node_z, &edge_z, 0.0, &workset);
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
            eps: 0.0,
            ..energy
        };
        let workset = ForceWorkSet::new(&state);
        for scale in [0.0, 0.5, 1.0, 2.0] {
            let (forces_v, forces_e) =
                compute_forces(&state, &energy, &node_z, &edge_z, scale, &workset);
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

        project_grouped_forces(&state, &mut forces_v, &mut forces_e);

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
