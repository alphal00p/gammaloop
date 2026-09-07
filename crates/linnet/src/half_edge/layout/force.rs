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
    pub z_spring: f64,
    pub z_spring_growth: f64,
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
    state.synchronize_grouped_coordinates();
    let mut rng = SmallRng::seed_from_u64(cfg.seed);
    let workset = ForceWorkSet::new(state);
    if workset.movable_nodes.is_empty() && workset.movable_edges.is_empty() {
        return;
    }

    // Break perfect symmetry (e.g., all x=0) so forces can separate axes.
    let jitter = 1e-3 * energy.spring_length * cfg.step.abs();
    if jitter > 0.0 {
        apply_initial_jitter(state, &mut rng, jitter, &workset);
    }
    let z_spread = 10.0 * energy.spring_length;
    // Entries without a directly movable planar axis have no z update either,
    // so keep them on the layout plane instead of stranding them at random z.
    let mut node_z = init_node_z(
        &mut rng,
        state.vertex_points.len().0,
        z_spread,
        &workset.movable_nodes,
    );
    let mut edge_z = init_edge_z(
        &mut rng,
        state.edge_points.len().0,
        z_spread,
        &workset.movable_edges,
    );

    // Virtual depth breaks symmetry, but its separation disappears in the drawing.
    // Finish on the exact plane with a fresh cooling schedule so projected
    // overlaps still have enough movement budget to relax.
    for planar in [false, true] {
        if planar {
            for (_, z) in node_z.iter_mut() {
                *z = 0.0;
            }
            for (_, z) in edge_z.iter_mut() {
                *z = 0.0;
            }
        }
        let mut step = cfg.step;
        'phase: for epoch in 0..cfg.epochs {
            let z_spring = if planar {
                0.0
            } else {
                cfg.z_spring * cfg.z_spring_growth.powi(epoch as i32)
            };
            for _ in 0..cfg.steps {
                let (mut forces_v, mut forces_e) =
                    compute_forces(state, energy, &node_z, &edge_z, z_spring, &workset);
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
                    let shift3 = clamp_shift3(forces_v[idx] * step, cfg.max_delta);
                    let shift2 = Vector2::new(shift3.x, shift3.y);
                    node_z[idx] += shift3.z;
                    apply_vertex_shift_with_groups(state, idx, shift2);
                    max_move = max_move.max(shift3.magnitude());
                }

                for &idx in &workset.movable_edges {
                    let shift3 = clamp_shift3(forces_e[idx] * step, cfg.max_delta);
                    let shift2 = Vector2::new(shift3.x, shift3.y);
                    edge_z[idx] += shift3.z;
                    apply_edge_shift_with_groups(state, idx, shift2);
                    max_move = max_move.max(shift3.magnitude());
                }

                if max_move < cfg.early_tol {
                    break 'phase;
                }
            }

            step *= cfg.cool;
            if step <= 0.0 {
                break;
            }
        }
    }
}

struct ForceWorkSet {
    movable_nodes: Vec<NodeIndex>,
    movable_edges: Vec<EdgeIndex>,
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
        let mut force_nodes = Vec::new();
        let mut force_node = NodeVec::with_capacity(n);
        for i in 0..n {
            let idx = NodeIndex(i);
            let constraints = state.graph[idx].point_constraint();
            let movable = can_shift_directly(constraints, LayoutPointIndex::Node(idx));
            if movable {
                movable_nodes.push(idx);
            }
            let receives_force = can_receive_force(constraints, LayoutPointIndex::Node(idx));
            if receives_force {
                force_nodes.push(idx);
            }
            force_node.push(receives_force);
        }

        let mut movable_edges = Vec::new();
        let mut force_edges = Vec::new();
        let mut force_edge = EdgeVec::with_capacity(m);
        for i in 0..m {
            let idx = EdgeIndex(i);
            let constraints = state.graph[idx].point_constraint();
            let movable = can_shift_directly(constraints, LayoutPointIndex::Edge(idx));
            if movable {
                movable_edges.push(idx);
            }
            let receives_force = can_receive_force(constraints, LayoutPointIndex::Edge(idx));
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
    z_spring: f64,
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
    // we only accumulate forces that can reach a planar degree of freedom.
    for &ni in &workset.force_nodes {
        let pi = point3_from_point(state.vertex_points[ni], node_z[ni]);
        for j in 0..n {
            if ni.0 == j {
                continue;
            }
            let nj = NodeIndex(j);
            let pj = point3_from_point(state.vertex_points[nj], node_z[nj]);
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
            let pi = point3_from_point(state.vertex_points[ni], node_z[ni]);
            for e in 0..m {
                let ei = EdgeIndex(e);
                let pe = point3_from_point(state.edge_points[ei], edge_z[ei]);
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
            let pi = point3_from_point(state.vertex_points[ni], node_z[ni]);
            for &ei in &workset.force_edges {
                let pe = point3_from_point(state.edge_points[ei], edge_z[ei]);
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
        let pi = point3_from_point(state.vertex_points[ni], node_z[ni]);
        let edges = &workset.incident_edges[ni];

        for &ei in edges {
            if !workset.force_node[ni] && !workset.force_edge[ei] {
                continue;
            }
            let pe = point3_from_point(state.edge_points[ei], edge_z[ei]);
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
                let pa = point3_from_point(state.edge_points[ea], edge_z[ea]);
                let pb = point3_from_point(state.edge_points[eb], edge_z[eb]);
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
                let pi = point3_from_point(state.edge_points[ei], edge_z[ei]);
                let pj = point3_from_point(state.edge_points[ej], edge_z[ej]);
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
    // introducing translational drift.
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
                point3_from_point(state.vertex_points[ni], node_z[ni]),
                energy.c_center,
            );
        }
    }

    if z_spring != 0.0 {
        for &ni in &workset.movable_nodes {
            forces_v[ni].z += -z_spring * node_z[ni];
        }
        for &ei in &workset.movable_edges {
            forces_e[ei].z += -z_spring * edge_z[ei];
        }
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

fn point3_from_point(p: Point2<f64>, z: f64) -> Point3<f64> {
    Point3::new(p.x, p.y, z)
}

fn center_gravity_force(point: Point3<f64>, c_center: f64) -> Vector3<f64> {
    (Point3::origin() - point) * c_center
}

fn init_node_z(rng: &mut impl Rng, len: usize, spread: f64, movable: &[NodeIndex]) -> NodeVec<f64> {
    let mut out = NodeVec::with_capacity(len);
    for _ in 0..len {
        out.push(0.0);
    }
    for &idx in movable {
        out[idx] = rng.gen_range(-spread..=spread);
    }
    out
}

fn init_edge_z(rng: &mut impl Rng, len: usize, spread: f64, movable: &[EdgeIndex]) -> EdgeVec<f64> {
    let mut out = EdgeVec::with_capacity(len);
    for _ in 0..len {
        out.push(0.0);
    }
    for &idx in movable {
        out[idx] = rng.gen_range(-spread..=spread);
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
    fn center_gravity_force_points_toward_origin() {
        let force = center_gravity_force(Point3::new(2.0, -3.0, 4.0), 0.5);

        assert_eq!(force, Vector3::new(-1.0, 1.5, -2.0));
    }

    #[test]
    fn planar_relaxation_restores_visible_spring_length_even_without_z_spring() {
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
        for z_spring in [0.0, 0.05, 2.0] {
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
                    z_spring,
                    z_spring_growth: 1.0,
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
        let node_z = init_node_z(&mut rng, 3, 10.0, &[NodeIndex(1)]);
        let edge_z = init_edge_z(&mut rng, 3, 10.0, &[EdgeIndex(2)]);

        assert_eq!(node_z[NodeIndex(0)], 0.0);
        assert_ne!(node_z[NodeIndex(1)], 0.0);
        assert_eq!(node_z[NodeIndex(2)], 0.0);
        assert_eq!(edge_z[EdgeIndex(0)], 0.0);
        assert_eq!(edge_z[EdgeIndex(1)], 0.0);
        assert_ne!(edge_z[EdgeIndex(2)], 0.0);
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
        let state = graph.new_layout_state(node_points, edge_points, 1.0, 0.0, false);

        let mut forces_v = NodeVec::new();
        forces_v.push(Vector3::new(1.0, 10.0, 100.0));
        forces_v.push(Vector3::new(2.0, 20.0, 200.0));
        let mut forces_e = EdgeVec::new();
        forces_e.push(Vector3::new(3.0, 30.0, 300.0));

        project_grouped_forces(&state, &mut forces_v, &mut forces_e);

        assert_eq!(forces_v[NodeIndex(0)], Vector3::new(4.0, 10.0, 100.0));
        assert_eq!(forces_v[NodeIndex(1)], Vector3::new(0.0, 0.0, 200.0));
        assert_eq!(forces_e[EdgeIndex(0)], Vector3::new(0.0, 0.0, 300.0));
    }
}
