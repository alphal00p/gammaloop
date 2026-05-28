use cgmath::{EuclideanSpace, InnerSpace, Point2, Point3, Vector2, Vector3, Zero};

use crate::half_edge::{
    involution::{EdgeIndex, EdgeVec},
    layout::spring::{
        apply_edge_shift_with_groups, apply_vertex_shift_with_groups, directional_force_shift,
        Constraint, HasPointConstraint, LayoutState, PointConstraint, Shiftable,
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
    E: Shiftable + HasPointConstraint,
    V: Shiftable + HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    let mut step = cfg.step;
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
    let mut node_z = init_node_z(&mut rng, state.vertex_points.len().0, z_spread);
    let mut edge_z = init_edge_z(&mut rng, state.edge_points.len().0, z_spread);

    for epoch in 0..cfg.epochs {
        let z_spring = cfg.z_spring * cfg.z_spring_growth.powi(epoch as i32);
        for _ in 0..cfg.steps {
            let (mut forces_v, mut forces_e) =
                compute_forces(state, energy, &node_z, &edge_z, z_spring, &workset);
            if state.directional_force != 0.0 {
                for &idx in &workset.movable_nodes {
                    let bias = directional_force_shift(
                        state.graph[idx].point_constraint(),
                        idx.0,
                        state.directional_force,
                    );
                    forces_v[idx] += Vector3::new(bias.x, bias.y, 0.0);
                }
                for &idx in &workset.movable_edges {
                    let bias = directional_force_shift(
                        state.graph[idx].point_constraint(),
                        idx.0,
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
                if apply_vertex_shift_with_groups(state, idx, shift2) {
                    max_move = max_move.max(shift3.magnitude());
                }
            }

            for &idx in &workset.movable_edges {
                let shift3 = clamp_shift3(forces_e[idx] * step, cfg.max_delta);
                let shift2 = Vector2::new(shift3.x, shift3.y);
                edge_z[idx] += shift3.z;
                if apply_edge_shift_with_groups(state, idx, shift2) {
                    max_move = max_move.max(shift3.magnitude());
                }
            }

            if max_move < cfg.early_tol {
                return;
            }
        }

        step *= cfg.cool;
        if step <= 0.0 {
            break;
        }
    }
}

struct ForceWorkSet {
    movable_nodes: Vec<NodeIndex>,
    movable_edges: Vec<EdgeIndex>,
    movable_node: NodeVec<bool>,
    movable_edge: EdgeVec<bool>,
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
        let mut movable_node = NodeVec::with_capacity(n);
        for i in 0..n {
            let idx = NodeIndex(i);
            let movable = can_shift_directly(state.graph[idx].point_constraint(), i);
            if movable {
                movable_nodes.push(idx);
            }
            movable_node.push(movable);
        }

        let mut movable_edges = Vec::new();
        let mut movable_edge = EdgeVec::with_capacity(m);
        for i in 0..m {
            let idx = EdgeIndex(i);
            let movable = can_shift_directly(state.graph[idx].point_constraint(), i);
            if movable {
                movable_edges.push(idx);
            }
            movable_edge.push(movable);
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
            movable_node,
            movable_edge,
            incident_edges,
            dangling_edges,
        }
    }
}

fn can_shift_directly(constraints: &PointConstraint, index: usize) -> bool {
    can_shift_axis(constraints.x, index) || can_shift_axis(constraints.y, index)
}

fn can_shift_axis(constraint: Constraint, index: usize) -> bool {
    match constraint {
        Constraint::Free => true,
        Constraint::Grouped(reference, _) => reference == index,
        Constraint::Fixed => false,
    }
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
    E: Shiftable + HasPointConstraint,
    V: Shiftable + HasPointConstraint,
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
    // we only accumulate forces for vertices that can move.
    for &ni in &workset.movable_nodes {
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
        for &ni in &workset.movable_nodes {
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
            for &ei in &workset.movable_edges {
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
            if !workset.movable_node[ni] && !workset.movable_edge[ei] {
                continue;
            }
            let pe = point3_from_point(state.edge_points[ei], edge_z[ei]);
            let d = pi - pe;
            let dist = d.magnitude();
            if dist <= 1e-9 {
                continue;
            }
            let dir = d / dist;
            let length = edge_spring_length(state, ei, energy.spring_length);
            let fmag = energy.k_spring * (length - dist);
            let f = dir * fmag;
            if workset.movable_node[ni] {
                forces_v[ni] += f;
            }
            if workset.movable_edge[ei] {
                forces_e[ei] -= f;
            }
        }

        for a in 0..edges.len() {
            for b in (a + 1)..edges.len() {
                let ea = edges[a];
                let eb = edges[b];
                if !workset.movable_edge[ea] && !workset.movable_edge[eb] {
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
                if workset.movable_edge[ea] {
                    forces_e[ea] += f;
                }
                if workset.movable_edge[eb] {
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
                if !workset.movable_edge[ei] && !workset.movable_edge[ej] {
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
                if workset.movable_edge[ei] {
                    forces_e[ei] += f;
                }
                if workset.movable_edge[ej] {
                    forces_e[ej] -= f;
                }
            }
        }
    }

    // Center gravity (if enabled).
    if energy.c_center != 0.0 {
        for &ni in &workset.movable_nodes {
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

    (forces_v, forces_e)
}

fn point3_from_point(p: Point2<f64>, z: f64) -> Point3<f64> {
    Point3::new(p.x, p.y, z)
}

fn center_gravity_force(point: Point3<f64>, c_center: f64) -> Vector3<f64> {
    (Point3::origin() - point) * c_center
}

fn edge_spring_length<'a, E, V, H, N>(
    state: &LayoutState<'a, E, V, H, N>,
    edge: EdgeIndex,
    base: f64,
) -> f64
where
    N: NodeStorageOps<NodeData = V> + Clone,
{
    let (_, pair) = &state.graph[&edge];
    match pair {
        crate::half_edge::involution::HedgePair::Unpaired { .. } => base * 2.0,
        _ => base,
    }
}

fn init_node_z(rng: &mut impl Rng, len: usize, spread: f64) -> NodeVec<f64> {
    let mut out = NodeVec::with_capacity(len);
    for _ in 0..len {
        out.push(rng.gen_range(-spread..=spread));
    }
    out
}

fn init_edge_z(rng: &mut impl Rng, len: usize, spread: f64) -> EdgeVec<f64> {
    let mut out = EdgeVec::with_capacity(len);
    for _ in 0..len {
        out.push(rng.gen_range(-spread..=spread));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn center_gravity_force_points_toward_origin() {
        let force = center_gravity_force(Point3::new(2.0, -3.0, 4.0), 0.5);

        assert_eq!(force, Vector3::new(-1.0, 1.5, -2.0));
    }
}
