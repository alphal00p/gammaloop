use cgmath::{EuclideanSpace, InnerSpace, Point2, Point3, Vector2, Vector3, Zero};

use crate::half_edge::{
    involution::{EdgeIndex, EdgeVec},
    layout::spring::{
        apply_edge_shift_with_groups, apply_vertex_shift_with_groups, directional_force_shift,
        HasPointConstraint, LayoutState, Shiftable, SpringChargeEnergy,
    },
    nodestore::NodeStorageOps,
    subgraph::SubSetLike,
    swap::Swap,
    NodeIndex, NodeVec,
};
use cgmath::MetricSpace;
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
    // Break perfect symmetry (e.g., all x=0) so forces can separate axes.
    let jitter = 1e-3 * energy.spring_length * cfg.step.abs();
    if jitter > 0.0 {
        apply_initial_jitter(state, &mut rng, jitter);
    }
    let z_spread = 10.0 * energy.spring_length;
    let mut node_z = init_node_z(&mut rng, state.vertex_points.len().0, z_spread);
    let mut edge_z = init_edge_z(&mut rng, state.edge_points.len().0, z_spread);

    for epoch in 0..cfg.epochs {
        let z_spring = cfg.z_spring * cfg.z_spring_growth.powi(epoch as i32);
        for _ in 0..cfg.steps {
            let (mut forces_v, mut forces_e) =
                compute_forces(state, energy, &node_z, &edge_z, z_spring);
            if state.directional_force != 0.0 {
                for i in 0..state.vertex_points.len().0 {
                    let idx = NodeIndex(i);
                    let bias = directional_force_shift(
                        state.graph[idx].point_constraint(),
                        idx.0,
                        state.directional_force,
                    );
                    forces_v[idx] += Vector3::new(bias.x, bias.y, 0.0);
                }
                for i in 0..state.edge_points.len().0 {
                    let idx = EdgeIndex(i);
                    let bias = directional_force_shift(
                        state.graph[idx].point_constraint(),
                        idx.0,
                        state.directional_force,
                    );
                    forces_e[idx] += Vector3::new(bias.x, bias.y, 0.0);
                }
            }

            let mut max_move: f64 = 0.0;

            for i in 0..state.vertex_points.len().0 {
                let idx = NodeIndex(i);
                let shift3 = clamp_shift3(forces_v[idx] * step, cfg.max_delta);
                let shift2 = Vector2::new(shift3.x, shift3.y);
                node_z[idx] += shift3.z;
                if apply_vertex_shift_with_groups(state, idx, shift2) {
                    max_move = max_move.max(shift3.magnitude());
                }
            }

            for i in 0..state.edge_points.len().0 {
                let idx = EdgeIndex(i);
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
) where
    E: Shiftable + HasPointConstraint,
    V: Shiftable + HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    for i in 0..state.vertex_points.len().0 {
        let idx = NodeIndex(i);
        let shift = Vector2::new(
            rng.gen_range(-jitter..=jitter),
            rng.gen_range(-jitter..=jitter),
        );
        if shift != Vector2::zero() {
            apply_vertex_shift_with_groups(state, idx, shift);
        }
    }

    for i in 0..state.edge_points.len().0 {
        let idx = EdgeIndex(i);
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

    // Vertex-vertex repulsion.
    for i in 0..n {
        let ni = NodeIndex(i);
        let pi = point3_from_point(state.vertex_points[ni], node_z[ni]);
        for j in (i + 1)..n {
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
            forces_v[nj] -= f;
        }
    }

    // Edge-vertex repulsion.
    if energy.c_ev != 0.0 {
        for i in 0..n {
            let ni = NodeIndex(i);
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
                forces_e[ei] -= f;
            }
        }
    }

    // Springs and local edge-edge repulsion around nodes.
    for i in 0..n {
        let ni = NodeIndex(i);
        let pi = point3_from_point(state.vertex_points[ni], node_z[ni]);
        let edges: Vec<EdgeIndex> = state
            .graph
            .iter_crown(ni)
            .map(|h| state.graph[&h])
            .collect();

        for &ei in &edges {
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
            forces_v[ni] += f;
            forces_e[ei] -= f;
        }

        for a in 0..edges.len() {
            for b in (a + 1)..edges.len() {
                let ea = edges[a];
                let eb = edges[b];
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
                forces_e[ea] += f;
                forces_e[eb] -= f;
            }
        }
    }

    // Dangling edge repulsion.
    if energy.dangling_charge != 0.0 {
        let ext_edges: Vec<EdgeIndex> =
            state.ext.included_iter().map(|h| state.graph[&h]).collect();

        for i in 0..ext_edges.len() {
            for j in (i + 1)..ext_edges.len() {
                let ei = ext_edges[i];
                let ej = ext_edges[j];
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
                forces_e[ei] += f;
                forces_e[ej] -= f;
            }
        }
    }

    // Center repulsion (if enabled).
    if energy.c_center != 0.0 {
        for i in 0..n {
            let ni = NodeIndex(i);
            let pi = point3_from_point(state.vertex_points[ni], node_z[ni]);
            let r = pi.distance(EuclideanSpace::origin());
            if r <= 1.0 {
                continue;
            }
            let dir = (pi - Point3::origin()) / r;
            let fmag = 0.5 * energy.c_center / (r + energy.eps).powi(2);
            forces_v[ni] += dir * fmag;
        }
    }

    if z_spring != 0.0 {
        for i in 0..n {
            let ni = NodeIndex(i);
            forces_v[ni].z += -z_spring * node_z[ni];
        }
        for i in 0..m {
            let ei = EdgeIndex(i);
            forces_e[ei].z += -z_spring * edge_z[ei];
        }
    }

    (forces_v, forces_e)
}

fn point3_from_point(p: Point2<f64>, z: f64) -> Point3<f64> {
    Point3::new(p.x, p.y, z)
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
