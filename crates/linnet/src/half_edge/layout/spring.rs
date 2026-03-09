use std::ops::IndexMut;

use cgmath::{EuclideanSpace, MetricSpace, Point2, Vector2, Zero};

use rand::{distributions::Uniform, prelude::Distribution, Rng};
use serde::{Deserialize, Serialize};

use crate::{
    half_edge::{
        involution::{EdgeIndex, EdgeVec, HedgePair},
        layout::simulatedanneale::{Energy, Neighbor},
        nodestore::NodeStorageOps,
        subgraph::{subset::SubSet, Inclusion, ModifySubSet, SuBitGraph, SubSetLike},
        swap::Swap,
        HedgeGraph, NodeIndex, NodeVec,
    },
    parser::GlobalData,
};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
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
pub enum ShiftDirection {
    Any,
    PositiveOnly,
    NegativeOnly,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum Constraint {
    Fixed,
    Free,
    Grouped(usize, ShiftDirection),
}

impl Default for Constraint {
    fn default() -> Self {
        Constraint::Free
    }
}

pub trait Shiftable {
    fn shift<I: From<usize> + PartialEq + Copy, R: IndexMut<I, Output = Point2<f64>>>(
        &self,
        shift: Vector2<f64>,
        index: I,
        values: &mut R,
    ) -> bool;
}

pub trait HasPointConstraint {
    fn point_constraint(&self) -> &PointConstraint;
}

impl HasPointConstraint for PointConstraint {
    fn point_constraint(&self) -> &PointConstraint {
        self
    }
}

fn apply_directional_shift(shift_val: f64, direction: ShiftDirection) -> f64 {
    match direction {
        ShiftDirection::Any | ShiftDirection::PositiveOnly | ShiftDirection::NegativeOnly => {
            shift_val
        }
    }
}

pub(crate) fn directional_force_shift(
    constraints: &PointConstraint,
    index: usize,
    magnitude: f64,
) -> Vector2<f64> {
    if magnitude == 0.0 {
        return Vector2::zero();
    }

    let force_x = match constraints.x {
        Constraint::Grouped(reference, ShiftDirection::PositiveOnly) if reference == index => {
            magnitude
        }
        Constraint::Grouped(reference, ShiftDirection::NegativeOnly) if reference == index => {
            -magnitude
        }
        _ => 0.0,
    };

    let force_y = match constraints.y {
        Constraint::Grouped(reference, ShiftDirection::PositiveOnly) if reference == index => {
            magnitude
        }
        Constraint::Grouped(reference, ShiftDirection::NegativeOnly) if reference == index => {
            -magnitude
        }
        _ => 0.0,
    };

    Vector2::from((force_x, force_y))
}

impl Shiftable for PointConstraint {
    fn shift<I: From<usize> + PartialEq + Copy, R: IndexMut<I, Output = Point2<f64>>>(
        &self,
        shift: Vector2<f64>,
        index: I,
        values: &mut R,
    ) -> bool {
        let mut changed = false;

        match (self.x, self.y) {
            (Constraint::Fixed, Constraint::Fixed) => {}
            (Constraint::Free, Constraint::Free) => {
                values[index] += shift;
                changed = true;
            }
            (Constraint::Fixed, Constraint::Free) => {
                values[index].y += shift.y;
                changed = true;
            }
            (Constraint::Free, Constraint::Fixed) => {
                values[index].x += shift.x;
                changed = true;
            }
            (Constraint::Grouped(r, dir), Constraint::Fixed) => {
                let i = r.into();
                if i != index {
                    values[index].x = values[i].x;
                } else {
                    let x_shift = apply_directional_shift(shift.x, dir);
                    values[index].x += x_shift;
                    changed = x_shift != 0.0;
                }
            }
            (Constraint::Grouped(r, dir), Constraint::Free) => {
                let i = r.into();
                if i != index {
                    values[index].x = values[i].x;
                    values[index].y += shift.y;
                    changed = true;
                } else {
                    let x_shift = apply_directional_shift(shift.x, dir);
                    values[index].x += x_shift;
                    values[index].y += shift.y;
                    changed = x_shift != 0.0 || shift.y != 0.0;
                }
            }

            (Constraint::Fixed, Constraint::Grouped(r, dir)) => {
                let i = r.into();
                if i != index {
                    values[index].y = values[i].y;
                } else {
                    let y_shift = apply_directional_shift(shift.y, dir);
                    values[index].y += y_shift;
                    changed = y_shift != 0.0;
                }
            }
            (Constraint::Free, Constraint::Grouped(r, dir)) => {
                let i = r.into();
                if i != index {
                    values[index].y = values[i].y;
                    values[index].x += shift.x;
                    changed = true;
                } else {
                    let y_shift = apply_directional_shift(shift.y, dir);
                    values[index].x += shift.x;
                    values[index].y += y_shift;
                    changed = shift.x != 0.0 || y_shift != 0.0;
                }
            }
            (Constraint::Grouped(xi, x_dir), Constraint::Grouped(yi, y_dir)) => {
                let ix = xi.into();
                let iy = yi.into();
                if ix != index && iy != index {
                    values[index].x = values[ix].x;
                    values[index].y = values[iy].y;
                } else if ix == index && iy != index {
                    let x_shift = apply_directional_shift(shift.x, x_dir);
                    values[index].x += x_shift;
                    values[index].y = values[iy].y;
                    changed = x_shift != 0.0;
                } else if ix != index && iy == index {
                    let y_shift = apply_directional_shift(shift.y, y_dir);
                    values[index].x = values[ix].x;
                    values[index].y += y_shift;
                    changed = y_shift != 0.0;
                } else {
                    let x_shift = apply_directional_shift(shift.x, x_dir);
                    let y_shift = apply_directional_shift(shift.y, y_dir);
                    values[index].x += x_shift;
                    values[index].y += y_shift;
                    changed = x_shift != 0.0 || y_shift != 0.0;
                }
            }
        }
        changed
    }
}

pub struct LayoutState<'a, E, V, H, N: NodeStorageOps<NodeData = V>> {
    pub graph: &'a HedgeGraph<E, V, H, N>,
    pub ext: SuBitGraph,
    pub vertex_points: NodeVec<Point2<f64>>,
    pub edge_points: EdgeVec<Point2<f64>>,
    pub delta: f64,
    pub directional_force: f64,
    // Tracks which node/edge entries were mutated during proposal generation so
    // the energy function can update only the affected terms.
    pub changed_nodes: SubSet<NodeIndex>,
    pub changed_edges: SubSet<EdgeIndex>,
    pub incremental: bool,
}

impl<'a, E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    pub fn new_layout_state(
        &self,
        vertex_points: NodeVec<Point2<f64>>,
        edge_points: EdgeVec<Point2<f64>>,
        delta: f64,
        directional_force: f64,
        incremental: bool,
    ) -> LayoutState<'_, E, V, H, N> {
        let ext = self.external_filter();
        let len_v = vertex_points.len().0;
        let len_e = edge_points.len().0;
        LayoutState {
            graph: self,
            ext,
            vertex_points,
            edge_points,
            delta,
            directional_force,
            changed_nodes: SubSet::empty(len_v),
            changed_edges: SubSet::empty(len_e),
            incremental,
        }
    }
}

impl<'a, E, V, H, N: NodeStorageOps<NodeData = V>> Clone for LayoutState<'a, E, V, H, N> {
    fn clone(&self) -> Self {
        LayoutState {
            graph: self.graph,
            ext: self.ext.clone(),
            vertex_points: self.vertex_points.clone(),
            edge_points: self.edge_points.clone(),
            delta: self.delta,
            directional_force: self.directional_force,
            changed_nodes: self.changed_nodes.clone(),
            changed_edges: self.changed_edges.clone(),
            incremental: self.incremental,
        }
    }
}

impl<'a, E, V, H, N: NodeStorageOps<NodeData = V>> LayoutState<'a, E, V, H, N> {
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

impl<'a, E: Shiftable, V: Shiftable, H, N: NodeStorageOps<NodeData = V> + Clone>
    Neighbor<LayoutState<'a, E, V, H, N>> for LayoutNeighbor
{
    fn propose(
        &self,
        s: &LayoutState<'a, E, V, H, N>,
        rng: &mut impl Rng,
        step: f64,
        _temp: f64,
    ) -> LayoutState<'a, E, V, H, N> {
        let mut st = s.clone();
        let n_v: NodeIndex = st.vertex_points.len();
        let n_e: EdgeIndex = st.edge_points.len();
        let step_range: Uniform<f64> = Uniform::try_from(-step..step).unwrap();

        let mut didnothing = true;
        while didnothing {
            match rng.gen_range(0..100) {
                0..=69 => {
                    // single-DOF
                    if rng.gen_bool(0.6) {
                        let v = NodeIndex(rng.gen_range(0..n_v.0));

                        let shift = LayoutNeighbor::axis_shift(&step_range, rng);
                        let changed = apply_vertex_shift(&mut st, v, shift);
                        didnothing = !changed;
                    } else {
                        let e = EdgeIndex(rng.gen_range(0..n_e.0));
                        let shift = LayoutNeighbor::axis_shift(&step_range, rng);
                        let changed = apply_edge_shift(&mut st, e, shift);
                        didnothing = !changed;
                    }
                }
                _ => {
                    // vertex block
                    let v = NodeIndex(rng.gen_range(0..n_v.0));

                    let shift = LayoutNeighbor::diagonal_shift(&step_range, rng, 0.6);

                    // Cache whether a vertex move succeeded so we can mark it.
                    let mut changed_any = apply_vertex_shift(&mut st, v, shift);

                    st.graph.iter_crown(v).for_each(|a| {
                        let index = st.graph[&a];

                        // Propagate to incident edge control points; any change gets recorded.
                        changed_any |= apply_edge_shift(&mut st, index, shift);
                    });

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
    E: Shiftable + HasPointConstraint,
    V: Shiftable + HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    fn propose(
        &self,
        s: &LayoutState<'a, E, V, H, N>,
        rng: &mut impl Rng,
        step: f64,
        _temp: f64,
    ) -> LayoutState<'a, E, V, H, N> {
        let mut st = s.clone();
        let n_v: NodeIndex = st.vertex_points.len();
        let n_e: EdgeIndex = st.edge_points.len();
        let step_range: Uniform<f64> = Uniform::try_from(-step..step).unwrap();

        let mut didnothing = true;
        while didnothing {
            match rng.gen_range(0..100) {
                0..=6 => {
                    // swap along pinned axis to allow reordering on fixed lines
                    let changed = if rng.gen_bool(0.6) {
                        swap_on_fixed_axis_nodes(&mut st, rng)
                    } else {
                        swap_on_fixed_axis_edges(&mut st, rng)
                    };
                    didnothing = !changed;
                }
                7..=69 => {
                    // single-DOF
                    if rng.gen_bool(0.6) {
                        let v = NodeIndex(rng.gen_range(0..n_v.0));

                        let mut shift = LayoutNeighbor::axis_shift(&step_range, rng);
                        let bias = directional_force_shift(
                            st.graph[v].point_constraint(),
                            v.0,
                            st.directional_force * step,
                        );
                        shift += bias;
                        let changed = apply_vertex_shift_with_groups(&mut st, v, shift);
                        didnothing = !changed;
                    } else {
                        let e = EdgeIndex(rng.gen_range(0..n_e.0));
                        let mut shift = LayoutNeighbor::axis_shift(&step_range, rng);
                        let bias = directional_force_shift(
                            st.graph[e].point_constraint(),
                            e.0,
                            st.directional_force * step,
                        );
                        shift += bias;
                        let changed = apply_edge_shift_with_groups(&mut st, e, shift);
                        didnothing = !changed;
                    }
                }
                _ => {
                    // vertex block
                    let v = NodeIndex(rng.gen_range(0..n_v.0));

                    let shift = LayoutNeighbor::diagonal_shift(&step_range, rng, 0.6);
                    let vertex_bias = directional_force_shift(
                        st.graph[v].point_constraint(),
                        v.0,
                        st.directional_force * step,
                    );

                    // Cache whether a vertex move succeeded so we can mark it.
                    let mut changed_any =
                        apply_vertex_shift_with_groups(&mut st, v, shift + vertex_bias);

                    st.graph.iter_crown(v).for_each(|a| {
                        let index = st.graph[&a];

                        // Propagate to incident edge control points; any change gets recorded.
                        let edge_bias = directional_force_shift(
                            st.graph[index].point_constraint(),
                            index.0,
                            st.directional_force * step,
                        );
                        let edge_shift = shift + edge_bias;
                        changed_any |= apply_edge_shift_with_groups(&mut st, index, edge_shift);
                    });

                    didnothing = !changed_any;
                }
            }
        }
        st
    }
}

pub(crate) fn apply_vertex_shift<
    'a,
    E: Shiftable,
    V: Shiftable,
    H,
    N: NodeStorageOps<NodeData = V> + Clone,
>(
    state: &mut LayoutState<'a, E, V, H, N>,
    idx: NodeIndex,
    shift: Vector2<f64>,
) -> bool {
    let changed = state.graph[idx].shift(shift, idx, &mut state.vertex_points);
    if changed {
        state.mark_node_changed(idx);
    }
    changed
}

pub(crate) fn apply_edge_shift<
    'a,
    E: Shiftable,
    V: Shiftable,
    H,
    N: NodeStorageOps<NodeData = V> + Clone,
>(
    state: &mut LayoutState<'a, E, V, H, N>,
    idx: EdgeIndex,
    shift: Vector2<f64>,
) -> bool {
    let changed = state.graph[idx].shift(shift, idx, &mut state.edge_points);
    if changed {
        state.mark_edge_changed(idx);
    }
    changed
}

fn is_group_reference(constraints: &PointConstraint, reference: usize) -> bool {
    matches!(constraints.x, Constraint::Grouped(r, _) if r == reference)
        || matches!(constraints.y, Constraint::Grouped(r, _) if r == reference)
}

fn propagate_grouped_nodes<'a, E, V, H, N>(
    state: &mut LayoutState<'a, E, V, H, N>,
    reference: NodeIndex,
) -> bool
where
    V: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    let graph = state.graph;
    let reference_point = state.vertex_points[reference];
    let reference_id = reference.0;
    let mut changed_any = false;
    let len = state.vertex_points.len().0;

    for i in 0..len {
        if i == reference_id {
            continue;
        }
        let idx = NodeIndex(i);
        let constraints = graph[idx].point_constraint();
        let mut changed = false;

        if matches!(constraints.x, Constraint::Grouped(r, _) if r == reference_id) {
            if state.vertex_points[idx].x != reference_point.x {
                state.vertex_points[idx].x = reference_point.x;
                changed = true;
            }
        }
        if matches!(constraints.y, Constraint::Grouped(r, _) if r == reference_id) {
            if state.vertex_points[idx].y != reference_point.y {
                state.vertex_points[idx].y = reference_point.y;
                changed = true;
            }
        }

        if changed {
            state.mark_node_changed(idx);
            changed_any = true;
        }
    }

    changed_any
}

fn propagate_grouped_edges<'a, E, V, H, N>(
    state: &mut LayoutState<'a, E, V, H, N>,
    reference: EdgeIndex,
) -> bool
where
    E: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    let graph = state.graph;
    let reference_point = state.edge_points[reference];
    let reference_id = reference.0;
    let mut changed_any = false;
    let len = state.edge_points.len().0;

    for i in 0..len {
        if i == reference_id {
            continue;
        }
        let idx = EdgeIndex(i);
        let constraints = graph[idx].point_constraint();
        let mut changed = false;

        if matches!(constraints.x, Constraint::Grouped(r, _) if r == reference_id) {
            if state.edge_points[idx].x != reference_point.x {
                state.edge_points[idx].x = reference_point.x;
                changed = true;
            }
        }
        if matches!(constraints.y, Constraint::Grouped(r, _) if r == reference_id) {
            if state.edge_points[idx].y != reference_point.y {
                state.edge_points[idx].y = reference_point.y;
                changed = true;
            }
        }

        if changed {
            state.mark_edge_changed(idx);
            changed_any = true;
        }
    }

    changed_any
}

pub(crate) fn apply_vertex_shift_with_groups<
    'a,
    E: Shiftable + HasPointConstraint,
    V: Shiftable + HasPointConstraint,
    H,
    N: NodeStorageOps<NodeData = V> + Clone,
>(
    state: &mut LayoutState<'a, E, V, H, N>,
    idx: NodeIndex,
    shift: Vector2<f64>,
) -> bool {
    let changed = state.graph[idx].shift(shift, idx, &mut state.vertex_points);
    if changed {
        state.mark_node_changed(idx);
        let constraints = state.graph[idx].point_constraint();
        if is_group_reference(constraints, idx.0) {
            return propagate_grouped_nodes(state, idx) || changed;
        }
    }
    changed
}

pub(crate) fn apply_edge_shift_with_groups<
    'a,
    E: Shiftable + HasPointConstraint,
    V: Shiftable + HasPointConstraint,
    H,
    N: NodeStorageOps<NodeData = V> + Clone,
>(
    state: &mut LayoutState<'a, E, V, H, N>,
    idx: EdgeIndex,
    shift: Vector2<f64>,
) -> bool {
    let changed = state.graph[idx].shift(shift, idx, &mut state.edge_points);
    if changed {
        state.mark_edge_changed(idx);
        let constraints = state.graph[idx].point_constraint();
        if is_group_reference(constraints, idx.0) {
            return propagate_grouped_edges(state, idx) || changed;
        }
    }
    changed
}

fn swap_on_fixed_axis_nodes<'a, E, V, H, N>(
    state: &mut LayoutState<'a, E, V, H, N>,
    rng: &mut impl Rng,
) -> bool
where
    V: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    fn is_pinned_axis(constraint: Constraint) -> bool {
        matches!(constraint, Constraint::Fixed | Constraint::Grouped(_, _))
    }

    fn is_free_axis(constraint: Constraint) -> bool {
        matches!(constraint, Constraint::Free)
    }

    let len = state.vertex_points.len().0;
    if len < 2 {
        return false;
    }

    let mut candidates: Vec<NodeIndex> = Vec::new();
    for i in 0..len {
        let idx = NodeIndex(i);
        let constraints = state.graph[idx].point_constraint();
        let free_x_fixed_y = is_free_axis(constraints.x) && is_pinned_axis(constraints.y);
        let fixed_x_free_y = is_pinned_axis(constraints.x) && is_free_axis(constraints.y);
        if free_x_fixed_y || fixed_x_free_y {
            candidates.push(idx);
        }
    }

    if candidates.len() < 2 {
        return false;
    }

    for _ in 0..8 {
        let a = candidates[rng.gen_range(0..candidates.len())];
        let b = candidates[rng.gen_range(0..candidates.len())];
        if a == b {
            continue;
        }
        let ca = state.graph[a].point_constraint();
        let cb = state.graph[b].point_constraint();
        let pa = state.vertex_points[a];
        let pb = state.vertex_points[b];
        let tol = 1e-9;

        let swap_x = is_pinned_axis(ca.y)
            && is_pinned_axis(cb.y)
            && is_free_axis(ca.x)
            && is_free_axis(cb.x)
            && (pa.y - pb.y).abs() <= tol;
        if swap_x {
            state.vertex_points[a].x = pb.x;
            state.vertex_points[b].x = pa.x;
            state.mark_node_changed(a);
            state.mark_node_changed(b);
            return true;
        }

        let swap_y = is_pinned_axis(ca.x)
            && is_pinned_axis(cb.x)
            && is_free_axis(ca.y)
            && is_free_axis(cb.y)
            && (pa.x - pb.x).abs() <= tol;
        if swap_y {
            state.vertex_points[a].y = pb.y;
            state.vertex_points[b].y = pa.y;
            state.mark_node_changed(a);
            state.mark_node_changed(b);
            return true;
        }
    }

    false
}

fn swap_on_fixed_axis_edges<'a, E, V, H, N>(
    state: &mut LayoutState<'a, E, V, H, N>,
    rng: &mut impl Rng,
) -> bool
where
    E: HasPointConstraint,
    N: NodeStorageOps<NodeData = V> + Clone,
{
    fn is_pinned_axis(constraint: Constraint) -> bool {
        matches!(constraint, Constraint::Fixed | Constraint::Grouped(_, _))
    }

    fn is_free_axis(constraint: Constraint) -> bool {
        matches!(constraint, Constraint::Free)
    }

    let len = state.edge_points.len().0;
    if len < 2 {
        return false;
    }

    let mut candidates: Vec<EdgeIndex> = Vec::new();
    for i in 0..len {
        let idx = EdgeIndex(i);
        let constraints = state.graph[idx].point_constraint();
        let free_x_fixed_y = is_free_axis(constraints.x) && is_pinned_axis(constraints.y);
        let fixed_x_free_y = is_pinned_axis(constraints.x) && is_free_axis(constraints.y);
        if free_x_fixed_y || fixed_x_free_y {
            candidates.push(idx);
        }
    }

    if candidates.len() < 2 {
        return false;
    }

    for _ in 0..8 {
        let a = candidates[rng.gen_range(0..candidates.len())];
        let b = candidates[rng.gen_range(0..candidates.len())];
        if a == b {
            continue;
        }
        let ca = state.graph[a].point_constraint();
        let cb = state.graph[b].point_constraint();
        let pa = state.edge_points[a];
        let pb = state.edge_points[b];
        let tol = 1e-9;

        let swap_x = is_pinned_axis(ca.y)
            && is_pinned_axis(cb.y)
            && is_free_axis(ca.x)
            && is_free_axis(cb.x)
            && (pa.y - pb.y).abs() <= tol;
        if swap_x {
            state.edge_points[a].x = pb.x;
            state.edge_points[b].x = pa.x;
            state.mark_edge_changed(a);
            state.mark_edge_changed(b);
            return true;
        }

        let swap_y = is_pinned_axis(ca.x)
            && is_pinned_axis(cb.x)
            && is_free_axis(ca.y)
            && is_free_axis(cb.y)
            && (pa.x - pb.x).abs() <= tol;
        if swap_y {
            state.edge_points[a].y = pb.y;
            state.edge_points[b].y = pa.y;
            state.mark_edge_changed(a);
            state.mark_edge_changed(b);
            return true;
        }
    }

    false
}

#[derive(Clone, Copy)]
pub struct SpringChargeEnergy {
    pub spring_length: f64,    // L
    pub k_spring: f64,         // 1.0
    pub c_vv: f64,             // vertex-vertex charge (≈ 0.14*L^3)
    pub dangling_charge: f64,  // dangling edge charge (≈ 0.14*L^3)
    pub c_ev: f64,             // edge-vertex (≈ 0.028*L^3)
    pub c_ee_local: f64,       // edge-edge local (≈ 0.014*L^3)
    pub c_center: f64,         // central pull (≈ 0.007*L^3)
    pub crossing_penalty: f64, // fixed penalty per crossing
    pub eps: f64,              // 1e-4
}

impl<'a, E, V, H, N: NodeStorageOps<NodeData = V> + Clone> Energy<LayoutState<'a, E, V, H, N>>
    for SpringChargeEnergy
{
    fn energy(
        &self,
        prev: Option<(&LayoutState<'a, E, V, H, N>, f64)>,
        next: &LayoutState<'a, E, V, H, N>,
    ) -> f64 {
        if !next.incremental {
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
    pub length_scale: f64,     // scales L: default 1.0
    pub k_spring: f64,         // spring stiffness: default 1.0
    pub beta: f64,             // vertex–vertex strength
    pub gamma_dangling: f64,   // dangling edge vs vertex–vertex
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
    #[cfg_attr(feature = "energy_trace", inline(never))]
    fn vv_term(&self, dist: f64) -> f64 {
        0.5 * self.c_vv / (dist + self.eps)
    }

    #[cfg_attr(feature = "energy_trace", inline(never))]
    fn ev_term(&self, dist: f64) -> f64 {
        self.c_ev / (dist + self.eps)
    }

    #[cfg_attr(feature = "energy_trace", inline(never))]
    fn spring_term(&self, dist: f64) -> f64 {
        let t = self.spring_length - dist;
        0.5 * self.k_spring * (t * t)
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
    fn center_term(&self, r: f64) -> f64 {
        0.5 * self.c_center / (r + self.eps)
    }

    #[cfg_attr(feature = "energy_trace", inline(never))]
    fn crossing_term(&self) -> f64 {
        self.crossing_penalty
    }

    fn edge_segments<'a, E, V, H, N>(
        s: &LayoutState<'a, E, V, H, N>,
        edge: EdgeIndex,
    ) -> Vec<(Point2<f64>, Point2<f64>, NodeIndex)>
    where
        N: NodeStorageOps<NodeData = V> + Clone,
    {
        let (_, pair) = &s.graph[&edge];
        let p = s.edge_points[edge];
        match *pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                let a = s.graph.node_id(source);
                let b = s.graph.node_id(sink);
                if a == b {
                    Vec::new()
                } else {
                    vec![(s.vertex_points[a], p, a), (s.vertex_points[b], p, b)]
                }
            }
            HedgePair::Unpaired { hedge, .. } => {
                let a = s.graph.node_id(hedge);
                vec![(s.vertex_points[a], p, a)]
            }
        }
    }

    fn edge_spring_length<'a, E, V, H, N>(
        s: &LayoutState<'a, E, V, H, N>,
        edge: EdgeIndex,
        base: f64,
    ) -> f64
    where
        N: NodeStorageOps<NodeData = V> + Clone,
    {
        let (_, pair) = &s.graph[&edge];
        match pair {
            HedgePair::Unpaired { .. } => base * 2.0,
            _ => base,
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
                if ai == aj {
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
            let np = s.vertex_points[ni];
            for j in (i + 1)..n {
                let nj = NodeIndex(j);
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
            let np = s.vertex_points[ni];
            if self.c_ev != 0.0 {
                for e in 0..m {
                    let ei = EdgeIndex(e);
                    let ep = s.edge_points[ei];
                    energy += self.ev_term(np.distance(ep));
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
            let np = s.vertex_points[ni];
            for e in s.graph.iter_crown(ni) {
                let ei = s.graph[&e];
                let ep = s.edge_points[ei];
                let length = Self::edge_spring_length(s, ei, self.spring_length);
                energy += self.spring_term_with_length(np.distance(ep), length);
                for e in s.graph.iter_crown(ni) {
                    let ej = s.graph[&e];
                    if ei == ej {
                        continue;
                    }
                    let ejp = s.edge_points[ej];
                    energy += self.ee_local_term(ejp.distance(ep));
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
            let np = s.vertex_points[ni];
            if self.c_center != 0.0 {
                let r = np.distance(EuclideanSpace::origin());
                if r > 1.0 {
                    energy += self.center_term(r);
                }
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
        #[cfg(feature = "energy_trace")]
        energy_trace::record_dangling(dangling_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let crossing_start = std::time::Instant::now();
        if self.crossing_penalty != 0.0 {
            let m = s.edge_points.len().0;
            for i in 0..m {
                let ei = EdgeIndex(i);
                for j in (i + 1)..m {
                    let ej = EdgeIndex(j);
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

    /// Recompute the energy contributions that touch the mutated nodes/edges.
    fn partial_energy<'a, E, V, H, N>(
        &self,
        s: &LayoutState<'a, E, V, H, N>,
        node_changes: &SubSet<NodeIndex>,
        edge_changes: &SubSet<EdgeIndex>,
    ) -> f64
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
            let node_changed = node_changes.includes(&ni);
            let np = s.vertex_points[ni];

            for j in (i + 1)..n {
                let nj = NodeIndex(j);
                if !(node_changed || node_changes.includes(&nj)) {
                    continue;
                }
                let vj = s.vertex_points[nj];
                // Vertex–vertex repulsion.
                energy += self.vv_term(np.distance(vj));
            }
        }
        #[cfg(feature = "energy_trace")]
        energy_trace::record_vv(vv_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let ev_start = std::time::Instant::now();
        for i in 0..n {
            let ni = NodeIndex(i);
            let node_changed = node_changes.includes(&ni);
            let np = s.vertex_points[ni];
            if self.c_ev != 0.0 {
                for e in 0..m {
                    let ei = EdgeIndex(e);
                    if !(node_changed || edge_changes.includes(&ei)) {
                        continue;
                    }
                    let ep = s.edge_points[ei];
                    // Edge–vertex repulsion.
                    energy += self.ev_term(np.distance(ep));
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
            let node_changed = node_changes.includes(&ni);
            let np = s.vertex_points[ni];

            let include_node = node_changed;
            for hedge in s.graph.iter_crown(ni) {
                let ei = s.graph[&hedge];
                let edge_changed = edge_changes.includes(&ei);
                if !(include_node || edge_changed) {
                    continue;
                }
                let ep = s.edge_points[ei];
                // Spring term for the edge control point.
                let length = Self::edge_spring_length(s, ei, self.spring_length);
                energy += self.spring_term_with_length(np.distance(ep), length);

                for other in s.graph.iter_crown(ni) {
                    let ej = s.graph[&other];
                    if ei == ej {
                        continue;
                    }
                    let other_edge_changed = edge_changes.includes(&ej);
                    if !(include_node || edge_changed || other_edge_changed) {
                        continue;
                    }
                    let ejp = s.edge_points[ej];
                    // Local edge–edge repulsion around `ni`.
                    energy += self.ee_local_term(ejp.distance(ep));
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
            let node_changed = node_changes.includes(&ni);
            if node_changed && self.c_center != 0.0 {
                let np = s.vertex_points[ni];
                let r = np.distance(EuclideanSpace::origin());
                if r > 1.0 {
                    // Central pull acts only on moved vertices.
                    energy += self.center_term(r);
                }
            }
        }
        #[cfg(feature = "energy_trace")]
        energy_trace::record_center(center_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let dangling_start = std::time::Instant::now();
        for hi in s.ext.included_iter() {
            let edge_i = s.graph[&hi];
            let edge_i_changed = edge_changes.includes(&edge_i);
            for hj in s.ext.included_iter() {
                if hi >= hj {
                    continue;
                }
                let edge_j = s.graph[&hj];
                let edge_j_changed = edge_changes.includes(&edge_j);
                if !(edge_i_changed || edge_j_changed) {
                    continue;
                }
                let pi = s.edge_points[edge_i];
                let pj = s.edge_points[edge_j];
                // Dangling charge interactions.
                energy += self.dangling_term(pi.distance(pj));
            }
        }
        #[cfg(feature = "energy_trace")]
        energy_trace::record_dangling(dangling_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let crossing_start = std::time::Instant::now();
        if self.crossing_penalty != 0.0 {
            let m = s.edge_points.len().0;
            for i in 0..m {
                let ei = EdgeIndex(i);
                for j in (i + 1)..m {
                    let ej = EdgeIndex(j);
                    if !(edge_changes.includes(&ei) || edge_changes.includes(&ej)) {
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
        energy_trace::record_partial(total_start.elapsed());
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
            let node_changed = node_changes.includes(&ni);
            let prev_np = prev.vertex_points[ni];
            let next_np = next.vertex_points[ni];
            for j in (i + 1)..n {
                let nj = NodeIndex(j);
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
                let node_changed = node_changes.includes(&ni);
                let prev_np = prev.vertex_points[ni];
                let next_np = next.vertex_points[ni];
                for e in 0..m {
                    let ei = EdgeIndex(e);
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
            let node_changed = node_changes.includes(&ni);
            let prev_np = prev.vertex_points[ni];
            let next_np = next.vertex_points[ni];

            let include_node = node_changed;
            for hedge in next.graph.iter_crown(ni) {
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

                for other in next.graph.iter_crown(ni) {
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
                if !node_changes.includes(&ni) {
                    continue;
                }
                let prev_np = prev.vertex_points[ni];
                let next_np = next.vertex_points[ni];
                let prev_r = prev_np.distance(EuclideanSpace::origin());
                let next_r = next_np.distance(EuclideanSpace::origin());
                if prev_r > 1.0 {
                    delta -= self.center_term(prev_r);
                }
                if next_r > 1.0 {
                    delta += self.center_term(next_r);
                }
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
        #[cfg(feature = "energy_trace")]
        energy_trace::record_dangling(dangling_start.elapsed());

        #[cfg(feature = "energy_trace")]
        let crossing_start = std::time::Instant::now();
        if self.crossing_penalty != 0.0 {
            let m = next.edge_points.len().0;
            for i in 0..m {
                let ei = EdgeIndex(i);
                for j in (i + 1)..m {
                    let ej = EdgeIndex(j);
                    if !(edge_changes.includes(&ei) || edge_changes.includes(&ej)) {
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
        let area = (viewport_w * viewport_h).max(1e-9);
        let spring_length = tune.length_scale * (area / (n_nodes.max(1) as f64)).sqrt();

        SpringChargeEnergy {
            spring_length,
            k_spring: tune.k_spring,
            c_vv: tune.beta * spring_length.powi(2),
            c_ev: tune.beta * tune.gamma_ev * spring_length.powi(2),
            c_ee_local: tune.beta * tune.gamma_ee * spring_length.powi(2),
            c_center: tune.beta * tune.g_center * spring_length.powi(2),
            dangling_charge: tune.gamma_dangling * tune.beta * spring_length.powi(2),
            crossing_penalty: tune.crossing_penalty,
            eps: tune.eps,
        }
    }
}
