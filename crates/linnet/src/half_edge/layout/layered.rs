use std::collections::BTreeMap;

use cgmath::Point2;

use crate::half_edge::{
    involution::{EdgeIndex, EdgeVec, HedgePair},
    nodestore::NodeStorageOps,
    subgraph::SubSetLike,
    HedgeGraph, NodeIndex, NodeVec,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LayeredProfile {
    Dot,
    Stable,
}

#[derive(Debug, Clone)]
pub struct LayeredConfig {
    pub profile: LayeredProfile,
    pub layer_gap: f64,
    pub node_gap: f64,
    pub edge_gap: f64,
    pub route_edge_weight: f64,
    pub route_label_width_scale: f64,
    pub route_label_width_cap: f64,
    pub sweeps: usize,
    pub roots: Vec<NodeIndex>,
    pub rank_same: Vec<Vec<NodeIndex>>,
    pub include_all_nodes: bool,
}

impl Default for LayeredConfig {
    fn default() -> Self {
        Self {
            profile: LayeredProfile::Dot,
            layer_gap: 1.2,
            node_gap: 0.9,
            edge_gap: 0.35,
            route_edge_weight: 0.15,
            route_label_width_scale: 1.0,
            route_label_width_cap: 2.0,
            sweeps: 8,
            roots: Vec::new(),
            rank_same: Vec::new(),
            include_all_nodes: true,
        }
    }
}

#[derive(Debug, Clone)]
pub struct LayeredGeometry {
    pub node_widths: NodeVec<f64>,
    pub node_heights: NodeVec<f64>,
    pub edge_label_widths: EdgeVec<f64>,
    pub edge_label_heights: EdgeVec<f64>,
    pub edge_minlens: EdgeVec<usize>,
    pub edge_weights: EdgeVec<f64>,
    pub edge_constrained: EdgeVec<bool>,
}

#[derive(Debug, Clone)]
pub struct LayeredOutput {
    pub node_positions: NodeVec<Option<Point2<f64>>>,
    pub edge_positions: EdgeVec<Option<Point2<f64>>>,
    pub edge_routes: EdgeVec<LayeredEdgeRoute>,
    pub ranks: NodeVec<Option<usize>>,
}

#[derive(Debug, Clone, Default)]
pub struct LayeredEdgeRoute {
    pub source: Vec<Point2<f64>>,
    pub sink: Vec<Point2<f64>>,
}

#[derive(Debug, Clone)]
struct LayoutEdge {
    edge: EdgeIndex,
    source: NodeIndex,
    sink: NodeIndex,
    rank_edge: bool,
}

#[derive(Debug, Clone)]
enum ItemKind {
    Real(NodeIndex),
    Dummy,
}

#[derive(Debug, Clone)]
struct Item {
    kind: ItemKind,
    rank: usize,
    width: f64,
    height: f64,
    order_key: f64,
    x: f64,
}

#[derive(Debug, Clone)]
struct SameRankEdge {
    edge: EdgeIndex,
    start_x: f64,
    end_x: f64,
    sink: NodeIndex,
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    pub fn layered_layout<S: SubSetLike>(
        &self,
        subgraph: &S,
        config: &LayeredConfig,
        geometry: &LayeredGeometry,
    ) -> LayeredOutput {
        let included = self.layered_included_nodes(subgraph, config.include_all_nodes);
        let ordered_nodes = self.layered_ordered_nodes(&included, &config.roots);
        let order = self.layered_order_positions(&ordered_nodes);
        let selected_edges = self.layered_selected_edges(subgraph, &included);
        let route_edges = self.layered_route_edges(&included, &selected_edges);
        let same_rank_group = self.layered_same_rank_groups(&included, &config.rank_same);
        let ranks = self.layered_ranks(
            &included,
            &ordered_nodes,
            &order,
            &selected_edges,
            &same_rank_group,
            geometry,
        );
        let layout_order = if config.profile == LayeredProfile::Dot {
            self.layered_hierarchy_order(&included, &ordered_nodes, &order, &selected_edges, &ranks)
        } else {
            order.iter().map(|&position| position as f64).collect()
        };

        let included_count = included.iter().filter(|(_, included)| **included).count();
        if included_count == 0 {
            return LayeredOutput {
                node_positions: self.new_nodevec(|_, _, _| None),
                edge_positions: self.new_edgevec(|_, _, _| None),
                edge_routes: (0..self.n_edges())
                    .map(|_| LayeredEdgeRoute::default())
                    .collect(),
                ranks: self.new_nodevec(|_, _, _| None),
            };
        }

        let mut layout = LayeredWorkspace::new(
            self,
            config,
            geometry,
            &included,
            &ordered_nodes,
            &layout_order,
            &route_edges,
            &ranks,
        );
        layout.order_items();
        layout.assign_initial_x();
        layout.straighten();
        let lane_offsets = layout.same_rank_lane_offsets();
        let rank_track_heights = layout.rank_track_heights(&lane_offsets);
        let rank_y = layout.rank_y_positions(&rank_track_heights);
        layout.into_output(&ranks, &route_edges, &lane_offsets, &rank_y)
    }

    fn layered_included_nodes<S: SubSetLike>(
        &self,
        subgraph: &S,
        include_all_nodes: bool,
    ) -> NodeVec<bool> {
        let mut included = self.new_nodevec(|_, _, _| include_all_nodes);
        if !include_all_nodes {
            for hedge in subgraph.included_iter() {
                included[self.node_id(hedge)] = true;
            }
        }
        included
    }

    fn layered_ordered_nodes(
        &self,
        included: &NodeVec<bool>,
        roots: &[NodeIndex],
    ) -> Vec<NodeIndex> {
        let mut seen = vec![false; self.n_nodes()];
        let mut nodes = Vec::new();

        for &root in roots {
            if root.0 < seen.len() && included[root] && !seen[root.0] {
                seen[root.0] = true;
                nodes.push(root);
            }
        }

        for (node, &is_included) in included.iter() {
            if is_included && !seen[node.0] {
                seen[node.0] = true;
                nodes.push(node);
            }
        }

        nodes
    }

    fn layered_order_positions(&self, ordered_nodes: &[NodeIndex]) -> Vec<usize> {
        let mut order = vec![usize::MAX; self.n_nodes()];
        for (position, &node) in ordered_nodes.iter().enumerate() {
            order[node.0] = position;
        }
        order
    }

    fn layered_selected_edges<S: SubSetLike>(
        &self,
        subgraph: &S,
        included: &NodeVec<bool>,
    ) -> Vec<LayoutEdge> {
        self.iter_edges_of(subgraph)
            .filter_map(|(pair, edge, _)| match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    let source = self.node_id(source);
                    let sink = self.node_id(sink);
                    (source != sink && included[source] && included[sink]).then_some(LayoutEdge {
                        edge,
                        source,
                        sink,
                        rank_edge: true,
                    })
                }
                HedgePair::Unpaired { .. } => None,
            })
            .collect()
    }

    fn layered_route_edges(
        &self,
        included: &NodeVec<bool>,
        selected_edges: &[LayoutEdge],
    ) -> Vec<LayoutEdge> {
        let mut rank_edge = (0..self.n_edges()).map(|_| false).collect::<EdgeVec<_>>();
        for edge in selected_edges {
            rank_edge[edge.edge] = true;
        }

        self.iter_edges()
            .filter_map(|(pair, edge, _)| match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    let source = self.node_id(source);
                    let sink = self.node_id(sink);
                    (source != sink && included[source] && included[sink]).then_some(LayoutEdge {
                        edge,
                        source,
                        sink,
                        rank_edge: rank_edge[edge],
                    })
                }
                HedgePair::Unpaired { .. } => None,
            })
            .collect()
    }

    fn layered_same_rank_groups(
        &self,
        included: &NodeVec<bool>,
        groups: &[Vec<NodeIndex>],
    ) -> NodeVec<Option<usize>> {
        let mut group_by_node = self.new_nodevec(|_, _, _| None);
        for (group_index, group) in groups.iter().enumerate() {
            for &node in group {
                if node.0 < self.n_nodes() && included[node] {
                    group_by_node[node] = Some(group_index);
                }
            }
        }
        group_by_node
    }

    fn layered_ranks(
        &self,
        included: &NodeVec<bool>,
        ordered_nodes: &[NodeIndex],
        order: &[usize],
        selected_edges: &[LayoutEdge],
        same_rank_group: &NodeVec<Option<usize>>,
        geometry: &LayeredGeometry,
    ) -> NodeVec<usize> {
        let mut ranks = self.new_nodevec(|_, _, _| 0usize);
        let mut dag_edges = selected_edges
            .iter()
            .filter_map(|edge| {
                if !geometry.edge_constrained[edge.edge] {
                    return None;
                }
                if same_rank_group[edge.source].is_some()
                    && same_rank_group[edge.source] == same_rank_group[edge.sink]
                {
                    return None;
                }
                let minlen = geometry.edge_minlens[edge.edge].max(1);
                if order[edge.source.0] <= order[edge.sink.0] {
                    Some((edge.source, edge.sink, minlen))
                } else {
                    Some((edge.sink, edge.source, minlen))
                }
            })
            .collect::<Vec<_>>();

        dag_edges.sort_by_key(|(source, sink, _)| (order[source.0], order[sink.0]));

        let max_iterations = ordered_nodes
            .len()
            .saturating_add(selected_edges.len())
            .max(1);
        for _ in 0..max_iterations {
            let mut changed = false;
            for &(source, sink, minlen) in &dag_edges {
                let target = ranks[source].saturating_add(minlen);
                if ranks[sink] < target {
                    ranks[sink] = target;
                    changed = true;
                }
            }

            for (group_index, group) in self.valid_rank_same_groups(included, same_rank_group) {
                let group_rank = group.iter().map(|&node| ranks[node]).max().unwrap_or(0);
                for node in group {
                    if same_rank_group[node] == Some(group_index) && ranks[node] != group_rank {
                        ranks[node] = group_rank;
                        changed = true;
                    }
                }
            }

            if !changed {
                break;
            }
        }

        ranks
    }

    fn valid_rank_same_groups(
        &self,
        included: &NodeVec<bool>,
        same_rank_group: &NodeVec<Option<usize>>,
    ) -> BTreeMap<usize, Vec<NodeIndex>> {
        let mut groups = BTreeMap::<usize, Vec<NodeIndex>>::new();
        for (node, &is_included) in included.iter() {
            if is_included {
                if let Some(group) = same_rank_group[node] {
                    groups.entry(group).or_default().push(node);
                }
            }
        }
        groups
    }

    fn layered_hierarchy_order(
        &self,
        included: &NodeVec<bool>,
        ordered_nodes: &[NodeIndex],
        order: &[usize],
        selected_edges: &[LayoutEdge],
        ranks: &NodeVec<usize>,
    ) -> Vec<f64> {
        let mut children = vec![Vec::<NodeIndex>::new(); self.n_nodes()];
        let mut has_parent = vec![false; self.n_nodes()];
        let mut forest_parent = (0..self.n_nodes()).collect::<Vec<_>>();
        let mut candidates = Vec::<(NodeIndex, NodeIndex)>::new();
        for edge in selected_edges {
            if ranks[edge.source] == ranks[edge.sink] {
                continue;
            }
            let (parent, child) = if ranks[edge.source] < ranks[edge.sink] {
                (edge.source, edge.sink)
            } else {
                (edge.sink, edge.source)
            };
            if included[parent] && included[child] {
                candidates.push((parent, child));
            }
        }

        candidates.sort_by_key(|(parent, child)| {
            (
                ranks[*parent],
                order[parent.0],
                ranks[*child],
                order[child.0],
            )
        });
        for (parent, child) in candidates {
            let parent_root = find_disjoint_root(&mut forest_parent, parent.0);
            let child_root = find_disjoint_root(&mut forest_parent, child.0);
            if parent_root == child_root {
                continue;
            }
            forest_parent[child_root] = parent_root;
            children[parent.0].push(child);
            has_parent[child.0] = true;
        }

        for node_children in &mut children {
            node_children.sort_by_key(|node| (ranks[*node], order[node.0]));
            node_children.dedup();
        }

        let mut seen = vec![false; self.n_nodes()];
        let mut hierarchy_order = vec![f64::INFINITY; self.n_nodes()];
        let mut next = 0usize;
        for &root in ordered_nodes {
            if !has_parent[root.0] {
                Self::assign_layered_hierarchy_order(
                    root,
                    &children,
                    &mut seen,
                    &mut hierarchy_order,
                    &mut next,
                );
            }
        }
        for &root in ordered_nodes {
            Self::assign_layered_hierarchy_order(
                root,
                &children,
                &mut seen,
                &mut hierarchy_order,
                &mut next,
            );
        }

        for (node, &is_included) in included.iter() {
            if is_included && !seen[node.0] {
                hierarchy_order[node.0] = next as f64;
                next += 1;
            }
        }

        hierarchy_order
    }

    fn assign_layered_hierarchy_order(
        node: NodeIndex,
        children: &[Vec<NodeIndex>],
        seen: &mut [bool],
        hierarchy_order: &mut [f64],
        next: &mut usize,
    ) {
        if seen[node.0] {
            return;
        }
        seen[node.0] = true;
        hierarchy_order[node.0] = *next as f64;
        *next += 1;
        for &child in &children[node.0] {
            Self::assign_layered_hierarchy_order(child, children, seen, hierarchy_order, next);
        }
    }
}

struct LayeredWorkspace<'a, E, V, H, N: NodeStorageOps<NodeData = V>> {
    graph: &'a HedgeGraph<E, V, H, N>,
    config: &'a LayeredConfig,
    geometry: &'a LayeredGeometry,
    layers: Vec<Vec<usize>>,
    items: Vec<Item>,
    node_items: Vec<Option<usize>>,
    edge_paths: EdgeVec<Vec<usize>>,
    edge_weights: EdgeVec<f64>,
    neighbors: Vec<Vec<(usize, f64)>>,
}

impl<'a, E, V, H, N: NodeStorageOps<NodeData = V>> LayeredWorkspace<'a, E, V, H, N> {
    #[allow(clippy::too_many_arguments)]
    fn new(
        graph: &'a HedgeGraph<E, V, H, N>,
        config: &'a LayeredConfig,
        geometry: &'a LayeredGeometry,
        included: &NodeVec<bool>,
        ordered_nodes: &[NodeIndex],
        order: &[f64],
        selected_edges: &[LayoutEdge],
        ranks: &NodeVec<usize>,
    ) -> Self {
        let max_rank = included
            .iter()
            .filter_map(|(node, included)| included.then_some(ranks[node]))
            .max()
            .unwrap_or(0);
        let mut layout = Self {
            graph,
            config,
            geometry,
            layers: vec![Vec::new(); max_rank + 1],
            items: Vec::new(),
            node_items: vec![None; graph.n_nodes()],
            edge_paths: (0..graph.n_edges()).map(|_| Vec::new()).collect(),
            edge_weights: (0..graph.n_edges()).map(|_| 0.0).collect(),
            neighbors: Vec::new(),
        };

        for &node in ordered_nodes {
            if !included[node] {
                continue;
            }
            let item = layout.push_item(Item {
                kind: ItemKind::Real(node),
                rank: ranks[node],
                width: geometry.node_widths[node],
                height: geometry.node_heights[node],
                order_key: order[node.0],
                x: 0.0,
            });
            layout.node_items[node.0] = Some(item);
        }

        for edge in selected_edges {
            let Some(source_item) = layout.node_items[edge.source.0] else {
                continue;
            };
            let Some(sink_item) = layout.node_items[edge.sink.0] else {
                continue;
            };

            let source_rank = ranks[edge.source];
            let sink_rank = ranks[edge.sink];
            let mut path = vec![source_item];
            let rank_delta = sink_rank as isize - source_rank as isize;
            let steps = rank_delta.unsigned_abs();
            if steps > 1 {
                let direction = rank_delta.signum();
                let middle_step = steps / 2;
                for step in 1..steps {
                    let rank = (source_rank as isize + direction * step as isize) as usize;
                    let is_middle = step == middle_step;
                    let label_width = geometry.edge_label_widths[edge.edge];
                    let label_height = geometry.edge_label_heights[edge.edge];
                    let width = if is_middle {
                        if edge.rank_edge {
                            label_width.max(config.edge_gap)
                        } else {
                            let width =
                                (label_width * config.route_label_width_scale).max(config.edge_gap);
                            if config.route_label_width_cap > 0.0 {
                                width.min(
                                    (config.node_gap * config.route_label_width_cap)
                                        .max(config.edge_gap),
                                )
                            } else {
                                width
                            }
                        }
                    } else {
                        config.edge_gap * 0.5
                    };
                    let height = if is_middle { label_height } else { 0.0 };
                    let order_key =
                        0.5 * (order[edge.source.0] + order[edge.sink.0]) + 0.001 * step as f64;
                    let dummy = layout.push_item(Item {
                        kind: ItemKind::Dummy,
                        rank,
                        width,
                        height,
                        order_key,
                        x: 0.0,
                    });
                    path.push(dummy);
                }
            }
            path.push(sink_item);
            layout.edge_paths[edge.edge] = path;
            layout.edge_weights[edge.edge] = if edge.rank_edge {
                geometry.edge_weights[edge.edge].max(0.0)
            } else {
                geometry.edge_weights[edge.edge].max(0.0) * config.route_edge_weight.max(0.0)
            };
        }

        layout.neighbors = vec![Vec::new(); layout.items.len()];
        for (edge, path) in layout.edge_paths.iter() {
            let weight = layout.edge_weights[edge];
            for segment in path.windows(2) {
                let a = segment[0];
                let b = segment[1];
                layout.neighbors[a].push((b, weight));
                layout.neighbors[b].push((a, weight));
            }
        }

        layout
    }

    fn push_item(&mut self, item: Item) -> usize {
        let id = self.items.len();
        let rank = item.rank;
        self.items.push(item);
        self.layers[rank].push(id);
        id
    }

    fn order_items(&mut self) {
        for layer in &mut self.layers {
            layer.sort_by(|&a, &b| {
                self.items[a]
                    .order_key
                    .partial_cmp(&self.items[b].order_key)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.cmp(&b))
            });
        }

        if self.config.profile == LayeredProfile::Dot {
            for _ in 0..self.config.sweeps {
                for rank in 1..self.layers.len() {
                    self.sort_layer_by_barycenter(rank, Direction::Up);
                }
                for rank in (0..self.layers.len().saturating_sub(1)).rev() {
                    self.sort_layer_by_barycenter(rank, Direction::Down);
                }
                self.transpose_to_reduce_crossings();
            }
        }
    }

    fn sort_layer_by_barycenter(&mut self, rank: usize, direction: Direction) {
        let positions = self.item_positions();
        let mut keyed = self.layers[rank]
            .iter()
            .enumerate()
            .map(|(old_position, &item)| {
                let mut weighted_sum = 0.0;
                let mut total_weight = 0.0;
                for &(neighbor, weight) in &self.neighbors[item] {
                    let use_neighbor = match direction {
                        Direction::Up => self.items[neighbor].rank < rank,
                        Direction::Down => self.items[neighbor].rank > rank,
                    };
                    if use_neighbor {
                        weighted_sum += positions[neighbor] as f64 * weight.max(0.05);
                        total_weight += weight.max(0.05);
                    }
                }
                let barycenter = if total_weight > 0.0 {
                    weighted_sum / total_weight
                } else {
                    old_position as f64
                };
                (barycenter, old_position, item)
            })
            .collect::<Vec<_>>();

        keyed.sort_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| a.1.cmp(&b.1))
        });
        self.layers[rank] = keyed
            .into_iter()
            .map(|(_, _, item)| item)
            .collect::<Vec<_>>();
    }

    fn item_positions(&self) -> Vec<usize> {
        let mut positions = vec![0usize; self.items.len()];
        for layer in &self.layers {
            for (position, &item) in layer.iter().enumerate() {
                positions[item] = position;
            }
        }
        positions
    }

    fn transpose_to_reduce_crossings(&mut self) {
        for _ in 0..self.config.sweeps.max(1) {
            let mut changed = false;
            for rank in 0..self.layers.len() {
                changed |= self.transpose_layer_to_reduce_crossings(rank);
            }
            for rank in (0..self.layers.len()).rev() {
                changed |= self.transpose_layer_to_reduce_crossings(rank);
            }
            if !changed {
                break;
            }
        }
    }

    fn transpose_layer_to_reduce_crossings(&mut self, rank: usize) -> bool {
        let mut changed = false;
        if self.layers[rank].len() < 2 {
            return changed;
        }

        for position in 0..self.layers[rank].len() - 1 {
            let before = self.crossing_score_around_rank(rank) + self.visible_edge_crossing_score();
            self.layers[rank].swap(position, position + 1);
            let after = self.crossing_score_around_rank(rank) + self.visible_edge_crossing_score();
            if after + 1e-9 < before {
                changed = true;
            } else {
                self.layers[rank].swap(position, position + 1);
            }
        }
        changed
    }

    fn crossing_score_around_rank(&self, rank: usize) -> f64 {
        let mut score = 0.0;
        if rank > 0 {
            score += self.crossing_score_between_ranks(rank - 1, rank);
        }
        if rank + 1 < self.layers.len() {
            score += self.crossing_score_between_ranks(rank, rank + 1);
        }
        score
    }

    fn crossing_score_between_ranks(&self, upper_rank: usize, lower_rank: usize) -> f64 {
        let positions = self.item_positions();
        let mut segments = Vec::<(usize, usize, f64)>::new();
        for &upper in &self.layers[upper_rank] {
            let upper_position = positions[upper];
            for &(neighbor, weight) in &self.neighbors[upper] {
                if self.items[neighbor].rank == lower_rank {
                    segments.push((upper_position, positions[neighbor], weight.max(0.05)));
                }
            }
        }

        let mut score = 0.0;
        for (left_index, &(left_upper, left_lower, left_weight)) in segments.iter().enumerate() {
            for &(right_upper, right_lower, right_weight) in &segments[left_index + 1..] {
                if left_upper == right_upper || left_lower == right_lower {
                    continue;
                }
                if (left_upper < right_upper && left_lower > right_lower)
                    || (left_upper > right_upper && left_lower < right_lower)
                {
                    score += left_weight * right_weight;
                }
            }
        }
        score
    }

    fn visible_edge_crossing_score(&self) -> f64 {
        let positions = self.item_positions();
        let mut segments = Vec::<VisibleSegment>::new();
        for (edge, path) in self.edge_paths.iter() {
            if path.len() < 2 {
                continue;
            }
            let source = path[0];
            let sink = *path.last().unwrap();
            if self.items[source].rank == self.items[sink].rank {
                continue;
            }
            segments.push(VisibleSegment {
                source,
                sink,
                source_rank: self.items[source].rank as f64,
                source_order: positions[source] as f64,
                sink_rank: self.items[sink].rank as f64,
                sink_order: positions[sink] as f64,
                weight: self.edge_weights[edge].max(0.05),
            });
        }

        let mut score = 0.0;
        for (left_index, left) in segments.iter().enumerate() {
            for right in &segments[left_index + 1..] {
                if left.source == right.source
                    || left.source == right.sink
                    || left.sink == right.source
                    || left.sink == right.sink
                {
                    continue;
                }
                if segments_cross(left, right) {
                    score += left.weight * right.weight;
                }
            }
        }
        score
    }

    fn assign_initial_x(&mut self) {
        for rank in 0..self.layers.len() {
            let mut x = 0.0;
            for (position, &item) in self.layers[rank].iter().enumerate() {
                if position > 0 {
                    let previous = self.layers[rank][position - 1];
                    x += 0.5 * self.items[previous].width
                        + self.config.node_gap
                        + 0.5 * self.items[item].width;
                }
                self.items[item].x = x;
            }
            self.center_layer(rank);
        }
    }

    fn straighten(&mut self) {
        self.solve_horizontal_constraints();

        if self.config.profile == LayeredProfile::Stable {
            self.offset_stable_dummy_runs();
        }
        self.center_all_layers();
    }

    fn solve_horizontal_constraints(&mut self) {
        let passes = match self.config.profile {
            LayeredProfile::Dot => self.config.sweeps.max(1),
            LayeredProfile::Stable => 4,
        };

        for _ in 0..passes {
            let mut targets = self
                .items
                .iter()
                .map(|item| (item.x, 1.0))
                .collect::<Vec<_>>();
            for (edge, path) in self.edge_paths.iter() {
                let weight = self.edge_weights[edge].max(0.05);
                for segment in path.windows(2) {
                    let a = segment[0];
                    let b = segment[1];
                    let ax = self.items[a].x;
                    let bx = self.items[b].x;
                    targets[a].0 += bx * weight;
                    targets[a].1 += weight;
                    targets[b].0 += ax * weight;
                    targets[b].1 += weight;
                }
            }

            for (item, (sum, weight)) in self.items.iter_mut().zip(targets) {
                let target = sum / weight;
                item.x = 0.5 * item.x + 0.5 * target;
            }

            for rank in 0..self.layers.len() {
                self.project_layer_constraints(rank);
            }
        }
    }

    fn project_layer_constraints(&mut self, rank: usize) {
        if self.layers[rank].is_empty() {
            return;
        }
        self.resolve_layer_overlaps(rank);
        self.center_layer(rank);
    }

    fn offset_stable_dummy_runs(&mut self) {
        for (_, path) in self.edge_paths.iter() {
            if path.len() <= 2 {
                continue;
            }
            let source = path[0];
            let sink = *path.last().unwrap();
            let source_x = self.items[source].x;
            let sink_x = self.items[sink].x;
            let offset = if self.items[source].rank <= self.items[sink].rank {
                -self.config.edge_gap
            } else {
                self.config.edge_gap
            };
            let target = source_x.min(sink_x) + 0.5 * (source_x - sink_x).abs() + offset;
            for &item in &path[1..path.len() - 1] {
                self.items[item].x = target;
            }
        }
        for rank in 0..self.layers.len() {
            self.resolve_layer_overlaps(rank);
        }
    }

    fn resolve_layer_overlaps(&mut self, rank: usize) {
        let Some((&first, rest)) = self.layers[rank].split_first() else {
            return;
        };
        let mut right = self.items[first].x + 0.5 * self.items[first].width;
        for &item in rest {
            let min_x = right + self.config.node_gap + 0.5 * self.items[item].width;
            if self.items[item].x < min_x {
                self.items[item].x = min_x;
            }
            right = self.items[item].x + 0.5 * self.items[item].width;
        }
    }

    fn center_layer(&mut self, rank: usize) {
        let Some((&first, rest)) = self.layers[rank].split_first() else {
            return;
        };
        let last = rest.last().copied().unwrap_or(first);
        let min_x = self.items[first].x - 0.5 * self.items[first].width;
        let max_x = self.items[last].x + 0.5 * self.items[last].width;
        let center = 0.5 * (min_x + max_x);
        for &item in &self.layers[rank] {
            self.items[item].x -= center;
        }
    }

    fn center_all_layers(&mut self) {
        let mut min_x = f64::INFINITY;
        let mut max_x = f64::NEG_INFINITY;
        for item in &self.items {
            min_x = min_x.min(item.x - 0.5 * item.width);
            max_x = max_x.max(item.x + 0.5 * item.width);
        }
        if !min_x.is_finite() || !max_x.is_finite() {
            return;
        }
        let center = 0.5 * (min_x + max_x);
        for item in &mut self.items {
            item.x -= center;
        }
    }

    fn same_rank_lane_offsets(&self) -> EdgeVec<Option<f64>> {
        let mut offsets: EdgeVec<Option<f64>> = (0..self.graph.n_edges()).map(|_| None).collect();
        let mut by_rank = BTreeMap::<usize, Vec<SameRankEdge>>::new();

        for (edge, path) in self.edge_paths.iter() {
            if path.len() != 2 {
                continue;
            }
            let source = path[0];
            let sink = path[1];
            if self.items[source].rank != self.items[sink].rank {
                continue;
            }
            let ItemKind::Real(sink_node) = self.items[sink].kind else {
                continue;
            };
            by_rank
                .entry(self.items[source].rank)
                .or_default()
                .push(SameRankEdge {
                    edge,
                    start_x: self.items[source].x,
                    end_x: self.items[sink].x,
                    sink: sink_node,
                });
        }

        for edges in by_rank.values_mut() {
            edges.sort_by(|a, b| {
                a.start_x
                    .partial_cmp(&b.start_x)
                    .unwrap_or(std::cmp::Ordering::Equal)
                    .then_with(|| a.edge.0.cmp(&b.edge.0))
            });
            let mut rightward = Vec::<Vec<SameRankEdge>>::new();
            let mut leftward = Vec::<Vec<SameRankEdge>>::new();
            for edge in edges.iter().cloned() {
                let tracks = if edge.end_x >= edge.start_x {
                    &mut rightward
                } else {
                    &mut leftward
                };
                Self::assign_same_rank_track(edge, tracks);
            }

            for (track_index, track) in rightward.iter().enumerate() {
                let offset = (track_index + 1) as f64 * self.config.edge_gap;
                for edge in track {
                    offsets[edge.edge] = Some(offset);
                }
            }
            for (track_index, track) in leftward.iter().enumerate() {
                let offset = -((track_index + 1) as f64) * self.config.edge_gap;
                for edge in track {
                    offsets[edge.edge] = Some(offset);
                }
            }
        }

        offsets
    }

    fn assign_same_rank_track(edge: SameRankEdge, tracks: &mut Vec<Vec<SameRankEdge>>) {
        let mut last_valid = None;
        for (track_index, track) in tracks.iter().enumerate().rev() {
            let mut overlaps = false;
            for other in track {
                if edge.sink == other.sink {
                    tracks[track_index].push(edge);
                    return;
                }
                if intervals_overlap(edge.start_x, edge.end_x, other.start_x, other.end_x) {
                    overlaps = true;
                    break;
                }
            }
            if overlaps {
                break;
            }
            last_valid = Some(track_index);
        }

        if let Some(track_index) = last_valid {
            tracks[track_index].push(edge);
        } else {
            tracks.push(vec![edge]);
        }
    }

    fn rank_track_heights(&self, lane_offsets: &EdgeVec<Option<f64>>) -> Vec<f64> {
        let mut heights = vec![0.0_f64; self.layers.len()];
        for (edge, offset) in lane_offsets.iter() {
            let Some(offset) = offset else {
                continue;
            };
            let label_height = self.geometry.edge_label_heights[edge];
            let height = offset.abs() + label_height + self.config.edge_gap;
            if let Some(rank) = self.edge_paths[edge]
                .first()
                .map(|&item| self.items[item].rank)
            {
                heights[rank] = heights[rank].max(height);
            }
        }
        heights
    }

    fn rank_y_positions(&self, rank_track_heights: &[f64]) -> Vec<f64> {
        let rank_heights = self
            .layers
            .iter()
            .map(|layer| {
                layer
                    .iter()
                    .map(|&item| self.items[item].height)
                    .fold(0.0, f64::max)
            })
            .collect::<Vec<_>>();
        let mut y = vec![0.0; self.layers.len()];
        for rank in 1..self.layers.len() {
            y[rank] = y[rank - 1]
                - (0.5 * rank_heights[rank - 1]
                    + self.config.layer_gap
                    + rank_track_heights[rank - 1]
                    + 0.5 * rank_heights[rank]);
        }
        y
    }

    fn into_output(
        self,
        ranks: &NodeVec<usize>,
        selected_edges: &[LayoutEdge],
        lane_offsets: &EdgeVec<Option<f64>>,
        rank_y: &[f64],
    ) -> LayeredOutput {
        let mut node_positions = self.graph.new_nodevec(|_, _, _| None);
        let mut rank_output = self.graph.new_nodevec(|_, _, _| None);
        for (node_index, item) in self.node_items.iter().enumerate() {
            let Some(item) = item else {
                continue;
            };
            let node = NodeIndex(node_index);
            node_positions[node] = Some(Point2::new(self.items[*item].x, rank_y[ranks[node]]));
            rank_output[node] = Some(ranks[node]);
        }

        let mut edge_positions: EdgeVec<Option<Point2<f64>>> =
            (0..self.graph.n_edges()).map(|_| None).collect();
        let mut edge_routes: EdgeVec<LayeredEdgeRoute> = (0..self.graph.n_edges())
            .map(|_| LayeredEdgeRoute::default())
            .collect();
        for edge in selected_edges {
            let path = &self.edge_paths[edge.edge];
            if path.len() < 2 {
                continue;
            }
            let source = path[0];
            let sink = *path.last().unwrap();
            let source_rank = self.items[source].rank;
            let sink_rank = self.items[sink].rank;
            let point = if source_rank == sink_rank {
                let y =
                    rank_y[source_rank] + lane_offsets[edge.edge].unwrap_or(self.config.edge_gap);
                Point2::new(0.5 * (self.items[source].x + self.items[sink].x), y)
            } else if path.len() > 2 {
                let middle = path[path.len() / 2];
                Point2::new(self.items[middle].x, rank_y[self.items[middle].rank])
            } else {
                Point2::new(
                    0.5 * (self.items[source].x + self.items[sink].x),
                    0.5 * (rank_y[source_rank] + rank_y[sink_rank]),
                )
            };
            edge_positions[edge.edge] = Some(point);
            edge_routes[edge.edge] = self.edge_route(path, rank_y);
        }

        LayeredOutput {
            node_positions,
            edge_positions,
            edge_routes,
            ranks: rank_output,
        }
    }

    fn edge_route(&self, path: &[usize], rank_y: &[f64]) -> LayeredEdgeRoute {
        if path.len() <= 3 {
            return LayeredEdgeRoute::default();
        }
        let middle = path.len() / 2;
        let point = |item: usize| Point2::new(self.items[item].x, rank_y[self.items[item].rank]);
        LayeredEdgeRoute {
            source: path[1..middle].iter().map(|&item| point(item)).collect(),
            sink: path[middle + 1..path.len() - 1]
                .iter()
                .rev()
                .map(|&item| point(item))
                .collect(),
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum Direction {
    Up,
    Down,
}

#[derive(Debug, Clone, Copy)]
struct VisibleSegment {
    source: usize,
    sink: usize,
    source_rank: f64,
    source_order: f64,
    sink_rank: f64,
    sink_order: f64,
    weight: f64,
}

fn segments_cross(left: &VisibleSegment, right: &VisibleSegment) -> bool {
    let left_start = (left.source_rank, left.source_order);
    let left_end = (left.sink_rank, left.sink_order);
    let right_start = (right.source_rank, right.source_order);
    let right_end = (right.sink_rank, right.sink_order);
    let left_side_start = orient(left_start, left_end, right_start);
    let left_side_end = orient(left_start, left_end, right_end);
    let right_side_start = orient(right_start, right_end, left_start);
    let right_side_end = orient(right_start, right_end, left_end);
    left_side_start * left_side_end < 0.0 && right_side_start * right_side_end < 0.0
}

fn orient(a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> f64 {
    (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
}

fn find_disjoint_root(parent: &mut [usize], index: usize) -> usize {
    let mut root = index;
    while parent[root] != root {
        root = parent[root];
    }
    let mut current = index;
    while parent[current] != current {
        let next = parent[current];
        parent[current] = root;
        current = next;
    }
    root
}

fn intervals_overlap(a_start: f64, a_end: f64, b_start: f64, b_end: f64) -> bool {
    let a_min = a_start.min(a_end);
    let a_max = a_start.max(a_end);
    let b_min = b_start.min(b_end);
    let b_max = b_start.max(b_end);
    a_max >= b_min && b_max >= a_min
}
