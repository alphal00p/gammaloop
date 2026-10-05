//! Directed strongly connected components over the native half-edge graph.
use std::collections::{BTreeSet, HashMap, HashSet};

use crate::half_edge::{
    HedgeGraph, NodeIndex, involution::Flow, nodestore::NodeStorageOps, subgraph::SubSetLike,
};

use super::DirectionBasis;

impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    /// Strongly connected components in reverse topological order: an edge
    /// from one component to another points to an earlier returned component.
    /// Nodes within each component are sorted. Among available sink components,
    /// the one with smallest node index comes first, making unrelated blocks
    /// deterministic. Isolated nodes are included; identification aliases that
    /// share incidence are visited once.
    ///
    /// The iterative Tarjan traversal uses native half-edge crowns and requires
    /// no recursive call stack or copied adjacency graph. Ordering the resulting
    /// components also traverses native incidence. Complexity is O(V+E+V log V).
    pub fn strongly_connected_components(&self, basis: DirectionBasis) -> Vec<Vec<NodeIndex>> {
        let full = self.full_filter();
        let nodes = self.iter_nodes_of(&full).map(|(node, _, _)| node).chain(
            self.iter_nodes()
                .filter_map(|(node, mut crown, _)| crown.next().is_none().then_some(node)),
        );
        self.strongly_connected_nodes(&full, nodes, basis)
    }

    /// Components among nodes touched by `subgraph`. Both half-edges must be
    /// selected for an internal edge to participate. Identity half-edges and
    /// split-boundary edges are ignored. The superficial basis ignores edges
    /// whose superficial orientation is undirected.
    pub fn strongly_connected_components_of<S: SubSetLike>(
        &self,
        subgraph: &S,
        basis: DirectionBasis,
    ) -> Vec<Vec<NodeIndex>> {
        self.strongly_connected_nodes(
            subgraph,
            self.iter_nodes_of(subgraph).map(|(node, _, _)| node),
            basis,
        )
    }

    fn strongly_connected_nodes<S: SubSetLike>(
        &self,
        subgraph: &S,
        nodes: impl IntoIterator<Item = NodeIndex>,
        basis: DirectionBasis,
    ) -> Vec<Vec<NodeIndex>> {
        let nodes = nodes.into_iter().collect::<BTreeSet<_>>();
        let mut indices = HashMap::<NodeIndex, usize>::new();
        let mut low = HashMap::<NodeIndex, usize>::new();
        let mut active = HashSet::new();
        let mut component_stack = Vec::new();
        let mut components = Vec::new();
        let successor = |hedge, direction| {
            let inverse = self.inv(hedge);
            (subgraph.includes(&hedge)
                && inverse != hedge
                && subgraph.includes(&inverse)
                && basis.flow(self, hedge) == Some(direction))
            .then(|| self.node_id(inverse))
        };
        for &start in &nodes {
            if indices.contains_key(&start) {
                continue;
            }
            let index = indices.len();
            indices.insert(start, index);
            low.insert(start, index);
            active.insert(start);
            component_stack.push(start);
            let mut walk = vec![(start, self.iter_crown(start))];
            while let Some((node, crown)) = walk.last_mut() {
                let current = *node;
                if let Some(child) = crown.find_map(|hedge| successor(hedge, Flow::Source)) {
                    if let Some(&index) = indices.get(&child) {
                        if active.contains(&child) {
                            *low.get_mut(&current).unwrap() = low[&current].min(index);
                        }
                    } else {
                        let index = indices.len();
                        indices.insert(child, index);
                        low.insert(child, index);
                        active.insert(child);
                        component_stack.push(child);
                        walk.push((child, self.iter_crown(child)));
                    }
                } else {
                    walk.pop();
                    if low[&current] == indices[&current] {
                        let mut component = Vec::new();
                        loop {
                            let node = component_stack.pop().expect("active component root");
                            active.remove(&node);
                            component.push(node);
                            if node == current {
                                break;
                            }
                        }
                        component.sort();
                        components.push(component);
                    }
                    if let Some((parent, _)) = walk.last() {
                        *low.get_mut(parent).unwrap() = low[parent].min(low[&current]);
                    }
                }
            }
        }
        let owners = components
            .iter()
            .enumerate()
            .flat_map(|(i, c)| c.iter().map(move |&node| (node, i)))
            .collect::<HashMap<_, _>>();
        let mut pending = vec![0usize; components.len()];
        for (i, component) in components.iter().enumerate() {
            for &node in component {
                pending[i] += self
                    .iter_crown(node)
                    .filter_map(|h| successor(h, Flow::Source))
                    .filter(|child| owners[child] != i)
                    .count();
            }
        }
        let mut ready = components
            .iter()
            .enumerate()
            .filter_map(|(i, c)| (pending[i] == 0).then_some((c[0], i)))
            .collect::<BTreeSet<_>>();
        let mut result = Vec::with_capacity(components.len());
        while let Some((_, index)) = ready.pop_first() {
            for &node in &components[index] {
                for parent in self
                    .iter_crown(node)
                    .filter_map(|h| successor(h, Flow::Sink))
                {
                    let owner = owners[&parent];
                    if owner != index {
                        pending[owner] -= 1;
                        if pending[owner] == 0 {
                            ready.insert((components[owner][0], owner));
                        }
                    }
                }
            }
            result.push(components[index].clone());
        }
        debug_assert_eq!(result.len(), components.len());
        result
    }
}

#[cfg(test)]
mod tests {
    use super::DirectionBasis;
    use crate::half_edge::{
        HedgeGraph,
        builder::HedgeGraphBuilder,
        involution::{Flow, Orientation},
        subgraph::{ModifySubSet, SuBitGraph, SubSetLike},
    };

    #[test]
    fn directed_components_match_mutual_reachability_exhaustively() {
        for mask in 0..(1usize << 9) {
            let mut builder = HedgeGraphBuilder::<(), ()>::new();
            let nodes = (0..3).map(|_| builder.add_node(())).collect::<Vec<_>>();
            for i in 0..3 {
                for j in 0..3 {
                    if mask & (1 << (3 * i + j)) != 0 {
                        builder.add_edge(nodes[i], nodes[j], (), true);
                    }
                }
            }
            let graph: HedgeGraph<(), ()> = builder.build();
            let components = graph.strongly_connected_components(DirectionBasis::Underlying);
            let owner = |node| components.iter().position(|c| c.contains(&node)).unwrap();
            for &i in &nodes {
                for &j in &nodes {
                    let forward = graph.is_reachable(i, j, DirectionBasis::Underlying);
                    let reverse = graph.is_reachable(j, i, DirectionBasis::Underlying);
                    assert_eq!(owner(i) == owner(j), forward && reverse);
                    if forward {
                        assert!(owner(j) <= owner(i));
                    }
                }
            }
            let mut remaining = (0..components.len()).collect::<std::collections::BTreeSet<_>>();
            for component in &components {
                let eligible = remaining
                    .iter()
                    .copied()
                    .filter(|&i| {
                        remaining.iter().all(|&j| {
                            i == j
                                || !graph.is_reachable(
                                    components[i][0],
                                    components[j][0],
                                    DirectionBasis::Underlying,
                                )
                        })
                    })
                    .min_by_key(|&i| components[i][0])
                    .unwrap();
                assert_eq!(&components[eligible], component);
                remaining.remove(&eligible);
            }
        }
    }

    #[test]
    fn direction_split_boundaries_isolated_nodes_and_parallel_edges() {
        let mut builder = HedgeGraphBuilder::<(), ()>::new();
        let a = builder.add_node(());
        let b = builder.add_node(());
        let c = builder.add_node(());
        let isolated = builder.add_node(());
        builder.add_edge(a, b, (), Orientation::Default);
        builder.add_edge(a, b, (), Orientation::Reversed);
        builder.add_edge(c, b, (), Orientation::Undirected);
        builder.add_external_edge(a, (), Orientation::Default, Flow::Source);
        let graph: HedgeGraph<(), ()> = builder.build();
        assert_eq!(
            graph.strongly_connected_components(DirectionBasis::Underlying),
            vec![vec![b], vec![a], vec![c], vec![isolated]]
        );
        assert_eq!(
            graph.strongly_connected_components(DirectionBasis::Superficial),
            vec![vec![a, b], vec![c], vec![isolated]]
        );
        let mut selection = SuBitGraph::empty(graph.n_hedges());
        let hedge = graph.iter_crown(a).find(|&h| graph.inv(h) != h).unwrap();
        selection.add(hedge);
        assert_eq!(
            graph.strongly_connected_components_of(&selection, DirectionBasis::Underlying),
            vec![vec![a]]
        );
    }

    #[test]
    fn long_directed_cycle_does_not_recurse_on_the_call_stack() {
        let mut builder = HedgeGraphBuilder::<(), ()>::new();
        let nodes = (0..20_000)
            .map(|_| builder.add_node(()))
            .collect::<Vec<_>>();
        for i in 0..nodes.len() {
            builder.add_edge(nodes[i], nodes[(i + 1) % nodes.len()], (), true);
        }
        let graph: HedgeGraph<(), ()> = builder.build();
        assert_eq!(
            graph.strongly_connected_components(DirectionBasis::Underlying),
            vec![nodes]
        );
    }

    #[test]
    fn identification_history_is_not_counted_as_an_extra_component() {
        use crate::{
            half_edge::nodestore::NodeStorageVec,
            tree::{Forest, child_vec::ChildVecStore},
        };
        let mut builder = HedgeGraphBuilder::<(), ()>::new();
        let first = builder.add_node(());
        let historical = builder.add_node(());
        let target = builder.add_node(());
        let isolated = builder.add_node(());
        builder.add_edge(first, target, (), true);
        builder.add_edge(historical, target, (), true);
        let mut vector: HedgeGraph<_, _, _, NodeStorageVec<()>> = builder.clone().build();
        let merged = vector.identify_nodes(&[first, historical], ());
        let result = vector.strongly_connected_components(DirectionBasis::Underlying);
        assert_eq!(result.len(), 3);
        assert_eq!(
            result
                .iter()
                .flatten()
                .copied()
                .collect::<std::collections::BTreeSet<_>>(),
            [merged, target, isolated]
                .into_iter()
                .collect::<std::collections::BTreeSet<_>>()
        );
        assert!(
            result.iter().position(|c| c.contains(&target)).unwrap()
                < result.iter().position(|c| c.contains(&merged)).unwrap()
        );
        let mut forest: HedgeGraph<_, _, _, Forest<(), ChildVecStore<()>>> = builder.build();
        let merged = forest.identify_nodes(&[first, historical], ());
        let result = forest.strongly_connected_components(DirectionBasis::Underlying);
        assert_eq!(result.len(), 3);
        assert_eq!(
            result
                .iter()
                .flatten()
                .copied()
                .collect::<std::collections::BTreeSet<_>>(),
            [merged, target, isolated]
                .into_iter()
                .collect::<std::collections::BTreeSet<_>>()
        );
        assert!(
            result.iter().position(|c| c.contains(&target)).unwrap()
                < result.iter().position(|c| c.contains(&merged)).unwrap()
        );
    }
}
