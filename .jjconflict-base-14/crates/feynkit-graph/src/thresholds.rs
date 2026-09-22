//! Initial-state attachment trees used to normalize physical threshold sides.

use linnet::half_edge::{
    HedgeGraph,
    involution::{EdgeIndex, HedgePair},
    nodestore::NodeStorageOps,
    subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetOps},
};

use crate::{EdgeId, FeynmanDiagram};

/// Shared topology cleanup for threshold partitions after initial-state sewing.
pub trait InitialStateTreeExt {
    /// Select the vertex crowns adjacent to loop-independent paired edges.
    /// Threshold sides are normalized during FeynKit finalization; runtime
    /// threshold queries reuse the same selection without changing those sides.
    ///
    /// Prefer the edge's source crown when it touches the initial cut; otherwise
    /// use its sink crown. This is GammaLoop's existing threshold normalization,
    /// including its deterministic source/sink choice and complete crown selection.
    fn initial_state_tree(
        &self,
        initial_state: &SuBitGraph,
        loop_independent: impl FnMut(EdgeIndex) -> bool,
    ) -> (SuBitGraph, Vec<EdgeIndex>);
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> InitialStateTreeExt for HedgeGraph<E, V, H, N> {
    fn initial_state_tree(
        &self,
        initial_state: &SuBitGraph,
        mut loop_independent: impl FnMut(EdgeIndex) -> bool,
    ) -> (SuBitGraph, Vec<EdgeIndex>) {
        let mut tree_like_edges = Vec::new();
        let selected = self.full_filter().subtract(initial_state);
        let mut result = self.empty_subgraph::<SuBitGraph>();
        for (pair, edge, _) in self.iter_edges_of(&selected) {
            if let HedgePair::Paired { source, sink } = pair
                && loop_independent(edge)
            {
                tree_like_edges.push(edge);
                let source_node = self.node_id(source);
                let sink_node = self.node_id(sink);
                let source_connects_initial_state = self.iter_crown(source_node).any(|hedge| {
                    let mut one_hedge = self.empty_subgraph::<SuBitGraph>();
                    one_hedge.add(hedge);
                    one_hedge.intersects(initial_state)
                });
                let initial_node = if source_connects_initial_state {
                    source_node
                } else {
                    sink_node
                };
                for hedge in self.iter_crown(initial_node) {
                    result.add(hedge);
                }
            }
        }
        (result, tree_like_edges)
    }
}

impl FeynmanDiagram {
    /// Native half-edge selection removed from topology threshold sides.
    ///
    /// Sewn external connections define the initial state. The edge list records
    /// the loop-independent paired edges whose attached crowns form the selection.
    pub fn initial_state_tree(&self) -> (SuBitGraph, Vec<EdgeId>) {
        let mut initial_state = self.graph.empty_subgraph::<SuBitGraph>();
        for (pair, _, edge) in self.graph.iter_edges() {
            if edge.data.external.is_some()
                && let HedgePair::Paired { source, sink } = pair
            {
                initial_state.add(source);
                initial_state.add(sink);
            }
        }
        let (tree, edges) = self.graph.initial_state_tree(&initial_state, |edge| {
            self.loop_momentum_basis.edge_signatures[&EdgeId(edge.0)]
                .loops
                .iter()
                .all(|sign| sign.is_zero())
        });
        (tree, edges.into_iter().map(|edge| EdgeId(edge.0)).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use linnet::half_edge::{
        builder::HedgeGraphBuilder, involution::Orientation, subgraph::SubSetLike,
    };

    #[test]
    fn initial_attachment_tree_keeps_the_historical_crown_choice() {
        let mut builder = HedgeGraphBuilder::new();
        let a = builder.add_node(());
        let b = builder.add_node(());
        let c = builder.add_node(());
        let d = builder.add_node(());
        builder.add_edge(a, b, (), Orientation::Default); // sewn initial edge
        builder.add_edge(a, c, (), Orientation::Default); // source touches initial state
        builder.add_edge(c, d, (), Orientation::Default); // sink fallback
        let graph: HedgeGraph<(), ()> = builder.into();
        let mut initial: SuBitGraph = graph.empty_subgraph();
        for (pair, edge, _) in graph.iter_edges() {
            if edge == EdgeIndex(0)
                && let HedgePair::Paired { source, sink } = pair
            {
                initial.add(source);
                initial.add(sink);
            }
        }
        let (tree, edges) = graph.initial_state_tree(&initial, |_| true);
        let mut expected: SuBitGraph = graph.empty_subgraph();
        for node in [a, d] {
            for hedge in graph.iter_crown(node) {
                expected.add(hedge);
            }
        }
        assert_eq!(tree, expected);
        assert_eq!(edges, vec![EdgeIndex(1), EdgeIndex(2)]);
        let (empty, edges) = graph.initial_state_tree(&initial, |_| false);
        assert!(empty.is_empty());
        assert!(edges.is_empty());
    }
}
