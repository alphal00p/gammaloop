use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::half_edge::{
    involution::Flow, nodestore::NodeStorageOps, subgraph::SubSetLike, HedgeGraph, NodeIndex,
};
use thiserror::Error;

use super::DirectionBasis;

#[derive(Error, Debug)]
pub enum TopoError {
    #[error("Not a DAG: {nodes_processed} nodes processed out of {total_nodes} total nodes. Remaining nodes with non-zero in-degrees: {remaining_nodes:?}")]
    NotDag {
        nodes_processed: usize,
        total_nodes: usize,
        remaining_nodes: Vec<(NodeIndex, usize)>,
    },
}

impl<E, V, H, N: NodeStorageOps<NodeData = V>> HedgeGraph<E, V, H, N> {
    /// Returns a topological ordering of the nodes touched by `subgraph`.
    ///
    /// An internal edge contributes only when both half-edges are selected.
    /// Identity and split-boundary edges therefore do not affect the ordering.
    pub fn topo_sort_kahn_of<S: SubSetLike>(
        &self,
        subgraph: &S,
        basis: DirectionBasis,
    ) -> Result<Vec<NodeIndex>, TopoError> {
        let nodes = self
            .iter_nodes_of(subgraph)
            .map(|(node, _, _)| node)
            .collect::<BTreeSet<_>>();
        self.topo_sort_kahn_nodes(subgraph, nodes, basis)
    }

    /// Returns a topological ordering in the selected direction basis.
    /// Identification-history aliases that share incidence are processed once.
    pub fn topo_sort_kahn(&self, basis: DirectionBasis) -> Result<Vec<NodeIndex>, TopoError> {
        let full = self.full_filter();
        let structural_nodes = self.iter_nodes_of(&full).map(|(node, _, _)| node).chain(
            self.iter_nodes()
                .filter_map(|(node, mut crown, _)| crown.next().is_none().then_some(node)),
        );
        self.topo_sort_kahn_nodes(&full, structural_nodes, basis)
    }

    fn topo_sort_kahn_nodes<S: SubSetLike>(
        &self,
        subgraph: &S,
        nodes: impl IntoIterator<Item = NodeIndex>,
        basis: DirectionBasis,
    ) -> Result<Vec<NodeIndex>, TopoError> {
        let nodes = nodes.into_iter().collect::<BTreeSet<_>>();
        let mut indegrees = nodes
            .iter()
            .map(|&node| {
                let degree = self
                    .iter_crown(node)
                    .filter(|&hedge| {
                        if !subgraph.includes(&hedge) {
                            return false;
                        }
                        let inverse = self.inv(hedge);
                        inverse != hedge
                            && subgraph.includes(&inverse)
                            && basis.flow(self, hedge) == Some(Flow::Sink)
                    })
                    .count();
                (node, degree)
            })
            .collect::<BTreeMap<_, _>>();

        let mut queue = indegrees
            .iter()
            .filter_map(|(&node, &degree)| (degree == 0).then_some(node))
            .collect::<VecDeque<_>>();
        let mut order = Vec::with_capacity(indegrees.len());
        while let Some(node) = queue.pop_front() {
            order.push(node);
            for hedge in self.iter_crown(node) {
                if !subgraph.includes(&hedge) {
                    continue;
                }
                let inverse = self.inv(hedge);
                if inverse == hedge
                    || !subgraph.includes(&inverse)
                    || basis.flow(self, hedge) != Some(Flow::Source)
                {
                    continue;
                }

                let target = self.node_id(inverse);
                let degree = indegrees
                    .get_mut(&target)
                    .expect("both endpoints of an active edge are selected");
                *degree -= 1;
                if *degree == 0 {
                    queue.push_back(target);
                }
            }
        }

        if order.len() != indegrees.len() {
            let remaining_nodes = indegrees
                .into_iter()
                .filter(|(_, degree)| *degree > 0)
                .collect();

            return Err(TopoError::NotDag {
                nodes_processed: order.len(),
                total_nodes: nodes.len(),
                remaining_nodes,
            });
        }
        Ok(order)
    }
}

#[cfg(test)]
mod test {
    use super::TopoError;
    use crate::{
        dot,
        half_edge::{
            algorithms::DirectionBasis,
            builder::HedgeGraphBuilder,
            involution::{Flow, HedgePair, Orientation},
            subgraph::{ModifySubSet, SuBitGraph, SubSetLike},
            HedgeGraph,
        },
        parser::DotGraph,
    };

    #[test]
    fn topo_sort_valid_dag() {
        let graph: DotGraph = dot!(
            digraph ToughKahnTest {
              rankdir=LR;
              node [shape=circle];

              // Many independent sources (queue tie stress)
              s0; s1; s2; s3; s4; s5; s6; s7; s8; s9;

              // Two mid-layer diamonds (lots of merges)
              a0; a1; a2; a3; a4; a5; a6; a7;
              b0; b1; b2; b3; b4; b5; b6; b7;

              // A dense-ish join layer (high indegree nodes)
              j0; j1; j2; j3; j4;

              // Long chain tail (propagation stress)
              t0; t1; t2; t3; t4; t5; t6; t7; t8; t9; t10; t11;

              // Sources feed into two diamond structures
              s0 -> a0; s1 -> a0; s2 -> a1; s3 -> a1; s4 -> a2; s5 -> a2; s6 -> a3; s7 -> a3;
              s8 -> b0; s9 -> b0;

              // Diamond 1: a0/a1 -> a4, a2/a3 -> a5, then merge to a6 -> a7
              a0 -> a4; a1 -> a4;
              a2 -> a5; a3 -> a5;
              a4 -> a6; a5 -> a6;
              a6 -> a7;

              // Diamond 2: b0 splits into b1,b2 then rejoins, then splits again
              b0 -> b1; b0 -> b2;
              b1 -> b3; b2 -> b3;
              b3 -> b4; b3 -> b5;
              b4 -> b6; b5 -> b6;
              b6 -> b7;

              // Cross edges to create many simultaneous indegree updates
              a4 -> b4;
              a5 -> b5;
              b3 -> a6;
              b6 -> a7;

              // Join layer: high indegree nodes receiving from multiple predecessors
              a7 -> j0; b7 -> j0; s0 -> j0; s9 -> j0;
              a7 -> j1; b7 -> j1; s1 -> j1; s8 -> j1;
              a6 -> j2; b6 -> j2; s2 -> j2; s7 -> j2;
              a5 -> j3; b5 -> j3; s3 -> j3; s6 -> j3;
              a4 -> j4; b4 -> j4; s4 -> j4; s5 -> j4;

              // Join layer to tail chain
              j0 -> t0; j1 -> t0; j2 -> t1; j3 -> t1; j4 -> t2;
              t0 -> t3; t1 -> t3; t2 -> t4;
              t3 -> t5 -> t6 -> t7 -> t8 -> t9 -> t10 -> t11;

              // Extra edges that preserve DAG but increase constraint density
              a0 -> t0;
              b0 -> t1;
              a2 -> t2;
              b2 -> t2;
              j2 -> t4;
              j3 -> t4;

              // NOTE: Removed the cycle component to make this a valid DAG
            }
        )
        .unwrap();

        let order = graph
            .graph
            .topo_sort_kahn(DirectionBasis::Underlying)
            .unwrap();
        insta::assert_ron_snapshot!(order);
    }

    #[test]
    fn topo_sort_with_cycle() {
        let graph: DotGraph = dot!(
            digraph {
                A
                B
                C
                D
                A -> B
                B -> C
                C -> D
                D -> B  // Creates a cycle B -> C -> D -> B
            }
        )
        .unwrap();

        let result = graph.graph.topo_sort_kahn(DirectionBasis::Underlying);
        assert!(result.is_err());

        if let Err(TopoError::NotDag {
            nodes_processed,
            total_nodes,
            remaining_nodes,
        }) = result
        {
            // Should have processed some nodes but not all due to the cycle
            assert!(nodes_processed < total_nodes);
            assert!(!remaining_nodes.is_empty());
            // The nodes in the cycle should have non-zero in-degrees
            for (_, degree) in remaining_nodes {
                assert!(degree > 0);
            }
        }
    }

    #[test]
    fn topo_sort_complex_with_cycle() {
        let graph: DotGraph = dot!(
            digraph ToughKahnTestWithCycle {
              rankdir=LR;
              node [shape=circle];

              // Many independent sources (queue tie stress)
              s0; s1; s2; s3; s4; s5; s6; s7; s8; s9;

              // Two mid-layer diamonds (lots of merges)
              a0; a1; a2; a3; a4; a5; a6; a7;
              b0; b1; b2; b3; b4; b5; b6; b7;

              // A dense-ish join layer (high indegree nodes)
              j0; j1; j2; j3; j4;

              // Long chain tail (propagation stress)
              t0; t1; t2; t3; t4; t5; t6; t7; t8; t9; t10; t11;

              // A cycle component (should trigger NotDag error)
              c0; c1; c2; c3;

              // Sources feed into two diamond structures
              s0 -> a0; s1 -> a0; s2 -> a1; s3 -> a1; s4 -> a2; s5 -> a2; s6 -> a3; s7 -> a3;
              s8 -> b0; s9 -> b0;

              // Diamond 1: a0/a1 -> a4, a2/a3 -> a5, then merge to a6 -> a7
              a0 -> a4; a1 -> a4;
              a2 -> a5; a3 -> a5;
              a4 -> a6; a5 -> a6;
              a6 -> a7;

              // Diamond 2: b0 splits into b1,b2 then rejoins, then splits again
              b0 -> b1; b0 -> b2;
              b1 -> b3; b2 -> b3;
              b3 -> b4; b3 -> b5;
              b4 -> b6; b5 -> b6;
              b6 -> b7;

              // Cross edges to create many simultaneous indegree updates
              a4 -> b4;
              a5 -> b5;
              b3 -> a6;
              b6 -> a7;

              // Join layer: high indegree nodes receiving from multiple predecessors
              a7 -> j0; b7 -> j0; s0 -> j0; s9 -> j0;
              a7 -> j1; b7 -> j1; s1 -> j1; s8 -> j1;
              a6 -> j2; b6 -> j2; s2 -> j2; s7 -> j2;
              a5 -> j3; b5 -> j3; s3 -> j3; s6 -> j3;
              a4 -> j4; b4 -> j4; s4 -> j4; s5 -> j4;

              // Join layer to tail chain
              j0 -> t0; j1 -> t0; j2 -> t1; j3 -> t1; j4 -> t2;
              t0 -> t3; t1 -> t3; t2 -> t4;
              t3 -> t5 -> t6 -> t7 -> t8 -> t9 -> t10 -> t11;

              // Extra edges that preserve DAG but increase constraint density
              a0 -> t0;
              b0 -> t1;
              a2 -> t2;
              b2 -> t2;
              j2 -> t4;
              j3 -> t4;

              // Cycle component (makes graph non-DAG)
              c0 -> c1 -> c2 -> c3 -> c0;
            }
        )
        .unwrap();

        let result = graph.graph.topo_sort_kahn(DirectionBasis::Underlying);
        assert!(result.is_err());

        if let Err(TopoError::NotDag {
            nodes_processed,
            total_nodes,
            remaining_nodes,
        }) = result
        {
            println!("Complex graph with cycle error details:");
            println!(
                "  Processed {} out of {} nodes",
                nodes_processed, total_nodes
            );
            println!("  Remaining nodes with non-zero in-degrees:");
            for (node, degree) in &remaining_nodes {
                println!("    {:?}: degree {}", node, degree);
            }

            // Should have processed most nodes except those in the cycle
            assert!(nodes_processed > 0);
            assert!(nodes_processed < total_nodes);
            assert!(!remaining_nodes.is_empty());
            // The cycle should contain exactly 4 nodes (c0, c1, c2, c3)
            assert_eq!(remaining_nodes.len(), 4);
            // All nodes in the cycle should have in-degree 1
            for (_, degree) in remaining_nodes {
                assert_eq!(degree, 1);
            }
        }
    }

    #[test]
    fn direction_basis_distinguishes_reversed_and_undirected_edges() {
        let mut reversed = HedgeGraphBuilder::<(), ()>::new();
        let left = reversed.add_node(());
        let right = reversed.add_node(());
        reversed.add_edge(left, right, (), Orientation::Reversed);
        let reversed: HedgeGraph<(), ()> = reversed.build();

        assert_eq!(
            reversed.topo_sort_kahn(DirectionBasis::Underlying).unwrap(),
            vec![left, right]
        );
        assert_eq!(
            reversed
                .topo_sort_kahn(DirectionBasis::Superficial)
                .unwrap(),
            vec![right, left]
        );

        let mut undirected = HedgeGraphBuilder::<(), ()>::new();
        let target = undirected.add_node(());
        let source = undirected.add_node(());
        undirected.add_edge(source, target, (), Orientation::Undirected);
        let undirected: HedgeGraph<(), ()> = undirected.build();

        assert_eq!(
            undirected
                .topo_sort_kahn(DirectionBasis::Underlying)
                .unwrap(),
            vec![source, target]
        );
        assert_eq!(
            undirected
                .topo_sort_kahn(DirectionBasis::Superficial)
                .unwrap(),
            vec![target, source]
        );
    }

    #[test]
    fn dangling_and_split_boundary_edges_do_not_add_indegree() {
        let mut builder = HedgeGraphBuilder::<(), ()>::new();
        let first = builder.add_node(());
        let middle = builder.add_node(());
        let last = builder.add_node(());
        builder.add_edge(first, middle, (), Orientation::Default);
        builder.add_edge(middle, last, (), Orientation::Default);
        builder.add_external_edge(middle, (), Orientation::Default, Flow::Sink);
        let graph: HedgeGraph<(), ()> = builder.build();

        assert_eq!(
            graph.topo_sort_kahn(DirectionBasis::Underlying).unwrap(),
            vec![first, middle, last]
        );

        let first_pair = graph.iter_edges().next().unwrap().0;
        let second_pair = graph.iter_edges().nth(1).unwrap().0;
        let HedgePair::Paired { sink, .. } = first_pair else {
            panic!("builder created a non-paired edge")
        };
        let mut selected = SuBitGraph::empty(graph.n_hedges());
        selected.add(sink);
        selected.add(second_pair);

        assert_eq!(
            graph
                .topo_sort_kahn_of(&selected, DirectionBasis::Underlying)
                .unwrap(),
            vec![middle, last]
        );
    }

    #[test]
    fn cycles_are_evaluated_in_the_selected_direction_basis() {
        let mut builder = HedgeGraphBuilder::<(), ()>::new();
        let first = builder.add_node(());
        let second = builder.add_node(());
        builder.add_edge(first, second, (), Orientation::Default);
        builder.add_edge(second, first, (), Orientation::Undirected);
        let graph: HedgeGraph<(), ()> = builder.build();

        assert!(graph.topo_sort_kahn(DirectionBasis::Underlying).is_err());
        assert_eq!(
            graph.topo_sort_kahn(DirectionBasis::Superficial).unwrap(),
            vec![first, second]
        );
    }

    #[test]
    fn full_sort_canonicalizes_identification_history() {
        use crate::{
            half_edge::nodestore::NodeStorageVec,
            tree::{child_vec::ChildVecStore, Forest},
        };

        let mut builder = HedgeGraphBuilder::<(), ()>::new();
        let first = builder.add_node(());
        let historical = builder.add_node(());
        let target = builder.add_node(());
        let isolated = builder.add_node(());
        builder.add_edge(first, target, (), Orientation::Default);
        builder.add_edge(historical, target, (), Orientation::Default);
        let mut vec_graph: HedgeGraph<_, _, _, NodeStorageVec<()>> = builder.clone().build();
        let merged = vec_graph.identify_nodes(&[first, historical], ());

        assert_eq!(
            vec_graph
                .topo_sort_kahn(DirectionBasis::Underlying)
                .unwrap(),
            vec![merged, isolated, target]
        );

        let mut forest_graph: HedgeGraph<_, _, _, Forest<(), ChildVecStore<()>>> = builder.build();
        let merged = forest_graph.identify_nodes(&[first, historical], ());
        assert_eq!(
            forest_graph
                .topo_sort_kahn(DirectionBasis::Underlying)
                .unwrap(),
            vec![merged, isolated, target]
        );
    }
}
