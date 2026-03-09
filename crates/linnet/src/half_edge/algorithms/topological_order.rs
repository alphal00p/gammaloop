use std::collections::VecDeque;

use crate::half_edge::{
    involution::Flow, nodestore::NodeStorageOps, swap::Swap, HedgeGraph, NodeIndex,
};
use thiserror::Error;

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
    pub fn topo_sort_kahn(&self) -> Result<Vec<NodeIndex>, TopoError> {
        let mut indeg = self
            .new_nodevec(|_i, neighs, _| neighs.filter(|h| self.flow(*h) == Flow::Sink).count());

        let mut q = VecDeque::new();
        for (i, d) in indeg.iter() {
            if *d == 0 {
                q.push_back(i);
            }
        }

        let mut order = Vec::with_capacity(indeg.len().0);
        while let Some(v) = q.pop_front() {
            order.push(v);
            for u in self.iter_crown(v) {
                if self.flow(u) == Flow::Sink {
                    continue;
                }
                let Some(u) = self.involved_node_id(u) else {
                    continue;
                };
                indeg[u] -= 1;
                if indeg[u] == 0 {
                    q.push_back(u);
                }
            }
        }

        if order.len() != indeg.len().0 {
            let remaining_nodes: Vec<(NodeIndex, usize)> = indeg
                .iter()
                .filter_map(|(node, &degree)| {
                    if degree > 0 {
                        Some((node, degree))
                    } else {
                        None
                    }
                })
                .collect();

            return Err(TopoError::NotDag {
                nodes_processed: order.len(),
                total_nodes: indeg.len().0,
                remaining_nodes,
            });
        }
        Ok(order)
    }
}

#[cfg(test)]
mod test {
    use super::TopoError;
    use crate::{dot, parser::DotGraph};

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

        let order = graph.graph.topo_sort_kahn().unwrap();
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

        let result = graph.graph.topo_sort_kahn();
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

        let result = graph.graph.topo_sort_kahn();
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
}
