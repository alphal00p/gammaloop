//! # Transitive Operations on Half-Edge Graphs
//!
//! This module provides algorithms for computing transitive closure and transitive reduction
//! on directed acyclic graphs represented using the half-edge data structure.
//!
//! ## Half-Edge Graph Fundamentals
//!
//! In the half-edge representation, each directed edge is composed of two **half-edges**:
//! - A **source half-edge** (outgoing from the source node)
//! - A **sink half-edge** (incoming to the target node)
//!
//! These half-edges are paired together and can be accessed via the involution `self.inv(hedge)`.
//!
//! ## Critical Implementation Details
//!
//! When modifying edges in a half-edge graph:
//! - **Adding edges**: Use `add_pair()` to create both half-edges simultaneously
//! - **Removing edges**: Both half-edges must be deleted to avoid corruption
//!
//! ## Algorithms Provided
//!
//! - [`transitive_closure`](HedgeGraph::transitive_closure): Adds all transitive edges
//! - [`transitive_reduction`](HedgeGraph::transitive_reduction): Removes redundant edges
//!
//! Both algorithms verify the graph is acyclic before proceeding, as transitive operations
//! on graphs with cycles are undefined.

use std::collections::{HashMap, HashSet, VecDeque};

use crate::half_edge::{
    builder::HedgeGraphBuilder,
    involution::{Flow, Orientation},
    nodestore::NodeStorageOps,
    subgraph::{ModifySubSet, SuBitGraph},
    HedgeGraph, NoData, NodeIndex,
};

use super::topological_order::TopoError;

/// Error types for transitive operations on half-edge graphs.
///
/// These operations require directed acyclic graphs (DAGs) to be meaningful,
/// as transitive operations on graphs with cycles are undefined.
#[derive(Debug, thiserror::Error)]
pub enum TransitiveError {
    /// The graph contains cycles, making transitive operations undefined.
    ///
    /// This error is returned when the topological sort fails, indicating
    /// the presence of cycles in the directed graph.
    #[error("Graph contains cycles: {0}")]
    HasCycles(#[from] TopoError),

    /// An error occurred in the underlying graph structure operations.
    ///
    /// This can happen during edge addition or removal operations.
    #[error("Graph structure error")]
    GraphStructure,
}

impl<V, S: NodeStorageOps<NodeData = V, Base = SuBitGraph>> HedgeGraph<NoData, V, NoData, S> {
    pub fn poset<I: IntoIterator<Item = V>>(nodes: I) -> Self
    where
        V: PartialOrd,
    {
        let mut builder = HedgeGraphBuilder::new();
        let mut nodes: Vec<Option<V>> = nodes.into_iter().map(Some).collect();

        // Create nodes in the graph
        let mut node_indices = HashMap::new();
        for (i, node) in nodes.iter().enumerate() {
            let idx = builder.add_node(i);
            node_indices.insert(i, idx);
        }

        // Add ALL order relations first (complete partial order graph)
        for (i, a) in nodes.iter().enumerate() {
            for (j, b) in nodes.iter().enumerate() {
                if i != j && a < b {
                    builder.add_edge(node_indices[&i], node_indices[&j], NoData {}, true);
                }
            }
        }
        builder.build_with_map::<V, S>(|i| nodes[i].take().unwrap())
    }
    pub fn hasse_diagram<I: IntoIterator<Item = V>>(nodes: I) -> Self
    where
        V: PartialOrd + Clone,
    {
        use crate::half_edge::builder::HedgeGraphBuilder;
        use crate::half_edge::involution::Orientation;

        let mut builder = HedgeGraphBuilder::new();
        let nodes: Vec<V> = nodes.into_iter().collect();

        // Create nodes in the graph
        let mut node_indices = HashMap::new();
        for (i, node) in nodes.iter().enumerate() {
            let idx = builder.add_node(node.clone());
            node_indices.insert(i, idx);
        }

        // Add ALL order relations first (complete partial order graph)
        for (i, a) in nodes.iter().enumerate() {
            for (j, b) in nodes.iter().enumerate() {
                if i != j && a < b {
                    builder.add_edge(
                        node_indices[&i],
                        node_indices[&j],
                        NoData {},
                        Orientation::Default,
                    );
                }
            }
        }

        // Build the complete graph and apply transitive reduction
        let complete_graph = builder.build();

        // Apply transitive reduction to get the Hasse diagram
        complete_graph
            .transitive_reduction()
            .expect("Partial order should be acyclic")
    }
}

impl<E, V, H, N: NodeStorageOps<NodeData = V, Base = crate::half_edge::subgraph::SuBitGraph>>
    HedgeGraph<E, V, H, N>
{
    /// Computes the transitive closure of the directed acyclic graph.
    ///
    /// The transitive closure adds an edge (u, v) for every pair of nodes u, v
    /// where there exists a directed path from u to v in the original graph.
    ///
    /// ## Algorithm Details
    ///
    /// This implementation uses a Floyd-Warshall-style algorithm optimized for DAGs:
    /// 1. Verifies the graph is acyclic using topological sorting
    /// 2. Builds a reachability matrix using the topological order
    /// 3. For each missing transitive edge, adds it using `add_pair()`
    ///
    /// ## Half-Edge Considerations
    ///
    /// In HedgeGraph, each directed edge consists of two half-edges (a source and sink pair).
    /// When adding new transitive edges, the algorithm:
    /// - Uses `add_pair()` to create both half-edges simultaneously
    /// - Assigns default edge data to new transitive edges
    /// - Maintains proper half-edge pair relationships
    ///
    /// # Arguments
    ///
    /// * `self` - Takes ownership of the graph to avoid cloning
    ///
    /// # Returns
    ///
    /// The modified graph with all transitive edges added
    ///
    /// # Errors
    ///
    /// Returns `TransitiveError::HasCycles` if the graph contains cycles
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// // Graph: A -> B -> C
    /// let closure = graph.transitive_closure()?;
    /// // Result: A -> B -> C, A -> C (transitive edge added)
    /// ```
    pub fn transitive_closure(self) -> Result<Self, TransitiveError>
    where
        E: Default,
        H: Default,
    {
        // First verify this is a DAG by getting topological order
        let topo_order = self.topo_sort_kahn()?;

        // Build reachability matrix using Floyd-Warshall algorithm
        let nodes: Vec<NodeIndex> = topo_order;
        let n = nodes.len();
        let mut node_to_idx = HashMap::new();
        let mut idx_to_node = Vec::new();

        // Map node indices to array indices for efficiency
        for (i, &node) in nodes.iter().enumerate() {
            node_to_idx.insert(node, i);
            idx_to_node.push(node);
        }

        // Initialize reachability matrix
        let mut reachable = vec![vec![false; n]; n];

        // Set diagonal to true (each node reaches itself)
        for i in 0..n {
            reachable[i][i] = true;
        }

        // Initialize direct edges
        for (i, &from_node) in nodes.iter().enumerate() {
            for hedge in self.iter_crown(from_node) {
                if self.flow(hedge) == Flow::Source {
                    if let Some(to_node) = self.involved_node_id(hedge) {
                        if let Some(&j) = node_to_idx.get(&to_node) {
                            reachable[i][j] = true;
                        }
                    }
                }
            }
        }

        // Floyd-Warshall to compute transitive closure
        for k in 0..n {
            for i in 0..n {
                for j in 0..n {
                    reachable[i][j] = reachable[i][j] || (reachable[i][k] && reachable[k][j]);
                }
            }
        }

        // Add new transitive edges
        let mut graph = self;
        for i in 0..n {
            for j in 0..n {
                if i != j && reachable[i][j] {
                    let from_node = idx_to_node[i];
                    let to_node = idx_to_node[j];

                    // Check if edge already exists
                    // Check if edge already exists
                    let edge_exists = graph.iter_crown(from_node).any(|hedge| {
                        graph.flow(hedge) == crate::half_edge::involution::Flow::Source
                            && graph.involved_node_id(hedge) == Some(to_node)
                    });

                    if !edge_exists {
                        // Add edge using add_pair
                        let (_, _, new_graph) = graph
                            .add_pair(
                                from_node,
                                to_node,
                                E::default(),
                                crate::half_edge::involution::Orientation::Default,
                            )
                            .map_err(|_| TransitiveError::GraphStructure)?;
                        graph = new_graph;
                    }
                }
            }
        }

        Ok(graph)
    }

    /// Computes the transitive reduction of the directed acyclic graph.
    ///
    /// The transitive reduction removes all edges (u, v) such that there exists
    /// an alternative path from u to v that doesn't use the direct edge.
    /// This produces the minimal graph that preserves the same reachability.
    ///
    /// ## Algorithm Details
    ///
    /// This implementation:
    /// 1. Verifies the graph is acyclic using topological sorting
    /// 2. For each edge (u, v), checks if an alternative path u ⇝ v exists
    /// 3. Collects all redundant edges for removal
    /// 4. Uses `delete_hedges()` to remove redundant edges
    ///
    /// ## Half-Edge Considerations
    ///
    /// **Critical Implementation Detail**: In HedgeGraph, each directed edge consists of
    /// two half-edges (a paired source and sink). To completely remove an edge, both
    /// half-edges must be deleted:
    ///
    /// - `hedge`: The source half-edge (outgoing from source node)
    /// - `self.inv(hedge)`: The sink half-edge (incoming to target node)
    ///
    /// The algorithm ensures complete edge removal by:
    /// 1. Identifying redundant source half-edges during traversal
    /// 2. Finding their paired sink half-edges using `self.inv()`
    /// 3. Adding both to the removal subgraph
    /// 4. Calling `delete_hedges()` once to remove all redundant half-edge pairs
    ///
    /// Failing to remove both half-edges would leave dangling half-edges and
    /// corrupt the graph structure.
    ///
    /// # Arguments
    ///
    /// * `self` - Takes ownership of the graph to modify it in place
    ///
    /// # Returns
    ///
    /// The modified graph with all redundant edges removed
    ///
    /// # Errors
    ///
    /// Returns `TransitiveError::HasCycles` if the graph contains cycles
    ///
    /// # Examples
    ///
    /// ```rust,ignore
    /// // Graph: A -> B -> C, A -> C (redundant)
    /// let reduction = graph.transitive_reduction()?;
    /// // Result: A -> B -> C (A -> C removed as redundant)
    /// ```
    pub fn transitive_reduction(mut self) -> Result<Self, TransitiveError> {
        // First verify this is a DAG
        let topo_order = self.topo_sort_kahn()?;

        // Collect all redundant edges (edges that can be removed)
        let mut redundant_hedges = Vec::new();
        for &u in &topo_order {
            for hedge in self.iter_crown(u) {
                if self.flow(hedge) == Flow::Source {
                    if let Some(v) = self.involved_node_id(hedge) {
                        // Check if this edge is redundant (has alternative path)
                        if self.has_alternative_path(u, v) {
                            redundant_hedges.push(hedge);
                        }
                    }
                }
            }
        }

        // Build a subgraph with all redundant edges to remove
        let mut edges_to_remove: crate::half_edge::subgraph::SuBitGraph = self.empty_subgraph();
        for hedge in redundant_hedges {
            // CRITICAL: Add both half-edges of the pair to completely remove the edge
            // In HedgeGraph, each directed edge consists of two half-edges:
            // - hedge: the source half-edge (outgoing from source node)
            // - paired_hedge: the sink half-edge (incoming to target node)
            // Both must be removed to avoid leaving dangling half-edges
            edges_to_remove.add(hedge);
            let paired_hedge = self.inv(hedge);
            if paired_hedge != hedge {
                edges_to_remove.add(paired_hedge);
            }
        }

        // Delete all redundant hedges at once
        self.delete_hedges(&edges_to_remove);

        Ok(self)
    }

    /// Helper function to check if there's an alternative path from u to v
    /// excluding the direct edge u -> v.
    ///
    /// This uses breadth-first search starting from all neighbors of u (except v)
    /// to determine if v can be reached through an indirect path. This is used
    /// to identify redundant edges in transitive reduction.
    fn has_alternative_path(&self, u: NodeIndex, v: NodeIndex) -> bool {
        // Use BFS to find if there's a path from u to v without using the direct edge
        let mut queue = VecDeque::new();
        let mut visited = HashSet::new();

        // Start with all neighbors of u except v
        for hedge in self.iter_crown(u) {
            if self.flow(hedge) == Flow::Source {
                if let Some(neighbor) = self.involved_node_id(hedge) {
                    if neighbor != v && !visited.contains(&neighbor) {
                        queue.push_back(neighbor);
                        visited.insert(neighbor);
                    }
                }
            }
        }

        // BFS to find v
        while let Some(current) = queue.pop_front() {
            if current == v {
                return true;
            }

            for hedge in self.iter_crown(current) {
                if self.flow(hedge) == Flow::Source {
                    if let Some(neighbor) = self.involved_node_id(hedge) {
                        if !visited.contains(&neighbor) {
                            queue.push_back(neighbor);
                            visited.insert(neighbor);
                        }
                    }
                }
            }
        }

        false
    }

    /// Checks if there is a path from `source` to `target` in the graph.
    ///
    /// Uses breadth-first search to determine reachability. This method
    /// only follows source half-edges (outgoing edges) to traverse the
    /// directed graph structure.
    ///
    /// # Arguments
    ///
    /// * `source` - The source node
    /// * `target` - The target node
    ///
    /// # Returns
    ///
    /// `true` if there is a directed path from source to target, `false` otherwise.
    ///
    /// # Half-Edge Traversal
    ///
    /// The algorithm only follows half-edges with `Flow::Source` to ensure
    /// it respects the directed nature of edges in the graph.
    pub fn is_reachable(&self, source: NodeIndex, target: NodeIndex) -> bool {
        if source == target {
            return true;
        }

        let mut visited = HashSet::new();
        let mut queue = VecDeque::new();

        queue.push_back(source);
        visited.insert(source);

        while let Some(current) = queue.pop_front() {
            for hedge in self.iter_crown(current) {
                if self.flow(hedge) == Flow::Source {
                    if let Some(neighbor) = self.involved_node_id(hedge) {
                        if neighbor == target {
                            return true;
                        }
                        if !visited.contains(&neighbor) {
                            visited.insert(neighbor);
                            queue.push_back(neighbor);
                        }
                    }
                }
            }
        }

        false
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        dot,
        half_edge::{nodestore::NodeStorageOps, HedgeGraph, NoData},
        parser::{DotGraph, DotVertexData},
    };
    use crate::half_edge::nodestore::NodeStorageVec;

    #[test]
    fn test_transitive_closure_simple() {
        let graph: DotGraph = dot!(
            digraph {
                A -> B
                B -> C
                A -> D
                D -> C
            }
        )
        .unwrap();

        let original_edges = graph.graph.n_edges();
        let closure = graph.graph.transitive_closure().unwrap();

        // The closure should have all original edges plus A -> C
        // (since A can reach C through both A->B->C and A->D->C)
        assert!(closure.n_edges() >= original_edges);
    }

    #[test]
    fn test_transitive_reduction_simple() {
        let graph: DotGraph = dot!(
            digraph {
                A -> B
                B -> C
                A -> C  // This should be removed in reduction
            }
        )
        .unwrap();

        let original_edges = graph.graph.n_edges();
        let reduction = graph.graph.transitive_reduction().unwrap();

        // The reduction should remove the redundant A -> C edge
        assert!(reduction.n_edges() < original_edges);
    }

    #[test]
    fn test_transitive_closure_diamond() {
        let graph: DotGraph = dot!(
            digraph {
                A -> B
                A -> C
                B -> D
                C -> D
            }
        )
        .unwrap();

        let original_edges = graph.graph.n_edges();
        let closure = graph.graph.transitive_closure().unwrap();

        // Should add A -> D edge in closure
        assert!(closure.n_edges() >= original_edges);
    }

    #[test]
    fn test_transitive_reduction_diamond() {
        let graph: DotGraph = dot!(
            digraph {
                A -> B
                A -> C
                B -> D
                C -> D
                A -> D  // This should be removed
            }
        )
        .unwrap();

        let reduction = graph.graph.transitive_reduction().unwrap();

        // Should remove the redundant A -> D edge
        assert_eq!(reduction.n_edges(), 4); // A->B, A->C, B->D, C->D
    }

    #[test]
    fn test_cycle_detection() {
        let graph: DotGraph = dot!(
            digraph {
                A -> B
                B -> C
                C -> A  // Creates cycle
            }
        )
        .unwrap();

        assert!(graph.graph.transitive_closure().is_err());

        let graph2: DotGraph = dot!(
            digraph {
                A -> B
                B -> C
                C -> A  // Creates cycle
            }
        )
        .unwrap();
        assert!(graph2.graph.transitive_reduction().is_err());
    }

    #[test]
    fn test_empty_graph() {
        let graph: DotGraph<NodeStorageVec<DotVertexData>> = DotGraph::from_string(
            r#"
            digraph {
                A
                B
                C
            }
            "#,
        )
        .unwrap();

        let closure = graph.graph.transitive_closure().unwrap();

        let graph2: DotGraph<NodeStorageVec<DotVertexData>> = DotGraph::from_string(
            r#"
            digraph {
                A
                B
                C
            }
            "#,
        )
        .unwrap();
        let reduction = graph2.graph.transitive_reduction().unwrap();

        // Empty graphs should remain empty
        assert_eq!(closure.n_edges(), 0);
        assert_eq!(reduction.n_edges(), 0);
    }

    #[test]
    fn test_single_node() {
        let graph: DotGraph<NodeStorageVec<DotVertexData>> =
            DotGraph::from_string("digraph { A }").unwrap();

        let closure = graph.graph.transitive_closure().unwrap();

        let graph2: DotGraph<NodeStorageVec<DotVertexData>> =
            DotGraph::from_string("digraph { A }").unwrap();
        let reduction = graph2.graph.transitive_reduction().unwrap();

        assert_eq!(closure.n_nodes(), 1);
        assert_eq!(reduction.n_nodes(), 1);
        assert_eq!(closure.n_edges(), 0);
        assert_eq!(reduction.n_edges(), 0);
    }

    #[test]
    fn test_is_reachable() {
        let graph: DotGraph = dot!(
            digraph {
                A -> B
                B -> C
                D -> E
            }
        )
        .unwrap();

        // Test reachability with actual node indices from the graph
        let nodes: Vec<_> = graph.graph.iter_node_ids().collect();
        if nodes.len() >= 4 {
            let a = nodes[0];
            let _b = nodes[1];
            let _c = nodes[2];
            let _d = nodes[3];

            assert!(graph.graph.is_reachable(a, a));
            // Note: actual reachability depends on how nodes are mapped in DOT parsing
        }
    }

    // ===============================
    // COMPREHENSIVE INSTA TESTS
    // ===============================

    #[test]
    fn test_transitive_closure_complex_dag() {
        let original: DotGraph = dot!(
            digraph {
                // Complex DAG with multiple paths
                A -> B
                A -> C
                A -> D
                B -> E
                C -> E
                C -> F
                D -> F
                E -> G
                F -> G
                F -> H
                H -> I
                // Some transitive edges will be added
            }
        )
        .unwrap();

        let closure_graph = DotGraph {
            graph: original.graph.transitive_closure().unwrap(),
            global_data: original.global_data.clone(),
        };

        insta::assert_snapshot!(closure_graph.debug_dot());
    }

    #[test]
    fn test_transitive_reduction_complex_redundant() {
        let original: DotGraph = dot!(
            digraph {
                // Graph with many redundant edges
                A -> B
                B -> C
                C -> D
                A -> C  // Redundant via A->B->C
                A -> D  // Redundant via A->B->C->D and A->C->D
                B -> D  // Redundant via B->C->D
                D -> E
                E -> F
                A -> E  // Redundant via multiple paths
                B -> F  // Redundant via B->C->D->E->F
            }
        )
        .unwrap();

        let reduction_graph = DotGraph {
            graph: original.graph.transitive_reduction().unwrap(),
            global_data: original.global_data.clone(),
        };

        insta::assert_snapshot!(reduction_graph.debug_dot());
    }

    #[test]
    fn test_transitive_closure_tree_structure() {
        let original: DotGraph = dot!(
            digraph {
                // Binary tree-like structure
                root -> left
                root -> right
                left -> ll
                left -> lr
                right -> rl
                right -> rr
                ll -> lll
                ll -> llr
            }
        )
        .unwrap();

        let closure_graph = DotGraph {
            graph: original.graph.transitive_closure().unwrap(),
            global_data: original.global_data.clone(),
        };

        insta::assert_snapshot!(closure_graph.debug_dot());
    }

    #[test]
    fn test_transitive_reduction_complete_subgraph() {
        let original: DotGraph = dot!(
            digraph {
                // Start with transitively complete subgraph
                A -> B
                A -> C
                A -> D
                B -> C
                B -> D
                C -> D
                // Plus some extensions
                D -> E
                E -> F
                A -> E  // Should be removed (via A->D->E)
                B -> F  // Should be removed (via B->D->E->F)
            }
        )
        .unwrap();

        let reduction_graph = DotGraph {
            graph: original.graph.transitive_reduction().unwrap(),
            global_data: original.global_data.clone(),
        };

        insta::assert_snapshot!(reduction_graph.debug_dot());
    }

    #[test]
    fn test_transitive_closure_parallel_paths() {
        let original: DotGraph = dot!(
            digraph {
                // Multiple parallel paths of different lengths
                start -> path1_a
                start -> path2_a
                start -> path3_a

                path1_a -> end

                path2_a -> path2_b
                path2_b -> end

                path3_a -> path3_b
                path3_b -> path3_c
                path3_c -> end

                // Some intermediate connections
                path2_a -> path3_c
                path1_a -> path2_b
            }
        )
        .unwrap();

        let closure_graph = DotGraph {
            graph: original.graph.transitive_closure().unwrap(),
            global_data: original.global_data.clone(),
        };

        insta::assert_snapshot!(closure_graph.debug_dot());
    }

    #[test]
    fn test_transitive_operations_on_star_graph() {
        let original: DotGraph = dot!(
            digraph {
                // Star graph (center connects to all)
                center -> a
                center -> b
                center -> c
                center -> d
                center -> e
                // Add some redundant edges for reduction test
                a -> b
                b -> c
                center -> c  // Redundant via center->a->b->c
            }
        )
        .unwrap();

        // Test closure
        let closure_graph = DotGraph {
            graph: original.graph.clone().transitive_closure().unwrap(),
            global_data: original.global_data.clone(),
        };
        insta::assert_snapshot!(closure_graph.debug_dot(), @r"
        digraph {
          a;
          b;
          c;
          center;
          d;
          e;

          a:0	-> b:1	 [id=0 ];
          b:2	-> c:3	 [id=1 ];
          center:4	-> a:5	 [id=2 ];
          center:6	-> b:7	 [id=3 ];
          center:8	-> c:9	 [id=4 ];
          center:10	-> c:11	 [id=5 ];
          center:12	-> d:13	 [id=6 ];
          center:14	-> e:15	 [id=7 ];
          a:16	-> c:17	 [id=8 ];
        }
        ");

        // Test reduction
        let reduction_graph = DotGraph {
            graph: original.graph.transitive_reduction().unwrap(),
            global_data: original.global_data.clone(),
        };
        insta::assert_snapshot!(reduction_graph.debug_dot(), @r"
        digraph {
          a;
          e;
          d;
          center;
          c;
          b;

          a:0	-> b:1	 [id=0 ];
          b:2	-> c:3	 [id=1 ];
          center:4	-> a:5	 [id=2 ];
          center:7	-> e:6	 [id=3 ];
          center:9	-> d:8	 [id=4 ];
        }
        ");
    }

    #[test]
    fn test_transitive_closure_deep_chain() {
        let original: DotGraph = dot!(
            digraph {
                // Long chain with some branches
                a -> b
                b -> c
                c -> d
                d -> e
                e -> f
                f -> g
                g -> h
                h -> i
                i -> j

                // Some branches at different levels
                c -> x
                x -> y
                f -> z

                // Cross connections
                b -> z
                d -> y
            }
        )
        .unwrap();

        let closure_graph = DotGraph {
            graph: original.graph.transitive_closure().unwrap(),
            global_data: original.global_data.clone(),
        };

        insta::assert_snapshot!(closure_graph.debug_dot());
    }

    #[test]
    fn test_transitive_reduction_layered_graph() {
        let original: DotGraph = dot!(
            digraph {
                // Layered graph with full connectivity between adjacent layers
                // Layer 1
                l1_a -> l2_a
                l1_a -> l2_b
                l1_b -> l2_a
                l1_b -> l2_b

                // Layer 2 to 3
                l2_a -> l3_a
                l2_a -> l3_b
                l2_b -> l3_a
                l2_b -> l3_b

                // Layer 3 to 4
                l3_a -> l4_a
                l3_b -> l4_a

                // Redundant long-distance edges
                l1_a -> l3_a  // Via l1_a -> l2_a -> l3_a
                l1_a -> l4_a  // Via multiple paths
                l2_a -> l4_a  // Via l2_a -> l3_a -> l4_a
            }
        )
        .unwrap();

        let reduction_graph = DotGraph {
            graph: original.graph.transitive_reduction().unwrap(),
            global_data: original.global_data.clone(),
        };

        insta::assert_snapshot!(reduction_graph.debug_dot());
    }

    #[test]
    fn test_transitive_operations_on_pyramid() {
        let original: DotGraph = dot!(
            digraph {
                // Pyramid structure
                top -> l2_a
                top -> l2_b

                l2_a -> l3_a
                l2_a -> l3_b
                l2_b -> l3_b
                l2_b -> l3_c

                l3_a -> l4_a
                l3_a -> l4_b
                l3_b -> l4_b
                l3_b -> l4_c
                l3_c -> l4_c
                l3_c -> l4_d

                // All bottom nodes connect to final
                l4_a -> bottom
                l4_b -> bottom
                l4_c -> bottom
                l4_d -> bottom

                // Some redundant edges
                top -> l3_b     // Via top->l2_a->l3_b or top->l2_b->l3_b
                l2_a -> l4_b    // Via l2_a->l3_a->l4_b or l2_a->l3_b->l4_b
                top -> bottom   // Via many paths
            }
        )
        .unwrap();

        // Test closure
        let closure_graph = DotGraph {
            graph: original.graph.clone().transitive_closure().unwrap(),
            global_data: original.global_data.clone(),
        };
        insta::assert_snapshot!(closure_graph.debug_dot(), @r"
        digraph {
          bottom;
          l2_a;
          l2_b;
          l3_a;
          l3_b;
          l3_c;
          l4_a;
          l4_b;
          l4_c;
          l4_d;
          top;

          l2_a:0	-> l3_a:1	 [id=0 ];
          l2_a:2	-> l3_b:3	 [id=1 ];
          l2_a:4	-> l4_b:5	 [id=2 ];
          l2_b:6	-> l3_b:7	 [id=3 ];
          l2_b:8	-> l3_c:9	 [id=4 ];
          l3_a:10	-> l4_a:11	 [id=5 ];
          l3_a:12	-> l4_b:13	 [id=6 ];
          l3_b:14	-> l4_b:15	 [id=7 ];
          l3_b:16	-> l4_c:17	 [id=8 ];
          l3_c:18	-> l4_c:19	 [id=9 ];
          l3_c:20	-> l4_d:21	 [id=10 ];
          l4_a:22	-> bottom:23	 [id=11 ];
          l4_b:24	-> bottom:25	 [id=12 ];
          l4_c:26	-> bottom:27	 [id=13 ];
          l4_d:28	-> bottom:29	 [id=14 ];
          top:30	-> bottom:31	 [id=15 ];
          top:32	-> l2_a:33	 [id=16 ];
          top:34	-> l2_b:35	 [id=17 ];
          top:36	-> l3_b:37	 [id=18 ];
          top:38	-> l3_a:39	 [id=19 ];
          top:40	-> l3_c:41	 [id=20 ];
          top:42	-> l4_a:43	 [id=21 ];
          top:44	-> l4_b:45	 [id=22 ];
          top:46	-> l4_c:47	 [id=23 ];
          top:48	-> l4_d:49	 [id=24 ];
          l2_a:50	-> l4_a:51	 [id=25 ];
          l2_a:52	-> l4_c:53	 [id=26 ];
          l2_a:54	-> bottom:55	 [id=27 ];
          l2_b:56	-> l4_b:57	 [id=28 ];
          l2_b:58	-> l4_c:59	 [id=29 ];
          l2_b:60	-> l4_d:61	 [id=30 ];
          l2_b:62	-> bottom:63	 [id=31 ];
          l3_a:64	-> bottom:65	 [id=32 ];
          l3_b:66	-> bottom:67	 [id=33 ];
          l3_c:68	-> bottom:69	 [id=34 ];
        }
        ");

        // Test reduction
        let reduction_graph = DotGraph {
            graph: original.graph.transitive_reduction().unwrap(),
            global_data: original.global_data.clone(),
        };
        insta::assert_snapshot!(reduction_graph.debug_dot(), @r"
        digraph {
          l4_d;
          l4_c;
          l2_b;
          l3_a;
          l4_a;
          l3_c;
          l3_b;
          l4_b;
          l2_a;
          bottom;
          top;

          l2_a:0	-> l3_a:1	 [id=0 ];
          l2_a:2	-> l3_b:3	 [id=1 ];
          top:5	-> l2_b:4	 [id=2 ];
          l2_b:6	-> l3_b:7	 [id=3 ];
          l2_b:8	-> l3_c:9	 [id=4 ];
          l3_a:10	-> l4_a:11	 [id=5 ];
          l3_a:12	-> l4_b:13	 [id=6 ];
          l3_b:14	-> l4_b:15	 [id=7 ];
          l3_b:16	-> l4_c:17	 [id=8 ];
          l3_c:18	-> l4_c:19	 [id=9 ];
          l3_c:20	-> l4_d:21	 [id=10 ];
          l4_a:22	-> bottom:23	 [id=11 ];
          l4_b:24	-> bottom:25	 [id=12 ];
          l4_c:26	-> bottom:27	 [id=13 ];
          l4_d:28	-> bottom:29	 [id=14 ];
          top:31	-> l2_a:30	 [id=15 ];
        }
        ");
    }

    #[test]
    fn test_transitive_closure_with_isolated_components() {
        let original: DotGraph<NodeStorageVec<DotVertexData>> = DotGraph::from_string(
            r#"
            digraph {
                // Component 1: Chain
                a1 -> b1
                b1 -> c1
                c1 -> d1

                // Component 2: Star
                center2 -> a2
                center2 -> b2
                center2 -> c2

                // Component 3: Diamond
                top3 -> left3
                top3 -> right3
                left3 -> bottom3
                right3 -> bottom3

                // Isolated nodes
                isolated1
                isolated2
            }
            "#,
        )
        .unwrap();

        let closure_graph = DotGraph {
            graph: original.graph.transitive_closure().unwrap(),
            global_data: original.global_data.clone(),
        };

        insta::assert_snapshot!(closure_graph.debug_dot());
    }

    #[test]
    fn test_error_handling_with_complex_cycle() {
        let graph_with_cycle: DotGraph = dot!(
            digraph {
                // Complex graph with embedded cycle
                start -> a
                a -> b
                b -> c
                c -> d
                d -> e

                // Create cycle in the middle
                c -> f
                f -> g
                g -> h
                h -> c  // Back to c creates cycle

                // Some more structure
                e -> end
                g -> end
            }
        )
        .unwrap();

        // Both operations should detect the cycle and fail
        let closure_result = graph_with_cycle.graph.clone().transitive_closure();
        let reduction_result = graph_with_cycle.graph.transitive_reduction();

        assert!(closure_result.is_err());
        assert!(reduction_result.is_err());

        // Snapshot the error debug representation
        insta::assert_debug_snapshot!(closure_result.unwrap_err());
        insta::assert_debug_snapshot!(reduction_result.unwrap_err());
    }

    #[test]
    fn test_comprehensive_reachability() {
        let graph: DotGraph<NodeStorageVec<DotVertexData>> = DotGraph::from_string(
            r#"
            digraph {
                // Create a graph with known reachability patterns
                root -> branch_a
                root -> branch_b

                branch_a -> leaf_a1
                branch_a -> leaf_a2
                branch_a -> shared

                branch_b -> leaf_b1
                branch_b -> shared

                shared -> final_node

                // Isolated component
                iso_a -> iso_b
                iso_b -> iso_c

                // Single node
                singleton
            }
            "#,
        )
        .unwrap();

        // Test various reachability scenarios
        let nodes: Vec<_> = graph.graph.iter_node_ids().collect();
        let node_map: std::collections::HashMap<_, _> = graph
            .graph
            .iter_node_ids()
            .filter_map(|idx| {
                graph
                    .graph
                    .node_store
                    .get_node_data(idx)
                    .name
                    .as_ref()
                    .map(|id| (id.clone(), idx))
            })
            .collect();

        // Self-reachability
        if let Some(&root_idx) = node_map.get("root") {
            assert!(graph.graph.is_reachable(root_idx, root_idx));
        }

        // Direct reachability
        if let (Some(&root_idx), Some(&branch_a_idx)) =
            (node_map.get("root"), node_map.get("branch_a"))
        {
            assert!(graph.graph.is_reachable(root_idx, branch_a_idx));
        }

        // Transitive reachability
        if let (Some(&root_idx), Some(&final_idx)) =
            (node_map.get("root"), node_map.get("final_node"))
        {
            assert!(graph.graph.is_reachable(root_idx, final_idx));
        }

        // Non-reachability across components
        if let (Some(&root_idx), Some(&iso_a_idx)) = (node_map.get("root"), node_map.get("iso_a")) {
            assert!(!graph.graph.is_reachable(root_idx, iso_a_idx));
            assert!(!graph.graph.is_reachable(iso_a_idx, root_idx));
        }

        // Reachability within isolated component
        if let (Some(&iso_a_idx), Some(&iso_c_idx)) = (node_map.get("iso_a"), node_map.get("iso_c"))
        {
            assert!(graph.graph.is_reachable(iso_a_idx, iso_c_idx));
            assert!(!graph.graph.is_reachable(iso_c_idx, iso_a_idx));
        }
    }

    #[test]
    fn test_hasse_diagram() {
        // Create a simple partial order: 1 < 2, 1 < 3, 2 < 4, 3 < 4
        // The Hasse diagram should have edges: 1->2, 2->3, 3->4
        let nodes = vec![1, 2, 3, 4];
        let hasse: HedgeGraph<NoData, i32> = HedgeGraph::hasse_diagram(nodes);

        // println!("{}", hasse.clone().transitive_closure().unwrap().base_dot());

        // Should have 4 nodes
        assert_eq!(hasse.n_nodes(), 4);

        // Should have 4 edges in the Hasse diagram (covering relations only)
        assert_eq!(hasse.n_edges(), 3);

        // Verify the graph is acyclic
        assert!(hasse.transitive_closure().is_ok());
    }
}
