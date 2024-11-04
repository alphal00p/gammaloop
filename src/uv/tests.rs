use std::{path::Path, time::Instant};

use crate::{
    graph::Graph,
    model::Model,
    numerator::UnInit,
    tests_from_pytest::{load_amplitude_output, load_generic_model},
    uv::{PoSet, UVGraph},
};

// fn self_energy_chain(num: usize) -> (Model, Graph<UnInit>) {
//     let model = load_generic_model("scalars");
// }

#[test]
#[allow(unused)]
fn lbl() {
    let (model, amplitude, _) = load_amplitude_output("TEST_AMPLITUDE_lbl_box/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    println!("{}", graph.bare_graph.dot());

    // graph.generate_uv();

    graph.generate_loop_momentum_bases();

    let uv_graph = UVGraph::from_graph(&graph.bare_graph);

    insta::assert_snapshot!("lbl_dot", uv_graph.0.base_dot());

    let lmb = graph.bare_graph.loop_momentum_basis.clone();

    let cycles = uv_graph.cycle_basis_from_lmb(&lmb);

    let all_cycles = uv_graph.0.read_tarjan();
    assert_eq!(all_cycles.len(), 1);

    for cycle in all_cycles {
        println!("{}", uv_graph.0.dot(&cycle));
    }

    insta::assert_ron_snapshot!("lbl_cycles", cycles);

    let uv_graph = UVGraph::from_graph(&graph.bare_graph);

    println!("tbt_dot{}", uv_graph.0.base_dot());

    let wood = uv_graph.wood();

    let structure = wood.unfold(&uv_graph);

    println!("{}", structure.show_structure(&wood, &uv_graph));
    println!("{}", structure.n_elements());
}

#[test]
#[allow(unused)]
fn tbt() {
    let (model, amplitude, _) =
        load_amplitude_output("TEST_AMPLITUDE_triangle_box_triangle_phys/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    // graph.generate_uv();

    graph.generate_loop_momentum_bases();

    let uv_graph = UVGraph::from_graph(&graph.bare_graph);

    println!("tbt_dot{}", uv_graph.0.base_dot());

    let wood = uv_graph.wood();

    let structure = wood.unfold(&uv_graph);

    println!("{}", structure.show_structure(&wood, &uv_graph));
    println!("{}", structure.n_elements());
}

use super::*;

use std::cmp::Ordering;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct TestNode(&'static str);

#[allow(clippy::non_canonical_partial_ord_impl)]
// Implement PartialOrd and Ord for TestNode to define the partial order
impl PartialOrd for TestNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self.0, other.0) {
            ("A", "A") => Some(Ordering::Equal),
            ("B", "B") => Some(Ordering::Equal),
            ("C", "C") => Some(Ordering::Equal),
            ("D", "D") => Some(Ordering::Equal),
            ("A", _) => Some(Ordering::Less), // A < B, A < C, A < D
            (_, "A") => Some(Ordering::Greater), // B > A, C > A, D > A
            ("B", "D") => Some(Ordering::Less), // B < D
            ("D", "B") => Some(Ordering::Greater), // D > B
            ("C", "D") => Some(Ordering::Less), // C < D
            ("D", "C") => Some(Ordering::Greater), // D > C
            _ => None,                        // Other pairs are incomparable
        }
    }
}

impl Ord for TestNode {
    fn cmp(&self, other: &Self) -> Ordering {
        // Since our partial order can be undefined for some pairs (None), we unwrap safely here
        // because we know our test cases only involve comparable nodes.
        self.partial_cmp(other).unwrap()
    }
}

#[test]
fn test_poset_topological_order() {
    // Define nodes
    let nodes = vec![TestNode("A"), TestNode("B"), TestNode("C"), TestNode("D")];

    // Build PoSet from iterator
    let poset = PoSet::from_iter(nodes.clone());

    // Expected topological order
    // Since A < B, A < C, B < D, C < D, the expected order is ["A", "B", "C", "D"]

    let expected_order = vec!["A", "B", "C", "D"];

    // Check that the nodes are in the expected order
    for (node, &expected_label) in poset.nodes.iter().zip(&expected_order) {
        assert_eq!(node.0, expected_label);
    }
}

#[test]
fn test_coverset_topological_order() {
    // Define nodes
    let nodes = vec![TestNode("A"), TestNode("B"), TestNode("C"), TestNode("D")];

    // Build PoSet from iterator
    let poset = PoSet::from_iter(nodes.clone());

    // Convert PoSet to CoverSet
    let coverset = poset.to_cover_set();

    // Expected topological order is ["A", "B", "C", "D"]
    let expected_order = vec!["A", "B", "C", "D"];

    // Check that the nodes are in the expected order
    for (node, &expected_label) in coverset.nodes.iter().zip(&expected_order) {
        assert_eq!(node.0, expected_label);
    }

    // Additionally, we can check that the covers are correct
    // For example, in the coverset, node A should cover nodes B and C
    // Node B and C should cover D

    // Get the indices of the nodes
    let a_index = coverset.nodes.iter().position(|n| n.0 == "A").unwrap();
    let b_index = coverset.nodes.iter().position(|n| n.0 == "B").unwrap();
    let c_index = coverset.nodes.iter().position(|n| n.0 == "C").unwrap();
    let d_index = coverset.nodes.iter().position(|n| n.0 == "D").unwrap();

    // Check that A covers B and C
    assert!(coverset.covers(a_index, b_index));
    assert!(coverset.covers(a_index, c_index));

    // Check that B and C cover D
    assert!(coverset.covers(b_index, d_index));
    assert!(coverset.covers(c_index, d_index));

    // Check that A does not directly cover D (since there is a node between them)
    assert!(!coverset.covers(a_index, d_index));
}

// #[test]
// fn test_triangle_graph_loop_count() {
//     // Create a triangle graph
//     let mut graph = BareGraph::new();
//     let v0 = graph.add_vertex(Vertex::new("v0"));
//     let v1 = graph.add_vertex(Vertex::new("v1"));
//     let v2 = graph.add_vertex(Vertex::new("v2"));

//     graph.add_edge(Edge::new("e0", v0, v1, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e1", v1, v2, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e2", v2, v0, EdgeType::Virtual));

//     // Generate UVGraph
//     let uv_graph = UVGraph::from_graph(&graph);

//     // Create SubGraph representing the full graph
//     let mut filter = bitvec![usize, Lsb0; 1; uv_graph.0.n_hedges()];
//     let subgraph = SubGraph::from(filter);

//     // Compute loop count
//     let loop_count = uv_graph.n_loops(&subgraph);
//     let basis_length = uv_graph.0.cycle_basis().len();

//     // Expected loop count is 1
//     assert_eq!(loop_count, 1);
//     assert_eq!(basis_length, 1);
// }

// #[test]
// fn test_square_graph_loop_count() {
//     // Create a square graph
//     let mut graph = BareGraph::new();
//     let v0 = graph.add_vertex(Vertex::new("v0"));
//     let v1 = graph.add_vertex(Vertex::new("v1"));
//     let v2 = graph.add_vertex(Vertex::new("v2"));
//     let v3 = graph.add_vertex(Vertex::new("v3"));

//     graph.add_edge(Edge::new("e0", v0, v1, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e1", v1, v2, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e2", v2, v3, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e3", v3, v0, EdgeType::Virtual));

//     // Generate UVGraph
//     let uv_graph = UVGraph::from_graph(&graph);

//     // Create SubGraph representing the full graph
//     let mut filter = bitvec![usize, Lsb0; 1; uv_graph.0.n_hedges()];
//     let subgraph = SubGraph::from(filter);

//     // Compute loop count
//     let loop_count = uv_graph.n_loops(&subgraph);
//     let basis_length = uv_graph.0.cycle_basis().len();

//     // Expected loop count is 1
//     assert_eq!(loop_count, 1);
//     assert_eq!(basis_length, 1);
// }

// #[test]
// fn test_figure_eight_graph_loop_count() {
//     // Create a figure-eight graph
//     let mut graph = BareGraph::new();
//     let v0 = graph.add_vertex(Vertex::new("v0"));
//     let v1 = graph.add_vertex(Vertex::new("v1"));
//     let v2 = graph.add_vertex(Vertex::new("v2"));

//     // First loop
//     graph.add_edge(Edge::new("e0", v0, v1, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e1", v1, v0, EdgeType::Virtual));

//     // Second loop
//     graph.add_edge(Edge::new("e2", v0, v2, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e3", v2, v0, EdgeType::Virtual));

//     // Generate UVGraph
//     let uv_graph = UVGraph::from_graph(&graph);

//     // Create SubGraph representing the full graph
//     let mut filter = bitvec![usize, Lsb0; 1; uv_graph.0.n_hedges()];
//     let subgraph = SubGraph::from(filter);

//     // Compute loop count
//     let loop_count = uv_graph.n_loops(&subgraph);
//     let basis_length = uv_graph.0.cycle_basis().len();

//     // Expected loop count is 2
//     assert_eq!(loop_count, 2);
//     assert_eq!(basis_length, 2);
// }

// #[test]
// fn test_complex_graph_loop_count() {
//     // Create a more complex graph
//     let mut graph = BareGraph::new();
//     let v0 = graph.add_vertex(Vertex::new("v0"));
//     let v1 = graph.add_vertex(Vertex::new("v1"));
//     let v2 = graph.add_vertex(Vertex::new("v2"));
//     let v3 = graph.add_vertex(Vertex::new("v3"));
//     let v4 = graph.add_vertex(Vertex::new("v4"));

//     graph.add_edge(Edge::new("e0", v0, v1, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e1", v1, v2, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e2", v2, v3, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e3", v3, v0, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e4", v1, v3, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e5", v3, v4, EdgeType::Virtual));
//     graph.add_edge(Edge::new("e6", v4, v1, EdgeType::Virtual));

//     // Generate UVGraph
//     let uv_graph = UVGraph::from_graph(&graph);

//     // Create SubGraph representing the full graph
//     let mut filter = bitvec![usize, Lsb0; 1; uv_graph.0.n_hedges()];
//     let subgraph = SubGraph::from(filter);

//     // Compute loop count
//     let loop_count = uv_graph.n_loops(&subgraph);
//     let basis_length = uv_graph.0.cycle_basis().len();

//     // Expected loop count is 3
//     // Calculated as E - N + 1 = 7 edges - 5 vertices + 1 = 3 loops
//     assert_eq!(loop_count, 3);
//     assert_eq!(basis_length, 3);
// }
