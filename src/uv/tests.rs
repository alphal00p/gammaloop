use std::time::Instant;

use ahash::HashMap;
use smartstring::SmartString;

use crate::{
    feyngen::diagram_generator::FeynGen,
    tests_from_pytest::{load_amplitude_output, load_generic_model},
    uv::{PoSet, UVGraph},
};

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

    // let all_cycles = uv_graph.0.read_tarjan();
    // assert_eq!(all_cycles.len(), 1);

    // for cycle in all_cycles {
    // println!("{}", uv_graph.0.dot(&cycle));
    // }

    // insta::assert_ron_snapshot!("lbl_cycles", cycles);

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

    println!("{}", wood.dot(&uv_graph));

    let mut ufold = wood.unfold_impl(&uv_graph);
    ufold.compute(&uv_graph);

    println!("unfolded : {}", ufold.show_structure(&uv_graph).unwrap());
    println!("graph: {}", ufold.graphs());

    // let structure = wood.unfold(&uv_graph);

    // println!("{}", structure.show_structure(&wood, &uv_graph));
    // println!("{}", structure.n_elements());
}

#[test]
#[allow(unused)]
fn bugblatter_forest() {
    let model = load_generic_model("sm");
    println!("{}", model.vertex_rules[0].name);
    let mut symbolica_graph = symbolica::graph::Graph::new();

    let l1 = symbolica_graph.add_node((0, "V_137".into()));
    let l2 = symbolica_graph.add_node((0, "V_137".into()));
    let l3 = symbolica_graph.add_node((0, "V_137".into()));
    let l4 = symbolica_graph.add_node((0, "V_137".into()));
    let l5 = symbolica_graph.add_node((0, "V_137".into()));
    let l6 = symbolica_graph.add_node((0, "V_137".into()));

    symbolica_graph.add_edge(l1, l2, true, "t");
    symbolica_graph.add_edge(l2, l3, true, "t");
    symbolica_graph.add_edge(l3, l1, true, "t");

    symbolica_graph.add_edge(l4, l5, true, "t");
    symbolica_graph.add_edge(l5, l6, true, "t");
    symbolica_graph.add_edge(l6, l4, true, "t");

    symbolica_graph.add_edge(l1, l4, true, "g");
    symbolica_graph.add_edge(l2, l5, true, "g");
    symbolica_graph.add_edge(l3, l6, true, "g");

    let bare_graph = BareGraph::from_symbolica_graph(
        &model,
        "bugblatter".into(),
        &symbolica_graph,
        "str".into(),
        vec![],
        None,
    )
    .unwrap();

    let uv_graph = UVGraph::from_graph(&bare_graph);

    println!("{}", uv_graph.0.base_dot());

    let wood = uv_graph.wood();

    assert_eq!(20, wood.n_spinneys());
    println!("{}", wood.dot(&uv_graph));
    println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold_impl(&uv_graph);
    assert_eq!(152, ufold.n_terms());
    ufold.compute(&uv_graph);

    println!("unfolded : {}", ufold.show_structure(&uv_graph).unwrap());
    println!("graph: {}", ufold.graphs());

    // let structure = wood.unfold(&uv_graph);

    // println!("{}", structure.show_structure(&wood, &uv_graph));
    // println!("{}", structure.n_elements());
}

#[test]
#[allow(unused)]
fn kaapo_triplering() {
    let model = load_generic_model("scalars");
    let mut symbolica_graph = symbolica::graph::Graph::new();

    let l1 = symbolica_graph.add_node((0, "V_4_SCALAR_0000".into()));
    let l2 = symbolica_graph.add_node((0, "V_4_SCALAR_0000".into()));
    let l3 = symbolica_graph.add_node((0, "V_4_SCALAR_0000".into()));

    symbolica_graph.add_edge(l1, l2, true, "scalar_0");
    symbolica_graph.add_edge(l2, l3, true, "scalar_0");
    symbolica_graph.add_edge(l3, l1, true, "scalar_0");
    symbolica_graph.add_edge(l1, l2, true, "scalar_0");
    symbolica_graph.add_edge(l2, l3, true, "scalar_0");
    symbolica_graph.add_edge(l3, l1, true, "scalar_0");

    let bare_graph = BareGraph::from_symbolica_graph(
        &model,
        "threeringscalar".into(),
        &symbolica_graph,
        "1".into(),
        vec![],
        None,
    )
    .unwrap();

    let uv_graph = UVGraph::from_graph(&bare_graph);

    // println!("{}", uv_graph.0.base_dot());

    let wood = uv_graph.wood();
    assert_eq!(26, wood.n_spinneys());

    // println!("{}", wood.dot(&uv_graph));
    // println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold_impl(&uv_graph);
    ufold.compute(&uv_graph);

    println!("unfolded : {}", ufold.show_structure(&uv_graph).unwrap());
    println!("graph: {}", ufold.graphs());

    assert_eq!(242, ufold.n_terms());
}

#[test]
#[allow(unused)]
fn kaapo_quintic_scalar() {
    let model = load_generic_model("scalars");
    let mut symbolica_graph = symbolica::graph::Graph::new();

    let l1 = symbolica_graph.add_node((0, "V_3_SCALAR_000".into()));
    let l2 = symbolica_graph.add_node((0, "V_4_SCALAR_0000".into()));
    let l3 = symbolica_graph.add_node((0, "V_5_SCALAR_00000".into()));

    symbolica_graph.add_edge(l1, l2, true, "scalar_0");
    symbolica_graph.add_edge(l2, l3, true, "scalar_0");
    symbolica_graph.add_edge(l3, l1, true, "scalar_0");
    symbolica_graph.add_edge(l3, l2, true, "scalar_0");
    symbolica_graph.add_edge(l2, l3, true, "scalar_0");
    symbolica_graph.add_edge(l3, l1, true, "scalar_0");

    let bare_graph = BareGraph::from_symbolica_graph(
        &model,
        "threeringscalar".into(),
        &symbolica_graph,
        "1".into(),
        vec![],
        None,
    )
    .unwrap();

    let uv_graph = UVGraph::from_graph(&bare_graph);

    println!("{}", uv_graph.0.base_dot());

    let wood = uv_graph.wood();

    assert_eq!(25, wood.n_spinneys());

    // println!("{}", wood.dot(&uv_graph));
    // println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold_impl(&uv_graph);
    ufold.compute(&uv_graph);

    println!("unfolded : {}", ufold.show_structure(&uv_graph).unwrap());
    println!("graph: {}", ufold.graphs());

    assert_eq!(248, ufold.n_terms());
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
