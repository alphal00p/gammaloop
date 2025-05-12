use std::{sync::Arc, time::Instant};

use ahash::HashMap;
use smartstring::SmartString;
use spenso::{
    network::parsing::ShadowedStructure,
    structure::{
        abstract_index::AbstractIndex,
        dimension::Dimension,
        representation::{Minkowski, RepName},
        NamedStructure, ToSymbolic,
    },
};
use symbolica::evaluate::{FunctionMap, OptimizationSettings};
use symbolica_community::physics::algebraic_simplification::metric::MetricSimplifier;

use crate::{
    feyngen::diagram_generator::{EdgeColor, FeynGen, NodeColorWithVertexRule},
    model::{ArcVertexRule, ColorStructure, VertexRule},
    tests_from_pytest::{load_amplitude_output, load_generic_model},
    uv::{PoSet, UVGraph},
};

pub fn spenso_lor(
    tag: i32,
    ind: impl Into<AbstractIndex>,
    dim: impl Into<Dimension>,
) -> ShadowedStructure {
    let mink = Minkowski {}.new_slot(dim, ind);
    NamedStructure::from_iter([mink], GS.emr_mom, Some(vec![Atom::new_num(tag)]))
}

pub fn spenso_lor_atom(tag: i32, ind: impl Into<AbstractIndex>, dim: impl Into<Dimension>) -> Atom {
    let mink = Minkowski {}.new_slot(dim, ind);
    // spenso_lor(tag, ind, dim).to_symbolic().unwrap()
    vec![mink].to_symbolic_with(GS.emr_mom, &[Atom::new_num(tag)])
}

#[test]
fn disconnect_forest_scalar() {
    let scalar_node = UVNode {
        dod: 0,
        num: Atom::new_num(1),
        color: None,
    };

    fn scalar_edge(eid: i32) -> UVEdge {
        let m2 = parse!("m^2").unwrap();
        UVEdge {
            og_edge: 1, // not needed
            dod: -2,
            num: Atom::new_num(1),
            den: spenso_lor_atom(eid, 1, GS.dim).npow(2).to_dots() + m2,
        }
    }

    let mut builder = HedgeGraphBuilder::new();

    let node1 = builder.add_node(scalar_node.clone());

    let node2 = builder.add_node(scalar_node.clone());

    let node3 = builder.add_node(scalar_node.clone());

    builder.add_edge(node1, node2, scalar_edge(0), false);

    builder.add_edge(node1, node3, scalar_edge(1), false);

    builder.add_edge(node2, node3, scalar_edge(2), false);
    builder.add_edge(node2, node3, scalar_edge(3), false);

    builder.add_external_edge(node1, scalar_edge(4), false, Flow::Sink);

    builder.add_external_edge(node2, scalar_edge(5), false, Flow::Source);

    let hedge_graph = builder.build();
    let mut uv_graph = UVGraph::from_hedge(hedge_graph);
    println!(
        "{}",
        uv_graph.dot_impl(
            &uv_graph.full_filter(),
            "",
            &|a| Some(a.den.to_string()),
            &|n| Some(n.num.to_string())
        )
    );

    let wood = uv_graph.wood();

    println!("{}", wood.dot(&uv_graph));
    println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold_impl(&uv_graph);
    // assert_eq!(152, ufold.n_terms());
    ufold.compute(&uv_graph);

    println!("unfolded : {}", ufold.show_structure(&uv_graph).unwrap());
    println!("graph: {}", ufold.graphs());

    let result = ufold.expr(&uv_graph).unwrap().0;

    println!("{}", ufold.structure_and_res(&uv_graph));
    println!("{:>}", result);

    let exp = result
        .replace(parse!("symbolica_community::dot(k_(x_),l_(y_))").unwrap())
        .with(parse!("k(x_,0)*k(y_,0)-k(x_,1)*k(y_,1)-k(x_,2)*k(y_,2)-k(x_,3)*k(y_,3)").unwrap());

    println!("{exp}");
    let mut fnmap = FunctionMap::new();

    fnmap.add_constant(Atom::new_var(symbol!("m")), (1.).into());
    let ev = exp
        .evaluator(
            &fnmap,
            &[
                parse!("k(1, 0)").unwrap(),
                parse!("k(1, 1)").unwrap(),
                parse!("k(1, 2)").unwrap(),
                parse!("k(1, 3)").unwrap(),
                parse!("k(3, 0)").unwrap(),
                parse!("k(3, 1)").unwrap(),
                parse!("k(3, 2)").unwrap(),
                parse!("k(3, 3)").unwrap(),
                parse!("k(5, 0)").unwrap(),
                parse!("k(5, 1)").unwrap(),
                parse!("k(5, 2)").unwrap(),
                parse!("k(5, 3)").unwrap(),
            ],
            OptimizationSettings::default(),
        )
        .unwrap();

    let mut ev2 = ev.map_coeff(&|x| x.to_f64());

    for t in (0..100_000).step_by(1000) {
        let r = ev2.evaluate_single(&[
            3.,
            43.,
            5.,
            6.5,
            t as f64 + 1.,
            t as f64 * 2. + 2.,
            t as f64 + 3.,
            t as f64 + 4.,
            7.,
            4.2,
            1.,
            2.4,
        ]);
        println!("{} {}", r, r.abs().log10());
    }
}

#[test]
fn manual() {
    let mut builder = HedgeGraphBuilder::new();

    let node1 = builder.add_node(UVNode {
        dod: 0,
        num: spenso_lor_atom(1, 2, GS.dim),
        color: None,
    });

    let node2 = builder.add_node(UVNode {
        dod: 1,
        num: parse!("anything").unwrap(),
        color: None,
    });

    let edge1 = UVEdge {
        og_edge: 1, // not needed
        dod: 1,
        num: spenso_lor_atom(1, 2, GS.dim),
        den: spenso_lor_atom(1, -2, GS.dim).npow(2).to_dots(),
    };

    builder.add_edge(node1, node2, edge1, false);

    let edge2 = UVEdge {
        og_edge: 1, // not needed
        dod: 1,
        num: spenso_lor_atom(2, 2, GS.dim),
        den: spenso_lor_atom(2, -1, GS.dim).npow(2).to_dots(),
    };

    builder.add_edge(node2, node1, edge2, false);

    let hedge_graph = builder.build();
    let cut_edges = hedge_graph.cycle_basis().1.tree_subgraph;

    let uv_graph = UVGraph {
        hedge_graph,
        cut_edges,
        lmb_replacement: vec![],
    };
    println!(
        "{}",
        uv_graph.dot_impl(
            &uv_graph.full_filter(),
            "",
            &|a| Some(a.den.to_string()),
            &|n| Some(n.num.to_string())
        )
    );

    let wood = uv_graph.wood();

    println!("{}", wood.dot(&uv_graph));
    println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold_impl(&uv_graph);
    // assert_eq!(152, ufold.n_terms());
    ufold.compute(&uv_graph);

    println!("unfolded : {}", ufold.show_structure(&uv_graph).unwrap());
    println!("graph: {}", ufold.graphs());

    println!("{}", ufold.expr(&uv_graph).unwrap());
}

#[test]
#[allow(unused)]
fn easy() {
    let model = load_generic_model("sm");
    let mut symbolica_graph = symbolica::graph::Graph::new();

    let dummy_external_vertex_rule = ArcVertexRule(Arc::new(VertexRule {
        name: "external".into(),
        couplings: vec![],
        lorentz_structures: vec![],
        particles: vec![],
        color_structures: ColorStructure::new(vec![]),
    }));

    let incoming = symbolica_graph.add_node(NodeColorWithVertexRule {
        external_tag: 1,
        vertex_rule: dummy_external_vertex_rule.clone(),
    });

    let outgoing = symbolica_graph.add_node(NodeColorWithVertexRule {
        external_tag: 2,
        vertex_rule: dummy_external_vertex_rule.clone(),
    });

    let tth = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_141"),
    };

    let t = EdgeColor::from_particle(model.get_particle("t"));
    let h = EdgeColor::from_particle(model.get_particle("H"));

    let l1 = symbolica_graph.add_node(tth.clone());
    let l2 = symbolica_graph.add_node(tth.clone());

    let e1 = symbolica_graph.add_edge(incoming, l1, false, h).unwrap();
    let e2 = symbolica_graph.add_edge(l2, outgoing, false, h).unwrap();
    symbolica_graph.add_edge(l1, l2, true, t);
    symbolica_graph.add_edge(l2, l1, true, t);

    let bare_graph = BareGraph::from_symbolica_graph(
        &model,
        "1l_prop".into(),
        &symbolica_graph,
        Atom::new_num(1),
        vec![((Some(1), Some(2)))],
        None,
    )
    .unwrap();

    let uv_graph = UVGraph::from_graph(&bare_graph);

    println!("{}", uv_graph.base_dot());

    let wood = uv_graph.wood();

    println!("{}", wood.dot(&uv_graph));
    println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold_impl(&uv_graph);
    // assert_eq!(152, ufold.n_terms());
    ufold.compute(&uv_graph);

    println!("unfolded : {}", ufold.show_structure(&uv_graph).unwrap());
    println!("graph: {}", ufold.graphs());

    // let structure = wood.unfold(&uv_graph);

    // println!("{}", structure.show_structure(&wood, &uv_graph));
    // println!("{}", structure.n_elements());
}

#[test]
#[allow(unused)]
fn lbl() {
    let (model, amplitude, _) = load_amplitude_output("TEST_AMPLITUDE_lbl_box/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    println!("{}", graph.bare_graph.dot());

    // graph.generate_uv();

    graph.generate_loop_momentum_bases();

    let uv_graph = UVGraph::from_graph(&graph.bare_graph);

    insta::assert_snapshot!("lbl_dot", uv_graph.base_dot());

    let lmb = graph.bare_graph.loop_momentum_basis.clone();

    // let cycles = uv_graph.cycle_basis_from_lmb(&lmb);

    // let all_cycles = uv_graph.read_tarjan();
    // assert_eq!(all_cycles.len(), 1);

    // for cycle in all_cycles {
    // println!("{}", uv_graph.dot(&cycle));
    // }

    // insta::assert_ron_snapshot!("lbl_cycles", cycles);

    let uv_graph = UVGraph::from_graph(&graph.bare_graph);

    println!("tbt_dot{}", uv_graph.base_dot());

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

    println!("tbt_dot{}", uv_graph.base_dot());

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
    // println!("{}", env!("CARGO_CRATE_NAME"));
    let model = load_generic_model("sm");
    println!("{}", model.vertex_rules[0].0.name);
    let mut symbolica_graph = symbolica::graph::Graph::new();

    let ttg = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_137"),
    };

    let t = EdgeColor::from_particle(model.get_particle("t"));
    let g = EdgeColor::from_particle(model.get_particle("g"));

    let l1 = symbolica_graph.add_node(ttg.clone());
    let l2 = symbolica_graph.add_node(ttg.clone());
    let l3 = symbolica_graph.add_node(ttg.clone());
    let l4 = symbolica_graph.add_node(ttg.clone());
    let l5 = symbolica_graph.add_node(ttg.clone());
    let l6 = symbolica_graph.add_node(ttg.clone());

    symbolica_graph.add_edge(l1, l2, true, t);
    symbolica_graph.add_edge(l2, l3, true, t);
    symbolica_graph.add_edge(l3, l1, true, t);

    symbolica_graph.add_edge(l4, l5, true, t);
    symbolica_graph.add_edge(l5, l6, true, t);
    symbolica_graph.add_edge(l6, l4, true, t);

    symbolica_graph.add_edge(l1, l4, true, g);
    symbolica_graph.add_edge(l2, l5, true, g);
    symbolica_graph.add_edge(l3, l6, true, g);

    let bare_graph = BareGraph::from_symbolica_graph(
        &model,
        "bugblatter".into(),
        &symbolica_graph,
        Atom::new_num(1),
        vec![],
        None,
    )
    .unwrap();

    let uv_graph = UVGraph::from_graph(&bare_graph);

    println!("{}", uv_graph.base_dot());

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

    let four_scalar = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_4_SCALAR_0000"),
    };

    let scalar = EdgeColor::from_particle(model.get_particle("scalar_0"));

    let l1 = symbolica_graph.add_node(four_scalar.clone());
    let l2 = symbolica_graph.add_node(four_scalar.clone());
    let l3 = symbolica_graph.add_node(four_scalar.clone());

    symbolica_graph.add_edge(l1, l2, true, scalar.clone());
    symbolica_graph.add_edge(l2, l3, true, scalar.clone());
    symbolica_graph.add_edge(l3, l1, true, scalar.clone());
    symbolica_graph.add_edge(l1, l2, true, scalar.clone());
    symbolica_graph.add_edge(l2, l3, true, scalar.clone());
    symbolica_graph.add_edge(l3, l1, true, scalar.clone());

    let bare_graph = BareGraph::from_symbolica_graph(
        &model,
        "threeringscalar".into(),
        &symbolica_graph,
        Atom::new_num(1),
        vec![],
        None,
    )
    .unwrap();

    let uv_graph = UVGraph::from_graph(&bare_graph);

    // println!("{}", uv_graph.base_dot());

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

    let three = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_3_SCALAR_000"),
    };
    let four = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_4_SCALAR_0000"),
    };
    let five = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_5_SCALAR_00000"),
    };

    let scalar = EdgeColor::from_particle(model.get_particle("scalar_0"));

    let l1 = symbolica_graph.add_node(three);
    let l2 = symbolica_graph.add_node(four);
    let l3 = symbolica_graph.add_node(five);

    symbolica_graph.add_edge(l1, l2, true, scalar.clone());
    symbolica_graph.add_edge(l2, l3, true, scalar.clone());
    symbolica_graph.add_edge(l3, l1, true, scalar.clone());
    symbolica_graph.add_edge(l3, l2, true, scalar.clone());
    symbolica_graph.add_edge(l2, l3, true, scalar.clone());
    symbolica_graph.add_edge(l3, l1, true, scalar.clone());

    let bare_graph = BareGraph::from_symbolica_graph(
        &model,
        "threeringscalar".into(),
        &symbolica_graph,
        Atom::new_num(1),
        vec![],
        None,
    )
    .unwrap();

    let uv_graph = UVGraph::from_graph(&bare_graph);

    println!("{}", uv_graph.base_dot());

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
