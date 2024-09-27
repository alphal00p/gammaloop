use std::time::Instant;

use itertools::Itertools;
use symbolica::printer::PrintOptions;

use crate::{
    tests_from_pytest::load_amplitude_output,
    uv::{PoSet, UVGraph},
};

use super::{BitFilter, Involution, NestingGraph, NestingGraphBuilder};

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

    // for cycle in cycles {
    //     println!("{:?}", cycle);
    //     println!("{}", uv_graph.0.dot(&cycle));
    // }

    // if let AtomView::Mul(mul) = numerator.as_view() {
    //     let net = SymbolicTensor::mul_to_tracking_network(mul).unwrap();

    //     println!("{}", net.dot());
    //     for (_, t) in net.graph.nodes {
    //         // println!("{}", t.structure());
    //     }
    // }

    for e in uv_graph
        .0
        .iter_internal_edge_data(&uv_graph.0.full_node().internal_graph)
    {
        println!("{}", e.num);
    }

    for (_, n) in uv_graph
        .0
        .iter_node_data(&uv_graph.0.full_node().internal_graph)
    {
        println!("{}", n.num);
    }

    println!("{}", uv_graph.dod(&uv_graph.0.full_node().internal_graph));
    println!("{}", uv_graph.wood().dot(&uv_graph));

    println!("{}", uv_graph.wood().unfold(&uv_graph));
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

    // let lmb = graph.bare_graph.loop_momentum_basis.clone();

    // let cycles = uv_graph.cycle_basis_from_lmb(&lmb);

    // let all_cycles = uv_graph.0.read_tarjan();
    // // assert_eq!(all_cycles.len(), 1);

    // for cycle in all_cycles {
    //     println!("{}", uv_graph.0.dot(&cycle));
    // }

    // // insta::assert_ron_snapshot!("lbl_cycles", cycles);

    // // for cycle in cycles {
    // //     println!("{:?}", cycle);
    // //     println!("{}", uv_graph.0.dot(&cycle));
    // // }

    // // if let AtomView::Mul(mul) = numerator.as_view() {
    // //     let net = SymbolicTensor::mul_to_tracking_network(mul).unwrap();

    // //     println!("{}", net.dot());
    // //     for (_, t) in net.graph.nodes {
    // //         // println!("{}", t.structure());
    // //     }
    // // }

    // for e in uv_graph
    //     .0
    //     .iter_internal_edge_data(&uv_graph.0.full_node().internal_graph)
    // {
    //     println!("{}", e.num);
    // }

    // for (_, n) in uv_graph
    //     .0
    //     .iter_node_data(&uv_graph.0.full_node().internal_graph)
    // {
    //     println!("{}", n.num);
    // }

    // println!("{}", uv_graph.dod(&uv_graph.0.full_node().internal_graph));

    // // println!("{}", pset.dot_structure());

    // // println!("{}", pset.remove_transitive_edges().dot_structure());

    // // uv_graph.wood().dot_spinneys(&uv_graph);

    // let mut wood = uv_graph.wood();

    // println!("{}", wood.dot(&uv_graph));

    // let shift = wood.poset.shift();

    // println!(
    //     "{}",
    //     wood.poset.nodes[wood.poset.minimum().unwrap()].dot_id(shift)
    // );

    // for (i, p) in wood.poset.bfs_paths().enumerate() {
    //     println!("path {}", i);
    //     for n in p {
    //         print!("{}->", wood.poset.nodes[n].dot_id(shift));
    //     }
    //     println!();
    // }

    // println!("Inverse");

    // wood.poset.invert();

    // println!("{}", wood.dot(&uv_graph));
    // for (i, p) in wood.poset.bfs_paths().enumerate() {
    //     println!("path {}", i);
    //     for n in p {
    //         print!("{}->", wood.poset.nodes[n].dot_id(shift));
    //     }
    //     println!();
    // }

    let wood = uv_graph.wood();

    let structure = wood.unfold(&uv_graph);
    println!("{}", structure.show_structure(&wood, &uv_graph));
}

#[test]
fn random_involution() {
    let inv = Involution::<(), ()>::random(10, 1);

    insta::assert_ron_snapshot!(inv);
}

#[test]
fn nine_loop() {
    let rand_graph = NestingGraph::<(), ()>::random(19, 30, 1);

    insta::assert_snapshot!("nine_loop_dot", rand_graph.base_dot());
    insta::assert_ron_snapshot!("nine_loop", rand_graph);
    let cycles = rand_graph.cycle_basis();

    for i in rand_graph.full_node().internal_graph.filter.iter_ones() {
        assert_eq!(
            9,
            rand_graph
                .paton_cycle_basis(&rand_graph.full_filter().into(), i)
                .unwrap()
                .len()
        );
    }

    let cycle_dots: Vec<String> = cycles.iter().map(|c| rand_graph.dot(c)).collect();

    let all_spinneys_other = rand_graph.all_spinneys();
    let all_cycles = rand_graph.all_cycles();

    assert_eq!(3361, all_spinneys_other.len());

    assert_eq!(all_cycles.len(), 109);

    insta::assert_toml_snapshot!("nine_loop_cycle_dots", cycle_dots);
}

#[test]
fn threeloop() {
    let mut builder: NestingGraphBuilder<(), ()> = NestingGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());

    builder.add_edge(a, b, ());
    builder.add_edge(b, a, ());
    builder.add_edge(a, b, ());

    builder.add_edge(b, c, ());
    builder.add_edge(c, d, ());
    builder.add_edge(d, a, ());

    let graph = builder.build();

    insta::assert_snapshot!("three_loop_dot", graph.base_dot());
    insta::assert_ron_snapshot!("three_loop", graph);

    for i in 0..graph.n_hedges() {
        assert_eq!(
            3,
            graph
                .paton_cycle_basis(&graph.full_filter().into(), i)
                .unwrap()
                .len()
        );
    }

    let cycles = graph.cycle_basis();

    assert_eq!(3, cycles.len());

    let all_cycles = graph.read_tarjan();

    assert_eq!(6, all_cycles.len());

    insta::assert_ron_snapshot!("three_loop_cycles", cycles);

    // let mut all_spinneys = graph.all_spinneys().into_iter().collect_vec();
    // all_spinneys.sort_by(|a, b| a.filter.cmp(&b.filter));

    // assert_eq!(graph.all_spinneys_alt().len(), all_spinneys.len());

    // // let posetSpinneys = all_spinneys.iter().collect::<PoSet<_>>();

    // // println!("{}", posetSpinneys.to_cover_set().dot());

    // insta::assert_toml_snapshot!(
    //     "three_loop_spinneys",
    //     all_spinneys
    //         .into_iter()
    //         .map(|s| graph.dot(&graph.nesting_node_from_subgraph(s)))
    //         .collect_vec()
    // );
}

#[test]
fn hairythreeloop() {
    let mut builder: NestingGraphBuilder<(), ()> = NestingGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());

    builder.add_edge(a, b, ());
    builder.add_edge(b, a, ());
    builder.add_edge(a, b, ());
    builder.add_external_edge(a, ());
    builder.add_external_edge(b, ());
    builder.add_external_edge(b, ());

    builder.add_edge(b, c, ());
    builder.add_edge(c, d, ());
    builder.add_edge(d, a, ());

    assert_eq!(builder.involution.len(), 15);
    let graph = builder.build();

    insta::assert_snapshot!("hairy_three_loop_dot", graph.base_dot());
    insta::assert_ron_snapshot!("hairy_three_loop", graph);
    insta::assert_snapshot!(
        "hairy_three_loop_dot_internal",
        graph.dot(&graph.full_node())
    );

    for i in graph.full_node().internal_graph.filter.iter_ones() {
        assert_eq!(
            3,
            graph
                .paton_cycle_basis(&graph.full_filter().into(), i)
                .unwrap()
                .len()
        );
    }

    let cycles = graph.cycle_basis();

    // for cycle in &cycles {
    //     println!("{}", graph.dot(cycle));
    // }

    insta::assert_ron_snapshot!("hairy_three_loop_cycles", cycles);
}

// #[test]
// fn poset() {
//     let a = (0..128)
//         .map(|i| BitFilter { data: i })
//         .collect::<PoSet<_>>();

//     println!("{}", a.to_cover_set().dot());
//     // println!("{}", a.dot());
// }

#[test]
fn cube() {
    let mut builder: NestingGraphBuilder<(), ()> = NestingGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());
    let e = builder.add_node(());
    let f = builder.add_node(());
    let g = builder.add_node(());
    let h = builder.add_node(());

    builder.add_edge(a, b, ());
    builder.add_edge(b, c, ());
    builder.add_edge(c, d, ());
    builder.add_edge(d, a, ());

    builder.add_edge(e, f, ());
    builder.add_edge(f, g, ());
    builder.add_edge(g, h, ());
    builder.add_edge(h, e, ());

    builder.add_edge(a, e, ());
    builder.add_edge(b, f, ());
    builder.add_edge(c, g, ());
    builder.add_edge(d, h, ());

    let graph = builder.build();

    insta::assert_snapshot!("cube_dot", graph.base_dot());

    // let mut all_spinneys = graph.all_spinneys().into_iter().collect_vec();
    // all_spinneys.sort_by(|a, b| a.filter.cmp(&b.filter));

    // assert_eq!(162, all_spinneys.len());
}

#[test]
fn alt_vs_pair() {
    for s in 0..100 {
        let rand_graph = NestingGraph::<(), ()>::random(10, 14, s);

        let before = Instant::now();
        let all_spinneys = rand_graph.all_spinneys();
        let after = before.elapsed();
        let before = Instant::now();
        let all_spinneys_alt = rand_graph.all_spinneys_alt();
        let after_alt = before.elapsed();
        // println!("{s} {after:?} {after_alt:?}");

        assert_eq!(
            all_spinneys.len(),
            all_spinneys_alt.len(),
            "{}",
            rand_graph.base_dot()
        );
    }
    // let rand_graph = NestingGraph::<(), ()>::random(6, 9, 8);

    // println!("{}", rand_graph.base_dot());

    // println!("loops {}", rand_graph.cycle_basis().len());

    // let all_spinneys_other = rand_graph.all_spinneys();

    // // println!("all spinneys read tarjan {}", all_spinneys.len());

    // println!("all spinneys {}", all_spinneys_other.len());
    // println!("all spinneys alt{}", rand_graph.all_spinneys_alt().len());
}

#[test]
#[should_panic]
fn read_tarjan_vs_cycle_space() {
    for s in 0..100 {
        let rand_graph = NestingGraph::<(), ()>::random(6, 9, s);

        let all_cycles = rand_graph.read_tarjan();
        let all_cycles_alt = rand_graph.all_cycles();

        assert_eq!(
            all_cycles.len(),
            all_cycles_alt.len(),
            "{} with seed {s}",
            rand_graph.base_dot()
        );
    }
}

// #[test]
// fn random_graph() {
//     let rand_graph = NestingGraph::<(), ()>::random(6, 9, 3);

//     println!(
//         "{} loop graph: \n {}",
//         rand_graph.cycle_basis().len(),
//         rand_graph.base_dot()
//     );

//     for c in rand_graph.all_cycles() {
//         println!(" {}", rand_graph.dot(&c));
//     }

//     for c in rand_graph.read_tarjan() {
//         println!("{}", rand_graph.dot(&c));
//     }
// }

#[test]
fn K33() {
    let mut builder: NestingGraphBuilder<(), ()> = NestingGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());
    let e = builder.add_node(());
    let f = builder.add_node(());

    builder.add_edge(a, d, ());
    builder.add_edge(a, e, ());
    builder.add_edge(a, f, ());

    builder.add_edge(b, d, ());
    builder.add_edge(b, e, ());
    builder.add_edge(b, f, ());

    builder.add_edge(c, d, ());
    builder.add_edge(c, e, ());
    builder.add_edge(c, f, ());

    let graph = builder.build();

    println!("{}", graph.dot(&graph.full_node()));

    for (s, v) in graph.all_spinneys() {
        println!("cyclotomatic_number: {}", graph.cyclotomatic_number(&s));
        println!(
            "paton_count_loops {}",
            graph
                .paton_count_loops(&s, s.filter.iter_ones().next().unwrap())
                .unwrap()
        );

        println!("paton_cycle_basislen {}", s.cycle_basis(&graph).len());
        println!("{}", graph.dot(&s.to_nesting_node(&graph)));
    }

    assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());
}

#[test]
fn petersen() {
    let mut builder: NestingGraphBuilder<(), ()> = NestingGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());
    let e = builder.add_node(());
    let f = builder.add_node(());
    let g = builder.add_node(());
    let h = builder.add_node(());
    let i = builder.add_node(());
    let j = builder.add_node(());

    builder.add_edge(a, b, ());
    builder.add_edge(a, f, ());
    builder.add_edge(a, e, ());

    builder.add_edge(b, c, ());
    builder.add_edge(b, g, ());

    builder.add_edge(c, d, ());
    builder.add_edge(c, h, ());

    builder.add_edge(d, e, ());
    builder.add_edge(d, i, ());

    builder.add_edge(e, j, ());

    builder.add_edge(f, h, ());
    builder.add_edge(f, i, ());

    builder.add_edge(g, i, ());
    builder.add_edge(g, j, ());

    builder.add_edge(h, j, ());

    let graph = builder.build();

    println!("{}", graph.base_dot());

    // assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());

    println!("loop count {}", graph.cycle_basis().len());
    println!("cycle count {}", graph.all_cycles().len());
    if let Some((s, v)) = graph
        .all_spinneys()
        .iter()
        .find(|(s, _)| graph.full_filter() == s.filter)
    {
        println!(
            "{}",
            graph.dot(&graph.nesting_node_from_subgraph(s.clone()))
        );
        for (ci, cj) in v {
            println!("{}", graph.dot(&ci));
            println!("{}", graph.dot(&cj.as_ref().unwrap()));
        }
    } else {
        println!("not found");
    }
}

#[test]
fn wagner_graph() {
    let mut builder: NestingGraphBuilder<(), ()> = NestingGraphBuilder::new();
    let n1 = builder.add_node(());
    let n2 = builder.add_node(());
    let n3 = builder.add_node(());
    let n4 = builder.add_node(());
    let n5 = builder.add_node(());
    let n6 = builder.add_node(());
    let n7 = builder.add_node(());
    let n8 = builder.add_node(());

    builder.add_edge(n1, n2, ());
    builder.add_edge(n1, n5, ());

    builder.add_edge(n2, n3, ());
    builder.add_edge(n2, n6, ());

    builder.add_edge(n3, n4, ());
    builder.add_edge(n3, n7, ());

    builder.add_edge(n4, n5, ());
    builder.add_edge(n4, n8, ());

    builder.add_edge(n5, n6, ());

    builder.add_edge(n6, n7, ());

    builder.add_edge(n7, n8, ());

    builder.add_edge(n8, n1, ());

    let graph = builder.build();

    println!("{}", graph.base_dot());

    // assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());

    println!("loop count {}", graph.cycle_basis().len());
    println!("cycle count {}", graph.all_cycles().len());
    if let Some((s, v)) = graph
        .all_spinneys()
        .iter()
        .find(|(s, _)| graph.full_filter() == s.filter)
    {
        println!(
            "{}",
            graph.dot(&graph.nesting_node_from_subgraph(s.clone()))
        );
        for (ci, cj) in v {
            println!("{}", graph.dot(&ci));
            println!("{}", graph.dot(&cj.as_ref().unwrap()));
        }
    } else {
        println!("not found");
    }
}

#[test]
fn flower_snark() {
    let mut builder: NestingGraphBuilder<(), ()> = NestingGraphBuilder::new();
    let n1 = builder.add_node(());
    let n2 = builder.add_node(());
    let n3 = builder.add_node(());
    let n4 = builder.add_node(());
    let n5 = builder.add_node(());
    let n6 = builder.add_node(());
    let n7 = builder.add_node(());
    let n8 = builder.add_node(());
    let n9 = builder.add_node(());
    let n10 = builder.add_node(());
    let n11 = builder.add_node(());
    let n12 = builder.add_node(());
    let n13 = builder.add_node(());
    let n14 = builder.add_node(());
    let n15 = builder.add_node(());
    let n16 = builder.add_node(());
    let n17 = builder.add_node(());
    let n18 = builder.add_node(());
    let n19 = builder.add_node(());
    let n20 = builder.add_node(());

    builder.add_edge(n1, n2, ());
    builder.add_edge(n2, n3, ());
    builder.add_edge(n3, n4, ());
    builder.add_edge(n4, n5, ());
    builder.add_edge(n5, n1, ());

    builder.add_edge(n6, n1, ()); // center
    builder.add_edge(n6, n7, ()); //next

    builder.add_edge(n7, n17, ()); //+10
    builder.add_edge(n7, n8, ()); //next

    builder.add_edge(n8, n13, ()); //+5
    builder.add_edge(n8, n9, ()); //next

    builder.add_edge(n9, n2, ()); //center
    builder.add_edge(n9, n10, ()); //next

    builder.add_edge(n10, n20, ()); //+10
    builder.add_edge(n10, n11, ()); //next

    builder.add_edge(n11, n16, ()); //+5
    builder.add_edge(n11, n12, ()); //next

    builder.add_edge(n12, n3, ()); //center
    builder.add_edge(n12, n13, ()); //next

    builder.add_edge(n13, n14, ()); //next

    builder.add_edge(n14, n19, ()); //+5
    builder.add_edge(n14, n15, ()); //next

    builder.add_edge(n15, n4, ()); //center
    builder.add_edge(n15, n16, ()); //next

    builder.add_edge(n16, n17, ()); //next

    builder.add_edge(n17, n18, ()); //next

    builder.add_edge(n18, n5, ()); //center
    builder.add_edge(n18, n19, ()); //next

    builder.add_edge(n19, n20, ()); //next

    builder.add_edge(n20, n6, ()); //next

    let graph = builder.build();

    println!("{}", graph.base_dot());

    // assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());

    println!("loop count {}", graph.cycle_basis().len());
    println!("cycle count {}", graph.all_cycles().len());
    println!(
        "loop count {}",
        graph.paton_count_loops(&graph.full_graph(), 0).unwrap()
    );
    if let Some((s, v)) = graph
        .all_spinneys()
        .iter()
        .find(|(s, _)| graph.full_filter() == s.filter)
    {
        println!(
            "{}",
            graph.dot(&graph.nesting_node_from_subgraph(s.clone()))
        );
        for (ci, cj) in v {
            println!("{}", graph.dot(&ci));
            println!("{}", graph.dot(&cj.as_ref().unwrap()));
        }
    } else {
        println!("not found");
    }
}

use super::*;

use std::cmp::Ordering;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct TestNode(&'static str);

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
