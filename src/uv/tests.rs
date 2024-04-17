use symbolica::{printer::PrintOptions, representations::Atom};

use crate::tests_from_pytest::load_amplitude_output;

use super::{Involution, NestingGraph, NestingGraphBuilder};

#[test]
#[allow(unused)]
fn lbl() {
    let (model, amplitude) = load_amplitude_output("./src/test_resources/lbl/");

    model
        .export_coupling_replacement_rules("./ignore/", PrintOptions::mathematica())
        .unwrap();
    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    graph.generate_uv();
    graph.generate_loop_momentum_bases();

    let uv_graph = graph.derived_data.uv.clone().unwrap();

    insta::assert_snapshot!("lbl_dot", uv_graph.0.base_dot());

    let lmb = &graph.derived_data.loop_momentum_bases.clone().unwrap()[0];

    let cycles = uv_graph.cycle_basis_from_lmb(lmb);

    let all_cycles = uv_graph.0.read_tarjan();
    assert_eq!(all_cycles.len(), 1);

    // for cycle in all_cycles {
    //     println!("{}", uv_graph.0.dot(&cycle));
    // }

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

    let all_cycles = rand_graph.read_tarjan();

    let all_spinneys = rand_graph.all_cycle_unions();

    println!("all spinneys {}", all_spinneys.len());

    assert_eq!(all_cycles.len(), 20);

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

    let all_cycles = graph.read_tarjan();

    assert_eq!(6, all_cycles.len());

    let all_spinneys = graph.all_cycle_unions();

    println!("all spinneys {}", all_spinneys.len());

    insta::assert_ron_snapshot!("three_loop_cycles", cycles);
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
