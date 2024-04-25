use std::time::Instant;

use itertools::Itertools;
use symbolica::{printer::PrintOptions, representations::Atom};

use crate::{
    tests_from_pytest::load_amplitude_output,
    uv::{PoSet, UVGraph},
};

use super::{BitFilter, Involution, NestingGraph, NestingGraphBuilder};

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

    let mut all_spinneys = graph.all_spinneys().into_iter().collect_vec();
    all_spinneys.sort_by(|a, b| a.filter.cmp(&b.filter));

    assert_eq!(graph.all_spinneys_alt().len(), all_spinneys.len());

    // let posetSpinneys = all_spinneys.iter().collect::<PoSet<_>>();

    // println!("{}", posetSpinneys.to_cover_set().dot());

    insta::assert_toml_snapshot!(
        "three_loop_spinneys",
        all_spinneys
            .into_iter()
            .map(|s| graph.dot(&graph.nesting_node_from_subgraph(s)))
            .collect_vec()
    );
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

    let mut all_spinneys = graph.all_spinneys().into_iter().collect_vec();
    all_spinneys.sort_by(|a, b| a.filter.cmp(&b.filter));

    assert_eq!(162, all_spinneys.len());
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
