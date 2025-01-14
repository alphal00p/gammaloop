#[test]
fn threeloop() {
    let mut builder: HedgeGraphBuilder<(), ()> = HedgeGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());

    builder.add_edge(a, b, (), false);
    builder.add_edge(b, a, (), false);
    builder.add_edge(a, b, (), false);

    builder.add_edge(b, c, (), false);
    builder.add_edge(c, d, (), false);
    builder.add_edge(d, a, (), false);

    let graph = builder.build();

    insta::assert_snapshot!("three_loop_dot", graph.base_dot());
    insta::assert_ron_snapshot!("three_loop", graph);

    for i in 0..graph.n_hedges() {
        assert_eq!(
            3,
            graph
                .paton_cycle_basis(&graph.full_graph(), graph.node_hairs(Hedge(i)), None)
                .unwrap()
                .0
                .len()
        );
    }

    let cycles = graph.cycle_basis().0;

    assert_eq!(3, cycles.len());

    let all_cycles = graph.all_cycles();

    assert_eq!(6, all_cycles.len());

    insta::assert_ron_snapshot!("three_loop_cycles", cycles);
}

#[test]
fn hairythreeloop() {
    let mut builder: HedgeGraphBuilder<(), ()> = HedgeGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());

    builder.add_edge(a, b, (), false);
    builder.add_edge(b, a, (), false);
    builder.add_edge(a, b, (), false);
    builder.add_external_edge(a, (), false, Flow::Sink);
    builder.add_external_edge(b, (), false, Flow::Sink);
    builder.add_external_edge(b, (), false, Flow::Sink);

    builder.add_edge(b, c, (), false);
    builder.add_edge(c, d, (), false);
    builder.add_edge(d, a, (), false);

    assert_eq!(builder.involution.len(), 15);
    let graph = builder.build();

    insta::assert_snapshot!("hairy_three_loop_dot", graph.base_dot());
    insta::assert_ron_snapshot!("hairy_three_loop", graph);
    insta::assert_snapshot!(
        "hairy_three_loop_dot_internal",
        graph.dot(&graph.full_node())
    );

    for i in graph.full_node().internal_graph.filter.included_iter() {
        assert_eq!(
            3,
            graph
                .paton_cycle_basis(&graph.full_graph(), graph.node_hairs(i), None)
                .unwrap()
                .0
                .len()
        );
    }

    let cycles = graph.cycle_basis().0;

    insta::assert_ron_snapshot!("hairy_three_loop_cycles", cycles);
}

#[test]
fn banana_cuts() {
    let mut builder: HedgeGraphBuilder<(), ()> = HedgeGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    builder.add_edge(a, b, (), false);
    builder.add_edge(a, b, (), false);
    builder.add_edge(a, b, (), false);

    let three_banana = builder.clone().build();

    assert_eq!(6, three_banana.non_cut_edges().len());
    builder.add_edge(a, b, (), false);

    let four_banana = builder.build();
    assert_eq!(14, four_banana.non_cut_edges().len());
}

#[test]
fn three_loop_fly() {
    let mut builder: HedgeGraphBuilder<(), ()> = HedgeGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());
    builder.add_edge(a, b, (), false);
    builder.add_edge(b, a, (), false);
    builder.add_edge(b, c, (), false);
    builder.add_edge(d, a, (), false);

    builder.add_edge(c, d, (), false);
    builder.add_edge(d, c, (), false);

    let fly = builder.clone().build();
    assert_eq!(32, fly.non_cut_edges().len());
}

#[test]
fn doubletriangle() {
    let mut builder: HedgeGraphBuilder<(), ()> = HedgeGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());
    builder.add_edge(a, b, (), true);
    builder.add_edge(b, c, (), true);
    builder.add_edge(d, a, (), true);

    builder.add_edge(c, d, (), true);
    builder.add_edge(a, c, (), true);

    let fly = builder.clone().build();

    for _c in fly.non_cut_edges() {
        // println!("{c:?}");
        //
        // sum += (2 ^ (c.count_ones() / 2)) / 2;

        // println!("{}", fly.dot(&c));
    }
    // println!("{sum}");
    assert_eq!(13, fly.non_cut_edges().len());
    // println!("{}", SignedCut::all_initial_state_cuts(&fly).len());

    // for c in SignedCut::all_initial_state_cuts(&fly) {
    //     println!("//{}", c.bare_signature(&fly));
    //     println!("{}", fly.dot(&c.cut_content));
    // }
}

#[test]
fn cube() {
    let mut builder: HedgeGraphBuilder<(), ()> = HedgeGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());
    let e = builder.add_node(());
    let f = builder.add_node(());
    let g = builder.add_node(());
    let h = builder.add_node(());

    builder.add_edge(a, b, (), false);
    builder.add_edge(b, c, (), false);
    builder.add_edge(c, d, (), false);
    builder.add_edge(d, a, (), false);

    builder.add_edge(e, f, (), false);
    builder.add_edge(f, g, (), false);
    builder.add_edge(g, h, (), false);
    builder.add_edge(h, e, (), false);

    builder.add_edge(a, e, (), false);
    builder.add_edge(b, f, (), false);
    builder.add_edge(c, g, (), false);
    builder.add_edge(d, h, (), false);

    let graph = builder.build();

    insta::assert_snapshot!("cube_dot", graph.base_dot());

    // let mut all_spinneys = graph.all_spinneys().into_iter().collect_vec();
    // all_spinneys.sort_by(|a, b| a.filter.cmp(&b.filter));

    // assert_eq!(162, all_spinneys.len());
}

#[test]
#[ignore]
fn alt_vs_pair() {
    for s in 0..100 {
        let rand_graph = HedgeGraph::<(), ()>::random(10, 14, s);

        let before = Instant::now();
        let all_spinneys = rand_graph.all_spinneys();
        let after = before.elapsed();
        let before = Instant::now();
        let all_spinneys_alt = rand_graph.all_spinneys_alt();
        let after_alt = before.elapsed();
        println!("{s} {after:?} {after_alt:?}");

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

// #[test]
// #[should_panic]
// fn read_tarjan_vs_cycle_space() {
//     for s in 0..100 {
//         let rand_graph = HedgeGraph::<(), ()>::random(6, 9, s);

//         let all_cycles = rand_graph.read_tarjan();
//         let all_cycles_alt = rand_graph.all_cycles();

//         assert_eq!(
//             all_cycles.len(),
//             all_cycles_alt.len(),
//             "{} with seed {s}",
//             rand_graph.base_dot()
//         );
//     }
// }

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
#[allow(non_snake_case)]
#[test]
fn K33() {
    let mut builder: HedgeGraphBuilder<(), ()> = HedgeGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());
    let e = builder.add_node(());
    let f = builder.add_node(());

    builder.add_edge(a, d, (), false);
    builder.add_edge(a, e, (), false);
    builder.add_edge(a, f, (), false);

    builder.add_edge(b, d, (), false);
    builder.add_edge(b, e, (), false);
    builder.add_edge(b, f, (), false);

    builder.add_edge(c, d, (), false);
    builder.add_edge(c, e, (), false);
    builder.add_edge(c, f, (), false);

    let graph = builder.build();
    graph.full_node();
    println!("built");

    println!("{}", graph.dot(&graph.full_node()));

    let t1 = TraversalTree::dfs(&graph, &graph.full_filter(), &graph.nodes[4], None);

    println!(
        "{}",
        graph.dot(&graph.nesting_node_from_subgraph(t1.tree.clone()))
    );
    println!("{:?}", t1.traversal);

    println!(
        "{}",
        t1.inv.print(
            &graph.full_filter(),
            &|a| match a {
                Parent::Root => Some("Root       \t |     \t|".to_string()),
                Parent::Hedge {
                    hedge_to_root,
                    traversal_order,
                } => Some(format!(
                    "parent {} \t | rank {} \t|",
                    hedge_to_root, traversal_order
                )),
                Parent::Unset => None,
            },
            &|_| None
        )
    );
    println!("{}", graph.dot(&graph.nesting_node_from_subgraph(t1.tree)));
    println!(
        "{}",
        graph
            .paton_count_loops(&graph.full_graph(), graph.node_hairs(Hedge(0)))
            .unwrap()
    );

    println!("{}", graph.cyclotomatic_number(&graph.full_graph()));

    let cycles = graph
        .paton_cycle_basis(&graph.full_graph(), &graph.nodes[4], None)
        .unwrap()
        .0;

    for c in cycles {
        println!(
            "{}",
            graph.dot(&graph.nesting_node_from_subgraph(c.internal_graph(&graph)))
        );
    }

    // assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());
}

#[test]
fn petersen() {
    let mut builder: HedgeGraphBuilder<(), ()> = HedgeGraphBuilder::new();
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

    builder.add_edge(a, b, (), false);
    builder.add_edge(a, f, (), false);
    builder.add_edge(a, e, (), false);

    builder.add_edge(b, c, (), false);
    builder.add_edge(b, g, (), false);

    builder.add_edge(c, d, (), false);
    builder.add_edge(c, h, (), false);

    builder.add_edge(d, e, (), false);
    builder.add_edge(d, i, (), false);

    builder.add_edge(e, j, (), false);

    builder.add_edge(f, h, (), false);
    builder.add_edge(f, i, (), false);

    builder.add_edge(g, i, (), false);
    builder.add_edge(g, j, (), false);

    builder.add_edge(h, j, (), false);

    let graph = builder.build();

    println!("{}", graph.base_dot());

    // assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());

    println!("loop count {}", graph.cycle_basis().0.len());
    println!("cycle count {}", graph.all_cycles().len());
    println!(
        "loop count alt {}",
        graph.cyclotomatic_number(&graph.full_node().internal_graph)
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
            println!("{}", graph.dot(&ci.to_hairy_subgraph(&graph)));
            println!(
                "{}",
                graph.dot(&cj.as_ref().unwrap().to_hairy_subgraph(&graph))
            );
        }
    } else {
        println!("not found");
    }
}

#[test]
fn wagner_graph() {
    let mut builder: HedgeGraphBuilder<(), ()> = HedgeGraphBuilder::new();
    let n1 = builder.add_node(());
    let n2 = builder.add_node(());
    let n3 = builder.add_node(());
    let n4 = builder.add_node(());
    let n5 = builder.add_node(());
    let n6 = builder.add_node(());
    let n7 = builder.add_node(());
    let n8 = builder.add_node(());

    builder.add_edge(n1, n2, (), false);
    builder.add_edge(n1, n5, (), false);

    builder.add_edge(n2, n3, (), false);
    builder.add_edge(n2, n6, (), false);

    builder.add_edge(n3, n4, (), false);
    builder.add_edge(n3, n7, (), false);

    builder.add_edge(n4, n5, (), false);
    builder.add_edge(n4, n8, (), false);

    builder.add_edge(n5, n6, (), false);

    builder.add_edge(n6, n7, (), false);

    builder.add_edge(n7, n8, (), false);

    builder.add_edge(n8, n1, (), false);

    let graph = builder.build();

    println!("{}", graph.base_dot());

    // assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());

    println!("loop count {}", graph.cycle_basis().0.len());
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
            println!("{}", graph.dot(&ci.to_hairy_subgraph(&graph)));
            println!(
                "{}",
                graph.dot(&cj.as_ref().unwrap().to_hairy_subgraph(&graph))
            );
        }
    } else {
        println!("not found");
    }
}

#[test]
fn flower_snark() {
    let mut builder: HedgeGraphBuilder<(), ()> = HedgeGraphBuilder::new();
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

    builder.add_edge(n1, n2, (), false);
    builder.add_edge(n2, n3, (), false);
    builder.add_edge(n3, n4, (), false);
    builder.add_edge(n4, n5, (), false);
    builder.add_edge(n5, n1, (), false);

    builder.add_edge(n6, n1, (), false); // center
    builder.add_edge(n6, n7, (), false); //next

    builder.add_edge(n7, n17, (), false); //+10
    builder.add_edge(n7, n8, (), false); //next

    builder.add_edge(n8, n13, (), false); //+5
    builder.add_edge(n8, n9, (), false); //next

    builder.add_edge(n9, n2, (), false); //center
    builder.add_edge(n9, n10, (), false); //next

    builder.add_edge(n10, n20, (), false); //+10
    builder.add_edge(n10, n11, (), false); //next

    builder.add_edge(n11, n16, (), false); //+5
    builder.add_edge(n11, n12, (), false); //next

    builder.add_edge(n12, n3, (), false); //center
    builder.add_edge(n12, n13, (), false); //next

    builder.add_edge(n13, n14, (), false); //next

    builder.add_edge(n14, n19, (), false); //+5
    builder.add_edge(n14, n15, (), false); //next

    builder.add_edge(n15, n4, (), false); //center
    builder.add_edge(n15, n16, (), false); //next

    builder.add_edge(n16, n17, (), false); //next

    builder.add_edge(n17, n18, (), false); //next

    builder.add_edge(n18, n5, (), false); //center
    builder.add_edge(n18, n19, (), false); //next

    builder.add_edge(n19, n20, (), false); //next

    builder.add_edge(n20, n6, (), false); //next

    let graph = builder.build();

    println!("{}", graph.base_dot());

    // assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());

    println!(
        "loop count {}",
        graph.cyclotomatic_number(&graph.full_graph())
    );
    println!("cycle count {}", graph.all_cycles().len());
    println!(
        "loop count {}",
        graph
            .paton_count_loops(&graph.full_graph(), graph.node_hairs(Hedge(0)))
            .unwrap()
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
            println!("{}", graph.dot(&ci.to_hairy_subgraph(&graph)));
            println!(
                "{}",
                graph.dot(&cj.as_ref().unwrap().to_hairy_subgraph(&graph))
            );
        }
    } else {
        println!("not found");
    }
}

#[test]
fn join() {
    let mut ab = HedgeGraphBuilder::<(), ()>::new();
    let v1 = ab.add_node(());
    let v2 = ab.add_node(());

    ab.add_edge(v1, v2, (), true);
    ab.add_edge(v2, v1, (), true);
    ab.add_external_edge(v1, (), true, Flow::Sink);
    ab.add_external_edge(v2, (), true, Flow::Source);

    let a = ab.build();
    println!("{}", a.base_dot());
    let b = a.clone();

    let c = a
        .join(b, |af, _, bf, _| af == -bf, |af, ad, _, _| (af, ad))
        .unwrap();

    println!("{}", c.base_dot());
}
use std::time::Instant;

use super::*;
