#[test]
fn threeloop() {
    let mut builder: HedgeGraphBuilder<(), ()> = HedgeGraphBuilder::new();
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
                .paton_cycle_basis(&graph.full_graph(), graph.get_incident_node_id(i))
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
                .paton_cycle_basis(&graph.full_graph(), graph.get_incident_node_id(i))
                .unwrap()
                .0
                .len()
        );
    }

    let cycles = graph.cycle_basis().0;

    insta::assert_ron_snapshot!("hairy_three_loop_cycles", cycles);
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

#[test]
#[should_panic]
fn read_tarjan_vs_cycle_space() {
    for s in 0..100 {
        let rand_graph = HedgeGraph::<(), ()>::random(6, 9, s);

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
    graph.full_node();
    println!("built");

    println!("{}", graph.dot(&graph.full_node()));

    let t1 = TraversalTree::dfs(
        &graph,
        &graph.full_filter(),
        graph.nodes.get_index(4).unwrap().0,
    );

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
                Parent::Hedge(p, q) => Some(format!("parent {} \t | rank {} \t|", p, q)),
                Parent::Unset => None,
            },
            &|_| None
        )
    );
    println!("{}", graph.dot(&graph.nesting_node_from_subgraph(t1.tree)));
    println!(
        "{}",
        graph
            .paton_count_loops(&graph.full_graph(), graph.get_incident_node_id(0))
            .unwrap()
    );

    println!("{}", graph.cyclotomatic_number(&graph.full_graph()));

    let cycles = graph
        .paton_cycle_basis(&graph.full_graph(), graph.nodes.get_index(4).unwrap().0)
        .unwrap()
        .0;

    for c in cycles {
        println!("{}", graph.dot(&graph.nesting_node_from_subgraph(c)));
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

    println!("loop count {}", graph.cycle_basis().0.len());
    println!("cycle count {}", graph.all_cycles().len());
    println!("cycle count alt {}", graph.read_tarjan().len());
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

    println!(
        "loop count {}",
        graph.cyclotomatic_number(&graph.full_graph())
    );
    println!("cycle count {}", graph.all_cycles().len());
    println!(
        "loop count {}",
        graph
            .paton_count_loops(&graph.full_graph(), graph.get_incident_node_id(0))
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

use std::time::Instant;

use super::*;
