#[test]
fn twobonds() {
    let two: DotGraph = dot!(digraph { 0 -> 1
        1 -> 2
        0 -> 3
        2 -> a ->b->c-> 3
        3 -> 0
        2 -> 1
    })
    .unwrap();

    let all_two_bonds = two.graph.all_bonds(&(2..3));
    insta::assert_snapshot!(all_two_bonds.len(),@"10");
    // for i in two.graph.all_bonds(&(2..3)) {
    //     println!("{}", two.dot(&i));
    // }

    two.dot(&two.graph.a_bond(&|c| c.n_included() == 3).unwrap());
}

#[test]
fn connected_tadpoles() {
    let two: DotGraph = dot!(
    digraph {

          0 -> 0
          1 -> 1
    })
    .unwrap();

    assert_eq!(two.connected_components(&two.full_filter()).len(), 2);

    let one: DotGraph = dot!(
    digraph {

          0:0 -> 0:1
          0:2-> 0:3
    })
    .unwrap();

    let mut single = one.full_filter();
    single.sub(Hedge(0));
    single.sub(Hedge(1));

    assert_eq!(one.count_connected_components(&single), 1);
}

#[test]
fn cycle_basis() {
    let two: DotGraph = dot!(
    digraph {
          ext0 [shape=none, label="" flow=sink];
          ext0 -> 2[dir=none color="red"];
          ext1 [shape=none, label="" flow=source];
          ext1 -> 3[dir=none color="blue"];
          2 -> 3[ dir=none color="red:blue;0.5"];
          3 -> 2[ dir=none color="red:blue;0.5"];
    })
    .unwrap();

    two.all_cycle_sym_diffs().unwrap();
    assert_eq!(two.cycle_basis().0.len(), 1);
}

#[test]
fn test_spanning_trees_of_tree() {
    let tree: DotGraph = dot!(
    digraph {
         1->2
         1->0
         2->5
          2 -> 3
          3 -> 4
    })
    .unwrap();
    let trees = tree.all_spanning_forests_of(&tree.full_filter());

    for t in &trees {
        println!("{}", tree.dot(t))
    }

    assert_eq!(trees.len(), 1);
}

#[test]
fn join_mut_simple() {
    let two: DotGraph = dot!(
    digraph {

      0 [label = "∑"];
      1 [label = "S:4"];
      ext0 [shape=none, label="" flow=sink];
      ext0 -> 0[dir=back color="red"];
      ext2 [shape=none, label="" flow=source];
      ext2 -> 0[dir=forward color="blue"];
      1 -> 0[ dir=forward color="red:blue;0.5"];
    })
    .unwrap();

    //with

    let one: DotGraph = dot!(digraph {
      node [shape=circle,height=0.1,label=""];  overlap="scale"; layout="neato";

      0 [label = "S:5"];
      ext0 [shape=none, label="" flow=sink];
      ext0 -> 0[dir=back color="red"];
    })
    .unwrap();

    let mut one = one.graph;

    one.join_mut(
        two.graph,
        |sf, _, of, _| {
            println!("{sf:?}vs{of:?}");
            sf == -of
        },
        |sf, sd, _, _| (sf, sd),
    )
    .unwrap();

    println!("{}", one.base_dot())
}

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

    let graph: HedgeGraph<(), (), NoData, NodeStorageVec<()>> = builder.build();

    insta::assert_snapshot!("three_loop_dot", graph.base_dot());
    #[cfg(feature = "serde")]
    insta::assert_ron_snapshot!("three_loop", graph);

    for i in 0..graph.n_hedges() {
        assert_eq!(
            3,
            graph
                .paton_cycle_basis(&graph.full_graph(), &graph.node_id(Hedge(i)), None)
                .unwrap()
                .0
                .len()
        );
    }

    let (cycles, tree) = graph.cycle_basis();
    assert_eq!(tree.covers(&graph.full_filter()), graph.full_filter());

    assert_eq!(3, cycles.len());

    let all_cycles = graph.all_cycles();

    assert_eq!(6, all_cycles.len());

    #[cfg(feature = "serde")]
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

    let n_h: Hedge = builder.involution.len();
    assert_eq!(n_h, Hedge(15));
    let graph: HedgeGraph<(), (), NoData, NodeStorageVec<()>> = builder.build();

    insta::assert_snapshot!("hairy_three_loop_dot", graph.base_dot());

    #[cfg(feature = "serde")]
    insta::assert_ron_snapshot!("hairy_three_loop", graph);
    insta::assert_snapshot!(
        "hairy_three_loop_dot_internal",
        graph.dot(&graph.full_node())
    );

    for i in graph.full_node().internal_graph.filter.included_iter() {
        assert_eq!(
            3,
            graph
                .paton_cycle_basis(&graph.full_graph(), &graph.node_id(i), None)
                .unwrap()
                .0
                .len()
        );
    }

    #[cfg(feature = "serde")]
    let cycles = graph.cycle_basis().0;

    #[cfg(feature = "serde")]
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

    let three_banana: HedgeGraph<(), ()> = builder.clone().build();

    assert_eq!(6, three_banana.non_cut_edges().len());
    builder.add_edge(a, b, (), false);

    let four_banana: HedgeGraph<(), ()> = builder.build();
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

    let fly: HedgeGraph<(), ()> = builder.clone().build();
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

    let fly: HedgeGraph<(), ()> = builder.clone().build();

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

    let graph: HedgeGraph<(), ()> = builder.build();

    insta::assert_snapshot!("cube_dot", graph.base_dot());

    // let mut all_spinneys = graph.all_spinneys().into_iter().collect_vec();
    // all_spinneys.sort_by(|a, b| a.filter.cmp(&b.filter));

    // assert_eq!(162, all_spinneys.len());
}

#[test]
#[ignore]
fn alt_vs_pair() {
    for s in 0..100 {
        let rand_graph = HedgeGraph::<(), (), (), NodeStorageVec<()>>::random(10, 14, s);

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

    let graph: HedgeGraph<(), ()> = builder.build();
    graph.full_node();
    println!("built");

    println!("{}", graph.dot(&graph.full_node()));

    let t1 = SimpleTraversalTree::depth_first_traverse(
        &graph,
        &graph.full_filter(),
        &NodeIndex(4),
        None,
    )
    .unwrap();

    println!(
        "{}",
        graph.dot(&graph.nesting_node_from_subgraph(t1.tree_subgraph(&graph).clone()))
    );
    // println!("{:?}", t1.traversal);

    // println!(
    //     "{}",
    //     t1.inv.print(
    //         &graph.full_filter(),
    //         &|a| match a {
    //             Parent::Root => Some("Root       \t |     \t|".to_string()),
    //             Parent::Hedge {
    //                 hedge_to_root,
    //                 traversal_order,
    //             } => Some(format!(
    //                 "parent {} \t | rank {} \t|",
    //                 hedge_to_root, traversal_order
    //             )),
    //             Parent::Unset => None,
    //         },
    //         &|_| None
    //     )
    // );
    println!(
        "{}",
        graph.dot(&graph.nesting_node_from_subgraph(t1.tree_subgraph(&graph)))
    );
    println!(
        "{}",
        graph
            .paton_count_loops(&graph.full_graph(), &graph.node_id(Hedge(0)))
            .unwrap()
    );

    println!("{}", graph.cyclotomatic_number(&graph.full_graph()));

    let cycles = graph
        .paton_cycle_basis(&graph.full_graph(), &NodeIndex(4), None)
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

    let graph: HedgeGraph<(), ()> = builder.build();

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

    let graph: HedgeGraph<(), ()> = builder.build();

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

    let graph: HedgeGraph<(), ()> = builder.build();

    println!("{}", graph.base_dot());

    // assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());
    assert_eq!(11, graph.cyclotomatic_number(&graph.full_graph()));
    println!("cycle count {}", graph.all_cycles().len());
    assert_eq!(
        11,
        graph
            .paton_count_loops(&graph.full_graph(), &graph.node_id(Hedge(0)))
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
    let mut ab = HedgeGraphBuilder::new();
    let v1 = ab.add_node("a");
    let v2 = ab.add_node("a");

    ab.add_edge(v1, v2, "ie", true);
    ab.add_edge(v2, v1, "ie", true);
    ab.add_external_edge(v1, "esink", true, Flow::Sink);
    ab.add_external_edge(v2, "esource", true, Flow::Source);

    let a: HedgeGraph<&'static str, &'static str> = ab.build();

    let mut ab = HedgeGraphBuilder::new();
    let v1 = ab.add_node("b");
    let v2 = ab.add_node("b");

    ab.add_edge(v1, v2, "if", true);
    ab.add_edge(v2, v1, "if", true);
    ab.add_external_edge(v1, "f", true, Flow::Sink);
    ab.add_external_edge(v2, "f", true, Flow::Source);

    let b: HedgeGraph<&'static str, &'static str> = ab.build();

    let mut c = a
        .clone()
        .join(b.clone(), |af, _, bf, _| af == -bf, |af, ad, _, _| (af, ad))
        .unwrap();

    let n_len: NodeIndex = c.node_store.len();
    assert_eq!(n_len, NodeIndex(4));
    let n_h: Hedge = c.node_store.len();
    assert_eq!(n_h, Hedge(12));
    assert!(c.node_store.check_and_set_nodes().is_ok());

    assert_snapshot!(c.dot_display(&c.full_filter()));

    let mut a = HedgeGraphBuilder::new();
    let n = a.add_node("a");
    a.add_external_edge(n, "e", true, Flow::Sink);
    let a: HedgeGraph<&'static str, &'static str> = a.build();

    let mut b = HedgeGraphBuilder::new();
    let n = b.add_node("b");
    b.add_external_edge(n, "f", true, Flow::Sink);
    let b: HedgeGraph<&'static str, &'static str> = b.build();
    let c = a
        .join(
            b,
            |_, _, _, _| true,
            |af, ad, bf, bd| {
                println!("af: {af:?}, ad: {ad:?}, bf: {bf:?}, bd: {bd:?}");
                (af, ad)
            },
        )
        .unwrap();

    assert_eq!(
        "e",
        *c.iter_edges_of(&c.full_filter()).next().unwrap().2.data
    );
}
use std::time::Instant;

use dot_parser::ast::CompassPt;
use insta::assert_snapshot;
use nodestore::NodeStorageVec;

use crate::{dot, half_edge::swap::Swap, parser::DotGraph};

use super::*;

#[test]
fn double_triangle() {
    let mut builder = HedgeGraphBuilder::new();
    let a = builder.add_node(());
    let b = builder.add_node(());
    let c = builder.add_node(());
    let d = builder.add_node(());

    builder.add_edge(a, b, (), true);
    builder.add_edge(b, c, (), true);
    builder.add_edge(c, d, (), true);
    builder.add_edge(d, a, (), true);
    builder.add_edge(b, d, (), true);

    let graph: HedgeGraph<(), (), ()> = builder.build();
    let cuts = graph.all_cuts_from_ids(&[a], &[c]);
    assert_eq!(cuts.len(), 4);
    for cut in cuts {
        assert!(!(cut.1.is_empty()))
    }
}

#[test]
fn double_self_loop() {
    let mut builder = HedgeGraphBuilder::new();
    let a = builder.add_node(());

    builder.add_edge(a, a, (), true);
    builder.add_edge(a, a, (), true);

    let graph: HedgeGraph<(), (), ()> = builder.build();
    let (loops, _) = graph.cycle_basis();
    assert_eq!(loops.len(), 2);
}

#[test]
fn self_energy_cut() {
    let mut epem_builder = HedgeGraphBuilder::new();
    let nodes = (0..4)
        .map(|_| epem_builder.add_node(()))
        .collect::<Vec<_>>();

    epem_builder.add_edge(nodes[0], nodes[1], (), true);
    epem_builder.add_edge(nodes[1], nodes[2], (), true);
    epem_builder.add_edge(nodes[1], nodes[2], (), true);
    epem_builder.add_edge(nodes[2], nodes[3], (), true);

    epem_builder.add_external_edge(nodes[0], (), true, Flow::Sink);
    epem_builder.add_external_edge(nodes[0], (), true, Flow::Sink);
    epem_builder.add_external_edge(nodes[3], (), true, Flow::Source);
    epem_builder.add_external_edge(nodes[3], (), true, Flow::Source);

    let epem: HedgeGraph<(), ()> = epem_builder.build();
    println!("{}", epem.dot(&epem.full_filter()));

    let cuts = epem.all_cuts_from_ids(&nodes[0..=0], &nodes[3..=3]);

    assert_eq!(cuts.len(), 3);
}

#[test]
fn double_pentagon_all_cuts() {
    let graph: DotGraph = dot!(
        digraph {
    node [shape=circle,height=0.1,label=""];  overlap="scale"; layout="neato";
    00 -> 07[ dir=none,label="a"];
    00 -> 12[ dir=forward,label="d"];
    01 -> 00[ dir=forward,label="d"];
    01 -> 03[ dir=none,label="a"];
    02 -> 01[ dir=forward,label="d"];
    02 -> 06[ dir=none,label="a"];
    03 -> 13[ dir=forward,label="d"];
    04 -> 03[ dir=forward,label="d"];
    04 -> 05[ dir=none,label="g"];
    05 -> 02[ dir=forward,label="d"];
    06 -> 07[ dir=forward,label="e-"];
    07 -> 11[ dir=forward,label="e-"];
    08 -> 06[ dir=forward,label="e-"];
    09 -> 04[ dir=forward,label="d"];
    10 -> 05[ dir=forward,label="d"];
    }
    )
    .unwrap();

    // println!(
    //     "{}",
    //     graph.dot_impl(&graph.full_filter(), "", &|a| None, &|n| Some(format!(
    //         "label=\"{}\"",
    //         n.id
    //     )))
    // );

    let cuts = graph.all_cuts(
        graph.combine_to_single_hedgenode(&[NodeIndex(10)]),
        graph.combine_to_single_hedgenode(&[NodeIndex(9)]),
    );

    assert_eq!(cuts.len(), 9);

    let cuts = graph.all_cuts(
        graph.combine_to_single_hedgenode(&[NodeIndex(10)]),
        graph.combine_to_single_hedgenode(&[NodeIndex(13)]),
    );

    assert_eq!(cuts.len(), 14);

    let cuts = graph.all_cuts(
        graph.combine_to_single_hedgenode(&[NodeIndex(1)]),
        graph.combine_to_single_hedgenode(&[NodeIndex(2)]),
    );

    assert_eq!(cuts.len(), 16);

    // for (l, c, r) in cuts {
    //     assert_eq!("//cut:\n{}", graph.dot(&c.reference));
    // }
    // let cuts = graph.all_cuts(NodeIndex(10), NodeIndex(13));

    // println!("All cuts: {}", cuts.len());
}

#[test]
fn box_all_cuts_multiple() {
    let graph: DotGraph = dot!(
        digraph G {
         e00->e01
         e01->e02
         e02->e03
         e03->e00
         e10->e00
         e11->e01
         e12->e02
         e13->e03
        }
    )
    .unwrap();

    // println!(
    //     "{}",
    //     graph.dot_impl(&graph.full_filter(), "", &|a| None, &|n| Some(format!(
    //         "label=\"{}\"",
    //         n.id
    //     )))
    // );

    let cuts =
        graph.all_cuts_from_ids(&[NodeIndex(4), NodeIndex(6)], &[NodeIndex(5), NodeIndex(7)]);

    // for (l, c, r) in &cuts {
    //     println!(
    //         "//cut:\n{}",
    //         graph.dot_impl(l, "start=2;\n", &|a| None, &|n| Some(format!(
    //             "label=\"{}\"",
    //             n.id
    //         )))
    //     );
    // }
    // let cuts = graph.all_cuts(NodeIndex(10), NodeIndex(13));

    assert_eq!(11, cuts.len());
}

#[test]
#[cfg(feature = "serde")]
fn self_energy_box() {
    let mut self_energy_builder: HedgeGraphBuilder<(), (), ()> = HedgeGraphBuilder::new();
    let nodes = (0..8)
        .map(|_| self_energy_builder.add_node(()))
        .collect::<Vec<_>>();

    self_energy_builder.add_edge(nodes[0], nodes[2], (), true);
    self_energy_builder.add_edge(nodes[2], nodes[1], (), true);
    self_energy_builder.add_edge(nodes[1], nodes[7], (), true);
    self_energy_builder.add_edge(nodes[7], nodes[5], (), true);
    self_energy_builder.add_edge(nodes[5], nodes[6], (), true);
    self_energy_builder.add_edge(nodes[6], nodes[0], (), true);
    self_energy_builder.add_edge(nodes[3], nodes[2], (), true);
    self_energy_builder.add_edge(nodes[4], nodes[5], (), true);
    self_energy_builder.add_edge(nodes[3], nodes[4], (), true);
    self_energy_builder.add_edge(nodes[4], nodes[3], (), true);

    self_energy_builder.add_external_edge(nodes[0], (), true, Flow::Sink);
    self_energy_builder.add_external_edge(nodes[1], (), true, Flow::Sink);
    self_energy_builder.add_external_edge(nodes[6], (), true, Flow::Source);
    self_energy_builder.add_external_edge(nodes[7], (), true, Flow::Source);

    let mut cut_to_look_for = vec![
        EdgeIndex::from(5),
        EdgeIndex::from(2),
        EdgeIndex::from(8),
        EdgeIndex::from(9),
    ];

    cut_to_look_for.sort();

    let self_energy: HedgeGraph<(), (), _> = self_energy_builder.build();

    let cuts = self_energy.all_cuts_from_ids(&[nodes[0]], &[nodes[7]]);

    for (left, cut, right) in cuts.iter() {
        let mut edges_in_cut: Vec<_> = self_energy
            .iter_edges_of(cut)
            .map(|(_, edge, _)| edge)
            .collect();

        edges_in_cut.sort();

        if cut_to_look_for == edges_in_cut {
            insta::assert_ron_snapshot!(left);

            insta::assert_ron_snapshot!(right);

            insta::assert_ron_snapshot!(cut.left);
        }
    }
}

#[test]
fn tadpoles() {
    let graph: DotGraph = dot!(
        digraph{
            a->b->c
            b->aa->d
            aa->ab
            ab->ab

            e->d->f
        }
    )
    .unwrap();

    let tads = graph.tadpoles(&[NodeIndex(0), NodeIndex(4), NodeIndex(6), NodeIndex(7)]);

    for t in tads {
        println!("//Tadpole: \n{}", graph.dot(&t));
    }
    println!("Graph: {}", graph.base_dot());
}

#[test]
fn extracting_network() {
    let mut graph: DotGraph = dot!(
        digraph {
          node	 [shape=circle,height=0.1,label=""];
          overlap = "scale";
          layout = "neato";

          0	 [label = "∑"];
          75	 [label = "∏"];
          1	 [label = "S:88"];
          9	 [label = "T:88"];
          3	 [label = "S:0"];
          4	 [label = "S:1"];
          5	 [label = "T:0"];
          6	 [label = "T:1"];
          7	 [label = "T:2"];
          8	 [label = "S:2"];
          2	 [label = "∏"];
          20	 [label = "T:87"];
          11	 [label = "S:3"];
          12	 [label = "S:4"];
          13	 [label = "S:5"];
          14	 [label = "T:3"];
          29	 [label = "T:86"];
          16	 [label = "T:4"];
          17	 [label = "T:5"];
          18	 [label = "T:6"];
          19	 [label = "S:6"];
          10	 [label = "∏"];
          38	 [label = "T:85"];
          22	 [label = "S:7"];
          23	 [label = "S:8"];
          24	 [label = "S:9"];
          25	 [label = "T:7"];
          26	 [label = "T:8"];
          27	 [label = "T:9"];
          28	 [label = "S:10"];
          15	 [label = "∏"];
          47	 [label = "T:84"];
          31	 [label = "S:11"];
          32	 [label = "S:12"];
          33	 [label = "S:13"];
          34	 [label = "T:10"];
          35	 [label = "T:11"];
          36	 [label = "T:12"];
          37	 [label = "S:14"];
          21	 [label = "∏"];
          56	 [label = "T:83"];
          40	 [label = "S:15"];
          41	 [label = "S:16"];
          42	 [label = "S:17"];
          43	 [label = "T:13"];
          44	 [label = "T:14"];
          45	 [label = "T:15"];
          46	 [label = "S:18"];
          30	 [label = "∏"];
          49	 [label = "S:19"];
          50	 [label = "S:20"];
          51	 [label = "S:21"];
          52	 [label = "T:16"];
          53	 [label = "T:17"];
          54	 [label = "T:18"];
          55	 [label = "S:22"];
          39	 [label = "∏"];
          65	 [label = "S:87"];
          74	 [label = "S:86"];
          58	 [label = "S:23"];
          59	 [label = "S:24"];
          60	 [label = "S:25"];
          61	 [label = "T:19"];
          62	 [label = "T:20"];
          63	 [label = "T:21"];
          64	 [label = "S:26"];
          48	 [label = "∏"];
          83	 [label = "S:85"];
          84	 [label = "∏"];
          67	 [label = "S:27"];
          68	 [label = "S:28"];
          69	 [label = "S:29"];
          70	 [label = "T:22"];
          71	 [label = "T:23"];
          72	 [label = "T:24"];
          73	 [label = "S:30"];
          57	 [label = "∏"];
          76	 [label = "S:31"];
          77	 [label = "S:32"];
          78	 [label = "S:33"];
          79	 [label = "T:25"];
          80	 [label = "T:26"];
          81	 [label = "T:27"];
          82	 [label = "S:34"];
          66	 [label = "∏"];
          85	 [label = "S:35"];
          86	 [label = "S:36"];
          87	 [label = "S:37"];
          88	 [label = "T:28"];
          89	 [label = "T:29"];
          90	 [label = "T:30"];
          91	 [label = "T:31"];
          92	 [label = "S:38"];
          94	 [label = "∏"];
          95	 [label = "S:39"];
          96	 [label = "S:40"];
          97	 [label = "S:41"];
          98	 [label = "T:32"];
          99	 [label = "T:33"];
          100	 [label = "T:34"];
          101	 [label = "T:35"];
          102	 [label = "S:42"];
          104	 [label = "∏"];
          105	 [label = "S:43"];
          106	 [label = "S:44"];
          107	 [label = "S:45"];
          108	 [label = "T:36"];
          109	 [label = "T:37"];
          110	 [label = "T:38"];
          111	 [label = "T:39"];
          112	 [label = "S:46"];
          114	 [label = "∏"];
          115	 [label = "S:47"];
          116	 [label = "S:48"];
          117	 [label = "S:49"];
          118	 [label = "S:50"];
          119	 [label = "T:40"];
          120	 [label = "T:41"];
          121	 [label = "T:42"];
          122	 [label = "S:51"];
          124	 [label = "∏"];
          93	 [label = "S:84"];
          103	 [label = "T:82"];
          125	 [label = "S:52"];
          126	 [label = "S:53"];
          127	 [label = "S:54"];
          128	 [label = "S:55"];
          129	 [label = "T:43"];
          130	 [label = "S:56"];
          131	 [label = "T:44"];
          132	 [label = "T:45"];
          133	 [label = "S:57"];
          135	 [label = "∏"];
          113	 [label = "T:81"];
          136	 [label = "S:58"];
          137	 [label = "S:59"];
          138	 [label = "S:60"];
          139	 [label = "T:46"];
          140	 [label = "T:47"];
          141	 [label = "T:48"];
          142	 [label = "T:49"];
          143	 [label = "T:50"];
          144	 [label = "S:61"];
          146	 [label = "∏"];
          123	 [label = "T:80"];
          147	 [label = "S:62"];
          148	 [label = "S:63"];
          149	 [label = "S:64"];
          150	 [label = "T:51"];
          151	 [label = "T:52"];
          152	 [label = "T:53"];
          153	 [label = "T:54"];
          154	 [label = "T:55"];
          155	 [label = "S:65"];
          157	 [label = "∏"];
          134	 [label = "T:79"];
          158	 [label = "S:66"];
          159	 [label = "S:67"];
          160	 [label = "S:68"];
          161	 [label = "T:56"];
          162	 [label = "T:57"];
          163	 [label = "T:58"];
          164	 [label = "T:59"];
          165	 [label = "T:60"];
          166	 [label = "S:69"];
          168	 [label = "∏"];
          145	 [label = "T:78"];
          169	 [label = "S:70"];
          170	 [label = "S:71"];
          171	 [label = "S:72"];
          172	 [label = "T:61"];
          173	 [label = "T:62"];
          174	 [label = "T:63"];
          175	 [label = "T:64"];
          176	 [label = "T:65"];
          177	 [label = "S:73"];
          178	 [label = "∏"];
          156	 [label = "T:77"];
          180	 [label = "S:74"];
          181	 [label = "S:75"];
          182	 [label = "S:76"];
          183	 [label = "T:66"];
          184	 [label = "T:67"];
          185	 [label = "T:68"];
          186	 [label = "T:69"];
          187	 [label = "T:70"];
          179	 [label = "∏"];
          190	 [label = "S:77"];
          191	 [label = "S:78"];
          192	 [label = "S:79"];
          193	 [label = "T:71"];
          194	 [label = "T:72"];
          195	 [label = "T:73"];
          196	 [label = "T:74"];
          197	 [label = "T:75"];
          198	 [label = "T:76"];
          199	 [label = "S:80"];
          189	 [label = "∏"];
          167	 [label = "S:83"];
          200	 [label = "∏"];
          188	 [label = "S:81"];
          ext0	 [style=invis];
            0:0	-> ext0	 [id=0 color="gray"];
            ext23	 [style=invis];
            0:23	-> ext23	 [id=23 dir=none color="gray" label="mink4|0"];
            ext46	 [style=invis];
            0:46	-> ext46	 [id=46 dir=none color="gray" label="mink4|1"];
            ext69	 [style=invis];
            0:69	-> ext69	 [id=69 dir=none color="gray" label="mink4|4"];
            75:92	-> 0:22	 [id=98  color="gray"];
            1:97:s	-> 84:355:s	 [id=70  color="red:blue;0.5"];
            3:101	-> 75:100	 [id=92  color="gray"];
            4:102	-> 75:99	 [id=97  color="gray"];
            5:103	-> 75:96	 [id=95  color="gray"];
            5:104	-> 0:45	 [id=22 dir=none  color="gray" label="mink4|0"];
            6:105	-> 75:95	 [id=94  color="gray"];
            6:106	-> 0:68	 [id=45 dir=none  color="gray" label="mink4|1"];
            7:107	-> 75:94	 [id=93  color="gray"];
            7:108	-> 0:91	 [id=68 dir=none  color="gray" label="mink4|4"];
            8:109	-> 75:93	 [id=91  color="gray"];
            2:110	-> 0:21	 [id=109  color="gray"];
            9:113:s	-> 84:331:s	 [id=166  color="red:blue;0.5"];
            20:116:s	-> 84:330:s	 [id=138  color="red:blue;0.5"];
            11:120	-> 2:119	 [id=105  color="gray"];
            12:121	-> 2:118	 [id=106  color="gray"];
            13:122	-> 2:117	 [id=104  color="gray"];
            14:123	-> 2:114	 [id=99  color="gray"];
            14:124	-> 0:44	 [id=21 dir=none  color="gray" label="mink4|0"];
            14:125	-> 0:67	 [id=44 dir=none  color="gray" label="mink4|1"];
            29:126:s	-> 0:70	 [id=47 dir=none  color="red:gray;0.5" label="mink4|4"];
            16:129	-> 2:128	 [id=108  color="gray"];
            16:130	-> 17:132	 [id=103 dir=none  color="gray" label="mink4|i"];
            17:131	-> 2:127	 [id=102  color="gray"];
            18:133	-> 2:112	 [id=100  color="gray"];
            18:134	-> 0:90	 [id=67 dir=none  color="gray" label="mink4|4"];
            19:135	-> 2:111	 [id=90  color="gray"];
            10:136	-> 0:20	 [id=117  color="gray"];
            29:141:s	-> 84:307:s	 [id=129  color="red:blue;0.5"];
            38:142:s	-> 9:98:s	 [id=161 dir=none  color="red:blue;0.5" label="mink4|26"];
            22:146	-> 10:145	 [id=116  color="gray"];
            23:147	-> 10:144	 [id=110  color="gray"];
            24:148	-> 10:143	 [id=115  color="gray"];
            25:149	-> 10:140	 [id=113  color="gray"];
            25:150	-> 0:43	 [id=20 dir=none  color="gray" label="mink4|0"];
            26:151	-> 10:139	 [id=112  color="gray"];
            26:152	-> 0:66	 [id=43 dir=none  color="gray" label="mink4|1"];
            27:153	-> 10:138	 [id=111  color="gray"];
            27:154	-> 0:89	 [id=66 dir=none  color="gray" label="mink4|4"];
            28:155	-> 10:137	 [id=89  color="gray"];
            15:156	-> 0:19	 [id=125  color="gray"];
            38:161:s	-> 84:306:s	 [id=146  color="red:blue;0.5"];
            47:162:s	-> 20:115:s	 [id=153 dir=none  color="red:blue;0.5" label="mink4|5"];
            31:166	-> 15:165	 [id=124  color="gray"];
            32:167	-> 15:164	 [id=118  color="gray"];
            33:168	-> 15:163	 [id=123  color="gray"];
            34:169	-> 15:160	 [id=121  color="gray"];
            34:170	-> 0:42	 [id=19 dir=none  color="gray" label="mink4|0"];
            35:171	-> 15:159	 [id=120  color="gray"];
            35:172	-> 0:65	 [id=42 dir=none  color="gray" label="mink4|1"];
            36:173	-> 15:158	 [id=119  color="gray"];
            36:174	-> 0:88	 [id=65 dir=none  color="gray" label="mink4|4"];
            37:175	-> 15:157	 [id=88  color="gray"];
            21:176	-> 0:18	 [id=133  color="gray"];
            47:181:s	-> 84:283:s	 [id=122  color="red:blue;0.5"];
            56:182:s	-> 0:47	 [id=24 dir=none  color="red:gray;0.5" label="mink4|1"];
            40:186	-> 21:185	 [id=132  color="gray"];
            41:187	-> 21:184	 [id=126  color="gray"];
            42:188	-> 21:183	 [id=131  color="gray"];
            43:189	-> 21:180	 [id=130  color="gray"];
            43:190	-> 0:64	 [id=41 dir=none  color="gray" label="mink4|1"];
            44:191	-> 21:179	 [id=128  color="gray"];
            44:192	-> 0:41	 [id=18 dir=none  color="gray" label="mink4|0"];
            45:193	-> 21:178	 [id=127  color="gray"];
            45:194	-> 0:87	 [id=64 dir=none  color="gray" label="mink4|4"];
            46:195	-> 21:177	 [id=87  color="gray"];
            30:196	-> 0:17	 [id=141  color="gray"];
            56:201:s	-> 0:24	 [id=1 dir=none  color="red:gray;0.5" label="mink4|0"];
            56:202:s	-> 84:282:s	 [id=176  color="red:blue;0.5"];
            49:206	-> 30:205	 [id=140  color="gray"];
            50:207	-> 30:204	 [id=134  color="gray"];
            51:208	-> 30:203	 [id=139  color="gray"];
            52:209	-> 30:200	 [id=137  color="gray"];
            52:210	-> 0:40	 [id=17 dir=none  color="gray" label="mink4|0"];
            53:211	-> 30:199	 [id=136  color="gray"];
            53:212	-> 0:63	 [id=40 dir=none  color="gray" label="mink4|1"];
            54:213	-> 30:198	 [id=135  color="gray"];
            54:214	-> 0:86	 [id=63 dir=none  color="gray" label="mink4|4"];
            55:215	-> 30:197	 [id=86  color="gray"];
            39:216	-> 0:16	 [id=149  color="gray"];
            65:221:s	-> 84:262:s	 [id=114  color="red:blue;0.5"];
            74:222:s	-> 84:261:s	 [id=101  color="red:blue;0.5"];
            58:226	-> 39:225	 [id=148  color="gray"];
            59:227	-> 39:224	 [id=142  color="gray"];
            60:228	-> 39:223	 [id=147  color="gray"];
            61:229	-> 39:220	 [id=145  color="gray"];
            61:230	-> 0:39	 [id=16 dir=none  color="gray" label="mink4|0"];
            62:231	-> 39:219	 [id=144  color="gray"];
            62:232	-> 0:62	 [id=39 dir=none  color="gray" label="mink4|1"];
            63:233	-> 39:218	 [id=143  color="gray"];
            63:234	-> 0:85	 [id=62 dir=none  color="gray" label="mink4|4"];
            64:235	-> 39:217	 [id=85  color="gray"];
            48:236	-> 0:15	 [id=157  color="gray"];
            83:241:s	-> 84:242:s	 [id=107  color="red:blue;0.5"];
            67:246	-> 48:245	 [id=156  color="gray"];
            68:247	-> 48:244	 [id=150  color="gray"];
            69:248	-> 48:243	 [id=155  color="gray"];
            70:249	-> 48:240	 [id=154  color="gray"];
            70:250	-> 0:61	 [id=38 dir=none  color="gray" label="mink4|1"];
            71:251	-> 48:239	 [id=152  color="gray"];
            71:252	-> 0:38	 [id=15 dir=none  color="gray" label="mink4|0"];
            72:253	-> 48:238	 [id=151  color="gray"];
            72:254	-> 0:84	 [id=61 dir=none  color="gray" label="mink4|4"];
            73:255	-> 48:237	 [id=84  color="gray"];
            57:256	-> 0:14	 [id=165  color="gray"];
            76:266	-> 57:265	 [id=164  color="gray"];
            77:267	-> 57:264	 [id=158  color="gray"];
            78:268	-> 57:263	 [id=163  color="gray"];
            79:269	-> 57:260	 [id=162  color="gray"];
            79:270	-> 0:60	 [id=37 dir=none  color="gray" label="mink4|1"];
            80:271	-> 57:259	 [id=160  color="gray"];
            80:272	-> 0:37	 [id=14 dir=none  color="gray" label="mink4|0"];
            81:273	-> 57:258	 [id=159  color="gray"];
            81:274	-> 0:83	 [id=60 dir=none  color="gray" label="mink4|4"];
            82:275	-> 57:257	 [id=83  color="gray"];
            66:276	-> 0:13	 [id=175  color="gray"];
            85:287	-> 66:286	 [id=173  color="gray"];
            86:288	-> 66:285	 [id=174  color="gray"];
            87:289	-> 66:284	 [id=172  color="gray"];
            88:290	-> 66:281	 [id=171  color="gray"];
            88:291	-> 0:36	 [id=13 dir=none  color="gray" label="mink4|0"];
            88:292	-> 0:59	 [id=36 dir=none  color="gray" label="mink4|1"];
            89:293	-> 66:280	 [id=170  color="gray"];
            89:294	-> 90:296	 [id=169 dir=none  color="gray" label="mink4|5"];
            90:295	-> 66:279	 [id=168  color="gray"];
            91:297	-> 66:278	 [id=167  color="gray"];
            91:298	-> 0:82	 [id=59 dir=none  color="gray" label="mink4|4"];
            92:299	-> 66:277	 [id=82  color="gray"];
            94:300	-> 0:12	 [id=185  color="gray"];
            95:311	-> 94:310	 [id=183  color="gray"];
            96:312	-> 94:309	 [id=184  color="gray"];
            97:313	-> 94:308	 [id=182  color="gray"];
            98:314	-> 94:305	 [id=181  color="gray"];
            98:315	-> 0:35	 [id=12 dir=none  color="gray" label="mink4|0"];
            98:316	-> 0:58	 [id=35 dir=none  color="gray" label="mink4|1"];
            99:317	-> 94:304	 [id=179  color="gray"];
            99:318	-> 101:322	 [id=178 dir=none  color="gray" label="mink4|5"];
            100:319	-> 94:303	 [id=180  color="gray"];
            100:320	-> 0:81	 [id=58 dir=none  color="gray" label="mink4|4"];
            101:321	-> 94:302	 [id=177  color="gray"];
            102:323	-> 94:301	 [id=81  color="gray"];
            104:324	-> 0:11	 [id=195  color="gray"];
            105:335	-> 104:334	 [id=193  color="gray"];
            106:336	-> 104:333	 [id=194  color="gray"];
            107:337	-> 104:332	 [id=192  color="gray"];
            108:338	-> 104:329	 [id=191  color="gray"];
            108:339	-> 0:34	 [id=11 dir=none  color="gray" label="mink4|0"];
            108:340	-> 0:57	 [id=34 dir=none  color="gray" label="mink4|1"];
            109:341	-> 104:328	 [id=190  color="gray"];
            109:342	-> 110:344	 [id=189 dir=none  color="gray" label="mink4|5"];
            110:343	-> 104:327	 [id=188  color="gray"];
            111:345	-> 104:326	 [id=187  color="gray"];
            111:346	-> 0:80	 [id=57 dir=none  color="gray" label="mink4|4"];
            112:347	-> 104:325	 [id=80  color="gray"];
            114:348	-> 0:10	 [id=204  color="gray"];
            84:356:s	-> 0:1	 [id=96  color="red:gray;0.5"];
            115:359	-> 114:358	 [id=196  color="gray"];
            116:360	-> 114:357	 [id=203  color="gray"];
            117:361	-> 114:354	 [id=201  color="gray"];
            118:362	-> 114:353	 [id=199  color="gray"];
            119:363	-> 114:352	 [id=200  color="gray"];
            119:364	-> 0:56	 [id=33 dir=none  color="gray" label="mink4|1"];
            120:365	-> 114:351	 [id=198  color="gray"];
            120:366	-> 0:33	 [id=10 dir=none  color="gray" label="mink4|0"];
            121:367	-> 114:350	 [id=197  color="gray"];
            121:368	-> 0:79	 [id=56 dir=none  color="gray" label="mink4|4"];
            122:369	-> 114:349	 [id=79  color="gray"];
            124:370	-> 0:9	 [id=214  color="gray"];
            93:377	-> 189:553	 [id=71  color="gray"];
            125:382	-> 124:381	 [id=213  color="gray"];
            126:383	-> 124:380	 [id=205  color="gray"];
            127:384	-> 124:379	 [id=212  color="gray"];
            128:385	-> 124:376	 [id=210  color="gray"];
            129:386	-> 124:375	 [id=209  color="gray"];
            129:387	-> 0:32	 [id=9 dir=none  color="gray" label="mink4|0"];
            130:388	-> 124:374	 [id=208  color="gray"];
            131:389	-> 124:373	 [id=207  color="gray"];
            131:390	-> 0:55	 [id=32 dir=none  color="gray" label="mink4|1"];
            132:391	-> 124:372	 [id=206  color="gray"];
            132:392	-> 0:78	 [id=55 dir=none  color="gray" label="mink4|4"];
            133:393	-> 124:371	 [id=78  color="gray"];
            135:394	-> 0:8	 [id=225  color="gray"];
            103:401	-> 189:554	 [id=283  color="gray"];
            136:406	-> 135:405	 [id=223  color="gray"];
            137:407	-> 135:404	 [id=224  color="gray"];
            138:408	-> 135:403	 [id=222  color="gray"];
            139:409	-> 135:400	 [id=219  color="gray"];
            139:410	-> 0:54	 [id=31 dir=none  color="gray" label="mink4|1"];
            140:411	-> 135:399	 [id=221  color="gray"];
            140:412	-> 143:418	 [id=217 dir=none  color="gray" label="mink4|26"];
            141:413	-> 135:398	 [id=218  color="gray"];
            141:414	-> 0:31	 [id=8 dir=none  color="gray" label="mink4|0"];
            142:415	-> 135:397	 [id=220  color="gray"];
            142:416	-> 0:77	 [id=54 dir=none  color="gray" label="mink4|4"];
            143:417	-> 135:396	 [id=216  color="gray"];
            144:419	-> 135:395	 [id=77  color="gray"];
            146:420	-> 0:7	 [id=236  color="gray"];
            113:427	-> 189:555	 [id=237  color="gray"];
            123:428	-> 0:71	 [id=48 dir=none  color="gray" label="mink4|4"];
            147:432	-> 146:431	 [id=234  color="gray"];
            148:433	-> 146:430	 [id=235  color="gray"];
            149:434	-> 146:429	 [id=233  color="gray"];
            150:435	-> 146:426	 [id=230  color="gray"];
            150:436	-> 0:53	 [id=30 dir=none  color="gray" label="mink4|1"];
            151:437	-> 146:425	 [id=232  color="gray"];
            151:438	-> 154:444	 [id=228 dir=none  color="gray" label="mink4|26"];
            152:439	-> 146:424	 [id=229  color="gray"];
            152:440	-> 0:30	 [id=7 dir=none  color="gray" label="mink4|0"];
            153:441	-> 146:423	 [id=231  color="gray"];
            153:442	-> 0:76	 [id=53 dir=none  color="gray" label="mink4|4"];
            154:443	-> 146:422	 [id=227  color="gray"];
            155:445	-> 146:421	 [id=76  color="gray"];
            157:446	-> 0:6	 [id=247  color="gray"];
            123:453	-> 189:556	 [id=248  color="gray"];
            134:454	-> 113:402	 [id=280 dir=none  color="gray" label="mink4|5"];
            158:458	-> 157:457	 [id=245  color="gray"];
            159:459	-> 157:456	 [id=246  color="gray"];
            160:460	-> 157:455	 [id=244  color="gray"];
            161:461	-> 157:452	 [id=243  color="gray"];
            161:462	-> 0:29	 [id=6 dir=none  color="gray" label="mink4|0"];
            162:463	-> 157:451	 [id=241  color="gray"];
            162:464	-> 165:470	 [id=239 dir=none  color="gray" label="mink4|26"];
            163:465	-> 157:450	 [id=240  color="gray"];
            163:466	-> 0:52	 [id=29 dir=none  color="gray" label="mink4|1"];
            164:467	-> 157:449	 [id=242  color="gray"];
            164:468	-> 0:75	 [id=52 dir=none  color="gray" label="mink4|4"];
            165:469	-> 157:448	 [id=238  color="gray"];
            166:471	-> 157:447	 [id=75  color="gray"];
            168:472	-> 0:5	 [id=258  color="gray"];
            134:479	-> 189:557	 [id=259  color="gray"];
            145:480	-> 103:378	 [id=284 dir=none  color="gray" label="mink4|26"];
            169:484	-> 168:483	 [id=256  color="gray"];
            170:485	-> 168:482	 [id=257  color="gray"];
            171:486	-> 168:481	 [id=255  color="gray"];
            172:487	-> 168:478	 [id=252  color="gray"];
            172:488	-> 0:51	 [id=28 dir=none  color="gray" label="mink4|1"];
            173:489	-> 168:477	 [id=254  color="gray"];
            173:490	-> 176:496	 [id=250 dir=none  color="gray" label="mink4|26"];
            174:491	-> 168:476	 [id=251  color="gray"];
            174:492	-> 0:28	 [id=5 dir=none  color="gray" label="mink4|0"];
            175:493	-> 168:475	 [id=253  color="gray"];
            175:494	-> 0:74	 [id=51 dir=none  color="gray" label="mink4|4"];
            176:495	-> 168:474	 [id=249  color="gray"];
            177:497	-> 168:473	 [id=74  color="gray"];
            178:498	-> 0:4	 [id=268  color="gray"];
            145:504	-> 189:558	 [id=226  color="gray"];
            156:505	-> 0:48	 [id=25 dir=none  color="gray" label="mink4|1"];
            180:509	-> 178:508	 [id=266  color="gray"];
            181:510	-> 178:507	 [id=267  color="gray"];
            182:511	-> 178:506	 [id=265  color="gray"];
            183:512	-> 178:503	 [id=262  color="gray"];
            183:513	-> 0:50	 [id=27 dir=none  color="gray" label="mink4|1"];
            184:514	-> 178:502	 [id=264  color="gray"];
            184:515	-> 187:521	 [id=260 dir=none  color="gray" label="mink4|26"];
            185:516	-> 178:501	 [id=261  color="gray"];
            185:517	-> 0:27	 [id=4 dir=none  color="gray" label="mink4|0"];
            186:518	-> 178:500	 [id=263  color="gray"];
            186:519	-> 0:73	 [id=50 dir=none  color="gray" label="mink4|4"];
            187:520	-> 178:499	 [id=73  color="gray"];
            179:522	-> 0:3	 [id=281  color="gray"];
            156:530	-> 0:25	 [id=2 dir=none  color="gray" label="mink4|0"];
            156:531	-> 189:559	 [id=282  color="gray"];
            190:535	-> 179:534	 [id=278  color="gray"];
            191:536	-> 179:533	 [id=279  color="gray"];
            192:537	-> 179:532	 [id=277  color="gray"];
            193:538	-> 179:529	 [id=269  color="gray"];
            193:539	-> 0:26	 [id=3 dir=none  color="gray" label="mink4|0"];
            193:540	-> 0:49	 [id=26 dir=none  color="gray" label="mink4|1"];
            194:541	-> 179:528	 [id=276  color="gray"];
            194:542	-> 197:548	 [id=272 dir=none  color="gray" label="mink4|5"];
            195:543	-> 179:527	 [id=273  color="gray"];
            195:544	-> 198:550	 [id=271 dir=none  color="gray" label="mink4|26"];
            196:545	-> 179:526	 [id=275  color="gray"];
            196:546	-> 0:72	 [id=49 dir=none  color="gray" label="mink4|4"];
            197:547	-> 179:525	 [id=274  color="gray"];
            198:549	-> 179:524	 [id=270  color="gray"];
            199:551	-> 179:523	 [id=72  color="gray"];
            189:552	-> 0:2	 [id=186  color="gray"];
            167:560	-> 189:562	 [id=215  color="gray"];
            200:561	-> 189:563	 [id=202  color="gray"];
            188:565	-> 189:564	 [id=211  color="gray"];
        }
    )
    .unwrap();

    // graph.iter_crown(NodeIndex(201)).for_each(|h| {
    //     println!("Hedge: {:?}", h);
    // });

    // return;
    println!(
        "{}",
        graph.dot_of(&graph.compass_subgraph::<SuBitGraph>(Some(CompassPt::S)))
    );

    let a = graph.clone().extract_nodes(
        [1, 9, 20, 29, 38, 47, 56, 65, 74, 83, 84]
            .into_iter()
            .map(NodeIndex),
        |a| a.map(Clone::clone),
        |a| a,
    );

    // let mut out = String::new();
    // a.dot_serialize_fmt(&mut out, (), &|a| a.clone(), &|a| a.clone(), &|a| a.clone())
    //     .unwrap();
    println!("{}", a.base_dot());

    let sub: SuBitGraph = graph.compass_subgraph::<SuBitGraph>(Some(CompassPt::S));
    let node: NodeIndex = graph.len();
    println!("{node}");
    let a = graph.extract(
        &sub,
        |a| a.map(Clone::clone),
        |a| a,
        |a| {
            println!("{a}");
            a.clone()
        },
        |a| a,
    );

    let mut out = String::new();
    a.dot_serialize_fmt(&mut out, (), &|a| a.clone(), &|a| a.clone(), &|a| a.clone())
        .unwrap();
    println!("{out}");
}

#[test]
fn extract_nodes() {
    let graph: DotGraph = dot!(digraph{
        A->B
        C->D
        B->C
        D->A
    })
    .unwrap();

    let a = graph.clone().extract_nodes(
        [1].into_iter().map(NodeIndex),
        |a| a.map(Clone::clone),
        |a| a,
    );

    println!("{a:#?}");

    // let mut out = String::new();
    // a.dot_serialize_fmt(&mut out, (), &|a| a.clone(), &|a| a.clone(), &|a| a.clone())
    //     .unwrap();
    println!("{}", a.base_dot());
}
