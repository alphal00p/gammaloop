use crate::{
    dot,
    half_edge::{
        builder::HedgeGraphBuilder,
        involution::{Flow, Hedge},
        nodestore::NodeStorageOps,
        subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike},
        HedgeGraph, NodeIndex,
    },
    parser::{DotGraph, DotVertexData},
    tree::{child_vec::ChildVecStore, Forest},
};

#[test]
fn extract_forest() {
    let mut aligned: DotGraph<Forest<DotVertexData, ChildVecStore<()>>> = dot!(
    digraph {
      ext4 [flow=sink];
      0 -> 1;
      2-> ext4;
      0 -> 2;
      0 -> 3[dir=none];
      1 -> 2;
      1 -> 1;
      1 -> 3;
      2 -> 3;
    })
    .unwrap();

    let mut subgraph: SuBitGraph = aligned.empty_subgraph();
    subgraph.add(Hedge(0));
    subgraph.add(Hedge(7));
    subgraph.add(Hedge(3));
    subgraph.add(Hedge(2));

    // for n in aligned.node_store.iter_node_ids() {
    //     aligned.node_store.nodes.set_node_data((), n);
    // }

    // println!("{}", aligned.node_store.nodes.debug_draw(|_| None));
    // aligned.node_store.swap(Hedge(1), Hedge(4));

    // println!("{}", aligned.node_store.nodes.debug_draw(|_| None));
    // aligned.node_store.swap(Hedge(1), Hedge(3));

    // println!("{}", aligned.node_store.nodes.debug_draw(|_| None));

    println!("{}", aligned.dot(&subgraph));

    // aligned.identify_nodes(&[NodeIndex(1), NodeIndex(2)], DotVertexData::empty());

    aligned.forget_identification_history();
    println!("{}", aligned.dot(&subgraph));

    let extracted = aligned.extract(
        &subgraph,
        |a| a.map(Clone::clone),
        |a| a,
        |a| a.clone(),
        |a| a,
    );
    // println!("{:?}", aligned.node_store.node_len());

    // println!("{:?}", extracted.node_store.node_len());

    println!("{}", extracted.base_dot());
    println!("{}", aligned.base_dot());
}

#[test]
fn extact_single_dangling() {
    let mut simple = HedgeGraphBuilder::new();
    let n1 = simple.add_node(());
    let n2 = simple.add_node(());
    simple.add_edge(n1, n2, (), false);
    simple.add_external_edge(n1, (), false, Flow::Sink);
    simple.add_edge(n1, n1, (), false);
    simple.add_external_edge(n2, (), false, Flow::Sink);
    let mut simple: HedgeGraph<(), (), ()> = simple.build();

    let mut single_hair: SuBitGraph = simple.empty_subgraph();
    if let Some(s) = simple.iter_edges().find(|a| a.0.is_unpaired()) {
        single_hair.add(s.0);
    }

    simple.extract(&single_hair, |a| a.map(Clone::clone), |a| a, |a| *a, |a| a);
}

#[test]
fn forest_storage_supports_empty_nodes() {
    type Store = Forest<(), ChildVecStore<()>>;

    let mut builder = HedgeGraphBuilder::<(), (), ()>::new();
    let n1 = builder.add_node(());
    let n2 = builder.add_node(());
    let empty = builder.add_node(());
    builder.add_edge(n1, n2, (), false);

    let graph: HedgeGraph<(), (), (), Store> = builder.build();
    assert_eq!(graph.iter_crown(empty).count(), 0);
    graph.node_store.check_nodes().unwrap();

    let (_, graph) = graph
        .add_dangling_edge(empty, (), Flow::Sink, false)
        .unwrap();
    assert_eq!(graph.iter_crown(empty).count(), 1);
    graph.node_store.check_nodes().unwrap();
}

#[test]
fn extract_buggy() {
    let mut aligned: DotGraph = dot!(
    digraph {
        ext0 [flow=sink];
        ext0 -> 0[dir=back];
        ext3 [flow=sink];
        ext3 -> 0[dir=none];
        ext6 [flow=sink];
        ext6 -> 0[dir=none];
        ext9 [flow=sink];
        ext9 -> 0[dir=none];
        1 -> 0[ dir=forward];
        2 -> 1[ dir=forward];
        2 -> 0[dir=none];
        2 -> 0[dir=none];
        3 -> 1[dir=forward];
        3 -> 0[dir=none];
        4 -> 1[ dir=forward];
        5 -> 0[ dir=forward];
        5 -> 0[ dir=none];
        5 -> 0[ dir=none];
        5 -> 0[ dir=none];
    })
    .unwrap();

    let mut subgraph: SuBitGraph = aligned.empty_subgraph();

    let nodes = [1, 2, 3, 4];

    for n in nodes {
        for h in aligned.iter_crown(NodeIndex(n)) {
            subgraph.add(h)
        }
    }
    // subgraph.add(Hedge(0));
    // subgraph.add(Hedge(7));
    // subgraph.add(Hedge(3));
    // subgraph.add(Hedge(2));

    // for n in aligned.node_store.iter_node_ids() {
    //     aligned.node_store.nodes.set_node_data((), n);
    // }

    // println!("{}", aligned.node_store.nodes.debug_draw(|_| None));
    // aligned.node_store.swap(Hedge(1), Hedge(4));

    // println!("{}", aligned.node_store.nodes.debug_draw(|_| None));
    // aligned.node_store.swap(Hedge(1), Hedge(3));

    // println!("{}", aligned.node_store.nodes.debug_draw(|_| None));

    println!("{}", aligned.dot(&subgraph));

    // aligned.identify_nodes(&[NodeIndex(1), NodeIndex(2)], DotVertexData::empty());

    // aligned.forget_identification_history();
    // println!("{}", aligned.dot(&subgraph));

    let extracted = aligned.extract(
        &subgraph,
        |a| a.map(Clone::clone),
        |a| a,
        |a| a.clone(),
        |a| a,
    );

    aligned.node_store.check_and_set_nodes().unwrap();
    // println!("{:?}", aligned.node_store.node_len());

    // println!("{:?}", extracted.node_store.node_len());

    println!("{}", extracted.base_dot());
    println!("{}", aligned.base_dot());
}

#[test]
fn extract_normal() {
    let mut aligned: DotGraph = dot!(
    digraph {
      ext4 [flow=sink];
      0 -> 1;
      2-> ext4;
      0 -> 2;
      0 -> 3[dir=none];
      1 -> 2;
      1 -> 1;
      1 -> 3;
      2 -> 3;
    })
    .unwrap();

    let mut subgraph: SuBitGraph = aligned.empty_subgraph();
    subgraph.add(Hedge(0));
    subgraph.add(Hedge(7));
    subgraph.add(Hedge(3));
    subgraph.add(Hedge(2));

    // for n in aligned.node_store.iter_node_ids() {
    //     aligned.node_store.nodes.set_node_data((), n);
    // }

    // println!("{}", aligned.node_store.nodes.debug_draw(|_| None));
    // aligned.node_store.swap(Hedge(1), Hedge(4));

    // println!("{}", aligned.node_store.nodes.debug_draw(|_| None));
    // aligned.node_store.swap(Hedge(1), Hedge(3));

    // println!("{}", aligned.node_store.nodes.debug_draw(|_| None));

    println!("{}", aligned.dot(&subgraph));

    let (_, s): (_, SuBitGraph) = aligned
        .identify_nodes_without_self_edges(&[NodeIndex(1), NodeIndex(2)], DotVertexData::empty());

    aligned.forget_identification_history();
    aligned.delete_hedges(&s);
    aligned.node_store.check_and_set_nodes().unwrap();

    println!("{}", aligned.dot(&subgraph));

    let extracted = aligned.extract(
        &subgraph,
        |a| a.map(Clone::clone),
        |a| a,
        |a| a.clone(),
        |a| a,
    );

    aligned.node_store.check_and_set_nodes().unwrap();
    // println!("{:?}", aligned.node_store.node_len());

    // println!("{:?}", extracted.node_store.node_len());

    println!("{}", extracted.base_dot());
    println!("{}", aligned.base_dot());
}

#[test]
fn identify_nodes_of_subgraph_marks_same_self_edges() {
    let aligned: DotGraph = dot!(
    digraph {
      ext4 [flow=sink];
      0 -> 1;
      2-> ext4;
      0 -> 2;
      0 -> 3[dir=none];
      1 -> 2;
      1 -> 1;
      1 -> 3;
      2 -> 3;
    })
    .unwrap();

    let nodes = [NodeIndex(1), NodeIndex(2)];
    let mut subgraph: SuBitGraph = aligned.empty_subgraph();
    for node in nodes {
        for hedge in aligned.iter_crown(node) {
            subgraph.add(hedge);
        }
    }

    let mut old = aligned.clone();
    let (_, old_self_edges): (_, SuBitGraph) =
        old.identify_nodes_without_self_edges(&nodes, DotVertexData::empty());

    let mut new = aligned.clone();
    let (_, new_self_edges) = new
        .identify_nodes_of_subgraph_without_self_edges::<_, SuBitGraph>(
            &subgraph,
            DotVertexData::empty(),
        )
        .unwrap();
    assert_eq!(
        old_self_edges.included_iter().collect::<Vec<_>>(),
        new_self_edges.included_iter().collect::<Vec<_>>()
    );

    let mut marked = aligned;
    let mut marked_self_edges: SuBitGraph = marked.empty_subgraph();
    marked
        .identify_nodes_of_subgraph_marking_self_edges(
            subgraph.included_iter(),
            |hedge| subgraph.includes(&hedge),
            DotVertexData::empty(),
            &mut marked_self_edges,
        )
        .unwrap();
    assert_eq!(
        old_self_edges.included_iter().collect::<Vec<_>>(),
        marked_self_edges.included_iter().collect::<Vec<_>>()
    );
}

#[test]
fn identify_compact_subgraph_preserves_hidden_edges_and_root_order() {
    let mut builder = HedgeGraphBuilder::<(), (), ()>::new();
    let n0 = builder.add_node(());
    let n1 = builder.add_node(());
    builder.add_edge(n1, n0, (), false);
    builder.add_edge(n0, n1, (), false);
    builder.add_edge(n1, n1, (), false);
    builder.add_external_edge(n0, (), false, Flow::Sink);
    let graph: HedgeGraph<(), (), (), Forest<(), ChildVecStore<()>>> = builder.build();

    for (visible, expected_marks, expected_root) in [
        (vec![0, 1, 2, 3, 4, 5, 6], vec![0, 1, 2, 3, 6], Some(n1)),
        (vec![0, 1, 2, 4, 5, 6], vec![0, 1, 6], Some(n1)),
        (vec![2, 3, 4, 5, 6], vec![2, 3, 6], Some(n0)),
        (vec![], vec![6], None),
    ] {
        let hedges = visible.into_iter().map(Hedge).collect::<Vec<_>>();
        let mut subgraph: SuBitGraph = graph.empty_subgraph();
        for &hedge in &hedges {
            subgraph.add(hedge);
        }

        let mut dense = graph.clone();
        let dense_result =
            dense.identify_nodes_of_subgraph_without_self_edges::<_, SuBitGraph>(&subgraph, ());
        assert_eq!(dense_result.as_ref().map(|(root, _)| *root), expected_root);

        let mut compact = graph.clone();
        let mut marked: SuBitGraph = compact.empty_subgraph();
        // Existing marks survive; hidden mates, existing loops and dangling
        // edges must not acquire new marks when the visible island collapses.
        marked.add(Hedge(6));
        let compact_root = compact.identify_nodes_of_subgraph_marking_self_edges(
            hedges.iter().copied(),
            |hedge| hedges.binary_search(&hedge).is_ok(),
            (),
            &mut marked,
        );
        assert_eq!(compact_root, expected_root);
        assert_eq!(
            marked.included_iter().collect::<Vec<_>>(),
            expected_marks.into_iter().map(Hedge).collect::<Vec<_>>()
        );
        if let Some((_, mut dense_marks)) = dense_result {
            dense_marks.add(Hedge(6));
            assert_eq!(
                marked.included_iter().collect::<Vec<_>>(),
                dense_marks.included_iter().collect::<Vec<_>>()
            );
        }
        for hedge in 0..7 {
            assert_eq!(dense.node_id(Hedge(hedge)), compact.node_id(Hedge(hedge)));
        }
    }
}

#[test]
fn orientation_hedges() {
    let mut single_node = HedgeGraphBuilder::new();
    let a = single_node.add_node(());
    single_node.add_external_edge(a, (), true, Flow::Source);
    single_node.add_external_edge(a, (), true, Flow::Sink);
    let aligned: HedgeGraph<(), (), ()> = single_node.build();

    println!("{}", aligned.base_dot())
}
