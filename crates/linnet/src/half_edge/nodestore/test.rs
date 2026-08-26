use crate::{
    dot,
    half_edge::{
        builder::HedgeGraphBuilder,
        involution::{Flow, Hedge},
        nodestore::{NodeStorageOps, NodeStorageVec},
        subgraph::{ModifySubSet, SuBitGraph, SubSetLike},
        HedgeGraph, HedgeGraphError, NodeIndex,
    },
    parser::{DotGraph, DotVertexData},
    tree::{child_vec::ChildVecStore, Forest},
};

#[test]
fn convert_node_store_preserves_graph_storage_and_isolated_nodes() {
    type ForestStore = Forest<&'static str, ChildVecStore<()>>;

    let mut builder = HedgeGraphBuilder::<&'static str, &'static str, &'static str>::new();
    let source = builder.add_node("source node");
    let sink = builder.add_node("sink node");
    let isolated = builder.add_node("isolated node");
    builder.add_edge(
        source.add_data("source half-edge"),
        sink.add_data("sink half-edge"),
        "internal edge",
        true,
    );
    builder.add_external_edge(
        source.add_data("external half-edge"),
        "external edge",
        false,
        Flow::Sink,
    );
    let graph: HedgeGraph<_, _, _, NodeStorageVec<_>> = builder.build();
    let before = graph
        .iter_hedges()
        .map(|(hedge, data)| {
            (
                hedge,
                *data,
                graph.node_id(hedge),
                graph.inv(hedge),
                graph.flow(hedge),
                graph[&hedge],
                graph[graph[&hedge]],
            )
        })
        .collect::<Vec<_>>();

    let forest: HedgeGraph<_, _, _, ForestStore> = graph.into_node_store().unwrap();
    let after = forest
        .iter_hedges()
        .map(|(hedge, data)| {
            (
                hedge,
                *data,
                forest.node_id(hedge),
                forest.inv(hedge),
                forest.flow(hedge),
                forest[&hedge],
                forest[forest[&hedge]],
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(after, before);
    assert_eq!(forest[source], "source node");
    assert_eq!(forest[sink], "sink node");
    assert_eq!(forest[isolated], "isolated node");
    assert_eq!(forest.iter_crown(isolated).count(), 0);

    let roundtrip: HedgeGraph<_, _, _, NodeStorageVec<_>> = forest.into_node_store().unwrap();
    assert_eq!(roundtrip.n_nodes(), 3);
    assert_eq!(roundtrip[source], "source node");
    assert_eq!(roundtrip[sink], "sink node");
    assert_eq!(roundtrip[isolated], "isolated node");
    assert_eq!(roundtrip.iter_crown(isolated).count(), 0);
    assert_eq!(roundtrip[Hedge(0)], "source half-edge");
    assert_eq!(roundtrip[Hedge(1)], "sink half-edge");
    assert_eq!(roundtrip[Hedge(2)], "external half-edge");
}

#[test]
fn converting_forest_history_requires_explicit_compaction() {
    type ForestStore = Forest<&'static str, ChildVecStore<()>>;

    let mut builder = HedgeGraphBuilder::<(), &'static str, ()>::new();
    let a = builder.add_node("a");
    let b = builder.add_node("b");
    let c = builder.add_node("c");
    builder.add_edge(a, b, (), false);
    builder.add_edge(b, c, (), false);
    let graph: HedgeGraph<_, _, _, ForestStore> = builder.build();

    let mut retained = graph.clone();
    retained.identify_nodes(&[a, b], "ab");
    let conversion: Result<HedgeGraph<_, _, _, NodeStorageVec<_>>, _> = retained.into_node_store();
    assert!(matches!(
        conversion,
        Err(HedgeGraphError::NodesDoNotPartition(_))
    ));

    let mut compacted = graph;
    compacted.identify_nodes(&[a, b], "ab");
    let discarded = compacted
        .forget_identification_history()
        .into_iter()
        .map(|(_, data)| data)
        .collect::<Vec<_>>();
    assert_eq!(discarded, vec!["b"]);

    let converted: HedgeGraph<_, _, _, NodeStorageVec<_>> = compacted.into_node_store().unwrap();
    let merged = converted.node_id(Hedge(0));
    assert_eq!(converted[merged], "ab");
    assert_eq!(converted.node_id(Hedge(1)), merged);
    assert_eq!(converted.node_id(Hedge(2)), merged);
    assert_eq!(converted[converted.node_id(Hedge(3))], "c");
}

#[test]
fn conversion_rejects_identified_isolated_forest_root() {
    type ForestStore = Forest<&'static str, ChildVecStore<()>>;

    let mut builder = HedgeGraphBuilder::<(), &'static str, ()>::new();
    let a = builder.add_node("a");
    let b = builder.add_node("b");
    let graph: HedgeGraph<_, _, _, ForestStore> = builder.build();

    let mut retained = graph.clone();
    retained.identify_nodes(&[a, b], "ab");
    let conversion: Result<HedgeGraph<_, _, _, NodeStorageVec<_>>, _> = retained.into_node_store();
    assert!(matches!(
        conversion,
        Err(HedgeGraphError::NodesDoNotPartition(_))
    ));

    let mut compacted = graph;
    compacted.identify_nodes(&[a, b], "ab");
    let discarded = compacted
        .forget_identification_history()
        .into_iter()
        .map(|(_, data)| data)
        .collect::<Vec<_>>();
    assert_eq!(discarded, vec!["b"]);

    let converted: HedgeGraph<_, _, _, NodeStorageVec<_>> = compacted.into_node_store().unwrap();
    assert_eq!(converted.n_nodes(), 1);
    assert_eq!(converted[NodeIndex(0)], "ab");
    assert_eq!(converted.iter_crown(NodeIndex(0)).count(), 0);
}

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
fn forest_extract_nodes_moves_selected_empty_nodes() {
    type Store = Forest<&'static str, ChildVecStore<()>>;

    let mut builder = HedgeGraphBuilder::<(), &'static str, ()>::new();
    let source = builder.add_node("source");
    let isolated = builder.add_node("isolated");
    let sink = builder.add_node("sink");
    builder.add_edge(source, sink, (), false);
    let mut graph: HedgeGraph<_, _, _, Store> = builder.build();

    let extracted = graph.extract_nodes([isolated], |edge| edge.map(|_| ()), |edge| edge);

    assert_eq!(graph.n_nodes(), 2);
    assert_eq!(graph.n_hedges(), 2);
    assert_eq!(graph[graph.node_id(Hedge(0))], "source");
    assert_eq!(graph[graph.node_id(Hedge(1))], "sink");
    graph.node_store.check_nodes().unwrap();

    assert_eq!(extracted.n_nodes(), 1);
    assert_eq!(extracted.n_hedges(), 0);
    assert_eq!(extracted[NodeIndex(0)], "isolated");
    assert_eq!(extracted.iter_crown(NodeIndex(0)).count(), 0);
    extracted.node_store.check_nodes().unwrap();
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
            &subgraph,
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
fn orientation_hedges() {
    let mut single_node = HedgeGraphBuilder::new();
    let a = single_node.add_node(());
    single_node.add_external_edge(a, (), true, Flow::Source);
    single_node.add_external_edge(a, (), true, Flow::Sink);
    let aligned: HedgeGraph<(), (), ()> = single_node.build();

    println!("{}", aligned.base_dot())
}
