use crate::{
    dot,
    half_edge::{
        builder::HedgeGraphBuilder,
        involution::{Flow, Hedge},
        nodestore::{NodeForest, NodeStorageOps, NodeStorageVec},
        subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike},
        swap::Swap,
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

#[test]
fn both_stores_report_absorbed_nodes() {
    fn check<N: NodeStorageOps<NodeData = ()>>() {
        let mut builder = HedgeGraphBuilder::<(), ()>::new();
        let [a, b, c, d] = [(); 4].map(|_| builder.add_node(()));
        builder.add_edge(a, b, (), true);
        let mut graph: HedgeGraph<(), (), _, N> = builder.build();
        let ab = graph.identify_nodes(&[a, b], ());
        let cd = graph.identify_nodes(&[c, d], ());

        let absorbed: Vec<_> = (0..graph.n_nodes())
            .map(NodeIndex)
            .filter(|&n| graph.node_store.is_absorbed(n))
            .collect();
        assert_eq!(absorbed.len(), 2);
        assert!(!absorbed.contains(&ab) && !absorbed.contains(&cd));
        // The absorbed isolated node is not an extra component.
        let components = graph.strongly_connected_components(
            crate::half_edge::algorithms::DirectionBasis::Underlying,
        );
        assert_eq!(components, [vec![ab], vec![cd]]);
    }
    check::<NodeStorageVec<()>>();
    check::<Forest<(), ChildVecStore<()>>>();
}

/// Drives every node store through the same random sequence of graph edits and
/// compares what does not depend on the store: node indices and survivors
/// differ, but partitions, payloads, crowns and counts must agree.
mod differential {
    use std::fmt::Debug;

    use crate::{
        half_edge::{
            builder::HedgeGraphBuilder,
            involution::{Flow, Hedge, Orientation},
            nodestore::{NodeForest, NodeStorageOps, NodeStorageVec},
            subgraph::{ModifySubSet, SuBitGraph, SubSetLike},
            HedgeGraph, NodeIndex,
        },
        tree::{child_pointer::ParentChildStore, child_vec::ChildVecStore, Forest},
    };

    type Graph<S> = HedgeGraph<usize, usize, usize, S>;

    trait Store: NodeStorageOps<NodeData = usize, Base = SuBitGraph> + Clone {}
    impl<S: NodeStorageOps<NodeData = usize, Base = SuBitGraph> + Clone> Store for S {}

    struct Random(u64);

    impl Random {
        fn below(&mut self, bound: usize) -> usize {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            (self.0 % bound.max(1) as u64) as usize
        }

        fn one_in(&mut self, n: usize) -> bool {
            self.below(n) == 0
        }
    }

    /// A random graph with unique node, edge and half-edge payloads from `next`.
    fn random_graph(
        random: &mut Random,
        next: &mut usize,
    ) -> HedgeGraphBuilder<usize, usize, usize> {
        let mut fresh = || {
            *next += 1;
            *next
        };
        let mut builder = HedgeGraphBuilder::new();
        let nodes: Vec<_> = (0..1 + random.below(7))
            .map(|_| builder.add_node(fresh()))
            .collect();
        for _ in 0..random.below(10) {
            let (a, b) = (
                nodes[random.below(nodes.len())],
                nodes[random.below(nodes.len())],
            );
            builder.add_edge(
                a.add_data(fresh()),
                b.add_data(fresh()),
                fresh(),
                random.one_in(2),
            );
        }
        for _ in 0..random.below(3) {
            let flow = if random.one_in(2) {
                Flow::Sink
            } else {
                Flow::Source
            };
            let node = nodes[random.below(nodes.len())];
            builder.add_external_edge(node.add_data(fresh()), fresh(), true, flow);
        }
        builder
    }

    /// One store's graph, with its index of every logical node and the node
    /// that each index was merged into by pending identifications.
    struct Lane<S: Store> {
        graph: Graph<S>,
        slots: Vec<NodeIndex>,
        merged_into: Vec<usize>,
        /// Whether `identify_nodes` accepts absorbed nodes.
        resolves_absorbed: bool,
    }

    /// What every store must agree on.
    #[derive(Debug, PartialEq)]
    struct Observation {
        hedge_payloads: Vec<usize>,
        involution: Vec<usize>,
        partition: Vec<usize>,
        node_payloads: Vec<usize>,
        crowns: Vec<Vec<usize>>,
        live_payloads: Vec<usize>,
        isolated_payloads: Vec<usize>,
    }

    /// Checks a store's consistency and observes its live nodes.
    ///
    /// Forest stores keep the subtree an absorbed root owned until history is
    /// forgotten, so their absorbed crowns overlap live ones and `check_nodes`
    /// holds for every store only without pending identifications.
    fn observe<S: Store>(graph: &Graph<S>, context: &str) -> Observation {
        graph.check().unwrap();
        let store = &graph.node_store;
        let pending = graph
            .iter_nodes()
            .any(|(node, _, _)| store.is_absorbed(node));
        if let Err(error) = store.check_nodes().map_err(|e| e.to_string()) {
            assert!(
                pending,
                "{context}: {error} in {}",
                std::any::type_name::<S>()
            );
        }
        let n_hedges = graph.n_hedges();
        let mut first_seen = vec![usize::MAX; graph.n_nodes()];
        let mut partition = Vec::with_capacity(n_hedges);
        let mut covered = 0;
        let (mut live_payloads, mut isolated_payloads) = (vec![], vec![]);
        for (node, crown, &payload) in graph.iter_nodes() {
            if store.is_absorbed(node) {
                continue;
            }
            let crown_len = crown.len();
            assert_eq!(crown_len, crown.clone().count());
            covered += crown_len;
            for hedge in crown {
                assert_eq!(graph.node_id(hedge), node);
            }
            live_payloads.push(payload);
            if crown_len == 0 {
                isolated_payloads.push(payload);
            }
        }
        assert_eq!(
            covered, n_hedges,
            "{context}: live crowns cover the half-edges"
        );
        for hedge in 0..n_hedges {
            let node = graph.node_id(Hedge(hedge)).0;
            assert!(!graph.node_store.is_absorbed(NodeIndex(node)));
            if first_seen[node] == usize::MAX {
                first_seen[node] = partition.iter().max().map_or(0, |&m: &usize| m + 1);
            }
            partition.push(first_seen[node]);
        }
        live_payloads.sort();
        isolated_payloads.sort();
        let crowns = (0..n_hedges)
            .map(|hedge| {
                let mut crown: Vec<_> = graph
                    .iter_crown(graph.node_id(Hedge(hedge)))
                    .map(|h| h.0)
                    .collect();
                crown.sort();
                crown
            })
            .collect();
        Observation {
            hedge_payloads: (0..n_hedges).map(|h| graph[Hedge(h)]).collect(),
            involution: (0..n_hedges).map(|h| graph.inv(Hedge(h)).0).collect(),
            node_payloads: (0..n_hedges)
                .map(|h| graph[graph.node_id(Hedge(h))])
                .collect(),
            partition,
            crowns,
            live_payloads,
            isolated_payloads,
        }
    }

    impl<S: Store> Lane<S> {
        fn new(builder: HedgeGraphBuilder<usize, usize, usize>, resolves_absorbed: bool) -> Self {
            let graph: Graph<S> = builder.build();
            let mut lane = Lane {
                graph,
                slots: vec![],
                merged_into: vec![],
                resolves_absorbed,
            };
            lane.renumbered();
            lane
        }

        fn live(&self, slot: usize) -> usize {
            let mut node = self.slots[slot].0;
            while self.merged_into[node] != node {
                node = self.merged_into[node];
            }
            node
        }

        /// Logical nodes follow payload order after a renumbering edit, which
        /// leaves no pending identification.
        fn renumbered(&mut self) {
            let mut nodes: Vec<_> = self.graph.iter_nodes().map(|(n, _, &p)| (p, n)).collect();
            nodes.sort();
            self.slots = nodes.into_iter().map(|(_, node)| node).collect();
            self.merged_into = (0..self.graph.n_nodes()).collect();
        }

        /// Appended nodes keep their order behind the existing ones.
        fn appended(&mut self) {
            let old = self.merged_into.len();
            self.slots
                .extend((old..self.graph.n_nodes()).map(NodeIndex));
            self.merged_into.extend(old..self.graph.n_nodes());
        }

        fn identify(&mut self, slots: &[usize], payload: usize) {
            let resolved: Vec<usize> = slots.iter().map(|&slot| self.live(slot)).collect();
            let requested: Vec<NodeIndex> = if self.resolves_absorbed {
                slots.iter().map(|&slot| self.slots[slot]).collect()
            } else {
                resolved.iter().copied().map(NodeIndex).collect()
            };
            let survivor = self.graph.identify_nodes(&requested, payload).0;
            assert!(resolved.contains(&survivor));
            for node in resolved {
                self.merged_into[node] = survivor;
                assert_eq!(
                    self.graph.node_store.is_absorbed(NodeIndex(node)),
                    node != survivor
                );
            }
        }

        /// Forgets and returns how many nodes were absorbed.
        fn forget(&mut self) -> usize {
            let removed = self
                .graph
                .forget_identification_history()
                .into_iter()
                .count();
            self.renumbered();
            removed
        }
    }

    struct Lanes {
        vec: Lane<NodeStorageVec<usize>>,
        parent_child: Lane<Forest<usize, ParentChildStore<()>>>,
        child_vec: Lane<Forest<usize, ChildVecStore<()>>>,
        ring: Lane<NodeForest<usize>>,
    }

    /// Runs `$body` on every lane and returns the four results.
    macro_rules! each {
        ($lanes:expr, |$lane:ident| $body:expr) => {
            [
                {
                    let $lane = &mut $lanes.vec;
                    $body
                },
                {
                    let $lane = &mut $lanes.parent_child;
                    $body
                },
                {
                    let $lane = &mut $lanes.child_vec;
                    $body
                },
                {
                    let $lane = &mut $lanes.ring;
                    $body
                },
            ]
        };
    }

    fn agree<T: PartialEq + Debug>(values: [T; 4], context: &str) -> T {
        let [vec, parent_child, child_vec, ring] = values;
        assert_eq!(vec, parent_child, "{context}: parent-child forest");
        assert_eq!(vec, child_vec, "{context}: child-vector forest");
        assert_eq!(vec, ring, "{context}: node forest");
        vec
    }

    fn random_subgraph(random: &mut Random, n_hedges: usize) -> SuBitGraph {
        let mut subgraph = SuBitGraph::empty(n_hedges);
        for hedge in 0..n_hedges {
            if random.one_in(3) {
                subgraph.add(Hedge(hedge));
            }
        }
        subgraph
    }

    impl Lanes {
        fn new(builder: HedgeGraphBuilder<usize, usize, usize>) -> Self {
            Lanes {
                vec: Lane::new(builder.clone(), false),
                parent_child: Lane::new(builder.clone(), true),
                child_vec: Lane::new(builder.clone(), true),
                ring: Lane::new(builder, true),
            }
        }

        fn n_slots(&self) -> usize {
            self.vec.slots.len()
        }

        fn forget(&mut self, context: &str) {
            // Survivors, and so the absorbed payloads, differ between stores.
            agree(each!(self, |lane| lane.forget()), context);
        }

        fn check(&mut self, context: &str) {
            agree(each!(self, |lane| observe(&lane.graph, context)), context);
            self.vec.graph.node_store.check_nodes().unwrap();
            // The node forest keeps absorbed crowns empty and live ones ascending
            // in every state.
            let ring = &self.ring.graph.node_store;
            ring.check_nodes().unwrap();
            let live = (0..ring.node_bound()).filter(|&n| !ring.is_absorbed(NodeIndex(n)));
            assert_eq!(ring.n_live(), live.count());
            for node in 0..ring.node_bound() {
                let crown: Vec<_> = ring.get_neighbor_iterator(NodeIndex(node)).collect();
                assert!(crown.is_sorted_by(|a, b| a < b), "{context}: {crown:?}");
                assert!(crown.is_empty() || !ring.is_absorbed(NodeIndex(node)));
            }
        }

        fn step(&mut self, random: &mut Random, next: &mut usize) -> String {
            let mut fresh = || {
                *next += 1;
                *next
            };
            let n_slots = self.n_slots();
            // An empty graph can only grow.
            let edit = if n_slots == 0 { 13 } else { random.below(16) };
            match edit {
                0..=4 => {
                    // Possibly repeated or already absorbed nodes.
                    let slots: Vec<_> = (0..1 + random.below(4))
                        .map(|_| random.below(n_slots))
                        .collect();
                    let payload = fresh();
                    each!(self, |lane| lane.identify(&slots, payload));
                    format!("identify {slots:?}")
                }
                5 => {
                    let slot = random.below(n_slots);
                    let (hedge, edge, flow) = (fresh(), fresh(), random.one_in(2));
                    let flow = if flow { Flow::Sink } else { Flow::Source };
                    each!(self, |lane| {
                        let node = NodeIndex(lane.live(slot));
                        let graph = lane.graph.clone();
                        lane.graph = graph
                            .add_dangling_edge(node.add_data(hedge), edge, flow, true)
                            .unwrap()
                            .1;
                    });
                    format!("dangling at {slot}")
                }
                6 | 7 => {
                    let (a, b) = (random.below(n_slots), random.below(n_slots));
                    let (ha, hb, edge) = (fresh(), fresh(), fresh());
                    each!(self, |lane| {
                        let (na, nb) = (NodeIndex(lane.live(a)), NodeIndex(lane.live(b)));
                        let graph = lane.graph.clone();
                        lane.graph = graph
                            .add_pair(na.add_data(ha), nb.add_data(hb), edge, Orientation::Default)
                            .unwrap()
                            .2;
                    });
                    format!("pair {a} {b}")
                }
                8 | 15 => {
                    self.forget("forget");
                    "forget".to_owned()
                }
                9 => {
                    self.forget("forget before delete");
                    let subgraph = random_subgraph(random, self.vec.graph.n_hedges());
                    each!(self, |lane| {
                        lane.graph.delete_hedges(&subgraph);
                        lane.renumbered();
                    });
                    format!("delete {subgraph:?}")
                }
                10 => {
                    self.forget("forget before extract");
                    let n_hedges = self.vec.graph.n_hedges();
                    if n_hedges == 0 {
                        return "extract nothing".to_owned();
                    }
                    // The bitset store cannot extract an empty selection.
                    let mut subgraph = random_subgraph(random, n_hedges);
                    subgraph.add(Hedge(random.below(n_hedges)));
                    let split_offset = 1_000_000;
                    let extracted = each!(self, |lane| {
                        let extracted = lane.graph.extract(
                            &subgraph,
                            |edge| edge.map(|data| *data),
                            |edge| edge,
                            |payload| payload + split_offset,
                            |payload| payload,
                        );
                        lane.renumbered();
                        observe(&extracted, "extracted")
                    });
                    agree(extracted, "extracted");
                    format!("extract {subgraph:?}")
                }
                11 => {
                    self.forget("forget before extract_nodes");
                    let slots: Vec<_> = (0..random.below(3))
                        .map(|_| random.below(self.n_slots()))
                        .collect();
                    let extracted = each!(self, |lane| {
                        let nodes: Vec<_> = slots.iter().map(|&slot| lane.slots[slot]).collect();
                        let extracted = lane.graph.extract_nodes(
                            nodes,
                            |edge| edge.map(|data| *data),
                            |edge| edge,
                        );
                        lane.renumbered();
                        observe(&extracted, "extracted")
                    });
                    agree(extracted, "extracted nodes");
                    format!("extract nodes {slots:?}")
                }
                12 => {
                    let n_hedges = self.vec.graph.n_hedges();
                    if n_hedges == 0 {
                        return "contract nothing".to_owned();
                    }
                    let mut subgraph = random_subgraph(random, n_hedges);
                    subgraph.add(Hedge(random.below(n_hedges)));
                    let payload = fresh();
                    each!(self, |lane| {
                        lane.graph.contract_subgraph(&subgraph, payload);
                        lane.renumbered();
                    });
                    format!("contract {subgraph:?}")
                }
                _ => {
                    let other = random_graph(random, next);
                    let join = random.one_in(2);
                    each!(self, |lane| {
                        let other: Graph<_> = other.clone().build();
                        if join {
                            let graph = lane.graph.clone();
                            lane.graph = graph
                                .join(
                                    other,
                                    |af, ad, bf, bd| af == -bf && ad.data % 2 == bd.data % 2,
                                    |af, ad, _, _| (af, ad),
                                )
                                .unwrap();
                        } else {
                            lane.graph.append_disconnected_mut(other).unwrap();
                        }
                        lane.appended();
                    });
                    format!("join {join}")
                }
            }
        }
    }

    #[test]
    fn node_stores_agree_on_random_edits() {
        for seed in 0..250u64 {
            let mut random = Random(0x9e37_79b9_7f4a_7c15 ^ (seed + 1).wrapping_mul(0x2545_f491));
            let mut next = 0;
            let mut lanes = Lanes::new(random_graph(&mut random, &mut next));
            lanes.check(&format!("seed {seed}: build"));
            let mut history = vec![];
            for _ in 0..24 {
                history.push(lanes.step(&mut random, &mut next));
                lanes.check(&format!("seed {seed}: {history:?}"));
            }
            lanes.forget(&format!("seed {seed}: final forget"));
            lanes.check(&format!("seed {seed}: final"));
        }
    }
}

#[test]
fn in_place_edge_insertion_preserves_graph_on_invalid_endpoints() {
    fn check<S>()
    where
        S: NodeStorageOps<NodeData = &'static str> + Clone + PartialEq + std::fmt::Debug,
    {
        let mut builder = HedgeGraphBuilder::<&str, &str, &str>::new();
        let a = builder.add_node("a");
        let b = builder.add_node("b");
        builder.add_edge(a.add_data("a port"), b.add_data("b port"), "existing", true);
        let mut graph: HedgeGraph<_, _, _, S> = builder.build();
        let before = graph.clone();
        let invalid = NodeIndex(graph.n_nodes());
        for (source, sink) in [(invalid, b), (a, invalid)] {
            assert!(matches!(
                graph.add_pair_mut(
                    source.add_data("new source"),
                    sink.add_data("new sink"),
                    "bad",
                    true
                ),
                Err(HedgeGraphError::NoNode)
            ));
            assert_eq!(graph, before);
        }
        assert!(matches!(
            graph.add_dangling_edge_mut(invalid.add_data("new port"), "bad", Flow::Sink, true),
            Err(HedgeGraphError::NoNode)
        ));
        assert_eq!(graph, before);
        graph.check().unwrap();
        graph.node_store.check_nodes().unwrap();
    }
    check::<NodeStorageVec<_>>();
    check::<Forest<_, ChildVecStore<()>>>();
    check::<Forest<_, crate::tree::child_pointer::ParentChildStore<()>>>();
    check::<NodeForest<_>>();
}

#[test]
fn in_place_edge_insertion_rejects_absorbed_endpoints_before_changing_the_graph() {
    fn check<S>()
    where
        S: NodeStorageOps<NodeData = &'static str> + Clone + PartialEq + std::fmt::Debug,
    {
        let mut builder = HedgeGraphBuilder::<&str, &str, &str>::new();
        let [a, b, c] = ["a", "b", "c"].map(|data| builder.add_node(data));
        builder.add_edge(a.add_data("a port"), b.add_data("b port"), "existing", true);
        let mut graph: HedgeGraph<_, _, _, S> = builder.build();
        let live = graph.identify_nodes(&[a, b], "ab");
        let absorbed = if live == a { b } else { a };
        assert!(graph.node_store.is_absorbed(absorbed));
        let before = graph.clone();
        for (source, sink) in [(c, absorbed), (absorbed, c), (live, absorbed)] {
            assert!(matches!(
                graph.add_pair_mut(
                    source.add_data("new source"),
                    sink.add_data("new sink"),
                    "bad",
                    true
                ),
                Err(HedgeGraphError::InvalidNode(node)) if node == absorbed
            ));
            assert_eq!(graph, before);
            assert_eq!(
                Swap::<Hedge>::len(&graph.node_store),
                Hedge(graph.n_hedges())
            );
        }
        assert!(matches!(
            graph.add_dangling_edge_mut(absorbed.add_data("new port"), "bad", Flow::Sink, true),
            Err(HedgeGraphError::InvalidNode(node)) if node == absorbed
        ));
        assert_eq!(graph, before);
        assert_eq!(
            Swap::<Hedge>::len(&graph.node_store),
            Hedge(graph.n_hedges())
        );
        graph.node_store.check_nodes().unwrap();
        // The rejected insertions leave a graph that still accepts live endpoints.
        let (source, sink) = graph
            .add_pair_mut(
                live.add_data("live source"),
                c.add_data("live sink"),
                "ok",
                true,
            )
            .unwrap();
        assert_eq!((graph.node_id(source), graph.node_id(sink)), (live, c));
        assert_eq!(
            Swap::<Hedge>::len(&graph.node_store),
            Hedge(graph.n_hedges())
        );
        graph.node_store.check_nodes().unwrap();
    }
    check::<NodeStorageVec<_>>();
    check::<Forest<_, ChildVecStore<()>>>();
    check::<Forest<_, crate::tree::child_pointer::ParentChildStore<()>>>();
    check::<NodeForest<_>>();
}

#[test]
fn in_place_edge_insertion_keeps_indices_and_incidence_payloads() {
    fn check<S: NodeStorageOps<NodeData = &'static str>>() {
        let mut builder = HedgeGraphBuilder::<&str, &str, &str>::new();
        let a = builder.add_node("a");
        let b = builder.add_node("b");
        let isolated = builder.add_node("isolated");
        builder.add_edge(a.add_data("a port"), b.add_data("b port"), "existing", true);
        let mut graph: HedgeGraph<_, _, _, S> = builder.build();
        let original_edge = graph[&Hedge(0)];
        let (source, sink) = graph
            .add_pair_mut(
                a.add_data("parallel source"),
                b.add_data("parallel sink"),
                "parallel",
                false,
            )
            .unwrap();
        let (loop_source, loop_sink) = graph
            .add_pair_mut(
                isolated.add_data("loop source"),
                isolated.add_data("loop sink"),
                "loop",
                true,
            )
            .unwrap();
        let dangling = graph
            .add_dangling_edge_mut(
                isolated.add_data("external port"),
                "external",
                Flow::Sink,
                true,
            )
            .unwrap();
        assert_eq!(graph.n_nodes(), 3);
        assert_eq!(graph.n_edges(), 4);
        assert_eq!(graph.n_hedges(), 7);
        assert_eq!(graph[&Hedge(0)], original_edge);
        assert_eq!(graph[original_edge], "existing");
        assert_eq!(graph.node_id(Hedge(0)), a);
        assert_eq!(graph.node_id(Hedge(1)), b);
        assert_eq!(graph.inv(source), sink);
        assert_eq!(graph.inv(loop_source), loop_sink);
        assert_eq!(graph.node_id(loop_source), isolated);
        assert_eq!(graph.node_id(loop_sink), isolated);
        assert_eq!(graph.inv(dangling), dangling);
        assert_eq!(graph.flow(dangling), Flow::Sink);
        assert_eq!(graph.iter_crown(isolated).count(), 3);
        let payloads: Vec<_> = graph.iter_hedges().map(|(_, data)| *data).collect();
        assert_eq!(
            payloads,
            [
                "a port",
                "b port",
                "parallel source",
                "parallel sink",
                "loop source",
                "loop sink",
                "external port"
            ]
        );
        graph.check().unwrap();
        graph.node_store.check_nodes().unwrap();
    }
    check::<NodeStorageVec<_>>();
    check::<Forest<_, ChildVecStore<()>>>();
    check::<Forest<_, crate::tree::child_pointer::ParentChildStore<()>>>();
    check::<NodeForest<_>>();
}
