use crate::{
    dot,
    half_edge::{
        builder::HedgeGraphBuilder,
        involution::{Flow, Hedge},
        nodestore::NodeStorageOps,
        subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike},
        swap::Swap,
        HedgeGraph, HedgeGraphError, NodeIndex,
    },
    parser::{DotGraph, DotVertexData},
};

#[test]
fn extract_forest() {
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
    let mut builder = HedgeGraphBuilder::<(), (), ()>::new();
    let n1 = builder.add_node(());
    let n2 = builder.add_node(());
    let empty = builder.add_node(());
    builder.add_edge(n1, n2, (), false);

    let graph: HedgeGraph<(), (), ()> = builder.build();
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
    let mut builder = HedgeGraphBuilder::<(), &'static str, ()>::new();
    let source = builder.add_node("source");
    let isolated = builder.add_node("isolated");
    let sink = builder.add_node("sink");
    builder.add_edge(source, sink, (), false);
    let mut graph: HedgeGraph<_, _, _> = builder.build();

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
    let graph: HedgeGraph<(), (), ()> = builder.build();

    for (visible, expected_marks, expected_root) in [
        (vec![0, 1, 2, 3, 4, 5, 6], vec![0, 1, 2, 3, 6], Some(n1)),
        (vec![0, 1, 2, 4, 5, 6], vec![0, 1, 6], Some(n1)),
        // The survivor is the node with the larger crown, not the first visited.
        (vec![2, 3, 4, 5, 6], vec![2, 3, 6], Some(n1)),
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

/// Contraction renumbers nodes but keeps every node it does not touch, isolated
/// ones included.
#[test]
fn contraction_keeps_unrelated_isolated_nodes() {
    let mut builder = HedgeGraphBuilder::new();
    let [a, b, c] = ["a", "b", "c"].map(|name| builder.add_node(name));
    builder.add_node("isolated d");
    builder.add_edge(a, b, (), true);
    builder.add_edge(b, c, (), true);
    let mut graph: HedgeGraph<(), &str> = builder.build();
    graph.add_node("isolated e");

    let mut edge: SuBitGraph = graph.empty_subgraph();
    edge.add(Hedge(0));
    edge.add(Hedge(1));
    graph.contract_subgraph(&edge, "ab");

    let mut data: Vec<_> = graph.iter_nodes().map(|(_, _, data)| *data).collect();
    data.sort();
    assert_eq!(data, ["ab", "c", "isolated d", "isolated e"]);
    graph.node_store.check_nodes().unwrap();
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
fn identification_reports_absorbed_nodes() {
    let mut builder = HedgeGraphBuilder::<(), ()>::new();
    let [a, b, c, d] = [(); 4].map(|_| builder.add_node(()));
    builder.add_edge(a, b, (), true);
    let mut graph: HedgeGraph<(), ()> = builder.build();
    let ab = graph.identify_nodes(&[a, b], ());
    let cd = graph.identify_nodes(&[c, d], ());

    let absorbed: Vec<_> = (0..graph.n_nodes())
        .map(NodeIndex)
        .filter(|&n| graph.node_store.is_absorbed(n))
        .collect();
    assert_eq!(absorbed.len(), 2);
    assert!(!absorbed.contains(&ab) && !absorbed.contains(&cd));
    // The absorbed isolated node is not an extra component.
    let components = graph
        .strongly_connected_components(crate::half_edge::algorithms::DirectionBasis::Underlying);
    assert_eq!(components, [vec![ab], vec![cd]]);
}

/// Drives the node store through random sequences of graph edits and checks
/// every state against a model of the node partition by payload: node indices
/// and survivors are the store's choice, but which half-edges share a node,
/// which nodes are live and how many nodes forgetting drops are not.
mod random_edits {
    use std::collections::{BTreeMap, BTreeSet};

    use crate::half_edge::{
        builder::HedgeGraphBuilder,
        involution::{Flow, Hedge, Orientation},
        nodestore::NodeStorageOps,
        subgraph::{ModifySubSet, SuBitGraph, SubSetLike},
        HedgeGraph, NodeIndex,
    };

    type Graph = HedgeGraph<usize, usize, usize>;

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

    fn random_subgraph(random: &mut Random, n_hedges: usize) -> SuBitGraph {
        let mut subgraph = SuBitGraph::empty(n_hedges);
        for hedge in 0..n_hedges {
            if random.one_in(3) {
                subgraph.add(Hedge(hedge));
            }
        }
        subgraph
    }

    /// The node partition by payload: the node of each half-edge and the live
    /// nodes, isolated ones included.
    #[derive(Debug, PartialEq)]
    struct Partition {
        node_of: BTreeMap<usize, usize>,
        live: BTreeSet<usize>,
    }

    impl Partition {
        /// Checks the consistency of `graph`, absorbed nodes and the ascending
        /// crowns included, and reads its partition.
        fn of(graph: &Graph, context: &str) -> Self {
            graph.check().unwrap();
            let store = &graph.node_store;
            if let Err(error) = store.check_nodes() {
                panic!("{context}: {error}");
            }
            let mut live = BTreeSet::new();
            let mut covered = 0;
            for (node, crown, &payload) in graph.iter_nodes() {
                let crown: Vec<_> = crown.collect();
                assert!(crown.is_sorted_by(|a, b| a < b), "{context}: {crown:?}");
                if store.is_absorbed(node) {
                    assert!(crown.is_empty(), "{context}: absorbed {node} has a crown");
                    continue;
                }
                covered += crown.len();
                for &hedge in &crown {
                    assert_eq!(graph.node_id(hedge), node, "{context}");
                }
                live.insert(payload);
            }
            assert_eq!(store.n_live(), live.len(), "{context}");
            assert_eq!(covered, graph.n_hedges(), "{context}: crowns cover");
            let node_of = (0..graph.n_hedges())
                .map(|hedge| (graph[Hedge(hedge)], graph[graph.node_id(Hedge(hedge))]))
                .collect();
            Partition { node_of, live }
        }

        /// Merges the live nodes `merged` into one carrying `payload`.
        fn identify(&mut self, merged: &BTreeSet<usize>, payload: usize) {
            for node in self.node_of.values_mut() {
                if merged.contains(node) {
                    *node = payload;
                }
            }
            self.live.retain(|node| !merged.contains(node));
            self.live.insert(payload);
        }

        /// Moves `hedges` out. Nodes that keep half-edges stay and leave a copy
        /// with payload `+ split_offset`; nodes that lose all of them leave.
        fn split_off(&mut self, hedges: &BTreeSet<usize>, split_offset: usize) -> Partition {
            let staying: BTreeSet<usize> = (self.node_of.iter())
                .filter(|(hedge, _)| !hedges.contains(hedge))
                .map(|(_, &node)| node)
                .collect();
            let mut moved = Partition {
                node_of: BTreeMap::new(),
                live: BTreeSet::new(),
            };
            for hedge in hedges {
                let node = self.node_of.remove(hedge).unwrap();
                let node = if staying.contains(&node) {
                    node + split_offset
                } else {
                    self.live.remove(&node);
                    node
                };
                moved.node_of.insert(*hedge, node);
                moved.live.insert(node);
            }
            moved
        }

        /// Moves the live `nodes` out with their half-edges.
        fn extract_nodes(&mut self, nodes: &BTreeSet<usize>) -> Partition {
            let (node_of, kept) = std::mem::take(&mut self.node_of)
                .into_iter()
                .partition(|(_, node)| nodes.contains(node));
            self.node_of = kept;
            self.live.retain(|node| !nodes.contains(node));
            Partition {
                node_of,
                live: nodes.clone(),
            }
        }

        fn append(&mut self, other: Partition) {
            self.node_of.extend(other.node_of);
            self.live.extend(other.live);
        }
    }

    /// The graph under test, its index of every logical node, the node that each
    /// index was merged into by pending identifications, and the model.
    struct Lane {
        graph: Graph,
        slots: Vec<NodeIndex>,
        merged_into: Vec<usize>,
        model: Partition,
        n_absorbed: usize,
    }

    impl Lane {
        fn new(builder: HedgeGraphBuilder<usize, usize, usize>) -> Self {
            let graph: Graph = builder.build();
            let model = Partition::of(&graph, "build");
            let mut lane = Lane {
                graph,
                slots: vec![],
                merged_into: vec![],
                model,
                n_absorbed: 0,
            };
            lane.renumbered();
            lane
        }

        fn n_slots(&self) -> usize {
            self.slots.len()
        }

        fn live(&self, slot: usize) -> usize {
            let mut node = self.slots[slot].0;
            while self.merged_into[node] != node {
                node = self.merged_into[node];
            }
            node
        }

        fn payload(&self, node: usize) -> usize {
            self.graph[NodeIndex(node)]
        }

        fn hedge_payloads(&self, subgraph: &SuBitGraph) -> BTreeSet<usize> {
            subgraph
                .included_iter()
                .map(|hedge| self.graph[hedge])
                .collect()
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

        /// Identifies possibly repeated or already absorbed nodes.
        fn identify(&mut self, slots: &[usize], payload: usize) {
            let resolved: BTreeSet<usize> = slots.iter().map(|&slot| self.live(slot)).collect();
            let merged = resolved.iter().map(|&node| self.payload(node)).collect();
            let requested: Vec<NodeIndex> = slots.iter().map(|&slot| self.slots[slot]).collect();
            let survivor = self.graph.identify_nodes(&requested, payload).0;
            assert!(resolved.contains(&survivor));
            for &node in &resolved {
                self.merged_into[node] = survivor;
                assert_eq!(
                    self.graph.node_store.is_absorbed(NodeIndex(node)),
                    node != survivor
                );
            }
            self.n_absorbed += resolved.len() - 1;
            self.model.identify(&merged, payload);
        }

        fn forget(&mut self) {
            let removed = self.graph.forget_identification_history().into_iter();
            assert_eq!(removed.count(), self.n_absorbed);
            self.n_absorbed = 0;
            self.renumbered();
        }

        fn check(&self, context: &str) {
            assert_eq!(Partition::of(&self.graph, context), self.model, "{context}");
        }

        fn step(&mut self, random: &mut Random, next: &mut usize) -> String {
            let mut fresh = || {
                *next += 1;
                *next
            };
            let n_slots = self.n_slots();
            // An empty graph can only grow.
            let edit = if n_slots == 0 { 13 } else { random.below(17) };
            match edit {
                16 => {
                    // Pending identifications neither block nor shift the new node.
                    let payload = fresh();
                    let node = self.graph.add_node(payload);
                    assert_eq!(node.0, self.merged_into.len());
                    self.model.live.insert(payload);
                    self.appended();
                    "add node".to_owned()
                }
                0..=4 => {
                    let slots: Vec<_> = (0..1 + random.below(4))
                        .map(|_| random.below(n_slots))
                        .collect();
                    self.identify(&slots, fresh());
                    format!("identify {slots:?}")
                }
                5 => {
                    let slot = random.below(n_slots);
                    let (hedge, edge) = (fresh(), fresh());
                    let flow = if random.one_in(2) {
                        Flow::Sink
                    } else {
                        Flow::Source
                    };
                    let node = self.live(slot);
                    (self.graph)
                        .add_dangling_edge_mut(NodeIndex(node).add_data(hedge), edge, flow, true)
                        .unwrap();
                    self.model.node_of.insert(hedge, self.payload(node));
                    format!("dangling at {slot}")
                }
                6 | 7 => {
                    let (a, b) = (random.below(n_slots), random.below(n_slots));
                    let (ha, hb, edge) = (fresh(), fresh(), fresh());
                    let (na, nb) = (self.live(a), self.live(b));
                    (self.graph)
                        .add_pair_mut(
                            NodeIndex(na).add_data(ha),
                            NodeIndex(nb).add_data(hb),
                            edge,
                            Orientation::Default,
                        )
                        .unwrap();
                    self.model.node_of.insert(ha, self.payload(na));
                    self.model.node_of.insert(hb, self.payload(nb));
                    format!("pair {a} {b}")
                }
                8 | 15 => {
                    self.forget();
                    "forget".to_owned()
                }
                9 => {
                    self.forget();
                    let subgraph = random_subgraph(random, self.graph.n_hedges());
                    let hedges = self.hedge_payloads(&subgraph);
                    self.graph.delete_hedges(&subgraph);
                    self.model.split_off(&hedges, 0);
                    self.renumbered();
                    format!("delete {subgraph:?}")
                }
                10 => {
                    self.forget();
                    let subgraph = random_subgraph(random, self.graph.n_hedges());
                    let hedges = self.hedge_payloads(&subgraph);
                    let split_offset = 1_000_000;
                    let extracted = self.graph.extract(
                        &subgraph,
                        |edge| edge.map(|data| *data),
                        |edge| edge,
                        |payload| payload + split_offset,
                        |payload| payload,
                    );
                    let expected = self.model.split_off(&hedges, split_offset);
                    assert_eq!(Partition::of(&extracted, "extracted"), expected);
                    self.renumbered();
                    format!("extract {subgraph:?}")
                }
                11 => {
                    self.forget();
                    let slots: Vec<_> = (0..random.below(3))
                        .map(|_| random.below(self.n_slots()))
                        .collect();
                    let nodes: BTreeSet<_> = slots.iter().map(|&slot| self.slots[slot]).collect();
                    let payloads = nodes.iter().map(|node| self.graph[*node]).collect();
                    let extracted =
                        self.graph
                            .extract_nodes(nodes, |edge| edge.map(|data| *data), |edge| edge);
                    let expected = self.model.extract_nodes(&payloads);
                    assert_eq!(Partition::of(&extracted, "extracted nodes"), expected);
                    self.renumbered();
                    format!("extract nodes {slots:?}")
                }
                12 => {
                    let n_hedges = self.graph.n_hedges();
                    if n_hedges == 0 {
                        return "contract nothing".to_owned();
                    }
                    let mut subgraph = random_subgraph(random, n_hedges);
                    subgraph.add(Hedge(random.below(n_hedges)));
                    let payload = fresh();
                    let touched = (subgraph.included_iter())
                        .map(|hedge| self.graph[self.graph.node_id(hedge)])
                        .collect();
                    let hedges = self.hedge_payloads(&subgraph);
                    self.graph.contract_subgraph(&subgraph, payload);
                    self.model.identify(&touched, payload);
                    self.model.split_off(&hedges, 0);
                    // Contraction forgets every pending identification.
                    self.n_absorbed = 0;
                    self.renumbered();
                    format!("contract {subgraph:?}")
                }
                _ => {
                    let other: Graph = random_graph(random, next).build();
                    self.model.append(Partition::of(&other, "other"));
                    if random.one_in(2) {
                        let graph =
                            std::mem::replace(&mut self.graph, HedgeGraphBuilder::new().build());
                        self.graph = graph
                            .join(
                                other,
                                |af, ad, bf, bd| af == -bf && ad.data % 2 == bd.data % 2,
                                |af, ad, _, _| (af, ad),
                            )
                            .unwrap();
                    } else {
                        self.graph.append_disconnected_mut(other).unwrap();
                    }
                    self.appended();
                    "join".to_owned()
                }
            }
        }
    }

    #[test]
    fn random_edits_keep_the_node_partition() {
        for seed in 0..250u64 {
            let mut random = Random(0x9e37_79b9_7f4a_7c15 ^ (seed + 1).wrapping_mul(0x2545_f491));
            let mut next = 0;
            let mut lane = Lane::new(random_graph(&mut random, &mut next));
            lane.check(&format!("seed {seed}: build"));
            let mut history = vec![];
            for _ in 0..24 {
                history.push(lane.step(&mut random, &mut next));
                lane.check(&format!("seed {seed}: {history:?}"));
            }
            lane.forget();
            lane.check(&format!("seed {seed}: final"));
        }
    }
}

#[test]
fn in_place_edge_insertion_preserves_graph_on_invalid_endpoints() {
    let mut builder = HedgeGraphBuilder::<&str, &str, &str>::new();
    let a = builder.add_node("a");
    let b = builder.add_node("b");
    builder.add_edge(a.add_data("a port"), b.add_data("b port"), "existing", true);
    let mut graph: HedgeGraph<_, _, _> = builder.build();
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

#[test]
fn in_place_edge_insertion_rejects_absorbed_endpoints_before_changing_the_graph() {
    let mut builder = HedgeGraphBuilder::<&str, &str, &str>::new();
    let [a, b, c] = ["a", "b", "c"].map(|data| builder.add_node(data));
    builder.add_edge(a.add_data("a port"), b.add_data("b port"), "existing", true);
    let mut graph: HedgeGraph<_, _, _> = builder.build();
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

#[test]
fn in_place_edge_insertion_keeps_indices_and_incidence_payloads() {
    let mut builder = HedgeGraphBuilder::<&str, &str, &str>::new();
    let a = builder.add_node("a");
    let b = builder.add_node("b");
    let isolated = builder.add_node("isolated");
    builder.add_edge(a.add_data("a port"), b.add_data("b port"), "existing", true);
    let mut graph: HedgeGraph<_, _, _> = builder.build();
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
