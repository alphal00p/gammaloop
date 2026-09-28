use linnet::half_edge::{
    builder::{HedgeData, HedgeGraphBuilder},
    involution::{Flow, Hedge, Orientation},
    EdgeAccessors, HedgeGraph, NodeIndex,
};
use std::sync::Mutex;

type Data = Mutex<usize>;

fn fixture(orientation: Orientation) -> HedgeGraphBuilder<Data, Data, Data> {
    let mut builder = HedgeGraphBuilder::new();
    let nodes = [0, 1, 2].map(|n| builder.add_node(Mutex::new(n)));
    let endpoint = |node, hedge| HedgeData {
        data: Mutex::new(hedge),
        is_in_subgraph: false,
        node: nodes[node],
    };
    builder.add_edge(endpoint(1, 0), endpoint(0, 1), Mutex::new(10), orientation);
    builder.add_edge(endpoint(0, 2), endpoint(2, 3), Mutex::new(20), orientation);
    builder.add_edge(endpoint(2, 4), endpoint(1, 5), Mutex::new(30), orientation);
    builder.add_edge(endpoint(0, 6), endpoint(0, 7), Mutex::new(40), orientation);
    builder.add_external_edge(endpoint(1, 8), Mutex::new(50), orientation, Flow::Sink);
    builder
}

#[test]
fn retention_preserves_moved_payloads_and_boundary_edge_orientations() {
    for orientation in [
        Orientation::Default,
        Orientation::Reversed,
        Orientation::Undirected,
    ] {
        for order in [vec![1, 0], vec![2], vec![2, 0, 1], vec![]] {
            let original: HedgeGraph<Data, Data, Data> = fixture(orientation).build();
            let order = order.into_iter().map(NodeIndex).collect::<Vec<_>>();
            let original_ref = &original;
            let expected_hedges = order
                .iter()
                .flat_map(|&node| {
                    (0..original_ref.n_hedges())
                        .map(Hedge)
                        .filter(move |&h| original_ref.node_id(h) == node)
                })
                .collect::<Vec<_>>();
            let (builder, mapping) = fixture(orientation).retain_nodes_in_order(&order);
            let actual: HedgeGraph<Data, Data, Data> = builder.build();
            actual.check().unwrap();
            assert_eq!(actual.n_nodes(), order.len());
            assert_eq!(actual.n_hedges(), expected_hedges.len());
            for (new, &old) in expected_hedges.iter().enumerate() {
                let new = Hedge(new);
                assert_eq!(mapping[old], Some(new));
                assert_eq!(*actual[new].lock().unwrap(), *original[old].lock().unwrap());
                assert_eq!(
                    *actual[[&new]].lock().unwrap(),
                    *original[[&old]].lock().unwrap()
                );
                assert_eq!(actual.inv(new), mapping[original.inv(old)].unwrap_or(new));
                assert_eq!(actual.flow(new), original.flow(old));
                assert_eq!(actual.orientation(new), original.orientation(old));
                assert_eq!(order[actual.node_id(new).0], original.node_id(old));
            }
            for (new, &old) in order.iter().enumerate() {
                assert_eq!(
                    *actual[NodeIndex(new)].lock().unwrap(),
                    *original[old].lock().unwrap()
                );
            }
            for old in (0..original.n_hedges()).map(Hedge) {
                assert_eq!(
                    mapping[old].is_some(),
                    order.contains(&original.node_id(old))
                );
            }
        }
    }
}

#[test]
#[should_panic(expected = "retained nodes must be unique")]
fn duplicate_retained_node_is_rejected() {
    let _ = fixture(Orientation::Default).retain_nodes_in_order(&[NodeIndex(0), NodeIndex(0)]);
}

#[test]
#[should_panic(expected = "index out of bounds")]
fn unknown_retained_node_is_rejected() {
    let _ = fixture(Orientation::Default).retain_nodes_in_order(&[NodeIndex(3)]);
}

#[test]
fn builder_connections_match_finished_connections_and_retained_boundaries() {
    use linnet::half_edge::involution::EdgeData;
    use std::{cell::RefCell, collections::BTreeMap};
    type Graph = HedgeGraph<Vec<usize>, usize, usize>;
    for orientation in [
        Orientation::Default,
        Orientation::Reversed,
        Orientation::Undirected,
    ] {
        for flow in [Flow::Source, Flow::Sink] {
            for order in [
                vec![NodeIndex(2), NodeIndex(0), NodeIndex(1)],
                vec![NodeIndex(1), NodeIndex(0)],
                vec![],
            ] {
                let fixture = || {
                    let mut builder = HedgeGraphBuilder::new();
                    let nodes = [0, 1, 2].map(|n| builder.add_node(n));
                    for (hedge, node) in [0, 1, 2, 0, 1, 2].into_iter().enumerate() {
                        builder.add_external_edge(
                            HedgeData {
                                data: hedge,
                                is_in_subgraph: false,
                                node: nodes[node],
                            },
                            vec![hedge],
                            orientation,
                            if hedge % 2 == 0 {
                                Flow::Source
                            } else {
                                Flow::Sink
                            },
                        );
                    }
                    builder
                };
                let (early_calls, late_calls) =
                    (RefCell::new(Vec::new()), RefCell::new(Vec::new()));
                let merge = |calls: &RefCell<Vec<_>>,
                             left_flow,
                             mut left: EdgeData<Vec<usize>>,
                             right_flow,
                             right: EdgeData<Vec<usize>>| {
                    calls
                        .borrow_mut()
                        .push((left_flow, left.clone(), right_flow, right.clone()));
                    left.data.extend(right.data);
                    (flow, left)
                };
                let mut early = fixture();
                let mut late: Graph = fixture().build();
                for (left, right) in [
                    (Hedge(4), Hedge(0)),
                    (Hedge(1), Hedge(5)),
                    (Hedge(2), Hedge(3)),
                ] {
                    early
                        .connect_identities(left, right, |lf, l, rf, r| {
                            merge(&early_calls, lf, l, rf, r)
                        })
                        .unwrap();
                    late.connect_identities(left, right, |lf, l, rf, r| {
                        merge(&late_calls, lf, l, rf, r)
                    });
                }
                assert_eq!(*early_calls.borrow(), *late_calls.borrow());
                let (early, mapping) = early.retain_nodes_in_order(&order);
                let actual: Graph = early.build();
                actual.check().unwrap();
                assert_eq!(actual.n_nodes(), order.len());
                let expected_hedges = order
                    .iter()
                    .flat_map(|&node| late.iter_crown(node))
                    .collect::<Vec<_>>();
                assert_eq!(actual.n_hedges(), expected_hedges.len());
                let mut edge_ids = BTreeMap::new();
                let mut inverse = BTreeMap::new();
                for (new, &old) in expected_hedges.iter().enumerate() {
                    let new = Hedge(new);
                    assert_eq!(mapping[old], Some(new));
                    assert_eq!(actual[new], late[old]);
                    assert_eq!(actual[[&new]], late[[&old]]);
                    assert_eq!(actual.inv(new), mapping[late.inv(old)].unwrap_or(new));
                    assert_eq!(actual.flow(new), late.flow(old));
                    assert_eq!(actual.orientation(new), late.orientation(old));
                    assert_eq!(order[actual.node_id(new).0], late.node_id(old));
                    if let Some(previous) = edge_ids.insert(late[&old], actual[&new]) {
                        assert_eq!(previous, actual[&new]);
                    }
                    if let Some(previous) = inverse.insert(actual[&new], late[&old]) {
                        assert_eq!(previous, late[&old]);
                    }
                }
                assert_eq!(inverse.len(), actual.iter_edges().count());
                for (new, &old) in order.iter().enumerate() {
                    assert_eq!(actual[NodeIndex(new)], late[old]);
                    assert_eq!(
                        actual.iter_crown(NodeIndex(new)).collect::<Vec<_>>(),
                        late.iter_crown(old)
                            .map(|h| mapping[h].unwrap())
                            .collect::<Vec<_>>()
                    );
                }
            }
        }
    }
}

#[test]
fn builder_connections_reject_already_paired_endpoints_without_calling_merge() {
    let mut builder = HedgeGraphBuilder::<usize, usize>::new();
    let a = builder.add_node(1);
    let b = builder.add_node(2);
    builder.add_edge(a, b, 7, Orientation::Default);
    builder.add_external_edge(a, 11, Orientation::Reversed, Flow::Source);
    assert!(builder
        .connect_identities(Hedge(0), Hedge(2), |_, _, _, _| panic!(
            "invalid seam must not merge"
        ))
        .is_err());
    let graph: HedgeGraph<usize, usize> = builder.build();
    graph.check().unwrap();
    assert_eq!(graph.n_hedges(), 3);
    assert_eq!(graph.inv(Hedge(0)), Hedge(1));
    assert_eq!(graph.inv(Hedge(2)), Hedge(2));
    assert_eq!(graph[[&Hedge(0)]], 7);
    assert_eq!(graph[[&Hedge(2)]], 11);
}
