use crate::{
    dot,
    half_edge::{
        builder::HedgeGraphBuilder,
        involution::{Flow, Hedge},
        nodestore::NodeStorageOps,
        subgraph::{ModifySubSet, SuBitGraph},
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
fn orientation_hedges() {
    let mut single_node = HedgeGraphBuilder::new();
    let a = single_node.add_node(());
    single_node.add_external_edge(a, (), true, Flow::Source);
    single_node.add_external_edge(a, (), true, Flow::Sink);
    let aligned: HedgeGraph<(), (), ()> = single_node.build();

    println!("{}", aligned.base_dot())
}
