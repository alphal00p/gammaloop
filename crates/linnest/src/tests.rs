use std::{collections::BTreeMap, fs};

use figment::providers::Serialized;
use figment::{Figment, Profile};
use linnet::half_edge::layout::spring::{Constraint, ShiftDirection};
use linnet::half_edge::swap::Swap;
use linnet::{dot, parser::set::DotGraphSet};

use crate::{
    layout_graph_bytes, layout_parsed_graphs_bytes, parse_dot_graphs_bytes, CBORTypstGraph,
    TreeInitCfg, TypstGraph, TypstOutput,
};

fn test_figment() -> Figment {
    Figment::from(Serialized::from(
        BTreeMap::<String, String>::new(),
        Profile::Default,
    ))
}

fn empty_config_bytes() -> Vec<u8> {
    let mut bytes = Vec::new();
    ciborium::ser::into_writer(&BTreeMap::<String, String>::new(), &mut bytes).unwrap();
    bytes
}

#[test]
fn test_parse_pass_returns_archived_dot_graphs() {
    let parsed =
        parse_dot_graphs_bytes(br#"digraph first { a -> b } digraph second { c -> d }"#).unwrap();

    let archived: &<DotGraphSet as rkyv::Archive>::Archived =
        unsafe { DotGraphSet::archived_from_bytes(&parsed) };

    assert_eq!(archived.global_data.len(), 2);
    let globals = archived.global_data.as_slice();
    assert_eq!(globals[0].name.as_str(), "first");
    assert_eq!(globals[1].name.as_str(), "second");
}

#[test]
fn test_split_parse_layout_matches_wrapper() {
    let input = br#"digraph archived { a [pin="x:1.0"]; a -> b }"#;
    let config = empty_config_bytes();

    let parsed = parse_dot_graphs_bytes(input).unwrap();
    let split = layout_parsed_graphs_bytes(&parsed, &config).unwrap();
    let wrapped = layout_graph_bytes(input, &config).unwrap();

    assert_eq!(split, wrapped);

    let output: TypstOutput = ciborium::de::from_reader(split.as_slice()).unwrap();
    assert_eq!(output.graphs.len(), 1);
}

#[test]
fn dot_cbor() {
    let figment = test_figment();
    let g = TypstGraph::from_dot(dot!(digraph{ a; a->b}).unwrap(), &figment);

    let _cbor = g.to_cbor();
}

#[test]
fn test_cbor_serialization() {
    let figment = test_figment();
    let g = TypstGraph::from_dot(dot!(digraph{ a; a->b; b->c}).unwrap(), &figment);
    let cbor = g.to_cbor();

    let test_path = "test_graph.cbor";

    cbor.serialize_to_file(test_path)
        .expect("Failed to serialize to file");

    let deserialized =
        CBORTypstGraph::deserialize_from_file(test_path).expect("Failed to deserialize from file");

    assert_eq!(cbor.nodes.len(), deserialized.nodes.len());
    assert_eq!(cbor.edges.len(), deserialized.edges.len());

    fs::remove_file(test_path).ok();
}

#[test]
fn test_typst_graph_convenience_serialization() {
    let figment = test_figment();
    let g = TypstGraph::from_dot(dot!(digraph{ a->b; b->c}).unwrap(), &figment);
    let test_path = "test_convenience_graph.cbor";

    g.serialize_to_file(test_path)
        .expect("Failed to serialize using convenience method");

    let deserialized =
        CBORTypstGraph::deserialize_from_file(test_path).expect("Failed to deserialize from file");

    assert_eq!(g.to_cbor().nodes.len(), deserialized.nodes.len());
    assert_eq!(g.to_cbor().edges.len(), deserialized.edges.len());

    fs::remove_file(test_path).ok();
}

#[test]
fn test_pin_parsing() {
    let figment = test_figment();
    let mut g = TypstGraph::from_dot(dot!( digraph dot_80_0_GL208 {

       steps=1
       step=0.4
       beta =13.1
       k_spring=20.3;
       g_center=0
       gamma_ee=0.3
       gamma_ev=0.01
       length_scale = 0.25
       node[
         eval="(stroke:blue,fill :black,
         radius:2pt,
         outset: -2pt)"
       ]

       edge[
         eval=top
       ]
       v0[pin="x:@initial,y:@p1", style=invis]
       v1[pin="x:@initial,y:@p2",style=invis]
       v2[pin="x:@final,y:@p1",style=invis]
       v3[pin="x:@final,y:@p2",style=invis]
       v0 -> v11 [eval=photon]
       v1 -> v10 [eval="(..photon,label:[$gamma$],label-side: left)", mom_eval="(label:[$p_1$],label-sep:0mm)"]
       v9 -> v2 [eval=photon]
       v8 -> v3 [eval=photon]
       v4 -> v10
       v10 -> v5
       v5 -> v11 [dir=back]
       v11 -> v4
       v4 -> v7 [eval=gluon]
       v5 -> v6 [eval=gluon]
       v6 -> v8
       v8 -> v7
       v7 -> v9
       v9 -> v6
   }).unwrap(), &figment);

    g.layout();
    println!("{}", g.to_dot_graph().debug_dot())
}

#[test]
fn test_parsing() {
    let figment = test_figment();
    let g = TypstGraph::from_dot(
        dot!(digraph qqx_aaa_pentagon {
            steps=600
            step=.2
            beta =3.1
            k_spring=15.3;
            g_center=0
            gamma_ee=0.3
            gamma_ev=0.01
            length_scale = 0.2

             node[
              eval="(stroke:blue,fill :black,
          radius:2pt,
          outset: -2pt)"
            ]

        exte0 [style=invis pin="x:-4"];
        v3:0 -> exte0;
        })
        .unwrap(),
        &figment,
    );

    println!("{}", g.to_dot_graph().debug_dot())
}

#[test]
fn test_directional_constraints() {
    let figment = test_figment();
    let typst_graph = TypstGraph::from_dot(
        dot!(digraph {
            tree_dy = 2.0
            tree_dx = 3.0
            a [pin="x:+@group1"]
            b [pin="x:-@group1"]
            c [pin="y:+@group2"]
            d [pin="y:-@group2"]
            e -> f [pin="x:+@edge_group"]
            g -> h [pin="x:-@edge_group"]
            a -> b
            c -> d
        })
        .unwrap(),
        &figment,
    );
    let cfg = TreeInitCfg { dy: 2.0, dx: 3.0 };
    let (node_positions, edge_positions) = typst_graph.new_positions(cfg);

    let mut group1_positive = None;
    let mut group1_negative = None;
    let mut group2_positive = None;
    let mut group2_negative = None;

    for (node_idx, _, node) in typst_graph.iter_nodes() {
        match &node.constraints.x {
            Constraint::Grouped(_, ShiftDirection::PositiveOnly) => {
                group1_positive = Some(node_positions[node_idx].x);
            }
            Constraint::Grouped(_, ShiftDirection::NegativeOnly) => {
                group1_negative = Some(node_positions[node_idx].x);
            }
            _ => {}
        }
        match &node.constraints.y {
            Constraint::Grouped(_, ShiftDirection::PositiveOnly) => {
                group2_positive = Some(node_positions[node_idx].y);
            }
            Constraint::Grouped(_, ShiftDirection::NegativeOnly) => {
                group2_negative = Some(node_positions[node_idx].y);
            }
            _ => {}
        }
    }

    if let Some(pos_x) = group1_positive {
        assert!(pos_x >= 0.0);
    }

    if let Some(neg_x) = group1_negative {
        assert!(neg_x <= 0.0);
    }

    if let Some(pos_y) = group2_positive {
        assert!(pos_y >= 0.0);
    }

    if let Some(neg_y) = group2_negative {
        assert!(neg_y <= 0.0);
    }

    if let (Some(pos_x), Some(neg_x)) = (group1_positive, group1_negative) {
        assert!(pos_x > neg_x);
    }

    if let (Some(pos_y), Some(neg_y)) = (group2_positive, group2_negative) {
        assert!(pos_y > neg_y);
    }

    let mut edge_group_positive = None;
    let mut edge_group_negative = None;

    for (_, eid, edge_data) in typst_graph.iter_edges() {
        let edge = edge_data.data;
        match &edge.constraints.x {
            Constraint::Grouped(_, ShiftDirection::PositiveOnly) => {
                edge_group_positive = Some(edge_positions[eid].x);
            }
            Constraint::Grouped(_, ShiftDirection::NegativeOnly) => {
                edge_group_negative = Some(edge_positions[eid].x);
            }
            _ => {}
        }
    }

    if let Some(pos_edge_x) = edge_group_positive {
        assert!(pos_edge_x >= 0.0);
    }

    if let Some(neg_edge_x) = edge_group_negative {
        assert!(neg_edge_x <= 0.0);
    }

    if let (Some(pos_edge_x), Some(neg_edge_x)) = (edge_group_positive, edge_group_negative) {
        assert!(pos_edge_x > neg_edge_x);
    }
}

#[test]
fn test_full() {
    let input = stringify!(
       digraph {
           steps=600
           step=0.4
           beta =13.1
           k_spring=3.3;
           g_center=0
           gamma_dangling=50
           gamma_ee=0.3
           gamma_ev=0.01
           length_scale = 0.25
           a->b
       }
    );

    let figment = Figment::from(Serialized::from(
        BTreeMap::<String, String>::new(),
        Profile::Default,
    ));

    let dots = DotGraphSet::from_string_with_figment(input, figment.clone())
        .map_err(|a| a.to_string())
        .unwrap()
        .into_iter();

    let mut graphs = Vec::new();
    for g in dots {
        let mut typst_graph = TypstGraph::from_dot(g, &figment);

        typst_graph.layout();
        graphs.push((
            typst_graph.to_cbor(),
            typst_graph.to_dot_graph().debug_dot(),
        ));
    }
}
