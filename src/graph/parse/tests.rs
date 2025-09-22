use insta::assert_snapshot;
use linnet::half_edge::involution::{EdgeIndex, HedgePair};
use log::info;
use spenso::{
    network::{library::DummyLibrary, parsing::ShadowedStructure, store::NetworkStore, Network},
    structure::{HasName, PermutedStructure},
    tensors::symbolic::SymbolicTensor,
};
use symbolica::atom::{Atom, FunctionBuilder};
use typed_index_collections::ti_vec;

use super::Graph;
use crate::{
    dot,
    feyngen::diagram_generator::{EdgeColor, NodeColorWithVertexRule},
    graph::{
        parse::{complete_group_parsing, IntoGraph},
        GraphGroup, LMBext,
    },
    initialisation::test_initialise,
    momentum_sample::LoopIndex,
    numerator::{aind::Aind, Numerator, UnInit},
    utils::test_utils::load_generic_model,
};

#[test]
fn symbolica_parse() {
    test_initialise().unwrap();
    let model = load_generic_model("sm");

    let mut a = symbolica::graph::Graph::new();

    let udwm = NodeColorWithVertexRule::from_particles(["d~", "W-", "u"], &model);
    let udwp = NodeColorWithVertexRule::from_particles(["d", "W+", "u~"], &model);
    let cswp = NodeColorWithVertexRule::from_particles(["s", "W+", "c~"], &model);
    let cswm = NodeColorWithVertexRule::from_particles(["s~", "W-", "c"], &model);

    let ext1 = NodeColorWithVertexRule {
        external_tag: 1,
        vertex_rule: model.vertex_rules[0].clone(),
    };
    let ext2 = NodeColorWithVertexRule {
        external_tag: 2,
        vertex_rule: model.vertex_rules[0].clone(),
    };
    let ext3 = NodeColorWithVertexRule {
        external_tag: 3,
        vertex_rule: model.vertex_rules[0].clone(),
    };
    let ext4 = NodeColorWithVertexRule {
        external_tag: 4,
        vertex_rule: model.vertex_rules[0].clone(),
    };

    let e1 = a.add_node(ext1.clone());
    let e2 = a.add_node(ext2.clone());
    let e3 = a.add_node(ext3.clone());
    let e4 = a.add_node(ext4.clone());
    let v1 = a.add_node(udwp.clone());
    let v2 = a.add_node(cswm.clone());
    let v3 = a.add_node(udwm.clone());
    let v4 = a.add_node(cswp.clone());

    let ed = EdgeColor::from_particle(model.get_particle("d"));
    let es = EdgeColor::from_particle(model.get_particle("s"));
    let ec = EdgeColor::from_particle(model.get_particle("c"));
    let ewp = EdgeColor::from_particle(model.get_particle("W+"));
    let ewm = EdgeColor::from_particle(model.get_particle("W-"));
    let eu = EdgeColor::from_particle(model.get_particle("u"));

    a.add_edge(e1, v1, true, ed).unwrap();
    a.add_edge(e2, v2, true, ec).unwrap();

    a.add_edge(v3, e3, true, ed).unwrap();
    a.add_edge(v4, e4, true, ec).unwrap();

    a.add_edge(v1, v3, true, eu).unwrap();
    a.add_edge(v2, v4, true, es).unwrap();
    a.add_edge(v1, v2, false, ewm).unwrap();
    a.add_edge(v3, v4, false, ewp).unwrap();

    println!("{}", a.to_dot());

    let g = Graph::from_symbolica_graph(
        &model,
        "test",
        &a,
        Atom::num(1),
        &[
            (Some(1), None),
            (Some(2), None),
            (None, Some(3)),
            (None, Some(4)),
        ],
    )
    .unwrap();

    println!("{}", g.dot_serialize());

    let g = Graph::from_symbolica_graph(
        &model,
        "test",
        &a,
        Atom::num(1),
        &[
            (Some(2), None),
            (Some(3), None),
            (None, Some(1)),
            (None, Some(4)),
        ],
    )
    .unwrap();

    println!("{}", g.dot_serialize());

    let g = Graph::from_symbolica_graph(
        &model,
        "test",
        &a,
        Atom::num(1),
        &[(Some(1), Some(3)), (Some(2), Some(4))],
    )
    .unwrap();

    println!("{}", g.dot_serialize());

    let mut a = symbolica::graph::Graph::new();

    let e1 = a.add_node(ext1.clone());
    let e2 = a.add_node(ext2.clone());
    let e3 = a.add_node(ext3.clone());
    let e4 = a.add_node(ext4.clone());
    let v1 = a.add_node(udwm.clone());
    let v2 = a.add_node(udwp.clone());
    let v3 = a.add_node(udwp.clone());
    let v4 = a.add_node(udwm.clone());

    a.add_edge(e1, v1, true, eu).unwrap();
    a.add_edge(e2, v2, true, ed).unwrap();

    a.add_edge(v3, e3, true, eu).unwrap();
    a.add_edge(v4, e4, true, ed).unwrap();

    a.add_edge(v1, v3, true, ed).unwrap();
    a.add_edge(v2, v4, true, eu).unwrap();
    a.add_edge(v1, v2, false, ewp).unwrap();
    a.add_edge(v3, v4, false, ewm).unwrap();

    let g = Graph::from_symbolica_graph(
        &model,
        "test",
        &a,
        Atom::num(1),
        &[(Some(1), Some(3)), (Some(2), Some(4))],
    )
    .unwrap();

    println!("{}", g.dot_serialize());

    // let g = Graph::from_symbolica_graph(
    //     &model,
    //     "test",
    //     &a,
    //     Atom::num(1),
    //     &[(Some(1), Some(4)), (Some(2), Some(3))],
    // )
    // .unwrap();

    // println!("{}", g.dot_serialize());
}

#[test]
fn test_load() {
    test_initialise().unwrap();
    info!("Loading graph");
    let graph: Vec<Graph> = dot!(digraph triangle {
        graph [
            overall_factor = 1;
            multiplicity_factor = 1;
        ]
        edge [
            pdg=1000
        ]
        ext [style=invis]
        ext -> v4
        ext -> v5
        v6 -> ext
        v5 -> v4 [lmb_index=0];
        v6 -> v5;
        v4 -> v6 ;
    }

    digraph triangle2 {
        graph [
            overall_factor = 1;
            multiplicity_factor = 1;
        ]
        edge [
            pdg=1000
        ]
        ext [style=invis]
        ext -> v4
        ext -> v5
        v6 -> ext [id = 2]
        v5 -> v4 [lmb_index=0];
        v6 -> v5;
        v4 -> v6 ;
    }
    ,"scalars")
    .unwrap();

    println!("{}", graph[0].dot_serialize());

    println!(
        "{}",
        graph[0].dot_lmb(&graph[0].underlying.full(), &graph[0].loop_momentum_basis)
    );

    println!("{}", graph[0].loop_momentum_basis);

    println!(
        "{}",
        graph[1].dot_lmb(&graph[1].underlying.full(), &graph[1].loop_momentum_basis)
    );
    println!("{}", graph[1].loop_momentum_basis);
}

#[test]
fn xs_glueing() {
    test_initialise().unwrap();
    let g: Graph = dot!(
        digraph{
          num = "1";
          overall_factor = "1";
          0[int_id=V_89];
          1[int_id=V_127];
          2[int_id=V_123];
          3[int_id=V_93];

          ext0	 [style=invis is_cut=0];
          ext1	 [style=invis is_cut=1];
          2	-> ext0	     [id=0  particle="d"];
          ext0	-> 0	 [id=1  particle="d"];
          3:2	-> ext1	 [id=2  particle="c"];
          ext1	-> 1:3	 [id=3  particle="c"];
          0:4	-> 2:5	 [id=4  particle="u"];
          1:6	-> 3:7	 [id=5  particle="s"];
          0:8	-> 1:9	 [id=6  particle="W-"];
          2:10	-> 3:11	 [id=7  particle="W+"];
        }
    )
    .unwrap();
    assert_snapshot!(g.debug_dot(),@r#"
    digraph {
      num = "1";
      overall_factor = "1";
      projector = "ubar(0,spenso::bis(4,hedge(0)))*ubar(1,spenso::bis(4,hedge(2)))*u(0,spenso::bis(4,hedge(1)))*u(1,spenso::bis(4,hedge(3)))";
      0[dod=0 int_id=V_89 num="UFO::GC_100*spenso::g(spenso::dind(spenso::cof(3,hedge(4))),spenso::cof(3,hedge(1)))*spenso::gamma(spenso::bis(4,hedge(4)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(8)))*spenso::projm(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(1)))"];
      1[dod=0 int_id=V_127 num="UFO::GC_45*spenso::g(spenso::dind(spenso::cof(3,hedge(6))),spenso::cof(3,hedge(3)))*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(9)))*spenso::projm(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3)))"];
      2[dod=0 int_id=V_123 num="UFO::GC_41*spenso::g(spenso::dind(spenso::cof(3,hedge(0))),spenso::cof(3,hedge(5)))*spenso::gamma(spenso::bis(4,hedge(0)),spenso::bis(4,vertex(2,1)),spenso::mink(4,hedge(10)))*spenso::projm(spenso::bis(4,vertex(2,1)),spenso::bis(4,hedge(5)))"];
      3[dod=0 int_id=V_93 num="UFO::GC_104*spenso::g(spenso::dind(spenso::cof(3,hedge(2))),spenso::cof(3,hedge(7)))*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(3,1)),spenso::mink(4,hedge(11)))*spenso::projm(spenso::bis(4,vertex(3,1)),spenso::bis(4,hedge(7)))"];

      2:3	-> 0:0	 [id=0 source=0 sink=1  dod=-2 is_cut=0 is_dummy=false lmb_rep="P(0,a___)" name=e0 num="spenso::g(spenso::dind(spenso::cof(3,hedge(1))),spenso::cof(3,hedge(0)))" particle="d"];
      3:2	-> 1:1	 [id=1 source=0 sink=1  dod=-2 is_cut=1 is_dummy=false lmb_rep="P(1,a___)" name=e1 num="spenso::g(spenso::dind(spenso::cof(3,hedge(3))),spenso::cof(3,hedge(2)))" particle="c"];
      2:10	-> 3:11	 [id=2 dir=none source=2 sink=2  dod=0 is_dummy=false lmb_rep="-P(0,a___)+K(0,a___)" name=e2 num="-spenso::g(spenso::mink(4,hedge(10)),spenso::mink(4,hedge(11)))+UFO::MW^-2*Q(2,spenso::mink(4,hedge(10)))*Q(2,spenso::mink(4,hedge(11)))" particle="W+"];
      0:8	-> 1:9	 [id=3 dir=none source=2 sink=2  dod=0 is_dummy=false lmb_rep="P(0,a___)-K(0,a___)" name=e3 num="-spenso::g(spenso::mink(4,hedge(8)),spenso::mink(4,hedge(9)))+UFO::MW^-2*Q(3,spenso::mink(4,hedge(8)))*Q(3,spenso::mink(4,hedge(9)))" particle="W-"];
      0:4	-> 2:5	 [id=4 source=0 sink=1  dod=-1 is_dummy=false lmb_id=0 lmb_rep="K(0,a___)" name=e4 num="Q(4,spenso::mink(4,edge(4,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(5))),spenso::cof(3,hedge(4)))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(4)),spenso::mink(4,edge(4,1)))" particle="u"];
      1:6	-> 3:7	 [id=5 source=0 sink=1  dod=-1 is_dummy=false lmb_rep="P(0,a___)+P(1,a___)-K(0,a___)" name=e5 num="Q(5,spenso::mink(4,edge(5,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(7))),spenso::cof(3,hedge(6)))*spenso::gamma(spenso::bis(4,hedge(7)),spenso::bis(4,hedge(6)),spenso::mink(4,edge(5,1)))" particle="s"];
    }
    "#);

    let g: Graph = dot!(
        digraph{
            a->b [particle= c];
            b->a [particle= d];
            b->a [particle= "W+"];
            a[int_id=V_90]
        }
    )
    .unwrap();
    assert_snapshot!(g.debug_dot(),@r#"
    digraph {
      num = "1";
      overall_factor = "1";
      projector = "1";
      a[dod=0 int_id=V_90 num="UFO::GC_103*spenso::g(spenso::dind(spenso::cof(3,hedge(0))),spenso::cof(3,hedge(3)))*spenso::gamma(spenso::bis(4,hedge(0)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(5)))*spenso::projm(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(3)))"];
      b[dod=0 int_id=V_126 num="UFO::GC_44*spenso::g(spenso::dind(spenso::cof(3,hedge(2))),spenso::cof(3,hedge(1)))*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(4)))*spenso::projm(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(1)))"];

      a:0	-> b:1	 [id=0 source=0 sink=1  dod=-1 is_dummy=false lmb_id=0 lmb_rep="K(0,a___)" name=e0 num="(-UFO::MC*spenso::g(spenso::bis(4,hedge(0)),spenso::bis(4,hedge(1)))+Q(0,spenso::mink(4,edge(0,1)))*spenso::gamma(spenso::bis(4,hedge(1)),spenso::bis(4,hedge(0)),spenso::mink(4,edge(0,1))))*spenso::g(spenso::dind(spenso::cof(3,hedge(1))),spenso::cof(3,hedge(0)))" particle="c"];
      b:2	-> a:3	 [id=1 source=0 sink=1  dod=-1 is_dummy=false lmb_id=1 lmb_rep="K(1,a___)" name=e1 num="Q(1,spenso::mink(4,edge(1,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(3))),spenso::cof(3,hedge(2)))*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(2)),spenso::mink(4,edge(1,1)))" particle="d"];
      b:4	-> a:5	 [id=2 dir=none source=2 sink=2  dod=0 is_dummy=false lmb_rep="K(0,a___)-K(1,a___)" name=e2 num="-spenso::g(spenso::mink(4,hedge(4)),spenso::mink(4,hedge(5)))+UFO::MW^-2*Q(2,spenso::mink(4,hedge(4)))*Q(2,spenso::mink(4,hedge(5)))" particle="W+"];
    }
    "#);
}

#[test]
fn vertex_normalization() {
    test_initialise().unwrap();
    let g: Graph = dot!(
        digraph GL1{
            num = "1";
            overall_factor = "AutG(1)^-1*InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)*AntiFermionSpinSumSign(1)";
            0[int_id=V_117];
            1[int_id=V_82];
            2[int_id=V_71];
            3[int_id=V_11];

            ext0 [style=invis];
            2:0 -> ext0 [id=0 dir=back is_cut=0 is_dummy=false particle="a"];
            ext1 [style=invis];
            ext1 -> 3:1 [id=1 is_cut=0 is_dummy=false particle="a"];
            0:2 -> 1:3 [id=2 is_dummy=false particle="c"];
            0:4 -> 2:5 [id=3 is_dummy=false particle="d~"];
            0:6 -> 3:7 [id=4 is_dummy=false particle="G-"];
            1:8 -> 2:9 [id=5 is_dummy=false particle="d"];
            1:10 -> 3:11 [id=6 is_dummy=false particle="G+"];
        }
    ).unwrap();
    assert_snapshot!(g.debug_dot(),@"");
}

#[test]
fn validate_verterx() {
    test_initialise().unwrap();
    let g: Graph = dot!(digraph GL1{
      //  0[dod=0 int_id=V_82 num="UFO::GC_16*spenso::g(spenso::dind(spenso::cof(3,hedge(2))),spenso::dind(spenso::cof(3,hedge(4))))*spenso::projp(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(4)))"];
      // 1[dod=0 int_id=V_117 num="UFO::GC_22*spenso::g(spenso::dind(spenso::cof(3,hedge(8))),spenso::cof(3,hedge(3)))*spenso::projm(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(3)))"];
      // 2[dod=0 int_id=V_71 num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(5)),spenso::cof(3,hedge(9)))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(9)),spenso::mink(4,hedge(0)))"];
      // 3[dod=0 int_id=V_11 num="UFO::GC_3*(UFO::P(spenso::mink(4,hedge(1)),2)-UFO::P(spenso::mink(4,hedge(1)),3))"];

      2:1   -> 3:0   [id=0 is_cut=0 particle="a"];
      1:8   -> 2:9   [id=1  particle="d"];
      0:2   -> 1:3   [id=2 particle="c"];
      0:4   -> 2:5   [id=3  particle="d~"];
      0:6   -> 3:7   [id=4  particle="G-"];
      1:10  -> 3:11  [id=5 particle="G+"];
    })
    .unwrap();

    assert_snapshot!(g.debug_dot(),@r#"
    digraph GL1{
      num = "1";
      overall_factor = "1";
      projector = "ϵ(0,spenso::mink(4,hedge(0)))*ϵbar(0,spenso::mink(4,hedge(1)))";
      0[dod=0 int_id=V_82 num="UFO::GC_16*spenso::g(spenso::dind(spenso::cof(3,hedge(2))),spenso::dind(spenso::cof(3,hedge(4))))*spenso::projp(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(4)))"];
      1[dod=0 int_id=V_117 num="UFO::GC_22*spenso::g(spenso::dind(spenso::cof(3,hedge(8))),spenso::cof(3,hedge(3)))*spenso::projm(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(3)))"];
      2[dod=0 int_id=V_71 num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(5)),spenso::cof(3,hedge(9)))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(9)),spenso::mink(4,hedge(1)))"];
      3[dod=1 int_id=V_11 num="UFO::GC_3*(Q(4,spenso::mink(4,hedge(0)))-Q(5,spenso::mink(4,hedge(0))))"];

      2:1	-> 3:0	 [id=0 dir=none source=2 sink=0  dod=-2 is_cut=0 is_dummy=false lmb_rep="P(0,a___)" name=e0 num="1" particle="a"];
      1:8	-> 2:9	 [id=1 source=0 sink=1  dod=-1 is_dummy=false lmb_rep="P(0,a___)-K(1,a___)" name=e1 num="Q(1,spenso::mink(4,edge(1,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(9))),spenso::cof(3,hedge(8)))*spenso::gamma(spenso::bis(4,hedge(9)),spenso::bis(4,hedge(8)),spenso::mink(4,edge(1,1)))" particle="d"];
      0:2	-> 1:3	 [id=2 source=0 sink=1  dod=-1 is_dummy=false lmb_id=0 lmb_rep="K(0,a___)" name=e2 num="(-UFO::MC*spenso::g(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(3)))+Q(2,spenso::mink(4,edge(2,1)))*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(2)),spenso::mink(4,edge(2,1))))*spenso::g(spenso::dind(spenso::cof(3,hedge(3))),spenso::cof(3,hedge(2)))" particle="c"];
      0:4	-> 2:5	 [id=3 dir=back source=1 sink=0  dod=-1 is_dummy=false lmb_id=1 lmb_rep="K(1,a___)" name=e3 num="Q(3,spenso::mink(4,edge(3,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(5))),spenso::cof(3,hedge(4)))*spenso::gamma(spenso::bis(4,hedge(4)),spenso::bis(4,hedge(5)),spenso::mink(4,edge(3,1)))" particle="d~"];
      0:6	-> 3:7	 [id=4 dir=none source=2 sink=1  dod=-2 is_dummy=false lmb_rep="-K(0,a___)-K(1,a___)" name=e4 num="1" particle="G-"];
      1:10	-> 3:11	 [id=5 dir=none source=2 sink=2  dod=-2 is_dummy=false lmb_rep="-P(0,a___)+K(0,a___)+K(1,a___)" name=e5 num="1" particle="G+"];
    }
    "#);
}

#[test]
fn test_ufo() {
    test_initialise().unwrap();
    let g: Graph = dot!(
        digraph GL0{
            num = "1";
            overall_factor = "AutG(1)^-1*InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)*AntiFermionSpinSumSign(1)";
            0[int_id=V_79];
        1[int_id=V_79];
        2[int_id=V_71];
        3[int_id=V_71];
        ext0 [style=invis];
        2:0-> ext0 [id=0 dir=back is_cut=0 is_dummy=false particle="a"];
        ext1 [style=invis];
        ext1-> 3:1 [id=1 is_cut=0 is_dummy=false particle="a"];
        0:2-> 1:3 [id=2  is_dummy=false particle="d"];
        0:4-> 1:5 [id=3  is_dummy=false particle="Z"];
        0:6-> 3:7 [id=4  is_dummy=false particle="d~"];
        1:8-> 2:9 [id=5  is_dummy=false particle="d"];
        2:10-> 3:11 [id=6  is_dummy=false particle="d"];
        }
    ).unwrap();

    assert_snapshot!(g.debug_dot(),@r#"
    digraph GL0{
      num = "1";
      overall_factor = "AutG(1)^-1*InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)*AntiFermionSpinSumSign(1)";
      projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵbar(0,spenso::mink(4,hedge(0)))";
      0[dod=0 int_id=V_79 num="(UFO::GC_58*(-2*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::projp(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(6)))+spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::projm(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(6))))+UFO::GC_50*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::projm(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(6))))*spenso::g(spenso::dind(spenso::cof(3,hedge(2))),spenso::dind(spenso::cof(3,hedge(6))))"];
      1[dod=0 int_id=V_79 num="(UFO::GC_58*(-2*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::projp(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3)))+spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::projm(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3))))+UFO::GC_50*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::projm(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3))))*spenso::g(spenso::dind(spenso::cof(3,hedge(8))),spenso::cof(3,hedge(3)))"];
      2[dod=0 int_id=V_71 num="UFO::GC_1*spenso::g(spenso::dind(spenso::cof(3,hedge(10))),spenso::cof(3,hedge(9)))*spenso::gamma(spenso::bis(4,hedge(10)),spenso::bis(4,hedge(9)),spenso::mink(4,hedge(0)))"];
      3[dod=0 int_id=V_71 num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(7)),spenso::cof(3,hedge(11)))*spenso::gamma(spenso::bis(4,hedge(7)),spenso::bis(4,hedge(11)),spenso::mink(4,hedge(1)))"];

      2:1	-> 3:0	 [id=0 dir=none source=2 sink=2  dod=-2 is_cut=0 is_dummy=false lmb_rep="P(0,a___)" name=e0 num="1" particle="a"];
      1:8	-> 2:9	 [id=1 source=0 sink=1  dod=-1 is_dummy=false lmb_rep="K(0,a___)+K(1,a___)" name=e1 num="Q(1,spenso::mink(4,edge(1,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(9))),spenso::cof(3,hedge(8)))*spenso::gamma(spenso::bis(4,hedge(9)),spenso::bis(4,hedge(8)),spenso::mink(4,edge(1,1)))" particle="d"];
      0:2	-> 1:3	 [id=2 source=0 sink=1  dod=-1 is_dummy=false lmb_id=0 lmb_rep="K(0,a___)" name=e2 num="Q(2,spenso::mink(4,edge(2,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(3))),spenso::cof(3,hedge(2)))*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(2)),spenso::mink(4,edge(2,1)))" particle="d"];
      0:4	-> 1:5	 [id=3 dir=none source=2 sink=2  dod=0 is_dummy=false lmb_id=1 lmb_rep="K(1,a___)" name=e3 num="-spenso::g(spenso::mink(4,hedge(4)),spenso::mink(4,hedge(5)))+UFO::MZ^-2*Q(3,spenso::mink(4,hedge(4)))*Q(3,spenso::mink(4,hedge(5)))" particle="Z"];
      0:6	-> 3:7	 [id=4 dir=back source=1 sink=0  dod=-1 is_dummy=false lmb_rep="-K(0,a___)-K(1,a___)" name=e4 num="Q(4,spenso::mink(4,edge(4,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(7))),spenso::cof(3,hedge(6)))*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,hedge(7)),spenso::mink(4,edge(4,1)))" particle="d~"];
      2:10	-> 3:11	 [id=5 source=0 sink=1  dod=-1 is_dummy=false lmb_rep="-P(0,a___)+K(0,a___)+K(1,a___)" name=e5 num="Q(5,spenso::mink(4,edge(5,1)))*spenso::g(spenso::dind(spenso::cof(3,hedge(11))),spenso::cof(3,hedge(10)))*spenso::gamma(spenso::bis(4,hedge(11)),spenso::bis(4,hedge(10)),spenso::mink(4,edge(5,1)))" particle="d"];
    }
    "#);
}
#[test]
fn test_loop_momentum_basis() {
    test_initialise().unwrap();
    let g: Graph = dot!(
        digraph triangle {
        graph [
        overall_factor = 1;
        multiplicity_factor = 1;
        ]
        edge [
        pdg=1000
        ]
        v4 [num=1]
        ext [style=invis]
        ext -> v4:0 [id=1 is_dummy=true]
        ext -> v5:1 [id=0]
        v6:2 -> ext [id=2]
        v5 -> v4 [lmb_index=0];
        v6 -> v5;
        v4 -> v6 ;
        },"scalars"
    )
    .unwrap();
    // println!("{}", g.loop_momentum_basis);

    insta::assert_ron_snapshot!(g.loop_momentum_basis);

    println!(
        "{}",
        g.dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
    );
}

#[test]
fn scalar_without_model() {
    test_initialise().unwrap();
    let g: Graph = dot!(
        digraph G{
            ext    [style=invis]
            node [num=1]
            ext -> A
            C -> A
            A -> D
            D -> B
            B -> C
            C -> D
            B -> ext
        }
    )
    .unwrap();

    println!("{}", g.dot_serialize())
}

#[test]
fn parse() {
    test_initialise().unwrap();
    let g: Graph = dot!(
        digraph G{
            ext    [style=invis]
            ext -> A   [particle=a]
            C -> A    [particle="b"]
            A -> B
            A -> D    [particle="b"]
            D -> B    [particle="b"]
            B -> C    [particle="b"]
            C -> D    [particle="g"]
            B -> ext   [particle=a]
        }
    )
    .unwrap();

    println!("{:#?}", g.loop_momentum_basis);

    println!(
        "{}",
        g.underlying
            .dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
    );

    println!("{}", g.dot_serialize());
    // return;
    let num = Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter());

    let expr = num.state.expr.as_view();

    let lib: DummyLibrary<SymbolicTensor<Aind>> = DummyLibrary::<_>::new();
    let net =
        Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, _, Aind>::try_from_view(expr, &lib)
            .unwrap();

    println!("{}", expr);
    println!(
        "{}",
        net.dot_display_impl(
            |a| a.to_string(),
            |_| None,
            |a| {
                if let Ok(a) =
                    PermutedStructure::<ShadowedStructure<Aind>>::try_from(a.expression.as_view())
                {
                    a.structure
                        .name()
                        .map(|s| {
                            if let Some(a) = a.structure.args() {
                                FunctionBuilder::new(s).add_args(&a).finish().to_string()
                            } else {
                                s.to_string()
                            }
                        })
                        .unwrap_or("".to_string())
                } else {
                    "".to_string()
                }
            }
        )
    );

    let num = num.color_simplify();

    println!("{}", num.state.expr);
    // println!("{}", num.state.expr.to_dots());
    // let num = num.gamma_simplify();

    // println!("{}", num.state.expr);
}
#[test]
fn parse_local() {
    test_initialise().unwrap();
    let graphs: Vec<Graph> = dot!(
        digraph G{
            edge [num="1" dod=-1000 pdg=1000]
            node [num="1" dod=-1000]
            ext     [style=invis]
            A [id=1]
            B [id=0]
            C
            D [id=2]
            ext -> A  [id=1]
            ext -> C
            C -> A    [num="P(eid,spenso::mink(4,eid(0)))*Q(eid,spenso::mink(4,eid(0)))"]
            A -> D    [num="P(eid,spenso::mink(4,0))"]
            D:1 -> B:2    [id=2]
            B -> ext
            C -> D
            B -> ext
        }

        digraph G2{
            ext [style=invis]
            A [num="1"]
            B [num="1"]
            C [num="1"]
            D [num="1"]
            ext -> A  [pdg=1000]
            ext -> C  [pdg=1000]
            C -> A    [pdg=1000, num="P(eid,spenso::mink(4,eid(0)))*Q(eid,spenso::mink(4,eid(0)))"]
            A -> D    [pdg=1000, num="K(eid,spenso::mink(4,0))"]
            D -> B    [pdg=1000, num="1"]
            B -> ext   [pdg=1000]
            C -> D    [pdg=1000, num="1"]
            B -> ext   [pdg=1000]
        }
        ,
        "scalars"
    )
    .unwrap();

    let g = &graphs[1];

    println!(
        "{}",
        g.underlying
            .dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
    );

    let num = Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter());

    let expr = num.state.expr.as_view();

    let lib: DummyLibrary<SymbolicTensor<Aind>> = DummyLibrary::<_>::new();
    let net =
        Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, _, Aind>::try_from_view(expr, &lib)
            .unwrap();

    println!("{}", expr);
    println!(
        "{}",
        net.dot_display_impl(
            |a| a.to_string(),
            |_| None,
            |a| {
                if let Ok(a) =
                    PermutedStructure::<ShadowedStructure<Aind>>::try_from(a.expression.as_view())
                {
                    a.structure
                        .name()
                        .map(|s| {
                            if let Some(a) = a.structure.args() {
                                FunctionBuilder::new(s).add_args(&a).finish().to_string()
                            } else {
                                s.to_string()
                            }
                        })
                        .unwrap_or("".to_string())
                } else {
                    "".to_string()
                }
            }
        )
    );

    let num = num.color_simplify();

    println!("{}", num.state.expr);
    // let num = num.gamma_simplify();

    // println!(
    //     "{}",
    //     num.state
    //         .colorless
    //         .get_owned_linear(FlatIndex::from(0))
    //         .unwrap()
    //         .to_dots()
    // );
}

#[test]
fn parse_lmbsetting_crossection() {
    test_initialise().unwrap();
    let g: Graph = dot!(
        digraph G{
            num = "1/2"
            overall_factor = "2*x"
            projector = "spenso::t"

            edge [num="1" pdg=1000]
            ext1 [style=invis is_cut=2]
            ext2 [style=invis is_cut=1]
            A [num="1" color_num="a"]
            B [num="1"]
            C [num="1"]
            D [num="1"]
            E [num="1"]
            ext1 -> A
            ext2 -> C
            C -> A    [ num="P(eid,spenso::mink(4,0))"]
            A -> D    [ num="P(eid,spenso::mink(4,0))"]
            D -> B    [ num="1"]
            B -> ext1
            E -> B    [lmb_id=0]
            E -> A    [lmb_id=1]
            E -> C
            C -> D    [ num="1"]
            B -> ext2
        },"scalars"
    )
    .unwrap();

    assert_snapshot!(g.debug_dot(),@r#"
    digraph G{
      num = "1/2";
      overall_factor = "2*x";
      projector = "spenso::t";
      A[dod=0 num="1"];
      B[dod=0 num="1"];
      C[dod=0 num="1"];
      D[dod=0 num="1"];
      E[dod=0 num="1"];

      B:3	-> C:0	 [id=0 dir=none source=1 sink=3  dod=-2 is_cut=0 is_dummy=false lmb_rep="P(0,a___)" name=e0 num="1" particle="scalar_0"];
      B:2	-> A:1	 [id=1 dir=none source=0 sink=3  dod=-2 is_cut=1 is_dummy=false lmb_rep="P(1,a___)" name=e1 num="1" particle="scalar_0"];
      E:14	-> C:15	 [id=2 dir=none source=2 sink=2  dod=-2 is_dummy=false lmb_rep="-P(0,a___)+K(0,a___)+K(1,a___)" name=e2 num="1" particle="scalar_0"];
      C:4	-> A:5	 [id=3 dir=none source=0 sink=1  dod=-2 is_dummy=false lmb_id=0 lmb_rep="K(0,a___)" name=e3 num="P(3,spenso::mink(4,0))" particle="scalar_0"];
      C:6	-> D:7	 [id=4 dir=none source=1 sink=1  dod=-2 is_dummy=false lmb_id=1 lmb_rep="K(1,a___)" name=e4 num="1" particle="scalar_0"];
      D:8	-> B:9	 [id=5 dir=none source=2 sink=2  dod=-2 is_dummy=false lmb_rep="K(1,a___)+K(2,a___)" name=e5 num="1" particle="scalar_0"];
      E:10	-> A:11	 [id=6 dir=none source=0 sink=2  dod=-2 is_dummy=false lmb_rep="-P(1,a___)-K(0,a___)+K(2,a___)" name=e6 num="1" particle="scalar_0"];
      E:12	-> B:13	 [id=7 dir=none source=1 sink=3  dod=-2 is_dummy=false lmb_rep="P(0,a___)+P(1,a___)-K(1,a___)-K(2,a___)" name=e7 num="1" particle="scalar_0"];
      A:16	-> D:17	 [id=8 dir=none source=0 sink=0  dod=-2 is_dummy=false lmb_id=2 lmb_rep="K(2,a___)" name=e8 num="P(0,spenso::mink(4,0))" particle="scalar_0"];
    }
    "#);
}
#[test]
fn parse_lmbsetting() {
    test_initialise().unwrap();
    let g: Graph = dot!(
        digraph G{
            graph [
                multiplicity_factor = "1/2"
                overall_factor = "2*x"
                color = "spenso::t"
                colorless = "spenso::gamma"
            ]
            ext [style=invis]
            A [num="1" color_num="a"]
            B [num="1"]
            C [num="1"]
            D [num="1"]
            E [num="1"]
            ext -> A  [pdg=1000]
            ext -> C  [pdg=1000]
            C -> A    [pdg=1000, num="P(eid,spenso::mink(4,0))"]
            A -> D    [pdg=1000, num="P(eid,spenso::mink(4,0))"]
            D -> B    [pdg=1000, num="1"]
            B -> ext   [pdg=1000]
            E -> B    [pdg=1000,lmb_id=0]
            E -> A    [pdg=1000,lmb_id=1]
            E -> C   [pdg=1000]
            C -> D    [pdg=1000, num="1"]
            B -> ext   [pdg=1000]
        },"scalars"
    )
    .unwrap();

    assert_snapshot!(g.debug_dot(),@r#"
    digraph G{
      num = "1";
      overall_factor = "2*x";
      projector = "1";
      A[dod=0 num="1"];
      B[dod=0 num="1"];
      C[dod=0 num="1"];
      D[dod=0 num="1"];
      E[dod=0 num="1"];

      A:0	-> D:1	 [id=0 dir=none source=0 sink=0  dod=-2 is_dummy=false lmb_id=0 lmb_rep="K(0,a___)" name=e0 num="P(0,spenso::mink(4,0))" particle="scalar_0"];
      ext1	 [style=invis];
      B:2	-> ext1	 [id=1 dir=none source=0 dod=-2 is_dummy=false lmb_rep="-P(1,a___)+P(2,a___)+P(3,a___)" name=e1 num="1𝑖" particle="scalar_0"];
      ext2	 [style=invis];
      B:3	-> ext2	 [id=2 dir=none source=1 dod=-2 is_dummy=false lmb_rep="P(1,a___)" name=e2 num="1𝑖" particle="scalar_0"];
      C:4	-> A:5	 [id=3 dir=none source=0 sink=1  dod=-2 is_dummy=false lmb_id=1 lmb_rep="K(1,a___)" name=e3 num="P(3,spenso::mink(4,0))" particle="scalar_0"];
      C:6	-> D:7	 [id=4 dir=none source=1 sink=1  dod=-2 is_dummy=false lmb_rep="P(2,a___)+P(3,a___)-K(0,a___)-K(2,a___)" name=e4 num="1" particle="scalar_0"];
      D:8	-> B:9	 [id=5 dir=none source=2 sink=2  dod=-2 is_dummy=false lmb_rep="P(2,a___)+P(3,a___)-K(2,a___)" name=e5 num="1" particle="scalar_0"];
      E:10	-> A:11	 [id=6 dir=none source=0 sink=2  dod=-2 is_dummy=false lmb_rep="-P(2,a___)+K(0,a___)-K(1,a___)" name=e6 num="1" particle="scalar_0"];
      E:12	-> B:13	 [id=7 dir=none source=1 sink=3  dod=-2 is_dummy=false lmb_id=2 lmb_rep="K(2,a___)" name=e7 num="1" particle="scalar_0"];
      E:14	-> C:15	 [id=8 dir=none source=2 sink=2  dod=-2 is_dummy=false lmb_rep="P(2,a___)-K(0,a___)+K(1,a___)-K(2,a___)" name=e8 num="1" particle="scalar_0"];
      ext9	 [style=invis];
      ext9	-> A:16	 [id=9 dir=none sink=3 dod=-2 is_dummy=false lmb_rep="P(2,a___)" name=e9 num="1𝑖" particle="scalar_0"];
      ext10	 [style=invis];
      ext10	-> C:17	 [id=10 dir=none sink=3 dod=-2 is_dummy=false lmb_rep="P(3,a___)" name=e10 num="1𝑖" particle="scalar_0"];
    }
    "#);

    println!(
        "{}",
        g.underlying
            .dot_lmb(&g.underlying.full(), &g.loop_momentum_basis)
    );

    let num = Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter());

    let expr = num.state.expr.as_view();

    let lib: DummyLibrary<SymbolicTensor<Aind>> = DummyLibrary::<_>::new();
    let net =
        Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, _, Aind>::try_from_view(expr, &lib)
            .unwrap();

    println!("{}", expr);
    println!(
        "{}",
        net.dot_display_impl(
            |a| a.to_string(),
            |_| None,
            |a| {
                if let Ok(a) =
                    PermutedStructure::<ShadowedStructure<Aind>>::try_from(a.expression.as_view())
                {
                    a.structure
                        .name()
                        .map(|s| {
                            if let Some(a) = a.structure.args() {
                                FunctionBuilder::new(s).add_args(&a).finish().to_string()
                            } else {
                                s.to_string()
                            }
                        })
                        .unwrap_or("".to_string())
                } else {
                    "".to_string()
                }
            }
        )
    );

    let num = num.color_simplify();

    println!("{}", num.state.expr);
    let num = num.gamma_simplify();

    println!("{}", num.state.expr);
}

#[test]
fn parse_and_build_edgevec() {
    test_initialise().unwrap();
    let graph_1: Graph = dot!(
                digraph triangle {
        graph [
            overall_factor = 1;
            multiplicity_factor = 1;
        ]
        edge [
            pdg=1000
        ]
        ext [style=invis]
        ext -> v4 [id=0]
        ext -> v5 [id=1]
        v6 -> ext
        v5 -> v4 [lmb_index=0];
        v6 -> v5;
        v4 -> v6 ;
    }, "scalars"
    )
    .unwrap();

    let graph_2: Graph = dot! (
                digraph triangle {
        graph [
            overall_factor = 1;
            multiplicity_factor = 1;
        ]
        edge [
            pdg=1000
        ]
        ext [style=invis]
        ext -> v4 [id=0]
        ext -> v5 [id=1]
        v6 -> ext [id=2]
        v5 -> v4 [lmb_index=0];
        v6 -> v5;
        v4 -> v6 ;
    }, "scalars"
    )
    .unwrap();

    let graph_1 = &graph_1.underlying;
    let graph_2 = &graph_2.underlying;

    let test_data = vec![true, true, true];

    let graph_1_test = graph_1.new_edgevec(|_, _, p| p.is_paired());

    let graph_2_test = graph_2.new_edgevec(|_, _, p| p.is_paired());

    let mut graph_1_iter = test_data.clone().into_iter();
    let mut graph_2_iter = test_data.clone().into_iter();

    let graph_1_test_2 = graph_1
        .new_edgevec_from_iter(graph_1.iter_edges().map(|(pair, _, _)| {
            if matches!(pair, HedgePair::Paired { .. }) {
                graph_1_iter.next().unwrap()
            } else {
                false
            }
        }))
        .unwrap();

    let graph_2_test_2 = graph_2
        .new_edgevec_from_iter(graph_2.iter_edges().map(|(pair, _, _)| {
            if matches!(pair, HedgePair::Paired { .. }) {
                graph_2_iter.next().unwrap()
            } else {
                false
            }
        }))
        .unwrap();

    assert_eq!(graph_1_test, graph_1_test_2);
    assert_eq!(graph_2_test, graph_2_test_2);
}

#[test]
fn test_group_parsing_1() {
    test_initialise().unwrap();
    let mut graphs: Vec<Graph> = dot! (
        digraph G1{
            graph [
                group_id=0
                is_group_master=true
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G2{
            graph [
                group_id=0
                is_group_master=false
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G3{
            graph [
                group_id=1
                is_group_master=true
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G4{
            graph [
                group_id=1
                is_group_master=false
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }


    )
    .unwrap();

    let groups = complete_group_parsing(&mut graphs).unwrap();
    let expected_result = ti_vec![
        GraphGroup {
            master: 0,
            remaining: vec![1],
        },
        GraphGroup {
            master: 2,
            remaining: vec![3],
        },
    ];

    assert_eq!(groups, expected_result);
}

#[test]
fn test_group_parsing_2() {
    test_initialise().unwrap();
    let mut graphs: Vec<Graph> = dot! (
        digraph G1{
            graph [
                group_id=0
                is_group_master=true
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G2{
            graph [
                group_id=0
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G3{
            graph [
                group_id=1

            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G4{
            graph [
                group_id=1
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }


    )
    .unwrap();

    let groups = complete_group_parsing(&mut graphs).unwrap();
    let expected_result = ti_vec![
        GraphGroup {
            master: 0,
            remaining: vec![1],
        },
        GraphGroup {
            master: 2,
            remaining: vec![3],
        },
    ];

    assert_eq!(groups, expected_result);
}

#[test]
fn test_group_parsing_3() {
    test_initialise().unwrap();
    let mut graphs: Vec<Graph> = dot! (
        digraph G1{
            graph [
                group_id=0
                is_group_master=true
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G2{
            graph [
                group_id=0
                is_group_master=false
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G3{
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G4{
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }


    )
    .unwrap();

    let groups = complete_group_parsing(&mut graphs).unwrap();
    let expected_result = ti_vec![
        GraphGroup {
            master: 0,
            remaining: vec![1],
        },
        GraphGroup {
            master: 2,
            remaining: vec![],
        },
        GraphGroup {
            master: 3,
            remaining: vec![],
        },
    ];

    assert_eq!(groups, expected_result);
}

#[test]
fn test_group_parsing_4() {
    test_initialise().unwrap();
    let mut graphs: Vec<Graph> = dot! (
        digraph G1{
            graph [
                group_id=0
                is_group_master=true
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G2{
            graph [
                group_id=0
                is_group_master=true
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G3{
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G4{
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }


    )
    .unwrap();

    let groups = complete_group_parsing(&mut graphs);
    assert!(groups.is_err());
}

#[test]
fn test_group_parsing_5() {
    test_initialise().unwrap();
    let mut graphs: Vec<Graph> = dot! (
        digraph G1{
            graph [
                group_id=1
                is_group_master=true
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G2{
            graph [
                group_id=1
                is_group_master=false
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G3{
            graph [
                group_id=2
                is_group_master=true
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }

        digraph G4{
            graph [
                group_id=2
                is_group_master=false
            ]
            ext    [style=invis]
            node [num=1]
            ext -> A
            A -> B
            A -> B
            B -> ext
        }


    )
    .unwrap();

    assert!(complete_group_parsing(&mut graphs).is_err());
}

#[test]
fn parse_triangle_lmb() {
    test_initialise().unwrap();

    let g: Graph = dot!(
        digraph massless_triangle_0 {
            ext [style=invis]
        ext -> v4:0 [pdg=1000, id=0 ];
        v5:1 -> ext [pdg=1000, id=1 ];
        v6:2 -> ext [pdg=1000, id=2 ];
        v4 -> v5 [pdg=1000, id=3, lmb_id=0];
        v5 -> v6 [pdg=1000, id=4];
        v6 -> v4 [pdg=1000, id=5];
    }
    ,"scalars"
        )
    .unwrap();

    let lmb = g.loop_momentum_basis.loop_edges;
    assert_eq!(lmb[LoopIndex::from(0)], EdgeIndex::from(3));

    let g: Graph = dot!(
        digraph massless_triangle_0 {
            ext [style=invis]
        ext -> v4:0 [pdg=1000, id=0 ];
        v5:1 -> ext [pdg=1000, id=1 ];
        v6:2 -> ext [pdg=1000, id=2 ];
        v4 -> v5 [pdg=1000, id=3];
        v5 -> v6 [pdg=1000, id=4, lmb_id=0];
        v6 -> v4 [pdg=1000, id=5];
    }
    ,"scalars"
        )
    .unwrap();

    let lmb = g.loop_momentum_basis.loop_edges;
    assert_eq!(lmb[LoopIndex::from(0)], EdgeIndex::from(4));

    let g: Graph = dot!(
        digraph massless_triangle_0 {
            ext [style=invis]
        ext -> v4:0 [pdg=1000, id=0 ];
        v5:1 -> ext [pdg=1000, id=1 ];
        v6:2 -> ext [pdg=1000, id=2 ];
        v4 -> v5 [pdg=1000, id=3];
        v5 -> v6 [pdg=1000, id=4];
        v6 -> v4 [pdg=1000, id=5, lmb_id=0];
    }
    ,"scalars"
        )
    .unwrap();

    let lmb = g.loop_momentum_basis.loop_edges;
    assert_eq!(lmb[LoopIndex::from(0)], EdgeIndex::from(5));
}
