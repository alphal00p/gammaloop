use insta::assert_snapshot;
use linnet::half_edge::{
    NodeIndex,
    involution::{EdgeIndex, HedgePair},
};
use spenso::{
    network::{
        Network,
        library::DummyLibrary,
        parsing::{ParseSettings, ShadowedStructure},
        store::NetworkStore,
    },
    shadowing::symbolica_utils::AtomCoreExt,
    structure::{HasName, PermutedStructure},
    tensors::symbolic::SymbolicTensor,
};
use symbolica::atom::{Atom, FunctionBuilder, Symbol};
use tracing::info;
use typed_index_collections::ti_vec;

use super::Graph;
use crate::{
    dot,
    graph::{
        GraphGroup, LMBext,
        parse::{IntoGraph, complete_group_parsing, string_utils::ToOrderedSimple},
    },
    initialisation::test_initialise,
    momentum::sample::LoopIndex,
    numerator::{Numerator, UnInit, aind::Aind},
    processes::DotExportSettings,
};

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

    println!("{}", graph[0].dot_serialize(&DotExportSettings::default()));

    println!(
        "{}",
        graph[0].dot_lmb_of(&graph[0].underlying.full(), &graph[0].loop_momentum_basis)
    );

    println!("{}", graph[0].loop_momentum_basis);

    println!(
        "{}",
        graph[1].dot_lmb_of(&graph[1].underlying.full(), &graph[1].loop_momentum_basis)
    );
    println!("{}", graph[1].loop_momentum_basis);
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
        g.dot_lmb_of(&g.underlying.full(), &g.loop_momentum_basis)
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

    println!("{}", g.dot_serialize(&DotExportSettings::default()))
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
            .dot_lmb_of(&g.underlying.full(), &g.loop_momentum_basis)
    );

    println!("{}", g.dot_serialize(&DotExportSettings::default()));
    // return;
    let num = Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter());

    let expr = num.state.expr.as_view();

    let lib: DummyLibrary<SymbolicTensor<Aind>> = DummyLibrary::<_>::new();
    let net =
        Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, _, Symbol, Aind>::try_from_view::<
            SymbolicTensor<Aind>,
            _,
        >(expr, &lib, &ParseSettings::default())
        .unwrap();

    println!("{}", expr);
    println!(
        "{}",
        net.dot_display_impl(
            ToString::to_string,
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
            },
            ToString::to_string,
        )
    );

    let num = num.color_simplify();

    println!("{}", num.state.expr);
    // println!("{}", num.state.expr.to_dots());
    // let num = num.gamma_simplify();

    // println!("{}", num.state.expr);
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
      overall_factor_evaluated = "2*x";
      projector = "spenso::t";
      A[num="1"];
      B[num="1"];
      C[num="1"];
      D[num="1"];
      E[num="1"];

      B:3	-> C:0	 [id=0 dir=none  is_cut="0" lmb_rep="P(0,a___)" num="1" particle="scalar_0"];
      B:2	-> A:1	 [id=1 dir=none  is_cut="1" lmb_rep="P(1,a___)" num="1" particle="scalar_0"];
      A:17	-> D:16	 [id=2 dir=none  lmb_rep="K(1,a___)+K(2,a___)+P(1,a___)" num="P(0,spenso::mink(4,0))" particle="scalar_0"];
      C:4	-> A:5	 [id=3 dir=none  lmb_id="2" lmb_rep="K(2,a___)" num="P(3,spenso::mink(4,0))" particle="scalar_0"];
      C:6	-> D:7	 [id=4 dir=none  lmb_rep="-1*K(0,a___)+-1*K(1,a___)+-1*K(2,a___)+P(0,a___)" num="1" particle="scalar_0"];
      D:8	-> B:9	 [id=5 dir=none  lmb_rep="-1*K(0,a___)+P(0,a___)+P(1,a___)" num="1" particle="scalar_0"];
      E:10	-> A:11	 [id=6 dir=none  lmb_id="1" lmb_rep="K(1,a___)" num="1" particle="scalar_0"];
      E:12	-> B:13	 [id=7 dir=none  lmb_id="0" lmb_rep="K(0,a___)" num="1" particle="scalar_0"];
      E:14	-> C:15	 [id=8 dir=none  lmb_rep="-1*K(0,a___)+-1*K(1,a___)" num="1" particle="scalar_0"];
    }
    "#);
}

#[test]
fn dod_rescales_only_internal_edge_qs() {
    test_initialise().unwrap();
    let g: Graph = dot!(
        digraph G {
            edge [pdg=1000]
            ext0 [style=invis]
            ext1 [style=invis]
            A [num="1"]
            B [num="1"]
            ext0 -> A [id=0]
            B -> ext1 [id=1]
            A -> B [id=2, num="Q(2,spenso::mink(4,edge(2,1)))*Q(0,spenso::mink(4,edge(2,1)))"]
            A -> B [id=3, num="1"]
        },
        "scalars"
    )
    .unwrap();

    assert_eq!(g.underlying[EdgeIndex::from(2)].dod.value, -1);
}

#[test]
fn lmb() {
    test_initialise().unwrap();

    let _g: Vec<Graph> = dot!(
        digraph G{
            edge [num="1"]
            ext    [style=invis]
            node [num="1"]
            ext -> A
            ext -> A
            ext -> A
            A -> B
            B -> ext
            B -> ext
            B -> ext
            B -> ext
        }

        digraph G{
            edge [num="1"]
            ext    [style=invis]
            node [num="1"]
            A ->C
            C -> D
            D -> A
            A -> B

            B->B1
            B1->B2
            B2->B
        }
    )
    .unwrap();
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

mod failing {
    use super::*;

    #[test]
    fn xs_parsing() {
        test_initialise().unwrap();
        let g: Graph = dot!(
            digraph{
              ext0	 [style=invis is_cut=0];
              ext1	 [style=invis is_cut=1];
              ext0	-> 0	 [dir=back particle="e+"];
              ext1	-> 0	 [particle="e-"];
                0	-> 1	 [particle="a"];
                1 -> 2	 [particle="d"];
                2 -> 4   [particle="g"];
                2 -> 3	 [particle="d"];
                3 -> 4	 [particle="d"];
                4 -> 1	 [particle="d"];
                3->5	 [particle="a"];
                5	-> ext0	 [dir=back particle="e+"];
                5	-> ext1	 [particle="e-"];
            }
        )
        .unwrap();
        assert_snapshot!(g.debug_dot(),@r#"
        digraph {
          num = "1";
          overall_factor = "1";
          overall_factor_evaluated = "1";
          projector = "u(1,spenso::bis(4,hedge(17)))*ubar(1,spenso::bis(4,hedge(15)))*v(0,spenso::bis(4,hedge(14)))*vbar(0,spenso::bis(4,hedge(16)))";
          0[dod="0" int_id="V_98" num="UFO::GC_3*spenso::gamma(spenso::bis(4,hedge(16)),spenso::bis(4,hedge(17)),spenso::mink(4,hedge(0)))"];
          1[dod="0" int_id="V_71" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(13)),spenso::dind(spenso::cof(3,hedge(2))))*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(13)),spenso::mink(4,hedge(1)))"];
          2[dod="0" int_id="V_74" num="UFO::GC_11*spenso::gamma(spenso::bis(4,hedge(4)),spenso::bis(4,hedge(3)),spenso::mink(4,hedge(6)))*spenso::t(spenso::coad(8,hedge(6)),spenso::cof(3,hedge(3)),spenso::dind(spenso::cof(3,hedge(4))))"];
          3[dod="0" int_id="V_71" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(5)),spenso::dind(spenso::cof(3,hedge(8))))*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(5)),spenso::mink(4,hedge(10)))"];
          4[dod="0" int_id="V_74" num="UFO::GC_11*spenso::gamma(spenso::bis(4,hedge(12)),spenso::bis(4,hedge(9)),spenso::mink(4,hedge(7)))*spenso::t(spenso::coad(8,hedge(7)),spenso::cof(3,hedge(9)),spenso::dind(spenso::cof(3,hedge(12))))"];
          5[dod="0" int_id="V_98" num="UFO::GC_3*spenso::gamma(spenso::bis(4,hedge(15)),spenso::bis(4,hedge(14)),spenso::mink(4,hedge(11)))"];

          5:14	-> 0:0	 [id=0 dir=back source="{ufo_order:1}" sink="{ufo_order:0}"  dod="-2" is_cut="0"  lmb_rep="P(0,a___)" name="e0" num="1" particle="e+"];
          5:15	-> 0:1	 [id=1 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-2" is_cut="1"  lmb_rep="P(1,a___)" name="e1" num="1" particle="e-"];
          2:4	-> 3:5	 [id=2 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_rep="-1*K(1,a___)+K(0,a___)" name="e2" num="Q(2,spenso::mink(4,edge(2,1)))*spenso::g(spenso::cof(3,hedge(4)),spenso::dind(spenso::cof(3,hedge(5))))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(4)),spenso::mink(4,edge(2,1)))" particle="d"];
          2:6	-> 4:7	 [id=3 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="-2"  lmb_id="1" lmb_rep="K(1,a___)" name="e3" num="spenso::g(spenso::coad(8,hedge(6)),spenso::coad(8,hedge(7)))*spenso::g(spenso::mink(4,hedge(6)),spenso::mink(4,hedge(7)))" particle="g"];
          3:8	-> 4:9	 [id=4 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_rep="-1*K(1,a___)+-1*P(0,a___)+-1*P(1,a___)+K(0,a___)" name="e4" num="Q(4,spenso::mink(4,edge(4,1)))*spenso::g(spenso::cof(3,hedge(8)),spenso::dind(spenso::cof(3,hedge(9))))*spenso::gamma(spenso::bis(4,hedge(9)),spenso::bis(4,hedge(8)),spenso::mink(4,edge(4,1)))" particle="d"];
          3:10	-> 5:11	 [id=5 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="-2"  lmb_rep="P(0,a___)+P(1,a___)" name="e5" num="-1*spenso::g(spenso::mink(4,hedge(10)),spenso::mink(4,hedge(11)))" particle="a"];
          4:12	-> 1:13	 [id=6 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_rep="-1*P(0,a___)+-1*P(1,a___)+K(0,a___)" name="e6" num="Q(6,spenso::mink(4,edge(6,1)))*spenso::g(spenso::cof(3,hedge(12)),spenso::dind(spenso::cof(3,hedge(13))))*spenso::gamma(spenso::bis(4,hedge(13)),spenso::bis(4,hedge(12)),spenso::mink(4,edge(6,1)))" particle="d"];
          0:16	-> 1:17	 [id=7 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="-2"  lmb_rep="P(0,a___)+P(1,a___)" name="e7" num="-1*spenso::g(spenso::mink(4,hedge(0)),spenso::mink(4,hedge(1)))" particle="a"];
          1:2	-> 2:3	 [id=8 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_id="0" lmb_rep="K(0,a___)" name="e8" num="Q(8,spenso::mink(4,edge(8,1)))*spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(3))))*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(2)),spenso::mink(4,edge(8,1)))" particle="d"];
        }
        "#);

        let g: Graph = dot!(
            digraph{
                edge [
                    pdg=1000
                ]
                ext0[style=invis is_cut=0];
                ext0->a
                a->b->c->d->a
                b->e
                c->e
                d->e
                c->ext0


            }
            ,"scalars")
        .unwrap();
        assert_snapshot!(g.debug_dot(),@r#"
        digraph {
          num = "1";
          overall_factor = "1";
          overall_factor_evaluated = "1";
          projector = "1";
          a[dod="0" int_id="V_3_SCALAR_000" num="UFO::SCALAR_COUPLING"];
          b[dod="0" int_id="V_3_SCALAR_000" num="UFO::SCALAR_COUPLING"];
          c[dod="0" int_id="V_4_SCALAR_0000" num="UFO::SCALAR_COUPLING"];
          d[dod="0" int_id="V_3_SCALAR_000" num="UFO::SCALAR_COUPLING"];
          e[dod="0" int_id="V_3_SCALAR_000" num="UFO::SCALAR_COUPLING"];

          c:10	-> a:0	 [id=0 dir=none source="{ufo_order:3}" sink="{ufo_order:2}"  dod="-2" is_cut="0"  lmb_rep="P(0,a___)" name="e0" num="1" particle="scalar_0"];
          b:2	-> c:3	 [id=1 dir=none source="{ufo_order:1}" sink="{ufo_order:0}"  dod="-2"  lmb_id="0" lmb_rep="K(0,a___)" name="e1" num="1" particle="scalar_0"];
          b:4	-> e:5	 [id=2 dir=none source="{ufo_order:2}" sink="{ufo_order:0}"  dod="-2"  lmb_rep="-1*K(0,a___)+K(2,a___)+P(0,a___)" name="e2" num="1" particle="scalar_0"];
          c:6	-> d:7	 [id=3 dir=none source="{ufo_order:1}" sink="{ufo_order:0}"  dod="-2"  lmb_rep="-1*K(1,a___)+-1*P(0,a___)+K(0,a___)" name="e3" num="1" particle="scalar_0"];
          c:8	-> e:9	 [id=4 dir=none source="{ufo_order:2}" sink="{ufo_order:1}"  dod="-2"  lmb_id="1" lmb_rep="K(1,a___)" name="e4" num="1" particle="scalar_0"];
          a:15	-> b:1	 [id=5 dir=none source="{ufo_order:0}" sink="{ufo_order:0}"  dod="-2"  lmb_rep="K(2,a___)+P(0,a___)" name="e5" num="1" particle="scalar_0"];
          d:11	-> a:12	 [id=6 dir=none source="{ufo_order:1}" sink="{ufo_order:1}"  dod="-2"  lmb_id="2" lmb_rep="K(2,a___)" name="e6" num="1" particle="scalar_0"];
          d:13	-> e:14	 [id=7 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="-2"  lmb_rep="-1*K(1,a___)+-1*K(2,a___)+-1*P(0,a___)+K(0,a___)" name="e7" num="1" particle="scalar_0"];
        }
        "#);

        let g: Graph = dot!(
            digraph GL06{
              num = "1";
              overall_factor = "(AutG(1))^(-1)*AntiFermionSpinSumSign(-1)*ExternalFermionOrderingSign(-1)*InternalFermionLoopSign(-1)*NumeratorIndependentSymmetryGrouping(2)";
              0[int_id=V_141];
              1[int_id=V_141];
              2[int_id=V_134];
              3[int_id=V_134];
              4[int_id=V_137];
              5[int_id=V_137];
              6[int_id=V_98];
              7[int_id=V_98];
              ext0 [style=invis];
              ext0-> 7:0 [id=0 dir=back is_cut=0 particle="e+"];
              ext1 [style=invis];
              ext1-> 7:1 [id=1 is_cut=1 particle="e-"];
              0:2-> 1:3 [id=2  particle="t"];
              0:4-> 1:5 [id=3 dir=none  particle="H"];
              5:6-> 0:7 [id=4  particle="t"];
              1:8-> 3:9 [id=5  particle="t"];
              4:10-> 2:11 [id=6  particle="t"];
              2:12-> 5:13 [id=7  particle="t"];
              2:14-> 7:15 [id=8 dir=none  particle="a"];
              3:16-> 4:17 [id=9  particle="t"];
              3:18-> 6:19 [id=10 dir=none  particle="a"];
              4:20-> 5:21 [id=11 dir=none  particle="g"];
              ext12 [style=invis];
              6:22-> ext12 [id=12 dir=back is_cut=0 particle="e+"];
              ext13 [style=invis];
              6:23-> ext13 [id=13 is_cut=1 particle="e-"];
            }
        )
        .unwrap();

        println!("{}", g.global_prefactor.projector);

        let g: Graph = dot!(
        digraph GL06{
            num = "-2";
                exte0    [style=invis, is_cut=0];
                exte0   -> 7:0   [id=0 particle="e+"];
                exte1    [style=invis, is_cut=1];
                exte1   -> 7:1   [id=1 particle="e-"];
                4:22    -> 5:23  [id=2 particle="g"];
                3:20    -> 6:21  [id=3 particle="a"];
                0:4     -> 1:5   [id=4 particle="t" lmb_id=0];
                0:6     -> 1:7   [id=5 particle="H" lmb_id=1];
                5:8     -> 0:9   [id=6 particle="t"];
                1:10    -> 3:11  [id=7 particle="t"];
                4:12    -> 2:13  [id=8 particle="t"];
                2:14    -> 5:15  [id=9 particle="t" lmb_id=2];
                2:16    -> 7:17  [id=10 particle="a"];
                3:18    -> 4:19  [id=11 particle="t"];

                6:3     -> exte0        [id=12 particle="e+"];

                6:2     -> exte1        [id=13 particle="e-"];
        })
        .unwrap();

        println!("{}", g.global_prefactor.projector);
    }

    #[test]
    fn massive_gluon() {
        test_initialise().unwrap();
        let _g: Graph = dot!(
        digraph ddxaaapentagon{
            // num = "-2";
            ext   [style=invis];
            ext	-> 0:0	 [id=0 particle="d"];
            ext	-> 1:1	 [id=1 particle="d~"];
            2:2	-> ext	 [id=2 particle="a"];
            3:3	-> ext	 [id=3 particle="a"];
            4:4	-> ext   [id=4 particle="a"];

            0	-> 2	 [id=5 particle="d"];
            2	-> 3	 [id=6 particle="d"];
            3	-> 4	 [id=7 particle="d"];
            4 -> 1     [id=8 particle="d"];
            // 1 -> 0     [id=9 lmb_id="0" particle="g"];
            // The option below would introduce an IR regulator
            1 -> 0     [id=9 lmb_id="0" particle="g" mass="1"];
        })
        .unwrap();
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
          overall_factor_evaluated = "1";
          projector = "u(0,spenso::bis(4,hedge(1)))*u(1,spenso::bis(4,hedge(3)))*ubar(0,spenso::bis(4,hedge(0)))*ubar(1,spenso::bis(4,hedge(2)))";
          0[dod="0" int_id="V_89" num="UFO::GC_100*spenso::g(spenso::cof(3,hedge(1)),spenso::dind(spenso::cof(3,hedge(4))))*spenso::gamma(spenso::bis(4,hedge(4)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(8)))*spenso::projm(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(1)))"];
          1[dod="0" int_id="V_127" num="UFO::GC_45*spenso::g(spenso::cof(3,hedge(3)),spenso::dind(spenso::cof(3,hedge(6))))*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(9)))*spenso::projm(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3)))"];
          2[dod="0" int_id="V_123" num="UFO::GC_41*spenso::g(spenso::cof(3,hedge(5)),spenso::dind(spenso::cof(3,hedge(0))))*spenso::gamma(spenso::bis(4,hedge(0)),spenso::bis(4,vertex(2,1)),spenso::mink(4,hedge(10)))*spenso::projm(spenso::bis(4,vertex(2,1)),spenso::bis(4,hedge(5)))"];
          3[dod="0" int_id="V_93" num="UFO::GC_104*spenso::g(spenso::cof(3,hedge(7)),spenso::dind(spenso::cof(3,hedge(2))))*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(3,1)),spenso::mink(4,hedge(11)))*spenso::projm(spenso::bis(4,vertex(3,1)),spenso::bis(4,hedge(7)))"];

          2:3	-> 0:0	 [id=0 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-2" is_cut="0"  lmb_rep="P(0,a___)" name="e0" num="spenso::g(spenso::cof(3,hedge(0)),spenso::dind(spenso::cof(3,hedge(1))))" particle="d"];
          3:2	-> 1:1	 [id=1 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-2" is_cut="1"  lmb_rep="P(1,a___)" name="e1" num="spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(3))))" particle="c"];
          0:8	-> 1:9	 [id=2 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="0"  lmb_rep="-1*K(0,a___)+P(0,a___)" name="e2" num="-1*spenso::g(spenso::mink(4,hedge(8)),spenso::mink(4,hedge(9)))+Q(2,spenso::mink(4,hedge(8)))*Q(2,spenso::mink(4,hedge(9)))*UFO::MW^(-2)" particle="W-"];
          2:10	-> 3:11	 [id=3 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="0"  lmb_rep="-1*P(0,a___)+K(0,a___)" name="e3" num="-1*spenso::g(spenso::mink(4,hedge(10)),spenso::mink(4,hedge(11)))+Q(3,spenso::mink(4,hedge(10)))*Q(3,spenso::mink(4,hedge(11)))*UFO::MW^(-2)" particle="W+"];
          0:4	-> 2:5	 [id=4 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_id="0" lmb_rep="K(0,a___)" name="e4" num="Q(4,spenso::mink(4,edge(4,1)))*spenso::g(spenso::cof(3,hedge(4)),spenso::dind(spenso::cof(3,hedge(5))))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(4)),spenso::mink(4,edge(4,1)))" particle="u"];
          1:6	-> 3:7	 [id=5 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_rep="-1*K(0,a___)+P(0,a___)+P(1,a___)" name="e5" num="Q(5,spenso::mink(4,edge(5,1)))*spenso::g(spenso::cof(3,hedge(6)),spenso::dind(spenso::cof(3,hedge(7))))*spenso::gamma(spenso::bis(4,hedge(7)),spenso::bis(4,hedge(6)),spenso::mink(4,edge(5,1)))" particle="s"];
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
          overall_factor_evaluated = "1";
          projector = "1";
          a[dod="0" int_id="V_90" num="UFO::GC_103*spenso::g(spenso::cof(3,hedge(3)),spenso::dind(spenso::cof(3,hedge(0))))*spenso::gamma(spenso::bis(4,hedge(0)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(5)))*spenso::projm(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(3)))"];
          b[dod="0" int_id="V_126" num="UFO::GC_44*spenso::g(spenso::cof(3,hedge(1)),spenso::dind(spenso::cof(3,hedge(2))))*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(4)))*spenso::projm(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(1)))"];

          a:0	-> b:1	 [id=0 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_id="0" lmb_rep="K(0,a___)" name="e0" num="(Q(0,spenso::mink(4,edge(0,1)))*spenso::gamma(spenso::bis(4,hedge(1)),spenso::bis(4,hedge(0)),spenso::mink(4,edge(0,1)))+UFO::MC*spenso::g(spenso::bis(4,hedge(0)),spenso::bis(4,hedge(1))))*spenso::g(spenso::cof(3,hedge(0)),spenso::dind(spenso::cof(3,hedge(1))))" particle="c"];
          b:2	-> a:3	 [id=1 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_id="1" lmb_rep="K(1,a___)" name="e1" num="Q(1,spenso::mink(4,edge(1,1)))*spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(3))))*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(2)),spenso::mink(4,edge(1,1)))" particle="d"];
          b:4	-> a:5	 [id=2 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="0"  lmb_rep="-1*K(1,a___)+K(0,a___)" name="e2" num="-1*spenso::g(spenso::mink(4,hedge(4)),spenso::mink(4,hedge(5)))+Q(2,spenso::mink(4,hedge(4)))*Q(2,spenso::mink(4,hedge(5)))*UFO::MW^(-2)" particle="W+"];
        }
        "#);
    }

    #[test]
    fn polarizations() {
        test_initialise().unwrap();
        let gs: Vec<Graph> = dot!(
            digraph g1{
               ext [style=invis];
                n [num=1]
                ext -> n [particle="d"];
            }
            digraph g2{
               ext [style=invis];
                n [num=1]
                n -> ext [particle="d"];
            }
            digraph g3{
               ext [style=invis];
                n [num=1]
                ext -> n [dir=back particle="d~"];
            }
            digraph g4{
               ext [style=invis];
                n [num=1]
                n -> ext [dir=back particle="d~"];
            }

            digraph g3p{
               ext [style=invis];
                n [num=1]
                ext -> n [ particle="d~"];
            }
            digraph g4p{
               ext [style=invis];
                n [num=1]
                n -> ext [particle="d~"];
            }
        )
        .unwrap();

        insta::assert_snapshot!(gs[0].global_prefactor.projector.to_bare_ordered_string(),@"u(0,bis(4,hedge(0)))");
        insta::assert_snapshot!(gs[1].global_prefactor.projector.to_bare_ordered_string(),@"ubar(0,bis(4,hedge(0)))");
        insta::assert_snapshot!(gs[2].global_prefactor.projector.to_bare_ordered_string(),@"vbar(0,bis(4,hedge(0)))");
        insta::assert_snapshot!(gs[3].global_prefactor.projector.to_bare_ordered_string(),@"v(0,bis(4,hedge(0)))");
        insta::assert_snapshot!(gs[4].global_prefactor.projector.to_bare_ordered_string(),@"vbar(0,bis(4,hedge(0)))");
        insta::assert_snapshot!(gs[5].global_prefactor.projector.to_bare_ordered_string(),@"v(0,bis(4,hedge(0)))");
    }

    #[test]
    fn vertex_rules() {
        test_initialise().unwrap();
        let gs: Vec<Graph> = dot!(
            digraph g1{
               ext [style=invis];
                ext -> n [particle="g"];
                ext -> n [particle="g"];
                ext -> n [particle="g"];
            }

            digraph g2{
               ext [style=invis];
                ext -> n [particle="d"];
                ext -> n [particle="d~"];
                ext -> n [particle="g"];
            }

            digraph g3{
               ext [style=invis];
                ext -> n [id = 0 particle="ghG"];
                ext -> n [id = 1 particle="ghG~"];
                ext -> n [id = 2 particle="g"];
            }

            digraph g4{
               ext [style=invis];
                ext -> n [id = 0 particle="ghG"];
                n -> ext [id = 1 particle="ghG"];
                ext -> n [id = 2 particle="g"];
            }

            digraph GL55{
              ext   [style=invis];
              0:0   -> ext   [id=0 dir=none particle="g"];
              0:1   -> ext   [id=1 dir=none particle="ghG~"];
              0:2   -> ext   [id=2 dir=none particle="ghG"];
            }

            digraph GL54{
              ext   [style=invis];
              0:0   -> ext  [id=0 dir=none  particle="g"];
              0:1   -> ext  [id=1 dir=none  particle="ghG"];
              0:2   -> ext  [id=2 dir=none  particle="ghG~"];
            }

        )
        .unwrap();

        insta::assert_snapshot!(gs[0].underlying[NodeIndex(0)].num.to_ordered_simple(),@"(-1*Q(0,mink(4,hedge(1)))*g(mink(4,hedge(0)),mink(4,hedge(2)))+-1*Q(1,mink(4,hedge(2)))*g(mink(4,hedge(0)),mink(4,hedge(1)))+-1*Q(2,mink(4,hedge(0)))*g(mink(4,hedge(1)),mink(4,hedge(2)))+Q(0,mink(4,hedge(2)))*g(mink(4,hedge(0)),mink(4,hedge(1)))+Q(1,mink(4,hedge(0)))*g(mink(4,hedge(1)),mink(4,hedge(2)))+Q(2,mink(4,hedge(1)))*g(mink(4,hedge(0)),mink(4,hedge(2))))*GC_10*f(coad(8,hedge(0)),coad(8,hedge(1)),coad(8,hedge(2)))");
        insta::assert_snapshot!(gs[1].underlying[NodeIndex(0)].num.to_ordered_simple(),@"GC_11*gamma(bis(4,hedge(1)),bis(4,hedge(0)),mink(4,hedge(2)))*t(coad(8,hedge(2)),cof(3,hedge(0)),dind(cof(3,hedge(1))))");
        insta::assert_snapshot!(gs[2].underlying[NodeIndex(0)].num.to_ordered_simple(),@"-1*GC_10*Q(1,mink(4,hedge(2)))*f(coad(8,hedge(2)),coad(8,hedge(0)),coad(8,hedge(1)))");
        insta::assert_snapshot!(gs[3].underlying[NodeIndex(0)].num.to_ordered_simple(),@"GC_10*Q(1,mink(4,hedge(1)))*f(coad(8,hedge(1)),coad(8,hedge(0)),coad(8,hedge(2)))");
        insta::assert_snapshot!(gs[4].underlying[NodeIndex(0)].num.to_ordered_simple(),@"GC_10*Q(2,mink(4,hedge(0)))*f(coad(8,hedge(0)),coad(8,hedge(1)),coad(8,hedge(2)))");
        insta::assert_snapshot!(gs[5].underlying[NodeIndex(0)].num.to_ordered_simple(),@"GC_10*Q(1,mink(4,hedge(0)))*f(coad(8,hedge(0)),coad(8,hedge(2)),coad(8,hedge(1)))");

        let gs: Vec<Graph> = dot!(
            digraph g1{
               ext [style=invis];
                ext -> n [particle="scalar_0"];
                ext -> n [particle="scalar_0"];
                ext -> n [particle="graviton"];
            },"scalar_gravity"
        )
        .unwrap();
        insta::assert_snapshot!(gs[0].underlying[NodeIndex(0)].num.to_ordered_simple(),@"(-1*Q(0,mink(4,vertex(0,1)))*Q(1,mink(4,vertex(0,1)))*g(mink(4,hedge(2)),mink(4,hedge(2,1)))+Q(0,mink(4,hedge(2)))*Q(1,mink(4,hedge(2,1)))+Q(0,mink(4,hedge(2,1)))*Q(1,mink(4,hedge(2))))*SST+-1*SSTmpart0*g(mink(4,hedge(2)),mink(4,hedge(2,1)))");
    }

    #[test]
    fn propagators() {
        test_initialise().unwrap();
        let gs: Vec<Graph> = dot!(
            digraph g1{
                node [num=1]
                a -> b[particle="g"];
            }

            digraph g2{
                node [num=1]
                a -> b [particle="d"];
            }

            digraph g3{
                node [num=1]
                a -> b [particle="Z"];
            }

        )
        .unwrap();

        insta::assert_snapshot!(gs[0].underlying[EdgeIndex(0)].num.to_ordered_simple(),@"g(coad(8,hedge(0)),coad(8,hedge(1)))*g(mink(4,hedge(0)),mink(4,hedge(1)))");
        insta::assert_snapshot!(gs[1].underlying[EdgeIndex(0)].num.to_ordered_simple(),@"Q(0,mink(4,edge(0,1)))*g(cof(3,hedge(0)),dind(cof(3,hedge(1))))*gamma(bis(4,hedge(1)),bis(4,hedge(0)),mink(4,edge(0,1)))");
        insta::assert_snapshot!(gs[2].underlying[EdgeIndex(0)].num.to_ordered_simple(),@"-1*g(mink(4,hedge(0)),mink(4,hedge(1)))+MZ^(-2)*Q(0,mink(4,hedge(0)))*Q(0,mink(4,hedge(1)))");
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
                2:0 -> ext0 [id=0 dir=back is_cut=0 particle="a"];
                ext1 [style=invis];
                ext1 -> 3:1 [id=1 is_cut=0 particle="a"];
                0:2 -> 1:3 [id=2 particle="c~"];
                0:4 -> 2:5 [id=3 particle="d"];
                0:6 -> 3:7 [id=4 particle="G+"];
                1:8 -> 2:9 [id=5 particle="d~"];
                1:10 -> 3:11 [id=6 particle="G-"];
            }
        ).unwrap();
        assert_snapshot!(g.debug_dot(),@r#"
        digraph GL1{
          num = "1";
          overall_factor = "(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)";
          overall_factor_evaluated = "-1";
          projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵbar(0,spenso::mink(4,hedge(0)))";
          0[dod="0" int_id="V_117" num="UFO::GC_22*spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(4))))*spenso::projm(spenso::bis(4,hedge(4)),spenso::bis(4,hedge(2)))"];
          1[dod="0" int_id="V_82" num="UFO::GC_16*spenso::g(spenso::cof(3,hedge(8)),spenso::dind(spenso::cof(3,hedge(3))))*spenso::projp(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(8)))"];
          2[dod="0" int_id="V_71" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(5)),spenso::dind(spenso::cof(3,hedge(9))))*spenso::gamma(spenso::bis(4,hedge(9)),spenso::bis(4,hedge(5)),spenso::mink(4,hedge(0)))"];
          3[dod="1" int_id="V_11" num="(-1*Q(4,spenso::mink(4,hedge(1)))+Q(1,spenso::mink(4,hedge(1))))*UFO::GC_3"];

          2:1	-> 3:0	 [id=0 dir=none source="{ufo_order:2}" sink="{ufo_order:0}"  dod="-2" is_cut="0"  lmb_rep="P(0,a___)" name="e0" num="1" particle="a"];
          1:10	-> 3:11	 [id=1 dir=none source="{ufo_order:2}" sink="{ufo_order:1}"  dod="-2"  lmb_rep="-1*K(1,a___)+-1*P(0,a___)" name="e1" num="1" particle="G-"];
          0:2	-> 1:3	 [id=2 dir=back source="{ufo_order:1}" sink="{ufo_order:0}"  dod="-1"  lmb_id="0" lmb_rep="K(0,a___)" name="e2" num="(-1*Q(2,spenso::mink(4,edge(2,1)))*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(3)),spenso::mink(4,edge(2,1)))+UFO::MC*spenso::g(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(3))))*spenso::g(spenso::cof(3,hedge(3)),spenso::dind(spenso::cof(3,hedge(2))))" particle="c~"];
          0:4	-> 2:5	 [id=3 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_rep="-1*K(0,a___)+-1*K(1,a___)" name="e3" num="Q(3,spenso::mink(4,edge(3,1)))*spenso::g(spenso::cof(3,hedge(4)),spenso::dind(spenso::cof(3,hedge(5))))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(4)),spenso::mink(4,edge(3,1)))" particle="d"];
          0:6	-> 3:7	 [id=4 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="-2"  lmb_id="1" lmb_rep="K(1,a___)" name="e4" num="1" particle="G+"];
          1:8	-> 2:9	 [id=5 dir=back source="{ufo_order:1}" sink="{ufo_order:0}"  dod="-1"  lmb_rep="K(0,a___)+K(1,a___)+P(0,a___)" name="e5" num="-1*Q(5,spenso::mink(4,edge(5,1)))*spenso::g(spenso::cof(3,hedge(9)),spenso::dind(spenso::cof(3,hedge(8))))*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(9)),spenso::mink(4,edge(5,1)))" particle="d~"];
        }
        "#);
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
          overall_factor_evaluated = "1";
          projector = "ϵ(0,spenso::mink(4,hedge(0)))*ϵbar(0,spenso::mink(4,hedge(1)))";
          0[dod="0" int_id="V_82" num="UFO::GC_16*spenso::g(spenso::cof(3,hedge(4)),spenso::dind(spenso::cof(3,hedge(2))))*spenso::projp(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(4)))"];
          1[dod="0" int_id="V_117" num="UFO::GC_22*spenso::g(spenso::cof(3,hedge(3)),spenso::dind(spenso::cof(3,hedge(8))))*spenso::projm(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(3)))"];
          2[dod="0" int_id="V_71" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(9)),spenso::dind(spenso::cof(3,hedge(5))))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(9)),spenso::mink(4,hedge(1)))"];
          3[dod="1" int_id="V_11" num="(-1*Q(5,spenso::mink(4,hedge(0)))+Q(4,spenso::mink(4,hedge(0))))*UFO::GC_3"];

          2:1	-> 3:0	 [id=0 dir=none source="{ufo_order:2}" sink="{ufo_order:0}"  dod="-2" is_cut="0"  lmb_rep="P(0,a___)" name="e0" num="1" particle="a"];
          1:8	-> 2:9	 [id=1 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_rep="K(0,a___)+K(1,a___)+P(0,a___)" name="e1" num="Q(1,spenso::mink(4,edge(1,1)))*spenso::g(spenso::cof(3,hedge(8)),spenso::dind(spenso::cof(3,hedge(9))))*spenso::gamma(spenso::bis(4,hedge(9)),spenso::bis(4,hedge(8)),spenso::mink(4,edge(1,1)))" particle="d"];
          0:2	-> 1:3	 [id=2 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_id="0" lmb_rep="K(0,a___)" name="e2" num="(Q(2,spenso::mink(4,edge(2,1)))*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(2)),spenso::mink(4,edge(2,1)))+UFO::MC*spenso::g(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(3))))*spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(3))))" particle="c"];
          0:4	-> 2:5	 [id=3 dir=back source="{ufo_order:1}" sink="{ufo_order:0}"  dod="-1"  lmb_rep="-1*K(0,a___)+-1*K(1,a___)" name="e3" num="-1*Q(3,spenso::mink(4,edge(3,1)))*spenso::g(spenso::cof(3,hedge(5)),spenso::dind(spenso::cof(3,hedge(4))))*spenso::gamma(spenso::bis(4,hedge(4)),spenso::bis(4,hedge(5)),spenso::mink(4,edge(3,1)))" particle="d~"];
          0:6	-> 3:7	 [id=4 dir=none source="{ufo_order:2}" sink="{ufo_order:1}"  dod="-2"  lmb_id="1" lmb_rep="K(1,a___)" name="e4" num="1" particle="G-"];
          1:10	-> 3:11	 [id=5 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="-2"  lmb_rep="-1*K(1,a___)+-1*P(0,a___)" name="e5" num="1" particle="G+"];
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
            2:0-> ext0 [id=0 dir=back is_cut=0  particle="a"];
            ext1 [style=invis];
            ext1-> 3:1 [id=1 is_cut=0  particle="a"];
            0:2-> 1:3 [id=2   particle="d"];
            0:4-> 1:5 [id=3   particle="Z"];
            0:6-> 3:7 [id=4   particle="d~"];
            1:8-> 2:9 [id=5   particle="d"];
            2:10-> 3:11 [id=6   particle="d"];
            }
        ).unwrap();

        assert_snapshot!(g.debug_dot(),@r#"
        digraph GL0{
          num = "1";
          overall_factor = "(AutG(1))^(-1)*AntiFermionSpinSumSign(1)*ExternalFermionOrderingSign(1)*InternalFermionLoopSign(-1)";
          overall_factor_evaluated = "-1";
          projector = "ϵ(0,spenso::mink(4,hedge(1)))*ϵbar(0,spenso::mink(4,hedge(0)))";
          0[dod="0" int_id="V_79" num="((-2*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::projp(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(6)))+spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::projm(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(6))))*UFO::GC_58+UFO::GC_50*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,vertex(0,1)),spenso::mink(4,hedge(4)))*spenso::projm(spenso::bis(4,vertex(0,1)),spenso::bis(4,hedge(6))))*spenso::g(spenso::cof(3,hedge(6)),spenso::dind(spenso::cof(3,hedge(2))))"];
          1[dod="0" int_id="V_79" num="((-2*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::projp(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3)))+spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::projm(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3))))*UFO::GC_58+UFO::GC_50*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,vertex(1,1)),spenso::mink(4,hedge(5)))*spenso::projm(spenso::bis(4,vertex(1,1)),spenso::bis(4,hedge(3))))*spenso::g(spenso::cof(3,hedge(3)),spenso::dind(spenso::cof(3,hedge(8))))"];
          2[dod="0" int_id="V_71" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(9)),spenso::dind(spenso::cof(3,hedge(10))))*spenso::gamma(spenso::bis(4,hedge(10)),spenso::bis(4,hedge(9)),spenso::mink(4,hedge(0)))"];
          3[dod="0" int_id="V_71" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(11)),spenso::dind(spenso::cof(3,hedge(7))))*spenso::gamma(spenso::bis(4,hedge(7)),spenso::bis(4,hedge(11)),spenso::mink(4,hedge(1)))"];

          2:1	-> 3:0	 [id=0 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="-2" is_cut="0"  lmb_rep="P(0,a___)" name="e0" num="1" particle="a"];
          2:10	-> 3:11	 [id=1 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_rep="-1*K(1,a___)+-1*P(0,a___)" name="e1" num="Q(1,spenso::mink(4,edge(1,1)))*spenso::g(spenso::cof(3,hedge(10)),spenso::dind(spenso::cof(3,hedge(11))))*spenso::gamma(spenso::bis(4,hedge(11)),spenso::bis(4,hedge(10)),spenso::mink(4,edge(1,1)))" particle="d"];
          0:2	-> 1:3	 [id=2 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_id="0" lmb_rep="K(0,a___)" name="e2" num="Q(2,spenso::mink(4,edge(2,1)))*spenso::g(spenso::cof(3,hedge(2)),spenso::dind(spenso::cof(3,hedge(3))))*spenso::gamma(spenso::bis(4,hedge(3)),spenso::bis(4,hedge(2)),spenso::mink(4,edge(2,1)))" particle="d"];
          0:4	-> 1:5	 [id=3 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="0"  lmb_rep="-1*K(0,a___)+-1*K(1,a___)" name="e3" num="-1*spenso::g(spenso::mink(4,hedge(4)),spenso::mink(4,hedge(5)))+Q(3,spenso::mink(4,hedge(4)))*Q(3,spenso::mink(4,hedge(5)))*UFO::MZ^(-2)" particle="Z"];
          0:6	-> 3:7	 [id=4 dir=back source="{ufo_order:1}" sink="{ufo_order:0}"  dod="-1"  lmb_id="1" lmb_rep="K(1,a___)" name="e4" num="-1*Q(4,spenso::mink(4,edge(4,1)))*spenso::g(spenso::cof(3,hedge(7)),spenso::dind(spenso::cof(3,hedge(6))))*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,hedge(7)),spenso::mink(4,edge(4,1)))" particle="d~"];
          1:8	-> 2:9	 [id=5 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_rep="-1*K(1,a___)" name="e5" num="Q(5,spenso::mink(4,edge(5,1)))*spenso::g(spenso::cof(3,hedge(8)),spenso::dind(spenso::cof(3,hedge(9))))*spenso::gamma(spenso::bis(4,hedge(9)),spenso::bis(4,hedge(8)),spenso::mink(4,edge(5,1)))" particle="d"];
        }
        "#);
    }

    #[test]
    fn test_pols_tree() {
        test_initialise().unwrap();
        let g: Graph = dot!(digraph qqx_aaa_tree_1 {
                    ext    [style=invis]
                    ext -> v1:1 [particle="d" id=1];
                    ext -> v3:2 [dir=back particle="d~" id=2];
                    v1:3 -> ext [particle="a" id=3];
                    v2:4 -> ext [particle="a" id=4];
                    v3:0 -> ext [particle="a" id=0];
                    v1:5 -> v2:6 [particle="d" id=5];
                    v2 -> v3 [particle="d" id=6];
        })
        .unwrap();

        assert_snapshot!(g.debug_dot(),@r#"
        digraph qqx_aaa_tree_1{
          num = "1";
          overall_factor = "1";
          overall_factor_evaluated = "1";
          projector = "u(1,spenso::bis(4,hedge(1)))*vbar(2,spenso::bis(4,hedge(2)))*ϵbar(0,spenso::mink(4,hedge(0)))*ϵbar(3,spenso::mink(4,hedge(3)))*ϵbar(4,spenso::mink(4,hedge(4)))";
          v1[dod="0" int_id="V_71" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(1)),spenso::dind(spenso::cof(3,hedge(5))))*spenso::gamma(spenso::bis(4,hedge(5)),spenso::bis(4,hedge(1)),spenso::mink(4,hedge(3)))"];
          v2[dod="0" int_id="V_71" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(6)),spenso::dind(spenso::cof(3,hedge(7))))*spenso::gamma(spenso::bis(4,hedge(7)),spenso::bis(4,hedge(6)),spenso::mink(4,hedge(4)))"];
          v3[dod="0" int_id="V_71" num="UFO::GC_1*spenso::g(spenso::cof(3,hedge(8)),spenso::dind(spenso::cof(3,hedge(2))))*spenso::gamma(spenso::bis(4,hedge(2)),spenso::bis(4,hedge(8)),spenso::mink(4,hedge(0)))"];

          ext0	 [style=invis];
          v3:0	-> ext0	 [id=0 dir=none source="{ufo_order:2}" dod="-2"  lmb_rep="P(0,a___)" name="e0" num="1𝑖" particle="a"];
          ext1	 [style=invis];
          ext1	-> v1:1	 [id=1 sink="{ufo_order:1}" dod="-2"  lmb_rep="P(1,a___)" name="e1" num="1𝑖" particle="d"];
          ext2	 [style=invis];
          ext2	-> v3:2	 [id=2 dir=back sink="{ufo_order:0}" dod="-2"  lmb_rep="P(2,a___)" name="e2" num="1𝑖" particle="d~"];
          ext3	 [style=invis];
          v1:3	-> ext3	 [id=3 dir=none source="{ufo_order:2}" dod="-2"  lmb_rep="P(3,a___)" name="e3" num="1𝑖" particle="a"];
          ext4	 [style=invis];
          v2:4	-> ext4	 [id=4 dir=none source="{ufo_order:2}" dod="-2"  lmb_rep="-1*P(0,a___)+-1*P(3,a___)+P(1,a___)+P(2,a___)" name="e4" num="1𝑖" particle="a"];
          v1:5	-> v2:6	 [id=5 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_rep="-1*P(3,a___)+P(1,a___)" name="e5" num="Q(5,spenso::mink(4,edge(5,1)))*spenso::g(spenso::cof(3,hedge(5)),spenso::dind(spenso::cof(3,hedge(6))))*spenso::gamma(spenso::bis(4,hedge(6)),spenso::bis(4,hedge(5)),spenso::mink(4,edge(5,1)))" particle="d"];
          v2:7	-> v3:8	 [id=6 source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-1"  lmb_rep="-1*P(2,a___)+P(0,a___)" name="e6" num="Q(6,spenso::mink(4,edge(6,1)))*spenso::g(spenso::cof(3,hedge(7)),spenso::dind(spenso::cof(3,hedge(8))))*spenso::gamma(spenso::bis(4,hedge(8)),spenso::bis(4,hedge(7)),spenso::mink(4,edge(6,1)))" particle="d"];
        }
        "#);
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
                .dot_lmb_of(&g.underlying.full(), &g.loop_momentum_basis)
        );

        let num = Numerator::<UnInit>::default().from_new_graph(g, &g.underlying.full_filter());

        let expr = num.state.expr.as_view();

        let lib: DummyLibrary<SymbolicTensor<Aind>> = DummyLibrary::<_>::new();
        let net =
            Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, _, Symbol, Aind>::try_from_view::<
                SymbolicTensor<Aind>,
                _,
            >(expr, &lib, &ParseSettings::default())
            .unwrap();

        println!("{}", expr);
        println!(
            "{}",
            net.dot_display_impl(
                ToString::to_string,
                |_| None,
                |a| {
                    if let Ok(a) = PermutedStructure::<ShadowedStructure<Aind>>::try_from(
                        a.expression.as_view(),
                    ) {
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
                },
                ToString::to_string
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
          overall_factor_evaluated = "2*x";
          projector = "1";
          A[dod="0" num="1"];
          B[dod="0" num="1"];
          C[dod="0" num="1"];
          D[dod="0" num="1"];
          E[dod="0" num="1"];

          A:0	-> D:1	 [id=0 dir=none source="{ufo_order:0}" sink="{ufo_order:0}"  dod="-2"  lmb_rep="K(1,a___)+K(2,a___)+P(2,a___)" name="e0" num="P(0,spenso::mink(4,0))" particle="scalar_0"];
          ext1	 [style=invis];
          B:2	-> ext1	 [id=1 dir=none source="{ufo_order:0}" dod="-2"  lmb_rep="P(0,a___)" name="e1" num="1𝑖" particle="scalar_0"];
          ext2	 [style=invis];
          B:3	-> ext2	 [id=2 dir=none source="{ufo_order:1}" dod="-2"  lmb_rep="P(1,a___)" name="e2" num="1𝑖" particle="scalar_0"];
          C:4	-> A:5	 [id=3 dir=none source="{ufo_order:0}" sink="{ufo_order:1}"  dod="-2"  lmb_id="2" lmb_rep="K(2,a___)" name="e3" num="P(3,spenso::mink(4,0))" particle="scalar_0"];
          C:6	-> D:7	 [id=4 dir=none source="{ufo_order:1}" sink="{ufo_order:1}"  dod="-2"  lmb_rep="-1*K(0,a___)+-1*K(1,a___)+-1*K(2,a___)+-1*P(2,a___)+P(0,a___)+P(1,a___)" name="e4" num="1" particle="scalar_0"];
          D:8	-> B:9	 [id=5 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="-2"  lmb_rep="-1*K(0,a___)+P(0,a___)+P(1,a___)" name="e5" num="1" particle="scalar_0"];
          E:10	-> A:11	 [id=6 dir=none source="{ufo_order:0}" sink="{ufo_order:2}"  dod="-2"  lmb_id="1" lmb_rep="K(1,a___)" name="e6" num="1" particle="scalar_0"];
          E:12	-> B:13	 [id=7 dir=none source="{ufo_order:1}" sink="{ufo_order:3}"  dod="-2"  lmb_id="0" lmb_rep="K(0,a___)" name="e7" num="1" particle="scalar_0"];
          E:14	-> C:15	 [id=8 dir=none source="{ufo_order:2}" sink="{ufo_order:2}"  dod="-2"  lmb_rep="-1*K(0,a___)+-1*K(1,a___)" name="e8" num="1" particle="scalar_0"];
          ext9	 [style=invis];
          ext9	-> A:16	 [id=9 dir=none sink="{ufo_order:3}" dod="-2"  lmb_rep="P(2,a___)" name="e9" num="1𝑖" particle="scalar_0"];
          ext10	 [style=invis];
          ext10	-> C:17	 [id=10 dir=none sink="{ufo_order:3}" dod="-2"  lmb_rep="-1*P(2,a___)+P(0,a___)+P(1,a___)" name="e10" num="1𝑖" particle="scalar_0"];
        }
        "#);

        println!(
            "{}",
            g.underlying
                .dot_lmb_of(&g.underlying.full(), &g.loop_momentum_basis)
        );

        let num = Numerator::<UnInit>::default().from_new_graph(&g, &g.underlying.full_filter());

        let expr = num.state.expr.as_view();

        let lib: DummyLibrary<SymbolicTensor<Aind>> = DummyLibrary::<_>::new();
        let net =
            Network::<NetworkStore<SymbolicTensor<Aind>, Atom>, _, Symbol, Aind>::try_from_view::<
                SymbolicTensor<Aind>,
                _,
            >(expr, &lib, &ParseSettings::default())
            .unwrap();

        println!("{}", expr);
        println!(
            "{}",
            net.dot_display_impl(
                ToString::to_string,
                |_| None,
                |a| {
                    if let Ok(a) = PermutedStructure::<ShadowedStructure<Aind>>::try_from(
                        a.expression.as_view(),
                    ) {
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
                },
                ToString::to_string
            )
        );

        let num = num.color_simplify();

        println!("{}", num.state.expr);
        let num = num.gamma_simplify();

        println!("{}", num.state.expr);
    }

    #[test]
    fn explicit_hedge_payload_round_trips_in_dot_export() {
        test_initialise().unwrap();
        let g: Graph = dot!(
            digraph payload_graph{
                ext_in [style=invis]
                ext_out [style=invis]
                A [num=1 dod=0]
                ext_in -> A [name=e_in num=1 dod=-2 sink="{ufo_order:0,dod:-2}"]
                A -> ext_out [name=e_out num=1 dod=-2 source="{ufo_order:1,dod:-2}"]
            }
        )
        .unwrap();

        let serialized = g.dot_serialize(&DotExportSettings::default());
        assert!(serialized.contains("ufo_order"));
        assert!(serialized.contains("dod"));
        assert!(serialized.contains("source=") || serialized.contains("sink="));
    }

    #[test]
    #[should_panic(expected = "vertex dod autogeneration is not implemented yet")]
    fn missing_vertex_dod_panics() {
        test_initialise().unwrap();
        let _: Graph = dot!(
            digraph missing_vertex_dod{
                ext_in [style=invis]
                ext_out [style=invis]
                A [num=1]
                ext_in -> A [name=e_in num=1 dod=-2 sink="{ufo_order:0,dod:-2}"]
                A -> ext_out [name=e_out num=1 dod=-2 source="{ufo_order:1,dod:-2}"]
            }
        )
        .unwrap();
    }

    #[test]
    #[should_panic(expected = "hedge dod autogeneration is not implemented yet")]
    fn missing_hedge_dod_panics() {
        test_initialise().unwrap();
        let _: Graph = dot!(
            digraph missing_hedge_dod{
                ext_in [style=invis]
                ext_out [style=invis]
                A [num=1 dod=0]
                ext_in -> A [name=e_in num=1 dod=-2 sink="{ufo_order:0}"]
                A -> ext_out [name=e_out num=1 dod=-2 source="{ufo_order:1}"]
            }
        )
        .unwrap();
    }
}
