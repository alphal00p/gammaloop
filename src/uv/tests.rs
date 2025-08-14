use crate::cff::expression::OrientationData;
use crate::momentum::SignOrZero;
use crate::momentum_sample::LoopIndex;
use crate::new_cs::Amplitude;
use crate::new_graph::edge::ParseEdge;
use crate::new_graph::global::ParseData;
use crate::new_graph::hedge_data::ParseHedge;
use crate::new_graph::parse::ParseGraph;
use crate::new_graph::vertex::ParseVertex;
use crate::new_graph::{LMBext, LoopMomentumBasis};
use crate::utils::W_;
use crate::uv::uv_graph::UVE;
use idenso::metric::MS;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, Orientation};

use linnet::half_edge::{builder::HedgeGraphBuilder, involution::Flow};
use symbolica::symbol;
use typed_index_collections::TiVec;

use crate::{
    cff::{cut_expression::SuperGraphOrientationID, expression::AmplitudeOrientationID},
    dot,
    feyngen::FeynGenFilters,
    inspect::inspect,
    integrands::Integrand,
    integrate::havana_integrate,
    momentum_sample::ExternalIndex,
    new_cs::{AmplitudeGraph, CrossSection, CrossSectionGraph, CutId, ProcessDefinition},
    new_graph::{parse::IntoGraph, ExternalConnection, FeynmanGraph, Graph},
    signature::LoopExtSignature,
    utils::test_utils::load_generic_model,
    utils::{external_energy_atom_from_index, F},
    GenerationSettings, RuntimeSettings,
};
use spenso::{algebra::algebraic_traits::IsZero, network::library::TensorLibraryData};
use symbolica::evaluate::{FunctionMap, OptimizationSettings};
use symbolica::{atom::AtomCore, function, parse};

#[test]
fn tri_uv_AMP() {
    let _ = env_logger::builder().is_test(true).try_init();
    let is_massless = true;
    let model = if is_massless {
        load_generic_model("sm_massless")
    } else {
        load_generic_model("sm")
    };

    let graph: Graph = dot!(
        digraph G{
            e        [style=invis]
            e -> A   [particle=H id=0]
            e -> B   [particle=H id=1]
            e -> C   [particle=H id=2]
            A -> B    [particle=H lmb_index=0]
            B -> C    [particle=H]
            C -> D    [particle=H]
        },&model
    )
    .unwrap();

    let mut amp = Amplitude::new("");
    amp.add_graph(graph);

    amp.preprocess(&model, &GenerationSettings::default()).unwrap();
    amp.build_integrand(RuntimeSettings::default(), &model).unwrap();
}

#[test]
fn tri_box_tri_LU() {
    let _ = env_logger::builder().is_test(true).try_init();
    let uv_dod = 1;
    let box_uv_dod = -1; // can be -1, 0, 1, 2
    let is_massless = false;
    let force_cut: Option<CutId> = None;
    //let force_cut: Option<CutId> = Some(CutId(3));

    // load the model and hack the masses, go through serializable model since arc is not mutable
    let model = if is_massless {
        load_generic_model("sm_massless")
    } else {
        load_generic_model("sm")
    };

    let mut underlying = HedgeGraphBuilder::new();

    let hhh = model.get_vertex_rule("V_9");
    let htt = model.get_vertex_rule("V_141");

    let hprop = model.get_propagator("H_propFeynman");
    let hp = model.get_particle("H");

    if is_massless {
        assert!(hp.0.mass.value.unwrap().is_zero());
    }

    let tprop = model.get_propagator("t_propFeynman");
    let tp = model.get_particle("t");

    if is_massless {
        assert!(tp.0.mass.value.unwrap().is_zero());
    }

    let n1 = underlying.add_node(ParseVertex::from(hhh.clone()).with_num(Atom::num(1)));
    let n2 = underlying.add_node(ParseVertex::from(hhh.clone()).with_num(Atom::num(1)));
    let n3 = underlying.add_node(ParseVertex::from(hhh.clone()).with_num(Atom::num(1)));
    let n4 = underlying.add_node(ParseVertex::from(htt.clone()).with_num(Atom::num(1)));
    let n5 = underlying.add_node(ParseVertex::from(htt.clone()).with_num(Atom::num(1)));
    let n6 = underlying.add_node(ParseVertex::from(htt.clone()).with_num(Atom::num(1)));

    underlying.add_edge(
        n1.add_data(ParseHedge::default()),
        n2.add_data(ParseHedge::default()),
        ParseEdge::new(hp.clone())
            .with_num(Atom::one())
            .with_lmb_id(LoopIndex(0)),
        false,
    );

    underlying.add_edge(
        n1.add_data(ParseHedge::default()),
        n3.add_data(ParseHedge::default()),
        ParseEdge::new(hp.clone()).with_num(Atom::one()),
        false,
    );

    underlying.add_edge(
        n2.add_data(ParseHedge::default()),
        n3.add_data(ParseHedge::default()),
        ParseEdge::new(hp.clone()).with_num(if box_uv_dod >= 1 {
            spenso_lor_atom(2, 30, GS.dim)
        } else {
            Atom::one()
        }),
        // Edge {
        //     name: "e2".into(),

        //     particle: hp.clone(),
        //     propagator: hprop.clone(),

        //     dod: if box_uv_dod >= 1 { -1 } else { -2 },
        //     num: if box_uv_dod >= 1 {
        //         spenso_lor_atom(2, 30, GS.dim)
        //     } else {
        //         Atom::one()
        //     },
        // },
        false,
    );

    underlying.add_edge(
        n2.add_data(ParseHedge::default()),
        n4.add_data(ParseHedge::default()),
        ParseEdge::new(hp.clone()).with_num(if box_uv_dod >= 0 && uv_dod >= 1 {
            spenso_lor_atom(3, 20, GS.dim)
        } else {
            Atom::one()
        }),
        false,
    );

    underlying.add_edge(
        n3.add_data(ParseHedge::default()),
        n5.add_data(ParseHedge::default()),
        ParseEdge::new(hp.clone())
            .with_num(if box_uv_dod == 2 {
                spenso_lor_atom(4, 30, GS.dim)
            } else {
                Atom::one()
            })
            .with_lmb_id(LoopIndex(1)),
        false,
    );

    underlying.add_edge(
        n4.add_data(ParseHedge::default()),
        n5.add_data(ParseHedge::default()),
        ParseEdge::new(tp.clone()).with_num(if uv_dod >= 0 {
            spenso_lor_atom(5, 10, GS.dim)
        } else {
            Atom::one()
        }),
        true,
    );

    underlying.add_edge(
        n6.add_data(ParseHedge::default()),
        n4.add_data(ParseHedge::default()),
        ParseEdge::new(tp.clone()).with_num(if uv_dod >= 1 {
            spenso_lor_atom(6, 20, GS.dim)
        } else {
            Atom::one()
        }),
        true,
    );

    underlying.add_edge(
        n5.add_data(ParseHedge::default()),
        n6.add_data(ParseHedge::default()),
        ParseEdge::new(tp.clone())
            .with_num(if uv_dod >= 0 {
                spenso_lor_atom(7, 10, GS.dim)
            } else {
                Atom::one()
            })
            .with_lmb_id(LoopIndex(2)),
        true,
    );

    underlying.add_external_edge(
        n1.add_data(ParseHedge::default()),
        ParseEdge::new(hp.clone()).with_num(if box_uv_dod == 1 {
            spenso_lor_atom(8, 30, GS.dim)
        } else {
            Atom::one()
        }),
        false,
        Flow::Sink,
    );

    underlying.add_external_edge(
        n6.add_data(ParseHedge::default()),
        ParseEdge::new(hp.clone()).with_num(if box_uv_dod == -1 && uv_dod >= 1 {
            spenso_lor_atom(9, 20, GS.dim)
        } else {
            Atom::one()
        }),
        false,
        Flow::Source,
    );

    let underlying = ParseGraph {
        graph: underlying.build(),
        global_data: ParseData::default(),
    };

    println!("{}", underlying.dot(&underlying.full_filter()));

    let graph = Graph::from_parsed(underlying, &model).unwrap();

    println!(
        "dot lmb:{}",
        graph
            .underlying
            .dot_lmb(&graph.underlying.full_filter(), &graph.loop_momentum_basis)
    );

    let hpdg = hp.pdg_code as i64;
    let tpdg = tp.pdg_code as i64;
    let mut cs: CrossSection = CrossSection::new("".into());
    cs.preprocess(
        &model,
        &ProcessDefinition {
            initial_pdgs: vec![hpdg],
            final_pdgs_lists: vec![
                vec![tpdg, -tpdg],
                vec![tpdg, -tpdg, hpdg],
                vec![hpdg, hpdg],
                vec![hpdg, hpdg, hpdg],
                vec![tpdg, -tpdg, hpdg, hpdg],
            ],
            n_unresolved: 0,
            unresolved_cut_content: HashSet::new(),
            amplitude_filters: FeynGenFilters(vec![]),
            cross_section_filters: FeynGenFilters(vec![]),
        },
    )
    .unwrap();
    cs.build_integrand(RuntimeSettings::default(), &model).unwrap();

    //println!("Final result: {:>}", sum.expand());
}

#[test]
fn double_triangle_LU() {
    let _ = env_logger::builder().is_test(true).try_init();
    let uv_dod = 1;

    // load the model and hack the masses, go through serializable model since arc is not mutable
    let model = load_generic_model("sm");

    let mut underlying = HedgeGraphBuilder::new();

    let hhh = model.get_vertex_rule("V_9");
    let htt = model.get_vertex_rule("V_141");

    let hprop = model.get_propagator("H_propFeynman");
    let hp = model.get_particle("H");

    let tprop = model.get_propagator("t_propFeynman");
    let tp = model.get_particle("t");

    let n1 = underlying.add_node(hhh.clone().into());
    let n2 = underlying.add_node(htt.clone().into());
    let n3 = underlying.add_node(htt.clone().into());
    let n4 = underlying.add_node(htt.clone().into());

    underlying.add_edge(
        n1.add_data(ParseHedge::default()),
        n2.add_data(ParseHedge::default()),
        ParseEdge::new(hp.clone()),
        false,
    );

    underlying.add_edge(
        n1.add_data(ParseHedge::default()),
        n3.add_data(ParseHedge::default()),
        ParseEdge::new(hp.clone()),
        false,
    );

    underlying.add_edge(
        n2.add_data(ParseHedge::default()),
        n3.add_data(ParseHedge::default()),
        ParseEdge::new(tp.clone()).with_num(if uv_dod >= 0 {
            spenso_lor_atom(2, 10, GS.dim)
        } else {
            Atom::one()
        }),
        true,
    );

    underlying.add_edge(
        n3.add_data(ParseHedge::default()),
        n4.add_data(ParseHedge::default()),
        ParseEdge::new(tp.clone()).with_num(if uv_dod >= 1 {
            spenso_lor_atom(3, 20, GS.dim)
        } else {
            Atom::one()
        }),
        true,
    );

    underlying.add_edge(
        n4.add_data(ParseHedge::default()),
        n2.add_data(ParseHedge::default()),
        ParseEdge::new(tp.clone()).with_num(if uv_dod >= 0 {
            spenso_lor_atom(4, 10, GS.dim)
        } else {
            Atom::one()
        }),
        true,
    );

    underlying.add_external_edge(
        n1.add_data(ParseHedge::default()),
        ParseEdge::new(hp.clone()),
        false,
        Flow::Sink,
    );

    underlying.add_external_edge(
        n4.add_data(ParseHedge::default()),
        ParseEdge::new(hp.clone()).with_num(if uv_dod >= 1 {
            spenso_lor_atom(6, 20, GS.dim)
        } else {
            Atom::one()
        }),
        false,
        Flow::Source,
    );

    let underlying = ParseGraph {
        graph: underlying.build(),
        global_data: ParseData::default(),
    };

    let mut loop_momentum_basis = LoopMomentumBasis {
        tree: None,
        loop_edges: vec![EdgeIndex::from(0), EdgeIndex::from(3)].into(),
        ext_edges: vec![EdgeIndex::from(5), EdgeIndex::from(6)].into(),
        edge_signatures: underlying.new_edgevec(|_, _, _| LoopExtSignature::from((vec![], vec![]))),
    };

    let graph = Graph::from_parsed(underlying, &model).unwrap();

    let mut cs: CrossSection = CrossSection::new("".into());
    cs.add_supergraph(graph);

    let hpdg = hp.pdg_code as i64;
    let tpdg = tp.pdg_code as i64;
    cs.preprocess(
        &model,
        &ProcessDefinition {
            initial_pdgs: vec![hpdg],
            final_pdgs_lists: vec![vec![tpdg, -tpdg], vec![tpdg, -tpdg, hpdg], vec![hpdg, hpdg]],
            n_unresolved: 0,
            unresolved_cut_content: HashSet::new(),
            amplitude_filters: FeynGenFilters(vec![]),
            cross_section_filters: FeynGenFilters(vec![]),
        },
    )
    .unwrap();

    cs.build_integrand(RuntimeSettings::default(), &model).unwrap();

    //println!("Final result: {:>}", sum.expand());
}

/*

#[test]
fn nested_bubble_soft_ct() {
    // let massive_particle = ArcParticle(Arc::new(Particle::))

    let scalar_node = UVNode {
        dod: 0,
        num: Atom::num(1),
        color: None,
    };

    fn scalar_edge(eid: i32) -> UVEdge {
        let model = load_generic_model("sm");
        let higgs = model.get_particle("H");
        let m2 = parse!(higgs.mass.name).npow(2);
        UVEdge {
            og_edge: 1, // not needed
            dod: -2,
            particle: higgs,
            num: Atom::num(1),
            den: spenso_lor_atom(eid, 1, GS.dim).npow(2).to_dots() - m2,
        }
    }

    fn scalar_edge_with_p(eid: i32, ind: impl Into<AbstractIndex>) -> UVEdge {
        let model = load_generic_model("sm");
        let higgs = model.get_particle("H");
        let m2 = parse!(higgs.mass.name).npow(2);
        UVEdge {
            og_edge: 1, // not needed
            dod: -1,
            particle: higgs,
            num: spenso_lor_atom(eid, ind, GS.dim),
            den: spenso_lor_atom(eid, 1, GS.dim).npow(2).to_dots() - m2,
        }
    }

    let mut builder = HedgeGraphBuilder::new();

    let node1 = builder.add_node(scalar_node.clone());
    let node2 = builder.add_node(scalar_node.clone());
    let node3 = builder.add_node(scalar_node.clone());
    let node4 = builder.add_node(scalar_node.clone());

    builder.add_edge(node1, node2, scalar_edge(0), true);
    builder.add_edge(node1, node4, scalar_edge(1), true);
    builder.add_edge(node2, node3, scalar_edge_with_p(2, 1), true);
    builder.add_edge(node2, node3, scalar_edge_with_p(3, 1), true);
    builder.add_edge(node3, node4, scalar_edge(4), true);
    builder.add_external_edge(node1, scalar_edge(5), true, Flow::Sink);
    builder.add_external_edge(node4, scalar_edge(6), true, Flow::Source);

    let hedge_graph = builder.build();
    let uv_graph = UVGraph::from_hedge(hedge_graph);
    println!(
        "{}",
        uv_graph.dot_impl(
            &uv_graph.full_filter(),
            "",
            &|a| Some(format!("label=\"{}/{}\"", a.num, a.den)),
            &|n| Some(format!("label=\"{}\"", n.num))
        )
    );

    let wood = uv_graph.wood(&uv_graph.full_graph());

    //println!("{}", wood.dot(&uv_graph));
    //println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.lmb);
    // assert_eq!(152, ufold.n_terms());
    // ufold.compute(&uv_graph);

    // println!("unfolded : {}", ufold.show_structure(&uv_graph).unwrap());
    // println!("graph: {}", ufold.graphs());

    let result = ufold.expr(&uv_graph, &uv_graph.full_graph()).unwrap().0;

    let result = result
        .replace(function!(GS.emr_mom, W_.x_, W_.y_))
        .with(W_.y_);

    // println!("{}", ufold.structure_and_res(&uv_graph));
    println!("{:>}", result);

    let t = symbol!("t");
    let series = Atom::var(t).npow(4)
        * result
            .replace(parse!("K(3)"))
            .with(parse!("t*K(3)"))
            .replace(parse!("spenso::dot(t*K(x_),y_)"))
            .repeat()
            .with(parse!("t*spenso::dot(K(x_),y_)"));

    println!("SERIESINPUT:\n{:>}", series);

    let s = series
        .replace(t)
        .with(Atom::var(t).npow(-1))
        .series(t, Atom::Zero, 0.into(), true)
        .unwrap();
    println!("Series: {:>}", s);
    println!("Correct UV cancellation if 0: {:>}", s.to_atom().expand());

    let exp = result
        .replace(parse!("spenso::dot(k_(x_),l_(y_))"))
        .with(parse!(
            "k(x_,0)*k(y_,0)-k(x_,1)*k(y_,1)-k(x_,2)*k(y_,2)-k(x_,3)*k(y_,3)"
        ));

    let mut fnmap = FunctionMap::new();

    fnmap.add_constant(
        Atom::var(symbol!("m")),
        symbolica::domains::float::Complex::new((1.).into(), (0.).into()),
    );
    fnmap.add_constant(
        Atom::var(symbol!("MH")),
        symbolica::domains::float::Complex::new((1.).into(), (0.).into()),
    );
    fnmap.add_constant(
        Atom::var(symbol!("mUV")),
        symbolica::domains::float::Complex::new((10.).into(), (0.).into()),
    );
    fnmap.add_constant(
        Atom::var(symbol!("ZERO")),
        symbolica::domains::float::Complex::new((0.).into(), (0.).into()),
    );
    let ev = exp
        .evaluator(
            &fnmap,
            &[
                parse!("k(1, 0)"),
                parse!("k(1, 1)"),
                parse!("k(1, 2)"),
                parse!("k(1, 3)"),
                parse!("k(3, 0)"),
                parse!("k(3, 1)"),
                parse!("k(3, 2)"),
                parse!("k(3, 3)"),
                parse!("k(5, 0)"),
                parse!("k(5, 1)"),
                parse!("k(5, 2)"),
                parse!("k(5, 3)"),
                parse!("k(6, 0)"),
                parse!("k(6, 1)"),
                parse!("k(6, 2)"),
                parse!("k(6, 3)"),
            ],
            OptimizationSettings::default(),
        )
        .unwrap();

    let mut ev2 = ev.map_coeff(&|x| x.re.to_f64());

    println!("Single limit:");
    for t in (0..100_000).step_by(5000) {
        let r = (t as f64).powf(4.)
            * ev2.evaluate_single(&[
                3.,
                43.,
                5.,
                6.5,
                t as f64 + 1.,
                t as f64 * 2. + 2.,
                t as f64 + 3.,
                t as f64 + 4.,
                7.,
                4.2,
                1.,
                2.4,
                5.,
                2.3,
                9.8,
                1.,
                2.,
            ]);
        println!("{} {}", r, r.abs().log10());
    }

    println!("Full graph limit:");
    for t in (0..100_000).step_by(5000) {
        let r = (t as f64).powf(8.)
            * ev2.evaluate_single(&[
                t as f64 + 3.,
                t as f64 + 43.,
                t as f64 * 3. + 5.,
                t as f64 + 6.5,
                t as f64 + 1.,
                t as f64 * 2. + 2.,
                t as f64 + 3.,
                t as f64 + 4.,
                7.,
                4.2,
                1.,
                2.4,
                2.3,
                9.8,
                1.,
                2.,
            ]);
        println!("{} {}", r, r.abs().log10());
    }
}

#[test]
fn nested_bubble_scalar_quad() {
    // let massive_particle = ArcParticle(Arc::new(Particle::))

    let scalar_node = UVNode {
        dod: 0,
        num: Atom::num(1),
        color: None,
    };

    fn scalar_edge(eid: i32) -> UVEdge {
        let model = load_generic_model("sm");
        let higgs = model.get_particle("H");
        let m2 = parse!(higgs.mass.name).npow(2);
        UVEdge {
            og_edge: 1, // not needed
            dod: -2,
            particle: higgs,
            num: Atom::num(1),
            den: spenso_lor_atom(eid, 1, GS.dim).npow(2).to_dots() - m2,
        }
    }

    fn scalar_edge_with_p(eid: i32, ind: impl Into<AbstractIndex>) -> UVEdge {
        let model = load_generic_model("sm");
        let higgs = model.get_particle("H");
        let m2 = parse!(higgs.mass.name).npow(2);
        UVEdge {
            og_edge: 1, // not needed
            dod: -1,
            particle: higgs,
            num: spenso_lor_atom(eid, ind, GS.dim),
            den: spenso_lor_atom(eid, 1, GS.dim).npow(2).to_dots() - m2,
        }
    }

    let mut builder = HedgeGraphBuilder::new();

    let node1 = builder.add_node(scalar_node.clone());
    let node2 = builder.add_node(scalar_node.clone());
    let node3 = builder.add_node(scalar_node.clone());

    builder.add_edge(node1, node2, scalar_edge(0), false);
    builder.add_edge(node1, node3, scalar_edge(1), false);
    builder.add_edge(node2, node3, scalar_edge_with_p(2, 1), false);
    builder.add_edge(node2, node3, scalar_edge_with_p(3, 1), false);
    builder.add_external_edge(node1, scalar_edge(4), false, Flow::Sink);
    builder.add_external_edge(node3, scalar_edge(5), false, Flow::Source);

    let hedge_graph = builder.build();
    let uv_graph = UVGraph::from_hedge(hedge_graph);
    println!(
        "{}",
        uv_graph.dot_impl(
            &uv_graph.full_filter(),
            "",
            &|a| Some(format!("label=\"{}/{}\"", a.num, a.den)),
            &|n| Some(format!("label=\"{}\"", n.num))
        )
    );

    let wood = uv_graph.wood(&uv_graph.full_graph());

    // println!("{}", wood.dot(&uv_graph));
    // println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.lmb);
    // assert_eq!(152, ufold.n_terms());
    // ufold.compute(&uv_graph);

    // println!("unfolded : {}", ufold.show_structure(&uv_graph).unwrap());
    // println!("graph: {}", ufold.graphs());

    let result = ufold.expr(&uv_graph, &uv_graph.full_graph()).unwrap().0;

    println!("{:>}", result);

    let result = result
        .replace(function!(GS.emr_mom, W_.x_, W_.y_))
        .with(W_.y_);

    // println!("{}", ufold.structure_and_res(&uv_graph));
    println!("{:>}", result);

    let t = symbol!("t");
    let series = Atom::var(t).npow(4)
        * result
            .replace(parse!("K(3)"))
            .with(parse!("t*K(3)"))
            .replace(parse!("spenso::dot(t*K(x_),y_)"))
            .repeat()
            .with(parse!("t*spenso::dot(K(x_),y_)"));

    println!("SERIESINPUT:\n{:>}", series);

    let s = series
        .replace(t)
        .with(Atom::var(t).npow(-1))
        .series(t, Atom::Zero, 0.into(), true)
        .unwrap();
    println!("Series: {:>}", s);
    println!(
        "Correct UV cancellation if 0: {:>}",
        s.to_atom().expand().factor()
    );

    panic!("STOP");

    let exp = result
        .replace(parse!("spenso::dot(k_(x_),l_(y_))"))
        .with(parse!(
            "k(x_,0)*k(y_,0)-k(x_,1)*k(y_,1)-k(x_,2)*k(y_,2)-k(x_,3)*k(y_,3)"
        ));

    let mut fnmap = FunctionMap::new();

    fnmap.add_constant(
        Atom::var(symbol!("m")),
        symbolica::domains::float::Complex::new((1.).into(), (0.).into()),
    );
    fnmap.add_constant(
        Atom::var(symbol!("MH")),
        symbolica::domains::float::Complex::new((1.).into(), (0.).into()),
    );
    fnmap.add_constant(
        Atom::var(symbol!("mUV")),
        symbolica::domains::float::Complex::new((10.).into(), (0.).into()),
    );
    let ev = exp
        .evaluator(
            &fnmap,
            &[
                parse!("k(1, 0)"),
                parse!("k(1, 1)"),
                parse!("k(1, 2)"),
                parse!("k(1, 3)"),
                parse!("k(3, 0)"),
                parse!("k(3, 1)"),
                parse!("k(3, 2)"),
                parse!("k(3, 3)"),
                parse!("k(5, 0)"),
                parse!("k(5, 1)"),
                parse!("k(5, 2)"),
                parse!("k(5, 3)"),
            ],
            OptimizationSettings::default(),
        )
        .unwrap();

    let mut ev2 = ev.map_coeff(&|x| x.re.to_f64());

    println!("Single limit:");
    for t in (0..100_000).step_by(5000) {
        let r = (t as f64).powf(4.)
            * ev2.evaluate_single(&[
                3.,
                43.,
                5.,
                6.5,
                t as f64 + 1.,
                t as f64 * 2. + 2.,
                t as f64 + 3.,
                t as f64 + 4.,
                7.,
                4.2,
                1.,
                2.4,
            ]);
        println!("{} {}", r, r.abs().log10());
    }

    println!("Full graph limit:");
    for t in (0..100_000).step_by(5000) {
        let r = (t as f64).powf(8.)
            * ev2.evaluate_single(&[
                t as f64 + 3.,
                t as f64 + 43.,
                t as f64 * 3. + 5.,
                t as f64 + 6.5,
                t as f64 + 1.,
                t as f64 * 2. + 2.,
                t as f64 + 3.,
                t as f64 + 4.,
                7.,
                4.2,
                1.,
                2.4,
            ]);
        println!("{} {}", r, r.abs().log10());
    }
}

#[test]
fn nested_bubble_scalar() {
    let scalar_node = UVNode {
        dod: 0,
        num: Atom::num(1),
        color: None,
    };

    fn scalar_edge(eid: i32) -> UVEdge {
        let model = load_generic_model("sm");
        let higgs = model.get_particle("H");
        let m2 = parse!(higgs.mass.name).npow(2);
        UVEdge {
            og_edge: 1, // not needed
            dod: -2,
            particle: higgs,
            num: Atom::num(1),
            den: spenso_lor_atom(eid, 1, GS.dim).npow(2).to_dots() - m2,
        }
    }

    let mut builder = HedgeGraphBuilder::new();

    let node1 = builder.add_node(scalar_node.clone());

    let node2 = builder.add_node(scalar_node.clone());

    let node3 = builder.add_node(scalar_node.clone());

    builder.add_edge(node1, node2, scalar_edge(0), false);

    builder.add_edge(node1, node3, scalar_edge(1), false);

    builder.add_edge(node2, node3, scalar_edge(2), false);
    builder.add_edge(node2, node3, scalar_edge(3), false);

    builder.add_external_edge(node1, scalar_edge(4), false, Flow::Sink);

    builder.add_external_edge(node3, scalar_edge(5), false, Flow::Source);

    let hedge_graph = builder.build();
    let uv_graph = UVGraph::from_hedge(hedge_graph);
    println!(
        "{}",
        uv_graph.dot_impl(
            &uv_graph.full_filter(),
            "",
            &|a| Some(a.den.to_string()),
            &|n| Some(n.num.to_string())
        )
    );

    let wood = uv_graph.wood(&uv_graph.full_graph());

    //println!("{}", wood.dot(&uv_graph));
    //println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.lmb);
    // assert_eq!(152, ufold.n_terms());
    // ufold.compute(&uv_graph);

    //println!("unfolded : {}", ufold.show_structure(&uv_graph).unwrap());
    //println!("graph: {}", ufold.graphs());

    let result = ufold.expr(&uv_graph, &uv_graph.full_graph()).unwrap().0;

    let result = result
        .replace(function!(GS.emr_mom, W_.x_, W_.y_))
        .with(W_.y_);

    println!(
        "{}",
        ufold.structure_and_res(&uv_graph, &uv_graph.full_graph())
    );
    println!("RESULT {:>}", result);

    let t = symbol!("t");
    let series = Atom::var(t).npow(4)
        * result
            .replace(parse!("K(3)"))
            .with(parse!("t*K(3)"))
            .replace(parse!("spenso::dot(t*K(x_),y_)"))
            .repeat()
            .with(parse!("t*spenso::dot(K(x_),y_)"));

    let s = series
        .replace(t)
        .with(Atom::var(t).npow(-1))
        .series(t, Atom::Zero, 0.into(), true)
        .unwrap();
    println!("Series: {}", s);
    println!("Correct UV cancellation if 0: {:>}", s.to_atom().expand());

    let exp = result
        .replace(parse!("spenso::dot(k_(x_),l_(y_))"))
        .with(parse!(
            "k(x_,0)*k(y_,0)-k(x_,1)*k(y_,1)-k(x_,2)*k(y_,2)-k(x_,3)*k(y_,3)"
        ));

    let mut fnmap = FunctionMap::new();

    fnmap.add_constant(
        Atom::var(symbol!("m")),
        symbolica::domains::float::Complex::new((1.).into(), (0.).into()),
    );
    fnmap.add_constant(
        Atom::var(symbol!("MH")),
        symbolica::domains::float::Complex::new((1.).into(), (0.).into()),
    );
    fnmap.add_constant(
        Atom::var(symbol!("mUV")),
        symbolica::domains::float::Complex::new((10.).into(), (0.).into()),
    );
    let ev = exp
        .evaluator(
            &fnmap,
            &[
                parse!("k(1, 0)"),
                parse!("k(1, 1)"),
                parse!("k(1, 2)"),
                parse!("k(1, 3)"),
                parse!("k(3, 0)"),
                parse!("k(3, 1)"),
                parse!("k(3, 2)"),
                parse!("k(3, 3)"),
                parse!("k(5, 0)"),
                parse!("k(5, 1)"),
                parse!("k(5, 2)"),
                parse!("k(5, 3)"),
            ],
            OptimizationSettings::default(),
        )
        .unwrap();

    let mut ev2 = ev.map_coeff(&|x| x.re.to_f64());

    for t in (0..100_000).step_by(5000) {
        let r = (t as f64).powf(4.)
            * ev2.evaluate_single(&[
                3.,
                43.,
                5.,
                6.5,
                t as f64 + 1.,
                t as f64 * 2. + 2.,
                t as f64 + 3.,
                t as f64 + 4.,
                7.,
                4.2,
                1.,
                2.4,
            ]);
        println!("{} {}", r, r.abs().log10());
    }
}

#[test]
fn disconnect_forest_scalar() {
    let scalar_node = UVNode {
        dod: 0,
        num: Atom::num(1),
        color: None,
    };

    fn scalar_edge(eid: i32) -> UVEdge {
        let model = load_generic_model("sm");
        let higgs = model.get_particle("H");
        let m2 = parse!(higgs.mass.name).npow(2);
        UVEdge {
            og_edge: 1, // not needed
            dod: -2,
            particle: higgs,
            num: Atom::num(1),
            den: spenso_lor_atom(eid, 1, GS.dim).npow(2).to_dots() - m2,
        }
    }

    let mut builder = HedgeGraphBuilder::new();

    let node1 = builder.add_node(scalar_node.clone());

    let node2 = builder.add_node(scalar_node.clone());

    let node3 = builder.add_node(scalar_node.clone());
    let node4 = builder.add_node(scalar_node.clone());

    builder.add_edge(node1, node2, scalar_edge(0), false);

    builder.add_edge(node1, node4, scalar_edge(1), false);
    builder.add_edge(node1, node4, scalar_edge(2), false);

    builder.add_edge(node2, node3, scalar_edge(3), false);
    builder.add_edge(node2, node3, scalar_edge(4), false);

    builder.add_edge(node3, node4, scalar_edge(5), false);

    builder.add_external_edge(node1, scalar_edge(6), false, Flow::Sink);

    builder.add_external_edge(node3, scalar_edge(7), false, Flow::Source);

    let hedge_graph = builder.build();
    let uv_graph = UVGraph::from_hedge(hedge_graph);
    println!(
        "{}",
        uv_graph.dot_impl(
            &uv_graph.full_filter(),
            "",
            &|a| Some(a.den.to_string()),
            &|n| Some(n.num.to_string())
        )
    );

    let wood = uv_graph.wood(&uv_graph.full_graph());

    println!("{}", wood.dot(&uv_graph));
    println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.lmb);
    // assert_eq!(152, ufold.n_terms());
    // ufold.compute(&uv_graph);

    println!(
        "unfolded : {}",
        ufold
            .show_structure(&uv_graph, &uv_graph.full_graph())
            .unwrap()
    );
    println!("graph: {}", ufold.graphs());

    let result = ufold.expr(&uv_graph, &uv_graph.full_graph()).unwrap().0;

    println!(
        "{}",
        ufold.structure_and_res(&uv_graph, &uv_graph.full_graph())
    );
    println!("{:>}", result);

    let exp = result
        .replace(parse!("spenso::dot(k_(x_),l_(y_))"))
        .with(parse!(
            "k(x_,0)*k(y_,0)-k(x_,1)*k(y_,1)-k(x_,2)*k(y_,2)-k(x_,3)*k(y_,3)"
        ));

    let mut fnmap = FunctionMap::new();

    fnmap.add_constant(
        Atom::var(symbol!("m")),
        symbolica::domains::float::Complex::new((1.).into(), (0.).into()),
    );
    let ev = exp
        .evaluator(
            &fnmap,
            &[
                parse!("k(2, 0)"),
                parse!("k(2, 1)"),
                parse!("k(2, 2)"),
                parse!("k(2, 3)"),
                parse!("k(4, 0)"),
                parse!("k(4, 1)"),
                parse!("k(4, 2)"),
                parse!("k(4, 3)"),
                parse!("k(5, 0)"),
                parse!("k(5, 1)"),
                parse!("k(5, 2)"),
                parse!("k(5, 3)"),
                parse!("k(7, 0)"),
                parse!("k(7, 1)"),
                parse!("k(7, 2)"),
                parse!("k(7, 3)"),
            ],
            OptimizationSettings::default(),
        )
        .unwrap();

    let mut ev2 = ev.map_coeff(&|x| x.re.to_f64());

    println!("Single limit");
    for t in (0..100_000).step_by(5000) {
        let r = (t as f64).powf(4.)
            * ev2.evaluate_single(&[
                3.,
                43.,
                5.,
                6.5,
                t as f64 + 1.,
                t as f64 * 2. + 2.,
                t as f64 + 3.,
                t as f64 + 4.,
                7.,
                4.2,
                1.,
                2.4,
                6.5,
                8.6,
                3.4,
                2.1,
            ]);
        println!("{} {}", r, r.abs().log10());
    }

    println!("Single disjoint limit");

    for t in (0..100_000).step_by(5000) {
        let r = (t as f64).powf(4.)
            * ev2.evaluate_single(&[
                t as f64 + 3.,
                t as f64 + 43.,
                t as f64 * 3. + 5.,
                t as f64 + 6.5,
                t as f64 + 1.,
                t as f64 * 2. + 2.,
                t as f64 + 3.,
                t as f64 + 4.,
                7.,
                4.2,
                1.,
                2.4,
                6.5,
                8.6,
                3.4,
                2.1,
            ]);
        println!("{} {}", r, r.abs().log10());
    }

    println!("Full graph limit");

    for t in (0..100_000).step_by(5000) {
        let r = (t as f64).powf(8.)
            * ev2.evaluate_single(&[
                t as f64 + 3.,
                t as f64 + 43.,
                t as f64 * 3. + 5.,
                t as f64 + 6.5,
                t as f64 + 1.,
                t as f64 * 2. + 2.,
                t as f64 + 3.,
                t as f64 + 4.,
                t as f64 * 5. + 7.,
                t as f64 + 4.2,
                t as f64 + 1.,
                t as f64 * 2. + 2.4,
                6.5,
                8.6,
                3.4,
                2.1,
            ]);
        println!("{} {}", r, r.abs().log10());
    }
}

    */

// #[test]
// #[allow(unused)]
// fn easy() {
//     let model = load_generic_model("sm");
//     let mut symbolica_graph = symbolica::graph::Graph::new();

//     let dummy_external_vertex_rule = ArcVertexRule(Arc::new(VertexRule {
//         name: "external".into(),
//         couplings: vec![],
//         lorentz_structures: vec![],
//         particles: vec![],
//         color_structures: ColorStructure::new(vec![]),
//     }));

//     let incoming = symbolica_graph.add_node(NodeColorWithVertexRule {
//         external_tag: 1,
//         vertex_rule: dummy_external_vertex_rule.clone(),
//     });

//     let outgoing = symbolica_graph.add_node(NodeColorWithVertexRule {
//         external_tag: 2,
//         vertex_rule: dummy_external_vertex_rule.clone(),
//     });

//     let tth = NodeColorWithVertexRule {
//         external_tag: 0,
//         vertex_rule: model.get_vertex_rule("V_141"),
//     };

//     let t = EdgeColor::from_particle(model.get_particle("t"));
//     let h = EdgeColor::from_particle(model.get_particle("H"));

//     let l1 = symbolica_graph.add_node(tth.clone());
//     let l2 = symbolica_graph.add_node(tth.clone());

//     let e1 = symbolica_graph.add_edge(incoming, l1, false, h).unwrap();
//     let e2 = symbolica_graph.add_edge(l2, outgoing, false, h).unwrap();
//     symbolica_graph.add_edge(l1, l2, true, t);
//     symbolica_graph.add_edge(l2, l1, true, t);

//     let bare_graph = BareGraph::from_symbolica_graph(
//         &model,
//         "1l_prop".into(),
//         &symbolica_graph,
//         Atom::num(1),
//         vec![((Some(1), Some(2)))],
//         None,
//     )
//     .unwrap();

//     let uv_graph = UVGraph::from_graph(&bare_graph);

//     println!("{}", uv_graph.base_dot());

//     let wood = uv_graph.wood(&uv_graph.full_graph());

//     println!("{}", wood.dot(&uv_graph));
//     println!("{}", wood.show_graphs(&uv_graph));

//     let mut ufold = wood.unfold(&uv_graph, &uv_graph.lmb);
//     // assert_eq!(152, ufold.n_terms());
//     // ufold.compute(&uv_graph);

//     println!(
//         "unfolded : {}",
//         ufold
//             .show_structure(&uv_graph, &uv_graph.full_graph())
//             .unwrap()
//     );
//     println!("graph: {}", ufold.graphs());

//     // let structure = wood.unfold(&uv_graph,&uv_graph.cut_edges);

//     // println!("{}", structure.show_structure(&wood, &uv_graph));
//     // println!("{}", structure.n_elements());
// }

// #[test]
// #[allow(unused)]
// fn tbt() {
//     let (model, amplitude, _) =
//         load_amplitude_output("TEST_AMPLITUDE_triangle_box_triangle_phys/GL_OUTPUT", true);

//     let mut graph = amplitude.amplitude_graphs[0].graph.clone();

//     // graph.generate_uv();

//     graph.generate_loop_momentum_bases();

//     let uv_graph = UVGraph::from_graph(&graph.bare_graph);

//     println!("tbt_dot{}", uv_graph.base_dot());

//     let wood = uv_graph.wood(&uv_graph.full_graph());

//     println!("{}", wood.dot(&uv_graph));

//     let mut ufold = wood.unfold(&uv_graph, &uv_graph.lmb);
//     // ufold.compute(&uv_graph);

//     println!(
//         "unfolded : {}",
//         ufold
//             .show_structure(&uv_graph, &uv_graph.full_graph())
//             .unwrap()
//     );
//     println!("graph: {}", ufold.graphs());

//     // let structure = wood.unfold(&uv_graph,&uv_graph.cut_edges);

//     // println!("{}", structure.show_structure(&wood, &uv_graph));
//     // println!("{}", structure.n_elements());
// }

// #[test]
// #[allow(unused)]
// fn bugblatter_forest() {
//     // println!("{}", env!("CARGO_CRATE_NAME"));
//     let model = load_generic_model("sm");
//     println!("{}", model.vertex_rules[0].0.name);
//     let mut symbolica_graph = symbolica::graph::Graph::new();

//     let ttg = NodeColorWithVertexRule {
//         external_tag: 0,
//         vertex_rule: model.get_vertex_rule("V_137"),
//     };

//     let t = EdgeColor::from_particle(model.get_particle("t"));
//     let g = EdgeColor::from_particle(model.get_particle("g"));

//     let l1 = symbolica_graph.add_node(ttg.clone());
//     let l2 = symbolica_graph.add_node(ttg.clone());
//     let l3 = symbolica_graph.add_node(ttg.clone());
//     let l4 = symbolica_graph.add_node(ttg.clone());
//     let l5 = symbolica_graph.add_node(ttg.clone());
//     let l6 = symbolica_graph.add_node(ttg.clone());

//     symbolica_graph.add_edge(l1, l2, true, t);
//     symbolica_graph.add_edge(l2, l3, true, t);
//     symbolica_graph.add_edge(l3, l1, true, t);

//     symbolica_graph.add_edge(l4, l5, true, t);
//     symbolica_graph.add_edge(l5, l6, true, t);
//     symbolica_graph.add_edge(l6, l4, true, t);

//     symbolica_graph.add_edge(l1, l4, true, g);
//     symbolica_graph.add_edge(l2, l5, true, g);
//     symbolica_graph.add_edge(l3, l6, true, g);

//     let bare_graph = BareGraph::from_symbolica_graph(
//         &model,
//         "bugblatter".into(),
//         &symbolica_graph,
//         Atom::num(1),
//         vec![],
//         None,
//     )
//     .unwrap();

//     let uv_graph = UVGraph::from_graph(&bare_graph);

//     println!("{}", uv_graph.base_dot());

//     let wood = uv_graph.wood(&uv_graph.full_graph());

//     assert_eq!(20, wood.n_spinneys());
//     println!("{}", wood.dot(&uv_graph));
//     println!("{}", wood.show_graphs(&uv_graph));

//     let mut ufold = wood.unfold(&uv_graph, &uv_graph.lmb);
//     assert_eq!(152, ufold.n_terms());
//     // ufold.compute(&uv_graph);

//     println!(
//         "unfolded : {}",
//         ufold
//             .show_structure(&uv_graph, &uv_graph.full_graph())
//             .unwrap()
//     );
//     println!("graph: {}", ufold.graphs());

//     // let structure = wood.unfold(&uv_graph,&uv_graph.cut_edges);

//     // println!("{}", structure.show_structure(&wood, &uv_graph));
//     // println!("{}", structure.n_elements());
// }

// #[test]
// #[allow(unused)]
// fn kaapo_triplering() {
//     let model = load_generic_model("scalars");
//     let mut symbolica_graph = symbolica::graph::Graph::new();

//     let four_scalar = NodeColorWithVertexRule {
//         external_tag: 0,
//         vertex_rule: model.get_vertex_rule("V_4_SCALAR_0000"),
//     };

//     let scalar = EdgeColor::from_particle(model.get_particle("scalar_0"));

//     let l1 = symbolica_graph.add_node(four_scalar.clone());
//     let l2 = symbolica_graph.add_node(four_scalar.clone());
//     let l3 = symbolica_graph.add_node(four_scalar.clone());

//     symbolica_graph.add_edge(l1, l2, true, scalar.clone());
//     symbolica_graph.add_edge(l2, l3, true, scalar.clone());
//     symbolica_graph.add_edge(l3, l1, true, scalar.clone());
//     symbolica_graph.add_edge(l1, l2, true, scalar.clone());
//     symbolica_graph.add_edge(l2, l3, true, scalar.clone());
//     symbolica_graph.add_edge(l3, l1, true, scalar.clone());

//     let bare_graph = BareGraph::from_symbolica_graph(
//         &model,
//         "threeringscalar".into(),
//         &symbolica_graph,
//         Atom::num(1),
//         vec![],
//         None,
//     )
//     .unwrap();

//     let uv_graph = UVGraph::from_graph(&bare_graph);

//     // println!("{}", uv_graph.base_dot());

//     let wood = uv_graph.wood(&uv_graph.full_graph());
//     assert_eq!(26, wood.n_spinneys());

//     // println!("{}", wood.dot(&uv_graph));
//     // println!("{}", wood.show_graphs(&uv_graph));

//     let mut ufold = wood.unfold(&uv_graph, &uv_graph.lmb);
//     // ufold.compute(&uv_graph);

//     println!(
//         "unfolded : {}",
//         ufold
//             .show_structure(&uv_graph, &uv_graph.full_graph())
//             .unwrap()
//     );
//     println!("graph: {}", ufold.graphs());

//     assert_eq!(242, ufold.n_terms());
// }

// #[test]
// #[allow(unused)]
// fn kaapo_quintic_scalar() {
//     let model = load_generic_model("scalars");
//     let mut symbolica_graph = symbolica::graph::Graph::new();

//     let three = NodeColorWithVertexRule {
//         external_tag: 0,
//         vertex_rule: model.get_vertex_rule("V_3_SCALAR_000"),
//     };
//     let four = NodeColorWithVertexRule {
//         external_tag: 0,
//         vertex_rule: model.get_vertex_rule("V_4_SCALAR_0000"),
//     };
//     let five = NodeColorWithVertexRule {
//         external_tag: 0,
//         vertex_rule: model.get_vertex_rule("V_5_SCALAR_00000"),
//     };

//     let scalar = EdgeColor::from_particle(model.get_particle("scalar_0"));

//     let l1 = symbolica_graph.add_node(three);
//     let l2 = symbolica_graph.add_node(four);
//     let l3 = symbolica_graph.add_node(five);

//     symbolica_graph.add_edge(l1, l2, true, scalar.clone());
//     symbolica_graph.add_edge(l2, l3, true, scalar.clone());
//     symbolica_graph.add_edge(l3, l1, true, scalar.clone());
//     symbolica_graph.add_edge(l3, l2, true, scalar.clone());
//     symbolica_graph.add_edge(l2, l3, true, scalar.clone());
//     symbolica_graph.add_edge(l3, l1, true, scalar.clone());

//     let bare_graph = BareGraph::from_symbolica_graph(
//         &model,
//         "threeringscalar".into(),
//         &symbolica_graph,
//         Atom::num(1),
//         vec![],
//         None,
//     )
//     .unwrap();

//     let uv_graph = UVGraph::from_graph(&bare_graph);

//     println!("{}", uv_graph.base_dot());

//     let wood = uv_graph.wood(&uv_graph.full_graph());

//     assert_eq!(25, wood.n_spinneys());

//     // println!("{}", wood.dot(&uv_graph));
//     // println!("{}", wood.show_graphs(&uv_graph));

//     let mut ufold = wood.unfold(&uv_graph, &uv_graph.lmb);
//     // ufold.compute(&uv_graph);

//     println!(
//         "unfolded : {}",
//         ufold
//             .show_structure(&uv_graph, &uv_graph.full_graph())
//             .unwrap()
//     );
//     println!("graph: {}", ufold.graphs());

//     assert_eq!(248, ufold.n_terms());
// }

use super::*;

use std::cmp::Ordering;
use std::collections::HashSet;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct TestNode(&'static str);

#[allow(clippy::non_canonical_partial_ord_impl)]
// Implement PartialOrd and Ord for TestNode to define the partial order
impl PartialOrd for TestNode {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self.0, other.0) {
            ("A", "A") => Some(Ordering::Equal),
            ("B", "B") => Some(Ordering::Equal),
            ("C", "C") => Some(Ordering::Equal),
            ("D", "D") => Some(Ordering::Equal),
            ("A", _) => Some(Ordering::Less), // A < B, A < C, A < D
            (_, "A") => Some(Ordering::Greater), // B > A, C > A, D > A
            ("B", "D") => Some(Ordering::Less), // B < D
            ("D", "B") => Some(Ordering::Greater), // D > B
            ("C", "D") => Some(Ordering::Less), // C < D
            ("D", "C") => Some(Ordering::Greater), // D > C
            _ => None,                        // Other pairs are incomparable
        }
    }
}

impl Ord for TestNode {
    fn cmp(&self, other: &Self) -> Ordering {
        // Since our partial order can be undefined for some pairs (None), we unwrap safely here
        // because we know our test cases only involve comparable nodes.
        self.partial_cmp(other).unwrap()
    }
}

const LU_TEST_SETTINGS: &'static str = "
General:
  amplitude_prefactor:
    im: 1.0
    re: 0.0
  debug: 0
  force_orientations: null
  joint_numerator_eval: true
  load_compiled_cff: true
  load_compiled_numerator: true
  load_compiled_separate_orientations: false
  use_ltd: false
Integrand:
  type: gamma_loop
Integrator:
  bin_number_evolution: null
  continuous_dim_learning_rate: 0.0
  discrete_dim_learning_rate: 0.0
  integrated_phase: real
  max_prob_ratio: 1000.0
  min_samples_for_update: 100
  n_bins: 16
  n_increase: 0
  n_max: 100000000
  n_start: 20000
  seed: 2
  show_max_wgt_info: false
  train_on_avg: false
Kinematics:
  e_cm: 800.0
  externals:
    data:
      helicities: []
      momenta:
      - - 800.0
        - 0.0
        - 0.0
        - 0.0
    type: constant
Observables: []
Selectors: []
Stability:
  levels:
  - escalate_for_large_weight_threshold: 0.9
    precision: Double
    required_precision_for_im: 1.0e-07
    required_precision_for_re: 1.0e-07
  - escalate_for_large_weight_threshold: -1.0
    precision: Quad
    required_precision_for_im: 1.0e-10
    required_precision_for_re: 1.0e-10
  rotate_numerator: false
  rotation_axis:
  - type: euler_angles
    alpha: 0.1
    beta: 0.2
    gamma: 0.3
sampling:
  sample_orientations: false
  sampling_type:
    b: 1.0
    mapping: linear
    mode: spherical
    subtype: default
  type: discrete_graph_sampling
subtraction:
  integrated_ct_settings:
    range:
      h_function_settings:
        enabled_dampening: true
        function: poly_exponential
        power: null
        sigma: 1.0
      type: infinite
  local_ct_settings:
    dampen_integrable_singularity:
      type: exponential
    uv_localisation:
      dynamic_width: false
      gaussian_width: 1.0
      sliver_width: 10.0
  overlap_settings:
    check_global_center: true
    force_global_center: null
    try_origin: false
    try_origin_all_lmbs: false";

#[test]
fn quick_test() {
    fn expose_mass_dimension(atom: &Atom) -> Atom {
        atom.replace(function!(GS.ose, W_.x___))
            .with(parse!("E"))
            .replace(function!(GS.dot, W_.x_, W_.y_))
            .with(parse!("-E^2"))
            .expand()
    }

    let edge_id = 4;

    let cff = parse!(
        "
        1 / 8 * (OSE(2) + OSE(3))
            ^ -1 * (OSE(2) + OSE(4))
            ^ -1 * OSE(2)
            ^ -1 * OSE(3)
            ^ -1 * OSE(4)
            ^ -1"
    );

    let num_2 = parse!("(-OSE(2) * OSE(4) + dot(Q3(3), Q3(3))) * OSE(3)^2");

    let mut expr = &num_2 * &cff;

    expr = expr
        .replace(function!(GS.ose, edge_id, W_.x___))
        .with(function!(GS.ose, edge_id, W_.x___, Atom::var(GS.rescale)))
        .derivative(GS.rescale)
        .replace(function!(Atom::DERIVATIVE, W_.x___, W_.x_))
        .with(Atom::var(W_.x_).npow(-1) / 2)
        .replace(function!(GS.ose, edge_id, W_.x___, Atom::var(GS.rescale)))
        .with(function!(GS.ose, edge_id, W_.x___));

    expr = expr.replace(function!(GS.ose, W_.x___)).with(parse!("E"));

    println!("expr: {}", expr.expand());
    expr = expose_mass_dimension(&expr);
    println!("expr: {}", expr);

    let new_expr = num_2 * cff;
    let mut alt_der = new_expr
        .replace(function!(GS.ose, edge_id))
        .with(parse!("OSE4"))
        .derivative(symbol!("OSE4"))
        / parse!("2*OSE4");

    alt_der = alt_der
        .replace(symbol!("OSE4"))
        .with(function!(GS.ose, edge_id));

    alt_der = alt_der
        .replace(function!(GS.ose, W_.x___))
        .with(parse!("E"));

    println!("alt_der: {}", alt_der.expand());
    alt_der = expose_mass_dimension(&alt_der);
    println!("alt_der: {}", alt_der);

    let test_expr = (parse!("1/(8*(OSE3+ OSE4))").derivative(symbol!("OSE4")) / parse!("2*OSE4"))
        .replace(parse!("OSE4"))
        .with(parse!("OSE3"))
        .expand();

    println!("test_expr: {}", test_expr);

    let test_expr_2 = (parse!("OSE4/(8*OSE4*(OSE3+ OSE4))").derivative(symbol!("OSE4"))
        / parse!("2*OSE4"))
    .replace(parse!("OSE4"))
    .with(parse!("OSE3"))
    .expand();

    println!("test_expr: {}", test_expr_2);
}
