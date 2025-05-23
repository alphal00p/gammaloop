use std::{sync::Arc, time::Instant};

use ahash::HashMap;
use linnet::half_edge::hedgevec::HedgeVec;
use nalgebra::LU;
use pathfinding::num_traits::real;
use smartstring::SmartString;
use spenso::{
    network::{library::TensorLibraryData, parsing::ShadowedStructure},
    structure::{
        abstract_index::AbstractIndex,
        dimension::Dimension,
        representation::{Minkowski, RepName},
        NamedStructure, ToSymbolic,
    },
};
use symbolica::evaluate::{FunctionMap, OptimizationSettings};
use symbolica_community::physics::algebraic_simplification::metric::MetricSimplifier;
use typed_index_collections::TiVec;

use crate::{
    cff::cut_expression::SuperGraphOrientationID,
    feyngen::{
        diagram_generator::{EdgeColor, FeynGen, NodeColorWithVertexRule},
        FeynGenFilters,
    },
    graph::InteractionVertexInfo,
    integrands::Integrand,
    integrate::{havana_integrate, UserData},
    model::{ArcVertexRule, ColorStructure, VertexRule},
    momentum_sample::ExternalIndex,
    new_cs::{CrossSection, CrossSectionGraph, CutId, ProcessDefinition},
    new_graph::{get_cff_inverse_energy_product_impl, Edge, ExternalConnection, Graph, Vertex},
    numerator::UnInit,
    signature::LoopExtSignature,
    tests_from_pytest::{load_amplitude_output, load_generic_model},
    uv::UVGraph,
    Settings,
};

#[test]
fn double_triangle_LU() {
    let _ = env_logger::builder().is_test(true).try_init();

    let with_log_uv = true;

    let model = load_generic_model("sm");
    let mut underlying = HedgeGraphBuilder::new();

    let hhh = VertexInfo::InteractonVertexInfo(InteractionVertexInfo {
        vertex_rule: model.get_vertex_rule("V_9"),
    });
    let htt = VertexInfo::InteractonVertexInfo(InteractionVertexInfo {
        vertex_rule: model.get_vertex_rule("V_141"),
    });

    let hprop = model.get_propagator("H_propFeynman");
    let hp = model.get_particle("H");

    let tprop = model.get_propagator("t_propFeynman");
    let tp = model.get_particle("t");

    let n1 = underlying.add_node(Vertex {
        name: "n1".into(),
        vertex_info: hhh,
        dod: 0,
        num: Atom::one(),
    });
    let n2 = underlying.add_node(Vertex {
        name: "n2".into(),
        vertex_info: htt.clone(),
        dod: 0,

        num: Atom::one(),
    });
    let n3 = underlying.add_node(Vertex {
        name: "n3".into(),
        vertex_info: htt.clone(),
        dod: 0,
        num: Atom::one(),
    });
    let n4 = underlying.add_node(Vertex {
        name: "n4".into(),
        vertex_info: htt.clone(),
        dod: 0,

        num: Atom::one(),
    });

    underlying.add_edge(
        n1,
        n2,
        Edge {
            name: "e0".into(),
            edge_type: EdgeType::Virtual,
            particle: hp.clone(),
            propagator: hprop.clone(),
            internal_index: vec![],
            dod: -2,
            num: Atom::one(),
        },
        false,
    );

    underlying.add_edge(
        n1,
        n3,
        Edge {
            name: "e1".into(),
            edge_type: EdgeType::Virtual,
            particle: hp.clone(),
            propagator: hprop.clone(),
            internal_index: vec![],
            dod: -2,
            num: Atom::one(),
        },
        false,
    );

    underlying.add_edge(
        n2,
        n3,
        Edge {
            name: "e2".into(),
            edge_type: EdgeType::Virtual,
            particle: tp.clone(),
            propagator: tprop.clone(),
            internal_index: vec![],
            dod: -2,
            num: Atom::one(),
        },
        true,
    );

    underlying.add_edge(
        n3,
        n4,
        Edge {
            name: "e3".into(),
            edge_type: EdgeType::Virtual,
            particle: tp.clone(),
            propagator: tprop.clone(),
            internal_index: vec![],
            dod: if with_log_uv { -1 } else { -2 },
            num: if with_log_uv {
                spenso_lor_atom(3, 1, GS.dim)
            } else {
                Atom::one()
            },
        },
        true,
    );

    underlying.add_edge(
        n4,
        n2,
        Edge {
            name: "e4".into(),
            edge_type: EdgeType::Virtual,
            particle: tp.clone(),
            propagator: tprop.clone(),
            internal_index: vec![],
            dod: if with_log_uv { -1 } else { -2 },
            num: if with_log_uv {
                spenso_lor_atom(4, 1, GS.dim)
            } else {
                Atom::one()
            },
        },
        true,
    );

    underlying.add_external_edge(
        n1,
        Edge {
            name: "q1".into(),
            edge_type: EdgeType::Incoming,
            particle: hp.clone(),
            propagator: hprop.clone(),
            internal_index: vec![],
            dod: 0,
            num: Atom::one(),
        },
        false,
        Flow::Sink,
    );

    underlying.add_external_edge(
        n4,
        Edge {
            name: "q2".into(),
            edge_type: EdgeType::Outgoing,
            particle: hp.clone(),
            propagator: hprop.clone(),
            internal_index: vec![],
            dod: 0,
            num: Atom::one(),
        },
        false,
        Flow::Source,
    );

    let underlying = underlying.build();

    let mut loop_momentum_basis = LoopMomentumBasis {
        tree: None,
        basis: vec![EdgeIndex::from(0), EdgeIndex::from(4)].into(),
        edge_signatures: underlying
            .new_hedgevec(|_, _, _| LoopExtSignature::from((vec![], vec![]))),
    };

    loop_momentum_basis
        .set_edge_signatures(&underlying)
        .unwrap();

    let graph = Graph {
        multiplicity: Atom::one(),
        name: "DT".into(),
        underlying,
        loop_momentum_basis,
        vertex_slots: vec![].into(),
        external_connections: None,
    };

    let mut cs: CrossSectionGraph<UnInit> = CrossSectionGraph::new(graph);

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

    let super_uv_graph = UVGraph::from_underlying(&cs.graph.underlying);
    let orientation_id = SuperGraphOrientationID(0);
    let supergraph_orientation_data = &cs
        .derived_data
        .cff_expression
        .as_ref()
        .unwrap()
        .orientation_data[orientation_id];

    let mut cut_atoms: TiVec<CutId, Atom> = TiVec::new();

    for (id, c) in cs.cuts.iter_enumerated() {
        let esurface_id = cs.cut_esurface_id_map[id];
        let cff_cut_expr = &cs
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .cut_expressions[id];

        if let Some((left_orientation, right_orientation)) =
            cff_cut_expr.orientation_map.get_lr_or(orientation_id)
        {
            let left_orientation_data =
                &cff_cut_expr.left_amplitude.orientations[left_orientation].data;
            let right_orientation_data =
                &cff_cut_expr.right_amplitude.orientations[right_orientation].data;

            let cut_mom_basis_id = cs.derived_data.esurface_data.as_ref().unwrap()[esurface_id]
                .as_ref()
                .unwrap()
                .cut_momentum_basis;
            let cut_lmb = &cs.derived_data.lmbs.as_ref().unwrap()[cut_mom_basis_id];

            let mut left_wood = super_uv_graph.wood(&c.left);

            let mut left_forest = left_wood.unfold(&super_uv_graph, &super_uv_graph.cut_edges);
            left_forest.compute(&super_uv_graph);
            left_forest.compute_cff(&super_uv_graph, left_orientation_data, &None);

            let mut right_forest = super_uv_graph
                .wood(&c.right)
                .unfold(&super_uv_graph, &super_uv_graph.cut_edges);
            right_forest.compute(&super_uv_graph);
            right_forest.compute_cff(&super_uv_graph, right_orientation_data, &None);

            println!("//left: \n{}", super_uv_graph.dot(&c.left));

            println!("//right: \n{}", super_uv_graph.dot(&c.right));
            let left_expr =
                left_forest.local_expr(&super_uv_graph, &c.left, &None, left_orientation_data);
            let right_expr =
                right_forest.local_expr(&super_uv_graph, &c.right, &None, right_orientation_data);

            let mut cut_res = cs.add_additional_factors_to_cff_atom(&(left_expr * right_expr), id);

            // add Feynman rules of cut edges
            for (_p, edge_index, d) in super_uv_graph.iter_edges(&c.cut.left) {
                let edge_id = usize::from(edge_index) as i64;
                let orientation = supergraph_orientation_data.orientation.clone();
                cut_res = (cut_res * &d.data.num)
                    .replace(function!(GS.emr_mom, edge_id, W_.y_))
                    .with_map(move |m| {
                        let index = m.get(W_.y_).unwrap().to_atom();

                        let sign = SignOrZero::from((&orientation[edge_index]).clone()) * 1;

                        function!(GS.ose, edge_id, index) * sign
                            + function!(GS.emr_vec, edge_id, index)
                    });
            }

            // add Feynman rules of external edges
            for (_p, edge_index, d) in super_uv_graph.iter_edges(&super_uv_graph.external_filter())
            {
                let edge_id = usize::from(edge_index) as i64;
                cut_res = (cut_res * &d.data.num)
                    .replace(function!(GS.emr_mom, edge_id, W_.y_))
                    .with_map(move |m| {
                        let index = m.get(W_.y_).unwrap().to_atom();

                        function!(GS.ose, edge_id, index) // NOTE: not OSE for external edge, understood as taking 0th part
                                + function!(GS.emr_vec, edge_id, index)
                    });
            }

            // contract all dot products, set all cross terms ose.q3 to 0
            // MS.dot is a 4d dot product
            cut_res = cut_res
                .expand()
                // .to_dots()//need to allow x___
                .replace(function!(GS.emr_vec, W_.x__, W_.y_).npow(2))
                .with(function!(
                    MS.dot,
                    function!(GS.emr_vec, W_.x__),
                    function!(GS.emr_vec, W_.x__)
                ))
                .replace(
                    function!(GS.emr_vec, W_.x__, W_.a_) * function!(GS.emr_vec, W_.y__, W_.a_),
                )
                .repeat()
                .with(function!(
                    MS.dot,
                    function!(GS.emr_vec, W_.x__),
                    function!(GS.emr_vec, W_.y__)
                ))
                .replace(function!(GS.ose, W_.y__, W_.x_).npow(2))
                .with(function!(GS.ose, W_.y__).npow(2))
                .replace(function!(GS.ose, W_.x__, W_.x_) * function!(GS.ose, W_.y__, W_.x_))
                .repeat()
                .with(function!(GS.ose, W_.x__) * function!(GS.ose, W_.y__))
                .replace(function!(GS.emr_vec, W_.x__, W_.a_) * function!(GS.ose, W_.y__, W_.a_))
                .with(Atom::Zero)
                .replace(function!(
                    MS.dot,
                    function!(GS.emr_vec, W_.x_),
                    function!(GS.ose, W_.y_)
                ))
                .with(Atom::Zero)
                .replace(function!(
                    MS.dot,
                    function!(GS.ose, W_.x_),
                    function!(GS.ose, W_.y_)
                ))
                .with(function!(GS.ose, W_.x_) * function!(GS.ose, W_.y_));

            // substitute all OSEs from subgraphs, they are in the form OSE(edge_id, momentum, mass)
            cut_res = cut_res
                .replace(function!(GS.ose, W_.x_, W_.y_, W_.z_))
                .with(
                    (-function!(
                        MS.dot,
                        function!(GS.emr_vec, W_.y_),
                        function!(GS.emr_vec, W_.y_)
                    ) + W_.z_)
                        .sqrt(),
                );

            // substitute all OSEs (add minus sign to cancel minus sign from 4d dot product)
            for (_p, edge_id, d) in super_uv_graph.iter_edges(&c.cut.left) {
                let mass2 = Atom::new_var(symbol!(d.data.particle.mass.name.as_str())).npow(2);
                cut_res = cut_res
                    .replace(function!(GS.ose, usize::from(edge_id) as i64))
                    .with(
                        (-function!(
                            MS.dot,
                            function!(GS.emr_vec, usize::from(edge_id) as i64),
                            function!(GS.emr_vec, usize::from(edge_id) as i64)
                        ) + mass2)
                            .sqrt(),
                    );
            }

            cut_res = cut_res
                .replace(function!(GS.emr_vec, function!(GS.emr_mom, W_.x_)))
                .with(function!(GS.emr_vec, W_.x_));

            cut_res = cut_res
                .replace(function!(
                    MS.dot,
                    function!(GS.emr_vec, W_.x_),
                    function!(GS.emr_vec, W_.y_)
                ))
                .with(
                    -(function!(GS.emr_vec, W_.x_, 1) * function!(GS.emr_vec, W_.y_, 1)
                        + function!(GS.emr_vec, W_.x_, 2) * function!(GS.emr_vec, W_.y_, 2)
                        + function!(GS.emr_vec, W_.x_, 3) * function!(GS.emr_vec, W_.y_, 3)),
                );

            let cut_res = cut_res.expand();
            println!("Cut {} result: {:>}", id, cut_res);

            if super_uv_graph.dod(&c.right) == 0 {
                // only check when this graph is a UV subgraph
                let t = symbol!("t");
                let series = Atom::new_var(t).npow(3)
                    * cut_res
                        .replace(parse!("Q3(3)").unwrap())
                        .with(parse!("t*Q3(3)").unwrap())
                        .replace(parse!("Q3(2)").unwrap())
                        .with(parse!("t*Q3(3)-Q3(1)").unwrap())
                        .replace(parse!("Q3(4)").unwrap())
                        .with(parse!("t*Q3(3)-Q3(6)").unwrap())
                        .replace(parse!("symbolica_community::dot(t*x_,y_)").unwrap())
                        .repeat()
                        .with(parse!("t*symbolica_community::dot(x_,y_)").unwrap());

                let s = series
                    .replace(t)
                    .with(Atom::new_var(t).npow(-1))
                    .series(t, Atom::Zero, 0.into(), true)
                    .unwrap();

                let mut r = s.to_atom().expand();
                r = r
                    .replace((-Atom::new_var(W_.x_)).pow(Atom::new_num((-5, 2))))
                    .with(Atom::new_var(W_.x_).pow(Atom::new_num((-5, 2))) * Atom::I)
                    .replace((-Atom::new_var(W_.x_)).pow(Atom::new_num((-3, 2))))
                    .with(-Atom::new_var(W_.x_).pow(Atom::new_num((-3, 2))) * Atom::I)
                    .expand(); // help Symbolica with cancellations and avoid bad simplification of (-1)^(-5/2)
                println!("Correct UV cancellation if 0: {:>}", r);
            }

            cut_atoms.push(cut_res);
        } else {
            cut_atoms.push(Atom::new());
        }
    }

    cs.derived_data.bare_cff_evaluators = None;
    cs.build_cut_evaluators(&model, Some(cut_atoms));

    let external_connection = ExternalConnection {
        incoming_index: ExternalIndex::from(0),
        outgoing_index: ExternalIndex::from(1),
    };

    let cs_struct = CrossSection {
        name: "LU_test".into(),
        supergraphs: vec![cs],
        external_connections: vec![external_connection],
        external_particles: vec![hp.clone(), hp.clone()],
        n_incmoming: 1,
    };

    let settings: Settings = serde_yaml::from_str(LU_TEST_SETTINGS).unwrap();
    let integrand = cs_struct.generate_integrand(settings.clone(), &model);

    let result = match integrand {
        Integrand::NewIntegrand(real_integrand) => havana_integrate(
            &settings,
            |set| real_integrand.user_data_generator(1, set),
            None,
            None,
            None,
        ),
        _ => unimplemented!(),
    };

    //println!("Final result: {:>}", sum.expand());
}

#[test]
fn nested_bubble_soft_ct() {
    // let massive_particle = ArcParticle(Arc::new(Particle::))

    let scalar_node = UVNode {
        dod: 0,
        num: Atom::new_num(1),
        color: None,
    };

    fn scalar_edge(eid: i32) -> UVEdge {
        let model = load_generic_model("sm");
        let higgs = model.get_particle("H");
        let m2 = parse!(higgs.mass.name).unwrap().npow(2);
        UVEdge {
            og_edge: 1, // not needed
            dod: -2,
            particle: higgs,
            num: Atom::new_num(1),
            den: spenso_lor_atom(eid, 1, GS.dim).npow(2).to_dots() - m2,
        }
    }

    fn scalar_edge_with_p(eid: i32, ind: impl Into<AbstractIndex>) -> UVEdge {
        let model = load_generic_model("sm");
        let higgs = model.get_particle("H");
        let m2 = parse!(higgs.mass.name).unwrap().npow(2);
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

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.cut_edges);
    // assert_eq!(152, ufold.n_terms());
    ufold.compute(&uv_graph);

    // println!("unfolded : {}", ufold.show_structure(&uv_graph).unwrap());
    // println!("graph: {}", ufold.graphs());

    let result = ufold.expr(&uv_graph, &uv_graph.full_graph()).unwrap().0;

    let result = result
        .replace(function!(GS.emr_mom, W_.x_, W_.y_))
        .with(W_.y_);

    // println!("{}", ufold.structure_and_res(&uv_graph));
    println!("{:>}", result);

    let t = symbol!("t");
    let series = Atom::new_var(t).npow(4)
        * result
            .replace(parse!("K(3)").unwrap())
            .with(parse!("t*K(3)").unwrap())
            .replace(parse!("symbolica_community::dot(t*K(x_),y_)").unwrap())
            .repeat()
            .with(parse!("t*symbolica_community::dot(K(x_),y_)").unwrap());

    println!("SERIESINPUT:\n{:>}", series);

    let s = series
        .replace(t)
        .with(Atom::new_var(t).npow(-1))
        .series(t, Atom::Zero, 0.into(), true)
        .unwrap();
    println!("Series: {:>}", s);
    println!("Correct UV cancellation if 0: {:>}", s.to_atom().expand());

    let exp = result
        .replace(parse!("symbolica_community::dot(k_(x_),l_(y_))").unwrap())
        .with(parse!("k(x_,0)*k(y_,0)-k(x_,1)*k(y_,1)-k(x_,2)*k(y_,2)-k(x_,3)*k(y_,3)").unwrap());

    let mut fnmap = FunctionMap::new();

    fnmap.add_constant(Atom::new_var(symbol!("m")), (1.).into());
    fnmap.add_constant(Atom::new_var(symbol!("MH")), (1.).into());
    fnmap.add_constant(Atom::new_var(symbol!("mUV")), (10.).into());
    fnmap.add_constant(Atom::new_var(symbol!("ZERO")), (0).into());
    let ev = exp
        .evaluator(
            &fnmap,
            &[
                parse!("k(1, 0)").unwrap(),
                parse!("k(1, 1)").unwrap(),
                parse!("k(1, 2)").unwrap(),
                parse!("k(1, 3)").unwrap(),
                parse!("k(3, 0)").unwrap(),
                parse!("k(3, 1)").unwrap(),
                parse!("k(3, 2)").unwrap(),
                parse!("k(3, 3)").unwrap(),
                parse!("k(5, 0)").unwrap(),
                parse!("k(5, 1)").unwrap(),
                parse!("k(5, 2)").unwrap(),
                parse!("k(5, 3)").unwrap(),
                parse!("k(6, 0)").unwrap(),
                parse!("k(6, 1)").unwrap(),
                parse!("k(6, 2)").unwrap(),
                parse!("k(6, 3)").unwrap(),
            ],
            OptimizationSettings::default(),
        )
        .unwrap();

    let mut ev2 = ev.map_coeff(&|x| x.to_f64());

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
        num: Atom::new_num(1),
        color: None,
    };

    fn scalar_edge(eid: i32) -> UVEdge {
        let model = load_generic_model("sm");
        let higgs = model.get_particle("H");
        let m2 = parse!(higgs.mass.name).unwrap().npow(2);
        UVEdge {
            og_edge: 1, // not needed
            dod: -2,
            particle: higgs,
            num: Atom::new_num(1),
            den: spenso_lor_atom(eid, 1, GS.dim).npow(2).to_dots() - m2,
        }
    }

    fn scalar_edge_with_p(eid: i32, ind: impl Into<AbstractIndex>) -> UVEdge {
        let model = load_generic_model("sm");
        let higgs = model.get_particle("H");
        let m2 = parse!(higgs.mass.name).unwrap().npow(2);
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

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.cut_edges);
    // assert_eq!(152, ufold.n_terms());
    ufold.compute(&uv_graph);

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
    let series = Atom::new_var(t).npow(4)
        * result
            .replace(parse!("K(3)").unwrap())
            .with(parse!("t*K(3)").unwrap())
            .replace(parse!("symbolica_community::dot(t*K(x_),y_)").unwrap())
            .repeat()
            .with(parse!("t*symbolica_community::dot(K(x_),y_)").unwrap());

    println!("SERIESINPUT:\n{:>}", series);

    let s = series
        .replace(t)
        .with(Atom::new_var(t).npow(-1))
        .series(t, Atom::Zero, 0.into(), true)
        .unwrap();
    println!("Series: {:>}", s);
    println!(
        "Correct UV cancellation if 0: {:>}",
        s.to_atom().expand().factor()
    );

    panic!("STOP");

    let exp = result
        .replace(parse!("symbolica_community::dot(k_(x_),l_(y_))").unwrap())
        .with(parse!("k(x_,0)*k(y_,0)-k(x_,1)*k(y_,1)-k(x_,2)*k(y_,2)-k(x_,3)*k(y_,3)").unwrap());

    let mut fnmap = FunctionMap::new();

    fnmap.add_constant(Atom::new_var(symbol!("m")), (1.).into());
    fnmap.add_constant(Atom::new_var(symbol!("MH")), (1.).into());
    fnmap.add_constant(Atom::new_var(symbol!("mUV")), (10.).into());
    let ev = exp
        .evaluator(
            &fnmap,
            &[
                parse!("k(1, 0)").unwrap(),
                parse!("k(1, 1)").unwrap(),
                parse!("k(1, 2)").unwrap(),
                parse!("k(1, 3)").unwrap(),
                parse!("k(3, 0)").unwrap(),
                parse!("k(3, 1)").unwrap(),
                parse!("k(3, 2)").unwrap(),
                parse!("k(3, 3)").unwrap(),
                parse!("k(5, 0)").unwrap(),
                parse!("k(5, 1)").unwrap(),
                parse!("k(5, 2)").unwrap(),
                parse!("k(5, 3)").unwrap(),
            ],
            OptimizationSettings::default(),
        )
        .unwrap();

    let mut ev2 = ev.map_coeff(&|x| x.to_f64());

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
        num: Atom::new_num(1),
        color: None,
    };

    fn scalar_edge(eid: i32) -> UVEdge {
        let model = load_generic_model("sm");
        let higgs = model.get_particle("H");
        let m2 = parse!(higgs.mass.name).unwrap().npow(2);
        UVEdge {
            og_edge: 1, // not needed
            dod: -2,
            particle: higgs,
            num: Atom::new_num(1),
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

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.cut_edges);
    // assert_eq!(152, ufold.n_terms());
    ufold.compute(&uv_graph);

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
    let series = Atom::new_var(t).npow(4)
        * result
            .replace(parse!("K(3)").unwrap())
            .with(parse!("t*K(3)").unwrap())
            .replace(parse!("symbolica_community::dot(t*K(x_),y_)").unwrap())
            .repeat()
            .with(parse!("t*symbolica_community::dot(K(x_),y_)").unwrap());

    let s = series
        .replace(t)
        .with(Atom::new_var(t).npow(-1))
        .series(t, Atom::Zero, 0.into(), true)
        .unwrap();
    println!("Series: {}", s);
    println!("Correct UV cancellation if 0: {:>}", s.to_atom().expand());

    let exp = result
        .replace(parse!("symbolica_community::dot(k_(x_),l_(y_))").unwrap())
        .with(parse!("k(x_,0)*k(y_,0)-k(x_,1)*k(y_,1)-k(x_,2)*k(y_,2)-k(x_,3)*k(y_,3)").unwrap());

    let mut fnmap = FunctionMap::new();

    fnmap.add_constant(Atom::new_var(symbol!("m")), (1.).into());
    fnmap.add_constant(Atom::new_var(symbol!("MH")), (1.).into());
    fnmap.add_constant(Atom::new_var(symbol!("mUV")), (10.).into());
    let ev = exp
        .evaluator(
            &fnmap,
            &[
                parse!("k(1, 0)").unwrap(),
                parse!("k(1, 1)").unwrap(),
                parse!("k(1, 2)").unwrap(),
                parse!("k(1, 3)").unwrap(),
                parse!("k(3, 0)").unwrap(),
                parse!("k(3, 1)").unwrap(),
                parse!("k(3, 2)").unwrap(),
                parse!("k(3, 3)").unwrap(),
                parse!("k(5, 0)").unwrap(),
                parse!("k(5, 1)").unwrap(),
                parse!("k(5, 2)").unwrap(),
                parse!("k(5, 3)").unwrap(),
            ],
            OptimizationSettings::default(),
        )
        .unwrap();

    let mut ev2 = ev.map_coeff(&|x| x.to_f64());

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
        num: Atom::new_num(1),
        color: None,
    };

    fn scalar_edge(eid: i32) -> UVEdge {
        let model = load_generic_model("sm");
        let higgs = model.get_particle("H");
        let m2 = parse!(higgs.mass.name).unwrap().npow(2);
        UVEdge {
            og_edge: 1, // not needed
            dod: -2,
            particle: higgs,
            num: Atom::new_num(1),
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

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.cut_edges);
    // assert_eq!(152, ufold.n_terms());
    ufold.compute(&uv_graph);

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
        .replace(parse!("symbolica_community::dot(k_(x_),l_(y_))").unwrap())
        .with(parse!("k(x_,0)*k(y_,0)-k(x_,1)*k(y_,1)-k(x_,2)*k(y_,2)-k(x_,3)*k(y_,3)").unwrap());

    let mut fnmap = FunctionMap::new();

    fnmap.add_constant(Atom::new_var(symbol!("m")), (1.).into());
    let ev = exp
        .evaluator(
            &fnmap,
            &[
                parse!("k(2, 0)").unwrap(),
                parse!("k(2, 1)").unwrap(),
                parse!("k(2, 2)").unwrap(),
                parse!("k(2, 3)").unwrap(),
                parse!("k(4, 0)").unwrap(),
                parse!("k(4, 1)").unwrap(),
                parse!("k(4, 2)").unwrap(),
                parse!("k(4, 3)").unwrap(),
                parse!("k(5, 0)").unwrap(),
                parse!("k(5, 1)").unwrap(),
                parse!("k(5, 2)").unwrap(),
                parse!("k(5, 3)").unwrap(),
                parse!("k(7, 0)").unwrap(),
                parse!("k(7, 1)").unwrap(),
                parse!("k(7, 2)").unwrap(),
                parse!("k(7, 3)").unwrap(),
            ],
            OptimizationSettings::default(),
        )
        .unwrap();

    let mut ev2 = ev.map_coeff(&|x| x.to_f64());

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

#[test]
#[allow(unused)]
fn easy() {
    let model = load_generic_model("sm");
    let mut symbolica_graph = symbolica::graph::Graph::new();

    let dummy_external_vertex_rule = ArcVertexRule(Arc::new(VertexRule {
        name: "external".into(),
        couplings: vec![],
        lorentz_structures: vec![],
        particles: vec![],
        color_structures: ColorStructure::new(vec![]),
    }));

    let incoming = symbolica_graph.add_node(NodeColorWithVertexRule {
        external_tag: 1,
        vertex_rule: dummy_external_vertex_rule.clone(),
    });

    let outgoing = symbolica_graph.add_node(NodeColorWithVertexRule {
        external_tag: 2,
        vertex_rule: dummy_external_vertex_rule.clone(),
    });

    let tth = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_141"),
    };

    let t = EdgeColor::from_particle(model.get_particle("t"));
    let h = EdgeColor::from_particle(model.get_particle("H"));

    let l1 = symbolica_graph.add_node(tth.clone());
    let l2 = symbolica_graph.add_node(tth.clone());

    let e1 = symbolica_graph.add_edge(incoming, l1, false, h).unwrap();
    let e2 = symbolica_graph.add_edge(l2, outgoing, false, h).unwrap();
    symbolica_graph.add_edge(l1, l2, true, t);
    symbolica_graph.add_edge(l2, l1, true, t);

    let bare_graph = BareGraph::from_symbolica_graph(
        &model,
        "1l_prop".into(),
        &symbolica_graph,
        Atom::new_num(1),
        vec![((Some(1), Some(2)))],
        None,
    )
    .unwrap();

    let uv_graph = UVGraph::from_graph(&bare_graph);

    println!("{}", uv_graph.base_dot());

    let wood = uv_graph.wood(&uv_graph.full_graph());

    println!("{}", wood.dot(&uv_graph));
    println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.cut_edges);
    // assert_eq!(152, ufold.n_terms());
    ufold.compute(&uv_graph);

    println!(
        "unfolded : {}",
        ufold
            .show_structure(&uv_graph, &uv_graph.full_graph())
            .unwrap()
    );
    println!("graph: {}", ufold.graphs());

    // let structure = wood.unfold(&uv_graph,&uv_graph.cut_edges);

    // println!("{}", structure.show_structure(&wood, &uv_graph));
    // println!("{}", structure.n_elements());
}

#[test]
#[allow(unused)]
fn tbt() {
    let (model, amplitude, _) =
        load_amplitude_output("TEST_AMPLITUDE_triangle_box_triangle_phys/GL_OUTPUT", true);

    let mut graph = amplitude.amplitude_graphs[0].graph.clone();

    // graph.generate_uv();

    graph.generate_loop_momentum_bases();

    let uv_graph = UVGraph::from_graph(&graph.bare_graph);

    println!("tbt_dot{}", uv_graph.base_dot());

    let wood = uv_graph.wood(&uv_graph.full_graph());

    println!("{}", wood.dot(&uv_graph));

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.cut_edges);
    ufold.compute(&uv_graph);

    println!(
        "unfolded : {}",
        ufold
            .show_structure(&uv_graph, &uv_graph.full_graph())
            .unwrap()
    );
    println!("graph: {}", ufold.graphs());

    // let structure = wood.unfold(&uv_graph,&uv_graph.cut_edges);

    // println!("{}", structure.show_structure(&wood, &uv_graph));
    // println!("{}", structure.n_elements());
}

#[test]
#[allow(unused)]
fn bugblatter_forest() {
    // println!("{}", env!("CARGO_CRATE_NAME"));
    let model = load_generic_model("sm");
    println!("{}", model.vertex_rules[0].0.name);
    let mut symbolica_graph = symbolica::graph::Graph::new();

    let ttg = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_137"),
    };

    let t = EdgeColor::from_particle(model.get_particle("t"));
    let g = EdgeColor::from_particle(model.get_particle("g"));

    let l1 = symbolica_graph.add_node(ttg.clone());
    let l2 = symbolica_graph.add_node(ttg.clone());
    let l3 = symbolica_graph.add_node(ttg.clone());
    let l4 = symbolica_graph.add_node(ttg.clone());
    let l5 = symbolica_graph.add_node(ttg.clone());
    let l6 = symbolica_graph.add_node(ttg.clone());

    symbolica_graph.add_edge(l1, l2, true, t);
    symbolica_graph.add_edge(l2, l3, true, t);
    symbolica_graph.add_edge(l3, l1, true, t);

    symbolica_graph.add_edge(l4, l5, true, t);
    symbolica_graph.add_edge(l5, l6, true, t);
    symbolica_graph.add_edge(l6, l4, true, t);

    symbolica_graph.add_edge(l1, l4, true, g);
    symbolica_graph.add_edge(l2, l5, true, g);
    symbolica_graph.add_edge(l3, l6, true, g);

    let bare_graph = BareGraph::from_symbolica_graph(
        &model,
        "bugblatter".into(),
        &symbolica_graph,
        Atom::new_num(1),
        vec![],
        None,
    )
    .unwrap();

    let uv_graph = UVGraph::from_graph(&bare_graph);

    println!("{}", uv_graph.base_dot());

    let wood = uv_graph.wood(&uv_graph.full_graph());

    assert_eq!(20, wood.n_spinneys());
    println!("{}", wood.dot(&uv_graph));
    println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.cut_edges);
    assert_eq!(152, ufold.n_terms());
    ufold.compute(&uv_graph);

    println!(
        "unfolded : {}",
        ufold
            .show_structure(&uv_graph, &uv_graph.full_graph())
            .unwrap()
    );
    println!("graph: {}", ufold.graphs());

    // let structure = wood.unfold(&uv_graph,&uv_graph.cut_edges);

    // println!("{}", structure.show_structure(&wood, &uv_graph));
    // println!("{}", structure.n_elements());
}

#[test]
#[allow(unused)]
fn kaapo_triplering() {
    let model = load_generic_model("scalars");
    let mut symbolica_graph = symbolica::graph::Graph::new();

    let four_scalar = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_4_SCALAR_0000"),
    };

    let scalar = EdgeColor::from_particle(model.get_particle("scalar_0"));

    let l1 = symbolica_graph.add_node(four_scalar.clone());
    let l2 = symbolica_graph.add_node(four_scalar.clone());
    let l3 = symbolica_graph.add_node(four_scalar.clone());

    symbolica_graph.add_edge(l1, l2, true, scalar.clone());
    symbolica_graph.add_edge(l2, l3, true, scalar.clone());
    symbolica_graph.add_edge(l3, l1, true, scalar.clone());
    symbolica_graph.add_edge(l1, l2, true, scalar.clone());
    symbolica_graph.add_edge(l2, l3, true, scalar.clone());
    symbolica_graph.add_edge(l3, l1, true, scalar.clone());

    let bare_graph = BareGraph::from_symbolica_graph(
        &model,
        "threeringscalar".into(),
        &symbolica_graph,
        Atom::new_num(1),
        vec![],
        None,
    )
    .unwrap();

    let uv_graph = UVGraph::from_graph(&bare_graph);

    // println!("{}", uv_graph.base_dot());

    let wood = uv_graph.wood(&uv_graph.full_graph());
    assert_eq!(26, wood.n_spinneys());

    // println!("{}", wood.dot(&uv_graph));
    // println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.cut_edges);
    ufold.compute(&uv_graph);

    println!(
        "unfolded : {}",
        ufold
            .show_structure(&uv_graph, &uv_graph.full_graph())
            .unwrap()
    );
    println!("graph: {}", ufold.graphs());

    assert_eq!(242, ufold.n_terms());
}

#[test]
#[allow(unused)]
fn kaapo_quintic_scalar() {
    let model = load_generic_model("scalars");
    let mut symbolica_graph = symbolica::graph::Graph::new();

    let three = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_3_SCALAR_000"),
    };
    let four = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_4_SCALAR_0000"),
    };
    let five = NodeColorWithVertexRule {
        external_tag: 0,
        vertex_rule: model.get_vertex_rule("V_5_SCALAR_00000"),
    };

    let scalar = EdgeColor::from_particle(model.get_particle("scalar_0"));

    let l1 = symbolica_graph.add_node(three);
    let l2 = symbolica_graph.add_node(four);
    let l3 = symbolica_graph.add_node(five);

    symbolica_graph.add_edge(l1, l2, true, scalar.clone());
    symbolica_graph.add_edge(l2, l3, true, scalar.clone());
    symbolica_graph.add_edge(l3, l1, true, scalar.clone());
    symbolica_graph.add_edge(l3, l2, true, scalar.clone());
    symbolica_graph.add_edge(l2, l3, true, scalar.clone());
    symbolica_graph.add_edge(l3, l1, true, scalar.clone());

    let bare_graph = BareGraph::from_symbolica_graph(
        &model,
        "threeringscalar".into(),
        &symbolica_graph,
        Atom::new_num(1),
        vec![],
        None,
    )
    .unwrap();

    let uv_graph = UVGraph::from_graph(&bare_graph);

    println!("{}", uv_graph.base_dot());

    let wood = uv_graph.wood(&uv_graph.full_graph());

    assert_eq!(25, wood.n_spinneys());

    // println!("{}", wood.dot(&uv_graph));
    // println!("{}", wood.show_graphs(&uv_graph));

    let mut ufold = wood.unfold(&uv_graph, &uv_graph.cut_edges);
    ufold.compute(&uv_graph);

    println!(
        "unfolded : {}",
        ufold
            .show_structure(&uv_graph, &uv_graph.full_graph())
            .unwrap()
    );
    println!("graph: {}", ufold.graphs());

    assert_eq!(248, ufold.n_terms());
}

use super::*;

use std::cmp::Ordering;

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
  continuous_dim_learning_rate: 0.3
  discrete_dim_learning_rate: 1.5
  integrated_phase: real
  max_prob_ratio: 1000.0
  min_samples_for_update: 100
  n_bins: 16
  n_increase: 0
  n_max: 100000000
  n_start: 2000000
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
  rotation_axis: []
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
