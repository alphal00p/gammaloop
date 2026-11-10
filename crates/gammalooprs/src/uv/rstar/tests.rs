use std::{collections::BTreeMap, sync::Arc};

use eyre::{Result, ensure, eyre};
use feynkit_generator::{
    GenerationFilter, GenerationOptions, NumeratorGrouping, Process, VertexSelector,
};
use idenso::{
    CookMode, CookSettings, IndexTooling,
    color::ColorSimplifySettings,
    representations::ColorAdjoint,
    tensor::{AlgebraContraction, AlgebraSettings, SymbolicTensor},
};
use linnet::half_edge::{
    NodeIndex,
    involution::{EdgeIndex, HedgePair},
};
use spenso::{
    network::{library::symbolic::ETS, tags::SPENSO_TAG},
    structure::{
        abstract_index::AbstractIndex,
        representation::{LibrarySlot, Minkowski, RepName},
        slot::IsAbstractSlot,
    },
};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};
use vakint::Vakint;

use crate::{
    debug_tags,
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    finalized_runtime_dot,
    graph::{Graph, parse::IntoFinalizedRuntimeGraph},
    initialisation::test_initialise,
    model::Model,
    utils::{GS, load_generic_model},
    uv::{
        ApproximationType, RenormalizationPrescriptionSettings, UVgenerationSettings,
        approx::integrated::Integrated, settings::VakintSettings as GammaVakintSettings,
    },
};

use super::RStarPlan;

pub(super) fn settings() -> UVgenerationSettings {
    UVgenerationSettings {
        renormalization_prescription: RenormalizationPrescriptionSettings {
            log_divergent: ApproximationType::PolePart,
            massive_power_divergent: ApproximationType::PolePart,
            massless_power_divergent: ApproximationType::PolePart,
            ..Default::default()
        },
        vakint: GammaVakintSettings {
            evaluation_methods: vec!["rustred".into(), "forcer".into()],
            normalization: "FMFTandMATAD".into(),
            additional_normalization: "1".into(),
            matad: crate::uv::settings::MATADSettings {
                direct_numerical_substition: false,
                substitute_hpls: false,
                ..Default::default()
            },
            ..Default::default()
        },
        ..Default::default()
    }
}

fn external_vertices(graph: &Graph) -> Vec<NodeIndex> {
    graph
        .iter_edges()
        .filter_map(|(pair, _, _)| match pair {
            HedgePair::Unpaired { hedge, .. } => Some(graph.node_id(hedge)),
            _ => None,
        })
        .collect()
}

pub(super) fn gluon_self_energy(
    model: &Model,
    loops: usize,
    vertices: usize,
    rule: &str,
    accept: impl Fn(&Graph) -> bool,
) -> Result<Graph> {
    let options = GenerationOptions::default()
        .with_loop_count(loops, loops)?
        .max_vertices(vertices)
        .allow_self_loops(false)
        .projector(Atom::one())
        .numerator_grouping(NumeratorGrouping::None)
        .with_graph_filter(GenerationFilter::VertexAllow(vec![VertexSelector::Name(
            rule.into(),
        )]))
        .with_graph_filter(GenerationFilter::LoopOneParticleIrreducible)
        .threads(1);
    let generated =
        Process::new([21_i64], [21_i64]).generate_diagrams(Arc::new(model.clone()), &options)?;
    for diagram in generated.diagrams {
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let external = external_vertices(&graph);
        if external.len() == 2 && external[0] != external[1] && accept(&graph) {
            return Ok(graph);
        }
    }
    Err(eyre!(
        "no {loops}-loop gluon self-energy with vertex {rule}"
    ))
}

pub(super) fn parallel_pair(graph: &Graph) -> Option<[EdgeIndex; 2]> {
    let mut endpoints = BTreeMap::new();
    for (pair, edge, _) in graph.iter_edges() {
        let HedgePair::Paired { source, sink } = pair else {
            continue;
        };
        let a = graph.node_id(source);
        let b = graph.node_id(sink);
        if let Some(previous) = endpoints.insert((a.min(b), a.max(b)), edge) {
            return Some([previous, edge]);
        }
    }
    None
}

pub(super) fn canonical_tensor(expression: &Atom) -> Result<Atom> {
    // Keep symbolic representation dimensions and all open tensor slots. The
    // production counterterm already performed its scalar Laurent projection.
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let tensor = SymbolicTensor::infer_indexed(cooking.try_cook(expression.as_view())?)?
        .undo_trace()?
        .undo_chain()?
        .simplify_algebra(&AlgebraSettings {
            color: Some(ColorSimplifySettings::default()),
            ..Default::default()
        })?
        .to_dots()?;
    Ok(cooking.uncook(tensor.expression().as_view()).expand())
}

pub(super) fn assert_same_counterterm(actual: &Atom, expected: &Atom) -> Result<()> {
    let actual = canonical_tensor(actual)?;
    let expected = canonical_tensor(expected)?;
    let difference = canonical_tensor(&(&actual - &expected))?;
    assert!(
        difference.together().is_zero(),
        "R*: {actual}; reference: {expected}; difference: {difference}"
    );
    Ok(())
}

pub(super) fn massive_oracle(
    plan: &RStarPlan,
    graph: &Graph,
    model: &Model,
    vakint: &Vakint,
    settings: &UVgenerationSettings,
) -> Result<Atom> {
    let mut settings = settings.clone();
    settings.vakint.evaluation_methods = vec!["matad".into()];
    let mut forest = plan.uv_forest(graph, &settings);
    forest.integrate(graph, model, vakint, &settings)?;
    Ok(forest
        .counterterm_for_support(&plan.support)?
        .replace(GS.m_uv_expansion)
        .with(GS.m_uv_vacuum))
}

#[test]
fn one_loop_cubic_gluon_tensor_matches_massive_r() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    let graph = gluon_self_energy(&model, 1, 2, "V_36", |_| true)?;
    let settings = settings();
    let vakint = Vakint::new()?;
    let marked = graph
        .iter_edges()
        .find_map(|(pair, edge, _)| pair.is_paired().then_some(edge))
        .unwrap();
    let plan = RStarPlan::new(&graph, &model, marked, &settings)?;
    assert_eq!(plan.taylor_order, 2);
    let actual = plan.counterterm(&graph, &model, &vakint, &settings)?;
    let expected = massive_oracle(&plan, &graph, &model, &vakint, &settings)?;
    assert!(!canonical_tensor(&expected)?.is_zero());
    assert_same_counterterm(&actual, &expected)
}

#[test]
fn projected_counterterm_contracts_only_the_root_after_taylor() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    let graph = gluon_self_energy(&model, 2, 4, "V_36", |graph| parallel_pair(graph).is_some())?;
    let settings = settings();
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&graph, &model, parallel_pair(&graph).unwrap()[0], &settings)?;
    assert!(!plan.counterterms.is_empty());
    assert!(!plan.ir.is_empty());
    let tensor = plan.counterterm(&graph, &model, &vakint, &settings)?;
    assert!(!tensor.is_zero());
    let cooking = CookSettings::indices()
        .with_mode(CookMode::ReversibleEncoding)
        .with_representation_payloads(true, true);
    let cooked = cooking.try_cook(tensor.as_view())?;
    let ports = cooked
        .list_dangling::<AbstractIndex>()?
        .iter()
        .map(|port| LibrarySlot::<AbstractIndex>::try_from(port.as_view()))
        .collect::<std::result::Result<Vec<_>, _>>()?;
    let lorentz = ports
        .iter()
        .filter(|port| port.rep_name() == (Minkowski {}).into())
        .collect::<Vec<_>>();
    let color = ports
        .iter()
        .filter(|port| port.rep_name() == (ColorAdjoint {}).into())
        .collect::<Vec<_>>();
    assert_eq!(lorentz.len(), 2);
    assert_eq!(color.len(), 2);
    assert_eq!(ports.len(), 4);
    let external = graph
        .iter_edges()
        .find_map(|(pair, edge, _)| (!pair.is_paired()).then_some(edge))
        .unwrap();
    let momentum = function!(
        GS.emr_mom,
        usize::from(external),
        Minkowski {}.new_rep(GS.dim).to_symbolic([])
    );
    let p_squared = function!(SPENSO_TAG.dot, &momentum, &momentum);
    let metric = ETS.metric(lorentz[0].to_atom(), lorentz[1].to_atom())
        * ETS.metric(color[0].to_atom(), color[1].to_atom());
    let normalization = p_squared * (Atom::var(GS.dim) - Atom::one());
    // Fix the original root ports. Applying this closure to a child would
    // leave unmatched ports. Its inverse p^2 is a scalar spectator of outer
    // Taylor and vacuum tensor reduction, and its 1/(d-1) tests evanescent
    // terms multiplying double poles. Keep the genuinely normalized map
    // throughout integration to exercise both tensor-bridge directions.
    let project = |expression: &Atom| -> Result<Atom> {
        if expression.is_zero() {
            return Ok(Atom::Zero);
        }
        let source = cooking.try_cook(expression.as_view())?;
        let contracted = SymbolicTensor::infer_indexed(source * &metric)?
            .undo_trace()?
            .undo_chain()?
            .simplify_algebra(&AlgebraSettings {
                gamma: Some(Default::default()),
                color: Some(ColorSimplifySettings::default()),
                contract: AlgebraContraction::Dots,
                ..Default::default()
            })?;
        ensure!(
            contracted
                .expression()
                .list_dangling::<AbstractIndex>()?
                .is_empty(),
            "the projector was applied to a boundary other than the root"
        );
        Ok(cooking.uncook(contracted.expression().as_view()) / &normalization)
    };
    let actual = plan.projected_counterterm(&graph, &model, &vakint, &settings, &project)?;
    let expected = Integrated::dimensionally_regularized(&project(&tensor)?)
        .series(GS.dim_epsilon, 0, -1)?
        .to_atom();
    assert_same_counterterm(&actual, &expected)
}

#[test]
fn two_loop_cubic_gluon_power_divergent_hard_moments_match_massive_r() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    // AB, AC, CD, CD, DB: the CD gluon bubble is quadratically divergent.
    // Marking CD requires its hard Taylor polynomial through degree two,
    // including tensor moments inserted coefficientwise into the quotient.
    let graph = gluon_self_energy(&model, 2, 4, "V_36", |graph| parallel_pair(graph).is_some())?;
    let bubble = parallel_pair(&graph).unwrap();
    let external = external_vertices(&graph);
    let outside = graph
        .iter_edges()
        .find_map(|(pair, edge, _)| {
            let HedgePair::Paired { source, sink } = pair else {
                return None;
            };
            (external.contains(&graph.node_id(source)) && external.contains(&graph.node_id(sink)))
                .then_some(edge)
        })
        .unwrap();
    let settings = settings();
    let vakint = Vakint::new()?;
    let inside = RStarPlan::new(&graph, &model, bubble[0], &settings)?;
    let outside = RStarPlan::new(&graph, &model, outside, &settings)?;
    assert_eq!(inside.taylor_order, 2);
    assert!(inside.ir.iter().any(|hard| hard.taylor_order == 2));
    assert!(
        inside
            .counterterms
            .iter()
            .any(|child| child.argument.loops == 1 && child.taylor_order == 2)
    );
    let inside_result = inside.counterterm(&graph, &model, &vakint, &settings)?;
    let outside_result = outside.counterterm(&graph, &model, &vakint, &settings)?;
    debug_tags!(#uv, #integrated, #rstar;
        stage = "two_loop_cubic_gluon_hard_line_comparison",
        log.inside_counterterm = canonical_tensor(&inside_result)?,
        log.outside_counterterm = canonical_tensor(&outside_result)?,
        "Two-loop cubic gluon counterterms before the hard-line comparison"
    );
    assert_same_counterterm(&inside_result, &outside_result)?;
    let expected = massive_oracle(&inside, &graph, &model, &vakint, &settings)?;
    assert!(!canonical_tensor(&expected)?.is_zero());
    debug_tags!(#uv, #integrated, #rstar;
        stage = "two_loop_cubic_gluon_massive_oracle",
        log.counterterm = canonical_tensor(&expected)?,
        "Two-loop cubic gluon massive R reference"
    );
    assert_same_counterterm(&inside_result, &expected)
}

#[test]
fn two_loop_quadratic_scalar_matches_massive_r() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("scalars");
    let mut graph: Graph = finalized_runtime_dot!(
        digraph sunrise {
            graph [projector=1];
            node [num=1];
            edge [particle="scalar_1" num=1];
            e [style=invis];
            e -> A:0 [id=3 sink="{ufo_order:0}"];
            B:1 -> e [id=4 source="{ufo_order:0}"];
            A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
            A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:2}"];
            A -> B [id=2 lmb_id=1 source="{ufo_order:3}" sink="{ufo_order:3}"];
        }, &model
    )?;
    let edges = graph
        .iter_edges()
        .map(|(_, edge, _)| edge)
        .collect::<Vec<_>>();
    for edge in edges {
        let data = &mut graph.underlying[edge];
        data.particle = data.particle.clone().override_mass(Some(Atom::Zero));
    }
    let settings = settings();
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&graph, &model, EdgeIndex(0), &settings)?;
    assert_eq!(plan.taylor_order, 2);
    assert!(!plan.ir.is_empty());
    let actual = plan.counterterm(&graph, &model, &vakint, &settings)?;
    let expected = massive_oracle(&plan, &graph, &model, &vakint, &settings)?;
    assert!(!expected.is_zero());
    assert_same_counterterm(&actual, &expected)
}

#[test]
fn embedded_massive_banana_counterterm_matches_standalone_polynomial() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("scalars");
    // Mixed orientations reproduce the five-line vacuum banana. Its three-loop
    // child has two boundary incidences of edge 4; its own nested children have
    // several closing lines. Routing just one endpoint corrupts their p² terms.
    let mut embedded: Graph = finalized_runtime_dot!(
        digraph G {
            graph [projector=0];
            node [num=1];
            edge [particle="scalar_1" num=1];
            B -> A [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
            B -> A [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
            A -> B [id=2 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
            A -> B [id=3 lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:3}"];
            B -> A [id=4 lmb_id=3 source="{ufo_order:4}" sink="{ufo_order:4}"];
        }, &model
    )?;
    let mut standalone: Graph = finalized_runtime_dot!(
        digraph G {
            graph [projector=0];
            node [num=1];
            edge [particle="scalar_1" num=1];
            e [style=invis];
            e -> A:0 [id=4 sink="{ufo_order:4}"];
            B:1 -> e [id=5 source="{ufo_order:4}"];
            A -> B [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
            A -> B [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
            A -> B [id=2 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
            A -> B [id=3 lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:3}"];
        }, &model
    )?;
    let mass = Atom::var(symbolica::symbol!("embedded_banana_mass"; Scalar));
    for graph in [&mut embedded, &mut standalone] {
        let edges = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        for edge in edges {
            graph.underlying[edge].particle = graph[edge]
                .particle
                .clone()
                .override_mass(Some(mass.clone()));
        }
    }
    let mut settings = settings();
    settings.generate_integrated = true;
    settings.add_marker = false;
    settings.vakint.evaluation_methods = vec!["rustred".into()];
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&embedded, &model, EdgeIndex(0), &settings)?;
    let child = plan
        .counterterms
        .iter()
        .find(|call| {
            call.argument.loops == 3
                && embedded
                    .iter_edges_of(&call.argument.support)
                    .all(|(_, edge, _)| edge != EdgeIndex(4))
        })
        .expect("the four-line child is UV divergent");
    let children = plan.evaluate_uv_children(&embedded, &model, &mass, &vakint, &settings)?;
    let actual = canonical_tensor(&children[&child.argument.support])?;
    let standalone_plan = RStarPlan::new(&standalone, &model, EdgeIndex(0), &settings)?;
    let mut forest = standalone_plan.uv_forest(&standalone, &settings);
    forest.integrate(&standalone, &model, &vakint, &settings)?;
    let independent = forest
        .counterterm_for_support(&standalone_plan.support)?
        .replace(GS.m_uv_expansion)
        .with(GS.m_uv_vacuum);
    let momentum = function!(GS.emr_mom, 4, Minkowski {}.new_rep(GS.dim).to_symbolic([]));
    let p2 = function!(SPENSO_TAG.dot, &momentum, &momentum);
    let epsilon = Atom::var(GS.dim_epsilon);
    let expected = (&mass.pow(2) * &p2 / Atom::num(2)
        - mass.pow(4) / Atom::num(2)
        - p2.pow(2) / Atom::num(36))
        / &epsilon
        + (-Atom::num(2) * mass.pow(2) * &p2 / Atom::num(3)
            + Atom::num(10) * mass.pow(4) / Atom::num(3))
            / epsilon.pow(2)
        - Atom::num(2) * mass.pow(4) / epsilon.pow(3);
    assert_same_counterterm(&independent, &expected)?;
    assert_same_counterterm(&actual, &expected)?;
    assert!(!actual.contains_symbol(GS.m_uv_vacuum));
    assert!(!actual.contains_symbol(GS.m_uv_expansion));
    Ok(())
}
