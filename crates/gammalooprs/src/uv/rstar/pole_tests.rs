//! Bare-pole references independent of the R* forest reconstruction.

use std::sync::Arc;

use eyre::Result;
use feynkit_generator::{
    GenerationFilter, GenerationOptions, NumeratorGrouping, Process, VertexSelector,
};
use linnet::half_edge::{
    involution::{EdgeIndex, HedgePair},
    subgraph::SuBitGraph,
};
use spenso::{
    network::tags::SPENSO_TAG,
    structure::representation::{Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    function,
    id::Replacement,
    transcendental::zeta,
};
use vakint::Vakint;

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    finalized_runtime_dot,
    graph::{Graph, LMBext, parse::IntoFinalizedRuntimeGraph},
    initialisation::test_initialise,
    utils::{GS, W_, load_generic_model},
    uv::{UVgenerationSettings, UltravioletGraph},
};

use super::{
    RStarPlan,
    tests::{assert_same_counterterm, canonical_tensor},
};

fn settings() -> UVgenerationSettings {
    let mut settings = super::tests::settings();
    // The normalized scalar bubble is (Q²/mu²)^(-eps)/eps. This keeps the
    // independent convolution references below free of measure conventions.
    settings.vakint.normalization = "GScheme".into();
    settings
}

fn make_massless(graph: &mut Graph) {
    let edges = graph
        .iter_edges()
        .map(|(_, edge, _)| edge)
        .collect::<Vec<_>>();
    for edge in edges {
        let data = &mut graph.underlying[edge];
        data.particle = data.particle.clone().override_mass(Some(Atom::Zero));
    }
}

fn external_square(graph: &Graph) -> Atom {
    let external = graph
        .iter_edges()
        .find_map(|(pair, edge, _)| (!pair.is_paired()).then_some(edge))
        .unwrap();
    let momentum = function!(
        GS.emr_mom,
        usize::from(external),
        Minkowski {}.new_rep(GS.dim).to_symbolic([])
    );
    function!(SPENSO_TAG.dot, &momentum, &momentum)
}

#[test]
fn two_loop_nested_bubble_bare_poles_include_the_physical_logarithm() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("scalars");
    // The external legs attach to A and C, so the AB bubble carries an
    // integrated outer momentum. Both outer lines have power one: unlike
    // placing both external legs on A/B, this graph has no physical IR pole.
    let mut graph: Graph = finalized_runtime_dot!(
        digraph G {
            graph [projector=1];
            node [num=1];
            edge [particle="scalar_1" num=1];
            e [style=invis];
            e -> A:0 [id=4 sink="{ufo_order:0}"];
            C:1 -> e [id=5 source="{ufo_order:0}"];
            A -> B [id=0 source="{ufo_order:1}" sink="{ufo_order:0}"];
            A -> B [id=1 lmb_id=0 source="{ufo_order:2}" sink="{ufo_order:1}"];
            A -> C [id=2 source="{ufo_order:3}" sink="{ufo_order:1}"];
            B -> C [id=3 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:2}"];
        }, &model
    )?;
    make_massless(&mut graph);
    let settings = settings();
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
    assert_eq!(graph.n_loops(&plan.support), 2);
    assert_eq!(graph.compute_dod(&plan.support), 0);
    assert!(!plan.ir.is_empty());

    let actual = plan.pole_part(&graph, &model, &vakint, &settings)?;
    let epsilon = Atom::var(GS.dim_epsilon);
    let logarithm = function!(
        Symbol::LOG,
        -external_square(&graph) / Atom::var(GS.mu_r_sq)
    );
    // Integrating the AB bubble first gives
    // G(1,1) G(1+eps,1)/(eps G(1,1))² × (Q²/mu²)^(-2eps).
    // Its first two coefficients are 1/(2eps²) + 1/(2eps).
    let expected = Atom::one() / (Atom::num(2) * epsilon.pow(2))
        + (Atom::num(1) / Atom::num(2) - logarithm) / &epsilon;
    assert!(actual.contains_symbol(Symbol::LOG), "bare poles: {actual}");
    assert_same_counterterm(&actual, &expected)
}

#[test]
fn overall_convergent_bubble_in_kite_has_a_nonlocal_subdivergence_pole() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("scalars");
    // Contracting the CE bubble leaves the finite five-line, two-loop kite.
    // The full graph is IR finite and overall UV convergent, while its CE
    // subintegral has a logarithmic UV pole. No other UV child is divergent.
    let mut graph: Graph = finalized_runtime_dot!(
        digraph G {
            graph [projector=1];
            node [num=1];
            edge [particle="scalar_1" num=1];
            e [style=invis];
            e -> A:0 [id=7 sink="{ufo_order:0}"];
            B:1 -> e [id=8 source="{ufo_order:0}"];
            C -> E [id=0 source="{ufo_order:0}" sink="{ufo_order:0}"];
            C -> E [id=1 lmb_id=0 source="{ufo_order:1}" sink="{ufo_order:1}"];
            A -> C [id=2 source="{ufo_order:1}" sink="{ufo_order:2}"];
            A -> D [id=3 source="{ufo_order:2}" sink="{ufo_order:0}"];
            B -> E [id=4 source="{ufo_order:1}" sink="{ufo_order:2}"];
            B -> D [id=5 lmb_id=1 source="{ufo_order:2}" sink="{ufo_order:1}"];
            C -> D [id=6 lmb_id=2 source="{ufo_order:3}" sink="{ufo_order:2}"];
        }, &model
    )?;
    make_massless(&mut graph);
    let settings = settings();
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&graph, &model, EdgeIndex::from(0), &settings)?;
    assert_eq!(graph.n_loops(&plan.support), 3);
    assert_eq!(graph.compute_dod(&plan.support), -2);
    assert_eq!(plan.counterterms.len(), 1);
    assert!(
        plan.counterterm(&graph, &model, &vakint, &settings)?
            .is_zero()
    );

    let actual = plan.pole_part(&graph, &model, &vakint, &settings)?;
    // The bubble residue is +1. The Euclidean kite is 6 zeta(3)/Q²;
    // seven Minkowski denominators supply an overall minus sign. Since
    // the quotient is finite, this pole has no scale logarithm or double pole.
    let expected =
        Atom::num(6) * function!(zeta(), 3) / (external_square(&graph) * Atom::var(GS.dim_epsilon));
    assert!(!canonical_tensor(&actual)?.is_zero());
    assert!(!actual.contains_symbol(Symbol::LOG), "bare poles: {actual}");
    assert_same_counterterm(&actual, &expected)
}

#[test]
fn one_loop_ghost_bare_pole_matches_the_exact_bubble_moment() -> Result<()> {
    test_initialise()?;
    let model = load_generic_model("sm");
    let options = GenerationOptions::default()
        .with_loop_count(1, 1)?
        .max_vertices(2)
        .allow_self_loops(false)
        .projector(Atom::one())
        .numerator_grouping(NumeratorGrouping::None)
        .with_graph_filter(GenerationFilter::VertexAllow(vec![VertexSelector::Name(
            "V_35".into(),
        )]))
        .with_graph_filter(GenerationFilter::LoopOneParticleIrreducible)
        .threads(1);
    let generated =
        Process::new(["ghG"], ["ghG"]).generate_diagrams(Arc::new(model.clone()), &options)?;
    let graph = generated.diagrams[0].to_gamma_loop_graph(None, true)?;
    let marked = graph
        .iter_edges()
        .find_map(|(pair, edge, _)| pair.is_paired().then_some(edge))
        .unwrap();
    let settings = settings();
    let vakint = Vakint::new()?;
    let plan = RStarPlan::new(&graph, &model, marked, &settings)?;

    // Derive the physical two-point routing directly from graph orientations,
    // independently of RStarIntegrand and pole_part's routing helper.
    let external = graph
        .iter_edges()
        .filter_map(|(pair, edge, _)| match pair {
            HedgePair::Unpaired { flow, .. } => Some((edge, flow)),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(external.len(), 2);
    let (p, p_flow) = external[0];
    let physical = external
        .iter()
        .map(|(edge, flow)| {
            let sign = if *edge == p || *flow != p_flow { 1 } else { -1 };
            Replacement::new(
                function!(GS.emr_mom, usize::from(*edge), W_.x___).to_pattern(),
                Atom::num(sign) * function!(GS.emr_mom, usize::from(p), W_.x___),
            )
        })
        .collect::<Vec<_>>();
    let basis = graph.lmb_of(&plan.support);
    assert_eq!(basis.loop_edges.len(), 1);
    let k = basis.loop_edges.raw[0];
    let other = graph
        .iter_edges_of(&plan.support)
        .find_map(|(_, edge, _)| (edge != k).then_some(edge))
        .unwrap();
    let a = basis
        .loop_atom::<Atom>(other, GS.emr_mom, &[], true)
        .replace(function!(GS.emr_mom, usize::from(k)))
        .with(Atom::one());
    let b = basis
        .ext_atom::<Atom>(other, GS.emr_mom, &[], true)
        .replace_multiple(&physical)
        .replace(function!(GS.emr_mom, usize::from(p)))
        .with(Atom::one());
    assert_eq!(a.pow(2), Atom::one());
    assert_eq!(b.pow(2), Atom::one());
    assert!(basis.ext_atom::<Atom>(k, GS.emr_mom, &[], true).is_zero());

    // With denominators k² (a k+b p)², symmetry gives <k>=-b p/(2a).
    // The genuine ghost-gluon vertices give a numerator linear in k, so
    // substituting this first moment integrates its complete color tensor,
    // coupling and external derivative exactly. The scalar bubble pole is 1/eps.
    let numerator = graph
        .numerator(&plan.support, &graph.empty_subgraph::<SuBitGraph>())
        .to_d_dim(GS.dim)
        .get_single_atom()
        .unwrap()
        .replace_multiple(graph.normal_emr_replacement(
            &plan.support,
            &basis,
            &[W_.x___],
            HedgePair::is_paired,
        ))
        .replace_multiple(&physical)
        .replace(function!(GS.emr_mom, usize::from(k), W_.x___))
        .with(-b / (Atom::num(2) * a) * function!(GS.emr_mom, usize::from(p), W_.x___));
    let expected = canonical_tensor(&numerator)? / Atom::var(GS.dim_epsilon);
    assert!(!expected.is_zero());
    let actual = plan.pole_part(&graph, &model, &vakint, &settings)?;
    let mut retains_other_external = false;
    actual.visitor(&mut |value| {
        if let AtomView::Fun(momentum) = value
            && momentum.get_symbol() == GS.emr_mom
            && momentum.get_nargs() > 0
            && usize::try_from(momentum.get(0)).ok() == Some(usize::from(external[1].0))
        {
            retains_other_external = true;
        }
        !retains_other_external
    });
    assert!(
        !retains_other_external,
        "bare pole has unrouted external momentum: {actual}"
    );
    assert_same_counterterm(&actual, &expected)
}
