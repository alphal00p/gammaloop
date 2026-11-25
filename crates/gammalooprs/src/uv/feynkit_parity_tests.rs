//! Parity between woods built from FeynKit diagrams and from their GammaLoop graphs.
//!
//! The runtime conversion keeps every half-edge and edge index, so subgraph
//! selections of a diagram and of its GammaLoop graph are directly comparable.
//! Both sides classify the diagram's internal propagators, FeynKit's UV domain;
//! GammaLoop's production domain selection (cuts and initial states) is not
//! exercised here.

use std::collections::BTreeSet;

use color_eyre::Result;
use feynkit_generator::FeynmanDiagramRulesExt;
use feynkit_graph::{
    EdgeId, FeynmanDiagram, LoopMomentumBasis as DiagramBasis, NumeratorInterface, VertexId,
    momentum_symbol,
};
use feynkit_model::Model;
use feynkit_uv::{
    ApproximationType, CTRenormalizationRule, RenormalizationPrescriptionSettings, Spinney,
    UvGraph, Wood,
};
use linnet::half_edge::involution::EdgeIndex;
use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};
use linnet::num_traits::SignOrZero;
use symbolica::atom::{Atom, AtomCore};
use symbolica::{function, symbol};

use crate::feyngen::feynkit::FeynmanDiagramGammaLoopExt;
use crate::graph::LoopMomentumBasis;
use crate::initialisation::test_initialise;

/// A sub-basis as its ordered loop edges and, for each internal edge of the
/// spinney, the signs of its loop coordinates.
type BasisSummary = (Vec<usize>, Vec<(usize, Vec<i8>)>);
/// A spinney as its half-edges, degree of divergence, scheme, component count
/// and sub-basis.
type SpinneySummary = (SuBitGraph, i32, ApproximationType, usize, BasisSummary);
/// A forest as its Foata levels of operation scopes.
type ForestSummary = Vec<Vec<SuBitGraph>>;

/// The classified spinneys, wood inclusions and forests of one graph, each as a
/// sorted list so that duplicates are compared too.
#[derive(Debug, PartialEq)]
struct Combinatorics {
    spinneys: Vec<SpinneySummary>,
    inclusions: Vec<(SuBitGraph, SuBitGraph, SuBitGraph)>,
    forests: Vec<ForestSummary>,
}

fn sign(sign: SignOrZero) -> i8 {
    match sign {
        SignOrZero::Plus => 1,
        SignOrZero::Minus => -1,
        SignOrZero::Zero => 0,
    }
}

fn gammaloop_basis(basis: &LoopMomentumBasis, edges: &[EdgeIndex]) -> BasisSummary {
    (
        basis.loop_edges.iter().map(|edge| edge.0).collect(),
        edges
            .iter()
            .map(|edge| {
                let loops = basis.edge_signatures[*edge]
                    .internal
                    .iter()
                    .map(|coefficient| sign(*coefficient))
                    .collect();
                (edge.0, loops)
            })
            .collect(),
    )
}

fn diagram_basis(basis: &DiagramBasis, edges: &[EdgeIndex]) -> BasisSummary {
    (
        basis.loop_edges.iter().map(|edge| edge.0).collect(),
        edges
            .iter()
            .map(|edge| {
                let loops = basis.edge_signatures[&EdgeId(edge.0)]
                    .loops
                    .iter()
                    .map(sign)
                    .collect();
                (edge.0, loops)
            })
            .collect(),
    )
}

impl Combinatorics {
    fn of<G: UvGraph>(
        graph: &G,
        domain: &SuBitGraph,
        model: &Model,
        prescription: &RenormalizationPrescriptionSettings,
        basis: &G::Basis,
        summarize_basis: impl Fn(&G::Basis, &[EdgeIndex]) -> BasisSummary,
    ) -> Result<Self>
    where
        G::Basis: PartialEq,
    {
        let spinneys = graph.classify_spinneys(domain, model, prescription, basis)?;
        let mut summaries = spinneys
            .iter()
            .filter(|spinney| !spinney.is_empty())
            .map(|spinney| Self::spinney_summary(graph, spinney, &summarize_basis))
            .collect::<Vec<_>>();
        summaries.sort();
        let wood = Wood::from_spinneys(spinneys, graph, basis)?;
        let mut inclusions = wood
            .inclusions()
            .map(|inclusion| {
                (
                    wood.spinney(inclusion.smaller).filter().clone(),
                    wood.spinney(inclusion.larger).filter().clone(),
                    inclusion.scope.clone(),
                )
            })
            .collect::<Vec<_>>();
        inclusions.sort();
        let mut forests = wood
            .unfold()
            .operations()
            .map(|(_, forest)| {
                forest
                    .levels()
                    .into_iter()
                    .map(|level| level.into_iter().cloned().collect())
                    .collect()
            })
            .collect::<Vec<_>>();
        forests.sort();
        Ok(Self {
            spinneys: summaries,
            inclusions,
            forests,
        })
    }

    fn spinney_summary<G: UvGraph>(
        graph: &G,
        spinney: &Spinney<G::Basis>,
        summarize_basis: impl Fn(&G::Basis, &[EdgeIndex]) -> BasisSummary,
    ) -> SpinneySummary {
        (
            spinney.filter().clone(),
            spinney.dod(),
            spinney.renormalization_scheme(),
            spinney.n_components(),
            summarize_basis(spinney.lmb(), &spinney.edges(graph)),
        )
    }
}

/// Cycle unions on which the two power countings disagree, as
/// `(subgraph label, FeynKit local, GammaLoop local, FeynKit, GammaLoop)`.
type DodDisagreement = (String, i32, i32, i32, i32);

struct Parity {
    feynkit: Combinatorics,
    gammaloop: Combinatorics,
    disagreements: Vec<DodDisagreement>,
}

/// Build the FeynKit and GammaLoop combinatorics of one diagram, together with
/// every cycle union whose degree of divergence differs between the two.
fn parity(
    diagram: &FeynmanDiagram,
    prescription: &RenormalizationPrescriptionSettings,
) -> Result<Parity> {
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    assert_eq!(
        diagram.hedge_graph().n_hedges(),
        graph.hedge_graph().n_hedges()
    );
    let domain = diagram.internal_subgraph();

    let mut disagreements = Vec::new();
    for cycle in diagram.all_cycle_unions(&domain) {
        let counts = (
            diagram.local_degree_of_divergence(&cycle.filter)?,
            graph.local_degree_of_divergence(&cycle.filter)?,
            diagram.degree_of_divergence(&cycle.filter)?,
            graph.degree_of_divergence(&cycle.filter)?,
        );
        if counts.0 != counts.1 || counts.2 != counts.3 {
            disagreements.push((cycle.string_label(), counts.0, counts.1, counts.2, counts.3));
        }
    }
    disagreements.sort();

    let feynkit = Combinatorics::of(
        diagram,
        &domain,
        diagram.model(),
        prescription,
        diagram.loop_momentum_basis(),
        diagram_basis,
    )?;
    let gammaloop = Combinatorics::of(
        &graph,
        &domain,
        &graph.model,
        prescription,
        &graph.loop_momentum_basis,
        gammaloop_basis,
    )?;
    Ok(Parity {
        feynkit,
        gammaloop,
        disagreements,
    })
}

fn assert_parity_with(
    diagram: &FeynmanDiagram,
    prescription: &RenormalizationPrescriptionSettings,
) -> Result<()> {
    let parity = parity(diagram, prescription)?;
    assert!(
        parity.disagreements.is_empty(),
        "power counting disagrees on {:#?}",
        parity.disagreements
    );
    assert_eq!(parity.feynkit, parity.gammaloop);
    Ok(())
}

fn ruled(model: Model, dot: &str) -> Result<FeynmanDiagram> {
    test_initialise()?;
    Ok(FeynmanDiagram::from_dot(model, dot)?.apply_feynman_rules()?)
}

fn assert_parity(model: Model, dot: &str) -> Result<()> {
    let diagram = ruled(model, dot)?;
    let prescription = RenormalizationPrescriptionSettings::default();
    assert_parity_with(&diagram, &prescription)?;
    assert!(!parity(&diagram, &prescription)?.feynkit.spinneys.is_empty());
    Ok(())
}

/// The energy component `momentum(edge, cind(0))` of an edge momentum.
fn energy(edge: EdgeId) -> Atom {
    function!(
        momentum_symbol(),
        edge.0,
        function!(symbol!("spenso::cind"), 0)
    )
}

fn first_edge(diagram: &FeynmanDiagram, external: bool) -> EdgeId {
    diagram
        .edges()
        .find(|(_, _, edge)| edge.external.is_some() == external)
        .expect("the fixture has such an edge")
        .0
}

fn vertex(diagram: &FeynmanDiagram, name: &str) -> VertexId {
    diagram
        .vertices()
        .find(|(_, vertex)| vertex.name == name)
        .expect("the fixture names this vertex")
        .0
}

const QUARTIC_BUBBLE: &str = r#"digraph bubble {
    edge [particle=phi];
    ext [style=invis];
    ext -> a; ext -> a; a -> b; a -> b; b -> ext; b -> ext;
}"#;

#[test]
fn scalar_basketball() -> Result<()> {
    assert_parity(
        Model::phi4(),
        r#"digraph basketball {
            edge [particle=phi];
            a -> b; a -> b; a -> b; a -> b;
        }"#,
    )
}

#[test]
fn scalar_nested_vertex() -> Result<()> {
    assert_parity(
        Model::phi4(),
        r#"digraph eye {
            edge [particle=phi];
            ext [style=invis];
            ext -> a; ext -> a; a -> b; a -> c; b -> c; b -> c; b -> ext; c -> ext;
        }"#,
    )
}

#[test]
fn scalar_disjoint_tadpoles() -> Result<()> {
    assert_parity(
        Model::phi3(),
        r#"digraph spectacles {
            edge [particle=phi];
            a -> a; a -> b; b -> b;
        }"#,
    )
}

#[test]
fn scalar_overlapping_sunrise() -> Result<()> {
    assert_parity(
        Model::phi4(),
        r#"digraph sunrise {
            edge [particle=phi];
            ext [style=invis];
            ext -> a; a -> b; a -> b; a -> b; b -> ext;
        }"#,
    )
}

#[test]
fn scalar_sewn_initial_state() -> Result<()> {
    // A forward-scattering cut sunrise: the sewn initial state is a paired
    // external edge, which neither side treats as a UV propagator.
    assert_parity(
        Model::phi4(),
        r#"digraph sewn_sunrise {
            graph [final_state="phi,phi,phi"];
            edge [particle=phi];
            ext [style=invis];
            ext -> a [is_cut=9]; b -> ext [is_cut=9];
            a -> b; a -> b; a -> b;
        }"#,
    )
}

#[test]
fn qcd_nested_quark_self_energy() -> Result<()> {
    assert_parity(
        Model::qcd(),
        r#"digraph rainbow {
            ext [style=invis];
            ext -> a [particle="d"];
            a -> b [particle="d"]; b -> c [particle="d"]; c -> d [particle="d"];
            d -> ext [particle="d"];
            a -> d [particle="g"]; b -> c [particle="g"];
        }"#,
    )
}

#[test]
fn qcd_three_gluon_vertices_with_external_momentum() -> Result<()> {
    assert_parity(
        Model::qcd(),
        r#"digraph gluon_bubble {
            edge [particle="g"];
            ext [style=invis];
            ext -> a; a -> b; a -> b; b -> ext;
        }"#,
    )
}

#[test]
fn qcd_fermion_loop() -> Result<()> {
    assert_parity(
        Model::qcd(),
        r#"digraph quark_loop {
            ext [style=invis];
            ext -> a [particle="g"];
            a -> b [particle="d"]; b -> a [particle="d"];
            b -> ext [particle="g"];
        }"#,
    )
}

#[test]
fn scheme_buckets_and_overrides() -> Result<()> {
    let rainbow = ruled(
        Model::qcd(),
        r#"digraph rainbow {
            ext [style=invis];
            ext -> a [particle="d"];
            a -> b [particle="d"]; b -> c [particle="d"]; c -> d [particle="d"];
            d -> ext [particle="d"];
            a -> d [particle="g"]; b -> c [particle="g"];
        }"#,
    )?;
    let domain = rainbow.internal_subgraph();
    let full = rainbow
        .classify_spinneys(
            &domain,
            rainbow.model(),
            &RenormalizationPrescriptionSettings::default(),
            rainbow.loop_momentum_basis(),
        )?
        .into_iter()
        .max_by_key(|spinney| spinney.loop_count())
        .expect("the rainbow has spinneys");
    // Both self-energies share the full graph's particle content, so the
    // override selects both; without it the massive-boundary bucket applies.
    let bucket = RenormalizationPrescriptionSettings {
        massive_power_divergent: ApproximationType::IR,
        ..RenormalizationPrescriptionSettings::default()
    };
    let overridden = RenormalizationPrescriptionSettings {
        overrides: vec![CTRenormalizationRule::new(
            rainbow.ct_identifier(full.filter(), rainbow.model())?,
            ApproximationType::PolePart,
        )],
        ..bucket.clone()
    };
    for (prescription, expected) in [
        (bucket, ApproximationType::IR),
        (overridden, ApproximationType::PolePart),
    ] {
        let parity = parity(&rainbow, &prescription)?;
        assert!(parity.disagreements.is_empty());
        assert_eq!(parity.feynkit, parity.gammaloop);
        let schemes = parity
            .feynkit
            .spinneys
            .iter()
            .map(|(_, _, scheme, _, _)| *scheme)
            .collect::<BTreeSet<_>>();
        assert_eq!(schemes, BTreeSet::from([expected]));
    }
    Ok(())
}

#[test]
fn external_energy_vertex_factor() -> Result<()> {
    // The reviewer's case: a factor of an external energy squared does not raise
    // the classifying degree on either side.
    let bubble = ruled(Model::phi4(), QUARTIC_BUBBLE)?;
    let factor = NumeratorInterface::admit(energy(first_edge(&bubble, true)).pow(2), None)?;
    let a = vertex(&bubble, "a");
    let bubble = bubble.with_numerator_fragments([(a, factor)], [] as [(EdgeId, _); 0])?;
    let prescription = RenormalizationPrescriptionSettings {
        massive_power_divergent: ApproximationType::Unsubtracted,
        ..RenormalizationPrescriptionSettings::default()
    };
    assert_parity_with(&bubble, &prescription)?;
    assert_eq!(bubble.degree_of_divergence(&bubble.internal_subgraph())?, 0);
    Ok(())
}

#[test]
fn derivative_coupling_on_an_external_leg() -> Result<()> {
    // A vertex factor linear in an external and in a loop momentum.
    let bubble = ruled(Model::phi4(), QUARTIC_BUBBLE)?;
    let factor = NumeratorInterface::admit(
        energy(first_edge(&bubble, true)) * energy(first_edge(&bubble, false)),
        None,
    )?;
    let a = vertex(&bubble, "a");
    let bubble = bubble.with_numerator_fragments([(a, factor)], [] as [(EdgeId, _); 0])?;
    assert_parity_with(&bubble, &RenormalizationPrescriptionSettings::default())?;
    assert_eq!(bubble.degree_of_divergence(&bubble.internal_subgraph())?, 1);
    Ok(())
}

/// An edge numerator factor that is an inverse quadratic in the edge's own
/// energy. Both hosts count it as a numerator factor next to one bare massless
/// quadratic per line; this does not certify propagator powers, for which
/// FeynKit has no diagram-level representation.
#[test]
fn inverse_quadratic_edge_factor() -> Result<()> {
    let bubble = ruled(Model::phi4(), QUARTIC_BUBBLE)?;
    let edge = first_edge(&bubble, false);
    let factor = NumeratorInterface::admit((energy(edge).pow(2) + Atom::num(1)).pow(-1), None)?;
    let bubble = bubble.with_numerator_fragments([] as [(VertexId, _); 0], [(edge, factor)])?;
    assert_parity_with(&bubble, &RenormalizationPrescriptionSettings::default())?;
    assert_eq!(
        bubble.degree_of_divergence(&bubble.internal_subgraph())?,
        -2
    );
    Ok(())
}

/// Expected difference, kept on purpose: two parallel bubble lines oriented
/// alike with the vertex factor `(q₁⁰ + q₂⁰)²`. Momentum conservation routes the
/// lines as `k` and `p − k`, so FeynKit sees the external energy `p⁰` alone and
/// classifies the bubble as logarithmic;
/// GammaLoop tags each line's hard momentum separately and counts the factor
/// with degree two. GammaLoop's convention is the reference for counterterms.
#[test]
fn expected_difference_tagged_hard_momenta() -> Result<()> {
    let bubble = ruled(Model::phi4(), QUARTIC_BUBBLE)?;
    let lines = bubble
        .edges()
        .filter(|(_, _, edge)| edge.external.is_none())
        .map(|(edge, _, _)| edge)
        .collect::<Vec<_>>();
    let [first, second] = lines.as_slice() else {
        panic!("the bubble has two lines");
    };
    let factor = NumeratorInterface::admit((energy(*first) + energy(*second)).pow(2), None)?;
    let a = vertex(&bubble, "a");
    let bubble = bubble.with_numerator_fragments([(a, factor)], [] as [(EdgeId, _); 0])?;
    let graph = bubble.to_gamma_loop_graph(None, true)?;
    let loops = bubble.internal_subgraph();
    assert_eq!(bubble.degree_of_divergence(&loops)?, 0);
    assert_eq!(graph.degree_of_divergence(&loops)?, 2);
    assert_eq!(
        bubble.local_degree_of_divergence(&loops)?,
        graph.local_degree_of_divergence(&loops)?
    );

    let power_divergences_unsubtracted = RenormalizationPrescriptionSettings {
        massive_power_divergent: ApproximationType::Unsubtracted,
        massless_power_divergent: ApproximationType::Unsubtracted,
        ..RenormalizationPrescriptionSettings::default()
    };
    let parity = parity(&bubble, &power_divergences_unsubtracted)?;
    assert_eq!(parity.feynkit.spinneys.len(), 1);
    assert!(parity.gammaloop.spinneys.is_empty());
    Ok(())
}

/// A local numerator that vanishes identically has degree zero on both hosts
/// and its spinney is kept. Recorded as current behaviour, not as a guarantee.
#[test]
fn vanishing_numerator_is_kept_with_degree_zero() -> Result<()> {
    let bubble = ruled(Model::phi4(), QUARTIC_BUBBLE)?;
    let a = vertex(&bubble, "a");
    let zero = NumeratorInterface::admit(Atom::Zero, None)?;
    let bubble = bubble.with_numerator_fragments([(a, zero)], [] as [(EdgeId, _); 0])?;
    let prescription = RenormalizationPrescriptionSettings::default();
    assert_parity_with(&bubble, &prescription)?;
    let loops = bubble.internal_subgraph();
    assert_eq!(bubble.degree_of_divergence(&loops)?, 0);
    assert_eq!(parity(&bubble, &prescription)?.feynkit.spinneys.len(), 1);
    Ok(())
}

/// Pruning and classifying degrees of both hosts on a scalar triangle with a
/// vertex factor that vanishes after routing, a literal zero factor, and a
/// squared energy of another line as a factor on edge 0. FeynKit prunes with at
/// least its classifying degree and keeps all three; GammaLoop's values are its
/// current behaviour, recorded, not changed: its stored-degree pruning is below
/// its classifier in the last two cases.
#[test]
fn pruning_and_classifying_degrees_on_vanishing_and_foreign_factors() -> Result<()> {
    test_initialise()?;
    let triangle = FeynmanDiagram::from_dot(
        Model::phi3(),
        r#"digraph triangle {
            ext [style=invis]; edge [particle=phi];
            a -> b [id=0]; a -> c [id=1]; b -> c [id=2];
            ext -> a [id=3]; b -> ext [id=4]; c -> ext [id=5];
        }"#,
    )?;
    let a = vertex(&triangle, "a");
    let vertex_factor = |factor| -> Result<FeynmanDiagram> {
        let fragment = NumeratorInterface::admit(factor, None)?;
        Ok(triangle
            .clone()
            .with_numerator_fragments([(a, fragment)], [] as [(EdgeId, _); 0])?)
    };
    let edge_factor = NumeratorInterface::admit(energy(EdgeId(1)).pow(2), None)?;
    let cases = [
        (
            "vanishes after routing",
            vertex_factor(energy(EdgeId(0)) + energy(EdgeId(1)) - energy(EdgeId(3)))?,
            (0, 0),
            (-1, -1),
        ),
        ("literal zero", vertex_factor(Atom::Zero)?, (0, 0), (-2, 0)),
        (
            "another line's energy on edge 0",
            triangle
                .clone()
                .with_numerator_fragments([], [(EdgeId(0), edge_factor)])?,
            (0, 0),
            (-2, 0),
        ),
    ];
    for (label, weighted, feynkit, gammaloop) in cases {
        let graph = weighted.to_gamma_loop_graph(None, true)?;
        let loops = weighted.internal_subgraph();
        assert_eq!(
            (
                weighted.local_degree_of_divergence(&loops)?,
                weighted.degree_of_divergence(&loops)?,
            ),
            feynkit,
            "{label}"
        );
        assert_eq!(
            (
                graph.local_degree_of_divergence(&loops)?,
                graph.degree_of_divergence(&loops)?,
            ),
            gammaloop,
            "{label}"
        );
        assert!(
            weighted
                .spinneys(&loops)?
                .iter()
                .any(|union| union.filter == loops),
            "{label}"
        );
    }
    Ok(())
}
