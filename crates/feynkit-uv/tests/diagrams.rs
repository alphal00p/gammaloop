//! Woods and forests built directly from FeynKit diagrams.

use std::collections::BTreeSet;
use std::error::Error;

use feynkit_generator::FeynmanDiagramRulesExt;
use feynkit_graph::{
    DiagramEdge, DiagramError, DiagramVertex, EdgeId, FeynmanDiagram, LoopMomentumBasis,
    NumeratorInterface, momentum_symbol,
};
use feynkit_model::Model;
use feynkit_uv::{
    ApproximationType, Counterterm, CountertermError, CountertermScheme, ForeignReason,
    ForestError, Forests, OperationSpinneyError, RenormalizationPrescriptionSettings, Spinney,
    SpinneyRelation, UvGraph, UvGraphError, Wood,
};
use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use linnet::half_edge::subgraph::{InternalSubGraph, SuBitGraph, SubSetLike};
use symbolica::atom::{Atom, AtomCore};
use symbolica::{function, symbol};

fn diagram(model: Model, dot: &str) -> FeynmanDiagram {
    FeynmanDiagram::from_dot(model, dot)
        .unwrap()
        .apply_feynman_rules()
        .unwrap()
}

fn spinneys(diagram: &FeynmanDiagram) -> Vec<Spinney<LoopMomentumBasis>> {
    diagram
        .classify_spinneys(
            &diagram.hedge_graph().full_filter(),
            diagram.model(),
            &RenormalizationPrescriptionSettings::default(),
            diagram.loop_momentum_basis(),
        )
        .unwrap()
}

/// Each nonempty spinney as its sorted internal edge ids, with its loop count and DOD.
/// Classification also returns the empty subgraph, which the wood merges with its root.
fn summary(diagram: &FeynmanDiagram) -> BTreeSet<(Vec<usize>, usize, i32)> {
    spinneys(diagram)
        .iter()
        .filter(|spinney| !spinney.is_empty())
        .map(|spinney| {
            let edges = spinney
                .edges(diagram)
                .into_iter()
                .map(|EdgeIndex(edge)| edge)
                .collect();
            (edges, spinney.loop_count(), spinney.dod())
        })
        .collect()
}

fn forests(diagram: &FeynmanDiagram) -> Forests<LoopMomentumBasis> {
    Wood::from_spinneys(spinneys(diagram), diagram, diagram.loop_momentum_basis())
        .unwrap()
        .unfold()
}

fn edges(diagram: &FeynmanDiagram, filter: &SuBitGraph) -> Vec<usize> {
    diagram
        .hedge_graph()
        .iter_edges_of(filter)
        .filter_map(|(pair, EdgeIndex(edge), _)| pair.is_paired().then_some(edge))
        .sorted()
        .collect()
}

#[test]
fn scalar_basketball_has_eleven_spinneys() {
    let basketball = diagram(
        Model::phi4(),
        r#"digraph basketball {
            a -> b [id=0, particle="phi"]; a -> b [id=1, particle="phi"];
            a -> b [id=2, particle="phi"]; a -> b [id=3, particle="phi"];
        }"#,
    );
    let spinneys = summary(&basketball);
    let by_loops = |loops: usize| spinneys.iter().filter(|(_, l, _)| *l == loops).count();
    assert_eq!(spinneys.len(), 11);
    assert_eq!((by_loops(1), by_loops(2), by_loops(3)), (6, 4, 1));
    // 4L - 2I in four dimensions for a quartic scalar.
    assert!(
        spinneys
            .iter()
            .all(|(edges, loops, dod)| { *dod == 4 * *loops as i32 - 2 * edges.len() as i32 })
    );

    let wood = Wood::from_spinneys(
        self::spinneys(&basketball),
        &basketball,
        basketball.loop_momentum_basis(),
    )
    .unwrap();
    assert_eq!(wood.n_spinneys(), 12);
    // Hasse diagram: root ⊂ 6 pairs, each pair ⊂ 2 triples, each triple ⊂ the full graph.
    assert_eq!(wood.cover_relation().len(), 6 + 12 + 4);
    assert_eq!(wood.max_loop_count(), 3);
    // Every inclusion goes from a smaller to a strictly larger spinney.
    for inclusion in wood.inclusions() {
        let smaller = wood.spinney(inclusion.smaller);
        let larger = wood.spinney(inclusion.larger);
        assert!(smaller.loop_count() < larger.loop_count());
        assert!(!inclusion.scope.is_empty());
    }

    let forests = wood.unfold();
    // Agrees with GammaLoop's hedge-poset basketball regression.
    assert_eq!(forests.n_forests(), 46);
    let root = forests.operation(forests.root());
    assert_eq!((root.depth(), root.sign(), root.covers()), (0, 1, None));
    for (node, forest) in forests.operations() {
        assert_eq!(
            forest.levels().iter().map(Vec::len).sum::<usize>(),
            forest.depth()
        );
        assert_eq!(forest.sign(), if forest.depth() % 2 == 0 { 1 } else { -1 });
        if node != forests.root() {
            assert_eq!(
                forest.covers().as_ref(),
                Some(forests.source_spinney(node).filter())
            );
        }
    }
    // The deepest forests are full chains: pair ⊂ triple ⊂ full graph.
    let max_depth = forests.operations().map(|(_, forest)| forest.depth()).max();
    assert_eq!(max_depth, Some(3));
    for edge in forests.edges() {
        assert_eq!(
            forests.operation(edge.child).depth(),
            forests.operation(edge.parent).depth() + 1
        );
    }
}

#[test]
fn nested_vertex_subgraph() {
    // Quartic four-point "eye": the b-c bubble sits inside the full graph.
    let eye = diagram(
        Model::phi4(),
        r#"digraph eye {
            ext [style=invis];
            ext -> a [id=0, particle="phi"]; ext -> a [id=1, particle="phi"];
            a -> b [id=2, particle="phi"]; a -> c [id=3, particle="phi"];
            b -> c [id=4, particle="phi"]; b -> c [id=5, particle="phi"];
            b -> ext [id=6, particle="phi"]; c -> ext [id=7, particle="phi"];
        }"#,
    );
    let spinneys = summary(&eye);
    assert_eq!(
        spinneys,
        BTreeSet::from([(vec![4, 5], 1, 0), (vec![2, 3, 4, 5], 2, 0)])
    );

    let forests = forests(&eye);
    // ∅, {inner}, {full}, {inner}·{full}.
    assert_eq!(forests.n_forests(), 4);
    let chains = forests
        .operations()
        .map(|(_, forest)| {
            forest
                .levels()
                .into_iter()
                .map(|level| level.into_iter().map(|scope| edges(&eye, scope)).collect())
                .collect::<Vec<Vec<_>>>()
        })
        .collect::<BTreeSet<_>>();
    assert!(chains.contains(&vec![vec![vec![4, 5]], vec![vec![2, 3, 4, 5]]]));
}

#[test]
fn disjoint_tadpoles_form_a_factorized_union() {
    // Cubic "spectacles": two tadpoles joined by a bridge.
    let spectacles = diagram(
        Model::phi3(),
        r#"digraph spectacles {
            a -> a [id=0, particle="phi"]; a -> b [id=1, particle="phi"]; b -> b [id=2, particle="phi"];
        }"#,
    );
    let spinneys = self::spinneys(&spectacles);
    assert_eq!(summary(&spectacles).len(), 3);
    let union = spinneys
        .iter()
        .find(|spinney| spinney.n_components() == 2)
        .expect("the two tadpoles form a disconnected union");
    assert_eq!(union.renormalization_scheme(), ApproximationType::MUV);
    assert_eq!(union.loop_count(), 2);
    assert_eq!(union.edges(&spectacles), vec![EdgeIndex(0), EdgeIndex(2)]);

    let forests = forests(&spectacles);
    // ∅, {t1}, {t2} and the commuting product {t1,t2}.
    assert_eq!(forests.n_forests(), 4);
    let union_nodes = forests
        .operations()
        .filter(|(node, _)| forests.is_disjoint_union(*node))
        .map(|(node, forest)| (node, forest.levels().len(), forest.depth()))
        .collect::<Vec<_>>();
    assert_eq!(union_nodes.len(), 1);
    let (union_node, levels, depth) = union_nodes[0];
    assert_eq!((levels, depth), (1, 2));
    assert_eq!(
        forests
            .disconnected_component_nodes(union_node)
            .unwrap()
            .len(),
        2
    );
}

#[test]
fn overlapping_sunrise_subgraphs() {
    let sunrise = diagram(
        Model::phi4(),
        r#"digraph sunrise {
            ext [style=invis];
            ext -> a [id=0, particle="phi"];
            a -> b [id=1, particle="phi"]; a -> b [id=2, particle="phi"]; a -> b [id=3, particle="phi"];
            b -> ext [id=4, particle="phi"];
        }"#,
    );
    let spinneys = summary(&sunrise);
    assert_eq!(
        spinneys,
        BTreeSet::from([
            (vec![1, 2], 1, 0),
            (vec![1, 3], 1, 0),
            (vec![2, 3], 1, 0),
            (vec![1, 2, 3], 2, 2),
        ])
    );
    let forests = forests(&sunrise);
    // The three bubbles overlap pairwise: no forest holds two of them.
    assert!(forests.operations().all(|(_, forest)| forest.depth() <= 2));
    assert_eq!(forests.n_forests(), 8);
}

#[test]
fn qcd_nested_quark_self_energy() {
    // Two-loop rainbow: the b-c quark self-energy is nested in the full graph.
    let rainbow = diagram(
        Model::qcd(),
        r#"digraph rainbow {
            ext [style=invis];
            ext -> a [id=0, particle="d"];
            a -> b [id=1, particle="d"]; b -> c [id=2, particle="d"]; c -> d [id=3, particle="d"];
            d -> ext [id=4, particle="d"];
            a -> d [id=5, particle="g"]; b -> c [id=6, particle="g"];
        }"#,
    );
    let spinneys = summary(&rainbow);
    assert_eq!(
        spinneys,
        BTreeSet::from([(vec![2, 6], 1, 1), (vec![1, 2, 3, 5, 6], 2, 1)])
    );
    let forests = forests(&rainbow);
    assert_eq!(forests.n_forests(), 4);

    let mut dot = String::new();
    forests
        .dot_serialize_with(&mut dot, |node| forests.operation(node).to_string())
        .unwrap();
    for attribute in ["label=", "foata=", "cover="] {
        assert!(dot.contains(attribute), "{dot}");
    }
}

#[test]
fn dummy_self_loops_never_become_spinneys() -> Result<(), Box<dyn Error>> {
    let dummy_loop = self_loops(true)?;

    for domain in [
        dummy_loop.hedge_graph().full_filter(),
        dummy_loop.internal_subgraph(),
    ] {
        let spinneys = dummy_loop.classify_spinneys(
            &domain,
            dummy_loop.model(),
            &RenormalizationPrescriptionSettings::default(),
            dummy_loop.loop_momentum_basis(),
        )?;
        let nonempty = spinneys
            .iter()
            .filter(|spinney| !spinney.is_empty())
            .collect::<Vec<_>>();
        assert_eq!(nonempty.len(), 1);
        let tadpole = nonempty[0];
        assert_eq!(tadpole.edges(&dummy_loop), vec![EdgeIndex(1)]);
        assert_eq!(tadpole.loop_count(), tadpole.lmb().loop_edges.len());
        assert_eq!(tadpole.loop_count(), 1);
    }
    Ok(())
}

#[test]
fn particle_queries_reject_a_foreign_model() {
    let bubble = diagram(
        Model::phi4(),
        r#"digraph bubble {
            edge [particle=phi];
            ext [style=invis];
            ext -> a; ext -> a; a -> b; a -> b; b -> ext; b -> ext;
        }"#,
    );
    let edge = EdgeIndex(0);
    assert_eq!(
        bubble.edge_pdg_code(edge, bubble.model()).unwrap(),
        Some(9_000_001)
    );
    let foreign = Model::qcd();
    assert!(matches!(
        bubble.edge_pdg_code(edge, &foreign),
        Err(DiagramError::ModelFingerprintMismatch { .. })
    ));
    assert!(matches!(
        bubble.is_massive_edge(edge, &foreign),
        Err(DiagramError::ModelFingerprintMismatch { .. })
    ));
    assert!(
        bubble
            .classify_spinneys(
                &bubble.hedge_graph().full_filter(),
                &foreign,
                &RenormalizationPrescriptionSettings::default(),
                bubble.loop_momentum_basis(),
            )
            .is_err()
    );
}

/// A one-loop scalar bubble whose vertex `a` carries the squared energy of an
/// external leg: a factor that does not grow with the loop momentum.
fn externally_weighted_bubble() -> Result<FeynmanDiagram, Box<dyn Error>> {
    let bubble = FeynmanDiagram::from_dot(
        Model::phi4(),
        r#"digraph bubble {
            ext [style=invis];
            edge [particle=phi];
            ext -> a; ext -> a; a -> b; a -> b; b -> ext; b -> ext;
        }"#,
    )?;
    let external = bubble
        .edges()
        .find(|(_, _, edge)| edge.external.is_some())
        .expect("the bubble has external legs")
        .0;
    let vertex = bubble
        .vertices()
        .find(|(_, vertex)| vertex.name == "a")
        .expect("the bubble has vertex a")
        .0;
    let energy = function!(
        momentum_symbol(),
        external.0,
        function!(symbol!("spenso::cind"), 0)
    );
    let fragment = NumeratorInterface::admit(energy.pow(2), None)?;
    Ok(bubble.with_numerator_fragments([(vertex, fragment)], [] as [(EdgeId, _); 0])?)
}

#[test]
fn external_numerator_factors_do_not_raise_the_classifying_degree() -> Result<(), Box<dyn Error>> {
    let bubble = externally_weighted_bubble()?;
    let loops = bubble.internal_subgraph();
    // The pruning degree scales every vertex momentum; the classifying degree
    // rescales only the loop momentum, as GammaLoop's runtime power counting does.
    assert_eq!(bubble.local_degree_of_divergence(&loops)?, 2);
    assert_eq!(bubble.degree_of_divergence(&loops)?, 0);

    let prescription = RenormalizationPrescriptionSettings {
        massive_power_divergent: ApproximationType::Unsubtracted,
        ..RenormalizationPrescriptionSettings::default()
    };
    let spinneys = bubble.classify_spinneys(
        &loops,
        bubble.model(),
        &prescription,
        bubble.loop_momentum_basis(),
    )?;
    let nonempty = spinneys
        .iter()
        .filter(|spinney| !spinney.is_empty())
        .collect::<Vec<_>>();
    assert_eq!(nonempty.len(), 1);
    assert_eq!(nonempty[0].renormalization_scheme(), ApproximationType::MUV);
    Ok(())
}

/// Structured forest DOT with FeynKit's own Typst-string labels for a tadpole
/// joined by a bridge to a bubble carrying a second tadpole: nested chains
/// (`CO;GG`) and factorised joins (`CO,3`).
const NESTED_FACTORISED_FOREST_DOT: &str = r##"
digraph {
  node	 [shape=circle,height=0.1,label=""];
  overlap = "scale";
  layout = "neato";
  start=2;
  0	 [label="\"∅\"" foata="" cover="0"];
  1	 [label="\"{CO}\"" foata="CO" cover="CO"];
  2	 [label="\"{3s}\"" foata="3s" cover="3s"];
  3	 [label="\"{GG}\"" foata="GG" cover="GG"];
  4	 [label="\"{3}\"" foata="3" cover="3"];
  5	 [label="\"{CO} · {GG}\"" foata="CO;GG" cover="GG"];
  6	 [label="\"{CO,3}\"" foata="CO,3" cover="CR"];
  7	 [label="\"{3s} · {GG}\"" foata="3s;GG" cover="GG"];
  8	 [label="\"{3s,3}\"" foata="3s,3" cover="3v"];
  9	 [label="\"{GG,3}\"" foata="GG,3" cover="GJ"];
  10	 [label="\"{CO} · {GG,3}\"" foata="CO;GG,3" cover="GJ"];
  11	 [label="\"{3s} · {GG,3}\"" foata="3s;GG,3" cover="GJ"];
  0:0:s	-> 1:1:s	 [id=0  color="red:blue;0.5"];
  0:2:s	-> 2:3:s	 [id=1  color="red:blue;0.5"];
  0:4:s	-> 3:5:s	 [id=2  color="red:blue;0.5"];
  0:6:s	-> 4:7:s	 [id=3  color="red:blue;0.5"];
  1:8:s	-> 5:9:s	 [id=4  color="red:blue;0.5"];
  1:10:s	-> 6:11:s	 [id=5  color="red:blue;0.5"];
  2:12:s	-> 7:13:s	 [id=6  color="red:blue;0.5"];
  2:14:s	-> 8:15:s	 [id=7  color="red:blue;0.5"];
  3:16:s	-> 9:17:s	 [id=8  color="red:blue;0.5"];
  4:18:s	-> 9:19:s	 [id=9  color="red:blue;0.5"];
  4:20:s	-> 8:21:s	 [id=10  color="red:blue;0.5"];
  4:22:s	-> 6:23:s	 [id=11  color="red:blue;0.5"];
  4:24:s	-> 10:25:s	 [id=12  color="red:blue;0.5"];
  5:26:s	-> 10:27:s	 [id=13  color="red:blue;0.5"];
  4:28:s	-> 11:29:s	 [id=14  color="red:blue;0.5"];
  7:30:s	-> 11:31:s	 [id=15  color="red:blue;0.5"];
}
"##;

#[test]
fn forest_dot_with_feynkit_labels_is_exact() {
    let graph = diagram(
        Model::phi_3_4(),
        r#"digraph nested_factorised {
            edge [particle=phi];
            a -> a; a -> b; b -> c; b -> c; c -> c;
        }"#,
    );
    assert_eq!(
        forests(&graph).dot_serialize(),
        NESTED_FACTORISED_FOREST_DOT
            .strip_prefix('\n')
            .expect("the fixture starts on its own line")
    );
}

/// The Petersen graph: cubic, six loops, without a nowhere-zero 4-flow, so it
/// is not the union of two cycle-space elements.
const PETERSEN_EDGES: &str = "
    v0 -> v1; v1 -> v2; v2 -> v3; v3 -> v4; v4 -> v0;
    v0 -> v5; v1 -> v6; v2 -> v7; v3 -> v8; v4 -> v9;
    v5 -> v7; v7 -> v9; v9 -> v6; v6 -> v8; v8 -> v5;";

/// Known limitation of the pruned enumeration, inherited from GammaLoop's
/// algorithm: a divergent union that is not the union of two cycle-space
/// elements is found only through retained intermediate unions. Weighting six
/// vertices of the scalar Petersen graph by one incident energy each makes the
/// full graph logarithmic while every proper union is convergent, so the full
/// graph is never reached. This test records today's answer and the complete
/// one; it must be revisited if the enumeration changes.
#[test]
fn known_limitation_petersen_full_graph_is_not_reached() -> Result<(), Box<dyn Error>> {
    let petersen = FeynmanDiagram::from_dot(
        Model::phi3(),
        &format!("digraph petersen {{ edge [particle=phi]; {PETERSEN_EDGES} }}"),
    )?;
    let weights = (0..6)
        .map(|index| {
            let vertex = petersen
                .vertices()
                .find(|(_, vertex)| vertex.name == format!("v{index}"))
                .expect("the Petersen graph names its vertices")
                .0;
            let edge = petersen
                .edges()
                .find(|(_, endpoints, _)| {
                    endpoints.source == Some(vertex) || endpoints.target == Some(vertex)
                })
                .expect("every vertex has an incident edge")
                .0;
            let energy = function!(
                momentum_symbol(),
                edge.0,
                function!(symbol!("spenso::cind"), 0)
            );
            Ok((vertex, NumeratorInterface::admit(energy, None)?))
        })
        .collect::<Result<Vec<_>, DiagramError>>()?;
    let petersen = petersen.with_numerator_fragments(weights, [] as [(EdgeId, _); 0])?;
    let domain = petersen.internal_subgraph();

    let unions = petersen.all_cycle_unions(&domain);
    assert_eq!(unions.len(), 594);
    let full = unions
        .iter()
        .find(|union| union.filter == domain)
        .expect("the full graph is a union of cycles");
    assert_eq!(petersen.local_degree_of_divergence(&domain)?, 0);
    assert_eq!(petersen.degree_of_divergence(&domain)?, 0);
    for union in unions
        .iter()
        .filter(|union| !union.is_empty() && union.filter != domain)
    {
        assert!(petersen.local_degree_of_divergence(&union.filter)? < 0);
    }

    // Complete answer: the full graph is the only divergent nonempty union.
    // Today's answer: only the empty subgraph is enumerated.
    let enumerated = petersen.spinneys(&domain)?;
    assert!(!enumerated.contains(full));
    assert_eq!(enumerated.len(), 1);
    Ok(())
}

#[test]
fn physical_petersen_vacua_lose_no_divergent_union() -> Result<(), Box<dyn Error>> {
    let gluons = diagram(
        Model::yang_mills(),
        &format!("digraph petersen {{ edge [particle=\"g\"]; {PETERSEN_EDGES} }}"),
    );
    // Two electron pentagons joined by five photons.
    let photons = diagram(
        Model::qed(),
        r#"digraph petersen {
            v0 -> v1 [particle="e-"]; v1 -> v2 [particle="e-"]; v2 -> v3 [particle="e-"];
            v3 -> v4 [particle="e-"]; v4 -> v0 [particle="e-"];
            v0 -> v5 [particle="a"]; v1 -> v6 [particle="a"]; v2 -> v7 [particle="a"];
            v3 -> v8 [particle="a"]; v4 -> v9 [particle="a"];
            v5 -> v7 [particle="e-"]; v7 -> v9 [particle="e-"]; v9 -> v6 [particle="e-"];
            v6 -> v8 [particle="e-"]; v8 -> v5 [particle="e-"];
        }"#,
    );
    for (vacuum, retained) in [(gluons, 116), (photons, 36)] {
        let domain = vacuum.internal_subgraph();
        assert_eq!(vacuum.degree_of_divergence(&domain)?, 4);
        let mut exhaustive = Vec::new();
        for union in vacuum.all_cycle_unions(&domain) {
            if vacuum.passes_pruning(&union.filter)? {
                exhaustive.push(union);
            }
        }
        assert_eq!(exhaustive.len(), retained);
        let enumerated = vacuum.spinneys(&domain)?;
        assert!(exhaustive.iter().all(|union| enumerated.contains(union)));
    }
    Ok(())
}

/// A quartic vertex with two self-loops; the first is a dummy attachment when
/// `dummy` is set.
fn self_loops(dummy: bool) -> Result<FeynmanDiagram, Box<dyn Error>> {
    let model = Model::phi4();
    let particle = model.particle_id("phi")?;
    let rule = model.vertex_rule_id("V_SCALAR")?;
    let mut builder = FeynmanDiagram::builder(model, "dummy-loop");
    let a = builder.add_vertex(DiagramVertex::interaction("a", rule));
    let mut first = DiagramEdge::new(particle, false);
    first.is_dummy = dummy;
    builder.add_edge(a, a, first)?;
    builder.add_edge(a, a, DiagramEdge::new(particle, false))?;
    Ok(builder.build()?)
}

#[test]
fn dummy_crown_does_not_change_the_scheme() -> Result<(), Box<dyn Error>> {
    let dummy_loop = self_loops(true)?;
    let full = dummy_loop.hedge_graph().full_filter();
    // The dummy edge forms the real tadpole's crown but carries no particle.
    let tadpole = dummy_loop.internal_subgraph();
    assert_eq!(
        dummy_loop.boundary_pdg_set(&tadpole, dummy_loop.model())?,
        BTreeSet::new()
    );
    assert!(!dummy_loop.has_massive_boundary_external(&tadpole, dummy_loop.model())?);
    let prescription = RenormalizationPrescriptionSettings {
        massive_power_divergent: ApproximationType::Unsubtracted,
        ..RenormalizationPrescriptionSettings::default()
    };
    let spinneys = dummy_loop.classify_spinneys(
        &full,
        dummy_loop.model(),
        &prescription,
        dummy_loop.loop_momentum_basis(),
    )?;
    let nonempty = spinneys
        .iter()
        .filter(|spinney| !spinney.is_empty())
        .collect::<Vec<_>>();
    assert_eq!(nonempty.len(), 1);
    assert_eq!(nonempty[0].edges(&dummy_loop), vec![EdgeIndex(1)]);
    assert_eq!(nonempty[0].renormalization_scheme(), ApproximationType::MUV);
    Ok(())
}

#[test]
fn spinney_constructors_reject_selections_outside_the_domain() -> Result<(), Box<dyn Error>> {
    let dummy_loop = self_loops(true)?;
    let graph = dummy_loop.hedge_graph();
    let both_loops = InternalSubGraph::cleaned_filter_optimist(graph.full_filter(), graph);
    assert!(matches!(
        Spinney::new(both_loops, &dummy_loop, dummy_loop.loop_momentum_basis()),
        Err(UvGraphError::OutsideDomain { .. })
    ));
    let tadpole = InternalSubGraph::cleaned_filter_optimist(dummy_loop.internal_subgraph(), graph);
    let spinney = Spinney::new(tadpole, &dummy_loop, dummy_loop.loop_momentum_basis())?
        .expect("the quartic tadpole is quadratically divergent");
    assert_eq!(spinney.loop_count(), spinney.lmb().loop_edges.len());
    Ok(())
}

#[test]
fn single_spinney_entry_points_validate_the_domain_first() -> Result<(), Box<dyn Error>> {
    let dummy_loop = self_loops(true)?;
    let graph = dummy_loop.hedge_graph();
    let both_loops = InternalSubGraph::cleaned_filter_optimist(graph.full_filter(), graph);
    let nothing_subtracted = RenormalizationPrescriptionSettings {
        log_divergent: ApproximationType::Unsubtracted,
        massive_power_divergent: ApproximationType::Unsubtracted,
        massless_power_divergent: ApproximationType::Unsubtracted,
        ..RenormalizationPrescriptionSettings::default()
    };
    assert!(matches!(
        dummy_loop.classify_spinney(
            both_loops.clone(),
            dummy_loop.model(),
            &nothing_subtracted,
            dummy_loop.loop_momentum_basis(),
        ),
        Err(UvGraphError::OutsideDomain { .. })
    ));
    assert!(matches!(
        Spinney::with_scheme(
            both_loops,
            &dummy_loop,
            dummy_loop.loop_momentum_basis(),
            ApproximationType::MUV,
            0,
        ),
        Err(UvGraphError::OutsideDomain { .. })
    ));
    Ok(())
}

#[test]
fn woods_reject_spinneys_of_another_graph() -> Result<(), Box<dyn Error>> {
    // The same half-edges are two real loops in one diagram and a real and a
    // dummy loop in the other.
    let real_loops = self_loops(false)?;
    let dummy_loop = self_loops(true)?;
    let graph = real_loops.hedge_graph();
    let both_loops = InternalSubGraph::cleaned_filter_optimist(graph.full_filter(), graph);
    let foreign = Spinney::new(both_loops, &real_loops, real_loops.loop_momentum_basis())?
        .expect("two quartic tadpoles are divergent");
    assert!(matches!(
        Wood::from_spinneys([foreign], &dummy_loop, dummy_loop.loop_momentum_basis()),
        Err(UvGraphError::OutsideDomain { .. })
    ));

    // A spinney of a differently sized diagram is foreign whatever its domain.
    let basketball = diagram(
        Model::phi4(),
        r#"digraph basketball { edge [particle=phi]; a -> b; a -> b; a -> b; a -> b; }"#,
    );
    let bubble = spinneys(&basketball)
        .into_iter()
        .find(|spinney| spinney.loop_count() == 1)
        .expect("the basketball has one-loop spinneys");
    assert!(matches!(
        Wood::from_spinneys([bubble], &real_loops, real_loops.loop_momentum_basis()),
        Err(UvGraphError::ForeignSpinney {
            reason: ForeignReason::Width,
            ..
        })
    ));
    Ok(())
}

/// FeynKit's pruning degree is never below its classifying degree, so a union it
/// classifies as divergent always passes pruning; it is enumerated here because a
/// one-loop union is reached directly (see the Petersen fixture for the separate
/// reachability limit). A scalar triangle with a vertex
/// factor that vanishes after routing, a literal zero factor, and a squared
/// energy of another line as a factor on edge 0.
#[test]
fn pruning_keeps_every_union_feynkit_classifies_as_divergent() -> Result<(), Box<dyn Error>> {
    let triangle = FeynmanDiagram::from_dot(
        Model::phi3(),
        r#"digraph triangle {
            ext [style=invis]; edge [particle=phi];
            a -> b [id=0]; a -> c [id=1]; b -> c [id=2];
            ext -> a [id=3]; b -> ext [id=4]; c -> ext [id=5];
        }"#,
    )?;
    let a = triangle
        .vertices()
        .find(|(_, vertex)| vertex.name == "a")
        .expect("the triangle names vertex a")
        .0;
    let energy = |edge: usize| {
        function!(
            momentum_symbol(),
            edge,
            function!(symbol!("spenso::cind"), 0)
        )
    };
    let vertex_factor = |factor| -> Result<FeynmanDiagram, Box<dyn Error>> {
        let fragment = NumeratorInterface::admit(factor, None)?;
        Ok(triangle
            .clone()
            .with_numerator_fragments([(a, fragment)], [] as [(EdgeId, _); 0])?)
    };
    let edge_factor = NumeratorInterface::admit(energy(1).pow(2), None)?;
    let cases = [
        (
            "vanishes after routing",
            vertex_factor(energy(0) + energy(1) - energy(3))?,
        ),
        ("literal zero", vertex_factor(Atom::Zero)?),
        (
            "another line's energy on edge 0",
            triangle
                .clone()
                .with_numerator_fragments([], [(EdgeId(0), edge_factor)])?,
        ),
    ];
    for (label, weighted) in cases {
        let loops = weighted.internal_subgraph();
        let classifying = weighted.degree_of_divergence(&loops)?;
        assert_eq!(classifying, 0, "{label}");
        assert!(
            weighted.local_degree_of_divergence(&loops)? >= classifying,
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

/// The energy component of edge `edge`'s momentum.
fn edge_energy(edge: usize) -> Atom {
    function!(
        momentum_symbol(),
        edge,
        function!(symbol!("spenso::cind"), 0)
    )
}

/// A cubic scalar bubble; `reversed` orients its second line the other way.
fn cubic_bubble(reversed: bool) -> Result<FeynmanDiagram, Box<dyn Error>> {
    let second = if reversed { "b -> a" } else { "a -> b" };
    Ok(FeynmanDiagram::from_dot(
        Model::phi3(),
        &format!(
            "digraph bubble {{ ext [style=invis]; edge [particle=phi]; \
             a -> b [id=0]; {second} [id=1]; ext -> a [id=2]; b -> ext [id=3]; }}"
        ),
    )?)
}

#[test]
fn pruning_degree_is_the_true_maximum() -> Result<(), Box<dyn Error>> {
    // Before: the superficial degree 0 was returned without the classifying 4.
    let bubble = cubic_bubble(false)?;
    let quartic = NumeratorInterface::admit(edge_energy(1).pow(4), None)?;
    let bubble = bubble.with_numerator_fragments([], [(EdgeId(0), quartic)])?;
    let loops = bubble.internal_subgraph();
    assert_eq!(bubble.local_degree_of_divergence(&loops)?, 4);
    assert_eq!(bubble.degree_of_divergence(&loops)?, 4);
    assert!(bubble.passes_pruning(&loops)?);

    // A negative superficial degree raised to the classifying one.
    let triangle = FeynmanDiagram::from_dot(
        Model::phi3(),
        r#"digraph triangle {
            ext [style=invis]; edge [particle=phi];
            a -> b [id=0]; a -> c [id=1]; b -> c [id=2];
            ext -> a [id=3]; b -> ext [id=4]; c -> ext [id=5];
        }"#,
    )?;
    let cubic = NumeratorInterface::admit(edge_energy(1).pow(3), None)?;
    let triangle = triangle.with_numerator_fragments([], [(EdgeId(0), cubic)])?;
    let loops = triangle.internal_subgraph();
    assert_eq!(triangle.local_degree_of_divergence(&loops)?, 1);
    assert_eq!(triangle.degree_of_divergence(&loops)?, 1);
    Ok(())
}

#[test]
fn algebraically_vanishing_numerator_is_kept_with_degree_zero() -> Result<(), Box<dyn Error>> {
    let bubble = cubic_bubble(false)?;
    let a = bubble
        .vertices()
        .find(|(_, vertex)| vertex.name == "a")
        .expect("the bubble names vertex a")
        .0;
    let (q0, q1) = (edge_energy(0), edge_energy(1));
    let identity = (&q0 + &q1).pow(2) - q0.pow(2) - Atom::num(2) * &q0 * &q1 - q1.pow(2);
    for (label, factor) in [
        ("algebraic identity", identity),
        ("literal zero", Atom::Zero),
    ] {
        let fragment = NumeratorInterface::admit(factor, None)?;
        let weighted = bubble
            .clone()
            .with_numerator_fragments([(a, fragment)], [] as [(EdgeId, _); 0])?;
        let loops = weighted.internal_subgraph();
        // Before: the identity gave pruning 2 and classifying −3.
        assert_eq!(weighted.degree_of_divergence(&loops)?, 0, "{label}");
        assert!(weighted.local_degree_of_divergence(&loops)? >= 0, "{label}");
        let kept = spinneys(&weighted)
            .into_iter()
            .any(|spinney| spinney.filter() == &loops);
        assert!(kept, "{label}");
    }
    Ok(())
}

#[test]
fn woods_reject_spinneys_routed_for_another_orientation() -> Result<(), Box<dyn Error>> {
    // Same selection, components and loop count; the second line is reversed,
    // so its loop signature differs.
    let aligned = cubic_bubble(false)?;
    let reversed = cubic_bubble(true)?;
    let graph = aligned.hedge_graph();
    let loops = InternalSubGraph::cleaned_filter_optimist(aligned.internal_subgraph(), graph);
    let foreign = Spinney::new(loops, &aligned, aligned.loop_momentum_basis())?
        .expect("the cubic bubble is logarithmic");
    assert!(matches!(
        Wood::from_spinneys([foreign.clone()], &reversed, reversed.loop_momentum_basis()),
        Err(UvGraphError::ForeignSpinney {
            reason: ForeignReason::SubBasis,
            ..
        })
    ));
    assert!(Wood::from_spinneys([foreign], &aligned, aligned.loop_momentum_basis()).is_ok());
    Ok(())
}

#[test]
fn direct_constructors_reject_foreign_inputs() -> Result<(), Box<dyn Error>> {
    let bubble = cubic_bubble(false)?;
    let small = self_loops(false)?;
    // A selection of another, smaller diagram.
    let small_loops =
        InternalSubGraph::cleaned_filter_optimist(small.internal_subgraph(), small.hedge_graph());
    assert!(matches!(
        Spinney::with_scheme(
            small_loops,
            &bubble,
            bubble.loop_momentum_basis(),
            ApproximationType::MUV,
            0,
        ),
        Err(UvGraphError::ForeignSpinney {
            reason: ForeignReason::Width,
            ..
        })
    ));
    // A valid selection with the parent basis of another, smaller diagram.
    let loops =
        InternalSubGraph::cleaned_filter_optimist(bubble.internal_subgraph(), bubble.hedge_graph());
    assert!(matches!(
        Spinney::new(loops, &bubble, small.loop_momentum_basis()),
        Err(UvGraphError::CompatibleBasis {
            source: DiagramError::InvalidLoopMomentumBasis(_),
            ..
        })
    ));
    Ok(())
}

#[test]
fn forest_operations_reject_a_smaller_diagram() -> Result<(), Box<dyn Error>> {
    let triangle = FeynmanDiagram::from_dot(
        Model::phi3(),
        r#"digraph triangle {
            ext [style=invis]; edge [particle=phi];
            a -> b [id=0]; a -> c [id=1]; b -> c [id=2];
            ext -> a [id=3]; b -> ext [id=4]; c -> ext [id=5];
        }"#,
    )?;
    let cubic = NumeratorInterface::admit(edge_energy(1).pow(3), None)?;
    let triangle = triangle.with_numerator_fragments([], [(EdgeId(0), cubic)])?;
    let forests = forests(&triangle);
    let smaller = self_loops(false)?;
    let (_, operation) = forests
        .operations()
        .find(|(_, operation)| operation.depth() > 0)
        .expect("the weighted triangle is divergent");
    assert!(matches!(
        operation.forest_node(&smaller, smaller.loop_momentum_basis(), 0),
        Err(OperationSpinneyError::Spinney {
            source: UvGraphError::ForeignSpinney {
                reason: ForeignReason::Width,
                ..
            },
            ..
        })
    ));
    let width = forests.graph().n_hedges();
    assert!(matches!(
        forests.compatible_topological_order(&SuBitGraph::empty(width + 1)),
        Err(ForestError::SelectionWidth { found, expected }) if found == width + 1 && expected == width
    ));
    Ok(())
}

/// The φ⁴ sunset has three loop momentum bases. A spinney classified from one
/// and put into a wood with another is accepted exactly when its routed
/// sub-basis coincides: never for the full graph, sometimes for a single loop.
#[test]
fn woods_accept_spinneys_of_another_basis_only_when_the_sub_basis_coincides()
-> Result<(), Box<dyn Error>> {
    let sunset = diagram(
        Model::phi4(),
        r#"digraph sunset {
            ext [style=invis]; edge [particle=phi];
            a -> b [id=0]; a -> b [id=1]; a -> b [id=2];
            ext -> a [id=3]; b -> ext [id=4];
        }"#,
    );
    let full = sunset.hedge_graph().full_filter();
    let bases = sunset.loop_momentum_bases_of(&full, 4)?;
    assert_eq!(bases.len(), 3);
    let mut accepted_one_loop = 0;
    for (classified, wood) in (0..bases.len()).cartesian_product(0..bases.len()) {
        if classified == wood {
            continue;
        }
        let spinneys = sunset.classify_spinneys(
            &full,
            sunset.model(),
            &RenormalizationPrescriptionSettings::default(),
            &bases[classified],
        )?;
        for spinney in spinneys.into_iter().filter(|spinney| !spinney.is_empty()) {
            let coincides =
                sunset.compatible_sub_basis(spinney.subgraph(), &bases[wood])? == *spinney.lmb();
            let loop_count = spinney.loop_count();
            let result = Wood::from_spinneys([spinney], &sunset, &bases[wood]);
            if loop_count == 2 {
                assert!(!coincides);
                assert!(matches!(
                    result,
                    Err(UvGraphError::ForeignSpinney {
                        reason: ForeignReason::SubBasis,
                        ..
                    })
                ));
            } else if coincides {
                assert!(result.is_ok());
                accepted_one_loop += 1;
            } else {
                assert!(matches!(
                    result,
                    Err(UvGraphError::ForeignSpinney {
                        reason: ForeignReason::SubBasis,
                        ..
                    })
                ));
            }
        }
    }
    assert!(accepted_one_loop > 0);
    Ok(())
}

/// Known limitation: `g_μν (q₀^μ q₃^ν − q₀^ν q₃^μ)` across two vertices is zero
/// only once the indices are contracted, which classification does not do. On a
/// triangle, where `q₃` is external, the routed degree is −1 and the triangle is
/// dropped although its integrand is zero; on a bubble, where both momenta are
/// loop lines, the routed degree is 1 and the bubble is kept. Contracting every
/// candidate was measured too costly (the gluon Petersen test exceeds 240 s).
#[test]
fn contraction_only_zeros() -> Result<(), Box<dyn Error>> {
    let mu = function!(symbol!("spenso::mink"), 4, symbol!("feynkit_uv_tests::mu"));
    let nu = function!(symbol!("spenso::mink"), 4, symbol!("feynkit_uv_tests::nu"));
    let metric = function!(symbol!("spenso::g"), &mu, &nu);
    let q = |edge: usize, index: &Atom| function!(momentum_symbol(), edge, index);
    let weighted =
        |diagram: FeynmanDiagram, other: usize| -> Result<FeynmanDiagram, Box<dyn Error>> {
            let vertex = |name: &str| {
                diagram
                    .vertices()
                    .find(|(_, vertex)| vertex.name == name)
                    .expect("the fixture names its vertices")
                    .0
            };
            let (a, b) = (vertex("a"), vertex("b"));
            let antisymmetric = q(0, &mu) * q(other, &nu) - q(0, &nu) * q(other, &mu);
            Ok(diagram.with_numerator_fragments(
                [
                    (a, NumeratorInterface::admit(metric.clone(), None)?),
                    (b, NumeratorInterface::admit(antisymmetric, None)?),
                ],
                [] as [(EdgeId, _); 0],
            )?)
        };
    let triangle = FeynmanDiagram::from_dot(
        Model::phi3(),
        r#"digraph triangle {
            ext [style=invis]; edge [particle=phi];
            a -> b [id=0]; a -> c [id=1]; b -> c [id=2];
            ext -> a [id=3]; b -> ext [id=4]; c -> ext [id=5];
        }"#,
    )?;
    for (label, diagram, degree, retained) in [
        ("triangle", weighted(triangle, 3)?, -1, false),
        ("bubble", weighted(cubic_bubble(false)?, 1)?, 1, true),
    ] {
        let loops = diagram.internal_subgraph();
        assert_eq!(diagram.degree_of_divergence(&loops)?, degree, "{label}");
        let kept = spinneys(&diagram)
            .into_iter()
            .any(|spinney| spinney.filter() == &loops);
        assert_eq!(kept, retained, "{label}");
    }
    Ok(())
}

#[test]
fn spinney_relations_classify_nesting_disjointness_and_overlap() {
    let sunrise = diagram(
        Model::phi4(),
        r#"digraph sunrise { edge [particle="phi"]; ext [style=invis];
            ext -> a; a -> b; a -> b; a -> b; b -> ext; }"#,
    );
    let sorted = spinneys(&sunrise)
        .into_iter()
        .filter(|spinney| !spinney.is_empty())
        .sorted_by_key(|spinney| (spinney.loop_count(), spinney.edges(&sunrise)))
        .collect::<Vec<_>>();
    let [bubble, other, _, full] = &sorted[..] else {
        panic!("the sunrise has three bubbles and the full graph");
    };
    assert_eq!(
        bubble.relation(full, &sunrise),
        SpinneyRelation::ContainedIn
    );
    assert_eq!(full.relation(bubble, &sunrise), SpinneyRelation::Contains);
    assert_eq!(bubble.relation(bubble, &sunrise), SpinneyRelation::Equal);
    assert_eq!(
        bubble.relation(other, &sunrise),
        SpinneyRelation::Overlapping
    );

    let spectacles = diagram(
        Model::phi3(),
        r#"digraph { a -> a [particle="phi"]; a -> b [particle="phi"]; b -> b [particle="phi"]; }"#,
    );
    let tadpoles = spinneys(&spectacles)
        .into_iter()
        .filter(|spinney| spinney.n_components() == 1 && !spinney.is_empty())
        .collect::<Vec<_>>();
    let relation = tadpoles[0].relation(&tadpoles[1], &spectacles);
    assert_eq!(relation, SpinneyRelation::Disjoint);
    assert!(relation.is_compatible());
}

#[test]
fn bubbles_sharing_a_vertex_touch_and_form_no_product_forest() {
    // The s-channel chain of two quartic bubbles: they share vertex b but no
    // line. Disjoint subgraphs of the forest formula share no vertex, so the
    // bubbles only nest into the full graph and never multiply.
    let chain = diagram(
        Model::phi4(),
        r#"digraph chain { edge [particle="phi"]; ext [style=invis]; ext -> a; ext -> a;
            a -> b; a -> b; b -> c; b -> c; c -> ext; c -> ext; }"#,
    );
    let all = spinneys(&chain);
    let bubbles = all
        .iter()
        .filter(|spinney| spinney.loop_count() == 1)
        .collect::<Vec<_>>();
    assert_eq!(bubbles.len(), 2);
    let relation = bubbles[0].relation(bubbles[1], &chain);
    assert_eq!(relation, SpinneyRelation::Touching);
    assert!(!relation.is_compatible());
    assert!(matches!(
        Counterterm::new(
            bubbles.iter().map(|spinney| (*spinney).clone()),
            CountertermScheme::MsBar,
            &chain,
            chain.loop_momentum_basis(),
        ),
        Err(CountertermError::Incompatible {
            relation: SpinneyRelation::Touching,
            ..
        })
    ));
    // ∅, the two bubbles, the full graph and each bubble nested in it.
    let forests = forests(&chain);
    assert_eq!(forests.n_forests(), 6);
    assert!(
        forests
            .operations()
            .all(|(_, forest)| forest.levels().iter().all(|level| level.len() == 1))
    );
}

#[test]
fn forest_transitions_chain_their_spinneys() {
    let basketball = diagram(
        Model::phi4(),
        r#"digraph basketball { edge [particle="phi"]; a -> b; a -> b; a -> b; a -> b; }"#,
    );
    let forests = forests(&basketball);
    let wood = forests.wood();
    for (node, forest) in forests.operations() {
        let transitions = forest.transitions(wood);
        assert_eq!(transitions.len(), forest.levels().len());
        let elements = forests.forest_spinneys(node);
        assert_eq!(elements.len(), forest.depth());
        let mut treated = vec![wood.root()];
        for (level, transition) in transitions.iter().flatten().enumerate() {
            // Each operation approximates a larger spinney after one already treated.
            assert!(treated.contains(&transition.smaller), "level {level}");
            assert_eq!(wood.spinney(transition.larger), elements[level]);
            assert!(wood.spinney(transition.smaller) < wood.spinney(transition.larger));
            treated.push(transition.larger);
        }
        // Forest elements are pairwise compatible.
        for (a, b) in elements.iter().tuple_combinations() {
            assert!(a.relation(*b, &basketball).is_compatible());
        }
    }
    let sorted = wood.nodes_by_edges(&basketball);
    assert_eq!(sorted[0], wood.root());
    assert_eq!(sorted.len(), wood.n_spinneys());
}

#[test]
fn counterterms_contract_compatible_spinneys_outermost_first() {
    let sunrise = diagram(
        Model::phi4(),
        r#"digraph sunrise { edge [particle="phi"]; ext [style=invis];
            ext -> a; a -> b; a -> b; a -> b; b -> ext; }"#,
    );
    let mut sorted = spinneys(&sunrise);
    sorted.sort_by_key(|spinney| (spinney.loop_count(), spinney.edges(&sunrise)));
    let (empty, bubble, other, full) = (
        sorted[0].clone(),
        sorted[1].clone(),
        sorted[2].clone(),
        sorted[4].clone(),
    );
    assert!(empty.is_empty());
    let counterterm = Counterterm::new(
        [bubble.clone(), empty.clone(), full.clone()],
        CountertermScheme::OnShell,
        &sunrise,
        sunrise.loop_momentum_basis(),
    )
    .unwrap();
    // A spinney selected twice is no forest.
    assert!(matches!(
        Counterterm::new(
            [bubble.clone(), full.clone(), bubble.clone()],
            CountertermScheme::OnShell,
            &sunrise,
            sunrise.loop_momentum_basis(),
        ),
        Err(CountertermError::Duplicate { spinney }) if spinney == "γ{0,1}"
    ));
    assert_eq!(counterterm.contracted(), [full.clone(), bubble.clone()]);
    assert_eq!(counterterm.scheme(), CountertermScheme::OnShell);
    assert_eq!(
        counterterm.with_scheme(CountertermScheme::MsBar).scheme(),
        CountertermScheme::MsBar
    );
    let basis = sunrise.loop_momentum_basis();
    assert!(matches!(
        Counterterm::new([empty], CountertermScheme::Unspecified, &sunrise, basis),
        Err(CountertermError::Empty)
    ));
    let incompatible = Counterterm::new([bubble, other], CountertermScheme::MsBar, &sunrise, basis);
    assert!(matches!(
        incompatible,
        Err(CountertermError::Incompatible {
            relation: SpinneyRelation::Overlapping,
            ..
        })
    ));
    // Errors name spinneys by their edges.
    assert_eq!(
        incompatible.unwrap_err().to_string(),
        "spinneys γ{0,1} and γ{0,2} are overlapping and cannot be contracted together"
    );

    let forests = forests(&sunrise);
    assert!(matches!(
        Counterterm::of_forest::<DiagramError>(&forests, forests.root(), CountertermScheme::MsBar),
        Err(CountertermError::Empty)
    ));
    for (node, forest) in forests
        .operations()
        .filter(|(_, forest)| forest.depth() > 0)
    {
        let counterterm =
            Counterterm::of_forest::<DiagramError>(&forests, node, CountertermScheme::MsBar)
                .unwrap();
        assert_eq!(counterterm.contracted().len(), forest.depth());
    }
    for name in ["MS-bar", "on-SHELL", "unspecified"] {
        let scheme: CountertermScheme = name.parse().unwrap();
        assert!(scheme.to_string().eq_ignore_ascii_case(name));
    }
}

/// The forests of `diagram` whose elements are a product of vertex-disjoint
/// spinneys, with those elements as sorted edge lists per level.
fn product_levels(diagram: &FeynmanDiagram) -> Vec<Vec<Vec<Vec<usize>>>> {
    let forests = forests(diagram);
    forests
        .operations()
        .filter(|(node, _)| forests.is_disjoint_union(*node))
        .map(|(node, _)| {
            forests
                .forest_levels(node)
                .into_iter()
                .map(|level| {
                    level
                        .into_iter()
                        .map(|spinney| edges(diagram, spinney.filter()))
                        .sorted()
                        .collect()
                })
                .collect()
        })
        .sorted()
        .collect()
}

#[test]
fn product_forests_list_each_disjoint_factor() {
    // Cubic spectacles: the tadpoles 0 and 2 multiply in one commuting level.
    let spectacles = diagram(
        Model::phi3(),
        r#"digraph { a -> a [id=0, particle="phi"]; a -> b [id=1, particle="phi"]; b -> b [id=2, particle="phi"]; }"#,
    );
    assert_eq!(product_levels(&spectacles), [[[vec![0], vec![2]]]]);

    // Two quartic bubbles joined by a line, and three in a chain.
    let pair = diagram(
        Model::phi4(),
        r#"digraph pair { edge [particle="phi"]; ext [style=invis];
            a -> b [id=0]; a -> b [id=1]; b -> c [id=2]; c -> d [id=3]; c -> d [id=4];
            ext -> a; ext -> a; b -> ext; c -> ext; d -> ext; d -> ext; }"#,
    );
    assert_eq!(product_levels(&pair), [[[vec![0, 1], vec![3, 4]]]]);
    let triple = diagram(
        Model::phi4(),
        r#"digraph triple { edge [particle="phi"]; ext [style=invis];
            a -> b [id=0]; a -> b [id=1]; b -> c [id=2]; c -> d [id=3]; c -> d [id=4];
            d -> e [id=5]; e -> f [id=6]; e -> f [id=7];
            ext -> a; ext -> a; b -> ext; c -> ext; d -> ext; e -> ext; f -> ext; f -> ext; }"#,
    );
    let products = product_levels(&triple);
    assert_eq!(products.len(), 4);
    assert!(products.contains(&vec![vec![vec![0, 1], vec![3, 4], vec![6, 7]]]));

    for diagram in [&spectacles, &pair, &triple] {
        let forests = forests(diagram);
        let basis = diagram.loop_momentum_basis();
        for (node, forest) in forests.operations().filter(|(_, f)| f.depth() > 0) {
            let elements = forests.forest_spinneys(node);
            // One connected element per operation, pairwise compatible.
            assert_eq!(elements.len(), forest.depth());
            assert!(elements.iter().all(|spinney| spinney.n_components() == 1));
            for (a, b) in elements.iter().tuple_combinations() {
                assert!(a.relation(*b, diagram).is_compatible());
            }
            // The forest's counterterm contracts what a hand-built one contracts.
            let of_forest =
                Counterterm::of_forest::<DiagramError>(&forests, node, CountertermScheme::MsBar)
                    .unwrap();
            let by_hand = Counterterm::new(
                elements.into_iter().cloned(),
                CountertermScheme::MsBar,
                diagram,
                basis,
            )
            .unwrap();
            assert_eq!(of_forest.contracted(), by_hand.contracted());
        }
    }
}

#[test]
fn counterterms_reject_spinneys_of_another_graph() {
    let basketball = diagram(
        Model::phi4(),
        r#"digraph basketball { edge [particle="phi"]; a -> b; a -> b; a -> b; a -> b; }"#,
    );
    let bubble = diagram(
        Model::phi4(),
        r#"digraph bubble { edge [particle="phi"]; ext [style=invis];
            ext -> a; ext -> a; a -> b; a -> b; b -> ext; b -> ext; }"#,
    );
    let full = spinneys(&basketball)
        .into_iter()
        .max_by_key(Spinney::loop_count)
        .unwrap();
    // Even a single spinney is checked against the destination graph. Both
    // diagrams have eight half-edges; the basketball's lines are not the
    // bubble's.
    let foreign = Counterterm::new(
        [full],
        CountertermScheme::MsBar,
        &bubble,
        bubble.loop_momentum_basis(),
    );
    assert!(
        matches!(foreign, Err(CountertermError::Foreign(_))),
        "{:?}",
        foreign.map(|counterterm| counterterm.contracted().len())
    );
    // A spinney of a smaller diagram has another width.
    let tadpole = diagram(
        Model::phi3(),
        r#"digraph { a -> a [particle="phi"]; a -> b [particle="phi"]; b -> b [particle="phi"]; }"#,
    );
    let small = spinneys(&tadpole)
        .into_iter()
        .find(|spinney| !spinney.is_empty())
        .unwrap();
    assert!(matches!(
        Counterterm::new(
            [small],
            CountertermScheme::MsBar,
            &bubble,
            bubble.loop_momentum_basis()
        ),
        Err(CountertermError::Foreign(UvGraphError::ForeignSpinney {
            reason: ForeignReason::Width,
            ..
        }))
    ));
}

#[test]
fn disconnected_spinneys_contract_as_their_components() {
    let spectacles = diagram(
        Model::phi3(),
        r#"digraph { a -> a [id=0, particle="phi"]; a -> b [id=1, particle="phi"]; b -> b [id=2, particle="phi"]; }"#,
    );
    let basis = spectacles.loop_momentum_basis();
    let classified = spinneys(&spectacles);
    let union = classified
        .iter()
        .find(|spinney| spinney.n_components() == 2)
        .expect("the two tadpoles form a disconnected union")
        .clone();
    let left = classified
        .iter()
        .find(|spinney| edges(&spectacles, spinney.filter()) == [0])
        .unwrap()
        .clone();
    // The union contracts as both tadpoles, as its product forest does.
    let by_union = Counterterm::new(
        [union.clone()],
        CountertermScheme::MsBar,
        &spectacles,
        basis,
    )
    .unwrap();
    let forests = forests(&spectacles);
    let (product, _) = forests
        .operations()
        .find(|(node, _)| forests.is_disjoint_union(*node))
        .unwrap();
    let of_forest =
        Counterterm::of_forest::<DiagramError>(&forests, product, CountertermScheme::MsBar)
            .unwrap();
    let labels = |counterterm: &Counterterm<LoopMomentumBasis>| {
        counterterm
            .contracted()
            .iter()
            .map(|spinney| edges(&spectacles, spinney.filter()))
            .collect::<Vec<_>>()
    };
    assert_eq!(labels(&by_union), labels(&of_forest));
    assert_eq!(by_union.contracted().len(), 2);
    // A tadpole together with the union holding it selects the tadpole twice.
    assert!(matches!(
        Counterterm::new([left, union], CountertermScheme::MsBar, &spectacles, basis),
        Err(CountertermError::Duplicate { .. })
    ));
}

#[test]
fn disjoint_products_close_disconnected_and_repeated_inputs() {
    // Three quartic tadpoles in a chain a–b–c; loops 0, 1 and 2 share no vertex.
    let chain = diagram(
        Model::phi4(),
        r#"digraph { edge [particle="phi"]; ext [style=invis];
            a -> a [id=0]; b -> b [id=1]; c -> c [id=2]; a -> b [id=3]; b -> c [id=4];
            ext -> a; c -> ext; }"#,
    );
    let basis = chain.loop_momentum_basis();
    let classified = spinneys(&chain);
    let classified_forests = Wood::from_spinneys(classified.clone(), &chain, basis)
        .unwrap()
        .unfold()
        .n_forests();
    assert_eq!(classified_forests, 8);
    let find = |wanted: &[usize]| {
        classified
            .iter()
            .find(|spinney| edges(&chain, spinney.filter()) == wanted)
            .unwrap()
            .clone()
    };
    // The single tadpoles close to the classified wood.
    let mut singles = vec![find(&[0]), find(&[1]), find(&[2])];
    let products = Wood::disjoint_products(&singles, &chain, basis).unwrap();
    assert_eq!(products.len(), 4);
    singles.extend(products);
    let closed = Wood::from_spinneys(singles, &chain, basis)
        .unwrap()
        .unfold();
    assert_eq!(closed.n_forests(), classified_forests);
    // A disconnected input X = {0, 1} multiplies with Y = {2}: ∅, X, Y and X·Y.
    let mut given = vec![find(&[0, 1]), find(&[2])];
    let products = Wood::disjoint_products(&given, &chain, basis).unwrap();
    assert_eq!(
        products
            .iter()
            .map(|spinney| edges(&chain, spinney.filter()))
            .collect::<Vec<_>>(),
        [vec![0, 1, 2]]
    );
    given.extend(products);
    let forests = Wood::from_spinneys(given, &chain, basis).unwrap().unfold();
    assert_eq!(forests.n_forests(), 4);
    // The product forest multiplies the two given factors with a positive sign.
    let (product, operation) = forests
        .operations()
        .find(|(node, _)| forests.is_disjoint_union(*node))
        .expect("X and Y multiply");
    assert_eq!(operation.sign(), 1);
    let elements: Vec<_> = forests
        .forest_spinneys(product)
        .into_iter()
        .map(|spinney| edges(&chain, spinney.filter()))
        .sorted()
        .collect();
    assert_eq!(elements, [vec![0, 1], vec![2]]);
    // Repeated inputs build each union once.
    let repeated: Vec<_> = std::iter::repeat_n(find(&[0]), 50)
        .chain(std::iter::repeat_n(find(&[2]), 50))
        .collect();
    let products = Wood::disjoint_products(&repeated, &chain, basis).unwrap();
    assert_eq!(products.len(), 1);
}
