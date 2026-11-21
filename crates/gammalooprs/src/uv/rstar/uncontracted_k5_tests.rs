//! U7c tests: the K5 insertion-surrogate records (kind `surrogate`,
//! [`super::uncontracted::k5`]).
//!
//! The ORACLE is the beta's verified fixtures (wh-fixtures k5-direct,
//! k5-surrogate, k5-stress): the record the emitter builds for the fixture's
//! site is compared field by field with the fixture's `input.txt` (atoms by
//! expanded difference, the label ignored), read back with vakint's frozen
//! reader (which validates it, R2 included), evaluated by vakint's entry and
//! finished as the hook finishes it; that value must equal the fixture's
//! CONTRACTED value (`value.txt`) through its `valid_through`. The surrogate
//! graphs come from this build's own `prepare` of the beta's prepared DOT;
//! the factors of the record's edges and vertices (the beta's colour branch)
//! are taken from the fixture, as the U7a oracle tests do ([`with_branch`]):
//! the factor lines then test the emitter's pass-through, the momenta, legs,
//! insertions, lines, cycles and ports its construction.
//!
//! Then one test per binding condition of the design (notes/wi/U7C_SITE.md;
//! beta's answers U7C_BETA_ANSWERS.md): the closed designed list, the loud
//! errors, the zero-comparison failures (diagram and run), the ledger frame,
//! the default kind list, the Δ-cache key, the off mode.
//!
//! k5-stress (five L5 3793 records) is vendored too (testdata/u7c/k5-stress-*,
//! sha256 equal to the fixture's VERIFIED record); the whole-graph L5 runs of
//! 22807 / 3793 are `#[ignore]`d (minutes; run with --run-ignored), the
//! multi-loop-hole graphs run by default.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use eyre::{Result, ensure, eyre};
use idenso::tensor::SymbolicTensor;
use linnet::half_edge::subgraph::SuBitGraph;
use symbolica::atom::{Atom, AtomCore};
use vakint::{
    Vakint,
    p_integral::factorized::{
        ClosureSpec, FactorizedKind, FactorizedNumerator, TaylorOrders, TaylorSpec,
        evaluate_factorized_ports,
    },
};

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt, graph::Graph, initialisation::test_initialise,
    numerator::symbolica_ext::NumeratorAtomExt, utils::GS,
    uv::approx::streamed::with_production_piece_limit,
};

use super::{
    DiagramLedger, PropagatorInsertionEvaluator, REFERENCE_MASS_SQUARED, RStarPlan,
    boundary::ClosedRootNormalization,
    expression::{RStarIntegrand, in_dressed_scope},
    plan::{join_edge_ids, uv_child_edge_ids},
    recursive_tests::{colour_factored, load, settings},
    top_mass::TopMassTerms,
    uncontracted::{
        self, Captured, IrHardSite, Kind, Mode, Outcome, Site, Switch,
        k5::{self, Provenance},
        with_capture, with_switch,
    },
    uncontracted_ledger,
    uncontracted_oracle_tests::{
        assume_defaults, atoms_equal, auxiliary_mass, equal_through, fixture_value, oracle_value,
        record, record_differences, valid_through,
    },
};

const L2_17_DOT: &str = include_str!("testdata/u7c/l2-d017.dot");
const L5_22807_DOT: &str = include_str!("testdata/u7c/l5-d22807.dot");
const L5_3793_DOT: &str = include_str!("testdata/u7c/l5-d3793.dot");
const K5_DIRECT: &str = include_str!("testdata/u7c/k5-direct.input.txt");
const K5_DIRECT_VALUE: &str = include_str!("testdata/u7c/k5-direct.value.txt");
const K5_SURROGATE: &str = include_str!("testdata/u7c/k5-surrogate.input.txt");
const K5_SURROGATE_VALUE: &str = include_str!("testdata/u7c/k5-surrogate.value.txt");
/// k5-stress r0051..r0055 (input, value), vendored from the beta's wh-fixtures.
const K5_STRESS: [(&str, &str, &str); 5] = [
    (
        "r0051",
        include_str!("testdata/u7c/k5-stress-r0051.input.txt"),
        include_str!("testdata/u7c/k5-stress-r0051.value.txt"),
    ),
    (
        "r0052",
        include_str!("testdata/u7c/k5-stress-r0052.input.txt"),
        include_str!("testdata/u7c/k5-stress-r0052.value.txt"),
    ),
    (
        "r0053",
        include_str!("testdata/u7c/k5-stress-r0053.input.txt"),
        include_str!("testdata/u7c/k5-stress-r0053.value.txt"),
    ),
    (
        "r0054",
        include_str!("testdata/u7c/k5-stress-r0054.input.txt"),
        include_str!("testdata/u7c/k5-stress-r0054.value.txt"),
    ),
    (
        "r0055",
        include_str!("testdata/u7c/k5-stress-r0055.input.txt"),
        include_str!("testdata/u7c/k5-stress-r0055.value.txt"),
    ),
];
/// L5 surrogates with holes of two or more loops (coordinator, 21 Nov):
/// (graph, prepared DOT).
const L5_MULTI_LOOP_HOLES: [(&str, &str); 5] = [
    ("143019", include_str!("testdata/u7c/l5-d143019.dot")),
    ("143190", include_str!("testdata/u7c/l5-d143190.dot")),
    ("163502", include_str!("testdata/u7c/l5-d163502.dot")),
    ("175270", include_str!("testdata/u7c/l5-d175270.dot")),
    ("139207", include_str!("testdata/u7c/l5-d139207.dot")),
];

fn check_surrogate() -> Switch {
    Switch::new(Mode::Check, &[Kind::Surrogate])
}

fn on_surrogate() -> Switch {
    Switch::new(Mode::On, &[Kind::Surrogate])
}

// ---------------------------------------------------------------------------
// Fixture plumbing
// ---------------------------------------------------------------------------

/// A surrogate of `dot` prepared by this build (switch `switch`, so that its
/// provenance is computed), and that provenance.
struct Prepared {
    surrogate: feynkit_graph::FeynmanDiagram,
    provenance: Arc<Provenance>,
}

fn prepare(dot: &str, vakint: &Vakint) -> Result<Prepared> {
    let original = colour_factored(&load(dot)?)?;
    let evaluator = PropagatorInsertionEvaluator::new(settings());
    with_switch(check_surrogate(), || -> Result<Prepared> {
        let prepared = evaluator
            .prepare(&original, vakint)?
            .ok_or_else(|| eyre!("no insertion selected"))?;
        let provenance = k5::provenance(
            &original,
            prepared
                .holes()
                .iter()
                .map(|hole| (hole.vertex.0, hole.loops, &hole.original.support)),
        )?
        .ok_or_else(|| eyre!("no provenance with the kind enabled"))?;
        Ok(Prepared {
            surrogate: prepared.surrogate().clone(),
            provenance,
        })
    })
}

/// The surrogate with the factors the fixture's record states for its own
/// edges and vertices (the beta's colour branch; an `ins:` vertex takes its
/// insertion's coefficient x factor).
fn with_branch(
    surrogate: &feynkit_graph::FeynmanDiagram,
    fixture: &FactorizedNumerator,
) -> Result<feynkit_graph::FeynmanDiagram> {
    let mut vertices = BTreeMap::new();
    for vertex in &fixture.vertices {
        if vertex.rule.starts_with("ct:") {
            continue;
        }
        let factor = if vertex.rule.starts_with("ins:") {
            let insertion = fixture
                .insertions
                .iter()
                .find(|insertion| insertion.vertex == vertex.id)
                .ok_or_else(|| eyre!("no insertion at the ins: vertex {}", vertex.id))?;
            &insertion.coefficient * &insertion.factor
        } else {
            vertex.factor.clone()
        };
        vertices.insert(vertex.id, factor.map_mink_dim(Atom::num(4)));
    }
    let edges = fixture
        .edges
        .iter()
        .filter(|edge| !edge.is_external())
        .map(|edge| (edge.id, edge.factor.map_mink_dim(Atom::num(4))))
        .collect::<BTreeMap<_, _>>();
    let vertex_fragments = surrogate
        .vertices()
        .filter_map(|(id, _)| vertices.get(&id.0).map(|f| (id, f.clone())))
        .map(|(id, f)| Ok((id, SymbolicTensor::infer(f)?)))
        .collect::<Result<Vec<_>>>()?;
    let edge_fragments = surrogate
        .edges()
        .filter_map(|(id, _, _)| edges.get(&id.0).map(|f| (id, f.clone())))
        .map(|(id, f)| Ok((id, SymbolicTensor::infer(f)?)))
        .collect::<Result<Vec<_>>>()?;
    ensure!(
        vertex_fragments.len() == vertices.len() && edge_fragments.len() == edges.len(),
        "the fixture's vertices / edges are not all in the prepared surrogate"
    );
    surrogate
        .clone()
        .with_numerator_fragments(vertex_fragments, edge_fragments)
        .map_err(|error| eyre!("{error}"))
}

/// `name=<value>` of a fixture label.
fn label_field(label: &str, name: &str) -> Option<String> {
    label
        .split_whitespace()
        .find_map(|token| token.strip_prefix(name))
        .map(str::to_owned)
}

/// The child supports of a fixture label (`children=[a,b]loops=..;[c]..`).
fn label_children(label: &str) -> BTreeSet<String> {
    let Some(children) = label_field(label, "children=") else {
        return BTreeSet::new();
    };
    children
        .split('[')
        .skip(1)
        .filter_map(|piece| piece.split(']').next())
        .map(str::to_owned)
        .collect()
}

/// One closed-root remainder of a surrogate above four loops, as the
/// fixture states it: the graph with the fixture's branch, the plan, the
/// term and the fixture's child counterterms by support.
struct RootCase {
    graph: Graph,
    model_diagram: feynkit_graph::FeynmanDiagram,
    plan: RStarPlan,
    term: usize,
    deltas: BTreeMap<SuBitGraph, Atom>,
    fixture: FactorizedNumerator,
    provenance: Arc<Provenance>,
}

fn root_case(prepared: &Prepared, input: &str) -> Result<RootCase> {
    let fixture = record(input)?;
    let diagram = with_branch(&prepared.surrogate, &fixture)?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let plan = RStarPlan::with_automatic_edge(&graph, diagram.model(), &settings())?;
    let FactorizedKind::Closure(spec) = &fixture.kind else {
        return Err(eyre!("the fixture is not a closure record"));
    };
    ensure!(
        usize::from(plan.marked_edge) == spec.edge,
        "the plan's marked edge {} is not the fixture's closing edge {}",
        usize::from(plan.marked_edge),
        spec.edge
    );
    let remaining =
        label_field(&fixture.label, "remaining=").ok_or_else(|| eyre!("no remaining="))?;
    let children = label_children(&fixture.label);
    let term = plan
        .ordinary_uv
        .iter()
        .position(|term| {
            join_edge_ids(&uv_child_edge_ids(&graph, &term.remainder.remaining)) == remaining
                && term
                    .remainder
                    .insertions
                    .iter()
                    .map(|slot| join_edge_ids(&uv_child_edge_ids(&graph, &slot.support)))
                    .collect::<BTreeSet<_>>()
                    == children
        })
        .ok_or_else(|| eyre!("no plan term for remaining={remaining} children={children:?}"))?;
    let mut deltas = BTreeMap::new();
    for slot in &plan.ordinary_uv[term].remainder.insertions {
        let rule = format!(
            "ct:edges-{}",
            uv_child_edge_ids(&graph, &slot.support)
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join("-")
        );
        let vertex = fixture
            .vertices
            .iter()
            .find(|vertex| vertex.rule == rule)
            .ok_or_else(|| eyre!("no {rule} vertex in the fixture"))?
            .id;
        let insertion = fixture
            .insertions
            .iter()
            .find(|insertion| insertion.vertex == vertex)
            .ok_or_else(|| eyre!("no insertion at {vertex}"))?;
        deltas.insert(
            slot.support.clone(),
            &insertion.coefficient * &insertion.factor,
        );
    }
    Ok(RootCase {
        graph,
        model_diagram: diagram,
        plan,
        term,
        deltas,
        fixture,
        provenance: prepared.provenance.clone(),
    })
}

/// The emitter's attempt on the case's remainder as the closed surrogate root
/// holds it (factor 1, σ the physical routing, the root metric closure), in
/// the surrogate context `context` (`None`: outside it); `ir_hard` makes it
/// a K3-shaped term.
fn root_attempt(
    case: &RootCase,
    vakint: &Vakint,
    context: Option<Arc<Provenance>>,
    ir_hard: bool,
) -> Result<uncontracted::Attempt> {
    let uv = settings();
    let vakint_settings = uv.vakint.true_settings();
    let model = case.model_diagram.model();
    let term = &case.plan.ordinary_uv[case.term];
    let closure = RStarIntegrand::root_closure(&case.graph, &case.plan.support)?;
    let (external, external_routing) = RStarIntegrand::physical_external_routing(
        &case.graph,
        &case.graph.external_filter::<SuBitGraph>(),
    )?;
    let closed_root = ClosedRootNormalization {
        external,
        external_routing,
        factor: Atom::num(1),
    };
    let auxiliary = auxiliary_mass();
    let insertions = term
        .counterterms
        .iter()
        .map(|child| {
            case.deltas
                .get(&child.argument.support)
                .cloned()
                .ok_or_else(|| eyre!("no Delta for a counterterm support"))
        })
        .collect::<Result<Vec<_>>>()?;
    let site = Site {
        graph: &case.graph,
        model,
        plan: &case.plan,
        term,
        children: &case.deltas,
        root_closure: Some(&closure),
        closed_root: Some(&closed_root),
        closure_data: None,
        root_projector: true,
        root_edges: &[],
        reference_mass_squared: *REFERENCE_MASS_SQUARED,
        auxiliary_mass_squared: &auxiliary,
        settings: &vakint_settings,
        vakint,
        uv_settings: &uv,
        top_mass: Some(TopMassTerms::Keep),
        epsilon: GS.dim_epsilon,
        today: "k5-site-test".into(),
        ir_hard: ir_hard.then_some(IrHardSite {
            correction: 0,
            term: 0,
            attempt: 0,
            hard: "test",
            degree: 1,
            depth: term.remainder.loops + 1,
            closure: None,
        }),
    };
    let run = || {
        uncontracted::attempt(&site, Kind::Surrogate, 7, Mode::On, || {
            RStarIntegrand::new_with_closure(
                &case.graph,
                model,
                &term.remainder,
                &insertions,
                case.plan.marked_edge,
                Atom::var(*REFERENCE_MASS_SQUARED),
                Some(&closure),
            )
        })
    };
    match context {
        Some(provenance) => k5::with_test_context(provenance, run),
        None => run(),
    }
}

/// Every check of a closed-root fixture: the record equals the fixture, its
/// value equals the fixture's contracted value, and so does the oracle value
/// of the record text (vakint's reader, the closure entry, finish_tensor).
fn assert_root_fixture(case: &RootCase, value_text: &str, vakint: &Vakint) -> Result<()> {
    let attempt = root_attempt(case, vakint, Some(case.provenance.clone()), false)?;
    let text = attempt
        .text
        .clone()
        .ok_or_else(|| eyre!("no record built: {:?}", attempt.outcome))?;
    let ours = record(&text)?;
    let differences = record_differences(&ours, &case.fixture);
    ensure!(
        differences.is_empty(),
        "record != fixture: {differences:#?}"
    );
    let value = attempt
        .value()
        .ok_or_else(|| eyre!("declined: {:?}", attempt.outcome))?;
    let expected = fixture_value(value_text)?;
    let order = valid_through(value_text)?;
    ensure!(
        equal_through(&value, &expected, order)?,
        "value != the fixture's contracted value through eps^{order}"
    );
    let oracle = oracle_value(vakint, &text, 7)?;
    ensure!(
        equal_through(&oracle, &expected, order)?,
        "oracle value != the fixture's contracted value"
    );
    Ok(())
}

/// The direct record of the L2 17 surrogate with the fixture's branch.
fn direct_case(vakint: &Vakint) -> Result<(Graph, feynkit_graph::FeynmanDiagram, Arc<Provenance>)> {
    let prepared = prepare(L2_17_DOT, vakint)?;
    let fixture = record(K5_DIRECT)?;
    let diagram = with_branch(&prepared.surrogate, &fixture)?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    Ok((graph, diagram, prepared.provenance))
}

/// `k5::direct` on `graph` inside a diagram ledger (when `ledger`), with
/// every attempt captured; the ledger's close result too.
fn direct_run(
    graph: &Graph,
    diagram: &feynkit_graph::FeynmanDiagram,
    provenance: Option<&Arc<Provenance>>,
    vakint: &Vakint,
    switch: Switch,
    ledger: bool,
) -> (Result<Atom>, Vec<Captured>, Option<Result<()>>) {
    let uv = settings();
    with_switch(switch, || {
        let opened = if ledger {
            DiagramLedger::open(&graph.name, "k5-test")
        } else {
            None
        };
        let (value, captured) = with_capture(|| {
            k5::direct(provenance, graph, diagram.model(), vakint, &uv, 2, || {
                super::integrate_propagator(graph, diagram.model(), vakint, &uv, 2)
            })
        });
        let closed = opened.map(|ledger| ledger.close(value.is_ok()).map(|_| ()));
        (value, captured, closed)
    })
}

/// A direct record's own fields (the kind, the provenance) and, through the
/// U7a comparator with the kinds blanked, all the others.
fn direct_differences(ours: &FactorizedNumerator, theirs: &FactorizedNumerator) -> Vec<String> {
    let mut out = Vec::new();
    if ours.kind != FactorizedKind::Direct || theirs.kind != FactorizedKind::Direct {
        out.push("kind is not direct on both sides".into());
    }
    let provenance = |r: &FactorizedNumerator| {
        r.record
            .as_ref()
            .map(|record| (record.original_cycles, record.depth))
    };
    if provenance(ours) != provenance(theirs) {
        out.push(format!(
            "direct provenance {:?} vs {:?}",
            provenance(ours),
            provenance(theirs)
        ));
    }
    let blank = ClosureSpec {
        edge: 0,
        mass_squared: Atom::Zero,
        auxiliary_mass_squared: Atom::Zero,
        power: 1,
        independent_externals: false,
        taylor: Some(TaylorSpec {
            degree: 0,
            scaled_edges: Vec::new(),
            compensate_mass: false,
            orders: TaylorOrders::All,
        }),
        softs: Vec::new(),
        original_cycles: None,
        depth: None,
        weights: Vec::new(),
    };
    let (mut a, mut b) = (ours.clone(), theirs.clone());
    a.kind = FactorizedKind::Closure(blank.clone());
    b.kind = FactorizedKind::Closure(blank);
    out.extend(record_differences(&a, &b));
    out
}

// ---------------------------------------------------------------------------
// The oracle: k5-direct, k5-surrogate, k5-stress
// ---------------------------------------------------------------------------

/// k5-direct (L2 17, <= 4 loops): the direct record of the surrogate equals
/// the fixture, and its value (vakint's ports entry, finish_tensor) equals
/// the fixture's contracted value through eps^-1; the hook in check mode
/// compares it with today's integrate_propagator (equal) and counts it in
/// the diagram's frame (site 1, expected 1, compared 1).
#[test]
fn k5_direct_record_and_value_equal_the_fixture() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let vakint = Vakint::new()?;
        let (graph, diagram, provenance) = direct_case(&vakint)?;
        let fixture = record(K5_DIRECT)?;
        let expected = fixture_value(K5_DIRECT_VALUE)?;
        let order = valid_through(K5_DIRECT_VALUE)?;
        ensure!(order == -1, "k5-direct is valid through eps^-1");
        uncontracted_ledger::take_lines();
        let (value, captured, closed) = direct_run(
            &graph,
            &diagram,
            Some(&provenance),
            &vakint,
            check_surrogate(),
            true,
        );
        let today = value?;
        closed.ok_or_else(|| eyre!("no ledger"))??;
        ensure!(captured.len() == 1, "one record: {captured:?}");
        let attempt = &captured[0];
        ensure!(attempt.kind == Kind::Surrogate, "kind surrogate");
        ensure!(attempt.equal == Some(true), "check: equal ({attempt:?})");
        let text = attempt.text.clone().ok_or_else(|| eyre!("no text"))?;
        let ours = record(&text)?;
        let differences = direct_differences(&ours, &fixture);
        ensure!(
            differences.is_empty(),
            "record != fixture: {differences:#?}"
        );
        let value = attempt.value.clone().ok_or_else(|| eyre!("no value"))?;
        ensure!(
            equal_through(&value, &expected, order)?,
            "ours != the fixture's contracted value through eps^-1:\n{value}\n{expected}"
        );
        ensure!(
            equal_through(&today, &expected, order)?,
            "today's value on the fixture branch != the fixture's (the test's branch)"
        );
        // the oracle: the text alone through vakint's reader and ports entry
        let mut vakint_settings = settings().vakint.true_settings();
        vakint_settings.number_of_terms_in_epsilon_expansion = 2;
        let oracle =
            evaluate_factorized_ports(&vakint, &vakint_settings, &ours, &ours.prefactor.total())
                .map_err(|error| eyre!("{error}"))?;
        let oracle = RStarIntegrand::finish_tensor(oracle.value, 0)?;
        ensure!(
            equal_through(&oracle, &expected, order)?,
            "oracle value != the fixture's"
        );
        let lines = uncontracted_ledger::take_lines();
        let disposition = lines
            .iter()
            .find(|line| {
                line.starts_with("UNCONTRACTED_DISPOSITION ") && line.contains(" kind=surrogate ")
            })
            .ok_or_else(|| eyre!("no surrogate disposition line: {lines:#?}"))?;
        ensure!(
            disposition.contains(" expected=1 ")
                && disposition.contains(" compared=1 ")
                && disposition.contains(" equal=1 ")
                && disposition.ends_with(" outcome=ok"),
            "{disposition}"
        );
        ensure!(
            lines
                .iter()
                .any(|line| line.starts_with("UNCONTRACTED_KIND_TOTAL ")
                    && line.contains(" kind=surrogate ")
                    && line.contains(" sites=1 ")
                    && line.contains(" compared=1 ")),
            "no kind total line: {lines:#?}"
        );
        // mode on: our value replaces today's, equal through eps^-1
        let (on, captured, closed) = direct_run(
            &graph,
            &diagram,
            Some(&provenance),
            &vakint,
            on_surrogate(),
            true,
        );
        closed.ok_or_else(|| eyre!("no ledger"))??;
        ensure!(captured.len() == 1 && captured[0].outcome == "uncontracted");
        ensure!(equal_through(&on?, &expected, order)?, "on != fixture");
        Ok(())
    })
}

/// k5-surrogate (L5 22807, closed surrogate root, seq 62: ct: [2,3,4,5,7,8]
/// and the `ins:` coefficient vertex 8): record and value equal the fixture.
#[test]
fn k5_surrogate_record_and_value_equal_the_fixture() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let vakint = Vakint::new()?;
        let prepared = prepare(L5_22807_DOT, &vakint)?;
        ensure!(
            prepared.provenance.original_cycles == (1, 0),
            "original cycles {:?}",
            prepared.provenance.original_cycles
        );
        let case = root_case(&prepared, K5_SURROGATE)?;
        assert_root_fixture(&case, K5_SURROGATE_VALUE, &vakint)
    })
}

/// k5-stress (L5 3793, skel1 x ins4: up to four `ins:` vertices, one ct:
/// child that contracts a hole in four of the five records): every record
/// and value equals its fixture; R2 holds (vakint's reader validates it).
#[test]
fn k5_stress_records_and_values_equal_the_fixtures() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let vakint = Vakint::new()?;
        let prepared = prepare(L5_3793_DOT, &vakint)?;
        let mut failures = Vec::new();
        let mut passed = 0;
        for (name, input, value) in K5_STRESS {
            let result = root_case(&prepared, input)
                .and_then(|case| assert_root_fixture(&case, value, &vakint));
            match result {
                Ok(()) => {
                    passed += 1;
                    eprintln!("K5_STRESS {name} PASS");
                }
                Err(error) => {
                    eprintln!("K5_STRESS {name} FAIL: {error:#}");
                    failures.push(format!("{name}: {error:#}"));
                }
            }
        }
        ensure!(failures.is_empty(), "{passed}/5 passed: {failures:#?}");
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// Designed declines and loud errors
// ---------------------------------------------------------------------------

/// The closed list: `surrogate_ir_hard` is designed; the K5 errors are not.
#[test]
fn k5_designed_codes_are_a_closed_list() {
    assert!(uncontracted::is_designed_code("surrogate_ir_hard"));
    assert!(uncontracted::is_designed_code("surrogate_split_chain"));
    for code in [
        "surrogate_chain",
        "surrogate_coefficient_mass",
        "surrogate_context_missing",
        "surrogate_hole_twice",
        "surrogate_insertion_eps",
        "ir_hard_guard:panic",
    ] {
        assert!(!uncontracted::is_designed_code(code), "{code}");
    }
}

/// A K3-shaped IR hard term inside a surrogate declines by design BEFORE its
/// record is built (no text), and is counted per kind: in check mode a
/// diagram whose surrogate records all decline fails (no exception named).
#[test]
fn k5_ir_hard_term_declines_before_any_record_and_is_counted() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let prepared = prepare(L5_22807_DOT, &vakint)?;
        let case = root_case(&prepared, K5_SURROGATE)?;
        let attempt = with_switch(check_surrogate(), || {
            root_attempt(&case, &vakint, Some(case.provenance.clone()), true)
        })?;
        ensure!(attempt.text.is_none(), "a record was built");
        match &attempt.outcome {
            Outcome::Declined(decline) => ensure!(
                decline.code == "surrogate_ir_hard" && decline.designed,
                "{decline:?}"
            ),
            other => return Err(eyre!("not declined: {other:?}")),
        }
        // counted per kind; all declined -> the diagram fails in check mode
        with_switch(check_surrogate(), || -> Result<()> {
            let ledger =
                DiagramLedger::open("k5-ir-hard", "k5-test").ok_or_else(|| eyre!("no ledger"))?;
            uncontracted_ledger::expect(Kind::Surrogate);
            uncontracted::report(&attempt, None);
            let closed = ledger.close(true);
            let error = closed
                .err()
                .ok_or_else(|| eyre!("an all-declined check passed"))?;
            ensure!(format!("{error}").contains("declined"), "{error}");
            Ok(())
        })
    })
}

/// A mass symbol inside an `ins:` coefficient is an ERROR (never the
/// designed `insertion_mass_symbol`, which stays a ct: child's decline).
#[test]
fn k5_mass_in_a_coefficient_is_an_error() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let prepared = prepare(L5_22807_DOT, &vakint)?;
        let mut case = root_case(&prepared, K5_SURROGATE)?;
        let node = case
            .graph
            .underlying
            .iter_nodes()
            .map(|(node, _, _)| node)
            .find(|node| usize::from(*node) == 8)
            .ok_or_else(|| eyre!("no vertex 8"))?;
        let massive = case.graph[node].num.value.clone() * Atom::var(*REFERENCE_MASS_SQUARED);
        case.graph.underlying[node].num.value = massive;
        let attempt = root_attempt(&case, &vakint, Some(case.provenance.clone()), false)?;
        match &attempt.outcome {
            Outcome::Error(decline) => ensure!(
                decline.code == "surrogate_coefficient_mass" && !decline.designed,
                "{decline:?}"
            ),
            other => return Err(eyre!("not an error: {other:?}")),
        }
        Ok(())
    })
}

/// A chain that is not the scalar bubble chain (a bubble edge with a factor
/// other than 1) and a surrogate record outside the surrogate context are
/// errors, not declines.
#[test]
fn k5_broken_chain_and_missing_context_are_errors() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let prepared = prepare(L5_22807_DOT, &vakint)?;
        let mut case = root_case(&prepared, K5_SURROGATE)?;
        let missing = root_attempt(&case, &vakint, None, false)?;
        match &missing.outcome {
            Outcome::Error(decline) => {
                ensure!(decline.code == "surrogate_context_missing", "{decline:?}")
            }
            other => return Err(eyre!("not an error: {other:?}")),
        }
        let bubble = case
            .graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .find(|edge| usize::from(*edge) == 14)
            .ok_or_else(|| eyre!("no edge 14"))?;
        case.graph.underlying[bubble].num.value = Atom::num(2);
        let broken = root_attempt(&case, &vakint, Some(case.provenance.clone()), false)?;
        match &broken.outcome {
            Outcome::Error(decline) => {
                ensure!(decline.code == "surrogate_chain", "{decline:?}")
            }
            other => return Err(eyre!("not an error: {other:?}")),
        }
        Ok(())
    })
}

/// An undesigned failure is loud at the direct hook: in mode on and in check
/// mode the call fails (no fallback to today's value), and the diagram's
/// line says error.
#[test]
fn k5_undesigned_failure_is_loud_at_the_hook() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let (mut graph, diagram, provenance) = direct_case(&vakint)?;
        let bubble = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .find(|edge| usize::from(*edge) == 5)
            .ok_or_else(|| eyre!("no edge 5"))?;
        graph.underlying[bubble].num.value = Atom::num(2);
        for switch in [on_surrogate(), check_surrogate()] {
            let (value, captured, closed) =
                direct_run(&graph, &diagram, Some(&provenance), &vakint, switch, true);
            let error = value
                .err()
                .ok_or_else(|| eyre!("{:?}: an undesigned failure passed", switch.mode))?;
            ensure!(
                format!("{error:#}").contains("outcome=error:surrogate_chain"),
                "{error:#}"
            );
            ensure!(
                captured.len() == 1 && captured[0].outcome == "error:surrogate_chain",
                "{captured:?}"
            );
            let _ = closed;
        }
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// Zero comparisons, the frame, the run's coverage
// ---------------------------------------------------------------------------

/// A surrogate site whose records never reached a hook (here: none at all)
/// fails its diagram in check mode and in mode on.
#[test]
fn k5_site_without_records_fails_the_diagram() -> Result<()> {
    for switch in [check_surrogate(), on_surrogate()] {
        with_switch(switch, || -> Result<()> {
            let ledger =
                DiagramLedger::open("k5-site", "k5-test").ok_or_else(|| eyre!("no ledger"))?;
            uncontracted_ledger::k5_site(Kind::Surrogate, "k5-site")?;
            let error = ledger
                .close(true)
                .err()
                .ok_or_else(|| eyre!("{:?}: a site without records passed", switch.mode))?;
            ensure!(
                format!("{error}").contains("insertion-surrogate evaluations reached no record"),
                "{error}"
            );
            Ok(())
        })?;
    }
    Ok(())
}

/// A surrogate evaluated outside every diagram ledger is an error (never an
/// implicit ledger): at the site, at the direct hook and at the R* scope;
/// inside one the frame is there.
#[test]
fn k5_no_frame_is_an_error() -> Result<()> {
    for switch in [check_surrogate(), on_surrogate()] {
        no_frame_is_an_error(switch)?;
    }
    Ok(())
}

fn no_frame_is_an_error(switch: Switch) -> Result<()> {
    let uv = settings();
    with_switch(switch, || -> Result<()> {
        ensure!(!uncontracted_ledger::in_diagram_ledger());
        ensure!(uncontracted_ledger::k5_site(Kind::Surrogate, "nowhere").is_err());
        let provenance = Arc::new(Provenance {
            original: "nowhere".into(),
            original_cycles: (0, 0),
            holes: Vec::new(),
        });
        let scoped = k5::scoped(Some(&provenance), "nowhere", true, &uv, || Ok(()));
        ensure!(
            scoped
                .err()
                .is_some_and(|e| format!("{e}").contains("outside every diagram ledger")),
            "scoped"
        );
        let ledger =
            DiagramLedger::open("k5-frame", "k5-test").ok_or_else(|| eyre!("no ledger"))?;
        ensure!(uncontracted_ledger::in_diagram_ledger());
        let inside = k5::scoped(Some(&provenance), "k5-frame", true, &uv, || {
            ensure!(
                uncontracted_ledger::in_diagram_ledger(),
                "no frame at the hook"
            );
            // the site's one record
            uncontracted_ledger::expect(Kind::Surrogate);
            Ok(())
        });
        drop(ledger);
        inside
    })
}

/// The direct hook outside every diagram ledger fails (with the kind
/// enabled), before any record.
#[test]
fn k5_direct_without_a_ledger_fails() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let (graph, diagram, provenance) = direct_case(&vakint)?;
        for switch in [check_surrogate(), on_surrogate()] {
            let (value, captured, _) =
                direct_run(&graph, &diagram, Some(&provenance), &vakint, switch, false);
            ensure!(
                value
                    .err()
                    .is_some_and(|e| format!("{e}").contains("outside every diagram ledger")),
                "{:?}: the direct hook ran outside a ledger",
                switch.mode
            );
            ensure!(captured.is_empty(), "a record was built: {captured:?}");
        }
        Ok(())
    })
}

/// The run's coverage: a check run in which the enabled kind compared
/// nothing fails; one comparison makes it pass; compared != equal fails; the
/// field names are the beta's (kind= compared= hit_compared= declined=).
#[test]
fn k5_run_with_zero_comparisons_fails() -> Result<()> {
    with_switch(check_surrogate(), || -> Result<()> {
        uncontracted_ledger::reset_totals();
        let ledger =
            DiagramLedger::open("k5-run-1", "k5-test").ok_or_else(|| eyre!("no ledger"))?;
        ledger.close(true)?;
        let error = uncontracted_ledger::run_coverage()
            .err()
            .ok_or_else(|| eyre!("a run that compared nothing passed"))?;
        ensure!(
            format!("{error}").contains("kind=surrogate: zero_comparisons"),
            "{error}"
        );
        let ledger =
            DiagramLedger::open("k5-run-2", "k5-test").ok_or_else(|| eyre!("no ledger"))?;
        uncontracted_ledger::expect(Kind::Surrogate);
        uncontracted_ledger::record_value(Kind::Surrogate, Mode::Check);
        ledger.close(true)?;
        let lines = uncontracted_ledger::run_coverage()?;
        ensure!(
            lines.len() == 1
                && lines[0].contains(" kind=surrogate ")
                && lines[0].contains(" compared=1 ")
                && lines[0].contains(" hit_compared=0 ")
                && lines[0].contains(" declined=0 ")
                && lines[0].ends_with(" outcome=ok"),
            "{lines:?}"
        );
        let ledger =
            DiagramLedger::open("k5-run-3", "k5-test").ok_or_else(|| eyre!("no ledger"))?;
        uncontracted_ledger::expect(Kind::Surrogate);
        uncontracted_ledger::record_compared_unequal(Kind::Surrogate);
        let _ = ledger.close(true);
        let error = uncontracted_ledger::run_coverage()
            .err()
            .ok_or_else(|| eyre!("compared != equal passed"))?;
        ensure!(
            format!("{error}").contains("kind=surrogate: unequal"),
            "{error}"
        );
        uncontracted_ledger::reset_totals();
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// The switch, the cache key, the off mode
// ---------------------------------------------------------------------------

/// `surrogate` is built but not in the unlisted default of on / check; it
/// must be listed.
#[test]
fn k5_surrogate_must_be_listed() {
    for value in ["on", "check", "1"] {
        let switch = Switch::parse(Some(value)).unwrap();
        assert!(!switch.enabled(Kind::Surrogate), "{value}");
        assert!(
            switch.enabled(Kind::Root) && switch.enabled(Kind::Child),
            "{value}"
        );
    }
    for value in ["on:surrogate", "check:root,child,surrogate", "check:e"] {
        let switch = Switch::parse(Some(value)).unwrap();
        assert!(switch.enabled(Kind::Surrogate), "{value}");
        assert!(switch.inert_kinds().is_empty(), "{value}");
    }
    assert!(Kind::Surrogate.implemented() && !Kind::Surrogate.in_default());
}

/// With `surrogate` enabled, a child looked up OUTSIDE the surrogate's R*
/// (its kind not enabled) gets the inert key, one inside does not: an inert
/// fill never serves a surrogate lookup.
#[test]
fn k5_inert_cache_key_with_surrogate_enabled() {
    with_switch(on_surrogate(), || {
        assert_eq!(
            uncontracted::cache_key_suffix().as_deref(),
            Some("#uncontracted=surrogate#uncontracted-inert=child")
        );
        in_dressed_scope(|| {
            assert_eq!(
                uncontracted::cache_key_suffix().as_deref(),
                Some("#uncontracted=surrogate")
            );
        });
    });
    with_switch(check_surrogate(), || {
        assert_eq!(
            uncontracted::cache_key_suffix().as_deref(),
            Some("#uncontracted-inert=child")
        );
    });
}

/// Off: the direct hook is today's call, value for value, and no provenance
/// is computed; the whole evaluator off equals the evaluator without the
/// switch (L2 17).
#[test]
fn k5_off_is_todays_path() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let original = colour_factored(&load(L2_17_DOT)?)?;
        let evaluator = PropagatorInsertionEvaluator::new(settings());
        let off = with_switch(Switch::OFF, || {
            evaluator.pole_part_closing_root(&original, &vakint, true)
        })?
        .ok_or_else(|| eyre!("no insertion"))?;
        let prepared = evaluator
            .prepare(&original, &vakint)?
            .ok_or_else(|| eyre!("no insertion"))?;
        // the pre-K5 computation, called directly: the surrogate's direct
        // propagator integral and the owner map (insertions.rs before U7c)
        let surrogate = prepared.surrogate().to_gamma_loop_graph(None, true)?;
        let direct = super::integrate_propagator(
            &surrogate,
            prepared.surrogate().model(),
            &vakint,
            &settings(),
            prepared.loop_count() as i64,
        )?;
        ensure!(
            off.poles == prepared.restore_atom(&direct)?,
            "off != today's integrate_propagator + restore_atom"
        );
        let none = with_switch(Switch::OFF, || {
            k5::provenance(
                &original,
                prepared
                    .holes()
                    .iter()
                    .map(|hole| (hole.vertex.0, hole.loops, &hole.original.support)),
            )
        })?;
        ensure!(none.is_none(), "a provenance computed with the switch off");
        Ok(())
    })
}

/// The whole evaluator on L2 17 (the direct record) in check and on mode,
/// inside a diagram ledger: compared and equal, on == off through eps^-1.
#[test]
fn k5_l2_17_check_and_on_equal_off() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let original = colour_factored(&load(L2_17_DOT)?)?;
        whole_graph_modes(&original, &vakint, 1)
    })
}

/// The whole evaluator (`pole_part_closing_root`) in off, check and on mode
/// inside a diagram ledger: every surrogate record compared and equal in
/// check (at least `min_records`), the diagrams' ledgers close ok, on == off
/// through eps^-1 of the bare poles.
fn whole_graph_modes(
    original: &feynkit_graph::FeynmanDiagram,
    vakint: &Vakint,
    min_records: usize,
) -> Result<()> {
    whole_graph_summary(original, vakint, min_records).map(|_| ())
}

/// [`whole_graph_modes`], returning the check run's surrogate outcomes by
/// name (`compared`, `declined:<code>`, ...).
fn whole_graph_summary(
    original: &feynkit_graph::FeynmanDiagram,
    vakint: &Vakint,
    min_records: usize,
) -> Result<BTreeMap<String, usize>> {
    let evaluator = PropagatorInsertionEvaluator::new(settings());
    let run = |switch: Switch| -> Result<(Atom, Vec<Captured>)> {
        with_switch(switch, || {
            let ledger = DiagramLedger::open(original.name(), "k5-whole");
            let (value, captured) =
                with_capture(|| evaluator.pole_part_closing_root(original, vakint, true));
            if let Some(ledger) = ledger {
                ledger.close(value.is_ok())?;
            }
            let value = value?.ok_or_else(|| eyre!("no insertion"))?;
            Ok((value.poles, captured))
        })
    };
    let (off, _) = run(Switch::OFF)?;
    let (check, captured) = run(check_surrogate())?;
    ensure!(off == check, "check mode must return today's value");
    let surrogate = captured
        .iter()
        .filter(|attempt| attempt.kind == Kind::Surrogate)
        .collect::<Vec<_>>();
    ensure!(
        surrogate.len() >= min_records,
        "{} surrogate records (at least {min_records} expected)",
        surrogate.len()
    );
    let mut summary = BTreeMap::<String, usize>::new();
    for attempt in &surrogate {
        ensure!(
            attempt.equal == Some(true)
                || attempt.outcome.starts_with("declined:")
                || attempt.outcome == "zero",
            "{} id {}: {} equal {:?}",
            attempt.graph,
            attempt.id,
            attempt.outcome,
            attempt.equal
        );
        let name = if attempt.equal == Some(true) {
            "compared".to_owned()
        } else {
            attempt.outcome.clone()
        };
        *summary.entry(name).or_default() += 1;
    }
    ensure!(
        summary.get("compared").copied().unwrap_or(0) > 0,
        "no surrogate record compared: {summary:?}"
    );
    let (on, captured) = run(on_surrogate())?;
    ensure!(
        captured
            .iter()
            .filter(|attempt| attempt.kind == Kind::Surrogate)
            .all(|attempt| !attempt.outcome.starts_with("error:")),
        "an error in mode on"
    );
    ensure!(
        equal_through(&on, &off, -1)? || atoms_equal(&on, &off),
        "on != off through eps^-1"
    );
    eprintln!(
        "K5_WHOLE graph={} surrogate_records={} {summary:?}",
        original.name(),
        surrogate.len()
    );
    Ok(summary)
}

/// L5 22807 (k5-surrogate's graph: the closed surrogate root, 9 remainders,
/// UV children inside the surrogate) through the whole evaluator.
#[test]
#[ignore = "about a minute (L5)"]
fn k5_l5_22807_check_and_on_equal_off() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let original = colour_factored(&load(L5_22807_DOT)?)?;
        whole_graph_modes(&original, &vakint, 2)
    })
}

/// L5 3793 (k5-stress's graph: four one-loop insertions) through the whole
/// evaluator.
#[test]
#[ignore = "about a minute (L5)"]
fn k5_l5_3793_check_and_on_equal_off() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let original = colour_factored(&load(L5_3793_DOT)?)?;
        whole_graph_modes(&original, &vakint, 2)
    })
}

/// The L5 surrogates whose holes have two or more loops (143019: two,
/// 143190: three, 163502 / 175270 / 139207: four) through the whole
/// evaluator in off, check and on mode: every record whose chain a UV child
/// boundary cuts is STATED (the anchor carries the hole; the far end's slot
/// a spectator port when the record lacks it) and compared equal, none
/// declined `surrogate_split_chain` (measured: no record of these five needs
/// it); on == off through eps^-1.
#[test]
fn k5_l5_multi_loop_holes_check_and_on_equal_off() -> Result<()> {
    multi_loop_holes(k5::HoleCut::Form)
}

/// The same five graphs with the DEFAULT hole-cut mode (`decline`: the
/// single-hedge `ins:` stays out of the G-D gate): every kept and cut hole is
/// declined by design (`surrogate_hole_cut`, today's path), counted per kind,
/// every other record compared and equal, on == off through eps^-1.
#[test]
fn k5_l5_multi_loop_holes_with_the_hole_cut_declined() -> Result<()> {
    multi_loop_holes(k5::HoleCut::Decline)
}

fn multi_loop_holes(mode: k5::HoleCut) -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let mut failures = Vec::new();
        let mut hole_cut_declines = 0;
        for (name, dot) in L5_MULTI_LOOP_HOLES {
            let original = colour_factored(&load(dot)?)?;
            let prepared = prepare(dot, &vakint)?;
            let holes = prepared
                .provenance
                .holes
                .iter()
                .map(|hole| hole.loops)
                .collect::<Vec<_>>();
            match k5::with_hole_cut(mode, || whole_graph_summary(&original, &vakint, 1)) {
                Ok(summary) if summary.contains_key("declined:surrogate_split_chain") => {
                    eprintln!(
                        "K5_MULTI {name} FAIL hole_loops={holes:?}: split chains declined {summary:?}"
                    );
                    failures.push(format!(
                        "{name}: a cut chain was declined, not stated: {summary:?}"
                    ));
                }
                Ok(summary)
                    if mode == k5::HoleCut::Form
                        && summary.contains_key("declined:surrogate_hole_cut") =>
                {
                    failures.push(format!(
                        "{name}: hole_cut declined in form mode: {summary:?}"
                    ));
                }
                Ok(summary) => {
                    hole_cut_declines += summary
                        .get("declined:surrogate_hole_cut")
                        .copied()
                        .unwrap_or(0);
                    eprintln!(
                        "K5_MULTI mode={mode:?} {name} PASS surrogate_loops={} hole_loops={holes:?} {summary:?}",
                        prepared.surrogate.loop_count()
                    )
                }
                Err(error) => {
                    eprintln!("K5_MULTI {name} FAIL hole_loops={holes:?}: {error:#}");
                    failures.push(format!("{name}: {error:#}"));
                }
            }
        }
        ensure!(failures.is_empty(), "{failures:#?}");
        if mode == k5::HoleCut::Decline {
            ensure!(
                hole_cut_declines > 0,
                "decline mode declined no kept and cut hole on the five graphs"
            );
        }
        Ok(())
    })
}

/// Per site: a surrogate whose R* adds no expected record of the kind (its
/// records missed the hook) fails at once, in check and in on, even when
/// another site of the same diagram recorded.
#[test]
fn k5_site_without_records_fails_per_site() -> Result<()> {
    let uv = settings();
    let provenance = Arc::new(Provenance {
        original: "k5-sites".into(),
        original_cycles: (0, 0),
        holes: Vec::new(),
    });
    for switch in [check_surrogate(), on_surrogate()] {
        with_switch(switch, || -> Result<()> {
            let ledger =
                DiagramLedger::open("k5-sites", "k5-test").ok_or_else(|| eyre!("no ledger"))?;
            k5::scoped(Some(&provenance), "k5-sites", true, &uv, || {
                uncontracted_ledger::expect(Kind::Surrogate);
                Ok(())
            })?;
            let missed = k5::scoped(Some(&provenance), "k5-sites", true, &uv, || Ok(()));
            ensure!(
                missed
                    .err()
                    .is_some_and(|e| format!("{e}").contains("surrogate_site_without_records")),
                "{:?}: a site without records passed beside one with",
                switch.mode
            );
            drop(ledger);
            Ok(())
        })?;
    }
    Ok(())
}

/// Settings under which every surrogate record would decline by design are
/// errors with the kind enabled (check and on): an open surrogate root, an
/// additional normalisation other than 1 (both at the R* scope and at the
/// direct hook); and the evaluator's open-root entry (`pole_part`, the
/// surrogate root not closed) on L5 3793 fails instead of declining every
/// root record.
#[test]
fn k5_settings_wide_declines_are_errors() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let provenance = Arc::new(Provenance {
            original: "k5-settings".into(),
            original_cycles: (0, 0),
            holes: Vec::new(),
        });
        let mut odd = settings();
        odd.vakint.additional_normalization = "2".into();
        let (graph, diagram, direct_provenance) = direct_case(&vakint)?;
        let original = colour_factored(&load(L5_3793_DOT)?)?;
        let evaluator = PropagatorInsertionEvaluator::new(settings());
        for switch in [check_surrogate(), on_surrogate()] {
            with_switch(switch, || -> Result<()> {
                let ledger = DiagramLedger::open("k5-settings", "k5-test")
                    .ok_or_else(|| eyre!("no ledger"))?;
                let open = k5::scoped(Some(&provenance), "k5-settings", false, &settings(), || {
                    Ok(())
                });
                ensure!(
                    open.err()
                        .is_some_and(|e| format!("{e}").contains("surrogate_settings")),
                    "open root"
                );
                let normalised =
                    k5::scoped(Some(&provenance), "k5-settings", true, &odd, || Ok(()));
                ensure!(
                    normalised
                        .err()
                        .is_some_and(|e| format!("{e}").contains("additional_normalization")),
                    "additional normalisation (R*)"
                );
                let direct = k5::direct(
                    Some(&direct_provenance),
                    &graph,
                    diagram.model(),
                    &vakint,
                    &odd,
                    2,
                    || Ok(Atom::Zero),
                );
                ensure!(
                    direct
                        .err()
                        .is_some_and(|e| format!("{e}").contains("additional_normalization")),
                    "additional normalisation (direct)"
                );
                drop(ledger);
                let ledger = DiagramLedger::open(original.name(), "k5-test")
                    .ok_or_else(|| eyre!("no ledger"))?;
                let pole_part = evaluator.pole_part(&original, &vakint);
                ensure!(
                    pole_part
                        .err()
                        .is_some_and(|e| format!("{e:#}").contains("surrogate_settings")),
                    "{:?}: the open surrogate root of L5 3793 did not fail",
                    switch.mode
                );
                drop(ledger);
                Ok(())
            })?;
        }
        Ok(())
    })
}

/// Correction 4: an insertion whose coefficient depends on eps other than
/// through sums, products and integer powers is an error (here log(eps) on
/// k5-surrogate's coefficient vertex); a direct record's coefficient vertex
/// with an eps pole is an error (its entry raises no depth).
#[test]
fn k5_insertion_eps_rules_are_errors() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let prepared = prepare(L5_22807_DOT, &vakint)?;
        let mut case = root_case(&prepared, K5_SURROGATE)?;
        let node = case
            .graph
            .underlying
            .iter_nodes()
            .map(|(node, _, _)| node)
            .find(|node| usize::from(*node) == 8)
            .ok_or_else(|| eyre!("no vertex 8"))?;
        let eps = Atom::var(GS.dim_epsilon);
        let logged = case.graph[node].num.value.clone()
            * symbolica::function!(symbolica::symbol!("log"), eps.clone());
        case.graph.underlying[node].num.value = logged;
        let attempt = root_attempt(&case, &vakint, Some(case.provenance.clone()), false)?;
        match &attempt.outcome {
            Outcome::Error(decline) => ensure!(
                decline.code == "surrogate_insertion_eps" && !decline.designed,
                "{decline:?}"
            ),
            other => return Err(eyre!("log(eps): not an error: {other:?}")),
        }
        let (mut graph, diagram, provenance) = direct_case(&vakint)?;
        let node = graph
            .underlying
            .iter_nodes()
            .map(|(node, _, _)| node)
            .find(|node| usize::from(*node) == 2)
            .ok_or_else(|| eyre!("no vertex 2"))?;
        let pole = graph[node].num.value.clone() / &eps;
        graph.underlying[node].num.value = pole;
        let (value, captured, _) = direct_run(
            &graph,
            &diagram,
            Some(&provenance),
            &vakint,
            on_surrogate(),
            true,
        );
        ensure!(
            value.is_err()
                && captured.len() == 1
                && captured[0].outcome == "error:surrogate_insertion_eps",
            "direct eps pole: {captured:?}"
        );
        Ok(())
    })
}

/// The record EMITTED from this build's own prepared surrogate of L2 17
/// (no factor taken from the fixture): its structure equals the fixture's
/// (edges, momenta, lines, vertices' rules, ports, the insertion's hedges and
/// cycles, the provenance), and its insertion factor (the K6 normalised
/// insertion this build's `prepare` computes) equals the fixture's up to the
/// colour branch's overall sign; check mode compares it with today's.
#[test]
fn k5_direct_emitted_from_the_prepared_surrogate() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let prepared = prepare(L2_17_DOT, &vakint)?;
        let graph = prepared.surrogate.to_gamma_loop_graph(None, true)?;
        let fixture = record(K5_DIRECT)?;
        let (value, captured, closed) = direct_run(
            &graph,
            &prepared.surrogate,
            Some(&prepared.provenance),
            &vakint,
            check_surrogate(),
            true,
        );
        let _today = value?;
        closed.ok_or_else(|| eyre!("no ledger"))??;
        ensure!(
            captured.len() == 1 && captured[0].equal == Some(true),
            "{captured:?}"
        );
        let ours = record(
            captured[0]
                .text
                .as_deref()
                .ok_or_else(|| eyre!("no text"))?,
        )?;
        let structure = direct_differences(&ours, &fixture)
            .into_iter()
            .filter(|difference| {
                !difference.contains(": factor")
                    && !difference.starts_with("vertex ")
                    && !difference.contains("coefficient x factor")
            })
            .collect::<Vec<_>>();
        ensure!(structure.is_empty(), "structure: {structure:#?}");
        let product = |r: &FactorizedNumerator| {
            r.insertions
                .iter()
                .map(|i| &i.coefficient * &i.factor)
                .fold(Atom::num(1), |a, b| a * b)
        };
        let (a, b) = (product(&ours), product(&fixture));
        ensure!(
            atoms_equal(&a, &b) || atoms_equal(&a, &(-b.clone())),
            "the prepared coefficient differs from the fixture's beyond a sign:\n{a}\n{b}"
        );
        eprintln!(
            "K5_DIRECT_EMITTED coefficient {}",
            if atoms_equal(&a, &b) {
                "equal"
            } else {
                "opposite sign (colour branch)"
            }
        );
        Ok(())
    })
}

/// The designed `surrogate_split_chain` (closed list) fires BEFORE any
/// record: a ct: child whose Delta carries the slot of a half-edge the
/// record has but the ct: vertex does not own (here k5-surrogate's site with
/// a metric of the far-end half-edge 22 of the hole at 8 injected into the
/// ct: child's Delta) is declined by design with no record text built.
#[test]
fn k5_split_chain_declines_before_any_record() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assert!(uncontracted::is_designed_code("surrogate_split_chain"));
        let vakint = Vakint::new()?;
        let prepared = prepare(L5_22807_DOT, &vakint)?;
        let mut case = root_case(&prepared, K5_SURROGATE)?;
        // edge 12 (vertex 9 -> 6, half-edges 22 / 23) is a remaining edge of
        // the record: 22 (at the chain end 9) is a half-edge of the record;
        // g(22, 23) in the ct: child's Delta (g(a, a) would be traced at once)
        ensure!(record(K5_SURROGATE)?.edges.iter().any(|edge| edge.id == 12));
        let far = symbolica::try_parse!(
            "spenso::{symmetric,linear,real,spenso::tensor}::g(spenso::{spenso::representation,spenso::self_dual}::mink(gammalooprs::{}::dim,gammalooprs::{spenso::index}::hedge(22,1)),spenso::{spenso::representation,spenso::self_dual}::mink(gammalooprs::{}::dim,gammalooprs::{spenso::index}::hedge(23,1)))"
        )
        .map_err(|error| eyre!("{error}"))?;
        for delta in case.deltas.values_mut() {
            *delta = &*delta * &far;
        }
        let attempt = root_attempt(&case, &vakint, Some(case.provenance.clone()), false)?;
        ensure!(attempt.text.is_none(), "a record was built");
        match &attempt.outcome {
            Outcome::Declined(decline) => ensure!(
                decline.code == "surrogate_split_chain" && decline.designed,
                "{decline:?}"
            ),
            other => return Err(eyre!("not declined: {other:?}")),
        }
        Ok(())
    })
}

/// The free hedge slots `hedge(h, k)` (Minkowski / bispinor) of an atom, as
/// canonical strings.
fn free_hedge_slots(atom: &Atom) -> BTreeSet<String> {
    let mut found = BTreeSet::new();
    atom.as_view().visitor(&mut |view| {
        if let symbolica::atom::AtomView::Fun(function) = view {
            let symbol = function.get_symbol();
            let name = symbol.get_name();
            if (name.ends_with("::mink") || name.ends_with("::bis"))
                && function.get_nargs() == 2
                && function
                    .iter()
                    .nth(1)
                    .is_some_and(|index| index.to_canonical_string().contains("::hedge("))
            {
                found.insert(view.to_canonical_string());
                return false;
            }
        }
        true
    });
    found
}

/// Spectator ports (beta's question (ii)): on the multi-loop-hole graphs the
/// records whose ports include a slot of a half-edge that is not a crown leg
/// (the coefficient's far-end slot, or a far edge's slot without the anchor)
/// return a value with exactly the free index slots of today's value for the
/// same remainder (the old path's child Delta), and the two values are equal
/// in check mode.
#[test]
fn k5_spectator_ports_have_todays_index_structure() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let evaluator = PropagatorInsertionEvaluator::new(settings());
        let mut spectators = 0;
        for (name, dot) in L5_MULTI_LOOP_HOLES {
            let original = colour_factored(&load(dot)?)?;
            let before = spectators;
            let captured = with_switch(check_surrogate(), || -> Result<Vec<Captured>> {
                let ledger = DiagramLedger::open(original.name(), "k5-spectator");
                // the single-hedge form (its spectator ports) is the explicit form mode
                let (value, captured) = k5::with_hole_cut(k5::HoleCut::Form, || {
                    with_capture(|| evaluator.pole_part_closing_root(&original, &vakint, true))
                });
                if let Some(ledger) = ledger {
                    ledger.close(value.is_ok())?;
                }
                value?;
                Ok(captured)
            })?;
            for attempt in captured.iter().filter(|a| a.kind == Kind::Surrogate) {
                let Some(text) = &attempt.text else { continue };
                let ours = record(text)?;
                let legs = ours
                    .edges
                    .iter()
                    .filter(|edge| edge.is_external())
                    .flat_map(|edge| [edge.source_hedge, edge.target_hedge])
                    .flatten()
                    .collect::<BTreeSet<_>>();
                let spectator = ours.ports.iter().any(|port| {
                    let text = port.to_canonical_string();
                    !legs.iter().any(|h| text.contains(&format!("::hedge({h},")))
                });
                if !spectator {
                    continue;
                }
                spectators += 1;
                let value = attempt.value.clone().ok_or_else(|| eyre!("no value"))?;
                let today = attempt
                    .contracted
                    .clone()
                    .ok_or_else(|| eyre!("not compared"))?;
                ensure!(attempt.equal == Some(true), "id {}: not equal", attempt.id);
                let (a, b) = (free_hedge_slots(&value), free_hedge_slots(&today));
                ensure!(
                    a == b,
                    "id {}: the record's value has the free slots {a:?}, today's {b:?}",
                    attempt.id
                );
                eprintln!(
                    "K5_SPECTATOR id={} ports={} free_slots={a:?} == today's",
                    attempt.id,
                    ours.ports.len()
                );
            }
            eprintln!("K5_SPECTATOR graph={name} records={}", spectators - before);
        }
        ensure!(
            spectators > 0,
            "no record with a spectator port on the five graphs"
        );
        eprintln!("K5_SPECTATOR records={spectators}");
        Ok(())
    })
}

/// The hole-cut switch: `decline` is the default (unset environment), and
/// the reader probe of this build's vakint (b6bb69835) accepts the
/// single-hedge form (the gated reader ccba5e0c8a does not: `auto` would
/// decline there); `surrogate_hole_cut` is in the closed list.
#[test]
fn k5_hole_cut_switch_defaults_to_decline() {
    assert!(std::env::var_os(k5::HOLE_CUT_VARIABLE).is_none());
    assert_eq!(k5::hole_cut(), k5::HoleCut::Decline);
    assert!(k5::reader_accepts_single_hedge());
    assert!(uncontracted::is_designed_code("surrogate_hole_cut"));
    assert_eq!(
        k5::with_hole_cut(k5::HoleCut::Form, k5::hole_cut),
        k5::HoleCut::Form
    );
}

/// Beta's point 6: the coefficient C_H carries the momentum of the excised
/// diagram's first external edge, which may be the chain's far edge e_out:
/// k5-surrogate's coefficient vertex 8 rewritten in Q(12) (= +Q(10): edge 10
/// enters the chain at 8, edge 12 leaves it at 9) is emitted in Q(10), the
/// anchor's own edge, and the record and value still equal the fixture.
#[test]
fn k5_coefficient_momentum_is_written_in_the_anchor_edge() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let prepared = prepare(L5_22807_DOT, &vakint)?;
        let mut case = root_case(&prepared, K5_SURROGATE)?;
        let node = case
            .graph
            .underlying
            .iter_nodes()
            .map(|(node, _, _)| node)
            .find(|node| usize::from(*node) == 8)
            .ok_or_else(|| eyre!("no vertex 8"))?;
        let q = |e: usize| symbolica::function!(GS.emr_mom, e, crate::utils::W_.x___);
        let far = case.graph[node]
            .num
            .value
            .replace(q(10).to_pattern())
            .with(q(12).to_pattern());
        ensure!(
            far.to_canonical_string().contains("::Q(12,")
                && !far.to_canonical_string().contains("::Q(10,"),
            "the rewrite to Q(12) did not take"
        );
        case.graph.underlying[node].num.value = far;
        let attempt = root_attempt(&case, &vakint, Some(case.provenance.clone()), false)?;
        let text = attempt
            .text
            .clone()
            .ok_or_else(|| eyre!("no record: {:?}", attempt.outcome))?;
        let ours = record(&text)?;
        let factor = ours
            .insertions
            .iter()
            .find(|insertion| insertion.vertex == 8)
            .ok_or_else(|| eyre!("no ins: at 8"))?
            .factor
            .to_canonical_string();
        ensure!(
            factor.contains("::Q(10,") && !factor.contains("::Q(12,"),
            "the coefficient is not written in Q(10)"
        );
        let differences = record_differences(&ours, &case.fixture);
        ensure!(differences.is_empty(), "{differences:#?}");
        let value = attempt.value().ok_or_else(|| eyre!("no value"))?;
        ensure!(equal_through(
            &value,
            &fixture_value(K5_SURROGATE_VALUE)?,
            valid_through(K5_SURROGATE_VALUE)?
        )?);
        Ok(())
    })
}

/// Beta's point 7: a cut child that is a run of whole bubbles through the
/// anchor (L5 143019's child {9,10}: superficial degree 0) carries the
/// Taylor degree of its local dod, here deg C_H = 2 (the ghost coefficient's
/// Q(e_in)^2), as today's outer Taylor step does.
#[test]
fn k5_cut_child_taylor_degree_is_the_local_dod() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let (_, dot) = L5_MULTI_LOOP_HOLES[0];
        let original = colour_factored(&load(dot)?)?;
        let evaluator = PropagatorInsertionEvaluator::new(settings());
        let captured = with_switch(check_surrogate(), || -> Result<Vec<Captured>> {
            let ledger = DiagramLedger::open(original.name(), "k5-taylor");
            let (value, captured) = k5::with_hole_cut(k5::HoleCut::Form, || {
                with_capture(|| evaluator.pole_part_closing_root(&original, &vakint, true))
            });
            if let Some(ledger) = ledger {
                ledger.close(value.is_ok())?;
            }
            value?;
            Ok(captured)
        })?;
        let child = captured
            .iter()
            .filter_map(|attempt| attempt.text.as_deref())
            .filter(|text| text.contains(" scope=surrogate:9,10 "))
            .map(record)
            .collect::<Result<Vec<_>>>()?;
        ensure!(!child.is_empty(), "no record of the child {{9,10}}");
        for value in &child {
            let FactorizedKind::Closure(spec) = &value.kind else {
                return Err(eyre!("not a closure record"));
            };
            let degree = spec.taylor.as_ref().map(|taylor| taylor.degree);
            ensure!(
                degree == Some(2),
                "child {{9,10}}: Taylor degree {degree:?}"
            );
        }
        Ok(())
    })
}

/// The v1.5 test cases (kept and cut holes, single-hedge `ins:` with a
/// spectator port; one directory per case with `input.txt` / `value.txt`)
/// under `GAMMALOOP_U7C_V15_DIR`: each record is read by vakint's reader,
/// evaluated by its closure entry at the record's depth and finished as the
/// hook finishes it; the value equals the case's reference (today's) value
/// through its valid_through.
#[test]
#[ignore = "reads GAMMALOOP_U7C_V15_DIR"]
fn k5_v15_cases_oracle() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let dir = std::path::PathBuf::from(
            std::env::var("GAMMALOOP_U7C_V15_DIR")
                .map_err(|_| eyre!("set GAMMALOOP_U7C_V15_DIR"))?,
        );
        let vakint = Vakint::new()?;
        let mut cases = std::fs::read_dir(&dir)?
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<std::io::Result<Vec<_>>>()?;
        cases.sort();
        let mut checked = 0;
        for case in cases.iter().filter(|path| path.join("input.txt").is_file()) {
            let text = std::fs::read_to_string(case.join("input.txt"))?;
            let value_text = std::fs::read_to_string(case.join("value.txt"))?;
            let ours = record(&text)?;
            let single = ours.insertions.iter().any(|insertion| {
                insertion.hedges.len() == 1
                    && ours
                        .vertices
                        .iter()
                        .any(|v| v.id == insertion.vertex && v.rule.starts_with("ins:"))
            });
            ensure!(single, "{}: no single-hedge ins:", case.display());
            let value = oracle_value(&vakint, &text, 1)?;
            let expected = fixture_value(&value_text)?;
            ensure!(
                equal_through(&value, &expected, valid_through(&value_text)?)?,
                "{}: the oracle value differs from the reference",
                case.display()
            );
            eprintln!("K5_V15_CASE {} PASS value={value}", case.display());
            checked += 1;
        }
        ensure!(checked > 0, "no case under {}", dir.display());
        Ok(())
    })
}
