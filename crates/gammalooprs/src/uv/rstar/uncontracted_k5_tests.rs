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
    DiagramLedger, PropagatorInsertionEvaluator, REFERENCE_MASS_SQUARED, RStarPlan, RunTally,
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
const L5_120364_DOT: &str = include_str!("testdata/u7c/l5-d120364.dot");
/// The variable that names the directory of the v1.5 cases (no default: the
/// cases live outside the repository; the deliverable's README names the
/// published directory).
const V15_DIR_VARIABLE: &str = "GAMMALOOP_U7C_V15_DIR";
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
/// The gluon surrogates of the kept + cut search whose remainders with a ct:
/// that contracts the anchor but keeps the far end failed with
/// IndexNotLocal before the locality-derived decline (2-loop and 3-loop holes).
const L5_ANCHOR_CONTRACTED: [(&str, &str); 12] = [
    ("59943", include_str!("testdata/u7c/l5-d59943.dot")),
    ("79937", include_str!("testdata/u7c/l5-d79937.dot")),
    ("44635", include_str!("testdata/u7c/l5-d44635.dot")),
    ("59942", include_str!("testdata/u7c/l5-d59942.dot")),
    ("62852", include_str!("testdata/u7c/l5-d62852.dot")),
    ("144283", include_str!("testdata/u7c/l5-d144283.dot")),
    ("59977", include_str!("testdata/u7c/l5-d59977.dot")),
    ("101342", include_str!("testdata/u7c/l5-d101342.dot")),
    ("102684", include_str!("testdata/u7c/l5-d102684.dot")),
    ("102671", include_str!("testdata/u7c/l5-d102671.dot")),
    ("59946", include_str!("testdata/u7c/l5-d59946.dot")),
    ("101213", include_str!("testdata/u7c/l5-d101213.dot")),
];
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

/// The run level is beta's RunTally (U7a fix round): a check run in which
/// the enabled kind `surrogate` compared nothing fails the run (one
/// diagram whose surrogate kind expected no record); one compared surrogate
/// record makes it pass.
#[test]
fn k5_run_with_zero_comparisons_fails() -> Result<()> {
    with_switch(check_surrogate(), || -> Result<()> {
        let mut tally = RunTally::new().ok_or_else(|| eyre!("no tally"))?;
        let ledger =
            DiagramLedger::open("k5-run-1", "k5-test").ok_or_else(|| eyre!("no ledger"))?;
        let evidence = ledger.close(true)?;
        tally.add("k5-run-1", evidence.as_ref(), true);
        let failures = tally.failures();
        ensure!(
            failures
                .iter()
                .any(|failure| failure.starts_with("kind surrogate: no record compared")),
            "a run that compared no surrogate record passed: {failures:?}"
        );
        let mut tally = RunTally::new().ok_or_else(|| eyre!("no tally"))?;
        let ledger =
            DiagramLedger::open("k5-run-2", "k5-test").ok_or_else(|| eyre!("no ledger"))?;
        uncontracted_ledger::expect(Kind::Surrogate);
        uncontracted_ledger::record_value(Kind::Surrogate, Mode::Check);
        let evidence = ledger.close(true)?;
        tally.add("k5-run-2", evidence.as_ref(), true);
        ensure!(
            tally.failures().is_empty(),
            "one compared surrogate record: {:?} ({})",
            tally.failures(),
            tally.fields()
        );
        ensure!(tally.fields().contains("surrogate:expected=1,compared=1"));
        Ok(())
    })
}

/// U7c delta pre-review N1 (BLOCKING): RunTally's diagram rule ("the
/// counterterm is nonzero, and no record of any kind was compared / taken")
/// holds for the DEFAULT kind lists (`check`, `on`) and for every list that
/// holds the default kinds, with or without `surrogate` (a kind built only
/// when listed). Base 292b5362a7 applied it to `check` / `on`; U7c d5070e583
/// had switched it off for every list without `surrogate`. A partial list
/// (not every default kind) keeps only the kind-level rule, as on the base.
#[test]
fn k5_tally_diagram_rule_holds_for_the_default_kind_lists() -> Result<()> {
    let run = |value: &str| -> Result<Vec<String>> {
        let switch = Switch::parse(Some(value)).map_err(|error| eyre!(error))?;
        with_switch(switch, || -> Result<Vec<String>> {
            let mut tally = RunTally::new().ok_or_else(|| eyre!("no tally"))?;
            // a diagram whose records all missed their hooks: nothing
            // expected in any kind, its counterterm nonzero
            let ledger =
                DiagramLedger::open("n1-empty", "k5-test").ok_or_else(|| eyre!("no ledger"))?;
            let evidence = ledger.close(true)?;
            tally.add("n1/000001", evidence.as_ref(), true);
            // another diagram compares (takes) one record of every enabled kind
            let ledger =
                DiagramLedger::open("n1-full", "k5-test").ok_or_else(|| eyre!("no ledger"))?;
            for kind in switch.enabled_kinds() {
                uncontracted_ledger::expect(kind);
                uncontracted_ledger::record_value(kind, switch.mode);
            }
            let evidence = ledger.close(true)?;
            tally.add("n1/000002", evidence.as_ref(), true);
            Ok(tally.failures())
        })
    };
    for value in [
        "check",
        "on",
        "1",
        "check:root,child,ir-hard",
        "on:root,child,ir-hard,quotient",
        "check:root,child,ir-hard,surrogate",
        "on:a,b,c,e",
    ] {
        let verb = if value.starts_with("check") {
            "compared"
        } else {
            "taken"
        };
        let failures = run(value)?;
        ensure!(
            failures
                == vec![format!(
                    "n1/000001: the counterterm is nonzero, and no record of any kind was {verb}"
                )],
            "{value}: {failures:#?}"
        );
    }
    for value in ["check:ir-hard", "on:root,surrogate", "check:surrogate"] {
        let failures = run(value)?;
        ensure!(
            failures.is_empty(),
            "{value} (a partial list): {failures:#?}"
        );
    }
    Ok(())
}

/// The base's run-tally test with `surrogate` LISTED in addition to the
/// default kinds (the variant U7c had edited the base's test into; the base's
/// own test stands unchanged): a nonzero diagram with nothing compared in any
/// kind fails the run, a surrogate kind that compared nothing fails it.
#[test]
fn k5_tally_with_surrogate_listed() -> Result<()> {
    test_initialise()?;
    let check = Switch::new(
        Mode::Check,
        &[Kind::Root, Kind::Child, Kind::IrHard, Kind::Surrogate],
    );
    let kind = |compared: usize| {
        serde_json::json!({"expected": compared, "attempted": compared, "uncontracted": compared,
                           "zero": 0, "equal": compared,
                           "classes": if compared > 0 { serde_json::json!({"compared": compared}) } else { serde_json::json!({}) },
                           "outcome": "ok"})
    };
    with_switch(check, || -> Result<()> {
        let identity = uncontracted_ledger::current_identity();
        let evidence = |graph: &str, root: usize, surrogate: usize| {
            serde_json::json!({
                "schema": uncontracted_ledger::EVIDENCE_SCHEMA, "identity": identity, "graph": graph,
                "diagram": "000001", "outcome": "ok", "inert": {},
                "kinds": {"root": kind(root), "child": kind(root), "ir-hard": kind(root),
                          "surrogate": kind(surrogate)},
            })
        };
        let mut tally = RunTally::new().ok_or_else(|| eyre!("no tally"))?;
        tally.add("g/000001", Some(&evidence("g", 1, 1)), true);
        ensure!(tally.failures().is_empty(), "{:#?}", tally.failures());
        tally.add("z/000009", Some(&evidence("z", 0, 0)), true);
        let failures = tally.failures();
        ensure!(
            failures
                == vec![
                    "z/000009: the counterterm is nonzero, and no record of any kind was compared"
                        .to_owned()
                ],
            "{failures:#?}"
        );
        let mut tally = RunTally::new().ok_or_else(|| eyre!("no tally"))?;
        tally.add("g/000001", Some(&evidence("g", 1, 0)), true);
        let failures = tally.failures();
        ensure!(
            failures.len() == 1 && failures[0].starts_with("kind surrogate: no record compared"),
            "{failures:#?}"
        );
        Ok(())
    })
}

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

/// The directory of the v1.5 cases from the value of [`V15_DIR_VARIABLE`]:
/// REQUIRED (there is no default directory; unset or empty is an error that
/// says what to set).
fn v15_dir(value: Option<String>) -> Result<std::path::PathBuf> {
    match value.filter(|value| !value.trim().is_empty()) {
        Some(value) => Ok(std::path::PathBuf::from(value)),
        None => Err(eyre!(
            "set {V15_DIR_VARIABLE} to the directory of the v1.5 cases: one sub-directory per case with input.txt and value.txt, or one per publication stamp with its cases below (the published directory: the deliverable's README)"
        )),
    }
}

/// The v1.5 oracle's directory is required (r4: r3's default, a host path,
/// is dropped): unset or empty is an error naming the variable and what it
/// must hold; a value is taken as given.
#[test]
fn k5_v15_oracle_requires_its_directory_variable() -> Result<()> {
    for unset in [None, Some(String::new()), Some("  ".to_owned())] {
        let error = v15_dir(unset.clone()).expect_err("no directory without the variable");
        let text = error.to_string();
        ensure!(
            text.starts_with(&format!(
                "set {V15_DIR_VARIABLE} to the directory of the v1.5 cases"
            )) && text.contains("input.txt and value.txt"),
            "{unset:?}: {text}"
        );
    }
    ensure!(v15_dir(Some("cases".into()))? == std::path::PathBuf::from("cases"));
    Ok(())
}

/// No host path in the K5 sources (r4; the r3 oracle's default directory was
/// one): no absolute path under the host's shared trees in uncontracted_k5.rs,
/// its tests, insertions.rs or the vendored K5 test data.
#[test]
fn k5_sources_hold_no_host_path() -> Result<()> {
    let rstar = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/uv/rstar");
    let mut files = vec![
        rstar.join("uncontracted_k5.rs"),
        rstar.join("uncontracted_k5_tests.rs"),
        rstar.join("insertions.rs"),
    ];
    for entry in std::fs::read_dir(rstar.join("testdata/u7c"))? {
        files.push(entry?.path());
    }
    // the needles are built, so that this file does not match itself
    let needles = ["common", "shared", "home", "tmp"].map(|root| format!("{}{root}/", '/'));
    let mut hits = Vec::new();
    for file in &files {
        let text = std::fs::read_to_string(file)?;
        for (n, line) in text.lines().enumerate() {
            if needles.iter().any(|needle| line.contains(needle.as_str())) {
                hits.push(format!("{}:{}: {line}", file.display(), n + 1));
            }
        }
    }
    ensure!(files.len() > 3, "the K5 test data was not found");
    ensure!(hits.is_empty(), "host paths: {hits:#?}");
    Ok(())
}

/// The published v1.5 test cases (kept and cut holes, single-hedge `ins:`
/// with a spectator port; one directory per case with `input.txt` /
/// `value.txt`) under the directory named by `GAMMALOOP_U7C_V15_DIR`
/// ([`V15_DIR_VARIABLE`], REQUIRED: unset is an error that says what to set;
/// the directory may also hold one sub-directory per publication stamp with
/// its cases below, as the published directory does): each record is read
/// by vakint's reader,
/// evaluated by its closure entry at the record's depth and finished as the
/// hook finishes it; the value equals the case's reference (today's) value
/// through its valid_through.
#[test]
#[ignore = "reads the v1.5 cases under GAMMALOOP_U7C_V15_DIR (required)"]
fn k5_v15_cases_oracle() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let dir = v15_dir(std::env::var(V15_DIR_VARIABLE).ok())?;
        let vakint = Vakint::new()?;
        // the cases: `dir` itself, its sub-directories and theirs (a stamp
        // directory, then its cases)
        let mut cases = Vec::new();
        let mut stack = vec![(dir.clone(), 0)];
        while let Some((path, depth)) = stack.pop() {
            if path.join("input.txt").is_file() {
                cases.push(path);
            } else if depth < 2 {
                for entry in std::fs::read_dir(&path)? {
                    let entry = entry?.path();
                    if entry.is_dir() {
                        stack.push((entry, depth + 1));
                    }
                }
            }
        }
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

/// Handed to beta as a FAILING test (patch 0010, ignored on 292b5362a7);
/// beta's follow-up tip 0484ddd9f6 (0c61cb1f1, item B: a setting-caused
/// decline of every record of one record class fails the diagram) makes it
/// PASS, so it runs by default since the rebase onto that tip. With the
/// surrogate root NOT closed (the evaluator's open-root entry `pole_part`,
/// as under GAMMALOOP_SURROGATE_CLOSED_ROOT=0) on L5 3793:
/// (1) with the K5 site's settings check (production): the site's error
///     `surrogate_settings` fires FIRST, at the surrogate's site, before any
///     record exists (none is captured), and the evaluation fails;
/// (2) with that check bypassed: every surrogate ROOT (K1-shaped) record
///     declines `open_root` while the surrogate's children (K2-shaped)
///     compare, and beta's ledger rule fails the diagram at its close
///     ("... setting-caused ..."). The two do not conflict: the K5 check is
///     earlier, beta's rule is the general one behind it.
#[test]
fn k5_beta_ledger_misses_settings_wide_root_declines() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let original = colour_factored(&load(L5_3793_DOT)?)?;
        let evaluator = PropagatorInsertionEvaluator::new(settings());
        with_switch(check_surrogate(), || -> Result<()> {
            // (1) the K5 site check fires first
            let ledger = DiagramLedger::open(original.name(), "k5-open-root-site")
                .ok_or_else(|| eyre!("no ledger"))?;
            let (value, captured) = with_capture(|| evaluator.pole_part(&original, &vakint));
            let Err(error) = value else {
                return Err(eyre!("the open surrogate root passed the K5 site check"));
            };
            ensure!(
                format!("{error:#}").contains("surrogate_settings"),
                "not the K5 site error: {error:#}"
            );
            ensure!(
                captured
                    .iter()
                    .all(|attempt| attempt.kind != Kind::Surrogate),
                "{} surrogate records exist before the site check",
                captured.len()
            );
            let _ = ledger.close(false);
            // (2) bypassed: beta's ledger rule
            let ledger = DiagramLedger::open(original.name(), "k5-open-root")
                .ok_or_else(|| eyre!("no ledger"))?;
            let (value, captured) = k5::without_settings_check(|| {
                with_capture(|| evaluator.pole_part(&original, &vakint))
            });
            let roots = captured
                .iter()
                .filter(|attempt| attempt.outcome == "declined:open_root")
                .count();
            let compared = captured
                .iter()
                .filter(|attempt| attempt.equal == Some(true))
                .count();
            eprintln!("K5_OPEN_ROOT roots_declined={roots} compared={compared}");
            ensure!(roots > 0, "no root declined open_root");
            let closed = ledger.close(value.is_ok());
            match closed {
                Err(error) => ensure!(
                    format!("{error:#}").contains("setting-caused"),
                    "the diagram failed for another reason: {error:#}"
                ),
                Ok(_) => {
                    return Err(eyre!(
                        "the diagram passed with every surrogate root declined by a setting ({roots} open_root, {compared} children compared)"
                    ));
                }
            }
            Ok(())
        })
    })
}

/// The 12 gluon surrogates that failed with IndexNotLocal (a ct: that
/// contracts the anchor while the record keeps the far end: C_H's far slot
/// non-local): through the whole evaluator in check / on / off and BOTH
/// hole-cut settings, no diagram fails, no record errors, the case declines
/// by design (`surrogate_split_chain`) before its record, every other record
/// compares equal, on == off through eps^-1.
#[test]
fn k5_anchor_contracted_far_kept_declines_by_design() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let mut failures = Vec::new();
        for mode in [k5::HoleCut::Decline, k5::HoleCut::Form] {
            for (name, dot) in L5_ANCHOR_CONTRACTED {
                let original = colour_factored(&load(dot)?)?;
                match k5::with_hole_cut(mode, || whole_graph_summary(&original, &vakint, 1)) {
                    Ok(summary) if summary.keys().any(|k| k.starts_with("error:")) => {
                        failures.push(format!("{mode:?} {name}: {summary:?}"))
                    }
                    Ok(summary) => {
                        eprintln!("K5_ANCHOR_CONTRACTED mode={mode:?} {name} PASS {summary:?}")
                    }
                    Err(error) => failures.push(format!("{mode:?} {name}: {error:#}")),
                }
            }
        }
        ensure!(failures.is_empty(), "{failures:#?}");
        Ok(())
    })
}

/// The scalar chains of a surrogate graph: per chain, its vertices
/// (`scalar_insertion_0..L` from the coefficient vertex) and the bubble edges
/// of each link, read from the graph's own vertex names and edges.
fn surrogate_chains(
    graph: &Graph,
    provenance: &Provenance,
) -> Vec<(Vec<usize>, Vec<BTreeSet<usize>>)> {
    use linnet::half_edge::involution::HedgePair;
    let mut name = BTreeMap::new();
    let mut ends = BTreeMap::new();
    for (pair, edge, _) in graph.iter_edges() {
        if let HedgePair::Paired { source, sink } = pair {
            let (s, t) = (graph.node_id(source), graph.node_id(sink));
            name.insert(usize::from(s), graph[s].name.value.clone());
            name.insert(usize::from(t), graph[t].name.value.clone());
            ends.insert(usize::from(edge), (usize::from(s), usize::from(t)));
        }
    }
    provenance
        .holes
        .iter()
        .map(|hole| {
            let mut vertices = vec![hole.vertex];
            let mut links = Vec::new();
            for i in 0..hole.loops {
                let current = vertices[i];
                let next = format!("scalar_insertion_{}", i + 1);
                let mut by = BTreeMap::<usize, BTreeSet<usize>>::new();
                for (edge, (a, b)) in &ends {
                    let other = if *a == current {
                        *b
                    } else if *b == current {
                        *a
                    } else {
                        continue;
                    };
                    if !vertices.contains(&other) && name.get(&other) == Some(&next) {
                        by.entry(other).or_default().insert(*edge);
                    }
                }
                let (v, edges) = by.into_iter().find(|(_, e)| e.len() == 2).expect("a link");
                vertices.push(v);
                links.push(edges);
            }
            (vertices, links)
        })
        .collect()
}

/// The class of a surrogate record against one chain (from its text): the
/// anchor kept / contracted / absent, the far end kept / contracted / absent.
fn hole_class(
    record: &FactorizedNumerator,
    chain: &(Vec<usize>, Vec<BTreeSet<usize>>),
) -> &'static str {
    let (vertices, links) = chain;
    let rule = |id: usize| {
        record
            .vertices
            .iter()
            .find(|v| v.id == id)
            .map(|v| v.rule.clone())
    };
    let ct_edges = record
        .vertices
        .iter()
        .filter_map(|v| v.rule.strip_prefix("ct:edges-"))
        .map(|edges| {
            edges
                .split('-')
                .filter_map(|e| e.parse::<usize>().ok())
                .collect::<BTreeSet<_>>()
        })
        .collect::<Vec<_>>();
    let contracted =
        |link: &BTreeSet<usize>| ct_edges.iter().any(|support| !support.is_disjoint(link));
    let anchor = match rule(vertices[0]) {
        Some(r) if r.starts_with("ins:") => "kept",
        _ if contracted(&links[0]) => "contracted",
        _ => "absent",
    };
    let last = *vertices.last().unwrap();
    let far = match rule(last) {
        Some(r) if r.starts_with("scalar_insertion_") => "kept",
        _ if contracted(links.last().unwrap()) => "contracted",
        _ => "absent",
    };
    match (anchor, far) {
        ("kept", "kept") => "whole",
        ("kept", "absent") => "kept+cut",
        ("kept", "contracted") => "kept,far-contracted",
        ("contracted", "contracted") => "contracted",
        ("contracted", _) => "contracted-cut",
        ("absent", _) => "no-hole",
        _ => "other",
    }
}

/// Every (child, ct) combination the R* forest produces for chains of
/// L = 2, 3, 4 bubbles (L5 143019, 143190, 163502 and the 12 gluon graphs of
/// 2- and 3-loop holes), in both hole-cut settings: each record named by
/// class (whole / kept+cut / kept,far-contracted / contracted / contracted-cut
/// / no-hole, or its designed decline), every EMITTED record read by vakint's
/// reader (its text) and compared equal, every other one declined by design
/// with no record text (or `zero`: a vanishing numerator, no record), none an
/// error; kept+cut emitted only in form mode.
#[test]
fn k5_every_child_ct_combination_of_cut_chains_is_stated_or_declined() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let graphs = L5_MULTI_LOOP_HOLES[..3]
            .iter()
            .chain(L5_ANCHOR_CONTRACTED.iter())
            .copied()
            .collect::<Vec<_>>();
        let evaluator = PropagatorInsertionEvaluator::new(settings());
        let mut failures = Vec::new();
        let mut classes = BTreeMap::<String, usize>::new();
        for mode in [k5::HoleCut::Decline, k5::HoleCut::Form] {
            for (name, dot) in &graphs {
                let original = colour_factored(&load(dot)?)?;
                let prepared = prepare(dot, &vakint)?;
                let graph = prepared.surrogate.to_gamma_loop_graph(None, true)?;
                let chains = surrogate_chains(&graph, &prepared.provenance);
                let captured = with_switch(check_surrogate(), || -> Result<Vec<Captured>> {
                    let ledger = DiagramLedger::open(original.name(), "k5-enum");
                    let (value, captured) = k5::with_hole_cut(mode, || {
                        with_capture(|| evaluator.pole_part_closing_root(&original, &vakint, true))
                    });
                    if let Some(ledger) = ledger {
                        ledger.close(value.is_ok())?;
                    }
                    value?;
                    Ok(captured)
                })?;
                for attempt in captured.iter().filter(|a| a.kind == Kind::Surrogate) {
                    let class = match &attempt.text {
                        Some(text) => {
                            let ours = record(text)?;
                            if attempt.equal != Some(true) {
                                failures.push(format!(
                                    "{mode:?} {name} id {}: emitted, not equal",
                                    attempt.id
                                ));
                            }
                            let class = chains
                                .iter()
                                .map(|chain| hole_class(&ours, chain))
                                .find(|class| *class != "no-hole")
                                .unwrap_or("no-hole");
                            if class == "kept+cut" && mode == k5::HoleCut::Decline {
                                failures.push(format!(
                                    "{name} id {}: kept+cut emitted in decline mode",
                                    attempt.id
                                ));
                            }
                            format!("emitted:{class}")
                        }
                        None if attempt.outcome.starts_with("declined:") => attempt.outcome.clone(),
                        // a vanishing numerator: no record, counted `zero` by the ledger
                        None if attempt.outcome == "zero" => "zero".to_string(),
                        None => {
                            failures.push(format!(
                                "{mode:?} {name} id {}: {}",
                                attempt.id, attempt.outcome
                            ));
                            attempt.outcome.clone()
                        }
                    };
                    *classes.entry(format!("{mode:?} {class}")).or_default() += 1;
                }
            }
        }
        for (class, n) in &classes {
            eprintln!("K5_ENUM {class} {n}");
        }
        ensure!(failures.is_empty(), "{failures:#?}");
        for class in [
            "Form emitted:kept+cut",
            "Form emitted:kept,far-contracted",
            "Decline declined:surrogate_hole_cut",
            "Decline declined:surrogate_split_chain",
        ] {
            ensure!(
                classes.contains_key(class),
                "no record of class {class}: {classes:?}"
            );
        }
        // the two settings differ exactly by kept+cut emitted <-> declined
        let in_mode = |mode: &str| {
            classes
                .iter()
                .filter_map(|(key, n)| {
                    key.strip_prefix(mode).map(|class| {
                        (
                            class.replace("emitted:kept+cut", "declined:surrogate_hole_cut"),
                            *n,
                        )
                    })
                })
                .collect::<BTreeMap<_, _>>()
        };
        ensure!(
            in_mode("Decline ") == in_mode("Form "),
            "the hole-cut settings differ beyond kept+cut: {classes:?}"
        );
        Ok(())
    })
}

/// The half-edges `h` of the bispinor slots `bis(.., hedge(h, k))` of an atom.
fn bispinor_hedges(atom: &Atom) -> BTreeSet<usize> {
    free_hedge_slots(atom)
        .iter()
        .filter(|slot| slot.contains("::bis("))
        .filter_map(|slot| {
            let tail = &slot[slot.find("::hedge(")? + "::hedge(".len()..];
            tail.split(',').next()?.trim().parse().ok()
        })
        .collect()
}

/// Beta's question (21 Nov, 17:0xZ): a CONTRACTED and CUT quark hole (a ct:
/// child holding the anchor, the far end's half-edge h_out absent from the
/// record) would need a bispinor port no projector of the record closes
/// (open format: OQ-6 / PORTS F5). It is NEVER emitted, in either hole-cut
/// setting, through the whole evaluator in check mode:
/// (a) L5 120364 (quark carrier, 3-loop hole; the 198-graph search's case
///     with a NON-ZERO anchor child Delta): the record declines
///     `surrogate_split_chain` (closed list) with the case named in its
///     detail, before any record text;
/// (b) L5 175270 (quark carrier, 4-loop hole): its contracted-cut records are
///     `zero` (the anchor child's Delta vanishes);
/// in both graphs and settings no emitted surrogate record has a ct:
/// insertion with a bispinor slot of a half-edge the record lacks, and every
/// attempt is compared equal, zero or declined by design.
#[test]
fn k5_contracted_cut_quark_holes_are_never_emitted() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let evaluator = PropagatorInsertionEvaluator::new(settings());
        let mut counts = BTreeMap::<String, usize>::new();
        for (name, dot) in [("120364", L5_120364_DOT), L5_MULTI_LOOP_HOLES[3]] {
            let original = colour_factored(&load(dot)?)?;
            for mode in [k5::HoleCut::Decline, k5::HoleCut::Form] {
                let captured = with_switch(check_surrogate(), || -> Result<Vec<Captured>> {
                    let ledger = DiagramLedger::open(original.name(), "k5-cut-quark");
                    let (value, captured) = k5::with_hole_cut(mode, || {
                        with_capture(|| evaluator.pole_part_closing_root(&original, &vakint, true))
                    });
                    if let Some(ledger) = ledger {
                        ledger.close(value.is_ok())?;
                    }
                    value?;
                    Ok(captured)
                })?;
                for attempt in captured.iter().filter(|a| a.kind == Kind::Surrogate) {
                    ensure!(
                        attempt.equal == Some(true)
                            || attempt.outcome == "zero"
                            || (attempt.outcome.starts_with("declined:")
                                && uncontracted::is_designed_code(
                                    &attempt.outcome["declined:".len()..]
                                )),
                        "{name} {mode:?} id {}: {} equal {:?} {:?}",
                        attempt.id,
                        attempt.outcome,
                        attempt.equal,
                        attempt.detail
                    );
                    let cut_quark = attempt.outcome == "declined:surrogate_split_chain"
                        && attempt
                            .detail
                            .as_deref()
                            .is_some_and(|detail| detail.contains("a contracted cut quark hole"));
                    if cut_quark {
                        ensure!(attempt.text.is_none(), "{name}: a record text was built");
                    }
                    let key = if cut_quark {
                        "declined:cut-quark".to_owned()
                    } else {
                        attempt.outcome.clone()
                    };
                    *counts.entry(format!("{name} {mode:?} {key}")).or_default() += 1;
                    let Some(text) = &attempt.text else { continue };
                    let ours = record(text)?;
                    let present = ours
                        .edges
                        .iter()
                        .flat_map(|edge| [edge.source_hedge, edge.target_hedge])
                        .flatten()
                        .collect::<BTreeSet<_>>();
                    for vertex in ours.vertices.iter().filter(|v| v.rule.starts_with("ct:")) {
                        for insertion in ours.insertions.iter().filter(|i| i.vertex == vertex.id) {
                            let slots =
                                bispinor_hedges(&(&insertion.coefficient * &insertion.factor));
                            let absent = slots.difference(&present).collect::<Vec<_>>();
                            ensure!(
                                absent.is_empty(),
                                "{name} {mode:?} id {}: the ct: insertion at {} carries the bispinor slots {absent:?} of half-edges the record lacks (an emitted contracted cut quark hole)",
                                attempt.id,
                                vertex.id
                            );
                        }
                    }
                }
            }
        }
        for (outcome, n) in &counts {
            eprintln!("K5_CUT_QUARK {outcome} {n}");
        }
        for mode in ["Decline", "Form"] {
            ensure!(
                counts.contains_key(&format!("120364 {mode} declined:cut-quark")),
                "120364 {mode}: no contracted cut quark decline: {counts:?}"
            );
            ensure!(
                counts.contains_key(&format!("175270 {mode} zero")),
                "175270 {mode}: no zero: {counts:?}"
            );
        }
        Ok(())
    })
}

/// U7c delta pre-review N3, after the gate's finding on 5e97a6eb5 (the mode
/// in the Delta-cache key split fills that r1 shared across kinds, so child
/// lookups missed): the hole-cut mode is resolved ONCE, at the switch's
/// initialisation, and stored immutably; it is in the resume identity (test
/// k5_hole_cut_mode_is_in_the_resume_identity) and NOT in the in-process
/// Delta-cache key, which is r1's again (k5_delta_cache_is_in_process_only
/// states why that is safe). Shown here:
/// - parsing at start-up (a typo is an error there, not a panic mid-run);
/// - a second resolution is an error, and changing the variable after
///   start-up changes nothing;
/// - structurally, in the production sources of uv/rstar (test files
///   excluded): the variable is read at exactly one place, inside the
///   `ENVIRONMENT` initialiser, through the one resolver; nothing else reads
///   the environment for it; the only override is `#[cfg(test)]`;
/// - the Delta-cache keys do not depend on the mode (inside and outside a
///   surrogate's R*, on and check, with and without `child`) and equal r1's.
#[test]
fn k5_hole_cut_mode_is_resolved_once_at_start_up_and_not_in_the_cache_key() -> Result<()> {
    use k5::{HoleCut, parse_hole_cut, with_hole_cut};
    // start-up parsing
    for (value, mode) in [
        (None, HoleCut::Decline),
        (Some(""), HoleCut::Decline),
        (Some("decline"), HoleCut::Decline),
        (Some(" form "), HoleCut::Form),
    ] {
        ensure!(parse_hole_cut(value, true) == Ok(mode), "{value:?}");
    }
    ensure!(parse_hole_cut(Some("auto"), true) == Ok(HoleCut::Form));
    ensure!(parse_hole_cut(Some("auto"), false) == Ok(HoleCut::Decline));
    let error = parse_hole_cut(Some("froms"), true).expect_err("a typo is refused");
    ensure!(error.contains(k5::HOLE_CUT_VARIABLE), "{error}");
    // the start-up (the test environment's switch is off: nothing resolved
    // there), one resolution, a second one refused, a later change of the
    // variable without effect
    uncontracted::validate()?;
    let before = k5::hole_cut();
    if k5::resolved_hole_cut().is_none() {
        ensure!(k5::resolve_hole_cut(None, false) == Ok(HoleCut::Decline));
    }
    let second = k5::resolve_hole_cut(Some("form"), true).expect_err("a second resolution");
    ensure!(second.contains("second resolution"), "{second}");
    ensure!(k5::resolved_hole_cut() == Some(HoleCut::Decline));
    // SAFETY: nextest runs every test in its own process; the variable is
    // removed again below
    unsafe { std::env::set_var(k5::HOLE_CUT_VARIABLE, "form") };
    let after = k5::hole_cut();
    unsafe { std::env::remove_var(k5::HOLE_CUT_VARIABLE) };
    ensure!(
        before == HoleCut::Decline && after == before,
        "{before:?} -> {after:?}"
    );
    // structure of the production sources
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/uv/rstar");
    let mut sources = Vec::new();
    for entry in std::fs::read_dir(&dir)? {
        let path = entry?.path();
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or("")
            .to_owned();
        if name.ends_with(".rs") && !name.ends_with("_tests.rs") {
            sources.push((name, std::fs::read_to_string(&path)?));
        }
    }
    let count = |needle: &str| {
        sources
            .iter()
            .map(|(_, text)| text.matches(needle).count())
            .sum::<usize>()
    };
    ensure!(
        count("\"GAMMALOOP_UNCONTRACTED_K5_HOLE_CUT\"") == 1,
        "the variable's name outside its constant"
    );
    ensure!(
        count("read(k5::HOLE_CUT_VARIABLE)") == 1,
        "not exactly one read of the variable"
    );
    ensure!(
        count("resolve_hole_cut(") == 2,
        "not exactly one call of the resolver"
    );
    ensure!(
        count("start_up_hole_cut(") == 2,
        "not exactly one start-up call"
    );
    ensure!(
        count("parse_hole_cut(") == 2,
        "parse_hole_cut called outside the resolver"
    );
    ensure!(
        count("RESOLVED_HOLE_CUT.set(") == 1,
        "a second writer of the mode"
    );
    for (name, text) in &sources {
        for line in text.lines() {
            ensure!(
                !(line.contains("HOLE_CUT")
                    && (line.contains("env::var") || line.contains("var_os"))),
                "{name}: a direct environment read of the hole-cut variable: {line}"
            );
        }
    }
    let (_, main) = sources
        .iter()
        .find(|(name, _)| name == "uncontracted.rs")
        .ok_or_else(|| eyre!("no uncontracted.rs"))?;
    let start = main
        .find("static ENVIRONMENT")
        .ok_or_else(|| eyre!("no ENVIRONMENT"))?;
    let end = start + main[start..].find("\n});").ok_or_else(|| eyre!("no end"))?;
    let read = main
        .find("read(k5::HOLE_CUT_VARIABLE)")
        .ok_or_else(|| eyre!("no read"))?;
    ensure!(
        start < read && read < end,
        "the read is not in the start-up initialiser"
    );
    let (_, k5_text) = sources
        .iter()
        .find(|(name, _)| name == "uncontracted_k5.rs")
        .ok_or_else(|| eyre!("no uncontracted_k5.rs"))?;
    let test_override = k5_text
        .find("static TEST_HOLE_CUT")
        .ok_or_else(|| eyre!("no override"))?;
    ensure!(
        k5_text[..test_override]
            .trim_end()
            .ends_with("thread_local! {")
            && k5_text[..test_override]
                .trim_end()
                .trim_end_matches("thread_local! {")
                .trim_end()
                .ends_with("#[cfg(test)]"),
        "the thread override is not #[cfg(test)]"
    );
    // the Delta-cache keys: independent of the mode, equal to r1's
    let key = |switch: Switch, mode: HoleCut, inside: bool| {
        with_switch(switch, || {
            with_hole_cut(mode, || {
                if inside {
                    in_dressed_scope(uncontracted::cache_key_suffix)
                } else {
                    uncontracted::cache_key_suffix()
                }
            })
        })
    };
    let gate =
        Switch::parse(Some("check:child,ir-hard,surrogate")).map_err(|error| eyre!(error))?;
    let gate_on =
        Switch::parse(Some("on:child,ir-hard,surrogate")).map_err(|error| eyre!(error))?;
    for (switch, inside, outside) in [
        (
            on_surrogate(),
            Some("#uncontracted=surrogate"),
            Some("#uncontracted=surrogate#uncontracted-inert=child"),
        ),
        (check_surrogate(), None, Some("#uncontracted-inert=child")),
        // with `child` enabled an ordinary child and a surrogate-scope fill
        // share the key (the cross-kind hits r1 had)
        (gate, None, None),
        (
            gate_on,
            Some("#uncontracted=child,ir-hard,surrogate"),
            Some("#uncontracted=child,ir-hard,surrogate"),
        ),
    ] {
        for mode in [HoleCut::Decline, HoleCut::Form] {
            ensure!(
                key(switch, mode, true).as_deref() == inside,
                "{switch:?} {mode:?} inside"
            );
            ensure!(
                key(switch, mode, false).as_deref() == outside,
                "{switch:?} {mode:?} outside"
            );
        }
    }
    Ok(())
}

/// The Delta cache never crosses a process, which is why a per-process
/// constant (the hole-cut mode) belongs in the resume identity and not in
/// its key. The type is `CountertermCache` (uv/counterterm_cache.rs:323-334):
/// `#[derive(Default)]`, a `Mutex<HashMap<Key, Entry>>` in memory, no serde,
/// no file I/O, and no static holds one. Beta (771c059d) builds ONE per
/// engine run (src/main.rs:3874-3880, `--delta-cache on`) or per integrator
/// (src/integrate.rs:782), shared by that process's worker threads through
/// an `Arc`; without one, a plan builds its own per call (plan.rs). What
/// crosses a process are the per-diagram checkpoints, bound by the resume
/// identity. Checked structurally here.
#[test]
fn k5_delta_cache_is_in_process_only() -> Result<()> {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/uv");
    let cache = std::fs::read_to_string(root.join("counterterm_cache.rs"))?;
    ensure!(
        cache.contains("#[derive(Default)]\npub struct CountertermCache {\n    entries: Mutex<HashMap<Key, Entry>>,"),
        "CountertermCache is not the plain in-memory map"
    );
    for needle in [
        "std::fs",
        "File::",
        "Deserialize",
        "OpenOptions",
        "BufWriter",
        "write_all",
    ] {
        ensure!(
            !cache.contains(needle),
            "counterterm_cache.rs uses {needle}"
        );
    }
    let mut stack = vec![root];
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir)? {
            let path = entry?.path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().is_some_and(|e| e == "rs")
                && !path
                    .file_name()
                    .and_then(|name| name.to_str())
                    .is_some_and(|name| name.ends_with("_tests.rs") || name == "tests.rs")
            {
                let text = std::fs::read_to_string(&path)?;
                for line in text.lines() {
                    ensure!(
                        !(line.contains("static") && line.contains("CountertermCache")),
                        "{}: a static holds a CountertermCache: {line}",
                        path.display()
                    );
                }
            }
        }
    }
    Ok(())
}

/// U7c delta pre-review N3 (resume identity): the hole-cut mode is part of
/// the resume identity when `surrogate` is enabled: a run resumed under the
/// other mode recomputes (its evidence's identity differs), the same mode
/// reuses it. Without `surrogate` the identity is the base's.
#[test]
fn k5_hole_cut_mode_is_in_the_resume_identity() -> Result<()> {
    use k5::{HoleCut, with_hole_cut};
    // identity
    let identity = |switch: Switch, mode: HoleCut| {
        with_switch(switch, || {
            with_hole_cut(mode, uncontracted_ledger::current_identity)
        })
    };
    let (decline, form) = (
        identity(check_surrogate(), HoleCut::Decline),
        identity(check_surrogate(), HoleCut::Form),
    );
    ensure!(
        decline.as_deref() == Some("check:surrogate:strict:hole_cut=decline")
            && form.as_deref() == Some("check:surrogate:strict:hole_cut=form"),
        "{decline:?} {form:?}"
    );
    let default = Switch::parse(Some("check")).map_err(|error| eyre!(error))?;
    ensure!(
        identity(default, HoleCut::Decline).as_deref() == Some("check:root,child,ir-hard:strict")
            && identity(default, HoleCut::Form) == identity(default, HoleCut::Decline),
        "the default identity changed"
    );
    // a resumed run under the other mode recomputes, under the same mode it
    // reuses the evidence
    let evidence = with_switch(check_surrogate(), || {
        with_hole_cut(HoleCut::Form, || -> Result<serde_json::Value> {
            let ledger =
                DiagramLedger::open("n3-graph", "000003").ok_or_else(|| eyre!("no ledger"))?;
            uncontracted_ledger::expect(Kind::Surrogate);
            uncontracted_ledger::record_value(Kind::Surrogate, Mode::Check);
            ledger.close(true)?.ok_or_else(|| eyre!("no evidence"))
        })
    })?;
    let resumed = |mode: HoleCut| {
        with_switch(check_surrogate(), || {
            with_hole_cut(mode, || {
                uncontracted_ledger::resumed("n3-graph", "000003", Some(&evidence))
            })
        })
    };
    ensure!(resumed(HoleCut::Form) == uncontracted_ledger::Resumed::Checked);
    match resumed(HoleCut::Decline) {
        uncontracted_ledger::Resumed::Unchecked(reason) => ensure!(
            reason.starts_with("identity:") && reason.contains("hole_cut=form"),
            "{reason}"
        ),
        other => {
            return Err(eyre!(
                "the decline run reused the form run's evidence: {other:?}"
            ));
        }
    }
    Ok(())
}

/// U7c delta pre-review N4: an open slot of a half-edge that is not a hole's
/// inner half-edge is an ERROR (`surrogate_stray_slot`), never a spectator
/// port, in mode on (no product check there): (a) a Minkowski slot pair of
/// absent half-edges in a ct: child's Delta, (b) the same pair on a
/// propagator outside every chain, (c) a bispinor slot of an absent
/// non-chain half-edge in a ct: Delta (not the designed contracted cut quark
/// decline, which is only a hole's far end). On every graph of the suite the
/// port set equals the one the record's structure implies
/// (`surrogate_ports` would be an error: covered by every whole-graph test).
#[test]
fn k5_stray_open_slots_are_errors_not_ports() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let prepared = prepare(L5_22807_DOT, &vakint)?;
        let parse = |text: &str| symbolica::try_parse!(text).map_err(|error| eyre!("{error}"));
        let mink = parse(
            "spenso::{symmetric,linear,real,spenso::tensor}::g(spenso::{spenso::representation,spenso::self_dual}::mink(gammalooprs::{}::dim,gammalooprs::{spenso::index}::hedge(98,1)),spenso::{spenso::representation,spenso::self_dual}::mink(gammalooprs::{}::dim,gammalooprs::{spenso::index}::hedge(99,1)))",
        )?;
        let bis = parse(
            "gammalooprs::{}::k5_spinor_probe(spenso::{spenso::representation,spenso::self_dual}::bis(4,gammalooprs::{spenso::index}::hedge(99,1)))",
        )?;
        let expect_error = |attempt: uncontracted::Attempt, what: &str| -> Result<()> {
            ensure!(attempt.text.is_none(), "{what}: a record was built");
            match &attempt.outcome {
                Outcome::Error(error) => ensure!(
                    error.code == "surrogate_stray_slot" && !error.designed,
                    "{what}: {error:?}"
                ),
                other => return Err(eyre!("{what}: not an error: {other:?}")),
            }
            Ok(())
        };
        for (what, probe) in [
            ("(a) mink in a ct: Delta", &mink),
            ("(c) bis in a ct: Delta", &bis),
        ] {
            let mut case = root_case(&prepared, K5_SURROGATE)?;
            ensure!(!case.deltas.is_empty());
            for delta in case.deltas.values_mut() {
                *delta = &*delta * probe;
            }
            for mode in [k5::HoleCut::Decline, k5::HoleCut::Form] {
                let attempt = k5::with_hole_cut(mode, || {
                    root_attempt(&case, &vakint, Some(case.provenance.clone()), false)
                })?;
                expect_error(attempt, what)?;
            }
        }
        // (b) a propagator outside every chain (not a bubble line, factor not 1)
        let case = root_case(&prepared, K5_SURROGATE)?;
        let done = std::rc::Rc::new(std::cell::Cell::new(false));
        let hook: uncontracted::InsertionHook = {
            let (done, mink) = (done.clone(), mink.clone());
            std::rc::Rc::new(move |kind: Kind, _node: usize, _delta: &mut Atom, edges: &mut [vakint::p_integral::factorized::FactorizedEdge]| {
                if kind != Kind::Surrogate || done.get() {
                    return;
                }
                if let Some(edge) = edges
                    .iter_mut()
                    .find(|edge| !edge.is_external() && edge.factor != Atom::num(1))
                {
                    edge.factor = &edge.factor * &mink;
                    done.set(true);
                }
            })
        };
        let attempt = uncontracted::with_insertion_hook(hook, || {
            root_attempt(&case, &vakint, Some(case.provenance.clone()), false)
        })?;
        ensure!(done.get(), "(b): no propagator mutated");
        expect_error(attempt, "(b) mink on a propagator")?;
        Ok(())
    })
}

/// The coordinator's decision (21 Nov, on beta's observation): a K5 variable
/// must not affect a run in which `surrogate` is not enabled. At start-up the
/// hole-cut variable is read and resolved ONLY with the kind enabled
/// (`k5::start_up_hole_cut`, the one call in uncontracted.rs `ENVIRONMENT`):
/// - a garbage value with `check:root,child,ir-hard` (and with the switch
///   off): the variable is not even read, nothing is resolved, the run
///   proceeds, the identity is the base's;
/// - a garbage value with `surrogate` listed: a start-up refusal naming the
///   variable, and nothing stored.
#[test]
fn k5_hole_cut_variable_is_ignored_without_the_kind() -> Result<()> {
    let before = k5::resolved_hole_cut();
    for value in ["off", "check:root,child,ir-hard", "on", "check"] {
        let switch = Switch::parse(Some(value)).map_err(|error| eyre!(error))?;
        let read = std::cell::Cell::new(false);
        let mode = k5::start_up_hole_cut(switch, || {
            read.set(true);
            Ok(Some("froms".to_owned()))
        });
        ensure!(mode == Ok(None), "{value}: {mode:?}");
        ensure!(
            !read.get(),
            "{value}: the variable was read without the kind"
        );
    }
    let default = Switch::parse(Some("check:root,child,ir-hard")).map_err(|error| eyre!(error))?;
    ensure!(
        with_switch(default, uncontracted_ledger::current_identity).as_deref()
            == Some("check:root,child,ir-hard:strict"),
        "the identity of a run without the kind changed"
    );
    for value in ["check:surrogate", "on:root,child,ir-hard,surrogate"] {
        let switch = Switch::parse(Some(value)).map_err(|error| eyre!(error))?;
        let refused = k5::start_up_hole_cut(switch, || Ok(Some("froms".to_owned())));
        ensure!(
            refused
                .as_ref()
                .is_err_and(|error| error.contains(k5::HOLE_CUT_VARIABLE)),
            "{value}: a garbage value was not refused: {refused:?}"
        );
    }
    ensure!(
        k5::resolved_hole_cut() == before,
        "a refused value was stored"
    );
    Ok(())
}

/// The direct record's record class (`Attempt::shape`, read by beta's
/// setting-caused-decline rule, uncontracted_ledger.rs `setting_wide_classes`
/// and RunTally): `"direct"`, its OWN class. `K1` means "a remainder of an
/// evaluation root" (uncontracted.rs `shape_of`), `K2` "a remainder of a
/// proper UV child", `K3` "an IR hard term"; the <= 4-loop surrogate record is
/// none of these (no R* plan, the format's `kind direct`), and the codes
/// that can decline it differ (no `open_root`). Labelled `K1`, compared
/// direct records would MASK a setting-caused decline of every surrogate
/// K1 root in the same diagram-kind (or run tally). Shown here: (a) with
/// `"direct"` such a decline fails the diagram; (b) the same records with the
/// direct ones labelled `K1` would pass; (c) the direct record of L2 17 is
/// counted in the class `direct` (compared).
#[test]
fn k5_direct_record_is_its_own_record_class() -> Result<()> {
    test_initialise()?;
    let close = |direct_shape: &str| {
        with_switch(check_surrogate(), || -> Result<bool> {
            let ledger =
                DiagramLedger::open("k5-shape", "000001").ok_or_else(|| eyre!("no ledger"))?;
            for _ in 0..3 {
                uncontracted_ledger::record_shaped(
                    Kind::Surrogate,
                    "K1",
                    Some("open_root"),
                    Mode::Check,
                );
            }
            uncontracted_ledger::record_shaped(Kind::Surrogate, direct_shape, None, Mode::Check);
            Ok(ledger.close(true).is_err())
        })
    };
    ensure!(
        close("direct")?,
        "(a) a setting-caused decline of every K1 root passed"
    );
    ensure!(
        !close("K1")?,
        "(b) the counterfactual did not mask (the test proves nothing)"
    );
    let vakint = Vakint::new()?;
    let original = colour_factored(&load(L2_17_DOT)?)?;
    let evaluator = PropagatorInsertionEvaluator::new(settings());
    let evidence = with_switch(check_surrogate(), || -> Result<serde_json::Value> {
        let ledger =
            DiagramLedger::open(original.name(), "000017").ok_or_else(|| eyre!("no ledger"))?;
        let value = evaluator.pole_part_closing_root(&original, &vakint, true);
        ledger
            .close(value.is_ok())?
            .ok_or_else(|| eyre!("no evidence"))
    })?;
    let shapes = &evidence["kinds"]["surrogate"]["shapes"];
    ensure!(
        shapes["direct"]["compared"] == serde_json::json!(1)
            && shapes.as_object().is_some_and(|classes| classes.len() == 1),
        "(c) the direct record's class: {shapes}"
    );
    Ok(())
}

/// Beta's guard on the DIRECT record (the coordinator's decision on beta's
/// suggestion (i)): `build_direct` calls `declared_momenta` after its
/// `from_text`, as `build` does for K1 / K2. The k5-direct case (L2 17), its
/// coefficient `-1/2 gs^2 Q(3).Q(3)` changed to `-1/2 gs^2 Q(3).Q(99)` (the
/// same degree; `Q(99)` is no edge of the record, so vakint's reader takes it
/// as a constant vector and no K5 rule sees it): our attempt is the ERROR
/// `undeclared_momentum` (never a decline) with no value, and the hook
/// fails; unchanged, the same case compares equal
/// (k5_direct_record_and_value_equal_the_fixture).
#[test]
fn k5_direct_record_with_an_undeclared_momentum_is_an_error() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let prepared = prepare(L2_17_DOT, &vakint)?;
        let text = K5_DIRECT
            .lines()
            .map(|line| {
                if line.starts_with("insertion_factor\t") {
                    let at = line.rfind("::Q(3,").expect("the coefficient's second Q(3)");
                    format!("{}::Q(99,{}", &line[..at], &line[at + "::Q(3,".len()..])
                } else {
                    line.to_owned()
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        ensure!(text.contains("::Q(99,"), "no probe in the coefficient");
        let fixture = record(&text)?;
        let diagram = with_branch(&prepared.surrogate, &fixture)?;
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let (value, captured, _) = direct_run(
            &graph,
            &diagram,
            Some(&prepared.provenance),
            &vakint,
            check_surrogate(),
            true,
        );
        let attempt = captured
            .iter()
            .find(|attempt| attempt.kind == Kind::Surrogate)
            .ok_or_else(|| eyre!("no surrogate attempt: {captured:?}"))?;
        ensure!(
            attempt.outcome == "error:undeclared_momentum",
            "not the guard's error: {} {:?}",
            attempt.outcome,
            attempt.detail
        );
        ensure!(attempt.value.is_none(), "an erroneous record has a value");
        ensure!(value.is_err(), "the direct hook passed an erroneous record");
        Ok(())
    })
}

/// The last rebase's item (3) (r4): `build_direct`, the second writer, runs
/// the pole check of every other record (`record_eps_pole`, as `build` does
/// since beta's 38c71f377d) BEFORE `to_text`. The k5-direct case (L2 17) with
/// the factor of its plain chain vertex 3 (`scalar_insertion_1`, factor 1)
/// changed to `1 + ε^-1`: a pole outside the insertions, which the gated
/// reader refuses as `NotScalar` (not a designed row: an error). Our attempt
/// is the DESIGNED decline `factor_eps_pole` (in the closed list), decided
/// before the record is written: no text, no value, the detail names the
/// field. Unchanged, the same case compares equal
/// (k5_direct_record_and_value_equal_the_fixture).
#[test]
fn k5_direct_record_with_an_eps_pole_outside_the_insertions_is_declined() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let prepared = prepare(L2_17_DOT, &vakint)?;
        let mut fixture = record(K5_DIRECT)?;
        let eps = Atom::var(GS.dim_epsilon);
        let vertex = fixture
            .vertices
            .iter_mut()
            .find(|vertex| vertex.id == 3 && vertex.rule == "scalar_insertion_1")
            .ok_or_else(|| eyre!("no plain chain vertex 3 in k5-direct"))?;
        vertex.factor = Atom::num(1) + eps.pow(Atom::num(-1));
        let diagram = with_branch(&prepared.surrogate, &fixture)?;
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let (_, captured, _) = direct_run(
            &graph,
            &diagram,
            Some(&prepared.provenance),
            &vakint,
            check_surrogate(),
            true,
        );
        let attempt = captured
            .iter()
            .find(|attempt| attempt.kind == Kind::Surrogate)
            .ok_or_else(|| eyre!("no surrogate attempt: {captured:?}"))?;
        ensure!(
            attempt.outcome == "declined:factor_eps_pole",
            "not the designed pole decline: {} {:?}",
            attempt.outcome,
            attempt.detail
        );
        ensure!(
            uncontracted::is_designed_code("factor_eps_pole"),
            "factor_eps_pole is not in the closed list"
        );
        ensure!(
            attempt
                .detail
                .as_deref()
                .is_some_and(|detail| detail.starts_with("vertex 3 factor: a pole of order 1")),
            "the detail does not name the field: {:?}",
            attempt.detail
        );
        ensure!(
            attempt.text.is_none() && attempt.value.is_none(),
            "a declined record was written or valued"
        );
        Ok(())
    })
}

/// Beta's suggestion (ii), the coordinator's decision: C_H is written in the
/// RECORD's name of the anchor's edge. When e_in is a paired crown edge of the
/// record whose sink is the hole's anchor, the record's line there is the
/// incoming leg `100000 + e_in` (uncontracted.rs `paired_leg_id`), and C_H must
/// name that leg (the same rename `crown_legs_at_node` makes in a ct: Delta),
/// or beta's guard refuses it (`insertion_momentum`). Constructed here at the
/// two functions `finish_with` calls (no graph of our sets has the case: the
/// dump scan finds 0 such anchors): (a) the anchor's leg of a paired crown
/// edge e_in = 5 is 100005, of a plain e_in its own id, any other owner an
/// error; (b) C_H in Q(e_in) becomes Q(100005), C_H in Q(e_out) becomes
/// s Q(100005), and with a plain leg the old rewrite is unchanged.
#[test]
fn k5_coefficient_names_the_anchor_leg_of_a_paired_crown_edge() -> Result<()> {
    use vakint::p_integral::factorized::{Endpoint, FactorizedEdge, LineKind};
    test_initialise()?;
    let edge = |id, source, target, source_hedge, target_hedge| FactorizedEdge {
        id,
        source,
        target,
        source_hedge,
        target_hedge,
        particle: "g".into(),
        line: LineKind::Vector,
        power: 1,
        mass_squared: Atom::Zero,
        momentum: vec![0, 1],
        factor: Atom::num(1),
    };
    // anchor vertex 2, its inner half-edge 9; e_in = 5 a paired crown edge:
    // outgoing leg 5 at vertex 7 (half-edge 8), incoming leg 100005 at 2 (9)
    let paired = vec![
        edge(5, Endpoint::Vertex(7), Endpoint::External, Some(8), None),
        edge(
            100005,
            Endpoint::External,
            Endpoint::Vertex(2),
            None,
            Some(9),
        ),
        edge(
            6,
            Endpoint::Vertex(2),
            Endpoint::Vertex(3),
            Some(10),
            Some(11),
        ),
    ];
    let plain = vec![
        edge(5, Endpoint::External, Endpoint::Vertex(2), None, Some(9)),
        edge(
            6,
            Endpoint::Vertex(2),
            Endpoint::Vertex(3),
            Some(10),
            Some(11),
        ),
    ];
    ensure!(
        k5::anchor_leg(&paired, 9, 5).ok() == Some(100005),
        "(a) paired"
    );
    ensure!(k5::anchor_leg(&plain, 9, 5).ok() == Some(5), "(a) plain");
    match k5::anchor_leg(&plain, 10, 5) {
        Err(decline) => ensure!(decline.code == "surrogate_chain" && !decline.designed),
        Ok(id) => return Err(eyre!("(a) a foreign owner {id} was accepted")),
    }
    let q = |e: usize| {
        symbolica::try_parse!(&format!(
            "gammalooprs::{{spenso::rank1,spenso::tensor,spenso::component-labels:subscript,spenso::tensor-label:q}}::Q({e},spenso::{{spenso::representation,spenso::self_dual}}::mink(4,gammalooprs::{{spenso::index}}::hedge(9,1)))"
        ))
        .map_err(|error| eyre!("{error}"))
    };
    // e_in = 5 enters the chain at the anchor (its sink), e_out = 12 leaves it
    let inner = [(5, true), (12, false)];
    ensure!(
        k5::in_anchor_edge(q(5)?, inner, 100005) == q(100005)?,
        "(b) Q(e_in)"
    );
    ensure!(
        k5::in_anchor_edge(q(12)?, inner, 100005) == q(100005)?,
        "(b) Q(e_out), s = +1"
    );
    ensure!(
        k5::in_anchor_edge(q(12)?, inner, 5) == q(5)?,
        "(b) plain leg: the old rewrite"
    );
    ensure!(
        k5::in_anchor_edge(q(5)?, inner, 5) == q(5)?,
        "(b) plain leg: Q(e_in) kept"
    );
    Ok(())
}
