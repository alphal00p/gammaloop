//! Oracle tests of the uncontracted K2 records (UV children,
//! `GAMMALOOP_UNCONTRACTED_NUMERATOR` kind `child`), written independently of
//! `uncontracted_tests`; the helpers are those of
//! [`super::uncontracted_oracle_tests`].
//!
//! The oracle, as for K1: the record the switch emits for a real child
//! remainder is read back with vakint's `from_text`, evaluated here with
//! `evaluate_factorized_closure` at the record's depth and finished with
//! `finish_tensor`; it must equal the contracted path's `value.txt` through
//! its `valid_through` order by expanded difference, with the isotropic
//! channel markers compared BY NAME (no relabelling). The emitted text must
//! equal the fixture `input.txt` field by field.
//!
//! Every graph run is a check run with `GAMMALOOP_UNCONTRACTED_STRICT` on, and
//! its `UNCONTRACTED_SUMMARY` lines are checked against the captured attempts:
//! attempted = captured, compared = attempted - declined, equal = compared,
//! no error, every decline code one the test lists. A test that reads data
//! from outside the tree is `#[ignore]` and fails when the data is missing.
//!
//! Data: the beta's k2-child, k2-child-ins and k1-root records as already
//! vendored in `testdata/u7a/` (byte-identical to the verified fixtures; the
//! sha256 are in `testdata/u7a/oracle/SHA256SUMS.txt`).

use std::collections::{BTreeMap, BTreeSet};

use eyre::{Result, ensure, eyre};
use idenso::{representations::Bispinor, tensor::SymbolicTensor};
use linnet::half_edge::subgraph::SuBitGraph;
use spenso::structure::{
    representation::{Minkowski, RepName},
    slot::{DummyAind, ParseableAind},
};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};
use vakint::{
    Vakint,
    p_integral::factorized::{FactorizedKind, FactorizedNumerator},
};

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    graph::Graph,
    initialisation::test_initialise,
    numerator::{aind::Aind, symbolica_ext::NumeratorAtomExt},
    utils::GS,
    uv::approx::streamed::with_production_piece_limit,
};

use super::{
    REFERENCE_MASS_SQUARED, RStarPlan,
    boundary::{BoundaryMode, with_boundary_mode},
    expression::RStarIntegrand,
    plan::{join_edge_ids, uv_child_edge_ids},
    recursive_tests::{colour_factored, load, settings},
    uncontracted::{
        self, Captured, ChildClosure, ClosureData, Kind, Mode, Outcome, Site, Switch, with_capture,
        with_strict, with_switch,
    },
    uncontracted_oracle_tests::{
        assume_defaults, atoms_equal, auxiliary_mass, branch, equal_through, fixture_value,
        fixtures_dir, key, metric_root, oracle_value, parse, read, record, record_differences,
        valid_through,
    },
};

const L2_22_DOT: &str = include_str!("testdata/u7a/l2-d22.dot");
const L3_36_DOT: &str = include_str!("testdata/u7a/l3-d036.dot");
const K1_ROOT: &str = include_str!("testdata/u7a/k1-root.input.txt");
const K2_CHILD: &str = include_str!("testdata/u7a/k2-child.input.txt");
const K2_CHILD_VALUE: &str = include_str!("testdata/u7a/k2-child.value.txt");
const K2_CHILD_INS: &str = include_str!("testdata/u7a/k2-child-ins.input.txt");
const K2_CHILD_INS_VALUE: &str = include_str!("testdata/u7a/k2-child-ins.value.txt");
const K2_CHILD_INS_BRANCH: &str = include_str!("testdata/u7a/k2-child-ins.branch.txt");

fn root_and_child() -> Switch {
    Switch::new(Mode::Check, &[Kind::Root, Kind::Child])
}

fn on_root_and_child() -> Switch {
    Switch::new(Mode::On, &[Kind::Root, Kind::Child])
}

// ---------------------------------------------------------------------------
// Branches
// ---------------------------------------------------------------------------

/// The colour-factored DOT with the given vertex and internal-edge factors
/// (the branch's own, as the record or the site dump states them).
fn branch_from(
    dot: &str,
    vertices: &BTreeMap<usize, Atom>,
    edges: &BTreeMap<usize, Atom>,
) -> Result<feynkit_graph::FeynmanDiagram> {
    let diagram = colour_factored(&load(dot)?)?;
    let vertex_fragments = diagram
        .vertices()
        .filter_map(|(id, _)| vertices.get(&id.0).map(|f| (id, f.clone())))
        .map(|(id, f)| Ok((id, SymbolicTensor::infer(f)?)))
        .collect::<Result<Vec<_>>>()?;
    let edge_fragments = diagram
        .edges()
        .filter_map(|(id, _, _)| edges.get(&id.0).map(|f| (id, f.clone())))
        .map(|(id, f)| Ok((id, SymbolicTensor::infer(f)?)))
        .collect::<Result<Vec<_>>>()?;
    diagram
        .with_numerator_fragments(vertex_fragments, edge_fragments)
        .map_err(|error| eyre!("{error}"))
}

/// The vertex and internal-edge factors stated by records (`ct:` vertices
/// skipped), in the DOT's mink(4).
fn factors_of(records: &[&FactorizedNumerator]) -> (BTreeMap<usize, Atom>, BTreeMap<usize, Atom>) {
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for value in records {
        for vertex in value.vertices.iter().filter(|v| !v.rule.starts_with("ct:")) {
            vertices.insert(vertex.id, vertex.factor.map_mink_dim(Atom::num(4)));
        }
        for edge in value.edges.iter().filter(|e| !e.is_external()) {
            edges.insert(edge.id, edge.factor.map_mink_dim(Atom::num(4)));
        }
    }
    (vertices, edges)
}

/// The factors of a site dump (`edge` / `vertex` lines with `factor=`).
fn factors_of_site(text: &str) -> Result<(BTreeMap<usize, Atom>, BTreeMap<usize, Atom>)> {
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for line in text.lines() {
        let mut fields = line.split('\t');
        let (Some(what), Some(id)) = (fields.next(), fields.next()) else {
            continue;
        };
        let Some(factor) = fields.find_map(|f| f.strip_prefix("factor=")) else {
            continue;
        };
        if factor == "contracted" {
            continue;
        }
        let id = id.parse::<usize>()?;
        let factor = parse(factor)?.map_mink_dim(Atom::num(4));
        match what {
            "vertex" => vertices.insert(id, factor),
            "edge" => edges.insert(id, factor),
            _ => None,
        };
    }
    ensure!(
        !vertices.is_empty() && !edges.is_empty(),
        "no factors in the site dump"
    );
    Ok((vertices, edges))
}

// ---------------------------------------------------------------------------
// One strict check run and its summary
// ---------------------------------------------------------------------------

/// The counts of one `UNCONTRACTED_SUMMARY` line.
#[derive(Debug, Default)]
struct Summary {
    attempted: usize,
    uncontracted: usize,
    zero: usize,
    compared: usize,
    equal: usize,
    declined: BTreeMap<String, usize>,
    error: String,
    outcome: String,
}

fn summary(line: &str) -> Result<Summary> {
    let field = |name: &str| {
        line.split_whitespace()
            .find_map(|token| token.strip_prefix(&format!("{name}=")))
            .ok_or_else(|| eyre!("no {name}= in {line}"))
    };
    let count = |name: &str| -> Result<usize> { Ok(field(name)?.parse()?) };
    let mut declined = BTreeMap::new();
    let list = field("declined")?;
    if list != "none" {
        for entry in list.split(',') {
            let (code, n) = entry
                .rsplit_once(':')
                .ok_or_else(|| eyre!("bad declined entry {entry}"))?;
            declined.insert(code.to_owned(), n.parse()?);
        }
    }
    Ok(Summary {
        attempted: count("attempted")?,
        uncontracted: count("uncontracted")?,
        zero: count("zero")?,
        compared: count("compared")?,
        equal: count("equal")?,
        declined,
        error: field("error")?.to_owned(),
        outcome: field("outcome")?.to_owned(),
    })
}

/// One diagram root (the beta's closed metric root) under `switch` with
/// STRICT on and the boundary path on: its value, every attempt and the
/// summary lines it printed.
fn strict_run(
    vakint: &Vakint,
    diagram: &feynkit_graph::FeynmanDiagram,
    switch: Switch,
) -> Result<(Atom, Vec<Captured>, Vec<String>)> {
    let uv = settings();
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    // a fresh plan and Δ cache per run: a cache filled by another mode would
    // serve the children of this run and hide their records
    let plan = RStarPlan::with_automatic_edge(&graph, model, &uv)?;
    let closed_root = metric_root(&graph);
    let callback = |tensor: &Atom| -> Result<Atom> {
        Ok(tensor.replace_multiple(&closed_root.external_routing) * &closed_root.factor)
    };
    let _ = uncontracted::take_summaries();
    let (value, captured) = with_switch(switch, || {
        with_strict(true, || {
            with_boundary_mode(BoundaryMode::On, false, || {
                with_capture(|| {
                    plan.projected_counterterm_closing_root_scalar(
                        &graph,
                        model,
                        vakint,
                        &uv,
                        &callback,
                        &closed_root,
                    )
                })
            })
        })
    });
    let summaries = uncontracted::take_summaries();
    Ok((value?, captured, summaries))
}

/// The accounting of one check run: the summary lines agree with the
/// captured attempts; nothing is an error; every decline has a listed code;
/// every compared record is equal. Returns (root, child) attempt counts.
fn assert_accounting(
    what: &str,
    captured: &[Captured],
    summaries: &[String],
    allowed: &[&str],
) -> Result<(usize, usize)> {
    ensure!(
        !summaries.is_empty(),
        "{what}: no UNCONTRACTED_SUMMARY line"
    );
    let mut total = Summary::default();
    for line in summaries {
        let s = summary(line)?;
        ensure!(s.error == "none", "{what}: errors in {line}");
        ensure!(s.outcome == "ok", "{what}: outcome in {line}");
        let declined = s.declined.values().sum::<usize>();
        ensure!(
            s.attempted == s.uncontracted + s.zero + declined,
            "{what}: attempted != uncontracted + zero + declined in {line}"
        );
        ensure!(
            s.compared == s.attempted - declined,
            "{what}: compared != attempted - declined in {line}"
        );
        ensure!(
            s.equal == s.compared,
            "{what}: not every compared record is equal: {line}"
        );
        for code in s.declined.keys() {
            ensure!(
                allowed.contains(&code.as_str()),
                "{what}: decline {code} is not one this test expects: {line}"
            );
        }
        total.attempted += s.attempted;
        total.compared += s.compared;
    }
    ensure!(
        total.attempted == captured.len(),
        "{what}: the summaries count {} attempts, {} were reported",
        total.attempted,
        captured.len()
    );
    let mut compared = 0;
    for attempt in captured {
        let outcome = attempt.outcome.as_str();
        let declined = outcome
            .strip_prefix("declined:")
            .is_some_and(|code| allowed.contains(&code));
        ensure!(
            outcome == "uncontracted" || outcome == "zero" || declined,
            "{what}: {} id {} outcome {outcome}",
            attempt.kind.as_str(),
            attempt.id
        );
        if attempt.equal.is_some() {
            ensure!(
                attempt.equal == Some(true),
                "{what}: id {} not equal",
                attempt.id
            );
            compared += 1;
        }
    }
    ensure!(
        compared == total.compared,
        "{what}: {compared} captured comparisons, the summaries say {}",
        total.compared
    );
    let roots = captured.iter().filter(|a| a.kind == Kind::Root).count();
    let children = captured.iter().filter(|a| a.kind == Kind::Child).count();
    Ok((roots, children))
}

/// One fixture record of `kind` against the captured attempts: the round
/// trip (fixed point of vakint's writer, equal to `input.txt` field by field),
/// the oracle value and the switch's value through `valid_through` (markers
/// by name), today's value too, the check equal.
fn compare_fixture(
    what: &str,
    kind: Kind,
    vakint: &Vakint,
    captured: &[Captured],
    input: &str,
    value_text: &str,
) -> Result<FactorizedNumerator> {
    let theirs = record(input)?;
    let wanted = key(&theirs);
    let mut found = Vec::new();
    for attempt in captured.iter().filter(|a| a.kind == kind) {
        let Some(text) = &attempt.text else { continue };
        let ours = record(text)?;
        if key(&ours) == wanted && ours.external_edges == theirs.external_edges {
            found.push((attempt, ours));
        }
    }
    ensure!(
        found.len() == 1,
        "{what}: {} captured {} records for {wanted:?} (outcomes {:?})",
        found.len(),
        kind.as_str(),
        captured
            .iter()
            .map(|a| format!("{}:{}:{}", a.kind.as_str(), a.id, a.outcome))
            .collect::<Vec<_>>()
    );
    let (attempt, ours) = found.pop().unwrap();
    let text = attempt.text.as_ref().unwrap();
    ensure!(
        attempt.outcome == "uncontracted",
        "{what}: outcome {}",
        attempt.outcome
    );
    ensure!(
        ours.to_text().map_err(|e| eyre!("{e}"))? == *text,
        "{what}: the emitted text is not a fixed point of vakint's writer"
    );
    let differences = record_differences(&ours, &theirs);
    ensure!(
        differences.is_empty(),
        "{what}: the emitted record differs from input.txt: {differences:#?}\n{text}"
    );
    let expected = fixture_value(value_text)?;
    let order = valid_through(value_text)?;
    let oracle = oracle_value(vakint, text, attempt.id)?;
    ensure!(
        equal_through(&oracle, &expected, order)?,
        "{what}: the oracle value differs from value.txt through eps^{order}\noracle {oracle}\nvalue  {expected}"
    );
    let value = attempt
        .value
        .as_ref()
        .ok_or_else(|| eyre!("{what}: no value"))?;
    ensure!(
        equal_through(value, &expected, order)?,
        "{what}: the switch's value differs from value.txt through eps^{order}"
    );
    let contracted = attempt
        .contracted
        .as_ref()
        .ok_or_else(|| eyre!("{what}: not compared (no check)"))?;
    ensure!(
        equal_through(contracted, &expected, order)?,
        "{what}: today's value differs from value.txt (not this branch)"
    );
    ensure!(
        attempt.equal == Some(true),
        "{what}: check equal={:?}",
        attempt.equal
    );
    eprintln!(
        "U7A_ORACLE_K2 {what} id={} insertions={} identical={:?} value_zero={}",
        attempt.id,
        ours.insertions.len(),
        attempt.identical,
        expected.is_zero()
    );
    Ok(ours)
}

/// A diagram in check (strict, accounted) and on (strict) mode: on == off.
fn check_and_on(
    what: &str,
    vakint: &Vakint,
    diagram: &feynkit_graph::FeynmanDiagram,
    allowed: &[&str],
) -> Result<Vec<Captured>> {
    let (off, nothing, _) = strict_run(vakint, diagram, Switch::OFF)?;
    ensure!(nothing.is_empty(), "{what}: off reports nothing");
    let (checked, captured, summaries) = strict_run(vakint, diagram, root_and_child())?;
    ensure!(checked == off, "{what}: check returns today's value");
    let (roots, children) = assert_accounting(what, &captured, &summaries, allowed)?;
    ensure!(children > 0, "{what}: no child record attempted");
    let (on, on_captured, on_summaries) = strict_run(vakint, diagram, on_root_and_child())?;
    ensure!(
        on_captured.iter().all(|a| a.outcome == "uncontracted"
            || a.outcome == "zero"
            || a.outcome
                .strip_prefix("declined:")
                .is_some_and(|c| allowed.contains(&c))),
        "{what}: on mode outcomes {:?}",
        on_captured
            .iter()
            .map(|a| a.outcome.clone())
            .collect::<Vec<_>>()
    );
    ensure!(
        on_summaries.iter().all(|l| l.contains(" error=none ")),
        "{what}: on mode errors {on_summaries:?}"
    );
    ensure!(
        atoms_equal(&on, &off),
        "{what}: on differs from off\noff {off}\non  {on}"
    );
    eprintln!(
        "U7A_ORACLE_K2 {what}: {roots} root + {children} child records, all accounted, on == off"
    );
    Ok(captured)
}

// ---------------------------------------------------------------------------
// Graph-driven K2 fixtures
// ---------------------------------------------------------------------------

/// k2-child: the isotropic closed child {2,3,5} of L2 22 (one loop, bare,
/// degree 0), from the graph with the branch's factors (k1-root + k2-child
/// records): record == input.txt, oracle value == value.txt (markers by
/// name), check equal, every record accounted, on == off.
#[test]
fn k2_child_from_the_graph() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let root = record(K1_ROOT)?;
        let child = record(K2_CHILD)?;
        let (vertices, edges) = factors_of(&[&root, &child]);
        let diagram = branch_from(L2_22_DOT, &vertices, &edges)?;
        let vakint = Vakint::new()?;
        let captured = check_and_on("L2 22", &vakint, &diagram, &[])?;
        let ours = compare_fixture(
            "k2-child",
            Kind::Child,
            &vakint,
            &captured,
            K2_CHILD,
            K2_CHILD_VALUE,
        )?;
        let FactorizedKind::Closure(spec) = &ours.kind else {
            return Err(eyre!("k2-child: not a closure record"));
        };
        ensure!(
            spec.independent_externals,
            "k2-child: externals=independent"
        );
        // the K1 record of the same run still equals its fixture
        compare_fixture(
            "k1-root (same run)",
            Kind::Root,
            &vakint,
            &captured,
            K1_ROOT,
            include_str!("testdata/u7a/k1-root.value.txt"),
        )?;
        Ok(())
    })
}

/// k2-child-ins: the open child {5,6,7,8} of L3 36, remainder {5,6} with the
/// child {7,8} contracted (a poled `ct:` insertion), from the graph with the
/// site's factors: record == input.txt, oracle value == value.txt, check
/// equal, every record accounted, on == off.
#[test]
fn k2_child_ins_from_the_graph() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let (vertices, edges) = factors_of_site(K2_CHILD_INS_BRANCH)?;
        let diagram = branch_from(L3_36_DOT, &vertices, &edges)?;
        let vakint = Vakint::new()?;
        let captured = check_and_on("L3 36", &vakint, &diagram, &[])?;
        let ours = compare_fixture(
            "k2-child-ins",
            Kind::Child,
            &vakint,
            &captured,
            K2_CHILD_INS,
            K2_CHILD_INS_VALUE,
        )?;
        ensure!(ours.insertions.len() == 1, "k2-child-ins: one insertion");
        ensure!(
            ours.insertion_pole_order(GS.dim_epsilon)
                .map_err(|e| eyre!("{e}"))?
                == 1,
            "k2-child-ins: the insertion's pole order is 1"
        );
        Ok(())
    })
}

/// Every K1 and K2 record of the k1-quick diagrams L4 40, 42 and 49 (the
/// branch's factors from each bare record): check strict and accounted
/// (L4-49's one-loop root remainder is the designed `one_loop_tree_cut`),
/// on == off, the 21 K1 fixtures still equal. Child records with two
/// insertions are counted and their pole order checked to be the sum.
#[test]
#[ignore = "reads k1-quick from GAMMALOOP_U7A_FIXTURES_DIR (about 15 min single core)"]
fn k2_children_of_l4_40_42_49() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let root = fixtures_dir()?.join("k1-quick");
        let mut fixtures = 0;
        let mut children_total = 0;
        let mut two_insertions = 0;
        let mut two_nonzero = 0;
        for (diagram, first) in [("L4-40", "rr"), ("L4-42", "rr"), ("L4-49", "")] {
            let base = root.join(diagram);
            let number = diagram.trim_start_matches("L4-");
            let dir = |i: usize| base.join(format!("{first}{i}"));
            let dot = read(&dir(1).join(format!("raw/diagram-0000{number}.dot")))?;
            let bare = record(&read(&dir(1).join("input.txt"))?)?;
            let branch = branch(&dot, &bare)?;
            let vakint = Vakint::new()?;
            let allowed: &[&str] = if diagram == "L4-49" {
                &["one_loop_tree_cut"]
            } else {
                &[]
            };
            let captured = check_and_on(diagram, &vakint, &branch, allowed)?;
            for i in 1..=7 {
                compare_fixture(
                    &format!("{diagram} {i}"),
                    Kind::Root,
                    &vakint,
                    &captured,
                    &read(&dir(i).join("input.txt"))?,
                    &read(&dir(i).join("value.txt"))?,
                )?;
                fixtures += 1;
            }
            for attempt in captured.iter().filter(|a| a.kind == Kind::Child) {
                children_total += 1;
                let Some(text) = &attempt.text else { continue };
                let value = record(text)?;
                let FactorizedKind::Closure(spec) = &value.kind else {
                    return Err(eyre!("{diagram}: a child record is not a closure"));
                };
                ensure!(
                    spec.independent_externals,
                    "{diagram}: child externals=independent"
                );
                let (fermion, ghost) = spec.original_cycles.ok_or_else(|| eyre!("no cycles"))?;
                let carried = value
                    .insertions
                    .iter()
                    .fold((0, 0), |(f, g), i| (f + i.cycles.0, g + i.cycles.1));
                ensure!(
                    carried.0 <= fermion && carried.1 <= ghost,
                    "{diagram}: child insertion cycles above original_cycles"
                );
                if value.insertions.len() >= 2 {
                    two_insertions += 1;
                    let mut sum = 0;
                    for insertion in &value.insertions {
                        let mut alone = value.clone();
                        alone.insertions = vec![insertion.clone()];
                        sum += alone
                            .insertion_pole_order(GS.dim_epsilon)
                            .map_err(|e| eyre!("{e}"))?;
                    }
                    let order = value
                        .insertion_pole_order(GS.dim_epsilon)
                        .map_err(|e| eyre!("{e}"))?;
                    ensure!(
                        order == sum,
                        "{diagram}: pole order {order}, the sum is {sum}"
                    );
                    let nonzero = attempt.value.as_ref().is_some_and(|v| !v.is_zero());
                    two_nonzero += usize::from(nonzero);
                    eprintln!(
                        "U7A_ORACLE_K2 {diagram} child id={} insertions={} pole_order={order} outcome={} equal={:?} value_zero={}",
                        attempt.id,
                        value.insertions.len(),
                        attempt.outcome,
                        attempt.equal,
                        !nonzero
                    );
                }
            }
        }
        ensure!(fixtures == 21, "{fixtures} K1 fixtures compared");
        ensure!(
            two_nonzero >= 1,
            "no child record with two poled insertions and a nonzero value ({two_insertions} with two)"
        );
        eprintln!(
            "U7A_ORACLE_K2 L4 40/42/49: {fixtures} K1 fixtures, {children_total} child records, {two_insertions} with two insertions ({two_nonzero} nonzero)"
        );
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// Site-driven: the k2-child-ins remainder through the emitter's entry
// ---------------------------------------------------------------------------

/// L3 36 with the site's factors, the child {5,6,7,8}'s own plan and its term
/// {5,6} with {7,8} contracted.
struct ChildCase {
    graph: Graph,
    diagram: feynkit_graph::FeynmanDiagram,
    support: SuBitGraph,
    term: usize,
    nested: SuBitGraph,
    delta: Atom,
    fixture: FactorizedNumerator,
}

fn child_case() -> Result<ChildCase> {
    let (vertices, edges) = factors_of_site(K2_CHILD_INS_BRANCH)?;
    let diagram = branch_from(L3_36_DOT, &vertices, &edges)?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let uv = settings();
    let root_plan = RStarPlan::with_automatic_edge(&graph, diagram.model(), &uv)?;
    let support = root_plan
        .counterterms
        .iter()
        .map(|call| call.argument.support.clone())
        .find(|support| uv_child_edge_ids(&graph, support) == [5, 6, 7, 8])
        .ok_or_else(|| eyre!("L3 36 has no child {{5,6,7,8}}"))?;
    let plan = RStarPlan::for_support(&graph, diagram.model(), &support, &uv)?;
    let term = plan
        .ordinary_uv
        .iter()
        .position(|term| {
            join_edge_ids(&uv_child_edge_ids(&graph, &term.remainder.remaining)) == "5,6"
                && term.remainder.insertions.len() == 1
                && uv_child_edge_ids(&graph, &term.remainder.insertions[0].support) == [7, 8]
        })
        .ok_or_else(|| eyre!("the child plan has no term {{5,6}} with {{7,8}} contracted"))?;
    let nested = plan.ordinary_uv[term].remainder.insertions[0]
        .support
        .clone();
    let fixture = record(K2_CHILD_INS)?;
    let insertion = fixture
        .insertions
        .first()
        .ok_or_else(|| eyre!("k2-child-ins has no insertion"))?;
    let delta = &insertion.coefficient * &insertion.factor;
    Ok(ChildCase {
        graph,
        diagram,
        support,
        term,
        nested,
        delta,
        fixture,
    })
}

/// The emitter's attempt on the child remainder in on mode, with the nested
/// child counterterm `delta`, a closure (`None`: the open child) and an
/// optional Taylor order.
fn child_attempt(
    case: &ChildCase,
    delta: &Atom,
    closure: Option<(Atom, ClosureData)>,
    taylor_order: Option<i32>,
) -> Result<uncontracted::Attempt> {
    let vakint = Vakint::new()?;
    let uv = settings();
    let vakint_settings = uv.vakint.true_settings();
    let model = case.diagram.model();
    let mut plan = RStarPlan::for_support(&case.graph, model, &case.support, &uv)?;
    if let Some(order) = taylor_order {
        plan.taylor_order = order;
    }
    let term = &plan.ordinary_uv[case.term];
    let children = BTreeMap::from([(case.nested.clone(), delta.clone())]);
    let insertions = term
        .counterterms
        .iter()
        .map(|call| {
            children
                .get(&call.argument.support)
                .cloned()
                .ok_or_else(|| eyre!("no Delta for a counterterm support"))
        })
        .collect::<Result<Vec<_>>>()?;
    // the crown of the child (the evaluation's root edges)
    let boundary = plan
        .root_boundary
        .clone()
        .ok_or_else(|| eyre!("a child plan has a root boundary"))?;
    let root_edges = case
        .graph
        .iter_edges_of(&boundary)
        .map(|(_, edge, _)| edge)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let auxiliary = auxiliary_mass();
    let (root_closure, closure_data) = match &closure {
        Some((tensor, data)) => (Some(tensor), Some(data)),
        None => (None, None),
    };
    let site = Site {
        graph: &case.graph,
        model,
        plan: &plan,
        term,
        children: &children,
        root_closure,
        closed_root: None,
        closure_data,
        root_projector: closure.is_some(),
        root_edges: &root_edges,
        reference_mass_squared: *REFERENCE_MASS_SQUARED,
        auxiliary_mass_squared: &auxiliary,
        settings: &vakint_settings,
        vakint: &vakint,
        uv_settings: &uv,
        top_mass: None,
        epsilon: GS.dim_epsilon,
        today: "site-test".into(),
        ir_hard: None,
    };
    uncontracted::attempt(&site, Kind::Child, 1, Mode::On, || {
        RStarIntegrand::new_with_closure(
            &case.graph,
            model,
            &term.remainder,
            &insertions,
            plan.marked_edge,
            Atom::var(*REFERENCE_MASS_SQUARED),
            root_closure,
        )
    })
}

fn declined(attempt: &uncontracted::Attempt) -> Option<String> {
    match &attempt.outcome {
        Outcome::Declined(decline) => Some(decline.code.clone()),
        _ => None,
    }
}

/// The child site with the fixture's nested Delta: the record equals the
/// fixture and its value equals value.txt (no in-process child computed).
#[test]
fn site_child_with_the_fixture_delta_equals_its_value() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = child_case()?;
        let attempt = child_attempt(&case, &case.delta, None, None)?;
        let text = attempt
            .text
            .clone()
            .ok_or_else(|| eyre!("no record built: {:?}", attempt.outcome))?;
        let differences = record_differences(&record(&text)?, &case.fixture);
        ensure!(
            differences.is_empty(),
            "k2-child-ins site: {differences:#?}"
        );
        let value = attempt
            .value()
            .ok_or_else(|| eyre!("declined: {:?}", attempt.outcome))?;
        let expected = fixture_value(K2_CHILD_INS_VALUE)?;
        let order = valid_through(K2_CHILD_INS_VALUE)?;
        ensure!(
            equal_through(&value, &expected, order)?,
            "k2-child-ins site: value"
        );
        let oracle = oracle_value(&Vakint::new()?, &text, 1)?;
        ensure!(
            equal_through(&oracle, &expected, order)?,
            "k2-child-ins site: oracle value"
        );
        Ok(())
    })
}

/// A mass symbol (M_aux, then R) in the nested child's Delta: declined
/// `insertion_mass_symbol` (designed), no value, no record.
#[test]
fn site_child_mass_symbol_in_its_insertion_declines() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = child_case()?;
        for mass in [auxiliary_mass(), Atom::var(*REFERENCE_MASS_SQUARED)] {
            let delta = &case.delta * (Atom::num(1) + &mass);
            let attempt = child_attempt(&case, &delta, None, None)?;
            ensure!(
                declined(&attempt).as_deref() == Some("insertion_mass_symbol"),
                "mass {mass}: outcome {:?}",
                attempt.outcome
            );
            ensure!(
                attempt.value().is_none() && attempt.text.is_none(),
                "no value, no record"
            );
        }
        Ok(())
    })
}

/// A child Taylor degree above 2: declined `vakint:taylor_degree` (designed).
#[test]
fn site_child_taylor_degree_three_declines() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = child_case()?;
        let attempt = child_attempt(&case, &case.delta, None, Some(3))?;
        ensure!(
            declined(&attempt).as_deref() == Some("vakint:taylor_degree"),
            "degree 3: outcome {:?}",
            attempt.outcome
        );
        ensure!(attempt.value().is_none(), "no value");
        Ok(())
    })
}

/// A child closed by a spinor chain that is not a closed pair (two open
/// bispinor slots): declined `vakint:open_spinor` (designed), no value.
#[test]
fn site_child_open_spinor_chain_declines() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = child_case()?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        let open = |h: u16| spin.to_symbolic([Aind::Hedge(h, 1).to_atom()]);
        let inner = spin.to_symbolic([Aind::new_dummy().to_atom()]);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let chain = function!(idenso::dirac::AGS.gamma, open(9001), &inner, &mu)
            * function!(idenso::dirac::AGS.gamma, &inner, open(9002), &mu);
        let data = ClosureData {
            kind: ChildClosure::DiracVertex,
            normalization: Atom::num(1),
            open: None,
        };
        let attempt = child_attempt(&case, &case.delta, Some((chain, data)), None)?;
        eprintln!("U7A_ORACLE_K2 open spinor chain: {:?}", attempt.outcome);
        ensure!(
            declined(&attempt).as_deref() == Some("vakint:open_spinor"),
            "open spinor chain: outcome {:?}",
            attempt.outcome
        );
        ensure!(attempt.value().is_none(), "no value");
        Ok(())
    })
}

/// A malformed nested Delta (a metric on half-edges that are not at its
/// vertex) is a bug, not a designed refusal: an ERROR outcome, no value.
#[test]
fn site_child_malformed_insertion_is_an_error() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = child_case()?;
        let stray = parse(
            "spenso::{symmetric,linear,real,spenso::tensor}::g(spenso::{spenso::representation,spenso::self_dual}::mink(gammalooprs::{}::dim,gammalooprs::{spenso::index}::hedge(9100,1)),spenso::{spenso::representation,spenso::self_dual}::mink(gammalooprs::{}::dim,gammalooprs::{spenso::index}::hedge(9101,1)))",
        )?;
        let attempt = child_attempt(&case, &(&case.delta * &stray), None, None)?;
        ensure!(
            attempt.is_error(),
            "malformed: outcome {:?}",
            attempt.outcome
        );
        ensure!(attempt.value().is_none(), "an error has no value");
        Ok(())
    })
}

/// A closed child whose closure comes without its closure data (the emitter
/// cannot state the normalisation): an ERROR outcome, never a silent
/// fallback.
#[test]
fn site_child_closure_without_its_data_is_an_error() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = child_case()?;
        let vakint = Vakint::new()?;
        let uv = settings();
        let vakint_settings = uv.vakint.true_settings();
        let model = case.diagram.model();
        let plan = RStarPlan::for_support(&case.graph, model, &case.support, &uv)?;
        let term = &plan.ordinary_uv[case.term];
        let children = BTreeMap::from([(case.nested.clone(), case.delta.clone())]);
        let boundary = plan.root_boundary.clone().unwrap();
        let root_edges = case
            .graph
            .iter_edges_of(&boundary)
            .map(|(_, edge, _)| edge)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let closure = Atom::num(1);
        let auxiliary = auxiliary_mass();
        let site = Site {
            graph: &case.graph,
            model,
            plan: &plan,
            term,
            children: &children,
            root_closure: Some(&closure),
            closed_root: None,
            closure_data: None,
            root_projector: true,
            root_edges: &root_edges,
            reference_mass_squared: *REFERENCE_MASS_SQUARED,
            auxiliary_mass_squared: &auxiliary,
            settings: &vakint_settings,
            vakint: &vakint,
            uv_settings: &uv,
            top_mass: None,
            epsilon: GS.dim_epsilon,
            today: "site-test".into(),
            ir_hard: None,
        };
        let insertions = vec![case.delta.clone()];
        let attempt = uncontracted::attempt(&site, Kind::Child, 1, Mode::On, || {
            RStarIntegrand::new_with_closure(
                &case.graph,
                model,
                &term.remainder,
                &insertions,
                plan.marked_edge,
                Atom::var(*REFERENCE_MASS_SQUARED),
                Some(&closure),
            )
        })?;
        ensure!(
            matches!(&attempt.outcome, Outcome::Error(d) if d.code == "closure_data_missing"),
            "outcome {:?}",
            attempt.outcome
        );
        Ok(())
    })
}
