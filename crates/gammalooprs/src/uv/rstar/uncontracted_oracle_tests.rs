//! Oracle tests of the uncontracted K1 records (`GAMMALOOP_UNCONTRACTED_NUMERATOR`,
//! [`super::uncontracted`]), written independently of `uncontracted_tests`.
//!
//! The oracle: the record the switch emits for a real remainder is read back
//! with vakint's `from_text`, evaluated here with `evaluate_factorized_closure`
//! and finished with `finish_tensor`; that value must equal the value the
//! contracted path gave for the same call (`value.txt`) through its
//! `valid_through` order, compared algebraically (expanded difference). The
//! emitted text must equal the fixture `input.txt` up to canonical order (every
//! line compared field by field, atoms by expanded difference).
//!
//! Three groups:
//! * graph-driven: a whole diagram root runs through R* in check, on and off
//!   mode; every captured root record is accounted for (no decline is
//!   allowed unless the test lists it);
//! * site-driven: one remainder of L4 diagram 40 is built through the
//!   emitter's own entry with chosen child counterterms, so that refusals and
//!   the momentum-conservation step are reached on a real record;
//! * the off mode: the switch off changes nothing.
//!
//! Data: `testdata/u7a/oracle/` (vendored, sha256 in its SHA256SUMS.txt). Tests that
//! read larger fixture sets take the directory from
//! `GAMMALOOP_U7A_FIXTURES_DIR` and are `#[ignore]`d.
//!
//! Settings: these tests go through R* with the defaults of this base, which
//! are the production defaults the fixtures were made with. Every switch
//! that changes the route of a root remainder must be unset; [`assume_defaults`]
//! fails otherwise.

use std::collections::{BTreeMap, BTreeSet};

use eyre::{Result, WrapErr, ensure, eyre};
use idenso::{representations::Bispinor, tensor::SymbolicTensor};
use linnet::half_edge::involution::{EdgeIndex, HedgePair};
use spenso::{
    network::tags::SPENSO_TAG,
    structure::{
        representation::{Minkowski, RepName},
        slot::{DummyAind, ParseableAind},
    },
};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    id::Replacement,
};
use vakint::{
    Vakint,
    p_integral::factorized::{FactorizedKind, FactorizedNumerator, evaluate_factorized_closure},
};

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    graph::Graph,
    initialisation::test_initialise,
    numerator::{aind::Aind, symbolica_ext::NumeratorAtomExt},
    utils::{GS, W_},
    uv::approx::streamed::with_production_piece_limit,
};

use super::{
    REFERENCE_MASS_SQUARED, RStarPlan,
    boundary::{BoundaryMode, ClosedRootNormalization, with_boundary_mode},
    expression::RStarIntegrand,
    plan::{join_edge_ids, uv_child_edge_ids},
    recursive_tests::{colour_factored, load, settings},
    top_mass::TopMassTerms,
    uncontracted::{self, Captured, Kind, Mode, Outcome, Site, Switch, with_capture, with_switch},
};

// ---------------------------------------------------------------------------
// Data and environment
// ---------------------------------------------------------------------------

const L4_40_DOT: &str = include_str!("testdata/u7a/oracle/l4-d40.dot");
const L4_40_RR1: &str = include_str!("testdata/u7a/oracle/l4-40-rr1.input.txt");
const L4_40_RR2: &str = include_str!("testdata/u7a/oracle/l4-40-rr2.input.txt");
const L4_40_RR2_VALUE: &str = include_str!("testdata/u7a/oracle/l4-40-rr2.value.txt");

/// The switches that change the route of a root remainder: the tests assume
/// the base's defaults (the production settings the fixtures were made with).
const ROUTE_SWITCHES: [&str; 12] = [
    "GAMMALOOP_UNCONTRACTED_NUMERATOR",
    "GAMMALOOP_UNCONTRACTED_DUMP",
    "GAMMALOOP_UV_CHILDREN",
    "GAMMALOOP_UV_CHILD_CLOSURE",
    "GAMMALOOP_UV_CHILD_TRUNCATE",
    "GAMMALOOP_UV_CHILD_EDGE_REMAINDER",
    "GAMMALOOP_REMAINDER_BOUNDARY",
    "GAMMALOOP_CHILD_BOUNDARY",
    "GAMMALOOP_ROOT_ROUTING",
    "GAMMALOOP_RSTAR_TERMS",
    "GAMMALOOP_EDGE_CONTRACTION",
    "PINT6_DUMP",
];

fn assume_defaults() -> Result<()> {
    let set = ROUTE_SWITCHES
        .iter()
        .filter(|name| std::env::var_os(name).is_some())
        .collect::<Vec<_>>();
    ensure!(
        set.is_empty(),
        "these tests assume the base's R* defaults; unset {set:?}"
    );
    Ok(())
}

/// The fixture directory of the ignored tests (`GAMMALOOP_U7A_FIXTURES_DIR`).
fn fixtures_dir() -> Result<std::path::PathBuf> {
    let dir = std::env::var("GAMMALOOP_U7A_FIXTURES_DIR").map_err(|_| {
        eyre!("set GAMMALOOP_U7A_FIXTURES_DIR to the wh-fixtures directory (k1-quick, k1-l5)")
    })?;
    let dir = std::path::PathBuf::from(dir);
    ensure!(dir.is_dir(), "{} is not a directory", dir.display());
    Ok(dir)
}

fn read(path: &std::path::Path) -> Result<String> {
    std::fs::read_to_string(path).wrap_err_with(|| format!("reading {}", path.display()))
}

fn parse(text: &str) -> Result<Atom> {
    symbolica::try_parse!(text).map_err(|error| eyre!("{error}: {text:.200}"))
}

fn auxiliary_mass() -> Atom {
    parse("gammalooprs::{scalar}::rstar_mass_squared").expect("the auxiliary mass symbol")
}

fn minkowski() -> Atom {
    Minkowski {}.new_rep(GS.dim).to_symbolic([])
}

fn record(text: &str) -> Result<FactorizedNumerator> {
    FactorizedNumerator::from_text(text).map_err(|error| eyre!("{error}"))
}

/// The `value` line of a value file.
fn fixture_value(text: &str) -> Result<Atom> {
    parse(
        text.lines()
            .find_map(|line| line.strip_prefix("value\t"))
            .ok_or_else(|| eyre!("no value line"))?,
    )
}

/// `n` of the value file's `valid_through eps^n` line (required).
fn valid_through(text: &str) -> Result<i64> {
    text.lines()
        .find_map(|line| line.strip_prefix("valid_through\t"))
        .and_then(|line| line.trim_start_matches("eps^").split_whitespace().next())
        .and_then(|order| order.parse().ok())
        .ok_or_else(|| eyre!("no valid_through line"))
}

/// `a == b` through ε^order: the expanded difference, truncated.
fn equal_through(a: &Atom, b: &Atom, order: i64) -> Result<bool> {
    let difference = (a - b).expand();
    if difference.is_zero() {
        return Ok(true);
    }
    let truncated = difference
        .series(GS.dim_epsilon, 0, order)?
        .to_atom()
        .expand();
    Ok(truncated.is_zero() || truncated.together().is_zero())
}

/// Two atoms agree: identical, or a zero expanded difference, or (dummy
/// labels) the tensor canonical form of the difference.
fn atoms_equal(a: &Atom, b: &Atom) -> bool {
    a == b || (a - b).expand().is_zero() || uncontracted::agree(a, b).unwrap_or(false)
}

// ---------------------------------------------------------------------------
// The record comparator (independent of the emitter's tests)
// ---------------------------------------------------------------------------

/// Every difference between an emitted record and a fixture, line by line
/// (the label is free text and ignored).
fn record_differences(ours: &FactorizedNumerator, theirs: &FactorizedNumerator) -> Vec<String> {
    let mut out = Vec::new();
    let mut want = |ok: bool, what: String| {
        if !ok {
            out.push(what);
        }
    };
    match (&ours.kind, &theirs.kind) {
        (FactorizedKind::Closure(a), FactorizedKind::Closure(b)) => {
            want(
                a.edge == b.edge,
                format!("closure edge {} vs {}", a.edge, b.edge),
            );
            want(
                atoms_equal(&a.mass_squared, &b.mass_squared),
                format!("closure mass {} vs {}", a.mass_squared, b.mass_squared),
            );
            want(
                atoms_equal(&a.auxiliary_mass_squared, &b.auxiliary_mass_squared),
                format!(
                    "auxiliary mass {} vs {}",
                    a.auxiliary_mass_squared, b.auxiliary_mass_squared
                ),
            );
            want(
                a.power == b.power,
                format!("closure power {} vs {}", a.power, b.power),
            );
            want(
                a.independent_externals == b.independent_externals,
                "externals".into(),
            );
            want(
                a.taylor == b.taylor,
                format!("taylor {:?} vs {:?}", a.taylor, b.taylor),
            );
            want(
                format!("{:?}", a.softs) == format!("{:?}", b.softs),
                "softs".into(),
            );
            want(
                a.original_cycles == b.original_cycles,
                format!(
                    "original_cycles {:?} vs {:?}",
                    a.original_cycles, b.original_cycles
                ),
            );
            want(
                a.depth == b.depth,
                format!("depth {:?} vs {:?}", a.depth, b.depth),
            );
            want(
                format!("{:?}", a.weights) == format!("{:?}", b.weights),
                "weights".into(),
            );
        }
        _ => want(false, "kind is not closure on both sides".into()),
    }
    want(ours.dimension == theirs.dimension, "dimension".into());
    want(
        ours.loop_edges == theirs.loop_edges,
        format!(
            "loop edges {:?} vs {:?}",
            ours.loop_edges, theirs.loop_edges
        ),
    );
    want(
        ours.external_edges == theirs.external_edges,
        format!(
            "external edges {:?} vs {:?}",
            ours.external_edges, theirs.external_edges
        ),
    );
    fn edges(
        r: &FactorizedNumerator,
    ) -> BTreeMap<usize, &vakint::p_integral::factorized::FactorizedEdge> {
        r.edges.iter().map(|e| (e.id, e)).collect()
    }
    let (a, b) = (edges(ours), edges(theirs));
    want(
        a.keys().eq(b.keys()),
        format!("edge ids {:?} vs {:?}", a.keys(), b.keys()),
    );
    for (id, x) in &a {
        let Some(y) = b.get(id) else { continue };
        want(
            x.source == y.source
                && x.target == y.target
                && x.source_hedge == y.source_hedge
                && x.target_hedge == y.target_hedge,
            format!("edge {id}: ends"),
        );
        want(
            x.particle == y.particle && x.line == y.line && x.power == y.power,
            format!("edge {id}: particle / line / power"),
        );
        want(
            x.momentum == y.momentum,
            format!("edge {id}: momentum {:?} vs {:?}", x.momentum, y.momentum),
        );
        want(
            atoms_equal(&x.mass_squared, &y.mass_squared),
            format!("edge {id}: mass"),
        );
        want(
            atoms_equal(&x.factor, &y.factor),
            format!("edge {id}: factor {} vs {}", x.factor, y.factor),
        );
    }
    fn vertices(
        r: &FactorizedNumerator,
    ) -> BTreeMap<usize, &vakint::p_integral::factorized::FactorizedVertex> {
        r.vertices.iter().map(|e| (e.id, e)).collect()
    }
    let (a, b) = (vertices(ours), vertices(theirs));
    want(
        a.keys().eq(b.keys()),
        format!("vertex ids {:?} vs {:?}", a.keys(), b.keys()),
    );
    for (id, x) in &a {
        let Some(y) = b.get(id) else { continue };
        want(
            x.rule == y.rule,
            format!("vertex {id}: rule {} vs {}", x.rule, y.rule),
        );
        want(
            atoms_equal(&x.factor, &y.factor),
            format!("vertex {id}: factor"),
        );
    }
    let canonical = |atoms: &[Atom]| {
        atoms
            .iter()
            .map(|a| a.to_canonical_string())
            .collect::<BTreeSet<_>>()
    };
    want(
        canonical(&ours.ports) == canonical(&theirs.ports),
        "ports".into(),
    );
    match (&ours.projector, &theirs.projector) {
        (Some(x), Some(y)) => {
            want(x.p2_power == y.p2_power, "projector p2_power".into());
            want(atoms_equal(&x.tensor, &y.tensor), "projector tensor".into());
        }
        (None, None) => {}
        _ => want(false, "projector present on one side only".into()),
    }
    let (x, y) = (&ours.prefactor, &theirs.prefactor);
    want(
        x.loop_sign == y.loop_sign && x.nf_power == y.nf_power,
        "prefactor loop_sign / nf_power".into(),
    );
    want(
        atoms_equal(&x.overall, &y.overall),
        "prefactor overall".into(),
    );
    want(
        atoms_equal(&x.scalars, &y.scalars),
        "prefactor scalars".into(),
    );
    want(atoms_equal(&x.colour, &y.colour), "prefactor colour".into());
    want(
        format!("{:?}", ours.colour_branch) == format!("{:?}", theirs.colour_branch),
        "colour branch".into(),
    );
    fn insertions(
        r: &FactorizedNumerator,
    ) -> BTreeMap<usize, &vakint::p_integral::factorized::Insertion> {
        r.insertions.iter().map(|e| (e.vertex, e)).collect()
    }
    let (a, b) = (insertions(ours), insertions(theirs));
    want(
        a.keys().eq(b.keys()),
        format!("insertion vertices {:?} vs {:?}", a.keys(), b.keys()),
    );
    for (v, x) in &a {
        let Some(y) = b.get(v) else { continue };
        want(x.hedges == y.hedges, format!("insertion {v}: hedges"));
        want(x.cycles == y.cycles, format!("insertion {v}: cycles"));
        want(
            atoms_equal(&(&x.coefficient * &x.factor), &(&y.coefficient * &y.factor)),
            format!("insertion {v}: coefficient x factor"),
        );
    }
    out
}

// ---------------------------------------------------------------------------
// The oracle
// ---------------------------------------------------------------------------

/// The oracle value of an emitted record text: vakint's reader, the closure
/// entry with the R* settings at the record's depth, `finish_tensor`.
fn oracle_value(vakint: &Vakint, text: &str, id: usize) -> Result<Atom> {
    let value = record(text)?;
    let FactorizedKind::Closure(spec) = &value.kind else {
        return Err(eyre!("not a closure record"));
    };
    let mut vakint_settings = settings().vakint.true_settings();
    vakint_settings.number_of_terms_in_epsilon_expansion =
        i64::try_from(spec.depth.ok_or_else(|| eyre!("no depth"))?)?;
    let prefactor = value.prefactor.total();
    let evaluated = evaluate_factorized_closure(
        vakint,
        &vakint_settings,
        std::slice::from_ref(&value),
        std::slice::from_ref(&prefactor),
    )
    .map_err(|error| eyre!("{error}"))?;
    let raw = evaluated
        .records
        .into_iter()
        .next()
        .ok_or_else(|| eyre!("no record value"))?
        .value;
    RStarIntegrand::finish_tensor(raw, id)
}

/// The physical external rules of a two-point graph: `Q(e, x___) -> ±Q(x0,
/// x___)` (momentum flowing through).
fn external_rules(graph: &Graph) -> (EdgeIndex, Vec<Replacement>) {
    let external = graph
        .iter_edges()
        .filter_map(|(pair, edge, _)| match pair {
            HedgePair::Unpaired { flow, .. } => Some((edge, flow)),
            _ => None,
        })
        .collect::<Vec<_>>();
    let (p, p_flow) = external[0];
    let rules = external
        .iter()
        .map(|(edge, flow)| {
            let sign = if *edge == p || *flow != p_flow { 1 } else { -1 };
            Replacement::new(
                function!(GS.emr_mom, usize::from(*edge), W_.x___).to_pattern(),
                (Atom::num(sign) * function!(GS.emr_mom, usize::from(p), W_.x___)).to_pattern(),
            )
        })
        .collect();
    (p, rules)
}

fn p_squared(x0: usize) -> Atom {
    function!(
        SPENSO_TAG.dot,
        function!(GS.emr_mom, x0, minkowski()),
        function!(GS.emr_mom, x0, minkowski())
    )
}

/// The closed metric root of the beta: `c = transverse / (𝑖 p² (d - 1))`
/// (graph prefactor 1 for these diagrams).
fn metric_root(graph: &Graph) -> ClosedRootNormalization {
    let (external, external_routing) = external_rules(graph);
    let factor = Atom::var(symbolica::symbol!("beta_projection::transverse"))
        * graph
            .global_prefactor
            .num
            .replace_multiple(&external_routing)
        / (Atom::i() * p_squared(usize::from(external)))
        / (Atom::var(GS.dim) - 1);
    ClosedRootNormalization {
        external,
        external_routing,
        factor,
    }
}

/// The colour branch of a fixture diagram: the colour-factored DOT with the
/// vertex and internal-edge factors the fixture's bare record states. The
/// beta weights the colour rows of a rule and its branch carries `+𝑖 g` on a
/// gluon propagator where the first colour row of the DOT rule gives `-𝑖 g`,
/// so the factors must come from the branch: the record's factor lines then
/// test the emitter's pass-through, not its source.
fn branch(dot: &str, bare: &FactorizedNumerator) -> Result<feynkit_graph::FeynmanDiagram> {
    let diagram = colour_factored(&load(dot)?)?;
    let vertices = bare
        .vertices
        .iter()
        .filter(|v| !v.rule.starts_with("ct:"))
        .map(|v| (v.id, v.factor.map_mink_dim(Atom::num(4))))
        .collect::<BTreeMap<_, _>>();
    let edges = bare
        .edges
        .iter()
        .filter(|e| !e.is_external())
        .map(|e| (e.id, e.factor.map_mink_dim(Atom::num(4))))
        .collect::<BTreeMap<_, _>>();
    let fragments = diagram
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
        .with_numerator_fragments(fragments, edge_fragments)
        .map_err(|error| eyre!("{error}"))
}

/// One diagram root through R* (the beta's closed metric root, the
/// remainder-boundary path on) under `switch`; every attempt captured.
fn root_run(
    diagram: &feynkit_graph::FeynmanDiagram,
    switch: Option<Switch>,
) -> Result<(Atom, Vec<Captured>)> {
    root_run_with(&Vakint::new()?, diagram, switch)
}

fn root_run_with(
    vakint: &Vakint,
    diagram: &feynkit_graph::FeynmanDiagram,
    switch: Option<Switch>,
) -> Result<(Atom, Vec<Captured>)> {
    let uv = settings();
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let plan = RStarPlan::with_automatic_edge(&graph, model, &uv)?;
    let closed_root = metric_root(&graph);
    let callback = |tensor: &Atom| -> Result<Atom> {
        Ok(tensor.replace_multiple(&closed_root.external_routing) * &closed_root.factor)
    };
    let run = || {
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
    };
    let (value, captured) = match switch {
        Some(switch) => with_switch(switch, run),
        None => run(),
    };
    Ok((value?, captured))
}

fn on_root() -> Switch {
    Switch::new(Mode::On, &[Kind::Root])
}

fn check_root() -> Switch {
    Switch::new(Mode::Check, &[Kind::Root])
}

/// The internal edge ids and `ct:` rules of a record.
fn key(value: &FactorizedNumerator) -> (Vec<usize>, Vec<String>) {
    let edges = value
        .edges
        .iter()
        .filter(|e| !e.is_external())
        .map(|e| e.id)
        .collect();
    let mut rules = value
        .vertices
        .iter()
        .filter(|v| v.rule.starts_with("ct:"))
        .map(|v| v.rule.clone())
        .collect::<Vec<_>>();
    rules.sort();
    (edges, rules)
}

/// Every root attempt must have been evaluated by vakint (or be a zero
/// insertion): a decline is not a pass. `allowed` lists the decline codes a
/// test expects.
fn assert_accounted(what: &str, captured: &[Captured], allowed: &[&str]) -> Result<usize> {
    let roots = captured
        .iter()
        .filter(|a| a.kind == Kind::Root)
        .collect::<Vec<_>>();
    ensure!(!roots.is_empty(), "{what}: no root record attempted");
    let unexpected = roots
        .iter()
        .filter(|a| {
            !(a.outcome == "uncontracted"
                || a.outcome == "zero"
                || allowed
                    .iter()
                    .any(|code| a.outcome == format!("declined:{code}")))
        })
        .map(|a| format!("id {}: {}", a.id, a.outcome))
        .collect::<Vec<_>>();
    ensure!(
        unexpected.is_empty(),
        "{what}: root records not evaluated (a decline is a fallback, not a pass): {unexpected:?}"
    );
    Ok(roots.len())
}

/// One fixture record against the captured attempts: the round trip (the
/// emitted text is a fixed point of vakint's writer and equals the fixture
/// up to canonical order), the oracle value and the captured value through
/// `valid_through`, the check.
fn compare_record(
    what: &str,
    vakint: &Vakint,
    captured: &[Captured],
    input: &str,
    value_text: &str,
) -> Result<FactorizedNumerator> {
    let theirs = record(input).wrap_err_with(|| what.to_owned())?;
    let wanted = key(&theirs);
    let mut found = Vec::new();
    for attempt in captured.iter().filter(|a| a.kind == Kind::Root) {
        let Some(text) = &attempt.text else { continue };
        let ours = record(text).wrap_err_with(|| what.to_owned())?;
        if key(&ours) == wanted {
            found.push((attempt, ours));
        }
    }
    ensure!(
        found.len() == 1,
        "{what}: {} captured records for {wanted:?}",
        found.len()
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
        "{what}: the emitted record differs from input.txt: {differences:#?}"
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
        .ok_or_else(|| eyre!("{what}: no captured value"))?;
    ensure!(
        equal_through(value, &expected, order)?,
        "{what}: the switch's value differs from value.txt through eps^{order}"
    );
    if let Some(contracted) = &attempt.contracted {
        ensure!(
            equal_through(contracted, &expected, order)?,
            "{what}: today's value differs from value.txt (the fixture is not this branch)"
        );
        ensure!(
            attempt.equal == Some(true),
            "{what}: check equal={:?}",
            attempt.equal
        );
    }
    eprintln!(
        "U7A_ORACLE {what} id={} insertions={} identical={:?} value_zero={}",
        attempt.id,
        ours.insertions.len(),
        attempt.identical,
        expected.is_zero()
    );
    Ok(ours)
}

/// One diagram's fixture records from the graph: off, check (every record
/// compared), on (the diagram's poles equal off's). `allowed`: the decline
/// codes the caller expects (each one a review finding).
fn diagram_from_graph(
    what: &str,
    dot: &str,
    bare: &str,
    records: &[(String, String, String)],
    allowed: &[&str],
    vakint: &Vakint,
) -> Result<Vec<FactorizedNumerator>> {
    let bare = record(bare)?;
    let diagram = branch(dot, &bare)?;
    let started = std::time::Instant::now();
    let (off, nothing) = root_run_with(vakint, &diagram, Some(Switch::OFF))?;
    ensure!(nothing.is_empty(), "{what}: off reports nothing");
    let (checked, captured) = root_run_with(vakint, &diagram, Some(check_root()))?;
    ensure!(checked == off, "{what}: check returns today's value");
    let attempted = assert_accounted(what, &captured, allowed)?;
    let mut out = Vec::new();
    for (name, input, value) in records {
        out.push(compare_record(
            &format!("{what} {name}"),
            vakint,
            &captured,
            input,
            value,
        )?);
    }
    let (on, on_captured) = root_run_with(vakint, &diagram, Some(on_root()))?;
    assert_accounted(what, &on_captured, allowed)?;
    ensure!(
        atoms_equal(&on, &off),
        "{what}: on differs from off\noff {off}\non  {on}"
    );
    eprintln!(
        "U7A_ORACLE {what}: {attempted} root records attempted, {} fixtures equal, on == off, {:.1} s",
        records.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(out)
}

// ---------------------------------------------------------------------------
// Graph-driven: the k1-quick diagrams the emitter's own tests do not use
// ---------------------------------------------------------------------------

/// k1-quick L4 40, 42 and 49 (seven boundary remainders each, 21 records):
/// emitted from the graph, record == input.txt, oracle value == value.txt,
/// check equal, on == off. L4-49 record 7 carries two poled insertions: its
/// value is nonzero and vakint's pole order is their SUM (R1). Every other
/// root record (the one-loop remainders, no fixture) must be evaluated too,
/// except L4-49's one-loop remainder, which vakint declines today (see
/// [`k1_quick_l4_49_one_loop_remainder_is_evaluated`]).
#[test]
#[ignore = "reads k1-quick from GAMMALOOP_U7A_FIXTURES_DIR (about 10 min single core)"]
fn k1_quick_l4_40_42_49_from_the_graph() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let root = fixtures_dir()?.join("k1-quick");
        for (diagram, dirs) in [
            (
                "L4-40",
                (1..=7).map(|i| format!("rr{i}")).collect::<Vec<_>>(),
            ),
            (
                "L4-42",
                (1..=7).map(|i| format!("rr{i}")).collect::<Vec<_>>(),
            ),
            ("L4-49", (1..=7).map(|i| format!("{i}")).collect::<Vec<_>>()),
        ] {
            let base = root.join(diagram);
            let number = diagram.trim_start_matches("L4-");
            let dot = read(
                &base
                    .join(&dirs[0])
                    .join(format!("raw/diagram-0000{number}.dot")),
            )?;
            let bare = read(&base.join(&dirs[0]).join("input.txt"))?;
            let records = dirs
                .iter()
                .map(|dir| {
                    Ok((
                        dir.clone(),
                        read(&base.join(dir).join("input.txt"))?,
                        read(&base.join(dir).join("value.txt"))?,
                    ))
                })
                .collect::<Result<Vec<_>>>()?;
            // L4-49: its one-loop remainder (no fixture) is declined by vakint,
            // `closure: every line is massive` (a finding, tested below); its
            // designed code (U7a review B1) is `vakint:all_massive`
            let allowed: &[&str] = if diagram == "L4-49" {
                &["vakint:all_massive"]
            } else {
                &[]
            };
            let emitted =
                diagram_from_graph(diagram, &dot, &bare, &records, allowed, &Vakint::new()?)?;
            ensure!(emitted.len() == 7, "{diagram}: seven records");
            if diagram == "L4-49" {
                let two = &emitted[6];
                ensure!(two.insertions.len() == 2, "L4-49/7: two insertions");
                let order = two
                    .insertion_pole_order(GS.dim_epsilon)
                    .map_err(|e| eyre!("{e}"))?;
                ensure!(order == 2, "L4-49/7: pole order {order}, the SUM is 2");
                ensure!(
                    !fixture_value(&records[6].2)?.is_zero(),
                    "L4-49/7: a nonzero value"
                );
            }
        }
        Ok(())
    })
}

/// k1-l5 L5-60 (15 K1 records of a five-loop root) from the graph. The 4-loop
/// p-integrals need the RustRed 1-4 loop package (`VAKINT_TEST_P_LIBRARY_1_4`,
/// loaded as the beta does); about 30 min single core, under 4 GB.
#[test]
#[ignore = "five loops: reads k1-l5 from GAMMALOOP_U7A_FIXTURES_DIR; set GAMMALOOP_U7A_L5=1"]
fn k1_l5_60_from_the_graph() -> Result<()> {
    if std::env::var("GAMMALOOP_U7A_L5").as_deref() != Ok("1") {
        return Err(eyre!("set GAMMALOOP_U7A_L5=1 to run the five-loop root"));
    }
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let base = fixtures_dir()?.join("k1-l5/L5-60");
        let mut dirs = std::fs::read_dir(&base)?
            .filter_map(|entry| entry.ok())
            .filter(|entry| entry.path().join("input.txt").is_file())
            .map(|entry| entry.file_name().to_string_lossy().into_owned())
            .collect::<Vec<_>>();
        dirs.sort_by_key(|dir| {
            dir.trim_start_matches(|c: char| !c.is_ascii_digit())
                .parse::<usize>()
                .unwrap_or(usize::MAX)
        });
        ensure!(dirs.len() == 15, "L5-60: {} records", dirs.len());
        let dot = std::fs::read_dir(base.join(&dirs[0]).join("raw"))?
            .filter_map(|entry| entry.ok())
            .map(|entry| entry.path())
            .find(|path| path.extension().is_some_and(|e| e == "dot"))
            .ok_or_else(|| eyre!("L5-60: no DOT under {}/raw", dirs[0]))?;
        let records = dirs
            .iter()
            .map(|dir| {
                Ok((
                    dir.clone(),
                    read(&base.join(dir).join("input.txt"))?,
                    read(&base.join(dir).join("value.txt"))?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        let bare = read(&base.join(&dirs[0]).join("input.txt"))?;
        let package = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4")
            .ok_or_else(|| eyre!("set VAKINT_TEST_P_LIBRARY_1_4 to the 1-4 loop p package"))?;
        let mut vakint = Vakint::new()?;
        vakint
            .load_p_library_directory(std::path::PathBuf::from(package), false)
            .map_err(|e| eyre!("{e}"))?;
        diagram_from_graph("L5-60", &read(&dot)?, &bare, &records, &[], &vakint)?;
        Ok(())
    })
}

/// Check == off and on == off on in-repo production diagrams the emitter's
/// own tests do not use (colour factored, the beta's closed metric root):
/// every root record is evaluated (none declines), check returns today's
/// value exactly, on gives the same poles, and the cycles rule R2 holds.
#[test]
fn k1_check_and_on_equal_off_on_more_graphs() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        for (name, dot) in [
            ("L2/6", include_str!("testdata/l2-d6.dot")),
            ("L2/7", include_str!("testdata/l2-d7.dot")),
            ("L3/166", include_str!("testdata/l3-d166.dot")),
            ("L3/241", include_str!("testdata/l3-d241.dot")),
        ] {
            let diagram = colour_factored(&load(dot)?)?;
            let (off, nothing) = root_run(&diagram, Some(Switch::OFF))?;
            ensure!(nothing.is_empty(), "{name}: off reports nothing");
            let (checked, captured) = root_run(&diagram, Some(check_root()))?;
            ensure!(checked == off, "{name}: check returns today's value");
            let attempted = assert_accounted(name, &captured, &[])?;
            for attempt in captured.iter().filter(|a| a.kind == Kind::Root) {
                ensure!(
                    attempt.outcome == "zero" || attempt.equal == Some(true),
                    "{name}: id {} check equal={:?}",
                    attempt.id,
                    attempt.equal
                );
                if let Some(text) = &attempt.text {
                    let value = record(text)?;
                    let FactorizedKind::Closure(spec) = &value.kind else {
                        return Err(eyre!("{name}: not a closure record"));
                    };
                    let (fermion, ghost) = spec
                        .original_cycles
                        .ok_or_else(|| eyre!("{name}: no original_cycles"))?;
                    let carried = value
                        .insertions
                        .iter()
                        .fold((0, 0), |(f, g), i| (f + i.cycles.0, g + i.cycles.1));
                    ensure!(
                        carried.0 <= fermion && carried.1 <= ghost,
                        "{name}: insertion cycles {carried:?} above original ({fermion}, {ghost})"
                    );
                }
            }
            let (on, on_captured) = root_run(&diagram, Some(on_root()))?;
            assert_accounted(name, &on_captured, &[])?;
            ensure!(
                atoms_equal(&on, &off),
                "{name}: on differs from off\noff {off}\non  {on}"
            );
            eprintln!("U7A_ORACLE {name}: {attempted} root records, check == off, on == off");
        }
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// The off mode
// ---------------------------------------------------------------------------

/// The switch off is today's behaviour: the same atom as a run that never
/// touches the switch (the environment's, unset), nothing captured, no
/// Δ-cache suffix.
#[test]
fn off_is_byte_identical_to_no_switch() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let diagram = colour_factored(&load(include_str!("testdata/l2-d8.dot"))?)?;
        let (plain, plain_captured) = root_run(&diagram, None)?;
        let (off, off_captured) = root_run(&diagram, Some(Switch::OFF))?;
        ensure!(plain == off, "off differs from a run without the switch");
        ensure!(
            plain_captured.is_empty() && off_captured.is_empty(),
            "nothing is reported"
        );
        ensure!(
            uncontracted::cache_key_suffix().is_none(),
            "no cache suffix (environment)"
        );
        ensure!(
            with_switch(Switch::OFF, uncontracted::cache_key_suffix).is_none(),
            "no cache suffix (off)"
        );
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// Site-driven: one real remainder through the emitter's entry
// ---------------------------------------------------------------------------

/// L4 40 with the beta's colour branch, its plan, and the remainder of the
/// fixture `rr2` (3 loops, one child [4,5,11]).
struct L440 {
    graph: Graph,
    diagram: feynkit_graph::FeynmanDiagram,
    plan: RStarPlan,
    term: usize,
    /// The fixture's child counterterms, by contracted support.
    deltas: Vec<(linnet::half_edge::subgraph::SuBitGraph, Atom, usize)>,
    fixture: FactorizedNumerator,
}

fn l4_40() -> Result<L440> {
    let bare = record(L4_40_RR1)?;
    let diagram = branch(L4_40_DOT, &bare)?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let plan = RStarPlan::with_automatic_edge(&graph, diagram.model(), &settings())?;
    let fixture = record(L4_40_RR2)?;
    let field = |name: &str| {
        fixture
            .label
            .split_whitespace()
            .find_map(|token| token.strip_prefix(name))
            .map(str::to_owned)
            .ok_or_else(|| eyre!("no {name} in the fixture label"))
    };
    let remaining = field("remaining=")?;
    let contracted = field("contracted=")?
        .split(';')
        .map(str::to_owned)
        .collect::<BTreeSet<_>>();
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
                    == contracted
        })
        .ok_or_else(|| eyre!("no plan term for remaining={remaining} contracted={contracted:?}"))?;
    let mut deltas = Vec::new();
    for slot in &plan.ordinary_uv[term].remainder.insertions {
        let edges = uv_child_edge_ids(&graph, &slot.support);
        let rule = format!(
            "ct:edges-{}",
            edges
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join("-")
        );
        let vertex = fixture
            .vertices
            .iter()
            .find(|v| v.rule == rule)
            .ok_or_else(|| eyre!("no {rule} vertex in the fixture"))?
            .id;
        let insertion = fixture
            .insertions
            .iter()
            .find(|i| i.vertex == vertex)
            .ok_or_else(|| eyre!("no insertion at {vertex}"))?;
        deltas.push((
            slot.support.clone(),
            &insertion.coefficient * &insertion.factor,
            vertex,
        ));
    }
    Ok(L440 {
        graph,
        diagram,
        plan,
        term,
        deltas,
        fixture,
    })
}

/// The attempt of the emitter's entry on the `rr2` remainder in on mode,
/// with the child counterterms `deltas` (by support), the root closure
/// `closure` (default: the metric root) and the plan's Taylor order
/// `taylor_order` (default: the plan's).
fn site_attempt(
    case: &L440,
    deltas: &BTreeMap<linnet::half_edge::subgraph::SuBitGraph, Atom>,
    closure: Option<Atom>,
    taylor_order: Option<i32>,
) -> Result<uncontracted::Attempt> {
    let vakint = Vakint::new()?;
    let uv = settings();
    let vakint_settings = uv.vakint.true_settings();
    let model = case.diagram.model();
    let mut plan = RStarPlan::with_automatic_edge(&case.graph, model, &uv)?;
    if let Some(order) = taylor_order {
        plan.taylor_order = order;
    }
    let term = &plan.ordinary_uv[case.term];
    let closure = match closure {
        Some(closure) => closure,
        None => RStarIntegrand::root_closure(&case.graph, &plan.support)?,
    };
    let closed_root = metric_root(&case.graph);
    let auxiliary = auxiliary_mass();
    let insertions = term
        .counterterms
        .iter()
        .map(|child| {
            deltas
                .get(&child.argument.support)
                .cloned()
                .ok_or_else(|| eyre!("no Delta for a counterterm support"))
        })
        .collect::<Result<Vec<_>>>()?;
    let site = Site {
        graph: &case.graph,
        model,
        plan: &plan,
        term,
        children: deltas,
        root_closure: Some(&closure),
        closed_root: Some(&closed_root),
        closure_data: None,
        root_projector: true,
        root_edges: &[],
        reference_mass_squared: *REFERENCE_MASS_SQUARED,
        auxiliary_mass_squared: &auxiliary,
        settings: &vakint_settings,
        vakint: &vakint,
        uv_settings: &uv,
        top_mass: Some(TopMassTerms::Keep),
        epsilon: GS.dim_epsilon,
        today: "site-test".into(),
    };
    uncontracted::attempt(&site, Kind::Root, 2, Mode::On, || {
        RStarIntegrand::new_with_closure(
            &case.graph,
            model,
            &term.remainder,
            &insertions,
            plan.marked_edge,
            Atom::var(*REFERENCE_MASS_SQUARED),
            Some(&closure),
        )
    })
}

fn deltas_of(case: &L440) -> BTreeMap<linnet::half_edge::subgraph::SuBitGraph, Atom> {
    case.deltas
        .iter()
        .map(|(support, delta, _)| (support.clone(), delta.clone()))
        .collect()
}

fn declined_code(attempt: &uncontracted::Attempt) -> Option<String> {
    match &attempt.outcome {
        Outcome::Declined(decline) => Some(decline.code.clone()),
        _ => None,
    }
}

/// The site with the fixture's own child counterterm: the record equals the
/// fixture and its value equals value.txt (no in-process child computed).
#[test]
fn site_rr2_with_the_fixture_delta_equals_its_value() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = l4_40()?;
        let attempt = site_attempt(&case, &deltas_of(&case), None, None)?;
        let text = attempt
            .text
            .clone()
            .ok_or_else(|| eyre!("no record built: {:?}", attempt.outcome))?;
        let differences = record_differences(&record(&text)?, &case.fixture);
        ensure!(differences.is_empty(), "rr2: {differences:#?}");
        let value = attempt
            .value()
            .ok_or_else(|| eyre!("rr2 declined: {:?}", attempt.outcome))?;
        let expected = fixture_value(L4_40_RR2_VALUE)?;
        let order = valid_through(L4_40_RR2_VALUE)?;
        ensure!(equal_through(&value, &expected, order)?, "rr2: value");
        let oracle = oracle_value(&Vakint::new()?, &text, 2)?;
        ensure!(
            equal_through(&oracle, &expected, order)?,
            "rr2: oracle value"
        );
        Ok(())
    })
}

/// A mass symbol inside a child counterterm (the auxiliary mass, then the
/// reference mass R) that momentum conservation cannot remove: declined
/// `insertion_mass_symbol`, no value (the hook then takes today's path).
#[test]
fn site_mass_symbol_in_an_insertion_declines() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = l4_40()?;
        for mass in [auxiliary_mass(), Atom::var(*REFERENCE_MASS_SQUARED)] {
            let mut deltas = deltas_of(&case);
            for delta in deltas.values_mut() {
                *delta = &*delta * (Atom::num(1) + &mass);
            }
            let attempt = site_attempt(&case, &deltas, None, None)?;
            ensure!(
                declined_code(&attempt).as_deref() == Some("insertion_mass_symbol"),
                "mass {mass}: outcome {:?}",
                attempt.outcome
            );
            ensure!(attempt.value().is_none(), "a declined record has no value");
            ensure!(
                attempt.text.is_none(),
                "no record text for a declined record"
            );
        }
        Ok(())
    })
}

/// `log(M_aux²) x (momentum-conservation zero)` added to the child
/// counterterm: the conservation step removes it, the record is evaluated,
/// and the value is the fixture's (the rewrite is exact on a real record).
#[test]
fn site_conservation_zero_is_removed_exactly() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = l4_40()?;
        let (_, _, vertex) = &case.deltas[0];
        // the edges at the ct: node and their signs (incoming +1)
        let mut zero = Atom::Zero;
        for edge in &case.fixture.edges {
            let mut sign = 0i64;
            if edge.target == vakint::p_integral::factorized::Endpoint::Vertex(*vertex) {
                sign += 1;
            }
            if edge.source == vakint::p_integral::factorized::Endpoint::Vertex(*vertex) {
                sign -= 1;
            }
            if sign != 0 {
                zero += Atom::num(sign) * function!(GS.emr_mom, edge.id, minkowski());
            }
        }
        let probe = case
            .fixture
            .external_edges
            .first()
            .copied()
            .flatten()
            .ok_or_else(|| eyre!("no external edge"))?;
        let zero = function!(
            SPENSO_TAG.dot,
            zero,
            function!(GS.emr_mom, probe, minkowski())
        );
        let log = parse("log(gammalooprs::{scalar}::rstar_mass_squared)")?;
        let mut deltas = deltas_of(&case);
        for delta in deltas.values_mut() {
            *delta = &*delta + &*delta * &log * &zero;
        }
        let attempt = site_attempt(&case, &deltas, None, None)?;
        let value = attempt
            .value()
            .ok_or_else(|| eyre!("declined: {:?}", attempt.outcome))?;
        let text = attempt.text.clone().unwrap_or_default();
        ensure!(
            !text
                .lines()
                .filter(|l| l.starts_with("insertion"))
                .any(|l| l.contains("rstar_mass_squared")),
            "the insertion still carries the auxiliary mass"
        );
        let expected = fixture_value(L4_40_RR2_VALUE)?;
        let order = valid_through(L4_40_RR2_VALUE)?;
        ensure!(
            equal_through(&value, &expected, order)?,
            "conserved record: value differs from value.txt"
        );
        Ok(())
    })
}

/// A Taylor degree above 2 (the format's limit): declined by vakint, the
/// designed refusal `vakint:taylor_degree` (Closure "Taylor degree", review
/// B1), no value.
#[test]
fn site_taylor_degree_three_declines() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = l4_40()?;
        let attempt = site_attempt(&case, &deltas_of(&case), None, Some(3))?;
        ensure!(
            declined_code(&attempt).as_deref() == Some("vakint:taylor_degree"),
            "degree 3: outcome {:?}",
            attempt.outcome
        );
        let Outcome::Declined(decline) = &attempt.outcome else {
            unreachable!()
        };
        ensure!(
            decline.detail.contains("Taylor degree 3"),
            "degree 3: detail {}",
            decline.detail
        );
        ensure!(attempt.value().is_none(), "a declined record has no value");
        Ok(())
    })
}

/// An open spinor chain in the root closure (two open bispinor slots, the
/// metric kept): declined by vakint, no value.
#[test]
fn site_open_spinor_chain_declines() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = l4_40()?;
        let metric = RStarIntegrand::root_closure(&case.graph, &case.plan.support)?;
        let mink = Minkowski {}.new_rep(GS.dim);
        let spin = Bispinor {}.new_rep(4);
        let open = |h: u16| spin.to_symbolic([Aind::Hedge(h, 1).to_atom()]);
        let inner = spin.to_symbolic([Aind::new_dummy().to_atom()]);
        let mu = mink.to_symbolic([Aind::new_dummy().to_atom()]);
        let chain = function!(idenso::dirac::AGS.gamma, open(9001), &inner, &mu)
            * function!(idenso::dirac::AGS.gamma, &inner, open(9002), &mu);
        let attempt = site_attempt(&case, &deltas_of(&case), Some(metric * chain), None)?;
        let code = declined_code(&attempt);
        eprintln!(
            "U7A_ORACLE open spinor chain: outcome {:?}",
            attempt.outcome
        );
        ensure!(
            code.as_deref().is_some_and(|c| c.starts_with("vakint:")),
            "open spinor chain: outcome {:?}",
            attempt.outcome
        );
        ensure!(attempt.value().is_none(), "a declined record has no value");
        Ok(())
    })
}

/// A child counterterm that is literally zero: the record is not built and
/// the value is 0.
#[test]
fn site_zero_insertion_is_zero() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = l4_40()?;
        let mut deltas = deltas_of(&case);
        for delta in deltas.values_mut() {
            *delta = Atom::Zero;
        }
        let attempt = site_attempt(&case, &deltas, None, None)?;
        ensure!(
            matches!(attempt.outcome, Outcome::Zero),
            "outcome {:?}",
            attempt.outcome
        );
        ensure!(attempt.value() == Some(Atom::Zero), "value 0");
        Ok(())
    })
}

/// FINDING (review B1, fixed): a malformed child counterterm (an open slot on
/// a half-edge that is not at its vertex) is an emitter bug, not a designed
/// refusal. It must not end as a silent per-record fallback (`declined`): it
/// is an ERROR outcome (`error:vakint:<Variant>`), which fails the diagram in
/// check mode and in a strict on run.
#[test]
fn site_malformed_insertion_is_not_a_silent_fallback() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let case = l4_40()?;
        // a metric on two half-edges that are not at the ct: vertex
        let stray = parse(
            "spenso::{symmetric,linear,real,spenso::tensor}::g(spenso::{spenso::representation,spenso::self_dual}::mink(gammalooprs::{}::dim,gammalooprs::{spenso::index}::hedge(9100,1)),spenso::{spenso::representation,spenso::self_dual}::mink(gammalooprs::{}::dim,gammalooprs::{spenso::index}::hedge(9101,1)))",
        )?;
        let mut deltas = deltas_of(&case);
        for delta in deltas.values_mut() {
            *delta = &*delta * &stray;
        }
        let result = site_attempt(&case, &deltas, None, None);
        match result {
            Err(_) => Ok(()),
            Ok(attempt) if attempt.is_error() => {
                eprintln!(
                    "U7A_ORACLE malformed insertion: outcome {:?}",
                    attempt.outcome
                );
                ensure!(attempt.value().is_none(), "an error has no value");
                Ok(())
            }
            Ok(attempt) => Err(eyre!(
                "a malformed record ends as {:?}: a fallback that hides the bug",
                attempt.outcome
            )),
        }
    })
}

/// FINDING (vakint side, review B3): the one-loop root remainder of L4
/// diagram 49 (no fixture: today's path takes it whole, `one_loop_tree_cut`)
/// is declined by vakint's closure entry, `closure: every line is massive:
/// not a closure` (a massive tadpole: only the closing line is left). Until
/// the entry serves it (or the format lists it), the emitter treats it as the
/// designed refusal `vakint:all_massive` (on the B1 list, counted in
/// UNCONTRACTED_SUMMARY); this test states the goal: the record evaluated.
#[test]
#[ignore = "vakint-side finding B3: the closure entry declines L4-49's one-loop massive tadpole remainder (vakint:all_massive); reads GAMMALOOP_U7A_FIXTURES_DIR"]
fn k1_quick_l4_49_one_loop_remainder_is_evaluated() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        assume_defaults()?;
        let base = fixtures_dir()?.join("k1-quick/L4-49/1");
        let bare = record(&read(&base.join("input.txt"))?)?;
        let diagram = branch(&read(&base.join("raw/diagram-000049.dot"))?, &bare)?;
        let (_, captured) = root_run(&diagram, Some(check_root()))?;
        assert_accounted("L4-49", &captured, &[])?;
        Ok(())
    })
}
