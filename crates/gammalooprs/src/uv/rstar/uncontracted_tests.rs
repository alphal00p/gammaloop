//! U7a (`GAMMALOOP_UNCONTRACTED_NUMERATOR`, [`super::uncontracted`]): the
//! switch, the Δ-cache suffix, the record pieces (insertion split, the Δ
//! momentum-conservation step), the comparator, and the K1 records of real
//! beta diagrams against the VERIFIED fixtures of the beta project
//! (results/consolidation/wh-fixtures): the emitted record equals the
//! fixture's `input.txt` up to canonical order (atoms by expanded
//! difference), its value equals `value.txt` through `valid_through`, and
//! check mode (both paths, pole parts compared) passes and returns today's
//! value.
//!
//! The fixture diagrams' colour branch is the fixture's: every vertex and
//! propagator factor that a fixture record states replaces the
//! colour-factored rule (`recursive_tests::colour_factored` keeps the first
//! colour row, the beta weights rows), so the in-process remainders are the
//! beta's. CAVEAT (U7a review I1): those injected factors are then read back
//! by the emitter, so "emitted record == input.txt" is circular for the edge
//! and vertex FACTORS (it tests momenta, half-edges, `ct:` nodes,
//! insertions, projector, prefactor, Taylor step, cycles and depth); the
//! values against `value.txt` stay independent. The factor extraction is
//! tested by the check-mode PRODUCT identity (the record's factors multiply
//! back to the site's numerator) on diagrams WITHOUT injection:
//! `k1_check_equals_off` (L2 0/8/12, L3 185) and `k1_check_equals_off_l4_*`
//! (L4 40/42/49/100), colour-factored as they are.
//!
//! Every run states its expectations by outcome (U7a review B2): the designed
//! declines it expects, by code and count; every other attempt evaluated by
//! vakint and, in check mode, compared equal; no error outcome (the runs are
//! strict, `GAMMALOOP_UNCONTRACTED_STRICT` on); mode on agrees with off on
//! the diagram's poles.
//!
//! Data: `testdata/u7a/` (k1-root, k1-root-rc and its r_c = 0 control, the
//! K2 fixtures and site values, the insertion-split records, the L3 016 Δ,
//! the DOTs). k1-quick and k1-l5 are read in place from
//! `GAMMALOOP_U7A_FIXTURES_DIR` (the beta's results/consolidation/wh-fixtures)
//! by `#[ignore]` tests, which FAIL when the variable is unset or its
//! directory is missing; k1-l5 also needs the beta's p library
//! (`GAMMALOOP_U7A_P_LIBRARY`, else `BETA_P_LIBRARY`).

use std::collections::{BTreeMap, BTreeSet};

use eyre::{Result, WrapErr, ensure, eyre};
use idenso::tensor::SymbolicTensor;
use linnet::half_edge::involution::{EdgeIndex, HedgePair};
use spenso::{
    network::tags::SPENSO_TAG,
    structure::representation::{Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    id::Replacement,
    symbol,
};
use vakint::{
    Vakint,
    p_integral::factorized::{
        Endpoint, FactorizedEdge, FactorizedKind, FactorizedNumerator, LineKind,
    },
};

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    graph::Graph,
    initialisation::test_initialise,
    numerator::symbolica_ext::NumeratorAtomExt,
    utils::{GS, W_},
    uv::approx::streamed::with_production_piece_limit,
};

use super::{
    RStarPlan, UvChildClosure, UvChildrenRoute,
    boundary::{BoundaryMode, ClosedRootNormalization, with_boundary_mode},
    recursive_tests::{colour_factored, load, settings},
    uncontracted::{
        self, Captured, DESIGNED_LOCAL, Decline, Kind, Mode, Switch, agree, conserve,
        is_unknown_master_text, split_insertion, split_root_factor, with_capture, with_strict,
        with_switch,
    },
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn parse(text: &str) -> Atom {
    symbolica::try_parse!(text).unwrap_or_else(|error| panic!("{error}: {text:.200}"))
}

fn minkowski() -> Atom {
    Minkowski {}.new_rep(GS.dim).to_symbolic([])
}

fn epsilon() -> symbolica::atom::Symbol {
    GS.dim_epsilon
}

/// The `value` line of a fixture value file.
fn fixture_value(text: &str) -> Atom {
    parse(
        text.lines()
            .find_map(|line| line.strip_prefix("value\t"))
            .expect("a value line"),
    )
}

/// The order `n` of a fixture's `valid_through eps^n` line (`-1` by default).
fn valid_through(text: &str) -> i64 {
    text.lines()
        .find_map(|line| line.strip_prefix("valid_through\t"))
        .and_then(|line| line.trim_start_matches("eps^").split_whitespace().next())
        .and_then(|order| order.parse().ok())
        .unwrap_or(-1)
}

/// `got == expected` through ε^order (the fixtures' identity rule).
fn equal_through(got: &Atom, expected: &Atom, order: i64) -> Result<bool> {
    let difference = (got - expected).expand();
    if difference.is_zero() {
        return Ok(true);
    }
    let truncated = difference.series(epsilon(), 0, order)?.to_atom().expand();
    Ok(truncated.is_zero() || truncated.together().is_zero())
}

/// The physical external rules of a two-point graph, as the beta builds
/// them: `(x0, Q(e, x___) -> ±Q(x0, x___))`.
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

/// The beta's closed metric-root normalisation (beta src/integrate.rs
/// `closed_root_scalar_normalization`): `c = transverse σ(prefactor) / (𝑖 p²
/// (d - 1))`.
fn beta_normalization(graph: &Graph) -> ClosedRootNormalization {
    let (external, external_routing) = external_rules(graph);
    let x0 = usize::from(external);
    let p2 = function!(
        SPENSO_TAG.dot,
        function!(GS.emr_mom, x0, minkowski()),
        function!(GS.emr_mom, x0, minkowski())
    );
    let factor = Atom::var(symbol!("beta_projection::transverse"))
        * graph
            .global_prefactor
            .num
            .replace_multiple(&external_routing)
        / (Atom::i() * p2)
        / (Atom::var(GS.dim) - 1);
    ClosedRootNormalization {
        external,
        external_routing,
        factor,
    }
}

/// The colour branch of a fixture: the colour-factored diagram with every
/// vertex and internal-edge factor that `records` state (mink(dim) back to
/// the DOT's mink(4)).
fn fixture_branch(
    dot: &str,
    records: &[&FactorizedNumerator],
) -> Result<feynkit_graph::FeynmanDiagram> {
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for record in records {
        for vertex in &record.vertices {
            if !vertex.rule.starts_with("ct:") {
                vertices.insert(vertex.id, vertex.factor.map_mink_dim(Atom::num(4)));
            }
        }
        for edge in &record.edges {
            if !edge.is_external() {
                edges.insert(edge.id, edge.factor.map_mink_dim(Atom::num(4)));
            }
        }
    }
    branch_with(dot, vertices, edges)
}

/// The vertex and internal-edge factors of a `wh-fixture.v0` site dump
/// (`vertex <id> rule=.. factor=<atom>`, `edge <id> .. factor=<atom>`;
/// `contracted` vertices skipped), mink(dim) back to the DOT's mink(4).
fn site_factors(text: &str) -> (BTreeMap<usize, Atom>, BTreeMap<usize, Atom>) {
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for line in text.lines() {
        let mut fields = line.split('\t');
        let (Some(key), Some(id)) = (fields.next(), fields.next()) else {
            continue;
        };
        let Some(factor) = fields.find_map(|field| field.strip_prefix("factor=")) else {
            continue;
        };
        if factor == "contracted" {
            continue;
        }
        let Ok(id) = id.parse::<usize>() else {
            continue;
        };
        let factor = parse(factor).map_mink_dim(Atom::num(4));
        match key {
            "vertex" => {
                vertices.insert(id, factor);
            }
            "edge" => {
                edges.insert(id, factor);
            }
            _ => {}
        }
    }
    (vertices, edges)
}

/// The colour-factored diagram with the given vertex and internal-edge
/// factors (the fixture's colour branch).
fn branch_with(
    dot: &str,
    vertices: BTreeMap<usize, Atom>,
    edges: BTreeMap<usize, Atom>,
) -> Result<feynkit_graph::FeynmanDiagram> {
    let diagram = colour_factored(&load(dot)?)?;
    let vertex_fragments = diagram
        .vertices()
        .filter_map(|(id, _)| vertices.get(&id.0).map(|factor| (id, factor.clone())))
        .map(|(id, factor)| Ok((id, SymbolicTensor::infer(factor)?)))
        .collect::<Result<Vec<_>>>()?;
    let edge_fragments = diagram
        .edges()
        .filter_map(|(id, _, _)| edges.get(&id.0).map(|factor| (id, factor.clone())))
        .map(|(id, factor)| Ok((id, SymbolicTensor::infer(factor)?)))
        .collect::<Result<Vec<_>>>()?;
    diagram
        .with_numerator_fragments(vertex_fragments, edge_fragments)
        .map_err(|error| eyre!("{error}"))
}

/// One diagram root (`projected_counterterm_closing_root_scalar`, the beta's
/// full-graph closing call) under `switch`, with the remainder-boundary path
/// on (production), and every attempt reported inside it.
fn root_run(
    diagram: &feynkit_graph::FeynmanDiagram,
    switch: Switch,
) -> Result<(Atom, Vec<Captured>, usize)> {
    root_run_with(diagram, switch, &Vakint::new()?)
}

/// [`root_run`] with the given vakint engine (e.g. one with a p library).
fn root_run_with(
    diagram: &feynkit_graph::FeynmanDiagram,
    switch: Switch,
    vakint: &Vakint,
) -> Result<(Atom, Vec<Captured>, usize)> {
    let settings = settings();
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
    let closed_root = beta_normalization(&graph);
    let callback = |tensor: &Atom| -> Result<Atom> {
        Ok(tensor.replace_multiple(&closed_root.external_routing) * &closed_root.factor)
    };
    let (value, captured) = with_switch(switch, || {
        with_strict(true, || {
            with_boundary_mode(BoundaryMode::On, false, || {
                with_capture(|| {
                    plan.projected_counterterm_closing_root_scalar(
                        &graph,
                        model,
                        vakint,
                        &settings,
                        &callback,
                        &closed_root,
                    )
                })
            })
        })
    });
    Ok((value?, captured, usize::from(plan.marked_edge)))
}

/// What a run's attempts must be (U7a review B2): exactly the designed
/// `declines` (code, count; codes as in `outcome=declined:<code>`), every
/// other attempt evaluated by vakint (or a zero insertion) and, `checked`,
/// compared equal; at least one record evaluated. An error outcome never
/// reaches here (a strict or check run fails on it).
fn expect_outcomes(
    what: &str,
    captured: &[Captured],
    declines: &[(&str, usize)],
    checked: bool,
) -> Result<()> {
    let mut census = BTreeMap::<String, usize>::new();
    for attempt in captured {
        eprintln!(
            "U7A_ATTEMPT {what} graph={} kind={} id={} outcome={} identical={:?} equal={:?} conserved={}",
            attempt.graph,
            attempt.kind.as_str(),
            attempt.id,
            attempt.outcome,
            attempt.identical,
            attempt.equal,
            attempt.conserved
        );
        match attempt.outcome.as_str() {
            "uncontracted" | "zero" => {
                if checked {
                    ensure!(
                        attempt.equal == Some(true),
                        "{what}: {} id {} was not compared equal ({:?})",
                        attempt.kind.as_str(),
                        attempt.id,
                        attempt.equal
                    );
                }
            }
            outcome => {
                let code = outcome.strip_prefix("declined:").ok_or_else(|| {
                    eyre!(
                        "{what}: {} id {} outcome {outcome}",
                        attempt.kind.as_str(),
                        attempt.id
                    )
                })?;
                *census.entry(code.to_owned()).or_default() += 1;
            }
        }
    }
    let expected = declines
        .iter()
        .filter(|(_, n)| *n > 0)
        .map(|(code, n)| ((*code).to_owned(), *n))
        .collect::<BTreeMap<_, _>>();
    ensure!(
        census == expected,
        "{what}: the declines {census:?} are not the expected {expected:?}"
    );
    ensure!(
        captured.iter().any(|a| a.outcome == "uncontracted"),
        "{what}: no record was evaluated by vakint"
    );
    Ok(())
}

/// One diagram under the switch with `kinds` enabled: check (the attempts
/// as [`expect_outcomes`]) and, `full`, also off (nothing reported; check
/// returns its value) and mode on (the same outcomes; the diagram's poles
/// equal off's). The attempts of the check run.
fn checked_run(
    what: &str,
    diagram: &feynkit_graph::FeynmanDiagram,
    vakint: &Vakint,
    kinds: &[Kind],
    declines: &[(&str, usize)],
    full: bool,
) -> Result<Vec<Captured>> {
    let (checked, captured, _) = root_run_with(diagram, Switch::new(Mode::Check, kinds), vakint)?;
    expect_outcomes(&format!("{what} check"), &captured, declines, true)?;
    if full {
        let (off, nothing, _) = root_run_with(diagram, Switch::OFF, vakint)?;
        ensure!(nothing.is_empty(), "{what}: mode off reports nothing");
        ensure!(
            checked == off,
            "{what}: check mode returns today's value\noff   {off}\ncheck {checked}"
        );
        let (value, attempts, _) = root_run_with(diagram, Switch::new(Mode::On, kinds), vakint)?;
        expect_outcomes(&format!("{what} on"), &attempts, declines, false)?;
        ensure!(
            atoms_agree(&value, &off),
            "{what}: mode on differs from off\noff {off}\non  {value}"
        );
    }
    Ok(captured)
}

/// The internal edge ids and the `ct:` rules of a record (its remainder).
fn record_key(record: &FactorizedNumerator) -> (Vec<usize>, Vec<String>) {
    let edges = record
        .edges
        .iter()
        .filter(|edge| !edge.is_external())
        .map(|edge| edge.id)
        .collect();
    let mut rules = record
        .vertices
        .iter()
        .filter(|vertex| vertex.rule.starts_with("ct:"))
        .map(|vertex| vertex.rule.clone())
        .collect::<Vec<_>>();
    rules.sort();
    (edges, rules)
}

fn atoms_agree(a: &Atom, b: &Atom) -> bool {
    a == b || (a - b).expand().is_zero() || agree(a, b).unwrap_or(false)
}

/// The comparator of U7a-DESIGN 11.2: the emitted record against a fixture
/// "up to canonical order" (the label ignored; structural lines equal; atoms
/// by expanded difference, then the tensor canonical form; an insertion by
/// its product coefficient × factor). The differences found, empty when they
/// agree.
fn record_differences(ours: &FactorizedNumerator, theirs: &FactorizedNumerator) -> Vec<String> {
    let mut out = Vec::new();
    let mut check = |ok: bool, what: String| {
        if !ok {
            out.push(what);
        }
    };
    match (&ours.kind, &theirs.kind) {
        (FactorizedKind::Closure(a), FactorizedKind::Closure(b)) => {
            check(
                a.edge == b.edge,
                format!("closure edge {} vs {}", a.edge, b.edge),
            );
            check(a.mass_squared == b.mass_squared, "closure mass".into());
            check(
                a.auxiliary_mass_squared == b.auxiliary_mass_squared,
                "auxiliary mass".into(),
            );
            check(a.power == b.power, "closure power".into());
            check(
                a.independent_externals == b.independent_externals,
                "externals".into(),
            );
            check(
                a.taylor == b.taylor,
                format!("taylor {:?} vs {:?}", a.taylor, b.taylor),
            );
            check(a.softs == b.softs, "softs".into());
            check(
                a.original_cycles == b.original_cycles,
                format!(
                    "original_cycles {:?} vs {:?}",
                    a.original_cycles, b.original_cycles
                ),
            );
            check(
                a.depth == b.depth,
                format!("depth {:?} vs {:?}", a.depth, b.depth),
            );
        }
        _ => check(false, "kind".into()),
    }
    check(
        ours.loop_edges == theirs.loop_edges,
        format!(
            "loop edges {:?} vs {:?}",
            ours.loop_edges, theirs.loop_edges
        ),
    );
    check(
        ours.external_edges == theirs.external_edges,
        format!(
            "external edges {:?} vs {:?}",
            ours.external_edges, theirs.external_edges
        ),
    );
    fn edges(r: &FactorizedNumerator) -> BTreeMap<usize, &FactorizedEdge> {
        r.edges.iter().map(|e| (e.id, e)).collect()
    }
    let (a, b) = (edges(ours), edges(theirs));
    check(
        a.keys().eq(b.keys()),
        format!(
            "edge ids {:?} vs {:?}",
            a.keys().collect::<Vec<_>>(),
            b.keys().collect::<Vec<_>>()
        ),
    );
    for (id, x) in &a {
        let Some(y) = b.get(id) else { continue };
        check(
            (x.source, x.target, x.source_hedge, x.target_hedge)
                == (y.source, y.target, y.source_hedge, y.target_hedge),
            format!("edge {id} ends"),
        );
        check(
            x.particle == y.particle,
            format!("edge {id} particle {} vs {}", x.particle, y.particle),
        );
        check(x.line == y.line, format!("edge {id} line"));
        check(x.power == y.power, format!("edge {id} power"));
        check(
            x.momentum == y.momentum,
            format!("edge {id} momentum {:?} vs {:?}", x.momentum, y.momentum),
        );
        check(x.mass_squared == y.mass_squared, format!("edge {id} mass"));
        check(
            atoms_agree(&x.factor, &y.factor),
            format!("edge {id} factor"),
        );
    }
    fn vertices(
        r: &FactorizedNumerator,
    ) -> BTreeMap<usize, &vakint::p_integral::factorized::FactorizedVertex> {
        r.vertices.iter().map(|v| (v.id, v)).collect()
    }
    let (a, b) = (vertices(ours), vertices(theirs));
    check(
        a.keys().eq(b.keys()),
        format!(
            "vertex ids {:?} vs {:?}",
            a.keys().collect::<Vec<_>>(),
            b.keys().collect::<Vec<_>>()
        ),
    );
    for (id, x) in &a {
        let Some(y) = b.get(id) else { continue };
        check(
            x.rule == y.rule,
            format!("vertex {id} rule {} vs {}", x.rule, y.rule),
        );
        check(
            atoms_agree(&x.factor, &y.factor),
            format!("vertex {id} factor"),
        );
    }
    let ports = |r: &FactorizedNumerator| {
        r.ports
            .iter()
            .map(|p| p.to_canonical_string())
            .collect::<BTreeSet<_>>()
    };
    check(ports(ours) == ports(theirs), "ports".into());
    match (&ours.projector, &theirs.projector) {
        (None, None) => {}
        (Some(x), Some(y)) => {
            check(
                x.p2_power == y.p2_power,
                format!("p2_power {} vs {}", x.p2_power, y.p2_power),
            );
            check(atoms_agree(&x.tensor, &y.tensor), "projector tensor".into());
        }
        _ => check(false, "projector presence".into()),
    }
    let (p, q) = (&ours.prefactor, &theirs.prefactor);
    check(
        (p.loop_sign, p.nf_power) == (q.loop_sign, q.nf_power),
        "prefactor integers".into(),
    );
    check(
        atoms_agree(&p.overall, &q.overall),
        "prefactor_overall".into(),
    );
    check(
        atoms_agree(&p.scalars, &q.scalars),
        format!("prefactor_scalars {} vs {}", p.scalars, q.scalars),
    );
    check(atoms_agree(&p.colour, &q.colour), "prefactor_colour".into());
    fn insertions(
        r: &FactorizedNumerator,
    ) -> BTreeMap<usize, &vakint::p_integral::factorized::Insertion> {
        r.insertions.iter().map(|i| (i.vertex, i)).collect()
    }
    let (a, b) = (insertions(ours), insertions(theirs));
    check(
        a.keys().eq(b.keys()),
        format!(
            "insertion vertices {:?} vs {:?}",
            a.keys().collect::<Vec<_>>(),
            b.keys().collect::<Vec<_>>()
        ),
    );
    for (vertex, x) in &a {
        let Some(y) = b.get(vertex) else { continue };
        check(
            x.hedges == y.hedges,
            format!("insertion {vertex} hedges {:?} vs {:?}", x.hedges, y.hedges),
        );
        check(x.cycles == y.cycles, format!("insertion {vertex} cycles"));
        check(
            atoms_agree(&(&x.coefficient * &x.factor), &(&y.coefficient * &y.factor)),
            format!("insertion {vertex} coefficient x factor"),
        );
    }
    out
}

/// One fixture record against the attempt captured for the same remainder:
/// the record (comparator), the value through `valid_through`, the check.
fn compare_fixture(what: &str, captured: &[Captured], input: &str, value_text: &str) -> Result<()> {
    compare_fixture_of(what, Kind::Root, captured, input, value_text)
}

/// [`compare_fixture`] for a record of `kind`.
fn compare_fixture_of(
    what: &str,
    kind: Kind,
    captured: &[Captured],
    input: &str,
    value_text: &str,
) -> Result<()> {
    let theirs = FactorizedNumerator::from_text(input).map_err(|e| eyre!("{what}: {e}"))?;
    let key = record_key(&theirs);
    let mut matches = Vec::new();
    for attempt in captured.iter().filter(|a| a.kind == kind) {
        let Some(text) = &attempt.text else { continue };
        let ours = FactorizedNumerator::from_text(text).map_err(|e| eyre!("{what}: {e}"))?;
        if record_key(&ours) == key {
            matches.push((attempt, ours));
        }
    }
    ensure!(
        matches.len() == 1,
        "{what}: {} captured records for the remainder {key:?} (outcomes {:?})",
        matches.len(),
        captured
            .iter()
            .map(|a| format!("{}:{}", a.id, a.outcome))
            .collect::<Vec<_>>()
    );
    let (attempt, ours) = &matches[0];
    // vakint's frozen writer and reader: a byte-identical round trip
    ensure!(
        ours.to_text().map_err(|e| eyre!("{e}"))? == *attempt.text.as_ref().unwrap(),
        "{what}: the record's text is not a fixed point of vakint's writer"
    );
    let differences = record_differences(ours, &theirs);
    ensure!(
        differences.is_empty(),
        "{what}: the emitted record differs from the fixture: {differences:#?}\n{}",
        attempt.text.as_ref().unwrap()
    );
    ensure!(
        attempt.outcome == "uncontracted",
        "{what}: outcome {}",
        attempt.outcome
    );
    let value = attempt
        .value
        .as_ref()
        .ok_or_else(|| eyre!("{what}: no value"))?;
    let expected = fixture_value(value_text);
    let order = valid_through(value_text);
    if let Some(contracted) = &attempt.contracted {
        ensure!(
            equal_through(contracted, &expected, order)?,
            "{what}: today's value differs from value.txt through eps^{order} (the fixture is not this branch)\ntoday  {contracted}\ntheirs {expected}"
        );
    }
    ensure!(
        equal_through(value, &expected, order)?,
        "{what}: the value differs from value.txt through eps^{order}\nours  {value}\ntheirs {expected}"
    );
    if let Some(equal) = attempt.equal {
        ensure!(equal, "{what}: check mode equal=false");
        eprintln!(
            "U7A_FIXTURE {what} id={} identical={} equal={equal}",
            attempt.id,
            attempt.identical.unwrap_or(false)
        );
    }
    Ok(())
}

/// The fixtures directory, from `GAMMALOOP_U7A_FIXTURES_DIR` only (U7a
/// review I2): an error when it is unset or not a directory, so an
/// `#[ignore]` test run on request never passes with nothing tested.
fn fixtures_dir() -> Result<std::path::PathBuf> {
    let dir = std::env::var("GAMMALOOP_U7A_FIXTURES_DIR").map_err(|_| {
        eyre!(
            "GAMMALOOP_U7A_FIXTURES_DIR is not set (the beta's results/consolidation/wh-fixtures)"
        )
    })?;
    let dir = std::path::PathBuf::from(dir);
    ensure!(
        dir.is_dir(),
        "GAMMALOOP_U7A_FIXTURES_DIR={} is not a directory",
        dir.display()
    );
    Ok(dir)
}

const ROOT_CHILD: [Kind; 2] = [Kind::Root, Kind::Child];

// ---------------------------------------------------------------------------
// The switch and the cache key
// ---------------------------------------------------------------------------

#[test]
fn switch_parses_every_form() {
    for off in [None, Some(""), Some("0"), Some("off"), Some(" off ")] {
        assert_eq!(Switch::parse(off), Ok(Switch::OFF), "{off:?}");
    }
    let on = Switch::parse(Some("on")).unwrap();
    assert_eq!(on.mode, Mode::On);
    assert!(on.enabled(Kind::Root) && on.enabled(Kind::Child));
    assert!(!on.enabled(Kind::IrHard) && !on.enabled(Kind::Quotient));
    assert_eq!(Switch::parse(Some("1")), Ok(on));
    let check = Switch::parse(Some("check")).unwrap();
    assert_eq!(check.mode, Mode::Check);
    assert_eq!(check.enabled_kinds(), vec![Kind::Root, Kind::Child]);
    assert_eq!(
        Switch::parse(Some("check:a,b")),
        Ok(check),
        "letter aliases"
    );
    let root = Switch::parse(Some("on:root, root")).unwrap();
    assert_eq!(root.enabled_kinds(), vec![Kind::Root], "duplicates ignored");
    let inert = Switch::parse(Some("on:ir-hard,quotient,surrogate,insertion")).unwrap();
    assert!(inert.enabled_kinds().is_empty());
    assert_eq!(
        inert.inert_kinds(),
        vec![
            Kind::IrHard,
            Kind::Quotient,
            Kind::Surrogate,
            Kind::Insertion
        ]
    );
    assert_eq!(
        inert.cache_key_suffix(),
        None,
        "nothing enabled, today's key"
    );
    let mixed = Switch::parse(Some("check:c,root")).unwrap();
    assert_eq!(mixed.enabled_kinds(), vec![Kind::Root]);
    assert_eq!(
        mixed.line(),
        "UNCONTRACTED_SWITCH mode=check kinds=root inert=ir-hard"
    );
    assert_eq!(
        on.line(),
        "UNCONTRACTED_SWITCH mode=on kinds=root,child inert=none"
    );
    for bad in [
        "yes",
        "on:",
        "on:bogus",
        "on:root,",
        "off:root",
        "check:root,,child",
        "checks",
    ] {
        let error = Switch::parse(Some(bad)).expect_err(bad);
        assert!(error.contains(uncontracted::VARIABLE), "{error}");
    }
}

#[test]
fn cache_suffix_only_when_enabled() {
    let key = || {
        RStarPlan::rstar_child_key(
            "settings".into(),
            UvChildrenRoute::Rstar { min_loops: 1 },
            UvChildClosure::Off,
            Some("#uv-child-truncate=graded@1".into()),
        )
    };
    let off = with_switch(Switch::OFF, key);
    assert!(!off.contains("#uncontracted"), "{off}");
    assert!(off.ends_with("#uv-child-truncate=graded@1"), "{off}");
    // mode on: the cacheable enabled kinds only (a root value never enters
    // the Δ cache; review M2)
    let on = with_switch(Switch::new(Mode::On, &[Kind::Root, Kind::Child]), key);
    assert_eq!(on, format!("{off}#uncontracted=child"));
    let root = with_switch(Switch::new(Mode::On, &[Kind::Root]), key);
    assert_eq!(root, off, "on:root leaves the key alone");
    // check mode: every value is today's, the key too
    let check = with_switch(Switch::new(Mode::Check, &[Kind::Child, Kind::Root]), key);
    assert_eq!(check, off, "check mode leaves the key alone");
    let inert = with_switch(Switch::new(Mode::On, &[Kind::Quotient]), key);
    assert_eq!(inert, off);
    assert_eq!(
        with_switch(Switch::OFF, uncontracted::cache_key_suffix),
        None
    );
}

/// The designed refusals are a closed list (U7a review B1): the emitter's
/// codes, and vakint's section-8 rows by variant AND detail prefix; every
/// other refusal (the entry's wrapped internal failures, the reader's
/// syntax errors, our own bugs) is an error.
#[test]
fn designed_refusals_are_a_closed_list() {
    use vakint::p_integral::factorized::FactorizedError as E;
    let closure = |detail: &str| E::Closure {
        detail: detail.into(),
    };
    for (error, code) in [
        (closure("Taylor degree 3 above 2"), "vakint:taylor_degree"),
        (
            closure(
                "the projector contains Q(0, ..): with a Taylor step it acts after the step and must be momentum-free",
            ),
            "vakint:projector_momentum",
        ),
        (
            E::Ports {
                detail: "a compact slash gamma(a, b, Q(x, mink(D))) in the projector".into(),
            },
            "vakint:compact_slash",
        ),
        (
            closure(
                "open slot bis(4,hedge(3,1)) is not a Lorentz slot (open bispinor ports are not supported; close the line with the projector)",
            ),
            "vakint:open_spinor",
        ),
        (
            E::Insertion {
                vertex: 3,
                detail: "contains a mass symbol gammalooprs::rstar_mass_squared".into(),
            },
            "vakint:insertion_mass_symbol",
        ),
        (
            closure("every line is massive: not a closure"),
            "vakint:all_massive",
        ),
        (
            E::Gamma5 {
                location: "edge 3".into(),
                tensor: "gamma5".into(),
            },
            "vakint:gamma5",
        ),
    ] {
        let decline = Decline::vakint(&error);
        assert!(decline.designed, "{error}");
        assert_eq!(decline.code, code, "{error}");
    }
    for (error, code) in [
        // the entry's internal failures (closure.rs 1212 / 1294 / 1317 / 1334)
        (closure("lowering: no key plan"), "vakint:Closure"),
        (
            closure("truncating the insertion product: x"),
            "vakint:Closure",
        ),
        (
            E::Insertion {
                vertex: 2,
                detail: "half-edge 9 is not in the graph".into(),
            },
            "vakint:Insertion",
        ),
        (
            E::Ports {
                detail: "port x does not occur in the factors".into(),
            },
            "vakint:Ports",
        ),
        (
            E::TextSyntax {
                line: 3,
                detail: "x".into(),
            },
            "vakint:TextSyntax",
        ),
        (
            E::MomentumNotConserved {
                vertex: 1,
                residual: vec![1],
            },
            "vakint:MomentumNotConserved",
        ),
        (
            E::Evaluation {
                detail: "vacuum: 2 massive lines among 3 lines".into(),
            },
            "vakint:Evaluation",
        ),
    ] {
        let decline = Decline::vakint(&error);
        assert!(!decline.designed, "{error}");
        assert_eq!(decline.code, code, "{error}");
    }
    for code in [
        "routing",
        "missing_child",
        "product_mismatch",
        "text",
        "taylor",
        "closure_data_missing",
        "normalisation_mismatch",
        "finish_tensor",
        "vakint:Empty",
    ] {
        assert!(
            !DESIGNED_LOCAL.contains(&code),
            "{code} is a bug, not designed"
        );
    }
    for code in ["open_root", "insertion_mass_symbol", "root_crown_paired"] {
        assert!(DESIGNED_LOCAL.contains(&code), "{code}");
    }
}

#[test]
fn unknown_master_text_is_recognised() {
    let unknown = vakint::UnknownMasterCoefficients {
        backend: "rustred".into(),
        requested_order: Some(0),
        coefficients: vec![vakint::UnknownMasterCoefficient {
            master: "p5_m1".into(),
            order: 1,
            first_affected_result_order: 0,
        }],
    };
    assert!(is_unknown_master_text(&unknown.to_string()), "{unknown}");
    // as vakint's factorised entry wraps it
    let wrapped = vakint::p_integral::factorized::FactorizedError::Evaluation {
        detail: unknown.to_string(),
    };
    assert!(is_unknown_master_text(&wrapped.to_string()));
    let missing = vakint::MissingMasterOrders {
        backend: "rustred".into(),
        requested_order: Some(0),
        masters: Vec::new(),
    };
    assert!(is_unknown_master_text(&missing.to_string()), "{missing}");
    assert!(!is_unknown_master_text(
        "closure: degree 3 exceeds the maximum Taylor degree"
    ));
}

// ---------------------------------------------------------------------------
// Record pieces
// ---------------------------------------------------------------------------

#[test]
fn split_matches_the_converter() {
    test_initialise().unwrap();
    for (what, text) in [
        (
            "k1-quick L4-42 rr2",
            include_str!("testdata/u7a/k1-quick-L4-42-rr2.input.txt"),
        ),
        (
            "k2-child-ins",
            include_str!("testdata/u7a/k2-child-ins.input.txt"),
        ),
    ] {
        let record = FactorizedNumerator::from_text(text).unwrap_or_else(|e| panic!("{what}: {e}"));
        assert!(!record.insertions.is_empty(), "{what}");
        for insertion in &record.insertions {
            let delta = (&insertion.coefficient * &insertion.factor).expand();
            let (coefficient, factor) = split_insertion(&delta);
            assert!(
                (&coefficient * &factor - &delta).expand().is_zero(),
                "{what}: coefficient x factor is not the counterterm"
            );
            let ratio = &coefficient / &insertion.coefficient;
            let ratio = ratio.to_canonical_string();
            assert!(
                ratio
                    .trim_start_matches('-')
                    .bytes()
                    .all(|b| b.is_ascii_digit() || b == b'/'),
                "{what}: coefficient {coefficient} is not a rational multiple of the converter's {}",
                insertion.coefficient
            );
            assert_ne!(
                coefficient,
                Atom::num(1),
                "{what}: the converter's split is non-trivial"
            );
        }
    }
    // a counterterm whose terms' scalars are not rational multiples: 1 x Δ
    let g = parse(
        "spenso::g(spenso::mink(gammalooprs::dim,gammalooprs::hedge(1,1)),spenso::mink(gammalooprs::dim,gammalooprs::hedge(2,1)))",
    );
    let delta = (Atom::var(symbol!("u7a_test_a")) + Atom::var(symbol!("u7a_test_b")) * &g).expand();
    let (coefficient, factor) = split_insertion(&delta);
    assert_eq!(coefficient, Atom::num(1));
    assert_eq!(factor, delta);
}

fn edge(id: usize, source: Endpoint, target: Endpoint, momentum: &[i64]) -> FactorizedEdge {
    FactorizedEdge {
        id,
        source,
        target,
        source_hedge: None,
        target_hedge: None,
        particle: "Q".into(),
        line: LineKind::Vector,
        power: 1,
        mass_squared: Atom::Zero,
        momentum: momentum.to_vec(),
        factor: Atom::num(1),
    }
}

#[test]
fn conservation_removes_zero_mass_logs() {
    test_initialise().unwrap();
    let q = |e: usize, h: usize| {
        function!(
            GS.emr_mom,
            e,
            function!(
                symbol!("spenso::mink"),
                GS.dim,
                function!(GS.hedgeaind, h, 1)
            )
        )
    };
    let mass = Atom::var(symbol!("rstar_mass_squared"; Scalar));
    let log = function!(symbolica::atom::Symbol::LOG, &mass);
    let g = function!(
        symbol!("spenso::g"),
        function!(
            symbol!("spenso::mink"),
            GS.dim,
            function!(GS.hedgeaind, 7, 1)
        ),
        function!(
            symbol!("spenso::mink"),
            GS.dim,
            function!(GS.hedgeaind, 8, 1)
        )
    );
    // node 5: edge 1 in (k1), edge 2 out (k1 - p1), edge 3 out (p1)
    let v = Endpoint::Vertex(5);
    let edges = [
        edge(1, Endpoint::Vertex(0), v, &[1, 0]),
        edge(2, v, Endpoint::Vertex(1), &[1, -1]),
        edge(3, v, Endpoint::External, &[0, 1]),
    ];
    let delta = (Atom::num(21) / 2048
        * &log
        * (q(1, 9) - q(2, 9) - q(3, 9))
        * Atom::var(symbol!("u7a_test_c"))
        + &g * Atom::var(symbol!("u7a_test_d")))
    .expand();
    assert!(delta.contains_symbol(mass.get_symbol().unwrap()));
    let rewritten = conserve(&delta, 5, &edges).expect("a conserving node");
    assert!(
        !rewritten.contains_symbol(mass.get_symbol().unwrap()),
        "the zero-mass log survives: {rewritten}"
    );
    assert!(
        (&rewritten - &g * Atom::var(symbol!("u7a_test_d")))
            .expand()
            .is_zero()
    );
    // a K2 crown vertex where independent legs meet: unchanged (None)
    let independent = [
        edge(1, Endpoint::External, v, &[0, 1, 0]),
        edge(2, Endpoint::External, v, &[0, 0, 1]),
        edge(3, v, Endpoint::Vertex(1), &[1, 0, 0]),
    ];
    assert!(conserve(&delta, 5, &independent).is_none());
}

/// The real case (WH-OQ-ANSWERS-ADDENDUM OQ-28): the Δ of the child
/// {4,5,6,7,8} at the root remainder {2,3} of L3_016 as the production run's
/// Δ cache spelled it (runs/L3 dump a_root-000718): `log(rstar_mass_squared)`
/// and `log(μᵣ²)` times `Q(1) - Q(2) - Q(3)`, zero by momentum conservation at
/// the ct: node (edges 2 and 3 in, leg 1 out; routes k, p - k, +p). The step
/// removes every mass symbol and leaves the routed Δ unchanged.
#[test]
fn conservation_removes_the_l3_016_zero() {
    test_initialise().unwrap();
    let delta = parse(include_str!("testdata/u7a/l3-016-delta.txt").trim());
    let mass = symbol!("rstar_mass_squared"; Scalar);
    assert!(
        delta.contains_symbol(mass),
        "the dumped Δ carries log(M_aux^2)"
    );
    let node = Endpoint::Vertex(3);
    let edges = [
        edge(1, node, Endpoint::External, &[0, 1]),
        edge(2, Endpoint::Vertex(0), node, &[1, 0]),
        edge(3, Endpoint::Vertex(0), node, &[-1, 1]),
    ];
    let rewritten = conserve(&delta, 3, &edges).expect("the ct: node conserves");
    assert!(
        !rewritten.contains_symbol(mass),
        "a mass symbol survives: {rewritten}"
    );
    // routed (Q(1) -> Q(0) by sigma, Q(3) -> Q(0) - Q(2)), both are the same
    let route = |atom: &Atom| {
        atom.replace(function!(GS.emr_mom, 1, W_.x___))
            .with(function!(GS.emr_mom, 0, W_.x___))
            .replace(function!(GS.emr_mom, 3, W_.x___))
            .with(function!(GS.emr_mom, 0, W_.x___) - function!(GS.emr_mom, 2, W_.x___))
            .expand()
    };
    assert!(
        (route(&delta) - route(&rewritten)).expand().is_zero(),
        "the step changed the routed Δ"
    );
    assert!(!route(&delta).contains_symbol(mass));
}

#[test]
fn root_factor_split() {
    test_initialise().unwrap();
    let p2 = function!(
        SPENSO_TAG.dot,
        function!(GS.emr_mom, 0, minkowski()),
        function!(GS.emr_mom, 0, minkowski())
    );
    let t = Atom::var(symbol!("beta_projection::transverse"));
    let c = &t / (Atom::i() * &p2) / (Atom::var(GS.dim) - 1);
    let (k, s) = split_root_factor(&c, 0).unwrap();
    assert_eq!(k, -1);
    let theirs = parse("(-1+gammalooprs::{}::dim)^(-1)*(-𝑖)*beta_projection::{}::transverse");
    assert!((&s - &theirs).expand().is_zero(), "{s}");
    let (k, s) = split_root_factor(&Atom::num(1), 0).unwrap();
    assert_eq!((k, s), (0, Atom::num(1)));
    let momentum = function!(GS.emr_mom, 0, minkowski());
    assert!(split_root_factor(&(&t * momentum), 0).is_err());
}

#[test]
fn comparator_catches_mutations() {
    test_initialise().unwrap();
    for (what, text) in [
        ("k1-root", include_str!("testdata/u7a/k1-root.input.txt")),
        (
            "k1-quick L4-42 rr2",
            include_str!("testdata/u7a/k1-quick-L4-42-rr2.input.txt"),
        ),
    ] {
        let record = FactorizedNumerator::from_text(text).unwrap();
        assert!(record_differences(&record, &record).is_empty(), "{what}");
        // a re-labelled copy still agrees (the label is ignored)
        let mut relabelled = record.clone();
        relabelled.label = "other".into();
        assert!(record_differences(&relabelled, &record).is_empty());
        let mutations: Vec<(&str, Box<dyn Fn(&mut FactorizedNumerator)>)> = vec![
            (
                "a wrong vector",
                Box::new(|r: &mut FactorizedNumerator| {
                    let e = r.edges.iter_mut().find(|e| !e.is_external()).unwrap();
                    e.momentum[0] += 1;
                }),
            ),
            (
                "a dropped (-i) in the scalars",
                Box::new(|r: &mut FactorizedNumerator| {
                    r.prefactor.scalars = &r.prefactor.scalars * Atom::i();
                }),
            ),
            (
                "a wrong p2_power",
                Box::new(|r: &mut FactorizedNumerator| {
                    r.projector.as_mut().unwrap().p2_power += 1;
                }),
            ),
            (
                "a wrong vertex factor",
                Box::new(|r: &mut FactorizedNumerator| {
                    let v = r
                        .vertices
                        .iter_mut()
                        .find(|v| !v.rule.starts_with("ct:"))
                        .unwrap();
                    v.factor = &v.factor * Atom::num(2);
                }),
            ),
        ];
        for (mutation, apply) in &mutations {
            let mut mutated = record.clone();
            apply(&mut mutated);
            assert!(
                !record_differences(&mutated, &record).is_empty(),
                "{what}: {mutation} is not caught"
            );
        }
        if !record.insertions.is_empty() {
            let mut removed = record.clone();
            removed.insertions.pop();
            assert!(
                !record_differences(&removed, &record).is_empty(),
                "{what}: a removed insertion"
            );
            let mut hedges = record.clone();
            hedges.insertions[0].hedges.pop();
            assert!(
                !record_differences(&hedges, &record).is_empty(),
                "{what}: wrong hedges"
            );
            let mut factor = record.clone();
            factor.insertions[0].coefficient = &factor.insertions[0].coefficient * Atom::num(3);
            assert!(
                !record_differences(&factor, &record).is_empty(),
                "{what}: wrong coefficient"
            );
        }
    }
}

// ---------------------------------------------------------------------------
// K1 records of real diagrams
// ---------------------------------------------------------------------------

/// One closed root with a single-record fixture: emitted record == the
/// fixture's input.txt, value == value.txt; check (every other record
/// compared, the declines stated) and on (== off) as [`checked_run`].
fn single_root_fixture(
    what: &str,
    dot: &str,
    input: &str,
    value_text: &str,
    declines: &[(&str, usize)],
) -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let theirs = FactorizedNumerator::from_text(input).map_err(|e| eyre!("{e}"))?;
        let diagram = fixture_branch(dot, &[&theirs])?;
        let (_, _, marked) = root_run(&diagram, Switch::OFF)?;
        let FactorizedKind::Closure(spec) = &theirs.kind else {
            return Err(eyre!("{what}: not a closure record"));
        };
        ensure!(
            marked == spec.edge,
            "{what}: marked edge {marked}, the fixture's closure edge {}",
            spec.edge
        );
        let captured = checked_run(what, &diagram, &Vakint::new()?, &ROOT_CHILD, declines, true)?;
        compare_fixture(what, &captured, input, value_text)
    })
}

#[test]
fn k1_root_fixture() -> Result<()> {
    single_root_fixture(
        "k1-root (L2 22)",
        include_str!("testdata/u7a/l2-d22.dot"),
        include_str!("testdata/u7a/k1-root.input.txt"),
        include_str!("testdata/u7a/k1-root.value.txt"),
        &[],
    )
}

#[test]
fn k1_root_rc_fixture() -> Result<()> {
    if std::env::var("GAMMALOOP_ROOT_ROUTING").is_ok_and(|v| v.trim() == "default") {
        return Err(eyre!(
            "k1-root-rc needs the num-est root routing (GAMMALOOP_ROOT_ROUTING unset)"
        ));
    }
    single_root_fixture(
        "k1-root-rc (L2 21, r_c = +p(1))",
        include_str!("testdata/u7a/l2-d21.dot"),
        include_str!("testdata/u7a/k1-root-rc.input.txt"),
        include_str!("testdata/u7a/k1-root-rc.value.txt"),
        &[],
    )
}

/// The r_c = 0 control of k1-root-rc: the default root routing, read once
/// per process, so run alone with `GAMMALOOP_ROOT_ROUTING=default`.
#[test]
#[ignore = "needs GAMMALOOP_ROOT_ROUTING=default in its own process"]
fn k1_root_rc0_fixture() -> Result<()> {
    if !std::env::var("GAMMALOOP_ROOT_ROUTING").is_ok_and(|v| v.trim() == "default") {
        return Err(eyre!("k1-root-rc0 needs GAMMALOOP_ROOT_ROUTING=default"));
    }
    single_root_fixture(
        "k1-root-rc0 (L2 21, r_c = 0)",
        include_str!("testdata/u7a/l2-d21.dot"),
        include_str!("testdata/u7a/k1-root-rc0.input.txt"),
        include_str!("testdata/u7a/k1-root-rc0.value.txt"),
        &[],
    )
}

/// The factor test (U7a review I1, I3): production diagrams colour factored
/// as they are, NO fixture factor injected, so the check mode's PRODUCT
/// identity tests the factor extraction; every root and child record is
/// compared equal except the stated declines (the one-loop massive tadpole
/// remainders vakint's closure entry refuses, `vakint:all_massive`); check
/// returns today's value; mode on == off on the diagram's poles.
fn equals_off(what: &str, dot: &str, declines: &[(&str, usize)]) -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let diagram = colour_factored(&load(dot)?)?;
        checked_run(what, &diagram, &Vakint::new()?, &ROOT_CHILD, declines, true)?;
        Ok(())
    })
}

#[test]
fn k1_check_equals_off() -> Result<()> {
    for (name, dot, declines) in [
        ("L2/0", include_str!("testdata/l2-d0.dot"), &[][..]),
        ("L2/8", include_str!("testdata/l2-d8.dot"), &[][..]),
        ("L2/12", include_str!("testdata/l2-d12.dot"), &[][..]),
        ("L3/185", include_str!("testdata/l3-d185.dot"), &[][..]),
    ] {
        equals_off(name, dot, declines)?;
    }
    Ok(())
}

#[test]
fn k1_check_equals_off_l4_40() -> Result<()> {
    equals_off("L4/40", include_str!("testdata/u7a/l4-d040.dot"), &[])
}

#[test]
fn k1_check_equals_off_l4_42() -> Result<()> {
    equals_off("L4/42", include_str!("testdata/u7a/l4-d042.dot"), &[])
}

#[test]
fn k1_check_equals_off_l4_49() -> Result<()> {
    equals_off(
        "L4/49",
        include_str!("testdata/u7a/l4-d049.dot"),
        &[("vakint:all_massive", 1)],
    )
}

#[test]
fn k1_check_equals_off_l4_100() -> Result<()> {
    // the first colour row: the one-loop tadpole remainder (id 8) has a zero
    // insertion (outcome zero)
    equals_off("L4/100", include_str!("testdata/u7a/l4-d100.dot"), &[])
}

/// The Δ momentum-conservation step on a real record (U7a review I4;
/// WH-RECORD-FORMAT section 8, OQ-28): L3 016 with the Δ of its child
/// {4,5,6,7,8} as the production run's Δ cache spelled it (`log(M_aux^2)`
/// and `log(μᵣ²)` times `Q(1) - Q(2) - Q(3)`, a momentum-conservation zero
/// at the `ct:` node), the other children as computed. The root remainder
/// {2,3} with that child contracted is emitted with `conserved=1` (no
/// `insertion_mass_symbol` decline), evaluated by vakint and, in check
/// mode, equal to today's value of the same Δ; mode on == off.
#[test]
fn conservation_on_l3_016() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let diagram = colour_factored(&load(include_str!("testdata/u7a/l3-d016.dot"))?)?;
        let settings = settings();
        let vakint = Vakint::new()?;
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let model = diagram.model();
        let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
        let reference_mass = Atom::var(*super::REFERENCE_MASS_SQUARED).sqrt();
        let mut children =
            plan.evaluate_uv_children(&graph, model, &reference_mass, &vakint, &settings)?;
        let support = children
            .keys()
            .find(|support| super::plan::uv_child_edge_ids(&graph, support) == [4, 5, 6, 7, 8])
            .cloned()
            .ok_or_else(|| eyre!("L3 016 has no child {{4,5,6,7,8}}"))?;
        let dumped = parse(include_str!("testdata/u7a/l3-016-delta.txt").trim());
        children.insert(support, dumped);
        let closed_root = beta_normalization(&graph);
        let callback = |tensor: &Atom| -> Result<Atom> {
            Ok(tensor.replace_multiple(&closed_root.external_routing) * &closed_root.factor)
        };
        let run = |switch: Switch| {
            let (value, captured) = with_switch(switch, || {
                with_strict(true, || {
                    with_boundary_mode(BoundaryMode::On, false, || {
                        with_capture(|| {
                            plan.projected_counterterm_closing_root_scalar_with_children(
                                &graph,
                                model,
                                &vakint,
                                &settings,
                                &callback,
                                &closed_root,
                                children.clone(),
                            )
                        })
                    })
                })
            });
            value.map(|value| (value, captured))
        };
        let (off, nothing) = run(Switch::OFF)?;
        ensure!(nothing.is_empty());
        let (checked, captured) = run(Switch::new(Mode::Check, &ROOT_CHILD))?;
        ensure!(checked == off, "check mode returns today's value");
        expect_outcomes("L3 016 check", &captured, &[], true)?;
        let conserved = captured
            .iter()
            .filter(|a| a.kind == Kind::Root && a.conserved > 0)
            .collect::<Vec<_>>();
        ensure!(
            conserved.len() == 1 && conserved[0].conserved == 1,
            "L3 016: the records with a conserved insertion: {:?}",
            captured
                .iter()
                .map(|a| (a.id, a.outcome.clone(), a.conserved))
                .collect::<Vec<_>>()
        );
        let text = conserved[0].text.as_ref().unwrap();
        let record = FactorizedNumerator::from_text(text).map_err(|e| eyre!("{e}"))?;
        ensure!(
            record_key(&record).1 == ["ct:edges-4-5-6-7-8"],
            "L3 016: the conserved record is not the remainder with {{4,5,6,7,8}} contracted\n{text}"
        );
        let mass = symbol!("rstar_mass_squared"; Scalar);
        ensure!(
            record
                .insertions
                .iter()
                .all(|i| !i.coefficient.contains_symbol(mass) && !i.factor.contains_symbol(mass)),
            "L3 016: a mass symbol survives in the record's insertion\n{text}"
        );
        let (on, attempts) = run(Switch::new(Mode::On, &ROOT_CHILD))?;
        expect_outcomes("L3 016 on", &attempts, &[], false)?;
        ensure!(
            attempts.iter().any(|a| a.conserved == 1),
            "L3 016 on: no conserved record"
        );
        ensure!(atoms_agree(&on, &off), "L3 016: mode on differs from off");
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// K2: UV children
// ---------------------------------------------------------------------------

/// k2-child: the isotropic closed child {2,3,5} of L2 22 (bare, one loop,
/// degree 0; a paired crown edge 4 = legs 4 and 100004), against
/// wh-fixtures/k2-child (record and value).
#[test]
fn k2_child_fixture() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let input = include_str!("testdata/u7a/k2-child.input.txt");
        let root = FactorizedNumerator::from_text(include_str!("testdata/u7a/k1-root.input.txt"))
            .map_err(|e| eyre!("{e}"))?;
        let child = FactorizedNumerator::from_text(input).map_err(|e| eyre!("{e}"))?;
        let diagram = fixture_branch(include_str!("testdata/u7a/l2-d22.dot"), &[&root, &child])?;
        let captured = checked_run(
            "k2-child",
            &diagram,
            &Vakint::new()?,
            &ROOT_CHILD,
            &[],
            true,
        )?;
        compare_fixture_of(
            "k2-child (L2 22, isotropic {2,3,5})",
            Kind::Child,
            &captured,
            input,
            include_str!("testdata/u7a/k2-child.value.txt"),
        )
    })
}

/// k2-child-ins: the open child {5,6,7,8} of L3 36, remainder {5,6} with
/// {7,8} contracted to the `ct:` vertex 4 (two split crown legs 3 and 4,
/// leg 3 on the merged vertex), against wh-fixtures/k2-child-ins.
#[test]
fn k2_child_ins_fixture() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let input = include_str!("testdata/u7a/k2-child-ins.input.txt");
        let (vertices, edges) = site_factors(include_str!("testdata/u7a/k2-child-ins.branch.txt"));
        let diagram = branch_with(include_str!("testdata/u7a/l3-d036.dot"), vertices, edges)?;
        let captured = checked_run(
            "k2-child-ins",
            &diagram,
            &Vakint::new()?,
            &ROOT_CHILD,
            &[],
            true,
        )?;
        compare_fixture_of(
            "k2-child-ins (L3 36, {5,6} with {7,8} contracted)",
            Kind::Child,
            &captured,
            input,
            include_str!("testdata/u7a/k2-child-ins.value.txt"),
        )
    })
}

/// `got == expected` through ε^order up to a relabelling of the isotropic
/// channel markers `uv_isotropic_channel_{1,2,3}` in `expected` (a site value
/// of another run orders a closed child's crown pairings by its own crown
/// order).
fn equal_through_channels(got: &Atom, expected: &Atom, order: i64) -> Result<bool> {
    let markers = [1, 2, 3]
        .map(|j| Atom::var(symbol!(&format!("gammalooprs::uv_isotropic_channel_{j}"); Scalar)));
    let staging = [1, 2, 3].map(|j| Atom::var(symbol!(&format!("u7a_test_channel_{j}"))));
    for permutation in [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ] {
        let mut relabelled = expected.clone();
        for (marker, stage) in markers.iter().zip(&staging) {
            relabelled = relabelled.replace(marker.clone()).with(stage.clone());
        }
        for (stage, target) in staging.iter().zip(permutation) {
            relabelled = relabelled
                .replace(stage.clone())
                .with(markers[target].clone());
        }
        if equal_through(got, &relabelled, order)? {
            return Ok(true);
        }
    }
    Ok(false)
}

/// The K2 shapes without a v3 fixture, against the raw site values
/// (wh-fixtures/b_uv_child: the contracted path's value at the same call,
/// through ε^-1): an open child (L2 22 {4,5}, degree 1), an isotropic child
/// (L2 25 {2,4}), a Dirac self-energy with its open ν (L2 26 {5,6}), a Dirac
/// vertex (L2 28 {4,5,6}) and a closed two-loop child through the child
/// boundary (L3 352 {2,3,4,5,6}). Each record is emitted, evaluated by
/// vakint and equal to the site value (the channel markers up to their
/// order); every other record is compared equal except the stated declines;
/// mode on agrees with off.
#[test]
fn k2_shapes_values() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let none: &[(&str, usize)] = &[];
        for (what, closure, dot, site, declines) in [
            (
                "open (L2 22 {4,5})",
                "none",
                include_str!("testdata/u7a/l2-d22.dot"),
                include_str!("testdata/u7a/k2-open.site.txt"),
                none,
            ),
            (
                "isotropic (L2 25 {2,4})",
                "isotropic",
                include_str!("testdata/u7a/l2-d25.dot"),
                include_str!("testdata/u7a/k2-isotropic.site.txt"),
                &[("vakint:all_massive", 1)][..],
            ),
            (
                "dirac self-energy (L2 26 {5,6})",
                "dirac_self_energy",
                include_str!("testdata/u7a/l2-d26.dot"),
                include_str!("testdata/u7a/k2-dirac-self-energy.site.txt"),
                none,
            ),
            (
                "dirac vertex (L2 28 {4,5,6})",
                "dirac_vertex",
                include_str!("testdata/u7a/l2-d28.dot"),
                include_str!("testdata/u7a/k2-dirac-vertex.site.txt"),
                none,
            ),
            (
                "closed child boundary (L3 352 {2,3,4,5,6})",
                "isotropic",
                include_str!("testdata/u7a/l3-d352.dot"),
                include_str!("testdata/u7a/k2-closed-child-boundary.site.txt"),
                none,
            ),
        ] {
            let (vertices, edges) = site_factors(site);
            let diagram = branch_with(dot, vertices, edges)?;
            let captured = checked_run(what, &diagram, &vakint, &ROOT_CHILD, declines, true)?;
            let field = |key: &str| {
                site.lines()
                    .find_map(|line| line.strip_prefix(key)?.strip_prefix('\t'))
                    .ok_or_else(|| eyre!("{what}: no {key} line"))
            };
            let remaining = field("remaining_edges")?
                .split(',')
                .map(|edge| edge.parse::<usize>())
                .collect::<std::result::Result<Vec<_>, _>>()?;
            let expected = parse(field("value")?);
            let mut found = 0;
            for attempt in captured.iter().filter(|a| a.kind == Kind::Child) {
                let Some(text) = &attempt.text else { continue };
                let ours = FactorizedNumerator::from_text(text).map_err(|e| eyre!("{e}"))?;
                if record_key(&ours) != (remaining.clone(), Vec::new())
                    || !ours.label.contains(&format!("closure={closure} "))
                {
                    continue;
                }
                found += 1;
                ensure!(
                    attempt.outcome == "uncontracted",
                    "{what}: outcome {}\n{text}",
                    attempt.outcome
                );
                let value = attempt.value.as_ref().ok_or_else(|| eyre!("no value"))?;
                if let Some(contracted) = &attempt.contracted {
                    ensure!(
                        equal_through_channels(contracted, &expected, -1)?,
                        "{what}: today's value differs from the site value (not this branch)\ntoday {contracted}\nsite  {expected}"
                    );
                }
                ensure!(
                    equal_through_channels(value, &expected, -1)?,
                    "{what}: the value differs from the site value through eps^-1\nours {value}\nsite {expected}\n{text}"
                );
                eprintln!(
                    "U7A_K2 {what} id={} identical={:?} equal={:?}",
                    attempt.id, attempt.identical, attempt.equal
                );
            }
            ensure!(
                found == 1,
                "{what}: {found} records of the remainder {remaining:?} closure={closure}"
            );
        }
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// The k1-quick and k1-l5 sets (read in place, #[ignore])
// ---------------------------------------------------------------------------

/// One VERIFIED K1 set read in place under `GAMMALOOP_U7A_FIXTURES_DIR`: per
/// colour branch (its first record is the bare root, which states every
/// vertex factor of the branch), one run of the diagram ([`checked_run`]:
/// check, and with `on` mode on == off; `declines` per branch); every record
/// of the branch emitted and equal to its `input.txt`, its value equal to
/// `value.txt`. The number of records compared.
fn k1_fixture_set(
    set: &str,
    dot: &str,
    branches: &[Vec<String>],
    p_library: bool,
    declines: &[(&str, usize)],
    on: bool,
) -> Result<usize> {
    let root = fixtures_dir()?.join(set);
    ensure!(
        root.is_dir(),
        "{set}: no such set under GAMMALOOP_U7A_FIXTURES_DIR ({})",
        root.display()
    );
    let read = |path: std::path::PathBuf| {
        std::fs::read_to_string(&path).wrap_err_with(|| format!("{}", path.display()))
    };
    let mut vakint = Vakint::new()?;
    if p_library {
        let library = std::env::var("GAMMALOOP_U7A_P_LIBRARY")
            .or_else(|_| std::env::var("BETA_P_LIBRARY"))
            .map_err(|_| {
                eyre!("{set} needs the beta's p library: GAMMALOOP_U7A_P_LIBRARY or BETA_P_LIBRARY")
            })?;
        ensure!(
            std::path::Path::new(&library).is_dir(),
            "{set}: the p library {library} is not a directory"
        );
        vakint.load_p_library_directory(&library, false)?;
    }
    let dot = read(root.join(dot))?;
    let mut records = 0;
    for branch in branches {
        let bare = read(root.join(&branch[0]).join("input.txt"))?;
        let bare = FactorizedNumerator::from_text(&bare).map_err(|e| eyre!("{e}"))?;
        let diagram = fixture_branch(&dot, &[&bare])?;
        let started = std::time::Instant::now();
        let captured = checked_run(
            &format!("{set} {}", branch[0]),
            &diagram,
            &vakint,
            &ROOT_CHILD,
            declines,
            on,
        )?;
        eprintln!(
            "U7A_SET {set} {} seconds={:.1} attempts={}",
            branch[0],
            started.elapsed().as_secs_f64(),
            captured.len()
        );
        for record in branch {
            let dir = root.join(record);
            compare_fixture(
                &format!("{set} {record}"),
                &captured,
                &read(dir.join("input.txt"))?,
                &read(dir.join("value.txt"))?,
            )?;
            records += 1;
        }
    }
    Ok(records)
}

fn numbered(prefix: &str, n: usize) -> Vec<String> {
    (1..=n).map(|i| format!("{prefix}{i}")).collect()
}

/// k1-quick L4-40: seven boundary remainders of one branch.
#[test]
#[ignore = "needs GAMMALOOP_U7A_FIXTURES_DIR (the beta's wh-fixtures)"]
fn k1_quick_l4_40_fixtures() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let n = k1_fixture_set(
            "k1-quick/L4-40",
            "rr1/raw/diagram-000040.dot",
            &[numbered("rr", 7)],
            false,
            &[],
            true,
        )?;
        ensure!(n == 7);
        Ok(())
    })
}

/// k1-quick L4-42: seven boundary remainders of one branch.
#[test]
#[ignore = "needs GAMMALOOP_U7A_FIXTURES_DIR (the beta's wh-fixtures)"]
fn k1_quick_l4_42_fixtures() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let n = k1_fixture_set(
            "k1-quick/L4-42",
            "rr1/raw/diagram-000042.dot",
            &[numbered("rr", 7)],
            false,
            &[],
            true,
        )?;
        ensure!(n == 7);
        Ok(())
    })
}

/// k1-quick L4-49: seven boundary remainders of one branch (directories
/// `1`..`7`).
#[test]
#[ignore = "needs GAMMALOOP_U7A_FIXTURES_DIR (the beta's wh-fixtures)"]
fn k1_quick_l4_49_fixtures() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let n = k1_fixture_set(
            "k1-quick/L4-49",
            "1/raw/diagram-000049.dot",
            &[numbered("", 7)],
            false,
            &[("vakint:all_massive", 1)],
            true,
        )?;
        ensure!(n == 7);
        Ok(())
    })
}

/// k1-quick L4-100: both root evaluations (colour branches b0 / b1), seven
/// boundary remainders each; per branch the one-loop massive tadpole
/// remainder (id 8) is vakint's `all_massive` refusal.
#[test]
#[ignore = "needs GAMMALOOP_U7A_FIXTURES_DIR (the beta's wh-fixtures)"]
fn k1_quick_l4_100_fixtures() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let n = k1_fixture_set(
            "k1-quick/L4-100",
            "b0-rr1/raw/diagram-000100.dot",
            &[numbered("b0-rr", 7), numbered("b1-rr", 7)],
            false,
            &[("vakint:all_massive", 1)],
            true,
        )?;
        ensure!(n == 14);
        Ok(())
    })
}

/// k1-l5 L5-60 (the cheap five-loop set): fifteen boundary remainders, the
/// four-loop cut through the beta's p library. Check mode only (about ten
/// minutes on one core).
#[test]
#[ignore = "needs GAMMALOOP_U7A_FIXTURES_DIR and the p library; about 10 minutes"]
fn k1_l5_60_fixtures() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let n = k1_fixture_set(
            "k1-l5/L5-60",
            "rr1/raw/diagram-000060.dot",
            &[numbered("rr", 15)],
            true,
            &[],
            false,
        )?;
        ensure!(n == 15);
        Ok(())
    })
}

/// k1-l5 L5-51 (heavy: about 20 minutes and 7 GB): six boundary remainders.
#[test]
#[ignore = "heavy: L5 51 (the bare root's reduction is about 330 s, the diagram 6 min+)"]
fn k1_l5_51_fixtures() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let n = k1_fixture_set(
            "k1-l5/L5-51",
            "rr1/raw/diagram-000051.dot",
            &[numbered("rr", 6)],
            true,
            &[],
            false,
        )?;
        ensure!(n == 6);
        Ok(())
    })
}

/// k1-l5 L5-34342 (heavy: a closed quark loop, 77 s contraction): the bare
/// root, the one boundary remainder.
#[test]
#[ignore = "heavy: L5 34342 (77 s contraction, 3.4 GB)"]
fn k1_l5_34342_fixtures() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let n = k1_fixture_set(
            "k1-l5/L5-34342",
            "rr1/raw/diagram-034342.dot",
            &[numbered("rr", 1)],
            true,
            &[],
            false,
        )?;
        ensure!(n == 1);
        Ok(())
    })
}
