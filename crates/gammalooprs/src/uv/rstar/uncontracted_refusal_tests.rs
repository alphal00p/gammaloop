//! U7a delta, lane B (what a record's outcome is): designed refusals against
//! errors, and the scopes that decide a record's kind (Astra U7a verdict
//! findings 2 and 3, `keep the rebase test narrow`; helper VH round 2 R2-I3,
//! R2-I4, R2-M4, R2-M5; U7-PLAN requirements 13 and 15). Rule: ONLY a named
//! designed refusal falls back to today's path; everything else is
//! `outcome=error`, loud in check and in on (strict is their only behaviour,
//! R2-I1); in the merged delta the ledger (lane A) reports it as the kind's
//! `error=<code>:n outcome=error`, and a designed refusal as `declined=<code>:n`.
//!
//! * A nonzero residual of the Δ momentum-conservation step in a component
//!   the record must conserve is a routing error, never
//!   `insertion_mass_symbol`: the signed mutation at the full hook (check and
//!   on fail), with its controls (a genuine mass dependence: the
//!   designed decline; independent crown legs: the record is conserved and
//!   evaluated).
//! * A panic caught inside the dressed scope (or the quotient scope) leaves
//!   no scope behind: the same thread's next root and child are classified
//!   and executed as root and child.
//! * `one_loop_tree_cut` is the emitter's own decline of the single
//!   propagator shape only (the closing line alone, a massive tadpole after
//!   the Taylor step); every other one-loop record stays on the new path.
//! * A vakint refusal is designed only by a row of `DESIGNED_VAKINT`, a
//!   subset of the frozen REFUSALS list (U7-PLAN requirements 13 and 14):
//!   four `Closure` rows and the `all_massive` backstop; every other
//!   `Closure` detail is an error; `Evaluation` and `Oracle` always are.
//! * The record builder's impossible states (`sigma`, `root_crown_paired`,
//!   `root_normalisation_momentum`) are errors.
//! * vakint's `open slot ` row needs `is not a Lorentz slot`.
//! * A K2 / K3 crown mismatch is its own error `crown` (not `routing`); the
//!   second leg of a paired crown edge has the record id
//!   `PAIRED_LEG_OFFSET + id`, every edge id below the offset.
//! * Requirement 15 on K1 / K2: the depth check passes the real ε^0-only K2
//!   records (U7a-EPS0-EVIDENCE) and fails a pole difference.

use std::{cell::RefCell, collections::BTreeMap, rc::Rc};

use eyre::{Result, ensure, eyre};
use linnet::half_edge::{involution::EdgeIndex, subgraph::SuBitGraph};
use spenso::{
    network::tags::SPENSO_TAG,
    structure::representation::{Minkowski, RepName},
};
use symbolica::{
    atom::{Atom, AtomCore},
    function, symbol,
};
use vakint::{
    Vakint,
    p_integral::factorized::{Endpoint, FactorizedEdge, FactorizedError, FactorizedNumerator},
};

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt, graph::Graph, initialisation::test_initialise,
    utils::GS, uv::approx::streamed::with_production_piece_limit,
};

use super::{
    REFERENCE_MASS_SQUARED, RStarPlan,
    boundary::{BoundaryMode, ClosedRootNormalization, with_boundary_mode},
    expression::RStarIntegrand,
    plan::{join_edge_ids, uv_child_edge_ids},
    recursive_tests::{colour_factored, load, settings},
    top_mass::TopMassTerms,
    uncontracted::{
        self, Captured, DESIGNED_LOCAL, Decline, InsertionHook, Kind, Mode, Outcome, Site, Switch,
        split_root_factor, with_capture, with_insertion_hook, with_strict, with_switch,
    },
    uncontracted_ledger as ledger,
    uncontracted_oracle_tests::{atoms_equal, auxiliary_mass, metric_root, parse},
};

const ROOT_CHILD: [Kind; 2] = [Kind::Root, Kind::Child];

/// The `<key>=` value of a ledger line.
fn field<'l>(line: &'l str, key: &str) -> Option<&'l str> {
    line.split_whitespace()
        .find_map(|token| token.strip_prefix(key)?.strip_prefix('='))
}

/// The one `UNCONTRACTED_DISPOSITION` line of `kind` and the one
/// `UNCONTRACTED_SUMMARY` line among the ledger lines `lines` (lane A's
/// per-diagram ledger: a root evaluation outside a diagram job opens its own).
fn disposition_and_summary(lines: &[String], kind: Kind) -> Result<(String, String)> {
    let of = |head: &str| {
        lines
            .iter()
            .filter(|line| line.starts_with(head))
            .filter(|line| {
                head != "UNCONTRACTED_DISPOSITION " || field(line, "kind") == Some(kind.as_str())
            })
            .cloned()
            .collect::<Vec<_>>()
    };
    let (dispositions, summaries) = (of("UNCONTRACTED_DISPOSITION "), of("UNCONTRACTED_SUMMARY "));
    ensure!(
        dispositions.len() == 1 && summaries.len() == 1,
        "{} {} disposition and {} summary lines: {lines:#?}",
        dispositions.len(),
        kind.as_str(),
        summaries.len()
    );
    Ok((dispositions[0].clone(), summaries[0].clone()))
}

fn check() -> Switch {
    Switch::new(Mode::Check, &ROOT_CHILD)
}

fn on() -> Switch {
    Switch::new(Mode::On, &ROOT_CHILD)
}

/// `rstar_mass_squared` (M_aux²), the mass symbol of the dumped L3 016 Δ.
fn auxiliary_symbol() -> symbolica::atom::Symbol {
    symbol!("rstar_mass_squared"; Scalar)
}

fn mink() -> Atom {
    Minkowski {}.new_rep(GS.dim).to_symbolic([])
}

/// The signed incidence of a record edge at `node` (in +1, out -1, 0 when
/// not incident or a self-loop).
fn incidence(edge: &FactorizedEdge, node: usize) -> i64 {
    i64::from(edge.target == Endpoint::Vertex(node))
        - i64::from(edge.source == Endpoint::Vertex(node))
}

/// The momentum residual of the record edges at `node`.
fn residual(edges: &[FactorizedEdge], node: usize) -> Vec<i64> {
    let width = edges.first().map_or(0, |edge| edge.momentum.len());
    let mut out = vec![0i64; width];
    for edge in edges {
        let sign = incidence(edge, node);
        for (slot, c) in out.iter_mut().zip(&edge.momentum) {
            *slot += sign * c;
        }
    }
    out
}

/// Every attempt evaluated (or a zero insertion) and, in check mode,
/// compared equal with its first difference above its order (requirement
/// 15), except the attempts whose outcome is in `declined`.
fn evaluated_except(what: &str, captured: &[Captured], declined: &[&str]) -> Result<()> {
    for attempt in captured {
        eprintln!(
            "U7AD_B {what} kind={} id={} loops={} outcome={} identical={:?} equal={:?} first={:?} conserved={}",
            attempt.kind.as_str(),
            attempt.id,
            attempt.loops,
            attempt.outcome,
            attempt.identical,
            attempt.equal,
            attempt.first_difference,
            attempt.conserved
        );
        if declined.contains(&attempt.outcome.as_str()) {
            continue;
        }
        ensure!(
            attempt.outcome == "uncontracted" || attempt.outcome == "zero",
            "{what}: {} id {} outcome {}",
            attempt.kind.as_str(),
            attempt.id,
            attempt.outcome
        );
        if let Some(equal) = attempt.equal {
            use uncontracted::FirstDifference;
            ensure!(
                equal,
                "{what}: {} id {} not equal",
                attempt.kind.as_str(),
                attempt.id
            );
            let recorded = match (attempt.identical, attempt.first_difference) {
                (Some(_), Some(FirstDifference::None)) => true,
                (Some(false), Some(FirstDifference::At(k))) => k > attempt.order,
                _ => false,
            };
            ensure!(
                recorded,
                "{what}: {} id {} first difference {:?} (order {})",
                attempt.kind.as_str(),
                attempt.id,
                attempt.first_difference,
                attempt.order
            );
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// U7a.2: a nonzero loop residual of the conservation step is a routing error
// ---------------------------------------------------------------------------

const L3_16_DOT: &str = include_str!("testdata/u7a/l3-d016.dot");
const L3_36_DOT: &str = include_str!("testdata/u7a/l3-d036.dot");
const L3_016_DELTA: &str = include_str!("testdata/u7a/l3-016-delta.txt");

/// L3 016 with the production Δ of its child {4,5,6,7,8} (`log(M_aux²)` and
/// `log(μᵣ²)` times `Q(1) - Q(2) - Q(3)`, a zero by momentum conservation at
/// the `ct:` node of the root remainder {2,3}), the other children as
/// computed (`conservation_on_l3_016`).
struct L3016 {
    diagram: feynkit_graph::FeynmanDiagram,
    graph: Graph,
    plan: RStarPlan,
    children: BTreeMap<SuBitGraph, Atom>,
    support: SuBitGraph,
    delta: Atom,
}

fn l3_016() -> Result<L3016> {
    let diagram = colour_factored(&load(L3_16_DOT)?)?;
    let settings = settings();
    let vakint = Vakint::new()?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
    let reference_mass = Atom::var(*REFERENCE_MASS_SQUARED).sqrt();
    let children = plan.evaluate_uv_children(&graph, model, &reference_mass, &vakint, &settings)?;
    let support = children
        .keys()
        .find(|support| uv_child_edge_ids(&graph, support) == [4, 5, 6, 7, 8])
        .cloned()
        .ok_or_else(|| eyre!("L3 016 has no child {{4,5,6,7,8}}"))?;
    let delta = parse(L3_016_DELTA.trim())?;
    ensure!(delta.contains_symbol(auxiliary_symbol()));
    Ok(L3016 {
        diagram,
        graph,
        plan,
        children,
        support,
        delta,
    })
}

/// L3 016's root under `switch` (strict) with `delta` as the Δ of
/// {4,5,6,7,8}: its value and every attempt.
fn l3_016_run(case: &L3016, delta: &Atom, switch: Switch) -> Result<(Result<Atom>, Vec<Captured>)> {
    let settings = settings();
    let vakint = Vakint::new()?;
    let model = case.diagram.model();
    let closed_root = metric_root(&case.graph);
    let callback = |tensor: &Atom| -> Result<Atom> {
        Ok(tensor.replace_multiple(&closed_root.external_routing) * &closed_root.factor)
    };
    let mut children = case.children.clone();
    children.insert(case.support.clone(), delta.clone());
    Ok(with_switch(switch, || {
        with_strict(true, || {
            with_boundary_mode(BoundaryMode::On, false, || {
                with_capture(|| {
                    case.plan
                        .projected_counterterm_closing_root_scalar_with_children(
                            &case.graph,
                            model,
                            &vakint,
                            &settings,
                            &callback,
                            &closed_root,
                            children,
                        )
                })
            })
        })
    }))
}

/// The step's classification of a residual: zero, a rewrite; nonzero in a
/// component the record must conserve (`conserved` leading components: the
/// loop ones, and the externals of a K1 record), the routing error; nonzero
/// only in an independent external component (K2 / K3 crown legs),
/// `Independent` (no relation, unchanged).
#[test]
fn conservation_classifies_the_residual() {
    test_initialise().unwrap();
    let q = |e: usize| function!(GS.emr_mom, e, mink());
    let delta = (function!(SPENSO_TAG.dot, q(3) - q(4), q(1))
        * function!(symbolica::atom::Symbol::LOG, Atom::var(auxiliary_symbol())))
    .expand();
    let edge = |id: usize, source: Endpoint, target: Endpoint, momentum: &[i64]| FactorizedEdge {
        id,
        source,
        target,
        source_hedge: None,
        target_hedge: None,
        particle: "g".into(),
        line: vakint::p_integral::factorized::LineKind::Vector,
        power: 1,
        mass_squared: Atom::Zero,
        momentum: momentum.to_vec(),
        factor: Atom::num(1),
    };
    let v = Endpoint::Vertex(5);
    // node 5: leg 1 in (p1), leg 2 in (p2), edge 3 in (k), edge 4 out
    let node = |out: &[i64]| {
        [
            edge(1, Endpoint::External, v, &[0, 1, 0]),
            edge(2, Endpoint::External, v, &[0, 0, 1]),
            edge(3, Endpoint::Vertex(0), v, &[1, 0, 0]),
            edge(4, v, Endpoint::Vertex(1), out),
        ]
    };
    let code = |result: std::result::Result<uncontracted::Conservation, Decline>| match result {
        Err(decline) => format!("error:{}:{}", decline.code, decline.designed),
        Ok(uncontracted::Conservation::Rewritten(_)) => "rewritten".to_owned(),
        Ok(other) => format!("{other:?}"),
    };
    for (out, conserved, expected) in [
        // k + p1 + p2 out: conserves in every component
        (&[1, 1, 1][..], 3, "rewritten"),
        (&[1, 1, 1][..], 1, "rewritten"),
        // k + p1 out: p2 left over, an independent external (K2 / K3) ...
        (&[1, 1, 0][..], 1, "Independent"),
        // ... but a conserved one in a K1 record: misrouted
        (&[1, 1, 0][..], 3, "error:routing:false"),
        // -k + p1 + p2 out (the emitted loop sign reversed): residual 2k
        (&[-1, 1, 1][..], 1, "error:routing:false"),
        (&[-1, 1, 1][..], 3, "error:routing:false"),
    ] {
        let got = code(uncontracted::conserve(&delta, 5, &node(out), conserved));
        assert_eq!(got, expected, "edge 4 out {out:?}, {conserved} conserved");
    }
    // no momentum of the node in the Δ: unchanged
    let free = Atom::var(symbol!("u7ad_b_free"));
    assert_eq!(
        code(uncontracted::conserve(&free, 5, &node(&[1, 1, 1]), 3)),
        "Absent"
    );
}

/// The signed mutation (Astra U7a finding 2) at the FULL hook: on the
/// mass-bearing insertion of L3 016's root remainder {2,3} (its `ct:` node:
/// leg 1 out with p, edge 2 in with k, edge 3 in with p - k, so the two loop
/// momenta agree), reverse ONLY the emitted loop coefficient of the edge
/// that carries k out of the node (edge 3: -k -> +k), the contracted source
/// and today's path untouched. The loop residual is then 2k: a misrouted
/// record. Check mode and mode on must FAIL the diagram with
/// `outcome=error:routing`; before the fix the residual made the step
/// return "no rewrite" and the record was the designed
/// `insertion_mass_symbol` (a silent fallback to today's value).
#[test]
fn loop_residual_at_the_full_hook_is_a_routing_error() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let case = l3_016()?;
        let flips = Rc::new(RefCell::new(Vec::<Vec<i64>>::new()));
        let hook: InsertionHook = {
            let flips = flips.clone();
            Rc::new(
                move |kind: Kind, node: usize, delta: &mut Atom, edges: &mut [FactorizedEdge]| {
                    if kind != Kind::Root || !delta.contains_symbol(auxiliary_symbol()) {
                        return;
                    }
                    // the K1 record: one external p(1) after the loop components
                    let loops = edges[0].momentum.len() - 1;
                    let Some(out) = edges.iter().position(|edge| {
                        let sign = incidence(edge, node);
                        !edge.is_external() && edge.momentum[..loops].iter().any(|c| sign * c < 0)
                    }) else {
                        return;
                    };
                    for c in &mut edges[out].momentum[..loops] {
                        *c = -*c;
                    }
                    flips.borrow_mut().push(residual(edges, node));
                },
            )
        };
        for (mode, switch) in [("check", check()), ("on", on())] {
            flips.borrow_mut().clear();
            let _ = ledger::take_lines();
            let (value, captured) =
                with_insertion_hook(hook.clone(), || l3_016_run(&case, &case.delta, switch))?;
            let lines = ledger::take_lines();
            let flipped = flips.borrow().clone();
            eprintln!("U7AD_B L3 016 {mode}: residuals after the flip {flipped:?}");
            ensure!(
                flipped.len() == 1 && flipped[0][0].abs() == 2,
                "L3 016 {mode}: the flip did not give a residual 2k: {flipped:?}"
            );
            let error = match value {
                Ok(_) => {
                    return Err(eyre!(
                        "L3 016 {mode}: the misrouted record did not fail the diagram; outcomes {:?}",
                        captured
                            .iter()
                            .map(|a| (a.kind.as_str(), a.id, a.outcome.clone()))
                            .collect::<Vec<_>>()
                    ));
                }
                Err(error) => format!("{error:#}"),
            };
            ensure!(
                error.contains("outcome=error:routing"),
                "L3 016 {mode}: not a routing error: {error}"
            );
            ensure!(
                captured
                    .iter()
                    .any(|a| a.kind == Kind::Root && a.outcome == "error:routing"),
                "L3 016 {mode}: no root attempt with outcome error:routing"
            );
            ensure!(
                !captured
                    .iter()
                    .any(|a| a.outcome == "declined:insertion_mass_symbol"),
                "L3 016 {mode}: the misrouted record counted as a mass refusal"
            );
            // the merged ledger (lane A): the routing error is the root
            // kind's `error=routing:1 outcome=error`, never a decline, and the
            // diagram's summary is `outcome=error`
            let (root, summary) = disposition_and_summary(&lines, Kind::Root)?;
            eprintln!("U7AD_MERGE L3 016 {mode}: {root}\n{summary}");
            for line in [&root, &summary] {
                ensure!(
                    field(line, "outcome") == Some("error")
                        && field(line, "error") == Some("routing:1")
                        && field(line, "declined") == Some("none"),
                    "L3 016 {mode}: the ledger does not report the routing error: {line}"
                );
            }
        }
        Ok(())
    })
}

/// Control 1 of the mutation: a GENUINE mass dependence on a loop-conserving
/// record (the production Δ times `1 + log(M_aux²)`: the conservation step
/// removes the zero, a `log(M_aux²)` stays): the designed
/// `insertion_mass_symbol`, the diagram passes check and on, on == off.
#[test]
fn genuine_mass_dependence_at_the_full_hook_declines() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let case = l3_016()?;
        let log = function!(symbolica::atom::Symbol::LOG, Atom::var(auxiliary_symbol()));
        let delta = (&case.delta * (Atom::num(1) + log)).expand();
        let (off, nothing) = l3_016_run(&case, &delta, Switch::OFF)?;
        let off = off?;
        ensure!(nothing.is_empty());
        let declined = ["declined:insertion_mass_symbol"];
        let _ = ledger::take_lines();
        let (checked, captured) = l3_016_run(&case, &delta, check())?;
        let lines = ledger::take_lines();
        let checked = checked?;
        ensure!(checked == off, "check returns today's value");
        // the merged ledger (lane A): a named designed decline, the other
        // root records compared, the diagram ok
        let (root, summary) = disposition_and_summary(&lines, Kind::Root)?;
        eprintln!("U7AD_MERGE L3 016 genuine mass check: {root}\n{summary}");
        ensure!(
            field(&root, "outcome") == Some("ok")
                && field(&root, "declined") == Some("insertion_mass_symbol:1")
                && field(&root, "error") == Some("none")
                && field(&root, "compared").is_some_and(|n| n != "0")
                && field(&summary, "outcome") == Some("ok"),
            "L3 016 genuine mass check: {root}\n{summary}"
        );
        evaluated_except("L3 016 genuine mass check", &captured, &declined)?;
        let refusals = captured
            .iter()
            .filter(|a| a.outcome == declined[0])
            .collect::<Vec<_>>();
        ensure!(
            refusals.len() == 1 && refusals[0].kind == Kind::Root,
            "L3 016: the genuine mass dependence is one designed root refusal: {:?}",
            refusals.iter().map(|a| a.id).collect::<Vec<_>>()
        );
        let (value, attempts) = l3_016_run(&case, &delta, on())?;
        let value = value?;
        evaluated_except("L3 016 genuine mass on", &attempts, &declined)?;
        ensure!(atoms_equal(&value, &off), "L3 016: on differs from off");
        Ok(())
    })
}

/// The R* of one proper UV child (`support`, its edge ids) of a diagram as
/// the parent's plan computes an `rstar`-route child, under `switch`
/// (strict, a fresh Δ cache): its Δ and every attempt.
fn child_run(
    diagram: &feynkit_graph::FeynmanDiagram,
    support: &[usize],
    switch: Switch,
) -> Result<(Result<Atom>, Vec<Captured>)> {
    let settings = settings();
    let vakint = Vakint::new()?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
    let child = plan
        .counterterms
        .iter()
        .map(|call| &call.argument.support)
        .find(|child| uv_child_edge_ids(&graph, child) == support)
        .ok_or_else(|| eyre!("{}: no UV child {support:?}", graph.name))?;
    let cache = std::sync::Arc::new(super::CountertermCache::default());
    Ok(with_switch(switch, || {
        with_strict(true, || {
            with_boundary_mode(BoundaryMode::On, false, || {
                with_capture(|| {
                    RStarPlan::rstar_child_counterterm(
                        &cache,
                        &graph,
                        child,
                        model,
                        &vakint,
                        &settings,
                        plan.uv_children,
                        plan.uv_child_closure,
                    )
                })
            })
        })
    }))
}

/// Control 2 of the mutation: INDEPENDENT external crown momenta (a K2
/// record, `externals=independent`: only the loop components must
/// conserve). The open child {5,6,7,8} of L3 36, remainder {5,6} with
/// {7,8} contracted at the `ct:` node that also carries crown leg 3: its Δ
/// gets `log(M_aux²) x dot(Σ σ_e Q(e), Q(leg))`, a zero by momentum
/// conservation at that node. The step removes it (`conserved=1`), the record
/// is evaluated and equal to today's value; on == off. No routing error.
#[test]
fn independent_crown_at_the_full_hook_conserves() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let diagram = colour_factored(&load(L3_36_DOT)?)?;
        let seen = Rc::new(RefCell::new(Vec::<(usize, Vec<i64>, bool)>::new()));
        let hook: InsertionHook = {
            let seen = seen.clone();
            Rc::new(
                move |kind: Kind, node: usize, delta: &mut Atom, edges: &mut [FactorizedEdge]| {
                    if kind != Kind::Child {
                        return;
                    }
                    let legs = edges
                        .iter()
                        .filter(|edge| edge.is_external() && incidence(edge, node) != 0)
                        .map(|edge| edge.id)
                        .collect::<Vec<_>>();
                    let rest = residual(edges, node);
                    let added = !legs.is_empty() && rest.iter().all(|c| *c == 0);
                    seen.borrow_mut().push((node, rest, !legs.is_empty()));
                    if !added {
                        return;
                    }
                    let mut zero = Atom::Zero;
                    for edge in edges.iter() {
                        let sign = incidence(edge, node);
                        if sign != 0 {
                            zero += Atom::num(sign) * function!(GS.emr_mom, edge.id, mink());
                        }
                    }
                    let zero =
                        function!(SPENSO_TAG.dot, zero, function!(GS.emr_mom, legs[0], mink()));
                    let log =
                        function!(symbolica::atom::Symbol::LOG, Atom::var(auxiliary_symbol()));
                    *delta = (&*delta + &*delta * &log * &zero).expand();
                },
            )
        };
        let support = [5, 6, 7, 8];
        let (off, nothing) = child_run(&diagram, &support, Switch::OFF)?;
        let off = off?;
        ensure!(nothing.is_empty());
        for (mode, switch) in [("check", check()), ("on", on())] {
            seen.borrow_mut().clear();
            let (value, captured) =
                with_insertion_hook(hook.clone(), || child_run(&diagram, &support, switch))?;
            eprintln!("U7AD_B L3 36 {mode}: ct: nodes seen {:?}", seen.borrow());
            let value = value?;
            evaluated_except(&format!("L3 36 child {mode}"), &captured, &[])?;
            ensure!(
                seen.borrow()
                    .iter()
                    .any(|(_, rest, leg)| *leg && rest.iter().all(|c| *c == 0)),
                "L3 36 {mode}: no ct: node with a crown leg conserves: {:?}",
                seen.borrow()
            );
            ensure!(
                captured.iter().any(|a| a.kind == Kind::Child
                    && a.conserved == 1
                    && a.outcome == "uncontracted"),
                "L3 36 {mode}: no conserved child record"
            );
            if mode == "check" {
                ensure!(value == off, "check returns today's Δ");
            } else {
                ensure!(atoms_equal(&value, &off), "L3 36: on differs from off");
            }
        }
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// U7a.3: the scopes restore on unwinding
// ---------------------------------------------------------------------------

const L2_22_DOT: &str = include_str!("testdata/u7a/l2-d22.dot");

/// The run of a diagram root (the beta's full-graph closing call) under
/// `switch`, strict, remainder boundary on: its value and every attempt.
fn root_run(
    diagram: &feynkit_graph::FeynmanDiagram,
    switch: Switch,
) -> Result<(Result<Atom>, Vec<Captured>)> {
    let settings = settings();
    let vakint = Vakint::new()?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
    let closed_root = metric_root(&graph);
    let callback = |tensor: &Atom| -> Result<Atom> {
        Ok(tensor.replace_multiple(&closed_root.external_routing) * &closed_root.factor)
    };
    Ok(with_switch(switch, || {
        with_strict(true, || {
            with_boundary_mode(BoundaryMode::On, false, || {
                with_capture(|| {
                    plan.projected_counterterm_closing_root_scalar(
                        &graph,
                        model,
                        &vakint,
                        &settings,
                        &callback,
                        &closed_root,
                    )
                })
            })
        })
    }))
}

/// Astra U7a finding 3: a panic inside an insertion surrogate's R* (the
/// dressed scope), caught on the same thread (beta catches a diagram's panic
/// and its worker goes on), must not leave the scope set: the next ordinary
/// root and child are classified root and child (not the inert surrogate),
/// and an ordinary diagram on that thread then emits and compares its root
/// and child records. The quotient scope likewise.
#[test]
fn a_caught_unwind_leaves_no_scope_behind() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let caught = std::panic::catch_unwind(|| {
            super::expression::in_dressed_scope(|| -> usize {
                panic!("U7a.3 control: a panic inside the dressed scope")
            })
        });
        ensure!(caught.is_err(), "the dressed scope's panic is caught");
        ensure!(
            uncontracted::classify(false) == Kind::Root
                && uncontracted::classify(true) == Kind::Child
                && uncontracted::classify_ir_hard() == Kind::IrHard,
            "after the unwind: root {:?}, child {:?}, ir-hard {:?}",
            uncontracted::classify(false),
            uncontracted::classify(true),
            uncontracted::classify_ir_hard()
        );
        let caught = with_switch(check(), || {
            std::panic::catch_unwind(|| {
                uncontracted::scoped(uncontracted::Scope::Quotient, || -> usize {
                    panic!("U7a.3 control: a panic inside the quotient scope")
                })
            })
        });
        ensure!(caught.is_err(), "the quotient scope's panic is caught");
        ensure!(
            uncontracted::classify(false) == Kind::Root
                && uncontracted::classify(true) == Kind::Child,
            "after the quotient unwind: root {:?}, child {:?}",
            uncontracted::classify(false),
            uncontracted::classify(true)
        );
        // executed as their kinds: an ordinary diagram on the same thread
        let diagram = colour_factored(&load(L2_22_DOT)?)?;
        let (value, captured) = root_run(&diagram, check())?;
        let _ = value?;
        evaluated_except("L2 22 after the unwinds", &captured, &[])?;
        for kind in ROOT_CHILD {
            ensure!(
                captured.iter().any(|a| a.kind == kind
                    && a.outcome == "uncontracted"
                    && a.equal == Some(true)),
                "L2 22 after the unwinds: no {} record emitted and compared",
                kind.as_str()
            );
        }
        ensure!(
            captured
                .iter()
                .all(|a| matches!(a.kind, Kind::Root | Kind::Child)),
            "L2 22 after the unwinds: a record of another kind"
        );
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// one_loop_tree_cut: the single-propagator shape only
// ---------------------------------------------------------------------------

const L2_25_DOT: &str = include_str!("testdata/u7a/l2-d25.dot");

/// `keep the rebase test narrow` (Astra U7a; U7-PLAN 14(d)): the emitter
/// declines with `one_loop_tree_cut` ONLY a remainder whose one internal
/// propagator is the closing line (after the Taylor step a massive tadpole,
/// which vakint's closure entry refuses as "every line is massive"; L2 25 has
/// one, as L4 49 id 8 has). Every other one-loop record (today's route names
/// every one-loop remainder `one_loop_tree_cut`: it integrates it whole;
/// 137 roots and 410 children at L3) stays on the new path: evaluated,
/// compared equal, with at least two internal propagators; L2 25 and L2 22
/// carry such root and child neighbours.
#[test]
fn one_loop_tree_cut_is_the_single_propagator_shape_only() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let mut neighbours = BTreeMap::<Kind, usize>::new();
        for (what, dot, shapes) in [("L2/25", L2_25_DOT, 1), ("L2/22", L2_22_DOT, 0)] {
            let diagram = colour_factored(&load(dot)?)?;
            let _ = ledger::take_lines();
            let (value, captured) = root_run(&diagram, check())?;
            let lines = ledger::take_lines();
            let _ = value?;
            evaluated_except(what, &captured, &["declined:one_loop_tree_cut"])?;
            // the merged ledger (lane A): the shape is a NAMED decline of the
            // root kind's disposition (no error, the neighbours compared, the
            // diagram passes check without an exception)
            let (root, summary) = disposition_and_summary(&lines, Kind::Root)?;
            eprintln!("U7AD_MERGE {what} check: {root}\n{summary}");
            let named = if shapes == 0 {
                "none".to_owned()
            } else {
                format!("one_loop_tree_cut:{shapes}")
            };
            ensure!(
                field(&root, "outcome") == Some("ok")
                    && field(&root, "declined") == Some(named.as_str())
                    && field(&root, "error") == Some("none")
                    && field(&root, "excepted") == Some("none")
                    && field(&root, "compared").is_some_and(|n| n != "0")
                    && field(&summary, "outcome") == Some("ok"),
                "{what} check: {root}\n{summary}"
            );
            let declined = captured
                .iter()
                .filter(|a| a.outcome.starts_with("declined:") || a.outcome.starts_with("error:"))
                .collect::<Vec<_>>();
            ensure!(
                declined.len() == shapes
                    && declined
                        .iter()
                        .all(|a| a.outcome == "declined:one_loop_tree_cut" && a.loops == 1),
                "{what}: the declines {:?} are not {shapes} single-propagator shape(s)",
                declined
                    .iter()
                    .map(|a| (a.kind.as_str(), a.id, a.loops, a.outcome.clone()))
                    .collect::<Vec<_>>()
            );
            for attempt in captured
                .iter()
                .filter(|a| a.loops == 1 && a.outcome == "uncontracted")
            {
                let text = attempt
                    .text
                    .as_ref()
                    .ok_or_else(|| eyre!("{what}: an evaluated record without text"))?;
                let record = FactorizedNumerator::from_text(text).map_err(|e| eyre!("{e}"))?;
                let internal = record.edges.iter().filter(|e| !e.is_external()).count();
                ensure!(
                    internal >= 2,
                    "{what}: an evaluated one-loop record with {internal} propagator(s)\n{text}"
                );
                *neighbours.entry(attempt.kind).or_default() += 1;
            }
        }
        eprintln!("U7AD_B one-loop neighbours evaluated: {neighbours:?}");
        ensure!(
            neighbours.get(&Kind::Root).is_some_and(|n| *n > 0)
                && neighbours.get(&Kind::Child).is_some_and(|n| *n > 0),
            "one-loop root and child neighbours on the new path: {neighbours:?}"
        );
        Ok(())
    })
}

/// The shape itself, on record edges: the closing line alone (a self-loop
/// at a `ct:` node, the legs external) is the single-propagator shape; a
/// one-loop remainder with a second propagator is not, nor one whose only
/// internal edge is not the closing line.
#[test]
fn single_propagator_is_the_closing_line_alone() {
    let edge = |id: usize, source: Endpoint, target: Endpoint| FactorizedEdge {
        id,
        source,
        target,
        source_hedge: None,
        target_hedge: None,
        particle: "g".into(),
        line: vakint::p_integral::factorized::LineKind::Vector,
        power: 1,
        mass_squared: Atom::Zero,
        momentum: vec![1, 0],
        factor: Atom::num(1),
    };
    let (ct, v) = (Endpoint::Vertex(3), Endpoint::Vertex(0));
    let legs = [
        edge(0, Endpoint::External, ct),
        edge(1, ct, Endpoint::External),
    ];
    let tadpole = [legs[0].clone(), legs[1].clone(), edge(8, ct, ct)];
    assert!(uncontracted::single_propagator(&tadpole, 8));
    assert!(!uncontracted::single_propagator(&tadpole, 7));
    let bubble = [
        legs[0].clone(),
        legs[1].clone(),
        edge(7, ct, v),
        edge(8, v, ct),
    ];
    assert!(!uncontracted::single_propagator(&bubble, 8));
    assert!(!uncontracted::single_propagator(&legs, 8));
}

// ---------------------------------------------------------------------------
// U7-PLAN requirements 13 and 14: the frozen REFUSALS list
// ---------------------------------------------------------------------------

/// U7-PLAN requirements 13 and 14: a vakint `Closure` refusal is designed
/// only by a row of `DESIGNED_VAKINT` (`Taylor degree `, `the projector
/// contains Q(`, the open-spinor port `open slot ` + `is not a Lorentz slot`,
/// and the backstop `every line is massive: not a closure`, which the
/// emitter declines first as `one_loop_tree_cut`); any other `Closure`
/// detail is an error, among them the entry's wrapped internal failures.
#[test]
fn closure_refusals_are_designed_only_by_a_frozen_row() {
    let closure = |detail: &str| FactorizedError::Closure {
        detail: detail.into(),
    };
    for (detail, code) in [
        ("Taylor degree 3 above 2", "vakint:taylor_degree"),
        (
            "the projector contains Q(0, ..): with a Taylor step it acts after the step and must be momentum-free",
            "vakint:projector_momentum",
        ),
        (
            "open slot bis(4,hedge(3,1)) is not a Lorentz slot (open bispinor ports are not supported; close the line with the projector)",
            "vakint:open_spinor",
        ),
        ("every line is massive: not a closure", "vakint:all_massive"),
    ] {
        let decline = Decline::vakint(&closure(detail));
        assert!(decline.designed, "{detail}: {decline:?}");
        assert_eq!(decline.code, code);
    }
    for detail in [
        "lowering: no key plan",
        "the lines do not form one closing line among massless ones",
        "truncating the insertion product: x",
        "the projector contains the momentum of internal edge 3",
        "open slot bis(4,hedge(3,1)): only Lorentz slots stay open",
    ] {
        let decline = Decline::vakint(&closure(detail));
        assert!(!decline.designed, "{detail}: {decline:?}");
        assert_eq!(decline.code, "vakint:Closure", "{detail}");
    }
    let prefixes = [
        "Taylor degree ",
        "the projector contains Q(",
        "open slot ",
        "every line is massive: not a closure",
    ];
    for row in uncontracted::DESIGNED_VAKINT
        .iter()
        .filter(|row| row.0 == "Closure")
    {
        assert!(
            prefixes.contains(&row.1),
            "the Closure row {:?} is not one of the four",
            row.1
        );
    }
}

/// Our designed vakint rows are a SUBSET of vakint's refusal table
/// (`crates/vakint/tests/fixtures/factorized_refusals.tsv`, read here; vakint's
/// `the_refusal_table_matches_the_source` keeps it equal to the source): each
/// row with a detail prefix lies in the literal start of a `refusal` row of
/// its variant (the template's text before its first `<>`), and the detail written
/// from that template is declined with the row's code; the rows matched by
/// name (Gamma5, LeviCivita: struct variants without a free-text detail) are
/// the only ones without a table row; no row is `Evaluation` or `Oracle`, and
/// those two are always errors.
#[test]
fn designed_vakint_rows_are_a_subset_of_the_frozen_refusals() {
    let table = include_str!("../../../../vakint/tests/fixtures/factorized_refusals.tsv");
    let refusals: Vec<(&str, &str)> = table
        .lines()
        .filter_map(|line| {
            let mut fields = line.split('\t');
            match (fields.next(), fields.next(), fields.next()) {
                (Some("refusal"), Some(variant), Some(template)) => Some((variant, template)),
                _ => None,
            }
        })
        .collect();
    assert!(refusals.len() > 100, "{} refusal rows", refusals.len());
    // does `prefix` lie in the literal start of `template` (its text before
    // the first `<>`)? A template that begins with `<>` starts no prefix.
    let starts = |template: &str, prefix: &str| {
        template
            .split("<>")
            .next()
            .unwrap_or("")
            .starts_with(prefix)
    };
    let mut by_name = Vec::new();
    for (variant, prefix, phrase, code) in uncontracted::DESIGNED_VAKINT {
        assert!(
            *variant != "Evaluation" && *variant != "Oracle",
            "{code}: {variant} is never designed"
        );
        if prefix.is_empty() {
            assert!(
                refusals.iter().all(|(v, _)| v != variant),
                "{code}: {variant} has table rows, match it by its detail"
            );
            by_name.push(*code);
            continue;
        }
        let rows: Vec<&str> = refusals
            .iter()
            .filter(|(v, template)| {
                v == variant && starts(template, prefix) && template.contains(phrase)
            })
            .map(|(_, template)| *template)
            .collect();
        assert!(
            !rows.is_empty(),
            "{code}: ({variant}, {prefix:?}) is in no table row"
        );
        for template in rows {
            let detail = template.replace("<>", "x");
            let error = match *variant {
                "Closure" => FactorizedError::Closure { detail },
                "Ports" => FactorizedError::Ports { detail },
                "Insertion" => FactorizedError::Insertion { vertex: 1, detail },
                other => panic!("{code}: no detail constructor for {other}"),
            };
            let decline = Decline::vakint(&error);
            assert!(decline.designed, "{code}: {template} -> {decline:?}");
            assert_eq!(decline.code, format!("vakint:{code}"), "{template}");
        }
    }
    assert_eq!(by_name, ["gamma5", "levi_civita"]);
    assert_eq!(uncontracted::DESIGNED_VAKINT.len(), 9);
    for (error, code) in [
        (
            FactorizedError::Evaluation {
                detail: "Taylor degree 3 above 2".into(),
            },
            "vakint:Evaluation",
        ),
        (
            FactorizedError::Oracle {
                detail: "every line is massive: not a closure".into(),
            },
            "vakint:Oracle",
        ),
    ] {
        let decline = Decline::vakint(&error);
        assert!(!decline.designed, "{error}: {decline:?}");
        assert_eq!(decline.code, code, "{error}");
    }
}

// ---------------------------------------------------------------------------
// R2-I3: the record builder's impossible states are errors
// ---------------------------------------------------------------------------

/// A K1 attempt on the bare root remainder of `dot` with the closed root
/// `closed_root` (default: the beta's metric root).
fn root_site_attempt(
    dot: &str,
    closed_root: impl FnOnce(&Graph) -> ClosedRootNormalization,
) -> Result<uncontracted::Attempt> {
    let diagram = colour_factored(&load(dot)?)?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let vakint = Vakint::new()?;
    let uv = settings();
    let vakint_settings = uv.vakint.true_settings();
    let plan = RStarPlan::with_automatic_edge(&graph, model, &uv)?;
    let term = plan
        .ordinary_uv
        .iter()
        .find(|term| term.remainder.insertions.is_empty() && term.remainder.loops > 0)
        .ok_or_else(|| eyre!("{}: no bare root remainder", graph.name))?;
    let closure = RStarIntegrand::root_closure(&graph, &plan.support)?;
    let closed_root = closed_root(&graph);
    let auxiliary = auxiliary_mass();
    let children = BTreeMap::new();
    let site = Site {
        graph: &graph,
        model,
        plan: &plan,
        term,
        children: &children,
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
        ir_hard: None,
    };
    uncontracted::attempt(&site, Kind::Root, 1, Mode::Check, || {
        RStarIntegrand::new_with_closure(
            &graph,
            model,
            &term.remainder,
            &[],
            plan.marked_edge,
            Atom::var(*REFERENCE_MASS_SQUARED),
            Some(&closure),
        )
    })
}

fn error_code(attempt: &uncontracted::Attempt) -> Option<&str> {
    match &attempt.outcome {
        Outcome::Error(decline) => Some(decline.code.as_str()),
        _ => None,
    }
}

/// R2-I3: `sigma`, `root_crown_paired` and `root_normalisation_momentum`
/// are states of the record builder that never occur on a correct site: a
/// regression that produced them must be an error, not a designed fallback.
/// Driven: a σ that does not fold the second root leg onto ±p(1) (L2 22's
/// bare root remainder, the external routing dropped) and a closed-root
/// factor that is not a scalar times a power of p².
#[test]
fn record_builder_states_are_errors() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        for code in ["sigma", "root_crown_paired", "root_normalisation_momentum"] {
            ensure!(
                !DESIGNED_LOCAL.contains(&code),
                "{code} is a designed refusal"
            );
        }
        let attempt = root_site_attempt(L2_22_DOT, |graph| {
            let mut closed_root = metric_root(graph);
            closed_root.external_routing.clear();
            closed_root
        })?;
        ensure!(
            error_code(&attempt) == Some("sigma"),
            "a σ that does not fold the legs: outcome {:?}",
            attempt.outcome
        );
        let q = |e: usize| function!(GS.emr_mom, e, mink());
        let factor = function!(SPENSO_TAG.dot, q(0), q(1)) * Atom::var(symbol!("u7ad_b_c"));
        let decline = split_root_factor(&factor, 0)
            .err()
            .ok_or_else(|| eyre!("the factor {factor} split"))?;
        ensure!(
            decline.code == "root_normalisation_momentum" && !decline.designed,
            "a closed-root factor with a momentum: {decline:?}"
        );
        // the control: the beta's own metric root is evaluated
        let attempt = root_site_attempt(L2_22_DOT, metric_root)?;
        ensure!(
            matches!(attempt.outcome, Outcome::Value(_)),
            "the metric root: outcome {:?}",
            attempt.outcome
        );
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// R2-I4: the open-spinor row needs its Lorentz phrase
// ---------------------------------------------------------------------------

/// R2-I4: vakint's `open slot ` row is the closure entry's open-spinor
/// refusal ("... is not a Lorentz slot (open bispinor ports are not
/// supported ...)"), not every detail with that prefix (e.g. the vacuum
/// entry's "open slot <port>: only Lorentz slots ... stay open").
#[test]
fn open_slot_needs_the_lorentz_phrase() {
    let closure = |detail: &str| FactorizedError::Closure {
        detail: detail.into(),
    };
    let designed = Decline::vakint(&closure(
        "open slot bis(4,hedge(3,1)) is not a Lorentz slot (open bispinor ports are not supported; close the line with the projector)",
    ));
    assert!(designed.designed, "{designed:?}");
    assert_eq!(designed.code, "vakint:open_spinor");
    for detail in [
        "open slot bis(4,hedge(3,1)): only Lorentz slots mink(D, index) stay open in the vacuum scope",
        "open slot bis(4,hedge(3,1))",
    ] {
        let decline = Decline::vakint(&closure(detail));
        assert!(!decline.designed, "{detail}: {decline:?}");
        assert_eq!(decline.code, "vakint:Closure", "{detail}");
    }
}

// ---------------------------------------------------------------------------
// helper VH round 2 minors: the crown code (R2-M5), the second leg id (R2-M4)
// ---------------------------------------------------------------------------

/// The K2 attempt on the remainder {5,6} (with {7,8} contracted) of L3 36's
/// open child {5,6,7,8}, colour factored, with the crown `root_edges`
/// (`None`: the child's own) and the nested Δ `delta`.
fn child_site_attempt(
    root_edges: Option<Vec<EdgeIndex>>,
    delta: Atom,
) -> Result<uncontracted::Attempt> {
    let diagram = colour_factored(&load(L3_36_DOT)?)?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let vakint = Vakint::new()?;
    let uv = settings();
    let vakint_settings = uv.vakint.true_settings();
    let root_plan = RStarPlan::with_automatic_edge(&graph, model, &uv)?;
    let support = root_plan
        .counterterms
        .iter()
        .map(|call| call.argument.support.clone())
        .find(|support| uv_child_edge_ids(&graph, support) == [5, 6, 7, 8])
        .ok_or_else(|| eyre!("L3 36 has no child {{5,6,7,8}}"))?;
    let plan = RStarPlan::for_support(&graph, model, &support, &uv)?;
    let term = plan
        .ordinary_uv
        .iter()
        .find(|term| {
            join_edge_ids(&uv_child_edge_ids(&graph, &term.remainder.remaining)) == "5,6"
                && term.remainder.insertions.len() == 1
        })
        .ok_or_else(|| eyre!("the child plan has no term {{5,6}} with {{7,8}} contracted"))?;
    let nested = term.remainder.insertions[0].support.clone();
    let children = BTreeMap::from([(nested, delta.clone())]);
    let boundary = plan
        .root_boundary
        .clone()
        .ok_or_else(|| eyre!("a child plan has a root boundary"))?;
    let crown = graph
        .iter_edges_of(&boundary)
        .map(|(_, edge, _)| edge)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let root_edges = root_edges.unwrap_or(crown);
    let auxiliary = auxiliary_mass();
    let site = Site {
        graph: &graph,
        model,
        plan: &plan,
        term,
        children: &children,
        root_closure: None,
        closed_root: None,
        closure_data: None,
        root_projector: false,
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
    let insertions = vec![delta];
    uncontracted::attempt(&site, Kind::Child, 1, Mode::Check, || {
        RStarIntegrand::new_with_closure(
            &graph,
            model,
            &term.remainder,
            &insertions,
            plan.marked_edge,
            Atom::var(*REFERENCE_MASS_SQUARED),
            None,
        )
    })
}

/// R2-M5: a K2 crown that is not the remainder's soft edges (here: one crown
/// edge dropped, then none) is its own error code `crown`, not the K1 code
/// `routing`.
#[test]
fn crown_mismatch_is_its_own_error() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        // the nested Δ of {7,8} as wh-fixtures/k2-child-ins states it
        let fixture =
            FactorizedNumerator::from_text(include_str!("testdata/u7a/k2-child-ins.input.txt"))
                .map_err(|e| eyre!("{e}"))?;
        let insertion = fixture
            .insertions
            .first()
            .ok_or_else(|| eyre!("k2-child-ins has no insertion"))?;
        let delta = &insertion.coefficient * &insertion.factor;
        let full = child_site_attempt(None, delta.clone())?;
        ensure!(
            matches!(full.outcome, Outcome::Value(_)),
            "the child's own crown: outcome {:?}",
            full.outcome
        );
        let diagram = colour_factored(&load(L3_36_DOT)?)?;
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let crown = graph
            .iter_edges()
            .map(|(_, edge, _)| edge)
            .filter(|edge| [3usize, 4].contains(&usize::from(*edge)))
            .collect::<Vec<_>>();
        ensure!(crown.len() == 2, "L3 36's crown edges 3 and 4");
        for root_edges in [vec![crown[0]], Vec::new()] {
            let attempt = child_site_attempt(Some(root_edges.clone()), delta.clone())?;
            ensure!(
                error_code(&attempt) == Some("crown"),
                "crown {root_edges:?}: outcome {:?}",
                attempt.outcome
            );
        }
        Ok(())
    })
}

/// R2-M4: the incoming leg of a crown edge with both ends in the scope has
/// the record id `PAIRED_LEG_OFFSET + id` (the VERIFIED k2-child states legs
/// 4 and 100004), valid only while every edge id of the graph stays below
/// the offset; a graph that reaches it is the `crown` error, not a silent
/// collision.
#[test]
fn paired_crown_leg_ids_stay_below_the_offset() {
    assert_eq!(uncontracted::PAIRED_LEG_OFFSET, 100_000);
    assert_eq!(uncontracted::paired_leg_id(4, 26), Ok(100_004));
    assert_eq!(uncontracted::paired_leg_id(4, 99_999), Ok(100_004));
    let decline = uncontracted::paired_leg_id(4, 100_000).unwrap_err();
    assert_eq!(decline.code, "crown");
    assert!(!decline.designed, "{decline:?}");
}

// ---------------------------------------------------------------------------
// Requirement 15 on K1 / K2: the real ε^0-only records
// ---------------------------------------------------------------------------

/// U7-PLAN requirement 15 for K1 / K2: the needed order of a K1 root or K2
/// child value is its poles (`valid_order(depth = loops, loops) = -1`: the
/// consumer adds it with coefficient one, then -K); a first difference
/// above it passes, at or below it fails. On two of U7a's 158 ε^0-only K2
/// records as the gated binary dumped them in check mode
/// (results/consolidation/U7a-EPS0-EVIDENCE/records, byte-identical copies:
/// L3 016 id 2 sha256 0b0ab27e..., L4 0042 id 4 sha256 95503f38...): the
/// poles are identical, today's value carries an ε^0 that vakint's (stopping
/// at ε^-1) does not; the check records `first_difference=eps^0` and passes.
/// A pole changed by `c/ε` fails at `eps^-1`.
#[test]
fn k1_k2_depth_check_on_the_eps0_only_records() -> Result<()> {
    test_initialise()?;
    let epsilon = GS.dim_epsilon;
    for (what, text) in [
        (
            "L3 016 id 2",
            include_str!("testdata/u7a/eps0-L3-016-id2.value.txt"),
        ),
        (
            "L4 0042 id 4",
            include_str!("testdata/u7a/eps0-L4-0042-id4.value.txt"),
        ),
    ] {
        let field = |key: &str| {
            text.lines()
                .find_map(|line| line.strip_prefix(key)?.strip_prefix('\t'))
                .ok_or_else(|| eyre!("{what}: no {key} line"))
        };
        ensure!(field("kind")? == "child" && field("equal")? == "true");
        let depth = field("depth")?.parse::<usize>()?;
        // K1 / K2: depth = the remainder's loops, so the order is the poles
        let order = uncontracted::valid_order(depth, depth);
        ensure!(order == -1, "{what}: the needed order eps^{order}");
        let today = parse(field("contracted")?)?;
        let ours = parse(field("uncontracted")?)?;
        let check = uncontracted::compare(&today, &ours, epsilon, order, 0.0, 0.0);
        eprintln!(
            "U7AD_B {what}: identical={} equal={} full={} first_difference={}",
            check.identical,
            check.equal,
            check.full,
            check.first_difference.as_text()
        );
        ensure!(
            !check.identical && check.equal && !check.full,
            "{what}: identical={} equal={} full={} error={:?}",
            check.identical,
            check.equal,
            check.full,
            check.error
        );
        ensure!(
            check.first_difference == uncontracted::FirstDifference::At(0),
            "{what}: first difference {}",
            check.first_difference.as_text()
        );
        // a pole changed: fails, at eps^-1
        let wrong = &ours + Atom::var(symbol!("u7ad_b_pole")) / Atom::var(epsilon);
        let failed = uncontracted::compare(&today, &wrong, epsilon, order, 0.0, 0.0);
        ensure!(
            !failed.equal && failed.first_difference == uncontracted::FirstDifference::At(-1),
            "{what}: a changed pole: equal={} first difference {}",
            failed.equal,
            failed.first_difference.as_text()
        );
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// U7a alone finding: a paired crown edge's incoming leg on a `ct:` node
// ---------------------------------------------------------------------------

const L5_101213_DOT: &str = include_str!("testdata/u7a/l5-d045537.dot");

/// U7a-ALONE-FINDING: L5 101213 (regrouped 45537) run alone, the 4-loop
/// child {2,3,4,5,6,8,9,10,11,12,13}. Its crown edge 7 has both ends in the
/// child (legs 7 out of vertex 1, 100007 into vertex 7); remainder record
/// id 2 contracts {12,13} onto the `ct:` vertex 6, which carries leg 100007
/// (half-edge 13). The Δ of {12,13} that this graph fills itself is written
/// with `Q(7, ..)`, the GRAPH's momentum of that boundary edge, which is not
/// a momentum of a leg on vertex 6: vakint refused the record
/// (`error:vakint:UnsupportedTensor`, "not a boundary momentum of the
/// insertion") and the diagram failed in check and on. In a set the Δ came
/// from 059943's fill, written in `Q(9)` and `Q(11)`, and compared. The
/// emitter writes the incoming leg's own momentum `Q(100007)` (equal to
/// `Q(7)`): every record is evaluated, check is today's value, on agrees.
#[test]
fn paired_crown_incoming_leg_on_a_ct_node() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let diagram = colour_factored(&load(L5_101213_DOT)?)?;
        let support = [2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 13];
        let (off, nothing) = child_run(&diagram, &support, Switch::OFF)?;
        let off = off?;
        ensure!(nothing.is_empty());
        for (mode, switch) in [("check", check()), ("on", on())] {
            let (value, captured) = child_run(&diagram, &support, switch)?;
            let value = value.map_err(|error| eyre!("L5 101213 child {mode}: {error:#}"))?;
            evaluated_except(&format!("L5 101213 child {mode}"), &captured, &[])?;
            let crown = captured
                .iter()
                .filter(|a| a.kind == Kind::Child && a.outcome == "uncontracted")
                .filter_map(|a| a.text.as_deref())
                .find(|text| {
                    text.contains("contracted=12,13 ")
                        && text.contains("\nedge\t100007\tsource=ext\ttarget=6\t")
                        && text.contains("\ninsertion\t6\thedges=13,17,21\t")
                });
            let Some(text) = crown else {
                return Err(eyre!(
                    "L5 101213 child {mode}: no evaluated record with {{12,13}} contracted onto vertex 6 next to leg 100007"
                ));
            };
            let factor = text
                .lines()
                .find(|line| line.starts_with("insertion_factor\t6\t"))
                .unwrap_or_default();
            eprintln!("U7A_ALONE L5 101213 {mode}: {factor}");
            ensure!(
                !factor.contains("::Q(7,"),
                "L5 101213 child {mode}: the insertion at vertex 6 keeps Q(7): {factor}"
            );
            if mode == "check" {
                ensure!(value == off, "check returns today's Δ");
            } else {
                ensure!(atoms_equal(&value, &off), "L5 101213: on differs from off");
            }
        }
        Ok(())
    })
}

/// The rename itself: only an incoming leg on the node renames, an edge
/// with both legs (or only its outgoing one) on the node keeps its id.
#[test]
fn crown_legs_at_node_renames_only_an_incoming_leg() {
    test_initialise().unwrap();
    let q = |e: usize| function!(GS.emr_mom, e, mink());
    let crowns = [
        uncontracted::PairedCrown {
            edge: 7,
            source_hedge: 12,
            incoming: 100_007,
            sink_hedge: 13,
        },
        uncontracted::PairedCrown {
            edge: 3,
            source_hedge: 4,
            incoming: 100_003,
            sink_hedge: 5,
        },
    ];
    let delta = q(7) + q(3) + q(9);
    // node 6 carries 13 (7's incoming leg) and 17; node 1 carries 12 (7's
    // outgoing leg); node 2 carries both legs of 3
    let ends = [(6, 13), (6, 17), (1, 12), (2, 4), (2, 5)];
    assert_eq!(
        uncontracted::crown_legs_at_node(&delta, 6, &ends, &crowns),
        q(100_007) + q(3) + q(9)
    );
    assert_eq!(
        uncontracted::crown_legs_at_node(&delta, 1, &ends, &crowns),
        delta
    );
    assert_eq!(
        uncontracted::crown_legs_at_node(&delta, 2, &ends, &crowns),
        delta
    );
}

/// IRH (iii) on the gated G-D: a field outside the insertions that is not
/// provably regular at ε = 0 is the designed decline `factor_eps_pole`
/// before the record is written. The pure test, every spelling: the
/// production fields (a vertex polynomial in tensors with `mink(d, ..)`, the
/// projector normalisations 1/(d-1), 1/d, 1/(d+2), a positive power of ε)
/// pass; ε^-1, a pole hidden in a sum, 1/(d-4), -2/(d-4), a tensor under a
/// negative power (Astra's tensor-hidden spelling) and Γ(ε) (an ε inside a
/// function argument) are refused.
#[test]
fn factor_eps_pole_reads_every_spelling() {
    let _ = test_initialise();
    let eps = Atom::var(GS.dim_epsilon);
    let d = Atom::var(GS.dim);
    let mu = Minkowski {}.new_rep(GS.dim).to_symbolic([Atom::var(symbol!("mu"))]);
    let nu = Minkowski {}.new_rep(GS.dim).to_symbolic([Atom::var(symbol!("nu"))]);
    let g = function!(symbol!("g_test"), mu.clone(), nu.clone());
    let q = function!(GS.emr_mom, 5, mu.clone());
    let vertex = (&q * Atom::var(symbol!("gs"))) + &g;
    let scalars = (&d - Atom::num(1)).pow(Atom::num(-1))
        * d.clone().pow(Atom::num(-1))
        * (&d + Atom::num(2)).pow(Atom::num(-1));
    for regular in [
        Atom::num(1),
        vertex.clone(),
        scalars.clone(),
        &scalars * &g,
        &vertex * eps.clone().pow(Atom::num(2)),
        &vertex * (Atom::num(1) + &eps),
    ] {
        assert_eq!(uncontracted::factor_eps_pole(&regular), None, "{regular}");
    }
    let pole = |atom: Atom, why: &str| {
        let found = uncontracted::factor_eps_pole(&atom);
        assert!(
            found.as_deref().is_some_and(|found| found.starts_with(why)),
            "{atom}: {found:?}"
        );
    };
    pole(&vertex * eps.clone().pow(Atom::num(-1)), "a pole of order 1");
    pole(&vertex + eps.clone().pow(Atom::num(-1)), "a pole of order 1");
    pole(&vertex * eps.clone().pow(Atom::num(-3)), "a pole of order 3");
    pole((&d - Atom::num(4)).pow(Atom::num(-1)), "a pole of order 1");
    pole(
        Atom::num(-2) * (&d - Atom::num(4)).pow(Atom::num(-1)) * &vertex,
        "a pole of order 1",
    );
    pole(
        Atom::num(-2) * (&g * &g - Atom::num(4)).pow(Atom::num(-1)),
        "a tensor under a negative power",
    );
    pole(
        function!(symbol!("Gamma"), eps.clone()) * &vertex,
        "an order in ε that is not provably",
    );
}

/// IRH (iii) end to end on a real K2 site (L3 36, child {5,6,7,8}, {7,8}
/// contracted): (1) the record as built has no such field; (2) the same
/// record with ε^-1 in a VERTEX factor (the `hidden` witness's shape) is
/// refused by vakint's gated validation with `NotScalar` "a pole in ε
/// outside an insertion", which is NOT designed (Decline::vakint gives an
/// error), and our pre-write check names it; (3) the builder declines such a
/// record BEFORE it is written: with ε^-1 on one edge factor and ε on another
/// (the product check of check mode still holds) the attempt is the
/// designed `declined:factor_eps_pole`, no record text, not an error.
#[test]
fn a_pole_outside_the_insertions_is_declined_before_the_record() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let fixture =
            FactorizedNumerator::from_text(include_str!("testdata/u7a/k2-child-ins.input.txt"))
                .map_err(|e| eyre!("{e}"))?;
        let insertion = fixture
            .insertions
            .first()
            .ok_or_else(|| eyre!("k2-child-ins has no insertion"))?;
        let delta = &insertion.coefficient * &insertion.factor;
        let eps = Atom::var(GS.dim_epsilon);
        // (1) the record as built
        let full = child_site_attempt(None, delta.clone())?;
        ensure!(
            matches!(full.outcome, Outcome::Value(_)),
            "unmodified: outcome {:?}",
            full.outcome
        );
        let text = full
            .text
            .clone()
            .ok_or_else(|| eyre!("a built record has its text"))?;
        let record = FactorizedNumerator::from_text(&text).map_err(|e| eyre!("{e}"))?;
        ensure!(
            uncontracted::record_eps_pole(&record).is_none(),
            "the built record: {:?}",
            uncontracted::record_eps_pole(&record)
        );
        // (2) ε^-1 in a vertex factor: vakint refuses it, not designed
        let mut hidden = record.clone();
        let vertex = hidden
            .vertices
            .iter_mut()
            .find(|vertex| !vertex.rule.starts_with("ct:"))
            .ok_or_else(|| eyre!("the record has a model vertex"))?;
        let id = vertex.id;
        vertex.factor = &vertex.factor * eps.clone().pow(Atom::num(-1));
        let ours = uncontracted::record_eps_pole(&hidden);
        ensure!(
            ours.as_deref()
                .is_some_and(|ours| ours.starts_with(&format!("vertex {id} factor: a pole of order 1"))),
            "ours: {ours:?}"
        );
        let refusal = hidden
            .validate()
            .err()
            .ok_or_else(|| eyre!("vakint accepted a pole in a vertex factor"))?;
        ensure!(
            matches!(&refusal, FactorizedError::NotScalar { detail, .. }
                if detail.starts_with("a pole in ε outside an insertion")),
            "vakint: {refusal:?}"
        );
        let decline = Decline::vakint(&refusal);
        ensure!(!decline.designed, "NotScalar must stay an error: {decline:?}");
        // (3) the builder declines it before writing
        let hook: InsertionHook = Rc::new(
            move |kind: Kind, _node: usize, _delta: &mut Atom, edges: &mut [FactorizedEdge]| {
                if kind != Kind::Child || edges.len() < 2 {
                    return;
                }
                let eps = Atom::var(GS.dim_epsilon);
                edges[0].factor = &edges[0].factor * eps.clone().pow(Atom::num(-1));
                edges[1].factor = &edges[1].factor * eps;
            },
        );
        let attempt = with_insertion_hook(hook, || child_site_attempt(None, delta.clone()))?;
        match &attempt.outcome {
            Outcome::Declined(decline) => {
                ensure!(decline.code == "factor_eps_pole", "{decline:?}");
                ensure!(
                    decline.detail.contains("factor: a pole of order 1"),
                    "{decline:?}"
                );
            }
            other => return Err(eyre!("expected declined:factor_eps_pole, got {other:?}")),
        }
        ensure!(attempt.text.is_none(), "a declined record is never written");
        ensure!(DESIGNED_LOCAL.contains(&"factor_eps_pole"));
        ensure!(!uncontracted::is_setting_code("factor_eps_pole"));
        Ok(())
    })
}
