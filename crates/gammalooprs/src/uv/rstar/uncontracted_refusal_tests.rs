//! U7a delta, lane B (what a record's outcome is): designed refusals against
//! errors, and the scopes that decide a record's kind (Astra U7a verdict
//! findings 2 and 3, `keep the rebase test narrow`; helper VH round 2 R2-I3,
//! R2-I4, R2-M4, R2-M5; U7-PLAN requirements 13 and 15). Rule: ONLY a named
//! designed refusal falls back to today's path; everything else is
//! `outcome=error`, loud in check and in strict on.
//!
//! * A nonzero residual of the Δ momentum-conservation step in a component
//!   the record must conserve is a routing error, never
//!   `insertion_mass_symbol`: the signed mutation at the full hook (check and
//!   strict on fail), with its controls (a genuine mass dependence: the
//!   designed decline; independent crown legs: the record is conserved and
//!   evaluated).
//! * A panic caught inside the dressed scope (or the quotient scope) leaves
//!   no scope behind: the same thread's next root and child are classified
//!   and executed as root and child.
//! * `one_loop_tree_cut` is the emitter's own decline of the single
//!   propagator shape only (the closing line alone, a massive tadpole after
//!   the Taylor step); every other one-loop record stays on the new path.
//! * A vakint `Closure` refusal without a WH-RECORD-FORMAT section-8 prefix
//!   is an error (requirement 13, interim until the rebase onto helper VH's
//!   frozen REFUSALS list), "every line is massive" among them.
//! * The record builder's impossible states (`sigma`, `root_crown_paired`,
//!   `root_normalisation_momentum`) are errors.

use std::{cell::RefCell, collections::BTreeMap, rc::Rc};

use eyre::{Result, ensure, eyre};
use linnet::half_edge::subgraph::SuBitGraph;
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
    plan::uv_child_edge_ids,
    recursive_tests::{colour_factored, load, settings},
    top_mass::TopMassTerms,
    uncontracted::{
        self, Captured, DESIGNED_LOCAL, Decline, InsertionHook, Kind, Mode, Outcome, Site, Switch,
        split_root_factor, with_capture, with_insertion_hook, with_strict, with_switch,
    },
    uncontracted_oracle_tests::{atoms_equal, auxiliary_mass, metric_root, parse},
};

const ROOT_CHILD: [Kind; 2] = [Kind::Root, Kind::Child];

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
/// record. Check mode and strict on must FAIL the diagram with
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
        for (mode, switch) in [("check", check()), ("strict on", on())] {
            flips.borrow_mut().clear();
            let (value, captured) =
                with_insertion_hook(hook.clone(), || l3_016_run(&case, &case.delta, switch))?;
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
        let (checked, captured) = l3_016_run(&case, &delta, check())?;
        let checked = checked?;
        ensure!(checked == off, "check returns today's value");
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
            let (value, captured) = root_run(&diagram, check())?;
            let _ = value?;
            evaluated_except(what, &captured, &["declined:one_loop_tree_cut"])?;
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
// U7-PLAN requirement 13 (interim): no section-8 prefix, no designed refusal
// ---------------------------------------------------------------------------

/// U7-PLAN requirement 13, interim rule (until the rebase onto helper VH's
/// frozen REFUSALS list): a vakint `Closure` refusal is designed only with a
/// WH-RECORD-FORMAT section-8 prefix (`Taylor degree `, `the projector
/// contains Q(`, the open-spinor port `open slot ` + `is not a Lorentz
/// slot`); any other `Closure` detail is an error, among them the entry's
/// wrapped internal failures and "every line is massive: not a closure" (no
/// section-8 row: the emitter declines that shape itself,
/// `one_loop_tree_cut`, so vakint never sees it on a correct record).
#[test]
fn closure_without_a_section_8_prefix_is_an_error() {
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
    ] {
        let decline = Decline::vakint(&closure(detail));
        assert!(decline.designed, "{detail}: {decline:?}");
        assert_eq!(decline.code, code);
    }
    for detail in [
        "every line is massive: not a closure",
        "lowering: no key plan",
        "the lines do not form one closing line among massless ones",
        "truncating the insertion product: x",
    ] {
        let decline = Decline::vakint(&closure(detail));
        assert!(!decline.designed, "{detail}: {decline:?}");
        assert_eq!(decline.code, "vakint:Closure", "{detail}");
    }
    let prefixes = ["Taylor degree ", "the projector contains Q(", "open slot "];
    for row in uncontracted::DESIGNED_VAKINT
        .iter()
        .filter(|row| row.0 == "Closure")
    {
        assert!(
            prefixes.contains(&row.1),
            "the Closure row {:?} has no section-8 prefix",
            row.1
        );
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
