//! Review of the U7a delta at 6c492c9ddc (zero-comparison,
//! refusal-integrity and values-and-off lenses): the reproducers of the
//! findings, as tests that FAIL on 6c492c9ddc and pass with the fixes. They
//! use only the interfaces of 6c492c9ddc, so this file compiles on it
//! unchanged (the before-run).
//!
//! * refusal-integrity 1: a child Δ filled inside an inert scope (an
//!   insertion surrogate's R*, the dressed scope) never serves an ordinary
//!   K2 child: the child's records are expected, compared (check) or taken
//!   (on) as with a fresh cache;
//! * refusal-integrity 3: a designed `insertion_mass_symbol` decline never
//!   hides a loop-routing error at another vertex of the record;
//! * zero-comparison 4 / refusal-integrity 2 (the ledger side): the close of
//!   a diagram whose ledger counts an error outcome fails, whether or not the
//!   caller propagated the error; so does mode on with a record without a
//!   disposition;
//! * zero-comparison 5: resumed evidence is verified again (graph and
//!   diagram key, the coverage rule recomputed from the saved counts);
//! * Astra U7a finding 3, the ledger side (zero-comparison 3): a diagram
//!   ledger opened inside a special scope (a scope left set) fails.

use std::{cell::RefCell, rc::Rc, sync::Arc};

use eyre::{Result, ensure, eyre};
use symbolica::{
    atom::{Atom, AtomCore},
    function, symbol,
};
use vakint::{
    Vakint,
    p_integral::factorized::{Endpoint, FactorizedEdge},
};

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt, initialisation::test_initialise,
    uv::approx::streamed::with_production_piece_limit,
};

use super::{
    CountertermCache, RStarPlan,
    boundary::{BoundaryMode, with_boundary_mode},
    plan::uv_child_edge_ids,
    recursive_tests::{colour_factored, load, settings},
    uncontracted::{
        self, Captured, InsertionHook, Kind, Mode, Switch, with_capture, with_insertion_hook,
        with_strict, with_switch,
    },
    uncontracted_ledger::{self as ledger, DiagramLedger, Resumed},
    uncontracted_oracle_tests::metric_root,
};

const ROOT_CHILD: [Kind; 2] = [Kind::Root, Kind::Child];
const L2_22_DOT: &str = include_str!("testdata/u7a/l2-d22.dot");
const L3_16_DOT: &str = include_str!("testdata/u7a/l3-d016.dot");
const L3_36_DOT: &str = include_str!("testdata/u7a/l3-d036.dot");

fn check() -> Switch {
    Switch::new(Mode::Check, &ROOT_CHILD)
}

fn on() -> Switch {
    Switch::new(Mode::On, &ROOT_CHILD)
}

fn field<'l>(line: &'l str, key: &str) -> Option<&'l str> {
    line.split_whitespace()
        .find_map(|token| token.strip_prefix(key)?.strip_prefix('='))
}

fn number(line: &str, key: &str) -> usize {
    field(line, key)
        .and_then(|value| value.parse().ok())
        .unwrap_or_else(|| panic!("no count {key}= in {line}"))
}

fn kind_line<'l>(lines: &'l [String], kind: Kind) -> Result<&'l str> {
    let found = lines
        .iter()
        .filter(|line| {
            line.starts_with("UNCONTRACTED_DISPOSITION ")
                && field(line, "kind") == Some(kind.as_str())
        })
        .collect::<Vec<_>>();
    ensure!(
        found.len() == 1,
        "{} {} lines: {lines:#?}",
        found.len(),
        kind.as_str()
    );
    Ok(found[0])
}

fn incidence(edge: &FactorizedEdge, node: usize) -> i64 {
    i64::from(edge.target == Endpoint::Vertex(node))
        - i64::from(edge.source == Endpoint::Vertex(node))
}

fn auxiliary_symbol() -> symbolica::atom::Symbol {
    symbol!("rstar_mass_squared"; Scalar)
}

/// The beta's full-graph closing root of `diagram` under `switch`.
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

/// refusal-integrity 3: the loop components of a root record's internal
/// edge AWAY from a mass-bearing `ct:` node reversed (a misrouted record),
/// with the node's Δ made genuinely mass dependent (× (1 + log M_aux²)): the
/// record must be the routing error, not the designed
/// `insertion_mass_symbol` (at 6c492c9ddc: declined, the diagram passed).
/// Control: the same reversal without the mass dependence is an error too.
#[test]
fn mass_decline_never_masks_a_routing_error_elsewhere() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        for (name, dot) in [("L3_016", L3_16_DOT), ("L3_036", L3_36_DOT)] {
            let diagram = colour_factored(&load(dot)?)?;
            for genuine in [false, true] {
                let done = Rc::new(RefCell::new(None::<(usize, usize)>));
                let hook: InsertionHook = {
                    let done = done.clone();
                    Rc::new(
                        move |kind: Kind,
                              node: usize,
                              delta: &mut Atom,
                              edges: &mut [FactorizedEdge]| {
                            if kind != Kind::Root || done.borrow().is_some() {
                                return;
                            }
                            let loops = edges[0].momentum.len() - 1;
                            let Some(at) = edges.iter().position(|edge| {
                                !edge.is_external()
                                    && incidence(edge, node) == 0
                                    && edge.source != edge.target
                                    && edge.momentum[..loops].iter().any(|c| *c != 0)
                            }) else {
                                return;
                            };
                            for c in &mut edges[at].momentum[..loops] {
                                *c = -*c;
                            }
                            if genuine {
                                let log = function!(
                                    symbolica::atom::Symbol::LOG,
                                    Atom::var(auxiliary_symbol())
                                );
                                *delta = (&*delta * (Atom::num(1) + log)).expand();
                            }
                            *done.borrow_mut() = Some((node, edges[at].id));
                        },
                    )
                };
                let (value, captured) =
                    with_insertion_hook(hook.clone(), || root_run(&diagram, check()))?;
                ensure!(
                    done.borrow().is_some(),
                    "{name}: the mutation found no edge"
                );
                let outcomes = captured
                    .iter()
                    .filter(|a| a.kind == Kind::Root && a.outcome != "uncontracted")
                    .map(|a| a.outcome.clone())
                    .collect::<Vec<_>>();
                ensure!(
                    value.is_err(),
                    "{name} genuine_mass={genuine}: a misrouted record passed the check ({outcomes:?})"
                );
                ensure!(
                    outcomes.iter().any(|o| o.starts_with("error:"))
                        && !outcomes.iter().any(|o| o.contains("insertion_mass_symbol")),
                    "{name} genuine_mass={genuine}: the misrouted record is not an error: {outcomes:?}"
                );
                if genuine {
                    ensure!(
                        outcomes.iter().any(|o| o == "error:routing"),
                        "{name}: the masked misrouting must be the routing error: {outcomes:?}"
                    );
                }
            }
        }
        Ok(())
    })
}

/// refusal-integrity 1: the child {5,6,7,8} of L3 036 filled once inside the
/// dressed scope (an insertion surrogate's R*, whose records are the inert
/// `surrogate` kind) and then looked up by an ordinary R* child of the same
/// shared cache: the ordinary lookup must count, attempt and compare (check)
/// or take (on) its records exactly as with a fresh cache (at 6c492c9ddc:
/// expected=0, no attempt, today's Δ served silently).
#[test]
fn inert_scope_fill_never_serves_an_ordinary_child() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let diagram = colour_factored(&load(L3_36_DOT)?)?;
        let settings = settings();
        let vakint = Vakint::new()?;
        let graph = diagram.to_gamma_loop_graph(None, true)?;
        let model = diagram.model();
        let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
        let child = plan
            .counterterms
            .iter()
            .map(|call| call.argument.support.clone())
            .find(|child| uv_child_edge_ids(&graph, child) == [5, 6, 7, 8])
            .ok_or_else(|| eyre!("no child 5,6,7,8"))?;
        let run = |cache: &Arc<CountertermCache>| {
            RStarPlan::rstar_child_counterterm(
                cache,
                &graph,
                &child,
                model,
                &vakint,
                &settings,
                plan.uv_children,
                plan.uv_child_closure,
            )
        };
        for switch in [check(), on()] {
            with_switch(switch, || {
                with_boundary_mode(BoundaryMode::On, false, || -> Result<()> {
                    let mut seen = Vec::new();
                    for shadowed in [false, true] {
                        let cache = Arc::new(CountertermCache::default());
                        if shadowed {
                            super::expression::in_dressed_scope(|| run(&cache))?;
                        }
                        let _ = ledger::take_lines();
                        let open = DiagramLedger::open("review_ri1", "000036")
                            .ok_or_else(|| eyre!("no ledger"))?;
                        let (value, captured) = with_capture(|| run(&cache));
                        let value = value?;
                        open.close(true)?;
                        let lines = ledger::take_lines();
                        let line = kind_line(&lines, Kind::Child)?.to_owned();
                        seen.push((value, captured.len(), line));
                    }
                    let (fresh, shadowed) = (&seen[0], &seen[1]);
                    ensure!(
                        fresh.0 == shadowed.0 || (&fresh.0 - &shadowed.0).expand().is_zero(),
                        "{}: the Δ differs after the inert fill",
                        switch.mode.as_str()
                    );
                    for key in ["expected", "attempted", "compared", "cache_hit"] {
                        ensure!(
                            number(&fresh.2, key) == number(&shadowed.2, key),
                            "{}: {key} differs after an inert-scope fill:\nfresh    {}\nshadowed {}",
                            switch.mode.as_str(),
                            fresh.2,
                            shadowed.2
                        );
                    }
                    ensure!(
                        fresh.1 == shadowed.1 && number(&shadowed.2, "expected") > 0,
                        "{}: attempts fresh {} shadowed {}: {}",
                        switch.mode.as_str(),
                        fresh.1,
                        shadowed.1,
                        shadowed.2
                    );
                    Ok(())
                })
            })?;
        }
        Ok(())
    })
}

/// zero-comparison 4 (and refusal-integrity 2, its ledger side): a diagram
/// ledger that counts an error outcome fails its close even when the caller
/// reports success (a guard swallowed the error); mode on with a record that
/// reached the hook but has no disposition fails likewise (at 6c492c9ddc
/// both closed Ok).
#[test]
fn close_fails_on_an_error_or_an_undisposed_record() -> Result<()> {
    test_initialise()?;
    for switch in [check(), on()] {
        with_switch(switch, || -> Result<()> {
            let _ = ledger::take_lines();
            let diagram =
                DiagramLedger::open("review_zc4", "000001").ok_or_else(|| eyre!("no ledger"))?;
            ledger::expect(Kind::Root);
            ledger::record_unknown_masters(Kind::Root);
            ensure!(
                diagram.close(true).is_err(),
                "{}: an error outcome closed ok",
                switch.mode.as_str()
            );
            let _ = ledger::take_lines();
            Ok(())
        })?;
    }
    with_switch(on(), || -> Result<()> {
        let diagram =
            DiagramLedger::open("review_zc4", "000002").ok_or_else(|| eyre!("no ledger"))?;
        ledger::expect(Kind::Child);
        ensure!(
            diagram.close(true).is_err(),
            "on: a record without a disposition closed ok"
        );
        let _ = ledger::take_lines();
        Ok(())
    })
}

/// zero-comparison 5: resumed evidence is verified again. A saved child
/// kind edited to `compared=0` (no class: its records undisposed) under the
/// same identity, and another diagram's evidence resumed as this one, are
/// reported unchecked (at 6c492c9ddc both were accepted as checked).
#[test]
fn resumed_evidence_is_verified_again() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let diagram = colour_factored(&load(L2_22_DOT)?)?;
        let graph = diagram.to_gamma_loop_graph(None, true)?.name;
        let evidence = with_switch(check(), || -> Result<serde_json::Value> {
            let open = DiagramLedger::open(&graph, "000022").ok_or_else(|| eyre!("no ledger"))?;
            let (value, _) = root_run(&diagram, check())?;
            let _ = value?;
            open.close(true)?
                .ok_or_else(|| eyre!("no evidence from a check ledger"))
        })?;
        let _ = ledger::take_lines();
        with_switch(check(), || -> Result<()> {
            ensure!(
                ledger::resumed(&graph, "000022", Some(&evidence)) == Resumed::Checked,
                "the untouched evidence is not checked"
            );
            let mut edited = evidence.clone();
            let child = &mut edited["kinds"]["child"];
            ensure!(
                child["expected"].as_u64().unwrap_or(0) > 0,
                "L2 22 child expects records: {child}"
            );
            child["compared"] = serde_json::json!(0);
            child["classes"] = serde_json::json!({});
            let edited_resumed = ledger::resumed(&graph, "000022", Some(&edited));
            ensure!(
                matches!(edited_resumed, Resumed::Unchecked(_)),
                "edited evidence accepted: {edited_resumed:?}"
            );
            let foreign = ledger::resumed(&graph, "000018", Some(&evidence));
            ensure!(
                matches!(foreign, Resumed::Unchecked(_)),
                "another diagram's evidence accepted: {foreign:?}"
            );
            let _ = ledger::take_lines();
            Ok(())
        })
    })
}

/// zero-comparison 3 (Astra U7a finding 3, the ledger's own check): a
/// diagram ledger the beta opens while this thread is still inside the
/// dressed scope (a scope left set) fails its close: every root record of
/// the diagram would be the inert `surrogate` kind and the kinds would
/// report ok with nothing expected (at 6c492c9ddc the close was ok).
#[test]
fn a_ledger_opened_inside_a_scope_fails() -> Result<()> {
    test_initialise()?;
    for switch in [check(), on()] {
        with_switch(switch, || -> Result<()> {
            let closed = super::expression::in_dressed_scope(|| -> Result<bool> {
                let open = DiagramLedger::open("review_zc3", "000003")
                    .ok_or_else(|| eyre!("no ledger"))?;
                Ok(open.close(true).is_err())
            })?;
            let _ = ledger::take_lines();
            ensure!(
                closed,
                "{}: a ledger opened inside the dressed scope closed ok",
                switch.mode.as_str()
            );
            Ok(())
        })?;
    }
    let _ = uncontracted::classify(false);
    Ok(())
}
