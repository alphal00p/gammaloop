//! `GAMMALOOP_IR_DIRAC_CLOSURE=check` with the K3 records of
//! `GAMMALOOP_UNCONTRACTED_NUMERATOR=on:ir-hard` (U7b; Astra verdict #10,
//! finding 2) on the UNMODIFIED L4 4966 (colour factored, no fixture factor
//! injected), inside the R* of its UV child {2,3,5,7,8,9,10,11} as the
//! production beta computes it: the IR correction with the Dirac-closed hard
//! part γ = {7,8,9,10,11} (ω = 1, slash) takes the child route, its closed
//! arm from the uncontracted records (valid through `ε^0`, the quotient's one
//! extra loop), its open arm from today's contracted integration (which
//! carries more orders).
//!
//! * Unchanged, with `GAMMALOOP_IR_HARD_CHILD_GUARD` off and on: every
//!   `IR_DIRAC_CHECK` is compared through the justified order `ε^0` and
//!   passes (the arms first differ at `ε^1`), nothing declines, and the
//!   child's Δ equals the contracted baseline's (uncontracted off).
//! * A controlled change of one closed coefficient at or below `ε^0`
//!   ([`super::ir_dirac::with_closed_mutation`]): the check fails at that
//!   order, and the diagram fails in both guard states (the guard never
//!   declines a failed comparison).
//! * A change only above `ε^0`: passes in both guard states, and the Δ is the
//!   unchanged one (the orders above the justified one never reach it).

use std::collections::BTreeMap;

use eyre::{Result, ensure, eyre};
use symbolica::atom::Atom;
use vakint::Vakint;

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt, initialisation::test_initialise,
    uv::approx::streamed::with_production_piece_limit,
};

use super::{
    CountertermCache, RStarPlan, UvChildClosure,
    boundary::{BoundaryMode, with_boundary_mode},
    ir_dirac::{self, IrDiracMode},
    recursive_tests::{colour_factored, load, settings},
    uncontracted::{self, Captured, Kind, Mode, Switch, with_capture, with_strict, with_switch},
};

const L4_4966: &str = include_str!("testdata/u7b/l4-d4966.dot");
/// The UV child of L4 4966 whose R* holds the Dirac-closed K3 site.
const CHILD: [usize; 8] = [2, 3, 5, 7, 8, 9, 10, 11];
/// The remainder loops of the K3 record with the `ct:edges-7-8` insertion
/// (term 1): the `all_rem1` group of `IR_DIRAC_CHECK`.
const REM1: usize = 1;

/// One R* of [`CHILD`] under `IR_DIRAC_CLOSURE=check`, strict.
struct Run {
    value: Result<Atom>,
    captured: Vec<Captured>,
    counts: BTreeMap<String, usize>,
}

impl Run {
    fn count(&self, key: &str) -> usize {
        self.counts.get(key).copied().unwrap_or(0)
    }

    /// Every `IR_DIRAC_CHECK` passed, at least one ran.
    fn checks_equal(&self) -> bool {
        self.count("check") > 0 && self.count("check") == self.count("check_equal")
    }

    fn value(&self) -> Result<&Atom> {
        self.value.as_ref().map_err(|error| eyre!("{error:#}"))
    }
}

fn child_run(
    vakint: &Vakint,
    switch: Switch,
    guard: bool,
    mutation: Option<(usize, i64)>,
) -> Result<Run> {
    let diagram = colour_factored(&load(L4_4966)?)?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let settings = settings();
    let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
    let child = plan
        .counterterms
        .iter()
        .map(|call| &call.argument.support)
        .find(|child| super::plan::uv_child_edge_ids(&graph, child) == CHILD)
        .ok_or_else(|| eyre!("{}: no UV child {CHILD:?}", graph.name))?;
    let cache = std::sync::Arc::new(super::CountertermCache::default());
    let run = || {
        with_switch(switch, || {
            with_strict(true, || {
                with_boundary_mode(BoundaryMode::On, false, || {
                    with_capture(|| {
                        RStarPlan::rstar_child_counterterm(
                            &cache,
                            &graph,
                            child,
                            model,
                            vakint,
                            &settings,
                            plan.uv_children,
                            plan.uv_child_closure,
                        )
                    })
                })
            })
        })
    };
    let _ = ir_dirac::take_counts();
    let (value, captured) = ir_dirac::with_ir_dirac_mode(IrDiracMode::Check, || {
        ir_dirac::with_child_guard(guard, || match mutation {
            Some((loops, order)) => ir_dirac::with_closed_mutation(loops, order, run),
            None => run(),
        })
    });
    let counts = ir_dirac::take_counts();
    eprintln!(
        "IR_DIRAC_TEST guard={guard} mutation={mutation:?} ok={} counts={counts:?} attempts={:?}",
        value.is_ok(),
        captured
            .iter()
            .map(|a| format!(
                "{}:{}:loops={}:order={}",
                a.kind.as_str(),
                a.outcome,
                a.loops,
                a.order
            ))
            .collect::<Vec<_>>()
    );
    if let Err(error) = &value {
        eprintln!("IR_DIRAC_TEST error: {error:#}");
    }
    Ok(Run {
        value,
        captured,
        counts,
    })
}

fn on() -> Switch {
    Switch::new(Mode::On, &[Kind::IrHard])
}

/// The contracted baseline (uncontracted off; both arms today's) and the
/// unchanged `on:ir-hard` run with `guard`: the baseline's checks pass; the
/// on run's closed arm is the uncontracted records' (a K3 record with a
/// value, valid through `ε^0`), its checks are compared through `ε^0` and
/// pass, and its Δ is the baseline's.
fn unchanged(vakint: &Vakint, guard: bool) -> Result<(Atom, Run)> {
    let baseline = child_run(vakint, Switch::OFF, false, None)?;
    let reference = baseline.value()?.clone();
    ensure!(baseline.checks_equal(), "baseline: {:?}", baseline.counts);
    let run = child_run(vakint, on(), guard, None)?;
    ensure!(
        run.captured
            .iter()
            .any(|a| a.kind == Kind::IrHard && a.outcome == "uncontracted" && a.loops == REM1),
        "no uncontracted K3 record of {REM1} loop: {:?}",
        run.captured.iter().map(|a| &a.outcome).collect::<Vec<_>>()
    );
    ensure!(
        run.captured
            .iter()
            .filter(|a| a.kind == Kind::IrHard)
            .all(|a| a.order == 0),
        "the K3 records' valid order"
    );
    let value = run.value()?;
    ensure!(
        run.checks_equal() && run.count("check_through:eps^0") > 0,
        "guard={guard}: an IR_DIRAC_CHECK failed or was not compared through eps^0: {:?}",
        run.counts
    );
    ensure!(
        uncontracted::agree(value, &reference)?,
        "guard={guard}: the on:ir-hard Δ differs from the contracted baseline\non       {value}\nbaseline {reference}"
    );
    Ok((reference, run))
}

/// A change at `ε^order` (at or below the justified `ε^0`) of the closed
/// `all_rem1` group fails the diagram with `guard`: a diagnostic failure
/// (no fallback took it), the check's first difference at `ε^order`.
fn in_range_fails(guard: bool) -> Result<()> {
    with_production_piece_limit(|| {
        test_initialise()?;
        let vakint = Vakint::new()?;
        for order in [0, -1] {
            let run = child_run(&vakint, on(), guard, Some((REM1, order)))?;
            ensure!(run.count("mutated") > 0, "order {order}: no change applied");
            let Err(error) = &run.value else {
                return Err(eyre!(
                    "guard={guard}: a closed coefficient changed at eps^{order} (in the valid range) passed: {:?}",
                    run.counts
                ));
            };
            let text = format!("{error:#}");
            ensure!(
                text.contains("IR_DIRAC_CHECK")
                    && text.contains(&format!("first_difference=eps^{order} "))
                    && text.contains("through=eps^0 "),
                "guard={guard} order {order}: not the failed comparison at eps^{order}: {text}"
            );
            ensure!(
                run.count(&format!("check_first:eps^{order}")) > 0,
                "guard={guard} order {order}: {:?}",
                run.counts
            );
        }
        Ok(())
    })
}

/// A change only above the justified order (from `ε^1`) passes with `guard`
/// and leaves the Δ unchanged.
fn above_range_passes(guard: bool) -> Result<()> {
    with_production_piece_limit(|| {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let (reference, _) = unchanged(&vakint, guard)?;
        let run = child_run(&vakint, on(), guard, Some((REM1, 1)))?;
        ensure!(run.count("mutated") > 0, "no change applied");
        let value = run.value()?;
        ensure!(
            run.checks_equal(),
            "guard={guard}: a change above eps^0 failed a check: {:?}",
            run.counts
        );
        ensure!(
            uncontracted::agree(value, &reference)?,
            "guard={guard}: a change above eps^0 reached the Δ\nchanged   {value}\nreference {reference}"
        );
        Ok(())
    })
}

#[test]
fn ir_dirac_check_k3_l4_4966_guard_off() -> Result<()> {
    with_production_piece_limit(|| {
        test_initialise()?;
        unchanged(&Vakint::new()?, false).map(|_| ())
    })
}

#[test]
fn ir_dirac_check_k3_l4_4966_guard_on() -> Result<()> {
    with_production_piece_limit(|| {
        test_initialise()?;
        unchanged(&Vakint::new()?, true).map(|_| ())
    })
}

#[test]
fn ir_dirac_check_k3_in_range_change_fails_guard_off() -> Result<()> {
    in_range_fails(false)
}

#[test]
fn ir_dirac_check_k3_in_range_change_fails_guard_on() -> Result<()> {
    in_range_fails(true)
}

#[test]
fn ir_dirac_check_k3_above_range_change_passes_guard_off() -> Result<()> {
    above_range_passes(false)
}

#[test]
fn ir_dirac_check_k3_above_range_change_passes_guard_on() -> Result<()> {
    above_range_passes(true)
}

/// The guard around the child route declines an ordinary error, a panic and
/// a designed decline as before, and returns a diagnostic failure (a failed
/// comparison, also under a context) as the error, its text unchanged.
#[test]
fn child_guard_returns_a_diagnostic_failure() {
    use eyre::WrapErr;
    let passed = ir_dirac::guarded_unless_diagnostic(None, None, || Ok(Ok(7_u8)));
    assert!(matches!(passed, Ok(Ok(7))));
    let declined: Result<std::result::Result<u8, String>> =
        ir_dirac::guarded_unless_diagnostic(None, None, || Ok(Err("depth:3".to_owned())));
    assert_eq!(declined.unwrap(), Err("depth:3".to_owned()));
    let ordinary: Result<std::result::Result<u8, String>> =
        ir_dirac::guarded_unless_diagnostic(None, None, || Err(eyre!("taylor failed here")));
    assert_eq!(
        ordinary.unwrap(),
        Err("error:taylor_failed_here".to_owned())
    );
    let panic: Result<std::result::Result<u8, String>> =
        ir_dirac::guarded_unless_diagnostic(None, None, || panic!("boom"));
    assert_eq!(panic.unwrap(), Err("panic:boom".to_owned()));
    let failed: Result<std::result::Result<u8, String>> =
        ir_dirac::guarded_unless_diagnostic(None, None, || {
            Err(ir_dirac::diagnostic_failure(eyre!(
                "IR_DIRAC_CHECK: differs"
            )))
            .wrap_err("while computing a child")
        });
    let error = failed.unwrap_err();
    assert!(ir_dirac::is_diagnostic_failure(&error));
    assert_eq!(
        format!("{error:#}"),
        "while computing a child: IR_DIRAC_CHECK: differs"
    );
}

/// L4 4966, its UV child [`CHILD`], the settings and the plan's children
/// route.
fn l4_4966_child<T>(
    f: impl FnOnce(
        &crate::graph::Graph,
        &linnet::half_edge::subgraph::SuBitGraph,
        &crate::model::Model,
        &crate::uv::UVgenerationSettings,
        super::UvChildrenRoute,
    ) -> Result<T>,
) -> Result<T> {
    let diagram = colour_factored(&load(L4_4966)?)?;
    let graph = diagram.to_gamma_loop_graph(None, true)?;
    let model = diagram.model();
    let settings = settings();
    let plan = RStarPlan::with_automatic_edge(&graph, model, &settings)?;
    let child = plan
        .counterterms
        .iter()
        .map(|call| &call.argument.support)
        .find(|child| super::plan::uv_child_edge_ids(&graph, child) == CHILD)
        .ok_or_else(|| eyre!("{}: no UV child {CHILD:?}", graph.name))?;
    f(&graph, child, model, &settings, plan.uv_children)
}

/// A failed comparison inside an R* computed through the Δ-cache (the
/// decorated quotients of `ir_outer_integrated`, inside the child guard)
/// stays a diagnostic failure across the cache: on the miss that computes it
/// (inside the guard: returned, not declined) and on a later hit of the
/// stored failure, its text unchanged; an ordinary failure stays ordinary
/// (declined by the guard).
#[test]
fn a_diagnostic_failure_survives_the_counterterm_cache() -> Result<()> {
    test_initialise()?;
    l4_4966_child(|graph, child, model, _, _| {
        let cache = CountertermCache::default();
        let computed = std::cell::Cell::new(0_usize);
        let lookup = |key: &str, diagnostic: bool| {
            cache.counterterm_with(graph, child, model, key.to_owned(), || {
                computed.set(computed.get() + 1);
                let error = eyre!("IR_DIRAC_CHECK: differs");
                Err(if diagnostic {
                    ir_dirac::diagnostic_failure(error)
                } else {
                    error
                })
            })
        };
        let miss =
            ir_dirac::guarded_unless_diagnostic(None, None, || lookup("diagnostic", true).map(Ok));
        let error = match miss {
            Err(error) => error,
            Ok(outcome) => {
                return Err(eyre!(
                    "a failed comparison through the cache was not returned by the child guard: {outcome:?}"
                ));
            }
        };
        ensure!(ir_dirac::is_diagnostic_failure(&error), "miss: {error:#}");
        ensure!(
            format!("{error:#}") == "IR_DIRAC_CHECK: differs",
            "miss: {error:#}"
        );
        let hit = lookup("diagnostic", false)
            .err()
            .ok_or_else(|| eyre!("hit: no error"))?;
        ensure!(computed.get() == 1, "the stored failure was recomputed");
        ensure!(ir_dirac::is_diagnostic_failure(&hit), "hit: {hit:#}");
        ensure!(
            format!("{hit:#}") == "IR_DIRAC_CHECK: differs",
            "hit: {hit:#}"
        );
        let ordinary =
            ir_dirac::guarded_unless_diagnostic(None, None, || lookup("ordinary", false).map(Ok))?;
        ensure!(
            ordinary == Err("error:IR_DIRAC_CHECK:_differs".to_owned()),
            "an ordinary failure: {ordinary:?}"
        );
        Ok(())
    })
}

/// `DIRAC_CLOSURE_CHECK` (`GAMMALOOP_UV_CHILD_CLOSURE=dirac-check`) under the
/// child guard, the children's Δ through the Δ-cache, as for a decorated
/// quotient's R* inside `ir_outer_integrated`: unchanged, the R* of
/// [`CHILD`] checks its Dirac-closed children and passes; with their closed Δ
/// changed ([`super::plan::with_dirac_closure_change`]) the failed comparison
/// is returned by the guard as a diagnostic failure, never declined.
#[test]
fn a_failed_dirac_closure_check_survives_the_child_guard() -> Result<()> {
    use std::sync::atomic::Ordering;
    with_production_piece_limit(|| {
        test_initialise()?;
        let vakint = Vakint::new()?;
        let closure = UvChildClosure::parse("isotropic,dirac-check")?;
        l4_4966_child(|graph, child, model, settings, route| {
            let run = |change: bool| {
                let cache = std::sync::Arc::new(CountertermCache::default());
                super::plan::with_dirac_closure_change(change, || {
                    ir_dirac::guarded_unless_diagnostic(None, None, || {
                        RStarPlan::rstar_child_counterterm(
                            &cache, graph, child, model, &vakint, settings, route, closure,
                        )
                        .map(Ok)
                    })
                })
            };
            let checks = || {
                super::plan::DIRAC_CHECKS
                    .each_ref()
                    .map(|count| count.load(Ordering::Relaxed))
            };
            let before = checks();
            let unchanged = run(false);
            let after = checks();
            ensure!(
                matches!(unchanged, Ok(Ok(_))),
                "unchanged: {:?}",
                unchanged.map_err(|error| format!("{error:#}"))
            );
            ensure!(
                after[0] > before[0] && after[2] - before[2] == after[0] - before[0],
                "no Dirac-closed child checked, or a check failed: {before:?} -> {after:?}"
            );
            let error = match run(true) {
                Err(error) => error,
                Ok(outcome) => {
                    return Err(eyre!(
                        "a failed DIRAC_CLOSURE_CHECK was declined or passed by the child guard: {outcome:?}"
                    ));
                }
            };
            let text = format!("{error:#}");
            ensure!(
                ir_dirac::is_diagnostic_failure(&error) && text.contains("DIRAC_CLOSURE_CHECK: "),
                "not the failed comparison: {text}"
            );
            Ok(())
        })
    })
}
