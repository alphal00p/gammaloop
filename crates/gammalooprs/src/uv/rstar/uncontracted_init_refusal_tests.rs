//! Helper VH's C7 (every exponent in a marked record an integer literal; no
//! exp, sqrt or log): the written statement's settings exceptions
//! (results/consolidation/format-c5-c15/c7-normalisation.md sections 7-8,
//! U7-PLAN item 33) are refused with `GAMMALOOP_UNCONTRACTED_NUMERATOR` on
//! or check, and a Δ that would carry a non-local or non-literal form into a
//! record's `insertion` line is the error `insertion_nonlocal`.
//!
//! These tests use only the interfaces of 292b5362a7 (the U7a delta after
//! its pre-review fixes), so this file compiles there unchanged and each
//! test FAILS there (the before-run) and passes with the refusals:
//!
//! * E2: `GAMMALOOP_IR_QUOTIENT` other than `rstar` (`massive`, `off`,
//!   `check`, `auto`) is refused at the switch's initialisation (the beta's
//!   start-up check); `rstar`, unset, and the switch off are accepted.
//! * E1: `uv.vakint.normalization` other than `MSbar` is an error of every
//!   record (`vakint_normalization`), in check and in on (at 292b5362a7 the
//!   records were built and evaluated under it); MSbar is unchanged.
//! * E3: a vakint evaluation method other than `rustred` and `forcer` is an
//!   error of every record (`vakint_evaluation_methods`), in check and in
//!   on; the beta's two orders are unchanged.
//!   `additional_normalization` other than "1" stays the designed decline
//!   `additional_normalization` (confirmed here; passes before and after).

use eyre::{Result, ensure, eyre};
use symbolica::atom::{Atom, AtomCore};
use vakint::Vakint;

use crate::{
    feyngen::feynkit::FeynmanDiagramGammaLoopExt,
    initialisation::test_initialise,
    uv::{UVgenerationSettings, approx::streamed::with_production_piece_limit},
};

use super::{
    RStarPlan,
    boundary::{BoundaryMode, with_boundary_mode},
    recursive_tests::{colour_factored, load, settings},
    uncontracted::{self, Captured, Kind, Mode, Switch, with_capture, with_switch},
    uncontracted_ledger::{self as ledger, DiagramLedger},
    uncontracted_ledger_tests::{L2_22, agree, disposition, field, lines_of},
    uncontracted_oracle_tests::metric_root,
};

/// Root records only: the first record a diagram builds is a root one.
const ROOT: [Kind; 1] = [Kind::Root];

/// One closed-root R* evaluation of `dot` under `switch` and the UV settings
/// `uv`, inside a diagram ledger: its value, every attempt and the ledger
/// lines.
struct Run {
    value: Result<Atom>,
    captured: Vec<Captured>,
    lines: Vec<String>,
}

fn run(dot: &str, uv: &UVgenerationSettings, switch: Switch) -> Result<Run> {
    let loaded = colour_factored(&load(dot)?)?;
    let graph = loaded.to_gamma_loop_graph(None, true)?;
    let model = loaded.model();
    let plan = RStarPlan::with_automatic_edge(&graph, model, uv)?;
    let closed_root = metric_root(&graph);
    let callback = |tensor: &Atom| -> Result<Atom> {
        Ok(tensor.replace_multiple(&closed_root.external_routing) * &closed_root.factor)
    };
    let vakint = Vakint::new()?;
    let _ = (ledger::take_lines(), ledger::take_timing_lines());
    let _ = uncontracted::take_summaries();
    let (value, captured) = with_switch(switch, || {
        with_boundary_mode(BoundaryMode::On, false, || {
            with_capture(|| {
                let opened = DiagramLedger::open(&graph.name, "000099");
                let value = plan.projected_counterterm_closing_root_scalar(
                    &graph,
                    model,
                    &vakint,
                    uv,
                    &callback,
                    &closed_root,
                );
                match opened {
                    None => value,
                    Some(opened) => {
                        let closed = opened.close(value.is_ok());
                        match (value, closed) {
                            (Ok(value), Ok(_)) => Ok(value),
                            (Err(error), _) | (Ok(_), Err(error)) => Err(error),
                        }
                    }
                }
            })
        })
    });
    let _ = uncontracted::take_summaries();
    Ok(Run {
        value,
        captured,
        lines: ledger::take_lines(),
    })
}

/// `run` failed its diagram with the error `code` of a root record: the
/// error text, the root kind's ledger line and the diagram's summary.
fn fails_with(run: &Run, code: &str, what: &str) -> Result<()> {
    let error = match &run.value {
        Ok(value) => {
            return Err(eyre!(
                "{what}: the diagram passed ({value}); outcomes {:?}",
                run.captured
                    .iter()
                    .map(|a| (a.kind.as_str(), a.id, a.outcome.clone()))
                    .collect::<Vec<_>>()
            ));
        }
        Err(error) => format!("{error:#}"),
    };
    ensure!(
        error.contains(&format!("outcome=error:{code}")),
        "{what}: not the error {code}: {error}"
    );
    ensure!(
        run.captured
            .iter()
            .any(|a| a.kind == Kind::Root && a.outcome == format!("error:{code}")),
        "{what}: no root attempt with outcome error:{code}"
    );
    let root = disposition(&run.lines, Kind::Root)?;
    ensure!(
        field(root, "outcome") == Some("error")
            && field(root, "error").is_some_and(|e| e.starts_with(&format!("{code}:")))
            && field(root, "declined") == Some("none"),
        "{what}: {root}"
    );
    let summaries = lines_of(&run.lines, "UNCONTRACTED_SUMMARY ");
    ensure!(
        summaries.len() == 1 && field(summaries[0], "outcome") == Some("error"),
        "{what}: {summaries:#?}"
    );
    Ok(())
}

/// Tests that set the switch's environment: only in a process of their own.
fn own_process() -> Result<()> {
    ensure!(
        std::env::var_os("NEXTEST").is_some(),
        "this test sets process environment variables: run it under cargo nextest (one process per test)"
    );
    Ok(())
}

/// The switch's initialisation in this process with `switch` and
/// `GAMMALOOP_IR_QUOTIENT=quotient` (`None`: unset).
fn initialise_with_quotient(switch: &str, quotient: Option<&str>) -> Result<Result<()>> {
    own_process()?;
    unsafe {
        std::env::set_var(uncontracted::VARIABLE, switch);
        match quotient {
            Some(value) => std::env::set_var("GAMMALOOP_IR_QUOTIENT", value),
            None => std::env::remove_var("GAMMALOOP_IR_QUOTIENT"),
        }
    }
    Ok(super::validate_uncontracted())
}

/// E2: refused in check and in on, naming the variable and its value.
fn quotient_is_refused(switch: &str, quotient: &str) -> Result<()> {
    let error = initialise_with_quotient(switch, Some(quotient))?
        .err()
        .ok_or_else(|| eyre!("{switch} initialised with GAMMALOOP_IR_QUOTIENT={quotient}"))?;
    let text = format!("{error}");
    ensure!(
        text.contains(&format!("GAMMALOOP_IR_QUOTIENT={quotient} is refused")),
        "not the quotient refusal: {text}"
    );
    // the switch is not read either (the first use panics, as for any
    // refused value)
    let read = std::panic::catch_unwind(uncontracted::switch);
    ensure!(read.is_err(), "the switch was read: {read:?}");
    Ok(())
}

#[test]
fn e2_ir_quotient_massive_is_refused_in_check() -> Result<()> {
    quotient_is_refused("check:root,child", "massive")
}

#[test]
fn e2_ir_quotient_massive_is_refused_in_on() -> Result<()> {
    quotient_is_refused("on:root,child", "massive")
}

#[test]
fn e2_ir_quotient_off_is_refused() -> Result<()> {
    quotient_is_refused("check:root,child", "off")
}

#[test]
fn e2_ir_quotient_check_is_refused() -> Result<()> {
    quotient_is_refused("on:root", "check")
}

/// `auto` is massive for every quotient of at most four loops.
#[test]
fn e2_ir_quotient_auto_is_refused() -> Result<()> {
    quotient_is_refused("check:root,child", "auto")
}

/// Controls (pass before and after): the production route, unset or
/// explicit, and the switch off with any route.
#[test]
fn e2_ir_quotient_rstar_is_accepted() -> Result<()> {
    initialise_with_quotient("check:root,child", Some("rstar"))?
}

#[test]
fn e2_ir_quotient_unset_is_accepted() -> Result<()> {
    initialise_with_quotient("on:root,child", None)?
}

#[test]
fn e2_ir_quotient_massive_with_the_switch_off_is_accepted() -> Result<()> {
    initialise_with_quotient("off", Some("massive"))?
}

/// E1: every normalisation but MSbar is an error of every record, in check
/// and in on; MSbar (the beta's) is the control, whose check compares and
/// whose on agrees with off.
#[test]
fn e1_a_normalisation_other_than_msbar_is_an_error() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let msbar = settings();
        let off = run(L2_22, &msbar, Switch::OFF)?.value?;
        for (mode, switch) in [
            ("check", Switch::new(Mode::Check, &ROOT)),
            ("on", Switch::new(Mode::On, &ROOT)),
        ] {
            let control = run(L2_22, &msbar, switch)?;
            ensure!(
                agree(control.value.as_ref().map_err(|e| eyre!("{e:#}"))?, &off),
                "{mode}: MSbar differs from off"
            );
            for normalization in ["GScheme", "FMFTandMATAD"] {
                let mut uv = settings();
                uv.vakint.normalization = normalization.into();
                fails_with(
                    &run(L2_22, &uv, switch)?,
                    "vakint_normalization",
                    &format!("{mode} {normalization:?}"),
                )?;
            }
        }
        Ok(())
    })
}

/// `additional_normalization` other than "1" stays the DESIGNED decline of
/// every record (not an error, not a settings refusal): each root record of
/// L2 22 is `declined:additional_normalization`, in check and in on. (The
/// diagram then fails by the all-declined rule unless an exception names it:
/// `uncontracted_prereview_tests`.)
#[test]
fn additional_normalization_stays_declined() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let mut uv = settings();
        uv.vakint.additional_normalization = "-1".into();
        for (mode, switch) in [
            ("check", Switch::new(Mode::Check, &ROOT)),
            ("on", Switch::new(Mode::On, &ROOT)),
        ] {
            let declined = run(L2_22, &uv, switch)?;
            let roots = declined
                .captured
                .iter()
                .filter(|a| a.kind == Kind::Root)
                .collect::<Vec<_>>();
            ensure!(
                !roots.is_empty()
                    && roots
                        .iter()
                        .all(|a| a.outcome == "declined:additional_normalization"),
                "{mode}: {:?}",
                roots.iter().map(|a| a.outcome.clone()).collect::<Vec<_>>()
            );
            let root = disposition(&declined.lines, Kind::Root)?;
            ensure!(
                field(root, "error") == Some("none")
                    && field(root, "declined")
                        .is_some_and(|d| d.starts_with("additional_normalization:")),
                "{mode}: {root}"
            );
        }
        Ok(())
    })
}

/// E3: an evaluation method outside {rustred, forcer} is an error of every
/// record (here after the two, so vakint itself would never reach it: the
/// setting alone is refused), in check and in on; the beta's order (RustRed
/// first) is the control (Forcer first: the settings test of
/// `uncontracted_ledger_tests`).
#[test]
fn e3_an_evaluation_method_outside_rustred_and_forcer_is_an_error() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let off = run(L2_22, &settings(), Switch::OFF)?.value?;
        for (mode, switch) in [
            ("check", Switch::new(Mode::Check, &ROOT)),
            ("on", Switch::new(Mode::On, &ROOT)),
        ] {
            for order in [["rustred", "forcer"]] {
                let mut uv = settings();
                uv.vakint.evaluation_methods = order.iter().map(|m| m.to_string()).collect();
                let control = run(L2_22, &uv, switch)?;
                ensure!(
                    agree(control.value.as_ref().map_err(|e| eyre!("{e:#}"))?, &off),
                    "{mode} {order:?}: differs from off"
                );
            }
            for extra in ["matad", "fmft", "alphaloop"] {
                let mut uv = settings();
                uv.vakint.evaluation_methods =
                    vec!["rustred".into(), "forcer".into(), extra.into()];
                fails_with(
                    &run(L2_22, &uv, switch)?,
                    "vakint_evaluation_methods",
                    &format!("{mode} {extra}"),
                )?;
            }
        }
        Ok(())
    })
}
