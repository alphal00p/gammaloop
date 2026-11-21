//! The coordinator's pre-review of the U7a delta at 86de2ec899
//! (`PREREVIEW-U7a-delta.md`): reproducers that use only the interfaces of
//! 86de2ec899 and two test hooks (a K3 value changed above its valid order,
//! a panic in the IR hard child route's outer step), so that each FAILS on
//! 86de2ec899 (with the hooks) and passes with the fixes:
//!
//! * finding 3: under `GAMMALOOP_IR_HARD_CHILD_GUARD=1` a panic or an error
//!   caught after the K3 records is an error (it was
//!   `declined:ir_hard_route:panic`, designed and exceptable);
//! * finding 4: in mode on a route decline caused by the new K3 values is an
//!   error (it fell back as `declined:ir_hard_route:<head>`), and check mode
//!   runs the route's value-dependent steps on its K3 values (it passed
//!   today's values only and never saw it);
//! * finding 2 (the ledger's half): a reused record without matching
//!   evidence asks for a recompute, never `outcome=unchecked`;
//! * finding 5: evidence whose comparisons did not all agree is not
//!   matching evidence;
//! * finding 7: mode on applies the all-declined rule.

use eyre::{Result, ensure, eyre};

use crate::{initialisation::test_initialise, uv::approx::streamed::with_production_piece_limit};

use super::{
    ir_dirac,
    uncontracted::{Kind, Mode, Switch, with_switch},
    uncontracted_ledger::{
        self as ledger, DiagramLedger, Exception, Exceptions, Resumed, Section, SectionEnd,
    },
    uncontracted_ledger_tests::{
        L2_22, L3_36, agree, count, disposition, field, graph_name, lines_of, run,
    },
};

const ALL_KINDS: [Kind; 3] = [Kind::Root, Kind::Child, Kind::IrHard];
const ROOT_CHILD: [Kind; 2] = [Kind::Root, Kind::Child];

/// Finding 3: the IR hard child route's codes and the guard's failures are
/// never designed, so an exception list naming one is refused at start-up;
/// the exhausted budget stays the one designed section code.
#[test]
fn exceptions_refuse_the_route_and_guard_failure_codes() {
    for code in [
        "ir_hard_route:panic",
        "ir_hard_route:error",
        "ir_hard_route:diagnostic",
        "ir_hard_route:moments",
        "ir_hard_route:depth",
        "ir_hard_guard:panic",
        "ir_hard_guard:error",
    ] {
        let text = format!("B_self_energy_L3_006\tir-hard\t{code}\twhy\n");
        assert!(
            ledger::parse_exceptions(&text).is_err(),
            "an exception naming {code} was accepted"
        );
    }
    assert!(
        ledger::parse_exceptions("B_self_energy_L3_006\tir-hard\tir_hard_guard:budget\twhy\n")
            .is_ok()
    );
}

/// One section with `records` K3 records evaluated (`taken` on, `compared`
/// check) and closed by `end`, in a diagram ledger: the ir-hard line and
/// whether the diagram's close passed.
fn section(mode: Mode, end: SectionEnd, records: usize) -> Result<(String, bool)> {
    with_switch(
        Switch::new(mode, &ALL_KINDS),
        || -> Result<(String, bool)> {
            let _ = ledger::take_lines();
            let diagram = DiagramLedger::open("prereview_section", "000003")
                .ok_or_else(|| eyre!("no ledger"))?;
            let section = Section::open().ok_or_else(|| eyre!("no section"))?;
            for _ in 0..records {
                ledger::expect(Kind::IrHard);
                ledger::record_value(Kind::IrHard, mode);
            }
            let _ = section.close(end);
            let closed = diagram.close(true);
            let lines = ledger::take_lines();
            Ok((
                disposition(&lines, Kind::IrHard)?.to_owned(),
                closed.is_ok(),
            ))
        },
    )
}

/// Finding 3 at the ledger: a panic (or an error) the guard caught AFTER the
/// K3 records ran fails the diagram in on and in check; at 86de2ec899 the
/// records became `declined:ir_hard_route:panic` (on, outcome ok) or stayed
/// `compared` (check, outcome ok).
#[test]
fn a_guard_failure_after_the_records_fails_the_diagram() -> Result<()> {
    test_initialise()?;
    for mode in [Mode::On, Mode::Check] {
        for reason in ["panic:index_out_of_bounds", "error:taylor_failed"] {
            let (line, closed) = section(mode, SectionEnd::GuardFailure(reason.into()), 2)?;
            let head = reason.split(':').next().unwrap_or_default();
            ensure!(
                field(&line, "error") == Some(format!("ir_hard_guard:{head}:2").as_str())
                    && field(&line, "declined") == Some("none")
                    && field(&line, "outcome") == Some("error")
                    && !closed,
                "{} {reason}: {line} (close ok: {closed})",
                mode.as_str()
            );
        }
    }
    Ok(())
}

/// Finding 4 at the ledger: in mode on, a route decline after the K3 values
/// were taken is the error `ir_hard_route:<head>` and fails the diagram (at
/// 86de2ec899: `declined:ir_hard_route:<head>`, outcome ok); check mode
/// keeps its comparisons (the route had today's values).
#[test]
fn a_route_decline_after_taken_values_fails_mode_on() -> Result<()> {
    test_initialise()?;
    for reason in ["moments:not_a_polynomial", "dirac_shape", "depth:3"] {
        let head = reason.split(':').next().unwrap_or_default();
        let (line, closed) = section(Mode::On, SectionEnd::RouteDeclined(reason.into()), 2)?;
        ensure!(
            field(&line, "error") == Some(format!("ir_hard_route:{head}:2").as_str())
                && field(&line, "declined") == Some("none")
                && field(&line, "outcome") == Some("error")
                && !closed,
            "on {reason}: {line} (close ok: {closed})"
        );
        let (line, closed) = section(Mode::Check, SectionEnd::RouteDeclined(reason.into()), 2)?;
        ensure!(
            count(&line, "compared") == 2 && field(&line, "outcome") == Some("ok") && closed,
            "check {reason}: {line}"
        );
    }
    Ok(())
}

/// The ir-hard disposition line of a run and its error text.
fn failed_ir_hard(run: &super::uncontracted_ledger_tests::Run) -> (String, String) {
    let error = match &run.value {
        Ok(_) => String::new(),
        Err(error) => format!("{error:#}"),
    };
    let line = disposition(&run.lines, Kind::IrHard)
        .map(str::to_owned)
        .unwrap_or_else(|error| format!("{error:#}"));
    (line, error)
}

/// Finding 3 at the full hook (L3 36, the IR hard child route,
/// `GAMMALOOP_IR_HARD_CHILD_GUARD=1`): a panic in the route's outer step,
/// after every K3 record of the correction ran, fails the diagram in on and
/// in check. At 86de2ec899 the guard declined the correction to the closed
/// path and the diagram passed (on: `declined=ir_hard_route:panic`, check:
/// the records `compared`, outcome ok).
#[test]
fn a_panic_after_the_k3_records_fails_on_and_check() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        for mode in [Mode::On, Mode::Check] {
            let run = ir_dirac::with_child_guard(true, || {
                ledger::with_outer_panic(|| {
                    run(L3_36, Switch::new(mode, &ALL_KINDS), None, Some("000036"))
                })
            })?;
            let (line, error) = failed_ir_hard(&run);
            ensure!(
                run.value.is_err()
                    && error.contains("GAMMALOOP_IR_HARD_CHILD_GUARD caught panic")
                    && field(&line, "error").is_some_and(|e| e.contains("ir_hard_guard:panic"))
                    && field(&line, "outcome") == Some("error"),
                "{}: the caught panic did not fail the diagram: {line}\n{error}",
                mode.as_str()
            );
        }
        Ok(())
    })
}

/// Finding 4 at the full hook (L3 36, mode on): every K3 value changed
/// ABOVE its valid order by a term the route's moment split refuses
/// (`with_k3_mutation`): the route declines after taking the values, which
/// fails the diagram with `error=ir_hard_route:moments`. Unchanged, the same
/// run passes and its K3 records are taken. At 86de2ec899 the mutated run
/// passed with `declined=ir_hard_route:moments` (the closed path computed
/// the correction).
#[test]
fn mode_on_fails_when_the_route_refuses_the_k3_values() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let on = Switch::new(Mode::On, &ALL_KINDS);
        let clean = run(L3_36, on, None, Some("000036"))?;
        let _ = clean.value?;
        let line = disposition(&clean.lines, Kind::IrHard)?;
        ensure!(
            field(line, "outcome") == Some("ok")
                && field(line, "error") == Some("none")
                && field(line, "declined") == Some("none"),
            "the unchanged on run: {line}"
        );
        let mutated = ledger::with_k3_mutation(|| run(L3_36, on, None, Some("000036")))?;
        let (line, error) = failed_ir_hard(&mutated);
        ensure!(
            mutated.value.is_err()
                && error.contains("declined after it took the uncontracted K3 values")
                && field(&line, "error").is_some_and(|e| e.starts_with("ir_hard_route:moments:"))
                && field(&line, "declined") == Some("none")
                && field(&line, "outcome") == Some("error"),
            "a route decline on the new values fell back: {line}\n{error}"
        );
        Ok(())
    })
}

/// Finding 4 in check mode (L3 36): the same change (above the valid order,
/// so every comparison agrees) is now seen by check mode, which runs the
/// route's value-dependent steps on its K3 values: `error=ir_hard_route:
/// moments:1`, the diagram fails. At 86de2ec899 check passed (the route had
/// today's values).
#[test]
fn check_mode_sees_the_route_refuse_the_k3_values() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let check = Switch::new(Mode::Check, &ALL_KINDS);
        let mutated = ledger::with_k3_mutation(|| run(L3_36, check, None, Some("000036")))?;
        let (line, error) = failed_ir_hard(&mutated);
        let records = mutated
            .captured
            .iter()
            .filter(|attempt| attempt.kind == Kind::IrHard)
            .count();
        ensure!(
            mutated.value.is_err()
                && error.contains("declines (moments:")
                && !error.contains("UNCONTRACTED_CHECK")
                && records > 0
                && count(&line, "compared") == records
                && count(&line, "equal") == records
                && field(&line, "error") == Some("ir_hard_route:moments:1")
                && field(&line, "outcome") == Some("error"),
            "check mode passed K3 values the route refuses ({records} records): {line}\n{error}"
        );
        Ok(())
    })
}

/// Finding 2 (the ledger's half): a reused record without matching
/// evidence prints one `UNCONTRACTED_RESUME ... outcome=recompute` line and
/// no disposition line (the beta recomputes the diagram, which prints its
/// own). At 86de2ec899 it printed `outcome=unchecked` lines, and the run
/// exited 0 with nothing compared.
#[test]
fn a_reused_record_without_evidence_is_recomputed() -> Result<()> {
    test_initialise()?;
    with_switch(Switch::new(Mode::Check, &ROOT_CHILD), || -> Result<()> {
        let _ = ledger::take_lines();
        ensure!(ledger::resumed("g", "000022", None) == Resumed::Unchecked("no_evidence".into()));
        let lines = ledger::take_lines();
        ensure!(
            lines.len() == 1
                && lines[0].starts_with("UNCONTRACTED_RESUME ")
                && field(&lines[0], "outcome") == Some("recompute")
                && field(&lines[0], "reason") == Some("no_evidence")
                && lines_of(&lines, "UNCONTRACTED_DISPOSITION ").is_empty(),
            "{lines:#?}"
        );
        Ok(())
    })
}

/// Finding 5: saved counts whose comparisons did not all agree (root
/// `equal` below `compared`; a child Δ-cache hit counted `hit:compared`
/// without an agreeing comparison) are not matching evidence. At 86de2ec899
/// both were resumed as checked (the verdict never read `equal`).
#[test]
fn evidence_with_disagreeing_comparisons_is_not_matching() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let check = Switch::new(Mode::Check, &ROOT_CHILD);
        let fresh = run(L2_22, check, None, Some("000022"))?;
        let _ = fresh.value?;
        let evidence = fresh.evidence.ok_or_else(|| eyre!("no evidence"))?;
        let graph = graph_name(L2_22)?;
        with_switch(check, || -> Result<()> {
            ensure!(ledger::resumed(&graph, "000022", Some(&evidence)) == Resumed::Checked);
            let mut unequal = evidence.clone();
            let equal = unequal["kinds"]["root"]["equal"].as_u64().unwrap_or(0);
            ensure!(equal > 0, "{evidence:#}");
            unequal["kinds"]["root"]["equal"] = serde_json::json!(equal - 1);
            let mut hit = evidence.clone();
            let child = &mut hit["kinds"]["child"];
            child["expected"] = serde_json::json!(child["expected"].as_u64().unwrap_or(0) + 1);
            child["classes"]["hit:compared"] = serde_json::json!(1);
            for (edited, kind) in [(unequal, "root"), (hit, "child")] {
                let resumed = ledger::resumed(&graph, "000022", Some(&edited));
                ensure!(
                    resumed == Resumed::Unchecked(format!("verdict:{kind}:error")),
                    "{kind}: {resumed:?}"
                );
            }
            let _ = ledger::take_lines();
            Ok(())
        })
    })
}

/// Finding 7: mode on applies the all-declined rule. Every root record of
/// L2 22 declined by a designed code (as a setting-wide decline does) fails
/// the diagram in mode on unless the run's exception list names it; named,
/// it passes on today's path. At 86de2ec899 mode on passed with nothing
/// taken (outcome ok).
#[test]
fn mode_on_fails_an_all_declined_kind_unless_named() -> Result<()> {
    with_production_piece_limit(|| -> Result<()> {
        test_initialise()?;
        let on = Switch::new(Mode::On, &ROOT_CHILD);
        let off = run(L2_22, Switch::OFF, None, None)?.value?;
        let declined = || {
            ledger::with_fault(Kind::Root, "additional_normalization", true, || {
                run(L2_22, on, None, Some("000022"))
            })
        };
        let failed = declined()?;
        let root = disposition(&failed.lines, Kind::Root)?;
        ensure!(
            failed.value.is_err()
                && field(root, "outcome") == Some("uncovered")
                && field(root, "declined")
                    .is_some_and(|d| d.starts_with("additional_normalization:")),
            "mode on passed an all-declined root: {root}"
        );
        let named = Exceptions {
            entries: vec![Exception {
                graph: graph_name(L2_22)?,
                kind: Kind::Root,
                code: "additional_normalization".into(),
                reason: "the test's setting-wide decline".into(),
            }],
            digest: Some("test".into()),
        };
        let excepted = ledger::with_exceptions(named, declined)?;
        let root = disposition(&excepted.lines, Kind::Root)?;
        ensure!(
            field(root, "outcome") == Some("excepted"),
            "the named decline: {root}"
        );
        ensure!(
            agree(&excepted.value?, &off),
            "the excepted run must take today's path"
        );
        Ok(())
    })
}
