//! IRH-DEPTH round 8 on the G-G stack: the uncontracted K3 records of the
//! IR hard parts (`GAMMALOOP_UNCONTRACTED_NUMERATOR`, kind `ir-hard`) under
//! the round-8 rules, on B_self_energy_L3_006 as the earlier rounds
//! ([`super`]), through the typed scalar root (the numerator prefactor
//! multiplied as beta does, `projected_counterterm_scalar`), K3 on and in
//! check, STRICT off and on, the uncontracted strict switch on (a K3 error
//! fails the run).
//!
//! * C1: the K3 record's factors are the graph's own numerator factors (the
//!   same objects the contracted path integrates; the weight is not in a
//!   record: the caller's root projection multiplies the normalised one).
//!   vakint substitutes `D -> 4 - 2ε` in the record's integrand before any
//!   expansion in ε (`p_integral::factorized`, closure.rs / evaluate.rs), and
//!   the record's pole check (`factor_eps_pole`) reads every field at
//!   d = 4 − 2ε. The hard factor −1 spelled `1/(1+(d−4)/ε)` gives the
//!   literal's value on the hard vertex and on a hard edge; the quotient
//!   form `ε/(ε+d−4)` on a hard factor keeps round 5's refusal (no record).
//! * The declines: a correction the child route declines before H(γ) (a
//!   declared k > 0: k3j1) builds no K3 record and the generic value is used.
//!   One declined AFTER H(γ) (C4, a row insertion that is not ε-free) keeps
//!   the stack's ledger rule (`uncontracted_ledger::Section`, review of the
//!   U7a delta, finding 4): in check mode the generic value is used and the
//!   records' comparisons stand; in mode on, where the route took the K3
//!   values, the decline is the error `error:ir_hard_route:..` (never a
//!   silent fallback after the new values were taken).

use eyre::{Result, ensure, eyre};
use symbolica::atom::{Atom, AtomCore};

use super::super::{
    closed::with_ir_hard,
    depth::{self, DepthCounts, DepthTest, RowHook},
    uncontracted::{Captured, Kind, Mode, Switch, with_capture, with_strict, with_switch},
};
use super::{Case, EPS, FIXED, L3_006, k, same};
use crate::utils::GS;

const DIM: &str = "gammalooprs::{}::dim";

/// The three spellings of −1 at d = 4 − 2ε (as round 8b's tests).
fn spellings() -> [(&'static str, String); 3] {
    [
        ("literal", "-1".to_owned()),
        ("cross", format!("1/(1+({DIM}-4)/{EPS})")),
        ("form", format!("{EPS}/({EPS}+{DIM}-4)")),
    ]
}

/// `dot` with the quoted `field` of the first line starting with `  {term} `
/// multiplied by `factor`.
fn times(dot: &str, term: &str, field: &str, factor: &str) -> String {
    let key = format!("{field}=\"");
    let mut count = 0;
    let out = dot
        .split_inclusive('\n')
        .map(|line| {
            if count == 0 && line.starts_with(&format!("  {term} ")) {
                let mut from = 0;
                while let Some(found) = line[from..].find(&key) {
                    let start = from + found;
                    let before = line[..start].chars().last();
                    if before.is_some_and(|c| c.is_ascii_lowercase() || c == '_') {
                        from = start + key.len();
                        continue;
                    }
                    let body = start + key.len();
                    let end = body + line[body..].find('"').expect("a closing quote");
                    count += 1;
                    return format!(
                        "{}({})*({factor}){}",
                        &line[..body],
                        &line[body..end],
                        &line[end..]
                    );
                }
                line.to_owned()
            } else {
                line.to_owned()
            }
        })
        .collect();
    assert_eq!(count, 1, "{term} {field}");
    out
}

/// `base` with the numerator of `term` and the aggregate graph numerator
/// multiplied by `factor` (the review's inputs).
fn spelled(base: &str, term: &str, factor: &str) -> String {
    let dot = times(base, term, "numerator", factor);
    times(&dot, "graph", "numerator", factor)
}

/// One run of the typed scalar root on `route`.
struct K3Run {
    value: Result<Atom>,
    counts: DepthCounts,
    /// The K3 attempts reported (records written).
    k3: Vec<Captured>,
}

fn typed_run(
    case: &Case,
    route: &'static str,
    mode: Mode,
    strict: bool,
    hook: Option<RowHook>,
) -> Result<K3Run> {
    let prefactor = case.graph.global_prefactor.num.clone();
    let test = DepthTest {
        strict: Some(strict),
        ..FIXED
    };
    let evaluate = || {
        with_ir_hard(route, || {
            case.plan.projected_counterterm_scalar(
                &case.graph,
                case.model(),
                &case.vakint,
                &case.settings,
                &prefactor,
            )
        })
    };
    let ((value, counts), captured) = with_switch(Switch::new(mode, &[Kind::IrHard]), || {
        with_strict(true, || {
            with_capture(|| {
                depth::with_depth_test(test, || {
                    depth::collect(|| match hook {
                        Some(hook) => depth::with_row_hook(hook, evaluate),
                        None => evaluate(),
                    })
                })
            })
        })
    });
    let k3 = captured
        .into_iter()
        .filter(|captured| captured.kind == Kind::IrHard)
        .collect::<Vec<_>>();
    eprintln!(
        "IRH_R8_K3 route={route} mode={mode:?} strict={strict} k3={} ok={} counts={counts:?}",
        k3.len(),
        value.is_ok()
    );
    Ok(K3Run { value, counts, k3 })
}

/// Today's (generic) value: the closed path, no record.
fn generic(case: &Case) -> Result<Atom> {
    typed_run(case, "closed", Mode::Off, false, None)?.value
}

/// C1 on the K3 records: every spelling of −1 on the hard vertex v0 and on
/// the hard edge v2 → v3 (both inside the K3 record) gives the literal's
/// value, which is −1 times the base's (vakint substitutes d = 4 − 2ε in the
/// record's integrand before it expands in ε).
#[test]
fn round8_k3_a_hard_spellings_equal_the_literal() -> Result<()> {
    let base = Case::new(L3_006, None)?;
    let base_value = generic(&base)?;
    drop(base);
    let mut failures = Vec::new();
    for term in ["v0", "v2 ->"] {
        let mut reference = None::<Atom>;
        for (name, factor) in spellings() {
            let case = Case::new(&spelled(L3_006, term, &factor), None)?;
            let today = generic(&case)?;
            match &reference {
                None => {
                    ensure!(
                        same(&today, &(-&base_value))?,
                        "{term}: the literal −1 is not −1 times the base"
                    );
                    reference = Some(today.clone());
                }
                Some(reference) => {
                    if !same(&today, reference)? {
                        failures.push(format!("{term} {name} generic: differs from the literal"));
                    }
                }
            }
            for mode in [Mode::On, Mode::Check] {
                for strict in [false, true] {
                    let cell = format!("{term} {name} K3 {mode:?} strict={strict}");
                    let run = typed_run(&case, "child", mode, strict, None)?;
                    let value = match &run.value {
                        Ok(value) => value,
                        Err(error) => {
                            failures.push(format!("{cell}: error {error:#}"));
                            continue;
                        }
                    };
                    if !same(value, reference.as_ref().expect("the literal"))? {
                        failures.push(format!("{cell}: differs from the literal"));
                    }
                    // The quotient form on a HARD factor keeps round 5's
                    // refusal (before H(γ): no record).
                    let refused = name == "form";
                    let disposition = if refused {
                        run.counts.taken == 0
                            && run.counts.reason == Some("symbol_in_denominator")
                            && run.k3.is_empty()
                    } else {
                        run.counts.taken >= 1
                            && run.counts.declined == 0
                            && run.counts.exceeded == 0
                            && !run.k3.is_empty()
                    };
                    if !disposition {
                        failures.push(format!(
                            "{cell}: disposition {:?}, {} K3 record(s)",
                            run.counts,
                            run.k3.len()
                        ));
                    }
                }
            }
        }
    }
    ensure!(failures.is_empty(), "{failures:#?}");
    Ok(())
}

/// A correction declined before H(γ) (k3j1: v1·ε⁻¹ declares k = 1) builds no
/// K3 record; the value is the generic one, K3 on and in check, STRICT off
/// and on.
#[test]
fn round8_k3_b_a_declined_correction_writes_no_record() -> Result<()> {
    let case = Case::new(&k(3, 1), None)?;
    let today = generic(&case)?;
    for mode in [Mode::On, Mode::Check] {
        for strict in [false, true] {
            let cell = format!("k3j1 K3 {mode:?} strict={strict}");
            let run = typed_run(&case, "child", mode, strict, None)?;
            let value = run.value.as_ref().map_err(|error| eyre!("{cell}: {error:#}"))?;
            ensure!(same(value, &today)?, "{cell}: differs from the generic value");
            ensure!(
                run.counts.taken == 0
                    && run.counts.declined >= 1
                    && run.counts.declared_order >= 1
                    && run.counts.reason == Some("declared_order"),
                "{cell}: {:?}",
                run.counts
            );
            ensure!(run.k3.is_empty(), "{cell}: {} K3 record(s)", run.k3.len());
        }
    }
    Ok(())
}

/// A correction declined AFTER H(γ) was integrated (C4: here every row
/// insertion is made ε-dependent by the test hook), on the stack's ledger
/// rule: in check mode the generic value is used and the K3 records keep
/// their comparisons; in mode on the route took the K3 values, so the
/// decline fails the diagram (`error:ir_hard_route:`, the U7a delta's
/// finding 4), never a silent fallback. The same run without the hook takes
/// the correction.
#[test]
fn round8_k3_c_a_late_decline_on_the_ledger_rule() -> Result<()> {
    let case = Case::new(L3_006, None)?;
    let today = generic(&case)?;
    for mode in [Mode::On, Mode::Check] {
        for strict in [false, true] {
            let cell = format!("L3_006 K3 {mode:?} strict={strict}");
            let taken = typed_run(&case, "child", mode, strict, None)?;
            let value = taken.value.as_ref().map_err(|error| eyre!("{cell}: {error:#}"))?;
            ensure!(same(value, &today)?, "{cell}: taken: differs from the generic value");
            ensure!(
                taken.counts.taken >= 1 && !taken.k3.is_empty(),
                "{cell}: taken: {:?}, {} K3 record(s)",
                taken.counts,
                taken.k3.len()
            );
            let hook: RowHook = Box::new(|_, rows| {
                for row in rows.iter_mut() {
                    row.1 = &row.1 * (Atom::num(1) + Atom::var(GS.dim_epsilon));
                }
            });
            let late = typed_run(&case, "child", mode, strict, Some(hook))?;
            ensure!(
                late.counts.taken == 0
                    && late.counts.declined >= 1
                    && late.counts.reason == Some("insertion_not_eps_free"),
                "{cell}: declined: {:?}",
                late.counts
            );
            match mode {
                Mode::Check => {
                    let value = late.value.as_ref().map_err(|error| eyre!("{cell}: {error:#}"))?;
                    ensure!(same(value, &today)?, "{cell}: declined: differs from the generic value");
                    ensure!(
                        !late.k3.is_empty(),
                        "{cell}: the compared K3 records are not reported"
                    );
                }
                _ => {
                    let error = late
                        .value
                        .as_ref()
                        .err()
                        .ok_or_else(|| eyre!("{cell}: a late decline after taken K3 values fell back silently"))?;
                    let text = format!("{error:#}");
                    ensure!(
                        text.contains("ir_hard_route") && text.contains("insertion_not_eps_free"),
                        "{cell}: {text}"
                    );
                }
            }
        }
    }
    Ok(())
}
