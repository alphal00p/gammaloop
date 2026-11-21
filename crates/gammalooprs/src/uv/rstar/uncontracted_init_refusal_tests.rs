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

use eyre::{Result, ensure, eyre};

use super::uncontracted;

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
