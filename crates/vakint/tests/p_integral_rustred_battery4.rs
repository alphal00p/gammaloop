//! Tier A of the Forcer battery at four loops (M3-V9; M3 plan, gate G7): RustRed
//! on M3's scope 1-4 p library against M3-V7's frozen Forcer fixtures,
//! through vakint's public API, without FORM (`support/p_battery4.rs`).
//!
//! Full profile only (M3 D-M3-10): `.config/nextest.toml` excludes this
//! binary from `test_vakint` by name; `just test_vakint` runs a few
//! four-loop cases (`src/p_integral/rustred_battery_tests.rs`). Every test
//! but the pin of the case counts needs `VAKINT_TEST_P_LIBRARY_1_4` (the
//! package of the PP-M3b pin) and reports a skip on stderr without it. The
//! evidence run (G7 Tier A at four loops), from the repository root, in a
//! profile without termination and without retries (the default profile
//! ends a test after 60 s); each case reports its seconds and peak RSS on
//! stderr ([`timed`]):
//!
//! ```text
//! VAKINT_TEST_P_LIBRARY_1_4=<package> cargo nextest run --cargo-profile ci-optim \
//!     --locked -p vakint --features vakint/feynkit-ingress -P test_gammaloop \
//!     --retries 0 --test p_integral_rustred_battery4 --no-fail-fast --no-capture
//! ```
//!
//! - The battery: the 54 cases of `p4l_battery.json` in parts of
//!   [`BATTERY_PART`], each in every mode through its frozen depth and at 11
//!   terms (M3 I-7, or the `Oep` failure on both backends), under the default
//!   arithmetic (`modular_polynomial`).
//! - The anchors: the 17 cases of `p4l_anchors.json`; together they reach
//!   every master of the four-loop basis in `Symbolic` and `Numerical` mode.
//! - The gate part of `p4l_grid.json`: its 55 sampled sectors in parts of
//!   [`GRID_PART`] (the 165 benchmark tops are Tier C's).
//! - `exact` and `modular_sampled`: one battery case per family and the
//!   anchors, every mode through the frozen depth.

#![cfg(feature = "feynkit-ingress")]

#[path = "support/p_backends.rs"]
mod p_backends;
#[path = "support/p_battery4.rs"]
mod p_battery4;
#[path = "support/p_fixtures.rs"]
mod p_fixtures;
#[path = "support/p_integrals.rs"]
mod p_integrals;
#[path = "support/p_library4.rs"]
mod p_library4;
#[path = "support/p_rustred.rs"]
mod p_rustred;

use std::collections::BTreeSet;

use p_battery4::*;
use vakint::p_integral::forcer_fixtures::Finalized;

/// At most this many battery cases per test.
const BATTERY_PART: usize = 9;
/// At most this many grid cases per test.
const GRID_PART: usize = 11;

/// The four-loop fixtures of Tier A, without an engine: the battery's 54
/// cases over 11 families, 26 of them `Oep` at 11 terms and 28 with an
/// 11-term value; 17 anchors, one per master of the basis but `no` and
/// `t105` (three-loop masters, which the four-loop tops carry); the grid's
/// gate part; the arithmetic subset of 11 + 17 cases.
#[test]
fn the_four_loop_cases_of_tier_a_are_pinned() {
    let battery = fixture4(BATTERY4, BATTERY4_CASES);
    let oep = (battery.cases.iter())
        .filter(|case| matches!(case.eleven_terms, Some(Finalized::MissingMasterOrders(_))))
        .count();
    let symbolic = (battery.cases.iter())
        .filter(|case| matches!(case.eleven_terms, Some(Finalized::Symbolic(_))))
        .count();
    assert_eq!(
        (oep, symbolic),
        (BATTERY4_OEP_CASES, BATTERY4_CASES - BATTERY4_OEP_CASES)
    );
    let anchored = fixture4(ANCHORS4, ANCHORS4_CASES)
        .cases
        .into_iter()
        .filter_map(|case| case.master)
        .collect::<BTreeSet<_>>();
    let expected = FOUR_LOOP_MASTERS
        .into_iter()
        .filter(|master| !["no", "t105"].contains(master))
        .map(String::from)
        .collect::<BTreeSet<_>>();
    assert_eq!(anchored, expected);
    assert_eq!(grid4_gate().len(), GRID4_GATE_CASES);
    assert_eq!(arithmetic_subset4().len(), 11 + ANCHORS4_CASES);
}

/// Part `part` of the four-loop battery: every case agrees with Forcer
/// through its frozen depth and at 11 terms.
fn assert_battery_part(part: usize, parts: usize) {
    let Some(engine) = p_battery4_engine(&format!("four-loop battery part {}", part + 1)) else {
        return;
    };
    let battery = fixture4(BATTERY4, BATTERY4_CASES);
    for case in part_of(&battery.cases, part, parts, BATTERY_PART) {
        timed(format!("{BATTERY4}: {}", case.name), || {
            assert_case_agrees(engine, BATTERY4, case, None);
            assert_eleven_terms(engine, BATTERY4, case);
        });
    }
}

macro_rules! battery_parts {
    ($($name:ident: $part:expr),* $(,)?) => {$(
        #[test]
        fn $name() {
            assert_battery_part($part, 6);
        }
    )*};
}

battery_parts!(
    the_four_loop_battery_part_1_agrees_with_forcer: 0,
    the_four_loop_battery_part_2_agrees_with_forcer: 1,
    the_four_loop_battery_part_3_agrees_with_forcer: 2,
    the_four_loop_battery_part_4_agrees_with_forcer: 3,
    the_four_loop_battery_part_5_agrees_with_forcer: 4,
    the_four_loop_battery_part_6_agrees_with_forcer: 5,
);

/// Every anchor agrees with Forcer through its depth, and the masters of
/// RustRed's values, each reached in `Symbolic` and `Numerical` mode by its
/// agreement, are every master of the four-loop basis.
#[test]
fn the_four_loop_anchors_agree_with_forcer_and_reach_every_master() {
    let Some(engine) = p_battery4_engine("the four-loop anchors") else {
        return;
    };
    let mut reached = BTreeSet::new();
    for case in fixture4(ANCHORS4, ANCHORS4_CASES).cases {
        let masters = timed(format!("{ANCHORS4}: {}", case.name), || {
            assert_case_agrees(engine, ANCHORS4, &case, None)
        });
        if let Some(master) = &case.master {
            assert!(masters.contains(master), "{}: {masters:?}", case.name);
        }
        reached.extend(masters);
    }
    let all = FOUR_LOOP_MASTERS
        .map(String::from)
        .into_iter()
        .collect::<BTreeSet<_>>();
    assert_eq!(reached, all);
}

/// Part `part` of the gate part of the four-loop grid.
fn assert_grid_part(part: usize, parts: usize) {
    let Some(engine) = p_battery4_engine(&format!("four-loop grid part {}", part + 1)) else {
        return;
    };
    for case in part_of(&grid4_gate(), part, parts, GRID_PART) {
        timed(format!("{GRID4}: {}", case.name), || {
            assert_case_agrees(engine, GRID4, case, None)
        });
    }
}

macro_rules! grid_parts {
    ($($name:ident: $part:expr),* $(,)?) => {$(
        #[test]
        fn $name() {
            assert_grid_part($part, 5);
        }
    )*};
}

grid_parts!(
    the_four_loop_grid_gate_part_1_agrees_with_forcer: 0,
    the_four_loop_grid_gate_part_2_agrees_with_forcer: 1,
    the_four_loop_grid_gate_part_3_agrees_with_forcer: 2,
    the_four_loop_grid_gate_part_4_agrees_with_forcer: 3,
    the_four_loop_grid_gate_part_5_agrees_with_forcer: 4,
);

/// [`arithmetic_subset4`] agrees with Forcer under `arithmetic`.
fn assert_subset_under(arithmetic: &str) {
    let Some(engine) = p_battery4_engine(&format!("the four-loop subset under {arithmetic}"))
    else {
        return;
    };
    for (file, case) in arithmetic_subset4() {
        timed(format!("{file}: {} ({arithmetic})", case.name), || {
            assert_case_agrees(engine, file, &case, Some(arithmetic))
        });
    }
}

/// One case per family and the anchors under `exact`.
#[test]
fn the_exact_arithmetic_agrees_with_forcer_on_the_four_loop_subset() {
    assert_subset_under(OTHER_ARITHMETICS[0]);
}

/// One case per family and the anchors under `modular_sampled`.
#[test]
fn the_sampled_arithmetic_agrees_with_forcer_on_the_four_loop_subset() {
    assert_subset_under(OTHER_ARITHMETICS[1]);
}
