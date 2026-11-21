//! The beta's REAL R* records WITH INSERTIONS through the closure entry,
//! against the value their contracted path returns at that call.
//!
//! Rules (shared by every real-record test of this series): a missing fixture
//! or an empty set is a FAILURE; the number of records compared is asserted
//! against the count written here; the comparison is exact, coefficient by
//! coefficient of the Laurent series in ε through the record's
//! `valid_through` (a SIGNED order: negative orders compare, they do not
//! hang) and no further, with a NONZERO reference through that order (a
//! record whose stated value is the exact atom `0` requires an exact zero). The
//! records of the repository are compared always; the sets that live outside
//! it are `#[ignore]` tests that take their root from ONE environment
//! variable, `VAKINT_TEST_WH_FIXTURES` (the beta's `wh-fixtures` directory),
//! and fail when it is unset or incomplete. Five-loop records need the 1-4
//! loop package of `VAKINT_TEST_P_LIBRARY_1_4` (required where needed).
use super::*;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
    RustRedMasterMode, Vakint, VakintSettings,
};
use std::path::{Path, PathBuf};
use symbolica::domains::rational::Rational;

/// The R* sites' settings (ε = `gammalooprs::ε`, μ_r² = `gammalooprs::μᵣ²`,
/// MSbar, RustRed only, FORM absent); the record's `depth` decides the terms.
fn site_settings() -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Symbolic,
            },
        )]),
        form_exe_path: "/nonexistent/helper-vh/form".into(),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        number_of_terms_in_epsilon_expansion: 3,
        use_dot_product_notation: true,
        n_cores: 1,
        epsilon_symbol: "gammalooprs::ε".into(),
        mu_r_sq_symbol: "gammalooprs::μᵣ²".into(),
        ..Default::default()
    }
}

fn value_of(value_text: &str) -> Atom {
    crate::Vakint::initialize_vakint_symbols();
    let line = value_text
        .lines()
        .find_map(|l| l.strip_prefix("value\t"))
        .expect("a value line");
    vk_parse!(line).unwrap()
}

/// The signed order of `valid_through eps^<n> ...`.
fn valid_through(value_text: &str) -> i64 {
    let line = value_text
        .lines()
        .find_map(|l| l.strip_prefix("valid_through\t"))
        .expect("a valid_through line");
    line.trim_start_matches("eps^")
        .split_whitespace()
        .next()
        .unwrap()
        .parse()
        .expect("an integer order")
}

/// The coefficients of `expression` in ε at the orders `<= order` (signed),
/// from its Laurent series truncated at `order`.
fn coefficients_through(expression: &Atom, order: i64) -> Vec<(i64, Atom)> {
    let epsilon = vk_symbol!("gammalooprs::ε");
    if expression.is_zero() {
        return Vec::new();
    }
    let series = expression
        .series(epsilon, Atom::Zero.as_view(), Rational::from(order))
        .unwrap();
    series
        .terms()
        .into_iter()
        .map(|(power, coefficient)| (power.numerator().to_i64().unwrap(), coefficient.expand()))
        .filter(|(n, c)| *n <= order && !c.is_zero())
        .collect()
}

/// Exact equality through `order` (signed) and no further; the reference
/// must be nonzero through `order`. Returns the number of nonzero orders of
/// the reference compared.
fn assert_equal_through(got: &Atom, expected: &Atom, order: i64, what: &str) -> usize {
    let reference = coefficients_through(expected, order);
    assert!(
        !reference.is_empty(),
        "{what}: the reference is zero through eps^{order}"
    );
    let difference = coefficients_through(&(got - expected).expand(), order);
    assert!(
        difference.is_empty(),
        "{what}: differs at eps^{:?}: {:?}\nentry {got}\nbeta  {expected}",
        difference.iter().map(|(n, _)| n).collect::<Vec<_>>(),
        difference.first().map(|(_, c)| c.to_string())
    );
    reference.len()
}

/// One record through the entry against its value (`engine` with the
/// package when the record needs it); returns the entry's seconds.
fn compare(vakint: &Vakint, input: &str, value_text: &str, what: &str) -> f64 {
    let value = FactorizedNumerator::from_text(input).unwrap_or_else(|e| panic!("{what}: {e}"));
    let again = FactorizedNumerator::from_text(&value.to_text().unwrap()).unwrap();
    assert_eq!(again, value, "{what}: round trip");
    let started = std::time::Instant::now();
    let got = evaluate_factorized_closure(
        vakint,
        &site_settings(),
        std::slice::from_ref(&value),
        &[value.prefactor.total()],
    )
    .unwrap_or_else(|e| panic!("{what}: {e}"));
    let seconds = started.elapsed().as_secs_f64();
    // a record whose stated value is the exact atom 0 (a remainder whose every
    // pattern integrates to zero): the entry must return exactly 0; any other
    // value must be nonzero through its order
    let stated_zero = value_text.lines().any(|l| l.trim_end() == "value\t0");
    let orders = if stated_zero {
        assert!(
            got.value.expand().is_zero(),
            "{what}: stated 0, entry {}",
            got.value
        );
        0
    } else {
        assert_equal_through(
            &got.value,
            &value_of(value_text),
            valid_through(value_text),
            what,
        )
    };
    eprintln!(
        "INSERTION_RECORD {what} insertions={} pole_order={} orders_compared={orders} seconds={seconds:.3}",
        value.insertions.len(),
        value
            .insertion_pole_order(vk_symbol!("gammalooprs::ε"))
            .unwrap()
    );
    seconds
}

fn plain_engine() -> Vakint {
    Vakint::new().unwrap()
}

/// The engine with the 1-4 loop package: REQUIRED (panics without it).
fn package_engine() -> Vakint {
    let directory = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4")
        .expect("five-loop records need VAKINT_TEST_P_LIBRARY_1_4");
    let mut vakint = Vakint::new().unwrap();
    let library = rustred_app::propagator::load_p_library(
        rustred_app::propagator::PLibrarySource::Directory(directory.into()),
        &rustred_app::propagator::PLibraryOptions {
            verify_descent: false,
        },
    )
    .unwrap();
    vakint.inject_p_library(library).unwrap();
    vakint
}

/// The root of the beta's fixtures: `VAKINT_TEST_WH_FIXTURES`, REQUIRED.
fn fixtures_root() -> PathBuf {
    let root = PathBuf::from(
        std::env::var_os("VAKINT_TEST_WH_FIXTURES")
            .expect("set VAKINT_TEST_WH_FIXTURES to the beta's wh-fixtures directory"),
    );
    assert!(root.is_dir(), "{}: not a directory", root.display());
    root
}

/// The record directories (holding input.txt and value.txt) under `dir`,
/// sorted; FAILS on a missing directory.
fn record_dirs(dir: &Path) -> Vec<PathBuf> {
    assert!(dir.is_dir(), "{}: missing", dir.display());
    let mut out = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.join("input.txt").is_file() && p.join("value.txt").is_file())
        .collect::<Vec<_>>();
    out.sort();
    out
}

/// Every record of `dirs` (or only those with insertions), compared; the
/// total number of records and the number compared are asserted.
fn compare_set(
    vakint: &Vakint,
    dirs: &[PathBuf],
    total: usize,
    with_insertions: Option<usize>,
    what: &str,
) {
    assert_eq!(
        dirs.len(),
        total,
        "{what}: {} records, expected {total}",
        dirs.len()
    );
    let mut compared = 0;
    for dir in dirs {
        let input = std::fs::read_to_string(dir.join("input.txt")).unwrap();
        let value_text = std::fs::read_to_string(dir.join("value.txt")).unwrap();
        let value = FactorizedNumerator::parse_text(&input)
            .unwrap_or_else(|e| panic!("{}: {e}", dir.display()));
        if with_insertions.is_some() && value.insertions.is_empty() {
            continue;
        }
        compare(vakint, &input, &value_text, &dir.display().to_string());
        compared += 1;
    }
    assert_eq!(
        compared,
        with_insertions.unwrap_or(total),
        "{what}: records compared"
    );
    eprintln!("INSERTION_RECORD_SET {what} total={total} compared={compared}");
}

#[test]
fn the_beta_k2_child_with_an_insertion_equals_its_contracted_value() {
    compare(
        &plain_engine(),
        include_str!("fixtures/wh-k2-child-ins.record.txt"),
        include_str!("fixtures/wh-k2-child-ins.value.txt"),
        "k2-child-ins",
    );
}

#[test]
fn the_beta_k4_quotient_equals_its_contracted_value() {
    compare(
        &plain_engine(),
        include_str!("fixtures/wh-k4-quotient.record.txt"),
        include_str!("fixtures/wh-k4-quotient.value.txt"),
        "k4-quotient",
    );
}

#[test]
fn the_beta_k5_surrogate_root_equals_its_contracted_value() {
    compare(
        &plain_engine(),
        include_str!("fixtures/wh-k5-surrogate.record.txt"),
        include_str!("fixtures/wh-k5-surrogate.value.txt"),
        "k5-surrogate",
    );
}

/// k5-stress (L5 graph 3793): 5 records, 4 insertion vertices each; r0051 has
/// five loops (the package is required).
#[test]
#[ignore = "needs VAKINT_TEST_WH_FIXTURES and VAKINT_TEST_P_LIBRARY_1_4"]
fn the_beta_k5_stress_records_equal_their_values() {
    let dirs = record_dirs(&fixtures_root().join("k5-stress"));
    compare_set(&package_engine(), &dirs, 5, None, "k5-stress");
}

/// k1-quick (L4 graphs 40 / 42 / 49 / 100): 35 records, the 30 with child
/// insertions compared here (the 5 bare roots are lane CLOSURE's test).
#[test]
#[ignore = "needs VAKINT_TEST_WH_FIXTURES"]
fn the_beta_k1_quick_records_with_insertions_equal_their_values() {
    let root = fixtures_root().join("k1-quick");
    let mut dirs = Vec::new();
    for graph in ["L4-40", "L4-42", "L4-49", "L4-100"] {
        dirs.extend(record_dirs(&root.join(graph)));
    }
    compare_set(&plain_engine(), &dirs, 35, Some(30), "k1-quick");
}

/// k1-l5/L5-60 (FIVE loops): 15 records, the 14 with child insertions.
#[test]
#[ignore = "needs VAKINT_TEST_WH_FIXTURES and VAKINT_TEST_P_LIBRARY_1_4"]
fn the_beta_k1_l5_60_records_with_insertions_equal_their_values() {
    let dirs = record_dirs(&fixtures_root().join("k1-l5").join("L5-60"));
    compare_set(&package_engine(), &dirs, 15, Some(14), "L5-60");
}

/// k1-l5/L5-51 (FIVE loops): 6 records, the 5 with child insertions.
#[test]
#[ignore = "needs VAKINT_TEST_WH_FIXTURES and VAKINT_TEST_P_LIBRARY_1_4"]
fn the_beta_k1_l5_51_records_with_insertions_equal_their_values() {
    let dirs = record_dirs(&fixtures_root().join("k1-l5").join("L5-51"));
    compare_set(&package_engine(), &dirs, 6, Some(5), "L5-51");
}

/// The comparator itself: a nonzero difference at a NEGATIVE order is
/// reported (not a hang), and a zero reference fails.
#[test]
fn the_comparator_reports_a_negative_order_and_refuses_a_zero_reference() {
    crate::Vakint::initialize_vakint_symbols();
    let name = "gammalooprs::ε";
    let eps = vk_parse!(name).unwrap();
    let a = Atom::num(3) / eps.pow(Atom::num(2)) + Atom::num(1) / &eps;
    let b = Atom::num(3) / eps.pow(Atom::num(2)) + Atom::num(2) / &eps;
    assert_eq!(assert_equal_through(&a, &a, -1, "same"), 2);
    let differs = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        assert_equal_through(&a, &b, -1, "differs")
    }));
    assert!(differs.is_err());
    let zero = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        assert_equal_through(&Atom::Zero, &Atom::Zero, -1, "zero")
    }));
    assert!(zero.is_err());
    // beyond the order nothing is compared
    let c = &a + &eps;
    assert_eq!(assert_equal_through(&c, &a, -1, "beyond"), 2);
}
