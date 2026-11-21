//! The beta's REAL R* records WITH INSERTIONS through the closure entry,
//! against the value their contracted path returns at that call.
//!
//! Through the SHARED fixture helper `tests/support/wh_fixtures.rs` (load or
//! fail, signed comparison through valid_through with a nonzero reference,
//! declared zeros, record counts). Rules: a missing fixture
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
use std::path::PathBuf;

#[path = "../../../tests/support/wh_fixtures.rs"]
mod wh_fixtures;

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
    let epsilon = vk_symbol!("gammalooprs::ε");
    let order = wh_fixtures::valid_through(value_text);
    let orders = if wh_fixtures::is_declared_zero(value_text) {
        wh_fixtures::assert_declared_zero(&got.value, epsilon, order, what)
    } else {
        crate::Vakint::initialize_vakint_symbols();
        let line = wh_fixtures::value_line(value_text);
        let expected = vk_parse!(line).unwrap();
        wh_fixtures::assert_equal_through(&got.value, &expected, epsilon, order, what)
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
    let directory = wh_fixtures::required_env(
        "VAKINT_TEST_P_LIBRARY_1_4",
        "five-loop records need the 1-4 loop package",
    );
    let mut vakint = Vakint::new().unwrap();
    let library = rustred_app::propagator::load_p_library(
        rustred_app::propagator::PLibrarySource::Directory(directory),
        &rustred_app::propagator::PLibraryOptions {
            verify_descent: false,
        },
    )
    .unwrap();
    vakint.inject_p_library(library).unwrap();
    vakint
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
        let input = wh_fixtures::read_fixture(&dir.join("input.txt"));
        let value_text = wh_fixtures::read_fixture(&dir.join("value.txt"));
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

/// The k5-surrogate record with insertions of the given coefficients: an
/// existing insertion's coefficient replaced, any other vertex turned into a
/// contracted child (ct:) with its own factor as the insertion's factor.
fn with_coefficients(coefficients: &[(usize, String)]) -> FactorizedNumerator {
    let mut value =
        FactorizedNumerator::from_text(include_str!("fixtures/wh-k5-surrogate.record.txt"))
            .unwrap();
    for (vertex, text) in coefficients {
        let coefficient = vk_parse!(text.as_str()).unwrap();
        if let Some(i) = value.insertions.iter().position(|x| x.vertex == *vertex) {
            value.insertions[i].coefficient = coefficient;
            continue;
        }
        let hedges = value
            .local_hedges(Location::Vertex(*vertex))
            .into_iter()
            .collect::<Vec<_>>();
        let v = value.vertices.iter_mut().find(|v| v.id == *vertex).unwrap();
        let factor = std::mem::replace(&mut v.factor, Atom::num(1));
        v.rule = format!("ct:child{vertex}");
        value.insertions.push(Insertion {
            vertex: *vertex,
            hedges,
            cycles: (0, 0),
            coefficient,
            factor,
        });
    }
    value.insertions.sort_by_key(|x| x.vertex);
    value
}

/// The pole order's arithmetic cannot wrap or panic: insertions whose powers
/// of ε add beyond the integers (two ε^-(2^63-1) and two ε^(2^63-1) on the
/// record's ε-free vertices with ε^-1 on its first insertion; three
/// ε^-(2^63-1); one, whose pole order does not fit the depth) are refused at
/// validation and by the closure entry, `Insertion` "an ε power out of range"
/// at the first such vertex (5). Without the bound the first two panicked on
/// the sum (and would wrap without overflow checks) and the third ran without
/// end on a negative depth. The defect-free twin (ε^-1 on the first insertion,
/// coefficient 1 elsewhere: pole order 3) equals ε^-1 times the record with
/// coefficient 1 evaluated one term deeper, truncated at ε^(depth - L - 1).
#[test]
fn the_pole_order_arithmetic_cannot_wrap() {
    let eps = "gammalooprs::ε";
    let big = "9223372036854775807";
    let vakint = plain_engine();
    for (what, coefficients) in [
        (
            "five",
            vec![
                (1, format!("{eps}^-1")),
                (5, format!("{eps}^-{big}")),
                (6, format!("{eps}^-{big}")),
                (7, format!("{eps}^{big}")),
                (9, format!("{eps}^{big}")),
            ],
        ),
        (
            "three",
            vec![
                (5, format!("{eps}^-{big}")),
                (6, format!("{eps}^-{big}")),
                (7, format!("{eps}^-{big}")),
            ],
        ),
        ("one", vec![(5, format!("{eps}^-{big}"))]),
    ] {
        let value = with_coefficients(&coefficients);
        for error in [
            value.validate().unwrap_err(),
            evaluate_factorized_closure(
                &vakint,
                &site_settings(),
                std::slice::from_ref(&value),
                &[value.prefactor.total()],
            )
            .unwrap_err(),
        ] {
            match &error {
                FactorizedError::Insertion { vertex: 5, detail } => {
                    assert!(
                        detail.starts_with("an ε power out of range"),
                        "{what}: {detail}"
                    )
                }
                other => {
                    panic!("{what}: expected the Insertion refusal at vertex 5, got {other:?}")
                }
            }
        }
    }
    let ones = [5, 6, 7, 9].map(|v| (v, "1".to_string()));
    let mut twin_coefficients = vec![(1, format!("{eps}^-1"))];
    twin_coefficients.extend(ones.iter().cloned());
    let twin = with_coefficients(&twin_coefficients);
    twin.validate().unwrap();
    let epsilon = vk_symbol!("gammalooprs::ε");
    assert_eq!(twin.insertion_pole_order(epsilon).unwrap(), 3);
    let got = evaluate_factorized_closure(
        &vakint,
        &site_settings(),
        std::slice::from_ref(&twin),
        &[twin.prefactor.total()],
    )
    .unwrap();
    let mut all_ones = vec![(1, "1".to_string())];
    all_ones.extend(ones.iter().cloned());
    let mut deeper = with_coefficients(&all_ones);
    let FactorizedKind::Closure(spec) = &mut deeper.kind else {
        unreachable!()
    };
    let depth = spec.depth.unwrap();
    spec.depth = Some(depth + 1);
    assert_eq!(deeper.insertion_pole_order(epsilon).unwrap(), 2);
    let bare = evaluate_factorized_closure(
        &vakint,
        &site_settings(),
        std::slice::from_ref(&deeper),
        &[deeper.prefactor.total()],
    )
    .unwrap();
    let target = depth as i64 - twin.loops() as i64 - 1;
    let inverse_text = format!("{eps}^-1");
    let inverse = vk_parse!(inverse_text.as_str()).unwrap();
    let reference = (inverse * &bare.value)
        .series(
            epsilon,
            Atom::Zero.as_view(),
            symbolica::domains::rational::Rational::from(target),
        )
        .unwrap()
        .to_atom()
        .expand();
    assert!(!reference.is_zero());
    assert!(
        (&got.value - &reference).expand().is_zero(),
        "twin {}\nreference {}",
        got.value,
        reference
    );
}

/// k5-stress (L5 graph 3793): 5 records, 4 insertion vertices each; r0051 has
/// five loops (the package is required).
#[test]
#[ignore = "needs VAKINT_TEST_WH_FIXTURES and VAKINT_TEST_P_LIBRARY_1_4"]
fn the_beta_k5_stress_records_equal_their_values() {
    let dirs = wh_fixtures::record_dirs(&wh_fixtures::fixture_dir("k5-stress"));
    compare_set(&package_engine(), &dirs, 5, None, "k5-stress");
}

/// k1-quick (L4 graphs 40 / 42 / 49 / 100): 35 records, the 30 with child
/// insertions compared here (the 5 bare roots are lane CLOSURE's test).
#[test]
#[ignore = "needs VAKINT_TEST_WH_FIXTURES"]
fn the_beta_k1_quick_records_with_insertions_equal_their_values() {
    let dirs = wh_fixtures::record_dirs(&wh_fixtures::fixture_dir("k1-quick"));
    compare_set(&plain_engine(), &dirs, 35, Some(30), "k1-quick");
}

/// k1-l5/L5-60 (FIVE loops): 15 records, the 14 with child insertions.
#[test]
#[ignore = "needs VAKINT_TEST_WH_FIXTURES and VAKINT_TEST_P_LIBRARY_1_4"]
fn the_beta_k1_l5_60_records_with_insertions_equal_their_values() {
    let dirs = wh_fixtures::record_dirs(&wh_fixtures::fixture_dir("k1-l5/L5-60"));
    compare_set(&package_engine(), &dirs, 15, Some(14), "L5-60");
}

/// k1-l5/L5-51 (FIVE loops): 6 records, the 5 with child insertions.
#[test]
#[ignore = "needs VAKINT_TEST_WH_FIXTURES and VAKINT_TEST_P_LIBRARY_1_4"]
fn the_beta_k1_l5_51_records_with_insertions_equal_their_values() {
    let dirs = wh_fixtures::record_dirs(&wh_fixtures::fixture_dir("k1-l5/L5-51"));
    compare_set(&package_engine(), &dirs, 6, Some(5), "L5-51");
}
