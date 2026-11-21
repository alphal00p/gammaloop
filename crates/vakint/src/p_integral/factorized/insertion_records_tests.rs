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
/// coefficient 1 evaluated one term deeper, truncated at ε^(depth - L - 1): a
/// CONSISTENCY check of the depth raise (both sides through the entry); the
/// independent checks are the closure tests'
/// `dimension_dependent_insertion_factors_equal_the_independent_reference`
/// and `a_cancelling_coefficient_raises_the_depth_by_its_exact_pole_order`.
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

/// The TOTAL depth raise is capped at 64 over all insertions of a record:
/// five insertions on the k5 surrogate record (its own factors add 2): ε^-13
/// on vertices 5, 6, 7, 9 and ε^-10 on 1 add up to exactly 64, accepted by
/// the depth raise (checked there, `insertion_depth_raise`, not evaluated 64
/// terms deeper); ε^-11 on 1 (65) and ε^-13 on all five (67) are refused, by
/// the depth raise and by the closure entry before any reduction (`Insertion`
/// at the first insertion's vertex, the total above 64). Each insertion alone
/// is within the per-insertion count.
#[test]
fn the_total_depth_raise_is_capped() {
    let eps = "gammalooprs::ε";
    let epsilon = vk_symbol!("gammalooprs::ε");
    let five = |first: i64, others: i64| {
        let mut coefficients = vec![(1, format!("{eps}^-{first}"))];
        coefficients.extend([5, 6, 7, 9].map(|v| (v, format!("{eps}^-{others}"))));
        with_coefficients(&coefficients)
    };
    let at_cap = five(10, 13);
    at_cap.validate().unwrap();
    assert_eq!(at_cap.insertion_depth_raise(epsilon).unwrap(), 64);
    let vakint = plain_engine();
    for (what, value, total) in [("65", five(11, 13), 65), ("67", five(13, 13), 67)] {
        value
            .validate()
            .unwrap_or_else(|e| panic!("{what}: each insertion is within the count: {e}"));
        let prefix = format!(
            "an ε power out of range: the insertions' pole orders add up to {total}, above 64"
        );
        for error in [
            value.insertion_depth_raise(epsilon).unwrap_err(),
            evaluate_factorized_closure(
                &vakint,
                &site_settings(),
                std::slice::from_ref(&value),
                &[value.prefactor.total()],
            )
            .unwrap_err(),
        ] {
            match &error {
                FactorizedError::Insertion { vertex: 1, detail } => {
                    assert!(detail.starts_with(&prefix), "{what}: {detail}")
                }
                other => {
                    panic!("{what}: expected the Insertion refusal at vertex 1, got {other:?}")
                }
            }
        }
    }
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

/// The depth of a marked value is a number of ε terms between 1 and
/// MARKED_MAX_TERMS (128) with the closure entry's insertion raise. The
/// gate's witness, the L2-22 K1 root (vendored k1-root, depth 2) with only
/// its depth set to usize::MAX, returned a zero double pole (the depth was
/// cast to -1); now it is refused (Closure, "the depth .. is out of range"),
/// by validation (text and typed records alike) and by the entry. Refused
/// too: i64::MAX and i64::MAX + 1 as record depths, 0, and 129 (checked at
/// validation, never evaluated); the settings' number of terms (record depth
/// none) -1, 0, i64::MIN (which overflowed in the truncation) and 129; with
/// an insertion of pole order 1 (k2-child-ins), depth 128 (raised to 129).
/// Accepted at validation: depth 128 (the cap). The ordinary control keeps
/// both poles: k1-root at depth 2 equals its value.
#[test]
fn an_out_of_range_depth_is_refused() {
    let vakint = plain_engine();
    let root = include_str!("fixtures/wh-k1-root.record.txt");
    let k2 = include_str!("fixtures/wh-k2-child-ins.record.txt");
    let with_depth = |text: &str, depth: Option<usize>| {
        let mut value = FactorizedNumerator::from_text(text).unwrap();
        let FactorizedKind::Closure(spec) = &mut value.kind else {
            unreachable!()
        };
        spec.depth = depth;
        value
    };
    let refused = |error: FactorizedError, depth: &str, what: &str| match &error {
        FactorizedError::Closure { detail } => assert!(
            detail.starts_with(&format!("the depth {depth} is out of range")),
            "{what}: {detail}"
        ),
        other => panic!("{what}: expected the Closure depth refusal, got {other:?}"),
    };
    let entry = |value: &FactorizedNumerator, settings: &VakintSettings| {
        evaluate_factorized_closure(
            &vakint,
            settings,
            std::slice::from_ref(value),
            &[value.prefactor.total()],
        )
    };
    let huge = usize::try_from(i64::MAX).unwrap();
    for depth in [usize::MAX, huge, huge + 1, 0, 129] {
        let value = with_depth(root, Some(depth));
        refused(
            value.validate().unwrap_err(),
            &depth.to_string(),
            "record depth",
        );
        refused(
            entry(&value, &site_settings()).unwrap_err(),
            &depth.to_string(),
            "record depth through the entry",
        );
        // the text form: the vendored text with only its depth line changed,
        // refused by the reader (which validates)
        let text = root.replace("\ndepth\t2\n", &format!("\ndepth\t{depth}\n"));
        assert_ne!(text, root);
        refused(
            FactorizedNumerator::from_text(&text).unwrap_err(),
            &depth.to_string(),
            "record depth from text",
        );
    }
    with_depth(root, Some(128))
        .validate()
        .expect("depth 128 is at the cap");
    for terms in [-1i64, 0, i64::MIN, 129] {
        let mut settings = site_settings();
        settings.number_of_terms_in_epsilon_expansion = terms;
        refused(
            entry(&with_depth(root, None), &settings).unwrap_err(),
            &terms.to_string(),
            "the settings' number of terms",
        );
    }
    refused(
        entry(&with_depth(k2, Some(128)), &site_settings()).unwrap_err(),
        "128",
        "depth 128 raised by the insertion's pole order 1",
    );
    compare(
        &vakint,
        root,
        include_str!("fixtures/wh-k1-root.value.txt"),
        "k1-root at depth 2 (control)",
    );
}

/// The R* site settings with RustRed master mode Retain (masters kept as
/// symbols: the value keeps its exact ε dependence, format section 7).
fn retain_settings() -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        ..site_settings()
    }
}

/// The record with its contracted (ct:) insertion at `vertex` moved onto the
/// vertex by hand: the vertex factor is coefficient x factor, no insertion.
fn insertion_moved_onto_its_vertex(
    value: &FactorizedNumerator,
    vertex: usize,
) -> FactorizedNumerator {
    let mut moved = value.clone();
    let i = moved
        .insertions
        .iter()
        .position(|x| x.vertex == vertex)
        .unwrap();
    let insertion = moved.insertions.remove(i);
    let v = moved.vertices.iter_mut().find(|x| x.id == vertex).unwrap();
    v.factor = (&insertion.coefficient * &insertion.factor).expand();
    v.rule = format!("V_moved{vertex}");
    moved
}

fn closure_value(vakint: &Vakint, settings: &VakintSettings, value: &FactorizedNumerator) -> Atom {
    evaluate_factorized_closure(
        vakint,
        settings,
        std::slice::from_ref(value),
        &[value.prefactor.total()],
    )
    .unwrap_or_else(|e| panic!("{}: {e}", value.label))
    .value
}

/// Retain keeps the EXACT ε dependence (format section 7) also with a poled
/// insertion. The K2 child of L3-036 (vendored k2-child-ins: one ct:
/// insertion at vertex 4, coefficient ∝ ε^-1, pole order 1, depth 1):
/// (1) linearity: its coefficient times 1/(1 + z ε), z a fresh scalar, gives
/// the retained value times 1/(1 + z ε) exactly (a finite series would be a
/// polynomial in z); (2) the exact product: the insertion moved onto its
/// vertex by hand gives the same value; (3) the value does not depend on the
/// depth. Symbolic control: the same record keeps its vendored value.
#[test]
fn retain_keeps_the_exact_product_of_a_poled_insertion() {
    let vakint = plain_engine();
    let retain = retain_settings();
    let k2 = FactorizedNumerator::from_text(include_str!("fixtures/wh-k2-child-ins.record.txt"))
        .unwrap();
    assert_eq!(
        k2.insertion_pole_order(vk_symbol!("gammalooprs::ε"))
            .unwrap(),
        1
    );
    let exact = closure_value(&vakint, &retain, &k2);
    let rational_text = "1/(1+witness::z*gammalooprs::ε)";
    let rational = vk_parse!(rational_text).unwrap();
    let mut twin = k2.clone();
    twin.insertions[0].coefficient = &twin.insertions[0].coefficient * &rational;
    let got = closure_value(&vakint, &retain, &twin);
    assert!(
        (&got - &(&exact * &rational)).expand().is_zero(),
        "linearity in Retain:\ntwin  {got}\nexact x 1/(1+z eps)  {}",
        &exact * &rational
    );
    let moved = closure_value(&vakint, &retain, &insertion_moved_onto_its_vertex(&k2, 4));
    assert!(
        (&exact - &moved).expand().is_zero(),
        "exact product in Retain:\ninsertion {exact}\nmoved     {moved}"
    );
    let mut deeper = k2.clone();
    if let FactorizedKind::Closure(spec) = &mut deeper.kind {
        spec.depth = Some(4);
    }
    assert!(
        (&exact - &closure_value(&vakint, &retain, &deeper))
            .expand()
            .is_zero(),
        "Retain must not depend on the depth"
    );
    compare(
        &vakint,
        include_str!("fixtures/wh-k2-child-ins.record.txt"),
        include_str!("fixtures/wh-k2-child-ins.value.txt"),
        "k2-child-ins (Symbolic control)",
    );
}

/// The same exact identities with TWO insertions (vendored k5-surrogate: a
/// ct: child at vertex 1 and a surrogate ins: vertex 8, pole order 2): the ct:
/// insertion moved onto its vertex gives the same retained value, and the
/// value does not depend on the depth.
#[test]
fn retain_keeps_the_exact_product_with_two_insertions() {
    let vakint = plain_engine();
    let retain = retain_settings();
    let k5 = FactorizedNumerator::from_text(include_str!("fixtures/wh-k5-surrogate.record.txt"))
        .unwrap();
    assert_eq!(
        k5.insertion_pole_order(vk_symbol!("gammalooprs::ε"))
            .unwrap(),
        2
    );
    let exact = closure_value(&vakint, &retain, &k5);
    let moved = closure_value(&vakint, &retain, &insertion_moved_onto_its_vertex(&k5, 1));
    assert!(
        (&exact - &moved).expand().is_zero(),
        "exact product in Retain (two insertions)"
    );
    let mut deeper = k5.clone();
    if let FactorizedKind::Closure(spec) = &mut deeper.kind {
        spec.depth = spec.depth.map(|d| d + 2);
    }
    assert!(
        (&exact - &closure_value(&vakint, &retain, &deeper))
            .expand()
            .is_zero(),
        "Retain must not depend on the depth (two insertions)"
    );
}
