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
/// The insertion at `vertex` moved onto its vertex by hand WITH ITS POLE
/// TAKEN OUT (a vertex factor may carry no pole in ε): the vertex factor
/// becomes coefficient × ε^n × factor, n the insertion's own pole order;
/// returns the moved value and n (the value with the insertion is ε^-n times
/// the moved one, exactly, in Retain).
fn insertion_moved_onto_its_vertex(
    value: &FactorizedNumerator,
    vertex: usize,
) -> (FactorizedNumerator, i64) {
    let eps_text = "gammalooprs::ε";
    let eps = vk_symbol!(eps_text);
    let mut moved = value.clone();
    let i = moved
        .insertions
        .iter()
        .position(|x| x.vertex == vertex)
        .unwrap();
    let insertion = moved.insertions.remove(i);
    let mut alone = value.clone();
    alone.insertions = vec![insertion.clone()];
    let n = i64::try_from(alone.insertion_pole_order(eps).unwrap()).unwrap();
    let v = moved.vertices.iter_mut().find(|x| x.id == vertex).unwrap();
    v.factor =
        (&insertion.coefficient * Atom::var(eps).pow(Atom::num(n)) * &insertion.factor).expand();
    v.rule = format!("V_moved{vertex}");
    (moved, n)
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
/// vertex by hand with its pole ε^-1 taken out gives the value times ε, exactly
/// (a vertex factor carries no pole); (3) the value does not depend on the
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
    let (moved, n) = insertion_moved_onto_its_vertex(&k2, 4);
    assert_eq!(n, 1);
    let moved = closure_value(&vakint, &retain, &moved);
    let eps = Atom::var(vk_symbol!("gammalooprs::ε"));
    assert!(
        (&exact - &(&moved * eps.pow(Atom::num(-n))))
            .expand()
            .is_zero(),
        "exact product in Retain:\ninsertion {exact}\nmoved     {moved} (times ε^-{n})"
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
/// insertion moved onto its vertex with its own pole taken out gives the
/// retained value times ε^n, exactly, and the value does not depend on the
/// depth.
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
    let (moved, n) = insertion_moved_onto_its_vertex(&k5, 1);
    let moved = closure_value(&vakint, &retain, &moved);
    let eps = Atom::var(vk_symbol!("gammalooprs::ε"));
    assert!(
        (&exact - &(&moved * eps.pow(Atom::num(-n))))
            .expand()
            .is_zero(),
        "exact product in Retain (two insertions, the moved one of pole order {n})"
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

/// An insertion UNDEFINED after D -> 4-2ε is refused; its two defined forms
/// agree. The gate's witness: the K2 child of L3-036 with its coefficient's
/// 1/ε replaced by 1/(D-4+2ε) (identically zero denominator; the entry
/// returned 0 before), refused (Insertion at vertex 4, "a denominator
/// vanishes identically after D -> 4-2ε") by validation and by the entry. The
/// two defined forms of the same pole, 1/ε and -2/(D-4) (equal after the
/// substitution), are accepted and give EQUAL values. The same with an
/// epsilon setting named user::ep (the coefficient written with it): the
/// undefined form is refused at the entry (validation recognises ε by name),
/// the two defined ones agree.
#[test]
fn an_undefined_insertion_is_refused_and_its_defined_twins_agree() {
    let vakint = plain_engine();
    let text = include_str!("fixtures/wh-k2-child-ins.record.txt");
    let mut user = site_settings();
    user.epsilon_symbol = "user::ep".into();
    let prefix = "a denominator vanishes identically after D -> 4-2ε";
    for (eps, settings) in [("gammalooprs::ε", site_settings()), ("user::ep", user)] {
        // the record with its coefficient's 1/ε (gammalooprs::ε) replaced by `pole`
        let with = |pole: &str| {
            let mut value = FactorizedNumerator::from_text(text).unwrap();
            let original_eps = "gammalooprs::ε";
            let replaced = format!("{original_eps}*{pole}");
            let factor = vk_parse!(replaced.as_str()).unwrap();
            value.insertions[0].coefficient = (&value.insertions[0].coefficient * &factor).expand();
            value
        };
        let undefined = with(&format!("1/(gammalooprs::{{}}::dim-4+2*{eps})"));
        if eps == "gammalooprs::ε" {
            match undefined.validate().unwrap_err() {
                FactorizedError::Insertion { vertex: 4, detail } => {
                    assert!(detail.starts_with(prefix), "{detail}")
                }
                other => panic!("expected the Insertion refusal, got {other:?}"),
            }
        }
        match evaluate_factorized_closure(
            &vakint,
            &settings,
            std::slice::from_ref(&undefined),
            &[undefined.prefactor.total()],
        )
        .unwrap_err()
        {
            FactorizedError::Insertion { vertex: 4, detail } => {
                assert!(detail.starts_with(prefix), "{eps}: {detail}")
            }
            other => panic!("{eps}: expected the Insertion refusal, got {other:?}"),
        }
        let one_over_eps = closure_value(&vakint, &settings, &with(&format!("1/{eps}")));
        let in_dimension = closure_value(&vakint, &settings, &with("-2/(gammalooprs::{}::dim-4)"));
        assert!(!one_over_eps.is_zero(), "{eps}");
        assert!(
            (&one_over_eps - &in_dimension).expand().is_zero(),
            "{eps}: 1/eps and -2/(D-4) must agree\n1/eps    {one_over_eps}\n-2/(D-4) {in_dimension}"
        );
    }
}

/// Every denominator of an insertion is tested for being identically zero,
/// whatever it depends on: with P(t) = (t-3)^2 - t^2 + 6t - 9 (identically 0
/// as a polynomial) on the K2 child of L3-036, the coefficient's 1/ε replaced
/// by 1/P(D), by 1/P(x) (x a caller scalar, no ε or dimension in it), the
/// same next to an explicit ε^-1, and P(x)^(-1/2): each refused (Insertion at
/// vertex 4, "a denominator vanishes identically"), by validation and by the
/// entry. The defined controls 1/(ε (P(t) + 1)) (P + 1 = 1) for t = D and
/// t = x give the original record's value exactly; 1/(x - 3), zero only for
/// one value of x, is defined and evaluates.
#[test]
fn every_identically_vanishing_denominator_is_refused() {
    let vakint = plain_engine();
    let text = include_str!("fixtures/wh-k2-child-ins.record.txt");
    let original = FactorizedNumerator::from_text(text).unwrap();
    let reference = closure_value(&vakint, &site_settings(), &original);
    // the record with its coefficient's 1/ε replaced by `replacement`
    let with = |replacement: &str| {
        let mut value = original.clone();
        // the coefficient's own factor c ε expanded (ε-free), the replacement
        // multiplied in as written (expanding it would turn P(t) into 0)
        let eps_text = "gammalooprs::ε";
        let eps = vk_parse!(eps_text).unwrap();
        let written = vk_parse!(replacement).unwrap();
        value.insertions[0].coefficient =
            (&value.insertions[0].coefficient * &eps).expand() * &written;
        value
    };
    let p = |t: &str| format!("(({t}-3)^2-{t}^2+6*{t}-9)");
    let d = "gammalooprs::{}::dim";
    let x = "witness::x";
    let prefix = "a denominator vanishes identically";
    for (what, replacement) in [
        ("1/P(D)", format!("1/{}", p(d))),
        ("1/P(x)", format!("1/{}", p(x))),
        (
            "1/P(D) next to an explicit ε^-1",
            format!("gammalooprs::ε^-1/{}", p(d)),
        ),
        (
            "1/P(x) next to an explicit ε^-1",
            format!("gammalooprs::ε^-1/{}", p(x)),
        ),
        ("P(x)^(-1/2)", format!("{}^(-1/2)", p(x))),
    ] {
        let value = with(&replacement);
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
                FactorizedError::Insertion { vertex: 4, detail } => {
                    assert!(detail.starts_with(prefix), "{what}: {detail}")
                }
                other => panic!("{what}: expected the Insertion refusal, got {other:?}"),
            }
        }
    }
    for t in [d, x] {
        let control = with(&format!("1/(gammalooprs::ε*({}+1))", p(t)));
        let got = closure_value(&vakint, &site_settings(), &control);
        assert!(
            (&got - &reference).expand().is_zero(),
            "1/(ε (P({t})+1)) must give the original value"
        );
    }
    let at_a_point = with(&format!("1/({x}-3)"));
    at_a_point.validate().expect("1/(x-3) is defined");
    assert!(!closure_value(&vakint, &site_settings(), &at_a_point).is_zero());
}

/// A vanishing base is refused under every exponent that is not a positive
/// rational number, and every exponent that is not an integer literal is
/// refused whatever its base: on the K2 child of L3-036, with P(t)
/// identically zero, for e = -sqrt(2), -π, P(y)-1 (identically -1, not a
/// literal), P(y) (0^0) and P(y)+2, ε^-1·P(x)^e is refused by the zero test
/// ("a denominator vanishes identically"), and its twin ε^-1·(P(x)+1)^e,
/// which the zero test passes, by the integer-exponent rule ("a non-integer
/// exponent"): `Insertion` at vertex 4, at validation and at the entry.
#[test]
fn a_vanishing_base_or_an_exponent_not_an_integer_literal_is_refused() {
    let vakint = plain_engine();
    let original =
        FactorizedNumerator::from_text(include_str!("fixtures/wh-k2-child-ins.record.txt"))
            .unwrap();
    let eps_text = "gammalooprs::ε";
    let eps = vk_parse!(eps_text).unwrap();
    // the record with its coefficient's 1/ε replaced by `replacement`, as written
    let with = |replacement: &str| {
        let mut value = original.clone();
        let written = vk_parse!(replacement).unwrap();
        value.insertions[0].coefficient =
            (&value.insertions[0].coefficient * &eps).expand() * &written;
        value
    };
    let p = |t: &str| format!("(({t}-3)^2-{t}^2+6*{t}-9)");
    let (x, y) = ("witness::x", "witness::y");
    for exponent in [
        "-2^(1/2)".to_string(),
        "-𝜋".to_string(),
        format!("{}-1", p(y)),
        p(y),
        format!("{}+2", p(y)),
    ] {
        for (base, prefix) in [
            (p(x), "a denominator vanishes identically"),
            (format!("({}+1)", p(x)), "a non-integer exponent"),
        ] {
            let value = with(&format!("gammalooprs::ε^-1*{base}^({exponent})"));
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
                    FactorizedError::Insertion { vertex: 4, detail } => {
                        assert!(detail.starts_with(prefix), "{base}^({exponent}): {detail}")
                    }
                    other => {
                        panic!("{base}^({exponent}): expected the Insertion refusal, got {other:?}")
                    }
                }
            }
        }
    }
}

/// EVERY CALLER FIELD outside the insertions is defined and finite at ε = 0
/// (D -> 4-2ε): on the K2 child of L3-036 (closure kind) and the real
/// direct-kind L3-185 record, each of an edge factor, a vertex factor, the
/// projector, prefactor_overall, prefactor_scalars and prefactor_colour
/// multiplied by 1/P(x) or 1/P(D) (P(t) = (t-3)^2-t^2+6t-9, identically
/// zero, written unexpanded) is refused by validation with "a denominator
/// vanishes identically" at that field (`NotScalar`), and multiplied by
/// 1/(D-4) with "a pole in ε outside an insertion". The finite forms pass:
/// (D^2-16)/(D-4) (the pole cancels in lowest terms; the normaliser keeps
/// the quotient) and the producer's prefactor bases 1/(D-1), 1/(D+2), 1/D.
#[test]
fn every_caller_field_is_defined_and_finite_at_epsilon_zero() {
    let closure =
        FactorizedNumerator::from_text(include_str!("fixtures/wh-k2-child-ins.record.txt"))
            .unwrap();
    let direct =
        FactorizedNumerator::from_text(include_str!("fixtures/l3-185-insertions.txt")).unwrap();
    let d = "gammalooprs::{}::dim";
    let p = |t: &str| format!("(({t}-3)^2-{t}^2+6*{t}-9)");
    let pa = |text: &str| vk_parse!(text).unwrap();
    // the first vertex without an insertion (an insertion's vertex factor is 1)
    let plain_vertex = |value: &FactorizedNumerator| {
        value
            .vertices
            .iter()
            .position(|v| value.insertions.iter().all(|i| i.vertex != v.id))
            .unwrap()
    };
    // the field, multiplied by `factor` as written (not expanded)
    let multiplied = |original: &FactorizedNumerator, field: &str, factor: &str| {
        let mut value = original.clone();
        let factor = pa(factor);
        let times = |atom: &mut Atom| *atom = &*atom * &factor;
        let vertex = plain_vertex(original);
        match field {
            "edge" => times(&mut value.edges[0].factor),
            "vertex" => times(&mut value.vertices[vertex].factor),
            "projector" => times(&mut value.projector.as_mut().unwrap().tensor),
            // (records without a projector skip this field below)
            "prefactor_overall" => times(&mut value.prefactor.overall),
            "prefactor_scalars" => times(&mut value.prefactor.scalars),
            "prefactor_colour" => times(&mut value.prefactor.colour),
            other => panic!("{other}"),
        }
        value
    };
    let location_of = |value: &FactorizedNumerator, field: &str| match field {
        "edge" => format!("edge {} factor", value.edges[0].id),
        "vertex" => format!("vertex {} factor", value.vertices[plain_vertex(value)].id),
        other => other.to_string(),
    };
    let fields = [
        "edge",
        "vertex",
        "projector",
        "prefactor_overall",
        "prefactor_scalars",
        "prefactor_colour",
    ];
    for original in [&closure, &direct] {
        original.validate().unwrap();
        for field in fields {
            // no projector; the closure kind's prefactor_overall is fixed
            // to 1 by its own rule (the original graph's weights)
            if (field == "projector" && original.projector.is_none())
                || (field == "prefactor_overall"
                    && matches!(original.kind, FactorizedKind::Closure(_)))
            {
                continue;
            }
            for (factor, prefix) in [
                (
                    format!("1/{}", p("witness::x")),
                    "a denominator vanishes identically",
                ),
                (format!("1/{}", p(d)), "a denominator vanishes identically"),
                (format!("1/({d}-4)"), "a pole in ε outside an insertion"),
            ] {
                let value = multiplied(original, field, &factor);
                match value.validate() {
                    Err(FactorizedError::NotScalar {
                        location: at,
                        detail,
                    }) => {
                        assert_eq!(at, location_of(original, field), "{factor} in {field}");
                        assert!(detail.starts_with(prefix), "{factor} in {field}: {detail}");
                    }
                    other => panic!(
                        "{}: {factor} in {field}: expected NotScalar, got {other:?}",
                        original.label
                    ),
                }
            }
            let finite = multiplied(original, field, &format!("({d}^2-16)/({d}-4)"));
            finite
                .validate()
                .unwrap_or_else(|e| panic!("(D^2-16)/(D-4) in {field}: {e}"));
        }
        for twin in [format!("1/({d}-1)"), format!("1/({d}+2)"), format!("1/{d}")] {
            multiplied(original, "prefactor_scalars", &twin)
                .validate()
                .unwrap_or_else(|e| panic!("{twin} in prefactor_scalars: {e}"));
        }
    }
    // through the closure entry: refused before any evaluation, and a twin
    // is evaluated
    let vakint = plain_engine();
    let bad = multiplied(&closure, "vertex", &format!("1/({d}-4)"));
    assert!(matches!(
        evaluate_factorized_closure(
            &vakint,
            &site_settings(),
            std::slice::from_ref(&bad),
            &[bad.prefactor.total()],
        ),
        Err(FactorizedError::NotScalar { .. })
    ));
    let twin = multiplied(&closure, "prefactor_scalars", &format!("1/({d}-1)"));
    assert!(!closure_value(&vakint, &site_settings(), &twin).is_zero());
    // a record WITHOUT insertions (the real k1-root): validation cannot know
    // its ε; the entry checks the fields with the settings' ε
    let root =
        FactorizedNumerator::from_text(include_str!("fixtures/wh-k1-root.record.txt")).unwrap();
    let poled = multiplied(&root, "vertex", "gammalooprs::ε^-1");
    poled.validate().unwrap();
    match evaluate_factorized_closure(
        &vakint,
        &site_settings(),
        std::slice::from_ref(&poled),
        &[poled.prefactor.total()],
    ) {
        Err(FactorizedError::NotScalar { location, detail }) => {
            assert_eq!(location, location_of(&root, "vertex"));
            assert!(
                detail.starts_with("a pole in ε outside an insertion"),
                "{detail}"
            );
        }
        other => panic!("ε^-1 in a vertex factor of k1-root: expected NotScalar, got {other:?}"),
    }
}

/// The real direct-kind L3-185 record carries `p2_power = -1` on its
/// projector: it validates, and the marked entry checks (prefactor argument,
/// ε setting, depth) accept it; the field checks act on the projector's
/// tensor, never on its power.
#[test]
fn a_real_record_with_p2_power_minus_one_is_accepted() {
    let value =
        FactorizedNumerator::from_text(include_str!("fixtures/l3-185-insertions.txt")).unwrap();
    assert_eq!(value.projector.as_ref().unwrap().p2_power, -1);
    value.validate().unwrap();
    let eps_text = "gammalooprs::ε";
    check_marked_entry(
        &value,
        &value.prefactor.total(),
        vk_symbol!(eps_text),
        Some(3),
    )
    .unwrap();
}

/// A v3 direct record's `depth` sets its number of terms (format section 7);
/// the direct, ports and vacuum entries take the settings' number: on the
/// real L3-185 record with depth 2, settings of 3 terms are refused ("the
/// depth 2 of the record differs from the settings' number of terms 3"), at
/// the entry check and through the direct entry; 2 terms pass the check, and
/// a record without a depth takes any number.
#[test]
fn a_record_depth_other_than_the_settings_terms_is_refused() {
    let mut value =
        FactorizedNumerator::from_text(include_str!("fixtures/l3-185-insertions.txt")).unwrap();
    let eps_text = "gammalooprs::ε";
    let eps = vk_symbol!(eps_text);
    check_marked_entry(&value, &value.prefactor.total(), eps, Some(3)).unwrap();
    value
        .record
        .get_or_insert_with(|| DirectRecord {
            original_cycles: None,
            depth: None,
            weights: Vec::new(),
        })
        .depth = Some(2);
    value.validate().unwrap();
    let detail = "the depth 2 of the record differs from the settings' number of terms 3";
    match check_marked_entry(&value, &value.prefactor.total(), eps, Some(3)) {
        Err(FactorizedError::Graph { detail: got }) => assert!(got.starts_with(detail), "{got}"),
        other => panic!("expected Graph, got {other:?}"),
    }
    check_marked_entry(&value, &value.prefactor.total(), eps, Some(2)).unwrap();
    match evaluate_factorized(
        &plain_engine(),
        &site_settings(),
        &value,
        &value.prefactor.total(),
    ) {
        Err(FactorizedError::Graph { detail: got }) => assert!(got.starts_with(detail), "{got}"),
        other => panic!("expected Graph through the direct entry, got {other:?}"),
    }
}
