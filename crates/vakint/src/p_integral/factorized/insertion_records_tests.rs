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
                // the exponents beyond 64 are refused by the atom rules'
                // bound before any count (the counts cannot wrap)
                FactorizedError::Insertion { vertex: 5, detail } => {
                    assert!(
                        detail.starts_with("a non-integer exponent in insertion 5")
                            && detail.contains("exceeds 64 in modulus"),
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
/// insertions compared here (the 5 bare roots are the closure tests').
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
/// MARKED_MAX_TERMS (32: far above the producer's 12; its largest depth +
/// raise is 5 + 4) with the closure entry's insertion raise. The
/// witness, the L2-22 K1 root (vendored k1-root, depth 2) with only
/// its depth set to usize::MAX, returned a zero double pole (the depth was
/// cast to -1); now it is refused (Closure, "the depth .. is out of range"),
/// by validation (text and typed records alike) and by the entry. Refused
/// too: i64::MAX and i64::MAX + 1 as record depths, 0, and 33 (checked at
/// validation, never evaluated); the settings' number of terms (record depth
/// none) -1, 0, i64::MIN (which overflowed in the truncation) and 33; with
/// an insertion of pole order 1 (k2-child-ins), depth 32 (raised to 33).
/// Accepted at validation: depth 32 (the cap). The ordinary control keeps
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
    for depth in [usize::MAX, huge, huge + 1, 0, 33] {
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
    with_depth(root, Some(32))
        .validate()
        .expect("depth 32 is at the cap");
    for terms in [-1i64, 0, i64::MIN, 33] {
        let mut settings = site_settings();
        settings.number_of_terms_in_epsilon_expansion = terms;
        refused(
            entry(&with_depth(root, None), &settings).unwrap_err(),
            &terms.to_string(),
            "the settings' number of terms",
        );
    }
    refused(
        entry(&with_depth(k2, Some(32)), &site_settings()).unwrap_err(),
        "32",
        "depth 32 raised by the insertion's pole order 1",
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
/// agree. The witness: the K2 child of L3-036 with its coefficient's
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

/// The data atom MZV is accepted exactly where vakint's evaluator is: on the
/// real K2 child of L3-036, its insertion coefficient times MZV(1) or
/// MZV(1, 2) (a FIRST index 1 diverges) is refused by validation and by the
/// entry (it is no data atom: the reserved-name refusal naming vakint::MZV and
/// the field), typed and through the text form, in Symbolic and Retain mode;
/// times MZV(2), MZV(2, 1) (a later 1 converges) or MZV(5, 3) it is accepted,
/// survives the text round trip, and the value is the original one times that
/// MZV exactly (both poles; in Retain the exact value times it).
#[test]
fn an_mzv_is_accepted_exactly_when_it_converges() {
    let vakint = plain_engine();
    let original =
        FactorizedNumerator::from_text(include_str!("fixtures/wh-k2-child-ins.record.txt"))
            .unwrap();
    let reference = closure_value(&vakint, &site_settings(), &original);
    let with = |mzv: &str| {
        let mut value = original.clone();
        value.insertions[0].coefficient =
            &value.insertions[0].coefficient * vk_parse!(mzv).unwrap();
        value
    };
    for bad in ["vakint::MZV(1)", "vakint::MZV(1,2)"] {
        let value = with(bad);
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
                FactorizedError::Graph { detail } | FactorizedError::Insertion { detail, .. } => {
                    assert!(
                        detail.contains("vakint::MZV")
                            && detail.contains("insertion 4 coefficient"),
                        "{bad}: {detail}"
                    )
                }
                other => panic!("{bad}: expected a refusal naming the atom, got {other:?}"),
            }
        }
    }
    // the text form: written and read back, the same refusal
    for bad in ["vakint::MZV(1)", "vakint::MZV(1,2)"] {
        let refused = match with(bad).to_text() {
            Err(_) => true,
            Ok(text) => FactorizedNumerator::from_text(&text)
                .is_err_and(|e| e.to_string().contains("vakint::MZV")),
        };
        assert!(refused, "{bad} (text): accepted");
    }
    let retain = retain_settings();
    let retained = closure_value(&vakint, &retain, &original);
    for good in ["vakint::MZV(2)", "vakint::MZV(2,1)", "vakint::MZV(5,3)"] {
        let mzv = vk_parse!(good).unwrap();
        let got = closure_value(&vakint, &site_settings(), &with(good));
        assert!(
            (&got - &(&reference * &mzv)).expand().is_zero(),
            "{good}: not linear"
        );
        // Retain: the exact value times the MZV
        let got = closure_value(&vakint, &retain, &with(good));
        assert!(
            (&got - &(&retained * &mzv)).expand().is_zero(),
            "{good} (Retain): not linear"
        );
        // the text form carries it unchanged
        let read = FactorizedNumerator::from_text(&with(good).to_text().unwrap()).unwrap();
        assert_eq!(read, with(good), "{good}: text round trip");
    }
    // Retain refuses the divergent ones too
    for bad in ["vakint::MZV(1)", "vakint::MZV(1,2)"] {
        let value = with(bad);
        assert!(
            evaluate_factorized_closure(
                &vakint,
                &retain,
                std::slice::from_ref(&value),
                &[value.prefactor.total()]
            )
            .is_err(),
            "{bad} (Retain): accepted"
        );
    }
}

/// SIGNED EXTREMES: every caller integer the validator bounds is compared
/// without a signed absolute value (|i64::MIN| does not exist): on the real
/// K2 child (insertion coefficient) and the real direct L3-185 record,
/// x^(i64::MIN), UNKNOWN(master, i64::MIN), MZV(i64::MIN), an index label
/// hedge(i64::MIN, 1) and the compact slash slash(i64::MIN, 1) in a slot of
/// a contracted momentum pair, an edge power and a momentum coefficient of
/// i64::MIN, and each with i64::MIN + 1 and i64::MAX, are refused with the
/// rule's text, never a panic.
#[test]
fn signed_extremes_are_refused_not_a_panic() {
    let k2 = FactorizedNumerator::from_text(include_str!("fixtures/wh-k2-child-ins.record.txt"))
        .unwrap();
    let direct =
        FactorizedNumerator::from_text(include_str!("fixtures/l3-185-insertions.txt")).unwrap();
    let d = "gammalooprs::{}::dim";
    for n in [i64::MIN, i64::MIN + 1, i64::MAX] {
        for (what, coefficient, text) in [
            ("x^n", format!("witness::x^({n})"), "exceeds 64 in modulus"),
            (
                "UNKNOWN(master, n)",
                format!("vakint::UNKNOWN(vakint::ggppss_int_1,{n})"),
                "vakint::UNKNOWN",
            ),
            ("MZV(n)", format!("vakint::MZV({n})"), "vakint::MZV"),
        ] {
            let mut v = k2.clone();
            v.insertions[0].coefficient =
                &v.insertions[0].coefficient * vk_parse!(&coefficient).unwrap();
            match v.validate() {
                Err(e) => assert!(e.to_string().contains(text), "{what} at {n}: {e}"),
                Ok(()) => panic!("{what} at {n}: accepted"),
            }
        }
        // an index label and the compact slash in a slot of a vertex factor
        // (a caller field: the shape rules)
        let plain = direct
            .vertices
            .iter()
            .position(|x| direct.insertions.iter().all(|i| i.vertex != x.id))
            .unwrap();
        for (what, label, text) in [
            (
                "hedge(n, 1)",
                format!("gammalooprs::hedge({n},1)"),
                "is not Id",
            ),
            (
                "slash(n, 1)",
                format!("vakint::slash({n},1)"),
                "is not SignedId",
            ),
        ] {
            let mut v = direct.clone();
            // a contracted pair of momenta (spenso contracts g(a, a) to D at
            // once, which would remove the label)
            let slot = format!("spenso::mink({d},{label})");
            v.vertices[plain].factor = &v.vertices[plain].factor
                * vk_parse!(&format!(
                    "gammalooprs::Q(8,{slot})*gammalooprs::Q(9,{slot})"
                ))
                .unwrap();
            match v.validate() {
                Err(e) => assert!(e.to_string().contains(text), "{what} at {n}: {e}"),
                Ok(()) => panic!("{what} at {n}: accepted"),
            }
        }
        let mut v = direct.clone();
        v.edges[0].power = n;
        match v.validate() {
            Err(e) => assert!(e.to_string().contains("power"), "edge power {n}: {e}"),
            Ok(()) => panic!("edge power {n}: accepted"),
        }
        let mut v = direct.clone();
        v.edges[0].momentum[0] = n;
        match v.validate() {
            Err(e) => assert!(
                e.to_string().contains("momentum coefficient"),
                "momentum {n}: {e}"
            ),
            Ok(()) => panic!("momentum coefficient {n}: accepted"),
        }
    }
    // the accepted endpoints of the exponent and the UNKNOWN order: +-64
    for n in [64, -64] {
        for coefficient in [
            format!("witness::x^({n})"),
            format!("vakint::UNKNOWN(vakint::ggppss_int_1,{n})"),
        ] {
            let mut v = k2.clone();
            v.insertions[0].coefficient =
                &v.insertions[0].coefficient * vk_parse!(&coefficient).unwrap();
            v.validate()
                .unwrap_or_else(|e| panic!("{coefficient}: {e}"));
        }
    }
}

/// The BOUNDS at the closure, direct and ports entries and at the public
/// input functions `vakint_input` and `contract_to_vakint_dots`, typed and
/// through the text form (the vacuum entry and its input functions:
/// vacuum::tests): the real K2 child (closure) and the real direct L3-185
/// record (the others) with a depth of
/// usize::MAX, an edge or a vertex id of i64::MAX + 1, an edge power or a
/// momentum coefficient of i64::MIN, 1/P(x) in a weight or a vertex factor:
/// each refused (no panic) with the bound's text.
#[test]
fn the_bounds_reach_every_entry_typed_and_through_text() {
    let vakint = plain_engine();
    let k2 = FactorizedNumerator::from_text(include_str!("fixtures/wh-k2-child-ins.record.txt"))
        .unwrap();
    let direct =
        FactorizedNumerator::from_text(include_str!("fixtures/l3-185-insertions.txt")).unwrap();
    let big = i64::MAX as usize + 1;
    let zero = "1/((witness::x-3)^2-witness::x^2+6*witness::x-9)";
    type Change = Box<dyn Fn(&mut FactorizedNumerator)>;
    let cases = || -> Vec<(&'static str, Change, &'static str)> {
        vec![
            (
                "depth usize::MAX",
                Box::new(move |v: &mut FactorizedNumerator| {
                    if let FactorizedKind::Closure(spec) = &mut v.kind {
                        spec.depth = Some(usize::MAX);
                    } else {
                        v.record
                            .get_or_insert_with(|| DirectRecord {
                                original_cycles: None,
                                depth: None,
                                weights: Vec::new(),
                            })
                            .depth = Some(usize::MAX);
                    }
                }),
                "depth",
            ),
            (
                "edge id i64::MAX + 1",
                Box::new(move |v: &mut FactorizedNumerator| {
                    let last = v.edges.len() - 1;
                    let old = v.edges[last].id;
                    v.edges[last].id = big;
                    for e in v.loop_edges.iter_mut().chain(v.external_edges.iter_mut()) {
                        if *e == Some(old) {
                            *e = Some(big);
                        }
                    }
                    if let FactorizedKind::Closure(spec) = &mut v.kind
                        && spec.edge == old
                    {
                        spec.edge = big;
                    }
                }),
                "ids",
            ),
            (
                "vertex id i64::MAX + 1",
                Box::new(move |v: &mut FactorizedNumerator| {
                    let last = v.vertices.len() - 1;
                    let old = v.vertices[last].id;
                    v.vertices[last].id = big;
                    for e in &mut v.edges {
                        for end in [&mut e.source, &mut e.target] {
                            if *end == Endpoint::Vertex(old) {
                                *end = Endpoint::Vertex(big);
                            }
                        }
                    }
                    for i in &mut v.insertions {
                        if i.vertex == old {
                            i.vertex = big;
                        }
                    }
                }),
                "ids",
            ),
            (
                "edge power i64::MIN",
                Box::new(|v: &mut FactorizedNumerator| v.edges[0].power = i64::MIN),
                "power",
            ),
            (
                "momentum coefficient i64::MIN",
                Box::new(|v: &mut FactorizedNumerator| v.edges[0].momentum[0] = i64::MIN),
                "momentum coefficient",
            ),
            (
                "1/P(x) in a weight",
                Box::new(move |v: &mut FactorizedNumerator| {
                    let weights = match &mut v.kind {
                        FactorizedKind::Closure(spec) => &mut spec.weights,
                        _ => {
                            &mut v
                                .record
                                .get_or_insert_with(|| DirectRecord {
                                    original_cycles: None,
                                    depth: None,
                                    weights: Vec::new(),
                                })
                                .weights
                        }
                    };
                    let bad = vk_parse!(zero).unwrap();
                    match weights.iter_mut().find(|(name, _, _)| name == "overall") {
                        Some((_, _, atom)) => *atom = &*atom * &bad,
                        None => weights.push(("overall".into(), false, bad)),
                    }
                }),
                "a denominator vanishes identically",
            ),
        ]
    };
    let settings = site_settings();
    // the public input functions and builders are entries of their own
    let entries: Vec<(&str, &FactorizedNumerator)> = vec![
        ("closure", &k2),
        ("direct", &direct),
        ("ports", &direct),
        ("vakint_input", &direct),
        ("contract_to_vakint_dots", &direct),
    ];
    for (entry, base) in entries {
        base.validate().unwrap();
        for (what, change, text) in cases() {
            let mut v = base.clone();
            change(&mut v);
            let run = |v: &FactorizedNumerator| -> std::result::Result<(), String> {
                match entry {
                    "closure" => evaluate_factorized_closure(
                        &vakint,
                        &settings,
                        std::slice::from_ref(v),
                        &[v.prefactor.total()],
                    )
                    .map(|_| ()),
                    "direct" => {
                        evaluate_factorized(&vakint, &settings, v, &v.prefactor.total()).map(|_| ())
                    }
                    "vakint_input" => vakint_input(&settings, v, &v.prefactor.total()).map(|_| ()),
                    "contract_to_vakint_dots" => contract_to_vakint_dots(v).map(|_| ()),
                    _ => evaluate_factorized_ports(&vakint, &settings, v, &v.prefactor.total())
                        .map(|_| ()),
                }
                .map_err(|e| e.to_string())
            };
            match run(&v) {
                Err(e) => assert!(e.contains(text), "{entry} {what} (typed): {e}"),
                Ok(()) => panic!("{entry} {what} (typed): accepted"),
            }
            let refused = match v.to_text() {
                Err(_) => true,
                Ok(written) => match FactorizedNumerator::from_text(&written) {
                    Err(_) => true,
                    Ok(read) => run(&read).is_err(),
                },
            };
            assert!(refused, "{entry} {what} (text): accepted");
        }
    }
}

/// EVERY CALLER-WRITTEN ATOM ([`FactorizedNumerator::record_atoms`], built
/// without `..`) outside the insertions gets the defined-and-finite check and
/// the grammar: on the real K2 child (closure kind, given the producer's
/// weights section: overall, loop_sign, colour; its scalars weight, applied by
/// vakint, must equal prefactor_scalars by its own rule) and the real direct
/// L3-185, each atom multiplied by 1/P(x) (P identically zero) and,
/// separately, by the caller function witness::f(witness::x) is refused by
/// validation (no panic; a zero atom takes the marker itself); where the
/// field's own structural rule does not refuse first (masses, ports, the
/// closure line, the closure kind's fixed overall, an insertion vertex's
/// factor 1), with "a denominator vanishes identically" and "a function head
/// that is not allowed", the field named. A WEIGHT multiplied by the unit
/// form P(x) + 1 is accepted and the value is unchanged (both poles).
#[test]
fn every_caller_written_atom_gets_the_finite_check_and_the_grammar() {
    let k2 = FactorizedNumerator::from_text(include_str!("fixtures/wh-k2-child-ins.record.txt"))
        .unwrap();
    let direct =
        FactorizedNumerator::from_text(include_str!("fixtures/l3-185-insertions.txt")).unwrap();
    let zero = parsed("1/((witness::x-3)^2-witness::x^2+6*witness::x-9)");
    let head = parsed("witness::f(witness::x)");
    let structural = |field: &str, closure: bool, base: &FactorizedNumerator| {
        base.insertions
            .iter()
            .any(|i| field == format!("vertex {} factor", i.vertex))
            || field.ends_with(" mass")
            || field == "ports"
            || field == "closure mass"
            || field == "taylor auxiliary mass"
            || (closure && field == "prefactor_overall")
    };
    // a weights section (overall 1, loop_sign -1: their prefactor
    // counterparts are 1 on this record)
    let mut k2 = k2;
    if let FactorizedKind::Closure(spec) = &mut k2.kind {
        spec.weights = vec![
            ("loop_sign".into(), false, parsed("-1")),
            ("overall".into(), false, parsed("1")),
        ];
    }
    let mut weights_seen = 0;
    for base in [&k2, &direct] {
        base.validate().unwrap();
        let closure = matches!(base.kind, FactorizedKind::Closure(_));
        let fields = base
            .record_atoms()
            .into_iter()
            .map(|(f, _)| f)
            .filter(|f| !f.starts_with("insertion "))
            .collect::<Vec<_>>();
        for (i, field) in fields.iter().enumerate() {
            weights_seen += usize::from(field.starts_with("weight "));
            for (marker, text) in [
                (&zero, "a denominator vanishes identically"),
                (&head, "a function head that is not allowed"),
            ] {
                let mut v = base.clone();
                {
                    let mut atoms = v
                        .record_atoms_mut()
                        .into_iter()
                        .filter(|(f, _)| !f.starts_with("insertion "))
                        .collect::<Vec<_>>();
                    let (_, atom) = &mut atoms[i];
                    // a zero field (a massless line's mass) takes the marker
                    // itself: a product with 0 is 0
                    **atom = if atom.is_zero() {
                        marker.clone()
                    } else {
                        &**atom * marker
                    };
                }
                let error = v
                    .validate()
                    .err()
                    .unwrap_or_else(|| panic!("{}: {field} x {marker}: accepted", base.label))
                    .to_string();
                if !structural(field, closure, base) {
                    assert!(
                        error.contains(text) && error.contains(field.as_str()),
                        "{}: {field} x {marker}: {error}",
                        base.label
                    );
                }
            }
        }
    }
    assert!(weights_seen >= 1, "no weight field was exercised");
    // a unit weight keeps the value
    let vakint = plain_engine();
    let reference = closure_value(&vakint, &site_settings(), &k2);
    let mut unit = k2.clone();
    if let FactorizedKind::Closure(spec) = &mut unit.kind {
        let (_, _, atom) = spec
            .weights
            .iter_mut()
            .find(|(n, _, _)| n == "overall")
            .unwrap();
        *atom = &*atom * parsed("((witness::x-3)^2-witness::x^2+6*witness::x-9)+1");
    }
    unit.validate().unwrap();
    let got = closure_value(&vakint, &site_settings(), &unit);
    assert!(
        (&got - &reference).expand().is_zero(),
        "a unit weight changed the value"
    );
}

/// The SHAPE rules run before any rule that computes, at every entry: a
/// TYPED marked value with the vertex factor 1/((1+UFO::x)^1000000000 -
/// UFO::y) (an exponent beyond 64, whose exact expansion does not end) is
/// refused within seconds by the closure (K2), direct and ports (L3-185)
/// entries ("exceeds 64 in modulus"); the vacuum entry: vacuum::tests.
#[test]
fn the_shape_rules_run_before_any_computing_rule_at_every_entry() {
    let vakint = plain_engine();
    let settings = site_settings();
    let k2 = FactorizedNumerator::from_text(include_str!("fixtures/wh-k2-child-ins.record.txt"))
        .unwrap();
    let direct =
        FactorizedNumerator::from_text(include_str!("fixtures/l3-185-insertions.txt")).unwrap();
    let huge = parsed("1/((1+UFO::x)^1000000000-UFO::y)");
    let with = |base: &FactorizedNumerator| {
        let mut v = base.clone();
        let plain = v
            .vertices
            .iter()
            .position(|x| base.insertions.iter().all(|i| i.vertex != x.id))
            .unwrap();
        v.vertices[plain].factor = &v.vertices[plain].factor * &huge;
        v
    };
    let (k2, direct) = (with(&k2), with(&direct));
    let started = std::time::Instant::now();
    let results = [
        (
            "closure",
            evaluate_factorized_closure(
                &vakint,
                &settings,
                std::slice::from_ref(&k2),
                &[k2.prefactor.total()],
            )
            .map(|_| ()),
        ),
        (
            "direct",
            evaluate_factorized(&vakint, &settings, &direct, &direct.prefactor.total()).map(|_| ()),
        ),
        (
            "ports",
            evaluate_factorized_ports(&vakint, &settings, &direct, &direct.prefactor.total())
                .map(|_| ()),
        ),
    ];
    for (entry, result) in results {
        match result {
            Err(e) => assert!(
                e.to_string().contains("exceeds 64 in modulus"),
                "{entry}: {e}"
            ),
            Ok(()) => panic!("{entry}: accepted"),
        }
    }
    assert!(started.elapsed().as_secs() < 10, "{:?}", started.elapsed());
}

/// An UNMARKED text with an unknown schema keeps the base's refusal text,
/// byte for byte (dea282615: "text line {line}: unknown schema {found:?}
/// (this reader knows {FACTORIZED_SCHEMA})"); only a v3 schema line written
/// inexactly names both schemas.
#[test]
fn an_unknown_schema_keeps_the_base_text() {
    let error =
        FactorizedNumerator::from_text("schema\tvakint-factorized/2\nlabel\tx\n").unwrap_err();
    assert_eq!(
        error.to_string(),
        "text line 1: unknown schema \"vakint-factorized/2\" (this reader knows pint6.factorized-numerator.v1)"
    );
    let error =
        FactorizedNumerator::from_text("schema\tpint6.factorized-numerator.v3 \nlabel\tx\n")
            .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("pint6.factorized-numerator.v1 and pint6.factorized-numerator.v3"),
        "{error}"
    );
}

/// `text` parsed (vk_parse! takes no literal here: the crate's `Result` alias).
fn parsed(text: &str) -> Atom {
    vk_parse!(text).unwrap()
}

/// INCONSISTENT COMBINATIONS: a closure value that also carries a direct
/// record (typed only; the text form cannot write both) is refused at
/// validation ("a closure value carries a direct record"); its record's
/// weights are walked anyway (record_atoms lists both weight sources).
#[test]
fn a_closure_value_with_a_direct_record_is_refused() {
    let mut v = FactorizedNumerator::from_text(include_str!("fixtures/wh-k2-child-ins.record.txt"))
        .unwrap();
    v.record = Some(DirectRecord {
        original_cycles: None,
        depth: None,
        weights: vec![("extra".into(), false, parsed("gammalooprs::ε^-1"))],
    });
    match v.validate() {
        Err(FactorizedError::Closure { detail }) => {
            assert!(
                detail.starts_with("a closure value carries a direct record"),
                "{detail}"
            )
        }
        other => panic!("expected the Closure refusal, got {other:?}"),
    }
    assert!(
        v.record_atoms()
            .iter()
            .any(|(field, _)| field == "weight extra")
    );
    assert!(
        v.record_atoms_mut()
            .iter()
            .any(|(field, _)| field == "weight extra")
    );
}

/// The producer's and vakint's evaluation orders are not mixed: the R*
/// sites' settings (RustRed in Symbolic), vakint's default order (RustRed
/// and Forcer both Numerical, then methods without a master mode) and the
/// producer's mapping (rustred and forcer both Symbolic, AlphaLoop / MATAD /
/// FMFT without a mode) pass the master-mode rule.
#[test]
fn the_producers_and_the_default_orders_have_one_master_mode() {
    check_master_modes(&site_settings()).unwrap();
    check_master_modes(&VakintSettings::default()).unwrap();
    let producer = VakintSettings {
        evaluation_order: crate::EvaluationOrder(vec![
            crate::EvaluationMethod::AlphaLoop(Default::default()),
            crate::EvaluationMethod::MATAD(Default::default()),
            crate::EvaluationMethod::FMFT(Default::default()),
            crate::EvaluationMethod::RustRed(RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Symbolic,
            }),
            crate::EvaluationMethod::Forcer(crate::ForcerOptions {
                master_mode: RustRedMasterMode::Symbolic,
            }),
        ]),
        ..site_settings()
    };
    check_master_modes(&producer).unwrap();
}

/// THE WEIGHTS of every marked kind (the record format's OQ-17 rule): on the
/// real direct K5 record (L2-17, prefactor_scalars 1), `weight scalars
/// applied_by=vakint 2` is refused by validation, through the text form and
/// at the ports entry in Symbolic, Numerical and Retain mode and at the direct
/// entry ("weight scalars applied_by=vakint is 2, prefactor_scalars is 1",
/// `Graph`); a weight applied by the caller with its prefactor counterpart not
/// 1 (overall with prefactor_overall 2) is refused; the twins (scalars = the
/// prefactor_scalars, overall 1) are accepted and leave the value unchanged.
/// The closure kind (the K2 child) refuses the same mismatch (`Closure`).
#[test]
fn the_weights_agree_with_the_prefactor_on_every_kind() {
    let vakint = plain_engine();
    let base =
        FactorizedNumerator::from_text(include_str!("fixtures/wh-k5-direct.record.txt")).unwrap();
    let depth = base.record.as_ref().unwrap().depth.unwrap() as i64;
    let with_weights = |v: &FactorizedNumerator, weights: Vec<(String, bool, Atom)>| {
        let mut v = v.clone();
        v.record.as_mut().unwrap().weights = weights;
        v
    };
    let mode = |m: RustRedMasterMode| VakintSettings {
        evaluation_order: crate::EvaluationOrder(vec![crate::EvaluationMethod::RustRed(
            RustRedEvaluationOptions { master_mode: m },
        )]),
        number_of_terms_in_epsilon_expansion: depth,
        ..site_settings()
    };
    let bad = with_weights(&base, vec![("scalars".into(), true, parsed("2"))]);
    let text = "weight scalars applied_by=vakint is 2, prefactor_scalars is 1";
    let graph = |e: FactorizedError, what: &str| match e {
        FactorizedError::Graph { detail } => assert!(detail.starts_with(text), "{what}: {detail}"),
        other => panic!("{what}: expected Graph, got {other:?}"),
    };
    graph(bad.validate().unwrap_err(), "validate");
    let written = bad.to_text().unwrap();
    assert!(written.contains("applied_by=vakint\t2"), "{written}");
    graph(
        FactorizedNumerator::from_text(&written).unwrap_err(),
        "from_text",
    );
    for m in [
        RustRedMasterMode::Symbolic,
        RustRedMasterMode::Numerical,
        RustRedMasterMode::Retain,
    ] {
        graph(
            evaluate_factorized_ports(&vakint, &mode(m), &bad, &bad.prefactor.total()).unwrap_err(),
            &format!("ports {m:?}"),
        );
    }
    graph(
        evaluate_factorized(
            &vakint,
            &mode(RustRedMasterMode::Symbolic),
            &bad,
            &bad.prefactor.total(),
        )
        .unwrap_err(),
        "direct",
    );
    let mut caller = with_weights(&base, vec![("overall".into(), false, parsed("1"))]);
    caller.prefactor.overall = parsed("2");
    match caller.validate() {
        Err(e) => assert!(
            e.to_string().contains("weight overall applied_by=beta"),
            "{e}"
        ),
        Ok(()) => panic!("a caller weight with prefactor_overall 2: accepted"),
    }
    // the twins
    let symbolic = mode(RustRedMasterMode::Symbolic);
    let reference = evaluate_factorized_ports(&vakint, &symbolic, &base, &base.prefactor.total())
        .unwrap()
        .value;
    for twin in [
        with_weights(
            &base,
            vec![("scalars".into(), true, base.prefactor.scalars.clone())],
        ),
        with_weights(&base, vec![("overall".into(), false, parsed("1"))]),
    ] {
        twin.validate().unwrap();
        let got = evaluate_factorized_ports(&vakint, &symbolic, &twin, &twin.prefactor.total())
            .unwrap()
            .value;
        assert!(
            (&got - &reference).expand().is_zero(),
            "a consistent weight changed the value"
        );
    }
    // no flavour symbol in prefactor_scalars, whatever the kind
    let mut nf = base.clone();
    nf.prefactor.scalars = parsed("UFO::nf");
    match nf.validate() {
        Err(e) => assert!(
            e.to_string().contains("a flavour symbol in the record"),
            "{e}"
        ),
        Ok(()) => panic!("nf in a direct record's prefactor_scalars: accepted"),
    }
    // the closure kind: unchanged
    let mut k2 =
        FactorizedNumerator::from_text(include_str!("fixtures/wh-k2-child-ins.record.txt"))
            .unwrap();
    if let FactorizedKind::Closure(spec) = &mut k2.kind {
        spec.weights = vec![("scalars".into(), true, parsed("2"))];
    }
    match k2.validate() {
        Err(FactorizedError::Closure { detail }) => {
            assert!(
                detail.starts_with("weight scalars applied_by=vakint is 2"),
                "{detail}"
            )
        }
        other => panic!("closure: expected Closure, got {other:?}"),
    }
}

/// TEXT ROUND TRIP: a closure value that carries a direct record (refused by
/// validation) cannot be written either: to_text refuses it (the writer would
/// drop the record and its text would be accepted).
#[test]
fn a_closure_value_with_a_direct_record_is_not_written() {
    let mut v = FactorizedNumerator::from_text(include_str!("fixtures/wh-k2-child-ins.record.txt"))
        .unwrap();
    v.record = Some(DirectRecord {
        original_cycles: None,
        depth: None,
        weights: vec![("extra".into(), false, parsed("gammalooprs::ε^-1"))],
    });
    match v.to_text() {
        Err(FactorizedError::Closure { detail }) => {
            assert!(
                detail.starts_with("a closure value carries a direct record"),
                "{detail}"
            )
        }
        other => panic!("expected the Closure refusal, got {other:?}"),
    }
}

/// An unknown schema line keeps the base's text unless it is the v3 name with
/// whitespace around it: "…v31" and "x…v3" get the base's text exactly.
#[test]
fn only_an_inexact_v3_schema_line_names_both_schemas() {
    for schema in [
        "pint6.factorized-numerator.v31",
        "xpint6.factorized-numerator.v3",
    ] {
        let error =
            FactorizedNumerator::from_text(&format!("schema\t{schema}\nlabel\tx\n")).unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "text line 1: unknown schema \"{schema}\" (this reader knows pint6.factorized-numerator.v1)"
            )
        );
    }
}

/// The exponent bound runs before the insertion rules that compute: on the
/// real K2 child, ((1+x)^1000000000 - y)^(-1) times the insertion coefficient,
/// or as the insertion factor's multiplier, is refused within a second by
/// validation, through the text form and at the closure entry ("exceeds 64 in
/// modulus"; the zero test of its base raised an integer point to that power
/// exactly before, about two minutes).
#[test]
fn a_huge_power_in_an_insertion_is_refused_before_any_computation() {
    let vakint = plain_engine();
    let k2 = FactorizedNumerator::from_text(include_str!("fixtures/wh-k2-child-ins.record.txt"))
        .unwrap();
    let huge = parsed("((1+UFO::x)^1000000000-UFO::y)^(-1)");
    for part in ["coefficient", "factor"] {
        let mut v = k2.clone();
        let atom = if part == "coefficient" {
            &mut v.insertions[0].coefficient
        } else {
            &mut v.insertions[0].factor
        };
        *atom = &*atom * &huge;
        let started = std::time::Instant::now();
        let refusals = [
            v.validate().map_err(|e| e.to_string()),
            match v.to_text() {
                Ok(text) => FactorizedNumerator::from_text(&text)
                    .map(|_| ())
                    .map_err(|e| e.to_string()),
                Err(e) => Err(e.to_string()),
            },
            evaluate_factorized_closure(
                &vakint,
                &site_settings(),
                std::slice::from_ref(&v),
                &[v.prefactor.total()],
            )
            .map(|_| ())
            .map_err(|e| e.to_string()),
        ];
        for (i, r) in refusals.iter().enumerate() {
            match r {
                Err(e) => assert!(e.contains("exceeds 64 in modulus"), "{part} #{i}: {e}"),
                Ok(()) => panic!("{part} #{i}: accepted"),
            }
        }
        assert!(
            started.elapsed().as_secs_f64() < 3.0,
            "{part}: {:?}",
            started.elapsed()
        );
    }
}
