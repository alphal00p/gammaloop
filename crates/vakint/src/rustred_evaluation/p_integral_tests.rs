//! The RustRed p backend against Forcer's frozen one-loop battery, with the
//! injected reference library `rustred::propagator::testing::
//! one_loop_reference_library()` (the massless bubble behind the real engine):
//! exact agreement in the tag basis, agreement after the Symbolic and
//! Numerical finalize, the Wick sign applied once, channel splitting, the
//! settings' arithmetic, the step cache, the direct entry points, the timing
//! span and the dispatch without p data.

use std::collections::{BTreeMap, BTreeSet};

use rustred::propagator::library::{PBasis, PReduction, PTerm};
use rustred::propagator::tag::GammaTag;
use symbolica::atom::{Atom, AtomCore};

use super::*;
use crate::p_integral::compare::tag_basis_difference;
use crate::p_integral::forcer::fixtures::{ForcerCase, ForcerFixture, registry_family};
use crate::utils::vakint_macros::vk_parse;
use crate::{
    EvaluationMethod, EvaluationOrder, ForcerOptions, LoopNormalizationFactor, OrderedCacheMode,
    OrderedCacheSettings, RustRedEvaluationOptions,
};

const P1L_BATTERY: &str = include_str!("../../tests/fixtures/forcer/p1l_battery.json");

fn battery() -> ForcerFixture {
    let fixture = ForcerFixture::parse(P1L_BATTERY).unwrap();
    assert_eq!(fixture.cases.len(), 45);
    fixture
}

fn case<'a>(fixture: &'a ForcerFixture, name: &str) -> &'a ForcerCase {
    fixture.cases.iter().find(|case| case.name == name).unwrap()
}

fn parse(text: &str) -> Atom {
    Vakint::initialize_vakint_symbols();
    vk_parse!(text).unwrap()
}

/// An engine reducing p-integrals with the reference library.
fn engine() -> Vakint {
    let mut engine = Vakint::new().unwrap();
    engine
        .inject_p_library(rustred::propagator::testing::one_loop_reference_library().unwrap())
        .unwrap();
    engine
}

/// RustRed alone in `mode`, G-scheme, through `ε^(terms-2)` at one loop,
/// exact `arithmetic` or the default.
fn settings(mode: RustRedMasterMode, arithmetic: Option<&str>) -> VakintSettings {
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions { master_mode: mode },
        )]),
        form_exe_path: "/rustred-p-integrals-must-not-call-form".into(),
        integral_normalization_factor: LoopNormalizationFactor::GScheme,
        number_of_terms_in_epsilon_expansion: 4,
        run_time_decimal_precision: 40,
        use_dot_product_notation: true,
        n_cores: 1,
        ..VakintSettings::default()
    };
    if let Some(arithmetic) = arithmetic {
        settings.ordered_weight_arithmetic = arithmetic.parse().unwrap();
    }
    settings
}

fn epsilon(settings: &VakintSettings) -> Symbol {
    vk_symbol!(&settings.epsilon_symbol)
}

/// `reduced` per scale power, summed.
fn by_scale_power(values: impl IntoIterator<Item = ReducedPropagatorValue>) -> BTreeMap<i64, Atom> {
    let mut sums = BTreeMap::<i64, Atom>::new();
    for value in values {
        *sums.entry(value.scale_power).or_default() += value.reduced;
    }
    sums.retain(|_, reduced| !reduced.expand().is_zero());
    sums
}

/// The RustRed reduced values of one fixture case, from its exact backend
/// request (the input Forcer reduced).
fn rustred_values(
    engine: &Vakint,
    settings: &VakintSettings,
    case: &ForcerCase,
) -> BTreeMap<i64, Atom> {
    let numerator = parse(&case.numerator);
    let request = PropagatorRequest {
        family: registry_family(&case.family).unwrap(),
        powers: &case.powers,
        canonical_numerator: numerator.as_view(),
        external: parse(&case.external),
    };
    let values = reduce_requests(engine, settings, &[request], None)
        .unwrap_or_else(|error| panic!("{}: {error}", case.name));
    assert!(values.len() <= 1, "{}: one external head", case.name);
    by_scale_power(values.into_values().flatten())
}

fn frozen_values(settings: &VakintSettings, case: &ForcerCase) -> BTreeMap<i64, Atom> {
    by_scale_power(case.frozen(epsilon(settings)).unwrap())
}

/// `rustred` and `forcer` agree exactly in the tag basis, power by power.
fn assert_tag_basis_equal(
    name: &str,
    rustred: &BTreeMap<i64, Atom>,
    forcer: &BTreeMap<i64, Atom>,
    epsilon: Symbol,
) {
    assert_eq!(
        rustred.keys().collect::<Vec<_>>(),
        forcer.keys().collect::<Vec<_>>(),
        "{name}: scale powers"
    );
    for (w, value) in rustred {
        let difference = tag_basis_difference(value, &forcer[w], epsilon).unwrap();
        assert!(
            difference.is_empty(),
            "{name}, w = {w}: RustRed {value} vs Forcer {}: {difference:?}",
            forcer[w]
        );
    }
}

/// All 45 cases of the frozen one-loop battery: RustRed's reduced value is
/// Forcer's frozen output exactly in the (Γ(1+mε) monomial × master) basis,
/// under the default (modular) and the exact arithmetic.
#[test]
fn the_one_loop_battery_matches_forcer_exactly_in_the_tag_basis() {
    let engine = engine();
    let fixture = battery();
    for arithmetic in [None, Some("exact")] {
        let settings = settings(RustRedMasterMode::Retain, arithmetic);
        for case in &fixture.cases {
            let rustred = rustred_values(&engine, &settings, case);
            let forcer = frozen_values(&settings, case);
            assert_eq!(rustred.is_empty(), case.outputs.is_empty(), "{}", case.name);
            assert_tag_basis_equal(&case.name, &rustred, &forcer, epsilon(&settings));
        }
    }
}

/// The battery through `Vakint::evaluate` (matching, batching, the p backend
/// and the shared finalizer) equals the frozen Forcer outputs finalized like
/// a Forcer evaluation, in Symbolic mode exactly and in Numerical mode to
/// forty digits.
#[test]
fn the_one_loop_battery_agrees_after_symbolic_and_numerical_finalize() {
    let engine = engine();
    let fixture = battery();
    let symbolic = settings(RustRedMasterMode::Symbolic, None);
    let numerical = settings(RustRedMasterMode::Numerical, None);
    for case in &fixture.cases {
        let input = case.expression().unwrap();
        let actual = engine.evaluate(&symbolic, input.as_view()).unwrap();
        // A Forcer evaluation returns its finalized value through simplify_real.
        let expected = simplify_real(
            case.finalized(&symbolic, RustRedMasterMode::Symbolic)
                .unwrap()
                .as_view(),
        );
        assert!(
            (&actual - &expected).expand().is_zero(),
            "{}: RustRed {actual} vs Forcer {expected}",
            case.name
        );

        let actual = engine.evaluate(&numerical, input.as_view()).unwrap();
        let expected = engine
            .resolve_rustred_numerical(
                crate::rustred_evaluation::combine_exact(
                    case.finalized(&numerical, RustRedMasterMode::Numerical)
                        .unwrap(),
                ),
                &numerical,
            )
            .unwrap();
        assert_numerically_equal(&engine, &numerical, &actual, &expected, &case.name);
    }
}

/// `actual` and `expected` at `Q² = μ² = 1` agree to `1e-30`.
fn assert_numerically_equal(
    engine: &Vakint,
    settings: &VakintSettings,
    actual: &Atom,
    expected: &Atom,
    name: &str,
) {
    let unit = |value: &Atom| -> Atom {
        let x = Atom::var(crate::symbols::S.x_);
        let p = function!(crate::symbols::S.p, &x);
        value
            .replace(function!(crate::symbols::S.dot, &p, &p).to_pattern())
            .with(Atom::num(-1).to_pattern())
            .replace(function!(Symbol::LOG, &x).to_pattern())
            .with(Atom::Zero.to_pattern())
            .expand()
    };
    let evaluate = |value: &Atom| {
        engine
            .numerical_evaluation(
                settings,
                unit(value).as_view(),
                &Default::default(),
                &Default::default(),
                None,
            )
            .unwrap_or_else(|error| panic!("{name}: {value}: {error}"))
            .0
    };
    let (actual, expected) = (evaluate(actual), evaluate(expected));
    let (equal, diagnostic) = actual.does_approx_match(&expected, None, 1e-30, 5.0);
    assert!(equal, "{name}: {diagnostic}");
}

/// The Wick sign `(-1)^(Σa+r)` is applied once: for an odd `Σa` (the bubble
/// (2,1)) and an odd rank (`k·q` for `q = p(3)` and `q = -p(3)`), RustRed's
/// value equals Forcer's sign-carrying frozen value and differs from its
/// negative; the plain bubble is the anchor `1/ε` with the empty tag.
#[test]
fn the_wick_sign_is_applied_once() {
    let engine = engine();
    let fixture = battery();
    let settings = settings(RustRedMasterMode::Retain, None);
    let epsilon = epsilon(&settings);
    for name in [
        "bubble_21_n0",
        "bubble_11_n1",
        "bubble_21_n1",
        "bubble_12_n3",
    ] {
        let case = case(&fixture, name);
        let rustred = rustred_values(&engine, &settings, case);
        let forcer = frozen_values(&settings, case);
        assert_tag_basis_equal(name, &rustred, &forcer, epsilon);
        for (w, value) in &rustred {
            assert!(
                !tag_basis_difference(value, &(-&forcer[w]), epsilon)
                    .unwrap()
                    .is_empty(),
                "{name}: the sign must not cancel"
            );
        }
    }
    let anchor = rustred_values(&engine, &settings, case(&fixture, "bubble_11_n0"));
    assert_eq!(anchor.len(), 1);
    let one_over_epsilon = Atom::var(epsilon).pow(Atom::num(-1));
    assert!(
        tag_basis_difference(&anchor[&0], &one_over_epsilon, epsilon)
            .unwrap()
            .is_empty(),
        "P1L(1,1) = 1/ε: {}",
        anchor[&0]
    );
}

/// One expression with p(1) and p(3) bubbles and a vacuum tadpole: the
/// p-integrals form separate channels (one finalize per external head and a
/// request per scale power) next to the vacuum batch, and the sum evaluates
/// like its terms evaluated apart and like Forcer's frozen values.
#[test]
fn mixed_heads_and_vacuum_are_split_into_channels() {
    let engine = engine();
    let fixture = battery();
    let settings = settings(RustRedMasterMode::Symbolic, None);
    let p1 = case(&fixture, "bubble_22_n4"); // external -p(1), rank 4
    let p3 = case(&fixture, "bubble_31_n1"); // external p(3), rank 1
    let p3b = case(&fixture, "bubble_12_n2"); // external -p(3), rank 2
    let numerators = [p1, p3, p3b].map(|case| parse(&case.numerator));
    let requests = [p1, p3, p3b]
        .iter()
        .zip(&numerators)
        .map(|(case, numerator)| PropagatorRequest {
            family: registry_family(&case.family).unwrap(),
            powers: &case.powers,
            canonical_numerator: numerator.as_view(),
            external: parse(&case.external),
        })
        .collect::<Vec<_>>();
    let values = reduce_requests(&engine, &settings, &requests, None).unwrap();
    assert_eq!(
        values
            .keys()
            .cloned()
            .collect::<std::collections::BTreeSet<_>>(),
        [(1, parse("p(1)")), (1, parse("p(3)"))]
            .into_iter()
            .collect::<std::collections::BTreeSet<_>>()
    );
    // Rank 1 and rank 2 at Σa = 4 and 3: scale powers -1 and 1 for p(3).
    let mut p3_powers = values[&(1, parse("p(3)"))]
        .iter()
        .map(|value| value.scale_power)
        .collect::<Vec<_>>();
    p3_powers.sort();
    assert_eq!(p3_powers, [-1, 1]);

    let vacuum = "topo(prop(1,edge(1,1),k(1),muvsq,2))";
    let terms = [
        p1.expression().unwrap(),
        p3.expression().unwrap(),
        p3b.expression().unwrap(),
        parse(vacuum),
    ];
    let sum = terms.iter().fold(Atom::Zero, |sum, term| sum + term);
    let together = engine.evaluate(&settings, sum.as_view()).unwrap();
    let apart = terms
        .iter()
        .map(|term| engine.evaluate(&settings, term.as_view()).unwrap())
        .fold(Atom::Zero, |sum, value| sum + value);
    assert!(
        (&together - &apart).expand().is_zero(),
        "{together}\n{apart}"
    );
    let forcer = [p1, p3, p3b]
        .iter()
        .map(|case| {
            simplify_real(
                case.finalized(&settings, RustRedMasterMode::Symbolic)
                    .unwrap()
                    .as_view(),
            )
        })
        .fold(Atom::Zero, |sum, value| sum + value)
        + engine.evaluate(&settings, parse(vacuum).as_view()).unwrap();
    assert!((&together - &forcer).expand().is_zero());
}

/// The p backend takes the Ordered weight arithmetic from the settings
/// (default `modular_polynomial`, the `Z_p(d)` replay; user directive
/// 2026-11-07): the request `reduce_requests` sends for a channel
/// (`p_reduction_request`) carries it, with the settings' scheduler, uncovered
/// keys rejected and the reduction onto the masters; the exact arithmetic
/// gives the same result.
#[test]
fn the_settings_arithmetic_is_used_and_exact_agrees() {
    let modular = settings(RustRedMasterMode::Retain, None);
    let exact = settings(RustRedMasterMode::Retain, Some("exact"));
    for (settings, arithmetic) in [(&modular, "modular_polynomial"), (&exact, "exact")] {
        let request = p_reduction_request(&request_options(settings), 4, Vec::new(), None);
        assert_eq!(request.arithmetic.name(), arithmetic);
        assert_eq!(
            request.arithmetic,
            settings.rustred_ordered_weight_arithmetic()
        );
        assert_eq!(
            request.scheduler,
            crate::rustred_evaluation::ordered_scheduler(settings)
        );
        assert_eq!(request.uncovered, PUncoveredPolicy::Reject);
        assert!(request.onto_masters && request.step_cache.is_none());
    }
    let engine = engine();
    let fixture = battery();
    for name in ["bubble_32_n4", "bubble_13_n3", "bubble_22_n2"] {
        let input = case(&fixture, name).expression().unwrap();
        assert_eq!(
            engine.evaluate(&modular, input.as_view()).unwrap(),
            engine.evaluate(&exact, input.as_view()).unwrap(),
            "{name}"
        );
    }
}

/// Five-loop p requests take the five-loop arithmetic (M5 batch 5, notes
/// arith-5l): under the default selection `modular_sampled` on the same
/// workers, while one to four loops keep the settings' `modular_polynomial`;
/// an explicit selection (also an explicit `modular_polynomial`) is sent at
/// every loop count.
#[test]
fn five_loop_p_requests_default_to_modular_sampled() {
    use crate::OrderedWeightArithmetic;

    let mut default = VakintSettings {
        n_cores: 3,
        ..settings(RustRedMasterMode::Retain, None)
    };
    default.ordered_weight_arithmetic = OrderedWeightArithmetic::default();
    let options = request_options(&default);
    for loops in 1..=4 {
        let request = p_reduction_request(&options, loops, Vec::new(), None);
        assert_eq!(request.arithmetic.name(), "modular_polynomial", "{loops}L");
        assert_eq!(request.arithmetic, default.rustred_ordered_weight_arithmetic());
    }
    let request = p_reduction_request(&options, 5, Vec::new(), None);
    assert_eq!(
        (request.arithmetic.name(), request.arithmetic.threads()),
        ("modular_sampled", 3)
    );
    for explicit in [
        OrderedWeightArithmetic::modular_polynomial(),
        OrderedWeightArithmetic::exact(),
    ] {
        let settings = VakintSettings {
            ordered_weight_arithmetic: explicit,
            ..default.clone()
        };
        let options = request_options(&settings);
        for loops in 1..=5 {
            let request = p_reduction_request(&options, loops, Vec::new(), None);
            assert_eq!(
                request.arithmetic,
                settings.rustred_ordered_weight_arithmetic(),
                "{loops}L"
            );
        }
    }
}

/// Non-rational coefficients (`muvsq`, `ε`) key the channels: parts of one
/// scale power with different coefficients are reduced apart and each is
/// multiplied by its own coefficient, so the reduced values are Forcer's
/// frozen values times the coefficients, before and after the Symbolic
/// finalize (review round 1, S2; the glue binary checks the reviewer's input
/// against a live Forcer).
#[test]
fn non_rational_coefficients_key_the_channels() {
    let engine = engine();
    let fixture = battery();
    let settings = settings(RustRedMasterMode::Symbolic, None);
    let epsilon = epsilon(&settings);
    // Two coefficients at w = 0 and two at w = -1, all at the head p(3).
    let terms = [
        ("muvsq+ε", "bubble_12_n1"), // dot(k(1),p(3)) P1L(1,2), w = 0
        ("muvsq", "bubble_21_n1"),   // dot(k(1),p(3)) P1L(2,1), w = 0
        ("ε", "bubble_22_n1"),       // dot(k(1),p(3)) P1L(2,2), w = -1
        ("muvsq^2", "bubble_21_n0"), // P1L(2,1), w = -1
    ];
    let external = parse("p(3)");
    let numerators = terms.map(|(coefficient, name)| {
        parse(&format!(
            "({coefficient})*({})",
            case(&fixture, name).numerator
        ))
    });
    let requests = terms
        .iter()
        .zip(&numerators)
        .map(|((_, name), numerator)| {
            let case = case(&fixture, name);
            PropagatorRequest {
                family: registry_family(&case.family).unwrap(),
                powers: &case.powers,
                canonical_numerator: numerator.as_view(),
                external: external.clone(),
            }
        })
        .collect::<Vec<_>>();
    let values = reduce_requests(&engine, &settings, &requests, None).unwrap();
    assert_eq!(values.keys().collect::<Vec<_>>(), [&(1, external.clone())]);
    assert_eq!(
        values[&(1, external.clone())].len(),
        4,
        "one value per channel"
    );
    let mut expected = BTreeMap::<i64, Atom>::new();
    let mut finalized = Atom::Zero;
    for (coefficient, name) in terms {
        let mut frozen = case(&fixture, name).frozen(epsilon).unwrap();
        for value in &mut frozen {
            value.reduced = &value.reduced * parse(coefficient);
            *expected.entry(value.scale_power).or_default() += &value.reduced;
        }
        finalized += simplify_real(
            finalize::finalize(
                &frozen,
                1,
                &external,
                &settings,
                RustRedMasterMode::Symbolic,
                "Forcer",
            )
            .unwrap()
            .as_view(),
        );
    }
    expected.retain(|_, value| !value.expand().is_zero());
    let rustred = by_scale_power(values.into_values().flatten());
    assert_tag_basis_equal("coefficients", &rustred, &expected, epsilon);

    let input = terms
        .iter()
        .map(|(coefficient, name)| {
            let case = case(&fixture, name);
            parse(&format!(
                "({coefficient})*({})*topo(P1L(p(3),{}))",
                case.numerator,
                case.powers
                    .iter()
                    .map(i64::to_string)
                    .collect::<Vec<_>>()
                    .join(",")
            ))
        })
        .fold(Atom::Zero, |sum, term| sum + term);
    let actual = engine.evaluate(&settings, input.as_view()).unwrap();
    assert!(
        (&actual - &finalized).expand().is_zero(),
        "RustRed {actual} vs Forcer {finalized}"
    );
}

/// A key the library cannot reduce names the vakint integrals of its
/// channel (at most eight, then a count) next to RustRed's owner class and
/// key; other RustRed errors map as for the vacuum programs.
#[test]
fn uncovered_keys_name_the_vakint_integrals() {
    use rustred::propagator::OwnerClassId;

    let integrals = (1..=10)
        .map(|power| format!("P1L(1,{power})"))
        .collect::<BTreeSet<_>>();
    let uncovered = PError::Uncovered {
        class: OwnerClassId(0x2a),
        key: IntegralKey::try_new(vec![1, 2]).unwrap(),
    };
    let RustRedEvaluationError::Reduction { detail } = reduction_error(uncovered, &integrals)
    else {
        panic!("an uncovered key is a reduction error");
    };
    assert!(detail.contains("does not cover [1, 2]"), "{detail}");
    assert!(detail.contains("P1L(1,1), P1L(1,10), P1L(1,2)"), "{detail}");
    assert!(
        detail.ends_with("2 more with the RustRed p library"),
        "{detail}"
    );
    let missing = PError::MissingOwner {
        class: OwnerClassId(7),
        ledger: 3,
    };
    let one = BTreeSet::from(["P2L(1,1,1,1,1)".to_owned()]);
    let error = reduction_error(missing, &one).to_string();
    assert!(error.contains("ledger entry 3"), "{error}");
    assert!(
        error.contains("P2L(1,1,1,1,1) with the RustRed p library"),
        "{error}"
    );
    let other = reduction_error(PError::MixedDimension, &one).to_string();
    assert!(!other.contains("P2L"), "{other}");
}

/// The direct entry points (`EvaluationMethod::evaluate_integral` and
/// `VakintTerm::evaluate_integral`) return the batch value of a p-integral.
#[test]
fn the_direct_entry_points_finalize_once() {
    let engine = engine();
    let fixture = battery();
    for mode in [RustRedMasterMode::Retain, RustRedMasterMode::Symbolic] {
        let settings = settings(mode, None);
        let case = case(&fixture, "bubble_21_n3");
        let input = case.expression().unwrap();
        let batched = engine
            .evaluate_integral(&settings, input.as_view())
            .unwrap();
        let mut term = crate::VakintTerm {
            integral: parse(&format!("topo(P1L({},2,1))", case.external)),
            numerator: parse(&case.numerator),
            vectors: Vec::new(),
        };
        let mut prepared = term.clone();
        assert!(prepared.normalize_propagator().unwrap());
        let matched = prepared
            .prepare_integral_specs(&engine, &settings, None)
            .unwrap();
        let direct = settings.evaluation_order.0[0]
            .evaluate_integral(&engine, &settings, prepared.numerator.as_view(), &matched)
            .unwrap();
        assert!((&direct - &batched).expand().is_zero(), "{mode}");
        term.evaluate_integral(&engine, &settings).unwrap();
        assert!((Atom::from(term) - &batched).expand().is_zero(), "{mode}");
    }
}

/// One step cache per library program, named by the textual rendering of the
/// library identity (an injected library with its own suffix): a read-write
/// directory is written under that name and a second engine attaches it
/// read-only and replays it, recomputing nothing, with the same result; on an
/// empty directory the same read-only engine misses.
#[test]
fn the_step_cache_directory_is_the_textual_library_identity() {
    let directory = tempfile::tempdir().unwrap();
    let fixture = battery();
    let input = case(&fixture, "bubble_32_n4").expression().unwrap();
    let mut writing = settings(RustRedMasterMode::Retain, None);
    writing.ordered_cache = OrderedCacheSettings::read_write(directory.path());
    let first = engine();
    let cold = first.evaluate(&writing, input.as_view()).unwrap();
    let handle = first.p_library().unwrap();
    let name = handle.cache.textual_tables.clone().unwrap();
    assert!(name.starts_with("p-library-1-1_"), "{name}");
    assert!(name.ends_with("_injected"), "{name}");
    assert!(directory.path().join(&name).is_dir(), "{name}");
    let statistics = first.ordered_cache_statistics();
    let program = statistics
        .programs
        .iter()
        .find(|program| program.program == "p-library-1-1")
        .unwrap();
    // No rule-table digest: zero, still 64 hexadecimal digits.
    assert_eq!(program.rule_tables, "0".repeat(64));
    assert_eq!(program.directory, Some(directory.path().join(&name)));
    assert_eq!(program.data_version, "injected");
    assert!(program.native.inserts > 0);
    assert!(program.native.misses > 0);
    let files = std::fs::read_dir(directory.path().join(&name))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect::<std::collections::BTreeSet<_>>();
    assert!(files.contains("ordered-step-cache.txt"), "{files:?}");

    let read_only = |root: &std::path::Path| {
        let mut reading = writing.clone();
        reading.ordered_cache = OrderedCacheSettings::read_only(root);
        let engine = engine();
        assert_eq!(engine.evaluate(&reading, input.as_view()).unwrap(), cold);
        let statistics = engine.ordered_cache_statistics();
        let program = statistics
            .programs
            .iter()
            .find(|program| program.program == "p-library-1-1")
            .expect("the p program is attached");
        assert_eq!(program.mode, OrderedCacheMode::ReadOnly);
        assert_eq!(program.directory, Some(root.join(&name)));
        (statistics.misses, statistics.inserts)
    };
    // Nothing is recomputed: the steps (and the tail table RustRed keeps
    // next to them) replay from the directory.
    assert_eq!(read_only(directory.path()), (0, 0));
    let empty = tempfile::tempdir().unwrap();
    let (misses, inserts) = read_only(empty.path());
    assert!(
        misses > 0 && inserts == 0,
        "{misses} misses, {inserts} inserts"
    );
    // Without a cache the result is the same.
    let mut off = writing;
    off.ordered_cache = OrderedCacheSettings::off();
    assert_eq!(engine().evaluate(&off, input.as_view()).unwrap(), cold);
}

/// The p reduction runs inside a reduction-math span (the benchmark's
/// "RustRed reduction only" line).
#[test]
fn the_p_reduction_is_a_timed_math_span() {
    use crate::rustred_evaluation::experimental::timing::capture_reduction_timing;

    let engine = engine();
    let fixture = battery();
    let input = case(&fixture, "bubble_22_n2").expression().unwrap();
    for arithmetic in [None, Some("exact")] {
        let settings = settings(RustRedMasterMode::Retain, arithmetic);
        let captured = capture_reduction_timing(|| engine.evaluate(&settings, input.as_view()));
        assert!(!captured.result.unwrap().is_zero());
        assert!(captured.timing.complete, "{:?}", captured.timing);
        assert!(captured.timing.interval_count > 0);
        assert_eq!(captured.timing.reduction_math, Some(captured.timing.gross));
    }
}

/// Residual terminals, untabulated masters and vanishing conditions fail
/// closed; a Forcer master renders as `PMaster(name)` with its tag.
#[test]
fn rendering_fails_closed_outside_the_forcer_basis() {
    let engine = engine();
    let handle = engine.p_library().unwrap();
    let context = handle.library.coefficient_context();
    let epsilon = epsilon(&VakintSettings::default());
    let reduction = |basis: PBasis, condition: Option<Coefficient>| PReduction {
        dimension: 1,
        terms: vec![PTerm {
            basis,
            tag: GammaTag::from_pairs([(1, 2), (-2, -1)]).unwrap(),
            coefficient: context.integer(3),
        }],
        nonzero_conditions: condition.into_iter().collect(),
        statistics: Default::default(),
    };
    let master = render(
        &handle,
        &reduction(PBasis::Master("t105".into()), None),
        epsilon,
    )
    .unwrap();
    assert_eq!(master, parse("3*PGamma(1)^2*PGamma(-2)^-1*PMaster(t105)"));
    // The coefficients of one tag and master are combined into one rational
    // function of ε: (d-2)/(d-4) + 1/(d-4) = (d-1)/(d-4) at d = 4-2ε.
    let limits = ExactAlgebraLimits::default();
    let d = context.parameter("d").unwrap();
    let over = |numerator: &Coefficient| {
        let denominator = context.try_sub(&d, &context.integer(4), limits).unwrap();
        context.try_div(numerator, &denominator, limits).unwrap()
    };
    let unit = |coefficient: Coefficient| PTerm {
        basis: PBasis::Unit,
        tag: GammaTag::one(),
        coefficient,
    };
    let split = PReduction {
        dimension: 1,
        terms: vec![
            unit(over(
                &context.try_sub(&d, &context.integer(2), limits).unwrap(),
            )),
            unit(over(&context.integer(1))),
        ],
        nonzero_conditions: Vec::new(),
        statistics: Default::default(),
    };
    assert_eq!(
        render(&handle, &split, epsilon).unwrap(),
        parse("(3-2*ε)/(-2*ε)").together()
    );
    for basis in [
        PBasis::Master("no7".into()),
        PBasis::Residual(rustred::propagator::PTerminalId(1)),
    ] {
        assert!(matches!(
            render(&handle, &reduction(basis, None), epsilon),
            Err(RustRedEvaluationError::MasterCapability { .. })
        ));
    }
    assert!(matches!(
        render(
            &handle,
            &reduction(PBasis::Unit, Some(context.zero())),
            epsilon
        ),
        Err(RustRedEvaluationError::InvalidMatchedFamily { .. })
    ));
}

/// An external-data root (vakint's vendored vacuum data) whose `p_integral/`
/// area lists a scope 1-4 library, as vakint-external-data lists the default
/// program's, that does not load (one placeholder `library.json`): what the
/// area lists decides the dispatch, the library's files only its evaluation.
fn listed_one_to_four_placeholder() -> tempfile::TempDir {
    use sha2::{Digest, Sha256};
    let root = tempfile::tempdir().unwrap();
    for (path, bytes) in super::super::data::vendored_files() {
        if !path.starts_with("p/") {
            let target = root.path().join(path);
            std::fs::create_dir_all(target.parent().unwrap()).unwrap();
            std::fs::write(target, bytes).unwrap();
        }
    }
    let library = root.path().join("p_integral/p-library-1-4");
    std::fs::create_dir_all(&library).unwrap();
    std::fs::write(library.join("library.json"), b"{}").unwrap();
    let area = serde_json::to_vec_pretty(&serde_json::json!({
        "schema": super::super::p_shipped_cache::AREA_SCHEMA,
        "libraries": [{
            "directory": "p-library-1-4",
            "schema": rustred::propagator::P_LIBRARY_SCHEMA,
            "scope": [1, 4],
            "identity": "p-library:scope=1-4:placeholder",
            "files": {"library.json": {"sha256": format!("{:x}", Sha256::digest(b"{}")), "bytes": 2}},
        }],
        "files": {},
    }))
    .unwrap();
    std::fs::write(root.path().join("p_integral/manifest.json"), &area).unwrap();
    let release_path = root.path().join("manifest.json");
    let mut release: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&release_path).unwrap()).unwrap();
    release["p_integral"] = serde_json::json!({
        "manifest": "p_integral/manifest.json",
        "manifest_sha256": format!("{:x}", Sha256::digest(&area)),
    });
    std::fs::write(&release_path, serde_json::to_vec_pretty(&release).unwrap()).unwrap();
    root
}

/// The default order sends a p-integral to RustRed exactly when the engine's
/// p library covers its loop count (plan M3 D-M3-7), whatever the p
/// backend's structural support: without a four-loop library (the built-in
/// one covers one to three loops) four loops keep their Forcer dispatch
/// (default order with a FORM 5); with an injected library only its scope
/// goes to RustRed; with a scope 1-4 library listed in the external data
/// (vakint-external-data's default program) one to four loops go to RustRed,
/// and a listed library that does not load fails them closed with its typed
/// error, never with a Forcer value.
#[test]
fn p_integrals_go_to_rustred_within_the_library_scope_else_to_forcer() {
    use std::io::Write;
    use std::os::unix::fs::PermissionsExt;

    let directory = tempfile::tempdir().unwrap();
    let form = directory.path().join("form");
    {
        let mut script = std::fs::File::create(&form).unwrap();
        writeln!(
            script,
            "#!/bin/sh\necho \"FORM 5.0.0 (vakint test fake) 64-bits\""
        )
        .unwrap();
    }
    std::fs::set_permissions(&form, std::fs::Permissions::from_mode(0o755)).unwrap();
    let settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![
            EvaluationMethod::RustRed(RustRedEvaluationOptions::default()),
            EvaluationMethod::Forcer(ForcerOptions::default()),
        ]),
        form_exe_path: form.to_str().unwrap().into(),
        ..VakintSettings::default()
    };
    let topology = |engine: &Vakint, input: &str| {
        engine
            .match_topologies_for_settings(&settings, parse(input).as_view())
            .unwrap()
            .unwrap()
            .canonical_topology
    };
    let bubble = "topo(P1L(p(1),1,1))";
    let kite = "topo(P2L(p(1),1,1,1,1,1))";
    let four_loops = "topo(P4L_LA4(p(1),1,1,1,1,1,1,1,1,1,1,1))";
    let bare = Vakint::new().unwrap();
    for input in [bubble, kite, four_loops] {
        let topology = topology(&bare, input);
        assert!(crate::rustred_evaluation::supports(
            &settings,
            &topology,
            &RustRedEvaluationOptions::default()
        ));
        let built_in = input != four_loops;
        assert_eq!(
            crate::rustred_evaluation::has_data(&bare, &settings, &topology),
            built_in
        );
        assert!(match bare.evaluation_method(&settings, &topology) {
            Some(EvaluationMethod::RustRed(_)) => built_in,
            Some(EvaluationMethod::Forcer(_)) => !built_in,
            _ => false,
        });
    }
    let injected = engine();
    assert!(matches!(
        injected.evaluation_method(&settings, &topology(&injected, bubble)),
        Some(EvaluationMethod::RustRed(_))
    ));
    assert!(matches!(
        injected.evaluation_method(&settings, &topology(&injected, kite)),
        Some(EvaluationMethod::Forcer(_))
    ));
    let root = listed_one_to_four_placeholder();
    let mut listed = Vakint::new().unwrap();
    listed.load_external_data(root.path(), 1).unwrap();
    for input in [bubble, kite, four_loops] {
        let topology = topology(&listed, input);
        assert!(crate::rustred_evaluation::has_data(
            &listed, &settings, &topology
        ));
        assert!(matches!(
            listed.evaluation_method(&settings, &topology),
            Some(EvaluationMethod::RustRed(_))
        ));
    }
    match listed.evaluate(&settings, parse(four_loops).as_view()) {
        Err(VakintError::RustRedEvaluation(error)) => {
            assert!(error.to_string().contains("p-library-1-4"), "{error}")
        }
        other => panic!("expected the listed library's error, got {other:?}"),
    }
}

/// RustRed is never the fallback method of a p-integral without a p library
/// for its loop count (it stays the fallback that reports missing vacuum
/// data): with RustRed alone and no FORM 5, the expression path reports
/// `ForcerUnavailable` and the term path `NoEvaluationMethodFound` for a
/// four-loop p-integral of an engine without a four-loop library (the
/// built-in library covers one to three loops; the scope 1-4 library comes
/// from vakint-external-data), as before the p backend existed.
#[test]
fn without_p_data_rustred_is_never_the_fallback_of_a_p_integral() {
    let settings = settings(RustRedMasterMode::Symbolic, None);
    let topology = |engine: &Vakint, input: &Atom| {
        engine
            .match_topologies_for_settings(&settings, input.as_view())
            .unwrap()
            .unwrap()
            .canonical_topology
    };
    let bare = Vakint::new().unwrap();
    let input = parse("topo(P4L_LA4(p(1),1,1,1,1,1,1,1,1,1,1,1))");
    assert!(
        bare.evaluation_method(&settings, &topology(&bare, &input))
            .is_none()
    );
    assert!(matches!(
        bare.evaluate(&settings, input.as_view()),
        Err(VakintError::ForcerUnavailable { .. })
    ));
    let mut term = crate::VakintTerm {
        integral: input,
        numerator: Atom::num(1),
        vectors: Vec::new(),
    };
    assert!(matches!(
        term.evaluate_integral(&bare, &settings),
        Err(VakintError::NoEvaluationMethodFound(..))
    ));
    let injected = engine();
    let bubble = parse("topo(P1L(p(1),1,1))");
    assert!(matches!(
        injected.evaluation_method(&settings, &topology(&injected, &bubble)),
        Some(EvaluationMethod::RustRed(_))
    ));
}

/// Every family of RustRed's families file has vakint's canonical lines as its
/// line slots (all 16, one to four loops), and random numerator polynomials
/// lower exactly onto lines and ISPs: at a random point of the scalar
/// products, `Σ c Π (k·k, k·q)^e = Σ c' Π D_s^(-key_s)` (`q·q = -1`).
#[test]
fn numerators_lower_exactly_onto_the_slots_of_every_family() {
    use crate::p_integral::P_FAMILIES;
    use rustred::propagator::families::{P_FAMILIES_V1_JSON, load_p_families};
    use rustred::propagator::presentation::Slot;

    let (file, presentations) = load_p_families(P_FAMILIES_V1_JSON).unwrap();
    assert_eq!(file.families.len(), P_FAMILIES.len());
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = |bound: u64| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 33) % bound
    };
    for family in &P_FAMILIES {
        let index = file
            .families
            .iter()
            .position(|entry| entry.name == family.name)
            .unwrap();
        let presentation = &presentations[index];
        let top =
            super::super::p_data::check_top_presentation(presentation, presentation.id, family)
                .unwrap_or_else(|error| panic!("{}: {error}", family.name));
        let loops = family.loops;
        // k_i·k_j for i <= j at `upper[i * loops + j]`, and k_i·q.
        let upper = (0..loops * loops)
            .map(|_| Rational::from((next(41) as i64 - 20, next(7) as i64 + 1)))
            .collect::<Vec<_>>();
        let loop_loop = |i: usize, j: usize| &upper[i.min(j) * loops + i.max(j)];
        let loop_external = (0..loops)
            .map(|_| Rational::from((next(41) as i64 - 20, next(5) as i64 + 1)))
            .collect::<Vec<_>>();
        let product = |a: &[i64], b: &[i64]| {
            let mut value = Rational::from(-a[loops] * b[loops]);
            for (i, (a_i, external)) in a[..loops].iter().zip(&loop_external).enumerate() {
                value += &(&Rational::from(a_i * b[loops] + a[loops] * b[i]) * external);
                for (j, b_j) in b[..loops].iter().enumerate() {
                    value += &(&Rational::from(a_i * b_j) * loop_loop(i, j));
                }
            }
            value
        };
        let slots = presentation
            .slots
            .iter()
            .map(|slot| match slot {
                Slot::Line(line) => product(&line.momentum, &line.momentum),
                Slot::Isp(isp) => product(&isp.left, &isp.right),
            })
            .collect::<Vec<_>>();
        let dot_value = |dot: PDot| {
            let (a, b) = dot_rows(dot, loops);
            product(&a, &b)
        };
        let dots = (1..=loops)
            .flat_map(|i| (i..=loops).map(move |j| PDot::Loops(i, j)))
            .chain((1..=loops).map(PDot::External))
            .collect::<Vec<_>>();
        let polynomial = (0..4)
            .map(|_| {
                let factors = (0..next(4))
                    .map(|_| (dots[next(dots.len() as u64) as usize], next(2) as i64 + 1))
                    .collect::<Vec<_>>();
                (
                    Fraction::new(next(19) as i128 - 9, next(3) as i128 + 1),
                    PMonomial(factors),
                )
            })
            .collect::<Vec<_>>();
        let mut direct = Rational::zero();
        for (value, monomial) in &polynomial {
            let mut term = fraction(*value);
            for &(dot, exponent) in &monomial.0 {
                for _ in 0..exponent {
                    term = &term * &dot_value(dot);
                }
            }
            direct += &term;
        }
        let lowered = lower_part(
            presentation,
            top,
            &vec![0; top.line_count],
            &polynomial,
            &mut BTreeMap::new(),
        )
        .unwrap();
        let mut via_slots = Rational::zero();
        for (key, value) in &lowered {
            let mut term = value.clone();
            for (slot, &power) in key.iter().enumerate() {
                assert!(power <= 0, "{}: {key:?}", family.name);
                for _ in 0..-power {
                    term = &term * &slots[slot];
                }
            }
            via_slots += &term;
        }
        assert_eq!(direct, via_slots, "{}: {polynomial:?}", family.name);
    }
}

/// What the p backend sends under each Ordered arithmetic (grant U6, review
/// R2-S1). `request_options` carries the settings' arithmetic on `n_cores`
/// workers and the settings' scheduler. A real request, with the inputs
/// `add_request` lowers from a case and built by `p_reduction_request` as
/// `reduce_requests` builds it, reduces under each arithmetic to Forcer's
/// frozen value, and so does `reduce_requests`. A modular option below
/// RustRed's minimum (`max_primes = 1`) in the settings fails
/// `reduce_requests`, which therefore sends the settings' arithmetic.
#[test]
fn every_arithmetic_is_sent_with_real_inputs() {
    use crate::{OrderedWeightArithmetic, RustRedModularCommon, RustRedOrderedWeightArithmetic};

    let engine = engine();
    let handle = engine.p_library().unwrap();
    let fixture = battery();
    let case = case(&fixture, "bubble_32_n4");
    let numerator = parse(&case.numerator);
    let request = PropagatorRequest {
        family: registry_family(&case.family).unwrap(),
        powers: &case.powers,
        canonical_numerator: numerator.as_view(),
        external: parse(&case.external),
    };
    let (mut forms, mut channels) = (BTreeMap::new(), Channels::new());
    add_request(&handle, &request, &mut forms, &mut channels).unwrap();
    assert!(!channels.is_empty());
    for arithmetic in ["exact", "modular_polynomial", "modular_sampled"] {
        let mut settings = settings(RustRedMasterMode::Retain, Some(arithmetic));
        for n_cores in [1, 3] {
            settings.n_cores = n_cores;
            let options = request_options(&settings);
            assert_eq!(options.arithmetic.name(), arithmetic);
            assert_eq!(
                options.arithmetic,
                settings.ordered_weight_arithmetic.resolve(n_cores)
            );
            let workers = if arithmetic == "exact" { 1 } else { n_cores };
            assert_eq!(options.arithmetic.threads(), workers, "{arithmetic}");
            assert_eq!(
                options.scheduler,
                crate::rustred_evaluation::ordered_scheduler(&settings)
            );
        }
        settings.n_cores = 1;
        let options = request_options(&settings);
        let epsilon = epsilon(&settings);
        let mut sent = BTreeMap::<i64, Atom>::new();
        for (channel, lowered) in &channels {
            let inputs = lowered
                .inputs
                .iter()
                .filter(|(_, value)| !value.is_zero())
                .map(|((presentation, powers), value)| {
                    let context = handle.library.coefficient_context();
                    let value = coefficient(context, value, options.limits.exact_algebra).unwrap();
                    (
                        *presentation,
                        IntegralKey::try_new(powers.clone()).unwrap(),
                        value,
                    )
                })
                .collect::<Vec<_>>();
            assert!(!inputs.is_empty());
            let request = p_reduction_request(&options, channel.loops, inputs, None);
            assert_eq!(request.arithmetic.name(), arithmetic);
            let reduction = handle.library.reduce(request).unwrap();
            *sent.entry(channel.scale_power).or_default() +=
                render(&handle, &reduction, epsilon).unwrap() * &channel.coefficient;
        }
        sent.retain(|_, reduced| !reduced.expand().is_zero());
        let forcer = frozen_values(&settings, case);
        assert_tag_basis_equal(arithmetic, &sent, &forcer, epsilon);
        let reduced = rustred_values(&engine, &settings, case);
        assert_tag_basis_equal(arithmetic, &reduced, &forcer, epsilon);
    }
    let below_minimum = RustRedModularCommon {
        max_primes: 1,
        ..RustRedModularCommon::default()
    };
    for invalid in [
        RustRedOrderedWeightArithmetic::ModularSampled(below_minimum, Default::default()),
        RustRedOrderedWeightArithmetic::ModularPolynomial(below_minimum, Default::default()),
    ] {
        let settings = VakintSettings {
            ordered_weight_arithmetic: OrderedWeightArithmetic::new(invalid),
            ..settings(RustRedMasterMode::Retain, None)
        };
        assert!(
            reduce_requests(&engine, &settings, std::slice::from_ref(&request), None).is_err(),
            "{}: max_primes = 1 must fail the sent request",
            invalid.name()
        );
    }
}
