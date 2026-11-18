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
        assert_eq!(
            request.arithmetic,
            default.rustred_ordered_weight_arithmetic()
        );
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
    let values = reduce_requests_unfolded(&engine, &settings, &requests, None, false).unwrap();
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
    // t162-export: a key of a boxed node outside its box names them too.
    let outside = PError::OutsideBox {
        class: OwnerClassId(9),
        key: IntegralKey::try_new(vec![3, 1, 1, 1, 1]).unwrap(),
        node: "KITE".to_owned(),
        definition: "box d1".to_owned(),
        violated: None,
        request: None,
    };
    let error = reduction_error(outside, &one).to_string();
    assert!(error.contains("(boxed node KITE)"), "{error}");
    assert!(
        error.contains("P2L(1,1,1,1,1) with the RustRed p library"),
        "{error}"
    );
}

/// vk5l WP4: RustRed's typed refusals of boxed and open nodes reach the
/// caller unchanged (`RustRedEvaluationError::PRefusal`, `p_refusal`): the
/// node, the key, the violated bound, the original RustRed request, its
/// reproduction command and the vakint integrals; a 20-slot (five-loop) key
/// also names the RustRed build the reproduction needs, a 5-slot key not.
#[test]
fn typed_refusals_pass_through_with_node_key_bound_request_and_reproduction() {
    use rustred::propagator::OwnerClassId;
    use rustred::propagator::error::{BoxViolation, PRequestContext, PRequestInput};
    let five = [1_i64; 14].into_iter().chain([0; 6]).collect::<Vec<_>>();
    let mut dotted = five.clone();
    dotted[6] = 6;
    let request = std::sync::Arc::new(PRequestContext {
        library: "p-library:scope=1-5:test".to_owned(),
        inputs: vec![PRequestInput {
            presentation: PresentationId(41),
            family: "P3e892c5a34c7100d_e15".to_owned(),
            powers: five.clone(),
        }],
        reproduce: "cargo run --release -p rustred-app --example p5l_library -- reduce '/pkg' \
                    p41:1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0"
            .to_owned(),
    });
    let integrals =
        BTreeSet::from(["P3e892c5a34c7100d_e15(1,1,1,1,1,1,1,1,1,1,1,1,1,1)".to_owned()]);
    let cases = [
        (
            PError::OutsideBox {
                class: OwnerClassId(9),
                key: IntegralKey::try_new(dotted.clone()).unwrap(),
                node: "T178".to_owned(),
                definition: "box B".to_owned(),
                violated: Some(Box::new(BoxViolation::LineCap {
                    slot: 6,
                    power: 6,
                    cap: 5,
                })),
                request: Some(request.clone()),
            },
            vec![
                "boxed node T178",
                "line slot 6 power 6 > cap 5",
                "outside its box B",
            ],
        ),
        (
            PError::OutsideCertifiedSet {
                class: OwnerClassId(9),
                key: IntegralKey::try_new(dotted.clone()).unwrap(),
                node: "T178".to_owned(),
                definition: "box B".to_owned(),
                violated: None,
                certified_keys: 335_306,
                request: Some(request.clone()),
            },
            vec![
                "boxed node T178",
                "not in the certified closure",
                "certified set (335306 keys)",
            ],
        ),
        (
            PError::NodeNotClosed {
                class: OwnerClassId(7),
                node: "T266".to_owned(),
                key: IntegralKey::try_new(five.clone()).unwrap(),
                status: "stuck: not closed".to_owned(),
                request: Some(request.clone()),
            },
            vec!["node T266", "status: stuck: not closed"],
        ),
    ];
    for (native, fragments) in cases {
        let error = reduction_error(native.clone(), &integrals);
        assert_eq!(error.p_refusal(), Some(&native));
        let text = crate::VakintError::from(error).to_string();
        for fragment in fragments.into_iter().chain([
            "request: P3e892c5a34c7100d_e15[p41] I(1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0)",
            "reproduce: cargo run --release -p rustred-app --example p5l_library -- reduce '/pkg'",
            "reducing the p-integrals P3e892c5a34c7100d_e15(1,1,1,1,1,1,1,1,1,1,1,1,1,1) with the \
             RustRed p library",
            "RUSTRED_RUNTIME_ARITIES=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,20",
        ]) {
            assert!(text.contains(fragment), "{fragment:?} not in {text}");
        }
    }
    let small = reduction_error(
        PError::NodeNotClosed {
            class: OwnerClassId(7),
            node: "T266".to_owned(),
            key: IntegralKey::try_new(vec![1, 1, 1, 1, 1]).unwrap(),
            status: "stuck".to_owned(),
            request: None,
        },
        &integrals,
    )
    .to_string();
    assert!(!small.contains("RUSTRED_RUNTIME_ARITIES"), "{small}");
}

/// t162-export V6: a library with a BOXED node (RustRed's toy
/// `testing::boxed_kite_library`: scope 1-2, P2L owned by the boxed kite
/// table, whose corner is a boxed terminal, an unknown integral). The corner
/// comes back as `PResidual(KITE,1,1,1,1,1)` (coefficient 1 in the G-scheme
/// reduction) under `master_mode=retain`; `symbolic` and `numerical` refuse
/// it (vk5l WP4: through the p master source, `NoCatalogue` naming it); a key
/// outside the box is a typed refusal naming RustRed's box and the vakint
/// integral.
#[test]
fn a_boxed_node_returns_its_terminals_as_presidual_and_refuses_numbers() {
    let mut engine = Vakint::new().unwrap();
    engine
        .inject_p_library(rustred::propagator::testing::boxed_kite_library().unwrap())
        .unwrap();
    let corner = parse("topo(P2L(p(1),1,1,1,1,1))");
    let residual = parse("PResidual(KITE,1,1,1,1,1)").to_string();
    let retain = settings(RustRedMasterMode::Retain, None);
    let kept = engine.evaluate(&retain, corner.as_view()).unwrap();
    assert_eq!(
        crate::p_integral::finalize::unevaluated_terminals(&kept),
        vec![residual.clone()],
        "{kept}"
    );
    // The reduced value itself is the terminal alone (tag basis, w = 0).
    let p2l = registry_family("P2L").unwrap();
    let values = reduce_requests(
        &engine,
        &retain,
        &[PropagatorRequest {
            family: p2l,
            powers: &[1, 1, 1, 1, 1],
            canonical_numerator: Atom::num(1).as_view(),
            external: parse("p(1)"),
        }],
        None,
    )
    .unwrap();
    let reduced = by_scale_power(values.into_values().flatten());
    assert_eq!(reduced.len(), 1, "{reduced:?}");
    // One scale power (w = 2L - Σa = -1 for the corner).
    let (_, value) = reduced.iter().next().unwrap();
    assert_eq!(
        crate::p_integral::finalize::unevaluated_terminals(value),
        vec![residual.clone()],
        "{reduced:?}"
    );
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let error = engine
            .evaluate(&settings(mode, None), corner.as_view())
            .unwrap_err()
            .to_string();
        assert!(error.contains("no p master catalogue"), "{mode}: {error}");
        assert!(
            error.contains("PResidual(KITE,1,1,1,1,1)"),
            "{mode}: {error}"
        );
        assert!(error.contains("master_mode=retain"), "{mode}: {error}");
    }
    let outside = engine
        .evaluate(&retain, parse("topo(P2L(p(1),3,1,1,1,1))").as_view())
        .unwrap_err();
    let crate::VakintError::RustRedEvaluation(native) = &outside else {
        panic!("not a RustRed error: {outside}");
    };
    assert!(
        matches!(native.p_refusal(), Some(PError::OutsideBox { node, .. }) if node == "KITE"),
        "{outside}"
    );
    let outside = outside.to_string();
    assert!(outside.contains("(boxed node KITE)"), "{outside}");
    assert!(outside.contains("outside its box"), "{outside}");
    assert!(outside.contains("P2L(3,1,1,1,1)"), "{outside}");
    assert!(outside.contains("reproduce:"), "{outside}");
}

/// vk5l WP4 (decision Q9): under `rustred_p_uncovered = Retain` and
/// `master_mode=retain` a key the boxed kite refuses stays in the result as
/// `PResidual(KITE, ..)` (the vacuum path's retained unresolved masters);
/// `Reject` (the default) and the other master modes keep the typed refusal.
/// The certified kite (box `d0`, certified set = corner + five dotted keys)
/// refuses a doubly dotted key too.
#[test]
fn the_retain_policy_keeps_refused_keys_as_presidual_in_retain_mode_only() {
    for library in [
        rustred::propagator::testing::boxed_kite_library().unwrap(),
        rustred::propagator::testing::boxed_kite_library_certified().unwrap(),
    ] {
        let mut engine = Vakint::new().unwrap();
        engine.inject_p_library(library).unwrap();
        let outside = parse("topo(P2L(p(1),3,1,1,1,1))");
        let mut retain = settings(RustRedMasterMode::Retain, None);
        assert!(engine.evaluate(&retain, outside.as_view()).is_err());
        retain.rustred_p_uncovered = crate::RustRedPUncoveredPolicy::Retain;
        let kept = engine.evaluate(&retain, outside.as_view()).unwrap();
        let residuals = crate::p_integral::finalize::unevaluated_terminals(&kept);
        assert_eq!(residuals.len(), 1, "{kept}");
        assert!(residuals[0].contains("KITE"), "{kept}");
        let mut symbolic = settings(RustRedMasterMode::Symbolic, None);
        symbolic.rustred_p_uncovered = crate::RustRedPUncoveredPolicy::Retain;
        let refused = engine
            .evaluate(&symbolic, outside.as_view())
            .unwrap_err()
            .to_string();
        assert!(refused.contains("boxed node KITE"), "{refused}");
    }
}

/// vk5l WP4: the p path writes its requests as the vacuum path does
/// (`VAKINT_DUMP_ORDERED_REQUEST`): the request (library identity, inputs,
/// vakint integrals, RustRed's reproduction command) and its result digest,
/// equal for two identical requests; a refused request records its error.
#[test]
fn a_p_request_and_its_result_are_dumped() {
    use super::super::request_dump::{PDumpRequest, PRequestDump};
    let directory = tempfile::tempdir().unwrap();
    let library = rustred::propagator::testing::boxed_kite_library().unwrap();
    let options = request_options(&settings(RustRedMasterMode::Retain, None));
    let one = library.coefficient_context().one();
    let integrals = BTreeSet::from(["P2L(1,1,1,1,1)".to_owned()]);
    let p2l = library.top_presentation("P2L").unwrap();
    let mut digests = Vec::new();
    for powers in [[1, 1, 1, 1, 1], [1, 1, 1, 1, 1], [3, 1, 1, 1, 1]] {
        let inputs = vec![(p2l, IntegralKey::try_new(powers).unwrap(), one.clone())];
        let dump = PRequestDump::write_to(
            directory.path(),
            PDumpRequest {
                library: &library,
                loops: 2,
                label: "test",
                arithmetic: options.arithmetic_for(2),
                uncovered: options.uncovered,
                master_mode: options.master_mode,
                inputs: &inputs,
                integrals: &integrals,
            },
        )
        .unwrap();
        let stem = dump.stem().to_owned();
        let result = library.reduce(p_reduction_request(&options, 2, inputs, None));
        dump.finish(&library, &result).unwrap();
        let read = |suffix: &str| -> serde_json::Value {
            serde_json::from_slice(
                &std::fs::read(directory.path().join(format!("{stem}.{suffix}.json"))).unwrap(),
            )
            .unwrap()
        };
        let request = read("request");
        assert_eq!(request["schema"], "vakint.p-request.v1");
        assert_eq!(request["inputs"][0]["family"], "P2L");
        assert_eq!(request["inputs"][0]["powers"], serde_json::json!(powers));
        assert_eq!(request["vakint_integrals"][0], "P2L(1,1,1,1,1)");
        assert!(
            request["reproduce"]
                .as_str()
                .unwrap()
                .contains("p5l_library -- reduce"),
            "{request}"
        );
        let result = read("result");
        digests.push((
            result["output"]["terms_blake3"].clone(),
            result["error"].clone(),
        ));
    }
    assert_eq!(digests[0], digests[1]);
    assert!(
        digests[0].0.is_string() && digests[0].1.is_null(),
        "{digests:?}"
    );
    assert!(digests[2].0.is_null(), "{digests:?}");
    assert!(
        digests[2].1.as_str().unwrap().contains("boxed node KITE"),
        "{digests:?}"
    );
    let index = std::fs::read_to_string(directory.path().join("index.jsonl")).unwrap();
    assert_eq!(index.lines().count(), 6);
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
        outside_certified: vec![],
        relation_cache: Default::default(),
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
        outside_certified: vec![],
        relation_cache: Default::default(),
    };
    assert_eq!(
        render(&handle, &split, epsilon).unwrap(),
        parse("(3-2*ε)/(-2*ε)").together()
    );
    assert!(matches!(
        render(
            &handle,
            &reduction(PBasis::Master("no7".into()), None),
            epsilon
        ),
        Err(RustRedEvaluationError::MasterCapability { .. })
    ));
    // t162-export: a residual renders as PResidual(..) when the library has
    // it; the reference library has none (terminal 1 does not exist).
    assert!(matches!(
        render(
            &handle,
            &reduction(PBasis::Residual(rustred::propagator::PTerminalId(1)), None),
            epsilon
        ),
        Err(RustRedEvaluationError::ExternalData { .. })
    ));
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

/// Every family of RustRed's families files has vakint's canonical lines as its
/// line slots (all 80: 16 of one to four loops, 64 of five loops, t162-export),
/// and random numerator polynomials
/// lower exactly onto lines and ISPs: at a random point of the scalar
/// products, `Σ c Π (k·k, k·q)^e = Σ c' Π D_s^(-key_s)` (`q·q = -1`).
#[test]
fn numerators_lower_exactly_onto_the_slots_of_every_family() {
    use crate::p_integral::P_FAMILIES;
    use rustred::propagator::families::{
        P_FAMILIES_5L_V1_JSON, P_FAMILIES_V1_JSON, load_p_families,
    };
    use rustred::propagator::presentation::Slot;

    let (file, presentations) = load_p_families(P_FAMILIES_V1_JSON).unwrap();
    let (five, five_presentations) = load_p_families(P_FAMILIES_5L_V1_JSON).unwrap();
    let names = (file.families.iter())
        .chain(&five.families)
        .map(|entry| entry.name.clone())
        .collect::<Vec<_>>();
    let presentations = (presentations.into_iter())
        .chain(five_presentations)
        .collect::<Vec<_>>();
    assert_eq!(names.len(), P_FAMILIES.len());
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut next = |bound: u64| {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (state >> 33) % bound
    };
    for family in &P_FAMILIES {
        let index = names.iter().position(|name| name == family.name).unwrap();
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

/// t162-export V6: five-loop smoke test on a scope 1-5 package that
/// `VAKINT_TEST_P_LIBRARY_1_5` names (skipped when unset): the package loads
/// as vakint loads a library (no descent certificate), passes the handle's
/// checks (tops of 1-4 complete, five-loop tops optional, residual terminals
/// accepted), and the short form `topo(TOP(p(1), a_1..a_14))` of the top
/// `VAKINT_TEST_P_5L_TOP` (default `P3dc41bbf2115e501_e15`, the s1 package's
/// top) with the powers `VAKINT_TEST_P_5L_POWERS` (comma separated, default
/// all 1) reduces with `master_mode=retain`; the result is printed with its
/// `PResidual` terminals and the wall time.
#[test]
fn five_loop_smoke_on_a_one_to_five_package() {
    let Some(directory) = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_5") else {
        eprintln!("five_loop_smoke: VAKINT_TEST_P_LIBRARY_1_5 unset, skipped");
        return;
    };
    let started = std::time::Instant::now();
    let library = rustred_app::propagator::load_p_library(
        rustred_app::propagator::PLibrarySource::Directory(directory.clone().into()),
        &rustred_app::propagator::PLibraryOptions {
            verify_descent: false,
        },
    )
    .unwrap_or_else(|error| panic!("{}: {error}", directory.display()));
    assert_eq!(*library.scope().end(), 5);
    let mut engine = Vakint::new().unwrap();
    engine.inject_p_library(library).unwrap();
    let loaded = started.elapsed().as_secs_f64();
    let top = std::env::var("VAKINT_TEST_P_5L_TOP")
        .unwrap_or_else(|_| "P3dc41bbf2115e501_e15".to_owned());
    let family = registry_family(&top).unwrap();
    assert_eq!(family.loops, 5);
    let handle = engine.p_library().unwrap();
    let index = crate::p_integral::P_FAMILIES
        .iter()
        .position(|entry| entry.name == top)
        .unwrap();
    handle
        .top(index)
        .unwrap_or_else(|error| panic!("{top}: {error}"));
    let powers = std::env::var("VAKINT_TEST_P_5L_POWERS").map_or_else(
        |_| vec!["1".to_owned(); family.line_count()],
        |powers| powers.split(',').map(str::to_owned).collect(),
    );
    assert_eq!(powers.len(), family.line_count());
    let input = parse(&format!("topo({top}(p(1),{}))", powers.join(",")));
    let started = std::time::Instant::now();
    let value = engine
        .evaluate(&settings(RustRedMasterMode::Retain, None), input.as_view())
        .unwrap_or_else(|error| panic!("{input}: {error}"));
    let residuals = crate::p_integral::finalize::unevaluated_terminals(&value);
    eprintln!(
        "five_loop_smoke: {input}: load {loaded:.1} s, evaluate {:.1} s, {} PResidual terminals {residuals:?}",
        started.elapsed().as_secs_f64(),
        residuals.len()
    );
    assert!(!value.is_zero(), "{input}");
}

/// vk5l WP4: the guards of the five-loop vacuum path on the p path at five
/// loops, through the engine (`Vakint::set_p_master_catalogue`, the hook
/// `p_master_catalogue` and `finalize_reduced`, the step after the RustRed
/// reduction of `evaluate_ordered_batch`) on WP5's synthetic five-loop
/// catalogue (`master_catalogue::tests::fixture`), with the vacuum tests'
/// shared helper `truncation_guard::expect_missing_orders`. Every refusal
/// names the five-loop vakint integral of the request. The guard "input
/// outside the box" is `a_boxed_node_returns_its_terminals_as_presidual_and_refuses_numbers`
/// and `typed_refusals_pass_through_with_node_key_bound_request_and_reproduction`
/// (and, on a real package, `five_loop_package_drill`).
mod five_loop_guards {
    use super::*;
    use crate::p_integral::master_catalogue::tests::fixture;
    use crate::p_integral::{PMasterCatalogue, PMasterCatalogueError, p_residual_symbol};
    use crate::truncation_guard::expect_missing_orders;

    const INTEGRAL: &str = "P3e892c5a34c7100d_e15(1,1,1,1,1,1,1,1,1,1,1,1,1,1)";

    fn engine(catalogue: bool) -> Vakint {
        let mut engine = Vakint::new().unwrap();
        engine.disable_default_p_master_catalogue();
        if catalogue {
            engine.set_p_master_catalogue(
                PMasterCatalogue::from_json("test catalogue", &fixture().to_string()).unwrap(),
            );
        }
        engine
    }

    fn residual(node: &str, key: &[i64]) -> Atom {
        let mut arguments = vec![Atom::var(vk_symbol!(node))];
        arguments.extend(key.iter().map(|&k| Atom::num(k)));
        p_residual_symbol().call(&arguments[..])
    }

    /// The p path after the reduction: `reduced` (scale power 0, external
    /// `p(1)`) of the five-loop request [`INTEGRAL`], finalized in `mode`
    /// through `terms` Laurent terms at `digits` digits.
    fn finalize(
        engine: &Vakint,
        reduced: Atom,
        mode: RustRedMasterMode,
        terms: i64,
        digits: u32,
    ) -> Result<Atom, crate::VakintError> {
        let mut settings = settings(mode, None);
        settings.number_of_terms_in_epsilon_expansion = terms;
        settings.run_time_decimal_precision = digits;
        let head = parse("p(1)");
        finalize_reduced(
            engine,
            &settings,
            BTreeMap::from([(
                (5, head.clone()),
                vec![ReducedPropagatorValue {
                    scale_power: 0,
                    reduced,
                    spectator: None,
                }],
            )]),
            &BTreeMap::from([((5, head), BTreeSet::from([INTEGRAL.to_owned()]))]),
            mode,
        )
    }

    fn catalogue_error(error: crate::VakintError) -> PMasterCatalogueError {
        match error {
            crate::VakintError::PMasterCatalogue(error) => error,
            other => panic!("expected a p master catalogue error, got: {other}"),
        }
    }

    /// Required number of coefficients (vacuum: `master_truncation_guard_tests`):
    /// through `ε^0` the known master suffices; one Laurent term more names
    /// the master, its first missing order and the five-loop request.
    #[test]
    fn required_number_of_coefficients() {
        let engine = engine(true);
        for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
            assert!(
                !finalize(&engine, residual("T900", &[1, 1, 0]), mode, 6, 40)
                    .unwrap()
                    .is_zero()
            );
            let missing = expect_missing_orders(
                finalize(&engine, residual("T900", &[1, 1, 0]), mode, 7, 40).unwrap_err(),
            );
            assert_eq!(missing.requested_order, Some(1), "{missing}");
            assert_eq!(missing.masters.len(), 1, "{missing}");
            assert_eq!(missing.masters[0].first_missing_order, 1, "{missing}");
            assert!(missing.backend.contains(INTEGRAL), "{missing}");
        }
    }

    /// Spurious poles after truncation error (vacuum:
    /// `spurious_poles_require_deeper_master_coefficients`).
    #[test]
    fn spurious_poles_after_truncation_error() {
        let engine = engine(true);
        let missing = expect_missing_orders(
            finalize(
                &engine,
                residual("T900", &[1, 1, 1]),
                RustRedMasterMode::Symbolic,
                6,
                40,
            )
            .unwrap_err(),
        );
        assert_eq!(missing.requested_order, Some(0), "{missing}");
        assert_eq!(missing.masters[0].first_affected_result_order, 0);
        assert!(missing.backend.contains(INTEGRAL), "{missing}");
    }

    /// Missing master (vacuum: `missing_coefficients_are_not_assumed_zero`):
    /// a master without coefficients, a terminal without a map row, and a
    /// five-loop terminal without any catalogue are refused by name; Retain
    /// keeps them.
    #[test]
    fn missing_master() {
        let with = engine(true);
        for (terminal, expected) in [
            (
                residual("T903", &[1]),
                "declared terminal has no computed coefficients",
            ),
            (residual("T901", &[2]), "no master map row yet"),
            (
                residual("T902", &[1]),
                "additional terminal of a boxed node",
            ),
        ] {
            assert!(
                !finalize(&with, terminal.clone(), RustRedMasterMode::Retain, 6, 40)
                    .unwrap()
                    .is_zero()
            );
            let error = catalogue_error(
                finalize(&with, terminal, RustRedMasterMode::Symbolic, 6, 40).unwrap_err(),
            )
            .to_string();
            assert!(
                error.contains(expected) && error.contains(INTEGRAL),
                "{error}"
            );
        }
        let without = engine(false);
        let error = catalogue_error(
            finalize(
                &without,
                residual("T900", &[1, 1, 0]),
                RustRedMasterMode::Numerical,
                6,
                40,
            )
            .unwrap_err(),
        );
        assert!(
            matches!(error, PMasterCatalogueError::NoCatalogue { .. }),
            "{error}"
        );
        assert!(error.to_string().contains(INTEGRAL), "{error}");
    }

    /// Required digits (vacuum: `check_coefficients` and `PrecisionWarning`):
    /// fewer validated digits than requested only warn, no validated digits
    /// are an error naming master and order.
    #[test]
    fn required_digits() {
        let engine = engine(true);
        let terminal = residual("T904", &[1]);
        assert!(
            !finalize(
                &engine,
                terminal.clone(),
                RustRedMasterMode::Numerical,
                6,
                20
            )
            .unwrap()
            .is_zero()
        );
        assert!(
            !finalize(
                &engine,
                terminal.clone(),
                RustRedMasterMode::Numerical,
                6,
                40
            )
            .unwrap()
            .is_zero()
        );
        let mut file = fixture();
        file["masters"][1]["coefficients"]["0"]
            .as_object_mut()
            .unwrap()
            .remove("significant_digits");
        let mut undigited = Vakint::new().unwrap();
        undigited.set_p_master_catalogue(
            PMasterCatalogue::from_json("test", &file.to_string()).unwrap(),
        );
        let error = catalogue_error(
            finalize(&undigited, terminal, RustRedMasterMode::Numerical, 6, 5).unwrap_err(),
        )
        .to_string();
        assert!(
            error.contains("syn-num")
                && error.contains("without validated significant digits")
                && error.contains(INTEGRAL),
            "{error}"
        );
    }

    /// One to four loops keep the Forcer table: the catalogue is a
    /// five-loop source only.
    #[test]
    fn the_catalogue_serves_five_loops_only() {
        let engine = engine(true);
        assert!(p_master_catalogue(&engine, 5).is_some());
        assert!(p_master_catalogue(&engine, 4).is_none());
    }
}

/// Number of summands of an expanded drill result (0 for zero).
fn drill_term_count(value: &Atom) -> usize {
    match value.as_view() {
        symbolica::atom::AtomView::Add(add) => add.get_nargs(),
        _ if value.is_zero() => 0,
        _ => 1,
    }
}

/// vk5l WP4 REAL TEST (env-gated, not in `test_vakint` unless set): a
/// scope 1-5 package assembled from existing tables, named by
/// `VAKINT_TEST_P5L_DRILL_PACKAGE` (loaded through `VAKINT_P_LIBRARY`'s code
/// path, `load_p_library_directory`, unverified), reduced through vakint's
/// public API. `VAKINT_TEST_P5L_DRILL` is a `;`-separated list of
/// `<label>=<expected>=<input>` with `<expected>` one of `ok` (Retain mode,
/// printed with its `PResidual` terminals), `outside-box`,
/// `outside-certified-set`, `outside` (either), `node-not-closed` (the typed
/// refusal, printed);
/// `VAKINT_TEST_P5L_DRILL_RETAIN=1` sets the Retain coverage policy.
/// `<expected>` = `zero`: the input is an exact identity (a sum, e.g. an IBP
/// relation); each summand is reduced on its own (must be nonzero), and the
/// sum of those reductions and the reduction of the whole input must vanish.

/// The largest absolute numerical factor of the terms of `values` (expanded):
/// the scale of a numerical identity check of [`five_loop_package_drill`].
fn drill_largest_coefficient(values: &[Atom]) -> f64 {
    let mut largest = 0.0f64;
    for value in values {
        let expanded = value.expand();
        let terms: Vec<symbolica::atom::AtomView> = match expanded.as_view() {
            symbolica::atom::AtomView::Add(sum) => sum.iter().collect(),
            other => vec![other],
        };
        for term in terms {
            let factors: Vec<symbolica::atom::AtomView> = match term {
                symbolica::atom::AtomView::Mul(product) => product.iter().collect(),
                other => vec![other],
            };
            for factor in factors {
                if let symbolica::atom::AtomView::Num(_) = factor {
                    let text = factor.to_string();
                    let parsed = text
                        .split('/')
                        .map(|part| {
                            part.trim()
                                .trim_end_matches(|c: char| c == '`' || c.is_ascii_digit() && false)
                                .parse::<f64>()
                        })
                        .collect::<Vec<_>>();
                    let magnitude = match parsed.as_slice() {
                        [Ok(a)] => a.abs(),
                        [Ok(a), Ok(b)] => (a / b).abs(),
                        _ => text
                            .split('`')
                            .next()
                            .and_then(|head| head.parse::<f64>().ok())
                            .map_or(0.0, f64::abs),
                    };
                    largest = largest.max(magnitude);
                }
            }
        }
    }
    largest
}

#[test]
fn five_loop_package_drill() {
    let Some(directory) = std::env::var_os("VAKINT_TEST_P5L_DRILL_PACKAGE") else {
        eprintln!("five_loop_package_drill: VAKINT_TEST_P5L_DRILL_PACKAGE unset, skipped");
        return;
    };
    let started = std::time::Instant::now();
    let mut engine = Vakint::new().unwrap();
    engine
        .load_p_library_directory(std::path::Path::new(&directory), false)
        .unwrap_or_else(|error| panic!("{}: {error}", directory.display()));
    eprintln!(
        "five_loop_package_drill: loaded {} in {:.1} s (scope {:?})",
        directory.display(),
        started.elapsed().as_secs_f64(),
        engine.p_library().unwrap().scope()
    );
    // vk5l WP4: VAKINT_TEST_P5L_DRILL_MODE=symbolic|numerical runs the
    // drill end to end through the p master catalogue (default retain).
    let drill_mode = match std::env::var("VAKINT_TEST_P5L_DRILL_MODE").as_deref() {
        Ok("symbolic") => RustRedMasterMode::Symbolic,
        Ok("numerical") => RustRedMasterMode::Numerical,
        _ => RustRedMasterMode::Retain,
    };
    let mut settings = settings(drill_mode, None);
    settings.number_of_terms_in_epsilon_expansion = 6;
    // lane vk5l-mastermap: the five-loop master basis of the drill's results
    // (VAKINT_TEST_P5L_DRILL_BASIS=pole_bounded; default original).
    if std::env::var("VAKINT_TEST_P5L_DRILL_BASIS").as_deref() == Ok("pole_bounded") {
        settings.p_master_basis = crate::p_integral::PMasterBasis::PoleBounded;
        eprintln!(
            "five_loop_package_drill: master basis {:?}",
            engine.p_master_basis_record(&settings).unwrap()
        );
    }
    settings.n_cores = std::env::var("VAKINT_TEST_P5L_DRILL_CORES")
        .ok()
        .and_then(|cores| cores.parse().ok())
        .unwrap_or(1);
    if std::env::var("VAKINT_TEST_P5L_DRILL_RETAIN").as_deref() == Ok("1") {
        settings.rustred_p_uncovered = crate::RustRedPUncoveredPolicy::Retain;
    }
    // Large inputs (physical requests) come from a file.
    let cases = match std::env::var_os("VAKINT_TEST_P5L_DRILL_FILE") {
        Some(file) => std::fs::read_to_string(file).unwrap(),
        None => std::env::var("VAKINT_TEST_P5L_DRILL").unwrap_or_default(),
    };
    if let Ok(epsilon) = std::env::var("VAKINT_TEST_P5L_DRILL_EPSILON") {
        settings.epsilon_symbol = epsilon;
    }
    for case in cases.split(';').filter(|case| !case.trim().is_empty()) {
        let mut parts = case.splitn(3, '=');
        let (label, expected, input) = (
            parts.next().unwrap().trim(),
            parts.next().unwrap().trim(),
            parts.next().unwrap().trim(),
        );
        let input = parse(input);
        let dispatch_before = engine.p_dispatch_counters();
        let started = std::time::Instant::now();
        let result = engine.evaluate(&settings, input.as_view());
        let seconds = started.elapsed().as_secs_f64();
        eprintln!(
            "DISPATCH {label}: {:?}",
            engine.p_dispatch_counters().since(&dispatch_before)
        );
        match (expected, result) {
            // An exact identity (e.g. an IBP relation with integer
            // coefficients): every summand reduced on its own must be
            // nonzero, and both the sum of those reductions and the
            // reduction of the whole input must vanish identically.
            ("zero", Ok(value)) => {
                let summands: Vec<Atom> = match input.as_view() {
                    symbolica::atom::AtomView::Add(add) => {
                        add.iter().map(|t| t.to_owned()).collect()
                    }
                    _ => vec![input.clone()],
                };
                let mut total = Atom::Zero;
                let mut nonzero = 0usize;
                let mut summand_values = Vec::<Atom>::new();
                for (index, summand) in summands.iter().enumerate() {
                    let started = std::time::Instant::now();
                    let reduced = engine
                        .evaluate(&settings, summand.as_view())
                        .unwrap_or_else(|error| panic!("{label} summand {index}: {error}"));
                    let reduced = reduced.expand();
                    let residuals = crate::p_integral::finalize::unevaluated_terminals(&reduced);
                    eprintln!(
                        "DRILL {label} summand {index}: {:.1} s, {} terms, {} PResidual terminals {residuals:?}\n  \
                         input {summand}\n  value {reduced}",
                        started.elapsed().as_secs_f64(),
                        drill_term_count(&reduced),
                        residuals.len()
                    );
                    // A summand may vanish on its own (e.g. a numerator odd
                    // under a symmetry of its sector); the identity needs at
                    // least two nonzero summands to be a check.
                    if reduced.is_zero() {
                        eprintln!("DRILL {label} summand {index}: reduces to 0 on its own");
                    } else {
                        nonzero += 1;
                    }
                    summand_values.push(reduced.clone());
                    total = total + reduced;
                }
                // Rational coefficients in eps (e.g. 1/(1+eps) from a
                // dimension-dependent IBP coefficient) cancel only over a
                // common denominator: `expand` alone leaves 15*eps/(1+eps)
                // next to -15*eps (the exact zero needs `together`).
                let total = total.expand().together().expand();
                let whole = value.expand().together().expand();
                eprintln!(
                    "DRILL {label}: ZERO CHECK whole input in {seconds:.1} s: {} terms ({whole}); \
                     sum of the {} separate reductions: {} terms ({total})",
                    drill_term_count(&whole),
                    summands.len(),
                    drill_term_count(&total)
                );
                assert!(
                    nonzero >= 2,
                    "{label}: fewer than two nonzero summands, no check"
                );
                if drill_mode == RustRedMasterMode::Numerical {
                    // Numerical masters: every coefficient of the sum must
                    // vanish to the working precision, relative to the
                    // largest summand coefficient.
                    let scale = drill_largest_coefficient(&summand_values);
                    let residue = drill_largest_coefficient(&[total.clone(), whole.clone()]);
                    eprintln!(
                        "DRILL {label}: numerical check: largest |coefficient| of the sums {residue:.3e}, \
                         of the summands {scale:.3e}"
                    );
                    assert!(
                        residue <= 1e-25 * scale.max(1.0),
                        "{label}: identity violated numerically"
                    );
                    eprintln!(
                        "DRILL {label}: IDENTITY HOLDS TO WORKING PRECISION ({} summands, {nonzero} nonzero)",
                        summands.len()
                    );
                    continue;
                } else {
                    assert!(
                        whole.is_zero() && total.is_zero(),
                        "{label}: identity violated"
                    );
                }
                eprintln!(
                    "DRILL {label}: IDENTITY HOLDS EXACTLY ({} summands, {nonzero} nonzero)",
                    summands.len()
                );
            }
            ("ok", Ok(value)) => {
                // vk5l WP4 showcase: VAKINT_TEST_P5L_DRILL_DUMP=<dir> writes the
                // value as Symbolica's canonical string (readable by Symbolica's
                // parser and its MathML viewer) and a JSON summary.
                if let Some(directory) = std::env::var_os("VAKINT_TEST_P5L_DRILL_DUMP") {
                    let directory = std::path::PathBuf::from(directory);
                    std::fs::create_dir_all(&directory).unwrap();
                    let mode = format!("{drill_mode:?}").to_lowercase();
                    let text = value.to_canonical_string();
                    std::fs::write(directory.join(format!("{label}.{mode}.sym")), &text).unwrap();
                    std::fs::write(
                        directory.join(format!("{label}.{mode}.json")),
                        serde_json::json!({
                            "label": label,
                            "mode": mode,
                            "seconds": seconds,
                            "terms": drill_term_count(&value.expand()),
                            "bytes": text.len(),
                            "input_bytes": input.to_canonical_string().len(),
                            "presidual_terminals": crate::p_integral::finalize::unevaluated_terminals(&value),
                        })
                        .to_string(),
                    )
                    .unwrap();
                }
                let residuals = crate::p_integral::finalize::unevaluated_terminals(&value);
                eprintln!(
                    "DRILL {label}: OK in {seconds:.1} s, {} terms, {} PResidual terminals {residuals:?}\n  \
                     input {input}\n  value {value}",
                    drill_term_count(&value.expand()),
                    residuals.len()
                );
                assert!(!value.is_zero(), "{label}: {input}");
            }
            (expected, Err(error)) if expected != "ok" => {
                let crate::VakintError::RustRedEvaluation(native) = &error else {
                    panic!("{label}: not a RustRed error: {error}");
                };
                let refusal = native
                    .p_refusal()
                    .unwrap_or_else(|| panic!("{label}: not a typed refusal: {error}"));
                let kind = match refusal {
                    PError::OutsideBox { .. } => "outside-box",
                    PError::OutsideCertifiedSet { .. } => "outside-certified-set",
                    PError::NodeNotClosed { .. } => "node-not-closed",
                    _ => "other",
                };
                eprintln!(
                    "DRILL {label}: {kind} in {seconds:.1} s\n  input {input}\n  ERROR {error}"
                );
                // `outside`: either refusal of a boxed node (the early box
                // check or the certified-set check at dispatch).
                let accepted = kind == expected
                    || (expected == "outside"
                        && matches!(kind, "outside-box" | "outside-certified-set"));
                assert!(
                    accepted,
                    "{label}: expected {expected}, got {kind}: {error}"
                );
            }
            (expected, result) => panic!("{label}: expected {expected}, got {result:?}"),
        }
    }
}

/// vk5l WP4, persistent relation cache, run as TWO processes of this binary
/// against one cache directory (environment-driven, skipped without it):
/// `VAKINT_TEST_P5L_RELATION_PACKAGE` = a package directory whose boxed node
/// carries a recipe bundle; `VAKINT_TEST_P5L_RELATION_INPUT` = an integral that
/// reaches a key solved by a recipe at run time (outside the certified set,
/// or a certified recipe key with no shipped relation);
/// `VAKINT_TEST_P5L_RELATION_DIR` = the cache directory (moved there with
/// `PRelationCacheLocation::Directory`); `VAKINT_TEST_P5L_RELATION_ROLE` =
/// `cold` (first process: at least one recipe solve, inserted into the cache,
/// value written to `<dir>/value.txt`) or `warm` (second process: the same
/// value, cache hits, no solve).
#[test]
fn p_relation_cache_two_processes() {
    let (Some(directory), Some(input), Some(cache), Ok(role)) = (
        std::env::var_os("VAKINT_TEST_P5L_RELATION_PACKAGE"),
        std::env::var("VAKINT_TEST_P5L_RELATION_INPUT").ok(),
        std::env::var_os("VAKINT_TEST_P5L_RELATION_DIR"),
        std::env::var("VAKINT_TEST_P5L_RELATION_ROLE"),
    ) else {
        eprintln!("p_relation_cache_two_processes: environment unset, skipped");
        return;
    };
    let cache = std::path::PathBuf::from(cache);
    let mut engine = Vakint::new().unwrap();
    engine
        .load_p_library_directory(std::path::Path::new(&directory), false)
        .unwrap_or_else(|error| panic!("{}: {error}", directory.to_string_lossy()));
    let mut settings = settings(RustRedMasterMode::Retain, None);
    settings.number_of_terms_in_epsilon_expansion = 6;
    settings.rustred_p_relation_cache.location =
        crate::PRelationCacheLocation::Directory(cache.clone());
    if let Some(shipped) = std::env::var_os("VAKINT_TEST_P5L_RELATION_SHIPPED") {
        settings.rustred_p_relation_cache.shipped = std::env::split_paths(&shipped).collect();
    }
    let dispatch_before = engine.p_dispatch_counters();
    let before = engine.p_relation_cache_statistics();
    let started = std::time::Instant::now();
    let value = engine
        .evaluate(&settings, parse(&input).as_view())
        .unwrap_or_else(|error| panic!("{input}: {error}"))
        .expand()
        .together()
        .expand();
    let after = engine.p_relation_cache_statistics();
    eprintln!(
        "DISPATCH relation-cache {role}: {:?}",
        engine.p_dispatch_counters().since(&dispatch_before)
    );
    let solves = after.recipes.solves - before.recipes.solves;
    let delta = after.cache.since(&before.cache);
    eprintln!(
        "RELATION-CACHE {role}: {:.1} s, {} terms, solves={solves} (outside {}), {}, directory {:?}, \
         bytes_on_disk={}",
        started.elapsed().as_secs_f64(),
        drill_term_count(&value),
        after.recipes.outside_solves - before.recipes.outside_solves,
        crate::rustred_evaluation::p_relations::describe(&delta),
        after.directory,
        after.bytes_on_disk
    );
    assert_eq!(after.directory.as_deref(), Some(cache.as_path()));
    let saved = cache.join("value.txt");
    match role.as_str() {
        "cold" => {
            assert!(
                solves > 0,
                "cold: no recipe solve; the input reaches no run-time key"
            );
            assert!(delta.inserts > 0, "cold: nothing inserted");
            assert!(after.bytes_on_disk > 0);
            std::fs::write(&saved, value.to_canonical_string()).unwrap();
        }
        "warm" => {
            assert_eq!(solves, 0, "warm: a recipe was solved again");
            assert!(delta.cache_hits > 0, "warm: no hit in the persistent cache");
            let cold = std::fs::read_to_string(&saved).unwrap();
            assert_eq!(
                value.to_canonical_string(),
                cold,
                "warm value differs from cold"
            );
            eprintln!("RELATION-CACHE warm: SERVED FROM DISK WITHOUT A SOLVE, value identical");
        }
        other => panic!("VAKINT_TEST_P5L_RELATION_ROLE {other}: cold or warm"),
    }
}

/// The records of a `VAKINT_DUMP_SPLIT_INPUTS` dump (`family`, `powers`,
/// `external` and `numerator` lines): family, powers, external, numerator.
fn request_records(
    text: &str,
) -> Vec<(
    &'static crate::p_integral::PFamilyData,
    Vec<i64>,
    Atom,
    Atom,
)> {
    text.split("\n---\n")
        .filter(|record| !record.trim().is_empty())
        .map(|record| {
            let field = |key: &str| {
                record
                    .lines()
                    .find_map(|line| line.strip_prefix(key)?.strip_prefix(' '))
                    .unwrap_or_else(|| panic!("a request record without {key}"))
            };
            (
                registry_family(field("family")).unwrap(),
                field("powers")
                    .split(',')
                    .map(|power| power.parse().unwrap())
                    .collect(),
                parse(field("external")),
                parse(field("numerator")),
            )
        })
        .collect()
}

/// The built-in p library (one to three loops), or the scope 1-4 package
/// named by `VAKINT_TEST_P_LIBRARY_1_4`.
fn library_engine() -> Vakint {
    let mut engine = Vakint::new().unwrap();
    if let Some(directory) = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4") {
        let library = rustred_app::propagator::load_p_library(
            rustred_app::propagator::PLibrarySource::Directory(directory.into()),
            &rustred_app::propagator::PLibraryOptions {
                verify_descent: false,
            },
        )
        .unwrap();
        engine.inject_p_library(library).unwrap();
    }
    engine
}

/// The lazy lowering of [`reduce_requests`] (one channel lowered just before
/// its reduction) gives the values, in the same order, and the errors of the
/// eager table of [`reduce_requests_eager`], request by request and on the
/// whole batch at once: the real five-loop requests of
/// `p_split_requests.txt` within the library's scope, and with
/// `VAKINT_P_REQUEST_CORPUS` a directory of dumped requests (one per file).
#[test]
fn lazy_channel_lowering_reproduces_the_eager_table() {
    let engine = library_engine();
    let settings = settings(RustRedMasterMode::Symbolic, None);
    let mut texts = vec![include_str!("../../tests/fixtures/p_split_requests.txt").to_string()];
    if let Some(directory) = std::env::var_os("VAKINT_P_REQUEST_CORPUS") {
        let mut paths = std::fs::read_dir(directory)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .collect::<Vec<_>>();
        paths.sort();
        texts.extend(
            paths
                .iter()
                .map(|path| std::fs::read_to_string(path).unwrap()),
        );
    }
    let records = texts
        .iter()
        .flat_map(|text| request_records(text))
        .filter(|(family, ..)| has_p_data(&engine, family.loops))
        .collect::<Vec<_>>();
    assert!(records.len() >= 5, "{} records in scope", records.len());
    let requests = records
        .iter()
        .map(|(family, powers, external, numerator)| PropagatorRequest {
            family,
            powers,
            canonical_numerator: numerator.as_view(),
            external: external.clone(),
        })
        .collect::<Vec<_>>();
    let same = |requests: &[PropagatorRequest]| {
        let lazy = reduce_requests_unfolded(&engine, &settings, requests, None, false);
        let eager = reduce_requests_eager(&engine, &settings, requests, None);
        match (&lazy, &eager) {
            (Ok(lazy), Ok(eager)) => {
                assert_eq!(lazy, eager);
                assert_eq!(format!("{lazy:?}"), format!("{eager:?}"));
            }
            _ => assert_eq!(
                lazy.as_ref().map_err(ToString::to_string).err(),
                eager.as_ref().map_err(ToString::to_string).err()
            ),
        }
        lazy.map_or(0, |values| values.values().map(Vec::len).sum::<usize>())
    };
    let mut values = 0;
    for request in &requests {
        values += same(std::slice::from_ref(request));
    }
    assert!(values > 0);
    assert!(same(&requests) > 0);
    // Batches whose requests share channels.
    let doubled = requests
        .iter()
        .chain(&requests)
        .cloned()
        .collect::<Vec<_>>();
    assert!(same(&doubled) > 0);
}

/// The values of one batch summed per (loops, head, scale power), the
/// spectator multiplied in.
fn summed(
    values: BTreeMap<(usize, Atom), Vec<ReducedPropagatorValue>>,
) -> BTreeMap<(usize, Atom, i64), Atom> {
    let mut sums = BTreeMap::<(usize, Atom, i64), Atom>::new();
    for ((loops, head), list) in values {
        for value in list {
            *sums
                .entry((loops, head.clone(), value.scale_power))
                .or_default() += value.value().as_ref();
        }
    }
    sums
}

/// Both fold modes give, per (loops, head, scale power), exactly the sum of
/// today's channel values, with and without factored spectators, request by
/// request and on whole batches: the real five-loop requests of
/// `p_split_requests.txt` within the library's scope (and the corpus of
/// `VAKINT_P_REQUEST_CORPUS`), and the one-loop requests with mixed
/// non-rational coefficients. A folded spectator is free of `ε` whenever
/// every channel coefficient is a polynomial in other symbols times
/// rational functions of `ε`.
#[test]
fn folded_channels_reproduce_the_unfolded_values() {
    let check = |engine: &Vakint,
                 settings: &VakintSettings,
                 requests: &[PropagatorRequest],
                 epsilon_free: bool| {
        let epsilon = epsilon(settings);
        for factored in [false, true] {
            let today =
                reduce_requests_unfolded(engine, settings, requests, None, factored).unwrap();
            let today = summed(today);
            for mode in [FoldMode::Basis, FoldMode::Coefficient] {
                let folded =
                    reduce_requests_folded(engine, settings, requests, None, mode, factored)
                        .unwrap();
                for value in folded.values().flatten() {
                    assert_eq!(value.spectator.is_some(), factored);
                    if epsilon_free && let Some(spectator) = &value.spectator {
                        assert!(!spectator.contains_symbol(epsilon), "{spectator}");
                    }
                }
                let folded = summed(folded);
                let keys = today.keys().chain(folded.keys()).collect::<BTreeSet<_>>();
                for key in keys {
                    let zero = Atom::Zero;
                    let difference = (today.get(key).unwrap_or(&zero)
                        - folded.get(key).unwrap_or(&zero))
                    .expand();
                    assert!(
                        difference.is_zero() || difference.together().expand().is_zero(),
                        "{mode:?}, factored {factored}, {key:?}: {difference}"
                    );
                }
            }
        }
    };

    let engine = library_engine();
    let settings = settings(RustRedMasterMode::Symbolic, None);
    let mut texts = vec![include_str!("../../tests/fixtures/p_split_requests.txt").to_string()];
    if let Some(directory) = std::env::var_os("VAKINT_P_REQUEST_CORPUS") {
        let mut paths = std::fs::read_dir(directory)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .collect::<Vec<_>>();
        paths.sort();
        texts.extend(
            paths
                .iter()
                .map(|path| std::fs::read_to_string(path).unwrap()),
        );
    }
    let records = texts
        .iter()
        .flat_map(|text| request_records(text))
        .filter(|(family, ..)| has_p_data(&engine, family.loops))
        .collect::<Vec<_>>();
    assert!(records.len() >= 5, "{} records in scope", records.len());
    let requests = records
        .iter()
        .map(|(family, powers, external, numerator)| PropagatorRequest {
            family,
            powers,
            canonical_numerator: numerator.as_view(),
            external: external.clone(),
        })
        .collect::<Vec<_>>();
    for request in &requests {
        check(&engine, &settings, std::slice::from_ref(request), false);
    }
    check(&engine, &settings, &requests, false);

    // One-loop requests whose coefficients mix spectator symbols and ε.
    let engine = engine_one_loop();
    let fixture = battery();
    let terms = [
        ("muvsq+ε", "bubble_12_n1"),
        ("muvsq*(1+ε)/(3-2*ε)+muvsq^2/(2-ε)", "bubble_21_n1"),
        ("ε/(1-ε)^2", "bubble_22_n1"),
        ("muvsq^2*ε", "bubble_21_n0"),
        ("𝑖*muvsq/(3-2*ε)", "bubble_12_n1"),
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
    for request in &requests {
        check(&engine, &settings, std::slice::from_ref(request), true);
    }
    check(&engine, &settings, &requests, true);
}

fn engine_one_loop() -> Vakint {
    engine()
}

/// The inline fold of `reduce_requests_folded` before helper VH lane refpath
/// split it into [`fold_requests`], kept verbatim (statistics aside) as the
/// reference of [`fold_requests_is_the_inline_fold`].
fn inline_fold_before_the_split(
    handle: &PLibraryHandle,
    settings: &VakintSettings,
    requests: &[PropagatorRequest],
) -> BTreeMap<Channel, FoldedInputs> {
    let epsilon = vk_symbol!(&settings.epsilon_symbol);
    let mut forms = BTreeMap::new();
    let folded_requests = requests
        .iter()
        .map(|request| fold_request(handle, settings, request, &mut forms, epsilon))
        .collect::<Result<Vec<_>, VakintError>>()
        .unwrap();
    let mut folded = BTreeMap::<Channel, FoldedInputs>::new();
    for request in &folded_requests {
        for (channel, by_basis) in &request.targets {
            let target = folded.entry(channel.clone()).or_default();
            target.integrals.insert(request.integral.clone());
            for (basis, inputs) in by_basis {
                match target.by_basis.get_mut(basis) {
                    None => {
                        target.by_basis.insert(basis.clone(), inputs.clone());
                    }
                    Some(sum) => {
                        for (key, value) in inputs {
                            *sum.entry(key.clone()).or_insert_with(Rational::zero) += value;
                        }
                    }
                }
            }
        }
    }
    folded
}

/// Helper VH lane refpath (reference path R2 of the numerator gate): the
/// split-off [`fold_requests`] returns exactly the folded channels of the
/// inline fold it replaced ([`inline_fold_before_the_split`]), request by
/// request and on whole batches (also with every request twice, so that
/// channels are shared), on the real requests of `p_split_requests.txt`
/// within the library's scope and on one-loop requests with mixed
/// spectators; and `reduce_requests_folded`, which now calls it, still gives
/// the values of the unfolded reduction (`folded_channels_reproduce_the_
/// unfolded_values` checks the values on the same inputs).
#[test]
fn fold_requests_is_the_inline_fold() {
    let check = |engine: &Vakint, settings: &VakintSettings, requests: &[PropagatorRequest]| {
        let handle = engine.p_library().unwrap();
        // Off and on the request memo: the second call of each pair is served
        // by the memo the first call filled.
        for _ in 0..2 {
            let split = fold_requests(&handle, settings, requests).unwrap();
            let inline = inline_fold_before_the_split(&handle, settings, requests);
            assert_eq!(split, inline);
            assert_eq!(format!("{split:?}"), format!("{inline:?}"));
            assert!(!split.is_empty());
        }
    };
    let settings = settings(RustRedMasterMode::Symbolic, None);
    let engine = library_engine();
    let text = include_str!("../../tests/fixtures/p_split_requests.txt");
    let records = request_records(text)
        .into_iter()
        .filter(|(family, ..)| has_p_data(&engine, family.loops))
        .collect::<Vec<_>>();
    assert!(records.len() >= 5, "{} records in scope", records.len());
    let requests = records
        .iter()
        .map(|(family, powers, external, numerator)| PropagatorRequest {
            family,
            powers,
            canonical_numerator: numerator.as_view(),
            external: external.clone(),
        })
        .collect::<Vec<_>>();
    for request in &requests {
        check(&engine, &settings, std::slice::from_ref(request));
    }
    check(&engine, &settings, &requests);
    let doubled = requests.iter().chain(&requests).cloned().collect::<Vec<_>>();
    check(&engine, &settings, &doubled);

    let engine = engine_one_loop();
    let fixture = battery();
    let external = parse("p(3)");
    let terms = [
        ("muvsq+ε", "bubble_12_n1"),
        ("𝑖*muvsq/(3-2*ε)", "bubble_12_n1"),
        ("ε/(1-ε)^2", "bubble_22_n1"),
    ];
    let numerators = terms.map(|(coefficient, name)| {
        parse(&format!("({coefficient})*({})", case(&fixture, name).numerator))
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
    check(&engine, &settings, &requests);
}

/// [`spectator_terms`] writes a coefficient as `Σ value · basis(ε) ·
/// monomial` with `basis` a rational function of `ε` alone and `monomial`
/// free of it; a factor that is not a rational function of `ε` alone (a Γ
/// function, a denominator mixing `ε` with another symbol) stays in the
/// monomial, so the decomposition is exact in every case.
#[test]
fn spectator_terms_split_off_the_rational_functions_of_epsilon() {
    let settings = settings(RustRedMasterMode::Symbolic, None);
    let epsilon = epsilon(&settings);
    let cases = [
        ("T*(1-ε)/(3-2*ε)+L*ε^2/((3-2*ε)*(2-ε))", true),
        ("2/3*𝑖*T*muvsq^2*(1+ε)/(3-2*ε)^2-5*𝑖*L*muvsq^2/(2-ε)", true),
        ("T*ε", true),
        ("T*exp(ε)/(3-2*ε)", false),
        ("T/(muvsq-ε)+L*(1+ε)", false),
    ];
    for (text, epsilon_free) in cases {
        let coefficient = parse(text);
        let terms = spectator_terms(&coefficient, epsilon);
        let mut sum = Atom::Zero;
        for (monomial, basis, value) in &terms {
            assert!(!value.is_zero());
            assert!(
                parts::rational(basis.as_view()).is_some() || basis.contains_symbol(epsilon),
                "{text}: basis {basis}"
            );
            if epsilon_free {
                assert!(
                    !monomial.contains_symbol(epsilon),
                    "{text}: monomial {monomial}"
                );
            }
            sum += monomial * basis * Atom::num(value.clone());
        }
        let difference = (&sum - &coefficient).expand();
        assert!(
            difference.is_zero() || difference.together().expand().is_zero(),
            "{text}: {sum}"
        );
    }
    // Two monomials, each with one basis function per distinct ε factor.
    let terms = spectator_terms(&parse("T*(1-ε)/(3-2*ε)+L/(3-2*ε)"), epsilon);
    let monomials = terms
        .iter()
        .map(|(m, ..)| m.clone())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        monomials,
        [parse("T"), parse("L")]
            .into_iter()
            .collect::<BTreeSet<_>>()
    );
}

/// vk5l WP4 like-for-like timing of `reduce_requests` on dumped requests
/// (`VAKINT_DUMP_SPLIT_INPUTS` records, e.g. the beta team's physical
/// requests), skipped unless `VAKINT_TEST_P_REQUEST_FILE` names a file: the
/// library of [`library_engine`] (`VAKINT_TEST_P_LIBRARY_1_4`), master mode
/// `VAKINT_TEST_P_REQUEST_MODE` (`symbolic` default, `retain`), one thread
/// unless `VAKINT_TEST_P_REQUEST_CORES`; prints the cold call with its
/// split / lower / prepare / library / render seconds, the finalization, and
/// a warm second call in the same process.
#[test]
fn p_request_file_timing() {
    let Some(file) = std::env::var_os("VAKINT_TEST_P_REQUEST_FILE") else {
        eprintln!("p_request_file_timing: VAKINT_TEST_P_REQUEST_FILE unset, skipped");
        return;
    };
    let started = std::time::Instant::now();
    let engine = library_engine();
    let mode = match std::env::var("VAKINT_TEST_P_REQUEST_MODE").as_deref() {
        Ok("retain") => RustRedMasterMode::Retain,
        _ => RustRedMasterMode::Symbolic,
    };
    let mut settings = settings(mode, None);
    settings.n_cores = std::env::var("VAKINT_TEST_P_REQUEST_CORES")
        .ok()
        .and_then(|cores| cores.parse().ok())
        .unwrap_or(1);
    let text = std::fs::read_to_string(&file).unwrap();
    // The beta team's requests write ε as gammalooprs' symbol, and their
    // settings name it (MSbar measure, five ε terms): as in production, so
    // that ε-rational coefficients fold into the inputs.
    if let Ok(epsilon) = std::env::var("VAKINT_TEST_P_REQUEST_EPSILON") {
        settings.epsilon_symbol = epsilon;
    } else if text.contains("gammalooprs::{}::ε") {
        settings.epsilon_symbol = "gammalooprs::ε".into();
        settings.integral_normalization_factor = crate::LoopNormalizationFactor::MSbar;
        settings.number_of_terms_in_epsilon_expansion = 5;
    }
    let records = request_records(&text);
    let requests = records
        .iter()
        .map(|(family, powers, external, numerator)| PropagatorRequest {
            family,
            powers,
            canonical_numerator: numerator.as_view(),
            external: external.clone(),
        })
        .collect::<Vec<_>>();
    let library_load = {
        let loading = std::time::Instant::now();
        engine.p_library().unwrap();
        loading.elapsed().as_secs_f64()
    };
    eprintln!(
        "P-TIMING setup {:.3} s (library {library_load:.3} s), {} requests, fold {:?}, mode {mode:?}, {} cores, epsilon {}",
        started.elapsed().as_secs_f64(),
        requests.len(),
        super::fold_mode(),
        settings.n_cores,
        settings.epsilon_symbol
    );
    let names = request_names(&requests).unwrap();
    for pass in ["cold", "warm"] {
        let calling = std::time::Instant::now();
        let reduced = reduce_requests(&engine, &settings, &requests, None).unwrap();
        let call = calling.elapsed().as_secs_f64();
        let stats = super::last_stats();
        let finalizing = std::time::Instant::now();
        let value = finalize_reduced(&engine, &settings, reduced, &names, mode).unwrap();
        let finalize = finalizing.elapsed().as_secs_f64();
        eprintln!(
            "P-TIMING {pass}: reduce_requests {call:.3} s = split {:.3} + lower/fold {:.3} + \
             coefficients {:.3} + library {:.3} + render {:.3} (channels {}, reductions {}, pairs {}, \
             memo hits {}); finalize {finalize:.3} s; {} terms",
            stats.split_seconds,
            stats.lower_seconds,
            stats.prepare_seconds,
            stats.library_seconds,
            stats.render_seconds,
            stats.channels,
            stats.reductions,
            stats.pairs,
            stats.memo_hits,
            drill_term_count(&value.expand())
        );
    }
}

/// vk5l WP4: a batch met again is served by the request memo (no split, no
/// lowering: every request a hit) and reduces to the same values; with the
/// memo bound at zero bytes nothing is kept.
#[test]
fn a_repeated_batch_reuses_the_request_memo() {
    let engine = library_engine();
    let settings = settings(RustRedMasterMode::Symbolic, None);
    let text = include_str!("../../tests/fixtures/p_split_requests.txt");
    let records = request_records(text)
        .into_iter()
        .filter(|(family, ..)| has_p_data(&engine, family.loops))
        .collect::<Vec<_>>();
    let requests = records
        .iter()
        .map(|(family, powers, external, numerator)| PropagatorRequest {
            family,
            powers,
            canonical_numerator: numerator.as_view(),
            external: external.clone(),
        })
        .collect::<Vec<_>>();
    assert!(!requests.is_empty());
    if super::fold_mode() == super::FoldMode::Off {
        return; // the memo serves the folded path
    }
    let first = reduce_requests(&engine, &settings, &requests, None).unwrap();
    assert_eq!(super::last_stats().memo_hits, 0);
    let second = reduce_requests(&engine, &settings, &requests, None).unwrap();
    assert_eq!(super::last_stats().memo_hits, requests.len());
    assert_eq!(super::last_stats().split_seconds, 0.0);
    assert_eq!(format!("{first:?}"), format!("{second:?}"));
    // A zero bound keeps nothing.
    let fresh = library_engine();
    let mut bounded = settings.clone();
    bounded.ordered_cache.max_resident_bytes = Some(0);
    reduce_requests(&fresh, &bounded, &requests, None).unwrap();
    reduce_requests(&fresh, &bounded, &requests, None).unwrap();
    assert_eq!(super::last_stats().memo_hits, 0);
}

/// vk5l WP4 for the p master catalogue (lane vk5l-mastermap): write every
/// terminal of the package `VAKINT_TEST_P_TERMINAL_PACKAGE` as vakint renders
/// it (`Vakint::p_terminal_listing`), one JSON object per line, to
/// `VAKINT_TEST_P_TERMINAL_OUT`; skipped unless both are set.
#[test]
fn p_terminal_export() {
    let (Some(directory), Some(out)) = (
        std::env::var_os("VAKINT_TEST_P_TERMINAL_PACKAGE"),
        std::env::var_os("VAKINT_TEST_P_TERMINAL_OUT"),
    ) else {
        eprintln!("p_terminal_export: environment unset, skipped");
        return;
    };
    let mut engine = Vakint::new().unwrap();
    engine
        .load_p_library_directory(std::path::Path::new(&directory), false)
        .unwrap();
    let records = engine.p_terminal_listing().unwrap();
    let mut text = String::new();
    for record in &records {
        text += &serde_json::json!({
            "schema": "vakint.p-terminal.v1",
            "id": record.id,
            "atom": record.atom,
            "node": record.node,
            "key": record.key,
            "boxed": record.boxed,
            "owner_class": record.owner_class,
            "lane": record.lane,
            "presentation": record.presentation,
            "loops": record.loops,
            "family": record.family,
            "slots": record.slots,
        })
        .to_string();
        text.push('\n');
    }
    std::fs::write(&out, text).unwrap();
    eprintln!(
        "p_terminal_export: {} terminals -> {}",
        records.len(),
        out.to_string_lossy()
    );
}

/// vk5l-physrun dashboard post-processing (env-gated):
/// `VAKINT_TEST_PP_RETAIN=<retain .sym>` -> `<VAKINT_TEST_PP_OUT>/result_exact_d_coef.txt`: the Retain value with
/// `gammalooprs::ε = (4 - d)/2`, grouped by master monomial (the product of its `PGamma(..)` / `gamma(..)`
/// factors) as `Σ coef(c) · master` (`c` with collect_num); `VAKINT_TEST_PP_SYMBOLIC=<symbolic .sym>` ->
/// `<OUT>/eps_orders_collected.json`: every order from `VAKINT_TEST_PP_MIN_ORDER` (default -6) to the highest
/// present, each `collect_num().collect_factors()`, an explicit "0" where it cancels.
#[test]
fn showcase_postprocess() {
    let Some(out) = std::env::var_os("VAKINT_TEST_PP_OUT") else {
        eprintln!("showcase_postprocess: environment unset, skipped");
        return;
    };
    let out = std::path::PathBuf::from(out);
    std::fs::create_dir_all(&out).unwrap();
    let _ = Vakint::new().unwrap();
    let epsilon = crate::utils::vakint_macros::vk_symbol!("gammalooprs::ε");
    if let Some(path) = std::env::var_os("VAKINT_TEST_PP_RETAIN") {
        let started = std::time::Instant::now();
        let value = parse(&std::fs::read_to_string(&path).unwrap());
        let d = parse("d");
        let rule = symbolica::id::Replacement::new(
            Atom::var(epsilon).to_pattern(),
            ((Atom::num(4) - &d) / Atom::num(2)).to_pattern(),
        );
        let in_d = value.replace_multiple(&[rule]).expand();
        let is_master = |factor: symbolica::atom::AtomView| -> bool {
            let base = match factor {
                symbolica::atom::AtomView::Pow(power) => power.get_base_exp().0,
                other => other,
            };
            match base {
                symbolica::atom::AtomView::Fun(function) => {
                    let symbol = function.get_symbol();
                    let name = symbol.get_name();
                    name.ends_with("::gamma")
                        || name == "gamma"
                        || name.ends_with("PGamma")
                        || name.ends_with("PMaster")
                        || name.ends_with("PResidual")
                        || name.ends_with("PG")
                }
                _ => false,
            }
        };
        let mut order = Vec::<Atom>::new();
        let mut groups = std::collections::HashMap::<Atom, Vec<Atom>>::new();
        let terms = match in_d.as_view() {
            symbolica::atom::AtomView::Add(sum) => {
                sum.iter().map(|t| t.to_owned()).collect::<Vec<_>>()
            }
            _ => vec![in_d.clone()],
        };
        for term in terms {
            let factors = match term.as_view() {
                symbolica::atom::AtomView::Mul(product) => {
                    product.iter().map(|f| f.to_owned()).collect::<Vec<_>>()
                }
                _ => vec![term.clone()],
            };
            let (masters, rest): (Vec<_>, Vec<_>) =
                factors.into_iter().partition(|f| is_master(f.as_view()));
            let master = if masters.is_empty() {
                Atom::num(1)
            } else {
                Atom::mul_many(masters.iter().map(|a| a.as_view()))
            };
            let coefficient = if rest.is_empty() {
                Atom::num(1)
            } else {
                Atom::mul_many(rest.iter().map(|a| a.as_view()))
            };
            if !groups.contains_key(&master) {
                order.push(master.clone());
            }
            groups.entry(master).or_default().push(coefficient);
        }
        let coef = crate::utils::vakint_macros::vk_symbol!("coef");
        let mut result = Atom::Zero;
        for master in &order {
            let coefficient =
                Atom::add_many(groups[master].iter().map(|a| a.as_view())).collect_num();
            result += symbolica::function!(coef, coefficient) * master;
        }
        let text = result.to_canonical_string();
        std::fs::write(out.join("result_exact_d_coef.txt"), &text).unwrap();
        // The same with each master monomial as a named symbol Master(i) and a
        // legend Master(i) -> its closed form (no closed form inline).
        let master_symbol = crate::utils::vakint_macros::vk_symbol!("Master");
        let mut named = Atom::Zero;
        let mut legend = Vec::new();
        for (index, master) in order.iter().enumerate() {
            let coefficient =
                Atom::add_many(groups[master].iter().map(|a| a.as_view())).collect_num();
            named += symbolica::function!(coef, coefficient)
                * symbolica::function!(master_symbol, Atom::num(index as i64 + 1));
            let mut pgamma = Vec::new();
            let mut gammas = Vec::new();
            for factor in match master.as_view() {
                symbolica::atom::AtomView::Mul(product) => {
                    product.iter().map(|f| f.to_owned()).collect::<Vec<_>>()
                }
                _ => vec![master.clone()],
            } {
                let text = factor.to_canonical_string();
                if text.contains("PGamma") {
                    pgamma.push(text)
                } else {
                    gammas.push(text)
                }
            }
            legend.push(serde_json::json!({"symbol": format!("vakint::Master({})", index + 1),
                "closed_form": master.to_canonical_string(), "pgamma_tags": pgamma, "gamma_functions": gammas}));
        }
        std::fs::write(
            out.join("result_exact_masters.txt"),
            named.to_canonical_string(),
        )
        .unwrap();
        std::fs::write(
            out.join("result_exact_masters_legend.json"),
            serde_json::to_string_pretty(&legend).unwrap(),
        )
        .unwrap();
        std::fs::write(
            out.join("result_exact_d_coef.json"),
            serde_json::json!({"masters": order.len(), "terms_in_d": groups.values().map(|g| g.len()).sum::<usize>(),
                               "bytes": text.len(), "seconds": started.elapsed().as_secs_f64(),
                               "masters_list": order.iter().map(|m| m.to_canonical_string()).collect::<Vec<_>>()}).to_string(),
        )
        .unwrap();
        eprintln!(
            "showcase_postprocess: {} master monomials, {} bytes",
            order.len(),
            text.len()
        );
    }
    if let Some(path) = std::env::var_os("VAKINT_TEST_PP_SYMBOLIC") {
        let value = parse(&std::fs::read_to_string(&path).unwrap()).expand();
        let mut orders = std::collections::BTreeMap::<i64, Atom>::new();
        for (power, coefficient) in value.coefficient_list::<i8>(&[Atom::var(epsilon)]) {
            let order = if power == Atom::num(1) {
                0
            } else if power == Atom::var(epsilon) {
                1
            } else {
                match power.as_view() {
                    symbolica::atom::AtomView::Pow(p) => {
                        crate::get_integer_from_atom(p.get_base_exp().1).unwrap()
                    }
                    _ => panic!("not a power of epsilon: {power}"),
                }
            };
            *orders.entry(order).or_insert(Atom::Zero) += coefficient;
        }
        let min = std::env::var("VAKINT_TEST_PP_MIN_ORDER")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(-6);
        let max = orders.keys().copied().max().unwrap_or(min);
        let mut json = serde_json::Map::new();
        for order in min.min(max)..=max {
            let text = match orders.get(&order) {
                Some(coefficient) => {
                    let collected = coefficient.expand().collect_num().collect_factors();
                    if collected.is_zero() {
                        "0".to_owned()
                    } else {
                        collected.to_canonical_string()
                    }
                }
                None => "0".to_owned(),
            };
            json.insert(order.to_string(), serde_json::Value::String(text));
        }
        std::fs::write(
            out.join("eps_orders_collected.json"),
            serde_json::to_string_pretty(&json).unwrap(),
        )
        .unwrap();
        eprintln!("showcase_postprocess: orders {min}..{max}");
    }
}

/// vk5l-physrun: the Laurent orders of a Symbolica canonical-string result
/// (`VAKINT_TEST_EPS_IN`, e.g. sum.symbolic.sym) in `gammalooprs::ε`, written
/// to `VAKINT_TEST_EPS_OUT` as JSON {"orders": {"<n>": "<canonical string>"},
/// "symbols": [non-ε symbols present]}.
#[test]
fn showcase_eps_orders() {
    let (Some(input), Some(output)) = (
        std::env::var_os("VAKINT_TEST_EPS_IN"),
        std::env::var_os("VAKINT_TEST_EPS_OUT"),
    ) else {
        eprintln!("showcase_eps_orders: environment unset, skipped");
        return;
    };
    let _ = Vakint::new().unwrap();
    let text = std::fs::read_to_string(&input).unwrap();
    let value = parse(&text).expand();
    let epsilon = crate::utils::vakint_macros::vk_symbol!("gammalooprs::ε");
    let mut orders = serde_json::Map::new();
    let mut symbols = std::collections::BTreeSet::new();
    for (power, coefficient) in value.coefficient_list::<i8>(&[Atom::var(epsilon)]) {
        let order = if power == Atom::num(1) {
            0
        } else if power == Atom::var(epsilon) {
            1
        } else {
            match power.as_view() {
                symbolica::atom::AtomView::Pow(p) => {
                    let (_, e) = p.get_base_exp();
                    crate::get_integer_from_atom(e).unwrap_or(i64::MAX)
                }
                _ => i64::MAX,
            }
        };
        coefficient.as_view().visitor(&mut |part| {
            if let symbolica::atom::AtomView::Var(v) = part {
                symbols.insert(v.get_symbol().get_name().to_string());
            }
            if let symbolica::atom::AtomView::Fun(f) = part {
                symbols.insert(format!("{}()", f.get_symbol().get_name()));
            }
            true
        });
        let key = order.to_string();
        let previous = orders.remove(&key).and_then(|v| v.as_str().map(parse));
        let total = match previous {
            Some(previous) => (previous + coefficient).expand(),
            None => coefficient,
        };
        orders.insert(key, serde_json::Value::String(total.to_canonical_string()));
    }
    let json = serde_json::json!({"orders": orders, "symbols": symbols});
    std::fs::write(&output, serde_json::to_string_pretty(&json).unwrap()).unwrap();
    eprintln!("showcase_eps_orders: {} orders", orders_len(&json));
}

fn orders_len(json: &serde_json::Value) -> usize {
    json["orders"].as_object().map(|o| o.len()).unwrap_or(0)
}

/// vk5l-physrun: a WHOLE GRAPH through [`MassiveClosure::evaluate_direct_sum`]:
/// all closure records of `VAKINT_TEST_CLOSURE_FILES` (`:`-separated) lowered
/// directly and summed with one key plan per closing line, reduced once per
/// mode (`VAKINT_TEST_CLOSURE_MODES`, default `retain`); writes
/// `sum.<mode>.sym` and `sum.<mode>.json` to `VAKINT_TEST_CLOSURE_OUT`. With
/// `VAKINT_TEST_CLOSURE_SUM_CHECK=<dir>` the per-closure results
/// `<dir>/<closure>/<closure>.<mode>.sym` are added and compared (exact).
/// vk4l-beta-slow tests: `topo(Π prop(id, edge, q, m, power))` with the power
/// of its `n`-th prop factor whose momentum contains a loop momentum raised by one.
fn raise_line_power(topology: &Atom, n: usize) -> Atom {
    let symbolica::atom::AtomView::Fun(topo) = topology.as_view() else {
        panic!("topology is not topo(..)");
    };
    let inner = topo.iter().next().unwrap();
    let mut factors = match inner {
        symbolica::atom::AtomView::Mul(product) => {
            product.iter().map(|f| f.to_owned()).collect::<Vec<_>>()
        }
        other => vec![other.to_owned()],
    };
    let candidates = factors
        .iter()
        .enumerate()
        .filter(|(_, factor)| {
            matches!(factor.as_view(), symbolica::atom::AtomView::Fun(f)
                if f.get_symbol() == crate::symbols::S.prop && factor.contains_symbol(crate::symbols::S.k))
        })
        .map(|(index, _)| index)
        .collect::<Vec<_>>();
    let index = candidates[n % candidates.len()];
    let symbolica::atom::AtomView::Fun(prop) = factors[index].as_view() else {
        unreachable!()
    };
    let mut args = prop.iter().map(|a| a.to_owned()).collect::<Vec<_>>();
    let last = args.len() - 1;
    args[last] = (&args[last] + Atom::num(1)).expand();
    factors[index] = symbolica::atom::FunctionBuilder::new(prop.get_symbol())
        .add_args(&args)
        .finish();
    symbolica::atom::FunctionBuilder::new(topo.get_symbol())
        .add_arg(Atom::mul_many(factors.iter()))
        .finish()
}

/// vk4l-beta-slow: one `pint6.massive-closure-request.v1` file (plain or
/// `.zst`) as (closure, direct input with the identity routing, settings).
fn closure_file_input(
    path: &std::path::Path,
) -> (
    crate::p_integral::MassiveClosure,
    crate::p_integral::DirectClosureInput,
    crate::VakintSettings,
) {
    let name = path.display().to_string();
    let text = if name.ends_with(".zst") {
        let raw = std::process::Command::new("zstd")
            .args(["-qdc"])
            .arg(path)
            .output()
            .unwrap();
        String::from_utf8(raw.stdout).unwrap()
    } else {
        std::fs::read_to_string(path).unwrap()
    };
    let field = |key: &str| {
        text.lines()
            .find_map(|line| line.strip_prefix(key)?.strip_prefix(' '))
            .unwrap_or_else(|| panic!("{name}: no {key}"))
            .to_owned()
    };
    let mut settings = settings(RustRedMasterMode::Retain, None);
    settings.epsilon_symbol = "gammalooprs::ε".into();
    settings.mu_r_sq_symbol = "gammalooprs::μᵣ²".into();
    settings.integral_normalization_factor = crate::LoopNormalizationFactor::MSbar;
    settings.number_of_terms_in_epsilon_expansion = field("epsilon_terms").parse().unwrap();
    settings.use_dot_product_notation = false;
    settings.n_cores = 1;
    let closure = crate::p_integral::MassiveClosure {
        momentum: parse(&field("closure_momentum")),
        mass_squared: parse(&field("closure_mass_squared")),
        power: field("closure_power").parse().unwrap(),
    };
    let input = parse(&field("input"));
    let (topology, numerator) = match input.as_view() {
        symbolica::atom::AtomView::Mul(product) => {
            let mut topology = Atom::num(1);
            let mut numerator = Atom::num(1);
            for factor in product.iter() {
                match factor {
                    symbolica::atom::AtomView::Fun(function)
                        if function.get_symbol() == crate::symbols::S.topo =>
                    {
                        topology = factor.to_owned()
                    }
                    _ => numerator = numerator * factor.to_owned(),
                }
            }
            (topology, numerator)
        }
        symbolica::atom::AtomView::Fun(_) => (input.clone(), Atom::num(1)),
        _ => panic!("{name}: input is not numerator * topo"),
    };
    let loops: i64 = field("p_loops").parse().unwrap();
    let mut routing = (1..=loops)
        .map(|i| (parse(&format!("k({i})")), parse(&format!("k({i})"))))
        .collect::<Vec<_>>();
    routing.push((closure.momentum.clone(), closure.momentum.clone()));
    (
        closure,
        crate::p_integral::DirectClosureInput {
            topology,
            routing,
            numerator,
        },
        settings,
    )
}

/// vk4l-beta-slow (R5 of the producer design): `pint6.remainder-boundary.v1`
/// records (header paths in `VAKINT_TEST_BOUNDARY_RECORDS`) read with the
/// streaming reader; every pattern with a `ref_closure` must lower to the
/// folded input of that v1 closure record (`closure-<pid>-<seq:06>.txt[.zst]`
/// beside the header); with `VAKINT_TEST_BOUNDARY_EVALUATE=1` the boundary
/// value must equal the sum of its reference closures (when every pattern has
/// one). Five-loop records need `VAKINT_TEST_P5L_DRILL_PACKAGE`.
#[test]
fn remainder_boundary_records() {
    let Some(files) = std::env::var_os("VAKINT_TEST_BOUNDARY_RECORDS") else {
        eprintln!("remainder_boundary_records: VAKINT_TEST_BOUNDARY_RECORDS unset, skipped");
        return;
    };
    let mut engine = Vakint::new().unwrap();
    if let Some(package) = std::env::var_os("VAKINT_TEST_P5L_DRILL_PACKAGE") {
        engine
            .load_p_library_directory(std::path::Path::new(&package), false)
            .unwrap();
        engine.p_library().unwrap();
    }
    let evaluate = std::env::var("VAKINT_TEST_BOUNDARY_EVALUATE").as_deref() == Ok("1");
    for header in std::env::split_paths(&files) {
        let name = header.file_name().unwrap().to_string_lossy().into_owned();
        let pid = name
            .strip_prefix("boundary-")
            .and_then(|rest| rest.split('-').next())
            .unwrap_or_else(|| panic!("{name}: not boundary-<pid>-<seq>.txt"))
            .to_owned();
        let reading = std::time::Instant::now();
        let record = crate::p_integral::read_boundary_record(&header).unwrap();
        let read_seconds = reading.elapsed().as_secs_f64();
        let directory = header.parent().unwrap();
        let reference_file = |seq: &str| -> std::path::PathBuf {
            let seq: usize = seq.parse().unwrap();
            let plain = directory.join(format!("closure-{pid}-{seq:06}.txt"));
            if plain.exists() {
                plain
            } else {
                directory.join(format!("closure-{pid}-{seq:06}.txt.zst"))
            }
        };
        let mut checked = 0;
        let mut references = Vec::new();
        for (index, reference) in record.ref_closures.iter().enumerate() {
            let Some(reference) = reference else { continue };
            let (closure, direct, settings) = closure_file_input(&reference_file(reference));
            let expected = closure
                .lower_direct(&engine, &settings, &direct)
                .unwrap()
                .unwrap_or_else(|| {
                    panic!("{name}: the direct path declines ref_closure {reference}")
                });
            let lowered = record
                .boundary
                .lower_pattern(&engine, &settings, index)
                .unwrap()
                .unwrap_or_else(|reason| panic!("{name}: pattern {index} declined: {reason}"));
            assert_eq!(
                lowered.fingerprint(),
                expected.fingerprint(),
                "{name}: pattern {index} against ref_closure {reference}"
            );
            references.push((closure, direct, settings));
            checked += 1;
        }
        eprintln!(
            "BOUNDARY_RECORD {name}: read {read_seconds:.1} s, {} blocks, {} patterns, {checked} lowered identically to their ref_closure",
            record.boundary.blocks.len(),
            record.boundary.patterns.len()
        );
        if evaluate && checked == record.boundary.patterns.len() && checked > 0 {
            let settings = references[0].2.clone();
            let (value, stats) = record
                .boundary
                .evaluate_with_stats(
                    &engine,
                    &settings,
                    crate::p_integral::BoundaryOptions { strict: true },
                )
                .unwrap();
            let inputs = references
                .into_iter()
                .map(|(closure, direct, _)| (closure, direct))
                .collect::<Vec<_>>();
            let expected =
                crate::p_integral::MassiveClosure::evaluate_direct_sum(&engine, &settings, &inputs)
                    .unwrap();
            assert_eq!(value, expected.expand(), "{name}: value");
            eprintln!("BOUNDARY_RECORD {name}: value identical; {stats:?}");
        }
    }
}

/// vk4l-beta-slow: the remainder boundary end to end on closure request files
/// (`pint6.massive-closure-request.v1`, plain or `.zst`; path list in
/// `VAKINT_TEST_BOUNDARY_FILES`; a five-loop package in
/// `VAKINT_TEST_P5L_DRILL_PACKAGE` when the closures are five-loop). For every
/// closure: (1) the boundary of the direct input lowers to exactly the folded
/// RustRed input of `lower_direct`; (2) the numerator split into two blocks
/// (two patterns on the same topology) lowers to the same key plan in one
/// class; (3) with `VAKINT_TEST_BOUNDARY_EVALUATE=1`, the boundary value equals
/// `evaluate_direct_sum`. Strict mode: a decline fails the test.
#[test]
fn remainder_boundary_closure_files() {
    let Some(files) = std::env::var_os("VAKINT_TEST_BOUNDARY_FILES") else {
        eprintln!("remainder_boundary_closure_files: VAKINT_TEST_BOUNDARY_FILES unset, skipped");
        return;
    };
    let mut engine = Vakint::new().unwrap();
    if let Some(package) = std::env::var_os("VAKINT_TEST_P5L_DRILL_PACKAGE") {
        engine
            .load_p_library_directory(std::path::Path::new(&package), false)
            .unwrap();
        engine.p_library().unwrap();
    }
    let evaluate = std::env::var("VAKINT_TEST_BOUNDARY_EVALUATE").as_deref() == Ok("1");
    let strict = crate::p_integral::BoundaryOptions { strict: true };
    let (mut checked, mut patterns_total) = (0usize, 0usize);
    for file in std::env::split_paths(&files) {
        let name = file.file_name().unwrap().to_string_lossy().into_owned();
        let text = if name.ends_with(".zst") {
            let raw = std::process::Command::new("zstd")
                .args(["-qdc"])
                .arg(&file)
                .output()
                .unwrap();
            String::from_utf8(raw.stdout).unwrap()
        } else {
            std::fs::read_to_string(&file).unwrap()
        };
        let field = |key: &str| {
            text.lines()
                .find_map(|line| line.strip_prefix(key)?.strip_prefix(' '))
                .unwrap_or_else(|| panic!("{name}: no {key}"))
                .to_owned()
        };
        let mut settings = settings(RustRedMasterMode::Retain, None);
        settings.epsilon_symbol = "gammalooprs::ε".into();
        settings.mu_r_sq_symbol = "gammalooprs::μᵣ²".into();
        settings.integral_normalization_factor = crate::LoopNormalizationFactor::MSbar;
        settings.number_of_terms_in_epsilon_expansion = field("epsilon_terms").parse().unwrap();
        settings.use_dot_product_notation = false;
        settings.n_cores = 1;
        let closure = crate::p_integral::MassiveClosure {
            momentum: parse(&field("closure_momentum")),
            mass_squared: parse(&field("closure_mass_squared")),
            power: field("closure_power").parse().unwrap(),
        };
        let input = parse(&field("input"));
        let symbolica::atom::AtomView::Mul(product) = input.as_view() else {
            panic!("{name}: input is not numerator * topo");
        };
        let mut topology = Atom::num(1);
        let mut numerator = Atom::num(1);
        for factor in product.iter() {
            match factor {
                symbolica::atom::AtomView::Fun(function)
                    if function.get_symbol() == crate::symbols::S.topo =>
                {
                    topology = factor.to_owned()
                }
                _ => numerator = numerator * factor.to_owned(),
            }
        }
        let loops: i64 = field("p_loops").parse().unwrap();
        let mut routing = (1..=loops)
            .map(|i| (parse(&format!("k({i})")), parse(&format!("k({i})"))))
            .collect::<Vec<_>>();
        routing.push((closure.momentum.clone(), closure.momentum.clone()));
        let direct = crate::p_integral::DirectClosureInput {
            topology,
            routing,
            numerator,
        };
        let started = std::time::Instant::now();
        let reference = closure
            .lower_direct(&engine, &settings, &direct)
            .unwrap()
            .unwrap_or_else(|| panic!("{name}: the direct path declines this closure"));
        let direct_seconds = started.elapsed().as_secs_f64();
        // (1) one block, one pattern.
        let boundary = crate::p_integral::RemainderBoundary::from_direct(&closure, &direct);
        let started = std::time::Instant::now();
        let (classes, declined) = boundary.lower(&engine, &settings, strict).unwrap();
        let boundary_seconds = started.elapsed().as_secs_f64();
        assert!(declined.is_empty());
        assert_eq!(classes.len(), 1, "{name}");
        assert_eq!(
            classes[0].fingerprint(),
            reference.fingerprint(),
            "{name}: (1)"
        );
        // (2) two blocks on two patterns of the same topology.
        let terms = match direct.numerator.as_view() {
            symbolica::atom::AtomView::Add(sum) => {
                sum.iter().map(|t| t.to_owned()).collect::<Vec<_>>()
            }
            other => vec![other.to_owned()],
        };
        let half = terms.len() / 2;
        let (first, second) = terms.split_at(half.max(1).min(terms.len()));
        let blocks = vec![Atom::add_many(first.iter()), Atom::add_many(second.iter())];
        let pattern = |block: usize| crate::p_integral::BoundaryPattern {
            topology: direct.topology.clone(),
            closure_power: closure.power,
            terms: vec![(block, Atom::num(1))],
        };
        let split = crate::p_integral::RemainderBoundary::new(
            name.clone(),
            closure.momentum.clone(),
            closure.mass_squared.clone(),
            direct.routing.clone(),
            blocks,
            vec![pattern(0), pattern(1)],
        )
        .unwrap();
        let (classes, declined) = split.lower(&engine, &settings, strict).unwrap();
        assert!(declined.is_empty());
        assert_eq!(classes.len(), 1, "{name}: same closing line, one key plan");
        assert_eq!(
            classes[0].fingerprint(),
            reference.fingerprint(),
            "{name}: (2)"
        );
        patterns_total += split.patterns.len();
        // (4) block reuse: three patterns over the two blocks, the third on the
        // same cut with the power of its first line raised by one; compared
        // with lowering every assembled pattern on its own and absorbing.
        let raised = raise_line_power(&direct.topology, 0);
        let mut reused = split.clone();
        reused.patterns.push(crate::p_integral::BoundaryPattern {
            topology: raised,
            closure_power: closure.power,
            terms: vec![(0, Atom::num(2)), (1, Atom::num(-1))],
        });
        let started = std::time::Instant::now();
        let (classes, declined, (reads, hits)) =
            reused.lower_counted(&engine, &settings, strict).unwrap();
        let reuse_seconds = started.elapsed().as_secs_f64();
        assert!(declined.is_empty());
        let started = std::time::Instant::now();
        let mut expected: Vec<crate::p_integral::LoweredClosure> = Vec::new();
        for pattern in &reused.patterns {
            let lowered = reused
                .closure(pattern)
                .lower_direct(&engine, &settings, &reused.direct_input(pattern))
                .unwrap()
                .unwrap();
            let mut lowered = Some(lowered);
            for class in expected.iter_mut() {
                match class.absorb(lowered.take().unwrap()) {
                    Ok(()) => break,
                    Err(back) => lowered = Some(back),
                }
            }
            if let Some(lowered) = lowered {
                expected.push(lowered);
            }
        }
        let separate_seconds = started.elapsed().as_secs_f64();
        assert_eq!(classes.len(), expected.len(), "{name}: (4) classes");
        for (class, expected) in classes.iter().zip(&expected) {
            assert_eq!(class.fingerprint(), expected.fingerprint(), "{name}: (4)");
        }
        eprintln!(
            "BOUNDARY {name} reuse: 3 patterns, {reads} block reads, {hits} reuses, {reuse_seconds:.3} s (separate {separate_seconds:.3} s), identical"
        );
        if evaluate {
            let started = std::time::Instant::now();
            let value = boundary.evaluate(&engine, &settings, strict).unwrap();
            let boundary_value_seconds = started.elapsed().as_secs_f64();
            let started = std::time::Instant::now();
            let expected = crate::p_integral::MassiveClosure::evaluate_direct_sum(
                &engine,
                &settings,
                &[(closure.clone(), direct.clone())],
            )
            .unwrap();
            assert_eq!(value, expected.expand(), "{name}: (3)");
            eprintln!(
                "BOUNDARY {name} evaluate {boundary_value_seconds:.3} s, direct sum {:.3} s, identical",
                started.elapsed().as_secs_f64()
            );
        }
        eprintln!(
            "BOUNDARY {name} lowered identically: direct {direct_seconds:.3} s, boundary {boundary_seconds:.3} s, {} numerator terms",
            terms.len()
        );
        checked += 1;
    }
    eprintln!("BOUNDARY {checked} closures checked ({patterns_total} split patterns)");
}

/// vk4l-beta-slow benchmark (env-gated): one closure file
/// (`VAKINT_TEST_SYNTH_FILE`, `.zst` or plain), its numerator multiplied by
/// `x(1) + ... + x(m)` (fresh loop-free spectators) and expanded to about
/// `VAKINT_TEST_SYNTH_MONOMIALS` monomials: one block shared by
/// `VAKINT_TEST_SYNTH_PATTERNS` patterns (line i's power raised by one in
/// pattern i). Prints the block size, `lower_counted` seconds, block reads /
/// reuses and the peak RSS. Run with and without `VAKINT_P_BOUNDARY_BLOCKS=off`.
#[test]
fn remainder_boundary_synthetic_blocks() {
    let Some(file) = std::env::var_os("VAKINT_TEST_SYNTH_FILE") else {
        eprintln!("remainder_boundary_synthetic_blocks: VAKINT_TEST_SYNTH_FILE unset, skipped");
        return;
    };
    let target: usize = std::env::var("VAKINT_TEST_SYNTH_MONOMIALS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(1_000_000);
    let patterns: usize = std::env::var("VAKINT_TEST_SYNTH_PATTERNS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(4);
    let peak = || -> u64 {
        std::fs::read_to_string("/proc/self/status")
            .ok()
            .and_then(|status| {
                status
                    .lines()
                    .find(|line| line.starts_with("VmHWM:"))
                    .and_then(|line| line.split_whitespace().nth(1)?.parse::<u64>().ok())
            })
            .unwrap_or(0)
    };
    let mut engine = Vakint::new().unwrap();
    if let Some(package) = std::env::var_os("VAKINT_TEST_P5L_DRILL_PACKAGE") {
        engine
            .load_p_library_directory(std::path::Path::new(&package), false)
            .unwrap();
        engine.p_library().unwrap();
    }
    let path = std::path::PathBuf::from(&file);
    let text = if path.to_string_lossy().ends_with(".zst") {
        let raw = std::process::Command::new("zstd")
            .args(["-qdc"])
            .arg(&path)
            .output()
            .unwrap();
        String::from_utf8(raw.stdout).unwrap()
    } else {
        std::fs::read_to_string(&path).unwrap()
    };
    let field = |key: &str| {
        text.lines()
            .find_map(|line| line.strip_prefix(key)?.strip_prefix(' '))
            .unwrap_or_else(|| panic!("no {key}"))
            .to_owned()
    };
    let mut settings = settings(RustRedMasterMode::Retain, None);
    settings.epsilon_symbol = "gammalooprs::ε".into();
    settings.mu_r_sq_symbol = "gammalooprs::μᵣ²".into();
    settings.integral_normalization_factor = crate::LoopNormalizationFactor::MSbar;
    settings.number_of_terms_in_epsilon_expansion = field("epsilon_terms").parse().unwrap();
    settings.use_dot_product_notation = false;
    settings.n_cores = 1;
    let closure = crate::p_integral::MassiveClosure {
        momentum: parse(&field("closure_momentum")),
        mass_squared: parse(&field("closure_mass_squared")),
        power: field("closure_power").parse().unwrap(),
    };
    let input = parse(&field("input"));
    let symbolica::atom::AtomView::Mul(product) = input.as_view() else {
        panic!("input is not numerator * topo");
    };
    let mut topology = Atom::num(1);
    let mut numerator = Atom::num(1);
    for factor in product.iter() {
        match factor {
            symbolica::atom::AtomView::Fun(function)
                if function.get_symbol() == crate::symbols::S.topo =>
            {
                topology = factor.to_owned()
            }
            _ => numerator = numerator * factor.to_owned(),
        }
    }
    let base_terms = match numerator.expand().as_view() {
        symbolica::atom::AtomView::Add(sum) => sum.get_nargs(),
        _ => 1,
    };
    let copies = target.div_ceil(base_terms.max(1)).max(1);
    let building = std::time::Instant::now();
    let spectators = Atom::add_many((1..=copies).map(|i| parse(&format!("vk4l_synth_x({i})"))));
    let block = (numerator * spectators).expand();
    let block_terms = match block.as_view() {
        symbolica::atom::AtomView::Add(sum) => sum.get_nargs(),
        _ => 1,
    };
    eprintln!(
        "SYNTH block {block_terms} monomials ({base_terms} x {copies}), {} MB, built in {:.1} s, peak {} MB",
        block.as_view().get_byte_size() >> 20,
        building.elapsed().as_secs_f64(),
        peak() >> 10
    );
    let loops: i64 = field("p_loops").parse().unwrap();
    let mut routing = (1..=loops)
        .map(|i| (parse(&format!("k({i})")), parse(&format!("k({i})"))))
        .collect::<Vec<_>>();
    routing.push((closure.momentum.clone(), closure.momentum.clone()));
    let raise = |line: usize| raise_line_power(&topology, line);
    let boundary = crate::p_integral::RemainderBoundary::new(
        "synthetic",
        closure.momentum.clone(),
        closure.mass_squared.clone(),
        routing,
        vec![block],
        (0..patterns)
            .map(|i| crate::p_integral::BoundaryPattern {
                topology: if i == 0 {
                    topology.clone()
                } else {
                    raise(i - 1)
                },
                closure_power: closure.power,
                terms: vec![(0, Atom::num(1))],
            })
            .collect(),
    )
    .unwrap();
    let started = std::time::Instant::now();
    let (classes, declined, (reads, hits)) = boundary
        .lower_counted(
            &engine,
            &settings,
            crate::p_integral::BoundaryOptions { strict: true },
        )
        .unwrap();
    eprintln!(
        "SYNTH lower {patterns} patterns: {:.1} s, {} key plans, {} declined, {reads} block reads, {hits} reuses, reuse={}, peak {} MB",
        started.elapsed().as_secs_f64(),
        classes.len(),
        declined.len(),
        std::env::var("VAKINT_P_BOUNDARY_BLOCKS").as_deref() != Ok("off"),
        peak() >> 10
    );
}

#[test]
fn massive_closure_showcase_sum() {
    let (Some(files), Some(out)) = (
        std::env::var_os("VAKINT_TEST_CLOSURE_FILES"),
        std::env::var_os("VAKINT_TEST_CLOSURE_OUT"),
    ) else {
        eprintln!("massive_closure_showcase_sum: environment unset, skipped");
        return;
    };
    let out = std::path::PathBuf::from(out);
    std::fs::create_dir_all(&out).unwrap();
    let mut engine = Vakint::new().unwrap();
    if let Some(package) = std::env::var_os("VAKINT_TEST_P5L_DRILL_PACKAGE") {
        let loading = std::time::Instant::now();
        engine
            .load_p_library_directory(std::path::Path::new(&package), false)
            .unwrap();
        engine.p_library().unwrap();
        eprintln!(
            "SHOWCASE_SUM package loaded in {:.1} s",
            loading.elapsed().as_secs_f64()
        );
    }
    let reading = std::time::Instant::now();
    let mut names = Vec::new();
    let mut inputs = Vec::new();
    let mut epsilon_terms = 0;
    for file in std::env::split_paths(&files) {
        let name = file
            .file_name()
            .unwrap()
            .to_string_lossy()
            .trim_end_matches(".zst")
            .trim_end_matches(".txt")
            .to_owned();
        let raw = std::process::Command::new("zstd")
            .args(["-qdc"])
            .arg(&file)
            .output()
            .unwrap();
        assert!(raw.status.success(), "{}: zstd failed", file.display());
        let text = String::from_utf8(raw.stdout).unwrap();
        let field = |key: &str| {
            text.lines()
                .find_map(|line| line.strip_prefix(key)?.strip_prefix(' '))
                .unwrap_or_else(|| panic!("{name}: no {key}"))
                .to_owned()
        };
        epsilon_terms = field("epsilon_terms").parse().unwrap();
        let closure = crate::p_integral::MassiveClosure {
            momentum: parse(&field("closure_momentum")),
            mass_squared: parse(&field("closure_mass_squared")),
            power: field("closure_power").parse().unwrap(),
        };
        let input = parse(&field("input"));
        let (topology, numerator) = match input.as_view() {
            symbolica::atom::AtomView::Mul(product) => {
                let mut topology = Atom::num(1);
                let mut numerator = Atom::num(1);
                for factor in product.iter() {
                    match factor {
                        symbolica::atom::AtomView::Fun(function)
                            if function.get_symbol() == crate::symbols::S.topo =>
                        {
                            topology = factor.to_owned()
                        }
                        _ => numerator = numerator * factor.to_owned(),
                    }
                }
                (topology, numerator)
            }
            other => panic!(
                "{name}: input is not numerator * topo: {:.100}",
                other.to_string()
            ),
        };
        let mut routing = (1..=5)
            .map(|i| (parse(&format!("k({i})")), parse(&format!("k({i})"))))
            .collect::<Vec<_>>();
        routing.push((closure.momentum.clone(), closure.momentum.clone()));
        names.push(name);
        inputs.push((
            closure,
            crate::p_integral::DirectClosureInput {
                topology,
                routing,
                numerator,
            },
        ));
    }
    let read_seconds = reading.elapsed().as_secs_f64();
    eprintln!(
        "SHOWCASE_SUM {} closures read in {read_seconds:.1} s",
        inputs.len()
    );
    let modes = std::env::var("VAKINT_TEST_CLOSURE_MODES").unwrap_or_else(|_| "retain".into());
    let make_settings = |mode: &str| {
        let master_mode = match mode {
            "retain" => RustRedMasterMode::Retain,
            "numerical" => RustRedMasterMode::Numerical,
            _ => RustRedMasterMode::Symbolic,
        };
        let arithmetic = std::env::var("VAKINT_TEST_CLOSURE_ARITHMETIC").ok();
        let mut settings = settings(master_mode, arithmetic.as_deref());
        // vk5l-physrun dry listing of uncovered boxed-node arrivals:
        // VAKINT_TEST_P_UNCOVERED=retain (PResidual symbols, Retain mode only),
        // VAKINT_TEST_P_OUTSIDE_RECIPES=<n> (0: no run-time recipe solve),
        // VAKINT_TEST_P_RULES_OUTSIDE=0 (no table rules outside the certified set).
        if std::env::var("VAKINT_TEST_P_UNCOVERED").as_deref() == Ok("retain") {
            settings.rustred_p_uncovered = crate::RustRedPUncoveredPolicy::Retain;
        }
        if let Some(recipes) = std::env::var("VAKINT_TEST_P_OUTSIDE_RECIPES")
            .ok()
            .and_then(|n| n.parse().ok())
        {
            settings.rustred_p_outside.max_recipes = recipes;
        }
        if std::env::var("VAKINT_TEST_P_RULES_OUTSIDE").as_deref() == Ok("0") {
            settings.rustred_p_outside.rules_outside = false;
        }
        settings.epsilon_symbol = "gammalooprs::ε".into();
        settings.mu_r_sq_symbol = "gammalooprs::μᵣ²".into();
        settings.integral_normalization_factor = crate::LoopNormalizationFactor::MSbar;
        settings.number_of_terms_in_epsilon_expansion = epsilon_terms;
        settings.run_time_decimal_precision = 100;
        settings.use_dot_product_notation = false;
        settings.n_cores = 1;
        settings
    };
    let mode_list = modes.split(',').collect::<Vec<_>>();
    let all_settings = mode_list
        .iter()
        .map(|mode| make_settings(mode))
        .collect::<Vec<_>>();
    // vk5l-physrun: several modes share ONE lowering and ONE reduction
    // (evaluate_direct_sum_modes) unless VAKINT_TEST_SUM_SHARED=0.
    let shared =
        mode_list.len() > 1 && std::env::var("VAKINT_TEST_SUM_SHARED").as_deref() != Ok("0");
    let mut values = Vec::new();
    if shared {
        let started = std::time::Instant::now();
        let refs = all_settings.iter().collect::<Vec<_>>();
        let results =
            crate::p_integral::MassiveClosure::evaluate_direct_sum_modes(&engine, &refs, &inputs)
                .unwrap();
        let seconds = started.elapsed().as_secs_f64();
        eprintln!(
            "SHOWCASE_SUM shared lowering + reduction for {} modes in {seconds:.1} s",
            mode_list.len()
        );
        values.extend(results.into_iter().map(|value| (value, seconds)));
    } else {
        for settings in &all_settings {
            let started = std::time::Instant::now();
            let value =
                crate::p_integral::MassiveClosure::evaluate_direct_sum(&engine, settings, &inputs)
                    .unwrap();
            values.push((value, started.elapsed().as_secs_f64()));
        }
    }
    for ((mode, settings), (value, seconds)) in
        mode_list.iter().copied().zip(&all_settings).zip(values)
    {
        let text = value.to_canonical_string();
        std::fs::write(out.join(format!("sum.{mode}.sym")), &text).unwrap();
        let mut check = serde_json::Value::Null;
        if let Some(dir) = std::env::var_os("VAKINT_TEST_CLOSURE_SUM_CHECK") {
            let dir = std::path::PathBuf::from(dir);
            let mut total = Atom::Zero;
            let mut found = 0;
            for name in &names {
                let path = dir.join(name).join(format!("{name}.{mode}.sym"));
                if let Ok(text) = std::fs::read_to_string(&path) {
                    total += parse(&text);
                    found += 1;
                }
            }
            let total = total.expand();
            let identical = total == value;
            let equal = identical || (&total - &value).expand().is_zero();
            eprintln!(
                "SHOWCASE_SUM check {mode}: {found}/{} per-closure results, identical={identical} equal={equal}",
                names.len()
            );
            check = serde_json::json!({"found": found, "identical": identical, "equal": equal});
        }
        let dispatch = engine.p_dispatch_counters();
        let load = std::fs::read_to_string("/proc/loadavg").unwrap_or_default();
        let summary = serde_json::json!({
            "dispatch_process_total": {"table": dispatch.table, "explicit": dispatch.explicit, "cached": dispatch.cached,
                                       "recipe_solved": dispatch.recipe_solved, "refused": dispatch.refused},
            "loadavg_at_end": load.trim(), "utc_end_unix": std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0),
            "closures": names, "mode": mode, "seconds": seconds, "read_seconds": read_seconds,
            "result_bytes": text.len(), "result_terms": drill_term_count(&value.expand()),
            "arithmetic": settings.ordered_weight_arithmetic.to_string(), "check": check,
        });
        std::fs::write(out.join(format!("sum.{mode}.json")), summary.to_string()).unwrap();
        eprintln!(
            "SHOWCASE_SUM {mode}: {} closures in {seconds:.1} s, {} terms",
            names.len(),
            drill_term_count(&value.expand())
        );
    }
}

/// vk5l WP4 end-to-end showcase: ready-made massive-closure requests
/// (`pint6.massive-closure-request.v1`, zstd, one per file, the beta team's
/// e2e-showcase input) through `MassiveClosure::evaluate_batch` as the R*
/// pipeline calls it. Environment: `VAKINT_TEST_CLOSURE_FILES` (`:`-separated
/// paths), `VAKINT_TEST_CLOSURE_OUT` (output directory),
/// `VAKINT_TEST_P5L_DRILL_PACKAGE` (the p package), `VAKINT_TEST_CLOSURE_MODES`
/// (`retain,symbolic` default). Per closure and mode: `<name>.<mode>.sym`
/// (Symbolica canonical string) and `<name>.<mode>.json` (seconds, peak RSS,
/// sizes, reduction and relation-cache counters).
#[test]
fn massive_closure_showcase() {
    let (Some(files), Some(out)) = (
        std::env::var_os("VAKINT_TEST_CLOSURE_FILES"),
        std::env::var_os("VAKINT_TEST_CLOSURE_OUT"),
    ) else {
        eprintln!("massive_closure_showcase: environment unset, skipped");
        return;
    };
    let out = std::path::PathBuf::from(out);
    std::fs::create_dir_all(&out).unwrap();
    let mut engine = Vakint::new().unwrap();
    if let Some(package) = std::env::var_os("VAKINT_TEST_P5L_DRILL_PACKAGE") {
        let loading = std::time::Instant::now();
        engine
            .load_p_library_directory(std::path::Path::new(&package), false)
            .unwrap();
        engine.p_library().unwrap();
        eprintln!(
            "SHOWCASE package loaded in {:.1} s",
            loading.elapsed().as_secs_f64()
        );
    }
    let peak = || -> u64 {
        std::fs::read_to_string("/proc/self/status")
            .ok()
            .and_then(|status| {
                status
                    .lines()
                    .find(|line| line.starts_with("VmHWM:"))
                    .and_then(|line| line.split_whitespace().nth(1)?.parse::<u64>().ok())
            })
            .unwrap_or(0)
            * 1024
    };
    let modes =
        std::env::var("VAKINT_TEST_CLOSURE_MODES").unwrap_or_else(|_| "retain,symbolic".into());
    for file in std::env::split_paths(&files) {
        let name = file
            .file_name()
            .unwrap()
            .to_string_lossy()
            .trim_end_matches(".zst")
            .trim_end_matches(".txt")
            .to_owned();
        let reading = std::time::Instant::now();
        let rss = || -> u64 {
            std::fs::read_to_string("/proc/self/status")
                .ok()
                .and_then(|status| {
                    status
                        .lines()
                        .find(|line| line.starts_with("VmRSS:"))
                        .and_then(|line| line.split_whitespace().nth(1)?.parse::<u64>().ok())
                })
                .unwrap_or(0)
                * 1024
        };
        let raw = std::process::Command::new("zstd")
            .args(["-qdc"])
            .arg(&file)
            .output()
            .unwrap();
        assert!(raw.status.success(), "{}: zstd failed", file.display());
        let text = String::from_utf8(raw.stdout).unwrap();
        let decompress_seconds = reading.elapsed().as_secs_f64();
        let decompress_rss = rss();
        let field = |key: &str| {
            text.lines()
                .find_map(|line| line.strip_prefix(key)?.strip_prefix(' '))
                .unwrap_or_else(|| panic!("{name}: no {key}"))
                .to_owned()
        };
        let epsilon_terms: i64 = field("epsilon_terms").parse().unwrap();
        let closure = crate::p_integral::MassiveClosure {
            momentum: parse(&field("closure_momentum")),
            mass_squared: parse(&field("closure_mass_squared")),
            power: field("closure_power").parse().unwrap(),
        };
        let input_text = field("input");
        let parsing = std::time::Instant::now();
        let input = parse(&input_text);
        let parse_seconds = parsing.elapsed().as_secs_f64();
        let parse_rss = rss();
        drop(text);
        let read_seconds = reading.elapsed().as_secs_f64();
        eprintln!(
            "SHOWCASE {name} ingest: decompress {decompress_seconds:.1} s, parse {parse_seconds:.1} s, \
             input {} terms, rss {:.2} GB",
            drill_term_count(&input),
            parse_rss as f64 / 1e9
        );
        for mode in modes.split(',') {
            let master_mode = match mode {
                "retain" => RustRedMasterMode::Retain,
                "numerical" => RustRedMasterMode::Numerical,
                _ => RustRedMasterMode::Symbolic,
            };
            // vk5l-physrun: VAKINT_TEST_CLOSURE_ARITHMETIC = RustRed ordered weight
            // arithmetic text (exact | modular_sampled[:..] | modular_polynomial[:..]).
            let arithmetic = std::env::var("VAKINT_TEST_CLOSURE_ARITHMETIC").ok();
            let mut settings = settings(master_mode, arithmetic.as_deref());
            // vk5l-physrun dry listing of uncovered boxed-node arrivals:
            // VAKINT_TEST_P_UNCOVERED=retain (PResidual symbols, Retain mode only),
            // VAKINT_TEST_P_OUTSIDE_RECIPES=<n> (0: no run-time recipe solve),
            // VAKINT_TEST_P_RULES_OUTSIDE=0 (no table rules outside the certified set).
            if std::env::var("VAKINT_TEST_P_UNCOVERED").as_deref() == Ok("retain") {
                settings.rustred_p_uncovered = crate::RustRedPUncoveredPolicy::Retain;
            }
            if let Some(recipes) = std::env::var("VAKINT_TEST_P_OUTSIDE_RECIPES")
                .ok()
                .and_then(|n| n.parse().ok())
            {
                settings.rustred_p_outside.max_recipes = recipes;
            }
            if std::env::var("VAKINT_TEST_P_RULES_OUTSIDE").as_deref() == Ok("0") {
                settings.rustred_p_outside.rules_outside = false;
            }
            settings.epsilon_symbol = "gammalooprs::ε".into();
            settings.mu_r_sq_symbol = "gammalooprs::μᵣ²".into();
            settings.integral_normalization_factor = crate::LoopNormalizationFactor::MSbar;
            settings.number_of_terms_in_epsilon_expansion = epsilon_terms;
            settings.run_time_decimal_precision = 100;
            settings.use_dot_product_notation = false;
            settings.n_cores = std::env::var("VAKINT_TEST_P5L_DRILL_CORES")
                .ok()
                .and_then(|cores| cores.parse().ok())
                .unwrap_or(1);
            let before = engine.p_relation_cache_statistics();
            let started = std::time::Instant::now();
            let result = if std::env::var("VAKINT_TEST_CLOSURE_DIRECT").as_deref() == Ok("1") {
                // The input is `numerator · topo(...)` in the cut momenta: the
                // direct lowering's input with the identity routing.
                let (topology, numerator) = match input.as_view() {
                    symbolica::atom::AtomView::Mul(product) => {
                        let mut topology = Atom::num(1);
                        let mut numerator = Atom::num(1);
                        for factor in product.iter() {
                            match factor {
                                symbolica::atom::AtomView::Fun(function)
                                    if function.get_symbol() == crate::symbols::S.topo =>
                                {
                                    topology = factor.to_owned()
                                }
                                _ => numerator = numerator * factor.to_owned(),
                            }
                        }
                        (topology, numerator)
                    }
                    other => panic!(
                        "{name}: input is not numerator * topo: {:.100}",
                        other.to_string()
                    ),
                };
                // vk5l-physrun: the records are already routed into the cut
                // momenta, external included: identity routing of k(1..5) AND of
                // the closing head (else "an unknown momentum in a scalar product").
                let mut routing = (1..=5)
                    .map(|i| (parse(&format!("k({i})")), parse(&format!("k({i})"))))
                    .collect::<Vec<_>>();
                routing.push((closure.momentum.clone(), closure.momentum.clone()));
                crate::p_integral::MassiveClosure::evaluate_direct_batch(
                    &engine,
                    &settings,
                    &[(
                        closure.clone(),
                        crate::p_integral::DirectClosureInput {
                            topology,
                            routing,
                            numerator,
                        },
                    )],
                )
            } else {
                crate::p_integral::MassiveClosure::evaluate_batch(
                    &engine,
                    &settings,
                    &[(closure.clone(), input.clone())],
                )
            };
            let seconds = started.elapsed().as_secs_f64();
            let stats = super::last_stats();
            let stages = crate::p_integral::massive::take_stage_seconds();
            let after = engine.p_relation_cache_statistics();
            let delta = after.cache.since(&before.cache);
            let (status, value_text, terms) = match &result {
                Ok(values) => {
                    let value = &values[0];
                    let text = value.to_canonical_string();
                    let terms = drill_term_count(&value.expand());
                    std::fs::write(out.join(format!("{name}.{mode}.sym")), &text).unwrap();
                    ("ok".to_owned(), text.len(), terms)
                }
                Err(error) => (format!("refused: {error}"), 0, 0),
            };
            let summary = serde_json::json!({
                "closure": name, "mode": mode, "status": status,
                "seconds": seconds, "read_seconds": read_seconds, "peak_rss_bytes": peak(),
                "ingest": {"decompress_s": decompress_seconds, "parse_s": parse_seconds,
                           "rss_after_decompress": decompress_rss, "rss_after_parse": parse_rss},
                "stages": stages.iter().map(|(stage, seconds)| (stage.to_string(), serde_json::json!(seconds)))
                    .collect::<serde_json::Map<_, _>>(),
                "rss_after": rss(),
                "input_bytes": input_text.len(), "result_bytes": value_text, "result_terms": terms,
                "epsilon_terms": epsilon_terms,
                "arithmetic": settings.ordered_weight_arithmetic.to_string(),
                "reduce": {"channels": stats.channels, "reductions": stats.reductions, "pairs": stats.pairs,
                           "split_s": stats.split_seconds, "lower_s": stats.lower_seconds,
                           "prepare_s": stats.prepare_seconds, "library_s": stats.library_seconds,
                           "render_s": stats.render_seconds, "memo_hits": stats.memo_hits},
                "relation_cache": {"shipped_hits": delta.shipped_hits, "cache_hits": delta.cache_hits,
                                   "misses": delta.misses, "inserts": delta.inserts},
                "recipes": {"solves": after.recipes.solves - before.recipes.solves,
                            "outside_solves": after.recipes.outside_solves - before.recipes.outside_solves,
                            "failures": after.recipes.failures - before.recipes.failures},
            });
            std::fs::write(out.join(format!("{name}.{mode}.json")), summary.to_string()).unwrap();
            eprintln!(
                "SHOWCASE {name} {mode}: {status:.200} in {seconds:.1} s, {terms} terms, peak {:.2} GB",
                peak() as f64 / 1e9
            );
        }
    }
}
