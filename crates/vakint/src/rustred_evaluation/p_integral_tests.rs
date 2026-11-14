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
    let mut settings = settings(RustRedMasterMode::Retain, None);
    settings.number_of_terms_in_epsilon_expansion = 6;
    settings.n_cores = std::env::var("VAKINT_TEST_P5L_DRILL_CORES")
        .ok()
        .and_then(|cores| cores.parse().ok())
        .unwrap_or(1);
    if std::env::var("VAKINT_TEST_P5L_DRILL_RETAIN").as_deref() == Ok("1") {
        settings.rustred_p_uncovered = crate::RustRedPUncoveredPolicy::Retain;
    }
    let cases = std::env::var("VAKINT_TEST_P5L_DRILL").unwrap_or_default();
    for case in cases.split(';').filter(|case| !case.trim().is_empty()) {
        let mut parts = case.splitn(3, '=');
        let (label, expected, input) = (
            parts.next().unwrap().trim(),
            parts.next().unwrap().trim(),
            parts.next().unwrap().trim(),
        );
        let input = parse(input);
        let started = std::time::Instant::now();
        let result = engine.evaluate(&settings, input.as_view());
        let seconds = started.elapsed().as_secs_f64();
        match (expected, result) {
            // An exact identity (e.g. an IBP relation with integer
            // coefficients): every summand reduced on its own must be
            // nonzero, and both the sum of those reductions and the
            // reduction of the whole input must vanish identically.
            ("zero", Ok(value)) => {
                let summands: Vec<Atom> = match input.as_view() {
                    symbolica::atom::AtomView::Add(add) => add.iter().map(|t| t.to_owned()).collect(),
                    _ => vec![input.clone()],
                };
                let mut total = Atom::Zero;
                let mut nonzero = 0usize;
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
                assert!(nonzero >= 2, "{label}: fewer than two nonzero summands, no check");
                assert!(whole.is_zero() && total.is_zero(), "{label}: identity violated");
                eprintln!(
                    "DRILL {label}: IDENTITY HOLDS EXACTLY ({} summands, {nonzero} nonzero)",
                    summands.len()
                );
            }
            ("ok", Ok(value)) => {
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
