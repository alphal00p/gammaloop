//! The RustRed p backend against Forcer's frozen one-loop battery, with the
//! injected reference library `rustred::propagator::testing::
//! one_loop_reference_library()` (the massless bubble behind the real engine):
//! exact agreement in the tag basis, agreement after the Symbolic and
//! Numerical finalize, the Wick sign applied once, channel splitting, the
//! settings' arithmetic, the step cache, the direct entry points, the timing
//! span and the dispatch without p data.

use std::collections::BTreeMap;

use rustred::propagator::library::{PBasis, PReduction, PTerm};
use rustred::propagator::tag::GammaTag;
use symbolica::atom::{Atom, AtomCore};

use super::*;
use crate::p_integral::compare::tag_basis_difference;
use crate::p_integral::forcer::fixtures::{ForcerCase, ForcerFixture, registry_family};
use crate::utils::vakint_macros::vk_parse;
use crate::{
    EvaluationMethod, EvaluationOrder, ForcerOptions, LoopNormalizationFactor,
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
            .collect()
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
/// 2026-11-07) and the exact arithmetic gives the same result.
#[test]
fn the_settings_arithmetic_is_used_and_exact_agrees() {
    let modular = settings(RustRedMasterMode::Retain, None);
    let exact = settings(RustRedMasterMode::Retain, Some("exact"));
    assert_eq!(
        request_options(&modular).arithmetic,
        modular.rustred_ordered_weight_arithmetic()
    );
    assert_eq!(
        request_options(&modular).arithmetic.name(),
        "modular_polynomial"
    );
    assert_eq!(request_options(&exact).arithmetic.name(), "exact");
    assert_eq!(
        request_options(&modular).scheduler,
        crate::rustred_evaluation::ordered_scheduler(&modular)
    );
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
/// library identity: a read-write directory is written under that name and a
/// second engine replays it read-only, recomputing nothing, with the same
/// result.
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
    assert!(directory.path().join(&name).is_dir(), "{name}");
    let statistics = first.ordered_cache_statistics();
    let program = statistics
        .programs
        .iter()
        .find(|program| program.program == "p-library-1-1")
        .unwrap();
    assert_eq!(program.rule_tables, name);
    assert_eq!(program.data_version, "injected");
    assert!(program.native.inserts > 0);

    let mut reading = writing.clone();
    reading.ordered_cache = OrderedCacheSettings::read_only(directory.path());
    let second = engine();
    assert_eq!(second.evaluate(&reading, input.as_view()).unwrap(), cold);
    // Nothing is recomputed: the steps (and the tail table RustRed keeps
    // next to them) replay from the directory, while the cold run above
    // computed and recorded its steps.
    let replayed = second.ordered_cache_statistics();
    assert_eq!((replayed.misses, replayed.inserts), (0, 0), "{replayed:?}");
    assert!(program.native.misses > 0);
    let files = std::fs::read_dir(directory.path().join(&name))
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        .collect::<std::collections::BTreeSet<_>>();
    assert!(files.contains("ordered-step-cache.txt"), "{files:?}");
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

/// Without p data, p-integrals keep their Forcer dispatch (default order with
/// a FORM 5), the selection does not depend on the p backend's structural
/// support; with an injected library covering the loop count RustRed is
/// chosen, and a p-integral outside its scope stays with Forcer.
#[test]
fn without_p_data_p_integrals_go_to_forcer() {
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
    let bare = Vakint::new().unwrap();
    for input in [bubble, kite] {
        let topology = topology(&bare, input);
        assert!(crate::rustred_evaluation::supports(
            &settings,
            &topology,
            &RustRedEvaluationOptions::default()
        ));
        assert!(!crate::rustred_evaluation::has_data(
            &bare, &settings, &topology
        ));
        assert!(matches!(
            bare.evaluation_method(&settings, &topology),
            Some(EvaluationMethod::Forcer(_))
        ));
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
}
