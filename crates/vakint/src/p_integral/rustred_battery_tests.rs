//! Tier A of the Forcer battery (M2-V5, M3-V9; plans M2 and M3, gate G7):
//! the RustRed p path against Forcer's frozen outputs, without FORM.
//!
//! Every case of the frozen batteries `tests/fixtures/forcer/p{1,2,3}l_battery.json`
//! (45 + 23 + 19 = 87) and of the anchors `p_anchors.json` (11) is checked
//! three ways, every case of the gate part of the grids `p{1,2,3}l_grid.json`
//! ([`in_gate`]) by check 1, the one- and two-loop grids also by checks 2 and
//! 3 (M3 I-6; Tier C runs every grid case through all three modes of the
//! public API):
//!
//! 1. exact: RustRed's reduced values (the request level below the finalizer,
//!    `rustred_evaluation::p_integral::reduce_requests`) equal Forcer's
//!    (`forcer_input` parts and `combine_parts` on the frozen outputs, i.e.
//!    `ForcerCase::frozen`), scale power by scale power, in the
//!    `Π Γ(1+mε)^e × Π PMaster` tag basis (`compare::tag_basis_difference`);
//! 2. symbolic: `Vakint::evaluate` with RustRed in `Symbolic` mode through
//!    `ε^(L+2)` equals the frozen value finalized like a Forcer evaluation;
//! 3. numerical: the same in `Numerical` mode at 40 digits, within `1e-30`
//!    at `p(i) = (0,0,0,1)` (`Q² = 1`) and `μ² = 2`, where the logarithms of
//!    `μ²/Q²` count.
//!
//! The Wick sign `(-1)^(Σa+r)` of the Forcer path is pinned to be visible on
//! the batteries ([`rustred_carries_the_wick_sign_of_odd_cases_once`]).
//!
//! All of it runs under the default Ordered arithmetic (`modular_polynomial`,
//! the `Z_p(d)` replay vakint users get; coordinator note 2026-11-07 09:40Z).
//! The other two arithmetics, `exact` and `modular_sampled` (M2 exit
//! criterion 2), repeat check 1 on the frozen batteries, the anchors and the
//! gate part of the grids ([`arithmetic_subset`]) and checks 2 and 3 on
//! [`FINALIZE_CASES`].
//!
//! The engine of the battery is [`battery_engine`]: `Vakint::new()` with its
//! built-in 1L-3L p library (M2-V4), or with the scope 1-4 library of M3
//! injected when `VAKINT_TEST_P_LIBRARY_1_4` names one ([`attached_library`]);
//! the battery covers the library's scope ([`battery_loops`]). Four loops
//! (M3-V9): on a 1-4 library, [`FAST_FOUR_LOOP_CASES`] run checks 1-3 through
//! the case's frozen depth (the deepest Laurent depth Forcer's value
//! finalizes at, `p4l_battery.json`), the 11-term value against Forcer's
//! frozen one (M3 I-7) or the expected `MissingMasterOrders` (`Oep`); the
//! complete four-loop battery, the anchors `p4l_anchors.json` and the gate
//! part of `p4l_grid.json` run in the full profile through the public API
//! (M3 D-M3-10). Without a four-loop library a four-loop case fails closed
//! on the RustRed path ([`rustred_covers_the_library_scope_and_no_more`]).
//! The full grids and the benchmark cells through vakint's public API are
//! Tier C (`tests/p_integral_rustred_full.rs`); live FORM comparisons are
//! Tier B (`tests/p_integral_rustred_tests.rs`).

use std::collections::{BTreeMap, HashMap};
use std::ops::RangeInclusive;
use std::sync::{Arc, OnceLock};

use rustred::propagator::library::PLibrary;
use rustred_app::propagator::{PLibraryOptions, PLibrarySource, load_p_library};
use symbolica::atom::{Atom, AtomCore, Symbol};
use symbolica::domains::float::Float;

use super::compare::tag_basis_difference;
use super::forcer::fixtures::{Finalized, ForcerCase, ForcerFixture, I7_TERMS, registry_family};
use super::forcer::{ForcerPart, combine_parts};
use super::{PropagatorRequest, ReducedPropagatorValue};
use crate::rustred_evaluation::p_integral::{has_p_data, reduce_requests};
use crate::rustred_evaluation::{RustRedEvaluationError, RustRedMasterMode, combine_exact};
use crate::utils::simplify_real;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, Momentum, RustRedEvaluationOptions,
    Vakint, VakintError, VakintSettings,
};

/// The environment variable naming M3's scope 1-4 p library package (set in
/// the M3 environment from the PP-M3b pin; `tests/support/p_library4.rs`).
const P_LIBRARY_1_4_VARIABLE: &str = "VAKINT_TEST_P_LIBRARY_1_4";

/// The p library attached to the battery: the package that
/// [`P_LIBRARY_1_4_VARIABLE`] names, loaded once per process like vakint
/// loads a library (without the descent certificate), or none. A named
/// package that does not load is a failure, never a fallback.
fn attached_library() -> Option<Arc<PLibrary>> {
    static LIBRARY: OnceLock<Option<Arc<PLibrary>>> = OnceLock::new();
    LIBRARY
        .get_or_init(|| {
            let directory = std::env::var_os(P_LIBRARY_1_4_VARIABLE)?;
            let library = load_p_library(
                PLibrarySource::Directory(directory.clone().into()),
                &PLibraryOptions {
                    verify_descent: false,
                },
            );
            Some(library.unwrap_or_else(|error| {
                panic!("{P_LIBRARY_1_4_VARIABLE}={}: {error}", directory.display())
            }))
        })
        .clone()
}

/// The engine of the battery: `Vakint::new()` with the [`attached_library`]
/// injected, else with its built-in p library (M2-V4).
fn battery_engine() -> Vakint {
    let mut engine = Vakint::new().unwrap();
    if let Some(library) = attached_library() {
        engine.inject_p_library(library).unwrap();
    }
    engine
}

/// The loop counts `engine` reduces with RustRed: the scope of its p library
/// (`1..=3` for the built-in library, `1..=4` with a scope 1-4 library).
fn battery_loops(engine: &Vakint) -> RangeInclusive<usize> {
    1..=(1..=4)
        .take_while(|&loops| has_p_data(engine, loops))
        .last()
        .unwrap_or(0)
}

/// A frozen fixture file and the number of cases it must hold.
struct Fixture {
    file: &'static str,
    json: &'static str,
    cases: usize,
}

const BATTERIES: [Fixture; 3] = [
    Fixture {
        file: "p1l_battery.json",
        json: include_str!("../../tests/fixtures/forcer/p1l_battery.json"),
        cases: 45,
    },
    Fixture {
        file: "p2l_battery.json",
        json: include_str!("../../tests/fixtures/forcer/p2l_battery.json"),
        cases: 23,
    },
    Fixture {
        file: "p3l_battery.json",
        json: include_str!("../../tests/fixtures/forcer/p3l_battery.json"),
        cases: 19,
    },
];

const GRIDS: [Fixture; 3] = [
    Fixture {
        file: "p1l_grid.json",
        json: include_str!("../../tests/fixtures/forcer/p1l_grid.json"),
        cases: 17,
    },
    Fixture {
        file: "p2l_grid.json",
        json: include_str!("../../tests/fixtures/forcer/p2l_grid.json"),
        cases: 39,
    },
    Fixture {
        file: "p3l_grid.json",
        json: include_str!("../../tests/fixtures/forcer/p3l_grid.json"),
        cases: 417,
    },
];

const ANCHORS: Fixture = Fixture {
    file: "p_anchors.json",
    json: include_str!("../../tests/fixtures/forcer/p_anchors.json"),
    cases: 11,
};

const FOUR_LOOP_BATTERY: Fixture = Fixture {
    file: "p4l_battery.json",
    json: include_str!("../../tests/fixtures/forcer/p4l_battery.json"),
    cases: 54,
};

impl Fixture {
    fn parse(&self) -> ForcerFixture {
        let fixture = ForcerFixture::parse(self.json)
            .unwrap_or_else(|error| panic!("{}: {error}", self.file));
        assert_eq!(fixture.cases.len(), self.cases, "{}", self.file);
        fixture
    }
}

/// The grid cases of the fast gate: all but the 350 labelled three-loop
/// sectors at rank 0 and 2 (`sector_P3L_*`), which Tier C checks through the
/// public API. The benchmark cells (`bench_*`), the one- and two-loop sectors
/// and every dotted sector are in.
fn in_gate(case: &ForcerCase) -> bool {
    !case.name.starts_with("sector_P3L_")
}

/// The cases of `loops` loops the battery checks: the frozen battery, the
/// gate part of the grid and the anchors of that loop count, as
/// `(file, case)`.
fn cases(loops: usize) -> Vec<(&'static str, ForcerCase)> {
    let mut cases = Vec::new();
    for (fixture, gate) in [
        (&BATTERIES[loops - 1], None),
        (&GRIDS[loops - 1], Some(in_gate as fn(&ForcerCase) -> bool)),
        (&ANCHORS, None),
    ] {
        for case in fixture.parse().cases {
            if case.loops().unwrap() == loops && gate.is_none_or(|gate| gate(&case)) {
                cases.push((fixture.file, case));
            }
        }
    }
    cases
}

fn parse(text: &str) -> Atom {
    Vakint::initialize_vakint_symbols();
    vk_parse!(text).unwrap()
}

/// The Laurent depth of the battery for `loops` loops: through `ε^(L+2)`,
/// `2L + 3` terms from `ε^-L`.
fn battery_terms(loops: usize) -> i64 {
    2 * loops as i64 + 3
}

/// The Laurent depth a case is checked through: its frozen depth (the
/// deepest at which Forcer's value finalizes, M3 fixtures), else
/// [`battery_terms`].
fn case_terms(case: &ForcerCase) -> i64 {
    case.depth
        .as_ref()
        .map_or_else(|| battery_terms(case.loops().unwrap()), |depth| depth.terms)
}

/// RustRed alone in `mode`, G-scheme, through `terms` Laurent terms (a FORM
/// call fails loudly), 40 digits, one core; the settings' default Ordered
/// arithmetic (`modular_polynomial`) or `arithmetic`.
fn settings(mode: RustRedMasterMode, terms: i64, arithmetic: Option<&str>) -> VakintSettings {
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions { master_mode: mode },
        )]),
        form_exe_path: "/rustred-battery-must-not-call-form".into(),
        integral_normalization_factor: LoopNormalizationFactor::GScheme,
        number_of_terms_in_epsilon_expansion: terms,
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

/// `reduced` per scale power, summed, zero sums dropped.
fn by_scale_power(values: impl IntoIterator<Item = ReducedPropagatorValue>) -> BTreeMap<i64, Atom> {
    let mut sums = BTreeMap::<i64, Atom>::new();
    for value in values {
        *sums.entry(value.scale_power).or_default() += value.reduced;
    }
    sums.retain(|_, reduced| !reduced.expand().is_zero());
    sums
}

/// The RustRed request of a fixture case: the input Forcer reduced.
fn rustred_reduced(
    engine: &Vakint,
    settings: &VakintSettings,
    case: &ForcerCase,
) -> Result<BTreeMap<i64, Atom>, VakintError> {
    let numerator = parse(&case.numerator);
    let request = PropagatorRequest {
        family: registry_family(&case.family)?,
        powers: &case.powers,
        canonical_numerator: numerator.as_view(),
        external: parse(&case.external),
    };
    let values = reduce_requests(engine, settings, &[request], None)?;
    assert!(values.len() <= 1, "{}: one external head", case.name);
    Ok(by_scale_power(values.into_values().flatten()))
}

/// Whether two reduced values agree exactly: the same scale powers, and a
/// zero difference in the tag basis at each; `Err` describes the first
/// difference.
fn reduced_agree(
    actual: &BTreeMap<i64, Atom>,
    expected: &BTreeMap<i64, Atom>,
    epsilon: Symbol,
) -> Result<(), String> {
    if actual.keys().ne(expected.keys()) {
        return Err(format!(
            "scale powers {:?} vs {:?}",
            actual.keys().collect::<Vec<_>>(),
            expected.keys().collect::<Vec<_>>()
        ));
    }
    for (w, value) in actual {
        let difference = tag_basis_difference(value, &expected[w], epsilon)
            .map_err(|error| format!("w = {w}: {error}"))?;
        if !difference.is_empty() {
            return Err(format!(
                "w = {w}: {value} vs {}: {difference:?}",
                expected[w]
            ));
        }
    }
    Ok(())
}

/// Check 1: RustRed's reduced values are Forcer's frozen ones, exactly in the
/// tag basis, scale power by scale power; a case with a value is nonzero.
fn assert_reduced_exactly(
    engine: &Vakint,
    settings: &VakintSettings,
    file: &str,
    case: &ForcerCase,
) {
    let name = format!("{file}: {}", case.name);
    let epsilon = epsilon(settings);
    let rustred =
        rustred_reduced(engine, settings, case).unwrap_or_else(|error| panic!("{name}: {error}"));
    let forcer = by_scale_power(case.frozen(epsilon).unwrap());
    assert_eq!(rustred.is_empty(), case.outputs.is_empty(), "{name}");
    reduced_agree(&rustred, &forcer, epsilon)
        .unwrap_or_else(|difference| panic!("{name}: RustRed vs Forcer: {difference}"));
}

/// Check 2: the case through `Vakint::evaluate` (matching, the p backend and
/// the shared finalizer) in `Symbolic` mode through `terms` Laurent terms is
/// its frozen value finalized like a Forcer evaluation (which returns
/// through `simplify_real`).
fn assert_symbolic(
    engine: &Vakint,
    terms: i64,
    arithmetic: Option<&str>,
    file: &str,
    case: &ForcerCase,
) {
    let settings = settings(RustRedMasterMode::Symbolic, terms, arithmetic);
    let input = case.expression().unwrap();
    let actual = engine
        .evaluate(&settings, input.as_view())
        .unwrap_or_else(|error| panic!("{file}: {}: {error}", case.name));
    let expected = simplify_real(
        case.finalized(&settings, RustRedMasterMode::Symbolic)
            .unwrap()
            .as_view(),
    );
    let difference = (&actual - &expected).expand();
    assert!(
        difference.is_zero(),
        "{file}: {}: RustRed {actual} vs Forcer {expected}: {difference}",
        case.name
    );
}

/// A point of the numerical checks: `p(1..4) = (0,0,0,1)`, so
/// `Q² = -p² = 1`, `μ²` and `mass = 3` (the free coefficient of the two-loop
/// battery case `kite_coefficients`), at the precision of `settings`.
type NumericalPoint = (
    HashMap<String, Float, ahash::RandomState>,
    HashMap<usize, Momentum, ahash::RandomState>,
);

fn numerical_point(settings: &VakintSettings, mu_squared: i64) -> NumericalPoint {
    let bits = settings.get_binary_precision();
    let float = |value: i64| Float::with_val(bits, value);
    let params = HashMap::from_iter([
        (settings.mu_r_sq_symbol.clone(), float(mu_squared)),
        ("mass".to_owned(), float(3)),
    ]);
    let externals = (1..=4)
        .map(|head| {
            (
                head,
                Momentum::Real((float(0), float(0), float(0), float(1))),
            )
        })
        .collect();
    (params, externals)
}

/// `actual` and `expected` evaluated to the digits of `settings` at
/// [`numerical_point`] with `μ² = mu_squared` agree within `1e-30`.
fn assert_numerically_equal(
    engine: &Vakint,
    settings: &VakintSettings,
    mu_squared: i64,
    name: &str,
    actual: &Atom,
    expected: &Atom,
) {
    let (params, externals) = numerical_point(settings, mu_squared);
    let at_point = |value: &Atom| {
        engine
            .numerical_evaluation(
                settings,
                value.as_view(),
                &params,
                &Default::default(),
                Some(&externals),
            )
            .unwrap_or_else(|error| panic!("{name}: {value}: {error}"))
            .0
    };
    let (actual, expected) = (at_point(actual), at_point(expected));
    let (equal, diagnostic) = actual.does_approx_match(&expected, None, 1e-30, 5.0);
    assert!(equal, "{name}: {diagnostic}");
}

/// Check 3: the same in `Numerical` mode through `terms` Laurent terms, both
/// sides evaluated to 40 digits at [`numerical_point`] with `μ² = 2` (the
/// logarithms of `μ²/Q²` count) and equal within `1e-30`.
fn assert_numerical(
    engine: &Vakint,
    terms: i64,
    arithmetic: Option<&str>,
    file: &str,
    case: &ForcerCase,
) {
    let name = format!("{file}: {}", case.name);
    let settings = settings(RustRedMasterMode::Numerical, terms, arithmetic);
    let input = case.expression().unwrap();
    let actual = engine
        .evaluate(&settings, input.as_view())
        .unwrap_or_else(|error| panic!("{name}: {error}"));
    let expected = engine
        .resolve_rustred_numerical(
            combine_exact(
                case.finalized(&settings, RustRedMasterMode::Numerical)
                    .unwrap(),
            ),
            &settings,
        )
        .unwrap();
    assert_numerically_equal(engine, &settings, 2, &name, &actual, &expected);
}

/// Check 1 on every case of `loops` loops ([`cases`]), whose number is
/// pinned (`battery + grid gate part + anchors`), and checks 2 and 3 on the
/// battery and anchor cases and, at one and two loops, on the grid cases too
/// (M3 I-6; the three-loop grid cases go through every mode in Tier C).
fn assert_loop_count_matches_forcer(loops: usize, expected_cases: usize) {
    let engine = battery_engine();
    assert!(
        battery_loops(&engine).contains(&loops),
        "no RustRed p library for {loops} loops"
    );
    let cases = cases(loops);
    assert_eq!(cases.len(), expected_cases, "{loops}-loop cases");
    let terms = battery_terms(loops);
    let retain = settings(RustRedMasterMode::Retain, terms, None);
    for (file, case) in &cases {
        assert_reduced_exactly(&engine, &retain, file, case);
        if *file != GRIDS[loops - 1].file || loops <= 2 {
            assert_symbolic(&engine, terms, None, file, case);
            assert_numerical(&engine, terms, None, file, case);
        }
    }
}

/// 1L: 45 battery cases, the 17 grid cases (the 15 benchmark cells and the
/// two labelled sectors) and the 3 one-loop anchors (`1/ε`, `k·p` for `±p`).
#[test]
fn one_loop_cases_match_forcer_exactly_symbolically_and_numerically() {
    assert_loop_count_matches_forcer(1, 45 + 17 + 3);
}

/// 2L: 23 battery cases (kite, dots, odd rank with `-p(3)`, sunset, pinched,
/// negative power), the 39 grid cases and the 3 kite anchors.
#[test]
fn two_loop_cases_match_forcer_exactly_symbolically_and_numerically() {
    assert_loop_count_matches_forcer(2, 23 + 39 + 3);
}

/// 3L: 19 battery cases (LA/BE/NO tops, ISPs, dots, pinched, contracted
/// 5-valent vertices), the 67 grid cases of the gate (45 benchmark cells,
/// 22 dotted sectors) and the 5 anchors (the master `no` and the parents of
/// the ε-line kite `t1star05`, including its `n = 0` chart).
#[test]
fn three_loop_cases_match_forcer_exactly_symbolically_and_numerically() {
    assert_loop_count_matches_forcer(3, 19 + 67 + 5);
}

/// The four-loop cases of `just test_vakint` (M3 D-M3-10; the complete
/// four-loop battery, the anchors and the gate part of the grid run in the
/// full profile): the LA4 top (insertions, carpets and the master `lala`,
/// `-70 ζ(7)` at `ε^0`), a dotted sector of the crossed master `cross` (a
/// Custom owner of the J2 family) and the top of NO1 (the master `no1`, a
/// node Forcer cannot check whole-family). The first and the last are `Oep`
/// at 11 terms (`lala`, `no1` known through `ε^5`), the second has an 11-term
/// value.
const FAST_FOUR_LOOP_CASES: [&str; 3] = ["la4_top", "cross_sector_dot", "no1_top"];

/// A four-loop battery case against its frozen Forcer data, on `engine`
/// under `arithmetic`: checks 1-3 through the case's frozen depth, and at
/// [`I7_TERMS`] terms either RustRed's `Numerical` value equals Forcer's
/// frozen `Symbolic` value (M3 I-7; unit-scale Laurent coefficients, at
/// `Q² = μ² = 1`) or the case fails closed on both backends with the
/// `MissingMasterOrders` (`Oep`) of the fixture.
fn assert_four_loop_case(engine: &Vakint, arithmetic: Option<&str>, case: &ForcerCase) {
    let file = match arithmetic {
        None => FOUR_LOOP_BATTERY.file.to_owned(),
        Some(arithmetic) => format!("{} ({arithmetic})", FOUR_LOOP_BATTERY.file),
    };
    let terms = case_terms(case);
    assert_reduced_exactly(
        engine,
        &settings(RustRedMasterMode::Retain, terms, arithmetic),
        &file,
        case,
    );
    assert_symbolic(engine, terms, arithmetic, &file, case);
    assert_numerical(engine, terms, arithmetic, &file, case);
    let name = format!("{file}: {} at {I7_TERMS} terms", case.name);
    match case.eleven_terms.as_ref().expect("an 11-term value") {
        Finalized::Symbolic(coefficients) => {
            let settings = settings(RustRedMasterMode::Numerical, I7_TERMS, arithmetic);
            let actual = engine
                .evaluate(&settings, case.expression().unwrap().as_view())
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            let e = Atom::var(epsilon(&settings));
            let loops = case.loops().unwrap() as i64;
            let expected = coefficients
                .iter()
                .zip(-loops..)
                .fold(Atom::Zero, |sum, (coefficient, power)| {
                    sum + parse(coefficient) * e.pow(Atom::num(power))
                });
            assert_numerically_equal(engine, &settings, 1, &name, &actual, &expected);
        }
        Finalized::MissingMasterOrders(expected) => {
            let expected = (expected.iter())
                .map(|m| {
                    (
                        m.master.as_str(),
                        m.first_missing_order,
                        m.first_affected_result_order,
                    )
                })
                .collect::<Vec<_>>();
            assert!(!expected.is_empty(), "{name}");
            for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
                let settings = settings(mode, I7_TERMS, arithmetic);
                match engine.evaluate(&settings, case.expression().unwrap().as_view()) {
                    Err(VakintError::MissingMasterOrders(missing)) => {
                        let found = (missing.masters.iter())
                            .map(|m| {
                                (
                                    m.master.as_str(),
                                    m.first_missing_order,
                                    m.first_affected_result_order,
                                )
                            })
                            .collect::<Vec<_>>();
                        assert_eq!(found, expected, "{name} ({mode})");
                    }
                    other => panic!("{name} ({mode}): expected MissingMasterOrders, got {other:?}"),
                }
            }
            assert_eq!(
                case.finalized_at(I7_TERMS).unwrap(),
                *case.eleven_terms.as_ref().unwrap(),
                "{name}: Forcer"
            );
        }
    }
}

/// 4L (M3-V9, on a scope 1-4 library): [`FAST_FOUR_LOOP_CASES`] through
/// [`assert_four_loop_case`]. Without a four-loop library the test reports
/// the skip on stderr; [`rustred_covers_the_library_scope_and_no_more`]
/// then pins that four loops fail closed.
#[test]
fn four_loop_cases_match_forcer_exactly_symbolically_and_numerically() {
    let engine = battery_engine();
    if !battery_loops(&engine).contains(&4) {
        eprintln!("no four-loop RustRed p library ({P_LIBRARY_1_4_VARIABLE} unset): skipped");
        return;
    }
    let battery = FOUR_LOOP_BATTERY.parse();
    for name in FAST_FOUR_LOOP_CASES {
        let case = (battery.cases.iter())
            .find(|case| case.name == name)
            .unwrap_or_else(|| panic!("{}: no case {name}", FOUR_LOOP_BATTERY.file));
        assert_four_loop_case(&engine, None, case);
    }
}

/// The arithmetics besides the default `modular_polynomial` (M2 exit
/// criterion 2).
const OTHER_ARITHMETICS: [&str; 2] = ["exact", "modular_sampled"];

/// The cases of `loops` loops (one to three) that check 1 repeats under
/// [`OTHER_ARITHMETICS`]: the frozen battery (every rank group, dots, ISPs,
/// both masters at three loops), the anchors and the gate part of the grid
/// (M3 I-6), i.e. every case of [`cases`].
fn arithmetic_subset(loops: usize) -> Vec<(&'static str, ForcerCase)> {
    cases(loops)
}

/// Check 1 on [`arithmetic_subset`] under `exact` and `modular_sampled`: every
/// arithmetic reproduces Forcer exactly, as the default does in the per-loop
/// tests. The subset sizes are pinned: battery + grid gate part + anchors
/// per loop count. The four-loop subset (one case per family and the
/// anchors) runs in the full profile.
#[test]
fn the_exact_and_sampled_arithmetics_reproduce_the_batteries_grids_and_anchors() {
    let engine = battery_engine();
    for (loops, expected) in [(1, 45 + 17 + 3), (2, 23 + 39 + 3), (3, 19 + 67 + 5)] {
        let cases = arithmetic_subset(loops);
        assert_eq!(cases.len(), expected, "{loops} loops");
        for arithmetic in OTHER_ARITHMETICS {
            let settings = settings(
                RustRedMasterMode::Retain,
                battery_terms(loops),
                Some(arithmetic),
            );
            assert_eq!(settings.ordered_weight_arithmetic.name(), arithmetic);
            for (file, case) in &cases {
                assert_reduced_exactly(&engine, &settings, &format!("{file} ({arithmetic})"), case);
            }
        }
    }
}

/// The cases per loop count for checks 2 and 3 under every arithmetic: odd
/// and even ranks with dots at one loop, the kite with a rank-3 and a dotted
/// rank-4 numerator at two loops (no master below three loops), at three
/// loops the masters `no` (`no_top`), `t105` (`la_top`) and both
/// (`no_dot`), and at four loops (on a scope 1-4 library) the dotted `cross`
/// sector of [`FAST_FOUR_LOOP_CASES`], through its frozen depth and at 11
/// terms ([`assert_four_loop_case`]).
const FINALIZE_CASES: [&[&str]; 4] = [
    &["bubble_32_n4", "bubble_13_n3", "bubble_21_n1"],
    &["kite", "kite_rank3", "kite_dots_rank4"],
    &["no_top", "la_top", "no_dot"],
    &["cross_sector_dot"],
];

/// Checks 2 and 3 (Symbolic through `ε^(L+2)`, Numerical within `1e-30` at
/// `μ² = 2`) on [`FINALIZE_CASES`] under `exact` and `modular_sampled`, for
/// every loop count of [`battery_loops`]; the default arithmetic finalizes
/// every case in the per-loop tests.
#[test]
fn the_exact_and_sampled_arithmetics_finalize_like_forcer() {
    let engine = battery_engine();
    for loops in battery_loops(&engine) {
        let fixture = BATTERIES.get(loops - 1).unwrap_or(&FOUR_LOOP_BATTERY);
        let battery = fixture.parse();
        for name in FINALIZE_CASES[loops - 1] {
            let case = battery
                .cases
                .iter()
                .find(|case| case.name == *name)
                .unwrap_or_else(|| panic!("{}: no case {name}", fixture.file));
            for arithmetic in OTHER_ARITHMETICS {
                if loops == 4 {
                    assert_four_loop_case(&engine, Some(arithmetic), case);
                    continue;
                }
                let file = format!("{} ({arithmetic})", fixture.file);
                assert_symbolic(&engine, battery_terms(loops), Some(arithmetic), &file, case);
                assert_numerical(&engine, battery_terms(loops), Some(arithmetic), &file, case);
            }
        }
    }
}

/// M3 I-7 on M2's three-loop `Oep` boundary (D-S2-2), the `P3L_NO` top that
/// `p4l_battery.json` freezes with its value at 11 terms: RustRed's
/// `Numerical` value at 11 terms equals Forcer's frozen `Symbolic` one
/// (unit-scale Laurent coefficients `ε^-3..ε^7`, at `Q² = μ² = 1`).
#[test]
fn the_three_loop_oep_boundary_at_eleven_terms_matches_forcer_numerically() {
    let engine = battery_engine();
    let fixture = FOUR_LOOP_BATTERY.parse();
    let [case] = fixture.boundary.as_slice() else {
        panic!("one boundary case, got {}", fixture.boundary.len());
    };
    assert_eq!(case.name, "p3l_no_top");
    let Some(Finalized::Symbolic(coefficients)) = &case.eleven_terms else {
        panic!("{}: no 11-term value", case.name);
    };
    assert_eq!(coefficients.len() as i64, I7_TERMS);
    let settings = settings(RustRedMasterMode::Numerical, I7_TERMS, None);
    let actual = engine
        .evaluate(&settings, case.expression().unwrap().as_view())
        .unwrap_or_else(|error| panic!("{}: {error}", case.name));
    let e = Atom::var(epsilon(&settings));
    let expected = coefficients
        .iter()
        .zip(-3..)
        .fold(Atom::Zero, |sum, (coefficient, power)| {
            sum + parse(coefficient) * e.pow(Atom::num(power))
        });
    assert_numerically_equal(&engine, &settings, 1, &case.name, &actual, &expected);
}

/// The frozen battery cases per loop count whose Wick sign is visible:
/// Forcer's frozen values with and without the sign differ.
const WICK_VISIBLE_CASES: [usize; 3] = [17, 11, 10];

/// The Wick sign `(-1)^(Σa+r)` belongs to the Forcer path (critique-phys SF3),
/// and RustRed's value equals Forcer's with the sign applied exactly once.
/// For every frozen battery case whose sign is visible (Forcer's frozen value
/// with the sign, `ForcerCase::frozen`, differs from `combine_parts` with every
/// sign `+1`, which is also the value with the sign applied twice), RustRed's
/// reduced value is the one with the sign. The number of such cases is
/// pinned per loop count, so that the batteries keep exercising the sign;
/// RustRed is checked on one to three loops (the four-loop cases carry the
/// sign in every comparison with their frozen value).
#[test]
fn rustred_carries_the_wick_sign_of_odd_cases_once() {
    let engine = battery_engine();
    let mut visible = Vec::new();
    for loops in 1..=3 {
        let settings = settings(RustRedMasterMode::Retain, battery_terms(loops), None);
        let epsilon = epsilon(&settings);
        let fixture = &BATTERIES[loops - 1];
        let mut count = 0;
        for case in fixture.parse().cases {
            let name = format!("{}: {}", fixture.file, case.name);
            let Some(input) = case.input().unwrap() else {
                continue;
            };
            let unsigned = input
                .parts
                .iter()
                .map(|part| ForcerPart {
                    sign: 1,
                    ..part.clone()
                })
                .collect::<Vec<_>>();
            let signed = by_scale_power(case.frozen(epsilon).unwrap());
            let unsigned =
                by_scale_power(combine_parts(&unsigned, &case.outputs, epsilon).unwrap());
            if reduced_agree(&signed, &unsigned, epsilon).is_ok() {
                continue;
            }
            count += 1;
            let rustred = rustred_reduced(&engine, &settings, &case)
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            reduced_agree(&rustred, &signed, epsilon).unwrap_or_else(|difference| {
                panic!("{name}: RustRed vs Forcer with its Wick sign: {difference}")
            });
        }
        visible.push(count);
    }
    assert_eq!(visible, WICK_VISIBLE_CASES, "Wick-visible battery cases");
}

/// The battery engine covers exactly its library's scope ([`battery_loops`]:
/// `1..=3` built in, `1..=4` with the M3 package), and an engine without a
/// four-loop library never reduces four loops: for the default engine
/// (`Vakint::new()`, whose built-in library stops at three loops until M3
/// switches the default), `has_p_data(4)` is false and every case of the
/// four-loop battery fails closed on the RustRed path (an unsupported
/// family, never a value).
#[test]
fn rustred_covers_the_library_scope_and_no_more() {
    let engine = battery_engine();
    let scope = battery_loops(&engine);
    let expected = if attached_library().is_some() { 4 } else { 3 };
    assert_eq!(scope, 1..=expected);
    for loops in 1..=4 {
        assert_eq!(
            has_p_data(&engine, loops),
            scope.contains(&loops),
            "{loops} loops"
        );
    }
    let default = Vakint::new().unwrap();
    assert_eq!(battery_loops(&default), 1..=3);
    assert!(!has_p_data(&default, 4));
    let settings = settings(RustRedMasterMode::Retain, battery_terms(4), None);
    for case in FOUR_LOOP_BATTERY.parse().cases {
        assert_eq!(case.loops().unwrap(), 4, "{}", case.name);
        match rustred_reduced(&default, &settings, &case) {
            Err(VakintError::RustRedEvaluation(
                RustRedEvaluationError::UnsupportedMatchedFamily { .. },
            )) => {}
            other => panic!(
                "{}: expected an unsupported family, got {other:?}",
                case.name
            ),
        }
    }
}

/// The gate part of the grids is what [`in_gate`] says: the three-loop grid
/// keeps its 45 benchmark cases and 22 dotted sectors and leaves its 350
/// labelled sectors to Tier C; the other grids are whole.
#[test]
fn the_gate_part_of_the_grids_is_pinned() {
    let counts = GRIDS
        .iter()
        .map(|fixture| fixture.parse().cases.iter().filter(|c| in_gate(c)).count())
        .collect::<Vec<_>>();
    assert_eq!(counts, [17, 39, 45 + 22]);
}
