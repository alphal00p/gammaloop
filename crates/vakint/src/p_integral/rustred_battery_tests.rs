//! Tier A of the Forcer battery (M2-V5; M2 plan, section 8, gate G7): the
//! RustRed p path against Forcer's frozen outputs, without FORM.
//!
//! Every case of the frozen batteries `tests/fixtures/forcer/p{1,2,3}l_battery.json`
//! (45 + 23 + 19 = 87) and of the anchors `p_anchors.json` (11) is checked
//! three ways, every case of the gate part of the grids `p{1,2,3}l_grid.json`
//! ([`in_gate`]) by check 1 (Tier C runs every grid case through all three
//! modes of the public API):
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
//! criterion 2), repeat check 1 on the frozen batteries and the anchors
//! ([`arithmetic_subset`]) and checks 2 and 3 on [`FINALIZE_CASES`]. Four
//! loops stay with Forcer: the engine has no p library for them, and a
//! four-loop case fails closed on the RustRed path.
//!
//! The engine of the battery is [`battery_engine`], `Vakint::new()` with
//! its built-in 1L-3L p library (M2-V4), covering [`BATTERY_LOOPS`]. The full
//! grids and the benchmark cells through vakint's public API are Tier C
//! (`tests/p_integral_rustred_full.rs`); live FORM comparisons are Tier B
//! (`tests/p_integral_rustred_tests.rs`).

use std::collections::{BTreeMap, HashMap};
use std::ops::RangeInclusive;

use symbolica::atom::{Atom, AtomCore, Symbol};
use symbolica::domains::float::Float;

use super::compare::tag_basis_difference;
use super::forcer::fixtures::{ForcerCase, ForcerFixture, registry_family};
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

/// The loop counts [`battery_engine`] reduces with RustRed: the scope of the
/// built-in p library.
const BATTERY_LOOPS: RangeInclusive<usize> = 1..=3;

/// The engine of the battery: `Vakint::new()`, whose built-in p library
/// (M2-V4) covers [`BATTERY_LOOPS`].
fn battery_engine() -> Vakint {
    Vakint::new().unwrap()
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
    cases: 10,
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

/// RustRed alone in `mode`, G-scheme, through `ε^(L+2)` for `loops` loops (a
/// FORM call fails loudly), 40 digits, one core; the settings' default
/// Ordered arithmetic (`modular_polynomial`) or `arithmetic`.
fn settings(mode: RustRedMasterMode, loops: usize, arithmetic: Option<&str>) -> VakintSettings {
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions { master_mode: mode },
        )]),
        form_exe_path: "/rustred-battery-must-not-call-form".into(),
        integral_normalization_factor: LoopNormalizationFactor::GScheme,
        number_of_terms_in_epsilon_expansion: 2 * loops as i64 + 3,
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
/// the shared finalizer) in `Symbolic` mode is its frozen value finalized like
/// a Forcer evaluation (which returns through `simplify_real`).
fn assert_symbolic(
    engine: &Vakint,
    loops: usize,
    arithmetic: Option<&str>,
    file: &str,
    case: &ForcerCase,
) {
    let settings = settings(RustRedMasterMode::Symbolic, loops, arithmetic);
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

/// The point of check 3: `p(1..4) = (0,0,0,1)`, so `Q² = -p² = 1`, `μ² = 2`
/// and `mass = 3` (the free coefficient of the two-loop battery case
/// `kite_coefficients`), at the precision of `settings`.
type NumericalPoint = (
    HashMap<String, Float, ahash::RandomState>,
    HashMap<usize, Momentum, ahash::RandomState>,
);

fn numerical_point(settings: &VakintSettings) -> NumericalPoint {
    let bits = settings.get_binary_precision();
    let float = |value: i64| Float::with_val(bits, value);
    let params = HashMap::from_iter([
        (settings.mu_r_sq_symbol.clone(), float(2)),
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

/// Check 3: the same in `Numerical` mode, both sides evaluated to 40 digits
/// at [`numerical_point`] (the logarithms of `μ²/Q²` count) and equal within
/// `1e-30`.
fn assert_numerical(
    engine: &Vakint,
    loops: usize,
    arithmetic: Option<&str>,
    file: &str,
    case: &ForcerCase,
) {
    let name = format!("{file}: {}", case.name);
    let settings = settings(RustRedMasterMode::Numerical, loops, arithmetic);
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
    let (params, externals) = numerical_point(&settings);
    let at_point = |value: &Atom| {
        engine
            .numerical_evaluation(
                &settings,
                value.as_view(),
                &params,
                &Default::default(),
                Some(&externals),
            )
            .unwrap_or_else(|error| panic!("{name}: {value}: {error}"))
            .0
    };
    let (actual, expected) = (at_point(&actual), at_point(&expected));
    let (equal, diagnostic) = actual.does_approx_match(&expected, None, 1e-30, 5.0);
    assert!(equal, "{name}: {diagnostic}");
}

/// Check 1 on every case of `loops` loops ([`cases`]), whose number is
/// pinned (`battery + grid gate part + anchors`), and checks 2 and 3 on the
/// battery and anchor cases (the grid cases go through every mode in Tier C).
fn assert_loop_count_matches_forcer(loops: usize, expected_cases: usize) {
    let engine = battery_engine();
    assert!(
        has_p_data(&engine, loops),
        "no RustRed p library for {loops} loops"
    );
    let cases = cases(loops);
    assert_eq!(cases.len(), expected_cases, "{loops}-loop cases");
    let retain = settings(RustRedMasterMode::Retain, loops, None);
    for (file, case) in &cases {
        assert_reduced_exactly(&engine, &retain, file, case);
        if *file != GRIDS[loops - 1].file {
            assert_symbolic(&engine, loops, None, file, case);
            assert_numerical(&engine, loops, None, file, case);
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

/// The arithmetics besides the default `modular_polynomial` (M2 exit
/// criterion 2).
const OTHER_ARITHMETICS: [&str; 2] = ["exact", "modular_sampled"];

/// The cases of `loops` loops that check 1 repeats under
/// [`OTHER_ARITHMETICS`]: the frozen battery (every rank group, dots, ISPs,
/// both masters at three loops) and the anchors.
fn arithmetic_subset(loops: usize) -> Vec<(&'static str, ForcerCase)> {
    let mut cases = Vec::new();
    for fixture in [&BATTERIES[loops - 1], &ANCHORS] {
        for case in fixture.parse().cases {
            if case.loops().unwrap() == loops {
                cases.push((fixture.file, case));
            }
        }
    }
    cases
}

/// Check 1 on [`arithmetic_subset`] under `exact` and `modular_sampled`: every
/// arithmetic reproduces Forcer exactly, as the default does in the per-loop
/// tests. The subset sizes are pinned: battery + anchors per loop count.
#[test]
fn the_exact_and_sampled_arithmetics_reproduce_the_batteries_and_anchors() {
    let engine = battery_engine();
    let anchors = [3, 3, 5];
    for loops in BATTERY_LOOPS {
        let cases = arithmetic_subset(loops);
        assert_eq!(
            cases.len(),
            BATTERIES[loops - 1].cases + anchors[loops - 1],
            "{loops} loops"
        );
        for arithmetic in OTHER_ARITHMETICS {
            let settings = settings(RustRedMasterMode::Retain, loops, Some(arithmetic));
            assert_eq!(settings.ordered_weight_arithmetic.name(), arithmetic);
            for (file, case) in &cases {
                assert_reduced_exactly(&engine, &settings, &format!("{file} ({arithmetic})"), case);
            }
        }
    }
}

/// Three cases per loop count for checks 2 and 3 under every arithmetic: odd
/// and even ranks with dots at one loop, the kite with a rank-3 and a dotted
/// rank-4 numerator at two loops (no master below three loops), and at three
/// loops the masters `no` (`no_top`), `t105` (`la_top`) and both (`no_dot`).
const FINALIZE_CASES: [[&str; 3]; 3] = [
    ["bubble_32_n4", "bubble_13_n3", "bubble_21_n1"],
    ["kite", "kite_rank3", "kite_dots_rank4"],
    ["no_top", "la_top", "no_dot"],
];

/// Checks 2 and 3 (Symbolic through `ε^(L+2)`, Numerical within `1e-30` at
/// `μ² = 2`) on [`FINALIZE_CASES`] under `exact` and `modular_sampled`; the
/// default arithmetic finalizes every case in the per-loop tests.
#[test]
fn the_exact_and_sampled_arithmetics_finalize_like_forcer() {
    let engine = battery_engine();
    for loops in BATTERY_LOOPS {
        let fixture = &BATTERIES[loops - 1];
        let battery = fixture.parse();
        for name in FINALIZE_CASES[loops - 1] {
            let case = battery
                .cases
                .iter()
                .find(|case| case.name == name)
                .unwrap_or_else(|| panic!("{}: no case {name}", fixture.file));
            for arithmetic in OTHER_ARITHMETICS {
                let file = format!("{} ({arithmetic})", fixture.file);
                assert_symbolic(&engine, loops, Some(arithmetic), &file, case);
                assert_numerical(&engine, loops, Some(arithmetic), &file, case);
            }
        }
    }
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
/// RustRed is checked on [`BATTERY_LOOPS`].
#[test]
fn rustred_carries_the_wick_sign_of_odd_cases_once() {
    let engine = battery_engine();
    let mut visible = Vec::new();
    for loops in 1..=3 {
        let settings = settings(RustRedMasterMode::Retain, loops, None);
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
            if BATTERY_LOOPS.contains(&loops) {
                let rustred = rustred_reduced(&engine, &settings, &case)
                    .unwrap_or_else(|error| panic!("{name}: {error}"));
                reduced_agree(&rustred, &signed, epsilon).unwrap_or_else(|difference| {
                    panic!("{name}: RustRed vs Forcer with its Wick sign: {difference}")
                });
            }
        }
        visible.push(count);
    }
    assert_eq!(visible, WICK_VISIBLE_CASES, "Wick-visible battery cases");
}

/// The engine covers exactly [`BATTERY_LOOPS`]; four loops never reach
/// RustRed: `has_p_data(4)` is false for the battery engine and for a bare
/// engine, and the four-loop battery fails closed on the RustRed path (an
/// unsupported family, never a value).
#[test]
fn rustred_covers_the_battery_loops_and_never_four_loops() {
    let engine = battery_engine();
    for loops in 1..=3 {
        assert_eq!(
            has_p_data(&engine, loops),
            BATTERY_LOOPS.contains(&loops),
            "{loops} loops"
        );
    }
    assert!(!has_p_data(&engine, 4));
    assert!(!has_p_data(&Vakint::new().unwrap(), 4));
    let settings = settings(RustRedMasterMode::Retain, 4, None);
    for case in FOUR_LOOP_BATTERY.parse().cases {
        assert_eq!(case.loops().unwrap(), 4, "{}", case.name);
        match rustred_reduced(&engine, &settings, &case) {
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
