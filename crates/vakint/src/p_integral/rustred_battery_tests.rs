//! Tier A of the Forcer battery (M2-V5; M2 plan, section 8, gate G7): the
//! RustRed p path against Forcer's frozen outputs, without FORM.
//!
//! Every case of the frozen batteries `tests/fixtures/forcer/p{1,2,3}l_battery.json`
//! (45 + 23 + 19 = 87), of the gate part of the grids
//! `p{1,2,3}l_grid.json` ([`in_gate`]) and of the anchors `p_anchors.json` is
//! checked three ways:
//!
//! 1. exact: RustRed's reduced values (the request level below the finalizer,
//!    `rustred_evaluation::p_integral::reduce_requests`) equal Forcer's
//!    (`forcer_input` parts and `combine_parts` on the frozen outputs, i.e.
//!    `ForcerCase::frozen`), scale power by scale power, in the
//!    `Π Γ(1+mε)^e × Π PMaster` tag basis (`compare::tag_basis_difference`);
//! 2. symbolic: `Vakint::evaluate` with RustRed in `Symbolic` mode through
//!    `ε^(L+2)` equals the frozen value finalized like a Forcer evaluation;
//! 3. numerical: the same in `Numerical` mode at 40 digits, within `1e-30`
//!    at `Q² = μ² = 1`.
//!
//! All of it runs under the default Ordered arithmetic (`modular_polynomial`,
//! the `Z_p(d)` replay vakint users get); the exact check of the frozen
//! batteries is repeated under the `exact` arithmetic (coordinator note
//! 2026-11-07 09:40Z). Four loops stay with Forcer: the engine has no p
//! library for them, and a four-loop case fails closed on the RustRed path.
//!
//! The engine of the battery is [`battery_engine`], covering
//! [`BATTERY_LOOPS`]. STAGE 1 (until M2-V4 ships the built-in 1L-3L library):
//! RustRed's one-loop reference library, injected; the two- and three-loop
//! tests are ignored. STAGE 2: the built-in library of `Vakint::new()`,
//! `BATTERY_LOOPS = 1..=3`, no test ignored. The full grids and the
//! benchmark cells through vakint's public API are Tier C
//! (`tests/p_integral_rustred_full.rs`); live FORM comparisons are Tier B
//! (`tests/p_integral_rustred_tests.rs`).

use std::collections::BTreeMap;
use std::ops::RangeInclusive;

use symbolica::atom::{Atom, AtomCore, Symbol};
use symbolica::function;

use super::compare::tag_basis_difference;
use super::forcer::fixtures::{ForcerCase, ForcerFixture, registry_family};
use super::{PropagatorRequest, ReducedPropagatorValue};
use crate::rustred_evaluation::p_integral::{has_p_data, reduce_requests};
use crate::rustred_evaluation::{RustRedEvaluationError, RustRedMasterMode, combine_exact};
use crate::utils::simplify_real;
use crate::utils::vakint_macros::{vk_parse, vk_symbol};
use crate::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions, Vakint,
    VakintError, VakintSettings,
};

/// The loop counts [`battery_engine`] reduces with RustRed (STAGE 1: the
/// one-loop reference library).
const BATTERY_LOOPS: RangeInclusive<usize> = 1..=1;

/// The engine of the battery: RustRed's one-loop reference library injected
/// (STAGE 1). STAGE 2 switches to `Vakint::new()` and its built-in library.
fn battery_engine() -> Vakint {
    let mut engine = Vakint::new().unwrap();
    engine
        .inject_p_library(rustred::propagator::testing::one_loop_reference_library().unwrap())
        .unwrap();
    engine
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
    assert_eq!(
        rustred.keys().collect::<Vec<_>>(),
        forcer.keys().collect::<Vec<_>>(),
        "{name}: scale powers"
    );
    for (w, value) in &rustred {
        let difference = tag_basis_difference(value, &forcer[w], epsilon)
            .unwrap_or_else(|error| panic!("{name}, w = {w}: {error}"));
        assert!(
            difference.is_empty(),
            "{name}, w = {w}: RustRed {value} vs Forcer {}: {difference:?}",
            forcer[w]
        );
    }
}

/// Check 2: the case through `Vakint::evaluate` (matching, the p backend and
/// the shared finalizer) in `Symbolic` mode is its frozen value finalized like
/// a Forcer evaluation (which returns through `simplify_real`).
fn assert_symbolic(engine: &Vakint, loops: usize, file: &str, case: &ForcerCase) {
    let settings = settings(RustRedMasterMode::Symbolic, loops, None);
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

/// Check 3: the same in `Numerical` mode, both sides evaluated at
/// `Q² = μ² = 1` to 40 digits and equal within `1e-30`.
fn assert_numerical(engine: &Vakint, loops: usize, file: &str, case: &ForcerCase) {
    let name = format!("{file}: {}", case.name);
    let settings = settings(RustRedMasterMode::Numerical, loops, None);
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
    let at_unit_scale = |value: &Atom| {
        let x = Atom::var(crate::symbols::S.x_);
        let p = function!(crate::symbols::S.p, &x);
        let unit = value
            .replace(function!(crate::symbols::S.dot, &p, &p).to_pattern())
            .with(Atom::num(-1).to_pattern())
            .replace(function!(Symbol::LOG, &x).to_pattern())
            .with(Atom::Zero.to_pattern())
            .expand();
        engine
            .numerical_evaluation(
                &settings,
                unit.as_view(),
                &Default::default(),
                &Default::default(),
                None,
            )
            .unwrap_or_else(|error| panic!("{name}: {unit}: {error}"))
            .0
    };
    let (actual, expected) = (at_unit_scale(&actual), at_unit_scale(&expected));
    let (equal, diagnostic) = actual.does_approx_match(&expected, None, 1e-30, 5.0);
    assert!(equal, "{name}: {diagnostic}");
}

/// The three checks on every case of `loops` loops ([`cases`]), whose number
/// is pinned: `battery + grid gate part + anchors`.
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
        assert_symbolic(&engine, loops, file, case);
        assert_numerical(&engine, loops, file, case);
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
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn two_loop_cases_match_forcer_exactly_symbolically_and_numerically() {
    assert_loop_count_matches_forcer(2, 23 + 39 + 3);
}

/// 3L: 19 battery cases (LA/BE/NO tops, ISPs, dots, pinched, contracted
/// 5-valent vertices), the 67 grid cases of the gate (45 benchmark cells,
/// 22 dotted sectors) and the 5 anchors (the master `no` and the parents of
/// the ε-line kite `t1star05`, including its `n = 0` chart).
#[test]
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn three_loop_cases_match_forcer_exactly_symbolically_and_numerically() {
    assert_loop_count_matches_forcer(3, 19 + 67 + 5);
}

/// The exact check of every frozen battery case the engine covers, repeated
/// under the `exact` Ordered arithmetic: the default `modular_polynomial` and
/// the exact replay both reproduce Forcer.
#[test]
fn the_exact_arithmetic_reproduces_the_batteries_too() {
    let engine = battery_engine();
    let mut checked = 0;
    for loops in BATTERY_LOOPS {
        let exact = settings(RustRedMasterMode::Retain, loops, Some("exact"));
        assert_eq!(exact.ordered_weight_arithmetic.name(), "exact");
        assert_eq!(
            settings(RustRedMasterMode::Retain, loops, None)
                .ordered_weight_arithmetic
                .name(),
            "modular_polynomial"
        );
        let fixture = &BATTERIES[loops - 1];
        for case in fixture.parse().cases {
            assert_reduced_exactly(&engine, &exact, fixture.file, &case);
            checked += 1;
        }
    }
    assert_eq!(
        checked,
        BATTERIES[..*BATTERY_LOOPS.end()]
            .iter()
            .map(|fixture| fixture.cases)
            .sum::<usize>()
    );
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
