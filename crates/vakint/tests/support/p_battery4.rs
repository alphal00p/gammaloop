//! Tier A of the Forcer battery at four loops (M3-V9; M3 plan, gate G7 and
//! D-M3-10): RustRed on M3's scope 1-4 p library against the frozen Forcer
//! fixtures of M3-V7, through vakint's public API, without FORM.
//!
//! The cases are every case of `tests/fixtures/forcer/p4l_battery.json`
//! (54, over the 11 four-loop families), the anchors `p4l_anchors.json` (17:
//! every master of Forcer's four-loop basis, `la4_top`, the NO3 and bubu
//! tops) and the gate part of `p4l_grid.json` ([`in_gate4`]: the 55 sampled
//! sectors; the 165 benchmark tops run in Tier C). Each agrees with its
//! frozen Forcer value ([`assert_case_agrees`], `support/p_backends.rs`
//! `agreement_on`): `Retain` exactly in the `Π Γ(1+mε)^e × Π PMaster` tag
//! basis, `Symbolic` exactly and `Numerical` at 40 digits within `1e-30`
//! (at `Q² = μ² = 1`), through the case's frozen depth ([`case_terms`]: the
//! deepest Laurent depth at which Forcer's value finalizes). A battery case
//! is also checked at 11 terms ([`assert_eleven_terms`]): RustRed's
//! `Numerical` value equals Forcer's frozen `Symbolic` value (M3 I-7), or
//! both backends fail closed with the case's `MissingMasterOrders` (`Oep`).
//!
//! The engine ([`p_battery4_engine`]) is `Vakint::new()` with the package
//! that `VAKINT_TEST_P_LIBRARY_1_4` names injected (`support/p_library4.rs`);
//! without it every four-loop check reports a skip on stderr.
//!
//! Include it after its sibling support files, whose items it uses:
//!
//! ```text
//! #[path = "support/p_integrals.rs"] mod p_integrals;
//! #[path = "support/p_rustred.rs"] mod p_rustred;
//! #[path = "support/p_backends.rs"] mod p_backends;
//! #[path = "support/p_fixtures.rs"] mod p_fixtures;
//! #[path = "support/p_library4.rs"] mod p_library4;
//! #[path = "support/p_battery4.rs"] mod p_battery4;
//! ```

#![allow(dead_code)]

use std::collections::BTreeSet;
use std::sync::LazyLock;

use symbolica::atom::{Atom, AtomCore};
use vakint::p_integral::compare::tag_basis;
use vakint::p_integral::forcer_fixtures::{Finalized, ForcerCase, ForcerFixture, I7_TERMS};
use vakint::utils::simplify_real;
use vakint::{
    LoopNormalizationFactor, RustRedMasterMode, Vakint, VakintError, VakintSettings, vakint_symbol,
};

use super::p_backends::{
    MODES, agreement_on, rustred_p_settings, rustred_p_settings_with_arithmetic, values_agree,
};
use super::p_fixtures::battery;
use super::p_library4::{
    P_LIBRARY_1_4_VARIABLE, PLibraryRoute, engine_with_p_library, p_library_1_4, p_library_scope,
};

/// The four-loop battery, its case count and its `Oep` cases at 11 terms.
pub const BATTERY4: &str = "p4l_battery.json";
pub const BATTERY4_CASES: usize = 54;
pub const BATTERY4_OEP_CASES: usize = 26;
/// The four-loop anchors and their count.
pub const ANCHORS4: &str = "p4l_anchors.json";
pub const ANCHORS4_CASES: usize = 17;
/// The four-loop grid, its case count and the size of its gate part.
pub const GRID4: &str = "p4l_grid.json";
pub const GRID4_CASES: usize = 220;
pub const GRID4_GATE_CASES: usize = 55;

/// The masters of Forcer's four-loop basis (vakint's `masters.rs`, M3-T17's
/// `P_BASIS_MASTERS_1_TO_4`), each reached by the four-loop checks.
pub const FOUR_LOOP_MASTERS: [&str; 16] = [
    "bebe", "cross", "fastar2", "haha", "lala", "no", "no1", "no2", "no6", "nono", "nostar5",
    "nostar6", "t105", "t124", "t145", "t155",
];

/// The other Ordered arithmetics (M3 G7: on [`arithmetic_subset4`]).
pub const OTHER_ARITHMETICS: [&str; 2] = ["exact", "modular_sampled"];

static P_BATTERY4_ENGINE: LazyLock<Option<Vakint>> = LazyLock::new(|| {
    let directory = p_library_1_4()?;
    assert_eq!(
        p_library_scope(&directory),
        [1, 4],
        "{P_LIBRARY_1_4_VARIABLE}={}",
        directory.display()
    );
    Some(engine_with_p_library(&directory, PLibraryRoute::Injected).engine)
});

/// The engine of the four-loop battery, shared by the tests of a binary:
/// `Vakint::new()` with M3's scope 1-4 package injected; `None` (a skip of
/// `what` reported on stderr) when `VAKINT_TEST_P_LIBRARY_1_4` is unset.
pub fn p_battery4_engine(what: &str) -> Option<&'static Vakint> {
    let engine = P_BATTERY4_ENGINE.as_ref();
    if engine.is_none() {
        eprintln!("{P_LIBRARY_1_4_VARIABLE} is unset: skipped {what}");
    }
    engine
}

/// A four-loop fixture file with its case count pinned.
pub fn fixture4(file: &str, cases: usize) -> ForcerFixture {
    let fixture = battery(file);
    assert_eq!(fixture.cases.len(), cases, "{file}");
    assert!(
        fixture.cases.iter().all(|case| case.loops().unwrap() == 4),
        "{file}: four loops"
    );
    fixture
}

/// The grid cases of Tier A: the sampled labelled sectors (rank 0 and 2) and
/// the dotted sectors; the benchmark tops (`bench_*`, up to rank 8 and two
/// dots) run through every mode in Tier C.
pub fn in_gate4(case: &ForcerCase) -> bool {
    !case.name.starts_with("bench_")
}

/// The gate part of the four-loop grid ([`in_gate4`]), its size pinned.
pub fn grid4_gate() -> Vec<ForcerCase> {
    let gate = fixture4(GRID4, GRID4_CASES)
        .cases
        .into_iter()
        .filter(in_gate4)
        .collect::<Vec<_>>();
    assert_eq!(gate.len(), GRID4_GATE_CASES, "{GRID4}: gate part");
    gate
}

/// The cases of the other arithmetics ([`OTHER_ARITHMETICS`]): one battery
/// case per family, its dotted case (`<family>_dot`, a reduction through the
/// family's owners), and every anchor (each master of the basis), as
/// `(file, case)`.
pub fn arithmetic_subset4() -> Vec<(&'static str, ForcerCase)> {
    let battery = fixture4(BATTERY4, BATTERY4_CASES);
    let families = (battery.cases.iter())
        .map(|case| case.family.clone())
        .collect::<BTreeSet<_>>();
    assert_eq!(families.len(), 11, "{BATTERY4}: families");
    let mut subset = families
        .iter()
        .map(|family| {
            let name = format!("{}_dot", family.trim_start_matches("P4L_").to_lowercase());
            let case = (battery.cases.iter())
                .find(|case| case.name == name && case.family == *family)
                .unwrap_or_else(|| panic!("{BATTERY4}: no case {name}"));
            (BATTERY4, case.clone())
        })
        .collect::<Vec<_>>();
    let anchors = fixture4(ANCHORS4, ANCHORS4_CASES).cases;
    subset.extend(anchors.into_iter().map(|case| (ANCHORS4, case)));
    subset
}

/// Part `part` of `parts` of `cases`, at most `size` cases each; the parts
/// cover `cases` exactly once.
pub fn part_of<T>(cases: &[T], part: usize, parts: usize, size: usize) -> &[T] {
    assert_eq!(cases.len().div_ceil(size), parts, "parts of {size}");
    &cases[part * size..cases.len().min((part + 1) * size)]
}

/// The Laurent depth a case is checked through: its frozen depth (the
/// deepest at which Forcer's value finalizes, `p4l_battery.json`,
/// `p4l_grid.json`), else that depth computed from the frozen outputs
/// (`ForcerCase::deepest`; the anchors store none).
pub fn case_terms(case: &ForcerCase) -> i64 {
    match &case.depth {
        Some(depth) => depth.terms,
        None => {
            (case.deepest())
                .unwrap_or_else(|error| panic!("{}: depth: {error}", case.name))
                .terms
        }
    }
}

/// RustRed alone in `mode`, G-scheme, `terms` Laurent terms, one core: the
/// default Ordered arithmetic, or `arithmetic`.
pub fn settings4(mode: RustRedMasterMode, terms: i64, arithmetic: Option<&str>) -> VakintSettings {
    match arithmetic {
        None => rustred_p_settings(mode, LoopNormalizationFactor::GScheme, terms),
        Some(arithmetic) => rustred_p_settings_with_arithmetic(
            mode,
            LoopNormalizationFactor::GScheme,
            terms,
            arithmetic,
            1,
        ),
    }
}

/// The Forcer side of a frozen case in `mode` through `terms`: the frozen
/// output finalized like a Forcer evaluation; in `Numerical` mode the exact
/// (`Symbolic`) value, which `values_agree` evaluates at 40 digits.
pub fn frozen_value(
    case: &ForcerCase,
    mode: RustRedMasterMode,
    terms: i64,
) -> Result<Atom, String> {
    let mode = match mode {
        RustRedMasterMode::Numerical => RustRedMasterMode::Symbolic,
        mode => mode,
    };
    case.finalized(&settings4(mode, terms, None), mode)
        .map(|value| simplify_real(value.as_view()))
        .map_err(|error| format!("Forcer fixture {}: {error}", case.name))
}

/// Runs `check` and reports its seconds and peak RSS on stderr as
/// `[tier A 4L] <what>: <s> s, peak <GB> GB` (the evidence of G7 records them
/// per case). The peak is the process's resident high-water mark, reset
/// before `check` (Linux `clear_refs`); elsewhere it reads `?`.
pub fn timed<T>(what: impl std::fmt::Display, check: impl FnOnce() -> T) -> T {
    let _ = std::fs::write("/proc/self/clear_refs", "5");
    let started = std::time::Instant::now();
    let value = check();
    let peak = std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|status| {
            let line = status
                .lines()
                .find_map(|line| line.strip_prefix("VmHWM:"))?;
            let kilobytes = line
                .trim()
                .trim_end_matches("kB")
                .trim()
                .parse::<f64>()
                .ok()?;
            Some(format!("{:.2}", kilobytes / 1e6))
        })
        .unwrap_or_else(|| "?".into());
    eprintln!(
        "[tier A 4L] {what}: {:.1} s, peak {peak} GB",
        started.elapsed().as_secs_f64()
    );
    value
}

/// The masters of a `Retain` value (its tag basis).
pub fn masters_of(value: &Atom) -> BTreeSet<String> {
    let epsilon = vakint_symbol!(&VakintSettings::default().epsilon_symbol);
    tag_basis(value, epsilon)
        .unwrap_or_else(|error| panic!("tag basis of {value}: {error}"))
        .into_keys()
        .flat_map(|(_, masters)| masters.into_iter().map(|(name, _)| name))
        .collect()
}

/// `case` of `file` agrees with its frozen Forcer value on `engine` under
/// `arithmetic` (`None`: the default) in every mode through [`case_terms`];
/// returns the masters RustRed's value carries, reached in `Symbolic` and
/// `Numerical` mode by this agreement.
pub fn assert_case_agrees(
    engine: &Vakint,
    file: &str,
    case: &ForcerCase,
    arithmetic: Option<&str>,
) -> BTreeSet<String> {
    let terms = case_terms(case);
    let values = agreement_on(
        engine,
        &case.expression().unwrap(),
        |mode| settings4(mode, terms, arithmetic),
        |mode| frozen_value(case, mode, terms),
    )
    .unwrap_or_else(|difference| {
        panic!(
            "{file}: {} ({terms} terms, {}):\n{difference}",
            case.name,
            arithmetic.unwrap_or("default arithmetic")
        )
    });
    assert_eq!(MODES[0], RustRedMasterMode::Retain);
    masters_of(&values[0])
}

/// The masters of a `MissingMasterOrders`, with the orders that identify the
/// failure.
fn missing_orders(missing: &vakint::MissingMasterOrders) -> Vec<(String, i64, i64)> {
    (missing.masters.iter())
        .map(|m| {
            (
                m.master.clone(),
                m.first_missing_order,
                m.first_affected_result_order,
            )
        })
        .collect()
}

/// A battery case at [`I7_TERMS`] terms on `engine`: with a frozen 11-term
/// value, RustRed's `Numerical` value equals it at 40 digits (M3 I-7: the
/// unit-scale Laurent coefficients of Forcer's `Symbolic` value, evaluated
/// like RustRed's at `Q² = μ² = 1`); for an `Oep` case, RustRed fails closed
/// in `Symbolic` and `Numerical` mode with the frozen `MissingMasterOrders`,
/// which Forcer's frozen value gives too. Returns whether it is an `Oep`
/// case.
pub fn assert_eleven_terms(engine: &Vakint, file: &str, case: &ForcerCase) -> bool {
    let name = format!("{file}: {} at {I7_TERMS} terms", case.name);
    let input = case.expression().unwrap();
    match case.eleven_terms.as_ref() {
        Some(Finalized::Symbolic(coefficients)) => {
            let settings = settings4(RustRedMasterMode::Numerical, I7_TERMS, None);
            let actual = engine
                .evaluate(&settings, input.as_view())
                .unwrap_or_else(|error| panic!("{name}: {error}"));
            let epsilon = Atom::var(vakint_symbol!(&settings.epsilon_symbol));
            let loops = case.loops().unwrap() as i64;
            assert_eq!(coefficients.len() as i64, I7_TERMS, "{name}");
            let expected = (coefficients.iter().zip(-loops..)).fold(
                Atom::Zero,
                |sum, (coefficient, power)| {
                    sum + super::p_integrals::parse(coefficient) * epsilon.pow(Atom::num(power))
                },
            );
            values_agree(&actual, &expected, RustRedMasterMode::Numerical)
                .unwrap_or_else(|difference| panic!("{name}: {difference}"));
            false
        }
        Some(Finalized::MissingMasterOrders(frozen)) => {
            let expected = (frozen.iter())
                .map(|m| {
                    (
                        m.master.clone(),
                        m.first_missing_order,
                        m.first_affected_result_order,
                    )
                })
                .collect::<Vec<_>>();
            assert!(!expected.is_empty(), "{name}");
            for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
                match engine.evaluate(&settings4(mode, I7_TERMS, None), input.as_view()) {
                    Err(VakintError::MissingMasterOrders(missing)) => {
                        assert_eq!(missing_orders(&missing), expected, "{name} ({mode})");
                    }
                    other => panic!("{name} ({mode}): expected MissingMasterOrders, got {other:?}"),
                }
            }
            match case.finalized(
                &settings4(RustRedMasterMode::Symbolic, I7_TERMS, None),
                RustRedMasterMode::Symbolic,
            ) {
                Err(VakintError::MissingMasterOrders(missing)) => {
                    assert_eq!(missing_orders(&missing), expected, "{name} (Forcer)");
                }
                other => panic!("{name} (Forcer): expected MissingMasterOrders, got {other:?}"),
            }
            true
        }
        None => panic!("{name}: no 11-term value"),
    }
}
