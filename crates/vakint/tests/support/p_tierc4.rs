//! Shared helpers of the four-loop Tiers B and C (M3-V9c; M3 plan, gate G7):
//! the RustRed p backend on M3's scope 1-4 package against Forcer, live or
//! frozen, through vakint's public API, and RustRed's routing and DOT
//! invariance at four loops.
//!
//! An agreement means what it means in `support/p_backends.rs` (`Retain`
//! exactly in the tag basis, `Symbolic` exactly, `Numerical` at 40 digits
//! within `1e-30` at `Q² = μ² = 1`), on an engine with the package that
//! `VAKINT_TEST_P_LIBRARY_1_4` names ([`p_library4_engine`]; without it the
//! four-loop checks report a skip on stderr). The RustRed side runs first,
//! so that a capped run ([`capped_agreement`]) spends no Forcer time on a
//! cell it does not finish.
//!
//! Include it after its sibling support files, whose items it uses:
//!
//! ```text
//! #[path = "support/p_integrals.rs"] mod p_integrals;
//! #[path = "support/p_rustred.rs"] mod p_rustred;
//! #[path = "support/p_backends.rs"] mod p_backends;
//! #[path = "support/p_library4.rs"] mod p_library4;
//! #[path = "support/p_tierc4.rs"] mod p_tierc4;
//! ```

#![allow(dead_code)]

use std::time::{Duration, Instant};

use symbolica::atom::{Atom, AtomCore};
use symbolica::id::Replacement;
use vakint::symbols::S;
use vakint::{LoopNormalizationFactor, RustRedMasterMode, Vakint, VakintError, VakintSettings};

use super::p_backends::{
    Disagreement, MODES, ROUTED_TERMS, battery_terms, forcer_side_value, rerouted,
    rustred_p_settings, values_agree,
};
use super::p_integrals::{VAKINT, canonical, from_dot, parse};
use super::p_library4::{
    P_LIBRARY_1_4_VARIABLE, PLibraryEngine, PLibraryRoute, engine_with_p_library, p_library_1_4,
    p_library_scope,
};

/// The eleven four-loop families (tops) of Forcer's basis.
pub const FOUR_LOOP_FAMILIES: [&str; 11] = [
    "P4L_LA4", "P4L_J1", "P4L_J2", "P4L_J4", "P4L_J5", "P4L_J6", "P4L_J11", "P4L_HAHA", "P4L_NO1",
    "P4L_NO2", "P4L_NO3",
];

/// A new engine on M3's scope 1-4 package, reached by `route` (a fresh
/// in-memory step cache per engine, so that a part of the grid never holds
/// the steps of another); `None`, a skip of `what` reported on stderr, when
/// `VAKINT_TEST_P_LIBRARY_1_4` is unset.
pub fn p_library4_engine(route: PLibraryRoute, what: &str) -> Option<PLibraryEngine> {
    let Some(directory) = p_library_1_4() else {
        eprintln!("{P_LIBRARY_1_4_VARIABLE} is unset: skipped {what}");
        return None;
    };
    assert_eq!(
        p_library_scope(&directory),
        [1, 4],
        "{P_LIBRARY_1_4_VARIABLE}={}",
        directory.display()
    );
    Some(engine_with_p_library(&directory, route))
}

/// Worker cores of a four-loop evaluation: 4, or fewer when this process may
/// use fewer. The full profile runs vakint tests one at a time.
pub fn four_loop_cores() -> usize {
    std::thread::available_parallelism().map_or(1, |cores| cores.get().min(4))
}

/// RustRed alone in `mode`, G-scheme, `terms` Laurent terms on `cores`
/// workers, under the default Ordered arithmetic or `arithmetic` (vakint's
/// text form, e.g. `exact`).
pub fn settings4(
    mode: RustRedMasterMode,
    terms: i64,
    arithmetic: Option<&str>,
    cores: usize,
) -> VakintSettings {
    let settings = VakintSettings {
        n_cores: cores,
        ..rustred_p_settings(mode, LoopNormalizationFactor::GScheme, terms)
    };
    match arithmetic {
        None => settings,
        Some(arithmetic) => VakintSettings {
            ordered_weight_arithmetic: arithmetic
                .parse()
                .unwrap_or_else(|error| panic!("arithmetic {arithmetic}: {error}")),
            ..settings
        },
    }
}

/// Why RustRed's side of an agreement has no values.
#[derive(Clone, Debug)]
pub enum RustRedSide {
    /// The wall-time budget ran out (`rustred_resources.max_elapsed`).
    OverBudget,
    /// An evaluation failed within the budget.
    Failed(String),
}

/// RustRed's value of `input` on `engine` in every mode of [`MODES`], with
/// the settings `rustred(mode)` and, given a `budget`, within it in total:
/// each evaluation runs with the remaining time as its
/// `rustred_resources.max_elapsed`, and a failure once the budget is spent is
/// [`RustRedSide::OverBudget`]. Returns the values and RustRed's seconds.
pub fn rustred_values(
    engine: &Vakint,
    input: &Atom,
    rustred: impl Fn(RustRedMasterMode) -> VakintSettings,
    budget: Option<Duration>,
) -> (Result<Vec<Atom>, RustRedSide>, f64) {
    let started = Instant::now();
    let mut values = Vec::new();
    for mode in MODES {
        let mut settings = rustred(mode);
        if let Some(budget) = budget {
            let Some(left) = budget.checked_sub(started.elapsed()) else {
                return (
                    Err(RustRedSide::OverBudget),
                    started.elapsed().as_secs_f64(),
                );
            };
            settings.rustred_resources.max_elapsed = Some(left);
        }
        match engine.evaluate(&settings, input.as_view()) {
            Ok(value) => values.push(value),
            Err(error) => {
                let seconds = started.elapsed().as_secs_f64();
                let spent = budget.is_some_and(|budget| started.elapsed() >= budget);
                let limited = matches!(error, VakintError::ResourceLimit(_));
                let side = if budget.is_some() && (spent || limited) {
                    RustRedSide::OverBudget
                } else {
                    RustRedSide::Failed(format!(
                        "RustRed ({mode}, {}): {error}",
                        settings.ordered_weight_arithmetic
                    ))
                };
                return (Err(side), seconds);
            }
        }
    }
    (Ok(values), started.elapsed().as_secs_f64())
}

/// RustRed's `values` (in the order of [`MODES`]) against the Forcer side
/// `forcer(mode)`: every mode agrees, and Forcer's value is nonzero (an
/// agreement of zeros tests nothing).
pub fn compare_with_forcer(
    values: &[Atom],
    forcer: impl Fn(RustRedMasterMode) -> Result<Atom, String>,
) -> Result<(), Disagreement> {
    for (mode, actual) in MODES.into_iter().zip(values) {
        let expected = forcer(mode).map_err(Disagreement::Failed)?;
        if expected.is_zero() {
            return Err(Disagreement::Failed(format!(
                "Forcer gives 0 ({mode}): the agreement tests nothing"
            )));
        }
        values_agree(actual, &expected, mode).map_err(|disagreement| match disagreement {
            Disagreement::Failed(error) => Disagreement::Failed(format!("{mode}: {error}")),
            Disagreement::Differs(difference) => Disagreement::Differs(format!(
                "{mode}: RustRed {actual}\nForcer {expected}\n{difference}"
            )),
        })?;
    }
    Ok(())
}

/// The agreement of `input` on `engine` (`support/p_backends.rs`
/// `agreement`, RustRed's side first): RustRed with `rustred(mode)` against
/// `forcer(mode)` in every mode. Returns RustRed's seconds.
pub fn agreement4(
    engine: &Vakint,
    input: &Atom,
    rustred: impl Fn(RustRedMasterMode) -> VakintSettings,
    forcer: impl Fn(RustRedMasterMode) -> Result<Atom, String>,
) -> Result<f64, Disagreement> {
    match rustred_values(engine, input, rustred, None) {
        (Ok(values), seconds) => compare_with_forcer(&values, forcer).map(|()| seconds),
        (Err(RustRedSide::Failed(error)), _) => Err(Disagreement::Failed(error)),
        (Err(RustRedSide::OverBudget), _) => unreachable!("no budget"),
    }
}

/// A capped agreement ([`capped_agreement`]).
#[derive(Debug)]
pub enum Capped {
    /// RustRed finished within the budget: the agreement and RustRed's
    /// seconds.
    Done(Result<(), Disagreement>, f64),
    /// RustRed did not finish within the budget (its seconds, at least the
    /// budget).
    OverBudget(f64),
}

/// [`agreement4`] with RustRed's side within `budget`; the Forcer side is
/// asked only when RustRed finished.
pub fn capped_agreement(
    engine: &Vakint,
    input: &Atom,
    rustred: impl Fn(RustRedMasterMode) -> VakintSettings,
    forcer: impl Fn(RustRedMasterMode) -> Result<Atom, String>,
    budget: Duration,
) -> Capped {
    match rustred_values(engine, input, rustred, Some(budget)) {
        (Ok(values), seconds) => Capped::Done(compare_with_forcer(&values, forcer), seconds),
        (Err(RustRedSide::OverBudget), seconds) => Capped::OverBudget(seconds),
        (Err(RustRedSide::Failed(error)), seconds) => {
            Capped::Done(Err(Disagreement::Failed(error)), seconds)
        }
    }
}

/// A frozen Forcer value that is an exact 0 (`p4l_grid.json` pins two): RustRed's
/// value on `engine` is an exact 0 in every mode.
pub fn assert_rustred_zero(
    engine: &Vakint,
    name: &str,
    input: &Atom,
    rustred: impl Fn(RustRedMasterMode) -> VakintSettings,
) {
    let values = rustred_values(engine, input, rustred, None)
        .0
        .unwrap_or_else(|side| panic!("{name}: {side:?}"));
    for (mode, value) in MODES.into_iter().zip(values) {
        assert!(
            value.expand().is_zero(),
            "{name} ({mode}): RustRed {value}, Forcer 0"
        );
    }
}

/// The Laurent depth of a live comparison of `input` (`loops` loops): the
/// depth of the battery, [`battery_terms`], when Forcer's `Symbolic` value
/// finalizes there; otherwise through `ε^(n-1)`, `n` the lowest order an
/// unknown master tail affects (`MissingMasterOrders`; as
/// `ForcerCase::deepest` does for frozen values). The Forcer run is
/// memoised per process, so the probe only finalizes again.
pub fn live_terms(form: &str, input: &Atom, loops: usize) -> i64 {
    let probe = |terms| forcer_side_value(&VAKINT, form, input, RustRedMasterMode::Symbolic, terms);
    let most = battery_terms(loops);
    match probe(most) {
        Ok(_) => most,
        Err(VakintError::MissingMasterOrders(missing)) => {
            let first = (missing.masters.iter())
                .map(|master| master.first_affected_result_order)
                .min()
                .unwrap_or_else(|| panic!("{input}: {missing}"));
            let terms = first + loops as i64;
            if let Err(error) = probe(terms) {
                panic!("{input} at {terms} terms: {error}");
            }
            terms
        }
        Err(error) => panic!("Forcer ({form}) on {input}: {error}"),
    }
}

/// The modes and depth of the invariance checks: the exact ones, through
/// `ε^0` at four loops (`support/p_backends.rs` [`ROUTED_TERMS`]).
pub const INVARIANCE_MODES: [RustRedMasterMode; 2] =
    [RustRedMasterMode::Retain, RustRedMasterMode::Symbolic];

/// RustRed's values of `input` on `engine` in [`INVARIANCE_MODES`].
pub fn invariance_values(engine: &Vakint, input: &Atom) -> Vec<Atom> {
    INVARIANCE_MODES
        .iter()
        .map(|&mode| {
            engine
                .evaluate(
                    &settings4(mode, ROUTED_TERMS, None, four_loop_cores()),
                    input.as_view(),
                )
                .unwrap_or_else(|error| panic!("{input} ({mode}): {error}"))
        })
        .collect()
}

/// `found` (the values of a variant of an integral) equals `expected` in
/// every mode of [`INVARIANCE_MODES`].
pub fn assert_same_values(found: &[Atom], expected: &[Atom], what: &str) {
    for ((mode, found), expected) in INVARIANCE_MODES.iter().zip(found).zip(expected) {
        values_agree(found, expected, *mode)
            .unwrap_or_else(|difference| panic!("{what} ({mode}): {difference}"));
    }
}

/// `input` with its external momentum reversed, `p(1) -> -p(1)` in the
/// numerator and every line (`q -> -q`): the same integral (substitute
/// `k -> -k`: every scalar product is even under both reversals).
pub fn reversed_head(input: &Atom) -> Atom {
    let p = S.p.call((1,));
    input.replace_multiple(&[Replacement::new(p.to_pattern(), (-p.clone()).to_pattern())])
}

/// Routing invariance of RustRed at four loops (M3 G7): `input` (a long form
/// `numerator*topo(prop(...))` with the external head `p(1)`) evaluates on
/// `engine` to the same nonzero values ([`INVARIANCE_MODES`]) as its
/// reversal `q -> -q` ([`reversed_head`]) and as its random unimodular
/// changes of loop momenta with line reversals (`support/p_backends.rs`
/// [`rerouted`], one per seed).
pub fn assert_routing_invariant(engine: &Vakint, input: &Atom, seeds: std::ops::Range<u64>) {
    let expected = invariance_values(engine, input);
    for (mode, value) in INVARIANCE_MODES.iter().zip(&expected) {
        assert!(!value.is_zero(), "{input} ({mode}) vanishes");
    }
    let reversed = reversed_head(input);
    assert_ne!(&reversed, input, "{input} has no head p(1)");
    assert_same_values(
        &invariance_values(engine, &reversed),
        &expected,
        &format!("q -> -q: {reversed}"),
    );
    for seed in seeds {
        let routed = rerouted(input, seed);
        assert_ne!(&routed, input, "seed {seed} does not reroute");
        assert_same_values(
            &invariance_values(engine, &routed),
            &expected,
            &format!("seed {seed}: {routed}"),
        );
    }
}

/// The DOT round trip at four loops (M3 G7): `input` (a numerator-free
/// p-integral) exported to FeynKit DOT and read back (`from_diagram`) is the
/// same integral (canonical form) with the same nonzero values on `engine`.
pub fn assert_dot_round_trip_value(engine: &Vakint, input: &Atom) {
    let dot = vakint::p_integral::to_dot(input.as_view(), "phi")
        .unwrap_or_else(|error| panic!("{input}: {error}"));
    let imported = from_dot(&dot);
    assert_eq!(canonical(&imported), canonical(input), "{input}\n{dot}");
    let expected = invariance_values(engine, input);
    for (mode, value) in INVARIANCE_MODES.iter().zip(&expected) {
        assert!(!value.is_zero(), "{input} ({mode}) vanishes");
    }
    assert_same_values(
        &invariance_values(engine, &imported),
        &expected,
        &format!("{input} through DOT\n{dot}"),
    );
}

/// The top of `family` with line powers `powers` in the canonical long form
/// (family routing, `k(1..=4)`, head `p(1)`).
pub fn long_form_top(family: &str, powers: &[i64]) -> Atom {
    let powers = powers.iter().map(i64::to_string).collect::<Vec<_>>();
    canonical(&parse(&format!(
        "topo({family}(p(1),{}))",
        powers.join(",")
    )))
}
