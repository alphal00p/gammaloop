//! Shared helpers of the RustRed p backend tests (M2-V3): engines that reduce
//! p-integrals with an injected RustRed p library, RustRed-only settings and
//! the comparison with the frozen Forcer fixtures.
//!
//! Include it with `#[path = "support/p_rustred.rs"] mod p_rustred;` (it needs
//! no other support file). Until vakint ships its built-in 1L-3L library
//! (M2-V4), the only library is RustRed's one-loop reference program
//! `rustred::propagator::testing::one_loop_reference_library()`: the massless
//! bubble as an exact closed form behind the real engine.

#![allow(dead_code)]

use std::sync::Arc;

use rustred::propagator::library::PLibrary;
use symbolica::atom::{Atom, AtomCore};
use vakint::p_integral::forcer_fixtures::ForcerCase;
use vakint::utils::simplify_real;
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedEvaluationOptions,
    RustRedMasterMode, Vakint, VakintSettings,
};

/// An engine whose p-integrals within the scope of `library` go to RustRed.
pub fn rustred_p_engine(library: Arc<PLibrary>) -> Vakint {
    let mut engine = Vakint::new().unwrap();
    engine.inject_p_library(library).unwrap();
    engine
}

/// [`rustred_p_engine`] with RustRed's one-loop reference library (scope 1).
pub fn reference_p_engine() -> Vakint {
    rustred_p_engine(rustred::propagator::testing::one_loop_reference_library().unwrap())
}

/// RustRed alone in `mode` (no FORM is ever called), `normalization`, through
/// `terms` Laurent terms, dot notation, one core, 40 digits, the default
/// Ordered arithmetic. The one RustRed p settings helper of the vakint tests
/// (M2-V5's `tests/support/p_backends.rs` re-exports it).
pub fn rustred_p_settings(
    mode: RustRedMasterMode,
    normalization: LoopNormalizationFactor,
    terms: i64,
) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions { master_mode: mode },
        )]),
        form_exe_path: "/rustred-p-integrals-must-not-call-form".into(),
        integral_normalization_factor: normalization,
        number_of_terms_in_epsilon_expansion: terms,
        run_time_decimal_precision: 40,
        use_dot_product_notation: true,
        n_cores: 1,
        ..VakintSettings::default()
    }
}

/// The evaluation of `input` with `settings`.
pub fn rustred_p_evaluate(engine: &Vakint, settings: &VakintSettings, input: &Atom) -> Atom {
    engine
        .evaluate(settings, input.as_view())
        .unwrap_or_else(|error| panic!("{input}: {error}"))
}

/// The case evaluated by `engine` (`settings` in `Symbolic` mode) equals its
/// frozen Forcer output finalized like a Forcer evaluation, exactly.
pub fn assert_rustred_p_symbolic_matches_frozen(
    engine: &Vakint,
    settings: &VakintSettings,
    case: &ForcerCase,
) {
    let actual = rustred_p_evaluate(engine, settings, &case.expression().unwrap());
    let expected = simplify_real(
        case.finalized(settings, RustRedMasterMode::Symbolic)
            .unwrap_or_else(|error| panic!("{}: {error}", case.name))
            .as_view(),
    );
    assert!(
        (&actual - &expected).expand().is_zero(),
        "{}: RustRed {actual} vs Forcer {expected}",
        case.name
    );
}
