//! The RustRed p backend through vakint's public API (M2-V3), with an injected
//! RustRed p library: the frozen one-loop Forcer battery and the dispatch next
//! to Forcer. The exact tag-basis battery, the Wick-sign, channel, cache and
//! arithmetic checks live with the backend (`rustred_evaluation::p_integral`).

#![cfg(feature = "feynkit-ingress")]

#[path = "support/p_integrals.rs"]
mod p_integrals;
#[path = "support/p_rustred.rs"]
mod p_rustred;

use p_integrals::{forcer_evaluate, forcer_form, forcer_settings, parse};
use p_rustred::{
    assert_rustred_p_symbolic_matches_frozen, reference_p_engine, rustred_p_evaluate,
    rustred_p_settings,
};
use symbolica::atom::AtomCore;
use vakint::p_integral::forcer_fixtures::ForcerFixture;
use vakint::{EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, RustRedMasterMode};

fn one_loop_battery() -> ForcerFixture {
    ForcerFixture::parse(include_str!("fixtures/forcer/p1l_battery.json")).unwrap()
}

/// Every case of the one-loop battery, evaluated by RustRed with the
/// injected reference library, is Forcer's frozen value.
#[test]
fn the_injected_reference_library_reproduces_the_one_loop_battery() {
    let engine = reference_p_engine();
    let settings = rustred_p_settings(RustRedMasterMode::Symbolic, 4);
    let battery = one_loop_battery();
    assert_eq!(battery.cases.len(), 45);
    for case in &battery.cases {
        assert_rustred_p_symbolic_matches_frozen(&engine, &settings, case);
    }
}

/// The default evaluation order with a live FORM 5: with the reference
/// library the bubbles of an expression go to RustRed and its kite (outside
/// the library's scope) to Forcer, and the whole expression, vacuum term
/// included, equals Forcer's evaluation. Skipped without FORM 5.
#[test]
fn the_default_order_splits_p_integrals_between_rustred_and_forcer() {
    let Some(form) = forcer_form() else {
        eprintln!("FORM 5 is unavailable: skipped the mixed RustRed/Forcer evaluation");
        return;
    };
    let forcer = forcer_settings(
        &form,
        RustRedMasterMode::Symbolic,
        LoopNormalizationFactor::GScheme,
        4,
    );
    let mut order = EvaluationOrder::available(&form, "python3");
    assert!(matches!(order.0[0], EvaluationMethod::RustRed(_)));
    for method in order.0.iter_mut() {
        match method {
            EvaluationMethod::Forcer(options) => options.master_mode = RustRedMasterMode::Symbolic,
            EvaluationMethod::RustRed(options) => options.master_mode = RustRedMasterMode::Symbolic,
            _ => {}
        }
    }
    let default = vakint::VakintSettings {
        evaluation_order: order,
        ..forcer.clone()
    };
    let bubble = parse("dot(k(1),p(2))^3*topo(P1L(-p(2),2,1)) + topo(P1L(p(1),1,1))");
    let kite = parse("topo(P2L(p(1),1,1,1,1,1))");
    let vacuum = parse("topo(prop(1,edge(1,1),k(1),muvsq,2))");
    let engine = reference_p_engine();
    let rustred_only = rustred_p_settings(RustRedMasterMode::Symbolic, 4);
    let bubble_value = rustred_p_evaluate(&engine, &rustred_only, &bubble);
    assert!(
        (&bubble_value - forcer_evaluate(&bubble, &forcer))
            .expand()
            .is_zero()
    );
    let sum = &bubble + &kite + &vacuum;
    let mixed = rustred_p_evaluate(&engine, &default, &sum);
    let expected = forcer_evaluate(&(&bubble + &kite), &forcer)
        + rustred_p_evaluate(&engine, &default, &vacuum);
    assert!(
        (&mixed - &expected).expand().is_zero(),
        "{mixed}\n{expected}"
    );
}
