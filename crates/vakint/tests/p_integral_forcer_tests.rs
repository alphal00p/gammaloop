//! The Forcer backend end to end (`docs/p_integrals.md`, sections 4, 5 and 7):
//! p-integrals written in any routing, as atoms, short forms or DOT graphs,
//! evaluated by vakint with the vendored Forcer.
//!
//! Live runs need FORM 5 (the dev shell and CI provide it); without it these
//! tests report a skip. The frozen batteries of `src/p_integral/forcer/tests.rs`
//! test the input generation, the parser and the finalizer without FORM.

#![cfg(feature = "feynkit-ingress")]

#[path = "support/p_integrals.rs"]
mod p_integrals;

use p_integrals::*;
use symbolica::atom::{Atom, AtomCore};
use vakint::{
    EvaluationOrder, LoopNormalizationFactor, RustRedMasterMode, VakintError, VakintSettings,
};

const BUBBLE: &str = "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))";

/// The two-loop kite of the interface example, in the user's routing.
const KITE: &str = "topo(prop(1,edge(1,3),k(1),0,1)*prop(2,edge(1,4),p(1)-k(1),0,1)\
    *prop(3,edge(3,4),k(2),0,1)*prop(4,edge(3,2),k(1)-k(2),0,1)*prop(5,edge(4,2),p(1)-k(1)+k(2),0,1))";

/// The kite with powers that break all of its symmetries.
const ASYMMETRIC_KITE: &str = "topo(prop(1,edge(1,3),k(1),0,1)*prop(2,edge(1,4),p(1)-k(1),0,2)\
    *prop(3,edge(3,4),k(2),0,1)*prop(4,edge(3,2),k(1)-k(2),0,3)*prop(5,edge(4,2),p(1)-k(1)+k(2),0,1))";

/// The compact FeynKit DOT of `docs/p_integrals.md` (section 3) and the same
/// graph as an atom: the canonical two-loop family with a dot on line 5.
const DOT_KITE: &str = "digraph p2l_kite {
  exte0 [style=invis]; exte1 [style=invis];
  exte0 -> v1 [id=0 particle=\"phi\"];
  v2 -> exte1 [id=1 particle=\"phi\"];
  v2 -> v4 [id=2 particle=\"phi\" lmb_id=0];
  v3 -> v4 [id=3 particle=\"phi\" lmb_id=1];
  v4 -> v1 [id=4 particle=\"phi\"];
  v3 -> v2 [id=5 particle=\"phi\"];
  v1 -> v3 [id=6 particle=\"phi\" power=2];
}";
const ATOM_KITE: &str = "topo(prop(2,edge(2,4),k(1),0,1)*prop(3,edge(3,4),k(2),0,1)\
    *prop(4,edge(4,1),k(1)+k(2),0,1)*prop(5,edge(3,2),k(1)+p(1),0,1)\
    *prop(6,edge(1,3),k(1)+k(2)+p(1),0,2))";

#[test]
fn the_bubble_is_one_over_epsilon_in_the_g_scheme() {
    assert_forcer_g_scheme(BUBBLE, 1, 1, &[Atom::num(1), Atom::Zero, Atom::Zero]);
}

#[test]
fn the_kite_is_six_zeta3() {
    // Five Minkowski propagators: (-1)^5 times the G-scheme value 6ζ3.
    assert_forcer_g_scheme(
        KITE,
        2,
        0,
        &[Atom::Zero, Atom::Zero, -Atom::num(6) * zeta(3)],
    );
}

#[test]
fn the_three_loop_tops_are_twenty_zeta5() {
    for top in ["P3L_LA", "P3L_NO"] {
        assert_forcer_g_scheme(
            &format!("topo({top}(p(1),1,1,1,1,1,1,1,1))"),
            3,
            0,
            &[Atom::Zero, Atom::Zero, Atom::Zero, Atom::num(20) * zeta(5)],
        );
    }
}

#[test]
fn the_four_loop_ladder_is_seventy_zeta7() {
    assert_forcer_g_scheme(
        "topo(P4L_LA4(p(1),1,1,1,1,1,1,1,1,1,1,1))",
        4,
        0,
        &[
            Atom::Zero,
            Atom::Zero,
            Atom::Zero,
            Atom::Zero,
            -Atom::num(70) * zeta(7),
        ],
    );
}

#[test]
fn forcer_evaluations_are_routing_invariant() {
    assert_forcer_routing_invariant(
        &format!("dot(k(1),p(1))*dot(k(2),p(1))*{ASYMMETRIC_KITE}"),
        0..2,
    );
}

#[test]
fn user_routings_evaluate_alike_at_one_three_and_four_loops() {
    assert_forcer_routing_invariant(
        "dot(k(1),p(1))^2*topo(prop(1,edge(1,2),k(1),0,2)*prop(2,edge(2,1),k(1)-p(1),0,1))",
        0..2,
    );
    assert_forcer_atom_routing_invariant(
        &(parse("dot(k(1),p(1))*dot(k(2),k(3))")
            * canonical(&parse("topo(P3L_BE(p(1),1,2,1,1,1,1,1,1))"))),
        0..2,
    );
    assert_forcer_atom_routing_invariant(
        &(parse("dot(k(1),p(1))") * canonical(&parse("topo(P4L_LA4(p(1),1,1,1,1,1,1,1,1,1,1,1))"))),
        0..1,
    );
}

#[test]
fn three_and_four_loop_dot_graphs_evaluate_like_their_atoms() {
    for top in [
        "topo(P3L_NO(p(1),1,1,1,1,1,1,1,2))",
        "topo(P4L_NO1(p(1),1,1,1,1,1,1,1,1,1,1,1))",
    ] {
        assert_forcer_dot_like_atom(&canonical(&parse(top)));
    }
}

#[test]
fn the_default_order_hands_p_integrals_to_forcer() {
    assert_default_order_dispatches_to_forcer("topo(prop(1,edge(1,1),k(1),muvsq,1))", KITE);
}

#[test]
fn dot_graphs_evaluate_like_their_atoms() {
    let Some(form) = forcer_form() else {
        eprintln!("FORM 5 is unavailable: skipped the live Forcer DOT evaluation");
        return;
    };
    let settings = forcer_settings(
        &form,
        RustRedMasterMode::Retain,
        LoopNormalizationFactor::pySecDec,
        4,
    );
    assert_eq!(
        forcer_evaluate(&from_dot(DOT_KITE), &settings),
        forcer_evaluate(&parse(ATOM_KITE), &settings)
    );
}

#[test]
fn numerical_masters_reproduce_the_kite() {
    let Some(form) = forcer_form() else {
        eprintln!("FORM 5 is unavailable: skipped the live numerical Forcer evaluation");
        return;
    };
    let settings = VakintSettings {
        run_time_decimal_precision: 40,
        ..forcer_settings(
            &form,
            RustRedMasterMode::Numerical,
            LoopNormalizationFactor::GScheme,
            3,
        )
    };
    let finite = at_unit_scale(&forcer_evaluate(&parse(KITE), &settings), 0, 0).remove(0);
    let six_zeta3 = parse("7.21234141895756571239842896906869994458992`40");
    let relative = ((finite + &six_zeta3) / six_zeta3).expand();
    let relative: f64 = relative
        .to_string()
        .split('`')
        .next()
        .unwrap()
        .parse()
        .unwrap();
    assert!(relative.abs() < 1e-30, "{relative}");
}

#[test]
fn forcer_without_form_five_is_reported() {
    let missing = "/vakint-forcer-tests/missing-form";
    let forcer_only = forcer_settings(
        missing,
        RustRedMasterMode::Symbolic,
        LoopNormalizationFactor::pySecDec,
        4,
    );
    // The default order without FORM 5 has no method for p-integrals.
    let default_order = VakintSettings {
        evaluation_order: EvaluationOrder::available(missing, "python3"),
        ..forcer_only.clone()
    };
    for settings in [forcer_only, default_order] {
        match vakint::Vakint::new()
            .unwrap()
            .evaluate(&settings, parse(KITE).as_view())
        {
            Err(VakintError::ForcerUnavailable { form_exe_path, .. }) => {
                assert_eq!(form_exe_path, missing)
            }
            other => panic!("expected ForcerUnavailable, got {other:?}"),
        }
    }
}
