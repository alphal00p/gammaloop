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

/// The two-loop sunset (three lines from the in to the out vertex).
const SUNSET: &str = "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(1,2),k(2),0,1)\
    *prop(3,edge(1,2),p(1)-k(1)-k(2),0,1))";

/// The kite from vertex 1 to vertex 5 followed by the sunset from 5 to 2.
const KITE_THEN_SUNSET: &str = "topo(prop(1,edge(1,3),k(1),0,1)*prop(2,edge(1,4),p(1)-k(1),0,1)\
    *prop(3,edge(3,4),k(2),0,1)*prop(4,edge(3,5),k(1)-k(2),0,1)*prop(5,edge(4,5),p(1)-k(1)+k(2),0,1)\
    *prop(6,edge(5,2),k(3),0,1)*prop(7,edge(5,2),k(4),0,1)*prop(8,edge(5,2),p(1)-k(3)-k(4),0,1))";

#[test]
fn a_kite_and_a_sunset_in_series_multiply() {
    // The blocks meet with five momenta at vertex 5: no Forcer topology has
    // such a vertex, so Forcer must be given the cubic family graph.
    let Some(form) = forcer_form() else {
        eprintln!("FORM 5 is unavailable: skipped the live kite-sunset evaluation");
        return;
    };
    let settings = forcer_settings(
        &form,
        RustRedMasterMode::Symbolic,
        LoopNormalizationFactor::GScheme,
        6,
    );
    let glued = forcer_evaluate(&parse(KITE_THEN_SUNSET), &settings);
    let blocks =
        forcer_evaluate(&parse(KITE), &settings) * forcer_evaluate(&parse(SUNSET), &settings);
    assert_eq!(
        at_unit_scale(&glued, -4, 0),
        at_unit_scale(&blocks.expand(), -4, 0)
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

/// Analytic identities between p-integrals with odd-rank numerators and with
/// the external momentum reversed (σ = -1), one to four loops (M1 review,
/// missing test 1): `k·q = (k² + q² - (k-q)²)/2` and its partial fractions.
mod identities {
    use super::*;

    /// The bubble with `q` written as `head` (`k - head` on line 2).
    fn bubble(head: &str) -> String {
        format!("topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-({head}),0,1))")
    }

    /// The kite of the interface example with line powers `a` and `q = head`.
    fn kite(a: [i64; 5], head: &str) -> String {
        format!(
            "topo(prop(1,edge(1,3),k(1),0,{})*prop(2,edge(1,4),{head}-k(1),0,{})\
             *prop(3,edge(3,4),k(2),0,{})*prop(4,edge(3,2),k(1)-k(2),0,{})\
             *prop(5,edge(4,2),{head}-k(1)+k(2),0,{}))",
            a[0], a[1], a[2], a[3], a[4]
        )
    }

    /// `family(±p(1), 1, ..., 1)` with the listed lines set to other powers.
    fn top(family: &str, head: &str, lines: usize, changes: &[(usize, i64)]) -> String {
        let mut powers = vec![1; lines];
        for &(line, power) in changes {
            powers[line - 1] = power;
        }
        let powers = powers.iter().map(i64::to_string).collect::<Vec<_>>();
        format!("topo({family}({head},{}))", powers.join(","))
    }

    #[test]
    fn one_loop_odd_and_even_rank() {
        let b = bubble("p(1)");
        assert_p_identity(
            &format!("dot(k(1),p(1))*{b}"),
            &format!("dot(p(1),p(1))/2*{b}"),
        );
        // (k·q)² = (k·q) q²/2 up to scaleless integrals.
        assert_p_identity(
            &format!("dot(k(1),p(1))^2*{b}"),
            &format!("dot(p(1),p(1))/2*dot(k(1),p(1))*{b}"),
        );
        assert_p_identity_with(
            &format!("dot(k(1),p(1))*{b}"),
            &format!("dot(p(1),p(1))/2*{b}"),
            LoopNormalizationFactor::MSbar,
        );
    }

    #[test]
    fn one_loop_with_the_external_momentum_reversed() {
        for head in ["-p(1)", "-p(3)"] {
            let b = bubble(head);
            let p = &head[1..];
            assert_p_identity(
                &format!("dot(k(1),{p})*{b}"),
                &format!("-dot({p},{p})/2*{b}"),
            );
        }
    }

    #[test]
    fn two_loop_partial_fractions() {
        let a = [1, 2, 1, 3, 1];
        let k = |a, head| kite(a, head);
        assert_p_identity(
            &format!("dot(k(1),p(1))*{}", k(a, "p(1)")),
            &format!(
                "1/2*{}+1/2*dot(p(1),p(1))*{}-1/2*{}",
                k([0, 2, 1, 3, 1], "p(1)"),
                k(a, "p(1)"),
                k([1, 1, 1, 3, 1], "p(1)")
            ),
        );
        assert_p_identity(
            &format!("dot(k(1),k(2))*{}", k(a, "p(1)")),
            &format!(
                "1/2*{}+1/2*{}-1/2*{}",
                k([0, 2, 1, 3, 1], "p(1)"),
                k([1, 2, 0, 3, 1], "p(1)"),
                k([1, 2, 1, 2, 1], "p(1)")
            ),
        );
        // Line 5 is q - k1 + k2: k2·q = [(q-k1+k2)² - (q-k1)² - k2²]/2 + k1·k2.
        let a5 = [1, 2, 1, 3, 2];
        assert_p_identity(
            &format!("dot(k(2),p(1))*{}", k(a5, "p(1)")),
            &format!(
                "1/2*{}-1/2*{}-1/2*{}+dot(k(1),k(2))*{}",
                k([1, 2, 1, 3, 1], "p(1)"),
                k([1, 1, 1, 3, 2], "p(1)"),
                k([1, 2, 0, 3, 2], "p(1)"),
                k(a5, "p(1)")
            ),
        );
    }

    #[test]
    fn two_loop_with_the_external_momentum_reversed() {
        let a = [1, 2, 1, 3, 1];
        assert_p_identity(
            &format!("dot(k(1),p(2))*{}", kite(a, "-p(2)")),
            &format!(
                "1/2*{}-1/2*dot(p(2),p(2))*{}-1/2*{}",
                kite([1, 1, 1, 3, 1], "-p(2)"),
                kite(a, "-p(2)"),
                kite([0, 2, 1, 3, 1], "-p(2)")
            ),
        );
    }

    #[test]
    fn three_and_four_loops_with_the_external_momentum_reversed() {
        let no = |changes: &[(usize, i64)]| top("P3L_NO", "-p(1)", 8, changes);
        assert_p_identity(
            &format!("dot(k(1),p(1))*{}", no(&[])),
            &format!(
                "1/2*{}+1/2*dot(p(1),p(1))*{}-1/2*{}",
                no(&[(1, 0)]),
                no(&[]),
                no(&[(5, 0)])
            ),
        );
        let la4 = |changes: &[(usize, i64)]| top("P4L_LA4", "-p(1)", 11, changes);
        assert_p_identity(
            &format!(
                "(dot(k(1),p(1))+dot(k(2),p(1))+dot(k(3),p(1))+dot(k(4),p(1)))*{}",
                la4(&[])
            ),
            &format!(
                "1/2*{}+1/2*dot(p(1),p(1))*{}-1/2*{}",
                la4(&[(9, 0)]),
                la4(&[]),
                la4(&[(11, 0)])
            ),
        );
    }

    #[test]
    fn index_notation_and_linear_dot_arguments_evaluate_alike() {
        let b = bubble("p(1)");
        assert_p_identity(
            &format!("k(1,mu)*p(1,mu)*{b}"),
            &format!("dot(k(1),p(1))*{b}"),
        );
        assert_p_identity(
            &format!("dot(k(1)+p(1),p(1))*{b}"),
            &format!("(dot(k(1),p(1))+dot(p(1),p(1)))*{b}"),
        );
        let k = kite([1, 2, 1, 3, 1], "p(1)");
        assert_p_identity(
            &format!("dot(k(1)+k(2),p(1))*{k}"),
            &format!("(dot(k(1),p(1))+dot(k(2),p(1)))*{k}"),
        );
    }
}

/// Evaluated normalizations of section 1, items 5 and 6 (M1 review, missing
/// tests 2 and 3).
mod normalization {
    use super::*;

    #[test]
    fn tree_lines_evaluate_as_numerator_factors() {
        for power in [1, 2] {
            assert_p_identity(
                &format!(
                    "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1)\
                     *prop(3,edge(3,1),p(1),0,{power}))"
                ),
                &format!("dot(p(1),p(1))^-{power}*{BUBBLE}"),
            );
        }
    }

    #[test]
    fn negative_short_form_powers_evaluate_as_numerators() {
        // Line 5 of P2L is k1 + k2 + p.
        assert_p_identity(
            "topo(P2L(p(1),1,1,1,1,-1))",
            "(dot(k(1),k(1))+dot(k(2),k(2))+dot(p(1),p(1))+2*dot(k(1),k(2))+2*dot(k(1),p(1))\
             +2*dot(k(2),p(1)))*topo(P2L(p(1),1,1,1,1,0))",
        );
    }

    #[test]
    fn negative_long_form_powers_evaluate_as_numerators() {
        // Line 3 of this kite carries k(2).
        let kite = |third: i64| {
            format!(
                "topo(prop(1,edge(1,3),k(1),0,1)*prop(2,edge(1,4),p(1)-k(1),0,1)\
                 *prop(3,edge(3,4),k(2),0,{third})*prop(4,edge(3,2),k(1)-k(2),0,1)\
                 *prop(5,edge(4,2),p(1)-k(1)+k(2),0,1))"
            )
        };
        for power in [1, 2] {
            assert_p_identity(
                &kite(-power),
                &format!("dot(k(2),k(2))^{power}*{}", kite(0)),
            );
        }
    }

    #[test]
    fn a_reversed_short_form_head_evaluates_alike() {
        // A scalar integral is even in q.
        assert_p_identity("topo(P2L(-p(1),1,2,1,3,1))", "topo(P2L(p(1),1,2,1,3,1))");
    }

    #[test]
    fn a_bubble_chain_through_a_cut_vertex_is_not_scaleless() {
        // Two bubbles in series through vertex 3: 1/ε² in the G-scheme.
        assert_forcer_g_scheme(
            "topo(prop(1,edge(1,3),k(1),0,1)*prop(2,edge(3,1),k(1)-p(1),0,1)\
             *prop(3,edge(3,2),k(2),0,1)*prop(4,edge(2,3),k(2)-p(1),0,1))",
            2,
            1,
            &[Atom::num(1), Atom::Zero, Atom::Zero, Atom::Zero],
        );
    }
}

/// The `Oep` guard through `Vakint::evaluate` at four loops (M1 review,
/// missing test 5): asking for more orders than the masters provide is a
/// `MissingMasterOrders` error, one order fewer passes. The boundary is
/// pinned (10 pass, 11 fail) so that a guard firing too early fails too.
#[test]
fn the_four_loop_oep_guard_fails_closed_through_evaluate() {
    let Some(form) = forcer_form() else {
        eprintln!("FORM 5 is unavailable: skipped the four-loop Oep guard");
        return;
    };
    let top = parse("topo(P4L_NO1(p(1),1,1,1,1,1,1,1,1,1,1,1))");
    for (mode, passing, failing) in [
        (RustRedMasterMode::Symbolic, 10, 11),
        (RustRedMasterMode::Numerical, 10, 11),
    ] {
        let evaluate = |terms| {
            vakint::Vakint::new().unwrap().evaluate(
                &forcer_settings(&form, mode, LoopNormalizationFactor::GScheme, terms),
                top.as_view(),
            )
        };
        let value = evaluate(passing).unwrap_or_else(|error| panic!("{mode:?}: {error}"));
        assert!(!value.to_string().contains("Oep"), "{mode:?}: {value}");
        match evaluate(failing) {
            Err(VakintError::MissingMasterOrders(_)) => {}
            other => panic!(
                "{mode:?} through {failing} terms: expected MissingMasterOrders, got {other:?}"
            ),
        }
    }
}

/// The G-scheme on vacuum integrals (M1 review, missing test 6).
mod g_scheme {
    use super::*;
    use vakint::EvaluationMethod;

    const TADPOLE: &str = "topo(prop(1,edge(1,1),k(1),muvsq,2))";

    /// ζ(2) to 40 digits.
    const ZETA2: &str = "1.644934066848226436472415166646025189219`40";

    /// The Laurent coefficients `ε^-1, ε^0, ε^1` of `value` at unit scales
    /// agree with `expected` to 25 digits.
    fn assert_unit_scale_values(value: &Atom, expected: [&str; 3]) {
        for (order, (found, expected)) in
            at_unit_scale(value, -1, 1).iter().zip(expected).enumerate()
        {
            let difference = (found - parse(expected)).expand();
            let difference: f64 = difference
                .to_string()
                .split('`')
                .next()
                .unwrap()
                .parse()
                .unwrap_or_else(|_| {
                    panic!(
                        "ε^{}: {difference} is not a number ({value})",
                        order as i64 - 1
                    )
                });
            assert!(
                difference.abs() < 1e-25,
                "ε^{}: {found} != {expected}",
                order as i64 - 1
            );
        }
    }

    #[test]
    fn the_massive_tadpole_without_form() {
        // Γ(ε)·Γ(2-2ε)/(Γ(1+ε)·Γ(1-ε)²)·(μ²/m²)^ε = 1/ε - 2 + ζ2·ε + ... at μ = m.
        let settings = VakintSettings {
            evaluation_order: EvaluationOrder::rustred_only(),
            integral_normalization_factor: LoopNormalizationFactor::GScheme,
            number_of_terms_in_epsilon_expansion: 3,
            form_exe_path: "/vakint-forcer-tests/missing-form".into(),
            ..VakintSettings::default()
        };
        let value = vakint::Vakint::new()
            .unwrap()
            .evaluate(&settings, parse(TADPOLE).as_view())
            .unwrap();
        assert_unit_scale_values(&value, ["1", "-2", ZETA2]);
    }

    #[test]
    fn a_vacuum_integral_plus_a_p_integral() {
        let Some(form) = forcer_form() else {
            eprintln!("FORM 5 is unavailable: skipped the mixed G-scheme evaluation");
            return;
        };
        let mut order = EvaluationOrder::available(&form, "python3");
        for method in order.0.iter_mut() {
            if let EvaluationMethod::Forcer(options) = method {
                options.master_mode = RustRedMasterMode::Numerical;
            }
        }
        let settings = VakintSettings {
            evaluation_order: order,
            run_time_decimal_precision: 32,
            ..forcer_settings(
                &form,
                RustRedMasterMode::Numerical,
                LoopNormalizationFactor::GScheme,
                3,
            )
        };
        // The tadpole above plus the bubble, 1/ε exactly in the G-scheme.
        let value = forcer_evaluate(&parse(&format!("{TADPOLE}+{BUBBLE}")), &settings);
        assert_unit_scale_values(&value, ["2", "-2", ZETA2]);
    }
}

/// A host whose FORM is older than 5 (M1 review F3): Forcer is optional in a
/// mixed evaluation order, never runs, and a p-integral without another
/// method reports `ForcerUnavailable`. No live FORM is needed.
#[test]
fn a_form_4_host_never_runs_forcer() {
    let form = fake_form("4.3.1");
    let vakint = vakint::Vakint::new().unwrap();
    let with_order = |evaluation_order| VakintSettings {
        evaluation_order,
        form_exe_path: form.path.clone(),
        ..VakintSettings::default()
    };
    // EvaluationOrder::all() no longer needs FORM 5; pySecDec may be missing.
    match vakint.validate_settings(&with_order(EvaluationOrder::all())) {
        Ok(()) | Err(VakintError::PySecDecUnavailable | VakintError::PySecDecVersion(_)) => {}
        Err(error) => panic!("EvaluationOrder::all() with FORM 4.3.1: {error}"),
    }
    let mut order = EvaluationOrder::all();
    order
        .0
        .retain(|method| !matches!(method, vakint::EvaluationMethod::PySecDec(_)));
    let settings = with_order(order);
    vakint.validate_settings(&settings).unwrap();
    // A vacuum integral evaluates (RustRed) ...
    let tadpole = vakint
        .evaluate(
            &settings,
            parse("topo(prop(1,edge(1,1),k(1),muvsq,2))").as_view(),
        )
        .unwrap();
    assert!(!tadpole.is_zero());
    // ... and a p-integral finds no method, because Forcer needs FORM 5.
    match vakint.evaluate(&settings, parse(KITE).as_view()) {
        Err(VakintError::ForcerUnavailable { reason, .. }) => {
            assert!(reason.contains("4.3.1"), "{reason}")
        }
        other => panic!("expected ForcerUnavailable, got {other:?}"),
    }
    // An order of Forcer alone still needs FORM 5 up front.
    match vakint.validate_settings(&with_order(EvaluationOrder::forcer_only())) {
        Err(VakintError::ForcerUnavailable { form_exe_path, .. }) => {
            assert_eq!(form_exe_path, form.path)
        }
        other => panic!("expected ForcerUnavailable, got {other:?}"),
    }
}
