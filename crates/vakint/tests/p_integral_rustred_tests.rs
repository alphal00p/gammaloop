//! Tier B of the Forcer battery (M2-V5; M2 plan, section 8, gate G7): the
//! RustRed p backend against a live Forcer (FORM 5), through vakint's public
//! API, exactly (`Retain` in the tag basis, `Symbolic`) and numerically
//! (`Numerical`, 40 digits) (`support/p_backends.rs`). Without FORM 5 the live
//! comparisons report a skip; the RustRed identities and reroutings need no
//! FORM.
//!
//! STAGE 1 (until M2-V4 ships the built-in 1L-3L library) runs the one-loop
//! tests against RustRed's one-loop reference library; the two- and
//! three-loop tests are ignored until STAGE 2.

#![cfg(feature = "feynkit-ingress")]

#[path = "support/p_integrals.rs"]
mod p_integrals;
#[path = "support/p_rustred.rs"]
mod p_rustred;
#[path = "support/p_backends.rs"]
mod p_backends;

use p_backends::*;
use p_integrals::{canonical, forcer_form, from_dot, parse};
use symbolica::atom::{Atom, AtomCore};
use vakint::RustRedMasterMode;
use vakint::p_integral::forcer_fixtures::nonzero_sectors;

/// The scalar numerators of rank 0 to 4 of the bubble grid, each with its own
/// coefficient (distinct primes), summed into one input per pair of powers.
const BUBBLE_NUMERATORS: [&str; 8] = [
    "1",
    "2*dot(k(1),p(1))",
    "3*dot(k(1),k(1))",
    "5*dot(k(1),p(1))^2",
    "7*dot(k(1),p(1))*dot(k(1),k(1))",
    "11*dot(k(1),p(1))^3",
    "13*dot(k(1),k(1))^2",
    "17*dot(k(1),p(1))^4",
];

/// `P1L(p(1), a, b)` times every numerator of [`BUBBLE_NUMERATORS`].
fn bubble_grid_input(a: i64, b: i64) -> String {
    format!("({})*topo(P1L(p(1),{a},{b}))", BUBBLE_NUMERATORS.join("+"))
}

/// The bubble with `q` written as `head` (`k - head` on line 2), long form.
fn bubble(head: &str) -> String {
    format!("topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-({head}),0,1))")
}

#[test]
fn the_bubble_grid_agrees_with_forcer() {
    for a in 1..=3 {
        for b in 1..=3 {
            assert_p_backends_agree(&bubble_grid_input(a, b), battery_terms(1));
        }
    }
}

/// Odd ranks with the head written `-p(3)`: the external momentum reversed
/// and renamed.
#[test]
fn odd_rank_bubbles_with_a_reversed_head_agree_with_forcer() {
    for input in [
        "dot(k(1),p(3))*topo(P1L(-p(3),1,1))",
        "dot(k(1),p(3))^3*topo(P1L(-p(3),2,1))",
        "dot(k(1),p(3))*dot(k(1),k(1))*topo(P1L(-p(3),1,3))",
    ] {
        assert_p_backends_agree(input, battery_terms(1));
    }
}

/// The Wick sign `(-1)^(Σa+r)` is applied once (critique-phys SF3): for an
/// odd rank (`Σa = 2`, `r = 1`) and an odd `Σa` (`(2,1)`, `r = 0`) RustRed
/// agrees with Forcer, and not with Forcer's negative.
#[test]
fn an_odd_bubble_carries_its_wick_sign_once() {
    let Some(form) = forcer_form() else {
        eprintln!("FORM 5 is unavailable: skipped the Wick sign of the odd bubbles");
        return;
    };
    for input in [
        "dot(k(1),p(1))*topo(P1L(p(1),1,1))",
        "topo(P1L(p(1),2,1))",
    ] {
        assert_p_backends_agree(input, battery_terms(1));
        assert_wick_sign_does_not_cancel(&form, input);
    }
}

/// RustRed's `Retain` value of `input` differs from minus Forcer's.
fn assert_wick_sign_does_not_cancel(form: &str, input: &str) {
    let input = parse(input);
    let mode = RustRedMasterMode::Retain;
    let rustred = rustred_value(&input, mode, battery_terms(1)).unwrap();
    let forcer = live_forcer_value(form, &input, mode, battery_terms(1)).unwrap();
    assert!(
        values_agree(&rustred, &(-forcer), mode).is_err(),
        "{input}: the sign must not cancel"
    );
}

/// One expression with p-integrals of the heads `p(1)` and `p(3)` (separate
/// channels) and a vacuum tadpole (RustRed's vacuum backend on both sides).
#[test]
fn mixed_heads_and_a_vacuum_tadpole_agree_with_forcer() {
    assert_p_backends_agree(
        "dot(k(1),p(1))^2*topo(P1L(p(1),2,1)) + dot(k(1),p(3))*topo(P1L(-p(3),1,2)) \
         + topo(prop(1,edge(1,1),k(1),muvsq,2))",
        battery_terms(1),
    );
}

/// DOT ingress: a dotted bubble exported to FeynKit DOT and read back
/// (`from_diagram`) evaluates alike with both backends.
#[test]
fn a_dot_bubble_agrees_with_forcer() {
    let atom = canonical(&parse("topo(P1L(p(1),2,1))"));
    let dot = vakint::p_integral::to_dot(atom.as_view(), "phi").unwrap();
    assert_p_backends_agree_atom(&from_dot(&dot), battery_terms(1));
}

/// A bubble in the user's routing agrees with Forcer, and its random
/// reroutings evaluate with RustRed exactly alike.
#[test]
fn a_rerouted_bubble_agrees_with_forcer() {
    assert_p_backends_agree_routed(
        "dot(k(1),p(1))^2*topo(prop(1,edge(1,2),k(1),0,2)*prop(2,edge(2,1),k(1)-p(1),0,1))",
        0..4,
    );
}

/// A one-loop p-integral involves no tabulated master, so no `Oep` guard
/// fires however deep the expansion (M1 review, missing test 5, under
/// RustRed): through `ε^14` both backends agree.
#[test]
fn a_one_loop_p_integral_expands_to_any_depth() {
    let input = "dot(k(1),p(1))^3*topo(P1L(p(1),2,3))";
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let value = rustred_value(&parse(input), mode, 16).unwrap();
        assert!(!value.to_string().contains("Oep"), "{value}");
    }
    assert_p_backends_agree(input, 16);
}

/// Analytic identities with odd ranks and the external momentum reversed
/// (M1 review, missing test 1), evaluated with RustRed alone.
mod rustred_identities {
    use super::*;

    #[test]
    fn one_loop_odd_and_even_rank() {
        let b = bubble("p(1)");
        assert_rustred_p_identity(
            &format!("dot(k(1),p(1))*{b}"),
            &format!("dot(p(1),p(1))/2*{b}"),
        );
        // (k·q)² = (k·q) q²/2 up to scaleless integrals.
        assert_rustred_p_identity(
            &format!("dot(k(1),p(1))^2*{b}"),
            &format!("dot(p(1),p(1))/2*dot(k(1),p(1))*{b}"),
        );
    }

    #[test]
    fn one_loop_with_the_external_momentum_reversed() {
        for head in ["-p(1)", "-p(3)"] {
            let b = bubble(head);
            let p = &head[1..];
            assert_rustred_p_identity(
                &format!("dot(k(1),{p})*{b}"),
                &format!("-dot({p},{p})/2*{b}"),
            );
        }
    }

    #[test]
    fn index_notation_and_linear_dot_arguments_evaluate_alike() {
        let b = bubble("p(1)");
        assert_rustred_p_identity(
            &format!("k(1,mu)*p(1,mu)*{b}"),
            &format!("dot(k(1),p(1))*{b}"),
        );
        assert_rustred_p_identity(
            &format!("dot(k(1)+p(1),p(1))*{b}"),
            &format!("(dot(k(1),p(1))+dot(p(1),p(1)))*{b}"),
        );
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

    #[test]
    #[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
    fn two_loop_partial_fractions() {
        let a = [1, 2, 1, 3, 1];
        let k = |a, head| kite(a, head);
        assert_rustred_p_identity(
            &format!("dot(k(1),p(1))*{}", k(a, "p(1)")),
            &format!(
                "1/2*{}+1/2*dot(p(1),p(1))*{}-1/2*{}",
                k([0, 2, 1, 3, 1], "p(1)"),
                k(a, "p(1)"),
                k([1, 1, 1, 3, 1], "p(1)")
            ),
        );
        assert_rustred_p_identity(
            &format!("dot(k(1),k(2))*{}", k(a, "p(1)")),
            &format!(
                "1/2*{}+1/2*{}-1/2*{}",
                k([0, 2, 1, 3, 1], "p(1)"),
                k([1, 2, 0, 3, 1], "p(1)"),
                k([1, 2, 1, 2, 1], "p(1)")
            ),
        );
    }

    #[test]
    #[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
    fn two_and_three_loops_with_the_external_momentum_reversed() {
        let a = [1, 2, 1, 3, 1];
        assert_rustred_p_identity(
            &format!("dot(k(1),p(2))*{}", kite(a, "-p(2)")),
            &format!(
                "1/2*{}-1/2*dot(p(2),p(2))*{}-1/2*{}",
                kite([1, 1, 1, 3, 1], "-p(2)"),
                kite(a, "-p(2)"),
                kite([0, 2, 1, 3, 1], "-p(2)")
            ),
        );
        let no = |changes: &[(usize, i64)]| {
            let mut powers = [1i64; 8];
            for &(line, power) in changes {
                powers[line - 1] = power;
            }
            let powers = powers.map(|power| power.to_string()).join(",");
            format!("topo(P3L_NO(-p(1),{powers}))")
        };
        assert_rustred_p_identity(
            &format!("dot(k(1),p(1))*{}", no(&[])),
            &format!(
                "1/2*{}+1/2*dot(p(1),p(1))*{}-1/2*{}",
                no(&[(1, 0)]),
                no(&[]),
                no(&[(5, 0)])
            ),
        );
    }
}

/// The kite and every non-scaleless labelled sector of the two-loop family,
/// each times a rank-0-to-4 numerator.
#[test]
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn the_kite_and_every_two_loop_sector_agree_with_forcer() {
    let numerator = "(3+5*dot(k(1),k(2))+7*dot(k(2),p(1))^2+11*dot(k(1),p(1))*dot(k(2),p(1))^3)";
    let sectors = nonzero_sectors("P2L").unwrap();
    assert!(sectors.contains(&vec![1; 5]), "{sectors:?}");
    for powers in sectors {
        let powers = powers.iter().map(i64::to_string).collect::<Vec<_>>();
        assert_p_backends_agree(
            &format!("{numerator}*topo(P2L(p(1),{}))", powers.join(",")),
            battery_terms(2),
        );
    }
}

/// The kite (`Σa = 5`, odd) carries its Wick sign once (critique-phys SF3).
#[test]
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn the_kite_carries_its_wick_sign_once() {
    let Some(form) = forcer_form() else {
        eprintln!("FORM 5 is unavailable: skipped the Wick sign of the kite");
        return;
    };
    let input = "topo(P2L(p(1),1,1,1,1,1))";
    assert_p_backends_agree(input, battery_terms(2));
    let atom = parse(input);
    let mode = RustRedMasterMode::Retain;
    let rustred = rustred_value(&atom, mode, battery_terms(2)).unwrap();
    let forcer = live_forcer_value(&form, &atom, mode, battery_terms(2)).unwrap();
    assert!(values_agree(&rustred, &(-forcer), mode).is_err());
}

/// The three-loop tops LA, BE and NO: plain, with a dot, and with numerators
/// that reach their ISPs.
#[test]
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn the_three_loop_tops_with_dots_and_isps_agree_with_forcer() {
    for family in ["P3L_LA", "P3L_BE", "P3L_NO"] {
        for (numerator, dotted) in [
            ("1", None),
            ("1", Some(2)),
            ("dot(k(1),k(2))*dot(k(2),k(3))", None),
            ("dot(k(1),p(1))*dot(k(3),p(1))+dot(k(1),k(3))", Some(5)),
        ] {
            let mut powers = [1i64; 8];
            if let Some(line) = dotted {
                powers[line - 1] = 2;
            }
            let powers = powers.map(|power| power.to_string()).join(",");
            assert_p_backends_agree(
                &format!("({numerator})*topo({family}(p(1),{powers}))"),
                battery_terms(3),
            );
        }
    }
}

/// Odd ranks with the head `-p(3)` at two and three loops.
#[test]
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn odd_rank_two_and_three_loop_integrals_with_a_reversed_head_agree_with_forcer() {
    for (input, loops) in [
        ("dot(k(1),p(3))*topo(P2L(-p(3),1,1,1,1,1))", 2),
        ("dot(k(2),p(3))^3*topo(P2L(-p(3),1,2,1,1,1))", 2),
        ("dot(k(1),p(3))*dot(k(2),k(3))*topo(P3L_BE(-p(3),1,1,1,1,1,1,1,1))", 3),
    ] {
        assert_p_backends_agree(input, battery_terms(loops));
    }
}

/// A kite, a bubble of another head and a vacuum tadpole in one expression.
#[test]
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn a_kite_with_mixed_heads_and_a_vacuum_tadpole_agrees_with_forcer() {
    assert_p_backends_agree(
        "dot(k(1),p(1))*topo(P2L(p(1),1,1,1,1,1)) + dot(k(1),p(3))^2*topo(P1L(-p(3),1,2)) \
         + topo(prop(1,edge(1,1),k(1),muvsq,2))",
        battery_terms(2),
    );
}

/// The compact FeynKit DOT kite of `docs/p_integrals.md` (section 3), with a
/// dot on one line.
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

/// DOT ingress at two and three loops: the DOT kite and an exported NO top.
#[test]
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn dot_graphs_agree_with_forcer() {
    assert_p_backends_agree_atom(&from_dot(DOT_KITE), battery_terms(2));
    let no = canonical(&parse("topo(P3L_NO(p(1),1,1,1,1,1,1,1,2))"));
    let dot = vakint::p_integral::to_dot(no.as_view(), "phi").unwrap();
    assert_p_backends_agree_atom(&from_dot(&dot), battery_terms(3));
}

/// A kite and a three-loop top in the user's routing, rerouted.
#[test]
#[ignore = "M2-V5 STAGE 2: needs the built-in 1L-3L p library (M2-V4)"]
fn rerouted_two_and_three_loop_integrals_agree_with_forcer() {
    assert_p_backends_agree_routed(
        "dot(k(1),p(1))*dot(k(2),p(1))*topo(prop(1,edge(1,3),k(1),0,1)*prop(2,edge(1,4),p(1)-k(1),0,2)\
         *prop(3,edge(3,4),k(2),0,1)*prop(4,edge(3,2),k(1)-k(2),0,3)*prop(5,edge(4,2),p(1)-k(1)+k(2),0,1))",
        0..2,
    );
    let be: Atom =
        parse("dot(k(1),p(1))*dot(k(2),k(3))") * canonical(&parse("topo(P3L_BE(p(1),1,2,1,1,1,1,1,1))"));
    assert_p_backends_agree_routed_atom(&be, 0..2);
}
