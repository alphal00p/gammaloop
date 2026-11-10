//! Tier B of the Forcer battery (M2-V5; M2 plan, section 8, gate G7): the
//! RustRed p backend against a live Forcer (FORM 5), through vakint's public
//! API, exactly (`Retain` in the tag basis, `Symbolic`) and numerically
//! (`Numerical`, 40 digits) (`support/p_backends.rs`). Without FORM 5 the live
//! comparisons report a skip, or fail with `VAKINT_REQUIRE_FORM5=1` (the gate
//! environment); the RustRed identities, the reroutings and the fail-closed
//! checks of the Forcer side need no FORM.
//!
//! The RustRed side is `Vakint::new()` with its built-in 1L-3L p library
//! (M2-V4); at four loops (M3-V9c, [`four_loops`]) an engine on M3's scope
//! 1-4 package listed in an external-data root (`VAKINT_TEST_P_LIBRARY_1_4`;
//! skipped without it), against a live Forcer that must be there.

#![cfg(feature = "feynkit-ingress")]

#[path = "support/p_backends.rs"]
mod p_backends;
#[path = "support/p_integrals.rs"]
mod p_integrals;
#[path = "support/p_library4.rs"]
mod p_library4;
#[path = "support/p_rustred.rs"]
mod p_rustred;
#[path = "support/p_tierc4.rs"]
mod p_tierc4;

use p_backends::*;
use p_integrals::{canonical, from_dot, parse};
use p_rustred::rustred_p_evaluate;
use symbolica::atom::Atom;
use vakint::p_integral::forcer_fixtures::nonzero_sectors;
use vakint::{
    LoopNormalizationFactor, OrderedCacheMode, OrderedCacheSettings, RustRedMasterMode,
    VakintError, VakintSettings,
};

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
        "dot(k(1),p(3))*dot(k(1),k(1))*topo(P1L(-p(3),2,3))",
    ] {
        assert_p_backends_agree(input, battery_terms(1));
    }
}

/// The Wick sign `(-1)^(Σa+r)` is applied once (critique-phys SF3): RustRed
/// agrees with Forcer for an odd rank (`Σa = 2`, `r = 1`) and an odd `Σa`
/// (`(2,1)`, `r = 0`), where a sign applied twice or never gives minus
/// Forcer's nonzero value. Tier A (`rustred_battery_tests.rs`) pins that the
/// sign is visible: Forcer's frozen values with and without it differ.
#[test]
fn an_odd_bubble_carries_its_wick_sign_once() {
    for input in ["dot(k(1),p(1))*topo(P1L(p(1),1,1))", "topo(P1L(p(1),2,1))"] {
        assert_p_backends_agree(input, battery_terms(1));
    }
}

/// One input per loop count for the concurrent path: an odd-rank dotted
/// bubble, a dotted kite with a rank-3 numerator, and the dotted NO top,
/// which reduces onto both masters `no` and `t105`.
const CONCURRENT_INPUTS: [&str; 3] = [
    "dot(k(1),p(1))^3*topo(P1L(p(1),2,3))",
    "dot(k(1),p(1))*dot(k(2),p(1))^2*topo(P2L(p(1),1,2,1,1,1))",
    "topo(P3L_NO(p(1),2,1,1,1,1,1,1,1))",
];

/// The concurrent path (M2 exit criterion 2; every other p test runs on one
/// core): with `n_cores = 4`, `modular_sampled` on four fixed workers and
/// `modular_polynomial` following `n_cores` agree with Forcer, at the
/// deepest loop count of [`P_BATTERY_LOOPS`].
#[test]
fn four_cores_under_both_modular_arithmetics_agree_with_forcer() {
    let loops = *P_BATTERY_LOOPS.end();
    for arithmetic in ["modular_sampled:threads=4", "modular_polynomial"] {
        assert_p_backends_agree_under(
            CONCURRENT_INPUTS[loops - 1],
            battery_terms(loops),
            arithmetic,
            4,
        );
    }
}

/// The pair of the step-cache test: the kite, then the kite with line 5
/// dotted, whose reduction passes through steps of the first.
const CACHE_PAIR: [&str; 2] = ["topo(P2L(p(1),1,1,1,1,1))", "topo(P2L(p(1),1,1,1,1,2))"];

/// The p library's Ordered step cache (M2 exit criterion 4, the vakint
/// half), FORM-free. The pair evaluated read-write into a directory under
/// `exact`, the second integral reusing steps of the first (step hits), is
/// replayed read-only by fresh engines under `modular_polynomial` and
/// `modular_sampled`, which compute nothing (step hits only: no miss, no
/// insert). Every value equals a cache-off engine's. The counts are printed
/// for the notes.
#[test]
fn the_p_step_cache_is_reused_and_replays_under_every_arithmetic() {
    let loops = 2;
    let pair = CACHE_PAIR.map(parse);
    let settings = |arithmetic, cache| VakintSettings {
        ordered_cache: cache,
        ..rustred_p_settings_with_arithmetic(
            RustRedMasterMode::Retain,
            LoopNormalizationFactor::GScheme,
            battery_terms(loops),
            arithmetic,
            1,
        )
    };
    let off = settings("exact", OrderedCacheSettings::off());
    let expected = pair.each_ref().map(|input| {
        let value = rustred_p_evaluate(&new_p_battery_engine(), &off, input);
        assert!(!value.is_zero(), "{input}");
        value
    });

    let directory = tempfile::tempdir().unwrap();
    let writing = settings("exact", OrderedCacheSettings::read_write(directory.path()));
    let writer = new_p_battery_engine();
    let mut before = PCacheCounts::default();
    for (index, (input, expected)) in pair.iter().zip(&expected).enumerate() {
        assert_eq!(&rustred_p_evaluate(&writer, &writing, input), expected);
        let counts = p_cache_counts(&writer, OrderedCacheMode::ReadWrite);
        let own = counts - before;
        eprintln!(
            "[p cache] exact read-write, integral {}: {own:?}",
            index + 1
        );
        assert!(own.misses > 0 && own.inserts > 0, "{input}: {own:?}");
        if index == 1 {
            assert!(own.hits > 0, "{input} reuses no step of {}", pair[0]);
        }
        before = counts;
    }
    let program = p_cache_program(&writer, OrderedCacheMode::ReadWrite);
    let mut files = std::fs::read_dir(program.directory.unwrap())
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    files.sort();
    eprintln!("[p cache] files written: {files:?}");
    assert!(!files.is_empty(), "the read-write cache wrote no file");

    for arithmetic in ["modular_polynomial", "modular_sampled"] {
        let reading = settings(
            arithmetic,
            OrderedCacheSettings::read_only(directory.path()),
        );
        let reader = new_p_battery_engine();
        for (input, expected) in pair.iter().zip(&expected) {
            let value = rustred_p_evaluate(&reader, &reading, input);
            assert_eq!(&value, expected, "{arithmetic}");
        }
        let counts = p_cache_counts(&reader, OrderedCacheMode::ReadOnly);
        eprintln!("[p cache] {arithmetic} read-only replay: {counts:?}");
        assert_eq!((counts.misses, counts.inserts), (0, 0), "{arithmetic}");
        assert!(counts.hits > 0, "{arithmetic}: {counts:?}");
    }
}

/// A FORM that does not run Forcer.
const BROKEN_FORM: &str = "/nonexistent/vakint-battery/form";

/// The Forcer side fails closed (review round 1, MF1): with a FORM path that
/// does not run Forcer, on the RustRed engine (whose p library evaluates
/// these p-integrals), it is `ForcerUnavailable` in every mode, never a
/// value, and an agreement against it is not decided. FORM-free.
#[test]
fn the_forcer_side_fails_closed_with_a_broken_form() {
    let terms = battery_terms(1);
    for input in [
        "dot(k(1),p(1))*topo(P1L(p(1),1,1))",
        "dot(k(1),p(1))^2*topo(P1L(p(1),2,1)) + topo(prop(1,edge(1,1),k(1),muvsq,2))",
    ] {
        let input = parse(input);
        for mode in MODES {
            let value = rustred_value(&input, mode, terms)
                .unwrap_or_else(|error| panic!("the engine reduces {input}: {error}"));
            assert!(!value.is_zero(), "{input} ({mode})");
            match forcer_side_value(p_battery_engine(), BROKEN_FORM, &input, mode, terms) {
                Err(VakintError::ForcerUnavailable { .. }) => {}
                other => panic!("{input} ({mode}): expected ForcerUnavailable, got {other:?}"),
            }
        }
        let broken = |mode| {
            forcer_side_value(p_battery_engine(), BROKEN_FORM, &input, mode, terms)
                .map_err(|error| error.to_string())
        };
        match agreement(&input, terms, broken) {
            Err(Disagreement::Failed(_)) => {}
            other => panic!("{input}: expected an undecided agreement, got {other:?}"),
        }
    }
}

/// The Forcer side splits a mixed input: the p-integral terms go to Forcer,
/// the vacuum terms to RustRed; a term without an integral, a product of
/// both kinds and an input without a p-integral are errors.
#[test]
fn the_forcer_side_splits_p_and_vacuum_terms() {
    let bubble = parse("dot(k(1),p(1))*topo(P1L(p(1),2,1))");
    let tadpole = parse("topo(prop(1,edge(1,1),k(1),muvsq,2))");
    let (p_terms, vacuum_terms) = p_and_vacuum_terms(&(&bubble + &tadpole)).unwrap();
    assert_eq!((p_terms, vacuum_terms), (bubble.clone(), tadpole.clone()));
    assert_eq!(
        p_and_vacuum_terms(&bubble).unwrap(),
        (bubble.clone(), Atom::Zero)
    );
    for input in [&tadpole, &(&bubble + &Atom::num(1)), &(&bubble * &tadpole)] {
        assert!(p_and_vacuum_terms(input).is_err(), "{input}");
    }
}

/// With `VAKINT_REQUIRE_FORM5=1` a missing FORM 5 fails instead of skipping.
#[test]
#[should_panic(expected = "VAKINT_REQUIRE_FORM5=1 requires it")]
fn a_missing_form5_fails_when_it_is_required() {
    form_or_skip(None, true, "a comparison");
}

/// Without the requirement, a missing FORM 5 is a skip, and a found one is
/// used.
#[test]
fn a_missing_form5_is_a_skip_unless_it_is_required() {
    assert_eq!(form_or_skip(None, false, "a comparison"), None);
    assert_eq!(
        form_or_skip(Some("form".into()), true, "a comparison"),
        Some("form".into())
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
/// fires at any depth (M1 review, missing test 5, under RustRed): through
/// `ε^9`, past the three-loop boundary of
/// [`the_three_loop_oep_guard_fails_closed_under_rustred`], both backends
/// agree.
#[test]
fn a_one_loop_p_integral_has_no_master_order_boundary() {
    let input = "dot(k(1),p(1))^3*topo(P1L(p(1),2,3))";
    for mode in [RustRedMasterMode::Symbolic, RustRedMasterMode::Numerical] {
        let value = rustred_value(&parse(input), mode, ONE_LOOP_DEPTH).unwrap();
        assert!(!value.to_string().contains("Oep"), "{value}");
    }
    assert_p_backends_agree(input, ONE_LOOP_DEPTH);
}

/// Laurent terms of [`a_one_loop_p_integral_has_no_master_order_boundary`].
const ONE_LOOP_DEPTH: i64 = 11;

/// The `Oep` guard under RustRed at three loops (M1 review, missing test 5):
/// `Master(no)` is tabulated through `ε^7` (`Masters.prc:354`, `ALARMno` at
/// `ep^8`), so the NO top passes through 11 terms (`ε^-3..ε^7`), equal to
/// Forcer's value, and one term more is `MissingMasterOrders` with RustRed.
/// The boundary is pinned so that a guard firing too early fails too.
#[test]
fn the_three_loop_oep_guard_fails_closed_under_rustred() {
    assert_p_backends_oep_boundary("topo(P3L_NO(p(1),1,1,1,1,1,1,1,1))", 11);
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

/// The kite (`Σa = 5`, odd) carries its Wick sign once (critique-phys SF3):
/// RustRed agrees with Forcer, whose nonzero value a doubled sign would
/// negate (Tier A pins the Wick-visible cases of the frozen batteries).
#[test]
fn the_kite_carries_its_wick_sign_once() {
    assert_p_backends_agree("topo(P2L(p(1),1,1,1,1,1))", battery_terms(2));
}

/// The three-loop tops LA, BE and NO: plain, with a dot, and with numerators
/// that reach their ISPs, summed into one input per family with distinct
/// coefficients (one Forcer run per family and mode).
#[test]
fn the_three_loop_tops_with_dots_and_isps_agree_with_forcer() {
    for family in ["P3L_LA", "P3L_BE", "P3L_NO"] {
        assert_p_backends_agree(
            &format!(
                "(2+3*dot(k(1),k(2))*dot(k(2),k(3)))*topo({family}(p(1),1,1,1,1,1,1,1,1)) \
                 + 5*topo({family}(p(1),1,2,1,1,1,1,1,1)) \
                 + (7*dot(k(1),p(1))*dot(k(3),p(1))+11*dot(k(1),k(3)))\
                 *topo({family}(p(1),1,1,1,1,2,1,1,1))"
            ),
            battery_terms(3),
        );
    }
}

/// Odd ranks with the head `-p(3)` at two and three loops.
#[test]
fn odd_rank_two_and_three_loop_integrals_with_a_reversed_head_agree_with_forcer() {
    for (input, loops) in [
        ("dot(k(1),p(3))*topo(P2L(-p(3),1,1,1,1,1))", 2),
        ("dot(k(2),p(3))^3*topo(P2L(-p(3),2,1,1,1,1))", 2),
        (
            "dot(k(1),p(3))*dot(k(2),k(3))*topo(P3L_BE(-p(3),1,1,1,1,1,1,1,1))",
            3,
        ),
    ] {
        assert_p_backends_agree(input, battery_terms(loops));
    }
}

/// A kite, a bubble of another head and a vacuum tadpole in one expression.
#[test]
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
fn dot_graphs_agree_with_forcer() {
    assert_p_backends_agree_atom(&from_dot(DOT_KITE), battery_terms(2));
    let no = canonical(&parse("topo(P3L_NO(p(1),1,1,1,1,1,1,1,2))"));
    let dot = vakint::p_integral::to_dot(no.as_view(), "phi").unwrap();
    assert_p_backends_agree_atom(&from_dot(&dot), battery_terms(3));
}

/// A kite and a three-loop top in the user's routing, rerouted.
#[test]
fn rerouted_two_and_three_loop_integrals_agree_with_forcer() {
    assert_p_backends_agree_routed(
        "dot(k(1),p(1))*dot(k(2),p(1))*topo(prop(1,edge(1,3),k(1),0,1)*prop(2,edge(1,4),p(1)-k(1),0,2)\
         *prop(3,edge(3,4),k(2),0,1)*prop(4,edge(3,2),k(1)-k(2),0,3)*prop(5,edge(4,2),p(1)-k(1)+k(2),0,1))",
        0..2,
    );
    let be: Atom = parse("dot(k(1),p(1))*dot(k(2),k(3))")
        * canonical(&parse("topo(P3L_BE(p(1),1,2,1,1,1,1,1,1))"));
    assert_p_backends_agree_routed_atom(&be, 0..2);
}

/// Tier B at four loops (M3-V9c; M3 plan, G7): RustRed on M3's scope 1-4
/// package, listed in an external-data root (the route that ships it, D-SIZE
/// (B)), against a live Forcer, through the deepest depth (at most `ε^6`) at
/// which Forcer's value finalizes. With the package, a missing FORM 5 fails.
mod four_loops {
    use super::p_backends::live_forcer_value;
    use super::p_integrals::{forcer_form, parse};
    use super::p_library4::PLibraryRoute;
    use super::p_tierc4::{agreement4, live_terms, p_library4_engine, resource_stop, settings4};
    use std::time::Duration;
    use symbolica::utils::resources::ResourceLimitError;
    use vakint::{RustRedEvaluationError, VakintError};

    /// `(input, arithmetic, cores)`: an odd rank on a dotted top with the head
    /// written `-p(3)`; a dotted top and a two-loop kite of another head in
    /// one expression under `modular_sampled` on four workers; a top with ISP
    /// numerators under `modular_polynomial` on four cores (the concurrent
    /// path at four loops). No vacuum term: at the four-loop depth a one-loop
    /// tadpole's `Symbolic` series runs for minutes (M3-V9c notes).
    const CASES: [(&str, Option<&str>, usize); 3] = [
        (
            "dot(k(1),p(3))*topo(P4L_HAHA(-p(3),1,1,1,1,1,1,1,2,1,1,1))",
            None,
            1,
        ),
        (
            "topo(P4L_NO3(p(1),1,1,1,1,1,1,1,1,1,1,2)) + dot(k(1),p(3))*topo(P2L(-p(3),1,1,1,1,1))",
            Some("modular_sampled:threads=4"),
            4,
        ),
        (
            "(2+3*dot(k(1),k(2))*dot(k(3),k(4)))*topo(P4L_J1(p(1),1,1,1,1,1,1,1,1,1,1,1))",
            Some("modular_polynomial"),
            4,
        ),
    ];

    #[test]
    fn four_loop_inputs_agree_with_a_live_forcer() {
        let what = "the four-loop Tier B";
        let Some(engine) = p_library4_engine(PLibraryRoute::Listed, what) else {
            return;
        };
        let form = forcer_form().expect("the four-loop Tier B needs FORM 5");
        for (input, arithmetic, cores) in CASES {
            let input = parse(input);
            let terms = live_terms(&form, &input, 4);
            let seconds = agreement4(
                &engine,
                &input,
                |mode| settings4(mode, terms, arithmetic, cores),
                |mode| live_forcer_value(&form, &input, mode, terms),
            )
            .unwrap_or_else(|difference| panic!("{input} ({terms} terms):\n{difference}"));
            eprintln!("[tier B 4L] {input}: {terms} terms, RustRed {seconds:.1} s");
        }
    }

    /// The four-loop tiers tell a resource stop from a failure by its cause
    /// (`p_tierc4::resource_stop`; review vk-tierc4 round 1): a typed time or
    /// RSS limit always counts; the untyped "malformed p library" form in
    /// which vakint's p path may report a stop counts only once the caller's
    /// budget is spent (not spent: a stop memoised by an earlier evaluation,
    /// a failure); every other error is a failure.
    #[test]
    fn resource_stops_are_told_from_failures_by_their_cause() {
        let elapsed = ResourceLimitError::Elapsed {
            elapsed: Duration::from_secs(5),
            limit: Duration::from_secs(5),
        };
        let memory = ResourceLimitError::Memory {
            rss_bytes: 2,
            additional_bytes: 1,
            limit_bytes: 2,
        };
        let unavailable = ResourceLimitError::MonitorUnavailable {
            resource: "RSS",
            message: "no /proc".into(),
        };
        let untyped = |detail: String| {
            VakintError::RustRedEvaluation(RustRedEvaluationError::Reduction {
                detail: format!("malformed p library: the route of sector 0x1ff: {detail}"),
            })
        };
        for spent in [false, true] {
            for limit in [&elapsed, &memory] {
                assert!(resource_stop(
                    &VakintError::ResourceLimit(limit.clone()),
                    spent
                ));
                let typed = RustRedEvaluationError::ResourceLimit(limit.clone());
                assert!(resource_stop(&VakintError::RustRedEvaluation(typed), spent));
                assert_eq!(resource_stop(&untyped(limit.to_string()), spent), spent);
            }
            assert!(!resource_stop(
                &VakintError::ResourceLimit(unavailable.clone()),
                spent
            ));
            assert!(!resource_stop(&untyped(unavailable.to_string()), spent));
            assert!(!resource_stop(&untyped("no route".into()), spent));
        }
    }
}
