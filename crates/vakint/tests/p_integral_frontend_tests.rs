//! The p-integral frontend (`docs/p_integrals.md`, sections 1-3): kind
//! inference, validation, exact normalization, scaleless zeros, matching onto
//! the 1L-4L registry for any routing, and the DOT round trip.
//!
//! Evaluation is not tested here: it belongs to the Forcer and RustRed
//! backends, which share the canonical forms checked below.

#![cfg(feature = "feynkit-ingress")]

#[path = "support/p_integrals.rs"]
mod p_integrals;

use p_integrals::*;

const BUBBLE: &str = "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,2))";

/// The two-loop kite with powers that break all of its symmetries.
const KITE: &str = "topo(prop(1,edge(1,3),k(1),0,1)*prop(2,edge(1,4),p(1)-k(1),0,2)\
    *prop(3,edge(3,4),k(2),0,1)*prop(4,edge(3,2),k(1)-k(2),0,3)*prop(5,edge(4,2),p(1)-k(1)+k(2),0,1))";

const FOUR_LOOP_TOPS: [&str; 11] = [
    "P4L_LA4", "P4L_J1", "P4L_J2", "P4L_J4", "P4L_J5", "P4L_J6", "P4L_J11", "P4L_HAHA", "P4L_NO1",
    "P4L_NO2", "P4L_NO3",
];

mod one_loop {
    use super::*;

    #[test]
    fn the_bubble_is_the_one_loop_family() {
        assert_p_family(BUBBLE, "P1L", 0);
        assert_p_same_canonical(
            BUBBLE,
            "topo(prop(7,edge(0,5),-k(2)+p(1),0,2)*prop(3,edge(0,5),k(2),0,1))",
        );
    }

    #[test]
    fn routing_invariance() {
        assert_p_reroutings_invariant(&format!("dot(k(1),p(1))^3*dot(k(1),k(1))*{BUBBLE}"), 0..12);
    }

    #[test]
    fn a_pinched_bubble_vanishes() {
        assert_p_scaleless("topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,0))");
        assert_p_scaleless(
            "dot(k(1),p(1))*topo(prop(1,edge(1,2),k(1),0,-1)*prop(2,edge(2,1),k(1)-p(1),0,2))",
        );
        assert_p_scaleless("topo(P1L(p(1),3,0))");
    }

    #[test]
    fn the_dot_round_trip_keeps_the_family() {
        assert_p_dot_round_trip(BUBBLE);
    }
}

mod two_loop {
    use super::*;

    #[test]
    fn the_kite_and_its_sectors_are_identified() {
        assert_p_family(KITE, "P2L", 0);
        // The spec example: the user's head p(3) and its own routing.
        assert_p_family(
            "dot(k(1),p(3))^2*topo(prop(1,edge(1,3),k(1),0,2)*prop(2,edge(1,4),p(3)-k(1),0,1)\
             *prop(3,edge(3,4),k(2),0,1)*prop(4,edge(3,2),k(1)-k(2),0,1)*prop(5,edge(4,2),p(3)-k(1)+k(2),0,1))",
            "P2L",
            0,
        );
        // Contracting the middle line (2) leaves two bubbles in series; an
        // outer line leaves a bubble with a one-loop insertion.
        assert_p_family("topo(P2L(p(1),1,0,1,1,1))", "P2L", 1);
        assert_p_family("topo(P2L(p(1),0,1,1,1,1))", "P2L", 1);
        assert_p_family("topo(P2L(p(1),1,1,1,1,0))", "P2L", 1);
    }

    #[test]
    fn routing_invariance() {
        assert_p_routing_invariant(
            &format!("dot(k(1),p(1))^2*dot(k(2),p(1))*dot(k(1),k(2))*{KITE}"),
            0..12,
        );
        assert_p_routing_invariant("dot(k(2),p(1))^2*topo(P2L(p(1),1,2,0,3,1))", 0..6);
    }

    #[test]
    fn the_user_head_is_kept() {
        assert_p_same_canonical(
            &format!("dot(k(1),p(3))*{}", KITE.replace("p(1)", "p(3)")),
            &format!(
                "dot(k(1),p(3))*{}",
                KITE.replace("p(1)", "p(3)").replace("prop(2,", "prop(12,")
            ),
        );
        let renamed = canonical(&parse(&KITE.replace("p(1)", "p(3)")));
        assert!(!renamed.to_string().contains("p(1)") && renamed.to_string().contains("p(3)"));
    }

    #[test]
    fn scaleless_sectors_vanish() {
        // A line pinched next to an external vertex closes a vacuum loop.
        assert_p_scaleless("topo(P2L(p(1),0,1,1,0,1))");
        assert_p_scaleless("topo(P2L(p(1),1,0,0,0,1))");
        assert_p_scaleless("dot(k(1),k(2))*topo(P2L(p(1),0,0,1,1,1))");
    }

    #[test]
    fn the_dot_round_trip_keeps_the_family() {
        assert_p_dot_round_trip(KITE);
        assert_p_dot_round_trip("topo(P2L(p(1),1,0,1,1,1))");
        assert_p_dot_round_trip("topo(P2L(p(1),0,1,1,1,1))");
        assert_p_dot_round_trip("topo(P2L(p(1),2,1,1,1,0))");
    }

    #[test]
    fn the_spec_dot_example_is_the_kite_with_one_dot() {
        let imported = from_dot(
            r#"digraph p2l_kite {
                exte0 [style=invis]; exte1 [style=invis];
                exte0 -> v1 [id=0 particle="phi"];
                v2 -> exte1 [id=1 particle="phi"];
                v2 -> v4 [id=2 particle="phi" lmb_id=0];
                v3 -> v4 [id=3 particle="phi" lmb_id=1];
                v4 -> v1 [id=4 particle="phi"];
                v3 -> v2 [id=5 particle="phi"];
                v1 -> v3 [id=6 particle="phi" power=2];
            }"#,
        );
        assert_eq!(
            canonical(&imported),
            canonical(&parse(
                "topo(prop(1,edge(1,3),k(1),0,2)*prop(2,edge(1,4),p(1)-k(1),0,1)\
                 *prop(3,edge(3,4),k(2),0,1)*prop(4,edge(3,2),k(1)-k(2),0,1)*prop(5,edge(4,2),p(1)-k(1)+k(2),0,1))"
            )),
        );
    }
}

mod three_loop {
    use super::*;

    #[test]
    fn the_three_tops_and_their_sectors_are_identified() {
        for family in ["P3L_LA", "P3L_BE", "P3L_NO"] {
            assert_p_family(&format!("topo({family}(p(1),1,1,1,1,1,1,1,1))"), family, 0);
        }
        // The first family in table order whose contraction is this graph.
        assert_p_family("topo(P3L_NO(p(1),1,1,1,1,1,1,1,0))", "P3L_BE", 1);
    }

    #[test]
    fn routing_invariance() {
        for family in ["P3L_LA", "P3L_BE", "P3L_NO"] {
            let top = distinct_power_top(family, 8, "dot(k(1),p(1))^2*dot(k(2),k(3))");
            assert_p_atom_routing_invariant(&top, 0..6);
        }
    }

    #[test]
    fn the_dot_round_trip_keeps_the_family() {
        for family in ["P3L_LA", "P3L_BE", "P3L_NO"] {
            assert_p_dot_round_trip(&format!("topo({family}(p(1),1,1,1,1,1,1,1,1))"));
        }
        assert_p_dot_round_trip("topo(P3L_NO(p(1),1,1,1,1,1,1,1,0))");
        assert_p_dot_round_trip("topo(P3L_LA(p(1),1,0,1,1,1,1,1,1))");
        assert_p_dot_round_trip("topo(P3L_BE(p(1),2,1,1,1,0,1,1,1))");
    }
}

mod four_loop {
    use super::*;

    #[test]
    fn the_eleven_tops_are_their_own_families() {
        for family in FOUR_LOOP_TOPS {
            assert_p_family(
                &format!("topo({family}(p(1),1,1,1,1,1,1,1,1,1,1,1))"),
                family,
                0,
            );
        }
    }

    #[test]
    fn routing_invariance() {
        for family in FOUR_LOOP_TOPS {
            let top = distinct_power_top(family, 11, "dot(k(1),p(1))*dot(k(2),k(4))");
            assert_p_atom_routing_invariant(&top, 0..3);
        }
    }

    #[test]
    fn the_dot_round_trip_keeps_the_family() {
        for family in FOUR_LOOP_TOPS {
            assert_p_dot_round_trip(&format!("topo({family}(p(1),1,1,1,1,1,1,1,1,1,1,1))"));
        }
        assert_p_dot_round_trip("topo(P4L_HAHA(p(1),1,1,1,1,0,1,1,1,1,1,1))");
        assert_p_dot_round_trip("topo(P4L_NO1(p(1),1,1,0,1,1,1,1,1,1,1,1))");
        assert_p_dot_round_trip("topo(P4L_LA4(p(1),1,1,1,1,1,1,1,1,0,0,1))");
    }
}

mod normalization {
    use super::*;

    #[test]
    fn a_tree_line_becomes_a_numerator_factor() {
        assert_p_same_canonical(
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,2)*prop(3,edge(3,1),p(1),0,2))",
            &format!("dot(p(1),p(1))^-2*{BUBBLE}"),
        );
    }

    #[test]
    fn lines_with_identical_quadrics_merge() {
        // Line 1 split by a two-valent vertex, written once against its flow.
        assert_p_same_canonical(
            "topo(prop(1,edge(1,3),k(1),0,1)*prop(4,edge(2,3),-k(1),0,2)*prop(2,edge(2,1),k(1)-p(1),0,2))",
            "topo(prop(1,edge(1,2),k(1),0,3)*prop(2,edge(2,1),k(1)-p(1),0,2))",
        );
    }

    #[test]
    fn a_negative_power_becomes_a_numerator() {
        assert_p_same_canonical(
            &KITE.replace("prop(3,edge(3,4),k(2),0,1)", "prop(3,edge(3,4),k(2),0,-2)"),
            &format!(
                "dot(k(2),k(2))^2*{}",
                KITE.replace("prop(3,edge(3,4),k(2),0,1)", "prop(3,edge(3,4),k(2),0,0)")
            ),
        );
    }

    #[test]
    fn the_short_form_is_the_family_long_form() {
        assert_p_same_canonical(
            "dot(k(1),p(1))*topo(P2L(-p(1),1,2,1,1,1))",
            "dot(k(1),p(1))*topo(prop(1,edge(2,4),k(1),0,1)*prop(2,edge(3,4),k(2),0,2)*prop(3,edge(4,1),k(1)+k(2),0,1)\
             *prop(4,edge(3,2),k(1)-p(1),0,1)*prop(5,edge(1,3),k(1)+k(2)-p(1),0,1))",
        );
    }

    #[test]
    fn massless_vacuum_integrals_vanish() {
        assert_p_scaleless("topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1),0,1))");
    }
}

mod errors {
    use super::*;

    #[test]
    fn invalid_p_integrals_are_rejected_with_their_reason() {
        assert_p_error(
            "topo(prop(1,edge(1,2),k(1),msq,1)*prop(2,edge(2,1),k(1)-p(1),0,1))",
            "p-integrals are massless; prop 1 has mass",
        );
        assert_p_error(
            "topo(prop(1,edge(1,2),k(1)+p(2),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))",
            "a p-integral has one external momentum; found",
        );
        assert_p_error(
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-2*p(1),0,1))",
            "momentum conservation",
        );
        assert_p_error(
            "topo(prop(1,edge(1,2),k(1),0,1/2)*prop(2,edge(2,1),k(1)-p(1),0,1))",
            "powers must be integers",
        );
        assert_p_error(
            "topo(prop(1,edge(1,2),k(1)+k(2),0,1)*prop(2,edge(2,1),k(1)+k(2)-p(1),0,1))",
            "the routing does not match the graph's loop count",
        );
        assert_p_error(
            &format!("dot(k(1),p(2))*{BUBBLE}"),
            "tensor p-integrals are not supported yet",
        );
        assert_p_error(
            &format!("k(1,mu)*{BUBBLE}"),
            "tensor p-integrals are not supported yet",
        );
    }

    #[test]
    fn five_loop_inputs_are_outside_the_registry() {
        assert_p_error(
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(1,2),k(2),0,1)*prop(3,edge(1,2),k(3),0,1)\
             *prop(4,edge(1,2),k(4),0,1)*prop(5,edge(1,2),k(5),0,1)\
             *prop(6,edge(2,1),k(1)+k(2)+k(3)+k(4)+k(5)-p(1),0,1))",
            "p-integral not covered by the 1L-4L registry",
        );
    }
}
