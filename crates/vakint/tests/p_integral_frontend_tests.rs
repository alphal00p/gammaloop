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

mod dot_ingress {
    use super::*;
    use vakint::diagram_integral::DiagramIntegralOptions;

    /// A dotted bubble whose two legs both point into the graph.
    const INCOMING_LEGS: &str = r#"digraph bubble {
        exte0 [style=invis]; exte1 [style=invis];
        exte0 -> v1 [id=0 particle="phi"];
        exte1 -> v2 [id=1 particle="phi"];
        v1 -> v2 [id=2 particle="phi" lmb_id=0];
        v1 -> v2 [id=3 particle="phi" power=2];
    }"#;

    #[test]
    fn the_first_leg_carries_p1_into_the_graph() {
        assert_eq!(
            canonical(&from_dot(INCOMING_LEGS)),
            canonical(&parse(
                "topo(prop(2,edge(1,2),k(1),0,1)*prop(3,edge(1,2),p(1)-k(1),0,2))"
            )),
        );
    }

    #[test]
    fn the_external_momentum_can_be_renamed_but_not_accompanied() {
        let renamed = from_dot_with(
            INCOMING_LEGS,
            &DiagramIntegralOptions {
                external_momentum: Some(parse("p(3)")),
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            canonical(&renamed),
            canonical(&parse(
                "topo(prop(2,edge(1,2),k(1),0,1)*prop(3,edge(1,2),p(3)-k(1),0,2))"
            )),
        );
        for options in [
            DiagramIntegralOptions {
                external_momenta: vec![parse("q")],
                ..Default::default()
            },
            DiagramIntegralOptions {
                external_momentum: Some(parse("p(1)+p(2)")),
                ..Default::default()
            },
        ] {
            assert!(from_dot_with(INCOMING_LEGS, &options).is_err());
        }
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
            FIVE_LOOP_SUNSET,
            "p-integral not covered by the 1L-4L registry",
        );
    }

    #[test]
    fn a_routing_that_is_not_unimodular_is_rejected() {
        // 2*k(1) would rescale the measure by 2^-d: never silently accepted.
        assert_p_malformed(
            "topo(prop(1,edge(1,2),2*k(1),0,1)*prop(2,edge(2,1),2*k(1)-p(1),0,1))",
            "no unimodular change of loop momenta",
        );
    }
}

/// The five-loop sunset: a p-integral outside the 1L-4L registry.
const FIVE_LOOP_SUNSET: &str = "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(1,2),k(2),0,1)\
    *prop(3,edge(1,2),k(3),0,1)*prop(4,edge(1,2),k(4),0,1)*prop(5,edge(1,2),k(5),0,1)\
    *prop(6,edge(2,1),k(1)+k(2)+k(3)+k(4)+k(5)-p(1),0,1))";

/// Section 1, "Out of domain" (M1 review F1, missing test 8).
mod out_of_domain {
    use super::*;

    #[test]
    fn out_of_domain_inputs_are_unknown_integrals() {
        assert_p_unknown(
            "topo(prop(1,edge(1,2),k(1),muvsq,1)*prop(2,edge(2,1),k(1)-p(1),muvsq,1))",
            "p-integrals are massless; prop 1 has mass muvsq",
        );
        assert_p_unknown(
            "dot(k(1),p(1))*topo(prop(1,edge(1,2),k(1),0,a)*prop(2,edge(2,1),k(1)-p(1),0,1))",
            "powers must be integers; prop 1 has power a",
        );
        assert_p_unknown(
            "topo(prop(1,edge(1,2),k(1),0,1/2)*prop(2,edge(2,1),k(1)-p(1),0,1))",
            "powers must be integers",
        );
        assert_p_unknown(
            "topo(prop(1,edge(1,2),k(1)+p(2),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))",
            "a p-integral has one external momentum; found",
        );
        // Several heads conserve momentum when every vertex balance is free of
        // k: here p(2)+3*p(1) enters at vertex 2 and leaves at vertex 1.
        assert_p_unknown(
            "topo(prop(1,edge(1,2),k(1)+p(2),0,1)*prop(2,edge(2,1),k(1)-3*p(1),0,1))",
            "a p-integral has one external momentum; found",
        );
        assert_p_unknown(
            FIVE_LOOP_SUNSET,
            "p-integral not covered by the 1L-4L registry",
        );
    }

    #[test]
    fn unknown_vacuum_integrals_are_unchanged() {
        let massive_sunset = "topo(prop(1,edge(1,2),k(1),muvsq,1)*prop(2,edge(2,1),k(1),muvsq,1)\
            *prop(3,edge(1,2),k(2),muvsq,1)*prop(4,edge(2,1),k(2),muvsq,1)*prop(5,edge(1,2),k(3),muvsq,1)\
            *prop(6,edge(2,1),k(1)+k(2)+k(3),muvsq,1))";
        assert_eq!(
            canonical(&parse(massive_sunset)),
            parse(&format!(
                "topo(UNKNOWN({}))",
                &massive_sunset[5..massive_sunset.len() - 1]
            )),
        );
    }

    #[test]
    fn malformed_inputs_are_errors_under_any_setting() {
        assert_p_malformed(
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-2*p(1),0,1))",
            "momentum conservation",
        );
        assert_p_malformed(
            "topo(prop(1,edge(1,2),k(1)+k(2),0,1)*prop(2,edge(2,1),k(1)+k(2)-p(1),0,1))",
            "the routing does not match the graph's loop count",
        );
        assert_p_malformed(
            "topo(prop(1,edge(1,2),k(1),0,1)*prop(1,edge(2,1),k(1)-p(1),0,1))",
            "prop ids must be distinct",
        );
        assert_p_malformed(
            "topo(prop(1,edge(1,-2),k(1),0,1)*prop(2,edge(-2,1),k(1)-p(1),0,1))",
            "vertex ids must be nonnegative",
        );
        assert_p_malformed(
            "topo(prop(1,edge(1,2),k(1)/2,0,1)*prop(2,edge(2,1),k(1)/2-p(1),0,1))",
            "is not an integer combination of k(j) and p(i)",
        );
        assert_p_malformed(
            &format!("dot(k(1),p(2))*{BUBBLE}"),
            "tensor p-integrals are not supported yet",
        );
        // A short form is no unknown integral: out of domain, it stays an error.
        assert_p_malformed("topo(P1L(p(1),1/2,1))", "powers must be integers");
    }

    /// Malformed wins: an input that is out of domain AND malformed is an error
    /// under any setting, never an unknown integral (M2-V1 review round 1).
    #[test]
    fn malformed_inputs_out_of_domain_as_well_are_errors_under_any_setting() {
        // Massive line and repeated prop ids.
        assert_p_malformed(
            "topo(prop(1,edge(1,2),k(1),muvsq,1)*prop(1,edge(2,1),k(1)-p(1),0,1))",
            "prop ids must be distinct",
        );
        // Massive line and two loop momenta on a one-cycle graph.
        assert_p_malformed(
            "topo(prop(1,edge(1,2),k(1)+k(2),muvsq,1)*prop(2,edge(2,1),k(1)+k(2)-p(1),0,1))",
            "the routing does not match the graph's loop count",
        );
        // Non-integer power and a loop-count mismatch.
        assert_p_malformed(
            "topo(prop(1,edge(1,2),k(1)+k(2),0,1/2)*prop(2,edge(2,1),k(1)+k(2)-p(1),0,1))",
            "the routing does not match the graph's loop count",
        );
        // Massive line and a routing of determinant 2.
        assert_p_malformed(
            "topo(prop(1,edge(1,2),2*k(1),muvsq,1)*prop(2,edge(2,1),2*k(1)-p(1),0,1))",
            "no unimodular change of loop momenta",
        );
        // Two heads and a loop momentum that is not conserved.
        assert_p_malformed(
            "topo(prop(1,edge(1,2),k(1)+p(2),0,1)*prop(2,edge(2,1),2*k(1)-3*p(1),0,1))",
            "momentum conservation",
        );
        // A graph outside the registry (five loops) with a routing of
        // determinant 2.
        assert_p_malformed(
            "topo(prop(1,edge(1,2),2*k(1),0,1)*prop(2,edge(1,2),k(2),0,1)*prop(3,edge(1,2),k(3),0,1)\
             *prop(4,edge(1,2),k(4),0,1)*prop(5,edge(1,2),k(5),0,1)\
             *prop(6,edge(2,1),2*k(1)+k(2)+k(3)+k(4)+k(5)-p(1),0,1))",
            "no unimodular change of loop momenta",
        );
    }
}

/// User numerators in the forms of section 1, item 7 (M1 review extras).
mod numerators {
    use super::*;

    const KITE_SYMMETRIC: &str = "topo(prop(1,edge(1,3),k(1),0,1)*prop(2,edge(1,4),p(1)-k(1),0,1)\
        *prop(3,edge(3,4),k(2),0,1)*prop(4,edge(3,2),k(1)-k(2),0,1)*prop(5,edge(4,2),p(1)-k(1)+k(2),0,1))";

    #[test]
    fn index_notation_is_a_scalar_product() {
        assert_p_same_canonical(
            &format!("k(1,mu)*p(1,mu)*{BUBBLE}"),
            &format!("dot(k(1),p(1))*{BUBBLE}"),
        );
        assert_p_same_canonical(
            &format!("k(1,mu)*p(1,mu)*k(2,nu)*k(1,nu)*{KITE}"),
            &format!("dot(k(1),p(1))*dot(k(1),k(2))*{KITE}"),
        );
    }

    #[test]
    fn scalar_products_of_combinations_are_linear() {
        for numerator in ["-dot(k(1),-p(1))", "1/2*dot(2*k(1),p(1))", "dot(p(1),k(1))"] {
            assert_p_same_canonical(
                &format!("{numerator}*{BUBBLE}"),
                &format!("dot(k(1),p(1))*{BUBBLE}"),
            );
        }
        assert_p_same_canonical(
            &format!("dot(k(1)+k(2),p(1))*{KITE_SYMMETRIC}"),
            &format!("(dot(k(1),p(1))+dot(k(2),p(1)))*{KITE_SYMMETRIC}"),
        );
    }

    #[test]
    fn a_scalar_product_of_a_sum_with_the_external_momentum_is_linear() {
        assert_p_same_canonical(
            &format!("dot(k(1)+p(1),p(1))*{BUBBLE}"),
            &format!("(dot(k(1),p(1))+dot(p(1),p(1)))*{BUBBLE}"),
        );
        // Another external vector contracted with a loop momentum is a tensor.
        assert_p_error(
            &format!("dot(k(1)+p(2),p(2))*{BUBBLE}"),
            "a loop momentum is contracted with p(2)",
        );
    }
}
