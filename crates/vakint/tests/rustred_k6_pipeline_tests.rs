//! Three-loop tensor/scalar integration groundwork for the shipped K6 program.
//!
//! Native tensor, existing scalar-oracle and RustRed peer checks are enabled
//! for the shipped reusable K6 artifact. A successful prepass alone still does
//! not establish scalar closure; the full native peer is a separate gate.

mod test_utils;

use symbolica::atom::AtomCore;
use test_utils::{
    RustRedParityPolicy, TensorPrepass, alphaloop_matad_lanes, alphaloop_matad_rustred_lanes,
    compare_evaluations,
};
use vakint::{
    TensorReductionMethod, Vakint, VakintSettings, externals_from_f64, params_from_f64,
    vakint_parse,
};

// These are input fixtures, not engine dispatch. Every input goes through
// Vakint's existing matcher, including its routing of the contracted graphs.
const THREE_LOOP_CLASSES: [(&str, &str); 5] = [
    ("I3L(", "topo(I3L(muvsq,1,1,1,1,1,1))"),
    ("I3L_pinch_6(", "topo(I3L_pinch_6(muvsq,1,1,1,1,1,0))"),
    ("I3L_pinch_3_6(", "topo(I3L_pinch_3_6(muvsq,1,1,0,1,1,0))"),
    ("I3L_pinch_1_6(", "topo(I3L_pinch_1_6(muvsq,0,1,1,1,1,0))"),
    (
        "I3L_pinch_1_3_6(",
        "topo(I3L_pinch_1_3_6(muvsq,0,1,0,1,1,0))",
    ),
];

// Retain the mixed-loop tensor and an even-loop rank-four term: the latter
// prevents an odd-loop zero from making a factorized routing check vacuous.
const NUMERATOR: &str = "(muvsq+2*mursq)*(k(1,1)*p(1,1)*k(2,2)*p(2,2)*k(3,3)*p(1,3)*k(1,4)*p(2,4)+k(1,1)*p(1,1)*k(1,2)*p(1,2)*k(2,3)*p(2,3)*k(2,4)*p(2,4))";

#[test]
fn feynkit_three_loop_class_prepasses_match_the_form_oracle() {
    let vakint = Vakint::new().unwrap();
    let native = VakintSettings {
        tensor_reduction_method: TensorReductionMethod::FeynKit,
        form_exe_path: "/this/path/must/not/be/invoked/by-k6-native-prepass".into(),
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    let oracle = VakintSettings {
        tensor_reduction_method: TensorReductionMethod::AlphaLoop,
        use_dot_product_notation: true,
        ..VakintSettings::default()
    };
    for (expected_class, topology) in THREE_LOOP_CLASSES {
        let input = vakint_parse!(&format!("({NUMERATOR})*{topology}")).unwrap();
        let canonical = vakint.to_canonical(&native, input.as_view(), true).unwrap();
        assert!(
            canonical.to_canonical_string().contains(expected_class),
            "fixture was not routed to {expected_class}: {canonical}"
        );
        let native_result = vakint.tensor_reduce(&native, canonical.as_view()).unwrap();
        assert!(
            !native_result.is_zero(),
            "vacuous tensor fixture for {expected_class}"
        );
        let oracle_result = vakint.tensor_reduce(&oracle, canonical.as_view()).unwrap();
        let difference = (Vakint::convert_to_dot_notation(native_result.as_view())
            - Vakint::convert_to_dot_notation(oracle_result.as_view()))
        .expand()
        .together();
        assert!(
            difference.is_zero(),
            "native/FORM tensor prepasses differ for {expected_class}: {difference}"
        );
    }
}

#[test]
fn feynkit_three_loop_class_scalar_oracles_agree() {
    for (expected_class, topology) in THREE_LOOP_CLASSES {
        eprintln!("three-loop tensor/scalar oracle groundwork: {expected_class}");
        let input = vakint_parse!(&format!("({NUMERATOR})*{topology}")).unwrap();
        compare_evaluations(
            VakintSettings {
                number_of_terms_in_epsilon_expansion: 4,
                run_time_decimal_precision: 32,
                ..VakintSettings::default()
            },
            &alphaloop_matad_lanes(),
            TensorPrepass::FeynKit,
            input.as_view(),
            params_from_f64(
                &[("muvsq".into(), 1.3), ("mursq".into(), 0.7)]
                    .into_iter()
                    .collect(),
                32,
            ),
            externals_from_f64(
                &[(1, (0.17, 0.4, 0.3, 0.12)), (2, (0.31, 0.2, 0.7, 0.43))]
                    .into_iter()
                    .collect(),
                32,
            ),
            1.0e-25,
            0.0,
            true,
        );
    }
}

#[test]
fn feynkit_rustred_three_loop_class_numerical_peers() {
    for (expected_class, topology) in THREE_LOOP_CLASSES {
        eprintln!("required FORM-less three-loop RustRed peer: {expected_class}");
        let input = vakint_parse!(&format!("({NUMERATOR})*{topology}")).unwrap();
        compare_evaluations(
            VakintSettings {
                number_of_terms_in_epsilon_expansion: 4,
                run_time_decimal_precision: 32,
                ..VakintSettings::default()
            },
            &alphaloop_matad_rustred_lanes(RustRedParityPolicy::NumericalOnly),
            TensorPrepass::FeynKit,
            input.as_view(),
            params_from_f64(
                &[("muvsq".into(), 1.3), ("mursq".into(), 0.7)]
                    .into_iter()
                    .collect(),
                32,
            ),
            externals_from_f64(
                &[(1, (0.17, 0.4, 0.3, 0.12)), (2, (0.31, 0.2, 0.7, 0.43))]
                    .into_iter()
                    .collect(),
                32,
            ),
            1.0e-25,
            0.0,
            true,
        );
    }
}
