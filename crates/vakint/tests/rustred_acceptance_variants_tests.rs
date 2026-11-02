//! Native peers for remaining legacy settings/input variants.
//! Existing tests are unchanged. Every three-loop peer remains explicitly
//! pending until a certified reusable K6 artifact is installed.

mod test_utils;

use test_utils::{
    External, ReferencePoint, RustRedParityPolicy, TensorPrepass, alphaloop_matad_rustred_lanes,
    analytic_matad_rustred_lanes, compare_exact_peers,
};
use vakint::{LoopNormalizationFactor, Vakint, VakintSettings, vakint_parse};

// Literal legacy input data shared only within this acceptance module.
const MERCEDES: &str = "topo(
    prop(1,edge(1,2),k(1),muvsq,1)*prop(2,edge(2,3),k(2),muvsq,1)*
    prop(3,edge(3,1),k(3),muvsq,1)*prop(4,edge(1,4),k(3)-k(1),muvsq,1)*
    prop(5,edge(2,4),k(1)-k(2),muvsq,1)*prop(6,edge(3,4),k(2)-k(3),muvsq,1))";
const RANK_FOUR: &str = "k(1,11)*k(2,11)*k(1,22)*k(2,22)
    +p(1,11)*k(3,11)*k(3,22)*p(2,22)
    +p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)";
const BASKETBALL: &str = "topo(
    prop(1,edge(2,1),k(1),muvsq,1)*prop(2,edge(2,1),k(2),muvsq,1)*
    prop(3,edge(1,2),k(3),muvsq,1)*prop(4,edge(1,2),k(1)+k(2)-k(3),muvsq,1))";
const SMALL_EXTERNALS: [External; 2] = [
    (1, (0.17 * 2.0, 0.4 * 3.0, 0.3 * 4.0, 0.12 * 5.0)),
    (2, (0.17 * 3.0, 0.4 * 4.0, 0.3 * 5.0, 0.12 * 6.0)),
];
const LARGE_EXTERNALS: [External; 2] = [
    (1, (17.0 * 2.0, 4.0 * 3.0, 3.0 * 4.0, 12.0 * 5.0)),
    (2, (17.0 * 3.0, 4.0 * 4.0, 3.0 * 5.0, 12.0 * 6.0)),
];

/// Unit mass at the unit and at the doubled renormalization scale.
const UNIT_SCALE: [(&str, f64); 2] = [("muvsq", 1.0), ("mursq", 1.0)];
const DOUBLED_SCALE: [(&str, f64); 2] = [("muvsq", 1.0), ("mursq", 2.0)];

/// `n` orders of the ε expansion, everything else at its default.
fn epsilon_terms(n: i64) -> VakintSettings {
    VakintSettings {
        number_of_terms_in_epsilon_expansion: n,
        ..VakintSettings::default()
    }
}

#[test]
fn rustred_decorated_one_loop_five_terms() {
    // Original: integral_evaluation_freeform_tests::
    // test_integrate_1l_decorated_indices_pysecdec. Retain decorated attributes,
    // custom scale name and epsilon^3 depth, but compare exact native peers
    // instead of reusing the original low-precision Monte Carlo references.
    Vakint::initialize_vakint_symbols();
    let input = vakint_parse!(
        "((user_space::{scalar,real}::AAsigma2+user_space::AAsigma)*
        vakint::k(1,user_space::mink4(4,11))*vakint::p(1,user_space::mink4(4,11))*
        vakint::k(1,user_space::mink4(4,12))*vakint::p(1,user_space::mink4(4,12)))*
        vakint::topo(vakint::prop(1,vakint::edge(1,1),vakint::k(1),
        user_space::{scalar,real}::AAmuvsq,2))"
    )
    .unwrap();
    compare_exact_peers(
        VakintSettings {
            number_of_terms_in_epsilon_expansion: 5,
            integral_normalization_factor: LoopNormalizationFactor::MSbar,
            mu_r_sq_symbol: "some_space::{real,scalar}::AAmursq".into(),
            ..VakintSettings::default()
        },
        &analytic_matad_rustred_lanes(RustRedParityPolicy::ExactMatadBasis),
        TensorPrepass::FeynKit,
        input.as_view(),
        &ReferencePoint {
            params: &[
                ("user_space::{scalar,real}::AAmuvsq", 1.0),
                ("some_space::{real,scalar}::AAmursq", 1.0),
                ("user_space::AAsigma", 0.5),
                ("user_space::{scalar,real}::AAsigma2", 0.5),
            ],
            externals: &SMALL_EXTERNALS[..1],
        },
        true,
    );
}

#[test]
fn rustred_three_loop_scalar_nonunit_scale() {
    // integral_comparison_vs_pysecdec_tests::test_integrate_3l_pysecdec.
    compare_exact_peers(
        VakintSettings::default(),
        &alphaloop_matad_rustred_lanes(RustRedParityPolicy::NumericalOnly),
        TensorPrepass::Skip,
        vakint_parse!(MERCEDES).unwrap().as_view(),
        &ReferencePoint {
            params: &DOUBLED_SCALE,
            externals: &LARGE_EXTERNALS,
        },
        true,
    );
}

#[test]
fn rustred_three_loop_rank_four_large_external_vectors() {
    // integral_comparison_vs_pysecdec_tests::test_integrate_3l_rank_4.
    compare_exact_peers(
        VakintSettings::default(),
        &alphaloop_matad_rustred_lanes(RustRedParityPolicy::NumericalOnly),
        TensorPrepass::FeynKit,
        vakint_parse!(&format!("({RANK_FOUR})*{MERCEDES}"))
            .unwrap()
            .as_view(),
        &ReferencePoint {
            params: &DOUBLED_SCALE,
            externals: &LARGE_EXTERNALS,
        },
        true,
    );
}

#[test]
fn rustred_three_loop_rank_four_five_terms() {
    // integral_comparison_vs_pysecdec_tests::test_integrate_3l_rank_4_matad.
    compare_exact_peers(
        epsilon_terms(5),
        &analytic_matad_rustred_lanes(RustRedParityPolicy::NumericalOnly),
        TensorPrepass::FeynKit,
        vakint_parse!(&format!("({RANK_FOUR})*{MERCEDES}"))
            .unwrap()
            .as_view(),
        &ReferencePoint {
            params: &DOUBLED_SCALE,
            externals: &SMALL_EXTERNALS,
        },
        true,
    );
}

#[test]
fn rustred_basketball_a_finite_part() {
    // Preserve the original pole-only test; this companion adds epsilon^0.
    compare_exact_peers(
        epsilon_terms(4),
        &alphaloop_matad_rustred_lanes(RustRedParityPolicy::NumericalOnly),
        TensorPrepass::FeynKit,
        vakint_parse!(&format!("k(1,1)^2*{BASKETBALL}"))
            .unwrap()
            .as_view(),
        &ReferencePoint {
            params: &UNIT_SCALE,
            externals: &SMALL_EXTERNALS,
        },
        true,
    );
}

#[test]
fn rustred_basketball_b_finite_part() {
    compare_exact_peers(
        epsilon_terms(4),
        &alphaloop_matad_rustred_lanes(RustRedParityPolicy::NumericalOnly),
        TensorPrepass::FeynKit,
        vakint_parse!(&format!("k(3,1)^2*{BASKETBALL}"))
            .unwrap()
            .as_view(),
        &ReferencePoint {
            params: &UNIT_SCALE,
            externals: &SMALL_EXTERNALS,
        },
        true,
    );
}
