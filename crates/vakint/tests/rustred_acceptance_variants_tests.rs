//! Native peers for remaining legacy settings/input variants; see
//! `RUSTRED_ACCEPTANCE.md` for the complete mapping and workload aliases.
//! Existing tests are unchanged. Every three-loop peer remains explicitly
//! pending until a certified reusable K6 artifact is installed.

mod test_utils;

use test_utils::{
    RustRedParityPolicy, TensorPrepass, alphaloop_matad_rustred_lanes,
    analytic_matad_rustred_lanes, compare_evaluations,
};
use vakint::{
    LoopNormalizationFactor, Vakint, VakintSettings, externals_from_f64, params_from_f64,
    vakint_parse,
};

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
const SMALL_EXTERNALS: [(usize, (f64, f64, f64, f64)); 2] = [
    (1, (0.17 * 2.0, 0.4 * 3.0, 0.3 * 4.0, 0.12 * 5.0)),
    (2, (0.17 * 3.0, 0.4 * 4.0, 0.3 * 5.0, 0.12 * 6.0)),
];
const LARGE_EXTERNALS: [(usize, (f64, f64, f64, f64)); 2] = [
    (1, (17.0 * 2.0, 4.0 * 3.0, 3.0 * 4.0, 12.0 * 5.0)),
    (2, (17.0 * 3.0, 4.0 * 4.0, 3.0 * 5.0, 12.0 * 6.0)),
];

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
    compare_evaluations(
        VakintSettings {
            number_of_terms_in_epsilon_expansion: 5,
            integral_normalization_factor: LoopNormalizationFactor::MSbar,
            mu_r_sq_symbol: "some_space::{real,scalar}::AAmursq".into(),
            run_time_decimal_precision: 32,
            ..VakintSettings::default()
        },
        &analytic_matad_rustred_lanes(RustRedParityPolicy::ExactMatadBasis),
        TensorPrepass::FeynKit,
        input.as_view(),
        params_from_f64(
            &[
                ("user_space::{scalar,real}::AAmuvsq".into(), 1.0),
                ("some_space::{real,scalar}::AAmursq".into(), 1.0),
                ("user_space::AAsigma".into(), 0.5),
                ("user_space::{scalar,real}::AAsigma2".into(), 0.5),
            ]
            .into_iter()
            .collect(),
            32,
        ),
        externals_from_f64(&SMALL_EXTERNALS[..1].iter().copied().collect(), 32),
        1.0e-25,
        0.0,
        true,
    );
}

#[test]
fn rustred_three_loop_scalar_nonunit_scale() {
    // integral_comparison_vs_pysecdec_tests::test_integrate_3l_pysecdec.
    let input = vakint_parse!(MERCEDES).unwrap();
    compare_evaluations(
        VakintSettings {
            run_time_decimal_precision: 32,
            ..VakintSettings::default()
        },
        &alphaloop_matad_rustred_lanes(RustRedParityPolicy::NumericalOnly),
        TensorPrepass::Skip,
        input.as_view(),
        params_from_f64(
            &[("muvsq".into(), 1.0), ("mursq".into(), 2.0)]
                .into_iter()
                .collect(),
            32,
        ),
        externals_from_f64(&LARGE_EXTERNALS.into_iter().collect(), 32),
        1.0e-25,
        0.0,
        true,
    );
}

#[test]
fn rustred_three_loop_rank_four_large_external_vectors() {
    // integral_comparison_vs_pysecdec_tests::test_integrate_3l_rank_4.
    let input = vakint_parse!(&format!("({RANK_FOUR})*{MERCEDES}")).unwrap();
    compare_evaluations(
        VakintSettings {
            run_time_decimal_precision: 32,
            ..VakintSettings::default()
        },
        &alphaloop_matad_rustred_lanes(RustRedParityPolicy::NumericalOnly),
        TensorPrepass::FeynKit,
        input.as_view(),
        params_from_f64(
            &[("muvsq".into(), 1.0), ("mursq".into(), 2.0)]
                .into_iter()
                .collect(),
            32,
        ),
        externals_from_f64(&LARGE_EXTERNALS.into_iter().collect(), 32),
        1.0e-25,
        0.0,
        true,
    );
}

#[test]
fn rustred_three_loop_rank_four_five_terms() {
    // integral_comparison_vs_pysecdec_tests::test_integrate_3l_rank_4_matad.
    let input = vakint_parse!(&format!("({RANK_FOUR})*{MERCEDES}")).unwrap();
    compare_evaluations(
        VakintSettings {
            number_of_terms_in_epsilon_expansion: 5,
            run_time_decimal_precision: 32,
            ..VakintSettings::default()
        },
        &analytic_matad_rustred_lanes(RustRedParityPolicy::NumericalOnly),
        TensorPrepass::FeynKit,
        input.as_view(),
        params_from_f64(
            &[("muvsq".into(), 1.0), ("mursq".into(), 2.0)]
                .into_iter()
                .collect(),
            32,
        ),
        externals_from_f64(&SMALL_EXTERNALS.into_iter().collect(), 32),
        1.0e-25,
        0.0,
        true,
    );
}

#[test]
fn rustred_basketball_a_finite_part() {
    // Preserve the original pole-only test; this companion adds epsilon^0.
    let input = vakint_parse!(&format!("k(1,1)^2*{BASKETBALL}")).unwrap();
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
            &[("muvsq".into(), 1.0), ("mursq".into(), 1.0)]
                .into_iter()
                .collect(),
            32,
        ),
        externals_from_f64(&SMALL_EXTERNALS.into_iter().collect(), 32),
        1.0e-25,
        0.0,
        true,
    );
}

#[test]
fn rustred_basketball_b_finite_part() {
    let input = vakint_parse!(&format!("k(3,1)^2*{BASKETBALL}")).unwrap();
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
            &[("muvsq".into(), 1.0), ("mursq".into(), 1.0)]
                .into_iter()
                .collect(),
            32,
        ),
        externals_from_f64(&SMALL_EXTERNALS.into_iter().collect(), 32),
        1.0e-25,
        0.0,
        true,
    );
}
