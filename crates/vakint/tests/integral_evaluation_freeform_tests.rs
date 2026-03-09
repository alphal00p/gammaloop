mod test_utils;
use vakint::{EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, PySecDecOptions};

use std::vec;

use symbolica::try_parse;
use vakint::{Vakint, VakintSettings};

use crate::test_utils::compare_vakint_evaluation_vs_reference;
use vakint::{externals_from_f64, params_from_f64};

const N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS: u32 = 32;

const N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS: u32 = 10;
// PySecDec QMC is often very optimistic
const MAX_PULL: f64 = 1.0e99;

#[test_log::test]
fn test_integrate_1l_decorated_indices_alphaloop() {
    #[rustfmt::skip]
    Vakint::initialize_vakint_symbols();
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 3, integral_normalization_factor: LoopNormalizationFactor::MSbar, mu_r_sq_symbol: "some_space::{positive,scalar}::DDmursq".into(),..VakintSettings::default()},
        EvaluationOrder::alphaloop_only(),
        try_parse!(
            "( (user_space::DDsigma(user_space::some_args) + user_space::{scalar}::DDsigma2(user_space::some_args) + user_space::{symmetric,real}::DDsigma3(user_space::some_args))*vakint::p(1,user_space::mink4(4,33))*vakint::p(2,user_space::mink4(4,33))*vakint::p(1,user_space::mink4(4,11))*vakint::p(2,user_space::mink4(4,22))+vakint::k(3,user_space::mink4(4,11))*vakint::k(3,user_space::mink4(4,22)) + vakint::k(3,user_space::mink4(4,77))*vakint::p(1,user_space::mink4(4,77)) ) \
             * vakint::topo( vakint::prop(9,vakint::edge(66,66),vakint::k(3),user_space::{real}::DDMUVsq,1 ))\
            ")
        .unwrap()
        .as_view(),
        params_from_f64(&[
            ("user_space::{real}::DDMUVsq".into(), 1.0), ("some_space::{positive,scalar}::DDmursq".into(), 1.0),
            ("vakint::g(user_space::mink4(4,22),user_space::mink4(4,11))".into(), 1.0),
            ("vakint::p(1,user_space::mink4(4,11))".into(), 1.0),
            ("vakint::p(2,user_space::mink4(4,22))".into(), 1.0),
            ("user_space::DDsigma(user_space::some_args)".into(), 0.0),
            ("user_space::{scalar}::DDsigma2(user_space::some_args)".into(), 0.5),
            ("user_space::{symmetric,real}::DDsigma3(user_space::some_args)".into(), 0.5),
            ].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-36.79980696990178527744327373291".into()),),
            (0,  ("0.0".into(), "99.70093613678094097571666372594".into()),),
            (1,  ("0.0".into(), "-183.0878169087836963770976657747".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0
    );
}

#[test_log::test]
fn test_integrate_1l_decorated_indices_matad() {
    #[rustfmt::skip]
    Vakint::initialize_vakint_symbols();
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 3, integral_normalization_factor: LoopNormalizationFactor::MSbar,..VakintSettings::default()},
        EvaluationOrder::matad_only(None),
        try_parse!(
            "( (user_space::CCsigma(user_space::some_args) + user_space::{scalar}::CCsigma2(user_space::some_args) + user_space::{real,symmetric}::CCsigma3(user_space::some_args))*vakint::p(1,user_space::mink4(4,33))*vakint::p(2,user_space::mink4(4,33))*vakint::p(1,user_space::mink4(4,11))*vakint::p(2,user_space::mink4(4,22))+vakint::k(3,user_space::mink4(4,11))*vakint::k(3,user_space::mink4(4,22)) + vakint::k(3,user_space::mink4(4,77))*vakint::p(1,user_space::mink4(4,77)) ) \
             * vakint::topo( vakint::prop(9,vakint::edge(66,66),vakint::k(3),user_space::MUVsq,1 ))\
            ")
        .unwrap()
        .as_view(),
        params_from_f64(&[
            ("user_space::MUVsq".into(), 1.0), ("vakint::mursq".into(), 1.0),
            ("vakint::g(user_space::mink4(4,22),user_space::mink4(4,11))".into(), 1.0),
            ("vakint::p(1,user_space::mink4(4,11))".into(), 1.0),
            ("vakint::p(2,user_space::mink4(4,22))".into(), 1.0),
            ("user_space::CCsigma(user_space::some_args)".into(), 0.5),
            ("user_space::{scalar}::CCsigma2(user_space::some_args)".into(), 0.5),
            ("user_space::{real,symmetric}::CCsigma3(user_space::some_args)".into(), 0.0),
            ].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-36.79980696990178527744327373291".into()),),
            (0,  ("0.0".into(), "99.70093613678094097571666372594".into()),),
            (1,  ("0.0".into(), "-183.0878169087836963770976657747".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0

    );
}

#[test_log::test]
fn test_integrate_1l_decorated_indices_fmft() {
    #[rustfmt::skip]
    Vakint::initialize_vakint_symbols();
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 5, integral_normalization_factor: LoopNormalizationFactor::MSbar, mu_r_sq_symbol: "some_space::{real,scalar}::BBmursq".into(),..VakintSettings::default()},
        EvaluationOrder::fmft_only(None),
        try_parse!(
            "( ( user_space::BBsigma(user_space::some_args) + user_space::{symmetric,scalar}::BBsigma2(user_space::{real}::some_args2) + user_space::{integer}::BBparam )*vakint::p(1,user_space::mink4(4,33))*vakint::p(2,user_space::mink4(4,33))*vakint::p(1,user_space::mink4(4,11))*vakint::p(2,user_space::mink4(4,22))+vakint::k(3,user_space::mink4(4,11))*vakint::k(3,user_space::mink4(4,22)) + vakint::k(3,user_space::mink4(4,77))*vakint::p(1,user_space::mink4(4,77)) ) \
             * vakint::topo( 
                  vakint::prop(9,vakint::edge(66,66),vakint::k(1),user_space::{real}::BBMUVsq,1 )
                * vakint::prop(9,vakint::edge(66,66),vakint::k(2),user_space::{real}::BBMUVsq,1 )
                * vakint::prop(9,vakint::edge(66,66),vakint::k(3),user_space::{real}::BBMUVsq,1 )
                * vakint::prop(9,vakint::edge(66,66),vakint::k(4),user_space::{real}::BBMUVsq,1 )
            )\
            ")
        .unwrap()
        .as_view(),
        params_from_f64(&[
            ("user_space::{real}::BBMUVsq".into(), 1.0), ("some_space::{real,scalar}::BBmursq".into(), 1.0),
            ("vakint::g(user_space::mink4(4,22),user_space::mink4(4,11))".into(), 1.0),
            ("vakint::p(1,user_space::mink4(4,11))".into(), 1.0),
            ("vakint::p(2,user_space::mink4(4,22))".into(), 1.0),
            ("user_space::BBsigma(user_space::some_args)".into(), 1.0),
            ("user_space::{symmetric,scalar}::BBsigma2(user_space::{real}::some_args2)".into(), 0.0),
            ("user_space::{integer}::BBparam".into(), 0.0),
            ].iter().cloned().collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=2)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS),
        vec![
            (-4, ("-35378.93674652074486878056225484".into(), "0.0".into()),),
            (-3,  ("379847.4112339445643392527721787".into(), "0.0".into()),),
            (-2,  ("-2225660.681830551493330196973382".into(), "0.0".into()),),
            (-1,  ("9310322.197027486227836524560952".into(), "0.0".into()),),
            (0,  ("-31010419.86210674125026829677739".into(), "0.0".into()),),
        ],
        N_DIGITS_ANLYTICAL_EVALUATION_FOR_TESTS, 1.0

    );
}

#[test_log::test]
fn test_integrate_1l_decorated_indices_pysecdec() {
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 5, integral_normalization_factor: LoopNormalizationFactor::MSbar, mu_r_sq_symbol: "some_space::{real,scalar}::AAmursq".into(), ..VakintSettings::default()},
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/test_integrate_1l_decorated_indices_pysecdec".into()) ,..PySecDecOptions::default() })]),
        try_parse!(
            "((user_space::{scalar,real}::AAsigma2+user_space::AAsigma)*vakint::k(1,user_space::mink4(4,11))*vakint::p(1,user_space::mink4(4,11))*vakint::k(1,user_space::mink4(4,12))*vakint::p(1,user_space::mink4(4,12)))*vakint::topo(\
                vakint::prop(1,vakint::edge(1,1),vakint::k(1),user_space::{scalar,real}::AAmuvsq,2)\
             )"
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[
            ("user_space::{scalar,real}::AAmuvsq".into(), 1.0), 
            ("some_space::{real,scalar}::AAmursq".into(), 1.0),
            ("user_space::AAsigma".into(), 0.5),
            ("user_space::{scalar,real}::AAsigma2".into(), 0.5)
            ].iter().cloned().collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(),  "-15.41829599538179739876641630989".into()),),
            (0,  ("0.0".into(),   "41.25556923066471696715797280407".into()),),
            (1,  ("0.0".into(),  "-75.58506810083682014451198579381".into()),),
            (2,  ("0.0".into(),   "104.8268973321405776943695037314".into()),),
            (3,  ("0.0".into(),  "-130.2125934660588794479583559489".into()),),
        ],
        N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS, MAX_PULL
    );
}

#[allow(dead_code)]
fn run_integral_evaluation_freeform_tests() {
    test_integrate_1l_decorated_indices_alphaloop();
    test_integrate_1l_decorated_indices_matad();
    test_integrate_1l_decorated_indices_fmft();
    test_integrate_1l_decorated_indices_pysecdec();
}
