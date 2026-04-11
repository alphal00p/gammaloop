mod test_utils;
#[allow(unused)]
use vakint::{
    AlphaLoopOptions, EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, MATADOptions,
    PySecDecOptions,
};

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
            (-1, ("0.0".into(), "-2.361163533305129016529279780695e-2".into()),),
            (0,  ("0.0".into(), "-2.282006358584552632588748950063e-2".into()),),
            (1,  ("0.0".into(), "-4.184406937890931251417541408739e-2".into()),),
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
            (-1, ("0.0".into(), "-2.361163533305129016529279780695e-2".into()),),
            (0,  ("0.0".into(), "-2.282006358584552632588748950063e-2".into()),),
            (1,  ("0.0".into(), "-4.184406937890931251417541408739e-2".into()),),
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
            (-4, ( "-5.996072606189216818994270605252e-9".into(), "0.0".into()),),
            (-3,  ("-2.378327420532500222026013157094e-8".into(), "0.0".into()),),
            (-2,  ("-7.878244126888092041321958350898e-8".into(), "0.0".into()),),
            (-1,  ("-1.860926787346530261093067881728e-7".into(), "0.0".into()),),
            (0,  ( "-3.997176154874809449646384405081e-7".into(), "0.0".into()),),
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
        //EvaluationOrder(vec![EvaluationMethod::MATAD(MATADOptions::default())]),
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
            (-1, ("0.0".into(),  "-9.892747066260199063e-3".into()),),
            (0,  ("0.0".into(),  "-9.892747066260199063e-3".into()),),
            (1,  ("0.0".into(),  "-1.802920539739716332e-2".into()),),
            (2,  ("0.0".into(),  "-1.406532376313407440e-2".into()),),
            (3,  ("0.0".into(),  "-2.008809563540125964e-2".into()),),
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
