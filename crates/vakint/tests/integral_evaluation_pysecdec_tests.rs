mod test_utils;
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, PySecDecOptions, VakintSettings,
    vakint_parse,
};

use std::{collections::HashMap, vec};

use vakint::{externals_from_f64, params_from_f64};

use crate::test_utils::compare_vakint_evaluation_vs_reference;

const N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS: u32 = 10;
// PySecDec QMC is often very optimistic
const MAX_PULL: f64 = 1.0e99;

#[test_log::test]
fn test_integrate_1l_simple() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{ number_of_terms_in_epsilon_expansion: 2, integral_normalization_factor: LoopNormalizationFactor::pySecDec, ..VakintSettings::default()},
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/pysecdec_test_integrate_1l_simple".into()) ,..PySecDecOptions::default() })]),
        vakint_parse!(
            "( 1 )*topo(\
                prop(1,edge(1,1),k(1),muvsq,1)\
             )"
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        HashMap::default(),
        vec![
            (-1, ("1.0".into(), "0.0".into()),),
            (0,  ("4.227843350984671393934879099176e-1".into(),  "0.0".into()),),
        ],
        N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS, MAX_PULL
    );
}

#[test_log::test]
fn test_integrate_1l_cross_product() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 5, integral_normalization_factor: LoopNormalizationFactor::MSbar, ..VakintSettings::default()},
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/test_integrate_1l_cross_product".into()) ,..PySecDecOptions::default() })]),
        vakint_parse!(
            "(k(1,11)*p(1,11)*k(1,12)*p(1,12))*topo(\
                prop(1,edge(1,1),k(1),muvsq,2)\
             )"
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-15.41829599538179739876641630989".into()),),
            (0,  ("0.0".into(),  "41.25556923066471696715797280407".into()),),
            (1,  ("0.0".into(), "-75.58506810083682014451198579381".into()),),
            (2,  ("0.0".into(),  "104.8268973321405776943695037314".into()),),
            (3,  ("0.0".into(),  "-130.2125934660588794479583559489".into()),),
        ],
        N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS, MAX_PULL
    );
}

#[test_log::test]
fn test_integrate_1l_cross_product_with_additional_symbols_numerator() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 5, integral_normalization_factor: LoopNormalizationFactor::MSbar, ..VakintSettings::default()},
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/test_integrate_1l_cross_product_additional_symbols_numerator".into()) ,..PySecDecOptions::default() })]),
        vakint_parse!(
            "(user_space::A*k(1,11)*p(1,11)*k(1,12)*p(1,12)+user_space::B)*topo(\
                prop(1,edge(1,1),k(1),muvsq,2)\
             )"
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0), ("user_space::A".into(), 3.0), ("user_space::B".into(), 4.0)].iter().cloned().collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
            .collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        vec![
            (-1, ("0.0".into(), "-6.776470381787957720961284930165".into()),),
            (0,  ("0.0".into(), "-21.34624897436485396509954629830".into()),),
            (1,  ("0.0".into(), "72.41426780477804749630966479154".into()),),
            (2,  ("0.0".into(), "-147.4626326303287005964359960851".into()),),
            (3,  ("0.0".into(), "211.1788648410998597813244933220".into()),),
        ],
        N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS, MAX_PULL
    );
}

#[test_log::test]
fn test_integrate_2l_different_masses() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar, allow_unknown_integrals: true, ..VakintSettings::default()},
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/test_integrate_2l_different_masses".into()) ,..PySecDecOptions::default() })]),
        vakint_parse!(
            "(1)*topo(\
                prop(1,edge(1,2),k(1),muvsqA,1)\
                *prop(2,edge(1,2),k(2),muvsqB,1)\
                *prop(3,edge(2,1),k(1)+k(2),muvsqC,1)\
            )"
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsqA".into(), 1.0), ("muvsqB".into(), 1.0), ("muvsqC".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        vec![
            (-2, ("-146.1136365510036558546604990331".into(), "0.0".into()),),
            (-1, ("635.8146971740286947808753047759".into(),  "0.0".into()),),
            (0,  ("-1646.531034471454483109201793220".into(), "0.0".into()),),
            (1,  ("2240.516116133454318298096346441".into(),  "0.0".into()),),
        ],
        N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS, MAX_PULL
    );
}

#[test_log::test]
fn test_integrate_3l_o_eps() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5, ..VakintSettings::default()},
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/test_integrate_3l_o_eps".into()) ,..PySecDecOptions::default() })]),
        vakint_parse!(
            "(1)*topo(\
                 prop(1,edge(1,2),k(1),muvsq,1)\
                *prop(2,edge(2,3),k(2),muvsq,1)\
                *prop(3,edge(3,1),k(3),muvsq,1)\
                *prop(4,edge(1,4),k(3)-k(1),muvsq,1)\
                *prop(5,edge(2,4),k(1)-k(2),muvsq,1)\
                *prop(6,edge(3,4),k(2)-k(3),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        externals_from_f64(
        &(1..=1)
            .map(|i| (i, (17.0*((i+1) as f64), 4.0*((i+2) as f64), 3.0*((i+3) as f64), 12.0*((i+4) as f64))))
            .collect(),
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
            vec![
                // These first two entries are obtained from the analytic expression
                (-1, ("0.0".into(), "-2311.289033520460340396770711738".into()),),
                ( 0, ("0.0".into(), "35134.99893627257345553503414002".into()),),
                // This last entry does not have an analytical expression (MATAD not implemented yet)
                ( 1, ("0.0".into(), "-287175.6919264755".into()),),
            ],
        N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS, MAX_PULL
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_h() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    let vakint_default_settings = VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        number_of_terms_in_epsilon_expansion: 5,
        ..VakintSettings::default()
    };
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        vakint_default_settings,
        EvaluationOrder(vec![EvaluationMethod::PySecDec(
            PySecDecOptions{ min_n_evals: 100_000, max_n_evals: 1_000_000, reuse_existing_output: Some("./tests_workspace/test_integrate_4l_h".into()), ..PySecDecOptions::default()} )]),
        vakint_parse!(
            "(1)*topo(\
                 prop(1,edge(5,1),k(1),muvsq,1)\
                *prop(2,edge(2,6),k(2),muvsq,1)\
                *prop(3,edge(6,5),k(3),muvsq,1)\
                *prop(4,edge(3,4),k(4),muvsq,1)\
                *prop(5,edge(4,5),k(1)-k(3),muvsq,1)\
                *prop(6,edge(6,3),k(2)-k(3),muvsq,1)\
                *prop(7,edge(4,1),k(3)-k(1)+k(4),muvsq,1)\
                *prop(8,edge(2,3),k(3)-k(2)+k(4),muvsq,1)\
                *prop(9,edge(1,2),k(3)+k(4),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            4),
        HashMap::default(),
        vec![
            (0,  ("-12799.53514305961548130719263292".into(), "0.0".into()),),
        ],
        4, MAX_PULL
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_PR9d_from_FG_pinch() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    let vakint_default_settings = VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        number_of_terms_in_epsilon_expansion: 6,
        ..VakintSettings::default()
    };
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        vakint_default_settings,
        EvaluationOrder(vec![EvaluationMethod::PySecDec(
            PySecDecOptions{ min_n_evals: 100_000, max_n_evals: 1_000_000, reuse_existing_output: Some("./tests_workspace/test_integrate_4l_PR9d_from_FG_pinch".into()), ..PySecDecOptions::default()} )]),
        vakint_parse!(
                "( 1 )*topo(\
                    prop(1,edge(5,3),k(1),muvsq,1)\
                    *prop(2,edge(3,4),k(2),muvsq,2)\
                    *prop(3,edge(4,5),k(3),muvsq,1)\
                    *prop(5,edge(5,1),k(4),muvsq,1)\
                    *prop(6,edge(4,1),k(2)-k(3),muvsq,1)\
                    *prop(7,edge(1,5),k(1)-k(3)+k(4),muvsq,1)\
                    *prop(8,edge(3,1),k(1)-k(2),muvsq,1)\
                )"
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            4),
        HashMap::default(),
        vec![
            (-4,  ("8.333333333333333333333333333333e-2".into(), "0.0".into()),),
            (-3,  ("3.333333333333333333333333333333e-1".into(), "0.0".into()),),
            (-2,  ("-3.144646082033583725553786166618e-1".into(), "0.0".into()),),
            (-1,  ("5.421352941798334340259377610275".into(), "0.0".into()),),
            ( 0,  ("-28.31064373017674207211847384976".into(), "0.0".into()),),
            // No analytical result for the O(ep) term
            ( 1,  ("154.6355".into(), "0.0".into()),),
        ],
        4, MAX_PULL
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_PR11d() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    let vakint_default_settings = VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        number_of_terms_in_epsilon_expansion: 5,
        ..VakintSettings::default()
    };
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        vakint_default_settings,
        EvaluationOrder(vec![EvaluationMethod::PySecDec(
            PySecDecOptions{ min_n_evals: 100_000, max_n_evals: 1_000_000, reuse_existing_output: Some("./tests_workspace/test_integrate_4l_PR11d".into()), ..PySecDecOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(\
                 prop(1,edge(1,2),k(1),muvsq,2)\
                *prop(2,edge(2,5),k(2),muvsq,1)\
                *prop(3,edge(3,4),k(3),muvsq,1)\
                *prop(4,edge(4,5),k(4),muvsq,1)\
                *prop(5,edge(2,3),k(1)-k(2),muvsq,1)\
                *prop(6,edge(4,1),k(3)-k(4),muvsq,1)\
                *prop(7,edge(5,3),k(2)+k(3)-k(1),muvsq,1)\
                *prop(8,edge(1,5),k(3)-k(4)-k(1),muvsq,1)\
            )"
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            4),
        HashMap::default(),
        vec![
            (0,  ("-2.906486288643112641819206002127".into(), "0.0".into()),),
        ],
        4, MAX_PULL
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_clover() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    let vakint_default_settings = VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        number_of_terms_in_epsilon_expansion: 5,
        ..VakintSettings::default()
    };
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        vakint_default_settings,
        EvaluationOrder(vec![EvaluationMethod::PySecDec(
            PySecDecOptions{ relative_precision: 1e-8, min_n_evals: 10_000_000, max_n_evals: 100_000_000, reuse_existing_output: Some("./tests_workspace/test_integrate_4l_clover".into()), ..PySecDecOptions::default()} )]),
        vakint_parse!(
            "( 1 )*topo(
              prop(1, edge(1, 1), k(1), muvsq, 1)*\
              prop(2, edge(1, 1), k(2), muvsq, 1)*\
              prop(3, edge(1, 1), k(3), muvsq, 1)*\
              prop(4, edge(1, 1), k(4), muvsq, 1)
          )"
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 1.0), ("mursq".into(), 1.0)].iter().cloned().collect(),
            4),
        HashMap::default(),
        vec![
            (-4,  ("1.000000000000000000000000000000".into(), "0.0".into()),),
            (-3,  ("4.000000000000000000000000000000".into(), "0.0".into()),),
            (-2,  ("13.28986813369645287294483033329".into(), "0.0".into()),),
            (-1,  ("31.55672999723968577791300378449".into(), "0.0".into()),),
            (0,   ("67.98165058904685502307905531744".into(), "0.0".into()),),
        ],
        4, MAX_PULL
    );
}

#[allow(non_snake_case)]
#[test_log::test]
fn test_integrate_4l_clover_with_numerator() {
    if test_utils::should_skip_pysecdec_tests() {
        return;
    }
    let vakint_default_settings = VakintSettings {
        integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
        number_of_terms_in_epsilon_expansion: 5,
        ..VakintSettings::default()
    };
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        vakint_default_settings,
        EvaluationOrder(vec![EvaluationMethod::PySecDec(
            PySecDecOptions{ relative_precision: 1e-8, min_n_evals: 10_000_000, max_n_evals: 100_000_000, reuse_existing_output: Some("./tests_workspace/test_integrate_4l_clover_with_numerator".into()), ..PySecDecOptions::default()} )]),
        vakint_parse!(
            "(
                user_space::A * k(1,11)*k(2,11)*k(1,22)*k(2,22)
              + user_space::B * p(1,11)*k(3,11)*k(3,22)*p(2,22)
              + user_space::C * p(1,11)*p(2,11)*(k(2,22)+k(1,22))*k(2,22)
           )*topo(
              prop(1, edge(1, 1), k(1), muvsq, 2)*\
              prop(2, edge(1, 1), k(2), muvsq, 1)*\
              prop(3, edge(1, 1), k(3), muvsq, 1)*\
              prop(4, edge(1, 1), k(4), muvsq, 1)
          )"
        )
        .unwrap()
        .as_view(),
        // Masses chosen equal on purpose here so as to have a reliable target analytical result
        params_from_f64(&[("muvsq".into(), 0.3), ("mursq".into(), 0.7), ("user_space::A".into(), 3.0), ("user_space::B".into(), 4.0), ("user_space::C".into(), 5.0)].iter().cloned().collect(),
            4),
        externals_from_f64(
            &(1..=2)
                .map(|i| (i, (0.17*((i+1) as f64), 0.4*((i+2) as f64), 0.3*((i+3) as f64), 0.12*((i+4) as f64))))
                .collect(),
                N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS),
        vec![
            (-4,  ("-1.897149599999999855007182247846e-1".into(), "0.0".into()),),
            (-3,  ("-1.495259819655131009380566817668".into(), "0.0".into()),),
            (-2,  ("-6.805240907875078933713325181389".into(), "0.0".into()),),
            (-1,  ("-22.56027900679456203938234477552".into(), "0.0".into()),),
            (0,   ("-60.49337040949871593265194938449".into(), "0.0".into()),),
        ],
        4, MAX_PULL
    );
}

#[allow(dead_code)]
fn run_integral_evaluation_pysecdec_tests() {
    // Convenience runner to execute all tests in this module.
    test_integrate_1l_simple();
    test_integrate_1l_cross_product();
    test_integrate_1l_cross_product_with_additional_symbols_numerator();
    test_integrate_2l_different_masses();
    test_integrate_3l_o_eps();
    test_integrate_4l_h();
    test_integrate_4l_PR9d_from_FG_pinch();
    test_integrate_4l_PR11d();
    test_integrate_4l_clover();
    test_integrate_4l_clover_with_numerator();
}
