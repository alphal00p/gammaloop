mod test_utils;
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, PySecDecOptions, VakintSettings,
    vakint_parse,
};
#[allow(unused)]
use vakint::{FMFTOptions, MATADOptions};

use std::{collections::HashMap, vec};

use symbolica::domains::float::{Complex, RealLike};
use vakint::{Vakint, externals_from_f64, params_from_complex_f64, params_from_f64};

use crate::test_utils::compare_vakint_evaluation_vs_reference;

const N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS: u32 = 10;
// PySecDec QMC is often very optimistic
const MAX_PULL: f64 = 1.0e99;

#[test_log::test]
#[ignore = "manual PySecDec validation"]
fn test_integrate_1l_simple() {
    test_utils::require_pysecdec_tests();
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{ number_of_terms_in_epsilon_expansion: 2, integral_normalization_factor: LoopNormalizationFactor::pySecDec, ..VakintSettings::default()},
        EvaluationOrder::pysecdec_only(None),
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
#[ignore = "manual PySecDec validation"]
fn test_integrate_1l_complex_and_signed_parameters() {
    test_utils::require_pysecdec_tests();
    let vakint = Vakint::new().unwrap();
    let euler_gamma = 0.577_215_664_901_532_9;
    let tadpole_coefficients = [
        1.0,
        1.0 - euler_gamma,
        1.0 - euler_gamma + euler_gamma * euler_gamma / 2.0 + std::f64::consts::PI.powi(2) / 12.0,
    ];
    for (coefficient, normalization, re, im, first_power, project_onto_tensor_integrals) in [
        ("user_space::phase", "1", 1.0, 2.0, -1, true),
        ("user_space::phase", "1", -3.0, 0.0, -1, true),
        ("user_space::phase/vakint::ε", "1", 1.0, 2.0, -2, true),
        ("user_space::phase", "eps^-1", 1.0, 2.0, -2, true),
        ("user_space::phase*vakint::ε^2", "eps^-1", 1.0, 2.0, 0, true),
        ("user_space::phase", "1", 1.0, 2.0, -1, false),
        ("user_space::phase", "1", -3.0, 0.0, -1, false),
        ("user_space::phase", "eps^-1", 1.0, 2.0, -2, false),
        (
            "user_space::phase*vakint::ε^2",
            "eps^-1",
            1.0,
            2.0,
            0,
            false,
        ),
    ] {
        let input = vakint_parse!(&format!(
            "({coefficient})*topo(prop(1,edge(1,1),k(1),muvsq,1))"
        ))
        .unwrap();
        let mut settings = VakintSettings {
            project_onto_tensor_integrals,
            number_of_terms_in_epsilon_expansion: 2,
            integral_normalization_factor: LoopNormalizationFactor::Custom(format!(
                "({normalization})*({})",
                LoopNormalizationFactor::pySecDec.to_expression()
            )),
            evaluation_order: EvaluationOrder::pysecdec_only(Some(PySecDecOptions {
                relative_precision: 1e-7,
                min_n_evals: 10_000,
                max_n_evals: 100_000,
                // Regenerate the package and its parameter declarations at each point.
                reuse_existing_output: None,
                ..PySecDecOptions::default()
            })),
            ..VakintSettings::default()
        };
        let parameters = params_from_complex_f64(
            &HashMap::from_iter([
                ("muvsq".into(), Complex::new(1.0, 0.0)),
                ("mursq".into(), Complex::new(1.0, 0.0)),
                ("user_space::phase".into(), Complex::new(re, im)),
            ]),
            settings.run_time_decimal_precision,
        );
        settings
            .evaluation_order
            .adjust(
                None,
                1e-7,
                &HashMap::default(),
                &parameters,
                &HashMap::default(),
            )
            .unwrap();
        let canonical = vakint
            .to_canonical(&settings, input.as_view(), true)
            .unwrap();
        let integral = vakint
            .evaluate_integral(&settings, canonical.as_view())
            .unwrap();
        let (result, error) = Vakint::full_numerical_evaluation(
            &settings,
            integral.as_view(),
            &HashMap::default(),
            &parameters,
            None,
        )
        .unwrap();
        assert_eq!(result.0.len(), (1 - first_power) as usize);
        // The same massive tadpole as test_integrate_1l_simple has residue one
        // and finite part 1-EulerGamma. Dividing by epsilon also requires its
        // order-epsilon coefficient to recover the requested finite part.
        // The same applies when the pole belongs to the normalization. A
        // quadratic evanescent numerator then requires coefficient order two
        // despite the requested answer ending at the finite term.
        // Check both components, including zero.
        for (power, coefficient) in (first_power..=0).zip(tadpole_coefficients) {
            let actual = result.get_epsilon_coefficient(power);
            for (component, value, expected) in [
                ("real", actual.re.to_f64(), re * coefficient),
                ("imaginary", actual.im.to_f64(), im * coefficient),
            ] {
                assert!(
                    (value - expected).abs() <= 1e-6 * (1.0 + expected.abs()),
                    "phase=({re},{im}), epsilon^{power} {component}: {value} != {expected}"
                );
            }
            if let Some(error) = &error {
                let uncertainty = error.get_epsilon_coefficient(power);
                for (component, value, expected) in [
                    ("real", uncertainty.re.to_f64(), re * coefficient),
                    ("imaginary", uncertainty.im.to_f64(), im * coefficient),
                ] {
                    assert!(
                        value.is_finite() && value.abs() <= 1e-6 * (1.0 + expected.abs()),
                        "phase=({re},{im}), epsilon^{power} {component} uncertainty: {value}"
                    );
                }
            }
        }
    }
}

#[test_log::test]
#[ignore = "manual PySecDec validation"]
fn test_integrate_1l_cross_product() {
    test_utils::require_pysecdec_tests();
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 5, integral_normalization_factor: LoopNormalizationFactor::MSbar, ..VakintSettings::default()},
        EvaluationOrder::pysecdec_only(None),
        //EvaluationOrder(vec![EvaluationMethod::MATAD(MATADOptions::default())]),
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
            (-1, ("0.0".into(), "-9.8927470662601990625262260437e-3".into()),),
            (0,  ("0.0".into(), "-9.8927470662601990625262260437e-3".into()),),
            (1,  ("0.0".into(), "-1.8029205397397163324058055878e-2".into()),),
            (2,  ("0.0".into(), "-1.4065323763134074397385120392e-2".into()),),
            (3,  ("0.0".into(), "-2.0088095635401259642094373703e-2".into()),),
        ],
        N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS, MAX_PULL
    );
}

#[test_log::test]
#[ignore = "manual PySecDec validation"]
fn test_integrate_1l_cross_product_with_additional_symbols_numerator() {
    test_utils::require_pysecdec_tests();
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{number_of_terms_in_epsilon_expansion: 5, integral_normalization_factor: LoopNormalizationFactor::MSbar, ..VakintSettings::default()},
        EvaluationOrder::pysecdec_only(None),
        //EvaluationOrder(vec![EvaluationMethod::MATAD(MATADOptions::default())]),
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
            (-1, ("0.0".into(), "-4.3479452906467486172914505005e-3".into()),),
            (0,  ("0.0".into(), "-2.9678241200599586591124534607e-2".into()),),
            (1,  ("0.0".into(), "-3.3254282865527784451842308044e-2".into()),),
            (2,  ("0.0".into(), "-5.2345456981129245832562446594e-2".into()),),
            (3,  ("0.0".into(), "-4.4843030036645359359681606293e-2".into()),),
        ],
        N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS, MAX_PULL
    );
}

#[test_log::test]
#[ignore = "manual PySecDec validation"]
fn test_integrate_2l_different_masses() {
    test_utils::require_pysecdec_tests();
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings{integral_normalization_factor: LoopNormalizationFactor::MSbar, allow_unknown_integrals: true, ..VakintSettings::default()},
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/test_integrate_2l_different_masses".into()) ,..PySecDecOptions::default() })]),
        //EvaluationOrder(vec![EvaluationMethod::MATAD(MATADOptions::default())]),
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
            (-2, ("-6.0152239775845828262390568852e-5".into(), "0.0".into()),),
            (-1, ("-1.8045671933464291214477270842e-4".into(),  "0.0".into()),),
            (0,  ("-3.7902087660768302157521247864e-4".into(), "0.0".into()),),
            (1,  ("-9.7081416197397629730403423309e-4".into(),  "0.0".into()),),
        ],
        N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS, MAX_PULL
    );
}

#[test_log::test]
#[ignore = "manual PySecDec validation"]
fn test_integrate_3l_o_eps() {
    test_utils::require_pysecdec_tests();
    #[rustfmt::skip]
    compare_vakint_evaluation_vs_reference(
        VakintSettings { integral_normalization_factor: LoopNormalizationFactor::MSbar, number_of_terms_in_epsilon_expansion: 5, ..VakintSettings::default()},
        EvaluationOrder(vec![EvaluationMethod::PySecDec(PySecDecOptions { reuse_existing_output: Some("./tests_workspace/test_integrate_3l_o_eps".into()) ,..PySecDecOptions::default() })]),
        //EvaluationOrder(vec![EvaluationMethod::MATAD(MATADOptions::default())]),
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
                (-1, ("0.0".into(), "-6.1051429657027990274129863910e-7".into()),),
                ( 0, ("0.0".into(), "+2.5484155391724767149965275403e-6".into()),),
                ( 1, ("0.0".into(), "-1.0634407259633205821030975460e-5".into()),),
            ],
            N_DIGITS_PYSECDEC_EVALUATION_FOR_TESTS, MAX_PULL
    );
}

#[allow(non_snake_case)]
#[test_log::test]
#[ignore = "manual PySecDec validation"]
fn test_integrate_4l_h() {
    test_utils::require_pysecdec_tests();
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
        //EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions::default())]),
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
            (0,  ("-2.169283452273432986058475530569e-9".into(), "0.0".into()),),
        ],
        4, MAX_PULL
    );
}

#[allow(non_snake_case)]
#[test_log::test]
#[ignore = "manual PySecDec validation"]
fn test_integrate_4l_PR9d_from_FG_pinch() {
    test_utils::require_pysecdec_tests();
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
#[ignore = "manual PySecDec validation"]
fn test_integrate_4l_PR11d() {
    test_utils::require_pysecdec_tests();
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
#[ignore = "manual PySecDec validation"]
fn test_integrate_4l_clover() {
    test_utils::require_pysecdec_tests();
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
#[ignore = "manual PySecDec validation"]
fn test_integrate_4l_clover_with_numerator() {
    test_utils::require_pysecdec_tests();
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
    test_integrate_1l_complex_and_signed_parameters();
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
