mod test_utils;

use std::collections::HashMap;
use test_utils::{compare_vakint_evaluation_vs_reference, should_skip_pysecdec_tests};
use vakint::{
    EvaluationOrder, FMFTOptions, LoopNormalizationFactor, MATADOptions, PySecDecOptions, Vakint,
    VakintError, VakintSettings, params_from_f64, vakint_parse,
};

#[test_log::test]
fn massive_tadpoles_in_various_dimensions() {
    if should_skip_pysecdec_tests() {
        return;
    }

    // With the mostly-minus metric and pySecDec normalization,
    // J_n = (-1)^n Gamma(n-D/2) (m^2)^(D/2-n) / Gamma(n).
    // At D=3 and m^2=4, J_1=4 sqrt(pi) and J_2=sqrt(pi)/2, with no epsilon pole.
    // Their logarithmic derivatives are 2-gamma-4 log(2) and -gamma-4 log(2).
    // The first case also checks that the input regulator keeps its meaning.
    for (dimension, power, numerator, normalization, imaginary, coefficients) in [
        (
            3,
            1,
            "1+delta",
            LoopNormalizationFactor::pySecDec,
            false,
            [
                "7.08981540362206410919266993336458073119",
                "-2.48004853220906458609302373132813387191",
            ],
        ),
        (
            3,
            2,
            "1",
            LoopNormalizationFactor::pySecDec,
            false,
            [
                "0.886226925452758013649083741670572591399",
                "-2.96868684288440711420887919142773450819",
            ],
        ),
        // For MSbar at mu^2=9: J_2=i/(16 pi) [1+eps (log(9)-4 log(2))].
        (
            3,
            2,
            "1",
            LoopNormalizationFactor::MSbar,
            true,
            [
                "0.0198943678864869169711104704215642952543",
                "-0.0114465059674054260994301977646579112840",
            ],
        ),
        // At D=0, J_2=1/16 with logarithmic derivative 1-gamma-log(4).
        (
            0,
            2,
            "1",
            LoopNormalizationFactor::pySecDec,
            false,
            ["0.0625", "-0.06021937662633896746506102081242222294957"],
        ),
        // At D=-2, J_2=1/32 with logarithmic derivative 3/2-gamma-log(4).
        (
            -2,
            2,
            "1",
            LoopNormalizationFactor::pySecDec,
            false,
            ["0.03125", "-0.01448468831316948373253051040621111147479"],
        ),
    ] {
        let integral = vakint_parse!(
            format!("({numerator})*topo(prop(1,edge(1,1),k(1),muvsq,{power}))").as_str()
        )
        .unwrap();
        let mut reference = vec![(-1, ("0".into(), "0".into()))];
        reference.extend(coefficients.into_iter().enumerate().map(|(order, value)| {
            (
                order as i64,
                if imaginary {
                    ("0".into(), value.into())
                } else {
                    (value.into(), "0".into())
                },
            )
        }));

        compare_vakint_evaluation_vs_reference(
            VakintSettings {
                dimension,
                epsilon_symbol: "vakint::delta".into(),
                number_of_terms_in_epsilon_expansion: 3,
                integral_normalization_factor: normalization,
                ..VakintSettings::default()
            },
            EvaluationOrder::pysecdec_only(Some(PySecDecOptions {
                max_n_evals: 100_000,
                ..PySecDecOptions::default()
            })),
            integral.as_view(),
            params_from_f64(
                &[("muvsq".into(), 4.0), ("mursq".into(), 9.0)]
                    .into_iter()
                    .collect(),
                10,
            ),
            HashMap::default(),
            reference,
            10,
            1.0e99,
        );
    }
}

#[test_log::test]
fn three_dimensional_msbar_matches_custom_measure() {
    let settings = VakintSettings {
        dimension: 3,
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        ..VakintSettings::default()
    };
    // Compare the built-in normalization against the unsimplified measure:
    // (2 pi)^(-D) [mu^2 exp(gamma)/(4 pi)]^eps per loop.
    let custom = LoopNormalizationFactor::Custom(
        "((2*pi)^(-(3-2*eps))*(exp(log_mu_sq)/(4*pi*exp(-EulerGamma)))^eps)^n_loops".into(),
    );
    let (_, _, msbar) = settings
        .integral_normalization_factor
        .validate(&settings)
        .unwrap();
    let (_, _, explicit_measure) = custom.validate(&settings).unwrap();
    let (matches, message) = msbar.does_approx_match(&explicit_measure, None, 1.0e-25, 0.0);
    assert!(matches, "{message}");
}

#[test_log::test]
fn unsupported_dimensions_fail_before_running_backends() {
    let vakint = Vakint::new().unwrap();
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder::analytic_only(),
        form_exe_path: "/vakint-test-missing-form".into(),
        python_exe_path: "/vakint-test-missing-python".into(),
        ..VakintSettings::default()
    };
    let integral = vakint_parse!("topo(prop(1,edge(1,1),k(1),muvsq,1))").unwrap();
    let specs = vakint
        .topologies
        .match_topologies_to_user_input(integral.as_view(), false)
        .unwrap()
        .unwrap();
    let numerator = vakint_parse!("1").unwrap();
    for dimension in [3, 0, -2] {
        settings.dimension = dimension;
        assert!(matches!(vakint.validate_settings(&settings),
            Err(VakintError::EvaluationError(message)) if message.contains("requires pySecDec")));
        for method in &settings.evaluation_order.0 {
            assert!(matches!(
                method.evaluate_integral(&vakint, &settings, numerator.as_view(), &specs),
                Err(VakintError::NoEvaluationMethodFound(_, _, d)) if d == dimension
            ));
        }
        assert!(
            matches!(vakint.matad_evaluate(&settings, numerator.as_view(), &specs, &MATADOptions::default()),
            Err(VakintError::MATADError(message)) if message.contains("d=4"))
        );
        assert!(
            matches!(vakint.fmft_evaluate(&settings, numerator.as_view(), &specs, &FMFTOptions::default()),
            Err(VakintError::FMFTError(message)) if message.contains("d=4"))
        );
    }
}

#[test_log::test]
fn three_dimensional_scalar_contractions_need_no_tensor_backend() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        dimension: 3,
        use_dot_product_notation: true,
        form_exe_path: "/vakint-test-missing-form".into(),
        ..VakintSettings::default()
    };
    for (numerator, expected) in [
        ("1+ε", "1+ε"),
        ("ε*k(1,1)^2", "ε*dot(k(1),k(1))"),
        ("ε*dot(k(1),k(1))", "ε*dot(k(1),k(1))"),
    ] {
        let input = vakint_parse!(format!("({numerator})*topo(I1L(muvsq,1))").as_str()).unwrap();
        let expected = vakint_parse!(format!("({expected})*topo(I1L(muvsq,1))").as_str()).unwrap();
        assert_eq!(
            vakint.tensor_reduce(&settings, input.as_view()).unwrap(),
            expected
        );
    }
    for numerator in ["k(1,1)*k(1,2)", "g(1,1)"] {
        let input = vakint_parse!(format!("({numerator})*topo(I1L(muvsq,1))").as_str()).unwrap();
        assert!(matches!(vakint.tensor_reduce(&settings, input.as_view()),
            Err(VakintError::InvalidNumerator(message)) if message.contains("Tensor reduction")));
    }
}

#[test_log::test]
fn four_dimensional_normalizations_are_unchanged() {
    let settings = VakintSettings::default();
    assert_eq!(settings.dimension, 4);
    for (normalization, expected) in [
        (
            LoopNormalizationFactor::pySecDec,
            "(𝑖*𝜋^((4-2*ε)/2))^(-n_loops)",
        ),
        (
            LoopNormalizationFactor::FMFTandMATAD,
            "(𝑖*𝜋^((4-2*ε)/2)*exp(-EulerGamma)^ε)^(-n_loops)",
        ),
        (
            LoopNormalizationFactor::MSbar,
            "(2*𝜋)^(-4*n_loops)*exp(log_mu_sq)^(ε*n_loops)*𝜋^(ε*n_loops)*exp(EulerGamma)^(ε*n_loops)",
        ),
        (
            LoopNormalizationFactor::Custom("(1+eps)^n_loops".into()),
            "(1+ε)^n_loops",
        ),
    ] {
        assert_eq!(
            normalization.to_atom(&settings).unwrap(),
            vakint_parse!(expected).unwrap()
        );
    }
}

#[test_log::test]
fn pysecdec_reuse_rejects_missing_or_mismatched_dimension() {
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        dimension: -2,
        python_exe_path: "/vakint-test-missing-python".into(),
        ..VakintSettings::default()
    };
    let directory =
        std::env::temp_dir().join(format!("vakint_dimension_reuse_{}", std::process::id()));
    std::fs::create_dir(&directory).unwrap();
    for metadata in [None, Some("4")] {
        if let Some(dimension) = metadata {
            std::fs::write(directory.join("dimension.txt"), dimension).unwrap();
        }
        assert!(matches!(vakint.run_pysecdec(&settings, &[], vec![], false,
            Some(directory.to_string_lossy().into_owned()), None),
            Err(VakintError::PySecDecError(message)) if message.contains("cached expansion dimension")));
    }
    std::fs::remove_dir_all(directory).unwrap();
}

#[test_log::test]
fn pysecdec_launcher_drains_stderr_and_reuses_matching_dimension() {
    let vakint = Vakint::new().unwrap();
    // The alarm also bounds this test if an unread stderr pipe is reintroduced.
    let script = r#"import signal
import sys
from pathlib import Path
signal.alarm(10)
print("stdout marker", flush=True)
sys.stderr.write("stderr marker\n" + "x" * (128 * 1024))
sys.stderr.flush()
Path("out.txt").write_text("1\n0")
"#;
    for dimension in [3, 0, -2] {
        let settings = VakintSettings {
            dimension,
            ..VakintSettings::default()
        };
        let directory = std::env::temp_dir().join(format!(
            "vakint_launcher_stderr_{}_{dimension}",
            std::process::id()
        ));
        for source in [
            script,
            "raise AssertionError(\"cached sources were overwritten\")",
        ] {
            let result = vakint
                .run_pysecdec(
                    &settings,
                    &[("run.py".into(), source.into())],
                    vec![],
                    false,
                    Some(directory.to_string_lossy().into_owned()),
                    None,
                )
                .unwrap();
            assert_eq!(result, ["1", "0"]);
            let log = std::fs::read_to_string(directory.join("follow_run.txt")).unwrap();
            assert!(log.contains("stdout marker") && log.contains("stderr marker"));
            assert!(log.len() > 128 * 1024);
        }
        std::fs::remove_dir_all(directory).unwrap();
    }
}
