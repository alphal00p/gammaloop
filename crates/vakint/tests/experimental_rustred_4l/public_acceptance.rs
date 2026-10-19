//! End-to-end tests through the public backend, without a test-local reducer.

use super::{acceptance_inputs, retained_parent_descriptor, test_utils};
use test_utils::{EvaluationTestInput, EvaluationTestLane, RustRedParityPolicy, TensorPrepass};
use vakint::{EvaluationMethod, EvaluationOrder, FMFTOptions};

#[path = "external_master_oracle.rs"]
mod external_master_oracle;

/// Frozen independent FORM/external-data reference, with neither an executable
/// oracle nor a FORM tensor prepass at test time. The historical FORM constants
/// have fewer digits; their original reference is still checked separately.
#[test]
#[ignore = "four-loop native numerator acceptance; explicit bounded execution required"]
fn h_rank_four_matches_reference_without_form() {
    let case = acceptance_inputs::cases()
        .into_iter()
        .find(|case| case.name == "test_integrate_4l_h_rank_4")
        .expect("retained H numerator reference");
    let mut settings = case.settings;
    settings.form_exe_path = "/this/path/must/not/be/invoked/by-rustred-acceptance".into();
    settings.tensor_reduction_method = vakint::TensorReductionMethod::FeynKit;
    assert_eq!(settings.run_time_decimal_precision, 32);
    let reference: serde_json::Value = serde_json::from_str(include_str!(
        "../inputs/h_rank_four_external_reference.json"
    ))
    .unwrap();
    assert_eq!(reference["case"].as_str(), Some(case.name));
    use sha2::Digest;
    assert_eq!(
        format!(
            "{:x}",
            sha2::Sha256::digest(include_bytes!("../inputs/h_rank_four_external_form.txt"))
        ),
        reference["raw_form_sha256"].as_str().unwrap()
    );
    let expected = serde_json::from_value(reference["coefficients"].clone()).unwrap();
    let masses = vakint::params_from_f64(
        &case
            .parameters
            .iter()
            .map(|(name, value)| ((*name).into(), *value))
            .collect(),
        settings.run_time_decimal_precision,
    );
    let external = vakint::externals_from_f64(
        &case.external_momenta.iter().copied().collect(),
        settings.run_time_decimal_precision,
    );
    test_utils::compare_vakint_evaluations_vs_reference(
        settings.clone(),
        &[EvaluationTestLane::rustred_scalar(
            EvaluationTestInput::TensorReduced,
            RustRedParityPolicy::NumericalOnly,
        )],
        TensorPrepass::FeynKit,
        case.input.as_view(),
        masses,
        external,
        expected,
        settings.run_time_decimal_precision,
        1.0,
    );
}

#[test]
#[ignore = "four-loop public RustRed/FMFT numerical inventory; explicit oracle executable required"]
fn all_fifteen_numerical_references_match_fmft_and_public_rustred() {
    let form = std::env::var("VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH").unwrap();
    let filter = std::env::var("VAKINT_4L_CANDIDATE_FAMILY_FILTER").ok();
    assert!(
        filter
            .as_deref()
            .is_none_or(|f| matches!(f, "H" | "FG" | "BMW" | "X"))
    );
    // Keep FORM's historical frozen references and numerical threshold. Some
    // external masters have more digits than the original FORM constants;
    // compare RustRed against FORM's independent exact reduction using that
    // same external scientific data, still at the original 1e-30 tolerance.
    let lanes = [EvaluationTestLane::oracle(
        EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions::default())]),
        EvaluationTestInput::TensorReduced,
    )];
    let mut passed = 0;
    for case in acceptance_inputs::cases() {
        let parent = match retained_parent_descriptor(&case.input) {
            acceptance_inputs::ParentDescriptor::H => "H",
            acceptance_inputs::ParentDescriptor::Fg => "FG",
            acceptance_inputs::ParentDescriptor::Bmw => "BMW",
            acceptance_inputs::ParentDescriptor::X => "X",
        };
        if filter.as_deref().is_some_and(|f| parent != f) {
            continue;
        }
        let mut settings = case.settings.clone();
        settings.form_exe_path = form.clone();
        let masses = vakint::params_from_f64(
            &case
                .parameters
                .iter()
                .map(|(name, value)| ((*name).into(), *value))
                .collect(),
            settings.run_time_decimal_precision,
        );
        let external = vakint::externals_from_f64(
            &case.external_momenta.iter().copied().collect(),
            settings.run_time_decimal_precision,
        );
        println!(
            "public inventory START {} (retained parent {parent})",
            case.name
        );
        test_utils::compare_vakint_evaluations_vs_reference(
            settings.clone(),
            &lanes,
            TensorPrepass::FeynKit,
            case.input.as_view(),
            masses,
            external,
            case.expected.clone(),
            settings.run_time_decimal_precision,
            1.0,
        );
        external_master_oracle::compare(&case, &form);
        println!("public inventory PASS {}", case.name);
        passed += 1;
    }
    assert!(passed > 0);
    if filter.is_none() {
        assert_eq!(passed, 15);
    }
}
