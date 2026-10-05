//! End-to-end tests through the public backend, without a test-local reducer.

use super::{acceptance_inputs, retained_parent_descriptor, test_utils};
use test_utils::{EvaluationTestInput, EvaluationTestLane, RustRedParityPolicy, TensorPrepass};
use vakint::{EvaluationMethod, EvaluationOrder, FMFTOptions};

/// Frozen 32-digit numerator reference, with neither an executable oracle nor a
/// FORM tensor prepass. This uses the same input, scales, external vectors and
/// tolerance as the existing FMFT/RustRed comparison above.
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
    let masses = vakint::params_from_f64(
        &case
            .parameters
            .iter()
            .map(|(name, value)| ((*name).into(), *value))
            .collect(),
        settings.run_time_decimal_precision,
    );
    let external = vakint::externals_from_f64(
        &case.external_momenta.into_iter().collect(),
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
        case.expected,
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
    let lanes = [
        EvaluationTestLane::oracle(
            EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions::default())]),
            EvaluationTestInput::TensorReduced,
        ),
        EvaluationTestLane::rustred_scalar(
            EvaluationTestInput::TensorReduced,
            RustRedParityPolicy::NumericalOnly,
        ),
    ];
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
        let mut settings = case.settings;
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
            &case.external_momenta.into_iter().collect(),
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
            case.expected,
            settings.run_time_decimal_precision,
            1.0,
        );
        println!("public inventory PASS {}", case.name);
        passed += 1;
    }
    assert!(passed > 0);
    if filter.is_none() {
        assert_eq!(passed, 15);
    }
}
