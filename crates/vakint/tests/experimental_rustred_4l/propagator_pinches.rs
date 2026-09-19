//! Independent numerator-cancellation and graph-contraction comparisons.
//! These finite tests do not certify arbitrary-index closure.

use super::acceptance_inputs::ParentDescriptor;
use super::input::{NumeratorPropagatorPinch, ParentInput};
use super::retained_parent_descriptor;
use super::test_utils::{
    self, EvaluationTestInput, EvaluationTestLane, RustRedParityPolicy, TensorPrepass,
};
use symbolica::atom::{AtomCore, AtomView};
use vakint::{
    EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor, Vakint, VakintSettings,
};

const PARENTS: [ParentDescriptor; 4] = [
    ParentDescriptor::H,
    ParentDescriptor::Fg,
    ParentDescriptor::Bmw,
    ParentDescriptor::X,
];

fn label(parent: ParentDescriptor) -> &'static str {
    match parent {
        ParentDescriptor::H => "H",
        ParentDescriptor::Fg => "FG",
        ParentDescriptor::Bmw => "BMW",
        ParentDescriptor::X => "X",
    }
}

fn cases(parent: ParentDescriptor) -> Vec<NumeratorPropagatorPinch> {
    let input = ParentInput::from_csv(parent.csv());
    // Each double contraction is a forest in all four input graphs. In
    // particular, do not cancel both FG edges 5 and 7: they carry its only
    // fourth-loop denominators and would give a vacuous scaleless zero.
    [
        ("D1", "P1", &[0][..]),
        ("D7", "P7", &[6][..]),
        ("D1D2", "P1P2", &[0, 1][..]),
        ("D5D6", "P5P6", &[4, 5][..]),
    ]
    .into_iter()
    .map(|(numerator, pinch, slots)| input.numerator_equals_propagators(slots, numerator, pinch))
    .collect()
}

#[test]
fn sixteen_inputs_are_expanded_and_have_registered_four_loop_witnesses() {
    Vakint::initialize_vakint_symbols();
    let mut count = 0;
    for parent in PARENTS {
        for pair in cases(parent) {
            // Both inputs go independently through the existing matcher.
            let first = retained_parent_descriptor(&pair.with_numerator);
            let second = retained_parent_descriptor(&pair.pinched);
            assert_eq!(first, parent);
            pair.with_numerator.visitor(&mut |node| {
                if let AtomView::Fun(fun) = node {
                    if fun.get_symbol() == vakint::symbols::S.dot {
                        for argument in fun.iter() {
                            assert!(matches!(argument, AtomView::Fun(f) if f.get_symbol() == vakint::symbols::S.k), "unexpanded scalar product: {node}");
                        }
                    }
                }
                true
            });
            println!(
                "{}/{} -> {} parent",
                label(parent),
                pair.pinch_name,
                label(second)
            );
            count += 1;
        }
    }
    assert_eq!(count, 16);
}

#[test]
#[ignore = "independent FMFT regression; requires explicit oracle executable"]
fn fg_mixed_momentum_cancellation_fmft() {
    test_utils::run_multi_lane_acceptance(|| {
        let mut vakint = test_utils::get_vakint(VakintSettings {
            form_exe_path: std::env::var("VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH").unwrap(),
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions::default())]),
            run_time_decimal_precision: 32,
            number_of_terms_in_epsilon_expansion: 5,
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            use_dot_product_notation: true,
            ..VakintSettings::default()
        });
        let pair = cases(ParentDescriptor::Fg).remove(1);
        let masses = vakint::params_from_f64(
            &[("muvsq".into(), 3.0), ("mursq".into(), 7.0)]
                .into_iter()
                .collect(),
            32,
        );
        let mut values = Vec::new();
        for (name, input) in [
            ("numerator", pair.with_numerator),
            ("pinched", pair.pinched),
        ] {
            let canonical = vakint.to_canonical(input.as_view(), true).unwrap();
            vakint.settings.tensor_reduction_method = vakint::TensorReductionMethod::FeynKit;
            let tensor = vakint.tensor_reduce(canonical.as_view()).unwrap();
            println!("{name}: input={input}; canonical={canonical}; tensor={tensor}");
            let value = vakint.evaluate_integral(tensor.as_view()).unwrap();
            let numerical = Vakint::full_numerical_evaluation(
                &vakint.settings,
                value.as_view(),
                &masses,
                &Default::default(),
                Some(&Default::default()),
            )
            .unwrap()
            .0;
            println!("{name}: {numerical}");
            values.push(numerical);
        }
        let (matches, message) = values[0].does_approx_match(&values[1], None, 1.0e-20, 0.0);
        assert!(matches, "{message}");
    });
}

#[test]
#[ignore = "release four-loop public RustRed acceptance and independent FMFT oracle"]
fn sixteen_numerator_propagator_pinches_match_fmft_and_rustred() {
    test_utils::run_multi_lane_acceptance(|| {
        Vakint::initialize_vakint_symbols();
        let form = std::env::var("VAKINT_4L_CANDIDATE_ORACLE_FORM_PATH")
            .expect("explicit FMFT oracle executable required");
        let filter = std::env::var("VAKINT_4L_CANDIDATE_FAMILY_FILTER").ok();
        assert!(
            filter
                .as_deref()
                .is_none_or(|f| PARENTS.iter().any(|p| label(*p) == f)),
            "unknown parent filter"
        );
        let mut passed = 0;
        for parent in PARENTS {
            if filter.as_deref().is_some_and(|f| label(parent) != f) {
                continue;
            }
            for pair in cases(parent) {
                let first = retained_parent_descriptor(&pair.with_numerator);
                let second = retained_parent_descriptor(&pair.pinched);
                let oracle = EvaluationTestLane::oracle(
                    EvaluationOrder(vec![EvaluationMethod::FMFT(FMFTOptions::default())]),
                    EvaluationTestInput::TensorReduced,
                );
                // The public backend loads shipped programs once. Neither
                // member of the pair receives a test-local reducer or solve.
                let native = EvaluationTestLane::rustred_scalar(
                    EvaluationTestInput::TensorReduced,
                    RustRedParityPolicy::NumericalOnly,
                );
                let settings = VakintSettings {
                    form_exe_path: form.clone(),
                    run_time_decimal_precision: 32,
                    number_of_terms_in_epsilon_expansion: 5,
                    integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
                    ..VakintSettings::default()
                };
                // Non-unit scales make a dropped mass term/power observable.
                let masses = vakint::params_from_f64(
                    &[("muvsq".into(), 3.0), ("mursq".into(), 7.0)]
                        .into_iter()
                        .collect(),
                    32,
                );
                println!(
                    "pinch START {}/{} == {} (owners {}/{})",
                    label(parent),
                    pair.numerator_name,
                    pair.pinch_name,
                    label(first),
                    label(second)
                );
                test_utils::compare_input_pair(
                    settings,
                    &[(oracle.clone(), oracle), (native.clone(), native)],
                    TensorPrepass::FeynKit,
                    pair.with_numerator.as_view(),
                    pair.pinched.as_view(),
                    masses,
                    Default::default(),
                    1.0e-20,
                    0.0,
                );
                println!(
                    "public pinch PASS {}/{} == {}",
                    label(parent),
                    pair.numerator_name,
                    pair.pinch_name
                );
                passed += 1;
            }
        }
        assert_eq!(passed, if filter.is_some() { 4 } else { 16 });
    });
}
