//! Keep the historical table-center reference distinct from its scientific
//! uncertainty. This test-only oracle declares the original stored centers as
//! exact rational points; the public backend still enforces every source cap.

use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    coefficient::CoefficientView,
    domains::float::{Complex, Float},
};
use vakint::{
    EvaluationOrder, FMFTOptions, NumericalEvaluationResult, TensorReductionMethod, Vakint,
    VakintSettings,
};

use std::collections::HashMap;
use vakint::Momentum;

fn frozen_centers(expression: Atom) -> Atom {
    let substitutions: Vec<_> = vakint::fmft_numerics::MASTERS_NUMERIC_SUBSTITUTIONS
        .iter()
        .map(|(source, (value, condition))| {
            let exact = value.replace_map(|part, _, out| {
                if let AtomView::Num(number) = part
                    && let CoefficientView::Float(real, imaginary) = number.get_coeff_view()
                {
                    **out = Atom::num(real.to_float().to_rational())
                        + Atom::i() * Atom::num(imaginary.to_float().to_rational());
                }
            });
            (source.clone(), exact, condition.clone())
        })
        .collect();
    let mut result = expression;
    result.repeat_map(Box::new(move |part: AtomView| {
        let mut value = part.to_owned();
        for (source, target, condition) in &substitutions {
            value = value
                .replace(source.to_pattern())
                .when(condition)
                .with(target.to_pattern());
        }
        value
    }));
    result
}

pub(super) fn compare_scalar(
    settings: &VakintSettings,
    scalar: &Atom,
    parameters: &HashMap<String, Float, ahash::RandomState>,
    external: &HashMap<usize, Momentum, ahash::RandomState>,
    expected_output: &[(i64, (String, String))],
    max_pull: f64,
) {
    let engine = Vakint::new().unwrap();
    let settings = VakintSettings {
        tensor_reduction_method: TensorReductionMethod::FeynKit,
        use_dot_product_notation: true,
        evaluation_order: EvaluationOrder::fmft_only(Some(FMFTOptions {
            susbstitute_masters: false,
            ..Default::default()
        })),
        ..settings.clone()
    };
    let symbolic = engine
        .evaluate_integral(&settings, scalar.as_view())
        .unwrap();
    let frozen = frozen_centers(symbolic);
    let working = VakintSettings {
        run_time_decimal_precision: 80,
        ..settings.clone()
    };
    let (central, _) = engine
        .numerical_evaluation(
            &working,
            frozen.as_view(),
            parameters,
            &Default::default(),
            Some(external),
        )
        .unwrap();
    let bits = working.get_binary_precision();
    let reference = NumericalEvaluationResult::from_values(
        expected_output
            .iter()
            .map(|(power, (real, imaginary))| {
                (
                    *power,
                    Complex::new(
                        Float::parse(real, Some(bits)).unwrap(),
                        Float::parse(imaginary, Some(bits)).unwrap(),
                    ),
                )
            })
            .collect(),
        80,
    );
    let (matches, detail) = central.does_approx_match(
        &reference,
        None,
        0.1_f64.powi(settings.run_time_decimal_precision as i32 - 2),
        max_pull,
    );
    assert!(
        matches,
        "historical FORM table-center reference: {detail}\n{central}"
    );

    let public_settings = VakintSettings {
        evaluation_order: EvaluationOrder::fmft_only(None),
        ..settings
    };
    let public = engine
        .evaluate_integral(&public_settings, scalar.as_view())
        .unwrap();
    let (tracked, _) = engine
        .numerical_evaluation(
            &public_settings,
            public.as_view(),
            parameters,
            &Default::default(),
            Some(external),
        )
        .unwrap();
    assert_eq!(
        tracked.coefficients().keys().collect::<Vec<_>>(),
        central.coefficients().keys().collect::<Vec<_>>()
    );
    for (power, actual) in tracked.coefficients() {
        let expected = &central.coefficients()[power];
        for (label, actual_value, actual_precision, expected_value, expected_precision) in [
            (
                "real",
                &actual.value().re,
                actual.precision().real(),
                &expected.value().re,
                expected.precision().real(),
            ),
            (
                "imaginary",
                &actual.value().im,
                actual.precision().imaginary(),
                &expected.value().im,
                expected.precision().imaginary(),
            ),
        ] {
            let difference = (actual_value.to_rational() - expected_value.to_rational()).abs();
            let bound = actual_precision.absolute_error().to_rational()
                + expected_precision.absolute_error().to_rational();
            assert!(
                difference <= bound,
                "epsilon^{power} {label}: the public tracked error does not contain the independently evaluated legacy table center"
            );
        }
    }
}
