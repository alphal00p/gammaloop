use super::{
    InspectFitStatus, InspectResult, UV_PROFILE_MAX_DOD, inspect_failure_reason,
    inspect_results_need_arbprec_retry, log_log_slope,
};
use crate::integrands::evaluation::{EvaluationMetaData, GenericEvaluationResult};
use crate::observables::events::{GenericEvent, GenericEventGroupList};
use crate::processes::CutId;
use crate::utils::{ArbPrec, F, QuadFloat};
use spenso::algebra::complex::Complex;
use symbolica::domains::float::DoubleFloat;

#[test]
fn nonfinite_profile_components_request_retry_without_hiding_invalid_points() {
    let one = F::<ArbPrec>::default().one();
    let scales = [1.0, 2.0, 4.0, 8.0, 16.0];
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        for imaginary in [false, true] {
            for select_cut in [false, true] {
                let value = if imaginary {
                    Complex::new(F(1.0), F(invalid))
                } else {
                    Complex::new(F(invalid), F(1.0))
                };
                let mut event = GenericEvent::default();
                event.cut_info.cut_id = 0;
                event.weight = value;
                let mut events = GenericEventGroupList::default();
                events.push_singleton(event);
                let result = GenericEvaluationResult {
                    reference_moments: None,
                    integrand_result: value,
                    absolute_integrand_result: None,
                    parameterization_jacobian: None,
                    integrator_weight: F(1.0),
                    event_groups: events,
                    evaluation_metadata: EvaluationMetaData::new_empty(),
                };
                let selected = [CutId(0)];
                let invalid_sample = InspectResult::from_evaluation(
                    result,
                    &one,
                    select_cut.then_some(selected.as_slice()),
                );
                assert!(invalid_sample.magnitude.is_none());

                // Four perfect power-law points cannot hide an invalid fifth point.
                let mut samples: Vec<_> = scales
                    .iter()
                    .map(|scale| InspectResult {
                        magnitude: Some(F::<ArbPrec>::from_f64(1.0 / scale)),
                    })
                    .collect();
                assert!(!inspect_results_need_arbprec_retry(&samples, &scales));
                let valid_sample = samples[2].clone();
                samples[2] = invalid_sample;
                assert!(inspect_results_need_arbprec_retry(&samples, &scales));
                assert!(log_log_slope(&samples, &scales).is_none());
                let status = InspectFitStatus::from_results(&samples, &scales);
                assert!(status.decay_slope.is_none());
                assert_eq!(
                    inspect_failure_reason(None, status, UV_PROFILE_MAX_DOD, true),
                    Some("missing_fit")
                );

                // Replacing the failed evaluation with a valid retry restores the fit.
                samples[2] = valid_sample;
                let fit = log_log_slope(&samples, &scales).unwrap();
                assert!((fit.slope + 1.0).abs() < 1.0e-12);
                assert!(!inspect_results_need_arbprec_retry(&samples, &scales));
            }
        }
    }
}

#[test]
fn invalid_weights_and_sanitized_zeros_cannot_certify_a_zero_tail() {
    let one = F::<ArbPrec>::default().one();
    for weight in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1.0] {
        let mut metadata = EvaluationMetaData::new_empty();
        metadata.is_nan = weight == 1.0;
        let result = GenericEvaluationResult {
            reference_moments: None,
            integrand_result: Complex::new_re(F(0.0)),
            absolute_integrand_result: None,
            parameterization_jacobian: None,
            integrator_weight: F(weight),
            event_groups: GenericEventGroupList::default(),
            evaluation_metadata: metadata,
        };
        let sample = InspectResult::from_evaluation(result, &one, None);
        assert!(sample.magnitude.is_none());
        let mut samples = vec![
            InspectResult {
                magnitude: Some(one.zero()),
            };
            5
        ];
        let scales = [1.0, 2.0, 4.0, 8.0, 16.0];
        assert_eq!(
            InspectFitStatus::from_results(&samples, &scales).decay_slope,
            Some(f64::NEG_INFINITY)
        );
        samples[2] = sample;
        assert!(inspect_results_need_arbprec_retry(&samples, &scales));
        assert!(log_log_slope(&samples, &scales).is_none());
        assert!(
            InspectFitStatus::from_results(&samples, &scales)
                .decay_slope
                .is_none()
        );
    }
}

#[test]
fn cut_projection_ignores_unselected_events_and_applies_each_weight_once() {
    let one = F::<ArbPrec>::default().one();
    let prefactor = one.from_usize(11);
    let mut events = GenericEventGroupList::default();
    for (cut_id, value) in [
        (0, Complex::new_re(F(2.0))),
        (1, Complex::new(F(f64::NAN), F(f64::INFINITY))),
    ] {
        let mut event = GenericEvent::default();
        event.cut_info.cut_id = cut_id;
        event.weight = value;
        events.push_singleton(event);
    }
    let result = GenericEvaluationResult {
        reference_moments: None,
        integrand_result: Complex::new_re(F(2.0)),
        absolute_integrand_result: None,
        parameterization_jacobian: None,
        integrator_weight: F(7.0),
        event_groups: events,
        evaluation_metadata: EvaluationMetaData::new_empty(),
    };
    let selected = InspectResult::from_evaluation(result.clone(), &prefactor, Some(&[CutId(0)]));
    assert_eq!(selected.magnitude, Some(one.from_usize(22)));
    let absent = InspectResult::from_evaluation(result.clone(), &prefactor, Some(&[CutId(2)]));
    assert_eq!(absent.magnitude, Some(one.zero()));
    let total = InspectResult::from_evaluation(result, &prefactor, None);
    assert_eq!(total.magnitude, Some(one.from_usize(154)));
}

#[test]
fn cut_projection_preserves_cancellation_before_weighting() {
    let one = F::<ArbPrec>::default().one();
    let prefactor = one.from_usize(10).powi(300);
    let residual = 2.0_f64.powi(-60);
    let mut events = GenericEventGroupList::default();
    for value in [1.0, residual, -1.0] {
        let mut event = GenericEvent::default();
        event.cut_info.cut_id = 0;
        event.weight = Complex::new_im(F(value));
        events.push_singleton(event);
    }
    let result = GenericEvaluationResult {
        reference_moments: None,
        integrand_result: Complex::new_re(F(0.0)),
        absolute_integrand_result: None,
        parameterization_jacobian: None,
        integrator_weight: F(7.0),
        event_groups: events,
        evaluation_metadata: EvaluationMetaData::new_empty(),
    };
    let selected = InspectResult::from_evaluation(result, &prefactor, Some(&[CutId(0)]));
    assert_eq!(
        selected.magnitude,
        Some(F(residual).to_arb_exact().unwrap() * prefactor)
    );
}

#[test]
fn separated_quad_limbs_preserve_residues_or_request_an_arb_retry() {
    let one = F::<ArbPrec>::default().one();
    let prefactor = one.from_usize(10).powi(300);
    for residual in [2.0_f64.powi(-400), f64::from_bits(1)] {
        let mut events = GenericEventGroupList::default();
        for value in [
            DoubleFloat::from_compensated_sum(1.0, residual),
            DoubleFloat::from(-1.0),
        ] {
            let mut event = GenericEvent::default();
            event.cut_info.cut_id = 0;
            event.weight = Complex::new_re(F(QuadFloat::from(value)));
            events.push_singleton(event);
        }
        let result = GenericEvaluationResult {
            reference_moments: None,
            integrand_result: Complex::new_re(F(QuadFloat::from(DoubleFloat::from(0.0)))),
            absolute_integrand_result: None,
            parameterization_jacobian: None,
            integrator_weight: F(QuadFloat::from(DoubleFloat::from(1.0))),
            event_groups: events,
            evaluation_metadata: EvaluationMetaData::new_empty(),
        };
        let selected = InspectResult::from_evaluation(result, &prefactor, Some(&[CutId(0)]));
        let expected = F(residual).to_arb_exact().unwrap() * &prefactor;
        if residual == f64::from_bits(1) {
            // Its exact source spans 1075 bits. Reject it locally instead of rounding
            // away the residue or aborting the profile before an Arb reevaluation.
            assert!(selected.magnitude.is_none());
            assert!(inspect_results_need_arbprec_retry(&[selected], &[1.0]));
            let mut event = GenericEvent::default();
            event.cut_info.cut_id = 0;
            event.weight = Complex::new_re(F(residual).to_arb_exact().unwrap());
            let mut events = GenericEventGroupList::default();
            events.push_singleton(event);
            let arb_result = GenericEvaluationResult {
                reference_moments: None,
                integrand_result: Complex::new_re(one.zero()),
                absolute_integrand_result: None,
                parameterization_jacobian: None,
                integrator_weight: one.clone(),
                event_groups: events,
                evaluation_metadata: EvaluationMetaData::new_empty(),
            };
            let retried = InspectResult::from_evaluation(arb_result, &prefactor, Some(&[CutId(0)]));
            assert_eq!(retried.magnitude, Some(expected));
        } else {
            assert_eq!(selected.magnitude, Some(expected));
        }
    }
}

#[test]
fn profile_magnitudes_preserve_extreme_double_components() {
    let one = F::<ArbPrec>::default().one();
    for scale in [1.0e200, 1.0e-200] {
        let result = GenericEvaluationResult {
            reference_moments: None,
            integrand_result: Complex::new(F(3.0 * scale), F(4.0 * scale)),
            absolute_integrand_result: None,
            parameterization_jacobian: None,
            integrator_weight: F(1.0),
            event_groups: GenericEventGroupList::default(),
            evaluation_metadata: EvaluationMetaData::new_empty(),
        };
        let sample = InspectResult::from_evaluation(result, &one, None);
        let magnitude = sample.magnitude.unwrap().into_ff64().0;
        assert!((magnitude / (5.0 * scale) - 1.0).abs() < 1.0e-15);
    }
}
