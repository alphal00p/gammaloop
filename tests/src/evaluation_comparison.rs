//! Physical event and threshold comparisons shared by integration-test binaries.

use gammalooprs::utils::F;
use itertools::Itertools;
use spenso::algebra::complex::Complex;

pub fn complex_distance(lhs: Complex<f64>, rhs: Complex<f64>) -> f64 {
    let delta_re = lhs.re - rhs.re;
    let delta_im = lhs.im - rhs.im;
    delta_re.hypot(delta_im)
}

pub fn assert_complex_approx_eq(
    actual: Complex<f64>,
    expected: Complex<f64>,
    context: impl std::fmt::Display,
) {
    let scale = actual
        .re
        .abs()
        .max(actual.im.abs())
        .max(expected.re.abs())
        .max(expected.im.abs());
    let relative_distance = if scale == 0.0 {
        0.0
    } else {
        complex_distance(actual / scale, expected / scale)
    };
    assert!(
        actual.re.is_finite()
            && actual.im.is_finite()
            && expected.re.is_finite()
            && expected.im.is_finite()
            && relative_distance <= 1.0e-10,
        "{context}: actual={actual}, expected={expected}, relative distance={relative_distance}, tolerance=1e-10"
    );
}

pub fn complex_ff64(value: &Complex<F<f64>>) -> Complex<f64> {
    Complex::new(value.re.0, value.im.0)
}

pub fn assert_f64_approx_eq(actual: f64, expected: f64, context: &str) {
    assert_complex_approx_eq(
        Complex::new(actual, 0.0),
        Complex::new(expected, 0.0),
        context,
    );
}

/// Exact source certificates for weights that vanish, expressed in physical
/// event units. Raw additional weights must first include their recorded event
/// normalization; decomposition values already include that normalization.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct CertifiedZeroWeights {
    pub original_cuts: std::collections::BTreeSet<(usize, usize)>,
    pub threshold_cuts: std::collections::BTreeSet<(usize, usize)>,
    pub threshold_components: std::collections::BTreeSet<(usize, usize, usize)>,
    pub all_cuts_zero: bool,
    pub dimensionless_tolerance: f64,
    pub e_cm: f64,
    pub energy_dimension: i32,
}

#[derive(Clone, Copy, PartialEq)]
pub enum EventWeightComparison<'a> {
    IndividualOrders,
    PhysicalResidues,
    CertifiedZeros {
        individual_orders: bool,
        weights: &'a CertifiedZeroWeights,
    },
}

pub fn assert_evaluation_outputs_match(
    actual: &gammalooprs::integrands::evaluation::EvaluationResultOutput,
    expected: &gammalooprs::integrands::evaluation::EvaluationResultOutput,
    context: &str,
    comparison: EventWeightComparison<'_>,
) {
    let (individual_orders, zeros) = match comparison {
        EventWeightComparison::IndividualOrders => (true, None),
        EventWeightComparison::PhysicalResidues => (false, None),
        EventWeightComparison::CertifiedZeros {
            individual_orders,
            weights,
        } => (individual_orders, Some(weights)),
    };
    let assert_weight = |mut actual: Complex<f64>,
                         mut expected: Complex<f64>,
                         context: &dyn std::fmt::Display,
                         certified_zero: bool,
                         remaining_factors: [Option<f64>; 2]| {
        if certified_zero {
            let zeros = zeros.expect("a zero comparison requires its exact source certificate");
            let physical_scale = zeros.e_cm.powi(zeros.energy_dimension);
            let tolerance = zeros.dimensionless_tolerance * physical_scale;
            assert!(
                zeros.dimensionless_tolerance.is_finite() && zeros.dimensionless_tolerance > 0.0
            );
            assert!(physical_scale.is_finite() && physical_scale > 0.0);
            assert!(tolerance.is_finite() && tolerance > 0.0);
            for (route, value) in [("actual", actual), ("expected", expected)] {
                assert!(
                    value.re.is_finite()
                        && value.im.is_finite()
                        && value.re.hypot(value.im) <= tolerance,
                    "{context}: {route} certified zero has value {value:e}, dimensionless tolerance={}, physical scale={}, bound={tolerance:e}",
                    zeros.dimensionless_tolerance,
                    physical_scale
                );
            }
        } else {
            // Match the ordinary reporting boundary: only completed physical
            // contributions below binary64's normal range may round to zero.
            // Bound each phase independently; a complex remaining factor can
            // promote either raw component into either physical phase.
            for (actual, expected) in [
                (&mut actual.re, &mut expected.re),
                (&mut actual.im, &mut expected.im),
            ] {
                let mut underflow = true;
                for (value, factor) in [*actual, *expected].into_iter().zip(remaining_factors) {
                    underflow &= factor.is_some_and(|factor| {
                        // Do not relax normal raw/bare values even if an outer
                        // factor suppresses them; only align reported subnormals.
                        let bound = value.abs() * factor.abs().max(1.0);
                        assert!(
                            value.is_finite() && factor.is_finite() && bound.is_finite(),
                            "{context}: nonfinite weight or completed contribution"
                        );
                        bound < f64::MIN_POSITIVE
                    });
                }
                if underflow {
                    *actual = 0.0;
                    *expected = 0.0;
                }
            }
            assert_complex_approx_eq(actual, expected, context);
        }
    };
    match (
        actual.parameterization_jacobian.as_ref(),
        expected.parameterization_jacobian.as_ref(),
    ) {
        (Some(actual), Some(expected)) => {
            assert_f64_approx_eq(actual.0, expected.0, &format!("{context}: jacobian"));
        }
        (None, None) => {}
        _ => panic!("{context}: jacobian presence differs"),
    }
    assert_f64_approx_eq(
        actual.integrator_weight.0,
        expected.integrator_weight.0,
        &format!("{context}: integrator weight"),
    );

    assert_eq!(
        actual.event_groups.len(),
        expected.event_groups.len(),
        "{context}: event-group count differs"
    );
    for (group_index, (actual_group, expected_group)) in actual
        .event_groups
        .iter()
        .zip(expected.event_groups.iter())
        .enumerate()
    {
        assert_eq!(
            actual_group.len(),
            expected_group.len(),
            "{context}: event count differs in group {group_index}"
        );
        let actual_events = actual_group
            .iter()
            .sorted_by_key(|event| {
                (
                    event.cut_info.graph_group_id,
                    event.cut_info.graph_id,
                    event.cut_info.cut_id,
                    event.cut_info.orientation_id,
                    event.cut_info.sampling_channel_id,
                )
            })
            .collect_vec();
        let expected_events = expected_group
            .iter()
            .sorted_by_key(|event| {
                (
                    event.cut_info.graph_group_id,
                    event.cut_info.graph_id,
                    event.cut_info.cut_id,
                    event.cut_info.orientation_id,
                    event.cut_info.sampling_channel_id,
                )
            })
            .collect_vec();

        for (event_index, (actual_event, expected_event)) in
            actual_events.iter().zip(expected_events.iter()).enumerate()
        {
            let event_context = format!(
                "{context}: group {group_index} event {event_index}; actual cut={:?}; expected cut={:?}; actual additional={:?}; expected additional={:?}",
                actual_event.cut_info,
                expected_event.cut_info,
                actual_event.additional_weights.weights,
                expected_event.additional_weights.weights
            );
            assert_eq!(
                actual_event.cut_info.cut_id, expected_event.cut_info.cut_id,
                "{event_context}: cut id differs"
            );
            assert_eq!(
                actual_event.cut_info.graph_id, expected_event.cut_info.graph_id,
                "{event_context}: graph id differs"
            );
            assert_eq!(
                actual_event.cut_info.graph_group_id, expected_event.cut_info.graph_group_id,
                "{event_context}: graph-group id differs"
            );
            assert_eq!(
                actual_event.cut_info.orientation_id, expected_event.cut_info.orientation_id,
                "{event_context}: orientation id differs"
            );
            assert_eq!(
                actual_event.cut_info.sampling_channel_id,
                expected_event.cut_info.sampling_channel_id,
                "{event_context}: sampling channel differs"
            );
            assert_eq!(
                actual_event.cut_info.sampling_channel_edge_ids,
                expected_event.cut_info.sampling_channel_edge_ids,
                "{event_context}: sampling channel edges differ"
            );
            assert_eq!(
                actual_event.cut_info.particle_pdgs, expected_event.cut_info.particle_pdgs,
                "{event_context}: observable particle identities differ"
            );
            for (actual_momenta, expected_momenta) in [
                (
                    actual_event.kinematic_configuration.0.as_slice(),
                    expected_event.kinematic_configuration.0.as_slice(),
                ),
                (
                    actual_event.kinematic_configuration.1.as_slice(),
                    expected_event.kinematic_configuration.1.as_slice(),
                ),
            ] {
                assert_eq!(actual_momenta.len(), expected_momenta.len());
                for (actual, expected) in actual_momenta.iter().zip(expected_momenta) {
                    for (actual, expected) in actual.into_iter().zip(expected) {
                        assert_f64_approx_eq(
                            actual.0,
                            expected.0,
                            &format!("{event_context}: observable momentum"),
                        );
                    }
                }
            }
            let cut_key = (actual_event.cut_info.graph_id, actual_event.cut_info.cut_id);
            let zero_original = zeros.is_some_and(|zeros| zeros.original_cuts.contains(&cut_key));
            let zero_thresholds =
                zeros.is_some_and(|zeros| zeros.threshold_cuts.contains(&cut_key));
            assert_weight(
                complex_ff64(&actual_event.weight),
                complex_ff64(&expected_event.weight),
                &format_args!("{event_context}: event weight"),
                zero_original && zero_thresholds,
                [Some(1.0); 2],
            );
            assert_eq!(
                actual_event.additional_weights.weights.keys().collect_vec(),
                expected_event
                    .additional_weights
                    .weights
                    .keys()
                    .collect_vec(),
                "{event_context}: additional-weight keys differ"
            );
            for (key, expected_weight) in &expected_event.additional_weights.weights {
                let actual_weight = &actual_event.additional_weights.weights[key];
                use gammalooprs::observables::events::AdditionalWeightKey;
                if *key == AdditionalWeightKey::FullMultiplicativeFactor {
                    assert_complex_approx_eq(
                        complex_ff64(actual_weight),
                        complex_ff64(expected_weight),
                        format_args!("{event_context}: additional weight {key:?}"),
                    );
                    continue;
                }
                let certified_zero = match key {
                    AdditionalWeightKey::Original => zero_original,
                    AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 } => {
                        zero_thresholds
                    }
                    _ => false,
                };
                // These additional weights precede event normalization. Transport
                // them to physical units before using an E_cm-based zero bound.
                let physical_weight =
                    |event: &gammalooprs::observables::events::GenericEvent<f64>, weight| {
                        if certified_zero {
                            complex_ff64(weight)
                                * complex_ff64(
                                    &event.additional_weights.weights
                                        [&AdditionalWeightKey::FullMultiplicativeFactor],
                                )
                        } else {
                            complex_ff64(weight)
                        }
                    };
                assert_weight(
                    physical_weight(actual_event, actual_weight),
                    physical_weight(expected_event, expected_weight),
                    &format_args!("{event_context}: additional weight {key:?}"),
                    certified_zero,
                    [actual_event, expected_event].map(|event| {
                        event
                            .additional_weights
                            .weights
                            .get(&AdditionalWeightKey::FullMultiplicativeFactor)
                            .and_then(|factor| {
                                // Both raw phases can reinforce one physical
                                // phase. An overflowing bound cannot certify
                                // underflow, even when the factor itself is finite.
                                let bound = factor.re.0.abs() + factor.im.0.abs();
                                let mixing = factor.re.0 != 0.0 && factor.im.0 != 0.0;
                                let weight = &event.additional_weights.weights[key];
                                let joint_bound = weight.re.0.abs().max(weight.im.0.abs()) * bound;
                                // For genuine phase mixing, decide from both
                                // ORIGINAL raw components before either is erased.
                                // Real or imaginary factors only scale/permute them.
                                (bound.is_finite() && (!mixing || joint_bound < f64::MIN_POSITIVE))
                                    .then_some(bound)
                            })
                    }),
                );
            }
            match (
                &actual_event.additional_weights.threshold_counterterms,
                &expected_event.additional_weights.threshold_counterterms,
            ) {
                (None, None) => {}
                (Some(actual), Some(expected)) => {
                    use gammalooprs::observables::events::{
                        GenericThresholdCountertermComponentWeight as Component,
                        ThresholdCountertermComponentOccurrence as Occurrence,
                    };
                    let assert_metadata =
                        |actual: &Component<f64>, expected: &Component<f64>, context: &str| {
                            assert_eq!(
                                actual.component_id, expected.component_id,
                                "{context}: component id differs"
                            );
                            assert_eq!(
                                actual.evaluation_skipped, expected.evaluation_skipped,
                                "{context}: skipped status differs"
                            );
                            assert_eq!(
                                actual.bare.is_some(),
                                expected.bare.is_some(),
                                "{context}: bare weight presence differs"
                            );
                            assert_eq!(
                                actual.multiplier_values.len(),
                                expected.multiplier_values.len(),
                                "{context}: multiplier count differs"
                            );
                            for (actual, expected) in actual
                                .multiplier_values
                                .iter()
                                .zip(&expected.multiplier_values)
                            {
                                assert_f64_approx_eq(
                                    actual.0,
                                    expected.0,
                                    &format!("{context}: multiplier"),
                                );
                            }
                            assert_f64_approx_eq(
                                actual.effective_multiplier.0,
                                expected.effective_multiplier.0,
                                &format!("{context}: effective multiplier"),
                            );
                        };
                    assert_weight(
                        complex_ff64(&actual.original),
                        complex_ff64(&expected.original),
                        &format_args!("{event_context}: threshold decomposition original"),
                        zero_original,
                        [Some(1.0); 2],
                    );
                    let is_raised = |component: &Component<f64>| match &component.occurrence {
                        Occurrence::LocalUnitarity {
                            left_threshold_order,
                            right_threshold_order,
                            lu_cut_order,
                            ..
                        } => [left_threshold_order, right_threshold_order, lu_cut_order]
                            .into_iter()
                            .any(|order| order.is_some_and(|order| order > 1)),
                        Occurrence::Amplitude { .. } => false,
                    };
                    // A vanishing physical residue may contain nonzero, compensating
                    // derivative-order coefficients. Keep within-CFF comparisons of
                    // those partial weights strict; only their physical sum is zero.
                    let compare_zero_components = !individual_orders
                        || !actual
                            .components
                            .iter()
                            .chain(&expected.components)
                            .any(is_raised);
                    let (actual_components, expected_components) = if !individual_orders {
                        // A graph-registry component fixes the threshold variants and
                        // L/I kind. Overlap groups and the occupied derivative axes fix
                        // its physical occurrence. Only the derivative orders may move
                        // between partial weights of that same residue.
                        let groups = |components: &[Component<f64>]| {
                            let mut groups = std::collections::BTreeMap::new();
                            for component in components {
                                let key = match &component.occurrence {
                                    Occurrence::Amplitude {
                                        raised_esurface_id,
                                        overlap_group,
                                    } => (
                                        component.component_id,
                                        Some((*raised_esurface_id, *overlap_group)),
                                        Vec::new(),
                                        [false; 3],
                                    ),
                                    Occurrence::LocalUnitarity {
                                        overlap_groups,
                                        left_threshold_order,
                                        right_threshold_order,
                                        lu_cut_order,
                                    } => (
                                        component.component_id,
                                        None,
                                        overlap_groups.to_vec(),
                                        [
                                            left_threshold_order.is_some(),
                                            right_threshold_order.is_some(),
                                            lu_cut_order.is_some(),
                                        ],
                                    ),
                                };
                                groups
                                    .entry(key)
                                    .or_insert_with(Vec::new)
                                    .push(component.clone());
                            }
                            groups
                        };
                        let actual_groups = groups(&actual.components);
                        let expected_groups = groups(&expected.components);
                        assert_eq!(
                            actual_groups.keys().collect_vec(),
                            expected_groups.keys().collect_vec(),
                            "{event_context}: physical threshold residues differ"
                        );
                        let mut actual_components = Vec::new();
                        let mut expected_components = Vec::new();
                        for (key, expected_group) in &expected_groups {
                            let actual_group = &actual_groups[key];
                            let raised = actual_group.iter().chain(expected_group).any(is_raised);
                            if !raised {
                                actual_components.extend(actual_group.iter().cloned());
                                expected_components.extend(expected_group.iter().cloned());
                                continue;
                            }
                            let reference = &expected_group[0];
                            let component_context =
                                format!("{event_context}: physical residue {key:?}");
                            for component in actual_group.iter().chain(expected_group) {
                                assert_metadata(component, reference, &component_context);
                            }
                            let sum = |components: &[Component<f64>]| {
                                let mut result = reference.clone();
                                if let Occurrence::LocalUnitarity {
                                    left_threshold_order,
                                    right_threshold_order,
                                    lu_cut_order,
                                    ..
                                } = &mut result.occurrence
                                {
                                    for order in
                                        [left_threshold_order, right_threshold_order, lu_cut_order]
                                    {
                                        *order = order.map(|_| 1);
                                    }
                                }
                                result.weighted = components
                                    .iter()
                                    .fold(Complex::new(F(0.0), F(0.0)), |sum, component| {
                                        sum + component.weighted
                                    });
                                result.bare = reference.bare.map(|_| {
                                    components
                                        .iter()
                                        .fold(Complex::new(F(0.0), F(0.0)), |sum, component| {
                                            sum + component.bare.unwrap()
                                        })
                                });
                                result
                            };
                            actual_components.push(sum(actual_group));
                            expected_components.push(sum(expected_group));
                        }
                        (actual_components, expected_components)
                    } else {
                        (actual.components.clone(), expected.components.clone())
                    };
                    assert_eq!(
                        actual_components.len(),
                        expected_components.len(),
                        "{event_context}: threshold component count differs"
                    );
                    for (actual, expected) in actual_components.iter().zip(&expected_components) {
                        let component_context = format!(
                            "{event_context}: threshold component {} {:?}",
                            expected.component_id, expected.occurrence
                        );
                        assert_metadata(actual, expected, &component_context);
                        assert_eq!(
                            actual.occurrence, expected.occurrence,
                            "{component_context}: occurrence differs"
                        );
                        let actual_component_id = actual.component_id;
                        let zero_component = compare_zero_components
                            && zeros.is_some_and(|zeros| {
                                zeros.threshold_components.contains(&(
                                    cut_key.0,
                                    cut_key.1,
                                    actual_component_id,
                                ))
                            });
                        let bare_factors = [actual, expected]
                            .map(|component| Some(component.effective_multiplier.0));
                        match (&actual.bare, &expected.bare) {
                            (None, None) => {}
                            (Some(actual), Some(expected)) => assert_weight(
                                complex_ff64(actual),
                                complex_ff64(expected),
                                &format_args!("{component_context}: bare weight"),
                                zero_component,
                                // Bare values already include event normalization.
                                // Preserve their own comparison and any promotion
                                // by the remaining dimensionless user multiplier.
                                bare_factors,
                            ),
                            _ => panic!("{component_context}: bare weight presence differs"),
                        }
                        assert_weight(
                            complex_ff64(&actual.weighted),
                            complex_ff64(&expected.weighted),
                            &format_args!("{component_context}: weighted contribution"),
                            zero_component,
                            [Some(1.0); 2],
                        );
                    }
                }
                _ => panic!("{event_context}: threshold decomposition presence differs"),
            }
        }
    }
    assert_weight(
        complex_ff64(&actual.integrand_result),
        complex_ff64(&expected.integrand_result),
        &format_args!("{context}: integrand result; actual={actual:?}; expected={expected:?}"),
        zeros.is_some_and(|zeros| zeros.all_cuts_zero),
        [actual, expected].map(|output| {
            Some(
                output
                    .parameterization_jacobian
                    .as_ref()
                    .map_or(1.0, |jacobian| jacobian.0)
                    * output.integrator_weight.0,
            )
        }),
    );
}
