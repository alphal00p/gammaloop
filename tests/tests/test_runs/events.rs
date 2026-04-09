use super::utils::*;
use super::*;

#[test]
#[serial]
fn amplitude_events_surface_threshold_counterterms_and_reproduce_weight() -> Result<()> {
    let test_name = "amplitude_events_surface_threshold_counterterms";
    let mut cli = setup_gg_hhh_threshold_amplitude_cli(test_name)?;
    let process_ref = ProcessRef::Unqualified("gg_hhh".to_string());
    let integrand_name = "1L".to_string();
    let (process_id, resolved_integrand_name) = cli
        .state
        .find_integrand_ref(Some(&process_ref), Some(&integrand_name))?;

    let points = arr2(&[[0.23, 0.41, 0.67]]);
    let discrete_dims = arr2(&[[0usize, 0usize, 0usize]]);
    let result = evaluate_sample(
        &mut cli.state,
        &EvaluateSamples {
            process_id: Some(process_id),
            integrand_name: Some(resolved_integrand_name),
            use_arb_prec: false,
            minimal_output: false,
            return_generated_events: None,
            momentum_space: false,
            points: points.view(),
            integrator_weights: None,
            discrete_dims: Some(discrete_dims.view()),
            graph_names: None,
            orientations: None,
        },
    )?;

    let evaluation = result.sample.evaluation;
    assert_eq!(evaluation.event_groups.len(), 1);
    assert_eq!(evaluation.event_groups[0].len(), 1);

    let event = &evaluation.event_groups[0][0];
    assert_eq!(event.cut_info.cut_id, 0);
    assert_eq!(event.cut_info.particle_pdgs.0.as_slice(), &[21, 21]);
    assert_eq!(event.cut_info.particle_pdgs.1.as_slice(), &[25, 25, 25]);
    assert_eq!(event.kinematic_configuration.0.len(), 2);
    assert_eq!(event.kinematic_configuration.1.len(), 3);

    let to_ff64 = |value: &Complex<F<f64>>| Complex::new(value.re.0, value.im.0);
    let weights = &event.additional_weights.weights;
    let original = to_ff64(
        weights
            .get(&gammalooprs::observables::events::AdditionalWeightKey::Original)
            .expect("amplitude events should expose the bare graph contribution"),
    );
    let full_factor = to_ff64(
        weights
            .get(&gammalooprs::observables::events::AdditionalWeightKey::FullMultiplicativeFactor)
            .expect("amplitude events should record the full multiplicative factor"),
    );
    let threshold_counterterms = weights
        .iter()
        .filter_map(|(key, value)| match key {
            gammalooprs::observables::events::AdditionalWeightKey::ThresholdCounterterm {
                subset_index,
            } => Some((*subset_index, to_ff64(value))),
            _ => None,
        })
        .sorted_by_key(|(subset_index, _)| *subset_index)
        .collect_vec();
    assert!(
        !threshold_counterterms.is_empty(),
        "physical gg->hhh amplitude events should expose threshold-counterterm entries",
    );
    for (expected_index, (subset_index, _)) in threshold_counterterms.iter().enumerate() {
        assert_eq!(*subset_index, expected_index);
    }

    let zero = Complex::new(0.0, 0.0);
    let threshold_sum = threshold_counterterms
        .iter()
        .fold(zero, |acc, (_, value)| acc + *value);
    let reconstructed_event_weight = (original + threshold_sum) * full_factor;
    let event_weight = to_ff64(&event.weight);
    assert!(
        complex_distance(event_weight, reconstructed_event_weight) < 1.0e-11,
        "event weight should match Original + sum(ThresholdCounterterm) times FullMultiplicativeFactor: {} vs {}",
        event_weight,
        reconstructed_event_weight,
    );

    let full_result_factor = evaluation
        .parameterization_jacobian
        .expect("x-space sample evaluation should surface a parameterization Jacobian")
        .0
        * evaluation.integrator_weight.0;
    let expected_total = Complex::new(
        evaluation.integrand_result.re.0 * full_result_factor,
        evaluation.integrand_result.im.0 * full_result_factor,
    );
    let summed_event_weights = evaluation
        .event_groups
        .iter()
        .flat_map(|event_group| event_group.iter())
        .fold(zero, |acc, event| acc + to_ff64(&event.weight));
    assert!(
        complex_distance(summed_event_weights, expected_total) < 1.0e-11,
        "summed amplitude event weights should reproduce the returned amplitude after the full-factor application: {} vs {}",
        summed_event_weights,
        expected_total,
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn amplitude_selectors_generate_internal_events_without_surfacing_them() -> Result<()> {
    let test_name = "amplitude_selectors_generate_internal_events_without_surfacing_them";
    let mut cli = setup_gg_hhh_threshold_amplitude_cli(test_name)?;
    cli.run_command(
        r#"set process -p gg_hhh -i 1L kv general.generate_events=false general.store_additional_weights_in_event=false"#,
    )?;
    cli.run_command(
        r#"set process -p gg_hhh -i 1L string '
[quantities.higgs_count]
type = "particle"
computation = "count"
pdgs = [25]

[selectors.reject_hhh]
quantity = "higgs_count"
selector = "discrete_range"
min = 4
max = 4
'"#,
    )?;

    let process_ref = ProcessRef::Unqualified("gg_hhh".to_string());
    let integrand_name = "1L".to_string();
    let (process_id, resolved_integrand_name) = cli
        .state
        .find_integrand_ref(Some(&process_ref), Some(&integrand_name))?;

    let points = arr2(&[[0.23, 0.41, 0.67]]);
    let discrete_dims = arr2(&[[0usize, 0usize, 0usize]]);
    let result = evaluate_sample(
        &mut cli.state,
        &EvaluateSamples {
            process_id: Some(process_id),
            integrand_name: Some(resolved_integrand_name),
            use_arb_prec: false,
            minimal_output: false,
            return_generated_events: None,
            momentum_space: false,
            points: points.view(),
            integrator_weights: None,
            discrete_dims: Some(discrete_dims.view()),
            graph_names: None,
            orientations: None,
        },
    )?;

    let evaluation = result.sample.evaluation;
    let metadata = evaluation
        .evaluation_metadata
        .as_ref()
        .expect("amplitude sample output should include evaluation metadata");
    assert!(evaluation.event_groups.is_empty());
    assert_eq!(evaluation.integrand_result.re.0, 0.0);
    assert_eq!(evaluation.integrand_result.im.0, 0.0);
    assert!(
        metadata.generated_event_count > 0,
        "selectors should force internal event generation for amplitudes",
    );
    assert_eq!(metadata.accepted_event_count, 0);
    assert!(result.observables.histograms.is_empty());

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn amplitude_observables_accumulate_without_returning_events() -> Result<()> {
    let test_name = "amplitude_observables_accumulate_without_returning_events";
    let mut cli = setup_gg_hhh_threshold_amplitude_cli(test_name)?;
    cli.run_command(
        r#"set process -p gg_hhh -i 1L kv general.generate_events=false general.store_additional_weights_in_event=false"#,
    )?;
    cli.run_command(
        r#"set process -p gg_hhh -i 1L string '
[quantities.integral]
type = "integral"

[quantities.graph_id]
type = "graph_id"

[quantities.orientation_id]
type = "orientation_id"

[quantities.lmb_channel_id]
type = "lmb_channel_id"

[observables.integral_hist]
quantity = "integral"
kind = "discrete"
domain = { type = "single_bin" }

[observables.graph_id_hist]
quantity = "graph_id"
kind = "discrete"
domain = { type = "explicit_range", min = 0, max = 64 }

[observables.orientation_id_hist]
quantity = "orientation_id"
kind = "discrete"
domain = { type = "explicit_range", min = 0, max = 64 }

[observables.lmb_channel_id_hist]
quantity = "lmb_channel_id"
kind = "discrete"
domain = { type = "explicit_range", min = 0, max = 64 }
'"#,
    )?;

    let process_ref = ProcessRef::Unqualified("gg_hhh".to_string());
    let integrand_name = "1L".to_string();
    let (process_id, resolved_integrand_name) = cli
        .state
        .find_integrand_ref(Some(&process_ref), Some(&integrand_name))?;

    let points = arr2(&[[0.23, 0.41, 0.67]]);
    let discrete_dims = arr2(&[[0usize, 0usize, 0usize]]);
    let result = evaluate_sample(
        &mut cli.state,
        &EvaluateSamples {
            process_id: Some(process_id),
            integrand_name: Some(resolved_integrand_name),
            use_arb_prec: false,
            minimal_output: false,
            return_generated_events: None,
            momentum_space: false,
            points: points.view(),
            integrator_weights: None,
            discrete_dims: Some(discrete_dims.view()),
            graph_names: None,
            orientations: None,
        },
    )?;

    let evaluation = result.sample.evaluation;
    let metadata = evaluation
        .evaluation_metadata
        .as_ref()
        .expect("amplitude sample output should include evaluation metadata");
    assert!(evaluation.event_groups.is_empty());
    assert!(
        metadata.generated_event_count > 0,
        "observables should force internal event generation for amplitudes",
    );
    assert!(
        metadata.accepted_event_count > 0,
        "observable processing should accept at least one amplitude event",
    );

    let histogram_names = result.observables.histograms.keys().cloned().collect_vec();
    assert_eq!(
        histogram_names,
        vec![
            "graph_id_hist".to_string(),
            "integral_hist".to_string(),
            "lmb_channel_id_hist".to_string(),
            "orientation_id_hist".to_string(),
        ]
    );
    for histogram in result.observables.histograms.values() {
        assert_eq!(histogram.sample_count, 1);
    }

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
