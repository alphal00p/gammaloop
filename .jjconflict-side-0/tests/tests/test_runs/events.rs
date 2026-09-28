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
    let amplitude_threshold_counterterms = weights
        .iter()
        .filter_map(|(key, value)| {
            match key {
            gammalooprs::observables::events::AdditionalWeightKey::AmplitudeThresholdCounterterm {
                esurface_id,
                overlap_group,
            } => Some(((*esurface_id, *overlap_group), to_ff64(value))),
            _ => None,
        }
        })
        .sorted_by_key(|((esurface_id, overlap_group), _)| (*esurface_id, *overlap_group))
        .collect_vec();
    assert!(
        !amplitude_threshold_counterterms.is_empty(),
        "physical gg->hhh amplitude events should expose threshold-counterterm entries",
    );
    let threshold_counterterm_labels = amplitude_threshold_counterterms
        .iter()
        .map(|((esurface_id, overlap_group), _)| (*esurface_id, *overlap_group))
        .collect_vec();
    assert_eq!(
        threshold_counterterm_labels.len(),
        threshold_counterterm_labels
            .iter()
            .copied()
            .unique()
            .count(),
        "amplitude threshold-counterterm weights should have unique (raised esurface id, overlap group) labels",
    );

    let zero = Complex::new(0.0, 0.0);
    let threshold_sum = amplitude_threshold_counterterms
        .iter()
        .fold(zero, |acc, (_, value)| acc + *value);
    let reconstructed_event_weight = (original + threshold_sum) * full_factor;
    let event_weight = to_ff64(&event.weight);
    assert!(
        complex_distance(event_weight, reconstructed_event_weight) < 1.0e-11,
        "event weight should match Original + sum(AmplitudeThresholdCounterterm) times FullMultiplicativeFactor: {} vs {}",
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
fn amplitude_threshold_counterterms_follow_lmb_channel_normalization() -> Result<()> {
    let test_name = "amplitude_threshold_counterterms_follow_lmb_channel_normalization";
    let mut cli = setup_gg_hhh_threshold_amplitude_cli(test_name)?;
    cli.run_command(
        "set process kv general.generate_events=true general.store_additional_weights_in_event=true stability.rotation_axis=[]",
    )?;

    let process_ref = ProcessRef::Unqualified("gg_hhh".to_string());
    let integrand_name = "1L".to_string();
    let (process_id, resolved_integrand_name) = cli
        .state
        .find_integrand_ref(Some(&process_ref), Some(&integrand_name))?;
    let info = cli
        .state
        .get_integrand_info(Some(&process_ref), Some(&integrand_name))?;
    let group = info
        .graph_groups
        .iter()
        .find(|group| {
            !group.loop_momentum_bases.is_empty() && !group.threshold_esurface_ids.is_empty()
        })
        .expect("expected an amplitude graph group with generated threshold counterterms");
    let group_id = group.group_id;
    let graph_name = group
        .graphs
        .iter()
        .find(|graph| graph.is_master)
        .expect("expected a master amplitude graph")
        .name
        .clone();
    let basis_ids = group
        .loop_momentum_bases
        .iter()
        .map(|basis| basis.basis_id)
        .take(2)
        .collect_vec();
    assert_eq!(
        basis_ids.len(),
        2,
        "amplitude regression fixture must provide two distinct LMBs"
    );

    let evaluate = |cli: &mut gammaloop_integration_tests::CLIState,
                    discrete_values: &[usize]|
     -> Result<gammalooprs::integrands::evaluation::SampleEvaluationResult> {
        let points = arr2(&[[0.23, 0.41, 0.67]]);
        let discrete_dims =
            ndarray::Array2::from_shape_vec((1, discrete_values.len()), discrete_values.to_vec())?;
        Ok(evaluate_sample(
            &mut cli.state,
            &EvaluateSamples {
                process_id: Some(process_id),
                integrand_name: Some(resolved_integrand_name.clone()),
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
        )?
        .sample)
    };
    let summarize = |result: &gammalooprs::integrands::evaluation::SampleEvaluationResult| {
        let mut event_weight = (0.0, 0.0);
        let mut original_weight = (0.0, 0.0);
        let mut threshold_weight = (0.0, 0.0);
        let mut threshold_magnitude = 0.0;
        for event in result
            .evaluation
            .event_groups
            .iter()
            .flat_map(|group| group.iter())
        {
            event_weight.0 += event.weight.re.0;
            event_weight.1 += event.weight.im.0;
            for (key, value) in &event.additional_weights.weights {
                match key {
                    gammalooprs::observables::events::AdditionalWeightKey::Original => {
                        original_weight.0 += value.re.0;
                        original_weight.1 += value.im.0;
                    }
                    gammalooprs::observables::events::AdditionalWeightKey::AmplitudeThresholdCounterterm {
                        ..
                    } => {
                        threshold_weight.0 += value.re.0;
                        threshold_weight.1 += value.im.0;
                        threshold_magnitude += value.re.0.abs() + value.im.0.abs();
                    }
                    _ => {}
                }
            }
        }
        (
            (
                result.evaluation.integrand_result.re.0,
                result.evaluation.integrand_result.im.0,
            ),
            event_weight,
            original_weight,
            threshold_weight,
            threshold_magnitude,
        )
    };
    let assert_pair_close = |actual: (f64, f64), expected: (f64, f64), label: &str| {
        let re_scale = expected.0.abs().max(actual.0.abs()).max(1.0e-30);
        let im_scale = expected.1.abs().max(actual.1.abs()).max(1.0e-30);
        assert!(
            (actual.0 - expected.0).abs() <= 1.0e-9 * re_scale,
            "{label} real mismatch: got {}, expected {}",
            actual.0,
            expected.0,
        );
        assert!(
            (actual.1 - expected.1).abs() <= 1.0e-9 * im_scale,
            "{label} imaginary mismatch: got {}, expected {}",
            actual.1,
            expected.1,
        );
    };

    let mut baseline_summaries = Vec::new();
    for basis_id in &basis_ids {
        cli.run_command(&format!(
            r#"set process string '
[sampling]
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
lmb_basis_ids = {{ "{graph_name}" = [{basis_id}] }}
'"#,
        ))?;
        baseline_summaries.push(summarize(&evaluate(&mut cli, &[group_id])?));
    }
    assert!(
        baseline_summaries.iter().any(|summary| summary.4 > 0.0),
        "amplitude regression fixture must evaluate a nonzero threshold counterterm"
    );

    for channel_weight in ["ose", "inverse_jacobian"] {
        cli.run_command(&format!(
            r#"set process string '
[sampling]
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = true
lmb_channels = "monte_carlo"
lmb_channel_weight = "{channel_weight}"
lmb_basis_ids = {{ "{graph_name}" = [{}, {}] }}
'"#,
            basis_ids[0], basis_ids[1],
        ))?;
        let mut monte_carlo_sum = ((0.0, 0.0), (0.0, 0.0), (0.0, 0.0), (0.0, 0.0));
        for (channel_id, baseline_summary) in baseline_summaries.iter().copied().enumerate() {
            let channel_summary = summarize(&evaluate(&mut cli, &[group_id, channel_id])?);
            let factor = if baseline_summary.2.0.abs() >= baseline_summary.2.1.abs() {
                assert!(baseline_summary.2.0.abs() > 1.0e-30);
                channel_summary.2.0 / baseline_summary.2.0
            } else {
                assert!(baseline_summary.2.1.abs() > 1.0e-30);
                channel_summary.2.1 / baseline_summary.2.1
            };
            assert!(
                factor.is_finite() && factor > 0.0 && factor < 1.0,
                "{channel_weight} amplitude channel {channel_id} produced invalid LMB weight {factor}"
            );
            let scale_pair = |pair: (f64, f64)| (pair.0 * factor, pair.1 * factor);
            assert_pair_close(
                channel_summary.0,
                scale_pair(baseline_summary.0),
                "amplitude channel integrand",
            );
            assert_pair_close(
                channel_summary.1,
                scale_pair(baseline_summary.1),
                "amplitude channel event weight",
            );
            assert_pair_close(
                channel_summary.2,
                scale_pair(baseline_summary.2),
                "amplitude channel original weight",
            );
            assert_pair_close(
                channel_summary.3,
                scale_pair(baseline_summary.3),
                "amplitude channel threshold-counterterm weight",
            );
            let expected_threshold_magnitude = baseline_summary.4 * factor;
            let threshold_scale = expected_threshold_magnitude
                .abs()
                .max(channel_summary.4.abs())
                .max(1.0e-30);
            assert!(
                (channel_summary.4 - expected_threshold_magnitude).abs()
                    <= 1.0e-9 * threshold_scale,
                "{channel_weight} amplitude channel {channel_id} threshold-counterterm magnitudes were not scaled by the LMB weight"
            );
            monte_carlo_sum.0.0 += channel_summary.0.0;
            monte_carlo_sum.0.1 += channel_summary.0.1;
            monte_carlo_sum.1.0 += channel_summary.1.0;
            monte_carlo_sum.1.1 += channel_summary.1.1;
            monte_carlo_sum.2.0 += channel_summary.2.0;
            monte_carlo_sum.2.1 += channel_summary.2.1;
            monte_carlo_sum.3.0 += channel_summary.3.0;
            monte_carlo_sum.3.1 += channel_summary.3.1;
        }

        cli.run_command(&format!(
            r#"set process string '
[sampling]
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = true
lmb_channels = "summed"
lmb_channel_weight = "{channel_weight}"
lmb_basis_ids = {{ "{graph_name}" = [{}, {}] }}
'"#,
            basis_ids[0], basis_ids[1],
        ))?;
        let summed_summary = summarize(&evaluate(&mut cli, &[group_id])?);
        assert_pair_close(
            summed_summary.0,
            monte_carlo_sum.0,
            "summed versus Monte Carlo amplitude",
        );
        assert_pair_close(
            summed_summary.1,
            monte_carlo_sum.1,
            "summed versus Monte Carlo amplitude event weight",
        );
        assert_pair_close(
            summed_summary.2,
            monte_carlo_sum.2,
            "summed versus Monte Carlo amplitude original weight",
        );
        assert_pair_close(
            summed_summary.3,
            monte_carlo_sum.3,
            "summed versus Monte Carlo amplitude threshold-counterterm weight",
        );
    }

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
