use std::{collections::BTreeMap, time::Duration};

use color_eyre::Result;
use gammaloop_api::{
    commands::evaluate_samples::{EvaluateSamples, EvaluateSamplesPrecise},
    integrand_info::IntegrandKind,
};
use gammaloop_integration_tests::{
    CLIState, default_momentum_space_point, default_xspace_point, setup_sm_differential_lu_cli,
};
use ndarray::{Array1, Array2};
use serial_test::serial;

fn configure_jet_quantities(cli: &mut CLIState) -> Result<()> {
    configure_jet_quantities_with_clustered_pdgs(cli, None)
}

fn configure_jet_quantities_with_clustered_pdgs(
    cli: &mut CLIState,
    clustered_pdgs: Option<&str>,
) -> Result<()> {
    let clustered_pdgs = clustered_pdgs
        .map(|clustered_pdgs| format!("clustered_pdgs = {clustered_pdgs}\n"))
        .unwrap_or_default();
    cli.run_command(&format!(
        r#"set process string '
[quantities.leading_jet_pt]
type = "jet"
quantity = "PT"
dR = 0.4
{clustered_pdgs}

[quantities.jet_count]
type = "jet"
computation = "count"
dR = 0.4
{clustered_pdgs}
'"#,
    ))
}

fn add_jet_observables(cli: &mut CLIState) -> Result<()> {
    cli.run_command(
        r#"set process string '
[observables.leading_jet_pt_hist]
quantity = "leading_jet_pt"
entry_selection = "leading_only"
x_min = 0.0
x_max = 1000.0
n_bins = 8

[observables.jet_count_hist]
quantity = "jet_count"
entry_selection = "all"
kind = "discrete"
domain = { type = "explicit_range", min = 0, max = 5 }
'"#,
    )
}

fn add_leading_jet_selector(cli: &mut CLIState, min: f64) -> Result<()> {
    cli.run_command(&format!(
        "set process string '\n[selectors.leading_jet_pt_cut]\nquantity = \"leading_jet_pt\"\nselector = \"value_range\"\nentry_selection = \"leading_only\"\nmin = {min}\n'"
    ))
}

fn update_leading_jet_selector(cli: &mut CLIState, min: f64) -> Result<()> {
    cli.run_command(&format!(
        "set process kv selectors.leading_jet_pt_cut.min={min}"
    ))
}

fn evaluate_x_samples(
    cli: &mut CLIState,
    points: &[Vec<f64>],
) -> Result<gammalooprs::integrands::evaluation::BatchSampleEvaluationResult> {
    evaluate_x_samples_with_weights(cli, points, None)
}

fn evaluate_x_samples_with_weights(
    cli: &mut CLIState,
    points: &[Vec<f64>],
    integrator_weights: Option<&[f64]>,
) -> Result<gammalooprs::integrands::evaluation::BatchSampleEvaluationResult> {
    evaluate_x_samples_with_weights_and_discrete_dims(cli, points, integrator_weights, None)
}

fn evaluate_x_samples_with_discrete_dims(
    cli: &mut CLIState,
    points: &[Vec<f64>],
    discrete_dims: &[Vec<usize>],
) -> Result<gammalooprs::integrands::evaluation::BatchSampleEvaluationResult> {
    evaluate_x_samples_with_weights_and_discrete_dims(cli, points, None, Some(discrete_dims))
}

fn evaluate_x_samples_with_weights_and_discrete_dims(
    cli: &mut CLIState,
    points: &[Vec<f64>],
    integrator_weights: Option<&[f64]>,
    discrete_dims: Option<&[Vec<usize>]>,
) -> Result<gammalooprs::integrands::evaluation::BatchSampleEvaluationResult> {
    let ncols = points.first().map(Vec::len).unwrap_or(0);
    let flat = points
        .iter()
        .flat_map(|point| point.iter().copied())
        .collect::<Vec<_>>();
    let array = Array2::from_shape_vec((points.len(), ncols), flat)?;
    let integrator_weights = integrator_weights.map(|weights| Array1::from_vec(weights.to_vec()));
    let discrete_dims = discrete_dims
        .map(|dims| {
            let ncols = dims.first().map(Vec::len).unwrap_or(0);
            let flat = dims
                .iter()
                .flat_map(|dim| dim.iter().copied())
                .collect::<Vec<_>>();
            Array2::from_shape_vec((dims.len(), ncols), flat)
        })
        .transpose()?;
    EvaluateSamples {
        process_id: None,
        integrand_name: None,
        use_arb_prec: false,
        minimal_output: false,
        return_generated_events: None,
        momentum_space: false,
        points: array.view(),
        integrator_weights: integrator_weights.as_ref().map(|weights| weights.view()),
        discrete_dims: discrete_dims.as_ref().map(|dims| dims.view()),
        graph_names: None,
        orientations: None,
    }
    .run(&mut cli.state)
}

fn evaluate_momentum_sample(
    cli: &mut CLIState,
    point: &[f64],
) -> Result<gammalooprs::integrands::evaluation::SingleSampleEvaluationResult> {
    evaluate_momentum_sample_with_weight(cli, point, None)
}

fn evaluate_momentum_sample_with_weight(
    cli: &mut CLIState,
    point: &[f64],
    integrator_weight: Option<f64>,
) -> Result<gammalooprs::integrands::evaluation::SingleSampleEvaluationResult> {
    let array = Array2::from_shape_vec((1, point.len()), point.to_vec())?;
    let integrator_weights = integrator_weight.map(|weight| Array1::from_vec(vec![weight]));
    let results = EvaluateSamples {
        process_id: None,
        integrand_name: None,
        use_arb_prec: false,
        minimal_output: false,
        return_generated_events: None,
        momentum_space: true,
        points: array.view(),
        integrator_weights: integrator_weights.as_ref().map(|weights| weights.view()),
        discrete_dims: None,
        graph_names: None,
        orientations: None,
    }
    .run(&mut cli.state)?;
    let mut samples = results.samples;
    Ok(
        gammalooprs::integrands::evaluation::SingleSampleEvaluationResult {
            sample: samples.remove(0),
            observables: results.observables,
        },
    )
}

fn evaluate_x_samples_minimal(
    cli: &mut CLIState,
    points: &[Vec<f64>],
) -> Result<gammalooprs::integrands::evaluation::BatchSampleEvaluationResult> {
    let ncols = points.first().map(Vec::len).unwrap_or(0);
    let flat = points
        .iter()
        .flat_map(|point| point.iter().copied())
        .collect::<Vec<_>>();
    let array = Array2::from_shape_vec((points.len(), ncols), flat)?;
    EvaluateSamples {
        process_id: None,
        integrand_name: None,
        use_arb_prec: false,
        minimal_output: true,
        return_generated_events: None,
        momentum_space: false,
        points: array.view(),
        integrator_weights: None,
        discrete_dims: None,
        graph_names: None,
        orientations: None,
    }
    .run(&mut cli.state)
}

fn evaluate_x_samples_precise(
    cli: &mut CLIState,
    points: &[Vec<f64>],
) -> Result<gammalooprs::integrands::evaluation::PreciseBatchSampleEvaluationResult> {
    let ncols = points.first().map(Vec::len).unwrap_or(0);
    let flat = points
        .iter()
        .flat_map(|point| point.iter().copied())
        .collect::<Vec<_>>();
    let array = Array2::from_shape_vec((points.len(), ncols), flat)?;
    EvaluateSamplesPrecise {
        process_id: None,
        integrand_name: None,
        use_arb_prec: true,
        minimal_output: false,
        return_generated_events: None,
        momentum_space: false,
        points: array.view(),
        integrator_weights: None,
        discrete_dims: None,
        graph_names: None,
        orientations: None,
    }
    .run(&mut cli.state)
}

fn assert_complex_scaled(
    scaled: &spenso::algebra::complex::Complex<gammalooprs::utils::F<f64>>,
    baseline: &spenso::algebra::complex::Complex<gammalooprs::utils::F<f64>>,
    factor: f64,
) {
    let expected_re = baseline.re.0 * factor;
    let expected_im = baseline.im.0 * factor;
    let re_scale = expected_re.abs().max(1.0);
    let im_scale = expected_im.abs().max(1.0);
    assert!(
        (scaled.re.0 - expected_re).abs() <= 1.0e-12 * re_scale,
        "real part mismatch: got {}, expected {}",
        scaled.re.0,
        expected_re
    );
    assert!(
        (scaled.im.0 - expected_im).abs() <= 1.0e-12 * im_scale,
        "imaginary part mismatch: got {}, expected {}",
        scaled.im.0,
        expected_im
    );
}

fn histogram_total_sum_weights(histogram: &gammalooprs::observables::HistogramSnapshot) -> f64 {
    histogram
        .bins
        .iter()
        .map(|bin| bin.sum_weights)
        .sum::<f64>()
        + histogram.underflow_bin.sum_weights
        + histogram.overflow_bin.sum_weights
}

fn metadata(
    result: &gammalooprs::integrands::evaluation::SampleEvaluationResult,
) -> &gammalooprs::integrands::evaluation::EvaluationMetaData {
    result
        .evaluation
        .evaluation_metadata
        .as_ref()
        .expect("full evaluate_sample(s) output should include evaluation metadata")
}

fn some_additional_weights_are_present(
    result: &gammalooprs::integrands::evaluation::SampleEvaluationResult,
) -> bool {
    result
        .evaluation
        .event_groups
        .iter()
        .flat_map(|group| group.iter())
        .any(|event| !event.additional_weights.weights.is_empty())
}

fn momentum_signature(
    momentum: &gammalooprs::momentum::FourMomentum<gammalooprs::utils::F<f64>>,
) -> String {
    format!(
        "{:.17e},{:.17e},{:.17e},{:.17e}",
        momentum.temporal.value.0,
        momentum.spatial.px.0,
        momentum.spatial.py.0,
        momentum.spatial.pz.0
    )
}

fn event_signature(event: &gammalooprs::observables::Event) -> String {
    let incoming = event
        .kinematic_configuration
        .0
        .iter()
        .map(momentum_signature)
        .collect::<Vec<_>>()
        .join("|");
    let outgoing = event
        .kinematic_configuration
        .1
        .iter()
        .map(momentum_signature)
        .collect::<Vec<_>>()
        .join("|");
    let lmb_channel = event
        .cut_info
        .lmb_channel_edge_ids
        .as_ref()
        .map(|edge_ids| {
            edge_ids
                .iter()
                .map(|edge_id| edge_id.to_string())
                .collect::<Vec<_>>()
                .join(",")
        })
        .unwrap_or_else(|| "none".to_string());

    format!(
        "graph={} cut={} lmb={} weight={:.17e},{:.17e} in={} out={}",
        event.cut_info.graph_id,
        event.cut_info.cut_id,
        lmb_channel,
        event.weight.re.0,
        event.weight.im.0,
        incoming,
        outgoing
    )
}

fn event_multiset(events: impl IntoIterator<Item = String>) -> BTreeMap<String, usize> {
    let mut counts = BTreeMap::new();
    for event in events {
        *counts.entry(event).or_insert(0) += 1;
    }
    counts
}

#[test]
#[serial]
fn lu_rust_get_integrand_info_reports_groups_orientations_lmbs_and_cuts() -> Result<()> {
    let cli = setup_sm_differential_lu_cli("lu_rust_get_integrand_info")?;
    let (process_id, integrand_name) = cli.state.find_integrand_ref(None, None)?;
    let info = cli.state.get_integrand_info(None, None)?;
    let process = &cli.state.process_list.processes[process_id];

    assert_eq!(info.process_id, process.definition.process_id);
    assert_eq!(info.process_name, process.definition.folder_name);
    assert_eq!(info.integrand_name, integrand_name);
    assert_eq!(info.kind, IntegrandKind::CrossSection);
    assert_eq!(info.graph_group_count, info.graph_groups.len());
    assert!(info.graph_count >= info.graph_group_count);
    assert!(info.record_size_bytes > 0);

    for group in &info.graph_groups {
        assert_eq!(
            group.graphs.iter().filter(|graph| graph.is_master).count(),
            1
        );
        assert!(!group.graphs.is_empty());

        assert!(!group.orientation_edge_ids.is_empty());

        for (expected_id, orientation) in group.orientations.iter().enumerate() {
            assert_eq!(orientation.orientation_id, expected_id);
            assert_eq!(
                orientation.signature.len(),
                group.orientation_edge_ids.len()
            );
            assert!(
                orientation
                    .signature
                    .iter()
                    .all(|sign| matches!(*sign, -1..=1))
            );
        }

        let mut seen_channel_ids = std::collections::BTreeSet::new();
        for (expected_id, basis) in group.loop_momentum_bases.iter().enumerate() {
            assert_eq!(basis.basis_id, expected_id);
            assert!(!basis.edge_ids.is_empty());
            if let Some(channel_id) = basis.channel_id {
                assert!(seen_channel_ids.insert(channel_id));
            }
        }
        assert_eq!(
            group
                .loop_momentum_bases
                .iter()
                .filter(|basis| basis.matches_generation_basis)
                .count(),
            1
        );

        for (expected_id, cut) in group.cuts.iter().enumerate() {
            assert_eq!(cut.cut_id, expected_id);
            assert!(!cut.edge_ids.is_empty());
            assert!(cut.raising_power >= 1);
        }
    }

    Ok(())
}

#[test]
#[serial]
fn lu_rust_generated_events_follow_graph_grouping_and_cut_ids() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_rust_generated_event_grouping")?;
    cli.run_command(
        "set process kv general.generate_events=true general.store_additional_weights_in_event=true",
    )?;

    let point = default_xspace_point(&cli)?;
    let mut results = evaluate_x_samples(&mut cli, std::slice::from_ref(&point))?;
    let result = results.samples.remove(0);
    let event_groups = &result.evaluation.event_groups;

    assert_eq!(metadata(&result).generated_event_count, 3);
    assert_eq!(metadata(&result).accepted_event_count, 3);
    assert_eq!(
        event_groups
            .iter()
            .map(|group| group.len())
            .collect::<Vec<_>>(),
        vec![1, 2]
    );

    let first_event = &event_groups[0][0];
    assert_eq!(first_event.cut_info.graph_id, 0);
    assert_eq!(first_event.cut_info.cut_id, 0);
    assert_eq!(
        first_event
            .cut_info
            .particle_pdgs
            .0
            .iter()
            .copied()
            .collect::<Vec<_>>(),
        vec![-11, 11]
    );

    let mut second_group_cut_ids = event_groups[1]
        .iter()
        .map(|event| (event.cut_info.graph_id, event.cut_info.cut_id))
        .collect::<Vec<_>>();
    second_group_cut_ids.sort();
    assert_eq!(second_group_cut_ids, vec![(1, 0), (1, 1)]);
    for event in event_groups[1].iter() {
        assert_eq!(
            event
                .cut_info
                .particle_pdgs
                .0
                .iter()
                .copied()
                .collect::<Vec<_>>(),
            vec![-11, 11]
        );
    }
    let formatted = result.to_string();
    assert!(formatted.contains("Kinematics"));
    assert!(formatted.contains("PDG"));
    assert!(formatted.contains("state"));
    assert!(!formatted.contains("incoming PDGs"));
    assert!(formatted.contains("Stability results"));

    Ok(())
}

#[test]
#[serial]
fn lu_rust_precise_evaluate_samples_returns_precision_tagged_results() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_rust_precise_evaluate_samples")?;
    configure_jet_quantities(&mut cli)?;
    add_jet_observables(&mut cli)?;

    let point = default_xspace_point(&cli)?;
    let results = evaluate_x_samples_precise(&mut cli, &[point])?;
    assert_eq!(results.samples.len(), 1);
    assert_eq!(results.observables.histograms.len(), 2);
    match &results.samples[0].evaluation {
        gammalooprs::integrands::evaluation::PreciseEvaluationResultOutput::Double(result) => {
            assert!(result.evaluation_metadata.is_some());
        }
        gammalooprs::integrands::evaluation::PreciseEvaluationResultOutput::Quad(result) => {
            assert!(result.evaluation_metadata.is_some());
        }
        gammalooprs::integrands::evaluation::PreciseEvaluationResultOutput::Arb(result) => {
            assert!(result.evaluation_metadata.is_some());
        }
    }

    Ok(())
}

#[test]
#[serial]
fn lu_rust_evaluate_samples_respect_event_generation_and_observables() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_rust_evaluate_samples_end_to_end")?;
    let point = default_xspace_point(&cli)?;

    let mut results = evaluate_x_samples(&mut cli, std::slice::from_ref(&point))?;
    let result = results.samples.remove(0);
    assert!(result.evaluation.event_groups.is_empty());
    assert!(results.observables.histograms.is_empty());
    assert_eq!(metadata(&result).generated_event_count, 0);
    assert_eq!(metadata(&result).accepted_event_count, 0);
    assert!(result.evaluation.parameterization_jacobian.is_some());

    cli.run_command(
        "set process kv general.generate_events=true general.store_additional_weights_in_event=true",
    )?;
    let mut results = evaluate_x_samples(&mut cli, std::slice::from_ref(&point))?;
    let result = results.samples.remove(0);
    assert!(!result.evaluation.event_groups.is_empty());
    assert!(
        metadata(&result).generated_event_count > 0,
        "expected generated events with general.generate_events=true"
    );
    assert_eq!(
        metadata(&result).accepted_event_count,
        metadata(&result).generated_event_count
    );
    assert!(some_additional_weights_are_present(&result));

    cli.run_command(
        "set process kv general.generate_events=false general.store_additional_weights_in_event=false",
    )?;
    configure_jet_quantities(&mut cli)?;
    add_leading_jet_selector(&mut cli, 1_000_000.0)?;
    let mut results = evaluate_x_samples(&mut cli, std::slice::from_ref(&point))?;
    let result = results.samples.remove(0);
    assert!(result.evaluation.event_groups.is_empty());
    assert!(results.observables.histograms.is_empty());
    assert!(
        metadata(&result).generated_event_count > 0,
        "selector-only mode should still generate temporary events"
    );
    assert_eq!(metadata(&result).accepted_event_count, 0);

    update_leading_jet_selector(&mut cli, 0.0)?;
    let mut results = evaluate_x_samples(&mut cli, std::slice::from_ref(&point))?;
    let result = results.samples.remove(0);
    assert!(result.evaluation.event_groups.is_empty());
    assert!(results.observables.histograms.is_empty());
    assert!(
        metadata(&result).generated_event_count > 0,
        "selector-only mode should still generate temporary events"
    );
    assert!(
        metadata(&result).accepted_event_count > 0,
        "accepted selector-only events should be counted even when not buffered"
    );
    assert_eq!(
        metadata(&result).accepted_event_count,
        metadata(&result).generated_event_count
    );

    add_jet_observables(&mut cli)?;
    let batch = vec![point.clone(), point];
    let results = evaluate_x_samples(&mut cli, &batch)?;
    assert_eq!(results.samples.len(), 2);
    for result in &results.samples {
        assert!(result.evaluation.event_groups.is_empty());
        assert!(
            metadata(result).generated_event_count > 0,
            "selector+observable mode should generate temporary events"
        );
        assert!(
            metadata(result).accepted_event_count > 0,
            "permissive selector should keep at least one event"
        );
        assert!(metadata(result).accepted_event_count <= metadata(result).generated_event_count);
        assert!(metadata(result).event_processing_time.as_nanos() > 0);
    }
    let histogram = results
        .observables
        .histograms
        .get("leading_jet_pt_hist")
        .expect("missing leading-jet histogram");
    assert_eq!(histogram.bins.len(), 8);
    assert_eq!(histogram.sample_count, 2);
    let jet_count_histogram = results
        .observables
        .histograms
        .get("jet_count_hist")
        .expect("missing jet-count histogram");
    assert_eq!(jet_count_histogram.bins.len(), 6);
    assert_eq!(jet_count_histogram.sample_count, 2);
    assert!(!jet_count_histogram.supports_misbinning_mitigation);

    Ok(())
}

#[test]
#[serial]
fn lu_rust_evaluate_samples_can_override_returned_events_without_forcing_internal_generation()
-> Result<()> {
    let mut cli = setup_sm_differential_lu_cli(
        "lu_rust_evaluate_samples_can_override_returned_events_without_forcing_internal_generation",
    )?;
    let point = default_xspace_point(&cli)?;
    let array = Array2::from_shape_vec((1, point.len()), point.clone())?;

    cli.run_command("set process kv general.generate_events=false")?;
    let mut results = EvaluateSamples {
        process_id: None,
        integrand_name: None,
        use_arb_prec: false,
        minimal_output: false,
        return_generated_events: Some(true),
        momentum_space: false,
        points: array.view(),
        integrator_weights: None,
        discrete_dims: None,
        graph_names: None,
        orientations: None,
    }
    .run(&mut cli.state)?;
    let result = results.samples.remove(0);
    assert!(
        !result.evaluation.event_groups.is_empty(),
        "return_generated_events=Some(true) should surface events even when the stored setting is false",
    );
    assert!(metadata(&result).generated_event_count > 0);
    assert_eq!(
        metadata(&result).accepted_event_count,
        metadata(&result).generated_event_count
    );

    cli.run_command("set process kv general.generate_events=true")?;
    let mut results = EvaluateSamples {
        process_id: None,
        integrand_name: None,
        use_arb_prec: false,
        minimal_output: false,
        return_generated_events: Some(false),
        momentum_space: false,
        points: array.view(),
        integrator_weights: None,
        discrete_dims: None,
        graph_names: None,
        orientations: None,
    }
    .run(&mut cli.state)?;
    let result = results.samples.remove(0);
    assert!(
        result.evaluation.event_groups.is_empty(),
        "return_generated_events=Some(false) should suppress surfaced events for this call",
    );
    assert_eq!(metadata(&result).generated_event_count, 0);
    assert_eq!(metadata(&result).accepted_event_count, 0);

    Ok(())
}

#[test]
#[serial]
fn lu_rust_xspace_integrator_weights_align_per_sample_and_scale_events() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_rust_xspace_integrator_weights_per_sample")?;
    cli.run_command(
        "set process kv general.generate_events=true general.store_additional_weights_in_event=true",
    )?;

    let point = default_xspace_point(&cli)?;
    let mut results =
        evaluate_x_samples_with_weights(&mut cli, &[point.clone(), point], Some(&[1.0, 3.0]))?;

    assert_eq!(results.samples.len(), 2);
    let first = results.samples.remove(0);
    let second = results.samples.remove(0);

    assert_eq!(first.evaluation.integrator_weight.0, 1.0);
    assert_eq!(second.evaluation.integrator_weight.0, 3.0);
    assert_eq!(
        first.evaluation.parameterization_jacobian,
        second.evaluation.parameterization_jacobian
    );
    assert_eq!(
        first.evaluation.integrand_result,
        second.evaluation.integrand_result
    );
    assert_eq!(
        first
            .evaluation
            .event_groups
            .iter()
            .map(|group| group.len())
            .collect::<Vec<_>>(),
        second
            .evaluation
            .event_groups
            .iter()
            .map(|group| group.len())
            .collect::<Vec<_>>()
    );

    for (first_group, second_group) in first
        .evaluation
        .event_groups
        .iter()
        .zip(second.evaluation.event_groups.iter())
    {
        for (first_event, second_event) in first_group.iter().zip(second_group.iter()) {
            assert_complex_scaled(&second_event.weight, &first_event.weight, 3.0);
        }
    }

    Ok(())
}

#[test]
#[serial]
fn lu_rust_xspace_integrator_weights_scale_observable_histograms() -> Result<()> {
    let mut unit_cli = setup_sm_differential_lu_cli("lu_rust_xspace_integrator_weights_unit")?;
    unit_cli.run_command("set process kv general.generate_events=true")?;
    configure_jet_quantities(&mut unit_cli)?;
    add_jet_observables(&mut unit_cli)?;

    let point = default_xspace_point(&unit_cli)?;
    let unit_results =
        evaluate_x_samples_with_weights(&mut unit_cli, std::slice::from_ref(&point), Some(&[1.0]))?;
    let unit_histogram = unit_results
        .observables
        .histograms
        .get("leading_jet_pt_hist")
        .expect("missing leading-jet histogram for unit-weight evaluation");

    let mut weighted_cli =
        setup_sm_differential_lu_cli("lu_rust_xspace_integrator_weights_weighted")?;
    weighted_cli.run_command("set process kv general.generate_events=true")?;
    configure_jet_quantities(&mut weighted_cli)?;
    add_jet_observables(&mut weighted_cli)?;

    let weighted_results = evaluate_x_samples_with_weights(
        &mut weighted_cli,
        std::slice::from_ref(&point),
        Some(&[3.0]),
    )?;
    let weighted_histogram = weighted_results
        .observables
        .histograms
        .get("leading_jet_pt_hist")
        .expect("missing leading-jet histogram for weighted evaluation");

    assert_eq!(unit_histogram.sample_count, 1);
    assert_eq!(weighted_histogram.sample_count, 1);
    assert_eq!(unit_histogram.statistics, weighted_histogram.statistics);
    let unit_total = histogram_total_sum_weights(unit_histogram);
    let weighted_total = histogram_total_sum_weights(weighted_histogram);
    let scale = unit_total.abs().max(1.0);
    assert!(
        (weighted_total - 3.0 * unit_total).abs() <= 1.0e-12 * scale,
        "histogram total weight mismatch: got {}, expected {}",
        weighted_total,
        3.0 * unit_total
    );

    Ok(())
}

#[test]
#[serial]
fn lu_rust_default_clustered_pdgs_match_explicit_massless_qcd_list() -> Result<()> {
    let mut default_cli = setup_sm_differential_lu_cli("lu_rust_default_clustered_pdgs")?;
    default_cli.run_command("set process kv general.generate_events=true")?;
    configure_jet_quantities(&mut default_cli)?;
    add_jet_observables(&mut default_cli)?;

    let point = default_xspace_point(&default_cli)?;
    let default_results = evaluate_x_samples(&mut default_cli, std::slice::from_ref(&point))?;

    let mut explicit_cli =
        setup_sm_differential_lu_cli("lu_rust_explicit_clustered_pdgs_default_match")?;
    explicit_cli.run_command("set process kv general.generate_events=true")?;
    configure_jet_quantities_with_clustered_pdgs(
        &mut explicit_cli,
        Some("[-4, -3, -2, -1, 1, 2, 3, 4, 21]"),
    )?;
    add_jet_observables(&mut explicit_cli)?;

    let explicit_results = evaluate_x_samples(&mut explicit_cli, std::slice::from_ref(&point))?;

    assert_eq!(default_results.observables, explicit_results.observables);

    Ok(())
}

#[test]
#[serial]
fn lu_rust_momentum_space_evaluate_sample_reports_no_parameterization_jacobian() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_rust_momentum_space_evaluate_sample")?;
    cli.run_command("set process kv general.generate_events=true")?;
    let point = default_momentum_space_point(&cli)?;
    let result = evaluate_momentum_sample(&mut cli, &point)?;

    assert!(result.sample.evaluation.parameterization_jacobian.is_none());
    assert!(metadata(&result.sample).event_processing_time >= Duration::ZERO);
    assert!(
        result.sample.evaluation.event_groups.len()
            <= metadata(&result.sample).generated_event_count
    );
    let formatted = result.to_string();
    assert!(formatted.contains("parameterization jacobian"));
    assert!(formatted.contains("None"));

    Ok(())
}

#[test]
#[serial]
fn lu_rust_momentum_space_integrator_weights_scale_events_without_jacobian() -> Result<()> {
    let mut cli =
        setup_sm_differential_lu_cli("lu_rust_momentum_space_integrator_weights_scale_events")?;
    cli.run_command("set process kv general.generate_events=true")?;
    let point = default_momentum_space_point(&cli)?;

    let unit = evaluate_momentum_sample_with_weight(&mut cli, &point, Some(1.0))?;
    let weighted = evaluate_momentum_sample_with_weight(&mut cli, &point, Some(2.5))?;

    assert!(unit.sample.evaluation.parameterization_jacobian.is_none());
    assert!(
        weighted
            .sample
            .evaluation
            .parameterization_jacobian
            .is_none()
    );
    assert_eq!(unit.sample.evaluation.integrator_weight.0, 1.0);
    assert_eq!(weighted.sample.evaluation.integrator_weight.0, 2.5);
    assert_eq!(
        unit.sample.evaluation.integrand_result,
        weighted.sample.evaluation.integrand_result
    );
    assert_eq!(
        unit.sample
            .evaluation
            .event_groups
            .iter()
            .map(|group| group.len())
            .collect::<Vec<_>>(),
        weighted
            .sample
            .evaluation
            .event_groups
            .iter()
            .map(|group| group.len())
            .collect::<Vec<_>>()
    );

    for (unit_group, weighted_group) in unit
        .sample
        .evaluation
        .event_groups
        .iter()
        .zip(weighted.sample.evaluation.event_groups.iter())
    {
        for (unit_event, weighted_event) in unit_group.iter().zip(weighted_group.iter()) {
            assert_complex_scaled(&weighted_event.weight, &unit_event.weight, 2.5);
        }
    }

    Ok(())
}

#[test]
#[serial]
fn lu_rust_explicit_lmb_multichanneling_groups_channel_events_and_tags_metadata() -> Result<()> {
    let point_group_label = "lu_rust_explicit_lmb_multichanneling_explicit";
    let mut explicit_cli = setup_sm_differential_lu_cli(point_group_label)?;
    explicit_cli.run_command(
        "set process kv general.generate_events=true general.store_additional_weights_in_event=true",
    )?;
    explicit_cli.run_command(
        r#"set process string '
[sampling]
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = true
lmb_channels = "summed"
'"#,
    )?;

    let (process_id, integrand_name) = explicit_cli.state.find_integrand_ref(None, None)?;
    let integrand = explicit_cli
        .state
        .process_list
        .get_integrand(process_id, &integrand_name)?
        .require_generated()?;
    let group_id = (0..integrand.graph_count())
        .filter_map(|graph_id| integrand.graph_name_by_id(graph_id))
        .filter_map(|graph_name| integrand.resolve_group_id_by_master_name(graph_name).ok())
        .find(|&group_id| integrand.group_channel_count(group_id).unwrap_or(0) > 1)
        .expect("expected at least one graph group with more than one LMB channel");
    let channel_count = integrand.group_channel_count(group_id).unwrap();
    let point = default_xspace_point(&explicit_cli)?;

    let mut explicit_results = evaluate_x_samples_with_discrete_dims(
        &mut explicit_cli,
        std::slice::from_ref(&point),
        &[vec![group_id.0]],
    )?;
    let explicit = explicit_results.samples.remove(0);
    assert_eq!(explicit.evaluation.event_groups.len(), 1);
    let explicit_events = &explicit.evaluation.event_groups[0];
    assert_eq!(
        metadata(&explicit).generated_event_count,
        metadata(&explicit).accepted_event_count
    );
    assert_eq!(
        metadata(&explicit).generated_event_count,
        explicit_events.len(),
        "expected one retained event per generated accepted event for explicit LMB summation"
    );
    let explicit_channel_tags = explicit_events
        .iter()
        .map(|event| {
            event
                .cut_info
                .lmb_channel_edge_ids
                .as_ref()
                .expect("explicit multi-channeling events should carry LMB channel metadata")
                .iter()
                .copied()
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let unique_channels = explicit_channel_tags
        .iter()
        .cloned()
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(
        unique_channels.len(),
        channel_count,
        "explicit multi-channeling should contribute events from every LMB channel into the same event group"
    );

    let mut per_channel_cli =
        setup_sm_differential_lu_cli("lu_rust_explicit_lmb_multichanneling_discrete_channels")?;
    per_channel_cli.run_command(
        "set process kv general.generate_events=true general.store_additional_weights_in_event=true",
    )?;
    per_channel_cli.run_command(
        r#"set process string '
[sampling]
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = true
lmb_channels = "monte_carlo"
'"#,
    )?;

    let mut discrete_event_signatures = Vec::new();
    for channel in 0..channel_count {
        let mut results = evaluate_x_samples_with_discrete_dims(
            &mut per_channel_cli,
            std::slice::from_ref(&point),
            &[vec![group_id.0, channel]],
        )?;
        let result = results.samples.remove(0);
        for group in result.evaluation.event_groups.iter() {
            for event in group.iter() {
                let lmb_channel_edge_ids = event
                    .cut_info
                    .lmb_channel_edge_ids
                    .as_ref()
                    .expect("discrete multi-channeling events should carry LMB channel metadata")
                    .iter()
                    .copied()
                    .collect::<Vec<_>>();
                assert!(
                    unique_channels.contains(&lmb_channel_edge_ids),
                    "channel tag should correspond to one explicit-sum LMB basis"
                );
                discrete_event_signatures.push(event_signature(event));
            }
        }
    }

    assert_eq!(
        event_multiset(explicit_events.iter().map(event_signature)),
        event_multiset(discrete_event_signatures),
        "explicit multi-channeling should return the union of the per-channel event sets, grouped into one graph-group event group"
    );

    Ok(())
}

#[test]
#[serial]
fn lu_rust_minimal_output_keeps_events_and_observables() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_rust_minimal_output")?;
    cli.run_command(
        "set process kv general.generate_events=true general.store_additional_weights_in_event=true",
    )?;
    configure_jet_quantities(&mut cli)?;
    add_jet_observables(&mut cli)?;

    let point = default_xspace_point(&cli)?;
    let mut results = evaluate_x_samples_minimal(&mut cli, std::slice::from_ref(&point))?;
    let result = results.samples.remove(0);

    assert!(result.evaluation.evaluation_metadata.is_none());
    assert!(!result.evaluation.event_groups.is_empty());
    assert!(
        results
            .observables
            .histograms
            .contains_key("leading_jet_pt_hist")
    );
    assert!(
        results
            .observables
            .histograms
            .contains_key("jet_count_hist")
    );

    let formatted = gammalooprs::integrands::evaluation::SingleSampleEvaluationResult {
        sample: result,
        observables: results.observables,
    }
    .to_string();
    assert!(!formatted.contains("Evaluation metadata"));
    assert!(!formatted.contains("Stability results"));
    assert!(formatted.contains("Observable snapshots"));

    Ok(())
}
