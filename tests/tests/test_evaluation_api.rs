use std::time::Duration;

use color_eyre::Result;
use gammaloop_api::commands::evaluate_samples::{EvaluateSamples, EvaluateSamplesPrecise};
use gammaloop_integration_tests::{
    CLIState, default_momentum_space_point, default_xspace_point, enable_lu_e2e_hack,
    setup_sm_differential_lu_cli,
};
use ndarray::Array2;
use serial_test::serial;

fn configure_leading_jet_quantity(cli: &mut CLIState) -> Result<()> {
    cli.run_command(
        r#"set process string '
[quantities.leading_jet_pt]
type = "jet_pt"
dR = 0.4
'"#,
    )
}

fn add_leading_jet_observable(cli: &mut CLIState) -> Result<()> {
    cli.run_command(
        r#"set process string '
[observables.leading_jet_pt_hist]
quantity = "leading_jet_pt"
entry_selection = "leading_only"
x_min = 0.0
x_max = 1000.0
n_bins = 8
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
) -> Result<Vec<gammalooprs::integrands::evaluation::SampleEvaluationResult>> {
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
        force_radius: false,
        minimal_output: false,
        momentum_space: false,
        points: array.view(),
        discrete_dims: None,
        graph_names: None,
        orientations: None,
    }
    .run(&mut cli.state)
}

fn evaluate_momentum_sample(
    cli: &mut CLIState,
    point: &[f64],
) -> Result<gammalooprs::integrands::evaluation::SampleEvaluationResult> {
    let array = Array2::from_shape_vec((1, point.len()), point.to_vec())?;
    let mut results = EvaluateSamples {
        process_id: None,
        integrand_name: None,
        use_arb_prec: false,
        force_radius: false,
        minimal_output: false,
        momentum_space: true,
        points: array.view(),
        discrete_dims: None,
        graph_names: None,
        orientations: None,
    }
    .run(&mut cli.state)?;
    Ok(results.remove(0))
}

fn evaluate_x_samples_minimal(
    cli: &mut CLIState,
    points: &[Vec<f64>],
) -> Result<Vec<gammalooprs::integrands::evaluation::SampleEvaluationResult>> {
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
        force_radius: false,
        minimal_output: true,
        momentum_space: false,
        points: array.view(),
        discrete_dims: None,
        graph_names: None,
        orientations: None,
    }
    .run(&mut cli.state)
}

fn evaluate_x_samples_precise(
    cli: &mut CLIState,
    points: &[Vec<f64>],
) -> Result<Vec<gammalooprs::integrands::evaluation::PreciseSampleEvaluationResult>> {
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
        force_radius: false,
        minimal_output: false,
        momentum_space: false,
        points: array.view(),
        discrete_dims: None,
        graph_names: None,
        orientations: None,
    }
    .run(&mut cli.state)
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

#[test]
#[serial]
fn lu_rust_generated_events_follow_graph_grouping_and_cut_ids() -> Result<()> {
    let _hack = enable_lu_e2e_hack();
    let mut cli = setup_sm_differential_lu_cli("lu_rust_generated_event_grouping")?;
    cli.run_command(
        "set process kv general.generate_events=true general.store_additional_weights_in_event=true",
    )?;

    let point = default_xspace_point(&cli)?;
    let mut results = evaluate_x_samples(&mut cli, std::slice::from_ref(&point))?;
    let result = results.remove(0);
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

    let second_group_cut_ids = event_groups[1]
        .iter()
        .map(|event| (event.cut_info.graph_id, event.cut_info.cut_id))
        .collect::<Vec<_>>();
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
    let _hack = enable_lu_e2e_hack();
    let mut cli = setup_sm_differential_lu_cli("lu_rust_precise_evaluate_samples")?;
    configure_leading_jet_quantity(&mut cli)?;
    add_leading_jet_observable(&mut cli)?;

    let point = default_xspace_point(&cli)?;
    let results = evaluate_x_samples_precise(&mut cli, &[point])?;
    assert_eq!(results.len(), 1);
    assert_eq!(results[0].observables.histograms.len(), 1);
    match &results[0].evaluation {
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
    let _hack = enable_lu_e2e_hack();
    let mut cli = setup_sm_differential_lu_cli("lu_rust_evaluate_samples_end_to_end")?;
    let point = default_xspace_point(&cli)?;

    let mut results = evaluate_x_samples(&mut cli, std::slice::from_ref(&point))?;
    let result = results.remove(0);
    assert!(result.evaluation.event_groups.is_empty());
    assert!(result.observables.histograms.is_empty());
    assert_eq!(metadata(&result).generated_event_count, 0);
    assert_eq!(metadata(&result).accepted_event_count, 0);
    assert!(result.evaluation.parameterization_jacobian.is_some());

    cli.run_command(
        "set process kv general.generate_events=true general.store_additional_weights_in_event=true",
    )?;
    let mut results = evaluate_x_samples(&mut cli, std::slice::from_ref(&point))?;
    let result = results.remove(0);
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
    configure_leading_jet_quantity(&mut cli)?;
    add_leading_jet_selector(&mut cli, 1_000_000.0)?;
    let mut results = evaluate_x_samples(&mut cli, std::slice::from_ref(&point))?;
    let result = results.remove(0);
    assert!(result.evaluation.event_groups.is_empty());
    assert!(result.observables.histograms.is_empty());
    assert!(
        metadata(&result).generated_event_count > 0,
        "selector-only mode should still generate temporary events"
    );
    assert_eq!(metadata(&result).accepted_event_count, 0);

    add_leading_jet_observable(&mut cli)?;
    update_leading_jet_selector(&mut cli, 0.0)?;
    let batch = vec![point.clone(), point];
    let results = evaluate_x_samples(&mut cli, &batch)?;
    assert_eq!(results.len(), 2);
    for result in &results {
        assert!(result.evaluation.event_groups.is_empty());
        assert!(
            result
                .observables
                .histograms
                .contains_key("leading_jet_pt_hist")
        );
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
        let histogram = result
            .observables
            .histograms
            .get("leading_jet_pt_hist")
            .expect("missing leading-jet histogram");
        assert_eq!(histogram.bins.len(), 8);
    }

    Ok(())
}

#[test]
#[serial]
fn lu_rust_momentum_space_evaluate_sample_reports_no_parameterization_jacobian() -> Result<()> {
    let _hack = enable_lu_e2e_hack();
    let mut cli = setup_sm_differential_lu_cli("lu_rust_momentum_space_evaluate_sample")?;
    cli.run_command("set process kv general.generate_events=true")?;
    let point = default_momentum_space_point(&cli)?;
    let result = evaluate_momentum_sample(&mut cli, &point)?;

    assert!(result.evaluation.parameterization_jacobian.is_none());
    assert!(metadata(&result).event_processing_time >= Duration::ZERO);
    assert!(result.evaluation.event_groups.len() <= metadata(&result).generated_event_count);
    let formatted = result.to_string();
    assert!(formatted.contains("parameterization jacobian"));
    assert!(formatted.contains("None"));

    Ok(())
}

#[test]
#[serial]
fn lu_rust_minimal_output_keeps_events_and_observables() -> Result<()> {
    let _hack = enable_lu_e2e_hack();
    let mut cli = setup_sm_differential_lu_cli("lu_rust_minimal_output")?;
    cli.run_command(
        "set process kv general.generate_events=true general.store_additional_weights_in_event=true",
    )?;
    configure_leading_jet_quantity(&mut cli)?;
    add_leading_jet_observable(&mut cli)?;

    let point = default_xspace_point(&cli)?;
    let mut results = evaluate_x_samples_minimal(&mut cli, std::slice::from_ref(&point))?;
    let result = results.remove(0);

    assert!(result.evaluation.evaluation_metadata.is_none());
    assert!(!result.evaluation.event_groups.is_empty());
    assert!(
        result
            .observables
            .histograms
            .contains_key("leading_jet_pt_hist")
    );

    let formatted = result.to_string();
    assert!(!formatted.contains("Evaluation metadata"));
    assert!(!formatted.contains("Stability results"));
    assert!(formatted.contains("Observable snapshots"));

    Ok(())
}
