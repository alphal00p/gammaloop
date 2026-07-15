use super::utils::*;
use super::*;
use gammaloop_integration_tests::{default_momentum_space_point, default_xspace_point};
use serde_json::json;

fn format_float_list(values: &[f64]) -> String {
    values
        .iter()
        .map(|value| format!("{value:.16e}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn format_axis(axis: &[f64]) -> String {
    axis.iter()
        .map(|value| format!("{value:.16e}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn first_axis(dimension: usize, scale: f64) -> Vec<f64> {
    let mut axis = vec![0.0; dimension];
    axis[0] = scale;
    axis
}

fn run_approach_command(
    cli: &mut gammaloop_integration_tests::CLIState,
    point: &[f64],
    axis: &[f64],
    n_cores: usize,
    output_path: &Path,
    extra_options: &str,
) -> Result<JsonValue> {
    cli.run_command(&format!(
        "approach -x {} --approach-axis {} --n-points 1 --n-cores {n_cores} --output-results {} {extra_options}",
        format_float_list(point),
        format_axis(axis),
        output_path.display(),
    ))?;
    Ok(serde_json::from_str(&std::fs::read_to_string(
        output_path,
    )?)?)
}

fn normalized_for_parallel_comparison(mut value: JsonValue) -> JsonValue {
    value["n_cores"] = JsonValue::from(0);
    for point in value["points"].as_array_mut().unwrap() {
        let Some(metadata) = point
            .get_mut("evaluation")
            .and_then(|evaluation| evaluation.get_mut("metadata"))
            .and_then(JsonValue::as_object_mut)
        else {
            continue;
        };
        for key in [
            "total_time_seconds",
            "parameterization_time_seconds",
            "integrand_evaluation_time_seconds",
            "evaluator_evaluation_time_seconds",
            "event_processing_time_seconds",
        ] {
            metadata.insert(key.to_string(), JsonValue::from(0.0));
        }
    }
    value
}

fn evaluated_points(value: &JsonValue) -> Vec<&JsonValue> {
    value["points"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|point| point["status"] == "evaluated")
        .collect()
}

fn ordered_point_indices(value: &JsonValue) -> Vec<usize> {
    value["points"]
        .as_array()
        .unwrap()
        .iter()
        .map(|point| point["index"].as_u64().unwrap() as usize)
        .collect()
}

fn multichannel_group(
    cli: &gammaloop_integration_tests::CLIState,
) -> Result<(usize, String, usize)> {
    let (process_id, integrand_name) = cli.state.find_integrand_ref(None, None)?;
    let integrand = cli
        .state
        .process_list
        .get_integrand(process_id, &integrand_name)?
        .require_generated()?;
    for graph_id in 0..integrand.graph_count() {
        let Some(graph_name) = integrand.graph_name_by_id(graph_id) else {
            continue;
        };
        let Ok(group_id) = integrand.resolve_group_id_by_master_name(graph_name) else {
            continue;
        };
        let channel_count = integrand.group_channel_count(group_id).unwrap_or(0);
        if channel_count > 1 {
            return Ok((group_id.0, graph_name.to_string(), channel_count));
        }
    }
    Err(color_eyre::eyre::eyre!(
        "expected at least one graph group with more than one LMB channel"
    ))
}

fn event_lmb_sample_ids(point: &JsonValue) -> std::collections::BTreeSet<usize> {
    point["evaluation"]["events"]
        .as_array()
        .unwrap()
        .iter()
        .filter_map(|event| event["lmb_sample_id"].as_u64().map(|id| id as usize))
        .collect()
}

fn contribution_lmb_sample_ids(point: &JsonValue) -> std::collections::BTreeSet<usize> {
    point["evaluation"]["contributions"]
        .as_array()
        .unwrap()
        .iter()
        .filter_map(|contribution| contribution["lmb_sample_id"].as_u64().map(|id| id as usize))
        .collect()
}

fn plot_fixture(process: &str, kind: &str, scale: f64) -> JsonValue {
    let mut points = Vec::new();
    let mut index = 0usize;
    for axis_index in 0..2 {
        for (axis_point_index, t) in [-1.0_f64, 0.0, 1.0].into_iter().enumerate() {
            let value = scale * (axis_index as f64 + 1.0) * (t + 1.5);
            let sign = if t < 0.0 { -1.0 } else { 1.0 };
            points.push(json!({
                "index": index,
                "axis_index": axis_index,
                "axis_point_index": axis_point_index,
                "t": t,
                "point": [0.5 + 0.01 * t, 0.25],
                "status": "evaluated",
                "evaluation": {
                    "integrand_result": {"re": value, "im": 0.1 * value},
                    "parameterization_jacobian": 2.0,
                    "integrator_weight": 1.0,
                    "total_weight": {"re": 2.0 * value, "im": 0.2 * value},
                    "event_weight_sum": {"re": sign * value, "im": 0.0},
                    "additional_weight_sums": {
                        "original": {"re": value, "im": 0.0},
                        "full_multiplicative_factor": {"re": 1.0, "im": 0.0},
                        "ct_8_3": {"re": 0.1 * sign * value, "im": 0.0}
                    },
                    "contributions": [{
                        "label": "event_weight graph=GL0 cut=0 edges=[1,2] orientation=0 lmb_sample=0",
                        "contribution": "event_weight",
                        "graph_id": axis_index,
                        "graph_name": format!("GL{axis_index}"),
                        "graph_group_id": axis_index,
                        "orientation_id": 0,
                        "cut_id": 0,
                        "cut_edges": if kind == "amplitude" { json!([]) } else { json!([1, 2]) },
                        "lmb_sample_id": axis_index,
                        "weight": {"re": sign * value, "im": 0.0}
                    }],
                    "events": [],
                    "metadata": {
                        "generated_event_count": 1,
                        "accepted_event_count": 1,
                        "is_nan": false,
                        "total_time_seconds": 0.0,
                        "parameterization_time_seconds": 0.0,
                        "integrand_evaluation_time_seconds": 0.0,
                        "evaluator_evaluation_time_seconds": 0.0,
                        "event_processing_time_seconds": 0.0
                    }
                }
            }));
            index += 1;
        }
    }

    json!({
        "schema_version": 1,
        "command": {
            "name": "approach",
            "use_arb_prec": false,
            "graph_id": null,
            "orientation_id": null,
            "discrete_dim": []
        },
        "process": {"id": 0, "name": process},
        "integrand": {"name": "default", "kind": kind},
        "space": "coordinate",
        "base_point": [0.5, 0.25],
        "axes": [[0.01, 0.0], [0.0, 0.02]],
        "spacing": {
            "kind": "linear",
            "n_points": 1,
            "min_abs_t": null,
            "t_values": [-1.0, 0.0, 1.0]
        },
        "n_cores": 1,
        "points_per_axis": 3,
        "evaluated_points": 6,
        "skipped_points": 0,
        "points": points
    })
}

fn command_succeeds(program: &str, args: &[&str]) -> bool {
    std::process::Command::new(program)
        .args(args)
        .status()
        .is_ok_and(|status| status.success())
}

#[test]
#[serial]
fn approach_cross_section_json_has_cut_metadata_and_parallel_matches_sequential() -> Result<()> {
    let test_name = "lu_approach_cross_section";
    let mut cli = setup_sm_differential_lu_cli(test_name)?;
    let point = default_xspace_point(&cli)?;
    let axis = first_axis(point.len(), 0.01);
    let output_dir = get_tests_workspace_path().join(test_name);
    let sequential_path = output_dir.join("approach_sequential.json");
    let parallel_path = output_dir.join("approach_parallel.json");

    let sequential = run_approach_command(&mut cli, &point, &axis, 1, &sequential_path, "")?;
    let parallel = run_approach_command(&mut cli, &point, &axis, 2, &parallel_path, "")?;

    assert_eq!(sequential["schema_version"], 1);
    assert_eq!(sequential["integrand"]["kind"], "cross_section");
    assert_eq!(sequential["space"], "coordinate");
    assert_eq!(sequential["evaluated_points"], 3);
    assert_eq!(sequential["skipped_points"], 0);

    let evaluated = evaluated_points(&sequential);
    assert_eq!(evaluated.len(), 3);
    let first_evaluation = &evaluated[0]["evaluation"];
    assert!(first_evaluation["parameterization_jacobian"].is_number());
    assert!(first_evaluation["total_weight"]["re"].is_number());
    assert!(first_evaluation["event_weight_sum"]["re"].is_number());
    assert!(
        first_evaluation["contributions"]
            .as_array()
            .unwrap()
            .iter()
            .any(|contribution| contribution["contribution"] == "event_weight")
    );
    assert!(
        first_evaluation["events"]
            .as_array()
            .unwrap()
            .iter()
            .any(|event| event["cut_edges"]
                .as_array()
                .is_some_and(|edges| !edges.is_empty())),
        "cross-section events should carry true cut-edge metadata",
    );

    assert_eq!(
        normalized_for_parallel_comparison(sequential),
        normalized_for_parallel_comparison(parallel),
        "parallel approach output should match sequential output modulo worker count and timing metadata",
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn approach_cross_section_momentum_space_uses_graph_orientation_selectors() -> Result<()> {
    let test_name = "lu_approach_cross_section_momentum";
    let mut cli = setup_sm_differential_lu_cli(test_name)?;
    let point = default_momentum_space_point(&cli)?;
    let axis = first_axis(point.len(), 0.001);
    let output_path = get_tests_workspace_path()
        .join(test_name)
        .join("approach_momentum.json");

    let result = run_approach_command(
        &mut cli,
        &point,
        &axis,
        1,
        &output_path,
        "--momentum-space --graph-id 0 --orientation-id 0",
    )?;

    assert_eq!(result["integrand"]["kind"], "cross_section");
    assert_eq!(result["space"], "momentum");
    assert_eq!(result["evaluated_points"], 3);
    assert_eq!(result["skipped_points"], 0);
    assert!(
        evaluated_points(&result)[0]["evaluation"]["events"]
            .as_array()
            .unwrap()
            .iter()
            .all(|event| event["orientation_id"] == 0),
        "momentum-space approach should honor explicit orientation selectors",
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn approach_amplitude_json_has_synthetic_cut_and_threshold_weights() -> Result<()> {
    let test_name = "gg_hhh_approach_amplitude";
    let mut cli = setup_gg_hhh_threshold_amplitude_cli(test_name)?;
    let point = vec![0.23, 0.41, 0.67];
    let axis = first_axis(point.len(), 0.01);
    let output_path = get_tests_workspace_path()
        .join(test_name)
        .join("approach_amplitude.json");

    let result = run_approach_command(
        &mut cli,
        &point,
        &axis,
        1,
        &output_path,
        "-p gg_hhh -i 1L --discrete-dim 0,0,0",
    )?;

    assert_eq!(result["integrand"]["kind"], "amplitude");
    let evaluated = evaluated_points(&result);
    assert_eq!(evaluated.len(), 3);
    let events = evaluated[0]["evaluation"]["events"].as_array().unwrap();
    assert!(!events.is_empty());
    assert!(events.iter().all(|event| event["cut_id"] == 0));
    assert!(
        events
            .iter()
            .all(|event| event["cut_edges"].as_array().unwrap().is_empty())
    );

    let additional_weight_sums = evaluated[0]["evaluation"]["additional_weight_sums"]
        .as_object()
        .unwrap();
    assert!(additional_weight_sums.contains_key("original"));
    assert!(additional_weight_sums.contains_key("full_multiplicative_factor"));
    assert!(
        additional_weight_sums
            .keys()
            .any(|label| label.starts_with("ct_") || label.starts_with("threshold_counterterm_")),
        "amplitude approach output should expose threshold-counterterm additional weights",
    );
    assert!(
        evaluated[0]["evaluation"]["contributions"]
            .as_array()
            .unwrap()
            .iter()
            .any(|contribution| contribution["label"]
                .as_str()
                .is_some_and(|label| label.contains("cut=0 edges=[]"))),
        "amplitude contributions should use the synthetic cut metadata",
    );

    let momentum_point = default_momentum_space_point(&cli)?;
    let momentum_axis = first_axis(momentum_point.len(), 0.001);
    let momentum_result = run_approach_command(
        &mut cli,
        &momentum_point,
        &momentum_axis,
        1,
        &get_tests_workspace_path()
            .join(test_name)
            .join("approach_amplitude_momentum.json"),
        "-p gg_hhh -i 1L --momentum-space --graph-id 0 --orientation-id 0",
    )?;
    assert_eq!(momentum_result["integrand"]["kind"], "amplitude");
    assert_eq!(momentum_result["space"], "momentum");
    assert_eq!(momentum_result["evaluated_points"], 3);
    assert_eq!(momentum_result["skipped_points"], 0);
    assert!(
        evaluated_points(&momentum_result)[0]["evaluation"]["events"]
            .as_array()
            .unwrap()
            .iter()
            .all(|event| event["cut_edges"].as_array().unwrap().is_empty()),
        "amplitude momentum-space approach should keep synthetic cut-edge metadata",
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn approach_sampling_modes_record_effective_lmb_sample_ids() -> Result<()> {
    let test_name = "lu_approach_lmb_sampling";
    let mut cli = setup_sm_differential_lu_cli(test_name)?;
    let point = default_xspace_point(&cli)?;
    let axis = first_axis(point.len(), 0.01);
    let output_dir = get_tests_workspace_path().join(test_name);

    let default_result = run_approach_command(
        &mut cli,
        &point,
        &axis,
        1,
        &output_dir.join("approach_default.json"),
        "",
    )?;
    assert_eq!(default_result["evaluated_points"], 3);
    assert_eq!(ordered_point_indices(&default_result), vec![0, 1, 2]);

    cli.run_command(
        r#"set process string '
[sampling]
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = true
lmb_channels = "summed"
'"#,
    )?;
    let (group_id, master_graph_name, channel_count) = multichannel_group(&cli)?;
    let summed_result = run_approach_command(
        &mut cli,
        &point,
        &axis,
        1,
        &output_dir.join("approach_lmb_summed.json"),
        &format!("--discrete-dim {group_id}"),
    )?;
    let summed_first_point = evaluated_points(&summed_result)[0];
    let summed_event_lmb_samples = event_lmb_sample_ids(summed_first_point);
    assert_eq!(
        summed_event_lmb_samples.len(),
        channel_count,
        "summed LMB multi-channeling should retain a distinct effective LMB sample id for every generated channel",
    );
    assert_eq!(
        contribution_lmb_sample_ids(summed_first_point),
        summed_event_lmb_samples,
        "per-contribution sums should stay distinct by lmb_sample_id in summed LMB mode",
    );

    cli.run_command(
        r#"set process string '
[sampling]
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = true
lmb_channels = "monte_carlo"
'"#,
    )?;
    let monte_carlo_result = run_approach_command(
        &mut cli,
        &point,
        &axis,
        1,
        &output_dir.join("approach_lmb_monte_carlo.json"),
        &format!("--discrete-dim {group_id},0"),
    )?;

    let event = evaluated_points(&monte_carlo_result)[0]["evaluation"]["events"]
        .as_array()
        .unwrap()
        .iter()
        .find(|event| event["lmb_channel_id"].is_number())
        .expect("discrete LMB sampling should record an event with lmb_channel_id");
    assert_eq!(event["lmb_channel_id"], 0);
    assert!(
        event["lmb_sample_id"].is_number(),
        "approach JSON should record the effective generated LMB basis id separately from the sampler channel id",
    );
    assert!(
        evaluated_points(&monte_carlo_result)[0]["evaluation"]["contributions"]
            .as_array()
            .unwrap()
            .iter()
            .any(|contribution| contribution["lmb_sample_id"].is_number()),
        "per-contribution sums should be keyed by lmb_sample_id",
    );

    cli.run_command(&format!(
        r#"set process string '
[sampling]
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = true
lmb_channels = "monte_carlo"
lmb_basis_ids = {{ {master_graph_name} = [1, 0] }}
'"#,
    ))?;
    let overridden_channel_zero = run_approach_command(
        &mut cli,
        &point,
        &axis,
        1,
        &output_dir.join("approach_lmb_override_channel_0.json"),
        &format!("--discrete-dim {group_id},0"),
    )?;
    let overridden_channel_one = run_approach_command(
        &mut cli,
        &point,
        &axis,
        1,
        &output_dir.join("approach_lmb_override_channel_1.json"),
        &format!("--discrete-dim {group_id},1"),
    )?;
    assert_eq!(
        event_lmb_sample_ids(evaluated_points(&overridden_channel_zero)[0]),
        std::collections::BTreeSet::from([1]),
        "lmb_sample_id should report the overridden generated basis id, not the channel slot",
    );
    assert_eq!(
        event_lmb_sample_ids(evaluated_points(&overridden_channel_one)[0]),
        std::collections::BTreeSet::from([0]),
        "lmb_sample_id should follow lmb_basis_ids overrides for each channel",
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn plot_approach_result_script_smoke() -> Result<()> {
    if !command_succeeds("python3", &["-c", "import matplotlib"]) {
        eprintln!("skipping plot smoke test because python3/matplotlib is unavailable");
        return Ok(());
    }
    if !command_succeeds("pdfinfo", &["-v"]) || !command_succeeds("pdftoppm", &["-v"]) {
        eprintln!("skipping plot smoke test because Poppler tools are unavailable");
        return Ok(());
    }

    let test_name = "plot_approach_result_script_smoke";
    let output_dir = get_tests_workspace_path().join(test_name);
    clean_test(&output_dir);
    std::fs::create_dir_all(&output_dir)?;

    let xs_path = output_dir.join("approach_xs.json");
    let amp_path = output_dir.join("approach_amp.json");
    std::fs::write(
        &xs_path,
        serde_json::to_vec_pretty(&plot_fixture("smoke_xs", "cross_section", 1.0))?,
    )?;
    std::fs::write(
        &amp_path,
        serde_json::to_vec_pretty(&plot_fixture("smoke_amp", "amplitude", 0.7))?,
    )?;

    let workspace_root = gammaloop_integration_tests::workspace_root();
    let script = workspace_root.join("assets/plot_approach_result.py");
    let default_pdf_path = output_dir.join("approach_default_pages.pdf");
    let default_status = std::process::Command::new("python3")
        .arg(&script)
        .arg(&xs_path)
        .arg(&amp_path)
        .arg("--output")
        .arg(&default_pdf_path)
        .arg("--include-contribution")
        .arg("event_weight|total_weight|ct_8_3")
        .arg("--exclude-contribution")
        .arg("full_multiplicative_factor")
        .status()?;
    assert!(
        default_status.success(),
        "default plot_approach_result.py invocation should succeed"
    );

    let default_pdfinfo_output = std::process::Command::new("pdfinfo")
        .arg(&default_pdf_path)
        .output()?;
    assert!(default_pdfinfo_output.status.success());
    let default_pdfinfo_text = String::from_utf8_lossy(&default_pdfinfo_output.stdout);
    let default_page_count = default_pdfinfo_text
        .lines()
        .find_map(|line| line.strip_prefix("Pages:"))
        .and_then(|value| value.trim().parse::<usize>().ok());
    assert_eq!(
        default_page_count,
        Some(4),
        "default plotting should produce one PDF page per result/axis: {default_pdfinfo_text}",
    );

    let pdf_path = output_dir.join("approach_combined.pdf");
    let status = std::process::Command::new("python3")
        .arg(script)
        .arg(&xs_path)
        .arg(&amp_path)
        .arg("--output")
        .arg(&pdf_path)
        .arg("--combine-plots")
        .arg("--combine-axes")
        .arg("--sum-lmb-samples-per-cut")
        .arg("--y-log-scale")
        .arg("--component")
        .arg("real")
        .arg("--include-contribution")
        .arg("event_weight|total_weight|ct_8_3")
        .arg("--exclude-contribution")
        .arg("full_multiplicative_factor")
        .status()?;
    assert!(status.success(), "plot_approach_result.py should succeed");

    let pdfinfo_output = std::process::Command::new("pdfinfo")
        .arg(&pdf_path)
        .output()?;
    assert!(pdfinfo_output.status.success());
    let pdfinfo_text = String::from_utf8_lossy(&pdfinfo_output.stdout);
    assert!(pdfinfo_text.contains("Pages:"));

    let rendered_prefix = output_dir.join("approach_combined_page");
    let render_status = std::process::Command::new("pdftoppm")
        .arg("-png")
        .arg("-singlefile")
        .arg(&pdf_path)
        .arg(&rendered_prefix)
        .status()?;
    assert!(render_status.success());
    assert!(rendered_prefix.with_extension("png").metadata()?.len() > 0);

    clean_test(&output_dir);
    Ok(())
}
