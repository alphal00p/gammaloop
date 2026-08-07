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
    let cut_group_id = (kind != "amplitude").then_some(0usize);
    let cut_edges = if kind == "amplitude" {
        Vec::new()
    } else {
        vec![1usize, 2]
    };
    let left_side = if kind == "amplitude" {
        "amplitude"
    } else {
        "left"
    };
    let right_side = if kind == "amplitude" {
        "amplitude"
    } else {
        "right"
    };
    let mut points = Vec::new();
    let mut index = 0usize;
    for axis_index in 0..2 {
        for (axis_point_index, t) in [-1.0_f64, -0.1, -0.01, -0.001, 0.0, 0.001, 0.01, 0.1, 1.0]
            .into_iter()
            .enumerate()
        {
            let value = scale * (axis_index as f64 + 1.0) * (t + 1.5);
            let sign = if t < 0.0 { -1.0 } else { 1.0 };
            let mut component_sums = vec![(0.0_f64, 0.0_f64, 0usize, 0usize); 14];
            let mut original_sum = (0.0_f64, 0.0_f64);
            let mut event_sum = (0.0_f64, 0.0_f64);
            let events = (0..2)
                .map(|event_index| {
                    let original = (
                        0.6 * sign * value,
                        0.15 * sign * value * (event_index as f64 + 1.0),
                    );
                    original_sum.0 += original.0;
                    original_sum.1 += original.1;
                    let components = (0..14)
                        .map(|component_id| {
                            let first_multiplier =
                                if event_index == 0 { 0.25 } else { 0.75 } + 0.01 * t.abs();
                            let second_multiplier = 0.4 + 0.1 * event_index as f64;
                            let multiplier_values = if component_id < 6 {
                                vec![first_multiplier]
                            } else {
                                vec![first_multiplier, second_multiplier]
                            };
                            let effective_multiplier = multiplier_values.iter().product::<f64>();
                            let evaluation_skipped =
                                component_id == 1 && event_index == 0 && t < 0.0;
                            let bare_re = -0.001
                                * (component_id as f64 + 1.0)
                                * (event_index as f64 + 1.0)
                                * sign
                                * value;
                            let bare_im = 0.25 * bare_re;
                            let (weighted_re, weighted_im) = if evaluation_skipped {
                                (0.0, 0.0)
                            } else {
                                (
                                    bare_re * effective_multiplier,
                                    bare_im * effective_multiplier,
                                )
                            };
                            component_sums[component_id].0 += weighted_re;
                            component_sums[component_id].1 += weighted_im;
                            component_sums[component_id].2 += 1;
                            component_sums[component_id].3 += usize::from(evaluation_skipped);
                            let pair = component_id >= 6;
                            json!({
                                "component_id": component_id,
                                "occurrence": {
                                    "kind": "local_unitarity",
                                    "overlap_groups": if pair { json!([0, 1]) } else { json!([0]) },
                                    "left_threshold_order": 1,
                                    "right_threshold_order": if pair { json!(2) } else { JsonValue::Null },
                                    "lu_cut_order": event_index
                                },
                                "multiplier_values": multiplier_values,
                                "effective_multiplier": effective_multiplier,
                                "bare": if evaluation_skipped {
                                    JsonValue::Null
                                } else {
                                    json!({"re": bare_re, "im": bare_im})
                                },
                                "weighted": {"re": weighted_re, "im": weighted_im},
                                "evaluation_skipped": evaluation_skipped
                            })
                        })
                        .collect::<Vec<_>>();
                    let counterterm_re = components
                        .iter()
                        .map(|component| component["weighted"]["re"].as_f64().unwrap())
                        .sum::<f64>();
                    let counterterm_im = components
                        .iter()
                        .map(|component| component["weighted"]["im"].as_f64().unwrap())
                        .sum::<f64>();
                    let event_weight = (original.0 + counterterm_re, original.1 + counterterm_im);
                    event_sum.0 += event_weight.0;
                    event_sum.1 += event_weight.1;
                    json!({
                        "event_group_index": 0,
                        "event_index": event_index,
                        "graph_id": 0,
                        "graph_name": "GL0",
                        "graph_group_id": 0,
                        "orientation_id": 0,
                        "cut_id": event_index,
                        "cut_edges": if kind == "amplitude" { json!([]) } else { json!([1, 2]) },
                        "lmb_channel_id": event_index,
                        "lmb_sample_id": event_index,
                        "weight": {"re": event_weight.0, "im": event_weight.1},
                        "full_multiplicative_factor": {"re": 2.0, "im": 0.0},
                        "display_only_weights": {
                            "original": {"re": 0.5 * original.0, "im": 0.5 * original.1}
                        },
                        "threshold_counterterms": {
                            "original": {"re": original.0, "im": original.1},
                            "components": components
                        }
                    })
                })
                .collect::<Vec<_>>();
            let counterterm_sum = component_sums.iter().fold((0.0, 0.0), |sum, component| {
                (sum.0 + component.0, sum.1 + component.1)
            });
            let summary_components = component_sums
                .iter()
                .enumerate()
                .map(
                    |(component_id, (re, im, occurrence_count, skipped_count))| {
                        json!({
                            "graph_id": 0,
                            "component_id": component_id,
                            "weighted": {"re": re, "im": im},
                            "occurrence_count": occurrence_count,
                            "skipped_count": skipped_count
                        })
                    },
                )
                .collect::<Vec<_>>();
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
                    "event_weight_sum": {"re": event_sum.0, "im": event_sum.1},
                    "additional_contribution_sums": {
                        "original": {"re": original_sum.0, "im": original_sum.1},
                        "ct_8_3": {"re": counterterm_sum.0, "im": counterterm_sum.1}
                    },
                    "contributions": [{
                        "label": "event_weight graph=GL0 cut=0 edges=[1,2] orientation=0 lmb_sample=0",
                        "contribution": "event_weight",
                        "graph_id": 0,
                        "graph_name": "GL0",
                        "graph_group_id": 0,
                        "orientation_id": 0,
                        "cut_id": 0,
                        "cut_edges": if kind == "amplitude" { json!([]) } else { json!([1, 2]) },
                        "lmb_sample_id": 0,
                        "weight": {"re": sign * value, "im": 0.0}
                    }],
                    "events": events,
                    "threshold_counterterm_summary": {
                        "original": {"re": original_sum.0, "im": original_sum.1},
                        "counterterm_sum": {"re": counterterm_sum.0, "im": counterterm_sum.1},
                        "reconstructed_total": {"re": event_sum.0, "im": event_sum.1},
                        "decomposed_event_sum": {"re": event_sum.0, "im": event_sum.1},
                        "closure": {"re": 0.0, "im": 0.0},
                        "components": summary_components
                    },
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
        "schema_version": 3,
        "command": {
            "name": "approach",
            "use_arb_prec": false,
            "skip_midpoint": false,
            "graph_id": null,
            "graph_name": null,
            "orientation_id": null,
            "discrete_dim": []
        },
        "process": {"id": 0, "name": process},
        "integrand": {
            "name": "default",
            "kind": kind,
            "threshold_counterterms": [{
                "graph_id": 0,
                "registry": {
                    "graph_name": "GL0",
                    "variants": [{
                        "variant_id": 0,
                        "name": "forced_1l",
                        "cut_group_id": cut_group_id,
                        "associations": [{
                            "cut_edges": cut_edges,
                            "threshold_edges": [3, 4]
                        }],
                        "side": left_side,
                        "resolved_subspace": [3]
                    }, {
                        "variant_id": 1,
                        "name": "default",
                        "cut_group_id": cut_group_id,
                        "associations": [{
                            "cut_edges": cut_edges,
                            "threshold_edges": [5, 6]
                        }],
                        "side": right_side,
                        "resolved_subspace": [5]
                    }, {
                        "variant_id": 2,
                        "name": "alternate",
                        "cut_group_id": cut_group_id,
                        "associations": [{
                            "cut_edges": cut_edges,
                            "threshold_edges": [7, 8]
                        }],
                        "side": right_side,
                        "resolved_subspace": [7]
                    }],
                    "components": [{
                        "component_id": 0, "cut_group_id": cut_group_id, "kind": "local", "variant_ids": [0]
                    }, {
                        "component_id": 1, "cut_group_id": cut_group_id, "kind": "integrated", "variant_ids": [0]
                    }, {
                        "component_id": 2, "cut_group_id": cut_group_id, "kind": "local", "variant_ids": [1]
                    }, {
                        "component_id": 3, "cut_group_id": cut_group_id, "kind": "integrated", "variant_ids": [1]
                    }, {
                        "component_id": 4, "cut_group_id": cut_group_id, "kind": "local", "variant_ids": [2]
                    }, {
                        "component_id": 5, "cut_group_id": cut_group_id, "kind": "integrated", "variant_ids": [2]
                    }, {
                        "component_id": 6, "cut_group_id": cut_group_id, "kind": "local_local", "variant_ids": [0, 1]
                    }, {
                        "component_id": 7, "cut_group_id": cut_group_id, "kind": "local_integrated", "variant_ids": [0, 1]
                    }, {
                        "component_id": 8, "cut_group_id": cut_group_id, "kind": "integrated_local", "variant_ids": [0, 1]
                    }, {
                        "component_id": 9, "cut_group_id": cut_group_id, "kind": "integrated_integrated", "variant_ids": [0, 1]
                    }, {
                        "component_id": 10, "cut_group_id": cut_group_id, "kind": "local_local", "variant_ids": [0, 2]
                    }, {
                        "component_id": 11, "cut_group_id": cut_group_id, "kind": "local_integrated", "variant_ids": [0, 2]
                    }, {
                        "component_id": 12, "cut_group_id": cut_group_id, "kind": "integrated_local", "variant_ids": [0, 2]
                    }, {
                        "component_id": 13, "cut_group_id": cut_group_id, "kind": "integrated_integrated", "variant_ids": [0, 2]
                    }]
                }
            }]
        },
        "space": "coordinate",
        "parameter_name": "lambda",
        "base_point": [0.5, 0.25],
        "axes": [{
            "label": "soft-threshold direction",
            "vector": [0.01, 0.0]
        }, {
            "label": "second direction",
            "vector": [0.0, 0.02]
        }],
        "spacing": {
            "kind": "logarithmic",
            "n_points": 4,
            "min_abs_t": 0.001,
            "max_abs_t": 1.0,
            "t_values": [-1.0, -0.1, -0.01, -0.001, 0.0, 0.001, 0.01, 0.1, 1.0]
        },
        "n_cores": 1,
        "points_per_axis": 9,
        "evaluated_points": 18,
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

    assert_eq!(sequential["schema_version"], 3);
    assert_eq!(sequential["parameter_name"], "t");
    assert_eq!(sequential["axes"][0]["label"], "axis 1");
    assert_eq!(
        sequential["axes"][0]["vector"],
        serde_json::to_value(&axis)?
    );
    assert_eq!(sequential["spacing"]["max_abs_t"], 1.0);
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
        "--momentum-space --graph-id 0 --orientation-id 0 --skip-midpoint \
         --parameter-name lambda --axis-label correlated-soft-threshold --max-abs-t 0.01",
    )?;

    assert_eq!(result["integrand"]["kind"], "cross_section");
    assert_eq!(result["space"], "momentum");
    assert_eq!(result["command"]["skip_midpoint"], true);
    assert_eq!(result["command"]["graph_name"], "GL0");
    assert_eq!(result["parameter_name"], "lambda");
    assert_eq!(result["axes"][0]["label"], "correlated-soft-threshold");
    assert_eq!(result["spacing"]["max_abs_t"], 0.01);
    assert_eq!(result["points_per_axis"], 2);
    assert_eq!(result["spacing"]["t_values"], json!([-0.01, 0.01]));
    assert_eq!(result["evaluated_points"], 2);
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

    let additional_contribution_sums = evaluated[0]["evaluation"]["additional_contribution_sums"]
        .as_object()
        .unwrap();
    assert!(additional_contribution_sums.contains_key("original"));
    assert!(!additional_contribution_sums.contains_key("full_multiplicative_factor"));
    assert!(
        additional_contribution_sums
            .keys()
            .any(|label| label.starts_with("ct_") || label.starts_with("threshold_counterterm_")),
        "amplitude approach output should expose threshold-counterterm additional weights",
    );
    let additive_re = additional_contribution_sums
        .values()
        .map(|weight| weight["re"].as_f64().unwrap())
        .sum::<f64>();
    let additive_im = additional_contribution_sums
        .values()
        .map(|weight| weight["im"].as_f64().unwrap())
        .sum::<f64>();
    let event_sum = &evaluated[0]["evaluation"]["event_weight_sum"];
    let scale = event_sum["re"]
        .as_f64()
        .unwrap()
        .abs()
        .max(event_sum["im"].as_f64().unwrap().abs())
        .max(1.0);
    assert!((additive_re - event_sum["re"].as_f64().unwrap()).abs() < 1.0e-12 * scale);
    assert!((additive_im - event_sum["im"].as_f64().unwrap()).abs() < 1.0e-12 * scale);
    assert!(events.iter().all(|event| {
        event["full_multiplicative_factor"].is_object()
            && event["display_only_weights"].is_object()
            && event.get("additional_weights").is_none()
    }));
    assert!(
        evaluated[0]["evaluation"]["contributions"]
            .as_array()
            .unwrap()
            .iter()
            .all(|contribution| contribution["contribution"] != "full_multiplicative_factor"),
        "multiplicative factors must not be presented as additive contributions",
    );
    if let Some(summary) = evaluated[0]["evaluation"]["threshold_counterterm_summary"].as_object() {
        assert!(summary["decomposed_event_sum"].is_object());
        assert!(summary["closure"].is_object());
        assert!(
            summary["components"]
                .as_array()
                .unwrap()
                .iter()
                .all(|component| {
                    component["occurrence_count"].is_u64() && component["skipped_count"].is_u64()
                })
        );
    }
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
    if !command_succeeds("pdfinfo", &["-v"])
        || !command_succeeds("pdftoppm", &["-v"])
        || !command_succeeds("pdftotext", &["-v"])
    {
        eprintln!("skipping plot smoke test because Poppler tools are unavailable");
        return Ok(());
    }

    let test_name = "plot_approach_result_script_smoke";
    let output_dir = get_tests_workspace_path().join(test_name);
    clean_test(&output_dir);
    std::fs::create_dir_all(&output_dir)?;

    let xs_path = output_dir.join("approach_xs.json");
    let amp_path = output_dir.join("approach_amp.json");
    let amp_threshold_path = output_dir.join("approach_amp_threshold.json");
    let amp_detailed_only_path = output_dir.join("approach_amp_detailed_only.json");
    let supplemental_path = output_dir.join("kinematic_series.json");
    let outdated_path = output_dir.join("approach_schema_v2.json");
    std::fs::write(
        &xs_path,
        serde_json::to_vec_pretty(&plot_fixture("smoke_xs", "cross_section", 1.0))?,
    )?;
    let mut outdated = plot_fixture("outdated", "cross_section", 1.0);
    outdated["schema_version"] = json!(2);
    std::fs::write(&outdated_path, serde_json::to_vec_pretty(&outdated)?)?;
    let amplitude = plot_fixture("smoke_amp", "amplitude", 0.7);
    std::fs::write(&amp_threshold_path, serde_json::to_vec_pretty(&amplitude)?)?;
    let mut amplitude_with_detailed_events_only = amplitude.clone();
    for point in amplitude_with_detailed_events_only["points"]
        .as_array_mut()
        .unwrap()
    {
        point["evaluation"]
            .as_object_mut()
            .unwrap()
            .remove("threshold_counterterm_summary");
    }
    std::fs::write(
        &amp_detailed_only_path,
        serde_json::to_vec_pretty(&amplitude_with_detailed_events_only)?,
    )?;
    let mut amplitude_without_threshold_details = amplitude;
    amplitude_without_threshold_details["integrand"]
        .as_object_mut()
        .unwrap()
        .remove("threshold_counterterms");
    for point in amplitude_without_threshold_details["points"]
        .as_array_mut()
        .unwrap()
    {
        point["evaluation"]
            .as_object_mut()
            .unwrap()
            .remove("threshold_counterterm_summary");
        for event in point["evaluation"]["events"].as_array_mut().unwrap() {
            event
                .as_object_mut()
                .unwrap()
                .remove("threshold_counterterms");
        }
    }
    std::fs::write(
        &amp_path,
        serde_json::to_vec_pretty(&amplitude_without_threshold_details)?,
    )?;
    let t_values = [-1.0_f64, -0.1, -0.01, -0.001, 0.0, 0.001, 0.01, 0.1, 1.0];
    std::fs::write(
        &supplemental_path,
        serde_json::to_vec_pretty(&json!({
            "schema_version": 1,
            "series": [{
                "label": "soft momentum",
                "axis_index": 0,
                "values": t_values.iter().enumerate().map(|(point_index, t)| {
                    let magnitude = t.abs();
                    json!({
                        "point_index": point_index,
                        "value": if magnitude > 0.1 {
                            json!(10.0 + magnitude)
                        } else {
                            json!(magnitude * magnitude)
                        }
                    })
                }).collect::<Vec<_>>()
            }, {
                "label": "diagnostic with gap",
                "axis_index": 0,
                "values": t_values.iter().enumerate().map(|(point_index, t)| json!({
                    "point_index": point_index,
                    "value": if point_index == 2 { JsonValue::Null } else { json!(t.abs()) }
                })).collect::<Vec<_>>()
            }]
        }))?,
    )?;

    let workspace_root = gammaloop_integration_tests::workspace_root();
    let script = workspace_root.join("assets/plot_approach_result.py");
    let outdated_output = std::process::Command::new("python3")
        .arg(&script)
        .arg(&outdated_path)
        .arg("--output")
        .arg(output_dir.join("outdated.pdf"))
        .output()?;
    assert!(!outdated_output.status.success());
    assert!(
        String::from_utf8_lossy(&outdated_output.stderr).contains("requires schema_version 3"),
        "the plotter must reject old approach schemas instead of guessing: {}",
        String::from_utf8_lossy(&outdated_output.stderr),
    );

    let default_pdf_path = output_dir.join("approach_default_pages.pdf");
    let default_status = std::process::Command::new("python3")
        .arg(&script)
        .arg(&xs_path)
        .arg(&amp_path)
        .arg("--output")
        .arg(&default_pdf_path)
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
    let per_page_pdfinfo = std::process::Command::new("pdfinfo")
        .args(["-f", "1", "-l", "4"])
        .arg(&default_pdf_path)
        .output()?;
    assert!(per_page_pdfinfo.status.success());
    let pdfinfo_text = String::from_utf8_lossy(&per_page_pdfinfo.stdout);
    let page_sizes = pdfinfo_text
        .lines()
        .filter(|line| line.starts_with("Page") && line.contains(" size:"))
        .filter_map(|line| line.split_once("size:").map(|(_, size)| size.trim()))
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(
        page_sizes.len(),
        1,
        "all plotter pages must retain one fixed physical size: {page_sizes:?}",
    );
    assert!(
        page_sizes
            .first()
            .is_some_and(|size| size.contains("972 x 590.4 pts")),
        "the fixed 13.5 by 8.2 inch page size should be preserved: {page_sizes:?}",
    );

    let pdf_path = output_dir.join("approach_combined.pdf");
    let status = std::process::Command::new("python3")
        .arg(&script)
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
        .arg("--t-branch")
        .arg("both")
        .arg("--x-log-scale")
        .arg("--branch-layout")
        .arg("split")
        .arg("--result-label")
        .arg("threshold CT on")
        .arg("--result-label")
        .arg("threshold CT off")
        .arg("--include-contribution")
        .arg("event_weight|total_weight|ct_8_3")
        .arg("--exclude-contribution")
        .arg("full_multiplicative_factor")
        .status()?;
    assert!(status.success(), "plot_approach_result.py should succeed");
    let combined_text = std::process::Command::new("pdftotext")
        .arg(&pdf_path)
        .arg("-")
        .output()?;
    assert!(combined_text.status.success());
    let combined_text = String::from_utf8_lossy(&combined_text.stdout);
    assert!(combined_text.contains("axis 1: soft-threshold direction"));
    assert!(combined_text.contains("axis 2: second direction"));

    let threshold_pdf_path = output_dir.join("approach_threshold_decomposition.pdf");
    let threshold_status = std::process::Command::new("python3")
        .arg(&script)
        .arg(&xs_path)
        .arg("--output")
        .arg(&threshold_pdf_path)
        .arg("--threshold-decomposition")
        .arg("all")
        .arg("--include-contribution")
        .arg("^threshold:")
        .arg("--t-branch")
        .arg("positive")
        .arg("--x-log-scale")
        .arg("--fit-log-slope")
        .arg("--hide-info-box")
        .arg("--axis-id")
        .arg("1")
        .arg("--title")
        .arg("Threshold decomposition smoke test")
        .status()?;
    assert!(
        threshold_status.success(),
        "threshold decomposition plot should succeed"
    );
    assert!(threshold_pdf_path.metadata()?.len() > 0);

    let supplemental_pdf_path = output_dir.join("approach_supplemental.pdf");
    let supplemental_output = std::process::Command::new("python3")
        .arg(&script)
        .arg(&xs_path)
        .arg("--output")
        .arg(&supplemental_pdf_path)
        .arg("--series-json")
        .arg(&supplemental_path)
        .arg("--include-contribution")
        .arg("^series:")
        .arg("--axis-id")
        .arg("1")
        .arg("--t-branch")
        .arg("both")
        .arg("--x-log-scale")
        .arg("--branch-layout")
        .arg("split")
        .arg("--fit-log-slope")
        .arg("--fit-points")
        .arg("3")
        .arg("--y-label")
        .arg("kinematic distance [GeV]")
        .output()?;
    assert!(
        supplemental_output.status.success(),
        "supplemental-series plot should succeed: {}",
        String::from_utf8_lossy(&supplemental_output.stderr),
    );
    assert!(
        String::from_utf8_lossy(&supplemental_output.stderr)
            .contains("Rendered non-finite/null value(s) as gaps"),
        "null supplemental values should be reported as plot gaps, never zeros",
    );
    let supplemental_text = std::process::Command::new("pdftotext")
        .arg(&supplemental_pdf_path)
        .arg("-")
        .output()?;
    assert!(supplemental_text.status.success());
    let supplemental_text = String::from_utf8_lossy(&supplemental_text.stdout);
    assert!(supplemental_text.contains("soft-threshold direction"));
    assert!(supplemental_text.contains("kinematic distance [GeV]"));
    assert!(supplemental_text.contains("lambda < 0"));
    assert!(supplemental_text.contains("p=+2.00"));
    assert!(supplemental_text.contains("R^2=1.000"));

    let threshold_report_path = output_dir.join("approach_threshold_report.pdf");
    let threshold_report_status = std::process::Command::new("python3")
        .arg(&script)
        .arg(&xs_path)
        .arg("--output")
        .arg(&threshold_report_path)
        .arg("--threshold-report")
        .arg("--axis-id")
        .arg("1")
        .arg("--t-branch")
        .arg("both")
        .arg("--x-log-scale")
        .arg("--branch-layout")
        .arg("split")
        .arg("--result-label")
        .arg("synthetic threshold report")
        .arg("--hide-info-box")
        .status()?;
    assert!(
        threshold_report_status.success(),
        "threshold report plot should succeed"
    );
    let threshold_report_pdfinfo = std::process::Command::new("pdfinfo")
        .arg(&threshold_report_path)
        .output()?;
    assert!(threshold_report_pdfinfo.status.success());
    let threshold_report_pdfinfo = String::from_utf8_lossy(&threshold_report_pdfinfo.stdout);
    let threshold_report_page_count = threshold_report_pdfinfo
        .lines()
        .find_map(|line| line.strip_prefix("Pages:"))
        .and_then(|value| value.trim().parse::<usize>().ok());
    assert_eq!(
        threshold_report_page_count,
        Some(3),
        "threshold report should contain closure, single-variant, and pair pages: {threshold_report_pdfinfo}",
    );
    let threshold_report_text = std::process::Command::new("pdftotext")
        .arg(&threshold_report_path)
        .arg("-")
        .output()?;
    assert!(threshold_report_text.status.success());
    let threshold_report_text = String::from_utf8_lossy(&threshold_report_text.stdout);
    assert!(threshold_report_text.contains("decomposed event sum"));
    assert!(threshold_report_text.contains("original O"));
    assert!(threshold_report_text.contains("weighted CT"));
    assert!(threshold_report_text.contains("reconstructed O+"));

    let filtered_report_path = output_dir.join("approach_filtered_bare_report.pdf");
    let filtered_report_output = std::process::Command::new("python3")
        .arg(&script)
        .arg(&xs_path)
        .arg("--output")
        .arg(&filtered_report_path)
        .arg("--threshold-report")
        .arg("--threshold-report-section")
        .arg("singles")
        .arg("--threshold-quantity")
        .arg("weighted")
        .arg("--threshold-quantity")
        .arg("bare")
        .arg("--include-contribution")
        .arg("forced_1l")
        .arg("--axis-id")
        .arg("1")
        .arg("--t-branch")
        .arg("both")
        .arg("--x-log-scale")
        .arg("--hide-info-box")
        .output()?;
    assert!(
        filtered_report_output.status.success(),
        "filtered weighted/bare report should succeed: {}",
        String::from_utf8_lossy(&filtered_report_output.stderr),
    );
    assert!(
        String::from_utf8_lossy(&filtered_report_output.stderr)
            .contains("missing because the component evaluation was skipped"),
        "missing bare values should produce explicit gaps and diagnostics",
    );
    let filtered_report_text = std::process::Command::new("pdftotext")
        .arg(&filtered_report_path)
        .arg("-")
        .output()?;
    assert!(filtered_report_text.status.success());
    let filtered_report_text = String::from_utf8_lossy(&filtered_report_text.stdout);
    assert!(filtered_report_text.contains("forced_1l"));
    assert!(filtered_report_text.contains("unmultiplied CT"));
    assert!(!filtered_report_text.contains("alternate"));

    let multiplier_report_path = output_dir.join("approach_multiplier_report.pdf");
    let multiplier_report_output = std::process::Command::new("python3")
        .arg(&script)
        .arg(&xs_path)
        .arg("--output")
        .arg(&multiplier_report_path)
        .arg("--threshold-report")
        .arg("--threshold-report-section")
        .arg("multipliers")
        .arg("--include-contribution")
        .arg("component c0 ")
        .arg("--axis-id")
        .arg("1")
        .arg("--t-branch")
        .arg("both")
        .arg("--x-log-scale")
        .arg("--branch-layout")
        .arg("split")
        .args(["--y-range", "-1", "2"])
        .arg("--hide-info-box")
        .output()?;
    assert!(
        multiplier_report_output.status.success(),
        "multiplier report should succeed: {}",
        String::from_utf8_lossy(&multiplier_report_output.stderr),
    );
    let multiplier_report_text = std::process::Command::new("pdftotext")
        .arg(&multiplier_report_path)
        .arg("-")
        .output()?;
    assert!(multiplier_report_text.status.success());
    let multiplier_report_text = String::from_utf8_lossy(&multiplier_report_text.stdout);
    assert!(multiplier_report_text.contains("eg0/e0"));
    assert!(multiplier_report_text.contains("eg0/e1"));
    assert!(multiplier_report_text.contains("f_eff"));
    assert!(multiplier_report_text.contains("multiplier value"));
    assert!(multiplier_report_text.contains("GL0"));
    assert!(multiplier_report_text.contains("cg0"));
    assert!(multiplier_report_text.contains("c0 L"));
    assert!(multiplier_report_text.contains("forced_1l"));

    let invalid_multiplier_range = std::process::Command::new("python3")
        .arg(&script)
        .arg(&xs_path)
        .arg("--output")
        .arg(output_dir.join("invalid_multiplier_range.pdf"))
        .arg("--threshold-report")
        .args(["--threshold-report-section", "multipliers"])
        .arg("--multiplier-log-scale")
        .args(["--y-range", "-1", "2"])
        .output()?;
    assert!(!invalid_multiplier_range.status.success());
    assert!(
        String::from_utf8_lossy(&invalid_multiplier_range.stderr)
            .contains("lower bound must be positive with logarithmic scale")
    );

    let amplitude_threshold_report_path = output_dir.join("approach_amp_threshold_report.pdf");
    let amplitude_threshold_report_status = std::process::Command::new("python3")
        .arg(&script)
        .arg(&amp_threshold_path)
        .arg("--output")
        .arg(&amplitude_threshold_report_path)
        .arg("--threshold-report")
        .arg("--component")
        .arg("imag")
        .arg("--axis-id")
        .arg("1")
        .status()?;
    assert!(
        amplitude_threshold_report_status.success(),
        "threshold reports should accept amplitude registries with null cut-group IDs",
    );
    let amplitude_threshold_report_text = std::process::Command::new("pdftotext")
        .arg(&amplitude_threshold_report_path)
        .arg("-")
        .output()?;
    assert!(amplitude_threshold_report_text.status.success());
    let amplitude_threshold_report_text =
        String::from_utf8_lossy(&amplitude_threshold_report_text.stdout);
    assert!(
        amplitude_threshold_report_text.contains("imag weight"),
        "--component imag should apply to every threshold-report section: {amplitude_threshold_report_text}",
    );

    let amplitude_threshold_decomposition_path =
        output_dir.join("approach_amp_threshold_decomposition.pdf");
    let amplitude_threshold_decomposition_status = std::process::Command::new("python3")
        .arg(&script)
        .arg(&amp_threshold_path)
        .arg("--output")
        .arg(&amplitude_threshold_decomposition_path)
        .arg("--threshold-decomposition")
        .arg("all")
        .arg("--include-contribution")
        .arg("^threshold:")
        .status()?;
    assert!(
        amplitude_threshold_decomposition_status.success(),
        "threshold decompositions should accept amplitude registries with null cut-group IDs",
    );

    let empty_threshold_report_path = output_dir.join("approach_empty_threshold_report.pdf");
    let empty_threshold_report_output = std::process::Command::new("python3")
        .arg(&script)
        .arg(&amp_path)
        .arg("--output")
        .arg(&empty_threshold_report_path)
        .arg("--threshold-report")
        .output()?;
    assert!(!empty_threshold_report_output.status.success());
    assert!(
        String::from_utf8_lossy(&empty_threshold_report_output.stderr)
            .contains("No requested threshold-counterterm data was found."),
        "threshold reports without decomposition data should fail clearly",
    );

    let missing_summary_output = std::process::Command::new("python3")
        .arg(&script)
        .arg(&amp_detailed_only_path)
        .arg("--output")
        .arg(output_dir.join("approach_missing_summary.pdf"))
        .arg("--threshold-report")
        .args(["--threshold-report-section", "summary"])
        .output()?;
    assert!(!missing_summary_output.status.success());
    assert!(
        String::from_utf8_lossy(&missing_summary_output.stderr)
            .contains("Summary pages require evaluation.threshold_counterterm_summary"),
        "detailed event components must not produce an empty aggregate-summary page",
    );

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
