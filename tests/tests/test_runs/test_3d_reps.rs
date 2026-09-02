use std::fs;

use gammalooprs::processes::ProcessCollection;
use three_dimensional_reps::{
    EnergyEdgeIndexMap, OrientationID, ThreeDExpression, ThreeDGraphSource,
    eval::{EvaluationInput, evaluate_expression},
};

use super::*;

#[test]
#[serial]
#[allow(clippy::type_complexity)]
fn diagnostic_shifted_double_poles_validate_and_build_with_exact_energy_bounds() -> Result<()> {
    let test_name = "threedreps_diagnostic_shifted_double_poles";
    let test_root = get_tests_workspace_path().join(test_name);
    let state_path = test_root.join("state");
    let workspace_path = test_root.join("workspace");
    clean_test(&test_root);

    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;

    let cases: [(&str, &str, &[(usize, usize)], &[&[usize]]); 4] = [
        (
            "qsq",
            "diagnostic_shifted_double_pole_qsq.dot",
            &[(2, 2)],
            &[&[0, 1]],
        ),
        (
            "affine",
            "diagnostic_shifted_double_pole_affine.dot",
            &[(2, 1), (4, 1)],
            &[&[0, 1]],
        ),
        (
            "factorized",
            "diagnostic_shifted_double_pole_factorized.dot",
            &[(2, 1), (3, 1)],
            &[&[0, 1]],
        ),
        (
            "post_t",
            "diagnostic_post_t_disconnected_poles.dot",
            &[(2, 2), (5, 3)],
            &[&[0, 1], &[3, 4, 5]],
        ),
    ];
    for (integrand, graph_file, _, _) in cases {
        cli.run_command(&format!(
            "import graphs ./tests/resources/graphs/uv_tests/{graph_file} -p diagnostic_shifted_double_poles -i {integrand}"
        ))?;
    }

    let mut qsq = None;
    let mut affine = None;
    let mut factorized = None;
    for (integrand, _, expected_bounds, expected_repeated_groups) in cases {
        let validate_path = test_root.join(format!("{integrand}_validate.json"));
        cli.run_command(&format!(
            "3Drep validate -p diagnostic_shifted_double_poles -i {integrand} -g 0 --json-out {}",
            validate_path.display()
        ))?;
        let validation = serde_json::from_str::<JsonValue>(&fs::read_to_string(validate_path)?)?;
        assert_eq!(validation["validation"]["ok"].as_bool(), Some(true));
        assert_eq!(
            validation["validation"]["repeated_groups"],
            serde_json::to_value(expected_repeated_groups)?
        );

        let build_path = test_root.join(format!("{integrand}_build.json"));
        cli.run_command(&format!(
            "3Drep build -p diagnostic_shifted_double_poles -i {integrand} -g 0 --workspace-path {} --json-out {} --no-pretty",
            workspace_path.display(),
            build_path.display()
        ))?;
        let artifact = serde_json::from_str::<JsonValue>(&fs::read_to_string(build_path)?)?;
        assert_eq!(artifact["validation"]["ok"].as_bool(), Some(true));
        assert_eq!(
            artifact["validation"]["repeated_groups"],
            serde_json::to_value(expected_repeated_groups)?
        );
        assert_eq!(
            artifact["energy_degree_bounds"],
            serde_json::to_value(expected_bounds)?
        );

        if matches!(integrand, "qsq" | "affine" | "factorized") {
            let process_ref =
                ProcessRef::Unqualified("diagnostic_shifted_double_poles".to_string());
            let integrand_name = integrand.to_string();
            let (process_id, integrand_name) = cli
                .state
                .find_integrand_ref(Some(&process_ref), Some(&integrand_name))?;
            let ProcessCollection::Amplitudes(amplitudes) =
                &cli.state.process_list.processes[process_id].collection
            else {
                panic!("the shifted-double-pole DOTs must import as amplitudes")
            };
            let graph = &amplitudes[&integrand_name].graphs[0].graph;
            let parsed = graph.to_three_d_parsed_graph()?;
            let edge_map = graph
                .energy_edge_index_map(&parsed)
                .expect("GammaLoop Graph sources expose their physical edge IDs");
            let local_edge_map = EnergyEdgeIndexMap {
                internal: edge_map
                    .internal
                    .into_iter()
                    .map(|(local, physical)| (physical, local))
                    .collect(),
                external: edge_map
                    .external
                    .into_iter()
                    .map(|(local, physical)| (physical, local))
                    .collect(),
                orientation_edge_count: parsed.internal_edges.len(),
            };
            let expression = serde_json::from_value::<ThreeDExpression<OrientationID>>(
                artifact["expression"].clone(),
            )?
            .remap_energy_edge_indices(&local_edge_map);
            match integrand {
                "qsq" => qsq = Some((parsed, expression)),
                "affine" => affine = Some((parsed, expression)),
                "factorized" => factorized = Some((parsed, expression)),
                _ => unreachable!(),
            }
        }
    }

    let (qsq_parsed, qsq_expression) = qsq.expect("the qsq CLI artifact was built");
    let (affine_parsed, affine_expression) = affine.expect("the affine CLI artifact was built");
    let (factorized_parsed, factorized_expression) =
        factorized.expect("the factorized CLI artifact was built");
    assert_eq!(
        qsq_parsed, affine_parsed,
        "the qsq and affine DOT numerators must share one exact parsed topology"
    );
    assert_eq!(
        qsq_parsed, factorized_parsed,
        "the qsq and locally factorized DOT numerators must share one exact parsed topology"
    );
    for seed in [17, 1337, 9100] {
        let input = EvaluationInput::deterministic(&qsq_parsed, seed, &Default::default(), None)?;
        let direct =
            evaluate_expression(&qsq_parsed, &qsq_expression, "edges[0][0]**2", &input)?.value;
        let shifted = evaluate_expression(
            &qsq_parsed,
            &affine_expression,
            "edges[0][0]*(edges[2][0]+ext[0][0])",
            &input,
        )?
        .value;
        let split = evaluate_expression(
            &qsq_parsed,
            &factorized_expression,
            "edges[0][0]*edges[1][0]",
            &input,
        )?
        .value;
        assert!(
            [direct, shifted, split]
                .into_iter()
                .all(|value| value.is_finite() && value != 0.0),
            "the literal CLI comparison must be finite and nonzero at seed {seed}: direct={direct:e}, shifted={shifted:e}, split={split:e}"
        );
        let scale = direct.abs().max(shifted.abs()).max(split.abs());
        for (comparison, left, right) in [
            ("qsq vs affine", direct, shifted),
            ("qsq vs factorized", direct, split),
            ("affine vs factorized", shifted, split),
        ] {
            assert!(
                (left - right).abs() <= 1.0e-11 * scale,
                "{comparison} CLI artifacts differ at seed {seed}: left={left:e}, right={right:e}"
            );
        }
    }

    clean_test(test_root);
    Ok(())
}

#[test]
#[serial]
fn cff_cli_validate_and_build_use_gammaloop_graph_state() -> Result<()> {
    let test_name = "threedreps_cff_cli_validate_build";
    let test_root = get_tests_workspace_path().join(test_name);
    let state_path = test_root.join("state");
    let workspace_path = test_root.join("workspace");
    let validate_path = test_root.join("validate.json");
    clean_test(&test_root);

    let mut cli = get_test_cli(None, &state_path, Some(test_name.to_string()), true)?;
    cli.run_command("import model scalars-default.json")?;
    cli.run_command(
        "import graphs ./tests/resources/graphs/scalar_box.dot -p threedreps_box_build -i default -o",
    )?;
    cli.run_command("set global kv global.generation.uniform_numerator_sampling_scale=none")?;

    cli.run_command(&format!(
        "3Drep validate -p threedreps_box_build -i default -g 0 --json-out {}",
        validate_path.display()
    ))?;
    let validation = serde_json::from_str::<JsonValue>(&fs::read_to_string(&validate_path)?)?;
    assert_eq!(validation["graph"]["n_internal_edges"].as_u64(), Some(4));
    assert_eq!(
        validation["graph"]["loop_names"].as_array().map(Vec::len),
        Some(1)
    );
    assert_eq!(validation["validation"]["ok"].as_bool(), Some(true));

    cli.run_command(&format!(
        "3Drep build -p threedreps_box_build -i default -g 0 --numerator-samples-normalization M_for_all --workspace-path {} --no-pretty",
        workspace_path.display()
    ))?;
    let pointer = fs::read_to_string(workspace_path.join("latest_oriented_expression_path.txt"))?;
    let expression_path = PathBuf::from(pointer.trim());
    let expression_path = if expression_path.is_absolute() {
        expression_path
    } else {
        env::current_dir()?.join(expression_path)
    };
    let artifact = serde_json::from_str::<JsonValue>(&fs::read_to_string(expression_path)?)?;
    assert_eq!(artifact["backend"].as_str(), Some("gammaloop-3Drep"));
    assert_eq!(artifact["family"].as_str(), Some("cff"));
    assert_eq!(artifact["validation"]["ok"].as_bool(), Some(true));
    assert_eq!(
        artifact["numerator_energy_support"]["monomials"],
        serde_json::json!([[]])
    );
    assert_eq!(
        artifact["numerator_sampling_scale_mode"].as_str(),
        Some("All")
    );
    assert_eq!(
        artifact["expression"]["orientations"]
            .as_array()
            .map(Vec::len),
        Some(14)
    );
    assert_eq!(
        cli.cli_settings
            .global
            .generation
            .uniform_numerator_sampling_scale,
        gammalooprs::settings::global::UniformNumeratorSamplingScale::None
    );

    let ltd_error = cli
        .run_command(
            "3Drep build -p threedreps_box_build -i default -g 0 --representation ltd --no-save-json --no-pretty",
        )
        .unwrap_err();
    assert!(
        format!("{ltd_error:?}").contains(
            "three-dimensional representation mode Ltd is not implemented; only CFF is currently supported"
        )
    );

    clean_test(test_root);
    Ok(())
}
