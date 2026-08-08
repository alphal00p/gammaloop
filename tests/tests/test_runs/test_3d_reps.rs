use std::fs;

use super::*;

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
