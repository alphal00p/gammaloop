use super::*;

#[test]
fn test_scalar_bubble_example_cli() -> Result<()> {
    let state_path = get_tests_workspace_path().join("scalar_bubble_example");
    let cli = get_example_cli(
        "scalar_topologies/bubble.toml",
        &["generate"],
        Some(state_path.clone()),
        None,
        true,
    )?;

    assert!(
        !cli.state.process_list.processes.is_empty(),
        "No processes were generated"
    );
    assert_eq!(
        cli.state.model.name, "scalars",
        "Expected scalars model to be loaded"
    );
    assert!(
        state_path.join("state_manifest.toml").exists(),
        "Expected saved state manifest in {}",
        state_path.display()
    );

    clean_test(&state_path);
    Ok(())
}

#[test]
fn test_aa_aa_example_card() -> Result<()> {
    let run_history = example_run_card("aa_aa/1L/aa_aa.toml")?;

    assert_eq!(
        run_history.command_blocks[0].commands[0]
            .raw_string
            .as_deref(),
        Some("import model sm-default.json")
    );
    assert_eq!(
        run_history.cli_settings.state.folder,
        PathBuf::from("./examples/cli/aa_aa/1L/state")
    );
    assert_eq!(run_history.default_runtime_settings.kinematics.e_cm, 300.0);
    Ok(())
}

#[test]
fn test_epem_a_tth_nlo_example_cli() -> Result<()> {
    let state_path = get_tests_workspace_path().join("epem_a_tth_nlo_example");
    let cli = get_example_cli(
        "epem_a_ttxh/NLO/epem_a_tth_NLO.toml",
        &["generate_diagrams"],
        Some(state_path.clone()),
        None,
        true,
    )?;

    assert_eq!(cli.state.model.name, "sm", "Expected SM model to be loaded");
    assert!(
        !cli.state.process_list.processes.is_empty(),
        "No processes were generated"
    );
    assert_eq!(
        cli.state.process_list.processes[0].definition.folder_name,
        "epem_a_tth"
    );
    assert!(
        state_path.join("justfile").exists(),
        "Expected dot export helper files to be created in {}",
        state_path.display()
    );

    clean_test(&state_path);
    Ok(())
}
