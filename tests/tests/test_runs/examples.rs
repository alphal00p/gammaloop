use super::*;

#[test]
fn test_simple_workflow_example_card() -> Result<()> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("tests crate should live in the workspace root");
    let mut run_history =
        RunHistory::load(workspace_root.join("examples/command_cards/simple_workflow.toml"))?;
    let state_path = get_tests_workspace_path().join("simple_workflow_example");

    clean_test(&state_path);

    let mut state = new_cli_for_test(&state_path, None);
    run_history.cli_settings.state.folder = state_path.clone();
    let mut cli_settings = run_history.cli_settings.clone();
    cli_settings.sync_settings()?;
    let mut default_runtime_settings = run_history.default_runtime_settings.clone();

    let _ = run_history.run(&mut state, &mut cli_settings, &mut default_runtime_settings)?;

    assert!(
        !state.process_list.processes.is_empty(),
        "No processes were generated"
    );
    assert_eq!(state.model.name, "sm", "Expected SM model to be loaded");
    assert!(
        state_path.join("justfile").exists(),
        "Expected dot export helper files to be created in {}",
        state_path.display()
    );

    clean_test(&state_path);
    Ok(())
}

#[test]
fn test_advanced_integration_example_card() -> Result<()> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("tests crate should live in the workspace root");
    let run_history =
        RunHistory::load(workspace_root.join("examples/command_cards/advanced_integration.toml"))?;

    assert_eq!(
        run_history.commands[0].raw_string.as_deref(),
        Some("import model sm-default")
    );
    assert_eq!(
        run_history.cli_settings.state.folder,
        PathBuf::from("./advanced_integration_state")
    );
    assert_eq!(run_history.default_runtime_settings.kinematics.e_cm, 500.0);
    Ok(())
}
