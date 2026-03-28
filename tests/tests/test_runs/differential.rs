use super::utils::*;
use super::*;

fn configure_differential_leading_jet_observable(
    cli: &mut gammaloop_integration_tests::CLIState,
) -> Result<()> {
    run_commands(
        cli,
        &[
            r#"set process string '
[quantities.leading_jet_pt]
type = "jet"
quantity = "PT"
dR = 0.4
'"#,
            r#"set process string '
[observables.leading_jet_pt_hist]
quantity = "leading_jet_pt"
entry_selection = "leading_only"
x_min = 0.0
x_max = 1000.0
n_bins = 8
'"#,
        ],
    )
}

fn configure_differential_leading_jet_selector(
    cli: &mut gammaloop_integration_tests::CLIState,
) -> Result<()> {
    cli.run_command(
        r#"set process string '
[selectors.leading_jet_pt_cut]
quantity = "leading_jet_pt"
selector = "value_range"
entry_selection = "leading_only"
min = 0.0
'"#,
    )
}

#[test]
#[serial]
fn lu_differential_integration_writes_json_observables() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_differential_integration_json")?;
    configure_differential_leading_jet_observable(&mut cli)?;
    configure_differential_leading_jet_selector(&mut cli)?;
    cli.run_command(
        "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=json",
    )?;

    let workspace = get_tests_workspace_path().join("lu_differential_integration_json/workspace");
    let integration_result = Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    assert!(single_slot_integral(&integration_result).neval > 0);
    assert!(workspace.join("integration_result.json").exists());
    let slot_workspace = selected_slot_workspace(&cli, &workspace, None, Some("default"))?;
    let final_file = slot_workspace.join("observables_final.json");
    assert!(final_file.exists());
    assert!(
        !slot_workspace
            .join("observables_final_iter_0001.json")
            .exists()
    );

    let final_bundle =
        gammalooprs::observables::ObservableSnapshotBundle::from_json_file(&final_file)?;
    assert!(final_bundle.histograms.contains_key("leading_jet_pt_hist"));

    Ok(())
}

#[test]
#[serial]
fn lu_differential_observables_without_selectors_still_fill_histograms() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_differential_integration_json")?;
    configure_differential_leading_jet_observable(&mut cli)?;
    cli.run_command(
        "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=json",
    )?;

    let workspace =
        get_tests_workspace_path().join("lu_differential_integration_json/no_selector_workspace");
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let slot_workspace = selected_slot_workspace(&cli, &workspace, None, Some("default"))?;
    let final_file = slot_workspace.join("observables_final.json");
    let final_bundle =
        gammalooprs::observables::ObservableSnapshotBundle::from_json_file(&final_file)?;
    let histogram = final_bundle
        .histograms
        .get("leading_jet_pt_hist")
        .expect("missing leading-jet histogram");
    assert!(histogram.sample_count > 0);

    Ok(())
}

#[test]
#[serial]
fn lu_differential_integration_cli_flag_writes_iteration_observables() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_differential_integration_json")?;
    configure_differential_leading_jet_observable(&mut cli)?;
    configure_differential_leading_jet_selector(&mut cli)?;
    cli.run_command(
        "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=json",
    )?;

    let workspace =
        get_tests_workspace_path().join("lu_differential_integration_json/cli_flag_workspace");
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        write_results_for_each_iteration: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let slot_workspace = selected_slot_workspace(&cli, &workspace, None, Some("default"))?;
    assert!(workspace.join("integration_result.json").exists());
    assert!(
        workspace
            .join("results/integration_result_iter_0001.json")
            .exists()
    );
    assert!(
        slot_workspace
            .join("observables_final_iter_0001.json")
            .exists()
    );
    assert!(slot_workspace.join("observables_final.json").exists());

    Ok(())
}

#[test]
#[serial]
fn lu_differential_integration_hwu_output_is_optional_and_single_file() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_differential_integration_hwu")?;
    cli.run_command(
        "set process kv integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=hwu",
    )?;

    let workspace_without_observables =
        get_tests_workspace_path().join("lu_differential_integration_hwu/without_observables");
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace_without_observables.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        workspace_without_observables
            .join("integration_result.json")
            .exists()
    );
    assert!(
        !selected_slot_workspace(&cli, &workspace_without_observables, None, Some("default"))?
            .join("observables_final.hwu")
            .exists()
    );

    configure_differential_leading_jet_observable(&mut cli)?;
    configure_differential_leading_jet_selector(&mut cli)?;
    let workspace_with_observables =
        get_tests_workspace_path().join("lu_differential_integration_hwu/with_observables");
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace_with_observables.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        workspace_with_observables
            .join("integration_result.json")
            .exists()
    );

    let slot_workspace =
        selected_slot_workspace(&cli, &workspace_with_observables, None, Some("default"))?;
    let hwu_file = slot_workspace.join("observables_final.hwu");
    assert!(hwu_file.exists());
    let hwu_contents = std::fs::read_to_string(hwu_file)?;
    assert!(hwu_contents.contains("<histogram>"));
    assert!(hwu_contents.contains("leading_jet_pt_hist"));
    assert!(
        !slot_workspace
            .join("observables_final_iter_0001.hwu")
            .exists()
    );

    Ok(())
}

#[test]
#[serial]
fn lu_differential_json_observables_resume_from_workspace() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_differential_integration_json")?;
    configure_differential_leading_jet_observable(&mut cli)?;
    configure_differential_leading_jet_selector(&mut cli)?;
    cli.run_command(
        "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=json",
    )?;

    let workspace =
        get_tests_workspace_path().join("lu_differential_integration_json/resume_workspace");
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let slot_workspace = selected_slot_workspace(&cli, &workspace, None, Some("default"))?;
    let (process_id, resolved_integrand_name) = cli
        .state
        .find_integrand_ref(None, Some(&"default".to_string()))?;
    let slot_meta = gammalooprs::integrate::SlotMeta {
        process_name: cli.state.process_list.processes[process_id]
            .definition
            .folder_name
            .clone(),
        integrand_name: resolved_integrand_name,
    };
    let final_file = slot_workspace.join("observables_final.json");
    let checkpoint_file =
        gammalooprs::integrate::latest_observable_resume_state_path(&workspace, &slot_meta);
    assert!(checkpoint_file.exists());
    let result_snapshot_path = workspace.join("integration_result.json");
    let result_before_resume = std::fs::read_to_string(&result_snapshot_path)?;
    let state_before_resume =
        std::fs::read(gammalooprs::integrate::workspace_state_path(&workspace))?;
    let before_resume =
        gammalooprs::observables::ObservableSnapshotBundle::from_json_file(&final_file)?;

    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: false,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let after_resume =
        gammalooprs::observables::ObservableSnapshotBundle::from_json_file(&final_file)?;
    assert_eq!(before_resume, after_resume);
    assert_eq!(
        result_before_resume,
        std::fs::read_to_string(result_snapshot_path)?
    );
    assert_eq!(
        state_before_resume,
        std::fs::read(gammalooprs::integrate::workspace_state_path(&workspace))?
    );

    Ok(())
}

#[test]
#[serial]
fn lu_differential_hwu_observables_resume_from_workspace() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_differential_integration_hwu")?;
    configure_differential_leading_jet_observable(&mut cli)?;
    configure_differential_leading_jet_selector(&mut cli)?;
    cli.run_command(
        "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=hwu",
    )?;

    let workspace =
        get_tests_workspace_path().join("lu_differential_integration_hwu/resume_workspace");
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let slot_workspace = selected_slot_workspace(&cli, &workspace, None, Some("default"))?;
    let final_file = slot_workspace.join("observables_final.hwu");
    let before_resume = std::fs::read_to_string(&final_file)?;
    let (process_id, resolved_integrand_name) = cli
        .state
        .find_integrand_ref(None, Some(&"default".to_string()))?;
    let slot_meta = gammalooprs::integrate::SlotMeta {
        process_name: cli.state.process_list.processes[process_id]
            .definition
            .folder_name
            .clone(),
        integrand_name: resolved_integrand_name,
    };
    let checkpoint_file =
        gammalooprs::integrate::latest_observable_resume_state_path(&workspace, &slot_meta);
    assert!(checkpoint_file.exists());
    let result_snapshot_path = workspace.join("integration_result.json");
    let result_before_resume = std::fs::read_to_string(&result_snapshot_path)?;
    let state_before_resume =
        std::fs::read(gammalooprs::integrate::workspace_state_path(&workspace))?;

    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: false,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let after_resume = std::fs::read_to_string(final_file)?;
    assert_eq!(before_resume, after_resume);
    assert_eq!(
        result_before_resume,
        std::fs::read_to_string(result_snapshot_path)?
    );
    assert_eq!(
        state_before_resume,
        std::fs::read(gammalooprs::integrate::workspace_state_path(&workspace))?
    );

    Ok(())
}

#[test]
#[serial]
fn lu_save_dot_silently_overwrites_existing_files() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_save_dot_overwrite")?;
    cli.run_command("save dot")?;
    cli.run_command("save dot")?;
    Ok(())
}
