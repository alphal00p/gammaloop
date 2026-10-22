use super::utils::{selected_slot_workspace, single_slot_integral};
use super::*;
use gammaloop_integration_tests::CLIState;
use gammaloop_integration_tests::SERIAL_TEST_LOCK;
use gammalooprs::integrands::Integrand;
use gammalooprs::integrate;
use gammalooprs::integrate::{
    ContributionSortMode, HavanaIntegrateRequest, IntegrationCheckpoint, IntegrationSlot,
    IntegrationStatusKind, IntegrationStatusPhaseDisplay, IntegrationStatusViewOptions,
    IntegrationWorkspace, IterationBatchingSettings, SamplingCorrelationMode,
    WorkspaceSnapshotControl,
};
use gammalooprs::observables::ObservableSnapshotBundle;
use std::fs;

fn configure_differential_leading_jet_observable(cli: &mut CLIState) -> Result<()> {
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

fn configure_differential_leading_jet_selector(cli: &mut CLIState) -> Result<()> {
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
fn lu_differential_observables_without_selectors_still_fill_histograms() -> Result<()> {
    let _serial_guard = SERIAL_TEST_LOCK.lock();
    let mut cli = setup_sm_differential_lu_cli(
        "lu_differential_observables_without_selectors_still_fill_histograms",
    )?;
    configure_differential_leading_jet_observable(&mut cli)?;
    cli.run_command(
        "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=[\"json\"]",
    )?;

    let workspace = cli.cli_settings.state.folder.join("no_selector_workspace");
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
    let final_bundle = ObservableSnapshotBundle::from_json_file(&final_file)?;
    let histogram = final_bundle
        .histograms
        .get("leading_jet_pt_hist")
        .expect("missing leading-jet histogram");
    assert!(histogram.sample_count > 0);

    Ok(())
}

#[test]
fn lu_differential_integration_cli_flag_writes_iteration_observables() -> Result<()> {
    let _serial_guard = SERIAL_TEST_LOCK.lock();
    let mut cli = setup_sm_differential_lu_cli(
        "lu_differential_integration_cli_flag_writes_iteration_observables",
    )?;
    configure_differential_leading_jet_observable(&mut cli)?;
    configure_differential_leading_jet_selector(&mut cli)?;
    cli.run_command(
        "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=[\"json\"]",
    )?;

    let workspace = cli.cli_settings.state.folder.join("cli_flag_workspace");
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
fn lu_differential_integration_hwu_output_is_optional_and_single_file() -> Result<()> {
    let _serial_guard = SERIAL_TEST_LOCK.lock();
    let mut cli = setup_sm_differential_lu_cli(
        "lu_differential_integration_hwu_output_is_optional_and_single_file",
    )?;
    cli.run_command(
        "set process kv integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=[\"hwu\"]",
    )?;

    let workspace_without_observables = cli.cli_settings.state.folder.join("without_observables");
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
    let workspace_with_observables = cli.cli_settings.state.folder.join("with_observables");
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
    let hwu_contents = fs::read_to_string(hwu_file)?;
    assert!(hwu_contents.contains("<histogram>"));
    assert!(hwu_contents.contains("leading_jet_pt_hist"));
    assert!(
        !slot_workspace
            .join("observables_final_iter_0001.hwu")
            .exists()
    );
    assert!(!slot_workspace.join("observables_final.json").exists());
    assert!(
        !slot_workspace
            .join("observables_final_iter_0001.json")
            .exists()
    );

    Ok(())
}

#[test]
fn lu_save_dot_silently_overwrites_existing_files() -> Result<()> {
    let _serial_guard = SERIAL_TEST_LOCK.lock();
    let mut cli = setup_sm_differential_lu_cli("lu_save_dot_overwrite")?;
    cli.run_command("save dot")?;
    cli.run_command("save dot")?;
    Ok(())
}

mod failing {
    use super::*;

    #[test]
    fn lu_differential_integration_writes_json_observables() -> Result<()> {
        let _serial_guard = SERIAL_TEST_LOCK.lock();
        let mut cli =
            setup_sm_differential_lu_cli("lu_differential_integration_writes_json_observables")?;
        configure_differential_leading_jet_observable(&mut cli)?;
        configure_differential_leading_jet_selector(&mut cli)?;
        cli.run_command(
            "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=[\"json\"]",
        )?;

        let workspace = cli.cli_settings.state.folder.join("workspace");
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

        let final_bundle = ObservableSnapshotBundle::from_json_file(&final_file)?;
        assert!(final_bundle.histograms.contains_key("leading_jet_pt_hist"));

        Ok(())
    }
}

fn differential_observables_resume_from_nonzero_checkpoint(
    format: &str,
    test_name: &str,
) -> Result<()> {
    let _serial_guard = SERIAL_TEST_LOCK.lock();
    let mut cli = setup_sm_differential_lu_cli(test_name)?;
    configure_differential_leading_jet_observable(&mut cli)?;
    configure_differential_leading_jet_selector(&mut cli)?;
    cli.run_command(&format!(
        "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=24 integrator.n_increase=0 integrator.observables_output.format=[\"{format}\"]"
    ))?;
    let full_workspace = cli
        .cli_settings
        .state
        .folder
        .join("uninterrupted_workspace");
    let mut command = Integrate {
        integrand_name: vec!["default".into()],
        n_cores: Some(1),
        restart: true,
        workspace_path: Some(full_workspace.clone()),
        batch_size: Some(4),
        batch_timing: 0.0,
        ..Default::default()
    };
    let full = command.run(&mut cli.state, &cli.cli_settings)?;
    let full_checkpoint = IntegrationCheckpoint::open(&full_workspace)?;
    assert_eq!(full_checkpoint.read_state()?.iter, 2);
    assert_eq!(full_checkpoint.read_state()?.num_points, 24);

    let workspace = cli.cli_settings.state.folder.join("resume_workspace");
    let manifest = full_checkpoint.manifest().workspace.clone();
    let meta = manifest.slots[0].clone();
    let (process_id, integrand_name) = cli
        .state
        .find_integrand_ref(None, Some(&meta.integrand_name))?;
    let model = cli
        .state
        .resolve_model_for_integrand(process_id, &integrand_name)?;
    let settings = full_checkpoint.read_slot_settings(0)?;
    let process = cli
        .state
        .process_list
        .get_integrand_mut(process_id, &integrand_name)?
        .clone();
    let mut writer = IntegrationWorkspace::acquire(&workspace, manifest)?;
    let stopped = integrate::havana_integrate(
        HavanaIntegrateRequest {
            slots: vec![IntegrationSlot::new(
                meta,
                settings.clone(),
                model,
                Integrand::ProcessIntegrand(Box::new(process)),
                None,
            )],
            sampling_correlation_mode: SamplingCorrelationMode::Correlated,
            n_cores: 1,
            state: None,
            workspace: Some(&mut writer),
            output_control: WorkspaceSnapshotControl::default(),
            batching: IterationBatchingSettings {
                batch_size: Some(4),
                batch_timing_seconds: 0.0,
                ..Default::default()
            },
            view_options: IntegrationStatusViewOptions {
                phase_display: IntegrationStatusPhaseDisplay::Both,
                training_phase_display: IntegrationStatusPhaseDisplay::Real,
                training_slot: 0,
                slot_training_phase_displays: vec![IntegrationStatusPhaseDisplay::Real],
                per_slot_training_phase: false,
                target_relative_accuracy: settings.integrator.target_relative_accuracy,
                target_absolute_accuracy: settings.integrator.target_absolute_accuracy,
                show_statistics: false,
                show_max_weight_details: false,
                show_top_discrete_grid: false,
                show_discrete_contributions_sum: false,
                contribution_sort: ContributionSortMode::Index,
                show_max_weight_info_for_discrete_bins: false,
            },
        },
        |status| {
            if status.kind() == IntegrationStatusKind::Iteration {
                gammalooprs::request_interrupt();
            }
            Ok(())
        },
    );
    gammalooprs::clear_interrupt_request();
    stopped?;
    drop(writer);
    let previous = IntegrationCheckpoint::open(&workspace)?;
    assert_eq!(previous.read_state()?.iter, 1);
    assert_eq!(previous.read_state()?.num_points, 12);
    let first_histogram = previous.read_observables(0)?.unwrap();
    assert_eq!(
        first_histogram.histograms["leading_jet_pt_hist"].sample_count,
        12
    );

    // Presentation files deliberately disagree: neither these nor absent latest JSON/HwU
    // files may influence authoritative restoration.
    fs::write(
        workspace.join("integration_result.json"),
        b"invalid presentation JSON",
    )?;
    let slot_workspace = selected_slot_workspace(&cli, &workspace, None, Some("default"))?;
    fs::remove_file(slot_workspace.join(format!("observables_final.{format}")))?;
    command.workspace_path = Some(workspace.clone());
    command.restart = false;
    let resumed = command.run(&mut cli.state, &cli.cli_settings)?;
    let resumed_checkpoint = IntegrationCheckpoint::open(&workspace)?;
    assert_eq!(resumed_checkpoint.read_state()?.iter, 2);
    assert_eq!(resumed_checkpoint.read_state()?.num_points, 24);
    let expected = full_checkpoint.read_observables(0)?.unwrap();
    let actual = resumed_checkpoint.read_observables(0)?.unwrap();
    assert_eq!(actual.histograms, expected.histograms);
    assert_eq!(actual.histograms["leading_jet_pt_hist"].sample_count, 24);
    assert_eq!(
        resumed.single_slot_observables().unwrap().histograms,
        expected.histograms
    );
    assert_eq!(
        serde_json::to_value(&resumed.result.single_slot().unwrap().integral)?,
        serde_json::to_value(&full.result.single_slot().unwrap().integral)?
    );
    let normalize_state = |checkpoint: &IntegrationCheckpoint| -> Result<JsonValue> {
        let mut value = serde_json::to_value(checkpoint.read_state()?)?;
        value["elapsed_seconds"] = 0.0.into();
        let remove_timings = |stats: &mut JsonValue| {
            for key in [
                "sum_integrand_evaluation_time",
                "sum_evaluator_evaluation_time",
                "sum_parameterization_time",
                "sum_event_time",
                "sum_integrator_overhead_time",
                "sum_total_evaluation_time",
            ] {
                stats.as_object_mut().unwrap().remove(key);
            }
        };
        // Preserve every evaluation/stability/event count, scalar moment and trained grid.
        // Only the explicit wall-clock duration fields are excluded from exact equality.
        remove_timings(&mut value["stats"]);
        for stats in value["slot_stats"].as_array_mut().unwrap() {
            remove_timings(stats);
        }
        Ok(value)
    };
    assert_eq!(
        normalize_state(&resumed_checkpoint)?,
        normalize_state(&full_checkpoint)?
    );
    if format == "json" {
        assert_eq!(
            ObservableSnapshotBundle::from_json_file(
                slot_workspace.join("observables_final.json")
            )?
            .histograms,
            expected.histograms
        );
    } else {
        let full_slot = selected_slot_workspace(&cli, &full_workspace, None, Some("default"))?;
        assert_eq!(
            fs::read(slot_workspace.join("observables_final.hwu"))?,
            fs::read(full_slot.join("observables_final.hwu"))?
        );
    }

    // Exhausted resume performs no publication, including with history enabled.
    let descriptor = fs::read(workspace.join("manifest.json"))?;
    let state_bytes = fs::read(resumed_checkpoint.state_path())?;
    command.write_results_for_each_iteration = true;
    let no_work = command.run(&mut cli.state, &cli.cli_settings)?;
    assert_eq!(fs::read(workspace.join("manifest.json"))?, descriptor);
    assert_eq!(fs::read(resumed_checkpoint.state_path())?, state_bytes);
    assert_eq!(
        no_work.single_slot_observables().unwrap().histograms,
        expected.histograms
    );
    let summary_root = cli.cli_settings.state.folder.join("summary_only_workspace");
    let summary_generation = summary_root
        .join("generations")
        .join(resumed_checkpoint.manifest().generation.as_str());
    fs::create_dir_all(summary_generation.join("slots/0"))?;
    let copies = [
        (
            workspace.join("manifest.json"),
            summary_root.join("manifest.json"),
        ),
        (
            resumed_checkpoint.state_path(),
            summary_generation.join("integration_state.bin"),
        ),
        (
            resumed_checkpoint.slot_settings_path(0),
            summary_generation.join("slots/0/settings.toml"),
        ),
        (
            resumed_checkpoint.observable_path(0).unwrap(),
            summary_generation.join("slots/0/observables.json"),
        ),
    ];
    let mut before_summary = Vec::new();
    for (source, destination) in &copies {
        fs::copy(source, destination)?;
        before_summary.push(fs::read(destination)?);
    }
    command.workspace_path = Some(summary_root.clone());
    command.show_summary_only = true;
    let summary = command.run(&mut cli.state, &cli.cli_settings)?;
    assert_eq!(
        summary.single_slot_observables().unwrap().histograms,
        expected.histograms
    );
    for ((_, path), before) in copies.iter().zip(before_summary) {
        assert_eq!(fs::read(path)?, before);
    }
    assert!(!summary_root.join(".snapshot.lock").exists());
    assert!(!summary_root.join("integration_result.json").exists());
    Ok(())
}

#[test]
fn lu_differential_json_observables_resume_from_workspace() -> Result<()> {
    differential_observables_resume_from_nonzero_checkpoint(
        "json",
        "lu_differential_json_observables_resume_from_workspace",
    )
}

#[test]
fn lu_differential_hwu_observables_resume_from_workspace() -> Result<()> {
    differential_observables_resume_from_nonzero_checkpoint(
        "hwu",
        "lu_differential_hwu_observables_resume_from_workspace",
    )
}
