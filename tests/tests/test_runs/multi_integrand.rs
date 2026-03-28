use super::utils::*;
use super::*;

#[test]
fn test_multi_integrand() -> Result<()> {
    let test_name = "test_multi_integrand";
    let mut cli = setup_scalar_topologies_cli(test_name)?;
    let shared_integrator_settings = "kv integrator.n_start=5000 integrator.n_max=10000 integrator.n_increase=5000 integrator.seed=1337";

    cli.run_command(&format!(
        "set process -p triangle -i scalar_tri {shared_integrator_settings}"
    ))?;
    let triangle_baseline = scalar_topology_integrate_command(
        test_name,
        "triangle_baseline",
        &[("triangle", "scalar_tri")],
        &[],
    )
    .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    let box_above_baseline = scalar_topology_integrate_command(
        test_name,
        "box_above_baseline",
        &[("box", "scalar_box")],
        &[],
    )
    .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command(SCALAR_BOX_BELOW_EXTERNALS)?;
    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    let box_below_baseline = scalar_topology_integrate_command(
        test_name,
        "box_below_baseline",
        &[("box", "scalar_box")],
        &[],
    )
    .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command(SCALAR_BOX_ABOVE_EXTERNALS)?;

    let triangle_target = single_slot_integral(&triangle_baseline).result;
    let box_above_target = single_slot_integral(&box_above_baseline).result;
    let box_below_target = single_slot_integral(&box_below_baseline).result;

    cli.run_command(&format!(
        "set process -p triangle -i scalar_tri {shared_integrator_settings}"
    ))?;
    let triangle_box_workspace = get_tests_workspace_path()
        .join(test_name)
        .join("triangle_and_box_correlated");
    scalar_topology_integrate_command(
        test_name,
        "triangle_and_box_correlated",
        &[("triangle", "scalar_tri"), ("box", "scalar_box")],
        &[
            ("triangle@scalar_tri", triangle_target),
            ("box@scalar_box", box_above_target),
        ],
    )
    .run(&mut cli.state, &cli.cli_settings)?;
    let triangle_box_result =
        load_integration_result(&triangle_box_workspace.join("integration_result.json"))?;

    assert_eq!(triangle_box_result.slots.len(), 2);
    let triangle_slot = triangle_box_result
        .slot("triangle@scalar_tri")
        .expect("triangle slot must be present");
    let box_slot = triangle_box_result
        .slot("box@scalar_box")
        .expect("box slot must be present");
    assert_eq!(triangle_slot.target, Some(triangle_target));
    assert_eq!(box_slot.target, Some(box_above_target));
    assert_eq!(triangle_slot.table_results.len(), 2);
    assert_eq!(box_slot.table_results.len(), 2);
    assert!(triangle_slot.integration_statistics.num_evals > 0);
    assert!(box_slot.integration_statistics.num_evals > 0);
    assert!(!triangle_slot.max_weight_info.is_empty());
    assert!(!box_slot.max_weight_info.is_empty());
    assert!(
        triangle_slot
            .integral
            .is_compatible_with_result(single_slot_integral(&triangle_baseline), 4),
        "triangle correlated result deviates from baseline beyond 4 sigma: {} vs {}",
        triangle_slot.integral,
        single_slot_integral(&triangle_baseline),
    );
    assert!(
        box_slot
            .integral
            .is_compatible_with_result(single_slot_integral(&box_above_baseline), 4),
        "box correlated result deviates from baseline beyond 4 sigma: {} vs {}",
        box_slot.integral,
        single_slot_integral(&box_above_baseline),
    );
    assert!(
        triangle_box_workspace
            .join("integrands/triangle@scalar_tri/settings.toml")
            .exists()
    );
    assert!(
        triangle_box_workspace
            .join("integrands/box@scalar_box/settings.toml")
            .exists()
    );

    cli.run_command(
        "duplicate integrand -p box -i scalar_box --output_process_name box_copy --output_integrand_name scalar_box_copy",
    )?;
    cli.run_command(SCALAR_BOX_COPY_BELOW_EXTERNALS)?;
    let output_process = ProcessRef::Unqualified("box_copy".to_string());
    let output_integrand = "scalar_box_copy".to_string();
    cli.state
        .find_integrand_ref(Some(&output_process), Some(&output_integrand))?;

    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    let box_pair_workspace = get_tests_workspace_path()
        .join(test_name)
        .join("box_above_and_below_correlated");
    scalar_topology_integrate_command(
        test_name,
        "box_above_and_below_correlated",
        &[("box", "scalar_box"), ("box_copy", "scalar_box_copy")],
        &[
            ("box@scalar_box", box_above_target),
            ("box_copy@scalar_box_copy", box_below_target),
        ],
    )
    .run(&mut cli.state, &cli.cli_settings)?;
    let box_pair_result =
        load_integration_result(&box_pair_workspace.join("integration_result.json"))?;

    assert_eq!(box_pair_result.slots.len(), 2);
    let box_above_slot = box_pair_result
        .slot("box@scalar_box")
        .expect("box slot must be present");
    let box_below_slot = box_pair_result
        .slot("box_copy@scalar_box_copy")
        .expect("box_copy slot must be present");
    assert_eq!(box_above_slot.target, Some(box_above_target));
    assert_eq!(box_below_slot.target, Some(box_below_target));
    assert_eq!(box_above_slot.table_results.len(), 2);
    assert_eq!(box_below_slot.table_results.len(), 2);
    assert!(box_above_slot.integration_statistics.num_evals > 0);
    assert!(box_below_slot.integration_statistics.num_evals > 0);
    assert!(!box_above_slot.max_weight_info.is_empty());
    assert!(!box_below_slot.max_weight_info.is_empty());
    assert!(
        box_above_slot
            .integral
            .is_compatible_with_result(single_slot_integral(&box_above_baseline), 4),
        "leading box correlated result deviates from baseline beyond 4 sigma: {} vs {}",
        box_above_slot.integral,
        single_slot_integral(&box_above_baseline),
    );
    assert!(
        box_below_slot
            .integral
            .is_compatible_with_result(single_slot_integral(&box_below_baseline), 4),
        "secondary box correlated result deviates from baseline beyond 4 sigma: {} vs {}",
        box_below_slot.integral,
        single_slot_integral(&box_below_baseline),
    );
    assert!(
        box_pair_workspace
            .join("integrands/box@scalar_box/settings.toml")
            .exists()
    );
    assert!(
        box_pair_workspace
            .join("integrands/box_copy@scalar_box_copy/settings.toml")
            .exists()
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn test_multi_integrand_batching_preserves_results() -> Result<()> {
    let test_name = "test_multi_integrand_batching_preserves_results";
    let mut cli = setup_scalar_topologies_cli(test_name)?;
    let shared_integrator_settings = "kv integrator.n_start=5000 integrator.n_max=10000 integrator.n_increase=5000 integrator.seed=1337";

    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    cli.run_command(SCALAR_BOX_BELOW_EXTERNALS)?;

    let mut single_batch = scalar_topology_integrate_command(
        test_name,
        "box_single_batch",
        &[("box", "scalar_box")],
        &[],
    );
    single_batch.batch_size = Some(10_000);
    single_batch.no_stream_iterations = true;
    single_batch.no_stream_updates = true;
    single_batch.run(&mut cli.state, &cli.cli_settings)?;

    let mut many_batches = scalar_topology_integrate_command(
        test_name,
        "box_many_batches",
        &[("box", "scalar_box")],
        &[],
    );
    many_batches.batch_size = Some(1);
    many_batches.no_stream_iterations = true;
    many_batches.no_stream_updates = true;
    many_batches.run(&mut cli.state, &cli.cli_settings)?;

    let single_batch_result = load_integration_result(
        &get_tests_workspace_path()
            .join(test_name)
            .join("box_single_batch")
            .join("integration_result.json"),
    )?;
    let many_batches_result = load_integration_result(
        &get_tests_workspace_path()
            .join(test_name)
            .join("box_many_batches")
            .join("integration_result.json"),
    )?;

    assert_integration_results_match_ignoring_timings(&single_batch_result, &many_batches_result);

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn test_multi_integrand_with_local_model_parameters() -> Result<()> {
    let test_name = "test_multi_integrand_with_local_model_parameters";
    let mut cli = setup_scalar_topologies_cli(test_name)?;
    let shared_integrator_settings = "kv integrator.n_start=5000 integrator.n_max=10000 integrator.n_increase=5000 integrator.seed=1337";

    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    let box_default_baseline = scalar_topology_integrate_command(
        test_name,
        "box_default_model",
        &[("box", "scalar_box")],
        &[],
    )
    .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command(
        "duplicate integrand -p box -i scalar_box --output_process_name box_copy --output_integrand_name scalar_box_copy",
    )?;
    cli.run_command("set model -p box_copy mass_scalar_2=1.0")?;
    cli.run_command(&format!(
        "set process -p box_copy -i scalar_box_copy {shared_integrator_settings}"
    ))?;

    let box_process = ProcessRef::Unqualified("box".to_string());
    let box_integrand = "scalar_box".to_string();
    let (box_process_id, box_integrand_name) = cli
        .state
        .find_integrand_ref(Some(&box_process), Some(&box_integrand))?;
    let box_card = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(
            box_process_id,
            &box_integrand_name,
        )?;
    assert_eq!(box_card.data.get("mass_scalar_2"), Some(&(F(2.0), F(0.0))));

    let box_copy_process = ProcessRef::Unqualified("box_copy".to_string());
    let box_copy_integrand = "scalar_box_copy".to_string();
    let (box_copy_process_id, box_copy_integrand_name) = cli
        .state
        .find_integrand_ref(Some(&box_copy_process), Some(&box_copy_integrand))?;
    let box_copy_card = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(
            box_copy_process_id,
            &box_copy_integrand_name,
        )?;
    assert_eq!(
        box_copy_card.data.get("mass_scalar_2"),
        Some(&(F(1.0), F(0.0)))
    );

    let box_local_baseline = scalar_topology_integrate_command(
        test_name,
        "box_local_model",
        &[("box_copy", "scalar_box_copy")],
        &[],
    )
    .run(&mut cli.state, &cli.cli_settings)?;

    let box_default_target = single_slot_integral(&box_default_baseline).result;
    let box_local_target = single_slot_integral(&box_local_baseline).result;

    let shared_workspace = get_tests_workspace_path()
        .join(test_name)
        .join("box_default_and_local_model");
    scalar_topology_integrate_command(
        test_name,
        "box_default_and_local_model",
        &[("box", "scalar_box"), ("box_copy", "scalar_box_copy")],
        &[
            ("box@scalar_box", box_default_target),
            ("box_copy@scalar_box_copy", box_local_target),
        ],
    )
    .run(&mut cli.state, &cli.cli_settings)?;
    let multi_result = load_integration_result(&shared_workspace.join("integration_result.json"))?;

    let box_default_slot = multi_result
        .slot("box@scalar_box")
        .expect("box slot must be present");
    let box_local_slot = multi_result
        .slot("box_copy@scalar_box_copy")
        .expect("box_copy slot must be present");
    assert!(
        box_default_slot
            .integral
            .is_compatible_with_result(single_slot_integral(&box_default_baseline), 4),
        "default-model box correlated result deviates from baseline beyond 4 sigma: {} vs {}",
        box_default_slot.integral,
        single_slot_integral(&box_default_baseline),
    );
    assert!(
        box_local_slot
            .integral
            .is_compatible_with_result(single_slot_integral(&box_local_baseline), 4),
        "local-model box correlated result deviates from baseline beyond 4 sigma: {} vs {}",
        box_local_slot.integral,
        single_slot_integral(&box_local_baseline),
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn test_inspect_uses_per_integrand_model_parameters() -> Result<()> {
    let test_name = "test_inspect_uses_per_integrand_model_parameters";
    let mut cli = setup_scalar_topologies_cli(test_name)?;

    cli.run_command(
        "duplicate integrand -p box -i scalar_box --output_process_name box_copy --output_integrand_name scalar_box_copy",
    )?;

    let inspect_point = vec![0.31, 0.52, 0.73];
    let (_, default_inspect) = Inspect {
        process: Some(ProcessRef::Unqualified("box".to_string())),
        integrand_name: Some("scalar_box".to_string()),
        point: inspect_point.clone(),
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;
    let (_, copied_inspect_before_override) = Inspect {
        process: Some(ProcessRef::Unqualified("box_copy".to_string())),
        integrand_name: Some("scalar_box_copy".to_string()),
        point: inspect_point.clone(),
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;

    assert!(
        complex_distance(default_inspect, copied_inspect_before_override) < 1e-12,
        "duplicated integrand should initially inspect identically: {} vs {}",
        default_inspect,
        copied_inspect_before_override,
    );

    cli.run_command("set model -p box_copy mass_scalar_2=1.0")?;

    let (_, copied_inspect_after_override) = Inspect {
        process: Some(ProcessRef::Unqualified("box_copy".to_string())),
        integrand_name: Some("scalar_box_copy".to_string()),
        point: inspect_point.clone(),
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;

    assert!(
        complex_distance(default_inspect, copied_inspect_after_override) > 1e-12,
        "per-integrand model override should change the inspect result: {} vs {}",
        default_inspect,
        copied_inspect_after_override,
    );

    let box_copy_process = ProcessRef::Unqualified("box_copy".to_string());
    let box_copy_integrand = "scalar_box_copy".to_string();
    let (box_copy_process_id, box_copy_integrand_name) = cli
        .state
        .find_integrand_ref(Some(&box_copy_process), Some(&box_copy_integrand))?;
    let box_copy_card = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(
            box_copy_process_id,
            &box_copy_integrand_name,
        )?;
    assert_eq!(
        box_copy_card.data.get("mass_scalar_2"),
        Some(&(F(1.0), F(0.0)))
    );

    cli.run_command("set process -p box_copy -i scalar_box_copy defaults")?;

    let box_copy_reset_card = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(
            box_copy_process_id,
            &box_copy_integrand_name,
        )?;
    assert_eq!(
        box_copy_reset_card.data.get("mass_scalar_2"),
        Some(&(F(2.0), F(0.0)))
    );

    cli.run_command(SCALAR_BOX_COPY_ABOVE_EXTERNALS)?;

    let (_, copied_inspect_after_defaults) = Inspect {
        process: Some(ProcessRef::Unqualified("box_copy".to_string())),
        integrand_name: Some("scalar_box_copy".to_string()),
        point: inspect_point,
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;

    assert!(
        complex_distance(default_inspect, copied_inspect_after_defaults) < 1e-12,
        "resetting process defaults should clear the local model override: {} vs {}",
        default_inspect,
        copied_inspect_after_defaults,
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn test_integration_workspace_model_mismatch_requires_restart() -> Result<()> {
    let test_name = "test_integration_workspace_model_mismatch_requires_restart";
    let mut cli = setup_scalar_topologies_cli(test_name)?;
    let shared_integrator_settings = "kv integrator.n_start=5000 integrator.n_max=10000 integrator.n_increase=5000 integrator.seed=1337";

    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    let workspace_name = "box_workspace_model_mismatch";
    scalar_topology_integrate_command(test_name, workspace_name, &[("box", "scalar_box")], &[])
        .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command("set model mass_scalar_2=1.0")?;
    let mut resume_command =
        scalar_topology_integrate_command(test_name, workspace_name, &[("box", "scalar_box")], &[]);
    resume_command.restart = false;
    let err = resume_command
        .run(&mut cli.state, &cli.cli_settings)
        .expect_err("workspace resume should fail on model-parameter mismatch");
    let err_text = format!("{err:?}");
    assert!(err_text.contains("Workspace effective model parameters do not match"));
    assert!(err_text.contains("box@scalar_box"));

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn test_integration_workspace_resume_preserves_current_effective_model_override() -> Result<()> {
    let test_name = "test_integration_workspace_resume_preserves_current_effective_model_override";
    let mut cli = setup_scalar_topologies_cli(test_name)?;
    let shared_integrator_settings = "kv integrator.n_start=5000 integrator.n_max=10000 integrator.n_increase=5000 integrator.seed=1337";

    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    let workspace_name = "box_workspace_effective_model_match";
    scalar_topology_integrate_command(test_name, workspace_name, &[("box", "scalar_box")], &[])
        .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command("set model mass_scalar_2=1.0")?;
    cli.run_command("set process -p box -i scalar_box kv model.mass_scalar_2=2.0")?;

    let mut resume_command =
        scalar_topology_integrate_command(test_name, workspace_name, &[("box", "scalar_box")], &[]);
    resume_command.restart = false;
    resume_command.run(&mut cli.state, &cli.cli_settings)?;

    let (process_id, integrand_name) = cli.state.find_integrand_ref(
        Some(&ProcessRef::Unqualified("box".to_string())),
        Some(&"scalar_box".to_string()),
    )?;
    let effective_card = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(process_id, &integrand_name)?;
    assert_eq!(
        effective_card.data.get("mass_scalar_2"),
        Some(&(F(2.0), F(0.0)))
    );

    let mut second_resume =
        scalar_topology_integrate_command(test_name, workspace_name, &[("box", "scalar_box")], &[]);
    second_resume.restart = false;
    second_resume.run(&mut cli.state, &cli.cli_settings)?;

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
