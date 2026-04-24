use super::*;
use gammaloop_integration_tests::run_card;

#[test]
fn massless_triangle_bulk_profile_passes() -> Result<()> {
    let state_path = get_tests_workspace_path().join("massless_triangle_bulk_profile");
    let run_history = run_card("generate_massless_triangle.toml")?;
    let setup_commands = run_history
        .commands
        .iter()
        .map(|command| {
            command
                .raw_string
                .as_deref()
                .expect("run-card commands loaded from strings should retain raw strings")
        })
        .take_while(|command| *command != "profile bulk")
        .collect_vec();

    assert!(
        run_history
            .commands
            .iter()
            .any(|command| command.raw_string.as_deref() == Some("profile bulk")),
        "expected tests/resources/run_cards/generate_massless_triangle.toml to contain the profile bulk command"
    );
    assert_eq!(setup_commands.last().copied(), Some("generate"));

    let mut cli = get_test_cli(None, &state_path, None, true)?;
    cli.cli_settings = run_history.cli_settings.clone();
    cli.cli_settings.state.folder = state_path.clone();
    cli.cli_settings.sync_settings()?;
    cli.default_runtime_settings = run_history.default_runtime_settings.clone();

    run_commands(&mut cli, &setup_commands)?;

    let profile_result = Profile::InfraRed(InfraRedProfile::default())
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_ir();

    assert!(profile_result.all_passed, "{profile_result}");

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn dotted_bubble_amp_bulk_profile_reports_threshold_component_scaling() -> Result<()> {
    let state_path = get_tests_workspace_path().join("dotted_bubble_amp_bulk_profile");
    let run_history = run_card("dotted_bubble_amp_generate.toml")?;
    let setup_commands = run_history
        .commands
        .iter()
        .map(|command| {
            command
                .raw_string
                .as_deref()
                .expect("run-card commands loaded from strings should retain raw strings")
        })
        .take_while(|command| *command != "profile bulk")
        .collect_vec();

    assert!(
        run_history
            .commands
            .iter()
            .any(|command| command.raw_string.as_deref() == Some("profile bulk")),
        "expected tests/resources/run_cards/dotted_bubble_amp_generate.toml to contain the profile bulk command"
    );
    assert_eq!(setup_commands.last().copied(), Some("generate"));

    let mut cli = get_test_cli(None, &state_path, None, true)?;
    cli.cli_settings = run_history.cli_settings.clone();
    cli.cli_settings.state.folder = state_path.clone();
    cli.cli_settings.sync_settings()?;
    cli.default_runtime_settings = run_history.default_runtime_settings.clone();

    run_commands(&mut cli, &setup_commands)?;

    let profile_result = Profile::InfraRed(InfraRedProfile::default())
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_ir();

    let graph_report = profile_result
        .results_per_graph
        .iter()
        .find(|report| report.graph_name == "dotted_bubble")
        .expect("expected dotted_bubble graph report");
    let threshold_report = graph_report
        .single_limit_reports
        .iter()
        .find(|report| report.limit_name == "T(t0)")
        .expect("expected dotted_bubble threshold report");

    assert!(
        !profile_result.all_passed,
        "expected dotted bubble bulk profile to remain failing:\n{profile_result}"
    );
    assert!(
        !threshold_report.passed,
        "expected dotted bubble threshold report to remain failing:\n{profile_result}"
    );

    let original_scaling = threshold_report
        .display_only_reports
        .iter()
        .find(|report| report.label == "original")
        .and_then(|report| report.scaling)
        .expect("expected original display-only profile scaling");
    let ct_scaling = threshold_report
        .display_only_reports
        .iter()
        .find(|report| report.label == "ct_0")
        .and_then(|report| report.scaling)
        .expect("expected ct_0 display-only profile scaling");

    assert!(
        (original_scaling + 2.0).abs() < 0.1,
        "expected original scaling near -2, got {original_scaling:+.4}:\n{profile_result}"
    );
    assert!(
        (ct_scaling + 2.0).abs() < 0.1,
        "expected ct_0 scaling near -2, got {ct_scaling:+.4}:\n{profile_result}"
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
