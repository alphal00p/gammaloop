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
fn dotted_bubble_amp_bulk_profile_passes() -> Result<()> {
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

    assert!(profile_result.all_passed, "{profile_result}");

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn double_dotted_bubble_amp_bulk_profile_passes() -> Result<()> {
    let state_path = get_tests_workspace_path().join("double_dotted_bubble_amp_bulk_profile");
    let run_history = run_card("double_dotted_bubble_amp_generate.toml")?;
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
        "expected tests/resources/run_cards/scalar_self_energy_amp_generate.toml to contain the profile bulk command"
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
fn scalar_self_energy_amp_bulk_profile_passes() -> Result<()> {
    let state_path = get_tests_workspace_path().join("scalar_self_energy_amp_bulk_profile");
    let run_history = run_card("scalar_self_energy_amp_generate.toml")?;
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
        "expected tests/resources/run_cards/scalar_self_energy_amp_generate.toml to contain the profile bulk command"
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
fn box_4e_bulk_profile_passes() -> Result<()> {
    let state_path = get_tests_workspace_path().join("box_4e_bulk_profile");
    let run_history = run_card("box_4e_generate.toml")?;
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
        "expected tests/resources/run_cards/box_4e_generate.toml to contain the profile bulk command"
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
