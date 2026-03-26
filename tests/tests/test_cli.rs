#![allow(dead_code)]
#![allow(unused_variables)]

use std::{
    fs,
    ops::{ControlFlow, Deref, DerefMut},
    path::PathBuf,
    sync::{Mutex, MutexGuard, OnceLock},
};

use color_eyre::Result;
use gammaloop_api::{
    CLISettings, OneShot, StateLoadOption,
    commands::Commands,
    state::{CommandHistory, CommandsBlock, RunHistory},
};
use gammaloop_integration_tests::{CLIState, clean_test, get_test_cli, get_tests_workspace_path};
use gammalooprs::{processes::ProcessCollection, settings::RuntimeSettings};
use serial_test::serial;

static TEMPLATE_CLI: OnceLock<Mutex<CLIState>> = OnceLock::new();

struct SharedCli<'a>(MutexGuard<'a, CLIState>);

impl Deref for SharedCli<'_> {
    type Target = CLIState;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for SharedCli<'_> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

fn cli_state_path(name: &str) -> PathBuf {
    get_tests_workspace_path().join("cli_refactor").join(name)
}

fn run_card_path(name: &str) -> PathBuf {
    get_tests_workspace_path()
        .join("cli_refactor")
        .join("run_cards")
        .join(name)
}

fn new_cli(name: &str) -> Result<SharedCli<'static>> {
    let template = TEMPLATE_CLI.get_or_init(|| {
        Mutex::new(
            get_test_cli(None, cli_state_path("_template"), None, true)
                .expect("template CLI should initialize"),
        )
    });
    let mut cli = SharedCli(template.lock().unwrap());
    let state_folder = cli_state_path(name);
    *cli = get_test_cli(None, &state_folder, None, true)?;
    Ok(cli)
}

fn command(raw: &str) -> CommandHistory {
    CommandHistory::from_raw_string(raw).unwrap()
}

fn block(name: &str, commands: &[&str]) -> CommandsBlock {
    CommandsBlock {
        name: name.to_string(),
        commands: commands.iter().map(|raw| command(raw)).collect(),
    }
}

fn history_strings(run_history: &RunHistory) -> Vec<String> {
    run_history
        .commands
        .iter()
        .map(|command| {
            command
                .raw_string
                .clone()
                .unwrap_or_else(|| format!("{:?}", command.command))
        })
        .collect()
}

fn duplicate_loaded_process(cli: &mut SharedCli<'_>, new_name: &str) {
    let mut cloned = cli.state.process_list.processes[0].clone();
    cloned.definition.folder_name = new_name.to_string();
    cloned.definition.process_id = cli.state.process_list.processes.len();
    cli.state.process_list.processes.push(cloned);
}

fn add_duplicate_integrand(
    cli: &mut SharedCli<'_>,
    process_index: usize,
    new_integrand_name: &str,
) {
    let process = &mut cli.state.process_list.processes[process_index];
    match &mut process.collection {
        ProcessCollection::Amplitudes(amplitudes) => {
            let source_name = amplitudes
                .keys()
                .next()
                .expect("fixture process should contain at least one amplitude")
                .clone();
            let mut duplicate = amplitudes
                .get(&source_name)
                .expect("fixture amplitude should exist")
                .clone();
            duplicate.name = new_integrand_name.to_string();
            amplitudes.insert(new_integrand_name.to_string(), duplicate);
        }
        ProcessCollection::CrossSections(cross_sections) => {
            let source_name = cross_sections
                .keys()
                .next()
                .expect("fixture process should contain at least one cross section")
                .clone();
            let mut duplicate = cross_sections
                .get(&source_name)
                .expect("fixture cross section should exist")
                .clone();
            duplicate.name = new_integrand_name.to_string();
            cross_sections.insert(new_integrand_name.to_string(), duplicate);
        }
    }
}

fn populate_generated_scalar_box_process(cli: &mut SharedCli<'_>) -> Result<()> {
    cli.run_command("import model scalars-default.json")?;
    cli.run_command(
        "generate amp scalar_0 scalar_0 > scalar_0 scalar_0 [{1}] --allowed-vertex-interactions V_3_SCALAR_022 -p box -i scalar_box --select-graphs GL0",
    )?;
    cli.run_command("generate")?;
    cli.save_state()?;
    Ok(())
}

fn run_without_arguments_is_a_noop() -> Result<()> {
    let mut cli = new_cli("run_without_arguments_is_a_noop")?;

    cli.run_command("run")?;

    assert!(cli.run_history.commands.is_empty());
    Ok(())
}

fn run_records_only_the_wrapper_command() -> Result<()> {
    let mut cli = new_cli("run_records_only_the_wrapper_command")?;
    cli.run_history.command_blocks = vec![
        block(
            "set_display",
            &["set global kv global.display_directive=warn"],
        ),
        block(
            "set_runtime",
            &["set default-runtime kv general.mu_r_sq=12.0"],
        ),
    ];

    cli.run_command(
        "run set_display set_runtime -c \"set global kv global.logfile_directive=error\"",
    )?;

    assert_eq!(cli.cli_settings.global.display_directive, "warn");
    assert_eq!(cli.default_runtime_settings.general.mu_r_sq, 12.0);
    assert_eq!(cli.cli_settings.global.logfile_directive, "error");
    assert_eq!(cli.run_history.commands.len(), 1);
    assert!(matches!(
        cli.run_history.commands[0].command,
        Commands::Run(_)
    ));
    Ok(())
}

fn run_prevalidation_is_all_or_nothing_for_inline_commands() -> Result<()> {
    let mut cli = new_cli("run_prevalidation_is_all_or_nothing_for_inline_commands")?;
    let original_display = cli.cli_settings.global.display_directive.clone();

    let err = cli
        .run_command("run -c \"set global kv global.display_directive=warn; no_such_command\"")
        .unwrap_err();
    let error_text = format!("{err:?}");

    assert!(error_text.contains("Failed to parse run -c command #2"));
    assert_eq!(cli.cli_settings.global.display_directive, original_display);
    assert!(cli.run_history.commands.is_empty());
    Ok(())
}

fn run_prevalidation_is_all_or_nothing_for_nested_block_failures() -> Result<()> {
    let mut cli = new_cli("run_prevalidation_is_all_or_nothing_for_nested_block_failures")?;
    let original_display = cli.cli_settings.global.display_directive.clone();
    cli.run_history.command_blocks = vec![
        block("first", &["set global kv global.display_directive=warn"]),
        block("broken", &["run missing_block"]),
    ];

    let err = cli.run_command("run first broken").unwrap_err();
    let error_text = format!("{err:?}");

    assert!(error_text.contains("command block 'broken' command #1"));
    assert!(error_text.contains("Unknown command block 'missing_block'"));
    assert_eq!(cli.cli_settings.global.display_directive, original_display);
    assert!(cli.run_history.commands.is_empty());
    Ok(())
}

fn nested_run_records_only_the_top_level_run() -> Result<()> {
    let mut cli = new_cli("nested_run_records_only_the_top_level_run")?;
    cli.run_history.command_blocks = vec![
        block("inner", &["set global kv global.display_directive=warn"]),
        block("outer", &["run inner"]),
    ];

    cli.run_command("run outer")?;

    assert_eq!(cli.cli_settings.global.display_directive, "warn");
    assert_eq!(history_strings(&cli.run_history), vec!["run outer"]);
    Ok(())
}

#[allow(clippy::needless_update)]
fn boot_run_history_merges_blocks_and_persists_commands_once() -> Result<()> {
    let mut cli = new_cli("boot_run_history_merges_blocks_and_persists_commands_once")?;
    let mut frozen_global = cli.cli_settings.global.clone();
    frozen_global.display_directive = "warn".into();
    let mut frozen_runtime = RuntimeSettings::default();
    frozen_runtime.general.mu_r_sq = 24.0;

    let boot_run_history = RunHistory {
        cli_settings: CLISettings {
            global: frozen_global.clone(),
            ..CLISettings::default()
        },
        default_runtime_settings: frozen_runtime.clone(),
        command_blocks: vec![block(
            "boot_block",
            &["set global kv global.display_directive=warn"],
        )],
        commands: vec![
            command("run boot_block"),
            command("set default-runtime kv general.mu_r_sq=24.0"),
        ],
        ..RunHistory::default()
    };

    let flow = cli.apply_boot_run_history(&boot_run_history)?;

    assert!(matches!(flow, ControlFlow::Continue(())));
    assert_eq!(cli.cli_settings.global.display_directive, "warn");
    assert_eq!(cli.default_runtime_settings.general.mu_r_sq, 24.0);
    assert_eq!(cli.run_history.command_blocks.len(), 1);
    assert_eq!(history_strings(&cli.run_history).len(), 2);

    cli.save_state()?;

    let persisted = RunHistory::load(cli.cli_settings.state.folder.join("run.toml"))?;
    assert_eq!(persisted.command_blocks.len(), 1);
    assert_eq!(persisted.commands.len(), 2);
    assert_eq!(persisted.cli_settings.global, frozen_global);
    assert_eq!(persisted.default_runtime_settings, frozen_runtime);
    Ok(())
}

fn boot_run_history_rejects_conflicting_block_redefinitions() -> Result<()> {
    let mut cli = new_cli("boot_run_history_rejects_conflicting_block_redefinitions")?;
    cli.run_history.command_blocks = vec![block(
        "shared",
        &["set global kv global.display_directive=warn"],
    )];
    let boot_run_history = RunHistory {
        command_blocks: vec![block(
            "shared",
            &["set global kv global.display_directive=error"],
        )],
        ..RunHistory::default()
    };

    let err = cli.apply_boot_run_history(&boot_run_history).unwrap_err();
    let error_text = format!("{err:?}");

    assert!(error_text.contains("redefines an existing block with different commands"));
    assert_eq!(cli.run_history.command_blocks.len(), 1);
    assert!(cli.run_history.commands.is_empty());
    Ok(())
}

fn boot_run_history_allows_conflicting_redefinitions_after_confirmation() -> Result<()> {
    let mut cli = new_cli("boot_run_history_allows_conflicting_redefinitions_after_confirmation")?;
    cli.run_history.command_blocks = vec![block(
        "shared",
        &["set global kv global.display_directive=warn"],
    )];
    let boot_run_history = RunHistory {
        command_blocks: vec![block(
            "shared",
            &["set global kv global.display_directive=error"],
        )],
        commands: vec![command("run shared")],
        ..RunHistory::default()
    };

    let flow = cli.apply_boot_run_history_with_conflict_resolution(&boot_run_history, true)?;

    assert!(matches!(flow, ControlFlow::Continue(())));
    assert_eq!(cli.run_history.command_blocks.len(), 1);
    assert_eq!(
        cli.run_history.command_blocks[0].commands[0]
            .raw_string
            .as_deref(),
        Some("set global kv global.display_directive=error")
    );
    assert_eq!(cli.cli_settings.global.display_directive, "error");
    Ok(())
}

fn boot_run_history_cancelled_conflicting_redefinitions_leave_state_untouched() -> Result<()> {
    let mut cli =
        new_cli("boot_run_history_cancelled_conflicting_redefinitions_leave_state_untouched")?;
    cli.run_history.command_blocks = vec![block(
        "shared",
        &["set global kv global.display_directive=warn"],
    )];
    let boot_run_history = RunHistory {
        command_blocks: vec![block(
            "shared",
            &["set global kv global.display_directive=error"],
        )],
        commands: vec![command("run shared")],
        ..RunHistory::default()
    };

    let flow = cli.apply_boot_run_history_with_conflict_resolution(&boot_run_history, false)?;

    assert!(matches!(flow, ControlFlow::Break(_)));
    assert_eq!(cli.run_history.command_blocks.len(), 1);
    assert_eq!(
        cli.run_history.command_blocks[0].commands[0]
            .raw_string
            .as_deref(),
        Some("set global kv global.display_directive=warn")
    );
    assert!(cli.run_history.commands.is_empty());
    Ok(())
}

fn boot_run_history_conflicting_redefinitions_error_in_read_only_mode() -> Result<()> {
    let mut cli = new_cli("boot_run_history_conflicting_redefinitions_error_in_read_only_mode")?;
    cli.cli_settings.session.read_only_state = true;
    cli.run_history.command_blocks = vec![block(
        "shared",
        &["set global kv global.display_directive=warn"],
    )];
    let boot_run_history = RunHistory {
        command_blocks: vec![block(
            "shared",
            &["set global kv global.display_directive=error"],
        )],
        ..RunHistory::default()
    };

    let err = cli
        .apply_boot_run_history_with_conflict_resolution(&boot_run_history, true)
        .unwrap_err();
    let error_text = format!("{err:?}");

    assert!(error_text.contains("shared"));
    assert!(error_text.contains("--read-only-state"));
    Ok(())
}

fn boot_run_history_allows_identical_semantic_redefinitions() -> Result<()> {
    let mut cli = new_cli("boot_run_history_allows_identical_semantic_redefinitions")?;
    cli.run_history.command_blocks = vec![block(
        "shared",
        &["set global kv global.display_directive=warn"],
    )];
    let boot_run_history = RunHistory {
        command_blocks: vec![CommandsBlock {
            name: "shared".to_string(),
            commands: vec![CommandHistory::new_with_raw(
                command("set global kv global.display_directive=warn").command,
                "set global kv   global.display_directive = warn".to_string(),
            )],
        }],
        ..RunHistory::default()
    };

    let flow = cli.apply_boot_run_history(&boot_run_history)?;

    assert!(matches!(flow, ControlFlow::Continue(())));
    assert_eq!(cli.run_history.command_blocks.len(), 1);
    assert!(cli.run_history.commands.is_empty());
    Ok(())
}

fn booting_existing_state_with_mismatched_frozen_settings_forces_read_only_and_warns() -> Result<()>
{
    let test_name =
        "booting_existing_state_with_mismatched_frozen_settings_forces_read_only_and_warns";
    let mut cli = new_cli(test_name)?;

    cli.run_history.cli_settings.global.display_directive = "info".into();
    cli.run_history.default_runtime_settings.general.mu_r_sq = 11.0;
    cli.save_state()?;

    let boot_card_path = run_card_path("boot_settings_mismatch.toml");
    if let Some(parent) = boot_card_path.parent() {
        fs::create_dir_all(parent)?;
    }

    let mut boot_global = CLISettings::default().global;
    boot_global.display_directive = "warn".into();
    let mut boot_runtime = RuntimeSettings::default();
    boot_runtime.general.mu_r_sq = 29.0;
    let boot_run_history = RunHistory {
        cli_settings: CLISettings {
            global: boot_global.clone(),
            ..CLISettings::default()
        },
        default_runtime_settings: boot_runtime.clone(),
        ..RunHistory::default()
    };
    fs::write(&boot_card_path, toml::to_string_pretty(&boot_run_history)?)?;

    let loaded = StateLoadOption {
        state_folder: Some(cli.cli_settings.state.folder.clone()),
        boot_commands_path: Some(boot_card_path),
        ..StateLoadOption::default()
    }
    .load()?;

    assert!(loaded.cli_settings.session.read_only_state);
    assert!(
        loaded
            .cli_settings
            .session
            .startup_warnings
            .iter()
            .any(|warning| warning.contains("forced into --read-only-state"))
    );
    assert!(
        loaded
            .cli_settings
            .session
            .startup_warnings
            .iter()
            .any(|warning| warning
                .contains("Line up cli_settings.global and default_runtime_settings"))
    );
    Ok(())
}

#[test]
#[serial]
fn state_load_option_clean_state_removes_saved_state_before_load() -> Result<()> {
    let test_name = "state_load_option_clean_state_removes_saved_state_before_load";
    let mut cli = new_cli(test_name)?;

    cli.run_history.commands.push(command("quit -n"));
    cli.save_state()?;

    let state_folder = cli.cli_settings.state.folder.clone();
    assert!(state_folder.exists());

    let loaded = StateLoadOption {
        state_folder: Some(state_folder.clone()),
        clean_state: true,
        ..StateLoadOption::default()
    }
    .load()?;

    assert!(loaded.run_history.commands.is_empty());
    assert!(loaded.state_load_summary.is_none());
    assert_eq!(loaded.cli_settings.state.folder, state_folder);
    assert!(!loaded.cli_settings.state.folder.exists());
    Ok(())
}

fn command_block_definition_mode_defers_execution_and_omits_history() -> Result<()> {
    let mut cli = new_cli("command_block_definition_mode_defers_execution_and_omits_history")?;
    let original_display = cli.cli_settings.global.display_directive.clone();

    cli.run_command("start_commands_block demo")?;
    cli.run_command("set global kv global.display_directive=warn")?;

    assert_eq!(cli.cli_settings.global.display_directive, original_display);
    assert!(cli.run_history.commands.is_empty());

    cli.run_command("finish_commands_block")?;

    assert_eq!(cli.run_history.command_blocks.len(), 1);
    assert_eq!(cli.run_history.command_blocks[0].name, "demo");
    assert!(cli.run_history.commands.is_empty());

    cli.run_command("run demo")?;

    assert_eq!(cli.cli_settings.global.display_directive, "warn");
    assert_eq!(history_strings(&cli.run_history), vec!["run demo"]);
    Ok(())
}

fn finish_commands_block_without_active_definition_fails() -> Result<()> {
    let mut cli = new_cli("finish_commands_block_without_active_definition_fails")?;

    let err = cli.run_command("finish_commands_block").unwrap_err();
    let error_text = format!("{err:?}");

    assert!(error_text.contains("No command block definition is currently being recorded"));
    Ok(())
}

fn start_commands_block_cannot_nest() -> Result<()> {
    let mut cli = new_cli("start_commands_block_cannot_nest")?;

    cli.run_command("start_commands_block first")?;
    let err = cli.run_command("start_commands_block second").unwrap_err();
    let error_text = format!("{err:?}");

    assert!(error_text.contains("Cannot start a new command block definition"));

    cli.run_command("finish_commands_block")?;

    assert_eq!(cli.run_history.command_blocks.len(), 1);
    assert_eq!(cli.run_history.command_blocks[0].name, "first");
    Ok(())
}

fn quit_during_block_definition_dismisses_the_pending_block() -> Result<()> {
    let mut cli = new_cli("quit_during_block_definition_dismisses_the_pending_block")?;

    cli.run_command("start_commands_block demo")?;
    cli.run_command("set global kv global.display_directive=warn")?;

    let flow = cli.run_command_flow("quit -o")?.flow;

    assert!(matches!(flow, ControlFlow::Continue(())));
    assert!(cli.run_history.command_blocks.is_empty());
    assert!(cli.run_history.commands.is_empty());

    let err = cli.run_command("finish_commands_block").unwrap_err();
    assert!(format!("{err:?}").contains("No command block definition is currently being recorded"));
    Ok(())
}

fn ctrl_c_dismisses_pending_block_definition() -> Result<()> {
    let mut cli = new_cli("ctrl_c_dismisses_pending_block_definition")?;

    cli.run_command("start_commands_block demo")?;
    cli.run_command("set global kv global.display_directive=warn")?;

    assert!(cli.dismiss_pending_commands_block("Ctrl-C"));
    assert!(cli.run_history.command_blocks.is_empty());
    assert!(cli.run_history.commands.is_empty());
    Ok(())
}

fn ctrl_d_dismisses_pending_block_definition() -> Result<()> {
    let mut cli = new_cli("ctrl_d_dismisses_pending_block_definition")?;

    cli.run_command("start_commands_block demo")?;
    cli.run_command("set global kv global.display_directive=warn")?;

    assert!(cli.dismiss_pending_commands_block("Ctrl-D"));
    assert!(cli.run_history.command_blocks.is_empty());
    assert!(cli.run_history.commands.is_empty());
    Ok(())
}

fn recursive_run_limit_is_enforced() -> Result<()> {
    let mut cli = new_cli("recursive_run_limit_is_enforced")?;
    cli.run_history.command_blocks = vec![block("loop", &["run loop"])];

    let err = cli.run_command("run loop").unwrap_err();
    let error_text = format!("{err:?}");

    assert!(error_text.contains("Maximum nested run depth of 100 reached"));
    assert!(cli.run_history.commands.is_empty());
    Ok(())
}

fn save_state_writes_global_settings_file() -> Result<()> {
    let mut cli = new_cli("save_state_writes_global_settings_file")?;

    cli.save_state()?;

    let global_settings_path = cli.cli_settings.state.folder.join("global_settings.toml");
    let global_settings_contents = fs::read_to_string(&global_settings_path)?;

    assert!(global_settings_path.exists());
    assert!(!global_settings_contents.contains("dummy"));
    Ok(())
}

fn reset_processes_with_process_selector_removes_only_that_process() -> Result<()> {
    let mut cli = new_cli("reset_processes_with_process_selector_removes_only_that_process")?;
    populate_generated_scalar_box_process(&mut cli)?;
    duplicate_loaded_process(&mut cli, "second_process");

    cli.run_command("reset processes -p second_process")?;

    assert_eq!(cli.state.process_list.processes.len(), 1);
    assert_eq!(
        cli.state.process_list.processes[0].definition.folder_name,
        "box"
    );
    Ok(())
}

fn reset_processes_with_integrand_selector_removes_only_that_integrand() -> Result<()> {
    let mut cli = new_cli("reset_processes_with_integrand_selector_removes_only_that_integrand")?;
    populate_generated_scalar_box_process(&mut cli)?;
    add_duplicate_integrand(&mut cli, 0, "virtual_copy");

    cli.run_command("reset processes -p box -i virtual_copy")?;

    let remaining = cli.state.process_list.processes[0]
        .collection
        .get_integrand_names();
    assert_eq!(remaining, vec!["scalar_box"]);
    Ok(())
}

fn reset_processes_removing_last_integrand_drops_empty_process() -> Result<()> {
    let mut cli = new_cli("reset_processes_removing_last_integrand_drops_empty_process")?;
    populate_generated_scalar_box_process(&mut cli)?;

    cli.run_command("reset processes -p box -i scalar_box")?;

    assert!(cli.state.process_list.processes.is_empty());
    Ok(())
}

fn reset_processes_rejects_integrand_without_process() -> Result<()> {
    let mut cli = new_cli("reset_processes_rejects_integrand_without_process")?;
    let err = cli.run_command("reset processes -i 1L").unwrap_err();
    assert!(format!("{err:?}").contains("--integrand-name requires --process"));
    Ok(())
}

#[test]
#[serial]
fn cli_stateful_workflow_behaviors() -> Result<()> {
    run_without_arguments_is_a_noop()?;
    run_records_only_the_wrapper_command()?;
    run_prevalidation_is_all_or_nothing_for_inline_commands()?;
    run_prevalidation_is_all_or_nothing_for_nested_block_failures()?;
    nested_run_records_only_the_top_level_run()?;
    boot_run_history_merges_blocks_and_persists_commands_once()?;
    boot_run_history_rejects_conflicting_block_redefinitions()?;
    boot_run_history_allows_conflicting_redefinitions_after_confirmation()?;
    boot_run_history_cancelled_conflicting_redefinitions_leave_state_untouched()?;
    boot_run_history_conflicting_redefinitions_error_in_read_only_mode()?;
    boot_run_history_allows_identical_semantic_redefinitions()?;
    booting_existing_state_with_mismatched_frozen_settings_forces_read_only_and_warns()?;
    command_block_definition_mode_defers_execution_and_omits_history()?;
    finish_commands_block_without_active_definition_fails()?;
    start_commands_block_cannot_nest()?;
    quit_during_block_definition_dismisses_the_pending_block()?;
    ctrl_c_dismisses_pending_block_definition()?;
    ctrl_d_dismisses_pending_block_definition()?;
    recursive_run_limit_is_enforced()?;
    reset_processes_with_process_selector_removes_only_that_process()?;
    reset_processes_with_integrand_selector_removes_only_that_integrand()?;
    reset_processes_removing_last_integrand_drops_empty_process()?;
    reset_processes_rejects_integrand_without_process()?;
    save_state_writes_global_settings_file()?;
    Ok(())
}

#[test]
#[serial]
fn run_history_load_reports_detailed_parse_errors() {
    let run_card_path = run_card_path("run_history_load_reports_detailed_parse_errors.toml");
    if let Some(parent) = run_card_path.parent() {
        fs::create_dir_all(parent).unwrap();
    }
    fs::write(&run_card_path, "commands = [\"no_such_command\"]\n").unwrap();

    let err = RunHistory::load(&run_card_path).unwrap_err();
    let error_text = format!("{err:?}");

    assert!(error_text.contains("Failed to parse top-level command #1 'no_such_command'"));
    assert!(error_text.contains("error:"));
}

#[test]
#[serial]
fn run_history_load_reports_block_command_context() {
    let run_card_path = run_card_path("run_history_load_reports_block_command_context.toml");
    if let Some(parent) = run_card_path.parent() {
        fs::create_dir_all(parent).unwrap();
    }
    fs::write(
        &run_card_path,
        r#"
[[command_blocks]]
name = "demo"
commands = ["no_such_command"]
"#,
    )
    .unwrap();

    let err = RunHistory::load(&run_card_path).unwrap_err();
    let error_text = format!("{err:?}");

    assert!(error_text.contains("command block 'demo' command #1"));
    assert!(error_text.contains("no_such_command"));
}

#[test]
#[serial]
fn existing_invalid_state_folder_fails_to_load_instead_of_falling_back() -> Result<()> {
    let state_path =
        cli_state_path("existing_invalid_state_folder_fails_to_load_instead_of_falling_back");
    clean_test(&state_path);
    fs::create_dir_all(&state_path)?;
    fs::write(state_path.join("state_manifest.toml"), "version = 999\n")?;

    let mut cli = OneShot::new_test(state_path.clone());
    let err = match cli.load() {
        Ok(_) => panic!("expected load to fail"),
        Err(err) => err,
    };
    let error_text = format!("{err:?}");

    assert!(error_text.contains("State version 999"));

    clean_test(&state_path);
    Ok(())
}
