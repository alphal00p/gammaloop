#![allow(dead_code)]
#![allow(unused_variables)]

use std::{
    collections::BTreeSet,
    env, fs,
    ops::{ControlFlow, Deref, DerefMut},
    path::{Path, PathBuf},
    sync::{Mutex, MutexGuard, OnceLock},
};

use color_eyre::Result;
use gammaloop_api::{
    CLISettings, OneShot, StateLoadOption,
    state::{CommandHistory, CommandsBlock, RunHistory},
};
use gammaloop_integration_tests::{CLIState, clean_test, get_test_cli, get_tests_workspace_path};
use gammalooprs::{processes::ProcessCollection, settings::RuntimeSettings};
use schemars::{JsonSchema, schema_for};
use serde_json::Value as JsonValue;
use serial_test::serial;
use toml::Value as TomlValue;

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

fn graph_fixture_path(name: &str) -> PathBuf {
    get_tests_workspace_path()
        .parent()
        .expect("tests workspace should have a parent")
        .join("resources")
        .join("graphs")
        .join(name)
}

struct CurrentDirGuard {
    previous: PathBuf,
}

impl CurrentDirGuard {
    fn enter(path: &Path) -> Result<Self> {
        let previous = env::current_dir()?;
        fs::create_dir_all(path)?;
        env::set_current_dir(path)?;
        Ok(Self { previous })
    }
}

impl Drop for CurrentDirGuard {
    fn drop(&mut self) {
        env::set_current_dir(&self.previous)
            .expect("current working directory should be restored after the test");
    }
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

fn run_records_executed_commands() -> Result<()> {
    let mut cli = new_cli("run_records_executed_commands")?;
    cli.run_history.command_blocks = vec![
        block(
            "set_display",
            &["set global kv global.display_directive=warn"],
        ),
        block("set_runtime", &["set default-runtime kv general.mu_r=12.0"]),
    ];

    cli.run_command(
        "run set_display set_runtime -c \"set global kv global.logfile_directive=error\"",
    )?;

    assert_eq!(cli.cli_settings.global.display_directive, "warn");
    assert_eq!(cli.default_runtime_settings.general.mu_r, 12.0);
    assert_eq!(cli.cli_settings.global.logfile_directive, "error");
    assert_eq!(
        history_strings(&cli.run_history),
        vec!["run set_display set_runtime -c \"set global kv global.logfile_directive=error\""]
    );
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

fn nested_run_records_executed_commands_once() -> Result<()> {
    let mut cli = new_cli("nested_run_records_executed_commands_once")?;
    cli.run_history.command_blocks = vec![
        block("inner", &["set global kv global.display_directive=warn"]),
        block("outer", &["run inner"]),
    ];

    cli.run_command("run outer")?;

    assert_eq!(cli.cli_settings.global.display_directive, "warn");
    assert_eq!(history_strings(&cli.run_history), vec!["run outer"]);
    Ok(())
}

fn run_command_block_substitutes_define_variables() -> Result<()> {
    let mut cli = new_cli("run_command_block_substitutes_define_variables")?;
    cli.run_history.command_blocks = vec![block(
        "set_display",
        &["set global kv global.display_directive={level}"],
    )];

    cli.run_command("run set_display -D level=warn")?;

    assert_eq!(cli.cli_settings.global.display_directive, "warn");
    assert_eq!(
        history_strings(&cli.run_history),
        vec!["run set_display -D level=warn"]
    );
    cli.save_state()?;
    let persisted = fs::read_to_string(cli.cli_settings.state.folder.join("run.toml"))?;
    assert!(persisted.contains("run set_display -D level=warn"));
    Ok(())
}

fn run_command_blocks_share_define_environment_and_allow_unused_keys() -> Result<()> {
    let mut cli = new_cli("run_command_blocks_share_define_environment_and_allow_unused_keys")?;
    cli.run_history.command_blocks = vec![
        block(
            "set_display",
            &["set global kv global.display_directive={display}"],
        ),
        block(
            "set_logfile",
            &["set global kv global.logfile_directive={logfile}"],
        ),
    ];

    cli.run_command("run set_display set_logfile -D display=warn -D logfile=error -D unused=ok")?;

    assert_eq!(cli.cli_settings.global.display_directive, "warn");
    assert_eq!(cli.cli_settings.global.logfile_directive, "error");
    assert_eq!(
        history_strings(&cli.run_history),
        vec!["run set_display set_logfile -D display=warn -D logfile=error -D unused=ok"]
    );
    Ok(())
}

fn nested_run_inherits_and_can_override_define_environment() -> Result<()> {
    let mut cli = new_cli("nested_run_inherits_and_can_override_define_environment")?;
    cli.run_history.command_blocks = vec![
        block("inner", &["set global kv global.display_directive={level}"]),
        block("outer_inherit", &["run inner"]),
        block("outer_override", &["run inner -D level=error"]),
    ];

    cli.run_command("run outer_inherit -D level=warn")?;
    assert_eq!(cli.cli_settings.global.display_directive, "warn");

    cli.run_command("run outer_override -D level=warn")?;
    assert_eq!(cli.cli_settings.global.display_directive, "error");
    Ok(())
}

fn run_define_placeholders_are_prevalidated_before_execution() -> Result<()> {
    let mut cli = new_cli("run_define_placeholders_are_prevalidated_before_execution")?;
    let original_display = cli.cli_settings.global.display_directive.clone();
    cli.run_history.command_blocks = vec![
        block("first", &["set global kv global.display_directive=warn"]),
        block(
            "broken",
            &["set global kv global.display_directive={level}"],
        ),
    ];

    let err = cli.run_command("run first broken").unwrap_err();
    let error_text = format!("{err:?}");

    assert!(error_text.contains("Missing command-block variable 'level'"));
    assert_eq!(cli.cli_settings.global.display_directive, original_display);
    assert!(cli.run_history.commands.is_empty());
    Ok(())
}

fn run_inline_commands_can_use_define_variables() -> Result<()> {
    let mut cli = new_cli("run_inline_commands_can_use_define_variables")?;

    cli.run_command("run -D level=warn -c 'set global kv global.display_directive={level}'")?;

    assert_eq!(cli.cli_settings.global.display_directive, "warn");
    assert_eq!(
        history_strings(&cli.run_history),
        vec!["run -D level=warn -c 'set global kv global.display_directive={level}'"]
    );
    Ok(())
}

#[allow(clippy::needless_update)]
fn boot_run_history_merges_blocks_and_persists_commands_once() -> Result<()> {
    let mut cli = new_cli("boot_run_history_merges_blocks_and_persists_commands_once")?;
    let mut frozen_global = cli.cli_settings.global.clone();
    frozen_global.display_directive = "warn".into();
    let mut frozen_runtime = RuntimeSettings::default();
    frozen_runtime.general.mu_r = 24.0;

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
            command("set default-runtime kv general.mu_r=24.0"),
        ],
        ..RunHistory::default()
    };

    let flow = cli.apply_boot_run_history(&boot_run_history)?;

    assert!(matches!(flow, ControlFlow::Continue(())));
    assert_eq!(cli.cli_settings.global.display_directive, "warn");
    assert_eq!(cli.default_runtime_settings.general.mu_r, 24.0);
    assert_eq!(cli.run_history.command_blocks.len(), 1);
    assert_eq!(
        history_strings(&cli.run_history),
        vec!["run boot_block", "set default-runtime kv general.mu_r=24.0",]
    );

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
    cli.run_history.default_runtime_settings.general.mu_r = 11.0;
    cli.save_state()?;

    let boot_card_path = run_card_path("boot_settings_mismatch.toml");
    if let Some(parent) = boot_card_path.parent() {
        fs::create_dir_all(parent)?;
    }

    let mut boot_global = CLISettings::default().global;
    boot_global.display_directive = "warn".into();
    let mut boot_runtime = RuntimeSettings::default();
    boot_runtime.general.mu_r = 29.0;
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

fn boot_mismatch_uses_conflicting_boot_blocks_read_only() -> Result<()> {
    let test_name = "boot_mismatch_uses_conflicting_boot_blocks_read_only";
    let mut cli = new_cli(test_name)?;

    cli.run_history.default_runtime_settings.general.mu_r = 11.0;
    cli.run_history.command_blocks = vec![block(
        "shared",
        &["set global kv global.display_directive=info"],
    )];
    cli.save_state()?;

    let boot_card_path = run_card_path("boot_settings_mismatch_conflicting_block.toml");
    if let Some(parent) = boot_card_path.parent() {
        fs::create_dir_all(parent)?;
    }

    let mut boot_runtime = RuntimeSettings::default();
    boot_runtime.general.mu_r = 29.0;
    let boot_run_history = RunHistory {
        default_runtime_settings: boot_runtime,
        command_blocks: vec![block(
            "shared",
            &["set global kv global.display_directive=warn"],
        )],
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
    assert_eq!(loaded.run_history.command_blocks.len(), 1);
    assert!(loaded.run_history.command_blocks[0].semantically_eq(&block(
        "shared",
        &["set global kv global.display_directive=warn"],
    )));
    assert!(
        loaded
            .cli_settings
            .session
            .startup_warnings
            .iter()
            .any(|warning| warning.contains("will use the boot card command block definitions"))
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

fn run_command_block_then_inline_quit_preserves_override_flag() -> Result<()> {
    let mut cli = new_cli("run_command_block_then_inline_quit_preserves_override_flag")?;
    cli.run_history.command_blocks.push(block(
        "demo",
        &["set global kv global.display_directive=warn"],
    ));

    let flow = cli.run_command_flow("run demo -c 'quit -o'")?.flow;
    let ControlFlow::Break(save_state) = flow else {
        panic!("expected run demo -c 'quit -o' to break the session");
    };

    assert_eq!(save_state.override_state, Some(true));
    assert_eq!(history_strings(&cli.run_history), vec!["run demo"]);
    assert_eq!(cli.cli_settings.global.display_directive, "warn");
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

#[test]
#[serial]
fn save_state_writes_exhaustive_default_settings_files() -> Result<()> {
    let mut cli = new_cli("save_state_writes_exhaustive_default_settings_files")?;

    cli.save_state()?;

    assert_settings_file_contains_schema_paths::<CLISettings>(
        &cli.cli_settings.state.folder.join("global_settings.toml"),
    )?;
    assert_settings_file_contains_schema_paths::<RuntimeSettings>(
        &cli.cli_settings
            .state
            .folder
            .join("default_runtime_settings.toml"),
    )?;

    Ok(())
}

fn assert_settings_file_contains_schema_paths<T: JsonSchema>(path: &Path) -> Result<()> {
    let contents = fs::read_to_string(path)?;
    let value: TomlValue = toml::from_str(&contents)?;
    let actual_paths = toml_paths(&value);
    let schema = serde_json::to_value(schema_for!(T))?;
    let mut expected_paths = BTreeSet::new();
    collect_schema_visible_paths(
        &schema,
        &schema,
        &mut Vec::new(),
        &mut expected_paths,
        &mut Vec::new(),
    );

    let missing = expected_paths
        .difference(&actual_paths)
        .cloned()
        .collect::<Vec<_>>();
    let missing_preview = missing
        .iter()
        .take(80)
        .map(String::as_str)
        .collect::<Vec<_>>()
        .join("\n");
    let omitted_count = missing.len().saturating_sub(80);
    let omitted_message = if omitted_count == 0 {
        String::new()
    } else {
        format!("\n... and {omitted_count} more")
    };

    assert!(
        missing.is_empty(),
        "{} is missing schema-backed default setting paths:\n{}",
        path.display(),
        format!("{missing_preview}{omitted_message}")
    );

    Ok(())
}

fn toml_paths(value: &TomlValue) -> BTreeSet<String> {
    let mut paths = BTreeSet::new();
    collect_toml_paths(value, &mut Vec::new(), &mut paths);
    paths
}

fn collect_toml_paths(value: &TomlValue, path: &mut Vec<String>, paths: &mut BTreeSet<String>) {
    if !path.is_empty() {
        paths.insert(path.join("."));
    }

    match value {
        TomlValue::Table(table) => {
            for (key, value) in table {
                path.push(key.clone());
                collect_toml_paths(value, path, paths);
                path.pop();
            }
        }
        TomlValue::Array(items) => {
            for item in items {
                if matches!(item, TomlValue::Table(_) | TomlValue::Array(_)) {
                    collect_toml_paths(item, path, paths);
                }
            }
        }
        _ => {}
    }
}

fn collect_schema_visible_paths(
    schema_root: &JsonValue,
    node: &JsonValue,
    path: &mut Vec<String>,
    paths: &mut BTreeSet<String>,
    ref_stack: &mut Vec<String>,
) {
    if let Some(reference) = node.get("$ref").and_then(JsonValue::as_str) {
        if ref_stack.iter().any(|existing| existing == reference) {
            if !path.is_empty() {
                paths.insert(path.join("."));
            }
            return;
        }
        let Some(resolved) = schema_ref_target(schema_root, reference) else {
            return;
        };
        ref_stack.push(reference.to_string());
        collect_schema_visible_paths(schema_root, resolved, path, paths, ref_stack);
        ref_stack.pop();
        return;
    }

    if schema_is_nullable(schema_root, node, &mut Vec::new()) {
        // TOML has no native null value, so default `None` fields cannot be
        // exhaustively represented without a field-specific sentinel encoding.
        return;
    }

    if let Some(properties) = node.get("properties").and_then(JsonValue::as_object) {
        for (key, child) in properties {
            path.push(key.clone());
            collect_schema_visible_paths(schema_root, child, path, paths, ref_stack);
            path.pop();
        }
        return;
    }

    if schema_has_keyword_variants(node) {
        for keyword in ["allOf", "anyOf", "oneOf"] {
            let Some(variants) = node.get(keyword).and_then(JsonValue::as_array) else {
                continue;
            };
            for variant in variants {
                collect_schema_visible_paths(schema_root, variant, path, paths, ref_stack);
            }
        }
        return;
    }

    if !path.is_empty() {
        paths.insert(path.join("."));
    }
}

fn schema_ref_target<'a>(schema_root: &'a JsonValue, reference: &str) -> Option<&'a JsonValue> {
    schema_root.pointer(reference.strip_prefix('#')?)
}

fn schema_is_nullable(
    schema_root: &JsonValue,
    node: &JsonValue,
    ref_stack: &mut Vec<String>,
) -> bool {
    if let Some(reference) = node.get("$ref").and_then(JsonValue::as_str) {
        if ref_stack.iter().any(|existing| existing == reference) {
            return false;
        }
        let Some(resolved) = schema_ref_target(schema_root, reference) else {
            return false;
        };
        ref_stack.push(reference.to_string());
        let nullable = schema_is_nullable(schema_root, resolved, ref_stack);
        ref_stack.pop();
        return nullable;
    }

    if schema_type_names(node).any(|type_name| type_name == "null") {
        return true;
    }

    ["allOf", "anyOf", "oneOf"].iter().any(|keyword| {
        node.get(keyword)
            .and_then(JsonValue::as_array)
            .is_some_and(|variants| {
                variants
                    .iter()
                    .any(|variant| schema_is_nullable(schema_root, variant, ref_stack))
            })
    })
}

fn schema_type_names(node: &JsonValue) -> impl Iterator<Item = &str> {
    node.get("type")
        .into_iter()
        .flat_map(|type_value| match type_value {
            JsonValue::String(type_name) => vec![type_name.as_str()],
            JsonValue::Array(type_names) => {
                type_names.iter().filter_map(JsonValue::as_str).collect()
            }
            _ => Vec::new(),
        })
}

fn schema_has_keyword_variants(node: &JsonValue) -> bool {
    ["allOf", "anyOf", "oneOf"].iter().any(|keyword| {
        node.get(keyword)
            .and_then(JsonValue::as_array)
            .is_some_and(|variants| !variants.is_empty())
    })
}

fn import_graphs_relative_path_prefers_active_state_root_parent() -> Result<()> {
    let root = cli_state_path("import_graphs_relative_path_prefers_active_state_root_parent");
    clean_test(&root);
    fs::create_dir_all(&root)?;
    let graph_name = "relative_state_root.dot";
    fs::copy(graph_fixture_path("scalar_box.dot"), root.join(graph_name))?;

    let mut cli = get_test_cli(None, root.join("state"), None, true)?;
    cli.run_command("import model scalars-default.json")?;
    cli.run_command(&format!("import graphs ./{graph_name}"))?;

    assert_eq!(cli.state.process_list.processes.len(), 1);
    assert_eq!(
        cli.state.process_list.processes[0].definition.folder_name,
        "relative_state_root"
    );
    Ok(())
}

fn import_graphs_relative_path_falls_back_to_current_working_directory() -> Result<()> {
    let root =
        cli_state_path("import_graphs_relative_path_falls_back_to_current_working_directory");
    clean_test(&root);
    fs::create_dir_all(&root)?;
    let cwd = root.join("cwd");
    let graph_name = "relative_cwd.dot";
    fs::create_dir_all(&cwd)?;
    fs::copy(graph_fixture_path("scalar_box.dot"), cwd.join(graph_name))?;
    let _cwd_guard = CurrentDirGuard::enter(&cwd)?;

    let mut cli = get_test_cli(None, root.join("state"), None, true)?;
    cli.run_command("import model scalars-default.json")?;
    cli.run_command(&format!("import graphs ./{graph_name}"))?;

    assert_eq!(cli.state.process_list.processes.len(), 1);
    assert_eq!(
        cli.state.process_list.processes[0].definition.folder_name,
        "relative_cwd"
    );
    Ok(())
}

fn import_graphs_relative_path_reports_lookup_locations() -> Result<()> {
    let root = cli_state_path("import_graphs_relative_path_reports_lookup_locations");
    clean_test(&root);
    fs::create_dir_all(&root)?;
    let cwd = root.join("cwd");
    let _cwd_guard = CurrentDirGuard::enter(&cwd)?;

    let mut cli = get_test_cli(None, root.join("state"), None, true)?;
    let missing_name = "missing_relative.dot";
    let err = cli
        .run_command(&format!("import graphs ./{missing_name}"))
        .unwrap_err();
    let error_text = format!("{err:?}");

    assert!(error_text.contains("Could not find graph file"));
    assert!(error_text.contains(&root.join(missing_name).display().to_string()));
    assert!(error_text.contains(&cwd.join(missing_name).display().to_string()));
    Ok(())
}

fn remove_processes_with_process_selector_removes_only_that_process() -> Result<()> {
    let mut cli = new_cli("remove_processes_with_process_selector_removes_only_that_process")?;
    populate_generated_scalar_box_process(&mut cli)?;
    duplicate_loaded_process(&mut cli, "second_process");

    cli.run_command("remove processes -p second_process")?;

    assert_eq!(cli.state.process_list.processes.len(), 1);
    assert_eq!(
        cli.state.process_list.processes[0].definition.folder_name,
        "box"
    );
    Ok(())
}

fn remove_processes_with_integrand_selector_removes_only_that_integrand() -> Result<()> {
    let mut cli = new_cli("remove_processes_with_integrand_selector_removes_only_that_integrand")?;
    populate_generated_scalar_box_process(&mut cli)?;
    add_duplicate_integrand(&mut cli, 0, "virtual_copy");

    cli.run_command("remove processes -p box -i virtual_copy")?;

    let remaining = cli.state.process_list.processes[0]
        .collection
        .get_integrand_names();
    assert_eq!(remaining, vec!["scalar_box"]);
    Ok(())
}

fn remove_processes_without_integrand_selector_drops_the_selected_process() -> Result<()> {
    let mut cli =
        new_cli("remove_processes_without_integrand_selector_drops_the_selected_process")?;
    populate_generated_scalar_box_process(&mut cli)?;
    add_duplicate_integrand(&mut cli, 0, "virtual_copy");

    cli.run_command("remove processes -p box")?;

    assert!(cli.state.process_list.processes.is_empty());
    Ok(())
}

fn remove_processes_removing_last_integrand_drops_empty_process() -> Result<()> {
    let mut cli = new_cli("remove_processes_removing_last_integrand_drops_empty_process")?;
    populate_generated_scalar_box_process(&mut cli)?;

    cli.run_command("remove processes -p box -i scalar_box")?;

    assert!(cli.state.process_list.processes.is_empty());
    Ok(())
}

fn remove_processes_without_selectors_clears_everything() -> Result<()> {
    let mut cli = new_cli("remove_processes_without_selectors_clears_everything")?;
    populate_generated_scalar_box_process(&mut cli)?;
    duplicate_loaded_process(&mut cli, "second_process");
    add_duplicate_integrand(&mut cli, 0, "virtual_copy");

    cli.run_command("remove processes")?;

    assert!(cli.state.process_list.processes.is_empty());
    Ok(())
}

fn remove_processes_ignores_missing_processes_and_integrands() -> Result<()> {
    let mut cli = new_cli("remove_processes_ignores_missing_processes_and_integrands")?;
    populate_generated_scalar_box_process(&mut cli)?;

    cli.run_command("remove processes -p does_not_exist")?;
    cli.run_command("remove processes -p box -i does_not_exist")?;

    assert_eq!(cli.state.process_list.processes.len(), 1);
    assert_eq!(
        cli.state.process_list.processes[0]
            .collection
            .get_integrand_names(),
        vec!["scalar_box"]
    );
    Ok(())
}

fn remove_processes_rejects_integrand_without_process() -> Result<()> {
    let mut cli = new_cli("remove_processes_rejects_integrand_without_process")?;
    let err = cli.run_command("remove processes -i 1L").unwrap_err();
    assert!(format!("{err:?}").contains("--integrand-name requires --process"));
    Ok(())
}

#[test]
#[serial]
fn cli_stateful_workflow_behaviors() -> Result<()> {
    run_without_arguments_is_a_noop()?;
    run_records_executed_commands()?;
    run_prevalidation_is_all_or_nothing_for_inline_commands()?;
    run_prevalidation_is_all_or_nothing_for_nested_block_failures()?;
    nested_run_records_executed_commands_once()?;
    run_command_block_substitutes_define_variables()?;
    run_command_blocks_share_define_environment_and_allow_unused_keys()?;
    nested_run_inherits_and_can_override_define_environment()?;
    run_define_placeholders_are_prevalidated_before_execution()?;
    run_inline_commands_can_use_define_variables()?;
    boot_run_history_merges_blocks_and_persists_commands_once()?;
    boot_run_history_rejects_conflicting_block_redefinitions()?;
    boot_run_history_allows_conflicting_redefinitions_after_confirmation()?;
    boot_run_history_cancelled_conflicting_redefinitions_leave_state_untouched()?;
    boot_run_history_conflicting_redefinitions_error_in_read_only_mode()?;
    boot_run_history_allows_identical_semantic_redefinitions()?;
    booting_existing_state_with_mismatched_frozen_settings_forces_read_only_and_warns()?;
    boot_mismatch_uses_conflicting_boot_blocks_read_only()?;
    command_block_definition_mode_defers_execution_and_omits_history()?;
    run_command_block_then_inline_quit_preserves_override_flag()?;
    finish_commands_block_without_active_definition_fails()?;
    start_commands_block_cannot_nest()?;
    quit_during_block_definition_dismisses_the_pending_block()?;
    ctrl_c_dismisses_pending_block_definition()?;
    ctrl_d_dismisses_pending_block_definition()?;
    recursive_run_limit_is_enforced()?;
    import_graphs_relative_path_prefers_active_state_root_parent()?;
    import_graphs_relative_path_falls_back_to_current_working_directory()?;
    import_graphs_relative_path_reports_lookup_locations()?;
    remove_processes_with_process_selector_removes_only_that_process()?;
    remove_processes_with_integrand_selector_removes_only_that_integrand()?;
    remove_processes_without_integrand_selector_drops_the_selected_process()?;
    remove_processes_removing_last_integrand_drops_empty_process()?;
    remove_processes_without_selectors_clears_everything()?;
    remove_processes_ignores_missing_processes_and_integrands()?;
    remove_processes_rejects_integrand_without_process()?;
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
fn run_history_load_accepts_command_block_templates_that_need_late_parsing() {
    let run_card_path = run_card_path(
        "run_history_load_accepts_command_block_templates_that_need_late_parsing.toml",
    );
    if let Some(parent) = run_card_path.parent() {
        fs::create_dir_all(parent).unwrap();
    }
    fs::write(
        &run_card_path,
        r#"
[[command_blocks]]
name = "bench_template"
commands = ["bench -s {samples} -p box -i default -c 1"]
"#,
    )
    .unwrap();

    let run_history = RunHistory::load(&run_card_path).unwrap();

    assert!(run_history.command_blocks[0].commands[0].is_template());
    assert_eq!(
        run_history
            .command_block_placeholder_names("bench_template")
            .into_iter()
            .collect::<Vec<_>>(),
        vec!["samples".to_string()]
    );
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
