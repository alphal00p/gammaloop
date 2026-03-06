//! Dedicated crate for integration tests that depend on `gammaloop-api`.
#![allow(dead_code)]
#![allow(unused_variables)]

use color_eyre::Result;

use gammaloop_api::{
    CLISettings,
    commands::{Commands, save::SaveState},
    state::{RunHistory, State, SyncSettings},
};

use gammalooprs::{initialisation::initialise, utils::test_utils::load_generic_model};
use std::{
    env,
    ops::{Deref, DerefMut},
    path::{Path, PathBuf},
    str::FromStr,
    sync::Once,
};
use tracing::{debug, warn};

static SET_WORKSPACE_CWD: Once = Once::new();

fn ensure_workspace_cwd() {
    SET_WORKSPACE_CWD.call_once(|| {
        let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
            .parent()
            .expect("integration-tests must live in <workspace>/tests")
            .to_path_buf();
        std::env::set_current_dir(&workspace_root).unwrap_or_else(|e| {
            panic!(
                "Failed to set working directory to workspace root '{}': {e}",
                workspace_root.display()
            )
        });
    });
}

pub fn run_card(resource_path: impl AsRef<Path>) -> Result<RunHistory> {
    ensure_workspace_cwd();
    initialise()?;

    RunHistory::load(PathBuf::from("./tests/resources/run_cards").join(resource_path))
}

pub const TESTS_WORKSPACE: &str = "./tests/workspace";

pub struct CLIState {
    pub state: State,
    pub cli_settings: CLISettings,
    pub default_runtime_settings: gammalooprs::settings::RuntimeSettings,
    pub run_history: RunHistory,
}

impl Deref for CLIState {
    type Target = State;

    fn deref(&self) -> &Self::Target {
        &self.state
    }
}

impl DerefMut for CLIState {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.state
    }
}

impl CLIState {
    pub fn run_command(&mut self, command: &str) -> Result<()> {
        let cmd = Commands::from_str(command)?;

        cmd.run(
            &mut self.state,
            &mut self.run_history,
            &mut self.cli_settings,
            &mut self.default_runtime_settings,
        )
        .map(|_| ())
    }
}

/// Yields a cli state for testing.
/// If a run card path is provided, the commands in the run card are executed first.
/// If clean is true, the state path is deleted before creating the cli state.
/// The resulting state will be saved in the state path specified
pub fn get_test_cli(
    run_card_path: Option<PathBuf>,
    state_path: impl AsRef<Path>,
    log_file_name: Option<String>,
    clean: bool,
) -> Result<CLIState> {
    ensure_workspace_cwd();
    if clean && state_path.as_ref().exists() {
        clean_test(state_path.as_ref());
    }
    let mut state = new_cli_for_test(state_path.as_ref(), log_file_name);

    let cmds = if let Some(rc_path) = run_card_path {
        let mut run_card_cmds = run_card(rc_path)?;
        let mut global_settings_for_running_cms = run_card_cmds.cli_settings.clone();
        let mut default_runtime_settings_for_running_cms =
            run_card_cmds.default_runtime_settings.clone();
        _ = run_card_cmds.run(
            &mut state,
            &mut global_settings_for_running_cms,
            &mut default_runtime_settings_for_running_cms,
        )?;
        run_card_cmds
    } else {
        RunHistory::default()
    };
    let mut global_settings = cmds.cli_settings.clone();
    global_settings.state_folder = state_path.as_ref().to_path_buf();
    global_settings.sync_settings()?;
    let default_runtime_settings = cmds.default_runtime_settings.clone();
    SaveState {
        override_state: Some(true),
        ..Default::default()
    }
    .save(
        &mut state,
        &cmds,
        &default_runtime_settings,
        &global_settings,
    )?;

    Ok(CLIState {
        state,
        cli_settings: global_settings,
        default_runtime_settings,
        run_history: cmds,
    })
}

pub fn new_cli_for_test(state_path: impl AsRef<Path>, log_file_name: Option<String>) -> State {
    debug!(
        "Using gammaloop state path: {}",
        state_path.as_ref().display()
    );
    let mut state = State::new(state_path.as_ref().to_path_buf(), log_file_name);
    state.model = load_generic_model("sm");

    state
}

pub fn clean_test(save_path: impl AsRef<Path>) {
    if env::var("GAMMALOOP_TESTS_NO_CLEAN_STATE").is_err() {
        match std::fs::remove_dir_all(save_path.as_ref()) {
            Ok(()) => debug!(
                "Gammaloop state folder '{}' deleted successfully",
                save_path.as_ref().display()
            ),
            Err(e) => debug!(
                "Error deleting Gammaloop state folder '{}'. Error: {}",
                save_path.as_ref().display(),
                e
            ),
        }
    } else {
        warn!(
            "Environment variable 'GAMMALOOP_TESTS_NO_CLEAN_STATE' is set so that the gammaloop state test folder '{}' is not cleaned.",
            save_path.as_ref().display()
        );
    }
}

pub fn get_tests_workspace_path() -> PathBuf {
    if let Ok(user_specified_state_path) = env::var("TESTS_GAMMALOOP_STATE_PATH") {
        if user_specified_state_path.eq_ignore_ascii_case("AUTO") {
            env::temp_dir().join("gammaloop_tests_workspace")
        } else {
            std::path::PathBuf::from(user_specified_state_path)
        }
    } else {
        PathBuf::from(TESTS_WORKSPACE)
    }
}
