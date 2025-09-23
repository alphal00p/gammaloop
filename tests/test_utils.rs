#![allow(dead_code)]
#![allow(unused_variables)]

use color_eyre::Result;

use gammaloop_api::{
    state::{RunHistory, State},
    Repl, Commands,
};

use gammalooprs::{initialisation::initialise, utils::test_utils::load_generic_model};
use std::{
    env,
    ops::{Deref, DerefMut},
    path::{Path, PathBuf},
    str::FromStr,
};
use tracing::{debug, warn};

pub(crate) fn run_card(resource_path: impl AsRef<Path>) -> Result<RunHistory> {
    initialise()?;

    RunHistory::load(PathBuf::from("./tests/resources/run_cards").join(resource_path))
}

pub(crate) const TESTS_WORKSPACE: &str = "./tests/workspace";

pub(crate) struct CLIState {
    pub state: State,
    pub cli: Repl,
    pub global_settings: gammalooprs::settings::GlobalSettings,
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
    pub(crate) fn run_command(&mut self, command: &str) -> Result<()> {
        let cmd = Commands::from_str(command)?;
        self.cli
            .run_command(
                cmd,
                &mut self.state,
                &mut self.run_history,
                &mut self.global_settings,
                &mut self.default_runtime_settings,
            )
            .map(|_| ())
    }
}

pub(crate) fn get_test_cli(
    run_card_path: Option<PathBuf>,
    state_path: impl AsRef<Path>,
    log_file_name: Option<String>,
    clean: bool,
) -> Result<CLIState> {
    if clean && state_path.as_ref().exists() {
        clean_test(state_path.as_ref());
    }
    let (mut cli, mut state) = new_cli_for_test(state_path, log_file_name);
    let cmds = if let Some(rc_path) = run_card_path {
        let mut run_card_cmds = run_card(rc_path)?;
        let mut global_settings_for_running_cms = run_card_cmds.global_settings.clone();
        let mut default_runtime_settings_for_running_cms =
            run_card_cmds.default_runtime_settings.clone();
        _ = run_card_cmds.run(
            &mut cli,
            &mut state,
            &mut global_settings_for_running_cms,
            &mut default_runtime_settings_for_running_cms,
        )?;
        run_card_cmds
    } else {
        RunHistory::default()
    };
    let global_settings = cmds.global_settings.clone();
    let default_runtime_settings = cmds.default_runtime_settings.clone();
    cli.save(
        &mut state,
        &cmds,
        &default_runtime_settings,
        &global_settings,
        None,
        false,
    )?;

    Ok(CLIState {
        state,
        cli,
        global_settings,
        default_runtime_settings,
        run_history: cmds,
    })
}

pub(crate) fn new_cli_for_test(
    state_path: impl AsRef<Path>,
    log_file_name: Option<String>,
) -> (Repl, State) {
    debug!(
        "Using gammaloop state path: {}",
        state_path.as_ref().display()
    );
    let mut state = State::new(state_path.as_ref().to_path_buf(), log_file_name);
    state.model = load_generic_model("sm");

    (state.new_test_cli(), state)
}

pub(crate) fn clean_test(save_path: impl AsRef<Path>) {
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
        warn!("Environment variable 'GAMMALOOP_TESTS_NO_CLEAN_STATE' is set so that the gammaloop state test folder '{}' is not cleaned.",save_path.as_ref().display());
    }
}

pub(crate) fn get_tests_workspace_path() -> PathBuf {
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
