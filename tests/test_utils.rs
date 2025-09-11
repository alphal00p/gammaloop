use color_eyre::Result;

use gammaloop_api::{
    state::{RunHistory, State},
    Cli, Commands,
};

use gammalooprs::{
    initialisation::initialise,
    utils::{test_utils::load_generic_model, F},
};
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
    pub cli: Cli,
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

pub(crate) fn run_run_card(
    run_card_path: impl AsRef<Path>,
    state_path: impl AsRef<Path>,
    log_file_name: Option<String>,
) -> Result<CLIState> {
    let (mut cli, mut state) = new_cli_for_test(state_path, log_file_name);
    let mut cmds: RunHistory = run_card(run_card_path)?;
    let mut global_settings = cmds.global_settings.clone();
    let mut default_runtime_settings = cmds.default_runtime_settings.clone();

    _ = cmds.run(
        &mut cli,
        &mut state,
        &mut global_settings,
        &mut default_runtime_settings,
    )?;

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
) -> (Cli, State) {
    debug!(
        "Using gammaloop state path: {}",
        state_path.as_ref().display()
    );
    let mut state = State::new(state_path.as_ref().to_path_buf(), log_file_name);
    state.model = load_generic_model("sm");

    (state.new_test_cli(), state)
}

pub(crate) fn clean_test(state: &State) {
    if env::var("GAMMALOOP_TESTS_NO_CLEAN_STATE").is_err() {
        match std::fs::remove_dir_all(state.save_path.clone()) {
            Ok(()) => debug!(
                "Gammaloop state folder '{}' deleted successfully",
                state.save_path.display()
            ),
            Err(e) => debug!(
                "Error deleting Gammaloop state folder '{}'. Error: {}",
                state.save_path.display(),
                e
            ),
        }
    } else {
        warn!("Environment variable 'GAMMALOOP_TESTS_NO_CLEAN_STATE' is set so that the gammaloop state test folder '{}' is not cleaned.",state.save_path.display());
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
