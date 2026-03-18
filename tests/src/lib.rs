//! Dedicated crate for integration tests that depend on `gammaloop-api`.
#![allow(dead_code)]
#![allow(unused_variables)]

use color_eyre::Result;

use gammaloop_api::{
    CLISettings,
    commands::CommandExecution,
    commands::save::SaveState,
    session::{CliSession, CliSessionState},
    state::{RunHistory, State, SyncSettings},
};

use gammalooprs::{
    initialisation::initialise, integrands::HasIntegrand, utils::test_utils::load_generic_model,
};
use std::{
    env,
    ffi::OsString,
    ops::{ControlFlow, Deref, DerefMut},
    path::{Path, PathBuf},
    sync::Once,
};
use tracing::{debug, warn};

static SET_WORKSPACE_CWD: Once = Once::new();

fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("integration-tests must live in <workspace>/tests")
        .to_path_buf()
}

fn ensure_workspace_cwd() {
    SET_WORKSPACE_CWD.call_once(|| {
        let workspace_root = workspace_root();
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

pub const TESTS_ARTIFACTS: &str = "tests/artifacts";

#[derive(Clone)]
pub struct CLIState {
    pub state: State,
    pub cli_settings: CLISettings,
    pub default_runtime_settings: gammalooprs::settings::RuntimeSettings,
    pub run_history: RunHistory,
    pub session_state: CliSessionState,
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
    fn with_session<T>(&mut self, f: impl FnOnce(&mut CliSession<'_>) -> Result<T>) -> Result<T> {
        let mut session = CliSession::new(
            &mut self.state,
            &mut self.run_history,
            &mut self.cli_settings,
            &mut self.default_runtime_settings,
            &mut self.session_state,
        );
        f(&mut session)
    }

    pub fn run_command_flow(&mut self, command: &str) -> Result<CommandExecution> {
        let command_history = gammaloop_api::state::CommandHistory::from_raw_string(command)?;
        self.with_session(|session| session.execute_command(command_history))
    }

    pub fn run_command(&mut self, command: &str) -> Result<()> {
        self.run_command_flow(command).map(|_| ())
    }

    pub fn apply_boot_run_history(
        &mut self,
        run_history: &RunHistory,
    ) -> Result<ControlFlow<SaveState>> {
        self.with_session(|session| session.apply_boot_run_history(run_history, run_history, false))
    }

    pub fn apply_boot_run_history_with_conflict_resolution(
        &mut self,
        run_history: &RunHistory,
        accept_conflicts: bool,
    ) -> Result<ControlFlow<SaveState>> {
        self.with_session(|session| {
            session.apply_boot_run_history_with_conflict_resolver(
                run_history,
                run_history,
                false,
                |_| Ok(accept_conflicts),
            )
        })
    }

    pub fn dismiss_pending_commands_block(&mut self, trigger: &str) -> bool {
        self.with_session(|session| Ok(session.dismiss_pending_commands_block(trigger)))
            .unwrap_or(false)
    }

    pub fn save_state(&mut self) -> Result<()> {
        SaveState {
            override_state: Some(true),
            ..Default::default()
        }
        .save(
            &mut self.state,
            &self.run_history,
            &self.default_runtime_settings,
            &self.cli_settings,
        )
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
    global_settings.state.folder = state_path.as_ref().to_path_buf();
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
        session_state: CliSessionState::default(),
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
            env::temp_dir().join("gammaloop_test_artifacts")
        } else {
            std::path::PathBuf::from(user_specified_state_path)
        }
    } else {
        workspace_root().join(TESTS_ARTIFACTS)
    }
}

pub struct EnvVarGuard {
    key: &'static str,
    previous: Option<OsString>,
}

impl EnvVarGuard {
    pub fn set(key: &'static str, value: &str) -> Self {
        let previous = env::var_os(key);
        unsafe {
            env::set_var(key, value);
        }
        Self { key, previous }
    }
}

impl Drop for EnvVarGuard {
    fn drop(&mut self) {
        match &self.previous {
            Some(value) => unsafe {
                env::set_var(self.key, value);
            },
            None => unsafe {
                env::remove_var(self.key);
            },
        }
    }
}

pub fn enable_lu_e2e_hack() -> EnvVarGuard {
    EnvVarGuard::set("GL_LU_E2E_HACK", "1")
}

pub fn new_test_artifact_dir(name: &str) -> Result<PathBuf> {
    let path = get_tests_workspace_path().join(name);
    if path.exists() {
        clean_test(&path);
    }
    std::fs::create_dir_all(&path)?;
    Ok(path)
}

pub fn run_commands(cli: &mut CLIState, commands: &[&str]) -> Result<()> {
    for command in commands {
        cli.run_command(command)?;
    }
    Ok(())
}

pub fn setup_sm_differential_lu_cli(test_name: &str) -> Result<CLIState> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join(test_name),
        Some(test_name.to_string()),
        true,
    )?;
    run_commands(
        &mut cli,
        &[
            "import model sm-default",
            r#"set default-runtime kv kinematics.externals='{"type":"constant","data":{"momenta":[[32.0,0.0,0.0,32.0],[32.0,0.0,0.0,-32.0]],"helicities":[1,1]}}'"#,
            "set default-runtime kv subtraction.disable_threshold_subtraction=true",
            "generate e+ e- > d d~ g | e- a d g QED^2==4 [{{2}} QCD=0] --numerator-grouping group_identical_graphs_up_to_sign --clear-existing-processes --only-diagrams",
            "generate",
        ],
    )?;
    Ok(cli)
}

pub fn default_xspace_point(cli: &CLIState) -> Result<Vec<f64>> {
    let (process_id, integrand_name) = cli.state.find_integrand_ref(None, None)?;
    let integrand = cli
        .state
        .process_list
        .get_integrand(process_id, &integrand_name)?;
    let n_dim = integrand.get_n_dim();
    let seed = [0.17, 0.31, 0.53, 0.23, 0.41, 0.67];
    Ok((0..n_dim).map(|index| seed[index % seed.len()]).collect())
}

pub fn default_momentum_space_point(cli: &CLIState) -> Result<Vec<f64>> {
    let (process_id, integrand_name) = cli.state.find_integrand_ref(None, None)?;
    let integrand = cli
        .state
        .process_list
        .get_integrand(process_id, &integrand_name)?;
    let n_dim = integrand.get_n_dim();
    let seed = [0.11, -0.07, 0.19, -0.13, 0.05, 0.29];
    Ok((0..n_dim).map(|index| seed[index % seed.len()]).collect())
}
