use color_eyre::Result;

use gammaloop_api::{
    inspect::Inspect,
    integrate::Integrate,
    state::{RunHistory, State},
    Cli, Commands,
};

use gammalooprs::{
    initialisation::initialise,
    status_info,
    utils::{test_utils::load_generic_model, F},
};
use insta::assert_snapshot;
use momtrop::assert_approx_eq;
use std::{
    env,
    ops::{ControlFlow, Deref, DerefMut},
    path::{Path, PathBuf},
    str::FromStr,
};
use tracing::{debug, warn};

fn run_card(resource_path: impl AsRef<Path>) -> Result<RunHistory> {
    initialise()?;

    RunHistory::load(PathBuf::from("./tests/resources/run_cards").join(resource_path))
}

const TESTS_WORKSPACE: &str = "./tests/workspace";

struct CLIState {
    state: State,
    cli: Cli,
    global_settings: gammalooprs::settings::GlobalSettings,
    default_runtime_settings: gammalooprs::settings::RuntimeSettings,
    run_history: RunHistory,
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
    fn run_command(&mut self, command: &str) -> Result<()> {
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

fn run_run_card(
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

fn new_cli_for_test(state_path: impl AsRef<Path>, log_file_name: Option<String>) -> (Cli, State) {
    debug!(
        "Using gammaloop state path: {}",
        state_path.as_ref().display()
    );
    let mut state = State::new(state_path.as_ref().to_path_buf(), log_file_name);
    state.model = load_generic_model("sm");

    (state.new_test_cli(), state)
}

fn clean_test(state: &State) {
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

fn get_tests_workspace_path() -> PathBuf {
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

#[test]
fn qqx_aaa_subtracted_nlo_amplitude_test() -> Result<()> {
    let mut state = run_run_card(
        "qqx_aaa_subtracted_nlo_amplitude.toml",
        get_tests_workspace_path().join("qqx_aaa_subtracted_nlo_amplitude"),
        None,
    )?;

    let a = Inspect {
        process_id: 0,
        integrand_name: "qqx_aaa_subtracted".to_string(),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        ..Default::default()
    }
    .run(&mut state)?;

    assert_snapshot!(format!("{a:.8e}"),@"(6.923469780323106e-4+-6.181987787694587e-4i)");

    Ok(())
}

#[test]
fn trees() -> Result<()> {
    let mut cli = run_run_card(
        "trees/qqx_aaa.toml",
        get_tests_workspace_path().join("qqx_aaa_tree"),
        None,
    )?;

    let a = Inspect {
        process_id: 0,
        integrand_name: "qqx_aaa".to_string(),
        ..Default::default()
    }
    .run(&mut cli)?;

    assert_snapshot!(format!("{a:.8e}"),@"(1.4727604164105595e-4+-1.1503139369130214e-3i)");

    clean_test(&cli);
    Ok(())
}

#[test]
fn photons_1l() -> Result<()> {
    let (mut state, _, _, _) = run_run_card("photons.toml", "./tests/photons")?;
    let inspect = Inspect {
        process_id: 0,
        process_name: "physical_1L_6photons".to_string(),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        ..Default::default()
    }
    .run(&mut state)?;

    assert_snapshot!(format!("{inspect:.8e}"),@"(-1.1779656944253678e-15+9.093348564643939e-16i)");

    let integrate = Integrate {
        process_id: 0,
        process_name: "physical_1L_6photons".to_string(),
        result_path: "./tests/photons/integration_workspace/integration_results.yaml".into(),
        workspace_path: "./tests/photons/integration_workspace/".into(),
        n_cores: 6,
        target: None,
        restart: true,
    };

    let integration_result = integrate.run(&mut state)?;
    Ok((()))
}

#[test]
fn test_grouped_subtraction() -> Result<()> {
    let mut cli = run_run_card(
        "test_grouped_subtraction.toml",
        get_tests_workspace_path().join("test_grouped_subtraction"),
        None,
    )?;

    let int1 = Integrate {
        process_id: 0,
        integrand_name: "no_group".to_string(),
        result_path: get_tests_workspace_path().join(
            "test_grouped_subtraction/integration_workspace_no_group/integration_results.yaml",
        ),
        n_cores: 1,
        workspace_path: get_tests_workspace_path()
            .join("test_grouped_subtraction/integration_workspace_no_group/"),
        target: None,
        restart: true,
    };
    let int2 = Integrate {
        process_id: 1,
        integrand_name: "group".to_string(),
        result_path: get_tests_workspace_path()
            .join("test_grouped_subtraction/integration_workspace_group/integration_results.yaml"),
        n_cores: 1,
        workspace_path: get_tests_workspace_path()
            .join("test_grouped_subtraction/integration_workspace_group/"),
        target: None,
        restart: true,
    };

    let integration_results_no_group = int1.run(&mut cli)?;
    let integration_results_group = int2.run(&mut cli)?;

    status_info!("No group result: {:#?}", integration_results_no_group);
    status_info!("Group result: {:#?}", integration_results_group);

    assert_approx_eq(
        &integration_results_group.result.re,
        &integration_results_no_group.result.re,
        &F(1e-1),
    );
    assert_approx_eq(
        &integration_results_group.result.im,
        &integration_results_no_group.result.im,
        &F(1e-1),
    );

    clean_test(&cli);

    Ok(())
}

#[test]
fn scalar_box() -> Result<()> {
    let mut cli = run_run_card(
        "scalar_box.toml",
        get_tests_workspace_path().join("scalar_box"),
        Some("scalar_box".to_string()),
    )?;

    let a = Inspect {
        process_id: 0,
        integrand_name: "scalar_box".to_string(),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        discrete_dim: vec![0],
        ..Default::default()
    }
    .run(&mut cli)?;

    assert_snapshot!(format!("{a:.8e}"),@"(1.3485885914334373e-3+0e0i)");

    cli.run_command("set process 0 scalar_box kv General.enable_cache=false")?;
    let integral_no_cache = Integrate {
        process_id: 0,
        integrand_name: "scalar_box".to_string(),
        result_path: get_tests_workspace_path()
            .join("scalar_box/integration_workspace/integration_results.toml"),
        n_cores: 1,
        workspace_path: get_tests_workspace_path().join("scalar_box/integration_workspace"),
        target: None,
        restart: true,
    }
    .run(&mut cli)?;

    cli.run_command("set process 0 scalar_box kv General.enable_cache=true")?;
    let integral_with_cache = Integrate {
        process_id: 0,
        integrand_name: "scalar_box".to_string(),
        result_path: get_tests_workspace_path()
            .join("scalar_box/integration_workspace/integration_results.toml"),
        n_cores: 1,
        workspace_path: get_tests_workspace_path().join("scalar_box/integration_workspace"),
        target: None,
        restart: true,
    }
    .run(&mut cli)?;

    status_info!("Integral result without caching: {:#?}", integral_no_cache);
    status_info!(
        "Integral result with caching.  : {:#?}",
        integral_with_cache
    );

    assert_approx_eq(
        &integral_no_cache.result.re,
        &integral_with_cache.result.re,
        &F(1e-4),
    );
    assert_approx_eq(
        &integral_no_cache.result.im,
        &integral_with_cache.result.im,
        &F(1e-4),
    );

    clean_test(&cli);

    Ok(())
}

#[test]
fn scalar_epem_ddx_generation() -> Result<()> {
    let mut cli = run_run_card(
        "sm_load.toml",
        get_tests_workspace_path().join("epem_ddx_generation"),
        Some("epem_ddx_generation".to_string()),
    )?;

    cli.run_command("generate amp e+ e- > d d~")?;

    clean_test(&cli);

    Ok(())
}
