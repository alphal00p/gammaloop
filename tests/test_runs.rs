use color_eyre::Result;
use gammaloop_api::{
    inspect::Inspect,
    integrate::Integrate,
    state::{RunHistory, State},
    Cli,
};
use gammalooprs::{
    initialisation::{initialise, test_initialise},
    settings::runtime::IntegrationResult,
    status_info,
    utils::{test_utils::load_generic_model, F},
};
use insta::assert_snapshot;
use momtrop::assert_approx_eq;
use std::{
    env,
    path::{Path, PathBuf},
};
use tracing::{debug, warn};

fn run_card(resource_path: impl AsRef<Path>) -> Result<RunHistory> {
    initialise()?;
    RunHistory::new(PathBuf::from("./tests/resources/runs/").join(resource_path))
}

fn run_run_card(
    run_card_path: impl AsRef<Path>,
    state_path: impl AsRef<Path>,
) -> Result<(
    State,
    Cli,
    gammalooprs::settings::GlobalSettings,
    gammalooprs::settings::RuntimeSettings,
)> {
    let mut cmds: RunHistory = run_card(run_card_path)?;
    let mut global_settings = cmds.global_settings.clone();
    let mut default_runtime_settings = cmds.default_runtime_settings.clone();
    let (mut cli, mut state) = new_cli_for_test(state_path);

    let _ = cmds.run(
        &mut cli,
        &mut state,
        &mut global_settings,
        &mut default_runtime_settings,
    )?;

    cli.save(
        &mut state,
        &mut cmds,
        &mut default_runtime_settings,
        &mut global_settings,
        None,
        false,
    )?;

    Ok((state, cli, global_settings, default_runtime_settings))
}

fn new_cli_for_test(test_path: impl AsRef<Path>) -> (Cli, State) {
    let state_path = if let Ok(user_specified_state_path) = env::var("TESTS_GAMMALOOP_STATE_PATH") {
        if user_specified_state_path.to_ascii_uppercase() == "AUTO" {
            env::temp_dir().join("tests_gammaloop_state")
        } else {
            std::path::PathBuf::from(user_specified_state_path)
        }
    } else {
        test_path.as_ref().join("tests_gammaloop_state")
    };
    debug!("Using gammaloop state path: {}", state_path.display());

    let mut state = State::new(state_path.clone(), None);
    state.model = load_generic_model("sm");

    (state.new_test_cli(), state)
}

fn clean_test(state: &State) {
    if let Err(_) = env::var("GAMMALOOP_TESTS_NO_CLEAN_STATE") {
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

#[test]
fn qqx_aaa_subtracted_nlo_amplitude_test() -> Result<()> {
    let (state, _, _, _) = run_run_card(
        "qqx_aaa_subtracted_nlo_amplitude.toml",
        "./tests/qqx_aaa_subtracted_nlo_amplitude",
    )?;

    Ok(())
}

#[test]
fn trees() -> Result<()> {
    let (mut state, _, _, _) = run_run_card("trees/qqx_aaa.toml", "./tests/qqx_aaa_tree")?;
    state.process_list.warm_up(&state.model)?;
    let a = Inspect {
        process_id: 0,
        process_name: "qqx_aaa".to_string(),
        ..Default::default()
    }
    .run(&mut state)?;

    assert_snapshot!(format!("{a:.8e}"),@"(1.4727604164105595e-4+-1.1503139369130214e-3i)");

    Ok(())
}

#[test]
fn test_grouped_subtraction() -> Result<()> {
    let (mut state, _, _, _) = run_run_card(
        "trees/test_grouped_subtraction.toml",
        "./tests/test_grouped_subtraction",
    )?;

    let int1 = Integrate{
        process_id: 0,
        process_name: "no_group".to_string(),
        result_path: "./tests/test_grouped_subtraction/integration_workspace_no_group/integration_results.yaml".into(),
        n_cores: 1,
        workspace_path: "./tests/test_grouped_subtraction/integration_workspace_no_group/".into(),
        target: None,
        restart: true,
    };
    let int2 = Integrate {
        process_id: 1,
        process_name: "group".to_string(),
        result_path:
            "./tests/test_grouped_subtraction/integration_workspace_group/integration_results.yaml"
                .into(),
        n_cores: 1,
        workspace_path: "./tests/test_grouped_subtraction/integration_workspace_group/".into(),
        target: None,
        restart: true,
    };

    let integration_results_no_group = int1.run(&mut state)?;
    let integration_results_group = int2.run(&mut state)?;

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

    clean_test(&state);

    Ok(())
}
