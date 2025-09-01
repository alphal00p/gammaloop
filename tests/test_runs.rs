use _gammaloop::{
    cli::{
        state::{RunHistory, State},
        Cli,
    },
    initialisation::test_initialise,
    settings::runtime::IntegrationResult,
    status_info,
    utils::{test_utils::load_generic_model, F},
};
use color_eyre::Result;
use momtrop::assert_approx_eq;
use std::{env, path::Path};
use tracing::{debug, warn};

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

    let mut state = State::new_test(state_path.clone());
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
    test_initialise()?;

    let mut cmds: RunHistory = RunHistory::from_file_yaml(
        "./tests/qqx_aaa_amplitude/qqx_aaa_subtracted_nlo_amplitude.yaml",
    )?;

    let (mut cli, mut state) = new_cli_for_test("./tests/qqx_aaa_amplitude");

    let _ = cmds.run(&mut cli, &mut state)?;

    clean_test(&state);

    Ok(())
}

#[test]
fn test_grouped_subtraction() -> Result<()> {
    test_initialise()?;

    let mut cmds: RunHistory = RunHistory::from_file_yaml(
        "./tests/test_grouped_subtraction/test_grouped_subtraction.yaml",
    )?;

    let (mut cli, mut state) = new_cli_for_test("./tests/test_grouped_subtraction");

    let _ = cmds.run(&mut cli, &mut state)?;

    let integration_results_no_group: IntegrationResult = serde_yaml::from_str(
        &std::fs::read_to_string(
            "./tests/test_grouped_subtraction/integration_workspace_no_group/integration_results.yaml",
        )
        .unwrap(),
    ).unwrap();

    let integration_results_group: IntegrationResult = serde_yaml::from_str(
        &std::fs::read_to_string(
            "./tests/test_grouped_subtraction/integration_workspace_group/integration_results.yaml",
        )
        .unwrap(),
    )
    .unwrap();

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
