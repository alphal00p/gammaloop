use _gammaloop::{
    cli::{
        state::{RunHistory, State},
        Cli,
    },
    initialisation::test_initialise,
    utils::test_utils::load_generic_model,
};
use color_eyre::Result;
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

    cmds.run(&mut cli, &mut state)?;

    clean_test(&state);

    Ok(())
}
