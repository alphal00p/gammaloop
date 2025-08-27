use _gammaloop::{
    cli::{
        state::{RunHistory, State},
        Cli,
    },
    initialisation::test_initialise,
    utils::test_utils::load_generic_model,
};
use color_eyre::Result;
use std::env;
use tracing::debug;

fn new_cli_for_test() -> (Cli, State) {
    let state_path = if let Ok(user_specified_state_path) = env::var("TESTS_GAMMALOOP_STATE_PATH") {
        std::path::PathBuf::from(user_specified_state_path)
    } else {
        env::temp_dir().join("tests_gammaloop_state")
    };
    debug!("Using gammaloop state path: {}", state_path.display());

    let mut state = State::new_test(state_path.clone());
    state.model = load_generic_model("sm");

    (state.new_test_cli(), state)
}

fn clean_test(state: &State) {
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
}

#[test]
fn qqx_aaa_subtracted_nlo_amplitude_test() -> Result<()> {
    test_initialise()?;

    let mut cmds: RunHistory = RunHistory::from_file_yaml(
        "./tests/qqx_aaa_amplitude/qqx_aaa_subtracted_nlo_amplitude.yaml",
    )?;

    let (mut cli, mut state) = new_cli_for_test();

    cmds.run(&mut cli, &mut state)?;

    clean_test(&state);

    Ok(())
}
