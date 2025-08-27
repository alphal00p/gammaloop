use _gammaloop::{
    cli::{
        state::{RunHistory, State},
        Cli,
    },
    initialisation::test_initialise,
    utils::test_utils::load_generic_model,
};
use color_eyre::Result;

#[test]
fn val_test() -> Result<()> {
    test_initialise()?;
    let mut run = RunHistory::from_file_yaml("./tests/babis_subtraction/valentin_run.yaml")?;
    let mut state = State::new_test("./tests/babis_subtraction/gammaloop_state".into());
    state.model = load_generic_model("sm");
    let mut cli = state.new_test_cli();

    run.run(&mut cli, &mut state)?;

    match std::fs::remove_dir_all("./tests/babis_subtraction/gammaloop_state") {
        Ok(()) => println!("Directory deleted successfully"),
        Err(e) => eprintln!("Error deleting directory: {}", e),
    }
    Ok(())
}
