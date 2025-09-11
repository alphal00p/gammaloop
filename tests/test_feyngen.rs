use color_eyre::Result;

mod test_utils;
use test_utils::{clean_test, get_tests_workspace_path, run_run_card};

#[test]
fn simple_epem_ddx_generation() -> Result<()> {
    let mut cli = run_run_card(
        "sm_load.toml",
        get_tests_workspace_path().join("epem_ddx_generation"),
        Some("epem_ddx_generation".to_string()),
    )?;

    cli.run_command("generate amp e+ e- > d d~")?;

    clean_test(&cli);

    Ok(())
}
