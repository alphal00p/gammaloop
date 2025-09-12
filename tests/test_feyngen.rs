use color_eyre::Result;

mod test_utils;
use test_utils::{clean_test, get_test_cli, get_tests_workspace_path};

#[test]
fn simple_epem_ddx_generation() -> Result<()> {
    let mut cli = get_test_cli(
        Some("sm_load.toml"),
        get_tests_workspace_path().join("epem_ddx_generation"),
        Some("epem_ddx_generation".to_string()),
    )?;

    cli.run_command("generate amp e+ e- > d d~")?;

    clean_test(&cli);

    Ok(())
}
