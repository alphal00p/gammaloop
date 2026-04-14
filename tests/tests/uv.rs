use color_eyre::Result;
use gammaloop_integration_tests::{get_test_cli, get_tests_workspace_path};
#[test]
fn soft_ct_se() -> Result<()> {
    let state = get_test_cli(
        Some("dgse.toml".into()),
        get_tests_workspace_path().join("dgse"),
        None,
        false,
    )?;
    Ok(())
}
