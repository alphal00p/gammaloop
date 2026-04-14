use color_eyre::Result;
use gammaloop_api::commands::{Profile, profile::InfraRedProfile};
use gammaloop_integration_tests::{get_test_cli, get_tests_workspace_path};

#[test]
fn soft_ct_se() -> Result<()> {
    let mut state = get_test_cli(
        Some("dgse.toml".into()),
        get_tests_workspace_path().join("dgse"),
        None,
        false,
    )?;

    let res = Profile::InfraRed(InfraRedProfile {
        select: Some("se S(e0)".into()),
        ..Default::default()
    })
    .run(&mut state.state, &state.cli_settings)?;

    assert!(res.unwrap_ir().all_passed);
    Ok(())
}
