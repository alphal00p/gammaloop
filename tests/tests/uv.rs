use color_eyre::Result;
use gammaloop_api::commands::{
    Profile,
    integrate::Integrate,
    profile::{InfraRedProfile, UltraVioletProfile},
};
use gammaloop_api::state::ProcessRef;
use gammaloop_integration_tests::{
    CLIState, clean_test, get_example_cli, get_test_cli, get_tests_workspace_path,
};
use gammalooprs::utils::F;
use spenso::algebra::complex::Complex;

fn run_bubble_below_threshold_integration(
    cli: &mut CLIState,
    test_name: &str,
    no_integrated_process: &str,
    integrated_process: &str,
    no_integrated_target: f64,
    integrated_target: f64,
) -> Result<()> {
    cli.run_command("run generate")?;
    cli.run_command(
        "set process kv integrator.target_relative_accuracy=0.005 \
         integrator.n_increase=0 integrator.n_start=50000 \
         integrator.n_max=1000000 integrator.seed=1337",
    )?;
    let integration_result = Integrate {
        process: vec![
            ProcessRef::Unqualified(no_integrated_process.to_string()),
            ProcessRef::Unqualified(integrated_process.to_string()),
        ],
        integrand_name: vec![
            "scalar_bubble_below_thres".to_string(),
            "scalar_bubble_below_thres".to_string(),
        ],
        workspace_path: Some(
            get_tests_workspace_path()
                .join(format!("{test_name}/integration_workspace_below_threshold")),
        ),
        n_cores: Some(1),
        target: vec![
            format!("{no_integrated_process}@scalar_bubble_below_thres={no_integrated_target},0.0"),
            format!("{integrated_process}@scalar_bubble_below_thres={integrated_target},0.0"),
        ],
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let no_integrated_key = format!("{no_integrated_process}@scalar_bubble_below_thres");
    let no_integrated_slot = integration_result
        .slot(&no_integrated_key)
        .expect("below-threshold no-integrated slot should exist");
    assert!(
        no_integrated_slot
            .integral
            .is_compatible_with_target(Complex::new_re(F(no_integrated_target)), 2),
        "Integration result {} not compatible with target {}",
        no_integrated_slot.integral,
        Complex::new_re(F(no_integrated_target))
    );

    let integrated_key = format!("{integrated_process}@scalar_bubble_below_thres");
    let integrated_slot = integration_result
        .slot(&integrated_key)
        .expect("below-threshold integrated slot should exist");
    assert!(
        integrated_slot
            .integral
            .is_compatible_with_target(Complex::new_re(F(integrated_target)), 2),
        "Integration result {} not compatible with target {}",
        integrated_slot.integral,
        Complex::new_re(F(integrated_target))
    );

    Ok(())
}

fn assert_bubble_uv_profile_passes(cli: &mut CLIState, process: &str) -> Result<()> {
    let res = Profile::UltraViolet(UltraVioletProfile {
        process: Some(ProcessRef::Unqualified(process.to_string())),
        integrand_name: Some("scalar_bubble_below_thres".to_string()),
        min_scale_exponent: 4.0,
        max_scale_exponent: 6.0,
        per_orientation: true,
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?;
    let uv = res.unwrap_uv();
    assert_eq!(uv.pass_fail(-0.9).failed, 0);
    Ok(())
}

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

    let res = Profile::UltraViolet(UltraVioletProfile {
        ..Default::default()
    })
    .run(&mut state.state, &state.cli_settings)?;

    let uv = res.unwrap_uv();
    assert_eq!(uv.pass_fail(-0.9).failed, 0);
    Ok(())
}

#[test]
fn epem_a_tth_nlo_uv() -> Result<()> {
    let state_path = get_tests_workspace_path().join("epem_a_tth_nlo_example");
    let mut cli = get_example_cli(
        "epem_a_ttxh/NLO/epem_a_tth_NLO.toml",
        &["generate_diagrams"],
        Some(state_path.clone()),
        None,
        true,
    )?;

    cli.run_command("run generate_diagrams generate_integrands")?;
    let res = Profile::UltraViolet(UltraVioletProfile {
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?;

    let uv = res.unwrap_uv();
    assert_eq!(uv.pass_fail(-0.9).failed, 0);
    clean_test(&state_path);
    Ok(())
}

#[test]
fn dod1_bubble_uv() -> Result<()> {
    let mut cli = get_test_cli(
        Some("dod1_bubble.toml".into()),
        get_tests_workspace_path().join("dod1_bubble"),
        None,
        false,
    )?;

    run_bubble_below_threshold_integration(
        &mut cli,
        "dod1_bubble",
        "bubble_dod1_no_integrated_UV",
        "bubble_dod1",
        7.358320108607984e-3,
        1.3514937842276269e-3,
    )?;
    assert_bubble_uv_profile_passes(&mut cli, "bubble_dod1")?;
    assert_bubble_uv_profile_passes(&mut cli, "bubble_dod1_no_integrated_UV")?;
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn dod2_bubble_uv() -> Result<()> {
    let mut cli = get_test_cli(
        Some("dod2_bubble.toml".into()),
        get_tests_workspace_path().join("dod2_bubble"),
        None,
        false,
    )?;

    run_bubble_below_threshold_integration(
        &mut cli,
        "dod2_bubble",
        "bubble_dod2_no_integrated_UV",
        "bubble_dod2",
        -0.5940830828502411,
        1.2143596454658382e-2,
    )?;
    assert_bubble_uv_profile_passes(&mut cli, "bubble_dod2")?;
    assert_bubble_uv_profile_passes(&mut cli, "bubble_dod2_no_integrated_UV")?;
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn dod0_bubble_uv() -> Result<()> {
    let mut cli = get_test_cli(
        Some("dod0_bubble.toml".into()),
        get_tests_workspace_path().join("dod0_bubble"),
        None,
        false,
    )?;

    run_bubble_below_threshold_integration(
        &mut cli,
        "dod0_bubble",
        "bubble_no_integrated_UV",
        "bubble",
        1.471664021721597e-2,
        2.7029875684552542e-3,
    )?;
    assert_bubble_uv_profile_passes(&mut cli, "bubble")?;
    assert_bubble_uv_profile_passes(&mut cli, "bubble_no_integrated_UV")?;
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
