use color_eyre::Result;
use gammaloop_api::commands::{
    Profile,
    integrate::Integrate,
    profile::{InfraRedProfile, UltraVioletProfile},
};
use gammaloop_api::state::ProcessRef;
use gammaloop_integration_tests::{
    clean_test, get_example_cli, get_test_cli, get_tests_workspace_path,
};
use gammalooprs::utils::F;
use spenso::algebra::complex::Complex;

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

    cli.run_command("run generate")?;
    cli.run_command(
        "set process kv integrator.target_relative_accuracy=0.005 \
         integrator.n_increase=0 integrator.n_start=50000 \
         integrator.n_max=1000000 integrator.seed=1337",
    )?;
    let integration_result = Integrate {
        process: vec![
            ProcessRef::Unqualified("bubble_dod1_no_integrated_UV".to_string()),
            ProcessRef::Unqualified("bubble_dod1".to_string()),
        ],
        integrand_name: vec![
            "scalar_bubble_below_thres".to_string(),
            "scalar_bubble_below_thres".to_string(),
        ],
        workspace_path: Some(
            get_tests_workspace_path().join("dod1_bubble/integration_workspace_below_threshold"),
        ),
        n_cores: Some(1),
        target: vec![
            "bubble_dod1_no_integrated_UV@scalar_bubble_below_thres=7.358320108607984e-3,0.0"
                .to_string(),
            "bubble_dod1@scalar_bubble_below_thres=1.3514937842276269e-3,0.0".to_string(),
        ],
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let no_integrated_slot = integration_result
        .slot("bubble_dod1_no_integrated_UV@scalar_bubble_below_thres")
        .expect("bubble_dod1_no_integrated_UV below-threshold slot should exist");
    assert!(
        no_integrated_slot
            .integral
            .is_compatible_with_target(Complex::new_re(F(7.358320108607984e-3)), 2),
        "Integration result {} not compatible with target {}",
        no_integrated_slot.integral,
        Complex::new_re(F(7.358320108607984e-3))
    );

    let integrated_slot = integration_result
        .slot("bubble_dod1@scalar_bubble_below_thres")
        .expect("bubble_dod1 below-threshold slot should exist");
    assert!(
        integrated_slot
            .integral
            .is_compatible_with_target(Complex::new_re(F(1.3514937842276269e-3)), 2),
        "Integration result {} not compatible with target {}",
        integrated_slot.integral,
        Complex::new_re(F(1.3514937842276269e-3))
    );

    let res = Profile::UltraViolet(UltraVioletProfile {
        process: Some(ProcessRef::Unqualified("bubble_dod1".to_string())),
        integrand_name: Some("scalar_bubble_below_thres".to_string()),
        min_scale_exponent: 3.0,
        max_scale_exponent: 5.0,
        per_orientation: true,
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?;
    let uv = res.unwrap_uv();
    assert_eq!(uv.pass_fail(-0.9).failed, 0);

    let res = Profile::UltraViolet(UltraVioletProfile {
        process: Some(ProcessRef::Unqualified(
            "bubble_dod1_no_integrated_UV".to_string(),
        )),
        integrand_name: Some("scalar_bubble_below_thres".to_string()),
        min_scale_exponent: 3.0,
        max_scale_exponent: 5.0,
        per_orientation: true,
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?;
    let uv = res.unwrap_uv();
    assert_eq!(uv.pass_fail(-0.9).failed, 0);
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

    cli.run_command("run generate integrate_bubble_below_threshold")?;
    let res = Profile::UltraViolet(UltraVioletProfile {
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?;
    let uv = res.unwrap_uv();
    assert_eq!(uv.pass_fail(-0.9).failed, 0);
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

    cli.run_command("run generate integrate_bubble_below_threshold")?;
    let res = Profile::UltraViolet(UltraVioletProfile {
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?;
    let uv = res.unwrap_uv();
    assert_eq!(uv.pass_fail(-0.9).failed, 0);
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
