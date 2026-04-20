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
use gammalooprs::settings::runtime::IntegralEstimate;
use gammalooprs::utils::F;
use spenso::algebra::complex::Complex;

struct BubbleBelowThresholdResults {
    no_integrated: IntegralEstimate,
    integrated: IntegralEstimate,
}

fn run_bubble_below_threshold_integration(
    cli: &mut CLIState,
    test_name: &str,
    no_integrated_process: &str,
    integrated_process: &str,
) -> Result<BubbleBelowThresholdResults> {
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
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let no_integrated_key = format!("{no_integrated_process}@scalar_bubble_below_thres");
    let no_integrated_slot = integration_result
        .slot(&no_integrated_key)
        .expect("below-threshold no-integrated slot should exist");

    let integrated_key = format!("{integrated_process}@scalar_bubble_below_thres");
    let integrated_slot = integration_result
        .slot(&integrated_key)
        .expect("below-threshold integrated slot should exist");

    Ok(BubbleBelowThresholdResults {
        no_integrated: no_integrated_slot.integral.clone(),
        integrated: integrated_slot.integral.clone(),
    })
}

fn assert_bubble_below_threshold_targets(
    results: &BubbleBelowThresholdResults,
    no_integrated_target: f64,
    integrated_target: f64,
) {
    assert!(
        results
            .no_integrated
            .is_compatible_with_target(Complex::new_re(F(no_integrated_target)), 2),
        "Integration result {} not compatible with target {}",
        results.no_integrated,
        Complex::new_re(F(no_integrated_target))
    );
    assert!(
        results
            .integrated
            .is_compatible_with_target(Complex::new_re(F(integrated_target)), 2),
        "Integration result {} not compatible with target {}",
        results.integrated,
        Complex::new_re(F(integrated_target))
    );
}

fn assert_bubble_uv_profile_passes(cli: &mut CLIState, process: &str) -> Result<()> {
    let res = Profile::UltraViolet(UltraVioletProfile {
        process: Some(ProcessRef::Unqualified(process.to_string())),
        integrand_name: Some("scalar_bubble_below_thres".to_string()),
        min_scale_exponent: 4.0,
        max_scale_exponent: 5.0,
        n_points: 100,
        per_orientation: true,
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?;
    let uv = res.unwrap_uv();
    assert_eq!(uv.pass_fail(-0.9).failed, 0);
    Ok(())
}

fn assert_integrated_result_is_muv_invariant(
    cli: &mut CLIState,
    test_name: &str,
    integrated_process: &str,
    baseline: &IntegralEstimate,
    original_m_uv: f64,
    alternate_m_uv: f64,
) -> Result<()> {
    cli.run_command(&format!("set process kv general.m_uv={alternate_m_uv}"))?;
    cli.run_command(
        "set process kv integrator.target_relative_accuracy=0.005 \
         integrator.n_increase=0 integrator.n_start=50000 \
         integrator.n_max=1000000 integrator.seed=1337",
    )?;
    let integration_result = Integrate {
        process: vec![ProcessRef::Unqualified(integrated_process.to_string())],
        integrand_name: vec!["scalar_bubble_below_thres".to_string()],
        workspace_path: Some(get_tests_workspace_path().join(format!(
            "{test_name}/integration_workspace_below_threshold_muv_shifted"
        ))),
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let integrated_key = format!("{integrated_process}@scalar_bubble_below_thres");
    let integrated_slot = integration_result
        .slot(&integrated_key)
        .expect("below-threshold integrated shifted-m_uv slot should exist");
    assert!(
        integrated_slot
            .integral
            .is_compatible_with_result(baseline, 2),
        "Integration result {} not compatible with baseline {}",
        integrated_slot.integral,
        baseline
    );

    cli.run_command(&format!("set process kv general.m_uv={original_m_uv}"))?;

    Ok(())
}

#[test]
fn dod1_bubble_uv() -> Result<()> {
    let mut cli = get_test_cli(
        Some("dod1_bubble.toml".into()),
        get_tests_workspace_path().join("dod1_bubble"),
        None,
        true,
    )?;
    cli.run_command("run generate")?;

    let baseline = run_bubble_below_threshold_integration(
        &mut cli,
        "dod1_bubble",
        "bubble_dod1_no_integrated_UV",
        "bubble_dod1",
    )?;
    assert_integrated_result_is_muv_invariant(
        &mut cli,
        "dod1_bubble",
        "bubble_dod1",
        &baseline.integrated,
        20.0,
        7.0,
    )?;
    assert_bubble_below_threshold_targets(
        &baseline,
        7.358320108607984e-3,
        1.351_493_784_227_627e-3,
    );
    assert_bubble_uv_profile_passes(&mut cli, "bubble_dod1")?;
    assert_bubble_uv_profile_passes(&mut cli, "bubble_dod1_no_integrated_UV")?;
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn dod0_bubble_uv() -> Result<()> {
    let mut cli = get_test_cli(
        Some("dod0_bubble.toml".into()),
        get_tests_workspace_path().join("dod0_bubble"),
        None,
        true,
    )?;
    cli.run_command("run generate")?;
    let baseline = run_bubble_below_threshold_integration(
        &mut cli,
        "dod0_bubble",
        "bubble_no_integrated_UV",
        "bubble",
    )?;
    assert_integrated_result_is_muv_invariant(
        &mut cli,
        "dod0_bubble",
        "bubble",
        &baseline.integrated,
        20.0,
        7.0,
    )?;
    assert_bubble_below_threshold_targets(&baseline, 1.471664021721597e-2, 2.7029875684552542e-3);
    assert_bubble_uv_profile_passes(&mut cli, "bubble")?;
    assert_bubble_uv_profile_passes(&mut cli, "bubble_no_integrated_UV")?;
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

mod slow {
    use super::*;

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
}

#[test]
fn dod2_bubble_uv() -> Result<()> {
    let mut cli = get_test_cli(
        Some("dod2_bubble.toml".into()),
        get_tests_workspace_path().join("dod2_bubble"),
        None,
        true,
    )?;
    cli.run_command("run generate")?;
    assert_bubble_uv_profile_passes(&mut cli, "bubble_dod2")?;
    assert_bubble_uv_profile_passes(&mut cli, "bubble_dod2_no_integrated_UV")?;
    let baseline = run_bubble_below_threshold_integration(
        &mut cli,
        "dod2_bubble",
        "bubble_dod2_no_integrated_UV",
        "bubble_dod2",
    )?;
    assert_integrated_result_is_muv_invariant(
        &mut cli,
        "dod2_bubble",
        "bubble_dod2",
        &baseline.integrated,
        20.0,
        7.0,
    )?;
    assert_bubble_below_threshold_targets(&baseline, -0.5940830828502411, 1.2143596454658382e-2);
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
mod failing {
    use super::*;

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
}
