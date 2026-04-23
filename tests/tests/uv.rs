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

struct IntegratedUvTargets {
    no_integrated: f64,
    integrated: f64,
}

struct IntegratedUvCase<'a> {
    run_card: &'a str,
    test_name: &'a str,
    no_integrated_process: &'a str,
    integrated_process: &'a str,
    integrand_name: &'a str,
    targets: Option<IntegratedUvTargets>,
    min_change_sigma: Option<f64>,
}

struct IntegratedUvResults {
    no_integrated: IntegralEstimate,
    integrated: IntegralEstimate,
}

struct IntegratedUvCaseResult {
    graph: String,
    uv_profile_passed: bool,
    muv_invariance_passed: bool,
    target_passed: Option<bool>,
    target_threshold: Option<&'static str>,
    ct_change_passed: Option<bool>,
    ct_change_sigma: Option<f64>,
    ct_change_threshold: Option<f64>,
    error: Option<String>,
}

impl IntegratedUvCaseResult {
    fn all_passed(&self) -> bool {
        self.error.is_none()
            && self.uv_profile_passed
            && self.muv_invariance_passed
            && self.target_passed.unwrap_or(true)
            && self.ct_change_passed.unwrap_or(true)
    }
}

fn set_fast_deterministic_integrator(cli: &mut CLIState) -> Result<()> {
    cli.run_command(
        "set process kv integrator.target_relative_accuracy=0.001 \
         integrator.n_increase=0 integrator.n_start=50000 \
         integrator.n_max=2000000 integrator.seed=1337",
    )
}

fn run_integrated_uv_integration(
    cli: &mut CLIState,
    test_name: &str,
    no_integrated_process: &str,
    integrated_process: &str,
    integrand_name: &str,
) -> Result<IntegratedUvResults> {
    set_fast_deterministic_integrator(cli)?;
    let integration_result = Integrate {
        process: vec![
            ProcessRef::Unqualified(no_integrated_process.to_string()),
            ProcessRef::Unqualified(integrated_process.to_string()),
        ],
        integrand_name: vec![integrand_name.to_string(), integrand_name.to_string()],
        workspace_path: Some(get_tests_workspace_path().join(format!(
            "{test_name}/integration_workspace_{integrand_name}"
        ))),
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let no_integrated_key = format!("{no_integrated_process}@{integrand_name}");
    let no_integrated_slot = integration_result
        .slot(&no_integrated_key)
        .expect("no-integrated slot should exist");

    let integrated_key = format!("{integrated_process}@{integrand_name}");
    let integrated_slot = integration_result
        .slot(&integrated_key)
        .expect("integrated slot should exist");

    Ok(IntegratedUvResults {
        no_integrated: no_integrated_slot.integral.clone(),
        integrated: integrated_slot.integral.clone(),
    })
}

fn integrated_uv_targets_pass(
    results: &IntegratedUvResults,
    no_integrated_target: f64,
    integrated_target: f64,
) -> bool {
    results
        .no_integrated
        .is_compatible_with_target(Complex::new_re(F(no_integrated_target)), 2)
        && results
            .integrated
            .is_compatible_with_target(Complex::new_re(F(integrated_target)), 2)
}

fn integrated_uv_result_change_sigma(results: &IntegratedUvResults) -> f64 {
    let delta_re = results.integrated.result.re.0 - results.no_integrated.result.re.0;
    let delta_im = results.integrated.result.im.0 - results.no_integrated.result.im.0;
    let delta_norm = delta_re.hypot(delta_im);
    let combined_error_norm = results
        .integrated
        .error
        .re
        .0
        .hypot(results.no_integrated.error.re.0)
        .hypot(
            results
                .integrated
                .error
                .im
                .0
                .hypot(results.no_integrated.error.im.0),
        );

    if combined_error_norm == 0.0 {
        if delta_norm == 0.0 {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        delta_norm / combined_error_norm
    }
}

fn integrated_uv_meaningfully_changes_result(
    results: &IntegratedUvResults,
    min_change_sigma: f64,
) -> bool {
    integrated_uv_result_change_sigma(results) >= min_change_sigma
}

fn integrated_uv_profile_passes(
    cli: &mut CLIState,
    process: &str,
    integrand_name: &str,
) -> Result<bool> {
    let res = Profile::UltraViolet(UltraVioletProfile {
        process: Some(ProcessRef::Unqualified(process.to_string())),
        integrand_name: Some(integrand_name.to_string()),
        min_scale_exponent: 4.0,
        max_scale_exponent: 5.0,
        n_points: 100,
        per_orientation: true,
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?;
    let uv = res.unwrap_uv();
    Ok(uv.pass_fail(-0.9).failed == 0)
}

fn integrated_result_is_muv_invariant(
    cli: &mut CLIState,
    test_name: &str,
    integrated_process: &str,
    integrand_name: &str,
    baseline: &IntegralEstimate,
    original_m_uv: f64,
    alternate_m_uv: f64,
) -> Result<bool> {
    cli.run_command(&format!("set process kv general.m_uv={alternate_m_uv}"))?;
    set_fast_deterministic_integrator(cli)?;
    let integration_result = Integrate {
        process: vec![ProcessRef::Unqualified(integrated_process.to_string())],
        integrand_name: vec![integrand_name.to_string()],
        workspace_path: Some(get_tests_workspace_path().join(format!(
            "{test_name}/integration_workspace_{integrand_name}_muv_shifted"
        ))),
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let integrated_key = format!("{integrated_process}@{integrand_name}");
    let integrated_slot = integration_result
        .slot(&integrated_key)
        .expect("integrated shifted-m_uv slot should exist");
    let passed = integrated_slot
        .integral
        .is_compatible_with_result(baseline, 2);

    cli.run_command(&format!("set process kv general.m_uv={original_m_uv}"))?;

    Ok(passed)
}

fn run_integrated_uv_case(case: &IntegratedUvCase<'_>) -> IntegratedUvCaseResult {
    let mut outcome = IntegratedUvCaseResult {
        graph: case.test_name.to_string(),
        uv_profile_passed: false,
        muv_invariance_passed: false,
        target_passed: case.targets.as_ref().map(|_| false),
        target_threshold: case.targets.as_ref().map(|_| "2sigma"),
        ct_change_passed: case.min_change_sigma.map(|_| false),
        ct_change_sigma: None,
        ct_change_threshold: case.min_change_sigma,
        error: None,
    };

    let cli_result = get_test_cli(
        Some(format!("{}.toml", case.run_card).into()),
        get_tests_workspace_path().join(case.test_name),
        None,
        true,
    );
    let mut cli = match cli_result {
        Ok(cli) => cli,
        Err(err) => {
            outcome.error = Some(format!("setup failed: {err:#}"));
            return outcome;
        }
    };

    let case_result = (|| -> Result<()> {
        cli.run_command("run generate")?;

        let integrated_profile =
            integrated_uv_profile_passes(&mut cli, case.integrated_process, case.integrand_name)?;
        let no_integrated_profile = integrated_uv_profile_passes(
            &mut cli,
            case.no_integrated_process,
            case.integrand_name,
        )?;
        outcome.uv_profile_passed = integrated_profile && no_integrated_profile;

        let baseline = run_integrated_uv_integration(
            &mut cli,
            case.test_name,
            case.no_integrated_process,
            case.integrated_process,
            case.integrand_name,
        )?;

        if let Some(min_change_sigma) = case.min_change_sigma {
            let change_sigma = integrated_uv_result_change_sigma(&baseline);
            outcome.ct_change_sigma = Some(change_sigma);
            outcome.ct_change_passed = Some(integrated_uv_meaningfully_changes_result(
                &baseline,
                min_change_sigma,
            ));
        }

        outcome.muv_invariance_passed = integrated_result_is_muv_invariant(
            &mut cli,
            case.test_name,
            case.integrated_process,
            case.integrand_name,
            &baseline.integrated,
            20.0,
            7.0,
        )?;

        if let Some(targets) = &case.targets {
            outcome.target_passed = Some(integrated_uv_targets_pass(
                &baseline,
                targets.no_integrated,
                targets.integrated,
            ));
        }

        Ok(())
    })();

    if let Err(err) = case_result {
        outcome.error = Some(format!("{err:#}"));
    }

    clean_test(&cli.cli_settings.state.folder);
    outcome
}

fn print_integrated_uv_summary(results: &[IntegratedUvCaseResult]) {
    println!(
        "{:<28} {:<6} {:<6} {:<8} {:<12} {:<8} {:<12} error",
        "graph", "uv", "m_uv", "target", "target_thr", "ct_diff", "ct_thr"
    );
    println!("{}", "-".repeat(100));
    for result in results {
        let target = match result.target_passed {
            Some(true) => "pass",
            Some(false) => "FAIL",
            None => "-",
        };
        let target_threshold = result.target_threshold.unwrap_or("-");
        let ct_diff = match result.ct_change_passed {
            Some(true) => "pass",
            Some(false) => "FAIL",
            None => "-",
        };
        let ct_threshold = match (result.ct_change_threshold, result.ct_change_sigma) {
            (Some(threshold), Some(sigma)) => format!(">={threshold:.1}σ ({sigma:.1}σ)"),
            (Some(threshold), None) => format!(">={threshold:.1}σ"),
            (None, _) => "-".to_string(),
        };
        let error = result.error.as_deref().unwrap_or("-");
        println!(
            "{:<28} {:<6} {:<6} {:<8} {:<12} {:<8} {:<12} {}",
            result.graph,
            if result.uv_profile_passed {
                "pass"
            } else {
                "FAIL"
            },
            if result.muv_invariance_passed {
                "pass"
            } else {
                "FAIL"
            },
            target,
            target_threshold,
            ct_diff,
            ct_threshold,
            error,
        );
    }
}

fn run_single_integrated_uv_case(case: &IntegratedUvCase<'_>) {
    let result = run_integrated_uv_case(case);
    print_integrated_uv_summary(std::slice::from_ref(&result));
    assert!(
        result.all_passed(),
        "Integrated UV case failed: {}",
        result.graph
    );
}

#[test]
fn dod1_bubble_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "dod1_bubble",
        test_name: "dod1_bubble",
        no_integrated_process: "bubble_dod1_no_integrated_UV",
        integrated_process: "bubble_dod1",
        integrand_name: "scalar_bubble_below_thres",
        targets: Some(IntegratedUvTargets {
            no_integrated: 7.358320108607984e-3,
            integrated: 1.351_493_784_227_627e-3,
        }),
        min_change_sigma: Some(5.0),
    });
}

#[test]
fn dod0_bubble_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "dod0_bubble",
        test_name: "dod0_bubble",
        no_integrated_process: "bubble_no_integrated_UV",
        integrated_process: "bubble",
        integrand_name: "scalar_bubble_below_thres",
        targets: Some(IntegratedUvTargets {
            no_integrated: 1.471664021721597e-2,
            integrated: 2.7029875684552542e-3,
        }),
        min_change_sigma: Some(5.0),
    });
}

#[test]
fn dod2_bubble_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "dod2_bubble",
        test_name: "dod2_bubble",
        no_integrated_process: "bubble_dod2_no_integrated_UV",
        integrated_process: "bubble_dod2",
        integrand_name: "scalar_bubble_below_thres",
        targets: Some(IntegratedUvTargets {
            no_integrated: -0.5940830828502411,
            integrated: 1.2143596454658382e-2,
        }),
        min_change_sigma: Some(5.0),
    });
}

mod slow {
    #[test]
    fn ad_ad_with_gluon_correction_uv() {
        run_single_integrated_uv_case(&IntegratedUvCase {
            run_card: "ad_ad_with_gluon_correction",
            test_name: "ad_ad_with_gluon_correction",
            no_integrated_process: "adad_no_integrated_UV",
            integrated_process: "adad",
            integrand_name: "adad_gluon",
            targets: None,
            min_change_sigma: Some(5.0),
        });
    }

    #[test]
    fn epem_a_bbx_amp_uv() {
        run_single_integrated_uv_case(&IntegratedUvCase {
            run_card: "epem_a_bbx_amp",
            test_name: "epem_a_bbx_amp",
            no_integrated_process: "epem_a_bbx_no_integrated_UV",
            integrated_process: "epem_a_bbx",
            integrand_name: "epem_a_bbx",
            targets: None,
            min_change_sigma: None,
        });
    }

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
