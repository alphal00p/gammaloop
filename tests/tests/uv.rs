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
    process: &'a str,
    integrand_name: &'a str,
    original_m_uv: f64,
    shifted_m_uv: f64,
    original_mu_r: f64,
    shifted_mu_r: f64,
    targets: Option<IntegratedUvTargets>,
    min_change_sigma: Option<f64>,
    min_mu_r_change_sigma: Option<f64>,
}

struct IntegratedUvResults {
    no_integrated: Option<IntegralEstimate>,
    integrated: IntegralEstimate,
    integrated_muv_shifted: Option<IntegralEstimate>,
    integrated_mur_shifted: Option<IntegralEstimate>,
}

struct IntegratedUvCaseResult {
    graph: String,
    uv_profile_passed: bool,
    muv_invariance_passed: bool,
    mur_dependence_passed: Option<bool>,
    target_passed: Option<bool>,
    target_threshold: Option<&'static str>,
    ct_change_passed: Option<bool>,
    ct_change_sigma: Option<f64>,
    ct_change_threshold: Option<f64>,
    mur_change_sigma: Option<f64>,
    mur_change_threshold: Option<f64>,
    error: Option<String>,
}

impl IntegratedUvCaseResult {
    fn all_passed(&self) -> bool {
        self.error.is_none()
            && self.uv_profile_passed
            && self.muv_invariance_passed
            && self.mur_dependence_passed.unwrap_or(true)
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

fn no_integrated_process_name(process: &str) -> String {
    format!("{process}_no_integrated_UV")
}

fn shifted_muv_integrand_name(integrand_name: &str) -> String {
    format!("{integrand_name}_muv_shifted")
}

fn shifted_mur_integrand_name(integrand_name: &str) -> String {
    format!("{integrand_name}_mur_shifted")
}

fn add_integrated_uv_scale_variants(cli: &mut CLIState, case: &IntegratedUvCase<'_>) -> Result<()> {
    if (case.original_m_uv - case.shifted_m_uv).abs() > f64::EPSILON {
        let shifted_name = shifted_muv_integrand_name(case.integrand_name);
        cli.run_command(&format!(
            "duplicate integrand -p {} -i {} --output_process_name {} --output_integrand_name {}",
            case.process, case.integrand_name, case.process, shifted_name
        ))?;
        cli.run_command(&format!(
            "set process -p {} -i {} kv general.m_uv={}",
            case.process, shifted_name, case.shifted_m_uv
        ))?;
    }

    if case.min_mu_r_change_sigma.is_some()
        && (case.original_mu_r - case.shifted_mu_r).abs() > f64::EPSILON
    {
        let shifted_name = shifted_mur_integrand_name(case.integrand_name);
        cli.run_command(&format!(
            "duplicate integrand -p {} -i {} --output_process_name {} --output_integrand_name {}",
            case.process, case.integrand_name, case.process, shifted_name
        ))?;
        cli.run_command(&format!(
            "set process -p {} -i {} kv general.mu_r={}",
            case.process, shifted_name, case.shifted_mu_r
        ))?;
    }

    Ok(())
}

fn process_integrand_exists(cli: &mut CLIState, process: &str, integrand_name: &str) -> bool {
    cli.state
        .find_integrand_ref(
            Some(&ProcessRef::Unqualified(process.to_string())),
            Some(&integrand_name.to_string()),
        )
        .is_ok()
}

fn run_integrated_uv_integration(
    cli: &mut CLIState,
    case: &IntegratedUvCase<'_>,
) -> Result<IntegratedUvResults> {
    set_fast_deterministic_integrator(cli)?;
    let no_integrated_process = no_integrated_process_name(case.process);
    let has_no_integrated =
        process_integrand_exists(cli, &no_integrated_process, case.integrand_name);
    let mut processes = vec![ProcessRef::Unqualified(case.process.to_string())];
    let mut integrand_names = vec![case.integrand_name.to_string()];

    if has_no_integrated {
        processes.push(ProcessRef::Unqualified(no_integrated_process.clone()));
        integrand_names.push(case.integrand_name.to_string());
    }

    let muv_shifted_name = ((case.original_m_uv - case.shifted_m_uv).abs() > f64::EPSILON)
        .then(|| shifted_muv_integrand_name(case.integrand_name));
    if let Some(name) = &muv_shifted_name {
        processes.push(ProcessRef::Unqualified(case.process.to_string()));
        integrand_names.push(name.clone());
    }

    let mur_shifted_name = (case.min_mu_r_change_sigma.is_some()
        && (case.original_mu_r - case.shifted_mu_r).abs() > f64::EPSILON)
        .then(|| shifted_mur_integrand_name(case.integrand_name));
    if let Some(name) = &mur_shifted_name {
        processes.push(ProcessRef::Unqualified(case.process.to_string()));
        integrand_names.push(name.clone());
    }

    let integration_result = Integrate {
        process: processes,
        integrand_name: integrand_names,
        workspace_path: Some(get_tests_workspace_path().join(format!(
            "{}/integration_workspace_{}",
            case.test_name, case.integrand_name
        ))),
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let integrated_key = format!("{}@{}", case.process, case.integrand_name);
    let integrated_slot = integration_result
        .slot(&integrated_key)
        .expect("integrated slot should exist");

    Ok(IntegratedUvResults {
        no_integrated: has_no_integrated.then(|| {
            integration_result
                .slot(&format!(
                    "{}@{}",
                    no_integrated_process, case.integrand_name
                ))
                .expect("no-integrated slot should exist")
                .integral
                .clone()
        }),
        integrated: integrated_slot.integral.clone(),
        integrated_muv_shifted: muv_shifted_name.map(|name| {
            integration_result
                .slot(&format!("{}@{}", case.process, name))
                .expect("integrated shifted-m_uv slot should exist")
                .integral
                .clone()
        }),
        integrated_mur_shifted: mur_shifted_name.map(|name| {
            integration_result
                .slot(&format!("{}@{}", case.process, name))
                .expect("integrated shifted-mu_r slot should exist")
                .integral
                .clone()
        }),
    })
}

fn integrated_uv_targets_pass(
    results: &IntegratedUvResults,
    no_integrated_target: f64,
    integrated_target: f64,
) -> bool {
    results
        .no_integrated
        .as_ref()
        .expect("no-integrated result should exist for target checks")
        .is_compatible_with_target(Complex::new_re(F(no_integrated_target)), 2)
        && results
            .integrated
            .is_compatible_with_target(Complex::new_re(F(integrated_target)), 2)
}

fn integral_estimate_change_sigma(lhs: &IntegralEstimate, rhs: &IntegralEstimate) -> f64 {
    let delta_re = lhs.result.re.0 - rhs.result.re.0;
    let delta_im = lhs.result.im.0 - rhs.result.im.0;
    let delta_norm = delta_re.hypot(delta_im);
    let combined_error_norm = lhs
        .error
        .re
        .0
        .hypot(rhs.error.re.0)
        .hypot(lhs.error.im.0.hypot(rhs.error.im.0));

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

fn integrated_uv_result_change_sigma(results: &IntegratedUvResults) -> f64 {
    integral_estimate_change_sigma(
        &results.integrated,
        results
            .no_integrated
            .as_ref()
            .expect("no-integrated result should exist for ct-diff checks"),
    )
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
    baseline: &IntegralEstimate,
    shifted: &IntegralEstimate,
) -> Result<bool> {
    Ok(shifted.is_compatible_with_result(baseline, 2))
}

fn integrated_result_mu_r_change_sigma(
    baseline: &IntegralEstimate,
    shifted: &IntegralEstimate,
) -> Result<f64> {
    Ok(integral_estimate_change_sigma(shifted, baseline))
}

fn run_integrated_uv_case(case: &IntegratedUvCase<'_>) -> IntegratedUvCaseResult {
    let mut outcome = IntegratedUvCaseResult {
        graph: case.test_name.to_string(),
        uv_profile_passed: false,
        muv_invariance_passed: false,
        mur_dependence_passed: case.min_mu_r_change_sigma.map(|_| false),
        target_passed: case.targets.as_ref().map(|_| false),
        target_threshold: case.targets.as_ref().map(|_| "2sigma"),
        ct_change_passed: case.min_change_sigma.map(|_| false),
        ct_change_sigma: None,
        ct_change_threshold: case.min_change_sigma,
        mur_change_sigma: None,
        mur_change_threshold: case.min_mu_r_change_sigma,
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
        add_integrated_uv_scale_variants(&mut cli, case)?;

        let no_integrated_process = no_integrated_process_name(case.process);
        let has_no_integrated =
            process_integrand_exists(&mut cli, &no_integrated_process, case.integrand_name);
        outcome.uv_profile_passed =
            integrated_uv_profile_passes(&mut cli, case.process, case.integrand_name)?
                && (!has_no_integrated
                    || integrated_uv_profile_passes(
                        &mut cli,
                        &no_integrated_process,
                        case.integrand_name,
                    )?);

        let baseline = run_integrated_uv_integration(&mut cli, case)?;

        if let Some(min_change_sigma) = case.min_change_sigma {
            if baseline.no_integrated.is_none() {
                return Err(eyre::eyre!(
                    "ct-diff check requested for {} but no {} process exists",
                    case.process,
                    no_integrated_process
                ));
            }
            let change_sigma = integrated_uv_result_change_sigma(&baseline);
            outcome.ct_change_sigma = Some(change_sigma);
            outcome.ct_change_passed = Some(integrated_uv_meaningfully_changes_result(
                &baseline,
                min_change_sigma,
            ));
        }

        outcome.muv_invariance_passed = integrated_result_is_muv_invariant(
            &baseline.integrated,
            baseline
                .integrated_muv_shifted
                .as_ref()
                .expect("shifted m_uv result should exist"),
        )?;

        if let Some(min_mu_r_change_sigma) = case.min_mu_r_change_sigma {
            let change_sigma = integrated_result_mu_r_change_sigma(
                &baseline.integrated,
                baseline
                    .integrated_mur_shifted
                    .as_ref()
                    .expect("shifted mu_r result should exist"),
            )?;
            outcome.mur_change_sigma = Some(change_sigma);
            outcome.mur_dependence_passed = Some(change_sigma >= min_mu_r_change_sigma);
        }

        if let Some(targets) = &case.targets {
            if baseline.no_integrated.is_none() {
                return Err(eyre::eyre!(
                    "target check requested for {} but no {} process exists",
                    case.process,
                    no_integrated_process
                ));
            }
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
        "{:<28} {:<6} {:<6} {:<8} {:<12} {:<8} {:<12} {:<8} {:<12} error",
        "graph", "uv", "m_uv", "mu_r", "mu_r_thr", "target", "target_thr", "ct_diff", "ct_thr"
    );
    println!("{}", "-".repeat(114));
    for result in results {
        let mu_r = match result.mur_dependence_passed {
            Some(true) => "pass",
            Some(false) => "FAIL",
            None => "-",
        };
        let mu_r_threshold = match (result.mur_change_threshold, result.mur_change_sigma) {
            (Some(threshold), Some(sigma)) => format!(">={threshold:.1}σ ({sigma:.1}σ)"),
            (Some(threshold), None) => format!(">={threshold:.1}σ"),
            (None, _) => "-".to_string(),
        };
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
            "{:<28} {:<6} {:<6} {:<8} {:<12} {:<8} {:<12} {:<8} {:<12} {}",
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
            mu_r,
            mu_r_threshold,
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
fn aa_aa_gl00_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/aa_aa_GL00",
        test_name: "aa_aa_GL00",
        process: "aa_aa",
        integrand_name: "2L",
        original_m_uv: 91.188,
        shifted_m_uv: 364.752,
        original_mu_r: 91.188,
        shifted_mu_r: 364.752,
        targets: None,
        min_change_sigma: Some(5.0),
        min_mu_r_change_sigma: Some(5.0),
    });
}
#[test]
fn epem_ttxh_gl00_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/epem_ttxh_GL00",
        test_name: "epem_ttxh_gl00",
        process: "epem_a_tth",
        integrand_name: "NLO",
        original_m_uv: 20.0,
        shifted_m_uv: 7.0,
        original_mu_r: 3.0,
        shifted_mu_r: 9.0,
        targets: None,
        min_change_sigma: Some(5.0),
        min_mu_r_change_sigma: Some(5.0),
    });
}

#[test]
fn dod0_bubble_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/dod0_bubble",
        test_name: "dod0_bubble",
        process: "bubble",
        integrand_name: "scalar_bubble_below_thres",
        original_m_uv: 20.0,
        shifted_m_uv: 7.0,
        original_mu_r: 3.0,
        shifted_mu_r: 9.0,
        targets: Some(IntegratedUvTargets {
            no_integrated: 1.4693e-2,
            integrated: 2.684e-3,
        }),
        min_change_sigma: Some(5.0),
        min_mu_r_change_sigma: Some(5.0),
    });
}

#[test]
fn dod1_bubble_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/dod1_bubble",
        test_name: "dod1_bubble",
        process: "bubble_dod1",
        integrand_name: "scalar_bubble_below_thres",
        original_m_uv: 20.0,
        shifted_m_uv: 7.0,
        original_mu_r: 3.0,
        shifted_mu_r: 9.0,
        targets: Some(IntegratedUvTargets {
            no_integrated: 7.358320108607984e-3,
            integrated: 1.351_493_784_227_627e-3,
        }),
        min_change_sigma: Some(5.0),
        min_mu_r_change_sigma: Some(5.0),
    });
}

#[test]
fn dod2_bubble_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/dod2_bubble",
        test_name: "dod2_bubble",
        process: "bubble_dod2",
        integrand_name: "scalar_bubble_below_thres",
        original_m_uv: 20.0,
        shifted_m_uv: 7.0,
        original_mu_r: 3.0,
        shifted_mu_r: 19.0,
        targets: Some(IntegratedUvTargets {
            no_integrated: -0.5940830828502411,
            integrated: 1.2143596454658382e-2,
        }),
        min_change_sigma: Some(5.0),
        min_mu_r_change_sigma: Some(5.0),
    });
}

#[test]
fn epem_a_bbx_amp_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/epem_a_bbx_amp",
        test_name: "epem_a_bbx_amp",
        process: "epem_a_bbx",
        integrand_name: "epem_a_bbx",
        original_m_uv: 20.0,
        shifted_m_uv: 7.0,
        original_mu_r: 20.0,
        shifted_mu_r: 9.0,
        targets: None,
        min_change_sigma: Some(5.0),
        min_mu_r_change_sigma: Some(5.0),
    });
}

mod slow {
    #[test]
    fn ad_ad_with_gluon_correction_uv() {
        run_single_integrated_uv_case(&IntegratedUvCase {
            run_card: "uv/ad_ad_with_gluon_correction",
            test_name: "ad_ad_with_gluon_correction",
            process: "adad",
            integrand_name: "adad_gluon",
            original_m_uv: 20.0,
            shifted_m_uv: 7.0,
            original_mu_r: 20.0,
            shifted_mu_r: 9.0,
            targets: None,
            min_change_sigma: Some(5.0),
            min_mu_r_change_sigma: Some(5.0),
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
