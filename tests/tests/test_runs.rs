#![allow(dead_code)]
#![allow(unused_variables)]

use color_eyre::Result;

use colored::{ColoredString, Colorize};
use eyre::Ok;
use gammaloop_api::{
    commands::{
        Profile, Renormalize, inspect::Inspect, integrate::Integrate, profile::UltraVioletProfile,
    },
    state::{ProcessRef, RunHistory, SyncSettings},
};

use gammalooprs::{
    model::UFOSymbol,
    settings::runtime::{IntegralEstimate, IntegrationResult as RuntimeIntegrationResult},
    utils::{self, ApproxEq, F},
};
use insta::assert_snapshot;
use itertools::Itertools;
use momtrop::assert_approx_eq;
use serde_json::Value as JsonValue;
use serial_test::serial;
use spenso::algebra::complex::Complex;
use std::path::{Path, PathBuf};
use std::{env, fmt::Write, time::Duration};
use symbolica::symbol;
use tabled::{Table, Tabled, settings::Style};
use tracing::info;

use gammaloop_integration_tests::{
    clean_test, get_test_cli, get_tests_workspace_path, new_cli_for_test, run_commands,
    setup_sm_differential_lu_cli,
};

fn single_slot_integral(result: &RuntimeIntegrationResult) -> &IntegralEstimate {
    &result
        .single_slot()
        .expect("expected a single integration-result slot")
        .integral
}

fn default_integrate_for(name: &str) -> Integrate {
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        n_cores: Some(1),
        workspace_path: Some(
            get_tests_workspace_path().join(format!("{name}/integration_workspace")),
        ),
        target: vec![],
        restart: true,
        ..Default::default()
    }
}

fn selected_slot_workspace(
    cli: &gammaloop_integration_tests::CLIState,
    workspace: &Path,
    process: Option<ProcessRef>,
    integrand_name: Option<&str>,
) -> Result<PathBuf> {
    let integrand_name = integrand_name.map(str::to_string);
    let (process_id, integrand_name) = cli
        .state
        .find_integrand_ref(process.as_ref(), integrand_name.as_ref())?;
    let process_name = cli.state.process_list.processes[process_id]
        .definition
        .folder_name
        .clone();
    Ok(workspace
        .join("integrands")
        .join(format!("{process_name}@{integrand_name}")))
}

const SCALAR_TRIANGLE_EXTERNALS: &str = r#"set process -p triangle -i scalar_tri string '
[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [1.0, 0.0, 0.0, 0.0],
    [0.5, 0.0, 0.0, -0.5],
    "dependent"
]
helicities = [0, 0, 0]
'"#;

const SCALAR_BOX_ABOVE_EXTERNALS: &str = r#"set process -p box -i scalar_box string '
[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [3.0, 0.0, 0.0, 3.0],
    [3.0, 0.0, 0.0, -3.0],
    [3.0, 0.0, 3.0, 0.0],
    "dependent"
]
helicities = [0, 0, 0, 0]
'"#;

const SCALAR_BOX_BELOW_EXTERNALS: &str = r#"set process -p box -i scalar_box string '
[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [0.3, 0.0, 0.0, 0.3],
    [0.3, 0.0, 0.0, -0.3],
    [0.3, 0.3, 0.0, 0.0],
    "dependent"
]
helicities = [0, 0, 0, 0]
'"#;

const SCALAR_BOX_COPY_ABOVE_EXTERNALS: &str = r#"set process -p box_copy -i scalar_box_copy string '
[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [3.0, 0.0, 0.0, 3.0],
    [3.0, 0.0, 0.0, -3.0],
    [3.0, 0.0, 3.0, 0.0],
    "dependent"
]
helicities = [0, 0, 0, 0]
'"#;

const SCALAR_BOX_COPY_BELOW_EXTERNALS: &str = r#"set process -p box_copy -i scalar_box_copy string '
[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [0.3, 0.0, 0.0, 0.3],
    [0.3, 0.0, 0.0, -0.3],
    [0.3, 0.3, 0.0, 0.0],
    "dependent"
]
helicities = [0, 0, 0, 0]
'"#;

fn setup_scalar_topologies_cli(test_name: &str) -> Result<gammaloop_integration_tests::CLIState> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join(test_name),
        Some(test_name.to_string()),
        true,
    )?;

    run_commands(
        &mut cli,
        &[
            "import model scalars-default.json",
            "reset processes",
            "generate amp scalar_1 > scalar_0 scalar_0 [{1}] --allowed-vertex-interactions V_3_SCALAR_022 V_3_SCALAR_122 -p triangle -i scalar_tri",
            "generate amp scalar_0 scalar_0 > scalar_0 scalar_0 [{1}] --allowed-vertex-interactions V_3_SCALAR_022 -p box -i scalar_box --select-graphs GL0",
            "generate",
            "set model mass_scalar_2=2.0",
            "set model mass_scalar_1=1.0",
            SCALAR_TRIANGLE_EXTERNALS,
            SCALAR_BOX_ABOVE_EXTERNALS,
        ],
    )?;

    Ok(cli)
}

fn scalar_topology_integrate_command(
    test_name: &str,
    workspace_name: &str,
    selections: &[(&str, &str)],
    targets: &[(&str, Complex<F<f64>>)],
) -> Integrate {
    let workspace = get_tests_workspace_path()
        .join(test_name)
        .join(workspace_name);
    Integrate {
        process: selections
            .iter()
            .map(|(process, _)| ProcessRef::Unqualified((*process).to_string()))
            .collect(),
        integrand_name: selections
            .iter()
            .map(|(_, integrand)| (*integrand).to_string())
            .collect(),
        workspace_path: Some(workspace),
        n_cores: Some(1),
        target: targets
            .iter()
            .map(|(key, target)| format!("{key}={},{}", target.re.0, target.im.0))
            .collect(),
        restart: true,
        ..Default::default()
    }
}

fn load_integration_result(path: &Path) -> Result<RuntimeIntegrationResult> {
    Ok(serde_json::from_str(&std::fs::read_to_string(path)?)?)
}

fn normalized_integration_result_json(result: &RuntimeIntegrationResult) -> JsonValue {
    let mut value = serde_json::to_value(result).expect("integration result must serialize");
    let Some(slots) = value.get_mut("slots").and_then(JsonValue::as_array_mut) else {
        return value;
    };
    for slot in slots {
        if let Some(statistics) = slot
            .get_mut("integration_statistics")
            .and_then(JsonValue::as_object_mut)
        {
            statistics.remove("average_total_time_seconds");
            statistics.remove("average_parameterization_time_seconds");
            statistics.remove("average_integrand_time_seconds");
            statistics.remove("average_evaluator_time_seconds");
            statistics.remove("average_observable_time_seconds");
            statistics.remove("average_integrator_time_seconds");
        }
    }
    value
}

fn assert_integration_results_match_ignoring_timings(
    lhs: &RuntimeIntegrationResult,
    rhs: &RuntimeIntegrationResult,
) {
    assert_eq!(
        normalized_integration_result_json(lhs),
        normalized_integration_result_json(rhs),
    );
}

fn complex_distance(lhs: Complex<f64>, rhs: Complex<f64>) -> f64 {
    let delta_re = lhs.re - rhs.re;
    let delta_im = lhs.im - rhs.im;
    (delta_re * delta_re + delta_im * delta_im).sqrt()
}

#[test]
fn oak() -> Result<()> {
    let state = get_test_cli(
        Some("generate_oak_diag.toml".into()),
        get_tests_workspace_path().join("generate_oak_diag"),
        None,
        false,
    )?;
    Ok(())
}

#[test]
fn addbar() -> Result<()> {
    let state = get_test_cli(
        Some("addbar_generate.toml".into()),
        get_tests_workspace_path().join("addbar_generate"),
        None,
        false,
    )?;
    Ok(())
}

#[test]
fn test_z_decay() -> Result<()> {
    let mut state = get_test_cli(
        Some("z_decay_test.toml".into()),
        get_tests_workspace_path().join("z_decay_test"),
        None,
        false,
    )?;

    let result = Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(
            get_tests_workspace_path().join("z_decay_test/integration_workspace/"),
        ),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut state.state, &state.cli_settings)?;

    assert_approx_eq(
        &single_slot_integral(&result).result.re,
        &F(-0.372),
        &F(1e-2),
    );
    Ok(())
}

#[test]
fn epem_ttxh_lo() -> Result<()> {
    let mut cli = get_test_cli(
        Some("epem_ttxh_lo.toml".into()),
        get_tests_workspace_path().join("epem_ttxh_lo"),
        None,
        false,
    )?;

    let integrate_command = Integrate {
        workspace_path: Some(get_tests_workspace_path().join("epem_ttxh_lo/integration_workspace")),
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    };

    let intergal = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        intergal.is_compatible_with_target(
            Complex::new(F(1.0185532594130467e-4), F(-2.7124366612352106e-7)),
            1
        ),
        "Not compatible: {intergal}",
    );
    Ok(())
}

#[test]
fn test_epem_dd_dt() -> Result<()> {
    let mut cli = get_test_cli(
        Some("test_epem_dd_dt.toml".into()),
        get_tests_workspace_path().join("test_epem_dd_dt"),
        None,
        false,
    )?;

    let (_, a) = Inspect {
        process: None,
        integrand_name: Some("default".to_string()),
        momentum_space: false,
        point: vec![
            0.3684805278343727,
            0.20242350926484026,
            0.48242524836619227,
            0.3919994435419629,
            0.14345675651437087,
            0.47075342347077187,
        ],
        ..Default::default()
    }
    .run(&mut cli)?;

    assert_snapshot!(format!("{a:.8e}"),@"(1.3677606162044735e-3+4.451464910288581e-4i)");

    let integrate_command = Integrate {
        workspace_path: Some(
            get_tests_workspace_path().join("test_epem_dd_dt/integration_workspace"),
        ),
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    };

    let intergal = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        intergal.is_compatible_with_target(
            Complex::new(F(1.0185532594130467e-4), F(-2.7124366612352106e-7)),
            1
        ),
        "Not compatible: {intergal}",
    );

    // todo add integration
    Ok(())
}
#[test]
fn test_pentabox_dario() -> Result<()> {
    let state = get_test_cli(
        Some("test_qqx_aaa_pentabox_generate.toml".into()),
        get_tests_workspace_path().join("test_qqx_aaa_pentabox_generate"),
        None,
        false,
    )?;

    // todo add integration
    Ok(())
}

#[test]
fn qqx_aaa_subtracted_nlo_amplitude_test() -> Result<()> {
    let mut state = get_test_cli(
        Some("qqx_aaa_subtracted_nlo_amplitude.toml".into()),
        get_tests_workspace_path().join("qqx_aaa_subtracted_nlo_amplitude"),
        None,
        false,
    )?;

    let (jac, a) = Inspect {
        process: None,
        integrand_name: Some("default".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        ..Default::default()
    }
    .run(&mut state)?;

    assert_snapshot!(format!("{:.8e}",a*jac.unwrap()),@"(6.923469780323106e-4+-6.181987787694587e-4i)");

    Ok(())
}

#[test]
fn trees() -> Result<()> {
    let mut cli = get_test_cli(
        Some("trees/qqx_aaa.toml".into()),
        get_tests_workspace_path().join("qqx_aaa_tree"),
        None,
        false,
    )?;

    let (_, a) = Inspect {
        process: None,
        integrand_name: Some("default".to_string()),
        ..Default::default()
    }
    .run(&mut cli)?;

    assert_snapshot!(format!("{a:.8e}"),@"(1.4727604164105595e-4+-1.1503139369130214e-3i)");

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn photons_1l_integrate() -> Result<()> {
    let mut cli = get_test_cli(
        Some("photons_eu.toml".into()),
        get_tests_workspace_path().join("photons_eu_integrate"),
        None,
        true,
    )?;

    // this can be moved to the run card once we have a set model param command

    cli.state
        .model_parameters
        .insert(UFOSymbol(symbol!("UFO::MT")), Complex::new_re(F(1500.0)));

    cli.state
        .model_parameters
        .insert(UFOSymbol(symbol!("UFO::aEWM1")), Complex::new_re(F(128.93)));

    cli.state
        .model
        .apply_param_card(&cli.state.model_parameters)?;

    let target = Complex::new(F(-1.22898408452706e-13), F(-3.94362534040412e-13));

    let integrate = Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        n_cores: Some(1),
        workspace_path: Some(
            get_tests_workspace_path().join("photons_eu_integrate/integration_workspace/"),
        ),
        target: vec![target.re.0.to_string(), target.im.0.to_string()],
        restart: true,
        ..Default::default()
    };

    let integration_result = integrate.run(&mut cli.state, &cli.cli_settings)?;

    let integration_result = single_slot_integral(&integration_result);
    assert_approx_eq(&integration_result.result.re, &target.re, &F(10e-2));
    assert_approx_eq(&integration_result.result.im, &target.im, &F(10e-2));

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn photons_1l_inspect() -> Result<()> {
    let mut cli = get_test_cli(
        Some("photons_eu.toml".into()),
        get_tests_workspace_path().join("photons_eu_inspect"),
        None,
        false,
    )?;

    // this can be moved to the run card once we have a set model param command

    cli.state
        .model_parameters
        .insert(UFOSymbol(symbol!("UFO::MT")), Complex::new_re(F(1500.0)));

    cli.state
        .model_parameters
        .insert(UFOSymbol(symbol!("UFO::aEWM1")), Complex::new_re(F(128.93)));

    cli.state
        .model
        .apply_param_card(&cli.state.model_parameters)?;

    let (_, inspect) = Inspect {
        process: None,
        integrand_name: Some("default".to_string()),
        point: vec![0.123, 0.3242, 0.4233],
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;

    println!("Inspect result: {inspect:.16e}");

    // The old test at a very bad value of e_cm, so I created a new value using the example card in the old main
    assert_snapshot!(format!("{inspect:.8e}"),@"(-4.236544183136419e-12+-3.710728958614228e-12i)");

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

#[test]
fn photons_phys_1l_inspect() -> Result<()> {
    let mut cli = get_test_cli(
        Some("generate_threshold_1L_6photons.toml".into()),
        get_tests_workspace_path().join("photons_phys_1l_inspect"),
        None,
        true,
    )?;
    let (_, inspect) = Inspect {
        process: None,
        integrand_name: Some("default".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;

    let target = Complex::new(-2.827365545920272e-10, -5.127347264133554e-10);
    assert_eq!(inspect, target);

    Ok(())
}

#[test]
fn photons_2l_inspect() -> Result<()> {
    let mut cli = get_test_cli(
        Some("photons_eu_2l.toml".into()),
        get_tests_workspace_path().join("photons_eu_2l_inspect"),
        None,
        false,
    )?;

    cli.state
        .model_parameters
        .insert(UFOSymbol(symbol!("UFO::MT")), Complex::new_re(F(1500.0)));

    cli.state
        .model_parameters
        .insert(UFOSymbol(symbol!("UFO::aEWM1")), Complex::new_re(F(128.93)));

    cli.state
        .model
        .apply_param_card(&cli.state.model_parameters)?;

    let (_, inspect) = Inspect {
        process: None,
        integrand_name: Some("default".to_string()),
        point: vec![0.123, 0.3242, 0.4233, 0.523, 0.314, 0.125],
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;

    println!("Inspect result: {inspect:.16e}");

    // wrong result
    assert_snapshot!(format!("{inspect:.8e}"));

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

#[test]
fn photonic_amplitudes() -> Result<()> {
    let mut test_failed = false;
    #[derive(Tabled)]
    struct ComparisonReport {
        gammaloop: ColoredString,
        reference: ColoredString,
    }

    impl ComparisonReport {
        fn display_nested(&self) -> String {
            let mut table = Table::new([self]);
            table
                .with(
                    Style::modern()
                        .remove_horizontals()
                        .remove_verticals()
                        .remove_frame(),
                )
                .with(tabled::settings::Rotate::Left);
            table.to_string()
        }

        fn empty() -> Self {
            Self {
                gammaloop: "N/A".red(),
                reference: "N/A".red(),
            }
        }

        fn only_gloop(gammaloop: ColoredString) -> Self {
            Self {
                gammaloop: gammaloop.blue(),
                reference: "N/A".red(),
            }
        }
    }

    #[derive(Tabled)]
    struct TestReportEntry {
        amplitude: String,
        #[tabled(display = "ComparisonReport::display_nested")]
        generation_time: ComparisonReport,
        #[tabled(display = "ComparisonReport::display_nested")]
        inspect: ComparisonReport,
        #[tabled(display = "ComparisonReport::display_nested")]
        integrated: ComparisonReport,
        #[tabled(display = "ComparisonReport::display_nested")]
        nvar: ComparisonReport,
        #[tabled(display = "ComparisonReport::display_nested")]
        sample_time: ComparisonReport,
    }

    impl TestReportEntry {
        fn default_with_name(name: String) -> Self {
            Self {
                amplitude: name,
                generation_time: ComparisonReport::empty(),
                inspect: ComparisonReport::empty(),
                integrated: ComparisonReport::empty(),
                nvar: ComparisonReport::empty(),
                sample_time: ComparisonReport::empty(),
            }
        }
    }

    struct BenchMarkData {
        run_card: String,
        state_path: String,
        amplitude: String,
        generation_time: Option<Duration>,
        inspect_point: Vec<f64>,
        inspect_target: Option<Complex<f64>>,
        integrated_target: Option<Complex<F<f64>>>,
        nvar_bench: Option<Complex<F<f64>>>,
        sample_time: Option<Duration>,
    }

    fn run_photonic_test(
        bench_data: BenchMarkData,
        test_failed: &mut bool,
    ) -> Result<TestReportEntry> {
        let mut cli = get_test_cli(
            Some(bench_data.run_card.clone().into()),
            bench_data.state_path.clone(),
            None,
            true,
        )?;

        let before_generation = std::time::Instant::now();
        cli.run_command("generate")?;
        let generation_time = before_generation.elapsed();

        let generation_time_comparison = if let Some(target_time) = bench_data.generation_time {
            let generation_time_string = if generation_time <= target_time {
                format!("{}s", generation_time.as_secs()).green()
            } else if generation_time <= target_time * 10 {
                format!("{}s", generation_time.as_secs()).yellow()
            } else {
                format!("{}s", generation_time.as_secs()).red()
            };

            ComparisonReport {
                gammaloop: generation_time_string,
                reference: format!("{}s", target_time.as_secs()).blue(),
            }
        } else {
            ComparisonReport::only_gloop(format!("{}s", generation_time.as_secs()).blue())
        };

        let (_, inspect_result) = Inspect {
            process: None,
            integrand_name: Some("default".to_string()),
            point: bench_data.inspect_point,
            momentum_space: false,
            ..Default::default()
        }
        .run(&mut cli)?;

        let inspect_result_comparison = if let Some(inspect_target) = bench_data.inspect_target {
            let inspect_result_f = Complex::new(F(inspect_result.re), F(inspect_result.im));
            let inspect_target_f = Complex::new(F(inspect_target.re), F(inspect_target.im));
            let inspect_result_string =
                if inspect_target_f.approx_eq(&inspect_result_f, &F(1.0e-13)) {
                    format!("{:.16e}", inspect_result).green()
                } else {
                    *test_failed = true;
                    format!("{:.16e}", inspect_result).red()
                };

            ComparisonReport {
                gammaloop: inspect_result_string,
                reference: format!("{:.16e}", inspect_target).blue(),
            }
        } else {
            ComparisonReport::only_gloop(format!("{:.16e}", inspect_result).blue())
        };

        let before_integration = std::time::Instant::now();
        let integrated_result = Integrate {
            process: vec![],
            integrand_name: vec!["default".to_string()],
            workspace_path: None,
            n_cores: Some(1),
            target: bench_data
                .integrated_target
                .map(|t| vec![t.re.0.to_string(), t.im.0.to_string()])
                .unwrap_or_default(),
            restart: true,
            ..Default::default()
        }
        .run(&mut cli.state, &cli.cli_settings)?;
        let integrated_result = single_slot_integral(&integrated_result);
        let integration_time = before_integration.elapsed();
        let sample_time = integration_time / (integrated_result.neval as u32);

        let real_part_string =
            utils::format_uncertainty(integrated_result.result.re, integrated_result.error.re);
        let imag_part_string =
            utils::format_uncertainty(integrated_result.result.im, integrated_result.error.im);

        let integrated_result_comparison =
            if let Some(integrated_target) = bench_data.integrated_target {
                let integrated_result_string =
                    if integrated_result.is_compatible_with_target(integrated_target, 3) {
                        format!("{} + {}i", real_part_string, imag_part_string).green()
                    } else {
                        *test_failed = true;
                        format!("{} + {}i", real_part_string, imag_part_string).red()
                    };

                ComparisonReport {
                    gammaloop: integrated_result_string,
                    reference: format!(
                        "{:.8e} + {:.8e}i",
                        integrated_target.re.0, integrated_target.im.0
                    )
                    .blue(),
                }
            } else {
                ComparisonReport::only_gloop(
                    format!("{} +  {}i", real_part_string, imag_part_string).blue(),
                )
            };

        let n_var_re = (integrated_result.error.re.0
            * integrated_result.error.re.0
            * integrated_result.neval as f64)
            / (integrated_result.result.re.0 * integrated_result.result.re.0);
        let n_var_im = (integrated_result.error.im.0
            * integrated_result.error.im.0
            * integrated_result.neval as f64)
            / (integrated_result.result.im.0 * integrated_result.result.im.0);

        let nvar_comparison = if let Some(nvar_target) = bench_data.nvar_bench {
            let nvar_string =
                if n_var_re < nvar_target.re.0 * 2.0 && n_var_im < nvar_target.im.0 * 2.0 {
                    format!("{:.2e} + {:.2e}", n_var_re, n_var_im).green()
                } else if n_var_re < nvar_target.re.0 * 4.0 && n_var_im < nvar_target.im.0 * 4.0 {
                    format!("{:.2e} + {:.2e}", n_var_re, n_var_im).yellow()
                } else {
                    *test_failed = true;
                    format!("{:.2e} + {:.2e}", n_var_re, n_var_im).red()
                };
            ComparisonReport {
                gammaloop: nvar_string,
                reference: format!("{:.2e} + {:.2e}", nvar_target.re.0, nvar_target.im.0).blue(),
            }
        } else {
            ComparisonReport::only_gloop(format!("{:.2e} + {:.2e}", n_var_re, n_var_im).blue())
        };

        let sample_time_comparison = if let Some(sample_time_bench) = bench_data.sample_time {
            let sample_time_string = if sample_time <= sample_time_bench {
                format!("{}μs", sample_time.as_micros()).green()
            } else if sample_time <= sample_time_bench * 10 {
                format!("{}μs", sample_time.as_micros()).yellow()
            } else {
                format!("{}μs", sample_time.as_micros()).red()
            };
            ComparisonReport {
                gammaloop: sample_time_string,
                reference: format!("{}μs", sample_time_bench.as_micros()).blue(),
            }
        } else {
            ComparisonReport::only_gloop(format!("{}μs", sample_time.as_micros()).blue())
        };

        Ok(TestReportEntry {
            amplitude: bench_data.amplitude,
            generation_time: generation_time_comparison,
            inspect: inspect_result_comparison,
            integrated: integrated_result_comparison,
            nvar: nvar_comparison,
            sample_time: sample_time_comparison,
        })
    }

    let mut test_reports: Vec<TestReportEntry> = Vec::new();

    let one_loop_eu = BenchMarkData {
        run_card: "photonic_amplitudes/1l_eu.toml".into(),
        state_path: get_tests_workspace_path()
            .join("photonic_amplitudes/1l_eu")
            .to_string_lossy()
            .into_owned(),
        amplitude: "1l_eu".into(),
        generation_time: Some(Duration::from_secs(7)),
        inspect_point: vec![0.123, 0.3242, 0.4233],
        inspect_target: Some(Complex::new(-4.236544183136417e-12, -3.710728958614226e-12)),
        integrated_target: Some(Complex::new(
            F(-1.22898408452706e-13),
            F(-3.94362534040412e-13),
        )),
        nvar_bench: None,
        sample_time: Some(Duration::from_micros(61)),
    };

    let one_loop_phys = BenchMarkData {
        run_card: "photonic_amplitudes/1l_phys.toml".into(),
        state_path: get_tests_workspace_path()
            .join("photonic_amplitudes/1l_phys")
            .to_string_lossy()
            .into_owned(),
        amplitude: "1l_phys".into(),
        inspect_point: vec![0.1, 0.2, 0.3],
        inspect_target: Some(Complex::new(4.660217572648287e-10, -6.496141401696065e-10)),
        integrated_target: Some(Complex::new(
            F(-9.277_595_006_874_547e-11),
            F(-3.683_945_762_498_705_4e-11),
        )),
        generation_time: Some(Duration::from_secs(7)),
        nvar_bench: None,
        sample_time: Some(Duration::from_micros(590)),
    };

    let two_loop_eu = BenchMarkData {
        run_card: "photonic_amplitudes/2l_eu.toml".into(),
        state_path: get_tests_workspace_path()
            .join("photonic_amplitudes/2l_eu")
            .to_string_lossy()
            .into_owned(),
        amplitude: "2l_eu".into(),
        inspect_point: vec![0.123, 0.3242, 0.4233, 0.523, 0.314, 0.125],
        inspect_target: None,
        integrated_target: Some(Complex::new(F(-1.8006e-15), F(-1.54335e-14))),
        generation_time: Some(Duration::from_secs(60)),
        nvar_bench: None,
        sample_time: Some(Duration::from_micros(1160)),
    };

    test_reports.push(run_photonic_test(one_loop_eu, &mut test_failed)?);
    test_reports.push(run_photonic_test(one_loop_phys, &mut test_failed)?);
    test_reports.push(run_photonic_test(two_loop_eu, &mut test_failed)?);
    test_reports.push(TestReportEntry::default_with_name("2l_phys".to_string()));
    test_reports.push(TestReportEntry::default_with_name("3l_eu".to_string()));
    test_reports.push(TestReportEntry::default_with_name("3l_phys".to_string()));

    let mut table = Table::new(test_reports);
    table.with(Style::modern());
    println!("test results: ");
    println!("{}", table);

    assert!(
        !test_failed,
        "Some photonic amplitude tests failed, see table for details",
    );

    Ok(())
}

#[test]
fn test_grouped_subtraction() -> Result<()> {
    let mut cli = get_test_cli(
        Some("test_grouped_subtraction.toml".into()),
        get_tests_workspace_path().join("test_grouped_subtraction"),
        None,
        false,
    )?;

    let int1 = Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(
            get_tests_workspace_path()
                .join("test_grouped_subtraction/integration_workspace_no_group/"),
        ),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    };
    let int2 = Integrate {
        process: vec![ProcessRef::Id(1)],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(
            get_tests_workspace_path()
                .join("test_grouped_subtraction/integration_workspace_group/"),
        ),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    };

    let integration_results_no_group = int1.run(&mut cli.state, &cli.cli_settings)?;
    let integration_results_group = int2.run(&mut cli.state, &cli.cli_settings)?;

    info!("No group result: {:#?}", integration_results_no_group);
    info!("Group result: {:#?}", integration_results_group);

    let integration_results_no_group = single_slot_integral(&integration_results_no_group);
    let integration_results_group = single_slot_integral(&integration_results_group);

    assert_approx_eq(
        &integration_results_group.result.re,
        &integration_results_no_group.result.re,
        &F(1e-1),
    );
    assert_approx_eq(
        &integration_results_group.result.im,
        &integration_results_no_group.result.im,
        &F(1e-1),
    );

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

#[test]
fn v_diag() -> Result<()> {
    let cli = get_test_cli(
        Some("v_diag.toml".into()),
        get_tests_workspace_path().join("v_diag"),
        Some("v_diag".to_string()),
        false,
    )?;
    Ok(())
}

#[test]
fn scalar_bubble() -> Result<()> {
    let mut cli = get_test_cli(
        Some("scalar_bubble.toml".into()),
        get_tests_workspace_path().join("scalar_bubble"),
        Some("scalar_bubble".to_string()),
        false,
    )?;

    let integrate_command = Integrate {
        ..default_integrate_for("scalar_bubble")
    };

    let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
        ..Default::default()
    });

    // from Kaapo: m=1 muv=5 2.03838e-02 m=2 muv=5 	1.16050e-02	 m=3 muv=5 6.46968e-03

    cli.run_command("set model mass_scalar_1=1.0")?;
    let res = profile_cmd
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
    assert_eq!(res.pass_fail(-0.9).failed, 0);
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(2.03838e-02)), 1));

    cli.run_command("set model mass_scalar_1=2.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(1.16050e-02)), 1));

    cli.run_command("set model mass_scalar_1=3.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(6.46968e-03)), 1));
    let renorm_command = Renormalize::default();

    let res = renorm_command.run(&mut cli.state, &cli.cli_settings)?;

    println!("{}", res[0]);

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

#[test]
fn scalar_sunrise() -> Result<()> {
    let mut cli = get_test_cli(
        Some("scalar_sunrise.toml".into()),
        get_tests_workspace_path().join("scalar_sunrise"),
        Some("scalar_sunrise".to_string()),
        false,
    )?;

    let integrate_command = Integrate {
        ..default_integrate_for("scalar_sunrise")
    };
    let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
        use_f128: false,
        max_scale_exponent: 6.0,
        min_scale_exponent: 1.0,
        ..Default::default()
    });
    // from Kaapo: m=1 muv=5 4.37688e-03 m=2 muv=5 	2.48100e-03	 m=3 muv=5 1.07231e-03
    cli.run_command("set model mass_scalar_1=1.0")?;
    let res = profile_cmd
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
    assert_eq!(res.pass_fail(-0.9).failed, 0);
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(-4.37688e-03)), 1) // .result
                                                                                         // .approx_eq(&Complex::new_re(F(-4.37688e-03)), &F(0.01))
    );

    // assert_snapshot!(format!("{integral_no_cache:.3}"),@"-4.359e-3");

    cli.run_command("set model mass_scalar_1=2.0")?;
    let res = profile_cmd
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
    assert_eq!(res.pass_fail(-0.9).failed, 0);
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    // assert_snapshot!(format!("{integral_no_cache:.3}"),@"-2.474e-3");
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(-2.48100e-03)), 1));

    // cli.run_command("set model mass_scalar_1=3.0")?;
    // let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    // assert_snapshot!(format!("{:.8e}",integral_no_cache.result),@"(-4.5184321377520566e-4+0e0i)");

    let renorm_command = Renormalize::default();

    let res = renorm_command.run(&mut cli.state, &cli.cli_settings)?;

    println!("{}", res[0]);
    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

#[test]
fn scalar_mercedes() -> Result<()> {
    let mut cli = get_test_cli(
        Some("scalar_mercedes.toml".into()),
        get_tests_workspace_path().join("scalar_mercedes"),
        Some("scalar_mercedes".to_string()),
        false,
    )?;

    let integrate_command = Integrate {
        ..default_integrate_for("scalar_mercedes")
    };

    let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
        max_scale_exponent: 4.0,
        use_f128: false,
        ..Default::default()
    });

    //5.89551e-06	3.35645e-06	1.87120e-06

    // from Kaapo: m=1 muv=5 5.89551e-06 m=2 muv=5 	3.35645e-06	 m=3 muv=5 1.87120e-06
    cli.run_command("set model mass_scalar_1=1.0")?;
    let res = profile_cmd
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
    assert_eq!(res.pass_fail(-0.9).failed, 0);
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(5.89551e-06)), 1),
        "Not compatible: {integral_no_cache}",
    );

    cli.run_command("set model mass_scalar_1=2.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(3.35645e-06)), 3),
        "Not compatible: {integral_no_cache}",
    );

    cli.run_command("set model mass_scalar_1=3.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(1.87120e-06)), 3),
        "Not compatible: {integral_no_cache}",
    );

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

#[test]
fn scalar_basketball() -> Result<()> {
    let mut cli = get_test_cli(
        Some("scalar_basketball.toml".into()),
        get_tests_workspace_path().join("scalar_basketball"),
        Some("scalar_basketball".to_string()),
        false,
    )?;

    let integrate_command = Integrate {
        ..default_integrate_for("scalar_basketball")
    };

    let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
        use_f128: false,
        max_scale_exponent: 4.,
        min_scale_exponent: 2.0,
        analyse_analytically: false,
        ..Default::default()
    });

    //1.47240e-03	7.15184e-04	2.27485e-04

    // from Kaapo: m=1 muv=5 1.47240e-03 m=2 muv=5 	7.15184e-04	 m=3 muv=5 2.27485e-04
    cli.run_command("set model mass_scalar_1=1.0")?;
    let res = profile_cmd
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
    assert_eq!(res.pass_fail(-0.9).failed, 0);
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(1.47240e-03)), 1),
        "Not compatible: {integral_no_cache}",
    );

    // cli.run_command("set model mass_scalar_1=2.0")?;
    // let res = profile_cmd.run(&mut cli.state, &cli.cli_settings)?;
    // assert_eq!(res.pass_fail(-0.9).failed, 0);
    // let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    // assert!(
    //     integral_no_cache.is_compatible_with_target(Complex::new_re(F(7.15184e-04)), 3),
    //     "Not compatible: {integral_no_cache}",
    // );

    // cli.run_command("set model mass_scalar_1=3.0")?;
    // let res = profile_cmd.run(&mut cli.state, &cli.cli_settings)?;
    // assert_eq!(res.pass_fail(-0.9).failed, 0);
    // let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    // assert!(
    //     integral_no_cache.is_compatible_with_target(Complex::new_re(F(2.27485e-04)), 3),
    //     "Not compatible: {integral_no_cache}",
    // );

    // clean_test(&cli.cli_settings.state.folder);

    Ok(())
}
#[test]
fn scalar_mercedes_with_extra_loop() -> Result<()> {
    let mut cli = get_test_cli(
        Some("scalar_mercedes_with_extra_loop.toml".into()),
        get_tests_workspace_path().join("scalar_mercedes_with_extra_loop"),
        Some("scalar_mercedes_with_extra_loop".to_string()),
        false,
    )?;

    let integrate_command = Integrate {
        ..default_integrate_for("scalar_mercedes_with_extra_loop")
    };

    let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
        max_scale_exponent: 7.,
        min_scale_exponent: 4.0,
        n_points: 15,
        use_f128: false,
        output_file: Some("uv_profile_extra_loop".into()),
        ..Default::default()
    });
    //2.90078e-06	1.59168e-06	6.86001e-07

    // from Kaapo: m=1 muv=5 2.90078e-06 m=2 muv=5 	1.59168e-06	 m=3 muv=5 6.86001e-07
    cli.run_command("set model mass_scalar_1=1.0")?;
    let res = profile_cmd
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
    assert_eq!(res.pass_fail(-0.9).failed, 0);

    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(2.90078e-06)), 1),
        "Not compatible: {integral_no_cache}",
    );

    cli.run_command("set model mass_scalar_1=2.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(1.59168e-06)), 3),
        "Not compatible: {integral_no_cache}",
    );

    cli.run_command("set model mass_scalar_1=3.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(6.86001e-07)), 3),
        "Not compatible: {integral_no_cache}",
    );

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}
#[test]
fn scalar_sunrise_inspect() -> Result<()> {
    symbolica::GLOBAL_SETTINGS
        .initialize_tracing
        .store(false, std::sync::atomic::Ordering::Relaxed);
    let mut cli = get_test_cli(
        Some("scalar_sunrise.toml".into()),
        get_tests_workspace_path().join("scalar_sunrise"),
        Some("scalar_sunrise".to_string()),
        false,
    )?;

    let point = [1., 1., 1., 2., 3., 4.];

    let point = [1., 1., 1., -3., -4., -5.];

    let point = vec![2., 3., 4., 1., 1., 1.];
    let mut ins = Inspect {
        point: point.clone(),
        momentum_space: true,
        discrete_dim: vec![0],
        ..Default::default()
    };

    // from Kaapo: m=1 muv=5 4.37688e-03 m=2 muv=5 	2.48100e-03	 m=3 muv=5 1.07231e-03
    cli.run_command("set model mass_scalar_1=1.0")?;

    fn string_with_prefactor(rs: &[Complex<f64>]) -> String {
        let mut out = String::new();
        let prefactor = -(2. * std::f64::consts::PI).powi(6);
        for r in rs {
            let re = r.re * prefactor;
            let im = r.im * prefactor;
            writeln!(&mut out, "{re:.5e}+i{im:.5e}").unwrap();
        }
        out
    }

    let (jac, rall_1) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, rall_10) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, rall_100) = ins.run(&mut cli.state)?;
    ins.point = point.clone();
    cli.run_command("set process -p 0 -i default kv general.additional_param_values=[1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]")?;
    let (jac, r1_1) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r1_10) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r1_100) = ins.run(&mut cli.state)?;
    ins.point = point.clone();
    cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0]")?;
    let (jac, r2_1) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r2_10) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r2_100) = ins.run(&mut cli.state)?;
    ins.point = point.clone();
    cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0]")?;
    let (jac, r3_1) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r3_10) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r3_100) = ins.run(&mut cli.state)?;
    ins.point = point.clone();
    cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0]")?;
    let (jac, r4_1) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r4_10) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r4_100) = ins.run(&mut cli.state)?;
    ins.point = point.clone();
    cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0]")?;
    let (jac, r5_1) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r5_10) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r5_100) = ins.run(&mut cli.state)?;
    ins.point = point.clone();
    cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0]")?;
    let (jac, r6_1) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r6_10) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r6_100) = ins.run(&mut cli.state)?;

    ins.point = point.clone();
    cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0]")?;
    let (jac, r7_1) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r7_10) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r7_100) = ins.run(&mut cli.state)?;
    ins.point = point.clone();
    cli.run_command("set process -p 0 -i default kv general.additional_param_values=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0]")?;
    let (jac, r8_1) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r8_10) = ins.run(&mut cli.state)?;
    ins.point = ins.point.iter().map(|a| a * 10.).collect_vec();
    let (jac, r8_100) = ins.run(&mut cli.state)?;

    insta::assert_snapshot!(string_with_prefactor(&[r1_1,r2_1,r3_1,r4_1,r5_1,r6_1,r7_1,r8_1,rall_1]),@"
    2.18603e-4+i-0.00000e0
    -4.41097e-5+i-0.00000e0
    -1.54032e-4+i-0.00000e0
    -1.57503e-4+i-0.00000e0
    -7.17631e-5+i-0.00000e0
    4.21936e-5+i-0.00000e0
    1.40322e-4+i-0.00000e0
    8.50437e-5+i-0.00000e0
    5.87544e-5+i-0.00000e0
    ");
    insta::assert_snapshot!(string_with_prefactor(&[r1_10,r2_10,r3_10,r4_10,r5_10,r6_10,r7_10,r8_10,rall_10]),@"
    2.66555e-8+i-0.00000e0
    -1.11736e-8+i-0.00000e0
    -3.96106e-7+i-0.00000e0
    -4.55447e-8+i-0.00000e0
    -2.65802e-8+i-0.00000e0
    1.11735e-8+i-0.00000e0
    3.96096e-7+i-0.00000e0
    4.54492e-8+i-0.00000e0
    -3.03929e-11+i-0.00000e0
    ");
    insta::assert_snapshot!(string_with_prefactor(&[r1_100,r2_100,r3_100,r4_100,r5_100,r6_100,r7_100,r8_100,rall_100]),@"
    2.67150e-12+i-0.00000e0
    -1.13180e-12+i-0.00000e0
    -4.46155e-11+i-0.00000e0
    -4.62050e-12+i-0.00000e0
    -2.67150e-12+i-0.00000e0
    1.13180e-12+i-0.00000e0
    4.46155e-11+i-0.00000e0
    4.62050e-12+i-0.00000e0
    -3.63173e-19+i-0.00000e0
    ");
    // clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

#[test]
fn scalar_box() -> Result<()> {
    let mut cli = get_test_cli(
        Some("scalar_box.toml".into()),
        get_tests_workspace_path().join("scalar_box"),
        Some("scalar_box".to_string()),
        false,
    )?;

    let (_, a) = Inspect {
        process: None,
        integrand_name: Some("default".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        discrete_dim: vec![0],
        ..Default::default()
    }
    .run(&mut cli)?;

    assert_snapshot!(format!("{a:.8e}"),@"(1.3485885914334373e-3+0e0i)");

    cli.run_command("set process -p 0 -i default kv general.enable_cache=false")?;
    let integral_no_cache = Integrate {
        ..default_integrate_for("scalar_box")
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command("set process -p 0 -i default kv general.enable_cache=true")?;
    let integral_with_cache = Integrate {
        ..default_integrate_for("scalar_box")
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    info!("Integral result without caching: {:#?}", integral_no_cache);
    info!(
        "Integral result with caching.  : {:#?}",
        integral_with_cache
    );

    let integral_no_cache = single_slot_integral(&integral_no_cache);
    let integral_with_cache = single_slot_integral(&integral_with_cache);

    assert_approx_eq(
        &integral_no_cache.result.re,
        &integral_with_cache.result.re,
        &F(1e-4),
    );
    assert_approx_eq(
        &integral_no_cache.result.im,
        &integral_with_cache.result.im,
        &F(1e-4),
    );

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

#[test]
fn test_multi_integrand() -> Result<()> {
    let test_name = "test_multi_integrand";
    let mut cli = setup_scalar_topologies_cli(test_name)?;
    let shared_integrator_settings = "kv integrator.n_start=5000 integrator.n_max=10000 integrator.n_increase=5000 integrator.seed=1337";

    cli.run_command(&format!(
        "set process -p triangle -i scalar_tri {shared_integrator_settings}"
    ))?;
    let triangle_baseline = scalar_topology_integrate_command(
        test_name,
        "triangle_baseline",
        &[("triangle", "scalar_tri")],
        &[],
    )
    .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    let box_above_baseline = scalar_topology_integrate_command(
        test_name,
        "box_above_baseline",
        &[("box", "scalar_box")],
        &[],
    )
    .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command(SCALAR_BOX_BELOW_EXTERNALS)?;
    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    let box_below_baseline = scalar_topology_integrate_command(
        test_name,
        "box_below_baseline",
        &[("box", "scalar_box")],
        &[],
    )
    .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command(SCALAR_BOX_ABOVE_EXTERNALS)?;

    let triangle_target = single_slot_integral(&triangle_baseline).result;
    let box_above_target = single_slot_integral(&box_above_baseline).result;
    let box_below_target = single_slot_integral(&box_below_baseline).result;

    cli.run_command(&format!(
        "set process -p triangle -i scalar_tri {shared_integrator_settings}"
    ))?;
    let triangle_box_workspace = get_tests_workspace_path()
        .join(test_name)
        .join("triangle_and_box_correlated");
    scalar_topology_integrate_command(
        test_name,
        "triangle_and_box_correlated",
        &[("triangle", "scalar_tri"), ("box", "scalar_box")],
        &[
            ("triangle@scalar_tri", triangle_target),
            ("box@scalar_box", box_above_target),
        ],
    )
    .run(&mut cli.state, &cli.cli_settings)?;
    let triangle_box_result =
        load_integration_result(&triangle_box_workspace.join("integration_result.json"))?;

    assert_eq!(triangle_box_result.slots.len(), 2);
    let triangle_slot = triangle_box_result
        .slot("triangle@scalar_tri")
        .expect("triangle slot must be present");
    let box_slot = triangle_box_result
        .slot("box@scalar_box")
        .expect("box slot must be present");
    assert_eq!(triangle_slot.target, Some(triangle_target));
    assert_eq!(box_slot.target, Some(box_above_target));
    assert_eq!(triangle_slot.table_results.len(), 2);
    assert_eq!(box_slot.table_results.len(), 2);
    assert!(triangle_slot.integration_statistics.num_evals > 0);
    assert!(box_slot.integration_statistics.num_evals > 0);
    assert!(!triangle_slot.max_weight_info.is_empty());
    assert!(!box_slot.max_weight_info.is_empty());
    assert!(
        triangle_slot
            .integral
            .is_compatible_with_result(single_slot_integral(&triangle_baseline), 4),
        "triangle correlated result deviates from baseline beyond 4 sigma: {} vs {}",
        triangle_slot.integral,
        single_slot_integral(&triangle_baseline),
    );
    assert!(
        box_slot
            .integral
            .is_compatible_with_result(single_slot_integral(&box_above_baseline), 4),
        "box correlated result deviates from baseline beyond 4 sigma: {} vs {}",
        box_slot.integral,
        single_slot_integral(&box_above_baseline),
    );
    assert!(
        triangle_box_workspace
            .join("integrands/triangle@scalar_tri/settings.toml")
            .exists()
    );
    assert!(
        triangle_box_workspace
            .join("integrands/box@scalar_box/settings.toml")
            .exists()
    );

    cli.run_command(
        "duplicate integrand -p box -i scalar_box --output_process_name box_copy --output_integrand_name scalar_box_copy",
    )?;
    cli.run_command(SCALAR_BOX_COPY_BELOW_EXTERNALS)?;
    let output_process = ProcessRef::Unqualified("box_copy".to_string());
    let output_integrand = "scalar_box_copy".to_string();
    cli.state
        .find_integrand_ref(Some(&output_process), Some(&output_integrand))?;

    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    let box_pair_workspace = get_tests_workspace_path()
        .join(test_name)
        .join("box_above_and_below_correlated");
    scalar_topology_integrate_command(
        test_name,
        "box_above_and_below_correlated",
        &[("box", "scalar_box"), ("box_copy", "scalar_box_copy")],
        &[
            ("box@scalar_box", box_above_target),
            ("box_copy@scalar_box_copy", box_below_target),
        ],
    )
    .run(&mut cli.state, &cli.cli_settings)?;
    let box_pair_result =
        load_integration_result(&box_pair_workspace.join("integration_result.json"))?;

    assert_eq!(box_pair_result.slots.len(), 2);
    let box_above_slot = box_pair_result
        .slot("box@scalar_box")
        .expect("box slot must be present");
    let box_below_slot = box_pair_result
        .slot("box_copy@scalar_box_copy")
        .expect("box_copy slot must be present");
    assert_eq!(box_above_slot.target, Some(box_above_target));
    assert_eq!(box_below_slot.target, Some(box_below_target));
    assert_eq!(box_above_slot.table_results.len(), 2);
    assert_eq!(box_below_slot.table_results.len(), 2);
    assert!(box_above_slot.integration_statistics.num_evals > 0);
    assert!(box_below_slot.integration_statistics.num_evals > 0);
    assert!(!box_above_slot.max_weight_info.is_empty());
    assert!(!box_below_slot.max_weight_info.is_empty());
    assert!(
        box_above_slot
            .integral
            .is_compatible_with_result(single_slot_integral(&box_above_baseline), 4),
        "leading box correlated result deviates from baseline beyond 4 sigma: {} vs {}",
        box_above_slot.integral,
        single_slot_integral(&box_above_baseline),
    );
    assert!(
        box_below_slot
            .integral
            .is_compatible_with_result(single_slot_integral(&box_below_baseline), 4),
        "secondary box correlated result deviates from baseline beyond 4 sigma: {} vs {}",
        box_below_slot.integral,
        single_slot_integral(&box_below_baseline),
    );
    assert!(
        box_pair_workspace
            .join("integrands/box@scalar_box/settings.toml")
            .exists()
    );
    assert!(
        box_pair_workspace
            .join("integrands/box_copy@scalar_box_copy/settings.toml")
            .exists()
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn test_multi_integrand_batching_preserves_results() -> Result<()> {
    let test_name = "test_multi_integrand_batching_preserves_results";
    let mut cli = setup_scalar_topologies_cli(test_name)?;
    let shared_integrator_settings = "kv integrator.n_start=5000 integrator.n_max=10000 integrator.n_increase=5000 integrator.seed=1337";

    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    cli.run_command(SCALAR_BOX_BELOW_EXTERNALS)?;

    let mut single_batch = scalar_topology_integrate_command(
        test_name,
        "box_single_batch",
        &[("box", "scalar_box")],
        &[],
    );
    single_batch.batch_size = Some(10_000);
    single_batch.no_stream_iterations = true;
    single_batch.no_stream_updates = true;
    single_batch.run(&mut cli.state, &cli.cli_settings)?;

    let mut many_batches = scalar_topology_integrate_command(
        test_name,
        "box_many_batches",
        &[("box", "scalar_box")],
        &[],
    );
    many_batches.batch_size = Some(1);
    many_batches.no_stream_iterations = true;
    many_batches.no_stream_updates = true;
    many_batches.run(&mut cli.state, &cli.cli_settings)?;

    let single_batch_result = load_integration_result(
        &get_tests_workspace_path()
            .join(test_name)
            .join("box_single_batch")
            .join("integration_result.json"),
    )?;
    let many_batches_result = load_integration_result(
        &get_tests_workspace_path()
            .join(test_name)
            .join("box_many_batches")
            .join("integration_result.json"),
    )?;

    assert_integration_results_match_ignoring_timings(&single_batch_result, &many_batches_result);

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn test_multi_integrand_with_local_model_parameters() -> Result<()> {
    let test_name = "test_multi_integrand_with_local_model_parameters";
    let mut cli = setup_scalar_topologies_cli(test_name)?;
    let shared_integrator_settings = "kv integrator.n_start=5000 integrator.n_max=10000 integrator.n_increase=5000 integrator.seed=1337";

    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    let box_default_baseline = scalar_topology_integrate_command(
        test_name,
        "box_default_model",
        &[("box", "scalar_box")],
        &[],
    )
    .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command(
        "duplicate integrand -p box -i scalar_box --output_process_name box_copy --output_integrand_name scalar_box_copy",
    )?;
    cli.run_command("set model -p box_copy mass_scalar_2=1.0")?;
    cli.run_command(&format!(
        "set process -p box_copy -i scalar_box_copy {shared_integrator_settings}"
    ))?;

    let box_process = ProcessRef::Unqualified("box".to_string());
    let box_integrand = "scalar_box".to_string();
    let (box_process_id, box_integrand_name) = cli
        .state
        .find_integrand_ref(Some(&box_process), Some(&box_integrand))?;
    let box_card = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(
            box_process_id,
            &box_integrand_name,
        )?;
    assert_eq!(box_card.data.get("mass_scalar_2"), Some(&(F(2.0), F(0.0))));

    let box_copy_process = ProcessRef::Unqualified("box_copy".to_string());
    let box_copy_integrand = "scalar_box_copy".to_string();
    let (box_copy_process_id, box_copy_integrand_name) = cli
        .state
        .find_integrand_ref(Some(&box_copy_process), Some(&box_copy_integrand))?;
    let box_copy_card = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(
            box_copy_process_id,
            &box_copy_integrand_name,
        )?;
    assert_eq!(
        box_copy_card.data.get("mass_scalar_2"),
        Some(&(F(1.0), F(0.0)))
    );

    let box_local_baseline = scalar_topology_integrate_command(
        test_name,
        "box_local_model",
        &[("box_copy", "scalar_box_copy")],
        &[],
    )
    .run(&mut cli.state, &cli.cli_settings)?;

    let box_default_target = single_slot_integral(&box_default_baseline).result;
    let box_local_target = single_slot_integral(&box_local_baseline).result;

    let shared_workspace = get_tests_workspace_path()
        .join(test_name)
        .join("box_default_and_local_model");
    scalar_topology_integrate_command(
        test_name,
        "box_default_and_local_model",
        &[("box", "scalar_box"), ("box_copy", "scalar_box_copy")],
        &[
            ("box@scalar_box", box_default_target),
            ("box_copy@scalar_box_copy", box_local_target),
        ],
    )
    .run(&mut cli.state, &cli.cli_settings)?;
    let multi_result = load_integration_result(&shared_workspace.join("integration_result.json"))?;

    let box_default_slot = multi_result
        .slot("box@scalar_box")
        .expect("box slot must be present");
    let box_local_slot = multi_result
        .slot("box_copy@scalar_box_copy")
        .expect("box_copy slot must be present");
    assert!(
        box_default_slot
            .integral
            .is_compatible_with_result(single_slot_integral(&box_default_baseline), 4),
        "default-model box correlated result deviates from baseline beyond 4 sigma: {} vs {}",
        box_default_slot.integral,
        single_slot_integral(&box_default_baseline),
    );
    assert!(
        box_local_slot
            .integral
            .is_compatible_with_result(single_slot_integral(&box_local_baseline), 4),
        "local-model box correlated result deviates from baseline beyond 4 sigma: {} vs {}",
        box_local_slot.integral,
        single_slot_integral(&box_local_baseline),
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn test_inspect_uses_per_integrand_model_parameters() -> Result<()> {
    let test_name = "test_inspect_uses_per_integrand_model_parameters";
    let mut cli = setup_scalar_topologies_cli(test_name)?;

    cli.run_command(
        "duplicate integrand -p box -i scalar_box --output_process_name box_copy --output_integrand_name scalar_box_copy",
    )?;

    let inspect_point = vec![0.31, 0.52, 0.73];
    let (_, default_inspect) = Inspect {
        process: Some(ProcessRef::Unqualified("box".to_string())),
        integrand_name: Some("scalar_box".to_string()),
        point: inspect_point.clone(),
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;
    let (_, copied_inspect_before_override) = Inspect {
        process: Some(ProcessRef::Unqualified("box_copy".to_string())),
        integrand_name: Some("scalar_box_copy".to_string()),
        point: inspect_point.clone(),
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;

    assert!(
        complex_distance(default_inspect, copied_inspect_before_override) < 1e-12,
        "duplicated integrand should initially inspect identically: {} vs {}",
        default_inspect,
        copied_inspect_before_override,
    );

    cli.run_command("set model -p box_copy mass_scalar_2=1.0")?;

    let (_, copied_inspect_after_override) = Inspect {
        process: Some(ProcessRef::Unqualified("box_copy".to_string())),
        integrand_name: Some("scalar_box_copy".to_string()),
        point: inspect_point.clone(),
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;

    assert!(
        complex_distance(default_inspect, copied_inspect_after_override) > 1e-12,
        "per-integrand model override should change the inspect result: {} vs {}",
        default_inspect,
        copied_inspect_after_override,
    );

    let box_copy_process = ProcessRef::Unqualified("box_copy".to_string());
    let box_copy_integrand = "scalar_box_copy".to_string();
    let (box_copy_process_id, box_copy_integrand_name) = cli
        .state
        .find_integrand_ref(Some(&box_copy_process), Some(&box_copy_integrand))?;
    let box_copy_card = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(
            box_copy_process_id,
            &box_copy_integrand_name,
        )?;
    assert_eq!(
        box_copy_card.data.get("mass_scalar_2"),
        Some(&(F(1.0), F(0.0)))
    );

    cli.run_command("set process -p box_copy -i scalar_box_copy defaults")?;

    let box_copy_reset_card = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(
            box_copy_process_id,
            &box_copy_integrand_name,
        )?;
    assert_eq!(
        box_copy_reset_card.data.get("mass_scalar_2"),
        Some(&(F(2.0), F(0.0)))
    );

    cli.run_command(SCALAR_BOX_COPY_ABOVE_EXTERNALS)?;

    let (_, copied_inspect_after_defaults) = Inspect {
        process: Some(ProcessRef::Unqualified("box_copy".to_string())),
        integrand_name: Some("scalar_box_copy".to_string()),
        point: inspect_point,
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;

    assert!(
        complex_distance(default_inspect, copied_inspect_after_defaults) < 1e-12,
        "resetting process defaults should clear the local model override: {} vs {}",
        default_inspect,
        copied_inspect_after_defaults,
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn test_integration_workspace_model_mismatch_requires_restart() -> Result<()> {
    let test_name = "test_integration_workspace_model_mismatch_requires_restart";
    let mut cli = setup_scalar_topologies_cli(test_name)?;
    let shared_integrator_settings = "kv integrator.n_start=5000 integrator.n_max=10000 integrator.n_increase=5000 integrator.seed=1337";

    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    let workspace_name = "box_workspace_model_mismatch";
    scalar_topology_integrate_command(test_name, workspace_name, &[("box", "scalar_box")], &[])
        .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command("set model mass_scalar_2=1.0")?;
    let mut resume_command =
        scalar_topology_integrate_command(test_name, workspace_name, &[("box", "scalar_box")], &[]);
    resume_command.restart = false;
    let err = resume_command
        .run(&mut cli.state, &cli.cli_settings)
        .expect_err("workspace resume should fail on model-parameter mismatch");
    let err_text = format!("{err:?}");
    assert!(err_text.contains("Workspace effective model parameters do not match"));
    assert!(err_text.contains("box@scalar_box"));

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn test_integration_workspace_resume_preserves_current_effective_model_override() -> Result<()> {
    let test_name = "test_integration_workspace_resume_preserves_current_effective_model_override";
    let mut cli = setup_scalar_topologies_cli(test_name)?;
    let shared_integrator_settings = "kv integrator.n_start=5000 integrator.n_max=10000 integrator.n_increase=5000 integrator.seed=1337";

    cli.run_command(&format!(
        "set process -p box -i scalar_box {shared_integrator_settings}"
    ))?;
    let workspace_name = "box_workspace_effective_model_match";
    scalar_topology_integrate_command(test_name, workspace_name, &[("box", "scalar_box")], &[])
        .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command("set model mass_scalar_2=1.0")?;
    cli.run_command("set process -p box -i scalar_box kv model.mass_scalar_2=2.0")?;

    let mut resume_command =
        scalar_topology_integrate_command(test_name, workspace_name, &[("box", "scalar_box")], &[]);
    resume_command.restart = false;
    resume_command.run(&mut cli.state, &cli.cli_settings)?;

    let (process_id, integrand_name) = cli.state.find_integrand_ref(
        Some(&ProcessRef::Unqualified("box".to_string())),
        Some(&"scalar_box".to_string()),
    )?;
    let effective_card = cli
        .state
        .resolve_effective_model_parameter_card_for_integrand(process_id, &integrand_name)?;
    assert_eq!(
        effective_card.data.get("mass_scalar_2"),
        Some(&(F(2.0), F(0.0)))
    );

    let mut second_resume =
        scalar_topology_integrate_command(test_name, workspace_name, &[("box", "scalar_box")], &[]);
    second_resume.restart = false;
    second_resume.run(&mut cli.state, &cli.cli_settings)?;

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn test_simple_workflow_example_card() -> Result<()> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("tests crate should live in the workspace root");
    let mut run_history =
        RunHistory::load(workspace_root.join("examples/command_cards/simple_workflow.toml"))?;
    let state_path = get_tests_workspace_path().join("simple_workflow_example");

    clean_test(&state_path);

    let mut state = new_cli_for_test(&state_path, None);
    run_history.cli_settings.state.folder = state_path.clone();
    let mut cli_settings = run_history.cli_settings.clone();
    cli_settings.sync_settings()?;
    let mut default_runtime_settings = run_history.default_runtime_settings.clone();

    let _ = run_history.run(&mut state, &mut cli_settings, &mut default_runtime_settings)?;

    assert!(
        !state.process_list.processes.is_empty(),
        "No processes were generated"
    );
    assert_eq!(state.model.name, "sm", "Expected SM model to be loaded");
    assert!(
        state_path.join("justfile").exists(),
        "Expected dot export helper files to be created in {}",
        state_path.display()
    );

    clean_test(&state_path);
    Ok(())
}

#[test]
fn test_advanced_integration_example_card() -> Result<()> {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("tests crate should live in the workspace root");
    let run_history =
        RunHistory::load(workspace_root.join("examples/command_cards/advanced_integration.toml"))?;

    assert_eq!(
        run_history.commands[0].raw_string.as_deref(),
        Some("import model sm-default")
    );
    assert_eq!(
        run_history.cli_settings.state.folder,
        PathBuf::from("./advanced_integration_state")
    );
    assert_eq!(run_history.default_runtime_settings.kinematics.e_cm, 500.0);
    Ok(())
}

#[test]
fn test_epem_tth_inspect_nlo_gl18() -> Result<()> {
    let mut cli = get_test_cli(
        Some("test_epem_tth_inspect_nlo_gl18.toml".into()),
        get_tests_workspace_path().join("test_epem_tth_inspect_nlo_gl18"),
        None,
        true,
    )
    .unwrap();

    let (_, inspect) = Inspect {
        process: None,
        integrand_name: None,
        point: vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.823, 0.214],
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)
    .unwrap();

    let target = Complex::new(-9.487984855932107e-6, 3.610476200052732e-5);
    assert_eq!(inspect, target);
    Ok(())
}

#[test]
fn test_qqx_aaa_ir_tree_unprocessed_inspect() -> Result<()> {
    let mut cli = get_test_cli(
        Some("generate_qqx_aaa_tree_unprocessed.toml".into()),
        get_tests_workspace_path().join("qqx_aaa_tree_unprocessed"),
        None,
        true,
    )
    .unwrap();

    let (_, inspect) = Inspect {
        process: None,
        integrand_name: None,
        point: vec![],
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)
    .unwrap();

    let target = Complex::new(0.00014727604164105595, -0.001150313936913021);
    assert_eq!(inspect, target);
    Ok(())
}

#[test]
fn test_qqx_aaa_ir_tree_user_numerator_unprocessed_with_momtrop_table_inspect() -> Result<()> {
    let mut cli = get_test_cli(
        Some("generate_qqx_aaa_tree_user_numerator_unprocessed_with_momtrop_table.toml".into()),
        get_tests_workspace_path()
            .join("qqx_aaa_tree_user_numerator_unprocessed_with_momtrop_table"),
        None,
        true,
    )
    .unwrap();

    let (_, inspect) = Inspect {
        process: None,
        integrand_name: None,
        point: vec![],
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)
    .unwrap();

    let target = Complex::new(1.47276041641056e-4, -1.1503139369130214e-3);
    assert_eq!(inspect, target);
    Ok(())
}

#[test]
fn test_qqx_aaa_ir_tree_user_numerator_inspect() -> Result<()> {
    let mut cli = get_test_cli(
        Some("generate_qqx_aaa_tree_user_numerator.toml".into()),
        get_tests_workspace_path().join("qqx_aaa_tree_user_numerator"),
        None,
        true,
    )
    .unwrap();

    let (_, inspect) = Inspect {
        process: None,
        integrand_name: None,
        point: vec![],
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)
    .unwrap();

    let target = Complex::new(1.47276041641056e-4, -1.1503139369130214e-3);
    assert_eq!(inspect, target);
    Ok(())
}

#[test]
fn epemttb_generate() -> Result<()> {
    let cli = get_test_cli(
        Some("epemttbar.toml".into()),
        get_tests_workspace_path().join("epemttbar"),
        None,
        true,
    )
    .unwrap();

    Ok(())
}

#[test]
fn test_qqx_aaa_ir_subtracted_inspect() -> Result<()> {
    let mut cli = get_test_cli(
        Some("generate_qqx_aaa_ir_subtracted_physical.toml".into()),
        get_tests_workspace_path().join("qqx_aaa_ir_subtracted_physical"),
        None,
        true,
    )
    .unwrap();

    let (_, inspect) = Inspect {
        process: None,
        integrand_name: None,
        point: vec![0.0, 0.10001, 0.10001],
        momentum_space: true,
        ..Default::default()
    }
    .run(&mut cli)
    .unwrap();

    let target = Complex::new(2.3159767780905335e-1, -1.8547720156633686e-4);
    assert_eq!(inspect, target);
    Ok(())
}
#[test]
fn test_integrate_dotted_bubble() -> Result<()> {
    let target = Complex::new(F(0.0), F(0.002871504663657095));
    let mut cli = get_test_cli(
        Some("dotted_bubble_generate.toml".into()),
        get_tests_workspace_path().join("dotted_bubble"),
        None,
        true,
    )?;

    let integrate_command = Integrate {
        workspace_path: Some(
            get_tests_workspace_path().join("dotted_bubble/integration_workspace"),
        ),
        n_cores: Some(1),

        restart: true,
        ..Default::default()
    };

    let integration_result = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integration_result.is_compatible_with_target(target, 2),
        "Integration result {integration_result} not compatible with target {target}"
    );

    Ok(())
}

#[test]
fn test_integrate_bubble_dot_1t1b() -> Result<()> {
    let target = Complex::new(F(-0.00041875868403258484), F(0.0));

    let mut cli = get_test_cli(
        Some("bubble_dot_1t1b_generate.toml".into()),
        get_tests_workspace_path().join("bubble_dot_1t1b"),
        None,
        true,
    )?;

    let integrate_command = Integrate {
        workspace_path: Some(
            get_tests_workspace_path().join("bubble_dot_1t1b/integration_workspace"),
        ),
        n_cores: Some(1),
        ..Default::default()
    };

    let integration_result = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integration_result.is_compatible_with_target(target, 2),
        "Integration result {integration_result} not compatible with target {target}"
    );

    Ok(())
}

#[test]
fn test_integrate_bubble_dot_2t0b() -> Result<()> {
    let target = Complex::new(F(-2.9911664985748502e-05), F(0.0));

    let mut cli = get_test_cli(
        Some("bubble_dot_2t0b_generate.toml".into()),
        get_tests_workspace_path().join("bubble_dot_2t0b"),
        None,
        true,
    )?;

    let integrate_command = Integrate {
        workspace_path: Some(
            get_tests_workspace_path().join("bubble_dot_2t0b/integration_workspace"),
        ),
        n_cores: Some(1),

        restart: true,
        ..Default::default()
    };

    let integration_result = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integration_result.is_compatible_with_target(target, 2),
        "Integration result {integration_result} not compatible with target {target}"
    );
    Ok(())
}

#[test]
fn test_mass_approach_scalar_self_energy() -> Result<()> {
    let mut cli = get_test_cli(
        Some("mass_approach_scalar_self_energy.toml".into()),
        get_tests_workspace_path().join("mass_approach_scalar_self_energy"),
        None,
        true,
    )?;

    let mass_values = [1.1, 1.05, 1.01, 1.001, 1.0001];
    let mut inspect_magnitudes = Vec::with_capacity(mass_values.len());

    for mass in mass_values {
        cli.run_command(&format!("set model mass_scalar_2={mass}"))?;

        let (_, inspect) = Inspect {
            process: None,
            integrand_name: Some("default".to_string()),
            point: vec![0.1, 0.2, 0.3, 0.3, 0.4, 0.5],
            momentum_space: false,
            use_arb_prec: true,
            ..Default::default()
        }
        .run(&mut cli)?;

        let magnitude = (inspect.re * inspect.re + inspect.im * inspect.im).sqrt();
        inspect_magnitudes.push(magnitude);
    }

    assert!(
        inspect_magnitudes.windows(2).all(|pair| pair[1] <= pair[0]),
        "Inspect magnitude is not monotonically decreasing as mass_scalar_2 approaches 1: {inspect_magnitudes:?}"
    );
    assert!(
        inspect_magnitudes.last().copied().unwrap_or(f64::INFINITY) < 1.0e-11,
        "Inspect did not approach zero closely enough near mass_scalar_2=1: {inspect_magnitudes:?}"
    );

    Ok(())
}

fn configure_differential_leading_jet_observable(
    cli: &mut gammaloop_integration_tests::CLIState,
) -> Result<()> {
    run_commands(
        cli,
        &[
            r#"set process string '
[quantities.leading_jet_pt]
type = "jet"
quantity = "PT"
dR = 0.4
'"#,
            r#"set process string '
[observables.leading_jet_pt_hist]
quantity = "leading_jet_pt"
entry_selection = "leading_only"
x_min = 0.0
x_max = 1000.0
n_bins = 8
'"#,
        ],
    )
}

fn configure_differential_leading_jet_selector(
    cli: &mut gammaloop_integration_tests::CLIState,
) -> Result<()> {
    cli.run_command(
        r#"set process string '
[selectors.leading_jet_pt_cut]
quantity = "leading_jet_pt"
selector = "value_range"
entry_selection = "leading_only"
min = 0.0
'"#,
    )
}

#[test]
#[serial]
fn lu_differential_integration_writes_json_observables() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_differential_integration_json")?;
    configure_differential_leading_jet_observable(&mut cli)?;
    configure_differential_leading_jet_selector(&mut cli)?;
    cli.run_command(
        "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=json",
    )?;

    let workspace = get_tests_workspace_path().join("lu_differential_integration_json/workspace");
    let integration_result = Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    assert!(single_slot_integral(&integration_result).neval > 0);
    assert!(workspace.join("integration_result.json").exists());
    let slot_workspace = selected_slot_workspace(&cli, &workspace, None, Some("default"))?;
    let final_file = slot_workspace.join("observables_final.json");
    assert!(final_file.exists());
    assert!(
        !slot_workspace
            .join("observables_final_iter_0001.json")
            .exists()
    );

    let final_bundle =
        gammalooprs::observables::ObservableSnapshotBundle::from_json_file(&final_file)?;
    assert!(final_bundle.histograms.contains_key("leading_jet_pt_hist"));

    Ok(())
}

#[test]
#[serial]
fn lu_differential_integration_cli_flag_writes_iteration_observables() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_differential_integration_json")?;
    configure_differential_leading_jet_observable(&mut cli)?;
    configure_differential_leading_jet_selector(&mut cli)?;
    cli.run_command(
        "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=json",
    )?;

    let workspace =
        get_tests_workspace_path().join("lu_differential_integration_json/cli_flag_workspace");
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        write_results_for_each_iteration: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let slot_workspace = selected_slot_workspace(&cli, &workspace, None, Some("default"))?;
    assert!(workspace.join("integration_result.json").exists());
    assert!(
        workspace
            .join("results/integration_result_iter_0001.json")
            .exists()
    );
    assert!(
        slot_workspace
            .join("observables_final_iter_0001.json")
            .exists()
    );
    assert!(slot_workspace.join("observables_final.json").exists());

    Ok(())
}

#[test]
#[serial]
fn lu_differential_integration_hwu_output_is_optional_and_single_file() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_differential_integration_hwu")?;
    cli.run_command(
        "set process kv integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=hwu",
    )?;

    let workspace_without_observables =
        get_tests_workspace_path().join("lu_differential_integration_hwu/without_observables");
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace_without_observables.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        workspace_without_observables
            .join("integration_result.json")
            .exists()
    );
    assert!(
        !selected_slot_workspace(&cli, &workspace_without_observables, None, Some("default"))?
            .join("observables_final.hwu")
            .exists()
    );

    configure_differential_leading_jet_observable(&mut cli)?;
    configure_differential_leading_jet_selector(&mut cli)?;
    let workspace_with_observables =
        get_tests_workspace_path().join("lu_differential_integration_hwu/with_observables");
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace_with_observables.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        workspace_with_observables
            .join("integration_result.json")
            .exists()
    );

    let slot_workspace =
        selected_slot_workspace(&cli, &workspace_with_observables, None, Some("default"))?;
    let hwu_file = slot_workspace.join("observables_final.hwu");
    assert!(hwu_file.exists());
    let hwu_contents = std::fs::read_to_string(hwu_file)?;
    assert!(hwu_contents.contains("<histogram>"));
    assert!(hwu_contents.contains("leading_jet_pt_hist"));
    assert!(
        !slot_workspace
            .join("observables_final_iter_0001.hwu")
            .exists()
    );

    Ok(())
}

#[test]
#[serial]
fn lu_differential_json_observables_resume_from_workspace() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_differential_integration_json")?;
    configure_differential_leading_jet_observable(&mut cli)?;
    configure_differential_leading_jet_selector(&mut cli)?;
    cli.run_command(
        "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=json",
    )?;

    let workspace =
        get_tests_workspace_path().join("lu_differential_integration_json/resume_workspace");
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let slot_workspace = selected_slot_workspace(&cli, &workspace, None, Some("default"))?;
    let (process_id, resolved_integrand_name) = cli
        .state
        .find_integrand_ref(None, Some(&"default".to_string()))?;
    let slot_meta = gammalooprs::integrate::SlotMeta {
        process_name: cli.state.process_list.processes[process_id]
            .definition
            .folder_name
            .clone(),
        integrand_name: resolved_integrand_name,
    };
    let final_file = slot_workspace.join("observables_final.json");
    let checkpoint_file =
        gammalooprs::integrate::latest_observable_resume_state_path(&workspace, &slot_meta);
    assert!(checkpoint_file.exists());
    let result_snapshot_path = workspace.join("integration_result.json");
    let result_before_resume = std::fs::read_to_string(&result_snapshot_path)?;
    let state_before_resume =
        std::fs::read(gammalooprs::integrate::workspace_state_path(&workspace))?;
    let before_resume =
        gammalooprs::observables::ObservableSnapshotBundle::from_json_file(&final_file)?;

    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: false,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let after_resume =
        gammalooprs::observables::ObservableSnapshotBundle::from_json_file(&final_file)?;
    assert_eq!(before_resume, after_resume);
    assert_eq!(
        result_before_resume,
        std::fs::read_to_string(result_snapshot_path)?
    );
    assert_eq!(
        state_before_resume,
        std::fs::read(gammalooprs::integrate::workspace_state_path(&workspace))?
    );

    Ok(())
}

#[test]
#[serial]
fn lu_differential_hwu_observables_resume_from_workspace() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_differential_integration_hwu")?;
    configure_differential_leading_jet_observable(&mut cli)?;
    configure_differential_leading_jet_selector(&mut cli)?;
    cli.run_command(
        "set process kv general.generate_events=false integrator.n_start=12 integrator.min_samples_for_update=12 integrator.n_max=12 integrator.n_increase=0 integrator.observables_output.format=hwu",
    )?;

    let workspace =
        get_tests_workspace_path().join("lu_differential_integration_hwu/resume_workspace");
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let slot_workspace = selected_slot_workspace(&cli, &workspace, None, Some("default"))?;
    let final_file = slot_workspace.join("observables_final.hwu");
    let before_resume = std::fs::read_to_string(&final_file)?;
    let (process_id, resolved_integrand_name) = cli
        .state
        .find_integrand_ref(None, Some(&"default".to_string()))?;
    let slot_meta = gammalooprs::integrate::SlotMeta {
        process_name: cli.state.process_list.processes[process_id]
            .definition
            .folder_name
            .clone(),
        integrand_name: resolved_integrand_name,
    };
    let checkpoint_file =
        gammalooprs::integrate::latest_observable_resume_state_path(&workspace, &slot_meta);
    assert!(checkpoint_file.exists());
    let result_snapshot_path = workspace.join("integration_result.json");
    let result_before_resume = std::fs::read_to_string(&result_snapshot_path)?;
    let state_before_resume =
        std::fs::read(gammalooprs::integrate::workspace_state_path(&workspace))?;

    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        workspace_path: Some(workspace.clone()),
        target: vec![],
        n_cores: Some(1),
        restart: false,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let after_resume = std::fs::read_to_string(final_file)?;
    assert_eq!(before_resume, after_resume);
    assert_eq!(
        result_before_resume,
        std::fs::read_to_string(result_snapshot_path)?
    );
    assert_eq!(
        state_before_resume,
        std::fs::read(gammalooprs::integrate::workspace_state_path(&workspace))?
    );

    Ok(())
}

#[test]
#[serial]
fn lu_save_dot_silently_overwrites_existing_files() -> Result<()> {
    let mut cli = setup_sm_differential_lu_cli("lu_save_dot_overwrite")?;
    cli.run_command("save dot")?;
    cli.run_command("save dot")?;
    Ok(())
}
