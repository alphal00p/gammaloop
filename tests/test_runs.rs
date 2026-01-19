#![allow(dead_code)]
#![allow(unused_variables)]

use color_eyre::Result;

use colored::{ColoredString, Colorize};
use eyre::Ok;
use gammaloop_api::commands::{inspect::Inspect, integrate::Integrate};

use gammalooprs::{
    GammaLoopContextContainer,
    model::UFOSymbol,
    numerator::aind::Aind,
    utils::{self, ApproxEq, F, FUN_LIB, GS, TENSORLIB, W_},
};
use insta::assert_snapshot;
use itertools::Itertools;
use libc::NF_IP_PRI_RAW_BEFORE_DEFRAG;
use momtrop::assert_approx_eq;
use serde_json::to_string;
use spenso::{
    algebra::complex::Complex,
    network::{Network, Sequential, SmallestDegree, store::NetworkStore},
    structure::{IndexlessNamedStructure, NamedStructure, representation::LibraryRep},
    tensors::{complex::RealOrComplexTensor, parametric::ParamOrConcrete},
};
use std::{env, fmt::Write, ops::Deref, time::Duration};
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    function,
    id::Pattern,
    symbol,
};
use tabled::{Table, Tabled, settings::Style};
use tracing::info;
use tracing_subscriber::fmt::format;

mod test_utils;
use test_utils::{clean_test, get_test_cli, get_tests_workspace_path};

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
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(
            get_tests_workspace_path()
                .join("z_decay_test/integration_workspace/integration_results.yaml"),
        ),
        workspace_path: Some(
            get_tests_workspace_path().join("z_decay_test/integration_workspace/"),
        ),
        target: None,
        n_cores: Some(1),
        restart: true,
    }
    .run(&mut state.state, &state.cli_settings)?;

    assert_approx_eq(&result.result.re, &F(-0.372), &F(1e-2));
    Ok(())
}

#[test]
fn test_epem_dd_dt() -> Result<()> {
    let state = get_test_cli(
        Some("test_epem_dd_dt.toml".into()),
        get_tests_workspace_path().join("test_epem_dd_dt"),
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
        process_id: Some(0),
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
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        ..Default::default()
    }
    .run(&mut cli)?;

    assert_snapshot!(format!("{a:.8e}"),@"(1.4727604164105595e-4+-1.1503139369130214e-3i)");

    clean_test(&cli.cli_settings.state_folder);
    Ok(())
}

#[test]
fn photons_1l_integrate() -> Result<()> {
    let mut cli = get_test_cli(
        Some("photons_eu.toml".into()),
        "./tests/photons_eu_integrate",
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
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(
            "./tests/photons_eu_integrate/integration_workspace/integration_results.yaml".into(),
        ),
        n_cores: Some(1),
        workspace_path: Some("./tests/photons_eu_integrate/integration_workspace/".into()),
        target: Some(vec![target.re.0, target.im.0]),
        restart: true,
    };

    let integration_result = integrate.run(&mut cli.state, &cli.cli_settings)?;

    assert_approx_eq(&integration_result.result.re, &target.re, &F(10e-2));
    assert_approx_eq(&integration_result.result.im, &target.im, &F(10e-2));

    clean_test(&cli.cli_settings.state_folder);
    Ok(())
}

#[test]
fn photons_1l_inspect() -> Result<()> {
    let mut cli = get_test_cli(
        Some("photons_eu.toml".into()),
        "./tests/photons_eu_inspect",
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
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        point: vec![0.123, 0.3242, 0.4233],
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;

    println!("Inspect result: {inspect:.16e}");

    // The old test at a very bad value of e_cm, so I created a new value using the example card in the old main
    assert_snapshot!(format!("{inspect:.8e}"),@"(-4.2365441831364175e-12+-3.71072895861423e-12i)");

    clean_test(&cli.cli_settings.state_folder);

    Ok(())
}

#[test]
fn photons_phys_1l_inspect() -> Result<()> {
    let mut cli = get_test_cli(
        Some("generate_threshold_1L_6photons.toml".into()),
        "./tests/photons_phys_1l_inspect",
        None,
        true,
    )?;
    let (_, inspect) = Inspect {
        process_id: Some(0),
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
        "./tests/photons_eu_2l_inspect",
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
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        point: vec![0.123, 0.3242, 0.4233, 0.523, 0.314, 0.125],
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)?;

    println!("Inspect result: {inspect:.16e}");

    // wrong result
    assert_snapshot!(format!("{inspect:.8e}"));

    clean_test(&cli.cli_settings.state_folder);

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
            table.with(
                Style::modern()
                    .remove_horizontals()
                    .remove_verticals()
                    .remove_frame(),
            );
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
        rsd: ComparisonReport,
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
                rsd: ComparisonReport::empty(),
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
        rsd_bench: Option<Complex<F<f64>>>,
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
            process_id: Some(0),
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
                if inspect_target_f.approx_eq(&inspect_result_f, &F(1.0e-14)) {
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
            process_id: Some(0),
            integrand_name: Some("default".to_string()),
            result_path: None,
            workspace_path: None,
            n_cores: Some(1),
            target: bench_data.integrated_target.map(|t| vec![t.re.0, t.im.0]),
            restart: true,
        }
        .run(&mut cli.state, &cli.cli_settings)?;
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
                        "{:.16e} + {:.16e}i",
                        integrated_target.re.0, integrated_target.im.0
                    )
                    .blue(),
                }
            } else {
                ComparisonReport::only_gloop(
                    format!("{} +  {}i", real_part_string, imag_part_string).blue(),
                )
            };

        let rsd_re = integrated_result.error.re.0 / integrated_result.result.re.0.abs();
        let rsd_im = integrated_result.error.im.0 / integrated_result.result.im.0.abs();

        let rsd_comparison = if let Some(rsd_target) = bench_data.rsd_bench {
            let rsd_string = if rsd_re < rsd_target.re.0 * 2.0 && rsd_im < rsd_target.im.0 * 2.0 {
                format!("{:.2e} + {:.2e}", rsd_re, rsd_im).green()
            } else if rsd_re < rsd_target.re.0 * 4.0 && rsd_im < rsd_target.im.0 * 4.0 {
                format!("{:.2e} + {:.2e}", rsd_re, rsd_im).yellow()
            } else {
                *test_failed = true;
                format!("{:.2e} + {:.2e}", rsd_re, rsd_im).red()
            };
            ComparisonReport {
                gammaloop: rsd_string,
                reference: format!("{:.2e} + {:.2e}", rsd_target.re.0, rsd_target.im.0).blue(),
            }
        } else {
            ComparisonReport::only_gloop(format!("{:.2e} + {:.2e}", rsd_re, rsd_im).blue())
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
            rsd: rsd_comparison,
            sample_time: sample_time_comparison,
        })
    }

    let mut test_reports: Vec<TestReportEntry> = Vec::new();

    let one_loop_eu = BenchMarkData {
        run_card: "photonic_amplitudes/1l_eu.toml".into(),
        state_path: "./tests/photonic_amplitudes/1l_eu".into(),
        amplitude: "1l_eu".into(),
        generation_time: Some(Duration::from_secs(7)),
        inspect_point: vec![0.123, 0.3242, 0.4233],
        inspect_target: Some(Complex::new(-4.236544183136417e-12, -3.710728958614226e-12)),
        integrated_target: Some(Complex::new(
            F(-1.22898408452706e-13),
            F(-3.94362534040412e-13),
        )),
        rsd_bench: None,
        sample_time: Some(Duration::from_micros(61)),
    };

    let one_loop_phys = BenchMarkData {
        run_card: "photonic_amplitudes/1l_phys.toml".into(),
        state_path: "./tests/photonic_amplitudes/1l_phys".into(),
        amplitude: "1l_phys".into(),
        inspect_point: vec![0.1, 0.2, 0.3],
        inspect_target: Some(Complex::new(4.660217572648287e-10, -6.496141401696065e-10)),
        integrated_target: Some(Complex::new(
            F(-9.27759500687454717e-11),
            F(-3.68394576249870544e-11),
        )),
        generation_time: Some(Duration::from_secs(7)),
        rsd_bench: None,
        sample_time: Some(Duration::from_micros(590)),
    };

    test_reports.push(run_photonic_test(one_loop_eu, &mut test_failed)?);
    test_reports.push(run_photonic_test(one_loop_phys, &mut test_failed)?);
    test_reports.push(TestReportEntry::default_with_name("2l_eu".to_string()));
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
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(get_tests_workspace_path().join(
            "test_grouped_subtraction/integration_workspace_no_group/integration_results.yaml",
        )),
        workspace_path: Some(
            get_tests_workspace_path()
                .join("test_grouped_subtraction/integration_workspace_no_group/"),
        ),
        target: None,
        n_cores: Some(1),
        restart: true,
    };
    let int2 =
        Integrate {
            process_id: Some(1),
            integrand_name: Some("default".to_string()),
            result_path: Some(get_tests_workspace_path().join(
                "test_grouped_subtraction/integration_workspace_group/integration_results.yaml",
            )),
            workspace_path: Some(
                get_tests_workspace_path()
                    .join("test_grouped_subtraction/integration_workspace_group/"),
            ),
            target: None,
            n_cores: Some(1),
            restart: true,
        };

    let integration_results_no_group = int1.run(&mut cli.state, &cli.cli_settings)?;
    let integration_results_group = int2.run(&mut cli.state, &cli.cli_settings)?;

    info!("No group result: {:#?}", integration_results_no_group);
    info!("Group result: {:#?}", integration_results_group);

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

    clean_test(&cli.cli_settings.state_folder);

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
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(
            get_tests_workspace_path()
                .join("scalar_bubble/integration_workspace/integration_results.toml"),
        ),
        workspace_path: Some(
            get_tests_workspace_path().join("scalar_bubble/integration_workspace"),
        ),
        n_cores: Some(1),
        target: None,
        restart: true,
    };

    // from Kaapo: m=1 muv=5 2.03838e-02 m=2 muv=5 	1.16050e-02	 m=3 muv=5 6.46968e-03

    cli.run_command("set model mass_scalar_1={re:1.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(2.03838e-02)), 1));

    cli.run_command("set model mass_scalar_1={re:2.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(1.16050e-02)), 1));

    cli.run_command("set model mass_scalar_1={re:3.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(6.46968e-03)), 1));

    clean_test(&cli.cli_settings.state_folder);

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
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(
            get_tests_workspace_path()
                .join("scalar_sunrise/integration_workspace/integration_results.toml"),
        ),
        workspace_path: Some(
            get_tests_workspace_path().join("scalar_sunrise/integration_workspace"),
        ),
        n_cores: Some(1),
        target: None,
        restart: true,
    };

    // from Kaapo: m=1 muv=5 4.37688e-03 m=2 muv=5 	2.48100e-03	 m=3 muv=5 1.07231e-03
    cli.run_command("set model mass_scalar_1={re:1.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(-4.37688e-03)), 1) // .result
                                                                                         // .approx_eq(&Complex::new_re(F(-4.37688e-03)), &F(0.01))
    );

    assert_snapshot!(format!("{integral_no_cache:.3}"),@"-4.359e-3");

    cli.run_command("set model mass_scalar_1={re:2.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert_snapshot!(format!("{integral_no_cache:.3}"),@"-2.474e-3");
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(-2.48100e-03)), 1));

    // cli.run_command("set model mass_scalar_1={re:3.0,im:0.0}")?;
    // let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    // assert_snapshot!(format!("{:.8e}",integral_no_cache.result),@"(-4.5184321377520566e-4+0e0i)");

    clean_test(&cli.cli_settings.state_folder);

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
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(
            get_tests_workspace_path()
                .join("scalar_mercedes/integration_workspace/integration_results.toml"),
        ),
        workspace_path: Some(
            get_tests_workspace_path().join("scalar_mercedes/integration_workspace"),
        ),
        n_cores: Some(1),
        target: None,
        restart: true,
    };

    //5.89551e-06	3.35645e-06	1.87120e-06

    // from Kaapo: m=1 muv=5 5.89551e-06 m=2 muv=5 	3.35645e-06	 m=3 muv=5 1.87120e-06
    cli.run_command("set model mass_scalar_1={re:1.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(5.89551e-06)), 1),
        "Not compatible: {integral_no_cache}",
    );

    cli.run_command("set model mass_scalar_1={re:2.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(3.35645e-06)), 3),
        "Not compatible: {integral_no_cache}",
    );

    cli.run_command("set model mass_scalar_1={re:3.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(1.87120e-06)), 3),
        "Not compatible: {integral_no_cache}",
    );

    clean_test(&cli.cli_settings.state_folder);

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
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(
            get_tests_workspace_path()
                .join("scalar_basketball/integration_workspace/integration_results.toml"),
        ),
        workspace_path: Some(
            get_tests_workspace_path().join("scalar_basketball/integration_workspace"),
        ),
        n_cores: Some(1),
        target: None,
        restart: true,
    };

    //1.47240e-03	7.15184e-04	2.27485e-04

    // from Kaapo: m=1 muv=5 1.47240e-03 m=2 muv=5 	7.15184e-04	 m=3 muv=5 2.27485e-04
    cli.run_command("set model mass_scalar_1={re:1.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(1.47240e-03)), 1),
        "Not compatible: {integral_no_cache}",
    );

    cli.run_command("set model mass_scalar_1={re:2.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(7.15184e-04)), 3),
        "Not compatible: {integral_no_cache}",
    );

    cli.run_command("set model mass_scalar_1={re:3.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(2.27485e-04)), 3),
        "Not compatible: {integral_no_cache}",
    );

    clean_test(&cli.cli_settings.state_folder);

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
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(get_tests_workspace_path().join(
            "scalar_mercedes_with_extra_loop/integration_workspace/integration_results.toml",
        )),
        workspace_path: Some(
            get_tests_workspace_path()
                .join("scalar_mercedes_with_extra_loop/integration_workspace"),
        ),
        n_cores: Some(1),
        target: None,
        restart: true,
    };

    //2.90078e-06	1.59168e-06	6.86001e-07

    // from Kaapo: m=1 muv=5 2.90078e-06 m=2 muv=5 	1.59168e-06	 m=3 muv=5 6.86001e-07
    cli.run_command("set model mass_scalar_1={re:1.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(2.90078e-06)), 1),
        "Not compatible: {integral_no_cache}",
    );

    cli.run_command("set model mass_scalar_1={re:2.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(1.59168e-06)), 3),
        "Not compatible: {integral_no_cache}",
    );

    cli.run_command("set model mass_scalar_1={re:3.0,im:0.0}")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(
        integral_no_cache.is_compatible_with_target(Complex::new_re(F(6.86001e-07)), 3),
        "Not compatible: {integral_no_cache}",
    );

    clean_test(&cli.cli_settings.state_folder);

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

    let point = vec![1., 1., 1., 2., 3., 4.];

    let point = vec![1., 1., 1., -3., -4., -5.];

    let point = vec![2., 3., 4., 1., 1., 1.];
    let mut ins = Inspect {
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        point: point.clone(),
        momentum_space: true,
        discrete_dim: vec![0],
        ..Default::default()
    };

    // from Kaapo: m=1 muv=5 4.37688e-03 m=2 muv=5 	2.48100e-03	 m=3 muv=5 1.07231e-03
    cli.run_command("set model mass_scalar_1={re:1.0,im:0.0}")?;

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
    // clean_test(&cli.cli_settings.state_folder);

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
        process_id: Some(0),
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
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(
            get_tests_workspace_path()
                .join("scalar_box/integration_workspace/integration_results.toml"),
        ),
        workspace_path: Some(get_tests_workspace_path().join("scalar_box/integration_workspace")),
        n_cores: Some(1),
        target: None,
        restart: true,
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    cli.run_command("set process -p 0 -i default kv general.enable_cache=true")?;
    let integral_with_cache = Integrate {
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(
            get_tests_workspace_path()
                .join("scalar_box/integration_workspace/integration_results.toml"),
        ),
        workspace_path: Some(get_tests_workspace_path().join("scalar_box/integration_workspace")),
        target: None,
        n_cores: Some(1),
        restart: true,
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    info!("Integral result without caching: {:#?}", integral_no_cache);
    info!(
        "Integral result with caching.  : {:#?}",
        integral_with_cache
    );

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

    clean_test(&cli.cli_settings.state_folder);

    Ok(())
}

#[test]
fn test_broken_network() -> Result<()> {
    let cli = get_test_cli(
        Some("photon_box.toml".into()),
        get_tests_workspace_path().join("test_broken_network"),
        Some("test_broken_network".to_string()),
        true,
    )?;

    let mut state_bin =
        std::fs::File::open("tests/resources/misc/product_esurface_10_symbolica_state.bin")?;

    let state = symbolica::state::State::import(&mut state_bin, None).unwrap();
    let model = cli.state.model.clone();

    let gammalooopcontainer = GammaLoopContextContainer {
        state_map: &state,
        model: &model,
    };

    let mut file = std::fs::File::open("tests/resources/misc/product_esurface_10.bin")?;

    let mut network: Network<
        NetworkStore<
            ParamOrConcrete<
                RealOrComplexTensor<F<f64>, NamedStructure<Symbol, Vec<Atom>, LibraryRep, Aind>>,
                NamedStructure<Symbol, Vec<Atom>, LibraryRep, Aind>,
            >,
            Atom,
        >,
        IndexlessNamedStructure<Symbol, Vec<Atom>, LibraryRep, Aind>,
        Symbol,
        Aind,
    > = bincode::decode_from_std_read_with_context(
        &mut file,
        bincode::config::standard(),
        gammalooopcontainer,
    )?;

    println!("{}", network.dot_pretty());

    network
        .execute::<Sequential, SmallestDegree, _, _, _>(
            TENSORLIB.read().unwrap().deref(),
            FUN_LIB.deref(),
        )
        .unwrap();

    let scalar: Atom = network.result_scalar().unwrap().into();
    let pattern = Pattern::from(function!(GS.ose, W_.x_) * function!(GS.sign, W_.x_));

    println!("res{scalar}");
    let num_ose_match = scalar.pattern_match(&pattern, None, None);
    //println!("num matches: {}", num_ose_match.count());

    let mut edge_5_in_expression = false;
    let mut edge_7_in_expression = false;

    for matched_exp in num_ose_match {
        println!("found:{}", pattern.replace_wildcards(&matched_exp));
        for (symbol, atom) in matched_exp.iter() {
            // println!("Matched {} to {atom}", symbol);
            if *atom == Atom::num(5) {
                edge_5_in_expression = true;
            }
            if *atom == Atom::num(7) {
                edge_7_in_expression = true;
            }
        }
    }

    assert!(edge_5_in_expression);
    assert!(edge_7_in_expression);

    //println!("Result: {scalar}");

    Ok(())
}

#[test]
fn test_simple_workflow_example_card() -> Result<()> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join("simple_workflow_example"),
        Some("simple_workflow_example".to_string()),
        true,
    )?;

    // Load and execute the simple workflow example card
    cli.run_command("run ../examples/command_cards/simple_workflow.toml")?;

    // Verify that processes were generated
    assert!(
        !cli.state.process_list.processes.is_empty(),
        "No processes were generated"
    );

    // Verify that the Standard Model was loaded
    assert_eq!(cli.state.model.name, "sm", "Expected SM model to be loaded");

    clean_test(&cli.cli_settings.state_folder);
    Ok(())
}

#[test]
fn test_advanced_integration_example_card() -> Result<()> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join("advanced_integration_example"),
        Some("advanced_integration_example".to_string()),
        true,
    )?;

    // Load and execute the advanced integration example card
    cli.run_command("run ../examples/command_cards/advanced_integration.toml")?;

    // Verify that processes were generated
    assert!(
        !cli.state.process_list.processes.is_empty(),
        "No processes were generated"
    );

    // Verify that the Standard Model was loaded
    assert_eq!(cli.state.model.name, "sm", "Expected SM model to be loaded");

    // Verify that settings were applied (check if we can access some expected configuration)
    // The example sets specific integrator settings
    let first_process = &cli.state.process_list.processes[0];
    match &first_process.collection {
        gammalooprs::processes::ProcessCollection::Amplitudes(amplitudes) => {
            // Check that amplitudes were generated
            assert!(!amplitudes.is_empty(), "No amplitudes were generated");
        }
        gammalooprs::processes::ProcessCollection::CrossSections(_) => {
            // Cross sections are also valid
        }
    }

    clean_test(&cli.cli_settings.state_folder);
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
        process_id: None,
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
        process_id: None,
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
        process_id: None,
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
        process_id: None,
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
        process_id: None,
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
fn quick_table_tryout() {
    #[derive(Tabled)]
    struct TopLevelThing {
        name: String,
        foo: String,
        #[tabled(display = "Nested::display_nested")]
        nested: Nested,
    }

    #[derive(Tabled)]
    struct Nested {
        nested_foo: String,
        nested_bar: String,
    }

    impl Nested {
        fn display_nested(&self) -> String {
            let mut table = Table::new([self]);
            table.with(
                Style::modern()
                    .remove_horizontals()
                    .remove_verticals()
                    .remove_frame(),
            );
            table.to_string()
        }
    }

    let things = vec![
        TopLevelThing {
            name: "Thing 1".to_string(),
            foo: "Foo 1".to_string(),
            nested: Nested {
                nested_foo: "Nested Foo 1".to_string(),
                nested_bar: "Nested Bar 1".to_string(),
            },
        },
        TopLevelThing {
            name: "Thing 2".to_string(),
            foo: "Foo 2".to_string(),
            nested: Nested {
                nested_foo: "Nested Foo 2".to_string(),
                nested_bar: "Nested Bar 2".to_string(),
            },
        },
    ];

    let table = Table::new(things).to_string();
    println!("{table}");
}
