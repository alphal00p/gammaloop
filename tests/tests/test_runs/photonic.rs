use super::utils::*;
use super::*;

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

        let real_part_string = gloop_utils::format_uncertainty(
            integrated_result.result.re,
            integrated_result.error.re,
        );
        let imag_part_string = gloop_utils::format_uncertainty(
            integrated_result.result.im,
            integrated_result.error.im,
        );

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
