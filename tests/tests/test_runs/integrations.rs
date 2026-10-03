use super::utils::*;
use super::*;
use gammaloop_api::commands::integrate::RendererOption;

#[test]
fn fermi_surface_one_loop_integration_converges_to_analytic_target() -> Result<()> {
    let name = "fermi_surface_one_loop_integration";
    let mut cli = get_test_cli(
        Some("fermi_surface_1l_integration.toml".into()),
        get_tests_workspace_path().join(name),
        Some(name.to_string()),
        true,
    )?;
    // For D=q0²-k²-m²+i0, the vacuum-subtracted double pole integrates to
    // I2=-i acosh(mu/m)/(8pi²). Since D^-3 = (1/2) d(D^-2)/d(m²),
    // I3=(1/2) dI2/d(m²)=i mu/(32pi² m² sqrt(mu²-m²)).
    // At m=3, mu=muB/3=5, pF=4 this fixes the absolute normalization and
    // phase independently of generated coefficients and h(t). Both delta
    // and delta' contribute: dropping delta' alone shifts the result by 10%.
    let expected = 5.0 / (32.0 * std::f64::consts::PI.powi(2) * 3.0_f64.powi(2) * 4.0);
    for (profile, sigma) in [("exponential", 1.0), ("poly_exponential", 1.3)] {
        cli.run_command(&format!(
            "set process kv h_function.function=\"{profile}\" h_function.sigma={sigma}"
        ))?;
        let mut errors = Vec::new();
        // Fixed seed and iteration size give identical initial samples in the
        // two fresh runs. Check decreasing uncertainty, not monotonic central
        // values: statistical fluctuations need not approach the target at
        // every iteration. Early stopping is disabled in the run card.
        for samples in [20_000, 100_000] {
            cli.run_command(&format!("set process kv integrator.n_max={samples}"))?;
            let started = std::time::Instant::now();
            let output = Integrate {
                renderer: RendererOption::Tabled,
                show_max_weight_info: false,
                no_stream_iterations: true,
                no_stream_updates: true,
                ..default_integrate_for(name)
            }
            .run(&mut cli.state, &cli.cli_settings)?;
            let estimate = single_slot_integral(&output);
            let (value, error) = (estimate.result.im.0, estimate.error.im.0);
            info!(
                profile,
                sigma,
                samples = estimate.neval,
                value,
                error,
                expected,
                elapsed = ?started.elapsed(),
                "One-loop Fermi-surface integration convergence"
            );
            assert_eq!(estimate.neval, samples, "the sample budget must be fixed");
            assert!(
                [value, error, estimate.result.re.0, estimate.error.re.0]
                    .into_iter()
                    .all(f64::is_finite)
                    && value > 0.0
                    && error > 0.0
                    && estimate.error.re.0 >= 0.0,
                "{profile}, samples={samples}: {estimate}"
            );
            assert!(
                estimate.result.re.0.abs() < 1.0e-12 * expected.abs()
                    && estimate.error.re.0 < 1.0e-12 * expected.abs(),
                "the real part must vanish: {estimate}"
            );
            assert!(
                (value - expected).abs() <= 5.0 * error,
                "{profile}, samples={samples}: {value:e} +/- {error:e}, target={expected:e}"
            );
            errors.push(error);
        }
        assert!(
            errors[1] < 0.005 * expected.abs() && errors[1] < 0.7 * errors[0],
            "{profile}: increased statistics must reduce uncertainty below 0.5%: {errors:?}"
        );
    }
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn fermi_surface_two_loop_integration_converges_in_nonaligned_basis() -> Result<()> {
    use gammalooprs::{momentum::SignOrZero, processes::ProcessCollection};

    let name = "fermi_surface_two_loop_integration";
    let mut cli = get_test_cli(
        Some("fermi_surface_2l_integration.toml".into()),
        get_tests_workspace_path().join(name),
        Some(name.to_string()),
        true,
    )?;
    let ProcessCollection::Amplitudes(amplitudes) = &cli.state.process_list.processes[0].collection
    else {
        panic!("expected the dotted-sunset amplitude")
    };
    let [sunset] = amplitudes["default"].graphs.as_slice() else {
        panic!("expected a single genuine two-loop graph")
    };
    let basis = &sunset.graph.loop_momentum_basis;
    assert_eq!(
        basis
            .loop_edges
            .iter()
            .map(|edge| edge.0)
            .collect::<Vec<_>>(),
        [2, 3],
        "the input basis must use the undotted fermion and the scalar chain"
    );
    for (edge, signature) in &basis.edge_signatures {
        if edge.0 < 2 {
            assert!(
                signature
                    .internal
                    .iter()
                    .all(|sign| *sign != SignOrZero::Zero),
                "the raised fermion route must involve both input loop momenta"
            );
        }
    }
    // Independent finite-density cut integrals, differentiated with respect
    // to the common fermion mass squared. The derivation and quadrature
    // convergence are recorded in fermi_surface_2l_integration.toml. Parameters are
    // m_b=3, mu_b=5, M_H=2; all UV subgraphs converge. Two loops and five
    // propagator factors give Wick phase i²(-1)^5=+1, so this target is real.
    let expected = -4.880_698_568_451e-6_f64;
    for (profile, sigma) in [("exponential", 1.0), ("poly_exponential", 1.3)] {
        cli.run_command(&format!(
            "set process kv h_function.function=\"{profile}\" h_function.sigma={sigma}"
        ))?;
        let mut errors = Vec::new();
        for samples in [20_000, 100_000] {
            cli.run_command(&format!("set process kv integrator.n_max={samples}"))?;
            let started = std::time::Instant::now();
            let output = Integrate {
                renderer: RendererOption::Tabled,
                show_max_weight_info: false,
                no_stream_iterations: true,
                no_stream_updates: true,
                ..default_integrate_for(name)
            }
            .run(&mut cli.state, &cli.cli_settings)?;
            let estimate = single_slot_integral(&output);
            let (value, error) = (estimate.result.re.0, estimate.error.re.0);
            info!(
                profile,
                sigma,
                samples = estimate.neval,
                value,
                error,
                expected,
                elapsed = ?started.elapsed(),
                "Non-aligned two-loop Fermi-surface integration convergence"
            );
            assert_eq!(estimate.neval, samples, "the sample budget must be fixed");
            assert!(
                [value, error, estimate.result.im.0, estimate.error.im.0]
                    .into_iter()
                    .all(f64::is_finite)
                    && value < 0.0
                    && error > 0.0
                    && estimate.error.im.0 >= 0.0,
                "{profile}, samples={samples}: {estimate}"
            );
            assert!(
                estimate.result.im.0.abs() < 1.0e-12 * expected.abs()
                    && estimate.error.im.0 < 1.0e-12 * expected.abs(),
                "the imaginary part must vanish: {estimate}"
            );
            assert!(
                (value - expected).abs() <= 5.0 * error,
                "{profile}, samples={samples}: {value:e} +/- {error:e}, target={expected:e}"
            );
            errors.push(error);
        }
        assert!(
            errors[1] < 0.01 * expected.abs() && errors[1] < 0.7 * errors[0],
            "{profile}: increased statistics must reduce uncertainty below 1%: {errors:?}"
        );
    }
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
fn v_diag() -> Result<()> {
    let cli = get_test_cli(
        Some("v_diag.toml".into()),
        get_tests_workspace_path().join("v_diag"),
        Some("v_diag".to_string()),
        true,
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

    cli.run_command("run import_graph")?;

    cli.run_command("run no_integrated")?;

    cli.run_command("generate")?;

    let integrate_command = Integrate {
        ..default_integrate_for("scalar_bubble")
    };

    let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
        ..Default::default()
    });

    // Kaapo reference magnitudes: m=1 muv=5 2.03838e-02 m=2 muv=5 	1.16050e-02	 m=3 muv=5 6.46968e-03

    cli.run_command("set model mass_scalar_1=1.0")?;
    let res = profile_cmd
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();
    assert_eq!(res.pass_fail(-0.9).failed, 0);
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-2.03838e-02)), 2));

    cli.run_command("set model mass_scalar_1=2.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-1.16050e-02)), 2));

    cli.run_command("set model mass_scalar_1=3.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-6.46968e-03)), 2));
    let renorm_command = Renormalize::default();

    let res = renorm_command.run(&mut cli.state, &cli.cli_settings)?;

    println!("{}", res[0]);

    clean_test(&cli.cli_settings.state.folder);

    Ok(())
}

mod failing {
    use super::*;

    #[test]
    fn test_grouped_subtraction() -> Result<()> {
        let mut cli = get_test_cli(
            Some("test_grouped_subtraction.toml".into()),
            get_tests_workspace_path().join("test_grouped_subtraction"),
            None,
            false,
        )?;

        let int1 = Integrate {
            process: vec![ProcessRef::Id(0)],
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
    fn vakint_compare_scalar_bubble() -> Result<()> {
        let mut cli = get_test_cli(
            Some("scalar_bubble.toml".into()),
            get_tests_workspace_path().join("scalar_bubble"),
            Some("scalar_bubble".to_string()),
            false,
        )?;

        cli.run_command("run import_graph")?;

        cli.run_command("run integrated")?;
        // These historical targets describe the local-only subtraction at m_uv=5.
        // Enabling the integrated addback instead gives the MS-bar result at mu_r;
        // its reference scale must be reconciled before this failing case is restored.

        cli.run_command("generate")?;

        let integrate_command = Integrate {
            ..default_integrate_for("scalar_bubble")
        };

        let profile_cmd = Profile::UltraViolet(UltraVioletProfile {
            ..Default::default()
        });

        // Kaapo reference magnitudes: m=1 muv=5 2.03838e-02 m=2 muv=5 	1.16050e-02	 m=3 muv=5 6.46968e-03

        cli.run_command("set model mass_scalar_1=1.0")?;
        let res = profile_cmd
            .run(&mut cli.state, &cli.cli_settings)?
            .unwrap_uv();
        assert_eq!(res.pass_fail(-0.9).failed, 0);
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-2.03838e-02)), 2)
        );

        cli.run_command("set model mass_scalar_1=2.0")?;
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-1.16050e-02)), 2)
        );

        cli.run_command("set model mass_scalar_1=3.0")?;
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-6.46968e-03)), 2)
        );
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
        // Four propagator i factors and three loops convert the historical
        // Euclidean targets below to the physical -i axis. Their subdivergence
        // prescription and p^2=1 magnitudes remain unverified; this phase migration
        // does not certify the old baselines for the current subtraction settings.
        cli.run_command("set model mass_scalar_1=1.0")?;
        let res = profile_cmd
            .run(&mut cli.state, &cli.cli_settings)?
            .unwrap_uv();
        assert_eq!(res.pass_fail(-0.9).failed, 0);
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-1.47240e-03)), 1),
            "Not compatible: {integral_no_cache}",
        );

        // cli.run_command("set model mass_scalar_1=2.0")?;
        // let res = profile_cmd.run(&mut cli.state, &cli.cli_settings)?;
        // assert_eq!(res.pass_fail(-0.9).failed, 0);
        // let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        // assert!(
        //     integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-7.15184e-04)), 3),
        //     "Not compatible: {integral_no_cache}",
        // );

        // cli.run_command("set model mass_scalar_1=3.0")?;
        // let res = profile_cmd.run(&mut cli.state, &cli.cli_settings)?;
        // assert_eq!(res.pass_fail(-0.9).failed, 0);
        // let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        // assert!(
        //     integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(-2.27485e-04)), 3),
        //     "Not compatible: {integral_no_cache}",
        // );

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
            ..Default::default()
        }
        .run(&mut cli)?;

        // At these subthreshold inputs the four-vertex/four-propagator scalar
        // box has positive imaginary phase. Preserve the historical magnitude
        // on that axis; it still needs an independent pointwise certificate.
        assert_snapshot!(format!("{a:.8e}"),@"(0e0+1.3485885914334373e-3i)");

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
}

mod slow {
    use super::*;

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

        // Kaapo magnitudes: m=1 muv=5 5.89551e-06 m=2 muv=5 	3.35645e-06	 m=3 muv=5 1.87120e-06
        // This primitive K4 vacuum graph has period 6*zeta(3), so local subtraction
        // gives i*6*zeta(3)*log(25/m^2)/(16*pi^2)^3. Vertex numerators are one;
        // six propagator i factors and three Wick rotations give +i.
        // Independent period: https://arxiv.org/abs/hep-ph/9504352.
        cli.run_command("set model mass_scalar_1=1.0")?;
        let res = profile_cmd
            .run(&mut cli.state, &cli.cli_settings)?
            .unwrap_uv();
        assert_eq!(res.pass_fail(-0.9).failed, 0);
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache
                .is_compatible_with_target(Complex::new(F(0.0), F(5.895509129899469e-06)), 1),
            "Not compatible: {integral_no_cache}",
        );

        cli.run_command("set model mass_scalar_1=2.0")?;
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache
                .is_compatible_with_target(Complex::new(F(0.0), F(3.3564515497441014e-06)), 3),
            "Not compatible: {integral_no_cache}",
        );

        cli.run_command("set model mass_scalar_1=3.0")?;
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache
                .is_compatible_with_target(Complex::new(F(0.0), F(1.8711980781814097e-06)), 3),
            "Not compatible: {integral_no_cache}",
        );

        clean_test(&cli.cli_settings.state.folder);

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
        // These historical targets have an unresolved input mismatch: the card
        // uses m_uv=1 and p^2=-12. Seven propagator i factors and four loops
        // convert the historical Euclidean targets to the physical +i axis;
        // this does not certify their signs or magnitudes for the current inputs.
        cli.run_command("set model mass_scalar_1=1.0")?;
        let res = profile_cmd
            .run(&mut cli.state, &cli.cli_settings)?
            .unwrap_uv();
        assert_eq!(res.pass_fail(-0.9).failed, 0);

        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(2.90078e-06)), 1),
            "Not compatible: {integral_no_cache}",
        );

        cli.run_command("set model mass_scalar_1=2.0")?;
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(1.59168e-06)), 3),
            "Not compatible: {integral_no_cache}",
        );

        cli.run_command("set model mass_scalar_1=3.0")?;
        let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
        assert!(
            integral_no_cache.is_compatible_with_target(Complex::new(F(0.0), F(6.86001e-07)), 3),
            "Not compatible: {integral_no_cache}",
        );

        clean_test(&cli.cli_settings.state.folder);

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
}
