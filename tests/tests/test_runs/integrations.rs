use super::utils::*;
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

    cli.run_command("run import_graph")?;

    cli.run_command("run no_integrated")?;

    cli.run_command("generate")?;

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
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(-2.03838e-02)), 1));

    cli.run_command("set model mass_scalar_1=2.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(-1.16050e-02)), 1));

    cli.run_command("set model mass_scalar_1=3.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(-6.46968e-03)), 1));
    let renorm_command = Renormalize::default();

    let res = renorm_command.run(&mut cli.state, &cli.cli_settings)?;

    println!("{}", res[0]);

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

    cli.run_command("generate")?;

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
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(-2.03838e-02)), 1));

    cli.run_command("set model mass_scalar_1=2.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(-1.16050e-02)), 1));

    cli.run_command("set model mass_scalar_1=3.0")?;
    let integral_no_cache = integrate_command.run(&mut cli.state, &cli.cli_settings)?;
    assert!(integral_no_cache.is_compatible_with_target(Complex::new_re(F(-6.46968e-03)), 1));
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
