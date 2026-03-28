use super::utils::*;
use super::*;

#[test]
fn inspect_x_space_reports_invalid_coordinate_count_cleanly() -> Result<()> {
    let mut cli = get_test_cli(
        Some("z_decay_test.toml".into()),
        get_tests_workspace_path().join("z_decay_invalid_inspect_point"),
        None,
        true,
    )?;

    let error = Inspect {
        process: Some(ProcessRef::Id(0)),
        integrand_name: Some("default".to_string()),
        point: vec![0.1, 0.2],
        momentum_space: false,
        discrete_dim: vec![0, 0],
        ..Default::default()
    }
    .run(&mut cli)
    .expect_err("inspect should reject an invalid x-space point length");

    let rendered = format!("{error:#}");
    assert!(
        rendered.contains("Expected 3 x-space coordinates for this integrand selection, got 2."),
        "{rendered}"
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn inspect_x_space_reports_missing_discrete_dimensions_cleanly() -> Result<()> {
    let mut cli =
        setup_gg_hhh_threshold_amplitude_cli("gg_hhh_inspect_missing_discrete_dimensions")?;

    let error = Inspect {
        process: Some(ProcessRef::Unqualified("gg_hhh".to_string())),
        integrand_name: Some("1L".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: false,
        ..Default::default()
    }
    .run(&mut cli)
    .expect_err("inspect should reject x-space points that omit required discrete dimensions");

    let rendered = format!("{error:#}");
    assert!(
        rendered.contains("This integrand uses discrete graph sampling")
            && rendered.contains("sample_orientations = true")
            && rendered.contains("sampling_type = discrete_multi_channeling")
            && rendered.contains(
                "requires 3 discrete dimensions [graph group, orientation, channel], but got 0."
            ),
        "{rendered}"
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
fn inspect_momentum_space_graph_id_is_sampling_agnostic() -> Result<()> {
    let mut cli = setup_scalar_topologies_cli("scalar_box_graph_id_inspect")?;

    cli.run_command(
        r#"set process -p box -i scalar_box string '
[sampling]
type = "discrete_graph_sampling"
[sampling.sampling_type]
subtype = "multi_channeling"
'"#,
    )?;

    let (_, graph_selected) = Inspect {
        process: Some(ProcessRef::Unqualified("box".to_string())),
        integrand_name: Some("scalar_box".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        graph_id: Some(0),
        ..Default::default()
    }
    .run(&mut cli)?;

    let (_, discrete_selected) = Inspect {
        process: Some(ProcessRef::Unqualified("box".to_string())),
        integrand_name: Some("scalar_box".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        discrete_dim: vec![0],
        ..Default::default()
    }
    .run(&mut cli)?;

    assert_eq!(graph_selected, discrete_selected);

    clean_test(&cli.cli_settings.state.folder);
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
