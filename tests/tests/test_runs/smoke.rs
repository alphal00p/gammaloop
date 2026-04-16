use super::utils::*;
use super::*;

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

mod slow {
    use super::*;

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
            workspace_path: Some(
                get_tests_workspace_path().join("epem_ttxh_lo/integration_workspace"),
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
}

mod failing {
    use super::*;

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
}
