use color_eyre::Result;

use gammaloop_api::{inspect::Inspect, integrate::Integrate};

use gammalooprs::{
    initialisation::initialise,
    model::UFOSymbol,
    numerator::ufo::UFOSymbols,
    status_info,
    utils::{test_utils::load_generic_model, F},
};
use insta::assert_snapshot;
use momtrop::assert_approx_eq;
use spenso::algebra::complex::Complex;
use std::{
    env,
    ops::{ControlFlow, Deref, DerefMut},
    path::{Path, PathBuf},
    str::FromStr,
};
use symbolica::symbol;
use tracing::{debug, warn};

mod test_utils;
use test_utils::{clean_test, get_test_cli, get_tests_workspace_path};

#[test]
fn qqx_aaa_subtracted_nlo_amplitude_test() -> Result<()> {
    let mut state = get_test_cli(
        Some("qqx_aaa_subtracted_nlo_amplitude.toml".into()),
        get_tests_workspace_path().join("qqx_aaa_subtracted_nlo_amplitude"),
        None,
    )?;

    let a = Inspect {
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        ..Default::default()
    }
    .run(&mut state)?;

    assert_snapshot!(format!("{a:.8e}"),@"(6.923469780323106e-4+-6.181987787694587e-4i)");

    Ok(())
}

#[test]
fn trees() -> Result<()> {
    let mut cli = get_test_cli(
        Some("trees/qqx_aaa.toml".into()),
        get_tests_workspace_path().join("qqx_aaa_tree"),
        None,
    )?;

    let a = Inspect {
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        ..Default::default()
    }
    .run(&mut cli)?;

    assert_snapshot!(format!("{a:.8e}"),@"(1.4727604164105595e-4+-1.1503139369130214e-3i)");

    clean_test(&cli);
    Ok(())
}

#[test]
fn photons_1l() -> Result<()> {
    let mut cli = get_test_cli(Some("photons.toml".into()), "./tests/photons", None)?;

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

    //let inspect = Inspect {
    //    process_id: 0,
    //    integrand_name: "physical_1L_6photons".to_string(),
    //    point: vec![0.1, 0.2, 0.3],
    //    momentum_space: true,
    //    ..Default::default()
    //}
    //.run(&mut cli)?;

    //assert_snapshot!(format!("{inspect:.8e}"),@"(-1.1779656944253678e-15+9.093348564643939e-16i)");

    let target = Complex::new(F(-1.22898408452706e-13), F(-3.94362534040412e-13));

    let integrate = Integrate {
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some("./tests/photons/integration_workspace/integration_results.yaml".into()),
        workspace_path: Some("./tests/photons/integration_workspace/".into()),
        n_cores: 6,
        target: Some(vec![target.re.0, target.im.0]),
        restart: true,
    };

    let integration_result = integrate.run(&mut cli)?;

    assert_approx_eq(&integration_result.result.re, &target.re, &F(10e-2));
    assert_approx_eq(&integration_result.result.im, &target.im, &F(10e-2));
    Ok(())
}

#[test]
fn test_grouped_subtraction() -> Result<()> {
    let mut cli = get_test_cli(
        Some("test_grouped_subtraction.toml".into()),
        get_tests_workspace_path().join("test_grouped_subtraction"),
        None,
    )?;

    let int1 = Integrate {
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(get_tests_workspace_path().join(
            "test_grouped_subtraction/integration_workspace_no_group/integration_results.yaml",
        )),
        n_cores: 1,
        workspace_path: Some(
            get_tests_workspace_path()
                .join("test_grouped_subtraction/integration_workspace_no_group/"),
        ),
        target: None,
        restart: true,
    };
    let int2 =
        Integrate {
            process_id: Some(1),
            integrand_name: Some("default".to_string()),
            result_path: Some(get_tests_workspace_path().join(
                "test_grouped_subtraction/integration_workspace_group/integration_results.yaml",
            )),
            n_cores: 1,
            workspace_path: Some(
                get_tests_workspace_path()
                    .join("test_grouped_subtraction/integration_workspace_group/"),
            ),
            target: None,
            restart: true,
        };

    let integration_results_no_group = int1.run(&mut cli)?;
    let integration_results_group = int2.run(&mut cli)?;

    status_info!("No group result: {:#?}", integration_results_no_group);
    status_info!("Group result: {:#?}", integration_results_group);

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

    clean_test(&cli);

    Ok(())
}

#[test]
fn scalar_box() -> Result<()> {
    let mut cli = get_test_cli(
        Some("scalar_box.toml".into()),
        get_tests_workspace_path().join("scalar_box"),
        Some("scalar_box".to_string()),
    )?;

    let a = Inspect {
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        point: vec![0.1, 0.2, 0.3],
        momentum_space: true,
        discrete_dim: vec![0],
        ..Default::default()
    }
    .run(&mut cli)?;

    assert_snapshot!(format!("{a:.8e}"),@"(1.3485885914334373e-3+0e0i)");

    cli.run_command("set process -p 0 -i default kv General.enable_cache=false")?;
    let integral_no_cache = Integrate {
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(
            get_tests_workspace_path()
                .join("scalar_box/integration_workspace/integration_results.toml"),
        ),
        n_cores: 1,
        workspace_path: Some(get_tests_workspace_path().join("scalar_box/integration_workspace")),
        target: None,
        restart: true,
    }
    .run(&mut cli)?;

    cli.run_command("set process -p 0 -i default kv General.enable_cache=true")?;
    let integral_with_cache = Integrate {
        process_id: Some(0),
        integrand_name: Some("default".to_string()),
        result_path: Some(
            get_tests_workspace_path()
                .join("scalar_box/integration_workspace/integration_results.toml"),
        ),
        n_cores: 1,
        workspace_path: Some(get_tests_workspace_path().join("scalar_box/integration_workspace")),
        target: None,
        restart: true,
    }
    .run(&mut cli)?;

    status_info!("Integral result without caching: {:#?}", integral_no_cache);
    status_info!(
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

    clean_test(&cli);

    Ok(())
}
