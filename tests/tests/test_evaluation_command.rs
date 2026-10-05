use std::{collections::HashMap, str::FromStr};

use color_eyre::Result;
use gammaloop_api::commands::Commands;
use gammaloop_integration_tests::{clean_test, get_test_cli, get_tests_workspace_path};
use gammalooprs::utils::{GS, vakint};
use symbolica::{
    atom::{Atom, AtomCore},
    parse,
    poly::series::SeriesDepth,
    symbol,
};
use vakint::NumericalEvaluationResult;
use which::which;

#[test]
fn evaluate_1l_scalar_vacuum() -> Result<()> {
    let mut cli = get_test_cli(
        Some("scalars_load.toml".into()),
        get_tests_workspace_path().join("evaluation_1L_scalar_vacuum"),
        Some("evaluation_1L_scalar_vacuum".to_string()),
        true,
    )?;
    let form_exe_path = which("form").expect("FORM executable could not be located.");
    cli.run_command("import graphs ./tests/resources/graphs/1l_vacuum.dot")?;
    cli.run_command(&format!(
        "set global kv global.generation.uv.vakint.run_time_decimal_precision=100 global.generation.uv.vakint.evaluation_methods=[\"alphaloop\",\"matad\"] global.generation.uv.vakint.form_exe_path={}",
        form_exe_path.display()
    ))?;
    // This fixture's reference uses the historical measure without (2*pi)^(-D).
    // Specify it rather than inheriting the current MSbar default.
    cli.run_command(
        "set global kv global.generation.uv.vakint.normalization='exp(eps*(log_mu_sq-log(4*pi)+EulerGamma))'",
    )?;

    // Independently expand the standard one-loop formula for q^2/(q^2-m^2)^4:
    // -i*pi^2/(3*m^2)*(1-eps/2)*exp(L*eps)*Gamma(1+eps)*exp(EulerGamma*eps),
    // where L=log(mu_r^2/(4*pi^2*m^2)). The coefficients retain the original
    // analytic reference without depending on Symbolica's printed ordering.
    let mass = symbol!("evaluation_reference_mass");
    let scale_squared = symbol!("evaluation_reference_scale_squared");
    let log_scale = symbol!("evaluation_reference_log_scale");
    let epsilon = symbol!("evaluation_reference_epsilon");
    let reference = parse!(
        "-1i*pi^2/(3*evaluation_reference_mass^2)*(1+(evaluation_reference_log_scale-1/2)*evaluation_reference_epsilon+(evaluation_reference_log_scale^2/2-evaluation_reference_log_scale/2+pi^2/12)*evaluation_reference_epsilon^2+(evaluation_reference_log_scale^3/6-evaluation_reference_log_scale^2/4+pi^2*evaluation_reference_log_scale/12-pi^2/24-zeta(3)/3)*evaluation_reference_epsilon^3)"
    )
    .replace(log_scale)
    .with(parse!(
        "log(evaluation_reference_scale_squared)-2*log(2)-2*log(pi)-2*log(evaluation_reference_mass)"
    ))
    .replace(epsilon)
    .with(GS.dim_epsilon);
    let settings = cli.cli_settings.global.generation.uv.vakint.true_settings();
    let vakint = vakint()?;
    let historical_reference = reference
        .replace(mass)
        .with(Atom::one())
        .replace(scale_squared)
        .with(Atom::num(12));
    let (historical_numerical, _) = vakint.numerical_evaluation(
        &settings,
        historical_reference.as_view(),
        &HashMap::default(),
        &HashMap::default(),
        None,
    )?;
    // The old fixture used mu_r^2=12. Keep its 100-digit coefficients as a
    // separate certificate; general.mu_r now stores mu_r, not its square.
    #[rustfmt::skip]
    let historical_target = NumericalEvaluationResult::from_vec(vec![
        (0, ("0.0".into(), "-3.28986813369645287294483033329205037843789980241359687547111645874001494080640174766725780123951741".into())),
        (1, ("0.0".into(), "5.56266525336352304035680257963623222447793296679948450388669066192714938626082105298641985966817556".into())),
        (2, ("0.0".into(), "-6.99738383886178490341873784327636031576356895962030139357161896913681610683215425206713129095599914".into())),
        (3, ("0.0".into(), "7.98563411114514294560056936573344248116767960794442757461392872305402199431688012106743992195327359".into())),
    ], &settings);
    let tolerance = 10.0_f64.powi(-((settings.run_time_decimal_precision - 4) as i32));
    let (matches, message) =
        historical_target.does_approx_match(&historical_numerical, None, tolerance, 1.0);
    assert!(matches, "historical analytic reference differs: {message}");

    for mass_value in [1, 2] {
        cli.run_command(&format!("set model mass_scalar_1={mass_value}.0"))?;
        for scale in [1, 12] {
            cli.run_command(&format!("set default-runtime kv general.mu_r={scale}.0"))?;
            for terms in [None, Some(4), Some(5)] {
                let degree = terms.unwrap_or(2) - 2;
                let expected = reference
                    .replace(mass)
                    .with(Atom::num(mass_value))
                    .replace(scale_squared)
                    .with(Atom::num(scale * scale))
                    .series(GS.dim_epsilon, 0, SeriesDepth::absolute(degree))?
                    .to_atom();
                let (expected_numerical, _) = vakint.numerical_evaluation(
                    &settings,
                    expected.as_view(),
                    &HashMap::default(),
                    &HashMap::default(),
                    None,
                )?;
                for numerical in [false, true] {
                    let command = format!(
                        "evaluate{}{}",
                        if numerical { " --numerical" } else { "" },
                        terms.map_or_else(String::new, |terms| format!(
                            " --n-epsilon-terms {terms}"
                        ))
                    );
                    let Commands::Evaluate(evaluate) = Commands::from_str(&command)? else {
                        unreachable!();
                    };
                    let evaluation = evaluate.run(
                        &mut cli.state,
                        &cli.cli_settings,
                        &cli.default_runtime_settings,
                    )?;
                    let numerical_result = if numerical {
                        NumericalEvaluationResult::from_atom(
                            evaluation.as_view(),
                            GS.dim_epsilon,
                            &settings,
                        )?
                    } else {
                        let physical = evaluation
                            .replace(symbol!("UFO::mass_scalar_1"))
                            .with(Atom::num(mass_value))
                            .replace(GS.mu_r_sq)
                            .with(Atom::num(scale * scale));
                        vakint
                            .numerical_evaluation(
                                &settings,
                                physical.as_view(),
                                &HashMap::default(),
                                &HashMap::default(),
                                None,
                            )?
                            .0
                    };
                    assert_eq!(
                        numerical_result.get_epsilon_coefficients().len(),
                        degree as usize + 1
                    );
                    let (matches, message) = expected_numerical.does_approx_match(
                        &numerical_result,
                        None,
                        tolerance,
                        1.0,
                    );
                    assert!(
                        matches,
                        "m={mass_value}, mu_r={scale}, {command}: {message}"
                    );
                }
            }
        }
    }
    // Keep the symbolic MATAD branch: its special functions must encode the
    // same coefficients before their numerical substitution is enabled.
    cli.run_command("set global kv global.generation.uv.vakint.matad.substitute_hpls=false global.generation.uv.vakint.matad.direct_numerical_substition=false")?;
    let Commands::Evaluate(evaluate) = Commands::from_str("evaluate --n-epsilon-terms 5")? else {
        unreachable!();
    };
    let evaluation = evaluate.run(
        &mut cli.state,
        &cli.cli_settings,
        &cli.default_runtime_settings,
    )?;
    let physical = evaluation
        .replace(symbol!("UFO::mass_scalar_1"))
        .with(Atom::num(2))
        .replace(GS.mu_r_sq)
        .with(Atom::num(144));
    let physical = vakint::matad::MATAD::with_settings(settings.clone())
        .substitute_poly_gamma(physical.as_view())?;
    let (numerical_result, _) = vakint.numerical_evaluation(
        &settings,
        physical.as_view(),
        &HashMap::default(),
        &HashMap::default(),
        None,
    )?;
    let expected = reference
        .replace(mass)
        .with(Atom::num(2))
        .replace(scale_squared)
        .with(Atom::num(144));
    let (expected_numerical, _) = vakint.numerical_evaluation(
        &settings,
        expected.as_view(),
        &HashMap::default(),
        &HashMap::default(),
        None,
    )?;
    let (matches, message) =
        expected_numerical.does_approx_match(&numerical_result, None, tolerance, 1.0);
    assert!(matches, "symbolic MATAD coefficients differ: {message}");
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
