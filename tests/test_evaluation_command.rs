use std::str::FromStr;

use color_eyre::Result;

mod test_utils;
use gammaloop_api::commands::Commands;
use gammalooprs::feyngen::diagram_generator::evaluate_overall_factor;
use gammalooprs::processes::ProcessCollection;
use symbolica::atom::{Atom, AtomCore};
use test_utils::{clean_test, get_test_cli, get_tests_workspace_path};
use tracing::debug;
use vakint::{vakint_symbol, NumericalEvaluationResult};

use crate::test_utils::CLIState;
use gammalooprs::utils::VAKINT;
use insta::assert_snapshot;
use which::which;

fn count_graphs_in_processes(cli: &CLIState) -> (usize, Atom) {
    assert_eq!(cli.state.process_list.processes.len(), 1);
    let process = &cli.state.process_list.processes[0];
    let integrand_names = process.collection.get_integrand_names();
    assert_eq!(integrand_names.len(), 1);
    let integrand_name = integrand_names[0];

    let graphs = match &process.collection {
        ProcessCollection::Amplitudes(amps) => {
            let amp = amps.get(integrand_name).unwrap();
            amp.graphs.iter().map(|g| &g.graph).collect::<Vec<_>>()
        }
        ProcessCollection::CrossSections(xss) => {
            let xs = xss.get(integrand_name).unwrap();
            xs.supergraphs.iter().map(|g| &g.graph).collect::<Vec<_>>()
        }
    };
    let n_graphs = graphs.len();
    let mut overall_factor_sum = Atom::Zero;
    for g in graphs {
        overall_factor_sum += &g.overall_factor;
    }
    overall_factor_sum = evaluate_overall_factor(overall_factor_sum.as_view());
    (n_graphs, overall_factor_sum)
}

fn split_before_flags(s: &str) -> (&str, &str) {
    match s.find("--") {
        Some(i) => {
            let (left, right) = s.split_at(i);
            (left.trim_end(), right) // X, Y
        }
        None => (s, ""), // no flags present
    }
}

fn generate_graphs_and_count(
    cli: &mut CLIState,
    generation_type: &str,
    process: &str,
) -> Result<(usize, Atom)> {
    let (process_input, process_options) = split_before_flags(process);
    let cmd = format!(
        "generate {} {} {} --clear-existing-processes --only-diagrams",
        generation_type, process_input, process_options
    );
    debug!("Running diagram generation command: {}", cmd);
    cli.run_command(&cmd)?;
    Ok(count_graphs_in_processes(cli))
}

fn feyngen_str(cli: &mut CLIState, generation_type: &str, process: &str) -> Result<String> {
    let (n_graphs, overall_factor_sum) = generate_graphs_and_count(cli, generation_type, process)?;
    Ok(format!(
        "{} | {}",
        n_graphs,
        overall_factor_sum.to_canonical_string()
    ))
}

fn compare_numerical_evaluation(
    cli: &mut CLIState,
    cmd: &str,
    target: &[(i64, (String, String))],
) -> Result<()> {
    let cmd = Commands::from_str(cmd)?;
    let eval = if let Commands::Evaluate(eval_cmd) = cmd {
        eval_cmd.run(
            &mut cli.state,
            &cli.cli_settings,
            &cli.default_runtime_settings,
        )?
    } else {
        panic!("Expected Evaluate command");
    };
    let vakint_settings = {
        let lock = VAKINT.read().unwrap();
        lock.as_ref()
            .unwrap_or_else(|| {
                panic!("VAKINT not initialized in CLI yet");
            })
            .settings
            .clone()
    };
    let num_eval = NumericalEvaluationResult::from_atom(
        eval.as_view(),
        vakint_symbol!(&vakint_settings.epsilon_symbol),
        &vakint_settings,
    )?;
    // println!("Numerical evaluation result:\n{}", num_eval);

    let target_eval = NumericalEvaluationResult::from_vec(target.into(), &vakint_settings);

    let (matches, match_msg) = target_eval.does_approx_match(
        &num_eval,
        None,
        10.0_f64.powi(-((vakint_settings.run_time_decimal_precision - 4) as i32)),
        1.0,
    );
    assert!(
        matches,
        "Numerical evaluation does not match target result:\n{}\nvs target\n{}\n{}",
        num_eval, target_eval, match_msg
    );

    Ok(())
}

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
    cli.run_command("set default-runtime kv general.mu_r_sq=12.0")?;
    cli.run_command("set model mass_scalar_1={re:1.0,im:0.0}")?;

    #[rustfmt::skip]
    cli.run_command(&format!(
        "{} {}",
        "set global kv",
        [ 
            ("run_time_decimal_precision", "100"),
            ("evaluation_methods", r#"["alphaloop","matad"]"#),
            ("form_exe_path", form_exe_path.to_str().unwrap(),)
        ]
        .iter()
        .map(|(k, v)| format!("global.generation.uv.vakint.{k}={v}"))
        .collect::<Vec<_>>()
        .join(" ")
    ))?;

    #[rustfmt::skip]
    compare_numerical_evaluation(
        &mut cli,
        "evaluate --numerical --n-epsilon-terms 4",
        &[
            (0, ("0.0".into(),"-3.28986813369645287294483033329205037843789980241359687547111645874001494080640174766725780123951741".into())),
            (1, ("0.0".into(),"5.56266525336352304035680257963623222447793296679948450388669066192714938626082105298641985966817556".into())),
            (2, ("0.0".into(),"-6.99738383886178490341873784327636031576356895962030139357161896913681610683215425206713129095599914".into())),
        ],
    )?;

    let cmd = Commands::from_str(&"evaluate --n-epsilon-terms 4")?;
    let eval = if let Commands::Evaluate(eval_cmd) = cmd {
        eval_cmd.run(
            &mut cli.state,
            &cli.cli_settings,
            &cli.default_runtime_settings,
        )?
    } else {
        panic!("Expected Evaluate command");
    };
    #[rustfmt::skip]
    assert_snapshot!(eval.to_canonical_string(),@"(((-1*log(𝜋)*vakint::𝑖*𝜋^2+log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::𝑖*𝜋^2)*-1/3*UFO::mass_scalar_1^(-2)+1/6*UFO::mass_scalar_1^(-2)*vakint::𝑖*𝜋^2)*-2*log(UFO::mass_scalar_1)+((log(1/4*vakint::μᵣ²*𝜋^(-1)))^2*1/2*vakint::𝑖*𝜋^2+(log(𝜋))^2*1/2*vakint::𝑖*𝜋^2+-1*log(1/4*vakint::μᵣ²*𝜋^(-1))*log(𝜋)*vakint::𝑖*𝜋^2)*-1/3*UFO::mass_scalar_1^(-2)+(-1*log(𝜋)*vakint::𝑖*𝜋^2+log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::𝑖*𝜋^2)*1/6*UFO::mass_scalar_1^(-2)+(log(UFO::mass_scalar_1))^2*-2/3*UFO::mass_scalar_1^(-2)*vakint::𝑖*𝜋^2+-1/36*UFO::mass_scalar_1^(-2)*vakint::𝑖*𝜋^4)*vakint::ε^2+((-1*log(𝜋)*vakint::𝑖*𝜋^2+log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::𝑖*𝜋^2)*-1/3*UFO::mass_scalar_1^(-2)+1/6*UFO::mass_scalar_1^(-2)*vakint::𝑖*𝜋^2+2/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::𝑖*𝜋^2)*vakint::ε+-1/3*UFO::mass_scalar_1^(-2)*vakint::𝑖*𝜋^2");

    #[rustfmt::skip]
    compare_numerical_evaluation(
        &mut cli,
        "evaluate --numerical --n-epsilon-terms 5",
        &[
            (0, ("0.0".into(),"-3.28986813369645287294483033329205037843789980241359687547111645874001494080640174766725780123951741".into())),
            (1, ("0.0".into(),"5.56266525336352304035680257963623222447793296679948450388669066192714938626082105298641985966817556".into())),
            (2, ("0.0".into(),"-6.99738383886178490341873784327636031576356895962030139357161896913681610683215425206713129095599914".into())),
            (3, ("0.0".into(),"7.98563411114514294560056936573344248116767960794442757461392872305402199431688012106743992195327359".into())),
        ],
    )?;

    cli.run_command("set global kv global.generation.uv.vakint.matad.substitute_hpls=false global.generation.uv.vakint.matad.direct_numerical_substition=false")?;

    let cmd = Commands::from_str(&"evaluate --n-epsilon-terms 5")?;
    let eval = if let Commands::Evaluate(eval_cmd) = cmd {
        eval_cmd.run(
            &mut cli.state,
            &cli.cli_settings,
            &cli.default_runtime_settings,
        )?
    } else {
        panic!("Expected Evaluate command");
    };
    #[rustfmt::skip]
    assert_snapshot!(eval.to_canonical_string(),@"((log(1/4*vakint::μᵣ²*𝜋^(-1)))^2*-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(UFO::mass_scalar_1))^2*-2/3*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(𝜋))^2*-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(vakint::PolyGamma(0,1))^2*-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+-1/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/3*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^2*vakint::𝑖*𝜋^2+-1/3*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε*vakint::𝑖*𝜋^2+-1/6*UFO::mass_scalar_1^(-2)*log(𝜋)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(1,1)*vakint::ε^3*vakint::𝑖*𝜋^2+-2/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*log(𝜋)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*log(𝜋)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/3*UFO::mass_scalar_1^(-2)*log(𝜋)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/3*UFO::mass_scalar_1^(-2)*log(𝜋)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+1/6*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+2/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+2/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+2/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2)*vakint::EulerGamma+(-1/6*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+1/12*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/6*UFO::mass_scalar_1^(-2)*log(𝜋)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2)*vakint::EulerGamma^2+(log(1/4*vakint::μᵣ²*𝜋^(-1)))^2*-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(1/4*vakint::μᵣ²*𝜋^(-1)))^2*-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+(log(1/4*vakint::μᵣ²*𝜋^(-1)))^2*1/12*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(1/4*vakint::μᵣ²*𝜋^(-1)))^2*1/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(1/4*vakint::μᵣ²*𝜋^(-1)))^2*1/6*UFO::mass_scalar_1^(-2)*log(𝜋)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(1/4*vakint::μᵣ²*𝜋^(-1)))^3*-1/18*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(UFO::mass_scalar_1))^2*-2/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(UFO::mass_scalar_1))^2*-2/3*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(UFO::mass_scalar_1))^2*-2/3*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+(log(UFO::mass_scalar_1))^2*1/3*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(UFO::mass_scalar_1))^2*2/3*UFO::mass_scalar_1^(-2)*log(𝜋)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(UFO::mass_scalar_1))^3*4/9*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(𝜋))^2*-1/6*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(𝜋))^2*-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(𝜋))^2*-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+(log(𝜋))^2*1/12*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(𝜋))^2*1/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(log(𝜋))^3*1/18*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(vakint::PolyGamma(0,1))^2*-1/6*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(vakint::PolyGamma(0,1))^2*-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+(vakint::PolyGamma(0,1))^2*1/12*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(vakint::PolyGamma(0,1))^2*1/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(vakint::PolyGamma(0,1))^2*1/6*UFO::mass_scalar_1^(-2)*log(𝜋)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+(vakint::PolyGamma(0,1))^3*-1/18*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/18*UFO::mass_scalar_1^(-2)*vakint::EulerGamma^3*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/18*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(2,1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^2*vakint::𝑖*𝜋^2+-1/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::Gamma(1)*vakint::ε*vakint::𝑖*𝜋^2+-1/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+-1/3*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε*vakint::𝑖*𝜋^2+-1/3*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::𝑖*𝜋^2+-1/6*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*log(𝜋)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/6*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::Gamma(1)*vakint::PolyGamma(1,1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/6*UFO::mass_scalar_1^(-2)*log(𝜋)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/6*UFO::mass_scalar_1^(-2)*log(𝜋)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::PolyGamma(1,1)*vakint::ε^3*vakint::𝑖*𝜋^2+-1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(1,1)*vakint::ε^2*vakint::𝑖*𝜋^2+-2/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*log(UFO::mass_scalar_1)*log(𝜋)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+-2/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*log(𝜋)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+-2/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*log(𝜋)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+1/12*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(1,1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*log(𝜋)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*log(𝜋)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+1/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*log(𝜋)*vakint::Gamma(1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::PolyGamma(1,1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/3*UFO::mass_scalar_1^(-2)*log(𝜋)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^2*vakint::𝑖*𝜋^2+1/3*UFO::mass_scalar_1^(-2)*log(𝜋)*vakint::Gamma(1)*vakint::ε*vakint::𝑖*𝜋^2+1/6*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/6*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+1/6*UFO::mass_scalar_1^(-2)*log(𝜋)*vakint::Gamma(1)*vakint::PolyGamma(1,1)*vakint::ε^3*vakint::𝑖*𝜋^2+1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^2*vakint::𝑖*𝜋^2+1/6*UFO::mass_scalar_1^(-2)*vakint::Gamma(1)*vakint::ε*vakint::𝑖*𝜋^2+2/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^3*vakint::𝑖*𝜋^2+2/3*UFO::mass_scalar_1^(-2)*log(1/4*vakint::μᵣ²*𝜋^(-1))*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::ε^2*vakint::𝑖*𝜋^2+2/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::PolyGamma(0,1)*vakint::ε^2*vakint::𝑖*𝜋^2+2/3*UFO::mass_scalar_1^(-2)*log(UFO::mass_scalar_1)*vakint::Gamma(1)*vakint::ε*vakint::𝑖*𝜋^2");

    clean_test(&cli.cli_settings.state_folder);

    Ok(())
}
