use super::*;
use gammaloop_api::commands::evaluate_samples::EvaluateSamplesPrecise;
use gammalooprs::{
    integrands::evaluation::PreciseEvaluationResultOutput,
    utils::{ArbPrec, F},
};
use ndarray::Array2;

const LOOP_X_POINT: [f64; 3] = [0.17, 0.31, 0.53];
const LOOP_MOMENTUM_POINT: [f64; 3] = [0.23, -0.17, 0.41];
const LOOP_MOMENTUM_POINT_Y_FLIPPED: [f64; 3] = [0.23, 0.17, 0.41];
const LOOP_MOMENTUM_POINT_NEGATED: [f64; 3] = [-0.23, 0.17, -0.41];

fn setup_adad_probe_cli(test_name: &str) -> Result<gammaloop_integration_tests::CLIState> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join(test_name),
        Some(test_name.to_string()),
        true,
    )?;
    run_commands(
        &mut cli,
        &[
            "import model sm-default",
            "set global kv global.generation.evaluator.iterative_orientation_optimization=false global.generation.evaluator.store_atom=true global.generation.evaluator.compile=false global.generation.evaluator.summed=false global.generation.evaluator.summed_function_map=true",
            "set global kv global.generation.threshold_subtraction.enable_thresholds=true global.generation.threshold_subtraction.check_esurface_at_generation=false",
            "set global kv global.generation.uv.generate_integrated=true",
            "import graphs ./tests/resources/graphs/uv_tests/ad_ad_1L_gluon.dot -p adad -i adad_gluon",
            "generate existing -p adad -i adad_gluon",
            r#"generate amp a d > a d [{0}] -p adad_tree -i adad_tree --global-prefactor-num "spenso::g(spenso::dind(spenso::cof(3,gammalooprs::hedge(1))),spenso::cof(3,gammalooprs::hedge(3)))""#,
        ],
    )?;
    cli.run_command(
        r#"set default-runtime string '
[general]
evaluator_method = "SummedFunctionMap"
enable_cache = false
debug_cache = false
generate_events = false
store_additional_weights_in_event = false
integral_unit = "none"
mu_r = 20.0
m_uv = 20.0

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
  [500.0, 0.0, 0.0, 500.0],
  [500.0, 0.0, 0.0, -500.0],
  [500.0, 500.0, 0.0, 0.0],
  "dependent"
]
helicities = [+1, +1, -1, -1]
'"#,
    )?;
    Ok(cli)
}

fn precise_probe(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
    point: &[f64],
    momentum_space: bool,
    graph_name: Option<&str>,
    orientation: Option<usize>,
) -> Result<PreciseEvaluationResultOutput> {
    let (process_id, integrand_name) = cli.state.find_integrand_ref(
        Some(&ProcessRef::Unqualified(process.to_string())),
        Some(&integrand.to_string()),
    )?;
    let points = Array2::from_shape_vec((1, point.len()), point.to_vec())?;
    Ok(EvaluateSamplesPrecise {
        process_id: Some(process_id),
        integrand_name: Some(integrand_name),
        use_arb_prec: true,
        minimal_output: false,
        return_generated_events: Some(true),
        momentum_space,
        points: points.view(),
        integrator_weights: None,
        discrete_dims: None,
        graph_names: graph_name.map(|name| vec![Some(name.to_string())]),
        orientations: orientation.map(|value| vec![Some(value)]),
    }
    .run(&mut cli.state)?
    .samples
    .into_iter()
    .next()
    .expect("expected a single precise sample")
    .evaluation)
}

fn arb_displayed_result(
    result: PreciseEvaluationResultOutput,
    momentum_space: bool,
) -> Result<Complex<F<ArbPrec>>> {
    match result {
        PreciseEvaluationResultOutput::Arb(result) => {
            if momentum_space {
                Ok(result.integrand_result)
            } else {
                let jacobian = result.parameterization_jacobian.ok_or_else(|| {
                    eyre::eyre!("x-space precise probe requires a parameterization jacobian")
                })?;
                Ok(result
                    .integrand_result
                    .map(|entry| entry * jacobian.clone()))
            }
        }
        _ => Err(eyre::eyre!("expected an Arb precise evaluation result")),
    }
}

#[test]
#[ignore = "investigation helper"]
#[serial]
fn adad_gluon_zero_probe() -> Result<()> {
    let mut cli = setup_adad_probe_cli("adad_gluon_zero_probe")?;

    {
        let loop_integrand = cli.state.process_list.get_integrand_mut(
            cli.state
                .find_integrand_ref(
                    Some(&ProcessRef::Unqualified("adad".to_string())),
                    Some(&"adad_gluon".to_string()),
                )?
                .0,
            "adad_gluon",
        )?;
        println!(
            "adad_gluon: graphs={}, graph0_name={:?}, graph0_orientations={:?}",
            loop_integrand.graph_count(),
            loop_integrand.graph_name_by_id(0),
            loop_integrand.graph_orientation_count(0)
        );
    }

    let summed_x = arb_displayed_result(
        precise_probe(
            &mut cli,
            "adad",
            "adad_gluon",
            &LOOP_X_POINT,
            false,
            None,
            None,
        )?,
        false,
    )?;
    println!("summed x-space @ {:?} = {}", LOOP_X_POINT, summed_x);

    let orientation0 = arb_displayed_result(
        precise_probe(
            &mut cli,
            "adad",
            "adad_gluon",
            &LOOP_MOMENTUM_POINT,
            true,
            Some("adad"),
            Some(0),
        )?,
        true,
    )?;
    println!(
        "graph adad orientation 0 momentum-space @ {:?} = {}",
        LOOP_MOMENTUM_POINT, orientation0
    );

    let orientation0_y_flipped = arb_displayed_result(
        precise_probe(
            &mut cli,
            "adad",
            "adad_gluon",
            &LOOP_MOMENTUM_POINT_Y_FLIPPED,
            true,
            Some("adad"),
            Some(0),
        )?,
        true,
    )?;
    println!(
        "graph adad orientation 0 momentum-space @ {:?} = {}",
        LOOP_MOMENTUM_POINT_Y_FLIPPED, orientation0_y_flipped
    );

    let orientation1 = arb_displayed_result(
        precise_probe(
            &mut cli,
            "adad",
            "adad_gluon",
            &LOOP_MOMENTUM_POINT,
            true,
            Some("adad"),
            Some(1),
        )?,
        true,
    )?;
    println!(
        "graph adad orientation 1 momentum-space @ {:?} = {}",
        LOOP_MOMENTUM_POINT, orientation1
    );

    let orientation1_y_flipped = arb_displayed_result(
        precise_probe(
            &mut cli,
            "adad",
            "adad_gluon",
            &LOOP_MOMENTUM_POINT_Y_FLIPPED,
            true,
            Some("adad"),
            Some(1),
        )?,
        true,
    )?;
    println!(
        "graph adad orientation 1 momentum-space @ {:?} = {}",
        LOOP_MOMENTUM_POINT_Y_FLIPPED, orientation1_y_flipped
    );

    let orientation1_negated = arb_displayed_result(
        precise_probe(
            &mut cli,
            "adad",
            "adad_gluon",
            &LOOP_MOMENTUM_POINT_NEGATED,
            true,
            Some("adad"),
            Some(1),
        )?,
        true,
    )?;
    println!(
        "graph adad orientation 1 momentum-space @ {:?} = {}",
        LOOP_MOMENTUM_POINT_NEGATED, orientation1_negated
    );

    let tree = arb_displayed_result(
        precise_probe(&mut cli, "adad_tree", "adad_tree", &[], true, None, None)?,
        true,
    )?;
    println!("tree amplitude momentum-space = {}", tree);

    let tree_gl0 = arb_displayed_result(
        precise_probe(
            &mut cli,
            "adad_tree",
            "adad_tree",
            &[],
            true,
            Some("GL0"),
            None,
        )?,
        true,
    )?;
    println!("tree graph GL0 momentum-space = {}", tree_gl0);

    let tree_gl1 = arb_displayed_result(
        precise_probe(
            &mut cli,
            "adad_tree",
            "adad_tree",
            &[],
            true,
            Some("GL1"),
            None,
        )?,
        true,
    )?;
    println!("tree graph GL1 momentum-space = {}", tree_gl1);

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
