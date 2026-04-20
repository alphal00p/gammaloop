use super::utils::{
    assert_complex_approx_eq, assert_complex_approx_eq_precise, decimal_complex, decimal_scalar,
};
use super::*;
use gammaloop_api::commands::evaluate_samples::{EvaluateSamplesPrecise, evaluate_sample_precise};
use gammalooprs::{
    integrands::evaluation::PreciseEvaluationResultOutput,
    utils::{ArbPrec, FloatLike, PrecisionUpgradable, f128},
};
use ndarray::Array2;
use symbolica::domains::float::Float as SymbolicaFloat;

const INSPECT_POINT: [f64; 3] = [
    9.999_986_763_137_556e-1,
    7.464_003_680_401_473e-1,
    6.115_729_727_208_097e-2,
];

fn bubble_no_integrated_inspect_f64_target() -> Complex<f64> {
    Complex::new(-1.580_657_072_640_913e-3, 0.0)
}

fn bubble_integrated_inspect_f64_target() -> Complex<f64> {
    Complex::new(1.585_638_298_408_533e-3, 0.0)
}

fn bubble_no_integrated_inspect_arb_target() -> Result<Complex<F<ArbPrec>>> {
    decimal_complex(
        "-1.58065707259153789178502163142117068526463969013510333564276531480561962745854241380044043667161360949437009469502210362095446423045005192306907814379082388080206628068736848837078345292618746320848779749755287907430238646092292488091093980397229936469812174411521638680368115422804978273334133533394619e-3",
        "0",
    )
}

fn bubble_integrated_inspect_arb_target() -> Result<Complex<F<ArbPrec>>> {
    decimal_complex(
        "1.58563829854623146817685961891979140882549385669771327706729197278504857216366149312468787585506886338192073840071485443029528929091987002980931067724774119021136369785653453585590970672064551013866039650085231428300325815259378884538340543405264269548949898960571949215728536882750006781538331536613484e-3",
        "0",
    )
}

fn bubble_no_integrated_inspect_quad_target() -> Result<Complex<F<f128>>> {
    decimal_complex("-1.58065707259153789178502133686503e-3", "0")
}

fn bubble_integrated_inspect_quad_target() -> Result<Complex<F<f128>>> {
    decimal_complex("1.58563829854623146817685992280347e-3", "0")
}

fn bubble_integrated_target() -> Complex<F<f64>> {
    Complex::new(F(2.7029875684552542e-3), F(0.0))
}

fn bubble_no_integrated_target() -> Complex<F<f64>> {
    Complex::new(F(1.471664021721597e-2), F(0.0))
}

fn bubble_runtime_block(process: &str, integrand: &str) -> String {
    format!(
        r#"set process -p {process} -i {integrand} string '
[general]
mu_r = 3.0
m_uv = 20.0

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [1.0, 0.0, 0.0, 0.0],
    "dependent"
]
helicities = [0, 0]

[sampling]
graphs = "monte_carlo"
orientations = "monte_carlo"
lmb_multichanneling = false
lmb_channels = "summed"
coordinate_system = "spherical"
mapping = "linear"
b = 1.0
'"#,
    )
}

fn setup_scalar_bubble_uv_cli(test_name: &str) -> Result<gammaloop_integration_tests::CLIState> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join(test_name),
        Some(test_name.to_string()),
        true,
    )?;

    run_commands(
        &mut cli,
        &[
            "import model scalars-default.json",
            "remove processes",
            "set global kv global.generation.evaluator.iterative_orientation_optimization=false global.generation.evaluator.store_atom=true global.generation.evaluator.compile=false global.generation.evaluator.summed=false global.generation.evaluator.summed_function_map=true",
            "set global kv global.generation.threshold_subtraction.enable_thresholds=true global.generation.threshold_subtraction.check_esurface_at_generation=false",
            "set global kv global.generation.uv.generate_integrated=true",
            "generate amp scalar_1 > scalar_1 [{1}] --allowed-vertex-interactions V_3_SCALAR_122 -p bubble -i scalar_bubble_below_thres",
            "set global kv global.generation.uv.generate_integrated=false",
            "generate amp scalar_1 > scalar_1 [{1}] --allowed-vertex-interactions V_3_SCALAR_122 -p bubble_no_integrated_UV -i scalar_bubble_below_thres",
            "set model mass_scalar_2=2.0",
            "set process kv integrator.target_relative_accuracy=0.00001 integrator.n_increase=0 integrator.n_start=1000000 integrator.n_max=1000000 integrator.seed=1337",
        ],
    )?;

    cli.run_command(&bubble_runtime_block("bubble", "scalar_bubble_below_thres"))?;
    cli.run_command(&bubble_runtime_block(
        "bubble_no_integrated_UV",
        "scalar_bubble_below_thres",
    ))?;

    Ok(cli)
}

fn inspect_bubble_process(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
    discrete_dim: Vec<usize>,
) -> Result<Complex<f64>> {
    let (_, inspect) = Inspect {
        process: Some(ProcessRef::Unqualified(process.to_string())),
        integrand_name: Some(integrand.to_string()),
        point: INSPECT_POINT.to_vec(),
        momentum_space: false,
        discrete_dim,
        ..Default::default()
    }
    .run(cli)?;
    Ok(inspect)
}

fn set_single_precision_level(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
    precision: &str,
) -> Result<()> {
    let thresholds = if precision == "Double" {
        "required_precision_for_re = 1.0e-10\nrequired_precision_for_im = 1.0e-10\nescalate_for_large_weight_threshold = 0.9"
    } else {
        "required_precision_for_re = 1.0e-30\nrequired_precision_for_im = 1.0e-30\nescalate_for_large_weight_threshold = -1.0"
    };
    cli.run_command(&format!(
        "set process -p {process} -i {integrand} string '\n[stability]\ncheck_on_norm = true\nescalate_if_exact_zero = false\nloop_momenta_norm_escalation_factor = -1.0\n\n[[stability.rotation_axis]]\ntype = \"x\"\n\n[[stability.levels]]\nprecision = \"{precision}\"\n{thresholds}\n'"
    ))
}

fn displayed_xspace_precise_result<T: FloatLike>(
    result: gammalooprs::integrands::evaluation::GenericEvaluationResultOutput<T>,
) -> Result<Complex<F<T>>> {
    let jacobian = result.parameterization_jacobian.ok_or_else(|| {
        eyre::eyre!("x-space precise inspect requires a parameterization jacobian")
    })?;
    Ok(result
        .integrand_result
        .map(|entry| entry * jacobian.clone()))
}

fn precise_inspect_bubble_process(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
    discrete_dim: Vec<usize>,
    use_arb_prec: bool,
) -> Result<PreciseEvaluationResultOutput> {
    let (process_id, integrand_name) = cli.state.find_integrand_ref(
        Some(&ProcessRef::Unqualified(process.to_string())),
        Some(&integrand.to_string()),
    )?;
    let points = Array2::from_shape_vec((1, INSPECT_POINT.len()), INSPECT_POINT.to_vec())?;
    let discrete_dims = Array2::from_shape_vec((1, discrete_dim.len()), discrete_dim)?;
    Ok(evaluate_sample_precise(
        &mut cli.state,
        &EvaluateSamplesPrecise {
            process_id: Some(process_id),
            integrand_name: Some(integrand_name),
            use_arb_prec,
            minimal_output: false,
            return_generated_events: Some(true),
            momentum_space: false,
            points: points.view(),
            integrator_weights: None,
            discrete_dims: Some(discrete_dims.view()),
            graph_names: None,
            orientations: None,
        },
    )?
    .sample
    .evaluation)
}

fn quad_result(result: PreciseEvaluationResultOutput) -> Result<Complex<F<f128>>> {
    match result {
        PreciseEvaluationResultOutput::Quad(result) => displayed_xspace_precise_result(result),
        _ => Err(eyre::eyre!("Expected a Quad precise evaluation result")),
    }
}

fn arb_result(result: PreciseEvaluationResultOutput) -> Result<Complex<F<ArbPrec>>> {
    match result {
        PreciseEvaluationResultOutput::Arb(result) => displayed_xspace_precise_result(result),
        _ => Err(eyre::eyre!("Expected an Arb precise evaluation result")),
    }
}

fn exact_f64_as_quad(value: Complex<f64>) -> Complex<F<f128>> {
    Complex::new(
        F(f128::from(SymbolicaFloat::from(value.re))),
        F(f128::from(SymbolicaFloat::from(value.im))),
    )
}

fn exact_quad_as_arb(value: &Complex<F<f128>>) -> Complex<F<ArbPrec>> {
    Complex::new(F(value.re.0.higher()), F(value.im.0.higher()))
}

fn quad_vs_f64_tolerance() -> F<f128> {
    F(f128::from(SymbolicaFloat::from(f64::EPSILON.sqrt())))
}

fn quad_benchmark_tolerance() -> Result<F<f128>> {
    decimal_scalar("1.0e-27")
}

fn arb_vs_quad_tolerance() -> Result<F<ArbPrec>> {
    decimal_scalar("2.0e-16")
}

mod important {
    use super::*;

    #[test]
    #[serial]
    fn scalar_bubble_inspect() -> Result<()> {
        let mut cli = setup_scalar_bubble_uv_cli("scalar_bubble_inspect")?;

        set_single_precision_level(
            &mut cli,
            "bubble_no_integrated_UV",
            "scalar_bubble_below_thres",
            "Double",
        )?;
        let bubble_no_integrated_f64 = inspect_bubble_process(
            &mut cli,
            "bubble_no_integrated_UV",
            "scalar_bubble_below_thres",
            vec![0, 0],
        )?;
        assert_complex_approx_eq(
            bubble_no_integrated_f64,
            bubble_no_integrated_inspect_f64_target(),
            "bubble_no_integrated_UV inspect f64 benchmark",
        );

        set_single_precision_level(&mut cli, "bubble", "scalar_bubble_below_thres", "Double")?;
        let bubble_integrated_f64 =
            inspect_bubble_process(&mut cli, "bubble", "scalar_bubble_below_thres", vec![0, 1])?;
        assert_complex_approx_eq(
            bubble_integrated_f64,
            bubble_integrated_inspect_f64_target(),
            "bubble integrated inspect f64 benchmark",
        );

        set_single_precision_level(
            &mut cli,
            "bubble_no_integrated_UV",
            "scalar_bubble_below_thres",
            "Quad",
        )?;
        let bubble_no_integrated_quad = quad_result(precise_inspect_bubble_process(
            &mut cli,
            "bubble_no_integrated_UV",
            "scalar_bubble_below_thres",
            vec![0, 0],
            false,
        )?)?;
        assert_complex_approx_eq_precise(
            &bubble_no_integrated_quad,
            &bubble_no_integrated_inspect_quad_target()?,
            &quad_benchmark_tolerance()?,
            "bubble_no_integrated_UV inspect Quad benchmark",
        );
        assert_complex_approx_eq_precise(
            &bubble_no_integrated_quad,
            &exact_f64_as_quad(bubble_no_integrated_f64),
            &quad_vs_f64_tolerance(),
            "bubble_no_integrated_UV inspect Quad compatible with f64",
        );

        set_single_precision_level(&mut cli, "bubble", "scalar_bubble_below_thres", "Quad")?;
        let bubble_integrated_quad = quad_result(precise_inspect_bubble_process(
            &mut cli,
            "bubble",
            "scalar_bubble_below_thres",
            vec![0, 1],
            false,
        )?)?;
        assert_complex_approx_eq_precise(
            &bubble_integrated_quad,
            &bubble_integrated_inspect_quad_target()?,
            &quad_benchmark_tolerance()?,
            "bubble integrated inspect Quad benchmark",
        );
        assert_complex_approx_eq_precise(
            &bubble_integrated_quad,
            &exact_f64_as_quad(bubble_integrated_f64),
            &quad_vs_f64_tolerance(),
            "bubble integrated inspect Quad compatible with f64",
        );

        set_single_precision_level(
            &mut cli,
            "bubble_no_integrated_UV",
            "scalar_bubble_below_thres",
            "Arb",
        )?;
        let bubble_no_integrated_arb = arb_result(precise_inspect_bubble_process(
            &mut cli,
            "bubble_no_integrated_UV",
            "scalar_bubble_below_thres",
            vec![0, 0],
            true,
        )?)?;
        assert_complex_approx_eq_precise(
            &bubble_no_integrated_arb,
            &bubble_no_integrated_inspect_arb_target()?,
            &decimal_scalar("1.0e-60")?,
            "bubble_no_integrated_UV inspect ArbPrec benchmark",
        );
        assert_complex_approx_eq_precise(
            &bubble_no_integrated_arb,
            &exact_quad_as_arb(&bubble_no_integrated_quad),
            &arb_vs_quad_tolerance()?,
            "bubble_no_integrated_UV inspect ArbPrec compatible with Quad",
        );

        set_single_precision_level(&mut cli, "bubble", "scalar_bubble_below_thres", "Arb")?;
        let bubble_integrated_arb = arb_result(precise_inspect_bubble_process(
            &mut cli,
            "bubble",
            "scalar_bubble_below_thres",
            vec![0, 1],
            true,
        )?)?;
        assert_complex_approx_eq_precise(
            &bubble_integrated_arb,
            &bubble_integrated_inspect_arb_target()?,
            &decimal_scalar("1.0e-60")?,
            "bubble integrated inspect ArbPrec benchmark",
        );
        assert_complex_approx_eq_precise(
            &bubble_integrated_arb,
            &exact_quad_as_arb(&bubble_integrated_quad),
            &arb_vs_quad_tolerance()?,
            "bubble integrated inspect ArbPrec compatible with Quad",
        );

        clean_test(&cli.cli_settings.state.folder);
        Ok(())
    }
}

#[test]
#[serial]
fn scalar_bubble_integrated() -> Result<()> {
    let mut cli = setup_scalar_bubble_uv_cli("scalar_bubble_integrated")?;

    let integration_result = Integrate {
        process: vec![
            ProcessRef::Unqualified("bubble_no_integrated_UV".to_string()),
            ProcessRef::Unqualified("bubble".to_string()),
        ],
        integrand_name: vec![
            "scalar_bubble_below_thres".to_string(),
            "scalar_bubble_below_thres".to_string(),
        ],
        workspace_path: Some(
            get_tests_workspace_path()
                .join("scalar_bubble_integrated")
                .join("integration_workspace"),
        ),
        n_cores: Some(1),
        target: vec![
            format!(
                "bubble_no_integrated_UV@scalar_bubble_below_thres={},{}",
                bubble_no_integrated_target().re.0,
                bubble_no_integrated_target().im.0
            ),
            format!(
                "bubble@scalar_bubble_below_thres={},{}",
                bubble_integrated_target().re.0,
                bubble_integrated_target().im.0
            ),
        ],
        restart: true,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let no_integrated_slot = integration_result
        .slot("bubble_no_integrated_UV@scalar_bubble_below_thres")
        .expect("bubble_no_integrated_UV slot should exist");
    let integrated_slot = integration_result
        .slot("bubble@scalar_bubble_below_thres")
        .expect("bubble slot should exist");

    assert!(
        no_integrated_slot
            .integral
            .is_compatible_with_target(bubble_no_integrated_target(), 2),
        "Integration result {} not compatible with target {}",
        no_integrated_slot.integral,
        bubble_no_integrated_target()
    );
    assert!(
        integrated_slot
            .integral
            .is_compatible_with_target(bubble_integrated_target(), 2),
        "Integration result {} not compatible with target {}",
        integrated_slot.integral,
        bubble_integrated_target()
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}
