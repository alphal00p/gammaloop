use super::utils::assert_complex_approx_eq;
use super::*;

const INSPECT_POINT: [f64; 3] = [
    9.9999867631375561e-01,
    7.4640036804014731e-01,
    6.1157297272080968e-02,
];

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

#[test]
#[serial]
fn scalar_bubble_inspect() -> Result<()> {
    let mut cli = setup_scalar_bubble_uv_cli("scalar_bubble_inspect")?;

    let bubble_no_integrated = inspect_bubble_process(
        &mut cli,
        "bubble_no_integrated_UV",
        "scalar_bubble_below_thres",
        vec![0, 0],
    )?;
    let bubble_integrated =
        inspect_bubble_process(&mut cli, "bubble", "scalar_bubble_below_thres", vec![0, 1])?;

    assert_complex_approx_eq(
        bubble_no_integrated,
        Complex::new(-1.5806570726409131e-3, 0.0),
        "bubble_no_integrated_UV inspect benchmark",
    );
    assert_complex_approx_eq(
        bubble_integrated,
        Complex::new(1.5856382985462317e-3, 0.0),
        "bubble integrated inspect benchmark",
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
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
