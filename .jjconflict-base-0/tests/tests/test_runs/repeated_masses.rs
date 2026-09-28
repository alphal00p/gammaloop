use super::utils::*;
use super::*;
use spenso::algebra::complex::FloatDerived;

const P_SQRT_5: f64 = 2.23606797749979;
const DELTA_SAMPLES: [f64; 3] = [1.0e-2, 1.0e-3, 1.0e-4];
const PROCESS_A: &str = "repeated_bubble_split_masses";
const PROCESS_B: &str = "repeated_bubble_equal_masses";
const INTEGRAND: &str = "default";

#[derive(Clone, Copy)]
enum MassSource {
    AdditionalParameters,
    ModelParameters,
}

fn setup_repeated_bubble_cli(
    test_name: &str,
    process: &str,
    graph_file: &str,
    mass_source: MassSource,
) -> Result<gammaloop_integration_tests::CLIState> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join(test_name),
        Some(test_name.to_string()),
        true,
    )?;
    let additional_param_values = match mass_source {
        MassSource::AdditionalParameters => "[1.0, 1.01]",
        MassSource::ModelParameters => "[]",
    };
    let mut commands = vec![
        "import model ./assets/models/json/scalars/scalars_2p_3p.json".to_string(),
        "remove processes".to_string(),
        "set global kv global.generation.explicit_orientation_sum_only=true global.generation.evaluator.compile=false global.generation.evaluator.summed=true global.generation.evaluator.summed_function_map=true global.generation.uv.subtract_uv=false global.generation.threshold_subtraction.enable_thresholds=false".to_string(),
        format!(
            r#"set default-runtime string '
[general]
evaluator_method = "Summed"
enable_cache = false
debug_cache = false
integral_unit = "none"
disable_flux_factor = true
mu_r = 3.0
m_uv = 20.0
additional_param_values = {additional_param_values}

[subtraction]
disable_threshold_subtraction = true

[sampling]
graphs = "summed"
orientations = "summed"
lmb_multichanneling = false
lmb_channels = "summed"
coordinate_system = "spherical"
mapping = "linear"
b = 1.0

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [{P_SQRT_5}, 0.0, 0.0, 0.0],
    "dependent"
]
helicities = [0]
'"#
        ),
    ];
    if matches!(mass_source, MassSource::ModelParameters) {
        commands.extend([
            "set model mass_scalar_1=1.0".to_string(),
            "set model mass_scalar_2=1.01".to_string(),
        ]);
    }
    commands.extend([
        format!(
            "import graphs ./tests/resources/graphs/{graph_file} -p {process} -i {INTEGRAND} -o"
        ),
        "generate".to_string(),
    ]);
    let command_refs = commands.iter().map(String::as_str).collect::<Vec<_>>();
    run_commands(&mut cli, &command_refs)?;
    Ok(cli)
}

fn set_mass_delta(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    delta: f64,
    mass_source: MassSource,
) -> Result<()> {
    match mass_source {
        MassSource::AdditionalParameters => cli.run_command(&format!(
            "set process -p {process} -i {INTEGRAND} kv general.additional_param_values=[1.0,{:.17}]",
            1.0 + delta
        )),
        MassSource::ModelParameters => {
            cli.run_command("set model mass_scalar_1=1.0")?;
            cli.run_command(&format!("set model mass_scalar_2={:.17}", 1.0 + delta))
        }
    }
}

#[test]
#[serial]
fn repeated_bubble_split_mass_limit_without_threshold_subtraction() -> Result<()> {
    let point = [0.34, 0.41, 0.23];
    let mut case_a = setup_repeated_bubble_cli(
        "repeated_bubble_split_mass_limit_no_threshold_a",
        PROCESS_A,
        "repeated_bubble_split_masses.dot",
        MassSource::AdditionalParameters,
    )?;
    let mut case_b = setup_repeated_bubble_cli(
        "repeated_bubble_split_mass_limit_no_threshold_b",
        PROCESS_B,
        "repeated_bubble_equal_masses.dot",
        MassSource::AdditionalParameters,
    )?;

    let mut rows = Vec::new();
    for delta in DELTA_SAMPLES {
        set_mass_delta(
            &mut case_a,
            PROCESS_A,
            delta,
            MassSource::AdditionalParameters,
        )?;
        set_mass_delta(
            &mut case_b,
            PROCESS_B,
            delta,
            MassSource::AdditionalParameters,
        )?;
        let value_a = inspect_xspace_process(&mut case_a, PROCESS_A, INTEGRAND, &point)?;
        let value_b = inspect_xspace_process(&mut case_b, PROCESS_B, INTEGRAND, &point)?;
        let abs_diff = complex_distance(value_a, value_b);
        let scale = value_a.norm().max(value_b.norm()).max(1.0);
        rows.push((delta, abs_diff, abs_diff / scale));
    }
    for pair in rows.windows(2) {
        assert!(
            pair[1].1 < pair[0].1,
            "threshold-disabled split-mass inspect: expected abs diff to decrease with Delta, got {rows:#?}"
        );
        assert!(
            pair[1].2 < pair[0].2,
            "threshold-disabled split-mass inspect: expected rel diff to decrease with Delta, got {rows:#?}"
        );
    }
    assert!(
        rows.last().unwrap().2 < 1.0e-3,
        "threshold-disabled split-mass inspect should be close at Delta=1e-4, got {rows:#?}"
    );

    clean_test(&case_a.cli_settings.state.folder);
    clean_test(&case_b.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn repeated_bubble_model_mass_limit_without_threshold_subtraction() -> Result<()> {
    let point = [0.34, 0.41, 0.23];
    let mut case_a = setup_repeated_bubble_cli(
        "repeated_bubble_model_mass_limit_no_threshold_a",
        PROCESS_A,
        "repeated_bubble_split_model_masses.dot",
        MassSource::ModelParameters,
    )?;
    let mut case_b = setup_repeated_bubble_cli(
        "repeated_bubble_model_mass_limit_no_threshold_b",
        PROCESS_B,
        "repeated_bubble_equal_model_masses.dot",
        MassSource::ModelParameters,
    )?;

    let mut rows = Vec::new();
    for delta in DELTA_SAMPLES {
        set_mass_delta(&mut case_a, PROCESS_A, delta, MassSource::ModelParameters)?;
        set_mass_delta(&mut case_b, PROCESS_B, delta, MassSource::ModelParameters)?;
        let value_a = inspect_xspace_process(&mut case_a, PROCESS_A, INTEGRAND, &point)?;
        let value_b = inspect_xspace_process(&mut case_b, PROCESS_B, INTEGRAND, &point)?;
        let abs_diff = complex_distance(value_a, value_b);
        let scale = value_a.norm().max(value_b.norm()).max(1.0);
        rows.push((delta, abs_diff, abs_diff / scale));
    }
    for pair in rows.windows(2) {
        assert!(
            pair[1].1 < pair[0].1,
            "threshold-disabled model-mass inspect: expected abs diff to decrease with Delta, got {rows:#?}"
        );
        assert!(
            pair[1].2 < pair[0].2,
            "threshold-disabled model-mass inspect: expected rel diff to decrease with Delta, got {rows:#?}"
        );
    }
    assert!(
        rows.last().unwrap().2 < 1.0e-3,
        "threshold-disabled model-mass inspect should be close at Delta=1e-4, got {rows:#?}"
    );

    clean_test(&case_a.cli_settings.state.folder);
    clean_test(&case_b.cli_settings.state.folder);
    Ok(())
}
