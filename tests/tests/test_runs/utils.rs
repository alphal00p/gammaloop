use super::*;
pub(super) use gammaloop_integration_tests::evaluation_comparison::{
    CertifiedZeroWeights, EventWeightComparison, assert_complex_approx_eq,
    assert_evaluation_outputs_match, complex_distance, complex_ff64,
};
use gammalooprs::{
    observables::{HistogramSnapshot, HistogramSnapshotKind},
    utils::FloatLike,
};
use symbolica::domains::float::Real;

pub(super) fn decimal_scalar<T>(value: &str) -> Result<F<T>>
where
    T: FloatLike + From<symbolica::domains::float::Float>,
{
    Ok(F(T::from(
        symbolica::domains::float::Float::parse(value, Some(T::new_zero().get_precision()))
            .map_err(|err| eyre::eyre!(err))?,
    )))
}

pub(super) fn decimal_complex<T>(re: &str, im: &str) -> Result<Complex<F<T>>>
where
    T: FloatLike + From<symbolica::domains::float::Float>,
{
    Ok(Complex::new(decimal_scalar(re)?, decimal_scalar(im)?))
}

pub(super) fn assert_complex_approx_eq_precise<T: FloatLike>(
    actual: &Complex<F<T>>,
    expected: &Complex<F<T>>,
    tolerance: &F<T>,
    context: &str,
) {
    let distance = (actual.clone() - expected.clone()).norm().re;
    let actual_norm = actual.norm().re;
    let expected_norm = expected.norm().re;
    let mut scale = if actual_norm > expected_norm {
        actual_norm
    } else {
        expected_norm
    };
    let one = tolerance.clone() / tolerance.clone();
    if scale < one {
        scale = one;
    }
    let scaled_tolerance = tolerance.clone() * scale.clone();
    assert!(
        distance <= scaled_tolerance,
        "{context}: actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}, relative_distance={:e}",
        distance / scale
    );
}

pub(super) fn discrete_histogram_bin_average_and_error(
    histogram: &HistogramSnapshot,
    bin_id: isize,
) -> Result<(f64, f64)> {
    if histogram.kind != HistogramSnapshotKind::Discrete {
        return Err(eyre::eyre!(
            "Expected a discrete histogram for bin lookup, got {:?}",
            histogram.kind
        ));
    }
    let bin = histogram
        .bins
        .iter()
        .find(|bin| bin.bin_id == Some(bin_id))
        .ok_or_else(|| eyre::eyre!("Discrete histogram bin_id {bin_id} is missing"))?;
    Ok((
        bin.average(histogram.sample_count),
        bin.error(histogram.sample_count),
    ))
}

pub(super) fn single_slot_integral(result: &RuntimeIntegrationResult) -> &IntegralEstimate {
    &result
        .single_slot()
        .expect("expected a single integration-result slot")
        .integral
}

pub(super) fn default_integrate_for(name: &str) -> Integrate {
    Integrate {
        process: vec![],
        integrand_name: vec!["default".to_string()],
        n_cores: Some(1),
        workspace_path: Some(
            get_tests_workspace_path().join(format!("{name}/integration_workspace")),
        ),
        target: vec![],
        restart: true,
        ..Default::default()
    }
}

pub(super) fn selected_slot_workspace(
    cli: &gammaloop_integration_tests::CLIState,
    workspace: &Path,
    process: Option<ProcessRef>,
    integrand_name: Option<&str>,
) -> Result<PathBuf> {
    let integrand_name = integrand_name.map(str::to_string);
    let (process_id, integrand_name) = cli
        .state
        .find_integrand_ref(process.as_ref(), integrand_name.as_ref())?;
    let process_name = cli.state.process_list.processes[process_id]
        .definition
        .folder_name
        .clone();
    Ok(workspace
        .join("integrands")
        .join(format!("{process_name}@{integrand_name}")))
}

pub(super) const SCALAR_TRIANGLE_EXTERNALS: &str = r#"set process -p triangle -i scalar_tri string '
[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [1.0, 0.0, 0.0, 0.0],
    [0.5, 0.0, 0.0, -0.5],
    "dependent"
]
helicities = [0, 0, 0]
'"#;

pub(super) const SCALAR_BOX_ABOVE_EXTERNALS: &str = r#"set process -p box -i scalar_box string '
[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [3.0, 0.0, 0.0, 3.0],
    [3.0, 0.0, 0.0, -3.0],
    [3.0, 0.0, 3.0, 0.0],
    "dependent"
]
helicities = [0, 0, 0, 0]
'"#;

pub(super) const SCALAR_BOX_BELOW_EXTERNALS: &str = r#"set process -p box -i scalar_box string '
[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [0.3, 0.0, 0.0, 0.3],
    [0.3, 0.0, 0.0, -0.3],
    [0.3, 0.3, 0.0, 0.0],
    "dependent"
]
helicities = [0, 0, 0, 0]
'"#;

pub(super) const SCALAR_BOX_COPY_ABOVE_EXTERNALS: &str = r#"set process -p box_copy -i scalar_box_copy string '
[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [3.0, 0.0, 0.0, 3.0],
    [3.0, 0.0, 0.0, -3.0],
    [3.0, 0.0, 3.0, 0.0],
    "dependent"
]
helicities = [0, 0, 0, 0]
'"#;

pub(super) const SCALAR_BOX_COPY_BELOW_EXTERNALS: &str = r#"set process -p box_copy -i scalar_box_copy string '
[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
    [0.3, 0.0, 0.0, 0.3],
    [0.3, 0.0, 0.0, -0.3],
    [0.3, 0.3, 0.0, 0.0],
    "dependent"
]
helicities = [0, 0, 0, 0]
'"#;

pub(super) fn setup_scalar_topologies_cli(
    test_name: &str,
) -> Result<gammaloop_integration_tests::CLIState> {
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
            "generate amp scalar_1 > scalar_0 scalar_0 [{1}] --allowed-vertex-interactions V_3_SCALAR_022 V_3_SCALAR_122 -p triangle -i scalar_tri",
            "generate amp scalar_0 scalar_0 > scalar_0 scalar_0 [{1}] --allowed-vertex-interactions V_3_SCALAR_022 -p box -i scalar_box --select-graphs GL0",
            "generate",
            "set model mass_scalar_2=2.0",
            "set model mass_scalar_1=1.0",
            SCALAR_TRIANGLE_EXTERNALS,
            SCALAR_BOX_ABOVE_EXTERNALS,
        ],
    )?;

    Ok(cli)
}

pub(super) fn setup_gg_hhh_threshold_amplitude_cli(
    test_name: &str,
    generation: Option<gammalooprs::settings::global::GenerationSettings>,
) -> Result<gammaloop_integration_tests::CLIState> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join(test_name),
        Some(test_name.to_string()),
        true,
    )?;

    if let Some(generation) = generation {
        cli.cli_settings.global.generation = generation;
    }
    let orientation_sampling = if cli
        .cli_settings
        .global
        .generation
        .explicit_orientation_sum_only
        || cli
            .cli_settings
            .global
            .generation
            .three_dimensional_representations
            .contains(&gammalooprs::settings::global::RepresentationMode::Ltd)
    {
        "summed"
    } else {
        "monte_carlo"
    };
    run_commands(
        &mut cli,
        &[
            "import model sm-default",
            "set global kv global.generation.evaluator.iterative_orientation_optimization=false global.generation.evaluator.store_atom=false global.generation.evaluator.compile=false global.generation.evaluator.summed=false global.generation.evaluator.summed_function_map=true",
            "set global kv global.generation.threshold_subtraction.enable_thresholds=true global.generation.threshold_subtraction.check_esurface_at_generation=true",
            &r#"set default-runtime string '
[general]
evaluator_method = "SingleParametric"
enable_cache = false
debug_cache = false
generate_events = true
store_additional_weights_in_event = true

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
  [500.0, 0.0, 0.0, 500.0],
  [500.0, 0.0, 0.0, -500.0],
  [
    438.5555662246945,
    155.3322001835378,
    348.0160396513587,
    -177.3773615718412
  ],
  [
    356.3696374921922,
    -16.802389008511,
    -318.7291102436005,
    97.48719163688098
  ],
  "dependent"
]
helicities = [1, 1, 0, 0, 0]

[sampling]
graphs = "monte_carlo"
orientations = "ORIENTATION_SAMPLING"
lmb_multichanneling = true
lmb_channels = "monte_carlo"

[stability]
rotation_axis = []

[subtraction]
disable_threshold_subtraction = false
'"#
            .replace("ORIENTATION_SAMPLING", orientation_sampling),
            "remove processes",
            r#"generate amp g g > h h h / u d c s b QED==3 [{1}]
                --only-diagrams
                --numerator-grouping only_detect_zeroes
                --select-graphs GL15
                --loop-momentum-bases GL15=8
                --global-prefactor-projector 'gammalooprs::ϵ(0,spenso::mink(4,gammalooprs::hedge(0)))
                                                * gammalooprs::ϵ(1,spenso::mink(4,gammalooprs::hedge(1)))
                                                * (1/8)*spenso::g(spenso::coad(8,gammalooprs::hedge(0)),spenso::coad(8,gammalooprs::hedge(1)))'
                -p gg_hhh
                -i 1L"#,
            "generate",
            "set model MT=173.0",
            "set model WT=0.0",
            "set model ymt=173.0",
        ],
    )?;

    Ok(cli)
}

pub(super) fn scalar_topology_integrate_command(
    test_name: &str,
    workspace_name: &str,
    selections: &[(&str, &str)],
    targets: &[(&str, Complex<F<f64>>)],
) -> Integrate {
    let workspace = get_tests_workspace_path()
        .join(test_name)
        .join(workspace_name);
    Integrate {
        process: selections
            .iter()
            .map(|(process, _)| ProcessRef::Unqualified((*process).to_string()))
            .collect(),
        integrand_name: selections
            .iter()
            .map(|(_, integrand)| (*integrand).to_string())
            .collect(),
        workspace_path: Some(workspace),
        n_cores: Some(1),
        target: targets
            .iter()
            .map(|(key, target)| format!("{key}={},{}", target.re.0, target.im.0))
            .collect(),
        restart: true,
        ..Default::default()
    }
}

pub(super) fn load_integration_result(path: &Path) -> Result<RuntimeIntegrationResult> {
    Ok(serde_json::from_str(&std::fs::read_to_string(path)?)?)
}

pub(super) fn normalized_integration_result_json(result: &RuntimeIntegrationResult) -> JsonValue {
    let mut value = serde_json::to_value(result).expect("integration result must serialize");
    let Some(slots) = value.get_mut("slots").and_then(JsonValue::as_array_mut) else {
        return value;
    };
    for slot in slots {
        if let Some(statistics) = slot
            .get_mut("integration_statistics")
            .and_then(JsonValue::as_object_mut)
        {
            statistics.remove("average_total_time_seconds");
            statistics.remove("average_parameterization_time_seconds");
            statistics.remove("average_integrand_time_seconds");
            statistics.remove("average_evaluator_time_seconds");
            statistics.remove("average_observable_time_seconds");
            statistics.remove("average_integrator_time_seconds");
        }
    }
    value
}

pub(super) fn assert_integration_results_match_ignoring_timings(
    lhs: &RuntimeIntegrationResult,
    rhs: &RuntimeIntegrationResult,
) {
    assert_eq!(
        normalized_integration_result_json(lhs),
        normalized_integration_result_json(rhs),
    );
}

pub(super) fn default_xspace_point_for(
    cli: &gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
) -> Result<Vec<f64>> {
    let (process_id, integrand_name) = cli.state.find_integrand_ref(
        Some(&ProcessRef::Unqualified(process.to_string())),
        Some(&integrand.to_string()),
    )?;
    let generated = cli
        .state
        .process_list
        .get_integrand(process_id, &integrand_name)?
        .require_generated()?;
    let n_dim = generated.get_n_dim();
    let seed = [0.17, 0.31, 0.53, 0.23, 0.41, 0.67];
    Ok((0..n_dim).map(|index| seed[index % seed.len()]).collect())
}

pub(super) fn nontrivial_xspace_point_for(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
) -> Result<Vec<f64>> {
    let base_point = default_xspace_point_for(cli, process, integrand)?;
    let candidate_seeds = [
        [0.17, 0.31, 0.53, 0.23, 0.41, 0.67],
        [0.73, 0.29, 0.61, 0.47, 0.19, 0.83],
        [0.11, 0.89, 0.37, 0.59, 0.71, 0.43],
        [0.64, 0.52, 0.28, 0.76, 0.34, 0.18],
    ];
    for seed in candidate_seeds {
        let point = (0..base_point.len())
            .map(|index| seed[index % seed.len()])
            .collect_vec();
        let value = inspect_xspace_process(cli, process, integrand, &point)?;
        let magnitude = (value.re * value.re + value.im * value.im).sqrt();
        if magnitude > 1.0e-16 {
            return Ok(point);
        }
    }
    Err(eyre::eyre!(
        "failed to find a non-trivial x-space point for {process}@{integrand} among the test seeds"
    ))
}

pub(super) fn set_process_incoming_helicities(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
    helicities: &str,
) -> Result<()> {
    cli.run_command(&format!(
        "set process -p {process} -i {integrand} string '\n[kinematics.externals.data]\nhelicities = {helicities}\n'"
    ))
}

pub(super) fn set_default_incoming_helicities(
    cli: &mut gammaloop_integration_tests::CLIState,
    helicities: &str,
) -> Result<()> {
    cli.run_command(&format!(
        "set default-runtime string '\n[kinematics.externals.data]\nhelicities = {helicities}\n'"
    ))
}

pub(super) fn inspect_xspace_process(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
    point: &[f64],
) -> Result<Complex<f64>> {
    let (_, inspect) = Inspect {
        process: Some(ProcessRef::Unqualified(process.to_string())),
        integrand_name: Some(integrand.to_string()),
        point: point.to_vec(),
        momentum_space: false,
        ..Default::default()
    }
    .run(cli)?;
    Ok(inspect)
}

pub(super) fn average_over_generated_helicity_processes(
    cli: &mut gammaloop_integration_tests::CLIState,
    processes: &[&str],
    explicit_integrand: &str,
    point: &[f64],
) -> Result<Complex<f64>> {
    let summed = processes.iter().try_fold(
        Complex::new(0.0, 0.0),
        |acc, process| -> Result<Complex<f64>> {
            Result::Ok(acc + inspect_xspace_process(cli, process, explicit_integrand, point)?)
        },
    )?;
    Ok(summed / processes.len() as f64)
}

pub(super) fn mink_dot(left_edge: usize, right_edge: usize, dummy_index: usize) -> String {
    format!(
        "gammalooprs::Q({left_edge},spenso::mink(4,{dummy_index}))*gammalooprs::Q({right_edge},spenso::mink(4,{dummy_index}))",
    )
}

pub(super) fn evaluate_xspace_process_with_events(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    integrand: &str,
    point: &[f64],
    discrete_dims: &[usize],
) -> Result<gammalooprs::integrands::evaluation::SingleSampleEvaluationResult> {
    let (process_id, resolved_integrand_name) = cli.state.find_integrand_ref(
        Some(&ProcessRef::Unqualified(process.to_string())),
        Some(&integrand.to_string()),
    )?;
    let points = ndarray::Array2::from_shape_vec((1, point.len()), point.to_vec())?;
    let discrete_dims =
        ndarray::Array2::from_shape_vec((1, discrete_dims.len()), discrete_dims.to_vec())?;
    evaluate_sample(
        &mut cli.state,
        &EvaluateSamples {
            process_id: Some(process_id),
            integrand_name: Some(resolved_integrand_name),
            use_arb_prec: false,
            minimal_output: false,
            return_generated_events: Some(true),
            momentum_space: false,
            points: points.view(),
            integrator_weights: None,
            discrete_dims: Some(discrete_dims.view()),
            graph_names: None,
            orientations: None,
        },
    )
}

pub(super) fn setup_epem_tth_spin_sum_cli(
    test_name: &str,
) -> Result<gammaloop_integration_tests::CLIState> {
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
            "set global kv global.generation.evaluator.iterative_orientation_optimization=false global.generation.evaluator.store_atom=false global.generation.evaluator.compile=false global.generation.evaluator.summed=false global.generation.evaluator.summed_function_map=true",
            "set global kv global.generation.threshold_subtraction.enable_thresholds=false global.generation.tropical_subgraph_table.disable_tropical_generation=true",
            r#"set default-runtime string '
[general]
evaluator_method = "SummedFunctionMap"
enable_cache = false
debug_cache = false
generate_events = false
store_additional_weights_in_event = false
integral_unit = "picobarn"

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
  [1000.0, 0.0, 0.0, 1000.0],
  [1000.0, 0.0, 0.0, -1000.0],
]
helicities = [1, 1]

[subtraction]
disable_threshold_subtraction = true
'"#,
            "set model MT=173.0",
            "set model WT=0.0",
            "set model ymt=173.0",
        ],
    )?;

    for (process, helicities) in [
        ("epem_a_tth_pp", "[1, 1]"),
        ("epem_a_tth_pm", "[1, -1]"),
        ("epem_a_tth_mp", "[-1, 1]"),
        ("epem_a_tth_mm", "[-1, -1]"),
        (
            "epem_a_tth_sum",
            r#"["summed_averaged", "summed_averaged"]"#,
        ),
    ] {
        set_default_incoming_helicities(&mut cli, helicities)?;
        cli.run_command(&format!(
            "generate xs e+ e- > t t~ h | e+ e- g t t~ h ghG ghG~ a QCD^2==0 QED^2==6 [{{{{2}}}} QCD=0] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --symmetrize-left-right-states true -p {process} -i LO"
        ))?;
        cli.run_command("generate")?;
    }

    Ok(cli)
}

pub(super) fn setup_aa_ttbar_spin_sum_cli(
    test_name: &str,
) -> Result<gammaloop_integration_tests::CLIState> {
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
            "set global kv global.generation.evaluator.iterative_orientation_optimization=false global.generation.evaluator.store_atom=false global.generation.evaluator.compile=false global.generation.evaluator.summed=false global.generation.evaluator.summed_function_map=true",
            "set global kv global.generation.threshold_subtraction.enable_thresholds=false global.generation.tropical_subgraph_table.disable_tropical_generation=true",
            r#"set default-runtime string '
[general]
evaluator_method = "SummedFunctionMap"
enable_cache = false
debug_cache = false
generate_events = false
store_additional_weights_in_event = false
integral_unit = "picobarn"

[kinematics.externals]
type = "constant"

[kinematics.externals.data]
momenta = [
  [500.0, 0.0, 0.0, 500.0],
  [500.0, 0.0, 0.0, -500.0],
]
helicities = [1, 1]

[subtraction]
disable_threshold_subtraction = true
'"#,
            "set model MT=173.0",
            "set model WT=0.0",
        ],
    )?;

    for (process, helicities) in [
        ("aa_ttbar_pp", "[1, 1]"),
        ("aa_ttbar_pm", "[1, -1]"),
        ("aa_ttbar_mp", "[-1, 1]"),
        ("aa_ttbar_mm", "[-1, -1]"),
        ("aa_ttbar_sum", r#"["summed_averaged", "summed_averaged"]"#),
    ] {
        set_default_incoming_helicities(&mut cli, helicities)?;
        cli.run_command(&format!(
            "generate xs a a > t t~ | a t t~ g ghG ghG~ QED^2==4 [{{{{1}}}} QCD=0] --numerator-grouping group_identical_graphs_up_to_scalar_rescaling --symmetrize-left-right-states true -p {process} -i LO"
        ))?;
        cli.run_command("generate")?;
    }

    Ok(cli)
}

#[test]
fn comparisons_reject_nonfinite_and_large_relative_errors() {
    use gammalooprs::integrands::evaluation::EvaluationResult;
    for (actual, expected) in [
        (Complex::new(f64::NAN, 0.0), Complex::new(1.0, 0.0)),
        (Complex::new(f64::INFINITY, 0.0), Complex::new(1.0, 0.0)),
        (Complex::new(1.0e-18, 0.0), Complex::new(2.0e-18, 0.0)),
        (
            Complex::new(1.7e308, 1.7e308),
            Complex::new(0.8e308, 1.7e308),
        ),
    ] {
        assert!(
            std::panic::catch_unwind(|| assert_complex_approx_eq(actual, expected, "invalid pair"))
                .is_err()
        );
        let mut actual_output = EvaluationResult::zero().into_output(true);
        actual_output.integrand_result = Complex::new(F(actual.re), F(actual.im));
        let mut expected_output = EvaluationResult::zero().into_output(true);
        expected_output.integrand_result = Complex::new(F(expected.re), F(expected.im));
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                assert_evaluation_outputs_match(
                    &actual_output,
                    &expected_output,
                    "invalid totals",
                    EventWeightComparison::IndividualOrders,
                );
            }))
            .is_err()
        );
    }
    for (actual, expected) in [(f64::INFINITY, 1.0), (f64::NAN, 1.0), (1.0e-18, 2.0e-18)] {
        for jacobian in [false, true] {
            let mut actual_output = EvaluationResult::zero().into_output(true);
            let mut expected_output = EvaluationResult::zero().into_output(true);
            if jacobian {
                actual_output.parameterization_jacobian = Some(F(actual));
                expected_output.parameterization_jacobian = Some(F(expected));
            } else {
                actual_output.integrator_weight = F(actual);
                expected_output.integrator_weight = F(expected);
            }
            assert!(
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    assert_evaluation_outputs_match(
                        &actual_output,
                        &expected_output,
                        "invalid real weight",
                        EventWeightComparison::IndividualOrders,
                    );
                }))
                .is_err()
            );
        }
    }
    for value in [0.0, 1.0e-300, 1.0e300, 1.7e308] {
        assert_complex_approx_eq(
            Complex::new(value, value),
            Complex::new(value, value),
            "identical finite values",
        );
    }
}

#[test]
fn rich_comparison_rejects_compensating_threshold_component_errors() {
    use gammalooprs::{
        integrands::evaluation::EvaluationResult,
        observables::events::{
            GenericEvent, GenericThresholdCountertermComponentWeight,
            GenericThresholdCountertermEventInfo, ThresholdCountertermComponentOccurrence,
        },
    };

    let mut expected = EvaluationResult::zero().into_output(false);
    let mut event = GenericEvent::default();
    event.additional_weights.threshold_counterterms = Some(GenericThresholdCountertermEventInfo {
        original: Complex::new(F(0.0), F(0.0)),
        components: [1.0, -1.0]
            .into_iter()
            .enumerate()
            .map(
                |(component_id, weight)| GenericThresholdCountertermComponentWeight {
                    component_id,
                    occurrence: ThresholdCountertermComponentOccurrence::LocalUnitarity {
                        overlap_groups: Default::default(),
                        left_threshold_order: Some(1),
                        right_threshold_order: None,
                        lu_cut_order: Some(1),
                    },
                    multiplier_values: vec![F(1.0)].into(),
                    effective_multiplier: F(1.0),
                    bare: Some(Complex::new(F(weight), F(0.0))),
                    weighted: Complex::new(F(weight), F(0.0)),
                    evaluation_skipped: false,
                },
            )
            .collect(),
    });
    expected.event_groups.push_singleton(event);
    let mut actual = expected.clone();
    let decomposition = actual.event_groups[0][0]
        .additional_weights
        .threshold_counterterms
        .as_mut()
        .unwrap();
    decomposition.components[0].weighted.re = F(2.0);
    decomposition.components[1].weighted.re = F(-2.0);
    assert_eq!(decomposition.total(), Complex::new(F(0.0), F(0.0)));
    for comparison in [
        EventWeightComparison::IndividualOrders,
        EventWeightComparison::PhysicalResidues,
    ] {
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                assert_evaluation_outputs_match(
                    &actual,
                    &expected,
                    "compensating simple-pole CT errors",
                    comparison,
                );
            }))
            .is_err()
        );
    }
}

#[test]
fn rich_comparison_sums_only_one_physical_raised_residue() {
    use gammalooprs::{
        integrands::evaluation::EvaluationResult,
        observables::events::{
            GenericEvent, GenericThresholdCountertermComponentWeight,
            GenericThresholdCountertermEventInfo,
            ThresholdCountertermComponentOccurrence as Occurrence,
        },
    };

    let mut expected = EvaluationResult::zero().into_output(false);
    let mut event = GenericEvent::default();
    event
        .kinematic_configuration
        .1
        .push([F(5.0), F(1.0), F(2.0), F(3.0)].into());
    event.additional_weights.threshold_counterterms = Some(GenericThresholdCountertermEventInfo {
        original: Complex::new(F(0.0), F(0.0)),
        components: (0..2)
            .cartesian_product(1..=2)
            .map(|(component_id, order)| {
                let weight = (component_id * 2 + order) as f64;
                GenericThresholdCountertermComponentWeight {
                    component_id,
                    occurrence: Occurrence::LocalUnitarity {
                        overlap_groups: vec![0].into(),
                        left_threshold_order: Some(1),
                        right_threshold_order: None,
                        lu_cut_order: Some(order),
                    },
                    multiplier_values: vec![F(2.0)].into(),
                    effective_multiplier: F(2.0),
                    bare: Some(Complex::new(F(weight), F(0.0))),
                    weighted: Complex::new(F(2.0 * weight), F(0.0)),
                    evaluation_skipped: false,
                }
            })
            .collect(),
    });
    event.weight = event
        .additional_weights
        .threshold_counterterms
        .as_ref()
        .unwrap()
        .total();
    expected.event_groups.push_singleton(event);
    let rejected = |actual: &_, expected: &_, context, comparison| {
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                assert_evaluation_outputs_match(actual, expected, context, comparison);
            }))
            .is_err(),
            "{context} was accepted"
        );
    };

    // A source-zero certificate describes each complete residue, not its
    // off-shell derivative coefficients. Both CFF partials can be nonzero.
    let mut zero_residues = expected.clone();
    let event = &mut zero_residues.event_groups[0][0];
    event.weight = Complex::new_re(F(0.0));
    for (component, weight) in event
        .additional_weights
        .threshold_counterterms
        .as_mut()
        .unwrap()
        .components
        .iter_mut()
        .zip([1.0, -1.0, 2.0, -2.0])
    {
        component.bare = Some(Complex::new_re(F(weight)));
        component.weighted = Complex::new_re(F(2.0 * weight));
    }
    let zeros = CertifiedZeroWeights {
        threshold_components: [(0, 0, 0), (0, 0, 1)].into(),
        dimensionless_tolerance: 1e-10,
        e_cm: 2.0,
        energy_dimension: 3,
        ..Default::default()
    };
    for individual_orders in [true, false] {
        let comparison = EventWeightComparison::CertifiedZeros {
            individual_orders,
            weights: &zeros,
        };
        assert_evaluation_outputs_match(
            &zero_residues,
            &zero_residues,
            "nonzero partials of zero physical residues",
            comparison,
        );
        let mut invalid = zero_residues.clone();
        let components = &mut invalid.event_groups[0][0]
            .additional_weights
            .threshold_counterterms
            .as_mut()
            .unwrap()
            .components;
        components[0].bare.as_mut().unwrap().re.0 += 1.0;
        components[2].bare.as_mut().unwrap().re.0 -= 1.0;
        rejected(
            &invalid,
            &zero_residues,
            "zero certificate cannot combine distinct residues",
            comparison,
        );
    }

    let mut redistributed = expected.clone();
    let components = &mut redistributed.event_groups[0][0]
        .additional_weights
        .threshold_counterterms
        .as_mut()
        .unwrap()
        .components;
    components[0].bare.as_mut().unwrap().re.0 += 3.0;
    components[0].weighted.re.0 += 6.0;
    components[1].bare.as_mut().unwrap().re.0 -= 3.0;
    components[1].weighted.re.0 -= 6.0;
    assert_evaluation_outputs_match(
        &redistributed,
        &expected,
        "same physical residue",
        EventWeightComparison::PhysicalResidues,
    );
    rejected(
        &redistributed,
        &expected,
        "CFF orders remain strict",
        EventWeightComparison::IndividualOrders,
    );

    // Missing partial rows and a different maximum order are permitted when
    // they describe the same complete physical residue and observable.
    let mut fewer_rows = expected.clone();
    let components = &mut fewer_rows.event_groups[0][0]
        .additional_weights
        .threshold_counterterms
        .as_mut()
        .unwrap()
        .components;
    components[0].bare = Some(Complex::new(F(3.0), F(0.0)));
    components[0].weighted = Complex::new(F(6.0), F(0.0));
    if let Occurrence::LocalUnitarity { lu_cut_order, .. } = &mut components[0].occurrence {
        *lu_cut_order = Some(3);
    }
    components.remove(1);
    assert_evaluation_outputs_match(
        &fewer_rows,
        &expected,
        "different partial decomposition",
        EventWeightComparison::PhysicalResidues,
    );

    for distinction in ["component", "overlap", "threshold side"] {
        let mut reference = expected.clone();
        let components = &mut reference.event_groups[0][0]
            .additional_weights
            .threshold_counterterms
            .as_mut()
            .unwrap()
            .components;
        for component in &mut components[2..] {
            if distinction != "component" {
                component.component_id = 0;
            }
            if let Occurrence::LocalUnitarity {
                overlap_groups,
                left_threshold_order,
                right_threshold_order,
                ..
            } = &mut component.occurrence
            {
                if distinction == "overlap" {
                    overlap_groups[0] = 1;
                } else if distinction == "threshold side" {
                    *left_threshold_order = None;
                    *right_threshold_order = Some(1);
                }
            }
        }
        let mut actual = reference.clone();
        let decomposition = actual.event_groups[0][0]
            .additional_weights
            .threshold_counterterms
            .as_mut()
            .unwrap();
        for (index, shift) in [(0, 1.0), (2, -1.0)] {
            decomposition.components[index].bare.as_mut().unwrap().re.0 += shift;
            decomposition.components[index].weighted.re.0 += 2.0 * shift;
        }
        assert_eq!(decomposition.total(), reference.event_groups[0][0].weight);
        rejected(
            &actual,
            &reference,
            distinction,
            EventWeightComparison::PhysicalResidues,
        );
    }

    let mut two_cuts = expected.clone();
    let mut second_cut = two_cuts.event_groups[0][0].clone();
    second_cut.cut_info.cut_id = 1;
    two_cuts.event_groups[0].push(second_cut);
    let mut shifted_cuts = two_cuts.clone();
    for (event, shift) in shifted_cuts.event_groups[0].iter_mut().zip([1.0, -1.0]) {
        let component = &mut event
            .additional_weights
            .threshold_counterterms
            .as_mut()
            .unwrap()
            .components[0];
        component.bare.as_mut().unwrap().re.0 += shift;
        component.weighted.re.0 += 2.0 * shift;
        event.weight.re.0 += 2.0 * shift;
    }
    assert_eq!(
        shifted_cuts.event_groups[0]
            .iter()
            .map(|event| event.weight.re.0)
            .sum::<f64>(),
        two_cuts.event_groups[0]
            .iter()
            .map(|event| event.weight.re.0)
            .sum::<f64>()
    );
    rejected(
        &shifted_cuts,
        &two_cuts,
        "compensating errors in different physical cuts",
        EventWeightComparison::PhysicalResidues,
    );

    for distinction in [
        "cut",
        "graph",
        "channel",
        "observable",
        "multiplier",
        "bare",
        "skipped",
    ] {
        let mut actual = redistributed.clone();
        let event = &mut actual.event_groups[0][0];
        match distinction {
            "cut" => event.cut_info.cut_id += 1,
            "graph" => event.cut_info.graph_id += 1,
            "channel" => event.cut_info.sampling_channel_id = Some(1),
            "observable" => {
                event.kinematic_configuration.1[0] = [F(6.0), F(1.0), F(2.0), F(3.0)].into()
            }
            "multiplier" => {
                event
                    .additional_weights
                    .threshold_counterterms
                    .as_mut()
                    .unwrap()
                    .components[0]
                    .multiplier_values[0] = F(3.0)
            }
            "bare" => {
                event
                    .additional_weights
                    .threshold_counterterms
                    .as_mut()
                    .unwrap()
                    .components[0]
                    .bare
                    .as_mut()
                    .unwrap()
                    .re
                    .0 += 1.0
            }
            "skipped" => {
                event
                    .additional_weights
                    .threshold_counterterms
                    .as_mut()
                    .unwrap()
                    .components[0]
                    .evaluation_skipped = true
            }
            _ => unreachable!(),
        }
        rejected(
            &actual,
            &expected,
            distinction,
            EventWeightComparison::PhysicalResidues,
        );
    }
}

#[test]
fn rich_comparison_bounds_component_underflow_after_normalization() {
    use gammalooprs::{
        integrands::evaluation::EvaluationResult,
        observables::events::{
            AdditionalWeightKey, GenericEvent, GenericThresholdCountertermComponentWeight,
            GenericThresholdCountertermEventInfo, ThresholdCountertermComponentOccurrence,
        },
    };

    let threshold = AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 };
    let factor_key = AdditionalWeightKey::FullMultiplicativeFactor;
    let mut expected = EvaluationResult::zero().into_output(false);
    expected.integrator_weight = F(1.0);
    let mut event = GenericEvent::default();
    for key in [AdditionalWeightKey::Original, threshold] {
        event
            .additional_weights
            .weights
            .insert(key, Complex::new_re(F(0.0)));
    }
    event
        .additional_weights
        .weights
        .insert(factor_key, Complex::new_re(F(1.0)));
    event.additional_weights.threshold_counterterms = Some(GenericThresholdCountertermEventInfo {
        original: Complex::new_re(F(0.0)),
        components: vec![GenericThresholdCountertermComponentWeight {
            component_id: 0,
            occurrence: ThresholdCountertermComponentOccurrence::LocalUnitarity {
                overlap_groups: Default::default(),
                left_threshold_order: Some(1),
                right_threshold_order: None,
                lu_cut_order: Some(1),
            },
            multiplier_values: vec![F(1.0)].into(),
            effective_multiplier: F(1.0),
            bare: Some(Complex::new_re(F(0.0))),
            weighted: Complex::new_re(F(0.0)),
            evaluation_skipped: false,
        }],
    });
    expected.event_groups.push_singleton(event);
    for comparison in [
        EventWeightComparison::IndividualOrders,
        EventWeightComparison::PhysicalResidues,
    ] {
        let rejected = |actual: &_, expected: &_, context| {
            assert!(
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    assert_evaluation_outputs_match(actual, expected, context, comparison);
                }))
                .is_err(),
                "{context} was accepted"
            );
        };
        // Actual GL43/GL47 reporting-boundary tails, including their remaining
        // event factors. The ordinary original weights are unchanged.
        for (tail, factor) in [
            (Complex::new(F(3.034e-321), F(7e-323)), 252.21491664745935),
            (Complex::new(F(1.84797234e-315), F(0.0)), 324.89999920591896),
        ] {
            let mut reference = expected.clone();
            reference.event_groups[0][0]
                .additional_weights
                .weights
                .insert(factor_key, Complex::new_re(F(factor)));
            let mut actual = reference.clone();
            actual.event_groups[0][0]
                .additional_weights
                .weights
                .insert(threshold, tail);
            assert_evaluation_outputs_match(
                &actual,
                &reference,
                "reported threshold tail",
                comparison,
            );
        }

        let tail = f64::MIN_POSITIVE / 4.0;
        for field in ["raw", "original", "bare", "weighted", "event", "total"] {
            for value in [tail, f64::NAN, f64::INFINITY] {
                let mut actual = expected.clone();
                let event = &mut actual.event_groups[0][0];
                let decomposition = event
                    .additional_weights
                    .threshold_counterterms
                    .as_mut()
                    .unwrap();
                let value = Complex::new(F(value), F(-value));
                match field {
                    "raw" => {
                        event.additional_weights.weights.insert(threshold, value);
                    }
                    "original" => decomposition.original = value,
                    "bare" => decomposition.components[0].bare = Some(value),
                    "weighted" => decomposition.components[0].weighted = value,
                    "event" => event.weight = value,
                    "total" => actual.integrand_result = value,
                    _ => unreachable!(),
                }
                if value.re.0.is_finite() {
                    assert_evaluation_outputs_match(&actual, &expected, field, comparison);
                } else {
                    rejected(&actual, &expected, field);
                }
            }
        }
        for factor in [
            Complex::new_re(F(8.0)),
            Complex::new(F(0.0), F(8.0)),
            Complex::new_re(F(f64::NAN)),
            Complex::new_re(F(f64::INFINITY)),
        ] {
            let mut reference = expected.clone();
            reference.event_groups[0][0]
                .additional_weights
                .weights
                .insert(factor_key, factor);
            let mut actual = reference.clone();
            actual.event_groups[0][0]
                .additional_weights
                .weights
                .insert(threshold, Complex::new_re(F(tail)));
            rejected(
                &actual,
                &reference,
                "remaining complex factor rescues a component or is nonfinite",
            );
        }
        for field in [
            "bare",
            "total",
            "unknown factor",
            "separate factor",
            "normal raw",
        ] {
            let mut reference = expected.clone();
            if field == "bare" {
                let component = &mut reference.event_groups[0][0]
                    .additional_weights
                    .threshold_counterterms
                    .as_mut()
                    .unwrap()
                    .components[0];
                component.effective_multiplier = F(8.0);
                component.multiplier_values[0] = F(8.0);
            } else if field == "total" {
                reference.parameterization_jacobian = Some(F(8.0));
            } else if field == "unknown factor" {
                reference.event_groups[0][0]
                    .additional_weights
                    .weights
                    .remove(&factor_key);
            } else if field == "separate factor" || field == "normal raw" {
                reference.event_groups[0][0]
                    .additional_weights
                    .weights
                    .insert(factor_key, Complex::new_re(F(0.0)));
            }
            let mut actual = reference.clone();
            let event = &mut actual.event_groups[0][0];
            match field {
                "bare" => {
                    event
                        .additional_weights
                        .threshold_counterterms
                        .as_mut()
                        .unwrap()
                        .components[0]
                        .bare = Some(Complex::new_re(F(tail)))
                }
                "total" => actual.integrand_result = Complex::new_re(F(tail)),
                "unknown factor" => {
                    event
                        .additional_weights
                        .weights
                        .insert(threshold, Complex::new_re(F(tail)));
                }
                "separate factor" => {
                    event
                        .additional_weights
                        .weights
                        .insert(factor_key, Complex::new_re(F(tail)));
                }
                "normal raw" => {
                    event
                        .additional_weights
                        .weights
                        .insert(threshold, Complex::new_re(F(1.0)));
                }
                _ => unreachable!(),
            }
            rejected(&actual, &reference, field);
        }
        for sign in [-1.0, 1.0] {
            let mut reference = expected.clone();
            reference.event_groups[0][0]
                .additional_weights
                .weights
                .insert(factor_key, Complex::new(F(1.0), F(sign)));
            let mut actual = reference.clone();
            actual.event_groups[0][0]
                .additional_weights
                .weights
                .insert(threshold, Complex::new(F(3.0 * tail), F(3.0 * tail)));
            rejected(
                &actual,
                &reference,
                "two subnormal raw phases reinforce one normal physical phase",
            );
            reference.event_groups[0][0]
                .additional_weights
                .weights
                .insert(threshold, Complex::new_re(F(3.0 * tail)));
            actual.event_groups[0][0]
                .additional_weights
                .weights
                .insert(threshold, Complex::new(F(3.0 * tail), F(1.2 * tail)));
            rejected(
                &actual,
                &reference,
                "erasing only one raw phase changes a normal physical phase",
            );
            reference.event_groups[0][0]
                .additional_weights
                .weights
                .insert(threshold, Complex::new_re(F(0.0)));
            actual.event_groups[0][0]
                .additional_weights
                .weights
                .insert(threshold, Complex::new(F(tail), F(tail)));
            assert_evaluation_outputs_match(
                &actual,
                &reference,
                "both mixed physical phases remain subnormal",
                comparison,
            );
        }
        let mut large_factor = expected.clone();
        large_factor.event_groups[0][0]
            .additional_weights
            .weights
            .insert(factor_key, Complex::new(F(f64::MAX), F(f64::MAX)));
        assert_evaluation_outputs_match(
            &large_factor,
            &large_factor,
            "finite factor with overflowing underflow bound",
            comparison,
        );
        for swap in [false, true] {
            let mut reference = expected.clone();
            let mut actual = expected.clone();
            for (output, normal, tiny) in [(&mut reference, 1.0, 0.0), (&mut actual, 2.0, tail)] {
                let weight = if swap {
                    Complex::new(F(tiny), F(normal))
                } else {
                    Complex::new(F(normal), F(tiny))
                };
                output.event_groups[0][0]
                    .additional_weights
                    .weights
                    .insert(threshold, weight);
            }
            rejected(
                &actual,
                &reference,
                "the other phase must retain its normal-scale comparison",
            );
        }
    }
}

#[test]
fn certified_zero_weights_follow_energy_dimensions_and_event_normalization() {
    use gammalooprs::{
        integrands::evaluation::EvaluationResult,
        observables::events::{
            AdditionalWeightKey, GenericEvent, GenericThresholdCountertermComponentWeight,
            GenericThresholdCountertermEventInfo, ThresholdCountertermComponentOccurrence,
        },
    };

    for dimension in [3, 5] {
        for energy in [1.0_f64, 2.0, 7.0] {
            let physical_scale = energy.powi(dimension);
            // Three spatial loop measures and one decay flux have dimension 8.
            // The raw additional weights carry the complementary dimension.
            let normalization = 1000.0 * energy.powi(8);
            let zeros = CertifiedZeroWeights {
                original_cuts: [(0, 0)].into(),
                threshold_cuts: [(0, 0)].into(),
                threshold_components: [(0, 0, 0)].into(),
                all_cuts_zero: true,
                dimensionless_tolerance: 1e-10,
                e_cm: energy,
                energy_dimension: dimension,
            };
            let mut expected = EvaluationResult::zero().into_output(false);
            let mut event = GenericEvent::default();
            for key in [
                AdditionalWeightKey::Original,
                AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 },
            ] {
                event
                    .additional_weights
                    .weights
                    .insert(key, Complex::new_re(F(0.0)));
            }
            event.additional_weights.weights.insert(
                AdditionalWeightKey::FullMultiplicativeFactor,
                Complex::new_re(F(normalization)),
            );
            event.additional_weights.threshold_counterterms =
                Some(GenericThresholdCountertermEventInfo {
                    original: Complex::new_re(F(0.0)),
                    components: vec![GenericThresholdCountertermComponentWeight {
                        component_id: 0,
                        occurrence: ThresholdCountertermComponentOccurrence::LocalUnitarity {
                            overlap_groups: Default::default(),
                            left_threshold_order: Some(1),
                            right_threshold_order: None,
                            lu_cut_order: Some(1),
                        },
                        multiplier_values: vec![F(1.0)].into(),
                        effective_multiplier: F(1.0),
                        bare: Some(Complex::new_re(F(0.0))),
                        weighted: Complex::new_re(F(0.0)),
                        evaluation_skipped: false,
                    }],
                });
            expected.event_groups.push_singleton(event);
            let mut actual = expected.clone();
            let event = &mut actual.event_groups[0][0];
            event.weight.re = F(5e-11 * physical_scale);
            event
                .additional_weights
                .weights
                .get_mut(&AdditionalWeightKey::Original)
                .unwrap()
                .re = F(2e-11 * physical_scale / normalization);
            event
                .additional_weights
                .weights
                .get_mut(&AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 })
                .unwrap()
                .re = F(3e-11 * physical_scale / normalization);
            let decomposition = event
                .additional_weights
                .threshold_counterterms
                .as_mut()
                .unwrap();
            decomposition.original.re = F(2e-11 * physical_scale);
            decomposition.components[0].bare.as_mut().unwrap().re = F(3e-11 * physical_scale);
            decomposition.components[0].weighted.re = F(3e-11 * physical_scale);
            actual.integrand_result.re = event.weight.re;
            for individual_orders in [true, false] {
                let comparison = EventWeightComparison::CertifiedZeros {
                    individual_orders,
                    weights: &zeros,
                };
                assert_evaluation_outputs_match(&actual, &expected, "dimensioned zero", comparison);
                for field in [
                    "raw",
                    "bare",
                    "weighted",
                    "original",
                    "event",
                    "total",
                    "multiplier",
                    "skip",
                ] {
                    let mut invalid = actual.clone();
                    let event = &mut invalid.event_groups[0][0];
                    let decomposition = event
                        .additional_weights
                        .threshold_counterterms
                        .as_mut()
                        .unwrap();
                    match field {
                        "raw" => {
                            event
                                .additional_weights
                                .weights
                                .get_mut(&AdditionalWeightKey::Original)
                                .unwrap()
                                .re = F(2e-10 * physical_scale / normalization)
                        }
                        "bare" => {
                            decomposition.components[0].bare.as_mut().unwrap().re =
                                F(2e-10 * physical_scale)
                        }
                        "weighted" => {
                            decomposition.components[0].weighted.re = F(2e-10 * physical_scale)
                        }
                        "original" => decomposition.original.re = F(2e-10 * physical_scale),
                        "event" => event.weight.re = F(2e-10 * physical_scale),
                        "total" => invalid.integrand_result.re = F(2e-10 * physical_scale),
                        "multiplier" => decomposition.components[0].effective_multiplier = F(2.0),
                        "skip" => decomposition.components[0].evaluation_skipped = true,
                        _ => unreachable!(),
                    }
                    assert!(
                        std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            assert_evaluation_outputs_match(&invalid, &expected, field, comparison);
                        }))
                        .is_err(),
                        "zero bound accepted invalid {field} at E_cm={energy}, dimension={dimension}"
                    );
                }
            }
        }
    }
}

#[test]
fn certified_zero_weights_do_not_cover_other_physical_residues() {
    use gammalooprs::{
        integrands::evaluation::EvaluationResult,
        observables::events::{AdditionalWeightKey, GenericEvent},
    };
    let zeros = CertifiedZeroWeights {
        original_cuts: [(0, 0)].into(),
        // This physical cut has a nonzero threshold: only its original is zero.
        dimensionless_tolerance: 1e-10,
        e_cm: 2.0,
        energy_dimension: 3,
        ..Default::default()
    };
    let mut expected = EvaluationResult::zero().into_output(false);
    for cut_id in [0, 1] {
        let mut event = GenericEvent::default();
        event.cut_info.cut_id = cut_id;
        event
            .additional_weights
            .weights
            .insert(AdditionalWeightKey::Original, Complex::new_re(F(0.0)));
        event.additional_weights.weights.insert(
            AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 },
            Complex::new_re(F(1e-12)),
        );
        event.additional_weights.weights.insert(
            AdditionalWeightKey::FullMultiplicativeFactor,
            Complex::new_re(F(1.0)),
        );
        expected.event_groups.push_singleton(event);
    }
    let comparison = EventWeightComparison::CertifiedZeros {
        individual_orders: false,
        weights: &zeros,
    };
    for field in ["other_original", "threshold", "event", "total", "cut_id"] {
        let mut invalid = expected.clone();
        match field {
            "other_original" => {
                invalid.event_groups[1][0]
                    .additional_weights
                    .weights
                    .get_mut(&AdditionalWeightKey::Original)
                    .unwrap()
                    .re = F(1e-12)
            }
            "threshold" => {
                invalid.event_groups[0][0]
                    .additional_weights
                    .weights
                    .get_mut(&AdditionalWeightKey::ThresholdCounterterm { subset_index: 0 })
                    .unwrap()
                    .re = F(2e-12)
            }
            "event" => invalid.event_groups[0][0].weight.re = F(1e-12),
            "total" => invalid.integrand_result.re = F(1e-12),
            "cut_id" => invalid.event_groups[1][0].cut_info.cut_id = 2,
            _ => unreachable!(),
        }
        assert!(
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                assert_evaluation_outputs_match(&invalid, &expected, field, comparison);
            }))
            .is_err(),
            "zero bound leaked to uncertified {field}"
        );
    }
}
