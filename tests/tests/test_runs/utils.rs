use super::*;
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
            "set global kv global.generation.threshold_subtraction.enable_thresholds=true global.generation.threshold_subtraction.check_esurface_at_generation=true",
            r#"set default-runtime string '
[general]
evaluator_method = "SummedFunctionMap"
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
orientations = "monte_carlo"
lmb_multichanneling = true
lmb_channels = "monte_carlo"

[stability]
rotation_axis = []

[subtraction]
disable_threshold_subtraction = false
'"#,
            "remove processes",
            r#"generate amp g g > h h h / u d c s b QED==3 [{1}]
                --only-diagrams
                --numerator-grouping only_detect_zeroes
                --select-graphs GL15
                --loop-momentum-bases GL15=8
                --global-prefactor-projector '1𝑖 * gammalooprs::ϵ(0,spenso::mink(4,gammalooprs::hedge(0)))
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

pub(super) fn complex_distance(lhs: Complex<f64>, rhs: Complex<f64>) -> f64 {
    let delta_re = lhs.re - rhs.re;
    let delta_im = lhs.im - rhs.im;
    (delta_re * delta_re + delta_im * delta_im).sqrt()
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

pub(super) fn assert_complex_approx_eq(
    actual: Complex<f64>,
    expected: Complex<f64>,
    context: &str,
) {
    let actual_norm = (actual.re * actual.re + actual.im * actual.im).sqrt();
    let expected_norm = (expected.re * expected.re + expected.im * expected.im).sqrt();
    let scale = actual_norm.max(expected_norm).max(1.0);
    let tolerance = 1.0e-10 * scale;
    assert!(
        complex_distance(actual, expected) <= tolerance,
        "{context}: actual={actual}, expected={expected}, tolerance={tolerance}"
    );
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

pub(super) fn complex_ff64(value: &Complex<F<f64>>) -> Complex<f64> {
    Complex::new(value.re.0, value.im.0)
}

pub(super) fn assert_f64_approx_eq(actual: f64, expected: f64, context: &str) {
    let scale = actual.abs().max(expected.abs()).max(1.0);
    let tolerance = 1.0e-10 * scale;
    assert!(
        (actual - expected).abs() <= tolerance,
        "{context}: actual={actual}, expected={expected}, tolerance={tolerance}"
    );
}

pub(super) fn assert_evaluation_outputs_match(
    actual: &gammalooprs::integrands::evaluation::EvaluationResultOutput,
    expected: &gammalooprs::integrands::evaluation::EvaluationResultOutput,
    context: &str,
) {
    match (
        actual.parameterization_jacobian.as_ref(),
        expected.parameterization_jacobian.as_ref(),
    ) {
        (Some(actual), Some(expected)) => {
            assert_f64_approx_eq(actual.0, expected.0, &format!("{context}: jacobian"));
        }
        (None, None) => {}
        _ => panic!("{context}: jacobian presence differs"),
    }
    assert_f64_approx_eq(
        actual.integrator_weight.0,
        expected.integrator_weight.0,
        &format!("{context}: integrator weight"),
    );

    assert_eq!(
        actual.event_groups.len(),
        expected.event_groups.len(),
        "{context}: event-group count differs"
    );
    for (group_index, (actual_group, expected_group)) in actual
        .event_groups
        .iter()
        .zip(expected.event_groups.iter())
        .enumerate()
    {
        assert_eq!(
            actual_group.len(),
            expected_group.len(),
            "{context}: event count differs in group {group_index}"
        );
        let actual_events = actual_group
            .iter()
            .sorted_by_key(|event| {
                (
                    event.cut_info.graph_group_id,
                    event.cut_info.graph_id,
                    event.cut_info.cut_id,
                    event.cut_info.orientation_id,
                    event.cut_info.lmb_channel_id,
                )
            })
            .collect_vec();
        let expected_events = expected_group
            .iter()
            .sorted_by_key(|event| {
                (
                    event.cut_info.graph_group_id,
                    event.cut_info.graph_id,
                    event.cut_info.cut_id,
                    event.cut_info.orientation_id,
                    event.cut_info.lmb_channel_id,
                )
            })
            .collect_vec();

        for (event_index, (actual_event, expected_event)) in
            actual_events.iter().zip(expected_events.iter()).enumerate()
        {
            let event_context = format!(
                "{context}: group {group_index} event {event_index}; actual cut={:?}; expected cut={:?}; actual additional={:?}; expected additional={:?}",
                actual_event.cut_info,
                expected_event.cut_info,
                actual_event.additional_weights.weights,
                expected_event.additional_weights.weights
            );
            assert_eq!(
                actual_event.cut_info.cut_id, expected_event.cut_info.cut_id,
                "{event_context}: cut id differs"
            );
            assert_eq!(
                actual_event.cut_info.graph_id, expected_event.cut_info.graph_id,
                "{event_context}: graph id differs"
            );
            assert_eq!(
                actual_event.cut_info.graph_group_id, expected_event.cut_info.graph_group_id,
                "{event_context}: graph-group id differs"
            );
            assert_eq!(
                actual_event.cut_info.orientation_id, expected_event.cut_info.orientation_id,
                "{event_context}: orientation id differs"
            );
            assert_complex_approx_eq(
                complex_ff64(&actual_event.weight),
                complex_ff64(&expected_event.weight),
                &format!("{event_context}: event weight"),
            );
            assert_eq!(
                actual_event.additional_weights.weights.keys().collect_vec(),
                expected_event
                    .additional_weights
                    .weights
                    .keys()
                    .collect_vec(),
                "{event_context}: additional-weight keys differ"
            );
            for (key, expected_weight) in &expected_event.additional_weights.weights {
                let actual_weight = &actual_event.additional_weights.weights[key];
                assert_complex_approx_eq(
                    complex_ff64(actual_weight),
                    complex_ff64(expected_weight),
                    &format!("{event_context}: additional weight {key:?}"),
                );
            }
        }
    }
    let actual_total = complex_ff64(&actual.integrand_result);
    let expected_total = complex_ff64(&expected.integrand_result);
    let actual_norm =
        (actual_total.re * actual_total.re + actual_total.im * actual_total.im).sqrt();
    let expected_norm =
        (expected_total.re * expected_total.re + expected_total.im * expected_total.im).sqrt();
    let scale = actual_norm.max(expected_norm).max(1.0);
    let tolerance = 1.0e-10 * scale;
    let distance = complex_distance(actual_total, expected_total);
    if distance > tolerance {
        let mut diagnostic = format!(
            "{context}: integrand result: actual={actual_total}, expected={expected_total}, tolerance={tolerance}; actual groups={}, expected groups={}",
            actual.event_groups.len(),
            expected.event_groups.len()
        );
        for (label, groups, total) in [
            ("actual", &actual.event_groups, actual_total),
            ("expected", &expected.event_groups, expected_total),
        ] {
            diagnostic.push_str(&format!("\n{label} integrand_result={total}"));
            for (group_index, group) in groups.iter().enumerate() {
                diagnostic.push_str(&format!(
                    "\n{label} group {group_index} events={}",
                    group.len()
                ));
                for (event_index, event) in group.iter().enumerate() {
                    diagnostic.push_str(&format!(
                        "\n{label} group {group_index} event {event_index} cut={:?} weight={} additional={:?}",
                        event.cut_info,
                        complex_ff64(&event.weight),
                        event.additional_weights.weights
                    ));
                }
            }
        }
        panic!("{diagnostic}");
    }
}

pub(super) fn assert_three_way_local_uv_outputs_match(
    cff_3d: &gammalooprs::integrands::evaluation::SingleSampleEvaluationResult,
    cff_4d: &gammalooprs::integrands::evaluation::SingleSampleEvaluationResult,
    ltd_4d: &gammalooprs::integrands::evaluation::SingleSampleEvaluationResult,
    context: &str,
) {
    assert_evaluation_outputs_match(
        &cff_4d.sample.evaluation,
        &cff_3d.sample.evaluation,
        &format!("{context}: CFF local-4D vs CFF local-3D"),
    );
    assert_evaluation_outputs_match(
        &ltd_4d.sample.evaluation,
        &cff_3d.sample.evaluation,
        &format!("{context}: LTD local-4D vs CFF local-3D"),
    );
}

pub(super) fn evaluation_has_threshold_counterterm(
    evaluation: &gammalooprs::integrands::evaluation::EvaluationResultOutput,
) -> bool {
    evaluation
        .event_groups
        .iter()
        .flat_map(|group| group.iter())
        .any(|event| {
            event.additional_weights.weights.keys().any(|key| {
                matches!(
                    key,
                    gammalooprs::observables::events::AdditionalWeightKey::ThresholdCounterterm {
                        ..
                    }
                )
            })
        })
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
