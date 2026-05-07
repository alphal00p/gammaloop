use gammaloop_api::commands::evaluate_samples::{EvaluateSamplesPrecise, evaluate_sample_precise};
use gammaloop_api::state::ProcessRef;
use gammalooprs::{
    integrands::evaluation::{GenericEvaluationResultOutput, PreciseEvaluationResultOutput},
    observables::{AdditionalWeightKey, GenericEvent},
    utils::{ArbPrec, FloatLike},
};
use ndarray::Array2;
use symbolica::domains::float::Float as SymbolicaFloat;

use super::utils::*;
use super::*;

const P_SQRT_5: f64 = 2.23606797749979;
const M1: f64 = 1.0;
const NO_THRESHOLD_DELTA1_SAMPLES: [f64; 3] = [1.0e-2, 1.0e-3, 1.0e-4];
const UNIT_LOCALIZER_THRESHOLD_DELTA2: f64 = 1.0e-1;
const UNIT_LOCALIZER_DELTA1_SAMPLES: [f64; 3] = [1.0e-3, 1.0e-5, 1.0e-7];
const RUNTIME_REQUIRED_PRECISION: f64 = 1.0e-30;
const RUNTIME_HALF_PRECISION_RELATIVE_TOLERANCE: f64 = 1.0e-15;
const OUTPUT_SIGNIFICANT_DIGITS: usize = 32;
const PROCESS_A: &str = "repeated_bubble_split_masses";
const PROCESS_B: &str = "repeated_bubble_equal_masses";
const INTEGRAND: &str = "default";

#[derive(Debug, Clone, Copy)]
enum MassUpdate {
    AdditionalParameters,
    ModelParameters,
}

#[derive(Debug, Clone, Copy)]
struct ComparisonRow {
    delta: f64,
    case_a: Complex<f64>,
    case_b: Complex<f64>,
    abs_diff: f64,
    rel_diff: f64,
}

#[derive(Debug, Clone)]
struct EventWeightSummary {
    total_weight: Complex<F<ArbPrec>>,
    original: Complex<F<ArbPrec>>,
    threshold_counterterms: Complex<F<ArbPrec>>,
    threshold_counterterm_entries: Vec<(String, Complex<F<ArbPrec>>)>,
    threshold_entry_count: usize,
}

#[derive(Debug, Clone)]
struct EventComparisonRow {
    k_abs: f64,
    delta1: f64,
    delta1_over_delta2: f64,
    delta2: f64,
    case_a: EventWeightSummary,
    case_b: EventWeightSummary,
}

fn threshold_k_abs() -> f64 {
    ((P_SQRT_5 * P_SQRT_5) / 4.0 - M1 * M1).sqrt()
}

fn setup_repeated_bubble_cli(
    test_name: &str,
    process: &str,
    graph_file: &str,
    threshold_subtraction: bool,
    uv_localisation_function: Option<&str>,
) -> Result<gammaloop_integration_tests::CLIState> {
    setup_repeated_bubble_cli_impl(
        test_name,
        process,
        graph_file,
        "./assets/models/json/scalars/scalars_2p_3p.json",
        "[1.0, 1.01]",
        &[],
        threshold_subtraction,
        uv_localisation_function,
    )
}

fn setup_repeated_bubble_model_cli(
    test_name: &str,
    process: &str,
    graph_file: &str,
    threshold_subtraction: bool,
    uv_localisation_function: Option<&str>,
) -> Result<gammaloop_integration_tests::CLIState> {
    let mut cli = setup_repeated_bubble_cli_impl(
        test_name,
        process,
        graph_file,
        "./assets/models/json/scalars/scalars_2p_3p.json",
        "[]",
        &[
            "set model mass_scalar_1=1.0",
            "set model mass_scalar_2=1.01",
        ],
        threshold_subtraction,
        uv_localisation_function,
    )?;
    set_mass_delta(&mut cli, process, 1.0e-2, MassUpdate::ModelParameters)?;
    Ok(cli)
}

fn setup_repeated_bubble_cli_impl(
    test_name: &str,
    process: &str,
    graph_file: &str,
    model_file: &str,
    additional_param_values: &str,
    model_setup_commands: &[&str],
    threshold_subtraction: bool,
    uv_localisation_function: Option<&str>,
) -> Result<gammaloop_integration_tests::CLIState> {
    let mut cli = get_test_cli(
        None,
        get_tests_workspace_path().join(test_name),
        Some(test_name.to_string()),
        true,
    )?;
    let uv_localisation_settings = uv_localisation_function
        .map(|function| {
            format!(
                r#"
[subtraction.local_ct_settings.uv_localisation]
function = "{function}"
"#
            )
        })
        .unwrap_or_default();

    let mut commands = vec![
        format!("import model {model_file}"),
        "remove processes".to_string(),
        format!(
            "set global kv global.3d_representation=CFF global.generation.explicit_orientation_sum_only=true global.generation.evaluator.compile=false global.generation.evaluator.summed=true global.generation.evaluator.summed_function_map=true global.generation.uv.subtract_uv=false global.generation.threshold_subtraction.enable_thresholds={threshold_subtraction}"
        ),
        format!(
            r#"set default-runtime string '
[general]
evaluator_method = "Summed"
enable_cache = false
debug_cache = false
generate_events = true
store_additional_weights_in_event = true
integral_unit = "none"
disable_flux_factor = true
mu_r = 3.0
m_uv = 20.0
additional_param_values = {additional_param_values}

[subtraction]
disable_threshold_subtraction = {}
{uv_localisation_settings}

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
'"#,
            !threshold_subtraction
        ),
    ];
    commands.extend(
        model_setup_commands
            .iter()
            .map(|command| (*command).to_string()),
    );
    commands.push(format!(
        "import graphs ./tests/resources/graphs/{graph_file} -p {process} -i {INTEGRAND} -o"
    ));
    commands.push("generate".to_string());
    let command_refs = commands.iter().map(String::as_str).collect::<Vec<_>>();
    run_commands(&mut cli, &command_refs)?;

    Ok(cli)
}

fn set_mass_delta(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    delta: f64,
    mass_update: MassUpdate,
) -> Result<()> {
    match mass_update {
        MassUpdate::AdditionalParameters => cli.run_command(&format!(
            "set process -p {process} -i {INTEGRAND} kv general.additional_param_values=[1.0,{:.17}]",
            1.0 + delta
        )),
        MassUpdate::ModelParameters => {
            cli.run_command("set model mass_scalar_1=1.0")?;
            cli.run_command(&format!("set model mass_scalar_2={:.17}", 1.0 + delta))
        }
    }
}

fn relative_diff(lhs: Complex<f64>, rhs: Complex<f64>) -> f64 {
    let scale = (lhs.re * lhs.re + lhs.im * lhs.im)
        .sqrt()
        .max((rhs.re * rhs.re + rhs.im * rhs.im).sqrt())
        .max(1.0);
    complex_distance(lhs, rhs) / scale
}

fn abs_arb(value: F<ArbPrec>) -> F<ArbPrec> {
    let zero = value.clone() - value.clone();
    if value < zero { -value } else { value }
}

fn imaginary_ratio_abs_minus_one(
    lhs: &Complex<F<ArbPrec>>,
    rhs: &Complex<F<ArbPrec>>,
) -> F<ArbPrec> {
    let lhs_im_abs = abs_arb(lhs.im.clone());
    let rhs_im_abs = abs_arb(rhs.im.clone());
    let zero = rhs_im_abs.clone() - rhs_im_abs.clone();
    if rhs_im_abs > zero {
        let one = rhs_im_abs.clone() / rhs_im_abs.clone();
        lhs_im_abs / rhs_im_abs - one
    } else {
        zero
    }
}

fn relative_diff_precise(lhs: &Complex<F<ArbPrec>>, rhs: &Complex<F<ArbPrec>>) -> F<ArbPrec> {
    let re_diff = abs_arb(lhs.re.clone() - rhs.re.clone());
    let im_diff = abs_arb(lhs.im.clone() - rhs.im.clone());
    let distance = if re_diff > im_diff { re_diff } else { im_diff };
    let lhs_re = abs_arb(lhs.re.clone());
    let lhs_im = abs_arb(lhs.im.clone());
    let rhs_re = abs_arb(rhs.re.clone());
    let rhs_im = abs_arb(rhs.im.clone());
    let lhs_scale = if lhs_re > lhs_im { lhs_re } else { lhs_im };
    let rhs_scale = if rhs_re > rhs_im { rhs_re } else { rhs_im };
    let scale = if lhs_scale > rhs_scale {
        lhs_scale
    } else {
        rhs_scale
    };
    let zero = scale.clone() - scale.clone();
    if scale > zero.clone() {
        distance / scale
    } else {
        distance
    }
}

fn assert_total_weights_match(
    model: &Complex<F<ArbPrec>>,
    additional_parameters: &Complex<F<ArbPrec>>,
    context: &str,
) {
    let rel_diff = relative_diff_precise(model, additional_parameters);
    let tolerance = F::<ArbPrec>::from_f64(RUNTIME_HALF_PRECISION_RELATIVE_TOLERANCE);
    assert!(
        rel_diff < tolerance,
        "{context}: total weight mismatch between model and additional-parameter setup: model={}, additional={}, rel_diff={}, tolerance={}",
        fmt_complex_32(model),
        fmt_complex_32(additional_parameters),
        fmt_scalar_32(&rel_diff),
        fmt_scalar_32(&tolerance)
    );
}

fn inspect_comparison_row(
    case_a: &mut gammaloop_integration_tests::CLIState,
    case_b: &mut gammaloop_integration_tests::CLIState,
    delta: f64,
    point: &[f64],
    mass_update: MassUpdate,
) -> Result<ComparisonRow> {
    set_mass_delta(case_a, PROCESS_A, delta, mass_update)?;
    set_mass_delta(case_b, PROCESS_B, delta, mass_update)?;
    let value_a = inspect_xspace_process(case_a, PROCESS_A, INTEGRAND, point)?;
    let value_b = inspect_xspace_process(case_b, PROCESS_B, INTEGRAND, point)?;
    Ok(ComparisonRow {
        delta,
        case_a: value_a,
        case_b: value_b,
        abs_diff: complex_distance(value_a, value_b),
        rel_diff: relative_diff(value_a, value_b),
    })
}

fn assert_diffs_decrease(rows: &[ComparisonRow], label: &str) {
    for pair in rows.windows(2) {
        assert!(
            pair[1].abs_diff < pair[0].abs_diff,
            "{label}: expected abs diff to decrease with Delta, got {rows:#?}"
        );
        assert!(
            pair[1].rel_diff < pair[0].rel_diff,
            "{label}: expected rel diff to decrease with Delta, got {rows:#?}"
        );
    }
}

fn zero_complex<T: FloatLike>() -> Complex<F<T>> {
    Complex::new(F(T::new_zero()), F(T::new_zero()))
}

fn one_complex<T: FloatLike>() -> Complex<F<T>> {
    let zero = T::new_zero();
    Complex::new(F(zero.one()), F(T::new_zero()))
}

fn scaled_additional_weight_precise<T: FloatLike>(
    event: &GenericEvent<T>,
    key: &AdditionalWeightKey,
) -> Option<Complex<F<T>>> {
    let weight = event.additional_weights.weights.get(key)?;
    let full_factor = event
        .additional_weights
        .weights
        .get(&AdditionalWeightKey::FullMultiplicativeFactor)
        .cloned()
        .unwrap_or_else(one_complex::<T>);
    Some(weight.clone() * full_factor)
}

fn event_weight_summary(evaluation: GenericEvaluationResultOutput<ArbPrec>) -> EventWeightSummary {
    let mut summary = EventWeightSummary {
        total_weight: evaluation.integrand_result,
        original: zero_complex::<ArbPrec>(),
        threshold_counterterms: zero_complex::<ArbPrec>(),
        threshold_counterterm_entries: Vec::new(),
        threshold_entry_count: 0,
    };

    for event in evaluation
        .event_groups
        .iter()
        .flat_map(|group| group.iter())
    {
        if let Some(original) =
            scaled_additional_weight_precise(event, &AdditionalWeightKey::Original)
        {
            summary.original = summary.original + original;
        }
        for (key, value) in &event.additional_weights.weights {
            if matches!(
                key,
                AdditionalWeightKey::ThresholdCounterterm { .. }
                    | AdditionalWeightKey::AmplitudeThresholdCounterterm { .. }
            ) {
                let full_factor = event
                    .additional_weights
                    .weights
                    .get(&AdditionalWeightKey::FullMultiplicativeFactor)
                    .cloned()
                    .unwrap_or_else(one_complex::<ArbPrec>);
                let scaled_value = value.clone() * full_factor;
                summary.threshold_counterterms =
                    summary.threshold_counterterms + scaled_value.clone();
                summary
                    .threshold_counterterm_entries
                    .push((key.to_string(), scaled_value));
                summary.threshold_entry_count += 1;
            }
        }
    }

    summary
}

fn set_integrand_precision(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    precision: &str,
) -> Result<()> {
    cli.run_command(&format!(
        "set process -p {process} -i {INTEGRAND} string '\n[stability]\ncheck_on_norm = true\nescalate_if_exact_zero = false\nloop_momenta_norm_escalation_factor = -1.0\n\n[[stability.rotation_axis]]\ntype = \"x\"\n\n[[stability.levels]]\nprecision = \"{precision}\"\nrequired_precision_for_re = {RUNTIME_REQUIRED_PRECISION:.1e}\nrequired_precision_for_im = {RUNTIME_REQUIRED_PRECISION:.1e}\nescalate_for_large_weight_threshold = -1.0\n'"
    ))
}

fn evaluate_momentum_process_with_events_precise(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    point: &[f64],
) -> Result<GenericEvaluationResultOutput<ArbPrec>> {
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
    let graph_name = cli
        .state
        .process_list
        .get_integrand(process_id, INTEGRAND)?
        .require_generated()?
        .graph_name_by_id(0)
        .expect("the repeated bubble fixture has one graph")
        .to_string();
    let points = Array2::from_shape_vec((1, point.len()), point.to_vec())?;
    let evaluation = evaluate_sample_precise(
        &mut cli.state,
        &EvaluateSamplesPrecise {
            process_id: Some(process_id),
            integrand_name: Some(INTEGRAND.to_string()),
            use_arb_prec: true,
            minimal_output: false,
            return_generated_events: Some(true),
            momentum_space: true,
            points: points.view(),
            integrator_weights: None,
            discrete_dims: None,
            graph_names: Some(vec![Some(graph_name)]),
            orientations: None,
        },
    )?
    .sample
    .evaluation;

    match evaluation {
        PreciseEvaluationResultOutput::Arb(result) => Ok(result),
        other => Err(eyre::eyre!(
            "expected arbitrary-precision evaluation, got {:?}",
            other.precision()
        )),
    }
}

fn event_summary_for_delta(
    cli: &mut gammaloop_integration_tests::CLIState,
    process: &str,
    delta1: f64,
    point: &[f64],
    mass_update: MassUpdate,
) -> Result<EventWeightSummary> {
    set_mass_delta(cli, process, delta1, mass_update)?;
    let result = evaluate_momentum_process_with_events_precise(cli, process, point)?;
    Ok(event_weight_summary(result))
}

fn event_comparison_row(
    case_a: &mut gammaloop_integration_tests::CLIState,
    case_b: &mut gammaloop_integration_tests::CLIState,
    k_abs: f64,
    delta1: f64,
    delta1_over_delta2: f64,
    delta2: f64,
    point: &[f64],
    mass_update: MassUpdate,
) -> Result<EventComparisonRow> {
    Ok(EventComparisonRow {
        k_abs,
        delta1,
        delta1_over_delta2,
        delta2,
        case_a: event_summary_for_delta(case_a, PROCESS_A, delta1, point, mass_update)?,
        case_b: event_summary_for_delta(case_b, PROCESS_B, delta1, point, mass_update)?,
    })
}

fn fmt_scalar_32(value: &F<ArbPrec>) -> String {
    let value = SymbolicaFloat::from(&value.0);
    format!("{:.*e}", OUTPUT_SIGNIFICANT_DIGITS, value)
}

fn fmt_complex_32(value: &Complex<F<ArbPrec>>) -> String {
    let re = SymbolicaFloat::from(&value.re.0);
    let im = SymbolicaFloat::from(&value.im.0);
    format!(
        "({:+.*e}{:+.*e}i)",
        OUTPUT_SIGNIFICANT_DIGITS, re, OUTPUT_SIGNIFICANT_DIGITS, im
    )
}

#[test]
#[serial]
fn repeated_bubble_split_mass_limit_without_threshold_subtraction() -> Result<()> {
    let point = [0.34, 0.41, 0.23];
    let mut case_a = setup_repeated_bubble_cli(
        "repeated_bubble_split_mass_limit_no_threshold_a",
        PROCESS_A,
        "repeated_bubble_split_masses.dot",
        false,
        None,
    )?;
    let mut case_b = setup_repeated_bubble_cli(
        "repeated_bubble_split_mass_limit_no_threshold_b",
        PROCESS_B,
        "repeated_bubble_equal_masses.dot",
        false,
        None,
    )?;

    let rows = NO_THRESHOLD_DELTA1_SAMPLES
        .iter()
        .map(|delta| {
            inspect_comparison_row(
                &mut case_a,
                &mut case_b,
                *delta,
                &point,
                MassUpdate::AdditionalParameters,
            )
        })
        .collect::<Result<Vec<_>>>()?;
    assert_diffs_decrease(&rows, "threshold-disabled split-mass inspect");
    assert!(
        rows.last().unwrap().rel_diff < 1.0e-3,
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
    let mut case_a = setup_repeated_bubble_model_cli(
        "repeated_bubble_model_mass_limit_no_threshold_a",
        PROCESS_A,
        "repeated_bubble_split_model_masses.dot",
        false,
        None,
    )?;
    let mut case_b = setup_repeated_bubble_model_cli(
        "repeated_bubble_model_mass_limit_no_threshold_b",
        PROCESS_B,
        "repeated_bubble_equal_model_masses.dot",
        false,
        None,
    )?;

    let rows = NO_THRESHOLD_DELTA1_SAMPLES
        .iter()
        .map(|delta| {
            inspect_comparison_row(
                &mut case_a,
                &mut case_b,
                *delta,
                &point,
                MassUpdate::ModelParameters,
            )
        })
        .collect::<Result<Vec<_>>>()?;
    assert_diffs_decrease(&rows, "threshold-disabled model-mass inspect");
    assert!(
        rows.last().unwrap().rel_diff < 1.0e-3,
        "threshold-disabled model-mass inspect should be close at Delta=1e-4, got {rows:#?}"
    );

    clean_test(&case_a.cli_settings.state.folder);
    clean_test(&case_b.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial]
fn repeated_bubble_unit_localizer_ct_sum_converges_and_matches_mass_sources() -> Result<()> {
    let mut model_case_a = setup_repeated_bubble_model_cli(
        "repeated_bubble_unit_localizer_model_threshold_a",
        PROCESS_A,
        "repeated_bubble_split_model_masses.dot",
        true,
        Some("unit"),
    )?;
    let mut model_case_b = setup_repeated_bubble_model_cli(
        "repeated_bubble_unit_localizer_model_threshold_b",
        PROCESS_B,
        "repeated_bubble_equal_model_masses.dot",
        true,
        Some("unit"),
    )?;
    let mut additional_case_a = setup_repeated_bubble_cli(
        "repeated_bubble_unit_localizer_additional_threshold_a",
        PROCESS_A,
        "repeated_bubble_split_masses.dot",
        true,
        Some("unit"),
    )?;
    let mut additional_case_b = setup_repeated_bubble_cli(
        "repeated_bubble_unit_localizer_additional_threshold_b",
        PROCESS_B,
        "repeated_bubble_equal_masses.dot",
        true,
        Some("unit"),
    )?;

    for (cli, process) in [
        (&mut model_case_a, PROCESS_A),
        (&mut model_case_b, PROCESS_B),
        (&mut additional_case_a, PROCESS_A),
        (&mut additional_case_b, PROCESS_B),
    ] {
        set_integrand_precision(cli, process, "Arb")?;
    }

    let k_abs = threshold_k_abs() + UNIT_LOCALIZER_THRESHOLD_DELTA2;
    let point = [k_abs, 0.0, 0.0];
    let mut model_rows = Vec::new();
    let mut additional_rows = Vec::new();

    for delta1 in UNIT_LOCALIZER_DELTA1_SAMPLES {
        model_rows.push(event_comparison_row(
            &mut model_case_a,
            &mut model_case_b,
            k_abs,
            delta1,
            delta1 / UNIT_LOCALIZER_THRESHOLD_DELTA2,
            UNIT_LOCALIZER_THRESHOLD_DELTA2,
            &point,
            MassUpdate::ModelParameters,
        )?);
        additional_rows.push(event_comparison_row(
            &mut additional_case_a,
            &mut additional_case_b,
            k_abs,
            delta1,
            delta1 / UNIT_LOCALIZER_THRESHOLD_DELTA2,
            UNIT_LOCALIZER_THRESHOLD_DELTA2,
            &point,
            MassUpdate::AdditionalParameters,
        )?);
    }

    for row in model_rows.iter().chain(additional_rows.iter()) {
        assert!(
            row.case_a.threshold_entry_count > 0 && row.case_b.threshold_entry_count > 0,
            "Delta={}: expected both cases to expose threshold counterterms: A={:#?}, B={:#?}",
            row.delta1,
            row.case_a,
            row.case_b,
        );
    }

    let x_values = model_rows
        .iter()
        .map(|row| {
            abs_arb(imaginary_ratio_abs_minus_one(
                &row.case_a.threshold_counterterms,
                &row.case_b.threshold_counterterms,
            ))
        })
        .collect::<Vec<_>>();
    let x_1_over_x_2 = x_values[0].clone() / x_values[1].clone();
    let x_2_over_x_3 = x_values[1].clone() / x_values[2].clone();

    assert!(
        x_1_over_x_2 > F::<ArbPrec>::from_f64(50.0),
        "expected CT-sum mismatch to scale down between Delta1={} and Delta1={}, got x1={}, x2={}, x1/x2={}",
        UNIT_LOCALIZER_DELTA1_SAMPLES[0],
        UNIT_LOCALIZER_DELTA1_SAMPLES[1],
        fmt_scalar_32(&x_values[0]),
        fmt_scalar_32(&x_values[1]),
        fmt_scalar_32(&x_1_over_x_2)
    );
    assert!(
        x_2_over_x_3 > F::<ArbPrec>::from_f64(70.0),
        "expected CT-sum mismatch to scale down between Delta1={} and Delta1={}, got x2={}, x3={}, x2/x3={}",
        UNIT_LOCALIZER_DELTA1_SAMPLES[1],
        UNIT_LOCALIZER_DELTA1_SAMPLES[2],
        fmt_scalar_32(&x_values[1]),
        fmt_scalar_32(&x_values[2]),
        fmt_scalar_32(&x_2_over_x_3)
    );

    for (model_row, additional_row) in model_rows.iter().zip(additional_rows.iter()) {
        assert_total_weights_match(
            &model_row.case_a.total_weight,
            &additional_row.case_a.total_weight,
            &format!("case A at Delta1={}", model_row.delta1),
        );
        assert_total_weights_match(
            &model_row.case_b.total_weight,
            &additional_row.case_b.total_weight,
            &format!("case B at Delta1={}", model_row.delta1),
        );
    }

    clean_test(&model_case_a.cli_settings.state.folder);
    clean_test(&model_case_b.cli_settings.state.folder);
    clean_test(&additional_case_a.cli_settings.state.folder);
    clean_test(&additional_case_b.cli_settings.state.folder);
    Ok(())
}
