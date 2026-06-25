use std::{collections::BTreeMap, fs, path::PathBuf};

use color_eyre::Result;
use gammaloop_api::commands::{
    Profile,
    evaluate_samples::{EvaluateSamples, evaluate_sample},
    integrate::Integrate,
    profile::{InfraRedProfile, UltraVioletProfile},
};
use gammaloop_api::state::ProcessRef;
use gammaloop_integration_tests::{
    CLIState, clean_test, get_example_cli, get_test_cli, get_tests_workspace_path, workspace_root,
};
use gammalooprs::integrands::{HasIntegrand, evaluation::EvaluationMetaData};
use gammalooprs::observables::events::AdditionalWeightKey;
use gammalooprs::settings::runtime::{IntegralEstimate, SlotIntegrationResult};
use gammalooprs::utils::F;
use ndarray::Array2;
use serde_json::{Map, Value, json};
use spenso::algebra::complex::Complex;
use tabled::{Table, Tabled};

const INSPECT_DEPENDENCE_ACCURACY_FACTOR: f64 = 1000.0;
const INSPECT_DEPENDENCE_SCALE_EXPONENTS: [f64; 5] =
    [0.0, 1.0, 2.0, 3.0, UV_PROFILE_MIN_SCALE_EXPONENT];
const INSPECT_DEPENDENCE_SEEDS: [[f64; 6]; 2] = [
    [0.11, -0.07, 0.19, -0.13, 0.05, 0.29],
    [0.23, 0.31, -0.17, 0.41, -0.29, 0.13],
];
const UV_PROFILE_MIN_SCALE_EXPONENT: f64 = 4.0;
const UV_PROFILE_MAX_SCALE_EXPONENT: f64 = 8.0;
const UV_PROFILE_N_POINTS: usize = 33;
const INTEGRATED_CT_RELATIVE_ERROR_LIMIT: f64 = 0.10;

#[derive(Clone, Copy)]
struct IntegratedUvIntegratorSettings {
    target_relative_accuracy: f64,
    n_start: usize,
    n_increase: usize,
    n_max: usize,
    n_cores: usize,
}

const DEFAULT_INTEGRATED_UV_INTEGRATOR: IntegratedUvIntegratorSettings =
    IntegratedUvIntegratorSettings {
        target_relative_accuracy: 0.001,
        n_start: 50_000,
        n_increase: 0,
        n_max: 400_000,
        n_cores: 1,
    };

const SUNRISE_INTEGRATED_UV_INTEGRATOR: IntegratedUvIntegratorSettings =
    IntegratedUvIntegratorSettings {
        target_relative_accuracy: INTEGRATED_CT_RELATIVE_ERROR_LIMIT,
        n_start: 50_000,
        n_increase: 0,
        n_max: 1_000_000,
        n_cores: 1,
    };

#[derive(Default)]
struct IntegratedUvTargets {
    integrated: Option<Complex<F<f64>>>,
}

impl IntegratedUvTargets {
    fn any(&self) -> bool {
        self.integrated.is_some()
    }
}

struct IntegratedUvCase<'a> {
    run_card: &'a str,
    test_name: &'a str,
    process: &'a str,
    integrand_name: &'a str,
    integrator: IntegratedUvIntegratorSettings,
    original_m_uv: f64,
    shifted_m_uv: f64,
    original_renormalization_localization_scale: f64,
    shifted_renormalization_localization_scale: f64,
    original_mu_r: f64,
    shifted_mu_r: f64,
    skip_uv_profile: bool,
    targets: IntegratedUvTargets,
    integrated_ct_relative_error_limit: Option<f64>,
    check_mu_r_dependence: bool,
}

struct IntegratedUvResults {
    integrated: IntegralEstimate,
    integrated_muv_shifted: Option<IntegralEstimate>,
    integrated_per_graph: Option<BTreeMap<String, IntegralEstimate>>,
    integrated_muv_shifted_per_graph: Option<BTreeMap<String, IntegralEstimate>>,
    integrated_rls_shifted_per_graph: Option<BTreeMap<String, IntegralEstimate>>,
}

struct CheckRow {
    check: &'static str,
    graph: String,
    passed: bool,
    sigma: Option<f64>,
    threshold: Option<String>,
}

struct IntegratedUvCaseResult {
    graph: String,
    uv_profile_passed: Option<bool>,
    muv_invariance_passed: bool,
    muv_rows: Vec<CheckRow>,
    rls_invariance_passed: bool,
    rls_rows: Vec<CheckRow>,
    mur_dependence_passed: Option<bool>,
    mur_rows: Vec<CheckRow>,
    target_passed: Option<bool>,
    target_threshold: Option<&'static str>,
    integrated_ct_accuracy_passed: Option<bool>,
    integrated_ct_accuracy_rows: Vec<CheckRow>,
    error: Option<String>,
}

impl IntegratedUvCaseResult {
    fn all_passed(&self) -> bool {
        self.error.is_none()
            && self.uv_profile_passed.unwrap_or(true)
            && self.muv_invariance_passed
            && self.rls_invariance_passed
            && self.mur_dependence_passed.unwrap_or(true)
            && self.target_passed.unwrap_or(true)
            && self.integrated_ct_accuracy_passed.unwrap_or(true)
    }
}

fn set_fast_deterministic_integrator(
    cli: &mut CLIState,
    settings: IntegratedUvIntegratorSettings,
) -> Result<()> {
    cli.run_command(&format!(
        "set process kv integrator.target_relative_accuracy={} \
         integrator.n_increase={} integrator.n_start={} \
         integrator.n_max={} integrator.seed=1337",
        settings.target_relative_accuracy, settings.n_increase, settings.n_start, settings.n_max
    ))
}

fn shifted_muv_integrand_name(integrand_name: &str) -> String {
    format!("{integrand_name}_muv_shifted")
}

fn shifted_mur_integrand_name(integrand_name: &str) -> String {
    format!("{integrand_name}_mur_shifted")
}

fn shifted_rls_integrand_name(integrand_name: &str) -> String {
    format!("{integrand_name}_rls_shifted")
}

fn set_original_integrated_uv_scales(
    cli: &mut CLIState,
    case: &IntegratedUvCase<'_>,
) -> Result<()> {
    cli.run_command(&format!(
        "set process -p {} -i {} kv general.m_uv={}",
        case.process, case.integrand_name, case.original_m_uv
    ))?;
    cli.run_command(&format!(
        "set process -p {} -i {} kv general.mu_r={}",
        case.process, case.integrand_name, case.original_mu_r
    ))?;
    cli.run_command(&format!(
        "set process -p {} -i {} kv general.renormalization_localization_scale={}",
        case.process, case.integrand_name, case.original_renormalization_localization_scale
    ))?;

    Ok(())
}

fn add_integrated_uv_scale_variants(cli: &mut CLIState, case: &IntegratedUvCase<'_>) -> Result<()> {
    if (case.original_m_uv - case.shifted_m_uv).abs() > f64::EPSILON {
        let shifted_name = shifted_muv_integrand_name(case.integrand_name);
        cli.run_command(&format!(
            "duplicate integrand -p {} -i {} --output_process_name {} --output_integrand_name {}",
            case.process, case.integrand_name, case.process, shifted_name
        ))?;
        cli.run_command(&format!(
            "set process -p {} -i {} kv general.m_uv={}",
            case.process, shifted_name, case.shifted_m_uv
        ))?;
    }

    let shifted_rls_name = shifted_rls_integrand_name(case.integrand_name);
    cli.run_command(&format!(
        "duplicate integrand -p {} -i {} --output_process_name {} --output_integrand_name {}",
        case.process, case.integrand_name, case.process, shifted_rls_name
    ))?;
    cli.run_command(&format!(
        "set process -p {} -i {} kv general.renormalization_localization_scale={}",
        case.process, shifted_rls_name, case.shifted_renormalization_localization_scale
    ))?;

    if case.check_mu_r_dependence && (case.original_mu_r - case.shifted_mu_r).abs() > f64::EPSILON {
        let shifted_name = shifted_mur_integrand_name(case.integrand_name);
        cli.run_command(&format!(
            "duplicate integrand -p {} -i {} --output_process_name {} --output_integrand_name {}",
            case.process, case.integrand_name, case.process, shifted_name
        ))?;
        cli.run_command(&format!(
            "set process -p {} -i {} kv general.mu_r={}",
            case.process, shifted_name, case.shifted_mu_r
        ))?;
    }

    Ok(())
}
fn graph_breakdown_estimates(
    slot: &SlotIntegrationResult,
) -> Option<BTreeMap<String, IntegralEstimate>> {
    let re_breakdown = slot
        .grid_breakdown
        .re
        .as_ref()
        .filter(|breakdown| breakdown.axis_label == "graph");
    let im_breakdown = slot
        .grid_breakdown
        .im
        .as_ref()
        .filter(|breakdown| breakdown.axis_label == "graph");

    if re_breakdown.is_none() && im_breakdown.is_none() {
        return None;
    }

    let mut estimates = BTreeMap::new();
    if let Some(breakdown) = re_breakdown {
        for entry in &breakdown.entries {
            let label = entry
                .bin_label
                .clone()
                .unwrap_or_else(|| entry.bin_index.to_string());
            estimates.insert(
                label,
                IntegralEstimate {
                    neval: entry.processed_samples,
                    real_zero: 0,
                    im_zero: 0,
                    result: Complex::new(entry.value, F(0.0)),
                    error: Complex::new(entry.error, F(0.0)),
                    real_chisq: entry.chi_sq,
                    im_chisq: F(0.0),
                },
            );
        }
    }

    if let Some(breakdown) = im_breakdown {
        for entry in &breakdown.entries {
            let label = entry
                .bin_label
                .clone()
                .unwrap_or_else(|| entry.bin_index.to_string());
            let estimate = estimates.entry(label).or_insert_with(|| IntegralEstimate {
                neval: entry.processed_samples,
                real_zero: 0,
                im_zero: 0,
                result: Complex::new(F(0.0), F(0.0)),
                error: Complex::new(F(0.0), F(0.0)),
                real_chisq: F(0.0),
                im_chisq: F(0.0),
            });
            estimate.neval = estimate.neval.max(entry.processed_samples);
            estimate.result.im = entry.value;
            estimate.error.im = entry.error;
            estimate.im_chisq = entry.chi_sq;
        }
    }

    Some(estimates)
}

fn graph_estimates_for_slot(
    cli: &CLIState,
    slot: &SlotIntegrationResult,
) -> Result<Option<BTreeMap<String, IntegralEstimate>>> {
    if let Some(estimates) = graph_breakdown_estimates(slot) {
        return Ok(Some(estimates));
    }

    let info = cli.state.get_integrand_info(
        Some(&ProcessRef::Unqualified(slot.process.clone())),
        Some(&slot.integrand.clone()),
    )?;
    let graph_names: Vec<_> = info
        .graph_groups
        .iter()
        .flat_map(|group| group.graphs.iter().map(|graph| graph.name.clone()))
        .collect();

    Ok((graph_names.len() == 1)
        .then(|| BTreeMap::from([(graph_names[0].clone(), slot.integral.clone())])))
}

fn integral_relative_error(estimate: &IntegralEstimate) -> f64 {
    let error_norm = estimate.error.re.0.hypot(estimate.error.im.0);
    let result_norm = estimate
        .result
        .re
        .0
        .hypot(estimate.result.im.0)
        .max(f64::MIN_POSITIVE);
    error_norm / result_norm
}

fn integrated_ct_relative_error_rows(results: &IntegratedUvResults, limit: f64) -> Vec<CheckRow> {
    let mut rows = Vec::new();
    for (check, graph, estimate) in [
        ("integrated.relerr", "integrated", Some(&results.integrated)),
        (
            "m_uv.integrated.relerr",
            "m_uv-shifted",
            results.integrated_muv_shifted.as_ref(),
        ),
    ] {
        if let Some(estimate) = estimate {
            let relative_error = integral_relative_error(estimate);
            rows.push(CheckRow {
                check,
                graph: graph.to_string(),
                passed: relative_error.is_finite() && relative_error <= limit,
                sigma: None,
                threshold: Some(format!(
                    "relative error {:.2}% <= {:.2}%",
                    100.0 * relative_error,
                    100.0 * limit
                )),
            });
        }
    }
    rows
}

fn run_integrated_uv_integration(
    cli: &mut CLIState,
    case: &IntegratedUvCase<'_>,
) -> Result<IntegratedUvResults> {
    let mut processes = vec![ProcessRef::Unqualified(case.process.to_string())];
    let mut integrand_names = vec![case.integrand_name.to_string()];

    let muv_shifted_name = ((case.original_m_uv - case.shifted_m_uv).abs() > f64::EPSILON)
        .then(|| shifted_muv_integrand_name(case.integrand_name));
    if let Some(name) = &muv_shifted_name {
        processes.push(ProcessRef::Unqualified(case.process.to_string()));
        integrand_names.push(name.clone());
    }

    let rls_shifted_name = shifted_rls_integrand_name(case.integrand_name);
    processes.push(ProcessRef::Unqualified(case.process.to_string()));
    integrand_names.push(rls_shifted_name.clone());

    set_fast_deterministic_integrator(cli, case.integrator)?;
    let integration_result = Integrate {
        process: processes,
        integrand_name: integrand_names,
        workspace_path: Some(get_tests_workspace_path().join(format!(
            "{}/integration_workspace_{}",
            case.test_name, case.integrand_name
        ))),
        n_cores: Some(case.integrator.n_cores),
        restart: true,
        renderer: gammaloop_api::commands::integrate::RendererOption::Tabled,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?;

    let integrated_key = format!("{}@{}", case.process, case.integrand_name);
    let integrated_slot = integration_result
        .slot(&integrated_key)
        .expect("integrated slot should exist");
    let integrated_per_graph = graph_estimates_for_slot(cli, integrated_slot)?;
    let integrated_muv_shifted = muv_shifted_name.as_ref().map(|name| {
        integration_result
            .slot(&format!("{}@{}", case.process, name))
            .expect("integrated shifted-m_uv slot should exist")
            .integral
            .clone()
    });
    let integrated_muv_shifted_per_graph = muv_shifted_name
        .as_ref()
        .map(|name| {
            graph_estimates_for_slot(
                cli,
                integration_result
                    .slot(&format!("{}@{}", case.process, name))
                    .expect("integrated shifted-m_uv slot should exist"),
            )
        })
        .transpose()?
        .flatten();
    let integrated_rls_shifted_per_graph = graph_estimates_for_slot(
        cli,
        integration_result
            .slot(&format!("{}@{}", case.process, rls_shifted_name))
            .expect("integrated shifted-rls slot should exist"),
    )?;

    Ok(IntegratedUvResults {
        integrated: integrated_slot.integral.clone(),
        integrated_muv_shifted,
        integrated_per_graph,
        integrated_muv_shifted_per_graph,
        integrated_rls_shifted_per_graph,
    })
}

fn required_stability_accuracy_floor(
    cli: &mut CLIState,
    process: &str,
    integrand_name: &str,
) -> Result<f64> {
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
    let integrand = cli
        .state
        .process_list
        .get_integrand(process_id, integrand_name)?
        .require_generated()?;

    Ok(integrand
        .get_settings()
        .stability
        .levels
        .first()
        .map(|level| {
            level
                .required_precision_for_re
                .max(level.required_precision_for_im)
        })
        .unwrap_or(f64::EPSILON))
}

fn stability_relative_accuracy(metadata: &EvaluationMetaData, accuracy_floor: f64) -> f64 {
    let stability_accuracy = metadata
        .stability_results
        .iter()
        .filter_map(|result| result.estimated_relative_accuracy.as_ref())
        .fold(0.0_f64, |max, accuracy| max.max(accuracy.0.abs()));
    let instability_error = metadata
        .relative_instability_error
        .re
        .0
        .abs()
        .max(metadata.relative_instability_error.im.0.abs());

    stability_accuracy
        .max(instability_error)
        .max(accuracy_floor.abs())
        .max(f64::EPSILON)
}

struct InspectProbePoint {
    point: Vec<f64>,
    scale_exponent: f64,
    seed_index: usize,
}

struct InspectProbeResult {
    probe: InspectProbePoint,
    relative_delta: f64,
    relative_accuracy: f64,
    required_relative_delta: f64,
}

impl InspectProbeResult {
    fn threshold_ratio(&self) -> f64 {
        let ratio = self.relative_delta / self.required_relative_delta.max(f64::MIN_POSITIVE);
        if ratio.is_finite() {
            ratio
        } else {
            f64::NEG_INFINITY
        }
    }
}

fn deterministic_uv_momentum_points(
    cli: &mut CLIState,
    process: &str,
    integrand_name: &str,
) -> Result<Vec<InspectProbePoint>> {
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
    let integrand = cli
        .state
        .process_list
        .get_integrand(process_id, integrand_name)?
        .require_generated()?;
    let e_cm = integrand.get_settings().kinematics.e_cm;
    let n_dim = integrand.get_n_dim();

    Ok(INSPECT_DEPENDENCE_SCALE_EXPONENTS
        .iter()
        .copied()
        .flat_map(|scale_exponent| {
            INSPECT_DEPENDENCE_SEEDS
                .iter()
                .enumerate()
                .map(move |(seed_index, seed)| {
                    let scale = 10_f64.powf(scale_exponent) * e_cm;
                    let point = (0..n_dim)
                        .map(|index| scale * seed[index % seed.len()])
                        .collect();
                    InspectProbePoint {
                        point,
                        scale_exponent,
                        seed_index,
                    }
                })
        })
        .collect())
}

struct InspectEvaluation {
    value: Complex<f64>,
    relative_accuracy: f64,
}

fn evaluate_summed_momentum_sample(
    cli: &mut CLIState,
    process: &str,
    integrand_name: &str,
    point: &[f64],
    accuracy_floor: f64,
) -> Result<InspectEvaluation> {
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
    let graph_name = cli
        .state
        .process_list
        .get_integrand(process_id, integrand_name)?
        .require_generated()?
        .graph_name_by_id(0)
        .ok_or_else(|| eyre::eyre!("Generated {integrand_name} integrand should have graph id 0"))?
        .to_string();
    let points = Array2::from_shape_vec((1, point.len()), point.to_vec())?;
    let result = evaluate_sample(
        &mut cli.state,
        &EvaluateSamples {
            process_id: Some(process_id),
            integrand_name: Some(integrand_name.to_string()),
            use_arb_prec: false,
            minimal_output: false,
            return_generated_events: Some(false),
            momentum_space: true,
            points: points.view(),
            integrator_weights: None,
            discrete_dims: None,
            graph_names: Some(vec![Some(graph_name)]),
            orientations: Some(vec![None]),
        },
    )?;
    let evaluation = result.sample.evaluation;
    let metadata = evaluation
        .evaluation_metadata
        .as_ref()
        .ok_or_else(|| eyre::eyre!("inspect evaluation should include metadata"))?;

    Ok(InspectEvaluation {
        value: evaluation.integrand_result.map(|entry| entry.0),
        relative_accuracy: stability_relative_accuracy(metadata, accuracy_floor),
    })
}

fn inspect_scale_dependence_row(
    cli: &mut CLIState,
    case: &IntegratedUvCase<'_>,
    baseline_name: &str,
    check: &'static str,
    shifted_name: &str,
    allow_absent_dependence: bool,
) -> Result<Vec<CheckRow>> {
    let accuracy_floor = required_stability_accuracy_floor(cli, case.process, baseline_name)?;
    let probes = deterministic_uv_momentum_points(cli, case.process, baseline_name)?;
    let n_probes = probes.len();
    let mut best: Option<InspectProbeResult> = None;
    for probe in probes {
        let baseline = evaluate_summed_momentum_sample(
            cli,
            case.process,
            baseline_name,
            &probe.point,
            accuracy_floor,
        )?;
        let shifted = evaluate_summed_momentum_sample(
            cli,
            case.process,
            shifted_name,
            &probe.point,
            accuracy_floor,
        )?;
        let delta_norm =
            (shifted.value.re - baseline.value.re).hypot(shifted.value.im - baseline.value.im);
        let scale = baseline
            .value
            .re
            .hypot(baseline.value.im)
            .max(shifted.value.re.hypot(shifted.value.im))
            .max(f64::MIN_POSITIVE);
        let relative_delta = delta_norm / scale;
        let relative_accuracy = baseline
            .relative_accuracy
            .max(shifted.relative_accuracy)
            .max(f64::EPSILON);
        let required_relative_delta = INSPECT_DEPENDENCE_ACCURACY_FACTOR * relative_accuracy;
        let probe_result = InspectProbeResult {
            probe,
            relative_delta,
            relative_accuracy,
            required_relative_delta,
        };
        if match &best {
            Some(current) => probe_result.threshold_ratio() > current.threshold_ratio(),
            None => true,
        } {
            best = Some(probe_result);
        }
    }

    let best = best.ok_or_else(|| eyre::eyre!("{check} should have at least one inspect probe"))?;
    let dependence_found =
        best.relative_delta.is_finite() && best.relative_delta >= best.required_relative_delta;
    let absent_dependence_allowed = allow_absent_dependence && best.relative_delta == 0.0;
    let passed = dependence_found || absent_dependence_allowed;
    let threshold = if absent_dependence_allowed {
        format!(
            "no inspect-level dependence across {n_probes} probes (allowed; integrated invariance checked separately)"
        )
    } else {
        format!(
            "best inspect relative delta {:.3e} >= {INSPECT_DEPENDENCE_ACCURACY_FACTOR:.0}*accuracy {:.3e} (accuracy={:.3e}, exponent={:+.1}, seed={}, probes={n_probes})",
            best.relative_delta,
            best.required_relative_delta,
            best.relative_accuracy,
            best.probe.scale_exponent,
            best.probe.seed_index,
        )
    };

    Ok(vec![CheckRow {
        check,
        graph: "inspect".to_string(),
        passed,
        sigma: None,
        threshold: Some(threshold),
    }])
}

fn integrated_uv_targets_pass(
    results: &IntegratedUvResults,
    targets: &IntegratedUvTargets,
) -> Result<bool> {
    if targets
        .integrated
        .as_ref()
        .is_some_and(|integrated_target| {
            !results
                .integrated
                .is_compatible_with_target(*integrated_target, 2)
        })
    {
        return Ok(false);
    }

    Ok(true)
}

fn integral_estimate_change_sigma(lhs: &IntegralEstimate, rhs: &IntegralEstimate) -> f64 {
    let delta_re = lhs.result.re.0 - rhs.result.re.0;
    let delta_im = lhs.result.im.0 - rhs.result.im.0;
    let delta_norm = delta_re.hypot(delta_im);
    let combined_error_norm = lhs
        .error
        .re
        .0
        .hypot(rhs.error.re.0)
        .hypot(lhs.error.im.0.hypot(rhs.error.im.0));

    if combined_error_norm == 0.0 {
        if delta_norm == 0.0 {
            0.0
        } else {
            f64::INFINITY
        }
    } else {
        delta_norm / combined_error_norm
    }
}

fn integrated_uv_profile_passes(
    cli: &mut CLIState,
    process: &str,
    integrand_name: &str,
) -> Result<bool> {
    let res = Profile::UltraViolet(UltraVioletProfile {
        process: Some(ProcessRef::Unqualified(process.to_string())),
        integrand_name: Some(integrand_name.to_string()),
        min_scale_exponent: UV_PROFILE_MIN_SCALE_EXPONENT,
        max_scale_exponent: UV_PROFILE_MAX_SCALE_EXPONENT,
        n_points: UV_PROFILE_N_POINTS,
        per_orientation: true,
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?;
    let uv = res.unwrap_uv();
    Ok(uv.pass_fail(-0.9).failed == 0)
}

fn integrated_result_scale_invariance_rows(
    results: &IntegratedUvResults,
    check: &'static str,
    shifted_per_graph: &Option<BTreeMap<String, IntegralEstimate>>,
) -> Result<Vec<CheckRow>> {
    match (&results.integrated_per_graph, shifted_per_graph) {
        (Some(baseline_graphs), Some(shifted_graphs)) if !baseline_graphs.is_empty() => {
            let mut rows = Vec::new();
            for (graph, baseline_graph) in baseline_graphs {
                let shifted_graph = shifted_graphs.get(graph).ok_or_else(|| {
                    eyre::eyre!(
                        "shifted {check} graph breakdown is missing graph '{}'",
                        graph
                    )
                })?;
                let sigma = integral_estimate_change_sigma(shifted_graph, baseline_graph);
                rows.push(CheckRow {
                    check,
                    graph: graph.clone(),
                    passed: sigma.is_finite() && sigma <= 2.0,
                    sigma: Some(sigma),
                    threshold: Some("compatible within 2sigma".to_string()),
                });
            }
            for graph in shifted_graphs.keys() {
                if !baseline_graphs.contains_key(graph) {
                    return Err(eyre::eyre!(
                        "baseline graph breakdown is missing shifted {check} graph '{}'",
                        graph
                    ));
                }
            }
            Ok(rows)
        }
        (None, None) => Err(eyre::eyre!(
            "{check} comparison requires graph breakdowns for baseline and shifted results"
        )),
        (Some(_), None) => Err(eyre::eyre!(
            "shifted {check} result has no graph breakdown to compare against baseline"
        )),
        (None, Some(_)) => Err(eyre::eyre!(
            "baseline result has no graph breakdown to compare against shifted {check}"
        )),
        (Some(_), Some(_)) => Err(eyre::eyre!(
            "{check} comparison requires a non-empty baseline graph breakdown"
        )),
    }
}

fn run_integrated_uv_case(case: &IntegratedUvCase<'_>) -> IntegratedUvCaseResult {
    let mut outcome = IntegratedUvCaseResult {
        graph: case.test_name.to_string(),
        uv_profile_passed: (!case.skip_uv_profile).then_some(false),
        muv_invariance_passed: false,
        muv_rows: Vec::new(),
        rls_invariance_passed: false,
        rls_rows: Vec::new(),
        mur_dependence_passed: case.check_mu_r_dependence.then_some(false),
        mur_rows: Vec::new(),
        target_passed: case.targets.any().then_some(false),
        target_threshold: case.targets.any().then_some("2sigma"),
        integrated_ct_accuracy_passed: case.integrated_ct_relative_error_limit.map(|_| false),
        integrated_ct_accuracy_rows: Vec::new(),
        error: None,
    };

    let cli_result = get_test_cli(
        Some(format!("{}.toml", case.run_card).into()),
        get_tests_workspace_path().join(case.test_name),
        None,
        true,
    );
    let mut cli = match cli_result {
        Ok(cli) => cli,
        Err(err) => {
            outcome.error = Some(format!("setup failed: {err:#}"));
            return outcome;
        }
    };

    let case_result = (|| -> Result<()> {
        cli.run_command("run generate")?;
        set_original_integrated_uv_scales(&mut cli, case)?;
        add_integrated_uv_scale_variants(&mut cli, case)?;

        if !case.skip_uv_profile {
            outcome.uv_profile_passed = Some(integrated_uv_profile_passes(
                &mut cli,
                case.process,
                case.integrand_name,
            )?);
        }

        let baseline = run_integrated_uv_integration(&mut cli, case)?;
        if let Some(limit) = case.integrated_ct_relative_error_limit {
            outcome.integrated_ct_accuracy_rows =
                integrated_ct_relative_error_rows(&baseline, limit);
            outcome.integrated_ct_accuracy_passed = Some(
                outcome
                    .integrated_ct_accuracy_rows
                    .iter()
                    .all(|row| row.passed),
            );
        }

        if (case.original_m_uv - case.shifted_m_uv).abs() > f64::EPSILON {
            let shifted_muv_name = shifted_muv_integrand_name(case.integrand_name);
            outcome.muv_rows = inspect_scale_dependence_row(
                &mut cli,
                case,
                case.integrand_name,
                "m_uv.inspect",
                &shifted_muv_name,
                false,
            )?;
        }
        outcome
            .muv_rows
            .extend(integrated_result_scale_invariance_rows(
                &baseline,
                "m_uv.integrated",
                &baseline.integrated_muv_shifted_per_graph,
            )?);
        outcome.muv_invariance_passed = outcome.muv_rows.iter().all(|row| row.passed);
        let shifted_rls_name = shifted_rls_integrand_name(case.integrand_name);
        outcome.rls_rows = inspect_scale_dependence_row(
            &mut cli,
            case,
            case.integrand_name,
            "rls.inspect",
            &shifted_rls_name,
            true,
        )?;
        outcome
            .rls_rows
            .extend(integrated_result_scale_invariance_rows(
                &baseline,
                "rls.integrated",
                &baseline.integrated_rls_shifted_per_graph,
            )?);
        outcome.rls_invariance_passed = outcome.rls_rows.iter().all(|row| row.passed);

        if case.check_mu_r_dependence {
            let shifted_mur_name = shifted_mur_integrand_name(case.integrand_name);
            outcome.mur_rows = inspect_scale_dependence_row(
                &mut cli,
                case,
                case.integrand_name,
                "mu_r.inspect",
                &shifted_mur_name,
                false,
            )?;
            outcome.mur_dependence_passed = Some(outcome.mur_rows.iter().all(|row| row.passed));
        }

        if case.targets.any() {
            outcome.target_passed = Some(integrated_uv_targets_pass(&baseline, &case.targets)?);
        }

        Ok(())
    })();

    if let Err(err) = case_result {
        outcome.error = Some(format!("{err:#}"));
    }

    // clean_test(&cli.cli_settings.state.folder);
    outcome
}

fn print_integrated_uv_summary(results: &[IntegratedUvCaseResult]) {
    #[derive(Tabled)]
    struct SummaryRow {
        case: String,
        check: String,
        graph: String,
        status: String,
        sigma: String,
        threshold: String,
        error: String,
    }

    fn status(value: Option<bool>) -> String {
        match value {
            Some(true) => "pass".to_string(),
            Some(false) => "FAIL".to_string(),
            None => "-".to_string(),
        }
    }

    fn sigma(value: Option<f64>) -> String {
        value.map_or_else(|| "-".to_string(), |value| format!("{value:.1}sigma"))
    }

    let mut rows = Vec::new();
    for result in results {
        let error = result.error.clone().unwrap_or_else(|| "-".to_string());
        rows.push(SummaryRow {
            case: result.graph.clone(),
            check: "uv_profile".to_string(),
            graph: "-".to_string(),
            status: status(result.uv_profile_passed),
            sigma: "-".to_string(),
            threshold: "-".to_string(),
            error: error.clone(),
        });
        rows.extend(result.muv_rows.iter().map(|row| SummaryRow {
            case: result.graph.clone(),
            check: row.check.to_string(),
            graph: row.graph.clone(),
            status: status(Some(row.passed)),
            sigma: sigma(row.sigma),
            threshold: row.threshold.clone().unwrap_or_else(|| "-".to_string()),
            error: error.clone(),
        }));
        if result.muv_rows.is_empty() {
            rows.push(SummaryRow {
                case: result.graph.clone(),
                check: "m_uv".to_string(),
                graph: "-".to_string(),
                status: status(Some(result.muv_invariance_passed)),
                sigma: "-".to_string(),
                threshold: "<=2σ".to_string(),
                error: error.clone(),
            });
        }
        rows.extend(result.rls_rows.iter().map(|row| SummaryRow {
            case: result.graph.clone(),
            check: row.check.to_string(),
            graph: row.graph.clone(),
            status: status(Some(row.passed)),
            sigma: sigma(row.sigma),
            threshold: row.threshold.clone().unwrap_or_else(|| "-".to_string()),
            error: error.clone(),
        }));
        if result.rls_rows.is_empty() {
            rows.push(SummaryRow {
                case: result.graph.clone(),
                check: "rls".to_string(),
                graph: "-".to_string(),
                status: status(Some(result.rls_invariance_passed)),
                sigma: "-".to_string(),
                threshold: "<=2σ".to_string(),
                error: error.clone(),
            });
        }
        rows.extend(result.mur_rows.iter().map(|row| SummaryRow {
            case: result.graph.clone(),
            check: row.check.to_string(),
            graph: row.graph.clone(),
            status: status(Some(row.passed)),
            sigma: sigma(row.sigma),
            threshold: row.threshold.clone().unwrap_or_else(|| "-".to_string()),
            error: error.clone(),
        }));
        if result.mur_rows.is_empty() {
            rows.push(SummaryRow {
                case: result.graph.clone(),
                check: "mu_r".to_string(),
                graph: "-".to_string(),
                status: status(result.mur_dependence_passed),
                sigma: "-".to_string(),
                threshold: result.mur_dependence_passed.map_or_else(
                    || "-".to_string(),
                    |_| {
                        format!(
                            "relative delta >= {INSPECT_DEPENDENCE_ACCURACY_FACTOR:.0}*accuracy"
                        )
                    },
                ),
                error: error.clone(),
            });
        }
        rows.extend(
            result
                .integrated_ct_accuracy_rows
                .iter()
                .map(|row| SummaryRow {
                    case: result.graph.clone(),
                    check: row.check.to_string(),
                    graph: row.graph.clone(),
                    status: status(Some(row.passed)),
                    sigma: sigma(row.sigma),
                    threshold: row.threshold.clone().unwrap_or_else(|| "-".to_string()),
                    error: error.clone(),
                }),
        );
        rows.push(SummaryRow {
            case: result.graph.clone(),
            check: "target".to_string(),
            graph: "-".to_string(),
            status: status(result.target_passed),
            sigma: "-".to_string(),
            threshold: result.target_threshold.unwrap_or("-").to_string(),
            error: error.clone(),
        });
    }

    println!("{}", Table::new(rows));
}

fn run_single_integrated_uv_case(case: &IntegratedUvCase<'_>) {
    let result = run_integrated_uv_case(case);
    print_integrated_uv_summary(std::slice::from_ref(&result));
    assert!(
        result.all_passed(),
        "Integrated UV case failed: {}",
        result.graph
    );
}

#[test]
fn dod0_bubble_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/dod0_bubble",
        test_name: "dod0_bubble",
        process: "bubble",
        integrand_name: "scalar_bubble_below_thres",
        integrator: DEFAULT_INTEGRATED_UV_INTEGRATOR,
        original_m_uv: 0.2,
        shifted_m_uv: 0.5,
        original_renormalization_localization_scale: 0.2,
        shifted_renormalization_localization_scale: 0.5,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: false,
        targets: IntegratedUvTargets {
            integrated: Some(Complex::new(F(0.0), F(0.01451155120018305))),
        },
        integrated_ct_relative_error_limit: None,
        check_mu_r_dependence: true,
    });
}

#[test]
fn dod1_bubble_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/dod1_bubble",
        test_name: "dod1_bubble",
        process: "bubble_dod1",
        integrand_name: "scalar_bubble_below_thres",
        integrator: DEFAULT_INTEGRATED_UV_INTEGRATOR,
        original_m_uv: 0.2,
        shifted_m_uv: 0.5,
        original_renormalization_localization_scale: 0.2,
        shifted_renormalization_localization_scale: 0.5,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: false,
        targets: IntegratedUvTargets {
            integrated: Some(Complex::new(F(0.0), F(0.003628430563793077))),
        },
        integrated_ct_relative_error_limit: None,
        check_mu_r_dependence: true,
    });
}

#[test]
fn dod2_bubble_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/dod2_bubble",
        test_name: "dod2_bubble",
        process: "bubble_dod2",
        integrand_name: "scalar_bubble_below_thres",
        integrator: DEFAULT_INTEGRATED_UV_INTEGRATOR,
        original_m_uv: 0.2,
        shifted_m_uv: 0.5,
        original_renormalization_localization_scale: 0.2,
        shifted_renormalization_localization_scale: 0.5,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: false,
        targets: IntegratedUvTargets {
            integrated: Some(Complex::new(F(0.0), F(0.01234018404957018))),
        },
        integrated_ct_relative_error_limit: None,
        check_mu_r_dependence: true,
    });
}

#[test]
fn se1l_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/se1l",
        test_name: "se1l",
        process: "se1l",
        integrand_name: "se1l",
        integrator: DEFAULT_INTEGRATED_UV_INTEGRATOR,
        original_m_uv: 0.2,
        shifted_m_uv: 0.5,
        original_renormalization_localization_scale: 0.2,
        shifted_renormalization_localization_scale: 0.5,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: false,
        targets: IntegratedUvTargets {
            integrated: Some(Complex::new(F(53.22768383202566), F(-63915.80926056468))),
        },
        integrated_ct_relative_error_limit: None,
        check_mu_r_dependence: true,
    });
}

#[test]
fn sunrise_scalar_1_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/sunrise_scalar_1",
        test_name: "sunrise_scalar_1",
        process: "sunrise_scalar_1",
        integrand_name: "scalar_sunrise",
        integrator: SUNRISE_INTEGRATED_UV_INTEGRATOR,
        original_m_uv: 0.2,
        shifted_m_uv: 0.5,
        original_renormalization_localization_scale: 0.2,
        shifted_renormalization_localization_scale: 0.5,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: false,
        targets: IntegratedUvTargets {
            integrated: Some(Complex::new(F(-0.000049944992589155517), F(0.0))),
        },
        integrated_ct_relative_error_limit: Some(INTEGRATED_CT_RELATIVE_ERROR_LIMIT),
        check_mu_r_dependence: true,
    });
}

#[test]
fn dotted_dt_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/dotted_dt",
        test_name: "dotted_dt",
        process: "dotted_dt",
        integrand_name: "dotted_dt",
        integrator: DEFAULT_INTEGRATED_UV_INTEGRATOR,
        original_m_uv: 0.2,
        shifted_m_uv: 0.5,
        original_renormalization_localization_scale: 0.2,
        shifted_renormalization_localization_scale: 0.5,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: false,
        targets: IntegratedUvTargets::default(),
        integrated_ct_relative_error_limit: None,
        check_mu_r_dependence: true,
    });
}

#[test]
fn epem_a_bbx_amp_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/epem_a_bbx_amp",
        test_name: "epem_a_bbx_amp",
        process: "epem_a_bbx",
        integrand_name: "epem_a_bbx",
        integrator: DEFAULT_INTEGRATED_UV_INTEGRATOR,
        original_m_uv: 0.2,
        shifted_m_uv: 0.5,
        original_renormalization_localization_scale: 0.2,
        shifted_renormalization_localization_scale: 0.5,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: false,
        targets: IntegratedUvTargets::default(),
        integrated_ct_relative_error_limit: None,
        check_mu_r_dependence: true,
    });
}

const AA_AA_2L_GRAPHS: &[&str] = &[
    "GL00", "GL01", "GL03", "GL05", "GL06", "GL08", "GL12", "GL14", "GL16", "GL17", "GL18",
];

const AA_AA_2L_UV_RICH_INSPECT: GraphUvRichInspectCase = GraphUvRichInspectCase {
    name: "aa_aa 2L",
    test_prefix: "aa_aa_2l",
    process_prefix: "aa_aa_2l_",
    integrand_name: "2L",
    run_card: "uv/aa_aa_2l_rich_inspect.toml",
    graph_block_prefix: "generate_",
    target_file: "aa_aa_2l_uv_inspect_event_targets.json",
    graphs: AA_AA_2L_GRAPHS,
};

#[derive(Clone, Copy)]
struct UvRichInspectMode {
    name: &'static str,
    process_suffix: &'static str,
}

const INTEGRATED_UV_RICH_INSPECT_MODE: UvRichInspectMode = UvRichInspectMode {
    name: "integrated",
    process_suffix: "",
};

const UV_RICH_INSPECT_MODES: &[UvRichInspectMode] = &[INTEGRATED_UV_RICH_INSPECT_MODE];

struct GraphUvRichInspectCase {
    name: &'static str,
    test_prefix: &'static str,
    process_prefix: &'static str,
    integrand_name: &'static str,
    run_card: &'static str,
    graph_block_prefix: &'static str,
    target_file: &'static str,
    graphs: &'static [&'static str],
}

impl GraphUvRichInspectCase {
    fn target_path(&self) -> PathBuf {
        workspace_root()
            .join("tests/resources/benchmarks")
            .join(self.target_file)
    }

    fn graph_test_name(&self, graph: &str) -> String {
        format!(
            "{}_{}_uv_profile_and_rich_inspect",
            self.test_prefix,
            graph.to_ascii_lowercase()
        )
    }

    fn process_name(&self, graph: &str, mode: UvRichInspectMode) -> String {
        format!(
            "{}{}{}",
            self.process_prefix,
            graph.to_ascii_lowercase(),
            mode.process_suffix
        )
    }

    fn graph_block_name(&self, graph: &str) -> String {
        format!("{}{}", self.graph_block_prefix, graph.to_ascii_lowercase())
    }

    fn setup_graph_cli(&self, graph: &str, test_name: &str) -> Result<CLIState> {
        let mut cli = get_test_cli(
            Some(self.run_card.into()),
            get_tests_workspace_path().join(test_name),
            Some(test_name.to_string()),
            true,
        )?;
        cli.run_command(&format!("run {}", self.graph_block_name(graph)))?;

        Ok(cli)
    }

    fn momentum_point(&self, cli: &CLIState, process: &str) -> Result<Vec<f64>> {
        let process_id = cli
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
        let integrand = cli
            .state
            .process_list
            .get_integrand(process_id, self.integrand_name)?
            .require_generated()?;
        let seed = [0.11, -0.07, 0.19, -0.13, 0.05, 0.29];
        Ok((0..integrand.get_n_dim())
            .map(|index| seed[index % seed.len()])
            .collect())
    }

    fn graph_name(&self, cli: &CLIState, process: &str) -> Result<String> {
        let process_id = cli
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified(process.to_string())))?;
        Ok(cli
            .state
            .process_list
            .get_integrand(process_id, self.integrand_name)?
            .require_generated()?
            .graph_name_by_id(0)
            .ok_or_else(|| eyre::eyre!("Generated {} graph should have graph id 0", self.name))?
            .to_string())
    }

    fn rich_inspect_record(
        &self,
        cli: &mut CLIState,
        graph: &str,
        mode: UvRichInspectMode,
        point: &[f64],
    ) -> Result<Value> {
        let process = self.process_name(graph, mode);
        let process_id = cli
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified(process.clone())))?;
        let graph_name = self.graph_name(cli, &process)?;
        let points = Array2::from_shape_vec((1, point.len()), point.to_vec())?;
        let result = evaluate_sample(
            &mut cli.state,
            &EvaluateSamples {
                process_id: Some(process_id),
                integrand_name: Some(self.integrand_name.to_string()),
                use_arb_prec: false,
                minimal_output: false,
                return_generated_events: Some(true),
                momentum_space: true,
                points: points.view(),
                integrator_weights: None,
                discrete_dims: None,
                graph_names: Some(vec![Some(graph_name)]),
                orientations: None,
            },
        )?;

        let evaluation = result.sample.evaluation;
        let metadata = evaluation
            .evaluation_metadata
            .as_ref()
            .expect("rich inspect evaluation should include metadata");
        let mut event_weight_sum = Complex::new(F(0.0), F(0.0));
        let mut additional_weights = BTreeMap::<String, Complex<F<f64>>>::new();
        let event_group_sizes = evaluation
            .event_groups
            .iter()
            .map(|group| {
                for event in group.iter() {
                    event_weight_sum += event.weight;
                    for (key, value) in &event.additional_weights.weights {
                        *additional_weights
                            .entry(additional_weight_key_label(*key))
                            .or_insert_with(|| Complex::new(F(0.0), F(0.0))) += *value;
                    }
                }
                group.len()
            })
            .collect::<Vec<_>>();

        let mut additional_weight_map = Map::new();
        for (key, value) in additional_weights {
            additional_weight_map.insert(key, complex_json(value));
        }

        Ok(json!({
            "graph": graph,
            "mode": mode.name,
            "result": complex_json(evaluation.integrand_result),
            "metadata": {
                "generated_event_count": metadata.generated_event_count,
                "accepted_event_count": metadata.accepted_event_count,
                "is_nan": metadata.is_nan,
            },
            "event_group_sizes": event_group_sizes,
            "event_weight_sum": complex_json(event_weight_sum),
            "additional_weights": Value::Object(additional_weight_map),
        }))
    }

    fn target_key(&self, record: &Value) -> Result<(String, String)> {
        let graph = record.get("graph").and_then(Value::as_str).ok_or_else(|| {
            eyre::eyre!(
                "{} target record is missing string field 'graph'",
                self.name
            )
        })?;
        let mode = record.get("mode").and_then(Value::as_str).ok_or_else(|| {
            eyre::eyre!("{} target record is missing string field 'mode'", self.name)
        })?;
        Ok((graph.to_string(), mode.to_string()))
    }

    fn load_targets(&self) -> Result<BTreeMap<(String, String), Value>> {
        let target_path = self.target_path();
        let contents = fs::read_to_string(&target_path)?;
        let records: Vec<Value> = serde_json::from_str(&contents)?;
        let mut targets = BTreeMap::new();
        for record in records {
            let key = self.target_key(&record)?;
            if targets.insert(key.clone(), record).is_some() {
                return Err(eyre::eyre!(
                    "Duplicate {} target record for graph={} mode={}",
                    self.name,
                    key.0,
                    key.1
                ));
            }
        }
        Ok(targets)
    }

    fn assert_target_matches(&self, record: &Value, targets: &BTreeMap<(String, String), Value>) {
        let key = self
            .target_key(record)
            .expect("actual target record should be keyed");
        let target = targets
            .get(&key)
            .unwrap_or_else(|| panic!("Missing {} inspect target for {key:?}", self.name));
        assert_json_approx_eq(record, target, &format!("{} {}", key.0, key.1));
    }

    fn run_graph_uv_profile_and_rich_inspect(&self, graph: &str) -> Result<()> {
        let test_name = self.graph_test_name(graph);
        let mut cli = self.setup_graph_cli(graph, &test_name)?;
        let targets = self.load_targets()?;
        let point_process = self.process_name(graph, INTEGRATED_UV_RICH_INSPECT_MODE);
        let point = self.momentum_point(&cli, &point_process)?;

        for mode in UV_RICH_INSPECT_MODES.iter().copied() {
            let process = self.process_name(graph, mode);
            assert!(
                integrated_uv_profile_passes(&mut cli, &process, self.integrand_name)?,
                "UV profile failed for case={} graph={graph} mode={}",
                self.name,
                mode.name
            );
            let record = self.rich_inspect_record(&mut cli, graph, mode, &point)?;
            self.assert_target_matches(&record, &targets);
        }

        clean_test(&cli.cli_settings.state.folder);
        Ok(())
    }

    fn collect_rich_inspect_records(&self) -> Result<Vec<Value>> {
        let mut records = Vec::new();
        for graph in self.graphs {
            let test_name = format!("{}_target_writer", self.graph_test_name(graph));
            let mut cli = self.setup_graph_cli(graph, &test_name)?;
            let point_process = self.process_name(graph, INTEGRATED_UV_RICH_INSPECT_MODE);
            let point = self.momentum_point(&cli, &point_process)?;
            for mode in UV_RICH_INSPECT_MODES.iter().copied() {
                records.push(self.rich_inspect_record(&mut cli, graph, mode, &point)?);
            }
            clean_test(&cli.cli_settings.state.folder);
        }
        Ok(records)
    }
}

fn complex_json(value: Complex<F<f64>>) -> Value {
    json!({
        "re": value.re.0,
        "im": value.im.0,
    })
}

fn additional_weight_key_label(key: AdditionalWeightKey) -> String {
    match key {
        AdditionalWeightKey::FullMultiplicativeFactor => "full_multiplicative_factor".to_string(),
        AdditionalWeightKey::Original => "original".to_string(),
        AdditionalWeightKey::ThresholdCounterterm { subset_index } => {
            format!("threshold_counterterm_{subset_index}")
        }
        AdditionalWeightKey::AmplitudeThresholdCounterterm {
            esurface_id,
            overlap_group,
        } => format!("ct_{esurface_id}_{overlap_group}"),
    }
}

fn assert_json_approx_eq(actual: &Value, expected: &Value, path: &str) {
    match (actual, expected) {
        (Value::Number(actual), Value::Number(expected)) => {
            let actual = actual
                .as_f64()
                .unwrap_or_else(|| panic!("{path}: actual number is not representable as f64"));
            let expected = expected
                .as_f64()
                .unwrap_or_else(|| panic!("{path}: expected number is not representable as f64"));
            let scale = actual.abs().max(expected.abs()).max(1.0);
            let tolerance = 1.0e-10 * scale;
            assert!(
                (actual - expected).abs() <= tolerance,
                "{path}: actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.17e}"
            );
        }
        (Value::Array(actual), Value::Array(expected)) => {
            assert_eq!(
                actual.len(),
                expected.len(),
                "{path}: array length mismatch; actual={actual:?}, expected={expected:?}"
            );
            for (index, (actual, expected)) in actual.iter().zip(expected).enumerate() {
                assert_json_approx_eq(actual, expected, &format!("{path}[{index}]"));
            }
        }
        (Value::Object(actual), Value::Object(expected)) => {
            let actual_keys = actual.keys().collect::<Vec<_>>();
            let expected_keys = expected.keys().collect::<Vec<_>>();
            assert_eq!(
                actual_keys, expected_keys,
                "{path}: object keys mismatch; actual={actual_keys:?}, expected={expected_keys:?}"
            );
            for (key, expected_value) in expected {
                let actual_value = actual
                    .get(key)
                    .unwrap_or_else(|| panic!("{path}.{key}: missing actual value"));
                assert_json_approx_eq(actual_value, expected_value, &format!("{path}.{key}"));
            }
        }
        _ => assert_eq!(actual, expected, "{path}: JSON value mismatch"),
    }
}

mod slow {
    use super::*;

    macro_rules! aa_aa_2l_uv_rich_inspect_tests {
        ($($name:ident => $graph:literal),+ $(,)?) => {
            $(
                #[test]
                #[serial_test::serial]
                fn $name() -> Result<()> {
                    AA_AA_2L_UV_RICH_INSPECT.run_graph_uv_profile_and_rich_inspect($graph)
                }
            )+
        };
    }

    aa_aa_2l_uv_rich_inspect_tests! {
        aa_aa_2l_gl00_uv_profile_and_rich_inspect => "GL00",
        aa_aa_2l_gl01_uv_profile_and_rich_inspect => "GL01",
        aa_aa_2l_gl03_uv_profile_and_rich_inspect => "GL03",
        aa_aa_2l_gl05_uv_profile_and_rich_inspect => "GL05",
        aa_aa_2l_gl06_uv_profile_and_rich_inspect => "GL06",
        aa_aa_2l_gl08_uv_profile_and_rich_inspect => "GL08",
        aa_aa_2l_gl12_uv_profile_and_rich_inspect => "GL12",
        aa_aa_2l_gl14_uv_profile_and_rich_inspect => "GL14",
        aa_aa_2l_gl16_uv_profile_and_rich_inspect => "GL16",
        aa_aa_2l_gl17_uv_profile_and_rich_inspect => "GL17",
        aa_aa_2l_gl18_uv_profile_and_rich_inspect => "GL18",
    }

    #[test]
    #[ignore = "target writer"]
    #[serial_test::serial]
    fn write_aa_aa_2l_uv_inspect_event_targets() -> Result<()> {
        let records = AA_AA_2L_UV_RICH_INSPECT.collect_rich_inspect_records()?;
        fs::write(
            AA_AA_2L_UV_RICH_INSPECT.target_path(),
            format!("{}\n", serde_json::to_string_pretty(&records)?),
        )?;
        Ok(())
    }

    #[test]
    fn aa_aa_gl00_uv() {
        run_single_integrated_uv_case(&IntegratedUvCase {
            run_card: "uv/aa_aa_GL00",
            test_name: "aa_aa_GL00",
            process: "aa_aa",
            integrand_name: "2L",
            integrator: DEFAULT_INTEGRATED_UV_INTEGRATOR,
            original_m_uv: 91.188,
            shifted_m_uv: 364.752,
            original_renormalization_localization_scale: 5.0,
            shifted_renormalization_localization_scale: 25.0,
            original_mu_r: 91.188,
            shifted_mu_r: 364.752,
            skip_uv_profile: false,
            targets: IntegratedUvTargets::default(),
            integrated_ct_relative_error_limit: None,
            check_mu_r_dependence: true,
        });
    }

    #[test]
    fn epem_ttxh_gl00_uv() {
        run_single_integrated_uv_case(&IntegratedUvCase {
            run_card: "uv/epem_ttxh_GL00",
            test_name: "epem_ttxh_gl00",
            process: "epem_a_tth",
            integrand_name: "NLO",
            integrator: DEFAULT_INTEGRATED_UV_INTEGRATOR,
            original_m_uv: 20.0,
            shifted_m_uv: 7.0,
            original_renormalization_localization_scale: 5.0,
            shifted_renormalization_localization_scale: 25.0,
            original_mu_r: 3.0,
            shifted_mu_r: 9.0,
            skip_uv_profile: false,
            targets: IntegratedUvTargets::default(),
            integrated_ct_relative_error_limit: None,
            check_mu_r_dependence: true,
        });
    }

    #[test]
    fn ad_ad_with_gluon_correction_uv() {
        run_single_integrated_uv_case(&IntegratedUvCase {
            run_card: "uv/ad_ad_with_gluon_correction",
            test_name: "ad_ad_with_gluon_correction",
            process: "adad",
            integrand_name: "adad_gluon",
            integrator: DEFAULT_INTEGRATED_UV_INTEGRATOR,
            original_m_uv: 20.0,
            shifted_m_uv: 7.0,
            original_renormalization_localization_scale: 5.0,
            shifted_renormalization_localization_scale: 25.0,
            original_mu_r: 20.0,
            shifted_mu_r: 9.0,
            skip_uv_profile: false,
            targets: IntegratedUvTargets::default(),
            integrated_ct_relative_error_limit: None,
            check_mu_r_dependence: true,
        });
    }

    #[test]
    fn soft_ct_se() -> Result<()> {
        let mut state = get_test_cli(
            Some("dgse.toml".into()),
            get_tests_workspace_path().join("dgse"),
            None,
            false,
        )?;

        let res = Profile::InfraRed(InfraRedProfile {
            select: Some("se S(e0)".into()),
            ..Default::default()
        })
        .run(&mut state.state, &state.cli_settings)?;

        assert!(res.unwrap_ir().all_passed);

        let res = Profile::UltraViolet(UltraVioletProfile {
            ..Default::default()
        })
        .run(&mut state.state, &state.cli_settings)?;

        let uv = res.unwrap_uv();
        assert_eq!(uv.pass_fail(-0.9).failed, 0);
        Ok(())
    }
}
#[test]
fn epem_a_ddx_xs_nlo_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/epem_a_ddx_xs_nlo",
        test_name: "epem_a_ddx_xs_nlo",
        process: "epem_a_ddx",
        integrand_name: "NLO",
        integrator: DEFAULT_INTEGRATED_UV_INTEGRATOR,
        original_m_uv: 20.0,
        shifted_m_uv: 7.0,
        original_renormalization_localization_scale: 5.0,
        shifted_renormalization_localization_scale: 25.0,
        original_mu_r: 3.0,
        shifted_mu_r: 9.0,
        skip_uv_profile: true,
        targets: IntegratedUvTargets {
            integrated: Some(Complex::new(F(0.0), F(1.163e-3))),
        },
        integrated_ct_relative_error_limit: None,
        check_mu_r_dependence: false,
    });
}

#[test]
fn epem_a_tth_nlo_uv() -> Result<()> {
    let state_path = get_tests_workspace_path().join("epem_a_tth_nlo_example");
    let mut cli = get_example_cli(
        "epem_a_ttxh/NLO/epem_a_tth_NLO.toml",
        &["generate_diagrams"],
        Some(state_path.clone()),
        None,
        true,
    )?;

    cli.run_command("run generate_diagrams generate_integrands")?;
    let res = Profile::UltraViolet(UltraVioletProfile {
        process: Some(ProcessRef::Unqualified("epem_a_tth".to_string())),
        integrand_name: Some("NLO".to_string()),
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?;

    let uv = res.unwrap_uv();
    assert_eq!(uv.pass_fail(-0.9).failed, 0);
    clean_test(&state_path);
    Ok(())
}
