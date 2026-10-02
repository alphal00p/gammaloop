use std::{
    collections::{BTreeMap, BTreeSet},
    env, fs,
    panic::{AssertUnwindSafe, catch_unwind},
    path::PathBuf,
    process::Command,
    time::{Duration, Instant},
};

use color_eyre::Result;
use gammaloop_api::commands::{
    Commands, Profile,
    evaluate_samples::{
        EvaluateSamples, EvaluateSamplesPrecise, evaluate_sample, evaluate_sample_precise,
    },
    integrate::Integrate,
    profile::{InfraRedProfile, UltraVioletProfile},
};
use gammaloop_api::{
    StateLoadOption,
    state::{CommandHistory, ProcessRef},
};
use gammaloop_integration_tests::{
    CLIState, clean_test, get_test_cli, get_tests_workspace_path, workspace_root,
};
use gammalooprs::integrands::{
    HasIntegrand,
    evaluation::{EvaluationMetaData, PreciseEvaluationResultOutput},
};
use gammalooprs::observables::events::AdditionalWeightKey;
use gammalooprs::processes::ProcessCollection;
use gammalooprs::settings::runtime::{IntegralEstimate, SlotIntegrationResult};
use gammalooprs::utils::{ArbPrec, F, FloatLike};
use gammalooprs::uv::{
    ApproximationType, CTIdentifier, CTRenormalizationRule, UVOrchestrator, UVProfileAnalysis,
    UltravioletGraph, settings::FinalIntegrandDimension,
};
use ndarray::Array2;
use serde_json::{Map, Value, json};
use spenso::algebra::{algebraic_traits::IsZero, complex::Complex};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::float::Real,
};
use tabled::{Table, Tabled};

const INSPECT_DEPENDENCE_ACCURACY_FACTOR: f64 = 1000.0;
const INSPECT_DEPENDENCE_SCALE_EXPONENTS: [f64; 5] = [
    0.0,
    1.0,
    2.0,
    UV_PROFILE_MIN_SCALE_EXPONENT,
    UV_PROFILE_MAX_SCALE_EXPONENT,
];
const INSPECT_DEPENDENCE_SEEDS: [[f64; 6]; 2] = [
    [0.11, -0.07, 0.19, -0.13, 0.05, 0.29],
    [0.23, 0.31, -0.17, 0.41, -0.29, 0.13],
];
const UV_PROFILE_MIN_SCALE_EXPONENT: f64 = 8.0;
const UV_PROFILE_MAX_SCALE_EXPONENT: f64 = 12.0;
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
        n_cores: 10,
    };

const EPEM_A_BBX_INTEGRATED_UV_INTEGRATOR: IntegratedUvIntegratorSettings =
    IntegratedUvIntegratorSettings {
        target_relative_accuracy: 0.01,
        n_start: 50_000,
        n_increase: 0,
        n_max: 400_000,
        n_cores: 1,
    };

const DOTTED_DT_INTEGRATED_UV_INTEGRATOR: IntegratedUvIntegratorSettings =
    IntegratedUvIntegratorSettings {
        target_relative_accuracy: 0.04,
        n_start: 10_000,
        n_increase: 0,
        n_max: 50_000,
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
    fn threshold_ratio(&self) -> Result<f64> {
        let ratio = self.relative_delta / self.required_relative_delta;
        if !self.relative_delta.is_finite()
            || self.relative_delta < 0.0
            || !self.relative_accuracy.is_finite()
            || self.relative_accuracy <= 0.0
            || !self.required_relative_delta.is_finite()
            || self.required_relative_delta <= 0.0
            || !ratio.is_finite()
        {
            return Err(eyre::eyre!(
                "invalid or non-finite inspect probe at exponent={:+.1}, seed={}: relative delta={}, accuracy={}, required delta={}",
                self.probe.scale_exponent,
                self.probe.seed_index,
                self.relative_delta,
                self.relative_accuracy,
                self.required_relative_delta,
            ));
        }
        Ok(ratio)
    }

    fn passes(&self, allow_absent_dependence: bool) -> Result<bool> {
        Ok(self.threshold_ratio()? >= 1.0
            || (allow_absent_dependence && self.relative_delta <= self.relative_accuracy))
    }
}

#[test]
fn inspect_probe_dependence_uses_reported_accuracy_and_rejects_nonfinite_probes() {
    let probe = |relative_delta, relative_accuracy| InspectProbeResult {
        probe: InspectProbePoint {
            point: vec![1.0],
            scale_exponent: 0.0,
            seed_index: 0,
        },
        relative_delta,
        relative_accuracy,
        required_relative_delta: INSPECT_DEPENDENCE_ACCURACY_FACTOR * relative_accuracy,
    };
    let accuracy = 1.0e-10;
    let threshold = INSPECT_DEPENDENCE_ACCURACY_FACTOR * accuracy;
    for (delta, allow_absence, expected) in [
        (0.0, true, true),
        (f64::EPSILON, true, true),
        (accuracy, true, true),
        (10.0 * accuracy, true, false),
        (0.0, false, false),
        (f64::EPSILON, false, false),
        (0.99 * threshold, false, false),
        (threshold, false, true),
        (threshold, true, true),
    ] {
        assert_eq!(
            probe(delta, accuracy).passes(allow_absence).unwrap(),
            expected
        );
    }
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(probe(invalid, accuracy).passes(true).is_err());
        assert!(probe(0.0, invalid).passes(true).is_err());
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

fn evaluate_momentum_sample(
    cli: &mut CLIState,
    process: &str,
    integrand_name: &str,
    point: &[f64],
    orientation: Option<usize>,
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
            orientations: Some(vec![orientation]),
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
    let mut best: Option<(f64, InspectProbeResult)> = None;
    for probe in probes {
        let baseline = evaluate_momentum_sample(
            cli,
            case.process,
            baseline_name,
            &probe.point,
            None,
            accuracy_floor,
        )?;
        let shifted = evaluate_momentum_sample(
            cli,
            case.process,
            shifted_name,
            &probe.point,
            None,
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
        // Validate every probe before ranking so a non-finite result cannot be
        // silently discarded in favor of a finite sample.
        let ratio = probe_result.threshold_ratio()?;
        if match &best {
            Some((current, _)) => ratio > *current,
            None => true,
        } {
            best = Some((ratio, probe_result));
        }
    }

    let (_, best) =
        best.ok_or_else(|| eyre::eyre!("{check} should have at least one inspect probe"))?;
    // The maximum delta/accuracy ratio bounds every validated probe. Allowed
    // absence is a numerical statement at the reported accuracy, not bit equality.
    let absent_dependence_allowed =
        allow_absent_dependence && best.relative_delta <= best.relative_accuracy;
    let passed = best.passes(allow_absent_dependence)?;
    let threshold = if absent_dependence_allowed {
        format!(
            "all {n_probes} inspect probes compatible with absence: best relative delta {:.3e} <= accuracy {:.3e} (allowed; integrated invariance checked separately)",
            best.relative_delta, best.relative_accuracy,
        )
    } else {
        format!(
            "best inspect relative delta {:.3e}; required >= {INSPECT_DEPENDENCE_ACCURACY_FACTOR:.0}*accuracy {:.3e} (accuracy={:.3e}, exponent={:+.1}, seed={}, probes={n_probes})",
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

fn squared_integral_estimate(estimate: &IntegralEstimate) -> IntegralEstimate {
    let re = estimate.result.re.0;
    let im = estimate.result.im.0;
    let err_re = estimate.error.re.0;
    let err_im = estimate.error.im.0;

    IntegralEstimate {
        neval: estimate.neval,
        real_zero: estimate.real_zero,
        im_zero: estimate.im_zero,
        result: Complex::new(F(re * re - im * im), F(2.0 * re * im)),
        error: Complex::new(
            F((2.0 * re * err_re).hypot(2.0 * im * err_im)),
            F((2.0 * im * err_re).hypot(2.0 * re * err_im)),
        ),
        real_chisq: estimate.real_chisq,
        im_chisq: estimate.im_chisq,
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

#[test]
fn orientation_local_dod_bubbles_are_uv_finite_per_orientation() -> Result<()> {
    for (run_card, test_name, process, expected_initial_dod, seed) in [
        ("uv/dod0_bubble", "dod0", "bubble", 0, 9300),
        ("uv/dod1_bubble", "dod1", "bubble_dod1", 1, 9301),
        ("uv/dod2_bubble", "dod2", "bubble_dod2", 2, 9302),
    ] {
        let mut cli = get_test_cli(
            Some(format!("{run_card}.toml").into()),
            get_tests_workspace_path()
                .join("orientation_local_dod_bubbles")
                .join(test_name),
            None,
            true,
        )?;
        cli.run_command("run generate")?;

        let generation = &cli.cli_settings.global.generation;
        assert!(!generation.explicit_orientation_sum_only);
        assert!(!generation.uv.local_uv_cts_from_expanded_4d_integrands);
        assert!(generation.uv.generate_integrated);
        assert!(generation.threshold_subtraction.enable_thresholds);

        let analysis = Profile::UltraViolet(UltraVioletProfile {
            process: Some(ProcessRef::Unqualified(process.to_string())),
            integrand_name: Some("scalar_bubble_below_thres".to_string()),
            min_scale_exponent: 4.0,
            max_scale_exponent: 8.0,
            n_points: 6,
            per_orientation: true,
            seed: Some(seed),
            uv_ray_directions: vec![1.0, 0.3, -0.2],
            uv_ray_norms: vec![3.0],
            ..Default::default()
        })
        .run(&mut cli.state, &cli.cli_settings)?
        .unwrap_uv();

        let subsets = analysis
            .graphs
            .iter()
            .flat_map(|graph| &graph.lmbs)
            .flat_map(|lmb| &lmb.subsets)
            .collect::<Vec<_>>();
        assert_eq!(subsets.len(), 1, "{test_name} selected UV limits");
        assert_eq!(subsets[0].initial_dod, expected_initial_dod);
        let fitted_orientations = subsets
            .iter()
            .flat_map(|subset| subset.per_orientation_inspect_entries.iter().flatten())
            .filter(|entry| entry.analysis.is_some())
            .count();
        let orientation_failures = analysis
            .pass_fail(-0.9)
            .failures
            .into_iter()
            .filter(|failure| failure.orientation_label.is_some())
            .collect::<Vec<_>>();

        clean_test(&cli.cli_settings.state.folder);
        assert!(
            fitted_orientations > 0,
            "{test_name} per-orientation UV profile was vacuous"
        );
        assert!(
            orientation_failures.is_empty(),
            "{test_name} has non-integrable active orientations: {orientation_failures:#?}"
        );
    }

    Ok(())
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
        cli.cli_settings.global.generation.uv.orchestrator = UVOrchestrator::Compare;
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
fn epem_a_ddx_nlo_raised_cff_generation() -> Result<()> {
    let mut cli = get_test_cli(
        Some("uv/epem_a_ddx_xs_nlo.toml".into()),
        get_tests_workspace_path().join("epem_a_ddx_nlo_raised_cff_generation"),
        None,
        true,
    )?;
    cli.cli_settings.global.generation.orientation_pattern = Default::default();
    cli.run_command("run generate")?;
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial_test::serial]
fn scalar_amplitudes_match_across_local_uv_routes() -> Result<()> {
    const AMPLITUDE_RUNTIME: &str = r#"
        [general]
        integral_unit = "none"
        disable_flux_factor = true
        mu_r = 3.0
        m_uv = 20.0

        [kinematics.externals]
        type = "constant"

        [kinematics.externals.data]
        momenta = [[4.0, 0.0, 0.0, 0.0], "dependent"]
        helicities = [0, 0]

        [sampling]
        graphs = "summed"
        orientations = "summed"
        lmb_multichanneling = true
        lmb_channel_weight = "ose"
        lmb_channels = "summed"

        [subtraction]
        disable_threshold_subtraction = true
    "#;
    let routes = [
        ("localized_local_3d", false, false),
        ("explicit_local_3d", true, false),
        ("projected_local_4d", true, true),
    ];
    let graphs = [
        (
            "dotted_bubble",
            "tests/resources/graphs/dotted_bubble_amp.dot",
        ),
        (
            "double_dotted_bubble",
            "tests/resources/graphs/double_dotted_bubble_amp.dot",
        ),
        (
            "spectacles",
            "tests/resources/graphs/uv_tests/scalar_spectacles_self_energy.dot",
        ),
        (
            "spectacles_overall_log",
            "tests/resources/graphs/uv_tests/scalar_spectacles_overall_log.dot",
        ),
        (
            "cubic_exact_uv_rewrite",
            "tests/resources/graphs/uv_tests/scalar_cubic_exact_uv_rewrite.dot",
        ),
        (
            "mirrored_cubic_exact_uv_rewrite",
            "tests/resources/graphs/uv_tests/scalar_mirrored_cubic_exact_uv_rewrite.dot",
        ),
    ];
    let total_started = Instant::now();
    let mut states = Vec::new();

    for (route, explicit_orientation_sum_only, project_local_4d) in routes {
        let test_root = get_tests_workspace_path().join(format!(
            "scalar_amplitudes_match_across_local_uv_routes_{route}"
        ));
        let mut cli = get_test_cli(None, &test_root, Some(route.to_string()), true)?;
        {
            let generation = &mut cli.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = explicit_orientation_sum_only;
            generation
                .tropical_subgraph_table
                .disable_tropical_generation = true;
            generation.evaluator.iterative_orientation_optimization = false;
            generation.evaluator.compile = false;
            generation.evaluator.store_atom = true;
            generation.threshold_subtraction.enable_thresholds = false;
            generation
                .threshold_subtraction
                .check_esurface_at_generation = false;
            generation.uv.softct = false;
            generation.uv.subtract_uv = true;
            generation.uv.generate_integrated = false;
            generation.uv.final_integrand = FinalIntegrandDimension::ThreeD;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d;
        }
        let generation_started = Instant::now();
        cli.run_command("import model ./assets/models/json/scalars/scalars_2p_3p.json")?;
        cli.run_command("set model mass_scalar_1=1.0")?;
        cli.run_command(&format!("set default-runtime string '{AMPLITUDE_RUNTIME}'"))?;
        for (name, path) in graphs {
            cli.run_command(&format!("import graphs ./{path} -p {name} -i amplitude"))?;
            cli.run_command(&format!("generate existing -p {name} -i amplitude"))?;
        }
        let generation_elapsed = generation_started.elapsed();
        println!("scalar local-UV route {route}: generation {generation_elapsed:?}");
        states.push((route, cli));
    }

    let cases: [(&str, &str, &[f64]); 10] = [
        ("dotted_bubble", "amplitude", &[0.17, -0.31, 0.53]),
        ("double_dotted_bubble", "amplitude", &[0.23, 0.41, -0.67]),
        (
            "spectacles",
            "amplitude",
            &[0.11, -0.29, 0.37, 0.59, -0.43, 0.71],
        ),
        (
            "spectacles_overall_log",
            "amplitude",
            &[0.11, -0.29, 0.37, 0.59, -0.43, 0.71, -0.17, 0.23, 0.41],
        ),
        (
            "cubic_exact_uv_rewrite",
            "amplitude",
            &[0.13, -0.19, 0.31, 0.43, -0.29, 0.61],
        ),
        (
            "cubic_exact_uv_rewrite",
            "amplitude",
            &[-0.37, 0.47, -0.11, 0.73, 0.17, -0.23],
        ),
        (
            "cubic_exact_uv_rewrite",
            "amplitude",
            &[1.3, -1.9, 3.1, 4.3, -2.9, 6.1],
        ),
        (
            "mirrored_cubic_exact_uv_rewrite",
            "amplitude",
            &[0.13, -0.19, 0.31, 0.43, -0.29, 0.61],
        ),
        (
            "mirrored_cubic_exact_uv_rewrite",
            "amplitude",
            &[-0.37, 0.47, -0.11, 0.73, 0.17, -0.23],
        ),
        (
            "mirrored_cubic_exact_uv_rewrite",
            "amplitude",
            &[1.3, -1.9, 3.1, 4.3, -2.9, 6.1],
        ),
    ];
    let mut evaluation_times = vec![Duration::ZERO; states.len()];
    let mut failures = Vec::new();
    for (process, integrand, point) in cases {
        let mut evaluations = Vec::with_capacity(states.len());
        for (index, (route, cli)) in states.iter_mut().enumerate() {
            let evaluation_started = Instant::now();
            let evaluation =
                evaluate_momentum_sample(cli, process, integrand, point, None, 1.0e-12)?;
            evaluation_times[index] += evaluation_started.elapsed();
            if !evaluation.value.re.is_finite() || !evaluation.value.im.is_finite() {
                failures.push(format!(
                    "{process} {route} evaluation is not finite: {:?}",
                    evaluation.value
                ));
            }
            if evaluation.value.re.hypot(evaluation.value.im) == 0.0 {
                failures.push(format!("{process} {route} evaluation is trivially zero"));
            }
            evaluations.push((*route, evaluation));
        }

        let (reference_route, reference) = evaluations
            .iter()
            .find(|(route, _)| *route == "explicit_local_3d")
            .expect("the scalar local-UV matrix must include its explicit local-3D reference");
        for (route, evaluation) in &evaluations {
            if route == reference_route {
                continue;
            }
            let delta = (evaluation.value.re - reference.value.re)
                .hypot(evaluation.value.im - reference.value.im);
            let scale = evaluation
                .value
                .re
                .hypot(evaluation.value.im)
                .max(reference.value.re.hypot(reference.value.im))
                .max(f64::MIN_POSITIVE);
            let tolerance = 1_000.0
                * evaluation
                    .relative_accuracy
                    .max(reference.relative_accuracy)
                    .max(1.0e-12);
            if delta / scale > tolerance {
                failures.push(format!(
                    "{process} {point:?} differs between {route} and {reference_route}: {route}={:?}, {reference_route}={:?}, relative delta={:.3e}, tolerance={tolerance:.3e}",
                    evaluation.value,
                    reference.value,
                    delta / scale,
                ));
            }
        }
    }

    let arb_started = Instant::now();
    let arb_process = "cubic_exact_uv_rewrite";
    let arb_integrand = "amplitude";
    let arb_point = [1.3, -1.9, 3.1, 4.3, -2.9, 6.1];
    let arb_points = Array2::from_shape_vec((1, arb_point.len()), arb_point.to_vec())?;
    let mut arb_values = Vec::new();
    for (route, cli) in &mut states {
        if !matches!(*route, "explicit_local_3d" | "projected_local_4d") {
            continue;
        }
        let process_id = cli
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified(arb_process.to_string())))?;
        let graph_name = cli
            .state
            .process_list
            .get_integrand(process_id, arb_integrand)?
            .require_generated()?
            .graph_name_by_id(0)
            .ok_or_else(|| {
                eyre::eyre!("Generated {arb_integrand} integrand should have graph id 0")
            })?
            .to_string();
        let result = evaluate_sample_precise(
            &mut cli.state,
            &EvaluateSamplesPrecise {
                process_id: Some(process_id),
                integrand_name: Some(arb_integrand.to_string()),
                use_arb_prec: true,
                minimal_output: true,
                return_generated_events: Some(false),
                momentum_space: true,
                points: arb_points.view(),
                integrator_weights: None,
                discrete_dims: None,
                graph_names: Some(vec![Some(graph_name)]),
                orientations: Some(vec![None]),
            },
        )?;
        let value = match result.sample.evaluation {
            PreciseEvaluationResultOutput::Arb(result) => result.integrand_result,
            evaluation => {
                return Err(eyre::eyre!(
                    "{route} precise local-UV comparison requested Arb but got {:?}",
                    evaluation.precision()
                ));
            }
        };
        arb_values.push((*route, value));
    }
    let (reference_route, reference) = arb_values
        .iter()
        .find(|(route, _)| *route == "explicit_local_3d")
        .expect("the precise scalar local-UV comparison must include its local-3D reference");
    let (projected_route, projected) = arb_values
        .iter()
        .find(|(route, _)| *route == "projected_local_4d")
        .expect("the precise scalar local-UV comparison must include its local-4D route");
    assert!(
        [reference, projected].into_iter().all(|value| {
            symbolica::domains::float::SingleFloat::is_finite(&value.re)
                && symbolica::domains::float::SingleFloat::is_finite(&value.im)
        }),
        "precise scalar local-UV route values must be finite: {reference:e}, {projected:e}",
    );
    let distance = (projected.clone() - reference.clone()).norm().re;
    let projected_norm = projected.norm().re;
    let reference_norm = reference.norm().re;
    let scale = if projected_norm > reference_norm {
        projected_norm
    } else {
        reference_norm
    };
    // Retaining one eighth of the requested precision accepts about 125 matching bits
    // at ArbPrec's 1000-bit precision, allowing substantial route-dependent bit loss.
    let tolerance = F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt();
    let relative_distance = if scale.is_zero() {
        distance.clone()
    } else {
        distance.clone() / scale
    };
    if relative_distance > tolerance.clone() {
        failures.push(format!(
            "{arb_process} differs between precise {projected_route} and {reference_route}: {projected_route}={projected:e}, {reference_route}={reference:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
        ));
    }
    println!(
        "scalar local-UV Arb-to-Arb route comparison: {:?}",
        arb_started.elapsed()
    );

    for ((route, cli), evaluation_elapsed) in states.iter().zip(evaluation_times) {
        println!("scalar local-UV route {route}: evaluation {evaluation_elapsed:?}");
        clean_test(&cli.cli_settings.state.folder);
    }
    println!(
        "scalar local-UV three-route acceptance: total {:?}",
        total_started.elapsed()
    );
    assert!(
        failures.is_empty(),
        "scalar local-UV route mismatches:\n{}",
        failures.join("\n")
    );
    Ok(())
}

fn scalar_spectacles_integrated_route_comparison(use_arb_prec: bool) -> Result<()> {
    let precision = if use_arb_prec { "arb" } else { "f64" };
    let routes = [
        ("localized_local_3d", false, false),
        ("explicit_local_3d", true, false),
        ("projected_local_4d", true, true),
    ];
    let total_started = Instant::now();
    let mut states = Vec::new();

    for (route, explicit_orientation_sum_only, project_local_4d) in routes {
        let test_root = get_tests_workspace_path().join(format!(
            "scalar_spectacles_integrated_matches_across_local_uv_routes_{precision}_{route}"
        ));
        let mut cli = get_test_cli(
            Some("uv/scalar_spectacles_self_energy.toml".into()),
            test_root,
            Some(format!("{precision}_{route}")),
            true,
        )?;
        {
            let generation = &mut cli.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = explicit_orientation_sum_only;
            generation.uv.generate_integrated = true;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d;
        }
        let generation_started = Instant::now();
        cli.run_command("run generate")?;
        cli.run_command(
            "import graphs ./tests/resources/graphs/uv_tests/scalar_spectacles_overall_log.dot -p spectacles_overall_log -i amplitude",
        )?;
        cli.run_command("generate existing -p spectacles_overall_log -i amplitude")?;
        println!(
            "integrated spectacles route {route}: generation {:?}",
            generation_started.elapsed()
        );
        states.push((route, cli));
    }

    let evaluation_started = Instant::now();
    let cases: [(&str, &str, &[f64]); 2] = [
        (
            "spectacles",
            "scalar_spectacles",
            &[0.11, -0.29, 0.37, 0.59, -0.43, 0.71],
        ),
        (
            "spectacles_overall_log",
            "amplitude",
            &[0.11, -0.29, 0.37, 0.59, -0.43, 0.71, -0.17, 0.23, 0.41],
        ),
    ];
    macro_rules! compare_routes {
        ($variant:ident, $tolerance:expr) => {{
            let tolerance = $tolerance;
            for (process, integrand, point) in cases {
                let points = Array2::from_shape_vec((1, point.len()), point.to_vec())?;
                let mut values = Vec::new();
                for (route, cli) in &mut states {
                    let process_id = cli.state.resolve_process_ref(Some(
                        &ProcessRef::Unqualified(process.to_string()),
                    ))?;
                    let graph_name = cli
                        .state
                        .process_list
                        .get_integrand(process_id, integrand)?
                        .require_generated()?
                        .graph_name_by_id(0)
                        .ok_or_else(|| {
                            eyre::eyre!("Generated {integrand} should have graph id 0")
                        })?
                        .to_string();
                    let result = evaluate_sample_precise(
                        &mut cli.state,
                        &EvaluateSamplesPrecise {
                            process_id: Some(process_id),
                            integrand_name: Some(integrand.to_string()),
                            use_arb_prec,
                            minimal_output: true,
                            return_generated_events: Some(false),
                            momentum_space: true,
                            points: points.view(),
                            integrator_weights: None,
                            discrete_dims: None,
                            graph_names: Some(vec![Some(graph_name)]),
                            orientations: Some(vec![None]),
                        },
                    )?;
                    let value = match result.sample.evaluation {
                        PreciseEvaluationResultOutput::$variant(result) => {
                            result.integrand_result
                        }
                        evaluation => {
                            return Err(eyre::eyre!(
                                "{route} integrated {process} comparison requested {precision} but got {:?}",
                                evaluation.precision()
                            ));
                        }
                    };
                    values.push((*route, value));
                }

                let (reference_route, reference) = values
                    .iter()
                    .find(|(route, _)| *route == "explicit_local_3d")
                    .expect("the integrated spectacles matrix must include explicit local 3D");
                assert!(
                    !reference.norm().re.is_zero(),
                    "integrated {process} is zero"
                );
                for (route, value) in &values {
                    if route == reference_route {
                        continue;
                    }
                    let distance = (value.clone() - reference.clone()).norm().re;
                    let value_norm = value.norm().re;
                    let reference_norm = reference.norm().re;
                    let scale = if value_norm > reference_norm {
                        value_norm
                    } else {
                        reference_norm.clone()
                    };
                    let relative_distance = distance / scale;
                    assert!(
                        relative_distance <= tolerance,
                        "integrated {process} differs between {route} and {reference_route}: {route}={value:e}, {reference_route}={reference:e}, relative delta={relative_distance:e}, tolerance={tolerance:e}"
                    );
                }
            }
        }};
    }
    if use_arb_prec {
        compare_routes!(Arb, F(ArbPrec::default().epsilon()).sqrt().sqrt().sqrt());
    } else {
        compare_routes!(Double, F(1.0e-14));
    }
    println!(
        "integrated spectacles three-route {precision} comparison: evaluation {:?}, total {:?}",
        evaluation_started.elapsed(),
        total_started.elapsed()
    );
    for (_, cli) in states {
        clean_test(&cli.cli_settings.state.folder);
    }
    Ok(())
}

#[test]
#[serial_test::serial]
fn scalar_spectacles_integrated_matches_across_local_uv_routes() -> Result<()> {
    scalar_spectacles_integrated_route_comparison(false)
}

#[test]
#[serial_test::serial]
#[ignore = "1000-bit Arb exposes a non-scaling roughly 4e-15 localized/explicit route delta"]
fn scalar_spectacles_integrated_arb_route_precision_diagnostic() -> Result<()> {
    scalar_spectacles_integrated_route_comparison(true)
}

#[test]
fn scalar_spectacles_integrated_uv_factorizes_over_bridge() -> Result<()> {
    let mut cli = get_test_cli(
        Some("uv/scalar_spectacles_self_energy.toml".into()),
        get_tests_workspace_path().join("scalar_spectacles_self_energy"),
        None,
        true,
    )?;
    cli.run_command("run generate")?;
    // The normalized localization density integrates to one at any positive scale.
    // Match its radial support to mu_r instead of the default scale 1000.
    cli.run_command("set process kv general.renormalization_localization_scale=3.0")?;
    let bubble_accuracy_floor =
        required_stability_accuracy_floor(&mut cli, "bubble", "scalar_bubble")?;
    let spectacles_accuracy_floor =
        required_stability_accuracy_floor(&mut cli, "spectacles", "scalar_spectacles")?;

    // The massless B-C propagator contributes i/s and the graph numerator -i/2,
    // so the spectacles graph is the product of both bubbles divided by 2s.
    let external_s = 1.0_f64;
    let bridge_factor = 1.0 / (2.0 * external_s);
    let bubble_left = evaluate_momentum_sample(
        &mut cli,
        "bubble",
        "scalar_bubble",
        &[0.11, -0.07, 0.19],
        None,
        bubble_accuracy_floor,
    )?;
    let bubble_right = evaluate_momentum_sample(
        &mut cli,
        "bubble",
        "scalar_bubble",
        &[0.23, 0.31, -0.17],
        None,
        bubble_accuracy_floor,
    )?;
    let spectacles_point = evaluate_momentum_sample(
        &mut cli,
        "spectacles",
        "scalar_spectacles",
        &[0.11, -0.07, 0.19, 0.23, 0.31, -0.17],
        None,
        spectacles_accuracy_floor,
    )?;
    let expected_point = Complex::new(
        bridge_factor
            * (bubble_left.value.re * bubble_right.value.re
                - bubble_left.value.im * bubble_right.value.im),
        bridge_factor
            * (bubble_left.value.re * bubble_right.value.im
                + bubble_left.value.im * bubble_right.value.re),
    );
    let point_delta = (spectacles_point.value.re - expected_point.re)
        .hypot(spectacles_point.value.im - expected_point.im);
    let point_scale = spectacles_point
        .value
        .re
        .hypot(spectacles_point.value.im)
        .max(expected_point.re.hypot(expected_point.im))
        .max(f64::MIN_POSITIVE);
    let point_accuracy = bubble_left
        .relative_accuracy
        .max(bubble_right.relative_accuracy)
        .max(spectacles_point.relative_accuracy)
        .max(f64::EPSILON);
    assert!(
        point_delta / point_scale <= INSPECT_DEPENDENCE_ACCURACY_FACTOR * point_accuracy,
        "expected pointwise spectacles = bubble_left*bubble_right/(2s); relative delta={:.3e}, accuracy={point_accuracy:.3e}",
        point_delta / point_scale,
    );

    cli.run_command("set process kv integrator.integrated_phase=imag")?;
    set_fast_deterministic_integrator(
        &mut cli,
        IntegratedUvIntegratorSettings {
            target_relative_accuracy: 1.0e-12,
            n_start: 50_000,
            n_increase: 50_000,
            n_max: 100_000,
            n_cores: 10,
        },
    )?;
    cli.run_command(
        "set process kv integrator.min_samples_for_update=50000 \
         integrator.discrete_dim_learning_rate=0.0 \
         integrator.continuous_dim_learning_rate=0.1",
    )?;
    // The bubble-derived target and spectacles estimate need independent errors.
    cli.run_command("set process kv integrator.seed=7331")?;

    let bubble = Integrate {
        process: vec![ProcessRef::Unqualified("bubble".to_string())],
        integrand_name: vec!["scalar_bubble".to_string()],
        workspace_path: Some(
            get_tests_workspace_path()
                .join("scalar_spectacles_self_energy/integration_workspace_bubble"),
        ),
        n_cores: Some(10),
        restart: true,
        renderer: gammaloop_api::commands::integrate::RendererOption::Tabled,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?
    .single_slot()
    .ok_or_else(|| eyre::eyre!("single bubble slot should exist"))?
    .integral
    .clone();

    cli.run_command("set process kv integrator.integrated_phase=real")?;
    set_fast_deterministic_integrator(
        &mut cli,
        IntegratedUvIntegratorSettings {
            target_relative_accuracy: 1.0e-12,
            n_start: 50_000,
            n_increase: 50_000,
            n_max: 100_000,
            n_cores: 10,
        },
    )?;
    cli.run_command("set process kv integrator.seed=1337")?;
    let spectacles = Integrate {
        process: vec![ProcessRef::Unqualified("spectacles".to_string())],
        integrand_name: vec!["scalar_spectacles".to_string()],
        workspace_path: Some(
            get_tests_workspace_path()
                .join("scalar_spectacles_self_energy/integration_workspace_spectacles"),
        ),
        n_cores: Some(10),
        restart: true,
        renderer: gammaloop_api::commands::integrate::RendererOption::Tabled,
        ..Default::default()
    }
    .run(&mut cli.state, &cli.cli_settings)?
    .single_slot()
    .ok_or_else(|| eyre::eyre!("single spectacles slot should exist"))?
    .integral
    .clone();
    let mut expected = squared_integral_estimate(&bubble);
    expected.result *= F(bridge_factor);
    expected.error *= F(bridge_factor.abs());

    println!("bubble: {bubble:#}");
    println!("spectacles: {spectacles:#}");
    println!("bubble squared over 2s: {expected:#}");

    for (name, estimate) in [
        ("bubble", &bubble),
        ("spectacles", &spectacles),
        ("bubble squared over 2s", &expected),
    ] {
        let relative_error = integral_relative_error(estimate);
        assert!(
            relative_error < 0.10,
            "{name} relative error must be below 10%; got {:.2}%",
            100.0 * relative_error,
        );
    }
    let delta_sigma = integral_estimate_change_sigma(&spectacles, &expected);
    assert!(
        delta_sigma <= 1.0,
        "expected spectacles self-energy = bubble^2/(2s) within 1sigma; delta={delta_sigma}sigma"
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
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
        shifted_renormalization_localization_scale: 1.0,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: true,
        targets: IntegratedUvTargets {
            integrated: Some(Complex::new(F(0.0), F(-0.01451155120018305))),
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
        shifted_renormalization_localization_scale: 1.0,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: true,
        targets: IntegratedUvTargets {
            integrated: Some(Complex::new(F(0.0), F(-0.003628430563793077))),
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
        shifted_renormalization_localization_scale: 1.0,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: true,
        targets: IntegratedUvTargets {
            integrated: Some(Complex::new(F(0.0), F(-0.01234018404957018))),
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
        shifted_renormalization_localization_scale: 1.0,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: false,
        targets: IntegratedUvTargets {
            integrated: Some(Complex::new(F(-13169.321086130927), F(-56907.637709852635))),
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
        shifted_renormalization_localization_scale: 1.0,
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
fn dotted_sunrise_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/dotted_sunrise",
        test_name: "dotted_sunrise",
        process: "dotted_sunrise",
        integrand_name: "dotted_sunrise",
        integrator: SUNRISE_INTEGRATED_UV_INTEGRATOR,
        original_m_uv: 0.2,
        shifted_m_uv: 0.5,
        original_renormalization_localization_scale: 0.2,
        shifted_renormalization_localization_scale: 1.0,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: false,
        targets: IntegratedUvTargets::default(),
        integrated_ct_relative_error_limit: None,
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
        integrator: DOTTED_DT_INTEGRATED_UV_INTEGRATOR,
        original_m_uv: 0.2,
        shifted_m_uv: 0.5,
        original_renormalization_localization_scale: 0.2,
        shifted_renormalization_localization_scale: 1.0,
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
        integrator: EPEM_A_BBX_INTEGRATED_UV_INTEGRATOR,
        original_m_uv: 0.2,
        shifted_m_uv: 0.5,
        original_renormalization_localization_scale: 0.2,
        shifted_renormalization_localization_scale: 1.0,
        original_mu_r: 0.2,
        shifted_mu_r: 0.8,
        skip_uv_profile: false,
        targets: IntegratedUvTargets::default(),
        integrated_ct_relative_error_limit: None,
        check_mu_r_dependence: true,
    });
}

#[test]
fn helicity_amplitude_norm_stability_uses_both_components() -> Result<()> {
    let mut cli = get_test_cli(
        Some("uv/epem_a_bbx_amp.toml".into()),
        get_tests_workspace_path().join("helicity_amplitude_norm_stability"),
        None,
        true,
    )?;
    cli.run_command("run generate")?;
    cli.run_command("set process kv general.m_uv=0.2 general.mu_r=0.2 general.renormalization_localization_scale=0.2")?;
    let points = Array2::from_shape_vec((1, 3), vec![0.31, 0.47, 0.63])?;
    let discrete_dims = Array2::from_shape_vec((1, 2), vec![0, 0])?;
    for multichanneling in [false, true] {
        cli.run_command(&format!(
            "set process kv sampling.sampling_multichanneling={multichanneling}"
        ))?;
        let mut primary_value: Option<Complex<F<f64>>> = None;
        for phase in ["real", "both", "imag"] {
            cli.run_command(&format!(
                "set process kv integrator.integrated_phase={phase}"
            ))?;
            let result = evaluate_sample(
                &mut cli.state,
                &EvaluateSamples {
                    process_id: Some(0),
                    integrand_name: Some("epem_a_bbx".into()),
                    use_arb_prec: false,
                    minimal_output: false,
                    return_generated_events: Some(false),
                    momentum_space: false,
                    points: points.view(),
                    integrator_weights: None,
                    discrete_dims: Some(discrete_dims.view()),
                    graph_names: None,
                    orientations: None,
                },
            )?;
            let evaluation = result.sample.evaluation;
            let levels = &evaluation.evaluation_metadata.unwrap().stability_results;
            // The returned primary point is independent of the component selected
            // for integration. The complex norm is invariant under helicity phases.
            if let Some(primary) = &primary_value {
                assert!((evaluation.integrand_result - *primary).norm_squared().0 < 1e-30);
            } else {
                primary_value = Some(evaluation.integrand_result);
            }
            assert_eq!(
                levels.len(),
                1,
                "the complex norm needs no precision rescue"
            );
            assert_eq!(
                levels[0].precision,
                gammalooprs::settings::runtime::Precision::Double
            );
            assert!(matches!(
                levels[0].status,
                gammalooprs::integrands::evaluation::StabilityStatus::Stable(2)
            ));
            assert!(levels[0].estimated_relative_accuracy.unwrap().0 < 1e-12);
        }
    }
    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

const AA_AA_2L_UV_RICH_INSPECT: GraphUvRichInspectCase = GraphUvRichInspectCase {
    name: "aa_aa 2L",
    test_prefix: "aa_aa_2l",
    process_prefix: "aa_aa_2l_",
    integrand_name: "2L",
    run_card: "uv/aa_aa_2l_rich_inspect.toml",
    graph_block_prefix: "generate_",
    target_file: "aa_aa_2l_uv_inspect_event_targets.json",
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

    fn run_graph_local_uv_route_comparison(
        &self,
        graph: &str,
        source: &str,
        integrand_name: &str,
        point: &[f64],
    ) -> Result<()> {
        use gammalooprs::{integrands::process::ProcessIntegrand, uv::UltravioletGraph};
        use linnet::half_edge::subgraph::SubSetLike;
        use std::collections::BTreeSet;

        let points = [
            point.to_vec(),
            point.iter().map(|value| 100.0 * value).collect(),
            [0.23, 0.31, -0.17, 0.41, -0.29, 0.13, -0.37, 0.47, 0.59][..point.len()].to_vec(),
        ];
        let mut results = Vec::new();
        for (route, explicit, direct) in [
            ("localized_3d", false, false),
            ("erased_3d", true, false),
            ("direct_4d", true, true),
        ] {
            let test_name = format!("aa_aa_{graph}_local_uv_routes_{route}");
            // The run card prepares the physical model and kinematics. Import
            // directly: its existing generate blocks enable integrated CTs.
            let mut cli = get_test_cli(
                Some(self.run_card.into()),
                get_tests_workspace_path().join(&test_name),
                Some(test_name),
                true,
            )?;
            let generation = &mut cli.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = explicit;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = direct;
            generation.uv.generate_integrated = false;
            generation.uv.orchestrator = UVOrchestrator::Compare;
            generation.threshold_subtraction.enable_thresholds = false;
            generation
                .threshold_subtraction
                .check_esurface_at_generation = false;
            let process = format!("aa_aa_{integrand_name}_{graph}_{route}");
            cli.run_command(&format!(
                "import graphs ./{source} -p {process} -i {integrand_name} -o"
            ))?;
            let started = Instant::now();
            cli.run_command(&format!(
                "generate existing -p {process} -i {integrand_name}"
            ))?;
            println!(
                "physical {graph} {route}: generation {:?}",
                started.elapsed()
            );

            let process_id = cli
                .state
                .resolve_process_ref(Some(&ProcessRef::Unqualified(process.clone())))?;
            let integrand = cli
                .state
                .process_list
                .get_integrand(process_id, integrand_name)?
                .require_generated()?;
            assert_eq!(integrand.get_n_dim(), point.len());
            if graph == "GL262" {
                let ProcessIntegrand::Amplitude(amplitude) = integrand else {
                    panic!("GL262 must be an amplitude")
                };
                let graph = &amplitude.data.graph_terms[0].graph;
                let regions = graph
                    .spinneys(&graph.full_filter())
                    .into_iter()
                    .filter(|region| !region.is_empty())
                    .map(|region| {
                        let mut edges = graph
                            .iter_edges_of(&region)
                            .map(|(_, edge, _)| usize::from(edge))
                            .collect::<Vec<_>>();
                        edges.sort_unstable();
                        (edges, graph.n_loops(&region), graph.local_dod(&region))
                    })
                    .collect::<BTreeSet<_>>();
                // The connected gluon-self-energy / top-self-energy / full-graph
                // chain gives eight conventional forests, including empty.
                // Production also admits the disconnected hexagon-plus-bubble
                // union by its aggregate DOD (-2 + 2 = 0), adding two nodes to
                // the unfolded forest. Keep this existing contribution covered.
                assert_eq!(
                    regions,
                    BTreeSet::from([
                        (vec![9, 10], 1, 2),
                        (vec![4, 5, 7, 9, 10], 2, 1),
                        (vec![4, 6, 8, 9, 10, 11, 12, 13], 2, 0),
                        ((4..14).collect(), 3, 0),
                    ]),
                    "GL262 must retain the physical DOD-2 gluon self-energy chain",
                );
            }

            let evaluations = points
                .iter()
                .map(|point| {
                    evaluate_momentum_sample(
                        &mut cli,
                        &process,
                        integrand_name,
                        point,
                        None,
                        1.0e-12,
                    )
                })
                .collect::<Result<Vec<_>>>()?;
            for evaluation in &evaluations {
                assert!(
                    evaluation.value.re.is_finite()
                        && evaluation.value.im.is_finite()
                        && evaluation.value.re.hypot(evaluation.value.im) > 0.0,
                    "{graph} {route} must evaluate to a finite, nonzero value",
                );
            }
            results.push((route, evaluations));
            clean_test(&cli.cli_settings.state.folder);
        }

        let reference = &results[0].1;
        for (route, evaluations) in &results[1..] {
            for ((point, expected), actual) in points.iter().zip(reference).zip(evaluations) {
                let delta = (actual.value.re - expected.value.re)
                    .hypot(actual.value.im - expected.value.im);
                let scale = actual
                    .value
                    .re
                    .hypot(actual.value.im)
                    .max(expected.value.re.hypot(expected.value.im));
                let tolerance = 1.0e-9;
                assert!(
                    delta / scale <= tolerance,
                    "{graph} {route} differs from localized 3D at {point:?}: actual={:?}, reference={:?}, relative delta={:.3e}, tolerance={tolerance:.3e}, reported accuracy actual={:.3e}, reference={:.3e}",
                    actual.value,
                    expected.value,
                    delta / scale,
                    actual.relative_accuracy,
                    expected.relative_accuracy,
                );
            }
        }
        Ok(())
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
        AdditionalWeightKey::AmplitudeThresholdCountertermVariant {
            variant_id,
            esurface_id,
            overlap_group,
        } => format!("ct_variant_{variant_id}_{esurface_id}_{overlap_group}"),
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

fn local_ct_cli(
    run_card: &str,
    test_name: &str,
    unsubtracted: bool,
    project_local_4d: bool,
    explicit_orientation_sum_only: bool,
) -> Result<CLIState> {
    let test_name = format!(
        "{test_name}_explicit_{explicit_orientation_sum_only}_project_local_4d_{project_local_4d}"
    );
    let mut cli = get_test_cli(
        Some(run_card.into()),
        get_tests_workspace_path().join(&test_name),
        Some(test_name),
        true,
    )?;
    let generation = &mut cli.cli_settings.global.generation;
    generation.explicit_orientation_sum_only = explicit_orientation_sum_only;
    if explicit_orientation_sum_only {
        generation.orientation_pattern = Default::default();
    }
    let uv = &mut generation.uv;
    uv.final_integrand = FinalIntegrandDimension::ThreeD;
    uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d;
    assert!(
        !uv.generate_integrated,
        "local profile fixtures intentionally isolate unintegrated counterterms"
    );
    if unsubtracted {
        let prescription = &mut uv.renormalization_prescription;
        prescription.log_divergent = ApproximationType::Unsubtracted;
        prescription.massive_power_divergent = ApproximationType::Unsubtracted;
        prescription.massless_power_divergent = ApproximationType::Unsubtracted;
        prescription.overrides.clear();
        // The legacy wood does not represent an entirely unsubtracted graph,
        // so Compare is meaningful only for the selected local counterterm.
        // Keep the bare side of each non-vacuity check on the reference path.
        if uv.orchestrator == UVOrchestrator::Compare {
            uv.orchestrator = UVOrchestrator::HedgePoset;
        }
    }
    cli.run_command("run generate")?;
    Ok(cli)
}

#[test]
#[serial_test::serial]
fn local_ir_integrated_generation_succeeds_across_local_uv_routes() -> Result<()> {
    for (name, external_pdg_sets, internal_pdgs, component_counts, momenta, require_nonzero) in [
        (
            "local_ir_quark_self_energy",
            &[[-1, 1]][..],
            &[1, 21][..],
            &[1][..],
            &[0.11, -0.29, 0.37][..],
            false,
        ),
        (
            "local_ir_gluon_self_energy",
            &[[-21, 21]][..],
            &[1][..],
            &[1][..],
            &[0.11, -0.29, 0.37][..],
            false,
        ),
        (
            "local_ir_dod2_scalar_spectacles",
            &[[-1000, 1001], [-1001, 1000]][..],
            &[1001][..],
            &[1, 1][..],
            &[0.11, -0.29, 0.37, 0.59, -0.43, 0.71][..],
            true,
        ),
        (
            "local_nested_soft_top_self_energy",
            &[[-6, 6]][..],
            &[6, 21][..],
            &[2][..],
            &[0.11, -0.29, 0.37, 0.59, -0.43, 0.71][..],
            true,
        ),
    ] {
        let mut values = Vec::new();
        let routes = [
            (false, false, true),
            (true, false, true),
            (true, true, true),
        ]
        .into_iter()
        .chain((name == "local_nested_soft_top_self_energy").then_some((true, false, false)));
        for (
            explicit_orientation_sum_only,
            project_local_4d,
            project_integrated_uv_cts_onto_tensor_integrals,
        ) in routes
        {
            let test_name = format!(
                "{name}_integrated_explicit_{explicit_orientation_sum_only}_project_local_4d_{project_local_4d}_project_integrated_{project_integrated_uv_cts_onto_tensor_integrals}"
            );
            let mut cli = get_test_cli(
                Some(format!("{name}.toml").into()),
                get_tests_workspace_path().join(&test_name),
                Some(test_name),
                true,
            )?;
            let generation = &mut cli.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = explicit_orientation_sum_only;
            generation.orientation_pattern = Default::default();
            generation.uv.final_integrand = FinalIntegrandDimension::ThreeD;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d;
            generation.uv.generate_integrated = true;
            generation
                .uv
                .project_integrated_uv_cts_onto_tensor_integrals =
                project_integrated_uv_cts_onto_tensor_integrals;
            cli.run_command("run generate")?;
            assert_generated_local_integrand(&cli, name, name, momenta.len())?;
            assert_eq!(
                selected_local_identifiers(&cli, name, name, ApproximationType::IR)?.len(),
                component_counts.iter().sum::<usize>(),
                "{name} must retain every expected integrated IR component",
            );
            for (external_pdgs, minimum_count) in external_pdg_sets.iter().zip(component_counts) {
                assert_selected_local_identifier(
                    &cli,
                    name,
                    name,
                    ApproximationType::IR,
                    external_pdgs.iter().copied(),
                    internal_pdgs.iter().copied(),
                    *minimum_count,
                )?;
            }

            // Force arbitrary precision for route equality; serialized values
            // are f64, so their rounding sets the comparison's accuracy floor.
            let process_id = cli
                .state
                .resolve_process_ref(Some(&ProcessRef::Unqualified(name.to_string())))?;
            let points = Array2::from_shape_vec((1, momenta.len()), momenta.to_vec())?;
            let result = evaluate_sample(
                &mut cli.state,
                &EvaluateSamples {
                    process_id: Some(process_id),
                    integrand_name: Some(name.to_string()),
                    use_arb_prec: true,
                    minimal_output: false,
                    return_generated_events: Some(false),
                    momentum_space: true,
                    points: points.view(),
                    integrator_weights: None,
                    discrete_dims: None,
                    graph_names: None,
                    orientations: Some(vec![None]),
                },
            )?;
            let evaluation = result.sample.evaluation;
            let metadata = evaluation.evaluation_metadata.as_ref().ok_or_else(|| {
                eyre::eyre!("integrated route comparison must report its numerical accuracy")
            })?;
            let value = InspectEvaluation {
                value: evaluation.integrand_result.map(|entry| entry.0),
                relative_accuracy: stability_relative_accuracy(metadata, f64::EPSILON),
            };
            assert!(
                value.value.re.is_finite() && value.value.im.is_finite(),
                "{name} integrated local IR generation must evaluate finitely for explicit={explicit_orientation_sum_only}, project_local_4d={project_local_4d}, project_integrated={project_integrated_uv_cts_onto_tensor_integrals}: {:?}",
                value.value,
            );
            assert!(
                !require_nonzero || value.value.re.hypot(value.value.im) > 0.0,
                "{name} integrated local IR generation must be nonzero for explicit={explicit_orientation_sum_only}, project_local_4d={project_local_4d}, project_integrated={project_integrated_uv_cts_onto_tensor_integrals}",
            );
            values.push((
                (
                    explicit_orientation_sum_only,
                    project_local_4d,
                    project_integrated_uv_cts_onto_tensor_integrals,
                ),
                value,
            ));
            clean_test(&cli.cli_settings.state.folder);
        }
        // External on-shell spinors can annihilate the quark self-energy's
        // finite addback. Its uncontracted nonzero atom has a separate unit oracle.
        // Compare complete physical sums with integration enabled in all routes.
        let reference = &values[0].1;
        for (
            (
                explicit_orientation_sum_only,
                project_local_4d,
                project_integrated_uv_cts_onto_tensor_integrals,
            ),
            value,
        ) in &values[1..]
        {
            let delta =
                (value.value.re - reference.value.re).hypot(value.value.im - reference.value.im);
            let scale = reference
                .value
                .re
                .hypot(reference.value.im)
                .max(value.value.re.hypot(value.value.im))
                .max(f64::MIN_POSITIVE);
            let relative_accuracy = reference
                .relative_accuracy
                .max(value.relative_accuracy)
                .max(f64::EPSILON);
            assert!(
                relative_accuracy.is_finite()
                    && INSPECT_DEPENDENCE_ACCURACY_FACTOR * relative_accuracy < 1.0
                    && delta <= INSPECT_DEPENDENCE_ACCURACY_FACTOR * relative_accuracy * scale,
                "{name} integrated local IR routes disagree for explicit={explicit_orientation_sum_only}, project_local_4d={project_local_4d}, project_integrated={project_integrated_uv_cts_onto_tensor_integrals}: localized direct={:?}, value={:?}, delta={delta:.3e}, scale={scale:.3e}, accuracy={relative_accuracy:.3e}",
                reference.value,
                value.value,
            );
        }
    }
    Ok(())
}

#[test]
#[serial_test::serial]
fn genuine_vacuum_tadpole_matches_all_local_uv_routes_and_vakint_inputs() -> Result<()> {
    const VACUUM_TADPOLE: &str = r#"
        digraph vacuum_tadpole {
            edge [num=1 mass=2]
            node [num=1]
            a -> a [id=0 lmb_id=0]
        }
    "#;
    const VACUUM_RUNTIME: &str = r#"
        [general]
        evaluator_method = "SingleParametric"
        m_uv = 1000.0
        renormalization_localization_scale = 1000.0
        mu_r = 1000.0

        [kinematics]
        e_cm = 64.0

        [kinematics.externals]
        type = "constant"

        [kinematics.externals.data]
        momenta = []
        helicities = []

        [sampling]
        graphs = "summed"
        orientations = "summed"

        [subtraction]
        disable_threshold_subtraction = true
    "#;
    let mut values = Vec::new();
    for project_integrated_uv_cts_onto_tensor_integrals in [false, true] {
        for (route, explicit_orientation_sum_only, project_local_4d) in [
            ("localized_direct_3d", false, false),
            ("explicit_direct_3d", true, false),
            ("projected_local_4d", true, true),
        ] {
            let test_name = format!(
                "genuine_vacuum_tadpole_{route}_project_integrated_{project_integrated_uv_cts_onto_tensor_integrals}"
            );
            let mut cli = get_test_cli(
                None,
                get_tests_workspace_path().join(&test_name),
                Some(test_name),
                true,
            )?;
            let generation = &mut cli.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = explicit_orientation_sum_only;
            generation.orientation_pattern = Default::default();
            generation
                .tropical_subgraph_table
                .disable_tropical_generation = true;
            generation.evaluator.iterative_orientation_optimization = false;
            generation.threshold_subtraction.enable_thresholds = false;
            generation.uv.generate_integrated = true;
            generation
                .uv
                .project_integrated_uv_cts_onto_tensor_integrals =
                project_integrated_uv_cts_onto_tensor_integrals;
            generation.uv.subtract_uv = true;
            generation.uv.final_integrand = FinalIntegrandDimension::ThreeD;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d;
            generation.uv.add_marker = false;
            generation.uv.inner_products = true;
            generation.uv.orchestrator = UVOrchestrator::HedgePoset;
            let prescription = &mut generation.uv.renormalization_prescription;
            prescription.log_divergent = ApproximationType::MUV;
            prescription.massive_power_divergent = ApproximationType::MUV;
            prescription.massless_power_divergent = ApproximationType::MUV;
            prescription.overrides.clear();

            cli.run_command("import model scalars-default.json")?;
            cli.run_command(&format!(
                r#"import graphs --inline-dot """{VACUUM_TADPOLE}""" -p vacuum -i tadpole"#
            ))?;
            cli.run_command(&format!("set default-runtime string '{VACUUM_RUNTIME}'"))?;
            cli.run_command("generate existing -p vacuum -i tadpole")?;
            assert_generated_local_integrand(&cli, "vacuum", "tadpole", 3)?;

            let spinneys = classified_local_spinneys(&cli, "vacuum", "tadpole")?;
            assert_eq!(
                spinneys.len(),
                2,
                "the vacuum tadpole must have its physical component and inert empty root: {spinneys:?}"
            );
            let spinney = spinneys
                .iter()
                .find(|spinney| spinney.edge_ids == [0])
                .expect("the vacuum tadpole must retain its physical component");
            assert!(spinney.identifier.external_pdg_set.is_empty());
            assert_eq!(spinney.edge_ids, [0]);
            assert_eq!(spinney.n_components, 1);
            assert_eq!(spinney.scheme, ApproximationType::MUV);
            assert_eq!(spinney.dod, 2);
            let root = spinneys
                .iter()
                .find(|spinney| spinney.edge_ids.is_empty())
                .expect("the vacuum tadpole must retain its inert empty root");
            assert!(root.identifier.external_pdg_set.is_empty());
            assert_eq!(root.identifier.internal_pdg_set, Some(BTreeSet::new()));
            assert_eq!(root.n_components, 0);
            assert_eq!(root.scheme, ApproximationType::MUV);
            assert_eq!(root.dod, 0);

            let value = evaluate_momentum_sample(
                &mut cli,
                "vacuum",
                "tadpole",
                &[0.11, -0.29, 0.37],
                None,
                f64::EPSILON,
            )?;
            assert!(
                value.value.re.is_finite() && value.value.im.is_finite(),
                "{route}, project_integrated={project_integrated_uv_cts_onto_tensor_integrals}: vacuum evaluation is not finite: {:?}",
                value.value,
            );
            assert!(
                value.value.re.hypot(value.value.im) > 0.0,
                "{route}, project_integrated={project_integrated_uv_cts_onto_tensor_integrals}: vacuum evaluation is zero",
            );
            values.push((
                (route, project_integrated_uv_cts_onto_tensor_integrals),
                value,
            ));
            clean_test(&cli.cli_settings.state.folder);
        }
    }

    let reference = &values
        .iter()
        .find(|((route, project_integrated), _)| {
            *route == "explicit_direct_3d" && *project_integrated
        })
        .expect("the vacuum matrix must contain its explicit direct reference")
        .1;
    assert!(
        reference.value.re.abs() <= 1.0e-12
            && (reference.value.im + 9.760529078735244e-4).abs() <= 1.0e-12,
        "the genuine-vacuum reference changed: {:?}",
        reference.value,
    );
    for ((route, project_integrated), value) in &values {
        let delta =
            (value.value.re - reference.value.re).hypot(value.value.im - reference.value.im);
        let scale = reference
            .value
            .re
            .hypot(reference.value.im)
            .max(value.value.re.hypot(value.value.im))
            .max(f64::MIN_POSITIVE);
        let relative_accuracy = reference
            .relative_accuracy
            .max(value.relative_accuracy)
            .max(f64::EPSILON);
        assert!(
            relative_accuracy.is_finite()
                && INSPECT_DEPENDENCE_ACCURACY_FACTOR * relative_accuracy < 1.0
                && delta <= INSPECT_DEPENDENCE_ACCURACY_FACTOR * relative_accuracy * scale,
            "{route}, project_integrated={project_integrated}: genuine-vacuum routes disagree: reference={:?}, value={:?}, delta={delta:.3e}, scale={scale:.3e}, accuracy={relative_accuracy:.3e}",
            reference.value,
            value.value,
        );
    }
    Ok(())
}

#[test]
#[serial_test::serial]
fn local_os_integrated_generation_is_rejected_by_run_generate() -> Result<()> {
    let mut cli = get_test_cli(
        Some(LOCAL_OS_CARD.into()),
        get_tests_workspace_path().join("local_os_integrated_generation_rejected"),
        Some("local_os_integrated_generation_rejected".to_string()),
        true,
    )?;
    cli.cli_settings.global.generation.uv.generate_integrated = true;

    let error = cli
        .run_command("run generate")
        .expect_err("the production generate command must reject integrated local OS CTs");
    let message = format!("{error:#}");
    assert!(
        message.contains("local on-shell/OS counterterms are local-only in phase 1")
            && message.contains("generate_integrated=false"),
        "the generation error lacks the phase-1 OS policy and corrective setting: {message}"
    );

    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(LOCAL_OS_NAME.to_string())))?;
    let resolved = cli
        .state
        .process_list
        .get_integrand(process_id, LOCAL_OS_NAME)?;
    assert!(
        resolved.integrand.is_none(),
        "a rejected production generation must not publish a partial integrand"
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial_test::serial]
fn local_ir_and_pole_part_wood_is_rejected_by_run_generate() -> Result<()> {
    const PROCESS: &str = "dgse";
    const INTEGRAND: &str = "dgse";
    let mut cli = get_test_cli(
        Some("dgse_local_ir.toml".into()),
        get_tests_workspace_path().join("local_ir_pole_part_wood_rejected"),
        Some("local_ir_pole_part_wood_rejected".to_string()),
        true,
    )?;
    cli.cli_settings
        .global
        .generation
        .uv
        .renormalization_prescription
        .overrides
        .push(CTRenormalizationRule::new(
            CTIdentifier::new(
                [-1, 1, 21].into_iter().collect(),
                Some([1, 21].into_iter().collect()),
            ),
            ApproximationType::PolePart,
        ));

    let error = cli
        .run_command("run generate")
        .expect_err("a local IR/PolePart wood must be rejected by production generation");
    let message = format!("{error:#}");
    assert!(
        message.contains("graph 'se' cut 0")
            && message.contains("local soft/IR and PolePart counterterms")
            && message.contains("cannot be combined in one UV wood")
            && message.contains("component"),
        "the production rejection lacks graph, cut, component, and scheme-policy context: {message}"
    );

    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(PROCESS.to_string())))?;
    let resolved = cli
        .state
        .process_list
        .get_integrand(process_id, INTEGRAND)?;
    assert!(
        resolved.integrand.is_none(),
        "a rejected IR/PolePart generation must not publish a partial integrand"
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial_test::serial]
fn local_ir_child_under_positive_degree_muv_parent_is_rejected_by_run_generate() -> Result<()> {
    const NAME: &str = "nested_positive_degree_hu";
    let mut cli = get_test_cli(
        Some("local_ir_child_muv_parent_rejected.toml".into()),
        get_tests_workspace_path().join("local_ir_child_muv_parent_rejected"),
        Some("local_ir_child_muv_parent_rejected".to_string()),
        true,
    )?;
    let uv = &cli.cli_settings.global.generation.uv;
    assert!(
        !uv.generate_integrated,
        "the policy fixture must be local-only"
    );
    assert_eq!(uv.orchestrator, UVOrchestrator::HedgePoset);
    assert_eq!(
        uv.renormalization_prescription.overrides,
        [CTRenormalizationRule::new(
            CTIdentifier::new(
                [-1001, 1001].into_iter().collect(),
                Some([1002].into_iter().collect()),
            ),
            ApproximationType::IR,
        )],
        "the run card must select only the scalar_2 child bubble as IR"
    );

    let error = cli.run_command("run generate").expect_err(
        "a positive-degree MUV parent over a positive-degree IR child must be rejected",
    );
    let message = format!("{error:#}");
    assert!(
        message.contains("graph 'nested_positive_degree_hu' cut 0")
            && message.contains("positive-degree MUV component 140 (d=1)")
            && message.contains("cannot contain local soft/IR component FU (d=1)")
            && message.contains("in phase 1")
            && message.contains("reduced-cograph subtraction is not local")
            && message.contains("assign IR to the containing component")
            && message.contains("deferred finite scheme-change/integrated policy"),
        "the production rejection lacks graph, cut, component, degree, or phase-policy context: {message}"
    );

    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(NAME.to_string())))?;
    let resolved = cli.state.process_list.get_integrand(process_id, NAME)?;
    assert!(
        resolved.integrand.is_none(),
        "a rejected positive-degree H/U generation must not publish a partial integrand"
    );

    clean_test(&cli.cli_settings.state.folder);
    Ok(())
}

#[derive(Debug)]
struct ClassifiedLocalSpinney {
    identifier: CTIdentifier,
    edge_ids: Vec<usize>,
    n_components: usize,
    scheme: ApproximationType,
    dod: i32,
}

fn classified_local_spinneys(
    cli: &CLIState,
    process_name: &str,
    integrand_name: &str,
) -> Result<Vec<ClassifiedLocalSpinney>> {
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(process_name.to_string())))?;
    let process = &cli.state.process_list.processes[process_id];
    let amplitudes = match &process.collection {
        ProcessCollection::Amplitudes(amplitudes) => amplitudes,
        ProcessCollection::CrossSections(_) => {
            return Err(eyre::eyre!(
                "local CT fixture '{process_name}' must be an amplitude"
            ));
        }
    };
    let amplitude = amplitudes
        .get(integrand_name)
        .ok_or_else(|| eyre::eyre!("missing local CT integrand '{integrand_name}'"))?;
    let settings = &cli.cli_settings.global.generation.uv;
    let mut classified = Vec::new();
    for amplitude_graph in &amplitude.graphs {
        let graph = &amplitude_graph.graph;
        classified.extend(
            graph
                .classified_spinneys(&graph.full_filter(), settings, &graph.loop_momentum_basis)?
                .into_iter()
                .map(|spinney| {
                    let mut edge_ids = graph
                        .iter_edges_of(spinney.filter())
                        .map(|(_, edge, _)| usize::from(edge))
                        .collect::<Vec<_>>();
                    edge_ids.sort_unstable();
                    ClassifiedLocalSpinney {
                        identifier: graph.ct_identifier(spinney.filter()),
                        edge_ids,
                        n_components: spinney.n_components(),
                        scheme: spinney.renormalization_scheme,
                        dod: spinney.dod,
                    }
                }),
        );
    }
    Ok(classified)
}

fn local_ct_canonical_routes(
    cli: &CLIState,
    process_name: &str,
    integrand_name: &str,
) -> Result<Vec<gammalooprs::graph::LoopMomentumBasis>> {
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(process_name.to_string())))?;
    let ProcessCollection::Amplitudes(amplitudes) =
        &cli.state.process_list.processes[process_id].collection
    else {
        return Err(eyre::eyre!(
            "local CT fixture '{process_name}' must be an amplitude"
        ));
    };
    Ok(amplitudes
        .get(integrand_name)
        .ok_or_else(|| eyre::eyre!("missing local CT integrand '{integrand_name}'"))?
        .graphs
        .iter()
        .map(|graph| graph.graph.loop_momentum_basis.clone())
        .collect())
}

fn assert_matched_local_routes(
    subtracted: &CLIState,
    bare: &CLIState,
    process_name: &str,
    integrand_name: &str,
) -> Result<()> {
    assert_eq!(
        local_ct_canonical_routes(subtracted, process_name, integrand_name)?,
        local_ct_canonical_routes(bare, process_name, integrand_name)?,
        "subtracted and bare {process_name}/{integrand_name} fixtures must use the same canonical routes"
    );
    Ok(())
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct UVForestTermIdentity {
    forest_index: usize,
    node_index: usize,
    node_key: String,
    term_index: usize,
    residue_index: gammalooprs::cff::CutCFFIndex,
}

impl From<&gammalooprs::uv::export::UVForestNodeTerm> for UVForestTermIdentity {
    fn from(term: &gammalooprs::uv::export::UVForestNodeTerm) -> Self {
        Self {
            forest_index: term.forest_index,
            node_index: term.node_index,
            node_key: term.node_key.clone(),
            term_index: term.term_index,
            residue_index: term.residue_index,
        }
    }
}

fn assert_unique_uv_forest_term_identities(
    export: &gammalooprs::uv::export::UVForestExport,
    context: &str,
) {
    assert_unique_uv_forest_term_identities_with_dimension(export, context, true);
}

fn assert_unique_uv_forest_term_identities_with_dimension(
    export: &gammalooprs::uv::export::UVForestExport,
    context: &str,
    is_three_dimensional: bool,
) {
    assert!(
        !export.node_terms.is_empty(),
        "{context} exported no final forest atoms"
    );
    let mut identities = BTreeSet::new();
    let mut combined_atoms = BTreeSet::new();
    for term in &export.node_terms {
        let identity = UVForestTermIdentity::from(term);
        assert!(
            identities.insert(identity.clone()),
            "{context} duplicated the full export identity {identity:?}"
        );
        assert!(
            combined_atoms.insert((
                identity.forest_index,
                identity.node_index,
                identity.node_key.clone(),
                identity.residue_index,
            )),
            "{context} emitted more than one final combined atom for forest {}, node {} ({}), residue {:?}",
            identity.forest_index,
            identity.node_index,
            identity.node_key,
            identity.residue_index,
        );
    }
    assert_eq!(
        identities.len(),
        export.node_terms.len(),
        "{context} must preserve one full export identity per final atom"
    );

    let mut canonical_component_routes = BTreeMap::<String, (Value, Value)>::new();
    let exported_nodes = export
        .node_terms
        .iter()
        .map(|term| (term.forest_index, term.node_key.as_str()))
        .collect::<BTreeSet<_>>();
    for term in &export.node_terms {
        let provenance = term
            .forest_provenance()
            .unwrap_or_else(|error| {
                panic!("{context} could not parse computed Hedge provenance: {error}")
            })
            .unwrap_or_else(|| panic!("{context} has no computed Hedge provenance"));
        assert_eq!(
            provenance.as_object().map(serde_json::Map::len),
            Some(1),
            "{context} forest provenance must contain only node-level data: {provenance}",
        );
        let local = &provenance["node"]["local"];
        let representation = local["representation"]
            .as_str()
            .unwrap_or_else(|| panic!("{context} has untagged local provenance: {local}"));
        assert_eq!(
            representation,
            if is_three_dimensional { "3d" } else { "4d" },
            "{context} provenance does not match its exported numerator"
        );
        match representation {
            "4d" => {
                for component in local["components"]
                    .as_array()
                    .unwrap_or_else(|| panic!("{context} has malformed component provenance"))
                {
                    let label = component["component"]
                        .as_str()
                        .unwrap_or_else(|| {
                            panic!("{context} has an unlabeled component: {component}")
                        })
                        .to_string();
                    let route = (
                        component["route_loop_edges"].clone(),
                        component["route_external_edges"].clone(),
                    );
                    if let Some(existing) =
                        canonical_component_routes.insert(label.clone(), route.clone())
                    {
                        assert_eq!(
                            existing, route,
                            "{context} exported incompatible canonical routes for component {label}",
                        );
                    }
                }
            }
            "3d" => {
                let paths = local["projection_paths"]
                    .as_array()
                    .unwrap_or_else(|| panic!("{context} has malformed 3D projection paths"));
                for path in paths {
                    let steps = path["steps"]
                        .as_array()
                        .unwrap_or_else(|| panic!("{context} has a malformed 3D path: {path}"));
                    assert!(
                        steps.windows(2).all(|pair| {
                            let child = pair[0]["topo_order"].as_u64().unwrap_or_else(|| {
                                panic!("{context} has no child topology order: {}", pair[0])
                            });
                            let parent = pair[1]["topo_order"].as_u64().unwrap_or_else(|| {
                                panic!("{context} has no parent topology order: {}", pair[1])
                            });
                            child <= parent
                        }),
                        "{context} must record each direct-3D path in child-to-parent replay order: {path}"
                    );
                    for step in steps {
                        let current = step["current_component"].as_str().unwrap_or_else(|| {
                            panic!("{context} has an unlabeled 3D projection step: {step}")
                        });
                        let rescaled = step["rescaled_subgraph"].as_str().unwrap_or_else(|| {
                            panic!("{context} has no rescaled component context: {step}")
                        });
                        let rescaling = step["rescaling"].as_str().unwrap_or_else(|| {
                            panic!("{context} has no 3D rescaling context: {step}")
                        });
                        let label = format!("{current}|{rescaled}|{rescaling}");
                        let route = (
                            step["route_loop_edges"].clone(),
                            step["route_external_edges"].clone(),
                        );
                        if let Some(existing) =
                            canonical_component_routes.insert(label.clone(), route.clone())
                        {
                            assert_eq!(
                                existing, route,
                                "{context} exported incompatible direct-3D routes for {label}",
                            );
                        }
                    }
                }
            }
            other => panic!("{context} has unsupported provenance representation {other:?}"),
        }
    }

    for term in &export.node_terms {
        let provenance = term
            .forest_provenance()
            .unwrap_or_else(|error| {
                panic!("{context} could not parse computed Hedge provenance: {error}")
            })
            .unwrap_or_else(|| panic!("{context} has no computed Hedge provenance"));
        let node = &provenance["node"];
        let parent_keys = node["parent_keys"]
            .as_array()
            .unwrap_or_else(|| panic!("{context} has malformed parent keys: {node}"));
        let local = &node["local"];
        let representation = local["representation"]
            .as_str()
            .unwrap_or_else(|| panic!("{context} has untagged local provenance: {local}"));
        for parent_key in parent_keys {
            let parent_key = parent_key
                .as_str()
                .unwrap_or_else(|| panic!("{context} has a non-string parent key: {node}"));
            assert!(
                exported_nodes.contains(&(term.forest_index, parent_key)),
                "{context} node {} refers to missing parent {parent_key:?} in forest {}",
                term.node_key,
                term.forest_index,
            );
        }
        if !is_three_dimensional {
            assert_eq!(
                term.residue_index,
                gammalooprs::cff::CutCFFIndex::new_all_none(),
                "{context} assigned a residue identity to a 4D forest atom"
            );
        }
        if representation == "4d" {
            let components = local["components"]
                .as_array()
                .unwrap_or_else(|| panic!("{context} has malformed components: {local}"));
            let branches = local["branches"]
                .as_array()
                .unwrap_or_else(|| panic!("{context} has malformed 4D branches: {local}"));
            let mut branch_names = BTreeSet::new();
            for branch in branches {
                let branch_name = branch["branch"]
                    .as_str()
                    .unwrap_or_else(|| panic!("{context} has an unnamed 4D branch: {branch}"));
                assert!(
                    branch["byte_size"].as_u64().is_some_and(|size| size > 0),
                    "{context} recorded a vacuous 4D branch: {branch}"
                );
                assert!(
                    branch_names.insert(branch_name),
                    "{context} duplicated a 4D branch for {}: {branches:?}",
                    term.file_name(),
                );
            }
            if !components.is_empty() {
                assert_eq!(
                    branches
                        .iter()
                        .filter(|branch| branch["branch"] == "Combined")
                        .count(),
                    1,
                    "{context} must expose exactly one completed 4D atom for an active node: {branches:?}",
                );
            }
            for component in components {
                let component_loop_edges =
                    component["route_loop_edges"].as_array().unwrap_or_else(|| {
                        panic!(
                            "{context} has malformed component loop-route provenance: {component}"
                        )
                    });
                assert!(
                    component["route_external_edges"].is_array()
                        && component["route_signatures"].is_array(),
                    "{context} has malformed component route provenance: {component}"
                );
                assert!(
                    !component_loop_edges.is_empty()
                        && component["route_signatures"]
                            .as_array()
                            .is_some_and(|signatures| !signatures.is_empty()),
                    "{context} has an active component without a canonical route: {component}"
                );
            }
            let has_positive_degree_soft_component = components.iter().any(|component| {
                component["scheme"] == "IR" && component["dod"].as_i64().is_some_and(|dod| dod > 0)
            });
            // A multi-parent node is a disconnected join. Its completed component
            // factors have already recorded their own U/S/US splits, so the join
            // itself correctly carries only one Combined branch.
            let is_disconnected_join = parent_keys.len() > 1;
            if has_positive_degree_soft_component && !is_disconnected_join {
                let expected = ["U", "S", "US", "Combined"]
                    .into_iter()
                    .collect::<BTreeSet<_>>();
                assert_eq!(
                    branch_names, expected,
                    "{context} did not expose the complete 4D U/S/US/H branch split"
                );
            }
        } else {
            let paths = local["projection_paths"]
                .as_array()
                .unwrap_or_else(|| panic!("{context} has malformed 3D projection paths: {local}"));
            for path in paths {
                let steps = path["steps"]
                    .as_array()
                    .unwrap_or_else(|| panic!("{context} has a malformed 3D path: {path}"));
                for step in steps {
                    assert!(
                        step["given_component"].is_string()
                            && step["active_subgraph"].is_string()
                            && step["rescaled_subgraph"].is_string()
                            && step["canonical_route_loop_edges"].is_array()
                            && step["canonical_route_external_edges"].is_array()
                            && step["canonical_route_signatures"].is_array()
                            && step["route_loop_edges"].is_array()
                            && step["route_external_edges"].is_array()
                            && step["route_signatures"].is_array(),
                        "{context} has incomplete direct-3D projection context: {step}"
                    );
                    let branches = step["conceptual_branches"].as_array().unwrap_or_else(|| {
                        panic!("{context} has malformed conceptual 3D branches: {step}")
                    });
                    let actual_branches = branches
                        .iter()
                        .map(|branch| {
                            (
                                branch["branch"].as_str().unwrap_or_else(|| {
                                    panic!("{context} has an unnamed 3D branch: {branch}")
                                }),
                                branch["coefficient"].as_i64().unwrap_or_else(|| {
                                    panic!("{context} has an unsigned 3D branch: {branch}")
                                }),
                            )
                        })
                        .collect::<Vec<_>>();
                    if step["scheme"] == "IR" && step["dod"].as_i64().is_some_and(|dod| dod > 0) {
                        assert_eq!(
                            actual_branches,
                            vec![("U", 1), ("S", 1), ("US", -1)],
                            "{context} did not expose the direct-3D U/S/US/H algebra"
                        );
                        assert_eq!(step["materialization"], "factorized_soft");
                    } else {
                        assert_eq!(actual_branches, vec![("U", 1)]);
                        assert_eq!(step["materialization"], "direct_u");
                    }
                }
            }
        }
    }
}

#[derive(Clone, Debug)]
struct CliUvProfileSummary {
    total: usize,
    resolved: usize,
    failed: usize,
    dod_failures: usize,
    orientation_total: usize,
    orientation_resolved: usize,
    orientation_failed: usize,
    orientation_dod_failures: usize,
    orientation_labels: BTreeSet<String>,
    orientation_label_counts: BTreeMap<String, usize>,
    failure_diagnostics: Value,
    profile_identity: Value,
}

fn selected_local_identifiers(
    cli: &CLIState,
    process_name: &str,
    integrand_name: &str,
    scheme: ApproximationType,
) -> Result<Vec<CTIdentifier>> {
    Ok(
        classified_local_spinneys(cli, process_name, integrand_name)?
            .into_iter()
            .filter(|spinney| spinney.n_components == 1 && spinney.scheme == scheme)
            .map(|spinney| spinney.identifier)
            .collect(),
    )
}

fn assert_selected_local_identifier(
    cli: &CLIState,
    process_name: &str,
    integrand_name: &str,
    scheme: ApproximationType,
    external_pdgs: impl IntoIterator<Item = isize>,
    internal_pdgs: impl IntoIterator<Item = isize>,
    minimum_count: usize,
) -> Result<()> {
    let expected = CTIdentifier::new(
        external_pdgs.into_iter().collect(),
        Some(internal_pdgs.into_iter().collect()),
    );
    let identifiers = selected_local_identifiers(cli, process_name, integrand_name, scheme)?;
    let count = identifiers
        .iter()
        .filter(|identifier| **identifier == expected)
        .count();
    assert!(
        count >= minimum_count,
        "expected at least {minimum_count} {scheme} component(s) matching {expected:?}; got {identifiers:?}"
    );
    Ok(())
}

fn assert_generated_local_integrand(
    cli: &CLIState,
    process_name: &str,
    integrand_name: &str,
    minimum_dimensions: usize,
) -> Result<()> {
    let process_id = cli
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified(process_name.to_string())))?;
    let integrand = cli
        .state
        .process_list
        .get_integrand(process_id, integrand_name)?
        .require_generated()?;
    assert!(
        integrand.get_n_dim() >= minimum_dimensions,
        "fixture '{process_name}' did not retain its loop integrations"
    );
    Ok(())
}

fn cli_uv_profile_pass_fail(
    cli: &mut CLIState,
    process: &str,
    integrand: &str,
) -> Result<CliUvProfileSummary> {
    let profiles_per_orientation = !cli
        .cli_settings
        .global
        .generation
        .explicit_orientation_sum_only;
    let profile_command = cli
        .run_history
        .command_blocks
        .iter()
        .find(|block| block.name == "profile_uv")
        .and_then(|block| block.commands.first())
        .and_then(|command| command.raw_string.as_deref())
        .ok_or_else(|| {
            eyre::eyre!("local CT run card must contain a non-empty 'profile_uv' command block")
        })?
        .to_string();
    assert!(
        profile_command.starts_with("profile ultra-violet "),
        "profile_uv must use the ultra-violet CLI syntax; got {profile_command:?}"
    );
    assert!(
        profile_command.contains(&format!("-p {process}"))
            && profile_command.contains(&format!("-i {integrand}")),
        "profile_uv does not target {process}/{integrand}: {profile_command}"
    );
    let command_tokens = profile_command.split_whitespace().collect::<Vec<_>>();
    let option_value = |name| {
        command_tokens
            .windows(2)
            .find_map(|pair| (pair[0] == name).then_some(pair[1]))
    };
    let min_scale_exponent = option_value("--min-scaling")
        .ok_or_else(|| eyre::eyre!("profile_uv has no --min-scaling: {profile_command}"))?
        .parse::<f64>()?;
    let max_scale_exponent = option_value("--max-scaling")
        .ok_or_else(|| eyre::eyre!("profile_uv has no --max-scaling: {profile_command}"))?
        .parse::<f64>()?;
    let n_points = option_value("--n-points")
        .ok_or_else(|| eyre::eyre!("profile_uv has no --n-points: {profile_command}"))?
        .parse::<usize>()?;
    let seed = option_value("--seed")
        .ok_or_else(|| eyre::eyre!("profile_uv has no --seed: {profile_command}"))?
        .parse::<u64>()?;
    assert!(
        min_scale_exponent < max_scale_exponent
            && n_points >= 2
            && seed == 1337
            && command_tokens.contains(&"--per-orientation"),
        "profile_uv must use an increasing per-orientation window with at least two points and seed 1337: {profile_command}"
    );

    // Explicit sums have no independently selectable production orientations.
    // Keep every scale, seed, LMB and numerical threshold of the original profile.
    let profile_command = if profiles_per_orientation {
        profile_command
    } else {
        profile_command.replace(" --per-orientation", "")
    };
    let output = cli.cli_settings.state.folder.join("local_ct_uv_profile");
    // Preserve the exhaustive LMB coverage used by the matched soft-CT controls.
    let Commands::Profile(profile @ Profile::UltraViolet(_)) =
        CommandHistory::from_raw_string(&format!(
            "{profile_command} --selected-limits all --output {}",
            output.display()
        ))?
        .command
    else {
        return Err(eyre::eyre!(
            "profile_uv must parse as an ultra-violet profile"
        ));
    };
    // The CLI handler writes and returns failed-limit reports; the outer
    // command dispatcher converts their verdict into an error. Bare controls
    // need those reports, while parser, evaluation and output errors still fail.
    profile.run(&mut cli.state, &cli.cli_settings)?;
    let profile_json = fs::read_to_string(output.join("uv_profile.json"))?;
    let profile: Value = serde_json::from_str(&profile_json)?;
    let typed_profile: UVProfileAnalysis = serde_json::from_str(&profile_json)?;
    let pass_fail = typed_profile.pass_fail(-0.9);
    let scales = profile["scales"]
        .as_array()
        .ok_or_else(|| eyre::eyre!("CLI UV profile JSON has no scales array"))?;
    assert_eq!(
        scales.len(),
        n_points,
        "profile_uv did not evaluate its requested number of scales"
    );
    assert_eq!(
        scales.first().and_then(Value::as_f64),
        Some(10.0_f64.powf(min_scale_exponent)),
    );
    assert_eq!(
        scales.last().and_then(Value::as_f64),
        Some(10.0_f64.powf(max_scale_exponent)),
    );

    let mut total = 0;
    let mut resolved = 0;
    let mut orientation_total = 0;
    let mut orientation_resolved = 0;
    let mut orientation_labels = BTreeSet::new();
    let mut orientation_label_counts = BTreeMap::new();
    for subset in profile["graphs"]
        .as_array()
        .into_iter()
        .flatten()
        .flat_map(|graph| graph["lmbs"].as_array().into_iter().flatten())
        .flat_map(|lmb| lmb["subsets"].as_array().into_iter().flatten())
    {
        total += 1;
        let orientation_entries = subset["per_orientation_inspect_entries"].as_array();
        if profiles_per_orientation {
            let orientation_entries = orientation_entries.ok_or_else(|| {
                eyre::eyre!(
                    "per-orientation CLI UV profile has no orientation entries for subset: {subset}"
                )
            })?;
            assert!(
                !orientation_entries.is_empty(),
                "per-orientation CLI UV profile selected no orientations for subset: {subset}"
            );
        } else {
            assert!(
                subset["per_orientation_inspect_entries"].is_null(),
                "summed CLI UV profile unexpectedly emitted production-orientation entries: {subset}"
            );
        }
        for entry in orientation_entries.into_iter().flatten() {
            orientation_total += 1;
            let orientation_label = entry["orientation_label"]
                .as_str()
                .ok_or_else(|| eyre::eyre!("CLI UV orientation entry has no label: {entry}"))?;
            assert!(
                !orientation_label.is_empty(),
                "CLI UV orientation entry has an empty label: {entry}"
            );
            orientation_labels.insert(orientation_label.to_string());
            *orientation_label_counts
                .entry(orientation_label.to_string())
                .or_default() += 1;
            let Some(inspect) = entry["analysis"].as_object() else {
                continue;
            };
            orientation_resolved += 1;
            let slope = inspect["result"]["slope"].as_f64().ok_or_else(|| {
                eyre::eyre!(
                    "resolved CLI UV orientation fit {orientation_label} has no finite slope"
                )
            })?;
            let r_squared = inspect["result"]["r_squared"].as_f64().ok_or_else(|| {
                eyre::eyre!(
                    "resolved CLI UV orientation fit {orientation_label} has no finite R-squared"
                )
            })?;
            assert!(
                slope.is_finite() && r_squared.is_finite(),
                "CLI UV orientation {orientation_label} produced a non-finite fit: slope={slope}, R-squared={r_squared}"
            );
        }
        let Some(inspect) = subset["analysis"]["inspect_level"].as_object() else {
            continue;
        };
        resolved += 1;
        let slope = inspect["result"]["slope"]
            .as_f64()
            .ok_or_else(|| eyre::eyre!("resolved CLI UV fit has no finite slope"))?;
        let r_squared = inspect["result"]["r_squared"]
            .as_f64()
            .ok_or_else(|| eyre::eyre!("resolved CLI UV fit has no finite R-squared"))?;
        assert!(
            slope.is_finite() && r_squared.is_finite(),
            "CLI UV profile produced a non-finite fit: slope={slope}, R-squared={r_squared}"
        );
    }
    if profiles_per_orientation {
        assert!(
            orientation_total > 0,
            "per-orientation CLI UV profile selected no orientation-level limits"
        );
    }
    assert_eq!(pass_fail.total, total + orientation_total);
    let dod_failures = pass_fail
        .failures
        .iter()
        .filter(|failure| failure.reason == "dod_exceeds_threshold")
        .count();
    let orientation_failed = pass_fail
        .failures
        .iter()
        .filter(|failure| failure.orientation_label.is_some())
        .count();
    let orientation_dod_failures = pass_fail
        .failures
        .iter()
        .filter(|failure| {
            failure.orientation_label.is_some() && failure.reason == "dod_exceeds_threshold"
        })
        .count();
    let failure_diagnostics = serde_json::to_value(&pass_fail.failures)?;
    let profile_identity = serde_json::json!({
        "scales": profile["scales"].clone(),
        "graphs": profile["graphs"]
            .as_array()
            .into_iter()
            .flatten()
            .map(|graph| serde_json::json!({
                "graph_index": graph["graph_index"].clone(),
                "lmbs": graph["lmbs"]
                    .as_array()
                    .into_iter()
                    .flatten()
                    .map(|lmb| serde_json::json!({
                        "lmb_index": lmb["lmb_index"].clone(),
                        "lmb_label": lmb["lmb_label"].clone(),
                        "subsets": lmb["subsets"]
                            .as_array()
                            .into_iter()
                            .flatten()
                            .map(|subset| serde_json::json!({
                                "subset_index": subset["subset_index"].clone(),
                                "free": subset["free"].clone(),
                                "fixed": subset["fixed"].clone(),
                                "initial_dod": subset["initial_dod"].clone(),
                                "orientation_labels": subset["per_orientation_inspect_entries"]
                                    .as_array()
                                    .into_iter()
                                    .flatten()
                                    .map(|entry| entry["orientation_label"].clone())
                                    .collect::<Vec<_>>(),
                            }))
                            .collect::<Vec<_>>(),
                    }))
                    .collect::<Vec<_>>(),
            }))
            .collect::<Vec<_>>(),
    });
    Ok(CliUvProfileSummary {
        total,
        resolved,
        failed: pass_fail.failed,
        dod_failures,
        orientation_total,
        orientation_resolved,
        orientation_failed,
        orientation_dod_failures,
        orientation_labels,
        orientation_label_counts,
        failure_diagnostics,
        profile_identity,
    })
}

#[derive(Clone, Debug, PartialEq)]
struct CliSoftProfileFit {
    scaling: f64,
    r_squared: f64,
    ray_fingerprint: Value,
}

fn cli_soft_profile_report(
    cli: &mut CLIState,
    process: &str,
    integrand: &str,
    per_orientation: bool,
) -> Result<Value> {
    let profile_command = cli
        .run_history
        .command_blocks
        .iter()
        .find(|block| block.name == "profile_soft")
        .and_then(|block| block.commands.first())
        .and_then(|command| command.raw_string.as_deref())
        .ok_or_else(|| {
            eyre::eyre!("top-bubble run card must contain a non-empty 'profile_soft' command block")
        })?
        .to_string();
    assert!(
        profile_command.starts_with("profile bulk ")
            && profile_command.contains(&format!("-p {process}"))
            && profile_command.contains(&format!("-i {integrand}"))
            && profile_command.contains("--select 'top_bubble_vertex S(e6)'")
            && profile_command.contains("--min-scaling=-2")
            && profile_command.contains("--max-scaling=-5")
            && profile_command.contains("--n-points 25")
            && profile_command.contains("--seed 1337"),
        "profile_soft does not encode the matched q_g ray: {profile_command}"
    );

    let output = cli.cli_settings.state.folder.join(format!(
        "{integrand}_soft{}_profile.json",
        if per_orientation {
            "_per_orientation"
        } else {
            ""
        }
    ));
    cli.run_command(&format!(
        "{profile_command} {}--output {}",
        if per_orientation {
            "--per-orientation "
        } else {
            ""
        },
        output.display()
    ))?;
    let report: Value = serde_json::from_str(&fs::read_to_string(output)?)?;
    assert_eq!(report["settings"]["n_points"], 25);
    assert_eq!(report["settings"]["min_scale_exponent"], -2.0);
    assert_eq!(report["settings"]["max_scale_exponent"], -5.0);
    assert_eq!(report["settings"]["seed"], 1337);
    assert_eq!(report["settings"]["select"], "top_bubble_vertex S(e6)");
    assert_eq!(report["settings"]["per_orientation"], per_orientation);
    let graphs = report["graphs"]
        .as_array()
        .ok_or_else(|| eyre::eyre!("profile_soft JSON has no graph reports"))?;
    assert_eq!(
        graphs.len(),
        1,
        "profile_soft must select exactly one graph"
    );
    assert_eq!(graphs[0]["graph_name"], "top_bubble_vertex");
    Ok(report)
}

fn parse_cli_soft_profile_fit(fit: &Value) -> Result<CliSoftProfileFit> {
    assert_eq!(fit["limit_name"], "S(e6)");
    let scaling = fit["scaling"]
        .as_f64()
        .ok_or_else(|| eyre::eyre!("profile_soft JSON has no finite scaling"))?;
    let r_squared = fit["r_squared"]
        .as_f64()
        .ok_or_else(|| eyre::eyre!("profile_soft JSON has no finite R-squared"))?;
    let ray_fingerprint = fit["ray_fingerprint"].clone();
    assert!(
        ray_fingerprint["lmb_edges"].is_array()
            && ray_fingerprint["digest"]
                .as_str()
                .is_some_and(|digest| digest.len() == 16),
        "profile_soft JSON has no routed-ray fingerprint: {ray_fingerprint}"
    );
    assert_eq!(
        ray_fingerprint["lmb_edges"],
        serde_json::json!([6, 8]),
        "the sampled q_g ray must use the fixture's canonical [e6,e8] route"
    );
    assert!(
        scaling.is_finite() && r_squared.is_finite(),
        "profile_soft returned a non-finite fit: scaling={scaling}, R-squared={r_squared}"
    );
    Ok(CliSoftProfileFit {
        scaling,
        r_squared,
        ray_fingerprint,
    })
}

fn cli_soft_profile_fit(
    cli: &mut CLIState,
    process: &str,
    integrand: &str,
) -> Result<CliSoftProfileFit> {
    let report = cli_soft_profile_report(cli, process, integrand, false)?;
    let reports = report["graphs"][0]["single_limit_reports"]
        .as_array()
        .ok_or_else(|| eyre::eyre!("profile_soft JSON has no limit reports"))?;
    assert_eq!(
        reports.len(),
        1,
        "profile_soft must resolve top_bubble_vertex S(e6) exactly once"
    );
    assert!(reports[0]["orientation_label"].is_null());
    parse_cli_soft_profile_fit(&reports[0])
}

fn cli_soft_profile_orientation_fits(
    cli: &mut CLIState,
    process: &str,
    integrand: &str,
) -> Result<BTreeMap<String, CliSoftProfileFit>> {
    let report = cli_soft_profile_report(cli, process, integrand, true)?;
    let reports = report["graphs"][0]["single_limit_reports"]
        .as_array()
        .ok_or_else(|| eyre::eyre!("profile_soft JSON has no limit reports"))?;
    let fits = reports
        .iter()
        .map(|report| {
            let orientation = report["orientation_label"]
                .as_str()
                .filter(|label| !label.is_empty())
                .ok_or_else(|| eyre::eyre!("per-orientation soft fit has no orientation label"))?
                .to_string();
            Ok((orientation, parse_cli_soft_profile_fit(report)?))
        })
        .collect::<Result<BTreeMap<_, _>>>()?;
    assert_eq!(
        fits.len(),
        reports.len(),
        "per-orientation soft profile contains duplicate orientation labels"
    );
    Ok(fits)
}

fn run_dgse_local_ir_profiles_are_non_vacuous(
    project_local_4d: bool,
    explicit_orientation_sum_only: bool,
) -> Result<()> {
    const BOUNDED_SOFT_SCALING: f64 = 3.0;
    const MINIMUM_SOFT_IMPROVEMENT: f64 = 1.0;
    const MINIMUM_R_SQUARED: f64 = 0.9;
    const EXPECTED_ORIENTATIONS: usize = 30;
    const EXPECTED_UV_LIMITS_PER_ORIENTATION: usize = 27;

    let mut subtracted = local_ct_cli(
        "dgse_local_ir.toml",
        "dgse_local_ir",
        false,
        project_local_4d,
        explicit_orientation_sum_only,
    )?;
    let classified = classified_local_spinneys(&subtracted, "dgse", "dgse")?;
    assert_selected_local_identifier(
        &subtracted,
        "dgse",
        "dgse",
        ApproximationType::IR,
        [-1, 1],
        [1, 21],
        1,
    )?;
    let process_id = subtracted
        .state
        .resolve_process_ref(Some(&ProcessRef::Unqualified("dgse".to_string())))?;
    let integrand_info = subtracted.state.get_integrand_info(
        Some(&ProcessRef::Unqualified("dgse".to_string())),
        Some(&"dgse".to_string()),
    )?;
    // Keep physical-signature uniqueness separate from the report's runtime slot labels.
    assert_eq!(
        integrand_info
            .graph_groups
            .iter()
            .flat_map(|group| &group.orientations)
            .map(|orientation| &orientation.signature)
            .collect::<BTreeSet<_>>()
            .len(),
        EXPECTED_ORIENTATIONS,
        "the unfiltered DGSE fixture must generate every expected acyclic orientation"
    );
    let generated_orientation_labels = integrand_info
        .graph_groups
        .iter()
        .flat_map(|group| &group.orientations)
        .map(|orientation| {
            let label = orientation
                .signature
                .iter()
                .map(|sign| match *sign {
                    -1 => '-',
                    0 => '0',
                    1 => '+',
                    _ => panic!("invalid generated orientation sign {sign}"),
                })
                .collect::<String>();
            format!("{label}|sigma({})", orientation.orientation_id)
        })
        .collect::<BTreeSet<_>>();
    assert_eq!(
        generated_orientation_labels.len(),
        EXPECTED_ORIENTATIONS,
        "the unfiltered DGSE fixture must generate every expected acyclic orientation"
    );
    let soft_export = subtracted.state.process_list.processes[process_id].export_uv_forest_graph(
        "dgse",
        0,
        &gammalooprs::uv::export::UVForestExportSettings { computed: true },
    )?;
    assert_unique_uv_forest_term_identities(&soft_export, "DGSE soft forest");
    let soft_forest_keys = soft_export
        .node_terms
        .iter()
        .map(|term| term.node_key.clone())
        .collect::<BTreeSet<_>>();

    let ir_settings = InfraRedProfile {
        process: Some(ProcessRef::Unqualified("dgse".to_string())),
        integrand_name: Some("dgse".to_string()),
        select: Some("se S(e0)".into()),
        seed: Some(1337),
        ..Default::default()
    };
    let subtracted_ir = Profile::InfraRed(ir_settings.clone())
        .run(&mut subtracted.state, &subtracted.cli_settings)?
        .unwrap_ir();
    let subtracted_ir_reports = subtracted_ir
        .results_per_graph
        .iter()
        .flat_map(|graph| &graph.single_limit_reports)
        .collect::<Vec<_>>();
    assert_eq!(
        subtracted_ir_reports.len(),
        1,
        "IR profile must resolve the requested S(e0) limit exactly once"
    );
    let subtracted_ir_report = subtracted_ir_reports[0];
    let subtracted_ir_r_squared = subtracted_ir_report.power_law_fit.r_squared();
    assert!(
        subtracted_ir_report.scaling.is_finite()
            && subtracted_ir_r_squared.is_finite()
            && subtracted_ir_r_squared >= MINIMUM_R_SQUARED,
        "subtracted DGSE soft profile must have a finite, resolved fit; scaling={}, R-squared={subtracted_ir_r_squared}",
        subtracted_ir_report.scaling,
    );
    let subtracted_ir_ray = subtracted_ir_report
        .ray_fingerprint
        .clone()
        .expect("subtracted DGSE soft profile must record its routed ray");
    assert!(
        !subtracted_ir_ray.lmb_edges.is_empty() && !subtracted_ir_ray.digest.is_empty(),
        "subtracted DGSE soft profile must have a non-empty routed-ray fingerprint: {subtracted_ir_ray:?}"
    );
    assert!(
        subtracted_ir_report.scaling > BOUNDED_SOFT_SCALING,
        "local IR subtraction must make the pointwise integrand bounded in S(e0); got scaling {}",
        subtracted_ir_report.scaling
    );
    let subtracted_ir_scaling = subtracted_ir_report.scaling;
    let subtracted_uv = cli_uv_profile_pass_fail(&mut subtracted, "dgse", "dgse")?;
    assert_eq!(
        subtracted_uv.total, EXPECTED_UV_LIMITS_PER_ORIENTATION,
        "the DGSE UV profile selected an unexpected hard-limit inventory"
    );
    assert_eq!(
        subtracted_uv.resolved, subtracted_uv.total,
        "every summed DGSE UV limit must resolve"
    );
    if !explicit_orientation_sum_only {
        assert_eq!(
            subtracted_uv.orientation_labels, generated_orientation_labels,
            "the per-orientation UV profile must cover exactly the generated DGSE orientations"
        );
        assert_eq!(
            subtracted_uv.orientation_total,
            EXPECTED_ORIENTATIONS * EXPECTED_UV_LIMITS_PER_ORIENTATION,
            "the DGSE UV profile must evaluate every hard limit in every orientation"
        );
        assert_eq!(
            subtracted_uv.orientation_label_counts,
            generated_orientation_labels
                .iter()
                .cloned()
                .map(|label| (label, EXPECTED_UV_LIMITS_PER_ORIENTATION))
                .collect::<BTreeMap<_, _>>(),
            "every generated DGSE orientation must occur once in each hard limit",
        );
        assert_eq!(
            subtracted_uv.orientation_resolved, subtracted_uv.orientation_total,
            "every orientation-local DGSE UV fit must resolve"
        );
    }
    if subtracted_uv.failed != 0 {
        let mut uv_only = local_ct_cli(
            "dgse_uv_only.toml",
            "dgse_uv_only_diagnostic",
            false,
            project_local_4d,
            explicit_orientation_sum_only,
        )?;
        let uv_process_id = uv_only
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified("dgse".to_string())))?;
        let uv_export = uv_only.state.process_list.processes[uv_process_id]
            .export_uv_forest_graph(
                "dgse",
                0,
                &gammalooprs::uv::export::UVForestExportSettings { computed: true },
            )?;
        assert_unique_uv_forest_term_identities(&uv_export, "DGSE ordinary-U forest");
        let uv_forest_keys = uv_export
            .node_terms
            .iter()
            .map(|term| term.node_key.clone())
            .collect::<BTreeSet<_>>();
        let dump = subtracted
            .cli_settings
            .state
            .folder
            .join("dgse_soft_uv_forest_comparison");
        fs::create_dir_all(&dump)?;
        for (prefix, export) in [("soft", &soft_export), ("uv", &uv_export)] {
            fs::write(
                dump.join(format!("{prefix}_forest.dot")),
                &export.forest_dot,
            )?;
            for term in &export.node_terms {
                fs::write(
                    dump.join(format!("{prefix}_{}", term.file_name())),
                    &term.dot,
                )?;
            }
        }
        let uv = cli_uv_profile_pass_fail(&mut uv_only, "dgse", "dgse")?;
        clean_test(&uv_only.cli_settings.state.folder);
        assert_eq!(
            subtracted_uv.failed,
            0,
            "local IR subtraction failed its UV profile; classified={classified:?}, soft forest={soft_forest_keys:?}, failures={}; matched ordinary-U control: total={}, resolved={}, orientation_total={}, orientation_resolved={}, failed={}, orientation_failed={}, forest={uv_forest_keys:?}",
            subtracted_uv.failure_diagnostics,
            uv.total,
            uv.resolved,
            uv.orientation_total,
            uv.orientation_resolved,
            uv.failed,
            uv.orientation_failed,
        );
    }

    let mut bare = local_ct_cli(
        "dgse_local_ir.toml",
        "dgse_local_ir_bare",
        true,
        project_local_4d,
        explicit_orientation_sum_only,
    )?;
    let bare_ir = Profile::InfraRed(ir_settings)
        .run(&mut bare.state, &bare.cli_settings)?
        .unwrap_ir();
    let bare_ir_reports = bare_ir
        .results_per_graph
        .iter()
        .flat_map(|graph| &graph.single_limit_reports)
        .collect::<Vec<_>>();
    assert_eq!(
        bare_ir_reports.len(),
        1,
        "bare IR profile must resolve the requested S(e0) limit exactly once"
    );
    let bare_ir_report = bare_ir_reports[0];
    let bare_ir_r_squared = bare_ir_report.power_law_fit.r_squared();
    assert!(
        bare_ir_report.scaling.is_finite()
            && bare_ir_r_squared.is_finite()
            && bare_ir_r_squared >= MINIMUM_R_SQUARED,
        "bare DGSE soft profile must have a finite, resolved fit; scaling={}, R-squared={bare_ir_r_squared}",
        bare_ir_report.scaling,
    );
    let bare_ir_ray = bare_ir_report
        .ray_fingerprint
        .as_ref()
        .expect("bare DGSE soft profile must record its routed ray");
    assert_eq!(
        &subtracted_ir_ray, bare_ir_ray,
        "subtracted and bare DGSE soft controls must evaluate the identical routed ray",
    );
    assert!(
        bare_ir_report.scaling <= BOUNDED_SOFT_SCALING,
        "bare dgse fixture must retain a pointwise soft singularity; got scaling {}",
        bare_ir_report.scaling
    );
    assert!(
        subtracted_ir_scaling - bare_ir_report.scaling >= MINIMUM_SOFT_IMPROVEMENT,
        "local IR subtraction must improve the matched DGSE S(e0) scaling by at least one power; bare={}, subtracted={subtracted_ir_scaling}",
        bare_ir_report.scaling,
    );
    let bare_uv = cli_uv_profile_pass_fail(&mut bare, "dgse", "dgse")?;
    assert!(bare_uv.total > 0, "bare UV profile selected no limits");
    assert!(
        bare_uv.resolved > 0
            && (explicit_orientation_sum_only
                || (bare_uv.orientation_total > 0 && bare_uv.orientation_resolved > 0)),
        "bare UV profile produced no resolved summed and per-orientation fits"
    );
    assert!(
        bare_uv.dod_failures > 0
            && (explicit_orientation_sum_only || (bare_uv.orientation_dod_failures > 0)),
        "bare dgse fixture must have a resolved per-orientation fit above the UV bound; summary={bare_uv:?}"
    );
    assert_eq!(
        subtracted_uv.profile_identity, bare_uv.profile_identity,
        "subtracted and bare DGSE controls must profile the same UV LMBs, subsets, and orientations",
    );

    clean_test(&subtracted.cli_settings.state.folder);
    clean_test(&bare.cli_settings.state.folder);
    Ok(())
}

#[test]
#[serial_test::serial]
fn dgse_local_ir_profiles_are_non_vacuous() -> Result<()> {
    for (explicit_orientation_sum_only, project_local_4d) in
        [(false, false), (true, false), (true, true)]
    {
        println!(
            "explicit_orientation_sum_only={explicit_orientation_sum_only}, project_local_4d={project_local_4d}"
        );
        let test = std::thread::Builder::new()
            .name("dgse-local-ir-profile".to_string())
            .stack_size(128 * 1024 * 1024)
            .spawn(move || {
                run_dgse_local_ir_profiles_are_non_vacuous(
                    project_local_4d,
                    explicit_orientation_sum_only,
                )
            })?;
        test.join()
            .map_err(|_| eyre::eyre!("DGSE local-IR profile thread panicked"))??;
    }
    Ok(())
}

#[test]
#[serial_test::serial]
fn dgse_local_ir_generation_orientation_filter_profiles_one_orientation() -> Result<()> {
    let project_local_4d = false;
    // Runtime production-orientation indices are exposed by the localized direct route.
    const UNFILTERED_ORIENTATION_INDEX: usize = 3;
    const PATTERN: &str = "(-,-,0,0,+,+,+,0,+)";
    const SIGNATURE: [i8; 9] = [-1, -1, 0, 0, 1, 1, 1, 0, 1];
    const LABEL: &str = "--00+++0+";
    const EXPECTED_UV_LIMITS: usize = 27;

    let test = std::thread::Builder::new()
        .name("dgse-local-ir-filtered-orientation-profile".to_string())
        .stack_size(128 * 1024 * 1024)
        .spawn(move || -> Result<()> {
            let test_name = format!("dgse_local_ir_generation_orientation_filter_project_local_4d_{project_local_4d}");
            let mut cli = get_test_cli(
                Some("dgse_local_ir.toml".into()),
                get_tests_workspace_path().join(&test_name),
                Some(test_name.to_string()),
                true,
            )?;
            let generation = &mut cli.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = false;
            generation.uv.final_integrand = FinalIntegrandDimension::ThreeD;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d;
            assert!(
                !cli.cli_settings.global.generation.uv.generate_integrated,
                "the focused local-IR orientation profile must remain local-only"
            );

            cli.run_command("run generate")?;
            let process = ProcessRef::Unqualified("dgse".to_string());
            let integrand_name = "dgse".to_string();
            let unfiltered_info = cli
                .state
                .get_integrand_info(Some(&process), Some(&integrand_name))?;
            assert!(
                !unfiltered_info.graph_groups.is_empty()
                    && unfiltered_info.graph_groups.iter().all(|group| {
                        group
                            .orientations
                            .get(UNFILTERED_ORIENTATION_INDEX)
                            .is_some_and(|orientation| orientation.signature == SIGNATURE)
                    }),
                "unfiltered DGSE orientation {UNFILTERED_ORIENTATION_INDEX} must be {LABEL}",
            );
            let point = deterministic_uv_momentum_points(&mut cli, "dgse", "dgse")?
                .remove(0)
                .point;
            let unfiltered_floor =
                required_stability_accuracy_floor(&mut cli, "dgse", "dgse")?;
            let unfiltered = evaluate_momentum_sample(
                &mut cli,
                "dgse",
                "dgse",
                &point,
                Some(UNFILTERED_ORIENTATION_INDEX),
                unfiltered_floor,
            )?;

            cli.run_command(&format!(
                "set global string '\n[global.generation.orientation_pattern]\npat = \"{PATTERN}\"\n'"
            ))?;
            cli.run_command("run generate")?;

            let info = cli
                .state
                .get_integrand_info(Some(&process), Some(&integrand_name))?;
            assert!(
                !info.graph_groups.is_empty()
                    && info.graph_groups.iter().all(|group| {
                        group.orientations.len() == 1
                            && group.orientations[0].orientation_id == 0
                            && group.orientations[0].signature == SIGNATURE
                    }),
                "the generation-time pattern must retain only {LABEL}: {:?}",
                info.graph_groups
                    .iter()
                    .map(|group| &group.orientations)
                    .collect::<Vec<_>>(),
            );
            let filtered_floor =
                required_stability_accuracy_floor(&mut cli, "dgse", "dgse")?;
            let filtered = evaluate_momentum_sample(
                &mut cli,
                "dgse",
                "dgse",
                &point,
                Some(0),
                filtered_floor,
            )?;
            let delta = (filtered.value.re - unfiltered.value.re)
                .hypot(filtered.value.im - unfiltered.value.im);
            let scale = filtered
                .value
                .re
                .hypot(filtered.value.im)
                .max(unfiltered.value.re.hypot(unfiltered.value.im))
                .max(f64::MIN_POSITIVE);
            let relative_delta = delta / scale;
            let relative_accuracy = filtered
                .relative_accuracy
                .max(unfiltered.relative_accuracy)
                .max(f64::EPSILON);
            assert!(
                filtered.value.re.is_finite()
                    && filtered.value.im.is_finite()
                    && unfiltered.value.re.is_finite()
                    && unfiltered.value.im.is_finite()
                    && INSPECT_DEPENDENCE_ACCURACY_FACTOR * relative_accuracy < 1.0
                    && relative_delta <= INSPECT_DEPENDENCE_ACCURACY_FACTOR * relative_accuracy,
                "filtered visible orientation 0 differs from unfiltered orientation {UNFILTERED_ORIENTATION_INDEX}: filtered={:?}, unfiltered={:?}, relative delta={relative_delta:.3e}, accuracy={relative_accuracy:.3e}",
                filtered.value,
                unfiltered.value,
            );

            // UV reports identify the visible runtime slot as well as its physical signs.
            let profile_label = format!(
                "{LABEL}|sigma({})",
                info.graph_groups[0].orientations[0].orientation_id
            );
            let profile = cli_uv_profile_pass_fail(&mut cli, "dgse", "dgse")?;
            assert_eq!(
                profile.orientation_labels,
                BTreeSet::from([profile_label.clone()]),
                "the filtered UV profile must contain only the generated orientation"
            );
            assert_eq!(
                profile.orientation_label_counts,
                BTreeMap::from([(profile_label, EXPECTED_UV_LIMITS)]),
                "the filtered orientation must occur once in every hard limit"
            );
            assert!(
                profile.total == EXPECTED_UV_LIMITS
                    && profile.resolved == profile.total
                    && profile.orientation_total == EXPECTED_UV_LIMITS
                    && profile.orientation_resolved == profile.orientation_total
                    && profile.failed == 0
                    && profile.orientation_failed == 0,
                "the focused DGSE orientation failed its UV subtraction: {profile:?}"
            );

            clean_test(&cli.cli_settings.state.folder);
            Ok(())
        })?;
    test.join()
        .map_err(|_| eyre::eyre!("filtered DGSE orientation profile thread panicked"))??;

    Ok(())
}

#[test]
#[serial_test::serial]
fn massless_quark_self_energy_local_ir_profile_is_non_vacuous() -> Result<()> {
    for (explicit_orientation_sum_only, project_local_4d) in
        [(false, false), (true, false), (true, true)]
    {
        println!(
            "explicit_orientation_sum_only={explicit_orientation_sum_only}, project_local_4d={project_local_4d}"
        );
        const CARD: &str = "local_ir_quark_self_energy.toml";
        const NAME: &str = "local_ir_quark_self_energy";

        let mut subtracted = local_ct_cli(
            CARD,
            NAME,
            false,
            project_local_4d,
            explicit_orientation_sum_only,
        )?;
        assert_selected_local_identifier(
            &subtracted,
            NAME,
            NAME,
            ApproximationType::IR,
            [-1, 1],
            [1, 21],
            1,
        )?;
        let subtracted_uv = cli_uv_profile_pass_fail(&mut subtracted, NAME, NAME)?;
        assert!(subtracted_uv.total > 0, "IR UV profile selected no limits");
        assert!(
            subtracted_uv.resolved > 0
                && (explicit_orientation_sum_only
                    || (subtracted_uv.orientation_total > 0
                        && subtracted_uv.orientation_resolved > 0)),
            "local IR subtraction produced no resolved summed and per-orientation UV fits"
        );
        assert_eq!(
            (subtracted_uv.failed, subtracted_uv.orientation_failed),
            (0, 0),
            "primitive local IR subtraction failed its UV profile"
        );

        let mut bare = local_ct_cli(
            CARD,
            "local_ir_quark_self_energy_bare",
            true,
            project_local_4d,
            explicit_orientation_sum_only,
        )?;
        assert_matched_local_routes(&subtracted, &bare, NAME, NAME)?;
        let bare_uv = cli_uv_profile_pass_fail(&mut bare, NAME, NAME)?;
        assert!(bare_uv.total > 0, "bare IR UV profile selected no limits");
        assert!(
            bare_uv.resolved > 0
                && (explicit_orientation_sum_only
                    || (bare_uv.orientation_total > 0 && bare_uv.orientation_resolved > 0)),
            "bare IR profile produced no resolved summed and per-orientation UV fits"
        );
        assert!(
            bare_uv.dod_failures > 0
                && (explicit_orientation_sum_only || (bare_uv.orientation_dod_failures > 0)),
            "bare massless quark self-energy must have a resolved per-orientation fit above the UV bound; summary={bare_uv:?}"
        );
        assert_eq!(
            subtracted_uv.profile_identity, bare_uv.profile_identity,
            "subtracted and bare massless-quark controls must profile the same UV LMBs, subsets, and orientations",
        );

        clean_test(&subtracted.cli_settings.state.folder);
        clean_test(&bare.cli_settings.state.folder);
    }
    Ok(())
}

#[test]
#[serial_test::serial]
fn massless_gluon_self_energy_local_ir_profile_is_non_vacuous() -> Result<()> {
    for (explicit_orientation_sum_only, project_local_4d) in
        [(false, false), (true, false), (true, true)]
    {
        println!(
            "explicit_orientation_sum_only={explicit_orientation_sum_only}, project_local_4d={project_local_4d}"
        );
        const CARD: &str = "local_ir_gluon_self_energy.toml";
        const NAME: &str = "local_ir_gluon_self_energy";

        let mut subtracted = local_ct_cli(
            CARD,
            NAME,
            false,
            project_local_4d,
            explicit_orientation_sum_only,
        )?;
        assert_selected_local_identifier(
            &subtracted,
            NAME,
            NAME,
            ApproximationType::IR,
            [-21, 21],
            [1],
            1,
        )?;
        let subtracted_uv = cli_uv_profile_pass_fail(&mut subtracted, NAME, NAME)?;
        assert!(
            subtracted_uv.total > 0
                && subtracted_uv.resolved > 0
                && subtracted_uv.failed == 0
                && (explicit_orientation_sum_only
                    || (subtracted_uv.orientation_total > 0
                        && subtracted_uv.orientation_resolved > 0
                        && subtracted_uv.orientation_failed == 0)),
            "massless-gluon H_2 subtraction must pass a non-vacuous UV profile: {subtracted_uv:?}"
        );

        let mut bare = local_ct_cli(
            CARD,
            "local_ir_gluon_self_energy_bare",
            true,
            project_local_4d,
            explicit_orientation_sum_only,
        )?;
        assert_matched_local_routes(&subtracted, &bare, NAME, NAME)?;
        let bare_uv = cli_uv_profile_pass_fail(&mut bare, NAME, NAME)?;
        assert!(
            bare_uv.total > 0
                && bare_uv.resolved > 0
                && bare_uv.dod_failures > 0
                && (explicit_orientation_sum_only
                    || (bare_uv.orientation_total > 0
                        && bare_uv.orientation_resolved > 0
                        && bare_uv.orientation_dod_failures > 0)),
            "bare massless gluon self-energy must have a resolved per-orientation fit above the UV bound: {bare_uv:?}"
        );
        assert_eq!(
            subtracted_uv.profile_identity, bare_uv.profile_identity,
            "subtracted and bare massless-gluon controls must profile the same UV LMBs, subsets, and orientations",
        );

        clean_test(&subtracted.cli_settings.state.folder);
        clean_test(&bare.cli_settings.state.folder);
    }
    Ok(())
}

const SOFT_CFF_PERSISTENCE_MODE: &str = "GAMMALOOP_SOFT_CFF_PERSISTENCE_MODE";
const SOFT_CFF_PERSISTENCE_STATE: &str = "GAMMALOOP_SOFT_CFF_PERSISTENCE_STATE";
const SOFT_CFF_PERSISTENCE_PROJECT_LOCAL_4D: &str =
    "GAMMALOOP_SOFT_CFF_PERSISTENCE_PROJECT_LOCAL_4D";

#[test]
#[ignore = "spawned by soft_cff_state_survives_a_fresh_process_reload"]
fn soft_cff_state_round_trip_child() -> Result<()> {
    const CARD: &str = "local_ir_gluon_self_energy.toml";
    const NAME: &str = "local_ir_gluon_self_energy";

    let Ok(mode) = env::var(SOFT_CFF_PERSISTENCE_MODE) else {
        return Ok(());
    };
    let state_path = env::var_os(SOFT_CFF_PERSISTENCE_STATE)
        .map(PathBuf::from)
        .ok_or_else(|| eyre::eyre!("missing {SOFT_CFF_PERSISTENCE_STATE} child state path"))?;

    let project_local_4d: bool = env::var(SOFT_CFF_PERSISTENCE_PROJECT_LOCAL_4D)?.parse()?;
    let worker = std::thread::Builder::new()
        .name(format!("soft-cff-persistence-{mode}"))
        .stack_size(128 * 1024 * 1024)
        .spawn(move || -> Result<()> {
            match mode.as_str() {
                "save" => {
                    let mut cli = get_test_cli(
                        Some(CARD.into()),
                        &state_path,
                        Some("soft_cff_persistence_save".to_string()),
                        true,
                    )?;
                    let generation = &mut cli.cli_settings.global.generation;
                    generation.explicit_orientation_sum_only = true;
                    generation.uv.final_integrand = FinalIntegrandDimension::ThreeD;
                    generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d;
                    assert!(!generation.uv.generate_integrated);
                    cli.run_command("run generate")?;
                    assert_generated_local_integrand(&cli, NAME, NAME, 3)?;
                    cli.save_state()?;
                }
                "load" => {
                    let mut loaded = StateLoadOption {
                        state_folder: Some(state_path),
                        ..Default::default()
                    }
                    .load()?;
                    assert_eq!(
                        loaded.cli_settings.global.generation.uv.local_uv_cts_from_expanded_4d_integrands,
                        project_local_4d,
                        "the fresh process must retain the requested local CT route",
                    );
                    assert!(loaded.cli_settings.global.generation.explicit_orientation_sum_only);
                    assert!(
                loaded.state_load_summary.is_some(),
                "the load child started a blank state instead of deserializing the saved soft CFF"
            );
                    loaded
                        .cli_session()
                        .execute_command(CommandHistory::from_raw_string(&format!(
                            "generate existing -p {NAME} -i {NAME}"
                        ))?)?;

                    let process_id = loaded
                        .state
                        .resolve_process_ref(Some(&ProcessRef::Unqualified(NAME.to_string())))?;
                    let integrand = loaded
                        .state
                        .process_list
                        .get_integrand(process_id, NAME)?
                        .require_generated()?;
                    let export = loaded.state.process_list.processes[process_id].export_uv_forest_graph(
NAME,
0,
                        &gammalooprs::uv::export::UVForestExportSettings {
                            computed: true,
                        },
                    )?;
                    assert_unique_uv_forest_term_identities(
                        &export,
                        "fresh-process soft-CFF state",
                    );

                    let n_dim = integrand.get_n_dim();
                    let graph_name = integrand
                        .graph_name_by_id(0)
                        .ok_or_else(|| eyre::eyre!("reloaded soft-CFF integrand has no graph 0"))?
                        .to_string();
                    let point = (0..n_dim)
                        .map(|index| 0.13 + 0.07 * index as f64)
                        .collect::<Vec<_>>();
                    let points = Array2::from_shape_vec((1, n_dim), point)?;
                    let result = evaluate_sample(
                        &mut loaded.state,
                        &EvaluateSamples {
                            process_id: Some(process_id),
                            integrand_name: Some(NAME.to_string()),
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
                    let value = evaluation.integrand_result.map(|entry| entry.0);
                    assert!(
                value.re.is_finite()
                    && value.im.is_finite()
                    && evaluation.evaluation_metadata.is_some(),
                "the freshly loaded and regenerated soft CFF did not evaluate normally: {value:?}"
            );
                }
                other => {
                    return Err(eyre::eyre!(
                        "unknown {SOFT_CFF_PERSISTENCE_MODE} child mode {other:?}"
                    ));
                }
            }
            Ok(())
        })?;
    worker
        .join()
        .map_err(|_| eyre::eyre!("soft-CFF persistence child worker panicked"))?
}

#[test]
#[serial_test::serial]
fn soft_cff_state_survives_a_fresh_process_reload() -> Result<()> {
    for project_local_4d in [false, true] {
        println!("project_local_4d={project_local_4d}");
        let state_path = get_tests_workspace_path().join(format!(
            "soft_cff_fresh_process_state_project_local_4d_{project_local_4d}"
        ));
        clean_test(&state_path);
        let outcome = (|| -> Result<()> {
            let executable = env::current_exe()?;
            for mode in ["save", "load"] {
                let output = Command::new(&executable)
                    .current_dir(workspace_root())
                    .args([
                        "--exact",
                        "soft_cff_state_round_trip_child",
                        "--ignored",
                        "--nocapture",
                    ])
                    .env(SOFT_CFF_PERSISTENCE_MODE, mode)
                    .env(SOFT_CFF_PERSISTENCE_STATE, &state_path)
                    .env(
                        SOFT_CFF_PERSISTENCE_PROJECT_LOCAL_4D,
                        project_local_4d.to_string(),
                    )
                    .output()?;
                if !output.status.success() {
                    return Err(eyre::eyre!(
                        "fresh-process soft-CFF {mode} child failed with {}\nstdout:\n{}\nstderr:\n{}",
                        output.status,
                        String::from_utf8_lossy(&output.stdout),
                        String::from_utf8_lossy(&output.stderr),
                    ));
                }
            }
            Ok(())
        })();
        clean_test(&state_path);
        outcome?;
    }
    Ok(())
}

const LOCAL_OS_CARD: &str = "local_os_top_self_energy.toml";
const LOCAL_OS_NAME: &str = "local_os_top_self_energy";

#[test]
#[serial_test::serial]
fn massive_top_self_energy_local_os_dispatch_is_deferred() {
    let state_path = get_tests_workspace_path().join(format!(
        "{LOCAL_OS_NAME}_explicit_false_project_local_4d_false"
    ));
    clean_test(&state_path);
    let panic = catch_unwind(AssertUnwindSafe(|| {
        let _ = local_ct_cli(LOCAL_OS_CARD, LOCAL_OS_NAME, false, false, false);
    }))
    .expect_err("the local OS fixture must reach the intentional deferred-dispatch panic");
    clean_test(&state_path);
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .unwrap_or("non-string panic payload");
    assert!(
        message.contains(
            "local on-shell counterterms are deferred until local counterterms can be derived from the 4D expanded representation"
        ),
        "unexpected local OS panic: {message}",
    );
}

#[test]
#[serial_test::serial]
fn paper_appendix_b1_os_child_reaches_deferred_dispatch() {
    const NAME: &str = "paper_appendix_b1_nested_gluon_self_energy_os";
    let state_path =
        get_tests_workspace_path().join(format!("{NAME}_explicit_false_project_local_4d_false"));
    clean_test(&state_path);
    let panic = catch_unwind(AssertUnwindSafe(|| {
        let _ = local_ct_cli(
            "paper_appendix_b1_nested_gluon_self_energy_os.toml",
            NAME,
            false,
            false,
            false,
        );
    }))
    .expect_err("the Appendix-B.1 OS fixture must reach the intentional deferred-dispatch panic");
    clean_test(&state_path);
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .unwrap_or("non-string panic payload");
    assert!(
        message.contains(
            "local on-shell counterterms are deferred until local counterterms can be derived from the 4D expanded representation"
        ),
        "unexpected Appendix-B.1 OS panic: {message}",
    );
}

#[test]
#[serial_test::serial]
fn local_ir_disconnected_completed_components_form_nonzero_union() -> Result<()> {
    for (explicit_orientation_sum_only, project_local_4d) in
        [(false, false), (true, false), (true, true)]
    {
        println!(
            "explicit_orientation_sum_only={explicit_orientation_sum_only}, project_local_4d={project_local_4d}"
        );
        let test = std::thread::Builder::new()
        .name("local-ir-disconnected-replay".to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(move || -> Result<()> {
            use gammalooprs::graph::parse::IntoGraph;

            const NAME: &str = "local_ir_dod2_scalar_spectacles";
            let mut spectacles =
                local_ct_cli("local_ir_dod2_scalar_spectacles.toml", NAME, false, project_local_4d, explicit_orientation_sum_only)?;
            assert_generated_local_integrand(&spectacles, NAME, NAME, 6)?;
            let classified = classified_local_spinneys(&spectacles, NAME, NAME)?;
            let connected = classified
                .iter()
                .filter(|spinney| spinney.n_components == 1)
                .collect::<Vec<_>>();
            assert!(
                connected.len() == 2
                    && connected.iter().all(|spinney| {
                        spinney.scheme == ApproximationType::IR && spinney.dod == 2
                    }),
                "dotted spectacles must select exactly two degree-2 soft components: {connected:?}"
            );
            let disconnected = classified
                .iter()
                .filter(|spinney| spinney.n_components > 1)
                .collect::<Vec<_>>();
            assert!(
                !disconnected.is_empty()
                    && disconnected
                        .iter()
                        .all(|spinney| spinney.scheme == ApproximationType::MUV),
                "disconnected unions must use only the neutral MUV tag: {disconnected:?}"
            );

            let process_id = spectacles
                .state
                .resolve_process_ref(Some(&ProcessRef::Unqualified(NAME.to_string())))?;
            let export = spectacles.state.process_list.processes[process_id].export_uv_forest_graph(
NAME,
0,
                &gammalooprs::uv::export::UVForestExportSettings {
                    computed: true,
                },
            )?;
            assert_unique_uv_forest_term_identities(
                &export,
                "disconnected spectacles forest",
            );
            let mut nodes = BTreeMap::<String, Atom>::new();
            for term in export.node_terms {
                let graph: gammalooprs::graph::Graph =
                    term.dot.as_str().into_graph(&spectacles.state.model)?;
                nodes
                    .entry(term.node_key)
                    .and_modify(|sum| *sum += &graph.global_prefactor.num)
                    .or_insert(graph.global_prefactor.num);
            }
            assert_eq!(
                nodes.len(),
                4,
                "fixture must expose a root, two completed components, and their union"
            );

            let root = nodes
                .remove("∅")
                .expect("the exported forest must contain its root node");
            let union_key = nodes
                .keys()
                .find(|key| key.contains(','))
                .cloned()
                .expect("the exported forest must contain a two-component union");
            let union = nodes
                .remove(&union_key)
                .expect("the identified union must remain in the node map");
            let components = nodes.into_values().collect::<Vec<_>>();
            assert_eq!(
                components.len(),
                2,
                "the disconnected union must have exactly two completed component nodes"
            );
            assert!(
                !root.is_zero()
                    && !union.is_zero()
                    && components.iter().all(|component| !component.is_zero()),
                "the disconnected-union assertion must not be satisfied by vanishing forest nodes"
            );

            // HedgePoset replays both connected-component operators with
            // independent local hard parameters. Its unit test checks the
            // resulting two-component CFF sector against both explicit replay
            // orders in one parent orientation. Here we additionally require
            // every exported family to be nonzero before checking the summed
            // and individual-orientation UV profiles below.

            let subtracted_uv = cli_uv_profile_pass_fail(&mut spectacles, NAME, NAME)?;
            let mut bare = local_ct_cli(
                "local_ir_dod2_scalar_spectacles.toml",
                "local_ir_dod2_scalar_spectacles_bare",
                true, project_local_4d, explicit_orientation_sum_only)?;
            let bare_uv = cli_uv_profile_pass_fail(&mut bare, NAME, NAME)?;
            let mut ordinary_uv = get_test_cli(
                Some("local_ir_dod2_scalar_spectacles.toml".into()),
                get_tests_workspace_path().join(format!("local_ir_dod2_scalar_spectacles_muv_explicit_{explicit_orientation_sum_only}_project_local_4d_{project_local_4d}")),
                Some(format!("local_ir_dod2_scalar_spectacles_muv_explicit_{explicit_orientation_sum_only}_project_local_4d_{project_local_4d}")),
                true,
            )?;
            let generation = &mut ordinary_uv.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = explicit_orientation_sum_only;
            generation.uv.final_integrand = FinalIntegrandDimension::ThreeD;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d;
            let ordinary_prescription = &mut ordinary_uv
                .cli_settings
                .global
                .generation
                .uv
                .renormalization_prescription;
            ordinary_prescription.log_divergent = ApproximationType::MUV;
            ordinary_prescription.massive_power_divergent = ApproximationType::MUV;
            ordinary_prescription.massless_power_divergent = ApproximationType::MUV;
            ordinary_prescription.overrides.clear();
            ordinary_uv.run_command("run generate")?;
            let ordinary_components = classified_local_spinneys(&ordinary_uv, NAME, NAME)?
                .into_iter()
                .filter(|spinney| spinney.n_components == 1)
                .map(|spinney| (spinney.edge_ids, spinney.dod, spinney.scheme))
                .collect::<Vec<_>>();
            assert!(
                ordinary_components.len() == 2
                    && ordinary_components.iter().all(|(_, dod, scheme)| {
                        *dod == 2 && *scheme == ApproximationType::MUV
                    }),
                "ordinary spectacles control must select exactly both d=2 components as MUV: {ordinary_components:?}",
            );
            let ordinary_uv_profile = cli_uv_profile_pass_fail(&mut ordinary_uv, NAME, NAME)?;
            assert_eq!(
                subtracted_uv.profile_identity, bare_uv.profile_identity,
                "subtracted and bare disconnected controls must profile the same UV LMBs, subsets, and orientations",
            );
            assert_eq!(
                subtracted_uv.profile_identity, ordinary_uv_profile.profile_identity,
                "soft-refined and ordinary-MUV spectacles controls must profile the same UV LMBs, subsets, and orientations",
            );
            assert!(
subtracted_uv.total > 0 && subtracted_uv.resolved > 0 && subtracted_uv.failed == 0
 && (explicit_orientation_sum_only || (subtracted_uv.orientation_total > 0 && subtracted_uv.orientation_resolved > 0 && subtracted_uv.orientation_failed == 0)),
                "completed disconnected soft components must pass a non-vacuous CLI UV profile: {subtracted_uv:?}"
            );
            assert!(
ordinary_uv_profile.total > 0 && ordinary_uv_profile.resolved > 0 && ordinary_uv_profile.failed == 0
 && (explicit_orientation_sum_only || (ordinary_uv_profile.orientation_total > 0 && ordinary_uv_profile.orientation_resolved > 0 && ordinary_uv_profile.orientation_failed == 0)),
                "ordinary-MUV disconnected components must pass the same non-vacuous CLI UV profile: {ordinary_uv_profile:?}"
            );
            assert!(
bare_uv.total > 0 && bare_uv.resolved > 0 && bare_uv.dod_failures > 0
 && (explicit_orientation_sum_only || (bare_uv.orientation_total > 0 && bare_uv.orientation_resolved > 0 && bare_uv.orientation_dod_failures > 0)),
                "bare disconnected fixture must have a resolved per-orientation fit above the UV bound: {bare_uv:?}"
            );

            clean_test(&spectacles.cli_settings.state.folder);
            clean_test(&bare.cli_settings.state.folder);
            clean_test(&ordinary_uv.cli_settings.state.folder);
            Ok(())
        })?;
        test.join()
            .map_err(|_| eyre::eyre!("disconnected replay test thread panicked"))??;
    }
    Ok(())
}

#[test]
#[serial_test::serial]
fn consecutive_soft_components_have_exact_forest_and_cli_profile() -> Result<()> {
    for (explicit_orientation_sum_only, project_local_4d) in
        [(false, false), (true, false), (true, true)]
    {
        println!(
            "explicit_orientation_sum_only={explicit_orientation_sum_only}, project_local_4d={project_local_4d}"
        );
        let test = std::thread::Builder::new()
        .name("consecutive-soft-components".to_string())
        .stack_size(64 * 1024 * 1024)
        .spawn(move || -> Result<()> {
            use gammalooprs::graph::parse::IntoGraph;

            const NAME: &str = "local_consecutive_soft_scalar_bubbles";
            const CARD: &str = "local_consecutive_soft_scalar_bubbles.toml";

            let mut subtracted = local_ct_cli(CARD, NAME, false, project_local_4d, explicit_orientation_sum_only)?;
            assert_eq!(
                subtracted.cli_settings.global.generation.uv.orchestrator,
                UVOrchestrator::HedgePoset
            );
            assert_generated_local_integrand(&subtracted, NAME, NAME, 6)?;

            let classified = classified_local_spinneys(&subtracted, NAME, NAME)?;
            let mut connected = classified
                .iter()
                .filter(|spinney| spinney.n_components == 1)
                .map(|spinney| {
                    (
                        spinney.edge_ids.clone(),
                        spinney.dod,
                        spinney.scheme,
                    )
                })
                .collect::<Vec<_>>();
            connected.sort_by(|left, right| left.0.cmp(&right.0));
            assert_eq!(
                connected,
                vec![
                    (vec![0, 1], 2, ApproximationType::IR),
                    (vec![3, 4], 2, ApproximationType::IR),
                ],
                "the connected carrier must contain exactly two consecutive degree-two soft bubbles"
            );
            let unions = classified
                .iter()
                .filter(|spinney| spinney.n_components == 2)
                .collect::<Vec<_>>();
            assert!(
                unions.len() == 1
                    && unions[0].edge_ids == [0, 1, 3, 4]
                    && unions[0].scheme == ApproximationType::MUV,
                "the two consecutive bubbles must have one neutral two-component union: {unions:?}"
            );
            assert_eq!(
                classified
                    .iter()
                    .filter(|spinney| !spinney.edge_ids.is_empty())
                    .count(),
                3,
                "the fixture must have only its two connected components and their compatible union: {classified:?}"
            );

            let process_id = subtracted
                .state
                .resolve_process_ref(Some(&ProcessRef::Unqualified(NAME.to_string())))?;
            let export = subtracted.state.process_list.processes[process_id].export_uv_forest_graph(
NAME,
0,
                &gammalooprs::uv::export::UVForestExportSettings {
                    computed: true,
                },
            )?;
            assert_unique_uv_forest_term_identities(&export, "consecutive soft forest");

            let mut atoms_by_identity = BTreeMap::new();
            for term in &export.node_terms {
                let graph: gammalooprs::graph::Graph =
                    term.dot.as_str().into_graph(&subtracted.state.model)?;
                let identity = UVForestTermIdentity::from(term);
                assert!(
                    atoms_by_identity
                        .insert(identity.clone(), graph.global_prefactor.num)
                        .is_none(),
                    "duplicate consecutive-soft forest identity {identity:?}"
                );
            }
            let forest_indices = atoms_by_identity
                .keys()
                .map(|identity| identity.forest_index)
                .collect::<BTreeSet<_>>();
            assert_eq!(
                forest_indices,
                [0].into_iter().collect::<BTreeSet<_>>(),
                "the consecutive fixture must export one HedgePoset wood"
            );
            let node_keys = atoms_by_identity
                .keys()
                .map(|identity| identity.node_key.clone())
                .collect::<BTreeSet<_>>();
            assert!(
                node_keys.len() == 4
                    && node_keys.contains("∅")
                    && node_keys.iter().filter(|key| key.contains(',')).count() == 1,
                "the consecutive forest must contain root, two components, and their union: {node_keys:?}"
            );
            for node_key in &node_keys {
                assert!(
                    atoms_by_identity.iter().any(|(identity, atom)| {
                        identity.node_key == *node_key && !atom.is_zero()
                    }),
                    "consecutive forest node {node_key} has no nonzero final atom"
                );
            }

            let subtracted_uv = cli_uv_profile_pass_fail(&mut subtracted, NAME, NAME)?;
            assert!(
subtracted_uv.total > 0 && subtracted_uv.resolved > 0 && subtracted_uv.failed == 0
 && (explicit_orientation_sum_only || (subtracted_uv.orientation_total > 0 && subtracted_uv.orientation_resolved > 0 && subtracted_uv.orientation_failed == 0)),
                "consecutive soft components must pass a non-vacuous CLI UV profile: {subtracted_uv:?}"
            );

            let mut bare = local_ct_cli(CARD, "local_consecutive_soft_scalar_bubbles_bare", true, project_local_4d, explicit_orientation_sum_only)?;
            let bare_classified = classified_local_spinneys(&bare, NAME, NAME)?;
            assert!(
                bare_classified.iter().all(|component| {
                    component.n_components == 0
                        && component.edge_ids.is_empty()
                        && component.scheme == ApproximationType::MUV
                        && component.dod == 0
                }),
                "the bare consecutive fixture must not retain counterterm components: {bare_classified:?}"
            );
            assert_matched_local_routes(&subtracted, &bare, NAME, NAME)?;
            let bare_uv = cli_uv_profile_pass_fail(&mut bare, NAME, NAME)?;
            assert!(
bare_uv.total > 0 && bare_uv.resolved > 0 && bare_uv.dod_failures > 0
 && (explicit_orientation_sum_only || (bare_uv.orientation_total > 0 && bare_uv.orientation_resolved > 0 && bare_uv.orientation_dod_failures > 0)),
                "bare consecutive fixture must have a resolved per-orientation fit above the UV bound: {bare_uv:?}"
            );
            assert_eq!(
                subtracted_uv.profile_identity, bare_uv.profile_identity,
                "subtracted and bare consecutive controls must profile the same UV LMBs, subsets, and orientations",
            );

            clean_test(&subtracted.cli_settings.state.folder);
            clean_test(&bare.cli_settings.state.folder);
            Ok(())
        })?;
        test.join()
            .map_err(|_| eyre::eyre!("consecutive soft-components test thread panicked"))??;
    }
    Ok(())
}

#[test]
#[serial_test::serial]
fn sunrise_pole_part_matches_muv_inspect() -> Result<()> {
    let process = "sunrise_scalar_1";
    let integrand_name = "scalar_sunrise";
    let setup = |test_name: &str, scheme| -> Result<CLIState> {
        let mut cli = get_test_cli(
            Some("uv/sunrise_scalar_1.toml".into()),
            get_tests_workspace_path().join(test_name),
            Some(test_name.to_string()),
            true,
        )?;
        let uv = &mut cli.cli_settings.global.generation.uv;
        uv.orchestrator = UVOrchestrator::Compare;
        let prescription = &mut uv.renormalization_prescription;
        prescription.log_divergent = scheme;
        prescription.massive_power_divergent = scheme;
        prescription.massless_power_divergent = scheme;
        cli.run_command("run generate")?;
        Ok(cli)
    };
    let mut muv = setup("sunrise_muv_inspect", ApproximationType::MUV)?;
    let point = deterministic_uv_momentum_points(&mut muv, process, integrand_name)?.remove(0);
    let mut pole_part = setup("sunrise_pole_part_inspect", ApproximationType::PolePart)?;

    let muv_value = evaluate_momentum_sample(
        &mut muv,
        process,
        integrand_name,
        &point.point,
        None,
        f64::EPSILON,
    )?
    .value;
    let pole_part_value = evaluate_momentum_sample(
        &mut pole_part,
        process,
        integrand_name,
        &point.point,
        None,
        f64::EPSILON,
    )?
    .value;

    clean_test(&muv.cli_settings.state.folder);
    clean_test(&pole_part.cli_settings.state.folder);

    let delta = (muv_value.re - pole_part_value.re).hypot(muv_value.im - pole_part_value.im);
    let scale = muv_value
        .re
        .hypot(muv_value.im)
        .max(pole_part_value.re.hypot(pole_part_value.im))
        .max(f64::MIN_POSITIVE);
    assert!(
        delta / scale <= 1.0e-10,
        "MUV and PolePart inspect values differ: MUV={muv_value:?}, PolePart={pole_part_value:?}, relative delta={}",
        delta / scale
    );
    Ok(())
}

mod slow {
    use super::*;

    fn run_double_triangle_soft_ir_case(
        card: &str,
        name: &str,
        external_pdgs: &[isize],
        internal_pdgs: &[isize],
        soft_dod: i32,
        project_local_4d: bool,
        explicit_orientation_sum_only: bool,
    ) -> Result<()> {
        use gammalooprs::cff::orientations::GraphOrientation;
        use gammalooprs::graph::FeynmanGraph;

        const POINT: [f64; 9] = [0.31, -0.27, 0.43, -0.22, 0.37, 0.19, 0.41, 0.16, -0.33];
        const EXPECTED_ORIENTATION: &str = "orientation_delta(0,0,-1,-1,-1,1,-1,-1,1,-1)";

        // Enumerating the parent CFF is cheap compared with constructing every
        // three-loop UV forest. Keep this probe ahead of local CT generation so
        // that the fixture's generation-time orientation can be audited against
        // the graph rather than copied from a generated evaluator.
        let probe_name = format!(
            "{name}_orientation_probe_explicit_{explicit_orientation_sum_only}_project_local_4d_{project_local_4d}"
        );
        let mut probe = get_test_cli(
            Some(card.into()),
            get_tests_workspace_path().join(&probe_name),
            Some(probe_name.clone()),
            true,
        )?;
        let generation = &mut probe.cli_settings.global.generation;
        generation.explicit_orientation_sum_only = explicit_orientation_sum_only;
        generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d;
        if explicit_orientation_sum_only {
            generation.orientation_pattern = Default::default();
        }
        let setup_commands = probe
            .run_history
            .command_blocks
            .iter()
            .find(|block| block.name == "generate")
            .ok_or_else(|| eyre::eyre!("{card} has no generate command block"))?
            .commands
            .iter()
            .filter_map(|command| command.raw_string.clone())
            .take_while(|command| !command.trim_start().starts_with("generate "))
            .collect::<Vec<_>>();
        for command in setup_commands {
            probe.run_command(&command)?;
        }
        let probe_process_id = probe
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified(name.to_string())))?;
        let mut parent_graph =
            match &probe.state.process_list.processes[probe_process_id].collection {
                ProcessCollection::Amplitudes(amplitudes) => amplitudes
                    .get(name)
                    .unwrap_or_else(|| panic!("missing imported amplitude {name}"))
                    .graphs[0]
                    .graph
                    .clone(),
                ProcessCollection::CrossSections(_) => {
                    return Err(eyre::eyre!("{name} must be an amplitude"));
                }
            };
        let parent_tree_edges = parent_graph
            .iter_edges_of(&parent_graph.tree_edges)
            .filter_map(|(pair, edge, _)| pair.is_paired().then_some(edge))
            .collect::<Vec<_>>();
        let options = parent_graph
            .production_cff_3d_expression_options(&probe.cli_settings.global.generation)?;
        let parent_cff = parent_graph.generate_3d_expression_for_integrand(
            &parent_tree_edges,
            &parent_graph.get_esurface_canonization(&parent_graph.loop_momentum_basis),
            &options,
            None,
        )?;
        let parent_orientations = parent_cff
            .expression
            .orientations
            .iter()
            .map(|orientation| &orientation.data.orientation)
            .collect::<Vec<_>>();
        let orientation_signatures = parent_orientations
            .iter()
            .filter(|orientation| !orientation.orientation_thetas().is_one())
            .map(|orientation| orientation.orientation_delta().to_string())
            .collect::<BTreeSet<_>>();
        assert!(
            orientation_signatures.len() > 1,
            "{name} must have a nontrivial set of valid nonzero parent-CFF orientations: {orientation_signatures:#?}"
        );
        assert_eq!(
            orientation_signatures.first().map(String::as_str),
            Some(EXPECTED_ORIENTATION),
            "the fixture must select the deterministic first nonzero parent-CFF orientation"
        );
        // Keep source-map multiplicities: equal direction signatures do not
        // identify maps with different numerator or mass sampling data.
        let mut source_orientations = parent_orientations
            .iter()
            .map(|orientation| orientation.orientation_delta().to_string())
            .collect::<Vec<_>>();
        source_orientations.sort();
        let expected_orientations = if explicit_orientation_sum_only {
            source_orientations.clone()
        } else {
            vec![EXPECTED_ORIENTATION.to_string()]
        };
        let mut configured_orientations = parent_orientations
            .iter()
            .filter(|&&orientation| {
                probe
                    .cli_settings
                    .global
                    .generation
                    .orientation_pattern
                    .filter(orientation)
            })
            .map(|orientation| orientation.orientation_delta().to_string())
            .collect::<Vec<_>>();
        configured_orientations.sort();
        assert_eq!(
            configured_orientations, expected_orientations,
            "the generation-time orientation pattern must retain exactly the independently audited parent-CFF inventory"
        );
        clean_test(&probe.cli_settings.state.folder);

        let mut subtracted = local_ct_cli(
            card,
            name,
            false,
            project_local_4d,
            explicit_orientation_sum_only,
        )?;
        let expected_identifier = CTIdentifier::new(
            external_pdgs.iter().copied().collect(),
            Some(internal_pdgs.iter().copied().collect()),
        );
        let classified = classified_local_spinneys(&subtracted, name, name)?;
        let selected = classified
            .iter()
            .filter(|spinney| {
                spinney.n_components == 1
                    && spinney.scheme == ApproximationType::IR
                    && spinney.identifier == expected_identifier
            })
            .collect::<Vec<_>>();
        assert!(
            selected.len() == 1 && selected[0].dod == soft_dod,
            "{name} must select exactly its degree-{soft_dod} central soft self-energy: {classified:?}"
        );
        let child_edges = &selected[0].edge_ids;
        assert!(
            classified.iter().any(|spinney| {
                spinney.n_components == 1
                    && spinney.scheme == ApproximationType::IR
                    && spinney.dod > 0
                    && spinney.edge_ids.len() > child_edges.len()
                    && child_edges
                        .iter()
                        .all(|edge| spinney.edge_ids.contains(edge))
            }),
            "{name} must use a complete soft wood for its positive-degree containing component: {classified:?}"
        );

        let process_id = subtracted
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified(name.to_string())))?;
        let integrand = subtracted
            .state
            .process_list
            .get_integrand(process_id, name)?
            .require_generated()?;
        assert_eq!(
            integrand.get_n_dim(),
            9,
            "{name} must retain its three loop integrations"
        );

        let mut bare = local_ct_cli(
            card,
            &format!("{name}_bare"),
            true,
            project_local_4d,
            explicit_orientation_sum_only,
        )?;
        assert_matched_local_routes(&subtracted, &bare, name, name)?;
        for (label, cli) in [("subtracted", &subtracted), ("bare", &bare)] {
            let process_id = cli
                .state
                .resolve_process_ref(Some(&ProcessRef::Unqualified(name.to_string())))?;
            let process = &cli.state.process_list.processes[process_id];
            let amplitude = match &process.collection {
                ProcessCollection::Amplitudes(amplitudes) => amplitudes
                    .get(name)
                    .unwrap_or_else(|| panic!("missing generated amplitude {name}")),
                ProcessCollection::CrossSections(_) => {
                    return Err(eyre::eyre!("{name} must be an amplitude"));
                }
            };
            let mut orientations = amplitude.graphs[0]
                .derived_data
                .cff_expression
                .as_ref()
                .expect("the generated amplitude must retain its CFF orientations")
                .expression
                .orientations
                .iter()
                .map(|orientation| orientation.data.orientation.orientation_delta().to_string())
                .collect::<Vec<_>>();
            orientations.sort();
            assert_eq!(
                orientations, source_orientations,
                "{name} {label} stored CFF must retain the full independently audited source inventory"
            );
            // Generation filters the evaluator's channels, while its stored
            // source CFF remains complete for local counterterm construction.
            let info = cli.state.get_integrand_info(
                Some(&ProcessRef::Unqualified(name.to_string())),
                Some(&name.to_string()),
            )?;
            assert!(!info.graph_groups.is_empty());
            for group in &info.graph_groups {
                let mut orientations = group
                    .orientations
                    .iter()
                    .map(|orientation| {
                        format!(
                            "orientation_delta({})",
                            orientation
                                .signature
                                .iter()
                                .map(ToString::to_string)
                                .collect::<Vec<_>>()
                                .join(",")
                        )
                    })
                    .collect::<Vec<_>>();
                orientations.sort();
                assert_eq!(
                    orientations, expected_orientations,
                    "{name} {label} evaluator must retain exactly the audited orientation inventory"
                );
            }
        }

        // Inspect the public, structure-only forest export without asking this
        // 3D fixture to build a second, unused 4D representation solely for
        // provenance metadata. The generated evaluator below is the
        // authoritative atom for the selected orientation representation;
        // focused operator tests independently retain U, S, and US.
        let structural_forest = subtracted.state.process_list.processes[process_id]
            .export_uv_forest_graph(
                name,
                0,
                &gammalooprs::uv::export::UVForestExportSettings { computed: false },
            )?;
        let structural_node_count = structural_forest.forest_dot.matches("foata=").count();
        assert!(
            structural_forest.node_terms.is_empty() && structural_node_count >= 4,
            "{name} must expose distinct root, child/parent, and nested forest families in its structure-only export; nodes={structural_node_count}"
        );

        let subtracted_floor = required_stability_accuracy_floor(&mut subtracted, name, name)?;
        let subtracted_value =
            evaluate_momentum_sample(&mut subtracted, name, name, &POINT, None, subtracted_floor)?;
        let bare_classified = classified_local_spinneys(&bare, name, name)?;
        assert!(
            bare_classified.iter().all(|component| {
                component.n_components == 0
                    && component.edge_ids.is_empty()
                    && component.scheme == ApproximationType::MUV
                    && component.dod == 0
            }),
            "the bare {name} fixture must not retain a counterterm component: {bare_classified:?}"
        );
        let bare_floor = required_stability_accuracy_floor(&mut bare, name, name)?;
        let bare_value = evaluate_momentum_sample(&mut bare, name, name, &POINT, None, bare_floor)?;

        let delta = (subtracted_value.value.re - bare_value.value.re)
            .hypot(subtracted_value.value.im - bare_value.value.im);
        let scale = subtracted_value
            .value
            .re
            .hypot(subtracted_value.value.im)
            .max(bare_value.value.re.hypot(bare_value.value.im))
            .max(f64::MIN_POSITIVE);
        let relative_delta = delta / scale;
        let relative_accuracy = subtracted_value
            .relative_accuracy
            .max(bare_value.relative_accuracy)
            .max(f64::EPSILON);
        assert!(
            subtracted_value.value.re.is_finite()
                && subtracted_value.value.im.is_finite()
                && bare_value.value.re.is_finite()
                && bare_value.value.im.is_finite()
                && INSPECT_DEPENDENCE_ACCURACY_FACTOR * relative_accuracy < 1.0,
            "{name} inspect values are not materially resolved: subtracted={:?}, bare={:?}, accuracy={relative_accuracy:.3e}",
            subtracted_value.value,
            bare_value.value,
        );
        assert!(
            relative_delta >= INSPECT_DEPENDENCE_ACCURACY_FACTOR * relative_accuracy,
            "{name} local counterterm is not visible above inspect accuracy: subtracted={:?}, bare={:?}, relative delta={relative_delta:.3e}, accuracy={relative_accuracy:.3e}",
            subtracted_value.value,
            bare_value.value,
        );

        let subtracted_uv = cli_uv_profile_pass_fail(&mut subtracted, name, name)?;
        assert!(
            subtracted_uv.total > 0
                && subtracted_uv.resolved > 0
                && subtracted_uv.failed == 0
                && (explicit_orientation_sum_only
                    || (subtracted_uv.orientation_total > 0
                        && subtracted_uv.orientation_resolved > 0
                        && subtracted_uv.orientation_failed == 0)),
            "{name} local counterterm must pass a non-vacuous CLI UV profile: {subtracted_uv:?}"
        );
        let bare_uv = cli_uv_profile_pass_fail(&mut bare, name, name)?;
        assert!(
            bare_uv.total > 0
                && bare_uv.resolved > 0
                && bare_uv.dod_failures > 0
                && (explicit_orientation_sum_only
                    || (bare_uv.orientation_total > 0
                        && bare_uv.orientation_resolved > 0
                        && bare_uv.orientation_dod_failures > 0)),
            "bare {name} must have a resolved per-orientation fit above the UV bound: {bare_uv:?}"
        );
        assert_eq!(
            subtracted_uv.profile_identity, bare_uv.profile_identity,
            "subtracted and bare {name} controls must profile the same UV LMBs, subsets, and orientations",
        );

        clean_test(&subtracted.cli_settings.state.folder);
        clean_test(&bare.cli_settings.state.folder);
        Ok(())
    }

    #[test]
    #[serial_test::serial]
    fn paper_figure_b1_massless_bubble_uses_a_complete_soft_wood() -> Result<()> {
        for (explicit_orientation_sum_only, project_local_4d) in
            [(false, false), (true, false), (true, true)]
        {
            println!(
                "explicit_orientation_sum_only={explicit_orientation_sum_only}, project_local_4d={project_local_4d}"
            );
            let test = std::thread::Builder::new()
                .name("local-ir-double-triangle".to_string())
                .stack_size(64 * 1024 * 1024)
                .spawn(move || {
                    run_double_triangle_soft_ir_case(
                        "paper_figure_b1_double_triangle_soft_ir.toml",
                        "paper_figure_b1_double_triangle_soft_ir",
                        &[-21, 21],
                        &[1],
                        2,
                        project_local_4d,
                        explicit_orientation_sum_only,
                    )
                })?;
            test.join()
                .map_err(|_| eyre::eyre!("local IR double-triangle test thread panicked"))??;
        }
        Ok(())
    }

    #[test]
    #[serial_test::serial]
    fn paper_figure_b1_double_triangle_is_a_non_vacuous_uv_fixture() -> Result<()> {
        for (explicit_orientation_sum_only, project_local_4d) in
            [(false, false), (true, false), (true, true)]
        {
            println!(
                "explicit_orientation_sum_only={explicit_orientation_sum_only}, project_local_4d={project_local_4d}"
            );
            const NAME: &str = "paper_figure_b1_double_triangle";
            let card = "paper_figure_b1_double_triangle.toml";
            let mut subtracted = local_ct_cli(
                card,
                NAME,
                false,
                project_local_4d,
                explicit_orientation_sum_only,
            )?;
            assert_generated_local_integrand(&subtracted, NAME, NAME, 6)?;
            assert!(
                classified_local_spinneys(&subtracted, NAME, NAME)?
                    .into_iter()
                    .any(|spinney| spinney.n_components == 1),
                "Figure-13 B.1 must expose at least one connected UV component"
            );
            let subtracted_uv = cli_uv_profile_pass_fail(&mut subtracted, NAME, NAME)?;
            assert!(
                subtracted_uv.total > 0
                    && subtracted_uv.resolved > 0
                    && subtracted_uv.failed == 0
                    && (explicit_orientation_sum_only
                        || (subtracted_uv.orientation_total > 0
                            && subtracted_uv.orientation_resolved > 0
                            && subtracted_uv.orientation_failed == 0)),
                "Figure-13 B.1 must pass a non-vacuous CLI UV profile: {subtracted_uv:?}"
            );

            let mut bare = local_ct_cli(
                card,
                "paper_figure_b1_double_triangle_bare",
                true,
                project_local_4d,
                explicit_orientation_sum_only,
            )?;
            let bare_uv = cli_uv_profile_pass_fail(&mut bare, NAME, NAME)?;
            assert!(
                bare_uv.total > 0
                    && bare_uv.resolved > 0
                    && bare_uv.dod_failures > 0
                    && (explicit_orientation_sum_only
                        || (bare_uv.orientation_total > 0
                            && bare_uv.orientation_resolved > 0
                            && bare_uv.orientation_dod_failures > 0)),
                "bare Figure-13 B.1 must have a resolved per-orientation fit above the UV bound: {bare_uv:?}"
            );
            assert_eq!(
                subtracted_uv.profile_identity, bare_uv.profile_identity,
                "subtracted and bare Figure-13 B.1 controls must profile the same UV LMBs, subsets, and orientations",
            );
            clean_test(&subtracted.cli_settings.state.folder);
            clean_test(&bare.cli_settings.state.folder);
        }
        Ok(())
    }

    #[test]
    #[serial_test::serial]
    fn local_nested_soft_fixture_generates() -> Result<()> {
        for (explicit_orientation_sum_only, project_local_4d) in
            [(false, false), (true, false), (true, true)]
        {
            println!(
                "explicit_orientation_sum_only={explicit_orientation_sum_only}, project_local_4d={project_local_4d}"
            );
            let test = std::thread::Builder::new()
            .name("local-nested-soft".to_string())
            .stack_size(128 * 1024 * 1024)
            .spawn(move || -> Result<()> {
        let mut nested_soft = local_ct_cli(
            "local_nested_soft_top_self_energy.toml",
            "local_nested_soft_top_self_energy",
            false, project_local_4d, explicit_orientation_sum_only)?;
        assert_eq!(
            nested_soft.cli_settings.global.generation.uv.orchestrator,
            UVOrchestrator::HedgePoset
        );
        assert_generated_local_integrand(
            &nested_soft,
            "local_nested_soft_top_self_energy",
            "local_nested_soft_top_self_energy",
            6,
        )?;
        assert_selected_local_identifier(
            &nested_soft,
            "local_nested_soft_top_self_energy",
            "local_nested_soft_top_self_energy",
            ApproximationType::IR,
            [-6, 6],
            [6, 21],
            2,
        )?;
        let classified = classified_local_spinneys(
            &nested_soft,
            "local_nested_soft_top_self_energy",
            "local_nested_soft_top_self_energy",
        )?;
        let connected = classified
            .iter()
            .filter(|spinney| spinney.n_components == 1)
            .collect::<Vec<_>>();
        let inert = classified
            .iter()
            .filter(|spinney| spinney.n_components == 0)
            .collect::<Vec<_>>();
        assert!(
            connected.len() == 2
                && classified.len() == connected.len() + inert.len()
                && inert.iter().all(|component| {
                    component.edge_ids.is_empty()
                        && component.scheme == ApproximationType::MUV
                        && component.dod == 0
                })
                && connected.iter().all(|component| {
                    component.scheme == ApproximationType::IR && component.dod == 1
                })
                && connected
                    .iter()
                    .any(|component| component.edge_ids == [3, 4])
                && connected
                    .iter()
                    .any(|component| component.edge_ids == [2, 3, 4, 5, 6]),
            "apart from its inert empty root, the nested-soft fixture must contain exactly its d=1 H child and d=1 H parent: {classified:?}",
        );
        let process_id = nested_soft
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified(
                "local_nested_soft_top_self_energy".to_string(),
            )))?;
        let structural_forest = nested_soft.state.process_list.processes[process_id].export_uv_forest_graph(
"local_nested_soft_top_self_energy",
0,
            &gammalooprs::uv::export::UVForestExportSettings { computed: false },
        )?;
        assert!(
            structural_forest.node_terms.is_empty(),
            "a structure-only nested H/H export must not construct local atoms",
        );
        assert_eq!(
            structural_forest.forest_dot.matches("foata=").count(),
            4,
            "the nested H/H forest must contain exactly the root, two singleton operations, and their child-to-parent chain",
        );
        let subtracted_uv = cli_uv_profile_pass_fail(
            &mut nested_soft,
            "local_nested_soft_top_self_energy",
            "local_nested_soft_top_self_energy",
        )?;
        assert!(
subtracted_uv.total > 0 && subtracted_uv.resolved > 0 && subtracted_uv.failed == 0
 && (explicit_orientation_sum_only || (subtracted_uv.orientation_total > 0 && subtracted_uv.orientation_resolved > 0 && subtracted_uv.orientation_failed == 0)),
            "nested soft fixture must pass a non-vacuous CLI UV profile: {subtracted_uv:?}"
        );

        let mut bare = local_ct_cli(
            "local_nested_soft_top_self_energy.toml",
            "local_nested_soft_top_self_energy_bare",
            true, project_local_4d, explicit_orientation_sum_only)?;
        let bare_uv = cli_uv_profile_pass_fail(
            &mut bare,
            "local_nested_soft_top_self_energy",
            "local_nested_soft_top_self_energy",
        )?;
        assert!(
bare_uv.total > 0 && bare_uv.resolved > 0 && bare_uv.dod_failures > 0
 && (explicit_orientation_sum_only || (bare_uv.orientation_total > 0 && bare_uv.orientation_resolved > 0 && bare_uv.orientation_dod_failures > 0)),
            "bare nested fixture must have a resolved per-orientation fit above the UV bound: {bare_uv:?}"
        );
        assert_eq!(
            subtracted_uv.profile_identity, bare_uv.profile_identity,
            "subtracted and bare nested-soft controls must profile the same UV LMBs, subsets, and orientations",
        );
        clean_test(&nested_soft.cli_settings.state.folder);
        clean_test(&bare.cli_settings.state.folder);
        Ok(())
            })?;
            test.join()
                .map_err(|_| eyre::eyre!("local nested soft test thread panicked"))??;
        }
        Ok(())
    }

    #[test]
    #[serial_test::serial]
    fn paper_appendix_b1_ir_specialization_generates_and_profiles() -> Result<()> {
        for (explicit_orientation_sum_only, project_local_4d) in
            [(false, false), (true, false), (true, true)]
        {
            println!(
                "explicit_orientation_sum_only={explicit_orientation_sum_only}, project_local_4d={project_local_4d}"
            );
            let test = std::thread::Builder::new()
            .name("paper-appendix-b1-soft-ir".to_string())
            .stack_size(128 * 1024 * 1024)
            .spawn(move || -> Result<()> {
                const NAME: &str = "paper_appendix_b1_nested_gluon_self_energy";
                let card = "paper_appendix_b1_nested_gluon_self_energy_ir.toml";
                let mut subtracted = local_ct_cli(card, NAME, false, project_local_4d, explicit_orientation_sum_only)?;
        assert_generated_local_integrand(&subtracted, NAME, "soft_ir", 6)?;
        let expected = CTIdentifier::new(
            [-21, 21].into_iter().collect(),
            Some([6, 21].into_iter().collect()),
        );
        let appendix_components = classified_local_spinneys(&subtracted, NAME, "soft_ir")?;
        let connected_components = appendix_components
            .iter()
            .filter(|spinney| spinney.n_components == 1)
            .collect::<Vec<_>>();
        let outer_soft_components = connected_components
            .iter()
            .copied()
            .filter(|spinney| {
                spinney.scheme == ApproximationType::IR
                    && spinney.identifier == expected
            })
            .collect::<Vec<_>>();
                assert!(
                    connected_components.len() == 3
                        && outer_soft_components.len() == 1
                        && outer_soft_components[0].edge_ids == [2, 3, 4, 5, 6]
                        && outer_soft_components[0].dod == 2
                        && connected_components.iter().any(|component| {
                            component.edge_ids == [3, 4]
                                && component.scheme == ApproximationType::MUV
                                && component.dod == 1
                        })
                        && connected_components.iter().any(|component| {
                            component.edge_ids == [2, 3, 5, 6]
                                && component.scheme == ApproximationType::MUV
                                && component.dod == 0
                        }),
                    "Appendix-B.1 must contain exactly gamma_1=U_1, gamma_2=U_0, and the outer H_2 component: {appendix_components:?}"
                );

                // Keep the two nested forest pairs independently observable. With the
                // outer component switched off, this is F_empty + F_child =
                // (1-U_child)I. It must already improve every child-hard ray; adding
                // the outer pair must not reintroduce that subdivergence.
                let mut child_only = get_test_cli(
                    Some(card.into()),
                    get_tests_workspace_path().join(format!("paper_appendix_b1_child_only_explicit_{explicit_orientation_sum_only}_project_local_4d_{project_local_4d}")),
                    Some(format!("paper_appendix_b1_child_only_explicit_{explicit_orientation_sum_only}_project_local_4d_{project_local_4d}")),
                    true,
                )?;
            let generation = &mut child_only.cli_settings.global.generation;
            generation.explicit_orientation_sum_only = explicit_orientation_sum_only;
            generation.uv.final_integrand = FinalIntegrandDimension::ThreeD;
            generation.uv.local_uv_cts_from_expanded_4d_integrands = project_local_4d;
                let outer_rule = child_only
                    .cli_settings
                    .global
                    .generation
                    .uv
                    .renormalization_prescription
                    .overrides
                    .iter_mut()
                    .find(|rule| rule.ct_identifier == expected)
                    .ok_or_else(|| {
                        eyre::eyre!("Appendix-B.1 child-only control has no outer CT rule")
                    })?;
                assert_eq!(outer_rule.prescription, ApproximationType::IR);
                outer_rule.prescription = ApproximationType::Unsubtracted;
                let gamma_two = CTIdentifier::new(
                    [-21, 21].into_iter().collect(),
                    Some([6].into_iter().collect()),
                );
                child_only
                    .cli_settings
                    .global
                    .generation
                    .uv
                    .renormalization_prescription
                    .overrides
                    .push(CTRenormalizationRule::new(
                        gamma_two,
                        ApproximationType::Unsubtracted,
                    ));
                child_only.run_command("run generate")?;
                let child_only_components = classified_local_spinneys(
                    &child_only,
                    NAME,
                    "soft_ir",
                )?
                .into_iter()
                .filter(|spinney| spinney.n_components == 1)
                .map(|spinney| (spinney.edge_ids, spinney.dod, spinney.scheme))
                .collect::<Vec<_>>();
                assert_eq!(
                    child_only_components,
                    vec![(vec![3, 4], 1, ApproximationType::MUV)],
                    "the Appendix-B.1 child-only control must retain exactly the ordinary d=1 gamma_1 component",
                );
                assert_matched_local_routes(&subtracted, &child_only, NAME, "soft_ir")?;
                let child_only_uv =
                    cli_uv_profile_pass_fail(&mut child_only, NAME, "soft_ir")?;
                assert!(
child_only_uv.total > 0 && child_only_uv.resolved > 0
 && (explicit_orientation_sum_only || (child_only_uv.orientation_total > 0 && child_only_uv.orientation_resolved > 0)),
                    "Appendix-B.1 child-only control produced no resolved summed and per-orientation UV limits: {child_only_uv:?}"
                );
                let child_only_report: Value = serde_json::from_str(&fs::read_to_string(
                    child_only
                        .cli_settings
                        .state
                        .folder
                        .join("local_ct_uv_profile/uv_profile.json"),
                )?)?;
                let child_hard_rays = child_only_report["graphs"]
                    .as_array()
                    .into_iter()
                    .flatten()
                    .flat_map(|graph| graph["lmbs"].as_array().into_iter().flatten())
                    .flat_map(|lmb| lmb["subsets"].as_array().into_iter().flatten())
                    .filter(|subset| subset["initial_dod"] == 1)
                    .map(|subset| {
                        let edges = |field: &str| {
                            subset[field]
                                .as_array()
                                .ok_or_else(|| {
                                    eyre::eyre!(
                                        "Appendix-B.1 child-hard control has no {field} edge list: {subset}"
                                    )
                                })?
                                .iter()
                                .map(|edge| {
                                    edge.as_u64().map(|edge| edge as usize).ok_or_else(|| {
                                        eyre::eyre!(
                                            "Appendix-B.1 child-hard control has a malformed {field} edge: {subset}"
                                        )
                                    })
                                })
                                .collect::<Result<Vec<_>>>()
                        };
                        let slope = subset["analysis"]["inspect_level"]["result"]["slope"]
                            .as_f64()
                            .ok_or_else(|| {
                                eyre::eyre!(
                                    "Appendix-B.1 child-hard control has no finite slope: {subset}"
                                )
                            })?;
                        Ok((edges("free")?, edges("fixed")?, slope))
                    })
                    .collect::<Result<Vec<_>>>()?;
                let observed_child_hard_pairs = child_hard_rays
                    .iter()
                    .map(|(free, fixed, _)| (free.clone(), fixed.clone()))
                    .collect::<BTreeSet<_>>();
                let expected_child_hard_pairs = [3, 4]
                    .into_iter()
                    .flat_map(|free| {
                        [2, 5, 6]
                            .into_iter()
                            .map(move |fixed| (vec![free], vec![fixed]))
                    })
                    .collect::<BTreeSet<_>>();
                assert_eq!(
                    child_hard_rays.len(),
                    expected_child_hard_pairs.len(),
                    "the Appendix-B.1 child-only profile must contain exactly six d=1 rays: {child_hard_rays:?}",
                );
                assert_eq!(
                    observed_child_hard_pairs, expected_child_hard_pairs,
                    "the Appendix-B.1 d=1 selector must resolve exactly the six gamma_1-hard rays with free e3/e4 and fixed outer representative e2/e5/e6",
                );
                assert!(
                    child_hard_rays
                        .iter()
                        .all(|(_, _, slope)| *slope < -0.9),
                    "the isolated (1-U_child) Appendix-B.1 forest pair did not remove its child subdivergence: {child_hard_rays:?}"
                );

                let process_id = subtracted
            .state
            .resolve_process_ref(Some(&ProcessRef::Unqualified(NAME.to_string())))?;
        let structural_forest = subtracted.state.process_list.processes[process_id].export_uv_forest_graph(
"soft_ir",
0,
            &gammalooprs::uv::export::UVForestExportSettings { computed: false },
        )?;
        assert!(
            structural_forest.node_terms.is_empty(),
            "a structure-only Appendix-B.1 export must not construct local atoms",
        );
        assert_eq!(
            structural_forest.forest_dot.matches("foata=").count(),
            6,
            "Appendix-B.1 must contain the root, three singleton operations, and the two allowed child-to-outer chains",
        );

        let computed_forest = subtracted.state.process_list.processes[process_id].export_uv_forest_graph(
"soft_ir",
0,
            &gammalooprs::uv::export::UVForestExportSettings { computed: true },
        )?;
        assert_unique_uv_forest_term_identities(
            &computed_forest,
            "Appendix-B.1 direct-3D forest",
        );
        let exported_nodes = computed_forest
            .node_terms
            .iter()
            .map(|term| (term.forest_index, term.node_key.clone()))
            .collect::<BTreeSet<_>>();
        assert_eq!(
            exported_nodes.len(),
            6,
            "the computed Appendix-B.1 forest must retain every structural node"
        );
        let mut saw_nested_path = false;
        let mut positive_degree_ir_steps = 0;
        for term in &computed_forest.node_terms {
            let provenance = term
                .forest_provenance()?
                .ok_or_else(|| {
                    eyre::eyre!(
                        "Appendix-B.1 computed term {} has no direct-3D provenance",
                        term.file_name()
                    )
                })?;
            let node = &provenance["node"];
            assert!(
                node.get("components").is_none()
                    && node.get("local_4d_branches").is_none(),
                "Appendix-B.1 direct-3D provenance must not expose legacy 4D fields: {node}"
            );
            let parent_keys = node["parent_keys"].as_array().unwrap_or_else(|| {
                panic!("Appendix-B.1 term has malformed parent provenance: {node}")
            });
            if term.node_key == "∅" {
                assert!(parent_keys.is_empty());
            } else {
                assert!(
                    !parent_keys.is_empty(),
                    "Appendix-B.1 non-root node {} has no forest parent",
                    term.node_key
                );
            }
            for parent_key in parent_keys {
                let parent_key = parent_key.as_str().unwrap_or_else(|| {
                    panic!("Appendix-B.1 has a non-string parent key: {node}")
                });
                assert!(
                    exported_nodes.contains(&(term.forest_index, parent_key.to_string())),
                    "Appendix-B.1 node {} refers to missing parent {parent_key:?}",
                    term.node_key
                );
            }

            let local = &node["local"];
            assert_eq!(
                local["representation"], "3d",
                "Appendix-B.1 computed CFF terms require native direct-3D provenance"
            );
            assert!(
                local.get("components").is_none()
                    && local.get("branches").is_none()
                    && local.get("local_4d_branches").is_none(),
                "Appendix-B.1 direct-3D payload must not force a 4D branch construction: {local}"
            );
            // Only direct 3D replays record projection paths. Node identities
            // and physical profiles remain checked in both routes.
            if !project_local_4d {
                let paths = local["projection_paths"].as_array().unwrap_or_else(|| {
                    panic!("Appendix-B.1 term has malformed direct-3D paths: {local}")
                });
                assert!(
                    !paths.is_empty(),
                    "Appendix-B.1 must retain at least the identity projection path"
                );
                if term.node_key == "∅" {
                    assert!(
                        paths.iter().all(|path| path["steps"]
                            .as_array()
                            .is_some_and(Vec::is_empty)),
                        "the Appendix-B.1 root must carry only its empty identity path: {paths:?}"
                    );
                } else {
                    assert!(
                        paths.iter().all(|path| path["steps"]
                            .as_array()
                            .is_some_and(|steps| !steps.is_empty())),
                        "each Appendix-B.1 non-root direct-3D path must contain an actual projection: {paths:?}"
                    );
                }

                for path in paths {
                    let steps = path["steps"].as_array().unwrap_or_else(|| {
                        panic!("Appendix-B.1 has a malformed projection path: {path}")
                    });
                    if steps.len() > 1 {
                        saw_nested_path = true;
                        assert!(
                            steps.windows(2).all(|pair| {
                                let child = pair[0]["topo_order"].as_u64().unwrap_or_else(|| {
                                    panic!("Appendix-B.1 child step has no topology order: {}", pair[0])
                                });
                                let parent = pair[1]["topo_order"].as_u64().unwrap_or_else(|| {
                                    panic!("Appendix-B.1 parent step has no topology order: {}", pair[1])
                                });
                                child < parent
                            }),
                            "Appendix-B.1 nested projection paths must be ordered child-to-parent: {path}"
                        );
                    }
                    for step in steps {
                        assert!(
                            step["current_component"].is_string()
                                && step["given_component"].is_string()
                                && step["active_subgraph"].is_string()
                                && step["rescaled_subgraph"].is_string()
                                && step["canonical_route_loop_edges"]
                                    .as_array()
                                    .is_some_and(|edges| !edges.is_empty())
                                && step["canonical_route_external_edges"].is_array()
                                && step["canonical_route_signatures"]
                                    .as_array()
                                    .is_some_and(|signatures| !signatures.is_empty())
                                && step["route_loop_edges"]
                                    .as_array()
                                    .is_some_and(|edges| !edges.is_empty())
                                && step["route_external_edges"].is_array()
                                && step["route_signatures"]
                                    .as_array()
                                    .is_some_and(|signatures| !signatures.is_empty()),
                            "Appendix-B.1 projection is missing its canonical or actual local LMB: {step}"
                        );
                        if step["scheme"] == "IR"
                            && step["dod"].as_i64().is_some_and(|dod| dod > 0)
                        {
                            positive_degree_ir_steps += 1;
                            let branches = step["conceptual_branches"]
                                .as_array()
                                .unwrap_or_else(|| {
                                    panic!("Appendix-B.1 IR step has malformed branches: {step}")
                                })
                                .iter()
                                .map(|branch| {
                                    (
                                        branch["branch"].as_str().unwrap_or_else(|| {
                                            panic!("Appendix-B.1 IR branch has no name: {branch}")
                                        }),
                                        branch["coefficient"].as_i64().unwrap_or_else(|| {
                                            panic!("Appendix-B.1 IR branch has no sign: {branch}")
                                        }),
                                    )
                                })
                                .collect::<Vec<_>>();
                            assert_eq!(step["dod"].as_i64(), Some(2));
                            assert_eq!(step["materialization"], "factorized_soft");
                            assert_eq!(branches, vec![("U", 1), ("S", 1), ("US", -1)]);
                        }
                    }
                }
            }
        }
        if !project_local_4d {
            assert!(
                saw_nested_path,
                "Appendix-B.1 computed forest must expose a nested child-to-parent projection path"
            );
            assert!(
                positive_degree_ir_steps > 0,
                "Appendix-B.1 computed forest must expose its outer factorized H_2 projection"
            );

        }

        let subtracted_uv = cli_uv_profile_pass_fail(&mut subtracted, NAME, "soft_ir")?;
        assert!(
subtracted_uv.total > 0 && subtracted_uv.resolved > 0 && subtracted_uv.failed == 0
 && (explicit_orientation_sum_only || (subtracted_uv.orientation_total > 0 && subtracted_uv.orientation_resolved > 0 && subtracted_uv.orientation_failed == 0)),
            "the Appendix-B.1 IR specialization must pass a non-vacuous CLI UV profile: {subtracted_uv:?}"
        );
        assert_eq!(
            subtracted_uv.profile_identity, child_only_uv.profile_identity,
            "the complete and child-only Appendix-B.1 forests must profile the same UV LMBs, subsets, and orientations",
        );

        let mut bare = local_ct_cli(card, "paper_appendix_b1_bare", true, project_local_4d, explicit_orientation_sum_only)?;
        let bare_uv = cli_uv_profile_pass_fail(&mut bare, NAME, "soft_ir")?;
        assert!(
bare_uv.total > 0 && bare_uv.resolved > 0 && bare_uv.dod_failures > 0
 && (explicit_orientation_sum_only || (bare_uv.orientation_total > 0 && bare_uv.orientation_resolved > 0 && bare_uv.orientation_dod_failures > 0)),
            "the bare Appendix-B.1 topology must have a resolved per-orientation fit above the UV bound: {bare_uv:?}"
        );
        assert_eq!(
            subtracted_uv.profile_identity, bare_uv.profile_identity,
            "subtracted and bare Appendix-B.1 controls must profile the same UV LMBs, subsets, and orientations",
        );
        clean_test(&subtracted.cli_settings.state.folder);
        clean_test(&child_only.cli_settings.state.folder);
        clean_test(&bare.cli_settings.state.folder);
        Ok(())
    })?;
            test.join()
                .map_err(|_| eyre::eyre!("Appendix-B.1 soft-IR test thread panicked"))??;
        }
        Ok(())
    }

    #[test]
    #[serial_test::serial]
    fn gamma_star_ddbar_top_bubble_child_only_has_two_power_soft_improvement() -> Result<()> {
        for (explicit_orientation_sum_only, project_local_4d) in
            [(false, false), (true, false), (true, true)]
        {
            println!(
                "explicit_orientation_sum_only={explicit_orientation_sum_only}, project_local_4d={project_local_4d}"
            );
            let test = std::thread::Builder::new()
            .name("gamma-star-ddbar-top-bubble-child-only".to_string())
            .stack_size(128 * 1024 * 1024)
            .spawn(move || -> Result<()> {
                const PROCESS: &str = "gamma_star_ddbar_top_bubble";
                let mut child_uv = local_ct_cli(
                    "gamma_star_ddbar_top_bubble_child_uv_only.toml",
                    "gamma_star_ddbar_top_bubble_child_uv_only",
                    false, project_local_4d, explicit_orientation_sum_only)?;
                let mut child_h = local_ct_cli(
                    "gamma_star_ddbar_top_bubble_child_soft_ir.toml",
                    "gamma_star_ddbar_top_bubble_child_soft_ir",
                    false, project_local_4d, explicit_orientation_sum_only)?;
                let expected = CTIdentifier::new(
                    [-21, 21].into_iter().collect(),
                    Some([6].into_iter().collect()),
                );

                for (cli, integrand, scheme) in [
                    (&child_uv, "child_uv_only", ApproximationType::MUV),
                    (&child_h, "child_soft_ir", ApproximationType::IR),
                ] {
                    assert_generated_local_integrand(cli, PROCESS, integrand, 6)?;
                    let classified = classified_local_spinneys(cli, PROCESS, integrand)?;
                    let selected = classified
                        .iter()
                        .filter(|spinney| spinney.n_components == 1)
                        .collect::<Vec<_>>();
                    assert!(
                        selected.len() == 1
                            && selected[0].scheme == scheme
                            && selected[0].dod == 2
                            && selected[0].edge_ids == [7, 8]
                            && selected[0].identifier == expected
                            && classified.iter().all(|spinney| {
                                spinney.n_components == 0 || std::ptr::eq(spinney, selected[0])
                            }),
                        "{integrand} must select exactly the massive-top child and no overall vertex counterterm: {classified:?}",
                    );
                }

                let uv_fit = cli_soft_profile_fit(&mut child_uv, PROCESS, "child_uv_only")?;
                let h_fit = cli_soft_profile_fit(&mut child_h, PROCESS, "child_soft_ir")?;
                assert_eq!(
                    uv_fit.ray_fingerprint, h_fit.ray_fingerprint,
                    "child-only U and H controls must evaluate the identical routed soft ray",
                );
                assert!(
                    uv_fit.r_squared >= 0.98 && h_fit.r_squared >= 0.98,
                    "child-only q_g fits must be resolved: U={uv_fit:?}, H={h_fit:?}",
                );
                assert!(
                    uv_fit.scaling < -1.0,
                    "the ordinary child counterterm must retain the spurious q_g enhancement: U={}",
                    uv_fit.scaling,
                );
                assert!(
                    h_fit.scaling > -0.5,
                    "the completed child H counterterm must restore marginal soft scaling: H={}",
                    h_fit.scaling,
                );
                assert!(
                    h_fit.scaling - uv_fit.scaling >= 1.5,
                    "the completed child H counterterm must improve the matched q_g ray by two powers: U={}, H={}",
                    uv_fit.scaling,
                    h_fit.scaling,
                );

                if !explicit_orientation_sum_only {
                let uv_orientations = cli_soft_profile_orientation_fits(
                    &mut child_uv,
                    PROCESS,
                    "child_uv_only",
                )?;
                let h_orientations = cli_soft_profile_orientation_fits(
                    &mut child_h,
                    PROCESS,
                    "child_soft_ir",
                )?;
                assert_eq!(
                    uv_orientations.keys().collect::<Vec<_>>(),
                    h_orientations.keys().collect::<Vec<_>>(),
                    "U and H must profile the same orientation labels",
                );
                assert_eq!(
                    h_orientations.len(),
                    30,
                    "the top-bubble fixture must retain all 30 acyclic orientations",
                );
                let mut divergent_u_orientations = 0;
                for (orientation, h_orientation) in &h_orientations {
                    let u_orientation = &uv_orientations[orientation];
                    assert_eq!(
                        h_orientation.ray_fingerprint, u_orientation.ray_fingerprint,
                        "orientation {orientation} used different U/H soft rays",
                    );
                    assert!(
                        h_orientation.scaling > -0.5,
                        "orientation {orientation} retains a soft enhancement after H: {h_orientation:?}",
                    );
                    if u_orientation.scaling < -1.0 {
                        divergent_u_orientations += 1;
                        assert!(
                            u_orientation.r_squared >= 0.98
                                && h_orientation.r_squared >= 0.98
                                && h_orientation.scaling - u_orientation.scaling >= 1.5,
                            "orientation {orientation} does not exhibit the resolved two-power U-to-H improvement: U={u_orientation:?}, H={h_orientation:?}",
                        );
                    }
                }
                assert_eq!(
                    divergent_u_orientations, 8,
                    "the matched U control must expose the eight orientation-local soft failures",
                );

                }
                clean_test(&child_uv.cli_settings.state.folder);
                clean_test(&child_h.cli_settings.state.folder);
                Ok(())
            })?;
            test.join()
                .map_err(|_| eyre::eyre!("top-bubble child-only profile thread panicked"))??;
        }
        Ok(())
    }

    #[test]
    #[serial_test::serial]
    fn gamma_star_ddbar_top_bubble_has_two_power_soft_improvement() -> Result<()> {
        for (explicit_orientation_sum_only, project_local_4d) in
            [(false, false), (true, false), (true, true)]
        {
            println!(
                "explicit_orientation_sum_only={explicit_orientation_sum_only}, project_local_4d={project_local_4d}"
            );
            let test = std::thread::Builder::new()
            .name("gamma-star-ddbar-top-bubble".to_string())
            .stack_size(128 * 1024 * 1024)
            .spawn(move || -> Result<()> {
                const PROCESS: &str = "gamma_star_ddbar_top_bubble";

                let mut bare = local_ct_cli(
                    "gamma_star_ddbar_top_bubble_bare.toml",
                    "gamma_star_ddbar_top_bubble_bare",
                    false, project_local_4d, explicit_orientation_sum_only)?;
                let mut uv_only = local_ct_cli(
                    "gamma_star_ddbar_top_bubble_uv_only.toml",
                    "gamma_star_ddbar_top_bubble_uv_only",
                    false, project_local_4d, explicit_orientation_sum_only)?;
                let mut soft_ir = local_ct_cli(
                    "gamma_star_ddbar_top_bubble_soft_ir.toml",
                    "gamma_star_ddbar_top_bubble_soft_ir",
                    false, project_local_4d, explicit_orientation_sum_only)?;

                for (cli, integrand) in [
                    (&bare, "bare"),
                    (&uv_only, "uv_only"),
                    (&soft_ir, "soft_ir"),
                ] {
                    assert_generated_local_integrand(cli, PROCESS, integrand, 6)?;
                    assert_eq!(
                        cli.default_runtime_settings.kinematics.e_cm, 300.0,
                        "{integrand} must use Q=300 GeV",
                    );
                    for (parameter, expected_value) in [("MT", 173.0), ("WT", 0.0)] {
                        let value = cli
                            .state
                            .model
                            .get_parameter(parameter)
                            .value
                            .as_ref()
                            .ok_or_else(|| {
                                eyre::eyre!(
                                    "{integrand} has no resolved model parameter {parameter}"
                                )
                            })?;
                        assert_eq!(
                            (value.re.0, value.im.0),
                            (expected_value, 0.0),
                            "{integrand} must use {parameter}={expected_value}",
                        );
                    }
                    let process_id = cli.state.resolve_process_ref(Some(
                        &ProcessRef::Unqualified(PROCESS.to_string()),
                    ))?;
                    let ProcessCollection::Amplitudes(amplitudes) =
                        &cli.state.process_list.processes[process_id].collection
                    else {
                        return Err(eyre::eyre!("top-bubble fixture must be an amplitude"));
                    };
                    let color_closure = amplitudes
                        .get(integrand)
                        .ok_or_else(|| eyre::eyre!("missing integrand {integrand}"))?
                        .graphs[0]
                        .graph
                        .global_prefactor
                        .num
                        .to_canonical_string();
                    assert!(
                        color_closure.contains("1/3")
                            && color_closure.contains("cof(3")
                            && color_closure.contains("dind"),
                        "{integrand} did not import the explicit delta_i^j/N_c color closure: {color_closure}",
                    );
                }
                let fixture_source = include_str!(
                    "../resources/graphs/gamma_star_ddbar_top_bubble.dot"
                );
                assert!(
                    fixture_source.contains(
                        "num = \"1/3*spenso::g(spenso::dind(spenso::cof(3,gammalooprs::hedge(0))),spenso::cof(3,gammalooprs::hedge(1)))\""
                    ),
                    "the top-bubble fixture must close the external fundamental string with delta_i^j/N_c",
                );
                let canonical_route = |cli: &CLIState,
                                       integrand_name: &str|
                 -> Result<gammalooprs::graph::LoopMomentumBasis> {
                    let process_id = cli.state.resolve_process_ref(Some(
                        &ProcessRef::Unqualified(PROCESS.to_string()),
                    ))?;
                    let ProcessCollection::Amplitudes(amplitudes) =
                        &cli.state.process_list.processes[process_id].collection
                    else {
                        return Err(eyre::eyre!("top-bubble fixture must be an amplitude"));
                    };
                    Ok(amplitudes
                        .get(integrand_name)
                        .ok_or_else(|| eyre::eyre!("missing integrand {integrand_name}"))?
                        .graphs[0]
                        .graph
                        .loop_momentum_basis
                        .clone())
                };
                let bare_route = canonical_route(&bare, "bare")?;
                assert_eq!(
                    bare_route
                        .loop_edges
                        .iter()
                        .copied()
                        .map(usize::from)
                        .collect::<Vec<_>>(),
                    vec![6, 8],
                    "the canonical route must use q_g=e6 and the top-loop momentum e8"
                );
                assert_eq!(
                    bare_route,
                    canonical_route(&uv_only, "uv_only")?,
                    "bare and ordinary-U variants must use the same canonical route"
                );
                assert_eq!(
                    bare_route,
                    canonical_route(&soft_ir, "soft_ir")?,
                    "bare and soft-refined variants must use the same canonical route"
                );
                let expected = CTIdentifier::new(
                    [-21, 21].into_iter().collect(),
                    Some([6].into_iter().collect()),
                );
                let bare_classified = classified_local_spinneys(&bare, PROCESS, "bare")?;
                assert!(
                    bare_classified.iter().all(|component| {
                        component.n_components == 0
                            && component.edge_ids.is_empty()
                            && component.scheme == ApproximationType::MUV
                            && component.dod == 0
                    }),
                    "the top-bubble bare control must contain no local counterterm components: {bare_classified:?}",
                );
                let uv_only_classified =
                    classified_local_spinneys(&uv_only, PROCESS, "uv_only")?;
                let uv_only_components = uv_only_classified
                    .iter()
                    .filter(|spinney| spinney.n_components == 1)
                    .collect::<Vec<_>>();
                let uv_only_inert = uv_only_classified
                    .iter()
                    .filter(|spinney| spinney.n_components == 0)
                    .collect::<Vec<_>>();
                assert!(
                    uv_only_components.len() == 2
                        && uv_only_classified.len()
                            == uv_only_components.len() + uv_only_inert.len()
                        && uv_only_inert.iter().all(|component| {
                            component.edge_ids.is_empty()
                                && component.scheme == ApproximationType::MUV
                                && component.dod == 0
                        })
                        && uv_only_components
                            .iter()
                            .all(|component| component.scheme == ApproximationType::MUV)
                        && uv_only_components.iter().any(|component| {
                            component.edge_ids == [7, 8]
                                && component.dod == 2
                                && component.identifier == expected
                        })
                        && uv_only_components.iter().any(|component| {
                            component.edge_ids == [3, 4, 5, 6, 7, 8]
                                && component.dod == 0
                        }),
                    "uv_only must contain exactly the ordinary-U top bubble and logarithmic overall vertex: {uv_only_components:?}",
                );
                let classified = classified_local_spinneys(&soft_ir, PROCESS, "soft_ir")?;
                let soft_components = classified
                    .iter()
                    .filter(|spinney| {
                        spinney.n_components == 1 && spinney.scheme == ApproximationType::IR
                    })
                    .collect::<Vec<_>>();
                let soft_inert = classified
                    .iter()
                    .filter(|spinney| spinney.n_components == 0)
                    .collect::<Vec<_>>();
                assert!(
                    soft_components.len() == 1
                        && soft_components[0].identifier == expected
                        && soft_components[0].edge_ids == [7, 8]
                        && soft_components[0].dod == 2,
                    "soft_ir must select exactly the degree-2 massive-top gluon self-energy: {soft_components:?}"
                );
                let child_edges = &soft_components[0].edge_ids;
                let containing_muv = classified
                    .iter()
                    .filter(|spinney| {
                        spinney.n_components == 1
                            && spinney.scheme == ApproximationType::MUV
                            && spinney.edge_ids.len() > child_edges.len()
                            && child_edges
                                .iter()
                                .all(|edge| spinney.edge_ids.contains(edge))
                    })
                    .collect::<Vec<_>>();
                assert!(
                    containing_muv.len() == 1
                        && containing_muv[0].dod == 0
                        && containing_muv[0].edge_ids == [3, 4, 5, 6, 7, 8]
                        && classified.len()
                            == soft_components.len() + containing_muv.len() + soft_inert.len()
                        && soft_inert.iter().all(|component| {
                            component.edge_ids.is_empty()
                                && component.scheme == ApproximationType::MUV
                                && component.dod == 0
                        }),
                    "the top-bubble IR child must be strictly nested in a logarithmic ordinary-U vertex component: child={:?}, containing={containing_muv:?}",
                    soft_components[0]
                );

                let process_id = soft_ir
                    .state
                    .resolve_process_ref(Some(&ProcessRef::Unqualified(PROCESS.to_string())))?;
                let structural_forest = soft_ir.state.process_list.processes[process_id].export_uv_forest_graph(
"soft_ir",
0,
                    &gammalooprs::uv::export::UVForestExportSettings { computed: false },
                )?;
                assert!(
                    structural_forest.node_terms.is_empty(),
                    "a structure-only top-bubble export must not construct local atoms",
                );
                assert_eq!(
                    structural_forest.forest_dot.matches("foata=").count(),
                    4,
                    "the top-bubble forest must contain exactly the root, H_2 child, U_0 parent, and child-to-parent chain",
                );

                let bare_uv = cli_uv_profile_pass_fail(&mut bare, PROCESS, "bare")?;
                assert!(
bare_uv.total > 0 && bare_uv.resolved > 0 && bare_uv.dod_failures > 0
 && (explicit_orientation_sum_only || (bare_uv.orientation_total > 0 && bare_uv.orientation_resolved > 0 && bare_uv.orientation_dod_failures > 0)),
                    "bare top-bubble vertex must have a resolved per-orientation fit above the UV bound: {bare_uv:?}"
                );
                let uv_only_uv = cli_uv_profile_pass_fail(&mut uv_only, PROCESS, "uv_only")?;
                assert!(
uv_only_uv.total > 0 && uv_only_uv.resolved > 0 && uv_only_uv.failed == 0
 && (explicit_orientation_sum_only || (uv_only_uv.orientation_total > 0 && uv_only_uv.orientation_resolved > 0 && uv_only_uv.orientation_failed == 0)),
                    "uv_only must pass a non-vacuous CLI UV profile: {uv_only_uv:?}"
                );
                let soft_ir_uv = cli_uv_profile_pass_fail(&mut soft_ir, PROCESS, "soft_ir")?;
                assert!(
soft_ir_uv.total > 0 && soft_ir_uv.resolved > 0 && soft_ir_uv.failed == 0
 && (explicit_orientation_sum_only || (soft_ir_uv.orientation_total > 0 && soft_ir_uv.orientation_resolved > 0 && soft_ir_uv.orientation_failed == 0)),
                    "soft_ir must pass a non-vacuous CLI UV profile: {soft_ir_uv:?}"
                );
                assert_eq!(
                    bare_uv.profile_identity, uv_only_uv.profile_identity,
                    "bare and ordinary-U variants must profile the same UV LMBs, subsets, and orientations",
                );
                assert_eq!(
                    bare_uv.profile_identity, soft_ir_uv.profile_identity,
                    "bare and soft-refined variants must profile the same UV LMBs, subsets, and orientations",
                );

                let bare_fit = cli_soft_profile_fit(&mut bare, PROCESS, "bare")?;
                let uv_fit = cli_soft_profile_fit(&mut uv_only, PROCESS, "uv_only")?;
                let soft_fit = cli_soft_profile_fit(&mut soft_ir, PROCESS, "soft_ir")?;
                assert_eq!(
                    bare_fit.ray_fingerprint, uv_fit.ray_fingerprint,
                    "bare and ordinary-U variants must evaluate the identical routed soft ray"
                );
                assert_eq!(
                    bare_fit.ray_fingerprint, soft_fit.ray_fingerprint,
                    "bare and soft-refined variants must evaluate the identical routed soft ray"
                );
                assert!(
                    [
                        bare_fit.r_squared,
                        uv_fit.r_squared,
                        soft_fit.r_squared,
                    ]
                        .into_iter()
                        .all(|r_squared| r_squared >= 0.98),
                    "top-bubble q_g fits must have R-squared >= 0.98: bare={}, U={}, H={}",
                    bare_fit.r_squared,
                    uv_fit.r_squared,
                    soft_fit.r_squared,
                );
                assert!(
                    uv_fit.scaling < -1.0,
                    "the ordinary-U vertex must retain the spurious q_g enhancement: U={}",
                    uv_fit.scaling,
                );
                assert!(
                    soft_fit.scaling > -0.5,
                    "the completed H counterterm must restore the physical marginal soft scaling; got {}",
                    soft_fit.scaling,
                );
                assert!(
                    soft_fit.scaling - uv_fit.scaling >= 1.5,
                    "H must improve the q_g scaling by approximately two powers: U={}, H={}",
                    uv_fit.scaling,
                    soft_fit.scaling,
                );

                if !explicit_orientation_sum_only {
                let uv_orientations =
                    cli_soft_profile_orientation_fits(&mut uv_only, PROCESS, "uv_only")?;
                let h_orientations =
                    cli_soft_profile_orientation_fits(&mut soft_ir, PROCESS, "soft_ir")?;
                assert_eq!(
                    uv_orientations.keys().collect::<Vec<_>>(),
                    h_orientations.keys().collect::<Vec<_>>(),
                    "the nested U and H woods must profile the same orientation labels",
                );
                assert_eq!(
                    h_orientations.len(),
                    30,
                    "the nested top-bubble fixture must retain all 30 acyclic orientations",
                );
                let mut divergent_u_orientations = 0;
                for (orientation, h_orientation) in &h_orientations {
                    let u_orientation = &uv_orientations[orientation];
                    assert_eq!(
                        h_orientation.ray_fingerprint, u_orientation.ray_fingerprint,
                        "nested orientation {orientation} used different U/H soft rays",
                    );
                    assert!(
                        h_orientation.scaling > -0.5,
                        "the outer U_0 destroyed the H child's soft improvement in orientation {orientation}: {h_orientation:?}",
                    );
                    if u_orientation.scaling < -1.0 {
                        divergent_u_orientations += 1;
                        assert!(
                            u_orientation.r_squared >= 0.98
                                && h_orientation.r_squared >= 0.98
                                && h_orientation.scaling - u_orientation.scaling >= 1.5,
                            "nested orientation {orientation} does not exhibit the resolved two-power U-to-H improvement: U={u_orientation:?}, H={h_orientation:?}",
                        );
                    }
                }
                assert_eq!(
                    divergent_u_orientations, 8,
                    "the nested U control must expose the same eight orientation-local soft failures as the isolated child",
                );

                }
                for cli in [&bare, &uv_only, &soft_ir] {
                    clean_test(&cli.cli_settings.state.folder);
                }
                Ok(())
            })?;
            test.join()
                .map_err(|_| eyre::eyre!("top-bubble profile test thread panicked"))??;
        }
        Ok(())
    }

    #[test]
    #[serial_test::serial]
    fn aa_aa_2l_gl00_matches_across_local_uv_routes() -> Result<()> {
        AA_AA_2L_UV_RICH_INSPECT.run_graph_local_uv_route_comparison(
            "GL00",
            "examples/cli/aa_aa/2L/graphs/GL00.dot",
            "2L",
            &[0.11, -0.07, 0.19, -0.13, 0.05, 0.29],
        )
    }

    #[test]
    #[serial_test::serial]
    fn aa_aa_2l_gl01_matches_across_local_uv_routes() -> Result<()> {
        AA_AA_2L_UV_RICH_INSPECT.run_graph_local_uv_route_comparison(
            "GL01",
            "examples/cli/aa_aa/2L/graphs/GL01.dot",
            "2L",
            &[0.11, -0.07, 0.19, -0.13, 0.05, 0.29],
        )
    }

    #[test]
    #[serial_test::serial]
    fn aa_aa_3l_gl262_matches_across_local_uv_routes() -> Result<()> {
        AA_AA_2L_UV_RICH_INSPECT.run_graph_local_uv_route_comparison(
            "GL262",
            "examples/cli/aa_aa/3L/graphs/processes/amplitudes/aa_aa/3L/GL262.dot",
            "3L",
            &[0.11, -0.07, 0.19, -0.13, 0.05, 0.29, 0.17, -0.23, 0.31],
        )
    }

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
}
