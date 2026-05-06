use std::collections::BTreeMap;

use color_eyre::Result;
use gammaloop_api::commands::{
    Profile,
    integrate::Integrate,
    profile::{InfraRedProfile, UltraVioletProfile},
};
use gammaloop_api::state::ProcessRef;
use gammaloop_integration_tests::{
    CLIState, clean_test, get_example_cli, get_test_cli, get_tests_workspace_path,
};
use gammalooprs::settings::runtime::{IntegralEstimate, SlotIntegrationResult};
use gammalooprs::utils::F;
use spenso::algebra::complex::Complex;
use tabled::{Table, Tabled};

#[derive(Default)]
struct IntegratedUvTargets {
    no_integrated: Option<Complex<F<f64>>>,
    integrated: Option<Complex<F<f64>>>,
}

impl IntegratedUvTargets {
    fn any(&self) -> bool {
        self.no_integrated.is_some() || self.integrated.is_some()
    }
}

struct IntegratedUvCase<'a> {
    run_card: &'a str,
    test_name: &'a str,
    process: &'a str,
    integrand_name: &'a str,
    original_m_uv: f64,
    shifted_m_uv: f64,
    original_mu_r: f64,
    shifted_mu_r: f64,
    skip_uv_profile: bool,
    targets: IntegratedUvTargets,
    min_change_sigma: Option<f64>,
    min_mu_r_change_sigma: Option<f64>,
}

struct IntegratedUvResults {
    no_integrated: Option<IntegralEstimate>,
    integrated: IntegralEstimate,
    integrated_per_graph: Option<BTreeMap<String, IntegralEstimate>>,
    integrated_muv_shifted_per_graph: Option<BTreeMap<String, IntegralEstimate>>,
    integrated_mur_shifted_per_graph: Option<BTreeMap<String, IntegralEstimate>>,
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
    mur_dependence_passed: Option<bool>,
    mur_rows: Vec<CheckRow>,
    target_passed: Option<bool>,
    target_threshold: Option<&'static str>,
    ct_change_passed: Option<bool>,
    ct_change_sigma: Option<f64>,
    ct_change_threshold: Option<f64>,
    mur_change_sigma: Option<f64>,
    mur_change_threshold: Option<f64>,
    error: Option<String>,
}

impl IntegratedUvCaseResult {
    fn all_passed(&self) -> bool {
        self.error.is_none()
            && self.uv_profile_passed.unwrap_or(true)
            && self.muv_invariance_passed
            && self.mur_dependence_passed.unwrap_or(true)
            && self.target_passed.unwrap_or(true)
            && self.ct_change_passed.unwrap_or(true)
    }
}

fn set_fast_deterministic_integrator(cli: &mut CLIState) -> Result<()> {
    cli.run_command(
        "set process kv integrator.target_relative_accuracy=0.001 \
         integrator.n_increase=0 integrator.n_start=50000 \
         integrator.n_max=200000 integrator.seed=1337",
    )
}

fn no_integrated_process_name(process: &str) -> String {
    format!("{process}_no_integrated_UV")
}

fn shifted_muv_integrand_name(integrand_name: &str) -> String {
    format!("{integrand_name}_muv_shifted")
}

fn shifted_mur_integrand_name(integrand_name: &str) -> String {
    format!("{integrand_name}_mur_shifted")
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

    if case.min_mu_r_change_sigma.is_some()
        && (case.original_mu_r - case.shifted_mu_r).abs() > f64::EPSILON
    {
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

fn process_integrand_exists(cli: &mut CLIState, process: &str, integrand_name: &str) -> bool {
    cli.state
        .find_integrand_ref(
            Some(&ProcessRef::Unqualified(process.to_string())),
            Some(&integrand_name.to_string()),
        )
        .is_ok()
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

fn run_integrated_uv_integration(
    cli: &mut CLIState,
    case: &IntegratedUvCase<'_>,
) -> Result<IntegratedUvResults> {
    set_fast_deterministic_integrator(cli)?;
    let no_integrated_process = no_integrated_process_name(case.process);
    let has_no_integrated =
        process_integrand_exists(cli, &no_integrated_process, case.integrand_name);
    let mut processes = vec![ProcessRef::Unqualified(case.process.to_string())];
    let mut integrand_names = vec![case.integrand_name.to_string()];

    if has_no_integrated {
        processes.push(ProcessRef::Unqualified(no_integrated_process.clone()));
        integrand_names.push(case.integrand_name.to_string());
    }

    let muv_shifted_name = ((case.original_m_uv - case.shifted_m_uv).abs() > f64::EPSILON)
        .then(|| shifted_muv_integrand_name(case.integrand_name));
    if let Some(name) = &muv_shifted_name {
        processes.push(ProcessRef::Unqualified(case.process.to_string()));
        integrand_names.push(name.clone());
    }

    let mur_shifted_name = (case.min_mu_r_change_sigma.is_some()
        && (case.original_mu_r - case.shifted_mu_r).abs() > f64::EPSILON)
        .then(|| shifted_mur_integrand_name(case.integrand_name));
    if let Some(name) = &mur_shifted_name {
        processes.push(ProcessRef::Unqualified(case.process.to_string()));
        integrand_names.push(name.clone());
    }

    let integration_result = Integrate {
        process: processes,
        integrand_name: integrand_names,
        workspace_path: Some(get_tests_workspace_path().join(format!(
            "{}/integration_workspace_{}",
            case.test_name, case.integrand_name
        ))),
        n_cores: Some(1),
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
    let integrated_mur_shifted_per_graph = mur_shifted_name
        .as_ref()
        .map(|name| {
            graph_estimates_for_slot(
                cli,
                integration_result
                    .slot(&format!("{}@{}", case.process, name))
                    .expect("integrated shifted-mu_r slot should exist"),
            )
        })
        .transpose()?
        .flatten();

    Ok(IntegratedUvResults {
        no_integrated: has_no_integrated.then(|| {
            integration_result
                .slot(&format!(
                    "{}@{}",
                    no_integrated_process, case.integrand_name
                ))
                .expect("no-integrated slot should exist")
                .integral
                .clone()
        }),
        integrated: integrated_slot.integral.clone(),
        integrated_per_graph,
        integrated_muv_shifted_per_graph,
        integrated_mur_shifted_per_graph,
    })
}

fn integrated_uv_targets_pass(
    results: &IntegratedUvResults,
    targets: &IntegratedUvTargets,
) -> Result<bool> {
    if let Some(no_integrated_target) = &targets.no_integrated {
        let Some(no_integrated_result) = results.no_integrated.as_ref() else {
            return Err(eyre::eyre!(
                "no-integrated result should exist for no-integrated target checks"
            ));
        };
        if !no_integrated_result.is_compatible_with_target(*no_integrated_target, 2) {
            return Ok(false);
        }
    }

    if let Some(integrated_target) = &targets.integrated
        && !results
            .integrated
            .is_compatible_with_target(*integrated_target, 2)
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

fn integrated_uv_result_change_sigma(results: &IntegratedUvResults) -> f64 {
    integral_estimate_change_sigma(
        &results.integrated,
        results
            .no_integrated
            .as_ref()
            .expect("no-integrated result should exist for ct-diff checks"),
    )
}

fn integrated_uv_meaningfully_changes_result(
    results: &IntegratedUvResults,
    min_change_sigma: f64,
) -> bool {
    integrated_uv_result_change_sigma(results) >= min_change_sigma
}

fn integrated_uv_profile_passes(
    cli: &mut CLIState,
    process: &str,
    integrand_name: &str,
) -> Result<bool> {
    let res = Profile::UltraViolet(UltraVioletProfile {
        process: Some(ProcessRef::Unqualified(process.to_string())),
        integrand_name: Some(integrand_name.to_string()),
        min_scale_exponent: 4.0,
        max_scale_exponent: 5.0,
        n_points: 100,
        per_orientation: true,
        ..Default::default()
    })
    .run(&mut cli.state, &cli.cli_settings)?;
    let uv = res.unwrap_uv();
    Ok(uv.pass_fail(-0.9).failed == 0)
}

fn integrated_result_muv_invariance_rows(results: &IntegratedUvResults) -> Result<Vec<CheckRow>> {
    match (
        &results.integrated_per_graph,
        &results.integrated_muv_shifted_per_graph,
    ) {
        (Some(baseline_graphs), Some(shifted_graphs)) if !baseline_graphs.is_empty() => {
            let mut rows = Vec::new();
            for (graph, baseline_graph) in baseline_graphs {
                let shifted_graph = shifted_graphs.get(graph).ok_or_else(|| {
                    eyre::eyre!("shifted m_uv graph breakdown is missing graph '{}'", graph)
                })?;
                rows.push(CheckRow {
                    check: "m_uv",
                    graph: graph.clone(),
                    passed: shifted_graph.is_compatible_with_result(baseline_graph, 2),
                    sigma: Some(integral_estimate_change_sigma(
                        shifted_graph,
                        baseline_graph,
                    )),
                    threshold: Some("compatible within 2sigma".to_string()),
                });
            }
            for graph in shifted_graphs.keys() {
                if !baseline_graphs.contains_key(graph) {
                    return Err(eyre::eyre!(
                        "baseline graph breakdown is missing shifted m_uv graph '{}'",
                        graph
                    ));
                }
            }
            Ok(rows)
        }
        (None, None) => Err(eyre::eyre!(
            "m_uv comparison requires graph breakdowns for baseline and shifted results"
        )),
        (Some(_), None) => Err(eyre::eyre!(
            "shifted m_uv result has no graph breakdown to compare against baseline"
        )),
        (None, Some(_)) => Err(eyre::eyre!(
            "baseline result has no graph breakdown to compare against shifted m_uv"
        )),
        (Some(_), Some(_)) => Err(eyre::eyre!(
            "m_uv comparison requires a non-empty baseline graph breakdown"
        )),
    }
}

fn integrated_result_mu_r_change_rows(
    results: &IntegratedUvResults,
    min_change_sigma: f64,
) -> Result<Vec<CheckRow>> {
    match (
        &results.integrated_per_graph,
        &results.integrated_mur_shifted_per_graph,
    ) {
        (Some(baseline_graphs), Some(shifted_graphs)) if !baseline_graphs.is_empty() => {
            let mut rows = Vec::new();
            for (graph, baseline_graph) in baseline_graphs {
                let shifted_graph = shifted_graphs.get(graph).ok_or_else(|| {
                    eyre::eyre!("shifted mu_r graph breakdown is missing graph '{}'", graph)
                })?;
                let sigma = integral_estimate_change_sigma(shifted_graph, baseline_graph);
                rows.push(CheckRow {
                    check: "mu_r",
                    graph: graph.clone(),
                    passed: sigma >= min_change_sigma,
                    sigma: Some(sigma),
                    threshold: Some(format!(">={min_change_sigma:.1}sigma")),
                });
            }
            for graph in shifted_graphs.keys() {
                if !baseline_graphs.contains_key(graph) {
                    return Err(eyre::eyre!(
                        "baseline graph breakdown is missing shifted mu_r graph '{}'",
                        graph
                    ));
                }
            }
            Ok(rows)
        }
        (None, None) => Err(eyre::eyre!(
            "mu_r comparison requires graph breakdowns for baseline and shifted results"
        )),
        (Some(_), None) => Err(eyre::eyre!(
            "shifted mu_r result has no graph breakdown to compare against baseline"
        )),
        (None, Some(_)) => Err(eyre::eyre!(
            "baseline result has no graph breakdown to compare against shifted mu_r"
        )),
        (Some(_), Some(_)) => Err(eyre::eyre!(
            "mu_r comparison requires a non-empty baseline graph breakdown"
        )),
    }
}

fn run_integrated_uv_case(case: &IntegratedUvCase<'_>) -> IntegratedUvCaseResult {
    let mut outcome = IntegratedUvCaseResult {
        graph: case.test_name.to_string(),
        uv_profile_passed: (!case.skip_uv_profile).then_some(false),
        muv_invariance_passed: false,
        muv_rows: Vec::new(),
        mur_dependence_passed: case.min_mu_r_change_sigma.map(|_| false),
        mur_rows: Vec::new(),
        target_passed: case.targets.any().then_some(false),
        target_threshold: case.targets.any().then_some("2sigma"),
        ct_change_passed: case.min_change_sigma.map(|_| false),
        ct_change_sigma: None,
        ct_change_threshold: case.min_change_sigma,
        mur_change_sigma: None,
        mur_change_threshold: case.min_mu_r_change_sigma,
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
        add_integrated_uv_scale_variants(&mut cli, case)?;

        let no_integrated_process = no_integrated_process_name(case.process);
        let has_no_integrated =
            process_integrand_exists(&mut cli, &no_integrated_process, case.integrand_name);
        if !case.skip_uv_profile {
            outcome.uv_profile_passed = Some(
                integrated_uv_profile_passes(&mut cli, case.process, case.integrand_name)?
                    && (!has_no_integrated
                        || integrated_uv_profile_passes(
                            &mut cli,
                            &no_integrated_process,
                            case.integrand_name,
                        )?),
            );
        }

        let baseline = run_integrated_uv_integration(&mut cli, case)?;

        if let Some(min_change_sigma) = case.min_change_sigma {
            if baseline.no_integrated.is_none() {
                return Err(eyre::eyre!(
                    "ct-diff check requested for {} but no {} process exists",
                    case.process,
                    no_integrated_process
                ));
            }
            let change_sigma = integrated_uv_result_change_sigma(&baseline);
            outcome.ct_change_sigma = Some(change_sigma);
            outcome.ct_change_passed = Some(integrated_uv_meaningfully_changes_result(
                &baseline,
                min_change_sigma,
            ));
        }

        outcome.muv_rows = integrated_result_muv_invariance_rows(&baseline)?;
        outcome.muv_invariance_passed = outcome.muv_rows.iter().all(|row| row.passed);

        if let Some(min_mu_r_change_sigma) = case.min_mu_r_change_sigma {
            outcome.mur_rows =
                integrated_result_mu_r_change_rows(&baseline, min_mu_r_change_sigma)?;
            let change_sigma = outcome
                .mur_rows
                .iter()
                .filter_map(|row| row.sigma)
                .fold(f64::INFINITY, f64::min);
            outcome.mur_change_sigma = Some(change_sigma);
            outcome.mur_dependence_passed = Some(outcome.mur_rows.iter().all(|row| row.passed));
        }

        if case.targets.any() {
            if case.targets.no_integrated.is_some() && baseline.no_integrated.is_none() {
                return Err(eyre::eyre!(
                    "no-integrated target check requested for {} but no {} process exists",
                    case.process,
                    no_integrated_process
                ));
            }
            outcome.target_passed = Some(integrated_uv_targets_pass(&baseline, &case.targets)?);
        }

        Ok(())
    })();

    if let Err(err) = case_result {
        outcome.error = Some(format!("{err:#}"));
    }

    clean_test(&cli.cli_settings.state.folder);
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
                sigma: sigma(result.mur_change_sigma),
                threshold: result
                    .mur_change_threshold
                    .map_or_else(|| "-".to_string(), |threshold| format!(">={threshold:.1}σ")),
                error: error.clone(),
            });
        }
        rows.push(SummaryRow {
            case: result.graph.clone(),
            check: "target".to_string(),
            graph: "-".to_string(),
            status: status(result.target_passed),
            sigma: "-".to_string(),
            threshold: result.target_threshold.unwrap_or("-").to_string(),
            error: error.clone(),
        });
        rows.push(SummaryRow {
            case: result.graph.clone(),
            check: "ct_diff".to_string(),
            graph: "-".to_string(),
            status: status(result.ct_change_passed),
            sigma: sigma(result.ct_change_sigma),
            threshold: result
                .ct_change_threshold
                .map_or_else(|| "-".to_string(), |threshold| format!(">={threshold:.1}σ")),
            error,
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
        original_m_uv: 20.0,
        shifted_m_uv: 7.0,
        original_mu_r: 3.0,
        shifted_mu_r: 9.0,
        skip_uv_profile: false,
        targets: IntegratedUvTargets {
            no_integrated: Some(Complex::new(F(1.4693e-2), F(0.0))),
            integrated: Some(Complex::new(F(2.684e-3), F(0.0))),
        },
        min_change_sigma: Some(5.0),
        min_mu_r_change_sigma: Some(5.0),
    });
}

#[test]
fn dod1_bubble_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/dod1_bubble",
        test_name: "dod1_bubble",
        process: "bubble_dod1",
        integrand_name: "scalar_bubble_below_thres",
        original_m_uv: 20.0,
        shifted_m_uv: 7.0,
        original_mu_r: 3.0,
        shifted_mu_r: 19.0,
        skip_uv_profile: false,
        targets: IntegratedUvTargets {
            no_integrated: Some(Complex::new(F(7.358320108607984e-3), F(0.0))),
            integrated: Some(Complex::new(F(1.351_493_784_227_627e-3), F(0.0))),
        },
        min_change_sigma: Some(5.0),
        min_mu_r_change_sigma: Some(5.0),
    });
}

#[test]
fn dod2_bubble_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/dod2_bubble",
        test_name: "dod2_bubble",
        process: "bubble_dod2",
        integrand_name: "scalar_bubble_below_thres",
        original_m_uv: 20.0,
        shifted_m_uv: 7.0,
        original_mu_r: 3.0,
        shifted_mu_r: 19.0,
        skip_uv_profile: false,
        targets: IntegratedUvTargets {
            no_integrated: Some(Complex::new(F(-0.5940830828502411), F(0.0))),
            integrated: Some(Complex::new(F(1.2143596454658382e-2), F(0.0))),
        },
        min_change_sigma: Some(5.0),
        min_mu_r_change_sigma: Some(5.0),
    });
}

#[test]
fn epem_a_bbx_amp_uv() {
    run_single_integrated_uv_case(&IntegratedUvCase {
        run_card: "uv/epem_a_bbx_amp",
        test_name: "epem_a_bbx_amp",
        process: "epem_a_bbx",
        integrand_name: "epem_a_bbx",
        original_m_uv: 20.0,
        shifted_m_uv: 7.0,
        original_mu_r: 20.0,
        shifted_mu_r: 9.0,
        skip_uv_profile: false,
        targets: IntegratedUvTargets::default(),
        min_change_sigma: Some(5.0),
        min_mu_r_change_sigma: Some(5.0),
    });
}

mod slow {
    use super::*;

    #[test]
    fn aa_aa_gl00_uv() {
        run_single_integrated_uv_case(&IntegratedUvCase {
            run_card: "uv/aa_aa_GL00",
            test_name: "aa_aa_GL00",
            process: "aa_aa",
            integrand_name: "2L",
            original_m_uv: 91.188,
            shifted_m_uv: 364.752,
            original_mu_r: 91.188,
            shifted_mu_r: 364.752,
            skip_uv_profile: false,
            targets: IntegratedUvTargets::default(),
            min_change_sigma: Some(5.0),
            min_mu_r_change_sigma: Some(5.0),
        });
    }

    #[test]
    fn epem_ttxh_gl00_uv() {
        run_single_integrated_uv_case(&IntegratedUvCase {
            run_card: "uv/epem_ttxh_GL00",
            test_name: "epem_ttxh_gl00",
            process: "epem_a_tth",
            integrand_name: "NLO",
            original_m_uv: 20.0,
            shifted_m_uv: 7.0,
            original_mu_r: 3.0,
            shifted_mu_r: 9.0,
            skip_uv_profile: false,
            targets: IntegratedUvTargets::default(),
            min_change_sigma: Some(5.0),
            min_mu_r_change_sigma: Some(5.0),
        });
    }

    #[test]
    fn ad_ad_with_gluon_correction_uv() {
        run_single_integrated_uv_case(&IntegratedUvCase {
            run_card: "uv/ad_ad_with_gluon_correction",
            test_name: "ad_ad_with_gluon_correction",
            process: "adad",
            integrand_name: "adad_gluon",
            original_m_uv: 20.0,
            shifted_m_uv: 7.0,
            original_mu_r: 20.0,
            shifted_mu_r: 9.0,
            skip_uv_profile: false,
            targets: IntegratedUvTargets::default(),
            min_change_sigma: Some(5.0),
            min_mu_r_change_sigma: Some(5.0),
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
mod failing {
    use super::*;

    #[test]
    fn epem_a_ddx_xs_nlo_uv() {
        run_single_integrated_uv_case(&IntegratedUvCase {
            run_card: "uv/epem_a_ddx_xs_nlo",
            test_name: "epem_a_ddx_xs_nlo",
            process: "epem_a_ddx",
            integrand_name: "NLO",
            original_m_uv: 20.0,
            shifted_m_uv: 7.0,
            original_mu_r: 3.0,
            shifted_mu_r: 9.0,
            skip_uv_profile: true,
            targets: IntegratedUvTargets {
                no_integrated: None,
                integrated: Some(Complex::new(F(0.0), F(1.163e-3))),
            },
            min_change_sigma: Some(5.0),
            min_mu_r_change_sigma: Some(5.0),
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
            ..Default::default()
        })
        .run(&mut cli.state, &cli.cli_settings)?;

        let uv = res.unwrap_uv();
        assert_eq!(uv.pass_fail(-0.9).failed, 0);
        clean_test(&state_path);
        Ok(())
    }
}
