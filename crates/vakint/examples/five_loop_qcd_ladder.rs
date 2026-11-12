//! Prepare the complete physical five-loop source, or reuse its verified scalar input.
//!
//! The saved zero-momentum result can be inspected or expanded through epsilon^5
//! without repeating IBPs using `five_loop_q0_coefficients.py`. Its master products
//! stay formal; see `docs/five_loop_q0_coefficients.typ` for the precision contract.
#[path = "../tests/fixtures/qcd_ladder/mod.rs"]
mod fixture;
#[path = "five_loop_qcd_ladder/resource_guard_admission.rs"]
mod resource_guard_admission;
use ::vakint;
use serde_json::json;
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::{BufWriter, Write},
    path::PathBuf,
    sync::Mutex,
    time::{Duration, Instant},
};
use symbolica::atom::{Atom, AtomCore, AtomView};
use vakint::{
    EvaluationMethod, EvaluationOrder, OfficialExternalDataOptions, RustRedEvaluationOptions,
    RustRedMasterMode, RustRedMonitorDisplay, RustRedMonitorOptions, RustRedMonitorScope, Vakint,
    VakintSettings,
};

const USAGE: &str = "usage: five_loop_qcd_ladder --work-dir PATH [--artifacts PATH] [--projection quadratic|zero-momentum] [--prepared-input PATH|--fresh-prepare] [--prepare-only|--dress-only] [--warm] [--monitor] [--progress-snapshot PATH] [--progress-interval-seconds N]\nThe reduction uses RustRed's ordered scheduler with the weight arithmetic of VAKINT_ORDERED_WEIGHT_ARITHMETIC (Vakint's default when unset, which five-loop requests run as modular_sampled) on VAKINT_N_CORES cores (Vakint's default when unset). --warm repeats the evaluation in the same engine and checks exact parity; there is no result cache.\nMonitoring defaults off. --monitor requires progress-tui and defaults its snapshot to work-dir/monitor.json. A snapshot path enables native JSON without a TUI; interval defaults to 1s and must be finite and nonnegative. Monitored timings include diagnostic observation overhead.";

struct StageLogger {
    phases: Mutex<BTreeMap<String, Vec<f64>>>,
}

impl StageLogger {
    fn stage_accounting() -> serde_json::Value {
        json!({
            "samples":"Every stages[name] entry is one completed span; repeated calls append samples. Sum samples within a phase, never all phases together.",
            "disjoint_evaluation_spans":[
                "frontend_setup","frontend_canonicalization","frontend_tensor_prepass",
                "scalar_input_split","frontend_admitted_witness_search","numerator_replacement",
                "ordered_batch_total","retained_finalization_total","public_result_packing"
            ],
            "nested_spans":{
                "ordered_batch_total":[
                    "parent_transport_and_scalar_lowering","coverage_mass_batching",
                    "lowering_guard_materialization","saved_rule_application",
                    "terminal_geometry_and_lower_factors","retained_guard_materialization",
                    "retained_coefficient_materialization","unprojected_terminal_materialization",
                    "retained_batch_mass_assembly","retained_batch_cleanup",
                    "retained_normalization_and_notation","detailed_result_packing"
                ],
                "retained_finalization_total":[
                    "retained_catalogue_exact_combination","retained_catalogue_classification"
                ]
            },
            "unassigned":"Public evaluation minus disjoint spans and ordered-batch/finalization totals minus their own children remain raw uninstrumented control, logging and caller cleanup time. Report each residual explicitly; no residual is a pure IBP estimate.",
            "saved_rule_application_scope":"The complete saved call includes entry transport, the ordered native reduction with weighted exact assembly, and saved terminal substitution.",
            "loading_scope":"Initial loading_stages are captured before evaluation. Preparation spans may recur inside frontend_setup when changed application options reload the package; products_prepared includes family_verified. Never add these overlapping spans to measured loading or evaluation wall totals.",
            "path_scope":"Successful and fallback batch spans append independently. Unprojected terminal materialization is a fallback path; cleanup is measured after existing batch locals naturally drop. Failed operations may omit completed spans."
        })
    }
}

static LOGGER: StageLogger = StageLogger {
    phases: Mutex::new(BTreeMap::new()),
};
static WARNINGS: Mutex<Vec<String>> = Mutex::new(Vec::new());
impl log::Log for StageLogger {
    fn enabled(&self, metadata: &log::Metadata<'_>) -> bool {
        metadata.level() <= log::Level::Warn || metadata.target() == "vakint::saved_application"
    }
    fn log(&self, record: &log::Record<'_>) {
        if !self.enabled(record.metadata()) {
            return;
        }
        let message = record.args().to_string();
        if record.level() <= log::Level::Warn {
            eprintln!("{message}");
            if let Ok(mut warnings) = WARNINGS.lock() {
                warnings.push(message.clone());
            }
        }
        let mut phase = None;
        let mut elapsed = None;
        for field in message.split_whitespace() {
            if let Some(value) = field.strip_prefix("phase=") {
                phase = Some(value);
            }
            if let Some(value) = field.strip_prefix("elapsed_seconds=") {
                elapsed = value.parse::<f64>().ok();
            }
        }
        if let (Some(phase), Some(elapsed)) = (phase, elapsed)
            && let Ok(mut stages) = self.phases.lock()
        {
            stages.entry(phase.into()).or_default().push(elapsed);
        }
    }
    fn flush(&self) {}
}

fn git(arguments: &[&str]) -> Option<String> {
    let output = std::process::Command::new("git")
        .args(arguments)
        .current_dir(concat!(env!("CARGO_MANIFEST_DIR"), "/../.."))
        .output()
        .ok()?;
    output
        .status
        .success()
        .then(|| String::from_utf8(output.stdout).ok())
        .flatten()
        .map(|value| value.trim().to_owned())
}

fn collector_helper_sources() -> BTreeMap<&'static str, String> {
    [(
        "crates/vakint/examples/five_loop_qcd_ladder/resource_guard_admission.rs",
        include_bytes!("five_loop_qcd_ladder/resource_guard_admission.rs").as_slice(),
    )]
    .into_iter()
    .map(|(path, bytes)| (path, format!("{:x}", Sha256::digest(bytes))))
    .collect()
}

fn main() -> fixture::Result<()> {
    let mut work = None;
    let mut artifacts = None;
    let mut prepared_input = None;
    let mut prepare_only = false;
    let mut fresh_prepare = false;
    let mut dress_only = false;
    let mut warm = false;
    let mut monitor = RustRedMonitorOptions::default();
    let mut projection = fixture::Projection::QuadraticTransverse;
    let mut arguments = std::env::args_os().skip(1);
    while let Some(argument) = arguments.next() {
        match argument.to_str() {
            Some("--help" | "-h") => {
                println!("{USAGE}");
                return Ok(());
            }
            Some("--work-dir") => {
                work = Some(PathBuf::from(
                    arguments.next().ok_or("missing work directory")?,
                ))
            }
            Some("--artifacts") => {
                artifacts = Some(PathBuf::from(
                    arguments.next().ok_or("missing artifact path")?,
                ))
            }
            Some("--prepared-input") => {
                prepared_input = Some(PathBuf::from(
                    arguments.next().ok_or("missing prepared-input directory")?,
                ))
            }
            Some("--prepare-only") => prepare_only = true,
            Some("--fresh-prepare") => fresh_prepare = true,
            Some("--dress-only") => dress_only = true,
            Some("--warm") => warm = true,
            Some("--monitor") => monitor.display = RustRedMonitorDisplay::Tui,
            Some("--progress-snapshot") => {
                monitor.progress_snapshot = Some(PathBuf::from(
                    arguments.next().ok_or("missing progress snapshot path")?,
                ));
            }
            Some("--progress-interval-seconds") => {
                let seconds = arguments
                    .next()
                    .ok_or("missing progress interval")?
                    .into_string()
                    .map_err(|_| "progress interval must be UTF-8")?
                    .parse::<f64>()?;
                monitor.interval = Duration::try_from_secs_f64(seconds)?;
            }

            Some("--projection") => {
                projection = match arguments
                    .next()
                    .and_then(|value| value.into_string().ok())
                    .as_deref()
                {
                    Some("quadratic") => fixture::Projection::QuadraticTransverse,
                    Some("zero-momentum") => fixture::Projection::ZeroMomentumMetric,
                    _ => return Err("projection must be quadratic or zero-momentum".into()),
                }
            }
            _ => return Err(USAGE.into()),
        }
    }
    if fresh_prepare && prepared_input.is_some() {
        return Err("--fresh-prepare cannot be combined with --prepared-input".into());
    }
    if dress_only && prepared_input.is_some() {
        return Err("--dress-only cannot be combined with --prepared-input".into());
    }
    let arithmetic = fixture::comparison::arithmetic_from_environment()?;
    let work = work.ok_or("--work-dir is required")?;
    if monitor.display == RustRedMonitorDisplay::Tui && monitor.progress_snapshot.is_none() {
        monitor.progress_snapshot = Some(work.join("monitor.json"));
    }
    monitor.validate()?;

    let _guard_admission = resource_guard_admission::admit_from_task_env()?;

    std::fs::create_dir_all(&work)?;
    let write_report = |report: &serde_json::Value| -> fixture::Result<()> {
        std::fs::write(
            work.join("benchmark.json"),
            serde_json::to_vec_pretty(report)?,
        )?;
        Ok(())
    };
    log::set_logger(&LOGGER).map_err(|error| error.to_string())?;
    log::set_max_level(log::LevelFilter::Debug);
    let operation = |mut scope: Option<&mut RustRedMonitorScope<'_, '_>>| -> fixture::Result<serde_json::Value> {
        if dress_only {
            if let Some(scope) = scope.as_deref_mut() { scope.phase("diagram_dressing"); }
            let started = Instant::now();
            let generated = fixture::dress(&work)?;
            return Ok(json!({"event":"native_qcd_ladder_dressed", "seconds":started.elapsed().as_secs_f64(),
                "diagrams":generated.diagrams.len(), "loops":generated.diagrams[0].loop_count()}));
        }
        let total_started = Instant::now();
        if let Some(scope) = scope.as_deref_mut() { scope.phase("prepared_input_selection"); }
        if prepared_input.is_none() && fixture::LadderFixture::has_prepared_input(&work)? {
            if fresh_prepare {
                return Err("--fresh-prepare requires a new --work-dir; preserve the completed prepared-input bundle".into());
            }
            prepared_input = Some(work.clone());
        }
        let fresh_native_preparation = prepared_input.is_none();
        let input_started = Instant::now();
        if let Some(scope) = scope.as_deref_mut() { scope.phase(if fresh_native_preparation { "scalar_input_preparation" } else { "prepared_input_loading" }); }
        let prepared = match &prepared_input {
            Some(directory) => fixture::LadderFixture::load_prepared(directory, projection)?,
            None => fixture::prepare(&work, projection)?,
        };
        let input_loading_seconds =
            (!fresh_native_preparation).then(|| input_started.elapsed().as_secs_f64());
        if let Some(scope) = scope.as_deref_mut() { scope.phase("prepared_input_authentication_and_report"); }
        let manifest_directory = prepared_input.as_ref().unwrap_or(&work);
        let manifest = match &prepared.prepared_input {
            Some(manifest) => manifest.clone(),
            None => prepared.save_prepared(&work)?,
        };
        let manifest_path =
            std::fs::canonicalize(manifest_directory.join(fixture::PREPARED_INPUT_MANIFEST))?;
        let sources: [(&str, &[u8]); 7] = [
            (
                "mod.rs",
                include_bytes!("../tests/fixtures/qcd_ladder/mod.rs"),
            ),
            (
                "physical.rs",
                include_bytes!("../tests/fixtures/qcd_ladder/physical.rs"),
            ),
            (
                "projector.rs",
                include_bytes!("../tests/fixtures/qcd_ladder/projector.rs"),
            ),
            (
                "vacuum.rs",
                include_bytes!("../tests/fixtures/qcd_ladder/vacuum.rs"),
            ),
            (
                "ladder.dot",
                include_bytes!("../tests/fixtures/qcd_ladder/ladder.dot"),
            ),
            ("qcd.json", include_bytes!("../tests/fixtures/qcd.json")),
            ("phi3.json", include_bytes!("../tests/fixtures/phi3.json")),
        ];
        let source_hashes = sources
            .into_iter()
            .map(|(name, bytes)| (name, format!("{:x}", Sha256::digest(bytes))))
            .collect::<BTreeMap<_, _>>();
        let mut report = json!({"schema":"vakint.qcd-ladder.benchmark.v1", "projection":prepared.projection,
            "reduction_strategy":"Ordered","ordered_weight_arithmetic":arithmetic.to_string(),
            "source_commit":git(&["rev-parse", "HEAD"]), "source_dirty":git(&["status", "--porcelain"]).map(|value| !value.is_empty()),
            "compiled_cohort":{
                "main_revision":option_env!("VAKINT_CACHE_BUILD_GIT_REVISION"),
                "native_revision":option_env!("VAKINT_CACHE_BUILD_RUSTRED_REVISION"),
                "kernel_revision":option_env!("VAKINT_CACHE_BUILD_SYMBOLICA_REVISION"),
                "main_source_snapshot_sha256":option_env!("VAKINT_CACHE_BUILD_MAIN_SOURCE_SNAPSHOT_SHA256"),
                "native_source_snapshot_sha256":option_env!("VAKINT_CACHE_BUILD_NATIVE_SOURCE_SNAPSHOT_SHA256"),
                "rustc":option_env!("VAKINT_CACHE_BUILD_RUSTC"),
                "profile":option_env!("VAKINT_CACHE_BUILD_PROFILE"),
            },
            "executable_sha256":format!("{:x}", Sha256::digest(std::fs::read(std::env::current_exe()?)?)),
            "collector_source_sha256":format!("{:x}", Sha256::digest(include_bytes!("five_loop_qcd_ladder.rs"))),
            "collector_helper_sources_sha256":collector_helper_sources(),
            "master_mode":"retain",
            "guard_policy":"canonical-full-sorted-deduplicated",
            "performance_claim":false,
            "warm_requested":warm,
            "input_canonical_sha256":format!("{:x}",Sha256::digest(prepared.input.to_canonical_string().as_bytes())),
            "fresh_native_preparation":fresh_native_preparation,
            "input_mode":if fresh_native_preparation {"fresh"} else {"prepared"},
            "fixture_sources_sha256":manifest.producer.fixture_sources_sha256,
            "consumer_fixture_sources_sha256":source_hashes,
            "prepared_input":{"manifest_path":manifest_path,
                "manifest_sha256":format!("{:x}",Sha256::digest(std::fs::read(&manifest_path)?)),
                "manifest":manifest},
            "input_loading_seconds":input_loading_seconds,
            "preparation":if fresh_native_preparation {serde_json::to_value(&prepared.timings)?} else {serde_json::Value::Null},
            "preparation_seconds":fresh_native_preparation.then_some(prepared.preparation_seconds),
            "end_to_end_scope":if fresh_native_preparation {"fresh native preparation through cold result serialization"} else {"verified prepared-input loading through cold result serialization; original preparation excluded"},
            "scalar_terms":prepared.scalar_terms,
            "preparation_policy":{"fixture_work_count_overrides":false,
                "optional_numerator_denominator_cancellation":false,
                "optional_scalar_expansion":false,"scalar_term_count":"stored representation",
                "native_denominator_basis_lowering":"required by saved-index application"},
            "contributions":prepared.contributions, "input_sha256":manifest.input.sha256,
            "input_bytes":manifest.input.bytes, "prepared_nonzero_conditions":prepared.nonzero_conditions.iter().map(AtomCore::to_canonical_string).collect::<Vec<_>>(),
            "persistent_cache":null,
            "reduction_completed":false});
        if monitor.enabled() {
            report["monitor"] = json!({
                "display": match monitor.display {
                    RustRedMonitorDisplay::Off => "off",
                    RustRedMonitorDisplay::Tui => "tui",
                    RustRedMonitorDisplay::Json => "json",
                },
                "progress_snapshot": monitor.progress_snapshot,
                "interval_seconds": monitor.interval.as_secs_f64(),
                "diagnostic_timing": true,
                "scope":"job", "api_completion":"separate from durable outputs",
            });
        }
        if let Some(scope) = scope.as_deref_mut() { scope.phase("durable_report_write"); }
        write_report(&report)?;
        if prepare_only {
            return Ok(report);
        }
        let mut settings = VakintSettings {
            evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
                RustRedEvaluationOptions {
                    master_mode: RustRedMasterMode::Retain,
                },
            )]),
            integral_normalization_factor: vakint::LoopNormalizationFactor::MSbar,
            form_exe_path: "/qcd-ladder-must-not-use-form".into(),
            use_dot_product_notation: true,
            ..Default::default()
        };
        // Match the ordinary five-loop evaluator's effective native policy, so
        // the recorded options describe the application that actually runs.
        settings.rustred_application_options.unresolved_policy =
            rustred::solver::UnresolvedPolicy::Preserve;
        settings.ordered_weight_arithmetic = arithmetic;
        settings.n_cores = fixture::comparison::n_cores_from_environment()?;
        // The arithmetic the five-loop requests run (the default selection
        // runs modular_sampled), next to the selection it comes from.
        report["ordered_weight_arithmetic"] =
            json!(settings.rustred_five_loop_ordered_weight_arithmetic().to_string());
        report["ordered_weight_arithmetic_selection"] = json!(arithmetic.to_string());
        report["conventions"] = json!({
            "integral_normalization":format!("{:?}",settings.integral_normalization_factor),
            "integral_normalization_expression":settings.get_integral_normalization_factor_atom()?.to_canonical_string(),
            "epsilon_symbol":settings.epsilon_symbol,
            "mass_convention":"one common denominator mass mUV; massless quark numerators",
            "application_options":format!("{:?}",settings.rustred_application_options),
            "scalar_limits":format!("{:?}",settings.rustred_scalar_numerator_limits),
            "matching_limits":format!("{:?}",settings.rustred_momentum_match_limits)});
        let started = Instant::now();
        if let Some(scope) = scope.as_deref_mut() { scope.phase("external_data_loading"); }
        let mut engine = Vakint::new()?;
        let loaded = engine.load_official_external_data(
            5,
            &OfficialExternalDataOptions {
                local_artifacts_path: artifacts,
                n_cores: settings.n_cores,
                ..Default::default()
            },
        )?;
        report["data_manifest_sha256"] = json!(format!(
            "{:x}",
            Sha256::digest(std::fs::read(loaded.path.join("manifest.json"))?)
        ));
        report["data"] = serde_json::to_value(loaded)?;
        report["data_loading_seconds"] = json!(started.elapsed().as_secs_f64());
        report["loading_stages"] = json!(*LOGGER.phases.lock().map_err(|error| error.to_string())?);
        let owners_after_data_load = engine.rustred_application_owner_preparation_statistics()?;
        let evaluation = (|| -> fixture::Result<()> {
            let owners_before = engine.rustred_application_owner_preparation_statistics()?;
            LOGGER
                .phases
                .lock()
                .map_err(|error| error.to_string())?
                .clear();
            report["execution_threads"] = json!({"math":1,
            "scope":"dedicated collector math; not measured OS thread count"});
            let started = Instant::now();
            let outcome = if let Some(scope) = scope.as_deref_mut() {
                engine.evaluate_rustred_symbolic_in_scope(&settings, prepared.input.as_view(), scope)
            } else if monitor.enabled() {
                engine.evaluate_rustred_symbolic_with_monitor(
                    &settings,
                    prepared.input.as_view(),
                    &monitor,
                )
            } else {
                engine.evaluate_rustred_symbolic(&settings, prepared.input.as_view())
            };
            report["reduction_seconds"] = json!(started.elapsed().as_secs_f64());
            report["end_to_end_seconds"] = json!(total_started.elapsed().as_secs_f64());
            let result = match outcome {
                Ok(result) => result,
                Err(error) => {
                    report["reduction_error"] = json!(error.to_string());
                    report["warnings"] = json!(WARNINGS.lock().ok().map(|warnings| warnings.clone()));
                    report["stages"] = json!(LOGGER.phases.lock().ok().map(|stages| stages.clone()));
                    if let Some(scope) = scope.as_deref_mut() { scope.phase("durable_report_write"); }
                    if let Err(receipt_error) = write_report(&report) {
                        eprintln!("could not write failed reduction receipt: {receipt_error}");
                    }
                    return Err(error.into());
                }
            };
            if let Some(scope) = scope.as_deref_mut() { scope.phase("evaluation_statistics"); }
            report["warnings"] = json!(*WARNINGS.lock().map_err(|error| error.to_string())?);
            report["stages"] = json!(*LOGGER.phases.lock().map_err(|error| error.to_string())?);
            report["witness_cache_statistics"] = json!(format!(
                "{:?}",
                engine.rustred_application_witness_cache_statistics()?
            ));
            let owners_after = engine.rustred_application_owner_preparation_statistics()?;
            report["owner_preparation"] = json!({
                "declared_owners":owners_after.declared_owners,
                "prepared_owners":owners_after.prepared_owners,
                "failed_owners":owners_after.failed_owners,
                "preparation_attempts":owners_after.preparation_attempts,
                "prepared_owners_after_data_load":owners_after_data_load.prepared_owners,
                "preparation_attempts_after_data_load":owners_after_data_load.preparation_attempts,
                "preparation_during_load_seconds":owners_after_data_load.preparation_time.as_secs_f64(),
                "preparation_during_reduction_seconds":owners_after.preparation_time.saturating_sub(owners_before.preparation_time).as_secs_f64()
            });
            report["native_timing_scope"] = json!({
                "saved_rule_application":"Native entry transport, lazy owner preparation, DAG reduction, exact sum, and terminal substitution",
                "owner_preparation":"Measured cumulative deferred-owner preparation counters at data load and reduction boundaries",
                "pure_ibp_seconds":null,"terminal_substitution_seconds":null,
                "separate_native_phase_durations_available":false,
                "phase_duration_estimates_by_subtraction":false,
                "stage_accounting":StageLogger::stage_accounting()
            });
            if let Some(scope) = scope.as_deref_mut() { scope.phase("result_catalogue_classification"); }
            let known = engine
                .master_symbols()?
                .into_iter()
                .map(|symbol| symbol.as_atom().clone())
                .collect::<BTreeSet<_>>();
            let mut masters = BTreeSet::new();
            result.value.visitor(&mut |part| {
            if matches!(part, AtomView::Fun(function) if function.get_symbol() == vakint::vakint_symbol!("RustRedMaster")) {
                masters.insert(part.to_owned());
                return false;
            }
            true
        });
            let unknown = masters
                .difference(&known)
                .map(AtomCore::to_canonical_string)
                .collect::<BTreeSet<_>>();
            let unknown_bytes = serde_json::to_vec_pretty(&unknown)?;
            std::fs::write(work.join("unknown-masters.json"), &unknown_bytes)?;
            report["catalogue_closure"] = json!({"closed":unknown.is_empty(),"master_count":masters.len(),
            "unknown_count":unknown.len(),"unknown_file":"unknown-masters.json",
            "unknown_sha256":format!("{:x}",Sha256::digest(&unknown_bytes)),"unknown_preview":unknown.iter().take(8).collect::<Vec<_>>()});
            if let Some(scope) = scope.as_deref_mut() { scope.expression_phase("result_expression_serialization", result.value.as_view()); }
            let value = result.value.to_canonical_string();
            report["result_sha256"] = json!(format!("{:x}", Sha256::digest(value.as_bytes())));
            report["result_bytes"] = json!(value.len());
            let d = rustred::algebra::CoefficientContext::try_new(["d"])?
                .parameter("d")
                .ok_or("missing native dimension")?
                .to_expression();
            let physical_d = Atom::num(4)
                - Atom::num(2) * Atom::var(vakint::vakint_symbol!(settings.epsilon_symbol.as_str()));
            if let Some(scope) = scope.as_deref_mut() { scope.phase("result_guard_combination"); }
            let complete_conditions = result
                .nonzero_conditions
                .iter()
                .cloned()
                .chain(
                    prepared
                        .nonzero_conditions
                        .iter()
                        .map(|condition| condition.replace(d.clone()).with(physical_d.clone())),
                )
                .collect::<BTreeSet<_>>();
            if complete_conditions.iter().any(AtomCore::is_zero) {
                return Err("a complete nonzero condition vanished".into());
            }
            let conditions_path = work.join("nonzero-conditions.symbolica");
            let mut conditions_output = BufWriter::new(std::fs::File::create(&conditions_path)?);
            let mut conditions_digest = Sha256::new();
            let mut conditions_preview = Vec::new();
            {
            let mut progress = scope.as_deref_mut().map(|scope| scope.work_phase("result_guard_serialization", complete_conditions.len(), "guards serialized"));
            for (index, condition) in complete_conditions.iter().enumerate() {
                let text = condition.to_canonical_string();
                conditions_output.write_all(text.as_bytes())?;
                conditions_output.write_all(b"\n")?;
                conditions_digest.update(text.as_bytes());
                conditions_digest.update(b"\n");
                if index < 3 {
                    let mut preview = text.chars().take(180).collect::<String>();
                    if text.chars().count() > 180 {
                        preview.push('…');
                    }
                    conditions_preview.push(preview);
                }
                if let Some(progress) = &mut progress { progress.completed(index + 1); }
            }
            }
            if let Some(scope) = scope.as_deref_mut() { scope.phase("result_guard_flush"); }
            conditions_output.flush()?;
            report["guards"] = json!({"count":complete_conditions.len(),
            "sha256":format!("{:x}",conditions_digest.finalize()),
            "file":"nonzero-conditions.symbolica", "preview":conditions_preview,
            "includes_preparation_conditions":true});
            if let Some(scope) = scope.as_deref_mut() { scope.phase("result_family_serialization"); }
            report["families"] = json!(
                result
                    .families
                    .iter()
                    .map(|(symbol, identity)| (symbol.to_canonical_string(), identity))
                    .collect::<BTreeMap<_, _>>()
            );
            if let Some(scope) = scope.as_deref_mut() { scope.phase("result_expression_write"); }
            std::fs::write(work.join("result.symbolica"), value)?;
            report["retained_result_available"] = json!(true);
            report["reduction_completed"] = json!(unknown.is_empty());
            // Include cold-result classification and artifact serialization, while
            // keeping the optional warm repetition outside the cold total.
            report["end_to_end_seconds"] = json!(total_started.elapsed().as_secs_f64());
            if warm {
                LOGGER
                    .phases
                    .lock()
                    .map_err(|error| error.to_string())?
                    .clear();
                let started = Instant::now();
                let outcome = if let Some(scope) = scope.as_deref_mut() {
                    engine.evaluate_rustred_symbolic_in_scope(&settings, prepared.input.as_view(), scope)
                } else if monitor.enabled() {
                    engine.evaluate_rustred_symbolic_with_monitor(
                        &settings,
                        prepared.input.as_view(),
                        &monitor,
                    )
                } else {
                    engine.evaluate_rustred_symbolic(&settings, prepared.input.as_view())
                };
                report["warm_reduction_seconds"] = json!(started.elapsed().as_secs_f64());
                let second = match outcome {
                    Ok(result) => result,
                    Err(error) => {
                        report["warm_reduction_error"] = json!(error.to_string());
                        if let Some(scope) = scope.as_deref_mut() { scope.phase("durable_report_write"); }
                        if let Err(receipt_error) = write_report(&report) {
                            eprintln!("could not write failed warm reduction receipt: {receipt_error}");
                        }
                        return Err(error.into());
                    }
                };
                report["warm_stages"] =
                    json!(*LOGGER.phases.lock().map_err(|error| error.to_string())?);
                if let Some(scope) = scope.as_deref_mut() { scope.phase("warm_exact_parity"); }
                if !(result.value - second.value).together().expand().is_zero()
                    || result.nonzero_conditions != second.nonzero_conditions
                    || result.families != second.families
                {
                    return Err("warm reduction failed exact result/guard/family parity".into());
                }
                report["warm_exact_parity"] = json!(true);
            }
            Ok(())
        })();
        if let Err(error) = evaluation {
            if let Some(scope) = scope.as_deref_mut() { scope.phase("durable_report_write"); }
            if let Err(receipt_error) = write_report(&report) {
                eprintln!("could not write terminal reduction receipt: {receipt_error}");
            }
            return Err(error);
        }
        if let Some(scope) = scope { scope.phase("final_durable_report"); }
        write_report(&report)?;
        Ok(report)
    };
    let output = RustRedMonitorScope::run(&monitor, |scope| operation(Some(scope)))?;
    // Deliver stdout after the TUI is restored; all durable artifacts are done.
    println!("{output}");
    Ok(())
}

#[cfg(test)]
mod stage_logger_tests {
    use super::*;
    use log::Log;

    #[test]
    fn physical_logger_keeps_repeated_nested_spans_separate_from_parent_totals() {
        let logger = StageLogger {
            phases: Mutex::new(BTreeMap::new()),
        };
        for (phase, seconds) in [
            ("frontend_setup", 1.0),
            ("coverage_mass_batching", 2.0),
            ("coverage_mass_batching", 3.0),
            ("retained_guard_materialization", 2.0),
            ("retained_guard_materialization", 3.0),
            ("retained_batch_cleanup", 0.5),
            ("ordered_batch_total", 12.0),
            ("retained_finalization_total", 5.0),
            ("retained_catalogue_exact_combination", 3.0),
            ("retained_catalogue_classification", 1.0),
            ("public_result_packing", 0.5),
        ] {
            logger.log(
                &log::Record::builder()
                    .level(log::Level::Debug)
                    .target("vakint::saved_application")
                    .args(format_args!("phase={phase} elapsed_seconds={seconds}"))
                    .build(),
            );
        }
        let phases = logger.phases.lock().unwrap();
        assert_eq!(phases["coverage_mass_batching"], [2.0, 3.0]);
        assert_eq!(phases["retained_guard_materialization"], [2.0, 3.0]);
        assert_eq!(phases["ordered_batch_total"], [12.0]);
        let accounting = StageLogger::stage_accounting();
        let disjoint = accounting["disjoint_evaluation_spans"]
            .as_array()
            .unwrap()
            .iter()
            .map(|phase| phase.as_str().unwrap())
            .collect::<BTreeSet<_>>();
        assert!(!disjoint.contains("retained_outer_exact_combination"));
        assert_eq!(
            accounting["nested_spans"]["retained_finalization_total"],
            json!([
                "retained_catalogue_exact_combination",
                "retained_catalogue_classification"
            ])
        );
        assert_eq!(
            disjoint.len(),
            accounting["disjoint_evaluation_spans"]
                .as_array()
                .unwrap()
                .len()
        );
        let mut nested = BTreeSet::new();
        for (parent, children) in accounting["nested_spans"].as_object().unwrap() {
            assert!(disjoint.contains(parent.as_str()));
            for child in children.as_array().unwrap() {
                let child = child.as_str().unwrap();
                assert!(!disjoint.contains(child));
                assert!(nested.insert(child));
            }
        }
        let accounted = disjoint
            .iter()
            .filter_map(|phase| phases.get(*phase))
            .flatten()
            .sum::<f64>();
        assert_eq!(accounted, 18.5);
        assert!(phases.values().flatten().sum::<f64>() > accounted);
        assert!(
            accounting["unassigned"]
                .as_str()
                .unwrap()
                .contains("explicitly")
        );
    }
}
