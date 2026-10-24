//! Rebuild and reduce the complete physical five-loop source in a fresh process.
#[path = "../tests/fixtures/qcd_ladder/mod.rs"]
mod fixture;
use ::vakint;
use serde_json::json;
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::{BufWriter, Write},
    path::PathBuf,
    sync::Mutex,
    time::Instant,
};
use symbolica::atom::{Atom, AtomCore, AtomView};
use vakint::{
    EvaluationMethod, EvaluationOrder, OfficialExternalDataOptions, RustRedEvaluationOptions,
    RustRedMasterMode, Vakint, VakintSettings,
};

struct StageLogger(Mutex<BTreeMap<String, Vec<f64>>>);
static LOGGER: StageLogger = StageLogger(Mutex::new(BTreeMap::new()));
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
            && let Ok(mut stages) = self.0.lock()
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

fn main() -> fixture::Result<()> {
    let mut work = None;
    let mut artifacts = None;
    let mut cache = None;
    let mut prepare_only = false;
    let mut dress_only = false;
    let mut warm = false;
    let mut projection = fixture::Projection::QuadraticTransverse;
    let mut arguments = std::env::args_os().skip(1);
    while let Some(argument) = arguments.next() {
        match argument.to_str() {
            Some("--work-dir") => work = Some(PathBuf::from(arguments.next().ok_or("missing work directory")?)),
            Some("--artifacts") => artifacts = Some(PathBuf::from(arguments.next().ok_or("missing artifact path")?)),
            Some("--cache") => cache = Some(PathBuf::from(arguments.next().ok_or("missing cache path")?)),
            Some("--prepare-only") => prepare_only = true,
            Some("--dress-only") => dress_only = true,
            Some("--warm") => warm = true,
            Some("--projection") => projection = match arguments.next().and_then(|value| value.into_string().ok()).as_deref() {
                Some("quadratic") => fixture::Projection::QuadraticTransverse,
                Some("zero-momentum") => fixture::Projection::ZeroMomentumMetric,
                _ => return Err("projection must be quadratic or zero-momentum".into()),
            },
            _ => return Err("usage: five_loop_qcd_ladder --work-dir PATH [--artifacts PATH] [--projection quadratic|zero-momentum] [--prepare-only|--dress-only] [--cache PATH] [--warm]".into()),
        }
    }
    let work = work.ok_or("--work-dir is required")?;
    std::fs::create_dir_all(&work)?;
    log::set_logger(&LOGGER).map_err(|error| error.to_string())?;
    log::set_max_level(log::LevelFilter::Debug);
    if dress_only {
        let started = Instant::now();
        let generated = fixture::dress(&work)?;
        println!(
            "{}",
            json!({"event":"native_qcd_ladder_dressed", "seconds":started.elapsed().as_secs_f64(),
            "diagrams":generated.diagrams.len(), "loops":generated.diagrams[0].loop_count()})
        );
        return Ok(());
    }
    let total_started = Instant::now();
    let prepared = fixture::prepare(&work, projection)?;
    let input = prepared.input.to_canonical_string();
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
        "source_commit":git(&["rev-parse", "HEAD"]), "source_dirty":git(&["status", "--porcelain"]).map(|value| !value.is_empty()),
        "executable_sha256":format!("{:x}", Sha256::digest(std::fs::read(std::env::current_exe()?)?)),
        "fresh_native_preparation":true, "fixture_sources_sha256":source_hashes, "preparation":prepared.timings,
        "preparation_seconds":prepared.preparation_seconds, "scalar_terms":prepared.scalar_terms,
        "preparation_policy":{"fixture_work_count_overrides":false,
            "optional_numerator_denominator_cancellation":false,
            "optional_scalar_expansion":false,"scalar_term_count":"stored representation",
            "native_denominator_basis_lowering":"required by saved-index application"},
        "contributions":prepared.contributions, "input_sha256":format!("{:x}", Sha256::digest(input.as_bytes())),
        "input_bytes":input.len(), "prepared_nonzero_conditions":prepared.nonzero_conditions.iter().map(AtomCore::to_canonical_string).collect::<Vec<_>>(),
        "persistent_cache":cache, "reduction_completed":false});
    std::fs::write(
        work.join("benchmark.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    if prepare_only {
        println!("{report}");
        return Ok(());
    }
    let mut settings = VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/qcd-ladder-must-not-use-form".into(),
        use_dot_product_notation: true,
        ..Default::default()
    };
    // Match the ordinary five-loop evaluator's effective native policy, so
    // the recorded options describe the application that actually runs.
    settings.rustred_application_options.unresolved_policy =
        rustred::solver::UnresolvedPolicy::Preserve;
    report["conventions"] = json!({
        "integral_normalization":format!("{:?}",settings.integral_normalization_factor),
        "integral_normalization_expression":settings.get_integral_normalization_factor_atom()?.to_canonical_string(),
        "epsilon_symbol":settings.epsilon_symbol,
        "mass_convention":"one common denominator mass mUV; massless quark numerators",
        "application_options":format!("{:?}",settings.rustred_application_options),
        "scalar_limits":format!("{:?}",settings.rustred_scalar_numerator_limits),
        "matching_limits":format!("{:?}",settings.rustred_momentum_match_limits)});
    let started = Instant::now();
    let mut engine = Vakint::new()?;
    let loaded = engine.load_official_external_data(
        5,
        &OfficialExternalDataOptions {
            local_artifacts_path: artifacts,
            load_prepared_caches: false,
            ..Default::default()
        },
    )?;
    let caches_before = (1..=5)
        .map(|loops| Ok((loops, engine.cache_statistics(loops)?.attached)))
        .collect::<fixture::Result<BTreeMap<_, _>>>()?;
    if caches_before.values().any(|attached| *attached) {
        return Err("baseline must never import a persistent cache".into());
    }
    report["data_manifest_sha256"] = json!(format!(
        "{:x}",
        Sha256::digest(std::fs::read(loaded.path.join("manifest.json"))?)
    ));
    report["data"] = serde_json::to_value(loaded)?;
    report["data_loading_seconds"] = json!(started.elapsed().as_secs_f64());
    report["loading_stages"] = json!(*LOGGER.0.lock().map_err(|error| error.to_string())?);
    let owners_after_data_load = engine.rustred_application_owner_preparation_statistics()?;
    if let Some(cache) = &cache {
        let started = Instant::now();
        engine.load_cache(5, cache, &Default::default())?;
        report["cache_loading_seconds"] = json!(started.elapsed().as_secs_f64());
        engine.clear_cache(5)?;
    }
    let caches_during_evaluation = (1..=5)
        .map(|loops| Ok((loops, engine.cache_statistics(loops)?.attached)))
        .collect::<fixture::Result<BTreeMap<_, _>>>()?;
    if caches_during_evaluation
        .iter()
        .any(|(&loops, &attached)| attached != (cache.is_some() && loops == 5))
    {
        return Err(
            "measurement requires no persistent backing or only the manual five-loop file".into(),
        );
    }
    report["persistent_cache_attached"] = json!(caches_during_evaluation[&5]);
    report["persistent_cache_attachments_before"] = json!(caches_before);
    report["persistent_cache_attachments_during_evaluation"] = json!(caches_during_evaluation);
    report["persistent_cache_policy"] = json!(if cache.is_some() {
        "manual-five-loop-only"
    } else {
        "none"
    });
    let owners_before = engine.rustred_application_owner_preparation_statistics()?;
    LOGGER.0.lock().map_err(|error| error.to_string())?.clear();
    let started = Instant::now();
    let outcome = engine.evaluate_rustred_symbolic(&settings, prepared.input.as_view());
    report["reduction_seconds"] = json!(started.elapsed().as_secs_f64());
    report["end_to_end_seconds"] = json!(total_started.elapsed().as_secs_f64());
    report["warnings"] = json!(*WARNINGS.lock().map_err(|error| error.to_string())?);
    report["stages"] = json!(*LOGGER.0.lock().map_err(|error| error.to_string())?);
    report["cache_statistics"] = json!(format!(
        "{:?}",
        engine.rustred_application_cache_statistics()?
    ));
    report["witness_cache_statistics"] = json!(format!(
        "{:?}",
        engine.rustred_application_witness_cache_statistics()?
    ));
    report["saved_application_cache"] = engine
        .cache_statistics(5)?
        .saved_application
        .map(|statistics| statistics.to_document())
        .unwrap_or(serde_json::Value::Null);
    let owners_after = engine.rustred_application_owner_preparation_statistics()?;
    report["owner_preparation"] = json!({
        "declared_owners":owners_after.declared_owners,
        "prepared_owners":owners_after.prepared_owners,
        "failed_owners":owners_after.failed_owners,
        "preparation_attempts":owners_after.preparation_attempts,
        "prepared_owners_after_data_load":owners_after_data_load.prepared_owners,
        "prepared_owners_after_cache_setup":owners_before.prepared_owners,
        "preparation_attempts_after_data_load":owners_after_data_load.preparation_attempts,
        "preparation_attempts_after_cache_setup":owners_before.preparation_attempts,
        "preparation_during_load_seconds":owners_after_data_load.preparation_time.as_secs_f64(),
        "preparation_during_cache_setup_seconds":owners_before.preparation_time.saturating_sub(owners_after_data_load.preparation_time).as_secs_f64(),
        "preparation_during_reduction_seconds":owners_after.preparation_time.saturating_sub(owners_before.preparation_time).as_secs_f64()
    });
    report["native_timing_scope"] = json!({
        "saved_rule_application":"Native entry transport, lazy owner preparation, DAG reduction, exact sum, and terminal substitution",
        "owner_preparation":"Measured cumulative deferred-owner preparation counters at data load, cache setup, and reduction boundaries",
        "pure_ibp_seconds":null,"terminal_substitution_seconds":null,
        "separate_native_phase_durations_available":false,
        "phase_duration_estimates_by_subtraction":false
    });
    let result = match outcome {
        Ok(result) => result,
        Err(error) => {
            report["reduction_error"] = json!(error.to_string());
            std::fs::write(
                work.join("benchmark.json"),
                serde_json::to_vec_pretty(&report)?,
            )?;
            return Err(error.into());
        }
    };
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
    let value = result.value.to_canonical_string();
    report["result_sha256"] = json!(format!("{:x}", Sha256::digest(value.as_bytes())));
    report["result_bytes"] = json!(value.len());
    let d = rustred::algebra::CoefficientContext::try_new(["d"])?
        .parameter("d")
        .ok_or("missing native dimension")?
        .to_expression();
    let physical_d = Atom::num(4)
        - Atom::num(2) * Atom::var(vakint::vakint_symbol!(settings.epsilon_symbol.as_str()));
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
    }
    conditions_output.flush()?;
    report["guards"] = json!({"count":complete_conditions.len(),
        "sha256":format!("{:x}",conditions_digest.finalize()),
        "file":"nonzero-conditions.symbolica", "preview":conditions_preview,
        "includes_preparation_conditions":true});
    report["families"] = json!(
        result
            .families
            .iter()
            .map(|(symbol, identity)| (symbol.to_canonical_string(), identity))
            .collect::<BTreeMap<_, _>>()
    );
    std::fs::write(work.join("result.symbolica"), value)?;
    report["retained_result_available"] = json!(true);
    report["reduction_completed"] = json!(unknown.is_empty());
    // Include cold-result classification and artifact serialization, while
    // keeping the optional warm repetition outside the cold total.
    report["end_to_end_seconds"] = json!(total_started.elapsed().as_secs_f64());
    if warm {
        LOGGER.0.lock().map_err(|error| error.to_string())?.clear();
        let started = Instant::now();
        let second = engine.evaluate_rustred_symbolic(&settings, prepared.input.as_view())?;
        report["warm_reduction_seconds"] = json!(started.elapsed().as_secs_f64());
        report["warm_stages"] = json!(*LOGGER.0.lock().map_err(|error| error.to_string())?);
        if !(result.value - second.value).together().expand().is_zero()
            || result.nonzero_conditions != second.nonzero_conditions
            || result.families != second.families
        {
            return Err("warm reduction failed exact result/guard/family parity".into());
        }
        report["warm_exact_parity"] = json!(true);
    }
    std::fs::write(
        work.join("benchmark.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    println!("{report}");
    Ok(())
}
