//! Prepare the complete physical five-loop source, or reuse its verified scalar input.
#[path = "../tests/fixtures/qcd_ladder/mod.rs"]
mod fixture;
use ::vakint;
use rustred::solver::CandidateCacheRepresentation;
use serde_json::json;
use sha2::{Digest, Sha256};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File, OpenOptions},
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    sync::{Mutex, mpsc},
    thread::{self, JoinHandle},
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};
use symbolica::atom::{Atom, AtomCore, AtomView};
use vakint::{
    EvaluationMethod, EvaluationOrder, OfficialExternalDataOptions, RustRedEvaluationOptions,
    RustRedMasterMode, RustRedMonitorDisplay, RustRedMonitorOptions, RustRedRunningCacheControl,
    RustRedRunningCacheOptions, Vakint, VakintSettings,
};

const USAGE: &str = "usage: five_loop_qcd_ladder --work-dir PATH [--artifacts PATH] [--projection quadratic|zero-momentum] [--coefficient-representation sparse|factorized] [--prepared-input PATH|--fresh-prepare] [--prepare-only|--dress-only] [--cache PATH] [--running-cache PATH] [--running-cache-capture-only] [--running-cache-no-read-delta] [--warm] [--monitor] [--progress-snapshot PATH] [--progress-interval-seconds N]\nMonitoring defaults off. --monitor requires progress-tui and defaults its snapshot to work-dir/monitor.json. A snapshot path enables native JSON without a TUI; interval defaults to 1s and must be finite and nonnegative. Monitored timings include diagnostic observation overhead.";

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

struct RunningTelemetry {
    stop: mpsc::Sender<()>,
    worker: Option<JoinHandle<serde_json::Value>>,
}

impl RunningTelemetry {
    fn start(
        control: RustRedRunningCacheControl,
        work: PathBuf,
        phase: &'static str,
    ) -> std::io::Result<Self> {
        let (stop, receiver) = mpsc::channel();
        let worker = thread::Builder::new()
            .name("running-cache-telemetry".into())
            .spawn(move || {
                let session = format!(
                    "{}-{phase}-{}",
                    std::process::id(),
                    SystemTime::now().duration_since(UNIX_EPOCH).unwrap_or_default().as_nanos()
                );
                let mut reports = 0_u64;
                let mut errors = 0_u64;
                let mut last_good_utc = None;
                let mut last_error = None;
                loop {
                    let utc = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .ok()
                        .map(|time| time.as_secs_f64());
                    let result = (|| -> Result<(), String> {
                        let native = control.report().map_err(|error| error.to_string())?;
                        let document = json!({
                            "schema":"vakint.running-cache.telemetry.v1",
                            "utc_unix_seconds":utc,"pid":std::process::id(),
                            "session":session,"phase":phase,"interval_seconds":10,
                            "native":native.to_document(),
                            "record_semantics":"accepted is enqueued; committed is acknowledged durable",
                        });
                        Self::write_atomic(&work.join("running-cache-live.json"), &document)
                            .map_err(|error| error.to_string())?;
                        last_good_utc = utc;
                        reports = reports.saturating_add(1);
                        Ok(())
                    })();
                    if let Err(error) = result {
                        errors = errors.saturating_add(1);
                        last_error = Some(error.clone());
                        eprintln!("running-cache telemetry failed: {error}");
                        // Preserve the last good live report. Reporting failure
                        // never touches the native writer or calculation.
                        let _ = Self::write_atomic(
                            &work.join("running-cache-telemetry-error.json"),
                            &json!({"utc_unix_seconds":utc,"pid":std::process::id(),
                                "session":session,"phase":phase,"error":error,
                                "last_good_utc_unix_seconds":last_good_utc}),
                        );
                    }
                    match receiver.recv_timeout(Duration::from_secs(10)) {
                        Ok(()) | Err(mpsc::RecvTimeoutError::Disconnected) => break,
                        Err(mpsc::RecvTimeoutError::Timeout) => {}
                    }
                }
                json!({"session":session,"phase":phase,"reports":reports,"errors":errors,
                    "last_good_utc_unix_seconds":last_good_utc,"last_error":last_error,
                    "stopped_and_joined":true})
            })?;
        Ok(Self {
            stop,
            worker: Some(worker),
        })
    }

    fn write_atomic(path: &Path, document: &serde_json::Value) -> std::io::Result<()> {
        let pending = path.with_extension(format!("{}.pending", std::process::id()));
        let mut file = File::create(&pending)?;
        serde_json::to_writer_pretty(&mut file, document).map_err(std::io::Error::other)?;
        file.sync_all()?;
        fs::rename(pending, path)?;
        File::open(path.parent().unwrap_or(Path::new(".")))?.sync_all()
    }

    fn stop(mut self) -> serde_json::Value {
        let _ = self.stop.send(());
        match self.worker.take().map(JoinHandle::join) {
            Some(Ok(report)) => report,
            Some(Err(_)) => json!({"stopped_and_joined":true,"thread_panicked":true}),
            None => serde_json::Value::Null,
        }
    }
}

impl Drop for RunningTelemetry {
    fn drop(&mut self) {
        let _ = self.stop.send(());
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
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
    let mut running_cache = None;
    let mut running_cache_capture_only = false;
    let mut running_cache_no_read_delta = false;
    let mut prepared_input = None;
    let mut prepare_only = false;
    let mut fresh_prepare = false;
    let mut dress_only = false;
    let mut warm = false;
    let mut monitor = RustRedMonitorOptions::default();
    let mut coefficient_representation = CandidateCacheRepresentation::Sparse;
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
            Some("--cache") => {
                cache = Some(PathBuf::from(arguments.next().ok_or("missing cache path")?))
            }
            Some("--running-cache") => {
                running_cache = Some(PathBuf::from(
                    arguments.next().ok_or("missing running-cache path")?,
                ))
            }
            Some("--running-cache-capture-only") => running_cache_capture_only = true,
            Some("--running-cache-no-read-delta") => running_cache_no_read_delta = true,
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

            Some("--coefficient-representation") => {
                coefficient_representation = match arguments
                    .next()
                    .and_then(|value| value.into_string().ok())
                    .as_deref()
                {
                    Some("sparse") => CandidateCacheRepresentation::Sparse,
                    Some("factorized") => CandidateCacheRepresentation::Factorized,
                    _ => {
                        return Err(
                            "coefficient representation must be sparse or factorized".into()
                        );
                    }
                }
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
    if running_cache_capture_only && running_cache.is_none() {
        return Err("--running-cache-capture-only requires --running-cache".into());
    }
    if running_cache_capture_only && cache.is_some() {
        return Err("capture-only running cache cannot be combined with --cache".into());
    }
    if running_cache_no_read_delta && running_cache.is_none() {
        return Err("--running-cache-no-read-delta requires --running-cache".into());
    }
    let representation_name = match coefficient_representation {
        CandidateCacheRepresentation::Sparse => "sparse",
        CandidateCacheRepresentation::Factorized => "factorized",
    };
    let work = work.ok_or("--work-dir is required")?;
    if monitor.display == RustRedMonitorDisplay::Tui && monitor.progress_snapshot.is_none() {
        monitor.progress_snapshot = Some(work.join("monitor.json"));
    }
    monitor.validate()?;

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
    if prepared_input.is_none() && fixture::LadderFixture::has_prepared_input(&work)? {
        if fresh_prepare {
            return Err("--fresh-prepare requires a new --work-dir; preserve the completed prepared-input bundle".into());
        }
        prepared_input = Some(work.clone());
    }
    let fresh_native_preparation = prepared_input.is_none();
    let input_started = Instant::now();
    let prepared = match &prepared_input {
        Some(directory) => fixture::LadderFixture::load_prepared(directory, projection)?,
        None => fixture::prepare(&work, projection)?,
    };
    let input_loading_seconds =
        (!fresh_native_preparation).then(|| input_started.elapsed().as_secs_f64());
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
        "coefficient_representation":representation_name,
        "source_commit":git(&["rev-parse", "HEAD"]), "source_dirty":git(&["status", "--porcelain"]).map(|value| !value.is_empty()),
        "executable_sha256":format!("{:x}", Sha256::digest(std::fs::read(std::env::current_exe()?)?)),
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
        "persistent_cache":cache, "running_cache_path":running_cache,
        "running_cache_capture_only":running_cache_capture_only,
        "running_cache_read_delta":!running_cache_no_read_delta,
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
        });
    }
    write_report(&report)?;
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
    settings.rustred_application_options.cache_representation = coefficient_representation;
    // The public loader prepares its default Sparse session. Selecting
    // Factorized reloads that session on the first evaluation; retain and
    // report this setup cost instead of introducing another loading path.
    report["representation_setup"] = json!({
        "initial_data_load":"sparse",
        "evaluation_requested":representation_name,
        "first_evaluation_reloads_application":coefficient_representation != CandidateCacheRepresentation::Sparse,
        "reload_timing_scope":"Any representation reload is included in reduction_seconds, not data_loading_seconds",
        "owner_preparation_counter_scope":"Before/after counters belong to different sessions when the representation reloads; their subtraction is not an isolated same-session preparation measurement"
    });
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
        .map(|loops| Ok((loops, engine.cache_statistics(loops)?.has_read_backing())))
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
    let running_read_backing = running_cache.is_some() && !running_cache_capture_only;
    let running_started = Instant::now();
    let running_control = running_cache
        .as_ref()
        .map(|path| {
            engine.attach_running_cache(
                5,
                path,
                &RustRedRunningCacheOptions {
                    read_backing: running_read_backing,
                    read_delta: !running_cache_no_read_delta,
                    ..Default::default()
                },
            )
        })
        .transpose()?;
    let evaluation = (|| -> fixture::Result<()> {
        report["running_cache_loading_seconds"] = json!(
            running_control
                .as_ref()
                .map(|_| running_started.elapsed().as_secs_f64())
        );
        report["running_cache_read_backing"] = json!(running_read_backing);
        report["running_cache_attached"] =
            json!(running_control.is_some() && engine.cache_statistics(5)?.attached);
        let caches_during_evaluation = (1..=5)
            .map(|loops| Ok((loops, engine.cache_statistics(loops)?.has_read_backing())))
            .collect::<fixture::Result<BTreeMap<_, _>>>()?;
        if caches_during_evaluation.iter().any(|(&loops, &attached)| {
            attached != ((cache.is_some() || running_read_backing) && loops == 5)
        }) {
            return Err(
                "measurement persistent backing differs from the declared five-loop cache policy"
                    .into(),
            );
        }
        report["persistent_cache_attached"] = json!(caches_during_evaluation[&5]);
        report["persistent_cache_attachments_before"] = json!(caches_before);
        report["persistent_cache_attachments_during_evaluation"] = json!(caches_during_evaluation);
        report["persistent_cache_policy"] =
            json!(if running_read_backing && running_cache_no_read_delta {
                "immutable-five-loop-base-only-with-running-capture"
            } else if running_read_backing && cache.is_some() {
                "manual-five-loop-base-and-running-delta"
            } else if running_read_backing {
                "running-five-loop-delta"
            } else if cache.is_some() {
                "manual-five-loop-only"
            } else {
                "none"
            });
        let flush_running_cache =
            |phase: &str| -> (serde_json::Value, Option<Box<dyn std::error::Error>>) {
                let Some(control) = &running_control else {
                    return (serde_json::Value::Null, None);
                };
                let started = Instant::now();
                let flushed = control.flush();
                let seconds = started.elapsed().as_secs_f64();
                let flush_complete = flushed.as_ref().is_ok_and(|report| {
                    report.pending_records == 0 && report.writer_error.is_none()
                });
                let document = match &flushed {
                    Ok(report) => Some(report.to_document()),
                    Err(_) => control.report().ok().map(|report| report.to_document()),
                };
                let flush_error = flushed.as_ref().err().map(ToString::to_string).or_else(|| {
                    (!flush_complete)
                        .then(|| "native flush did not commit all accepted records".to_owned())
                });
                let mut persistence_error: Option<Box<dyn std::error::Error>> = flushed
                    .err()
                    .map(|error| Box::new(error) as Box<dyn std::error::Error>);
                if persistence_error.is_none() && !flush_complete {
                    persistence_error =
                        Some("native flush did not commit all accepted records".into());
                }
                let mut snapshot_seconds = None;
                let mut snapshot_sidecar = None;
                let mut snapshot_descriptor = None;
                let mut snapshot_error = None;
                let snapshot = if flush_complete {
                    let started = Instant::now();
                    let frozen = (|| -> fixture::Result<serde_json::Value> {
                        let delta = fs::canonicalize(
                            running_cache
                                .as_ref()
                                .ok_or("running cache has no delta path")?,
                        )?;
                        let directory = fs::canonicalize(&work)?.join("running-cache-snapshots");
                        if directory.starts_with(&delta) {
                            return Err("snapshot must be outside the running delta".into());
                        }
                        fs::create_dir_all(&directory)?;
                        let directory = fs::canonicalize(directory)?;
                        if directory.starts_with(&delta) {
                            return Err(
                                "snapshot directory resolves inside the running delta".into()
                            );
                        }
                        let stamp = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
                        let sidecar = directory.join(format!(
                            "ladder-{phase}-{}-{stamp}.rrindex",
                            std::process::id()
                        ));
                        let descriptor = sidecar.with_extension("json");
                        if sidecar.exists() || descriptor.exists() {
                            return Err("preserve the existing running-cache snapshot".into());
                        }
                        snapshot_sidecar = Some(sidecar.clone());
                        snapshot_descriptor = Some(descriptor.clone());
                        let document = control.snapshot(&sidecar)?.to_document();
                        let mut file = OpenOptions::new()
                            .write(true)
                            .create_new(true)
                            .open(&descriptor)?;
                        serde_json::to_writer_pretty(&mut file, &document)?;
                        file.sync_all()?;
                        File::open(directory)?.sync_all()?;
                        Ok(document)
                    })();
                    snapshot_seconds = Some(started.elapsed().as_secs_f64());
                    match frozen {
                        Ok(document) => Some(document),
                        Err(error) => {
                            snapshot_error = Some(error.to_string());
                            persistence_error = Some(error);
                            None
                        }
                    }
                } else {
                    None
                };
                let mut receipt = json!({
                    "path":running_cache,"read_backing":running_read_backing,
                    "phase":phase,"read_delta":!running_cache_no_read_delta,
                    "flush_seconds":seconds,"report":document,
                    "flush_error":flush_error,"flush_complete":flush_complete,
                    "snapshot":snapshot,"snapshot_error":snapshot_error,"snapshot_seconds":snapshot_seconds,
                    "snapshot_sidecar":snapshot_sidecar,"snapshot_descriptor":snapshot_descriptor,
                });
                let saved = (|| -> fixture::Result<()> {
                    let stamp = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
                    let temporary = work.join(format!(
                        "running-cache-{phase}-{}-{stamp}.json.pending",
                        std::process::id()
                    ));
                    let mut file = OpenOptions::new()
                        .write(true)
                        .create_new(true)
                        .open(&temporary)?;
                    serde_json::to_writer_pretty(&mut file, &receipt)?;
                    file.sync_all()?;
                    fs::rename(&temporary, work.join("running-cache.json"))?;
                    File::open(&work)?.sync_all()?;
                    Ok(())
                })();
                if let Err(error) = saved {
                    eprintln!("could not write running-cache receipt: {error}");
                    receipt["receipt_error"] = json!(error.to_string());
                    if persistence_error.is_none() {
                        persistence_error = Some(error);
                    }
                }
                (receipt, persistence_error)
            };
        let start_telemetry = |phase| {
            let Some(control) = &running_control else {
                return (None, None);
            };
            match RunningTelemetry::start(control.clone(), work.clone(), phase) {
                Ok(telemetry) => (Some(telemetry), None),
                Err(error) => {
                    eprintln!("could not start running-cache telemetry: {error}");
                    let _ = RunningTelemetry::write_atomic(
                        &work.join("running-cache-telemetry-error.json"),
                        &json!({"pid":std::process::id(),"phase":phase,
                        "utc_unix_seconds":SystemTime::now().duration_since(UNIX_EPOCH).ok().map(|time|time.as_secs_f64()),
                        "error":error.to_string(),"reporter_started":false}),
                    );
                    (None, Some(error.to_string()))
                }
            }
        };
        let owners_before = engine.rustred_application_owner_preparation_statistics()?;
        LOGGER.0.lock().map_err(|error| error.to_string())?.clear();
        let (telemetry, telemetry_error) = start_telemetry("cold");
        report["running_cache_telemetry_start_error"] = json!(telemetry_error);
        report["execution_threads"] = json!({"math":1,
        "native_io_writer":usize::from(running_control.is_some()),
        "telemetry":usize::from(telemetry.is_some()),
        "scope":"dedicated collector math, native writer and telemetry; not measured OS thread count"});
        let started = Instant::now();
        let outcome = if monitor.enabled() {
            engine.evaluate_rustred_symbolic_with_monitor(
                &settings,
                prepared.input.as_view(),
                &monitor,
            )
        } else {
            engine.evaluate_rustred_symbolic(&settings, prepared.input.as_view())
        };
        report["reduction_seconds"] = json!(started.elapsed().as_secs_f64());
        report["running_cache_telemetry"] = telemetry
            .map(RunningTelemetry::stop)
            .unwrap_or(serde_json::Value::Null);
        // The control does not borrow the reducer. Flush completed native nodes
        // before fallible telemetry or error propagation can leave this scope.
        let (running_receipt, running_flush_error) = flush_running_cache("cold");
        report["running_cache"] = running_receipt;
        report["end_to_end_seconds"] = json!(total_started.elapsed().as_secs_f64());
        let result = match outcome {
            Ok(result) => result,
            Err(error) => {
                report["reduction_error"] = json!(error.to_string());
                report["warnings"] = json!(WARNINGS.lock().ok().map(|warnings| warnings.clone()));
                report["stages"] = json!(LOGGER.0.lock().ok().map(|stages| stages.clone()));
                if let Err(receipt_error) = write_report(&report) {
                    eprintln!("could not write failed reduction receipt: {receipt_error}");
                }
                return Err(error.into());
            }
        };
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
        if let Some(error) = running_flush_error {
            report["running_cache_capture_failed"] = json!(true);
            if let Err(receipt_error) = write_report(&report) {
                eprintln!("could not write failed running-cache receipt: {receipt_error}");
            }
            return Err(error);
        }
        if warm {
            LOGGER.0.lock().map_err(|error| error.to_string())?.clear();
            let (telemetry, telemetry_error) = start_telemetry("warm");
            report["warm_running_cache_telemetry_start_error"] = json!(telemetry_error);
            let started = Instant::now();
            let outcome = if monitor.enabled() {
                engine.evaluate_rustred_symbolic_with_monitor(
                    &settings,
                    prepared.input.as_view(),
                    &monitor,
                )
            } else {
                engine.evaluate_rustred_symbolic(&settings, prepared.input.as_view())
            };
            report["warm_reduction_seconds"] = json!(started.elapsed().as_secs_f64());
            report["warm_running_cache_telemetry"] = telemetry
                .map(RunningTelemetry::stop)
                .unwrap_or(serde_json::Value::Null);
            let (running_receipt, running_flush_error) = flush_running_cache("warm");
            report["warm_running_cache"] = running_receipt;
            let second = match outcome {
                Ok(result) => result,
                Err(error) => {
                    report["warm_reduction_error"] = json!(error.to_string());
                    if let Err(receipt_error) = write_report(&report) {
                        eprintln!("could not write failed warm reduction receipt: {receipt_error}");
                    }
                    return Err(error.into());
                }
            };
            if let Some(error) = running_flush_error {
                report["running_cache_capture_failed"] = json!(true);
                if let Err(receipt_error) = write_report(&report) {
                    eprintln!("could not write failed running-cache receipt: {receipt_error}");
                }
                return Err(error);
            }
            report["warm_stages"] = json!(*LOGGER.0.lock().map_err(|error| error.to_string())?);
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
    // Close exactly once after the final cold/warm durability barrier. Even an
    // ordinary calculation/reporting error must release the physical writer.
    let close_error = running_control.as_ref().and_then(|control| {
        let started = Instant::now();
        let closed = control.close();
        let document = match &closed {
            Ok(closed) => Some(closed.to_document()),
            Err(_) => control.report().ok().map(|report| report.to_document()),
        };
        report["running_cache_close"] = json!({"seconds":started.elapsed().as_secs_f64(),
            "report":document,"error":closed.as_ref().err().map(ToString::to_string)});
        closed.err()
    });
    let terminal_receipt_error = running_control.as_ref().and_then(|_| {
        let mut terminal = if warm && !report["warm_running_cache"].is_null() {
            report["warm_running_cache"].clone()
        } else {
            report["running_cache"].clone()
        };
        terminal["flush_report"] = terminal["report"].clone();
        terminal["report"] = report["running_cache_close"]["report"].clone();
        terminal["close"] = report["running_cache_close"].clone();
        terminal["terminal"] = json!(true);
        terminal["utc_unix_seconds"] = json!(
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .ok()
                .map(|time| time.as_secs_f64())
        );
        RunningTelemetry::write_atomic(&work.join("running-cache.json"), &terminal).err()
    });
    if let Some(error) = &terminal_receipt_error {
        report["running_cache_terminal_receipt_error"] = json!(error.to_string());
        eprintln!("could not write terminal running-cache receipt: {error}");
    }
    if let Err(error) = evaluation {
        if let Err(receipt_error) = write_report(&report) {
            eprintln!("could not write terminal reduction receipt: {receipt_error}");
        }
        return Err(error);
    }
    if let Some(error) = close_error {
        report["running_cache_capture_failed"] = json!(true);
        if let Err(receipt_error) = write_report(&report) {
            eprintln!("could not write terminal running-cache receipt: {receipt_error}");
        }
        return Err(error.into());
    }
    if let Some(error) = terminal_receipt_error {
        return Err(error.into());
    }
    write_report(&report)?;
    println!("{report}");
    Ok(())
}
