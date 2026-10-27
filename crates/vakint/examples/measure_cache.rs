//! Fresh-process public-evaluator measurements of a frozen cache suite.
#[path = "build_cache/cache_file_args.rs"]
mod cache_file_args;
#[path = "build_cache/cache_hit_capture.rs"]
mod cache_hit_capture;
#[path = "build_cache/cache_identity.rs"]
mod cache_identity;
#[path = "build_cache/owner_proof.rs"]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "Shared mapper creates frontier steps; measurement only validates saved proofs"
    )
)]
mod owner_proof;
#[path = "build_cache/resource_guard_admission.rs"]
mod resource_guard_admission;
#[path = "build_cache/workload.rs"]
mod workload;

use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File, OpenOptions},
    io::Write,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::Mutex,
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};

use cache_file_args::CacheFileArgs;
use cache_hit_capture::{CacheHitDiagnostic, ClassificationLogging};
use cache_identity::{CacheIdentity, digest};
use eyre::{Context, Result, bail, ensure};
use owner_proof::OwnerProof;
use owner_proof::validate_hit_diagnostic;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use symbolica::atom::{AtomCore, AtomView};
use vakint::{
    EvaluationMethod, EvaluationOrder, OfficialExternalDataOptions, RustRedEvaluationOptions,
    RustRedMasterMode, RustRedMonitorDisplay, RustRedMonitorOptions, RustRedRunningCacheOptions,
    RustRedSnapshotViewer, Vakint, VakintSettings,
};
use workload::{Root, Workload};

#[derive(Default)]
struct Stages {
    seconds: BTreeMap<String, Vec<f64>>,
    lowered: Vec<String>,
    warnings: Vec<String>,
    cache_hits: CacheHitDiagnostic,
}

impl Stages {
    fn authenticate(
        &self,
        lowered: &[Vec<i64>],
        protected: &[Vec<i64>],
        provisional: bool,
        proof: Option<&OwnerProof>,
        case_id: &str,
    ) -> Result<bool> {
        ensure!(
            self.cache_hits.input_roots_error.is_none(),
            "invalid input-root telemetry: {:?}",
            self.cache_hits.input_roots_error
        );
        if !self.cache_hits.input_roots_seen || lowered.is_empty() || protected.is_empty() {
            ensure!(
                provisional,
                "native key authentication unavailable: require complete telemetry and frozen protected_roots, or use --provisional without final claims"
            );
            return Ok(false);
        }
        let lowered = lowered.iter().cloned().collect::<BTreeSet<_>>();
        ensure!(
            self.cache_hits.input_roots == lowered,
            "public input lowered input roots differ from frozen lowered_roots"
        );
        if let Some(proof) = proof {
            proof.observe(case_id, &self.cache_hits.input_roots, protected)?;
            Ok(true)
        } else {
            ensure!(
                provisional,
                "native owner authentication requires a validated pinned proof"
            );
            Ok(false)
        }
    }
}

struct StageLogger(Mutex<Stages>);
static LOGGER: StageLogger = StageLogger(Mutex::new(Stages {
    seconds: BTreeMap::new(),
    lowered: Vec::new(),
    warnings: Vec::new(),
    cache_hits: CacheHitDiagnostic::new(),
}));

impl StageLogger {
    fn stages(&self) -> std::sync::MutexGuard<'_, Stages> {
        self.0.lock().unwrap_or_else(|poison| {
            let mut stages = poison.into_inner();
            stages.cache_hits.hit_error = Some("cache-hit logger mutex was poisoned".into());
            stages
        })
    }
}

impl log::Log for StageLogger {
    fn enabled(&self, metadata: &log::Metadata<'_>) -> bool {
        if metadata.target() == "vakint::benchmark_cache_hits" {
            return cfg!(feature = "experimental-rustred") && ClassificationLogging::enabled();
        }
        metadata.level() <= log::Level::Warn
            || matches!(
                metadata.target(),
                "vakint::saved_application" | "vakint::cache_workload"
            )
    }

    fn log(&self, record: &log::Record<'_>) {
        if !self.enabled(record.metadata()) {
            return;
        }
        let message = record.args().to_string();
        {
            let mut stages = self.stages();
            if record.target() == "vakint::benchmark_cache_hits" {
                if stages.cache_hits.hit_error.is_none()
                    && let Err(error) = stages.cache_hits.capture_hit_event(&message)
                {
                    stages.cache_hits.hit_error = Some(error.to_string());
                }
                return;
            }
            if record.target() == "vakint::cache_workload" {
                if let Err(error) = stages.cache_hits.capture_input_roots(&message) {
                    stages.cache_hits.input_roots_error = Some(error.to_string());
                }
                return;
            }
            if record.level() <= log::Level::Warn {
                eprintln!("{message}");
                stages.warnings.push(message.clone());
            }
            if message.starts_with("phase=scalar_lowered ") {
                stages.lowered.push(message.clone());
            }
            let phase = message
                .split_whitespace()
                .find_map(|field| field.strip_prefix("phase="));
            let elapsed = message.split_whitespace().find_map(|field| {
                field
                    .strip_prefix("elapsed_seconds=")
                    .and_then(|value| value.parse::<f64>().ok())
            });
            if let (Some(phase), Some(elapsed)) = (phase, elapsed) {
                stages
                    .seconds
                    .entry(phase.into())
                    .or_default()
                    .push(elapsed);
            }
        }
    }

    fn flush(&self) {}
}

struct Options {
    artifacts: PathBuf,
    workload: PathBuf,
    source_census: PathBuf,
    cache: Option<PathBuf>,
    cache_file_options: CacheFileArgs,
    running_cache: Option<PathBuf>,
    running_cache_capture_only: bool,
    running_cache_no_read_delta: bool,
    work_dir: PathBuf,
    suite: String,
    repetitions: usize,
    max_case_seconds: Option<u64>,
    case_id: Option<String>,
    variant: String,
    worker: bool,
    sample: Option<usize>,
    worker_deadline: Option<f64>,
    provisional: bool,
    only_case: Option<String>,
    monitor: RustRedMonitorOptions,
    classify_cache_hits: bool,
}

impl Options {
    fn case_deadline(&self, now: f64) -> Result<Option<f64>> {
        ensure!(now.is_finite(), "invalid supervisor clock");
        self.max_case_seconds
            .map(|seconds| {
                ensure!(seconds > 0, "positive diagnostic deadline required");
                Ok(now + seconds as f64)
            })
            .transpose()
    }

    fn worker_allowance(&self, now: f64) -> Result<Option<Duration>> {
        match (self.max_case_seconds, self.worker_deadline) {
            (None, None) => Ok(None),
            (Some(seconds), Some(deadline)) => {
                let remaining = (deadline - now).min(seconds as f64);
                ensure!(
                    seconds > 0 && deadline.is_finite() && now.is_finite() && remaining > 0.0,
                    "worker diagnostic deadline has expired or is invalid"
                );
                Ok(Some(Duration::try_from_secs_f64(remaining)?))
            }
            _ => bail!("worker deadline must match the explicit diagnostic budget"),
        }
    }

    fn deadline_exhausted(
        &self,
        deadline: Option<f64>,
        elapsed: Duration,
        now: f64,
    ) -> Result<bool> {
        match (self.max_case_seconds, deadline) {
            (None, None) => Ok(false),
            (Some(seconds), Some(deadline)) => {
                ensure!(
                    seconds > 0 && deadline.is_finite() && now.is_finite(),
                    "invalid diagnostic deadline"
                );
                Ok(elapsed >= Duration::from_secs(seconds) || now >= deadline)
            }
            _ => bail!("supervisor deadline must match the explicit diagnostic budget"),
        }
    }

    fn interrupted_supervision(&self, active: &Value, now: f64) -> Result<Value> {
        if self.max_case_seconds.is_none() {
            bail!(
                "an unbounded worker may still be active; preserve the interrupted attempt \
                 until its terminal completion is independently verified"
            );
        }
        let deadline = active["deadline_unix_seconds"]
            .as_f64()
            .ok_or_else(|| eyre::eyre!("invalid active diagnostic deadline"))?;
        ensure!(
            deadline.is_finite() && now.is_finite() && now >= deadline,
            "an interrupted supervisor's worker may still be active; preserve its diagnostic deadline"
        );
        Ok(
            json!({"success":false,"interrupted":true,"deadline_exhausted":true,
            "error":"interrupted attempt retained without retry"}),
        )
    }

    fn arm_worker_watchdog(&self) -> Result<Value> {
        let active: Value = serde_json::from_slice(&fs::read(self.work_dir.join("running.json"))?)?;
        ensure!(
            active.get("deadline_unix_seconds") == Some(&json!(self.worker_deadline))
                && active["case_id"].as_str() == self.case_id.as_deref()
                && active["variant"] == self.variant
                && active["sample"].as_u64() == self.sample.map(|sample| sample as u64),
            "worker has no matching durable active case"
        );
        let now = SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs_f64();
        if let Some(allowance) = self.worker_allowance(now)? {
            let deadline = self
                .worker_deadline
                .ok_or_else(|| eyre::eyre!("bounded worker requires a supervisor deadline"))?;
            std::thread::Builder::new()
                .name("measurement-deadline".into())
                .spawn(move || {
                    let started = Instant::now();
                    loop {
                        let remaining = allowance.saturating_sub(started.elapsed());
                        let expired = SystemTime::now()
                            .duration_since(UNIX_EPOCH)
                            .map_or(true, |now| now.as_secs_f64() >= deadline);
                        if remaining.is_zero() || expired {
                            std::process::exit(124);
                        }
                        std::thread::sleep(remaining.min(Duration::from_secs(1)));
                    }
                })
                .context("start independent diagnostic measurement watchdog")?;
        }
        Ok(active["identity"].clone())
    }

    fn parse() -> Result<Option<Self>> {
        Self::parse_args(std::env::args().skip(1))
    }

    fn parse_args(mut args: impl Iterator<Item = String>) -> Result<Option<Self>> {
        let mut artifacts = None;
        let mut workload = None;
        let mut source_census = None;
        let mut cache = None;
        let mut cache_file_options = CacheFileArgs::default();
        let mut running_cache = None;
        let mut running_cache_capture_only = false;
        let mut running_cache_no_read_delta = false;
        let mut work_dir = None;
        let mut suite = "tuning".to_owned();
        let mut repetitions = 3;
        let mut max_case_seconds: Option<u64> = None;
        let mut case_id = None;
        let mut variant = "uncached".to_owned();
        let mut worker = false;
        let mut sample = None;
        let mut worker_deadline: Option<f64> = None;
        let mut provisional = false;
        let mut only_case = None;
        let mut monitor = RustRedMonitorOptions::default();
        let mut classify_cache_hits = false;
        while let Some(arg) = args.next() {
            if matches!(arg.as_str(), "--help" | "-h") {
                println!(
                    "measure_cache --artifacts PATH --workload JSON --source-census FILE --work-dir PATH\n\
                     [--cache FILE.gz|NATIVE_DIRECTORY] [--suite tuning|held-out] [--repetitions N] [--max-case-seconds N]\n\
                     [--cache-max-compressed-bytes N] [--cache-max-decoded-bytes N] [--cache-max-index-bytes N]\n\
                     [--cache-decoded-directory PATH] uses file-backed decoded storage in an existing directory.\n\
                     File grants require --cache; omitted grants keep compact defaults and are never raised automatically.\n\
                     [--running-cache PATH] [--running-cache-capture-only] persist completed native records.\n\
                     PATH is a capture root; each case/sample/variant owns a separate native delta.\n\
                     Capture-only disables persistent reads and cannot be combined with --cache.\n\
                     Held-out measurements forbid --running-cache and read --cache without modifying it.\n\
                     [--running-cache-no-read-delta] reads immutable bases only; requires --running-cache.\n\
                     [--provisional] permits legacy functional pilots without native key authentication.\n\
                     [--only-case ID] selects one frozen case and requires --provisional.\n\
                     Every case/sample/variant runs in a fresh process. Default: three samples, no deadline.\n\
                     --max-case-seconds N opts into a positive diagnostic budget for the whole worker.\n\
                     Omit --cache for an uncached baseline; final cache comparisons require three samples.\n\
                     [--monitor] [--progress-snapshot PATH] [--progress-interval-seconds N]\n\
                     Monitoring defaults off; --monitor needs progress-tui and uses worker work-dir/monitor.json.\n\
                     A snapshot path enables native JSON without a TUI; interval defaults to 1s.\n\
                     Monitored timings include diagnostic observation overhead.\n\
                     [--classify-cache-hits] records a lossless diagnostic for offline owner mapping.\n\
                     Requires experimental-rustred, --cache and --repetitions 1; runs only the cached variant.\n\
                     Forbids presentation, running capture and performance claims.\n\
                     Run with fixed affinity in a serialized measurement window."
                );
                return Ok(None);
            }
            if arg == "--worker" {
                worker = true;
                continue;
            }
            if arg == "--provisional" {
                provisional = true;
                continue;
            }
            if arg == "--running-cache-capture-only" {
                running_cache_capture_only = true;
                continue;
            }
            if arg == "--running-cache-no-read-delta" {
                running_cache_no_read_delta = true;
                continue;
            }
            if arg == "--monitor" {
                monitor.display = RustRedMonitorDisplay::Tui;
                continue;
            }
            if arg == "--classify-cache-hits" {
                classify_cache_hits = true;
                continue;
            }
            let value = args
                .next()
                .ok_or_else(|| eyre::eyre!("missing value for {arg}"))?;
            match arg.as_str() {
                "--artifacts" => artifacts = Some(value.into()),
                "--workload" => workload = Some(value.into()),
                "--source-census" => source_census = Some(value.into()),
                "--cache" => cache = Some(value.into()),
                flag if CacheFileArgs::recognizes(flag) => {
                    cache_file_options.set(flag, value.into())?;
                }
                "--running-cache" => running_cache = Some(value.into()),
                "--work-dir" => work_dir = Some(value.into()),
                "--suite" => suite = value,
                "--repetitions" => repetitions = value.parse()?,
                "--max-case-seconds" => max_case_seconds = Some(value.parse()?),
                "--case-id" => case_id = Some(value),
                "--only-case" => only_case = Some(value),
                "--variant" => variant = value,
                "--sample" => sample = Some(value.parse()?),
                "--worker-deadline" => worker_deadline = Some(value.parse()?),
                "--progress-snapshot" => monitor.progress_snapshot = Some(value.into()),
                "--progress-interval-seconds" => {
                    monitor.interval = Duration::try_from_secs_f64(value.parse()?)?;
                }
                _ => bail!("unknown argument {arg}"),
            }
        }
        cache_file_options.normalize(cache.as_deref())?;
        ensure!(
            matches!(suite.as_str(), "tuning" | "held-out"),
            "unknown suite"
        );
        ensure!(
            repetitions > 0 && max_case_seconds.is_none_or(|seconds| seconds > 0),
            "positive repetitions and an optional positive diagnostic deadline required"
        );
        ensure!(
            matches!(variant.as_str(), "uncached" | "cached"),
            "unknown variant"
        );
        ensure!(
            variant != "cached" || cache.is_some(),
            "cached worker needs --cache"
        );
        ensure!(
            suite != "held-out" || running_cache.is_none(),
            "held-out measurements require a read-only --cache and forbid --running-cache"
        );
        ensure!(
            !running_cache_capture_only || (running_cache.is_some() && cache.is_none()),
            "--running-cache-capture-only requires --running-cache and forbids --cache"
        );
        ensure!(
            !running_cache_no_read_delta || running_cache.is_some(),
            "--running-cache-no-read-delta requires --running-cache"
        );
        ensure!(
            only_case.is_none() || provisional,
            "--only-case requires --provisional; final comparisons use the complete frozen suite"
        );
        #[cfg(feature = "experimental-rustred")]
        {
            monitor.cache_hit_diagnostic = classify_cache_hits;
        }
        ensure!(
            !classify_cache_hits || cfg!(feature = "experimental-rustred"),
            "--classify-cache-hits requires experimental-rustred"
        );
        ensure!(
            !classify_cache_hits || !monitor.enabled(),
            "lossless cache-hit classification cannot be combined with presentation"
        );
        ensure!(
            !classify_cache_hits
                || (cache.is_some() && repetitions == 1 && running_cache.is_none()),
            "classification requires one cached diagnostic and forbids running capture"
        );
        ensure!(
            !classify_cache_hits || !worker || variant == "cached",
            "a classification worker must use the cached variant"
        );
        ensure!(
            worker_deadline.is_none_or(|deadline| deadline.is_finite()),
            "finite worker diagnostic deadline required"
        );
        ensure!(
            !worker || max_case_seconds.is_some() == worker_deadline.is_some(),
            "worker deadline must match the explicit diagnostic budget"
        );
        monitor.validate()?;
        Ok(Some(Self {
            artifacts: artifacts.ok_or_else(|| eyre::eyre!("--artifacts is required"))?,
            workload: workload.ok_or_else(|| eyre::eyre!("--workload is required"))?,
            source_census: source_census
                .ok_or_else(|| eyre::eyre!("--source-census is required"))?,
            work_dir: work_dir.ok_or_else(|| eyre::eyre!("--work-dir is required"))?,
            cache,
            cache_file_options,
            running_cache,
            running_cache_capture_only,
            running_cache_no_read_delta,
            suite,
            repetitions,
            max_case_seconds,
            case_id,
            variant,
            worker,
            sample,
            worker_deadline,
            provisional,
            only_case,
            monitor,
            classify_cache_hits,
        }))
    }

    fn append_cache_file_arguments(&self, command: &mut Command) {
        let Some(options) = self.cache_file_options.load_options().file_options else {
            return;
        };
        for (flag, bytes) in [
            (
                "--cache-max-compressed-bytes",
                options.limits.max_compressed_bytes,
            ),
            (
                "--cache-max-decoded-bytes",
                options.limits.max_decoded_bytes(),
            ),
            (
                "--cache-max-index-bytes",
                options.limits.chunk.max_index_owned_bytes,
            ),
        ] {
            command.arg(flag).arg(bytes.to_string());
        }
        if let Some(directory) = options.decoded_directory {
            command.arg("--cache-decoded-directory").arg(directory);
        }
    }

    fn identity(&self) -> Result<Value> {
        let workload: Workload = serde_json::from_slice(&fs::read(&self.workload)?)?;
        let mut identity = json!({
            "workload_sha256":digest(&self.workload)?,
            "load_prepared_caches":false,
            "cached_variant_persistent_loops":[5],
            "source_census_sha256":digest(&self.source_census)?,
            "data_manifest_sha256":digest(&self.artifacts.join("manifest.json"))?,
            "cache_sha256":null,
            "cache_file_options":self.cache_file_options.document()?,
            "executable_sha256":digest(&std::env::current_exe()?)?,
            "compiled_driver_sha256":format!("{:x}",Sha256::digest(include_bytes!("measure_cache.rs"))),
            "compiled_workload_schema_sha256":format!("{:x}",Sha256::digest(include_bytes!("build_cache/workload.rs"))),
            "collector_helper_sources_sha256":{
                "crates/vakint/examples/build_cache/cache_identity.rs":format!("{:x}",Sha256::digest(include_bytes!("build_cache/cache_identity.rs"))),
                "crates/vakint/examples/build_cache/cache_hit_diagnostic.rs":format!("{:x}",Sha256::digest(include_bytes!("build_cache/cache_hit_diagnostic.rs"))),
                "crates/vakint/examples/build_cache/cache_hit_capture.rs":format!("{:x}",Sha256::digest(include_bytes!("build_cache/cache_hit_capture.rs"))),
            },
            "compiled_owner_proof_sha256":format!("{:x}",Sha256::digest(include_bytes!("build_cache/owner_proof.rs"))),
            "compiled_cache_file_args_sha256":format!("{:x}",Sha256::digest(include_bytes!("build_cache/cache_file_args.rs"))),
            "owner_proof":OwnerProof::pinned_identity(&workload,&self.workload,&self.artifacts)?,
            "compiled_cargo_lock_sha256":format!("{:x}",Sha256::digest(include_bytes!("../../../Cargo.lock"))),
            "build_rustc":option_env!("VAKINT_CACHE_BUILD_RUSTC"),
            "build_profile":option_env!("VAKINT_CACHE_BUILD_PROFILE"),
            "build_git_revision":option_env!("VAKINT_CACHE_BUILD_GIT_REVISION"),
            "suite":self.suite,"repetitions":self.repetitions,"max_case_seconds":self.max_case_seconds,
            "provisional":self.provisional,"only_case":self.only_case,
            "cache_hit_diagnostic":self.classify_cache_hits,
            "cpu_affinity":fs::read_to_string("/proc/self/status").ok().and_then(|status|status.lines().find_map(|line|line.strip_prefix("Cpus_allowed_list:").map(str::trim).map(str::to_owned))),
            "rayon_num_threads":std::env::var("RAYON_NUM_THREADS").ok(),
        });
        if let Some(path) = &self.cache {
            match CacheIdentity::read(path)? {
                CacheIdentity::File(sha256) => identity["cache_sha256"] = sha256.into(),
                CacheIdentity::Native(indexes) => identity["cache_native_index_identity"] = indexes,
            }
        }
        if self.monitor.enabled() {
            identity["monitor"] = json!({
                "display": match self.monitor.display {
                    RustRedMonitorDisplay::Off => "off",
                    RustRedMonitorDisplay::Tui => "tui",
                    RustRedMonitorDisplay::Json => "json",
                },
                "progress_snapshot": self.monitor.progress_snapshot,
                "default_snapshot": "worker-work-dir/monitor.json",
                "interval_seconds": self.monitor.interval.as_secs_f64(),
                "diagnostic_timing": true,
            });
        }
        if self.classify_cache_hits {
            identity["cache_hit_diagnostic_scope"] = json!({
                "capability":"vakint.initial-routing-frontier-hit-proof.v1",
                "interval_ns":0,"presentation":false,
                "requires_offline_recursive_mapper":true,
                "performance_claim_allowed":false,
            });
        }
        if let Some(path) = &self.running_cache {
            identity["running_cache"] = json!({
                "path":path,"read_backing":!self.running_cache_capture_only,
                "capture_only":self.running_cache_capture_only,
                "read_delta":!self.running_cache_no_read_delta,
                "delta_scope":"one native delta per case, sample and variant",
            });
        }
        Ok(identity)
    }
}

/// Index bytes commit to native bindings, generations, ordered runs and chunk
/// hashes. The native loader validates their framing and any payload it uses;
/// this metadata identity does not reread or certify unused chunk payloads.

fn write_json(path: &Path, value: &Value) -> Result<()> {
    write_bytes(path, &serde_json::to_vec_pretty(value)?)
}

fn write_bytes(path: &Path, bytes: &[u8]) -> Result<()> {
    let temporary = path.with_extension("pending");
    let mut file = File::create(&temporary)?;
    file.write_all(bytes)?;
    file.sync_all()?;
    fs::rename(temporary, path)?;
    File::open(path.parent().unwrap_or(Path::new(".")))?.sync_all()?;
    Ok(())
}

fn filter_roots<'a>(roots: &'a [Root], only_case: Option<&str>) -> Result<Vec<&'a Root>> {
    let selected = roots
        .iter()
        .filter(|root| only_case.is_none_or(|id| root.id == id))
        .collect::<Vec<_>>();
    ensure!(
        !selected.is_empty(),
        "--only-case is absent from the frozen selected suite"
    );
    Ok(selected)
}

fn validate_receipt(
    receipt: &Value,
    identity: &Value,
    case_id: &str,
    variant: &str,
    sample: usize,
) -> Result<()> {
    ensure!(
        &receipt["identity"] == identity
            && receipt["case_id"] == case_id
            && receipt["variant"] == variant
            && receipt["sample"].as_u64() == Some(sample as u64),
        "worker receipt belongs to another identity, case, variant or sample"
    );
    Ok(())
}

fn complete_samples(receipts: &[&Value], variants: &[&str], repetitions: usize) -> bool {
    !receipts.is_empty()
        && receipts.len() == variants.len().saturating_mul(repetitions)
        && variants.iter().all(|variant| {
            let selected = receipts
                .iter()
                .filter(|receipt| receipt["variant"] == *variant)
                .collect::<Vec<_>>();
            selected.len() == repetitions
                && selected.iter().all(|receipt| {
                    receipt["sample"]
                        .as_u64()
                        .is_some_and(|sample| sample < repetitions as u64)
                })
                && selected
                    .iter()
                    .filter_map(|receipt| receipt["sample"].as_u64())
                    .collect::<BTreeSet<_>>()
                    .len()
                    == repetitions
        })
}

fn validate_catalogue_closure(receipt: &Value, directory: &Path) -> Result<()> {
    if receipt["success"] != true {
        return Ok(());
    }
    let closure = &receipt["catalogue_closure"];
    let unknown: BTreeSet<String> = serde_json::from_value(closure["unknown_masters"].clone())?;
    let recorded: BTreeSet<String> =
        serde_json::from_slice(&fs::read(directory.join("unknown-masters.json"))?)?;
    ensure!(
        recorded == unknown
            && closure["unknown_count"].as_u64() == Some(unknown.len() as u64)
            && closure["master_count"]
                .as_u64()
                .is_some_and(|count| count >= unknown.len() as u64)
            && closure["closed"].as_bool() == Some(unknown.is_empty())
            && closure["unknown_sha256"] == digest(&directory.join("unknown-masters.json"))?,
        "saved catalogue classification differs from its full unknown-master inventory"
    );
    Ok(())
}

fn persistent_cache_attachments(engine: &Vakint) -> Result<BTreeMap<usize, bool>> {
    (1..=5)
        .map(|loops| Ok((loops, engine.cache_statistics(loops)?.has_read_backing())))
        .collect()
}

fn measure(options: &Options, root: &Root, identity: &Value, workload: &Workload) -> Result<()> {
    log::set_logger(&LOGGER).map_err(|error| eyre::eyre!(error.to_string()))?;
    log::set_max_level(log::LevelFilter::Debug);
    let input = vakint::vakint_parse!(
        root.public_input
            .as_deref()
            .ok_or_else(|| eyre::eyre!("missing public input"))?
    )
    .map_err(|error| eyre::eyre!(error))?;
    let settings = VakintSettings {
        form_exe_path: "/cache-measurement-must-not-use-form".into(),
        use_dot_product_notation: true,
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        ..Default::default()
    };
    let started = Instant::now();
    let mut engine = Vakint::new()?;
    let data_report = engine.load_official_external_data(
        5,
        &OfficialExternalDataOptions {
            local_artifacts_path: Some(options.artifacts.clone()),
            load_prepared_caches: false,
            ..Default::default()
        },
    )?;
    let caches_before = persistent_cache_attachments(&engine)?;
    ensure!(
        caches_before.values().all(|attached| !attached),
        "baseline must never import a persistent cache"
    );
    let data_loading_seconds = started.elapsed().as_secs_f64();
    let cache_started = Instant::now();
    if options.variant == "cached" {
        engine.load_cache(
            5,
            options.cache.as_ref().unwrap(),
            &options.cache_file_options.load_options(),
        )?;
        engine.clear_cache(5)?;
    }
    let cache_loading_seconds = cache_started.elapsed().as_secs_f64();
    let running_cache_path = options
        .running_cache
        .as_ref()
        .map(|root_path| {
            let sample = options
                .sample
                .ok_or_else(|| eyre::eyre!("running-cache worker requires a sample number"))?;
            Ok::<_, eyre::Report>(root_path.join(format!(
                "case-{:x}-sample-{sample}-{}",
                Sha256::digest(root.id.as_bytes()),
                options.variant
            )))
        })
        .transpose()?;
    let running_read_backing = running_cache_path.is_some() && !options.running_cache_capture_only;
    let running_cache_started = Instant::now();
    let running_cache = running_cache_path
        .as_ref()
        .map(|path| {
            ensure!(
                !path.exists(),
                "preserve the existing worker running cache; choose a fresh --running-cache root"
            );
            engine
                .attach_running_cache(
                    5,
                    path,
                    &RustRedRunningCacheOptions {
                        read_backing: running_read_backing,
                        read_delta: !options.running_cache_no_read_delta,
                        ..Default::default()
                    },
                )
                .map_err(eyre::Report::from)
        })
        .transpose()?;
    let running_cache_loading_seconds = running_cache_started.elapsed().as_secs_f64();
    let mut running_cache_receipt = None;
    let measured = (|| -> Result<()> {
        let running_cache_attached =
            running_cache.is_some() && engine.cache_statistics(5)?.attached;
        let caches_during_evaluation = persistent_cache_attachments(&engine)?;
        ensure!(
            caches_during_evaluation.iter().all(|(&loops, &attached)| {
                attached == ((options.variant == "cached" || running_read_backing) && loops == 5)
            }),
            "measurement persistent backing differs from the declared five-loop cache policy"
        );
        let proof_started = Instant::now();
        let binding = engine
            .cache_statistics(5)?
            .saved_application
            .ok_or_else(|| eyre::eyre!("missing loaded application binding"))?
            .to_document()["binding"]
            .clone();
        let proof = OwnerProof::load(workload, &options.workload, &options.artifacts, &binding)?;
        ensure!(
            options.provisional || proof.is_some(),
            "final measurement requires validated owner proof"
        );
        if let Some(proof) = &proof {
            proof.observe_input(
                &root.id,
                &format!(
                    "{:x}",
                    Sha256::digest(input.to_canonical_string().as_bytes())
                ),
            )?;
        }
        let owner_proof_seconds = proof_started.elapsed().as_secs_f64();
        let owners_before = engine.rustred_application_owner_preparation_statistics()?;
        *LOGGER.stages() = Stages::default();
        let started = Instant::now();
        let outcome = {
            let _logging = ClassificationLogging::enter(options.classify_cache_hits);
            if options.monitor.enabled() || options.classify_cache_hits {
                let mut monitor = options.monitor.clone();
                if monitor.display == RustRedMonitorDisplay::Tui
                    && monitor.progress_snapshot.is_none()
                {
                    monitor.progress_snapshot = Some(options.work_dir.join("monitor.json"));
                }
                engine.evaluate_rustred_symbolic_with_monitor(&settings, input.as_view(), &monitor)
            } else {
                engine.evaluate_rustred_symbolic(&settings, input.as_view())
            }
        };
        let evaluation_seconds = started.elapsed().as_secs_f64();
        // The native flush is the durability acknowledgement on success and ordinary
        // reduction errors. Killed workers rely on the native periodic checkpoints.
        running_cache_receipt = running_cache.as_ref().map(|control| {
        let started = Instant::now();
        let flushed = control.flush();
        let flush_seconds = started.elapsed().as_secs_f64();
        let flush_complete = flushed
            .as_ref()
            .is_ok_and(|report| report.pending_records == 0 && report.writer_error.is_none());
        let (report, error) = match flushed {
            Ok(report) => (
                Some(report.to_document()),
                (!flush_complete)
                    .then(|| "native flush did not commit all accepted records".to_owned()),
            ),
            Err(error) => {
                eprintln!("running-cache flush failed: {error}");
                (
                    control.report().ok().map(|report| report.to_document()),
                    Some(error.to_string()),
                )
            }
        };
        let mut snapshot_seconds = None;
        let mut snapshot_sidecar = None;
        let mut snapshot_descriptor = None;
        let (snapshot, snapshot_error) = if flush_complete {
            let started = Instant::now();
            let frozen = (|| -> Result<Value> {
                let delta = fs::canonicalize(
                    running_cache_path
                        .as_ref()
                        .ok_or_else(|| eyre::eyre!("running cache has no delta path"))?,
                )?;
                let directory =
                    fs::canonicalize(&options.work_dir)?.join("running-cache-snapshots");
                ensure!(
                    !directory.starts_with(&delta),
                    "snapshot must be outside the running delta"
                );
                fs::create_dir_all(&directory)?;
                let directory = fs::canonicalize(directory)?;
                ensure!(
                    !directory.starts_with(&delta),
                    "snapshot directory resolves inside the running delta"
                );
                let stamp = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
                let sidecar =
                    directory.join(format!("measure-{}-{stamp}.rrindex", std::process::id()));
                let descriptor = sidecar.with_extension("json");
                ensure!(
                    !sidecar.exists() && !descriptor.exists(),
                    "preserve the existing running-cache snapshot"
                );
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
                Ok(snapshot) => (Some(snapshot), None),
                Err(error) => {
                    eprintln!("running-cache snapshot failed: {error}");
                    (None, Some(error.to_string()))
                }
            }
        } else {
            (None, None)
        };
        json!({
            "path":running_cache_path,"read_backing":running_read_backing,
            "capture_only":options.running_cache_capture_only,
            "read_delta":!options.running_cache_no_read_delta,
            "loading_seconds":running_cache_loading_seconds,"flush_seconds":flush_seconds,
            "report":report,"flush_error":error,"flush_complete":flush_complete,
            "snapshot":snapshot,"snapshot_error":snapshot_error,"snapshot_seconds":snapshot_seconds,
            "snapshot_sidecar":snapshot_sidecar,"snapshot_descriptor":snapshot_descriptor,
            "reduction_error":outcome.as_ref().err().map(ToString::to_string),
        })
    });
        let reduction_error_message = outcome.as_ref().err().map(ToString::to_string);
        let (outcome, mut reduction_error): (_, Option<eyre::Report>) = match outcome {
            Ok(result) => (Some(result), None),
            Err(error) => (None, Some(error.into())),
        };
        let mut preserve_reduction_error = |secondary: eyre::Report| {
            if running_cache_receipt.is_some()
                && let Some(primary) = reduction_error.take()
            {
                eprintln!("measurement reporting error after reduction failed: {secondary}");
                primary
            } else {
                secondary
            }
        };
        if let Some(receipt) = &running_cache_receipt {
            write_json(&options.work_dir.join("running-cache.json"), receipt)
                .map_err(&mut preserve_reduction_error)?;
        }
        let owners_after = engine
            .rustred_application_owner_preparation_statistics()
            .map_err(eyre::Report::from)
            .map_err(&mut preserve_reduction_error)?;
        let owner_seconds = owners_after
            .preparation_time
            .saturating_sub(owners_before.preparation_time)
            .as_secs_f64();
        let stages = LOGGER.stages();
        let saved_seconds = stages
            .seconds
            .get("saved_rule_application")
            .map_or(0.0, |samples| samples.iter().sum::<f64>());
        let statistics = engine
            .cache_statistics(5)
            .map_err(eyre::Report::from)
            .map_err(&mut preserve_reduction_error)?
            .saved_application
            .map(|s| s.to_document());
        let mut receipt = json!({
            "schema":"vakint.five-loop-cache-measurement.v1","identity":identity,
            "case_id":root.id,"variant":options.variant,"sample":options.sample,"input_sha256":format!("{:x}",Sha256::digest(input.to_canonical_string().as_bytes())),
            "declared_common_powers":root.powers,"rationale":root.rationale,
            "frozen_lowered_roots":root.lowered_roots,"frozen_protected_roots":root.protected_roots,
            "actual_lowered_roots":stages.cache_hits.input_roots,
            "provisional":options.provisional,
            "owner_proof":proof.as_ref().map(OwnerProof::identity),"owner_proof_validation_seconds":owner_proof_seconds,
            "owner_root_correlations":proof.as_ref().map(OwnerProof::correlations),
            "fresh_process":true,"data_loading_seconds":data_loading_seconds,"cache_loading_seconds":cache_loading_seconds,
            "external_data":data_report,"persistent_cache_attachments_before":caches_before,
            "persistent_cache_attachments_during_evaluation":caches_during_evaluation,
            "persistent_cache_policy":match (options.variant=="cached", running_read_backing) {
                (true,true)=>"manual-five-loop-and-running-cache",
                (false,true)=>"running-cache",
                (true,false)=>"manual-five-loop-only",
                (false,false)=>"none",
            },
            "evaluation_seconds":evaluation_seconds,"saved_application_including_owner_seconds":saved_seconds,
            "owner_preparation_seconds":owner_seconds,"ibp_excluding_owner_seconds":(saved_seconds-owner_seconds).max(0.0),
            // Retain the existing field name, but never label this derived interval
            // pure IBP: the native call also transports and combines exact terms.
            "native_timing_scope":{
                "saved_application_including_owner_seconds":"Native entry transport, lazy owner preparation, DAG reduction, exact sum, and terminal substitution",
                "ibp_excluding_owner_seconds":"Saved application interval minus actual deferred-owner preparation; entry transport, exact combination and terminal substitution remain included",
                "pure_ibp_seconds":null,"terminal_substitution_seconds":null,
                "separate_native_phase_durations_available":false
            },
            "owner_preparation_policy":"deferred on demand; actual stat delta, no common prewarm",
            "end_to_end_scope":"supervisor worker spawn through exit, including loading, evaluation and output serialization",
            "owner_preparation_during_load_seconds":owners_before.preparation_time.as_secs_f64(),
            "stages":stages.seconds,"lowering_receipts":stages.lowered,"warnings":stages.warnings,"cache_statistics":statistics,
        });
        if let Some(running_cache) = &running_cache_receipt {
            receipt["running_cache"] = running_cache.clone();
            receipt["running_cache_attached"] = running_cache_attached.into();
        }
        if options.classify_cache_hits {
            receipt["cache_hit_diagnostic"] = stages
                .cache_hits
                .save_hit_diagnostic(&options.work_dir, outcome.is_some());
            receipt["performance_claim_allowed"] = false.into();
        }
        if let Some(proof) = &proof {
            let evidence = proof.initial_routing_frontier(&root.id).unwrap_or_else(
                |error| json!({"complete":false,"scope":"unclassified","reason":error.to_string()}),
            );
            if !evidence.is_null() {
                receipt["owner_frontier_evidence"] = json!({
                    "scope":"Bound saved diagnostic only; not classification of this request or its timings",
                    "evidence":evidence,
                });
            }
        }
        let authentication = stages.authenticate(
            &root.lowered_roots,
            &root.protected_roots,
            options.provisional,
            proof.as_ref(),
            &root.id,
        );
        receipt["native_keys_authenticated"] = authentication
            .as_ref()
            .is_ok_and(|authenticated| *authenticated)
            .into();
        receipt["native_key_authentication"] = if !stages.cache_hits.input_roots_seen {
            "unavailable"
        } else if authentication
            .as_ref()
            .is_ok_and(|authenticated| *authenticated)
        {
            "verified-lowered-inputs-with-frozen-owner-protection"
        } else {
            "unverified"
        }
        .into();
        match outcome {
            Some(result) => {
                // Inspect the final combined expression against the instance's
                // declared identities. Retain success does not imply catalogue
                // closure, an independent master basis, or numerical availability.
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
                let unknown_path = options.work_dir.join("unknown-masters.json");
                write_json(&unknown_path, &json!(unknown))?;
                receipt["catalogue_closure"] = json!({
                    "closed":unknown.is_empty(),"master_count":masters.len(),
                    "unknown_count":unknown.len(),"unknown_masters":unknown,
                    "unknown_file":"unknown-masters.json","unknown_sha256":digest(&unknown_path)?,
                    "scope":"Final combined retained master identities versus loaded catalogue declarations",
                    "independent_basis_claim":false,"numerical_evaluation_claim":false
                });
                write_bytes(
                    &options.work_dir.join("value.symbolica"),
                    result.value.to_canonical_string().as_bytes(),
                )?;
                let guards = result
                    .nonzero_conditions
                    .iter()
                    .map(AtomCore::to_canonical_string)
                    .collect::<BTreeSet<_>>();
                write_json(&options.work_dir.join("guards.json"), &json!(guards))?;
                receipt["success"] = true.into();
                receipt["retained_result_available"] = true.into();
                receipt["value_sha256"] = digest(&options.work_dir.join("value.symbolica"))?.into();
                receipt["guards_sha256"] = digest(&options.work_dir.join("guards.json"))?.into();
            }
            None => {
                receipt["success"] = false.into();
                receipt["retained_result_available"] = false.into();
                receipt["catalogue_closure"] = Value::Null;
                receipt["error"] = json!(reduction_error_message);
            }
        }
        if let Err(error) = authentication {
            receipt["success"] = false.into();
            if running_cache_receipt.is_some() && !receipt["error"].is_null() {
                receipt["native_key_authentication_error"] = error.to_string().into();
            } else {
                receipt["error"] = error.to_string().into();
            }
        }
        if let Some(running_cache) = &running_cache_receipt {
            for phase in ["flush", "snapshot"] {
                if let Some(error) = running_cache[format!("{phase}_error")].as_str() {
                    receipt["success"] = false.into();
                    if receipt["error"].is_null() {
                        receipt["error"] = format!("running-cache {phase} failed: {error}").into();
                    }
                }
            }
        }
        if &options.identity().map_err(&mut preserve_reduction_error)? != identity {
            return Err(preserve_reduction_error(eyre::eyre!(
                "measurement inputs changed during evaluation"
            )));
        }
        write_json(&options.work_dir.join("worker.json"), &receipt)
            .map_err(&mut preserve_reduction_error)?;
        if running_cache_receipt.is_some()
            && let Some(error) = reduction_error
        {
            return Err(error);
        }
        ensure!(
            receipt["success"] == true,
            "public evaluation failed; inspect measurement.json"
        );
        Ok(())
    })();
    // Finalize the shared native writer even when validation or output fails.
    // The original evaluation error remains primary if persistence also fails.
    let close_error = running_cache.as_ref().and_then(|control| {
        let started = Instant::now();
        let closed = control.close();
        let document = match &closed {
            Ok(report) => Some(report.to_document()),
            Err(_) => control.report().ok().map(|report| report.to_document()),
        };
        let receipt = running_cache_receipt.get_or_insert_with(|| {
            json!({
                "path":running_cache_path,"read_backing":running_read_backing,
                "capture_only":options.running_cache_capture_only,
                "read_delta":!options.running_cache_no_read_delta,
                "loading_seconds":running_cache_loading_seconds,
            })
        });
        receipt["flush_report"] = receipt["report"].clone();
        receipt["report"] = document.clone().into();
        receipt["close"] = json!({"seconds":started.elapsed().as_secs_f64(),
            "report":document,"error":closed.as_ref().err().map(ToString::to_string)});
        receipt["terminal"] = true.into();
        closed.err()
    });
    let terminal_receipt = (|| -> Result<()> {
        if let Some(receipt) = &running_cache_receipt {
            write_json(&options.work_dir.join("running-cache.json"), receipt)?;
            let worker_path = options.work_dir.join("worker.json");
            if worker_path.exists() {
                let mut worker: Value = serde_json::from_slice(&fs::read(&worker_path)?)?;
                worker["running_cache"] = receipt.clone();
                if let Some(error) = &close_error {
                    worker["success"] = false.into();
                    if worker["error"].is_null() {
                        worker["error"] = json!(format!("running-cache close failed: {error}"));
                    }
                }
                write_json(&worker_path, &worker)?;
            }
        }
        Ok(())
    })();
    if let Err(error) = measured {
        if let Err(secondary) = terminal_receipt {
            eprintln!("running-cache terminal receipt failed: {secondary}");
        }
        return Err(error);
    }
    if let Some(error) = close_error {
        if let Err(secondary) = terminal_receipt {
            eprintln!("running-cache terminal receipt failed: {secondary}");
        }
        return Err(error.into());
    }
    terminal_receipt
}

fn main() -> Result<()> {
    let Some(options) = Options::parse()? else {
        return Ok(());
    };
    if options.worker {
        let _guard_admission = resource_guard_admission::admit_from_task_env()?;
    }
    let expected_identity = if options.worker {
        Some(options.arm_worker_watchdog()?)
    } else {
        ensure!(
            options.worker_deadline.is_none()
                && options.sample.is_none()
                && options.case_id.is_none(),
            "case arguments are reserved for supervised workers"
        );
        None
    };
    let workload_bytes = fs::read(&options.workload)?;
    let workload_hash = format!("{:x}", Sha256::digest(&workload_bytes));
    let workload: Workload = serde_json::from_slice(&workload_bytes)?;
    if !options.provisional {
        workload.validate_final_selection()?;
    } else {
        workload.validate()?;
    }
    ensure!(
        digest(&options.artifacts.join("manifest.json"))? == workload.data_manifest_sha256,
        "data differs from frozen workload"
    );
    ensure!(
        digest(&options.source_census)? == workload.source_census_sha256,
        "census differs from frozen workload"
    );
    let all_roots = if options.suite == "tuning" {
        &workload.tuning
    } else {
        &workload.held_out
    };
    let roots = filter_roots(all_roots, options.only_case.as_deref())?;
    ensure!(
        options.provisional
            || roots
                .iter()
                .all(|root| !root.lowered_roots.is_empty() && !root.protected_roots.is_empty()),
        "authenticated measurements require frozen lowered_roots and protected_roots; incomplete mapping is only permitted with --provisional"
    );
    let identity = options.identity()?;
    ensure!(
        identity["workload_sha256"] == workload_hash
            && identity["data_manifest_sha256"] == workload.data_manifest_sha256
            && identity["source_census_sha256"] == workload.source_census_sha256,
        "measurement inputs changed during preflight"
    );
    if options.worker {
        ensure!(
            expected_identity.as_ref() == Some(&identity),
            "worker inputs differ from the durable running identity"
        );
        let id = options
            .case_id
            .as_ref()
            .ok_or_else(|| eyre::eyre!("worker needs case id"))?;
        return measure(
            &options,
            roots
                .iter()
                .find(|root| &root.id == id)
                .ok_or_else(|| eyre::eyre!("case is outside frozen suite"))?,
            &identity,
            &workload,
        );
    }
    fs::create_dir_all(&options.work_dir)?;
    let lock = OpenOptions::new()
        .create(true)
        .truncate(false)
        .write(true)
        .open(options.work_dir.join("supervisor.lock"))?;
    lock.try_lock()
        .context("another suite supervisor is active")?;
    let identity_path = options.work_dir.join("identity.json");
    if identity_path.exists() {
        ensure!(
            serde_json::from_slice::<Value>(&fs::read(&identity_path)?)? == identity,
            "work directory has a different measurement identity"
        );
    } else {
        write_json(&identity_path, &identity)?;
    }
    let variants: &[&str] = if options.classify_cache_hits {
        &["cached"]
    } else if options.cache.is_some() {
        &["uncached", "cached"]
    } else {
        &["uncached"]
    };
    let mut measurements = Vec::new();
    for (case_index, root) in roots.iter().enumerate() {
        for sample in 0..options.repetitions {
            for variant in variants {
                let directory = options
                    .work_dir
                    .join(format!("case-{case_index:03}-{sample:02}-{variant}"));
                fs::create_dir_all(&directory)?;
                let receipt_path = directory.join("measurement.json");
                if !receipt_path.exists() {
                    let active_path = directory.join("running.json");
                    let now = SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs_f64();
                    let supervision = if active_path.exists() {
                        let active: Value = serde_json::from_slice(&fs::read(&active_path)?)?;
                        validate_receipt(&active, &identity, &root.id, variant, sample)?;
                        options.interrupted_supervision(&active, now)?
                    } else {
                        let started = Instant::now();
                        let wall_started = SystemTime::now();
                        let deadline = options.case_deadline(now)?;
                        write_json(
                            &active_path,
                            &json!({
                                "identity":identity,"case_id":root.id,"variant":variant,"sample":sample,
                                "deadline_unix_seconds":deadline,
                            }),
                        )?;
                        let log = File::create(directory.join("worker.log"))?;
                        let mut viewer = if options.monitor.display == RustRedMonitorDisplay::Tui {
                            let path = options
                                .monitor
                                .progress_snapshot
                                .clone()
                                .unwrap_or_else(|| directory.join("monitor.json"));
                            Some(RustRedSnapshotViewer::new(&path, options.monitor.no_color))
                        } else {
                            None
                        };
                        let mut command = Command::new(std::env::current_exe()?);
                        command
                            .args(std::env::args_os().skip(1))
                            .arg("--worker")
                            .arg("--case-id")
                            .arg(&root.id)
                            .arg("--variant")
                            .arg(variant)
                            .arg("--sample")
                            .arg(sample.to_string());
                        if let Some(deadline) = deadline {
                            command.arg("--worker-deadline").arg(deadline.to_string());
                        }
                        options.append_cache_file_arguments(&mut command);
                        let launched = command
                            .arg("--work-dir")
                            .arg(&directory)
                            .stdin(Stdio::null())
                            .stdout(log.try_clone()?)
                            .stderr(log)
                            .spawn();
                        match launched {
                            Ok(mut child) => {
                                let child_pid = child.id();
                                let supervisor_pid = std::process::id();
                                if let Some(viewer) = viewer.as_mut() {
                                    viewer.bind_worker_pid(child_pid);
                                }
                                let (status, expired) = loop {
                                    let expired = options.deadline_exhausted(
                                        deadline,
                                        started.elapsed(),
                                        SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs_f64(),
                                    )?;
                                    if let Some(status) = child.try_wait()? {
                                        break (status, expired);
                                    }
                                    if expired {
                                        child.kill()?;
                                        break (child.wait()?, true);
                                    }
                                    if let Some(viewer) = viewer.as_mut() {
                                        viewer.refresh();
                                    }
                                    std::thread::sleep(Duration::from_millis(200));
                                };
                                #[cfg(feature = "progress-tui")]
                                drop(viewer);
                                let wall_finished = SystemTime::now();
                                json!({"success":status.success() && !expired,
                                    "exit_code":status.code(),"exit_status":status.to_string(),
                                    "deadline_exhausted":expired || (options.max_case_seconds.is_some() && status.code()==Some(124)),
                                    "elapsed_seconds":started.elapsed().as_secs_f64(),
                                    "child_pid":child_pid,"supervisor_pid":supervisor_pid,
                                    "started_unix_seconds":wall_started.duration_since(UNIX_EPOCH)?.as_secs_f64(),
                                    "finished_unix_seconds":wall_finished.duration_since(UNIX_EPOCH)?.as_secs_f64(),
                                    "interval_scope":"supervisor pre-spawn setup through supervisor-observed worker completion"})
                            }
                            Err(error) => {
                                #[cfg(feature = "progress-tui")]
                                drop(viewer);
                                json!({"success":false,"launch_error":error.to_string(),
                                    "elapsed_seconds":started.elapsed().as_secs_f64()})
                            }
                        }
                    };
                    let worker_path = directory.join("worker.json");
                    let mut receipt = if worker_path.exists() {
                        let receipt: Value = serde_json::from_slice(&fs::read(worker_path)?)?;
                        validate_receipt(&receipt, &identity, &root.id, variant, sample)?;
                        receipt
                    } else {
                        json!({"identity":identity,"case_id":root.id,"variant":variant,"sample":sample,
                            "success":false,"error":"worker produced no receipt; retained without retry; inspect worker.log"})
                    };
                    receipt["worker_success"] = receipt["success"].clone();
                    receipt["success"] =
                        (receipt["success"] == true && supervision["success"] == true).into();
                    if supervision["success"] != true {
                        receipt["supervisor_error"] = "worker did not exit successfully; retained without retry; inspect supervision and worker.log".into();
                        if receipt["error"].is_null() {
                            receipt["error"] = receipt["supervisor_error"].clone();
                        }
                    }
                    receipt["end_to_end_seconds"] = supervision["elapsed_seconds"].clone();
                    receipt["supervision"] = supervision;
                    write_json(&receipt_path, &receipt)?;
                }
                let receipt: Value = serde_json::from_slice(&fs::read(&receipt_path)?)?;
                owner_proof::validate_completed_outputs(&receipt, &directory, root)?;
                validate_catalogue_closure(&receipt, &directory)?;
                validate_receipt(&receipt, &identity, &root.id, variant, sample)?;
                measurements.push(receipt);
            }
        }
    }
    let mut cases = Vec::new();
    for root in roots {
        let selected = measurements
            .iter()
            .filter(|receipt| receipt["case_id"] == root.id)
            .collect::<Vec<_>>();
        let successful = complete_samples(&selected, variants, options.repetitions)
            && selected.iter().all(|receipt| receipt["success"] == true);
        let exact = !options.classify_cache_hits
            && successful
            && selected.iter().all(|receipt| {
                receipt["value_sha256"] == selected[0]["value_sha256"]
                    && receipt["guards_sha256"] == selected[0]["guards_sha256"]
            });
        let mut row =
            json!({"id":root.id,"successful":successful,"exact_output_and_guard_parity":exact});
        if options.classify_cache_hits {
            // One cached diagnostic cannot compare its output with itself.
            // The mapper binds a prior accepted same-cohort uncached receipt.
            row["exact_output_and_guard_parity"] = Value::Null;
            row["hit_stream_complete"] = selected
                .iter()
                .all(|receipt| receipt["cache_hit_diagnostic"]["stream_complete"] == true)
                .into();
        }
        // Missing or failed samples cannot establish closure. A completed
        // retained expression may still contain correctly reported unknowns.
        row["catalogue_closure"] = if successful {
            let unknown = selected
                .iter()
                .map(|receipt| {
                    serde_json::from_value::<BTreeSet<String>>(
                        receipt["catalogue_closure"]["unknown_masters"].clone(),
                    )
                })
                .collect::<Result<Vec<_>, _>>()?
                .into_iter()
                .flatten()
                .collect::<BTreeSet<_>>();
            json!({"closed":unknown.is_empty(),"unknown_count":unknown.len(),
                "unknown_masters":unknown,"numerical_evaluation_claim":false})
        } else {
            Value::Null
        };
        row["native_keys_authenticated"] = (successful
            && selected
                .iter()
                .all(|receipt| receipt["native_keys_authenticated"] == true))
        .into();
        let root_hits = selected
            .iter()
            .filter(|receipt| receipt["variant"] == "cached")
            .filter_map(|receipt| {
                receipt["cache_statistics"]["persistent"]["requested_root_hits"].as_u64()
            })
            .sum::<u64>();
        let intermediate_hits = selected
            .iter()
            .filter(|receipt| receipt["variant"] == "cached")
            .filter_map(|receipt| {
                receipt["cache_statistics"]["persistent"]["intermediate_hits"].as_u64()
            })
            .sum::<u64>();
        row["requested_root_hits"] = root_hits.into();
        row["intermediate_hits"] = intermediate_hits.into();
        row["diagnostic_only"] = options.classify_cache_hits.into();
        row["backing_reuse"] = if root_hits > 0 {
            "incidental-requested-root-hit"
        } else if intermediate_hits > 0 {
            // The legacy Native counter includes Apply records for requested
            // roots. Only a complete recursive frontier proof can separate them.
            "root-scope-unclassified"
        } else {
            "no-backing-hit"
        }
        .into();
        for variant in variants {
            let mut samples = selected
                .iter()
                .filter(|receipt| receipt["variant"] == *variant)
                .filter_map(|receipt| receipt["ibp_excluding_owner_seconds"].as_f64())
                .collect::<Vec<_>>();
            samples.sort_by(f64::total_cmp);
            if successful
                && !options.classify_cache_hits
                && samples.len() == options.repetitions
                && samples
                    .iter()
                    .all(|sample| sample.is_finite() && *sample >= 0.0)
            {
                let middle = samples.len() / 2;
                let median = if samples.len() % 2 == 0 {
                    (samples[middle - 1] + samples[middle]) / 2.0
                } else {
                    samples[middle]
                };
                row[*variant] = json!({"median_seconds":median,"min_seconds":samples[0],"max_seconds":samples[samples.len()-1],"samples_seconds":samples});
            }
        }
        if exact
            && options.repetitions >= 3
            && options.cache.is_some()
            && let (Some(uncached), Some(cached)) = (
                row["uncached"]["median_seconds"].as_f64(),
                row["cached"]["median_seconds"].as_f64(),
            )
            && cached > 0.0
        {
            row["uncached_over_cached_ratio"] = (uncached / cached).into();
        }
        cases.push(row);
    }
    let complete = !cases.is_empty()
        && cases.iter().all(|case| {
            if options.classify_cache_hits {
                case["successful"] == true && case["hit_stream_complete"] == true
            } else {
                case["exact_output_and_guard_parity"] == true
            }
        });
    let ratios = cases
        .iter()
        .filter_map(|case| case["uncached_over_cached_ratio"].as_f64())
        .filter(|ratio| ratio.is_finite() && *ratio > 0.0)
        .collect::<Vec<_>>();
    let geometric_mean = if complete && ratios.len() == cases.len() {
        Some((ratios.iter().map(|ratio| ratio.ln()).sum::<f64>() / ratios.len() as f64).exp())
    } else {
        None
    };
    let authenticated = complete
        && !options.provisional
        && cases
            .iter()
            .all(|case| case["native_keys_authenticated"] == true);
    let provenance_complete = ["build_rustc", "build_profile", "build_git_revision"]
        .iter()
        .all(|field| {
            identity[field]
                .as_str()
                .is_some_and(|value| !value.trim().is_empty())
        })
        && identity["cpu_affinity"].is_string();
    write_json(
        &options.work_dir.join("summary.json"),
        &json!({"schema":"vakint.five-loop-cache-suite.v1","identity":identity,"complete":complete,"cases":cases,"measurements":measurements,"geometric_mean_uncached_over_cached_ratio":geometric_mean,"provisional":options.provisional,"diagnostic_only":options.classify_cache_hits,"native_keys_authenticated":authenticated,"execution_provenance_complete":provenance_complete,"final_held_out_comparison":!options.classify_cache_hits && authenticated && provenance_complete && options.suite=="held-out" && options.cache.is_some() && options.repetitions>=3,"selection_allowed":!options.classify_cache_hits && authenticated && provenance_complete && options.suite=="tuning" && options.cache.is_some() && options.repetitions>=3}),
    )?;
    ensure!(
        complete,
        "suite has failures or exact parity differences; inspect summary.json"
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parsed_options(extra: &[&str]) -> Result<Options> {
        let mut args = vec![
            "--artifacts",
            "data",
            "--workload",
            "workload.json",
            "--source-census",
            "census.json",
            "--work-dir",
            "measurement",
        ];
        args.extend_from_slice(extra);
        Options::parse_args(args.into_iter().map(str::to_owned))?
            .ok_or_else(|| eyre::eyre!("unexpected help request"))
    }

    #[test]
    fn cache_file_cli_rejects_orphans_and_preserves_explicit_zero_grants() -> Result<()> {
        for (flag, value) in [
            ("--cache-max-compressed-bytes", "0"),
            ("--cache-max-decoded-bytes", "1048576"),
            ("--cache-max-index-bytes", "0"),
            ("--cache-decoded-directory", "unopened-directory"),
        ] {
            let error = parsed_options(&[flag, value]).err().unwrap();
            assert!(error.to_string().contains("require --cache"));
        }
        let options = parsed_options(&[
            "--cache",
            "cache.rrcache.gz",
            "--cache-max-compressed-bytes",
            "0",
            "--cache-max-index-bytes",
            "0",
        ])?;
        let policy = options.cache_file_options.document()?;
        assert_eq!(policy["max_compressed_bytes"], 0);
        assert_eq!(policy["max_index_bytes"], 0);
        assert!(parsed_options(&["--cache", "cache.gz", "--cache-max-records", "1"]).is_err());
        Ok(())
    }

    #[test]
    fn cache_file_worker_arguments_preserve_policy_and_resume_identity() -> Result<()> {
        let scratch = tempfile::tempdir()?;
        fs::create_dir(scratch.path().join("decoded"))?;
        let directory = scratch.path().join("decoded/../decoded");
        let args = [
            "--artifacts",
            "data",
            "--workload",
            "workload.json",
            "--source-census",
            "census.json",
            "--work-dir",
            "measurement",
            "--cache",
            "cache.rrcache.gz",
            "--cache-max-compressed-bytes",
            "0",
            "--cache-max-decoded-bytes",
            "2147483648",
            "--cache-max-index-bytes",
            "0",
            "--cache-decoded-directory",
            directory.to_str().unwrap(),
        ]
        .map(str::to_owned);
        let parent = Options::parse_args(args.clone().into_iter())?.unwrap();
        let mut command = Command::new("unlaunched-collector");
        command.args(args).args([
            "--worker",
            "--case-id",
            "case",
            "--variant",
            "cached",
            "--sample",
            "0",
        ]);
        parent.append_cache_file_arguments(&mut command);
        let worker = Options::parse_args(
            command
                .get_args()
                .map(|arg| arg.to_str().unwrap().to_owned()),
        )?
        .unwrap();
        assert_eq!(parent.cache_file_options, worker.cache_file_options);
        let policy = parent.cache_file_options.document()?;
        assert_eq!(policy, worker.cache_file_options.document()?);
        assert_eq!(
            policy["decoded_directory"],
            serde_json::to_value(scratch.path().join("decoded").canonicalize()?)?
        );
        let identity = json!({"cache_file_options":policy});
        let receipt = json!({"identity":identity,"case_id":"case","variant":"cached","sample":0});
        validate_receipt(&receipt, &identity, "case", "cached", 0)?;
        let mut changed = identity.clone();
        changed["cache_file_options"]["max_index_bytes"] = 1.into();
        assert!(validate_receipt(&receipt, &changed, "case", "cached", 0).is_err());
        let mut omitted = Command::new("unlaunched-collector");
        parsed_options(&[])?.append_cache_file_arguments(&mut omitted);
        assert_eq!(omitted.get_args().count(), 0);
        Ok(())
    }

    #[test]
    fn collector_defaults_to_unbounded_non_diagnostic_execution() -> Result<()> {
        let options = parsed_options(&[])?;
        assert_eq!(options.repetitions, 3);
        assert_eq!(options.max_case_seconds, None);
        assert_eq!(options.case_deadline(100.0)?, None);
        assert_eq!(options.worker_allowance(100.0)?, None);
        assert!(!options.classify_cache_hits);
        assert!(!options.monitor.enabled());
        #[cfg(feature = "experimental-rustred")]
        assert!(!options.monitor.cache_hit_diagnostic);
        assert!(!options.deadline_exhausted(None, Duration::from_secs(u64::MAX), f64::MAX)?);
        Ok(())
    }

    #[test]
    fn explicit_positive_budget_bounds_both_clocks_and_worker_allowance() -> Result<()> {
        let mut options = parsed_options(&["--max-case-seconds", "30"])?;
        assert_eq!(options.max_case_seconds, Some(30));
        let deadline = options.case_deadline(100.0)?;
        assert_eq!(deadline, Some(130.0));
        options.worker_deadline = deadline;
        assert_eq!(
            options.worker_allowance(110.0)?,
            Some(Duration::from_secs(20))
        );
        assert_eq!(
            options.worker_allowance(90.0)?,
            Some(Duration::from_secs(30))
        );
        assert!(!options.deadline_exhausted(deadline, Duration::from_secs(29), 129.0)?);
        assert!(options.deadline_exhausted(deadline, Duration::from_secs(30), 100.0)?);
        assert!(options.deadline_exhausted(deadline, Duration::ZERO, 130.0)?);
        assert!(options.worker_allowance(130.0).is_err());
        Ok(())
    }

    #[test]
    fn cli_rejects_invalid_budget_and_inconsistent_worker_deadlines() {
        for args in [
            vec!["--max-case-seconds", "0"],
            vec!["--max-case-seconds", "-1"],
            vec!["--max-case-seconds", "NaN"],
            vec!["--max-case-seconds"],
            vec!["--worker", "--max-case-seconds", "30"],
            vec!["--worker", "--worker-deadline", "130"],
            vec![
                "--worker",
                "--max-case-seconds",
                "30",
                "--worker-deadline",
                "NaN",
            ],
            vec![
                "--worker",
                "--max-case-seconds",
                "30",
                "--worker-deadline",
                "inf",
            ],
        ] {
            assert!(
                parsed_options(&args).is_err(),
                "accepted invalid arguments {args:?}"
            );
        }
    }

    #[test]
    fn deadline_policy_rejects_unpaired_or_nonfinite_state() -> Result<()> {
        let mut options = parsed_options(&[])?;
        options.worker_deadline = Some(130.0);
        assert!(options.worker_allowance(100.0).is_err());
        assert!(
            options
                .deadline_exhausted(Some(130.0), Duration::ZERO, 100.0)
                .is_err()
        );
        options.max_case_seconds = Some(30);
        options.worker_deadline = None;
        assert!(options.worker_allowance(100.0).is_err());
        assert!(
            options
                .deadline_exhausted(None, Duration::ZERO, 100.0)
                .is_err()
        );
        assert!(
            options
                .deadline_exhausted(Some(f64::NAN), Duration::ZERO, 100.0)
                .is_err()
        );
        assert!(options.case_deadline(f64::INFINITY).is_err());
        Ok(())
    }

    #[test]
    fn interrupted_unbounded_attempt_never_infers_expiry_or_retries() -> Result<()> {
        let options = parsed_options(&[])?;
        for active in [
            json!({"deadline_unix_seconds":null}),
            json!({"deadline_unix_seconds":1}),
        ] {
            let error = options
                .interrupted_supervision(&active, f64::MAX)
                .err()
                .unwrap();
            assert!(
                error
                    .to_string()
                    .contains("terminal completion is independently verified")
            );
        }
        Ok(())
    }

    #[test]
    fn interrupted_bounded_attempt_preserves_existing_no_retry_failure() -> Result<()> {
        let options = parsed_options(&["--max-case-seconds", "30"])?;
        let active = json!({"deadline_unix_seconds":130.0});
        assert!(options.interrupted_supervision(&active, 129.0).is_err());
        let supervision = options.interrupted_supervision(&active, 130.0)?;
        assert_eq!(supervision["success"], false);
        assert_eq!(supervision["interrupted"], true);
        assert_eq!(supervision["deadline_exhausted"], true);
        assert_eq!(
            supervision["error"],
            "interrupted attempt retained without retry"
        );
        assert!(options.interrupted_supervision(&json!({}), 200.0).is_err());
        Ok(())
    }

    #[test]
    fn unbounded_worker_requires_matching_durable_active_identity() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let mut options = parsed_options(&["--worker", "--case-id", "case", "--sample", "0"])?;
        options.work_dir = directory.path().into();
        let mut active = json!({"identity":{"frozen":"identity"},"case_id":"case",
            "variant":"uncached","sample":0,"deadline_unix_seconds":null});
        write_json(&directory.path().join("running.json"), &active)?;
        assert_eq!(options.arm_worker_watchdog()?, active["identity"]);
        active
            .as_object_mut()
            .unwrap()
            .remove("deadline_unix_seconds");
        write_json(&directory.path().join("running.json"), &active)?;
        assert!(options.arm_worker_watchdog().is_err());
        active["deadline_unix_seconds"] = Value::Null;
        active["case_id"] = "different".into();
        write_json(&directory.path().join("running.json"), &active)?;
        assert!(options.arm_worker_watchdog().is_err());
        Ok(())
    }

    #[cfg(feature = "experimental-rustred")]
    #[test]
    fn classification_cli_explicitly_opts_into_lossless_monitor_events() -> Result<()> {
        let options = parsed_options(&[
            "--classify-cache-hits",
            "--cache",
            "native-cache",
            "--repetitions",
            "1",
        ])?;
        assert!(options.classify_cache_hits);
        assert!(options.monitor.cache_hit_diagnostic);
        assert!(!options.monitor.enabled());
        assert_eq!(options.max_case_seconds, None);
        for args in [
            vec!["--classify-cache-hits", "--repetitions", "1"],
            vec!["--classify-cache-hits", "--cache", "native-cache"],
            vec![
                "--classify-cache-hits",
                "--cache",
                "native-cache",
                "--repetitions",
                "1",
                "--progress-snapshot",
                "progress.json",
            ],
            vec![
                "--classify-cache-hits",
                "--cache",
                "native-cache",
                "--repetitions",
                "1",
                "--running-cache",
                "delta",
            ],
        ] {
            assert!(
                parsed_options(&args).is_err(),
                "accepted invalid diagnostic arguments {args:?}"
            );
        }
        Ok(())
    }

    #[cfg(not(feature = "experimental-rustred"))]
    #[test]
    fn classification_cli_rejects_builds_without_experimental_diagnostics() {
        assert!(
            parsed_options(&[
                "--classify-cache-hits",
                "--cache",
                "native-cache",
                "--repetitions",
                "1"
            ])
            .is_err()
        );
    }

    #[test]
    fn lossless_hit_buffer_preserves_multiplicity_and_dispatch_boundaries() -> Result<()> {
        let document = owner_proof::hit_diagnostic_fixture();
        let mut stages = Stages::default();
        stages
            .cache_hits
            .capture_input_roots(&json!({"phase":"input_roots","roots":[]}).to_string())?;
        for event in document["events"].as_array().unwrap() {
            stages.cache_hits.capture_hit_event(&event.to_string())?;
        }
        assert_eq!(stages.cache_hits.hit_diagnostic(true), {
            let mut expected = document;
            expected["diagnostic_only"] = true.into();
            expected["timing_scope"] =
                "Lossless observation overhead; not an unobserved performance sample".into();
            expected
        });
        assert!(
            owner_proof::validate_hit_diagnostic(&stages.cache_hits.hit_diagnostic(true)).is_ok()
        );
        let hit = stages
            .cache_hits
            .hit_events
            .iter()
            .find(|event| event["event"] == "backing_hit")
            .unwrap()
            .clone();
        stages.cache_hits.capture_hit_event(&hit.to_string())?;
        // Duplicate delivery is retained and fails reconciliation, not deduplicated away.
        assert_eq!(
            stages
                .cache_hits
                .hit_events
                .iter()
                .filter(|event| event["event"] == "backing_hit")
                .count(),
            2
        );
        assert!(
            owner_proof::validate_hit_diagnostic(&stages.cache_hits.hit_diagnostic(true)).is_err()
        );
        Ok(())
    }

    #[test]
    fn missing_dispatch_and_cancelled_reduction_never_become_complete() -> Result<()> {
        let document = owner_proof::hit_diagnostic_fixture();
        let mut stages = Stages::default();
        stages.cache_hits.hit_events = document["events"].as_array().unwrap().clone();
        assert!(
            owner_proof::validate_hit_diagnostic(&stages.cache_hits.hit_diagnostic(true)).is_err()
        );
        stages
            .cache_hits
            .capture_input_roots(&json!({"phase":"input_roots","roots":[]}).to_string())?;
        assert!(
            owner_proof::validate_hit_diagnostic(&stages.cache_hits.hit_diagnostic(false)).is_err()
        );
        stages.cache_hits.hit_error = Some("a telemetry record could not be retained".into());
        assert!(
            owner_proof::validate_hit_diagnostic(&stages.cache_hits.hit_diagnostic(true)).is_err()
        );
        Ok(())
    }

    #[test]
    fn collected_seeds_are_not_reconstructed_from_the_unweighted_input_union() -> Result<()> {
        let document = owner_proof::hit_diagnostic_fixture();
        let first = vec![1_i64; 15];
        let second = vec![2_i64; 15];
        let mut stages = Stages::default();
        stages.cache_hits.capture_input_roots(
            &json!({"phase":"input_roots","roots":[first.clone(),second]}).to_string(),
        )?;
        stages.cache_hits.hit_events = document["events"].as_array().unwrap().clone();
        let observed =
            owner_proof::validate_hit_diagnostic(&stages.cache_hits.hit_diagnostic(true))?;
        // The raw logger cannot encode coefficient cancellation. Only the
        // native collected-root stream supplies the surviving combined seeds.
        assert_eq!(stages.cache_hits.input_roots.len(), 2);
        assert_eq!(
            observed["dispatches"][0]["seeds"].as_array().unwrap().len(),
            1
        );
        assert_eq!(
            observed["dispatches"][0]["seeds"][0]["powers"],
            json!(first)
        );
        stages.cache_hits.input_dispatch_count = u64::MAX;
        assert!(
            stages
                .cache_hits
                .capture_input_roots(&json!({"phase":"input_roots","roots":[]}).to_string())
                .is_err()
        );
        assert_eq!(stages.cache_hits.input_dispatch_count, u64::MAX);
        Ok(())
    }

    #[test]
    fn hit_evidence_commit_failure_keeps_classification_unavailable() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let mut stages = Stages::default();
        let document = owner_proof::hit_diagnostic_fixture();
        stages.cache_hits.hit_events = document["events"].as_array().unwrap().clone();
        stages.cache_hits.input_dispatch_count = 1;
        let committed = stages
            .cache_hits
            .save_hit_diagnostic(directory.path(), true);
        assert_eq!(committed["stream_complete"], true);
        assert_eq!(committed["classification_available"], false);
        assert_eq!(
            committed["sha256"],
            digest(&directory.path().join("cache-hit-events.json"))?
        );
        let blocker = directory.path().join("not-a-directory");
        fs::write(&blocker, b"preserve prior data")?;
        let failed = stages.cache_hits.save_hit_diagnostic(&blocker, true);
        assert_eq!(failed["stream_complete"], false);
        assert_eq!(failed["classification_available"], false);
        assert!(failed["sha256"].is_null());
        assert!(
            failed["error"]
                .as_str()
                .unwrap()
                .contains("could not be committed")
        );
        assert_eq!(fs::read(&blocker)?, b"preserve prior data");
        assert_eq!(stages.cache_hits.hit_diagnostic(true)["success"], true);
        Ok(())
    }

    #[test]
    fn classification_log_admission_is_released_on_scope_exit() {
        assert!(!ClassificationLogging::enabled());
        {
            let _guard = ClassificationLogging::enter(true);
            assert!(ClassificationLogging::enabled());
        }
        assert!(!ClassificationLogging::enabled());
    }

    #[test]
    fn cache_file_identity_preserves_the_existing_digest() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let path = directory.path().join("cache.rrcache.gz");
        fs::write(&path, b"existing opaque single-file payload")?;
        assert_eq!(
            CacheIdentity::read(&path)?,
            CacheIdentity::File(digest(&path)?)
        );
        Ok(())
    }

    #[cfg(unix)]
    #[test]
    fn cache_file_identity_preserves_existing_symlink_support() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let path = directory.path().join("cache.rrcache.gz");
        let link = directory.path().join("cache-link.rrcache.gz");
        fs::write(&path, b"existing opaque single-file payload")?;
        std::os::unix::fs::symlink(&path, &link)?;
        assert_eq!(CacheIdentity::read(&path)?, CacheIdentity::read(&link)?);
        Ok(())
    }

    #[test]
    fn native_pack_identity_reads_only_the_index_and_detects_changes() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let index = directory.path().join("cache.rrindex");
        fs::write(&index, b"opaque native index generation one")?;
        fs::create_dir(directory.path().join("chunks"))?;
        fs::write(
            directory.path().join("chunks/unused"),
            b"not decoded or hashed",
        )?;
        let before = CacheIdentity::read(directory.path())?;
        let CacheIdentity::Native(document) = &before else {
            panic!("expected native identity")
        };
        assert_eq!(document["layout"], "pack");
        assert_eq!(document["indexes"].as_object().unwrap().len(), 1);
        fs::write(directory.path().join("chunks/unused"), b"still unused")?;
        assert_eq!(before, CacheIdentity::read(directory.path())?);
        fs::write(&index, b"opaque native index generation two")?;
        assert_ne!(before, CacheIdentity::read(directory.path())?);
        Ok(())
    }

    #[test]
    fn native_bundle_identity_is_sorted_relocatable_and_binds_child_indexes() -> Result<()> {
        let directories = [tempfile::tempdir()?, tempfile::tempdir()?];
        let children = ["a".repeat(64), "b".repeat(64)];
        for (number, directory) in directories.iter().enumerate() {
            fs::write(
                directory.path().join("bundle.rrindex"),
                b"opaque namespace index",
            )?;
            for child in if number == 0 {
                children.iter().collect::<Vec<_>>()
            } else {
                children.iter().rev().collect()
            } {
                let path = directory.path().join("bindings").join(child);
                fs::create_dir_all(&path)?;
                fs::write(path.join("cache.rrindex"), child)?;
            }
        }
        let before = CacheIdentity::read(directories[0].path())?;
        assert_eq!(before, CacheIdentity::read(directories[1].path())?);
        let CacheIdentity::Native(document) = &before else {
            panic!("expected native identity")
        };
        assert_eq!(document["layout"], "bundle");
        assert_eq!(document["indexes"].as_object().unwrap().len(), 3);
        fs::write(
            directories[0]
                .path()
                .join("bindings")
                .join(&children[0])
                .join("cache.rrindex"),
            b"appended child generation",
        )?;
        assert_ne!(before, CacheIdentity::read(directories[0].path())?);
        Ok(())
    }

    #[test]
    fn native_identity_rejects_incomplete_or_ambiguous_layouts() -> Result<()> {
        let directory = tempfile::tempdir()?;
        assert!(CacheIdentity::read(directory.path()).is_err());
        fs::write(directory.path().join("bundle.rrindex"), b"opaque")?;
        assert!(CacheIdentity::read(directory.path()).is_err());
        let bindings = directory.path().join("bindings");
        fs::create_dir(&bindings)?;
        assert!(CacheIdentity::read(directory.path()).is_err());
        let invalid = bindings.join("A".repeat(64));
        fs::create_dir(&invalid)?;
        fs::write(invalid.join("cache.rrindex"), b"opaque")?;
        assert!(CacheIdentity::read(directory.path()).is_err());
        fs::remove_dir_all(invalid)?;
        let child = bindings.join("a".repeat(64));
        fs::create_dir(&child)?;
        assert!(CacheIdentity::read(directory.path()).is_err());
        fs::write(child.join("cache.rrindex"), b"opaque")?;
        assert!(CacheIdentity::read(directory.path()).is_ok());
        fs::write(directory.path().join("cache.rrindex"), b"ambiguous")?;
        assert!(CacheIdentity::read(directory.path()).is_err());
        Ok(())
    }

    #[cfg(unix)]
    #[test]
    fn native_identity_refuses_symlinked_layout_entries_without_fallback() -> Result<()> {
        use std::os::unix::fs::symlink;
        let directory = tempfile::tempdir()?;
        let outside = tempfile::tempdir()?;
        fs::write(outside.path().join("cache.rrindex"), b"opaque")?;
        let root_link = directory.path().join("root-link");
        symlink(outside.path(), &root_link)?;
        assert!(CacheIdentity::read(&root_link).is_err());
        let pack = directory.path().join("cache.rrindex");
        symlink(outside.path().join("cache.rrindex"), &pack)?;
        assert!(CacheIdentity::read(directory.path()).is_err());
        fs::remove_file(&pack)?;
        fs::write(&pack, b"opaque")?;
        symlink(
            directory.path().join("missing-index"),
            directory.path().join("bundle.rrindex"),
        )?;
        assert!(CacheIdentity::read(directory.path()).is_err());
        fs::remove_file(&pack)?;
        fs::remove_file(directory.path().join("bundle.rrindex"))?;
        fs::write(directory.path().join("bundle.rrindex"), b"opaque")?;
        symlink(outside.path(), directory.path().join("bindings"))?;
        assert!(CacheIdentity::read(directory.path()).is_err());
        fs::remove_file(directory.path().join("bindings"))?;
        fs::create_dir(directory.path().join("bindings"))?;
        let child = directory.path().join("bindings").join("a".repeat(64));
        symlink(outside.path(), &child)?;
        assert!(CacheIdentity::read(directory.path()).is_err());
        fs::remove_file(&child)?;
        fs::create_dir(&child)?;
        symlink(
            outside.path().join("cache.rrindex"),
            child.join("cache.rrindex"),
        )?;
        assert!(CacheIdentity::read(directory.path()).is_err());
        Ok(())
    }

    #[test]
    fn retained_unknowns_remain_successful_but_cannot_claim_catalogue_closure() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let unknown = json!(["vakint::RustRedMaster(vakint::rustred_family_unknown,2)"]);
        write_json(&directory.path().join("unknown-masters.json"), &unknown)?;
        let mut receipt = json!({"success":true,"catalogue_closure":{
            "closed":false,"master_count":1,"unknown_count":1,"unknown_masters":unknown,
            "unknown_sha256":digest(&directory.path().join("unknown-masters.json"))?}});
        validate_catalogue_closure(&receipt, directory.path())?;
        assert_eq!(receipt["success"], true);
        receipt["catalogue_closure"]["closed"] = true.into();
        assert!(validate_catalogue_closure(&receipt, directory.path()).is_err());
        receipt["catalogue_closure"]["closed"] = false.into();
        write_json(&directory.path().join("unknown-masters.json"), &json!([]))?;
        assert!(validate_catalogue_closure(&receipt, directory.path()).is_err());
        Ok(())
    }

    #[test]
    fn complete_cases_require_every_sample_and_variant() {
        let variants = ["uncached", "cached"];
        assert!(!complete_samples(&[], &variants, 1));
        let a = json!({"variant":"uncached","sample":0});
        let b = json!({"variant":"cached","sample":0});
        assert!(complete_samples(&[&a, &b], &variants, 1));
        assert!(!complete_samples(&[&a], &variants, 1));
        assert!(!complete_samples(&[&a, &b, &a, &b], &variants, 2));
        let c = json!({"variant":"uncached","sample":1});
        let d = json!({"variant":"cached","sample":1});
        assert!(complete_samples(&[&a, &b, &c, &d], &variants, 2));
    }

    #[test]
    fn resumed_receipts_must_match_the_case_variant_and_sample() {
        let identity = json!({"workload_sha256":"frozen"});
        let receipt = json!({"identity":identity,"case_id":"case-a","variant":"cached","sample":1});
        validate_receipt(&receipt, &identity, "case-a", "cached", 1).unwrap();
        assert!(validate_receipt(&receipt, &identity, "case-b", "cached", 1).is_err());
        assert!(validate_receipt(&receipt, &identity, "case-a", "uncached", 1).is_err());
        assert!(validate_receipt(&receipt, &identity, "case-a", "cached", 0).is_err());
    }

    #[test]
    fn public_input_requires_the_complete_frozen_lowered_union_and_owner_proof() {
        let workload = workload::tests::selected();
        let root = &workload.held_out[0];
        let first = root.lowered_roots[0].clone();
        let second = workload.held_out[1].lowered_roots[0].clone();
        let proof = owner_proof::tests::fixture();
        let mut stages = Stages::default();
        assert!(
            stages
                .authenticate(
                    &[first.clone()],
                    &[first.clone()],
                    false,
                    Some(&proof),
                    &root.id
                )
                .is_err()
        );
        assert!(!stages.authenticate(&[], &[], true, None, &root.id).unwrap());
        stages
            .cache_hits
            .capture_input_roots(
                &json!({"phase":"input_roots","roots":[first.clone()]}).to_string(),
            )
            .unwrap();
        assert!(
            stages
                .authenticate(&[second], &[first.clone()], false, Some(&proof), &root.id)
                .is_err()
        );
        assert!(
            stages
                .authenticate(&[first.clone()], &[first.clone()], false, None, &root.id)
                .is_err()
        );
        assert!(
            !stages
                .authenticate(&[first.clone()], &[first.clone()], true, None, &root.id)
                .unwrap()
        );
        assert!(
            stages
                .authenticate(&[first.clone()], &[first], false, Some(&proof), &root.id)
                .unwrap()
        );
    }
}
