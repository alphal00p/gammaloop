//! Fresh-process public-evaluator measurements of a frozen cache suite.
#[path = "build_cache/workload.rs"]
mod workload;

use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::Mutex,
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};

use eyre::{Context, Result, bail, ensure};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use symbolica::atom::AtomCore;
use vakint::{
    EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions, RustRedMasterMode, Vakint,
    VakintSettings,
};
use workload::{Root, Workload};

#[derive(Default)]
struct Stages {
    seconds: BTreeMap<String, Vec<f64>>,
    lowered: Vec<String>,
    warnings: Vec<String>,
    input_roots: BTreeSet<Vec<i64>>,
    input_roots_seen: bool,
    input_roots_error: Option<String>,
}

impl Stages {
    fn capture_input_roots(&mut self, message: &str) -> Result<()> {
        let message: Value = serde_json::from_str(message)?;
        ensure!(
            message["phase"] == "input_roots",
            "unknown input-root telemetry phase"
        );
        let roots: Vec<Vec<i64>> = serde_json::from_value(message["roots"].clone())?;
        ensure!(
            roots.iter().all(|root| root.len() == 15),
            "invalid input-root telemetry arity"
        );
        self.input_roots_seen = true;
        self.input_roots.extend(roots);
        Ok(())
    }

    fn authenticate(
        &self,
        lowered: &[Vec<i64>],
        protected: &[Vec<i64>],
        provisional: bool,
    ) -> Result<bool> {
        ensure!(
            self.input_roots_error.is_none(),
            "invalid input-root telemetry: {:?}",
            self.input_roots_error
        );
        if !self.input_roots_seen || lowered.is_empty() || protected.is_empty() {
            ensure!(
                provisional,
                "native key authentication unavailable: require complete telemetry and frozen protected_roots, or use --provisional without final claims"
            );
            return Ok(false);
        }
        let lowered = lowered.iter().cloned().collect::<BTreeSet<_>>();
        ensure!(
            self.input_roots == lowered,
            "public input lowered input roots differ from frozen lowered_roots"
        );
        Ok(true)
    }
}

struct StageLogger(Mutex<Stages>);
static LOGGER: StageLogger = StageLogger(Mutex::new(Stages {
    seconds: BTreeMap::new(),
    lowered: Vec::new(),
    warnings: Vec::new(),
    input_roots: BTreeSet::new(),
    input_roots_seen: false,
    input_roots_error: None,
}));

impl log::Log for StageLogger {
    fn enabled(&self, metadata: &log::Metadata<'_>) -> bool {
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
        if let Ok(mut stages) = self.0.lock() {
            if record.target() == "vakint::cache_workload" {
                if let Err(error) = stages.capture_input_roots(&message) {
                    stages.input_roots_error = Some(error.to_string());
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
    work_dir: PathBuf,
    suite: String,
    repetitions: usize,
    max_case_seconds: u64,
    case_id: Option<String>,
    variant: String,
    worker: bool,
    sample: Option<usize>,
    worker_deadline: Option<f64>,
    provisional: bool,
}

impl Options {
    fn arm_worker_watchdog(&self) -> Result<Value> {
        let deadline = self
            .worker_deadline
            .ok_or_else(|| eyre::eyre!("worker requires an absolute supervisor deadline"))?;
        let now = SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs_f64();
        let remaining = (deadline - now).min(self.max_case_seconds as f64);
        ensure!(
            deadline.is_finite() && remaining.is_finite() && remaining > 0.0,
            "worker deadline has expired or is invalid"
        );
        std::thread::Builder::new()
            .name("measurement-deadline".into())
            .spawn(move || {
                let started = Instant::now();
                let allowance = Duration::from_secs_f64(remaining);
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
            .context("start independent measurement watchdog")?;
        let active: Value = serde_json::from_slice(&fs::read(self.work_dir.join("running.json"))?)?;
        ensure!(
            active["deadline_unix_seconds"].as_f64() == Some(deadline)
                && active["case_id"].as_str() == self.case_id.as_deref()
                && active["variant"] == self.variant
                && active["sample"].as_u64() == self.sample.map(|sample| sample as u64),
            "worker has no matching durable active case"
        );
        Ok(active["identity"].clone())
    }

    fn parse() -> Result<Option<Self>> {
        let mut artifacts = None;
        let mut workload = None;
        let mut source_census = None;
        let mut cache = None;
        let mut work_dir = None;
        let mut suite = "tuning".to_owned();
        let mut repetitions = 3;
        let mut max_case_seconds = 300;
        let mut case_id = None;
        let mut variant = "uncached".to_owned();
        let mut worker = false;
        let mut sample = None;
        let mut worker_deadline = None;
        let mut provisional = false;
        let mut args = std::env::args().skip(1);
        while let Some(arg) = args.next() {
            if matches!(arg.as_str(), "--help" | "-h") {
                println!(
                    "measure_cache --artifacts PATH --workload JSON --source-census FILE --work-dir PATH\n\
                     [--cache FILE.gz] [--suite tuning|held-out] [--repetitions N] [--max-case-seconds N]\n\
                     [--provisional] permits legacy functional pilots without native key authentication.\n\
                     Every case/sample/variant runs in a fresh process. Default: three samples, 300s each.\n\
                     Omit --cache for an uncached baseline; final cache comparisons require three samples.\n\
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
            let value = args
                .next()
                .ok_or_else(|| eyre::eyre!("missing value for {arg}"))?;
            match arg.as_str() {
                "--artifacts" => artifacts = Some(value.into()),
                "--workload" => workload = Some(value.into()),
                "--source-census" => source_census = Some(value.into()),
                "--cache" => cache = Some(value.into()),
                "--work-dir" => work_dir = Some(value.into()),
                "--suite" => suite = value,
                "--repetitions" => repetitions = value.parse()?,
                "--max-case-seconds" => max_case_seconds = value.parse()?,
                "--case-id" => case_id = Some(value),
                "--variant" => variant = value,
                "--sample" => sample = Some(value.parse()?),
                "--worker-deadline" => worker_deadline = Some(value.parse()?),
                _ => bail!("unknown argument {arg}"),
            }
        }
        ensure!(
            matches!(suite.as_str(), "tuning" | "held-out"),
            "unknown suite"
        );
        ensure!(
            repetitions > 0 && max_case_seconds > 0,
            "positive repetition/deadline required"
        );
        ensure!(
            matches!(variant.as_str(), "uncached" | "cached"),
            "unknown variant"
        );
        ensure!(
            variant != "cached" || cache.is_some(),
            "cached worker needs --cache"
        );
        Ok(Some(Self {
            artifacts: artifacts.ok_or_else(|| eyre::eyre!("--artifacts is required"))?,
            workload: workload.ok_or_else(|| eyre::eyre!("--workload is required"))?,
            source_census: source_census
                .ok_or_else(|| eyre::eyre!("--source-census is required"))?,
            work_dir: work_dir.ok_or_else(|| eyre::eyre!("--work-dir is required"))?,
            cache,
            suite,
            repetitions,
            max_case_seconds,
            case_id,
            variant,
            worker,
            sample,
            worker_deadline,
            provisional,
        }))
    }

    fn identity(&self) -> Result<Value> {
        Ok(json!({
            "workload_sha256":digest(&self.workload)?,
            "source_census_sha256":digest(&self.source_census)?,
            "data_manifest_sha256":digest(&self.artifacts.join("manifest.json"))?,
            "cache_sha256":self.cache.as_ref().map(|path|digest(path)).transpose()?,
            "executable_sha256":digest(&std::env::current_exe()?)?,
            "compiled_driver_sha256":format!("{:x}",Sha256::digest(include_bytes!("measure_cache.rs"))),
            "compiled_workload_schema_sha256":format!("{:x}",Sha256::digest(include_bytes!("build_cache/workload.rs"))),
            "compiled_cargo_lock_sha256":format!("{:x}",Sha256::digest(include_bytes!("../../../Cargo.lock"))),
            "build_rustc":option_env!("VAKINT_CACHE_BUILD_RUSTC"),
            "build_profile":option_env!("VAKINT_CACHE_BUILD_PROFILE"),
            "build_git_revision":option_env!("VAKINT_CACHE_BUILD_GIT_REVISION"),
            "suite":self.suite,"repetitions":self.repetitions,"max_case_seconds":self.max_case_seconds,
            "provisional":self.provisional,
            "cpu_affinity":fs::read_to_string("/proc/self/status").ok().and_then(|status|status.lines().find_map(|line|line.strip_prefix("Cpus_allowed_list:").map(str::trim).map(str::to_owned))),
            "rayon_num_threads":std::env::var("RAYON_NUM_THREADS").ok(),
        }))
    }
}

fn digest(path: &Path) -> Result<String> {
    let mut input = File::open(path)?;
    let mut digest = Sha256::new();
    let mut block = [0; 65536];
    loop {
        let count = input.read(&mut block)?;
        if count == 0 {
            break;
        }
        digest.update(&block[..count]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

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

fn measure(options: &Options, root: &Root, identity: &Value) -> Result<()> {
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
    engine.load_external_data(&options.artifacts, 5)?;
    engine.detach_cache(5)?;
    let data_loading_seconds = started.elapsed().as_secs_f64();
    let cache_started = Instant::now();
    if options.variant == "cached" {
        engine.load_cache(5, options.cache.as_ref().unwrap(), &Default::default())?;
        engine.clear_cache(5)?;
    }
    let cache_loading_seconds = cache_started.elapsed().as_secs_f64();
    let owners_before = engine.rustred_application_owner_preparation_statistics()?;
    *LOGGER
        .0
        .lock()
        .map_err(|error| eyre::eyre!(error.to_string()))? = Stages::default();
    let started = Instant::now();
    let outcome = engine.evaluate_rustred_symbolic(&settings, input.as_view());
    let evaluation_seconds = started.elapsed().as_secs_f64();
    let owners_after = engine.rustred_application_owner_preparation_statistics()?;
    let owner_seconds = owners_after
        .preparation_time
        .saturating_sub(owners_before.preparation_time)
        .as_secs_f64();
    let stages = LOGGER
        .0
        .lock()
        .map_err(|error| eyre::eyre!(error.to_string()))?;
    let saved_seconds = stages
        .seconds
        .get("saved_rule_application")
        .map_or(0.0, |samples| samples.iter().sum::<f64>());
    let statistics = engine
        .cache_statistics(5)?
        .saved_application
        .map(|s| s.to_document());
    let mut receipt = json!({
        "schema":"vakint.five-loop-cache-measurement.v1","identity":identity,
        "case_id":root.id,"variant":options.variant,"sample":options.sample,"input_sha256":format!("{:x}",Sha256::digest(input.to_canonical_string().as_bytes())),
        "declared_common_powers":root.powers,"rationale":root.rationale,
        "frozen_lowered_roots":root.lowered_roots,"frozen_protected_roots":root.protected_roots,
        "actual_lowered_roots":stages.input_roots,
        "provisional":options.provisional,
        "fresh_process":true,"data_loading_seconds":data_loading_seconds,"cache_loading_seconds":cache_loading_seconds,
        "evaluation_seconds":evaluation_seconds,"saved_application_including_owner_seconds":saved_seconds,
        "owner_preparation_seconds":owner_seconds,"ibp_excluding_owner_seconds":(saved_seconds-owner_seconds).max(0.0),
        "owner_preparation_during_load_seconds":owners_before.preparation_time.as_secs_f64(),
        "stages":stages.seconds,"lowering_receipts":stages.lowered,"warnings":stages.warnings,"cache_statistics":statistics,
    });
    let authentication = stages.authenticate(
        &root.lowered_roots,
        &root.protected_roots,
        options.provisional,
    );
    receipt["native_keys_authenticated"] = authentication
        .as_ref()
        .is_ok_and(|authenticated| *authenticated)
        .into();
    receipt["native_key_authentication"] = if !stages.input_roots_seen {
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
        Ok(result) => {
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
            receipt["value_sha256"] = digest(&options.work_dir.join("value.symbolica"))?.into();
            receipt["guards_sha256"] = digest(&options.work_dir.join("guards.json"))?.into();
        }
        Err(error) => {
            receipt["success"] = false.into();
            receipt["error"] = error.to_string().into();
        }
    }
    if let Err(error) = authentication {
        receipt["success"] = false.into();
        receipt["error"] = error.to_string().into();
    }
    ensure!(
        &options.identity()? == identity,
        "measurement inputs changed during evaluation"
    );
    write_json(&options.work_dir.join("worker.json"), &receipt)?;
    ensure!(
        receipt["success"] == true,
        "public evaluation failed; inspect measurement.json"
    );
    Ok(())
}

fn main() -> Result<()> {
    let Some(options) = Options::parse()? else {
        return Ok(());
    };
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
    let roots = if options.suite == "tuning" {
        &workload.tuning
    } else {
        &workload.held_out
    };
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
    let variants: &[&str] = if options.cache.is_some() {
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
                        let deadline = active["deadline_unix_seconds"]
                            .as_f64()
                            .ok_or_else(|| eyre::eyre!("invalid active worker deadline"))?;
                        ensure!(
                            deadline.is_finite() && now >= deadline,
                            "an interrupted supervisor's worker may still be active; preserve its deadline"
                        );
                        json!({"success":false,"interrupted":true,"deadline_exhausted":true,
                            "error":"interrupted attempt retained without retry"})
                    } else {
                        let started = Instant::now();
                        let deadline = now + options.max_case_seconds as f64;
                        write_json(
                            &active_path,
                            &json!({
                                "identity":identity,"case_id":root.id,"variant":variant,"sample":sample,
                                "deadline_unix_seconds":deadline,
                            }),
                        )?;
                        let log = File::create(directory.join("worker.log"))?;
                        let launched = Command::new(std::env::current_exe()?)
                            .args(std::env::args_os().skip(1))
                            .arg("--worker")
                            .arg("--case-id")
                            .arg(&root.id)
                            .arg("--variant")
                            .arg(variant)
                            .arg("--sample")
                            .arg(sample.to_string())
                            .arg("--worker-deadline")
                            .arg(deadline.to_string())
                            .arg("--work-dir")
                            .arg(&directory)
                            .stdin(Stdio::null())
                            .stdout(log.try_clone()?)
                            .stderr(log)
                            .spawn();
                        match launched {
                            Ok(mut child) => {
                                let (status, expired) = loop {
                                    let expired = started.elapsed()
                                        >= Duration::from_secs(options.max_case_seconds)
                                        || SystemTime::now()
                                            .duration_since(UNIX_EPOCH)?
                                            .as_secs_f64()
                                            >= deadline;
                                    if let Some(status) = child.try_wait()? {
                                        break (status, expired);
                                    }
                                    if expired {
                                        child.kill()?;
                                        break (child.wait()?, true);
                                    }
                                    std::thread::sleep(Duration::from_millis(200));
                                };
                                json!({"success":status.success() && !expired,
                                    "exit_code":status.code(),"exit_status":status.to_string(),
                                    "deadline_exhausted":expired || status.code()==Some(124),
                                    "elapsed_seconds":started.elapsed().as_secs_f64()})
                            }
                            Err(error) => json!({"success":false,"launch_error":error.to_string(),
                                "elapsed_seconds":started.elapsed().as_secs_f64()}),
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
                        receipt["supervisor_error"] = "worker did not exit successfully within its deadline; retained without retry; inspect supervision and worker.log".into();
                        if receipt["error"].is_null() {
                            receipt["error"] = receipt["supervisor_error"].clone();
                        }
                    }
                    receipt["supervision"] = supervision;
                    write_json(&receipt_path, &receipt)?;
                }
                let receipt: Value = serde_json::from_slice(&fs::read(&receipt_path)?)?;
                validate_receipt(&receipt, &identity, &root.id, variant, sample)?;
                if receipt["success"] == true {
                    ensure!(
                        receipt["supervision"]["success"] == true,
                        "successful receipt lacks successful supervision"
                    );
                    ensure!(
                        receipt["value_sha256"] == digest(&directory.join("value.symbolica"))?
                            && receipt["guards_sha256"] == digest(&directory.join("guards.json"))?,
                        "saved measurement value or guards differ from their receipt"
                    );
                }
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
        let exact = successful
            && selected.iter().all(|receipt| {
                receipt["value_sha256"] == selected[0]["value_sha256"]
                    && receipt["guards_sha256"] == selected[0]["guards_sha256"]
            });
        let mut row =
            json!({"id":root.id,"successful":successful,"exact_output_and_guard_parity":exact});
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
        row["backing_reuse"] = if root_hits > 0 {
            "incidental-requested-root-hit"
        } else if intermediate_hits > 0 {
            "intermediate-only"
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
        && cases
            .iter()
            .all(|case| case["exact_output_and_guard_parity"] == true);
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
        &json!({"schema":"vakint.five-loop-cache-suite.v1","identity":identity,"complete":complete,"cases":cases,"measurements":measurements,"geometric_mean_uncached_over_cached_ratio":geometric_mean,"provisional":options.provisional,"native_keys_authenticated":authenticated,"execution_provenance_complete":provenance_complete,"final_held_out_comparison":authenticated && provenance_complete && options.suite=="held-out" && options.cache.is_some() && options.repetitions>=3,"selection_allowed":authenticated && provenance_complete && options.suite=="tuning" && options.cache.is_some() && options.repetitions>=3}),
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
    fn public_input_requires_the_complete_frozen_lowered_union() {
        let first = vec![1; 15];
        let mut second = first.clone();
        second[0] = 2;
        let mut stages = Stages::default();
        assert!(
            stages
                .authenticate(&[first.clone()], &[first.clone()], false)
                .is_err()
        );
        assert!(!stages.authenticate(&[], &[], true).unwrap());
        stages
            .capture_input_roots(
                &json!({"phase":"input_roots","roots":[first.clone()]}).to_string(),
            )
            .unwrap();
        assert!(
            stages
                .authenticate(&[second.clone()], &[first.clone()], false)
                .is_err()
        );
        assert!(
            stages
                .authenticate(&[second.clone()], &[first.clone()], true)
                .is_err()
        );
        stages
            .capture_input_roots(
                &json!({"phase":"input_roots","roots":[second.clone(), first.clone()]}).to_string(),
            )
            .unwrap();
        assert!(
            stages
                .authenticate(&[first.clone(), second], &[first], false)
                .unwrap()
        );
    }
}
