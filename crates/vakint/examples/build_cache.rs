//! Reproducible, explicitly selected five-loop cache preparation.
//!
//! Run `cargo run --release -p vakint --example build_cache -- --help`.
//! A supervisor and an independent worker watchdog enforce the wall deadline.
//! Interrupted checkpoints are retained, but never called complete.
use std::{
    fs::{self, File, OpenOptions},
    io::{Read, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::atomic::AtomicBool,
    time::{Duration, Instant, SystemTime, UNIX_EPOCH},
};

use eyre::{Context, Result, bail, ensure};
use rustred::family::IntegralKey;
use rustred_app::CachePreparationStatus;
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use symbolica::utils::resources::{ResourceBudget, ResourceLimits, checkpoint};
use vakint::{RustRedCachePreparationOptions, RustRedCacheRoot, RustRedCacheWorkload, Vakint};

const MIB: u64 = 1024 * 1024;

#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "Shared proof validation and observation serve measure_cache; builder uses identity and loading"
    )
)]
#[path = "build_cache/owner_proof.rs"]
mod owner_proof;
#[path = "build_cache/workload.rs"]
mod workload;
use owner_proof::OwnerProof;
use workload::Workload;

struct Options {
    artifacts: PathBuf,
    workload: PathBuf,
    source_census: PathBuf,
    output: PathBuf,
    work_dir: PathBuf,
    phase: String,
    prepare_seconds: u64,
    max_seconds: u64,
    retained_coefficient_mib: Option<u64>,
    max_process_rss_mib: Option<u64>,
    max_elapsed_seconds: Option<u64>,
    worker: bool,
    worker_budget: Option<Attempt>,
}

impl Options {
    fn validate_budget(phase: &str, prepare_seconds: u64, max_seconds: u64) -> Result<()> {
        ensure!(matches!(phase, "pilot" | "final"), "invalid phase");
        ensure!(
            phase != "pilot" || max_seconds <= 10800,
            "pilot time is limited to three cumulative hours"
        );
        ensure!(
            prepare_seconds > 0 && prepare_seconds < max_seconds && max_seconds <= 36000,
            "require 0 < preparation budget < hard budget <= 36000 seconds"
        );
        ensure!(
            phase != "final" || (prepare_seconds <= 32400 && max_seconds - prepare_seconds >= 3600),
            "reserve at least one hour of the ten-hour final budget for finalization"
        );
        Ok(())
    }

    fn identity(&self, workload_hash: &str, workload: &Workload) -> Result<Value> {
        let source = Path::new(env!("CARGO_MANIFEST_DIR"));
        let repository = source
            .parent()
            .and_then(Path::parent)
            .ok_or_else(|| eyre::eyre!("cannot locate workspace source provenance"))?;
        let revision = Command::new("git")
            .args(["rev-parse", "HEAD"])
            .current_dir(repository)
            .output()?;
        let changes = Command::new("git")
            .args(["diff", "--binary", "HEAD"])
            .current_dir(repository)
            .output()?;
        let untracked = Command::new("git")
            .args(["ls-files", "--others", "--exclude-standard", "-z"])
            .current_dir(repository)
            .output()?;
        ensure!(
            revision.status.success() && changes.status.success() && untracked.status.success(),
            "cannot record runtime source provenance"
        );
        let mut untracked_hashes = serde_json::Map::new();
        for path in String::from_utf8(untracked.stdout)?.split('\0') {
            let extension = Path::new(path).extension().and_then(|value| value.to_str());
            if matches!(
                extension,
                Some("rs" | "py" | "toml" | "lock" | "nix" | "typ")
            ) && fs::symlink_metadata(repository.join(path))?.is_file()
            {
                untracked_hashes.insert(path.to_owned(), digest(&repository.join(path))?.into());
            }
        }
        let compiler = option_env!("VAKINT_CACHE_BUILD_RUSTC");
        let profile = option_env!("VAKINT_CACHE_BUILD_PROFILE");
        let built_revision = option_env!("VAKINT_CACHE_BUILD_GIT_REVISION");
        ensure!(
            self.phase != "final"
                || [compiler, profile, built_revision]
                    .iter()
                    .all(|value| value.is_some_and(|value| !value.trim().is_empty())),
            "final builds require compile-time VAKINT_CACHE_BUILD_RUSTC, VAKINT_CACHE_BUILD_PROFILE and VAKINT_CACHE_BUILD_GIT_REVISION provenance"
        );
        #[cfg(target_os = "linux")]
        let affinity = fs::read_to_string("/proc/self/status")?
            .lines()
            .find_map(|line| line.strip_prefix("Cpus_allowed_list:").map(str::trim))
            .ok_or_else(|| eyre::eyre!("missing CPU affinity provenance"))?
            .to_owned();
        #[cfg(not(target_os = "linux"))]
        let affinity = "unavailable on this platform";
        let worker_environment = [
            "RAYON_NUM_THREADS",
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
        ]
        .map(|name| (name, std::env::var(name).ok()))
        .into_iter()
        .collect::<std::collections::BTreeMap<_, _>>();
        Ok(json!({
            "schema": "vakint.five-loop-cache-build-identity.v3",
            "workload_sha256": workload_hash,
            "source_census": self.source_census.canonicalize()?,
            "source_census_sha256": workload.source_census_sha256,
            "data_manifest_sha256": workload.data_manifest_sha256,
            "output": self.output, "phase": self.phase,
            "prepare_seconds": self.prepare_seconds, "max_seconds": self.max_seconds,
            "retained_coefficient_mib": self.retained_coefficient_mib,
            "resources": {
                "max_process_rss_mib": self.max_process_rss_mib,
                "max_elapsed_seconds": self.max_elapsed_seconds,
                "hard_cumulative_seconds": self.max_seconds,
            },
            "owner_proof": OwnerProof::pinned_identity(workload, &self.workload, &self.artifacts)?,
            "build": {
                "executable_sha256": digest(&std::env::current_exe()?)?,
                "driver_sha256": format!("{:x}", Sha256::digest(include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/examples/build_cache.rs")))),
                "cargo_lock_sha256": format!("{:x}", Sha256::digest(include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../Cargo.lock")))),
                "workload_schema_sha256": format!("{:x}", Sha256::digest(include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/examples/build_cache/workload.rs")))),
                "owner_proof_validator_sha256": format!("{:x}", Sha256::digest(include_bytes!(concat!(env!("CARGO_MANIFEST_DIR"), "/examples/build_cache/owner_proof.rs")))),
                "rustc": compiler, "profile": profile, "git_revision": built_revision,
                "debug_assertions": cfg!(debug_assertions),
                "target_arch": std::env::consts::ARCH, "target_os": std::env::consts::OS,
            },
            // Runtime files are recorded separately: their contents are not
            // evidence of which source produced an existing executable.
            "runtime_source": {
                "git_revision": String::from_utf8(revision.stdout)?.trim(),
                "tracked_diff_sha256": format!("{:x}", Sha256::digest(changes.stdout)),
                "untracked_files": untracked_hashes,
                "cargo_lock_sha256": digest(&repository.join("Cargo.lock"))?,
            },
            "execution": {
                "cpu_affinity": affinity,
                "available_parallelism": std::thread::available_parallelism()?.get(),
                "worker_environment": worker_environment,
            },
        }))
    }

    fn parse() -> Result<Option<Self>> {
        let mut artifacts = None;
        let mut workload = None;
        let mut source_census = None;
        let mut output = None;
        let mut work_dir = None;
        let mut phase = "pilot".to_owned();
        let mut prepare_seconds = 840;
        let mut max_seconds = 900;
        let mut retained_coefficient_mib: Option<u64> = None;
        let mut max_process_rss_mib: Option<u64> = None;
        let mut max_elapsed_seconds: Option<u64> = None;
        let mut worker = false;
        let mut worker_budget = None;
        let mut args = std::env::args().skip(1);
        while let Some(arg) = args.next() {
            if matches!(arg.as_str(), "--help" | "-h") {
                println!(
                    "build_cache --artifacts PATH --workload JSON --source-census JSON --output FILE.gz --work-dir PATH\n\
                    --phase pilot|final       Label the measured run (default: pilot)\n\
                    --prepare-seconds N       Cooperative native preparation budget (default: 840)\n\
                    --max-seconds N           Hard process deadline, including finalization (default: 900)\n\
                    --retained-coefficient-mib N Native coefficient capacity (default: unbounded; alias: --ram-mib)\n\
                    --max-process-rss-mib N   Whole worker RSS budget (default: unbounded)\n\
                    --max-elapsed-seconds N   Tighter cumulative worker time budget (default: hard deadline)\n\
                    Final production uses --phase final --prepare-seconds 32400 --max-seconds 36000.\n\
                    Freeze training/tuning/held-out inputs first. This command does not select roots.\n\
                    Run in an optimized build with fixed affinity; all work and receipts stay in --work-dir.\n\
                    Record build metadata using VAKINT_CACHE_BUILD_RUSTC, VAKINT_CACHE_BUILD_PROFILE,\n\
                    and VAKINT_CACHE_BUILD_GIT_REVISION when compiling the driver."
                );
                return Ok(None);
            }
            if arg == "--worker" {
                worker = true;
                continue;
            }
            let value = args
                .next()
                .ok_or_else(|| eyre::eyre!("missing value for {arg}"))?;
            match arg.as_str() {
                "--artifacts" => artifacts = Some(value.into()),
                "--workload" => workload = Some(value.into()),
                "--source-census" => source_census = Some(value.into()),
                "--output" => output = Some(value.into()),
                "--work-dir" => work_dir = Some(value.into()),
                "--phase" => phase = value,
                "--prepare-seconds" => prepare_seconds = value.parse()?,
                "--max-seconds" => max_seconds = value.parse()?,
                "--retained-coefficient-mib" | "--ram-mib" => {
                    retained_coefficient_mib = Some(value.parse()?);
                }
                "--max-process-rss-mib" => max_process_rss_mib = Some(value.parse()?),
                "--max-elapsed-seconds" => max_elapsed_seconds = Some(value.parse()?),
                "--worker-budget" => worker_budget = Some(serde_json::from_str(&value)?),
                _ => bail!("unknown argument {arg}"),
            }
        }
        Self::validate_budget(&phase, prepare_seconds, max_seconds)?;
        for mib in [retained_coefficient_mib, max_process_rss_mib]
            .into_iter()
            .flatten()
        {
            ensure!(
                mib > 0 && mib.checked_mul(MIB).is_some(),
                "invalid memory allowance"
            );
        }
        ensure!(
            max_elapsed_seconds != Some(0),
            "elapsed allowance must be positive"
        );
        Ok(Some(Self {
            artifacts: artifacts.ok_or_else(|| eyre::eyre!("--artifacts is required"))?,
            workload: workload.ok_or_else(|| eyre::eyre!("--workload is required"))?,
            source_census: source_census
                .ok_or_else(|| eyre::eyre!("--source-census is required"))?,
            output: output.ok_or_else(|| eyre::eyre!("--output is required"))?,
            work_dir: work_dir.ok_or_else(|| eyre::eyre!("--work-dir is required"))?,
            phase,
            prepare_seconds,
            max_seconds,
            retained_coefficient_mib,
            max_process_rss_mib,
            max_elapsed_seconds,
            worker,
            worker_budget,
        }))
    }
}

/// All earlier elapsed time is charged against preparation, conservatively
/// retaining the finalization reserve after an ordinary resumed attempt.
#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Attempt {
    previous_elapsed_seconds: f64,
    started_unix_seconds: f64,
    prepare_seconds: f64,
    hard_seconds: f64,
    hard_deadline_unix_seconds: f64,
}

impl Attempt {
    fn new(previous: f64, started: f64, prepare_seconds: u64, max_seconds: u64) -> Result<Self> {
        ensure!(
            previous.is_finite() && previous >= 0.0 && started.is_finite() && started >= 0.0,
            "invalid previous elapsed time or attempt start"
        );
        let prepare_seconds = prepare_seconds as f64 - previous;
        let hard_seconds = max_seconds as f64 - previous;
        ensure!(
            prepare_seconds >= 1.0 && hard_seconds > prepare_seconds,
            "the cumulative preparation budget is exhausted; preserve this incomplete candidate"
        );
        Ok(Self {
            previous_elapsed_seconds: previous,
            started_unix_seconds: started,
            prepare_seconds,
            hard_seconds,
            hard_deadline_unix_seconds: started + hard_seconds,
        })
    }

    fn validate(&self, prepare_seconds: u64, max_seconds: u64) -> Result<()> {
        let expected = Self::new(
            self.previous_elapsed_seconds,
            self.started_unix_seconds,
            prepare_seconds,
            max_seconds,
        )?;
        // The command line and durable JSON round-trip floating-point seconds.
        // Permit their sub-microsecond rounding without extending a budget.
        ensure!(
            (self.prepare_seconds - expected.prepare_seconds).abs() < 1e-6
                && (self.hard_seconds - expected.hard_seconds).abs() < 1e-6
                && (self.hard_deadline_unix_seconds - expected.hard_deadline_unix_seconds).abs()
                    < 1e-6,
            "worker budget differs from the cumulative build limits"
        );
        Ok(())
    }

    fn preparation_allowance(&self, now: f64, worker_elapsed: Duration) -> Result<u64> {
        let elapsed = (now - self.started_unix_seconds)
            .max(worker_elapsed.as_secs_f64())
            .max(0.0);
        let remaining = (self.prepare_seconds - elapsed).floor();
        ensure!(
            remaining >= 1.0,
            "loading exhausted the remaining cumulative preparation allowance"
        );
        Ok(remaining as u64)
    }

    fn resource_allowance(
        &self,
        now: f64,
        worker_elapsed: Duration,
        max_elapsed_seconds: Option<u64>,
    ) -> Result<Duration> {
        let allowed = max_elapsed_seconds
            .map(|limit| limit as f64 - self.previous_elapsed_seconds)
            .unwrap_or(self.hard_seconds)
            .min(self.hard_seconds);
        let elapsed = (now - self.started_unix_seconds)
            .max(worker_elapsed.as_secs_f64())
            .max(0.0);
        let remaining = allowed - elapsed;
        ensure!(
            remaining > 0.0,
            "the cumulative worker resource time budget is exhausted"
        );
        Ok(Duration::try_from_secs_f64(remaining)?)
    }

    fn arm_watchdog(&self) -> Result<()> {
        let seconds = (self.hard_deadline_unix_seconds - unix_seconds()?).min(self.hard_seconds);
        ensure!(seconds > 0.0, "the worker hard deadline has expired");
        let allowance = Duration::from_secs_f64(seconds);
        let deadline = self.hard_deadline_unix_seconds;
        std::thread::Builder::new()
            .name("cache-deadline".into())
            .spawn(move || {
                let started = Instant::now();
                loop {
                    let remaining = allowance.saturating_sub(started.elapsed());
                    if remaining.is_zero() || unix_seconds().map_or(true, |now| now >= deadline) {
                        eprintln!("cache worker exhausted its hard wall-clock deadline");
                        std::process::exit(124);
                    }
                    std::thread::sleep(remaining.min(Duration::from_secs(1)));
                }
            })
            .context("start independent worker deadline watchdog")?;
        Ok(())
    }

    fn interrupted_elapsed(&self, now: f64) -> Result<f64> {
        ensure!(
            now >= self.hard_deadline_unix_seconds,
            "an interrupted worker may still be active until UNIX time {}; wait for its watchdog deadline",
            self.hard_deadline_unix_seconds
        );
        // A dead supervisor cannot prove when its worker stopped. Charge the
        // complete allowance rather than granting an unaccounted retry.
        Ok(self.previous_elapsed_seconds + self.hard_seconds)
    }
}

fn unix_seconds() -> Result<f64> {
    Ok(SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs_f64())
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

fn write_json(path: &Path, document: &Value) -> Result<()> {
    let temporary = path.with_extension("json.pending");
    let mut file = File::create(&temporary)?;
    file.write_all(&serde_json::to_vec_pretty(document)?)?;
    file.sync_all()?;
    fs::rename(temporary, path)?;
    File::open(path.parent().unwrap_or(Path::new(".")))?.sync_all()?;
    Ok(())
}

fn build_cache(
    options: &Options,
    workload: &Workload,
    workload_hash: &str,
    attempt: &Attempt,
    started: Instant,
) -> Result<()> {
    let identity: Value =
        serde_json::from_slice(&fs::read(options.work_dir.join("identity.json"))?)?;
    ensure!(
        identity == options.identity(workload_hash, workload)?,
        "workload, census, source or execution settings changed after supervisor preflight"
    );
    checkpoint()?;
    let mut engine = Vakint::new()?;
    engine.load_external_data(&options.artifacts, 5)?;
    let binding = engine
        .cache_statistics(5)?
        .saved_application
        .ok_or_else(|| eyre::eyre!("missing loaded application binding"))?
        .to_document()["binding"]
        .clone();
    let proof = OwnerProof::load(workload, &options.workload, &options.artifacts, &binding)?;
    ensure!(
        options.phase != "final" || proof.is_some(),
        "final preparation requires validated owner proof"
    );
    let loading_seconds = started.elapsed().as_secs_f64();
    let mut preparation = RustRedCachePreparationOptions {
        checkpoint_path: Some(options.work_dir.join("checkpoints")),
        ..Default::default()
    };
    preparation.storage.time_limit_seconds =
        attempt.preparation_allowance(unix_seconds()?, started.elapsed())?;
    preparation.storage.ram_bytes = options.retained_coefficient_mib.map(|mib| mib * MIB);
    // The final gzip target is checked below; scratch checkpoints are not
    // constrained by the final representation size.
    preparation.storage.disk_bytes = None;
    preparation.storage.max_roots = workload.training.len();
    preparation.workload = RustRedCacheWorkload::Explicit(
        workload
            .training
            .iter()
            .map(|root| {
                Ok(RustRedCacheRoot {
                    namespace: 0,
                    integral: IntegralKey::try_new(root.powers.clone())?,
                })
            })
            .collect::<Result<Vec<_>>>()?,
    );
    let mut progress = OpenOptions::new()
        .create(true)
        .append(true)
        .open(options.work_dir.join("progress.jsonl"))?;
    let report = engine.prepare_cache(
        5,
        &options.output,
        &preparation,
        &AtomicBool::new(false),
        |report| {
            serde_json::to_writer(&mut progress, &report.to_document())
                .and_then(|()| writeln!(progress).map_err(serde_json::Error::io))
                .map_err(|error| {
                    rustred_app::AppError::new(
                        rustred_app::AppErrorKind::Serialization,
                        error.to_string(),
                    )
                })?;
            progress.flush().map_err(|error| {
                rustred_app::AppError::new(
                    rustred_app::AppErrorKind::Serialization,
                    error.to_string(),
                )
            })
        },
    )?;
    ensure!(
        digest(&options.workload)? == workload_hash,
        "workload changed during preparation; preserve the candidate and investigate"
    );
    ensure!(
        digest(&options.artifacts.join("manifest.json"))? == workload.data_manifest_sha256,
        "external data changed during preparation; preserve the candidate and investigate"
    );
    ensure!(
        digest(&options.source_census)? == workload.source_census_sha256,
        "source census changed during preparation; preserve the candidate and investigate"
    );
    let complete = report.status == CachePreparationStatus::Complete;
    let mut receipt = json!({
        "schema": "vakint.five-loop-cache-build.v2", "phase": options.phase,
        "workload_sha256": workload_hash, "data_manifest_sha256": workload.data_manifest_sha256,
        "source_census_sha256": workload.source_census_sha256,
        "loading_seconds": loading_seconds, "elapsed_seconds": started.elapsed().as_secs_f64(),
        "complete": complete, "preparation": report.to_document(),
        "selected_training_roots": workload.training.len(), "frozen_held_out_inputs": workload.held_out.len(),
        "source_census": options.source_census, "attempt": attempt,
        "publishable": false,
        "owner_proof":proof.as_ref().map(OwnerProof::identity),
        "owner_root_correlations":proof.as_ref().map(OwnerProof::correlations),
    });
    if complete {
        let bytes = fs::metadata(&options.output)?.len();
        receipt["compressed_bytes"] = bytes.into();
        receipt["cache_sha256"] = digest(&options.output)?.into();
        receipt["within_50_mib_target"] = (bytes <= 50 * MIB).into();
        // A larger candidate requires a separately reviewed measured exception;
        // this baseline builder never silently relaxes the selected size target.
    }
    // Publication is decided only after the supervisor verifies cumulative
    // elapsed time and the worker's successful exit.
    write_json(&options.work_dir.join("worker.json"), &receipt)?;
    ensure!(
        complete,
        "cache preparation was incomplete; inspect its durable progress and report"
    );
    ensure!(
        receipt["within_50_mib_target"] == true,
        "completed cache exceeds 50 MiB; retain it as a candidate, not a publishable selection"
    );
    Ok(())
}

fn main() -> Result<()> {
    let Some(options) = Options::parse()? else {
        return Ok(());
    };
    let worker_started = Instant::now();
    if options.worker {
        let attempt = options
            .worker_budget
            .as_ref()
            .ok_or_else(|| eyre::eyre!("worker requires supervisor budget"))?;
        attempt.validate(options.prepare_seconds, options.max_seconds)?;
        attempt.arm_watchdog()?;
        let limits = ResourceLimits {
            max_process_rss_bytes: options.max_process_rss_mib.map(|mib| mib * MIB),
            max_elapsed: Some(attempt.resource_allowance(
                unix_seconds()?,
                worker_started.elapsed(),
                options.max_elapsed_seconds,
            )?),
        };
        // One inherited scope spans loading, native preparation and compaction.
        // Native publication checks its budget before commit; no blanket
        // post-return check can turn its committed report into a ResourceError.
        return ResourceBudget::new(limits).run(|| {
            checkpoint()?;
            run(options, worker_started)
        });
    }
    run(options, worker_started)
}

fn run(options: Options, worker_started: Instant) -> Result<()> {
    if options.worker {
        let attempt = options
            .worker_budget
            .as_ref()
            .ok_or_else(|| eyre::eyre!("worker requires supervisor budget"))?;
        let active: Value =
            serde_json::from_slice(&fs::read(options.work_dir.join("supervisor.json"))?)?;
        ensure!(
            active["state"] == "active" && active["attempt"] == serde_json::to_value(attempt)?,
            "worker has no matching durable active attempt"
        );
    } else {
        ensure!(
            options.worker_budget.is_none(),
            "worker budget is reserved for supervisor invocations"
        );
    }
    let workload_bytes = fs::read(&options.workload)?;
    let workload_hash = format!("{:x}", Sha256::digest(&workload_bytes));
    let workload: Workload = serde_json::from_slice(&workload_bytes)?;
    if options.phase == "final" {
        workload.validate_final_selection()?;
    } else {
        workload.validate()?;
    }
    checkpoint()?;
    ensure!(
        digest(&options.artifacts.join("manifest.json"))? == workload.data_manifest_sha256,
        "external data manifest differs from the frozen workload"
    );
    ensure!(
        digest(&options.source_census)? == workload.source_census_sha256,
        "source census differs from the frozen workload"
    );
    if options.worker {
        return build_cache(
            &options,
            &workload,
            &workload_hash,
            options.worker_budget.as_ref().unwrap(),
            worker_started,
        );
    }
    fs::create_dir_all(&options.work_dir)?;
    let supervisor_lock = OpenOptions::new()
        .read(true)
        .write(true)
        .create(true)
        .truncate(false)
        .open(options.work_dir.join("supervisor.lock"))?;
    supervisor_lock.try_lock().map_err(|error| {
        eyre::eyre!("cannot lock cache work directory; another supervisor may be active: {error}")
    })?;
    ensure!(
        !options.output.exists(),
        "preserve the existing output; choose a new candidate file"
    );
    let identity = options.identity(&workload_hash, &workload)?;
    let identity_path = options.work_dir.join("identity.json");
    if identity_path.exists() {
        ensure!(
            serde_json::from_slice::<Value>(&fs::read(&identity_path)?)? == identity,
            "work directory belongs to a different workload, census, build, budget or execution setup"
        );
    } else {
        write_json(&identity_path, &identity)?;
    }
    let supervisor_path = options.work_dir.join("supervisor.json");
    let previous_seconds = if supervisor_path.exists() {
        let previous: Value = serde_json::from_slice(&fs::read(&supervisor_path)?)?;
        if previous["state"] == "active" {
            let interrupted: Attempt = serde_json::from_value(previous["attempt"].clone())?;
            interrupted.validate(options.prepare_seconds, options.max_seconds)?;
            let elapsed = interrupted.interrupted_elapsed(unix_seconds()?)?;
            write_json(
                &supervisor_path,
                &json!({
                    "state": "interrupted", "success": false,
                    "cumulative_elapsed_seconds": elapsed,
                    "attempt": interrupted,
                }),
            )?;
            elapsed
        } else {
            ensure!(
                matches!(previous["state"].as_str(), Some("finished" | "interrupted")),
                "unknown supervisor receipt; do not reset an unaccounted build"
            );
            previous["cumulative_elapsed_seconds"]
                .as_f64()
                .ok_or_else(|| {
                    eyre::eyre!("previous supervisor receipt has no cumulative elapsed time")
                })?
        }
    } else {
        0.0
    };
    let log = OpenOptions::new()
        .create(true)
        .append(true)
        .open(options.work_dir.join("worker.log"))?;
    let started = Instant::now();
    let attempt = Attempt::new(
        previous_seconds,
        unix_seconds()?,
        options.prepare_seconds,
        options.max_seconds,
    )?;
    write_json(
        &supervisor_path,
        &json!({
            "state": "active", "attempt": attempt,
            "cumulative_elapsed_seconds": previous_seconds,
            "build_identity": identity,
        }),
    )?;
    let mut child = Command::new(std::env::current_exe()?)
        .args(std::env::args_os().skip(1))
        .arg("--worker")
        .arg("--worker-budget")
        .arg(serde_json::to_string(&attempt)?)
        .stdin(Stdio::null())
        .stdout(log.try_clone()?)
        .stderr(log)
        .spawn()
        .context("start cache worker")?;
    let (success, hard_deadline) = loop {
        let expired = started.elapsed() >= Duration::from_secs_f64(attempt.hard_seconds)
            || unix_seconds()? >= attempt.hard_deadline_unix_seconds;
        if let Some(status) = child.try_wait()? {
            break (status.success() && !expired, expired);
        }
        if expired {
            child.kill()?;
            child.wait()?;
            break (false, true);
        }
        std::thread::sleep(Duration::from_millis(200));
    };
    let elapsed_seconds = started.elapsed().as_secs_f64();
    let supervision = json!({
        "state": "finished", "success": success, "hard_deadline_exhausted": hard_deadline,
        "elapsed_seconds": elapsed_seconds, "cumulative_elapsed_seconds": previous_seconds + elapsed_seconds,
        "hard_budget_seconds": options.max_seconds,
        "executable_sha256": identity["build"]["executable_sha256"],
        "workload_sha256": workload_hash, "attempt": attempt,
        "build_identity": identity,
        "command_line": std::env::args().skip(1).collect::<Vec<_>>(),
    });
    write_json(&supervisor_path, &supervision)?;
    let mut attempts = OpenOptions::new()
        .create(true)
        .append(true)
        .open(options.work_dir.join("attempts.jsonl"))?;
    serde_json::to_writer(&mut attempts, &supervision)?;
    writeln!(attempts)?;
    attempts.sync_all()?;
    let worker_receipt = options.work_dir.join("worker.json");
    let mut publishable = false;
    if worker_receipt.exists() {
        let mut receipt: Value = serde_json::from_slice(&fs::read(&worker_receipt)?)?;
        ensure!(
            receipt["attempt"] == serde_json::to_value(&attempt)?,
            "worker receipt belongs to another attempt; preserve the candidate"
        );
        receipt["last_attempt_elapsed_seconds"] = receipt["elapsed_seconds"].clone();
        let cumulative_seconds = previous_seconds + started.elapsed().as_secs_f64();
        publishable = success
            && cumulative_seconds <= options.max_seconds as f64
            && receipt["complete"] == true
            && receipt["within_50_mib_target"] == true;
        receipt["elapsed_seconds"] = cumulative_seconds.into();
        receipt["executable_sha256"] = supervision["executable_sha256"].clone();
        receipt["build_identity"] = identity;
        receipt["supervision"] = supervision;
        receipt["publishable"] = publishable.into();
        write_json(&options.work_dir.join("build.json"), &receipt)?;
    }
    ensure!(
        publishable,
        "cache worker did not complete; see {}",
        options.work_dir.join("worker.log").display()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resumed_preparation_charges_loading_and_preserves_finalization_reserve() {
        Options::validate_budget("final", 32400, 36000).unwrap();
        assert!(Options::validate_budget("final", 32400, 35000).is_err());
        assert!(Options::validate_budget("final", 1, 3600).is_err());
        assert!(Options::validate_budget("final", 32401, 36000).is_err());
        let resumed = Attempt::new(7200.0, 1000.0, 32400, 36000).unwrap();
        assert_eq!(resumed.hard_seconds - resumed.prepare_seconds, 3600.0);
        assert_eq!(
            resumed
                .preparation_allowance(1100.75, Duration::from_secs(100))
                .unwrap(),
            25099,
        );
        // A backwards wall-clock step cannot restore loading time measured by
        // the worker's monotonic clock.
        assert_eq!(
            resumed
                .preparation_allowance(900.0, Duration::from_secs(100))
                .unwrap(),
            25100,
        );
        assert!(
            resumed
                .preparation_allowance(26200.0, Duration::ZERO)
                .is_err()
        );
    }

    #[test]
    fn one_outer_resource_clock_charges_previous_attempts_and_loading() {
        let attempt = Attempt::new(100.0, 1000.0, 840, 900).unwrap();
        assert_eq!(
            attempt
                .resource_allowance(1010.0, Duration::from_secs(5), None)
                .unwrap(),
            Duration::from_secs(790),
        );
        assert_eq!(
            attempt
                .resource_allowance(900.0, Duration::from_secs(10), Some(600))
                .unwrap(),
            Duration::from_secs(490),
        );
        assert!(
            attempt
                .resource_allowance(1500.0, Duration::ZERO, Some(600))
                .is_err()
        );
        assert!(
            attempt
                .resource_allowance(1000.0, Duration::ZERO, Some(100))
                .is_err()
        );
    }

    #[test]
    fn an_interrupted_supervisor_never_grants_an_unaccounted_retry() {
        let attempt = Attempt::new(7200.0, 1000.0, 32400, 36000).unwrap();
        assert!(attempt.interrupted_elapsed(29799.0).is_err());
        let charged = attempt.interrupted_elapsed(29800.0).unwrap();
        assert_eq!(charged, 36000.0);
        assert!(Attempt::new(charged, 30000.0, 32400, 36000).is_err());
        assert!(Attempt::new(f64::NAN, 1000.0, 32400, 36000).is_err());
        assert!(Attempt::new(-1.0, 1000.0, 32400, 36000).is_err());
    }

    #[test]
    fn worker_watchdog_enforces_deadline_without_a_supervisor() {
        const CHILD: &str = "VAKINT_CACHE_WATCHDOG_TEST_CHILD";
        if std::env::var_os(CHILD).is_some() {
            Attempt::new(0.0, unix_seconds().unwrap(), 1, 2)
                .unwrap()
                .arm_watchdog()
                .unwrap();
            std::thread::sleep(Duration::from_secs(3));
            panic!("worker outlived its watchdog deadline");
        }
        let status = Command::new(std::env::current_exe().unwrap())
            .args([
                "--exact",
                "tests::worker_watchdog_enforces_deadline_without_a_supervisor",
            ])
            .env(CHILD, "1")
            .status()
            .unwrap();
        assert_eq!(status.code(), Some(124));
    }
}
