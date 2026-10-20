//! Refresh docs/benchmark/data/results.json with independent first-call timings.
//!
//! Run with `cargo run --release -p vakint --example benchmark_reduction -- --help`. For one CPU,
//! prefix Cargo (or the compiled executable) with `taskset -c CPU` on Linux.
use std::{
    fs::{self, File},
    io::Read,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    time::Instant,
};

use eyre::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use symbolica::{
    atom::{Atom, AtomCore},
    id::Replacement,
};
use vakint::{
    EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor,
    OfficialExternalDataOptions, RustRedEvaluationOptions, RustRedMasterMode, Vakint,
    VakintSettings, vakint_parse, vakint_symbol,
};

const CASES: &str = include_str!("../docs/benchmark/cases.json");
const SCOPE: &str = "First public evaluate_integral call in a fresh process; retained masters, \
    common symbolic mUVsq and FMFTandMATAD normalization. Excludes engine/data loading, cache \
    attachment, input parsing and validation. FMFT includes version check, process launch, \
    reduction, I/O and output parsing. Cache chunks are cold in-process; OS pages are not flushed.";

#[derive(Deserialize)]
struct Case {
    case: String,
    family: String,
    power: i64,
    powers: Vec<i64>,
    kind: String,
    input: String,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Backend {
    Fmft,
    Cached,
}

impl Backend {
    fn name(self) -> &'static str {
        match self {
            Self::Fmft => "fmft",
            Self::Cached => "cached",
        }
    }
}

struct Options {
    output: PathBuf,
    artifacts: Option<PathBuf>,
    version: String,
    form: PathBuf,
    repetitions: usize,
    work_dir: PathBuf,
    child: Option<usize>,
    backend: Backend,
    repetition: usize,
}

impl Options {
    fn parse() -> Result<Option<Self>> {
        let mut result = Self {
            output: Path::new(env!("CARGO_MANIFEST_DIR")).join("docs/benchmark/data/results.json"),
            artifacts: None,
            version: "latest".into(),
            form: "form".into(),
            repetitions: 3,
            work_dir: std::env::current_dir()?.join("workspace/benchmark"),
            child: None,
            backend: Backend::Fmft,
            repetition: 0,
        };
        let mut args = std::env::args().skip(1);
        while let Some(arg) = args.next() {
            if arg == "--help" || arg == "-h" {
                println!(
                    "cargo run --release --locked -p vakint --example benchmark_reduction -- [OPTIONS]\n\
                     --output PATH       JSON destination (default: crate docs/benchmark/data/results.json)\n\
                     --artifacts PATH    Data checkout or installation base (default: current directory)\n\
                     --version VERSION   Official data identity (default: latest, resolved once)\n\
                     --form PATH         FORM executable (default: form on PATH)\n\
                     --repetitions N     Fresh processes per case/backend (default: 3)\n\
                     --work-dir PATH     Scratch/log directory (default: cwd/workspace/benchmark)\n\
                     Use an optimized Cargo profile; optionally pin the entire command to one CPU."
                );
                return Ok(None);
            }
            let value = args
                .next()
                .ok_or_else(|| eyre::eyre!("missing value for {arg}"))?;
            match arg.as_str() {
                "--output" => result.output = value.into(),
                "--artifacts" => result.artifacts = Some(value.into()),
                "--version" => result.version = value,
                "--form" => result.form = value.into(),
                "--repetitions" => result.repetitions = value.parse()?,
                "--work-dir" => result.work_dir = value.into(),
                "--child" => result.child = Some(value.parse()?),
                "--repetition" => result.repetition = value.parse()?,
                "--mode" => {
                    result.backend = match value.as_str() {
                        "fmft" => Backend::Fmft,
                        "cached" => Backend::Cached,
                        _ => bail!("unknown backend {value}"),
                    }
                }
                _ => bail!("unknown option {arg}; see --help"),
            }
        }
        ensure!(result.repetitions > 0, "repetitions must be positive");
        Ok(Some(result))
    }

    fn child(&self, cases: &[Case], index: usize) -> Result<()> {
        let case = cases
            .get(index)
            .ok_or_else(|| eyre::eyre!("invalid case index {index}"))?;
        let start = Instant::now();
        let mut engine = Vakint::new()?;
        if self.backend == Backend::Cached {
            let report = engine.load_official_external_data(
                4,
                &OfficialExternalDataOptions {
                    version: self.version.clone(),
                    local_artifacts_path: self.artifacts.clone(),
                    ..Default::default()
                },
            )?;
            ensure!(report.version == self.version, "data version changed");
            ensure!(engine.cache_statistics(4)?.attached, "4L cache is absent");
        }
        let loader_seconds = start.elapsed().as_secs_f64();
        let input = vakint_parse!(&case.input)?;
        let settings = VakintSettings {
            evaluation_order: match self.backend {
                Backend::Fmft => EvaluationOrder::fmft_only(Some(FMFTOptions {
                    expand_masters: false,
                    susbstitute_masters: false,
                })),
                Backend::Cached => {
                    EvaluationOrder(vec![EvaluationMethod::RustRed(RustRedEvaluationOptions {
                        master_mode: RustRedMasterMode::Retain,
                    })])
                }
            },
            integral_normalization_factor: LoopNormalizationFactor::FMFTandMATAD,
            form_exe_path: if self.backend == Backend::Fmft {
                self.form.to_string_lossy().into_owned()
            } else {
                "/cached-benchmark-must-not-use-form".into()
            },
            temporary_directory: Some(self.work_dir.to_string_lossy().into_owned()),
            use_dot_product_notation: true,
            ..Default::default()
        };
        let start = Instant::now();
        let value = engine.evaluate_integral(&settings, input.as_view())?;
        let reduction_seconds = start.elapsed().as_secs_f64();
        let stats = engine.cache_statistics(4)?;
        if self.backend == Backend::Cached {
            if case.kind == "base" {
                ensure!(stats.native.backing_hits > 0, "base root missed cache");
                ensure!(
                    stats.native.work.rule_applications() == 0,
                    "base root performed native reduction"
                );
            } else {
                ensure!(stats.native.backing_hits == 0, "held-out root was cached");
                ensure!(
                    stats.native.descendant_backing_hits > 0,
                    "held-out reduction did not reuse intermediate cache entries"
                );
            }
        }
        Json::write(
            &self.output,
            &json!({
                "case": case.case, "family": case.family, "power": case.power,
                "powers": case.powers, "kind": case.kind, "mode": self.backend.name(),
                "repetition": self.repetition, "input": input.to_canonical_string(),
                "output": value.to_canonical_string(), "loader_seconds": loader_seconds,
                "reduction_seconds": reduction_seconds,
                "prepared_entry_hits": stats.native.backing_hits,
                "prepared_descendant_hits": stats.native.descendant_backing_hits,
                "native_rule_applications": stats.native.work.rule_applications(),
            }),
        )
    }

    fn run(&self, cases: &[Case]) -> Result<()> {
        ensure!(
            !cfg!(debug_assertions),
            "use an optimized build: cargo run --release -p vakint --example benchmark_reduction"
        );
        let started = Instant::now();
        fs::create_dir_all(&self.work_dir)?;
        let scratch = tempfile::Builder::new()
            .prefix("reduction-")
            .tempdir_in(&self.work_dir)?
            .keep();
        eprintln!("Benchmark logs: {}", scratch.display());
        let form = Provenance::executable(&self.form)?;
        let executable = std::env::current_exe()?;
        let mut engine = Vakint::new()?;
        let report = engine.load_official_external_data(
            4,
            &OfficialExternalDataOptions {
                version: self.version.clone(),
                local_artifacts_path: self.artifacts.clone(),
                ..Default::default()
            },
        )?;
        ensure!(engine.cache_statistics(4)?.attached, "4L cache is absent");
        drop(engine);
        let root_manifest_path = report.path.join("manifest.json");
        let root_manifest_hash = Provenance::sha256(&root_manifest_path)?;
        let manifest_path = report.path.join("4L/manifest.json");
        let manifest_hash = Provenance::sha256(&manifest_path)?;
        let mut runs = Vec::new();
        for repetition in 0..self.repetitions {
            for (index, case) in cases.iter().enumerate() {
                let modes = if repetition % 2 == 0 {
                    [Backend::Fmft, Backend::Cached]
                } else {
                    [Backend::Cached, Backend::Fmft]
                };
                for mode in modes {
                    let label = format!("{index}-{}-{repetition}", mode.name());
                    let output = scratch.join(format!("{label}.json"));
                    let log = File::create(scratch.join(format!("{label}.log")))?;
                    let start = Instant::now();
                    let status = Command::new(&executable)
                        .arg("--child")
                        .arg(index.to_string())
                        .arg("--mode")
                        .arg(mode.name())
                        .arg("--repetition")
                        .arg(repetition.to_string())
                        .arg("--output")
                        .arg(&output)
                        .arg("--artifacts")
                        .arg(&report.path)
                        .arg("--version")
                        .arg(&report.version)
                        .arg("--form")
                        .arg(&form)
                        .arg("--work-dir")
                        .arg(&scratch)
                        .env("RAYON_NUM_THREADS", "1")
                        .env("OMP_NUM_THREADS", "1")
                        .env("MKL_NUM_THREADS", "1")
                        .stdout(Stdio::from(log.try_clone()?))
                        .stderr(Stdio::from(log))
                        .status()?;
                    ensure!(
                        status.success(),
                        "{label} failed ({status}); inspect {}",
                        scratch.display()
                    );
                    let mut record: Value = Json::read(&output)?;
                    record["process_seconds"] = json!(start.elapsed().as_secs_f64());
                    record["process_exit_code"] = json!(0);
                    eprintln!(
                        "{} {} {} repetition {}: {} s",
                        case.family,
                        case.kind,
                        mode.name(),
                        repetition + 1,
                        record["reduction_seconds"]
                    );
                    runs.push(record);
                    Json::write(&scratch.join("progress.json"), &runs)?;
                }
            }
        }
        ensure!(
            Provenance::sha256(&root_manifest_path)? == root_manifest_hash,
            "external release manifest changed during measurement"
        );
        ensure!(
            Provenance::sha256(&manifest_path)? == manifest_hash,
            "external data changed during measurement"
        );
        let parity = Validation::exact_parity(cases, &runs, &manifest_path)?;
        let summary = cases
            .iter()
            .map(|case| {
                let fmft = Measurements::summarize(case, "fmft", &runs)?;
                let cached = Measurements::summarize(case, "cached", &runs)?;
                Ok(json!({
                    "family": case.family, "power": case.power, "kind": case.kind,
                    "fmft_over_cached": fmft.median_seconds / cached.median_seconds,
                    "fmft": fmft, "cached": cached,
                }))
            })
            .collect::<Result<Vec<_>>>()?;
        Json::write(
            &self.output,
            &json!({
                "completed_utc": Provenance::command("date", &["-u", "+%Y-%m-%dT%H:%M:%SZ"])?,
                "elapsed_seconds": started.elapsed().as_secs_f64(),
                "validation_passed": true, "runs": runs.len(),
                "protocol": Provenance::collect(self, &form, &executable, &report.version, &root_manifest_hash, &manifest_hash)?,
                "summary": summary, "exact_parity": parity,
                "cached_loader_seconds": runs.iter().filter(|r| r["mode"] == "cached")
                    .map(|r| &r["loader_seconds"]).collect::<Vec<_>>(),
                "raw_runs": runs,
            }),
        )?;
        eprintln!("Validated results: {}", self.output.display());
        Ok(())
    }
}

struct Validation;

impl Validation {
    fn exact_parity(cases: &[Case], runs: &[Value], manifest_path: &Path) -> Result<Value> {
        #[derive(Deserialize)]
        struct Master {
            family_sha256: String,
            powers: Vec<i64>,
            exact_basis_expression: String,
        }
        #[derive(Deserialize)]
        struct Manifest {
            masters: Vec<Master>,
        }
        let manifest: Manifest = Json::read(manifest_path)?;
        let mut aliases = Vec::new();
        for master in manifest.masters {
            let powers = master
                .powers
                .iter()
                .map(ToString::to_string)
                .collect::<Vec<_>>()
                .join(",");
            let source = format!(
                "RustRedMaster(rustred_family_{},{powers})",
                master.family_sha256
            );
            let basis = vakint_parse!(&master.exact_basis_expression)?
                .replace(vakint_parse!("ep")?.to_pattern())
                .with(Atom::var(vakint_symbol!(
                    &VakintSettings::default().epsilon_symbol
                )));
            aliases.push(Replacement::new(vakint_parse!(&source)?, basis));
        }
        for case in cases {
            let selected = runs
                .iter()
                .filter(|r| r["case"] == case.case)
                .collect::<Vec<_>>();
            ensure!(!selected.is_empty(), "no results for {}", case.case);
            ensure!(
                selected.iter().all(|r| r["input"] == selected[0]["input"]),
                "input differs across measurements for {}",
                case.case
            );
            let get = |mode: &str| -> Result<Atom> {
                let records = selected
                    .iter()
                    .filter(|r| r["mode"] == mode)
                    .collect::<Vec<_>>();
                ensure!(!records.is_empty(), "missing {mode} result");
                ensure!(
                    records.iter().all(|r| r["output"] == records[0]["output"]),
                    "repetitions disagree for {} {mode}",
                    case.case
                );
                let output = records[0]["output"]
                    .as_str()
                    .ok_or_else(|| eyre::eyre!("missing output expression"))?;
                Ok(vakint_parse!(output)?
                    .replace(vakint_parse!("mUVsq")?.to_pattern())
                    .with(Atom::num(1))
                    .replace(vakint_parse!(&VakintSettings::default().mu_r_sq_symbol)?.to_pattern())
                    .with(Atom::num(1)))
            };
            let form = get("fmft")?;
            let native = get("cached")?.replace_multiple(&aliases);
            ensure!(
                !native.contains_symbol(vakint_symbol!("RustRedMaster")),
                "unmapped master for {}",
                case.case
            );
            ensure!(
                (native - form).expand().together().is_zero(),
                "exact backend mismatch for {}",
                case.case
            );
        }
        Ok(json!({
            "passed": cases.iter().map(|c| &c.case).collect::<Vec<_>>(),
            "method": "Exact rational equality at unit mass and renormalization scale in authenticated FMFT master basis; no epsilon truncation or floating tolerance",
        }))
    }
}

#[derive(Serialize)]
struct Measurements {
    median_seconds: f64,
    min_seconds: f64,
    max_seconds: f64,
    samples_seconds: Vec<f64>,
}

impl Measurements {
    fn summarize(case: &Case, mode: &str, runs: &[Value]) -> Result<Self> {
        let samples_seconds = runs
            .iter()
            .filter(|r| r["case"] == case.case && r["mode"] == mode)
            .map(|r| {
                r["reduction_seconds"]
                    .as_f64()
                    .filter(|v| v.is_finite() && *v > 0.0)
                    .ok_or_else(|| eyre::eyre!("invalid reduction duration"))
            })
            .collect::<Result<Vec<_>>>()?;
        ensure!(!samples_seconds.is_empty(), "no samples");
        let mut sorted = samples_seconds.clone();
        sorted.sort_by(f64::total_cmp);
        let median_seconds = (sorted[(sorted.len() - 1) / 2] + sorted[sorted.len() / 2]) / 2.0;
        Ok(Self {
            median_seconds,
            min_seconds: sorted[0],
            max_seconds: sorted[sorted.len() - 1],
            samples_seconds,
        })
    }
}

struct Provenance;

impl Provenance {
    fn command(program: &str, args: &[&str]) -> Result<String> {
        let output = Command::new(program).args(args).output()?;
        ensure!(output.status.success(), "{program} {args:?} failed");
        Ok(String::from_utf8(output.stdout)?.trim().to_owned())
    }

    fn executable(path: &Path) -> Result<PathBuf> {
        if path.is_file() {
            return Ok(path.canonicalize()?);
        }
        if path.components().count() == 1 {
            for directory in std::env::split_paths(&std::env::var_os("PATH").unwrap_or_default()) {
                let candidate = directory.join(path);
                if candidate.is_file() {
                    return Ok(candidate.canonicalize()?);
                }
            }
        }
        bail!("FORM executable not found: {}", path.display())
    }

    fn sha256(path: &Path) -> Result<String> {
        let mut file = File::open(path)?;
        let mut hasher = Sha256::new();
        let mut buffer = [0; 64 * 1024];
        loop {
            let count = file.read(&mut buffer)?;
            if count == 0 {
                break;
            }
            hasher.update(&buffer[..count]);
        }
        Ok(format!("{:x}", hasher.finalize()))
    }

    fn collect(
        options: &Options,
        form: &Path,
        executable: &Path,
        version: &str,
        root_manifest_hash: &str,
        manifest_hash: &str,
    ) -> Result<Value> {
        let cpuinfo = fs::read_to_string("/proc/cpuinfo").unwrap_or_default();
        let status = fs::read_to_string("/proc/self/status").unwrap_or_default();
        let crate_path = env!("CARGO_MANIFEST_DIR");
        Ok(json!({
            "source_commit": Self::command("git", &["-C", crate_path, "rev-parse", "HEAD"])?,
            "source_dirty": !Self::command("git", &["-C", crate_path, "status", "--porcelain"])?.is_empty(),
            "build_profile": "release (default documented invocation)",
            "configured_profile": {"opt_level": 3, "lto": "fat"},
            "debug_assertions": cfg!(debug_assertions),
            "binary_sha256": Self::sha256(executable)?,
            "bench_source_sha256": format!("{:x}", Sha256::digest(include_bytes!("benchmark_reduction.rs"))),
            "cases_sha256": format!("{:x}", Sha256::digest(CASES.as_bytes())),
            "form_sha256": Self::sha256(form)?,
            "form_version": Self::command(&form.to_string_lossy(), &["-v"])?,
            "rustc": Self::command("rustc", &["--version"])?,
            "platform": format!("{}-{}", std::env::consts::ARCH, std::env::consts::OS),
            "cpu_model": cpuinfo.lines().find_map(|line| line.strip_prefix("model name").and_then(|s| s.split_once(':')).map(|(_, v)| v.trim())),
            "cpu_affinity": status.lines().find_map(|line| line.strip_prefix("Cpus_allowed_list:")).map(str::trim),
            "data_version": version, "data_manifest_sha256": manifest_hash,
            "release_manifest_sha256": root_manifest_hash,
            "repetitions": options.repetitions,
            "mode_order": "alternating between repetitions; sequential, never overlapping",
            "scope": SCOPE,
        }))
    }
}

struct Json;

impl Json {
    fn read<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T> {
        serde_json::from_reader(File::open(path)?)
            .with_context(|| format!("reading {}", path.display()))
    }

    fn write(path: &Path, value: &impl Serialize) -> Result<()> {
        if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
            fs::create_dir_all(parent)?;
        }
        let mut text = serde_json::to_vec_pretty(value)?;
        text.push(b'\n');
        fs::write(path, text)?;
        Ok(())
    }
}

fn main() -> Result<()> {
    let Some(options) = Options::parse()? else {
        return Ok(());
    };
    let cases: Vec<Case> = serde_json::from_str(CASES)?;
    if let Some(index) = options.child {
        options.child(&cases, index)
    } else {
        options.run(&cases)
    }
}
