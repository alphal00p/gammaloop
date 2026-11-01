//! Prepare and check the Ordered step caches shipped in vakint-external-data.
//! Build with `--release --features feynkit-ingress`.
//!
//! The shipped workload is deterministic: the zero-momentum and quadratic QCD
//! ladder rungs at one to four loops (the inputs of
//! `tests/ordered_reduction_rungs.rs`), the four-loop benchmark integrals of
//! `tests/fixtures/benchmark_reduction.json`, and the five-loop benchmarks
//! `5L/benchmarks/{zero-momentum,quadratic}` of the data. Each case runs
//! through the same public calls as the rung tests, cold (shipped caches
//! ignored), recording into one engine cache directory; the caches depend only
//! on the recorded steps, so the export is reproducible byte for byte.
//!
//! ```text
//! ordered_cache_maintenance prepare --artifacts DATA --work WORK --output EXPORT
//!     [--loops 1,2,3,4,5] [--five-loop zero,quadratic]
//!     [--vakint-revision REV] [--rustred-revision REV]
//!   cold workload into WORK/cache-root, then export the caches of the programs
//!   loaded from DATA to EXPORT/<L>L/ordered-cache{/,.json}; stage EXPORT into
//!   vakint-external-data with its tooling (see its README).
//! ordered_cache_maintenance check --artifacts DATA [--loops 1,2,3,4,5]
//!   for every program loaded from DATA: its cache directory name, rule tables
//!   and stamp under this build, and whether the cache DATA ships for it is
//!   accepted (exit status 1 if any shipped cache is rejected).
//! ```
//! `VAKINT_ORDERED_WEIGHT_ARITHMETIC` selects the weight arithmetic (records do
//! not depend on it) and `VAKINT_ORDERED_CACHE_MAX_RESIDENT_BYTES` a resident
//! ceiling per program during `prepare` (records that do not fit are not
//! kept: a way to prepare a trimmed cache of the records computed first).
use ::vakint;
#[path = "../tests/fixtures/qcd_ladder/mod.rs"]
#[allow(dead_code)]
mod fixture;

use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};
use symbolica::atom::AtomCore;
use vakint::{
    EvaluationMethod, EvaluationOrder, LoopNormalizationFactor, OfficialExternalDataOptions,
    OrderedWeightArithmetic, RustRedEvaluationOptions, RustRedMasterMode, RustRedMonitorOptions,
    ShippedOrderedCacheStatus, Vakint, VakintSettings, vakint_parse,
};

/// The four-loop single-integral benchmark (prepared-cache roots and held-out
/// powers of the former depth-first gate), evaluated after the rungs.
const FOUR_LOOP_INTEGRALS: &str = include_str!("../tests/fixtures/benchmark_reduction.json");

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

struct Options {
    artifacts: PathBuf,
    work: Option<PathBuf>,
    output: Option<PathBuf>,
    loops: Vec<usize>,
    five_loop: Vec<String>,
    vakint_revision: String,
    rustred_revision: String,
}

fn parse() -> fixture::Result<(String, Options)> {
    let mut arguments = std::env::args().skip(1);
    let command = arguments
        .next()
        .ok_or("usage: ordered_cache_maintenance prepare|check [options]")?;
    let mut options = Options {
        artifacts: PathBuf::new(),
        work: None,
        output: None,
        loops: vec![1, 2, 3, 4, 5],
        five_loop: vec!["zero".into()],
        vakint_revision: "unknown".into(),
        rustred_revision: "unknown".into(),
    };
    while let Some(argument) = arguments.next() {
        let mut value = || {
            arguments
                .next()
                .ok_or_else(|| format!("{argument} needs a value"))
        };
        match argument.as_str() {
            "--artifacts" => options.artifacts = value()?.into(),
            "--work" => options.work = Some(value()?.into()),
            "--output" => options.output = Some(value()?.into()),
            "--loops" => {
                options.loops = value()?
                    .split(',')
                    .map(str::parse)
                    .collect::<Result<_, _>>()?;
                if options.loops.iter().any(|loops| !(1..=5).contains(loops)) {
                    return Err("--loops takes loop counts 1 through 5".into());
                }
            }
            "--five-loop" => {
                options.five_loop = value()?.split(',').map(str::to_owned).collect();
                if options
                    .five_loop
                    .iter()
                    .any(|projection| !matches!(projection.as_str(), "zero" | "quadratic"))
                {
                    return Err("--five-loop takes zero and/or quadratic".into());
                }
            }
            "--vakint-revision" => options.vakint_revision = value()?,
            "--rustred-revision" => options.rustred_revision = value()?,
            other => return Err(format!("unknown argument {other}").into()),
        }
    }
    if options.artifacts.as_os_str().is_empty() {
        return Err("--artifacts DATA (an external-data directory) is required".into());
    }
    Ok((command, options))
}

/// One workload case: its prepared input and description.
struct Case {
    name: String,
    prepared: fixture::LadderFixture,
    row: serde_json::Value,
}

fn workload(options: &Options, inputs: &Path) -> fixture::Result<Vec<Case>> {
    let mut cases = Vec::new();
    for &loops in &options.loops {
        if loops == 5 {
            for projection in &options.five_loop {
                let (kind, label, directory) = match projection.as_str() {
                    "zero" => (
                        fixture::Projection::ZeroMomentumMetric,
                        "zero-momentum",
                        "5L/benchmarks/zero-momentum",
                    ),
                    _ => (
                        fixture::Projection::QuadraticTransverse,
                        "quadratic",
                        "5L/benchmarks/quadratic",
                    ),
                };
                let path = options.artifacts.join(directory);
                let prepared = fixture::LadderFixture::load_prepared_at_loop_count(&path, kind, 5)?;
                let manifest = std::fs::read(path.join("prepared-input.json"))?;
                let name = format!("5L-{label}");
                cases.push(Case {
                    row: serde_json::json!({
                        "name": name, "kind": "five-loop-benchmark", "loop_count": 5,
                        "projection": label, "prepared_input": format!("{directory}/prepared-input.json"),
                        "prepared_input_sha256": sha256(&manifest),
                        "input_sha256": sha256(prepared.input.to_canonical_string().as_bytes()),
                    }),
                    name,
                    prepared,
                });
            }
            continue;
        }
        for (projection, label) in [
            (fixture::Projection::ZeroMomentumMetric, "zero-momentum"),
            (fixture::Projection::QuadraticTransverse, "quadratic"),
        ] {
            let name = format!("{loops}L-{label}");
            let prepared = fixture::comparison::prepared(&inputs.join(&name), projection, loops)?;
            cases.push(Case {
                row: serde_json::json!({
                    "name": name, "kind": "qcd-ladder-rung", "loop_count": loops,
                    "projection": label,
                    "source": "crates/vakint/tests/fixtures/qcd_ladder (tests/ordered_reduction_rungs.rs)",
                    "input_sha256": sha256(prepared.input.to_canonical_string().as_bytes()),
                }),
                name,
                prepared,
            });
        }
    }
    Ok(cases)
}

/// Run the workload cold into `work/cache-root`; returns the root and the
/// workload rows with the retained results' digests.
fn run_cold(options: &Options, work: &Path) -> fixture::Result<(PathBuf, Vec<serde_json::Value>)> {
    let root = work.join("cache-root");
    if root.exists() {
        return Err(format!("{} exists; use a fresh --work directory", root.display()).into());
    }
    std::fs::create_dir_all(&root)?;
    // The fixture reads its cache settings from the environment.
    // SAFETY: single-threaded before any evaluation starts.
    unsafe {
        std::env::set_var("VAKINT_ORDERED_CACHE", "read-write");
        std::env::set_var("VAKINT_ORDERED_CACHE_DIR", &root);
        std::env::set_var("VAKINT_ORDERED_CACHE_SHIPPED", "off");
    }
    let arithmetic = fixture::comparison::arithmetic_from_environment()?;
    let mut rows = Vec::new();
    for case in workload(options, &work.join("inputs"))? {
        let started = std::time::Instant::now();
        let completed = fixture::comparison::run(
            &case.prepared,
            &work.join("runs").join(&case.name),
            &fixture::comparison::RunOptions {
                arithmetic,
                artifacts: Some(options.artifacts.clone()),
                monitor: RustRedMonitorOptions::default(),
                numerical_reference: false,
            },
        )?;
        let cache = &completed.report["ordered_step_cache"];
        eprintln!(
            "{}: {:.1} s, retained {}, cache records {} inserts {} hits {}",
            case.name,
            started.elapsed().as_secs_f64(),
            completed.report["retained"]["sha256"],
            cache["records"],
            cache["inserts"],
            cache["hits"]
        );
        let mut row = case.row;
        row["retained_sha256"] = completed.report["retained"]["sha256"].clone();
        if let Some(guards) = completed.report["guards"]["payload"]["sha256"].as_str() {
            row["nonzero_conditions_sha256"] = guards.into();
        }
        rows.push(row);
        if case.name == "4L-quadratic" {
            rows.extend(four_loop_integrals(options, arithmetic)?);
        }
    }
    Ok((root, rows))
}

fn retained_settings(arithmetic: OrderedWeightArithmetic) -> fixture::Result<VakintSettings> {
    Ok(VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        use_dot_product_notation: true,
        form_exe_path: "/ordered-cache-maintenance-must-not-use-form".into(),
        number_of_terms_in_epsilon_expansion: 5,
        ordered_cache: fixture::comparison::ordered_cache_from_environment()?,
        n_cores: fixture::comparison::n_cores_from_environment()?,
        ordered_weight_arithmetic: arithmetic,
        ..Default::default()
    })
}

/// The four-loop benchmark integrals, one engine, retained masters.
fn four_loop_integrals(
    options: &Options,
    arithmetic: OrderedWeightArithmetic,
) -> fixture::Result<Vec<serde_json::Value>> {
    let engine = engine(options, 4)?;
    let settings = retained_settings(arithmetic)?;
    let cases: Vec<serde_json::Value> = serde_json::from_str(FOUR_LOOP_INTEGRALS)?;
    let mut rows = Vec::new();
    for case in cases {
        let input = case["input"]
            .as_str()
            .ok_or("benchmark case without input")?;
        let started = std::time::Instant::now();
        let value = engine.evaluate(&settings, vakint_parse!(input)?.as_view())?;
        let retained = sha256(value.to_canonical_string().as_bytes());
        let name = format!(
            "4L-integral-{}-p{}",
            case["family"].as_str().unwrap_or("?"),
            case["power"]
        );
        eprintln!(
            "{name}: {:.1} s, retained {retained}",
            started.elapsed().as_secs_f64()
        );
        rows.push(serde_json::json!({
            "name": name, "kind": "integral", "loop_count": 4,
            "source": format!("crates/vakint/tests/fixtures/benchmark_reduction.json case {}", case["case"]),
            "input_sha256": sha256(input.as_bytes()), "retained_sha256": retained,
        }));
    }
    Ok(rows)
}

/// An engine with the data of `options` loaded up to `loops` (and the
/// five-loop saved application at five loops), so that it lists every
/// program the shipped caches can belong to.
fn engine(options: &Options, loops: usize) -> fixture::Result<Vakint> {
    let mut engine = Vakint::new()?;
    engine.load_official_external_data(
        loops,
        &OfficialExternalDataOptions {
            local_artifacts_path: Some(options.artifacts.clone()),
            ..Default::default()
        },
    )?;
    if loops == 5 {
        let mut settings = VakintSettings::default();
        engine =
            engine.with_rustred_application(&mut settings, options.artifacts.join("5L/ibps"))?;
    }
    Ok(engine)
}

fn main() -> fixture::Result<()> {
    let (command, options) = parse()?;
    let maximum = *options.loops.iter().max().ok_or("--loops is empty")?;
    match command.as_str() {
        "prepare" => {
            let work = options.work.clone().ok_or("prepare needs --work")?;
            let output = options.output.clone().ok_or("prepare needs --output")?;
            let (root, rows) = run_cold(&options, &work)?;
            let producer = vakint::ordered_cache_producer(
                &options.vakint_revision,
                &options.rustred_revision,
                &format!(
                    "ordered_cache_maintenance prepare --loops {} --five-loop {} \
                     (VAKINT_ORDERED_WEIGHT_ARITHMETIC={}, \
                     VAKINT_ORDERED_CACHE_MAX_RESIDENT_BYTES={})",
                    options
                        .loops
                        .iter()
                        .map(usize::to_string)
                        .collect::<Vec<_>>()
                        .join(","),
                    options.five_loop.join(","),
                    std::env::var("VAKINT_ORDERED_WEIGHT_ARITHMETIC").unwrap_or("default".into()),
                    std::env::var("VAKINT_ORDERED_CACHE_MAX_RESIDENT_BYTES")
                        .unwrap_or("default".into()),
                ),
            );
            let report =
                engine(&options, maximum)?.export_ordered_caches(&root, &output, producer, rows)?;
            println!("{}", serde_json::to_string_pretty(&report)?);
        }
        "check" => {
            let checks = engine(&options, maximum)?.shipped_ordered_caches()?;
            println!("{}", serde_json::to_string_pretty(&checks)?);
            if checks.iter().any(|check| {
                check
                    .shipped
                    .as_ref()
                    .is_some_and(|shipped| shipped.status == ShippedOrderedCacheStatus::Rejected)
            }) {
                return Err("some shipped Ordered caches are rejected by this build".into());
            }
        }
        other => return Err(format!("unknown command {other}").into()),
    }
    Ok(())
}
