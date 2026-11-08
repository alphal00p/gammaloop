//! Prepare, check and replay the light 1L-3L p cache that vakint-external-data
//! ships in its `p_integral/` area (branch `p_integral`; M2-V6). Build with
//! `cargo build --release -p vakint --example p_cache_maintenance`.
//!
//! The cache is the complete Ordered step cache of Vakint's built-in 1L-3L p
//! library for one fixed workload (`tests/support/p_cache.rs`: the 1L-3L
//! benchmark grid, the 2L/3L physics cells, the Forcer batteries and the
//! anchors), recorded cold under each of the three weight arithmetics in turn
//! (exact first): a modular arithmetic visits a few steps exact does not (two
//! 2L benchmark cells, 9 steps), so the cache holds the union and serves all
//! three without a miss.
//!
//! ```text
//! p_cache_maintenance prepare --output WORK [--n-cores N]
//!     [--vakint-revision REV] [--rustred-revision REV]
//!   the workload cold (shipped caches ignored) under exact, modular_polynomial
//!   and modular_sampled into one read-write WORK/cache-root,
//!   then export_ordered_caches into WORK/export (p_integral/ordered-cache/,
//!   p_integral/ordered-cache.json); WORK must not exist. Stage it with
//!   `python3 tools/p_integral_area.py stage --export WORK/export` in
//!   vakint-external-data.
//! p_cache_maintenance check --artifacts DATA
//!   loads DATA (1L-3L) from its directory and the p library (one cache-off
//!   evaluation) and prints the p row of shipped_ordered_caches(); exit status
//!   0 iff that cache is used in place from DATA.
//! p_cache_maintenance replay --artifacts DATA --arithmetic LIST [--n-cores N]
//!     [--report FILE]
//!   per arithmetic (exact, modular_polynomial, modular_sampled): the cache-off
//!   value of every workload input, then the same inputs read-only in place
//!   from DATA (no cache directory). Fails unless the p cache is in place,
//!   every input has 0 misses, the workload has hits and every value equals
//!   its cache-off value. Prints (and writes) hits, misses and seconds.
//! ```
#[path = "../tests/support/p_cache.rs"]
mod p_cache;

use std::path::{Path, PathBuf};
use std::time::Instant;

use symbolica::atom::Atom;
use vakint::{
    EvaluationMethod, EvaluationOrder, OrderedCacheMode, OrderedCacheSettings,
    OrderedWeightArithmetic, RustRedEvaluationOptions, RustRedMasterMode,
    ShippedOrderedCacheStatus, Vakint, VakintSettings, vakint_parse,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

const P_LABEL: &str = "p-library-1-3";
/// The weight arithmetics prepare records and replay is run under.
const ARITHMETICS: [&str; 3] = ["exact", "modular_polynomial", "modular_sampled"];
const BUBBLE: &str = "topo(prop(1,edge(1,2),k(1),0,1)*prop(2,edge(2,1),k(1)-p(1),0,1))";

struct Options {
    artifacts: Option<PathBuf>,
    output: Option<PathBuf>,
    report: Option<PathBuf>,
    arithmetic: Vec<String>,
    n_cores: usize,
    vakint_revision: String,
    rustred_revision: String,
}

fn parse() -> Result<(String, Options)> {
    let mut arguments = std::env::args().skip(1);
    let command = arguments
        .next()
        .ok_or("usage: p_cache_maintenance prepare|check|replay [options]")?;
    let mut options = Options {
        artifacts: None,
        output: None,
        report: None,
        arithmetic: Vec::new(),
        n_cores: 4,
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
            "--artifacts" => options.artifacts = Some(value()?.into()),
            "--output" => options.output = Some(value()?.into()),
            "--report" => options.report = Some(value()?.into()),
            "--arithmetic" => {
                options.arithmetic = value()?.split(',').map(str::to_owned).collect();
            }
            "--n-cores" => options.n_cores = value()?.parse()?,
            "--vakint-revision" => options.vakint_revision = value()?,
            "--rustred-revision" => options.rustred_revision = value()?,
            other => return Err(format!("unknown argument {other}").into()),
        }
    }
    Ok((command, options))
}

/// RustRed alone, retained masters, `arithmetic`, `cache`.
fn settings(
    options: &Options,
    arithmetic: OrderedWeightArithmetic,
    ordered_cache: OrderedCacheSettings,
) -> VakintSettings {
    VakintSettings {
        evaluation_order: EvaluationOrder(vec![EvaluationMethod::RustRed(
            RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            },
        )]),
        form_exe_path: "/p-cache-maintenance-must-not-use-form".into(),
        use_dot_product_notation: true,
        n_cores: options.n_cores,
        ordered_weight_arithmetic: arithmetic,
        ordered_cache,
        ..Default::default()
    }
}

/// An engine with the 1L-3L data of `artifacts` loaded from its directory.
fn engine(artifacts: &Path) -> Result<Vakint> {
    let mut engine = Vakint::new()?;
    engine.load_external_data(artifacts, 3)?;
    Ok(engine)
}

fn p_counters(engine: &Vakint) -> (u64, u64, u64) {
    engine
        .ordered_cache_statistics()
        .programs
        .iter()
        .find(|program| program.program == P_LABEL)
        .map_or((0, 0, 0), |program| {
            (
                program.native.hits,
                program.native.misses,
                program.native.inserts,
            )
        })
}

fn prepare(options: &Options) -> Result<()> {
    let work = options
        .output
        .as_ref()
        .ok_or("prepare needs --output WORK")?;
    if work.exists() {
        return Err(format!("{} exists; use a fresh WORK directory", work.display()).into());
    }
    let root = work.join("cache-root");
    std::fs::create_dir_all(&root)?;
    let engine = Vakint::new()?;
    let workload = p_cache::p_cache_workload()?;
    let mut rows: Vec<serde_json::Value> = workload.iter().map(|input| input.row()).collect();
    for name in ARITHMETICS {
        let settings = settings(
            options,
            name.parse()?,
            OrderedCacheSettings {
                shipped: false,
                ..OrderedCacheSettings::read_write(&root)
            },
        );
        for (input, row) in workload.iter().zip(&mut rows) {
            let before = p_counters(&engine);
            let started = Instant::now();
            drop(engine.evaluate(&settings, input.input.as_view())?);
            let seconds = started.elapsed().as_secs_f64();
            let inserts = p_counters(&engine).2 - before.2;
            row["inserts"][name] = inserts.into();
            row["seconds"][name] = seconds.into();
            eprintln!("{name} {}: {seconds:.2} s, {inserts} inserts", input.name);
        }
    }
    let producer = {
        let mut producer = vakint::ordered_cache_producer(
            &options.vakint_revision,
            &options.rustred_revision,
            &format!(
                "p_cache_maintenance prepare --n-cores {} (exact, modular_polynomial, \
                 modular_sampled; built-in library)",
                options.n_cores
            ),
        );
        producer["p"] = serde_json::json!({
            "persistence": "rustred ordered step cache: ordered-step-cache.txt + RR-ORDERED-STEP1 \
                            shards (+ ordered-tail.rrtail)",
            "arithmetics": ARITHMETICS,
            "adjoint_tail": std::env::var("RUSTRED_ORDERED_ADJOINT_TAIL").unwrap_or("default".into()),
            "n_cores": options.n_cores,
        });
        producer
    };
    let report = engine.export_ordered_caches(&root, &work.join("export"), producer, rows)?;
    if report.p_integral.is_none() {
        return Err("the workload recorded no p cache".into());
    }
    println!("{}", serde_json::to_string_pretty(&report)?);
    Ok(())
}

fn check(options: &Options) -> Result<()> {
    let artifacts = options
        .artifacts
        .as_ref()
        .ok_or("check needs --artifacts")?;
    let engine = engine(artifacts)?;
    let load = settings(
        options,
        OrderedWeightArithmetic::default(),
        OrderedCacheSettings::off(),
    );
    drop(engine.evaluate(&load, vakint_parse!(BUBBLE)?.as_view())?);
    let checks = engine.shipped_ordered_caches()?;
    let row = checks
        .iter()
        .find(|check| check.program == P_LABEL)
        .ok_or("the p program is not listed")?;
    println!("{}", serde_json::to_string_pretty(row)?);
    let data = artifacts.canonicalize()?;
    match &row.shipped {
        Some(shipped)
            if shipped.status == ShippedOrderedCacheStatus::InPlace
                && shipped.source.starts_with(&data) =>
        {
            Ok(())
        }
        _ => Err(format!(
            "the p cache of {} is not used in place",
            artifacts.display()
        )
        .into()),
    }
}

fn replay(options: &Options) -> Result<()> {
    let artifacts = options
        .artifacts
        .as_ref()
        .ok_or("replay needs --artifacts")?;
    if options.arithmetic.is_empty() {
        return Err("replay needs --arithmetic LIST".into());
    }
    let workload = p_cache::p_cache_workload()?;
    let mut failures = Vec::new();
    let mut summary = Vec::new();
    for name in &options.arithmetic {
        let arithmetic: OrderedWeightArithmetic = name.parse()?;
        let cold_engine = engine(artifacts)?;
        let off = settings(options, arithmetic, OrderedCacheSettings::off());
        let started = Instant::now();
        let cold: Vec<Atom> = workload
            .iter()
            .map(|input| cold_engine.evaluate(&off, input.input.as_view()))
            .collect::<std::result::Result<_, _>>()?;
        let cold_seconds = started.elapsed().as_secs_f64();
        drop(cold_engine);

        let warm_engine = engine(artifacts)?;
        if !warm_engine.ordered_cache_statistics().programs.is_empty() {
            return Err("a fresh engine has attached caches".into());
        }
        let read_only = settings(
            options,
            arithmetic,
            OrderedCacheSettings {
                mode: OrderedCacheMode::ReadOnly,
                ..OrderedCacheSettings::default()
            },
        );
        let started = Instant::now();
        let (mut hits, mut misses, mut without_steps) = (0, 0, 0);
        for (input, cold) in workload.iter().zip(&cold) {
            let before = p_counters(&warm_engine);
            let value = warm_engine.evaluate(&read_only, input.input.as_view())?;
            let after = p_counters(&warm_engine);
            let (input_hits, input_misses) = (after.0 - before.0, after.1 - before.1);
            hits += input_hits;
            misses += input_misses;
            if input_hits + input_misses == 0 {
                without_steps += 1;
            }
            if input_misses != 0 {
                failures.push(format!("{name} {}: {input_misses} misses", input.name));
            }
            if value != *cold {
                failures.push(format!(
                    "{name} {}: value differs from cache-off",
                    input.name
                ));
            }
        }
        let warm_seconds = started.elapsed().as_secs_f64();
        let statistics = warm_engine.ordered_cache_statistics();
        let status = statistics
            .programs
            .iter()
            .find(|program| program.program == P_LABEL)
            .and_then(|program| program.shipped.as_ref())
            .map(|shipped| shipped.status.as_str());
        if status != Some("in-place") {
            failures.push(format!("{name}: the p cache is {status:?}, not in-place"));
        }
        if hits == 0 {
            failures.push(format!("{name}: no hits"));
        }
        let row = serde_json::json!({
            "arithmetic": name, "inputs": workload.len(), "hits": hits, "misses": misses,
            "inputs_without_steps": without_steps, "status": status,
            "cold_seconds": cold_seconds, "warm_seconds": warm_seconds,
        });
        eprintln!("{row}");
        summary.push(row);
    }
    let report = serde_json::json!({
        "verdict": if failures.is_empty() { "pass" } else { "fail" },
        "artifacts": artifacts, "n_cores": options.n_cores,
        "arithmetics": summary, "failures": failures,
    });
    let text = serde_json::to_string_pretty(&report)?;
    println!("{text}");
    if let Some(path) = &options.report {
        std::fs::write(path, text + "\n")?;
    }
    if failures.is_empty() {
        Ok(())
    } else {
        Err(format!("{} replay failures", failures.len()).into())
    }
}

fn main() -> Result<()> {
    let (command, options) = parse()?;
    match command.as_str() {
        "prepare" => prepare(&options),
        "check" => check(&options),
        "replay" => replay(&options),
        other => Err(format!("unknown command {other}").into()),
    }
}
