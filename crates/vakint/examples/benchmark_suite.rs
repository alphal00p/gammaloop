//! Measurement harness of `docs/benchmark.typ`: one measurement per fresh process.
//!
//! `docs/benchmark/make-benchmark.sh DATA_ROOT THIS_BINARY OUT_DIR` drives it
//! (repetitions, CPU pinning, comparisons and the report live in
//! `docs/benchmark/tools/run_benchmarks.py`). Build it with
//! `cargo build --release -p vakint --features feynkit-ingress,experimental-rustred --example benchmark_suite`.
//!
//! Subcommands (all write one JSON object to `--out`):
//! * `prepare-rung --loops L --projection zero-momentum|quadratic --work-dir DIR --out FILE.json`:
//!   the massless-quark QCD ladder (gluon self-energy with L-1 gluon rungs) of
//!   `tests/fixtures/qcd_ladder`, prepared natively (or reused from DIR); the
//!   evaluation input binds the native dimension to 4-2ε and expands the
//!   model couplings at G = 1 (identical input for every backend).
//! * `prepare-graph --data-root ROOT --loops L --master NAME [--powers p,...] --out FILE.json`:
//!   a master graph of the data root through native FeynKit ingress, with
//!   optional raised (dotted) common-basis powers.
//! * `measure --data-root ROOT --input FILE.json --backend rustred|alphaloop|matad|fmft
//!   --scope reduction|full [--cache cold|warm|default] [--form PATH] [--n-cores N]
//!   [--epsilon-terms K] [--precision DIGITS] --out FILE.json`.
//!   `reduction`: masters retained (RustRed Retain; FORM backends without
//!   master expansion/substitution). `full`: the complete public call
//!   (RustRed Numerical, FORM backends with their defaults), followed by an
//!   untimed numerical evaluation at unit mass and scale for the cross-backend
//!   comparison. Cache (RustRed): `cold` = no step cache, first call;
//!   `warm` = in-memory step cache, second identical call in the same process;
//!   `default` = engine defaults (shipped step caches on), first call.
//! * `identity --data-root ROOT --out FILE.json`: build and data identities.
use ::vakint;
#[path = "../tests/fixtures/qcd_ladder/mod.rs"]
#[allow(dead_code)]
mod fixture;
#[path = "external_data_acceptance/dot.rs"]
#[allow(dead_code)]
mod dot;

use std::{
    collections::BTreeMap,
    fs,
    io::Read,
    path::{Path, PathBuf},
    time::Instant,
};

use ahash::HashMap;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore, AtomView};
use vakint::{
    AlphaLoopOptions, EvaluationMethod, EvaluationOrder, FMFTOptions, LoopNormalizationFactor,
    MATADOptions, OrderedCacheMode, OrderedCacheSettings, OrderedWeightArithmetic,
    RustRedEvaluationOptions, RustRedMasterMode, Vakint, VakintSettings, vakint_parse,
    vakint_symbol,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

const SOURCE: &str = include_str!("benchmark_suite.rs");

struct Args(BTreeMap<String, String>);

impl Args {
    fn parse(mut raw: impl Iterator<Item = String>) -> Result<(String, Self)> {
        let command = raw.next().ok_or("missing subcommand; see the module docs")?;
        let mut map = BTreeMap::new();
        while let Some(key) = raw.next() {
            let key = key
                .strip_prefix("--")
                .ok_or_else(|| format!("unexpected argument {key}"))?
                .to_owned();
            let value = raw.next().ok_or_else(|| format!("--{key} needs a value"))?;
            map.insert(key, value);
        }
        Ok((command, Self(map)))
    }
    fn get(&self, key: &str) -> Result<&str> {
        self.0
            .get(key)
            .map(String::as_str)
            .ok_or_else(|| format!("missing --{key}").into())
    }
    fn opt(&self, key: &str) -> Option<&str> {
        self.0.get(key).map(String::as_str)
    }
    fn path(&self, key: &str) -> Result<PathBuf> {
        Ok(PathBuf::from(self.get(key)?))
    }
}

fn sha256_bytes(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn sha256_file(path: &Path) -> Result<String> {
    let mut file = fs::File::open(path)?;
    let mut digest = Sha256::new();
    let mut buffer = vec![0_u8; 1 << 16];
    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        digest.update(&buffer[..n]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

fn write_json(path: &Path, value: &Value) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let tmp = path.with_extension("json.tmp");
    fs::write(&tmp, serde_json::to_vec_pretty(value)?)?;
    fs::rename(tmp, path)?;
    Ok(())
}

fn peak_rss_bytes() -> Option<u64> {
    let status = fs::read_to_string("/proc/self/status").ok()?;
    status.lines().find_map(|line| {
        line.strip_prefix("VmHWM:")
            .and_then(|v| v.trim().strip_suffix("kB"))
            .and_then(|v| v.trim().parse::<u64>().ok())
            .map(|kb| kb * 1024)
    })
}

fn allowed_cpus() -> Option<String> {
    let status = fs::read_to_string("/proc/self/status").ok()?;
    status.lines().find_map(|line| {
        line.strip_prefix("Cpus_allowed_list:")
            .map(|v| v.trim().to_owned())
    })
}

fn parse_atom(text: &str) -> Result<Atom> {
    vakint_parse!(text).map_err(|e| format!("parse error: {e}").into())
}

fn epsilon(settings: &VakintSettings) -> Atom {
    Atom::var(vakint_symbol!(&settings.epsilon_symbol))
}

/// Count the distinct retained masters (RustRed or FORM-backend symbols).
fn master_census(value: AtomView) -> (usize, usize) {
    let mut rustred = std::collections::BTreeSet::new();
    let mut other = std::collections::BTreeSet::new();
    value.visitor(&mut |part| {
        if let AtomView::Fun(f) = part {
            let name = f.get_symbol().get_name().to_owned();
            if name.ends_with("RustRedMaster") {
                rustred.insert(part.to_canonical_string());
                return false;
            }
        }
        if let AtomView::Var(v) = part {
            let symbol = v.get_symbol();
            let name = symbol.get_name().to_owned();
            if name.contains("::PR") || name.contains("::mi") {
                other.insert(name);
            }
        }
        true
    });
    (rustred.len(), other.len())
}

// ------------------------------------------------------------- preparation --

fn prepare_rung(args: &Args) -> Result<()> {
    let loops: usize = args.get("loops")?.parse()?;
    let projection = match args.get("projection")? {
        "zero-momentum" => fixture::Projection::ZeroMomentumMetric,
        "quadratic" => fixture::Projection::QuadraticTransverse,
        other => return Err(format!("unknown projection {other}").into()),
    };
    let work = args.path("work-dir")?;
    let started = Instant::now();
    let prepared = fixture::comparison::prepared(&work, projection, loops)?;
    let preparation_seconds = started.elapsed().as_secs_f64();
    let settings = VakintSettings::default();
    // The fixture's native trace dimension, bound to 4-2ε exactly as the rung tests do.
    let native_d = rustred::algebra::CoefficientContext::try_new(["d"])?
        .parameter("d")
        .ok_or("missing native dimension")?
        .to_expression();
    let physical_d = Atom::num(4) - Atom::num(2) * epsilon(&settings);
    let model = feynkit_model::Model::from_json(include_str!("../tests/fixtures/qcd.json"))?;
    let input = prepared.input.replace(native_d).with(physical_d);
    let input = model
        .expand_couplings(&input)
        .replace(Atom::var(symbolica::symbol!("UFO::G")))
        .with(Atom::one())
        .expand();
    let text = input.to_canonical_string();
    write_json(
        &args.path("out")?,
        &json!({
            "schema": "vakint.benchmark-input.v1",
            "kind": "qcd-ladder-rung",
            "loop_count": loops,
            "projection": args.get("projection")?,
            "input": text,
            "input_sha256": sha256_bytes(text.as_bytes()),
            "fixture_scalar_terms": prepared.scalar_terms,
            "fixture_contributions": prepared.contributions,
            "fixture_preparation_seconds": prepared.preparation_seconds,
            "preparation_seconds": preparation_seconds,
            "substitutions": {"dimension": "rustred::{}::d -> 4-2ε", "couplings": "native QCD model couplings expanded, UFO::G = 1"},
            "numerical_point": ["vakint::{}::mUV", "vakint::{}::mUVsq"],
        }),
    )
}

fn prepare_graph(args: &Args) -> Result<()> {
    let root = args.path("data-root")?;
    let loops: usize = args.get("loops")?.parse()?;
    let pack = root.join(format!("{loops}L"));
    let manifest: Value = serde_json::from_str(&fs::read_to_string(pack.join("manifest.json"))?)?;
    let name = args.get("master")?;
    let master = manifest["masters"]
        .as_array()
        .ok_or("manifest without masters")?
        .iter()
        .find(|m| m["name"] == name)
        .ok_or_else(|| format!("no master {name}"))?
        .clone();
    let powers = args
        .opt("powers")
        .map(|text| {
            text.split(',')
                .map(|p| p.trim().parse::<i64>())
                .collect::<std::result::Result<Vec<_>, _>>()
        })
        .transpose()?;
    let input = dot::input_powers(&pack, &master, powers.as_deref())?;
    let text = input.to_canonical_string();
    write_json(
        &args.path("out")?,
        &json!({
            "schema": "vakint.benchmark-input.v1",
            "kind": "master-graph",
            "loop_count": loops,
            "master": name,
            "master_powers": master["powers"],
            "powers": powers,
            "input": text,
            "input_sha256": sha256_bytes(text.as_bytes()),
            "numerical_point": ["vakint::{}::mUV", "vakint::{}::mUVsq"],
        }),
    )
}

// ------------------------------------------------------------- measurement --

fn backend_order(backend: &str, scope: &str) -> Result<EvaluationOrder> {
    let full = scope == "full";
    Ok(match backend {
        "rustred" => EvaluationOrder(vec![EvaluationMethod::RustRed(RustRedEvaluationOptions {
            master_mode: if full {
                RustRedMasterMode::Numerical
            } else {
                RustRedMasterMode::Retain
            },
        })]),
        "alphaloop" => EvaluationOrder(vec![EvaluationMethod::AlphaLoop(AlphaLoopOptions {
            susbstitute_masters: full,
        })]),
        "matad" => EvaluationOrder::matad_only(Some(if full {
            MATADOptions::default()
        } else {
            MATADOptions {
                expand_masters: false,
                susbstitute_masters: false,
                substitute_hpls: false,
                direct_numerical_substition: false,
            }
        })),
        "fmft" => EvaluationOrder::fmft_only(Some(if full {
            FMFTOptions::default()
        } else {
            FMFTOptions {
                expand_masters: false,
                susbstitute_masters: false,
            }
        })),
        other => return Err(format!("unknown backend {other}").into()),
    })
}

fn measure(args: &Args) -> Result<()> {
    let process_start = Instant::now();
    let root = args.path("data-root")?;
    let input_file = args.path("input")?;
    let input_record: Value = serde_json::from_str(&fs::read_to_string(&input_file)?)?;
    let loops = input_record["loop_count"]
        .as_u64()
        .ok_or("input without loop_count")? as usize;
    let backend = args.get("backend")?.to_owned();
    let scope = args.get("scope")?.to_owned();
    let cache = args.opt("cache").unwrap_or("cold").to_owned();
    let n_cores: usize = args.opt("n-cores").unwrap_or("1").parse()?;
    let terms: i64 = match args.opt("epsilon-terms") {
        Some(t) => t.parse()?,
        None => loops as i64 + 1,
    };
    let precision: u32 = args.opt("precision").unwrap_or("32").parse()?;
    let form = args.opt("form").unwrap_or("form").to_owned();
    let work = args.path("out")?.with_extension("work");
    fs::create_dir_all(&work)?;

    let mut settings = VakintSettings {
        evaluation_order: backend_order(&backend, &scope)?,
        integral_normalization_factor: LoopNormalizationFactor::MSbar,
        number_of_terms_in_epsilon_expansion: terms,
        use_dot_product_notation: true,
        form_exe_path: if backend == "rustred" {
            "/benchmark-rustred-must-not-use-form".into()
        } else {
            form.clone()
        },
        temporary_directory: Some(work.to_string_lossy().into_owned()),
        run_time_decimal_precision: precision,
        ..Default::default()
    };
    settings.n_cores = n_cores;
    if let Some(arithmetic) = args.opt("arithmetic") {
        settings.ordered_weight_arithmetic = arithmetic.parse::<OrderedWeightArithmetic>()?;
    }
    settings.ordered_cache = match cache.as_str() {
        "cold" => OrderedCacheSettings {
            mode: OrderedCacheMode::Off,
            shipped: false,
            ..Default::default()
        },
        "warm" => OrderedCacheSettings {
            shipped: false,
            ..Default::default()
        },
        "default" => OrderedCacheSettings::default(),
        other => return Err(format!("unknown cache mode {other}").into()),
    };

    // Engine and data (excluded from the call time).
    let loading = Instant::now();
    let mut engine = Vakint::new()?;
    if backend == "rustred" {
        engine.load_external_data(&root, loops)?;
    }
    let loader_seconds = loading.elapsed().as_secs_f64();
    let parsing = Instant::now();
    let input = parse_atom(
        input_record["input"]
            .as_str()
            .ok_or("input record without input")?,
    )?;
    let parse_seconds = parsing.elapsed().as_secs_f64();

    let call = |engine: &Vakint| -> Result<(Atom, f64, Value)> {
        let started = Instant::now();
        #[cfg(feature = "experimental-rustred")]
        let (result, math) = {
            let capture = vakint::rustred_evaluation::experimental::capture_reduction_timing(|| {
                engine.evaluate(&settings, input.as_view())
            });
            let t = capture.timing;
            (
                capture.result,
                json!({
                    "reduction_math_seconds": t.reduction_math.map(|d| d.as_secs_f64()),
                    "gross_seconds": t.gross.as_secs_f64(),
                    "owner_preparation_seconds": t.measured_owner_preparation.as_secs_f64(),
                    "complete": t.complete,
                    "intervals": t.interval_count,
                    "incompleteness_reason": t.incompleteness_reason,
                }),
            )
        };
        #[cfg(not(feature = "experimental-rustred"))]
        let (result, math) = (engine.evaluate(&settings, input.as_view()), Value::Null);
        let seconds = started.elapsed().as_secs_f64();
        Ok((result?, seconds, math))
    };

    let (first, first_seconds, first_math) = call(&engine)?;
    let mut timed = (first.clone(), first_seconds, first_math.clone());
    let mut first_record = Value::Null;
    if cache == "warm" {
        first_record = json!({"seconds": first_seconds, "capture": first_math});
        let second = call(&engine)?;
        if second.0 != first {
            return Err("the warm call returned a different expression".into());
        }
        timed = second;
    }
    let (value, call_seconds, math) = timed;
    let rss_after_call = peak_rss_bytes();
    let canonical = value.to_canonical_string();
    let (rustred_masters, form_masters) = master_census(value.as_view());

    // Numerical evaluation (untimed for the comparison) of complete results.
    let mut numerical = Value::Null;
    if scope == "full" {
        let started = Instant::now();
        let mut evaluated = value.clone();
        for name in ["mUV", "mUVsq", "muvsq"] {
            evaluated = evaluated
                .replace(Atom::var(vakint_symbol!(name)))
                .with(Atom::one());
        }
        // master graphs of the data root carry the model's mass symbol
        evaluated = evaluated
            .replace(Atom::var(symbolica::symbol!("UFO::mUV")))
            .with(Atom::one());
        evaluated = evaluated
            .replace(Atom::var(vakint_symbol!(&settings.mu_r_sq_symbol)))
            .with(Atom::one());
        let outcome = engine.numerical_evaluation(
            &settings,
            evaluated.as_view(),
            &HashMap::default(),
            &HashMap::default(),
            None,
        );
        numerical = match outcome {
            Ok((result, _)) => json!({
                "seconds": started.elapsed().as_secs_f64(),
                "point": {"mUV": 1, "mUVsq": 1, "UFO::mUV": 1, "mu_r_squared": 1},
                "coefficients": result.coefficients().iter().map(|(power, c)| {
                    let p = c.precision();
                    json!({
                        "epsilon_power": power,
                        "real": c.value().re.to_string(),
                        "imaginary": c.value().im.to_string(),
                        "real_significant_digits": p.real().significant_digits(),
                        "imaginary_significant_digits": p.imaginary().significant_digits(),
                        "real_exact_zero": p.real().exact_zero(),
                        "imaginary_exact_zero": p.imaginary().exact_zero(),
                    })
                }).collect::<Vec<_>>(),
                "warnings": result.warnings().iter().map(ToString::to_string).collect::<Vec<_>>(),
            }),
            Err(error) => json!({"error": error.to_string()}),
        };
    }

    let statistics = if backend == "rustred" {
        let s = engine.ordered_cache_statistics();
        json!({"hits": s.hits, "misses": s.misses, "records": s.records})
    } else {
        Value::Null
    };
    write_json(
        &args.path("out")?,
        &json!({
            "schema": "vakint.benchmark-measurement.v1",
            "input_sha256": input_record["input_sha256"],
            "input_file_sha256": sha256_file(&input_file)?,
            "loop_count": loops,
            "backend": backend,
            "scope": scope,
            "cache": cache,
            "n_cores": n_cores,
            "ordered_weight_arithmetic": settings.ordered_weight_arithmetic.to_string(),
            "epsilon_terms": terms,
            "normalization": "MSbar",
            "loader_seconds": loader_seconds,
            "parse_seconds": parse_seconds,
            "call_seconds": call_seconds,
            "first_call": first_record,
            "reduction_capture": math,
            "process_seconds": process_start.elapsed().as_secs_f64(),
            "peak_rss_bytes_after_call": rss_after_call,
            "peak_rss_bytes": peak_rss_bytes(),
            "output_sha256": sha256_bytes(canonical.as_bytes()),
            "output_bytes": canonical.len(),
            "retained_rustred_masters": rustred_masters,
            "retained_form_symbols": form_masters,
            "output": if canonical.len() <= 200_000 { Value::from(canonical) } else { Value::Null },
            "ordered_step_cache": statistics,
            "numerical": numerical,
            "allowed_cpus": allowed_cpus(),
            "debug_assertions": cfg!(debug_assertions),
        }),
    )?;
    let _ = fs::remove_dir_all(&work);
    Ok(())
}

fn identity(args: &Args) -> Result<()> {
    let root = args.path("data-root")?;
    let mut loops = BTreeMap::new();
    for l in 1..=5 {
        let manifest = root.join(format!("{l}L/manifest.json"));
        if manifest.exists() {
            loops.insert(format!("{l}L"), sha256_file(&manifest)?);
        }
    }
    let exe = std::env::current_exe()?;
    write_json(
        &args.path("out")?,
        &json!({
            "schema": "vakint.benchmark-identity.v1",
            "vakint_version": env!("CARGO_PKG_VERSION"),
            "harness_source_sha256": sha256_bytes(SOURCE.as_bytes()),
            "executable_sha256": sha256_file(&exe)?,
            "debug_assertions": cfg!(debug_assertions),
            "experimental_rustred": cfg!(feature = "experimental-rustred"),
            "data_root": root,
            "root_manifest_sha256": sha256_file(&root.join("manifest.json")).ok(),
            "loop_manifest_sha256": loops,
            "default_n_cores": vakint::default_n_cores(),
            "default_ordered_weight_arithmetic": vakint::DEFAULT_ORDERED_WEIGHT_ARITHMETIC,
        }),
    )
}

fn main() -> Result<()> {
    let (command, args) = Args::parse(std::env::args().skip(1))?;
    match command.as_str() {
        "prepare-rung" => prepare_rung(&args),
        "prepare-graph" => prepare_graph(&args),
        "measure" => measure(&args),
        "identity" => identity(&args),
        "--help" | "-h" | "help" => {
            println!("{}", SOURCE.lines().take_while(|l| l.starts_with("//!")).map(|l| l.trim_start_matches("//!").trim_start()).collect::<Vec<_>>().join("\n"));
            Ok(())
        }
        other => Err(format!("unknown subcommand {other}").into()),
    }
}
