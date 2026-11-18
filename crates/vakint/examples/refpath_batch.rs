//! Batch driver of the reference path R2 (helper VH, numerator assignment,
//! lane `l5-stage1`; `vakint::refpath`): the PER-GRAPH results of a list of
//! background-field `B -> B` DOT graphs (e.g. one complete n_f class) through
//! the production pieces only. Class sums are NOT made here: they are made
//! from these per-graph records by the canonical hand-off script.
//!
//! Per graph: the input type on the DOT text
//! (`FactorizedNumerator::from_feynkit_dot`: the graph weight `overall_factor`
//! required and checked; colour split in invariants mode; scalars;
//! `--loop-sign-in-overall` for the reduced five-loop set), refpath's idenso
//! contraction of every colour branch (`refpath::factorized_branches` +
//! `refpath::contract`), the FeynKit ingress (`VakintExpression::from_diagram`,
//! `D -> 4 - 2ε`), and one `Vakint::evaluate` per master mode on a fresh
//! engine (direct lowering with the zero-sector pruning of
//! `VAKINT_P_ZERO_PRUNE`, one RustRed reduction per folded channel; Retain then
//! Symbolic: the request memo serves the second), observed by `p_gate`.
//!
//! ```text
//! refpath_batch <dir with model.json + diagram-NNNNNN.dot> --out <dir>
//!     (--graphs 1,2,7-9 | --graph-list <file with one number per line>)
//!     [--class <name>]           class label of the records (default nf<closed quark loops>)
//!     [--p-library <dir>]        RustRed p package (default $VAKINT_TEST_P_LIBRARY_1_4)
//!     [--projection metric|transverse|longitudinal|both]   (default metric)
//!     [--extra-terms <k>]        Laurent window eps^-L .. eps^(k-1) (default 0: the poles)
//!     [--modes retain,symbolic]  master modes, in this order (default both)
//!     [--loop-sign-in-overall]   FromFeynkitOptions::loop_sign_in_overall
//!     [--contract-only]          stop at the vakint input (no reduction)
//!     [--shared-engine]          one engine for all graphs of the run (warm caches; default: a fresh one per graph)
//!     [--no-gate]                no per-graph gate dump
//!     [--stop-file <path>]       stop between graphs while it exists (repeatable; exit 3)
//! ```
//!
//! REUSE IS BOUND TO THE CONFIGURATION: every output lives in
//! `<out>/cfg-<digest>/`, the digest of `config.json` = the binary's sha256,
//! `model.json`'s sha256, projection, colour treatment, weight source and
//! input-sign option, the package (`library.json` sha256) and the master
//! catalogue (file sha256), the eps window, the master modes and every other
//! `VAKINT_*` / `RUSTRED_*` setting that can change values (pruning mode
//! included). Another configuration is another directory (a miss, never a
//! hit). A per-graph record also binds the sha256 of its DOT file: a record of
//! another configuration or of a changed DOT is an ERROR.
//!
//! Per requested (graph, projection) the record `graphs/<id>.json` holds a
//! VALUE (per mode: value file, its sha256, for Symbolic the declared eps window,
//! the Laurent support found (which must lie inside it) and the value inline;
//! the Retain value keeps the reduced master-basis expression) or a
//! PROVED_ZERO entry with its reason. A failure is a FAILED record (never
//! retried by a later run: delete the record to retry), never a missing file
//! or a zero; a graph interrupted in an earlier run (`.started` without a
//! record) is recorded FAILED. Duplicate ids and `RUSTRED_RESUME_DIR` are
//! refused. `summary.json` says whether the request is complete.
use std::collections::BTreeMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use feynkit_kinematics::Kinematics;
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::symbol;
use vakint::diagram_integral::DiagramIntegralOptions;
use vakint::p_integral::FactorizedNumerator;
use vakint::p_integral::factorized::{FromFeynkitOptions, ProjectorChoice};
use vakint::refpath::{self, Projection};
use vakint::rustred_evaluation::p_gate;
use vakint::{
    EvaluationMethod, EvaluationOrder, RustRedEvaluationOptions, RustRedMasterMode, Vakint,
    VakintExpression,
};

const GRAPH_SCHEMA: &str = "refpath-batch.graph.v1";

fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect()
}

fn sha256_file(path: &Path) -> std::io::Result<String> {
    Ok(sha256_hex(&std::fs::read(path)?))
}

/// The CONTENT digest of a directory (lane handoff-fix's package_digest.py):
/// sha256 of the concatenated lines `"{relative path}\t{size}\t{sha256 of the
/// file}\n"` of every regular file below `root` (symbolic links to files
/// followed, a symbolic link to a directory is an error), sorted; with the
/// number of files and bytes.
fn content_digest(root: &Path) -> eyre::Result<(String, usize, u64)> {
    use std::io::Read as _;
    let (mut rows, mut total, mut stack) = (Vec::new(), 0u64, vec![root.to_path_buf()]);
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir)? {
            let path = entry?.path();
            let meta = std::fs::metadata(&path)?;
            if meta.is_dir() {
                eyre::ensure!(
                    !std::fs::symlink_metadata(&path)?.file_type().is_symlink(),
                    "{}: a symbolic link to a directory",
                    path.display()
                );
                stack.push(path);
                continue;
            }
            let (mut hasher, mut file, mut buffer) = (
                Sha256::new(),
                std::fs::File::open(&path)?,
                vec![0u8; 1 << 22],
            );
            loop {
                let n = file.read(&mut buffer)?;
                if n == 0 {
                    break;
                }
                hasher.update(&buffer[..n]);
            }
            let digest: String = hasher
                .finalize()
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect();
            let relative = path
                .strip_prefix(root)?
                .components()
                .map(|c| c.as_os_str().to_string_lossy().into_owned())
                .collect::<Vec<_>>()
                .join("/");
            total += meta.len();
            rows.push(format!("{relative}\t{}\t{digest}\n", meta.len()));
        }
    }
    rows.sort();
    Ok((sha256_hex(rows.concat().as_bytes()), rows.len(), total))
}

/// The canonical JSON text of the digest: object keys sorted, separators `,`
/// and `:`, no whitespace, strings as serde_json escapes them (= Python's
/// `json.dumps(v, sort_keys=True, separators=(',', ':'), ensure_ascii=False)`
/// for this content).
fn canonical_json(value: &Value) -> String {
    match value {
        Value::Object(map) => {
            let mut keys = map.keys().collect::<Vec<_>>();
            keys.sort();
            let items = keys
                .iter()
                .map(|k| {
                    format!(
                        "{}:{}",
                        serde_json::to_string(k).unwrap(),
                        canonical_json(&map[k.as_str()])
                    )
                })
                .collect::<Vec<_>>();
            format!("{{{}}}", items.join(","))
        }
        Value::Array(items) => format!(
            "[{}]",
            items
                .iter()
                .map(canonical_json)
                .collect::<Vec<_>>()
                .join(",")
        ),
        other => serde_json::to_string(other).unwrap(),
    }
}

fn parse_graphs(text: &str) -> Vec<usize> {
    let mut out = Vec::new();
    for part in text.split(',').filter(|part| !part.is_empty()) {
        match part.split_once('-') {
            Some((a, b)) => out.extend(a.parse::<usize>().unwrap()..=b.parse::<usize>().unwrap()),
            None => out.push(part.parse().unwrap()),
        }
    }
    out
}

fn cpu_seconds() -> f64 {
    // utime + stime of this process (clock ticks of 100 Hz).
    let stat = std::fs::read_to_string("/proc/self/stat").unwrap_or_default();
    let fields = stat
        .rsplit_once(')')
        .map(|(_, rest)| rest.split_whitespace().collect::<Vec<_>>())
        .unwrap_or_default();
    let tick = |i: usize| {
        fields
            .get(i)
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(0.0)
    };
    (tick(11) + tick(12)) / 100.0
}

/// Write `text` to `path` through a temporary file (unique per process) and a
/// rename.
fn write_atomic(path: &Path, text: &str) -> std::io::Result<()> {
    let tmp = path.with_file_name(format!(
        "{}.{}.tmp",
        path.file_name().unwrap().to_string_lossy(),
        std::process::id()
    ));
    std::fs::write(&tmp, text)?;
    std::fs::rename(tmp, path)
}

fn projection_names(projection: Projection) -> Vec<&'static str> {
    match projection {
        Projection::Both => vec!["transverse", "longitudinal"],
        other => vec![other.name()],
    }
}

fn mode_from(text: &str) -> eyre::Result<RustRedMasterMode> {
    match text {
        "retain" => Ok(RustRedMasterMode::Retain),
        "symbolic" => Ok(RustRedMasterMode::Symbolic),
        "numerical" => Ok(RustRedMasterMode::Numerical),
        other => eyre::bail!("unknown master mode {other}"),
    }
}

/// Everything a per-graph result depends on besides the graph itself.
#[derive(Clone, Debug)]
struct RunSettings {
    binary_sha256: String,
    model_sha256: String,
    projection: Projection,
    loop_sign_in_overall: bool,
    library: Option<PathBuf>,
    library_sha256: Option<String>,
    /// [`content_digest`] of the package directory (every file's bytes).
    library_content_sha256: Option<String>,
    catalogue_sha256: Option<String>,
    /// Where the five-loop master catalogue comes from (`VAKINT_P_MASTER_CATALOGUE`,
    /// the package's `p-master-catalogue.json`, or vakint's built-in one) and its
    /// data version.
    catalogue_source: String,
    catalogue_version: Option<String>,
    /// sha256 of the catalogue file in use (variable, package or built-in);
    /// `catalogue_sha256` keeps the reader's meaning: the variable's file only.
    catalogue_in_use_sha256: Option<String>,
    extra_terms: i64,
    modes: Vec<RustRedMasterMode>,
    contract_only: bool,
    /// One engine for every graph of the run (warm request memo / step
    /// cache) instead of a fresh one per graph; values are exact either way.
    shared_engine: bool,
    /// `VAKINT_*` / `RUSTRED_*` variables that can change values (all of them
    /// except output paths and diagnostics), as set.
    environment: BTreeMap<String, String>,
    /// For each bound variable whose value names an existing file or
    /// directory: the sha256 of the file / the [`content_digest`] of the directory.
    environment_content: BTreeMap<String, String>,
}

/// The variables of the environment that are bound into the configuration.
fn bound_environment() -> BTreeMap<String, String> {
    // output paths and diagnostics (the trace spill directory is scratch:
    // RustRed writes and removes its own files there; values do not depend on it)
    const NOT_BOUND: [&str; 4] = [
        "RUSTRED_ORDERED_MODULAR_STATS",
        "RUSTRED_TRACE_SPILL_DIR",
        "VAKINT_P_PHASES",
        "VAKINT_LOWERING_CENSUS",
    ];
    std::env::vars()
        .filter(|(name, _)| {
            (name.starts_with("VAKINT_") || name.starts_with("RUSTRED_"))
                && !name.starts_with("VAKINT_TEST_")
                && !NOT_BOUND.contains(&name.as_str())
        })
        .collect()
}

/// The configuration record (`config.json`) and its digest.
fn config(settings: &RunSettings) -> (Value, String) {
    let value = json!({
        "schema": "refpath-batch.config.v1",
        "binary_sha256": settings.binary_sha256,
        "model_sha256": settings.model_sha256,
        "projection": settings.projection.name(),
        "colour": "invariants (refpath symbols CA CF TF NA nf; higher invariants named)",
        "dimension": refpath::SYMBOLS.d.get_name(),
        "overall_factor": "DOT attribute, required and checked (OverallFactorSource::Dot)",
        "loop_sign_in_overall": settings.loop_sign_in_overall,
        "p_library": settings.library.as_ref().map(|p| p.display().to_string()),
        "p_library_json_sha256": settings.library_sha256,
        "p_library_content_sha256": settings.library_content_sha256,
        "master_catalogue_sha256": settings.catalogue_sha256,
        "master_catalogue_source": settings.catalogue_source,
        "master_catalogue_in_use_sha256": settings.catalogue_in_use_sha256,
        "master_catalogue_data_version": settings.catalogue_version,
        "eps_window": format!("eps^-L .. eps^{}", settings.extra_terms - 1),
        "extra_terms": settings.extra_terms,
        "zero_prune": std::env::var("VAKINT_P_ZERO_PRUNE").ok(),
        "modes": settings.modes.iter().map(|m| m.to_string()).collect::<Vec<_>>(),
        "contract_only": settings.contract_only,
        "shared_engine": settings.shared_engine,
        "environment": settings.environment,
        "environment_content_sha256": settings.environment_content,
    });
    let digest = sha256_hex(canonical_json(&value).as_bytes());
    (value, digest)
}

/// vakint's ε: the settings' name is already namespaced (`vakint::ε`).
fn epsilon_symbol() -> Symbol {
    let name = refpath::vakint_settings(1).epsilon_symbol;
    if name.contains("::") {
        symbol!(&name)
    } else {
        symbol!(&format!("vakint::{name}"))
    }
}

/// Variables of a Symbolic value outside its basis (ε, μ², π and this
/// chain's spectators `refpath::*`): a foreign symbol (e.g. another ε or an
/// unreplaced D) means a wrong value, not a valid one.
fn foreign_symbols(value: &Atom, epsilon: Symbol) -> Vec<String> {
    let mut out = std::collections::BTreeSet::new();
    value.visitor(&mut |part| {
        if let AtomView::Var(v) = part {
            let symbol = v.get_symbol();
            let name = symbol.get_name();
            if symbol != epsilon
                && symbol != Symbol::PI
                && !name.starts_with("refpath::")
                && name != "vakint::mursq"
            {
                out.insert(name.to_string());
            }
        }
        true
    });
    out.into_iter().collect()
}

/// The eps powers present in a value (`log(-p²/μ²)` counts as eps-free).
fn laurent_support(value: &Atom, epsilon: Symbol) -> Vec<i64> {
    let expanded = value.expand();
    let terms = match expanded.as_view() {
        AtomView::Add(add) => add.iter().map(|t| t.to_owned()).collect::<Vec<_>>(),
        view if view.is_zero() => Vec::new(),
        view => vec![view.to_owned()],
    };
    let mut orders = std::collections::BTreeSet::new();
    for term in terms {
        let factors = match term.as_view() {
            AtomView::Mul(mul) => mul.iter().map(|f| f.to_owned()).collect::<Vec<_>>(),
            view => vec![view.to_owned()],
        };
        let mut order = 0i64;
        for factor in factors {
            match factor.as_view() {
                AtomView::Var(v) if v.get_symbol() == epsilon => order += 1,
                AtomView::Pow(pow) => {
                    let (base, exp) = pow.get_base_exp();
                    if matches!(base, AtomView::Var(v) if v.get_symbol() == epsilon) {
                        order += i64::try_from(exp).unwrap_or(i64::MIN / 2);
                    }
                }
                _ => {}
            }
        }
        orders.insert(order);
    }
    orders.into_iter().collect()
}

/// One graph: the input type -> contraction -> vakint inputs.
struct Prepared {
    record: Value,
    loops: usize,
    /// per projection: the vakint input, or why it is zero
    inputs: Vec<(&'static str, Result<Atom, String>)>,
}

fn prepare(
    model: &Arc<feynkit_model::Model>,
    dot: &str,
    settings: &RunSettings,
    epsilon: Symbol,
    seconds: &mut BTreeMap<&'static str, f64>,
) -> eyre::Result<Prepared> {
    let s = &*refpath::SYMBOLS;
    let names = projection_names(settings.projection);
    let zero_all = |record: Value, loops: usize, reason: &str| Prepared {
        record,
        loops,
        inputs: names
            .iter()
            .map(|name| (*name, Err(reason.to_string())))
            .collect(),
    };
    let started = Instant::now();
    let options = FromFeynkitOptions {
        dimension: s.d,
        projector: ProjectorChoice::None,
        loop_sign_in_overall: settings.loop_sign_in_overall,
        ..FromFeynkitOptions::default()
    };
    // The input type on the DOT text: parse, weight required and checked,
    // colour split, scalars, validation.
    let (diagram, values) = FactorizedNumerator::from_feynkit_dot(dot, model, &options)?;
    seconds.insert("parse_input_type", started.elapsed().as_secs_f64());
    eyre::ensure!(
        diagram.projector().expression().is_one()
            && diagram.numerator_prefactor().expression().is_one(),
        "the reference path expects projector 1 and numerator_prefactor 1"
    );
    let loops = diagram.loop_count();
    let quark_loops = refpath::closed_quark_loops(&diagram)?;
    let overall =
        feynkit_graph::expressions::evaluate_overall_factor(diagram.overall_factor().as_view());
    let mut record = json!({
        "name": diagram.name(), "loops": loops, "quark_loops": quark_loops,
        "overall_factor_dot": diagram.overall_factor().to_canonical_string(),
        "overall_factor": overall.to_canonical_string(),
        "colour_branches": values.len(),
    });
    if refpath::scaleless_contact(&diagram) {
        let reason = "scaleless: no internal edge carries the external momentum";
        return Ok(zero_all(record, loops, reason));
    }
    if values.is_empty() {
        let reason = "zero colour factor: every colour branch vanishes";
        return Ok(zero_all(record, loops, reason));
    }
    // The input type's weight equals the chain's reading of the graph.
    let weight = refpath::factorized_weight(&values)?;
    let chain_weight = &overall * Atom::var(s.nf).pow(Atom::num(quark_loops as i64));
    eyre::ensure!(
        (&weight - &chain_weight).expand().is_zero(),
        "input-type weight {weight} differs from the graph's {chain_weight}"
    );
    record["weight"] = json!(weight.to_canonical_string());
    record["loop_sign"] = json!(values[0].prefactor.loop_sign);
    record["nf_power"] = json!(values[0].prefactor.nf_power);
    let branches = refpath::factorized_branches(&values)?;
    record["branch_colours"] = json!(
        values
            .iter()
            .map(|v| refpath::factorized_colour(v).to_canonical_string())
            .collect::<Vec<_>>()
    );
    record["branch_coefficients"] = json!(
        branches
            .iter()
            .map(|b| b.coefficient.to_canonical_string())
            .collect::<Vec<_>>()
    );

    let dimension = Atom::var(s.d);
    let kinematics = Kinematics::in_dimension(&dimension)?;
    let family = diagram.propagator_family(&kinematics)?;
    let [momentum] = family.external_momenta() else {
        eyre::bail!("a two-point graph needs one independent external momentum");
    };
    let momentum = momentum.clone();
    let q_squared = Kinematics::in_dimension(&dimension)?
        .with_momenta([momentum.clone()])?
        .scalar_product(&momentum, &momentum)?;
    let started = Instant::now();
    let numerators = refpath::contract(
        &diagram,
        &branches,
        &momentum,
        &q_squared,
        settings.projection,
    )?;
    seconds.insert("contract", started.elapsed().as_secs_f64());

    let started = Instant::now();
    let in_epsilon = Atom::num(4) - Atom::num(2) * Atom::var(epsilon);
    let mut inputs = Vec::new();
    for (name, numerator) in numerators {
        let numerator = numerator * &weight;
        if numerator.expand().is_zero() {
            inputs.push((name, Err("the contracted numerator is 0".to_string())));
            continue;
        }
        let input = VakintExpression::from_diagram(
            &diagram,
            &family,
            &numerator,
            &DiagramIntegralOptions::default(),
        )?;
        let input = Atom::from(input)
            .replace(Atom::var(s.d))
            .with(in_epsilon.clone());
        inputs.push((name, Ok(input)));
    }
    seconds.insert("ingress", started.elapsed().as_secs_f64());
    Ok(Prepared {
        record,
        loops,
        inputs,
    })
}

/// RustRed's per-reduction statistics lines appended to
/// `RUSTRED_ORDERED_MODULAR_STATS` after byte `from`: (file length, summary).
fn rustred_stats_since(from: u64) -> (u64, Value) {
    let Some(path) = std::env::var_os("RUSTRED_ORDERED_MODULAR_STATS") else {
        return (0, json!(null));
    };
    let text = std::fs::read(&path).unwrap_or_default();
    let end = text.len() as u64;
    let new = String::from_utf8_lossy(&text[(from.min(end) as usize)..]).to_string();
    let (mut reductions, mut ops, mut record, mut replay) = (0, 0u64, 0.0, 0.0);
    for line in new.lines() {
        if let Ok(d) = serde_json::from_str::<Value>(line)
            && d.get("arithmetic").is_some()
        {
            reductions += 1;
            ops += d["ops"].as_u64().unwrap_or(0);
            record += d["record_seconds"].as_f64().unwrap_or(0.0);
            replay += d["replay_seconds"].as_f64().unwrap_or(0.0);
        }
    }
    let summary = json!({"reductions": reductions, "ops": ops,
        "record_seconds": record, "replay_seconds": replay});
    (end, summary)
}

/// Where the outputs of one graph go.
struct Outputs<'a> {
    values: &'a Path,
    gate: Option<&'a Path>,
    id: &'a str,
}

/// Evaluate the vakint inputs of one graph: per projection, per mode.
fn evaluate(
    template: Option<&Vakint>,
    prepared: &Prepared,
    settings: &RunSettings,
    epsilon: Symbol,
    outputs: &Outputs<'_>,
) -> eyre::Result<Vec<Value>> {
    let loops = prepared.loops as i64;
    let window = (-loops, settings.extra_terms - 1);
    let id = outputs.id;
    let mut entries = Vec::new();
    for (name, input) in &prepared.inputs {
        let input = match input {
            Err(reason) => {
                entries
                    .push(json!({"projection": name, "status": "PROVED_ZERO", "reason": reason}));
                continue;
            }
            Ok(input) => input,
        };
        if settings.contract_only {
            let text = input.to_canonical_string();
            let file = format!("{id}.{name}.input");
            write_atomic(&outputs.values.join(&file), &text)?;
            entries.push(json!({"projection": name, "status": "INPUT", "file": file,
                "sha256": sha256_hex(text.as_bytes()), "bytes": text.len()}));
            continue;
        }
        let template = template.ok_or_else(|| eyre::eyre!("no engine"))?;
        let fresh;
        let engine = if settings.shared_engine {
            template
        } else {
            fresh = refpath::fresh_engine(template)?;
            &fresh
        };
        let mut modes = serde_json::Map::new();
        for (index, mode) in settings.modes.iter().enumerate() {
            let mut vakint_settings = refpath::vakint_settings(loops + settings.extra_terms);
            vakint_settings.evaluation_order =
                EvaluationOrder(vec![EvaluationMethod::RustRed(RustRedEvaluationOptions {
                    master_mode: *mode,
                })]);
            let (before, _) = rustred_stats_since(u64::MAX);
            let cpu_start = cpu_seconds();
            let started = Instant::now();
            let (value, observation) =
                p_gate::observe(|| engine.evaluate(&vakint_settings, input.as_view()));
            let wall = started.elapsed().as_secs_f64();
            let cpu = cpu_seconds() - cpu_start;
            let value = value.map_err(|e| eyre::eyre!("{name} {mode}: {e}"))?;
            let reduce = observation.prepare_seconds + observation.library_seconds;
            let finalise = observation.render_seconds + observation.finalize_seconds;
            let (_, rustred) = rustred_stats_since(before);
            let text = value.to_canonical_string();
            let file = format!("{id}.{name}.{mode}.txt");
            write_atomic(&outputs.values.join(&file), &text)?;
            let mut entry = json!({
                "file": file, "sha256": sha256_hex(text.as_bytes()), "bytes": text.len(),
                "wall": wall, "cpu": cpu, "lower": wall - reduce - finalise - observation.observer_seconds,
                "reduce": reduce, "finalise": finalise, "reductions": observation.reductions,
                "pairs": observation.pairs, "folded_channels": observation.folded,
                "request_memo_hits": observation.memo_hits, "direct_lowered": observation.direct_lowered,
                "direct_declined": observation.direct_declined, "rustred": rustred,
                "step_memo": "not exposed by rustred d5b63a659 (only with its opt-in fwd-profile feature)",
            });
            if *mode == RustRedMasterMode::Retain {
                entry["window"] = json!("exact in eps: the reduced master-basis expression");
            } else {
                let foreign = foreign_symbols(&value, epsilon);
                eyre::ensure!(
                    foreign.is_empty(),
                    "{name} {mode}: symbols foreign to the value's basis: {foreign:?}"
                );
                let support = laurent_support(&value, epsilon);
                entry["window"] = json!([window.0, window.1]);
                entry["support"] = json!(support);
                eyre::ensure!(
                    support.iter().all(|k| *k >= window.0 && *k <= window.1),
                    "{name} {mode}: Laurent support {support:?} outside the declared window {window:?}"
                );
                entry["oep"] = json!(text.contains("Oep("));
                entry["residual"] = json!(text.contains("PResidual("));
                entry["value"] = json!(text);
            }
            if index == 0 {
                let (view, count) = observation.view_without_zero_sectors();
                entry["keys"] = json!(count.total);
                entry["keys_zero_sector"] = json!(count.zero);
                entry["keys_census_unknown"] = json!(count.unknown);
                entry["keys_kept"] = json!(view.len());
                if observation.reductions == 0 {
                    entry["zero_by"] = json!(format!(
                        "no reduction: {} keys, {} in zero sectors (pruned)",
                        count.total, count.zero
                    ));
                }
                if let Some(gate) = outputs.gate {
                    let header = BTreeMap::from([
                        ("graph".to_string(), id.to_string()),
                        ("projection".to_string(), name.to_string()),
                        ("loops".to_string(), loops.to_string()),
                        ("zero_sectors_dropped".to_string(), count.zero.to_string()),
                        ("keys_before_drop".to_string(), count.total.to_string()),
                        ("census_unknown_kept".to_string(), count.unknown.to_string()),
                        (
                            "producer".to_string(),
                            "vakint refpath_batch (helper VH lane l5-stage1)".to_string(),
                        ),
                    ]);
                    write_atomic(
                        &gate.join(format!("{id}.{name}.gate")),
                        &p_gate::write_gate_view(&view, &header),
                    )?;
                }
            }
            modes.insert(mode.to_string(), entry);
        }
        entries.push(json!({"projection": name, "status": "VALUE", "modes": modes}));
    }
    Ok(entries)
}

/// The five-loop master catalogue vakint will use, as vakint chooses it
/// (`Vakint::p_master_catalogue`): `VAKINT_P_MASTER_CATALOGUE` (a file or `off`),
/// else the package's `p-master-catalogue.json`, else the built-in one compiled
/// into this binary: (source, sha256 of the file's bytes, its data version).
fn catalogue_identity(library: Option<&Path>) -> eyre::Result<CatalogueIdentity> {
    const BUILT_IN: &[u8] = include_bytes!("../data/p_master/p-master-catalogue-5l.json");
    let version = |bytes: &[u8]| {
        serde_json::from_slice::<Value>(bytes)
            .ok()
            .and_then(|v| v["data_version"].as_str().map(str::to_string))
    };
    if let Some(path) = std::env::var_os("VAKINT_P_MASTER_CATALOGUE") {
        if path == "off" {
            return Ok(CatalogueIdentity {
                source: "VAKINT_P_MASTER_CATALOGUE=off".into(),
                variable_sha256: Some("off".into()),
                in_use_sha256: None,
                version: None,
            });
        }
        let bytes = std::fs::read(&path)?;
        return Ok(CatalogueIdentity {
            source: format!("VAKINT_P_MASTER_CATALOGUE={}", Path::new(&path).display()),
            variable_sha256: Some(sha256_hex(&bytes)),
            in_use_sha256: Some(sha256_hex(&bytes)),
            version: version(&bytes),
        });
    }
    if let Some(file) = library
        .map(|dir| dir.join("p-master-catalogue.json"))
        .filter(|f| f.is_file())
    {
        let bytes = std::fs::read(&file)?;
        return Ok(CatalogueIdentity {
            source: format!("package {}", file.display()),
            variable_sha256: None,
            in_use_sha256: Some(sha256_hex(&bytes)),
            version: version(&bytes),
        });
    }
    Ok(CatalogueIdentity {
        source: "built-in (crates/vakint/data/p_master/p-master-catalogue-5l.json)".into(),
        variable_sha256: None,
        in_use_sha256: Some(sha256_hex(BUILT_IN)),
        version: version(BUILT_IN),
    })
}

/// The master catalogue vakint uses, see [`catalogue_identity`].
struct CatalogueIdentity {
    source: String,
    /// sha256 of the file of `VAKINT_P_MASTER_CATALOGUE` ("off" when off;
    /// `None` when unset): the reader re-hashes that file.
    variable_sha256: Option<String>,
    in_use_sha256: Option<String>,
    version: Option<String>,
}

fn main() -> eyre::Result<()> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    let mut directory = None;
    let mut out = None;
    let mut graphs = None;
    let mut list_file = None;
    let mut class = None;
    let mut projection = Projection::Metric;
    let mut extra_terms = 0i64;
    let mut modes = vec![RustRedMasterMode::Retain, RustRedMasterMode::Symbolic];
    let mut library = std::env::var_os("VAKINT_TEST_P_LIBRARY_1_4").map(PathBuf::from);
    let mut loop_sign_in_overall = false;
    let mut contract_only = false;
    let mut shared_engine = false;
    let mut gate = true;
    let mut stop_files = Vec::new();
    let mut index = 0;
    while index < args.len() {
        let value = || args.get(index + 1).cloned().expect("missing value");
        match args[index].as_str() {
            "--out" => {
                out = Some(PathBuf::from(value()));
                index += 1;
            }
            "--graphs" => {
                graphs = Some(parse_graphs(&value()));
                index += 1;
            }
            "--graph-list" => {
                let path = PathBuf::from(value());
                let text = std::fs::read_to_string(&path)?;
                graphs = Some(
                    text.lines()
                        .filter(|line| !line.trim().is_empty() && !line.starts_with('#'))
                        .map(|line| line.trim().parse().unwrap())
                        .collect(),
                );
                list_file = Some(path);
                index += 1;
            }
            "--class" => {
                class = Some(value());
                index += 1;
            }
            "--projection" => {
                projection = value().parse().map_err(|e: String| eyre::eyre!(e))?;
                index += 1;
            }
            "--extra-terms" => {
                extra_terms = value().parse()?;
                index += 1;
            }
            "--modes" => {
                modes = value()
                    .split(',')
                    .map(mode_from)
                    .collect::<eyre::Result<_>>()?;
                index += 1;
            }
            "--p-library" => {
                library = Some(PathBuf::from(value()));
                index += 1;
            }
            "--stop-file" => {
                stop_files.push(PathBuf::from(value()));
                index += 1;
            }
            "--loop-sign-in-overall" => loop_sign_in_overall = true,
            "--contract-only" => contract_only = true,
            "--shared-engine" => shared_engine = true,
            "--no-gate" => gate = false,
            other if directory.is_none() && !other.starts_with("--") => {
                directory = Some(PathBuf::from(other))
            }
            other => eyre::bail!("unknown argument {other}"),
        }
        index += 1;
    }
    eyre::ensure!(
        std::env::var_os("RUSTRED_RESUME_DIR").is_none(),
        "RUSTRED_RESUME_DIR is set: refused"
    );
    let directory = directory.ok_or_else(|| eyre::eyre!("no DOT directory"))?;
    let out = out.ok_or_else(|| eyre::eyre!("no --out"))?;
    let graphs = graphs.ok_or_else(|| eyre::eyre!("no --graphs / --graph-list"))?;
    let catalogue = catalogue_identity(library.as_deref())?;
    let started = Instant::now();
    let library_content_sha256 = library
        .as_deref()
        .map(content_digest)
        .transpose()?
        .map(|(digest, files, bytes)| {
            eprintln!(
                "[refpath_batch] package content digest {digest}: {files} files, {bytes} bytes, {:.1} s",
                started.elapsed().as_secs_f64()
            );
            digest
        });
    let environment = bound_environment();
    let mut environment_content = BTreeMap::new();
    for (name, value) in &environment {
        let path = Path::new(value);
        if path.is_file() {
            environment_content.insert(name.clone(), sha256_file(path)?);
        } else if path.is_dir() {
            environment_content.insert(name.clone(), content_digest(path)?.0);
        }
    }
    let settings = RunSettings {
        binary_sha256: sha256_file(&std::env::current_exe()?)?,
        model_sha256: sha256_file(&directory.join("model.json"))?,
        projection,
        loop_sign_in_overall,
        library_sha256: library
            .as_ref()
            .map(|dir| sha256_file(&dir.join("library.json")))
            .transpose()?,
        library,
        library_content_sha256,
        catalogue_sha256: catalogue.variable_sha256,
        catalogue_source: catalogue.source,
        catalogue_version: catalogue.version,
        catalogue_in_use_sha256: catalogue.in_use_sha256,
        extra_terms,
        modes,
        contract_only,
        shared_engine,
        environment,
        environment_content,
    };
    let request = Request {
        directory: &directory,
        out: &out,
        graphs: &graphs,
        list_file: list_file.as_deref(),
        class: class.as_deref(),
        gate,
        stop_files: &stop_files,
    };
    let summary = run(&request, &settings)?;
    eprintln!("[refpath_batch] {summary}");
    if summary["failed"].as_u64().unwrap_or(0) > 0 {
        std::process::exit(2);
    }
    if summary["complete"] != true {
        std::process::exit(3);
    }
    Ok(())
}

/// What is asked of one run.
struct Request<'a> {
    directory: &'a Path,
    out: &'a Path,
    graphs: &'a [usize],
    list_file: Option<&'a Path>,
    class: Option<&'a str>,
    gate: bool,
    stop_files: &'a [PathBuf],
}

/// The batch: per graph, reuse a record of THIS configuration and DOT, or
/// compute one; returns the summary (also written to `summary.json`).
fn run(request: &Request<'_>, settings: &RunSettings) -> eyre::Result<Value> {
    let graphs = request.graphs;
    let mut ids = graphs.to_vec();
    ids.sort_unstable();
    let distinct = ids.len();
    ids.dedup();
    eyre::ensure!(
        ids.len() == distinct,
        "the graph list repeats {} id(s): refused",
        distinct - ids.len()
    );
    let (config_value, digest) = config(settings);
    let cfg_dir = request.out.join(format!("cfg-{}", &digest[..16]));
    for sub in ["graphs", "values", "gate"] {
        std::fs::create_dir_all(cfg_dir.join(sub))?;
    }
    let config_path = cfg_dir.join("config.json");
    let config_text = serde_json::to_string_pretty(&config_value)?;
    if config_path.exists() {
        eyre::ensure!(
            std::fs::read_to_string(&config_path)? == config_text,
            "{} holds another configuration: refused",
            config_path.display()
        );
    } else {
        write_atomic(&config_path, &config_text)?;
    }
    let model = refpath::load_model(request.directory)?;
    let epsilon = epsilon_symbol();
    let mut log = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(cfg_dir.join("graphs.jsonl"))?;
    let mut template = None;
    let started_run = Instant::now();
    let (mut reused, mut computed, mut failed, mut stopped) = (0usize, 0usize, Vec::new(), false);
    // per-stage totals of the graphs computed in this run (one DOT set)
    let mut totals = BTreeMap::<String, f64>::new();
    let mut peak = 0u64;
    for (position, &graph) in graphs.iter().enumerate() {
        let id = format!("{graph:06}");
        let dot_path = request.directory.join(format!("diagram-{id}.dot"));
        let dot = std::fs::read_to_string(&dot_path)?;
        let dot_sha256 = sha256_hex(dot.as_bytes());
        let record_path = cfg_dir.join("graphs").join(format!("{id}.json"));
        let started_path = cfg_dir.join("graphs").join(format!("{id}.started"));
        if record_path.exists() {
            let record: Value = serde_json::from_str(&std::fs::read_to_string(&record_path)?)?;
            eyre::ensure!(
                record["schema"] == GRAPH_SCHEMA
                    && record["cfg_digest"] == digest.as_str()
                    && record["dot_sha256"] == dot_sha256.as_str(),
                "{}: a record of another configuration or DOT file (dot sha256 {} vs {dot_sha256}): refused",
                record_path.display(),
                record["dot_sha256"]
            );
            if record["status"] == "FAILED" {
                failed.push(graph);
            } else {
                reused += 1;
            }
            continue;
        }
        if let Some(stop) = request.stop_files.iter().find(|path| path.exists()) {
            eprintln!(
                "[refpath_batch] stop file {} exists: stopping before graph {id}",
                stop.display()
            );
            stopped = true;
            break;
        }
        let mut record = json!({
            "schema": GRAPH_SCHEMA, "graph": graph, "id": id, "dot_file": dot_path.display().to_string(),
            "dot_sha256": dot_sha256, "model_sha256": settings.model_sha256, "cfg_digest": digest,
            "binary_sha256": settings.binary_sha256, "projection": settings.projection.name(),
            "master_catalogue_sha256": settings.catalogue_sha256,
            "master_catalogue_in_use_sha256": settings.catalogue_in_use_sha256,
            "master_catalogue_data_version": settings.catalogue_version,
            "colour_branch": "all: the record is the sum over the graph's non-zero colour branches",
        });
        if started_path.exists() {
            record["status"] = json!("FAILED");
            record["error"] = json!(
                "interrupted during an earlier run (crash, kill or out of memory): not retried"
            );
            write_atomic(&record_path, &record.to_string())?;
            std::fs::remove_file(&started_path)?;
            failed.push(graph);
            continue;
        }
        std::fs::write(&started_path, std::process::id().to_string())?;
        if template.is_none() && !settings.contract_only {
            let mut engine = Vakint::new()?;
            if let Some(library) = &settings.library {
                engine.load_p_library_directory(library, false)?;
            }
            template = Some(engine);
        }
        refpath::start_peak_rss();
        let cpu_start = cpu_seconds();
        let started = Instant::now();
        let mut seconds = BTreeMap::new();
        let gate_dir = cfg_dir.join("gate");
        let values_dir = cfg_dir.join("values");
        let outputs = Outputs {
            values: &values_dir,
            gate: request.gate.then_some(gate_dir.as_path()),
            id: &id,
        };
        let result = prepare(&model, &dot, settings, epsilon, &mut seconds).and_then(|prepared| {
            let entries = evaluate(template.as_ref(), &prepared, settings, epsilon, &outputs)?;
            Ok((prepared.record, entries))
        });
        record["seconds"] = json!(seconds);
        record["wall"] = json!(started.elapsed().as_secs_f64());
        record["cpu"] = json!(cpu_seconds() - cpu_start);
        let rss = refpath::peak_rss_bytes().unwrap_or(0);
        record["peak_rss_bytes"] = json!(rss);
        peak = peak.max(rss);
        match result {
            Ok((graph_record, entries)) => {
                for (key, value) in graph_record.as_object().unwrap() {
                    record[key] = value.clone();
                }
                record["class"] = json!(request.class.map_or_else(
                    || format!("nf{}", graph_record["quark_loops"]),
                    str::to_string
                ));
                for (stage, value) in &seconds {
                    *totals.entry((*stage).to_string()).or_default() += value;
                }
                for entry in &entries {
                    for (mode, m) in entry["modes"].as_object().into_iter().flatten() {
                        for stage in ["lower", "reduce", "finalise", "cpu"] {
                            *totals.entry(format!("{mode}.{stage}")).or_default() +=
                                m[stage].as_f64().unwrap_or(0.0);
                        }
                    }
                }
                record["projections"] = json!(entries);
                record["status"] = json!("OK");
                computed += 1;
            }
            Err(error) => {
                record["status"] = json!("FAILED");
                record["error"] = json!(format!("{error:?}"));
                failed.push(graph);
            }
        }
        *totals.entry("cpu".to_string()).or_default() += record["cpu"].as_f64().unwrap_or(0.0);
        write_atomic(&record_path, &record.to_string())?;
        std::fs::remove_file(&started_path)?;
        let line = json!({"graph": graph, "status": record["status"], "wall": record["wall"],
            "cpu": record["cpu"], "peak_rss_bytes": record["peak_rss_bytes"], "seconds": record["seconds"]});
        log.write_all(format!("{line}\n").as_bytes())?;
        eprintln!("[refpath_batch] {}/{} {line}", position + 1, graphs.len());
        let elapsed = started_run.elapsed().as_secs_f64();
        let per_graph = if computed > 0 {
            elapsed / computed as f64
        } else {
            0.0
        };
        let progress = json!({
            "updated_unix": std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)?.as_secs(),
            "listed": graphs.len(), "reused": reused, "computed": computed, "failed": failed.len(),
            "remaining": graphs.len() - position - 1, "elapsed_seconds": elapsed,
            "seconds_per_computed_graph": per_graph,
            "stage_totals_seconds": totals, "peak_rss_bytes": peak,
        });
        write_atomic(&cfg_dir.join("progress.json"), &progress.to_string())?;
    }
    let summary = json!({
        "cfg_dir": cfg_dir.display().to_string(), "cfg_digest": digest,
        "dot_directory": request.directory.display().to_string(), "listed": graphs.len(),
        "list_file": request.list_file.map(|p| p.display().to_string()),
        "list_file_sha256": request.list_file.map(sha256_file).transpose()?,
        "ids_sha256": sha256_hex(ids.iter().map(|g| format!("{g}\n")).collect::<String>().as_bytes()),
        "reused": reused, "computed": computed, "failed": failed.len(), "failed_graphs": failed,
        "stopped_by_stop_file": stopped,
        "complete": !stopped && failed.is_empty() && reused + computed == graphs.len(),
        "stage_totals_seconds_this_run": totals, "peak_rss_bytes_this_run": peak,
    });
    write_atomic(&cfg_dir.join("summary.json"), &summary.to_string())?;
    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn settings(projection: Projection) -> RunSettings {
        RunSettings {
            binary_sha256: "b".into(),
            model_sha256: "m".into(),
            projection,
            loop_sign_in_overall: false,
            library: None,
            library_sha256: None,
            library_content_sha256: None,
            catalogue_sha256: None,
            catalogue_source: "none".into(),
            catalogue_version: None,
            catalogue_in_use_sha256: None,
            extra_terms: 0,
            modes: vec![RustRedMasterMode::Retain, RustRedMasterMode::Symbolic],
            contract_only: true,
            shared_engine: false,
            environment: BTreeMap::new(),
            environment_content: BTreeMap::new(),
        }
    }

    /// Every bound field changes the configuration digest (another
    /// configuration = another directory: a miss, never a hit).
    #[test]
    fn every_bound_setting_changes_the_configuration() {
        let base = settings(Projection::Metric);
        let digest = |s: &RunSettings| config(s).1;
        let d0 = digest(&base);
        let variants: Vec<(&str, RunSettings)> = vec![
            ("projection", settings(Projection::Transverse)),
            (
                "binary",
                RunSettings {
                    binary_sha256: "b2".into(),
                    ..base.clone()
                },
            ),
            (
                "model",
                RunSettings {
                    model_sha256: "m2".into(),
                    ..base.clone()
                },
            ),
            (
                "input sign",
                RunSettings {
                    loop_sign_in_overall: true,
                    ..base.clone()
                },
            ),
            (
                "package",
                RunSettings {
                    library_sha256: Some("p".into()),
                    ..base.clone()
                },
            ),
            (
                "package content",
                RunSettings {
                    library_content_sha256: Some("pc".into()),
                    ..base.clone()
                },
            ),
            (
                "content of a bound path variable",
                RunSettings {
                    environment_content: BTreeMap::from([("X".into(), "h".into())]),
                    ..base.clone()
                },
            ),
            (
                "catalogue",
                RunSettings {
                    catalogue_sha256: Some("c".into()),
                    ..base.clone()
                },
            ),
            (
                "eps depth",
                RunSettings {
                    extra_terms: 1,
                    ..base.clone()
                },
            ),
            (
                "modes",
                RunSettings {
                    modes: vec![RustRedMasterMode::Symbolic],
                    ..base.clone()
                },
            ),
            (
                "contract only",
                RunSettings {
                    contract_only: false,
                    ..base.clone()
                },
            ),
            (
                "shared engine",
                RunSettings {
                    shared_engine: true,
                    ..base.clone()
                },
            ),
            (
                "environment (pruning mode VAKINT_P_ZERO_PRUNE and the like)",
                RunSettings {
                    environment: BTreeMap::from([("VAKINT_P_ZERO_PRUNE".into(), "early".into())]),
                    ..base.clone()
                },
            ),
        ];
        for (what, variant) in variants {
            assert_ne!(digest(&variant), d0, "{what} is not bound");
        }
        assert_eq!(digest(&base.clone()), d0);
    }

    fn l1() -> Option<PathBuf> {
        let dir = PathBuf::from(
            std::env::var("VAKINT_TEST_FACTORIZED_DOT_ROOT")
                .unwrap_or_else(|_| "/common/ben/beta/results/prepared".into()),
        )
        .join("L1");
        if dir.join("model.json").is_file() {
            Some(dir)
        } else {
            eprintln!("SKIP: no L1 DOTs under {}", dir.display());
            None
        }
    }

    /// (projection, status, sha256 of the stored input) per projection.
    fn entries(cfg: &Path, graph: usize) -> Vec<(String, String, String)> {
        let record: Value = serde_json::from_str(
            &std::fs::read_to_string(cfg.join("graphs").join(format!("{graph:06}.json"))).unwrap(),
        )
        .unwrap();
        assert_eq!(record["status"], "OK", "{record}");
        record["projections"]
            .as_array()
            .unwrap()
            .iter()
            .map(|p| {
                let text = |key: &str| p[key].as_str().unwrap_or("").to_string();
                (text("projection"), text("status"), text("sha256"))
            })
            .collect()
    }

    fn request<'a>(dir: &'a Path, out: &'a Path, graphs: &'a [usize]) -> Request<'a> {
        Request {
            directory: dir,
            out,
            graphs,
            list_file: None,
            class: None,
            gate: false,
            stop_files: &[],
        }
    }

    /// ACCEPTANCE (reviewer's control + cold vs warm): metric on graphs {0,1},
    /// then transverse on {0,1,2} and metric on {0,1,2} in the SAME output
    /// directory (warm) against cold runs in fresh directories: the warm records
    /// equal the cold ones entry by entry (status and input digest per
    /// projection), every requested (graph, projection) has an entry, nothing
    /// of another configuration is reused; a record planted from another
    /// configuration, a changed DOT and a repeated id are refused.
    #[test]
    fn warm_runs_with_a_changed_projection_and_graph_set_equal_cold_runs() {
        let Some(dir) = l1() else { return };
        let root = std::env::temp_dir().join(format!("refpath_batch_test_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let (warm, cold) = (root.join("warm"), root.join("cold"));
        let metric = settings(Projection::Metric);
        let transverse = settings(Projection::Transverse);
        let first = run(&request(&dir, &warm, &[0, 1]), &metric).unwrap();
        assert_eq!(
            (first["computed"].as_u64(), first["complete"].as_bool()),
            (Some(2), Some(true))
        );
        let all = [0usize, 1, 2];
        let t_warm = run(&request(&dir, &warm, &all), &transverse).unwrap();
        assert_eq!(
            t_warm["computed"], 3,
            "transverse reused a metric record: {t_warm}"
        );
        let m_warm = run(&request(&dir, &warm, &all), &metric).unwrap();
        assert_eq!(
            (m_warm["reused"].as_u64(), m_warm["computed"].as_u64()),
            (Some(2), Some(1))
        );
        for (settings, warm_summary, name) in [
            (&transverse, &t_warm, "transverse"),
            (&metric, &m_warm, "metric"),
        ] {
            let cold_summary = run(&request(&dir, &cold.join(name), &all), settings).unwrap();
            let (w, c) = (
                PathBuf::from(warm_summary["cfg_dir"].as_str().unwrap()),
                PathBuf::from(cold_summary["cfg_dir"].as_str().unwrap()),
            );
            for graph in all {
                let warm_entries = entries(&w, graph);
                assert_eq!(warm_entries, entries(&c, graph), "graph {graph} {name}");
                assert_eq!(warm_entries.len(), 1);
                assert_eq!(warm_entries[0].0, name);
                assert!(["INPUT", "PROVED_ZERO"].contains(&warm_entries[0].1.as_str()));
            }
        }
        // A metric record planted under the transverse configuration: refused.
        let (m_dir, t_dir) = (
            PathBuf::from(m_warm["cfg_dir"].as_str().unwrap()),
            PathBuf::from(t_warm["cfg_dir"].as_str().unwrap()),
        );
        std::fs::copy(
            m_dir.join("graphs/000001.json"),
            t_dir.join("graphs/000001.json"),
        )
        .unwrap();
        let error = run(&request(&dir, &warm, &all), &transverse).unwrap_err();
        assert!(
            error.to_string().contains("another configuration"),
            "{error}"
        );
        // A repeated id: refused.
        let error = run(&request(&dir, &cold.join("dup"), &[0, 0]), &metric).unwrap_err();
        assert!(error.to_string().contains("repeats"), "{error}");
        // A changed DOT under the same name: refused.
        let copy = root.join("dots");
        std::fs::create_dir_all(&copy).unwrap();
        for name in ["model.json", "diagram-000000.dot"] {
            std::fs::copy(dir.join(name), copy.join(name)).unwrap();
        }
        let out = root.join("changed");
        assert_eq!(
            run(&request(&copy, &out, &[0]), &metric).unwrap()["computed"],
            1
        );
        let text = std::fs::read_to_string(copy.join("diagram-000000.dot")).unwrap();
        std::fs::write(copy.join("diagram-000000.dot"), format!("{text}\n")).unwrap();
        let error = run(&request(&copy, &out, &[0]), &metric).unwrap_err();
        assert!(error.to_string().contains("DOT file"), "{error}");
        let _ = std::fs::remove_dir_all(&root);
    }

    /// The driver's ε is vakint's (`vakint::ε`, as in every value), and a
    /// value with another ε-like or D symbol is flagged.
    #[test]
    fn the_epsilon_is_vakints_and_foreign_symbols_are_flagged() {
        let e = epsilon_symbol();
        assert_eq!(e.get_name(), "vakint::ε");
        assert_eq!(e, symbol!("vakint::ε"));
        let good = symbolica::try_parse!(
            "refpath::CA*vakint::ε^-1+refpath::nf*log(vakint::mursq)+𝜋^2*vakint::ε"
        )
        .unwrap();
        assert!(foreign_symbols(&good, e).is_empty());
        let bad = symbolica::try_parse!("refpath::CA/(3-2*vakint::vakint::ε)").unwrap();
        assert_eq!(
            foreign_symbols(&bad, e),
            vec!["vakint::vakint::ε".to_string()]
        );
    }

    /// The digest's canonical form = the reader's Python definition
    /// (`json.dumps(config, sort_keys=True, separators=(',', ':'),
    /// ensure_ascii=False)`): build 8's config.json (keys here in another
    /// order) hashes to the digest the reader recomputed for it.
    #[test]
    fn the_canonical_form_reproduces_build_8s_digest() {
        let text = r#"{"binary_sha256":"34cf5dbf6c1e97766a9459d0ff3e17b38e8c89cce5db10922a2ba8db67a8841a","colour":"invariants (refpath symbols CA CF TF NA nf; higher invariants named)","contract_only":false,"dimension":"vakint_gate::D","environment":{"RUSTRED_P_LIBRARY_1_4":"/shared/vhirschi/t3code/vakint_dev/p_integral_dev/artefacts/m5/batch6/p-library-1-4","VAKINT_P_MASTER_CATALOGUE":"/shared/vhirschi/t3code/vakint_dev/p_integral_dev/artefacts/vk5l/ggppss_catalogue.v3.json","VAKINT_P_ZERO_PRUNE":"early","VAKINT_REQUIRE_FORM5":"1"},"eps_window":"eps^-L .. eps^-1","extra_terms":0,"loop_sign_in_overall":true,"master_catalogue_sha256":"b787f265e70ff4eb82e72d3655ac14c56c83c2849b4736f96df03f3f4065161d","model_sha256":"3e73d91ce3889fce13785bf2f9c30651b8c9e71ac9821d593471e3563556fbd9","modes":["retain","symbolic"],"overall_factor":"DOT attribute, required and checked (OverallFactorSource::Dot)","p_library":"/common/build-cache/vhirschi/p_integral_dev/vk5l-wp6/pkg-1-5/s12349d-pacc/lib-main","p_library_json_sha256":"600bb00a0bf0392920cce83a5f6b1726e9061cf090cf412ffed97db301a73e9e","projection":"metric","schema":"refpath-batch.config.v1","shared_engine":false,"zero_prune":"early"}"#;
        let mut value: Value = serde_json::from_str(text).unwrap();
        // reorder the keys: the canonical form must not depend on it
        let object = value.as_object_mut().unwrap();
        let schema = object.remove("schema").unwrap();
        object.insert("schema".into(), schema);
        assert_eq!(canonical_json(&value), text);
        assert_eq!(
            sha256_hex(canonical_json(&value).as_bytes()),
            "1f9be719271ef95bc0aa7f9982efcd2d8fa880f2717e5b61b1cd1d952f4c9c90"
        );
    }

    /// The content digest is sha256 over sorted "path\tsize\tsha256\n" rows of
    /// every file (written out here independently of the implementation).
    #[test]
    fn the_content_digest_follows_the_readers_definition() {
        let root =
            std::env::temp_dir().join(format!("refpath_batch_digest_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        std::fs::create_dir_all(root.join("b")).unwrap();
        std::fs::write(root.join("a.json"), "{}").unwrap();
        std::fs::write(root.join("b/c.bin"), [1u8, 2, 3]).unwrap();
        let rows = format!(
            "a.json\t2\t{}\nb/c.bin\t3\t{}\n",
            sha256_hex(b"{}"),
            sha256_hex(&[1u8, 2, 3])
        );
        let (digest, files, bytes) = content_digest(&root).unwrap();
        assert_eq!((digest, files, bytes), (sha256_hex(rows.as_bytes()), 2, 5));
        #[cfg(unix)]
        {
            std::os::unix::fs::symlink(root.join("b"), root.join("link")).unwrap();
            assert!(content_digest(&root).is_err());
        }
        let _ = std::fs::remove_dir_all(&root);
    }

    /// The Laurent support reader: eps powers of every term, logs eps-free.
    #[test]
    fn laurent_support_reads_every_order() {
        let e = symbol!("vakint::ε");
        let value =
            symbolica::try_parse!("vakint::ε^-2*3+vakint::ε^-1*log(2)+5+7*vakint::ε^2*vakint::x")
                .unwrap();
        assert_eq!(laurent_support(&value, e), vec![-2, -1, 0, 2]);
        assert_eq!(laurent_support(&Atom::Zero, e), Vec::<i64>::new());
    }
}
