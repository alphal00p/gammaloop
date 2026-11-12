//! The append-only results file of the p benchmark,
//! `docs/benchmark/data/p_results.json` (schema `vakint.p-benchmark.v1`): runs
//! with their protocol, each holding the raw measurements in the order they
//! were taken. Nothing recorded is ever rewritten; the report reads the
//! measurements of every run of the same source as the selected run, so that
//! a resumed run continues the interrupted one, and those of the source of
//! the cited run (`cited_run`) for what the selected source does not measure.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use eyre::{Context, Result, bail, ensure};
use serde_json::{Value, json};
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    symbol,
};

pub const SCHEMA: &str = "vakint.p-benchmark.v1";
pub const EXACT_SCHEMA: &str = "vakint.p-exact.v1";

/// What is measured: the sources and the build. `gammaloop_sha` is the last
/// commit outside `crates/vakint/docs/benchmark`, `rustred_sha` the commit
/// pinned in `Cargo.lock`, `p_library` the identity of the RustRed p library
/// the children reduce with (protocol `p_library.identity`: the built-in one,
/// or the one an external-data root lists). A dirty tree is another source
/// than its clean sha.
#[derive(Clone, Copy, Debug)]
pub struct Source<'a> {
    pub gammaloop_sha: &'a str,
    pub rustred_sha: &'a str,
    pub build_profile: &'a str,
    pub source_dirty: bool,
    pub p_library: &'a str,
}

/// The arithmetic of a Forcer measurement, which has none to choose.
pub const NO_ARITHMETIC: &str = "none";

/// What one measurement measures. A recorded key is never measured again for
/// the same source, unless a capped record leaves room for a longer cap.
/// `arithmetic` is RustRed's Ordered weight arithmetic by mode name (`exact`,
/// `modular_polynomial`, `modular_sampled`), [`NO_ARITHMETIC`] for Forcer.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Key {
    pub cell: String,
    pub backend: String,
    pub arithmetic: String,
    pub threads: usize,
    pub state: String,
    pub repetition: usize,
}

impl Key {
    fn of(measurement: &Value) -> Result<Self> {
        let text = |name: &str| {
            measurement[name]
                .as_str()
                .map(str::to_owned)
                .ok_or_else(|| eyre::eyre!("measurement without {name}: {measurement}"))
        };
        let number = |name: &str| {
            measurement[name]
                .as_u64()
                .map(|value| value as usize)
                .ok_or_else(|| eyre::eyre!("measurement without {name}: {measurement}"))
        };
        Ok(Self {
            cell: text("cell")?,
            backend: text("backend")?,
            arithmetic: text("arithmetic")?,
            threads: number("threads")?,
            state: text("state")?,
            repetition: number("repetition")?,
        })
    }
}

/// The recorded document at `path`, or an empty one when there is no file.
pub fn read(path: &Path) -> Result<Value> {
    if !path.exists() {
        return Ok(json!({"schema": SCHEMA, "selected_run": null, "runs": []}));
    }
    let document: Value = serde_json::from_reader(File::open(path)?)
        .with_context(|| format!("reading {}", path.display()))?;
    ensure!(
        document["schema"] == SCHEMA && document["runs"].is_array(),
        "{} is not a {SCHEMA} document",
        path.display()
    );
    Ok(document)
}

fn write(path: &Path, document: &Value) -> Result<()> {
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent)?;
    let mut temporary = tempfile::NamedTempFile::new_in(parent)?;
    serde_json::to_writer_pretty(&mut temporary, document)?;
    temporary.write_all(b"\n")?;
    temporary.persist(path)?;
    Ok(())
}

/// Record a new run with its `protocol` and, on request, select it.
pub fn begin_run(path: &Path, id: &str, protocol: &Value, select: bool) -> Result<()> {
    let mut document = read(path)?;
    let runs = document["runs"].as_array_mut().expect("checked by read");
    if runs.iter().any(|run| run["id"] == id) {
        bail!("run {id} is already recorded");
    }
    runs.push(json!({"id": id, "protocol": protocol, "measurements": []}));
    if select {
        document["selected_run"] = json!(id);
    }
    write(path, &document)
}

/// Cite the recorded run `id` in the report: the measurements of its source
/// fill the cells the selected source does not measure.
pub fn cite_run(path: &Path, id: &str) -> Result<()> {
    let mut document = read(path)?;
    let runs = document["runs"].as_array().expect("checked by read");
    ensure!(
        runs.iter().any(|run| run["id"] == id),
        "no run {id} to cite in {}",
        path.display()
    );
    document["cited_run"] = json!(id);
    write(path, &document)
}

/// Append one measurement to the run `id` and save at once: an interrupted
/// harness loses at most the measurement in flight.
pub fn append(path: &Path, id: &str, measurement: Value) -> Result<()> {
    let key = Key::of(&measurement)?;
    let mut document = read(path)?;
    let runs = document["runs"].as_array_mut().expect("checked by read");
    let Some(run) = runs.iter_mut().find(|run| run["id"] == id) else {
        bail!("run {id} is not recorded");
    };
    let measurements = run["measurements"]
        .as_array_mut()
        .ok_or_else(|| eyre::eyre!("run {id} has no measurement list"))?;
    ensure!(
        !measurements
            .iter()
            .any(|recorded| Key::of(recorded).is_ok_and(|recorded| recorded == key)),
        "{key:?} is already recorded in run {id}"
    );
    measurements.push(measurement);
    write(path, &document)
}

/// The measurements of `source`.
fn of_source<'a>(document: &'a Value, source: Source<'a>) -> impl Iterator<Item = &'a Value> {
    document["runs"]
        .as_array()
        .into_iter()
        .flatten()
        .filter(move |run| {
            let protocol = &run["protocol"];
            protocol["gammaloop_sha"] == source.gammaloop_sha
                && protocol["rustred_sha"] == source.rustred_sha
                && protocol["p_library"]["identity"] == source.p_library
                && protocol["build_profile"] == source.build_profile
                && protocol["source_dirty"] == source.source_dirty
        })
        .flat_map(|run| run["measurements"].as_array().into_iter().flatten())
}

/// Whether a measurement needs no new attempt under the cap `cap_seconds`: a
/// measurement that hit a cap, or was skipped for that reason, proves nothing
/// about a longer one.
fn is_final(measurement: &Value, cap_seconds: Option<f64>) -> bool {
    measurement["status"] == "ok"
        || cap_seconds
            .zip(measurement["cap_seconds"].as_f64())
            .is_some_and(|(cap, recorded)| recorded >= cap)
}

/// The keys already measured for `source` that a run under `cap_seconds` skips.
pub fn recorded(
    document: &Value,
    source: Source,
    cap_seconds: Option<f64>,
) -> Result<BTreeSet<Key>> {
    of_source(document, source)
        .filter(|measurement| is_final(measurement, cap_seconds))
        .map(Key::of)
        .collect()
}

/// The result digest recorded per (cell, backend) for `source`.
pub fn outputs(document: &Value, source: Source) -> BTreeMap<(String, String), String> {
    of_source(document, source)
        .filter_map(|measurement| {
            let cell = measurement["cell"].as_str()?.to_owned();
            let backend = measurement["backend"].as_str()?.to_owned();
            let digest = measurement["output_sha256"].as_str()?.to_owned();
            Some(((cell, backend), digest))
        })
        .collect()
}

/// A RustRed child must have reduced with the library of the protocol
/// (`p_library`, see [`Source`]): every p program its cache counters list
/// (`p-library-<a>-<b>`) is that library's directory, loaded from its origin
/// (`built-in` or `external-data`), and there is one.
pub fn check_p_program(report: &Value, p_library: &Value) -> Result<()> {
    let programs = report["cache"]["programs"]
        .as_array()
        .into_iter()
        .flatten()
        .filter(|program| {
            program["program"]
                .as_str()
                .is_some_and(|name| name.starts_with("p-library-"))
        })
        .collect::<Vec<_>>();
    ensure!(
        !programs.is_empty(),
        "the child reports no p library program"
    );
    for program in programs {
        ensure!(
            program["program"] == p_library["directory"]
                && program["data_version"] == p_library["origin"],
            "the child reduced with {} ({}), the protocol names {} ({})",
            program["program"],
            program["data_version"],
            p_library["directory"],
            p_library["origin"]
        );
    }
    Ok(())
}

/// An exact-agreement file marks what was measured only for the same sources;
/// the stub of the repository marks nothing.
pub fn check_exact(exact: &Value, source: Source) -> Result<()> {
    ensure!(
        exact["schema"] == EXACT_SCHEMA && exact["cells"].is_array(),
        "not a {EXACT_SCHEMA} document"
    );
    if exact["stub"] != true {
        for (name, sha) in [
            ("gammaloop_sha", source.gammaloop_sha),
            ("rustred_sha", source.rustred_sha),
        ] {
            ensure!(
                exact[name] == sha,
                "the exact-agreement file has {name} {}, this source {sha}",
                exact[name]
            );
        }
    }
    Ok(())
}

/// Whether this binary may record measurements of `source`. `stamp` is the
/// source it was built from (`just vakint-p-benchmark-build`, an optional
/// compile-time variable): the tree alone cannot tell an old binary from a
/// current one. A release measurement is meant for the report, so it needs the
/// stamp, a clean tree and the real measuring child.
pub fn check_build(source: Source, stamp: Option<&str>, stub_child: bool) -> Result<()> {
    if source.build_profile == "release" {
        ensure!(
            !source.source_dirty,
            "a release measurement needs a clean tree: commit, or keep the change under \
            crates/vakint/docs/benchmark"
        );
        ensure!(
            !stub_child,
            "--child-program does not measure vakint: refused for a release measurement"
        );
        ensure!(
            stamp.is_some(),
            "a release measurement needs a stamped binary: build it with \
            `just vakint-p-benchmark-build`"
        );
    }
    if let Some(stamp) = stamp {
        let dirty = if source.source_dirty { "-dirty" } else { "" };
        let tree = format!("{}{dirty}", source.gammaloop_sha);
        ensure!(
            stamp == tree,
            "this binary is built from {stamp}, the tree is at {tree}: rebuild with \
            `just vakint-p-benchmark-build`"
        );
    }
    Ok(())
}

/// The seconds of the reduction sub-line. A capture without a math span
/// reports zero, which is no measurement.
pub fn reduction_seconds(complete: bool, intervals: usize, math: Option<f64>) -> Option<f64> {
    math.filter(|_| complete && intervals > 0)
}

/// The integrals of an input: its expanded terms up to a coefficient free of
/// loop momenta (the projector's powers of epsilon and of the external scale).
pub fn integral_count(input: &Atom) -> usize {
    let momentum = symbol!("vakint::k");
    let expanded = input.expand();
    let terms = match expanded.as_view() {
        AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
        term => vec![term],
    };
    let integrals = terms.into_iter().map(|term| {
        let factors = match term {
            AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
            factor => vec![factor],
        };
        let mut integral = factors
            .into_iter()
            .filter(|factor| factor.contains_symbol(momentum))
            .map(|factor| factor.to_canonical_string())
            .collect::<Vec<_>>();
        integral.sort();
        integral
    });
    integrals.collect::<BTreeSet<_>>().len()
}
