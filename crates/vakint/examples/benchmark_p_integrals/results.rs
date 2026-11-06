//! The append-only results file of benchmark v1,
//! `docs/benchmark/data/p_results.json` (schema `vakint.p-benchmark.v1`): runs
//! with their protocol, each holding the raw measurements in the order they
//! were taken. Nothing recorded is ever rewritten; the report reads the
//! measurements of every run with the shas of the selected run, so that a
//! resumed run continues the interrupted one.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

use eyre::{Context, Result, bail, ensure};
use serde_json::{Value, json};

pub const SCHEMA: &str = "vakint.p-benchmark.v1";

/// What one measurement measures. A recorded key is never measured again for
/// the same shas and profile.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Key {
    pub cell: String,
    pub backend: String,
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

/// The measurements of this source and build: every run with the same
/// gammaloop and rustred shas and the same build profile.
fn of_source<'a>(
    document: &'a Value,
    gammaloop_sha: &'a str,
    rustred_sha: &'a str,
    build_profile: &'a str,
) -> impl Iterator<Item = &'a Value> {
    document["runs"]
        .as_array()
        .into_iter()
        .flatten()
        .filter(move |run| {
            let protocol = &run["protocol"];
            protocol["gammaloop_sha"] == gammaloop_sha
                && protocol["rustred_sha"] == rustred_sha
                && protocol["build_profile"] == build_profile
        })
        .flat_map(|run| run["measurements"].as_array().into_iter().flatten())
}

/// The keys already measured for this source and build.
pub fn recorded(
    document: &Value,
    gammaloop_sha: &str,
    rustred_sha: &str,
    build_profile: &str,
) -> Result<BTreeSet<Key>> {
    of_source(document, gammaloop_sha, rustred_sha, build_profile)
        .map(Key::of)
        .collect()
}

/// The result digest recorded per (cell, backend) for this source and build.
pub fn outputs(
    document: &Value,
    gammaloop_sha: &str,
    rustred_sha: &str,
    build_profile: &str,
) -> BTreeMap<(String, String), String> {
    of_source(document, gammaloop_sha, rustred_sha, build_profile)
        .filter_map(|measurement| {
            let cell = measurement["cell"].as_str()?.to_owned();
            let backend = measurement["backend"].as_str()?.to_owned();
            let digest = measurement["output_sha256"].as_str()?.to_owned();
            Some(((cell, backend), digest))
        })
        .collect()
}
