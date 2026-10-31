//! Opt-in capture of the native RustRed reduction requests Vakint issues.
//!
//! Setting `VAKINT_DUMP_ORDERED_REQUEST=<directory>` writes, for every native
//! reduction request (one per compatible mass/spectator channel):
//!
//! * `<stem>.request.json` before the request runs: loop count, the program
//!   identity needed to load the same native program from the external data,
//!   strategy (always Ordered), weight arithmetic, master mode, channel label, the exact
//!   input keys and
//!   coefficients (Symbolica canonical strings of rational functions in the
//!   family dimension) and the caller-side nonzero conditions;
//! * `<stem>.phases.jsonl` while it runs: the native ordered counters at every
//!   phase end, flushed per line, when no caller observer is attached;
//! * `<stem>.result.json` afterwards: BLAKE3 digests of the exact native output
//!   (terms, and the guard set) plus the final counters.
//!
//! Both also record the request's cumulative Ordered step-cache hits and
//! misses (`step_cache`), apart from the structural counters (`progress`),
//! which do not depend on the cache.
//!
//! `index.jsonl` in the same directory lists every written file. With
//! `VAKINT_DUMP_ORDERED_REQUEST_ONLY=1` the native reductions are skipped after
//! their requests are written and the evaluation fails with an explicit error.
//! RustRed's `ordered_replay` example replays these files without Vakint.
//! Without the variables nothing happens beyond one cached environment lookup
//! per native request.
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::{Duration, Instant};

use rustred::algebra::Coefficient;
use rustred::family::IntegralKey;
use rustred::solver::{
    CandidateOrderedReductionProgress, CandidateRoutedApplicationEnvelope,
    CandidateRoutedApplicationObserver, CandidateRoutedObservationOptions,
    CandidateRoutedObserverControl,
};
use serde_json::{Value, json};
use symbolica::atom::{Atom, AtomCore};

use super::RustRedEvaluationError;

pub(super) const DIRECTORY_VARIABLE: &str = "VAKINT_DUMP_ORDERED_REQUEST";
pub(super) const ONLY_VARIABLE: &str = "VAKINT_DUMP_ORDERED_REQUEST_ONLY";
const REQUEST_SCHEMA: &str = "vakint.ordered-native-request.v1";
const PHASE_SCHEMA: &str = "vakint.ordered-native-phase.v1";
const RESULT_SCHEMA: &str = "vakint.ordered-native-result.v1";
/// Identifies the canonical byte layout hashed by [`output_digest`].
pub(super) const OUTPUT_DIGEST_SCHEME: &str = "vakint.ordered-native-output.v1";

struct DumpSettings {
    directory: PathBuf,
    only: bool,
}

fn settings() -> Option<&'static DumpSettings> {
    static SETTINGS: OnceLock<Option<DumpSettings>> = OnceLock::new();
    SETTINGS
        .get_or_init(|| {
            let directory =
                std::env::var_os(DIRECTORY_VARIABLE).filter(|value| !value.is_empty())?;
            let only = std::env::var(ONLY_VARIABLE).is_ok_and(|value| value == "1");
            Some(DumpSettings {
                directory: directory.into(),
                only,
            })
        })
        .as_ref()
}

fn failure(detail: impl std::fmt::Display) -> RustRedEvaluationError {
    RustRedEvaluationError::Reduction {
        detail: format!("ordered request dump: {detail}"),
    }
}

/// Native program identity: enough to load the same program from external data.
pub(super) enum DumpProgram<'a> {
    /// 1-3L sealed `ClosedArtifact`.
    ClosedArtifact {
        family_fingerprint: &'a str,
        artifact_sha256: [u8; 32],
        data: &'a super::data::LoopData,
    },
    /// 4L shipped candidate program.
    CandidateProgram {
        label: &'a str,
        index: usize,
        family_fingerprint: &'a str,
        data: &'a super::data::LoopData,
    },
    /// 5L saved Route/Apply application.
    SavedApplication {
        path: &'a Path,
        family_fingerprint: &'a str,
        options: &'a rustred_app::SavedApplicationOptions,
    },
}

impl DumpProgram<'_> {
    fn describe(&self) -> Value {
        let files = |data: &super::data::LoopData, names: &[String]| -> Value {
            names
                .iter()
                .map(|name| {
                    let digest = data.manifest.files.get(name);
                    let path = data.path(name).ok().flatten();
                    (
                        name.clone(),
                        json!({
                            "sha256": digest.map(|digest| &digest.sha256),
                            "bytes": digest.map(|digest| digest.bytes),
                            "path": path,
                        }),
                    )
                })
                .collect::<serde_json::Map<_, _>>()
                .into()
        };
        match self {
            Self::ClosedArtifact {
                family_fingerprint,
                artifact_sha256,
                data,
            } => json!({
                "kind": "closed-artifact",
                "family_fingerprint": family_fingerprint,
                "artifact_sha256": hex(artifact_sha256),
                "data_version": data.manifest.version,
                "files": files(data, &["ibps/rules.rrbin".into()]),
            }),
            Self::CandidateProgram {
                label,
                index,
                family_fingerprint,
                data,
            } => {
                let name = label.trim_start_matches("four-loop-");
                json!({
                    "kind": "candidate-program",
                    "label": label,
                    "index": index,
                    "family_fingerprint": family_fingerprint,
                    "data_version": data.manifest.version,
                    "files": files(data, &[
                        format!("ibps/{name}.csv"),
                        format!("ibps/{name}.candidates.rrbin.gz"),
                        format!("ibps/{name}.rrcat.bin.gz"),
                        format!("ibps/{name}.rrnorm.bin"),
                    ]),
                })
            }
            Self::SavedApplication {
                path,
                family_fingerprint,
                options,
            } => {
                let replayable = rustred_app::SavedApplicationOptions {
                    reduction_strategy: rustred::reduction::ReductionStrategy::Ordered,
                    unresolved_policy: rustred::solver::UnresolvedPolicy::Preserve,
                    ..Default::default()
                };
                json!({
                    "kind": "saved-application",
                    "path": path,
                    "family_fingerprint": family_fingerprint,
                    "unresolved_policy": format!("{:?}", options.unresolved_policy),
                    "entry_coverage_policy": format!("{:?}", options.entry_coverage_policy),
                    "preparation_workers": options.preparation_workers,
                    "options_are_defaults_with_preserve": **options == replayable,
                    "operational_options": format!("{options:?}"),
                })
            }
        }
    }
}

/// Caller label of one native request.
pub(super) struct DumpChannel<'a> {
    pub(super) index: usize,
    pub(super) total: usize,
    pub(super) mass_squared: &'a Atom,
    pub(super) spectator: &'a Atom,
    pub(super) exponent: i64,
}

pub(super) struct DumpRequest<'a> {
    pub(super) loop_count: usize,
    pub(super) program: DumpProgram<'a>,
    /// Recorded for information: every arithmetic yields the same exact
    /// output, but modular weights report a different monomial census.
    pub(super) weight_arithmetic: crate::RustRedOrderedWeightArithmetic,
    pub(super) master_mode: Option<crate::RustRedMasterMode>,
    pub(super) channel: DumpChannel<'a>,
    pub(super) limits: rustred::reduction::ReductionLimits,
    pub(super) coefficient_variables: Vec<String>,
    pub(super) inputs: &'a [(IntegralKey, Coefficient)],
    /// Caller-side conditions; RustRed receives only `inputs`.
    pub(super) caller_conditions: Vec<String>,
}

/// One written request. Finish it with the exact native output.
pub(super) struct RequestDump {
    directory: PathBuf,
    stem: String,
    only: bool,
    started: Instant,
    phases: Option<PhaseLog>,
    last_progress: Option<CandidateOrderedReductionProgress>,
}

struct PhaseLog {
    writer: BufWriter<File>,
    recorded: usize,
    error: Option<std::io::Error>,
}

impl RequestDump {
    /// Write the request if dumping is enabled; zero-cost otherwise.
    pub(super) fn begin(
        request: DumpRequest<'_>,
    ) -> Result<Option<RequestDump>, RustRedEvaluationError> {
        settings()
            .map(|settings| Self::write(settings, request))
            .transpose()
    }

    fn write(
        settings: &DumpSettings,
        request: DumpRequest<'_>,
    ) -> Result<RequestDump, RustRedEvaluationError> {
        static SEQUENCE: AtomicUsize = AtomicUsize::new(0);
        let sequence = SEQUENCE.fetch_add(1, Ordering::Relaxed);
        std::fs::create_dir_all(&settings.directory).map_err(failure)?;
        let stem = format!(
            "L{}-{}-{sequence:04}",
            request.loop_count,
            std::process::id()
        );
        let mut hasher = blake3::Hasher::new();
        let inputs = request
            .inputs
            .iter()
            .map(|(key, coefficient)| {
                let coefficient = canonical(coefficient);
                hash_key(&mut hasher, key);
                hasher.update(coefficient.as_bytes());
                hasher.update(b"\n");
                json!({"powers": key.powers(), "coefficient": coefficient})
            })
            .collect::<Vec<_>>();
        let document = json!({
            "schema": REQUEST_SCHEMA,
            "sequence": sequence,
            "process_id": std::process::id(),
            "loop_count": request.loop_count,
            "program": request.program.describe(),
            // Vakint reduces only with the ordered scheduler; replay expects it.
            "reduction_strategy": "Ordered",
            "ordered_weight_arithmetic": format!("{:?}", request.weight_arithmetic),
            "master_mode": request.master_mode.map(|mode| format!("{mode:?}")),
            "channel": {
                "index": request.channel.index,
                "total": request.channel.total,
                "mass_squared": request.channel.mass_squared.to_canonical_string(),
                "spectator": request.channel.spectator.to_canonical_string(),
                "mass_exponent": request.channel.exponent,
            },
            "reduction_limits_are_default":
                request.limits == rustred::reduction::ReductionLimits::default(),
            "reduction_limits": format!("{:?}", request.limits),
            "coefficient_variables": request.coefficient_variables,
            "input_terms": inputs.len(),
            "input_digest": hasher.finalize().to_hex().to_string(),
            "inputs": inputs,
            "caller_nonzero_conditions": request.caller_conditions,
            "skipped_by_request_only": settings.only,
        });
        let file = format!("{stem}.request.json");
        write_json(&settings.directory.join(&file), &document)?;
        append_index(
            &settings.directory,
            &json!({"event": "request", "file": file, "loop_count": request.loop_count,
                "input_terms": request.inputs.len(), "channel": request.channel.index,
                "channels": request.channel.total}),
        )?;
        log::info!(
            "ordered request dump: wrote {} ({} input terms)",
            settings.directory.join(&file).display(),
            request.inputs.len()
        );
        Ok(RequestDump {
            directory: settings.directory.clone(),
            stem,
            only: settings.only,
            started: Instant::now(),
            phases: None,
            last_progress: None,
        })
    }

    /// The native reduction must be skipped after this request was written.
    pub(super) fn only(&self) -> bool {
        self.only
    }

    fn record_phase(&mut self, progress: &CandidateOrderedReductionProgress) {
        let elapsed = self.started.elapsed();
        self.last_progress = Some(*progress);
        let Some(log) = self.phases.as_mut() else {
            return;
        };
        if log.error.is_some() || progress.completed_phases <= log.recorded {
            return;
        }
        let record = json!({
            "schema": PHASE_SCHEMA,
            "phase": progress.completed_phases - 1,
            // Earlier ends were not delivered separately; never expected.
            "missed_phase_ends": progress.completed_phases - 1 - log.recorded,
            "elapsed_seconds": elapsed.as_secs_f64(),
            "progress": progress_json(progress),
            "step_cache": step_cache_json(progress),
        });
        log.recorded = progress.completed_phases;
        if let Err(error) = writeln!(log.writer, "{record}").and_then(|()| log.writer.flush()) {
            log.error = Some(error);
        }
    }

    /// Record the exact output digests and close the dump. `progress` is the
    /// caller observer's final census; the dump's own recorder supplies it otherwise.
    pub(super) fn finish<'k>(
        self,
        terms: impl IntoIterator<Item = (&'k str, &'k IntegralKey, &'k Coefficient)>,
        nonzero_conditions: &[Coefficient],
        progress: Option<CandidateOrderedReductionProgress>,
    ) -> Result<(), RustRedEvaluationError> {
        let seconds = self.started.elapsed().as_secs_f64();
        if let Some(log) = self.phases
            && let Some(error) = log.error
        {
            return Err(failure(error));
        }
        let digest = output_digest(terms, nonzero_conditions);
        let file = format!("{}.result.json", self.stem);
        write_json(
            &self.directory.join(&file),
            &json!({
                "schema": RESULT_SCHEMA,
                "request": format!("{}.request.json", self.stem),
                "wall_seconds": seconds,
                "output": digest,
                "progress": progress.or(self.last_progress).as_ref().map(progress_json),
                "step_cache": progress.or(self.last_progress).as_ref().map(step_cache_json),
            }),
        )?;
        append_index(
            &self.directory,
            &json!({"event": "result", "file": file, "wall_seconds": seconds,
                "terms_blake3": digest["terms_blake3"]}),
        )
    }
}

/// Reduce with the caller's observer, or record phase ends with a minimal
/// observer while dumping a request that has none. Observation never changes
/// the native mathematics; it only reads safe-boundary counters.
pub(super) fn observe_native<T>(
    dump: Option<&mut RequestDump>,
    observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    reduce: impl FnOnce(Option<&mut CandidateRoutedApplicationObserver<'_>>) -> T,
) -> Result<T, RustRedEvaluationError> {
    let Some(dump) = dump.filter(|_| observation.is_none()) else {
        return Ok(reduce(observation));
    };
    let file = format!("{}.phases.jsonl", dump.stem);
    dump.phases = Some(PhaseLog {
        writer: BufWriter::new(File::create(dump.directory.join(&file)).map_err(failure)?),
        recorded: 0,
        error: None,
    });
    append_index(&dump.directory, &json!({"event": "phases", "file": file}))?;
    let mut callback = |envelope: &CandidateRoutedApplicationEnvelope<'_>| {
        if let Some(progress) = &envelope.ordered {
            dump.record_phase(progress);
        }
        Ok(CandidateRoutedObserverControl::Continue)
    };
    let mut observer = CandidateRoutedApplicationObserver::new(
        CandidateRoutedObservationOptions {
            // Phase ends are always delivered; nothing else is sampled.
            interval: Duration::MAX,
            max_frontier_frames: 0,
            max_path_nodes: 0,
        },
        &mut callback,
    );
    Ok(reduce(Some(&mut observer)))
}

fn canonical(coefficient: &Coefficient) -> String {
    coefficient.to_expression().to_canonical_string()
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn hash_key(hasher: &mut blake3::Hasher, key: &IntegralKey) {
    let powers = key
        .powers()
        .iter()
        .map(i64::to_string)
        .collect::<Vec<_>>()
        .join(",");
    hasher.update(powers.as_bytes());
    hasher.update(b"\n");
}

/// Digest of an exact native output, in the layout `OUTPUT_DIGEST_SCHEME`:
/// terms in native order as `family\npowers\ncoefficient\n`; guards as the
/// sorted set of canonical strings, each followed by `\n`. A valid guard
/// superset changes only the guard digest.
pub(super) fn output_digest<'k>(
    terms: impl IntoIterator<Item = (&'k str, &'k IntegralKey, &'k Coefficient)>,
    nonzero_conditions: &[Coefficient],
) -> Value {
    let mut hasher = blake3::Hasher::new();
    let mut count = 0usize;
    for (family, key, coefficient) in terms {
        hasher.update(family.as_bytes());
        hasher.update(b"\n");
        hash_key(&mut hasher, key);
        hasher.update(canonical(coefficient).as_bytes());
        hasher.update(b"\n");
        count += 1;
    }
    let guards = nonzero_conditions
        .iter()
        .map(canonical)
        .collect::<BTreeSet<_>>();
    let mut guard_hasher = blake3::Hasher::new();
    for guard in &guards {
        guard_hasher.update(guard.as_bytes());
        guard_hasher.update(b"\n");
    }
    json!({
        "scheme": OUTPUT_DIGEST_SCHEME,
        "terms": count,
        "terms_blake3": hasher.finalize().to_hex().to_string(),
        "nonzero_conditions": nonzero_conditions.len(),
        "distinct_nonzero_conditions": guards.len(),
        "nonzero_conditions_blake3": guard_hasher.finalize().to_hex().to_string(),
    })
}

/// Cumulative Ordered step-cache counters of the request, kept apart from
/// the structural counters (which do not depend on the cache).
fn step_cache_json(progress: &CandidateOrderedReductionProgress) -> Value {
    json!({
        "hits": progress.step_cache_hits,
        "misses": progress.step_cache_misses,
    })
}

fn progress_json(progress: &CandidateOrderedReductionProgress) -> Value {
    json!({
        "completed_phases": progress.completed_phases,
        "total_phases": progress.total_phases,
        "input_keys": progress.input_keys,
        "processed_keys": progress.processed_keys,
        "emitted_terms": progress.emitted_terms,
        "merged_terms": progress.merged_terms,
        "cancelled_terms": progress.cancelled_terms,
        "terminal_merges": progress.terminal_merges,
        "terminal_cancellations": progress.terminal_cancellations,
        "frontier_keys": progress.frontier_keys,
        "peak_frontier_keys": progress.peak_frontier_keys,
        "terminal_keys": progress.terminal_keys,
        "frontier_coefficient_monomials": progress.frontier_coefficient_monomials,
        "peak_frontier_coefficient_monomials": progress.peak_frontier_coefficient_monomials,
        "terminal_coefficient_monomials": progress.terminal_coefficient_monomials,
        "peak_terminal_coefficient_monomials": progress.peak_terminal_coefficient_monomials,
    })
}

fn write_json(path: &Path, document: &Value) -> Result<(), RustRedEvaluationError> {
    let mut writer = BufWriter::new(File::create(path).map_err(failure)?);
    serde_json::to_writer(&mut writer, document).map_err(failure)?;
    writer.write_all(b"\n").map_err(failure)?;
    writer.flush().map_err(failure)
}

fn append_index(directory: &Path, record: &Value) -> Result<(), RustRedEvaluationError> {
    let mut file = std::fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(directory.join("index.jsonl"))
        .map_err(failure)?;
    writeln!(file, "{record}").map_err(failure)
}

/// Master mode of the first RustRed backend in the evaluation order.
pub(super) fn master_mode(settings: &crate::VakintSettings) -> Option<crate::RustRedMasterMode> {
    settings
        .evaluation_order
        .0
        .iter()
        .find_map(|method| match method {
            crate::EvaluationMethod::RustRed(options) => Some(options.master_mode),
            _ => None,
        })
}

/// Error returned after all requests of a batch were written and skipped.
pub(super) fn skipped(count: usize) -> RustRedEvaluationError {
    failure(format!(
        "{count} native request(s) written and skipped because {ONLY_VARIABLE}=1"
    ))
}

/// Terms of a single-family native output, for [`RequestDump::finish`].
pub(super) fn family_terms<'k>(
    family: &'k str,
    terms: &'k BTreeMap<IntegralKey, Coefficient>,
) -> impl Iterator<Item = (&'k str, &'k IntegralKey, &'k Coefficient)> {
    terms
        .iter()
        .map(move |(key, coefficient)| (family, key, coefficient))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rustred::algebra::CoefficientContext;

    #[test]
    fn request_and_result_record_exact_canonical_inputs_and_output_digests() {
        let directory = tempfile::tempdir().unwrap();
        let settings = DumpSettings {
            directory: directory.path().to_owned(),
            only: true,
        };
        let context = CoefficientContext::try_new(["d"]).unwrap();
        let d = context.parameter("d").unwrap();
        let weight = context
            .try_div(&context.one(), &d, Default::default())
            .unwrap();
        let key = IntegralKey::try_new([2, -1]).unwrap();
        let inputs = vec![(key.clone(), weight.clone())];
        let options = rustred_app::SavedApplicationOptions {
            reduction_strategy: rustred::reduction::ReductionStrategy::Ordered,
            unresolved_policy: rustred::solver::UnresolvedPolicy::Preserve,
            ..Default::default()
        };
        let mass = Atom::num(1);
        let dump = RequestDump::write(
            &settings,
            DumpRequest {
                loop_count: 5,
                program: DumpProgram::SavedApplication {
                    path: Path::new("/data/5L/ibps"),
                    family_fingerprint: "family",
                    options: &options,
                },
                weight_arithmetic: options.ordered_weight_arithmetic,
                master_mode: Some(crate::RustRedMasterMode::Retain),
                channel: DumpChannel {
                    index: 0,
                    total: 1,
                    mass_squared: &mass,
                    spectator: &mass,
                    exponent: -3,
                },
                limits: Default::default(),
                coefficient_variables: vec!["d".into()],
                inputs: &inputs,
                caller_conditions: vec![],
            },
        )
        .unwrap();
        assert!(dump.only());
        let stem = dump.stem.clone();
        let request: Value = serde_json::from_slice(
            &std::fs::read(directory.path().join(format!("{stem}.request.json"))).unwrap(),
        )
        .unwrap();
        assert_eq!(request["schema"], REQUEST_SCHEMA);
        assert_eq!(request["ordered_weight_arithmetic"], "Exact");
        assert_eq!(
            request["program"]["options_are_defaults_with_preserve"],
            true
        );
        assert_eq!(request["reduction_limits_are_default"], true);
        assert_eq!(request["inputs"][0]["powers"], json!([2, -1]));
        let text = request["inputs"][0]["coefficient"].as_str().unwrap();
        assert_eq!(text, canonical(&weight));
        let expected = blake3::hash(format!("2,-1\n{text}\n").as_bytes());
        assert_eq!(request["input_digest"], expected.to_hex().to_string());

        let terms = BTreeMap::from([(key, weight.clone())]);
        dump.finish(family_terms("family", &terms), &[d.clone(), d], None)
            .unwrap();
        let result: Value = serde_json::from_slice(
            &std::fs::read(directory.path().join(format!("{stem}.result.json"))).unwrap(),
        )
        .unwrap();
        let output = &result["output"];
        assert_eq!(output["scheme"], OUTPUT_DIGEST_SCHEME);
        assert_eq!(
            (
                output["terms"].as_u64(),
                output["nonzero_conditions"].as_u64()
            ),
            (Some(1), Some(2))
        );
        assert_eq!(output["distinct_nonzero_conditions"], 1);
        let terms = blake3::hash(format!("family\n2,-1\n{text}\n").as_bytes());
        assert_eq!(output["terms_blake3"], terms.to_hex().to_string());
        let index = std::fs::read_to_string(directory.path().join("index.jsonl")).unwrap();
        assert_eq!(index.lines().count(), 2);
    }
}
