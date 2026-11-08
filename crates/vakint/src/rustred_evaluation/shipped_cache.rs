//! Ordered step caches shipped with the external data.
//!
//! A loop pack may ship, for each of its reduction programs, the RustRed
//! Ordered step-cache directory that a cold run of a deterministic workload
//! wrote: `ordered-cache/<label>.<tables>.<stamp>/` with the RustRed manifest,
//! the shards it lists and optionally the tail table. The directory name is
//! the one an engine cache directory uses for the same program, rule tables
//! and binding stamp (see [`super::cache`]). `ordered-cache.json` (schema
//! [`METADATA_SCHEMA`]) describes the programs and the workload, and the loop
//! manifest lists every file with its sha256 like any other payload.
//!
//! Validity is decided by versions only, as for every Ordered cache. A
//! shipped cache is used iff it was recorded for the program's rule tables
//! under this build's binding stamp (RustRed's step-cache schema, its
//! step-semantics version and the binding's identity fields; a tail table must
//! carry this build's tail schema and semantics) and its files match the loop
//! manifest's sha256. Otherwise vakint prints one line on stderr and the
//! program runs without it (cold); a shipped cache never fails an evaluation
//! and never changes a result. Only such verdicts about the data are kept for
//! the engine: when the request's resource limits stop the verification
//! (reading the files checks them), that request fails with the resource
//! error and the next one verifies again.
//!
//! Without a cache directory (the default in-memory mode) an accepted cache
//! is opened in place: its shards are read lazily, new records stay in memory
//! and nothing is ever written into the data directory. With a read-write
//! cache directory its files are copied into the program's directory there
//! when that directory holds no cache yet; read-only mode never writes and
//! opens it in place instead.
//!
//! The maintenance example `ordered_cache_maintenance` prepares the shipped
//! caches ([`crate::Vakint::export_ordered_caches`]) and checks a data
//! directory against a build ([`crate::Vakint::shipped_ordered_caches`]).
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
};

use rustred::reduction::{
    ORDERED_STEP_CACHE_SCHEMA, ORDERED_STEP_SEMANTICS, ORDERED_TAIL_SCHEMA, ORDERED_TAIL_SEMANTICS,
};
use serde::{Deserialize, Serialize};

use super::{
    RustRedEvaluationError,
    cache::{OrderedCacheProgram, hex, ordered_cache_failure, program_directory_name},
    data::{LoopData, RustRedData},
};
use crate::{Vakint, VakintError};

/// Schema of a loop pack's `ordered-cache.json`.
pub const METADATA_SCHEMA: &str = "vakint.external-data.ordered-cache.v2";
/// Directory of the shipped caches inside a loop pack.
pub(crate) const PACK_DIRECTORY: &str = "ordered-cache";
/// Metadata file of the shipped caches inside a loop pack.
pub(crate) const METADATA_FILE: &str = "ordered-cache.json";
/// RustRed's step-cache manifest inside a cache directory.
pub(super) const RUSTRED_MANIFEST: &str = "ordered-step-cache.txt";
/// RustRed's tail table inside a cache directory.
const TAIL_FILE: &str = "ordered-tail.rrtail";
const SHARD_DIRECTORY: &str = "shards";

/// The loop-manifest section of shipped Ordered caches.
#[derive(Clone, Debug, PartialEq, Eq, Deserialize, Serialize)]
pub(crate) struct OrderedCacheSection {
    pub path: String,
    pub metadata: String,
}

/// One program cache of a pack's `ordered-cache.json`.
#[derive(Clone, Debug, PartialEq, Eq, Deserialize, Serialize)]
pub struct ShippedOrderedCacheProgram {
    /// Vakint program label, e.g. `4L-four-loop-h`.
    pub label: String,
    /// `ordered-cache/<label>.<tables>.<stamp>`, relative to the loop pack.
    pub directory: String,
    /// Hexadecimal identity of the program's rule tables.
    pub rule_tables: String,
    /// Hexadecimal RustRed binding stamp the records were written under.
    pub stamp: String,
    /// RustRed adapter (`closed-artifact`, `candidate-program`,
    /// `saved-application`).
    pub adapter: String,
    pub records: usize,
    pub shards: usize,
    /// Bytes of the shard files.
    pub bytes: u64,
    /// Bytes of the tail table (0: none shipped).
    pub tail_bytes: u64,
    /// Whether the directory holds every record (and the tail table) of the
    /// workload, or a trimmed selection.
    pub complete: bool,
    /// How a trimmed cache was selected from the complete one.
    #[serde(default)]
    pub trim: Option<serde_json::Value>,
}

/// A pack's `ordered-cache.json`.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ShippedOrderedCacheMetadata {
    pub schema: String,
    pub loop_count: usize,
    /// Producing vakint/RustRed revisions, versions and command (reported,
    /// never part of the decision).
    #[serde(default)]
    pub producer: serde_json::Value,
    /// The deterministic workload whose steps the caches record.
    #[serde(default)]
    pub workload: Vec<serde_json::Value>,
    pub programs: Vec<ShippedOrderedCacheProgram>,
}

/// How one program uses the cache its data ships.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum ShippedOrderedCacheStatus {
    /// Opened in place without a cache directory: replayed, never written.
    InPlace,
    /// Copied into the program's directory of the cache directory.
    Installed,
    /// The program's directory of the cache directory already held a cache,
    /// which is used (and kept) instead.
    AlreadyPresent,
    /// Refused (another version, other rule tables, corrupt or inconsistent
    /// files); the program runs without it.
    Rejected,
}

impl ShippedOrderedCacheStatus {
    /// The canonical text: `in-place`, `installed`, `already-present` or
    /// `rejected`.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::InPlace => "in-place",
            Self::Installed => "installed",
            Self::AlreadyPresent => "already-present",
            Self::Rejected => "rejected",
        }
    }
}

/// The shipped cache of one program and what became of it.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct ShippedOrderedCacheUse {
    pub status: ShippedOrderedCacheStatus,
    /// The shipped directory, `<L>L/ordered-cache/<label>.<tables>.<stamp>`
    /// inside the data.
    pub source: PathBuf,
    /// Records and bytes (shards and tail) declared by `ordered-cache.json`.
    pub records: usize,
    pub bytes: u64,
    pub complete: bool,
    /// Why the cache was rejected, or why it could not be installed.
    pub detail: Option<String>,
}

/// One loaded program, its cache identity and its shipped cache.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct ShippedOrderedCacheCheck {
    pub loop_count: usize,
    pub program: String,
    /// The program's cache directory name, in a pack's `ordered-cache/` and in
    /// an engine cache directory alike.
    pub directory: String,
    pub rule_tables: String,
    pub stamp: String,
    /// `None` when the data ships no cache for the program.
    pub shipped: Option<ShippedOrderedCacheUse>,
}

/// A shipped cache that passed every check: the files of `source` (paths
/// relative to it, manifest last) are the ones the loop manifest lists.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct AcceptedOrderedCache {
    pub source: PathBuf,
    pub files: Vec<String>,
    pub entry: ShippedOrderedCacheProgram,
}

impl AcceptedOrderedCache {
    pub(crate) fn report(
        &self,
        status: ShippedOrderedCacheStatus,
        detail: Option<String>,
    ) -> ShippedOrderedCacheUse {
        ShippedOrderedCacheUse {
            status,
            source: self.source.clone(),
            records: self.entry.records,
            bytes: self.entry.bytes + self.entry.tail_bytes,
            complete: self.entry.complete,
            detail,
        }
    }

    /// Copy the cache into `destination` (a program directory of an engine
    /// cache directory) unless it already holds a cache. Files are staged
    /// next to it and moved in, the RustRed manifest last, so an interrupted
    /// copy never leaves a partial cache behind a manifest.
    pub(crate) fn install(&self, destination: &Path) -> ShippedOrderedCacheUse {
        if destination.join(RUSTRED_MANIFEST).exists() {
            return self.report(ShippedOrderedCacheStatus::AlreadyPresent, None);
        }
        match self.copy_into(destination) {
            Ok(()) => self.report(ShippedOrderedCacheStatus::Installed, None),
            Err(error) => {
                let detail = format!("could not copy it into {}: {error}", destination.display());
                loud_rejection(&self.source, &self.entry.label, &detail);
                self.report(ShippedOrderedCacheStatus::Rejected, Some(detail))
            }
        }
    }

    fn copy_into(&self, destination: &Path) -> Result<(), String> {
        let io = |path: &Path, error: std::io::Error| format!("{}: {error}", path.display());
        let parent = destination
            .parent()
            .ok_or("the cache directory has no parent")?;
        fs::create_dir_all(parent).map_err(|error| io(parent, error))?;
        let staging = tempfile::Builder::new()
            .prefix(".shipped-ordered-cache-")
            .tempdir_in(parent)
            .map_err(|error| io(parent, error))?;
        for file in &self.files {
            let target = staging.path().join(file);
            if let Some(directory) = target.parent() {
                fs::create_dir_all(directory).map_err(|error| io(directory, error))?;
            }
            fs::copy(self.source.join(file), &target).map_err(|error| io(&target, error))?;
        }
        for file in &self.files {
            let target = destination.join(file);
            if let Some(directory) = target.parent() {
                fs::create_dir_all(directory).map_err(|error| io(directory, error))?;
            }
            fs::rename(staging.path().join(file), &target).map_err(|error| io(&target, error))?;
        }
        Ok(())
    }
}

/// The shipped cache of one program, after verification.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum ShippedOrderedCache {
    /// The data ships no cache for the program.
    None,
    Accepted(AcceptedOrderedCache),
    /// Refused; the stderr line has been printed.
    Rejected(ShippedOrderedCacheUse),
}

/// Why a shipped cache is not used.
enum Refusal {
    /// A verdict about the data (version, rule tables, file set, sha256):
    /// reported once and kept.
    Data(String),
    /// The request's resource limits stopped the verification: the request's
    /// error, never a verdict about the data.
    Budget(RustRedEvaluationError),
}

impl From<String> for Refusal {
    fn from(detail: String) -> Self {
        Self::Data(detail)
    }
}

/// A failed read of the pack file `file`.
fn read_refusal(file: &str, error: RustRedEvaluationError) -> Refusal {
    match error {
        RustRedEvaluationError::ResourceLimit(_) => Refusal::Budget(error),
        error => Refusal::Data(format!("{file}: {error}")),
    }
}

impl ShippedOrderedCache {
    /// Look up and verify the cache the loaded data ship for `program`. A
    /// refusal prints one line on stderr. Fails only when the request's
    /// resource limits stop the verification.
    pub(crate) fn resolve(
        data: &RustRedData,
        program: &OrderedCacheProgram,
    ) -> Result<Self, RustRedEvaluationError> {
        if program.textual_tables.is_some() {
            return super::p_shipped_cache::resolve(data, program);
        }
        let Some(pack) = data.loops.get(&program.loop_count) else {
            return Ok(Self::None);
        };
        let Some(section) = &pack.manifest.ordered_cache else {
            return Ok(Self::None);
        };
        let name = program_directory_name(program);
        let reject =
            |source: PathBuf, refusal: Refusal, entry: Option<&ShippedOrderedCacheProgram>| {
                let detail = match refusal {
                    Refusal::Data(detail) => detail,
                    Refusal::Budget(error) => return Err(error),
                };
                loud_rejection(&source, &program.label, &detail);
                Ok(Self::Rejected(ShippedOrderedCacheUse {
                    status: ShippedOrderedCacheStatus::Rejected,
                    source,
                    records: entry.map_or(0, |entry| entry.records),
                    bytes: entry.map_or(0, |entry| entry.bytes + entry.tail_bytes),
                    complete: entry.is_some_and(|entry| entry.complete),
                    detail: Some(detail),
                }))
            };
        let pack_path = |relative: &str| {
            pack.path(relative)
                .ok()
                .flatten()
                .unwrap_or_else(|| PathBuf::from(format!("{}L/{relative}", program.loop_count)))
        };
        let metadata = match metadata(pack, section) {
            Ok(metadata) => metadata,
            Err(refusal) => return reject(pack_path(&section.metadata), refusal, None),
        };
        let Some(entry) = metadata
            .programs
            .iter()
            .find(|entry| entry.label == program.label)
        else {
            return Ok(Self::None);
        };
        let source = match pack.path(&entry.directory) {
            Ok(Some(source)) => source,
            // Vendored or locally overridden data have no directory to read.
            Ok(None) => return Ok(Self::None),
            Err(error) => {
                return reject(
                    pack_path(&section.metadata),
                    Refusal::Data(error.to_string()),
                    Some(entry),
                );
            }
        };
        match verify(pack, program, entry, &name) {
            Ok(files) => Ok(Self::Accepted(AcceptedOrderedCache {
                source,
                files,
                entry: entry.clone(),
            })),
            Err(refusal) => reject(source, refusal, Some(entry)),
        }
    }
}

fn loud_rejection(source: &Path, label: &str, detail: &str) {
    eprintln!(
        "vakint: shipped Ordered step cache {} not used ({detail}); {label} runs without it",
        source.display()
    );
}

/// Read and check a pack's `ordered-cache.json` (sha256 checked by the pack).
fn metadata(
    pack: &LoopData,
    section: &OrderedCacheSection,
) -> Result<ShippedOrderedCacheMetadata, Refusal> {
    if section.path != PACK_DIRECTORY || section.metadata != METADATA_FILE {
        return Err(Refusal::Data(format!(
            "unsupported ordered_cache section {{path: {:?}, metadata: {:?}}}",
            section.path, section.metadata
        )));
    }
    let bytes = pack
        .read(METADATA_FILE)
        .map_err(|error| read_refusal(METADATA_FILE, error))?;
    let metadata: ShippedOrderedCacheMetadata =
        serde_json::from_slice(&bytes).map_err(|error| format!("{METADATA_FILE}: {error}"))?;
    if metadata.schema != METADATA_SCHEMA {
        return Err(Refusal::Data(format!(
            "{METADATA_FILE} has schema {:?}; this vakint reads {METADATA_SCHEMA:?}",
            metadata.schema
        )));
    }
    if metadata.loop_count != pack.manifest.loop_count {
        return Err(Refusal::Data(format!(
            "{METADATA_FILE} describes {} loops",
            metadata.loop_count
        )));
    }
    Ok(metadata)
}

/// The fields of a RustRed step-cache manifest that decide its use.
struct RustRedManifest {
    schema: String,
    semantics: String,
    adapter: String,
    stamp: String,
    /// `(file, records, bytes)` of every shard line.
    shards: Vec<(String, usize, u64)>,
}

fn parse_manifest(text: &str) -> Result<RustRedManifest, String> {
    let mut manifest = RustRedManifest {
        schema: String::new(),
        semantics: String::new(),
        adapter: String::new(),
        stamp: String::new(),
        shards: Vec::new(),
    };
    let mut ended = false;
    for line in text.lines() {
        let (key, value) = line.split_once(' ').unwrap_or((line, ""));
        match key {
            "schema" => manifest.schema = value.to_owned(),
            "semantics" => manifest.semantics = value.to_owned(),
            "adapter" => manifest.adapter = value.to_owned(),
            "stamp" => manifest.stamp = value.to_owned(),
            "shard" => {
                // shard <stage> <support> <class> <file> <records> <bytes> <digest>
                let fields: Vec<&str> = value.split(' ').collect();
                let [_, _, _, file, records, bytes, _] = fields[..] else {
                    return Err(format!("malformed shard line {line:?}"));
                };
                manifest.shards.push((
                    file.to_owned(),
                    records.parse().map_err(|_| format!("malformed {line:?}"))?,
                    bytes.parse().map_err(|_| format!("malformed {line:?}"))?,
                ));
            }
            "end" => ended = true,
            _ => {}
        }
    }
    if !ended || manifest.stamp.is_empty() {
        return Err(format!("truncated {RUSTRED_MANIFEST}"));
    }
    Ok(manifest)
}

/// The version header of a RustRed tail table: magic (16 bytes), the header
/// length (u64 little endian) and the header text, whose first lines are the
/// file schema and the tail semantics version.
fn check_tail_header(bytes: &[u8]) -> Result<(), String> {
    let header = bytes
        .get(16..24)
        .and_then(|length| usize::try_from(u64::from_le_bytes(length.try_into().ok()?)).ok())
        .and_then(|length| bytes.get(24..24usize.checked_add(length)?))
        .ok_or("truncated tail table header")?;
    let header = String::from_utf8_lossy(header);
    let mut lines = header.lines();
    let schema = lines.next().unwrap_or_default();
    let semantics = lines.next().unwrap_or_default();
    if schema != format!("schema {ORDERED_TAIL_SCHEMA}")
        || semantics != format!("semantics {ORDERED_TAIL_SEMANTICS}")
    {
        return Err(format!(
            "tail table version {schema:?} / {semantics:?}; this build reads \
             {ORDERED_TAIL_SCHEMA} / {ORDERED_TAIL_SEMANTICS}"
        ));
    }
    Ok(())
}

/// Check one shipped program cache; returns its files (relative to its
/// directory, RustRed manifest last).
fn verify(
    pack: &LoopData,
    program: &OrderedCacheProgram,
    entry: &ShippedOrderedCacheProgram,
    name: &str,
) -> Result<Vec<String>, Refusal> {
    let stamp = hex(&program.binding.stamp());
    let tables = hex(&program.rule_tables);
    if entry.rule_tables != tables {
        return Err(Refusal::Data(format!(
            "it was prepared for other rule tables ({}..., these data have {}...)",
            entry.rule_tables.get(..16).unwrap_or(&entry.rule_tables),
            &tables[..16]
        )));
    }
    if entry.stamp != stamp {
        return Err(Refusal::Data(format!(
            "it was recorded under another cache version (stamp {}..., this build {}...: \
             {ORDERED_STEP_CACHE_SCHEMA}, {ORDERED_STEP_SEMANTICS})",
            entry.stamp.get(..16).unwrap_or(&entry.stamp),
            &stamp[..16]
        )));
    }
    if entry.directory != format!("{PACK_DIRECTORY}/{name}") {
        return Err(Refusal::Data(format!(
            "its directory {:?} is not {PACK_DIRECTORY}/{name}",
            entry.directory
        )));
    }
    let prefix = format!("{}/", entry.directory);
    let listed: BTreeSet<String> = pack
        .manifest
        .files
        .range(prefix.clone()..)
        .map_while(|(path, _)| path.strip_prefix(&prefix).map(str::to_owned))
        .collect();
    let read = |file: &str| {
        pack.read(&format!("{prefix}{file}"))
            .map_err(|error| read_refusal(file, error))
    };
    if !listed.contains(RUSTRED_MANIFEST) {
        return Err(Refusal::Data(format!(
            "the loop manifest lists no {RUSTRED_MANIFEST}"
        )));
    }
    let manifest = parse_manifest(&String::from_utf8_lossy(&read(RUSTRED_MANIFEST)?))?;
    if manifest.schema != ORDERED_STEP_CACHE_SCHEMA || manifest.semantics != ORDERED_STEP_SEMANTICS
    {
        return Err(Refusal::Data(format!(
            "it was written as {} / {}; this build reads {ORDERED_STEP_CACHE_SCHEMA} / \
             {ORDERED_STEP_SEMANTICS}",
            manifest.schema, manifest.semantics
        )));
    }
    if manifest.stamp != stamp {
        return Err(Refusal::Data(format!(
            "its {RUSTRED_MANIFEST} carries stamp {}..., this build {}...",
            manifest.stamp.get(..16).unwrap_or(&manifest.stamp),
            &stamp[..16]
        )));
    }
    let shards: BTreeSet<String> = manifest
        .shards
        .iter()
        .map(|(file, _, _)| file.clone())
        .collect();
    let mut expected = shards.clone();
    expected.insert(RUSTRED_MANIFEST.to_owned());
    let tail = listed.contains(TAIL_FILE);
    if tail {
        expected.insert(TAIL_FILE.to_owned());
    }
    if listed != expected
        || shards
            .iter()
            .any(|file| !file.starts_with(&format!("{SHARD_DIRECTORY}/")))
    {
        return Err(Refusal::Data(format!(
            "its files differ from the shards its {RUSTRED_MANIFEST} lists"
        )));
    }
    for file in &shards {
        read(file)?;
    }
    if tail {
        check_tail_header(&read(TAIL_FILE)?)?;
    }
    let mut files: Vec<String> = shards.into_iter().collect();
    if tail {
        files.push(TAIL_FILE.to_owned());
    }
    files.push(RUSTRED_MANIFEST.to_owned());
    Ok(files)
}

/// Result of [`Vakint::export_ordered_caches`].
#[derive(Clone, Debug, Default, Serialize)]
pub struct OrderedCacheExportReport {
    /// Loop count -> exported programs.
    pub packs: BTreeMap<usize, Vec<ShippedOrderedCacheProgram>>,
    /// Program caches of this engine with no directory in the cache root.
    pub missing: Vec<String>,
    /// The exported cache of the RustRed p library, written as
    /// `p_integral/ordered-cache/<directory>/` and
    /// `p_integral/ordered-cache.json` (see `p_shipped_cache`).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub p_integral: Option<ShippedOrderedCacheProgram>,
}

/// The producer description of exports by this build.
pub fn ordered_cache_producer(
    vakint_revision: &str,
    rustred_revision: &str,
    command: &str,
) -> serde_json::Value {
    serde_json::json!({
        "vakint": {"version": env!("CARGO_PKG_VERSION"), "revision": vakint_revision},
        "rustred": {
            "revision": rustred_revision,
            "ordered_step_cache_schema": ORDERED_STEP_CACHE_SCHEMA,
            "ordered_step_semantics": ORDERED_STEP_SEMANTICS,
            "ordered_tail_schema": ORDERED_TAIL_SCHEMA,
            "ordered_tail_semantics": ORDERED_TAIL_SEMANTICS,
        },
        "command": command,
    })
}

impl Vakint {
    /// Every loaded program (the five-loop application once loaded) with its
    /// cache directory name, rule tables, binding stamp and the verified
    /// state of the cache its data ship: [`ShippedOrderedCacheStatus::InPlace`]
    /// when the cache would be used, `Rejected` (with the reason; the stderr
    /// line is printed) otherwise. Nothing is attached or copied. Fails when
    /// the active resource limits stop the verification.
    pub fn shipped_ordered_caches(&self) -> Result<Vec<ShippedOrderedCacheCheck>, VakintError> {
        let data = self.rustred_data()?;
        let mut out = Vec::new();
        for program in self.ordered_cache_programs()? {
            let shipped = match ShippedOrderedCache::resolve(&data, &program)? {
                ShippedOrderedCache::None => None,
                ShippedOrderedCache::Accepted(accepted) => {
                    Some(accepted.report(ShippedOrderedCacheStatus::InPlace, None))
                }
                ShippedOrderedCache::Rejected(report) => Some(report),
            };
            out.push(ShippedOrderedCacheCheck {
                loop_count: program.loop_count,
                directory: program_directory_name(&program),
                rule_tables: hex(&program.rule_tables),
                stamp: hex(&program.binding.stamp()),
                program: program.label,
                shipped,
            });
        }
        Ok(out)
    }

    /// Export the caches that an engine cache root `cache_root` (the
    /// `directory` of [`crate::OrderedCacheSettings`]) holds for this engine's
    /// loaded programs into the external-data layout under `output`:
    /// `<L>L/ordered-cache/<label>.<tables>.<stamp>/` (the RustRed manifest,
    /// the shards it lists and the tail table, byte for byte) and
    /// `<L>L/ordered-cache.json` with `producer` and `workload`. `output` must
    /// not exist. The files depend only on the recorded steps, so the same
    /// workload exports the same bytes.
    pub fn export_ordered_caches(
        &self,
        cache_root: &Path,
        output: &Path,
        producer: serde_json::Value,
        workload: Vec<serde_json::Value>,
    ) -> Result<OrderedCacheExportReport, VakintError> {
        let failure = |detail: String| -> VakintError { ordered_cache_failure(detail).into() };
        if output.exists() {
            return Err(failure(format!("{} already exists", output.display())));
        }
        let mut report = OrderedCacheExportReport::default();
        for program in self.ordered_cache_programs()? {
            if program.textual_tables.is_some() {
                match super::p_shipped_cache::export(
                    cache_root, output, &program, &producer, &workload,
                )? {
                    Some(entry) => report.p_integral = Some(entry),
                    None => report.missing.push(program_directory_name(&program)),
                }
                continue;
            }
            let name = program_directory_name(&program);
            let source = cache_root.join(&name);
            let manifest_path = source.join(RUSTRED_MANIFEST);
            if !manifest_path.is_file() {
                report.missing.push(name);
                continue;
            }
            let text = fs::read_to_string(&manifest_path)
                .map_err(|error| failure(format!("{}: {error}", manifest_path.display())))?;
            let manifest = parse_manifest(&text)
                .map_err(|error| failure(format!("{}: {error}", manifest_path.display())))?;
            if manifest.stamp != hex(&program.binding.stamp()) {
                return Err(failure(format!(
                    "{name}: the manifest carries another stamp"
                )));
            }
            let directory = format!("{PACK_DIRECTORY}/{name}");
            let target = output
                .join(format!("{}L", program.loop_count))
                .join(&directory);
            let mut files: Vec<String> = manifest
                .shards
                .iter()
                .map(|(file, _, _)| file.clone())
                .collect();
            let tail_bytes = match fs::metadata(source.join(TAIL_FILE)) {
                Ok(metadata) => {
                    files.push(TAIL_FILE.to_owned());
                    metadata.len()
                }
                Err(_) => 0,
            };
            files.push(RUSTRED_MANIFEST.to_owned());
            for file in &files {
                let to = target.join(file);
                if let Some(parent) = to.parent() {
                    fs::create_dir_all(parent)
                        .map_err(|error| failure(format!("{}: {error}", parent.display())))?;
                }
                fs::copy(source.join(file), &to)
                    .map_err(|error| failure(format!("{name}/{file}: {error}")))?;
            }
            report
                .packs
                .entry(program.loop_count)
                .or_default()
                .push(ShippedOrderedCacheProgram {
                    label: program.label.clone(),
                    directory,
                    rule_tables: hex(&program.rule_tables),
                    stamp: manifest.stamp.clone(),
                    adapter: manifest.adapter.clone(),
                    records: manifest.shards.iter().map(|(_, records, _)| records).sum(),
                    shards: manifest.shards.len(),
                    bytes: manifest.shards.iter().map(|(_, _, bytes)| bytes).sum(),
                    tail_bytes,
                    complete: true,
                    trim: None,
                });
        }
        if report.packs.is_empty() && report.p_integral.is_none() {
            return Err(failure(format!(
                "{} holds no Ordered step cache of the loaded programs",
                cache_root.display()
            )));
        }
        for (&loops, programs) in &mut report.packs {
            programs.sort_by(|left, right| left.label.cmp(&right.label));
            let metadata = ShippedOrderedCacheMetadata {
                schema: METADATA_SCHEMA.to_owned(),
                loop_count: loops,
                producer: producer.clone(),
                workload: workload
                    .iter()
                    .filter(|row| {
                        row["loop_count"]
                            .as_u64()
                            .is_none_or(|count| count == loops as u64)
                    })
                    .cloned()
                    .collect(),
                programs: programs.clone(),
            };
            let path = output.join(format!("{loops}L")).join(METADATA_FILE);
            let mut text = serde_json::to_string_pretty(&metadata)
                .map_err(|error| failure(error.to_string()))?;
            text.push('\n');
            fs::write(&path, text)
                .map_err(|error| failure(format!("{}: {error}", path.display())))?;
        }
        Ok(report)
    }
}

#[cfg(test)]
#[path = "shipped_cache_tests.rs"]
mod tests;
