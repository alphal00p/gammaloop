//! The shipped Ordered step cache of the RustRed p library (M2-V6).
//!
//! External data read from a directory may declare a top-level `p_integral/`
//! area through the root-manifest key `p_integral` (vakint-external-data
//! branch `p_integral`):
//!
//! ```text
//! manifest.json            ... + "p_integral": {"manifest": "p_integral/manifest.json",
//!                                               "manifest_sha256": "<hex>"}
//! p_integral/manifest.json {"schema": "vakint.external-data.p-integral.v1", "libraries": [],
//!                           "ordered_cache": {"path": "ordered-cache", "metadata": "ordered-cache.json"},
//!                           "files": {"<path below p_integral/>": {"sha256": ..., "bytes": ...}}}
//! p_integral/ordered-cache.json            schema vakint.external-data.p-ordered-cache.v1
//! p_integral/ordered-cache/<directory>/    ordered-step-cache.txt, shards/*, [ordered-tail.rrtail]
//! ```
//!
//! `<directory>` is the textual cache directory of the library program
//! ([`super::p_data::p_cache_directory_name`]), the one an engine cache
//! directory uses. Validity is decided by versions only, as for the vacuum
//! caches of a loop pack ([`super::shipped_cache`]): the cache of the running
//! program is used iff its metadata entry carries the program's label,
//! directory, library identity and binding stamp, its RustRed manifest carries
//! this build's step-cache schema, semantics and stamp, a tail table carries
//! this build's tail versions, and every file of the cache directory is listed
//! by the area manifest with matching sha256 and bytes (integrity). Otherwise
//! vakint prints one line on stderr and the program runs without it (cold);
//! only a resource limit of the request propagates.
//!
//! No root key, a root key whose area manifest is absent (a packed or managed
//! install keeps the key), an area without an Ordered cache or without an
//! entry for the program: no shipped cache, silently. An injected library
//! never uses one, and vendored data have no directory to read.
//!
//! Once accepted, the cache is used exactly like a loop pack's shipped cache
//! (opened in place, installed into a read-write cache directory, replayed
//! read-only; see [`super::cache`]). [`export`] writes the area layout below
//! an export directory for [`Vakint::export_ordered_caches`].
//!
//! The area may also ship RustRed p libraries (M3-V8; `tools/p_integral_area.py
//! stage-library`), each listed in `libraries`:
//!
//! ```text
//! {"directory": "p-library-<a>-<b>", "schema": <its P_LIBRARY_SCHEMA>, "scope": [a, b],
//!  "identity": <its PLibrary identity>, "files": {<path below the directory>: {"sha256", "bytes"}}}
//! p_integral/p-library-<a>-<b>/            the `rustred p-package` output
//! ```
//!
//! [`listed_library`] adopts only a listed library whose schema is this
//! RustRed's `P_LIBRARY_SCHEMA`. A listed library of another schema, a
//! `p-library-<a>-<b>` directory the manifest does not list, an area manifest
//! of another area schema or one that cannot be read are ignored with one
//! stderr line each (once per area manifest and process). Only a listed
//! library of the supported schema that fails verification is an error (the
//! engine's p data are broken): an area manifest that does not match the root
//! key, several such libraries, a malformed entry, a missing directory, or,
//! when the library is loaded ([`verify_listed_files`]), a file inventory or
//! sha256 that differs from the entry. The area schema stays
//! `vakint.external-data.p-integral.v1`: `libraries` and the library
//! directories are additive, and the cache reader above ignores both.
use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    ops::RangeInclusive,
    path::{Path, PathBuf},
    sync::{Mutex, PoisonError},
};

use rustred::propagator::P_LIBRARY_SCHEMA;
use rustred::reduction::{
    ORDERED_STEP_CACHE_SCHEMA, ORDERED_STEP_SEMANTICS, ORDERED_TAIL_SCHEMA, ORDERED_TAIL_SEMANTICS,
};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use symbolica::utils::resources;

use super::{
    RustRedEvaluationError,
    cache::{OrderedCacheProgram, hex, ordered_cache_failure, program_directory_name},
    data::{FileDigest, RustRedData, relative_path},
    p_data::directory_scope,
    shipped_cache::{
        AcceptedOrderedCache, OrderedCacheSection, RUSTRED_MANIFEST, ShippedOrderedCache,
        ShippedOrderedCacheProgram, ShippedOrderedCacheStatus, ShippedOrderedCacheUse,
    },
};
use crate::VakintError;

/// Root-manifest key declaring the area.
pub(crate) const ROOT_KEY: &str = "p_integral";
/// The area directory under the external-data root.
pub(crate) const AREA_DIRECTORY: &str = "p_integral";
/// Schema of `p_integral/manifest.json`.
pub(crate) const AREA_SCHEMA: &str = "vakint.external-data.p-integral.v1";
/// Schema of `p_integral/ordered-cache.json`.
pub(crate) const METADATA_SCHEMA: &str = "vakint.external-data.p-ordered-cache.v1";
const AREA_MANIFEST: &str = "manifest.json";
const CACHE_DIRECTORY: &str = "ordered-cache";
const METADATA_FILE: &str = "ordered-cache.json";
const TAIL_FILE: &str = "ordered-tail.rrtail";
const SHARD_DIRECTORY: &str = "shards";
/// The data version of an injected library's program
/// (`p_data::PLibraryOrigin::Injected`).
const INJECTED: &str = "injected";

/// The root-manifest key.
#[derive(Deserialize)]
struct RootKey {
    manifest: String,
    manifest_sha256: String,
}

/// `p_integral/manifest.json`, as far as the cache reader goes.
#[derive(Deserialize)]
struct AreaManifest {
    schema: String,
    // `libraries` lists p libraries shipped in the area ([`listed_library`]);
    // the cache reader ignores it.
    #[serde(default)]
    ordered_cache: Option<OrderedCacheSection>,
    files: BTreeMap<String, FileDigest>,
}

/// One program of `p_integral/ordered-cache.json`.
#[derive(Clone, Debug, Deserialize, Serialize)]
pub(crate) struct PCacheProgram {
    /// `p-library-<a>-<b>`.
    pub label: String,
    pub scope: [u8; 2],
    /// `ordered-cache/<textual directory name>`, relative to the area.
    pub directory: String,
    /// `PLibrary::identity()` of the program.
    pub identity: String,
    /// The `P_*` versions of the producing build (provenance).
    #[serde(default)]
    pub p_versions: BTreeMap<String, String>,
    pub stamp: String,
    pub adapter: String,
    pub records: usize,
    pub shards: usize,
    pub bytes: u64,
    pub tail_bytes: u64,
    pub complete: bool,
    #[serde(default)]
    pub trim: Option<serde_json::Value>,
}

/// `p_integral/ordered-cache.json`.
#[derive(Deserialize, Serialize)]
struct PCacheMetadata {
    schema: String,
    #[serde(default)]
    producer: serde_json::Value,
    #[serde(default)]
    workload: Vec<serde_json::Value>,
    programs: Vec<PCacheProgram>,
}

impl PCacheProgram {
    /// The vacuum-shaped entry of reports: zero rule tables (the p program is
    /// identified by its textual directory).
    fn shipped_entry(&self) -> ShippedOrderedCacheProgram {
        ShippedOrderedCacheProgram {
            label: self.label.clone(),
            directory: self.directory.clone(),
            rule_tables: "0".repeat(64),
            stamp: self.stamp.clone(),
            adapter: self.adapter.clone(),
            records: self.records,
            shards: self.shards,
            bytes: self.bytes,
            tail_bytes: self.tail_bytes,
            complete: self.complete,
            trim: self.trim.clone(),
        }
    }
}

/// Why the area is not used.
enum Refusal {
    /// A verdict about the data: one stderr line, kept.
    Data(String),
    /// The request's resource limits stopped the verification.
    Budget(RustRedEvaluationError),
}

impl From<String> for Refusal {
    fn from(detail: String) -> Self {
        Self::Data(detail)
    }
}

/// A failed read of `file`: only a resource limit is the request's error.
fn read_refusal(file: &str, error: RustRedEvaluationError) -> Refusal {
    match error {
        RustRedEvaluationError::ResourceLimit(_) => Refusal::Budget(error),
        error => Refusal::Data(format!("{file}: {error}")),
    }
}

/// The external-data root of `data`, if they were read from a directory
/// (as `Vakint::external_data_root` finds it).
fn data_root(data: &RustRedData) -> Option<PathBuf> {
    data.loops.values().find_map(|pack| {
        let manifest = pack.path("manifest.json").ok().flatten()?;
        Some(manifest.parent()?.parent()?.to_owned())
    })
}

/// `root/relative`, refusing symlinks on the way (as the loop packs do);
/// `Ok(None)` when it does not exist.
fn read_file(root: &Path, relative: &str) -> Result<Option<Vec<u8>>, RustRedEvaluationError> {
    crate::execution::checkpoint()?;
    let mut current = root.to_owned();
    for part in relative_path(relative)?.components() {
        current.push(part.as_os_str());
        match fs::symlink_metadata(&current) {
            Ok(metadata) if metadata.file_type().is_symlink() => {
                return Err(super::data::failure(format!(
                    "symlink in {}",
                    current.display()
                )));
            }
            Ok(_) => {}
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
            Err(error) => return Err(super::data::failure(error)),
        }
    }
    resources::check_allocation(fs::metadata(&current).map_err(super::data::failure)?.len())?;
    fs::read(&current).map(Some).map_err(super::data::failure)
}

fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

/// The area of one external-data root, as far as the root key and the area
/// manifest go: `Ok(None)` without a root key or without an area manifest.
struct Area {
    root: PathBuf,
    manifest_sha256: String,
}

fn area_of(root: &Path) -> Result<Option<Area>, Refusal> {
    let bytes = read_file(root, "manifest.json")
        .map_err(|error| read_refusal("manifest.json", error))?
        .ok_or_else(|| format!("{} has no manifest.json", root.display()))?;
    let manifest: serde_json::Value =
        serde_json::from_slice(&bytes).map_err(|error| format!("manifest.json: {error}"))?;
    let Some(key) = manifest.get(ROOT_KEY) else {
        return Ok(None);
    };
    let key: RootKey = serde_json::from_value(key.clone())
        .map_err(|error| format!("manifest.json key {ROOT_KEY}: {error}"))?;
    if key.manifest != format!("{AREA_DIRECTORY}/{AREA_MANIFEST}") {
        return Err(Refusal::Data(format!(
            "manifest.json key {ROOT_KEY} names {:?}",
            key.manifest
        )));
    }
    if fs::symlink_metadata(root.join(&key.manifest))
        .is_err_and(|error| error.kind() == std::io::ErrorKind::NotFound)
    {
        return Ok(None);
    }
    Ok(Some(Area {
        root: root.to_owned(),
        manifest_sha256: key.manifest_sha256,
    }))
}

impl Area {
    /// The program's cache in the area: `Ok(None)` when the area ships none.
    fn verify(
        &self,
        program: &OrderedCacheProgram,
    ) -> Result<Option<AcceptedOrderedCache>, Refusal> {
        let manifest_path = format!("{AREA_DIRECTORY}/{AREA_MANIFEST}");
        let bytes = read_file(&self.root, &manifest_path)
            .map_err(|error| read_refusal(&manifest_path, error))?
            .ok_or_else(|| format!("{manifest_path} disappeared"))?;
        if sha256(&bytes) != self.manifest_sha256 {
            return Err(Refusal::Data(format!(
                "{manifest_path} does not match the root manifest's sha256"
            )));
        }
        let area: AreaManifest =
            serde_json::from_slice(&bytes).map_err(|error| format!("{manifest_path}: {error}"))?;
        if area.schema != AREA_SCHEMA {
            return Err(Refusal::Data(format!(
                "{manifest_path} has schema {:?}; this vakint reads {AREA_SCHEMA:?}",
                area.schema
            )));
        }
        let Some(section) = &area.ordered_cache else {
            return Ok(None);
        };
        if section.path != CACHE_DIRECTORY || section.metadata != METADATA_FILE {
            return Err(Refusal::Data(format!(
                "unsupported ordered_cache section {{path: {:?}, metadata: {:?}}}",
                section.path, section.metadata
            )));
        }
        // Every read file is listed in the area manifest with its sha256 and size.
        let read = |file: &str| -> Result<Vec<u8>, Refusal> {
            let digest = area
                .files
                .get(file)
                .ok_or_else(|| format!("{manifest_path} does not list {file}"))?;
            let path = format!("{AREA_DIRECTORY}/{file}");
            let bytes = read_file(&self.root, &path)
                .map_err(|error| read_refusal(&path, error))?
                .ok_or_else(|| format!("{path} is missing"))?;
            if bytes.len() as u64 != digest.bytes || sha256(&bytes) != digest.sha256 {
                return Err(Refusal::Data(format!(
                    "{path} does not match its sha256 and size in {manifest_path}"
                )));
            }
            Ok(bytes)
        };
        let metadata: PCacheMetadata = serde_json::from_slice(&read(METADATA_FILE)?)
            .map_err(|error| format!("{METADATA_FILE}: {error}"))?;
        if metadata.schema != METADATA_SCHEMA {
            return Err(Refusal::Data(format!(
                "{METADATA_FILE} has schema {:?}; this vakint reads {METADATA_SCHEMA:?}",
                metadata.schema
            )));
        }
        let Some(entry) = metadata
            .programs
            .iter()
            .find(|entry| entry.label == program.label)
        else {
            return Ok(None);
        };
        let name = program_directory_name(program);
        let stamp = hex(&program.binding.stamp());
        if entry.directory != format!("{CACHE_DIRECTORY}/{name}") {
            return Err(Refusal::Data(format!(
                "its directory {:?} is not {CACHE_DIRECTORY}/{name}",
                entry.directory
            )));
        }
        if entry.identity != program.binding.family_fingerprint {
            return Err(Refusal::Data(format!(
                "it was prepared for the library {:?}, this build runs {:?}",
                entry.identity, program.binding.family_fingerprint
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
        let prefix = format!("{}/", entry.directory);
        let listed: BTreeSet<String> = area
            .files
            .keys()
            .filter_map(|path| path.strip_prefix(&prefix).map(str::to_owned))
            .collect();
        let read_cache = |file: &str| read(&format!("{prefix}{file}"));
        let manifest = parse_manifest(&String::from_utf8_lossy(&read_cache(RUSTRED_MANIFEST)?))?;
        if manifest.versions.0 != ORDERED_STEP_CACHE_SCHEMA
            || manifest.versions.1 != ORDERED_STEP_SEMANTICS
            || manifest.stamp != stamp
        {
            return Err(Refusal::Data(format!(
                "its {RUSTRED_MANIFEST} was written as {} / {} (stamp {}...); this build \
                 reads {ORDERED_STEP_CACHE_SCHEMA} / {ORDERED_STEP_SEMANTICS} (stamp {}...)",
                manifest.versions.0,
                manifest.versions.1,
                manifest.stamp.get(..16).unwrap_or(&manifest.stamp),
                &stamp[..16]
            )));
        }
        let tail = listed.contains(TAIL_FILE);
        let mut expected: BTreeSet<String> = manifest.shards.iter().cloned().collect();
        expected.insert(RUSTRED_MANIFEST.to_owned());
        if tail {
            expected.insert(TAIL_FILE.to_owned());
        }
        if listed != expected
            || manifest
                .shards
                .iter()
                .any(|file| !file.starts_with(&format!("{SHARD_DIRECTORY}/")))
        {
            return Err(Refusal::Data(format!(
                "its files differ from the shards its {RUSTRED_MANIFEST} lists"
            )));
        }
        let mut files = Vec::new();
        for file in &manifest.shards {
            read_cache(file)?;
            files.push(file.clone());
        }
        if tail {
            check_tail_header(&read_cache(TAIL_FILE)?)?;
            files.push(TAIL_FILE.to_owned());
        }
        files.push(RUSTRED_MANIFEST.to_owned());
        Ok(Some(AcceptedOrderedCache {
            source: self.root.join(AREA_DIRECTORY).join(&entry.directory),
            files,
            entry: entry.shipped_entry(),
        }))
    }
}

/// One file of a listed library, by its path below the library directory.
#[derive(Clone, Debug, PartialEq, Eq, Deserialize)]
pub(crate) struct ListedFile {
    pub sha256: String,
    pub bytes: u64,
}

/// A `libraries` entry of the supported schema (module documentation).
#[derive(Deserialize)]
struct LibraryEntry {
    directory: String,
    scope: [u8; 2],
    identity: String,
    files: BTreeMap<String, ListedFile>,
}

/// A library the area lists with this RustRed's `P_LIBRARY_SCHEMA`, checked
/// as far as its entry goes ([`verify_listed_files`] checks the files).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct ListedLibrary {
    /// `<root>/p_integral/p-library-<a>-<b>`.
    pub path: PathBuf,
    pub scope: RangeInclusive<u8>,
    /// The `PLibrary::identity()` the entry records.
    pub identity: String,
    pub files: BTreeMap<String, ListedFile>,
}

/// What the area manifest of one root lists.
#[derive(Debug)]
pub(crate) struct ListedResolution {
    /// `Ok(None)` without a listed library of the supported schema; `Err` for
    /// one that cannot be used (the detail of the engine's error).
    pub library: Result<Option<ListedLibrary>, String>,
    /// One line per ignored library or unreadable area manifest.
    pub notes: Vec<String>,
}

/// The listed library of `root`'s area, given the area manifest sha256 the
/// root key declares (module documentation). Reads only the area manifest and
/// the area directory listing; fails only when the request's resource limits
/// stop the read.
pub(crate) fn resolve_listed(
    root: &Path,
    declared_sha256: &str,
) -> Result<ListedResolution, RustRedEvaluationError> {
    let manifest_path = format!("{AREA_DIRECTORY}/{AREA_MANIFEST}");
    let ignored = |detail: String| ListedResolution {
        library: Ok(None),
        notes: vec![format!("{detail}; its p libraries are ignored")],
    };
    let bytes = match read_file(root, &manifest_path) {
        Ok(Some(bytes)) => bytes,
        Ok(None) => {
            return Ok(ListedResolution {
                library: Ok(None),
                notes: Vec::new(),
            });
        }
        Err(error @ RustRedEvaluationError::ResourceLimit(_)) => return Err(error),
        Err(error) => return Ok(ignored(format!("{manifest_path}: {error}"))),
    };
    let area: serde_json::Value = match serde_json::from_slice(&bytes) {
        Ok(area) => area,
        Err(error) => return Ok(ignored(format!("{manifest_path}: {error}"))),
    };
    let schema = area.get("schema").and_then(serde_json::Value::as_str);
    if schema != Some(AREA_SCHEMA) {
        return Ok(ignored(format!(
            "{manifest_path} has the schema {schema:?}; this vakint reads {AREA_SCHEMA:?}"
        )));
    }
    let entries = match area.get("libraries") {
        None => Vec::new(),
        Some(serde_json::Value::Array(entries)) => entries.clone(),
        Some(_) => return Ok(ignored(format!("{manifest_path}: libraries is not a list"))),
    };
    let mut notes = Vec::new();
    let mut supported = Vec::new();
    let mut listed = BTreeSet::new();
    for entry in entries {
        let directory = entry
            .get("directory")
            .and_then(serde_json::Value::as_str)
            .unwrap_or("(no directory)")
            .to_owned();
        match entry.get("schema").and_then(serde_json::Value::as_str) {
            Some(schema) if schema == P_LIBRARY_SCHEMA => supported.push(entry),
            schema => notes.push(format!(
                "{AREA_DIRECTORY}/{directory} has the library schema {schema:?}; this vakint \
                 reads {P_LIBRARY_SCHEMA:?}; ignored"
            )),
        }
        listed.insert(directory);
    }
    let area_directory = root.join(AREA_DIRECTORY);
    if let Ok(read) = fs::read_dir(&area_directory) {
        let mut unlisted = read
            .filter_map(|entry| entry.ok()?.file_name().into_string().ok())
            .filter(|name| directory_scope(name).is_some() && !listed.contains(name))
            .collect::<Vec<_>>();
        unlisted.sort();
        notes.extend(unlisted.into_iter().map(|name| {
            format!("{AREA_DIRECTORY}/{name} is not listed in {manifest_path}; ignored")
        }));
    }
    let library = match supported.as_slice() {
        [] => Ok(None),
        _ if sha256(&bytes) != declared_sha256 => Err(format!(
            "{manifest_path} does not match the root manifest's sha256"
        )),
        [entry] => listed_entry(&area_directory, entry.clone()).map(Some),
        _ => Err(format!(
            "{manifest_path} lists {} p libraries of the schema {P_LIBRARY_SCHEMA:?}; one \
             engine uses one library program",
            supported.len()
        )),
    };
    Ok(ListedResolution { library, notes })
}

/// The entry `value` of the supported schema under `area`: its directory
/// names its scope, its identity is of that scope, its file paths are plain
/// relative paths and the directory exists (no symlink).
fn listed_entry(area: &Path, value: serde_json::Value) -> Result<ListedLibrary, String> {
    let entry: LibraryEntry = serde_json::from_value(value)
        .map_err(|error| format!("a {P_LIBRARY_SCHEMA:?} entry of libraries: {error}"))?;
    let name = &entry.directory;
    let scope = directory_scope(name)
        .filter(|scope| [*scope.start(), *scope.end()] == entry.scope)
        .ok_or_else(|| {
            format!(
                "the library directory {name:?} does not declare the scope {:?}",
                entry.scope
            )
        })?;
    let [start, end] = entry.scope;
    if !entry
        .identity
        .starts_with(&format!("p-library:scope={start}-{end}:"))
    {
        return Err(format!(
            "{name}: the identity {:?} is not of the scope {start}-{end}",
            entry.identity
        ));
    }
    if entry.files.is_empty() {
        return Err(format!("{name}: no files are listed"));
    }
    for file in entry.files.keys() {
        relative_path(file).map_err(|error| format!("{name}: {error}"))?;
    }
    let path = area.join(name);
    if !fs::symlink_metadata(&path).is_ok_and(|metadata| metadata.is_dir()) {
        return Err(format!(
            "the listed library directory {} is missing (symlinks are refused)",
            path.display()
        ));
    }
    Ok(ListedLibrary {
        path,
        scope,
        identity: entry.identity,
        files: entry.files,
    })
}

/// Listed-library verdicts by external-data root and the area manifest sha256
/// its root key declares (a restaged area declares another one): resolved
/// once, their notes printed once.
type ListedVerdict = Result<Option<ListedLibrary>, String>;

fn listed_memo() -> &'static Mutex<BTreeMap<(PathBuf, String), ListedVerdict>> {
    static MEMO: Mutex<BTreeMap<(PathBuf, String), ListedVerdict>> = Mutex::new(BTreeMap::new());
    &MEMO
}

/// The library the `p_integral/` area of the external-data `root` lists, if
/// any ([`resolve_listed`], memoised; each note is printed once on stderr).
/// An error is a listed library of the supported schema that cannot be used.
pub(crate) fn listed_library(root: &Path) -> Result<Option<ListedLibrary>, RustRedEvaluationError> {
    let declared = match area_of(root) {
        Ok(Some(area)) => area.manifest_sha256,
        Ok(None) => return Ok(None),
        Err(Refusal::Budget(error)) => return Err(error),
        // The cache reader reports a malformed root key; no library is listed.
        Err(Refusal::Data(_)) => return Ok(None),
    };
    let key = (root.to_owned(), declared);
    if let Some(known) = listed_memo()
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .get(&key)
    {
        return known.clone().map_err(super::data::failure);
    }
    let resolution = resolve_listed(root, &key.1)?;
    for note in &resolution.notes {
        eprintln!("vakint: external data {}: {note}", root.display());
    }
    listed_memo()
        .lock()
        .unwrap_or_else(PoisonError::into_inner)
        .insert(key, resolution.library.clone());
    resolution.library.map_err(super::data::failure)
}

/// The files of `library` are exactly its entry's, each with its sha256 and
/// size; no symlink below its directory.
pub(crate) fn verify_listed_files(library: &ListedLibrary) -> Result<(), RustRedEvaluationError> {
    fn visit(
        root: &Path,
        directory: &Path,
        found: &mut BTreeSet<String>,
    ) -> Result<(), RustRedEvaluationError> {
        for entry in fs::read_dir(directory).map_err(super::data::failure)? {
            let path = entry.map_err(super::data::failure)?.path();
            let metadata = fs::symlink_metadata(&path).map_err(super::data::failure)?;
            if metadata.is_dir() {
                visit(root, &path, found)?;
            } else if metadata.is_file() {
                let relative = path.strip_prefix(root).map_err(super::data::failure)?;
                found.insert(relative.to_string_lossy().replace('\\', "/"));
            } else {
                return Err(super::data::failure(format!(
                    "{} is not a plain file (symlinks are refused)",
                    path.display()
                )));
            }
        }
        Ok(())
    }
    let mut found = BTreeSet::new();
    visit(&library.path, &library.path, &mut found)?;
    let listed: BTreeSet<String> = library.files.keys().cloned().collect();
    if found != listed {
        let first = |set: BTreeSet<&String>| set.into_iter().take(3).cloned().collect::<Vec<_>>();
        return Err(super::data::failure(format!(
            "{}: the files differ from its {AREA_DIRECTORY}/{AREA_MANIFEST} entry (unlisted {:?}, \
             missing {:?})",
            library.path.display(),
            first(found.difference(&listed).collect()),
            first(listed.difference(&found).collect()),
        )));
    }
    for (name, digest) in &library.files {
        let bytes = read_file(&library.path, name)?
            .ok_or_else(|| super::data::failure(format!("{name} disappeared")))?;
        if bytes.len() as u64 != digest.bytes || sha256(&bytes) != digest.sha256 {
            return Err(super::data::failure(format!(
                "{}/{name} does not match its sha256 and size in {AREA_DIRECTORY}/{AREA_MANIFEST}",
                library.path.display()
            )));
        }
    }
    Ok(())
}

/// The fields of a RustRed step-cache manifest that decide its use.
struct RustRedManifest {
    /// Schema and step semantics.
    versions: (String, String),
    stamp: String,
    adapter: String,
    /// `(file, records, bytes)` of every shard line.
    shard_rows: Vec<(String, usize, u64)>,
    shards: Vec<String>,
}

fn parse_manifest(text: &str) -> Result<RustRedManifest, String> {
    let mut manifest = RustRedManifest {
        versions: (String::new(), String::new()),
        stamp: String::new(),
        adapter: String::new(),
        shard_rows: Vec::new(),
        shards: Vec::new(),
    };
    let mut ended = false;
    for line in text.lines() {
        let (key, value) = line.split_once(' ').unwrap_or((line, ""));
        match key {
            "schema" => manifest.versions.0 = value.to_owned(),
            "semantics" => manifest.versions.1 = value.to_owned(),
            "adapter" => manifest.adapter = value.to_owned(),
            "stamp" => manifest.stamp = value.to_owned(),
            "shard" => {
                // shard <stage> <support> <class> <file> <records> <bytes> <digest>
                let fields: Vec<&str> = value.split(' ').collect();
                let [_, _, _, file, records, bytes, _] = fields[..] else {
                    return Err(format!("malformed shard line {line:?}"));
                };
                let number = || format!("malformed shard line {line:?}");
                manifest.shard_rows.push((
                    file.to_owned(),
                    records.parse().map_err(|_| number())?,
                    bytes.parse().map_err(|_| number())?,
                ));
                manifest.shards.push(file.to_owned());
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
    let (schema, semantics) = (
        lines.next().unwrap_or_default(),
        lines.next().unwrap_or_default(),
    );
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

fn loud_rejection(source: &Path, label: &str, detail: &str) {
    eprintln!(
        "vakint: shipped p cache {} not used ({detail}); {label} runs without it",
        source.display()
    );
}

/// Look up and verify the cache the loaded data ship for the p `program`
/// (unmemoised; [`ShippedOrderedCache::resolve`] of a textual program). A
/// refusal prints one line on stderr. Fails only when the request's resource
/// limits stop the verification.
pub(crate) fn resolve(
    data: &RustRedData,
    program: &OrderedCacheProgram,
) -> Result<ShippedOrderedCache, RustRedEvaluationError> {
    resolve_memo(data, program, false)
}

/// [`resolve`], memoised per external-data root, area manifest (its sha256
/// declared by the root key) and program: the shipped p cache of
/// [`super::cache`]'s `cache()`. Only verdicts about the data are kept.
pub(crate) fn shipped(
    data: &RustRedData,
    program: &OrderedCacheProgram,
) -> Result<ShippedOrderedCache, RustRedEvaluationError> {
    resolve_memo(data, program, true)
}

type MemoKey = (PathBuf, String, String, String);

fn memo() -> &'static Mutex<BTreeMap<MemoKey, ShippedOrderedCache>> {
    static MEMO: Mutex<BTreeMap<MemoKey, ShippedOrderedCache>> = Mutex::new(BTreeMap::new());
    &MEMO
}

fn resolve_memo(
    data: &RustRedData,
    program: &OrderedCacheProgram,
    memoise: bool,
) -> Result<ShippedOrderedCache, RustRedEvaluationError> {
    if program.textual_tables.is_none() || program.data_version == INJECTED {
        return Ok(ShippedOrderedCache::None);
    }
    let Some(root) = data_root(data) else {
        return Ok(ShippedOrderedCache::None);
    };
    let reject = |source: PathBuf, refusal: Refusal| match refusal {
        Refusal::Budget(error) => Err(error),
        Refusal::Data(detail) => {
            loud_rejection(&source, &program.label, &detail);
            Ok(ShippedOrderedCache::Rejected(ShippedOrderedCacheUse {
                status: ShippedOrderedCacheStatus::Rejected,
                source,
                records: 0,
                bytes: 0,
                complete: false,
                detail: Some(detail),
            }))
        }
    };
    let area_path = root.join(AREA_DIRECTORY);
    let area = match area_of(&root) {
        Ok(Some(area)) => area,
        Ok(None) => return Ok(ShippedOrderedCache::None),
        Err(refusal) => return reject(area_path, refusal),
    };
    let key = (
        root,
        area.manifest_sha256.clone(),
        program_directory_name(program),
        hex(&program.binding.stamp()),
    );
    if memoise
        && let Some(known) = memo()
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .get(&key)
    {
        return Ok(known.clone());
    }
    let resolved = match area.verify(program) {
        Ok(Some(accepted)) => ShippedOrderedCache::Accepted(accepted),
        Ok(None) => ShippedOrderedCache::None,
        Err(refusal) => reject(area_path, refusal)?,
    };
    if memoise {
        memo()
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .insert(key, resolved.clone());
    }
    Ok(resolved)
}

/// `[a, b]` of a label `p-library-<a>-<b>`.
fn label_scope(label: &str) -> Option<[u8; 2]> {
    let (start, end) = label.strip_prefix("p-library-")?.split_once('-')?;
    Some([start.parse().ok()?, end.parse().ok()?])
}

/// Export the cache an engine cache root holds for the p `program` into
/// `output/p_integral/`: `ordered-cache/<directory>/` (the RustRed manifest,
/// its shards and the tail table, byte for byte) and `ordered-cache.json`
/// with `producer` and the `workload` rows within the program's scope.
/// `Ok(None)` when the root holds no cache of the program; an injected
/// library's cache is never exported. Staging into a data checkout (the area
/// manifest and the root key) is `tools/p_integral_area.py stage` of
/// vakint-external-data.
pub(crate) fn export(
    cache_root: &Path,
    output: &Path,
    program: &OrderedCacheProgram,
    producer: &serde_json::Value,
    workload: &[serde_json::Value],
) -> Result<Option<ShippedOrderedCacheProgram>, VakintError> {
    let failure = |detail: String| -> VakintError { ordered_cache_failure(detail).into() };
    if program.data_version == INJECTED {
        return Ok(None);
    }
    let name = program_directory_name(program);
    let source = cache_root.join(&name);
    let manifest_path = source.join(RUSTRED_MANIFEST);
    if !manifest_path.is_file() {
        return Ok(None);
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
    let scope = label_scope(&program.label)
        .ok_or_else(|| failure(format!("unexpected p program label {}", program.label)))?;
    let area = output.join(AREA_DIRECTORY);
    let directory = format!("{CACHE_DIRECTORY}/{name}");
    let mut files = manifest.shards.clone();
    let tail_bytes = match fs::metadata(source.join(TAIL_FILE)) {
        Ok(metadata) => {
            files.push(TAIL_FILE.to_owned());
            metadata.len()
        }
        Err(_) => 0,
    };
    files.push(RUSTRED_MANIFEST.to_owned());
    for file in &files {
        let to = area.join(&directory).join(file);
        if let Some(parent) = to.parent() {
            fs::create_dir_all(parent)
                .map_err(|error| failure(format!("{}: {error}", parent.display())))?;
        }
        fs::copy(source.join(file), &to)
            .map_err(|error| failure(format!("{name}/{file}: {error}")))?;
    }
    let entry = PCacheProgram {
        label: program.label.clone(),
        scope,
        directory,
        identity: program.binding.family_fingerprint.clone(),
        p_versions: rustred::propagator::p_versions()
            .into_iter()
            .map(|(name, value)| (name.to_owned(), value))
            .collect(),
        stamp: manifest.stamp.clone(),
        adapter: manifest.adapter.clone(),
        records: manifest
            .shard_rows
            .iter()
            .map(|(_, records, _)| records)
            .sum(),
        shards: manifest.shard_rows.len(),
        bytes: manifest.shard_rows.iter().map(|(_, _, bytes)| bytes).sum(),
        tail_bytes,
        complete: true,
        trim: None,
    };
    let metadata = PCacheMetadata {
        schema: METADATA_SCHEMA.to_owned(),
        producer: producer.clone(),
        workload: workload
            .iter()
            .filter(|row| {
                row["loop_count"]
                    .as_u64()
                    .is_none_or(|count| (scope[0] as u64..=scope[1] as u64).contains(&count))
            })
            .cloned()
            .collect(),
        programs: vec![entry.clone()],
    };
    let path = area.join(METADATA_FILE);
    let mut text =
        serde_json::to_string_pretty(&metadata).map_err(|error| failure(error.to_string()))?;
    text.push('\n');
    fs::write(&path, text).map_err(|error| failure(format!("{}: {error}", path.display())))?;
    Ok(Some(entry.shipped_entry()))
}
