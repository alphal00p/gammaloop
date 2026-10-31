//! Persistent Ordered step cache.
//!
//! RustRed's ordered scheduler expands every operational identity with one
//! checked step: the ordered children with their exact step coefficients, the
//! emitted terminal coefficients and the step's nonzero conditions. That step
//! is a deterministic function of the identity and the immutable program, so
//! RustRed's `OrderedStepCache` records it once and replays it on every later
//! expansion of the same identity. A replayed step is identical to the checked
//! one, so the scheduler performs the same merges, cancellations and guard
//! interning as without a cache, in every weight arithmetic and every scheduler
//! pass. The cache never stores frontier or terminal weights, zero decisions or
//! final answers.
//!
//! Vakint keeps one cache per reduction program (the closed one- to three-loop
//! artifacts, the four four-loop candidate programs and the five-loop saved
//! application). [`OrderedCacheSettings`] selects how evaluations use them:
//! off, read-only or read-write (live population), in memory or backed by a
//! directory, with a resident-byte ceiling and an autosave policy.
//!
//! Validity is decided by versions and cheap identities only, never by digests
//! of engine sources. RustRed's binding stamp covers its cache schema version,
//! its hand-bumped step-semantics version and cheap identity fields (adapter,
//! family fingerprint, coefficient variables, unresolved policy). Vakint adds
//! the identity of the program's rule tables, read from digests the loader has
//! already verified (see [`OrderedCacheProgram::rule_tables`]); no data file is
//! hashed for the cache. External-data checkouts change their tables without
//! changing their release version, so the release version is reported but
//! never part of the key. A cache directory holds one subdirectory per
//! program, rule tables and stamp, `<loops>L-<program>.<tables>.<stamp>`: new
//! rule tables (a data update or a local override) or a RustRed semantics
//! bump use a new subdirectory, and RustRed rejects a file whose stamp
//! differs. A data update that leaves a program's tables unchanged keeps its
//! records.
//!
//! The depth-first completed-result caches (prepared bundles, running caches
//! and their transport) were removed together with depth-first reduction: a
//! completed terminal image replaces a frontier key and would skip every
//! grouping and cancellation below it.
use std::{
    collections::{BTreeMap, BTreeSet},
    fmt, fs,
    path::{Path, PathBuf},
    str::FromStr,
    sync::{Arc, Mutex, PoisonError},
    time::Instant,
};

use rustred::reduction::{OrderedStepCache, OrderedStepCacheBinding, OrderedStepCacheOptions};
use sha2::{Digest, Sha256};
use symbolica::atom::AtomView;

use super::{
    RustRedEvaluationError, RustRedEvaluationOptions, RustRedMasterMode, data::RustRedData,
};
use crate::{EvaluationMethod, EvaluationOrder, Vakint, VakintError, VakintSettings};

pub use rustred::reduction::{
    OrderedStepCacheSaveReport as RustRedOrderedStepCacheSaveReport,
    OrderedStepCacheShardStatistics as RustRedOrderedStepCacheShardStatistics,
    OrderedStepCacheStatistics as RustRedOrderedStepCacheStatistics,
};

/// Default ceiling on the estimated resident bytes of each program's Ordered
/// step cache. The complete five-loop zero-momentum and quadratic records are
/// estimated at about 2.5 GB; records beyond the ceiling are recomputed.
pub const DEFAULT_ORDERED_CACHE_MAX_RESIDENT_BYTES: usize = 8 << 30;

/// How evaluations use the Ordered step cache.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum OrderedCacheMode {
    /// No cache: every expanded identity runs its checked step.
    Off,
    /// Replay the records of the cache directory; record nothing new.
    ReadOnly,
    /// Replay existing records and record the step of every other expanded
    /// identity (live population).
    #[default]
    ReadWrite,
}

impl OrderedCacheMode {
    /// The canonical text: `off`, `read-only` or `read-write`.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::ReadOnly => "read-only",
            Self::ReadWrite => "read-write",
        }
    }
}

impl fmt::Display for OrderedCacheMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for OrderedCacheMode {
    type Err = VakintError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().replace('_', "-").as_str() {
            "off" => Ok(Self::Off),
            "read-only" | "readonly" => Ok(Self::ReadOnly),
            "read-write" | "readwrite" => Ok(Self::ReadWrite),
            _ => Err(invalid_setting(format!(
                "unknown Ordered cache mode {value:?}; expected 'off', 'read-only' or 'read-write'"
            ))),
        }
    }
}

/// When recorded steps of a directory-backed read-write cache are written.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub enum OrderedCacheAutosave {
    /// Only [`Vakint::save_ordered_cache`] writes records.
    Never,
    /// At the end of every public evaluation call that recorded new steps,
    /// including a failed one (its records are valid steps).
    #[default]
    AfterEvaluation,
    /// When the last handle of the engine (all its clones) is dropped.
    OnDrop,
}

impl OrderedCacheAutosave {
    /// The canonical text: `never`, `after-evaluation` or `on-drop`.
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Never => "never",
            Self::AfterEvaluation => "after-evaluation",
            Self::OnDrop => "on-drop",
        }
    }
}

impl fmt::Display for OrderedCacheAutosave {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for OrderedCacheAutosave {
    type Err = VakintError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().replace('_', "-").as_str() {
            "never" => Ok(Self::Never),
            "after-evaluation" => Ok(Self::AfterEvaluation),
            "on-drop" => Ok(Self::OnDrop),
            _ => Err(invalid_setting(format!(
                "unknown Ordered cache autosave policy {value:?}; expected 'never', 'after-evaluation' or 'on-drop'"
            ))),
        }
    }
}

/// Use of the persistent Ordered step cache by the evaluations of one engine.
///
/// The default records steps in memory for the lifetime of the engine (and its
/// clones), so later channels and requests replay the steps of earlier ones.
/// With a `directory`, caches are read from it lazily (one shard file at a
/// time) and, in read-write mode, written back according to `autosave`. A
/// directory has a single writer: do not let two engines or processes write
/// the same directory concurrently (records of one of them would be lost,
/// never mixed). Changing these settings between calls reopens the caches;
/// unsaved records of the previous configuration are saved first unless its
/// autosave policy is [`OrderedCacheAutosave::Never`], in which case they are
/// dropped. Results never depend on these settings.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct OrderedCacheSettings {
    pub mode: OrderedCacheMode,
    /// Root directory; `None` keeps records in memory only. Read-only mode
    /// replays saved records, so without a directory it uses no cache.
    pub directory: Option<PathBuf>,
    /// Ceiling on the estimated resident bytes of each program's cache;
    /// records and shards beyond it are not kept (their steps are
    /// recomputed). `None` is bounded only by
    /// [`VakintSettings::rustred_resources`].
    pub max_resident_bytes: Option<usize>,
    pub autosave: OrderedCacheAutosave,
}

impl Default for OrderedCacheSettings {
    fn default() -> Self {
        Self {
            mode: OrderedCacheMode::ReadWrite,
            directory: None,
            max_resident_bytes: Some(DEFAULT_ORDERED_CACHE_MAX_RESIDENT_BYTES),
            autosave: OrderedCacheAutosave::AfterEvaluation,
        }
    }
}

impl OrderedCacheSettings {
    /// No Ordered step cache.
    pub fn off() -> Self {
        Self {
            mode: OrderedCacheMode::Off,
            ..Self::default()
        }
    }

    /// A read-write cache in `directory`, saved after every evaluation.
    pub fn read_write(directory: impl Into<PathBuf>) -> Self {
        Self {
            directory: Some(directory.into()),
            ..Self::default()
        }
    }

    /// Replay the records saved in `directory` without recording new ones.
    pub fn read_only(directory: impl Into<PathBuf>) -> Self {
        Self {
            mode: OrderedCacheMode::ReadOnly,
            directory: Some(directory.into()),
            ..Self::default()
        }
    }
}

/// Counters of one program's attached cache.
#[derive(Clone, Debug, PartialEq)]
pub struct OrderedCacheProgramStatistics {
    pub loop_count: usize,
    /// Program label, e.g. `3L-unit-mass-vacuum-k6` or `4L-four-loop-h`.
    pub program: String,
    /// Version of the data the program was loaded from (file-name safe;
    /// reported only, not part of the cache identity).
    pub data_version: String,
    /// Hexadecimal identity of the program's rule tables (part of the cache
    /// identity, see the module documentation).
    pub rule_tables: String,
    /// Hexadecimal binding stamp of the cache.
    pub stamp: String,
    /// The program's cache directory, if the cache is backed by one.
    pub directory: Option<PathBuf>,
    pub mode: OrderedCacheMode,
    pub autosave: OrderedCacheAutosave,
    /// Cumulative native counters (hits, misses, inserts, records, bytes,
    /// per-shard states) since the cache was attached.
    pub native: RustRedOrderedStepCacheStatistics,
}

/// Counters of every attached cache of an engine.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OrderedCacheStatistics {
    /// One entry per attached program cache, sorted by loop count and label.
    pub programs: Vec<OrderedCacheProgramStatistics>,
    pub records: usize,
    pub hits: u64,
    pub misses: u64,
    pub inserts: u64,
    pub resident_bytes: usize,
    pub file_bytes: u64,
    /// The last failed automatic save, if any (automatic saves never fail an
    /// evaluation; an explicit save reports its error).
    pub last_autosave_error: Option<String>,
}

/// Result of one program's save.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct OrderedCacheProgramSaveReport {
    pub loop_count: usize,
    pub program: String,
    pub native: RustRedOrderedStepCacheSaveReport,
}

/// Result of [`Vakint::save_ordered_cache`]: only directory-backed read-write
/// caches with unsaved records are written.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct OrderedCacheSaveReport {
    pub programs: Vec<OrderedCacheProgramSaveReport>,
    pub shards_written: usize,
    pub records_written: usize,
    pub bytes_written: u64,
}

/// Result of [`Vakint::warm_ordered_cache`].
#[derive(Clone, Debug, Default, PartialEq)]
pub struct OrderedCacheWarmReport {
    pub integrals: usize,
    pub seconds: f64,
    /// Replayed and recorded steps of the warming evaluations.
    pub hits: u64,
    pub misses: u64,
    pub inserts: u64,
    /// Present when the warmed caches have a directory.
    pub save: Option<OrderedCacheSaveReport>,
    pub statistics: OrderedCacheStatistics,
}

/// One program whose Ordered steps a cache records.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct OrderedCacheProgram {
    pub loop_count: usize,
    /// File-name safe label, unique per program within a loop count.
    pub label: String,
    /// File-name safe version of the data the program was loaded from
    /// (reported only: checkouts do not bump it when their tables change).
    pub data_version: String,
    /// Identity of the program's rule tables, from digests verified on load:
    /// the sha256 of a closed program's `rules.rrbin`, a digest of the pack
    /// manifest's declared sha256 of a four-loop program's inputs, and the
    /// saved-application identity of a five-loop package (its
    /// `application.json`, which pins every payload, with the pack's families
    /// and masters).
    pub rule_tables: [u8; 32],
    pub binding: OrderedStepCacheBinding,
}

/// Registry key: program label, rule tables and binding stamp.
type ProgramKey = (String, [u8; 32], [u8; 32]);

/// The rule-table identity of a program read from several data files: a
/// digest of their paths and declared (verified) sha256 values, in order.
pub(super) fn declared_rule_tables(
    files: impl IntoIterator<Item = (impl AsRef<str>, impl AsRef<str>)>,
) -> [u8; 32] {
    let mut digest = Sha256::new();
    digest.update(b"vakint-ordered-cache-rule-tables-v1\0");
    for (path, sha256) in files {
        for field in [path.as_ref(), sha256.as_ref()] {
            digest.update((field.len() as u64).to_le_bytes());
            digest.update(field.as_bytes());
        }
    }
    digest.finalize().into()
}

/// At most 64 ASCII letters, digits, `-`, `.` and `_` (others become `-`).
fn file_name_safe(text: &str, lowercase: bool) -> String {
    text.chars()
        .map(|character| {
            if character.is_ascii_alphanumeric() || matches!(character, '-' | '.' | '_') {
                if lowercase {
                    character.to_ascii_lowercase()
                } else {
                    character
                }
            } else {
                '-'
            }
        })
        .take(64)
        .collect()
}

impl OrderedCacheProgram {
    pub(crate) fn new(
        loop_count: usize,
        name: &str,
        data_version: &str,
        rule_tables: [u8; 32],
        binding: OrderedStepCacheBinding,
    ) -> Self {
        Self {
            loop_count,
            label: format!(
                "{loop_count}L-{}",
                file_name_safe(&name.replace(['.', '_'], "-"), true)
            ),
            data_version: file_name_safe(data_version, false),
            rule_tables,
            binding,
        }
    }

    /// Registry key: one cache per program, rule tables and binding stamp.
    fn key(&self) -> ProgramKey {
        (self.label.clone(), self.rule_tables, self.binding.stamp())
    }
}

/// How a program cache was opened. Changing any of these reopens it.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Configuration {
    directory: Option<PathBuf>,
    populate: bool,
    max_resident_bytes: Option<usize>,
}

struct Attached {
    loop_count: usize,
    label: String,
    data_version: String,
    rule_tables: [u8; 32],
    stamp: [u8; 32],
    configuration: Configuration,
    autosave: OrderedCacheAutosave,
    cache: Arc<OrderedStepCache>,
}

impl Attached {
    fn mode(&self) -> OrderedCacheMode {
        if self.configuration.populate {
            OrderedCacheMode::ReadWrite
        } else {
            OrderedCacheMode::ReadOnly
        }
    }

    /// Whether a save would write anything.
    fn has_unsaved_records(&self) -> bool {
        self.configuration.populate
            && self.configuration.directory.is_some()
            && self
                .cache
                .statistics()
                .shards
                .iter()
                .any(|shard| shard.dirty)
    }

    fn save(&self) -> Result<Option<OrderedCacheProgramSaveReport>, RustRedEvaluationError> {
        if !self.has_unsaved_records() {
            return Ok(None);
        }
        let native = self
            .cache
            .save()
            .map_err(|error| ordered_cache_failure(format!("{}: {error}", self.label)))?;
        Ok(Some(OrderedCacheProgramSaveReport {
            loop_count: self.loop_count,
            program: self.label.clone(),
            native,
        }))
    }

    fn statistics(&self) -> OrderedCacheProgramStatistics {
        OrderedCacheProgramStatistics {
            loop_count: self.loop_count,
            program: self.label.clone(),
            data_version: self.data_version.clone(),
            rule_tables: hex(&self.rule_tables),
            stamp: hex(&self.stamp),
            directory: self.configuration.directory.clone(),
            mode: self.mode(),
            autosave: self.autosave,
            native: self.cache.statistics(),
        }
    }
}

/// The Ordered step caches of one engine, shared by its clones. Caches are
/// keyed by program, rule tables and binding stamp, so one program has one
/// cache.
#[derive(Default)]
pub(crate) struct OrderedCacheRegistry {
    attached: Mutex<BTreeMap<ProgramKey, Attached>>,
    last_autosave_error: Mutex<Option<String>>,
}

impl fmt::Debug for OrderedCacheRegistry {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("OrderedCacheRegistry")
            .field("attached", &self.lock().len())
            .finish()
    }
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

pub(super) fn ordered_cache_failure(detail: impl Into<String>) -> RustRedEvaluationError {
    RustRedEvaluationError::OrderedStepCache {
        detail: detail.into(),
    }
}

fn invalid_setting(detail: String) -> VakintError {
    ordered_cache_failure(detail).into()
}

#[cfg(test)]
#[path = "ordered_cache_tests.rs"]
mod tests;

impl OrderedCacheRegistry {
    fn lock(&self) -> std::sync::MutexGuard<'_, BTreeMap<ProgramKey, Attached>> {
        self.attached.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// The program's cache under `settings`, opened (or reopened after a
    /// settings change) on first use. `None` when the cache is off.
    pub(crate) fn cache(
        &self,
        settings: &OrderedCacheSettings,
        program: &OrderedCacheProgram,
    ) -> Result<Option<Arc<OrderedStepCache>>, RustRedEvaluationError> {
        // Read-only replay needs records saved in a directory.
        if settings.mode == OrderedCacheMode::Off
            || (settings.mode == OrderedCacheMode::ReadOnly && settings.directory.is_none())
        {
            return Ok(None);
        }
        let key = program.key();
        let stamp = key.2;
        let configuration = Configuration {
            directory: settings
                .directory
                .as_ref()
                .map(|root| program_directory(root, program, &stamp)),
            populate: settings.mode == OrderedCacheMode::ReadWrite,
            max_resident_bytes: settings.max_resident_bytes,
        };
        let mut attached = self.lock();
        if let Some(current) = attached.get_mut(&key) {
            if current.configuration == configuration {
                current.autosave = settings.autosave;
                current.data_version.clone_from(&program.data_version);
                return Ok(Some(Arc::clone(&current.cache)));
            }
            // An automatic save: a failure is reported, not raised.
            if current.autosave != OrderedCacheAutosave::Never
                && let Err(error) = current.save()
            {
                self.autosave_failed(error);
            }
        }
        let options = OrderedStepCacheOptions {
            populate: configuration.populate,
            max_resident_bytes: configuration.max_resident_bytes,
            ..OrderedStepCacheOptions::default()
        };
        let cache = match &configuration.directory {
            Some(directory) => OrderedStepCache::open(directory, program.binding.clone(), options)
                .map_err(|error| {
                    ordered_cache_failure(format!(
                        "{}: {error}; remove the directory to rebuild it",
                        program.label
                    ))
                })?,
            None => OrderedStepCache::in_memory(program.binding.clone(), options),
        };
        attached.insert(
            key,
            Attached {
                loop_count: program.loop_count,
                label: program.label.clone(),
                data_version: program.data_version.clone(),
                rule_tables: program.rule_tables,
                stamp,
                configuration,
                autosave: settings.autosave,
                cache: Arc::clone(&cache),
            },
        );
        Ok(Some(cache))
    }

    pub(crate) fn statistics(&self) -> OrderedCacheStatistics {
        let attached = self.lock();
        let mut programs: Vec<_> = attached.values().map(Attached::statistics).collect();
        drop(attached);
        programs.sort_by(|left, right| {
            (
                left.loop_count,
                &left.program,
                &left.rule_tables,
                &left.stamp,
            )
                .cmp(&(
                    right.loop_count,
                    &right.program,
                    &right.rule_tables,
                    &right.stamp,
                ))
        });
        let mut statistics = OrderedCacheStatistics {
            last_autosave_error: self
                .last_autosave_error
                .lock()
                .unwrap_or_else(PoisonError::into_inner)
                .clone(),
            ..Default::default()
        };
        for program in &programs {
            statistics.records += program.native.records;
            statistics.hits += program.native.hits;
            statistics.misses += program.native.misses;
            statistics.inserts += program.native.inserts;
            statistics.resident_bytes += program.native.resident_bytes;
            statistics.file_bytes += program.native.file_bytes;
        }
        statistics.programs = programs;
        statistics
    }

    /// Save every directory-backed read-write cache with unsaved records
    /// (optionally only those with the given autosave policy).
    fn save(
        &self,
        policy: Option<OrderedCacheAutosave>,
    ) -> Result<OrderedCacheSaveReport, RustRedEvaluationError> {
        let attached = self.lock();
        let mut report = OrderedCacheSaveReport::default();
        for program in attached.values() {
            if policy.is_some_and(|policy| program.autosave != policy) {
                continue;
            }
            if let Some(saved) = program.save()? {
                report.shards_written += saved.native.shards_written;
                report.records_written += saved.native.records_written;
                report.bytes_written += saved.native.bytes_written;
                report.programs.push(saved);
            }
        }
        Ok(report)
    }

    /// The automatic save at the end of a public evaluation call. A failure
    /// is logged and kept for the statistics; it never fails the evaluation.
    pub(crate) fn after_evaluation(&self) {
        if let Err(error) = self.save(Some(OrderedCacheAutosave::AfterEvaluation)) {
            self.autosave_failed(error);
        }
    }

    fn autosave_failed(&self, error: RustRedEvaluationError) {
        log::warn!("Vakint could not save its Ordered step cache: {error}");
        *self
            .last_autosave_error
            .lock()
            .unwrap_or_else(PoisonError::into_inner) = Some(error.to_string());
    }

    fn clear(&self) {
        self.lock().clear();
    }

    /// Detach the caches of every program not in `programs` (superseded by a
    /// data or five-loop application change), so their resident records are
    /// released. Records of a directory-backed cache are saved first unless
    /// its autosave policy is `Never`; a failed save is reported, not raised.
    pub(crate) fn detach_except(&self, programs: &[OrderedCacheProgram]) {
        let keep: BTreeSet<ProgramKey> = programs.iter().map(OrderedCacheProgram::key).collect();
        let detached: Vec<Attached> = {
            let mut attached = self.lock();
            let superseded: Vec<ProgramKey> = attached
                .keys()
                .filter(|key| !keep.contains(*key))
                .cloned()
                .collect();
            superseded
                .iter()
                .filter_map(|key| attached.remove(key))
                .collect()
        };
        for program in detached {
            if program.autosave != OrderedCacheAutosave::Never
                && let Err(error) = program.save()
            {
                self.autosave_failed(error);
            }
        }
    }
}

impl Drop for OrderedCacheRegistry {
    fn drop(&mut self) {
        let attached = self
            .attached
            .get_mut()
            .unwrap_or_else(PoisonError::into_inner);
        for program in attached.values() {
            if program.autosave == OrderedCacheAutosave::Never {
                continue;
            }
            if let Err(error) = program.save() {
                log::warn!("Vakint could not save its Ordered step cache on drop: {error}");
            }
        }
    }
}

/// `<root>/<label>.<first 16 hex digits of the rule tables>.<first 16 hex
/// digits of the stamp>`.
fn program_directory(root: &Path, program: &OrderedCacheProgram, stamp: &[u8; 32]) -> PathBuf {
    root.join(format!(
        "{}.{}.{}",
        program.label,
        hex(&program.rule_tables[..8]),
        hex(&stamp[..8])
    ))
}

impl Vakint {
    /// Attach the Ordered step caches of every loaded program now, as
    /// evaluations would on first use, and return their statistics. Opening
    /// a directory reads and checks only its manifests; shard files are read
    /// at their first lookup. A five-loop saved application is covered once
    /// it is loaded. Fails when a directory holds a cache with another
    /// binding (RustRed's stamp check: another cache schema or step-semantics
    /// version, family, variables or policy) or cannot be read.
    pub fn load_ordered_cache(
        &self,
        settings: &VakintSettings,
    ) -> Result<OrderedCacheStatistics, VakintError> {
        for program in self.ordered_cache_programs()? {
            self.ordered_cache
                .cache(&settings.ordered_cache, &program)?;
        }
        Ok(self.ordered_cache_statistics())
    }

    /// Write the unsaved records of every directory-backed read-write cache.
    /// Shard files are content-addressed and the manifest is replaced
    /// atomically, so an interrupted save keeps the previous cache valid.
    pub fn save_ordered_cache(&self) -> Result<OrderedCacheSaveReport, VakintError> {
        Ok(self.ordered_cache.save(None)?)
    }

    /// Hits, misses, inserts, records and bytes of every attached cache.
    pub fn ordered_cache_statistics(&self) -> OrderedCacheStatistics {
        self.ordered_cache.statistics()
    }

    /// Detach every cache and drop its resident records. Files are kept; the
    /// next evaluation reopens its directory. Unsaved records are discarded.
    pub fn clear_ordered_cache(&self) {
        self.ordered_cache.clear();
    }

    /// Populate the read-write caches by reducing every integral expression
    /// with RustRed (retained masters; results are discarded), then save
    /// them when they have a directory. Each expression is one evaluation, so
    /// its own cancellations apply. Read-only and disabled cache modes are
    /// treated as read-write here. The first failing evaluation is returned
    /// as an error; records of earlier ones are kept.
    pub fn warm_ordered_cache<'a>(
        &self,
        settings: &VakintSettings,
        integrals: impl IntoIterator<Item = AtomView<'a>>,
    ) -> Result<OrderedCacheWarmReport, VakintError> {
        let started = Instant::now();
        let mut settings = settings.clone();
        settings.ordered_cache.mode = OrderedCacheMode::ReadWrite;
        settings.evaluation_order =
            EvaluationOrder(vec![EvaluationMethod::RustRed(RustRedEvaluationOptions {
                master_mode: RustRedMasterMode::Retain,
            })]);
        let before = self.ordered_cache_statistics();
        let mut report = OrderedCacheWarmReport::default();
        for integral in integrals {
            // Only the recorded steps matter; the retained value is dropped.
            drop(self.evaluate_observed(&settings, integral, None, None)?);
            report.integrals += 1;
        }
        report.save = settings
            .ordered_cache
            .directory
            .is_some()
            .then(|| self.save_ordered_cache())
            .transpose()?;
        report.statistics = self.ordered_cache_statistics();
        report.hits = report.statistics.hits.saturating_sub(before.hits);
        report.misses = report.statistics.misses.saturating_sub(before.misses);
        report.inserts = report.statistics.inserts.saturating_sub(before.inserts);
        report.seconds = started.elapsed().as_secs_f64();
        Ok(report)
    }

    /// Run the automatic save after a public evaluation call.
    pub(crate) fn autosave_ordered_cache<T>(&self, result: T) -> T {
        self.ordered_cache.after_evaluation();
        result
    }

    /// The Ordered step cache of one program under these settings.
    pub(crate) fn ordered_step_cache(
        &self,
        settings: &VakintSettings,
        program: &OrderedCacheProgram,
    ) -> Result<Option<Arc<OrderedStepCache>>, RustRedEvaluationError> {
        self.ordered_cache.cache(&settings.ordered_cache, program)
    }

    /// Detach the caches of programs this engine no longer has, after its
    /// data were replaced. Nothing is detached if the programs cannot be
    /// listed.
    pub(super) fn detach_superseded_ordered_caches(&self) {
        if let Ok(programs) = self.ordered_cache_programs() {
            self.ordered_cache.detach_except(&programs);
        }
    }

    /// Every program currently loaded in this engine.
    fn ordered_cache_programs(&self) -> Result<Vec<OrderedCacheProgram>, VakintError> {
        let data = self.rustred_data()?;
        let mut programs = data.ordered_cache_programs()?;
        if let Ok(application) = self.loaded_application() {
            programs.push(application.ordered_cache_program().clone());
        }
        Ok(programs)
    }
}

pub(super) fn operation(error: rustred_app::AppError) -> RustRedEvaluationError {
    RustRedEvaluationError::CacheOperation { source: error }
}

fn failure(error: impl std::fmt::Display) -> VakintError {
    operation(rustred_app::AppError::new(
        rustred_app::AppErrorKind::Input,
        error.to_string(),
    ))
    .into()
}

/// An application identity deliberately excludes cache inventory, release names
/// and numerical master coefficients: none changes a raw Route/Apply record.
pub(super) fn saved_application_digest(
    data: &RustRedData,
    path: &Path,
) -> Result<[u8; 32], VakintError> {
    let manifest = if path.is_dir() {
        path.join("application.json")
    } else {
        path.to_owned()
    };
    let application = fs::read(manifest).map_err(failure)?;
    let mut digest = Sha256::new();
    digest.update(b"vakint-saved-application-cache-v1\0");
    digest.update((application.len() as u64).to_le_bytes());
    digest.update(application);
    if let Some(pack) = data.loops.get(&5) {
        fn field(digest: &mut Sha256, bytes: &[u8]) {
            digest.update((bytes.len() as u64).to_le_bytes());
            digest.update(bytes);
        }
        field(&mut digest, pack.manifest.dimension.as_bytes());
        field(&mut digest, pack.manifest.normalization.as_bytes());
        let mut families = pack
            .manifest
            .families
            .iter()
            .map(|family| {
                family["fingerprint"]
                    .as_str()
                    .ok_or_else(|| failure("family fingerprint is absent"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        families.sort_unstable();
        digest.update((families.len() as u64).to_le_bytes());
        for family in families {
            field(&mut digest, family.as_bytes());
        }
        let mut masters = pack
            .manifest
            .masters
            .iter()
            .map(|master| {
                let family = master["family_sha256"]
                    .as_str()
                    .ok_or_else(|| failure("master family identity is absent"))?;
                let powers = master["powers"]
                    .as_array()
                    .ok_or_else(|| failure("master powers are absent"))?
                    .iter()
                    .map(|power| {
                        power
                            .as_i64()
                            .ok_or_else(|| failure("invalid master power"))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                Ok((family, powers))
            })
            .collect::<Result<Vec<_>, VakintError>>()?;
        masters.sort_unstable();
        digest.update((masters.len() as u64).to_le_bytes());
        for (family, powers) in masters {
            field(&mut digest, family.as_bytes());
            digest.update((powers.len() as u64).to_le_bytes());
            for power in powers {
                digest.update(power.to_le_bytes());
            }
        }
    }
    Ok(digest.finalize().into())
}
