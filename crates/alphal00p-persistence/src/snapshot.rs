use crate::{OutputPath, OutputPathError};
use same_file::Handle;
use serde::{Deserialize, Serialize};
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::fs::{self, File, TryLockError};
use std::io::{Error as IoError, ErrorKind, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::str::FromStr;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, SystemTimeError, UNIX_EPOCH};
use thiserror::Error;

pub(crate) const CREATE_ATTEMPTS: usize = 128;

/// Directory containing immutable payload generations. See the [crate example](crate).
pub const GENERATIONS_DIRECTORY: &str = "generations";
/// Persistent advisory writer lock. See [SnapshotWriter::acquire].
pub const SNAPSHOT_LOCK_FILE: &str = ".snapshot.lock";
static NEXT_GENERATION: AtomicU64 = AtomicU64::new(0);

/// A single component in the format `gen-<8 hex pid>-<32 hex nanoseconds>-<16 hex counter>`.
///
/// Hexadecimal digits are lowercase. Only the writer generates new identities; parsing and Serde
/// decoding validate references read from manifests. See the [crate example](crate).
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(try_from = "String")]
pub struct GenerationId(String);

impl GenerationId {
    /// Borrow the validated component, as shown in the [crate example](crate).
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub(crate) fn generate() -> Result<Self, SystemTimeError> {
        let nanos = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
        let counter = NEXT_GENERATION.fetch_add(1, Ordering::Relaxed);
        Ok(Self(format!(
            "gen-{:08x}-{nanos:032x}-{counter:016x}",
            process::id()
        )))
    }
}

impl Display for GenerationId {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> FmtResult {
        self.0.fmt(formatter)
    }
}

impl TryFrom<String> for GenerationId {
    type Error = SnapshotError;

    fn try_from(value: String) -> Result<Self, Self::Error> {
        let mut parts = value.split('-');
        let valid = parts.next() == Some("gen")
            && [8, 32, 16].into_iter().all(|length| {
                parts.next().is_some_and(|part| {
                    part.len() == length
                        && part
                            .bytes()
                            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
                })
            })
            && parts.next().is_none();
        if !valid {
            return Err(SnapshotError::InvalidGeneration { value });
        }
        Ok(Self(value))
    }
}

impl FromStr for GenerationId {
    type Err = SnapshotError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        Self::try_from(value.to_owned())
    }
}

/// A validated, caller-provided commit-pointer filename. See the [crate example](crate).
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct ManifestName<'name>(&'name str);

impl<'name> ManifestName<'name> {
    /// Validate one filename, excluding reserved snapshot entries. See the [crate example](crate).
    ///
    /// # Errors
    /// Rejects invalid components and reserved lock and generation entries.
    pub fn new(name: &'name str) -> Result<Self, SnapshotError> {
        OutputPath::validate_component(name)?;
        if matches!(name, SNAPSHOT_LOCK_FILE | GENERATIONS_DIRECTORY) {
            return Err(SnapshotError::ReservedManifest {
                name: name.to_owned(),
            });
        }
        Ok(Self(name))
    }

    /// Borrow the exact commit-pointer filename. See the [crate example](crate).
    #[must_use]
    pub fn as_str(self) -> &'name str {
        self.0
    }

    fn is_temporary(self, name: &str) -> bool {
        name.strip_prefix('.')
            .and_then(|name| name.strip_prefix(self.0))
            .and_then(|name| name.strip_prefix('.'))
            .and_then(|name| name.strip_suffix(".tmp"))
            .is_some_and(|id| id.parse::<GenerationId>().is_ok())
    }
}

impl Display for ManifestName<'_> {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> FmtResult {
        self.0.fmt(formatter)
    }
}

/// Publication failures, retaining I/O sources and distinguishing immediate writer contention.
///
/// Only [`Self::Committed`] follows the publication commit point. The [crate example](crate) shows
/// propagation; [`SnapshotWriter::retire_root`] demonstrates explicit retirement.
#[derive(Debug, Error)]
pub enum SnapshotError {
    /// A filesystem operation failed, including unsupported locking or replacement.
    #[error("cannot {operation} at {}", path.display())]
    Io {
        /// The attempted operation.
        operation: &'static str,
        /// The affected entry.
        path: PathBuf,
        /// The original filesystem error.
        #[source]
        source: IoError,
    },
    /// Another cooperating writer holds the lock; retry after it releases ownership.
    #[error("snapshot writer is busy at {}", path.display())]
    Busy {
        /// The persistent lock entry.
        path: PathBuf,
    },
    /// The acquired lock no longer belongs to the root at the requested path.
    #[error(
        "snapshot root disappeared or its writer lock was replaced at {}; retry the operation",
        path.display()
    )]
    StaleLock {
        /// The lock path that must still name the acquired file.
        path: PathBuf,
    },
    /// A manifest supplied a component outside the internally generated format.
    #[error("invalid snapshot generation {value:?}")]
    InvalidGeneration {
        /// The rejected component.
        value: String,
    },
    /// A caller-selected manifest would replace a reserved snapshot entry.
    #[error("snapshot manifest name {name:?} is reserved")]
    ReservedManifest {
        /// The rejected filename.
        name: String,
    },
    /// Shared path validation or filesystem synchronization failed before commit.
    #[error(transparent)]
    Path(#[from] OutputPathError),
    /// The new pointer is visible, but its directory synchronization failed.
    ///
    /// Retain this generation and prior generations; do not run restart cleanup or assume rollback.
    #[error("snapshot {} committed but durability is unconfirmed", generation.id())]
    Committed {
        /// The generation now selected by the visible manifest.
        generation: PublishedGeneration,
        /// The failed directory synchronization.
        #[source]
        source: OutputPathError,
    },
    /// The root was detached, but synchronizing its parent failed.
    #[error("snapshot root retired to {}; parent durability is unconfirmed", path.display())]
    Retired {
        /// The detached sibling retained for the next explicit cleanup.
        path: PathBuf,
        /// The failed parent synchronization.
        #[source]
        source: OutputPathError,
    },
    /// The detached tree was removed, but synchronizing its parent failed.
    #[error(
        "retired snapshot tree removed at {}; parent durability is unconfirmed",
        path.display()
    )]
    Removed {
        /// The removed detached sibling, which must not be treated as a retained backup.
        path: PathBuf,
        /// The failed parent synchronization.
        #[source]
        source: OutputPathError,
    },
    /// The system clock cannot supply a timestamp relative to the Unix epoch.
    #[error("cannot generate a snapshot identity")]
    Clock(#[from] SystemTimeError),
}

impl SnapshotError {
    /// Whether a cooperating writer currently owns the lock. See [`SnapshotWriter::acquire`].
    #[must_use]
    pub fn is_busy(&self) -> bool {
        matches!(self, Self::Busy { .. })
    }

    fn io(operation: &'static str, path: &Path, source: IoError) -> Self {
        Self::Io {
            operation,
            path: path.to_owned(),
            source,
        }
    }
}

/// Exclusive owner of a logical snapshot root. Dropping it releases the OS lock.
///
/// An active snapshot root retains its lock inode. Explicit [`Self::retire_root`] detaches the
/// entire root before removal deletes that detached inode. See the [crate example](crate).
#[derive(Debug)]
pub struct SnapshotWriter {
    root: OutputPath,
    lock: File,
}

impl SnapshotWriter {
    /// Create an authorized root if needed and acquire its persistent `.snapshot.lock`.
    ///
    /// The returned root is absolute, with root symlinks resolved. This performs no application
    /// authorization: callers must check their destination policy before calling it and recheck
    /// under the acquired lock. See the [crate example](crate).
    /// Under the acquired lock, removes stale temporary files belonging to the supplied manifests.
    /// On Unix, synchronizes directories created here and the root. Parent synchronization is
    /// attempted on every acquisition, tolerating permission and unsupported-operation errors
    /// so roots work beneath search-only parents. Other parent I/O errors remain fatal on retry.
    /// Neither payload generations nor manifests are removed.
    ///
    /// # Errors
    /// Returns [`SnapshotError::Busy`] immediately on contention. Invalid entry types, symlinks,
    /// I/O errors, and unsupported filesystem operations are errors. A failed acquisition does not
    /// publish a manifest, but can have removed abandoned manifest temporary files.
    /// [`SnapshotError::StaleLock`] rejects an opener delayed across retirement or replacement;
    /// it never cleans temporary files in the replacement root.
    pub fn acquire(root: &Path, manifests: &[ManifestName<'_>]) -> Result<Self, SnapshotError> {
        let root = OutputPath::new(root)?;
        root.ensure_durable_directory()?;
        let lock_path = root.join(SNAPSHOT_LOCK_FILE)?;
        lock_path.check_regular_or_missing()?;
        let lock = File::options()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(lock_path.as_path())
            .map_err(|source| SnapshotError::io("open writer lock", lock_path.as_path(), source))?;
        Self::finish_acquire(root, lock, manifests)
    }

    fn finish_acquire(
        root: OutputPath,
        lock: File,
        manifests: &[ManifestName<'_>],
    ) -> Result<Self, SnapshotError> {
        let lock_path = root.as_path().join(SNAPSHOT_LOCK_FILE);
        Self::lock_file(&lock, &lock_path)?;
        let writer = Self { root, lock };
        writer.verify_current_lock()?;
        writer.remove_manifest_temporaries(manifests)?;
        Ok(writer)
    }

    fn lock_file(lock: &File, path: &Path) -> Result<(), SnapshotError> {
        match lock.try_lock() {
            Ok(()) => Ok(()),
            Err(TryLockError::WouldBlock) => Err(SnapshotError::Busy {
                path: path.to_owned(),
            }),
            Err(TryLockError::Error(source)) => {
                Err(SnapshotError::io("lock snapshot root", path, source))
            }
        }
    }

    fn verify_current_lock(&self) -> Result<(), SnapshotError> {
        let lock_path = self.root.join(SNAPSHOT_LOCK_FILE)?;
        let stale = || SnapshotError::StaleLock {
            path: lock_path.as_path().to_owned(),
        };
        if lock_path.metadata()?.is_none() {
            return Err(stale());
        }
        lock_path.check_regular_or_missing()?;
        Self::verify_lock(&self.lock, lock_path.as_path())
    }

    fn verify_lock(lock: &File, path: &Path) -> Result<(), SnapshotError> {
        let stale = || SnapshotError::StaleLock {
            path: path.to_owned(),
        };
        // Verification must follow try_lock: an opener can retain the old descriptor while
        // retirement holds the lock, then acquire it after the original path has been reused.
        let acquired = lock
            .try_clone()
            .and_then(Handle::from_file)
            .map_err(|source| SnapshotError::io("identify acquired writer lock", path, source))?;
        let current = Handle::from_path(path).map_err(|source| {
            if source.kind() == ErrorKind::NotFound {
                stale()
            } else {
                SnapshotError::io("identify current writer lock", path, source)
            }
        })?;
        if acquired != current {
            return Err(stale());
        }
        Ok(())
    }

    /// Try-lock existing nested writers and atomically detach this root for removal.
    ///
    /// The caller authorizes deletion of the entire tree and must not start new nested
    /// snapshot writers during cleanup. Existing nested locks are acquired without following
    /// symlinks.
    /// The root is renamed to a unique hidden sibling, named
    /// `.snapshot-retired-<root-name-hash>-<generation>`. All acquired locks remain held until
    /// [`RetiredSnapshotRoot::remove`] completes or the returned owner is dropped.
    /// A new writer can initialize the original path independently after detachment.
    /// A drop, crash, or failed removal leaves the sibling for [`Self::remove_retired_roots`]
    /// during the next explicit cleanup.
    ///
    /// ```
    /// use alphal00p_persistence::SnapshotWriter;
    /// use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// # let temporary = tempfile::tempdir()?;
    /// let root = temporary.path().join("authorized-snapshot");
    /// SnapshotWriter::remove_retired_roots(&root)?;
    /// let writer = SnapshotWriter::acquire(&root, &[])?;
    /// writer.retire_root()?.remove()?;
    /// assert!(!root.exists());
    /// # Ok(()) }
    /// ```
    ///
    /// # Errors
    /// Refuses busy nested writers, stale locks, and inspection or rename failures. A failed
    /// parent synchronization after rename leaves a detached sibling for the next cleanup.
    pub fn retire_root(mut self) -> Result<RetiredSnapshotRoot, SnapshotError> {
        self.verify_current_lock()?;
        self.root.check_directory()?;
        let nested_locks = self.acquire_nested_locks()?;
        self.verify_current_lock()?;
        let parent_path = self
            .root()
            .parent()
            .ok_or_else(|| OutputPathError::InvalidRoot(self.root().to_owned()))?;
        let parent = OutputPath::new(parent_path)?;
        let prefix = Self::retired_prefix(self.root())?;
        for attempt in 0..CREATE_ATTEMPTS {
            let name = format!("{prefix}{}", GenerationId::generate()?);
            let directory = parent.join(&name)?;
            if directory.metadata()?.is_some() {
                if attempt + 1 < CREATE_ATTEMPTS {
                    continue;
                }
                return Err(SnapshotError::io(
                    "select unused retired snapshot path",
                    directory.as_path(),
                    IoError::from(ErrorKind::AlreadyExists),
                ));
            }
            // The generated identity is private to this cleanup. As with checked output paths,
            // callers exclude hostile concurrent replacement of the destination.
            fs::rename(self.root(), directory.as_path())
                .map_err(|source| SnapshotError::io("retire snapshot root", self.root(), source))?;
            self.root = directory;
            let retired = RetiredSnapshotRoot {
                writer: self,
                _nested_locks: nested_locks,
            };
            parent
                .sync_directory()
                .map_err(|source| SnapshotError::Retired {
                    path: retired.writer.root().to_owned(),
                    source,
                })?;
            return Ok(retired);
        }
        unreachable!("the last reservation attempt returns its error")
    }

    fn retired_prefix(root: &Path) -> Result<String, SnapshotError> {
        let hash = OutputPath::filename_hash(root)?;
        Ok(format!(".snapshot-retired-{hash}-"))
    }

    /// Remove abandoned detached siblings belonging to this root during explicit cleanup.
    ///
    /// This never creates or removes the requested root. It also works when that root is absent;
    /// an unreadable parent is skipped when the root is missing.
    /// Only directories with this root's name hash and a valid generation suffix are inspected;
    /// symlinks, unrelated entries, and siblings with a busy root or nested writer are left alone.
    /// Each reclaimed sibling is removed while its root and existing nested locks remain held.
    /// Partial removals, including a missing lock file, can be retried on the next cleanup.
    /// See [`Self::retire_root`] for an example.
    ///
    /// # Errors
    /// Returns path-validation, lock, inspection, removal, or synchronization failures other
    /// than contention. Callers must authorize cleanup and exclude new nested snapshot writers.
    pub fn remove_retired_roots(root: &Path) -> Result<(), SnapshotError> {
        let root = OutputPath::new(root)?;
        let prefix = Self::retired_prefix(root.as_path())?;
        let parent_path = root
            .as_path()
            .parent()
            .ok_or_else(|| OutputPathError::InvalidRoot(root.as_path().to_owned()))?;
        let entries = match fs::read_dir(parent_path) {
            Ok(entries) => entries,
            Err(source) if source.kind() == ErrorKind::NotFound => return Ok(()),
            Err(source)
                if source.kind() == ErrorKind::PermissionDenied && root.metadata()?.is_none() =>
            {
                return Ok(());
            }
            Err(source) => {
                return Err(SnapshotError::io(
                    "list retired snapshots",
                    parent_path,
                    source,
                ));
            }
        };
        for entry in entries {
            let entry = entry.map_err(|source| {
                SnapshotError::io("inspect retired snapshot", parent_path, source)
            })?;
            let name = entry.file_name();
            let matches_root = name
                .to_str()
                .and_then(|name| name.strip_prefix(&prefix))
                .is_some_and(|id| id.parse::<GenerationId>().is_ok());
            if !matches_root {
                continue;
            }
            let path = entry.path();
            let kind = entry.file_type().map_err(|source| {
                SnapshotError::io("inspect retired snapshot type", &path, source)
            })?;
            if !kind.is_dir() {
                continue;
            }
            Self::remove_retired_root(&path)?;
        }
        Ok(())
    }

    fn remove_retired_root(path: &Path) -> Result<(), SnapshotError> {
        let writer = match Self::acquire(path, &[]) {
            Ok(writer) => writer,
            Err(SnapshotError::Busy { .. }) => return Ok(()),
            Err(error) => return Err(error),
        };
        // Acquisition resolves a user-selected root once. A sibling swapped for a symlink
        // after its directory type check must never authorize deletion of the link's target.
        if writer.root() != path {
            return Ok(());
        }
        let nested_locks = match writer.acquire_nested_locks() {
            Ok(locks) => locks,
            Err(SnapshotError::Busy { .. }) => return Ok(()),
            Err(error) => return Err(error),
        };
        RetiredSnapshotRoot {
            writer,
            _nested_locks: nested_locks,
        }
        .remove()
    }

    fn acquire_nested_locks(&self) -> Result<Vec<File>, SnapshotError> {
        // Inspect only existing regular locks; never follow symlinks or initialize descendants.
        let mut locks = Vec::new();
        let mut pending = vec![self.root().to_owned()];
        while let Some(directory) = pending.pop() {
            let entries = fs::read_dir(&directory).map_err(|source| {
                SnapshotError::io("list nested writer locks", &directory, source)
            })?;
            for entry in entries {
                let entry = entry.map_err(|source| {
                    SnapshotError::io("inspect nested writer entry", &directory, source)
                })?;
                let path = entry.path();
                let kind = entry.file_type().map_err(|source| {
                    SnapshotError::io("inspect nested writer type", &path, source)
                })?;
                if kind.is_dir() {
                    pending.push(path);
                } else if kind.is_file()
                    && entry.file_name() == SNAPSHOT_LOCK_FILE
                    && directory != self.root()
                {
                    let lock =
                        File::options()
                            .read(true)
                            .write(true)
                            .open(&path)
                            .map_err(|source| {
                                SnapshotError::io("open nested writer lock", &path, source)
                            })?;
                    Self::lock_file(&lock, &path)?;
                    Self::verify_lock(&lock, &path)?;
                    locks.push(lock);
                }
            }
        }
        Ok(locks)
    }

    fn remove_manifest_temporaries(
        &self,
        manifests: &[ManifestName<'_>],
    ) -> Result<(), SnapshotError> {
        let mut temporaries = Vec::new();
        for entry in fs::read_dir(self.root()).map_err(|source| {
            SnapshotError::io("list abandoned manifest files", self.root(), source)
        })? {
            let entry = entry.map_err(|source| {
                SnapshotError::io("inspect abandoned manifest file", self.root(), source)
            })?;
            let name = entry.file_name();
            let Some(name) = name.to_str() else {
                continue;
            };
            if manifests.iter().any(|manifest| manifest.is_temporary(name)) {
                let temporary = self.root.join(name)?;
                temporary.check_regular_or_missing()?;
                temporaries.push(temporary);
            }
        }
        for temporary in temporaries {
            temporary.remove_tree()?;
        }
        Ok(())
    }

    /// Borrow the resolved logical root; see the [crate example](crate).
    #[must_use]
    pub fn root(&self) -> &Path {
        self.root.as_path()
    }

    /// Borrow the checked root to join descendants without authorizing a new boundary.
    ///
    /// ```
    /// use alphal00p_persistence::SnapshotWriter;
    /// use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// # let temporary = tempfile::tempdir()?;
    /// let writer = SnapshotWriter::acquire(&temporary.path().join("state"), &[])?;
    /// let report = writer.checked_root().join_relative("reports/result.json".as_ref())?;
    /// assert!(report.as_path().starts_with(writer.root()));
    /// # Ok(()) }
    /// ```
    #[must_use]
    pub fn checked_root(&self) -> &OutputPath {
        &self.root
    }

    /// Exclusively create a fresh unpublished generation, retrying identity collisions.
    ///
    /// Its mutable borrow keeps the writer locked throughout staging and publication. See the
    /// [crate example](crate).
    ///
    /// # Errors
    /// Returns entry-validation, clock, or I/O errors without touching published generations.
    pub fn begin(&mut self) -> Result<StagedGeneration<'_>, SnapshotError> {
        self.root.check_directory()?;
        let generations = self.root.join(GENERATIONS_DIRECTORY)?;
        generations.ensure_directory()?;
        for attempt in 0..CREATE_ATTEMPTS {
            let id = GenerationId::generate()?;
            let path = generations.join(id.as_str())?.as_path().to_owned();
            match fs::create_dir(&path) {
                Ok(()) => {
                    return Ok(StagedGeneration {
                        writer: self,
                        generation: Some(PublishedGeneration { id, path }),
                        manifest_temp: None,
                        #[cfg(test)]
                        fault: None,
                        #[cfg(test)]
                        barrier: None,
                    });
                }
                Err(source)
                    if source.kind() == ErrorKind::AlreadyExists
                        && attempt + 1 < CREATE_ATTEMPTS => {}
                Err(source) => {
                    return Err(SnapshotError::io("create generation", &path, source));
                }
            }
        }
        unreachable!("the last creation attempt returns its error")
    }

    /// Remove unlisted generations while this writer holds the root lock.
    ///
    /// The caller supplies every current manifest generation and any coordinated reader pins.
    /// Ordinary publication retains everything; explicitly requested cleanup may invalidate older
    /// reader pins. This crate does not parse domain manifests and cannot infer their references.
    /// Non-generation names are ignored. All generation trees and keep IDs are checked before
    /// deleting anything; I/O failures during removal can still leave a partially pruned set.
    /// Returns the removed IDs in sorted order, for reporting by the caller.
    ///
    /// ```no_run
    /// use alphal00p_persistence::{GenerationId, ManifestName, SnapshotWriter};
    /// use std::error::Error;
    /// # fn prune(current: &GenerationId, previous: &GenerationId) -> Result<(), Box<dyn Error>> {
    /// let name = ManifestName::new("checkpoint.json")?;
    /// let mut writer = SnapshotWriter::acquire("checkpoints".as_ref(), &[name])?;
    /// writer.prune_generations([current, previous])?;
    /// # Ok(()) }
    /// ```
    ///
    /// # Errors
    /// Rejects missing keep generations and symlinks or special entries inside generation trees,
    /// and inspection/removal/synchronization failures. Published manifests are never modified.
    pub fn prune_generations<'a>(
        &mut self,
        keep: impl IntoIterator<Item = &'a GenerationId>,
    ) -> Result<Vec<GenerationId>, SnapshotError> {
        let keep: Vec<_> = keep.into_iter().collect();
        let generations = self.root.join(GENERATIONS_DIRECTORY)?;
        if generations.metadata()?.is_none() && keep.is_empty() {
            return Ok(Vec::new());
        }
        generations.check_directory()?;
        for id in &keep {
            generations.join(id.as_str())?.check_directory()?;
        }
        let mut remove = Vec::new();
        for entry in fs::read_dir(generations.as_path()).map_err(|source| {
            SnapshotError::io("list generations", generations.as_path(), source)
        })? {
            let entry = entry.map_err(|source| {
                SnapshotError::io("read generation entry", generations.as_path(), source)
            })?;
            let name = entry.file_name().to_string_lossy().into_owned();
            let Ok(id) = name.parse::<GenerationId>() else {
                continue;
            };
            let path = generations.join(id.as_str())?;
            path.check_directory()?;
            path.check_tree()?;
            if !keep.contains(&&id) {
                remove.push((id, path));
            }
        }
        remove.sort_by(|(left, _), (right, _)| left.cmp(right));
        let mut removed = Vec::with_capacity(remove.len());
        for (id, path) in remove {
            path.remove_tree()?;
            removed.push(id);
        }
        generations.sync_directory()?;
        Ok(removed)
    }

    /// Recognize a missing root or only unpublished runtime artifacts, without writes.
    ///
    /// Empty roots, caller-approved directories, regular lock files, valid generation directories,
    /// and temporary files for supplied manifests are scratch. A manifest or unrelated entry is
    /// not scratch. Approved directories are checked recursively.
    /// Temporary manifest bytes need not be complete: an interrupted first publication may leave
    /// a partially written file. This is classification, not destination authorization.
    ///
    /// ```no_run
    /// use alphal00p_persistence::{ManifestName, SnapshotWriter};
    /// use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// let name = ManifestName::new("session.toml")?;
    /// let scratch = SnapshotWriter::is_scratch("session".as_ref(), &[name], &["logs"])?;
    /// # let _ = scratch;
    /// # Ok(()) }
    /// ```
    ///
    /// # Errors
    /// Returns path-validation and inspection errors for symlinks, special entries, or malformed
    /// generation directories. Never creates a root or acquires a writer lock.
    pub fn is_scratch(
        root: &Path,
        manifests: &[ManifestName<'_>],
        scratch_directories: &[&str],
    ) -> Result<bool, SnapshotError> {
        for name in scratch_directories {
            OutputPath::validate_component(name)?;
        }
        let root = OutputPath::new(root)?;
        if root.metadata()?.is_none() {
            return Ok(true);
        }
        root.check_directory()?;
        for manifest in manifests {
            if root.join(manifest.as_str())?.metadata()?.is_some() {
                return Ok(false);
            }
        }
        for entry in fs::read_dir(root.as_path())
            .map_err(|source| SnapshotError::io("inspect scratch root", root.as_path(), source))?
        {
            let entry = entry.map_err(|source| {
                SnapshotError::io("inspect scratch entry", root.as_path(), source)
            })?;
            let name = entry.file_name();
            let Some(name) = name.to_str() else {
                return Ok(false);
            };
            let path = root.join(name)?;
            match name {
                SNAPSHOT_LOCK_FILE => path.check_regular_or_missing()?,
                GENERATIONS_DIRECTORY => {
                    path.check_directory()?;
                    for generation in fs::read_dir(path.as_path()).map_err(|source| {
                        SnapshotError::io("inspect scratch generations", path.as_path(), source)
                    })? {
                        let generation = generation.map_err(|source| {
                            SnapshotError::io("inspect scratch generation", path.as_path(), source)
                        })?;
                        let id: GenerationId = generation.file_name().to_string_lossy().parse()?;
                        let generation = path.join(id.as_str())?;
                        generation.check_directory()?;
                        generation.check_tree()?;
                    }
                }
                _ if manifests.iter().any(|manifest| manifest.is_temporary(name)) => {
                    path.check_regular_or_missing()?;
                }
                _ if scratch_directories.contains(&name) => {
                    path.check_directory()?;
                    path.check_tree()?;
                }
                _ => return Ok(false),
            }
        }
        Ok(true)
    }
}

/// Owns a detached snapshot tree and all its writer locks until explicit removal or drop.
///
/// Dropping this owner leaves the tree intact. See [`SnapshotWriter::retire_root`] for an example.
#[derive(Debug)]
pub struct RetiredSnapshotRoot {
    writer: SnapshotWriter,
    _nested_locks: Vec<File>,
}

impl RetiredSnapshotRoot {
    /// Remove only the detached tree, holding all locks until removal completes.
    ///
    /// The original root may already belong to a new writer; it is never inspected or modified.
    /// Recursive removal unlinks descendant symlinks and special entries without following them.
    /// Synchronizes the detached directory's parent on Unix. See [`SnapshotWriter::retire_root`].
    ///
    /// # Errors
    /// Returns stale-lock, validation, removal, or synchronization errors with path context.
    /// Failed removal can leave a partial tree for [`SnapshotWriter::remove_retired_roots`].
    /// Failed synchronization can leave no tree but uncertain durability; this never restores
    /// the original root.
    pub fn remove(self) -> Result<(), SnapshotError> {
        self.writer.verify_current_lock()?;
        self.writer.root.check_directory()?;
        let root = self.writer.root();
        self.remove_contents()?;
        // Keep the root lock until all payload deletion has finished. A concurrent cleaner
        // cannot recreate it and gain admission during the potentially long recursive removal.
        let lock = root.join(SNAPSHOT_LOCK_FILE);
        fs::remove_file(&lock)
            .map_err(|source| SnapshotError::io("remove retired snapshot lock", &lock, source))?;
        fs::remove_dir(root)
            .map_err(|source| SnapshotError::io("remove retired snapshot tree", root, source))?;
        let parent = root
            .parent()
            .expect("a detached root has a sibling's parent");
        OutputPath::new(parent)
            .and_then(|parent| parent.sync_directory())
            .map_err(|source| SnapshotError::Removed {
                path: root.to_owned(),
                source,
            })?;
        Ok(())
    }

    fn remove_contents(&self) -> Result<(), SnapshotError> {
        let root = self.writer.root();
        let entries = fs::read_dir(root)
            .map_err(|source| SnapshotError::io("list retired snapshot entries", root, source))?;
        for entry in entries {
            let entry = entry.map_err(|source| {
                SnapshotError::io("inspect retired snapshot entry", root, source)
            })?;
            if entry.file_name() == SNAPSHOT_LOCK_FILE {
                continue;
            }
            let path = entry.path();
            let kind = entry.file_type().map_err(|source| {
                SnapshotError::io("inspect retired snapshot entry type", &path, source)
            })?;
            let removal = if kind.is_dir() || kind.is_symlink() {
                fs::remove_dir_all(&path)
            } else {
                fs::remove_file(&path)
            };
            removal.map_err(|source| {
                SnapshotError::io("remove retired snapshot entry", &path, source)
            })?;
        }
        Ok(())
    }
}

/// Owns one unpublished generation and its manifest temporary file while borrowing the writer.
///
/// Dropping it removes only those unpublished entries, on a best effort basis. See the
/// [crate example](crate).
#[derive(Debug)]
pub struct StagedGeneration<'writer> {
    writer: &'writer mut SnapshotWriter,
    generation: Option<PublishedGeneration>,
    manifest_temp: Option<PathBuf>,
    #[cfg(test)]
    fault: Option<PublicationStep>,
    #[cfg(test)]
    barrier: Option<PublicationBarrier>,
}

impl StagedGeneration<'_> {
    /// Borrow the fresh payload directory; see the [crate example](crate).
    #[must_use]
    pub fn path(&self) -> &Path {
        self.generation().path()
    }

    /// Borrow the internally generated identity for the manifest; see the [crate example](crate).
    #[must_use]
    pub fn id(&self) -> &GenerationId {
        self.generation().id()
    }

    /// Synchronize all regular payload files, then write and synchronize a unique manifest
    /// temporary file, close it, and atomically replace the pointer.
    ///
    /// The rename is the commit point. Immediately afterward, cleanup is disarmed and this method
    /// synchronizes the root directory on Unix before returning. It never unlinks an old
    /// manifest first and never deletes prior committed generations. See the
    /// [crate example](crate).
    ///
    /// # Errors
    /// All errors except [`SnapshotError::Committed`] precede the commit point, leaving the old
    /// manifest untouched. A committed error preserves the new generation and means its visible
    /// pointer could not be confirmed durable; do not delete prior state. Symlinks, special
    /// payload entries, failed writes/flushes/synchronization, unsupported operations, and failed
    /// replacement are reported with context. Callers must finish and check their writes and
    /// application-buffer flushes before invoking this method, and stop modifying the generation.
    /// Payload files must remain writable on Windows. On Unix, files and directories synchronize
    /// bottom up before commit, then the manifest directory is synchronized after replacement.
    /// Other platforms cannot promise directory durability. See the [crate contract](crate)
    /// for synchronization cost.
    pub fn publish(
        mut self,
        manifest: &[u8],
        name: ManifestName<'_>,
    ) -> Result<PublishedGeneration, SnapshotError> {
        self.writer.root.check_directory()?;
        self.writer
            .root
            .join(GENERATIONS_DIRECTORY)?
            .check_directory()?;
        self.output_path()?.check_directory()?;
        self.complete_payload()?;
        self.writer
            .root
            .join(GENERATIONS_DIRECTORY)?
            .sync_directory()?;
        let manifest_output = self.writer.root.join(name.as_str())?;
        let manifest_path = manifest_output.as_path();
        manifest_output.check_regular_or_missing()?;
        let (temp, mut file) = self.create_manifest_temp(name)?;
        self.manifest_temp = Some(temp);
        let temp = self
            .manifest_temp
            .as_deref()
            .expect("the exclusively created manifest temporary file is owned until commit");
        #[cfg(test)]
        if self.fault == Some(PublicationStep::ManifestWrite) {
            file.write_all(&manifest[..manifest.len() / 2])
                .map_err(|source| SnapshotError::io("write partial test manifest", temp, source))?;
            self.injected_error(PublicationStep::ManifestWrite)
                .map_err(|source| {
                    SnapshotError::io("write manifest temporary file", temp, source)
                })?;
        }
        file.write_all(manifest)
            .map_err(|source| SnapshotError::io("write manifest temporary file", temp, source))?;
        file.flush()
            .map_err(|source| SnapshotError::io("flush manifest temporary file", temp, source))?;
        let completion = file.sync_all();
        #[cfg(test)]
        let completion =
            completion.and_then(|()| self.injected_error(PublicationStep::ManifestCompletion));
        completion.map_err(|source| {
            SnapshotError::io("complete manifest temporary file", temp, source)
        })?;
        drop(file);
        self.writer.root.sync_directory()?;
        manifest_output.check_regular_or_missing()?;
        #[cfg(test)]
        self.injected_error(PublicationStep::ManifestRename)
            .map_err(|source| SnapshotError::io("publish manifest", manifest_path, source))?;
        #[cfg(test)]
        if let Some(barrier) = &self.barrier {
            barrier.wait(PublicationStep::ManifestRename);
        }
        fs::rename(temp, manifest_path)
            .map_err(|source| SnapshotError::io("publish manifest", manifest_path, source))?;
        let published = self
            .generation
            .take()
            .expect("a staged generation is present until its successful publication");
        self.manifest_temp = None;
        #[cfg(test)]
        if let Some(barrier) = &self.barrier {
            barrier.wait(PublicationStep::Committed);
        }
        let completion = self.writer.root.sync_directory();
        #[cfg(test)]
        let completion = completion.and_then(|()| {
            self.injected_error(PublicationStep::DirectoryCompletion)
                .map_err(|source| OutputPathError::io(self.writer.root(), source))
        });
        match completion {
            Ok(()) => Ok(published),
            Err(source) => Err(SnapshotError::Committed {
                generation: published,
                source,
            }),
        }
    }

    /// Retain the writer's managed boundary when accessing this generation's descendants.
    /// See [`OutputPath::join_relative`] for checked descendant construction.
    ///
    /// # Errors
    /// Rejects symlinks or inspection failures in the managed path.
    pub fn output_path(&self) -> Result<OutputPath, SnapshotError> {
        Ok(self
            .writer
            .root
            .join(GENERATIONS_DIRECTORY)?
            .join(self.id().as_str())?)
    }

    fn complete_payload(&self) -> Result<(), SnapshotError> {
        self.output_path()?.sync_tree()?;
        #[cfg(test)]
        self.injected_error(PublicationStep::PayloadCompletion)
            .map_err(|source| SnapshotError::io("complete payload file", self.path(), source))?;
        Ok(())
    }

    fn generation(&self) -> &PublishedGeneration {
        self.generation
            .as_ref()
            .expect("a staged generation is present until its successful publication")
    }

    fn create_manifest_temp(
        &self,
        name: ManifestName<'_>,
    ) -> Result<(PathBuf, File), SnapshotError> {
        for attempt in 0..CREATE_ATTEMPTS {
            let nonce = GenerationId::generate()?;
            let path = self
                .writer
                .root()
                .join(format!(".{}.{nonce}.tmp", name.as_str()));
            match File::create_new(&path) {
                Ok(file) => return Ok((path, file)),
                Err(source)
                    if source.kind() == ErrorKind::AlreadyExists
                        && attempt + 1 < CREATE_ATTEMPTS => {}
                Err(source) => {
                    return Err(SnapshotError::io(
                        "create manifest temporary file",
                        &path,
                        source,
                    ));
                }
            }
        }
        unreachable!("the last creation attempt returns its error")
    }

    #[cfg(test)]
    fn injected_error(&self, step: PublicationStep) -> Result<(), IoError> {
        if self.fault == Some(step) {
            return Err(IoError::other("test publication fault"));
        }
        Ok(())
    }
}

impl Drop for StagedGeneration<'_> {
    fn drop(&mut self) {
        if self.writer.root.check_directory().is_err() {
            return;
        }
        if let Some(temp) = &self.manifest_temp {
            let _ = fs::remove_file(temp);
        }
        if self.generation.is_some()
            && let Ok(path) = self.output_path()
            && path.check_directory().is_ok()
        {
            // This owned directory may contain invalid symlinks; remove_dir_all unlinks them.
            let _ = fs::remove_dir_all(path.as_path());
        }
    }
}

/// A committed generation, retained independently of the writer and future publications.
///
/// Dropping this descriptor does not delete its files. See the [crate example](crate).
#[derive(Clone, Debug, Eq, Hash, PartialEq)]
pub struct PublishedGeneration {
    id: GenerationId,
    path: PathBuf,
}

impl PublishedGeneration {
    /// Borrow the committed identity; see the [crate example](crate).
    #[must_use]
    pub fn id(&self) -> &GenerationId {
        &self.id
    }

    /// Borrow the retained payload directory; see the [crate example](crate).
    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }
}

#[cfg(test)]
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum PublicationStep {
    PayloadCompletion,
    ManifestWrite,
    ManifestCompletion,
    ManifestRename,
    Committed,
    DirectoryCompletion,
}

#[cfg(test)]
#[derive(Debug)]
struct PublicationBarrier {
    step: PublicationStep,
    ready: PathBuf,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::DirectoryPublication;
    use crate::test_support::TestDirectory;
    use std::env;
    use std::error::Error as StdError;
    #[cfg(unix)]
    use std::ffi::OsStr;
    #[cfg(unix)]
    use std::os::unix::ffi::OsStrExt;
    #[cfg(unix)]
    use std::os::unix::fs as unix_fs;
    #[cfg(unix)]
    use std::os::unix::fs::{MetadataExt, PermissionsExt};
    #[cfg(unix)]
    use std::os::unix::net::UnixListener;
    #[cfg(unix)]
    use std::os::unix::process::ExitStatusExt;
    use std::process::{Child, Command, ExitStatus};
    use std::thread;
    use std::time::{Duration, Instant};

    const STATE_MANIFEST_FILE: &str = "state_manifest.toml";
    const INTEGRATION_MANIFEST_FILE: &str = "manifest.json";
    const STATE: ManifestName<'static> = ManifestName(STATE_MANIFEST_FILE);
    const INTEGRATION: ManifestName<'static> = ManifestName(INTEGRATION_MANIFEST_FILE);
    const MANIFESTS: &[ManifestName<'static>] = &[STATE, INTEGRATION];

    const CHILD_TIMEOUT: Duration = Duration::from_secs(15);
    const POLL_INTERVAL: Duration = Duration::from_millis(10);

    impl PublicationBarrier {
        pub(super) fn wait(&self, step: PublicationStep) {
            if self.step != step {
                return;
            }
            fs::write(&self.ready, b"ready").unwrap();
            // A crashed parent must not leave a permanent worker behind.
            let deadline = Instant::now() + Duration::from_secs(30);
            while Instant::now() < deadline {
                thread::sleep(POLL_INTERVAL);
            }
            panic!("publication barrier was not interrupted within 30 seconds");
        }
    }

    struct ChildProcess {
        child: Child,
        control: PathBuf,
    }

    impl ChildProcess {
        fn spawn(fixture: &Fixture, action: &str, expected: &str, directory: &Path) -> Self {
            let control = fixture.0.0.join(action);
            fs::create_dir(&control).unwrap();
            let log = File::create_new(control.join("child.log")).unwrap();
            // Both streams deliberately share one descriptor for diagnostics on assertion failure.
            let stderr = log.try_clone().unwrap();
            let child = Command::new(env::current_exe().unwrap())
                .current_dir(directory)
                .args([
                    "--exact",
                    "snapshot::tests::subprocess_entry",
                    "--ignored",
                    "--nocapture",
                ])
                .env("GAMMALOOP_SNAPSHOT_TEST_ACTION", action)
                .env("GAMMALOOP_SNAPSHOT_TEST_ROOT", fixture.root())
                .env("GAMMALOOP_SNAPSHOT_TEST_CONTROL", &control)
                .env("GAMMALOOP_SNAPSHOT_TEST_EXPECTED", expected)
                .stdout(log)
                .stderr(stderr)
                .spawn()
                .unwrap();
            Self { child, control }
        }

        fn log(&self) -> String {
            fs::read_to_string(self.control.join("child.log")).unwrap_or_default()
        }

        fn wait_ready(&mut self) {
            let deadline = Instant::now() + CHILD_TIMEOUT;
            while !self.control.join("ready").is_file() {
                assert!(
                    self.child.try_wait().unwrap().is_none(),
                    "writer exited before its barrier: {}",
                    self.log()
                );
                assert!(
                    Instant::now() < deadline,
                    "writer barrier timeout: {}",
                    self.log()
                );
                thread::sleep(POLL_INTERVAL);
            }
        }

        fn wait(&mut self) -> ExitStatus {
            let deadline = Instant::now() + CHILD_TIMEOUT;
            loop {
                if let Some(status) = self.child.try_wait().unwrap() {
                    return status;
                }
                assert!(Instant::now() < deadline, "child timeout: {}", self.log());
                thread::sleep(POLL_INTERVAL);
            }
        }

        fn assert_success(&mut self) {
            let status = self.wait();
            assert!(status.success(), "child failed ({status}): {}", self.log());
            // An unmatched libtest filter also exits successfully; require the protocol receipt.
            assert_eq!(fs::read(self.control.join("done")).unwrap(), b"ok");
        }

        fn kill(&mut self) {
            self.child.kill().unwrap();
            let status = self.wait();
            assert!(
                !status.success(),
                "the interrupted writer exited successfully"
            );
            #[cfg(unix)]
            assert_eq!(
                status.signal(),
                Some(9),
                "the writer must be killed without Drop"
            );
        }
    }

    impl Drop for ChildProcess {
        fn drop(&mut self) {
            // Always reap the child, including assertion failures and timeouts in the parent.
            let _ = self.child.kill();
            let _ = self.child.wait();
        }
    }

    struct Fixture(TestDirectory);

    #[cfg(unix)]
    struct SearchOnlyDirectory<'path> {
        path: &'path Path,
        original_mode: u32,
    }

    #[cfg(unix)]
    impl<'path> SearchOnlyDirectory<'path> {
        fn new(path: &'path Path) -> Option<Self> {
            let restriction = Self {
                path,
                original_mode: fs::metadata(path).unwrap().permissions().mode(),
            };
            fs::set_permissions(path, PermissionsExt::from_mode(0o111)).unwrap();
            match File::open(path) {
                Ok(_) => {
                    eprintln!("skipping permission test: mode 0111 remains readable for this user");
                    return None;
                }
                Err(error) => assert_eq!(error.kind(), ErrorKind::PermissionDenied),
            }
            Some(restriction)
        }
    }

    #[cfg(unix)]
    impl Drop for SearchOnlyDirectory<'_> {
        fn drop(&mut self) {
            fs::set_permissions(self.path, PermissionsExt::from_mode(self.original_mode)).unwrap();
        }
    }

    struct DelayedWriter {
        root: OutputPath,
        lock: File,
    }

    impl DelayedWriter {
        fn open(root: &Path) -> Self {
            let root = OutputPath::new(root).unwrap();
            let lock = File::options()
                .read(true)
                .write(true)
                .open(root.as_path().join(SNAPSHOT_LOCK_FILE))
                .unwrap();
            Self { root, lock }
        }

        fn assert_stale(self) {
            let expected = self.root.as_path().join(SNAPSHOT_LOCK_FILE);
            let error =
                SnapshotWriter::finish_acquire(self.root, self.lock, MANIFESTS).unwrap_err();
            assert!(matches!(error, SnapshotError::StaleLock { path } if path == expected));
        }
    }

    impl Fixture {
        fn new() -> Self {
            Self(TestDirectory::new())
        }

        fn root(&self) -> PathBuf {
            self.0.0.join("snapshot")
        }

        fn publish(
            writer: &mut SnapshotWriter,
            marker: &[u8],
            name: ManifestName<'_>,
        ) -> PublishedGeneration {
            let stage = writer.begin().unwrap();
            fs::write(stage.path().join("payload"), marker).unwrap();
            stage.publish(marker, name).unwrap()
        }

        fn entries(path: &Path) -> Vec<PathBuf> {
            let mut entries: Vec<_> = fs::read_dir(path)
                .unwrap()
                .map(|entry| entry.unwrap().path())
                .collect();
            entries.sort();
            entries
        }

        fn publish_pointer(
            writer: &mut SnapshotWriter,
            marker: &[u8],
            barrier: Option<PublicationBarrier>,
        ) -> PublishedGeneration {
            let mut stage = writer.begin().unwrap();
            fs::write(stage.path().join("payload"), marker).unwrap();
            fs::create_dir(stage.path().join("nested")).unwrap();
            fs::write(stage.path().join("nested/member"), marker).unwrap();
            let manifest = stage.id().to_string();
            stage.barrier = barrier;
            stage.publish(manifest.as_bytes(), INTEGRATION).unwrap()
        }

        fn read_pointer(root: &Path) -> GenerationId {
            let manifest = fs::read_to_string(root.join("manifest.json")).unwrap();
            manifest.parse().unwrap()
        }

        fn assert_payload(path: &Path, expected: &[u8]) {
            assert_eq!(fs::read(path.join("payload")).unwrap(), expected);
            assert_eq!(fs::read(path.join("nested/member")).unwrap(), expected);
        }

        fn interrupted_publication(step: PublicationStep) {
            let action = match step {
                PublicationStep::ManifestRename => "before",
                PublicationStep::Committed => "after",
                _ => panic!("only pre-rename and post-commit barriers are valid"),
            };
            let committed = step == PublicationStep::Committed;
            let fixture = Self::new();
            let root = fixture.root();
            let mut writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
            let first = Self::publish_pointer(&mut writer, b"A", None);
            drop(writer);
            let lock_path = root.join(".snapshot.lock");
            #[cfg(unix)]
            let original_lock = fs::metadata(&lock_path).unwrap();

            let mut interrupted = ChildProcess::spawn(&fixture, action, "", &fixture.0.0);
            interrupted.wait_ready();
            let generations = Self::entries(&root.join("generations"));
            assert_eq!(generations.len(), 2);
            let second = generations
                .into_iter()
                .find(|path| path != first.path())
                .unwrap();
            let second_id: GenerationId = second
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .parse()
                .unwrap();
            let expected = if committed { "B" } else { "A" };
            let expected_id = if committed { &second_id } else { first.id() };
            assert_eq!(&Self::read_pointer(&root), expected_id);
            Self::assert_payload(first.path(), b"A");
            Self::assert_payload(&second, b"B");
            let temporaries: Vec<_> = Self::entries(&root)
                .into_iter()
                .filter(|path| {
                    path.file_name()
                        .unwrap()
                        .to_str()
                        .unwrap()
                        .starts_with(".manifest.json.")
                })
                .collect();
            assert_eq!(temporaries.len(), usize::from(!committed));

            let mut contender = ChildProcess::spawn(&fixture, "busy", "", &fixture.0.0);
            contender.assert_success();
            interrupted.kill();
            // SIGKILL skips all destructors: unpublished payloads and the pre-rename temp survive.
            assert_eq!(&Self::read_pointer(&root), expected_id);
            Self::assert_payload(first.path(), b"A");
            Self::assert_payload(&second, b"B");
            assert!(temporaries.iter().all(|path| path.is_file()));
            assert!(fs::symlink_metadata(&lock_path).unwrap().is_file());

            let mut recovery = ChildProcess::spawn(&fixture, "recover", expected, &fixture.0.0);
            recovery.assert_success();
            let third_id = Self::read_pointer(&root);
            assert_ne!(&third_id, first.id());
            assert_ne!(third_id, second_id);
            Self::assert_payload(&root.join("generations").join(third_id.as_str()), b"C");
            Self::assert_payload(first.path(), b"A");
            Self::assert_payload(&second, b"B");
            assert!(temporaries.iter().all(|path| !path.exists()));
            assert_eq!(Self::entries(&root.join("generations")).len(), 3);
            #[cfg(unix)]
            {
                let retained_lock = fs::metadata(lock_path).unwrap();
                assert_eq!(
                    (retained_lock.dev(), retained_lock.ino()),
                    (original_lock.dev(), original_lock.ino()),
                    "publication and recovery must retain the same lock inode"
                );
            }
        }
    }

    #[test]
    #[ignore = "subprocess protocol entry; invoked by the bounded interruption tests"]
    fn subprocess_entry() {
        // Dispatch exists only in the test binary; production has no environment fault switch.
        let Ok(action) = env::var("GAMMALOOP_SNAPSHOT_TEST_ACTION") else {
            return;
        };
        let root = PathBuf::from(env::var_os("GAMMALOOP_SNAPSHOT_TEST_ROOT").unwrap());
        let control = PathBuf::from(env::var_os("GAMMALOOP_SNAPSHOT_TEST_CONTROL").unwrap());
        match action.as_str() {
            "before" | "after" => {
                let mut writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
                let step = if action == "before" {
                    PublicationStep::ManifestRename
                } else {
                    PublicationStep::Committed
                };
                Fixture::publish_pointer(
                    &mut writer,
                    b"B",
                    Some(PublicationBarrier {
                        step,
                        ready: control.join("ready"),
                    }),
                );
                panic!("the writer must be interrupted at its publication barrier");
            }
            "busy" => {
                assert!(matches!(
                    SnapshotWriter::acquire(&root, MANIFESTS),
                    Err(SnapshotError::Busy { .. })
                ));
            }
            "recover" => {
                let mut writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
                let expected = env::var("GAMMALOOP_SNAPSHOT_TEST_EXPECTED").unwrap();
                let pinned = Fixture::read_pointer(&root);
                Fixture::assert_payload(
                    &root.join("generations").join(pinned.as_str()),
                    expected.as_bytes(),
                );
                Fixture::publish_pointer(&mut writer, b"C", None);
            }
            "retire" => {
                let writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
                let _retired = writer.retire_root().unwrap();
                PublicationBarrier {
                    step: PublicationStep::Committed,
                    ready: control.join("ready"),
                }
                .wait(PublicationStep::Committed);
                panic!("the cleaner must be interrupted before removal");
            }
            #[cfg(unix)]
            "socket" => {
                drop(UnixListener::bind("invalid").unwrap());
            }
            _ => panic!("unknown snapshot subprocess action {action:?}"),
        }
        fs::write(control.join("done"), b"ok").unwrap();
    }

    #[test]
    fn killing_writer_before_commit_preserves_a_and_releases_os_lock() {
        Fixture::interrupted_publication(PublicationStep::ManifestRename);
    }

    #[test]
    fn killing_writer_after_commit_retains_a_and_b_and_releases_os_lock() {
        Fixture::interrupted_publication(PublicationStep::Committed);
    }

    #[test]
    fn killing_cleaner_after_detachment_allows_next_cleanup_to_reclaim_its_tree() {
        let fixture = Fixture::new();
        let root = fixture.root();
        let mut writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        Fixture::publish(&mut writer, b"old state", STATE);
        drop(writer);
        drop(SnapshotWriter::acquire(&root.join("nested"), MANIFESTS).unwrap());
        let mut cleaner = ChildProcess::spawn(&fixture, "retire", "", &fixture.0.0);
        cleaner.wait_ready();
        assert!(!root.exists());
        let prefix = SnapshotWriter::retired_prefix(&root).unwrap();
        let siblings: Vec<_> = Fixture::entries(&fixture.0.0)
            .into_iter()
            .filter(|path| {
                path.file_name()
                    .unwrap()
                    .to_str()
                    .unwrap()
                    .starts_with(&prefix)
            })
            .collect();
        assert_eq!(siblings.len(), 1);
        SnapshotWriter::remove_retired_roots(&root).unwrap();
        assert_eq!(
            fs::read(siblings[0].join(STATE_MANIFEST_FILE)).unwrap(),
            b"old state"
        );
        cleaner.kill();
        SnapshotWriter::remove_retired_roots(&root).unwrap();
        assert!(!siblings[0].exists());
        assert!(!root.exists());
    }

    #[test]
    fn first_and_repeated_publication_retain_pinned_generations() {
        for name in [STATE, INTEGRATION] {
            let fixture = Fixture::new();
            let mut writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
            let manifest = writer.root().join(name.as_str());
            assert!(!manifest.exists());
            let first = Fixture::publish(&mut writer, b"first", name);
            assert_eq!(fs::read(&manifest).unwrap(), b"first");
            let second = Fixture::publish(&mut writer, b"second", name);
            assert_eq!(fs::read(&manifest).unwrap(), b"second");
            assert_ne!(first.id(), second.id());
            drop(writer);
            assert_eq!(fs::read(first.path().join("payload")).unwrap(), b"first");
            assert_eq!(fs::read(second.path().join("payload")).unwrap(), b"second");
            assert_eq!(
                Fixture::entries(&fixture.root().join("generations")).len(),
                2
            );
            assert_eq!(Fixture::entries(&fixture.root()).len(), 3);
            assert!(fixture.root().join(".snapshot.lock").is_file());
        }
    }

    #[test]
    fn precommit_failures_preserve_old_manifest_and_clean_owned_entries() {
        for fault in [
            PublicationStep::PayloadCompletion,
            PublicationStep::ManifestWrite,
            PublicationStep::ManifestCompletion,
            PublicationStep::ManifestRename,
        ] {
            let fixture = Fixture::new();
            let mut writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
            let first = Fixture::publish(&mut writer, b"old", STATE);
            let root_entries = Fixture::entries(writer.root());
            let mut stage = writer.begin().unwrap();
            let abandoned = stage.path().to_owned();
            fs::write(stage.path().join("payload"), b"new").unwrap();
            stage.fault = Some(fault);
            let error = stage.publish(b"new", STATE).unwrap_err();
            assert!(matches!(error, SnapshotError::Io { .. }));
            assert!(error.source().is_some());
            if let SnapshotError::Io { operation, .. } = &error {
                match fault {
                    PublicationStep::PayloadCompletion => {
                        assert_eq!(*operation, "complete payload file")
                    }
                    PublicationStep::ManifestCompletion => {
                        assert_eq!(*operation, "complete manifest temporary file")
                    }
                    _ => {}
                }
            }
            assert_eq!(
                fs::read(writer.root().join("state_manifest.toml")).unwrap(),
                b"old"
            );
            assert_eq!(fs::read(first.path().join("payload")).unwrap(), b"old");
            assert!(!abandoned.exists());
            assert_eq!(Fixture::entries(writer.root()), root_entries);
            assert_eq!(
                Fixture::entries(&writer.root().join("generations")),
                [first.path()]
            );
        }
    }

    #[test]
    fn failed_first_publication_leaves_no_manifest_or_generation() {
        let fixture = Fixture::new();
        let mut writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        let mut stage = writer.begin().unwrap();
        stage.fault = Some(PublicationStep::ManifestRename);
        assert!(stage.publish(b"first", INTEGRATION).is_err());
        assert!(!writer.root().join("manifest.json").exists());
        assert!(Fixture::entries(&writer.root().join("generations")).is_empty());
        assert_eq!(Fixture::entries(writer.root()).len(), 2);
    }

    #[test]
    fn dropping_stage_removes_only_its_unpublished_payload() {
        let fixture = Fixture::new();
        let mut writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        let committed = Fixture::publish(&mut writer, b"old", INTEGRATION);
        let foreign = writer.root().join(".unrelated.tmp");
        fs::write(&foreign, b"another owner").unwrap();
        let orphan = writer.root().join("generations").join("unreferenced");
        fs::create_dir(&orphan).unwrap();
        let stage = writer.begin().unwrap();
        let abandoned = stage.path().to_owned();
        fs::create_dir(stage.path().join("nested")).unwrap();
        fs::write(stage.path().join("nested/payload"), b"uncommitted").unwrap();
        drop(stage);
        assert!(!abandoned.exists());
        assert!(committed.path().is_dir());
        assert!(orphan.is_dir());
        assert_eq!(fs::read(foreign).unwrap(), b"another owner");
        assert_eq!(
            fs::read(writer.root().join("manifest.json")).unwrap(),
            b"old"
        );
        Fixture::publish(&mut writer, b"next", INTEGRATION);
    }

    #[test]
    fn generation_parsing_and_serde_reject_escaping_and_malformed_components() {
        let generated = GenerationId::generate().unwrap();
        let parsed: GenerationId = generated.as_str().parse().unwrap();
        assert_eq!(parsed, generated);
        assert_eq!(
            GenerationId::try_from(generated.to_string()).unwrap(),
            generated
        );
        let encoded = serde_json::to_string(&generated).unwrap();
        assert_eq!(
            serde_json::from_str::<GenerationId>(&encoded).unwrap(),
            generated
        );
        for invalid in [
            "",
            ".",
            "..",
            "../generation",
            "/generation",
            "a/b",
            "a\\b",
            "C:escape",
            "generation",
            "gen-1-2-3",
            "gen-00000001-00000000000000000000000000000000-000000000000000A",
            "gen-00000001-00000000000000000000000000000000-0000000000000000-extra",
            "gen-00000001-00000000000000000000000000000000-000000000000000/",
        ] {
            assert!(
                matches!(
                    invalid.parse::<GenerationId>(),
                    Err(SnapshotError::InvalidGeneration { .. })
                ),
                "{invalid:?}"
            );
            let encoded = serde_json::to_string(invalid).unwrap();
            assert!(
                serde_json::from_str::<GenerationId>(&encoded).is_err(),
                "{invalid:?}"
            );
        }
    }

    #[test]
    fn caller_manifest_names_validate_and_display_without_domain_variants() {
        let filename = String::from("application-checkpoint.json");
        let manifest = ManifestName::new(&filename).unwrap();
        assert_eq!(manifest.as_str(), filename);
        assert_eq!(manifest.to_string(), filename);
        for invalid in ["", ".", "..", "../escape", "/absolute", "a/b", "a\\b"] {
            assert!(matches!(
                ManifestName::new(invalid),
                Err(SnapshotError::Path(OutputPathError::InvalidComponent(_)))
            ));
        }
        for reserved in [SNAPSHOT_LOCK_FILE, GENERATIONS_DIRECTORY] {
            assert!(matches!(
                ManifestName::new(reserved),
                Err(SnapshotError::ReservedManifest { .. })
            ));
        }
    }

    #[test]
    fn acquisition_removes_only_known_manifest_temporaries_under_the_lock() {
        let fixture = Fixture::new();
        let manifest = ManifestName::new("application.json").unwrap();
        let writer = SnapshotWriter::acquire(&fixture.root(), &[manifest]).unwrap();
        let id = GenerationId::generate().unwrap();
        let abandoned = writer.root().join(format!(".{manifest}.{id}.tmp"));
        let other_manifest = writer.root().join(format!(".other.json.{id}.tmp"));
        let malformed = writer
            .root()
            .join(format!(".{manifest}.not-a-generation.tmp"));
        for path in [&abandoned, &other_manifest, &malformed] {
            fs::write(path, b"partial").unwrap();
        }
        assert!(matches!(
            SnapshotWriter::acquire(&fixture.root(), &[manifest]),
            Err(SnapshotError::Busy { .. })
        ));
        assert!(
            abandoned.is_file(),
            "a busy acquisition cannot remove a writer's file"
        );
        drop(writer);
        let _writer = SnapshotWriter::acquire(&fixture.root(), &[manifest]).unwrap();
        assert!(!abandoned.exists());
        assert_eq!(fs::read(other_manifest).unwrap(), b"partial");
        assert_eq!(fs::read(malformed).unwrap(), b"partial");
    }

    #[test]
    fn retirement_refuses_busy_nested_writers_without_changing_the_root() {
        let fixture = Fixture::new();
        let root = fixture.root();
        let writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        let nested_root = root.join("new/nested");
        let mut nested = SnapshotWriter::acquire(&nested_root, MANIFESTS).unwrap();
        let before = Fixture::entries(&root);
        let error = writer.retire_root().unwrap_err();
        let expected_lock = nested_root.join(SNAPSHOT_LOCK_FILE);
        assert!(error.is_busy());
        assert_eq!(
            error.to_string(),
            format!("snapshot writer is busy at {}", expected_lock.display())
        );
        assert!(matches!(error, SnapshotError::Busy { path } if path == expected_lock));
        assert_eq!(Fixture::entries(&root), before);
        assert_eq!(Fixture::entries(&fixture.0.0), [root]);
        Fixture::publish(&mut nested, b"still owned", STATE);
    }

    #[test]
    fn cleaner_nested_lock_acquisition_preserves_temporaries_until_removal() {
        let fixture = Fixture::new();
        let root = fixture.root();
        let writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        let nested_root = root.join("nested");
        drop(SnapshotWriter::acquire(&nested_root, MANIFESTS).unwrap());
        let temporary = format!(".{STATE}.{}.tmp", GenerationId::generate().unwrap());
        fs::write(nested_root.join(&temporary), b"unpublished").unwrap();
        let retired = writer.retire_root().unwrap();
        let detached_nested = retired.writer.root().join("nested");
        assert_eq!(
            fs::read(detached_nested.join(temporary)).unwrap(),
            b"unpublished"
        );
        // There is no container, extra lock, or marker around the detached root.
        assert_eq!(Fixture::entries(retired.writer.root()).len(), 2);
        assert!(retired.writer.root().join(SNAPSHOT_LOCK_FILE).is_file());
        retired.remove().unwrap();
        assert!(Fixture::entries(&fixture.0.0).is_empty());
    }

    #[test]
    fn retirement_removal_refuses_replaced_detached_lock() {
        let fixture = Fixture::new();
        let writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        let retired = writer.retire_root().unwrap();
        let retained = retired.writer.root().to_owned();
        let lock = retained.join(SNAPSHOT_LOCK_FILE);
        fs::remove_file(&lock).unwrap();
        fs::write(&lock, b"replacement lock").unwrap();
        let error = retired.remove().unwrap_err();
        assert!(matches!(error, SnapshotError::StaleLock { path } if path == lock));
        assert_eq!(fs::read(lock).unwrap(), b"replacement lock");
        assert!(retained.is_dir());
    }

    #[test]
    fn delayed_root_and_nested_openers_refuse_a_removed_root() {
        let fixture = Fixture::new();
        let root = fixture.root();
        let nested_root = root.join("nested");
        let writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        let nested = SnapshotWriter::acquire(&nested_root, MANIFESTS).unwrap();
        let delayed = [
            DelayedWriter::open(&root),
            DelayedWriter::open(&nested_root),
        ];

        drop(nested);
        let retired = writer.retire_root().unwrap();
        assert!(!root.exists());
        for opener in &delayed {
            assert!(matches!(
                opener.lock.try_lock(),
                Err(TryLockError::WouldBlock)
            ));
        }
        retired.remove().unwrap();
        for opener in delayed {
            opener.assert_stale();
        }
        assert!(!root.exists(), "a stale opener must not recreate the root");
        assert!(Fixture::entries(&fixture.0.0).is_empty());
    }

    #[test]
    fn retired_cleanup_and_stale_openers_leave_replacement_writers_untouched() {
        let fixture = Fixture::new();
        let root = fixture.root();
        let nested_root = root.join("nested");
        let writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        let nested = SnapshotWriter::acquire(&nested_root, MANIFESTS).unwrap();
        let delayed = [
            DelayedWriter::open(&root),
            DelayedWriter::open(&nested_root),
        ];
        drop(nested);
        let retired = writer.retire_root().unwrap();
        let retired_path = retired.writer.root().to_owned();

        let mut replacement = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        let mut replacement_nested = SnapshotWriter::acquire(&nested_root, MANIFESTS).unwrap();
        let published = Fixture::publish(&mut replacement, b"new state", STATE);
        let published_nested = Fixture::publish(&mut replacement_nested, b"new nested", STATE);
        let temporary = format!(".{STATE}.{}.tmp", GenerationId::generate().unwrap());
        for root in [&root, &nested_root] {
            fs::write(root.join(&temporary), b"live writer temporary").unwrap();
        }

        retired.remove().unwrap();
        assert!(!retired_path.exists());
        for opener in delayed {
            opener.assert_stale();
        }
        for root in [&root, &nested_root] {
            assert_eq!(
                fs::read(root.join(&temporary)).unwrap(),
                b"live writer temporary"
            );
            assert!(matches!(
                SnapshotWriter::acquire(root, MANIFESTS),
                Err(SnapshotError::Busy { .. })
            ));
        }
        assert_eq!(
            fs::read(root.join(STATE_MANIFEST_FILE)).unwrap(),
            b"new state"
        );
        assert_eq!(
            fs::read(nested_root.join(STATE_MANIFEST_FILE)).unwrap(),
            b"new nested"
        );
        assert_eq!(
            fs::read(published.path().join("payload")).unwrap(),
            b"new state"
        );
        assert_eq!(
            fs::read(published_nested.path().join("payload")).unwrap(),
            b"new nested"
        );
        Fixture::publish(&mut replacement, b"next state", STATE);
        Fixture::publish(&mut replacement_nested, b"next nested", STATE);
    }

    #[test]
    fn dropping_retired_owner_preserves_tree_and_next_cleanup_removes_it() {
        let fixture = Fixture::new();
        let root = fixture.root();
        let mut writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        Fixture::publish(&mut writer, b"retained state", STATE);
        drop(SnapshotWriter::acquire(&root.join("nested"), MANIFESTS).unwrap());
        let retired = writer.retire_root().unwrap();
        let retired_root = retired.writer.root().to_owned();
        let retired_nested = retired_root.join("nested");
        for path in [&retired_root, &retired_nested] {
            assert!(matches!(
                SnapshotWriter::acquire(path, MANIFESTS),
                Err(SnapshotError::Busy { .. })
            ));
        }
        SnapshotWriter::remove_retired_roots(&root).unwrap();
        assert!(
            retired_root.exists(),
            "a live cleaner retains its detached tree"
        );
        drop(retired);
        assert!(!root.exists());
        assert_eq!(
            fs::read(retired_root.join(STATE_MANIFEST_FILE)).unwrap(),
            b"retained state"
        );
        // A busy descendant also prevents sweeping, even after the root owner is gone.
        let nested = SnapshotWriter::acquire(&retired_nested, MANIFESTS).unwrap();
        SnapshotWriter::remove_retired_roots(&root).unwrap();
        assert!(retired_root.exists());
        drop(nested);
        SnapshotWriter::remove_retired_roots(&root).unwrap();
        assert!(!root.exists(), "sweeping must not recreate the live root");
        assert!(Fixture::entries(&fixture.0.0).is_empty());
    }

    #[test]
    fn retired_sweep_preserves_other_states_unrelated_entries_and_live_replacements() {
        let fixture = Fixture::new();
        let root = fixture.root();
        let other_root = fixture.0.0.join("other-state");
        let mut retired_paths = Vec::new();
        for path in [&root, &root, &other_root] {
            let writer = SnapshotWriter::acquire(path, MANIFESTS).unwrap();
            let retired = writer.retire_root().unwrap();
            retired_paths.push(retired.writer.root().to_owned());
            drop(retired);
        }
        // Model an interrupted removal that already unlinked its root lock.
        fs::remove_file(retired_paths[0].join(SNAPSHOT_LOCK_FILE)).unwrap();
        let malformed = fixture.0.0.join(format!(
            "{}not-a-generation",
            SnapshotWriter::retired_prefix(&root).unwrap()
        ));
        fs::create_dir(&malformed).unwrap();
        fs::write(malformed.join("sentinel"), b"unrelated").unwrap();
        let mut replacement = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        let published = Fixture::publish(&mut replacement, b"replacement state", STATE);
        SnapshotWriter::remove_retired_roots(&root).unwrap();
        assert!(!retired_paths[0].exists());
        assert!(!retired_paths[1].exists());
        assert!(retired_paths[2].exists());
        assert_eq!(fs::read(malformed.join("sentinel")).unwrap(), b"unrelated");
        assert_eq!(
            fs::read(published.path().join("payload")).unwrap(),
            b"replacement state"
        );
        assert!(matches!(
            SnapshotWriter::acquire(&root, MANIFESTS),
            Err(SnapshotError::Busy { .. })
        ));
        SnapshotWriter::remove_retired_roots(&other_root).unwrap();
        assert!(!retired_paths[2].exists());
        assert!(!other_root.exists());
    }

    #[test]
    fn retired_sweep_leaves_missing_parents_absent() {
        let fixture = Fixture::new();
        let root = fixture.0.0.join("missing/state");
        SnapshotWriter::remove_retired_roots(&root).unwrap();
        assert!(Fixture::entries(&fixture.0.0).is_empty());
    }

    #[cfg(unix)]
    #[test]
    fn retired_sweep_skips_an_unreadable_parent_only_when_the_root_is_absent() {
        let fixture = Fixture::new();
        let parent = fixture.0.0.join("search-only");
        let present = parent.join("present");
        fs::create_dir_all(&present).unwrap();
        let Some(_restriction) = SearchOnlyDirectory::new(&parent) else {
            return;
        };
        let missing = parent.join("missing");
        SnapshotWriter::remove_retired_roots(&missing).unwrap();
        assert!(!missing.exists());
        let error = SnapshotWriter::remove_retired_roots(&present).unwrap_err();
        assert!(
            matches!(error, SnapshotError::Io { operation, path, source }
            if operation == "list retired snapshots" && path == parent
                && source.kind() == ErrorKind::PermissionDenied)
        );
    }

    #[cfg(unix)]
    #[test]
    fn retired_sweep_does_not_delete_a_sibling_swapped_for_a_symlink() {
        let fixture = Fixture::new();
        let root = fixture.root();
        let retired = SnapshotWriter::acquire(&root, &[])
            .unwrap()
            .retire_root()
            .unwrap();
        let sibling = retired.writer.root().to_owned();
        drop(retired);
        assert!(fs::symlink_metadata(&sibling).unwrap().is_dir());
        let displaced = fixture.0.0.join("displaced");
        fs::rename(&sibling, &displaced).unwrap();
        let external = fixture.0.0.join("external");
        drop(SnapshotWriter::acquire(&external, &[]).unwrap());
        fs::write(external.join("payload"), b"external data").unwrap();
        unix_fs::symlink(&external, &sibling).unwrap();

        SnapshotWriter::remove_retired_root(&sibling).unwrap();

        assert!(fs::symlink_metadata(&sibling).unwrap().is_symlink());
        assert_eq!(
            fs::read(external.join("payload")).unwrap(),
            b"external data"
        );
        assert!(displaced.join(SNAPSHOT_LOCK_FILE).is_file());
    }

    #[test]
    fn retired_payload_removal_keeps_the_root_lock_until_the_tree_is_empty() {
        let fixture = Fixture::new();
        let root = fixture.root();
        let writer = SnapshotWriter::acquire(&root, &[]).unwrap();
        fs::create_dir(root.join("payload")).unwrap();
        fs::write(root.join("payload/data"), b"old snapshot").unwrap();
        let retired = writer.retire_root().unwrap();
        let sibling = retired.writer.root().to_owned();
        let lock = sibling.join(SNAPSHOT_LOCK_FILE);
        let original_lock = Handle::from_path(&lock).unwrap();

        retired.remove_contents().unwrap();

        assert_eq!(Fixture::entries(&sibling), vec![lock.to_owned()]);
        assert_eq!(original_lock, Handle::from_path(&lock).unwrap());
        let busy = SnapshotWriter::acquire(&sibling, &[]).unwrap_err();
        assert!(busy.is_busy());
        assert_eq!(
            busy.to_string(),
            format!("snapshot writer is busy at {}", lock.display())
        );
        let invalid = "not-a-generation".parse::<GenerationId>().unwrap_err();
        assert!(!invalid.is_busy());
        SnapshotWriter::remove_retired_roots(&root).unwrap();
        assert!(lock.is_file());
        drop(original_lock);
        retired.remove().unwrap();
        assert!(!sibling.exists());
    }

    #[cfg(unix)]
    #[test]
    fn retired_sweep_ignores_symlink_and_regular_file_siblings() {
        let fixture = Fixture::new();
        let root = fixture.root();
        let external = fixture.0.0.join("external");
        fs::create_dir(&external).unwrap();
        fs::write(external.join("sentinel"), b"preserve").unwrap();
        let prefix = SnapshotWriter::retired_prefix(&root).unwrap();
        let symlink = fixture
            .0
            .0
            .join(format!("{prefix}{}", GenerationId::generate().unwrap()));
        unix_fs::symlink(&external, &symlink).unwrap();
        // Unix socket paths are length-limited, so use a regular file as a non-directory sibling.
        let file = fixture
            .0
            .0
            .join(format!("{prefix}{}", GenerationId::generate().unwrap()));
        fs::write(&file, b"unrelated file").unwrap();
        SnapshotWriter::remove_retired_roots(&root).unwrap();
        assert!(
            fs::symlink_metadata(symlink)
                .unwrap()
                .file_type()
                .is_symlink()
        );
        assert_eq!(fs::read(external.join("sentinel")).unwrap(), b"preserve");
        assert_eq!(fs::read(file).unwrap(), b"unrelated file");
        assert!(!root.exists());
    }

    #[cfg(unix)]
    #[test]
    fn retirement_locks_and_removes_non_utf8_descendants() {
        let fixture = Fixture::new();
        let root = fixture.root();
        let writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        let nested = root.join(OsStr::from_bytes(b"nested-\xff"));
        fs::create_dir(&nested).unwrap();
        let lock_path = nested.join(SNAPSHOT_LOCK_FILE);
        let lock = File::create_new(&lock_path).unwrap();
        lock.try_lock().unwrap();
        let error = writer.retire_root().unwrap_err();
        assert!(matches!(error, SnapshotError::Busy { path } if path == lock_path));
        assert!(nested.is_dir());
        drop(lock);
        let writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        writer.retire_root().unwrap().remove().unwrap();
        assert!(Fixture::entries(&fixture.0.0).is_empty());
    }

    #[cfg(unix)]
    #[test]
    fn retirement_unlinks_symlinks_and_special_entries_without_following_them() {
        let fixture = Fixture::new();
        let root = fixture.root();
        let external = fixture.0.0.join("external");
        let external_writer = SnapshotWriter::acquire(&external, MANIFESTS).unwrap();
        let target = external.join("sentinel");
        fs::write(&target, b"preserve").unwrap();
        let writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        unix_fs::symlink(&external, root.join("directory-link")).unwrap();
        unix_fs::symlink(&target, root.join("file-link")).unwrap();
        fs::create_dir(root.join("nested")).unwrap();
        unix_fs::symlink(
            external.join(SNAPSHOT_LOCK_FILE),
            root.join("nested").join(SNAPSHOT_LOCK_FILE),
        )
        .unwrap();
        let socket = UnixListener::bind(root.join("socket")).unwrap();
        let retired = writer.retire_root().unwrap();
        let retired_path = retired.writer.root().to_owned();
        unix_fs::symlink(&target, retired_path.join("late-link")).unwrap();
        let _replacement = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        retired.remove().unwrap();
        assert!(!retired_path.exists());
        assert!(root.join(SNAPSHOT_LOCK_FILE).is_file());
        assert_eq!(fs::read(target).unwrap(), b"preserve");
        assert!(matches!(
            SnapshotWriter::acquire(external_writer.root(), MANIFESTS),
            Err(SnapshotError::Busy { .. })
        ));
        drop(socket);
    }

    #[cfg(unix)]
    #[test]
    fn temporary_cleanup_rejects_symlinks_before_removing_any_files() {
        let fixture = Fixture::new();
        fs::create_dir(fixture.root()).unwrap();
        let target = fixture.0.0.join("sentinel");
        fs::write(&target, b"preserve").unwrap();
        let first = fixture.root().join(format!(
            ".{STATE}.{}.tmp",
            GenerationId::generate().unwrap()
        ));
        let second = fixture.root().join(format!(
            ".{STATE}.{}.tmp",
            GenerationId::generate().unwrap()
        ));
        fs::write(&first, b"abandoned").unwrap();
        unix_fs::symlink(&target, &second).unwrap();
        assert!(matches!(
            SnapshotWriter::acquire(&fixture.root(), MANIFESTS),
            Err(SnapshotError::Path(OutputPathError::Symlink(_)))
        ));
        assert_eq!(fs::read(first).unwrap(), b"abandoned");
        assert_eq!(fs::read(target).unwrap(), b"preserve");
        assert!(fs::symlink_metadata(second).unwrap().is_symlink());
    }

    #[test]
    fn scratch_directories_and_manifest_temporary_names_are_caller_supplied() {
        let fixture = Fixture::new();
        let manifest = ManifestName::new("application.json").unwrap();
        fs::create_dir_all(fixture.root().join("diagnostics")).unwrap();
        fs::write(fixture.root().join("diagnostics/message"), b"scratch").unwrap();
        assert!(!SnapshotWriter::is_scratch(&fixture.root(), &[manifest], &[]).unwrap());
        assert!(
            SnapshotWriter::is_scratch(&fixture.root(), &[manifest], &["diagnostics"]).unwrap()
        );
        let temporary = fixture.root().join(format!(
            ".{manifest}.{}.tmp",
            GenerationId::generate().unwrap()
        ));
        fs::write(&temporary, b"partial").unwrap();
        assert!(
            SnapshotWriter::is_scratch(&fixture.root(), &[manifest], &["diagnostics"]).unwrap()
        );
        assert!(!SnapshotWriter::is_scratch(&fixture.root(), &[], &["diagnostics"]).unwrap());
        assert_eq!(fs::read(&temporary).unwrap(), b"partial");
        fs::write(fixture.root().join(manifest.as_str()), b"pointer").unwrap();
        assert!(
            !SnapshotWriter::is_scratch(&fixture.root(), &[manifest], &["diagnostics"]).unwrap()
        );
    }

    #[test]
    fn writer_contention_is_immediate_and_drop_releases_persistent_lock() {
        let fixture = Fixture::new();
        let writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        assert!(matches!(
            SnapshotWriter::acquire(&fixture.root(), MANIFESTS),
            Err(SnapshotError::Busy { .. })
        ));
        drop(writer);
        assert!(fixture.root().join(".snapshot.lock").is_file());
        let writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        assert!(matches!(
            SnapshotWriter::acquire(writer.root(), MANIFESTS),
            Err(SnapshotError::Busy { .. })
        ));
        drop(writer);
    }

    #[test]
    fn invalid_manifest_target_is_preserved() {
        let fixture = Fixture::new();
        let mut writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        let manifest = writer.root().join("manifest.json");
        fs::create_dir(&manifest).unwrap();
        fs::write(manifest.join("foreign"), b"preserve").unwrap();
        let stage = writer.begin().unwrap();
        assert!(matches!(
            stage.publish(b"pointer", INTEGRATION),
            Err(SnapshotError::Path(OutputPathError::NotFile(_)))
        ));
        assert_eq!(fs::read(manifest.join("foreign")).unwrap(), b"preserve");
        assert!(Fixture::entries(&writer.root().join("generations")).is_empty());
    }

    #[cfg(unix)]
    #[test]
    fn payload_symlinks_and_special_entries_preserve_the_previous_commit() {
        let fixture = Fixture::new();
        let external = fixture.0.0.join("external");
        fs::create_dir(&external).unwrap();
        fs::write(external.join("sentinel"), b"preserve").unwrap();
        let mut writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        let first = Fixture::publish(&mut writer, b"old", STATE);
        for kind in ["file-link", "directory-link", "socket"] {
            let stage = writer.begin().unwrap();
            let abandoned = stage.path().to_owned();
            fs::create_dir(stage.path().join("nested")).unwrap();
            let invalid = stage.path().join("nested/invalid");
            match kind {
                "file-link" => unix_fs::symlink(external.join("sentinel"), invalid).unwrap(),
                "directory-link" => unix_fs::symlink(&external, invalid).unwrap(),
                "socket" => {
                    // Relative binding in a child avoids SUN_LEN and process-wide cwd mutation.
                    let mut child =
                        ChildProcess::spawn(&fixture, "socket", "", &stage.path().join("nested"));
                    child.assert_success();
                }
                _ => unreachable!(),
            }
            assert!(matches!(
                stage.publish(b"new", STATE),
                Err(SnapshotError::Path(
                    OutputPathError::Symlink(_) | OutputPathError::SpecialEntry(_)
                ))
            ));
            assert!(!abandoned.exists());
            assert_eq!(
                fs::read(writer.root().join("state_manifest.toml")).unwrap(),
                b"old"
            );
            assert_eq!(fs::read(first.path().join("payload")).unwrap(), b"old");
            assert_eq!(fs::read(external.join("sentinel")).unwrap(), b"preserve");
            assert_eq!(Fixture::entries(&external).len(), 1);
        }
    }

    #[cfg(unix)]
    #[test]
    fn selected_root_symlinks_resolve_but_descendant_symlinks_are_rejected() {
        let fixture = Fixture::new();
        let external = fixture.0.0.join("external");
        fs::create_dir(&external).unwrap();
        fs::write(external.join("sentinel"), b"preserve").unwrap();
        unix_fs::symlink(&external, fixture.root()).unwrap();
        for root in [
            fixture.root(),
            fixture.root().join("."),
            fixture.root().join(""),
        ] {
            let writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
            assert_eq!(writer.root(), external);
        }
        assert_eq!(fs::read(external.join("sentinel")).unwrap(), b"preserve");
        assert!(external.join(SNAPSHOT_LOCK_FILE).is_file());
        fs::remove_file(fixture.root()).unwrap();
        fs::create_dir(fixture.root()).unwrap();

        let lock = fixture.root().join(".snapshot.lock");
        unix_fs::symlink(external.join("sentinel"), &lock).unwrap();
        assert!(matches!(
            SnapshotWriter::acquire(&fixture.root(), MANIFESTS),
            Err(SnapshotError::Path(OutputPathError::Symlink(_)))
        ));
        fs::remove_file(lock).unwrap();

        let mut writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        let generations = writer.root().join("generations");
        unix_fs::symlink(&external, &generations).unwrap();
        assert!(matches!(
            writer.begin(),
            Err(SnapshotError::Path(OutputPathError::Symlink(_)))
        ));
        fs::remove_file(&generations).unwrap();

        let manifest = writer.root().join("state_manifest.toml");
        unix_fs::symlink(external.join("sentinel"), &manifest).unwrap();
        let stage = writer.begin().unwrap();
        assert!(matches!(
            stage.publish(b"new", STATE),
            Err(SnapshotError::Path(OutputPathError::Symlink(_)))
        ));
        assert_eq!(fs::read(external.join("sentinel")).unwrap(), b"preserve");
        assert_eq!(Fixture::entries(&external).len(), 2);
        assert!(fs::symlink_metadata(manifest).unwrap().is_symlink());
    }
    #[test]
    fn committed_sync_failure_retains_visible_and_previous_generations() {
        let fixture = Fixture::new();
        let mut writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        let old = Fixture::publish(&mut writer, b"old", STATE);
        let mut stage = writer.begin().unwrap();
        fs::write(stage.path().join("payload"), b"new").unwrap();
        let expected_id = stage.id().to_string();
        stage.fault = Some(PublicationStep::DirectoryCompletion);
        let error = stage.publish(expected_id.as_bytes(), STATE).unwrap_err();
        let SnapshotError::Committed { generation, source } = error else {
            panic!("post-commit failure must identify the visible generation: {error}");
        };
        assert_eq!(generation.id().as_str(), expected_id);
        assert!(source.source().unwrap().is::<IoError>());
        assert_eq!(
            fs::read_to_string(writer.root().join(STATE_MANIFEST_FILE)).unwrap(),
            expected_id
        );
        assert_eq!(fs::read(generation.path().join("payload")).unwrap(), b"new");
        assert_eq!(fs::read(old.path().join("payload")).unwrap(), b"old");
        assert_eq!(
            Fixture::entries(&writer.root().join(GENERATIONS_DIRECTORY)).len(),
            2
        );
    }

    #[test]
    fn pruning_keeps_current_previous_other_manifest_and_pinned_readers() {
        let fixture = Fixture::new();
        let mut writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        let state = Fixture::publish(&mut writer, b"state", STATE);
        let pinned = Fixture::publish(&mut writer, b"pinned reader", INTEGRATION);
        let unneeded = Fixture::publish(&mut writer, b"unneeded", INTEGRATION);
        let previous = Fixture::publish(&mut writer, b"previous", INTEGRATION);
        let current = Fixture::publish(&mut writer, b"current", INTEGRATION);
        writer
            .prune_generations([state.id(), pinned.id(), previous.id(), current.id()])
            .unwrap();
        assert!(!unneeded.path().exists());
        for generation in [&state, &pinned, &previous, &current] {
            assert!(generation.path().is_dir());
        }
        assert_eq!(
            fs::read(writer.root().join(STATE_MANIFEST_FILE)).unwrap(),
            b"state"
        );
        assert_eq!(
            fs::read(writer.root().join(INTEGRATION_MANIFEST_FILE)).unwrap(),
            b"current"
        );
        writer
            .prune_generations([state.id(), previous.id(), current.id()])
            .unwrap();
        assert!(!pinned.path().exists());
    }

    #[test]
    fn pruning_ignores_strays_and_checks_missing_keep_before_removal() {
        let fixture = Fixture::new();
        let mut writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        let old = Fixture::publish(&mut writer, b"old", STATE);
        let current = Fixture::publish(&mut writer, b"current", STATE);
        let generations = writer.root().join(GENERATIONS_DIRECTORY);
        let stray_directory = generations.join("unknown");
        let stray_file = generations.join(".DS_Store");
        fs::create_dir(&stray_directory).unwrap();
        fs::write(&stray_file, b"unmanaged").unwrap();
        let missing = GenerationId::generate().unwrap();
        assert!(matches!(
            writer.prune_generations([current.id(), &missing]),
            Err(SnapshotError::Path(OutputPathError::NotDirectory(_)))
        ));
        assert!(old.path().is_dir());
        assert!(current.path().is_dir());
        let removed = writer.prune_generations([current.id()]).unwrap();
        assert_eq!(removed.len(), 1);
        assert_eq!(&removed[0], old.id());
        assert!(!old.path().exists());
        assert!(current.path().is_dir());
        assert!(stray_directory.is_dir());
        assert_eq!(fs::read(stray_file).unwrap(), b"unmanaged");
        assert!(writer.prune_generations([current.id()]).unwrap().is_empty());
    }

    #[cfg(unix)]
    #[test]
    fn pruning_rejects_nested_symlinks_without_removing_any_generation() {
        let fixture = Fixture::new();
        let mut writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        let old = Fixture::publish(&mut writer, b"old", STATE);
        let current = Fixture::publish(&mut writer, b"current", STATE);
        let external = fixture.0.0.join("external");
        fs::write(&external, b"preserve").unwrap();
        unix_fs::symlink(&external, old.path().join("link")).unwrap();
        assert!(matches!(
            writer.prune_generations([current.id()]),
            Err(SnapshotError::Path(OutputPathError::Symlink(_)))
        ));
        assert!(old.path().is_dir());
        assert!(current.path().is_dir());
        assert_eq!(fs::read(external).unwrap(), b"preserve");
    }

    #[test]
    fn scratch_classification_accepts_interrupted_first_save_without_writes() {
        let fixture = Fixture::new();
        assert!(SnapshotWriter::is_scratch(&fixture.root(), MANIFESTS, &["logs"]).unwrap());
        assert!(!fixture.root().exists());
        fs::create_dir_all(fixture.root().join("logs")).unwrap();
        fs::write(fixture.root().join("logs/run.jsonl"), b"log").unwrap();
        assert!(SnapshotWriter::is_scratch(&fixture.root(), MANIFESTS, &["logs"]).unwrap());
        assert!(!fixture.root().join(SNAPSHOT_LOCK_FILE).exists());
        let mut writer = SnapshotWriter::acquire(&fixture.root(), MANIFESTS).unwrap();
        // Simulate SIGKILL leaving an otherwise valid, unpublished generation.
        let orphan = GenerationId::generate().unwrap();
        let generation = writer
            .root()
            .join(GENERATIONS_DIRECTORY)
            .join(orphan.as_str());
        fs::create_dir_all(&generation).unwrap();
        fs::write(generation.join("partial"), b"payload").unwrap();
        let temporary = writer
            .root()
            .join(format!(".{}.{}.tmp", STATE_MANIFEST_FILE, orphan));
        fs::write(&temporary, b"incomplete manifest").unwrap();
        let before = Fixture::entries(writer.root());
        assert!(SnapshotWriter::is_scratch(writer.root(), MANIFESTS, &["logs"]).unwrap());
        assert_eq!(Fixture::entries(writer.root()), before);
        assert_eq!(fs::read(temporary).unwrap(), b"incomplete manifest");
        let unrelated = writer.root().join("user-data");
        fs::write(&unrelated, b"preserve").unwrap();
        assert!(!SnapshotWriter::is_scratch(writer.root(), MANIFESTS, &["logs"]).unwrap());
        fs::remove_file(unrelated).unwrap();
        Fixture::publish(&mut writer, b"committed", STATE);
        assert!(!SnapshotWriter::is_scratch(writer.root(), MANIFESTS, &["logs"]).unwrap());
    }

    #[cfg(unix)]
    #[test]
    fn scratch_classification_rejects_managed_symlinks_and_supports_system_ancestors() {
        let fixture = Fixture::new();
        let private = fixture.0.0.join("private");
        fs::create_dir(&private).unwrap();
        let alias = fixture.0.0.join("var");
        unix_fs::symlink(&private, &alias).unwrap();
        let root = alias.join("managed");
        assert!(SnapshotWriter::is_scratch(&root, MANIFESTS, &["logs"]).unwrap());
        assert!(!private.join("managed").exists());
        let writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        assert_eq!(
            writer.root(),
            fs::canonicalize(&private).unwrap().join("managed")
        );
        drop(writer);
        let external = fixture.0.0.join("external");
        fs::create_dir(&external).unwrap();
        unix_fs::symlink(&external, root.join("logs")).unwrap();
        assert!(matches!(
            SnapshotWriter::is_scratch(&root, MANIFESTS, &["logs"]),
            Err(SnapshotError::Path(OutputPathError::Symlink(_)))
        ));
        assert_eq!(Fixture::entries(&external).len(), 0);
    }
    #[cfg(unix)]
    #[test]
    fn search_only_higher_ancestor_allows_snapshot_and_directory_publication() {
        let fixture = Fixture::new();
        let gate = fixture.0.0.join("gate");
        let work = gate.join("work");
        fs::create_dir_all(&work).unwrap();
        fs::set_permissions(&work, PermissionsExt::from_mode(0o700)).unwrap();
        let root = work.join("snapshot");
        let mut writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        let old = Fixture::publish(&mut writer, b"old checkpoint", STATE);
        drop(writer);
        let site = work.join("site");
        fs::create_dir(&site).unwrap();
        fs::write(site.join("index.html"), b"old site").unwrap();
        let Some(restriction) = SearchOnlyDirectory::new(&gate) else {
            return;
        };
        let snapshot = SnapshotWriter::acquire(&root, MANIFESTS);
        let directory = DirectoryPublication::new(&site);
        // Restore permissions before assertions so failure never obstructs fixture cleanup.
        drop(restriction);
        let mut writer = snapshot.unwrap();
        let publication = directory.unwrap();
        assert_eq!(
            fs::read(root.join(STATE_MANIFEST_FILE)).unwrap(),
            b"old checkpoint"
        );
        assert_eq!(fs::read(site.join("index.html")).unwrap(), b"old site");
        let next = Fixture::publish(&mut writer, b"next checkpoint", STATE);
        assert_eq!(
            fs::read(next.path().join("payload")).unwrap(),
            b"next checkpoint"
        );
        assert_eq!(
            fs::read(old.path().join("payload")).unwrap(),
            b"old checkpoint"
        );
        fs::create_dir(publication.staging_path()).unwrap();
        fs::write(publication.staging_path().join("index.html"), b"next site").unwrap();
        publication.publish().unwrap();
        assert_eq!(fs::read(site.join("index.html")).unwrap(), b"next site");
    }

    #[cfg(unix)]
    #[test]
    fn existing_snapshot_beneath_search_only_parent_can_acquire_and_publish() {
        let fixture = Fixture::new();
        let parent = fixture.0.0.join("search-only");
        let root = parent.join("snapshot");
        let mut writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        let first = Fixture::publish(&mut writer, b"first", STATE);
        drop(writer);
        let Some(_restriction) = SearchOnlyDirectory::new(&parent) else {
            return;
        };
        let mut writer = SnapshotWriter::acquire(&root, MANIFESTS).unwrap();
        let next = Fixture::publish(&mut writer, b"next", STATE);
        assert_eq!(fs::read(first.path().join("payload")).unwrap(), b"first");
        assert_eq!(fs::read(next.path().join("payload")).unwrap(), b"next");
        assert_eq!(fs::read(root.join(STATE_MANIFEST_FILE)).unwrap(), b"next");
    }
}
