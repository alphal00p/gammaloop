//! Publish immutable generations by replacing one manifest under an advisory writer lock.
//!
//! Callers own destination authorization, payload serialization, and manifest validation. Check
//! read-only and overwrite policies **before** acquiring a writer, then recheck under its lock.
//! Payload writers must check their writes, finish and flush all application buffers, and stop
//! modifying the generation before calling `publish`. Publication opens and synchronizes every
//! regular payload file, then writes and synchronizes the manifest before replacing its pointer.
//! Symlinks and special payload entries are rejected. Files must remain writable for
//! synchronization.
//!
//! The root itself, `generations`, the lock, and an existing manifest must not be symlinks.
//! Ancestor symlinks are resolved once when acquiring the root. Callers must keep those directories
//! under their control: this is a cooperating-writer contract, not protection against concurrent
//! hostile directory replacement. Readers pin one manifest and retain its generation throughout
//! their operation. Published generations are retained indefinitely; no garbage collection occurs.
//! The filesystem must support locking, file synchronization, and atomic replacement. Checked
//! `File::sync_all` completion surfaces delayed write errors reported by the filesystem before
//! commit; `File::flush` and dropping a handle alone do not establish that completion. Successful
//! publication provides atomic visibility across ordinary I/O failures and process interruption,
//! assuming callers have checked their writes and flushes. Synchronization costs a filesystem walk,
//! an open, and a synchronization operation per payload file, plus manifest synchronization; it may
//! wait for dirty data to reach storage. Directories are not synchronized, so this is **not** a
//! power-loss durability guarantee.
//!
//! # Example
//!
//! Both domain transactions use the same owner. These minimal manifests illustrate the pointer;
//! real callers serialize their complete domain descriptor.
//!
//! ```
//! use gammalooprs::utils::snapshot::{GenerationId, ManifestName, SnapshotWriter};
//! use std::env;
//! use std::error::Error;
//! use std::fs;
//! use std::process;
//! use std::time::{SystemTime, UNIX_EPOCH};
//!
//! # fn main() -> Result<(), Box<dyn Error>> {
//! let nonce = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
//! let root = env::temp_dir().join(format!("snapshot-example-{}-{nonce}", process::id()));
//! // Application policy has authorized this fresh destination.
//! let mut writer = SnapshotWriter::acquire(&root)?;
//! for name in [ManifestName::State, ManifestName::Integration] {
//!     let staged = writer.begin()?;
//!     fs::write(staged.path().join("payload.bin"), b"complete payload")?;
//!     let id: GenerationId = staged.id().as_str().parse()?;
//!     let manifest = match name {
//!         ManifestName::State => format!("generation = \"{id}\"\n"),
//!         ManifestName::Integration => format!("{{\"generation\":\"{id}\"}}\n"),
//!     };
//!     let published = staged.publish(manifest.as_bytes(), name)?;
//!     assert_eq!(published.id(), &id);
//!     assert_eq!(published.path(), writer.root().join("generations").join(id.as_str()));
//!     assert_eq!(fs::read(published.path().join("payload.bin"))?, b"complete payload");
//! }
//! drop(writer); // Releases the lock; the persistent lock file remains.
//! fs::remove_dir_all(root)?;
//! # Ok(())
//! # }
//! ```

use serde::{Deserialize, Serialize};
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::fs;
use std::fs::{File, TryLockError};
use std::io::{Error as IoError, ErrorKind, Write};
use std::path::{Path, PathBuf};
use std::process;
use std::str::FromStr;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::{SystemTime, SystemTimeError, UNIX_EPOCH};
use thiserror::Error;

const CREATE_ATTEMPTS: usize = 128;
static NEXT_GENERATION: AtomicU64 = AtomicU64::new(0);

/// A single component in the format `gen-<8 hex pid>-<32 hex nanoseconds>-<16 hex counter>`.
///
/// Hexadecimal digits are lowercase. Only the writer generates new identities; parsing and Serde
/// decoding validate references read from manifests. See the [module example](self).
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(try_from = "String")]
pub struct GenerationId(String);

impl GenerationId {
    /// Borrow the validated component, as shown in the [module example](self).
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }

    fn generate() -> Result<Self, SnapshotError> {
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

/// The two supported logical commit pointers. See the [module example](self).
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum ManifestName {
    /// Saved session pointer, exactly `state_manifest.toml`.
    State,
    /// Completed integration checkpoint pointer, exactly `manifest.json`.
    Integration,
}

impl ManifestName {
    fn as_str(self) -> &'static str {
        match self {
            Self::State => "state_manifest.toml",
            Self::Integration => "manifest.json",
        }
    }
}

/// Publication failures, retaining I/O sources and distinguishing immediate writer contention.
///
/// Every error from publication precedes its commit point. Cleanup is best effort and never
/// replaces the primary error. The [module example](self) demonstrates propagation.
#[derive(Debug, Error)]
pub enum SnapshotError {
    /// A filesystem operation failed, including unsupported locking or replacement.
    #[error("cannot {operation} at {}: {source}", path.display())]
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
    /// A manifest supplied a component outside the internally generated format.
    #[error("invalid snapshot generation {value:?}")]
    InvalidGeneration {
        /// The rejected component.
        value: String,
    },
    /// An existing entry violates the directory, regular-file, or symlink policy.
    #[error("expected {expected} at {}", path.display())]
    InvalidPath {
        /// The rejected entry.
        path: PathBuf,
        /// The required entry type.
        expected: &'static str,
    },
    /// The system clock cannot supply a timestamp relative to the Unix epoch.
    #[error("cannot generate a snapshot identity: {0}")]
    Clock(#[from] SystemTimeError),
}

impl SnapshotError {
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
/// The lock file is never removed. See the [module example](self).
#[derive(Debug)]
pub struct SnapshotWriter {
    root: PathBuf,
    _lock: File,
}

impl SnapshotWriter {
    /// Create an authorized root if needed and acquire its persistent `.snapshot.lock`.
    ///
    /// The returned root is absolute, with ancestors resolved. This method performs no application
    /// authorization: callers must check their destination policy before calling it and recheck
    /// under the acquired lock. See the [module example](self).
    ///
    /// # Errors
    /// Returns [`SnapshotError::Busy`] immediately on contention. Invalid entry types, symlinks,
    /// I/O errors, and unsupported filesystem operations are errors.
    pub fn acquire(root: &Path) -> Result<Self, SnapshotError> {
        // Remove trailing separators and `.` so they cannot hide a symlink at the root entry.
        let root: PathBuf = root.components().collect();
        Self::ensure_directory(&root)?;
        let root = fs::canonicalize(&root)
            .map_err(|source| SnapshotError::io("resolve snapshot root", &root, source))?;
        let lock_path = root.join(".snapshot.lock");
        Self::check_regular_or_missing(&lock_path)?;
        let lock = File::options()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(&lock_path)
            .map_err(|source| SnapshotError::io("open writer lock", &lock_path, source))?;
        match lock.try_lock() {
            Ok(()) => {}
            Err(TryLockError::WouldBlock) => {
                return Err(SnapshotError::Busy { path: lock_path });
            }
            Err(TryLockError::Error(source)) => {
                return Err(SnapshotError::io("lock snapshot root", &lock_path, source));
            }
        }
        Self::check_regular_or_missing(&lock_path)?;
        Ok(Self { root, _lock: lock })
    }

    /// Borrow the resolved logical root; see the [module example](self).
    #[must_use]
    pub fn root(&self) -> &Path {
        &self.root
    }

    /// Exclusively create a fresh unpublished generation, retrying identity collisions.
    ///
    /// Its mutable borrow keeps the writer locked throughout staging and publication. See the
    /// [module example](self).
    ///
    /// # Errors
    /// Returns entry-validation, clock, or I/O errors without touching published generations.
    pub fn begin(&mut self) -> Result<StagedGeneration<'_>, SnapshotError> {
        Self::check_directory(&self.root)?;
        let generations = self.root.join("generations");
        Self::ensure_directory(&generations)?;
        for attempt in 0..CREATE_ATTEMPTS {
            let id = GenerationId::generate()?;
            let path = generations.join(id.as_str());
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

    fn ensure_directory(path: &Path) -> Result<(), SnapshotError> {
        match fs::symlink_metadata(path) {
            Ok(_) => Self::check_directory(path),
            Err(source) if source.kind() == ErrorKind::NotFound => {
                fs::create_dir_all(path).map_err(|source| {
                    SnapshotError::io("create snapshot directory", path, source)
                })?;
                Self::check_directory(path)
            }
            Err(source) => Err(SnapshotError::io(
                "inspect snapshot directory",
                path,
                source,
            )),
        }
    }

    fn check_directory(path: &Path) -> Result<(), SnapshotError> {
        let metadata = fs::symlink_metadata(path)
            .map_err(|source| SnapshotError::io("inspect snapshot directory", path, source))?;
        if !metadata.is_dir() {
            return Err(SnapshotError::InvalidPath {
                path: path.to_owned(),
                expected: "a directory without a symlink",
            });
        }
        Ok(())
    }

    fn check_regular_or_missing(path: &Path) -> Result<(), SnapshotError> {
        match fs::symlink_metadata(path) {
            Ok(metadata) if metadata.is_file() => Ok(()),
            Ok(_) => Err(SnapshotError::InvalidPath {
                path: path.to_owned(),
                expected: "a regular file without a symlink, or no entry",
            }),
            Err(source) if source.kind() == ErrorKind::NotFound => Ok(()),
            Err(source) => Err(SnapshotError::io("inspect snapshot file", path, source)),
        }
    }
}

/// Owns one unpublished generation and its manifest temporary file while borrowing the writer.
///
/// Dropping it removes only those unpublished entries, on a best effort basis. See the
/// [module example](self).
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
    /// Borrow the fresh payload directory; see the [module example](self).
    #[must_use]
    pub fn path(&self) -> &Path {
        self.generation().path()
    }

    /// Borrow the internally generated identity for the manifest; see the [module example](self).
    #[must_use]
    pub fn id(&self) -> &GenerationId {
        self.generation().id()
    }

    /// Synchronize all regular payload files, then write and synchronize a unique manifest
    /// temporary file, close it, and atomically replace the pointer.
    ///
    /// The rename is the commit point. Immediately afterward, cleanup is disarmed and this method
    /// returns the committed generation, with no further fallible work. It never unlinks an old
    /// manifest first and never deletes prior committed generations. See the
    /// [module example](self).
    ///
    /// # Errors
    /// All errors precede the commit point, leaving the old manifest untouched. Symlinks, special
    /// payload entries, failed writes/flushes/synchronization, unsupported operations, and failed
    /// replacement are reported with context. Callers must finish and check their writes and
    /// application-buffer flushes before invoking this method, and stop modifying the generation.
    /// Payload files must remain writable. This method synchronizes files, not directories, and
    /// does not promise power-loss durability. See the [module contract](self) for the I/O cost.
    pub fn publish(
        mut self,
        manifest: &[u8],
        name: ManifestName,
    ) -> Result<PublishedGeneration, SnapshotError> {
        SnapshotWriter::check_directory(self.writer.root())?;
        SnapshotWriter::check_directory(&self.writer.root().join("generations"))?;
        SnapshotWriter::check_directory(self.path())?;
        self.complete_payload()?;
        let manifest_path = self.writer.root().join(name.as_str());
        SnapshotWriter::check_regular_or_missing(&manifest_path)?;
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
        SnapshotWriter::check_regular_or_missing(&manifest_path)?;
        #[cfg(test)]
        self.injected_error(PublicationStep::ManifestRename)
            .map_err(|source| SnapshotError::io("publish manifest", &manifest_path, source))?;
        #[cfg(test)]
        if let Some(barrier) = &self.barrier {
            barrier.wait(PublicationStep::ManifestRename);
        }
        fs::rename(temp, &manifest_path)
            .map_err(|source| SnapshotError::io("publish manifest", &manifest_path, source))?;
        let published = self
            .generation
            .take()
            .expect("a staged generation is present until its successful publication");
        self.manifest_temp = None;
        #[cfg(test)]
        if let Some(barrier) = &self.barrier {
            barrier.wait(PublicationStep::Committed);
        }
        Ok(published)
    }

    fn complete_payload(&self) -> Result<(), SnapshotError> {
        let mut directories = vec![self.path().to_owned()];
        while let Some(directory) = directories.pop() {
            let entries = fs::read_dir(&directory).map_err(|source| {
                SnapshotError::io("read payload directory", &directory, source)
            })?;
            for entry in entries {
                let entry = entry.map_err(|source| {
                    SnapshotError::io("read payload directory entry", &directory, source)
                })?;
                let path = entry.path();
                let kind = entry
                    .file_type()
                    .map_err(|source| SnapshotError::io("inspect payload entry", &path, source))?;
                if kind.is_dir() {
                    directories.push(path);
                } else if kind.is_file() {
                    self.complete_payload_file(&path)?;
                } else {
                    return Err(SnapshotError::InvalidPath {
                        path,
                        expected: "a regular payload file or directory without symlinks",
                    });
                }
            }
        }
        Ok(())
    }

    fn complete_payload_file(&self, path: &Path) -> Result<(), SnapshotError> {
        // Windows synchronization requires a writable handle; opening does not truncate the file.
        let file = File::options()
            .write(true)
            .open(path)
            .map_err(|source| SnapshotError::io("open payload for completion", path, source))?;
        let completion = file.sync_all();
        #[cfg(test)]
        let completion =
            completion.and_then(|()| self.injected_error(PublicationStep::PayloadCompletion));
        completion.map_err(|source| SnapshotError::io("complete payload file", path, source))
    }

    fn generation(&self) -> &PublishedGeneration {
        self.generation
            .as_ref()
            .expect("a staged generation is present until its successful publication")
    }

    fn create_manifest_temp(&self, name: ManifestName) -> Result<(PathBuf, File), SnapshotError> {
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
        if let Some(temp) = &self.manifest_temp {
            let _ = fs::remove_file(temp);
        }
        if let Some(generation) = &self.generation {
            let _ = fs::remove_dir_all(generation.path());
        }
    }
}

/// A committed generation, retained independently of the writer and future publications.
///
/// Dropping this descriptor does not delete its files. See the [module example](self).
#[derive(Debug, Eq, PartialEq)]
pub struct PublishedGeneration {
    id: GenerationId,
    path: PathBuf,
}

impl PublishedGeneration {
    /// Borrow the committed identity; see the [module example](self).
    #[must_use]
    pub fn id(&self) -> &GenerationId {
        &self.id
    }

    /// Borrow the retained payload directory; see the [module example](self).
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
    use std::env;
    use std::error::Error as StdError;
    #[cfg(unix)]
    use std::os::unix::fs as unix_fs;
    #[cfg(unix)]
    use std::os::unix::fs::MetadataExt;
    #[cfg(unix)]
    use std::os::unix::net::UnixListener;
    #[cfg(unix)]
    use std::os::unix::process::ExitStatusExt;
    use std::process::{Child, Command, ExitStatus};
    use std::thread;
    use std::time::{Duration, Instant};

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
            let control = fixture.0.join(action);
            fs::create_dir(&control).unwrap();
            let log = File::create_new(control.join("child.log")).unwrap();
            // Both streams deliberately share one descriptor for diagnostics on assertion failure.
            let stderr = log.try_clone().unwrap();
            let child = Command::new(env::current_exe().unwrap())
                .current_dir(directory)
                .args([
                    "--exact",
                    "utils::snapshot::tests::subprocess_entry",
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

    struct Fixture(PathBuf);

    impl Fixture {
        fn new() -> Self {
            let parent = env::temp_dir();
            let generation = GenerationId::generate().unwrap();
            let path = parent.join(format!("gammaloop-snapshot-test-{generation}"));
            fs::create_dir(&path).unwrap();
            Self(path)
        }

        fn root(&self) -> PathBuf {
            self.0.join("snapshot")
        }

        fn publish(
            writer: &mut SnapshotWriter,
            marker: &[u8],
            name: ManifestName,
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
            stage
                .publish(manifest.as_bytes(), ManifestName::Integration)
                .unwrap()
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
            let mut writer = SnapshotWriter::acquire(&root).unwrap();
            let first = Self::publish_pointer(&mut writer, b"A", None);
            drop(writer);
            let lock_path = root.join(".snapshot.lock");
            #[cfg(unix)]
            let original_lock = fs::metadata(&lock_path).unwrap();

            let mut interrupted = ChildProcess::spawn(&fixture, action, "", &fixture.0);
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

            let mut contender = ChildProcess::spawn(&fixture, "busy", "", &fixture.0);
            contender.assert_success();
            interrupted.kill();
            // SIGKILL skips all destructors: unpublished payloads and the pre-rename temp survive.
            assert_eq!(&Self::read_pointer(&root), expected_id);
            Self::assert_payload(first.path(), b"A");
            Self::assert_payload(&second, b"B");
            assert!(temporaries.iter().all(|path| path.is_file()));
            assert!(fs::symlink_metadata(&lock_path).unwrap().is_file());

            let mut recovery = ChildProcess::spawn(&fixture, "recover", expected, &fixture.0);
            recovery.assert_success();
            let third_id = Self::read_pointer(&root);
            assert_ne!(&third_id, first.id());
            assert_ne!(third_id, second_id);
            Self::assert_payload(&root.join("generations").join(third_id.as_str()), b"C");
            Self::assert_payload(first.path(), b"A");
            Self::assert_payload(&second, b"B");
            assert!(temporaries.iter().all(|path| path.is_file()));
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

    impl Drop for Fixture {
        fn drop(&mut self) {
            fs::remove_dir_all(&self.0).unwrap();
        }
    }

    #[test]
    #[ignore = "subprocess protocol entry; invoked by the bounded interruption tests"]
    fn subprocess_entry() {
        // Dispatch exists only in the test binary; production has no environment fault switch.
        let action = env::var("GAMMALOOP_SNAPSHOT_TEST_ACTION").unwrap();
        let root = PathBuf::from(env::var_os("GAMMALOOP_SNAPSHOT_TEST_ROOT").unwrap());
        let control = PathBuf::from(env::var_os("GAMMALOOP_SNAPSHOT_TEST_CONTROL").unwrap());
        match action.as_str() {
            "before" | "after" => {
                let mut writer = SnapshotWriter::acquire(&root).unwrap();
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
                    SnapshotWriter::acquire(&root),
                    Err(SnapshotError::Busy { .. })
                ));
            }
            "recover" => {
                let mut writer = SnapshotWriter::acquire(&root).unwrap();
                let expected = env::var("GAMMALOOP_SNAPSHOT_TEST_EXPECTED").unwrap();
                let pinned = Fixture::read_pointer(&root);
                Fixture::assert_payload(
                    &root.join("generations").join(pinned.as_str()),
                    expected.as_bytes(),
                );
                Fixture::publish_pointer(&mut writer, b"C", None);
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
    fn first_and_repeated_publication_retain_pinned_generations() {
        for name in [ManifestName::State, ManifestName::Integration] {
            let fixture = Fixture::new();
            let mut writer = SnapshotWriter::acquire(&fixture.root()).unwrap();
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
            let mut writer = SnapshotWriter::acquire(&fixture.root()).unwrap();
            let first = Fixture::publish(&mut writer, b"old", ManifestName::State);
            let root_entries = Fixture::entries(writer.root());
            let mut stage = writer.begin().unwrap();
            let abandoned = stage.path().to_owned();
            fs::write(stage.path().join("payload"), b"new").unwrap();
            stage.fault = Some(fault);
            let error = stage.publish(b"new", ManifestName::State).unwrap_err();
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
        let mut writer = SnapshotWriter::acquire(&fixture.root()).unwrap();
        let mut stage = writer.begin().unwrap();
        stage.fault = Some(PublicationStep::ManifestRename);
        assert!(stage.publish(b"first", ManifestName::Integration).is_err());
        assert!(!writer.root().join("manifest.json").exists());
        assert!(Fixture::entries(&writer.root().join("generations")).is_empty());
        assert_eq!(Fixture::entries(writer.root()).len(), 2);
    }

    #[test]
    fn dropping_stage_removes_only_its_unpublished_payload() {
        let fixture = Fixture::new();
        let mut writer = SnapshotWriter::acquire(&fixture.root()).unwrap();
        let committed = Fixture::publish(&mut writer, b"old", ManifestName::Integration);
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
        Fixture::publish(&mut writer, b"next", ManifestName::Integration);
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
    fn writer_contention_is_immediate_and_drop_releases_persistent_lock() {
        let fixture = Fixture::new();
        let writer = SnapshotWriter::acquire(&fixture.root()).unwrap();
        assert!(matches!(
            SnapshotWriter::acquire(&fixture.root()),
            Err(SnapshotError::Busy { .. })
        ));
        drop(writer);
        assert!(fixture.root().join(".snapshot.lock").is_file());
        let writer = SnapshotWriter::acquire(&fixture.root()).unwrap();
        assert!(matches!(
            SnapshotWriter::acquire(writer.root()),
            Err(SnapshotError::Busy { .. })
        ));
        drop(writer);
    }

    #[test]
    fn invalid_manifest_target_is_preserved() {
        let fixture = Fixture::new();
        let mut writer = SnapshotWriter::acquire(&fixture.root()).unwrap();
        let manifest = writer.root().join("manifest.json");
        fs::create_dir(&manifest).unwrap();
        fs::write(manifest.join("foreign"), b"preserve").unwrap();
        let stage = writer.begin().unwrap();
        assert!(matches!(
            stage.publish(b"pointer", ManifestName::Integration),
            Err(SnapshotError::InvalidPath { .. })
        ));
        assert_eq!(fs::read(manifest.join("foreign")).unwrap(), b"preserve");
        assert!(Fixture::entries(&writer.root().join("generations")).is_empty());
    }

    #[cfg(unix)]
    #[test]
    fn payload_symlinks_and_special_entries_preserve_the_previous_commit() {
        let fixture = Fixture::new();
        let external = fixture.0.join("external");
        fs::create_dir(&external).unwrap();
        fs::write(external.join("sentinel"), b"preserve").unwrap();
        let mut writer = SnapshotWriter::acquire(&fixture.root()).unwrap();
        let first = Fixture::publish(&mut writer, b"old", ManifestName::State);
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
                stage.publish(b"new", ManifestName::State),
                Err(SnapshotError::InvalidPath { .. })
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
    fn symlink_entries_are_rejected_without_touching_their_targets() {
        let fixture = Fixture::new();
        let external = fixture.0.join("external");
        fs::create_dir(&external).unwrap();
        fs::write(external.join("sentinel"), b"preserve").unwrap();
        unix_fs::symlink(&external, fixture.root()).unwrap();
        for root in [
            fixture.root(),
            fixture.root().join("."),
            fixture.root().join(""),
        ] {
            assert!(matches!(
                SnapshotWriter::acquire(&root),
                Err(SnapshotError::InvalidPath { .. })
            ));
        }
        assert_eq!(Fixture::entries(&external).len(), 1);
        fs::remove_file(fixture.root()).unwrap();
        fs::create_dir(fixture.root()).unwrap();

        let lock = fixture.root().join(".snapshot.lock");
        unix_fs::symlink(external.join("sentinel"), &lock).unwrap();
        assert!(matches!(
            SnapshotWriter::acquire(&fixture.root()),
            Err(SnapshotError::InvalidPath { .. })
        ));
        fs::remove_file(lock).unwrap();

        let mut writer = SnapshotWriter::acquire(&fixture.root()).unwrap();
        let generations = writer.root().join("generations");
        unix_fs::symlink(&external, &generations).unwrap();
        assert!(matches!(
            writer.begin(),
            Err(SnapshotError::InvalidPath { .. })
        ));
        fs::remove_file(&generations).unwrap();

        let manifest = writer.root().join("state_manifest.toml");
        unix_fs::symlink(external.join("sentinel"), &manifest).unwrap();
        let stage = writer.begin().unwrap();
        assert!(matches!(
            stage.publish(b"new", ManifestName::State),
            Err(SnapshotError::InvalidPath { .. })
        ));
        assert_eq!(fs::read(external.join("sentinel")).unwrap(), b"preserve");
        assert_eq!(Fixture::entries(&external).len(), 1);
        assert!(fs::symlink_metadata(manifest).unwrap().is_symlink());
    }
}
