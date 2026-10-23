//! Complete-directory replacement with explicit rollback and recovery artifacts.

use std::error::Error as StdError;
use std::fs::{self, File, TryLockError};
#[cfg(test)]
use std::io::Error as IoError;
use std::io::ErrorKind;
use std::path::{Path, PathBuf};
use std::time::SystemTimeError;

use crate::snapshot::CREATE_ATTEMPTS;
use crate::{GenerationId, OutputPath, OutputPathError};
use thiserror::Error;

const PUBLICATION_LOCK_PREFIX: &str = ".publication-lock-";
const PUBLICATION_PREFIX: &str = ".publication-";

/// Failures before, during, or after directory replacement. See [`DirectoryPublication`].
#[derive(Debug, Error)]
pub enum DirectoryPublicationError {
    /// Path validation, recovery, or staging failed.
    #[error(transparent)]
    Path(#[from] OutputPathError),
    /// A unique staging identity could not be generated.
    #[error("cannot generate a directory publication identity")]
    Clock(#[from] SystemTimeError),
    /// Another cooperating publisher owns this destination.
    #[error("directory publisher is busy at {path}")]
    Busy {
        /// The persistent destination-specific lock.
        path: PathBuf,
    },
    /// The requested destination has a reserved publication-lock filename.
    #[error("directory destination {0} is reserved for publication locking")]
    ReservedDestination(PathBuf),
    /// Recovery cannot choose a previous tree safely; all journals remain untouched.
    #[error(
        "cannot recover directory {destination}: {reason}; \
         inspect and resolve retained publication journals before retrying: {journals:?}"
    )]
    RecoveryRequired {
        /// The requested destination.
        destination: PathBuf,
        /// The matching journals, including their retained backups and candidates.
        journals: Vec<PathBuf>,
        /// Why automatic recovery cannot select a tree.
        reason: &'static str,
    },
    /// Publication failed; the previous destination was restored, or none existed.
    #[error("failed to publish {destination}")]
    Publish {
        /// The requested destination.
        destination: PathBuf,
        /// The original publication failure.
        #[source]
        source: OutputPathError,
    },
    /// Publication and rollback both failed. The complete old tree remains at `backup`.
    #[error(
        "failed to publish {destination}: {primary}; restoration also failed; \
         previous tree retained at {backup}",
        primary = Self::primary_chain(.primary)
    )]
    Restore {
        /// The requested destination.
        destination: PathBuf,
        /// The complete previous tree; automatic cleanup is disabled.
        backup: PathBuf,
        /// The original publication failure, preserved separately from restoration.
        primary: OutputPathError,
        /// The original restoration failure.
        #[source]
        source: OutputPathError,
    },
    /// Rollback restored the destination, but its directory synchronization failed.
    #[error(
        "failed to publish {destination}: {primary}; restored tree durability is unconfirmed",
        primary = Self::primary_chain(.primary)
    )]
    RestoreDurability {
        /// The restored destination.
        destination: PathBuf,
        /// The original publication failure.
        primary: OutputPathError,
        /// The restoration synchronization failure.
        #[source]
        source: OutputPathError,
    },
    /// Publication is visible; synchronization or subsequent cleanup failed.
    #[error(
        "directory {destination} committed but synchronization or cleanup failed; \
         retained backup: {backup:?}"
    )]
    Committed {
        /// The destination now containing the complete candidate.
        destination: PathBuf,
        /// The retained backup, if it has not yet been removed. A cleanup failure may have
        /// partially removed its contents; a synchronization failure leaves it complete.
        backup: Option<PathBuf>,
        /// The post-commit synchronization or cleanup failure.
        #[source]
        source: OutputPathError,
    },
}

impl DirectoryPublicationError {
    fn primary_chain(primary: &OutputPathError) -> String {
        let mut message = primary.to_string();
        let mut source = primary.source();
        while let Some(error) = source {
            message.push_str(": ");
            message.push_str(&error.to_string());
            source = error.source();
        }
        message
    }
}

/// Owns a private sibling candidate and restores the previous tree if publication fails.
///
/// Callers authorize the destination, then construct this owner before inspecting or copying the
/// previous tree. Construction locks this destination and recovers its matching stale journals.
/// Keep the owner alive while building the complete tree at [`Self::staging_path`], flush
/// application buffers, then call [`Self::publish`]. Both trees must contain only regular
/// files and directories. Windows synchronization requires writable files; Unix permits read-only
/// files. A persistent sibling `.publication-lock-<digest>` excludes concurrent publication to the
/// same destination; distinct destinations have independent locks. The digest is the full lowercase
/// BLAKE3 hash of the resolved destination filename's encoded bytes, keeping names bounded even for
/// long or non-UTF-8 filenames. Lock filenames are reserved and lock files are never removed.
/// Callers must retain one destination spelling across publication and recovery: alternate case or
/// Unicode spellings on case-insensitive filesystems are not normalized while the leaf is absent.
///
/// Portable replacement uses two renames: old destination to backup, then candidate to
/// destination. Each rename is atomic, so readers see complete trees, but there is a brief interval
/// with no destination. A crash in that interval leaves a discoverable sibling
/// `.publication-<digest>-<generation>/backup`, with an exact resolved `destination` record.
/// Construction restores a single matching backup when the destination is absent. Multiple backups
/// or an existing destination alongside a backup require manual recovery; neither tree is deleted.
/// Recovery filters journals by this destination's prefix, so unrelated publishers can prepare
/// private trees concurrently. Journals without backups are discarded, including interrupted
/// destination records. Journals for other destinations are left untouched. Use manifest snapshots
/// when readers require an uninterrupted commit pointer.
/// On Unix all candidate files/directories and both rename parents are synchronized. The old
/// backup is removed only after successful post-commit directory synchronization. Other platforms
/// retain the same ordinary-failure rollback without directory-durability guarantees.
///
/// Dropping an unpublished owner cleans only its staging. Failed restoration or post-commit sync
/// disables cleanup and returns the retained backup location in a typed error.
///
/// ```
/// use alphal00p_persistence::DirectoryPublication;
/// use std::error::Error;
/// use std::fs;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// let temporary = tempfile::tempdir()?;
/// let destination = temporary.path().join("site");
/// let publication = DirectoryPublication::new(&destination)?;
/// assert!(publication.destination().is_absolute());
/// fs::create_dir(publication.staging_path())?;
/// fs::write(publication.staging_path().join("index.html"), b"complete site")?;
/// publication.publish()?;
/// assert_eq!(fs::read(destination.join("index.html"))?, b"complete site");
/// # Ok(()) }
/// ```
#[derive(Debug)]
pub struct DirectoryPublication {
    destination: OutputPath,
    parent: OutputPath,
    transaction: OutputPath,
    staged: OutputPath,
    backup: OutputPath,
    preserve: bool,
    _lock: File,
    #[cfg(test)]
    fault: Option<PublicationFault>,
}

impl DirectoryPublication {
    /// Lock this destination, recover it, and create staging. See [`Self`].
    ///
    /// The returned staging path does not yet exist, so a serializer can create it exclusively.
    /// Missing destination ancestors are created. A stale backup may be restored before this
    /// returns, so inspect and copy [`Self::destination`] only after construction. On Unix, durable
    /// initialization synchronizes newly created directories, the transaction, and the
    /// destination's parent; it does not synchronize preexisting higher ancestors.
    ///
    /// # Errors
    /// Returns [`DirectoryPublicationError::Busy`] immediately on contention,
    /// [`DirectoryPublicationError::RecoveryRequired`] for ambiguous retained backups. Rejects
    /// symlinks within managed trees, non-directory destinations, filesystem-root destinations,
    /// clock failures, and filesystem errors. Recovery can restore the previous tree before a
    /// later synchronization or staging error; it never removes the only previous tree.
    pub fn new(destination: &Path) -> Result<Self, DirectoryPublicationError> {
        let destination = OutputPath::new(destination)?;
        let name = destination
            .as_path()
            .file_name()
            .ok_or_else(|| OutputPathError::InvalidRoot(destination.as_path().to_owned()))?;
        if name.to_str().is_some_and(Self::is_lock_filename) {
            return Err(DirectoryPublicationError::ReservedDestination(
                destination.as_path().to_owned(),
            ));
        }
        let parent = destination
            .as_path()
            .parent()
            .ok_or_else(|| OutputPathError::InvalidRoot(destination.as_path().to_owned()))?;
        let parent = OutputPath::new(parent)?;
        parent.ensure_directory_synced()?;
        let identity = OutputPath::filename_hash(destination.as_path())?;
        let lock_path = parent.join(&format!("{PUBLICATION_LOCK_PREFIX}{identity}"))?;
        lock_path.check_regular_or_missing()?;
        let lock = File::options()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(lock_path.as_path())
            .map_err(|source| OutputPathError::io(lock_path.as_path(), source))?;
        match lock.try_lock() {
            Ok(()) => {}
            Err(TryLockError::WouldBlock) => {
                return Err(DirectoryPublicationError::Busy {
                    path: lock_path.as_path().to_owned(),
                });
            }
            Err(TryLockError::Error(source)) => {
                return Err(OutputPathError::io(lock_path.as_path(), source).into());
            }
        }
        lock_path.check_regular_or_missing()?;
        if destination.metadata()?.is_some() {
            destination.check_directory()?;
        }
        destination.check_tree()?;
        let journal_prefix = format!("{PUBLICATION_PREFIX}{identity}-");
        Self::recover(&parent, &destination, &journal_prefix)?;
        for attempt in 0..CREATE_ATTEMPTS {
            let id = GenerationId::generate()?;
            let transaction = parent.join(&format!("{journal_prefix}{id}"))?;
            match transaction.create_directory() {
                Ok(()) => {
                    let publication = Self {
                        staged: transaction.join("candidate")?,
                        backup: transaction.join("backup")?,
                        destination,
                        parent,
                        transaction,
                        preserve: false,
                        _lock: lock,
                        #[cfg(test)]
                        fault: None,
                    };
                    // Recovery compares exact bytes; it never follows a path read from this record.
                    publication.transaction.join("destination")?.write(
                        publication
                            .destination
                            .as_path()
                            .as_os_str()
                            .as_encoded_bytes(),
                    )?;
                    publication.transaction.sync_tree()?;
                    publication.parent.sync_directory()?;
                    return Ok(publication);
                }
                Err(OutputPathError::Io { source, .. })
                    if source.kind() == ErrorKind::AlreadyExists
                        && attempt + 1 < CREATE_ATTEMPTS => {}
                Err(source) => return Err(source.into()),
            }
        }
        unreachable!("the last staging creation attempt returns its error")
    }

    fn is_lock_filename(name: &str) -> bool {
        name.strip_prefix(PUBLICATION_LOCK_PREFIX)
            .is_some_and(|digest| {
                digest.len() == 64
                    && digest
                        .bytes()
                        .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
            })
    }

    fn recover(
        parent: &OutputPath,
        destination: &OutputPath,
        journal_prefix: &str,
    ) -> Result<(), DirectoryPublicationError> {
        // The caller holds this destination's lock throughout validation and recovery. Filter by
        // its digest before accessing entries that other publishers may mutate or keep private.
        let mut journals = Vec::new();
        let mut backups = Vec::new();
        for entry in fs::read_dir(parent.as_path())
            .map_err(|source| OutputPathError::io(parent.as_path(), source))?
        {
            let entry = entry.map_err(|source| OutputPathError::io(parent.as_path(), source))?;
            let name = entry.file_name();
            let Some(name) = name.to_str() else {
                continue;
            };
            if !name
                .strip_prefix(journal_prefix)
                .is_some_and(|id| id.parse::<GenerationId>().is_ok())
            {
                continue;
            }
            let journal = parent.join(name)?;
            if journal.as_path() == destination.as_path() {
                continue;
            }
            if !journal
                .metadata()?
                .is_some_and(|metadata| metadata.is_dir())
            {
                continue;
            }
            let backup = journal.join("backup")?;
            if backup.metadata()?.is_some() {
                // Backups contain prior published data, so validate their complete journals
                // before restoring or deleting anything. Backup-less candidates are disposable.
                journal.check_tree()?;
                backup.check_directory()?;
                let record = journal.join("destination")?;
                record.check_regular_or_missing()?;
                if record.metadata()?.is_none() {
                    return Err(DirectoryPublicationError::RecoveryRequired {
                        destination: destination.as_path().to_owned(),
                        journals: vec![journal.as_path().to_owned()],
                        reason: "a retained backup has no destination record",
                    });
                }
                let recorded = fs::read(record.as_path())
                    .map_err(|source| OutputPathError::io(record.as_path(), source))?;
                if recorded != destination.as_path().as_os_str().as_encoded_bytes() {
                    return Err(DirectoryPublicationError::RecoveryRequired {
                        destination: destination.as_path().to_owned(),
                        journals: vec![journal.as_path().to_owned()],
                        reason: "a matching journal has an inconsistent destination record",
                    });
                }
                backups.push(backup);
            }
            // Without a backup this journal cannot contain the previous tree, even if a crash
            // interrupted writing its destination record. Reclaim it with the other candidates.
            journals.push(journal);
        }
        let reason = if backups.len() > 1 {
            Some("multiple previous trees exist and their publication order is unknown")
        } else if !backups.is_empty() && destination.metadata()?.is_some() {
            Some("both a destination and a previous tree exist; neither can be discarded safely")
        } else {
            None
        };
        if let Some(reason) = reason {
            return Err(DirectoryPublicationError::RecoveryRequired {
                destination: destination.as_path().to_owned(),
                journals: journals
                    .into_iter()
                    .map(|journal| journal.as_path().to_owned())
                    .collect(),
                reason,
            });
        }
        if let Some(backup) = backups.first() {
            Self::rename(backup, destination)?;
        }
        if !journals.is_empty() {
            // A restored tree must be durable before removing its recovery evidence. Repeating
            // these syncs also completes an earlier recovery interrupted after its rename.
            for journal in &journals {
                journal.sync_directory()?;
            }
            parent.sync_directory()?;
            for journal in journals {
                // Keep the destination record until the potentially large abandoned candidate
                // has been removed, so interrupted cleanup retains its destination context.
                Self::remove_candidate(&journal)?;
                journal.sync_directory()?;
                fs::remove_dir_all(journal.as_path())
                    .map_err(|source| OutputPathError::io(journal.as_path(), source))?;
            }
            parent.sync_directory()?;
        }
        Ok(())
    }

    fn remove_candidate(journal: &OutputPath) -> Result<(), OutputPathError> {
        journal.check_directory()?;
        let candidate = journal.as_path().join("candidate");
        // This exclusively owned candidate cannot be a backup. Unlink symlinks and special
        // entries without following them, just as when an unpublished transaction is dropped.
        match fs::symlink_metadata(&candidate) {
            Ok(metadata) if metadata.is_dir() || metadata.is_symlink() => {
                fs::remove_dir_all(&candidate)
            }
            Ok(_) => fs::remove_file(&candidate),
            Err(error) if error.kind() == ErrorKind::NotFound => Ok(()),
            Err(error) => Err(error),
        }
        .map_err(|source| OutputPathError::io(&candidate, source))
    }

    /// Borrow the absent candidate path to materialize before publishing. See [`Self`].
    #[must_use]
    pub fn staging_path(&self) -> &Path {
        self.staged.as_path()
    }

    /// Borrow the recovered, resolved destination while holding its publication lock. See [`Self`].
    #[must_use]
    pub fn destination(&self) -> &Path {
        self.destination.as_path()
    }

    /// Synchronize and publish the candidate, restoring the previous tree on rename failure.
    ///
    /// See [`Self`] for the two-rename visibility and crash-recovery contract.
    ///
    /// # Errors
    /// Returns validation errors before moving the old tree. [`DirectoryPublicationError::Publish`]
    /// means the old tree was restored (or absent initially).
    /// [`DirectoryPublicationError::Restore`] retains both failures and the complete backup.
    /// A committed error means the complete new destination is already visible; it must not be
    /// treated as a rollback.
    pub fn publish(mut self) -> Result<(), DirectoryPublicationError> {
        self.staged.check_directory()?;
        self.staged.sync_tree()?;
        self.transaction.sync_directory()?;
        let had_destination = self.destination.metadata()?.is_some();
        if had_destination {
            self.destination.check_directory()?;
            self.destination.check_tree()?;
            Self::rename(&self.destination, &self.backup)?;
        }
        // From this point onward no destructor may remove the only copy of the old tree.
        self.preserve = had_destination;
        let preparation = self.sync_parents();
        let publication = preparation.and_then(|()| {
            #[cfg(test)]
            self.injected_error(PublicationFault::Publish)?;
            Self::rename(&self.staged, &self.destination)
        });
        if let Err(primary) = publication {
            if had_destination {
                let restoration = self.restore();
                if let Err(source) = restoration {
                    return Err(DirectoryPublicationError::Restore {
                        destination: self.destination().to_owned(),
                        backup: self.backup.as_path().to_owned(),
                        primary,
                        source,
                    });
                }
                if let Err(source) = self.sync_parents() {
                    return Err(DirectoryPublicationError::RestoreDurability {
                        destination: self.destination().to_owned(),
                        primary,
                        source,
                    });
                }
            }
            self.preserve = false;
            return Err(DirectoryPublicationError::Publish {
                destination: self.destination().to_owned(),
                source: primary,
            });
        }
        self.preserve = true;
        let completion = self.sync_parents();
        #[cfg(test)]
        let completion = completion.and_then(|()| self.injected_error(PublicationFault::Committed));
        if let Err(source) = completion {
            let backup = if had_destination {
                BackupState::Retained
            } else {
                BackupState::Absent
            };
            return Err(self.committed_error(backup, source));
        }
        if had_destination && let Err(source) = self.backup.remove_tree() {
            return Err(self.committed_error(BackupState::Retained, source));
        }
        if let Err(source) = self
            .transaction
            .remove_tree()
            .and_then(|()| self.parent.sync_directory())
        {
            return Err(self.committed_error(BackupState::Absent, source));
        }
        self.preserve = false;
        Ok(())
    }

    fn rename(source: &OutputPath, destination: &OutputPath) -> Result<(), OutputPathError> {
        source.check()?;
        destination.check()?;
        fs::rename(source.as_path(), destination.as_path())
            .map_err(|source| OutputPathError::io(destination.as_path(), source))
    }

    fn sync_parents(&self) -> Result<(), OutputPathError> {
        self.transaction.sync_directory()?;
        self.parent.sync_directory()
    }

    fn restore(&self) -> Result<(), OutputPathError> {
        #[cfg(test)]
        self.injected_error(PublicationFault::Restore)?;
        Self::rename(&self.backup, &self.destination)
    }

    fn committed_error(
        &self,
        backup: BackupState,
        source: OutputPathError,
    ) -> DirectoryPublicationError {
        DirectoryPublicationError::Committed {
            destination: self.destination().to_owned(),
            backup: match backup {
                BackupState::Absent => None,
                BackupState::Retained => Some(self.backup.as_path().to_owned()),
            },
            source,
        }
    }

    #[cfg(test)]
    fn injected_error(&self, step: PublicationFault) -> Result<(), OutputPathError> {
        let fault = self.fault == Some(step)
            || (self.fault == Some(PublicationFault::Restore) && step == PublicationFault::Publish);
        if fault {
            return Err(OutputPathError::io(
                self.destination(),
                IoError::other(format!("injected {step:?} failure")),
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BackupState {
    Absent,
    Retained,
}

impl Drop for DirectoryPublication {
    fn drop(&mut self) {
        if !self.preserve && self.transaction.check_directory().is_ok() {
            // Keep the record if candidate removal fails or is interrupted.
            if Self::remove_candidate(&self.transaction).is_ok() {
                let _ = fs::remove_dir_all(self.transaction.as_path());
            }
        }
    }
}

#[cfg(test)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PublicationFault {
    Publish,
    Restore,
    Committed,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::TestDirectory;
    #[cfg(unix)]
    use std::ffi::OsString;
    #[cfg(unix)]
    use std::os::unix::ffi::OsStringExt;
    #[cfg(unix)]
    use std::os::unix::fs::{self as unix_fs, PermissionsExt};
    #[cfg(unix)]
    use std::os::unix::net::UnixListener;

    const SNAPSHOT: &str = "products/linnet/snapshots/v1/index.html";

    fn write_tree(path: &Path, index: &[u8], snapshot: &[u8]) {
        fs::create_dir_all(path.join("products/linnet/snapshots/v1")).unwrap();
        fs::write(path.join("index.html"), index).unwrap();
        fs::write(path.join(SNAPSHOT), snapshot).unwrap();
    }

    fn assert_tree(path: &Path, index: &[u8], snapshot: &[u8]) {
        assert_eq!(fs::read(path.join("index.html")).unwrap(), index);
        assert_eq!(fs::read(path.join(SNAPSHOT)).unwrap(), snapshot);
    }

    fn lock_path(destination: &Path) -> PathBuf {
        let name = destination.file_name().unwrap();
        let digest = blake3::hash(name.as_encoded_bytes()).to_hex();
        destination
            .parent()
            .unwrap()
            .join(format!("{PUBLICATION_LOCK_PREFIX}{digest}"))
    }

    fn stale_journal(parent: &Path, destination: &Path) -> PathBuf {
        let id = GenerationId::generate().unwrap();
        let name = destination.file_name().unwrap();
        let digest = blake3::hash(name.as_encoded_bytes()).to_hex();
        let journal = parent.join(format!("{PUBLICATION_PREFIX}{digest}-{id}"));
        fs::create_dir(&journal).unwrap();
        fs::write(
            journal.join("destination"),
            destination.as_os_str().as_encoded_bytes(),
        )
        .unwrap();
        journal
    }

    #[test]
    fn publication_lock_covers_preparation_and_is_released_on_drop() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        write_tree(&destination, b"old index", b"immutable old snapshot");
        let publication = DirectoryPublication::new(&destination).unwrap();
        let transaction = publication.transaction.as_path().to_owned();
        write_tree(publication.staging_path(), b"new index", b"new snapshot");
        let error = DirectoryPublication::new(&destination).unwrap_err();
        let DirectoryPublicationError::Busy { path } = error else {
            panic!("expected publication lock contention: {error}");
        };
        assert_eq!(path, lock_path(&destination));
        assert_tree(&destination, b"old index", b"immutable old snapshot");
        assert_tree(publication.staging_path(), b"new index", b"new snapshot");
        drop(publication);
        assert!(!transaction.exists());
        let next = DirectoryPublication::new(&destination).unwrap();
        assert_tree(next.destination(), b"old index", b"immutable old snapshot");
    }

    #[test]
    fn distinct_destinations_publish_concurrently_with_independent_persistent_locks() {
        let temporary = TestDirectory::new();
        let first = temporary.0.join("alice-site");
        let second = temporary.0.join("bob-site");
        let alice = DirectoryPublication::new(&first).unwrap();
        let bob = DirectoryPublication::new(&second).unwrap();
        assert_ne!(lock_path(&first), lock_path(&second));
        for destination in [&first, &second] {
            assert!(lock_path(destination).is_file());
            let error = DirectoryPublication::new(destination).unwrap_err();
            let DirectoryPublicationError::Busy { path } = error else {
                panic!("expected this destination's lock to remain held: {error}");
            };
            assert_eq!(path, lock_path(destination));
        }
        write_tree(alice.staging_path(), b"Alice index", b"Alice snapshot");
        write_tree(bob.staging_path(), b"Bob index", b"Bob snapshot");
        alice.publish().unwrap();
        assert!(matches!(
            DirectoryPublication::new(&second),
            Err(DirectoryPublicationError::Busy { .. })
        ));
        bob.publish().unwrap();
        assert_tree(&first, b"Alice index", b"Alice snapshot");
        assert_tree(&second, b"Bob index", b"Bob snapshot");
        assert_eq!(fs::read_dir(&temporary.0).unwrap().count(), 4);
    }

    #[cfg(unix)]
    #[test]
    fn private_foreign_lock_and_scoped_journal_do_not_block_a_sibling_destination() {
        let temporary = TestDirectory::new();
        let first = temporary.0.join("alice-site");
        let second = temporary.0.join("bob-site");
        let alice = DirectoryPublication::new(&first).unwrap();
        let journal = alice.transaction.as_path();
        let lock = lock_path(&first);
        let lock_permissions = fs::metadata(&lock).unwrap().permissions();
        let journal_permissions = fs::metadata(journal).unwrap().permissions();
        // Model artifacts inaccessible to another user without requiring a privileged uid switch.
        fs::set_permissions(&lock, PermissionsExt::from_mode(0o000)).unwrap();
        fs::set_permissions(journal, PermissionsExt::from_mode(0o000)).unwrap();
        let result = DirectoryPublication::new(&second);
        // Restore before any assertion so failing tests still clean their private directories.
        fs::set_permissions(&lock, lock_permissions).unwrap();
        fs::set_permissions(journal, journal_permissions).unwrap();
        let bob = result.unwrap();
        assert!(alice.transaction.as_path().is_dir());
        assert!(bob.transaction.as_path().is_dir());
        assert_ne!(alice.transaction, bob.transaction);
        // A foreign invalid record must also be skipped even when permissions are bypassed.
        drop(bob);
        fs::remove_file(journal.join("destination")).unwrap();
        unix_fs::symlink(temporary.0.join("absent"), journal.join("destination")).unwrap();
        let bob = DirectoryPublication::new(&second).unwrap();
        assert!(bob.transaction.as_path().is_dir());
        fs::remove_file(journal.join("destination")).unwrap();
    }

    #[cfg(unix)]
    #[test]
    fn lock_identity_handles_long_and_non_utf8_destination_names() {
        let temporary = TestDirectory::new();
        for name in [
            OsString::from("x".repeat(255)),
            OsString::from_vec(vec![b'x', 0xff]),
        ] {
            let destination = temporary.0.join(name);
            let publication = DirectoryPublication::new(&destination).unwrap();
            let lock = lock_path(&destination);
            assert_eq!(lock.file_name().unwrap().as_encoded_bytes().len(), 82);
            assert!(lock.is_file());
            assert!(publication.transaction.as_path().file_name().unwrap().len() < 255);
            let error = DirectoryPublication::new(&destination).unwrap_err();
            let DirectoryPublicationError::Busy { path } = error else {
                panic!("expected stable destination identity: {error}");
            };
            assert_eq!(path, lock);
            drop(publication);
            let again = DirectoryPublication::new(&destination).unwrap();
            assert!(lock.is_file());
            drop(again);
        }
    }

    #[test]
    fn matching_journal_with_a_missing_or_inconsistent_record_preserves_its_backup() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        let journal = stale_journal(&temporary.0, &destination);
        write_tree(
            &journal.join("backup"),
            b"previous index",
            b"previous snapshot",
        );
        for contents in [None, Some(b"wrong destination".as_slice())] {
            fs::remove_file(journal.join("destination")).unwrap();
            if let Some(contents) = contents {
                fs::write(journal.join("destination"), contents).unwrap();
            }
            assert!(matches!(
                DirectoryPublication::new(&destination),
                Err(DirectoryPublicationError::RecoveryRequired { .. })
            ));
            assert!(!destination.exists());
            assert_tree(
                &journal.join("backup"),
                b"previous index",
                b"previous snapshot",
            );
            if contents.is_none() {
                fs::write(journal.join("destination"), b"reset fixture record").unwrap();
            }
        }
    }

    #[test]
    fn recovery_discards_incomplete_records_without_backups() {
        for destination_exists in [false, true] {
            for contents in [None, Some(b"truncated".as_slice())] {
                let temporary = TestDirectory::new();
                let destination = temporary.0.join("site");
                if destination_exists {
                    write_tree(&destination, b"previous index", b"previous snapshot");
                }
                let journal = stale_journal(&temporary.0, &destination);
                write_tree(
                    &journal.join("candidate"),
                    b"partial index",
                    b"partial snapshot",
                );
                fs::remove_file(journal.join("destination")).unwrap();
                if let Some(contents) = contents {
                    fs::write(journal.join("destination"), contents).unwrap();
                }

                let publication = DirectoryPublication::new(&destination).unwrap();

                assert!(!journal.exists());
                assert!(!publication.staging_path().exists());
                if destination_exists {
                    assert_tree(
                        publication.destination(),
                        b"previous index",
                        b"previous snapshot",
                    );
                } else {
                    assert!(!publication.destination().exists());
                }
            }
        }
    }

    #[cfg(unix)]
    #[test]
    fn recovery_and_drop_unlink_disposable_symlinks_and_special_entries() {
        for recover in [false, true] {
            for candidate_is_link in [false, true] {
                let temporary = TestDirectory::new();
                let destination = temporary.0.join("site");
                write_tree(&destination, b"old index", b"old snapshot");
                let external = temporary.0.join("external");
                fs::create_dir(&external).unwrap();
                fs::write(external.join("sentinel"), b"outside").unwrap();
                let mut abandoned = DirectoryPublication::new(&destination).unwrap();
                let journal = abandoned.transaction.as_path().to_owned();
                if candidate_is_link {
                    unix_fs::symlink(&external, abandoned.staging_path()).unwrap();
                } else {
                    fs::create_dir(abandoned.staging_path()).unwrap();
                    unix_fs::symlink(&external, abandoned.staging_path().join("link")).unwrap();
                    unix_fs::symlink(
                        external.join("sentinel"),
                        abandoned.staging_path().join("file-link"),
                    )
                    .unwrap();
                    // Bind at a short path: the journal prefix exceeds Unix socket path limits.
                    let socket_path = temporary.0.join("socket");
                    let socket = UnixListener::bind(&socket_path).unwrap();
                    fs::rename(&socket_path, abandoned.staging_path().join("socket")).unwrap();
                    drop(socket);
                }
                fs::remove_file(journal.join("destination")).unwrap();
                unix_fs::symlink(external.join("sentinel"), journal.join("destination")).unwrap();
                abandoned.preserve = recover;
                drop(abandoned);

                if recover {
                    assert!(journal.is_dir());
                    let publication = DirectoryPublication::new(&destination).unwrap();
                    assert!(!journal.exists());
                    drop(publication);
                }
                assert!(!journal.exists());
                assert_tree(&destination, b"old index", b"old snapshot");
                assert_eq!(fs::read(external.join("sentinel")).unwrap(), b"outside");
            }
        }
    }

    #[test]
    fn interrupted_between_renames_restores_before_caller_copies_old_output() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        write_tree(&destination, b"old index", b"immutable old snapshot");
        let mut interrupted = DirectoryPublication::new(&destination).unwrap();
        write_tree(
            interrupted.staging_path(),
            b"unpublished index",
            b"unpublished snapshot",
        );
        DirectoryPublication::rename(&interrupted.destination, &interrupted.backup).unwrap();
        interrupted.sync_parents().unwrap();
        interrupted.preserve = true;
        let journal = interrupted.transaction.as_path().to_owned();
        drop(interrupted);
        assert!(!destination.exists());
        assert_tree(
            &journal.join("backup"),
            b"old index",
            b"immutable old snapshot",
        );

        let recovered = DirectoryPublication::new(&destination).unwrap();
        assert_tree(
            recovered.destination(),
            b"old index",
            b"immutable old snapshot",
        );
        assert!(!journal.exists());
        let snapshot = fs::read(recovered.destination().join(SNAPSHOT)).unwrap();
        write_tree(recovered.staging_path(), b"next index", &snapshot);
        recovered.publish().unwrap();
        assert_tree(&destination, b"next index", b"immutable old snapshot");
    }

    #[test]
    fn recovery_removes_matching_abandoned_staging_and_preserves_other_destinations() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        write_tree(&destination, b"old index", b"immutable old snapshot");
        let abandoned = stale_journal(&temporary.0, &destination);
        write_tree(
            &abandoned.join("candidate"),
            b"abandoned index",
            b"abandoned snapshot",
        );
        let foreign = stale_journal(&temporary.0, &temporary.0.join("other"));
        write_tree(
            &foreign.join("backup"),
            b"foreign index",
            b"foreign snapshot",
        );
        let unrecognized = temporary.0.join(".publication-user-data");
        fs::create_dir(&unrecognized).unwrap();
        fs::write(
            unrecognized.join("destination"),
            destination.as_os_str().as_encoded_bytes(),
        )
        .unwrap();
        write_tree(
            &unrecognized.join("backup"),
            b"user index",
            b"user snapshot",
        );
        let publication = DirectoryPublication::new(&destination).unwrap();
        assert!(!abandoned.exists());
        assert_tree(
            publication.destination(),
            b"old index",
            b"immutable old snapshot",
        );
        assert_tree(
            &foreign.join("backup"),
            b"foreign index",
            b"foreign snapshot",
        );
        assert_tree(
            &unrecognized.join("backup"),
            b"user index",
            b"user snapshot",
        );
        assert_eq!(
            fs::read(foreign.join("destination")).unwrap(),
            temporary.0.join("other").as_os_str().as_encoded_bytes()
        );
    }

    #[test]
    fn recovery_refuses_multiple_backups_before_mutating_any_journal() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        let first = stale_journal(&temporary.0, &destination);
        let second = stale_journal(&temporary.0, &destination);
        let abandoned = stale_journal(&temporary.0, &destination);
        write_tree(&first.join("backup"), b"first index", b"first snapshot");
        write_tree(&second.join("backup"), b"second index", b"second snapshot");
        write_tree(
            &abandoned.join("candidate"),
            b"candidate index",
            b"candidate snapshot",
        );
        let error = DirectoryPublication::new(&destination).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("inspect and resolve retained publication journals")
        );
        let DirectoryPublicationError::RecoveryRequired {
            journals, reason, ..
        } = error
        else {
            panic!("expected ambiguous recovery error: {error}");
        };
        assert!(reason.contains("multiple previous trees"));
        assert_eq!(journals.len(), 3);
        assert!(journals.contains(&first));
        assert!(journals.contains(&second));
        assert!(journals.contains(&abandoned));
        assert!(!destination.exists());
        assert_tree(&first.join("backup"), b"first index", b"first snapshot");
        assert_tree(&second.join("backup"), b"second index", b"second snapshot");
        assert_tree(
            &abandoned.join("candidate"),
            b"candidate index",
            b"candidate snapshot",
        );
    }

    #[test]
    fn recovery_refuses_a_regular_file_backup_without_discarding_it() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        let journal = stale_journal(&temporary.0, &destination);
        let backup = journal.join("backup");
        fs::write(&backup, b"unexpected backup contents").unwrap();
        let error = DirectoryPublication::new(&destination).unwrap_err();
        let DirectoryPublicationError::Path(OutputPathError::NotDirectory(path)) = error else {
            panic!("expected an invalid backup directory error: {error}");
        };
        assert_eq!(path, backup);
        assert_eq!(fs::read(&backup).unwrap(), b"unexpected backup contents");
        assert!(!destination.exists());
    }

    #[test]
    fn destination_with_a_journal_name_is_never_recovered_as_staging() {
        let temporary = TestDirectory::new();
        let id = GenerationId::generate().unwrap();
        for name in [
            ".publication-archive".to_owned(),
            format!("{PUBLICATION_PREFIX}{id}"),
        ] {
            let destination = temporary.0.join(name);
            write_tree(&destination, b"old index", b"immutable old snapshot");
            fs::write(
                destination.join("destination"),
                destination.as_os_str().as_encoded_bytes(),
            )
            .unwrap();
            let publication = DirectoryPublication::new(&destination).unwrap();
            assert_tree(
                publication.destination(),
                b"old index",
                b"immutable old snapshot",
            );
            assert_eq!(
                fs::read(destination.join("destination")).unwrap(),
                destination.as_os_str().as_encoded_bytes()
            );
        }
    }

    #[test]
    fn reserved_lock_destination_is_rejected_without_creating_it() {
        let temporary = TestDirectory::new();
        let destination = lock_path(&temporary.0.join("site"));
        let error = DirectoryPublication::new(&destination).unwrap_err();
        let DirectoryPublicationError::ReservedDestination(path) = error else {
            panic!("expected the reserved destination error: {error}");
        };
        assert_eq!(path, destination);
        assert!(!destination.exists());
        assert_eq!(fs::read_dir(&temporary.0).unwrap().count(), 0);
    }

    #[test]
    fn abandoned_first_publication_is_never_promoted_by_recovery() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        let partial = stale_journal(&temporary.0, &destination);
        let empty = stale_journal(&temporary.0, &destination);
        fs::create_dir(partial.join("candidate")).unwrap();
        fs::write(partial.join("candidate/index.html"), b"partial candidate").unwrap();
        let publication = DirectoryPublication::new(&destination).unwrap();
        assert!(!publication.destination().exists());
        assert!(!publication.staging_path().exists());
        assert!(!partial.exists());
        assert!(!empty.exists());
    }

    #[cfg(unix)]
    #[test]
    fn contention_is_reported_before_inspecting_the_previous_tree() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        write_tree(&destination, b"old index", b"immutable old snapshot");
        let publication = DirectoryPublication::new(&destination).unwrap();
        let external = temporary.0.join("external");
        fs::write(&external, b"external contents").unwrap();
        let link = destination.join("linked-file");
        unix_fs::symlink(&external, &link).unwrap();
        assert!(matches!(
            DirectoryPublication::new(&destination),
            Err(DirectoryPublicationError::Busy { .. })
        ));
        drop(publication);
        let error = DirectoryPublication::new(&destination).unwrap_err();
        let DirectoryPublicationError::Path(OutputPathError::Symlink(path)) = error else {
            panic!("expected a managed-tree symlink error: {error}");
        };
        assert_eq!(path, link);
        assert_eq!(fs::read(&external).unwrap(), b"external contents");
        assert_tree(&destination, b"old index", b"immutable old snapshot");
    }

    #[cfg(unix)]
    #[test]
    fn recovery_refuses_symlinks_in_a_backup_and_preserves_external_contents() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        let journal = stale_journal(&temporary.0, &destination);
        let backup = journal.join("backup");
        write_tree(&backup, b"old index", b"immutable old snapshot");
        let external = temporary.0.join("external");
        fs::write(&external, b"external contents").unwrap();
        let link = backup.join("linked-file");
        unix_fs::symlink(&external, &link).unwrap();
        let error = DirectoryPublication::new(&destination).unwrap_err();
        let DirectoryPublicationError::Path(OutputPathError::Symlink(path)) = error else {
            panic!("expected a backup symlink error: {error}");
        };
        assert_eq!(path, link);
        assert!(!destination.exists());
        assert_tree(&backup, b"old index", b"immutable old snapshot");
        assert_eq!(fs::read(&external).unwrap(), b"external contents");
    }

    #[test]
    fn publication_replaces_complete_tree_and_cleans_only_owned_staging() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        write_tree(&destination, b"old index", b"immutable old snapshot");
        let foreign = temporary.0.join(".publication-foreign");
        fs::create_dir(&foreign).unwrap();
        let publication = DirectoryPublication::new(&destination).unwrap();
        assert!(!publication.staging_path().exists());
        write_tree(
            publication.staging_path(),
            b"new index",
            b"retained old snapshot",
        );
        assert_tree(&destination, b"old index", b"immutable old snapshot");
        publication.publish().unwrap();
        assert_tree(&destination, b"new index", b"retained old snapshot");
        assert_eq!(fs::read_dir(&temporary.0).unwrap().count(), 3);
        assert!(lock_path(&destination).is_file());
        assert!(foreign.is_dir());
    }

    #[test]
    fn publication_failure_restores_index_and_nested_immutable_snapshot() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        write_tree(&destination, b"old index", b"immutable old snapshot");
        let mut publication = DirectoryPublication::new(&destination).unwrap();
        write_tree(publication.staging_path(), b"new index", b"new snapshot");
        publication.fault = Some(PublicationFault::Publish);
        let error = publication.publish().unwrap_err();
        assert!(matches!(error, DirectoryPublicationError::Publish { .. }));
        assert!(error.source().is_some());
        assert_tree(&destination, b"old index", b"immutable old snapshot");
        assert_eq!(fs::read_dir(&temporary.0).unwrap().count(), 2);
        assert!(lock_path(&destination).is_file());
    }

    #[test]
    fn failed_restoration_preserves_both_sources_and_complete_recoverable_backup() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        write_tree(&destination, b"old index", b"immutable old snapshot");
        let mut publication = DirectoryPublication::new(&destination).unwrap();
        write_tree(publication.staging_path(), b"new index", b"new snapshot");
        publication.fault = Some(PublicationFault::Restore);
        let error = publication.publish().unwrap_err();
        assert!(error.to_string().contains("injected Publish failure"));
        assert!(
            error
                .to_string()
                .contains(&destination.display().to_string())
        );
        let DirectoryPublicationError::Restore {
            backup,
            primary,
            source,
            ..
        } = error
        else {
            panic!("expected both publication and restoration failures: {error}");
        };
        assert!(matches!(primary, OutputPathError::Io { .. }));
        assert!(matches!(source, OutputPathError::Io { .. }));
        assert!(
            primary
                .source()
                .unwrap()
                .to_string()
                .contains("injected Publish failure")
        );
        assert!(
            source
                .source()
                .unwrap()
                .to_string()
                .contains("injected Restore failure")
        );
        assert!(primary.source().unwrap().is::<IoError>());
        assert!(source.source().unwrap().is::<IoError>());
        assert!(!destination.exists());
        assert_tree(&backup, b"old index", b"immutable old snapshot");
        assert_eq!(
            fs::read(backup.parent().unwrap().join("destination")).unwrap(),
            destination.as_os_str().as_encoded_bytes()
        );
        let recovered = DirectoryPublication::new(&destination).unwrap();
        assert_tree(
            recovered.destination(),
            b"old index",
            b"immutable old snapshot",
        );
        assert!(!backup.parent().unwrap().exists());
    }

    #[test]
    fn restoration_durability_error_displays_the_primary_filesystem_cause() {
        let destination = PathBuf::from("site");
        let error = DirectoryPublicationError::RestoreDurability {
            primary: OutputPathError::io(&destination, IoError::other("publication rename failed")),
            source: OutputPathError::io(&destination, IoError::other("restored tree sync failed")),
            destination,
        };
        let message = error.to_string();
        assert!(message.contains("output I/O at site: publication rename failed"));
        assert!(message.contains("restored tree durability is unconfirmed"));
        assert_eq!(
            error.source().unwrap().source().unwrap().to_string(),
            "restored tree sync failed"
        );
    }

    #[test]
    fn postcommit_failure_preserves_new_destination_and_previous_backup() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        write_tree(&destination, b"old index", b"immutable old snapshot");
        let mut publication = DirectoryPublication::new(&destination).unwrap();
        write_tree(publication.staging_path(), b"new index", b"new snapshot");
        publication.fault = Some(PublicationFault::Committed);
        let error = publication.publish().unwrap_err();
        let DirectoryPublicationError::Committed {
            backup: Some(backup),
            source,
            ..
        } = error
        else {
            panic!("expected committed error with retained backup: {error}");
        };
        assert!(source.source().unwrap().is::<IoError>());
        assert_tree(&destination, b"new index", b"new snapshot");
        assert_tree(&backup, b"old index", b"immutable old snapshot");
        let error = DirectoryPublication::new(&destination).unwrap_err();
        let DirectoryPublicationError::RecoveryRequired {
            journals, reason, ..
        } = error
        else {
            panic!("expected a retained-backup recovery error: {error}");
        };
        assert!(reason.contains("both a destination and a previous tree"));
        assert_eq!(journals, vec![backup.parent().unwrap().to_owned()]);
        assert_tree(&destination, b"new index", b"new snapshot");
        assert_tree(&backup, b"old index", b"immutable old snapshot");
        // A cleanup failure can leave only part of the previous tree. Never assume such a
        // backup is complete or delete it merely because a destination is already visible.
        fs::remove_file(backup.join("index.html")).unwrap();
        assert!(matches!(
            DirectoryPublication::new(&destination),
            Err(DirectoryPublicationError::RecoveryRequired { .. })
        ));
        assert_tree(&destination, b"new index", b"new snapshot");
        assert_eq!(
            fs::read(backup.join(SNAPSHOT)).unwrap(),
            b"immutable old snapshot"
        );
        assert!(!backup.join("index.html").exists());
    }

    #[test]
    fn absent_and_abandoned_destinations_leave_no_partial_tree() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        let publication = DirectoryPublication::new(&destination).unwrap();
        write_tree(publication.staging_path(), b"new index", b"new snapshot");
        drop(publication);
        assert!(!destination.exists());
        assert_eq!(fs::read_dir(&temporary.0).unwrap().count(), 1);
        assert!(lock_path(&destination).is_file());
        let mut publication = DirectoryPublication::new(&destination).unwrap();
        write_tree(publication.staging_path(), b"new index", b"new snapshot");
        publication.fault = Some(PublicationFault::Publish);
        assert!(matches!(
            publication.publish(),
            Err(DirectoryPublicationError::Publish { .. })
        ));
        assert!(!destination.exists());
        assert_eq!(fs::read_dir(&temporary.0).unwrap().count(), 1);
        assert!(lock_path(&destination).is_file());
        let publication = DirectoryPublication::new(&destination).unwrap();
        write_tree(publication.staging_path(), b"new index", b"new snapshot");
        publication.publish().unwrap();
        assert_tree(&destination, b"new index", b"new snapshot");
    }

    #[cfg(unix)]
    #[test]
    fn publication_accepts_readonly_assets_and_system_ancestor_alias() {
        let temporary = TestDirectory::new();
        let parent = temporary.0.join("private");
        fs::create_dir(&parent).unwrap();
        let alias = temporary.0.join("var");
        unix_fs::symlink(&parent, &alias).unwrap();
        let publication = DirectoryPublication::new(&alias.join("site")).unwrap();
        write_tree(
            publication.staging_path(),
            b"new index",
            b"immutable snapshot",
        );
        let asset = publication.staging_path().join(SNAPSHOT);
        fs::set_permissions(&asset, PermissionsExt::from_mode(0o444)).unwrap();
        publication.publish().unwrap();
        assert_tree(&parent.join("site"), b"new index", b"immutable snapshot");
        assert_eq!(
            fs::metadata(parent.join("site").join(SNAPSHOT))
                .unwrap()
                .permissions()
                .mode()
                & 0o777,
            0o444
        );
    }

    #[cfg(unix)]
    #[test]
    fn publication_refuses_symlinks_without_touching_old_or_external_trees() {
        let temporary = TestDirectory::new();
        let destination = temporary.0.join("site");
        write_tree(&destination, b"old index", b"immutable old snapshot");
        let external = temporary.0.join("external");
        write_tree(&external, b"external index", b"external snapshot");
        let publication = DirectoryPublication::new(&destination).unwrap();
        write_tree(publication.staging_path(), b"new index", b"new snapshot");
        unix_fs::symlink(&external, publication.staging_path().join("escape")).unwrap();
        assert!(matches!(
            publication.publish(),
            Err(DirectoryPublicationError::Path(OutputPathError::Symlink(_)))
        ));
        assert_tree(&destination, b"old index", b"immutable old snapshot");
        assert_tree(&external, b"external index", b"external snapshot");
        let alias = temporary.0.join("alias");
        unix_fs::symlink(&destination, &alias).unwrap();
        let publication = DirectoryPublication::new(&alias).unwrap();
        assert_eq!(publication.destination(), destination);
        assert_tree(
            publication.destination(),
            b"old index",
            b"immutable old snapshot",
        );
    }
}
