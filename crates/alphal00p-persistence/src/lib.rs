//! Checked filesystem paths and publication transactions for application-owned data.
//!
//! [OutputPath] resolves a user-selected root (including a root symlink) once, retains that
//! boundary through joins, and rejects symlinks and special entries within the managed tree.
//! Ordinary directory creation does not synchronize storage. Construct checked paths at managed
//! roots and pass them down instead of authorizing new roots at payload leaves. Callers exclude
//! concurrent hostile directory replacement; validation is not descriptor-relative race protection.
//!
//! [SnapshotWriter] publishes immutable generations by replacing one caller-named manifest under
//! an advisory writer lock. Callers own destination authorization, payload serialization, and
//! manifest validation. Check overwrite/read-only policies before acquisition and again under the
//! lock. Finish and check writes and buffer flushes before publication; do not modify the payload
//! afterward. Publication synchronizes every regular payload file and the manifest. Generations
//! remain until explicit writer-locked pruning; callers supply the full retention set. Readers
//! pinned to older generations remain valid across ordinary publication, but explicit cleanup can
//! remove those generations. Manifest names and allowed scratch directories belong to the caller.
//! Explicit root retirement locks existing nested writers, then atomically detaches the tree.
//! Callers must prevent new nested snapshot writers from starting concurrently with retirement.
//! [RetiredSnapshotRoot] retains those locks through removal, which
//! unlinks descendant symlinks without following them and leaves replacement roots untouched.
//! [SnapshotWriter::remove_retired_roots] retries interrupted removals for one selected snapshot
//! root.
//!
//! [DirectoryPublication] stages a replacement directory beside its destination. Acquisition
//! locks its destination and recovers or refuses stale recovery journals before callers inspect
//! the prior destination. Keep the transaction alive while copying retained content into its stage.
//! Publication uses two renames with rollback and a brief destination gap; see its contract.
//!
//! The filesystem must support locking, synchronization, and atomic replacement. Checked
//! synchronization reports delayed write errors before commit. Unix initialization synchronizes
//! directories created by that call plus the managed root. Snapshot acquisition also attempts
//! parent synchronization on every call, tolerating permission and unsupported-operation errors
//! at that parent; real I/O failures remain fatal on retry. This supports search-only parents and
//! filesystem mount boundaries, with parent-entry durability then left to the caller. Preexisting
//! higher ancestors are never synchronized. If `/existing/new/root` is created, the `new` entry in
//! `/existing` is not synchronized and can be lost after power failure. Prepare durable ancestors
//! first.
//! Publication synchronizes payload files and directories from children to parents, at a cost of
//! one walk, open, and synchronization per entry. Windows files must remain writable; Unix accepts
//! read-only payload files. Other platforms provide atomic visibility without directory-durability
//! guarantees. Hardware and filesystem guarantees still apply. [SnapshotError::Committed] means
//! the new pointer is visible but durability is unconfirmed: retain both new and prior generations.
//!
//! # Example
//!
//! A caller supplies its own manifest filename and serializes its domain descriptor.
//!
//! ```
//! use alphal00p_persistence::{ManifestName, SnapshotWriter, GENERATIONS_DIRECTORY};
//! use std::error::Error;
//! use std::fs;
//!
//! # fn main() -> Result<(), Box<dyn Error>> {
//! # let temporary = tempfile::tempdir()?;
//! # let root = temporary.path().join("snapshot");
//! // Application policy has authorized this fresh destination.
//! let name = ManifestName::new("checkpoint.toml")?;
//! let mut writer = SnapshotWriter::acquire(&root, &[name])?;
//! let staged = writer.begin()?;
//! staged.output_path()?.join("payload.bin")?.write(b"complete payload")?;
//! let manifest = format!("generation = \"{}\"\n", staged.id());
//! let published = staged.publish(manifest.as_bytes(), name)?;
//! assert_eq!(fs::read(writer.root().join(name.as_str()))?, manifest.as_bytes());
//! let expected = writer.root().join(GENERATIONS_DIRECTORY).join(published.id().as_str());
//! assert_eq!(published.path(), expected);
//! assert_eq!(fs::read(published.path().join("payload.bin"))?, b"complete payload");
//! drop(writer); // Releases the lock; the persistent lock file remains.
//! # Ok(())
//! # }
//! ```

#![forbid(unsafe_code)]

mod output;
mod publication;
mod snapshot;

pub use output::{OutputPath, OutputPathError};
pub use publication::{DirectoryPublication, DirectoryPublicationError};
pub use snapshot::{
    GENERATIONS_DIRECTORY, GenerationId, ManifestName, PublishedGeneration, RetiredSnapshotRoot,
    SNAPSHOT_LOCK_FILE, SnapshotError, SnapshotWriter, StagedGeneration,
};

#[cfg(test)]
mod test_support;
