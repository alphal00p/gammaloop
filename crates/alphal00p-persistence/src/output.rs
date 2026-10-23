//! Filesystem validation shared by snapshots and directory publication.

#[cfg(test)]
use std::cell::RefCell;
use std::env;
use std::fs::{self, File, Metadata};
use std::io::{Error as IoError, ErrorKind, Write};
use std::path::{Component, Path, PathBuf};
#[cfg(test)]
use std::thread_local;

use blake3::Hash;
use thiserror::Error;

#[cfg(test)]
thread_local! {
    static DIRECTORY_SYNC_FAULT: RefCell<Option<(PathBuf, ErrorKind)>> =
        const { RefCell::new(None) };
    static DIRECTORY_SYNC_TRACE: RefCell<Option<Vec<PathBuf>>> = const { RefCell::new(None) };
}

/// A rejected path or failed filesystem operation. See the [`OutputPath`] example.
#[derive(Debug, Error)]
pub enum OutputPathError {
    /// Names must be a single nonempty filename without separators or dot components.
    #[error("invalid output component {0:?}: expected one nonempty filename")]
    InvalidComponent(String),
    /// Relative paths must contain only nonempty normal filename components.
    #[error("invalid relative output path {0:?}")]
    InvalidRelative(PathBuf),
    /// The requested publication root has no usable parent.
    #[error("output root has no usable parent: {0}")]
    InvalidRoot(PathBuf),
    /// Symlinks are forbidden below the resolved managed root.
    #[error("refusing symlink in output path {0}")]
    Symlink(PathBuf),
    /// The entry must be a regular file.
    #[error("expected a regular file at {0}")]
    NotFile(PathBuf),
    /// The entry must be a directory.
    #[error("expected a directory at {0}")]
    NotDirectory(PathBuf),
    /// Payload trees may contain only directories and regular files.
    #[error("expected a regular file or directory at {0}")]
    SpecialEntry(PathBuf),
    /// An operating-system failure, with its original source.
    #[error("output I/O at {path}")]
    Io {
        /// The affected entry.
        path: PathBuf,
        /// The original filesystem failure.
        #[source]
        source: IoError,
    },
}

impl OutputPathError {
    pub(crate) fn io(path: &Path, source: IoError) -> Self {
        Self::Io {
            path: path.to_owned(),
            source,
        }
    }
}

/// A path retaining the boundary of its user-selected managed root.
///
/// The supplied root and its ancestors are resolved once, permitting user-selected symlinks
/// and system paths such as macOS `/var` and `/tmp`. Descendant symlinks remain forbidden.
/// [`Self::join`] preserves this boundary and validates names without changing them. Constructors
/// and checks never write. Construct this at the managed root, then join descendants;
/// constructing a new root at a payload leaf would deliberately authorize its ancestors.
/// The caller must exclude concurrent hostile directory replacement; these checks are not
/// descriptor-relative protection against symlink races.
///
/// ```no_run
/// use alphal00p_persistence::OutputPath;
/// use std::error::Error;
/// # fn main() -> Result<(), Box<dyn Error>> {
/// let root = OutputPath::new("export".as_ref())?;
/// root.ensure_directory()?;
/// let child = root.join("process")?;
/// child.create_directory()?;
/// child.join("graph.dot")?.create_file()?;
/// child.join("template.typ")?.write(b"complete template")?;
/// child.check_tree()?;
/// # Ok(()) }
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct OutputPath {
    root: PathBuf,
    path: PathBuf,
}

impl OutputPath {
    pub(crate) fn filename_hash(path: &Path) -> Result<Hash, OutputPathError> {
        let name = path
            .file_name()
            .ok_or_else(|| OutputPathError::InvalidRoot(path.to_owned()))?;
        Ok(blake3::hash(name.as_encoded_bytes()))
    }

    /// Resolve the user-selected root once without creating anything. See [`Self`].
    ///
    /// # Errors
    /// Rejects dangling symlinks and I/O failures resolving the root or its ancestors. Parent
    /// components in the user-selected root follow filesystem semantics; child names stay strict.
    pub fn new(path: &Path) -> Result<Self, OutputPathError> {
        let absolute = if path.is_absolute() {
            path.to_owned()
        } else {
            env::current_dir()
                .map_err(|source| OutputPathError::io(path, source))?
                .join(path)
        };
        let root = Self::resolve_ancestor(&absolute)?;
        let output = Self {
            path: root.to_owned(),
            root,
        };
        output.check()?;
        Ok(output)
    }

    fn resolve_ancestor(path: &Path) -> Result<PathBuf, OutputPathError> {
        let mut resolved = PathBuf::new();
        for component in path.components() {
            match component {
                Component::Prefix(_) | Component::RootDir => resolved.push(component),
                Component::CurDir => {}
                Component::ParentDir => {
                    resolved.pop();
                }
                Component::Normal(name) => {
                    let candidate = resolved.join(name);
                    match fs::canonicalize(&candidate) {
                        Ok(path) => resolved = path,
                        Err(error) if error.kind() == ErrorKind::NotFound => {
                            // A dangling symlink is not a missing directory. Recheck each component
                            // even after a missing/.. pair: the next component may exist as a link.
                            match fs::symlink_metadata(&candidate) {
                                Ok(metadata) if metadata.file_type().is_symlink() => {
                                    return Err(OutputPathError::Symlink(candidate));
                                }
                                Ok(_) => return Err(OutputPathError::io(&candidate, error)),
                                Err(source) if source.kind() == ErrorKind::NotFound => {
                                    resolved = candidate
                                }
                                Err(source) => return Err(OutputPathError::io(&candidate, source)),
                            }
                        }
                        Err(source) => return Err(OutputPathError::io(&candidate, source)),
                    }
                }
            }
        }
        Ok(resolved)
    }

    /// Reject separators, absolute paths, dot components, and empty names.
    ///
    /// ```
    /// use alphal00p_persistence::{OutputPath, OutputPathError};
    /// # fn main() -> Result<(), OutputPathError> {
    /// OutputPath::validate_component("payload.bin")?;
    /// assert!(OutputPath::validate_component("../payload.bin").is_err());
    /// # Ok(()) }
    /// ```
    ///
    /// # Errors
    /// Returns [`OutputPathError::InvalidComponent`] for a name outside this contract.
    pub fn validate_component(name: &str) -> Result<(), OutputPathError> {
        if name.is_empty()
            || name.contains(['/', '\\', '\0'])
            || !matches!(
                Path::new(name).components().next(),
                Some(Component::Normal(_))
            )
            || Path::new(name).components().count() != 1
        {
            return Err(OutputPathError::InvalidComponent(name.to_owned()));
        }
        Ok(())
    }

    /// Append one unchanged component, retaining the original boundary. See [`Self`].
    ///
    /// # Errors
    /// Rejects invalid names, symlinks in the managed path, and inspection failures.
    pub fn join(&self, name: &str) -> Result<Self, OutputPathError> {
        Self::validate_component(name)?;
        let child = Self {
            root: self.root.to_owned(),
            path: self.path.join(name),
        };
        child.check()?;
        Ok(child)
    }

    /// Append a relative path without changing the managed boundary. See [`Self`].
    ///
    /// Unlike [`Self::join`], this accepts multiple components, such as `process/integrand.bin`.
    /// Every component must satisfy [`Self::validate_component`]; empty, dot, parent, and absolute
    /// paths are rejected rather than normalized. No entries are created.
    ///
    /// ```
    /// # use alphal00p_persistence::OutputPath;
    /// # use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// let root = OutputPath::new("export".as_ref())?;
    /// let payload = root.join_relative("process/integrand.bin".as_ref())?;
    /// assert!(payload.as_path().starts_with(root.as_path()));
    /// assert!(root.join_relative("../escape".as_ref()).is_err());
    /// # Ok(()) }
    /// ```
    ///
    /// # Errors
    /// Rejects invalid relative paths, non-UTF-8 names, symlinks, and inspection failures.
    pub fn join_relative(&self, relative: &Path) -> Result<Self, OutputPathError> {
        let invalid = || OutputPathError::InvalidRelative(relative.to_owned());
        let text = relative.to_str().ok_or_else(invalid)?;
        if text.is_empty() || relative.is_absolute() {
            return Err(invalid());
        }
        // Path::components normalizes internal `.` and duplicate separators; validate the raw
        // spelling first so malformed persisted members cannot acquire a different meaning.
        for name in text.split('/') {
            Self::validate_component(name).map_err(|_| invalid())?;
        }
        let child = Self {
            root: self.root.to_owned(),
            path: self.path.join(relative),
        };
        child.check()?;
        Ok(child)
    }

    /// Borrow the absolute path with resolved ancestors. See [`Self`].
    #[must_use]
    pub fn as_path(&self) -> &Path {
        &self.path
    }

    /// Recheck every existing component from the retained root, without writes. See [`Self`].
    ///
    /// # Errors
    /// Rejects symlinks and filesystem inspection failures.
    pub fn check(&self) -> Result<(), OutputPathError> {
        let relative = self
            .path
            .strip_prefix(&self.root)
            .expect("joined output paths remain within their retained root");
        let mut current = self.root.to_owned();
        Self::metadata_at(&current)?;
        for component in relative.components() {
            current.push(component);
            Self::metadata_at(&current)?;
        }
        Ok(())
    }

    fn metadata_at(path: &Path) -> Result<Option<Metadata>, OutputPathError> {
        match fs::symlink_metadata(path) {
            Ok(metadata) if metadata.file_type().is_symlink() => {
                Err(OutputPathError::Symlink(path.to_owned()))
            }
            Ok(metadata) => Ok(Some(metadata)),
            Err(error) if error.kind() == ErrorKind::NotFound => Ok(None),
            Err(source) => Err(OutputPathError::io(path, source)),
        }
    }

    pub(crate) fn metadata(&self) -> Result<Option<Metadata>, OutputPathError> {
        self.check()?;
        Self::metadata_at(&self.path)
    }

    pub(crate) fn check_directory(&self) -> Result<(), OutputPathError> {
        if !self.metadata()?.is_some_and(|metadata| metadata.is_dir()) {
            return Err(OutputPathError::NotDirectory(self.path.to_owned()));
        }
        Ok(())
    }

    pub(crate) fn check_regular_or_missing(&self) -> Result<(), OutputPathError> {
        if self.metadata()?.is_some_and(|metadata| !metadata.is_file()) {
            return Err(OutputPathError::NotFile(self.path.to_owned()));
        }
        Ok(())
    }

    /// Check the entire subtree before handing it to a serializer. Never writes. See [`Self`].
    ///
    /// # Errors
    /// Rejects symlinks, special entries, and filesystem inspection failures.
    pub fn check_tree(&self) -> Result<(), OutputPathError> {
        self.tree_entries().map(|_| ())
    }

    fn tree_entries(&self) -> Result<Vec<PathBuf>, OutputPathError> {
        self.check()?;
        let mut pending = vec![self.path.to_owned()];
        let mut entries = Vec::new();
        while let Some(path) = pending.pop() {
            let Some(metadata) = Self::metadata_at(&path)? else {
                continue;
            };
            if metadata.is_dir() {
                for entry in
                    fs::read_dir(&path).map_err(|source| OutputPathError::io(&path, source))?
                {
                    pending.push(
                        entry
                            .map_err(|source| OutputPathError::io(&path, source))?
                            .path(),
                    );
                }
            } else if !metadata.is_file() {
                return Err(OutputPathError::SpecialEntry(path));
            }
            entries.push(path);
        }
        Ok(entries)
    }

    /// Synchronize a checked regular file or subtree, visiting children before their directories.
    ///
    /// Callers must finish and flush application buffers, authorize this operation, and exclude
    /// concurrent modifications. Unix synchronization accepts read-only files and includes
    /// directories; Windows requires writable files. Other platforms synchronize files without
    /// directory-durability guarantees. A missing target is a no-op; no entries are created.
    /// This does not synchronize the target's containing directory or publish a commit pointer.
    /// See [`Self::remove_tree`] for an example using both tree and parent synchronization.
    ///
    /// # Errors
    /// Rejects symlinks and special entries before synchronization. Returns inspection, file-open,
    /// or synchronization errors with their original sources. An I/O failure can occur after
    /// earlier entries have been synchronized; there is no rollback.
    pub fn sync_tree(&self) -> Result<(), OutputPathError> {
        // Reverse preorder synchronizes children before their directory entries.
        for path in self.tree_entries()?.into_iter().rev() {
            let Some(metadata) = Self::metadata_at(&path)? else {
                return Err(OutputPathError::io(
                    &path,
                    IoError::from(ErrorKind::NotFound),
                ));
            };
            if metadata.is_dir() {
                Self::sync_directory_at(&path)?;
            } else {
                // Windows requires writable handles; Unix can synchronize read-only artifacts.
                #[cfg(windows)]
                let file = File::options().write(true).open(&path);
                #[cfg(not(windows))]
                let file = File::open(&path);
                file.and_then(|file| file.sync_all())
                    .map_err(|source| OutputPathError::io(&path, source))?;
            }
        }
        Ok(())
    }

    /// Synchronize this existing directory's entries on Unix, without traversing its children.
    ///
    /// Callers authorize the operation and exclude concurrent modifications. Use this on the
    /// containing directory after creating, renaming, or removing an entry; child contents need
    /// their own synchronization. Untouched ancestors are not opened. On other platforms this
    /// validates the directory without providing directory-durability guarantees.
    /// See [`Self::remove_tree`] for a complete example.
    ///
    /// # Errors
    /// Rejects missing or non-directory targets, symlinks in the managed path, and filesystem
    /// inspection or synchronization failures. No directories are created.
    pub fn sync_directory(&self) -> Result<(), OutputPathError> {
        self.check_directory()?;
        Self::sync_directory_at(&self.path)
    }

    fn sync_directory_at(path: &Path) -> Result<(), OutputPathError> {
        #[cfg(test)]
        DIRECTORY_SYNC_TRACE.with(|trace| {
            if let Some(paths) = trace.borrow_mut().as_mut() {
                paths.push(path.to_owned());
            }
        });
        #[cfg(test)]
        if let Some(kind) = DIRECTORY_SYNC_FAULT.with(|fault| {
            fault
                .borrow()
                .as_ref()
                .and_then(|(target, kind)| (target == path).then_some(*kind))
        }) {
            return Err(OutputPathError::io(path, IoError::from(kind)));
        }
        #[cfg(unix)]
        File::open(path)
            .and_then(|file| file.sync_all())
            .map_err(|source| OutputPathError::io(path, source))?;
        #[cfg(not(unix))]
        let _ = path;
        Ok(())
    }

    /// Create checked missing directories without synchronizing storage. See [`Self`].
    ///
    /// Ordinary exports and presentation writes need only checked creation. Snapshot and
    /// publication roots additionally synchronize their initialization boundary.
    ///
    /// # Errors
    /// Rejects symlinks, non-directory entries, and creation or inspection failures.
    pub fn ensure_directory(&self) -> Result<(), OutputPathError> {
        self.check()?;
        fs::create_dir_all(&self.path).map_err(|source| OutputPathError::io(&self.path, source))?;
        self.check_directory()
    }

    /// Create directories and synchronize only the initialization boundary on Unix.
    ///
    /// Synchronizes this directory and directories created by this call. Its immediate parent
    /// is synchronized on every acquisition, so retrying also completes a failed parent sync.
    /// Parent synchronization tolerates permission and unsupported-operation errors: a managed
    /// root may live beneath a search-only parent or a filesystem mount boundary. That parent's
    /// entry durability then depends on the caller. Preexisting higher ancestors are never
    /// opened for synchronization. Ordinary writes use [`Self::ensure_directory`].
    /// Creating `/existing/new/root` synchronizes `root` and `new`, but not `/existing`.
    /// The new `new` entry in `/existing` is therefore not guaranteed durable after power loss;
    /// callers needing that guarantee must prepare durable ancestors before initialization.
    ///
    /// # Errors
    /// Rejects symlinks, non-directory entries, and creation or synchronization failures.
    /// Created directories can remain visible on error. Retrying synchronizes the target and
    /// immediate parent, but cannot guarantee durability of previously created higher ancestors.
    pub(crate) fn ensure_durable_directory(&self) -> Result<(), OutputPathError> {
        let parent_synced = self.ensure_directory_synced()?;
        if !parent_synced && let Some(parent) = self.path.parent() {
            Self::sync_initialization_parent(parent)?;
        }
        Ok(())
    }

    // Directory publication initializes its destination's parent before creating a sibling stage.
    // Synchronize that parent and newly created ancestors, without opening its existing parent.
    // Returns whether this call already synchronized the target's immediate parent.
    pub(crate) fn ensure_directory_synced(&self) -> Result<bool, OutputPathError> {
        self.check()?;
        let mut missing = Vec::new();
        for path in self.path.ancestors() {
            match Self::metadata_at(path)? {
                Some(metadata) if metadata.is_dir() => break,
                Some(_) => return Err(OutputPathError::NotDirectory(path.to_owned())),
                None => missing.push(path),
            }
        }
        let mut created = Vec::new();
        for path in missing.into_iter().rev() {
            match fs::create_dir(path) {
                Ok(()) => created.push(path),
                Err(source) if source.kind() == ErrorKind::AlreadyExists => {
                    if !Self::metadata_at(path)?.is_some_and(|metadata| metadata.is_dir()) {
                        return Err(OutputPathError::NotDirectory(path.to_owned()));
                    }
                }
                Err(source) => return Err(OutputPathError::io(path, source)),
            }
        }
        self.check_directory()?;
        Self::sync_directory_at(&self.path)?;
        let parent = self.path.parent();
        let parent_synced = parent.is_some_and(|parent| created.contains(&parent));
        for path in created.into_iter().rev() {
            if Some(path) == parent {
                Self::sync_initialization_parent(path)?;
            } else if path != self.path {
                Self::sync_directory_at(path)?;
            }
        }
        Ok(parent_synced)
    }

    fn sync_initialization_parent(path: &Path) -> Result<(), OutputPathError> {
        match Self::sync_directory_at(path) {
            // EACCES/EPERM, EINVAL and ENOTSUP describe parent access or filesystem capabilities,
            // not delayed payload write failures. The managed directory itself remains strict.
            Err(OutputPathError::Io { source, .. })
                if matches!(
                    source.kind(),
                    ErrorKind::PermissionDenied | ErrorKind::InvalidInput | ErrorKind::Unsupported
                ) =>
            {
                Ok(())
            }
            // Unix EBADF can also mean the filesystem does not support fsync on directories.
            #[cfg(unix)]
            Err(OutputPathError::Io { source, .. }) if source.raw_os_error() == Some(9) => Ok(()),
            result => result,
        }
    }

    /// Exclusively create a staging directory. See [`Self`].
    ///
    /// # Errors
    /// Rejects symlinks, existing entries, missing parents, and filesystem failures.
    pub fn create_directory(&self) -> Result<(), OutputPathError> {
        self.check()?;
        fs::create_dir(&self.path).map_err(|source| OutputPathError::io(&self.path, source))
    }

    /// Replace a regular file with a fresh entry, preserving existing hard links. See [`Self`].
    ///
    /// # Errors
    /// Rejects symlinks, directories, special entries, and filesystem failures.
    pub fn create_file(&self) -> Result<File, OutputPathError> {
        self.check_regular_or_missing()?;
        if self.metadata()?.is_some() {
            fs::remove_file(&self.path)
                .map_err(|source| OutputPathError::io(&self.path, source))?;
        }
        File::create_new(&self.path).map_err(|source| OutputPathError::io(&self.path, source))
    }

    /// Write and flush a checked regular file without synchronizing storage. See [`Self`].
    ///
    /// # Errors
    /// Returns path-validation, file-creation, write, or flush errors.
    pub fn write(&self, contents: &[u8]) -> Result<(), OutputPathError> {
        let mut file = self.create_file()?;
        file.write_all(contents)
            .and_then(|()| file.flush())
            .map_err(|source| OutputPathError::io(&self.path, source))
    }

    /// Remove an authorized regular file or directory subtree after checking every entry.
    ///
    /// The caller must own the entries being removed, authorize their deletion, and exclude
    /// concurrent modifications. No lock is acquired. A missing target is a no-op. This does not
    /// synchronize the containing directory: retain its checked path and synchronize it separately
    /// when removal must be durable. See [`Self`] for the retained managed-root boundary.
    ///
    /// ```no_run
    /// use alphal00p_persistence::OutputPath;
    /// use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// // Application policy authorizes this root and its exclusively owned staging directory.
    /// let root = OutputPath::new("authorized-output".as_ref())?;
    /// root.ensure_directory()?;
    /// let staging = root.join("abandoned-stage")?;
    /// staging.create_directory()?;
    /// staging.join("payload.bin")?.write(b"flushed payload")?;
    /// staging.sync_tree()?;
    /// root.sync_directory()?;
    /// staging.remove_tree()?;
    /// root.sync_directory()?;
    /// # Ok(()) }
    /// ```
    ///
    /// # Errors
    /// Rejects symlinks and special entries anywhere in the subtree before deleting anything.
    /// Returns inspection or removal errors with their original sources. An I/O failure during
    /// deletion can leave a partially removed subtree; the operation does not roll back.
    pub fn remove_tree(&self) -> Result<(), OutputPathError> {
        self.check_tree()?;
        match self.metadata()? {
            Some(metadata) if metadata.is_dir() => fs::remove_dir_all(&self.path),
            Some(_) => fs::remove_file(&self.path),
            None => return Ok(()),
        }
        .map_err(|source| OutputPathError::io(&self.path, source))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_support::TestDirectory;
    use crate::{DirectoryPublication, GenerationId, ManifestName, SnapshotError, SnapshotWriter};
    #[cfg(unix)]
    use std::os::unix::fs as unix_fs;

    #[test]
    fn reads_preserve_missing_roots_and_joins_retain_names() {
        let temporary = TestDirectory::new();
        let missing = temporary.0.join("absent/managed");
        let root = OutputPath::new(&missing).unwrap();
        let child = root
            .join("a name with spaces and λ")
            .unwrap()
            .join("result.dat")
            .unwrap();
        child.check_tree().unwrap();
        assert_eq!(
            child.as_path(),
            missing.join("a name with spaces and λ/result.dat")
        );
        assert!(!temporary.0.join("absent").exists());
        for bad in [
            "",
            ".",
            "..",
            "../escape",
            "/absolute",
            "a/b",
            "a\\b",
            "a\0b",
        ] {
            assert!(
                matches!(root.join(bad), Err(OutputPathError::InvalidComponent(_))),
                "{bad:?}"
            );
        }
        assert_eq!(
            OutputPath::new(&missing.join("..")).unwrap().as_path(),
            temporary.0.join("absent")
        );
    }

    #[test]
    fn regular_file_replacement_preserves_hard_links_and_directory_targets() {
        let temporary = TestDirectory::new();
        let root = OutputPath::new(&temporary.0).unwrap();
        let file = root.join("file").unwrap();
        file.write(b"old").unwrap();
        fs::hard_link(file.as_path(), temporary.0.join("alias")).unwrap();
        file.write(b"new").unwrap();
        assert_eq!(fs::read(temporary.0.join("alias")).unwrap(), b"old");
        let directory = root.join("directory").unwrap();
        directory.create_directory().unwrap();
        directory
            .join("sentinel")
            .unwrap()
            .write(b"preserve")
            .unwrap();
        assert!(matches!(
            directory.create_file(),
            Err(OutputPathError::NotFile(_))
        ));
        assert_eq!(
            fs::read(directory.as_path().join("sentinel")).unwrap(),
            b"preserve"
        );
    }

    #[cfg(unix)]
    #[test]
    fn system_ancestor_symlink_resolves_once_but_managed_symlinks_are_rejected() {
        let temporary = TestDirectory::new();
        let real = temporary.0.join("private");
        let alias = temporary.0.join("var");
        fs::create_dir(&real).unwrap();
        unix_fs::symlink(&real, &alias).unwrap();
        let root = OutputPath::new(&alias.join("managed")).unwrap();
        assert_eq!(root.as_path(), real.join("managed"));
        root.ensure_directory().unwrap();
        let child = root.join("child").unwrap();
        child.create_directory().unwrap();
        let leaf = child.join("payload").unwrap();
        // The original spelling ceases to exist: retained paths must not resolve it again.
        fs::remove_file(&alias).unwrap();
        leaf.write(b"payload").unwrap();
        let external = temporary.0.join("external");
        fs::rename(child.as_path(), &external).unwrap();
        unix_fs::symlink(&external, child.as_path()).unwrap();
        assert!(matches!(
            leaf.write(b"bad"),
            Err(OutputPathError::Symlink(_))
        ));
        assert!(matches!(
            root.join("child"),
            Err(OutputPathError::Symlink(_))
        ));
        assert!(matches!(
            root.check_tree(),
            Err(OutputPathError::Symlink(_))
        ));
        assert_eq!(fs::read(external.join("payload")).unwrap(), b"payload");
        fs::remove_file(child.as_path()).unwrap();
        unix_fs::symlink(temporary.0.join("does-not-exist"), child.as_path()).unwrap();
        assert!(matches!(
            root.join("child"),
            Err(OutputPathError::Symlink(_))
        ));
        for path in [child.as_path().to_owned(), child.as_path().join(".")] {
            assert!(matches!(
                OutputPath::new(&path),
                Err(OutputPathError::Symlink(_))
            ));
        }
    }
    #[test]
    fn relative_parent_root_resolves_without_writing() {
        let relative = Path::new("../").join(format!("root-{}", GenerationId::generate().unwrap()));
        let checked = OutputPath::new(&relative).unwrap();
        let cwd = env::current_dir().unwrap();
        assert_eq!(
            checked.as_path(),
            cwd.parent().unwrap().join(relative.file_name().unwrap())
        );
        assert!(!checked.as_path().exists());
        assert!(matches!(
            checked.join(".."),
            Err(OutputPathError::InvalidComponent(_))
        ));
    }

    #[cfg(unix)]
    #[test]
    fn parent_components_and_selected_roots_resolve_symlinks_once() {
        let temporary = TestDirectory::new();
        let private = temporary.0.join("private");
        fs::create_dir_all(private.join("nested")).unwrap();
        let alias = temporary.0.join("alias");
        unix_fs::symlink(private.join("nested"), &alias).unwrap();
        let path = OutputPath::new(&alias.join("../managed")).unwrap();
        assert_eq!(path.as_path(), private.join("managed"));
        assert!(!path.as_path().exists());
        unix_fs::symlink(private.join("nested"), private.join("managed")).unwrap();
        let root = OutputPath::new(&alias.join("../managed/.")).unwrap();
        assert_eq!(root.as_path(), private.join("nested"));
        assert_eq!(
            OutputPath::new(&alias.join(".")).unwrap().as_path(),
            root.as_path()
        );
        fs::remove_file(&alias).unwrap();
        root.join("payload")
            .unwrap()
            .write(b"resolved once")
            .unwrap();
        unix_fs::symlink(&private, root.as_path().join("link")).unwrap();
        assert!(matches!(
            root.join_relative("link/payload".as_ref()),
            Err(OutputPathError::Symlink(_))
        ));
    }
    #[cfg(unix)]
    #[test]
    fn cancelled_missing_components_still_resolve_subsequent_existing_symlinks() {
        let temporary = TestDirectory::new();
        let private = temporary.0.join("private");
        fs::create_dir_all(private.join("nested")).unwrap();
        let alias = temporary.0.join("alias");
        unix_fs::symlink(private.join("nested"), &alias).unwrap();
        let child = OutputPath::new(&temporary.0.join("missing/../alias/managed")).unwrap();
        assert_eq!(child.as_path(), private.join("nested/managed"));
        let parent = OutputPath::new(&temporary.0.join("missing/../alias/../state")).unwrap();
        assert_eq!(parent.as_path(), private.join("state"));
        fs::remove_file(alias).unwrap();
        child.ensure_directory().unwrap();
        parent.ensure_directory().unwrap();
        child
            .join("payload")
            .unwrap()
            .write(b"retained resolved root")
            .unwrap();
        assert_eq!(
            fs::read(private.join("nested/managed/payload")).unwrap(),
            b"retained resolved root"
        );
        assert!(!temporary.0.join("missing").exists());
        assert!(!temporary.0.join("state").exists());
    }
    struct DirectorySyncFault;

    impl DirectorySyncFault {
        fn new(path: &Path) -> Self {
            Self::with_kind(path, ErrorKind::PermissionDenied)
        }

        fn with_kind(path: &Path, kind: ErrorKind) -> Self {
            DIRECTORY_SYNC_FAULT.with(|fault| {
                assert!(fault.replace(Some((path.to_owned(), kind))).is_none());
            });
            Self
        }
    }

    impl Drop for DirectorySyncFault {
        fn drop(&mut self) {
            DIRECTORY_SYNC_FAULT.with(|fault| {
                fault.replace(None);
            });
        }
    }

    struct DirectorySyncTrace;

    impl DirectorySyncTrace {
        fn new() -> Self {
            DIRECTORY_SYNC_TRACE.with(|trace| {
                assert!(trace.replace(Some(Vec::new())).is_none());
            });
            Self
        }

        fn count(&self, path: &Path) -> usize {
            DIRECTORY_SYNC_TRACE.with(|trace| {
                trace
                    .borrow()
                    .as_ref()
                    .unwrap()
                    .iter()
                    .filter(|item| *item == path)
                    .count()
            })
        }
    }

    impl Drop for DirectorySyncTrace {
        fn drop(&mut self) {
            DIRECTORY_SYNC_TRACE.with(|trace| {
                trace.replace(None);
            });
        }
    }

    #[test]
    fn durable_initialization_requires_root_sync_for_existing_directories() {
        let temporary = TestDirectory::new();
        let root = temporary.0.join("parent/snapshot");
        fs::create_dir_all(&root).unwrap();
        fs::write(temporary.0.join("sentinel"), b"prior data").unwrap();
        let checked = OutputPath::new(&root).unwrap();
        let manifest = ManifestName::new("checkpoint.toml").unwrap();
        let _fault = DirectorySyncFault::new(&root);
        // Ordinary output creation has no synchronization prerequisite.
        checked.ensure_directory().unwrap();
        let error = checked.ensure_durable_directory().unwrap_err();
        assert!(matches!(error, OutputPathError::Io { path, source }
            if path == root && source.kind() == ErrorKind::PermissionDenied));
        let error = SnapshotWriter::acquire(&root, &[manifest]).unwrap_err();
        assert!(
            matches!(error, SnapshotError::Path(OutputPathError::Io { path, source })
            if path == root && source.kind() == ErrorKind::PermissionDenied)
        );
        assert!(!root.join(manifest.as_str()).exists());
        assert!(!root.join(".snapshot.lock").exists());
        assert_eq!(
            fs::read(temporary.0.join("sentinel")).unwrap(),
            b"prior data"
        );
    }

    #[test]
    fn existing_roots_do_not_require_parent_or_higher_ancestor_sync() {
        let temporary = TestDirectory::new();
        let root = temporary.0.join("parent/snapshot");
        fs::create_dir_all(&root).unwrap();
        let checked = OutputPath::new(&root).unwrap();
        let manifest = ManifestName::new("checkpoint.toml").unwrap();
        for ancestor in [root.parent().unwrap(), temporary.0.as_path()] {
            let _fault = DirectorySyncFault::new(ancestor);
            checked.ensure_durable_directory().unwrap();
            let mut writer = SnapshotWriter::acquire(&root, &[manifest]).unwrap();
            let stage = writer.begin().unwrap();
            fs::write(stage.path().join("payload"), b"complete").unwrap();
            let published = stage.publish(b"pointer", manifest).unwrap();
            assert_eq!(
                fs::read(published.path().join("payload")).unwrap(),
                b"complete"
            );
        }
    }

    #[test]
    fn parent_access_and_capability_errors_allow_initial_acquisition_and_retry() {
        for kind in [
            ErrorKind::PermissionDenied,
            ErrorKind::InvalidInput,
            ErrorKind::Unsupported,
        ] {
            for parent_exists in [false, true] {
                let temporary = TestDirectory::new();
                let parent = temporary.0.join("parent");
                if parent_exists {
                    fs::create_dir(&parent).unwrap();
                }
                let root = parent.join("snapshot");
                let _fault = DirectorySyncFault::with_kind(&parent, kind);
                let trace = DirectorySyncTrace::new();
                for attempt in 1..=2 {
                    drop(SnapshotWriter::acquire(&root, &[]).unwrap());
                    assert_eq!(trace.count(&parent), attempt);
                    assert_eq!(trace.count(&root), attempt);
                }
            }
        }
    }

    #[test]
    fn parent_io_failure_remains_fatal_on_retry() {
        for parent_exists in [false, true] {
            let temporary = TestDirectory::new();
            let parent = temporary.0.join("parent");
            if parent_exists {
                fs::create_dir(&parent).unwrap();
            }
            let root = parent.join("snapshot");
            let fault = DirectorySyncFault::with_kind(&parent, ErrorKind::Other);
            for _ in 0..2 {
                let error = SnapshotWriter::acquire(&root, &[]).unwrap_err();
                assert!(
                    matches!(error, SnapshotError::Path(OutputPathError::Io { path, source })
                    if path == parent && source.kind() == ErrorKind::Other)
                );
                assert!(root.is_dir());
                assert!(!root.join(".snapshot.lock").exists());
            }
            drop(fault);
            drop(SnapshotWriter::acquire(&root, &[]).unwrap());
        }
    }

    #[test]
    fn writer_acquisition_syncs_root_and_parent_once() {
        for preexisting in [false, true] {
            let temporary = TestDirectory::new();
            let root = temporary.0.join("snapshot");
            if preexisting {
                fs::create_dir(&root).unwrap();
            }
            let trace = DirectorySyncTrace::new();
            let _writer = SnapshotWriter::acquire(&root, &[]).unwrap();
            assert_eq!(trace.count(&root), 1);
            assert_eq!(trace.count(&temporary.0), 1);
        }
    }

    #[test]
    fn newly_created_parent_is_synchronized_once() {
        let temporary = TestDirectory::new();
        let parent = temporary.0.join("new/nested");
        let root = parent.join("snapshot");
        let trace = DirectorySyncTrace::new();
        let _writer = SnapshotWriter::acquire(&root, &[]).unwrap();
        assert_eq!(trace.count(&root), 1);
        assert_eq!(trace.count(&parent), 1);
        assert_eq!(trace.count(&temporary.0.join("new")), 1);
        assert_eq!(trace.count(&temporary.0), 0);
    }

    #[test]
    fn retirement_sync_failure_reports_the_detached_tree() {
        let temporary = TestDirectory::new();
        let root = temporary.0.join("snapshot");
        let writer = SnapshotWriter::acquire(&root, &[]).unwrap();
        fs::write(root.join("payload"), b"retained").unwrap();
        let fault = DirectorySyncFault::new(&temporary.0);
        let error = writer.retire_root().unwrap_err();
        let SnapshotError::Retired { path, source } = &error else {
            panic!("expected retirement durability error, got {error:?}");
        };
        assert!(!root.exists());
        assert_eq!(fs::read(path.join("payload")).unwrap(), b"retained");
        assert!(error.to_string().contains(&path.display().to_string()));
        assert!(matches!(source, OutputPathError::Io { path, .. } if path == &temporary.0));
        drop(fault);
        SnapshotWriter::remove_retired_roots(&root).unwrap();
        assert!(!path.exists());
    }

    #[test]
    fn removal_sync_failure_reports_that_the_detached_tree_is_gone() {
        let temporary = TestDirectory::new();
        let root = temporary.0.join("snapshot");
        let writer = SnapshotWriter::acquire(&root, &[]).unwrap();
        let retired = writer.retire_root().unwrap();
        let _fault = DirectorySyncFault::new(&temporary.0);
        let error = retired.remove().unwrap_err();
        let SnapshotError::Removed { path, source } = &error else {
            panic!("expected removal durability error, got {error:?}");
        };
        assert!(!root.exists());
        assert!(!path.exists());
        let message = error.to_string();
        assert!(message.contains("tree removed"));
        assert!(message.contains(&path.display().to_string()));
        assert!(matches!(source, OutputPathError::Io { path, .. } if path == &temporary.0));
    }

    #[test]
    fn durable_initialization_syncs_new_intermediate_directories() {
        let temporary = TestDirectory::new();
        let root = OutputPath::new(&temporary.0.join("new/nested/snapshot")).unwrap();
        let created_ancestor = temporary.0.join("new");
        let _fault = DirectorySyncFault::new(&created_ancestor);
        let error = root.ensure_durable_directory().unwrap_err();
        assert!(matches!(error, OutputPathError::Io { path, source }
            if path == created_ancestor && source.kind() == ErrorKind::PermissionDenied));
        assert!(root.as_path().is_dir());
    }

    #[test]
    fn beginning_generations_does_not_resynchronize_the_roots_parent() {
        let temporary = TestDirectory::new();
        let root = temporary.0.join("snapshot");
        let manifest = ManifestName::new("checkpoint.toml").unwrap();
        let mut writer = SnapshotWriter::acquire(&root, &[manifest]).unwrap();
        let _fault = DirectorySyncFault::new(&temporary.0);
        for contents in [b"first", b"later"] {
            let stage = writer.begin().unwrap();
            stage
                .output_path()
                .unwrap()
                .join("payload")
                .unwrap()
                .write(contents)
                .unwrap();
            let generation = stage.publish(contents, manifest).unwrap();
            assert_eq!(
                fs::read(generation.path().join("payload")).unwrap(),
                contents
            );
        }
    }

    #[test]
    fn directory_publication_does_not_synchronize_preexisting_higher_ancestors() {
        let temporary = TestDirectory::new();
        let work = temporary.0.join("work");
        fs::create_dir(&work).unwrap();
        let _fault = DirectorySyncFault::new(&temporary.0);
        let publication = DirectoryPublication::new(&work.join("site")).unwrap();
        fs::create_dir(publication.staging_path()).unwrap();
        fs::write(publication.staging_path().join("index.html"), b"complete").unwrap();
        publication.publish().unwrap();
        assert_eq!(fs::read(work.join("site/index.html")).unwrap(), b"complete");
    }

    #[test]
    fn joining_relative_paths_checks_spelling_and_preserves_boundary() {
        let temporary = TestDirectory::new();
        let root = OutputPath::new(&temporary.0).unwrap();
        let child = root
            .join_relative("name with spaces/λ/data.bin".as_ref())
            .unwrap();
        assert_eq!(child.root, root.root);
        assert_eq!(
            child.as_path(),
            temporary.0.join("name with spaces/λ/data.bin")
        );
        for invalid in ["", ".", "..", "a/../b", "a/./b", "a//b", "a/", "/a", "a\\b"] {
            let error = root.join_relative(invalid.as_ref()).unwrap_err();
            assert!(
                matches!(error, OutputPathError::InvalidRelative(_)),
                "{invalid:?}"
            );
        }
        assert!(fs::read_dir(&temporary.0).unwrap().next().is_none());
    }
}
