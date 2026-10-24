use std::env;
use std::fs;
use std::io::{Error as IoError, ErrorKind};
use std::path::{Component, Path, PathBuf};
use std::sync::OnceLock;
use thiserror::Error;

/// A required file or directory relative to a workspace root.
///
/// See [`WorkspaceRoot::discover`] for an example combining both marker kinds.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum WorkspaceMarker {
    /// A regular file, following symbolic links.
    File(&'static str),
    /// A directory, following symbolic links.
    Directory(&'static str),
}

impl WorkspaceMarker {
    fn path(&self) -> &Path {
        let (Self::File(path) | Self::Directory(path)) = self;
        Path::new(path)
    }

    fn matches(&self, root: &Path) -> Result<bool, WorkspaceError> {
        let path = root.join(self.path());
        match fs::metadata(&path) {
            Ok(metadata) => Ok(match self {
                Self::File(_) => metadata.is_file(),
                Self::Directory(_) => metadata.is_dir(),
            }),
            Err(source)
                if matches!(
                    source.kind(),
                    ErrorKind::NotFound | ErrorKind::NotADirectory
                ) =>
            {
                Ok(false)
            }
            Err(source) => Err(WorkspaceError::Inspect { path, source }),
        }
    }
}

/// A canonical directory containing `Cargo.toml` and the requested workspace markers.
///
/// Discovery never creates or modifies files. See [`Self::discover`] for an example.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct WorkspaceRoot(PathBuf);

const TEST_WORKSPACE_MARKERS: &[WorkspaceMarker] = &[
    WorkspaceMarker::Directory("tests/resources"),
    WorkspaceMarker::Directory("examples/cli"),
];

/// Borrows the shared GammaLoop test root, discovering it once on first use.
///
/// Requires `Cargo.toml`, `tests/resources/`, and `examples/cli/`. Discovery follows
/// [`WorkspaceRoot::discover_runtime`], with the caller's compile-time manifest
/// directory as the final fallback. The first result, including an error, is
/// retained across callers even if runtime paths later change.
///
/// # Panics
///
/// Panics with the discovery error if no usable test workspace can be found.
///
/// ```no_run
/// use alphal00p_workspace as workspace;
/// use std::path::Path;
///
/// let root = workspace::test_workspace_root(Path::new(env!("CARGO_MANIFEST_DIR")));
/// assert!(root.join("tests/resources").is_dir());
/// assert!(root.join("examples/cli").is_dir());
/// ```
#[must_use]
pub fn test_workspace_root(compiled_manifest_dir: &Path) -> &'static Path {
    static ROOT: OnceLock<Result<WorkspaceRoot, WorkspaceError>> = OnceLock::new();

    let root = ROOT.get_or_init(|| {
        WorkspaceRoot::discover_runtime(compiled_manifest_dir, TEST_WORKSPACE_MARKERS)
    });
    root.as_ref()
        .expect("GammaLoop tests require a runtime workspace containing their fixtures")
        .as_path()
}

/// Why workspace discovery failed.
///
/// I/O failures retain their source and the affected path. See
/// [`WorkspaceRoot::discover`] for an example propagating this error.
#[derive(Debug, Error)]
pub enum WorkspaceError {
    /// The current directory could not be read when no runtime manifest was supplied.
    #[error("cannot read the current directory for workspace discovery")]
    CurrentDirectory {
        /// Underlying operating-system failure.
        #[source]
        source: IoError,
    },
    /// A marker was empty, absolute, or contained unsupported path components.
    #[error("workspace marker must be a nonempty relative path without traversal: {path}")]
    InvalidMarker {
        /// Rejected caller-supplied marker.
        path: PathBuf,
    },
    /// A discovery path could not be canonicalized or inspected.
    #[error("cannot inspect workspace discovery path {path}")]
    Inspect {
        /// Path whose inspection failed.
        path: PathBuf,
        /// Underlying filesystem failure.
        #[source]
        source: IoError,
    },
    /// The requested start exists but is not a directory.
    #[error("workspace discovery must start from a directory: {path}")]
    NotDirectory {
        /// Canonical start path.
        path: PathBuf,
    },
    /// No ancestor of any searched directory matched all markers.
    #[error("no workspace containing Cargo.toml and the requested markers above {starts:?}")]
    NotFound {
        /// Searched starting directories in selection order.
        starts: Vec<PathBuf>,
    },
}

impl WorkspaceRoot {
    /// Finds the nearest matching ancestor of an existing starting directory.
    ///
    /// Every root must contain a `Cargo.toml` file and all supplied markers. Markers
    /// must be nonempty relative paths without traversal. The start is canonicalized
    /// before searching, so a nonexistent start cannot accidentally select an existing
    /// ancestor. No other runtime paths are searched.
    ///
    /// # Errors
    ///
    /// Returns [`WorkspaceError`] for invalid markers, an inaccessible or nondirectory
    /// start, an I/O failure inspecting a candidate, or no matching ancestor.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use alphal00p_workspace::{WorkspaceError, WorkspaceMarker, WorkspaceRoot};
    /// use std::path::Path;
    ///
    /// let root = WorkspaceRoot::discover(
    ///     Path::new("crates/alphal00p-workspace"),
    ///     &[
    ///         WorkspaceMarker::File("docs/products/registry.toml"),
    ///         WorkspaceMarker::Directory("tests/resources"),
    ///     ],
    /// )?;
    /// assert!(root.as_path().join("Cargo.toml").is_file());
    /// let owned_path = root.into_path_buf();
    /// assert!(owned_path.is_absolute());
    /// # Ok::<(), WorkspaceError>(())
    /// ```
    pub fn discover(start: &Path, markers: &[WorkspaceMarker]) -> Result<Self, WorkspaceError> {
        for marker in markers {
            let path = marker.path();
            if path.as_os_str().is_empty()
                || !path
                    .components()
                    .all(|part| matches!(part, Component::Normal(_)))
            {
                return Err(WorkspaceError::InvalidMarker {
                    path: path.to_path_buf(),
                });
            }
        }
        let start = start
            .canonicalize()
            .map_err(|source| WorkspaceError::Inspect {
                path: start.to_path_buf(),
                source,
            })?;
        let metadata = fs::metadata(&start).map_err(|source| WorkspaceError::Inspect {
            path: start.to_path_buf(),
            source,
        })?;
        if !metadata.is_dir() {
            return Err(WorkspaceError::NotDirectory { path: start });
        }
        for ancestor in start.ancestors() {
            if !WorkspaceMarker::File("Cargo.toml").matches(ancestor)? {
                continue;
            }
            let mut matches = true;
            for marker in markers {
                if !marker.matches(ancestor)? {
                    matches = false;
                    break;
                }
            }
            if matches {
                return Ok(Self(ancestor.to_path_buf()));
            }
        }
        Err(WorkspaceError::NotFound {
            starts: vec![start],
        })
    }

    /// Discovers from `CARGO_MANIFEST_DIR`, or the current directory when unset.
    ///
    /// A supplied manifest directory is authoritative, even when invalid. Unlike
    /// [`Self::discover_runtime`], this never searches executable or compile-time paths.
    /// Repository tools must stay in the selected tree: those fallbacks could silently
    /// read a different checkout from the one the caller intends to inspect or build.
    ///
    /// # Errors
    ///
    /// Returns [`WorkspaceError::CurrentDirectory`] if no manifest directory is set
    /// and the current directory cannot be read, or an error from [`Self::discover`].
    ///
    /// ```no_run
    /// use alphal00p_workspace::{WorkspaceError, WorkspaceMarker, WorkspaceRoot};
    ///
    /// let root = WorkspaceRoot::discover_from_manifest_or_cwd(&[
    ///     WorkspaceMarker::File("docs/products/registry.toml"),
    /// ])?;
    /// assert!(root.as_path().join("docs/products/registry.toml").is_file());
    /// # Ok::<(), WorkspaceError>(())
    /// ```
    pub fn discover_from_manifest_or_cwd(
        markers: &[WorkspaceMarker],
    ) -> Result<Self, WorkspaceError> {
        let runtime_manifest = env::var_os("CARGO_MANIFEST_DIR").map(PathBuf::from);
        Self::discover_manifest_or_cwd(runtime_manifest.as_deref(), env::current_dir, markers)
    }

    /// Discovers a test workspace, honoring Cargo/Nextest runtime path remapping.
    ///
    /// An explicitly set `CARGO_MANIFEST_DIR` is authoritative, even when invalid.
    /// Otherwise, searches the current directory, executable's directory, then the
    /// caller's compile-time manifest directory. Unavailable current paths and missing
    /// compile-time paths are skipped, supporting binaries relocated from build machines.
    /// Only test discovery uses these fallbacks because archived binaries may run outside
    /// their source tree; repository tools use [`Self::discover_from_manifest_or_cwd`].
    ///
    /// # Errors
    ///
    /// Returns [`WorkspaceError`] if the selected runtime manifest directory is invalid,
    /// no candidate contains the requested workspace, markers are invalid, or inspection
    /// fails for a reason other than a missing optional candidate.
    ///
    /// # Examples
    ///
    /// ```no_run
    /// use alphal00p_workspace::{WorkspaceError, WorkspaceMarker, WorkspaceRoot};
    /// use std::path::Path;
    ///
    /// let root = WorkspaceRoot::discover_runtime(
    ///     Path::new(env!("CARGO_MANIFEST_DIR")),
    ///     &[WorkspaceMarker::Directory("tests/resources")],
    /// )?;
    /// assert!(root.as_path().join("tests/resources").is_dir());
    /// # Ok::<(), WorkspaceError>(())
    /// ```
    pub fn discover_runtime(
        compiled_manifest_dir: &Path,
        markers: &[WorkspaceMarker],
    ) -> Result<Self, WorkspaceError> {
        let runtime_manifest = env::var_os("CARGO_MANIFEST_DIR").map(PathBuf::from);
        let current_dir = env::current_dir().ok();
        let executable = env::current_exe().ok();
        let candidates = [
            current_dir.as_deref(),
            executable.as_deref().and_then(Path::parent),
            Some(compiled_manifest_dir),
        ];
        Self::discover_candidates(
            runtime_manifest.as_deref(),
            candidates.into_iter().flatten(),
            markers,
        )
    }

    /// Borrows the canonical root directory; see [`Self::discover`] for an example.
    #[must_use]
    pub fn as_path(&self) -> &Path {
        &self.0
    }

    /// Takes ownership of the canonical root; see [`Self::discover`] for an example.
    #[must_use]
    pub fn into_path_buf(self) -> PathBuf {
        self.0
    }

    fn discover_manifest_or_cwd(
        runtime_manifest: Option<&Path>,
        current_dir: impl FnOnce() -> Result<PathBuf, IoError>,
        markers: &[WorkspaceMarker],
    ) -> Result<Self, WorkspaceError> {
        if let Some(start) = runtime_manifest {
            return Self::discover(start, markers);
        }
        let start = current_dir().map_err(|source| WorkspaceError::CurrentDirectory { source })?;
        Self::discover(&start, markers)
    }

    fn discover_candidates<'a>(
        runtime_manifest: Option<&Path>,
        candidates: impl IntoIterator<Item = &'a Path>,
        markers: &[WorkspaceMarker],
    ) -> Result<Self, WorkspaceError> {
        if let Some(start) = runtime_manifest {
            return Self::discover(start, markers);
        }
        let mut starts = Vec::new();
        for start in candidates {
            match Self::discover(start, markers) {
                Ok(root) => return Ok(root),
                Err(WorkspaceError::NotFound { .. }) => {}
                Err(WorkspaceError::Inspect { source, .. })
                    if source.kind() == ErrorKind::NotFound => {}
                Err(error) => return Err(error),
            }
            starts.push(start.to_path_buf());
        }
        Err(WorkspaceError::NotFound { starts })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error as _;
    use tempfile::{self, TempDir};

    const MARKERS: &[WorkspaceMarker] = &[
        WorkspaceMarker::File("docs/products/registry.toml"),
        WorkspaceMarker::Directory("tests/resources"),
    ];

    struct Fixture(TempDir);

    impl Fixture {
        fn new() -> Self {
            Self(tempfile::tempdir().unwrap())
        }

        fn workspace(&self, relative: &str) -> PathBuf {
            let root = self.0.path().join(relative);
            fs::create_dir_all(root.join("docs/products")).unwrap();
            fs::create_dir_all(root.join("tests/resources")).unwrap();
            fs::write(root.join("Cargo.toml"), "[workspace]\n").unwrap();
            fs::write(root.join("docs/products/registry.toml"), "").unwrap();
            root.canonicalize().unwrap()
        }
    }

    #[test]
    fn discovers_nearest_matching_root() {
        let fixture = Fixture::new();
        let outer = fixture.workspace("outer");
        let inner = fixture.workspace("outer/inner");
        let start = inner.join("tests/resources");
        assert_eq!(
            WorkspaceRoot::discover(&start, MARKERS).unwrap().as_path(),
            inner
        );
        fs::remove_file(inner.join("docs/products/registry.toml")).unwrap();
        assert_eq!(
            WorkspaceRoot::discover(&start, MARKERS)
                .unwrap()
                .into_path_buf(),
            outer
        );
    }

    #[test]
    fn rejects_wrong_marker_types_and_missing_workspace() {
        let fixture = Fixture::new();
        assert!(matches!(
            WorkspaceRoot::discover(fixture.0.path(), MARKERS),
            Err(WorkspaceError::NotFound { .. })
        ));
        let root = fixture.workspace("workspace");
        let marker = root.join("docs/products/registry.toml");
        fs::remove_file(&marker).unwrap();
        fs::create_dir(&marker).unwrap();
        assert!(matches!(
            WorkspaceRoot::discover(&root, MARKERS),
            Err(WorkspaceError::NotFound { .. })
        ));
        fs::remove_dir(&marker).unwrap();
        fs::write(marker, "").unwrap();
        fs::remove_dir(root.join("tests/resources")).unwrap();
        fs::write(root.join("tests/resources"), "").unwrap();
        assert!(matches!(
            WorkspaceRoot::discover(&root, MARKERS),
            Err(WorkspaceError::NotFound { .. })
        ));
    }

    #[test]
    fn shared_test_markers_require_the_manifest_and_both_fixture_directories() {
        let fixture = Fixture::new();
        let root = fixture.workspace("workspace");
        let cli = root.join("examples/cli");
        fs::create_dir_all(cli.parent().unwrap()).unwrap();
        fs::write(&cli, "not a directory").unwrap();
        assert!(matches!(
            WorkspaceRoot::discover(&root, TEST_WORKSPACE_MARKERS),
            Err(WorkspaceError::NotFound { .. })
        ));
        fs::remove_file(&cli).unwrap();
        fs::create_dir(&cli).unwrap();
        assert_eq!(
            WorkspaceRoot::discover(&root, TEST_WORKSPACE_MARKERS)
                .unwrap()
                .as_path(),
            root
        );

        for relative in ["Cargo.toml", "tests/resources", "examples/cli"] {
            let member = root.join(relative);
            let moved = root.join("temporarily-removed");
            fs::rename(&member, &moved).unwrap();
            assert!(
                matches!(
                    WorkspaceRoot::discover(&root, TEST_WORKSPACE_MARKERS),
                    Err(WorkspaceError::NotFound { .. })
                ),
                "test discovery must require {relative}"
            );
            fs::rename(moved, member).unwrap();
        }
    }

    #[test]
    fn rejects_missing_or_nondirectory_requested_start() {
        let fixture = Fixture::new();
        let root = fixture.workspace("workspace");
        let error = WorkspaceRoot::discover(&root.join("missing"), MARKERS).unwrap_err();
        assert!(error.source().is_some());
        let WorkspaceError::Inspect { path, source } = error else {
            panic!("missing discovery start should retain its inspection error");
        };
        assert_eq!(path, root.join("missing"));
        assert_eq!(source.kind(), ErrorKind::NotFound);
        assert!(matches!(
            WorkspaceRoot::discover(&root.join("Cargo.toml"), MARKERS),
            Err(WorkspaceError::NotDirectory { .. })
        ));
    }

    #[test]
    fn rejects_markers_outside_the_root() {
        let fixture = Fixture::new();
        for marker in ["", ".", "../Cargo.toml", "/Cargo.toml"] {
            assert!(matches!(
                WorkspaceRoot::discover(fixture.0.path(), &[WorkspaceMarker::File(marker)]),
                Err(WorkspaceError::InvalidMarker { .. })
            ));
        }
    }

    #[test]
    fn manifest_or_cwd_never_falls_back_from_a_selected_manifest() {
        let fixture = Fixture::new();
        let root = fixture.workspace("selected");
        let discover = |start: &Path| {
            WorkspaceRoot::discover_manifest_or_cwd(
                Some(start),
                || panic!("an explicit manifest must not read cwd"),
                MARKERS,
            )
        };
        assert_eq!(
            discover(&root.join("tests/resources")).unwrap().as_path(),
            root
        );
        let missing = root.join("missing");
        let WorkspaceError::Inspect { path, source } = discover(&missing).unwrap_err() else {
            panic!("a missing manifest directory should report inspection failure");
        };
        assert_eq!(path, missing);
        assert_eq!(source.kind(), ErrorKind::NotFound);
        let manifest = root.join("Cargo.toml");
        let error = discover(&manifest).unwrap_err();
        assert!(matches!(error, WorkspaceError::NotDirectory { path } if path == manifest));
        let error = discover(fixture.0.path()).unwrap_err();
        let WorkspaceError::NotFound { starts } = error else {
            panic!("a markerless manifest directory should report no matching workspace");
        };
        assert_eq!(starts, [fixture.0.path().canonicalize().unwrap()]);
    }

    #[test]
    fn manifest_or_cwd_searches_only_the_current_directory_when_unset() {
        let fixture = Fixture::new();
        let root = fixture.workspace("workspace");
        let discovered = WorkspaceRoot::discover_manifest_or_cwd(
            None,
            || Ok(root.join("tests/resources")),
            MARKERS,
        )
        .unwrap();
        assert_eq!(discovered.as_path(), root);

        let error = WorkspaceRoot::discover_manifest_or_cwd(
            None,
            || Ok(fixture.0.path().to_path_buf()),
            MARKERS,
        )
        .unwrap_err();
        let WorkspaceError::NotFound { starts } = error else {
            panic!("a markerless cwd must not fall back to the compile-time workspace");
        };
        assert_eq!(starts, [fixture.0.path().canonicalize().unwrap()]);
    }

    #[test]
    fn manifest_or_cwd_preserves_current_directory_errors() {
        let error = WorkspaceRoot::discover_manifest_or_cwd(
            None,
            || {
                Err(IoError::new(
                    ErrorKind::PermissionDenied,
                    "cwd is inaccessible",
                ))
            },
            MARKERS,
        )
        .unwrap_err();
        assert!(error.source().is_some());
        let WorkspaceError::CurrentDirectory { source } = error else {
            panic!("reading cwd should retain its distinct error type");
        };
        assert_eq!(source.kind(), ErrorKind::PermissionDenied);
        assert_eq!(source.to_string(), "cwd is inaccessible");
    }

    #[test]
    fn runtime_manifest_is_authoritative_even_when_invalid() {
        let fixture = Fixture::new();
        let selected = fixture.workspace("selected");
        let other = fixture.workspace("other");
        let discover = |start: &Path| {
            WorkspaceRoot::discover_candidates(Some(start), [other.as_path()], MARKERS)
        };
        assert_eq!(discover(&selected).unwrap().as_path(), selected);
        assert!(matches!(
            discover(&selected.join("missing")),
            Err(WorkspaceError::Inspect { .. })
        ));
        assert!(matches!(
            discover(fixture.0.path()),
            Err(WorkspaceError::NotFound { .. })
        ));
    }

    #[test]
    fn portable_runtime_candidates_support_current_executable_and_compile_paths() {
        let fixture = Fixture::new();
        let root = fixture.workspace("workspace");
        let binary_dir = root.join("target/debug/deps");
        fs::create_dir_all(&binary_dir).unwrap();
        let crate_dir = root.join("crates/example");
        fs::create_dir_all(&crate_dir).unwrap();
        let missing = fixture.0.path().join("remapped/compile/path");
        for candidates in [
            [root.as_path(), missing.as_path(), missing.as_path()],
            [fixture.0.path(), binary_dir.as_path(), missing.as_path()],
            [fixture.0.path(), missing.as_path(), crate_dir.as_path()],
        ] {
            assert_eq!(
                WorkspaceRoot::discover_candidates(None, candidates, MARKERS)
                    .unwrap()
                    .as_path(),
                root
            );
        }
    }

    #[test]
    fn runtime_candidates_prefer_current_then_executable_then_compile_directory() {
        let fixture = Fixture::new();
        let current = fixture.workspace("current");
        let executable = fixture.workspace("executable");
        let compiled = fixture.workspace("compiled");
        let missing = fixture.0.path().join("missing");
        for (candidates, expected) in [
            (
                [current.as_path(), executable.as_path(), compiled.as_path()],
                &current,
            ),
            (
                [missing.as_path(), executable.as_path(), compiled.as_path()],
                &executable,
            ),
            (
                [missing.as_path(), missing.as_path(), compiled.as_path()],
                &compiled,
            ),
        ] {
            let root = WorkspaceRoot::discover_candidates(None, candidates, MARKERS).unwrap();
            assert_eq!(root.as_path(), expected);
        }
    }

    #[test]
    fn runtime_candidates_propagate_nondirectory_errors_before_valid_fallbacks() {
        let fixture = Fixture::new();
        let root = fixture.workspace("workspace");
        let manifest = root.join("Cargo.toml");
        let error =
            WorkspaceRoot::discover_candidates(None, [manifest.as_path(), root.as_path()], MARKERS)
                .unwrap_err();
        assert!(matches!(error, WorkspaceError::NotDirectory { path } if path == manifest));
    }
}
