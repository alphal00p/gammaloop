//! Read-only discovery of development workspaces from explicit paths or test runtime paths.
//!
//! See [`WorkspaceRoot::discover`] for explicit discovery and
//! [`WorkspaceRoot::discover_from_manifest_or_cwd`] for repository tools.
//! [`test_workspace_root`] shares cached fixture discovery across GammaLoop tests.

mod discovery;

pub use discovery::{WorkspaceError, WorkspaceMarker, WorkspaceRoot, test_workspace_root};
