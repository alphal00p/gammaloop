//! Explicit reclamation of superseded immutable generations.
//!
//! For example, `cleanup state --path saved-state` retains its committed snapshot;
//! `cleanup integration --workspace-path integration-workspace` retains its checkpoint.

use super::Integrate;
use crate::state::State;
use crate::CLISettings;
use alphal00p_persistence::SnapshotError;
use clap::{Subcommand, ValueHint};
use color_eyre::{Report, Result, Section};
use gammalooprs::integrate::{IntegrationCheckpointError, IntegrationWorkspace};
use gammalooprs::INTEGRATION_MANIFEST_FILE;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use std::io::ErrorKind;
use std::path::PathBuf;
use tracing::info;

/// Remove superseded generations while retaining the validated current snapshot.
///
/// Ordinary saves and integration iterations retain every generation until explicitly cleaned up.
///
/// ```text
/// cleanup state --path saved-state
/// cleanup integration -w integration-workspace
/// cleanup integration
/// ```
/// Without a workspace path, integration cleanup uses the same default as `integrate`.
#[derive(Subcommand, Debug, Serialize, Deserialize, Clone, JsonSchema, PartialEq)]
pub enum Cleanup {
    /// Prune older saved-state generations; keep the current committed state.
    State {
        /// Saved-state folder; defaults to the active state folder.
        #[arg(short, long, value_hint = ValueHint::DirPath)]
        path: Option<PathBuf>,
    },
    /// Prune older integration generations; keep the last committed iteration.
    Integration {
        /// Integration workspace; defaults to the same workspace as `integrate`.
        #[arg(short = 'w', long, value_hint = ValueHint::DirPath)]
        workspace_path: Option<PathBuf>,
    },
}

impl Cleanup {
    /// Acquire the destination's writer lock, validate its snapshot, and prune it.
    ///
    /// See [`Self`] for command syntax. Every removed generation
    /// is reported; presentation files and unrelated entries are retained.
    ///
    /// # Errors
    /// Fails if the destination lies in a read-only active state, another writer
    /// owns it, its current snapshot is invalid, or filesystem operations fail.
    pub fn run(self, settings: &CLISettings) -> Result<()> {
        let (root, removed) = match self {
            Self::State { path } => {
                let root = path.unwrap_or_else(|| settings.state.folder.clone());
                settings.ensure_write_target_outside_active_state(
                    &root,
                    "clean up saved-state generations",
                )?;
                let removed = State::cleanup_generations(&root).map_err(|error| {
                    let busy = error.is_busy();
                    Self::diagnostic(Report::new(error), busy)
                })?;
                (root, removed)
            }
            Self::Integration { workspace_path } => {
                let workspace_path =
                    workspace_path.unwrap_or_else(|| Integrate::default_workspace_path(settings));
                settings.ensure_write_target_outside_active_state(
                    &workspace_path,
                    "clean up integration generations",
                )?;
                let removed = IntegrationWorkspace::cleanup_superseded_generations(&workspace_path)
                    .map_err(|error| {
                        let missing_workspace = matches!(
                            &error,
                            IntegrationCheckpointError::Io { path, source }
                                if path.ends_with(INTEGRATION_MANIFEST_FILE)
                                    && source.kind() == ErrorKind::NotFound
                        );
                        let busy = error.is_busy();
                        let report = Self::diagnostic(Report::new(error), busy);
                        if missing_workspace {
                            report.wrap_err(format!(
                                "No integration workspace at '{}'",
                                workspace_path.display()
                            ))
                        } else {
                            report
                        }
                    })?;
                (workspace_path, removed)
            }
        };
        for generation in &removed {
            info!("Removed generation {generation} from '{}'", root.display());
        }
        info!(
            "Removed {} superseded generation(s) from '{}'; retained the current generation",
            removed.len(),
            root.display()
        );
        Ok(())
    }

    /// Attach the same writer-contention guidance to whole-state removal.
    pub(crate) fn snapshot_diagnostic(error: SnapshotError) -> Report {
        let busy = error.is_busy();
        Self::diagnostic(Report::new(error), busy)
    }

    fn diagnostic(error: Report, busy: bool) -> Report {
        if busy {
            error.suggestion("Another GammaLoop process is writing; retry when it finishes.")
        } else {
            error
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Cleanup;
    use crate::commands::Commands;
    use crate::state::StateSnapshotError;
    use crate::CLISettings;
    use alphal00p_persistence::{SnapshotWriter, SNAPSHOT_LOCK_FILE};
    use color_eyre::Report;
    use gammalooprs::initialisation;
    use gammalooprs::integrate::IntegrationCheckpointError;
    use std::io::ErrorKind;
    use std::path::PathBuf;

    #[test]
    fn cleanup_commands_parse_both_targets() {
        assert_eq!(
            "cleanup state --path saved".parse::<Commands>().unwrap(),
            Commands::Cleanup(Cleanup::State {
                path: Some(PathBuf::from("saved")),
            })
        );
        assert_eq!(
            "cleanup integration --workspace-path run"
                .parse::<Commands>()
                .unwrap(),
            Commands::Cleanup(Cleanup::Integration {
                workspace_path: Some(PathBuf::from("run")),
            })
        );
        assert_eq!(
            "cleanup integration -w run".parse::<Commands>().unwrap(),
            Commands::Cleanup(Cleanup::Integration {
                workspace_path: Some(PathBuf::from("run")),
            })
        );
        assert_eq!(
            "cleanup integration".parse::<Commands>().unwrap(),
            Commands::Cleanup(Cleanup::Integration {
                workspace_path: None
            })
        );
    }

    #[test]
    fn cleanup_integration_uses_the_active_state_workspace_default() {
        let temporary = tempfile::tempdir().unwrap();
        let mut settings = CLISettings::default();
        settings.state.folder = temporary.path().join("saved");
        let error = Cleanup::Integration {
            workspace_path: None,
        }
        .run(&settings)
        .unwrap_err();
        let workspace = settings.state.folder.join("integration_workspace");
        let expected = workspace.join("manifest.json");
        assert_eq!(
            error.to_string(),
            format!("No integration workspace at '{}'", workspace.display())
        );
        assert!(matches!(
            error.downcast_ref::<IntegrationCheckpointError>(),
            Some(IntegrationCheckpointError::Io { path, source })
                if path == &expected && source.kind() == ErrorKind::NotFound
        ));
        assert!(!format!("{error:?}").contains("retry when it finishes"));
        assert!(!settings.state.folder.exists());
    }

    #[test]
    fn cleanup_busy_diagnostics_preserve_sources_and_suggest_retrying() {
        initialisation::test_initialise().unwrap();
        let temporary = tempfile::tempdir().unwrap();
        let _writer = SnapshotWriter::acquire(temporary.path(), &[]).unwrap();
        let state_error = StateSnapshotError::from(Box::new(
            SnapshotWriter::acquire(temporary.path(), &[]).unwrap_err(),
        ));
        let integration_error = IntegrationCheckpointError::from(
            SnapshotWriter::acquire(temporary.path(), &[]).unwrap_err(),
        );
        let state_busy = state_error.is_busy();
        let integration_busy = integration_error.is_busy();
        for (cause, busy) in [
            (Report::new(state_error), state_busy),
            (Report::new(integration_error), integration_busy),
        ] {
            assert!(busy);
            let message = cause.to_string();
            let error = Cleanup::diagnostic(cause, busy);
            assert_eq!(error.to_string(), message);
            assert!(message.contains(SNAPSHOT_LOCK_FILE), "{message}");
            let diagnostic = format!("{error:?}");
            assert_eq!(
                diagnostic.matches("retry when it finishes").count(),
                1,
                "{diagnostic}"
            );
            assert!(
                diagnostic.contains("Another GammaLoop process is writing"),
                "{diagnostic}"
            );
            assert!(
                error.downcast_ref::<StateSnapshotError>().is_some()
                    || error.downcast_ref::<IntegrationCheckpointError>().is_some()
            );
        }
    }

    #[test]
    fn cleanup_honors_read_only_active_state_before_creating_a_lock() {
        let temporary = tempfile::tempdir().unwrap();
        let mut settings = CLISettings::default();
        settings.state.folder = temporary.path().join("saved");
        settings.session.read_only_state = true;
        for command in [
            Cleanup::State { path: None },
            Cleanup::Integration {
                workspace_path: Some(settings.state.folder.join("integration_workspace")),
            },
        ] {
            let error = command.run(&settings).unwrap_err();
            assert!(error.to_string().contains("--read-only-state"));
            assert!(!settings.state.folder.exists());
        }
    }
}
