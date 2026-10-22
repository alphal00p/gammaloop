use std::error::Error as StdError;
use std::io::Error as IoError;
use std::path::{Path, PathBuf};

use bincode::error::DecodeError;
use eyre::Report;
use feynkit_generator::GenerationType;
use thiserror::Error;
use toml::de::Error as SettingsDecodeError;

use super::{OutputPathError, ProcessLoadIntegrandSelector};

/// Actionable failures at the indexed process snapshot boundary.
///
/// Membership failures are distinct from path rejection, definition decoding,
/// and nested payload failures. Nested serializers retain their original error
/// objects and source chains. See [`super::ProcessSnapshotIndex::validate`] for
/// an example of matching a validation failure.
#[derive(Debug, Error)]
pub enum ProcessSnapshotError {
    #[error(transparent)]
    Output(#[from] OutputPathError),
    #[error("Duplicate saved process {process:?}")]
    DuplicateProcess { process: String },
    #[error("Duplicate or unordered saved child {child:?} after {previous:?} in {process:?}")]
    InvalidChildOrder {
        process: String,
        previous: String,
        child: String,
    },
    #[error("Process {process:?} has ID {actual}, expected {expected}")]
    InvalidProcessId {
        process: String,
        expected: usize,
        actual: usize,
    },
    #[error("Process {process:?} has kind {actual}, expected {expected}")]
    ProcessKindMismatch {
        process: String,
        expected: GenerationType,
        actual: GenerationType,
    },
    #[error("Saved process name {actual:?} differs from indexed name {expected:?}")]
    ProcessNameMismatch { expected: String, actual: String },
    #[error("Child name {actual:?} differs from {expected:?} in process {process:?}")]
    ChildNameMismatch {
        process: String,
        expected: String,
        actual: String,
    },
    #[error("Generated {kind} payload disagrees with child {child:?} in process {process:?}")]
    GeneratedChildMismatch {
        process: String,
        child: String,
        kind: GenerationType,
    },
    #[error("Cannot load requested process {process:?}; available processes: {available:?}")]
    UnknownProcess {
        process: String,
        available: Vec<String>,
    },
    #[error("Cannot load requested integrand {selector:?}")]
    UnknownIntegrand {
        selector: ProcessLoadIntegrandSelector,
    },
    #[error("Could not read process snapshot member {path}")]
    Read {
        path: PathBuf,
        #[source]
        source: IoError,
    },
    #[error("Could not decode process definition {path}")]
    DefinitionDecode {
        path: PathBuf,
        #[source]
        source: DecodeError,
    },
    #[error("Could not decode process settings history {path}")]
    SettingsDecode {
        path: PathBuf,
        #[source]
        source: SettingsDecodeError,
    },
    #[error("Could not {operation} for process {process:?} at {path}")]
    Payload {
        operation: &'static str,
        process: String,
        path: PathBuf,
        #[source]
        source: Box<dyn StdError + Send + Sync>,
    },
}

impl ProcessSnapshotError {
    pub(super) fn payload(
        operation: &'static str,
        process: &str,
        path: &Path,
        source: Report,
    ) -> Self {
        Self::Payload {
            operation,
            process: process.to_owned(),
            path: path.to_owned(),
            // eyre 0.6 boxes an opaque ErrorImpl whose source skips its root.
            // Add our operation context so the original native root remains exposed.
            source: source.wrap_err(operation).into(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::ErrorKind;

    fn native_source<E: StdError + 'static>(error: &ProcessSnapshotError) -> &E {
        let mut source = error.source();
        while let Some(cause) = source {
            if let Some(native) = cause.downcast_ref::<E>() {
                return native;
            }
            source = cause.source();
        }
        panic!("Native error absent from source chain: {error:?}");
    }

    #[test]
    fn process_snapshot_payload_exposes_bare_io() {
        let report = Report::new(IoError::from(ErrorKind::PermissionDenied));
        let error =
            ProcessSnapshotError::payload("write payload", "process", Path::new("stage"), report);
        assert_eq!(
            native_source::<IoError>(&error).kind(),
            ErrorKind::PermissionDenied
        );
    }

    #[test]
    fn process_snapshot_payload_exposes_contextual_io() {
        let report = Report::new(IoError::from(ErrorKind::NotFound)).wrap_err("read child");
        let error =
            ProcessSnapshotError::payload("load payload", "process", Path::new("stage"), report);
        assert_eq!(native_source::<IoError>(&error).kind(), ErrorKind::NotFound);
    }

    #[test]
    fn process_snapshot_payload_exposes_contextual_decode() {
        let report = Report::new(DecodeError::UnexpectedEnd { additional: 5 })
            .wrap_err("decode child metadata");
        let error =
            ProcessSnapshotError::payload("load payload", "process", Path::new("stage"), report);
        assert!(matches!(
            native_source::<DecodeError>(&error),
            DecodeError::UnexpectedEnd { additional: 5 }
        ));
    }
}
