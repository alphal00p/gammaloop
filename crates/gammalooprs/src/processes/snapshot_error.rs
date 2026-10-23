use std::io::Error as IoError;
use std::path::{Path, PathBuf};

use alphal00p_persistence::OutputPathError;
use bincode::error::{DecodeError, EncodeError};
use feynkit_generator::GenerationType;
use thiserror::Error;
use toml::de::Error as SettingsDecodeError;
use toml::ser::Error as SettingsEncodeError;

use super::ProcessLoadIntegrandSelector;

/// Actionable failures at the indexed process snapshot boundary.
///
/// Membership failures are distinct from path rejection, definition decoding,
/// and nested payload failures. Nested serializers retain their original error
/// objects and source chains. See [`super::ProcessSnapshotIndex::validate`] for
/// an example of matching a validation failure.
#[derive(Debug, Error)]
pub enum ProcessSnapshotError {
    /// A managed snapshot path was rejected or could not be inspected or created.
    #[error(transparent)]
    Output(#[from] OutputPathError),
    /// Two indexed processes share a folder name, including across process kinds.
    #[error("Duplicate saved process {process:?}")]
    DuplicateProcess {
        /// Folder name repeated in the snapshot index.
        process: String,
    },
    /// Child names are not strictly increasing; duplicate names are also invalid.
    #[error("Duplicate or unordered saved child {child:?} after {previous:?} in {process:?}")]
    InvalidChildOrder {
        /// Folder name of the process owning the children.
        process: String,
        /// Child name immediately preceding the invalid entry.
        previous: String,
        /// Name that must sort strictly after `previous`.
        child: String,
    },
    /// A process ID differs from its position in the ordered snapshot membership.
    #[error("Process {process:?} has ID {actual}, expected {expected}")]
    InvalidProcessId {
        /// Folder name of the process with the inconsistent ID.
        process: String,
        /// Zero-based position in the process list or snapshot index.
        expected: usize,
        /// ID recorded in the process definition.
        actual: usize,
    },
    /// The process definition disagrees with its collection or indexed kind.
    #[error("Process {process:?} has kind {actual}, expected {expected}")]
    ProcessKindMismatch {
        /// Folder name of the inconsistent process.
        process: String,
        /// Kind required by the in-memory collection or snapshot index.
        expected: GenerationType,
        /// Kind recorded in the process definition.
        actual: GenerationType,
    },
    /// The decoded definition names a different process than its index entry.
    #[error("Saved process name {actual:?} differs from indexed name {expected:?}")]
    ProcessNameMismatch {
        /// Folder name required by the snapshot index.
        expected: String,
        /// Folder name recorded in the decoded process definition.
        actual: String,
    },
    /// A child's metadata name differs from its collection key or index entry.
    #[error("Child name {actual:?} differs from {expected:?} in process {process:?}")]
    ChildNameMismatch {
        /// Folder name of the process owning the child.
        process: String,
        /// Name required by the collection key or snapshot index.
        expected: String,
        /// Name recorded in the child metadata.
        actual: String,
    },
    /// A generated integrand disagrees with its owning child's name or process kind.
    #[error("Generated {kind} payload disagrees with child {child:?} in process {process:?}")]
    GeneratedChildMismatch {
        /// Folder name of the process owning the generated integrand.
        process: String,
        /// Name required by the owning child.
        child: String,
        /// Integrand kind required by the owning process collection.
        kind: GenerationType,
    },
    /// A requested process name is absent from the snapshot index.
    #[error("Cannot load requested process {process:?}; available processes: {available:?}")]
    UnknownProcess {
        /// Requested process folder name.
        process: String,
        /// Process folder names present in the snapshot index.
        available: Vec<String>,
    },
    /// No indexed child matches the requested qualified or unqualified selector.
    #[error("Cannot load requested integrand {selector:?}")]
    UnknownIntegrand {
        /// Selector for which the snapshot contains no matching integrand.
        selector: ProcessLoadIntegrandSelector,
    },
    /// An indexed process definition or settings-history file could not be read.
    #[error("Could not read process snapshot member {path}")]
    Read {
        /// Snapshot member whose read failed.
        path: PathBuf,
        /// Original filesystem or text-decoding failure.
        #[source]
        source: IoError,
    },
    /// A saved process definition could not be decoded from its binary representation.
    #[error("Could not decode process definition {path}")]
    DefinitionDecode {
        /// Definition file whose bytes could not be decoded.
        path: PathBuf,
        /// Original binary decoder failure.
        #[source]
        source: DecodeError,
    },
    /// Saved process settings history could not be decoded as TOML.
    #[error("Could not decode process settings history {path}")]
    SettingsDecode {
        /// Settings-history file whose contents could not be decoded.
        path: PathBuf,
        /// Original TOML decoder failure.
        #[source]
        source: SettingsDecodeError,
    },
    /// A process payload operation failed, retaining its concrete nested error.
    #[error(transparent)]
    Payload(Box<ProcessPayloadFailure>),
}

impl ProcessSnapshotError {
    pub(super) fn payload(
        operation: &'static str,
        process: &str,
        path: &Path,
        source: ProcessPayloadError,
    ) -> Self {
        Self::Payload(Box::new(ProcessPayloadFailure {
            operation,
            process: process.to_owned(),
            path: path.to_owned(),
            source,
        }))
    }
}

/// A payload failure with the owning process and operation context.
///
/// [`ProcessSnapshotError::Payload`] boxes this context and transparently exposes
/// its native [`ProcessPayloadError`] source. Boxing the context keeps the outer
/// enum small without hiding the concrete source behind a boxed error type.
///
/// ```
/// use gammalooprs::processes::{
///     ProcessPayloadError, ProcessPayloadFailure, ProcessSnapshotError,
/// };
/// use std::error::Error;
/// use std::io::{Error as IoError, ErrorKind};
/// let error = ProcessSnapshotError::Payload(Box::new(ProcessPayloadFailure {
///     operation: "read payload",
///     process: "example".into(),
///     path: "example/LO".into(),
///     source: ProcessPayloadError::Read {
///         path: "example/LO/amp.bin".into(),
///         source: IoError::from(ErrorKind::NotFound),
///     },
/// }));
/// assert!(error.source().unwrap().is::<ProcessPayloadError>());
/// ```
#[derive(Debug, Error)]
#[error("Could not {operation} for process {process:?} at {path}")]
pub struct ProcessPayloadFailure {
    /// Action being performed when the payload operation failed.
    pub operation: &'static str,
    /// Folder name of the process owning the payload.
    pub process: String,
    /// Directory at the process or child operation boundary.
    pub path: PathBuf,
    /// Native payload failure, including its more specific affected path.
    #[source]
    pub source: ProcessPayloadError,
}

/// Native failures while reading or writing a process payload.
///
/// Each failure retains the exact payload path and original source. These are
/// exposed through [`ProcessSnapshotError::Payload`].
///
/// ```
/// use gammalooprs::processes::ProcessPayloadError;
/// use std::error::Error;
/// use std::io::{Error as IoError, ErrorKind};
/// let error = ProcessPayloadError::Read {
///     path: "integrand/settings.toml".into(),
///     source: IoError::from(ErrorKind::NotFound),
/// };
/// let source = error.source().unwrap().downcast_ref::<IoError>().unwrap();
/// assert_eq!(source.kind(), ErrorKind::NotFound);
/// ```
#[derive(Debug, Error)]
pub enum ProcessPayloadError {
    /// A managed payload path was rejected or could not be inspected or created.
    #[error(transparent)]
    Output(#[from] OutputPathError),
    /// A binary payload or pinned runtime-settings file could not be read.
    #[error("Could not read process payload {path}")]
    Read {
        /// Payload file whose read failed.
        path: PathBuf,
        /// Original filesystem or text-decoding failure.
        #[source]
        source: IoError,
    },
    /// Creating, writing, or flushing a payload entry failed.
    #[error("Could not write process payload {path}")]
    Write {
        /// Payload entry whose creation, write, or flush failed.
        path: PathBuf,
        /// Original filesystem failure.
        #[source]
        source: IoError,
    },
    /// Child metadata or generated integrand data could not be decoded from binary.
    #[error("Could not decode process payload {path}")]
    Decode {
        /// Binary payload file whose bytes could not be decoded.
        path: PathBuf,
        /// Original binary decoder failure.
        #[source]
        source: DecodeError,
    },
    /// A process definition, child metadata, or generated integrand could not be encoded.
    #[error("Could not encode process payload {path}")]
    Encode {
        /// Intended destination of the binary payload.
        path: PathBuf,
        /// Original binary encoder failure.
        #[source]
        source: EncodeError,
    },
    /// A generated integrand's pinned runtime settings could not be decoded as TOML.
    #[error("Could not decode process settings {path}")]
    SettingsDecode {
        /// Pinned settings file whose contents could not be decoded.
        path: PathBuf,
        /// Original TOML decoder failure.
        #[source]
        source: SettingsDecodeError,
    },
    /// Process settings history or pinned runtime settings could not be encoded as TOML.
    #[error("Could not encode process settings {path}")]
    SettingsEncode {
        /// Intended destination of the serialized settings.
        path: PathBuf,
        /// Original TOML encoder failure.
        #[source]
        source: SettingsEncodeError,
    },
    /// Child metadata embeds an integrand, which must be stored in a separate payload.
    #[error("Saved child metadata {path} unexpectedly contains an integrand")]
    MetadataIntegrand {
        /// Child metadata file containing the unexpected embedded integrand.
        path: PathBuf,
    },
    /// A decoded generated integrand names a different child than its metadata does.
    #[error("Generated payload {path} has name {actual:?}, expected {expected:?}")]
    GeneratedNameMismatch {
        /// Directory containing the generated integrand payload.
        path: PathBuf,
        /// Child name required by the decoded metadata.
        expected: String,
        /// Name recorded in the generated integrand data.
        actual: String,
    },
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error as StdError;
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
    fn process_snapshot_payload_exposes_write_io() {
        let payload = ProcessPayloadError::Write {
            path: "stage/amp.bin".into(),
            source: IoError::from(ErrorKind::PermissionDenied),
        };
        let error =
            ProcessSnapshotError::payload("write payload", "process", Path::new("stage"), payload);
        assert_eq!(
            native_source::<IoError>(&error).kind(),
            ErrorKind::PermissionDenied
        );
        assert!(matches!(
            error.source().unwrap().downcast_ref::<ProcessPayloadError>(),
            Some(ProcessPayloadError::Write { path, .. }) if path == Path::new("stage/amp.bin")
        ));
    }

    #[test]
    fn process_snapshot_payload_exposes_read_io() {
        let payload = ProcessPayloadError::Read {
            path: "stage/amp.bin".into(),
            source: IoError::from(ErrorKind::NotFound),
        };
        let error =
            ProcessSnapshotError::payload("load payload", "process", Path::new("stage"), payload);
        assert_eq!(native_source::<IoError>(&error).kind(), ErrorKind::NotFound);
    }

    #[test]
    fn process_snapshot_payload_exposes_decode() {
        let payload = ProcessPayloadError::Decode {
            path: "stage/amp.bin".into(),
            source: DecodeError::UnexpectedEnd { additional: 5 },
        };
        let error =
            ProcessSnapshotError::payload("load payload", "process", Path::new("stage"), payload);
        assert!(matches!(
            native_source::<DecodeError>(&error),
            DecodeError::UnexpectedEnd { additional: 5 }
        ));
    }
}
