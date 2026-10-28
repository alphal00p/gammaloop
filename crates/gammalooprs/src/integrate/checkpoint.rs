//! One manifest commits the scalar/grid state, settings, and observables of a completed iteration.
//!
//! Readers select an immutable generation; writers retain their OS lock throughout a run.
//! Ordinary commits retain every generation, so readers remain pinned across later publications.
//! Explicit cleanup and a successful restart may discard old generations; coordinate those
//! operations with readers. Failed restart preparation preserves the old run. See
//! [`IntegrationWorkspace`].

use super::{
    INTEGRANDS_DIRECTORY, INTEGRATION_RESULT_FILE, INTERRUPTED_DIRECTORY, IntegrationSlot,
    IntegrationState, IntegrationWorkspaceManifest, NUMERICAL_STABILITY_DIRECTORY,
    PRESENTATION_TEMP_PREFIX, PRESENTATION_TEMP_SUFFIX, RESULTS_DIRECTORY, SlotMeta,
};
use crate::observables::ObservableSnapshotBundle;
use crate::settings::RuntimeSettings;
use crate::utils::F;
use crate::{INTEGRATION_MANIFEST_FILE, STATE_MANIFEST_FILE, integrate};
use alphal00p_persistence::{
    GENERATIONS_DIRECTORY, GenerationId, ManifestName, OutputPath, OutputPathError, SnapshotError,
    SnapshotWriter,
};
use bincode::config;
use bincode::error::{DecodeError, EncodeError};
use serde::{Deserialize, Serialize};
use serde_json::Error as JsonError;
use spenso::algebra::complex::Complex;
use std::collections::HashSet;
use std::fmt::{Display, Formatter, Result as FmtResult};
use std::fs;
use std::fs::File;
use std::io::{Error as IoError, ErrorKind, Write};
use std::path::{Path, PathBuf};
use std::string::FromUtf8Error;
use thiserror::Error;
use toml::de::Error as SettingsDecodeError;
use toml::ser::Error as SettingsEncodeError;
use tracing::warn;

#[cfg(test)]
use std::env;
#[cfg(test)]
use std::process;
#[cfg(test)]
use std::sync::atomic::{AtomicU64, Ordering};

#[cfg(test)]
pub(super) struct CheckpointTestDirectory(PathBuf);

#[cfg(test)]
impl CheckpointTestDirectory {
    pub(super) fn new() -> Self {
        static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(0);
        for _ in 0..128 {
            let sequence = NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed);
            let path = env::temp_dir().join(format!(
                "gammaloop-checkpoint-{}-{sequence:x}",
                process::id()
            ));
            match fs::create_dir(&path) {
                Ok(()) => return Self(path),
                Err(source) if source.kind() == ErrorKind::AlreadyExists => continue,
                Err(source) => panic!("cannot create checkpoint test directory: {source}"),
            }
        }
        panic!("checkpoint test directory collision limit exceeded");
    }

    pub(super) fn path(&self) -> &Path {
        &self.0
    }
}

#[cfg(test)]
impl Drop for CheckpointTestDirectory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Failures before publication leave the last committed iteration authoritative.
/// See the [`IntegrationWorkspace`] example for propagation at the API boundary.
#[derive(Debug, Error)]
pub enum IntegrationCheckpointError {
    /// Lock acquisition, generation publication, or pruning failed in the shared owner.
    #[error(transparent)]
    Snapshot(#[from] SnapshotError),
    /// A managed path violates the shared filesystem policy.
    #[error(transparent)]
    Output(#[from] OutputPathError),
    /// A checkpoint file could not be read, written, renamed, or synchronized.
    #[error("Integration checkpoint I/O failed at {path}")]
    Io {
        /// The affected member or journal entry.
        path: PathBuf,
        /// The original filesystem failure.
        source: IoError,
    },
    /// A payload referenced by the committed manifest is missing.
    #[error("Integration checkpoint member is missing at {path}")]
    MissingMember {
        /// The absent scalar, settings, or observable payload.
        path: PathBuf,
        /// The original missing-file failure.
        source: IoError,
    },
    /// A JSON member or restart journal could not be encoded or decoded.
    #[error("Invalid integration checkpoint JSON at {path}")]
    Json {
        /// The affected JSON member.
        path: PathBuf,
        /// The original serialization failure.
        source: JsonError,
    },
    /// A recognizable flat manifest predates immutable checkpoint generations.
    #[error("Integration workspace at {path} was created by an older GammaLoop version")]
    LegacyLayout {
        /// The workspace root, which this version will not migrate or overwrite.
        path: PathBuf,
    },
    /// The manifest could not be decoded into a checkpoint envelope.
    #[error("Invalid integration checkpoint manifest at {path}")]
    ManifestDecode {
        /// The manifest that failed validation before any state decoding.
        path: PathBuf,
        /// The original JSON decoding failure.
        source: JsonError,
    },
    /// A slot's pinned TOML settings could not be decoded.
    #[error("Invalid integration checkpoint settings at {path}")]
    SettingsDecode {
        /// The slot settings member.
        path: PathBuf,
        /// The original TOML decoding failure.
        source: SettingsDecodeError,
    },
    /// A slot's pinned settings are not UTF-8 text.
    #[error("Invalid integration checkpoint settings encoding at {path}")]
    SettingsEncoding {
        /// The settings member containing invalid text.
        path: PathBuf,
        /// The original text-decoding failure.
        source: FromUtf8Error,
    },
    /// Candidate slot settings could not be serialized before publication.
    #[error("Could not encode integration checkpoint settings at {path}")]
    SettingsEncode {
        /// The candidate settings member.
        path: PathBuf,
        /// The original TOML encoding failure.
        source: SettingsEncodeError,
    },
    /// The pinned scalar/grid payload could not be decoded.
    #[error("Could not decode integration checkpoint state at {path}")]
    StateDecode {
        /// The scalar/grid member.
        path: PathBuf,
        /// The original binary decoding failure.
        source: DecodeError,
    },
    /// Candidate scalar/grid state could not be serialized before publication.
    #[error("Could not encode integration checkpoint state at {path}")]
    StateEncode {
        /// The candidate scalar/grid member.
        path: PathBuf,
        /// The original binary encoding failure.
        source: EncodeError,
    },
    /// The stored version is absent or differs from the implemented format.
    #[error(
        "Integration workspace version {found} is unsupported (expected {expected}; \
         zero denotes an unversioned workspace)"
    )]
    UnsupportedVersion {
        /// The stored version, or zero for an absent version tag.
        found: u32,
        /// The version implemented by this build.
        expected: u32,
    },
    /// Slot identities are empty or repeated, or per-slot configuration lengths disagree.
    #[error("Workspace slot identities or configuration lengths disagree ({slots} slots)")]
    WorkspaceSlots {
        /// The declared slot count.
        slots: usize,
    },
    /// Zero iterations or zero samples cannot form an authoritative checkpoint.
    #[error(
        "Checkpoint must describe a completed iteration \
         (iteration {iteration}, samples {num_points})"
    )]
    IncompleteIteration {
        /// The candidate completed-iteration count.
        iteration: usize,
        /// The candidate completed-sample count.
        num_points: usize,
    },
    /// Observable membership contains an out-of-range, repeated, or unsorted slot index.
    #[error("Checkpoint observable indices must be sorted, unique, and within {slots} slots")]
    ObservableMembership {
        /// The number of slots available to the membership list.
        slots: usize,
    },
    /// Scalar, grid, or accumulator identity differs from its committed envelope.
    #[error("Scalar/grid state at {path} disagrees with its completed-iteration descriptor")]
    StateDescriptor {
        /// The inconsistent scalar/grid member.
        path: PathBuf,
    },
    /// Binary decoding did not consume the complete scalar/grid member.
    #[error("Scalar/grid state at {path} contains {remaining} trailing bytes")]
    TrailingState {
        /// The member containing trailing bytes.
        path: PathBuf,
        /// The number of unconsumed bytes.
        remaining: usize,
    },
    /// A caller requested a slot outside the ordered manifest slots.
    #[error("Checkpoint slot index {slot} is outside its {slots} slots")]
    SlotOutOfBounds {
        /// The requested slot.
        slot: SlotIndex,
        /// The number of available slots.
        slots: usize,
    },
    /// A settings or observable envelope names a different slot, generation, or iteration.
    #[error(
        "Slot {slot} member at {path} disagrees with its \
         generation/iteration/configuration descriptor"
    )]
    MemberDescriptor {
        /// The expected ordered slot.
        slot: SlotIndex,
        /// The mismatched member.
        path: PathBuf,
    },
    /// An observable payload exists despite the manifest declaring its absence.
    #[error("Undeclared observable bundle for slot {slot} at {path}")]
    UndeclaredObservables {
        /// The undeclared member's slot.
        slot: SlotIndex,
        /// The unexpected observable file.
        path: PathBuf,
    },
    /// Histogram membership or bin layout differs from pinned slot settings.
    #[error("Slot {slot} observable configuration or membership disagrees with its settings")]
    ObservableConfiguration {
        /// The inconsistent histogram owner.
        slot: SlotIndex,
    },
    /// An expected payload is a directory or another non-file entry.
    #[error("Checkpoint member at {path} is not a regular file")]
    NotRegularFile {
        /// The rejected member.
        path: PathBuf,
    },
    /// Worker targets or their count differ from the workspace's requested targets.
    #[error("Integration slot targets disagree with the workspace configuration")]
    SlotTargets,
    /// Worker slot identity, bundle count, or fingerprint count differs from the workspace.
    #[error("Integration slots disagree with the workspace configuration")]
    CommitSlots,
    /// A normal commit attempts to replay an iteration instead of advancing the run.
    #[error(
        "Completed iteration {iteration} does not advance previous iteration {previous}; \
         prepare a restart before replacing a run"
    )]
    IterationOrder {
        /// The latest committed iteration.
        previous: usize,
        /// The proposed iteration.
        iteration: usize,
    },
    /// Outputs from an earlier attempt exist without a readable committed checkpoint.
    #[error("Integration outputs at {path} exist without a committed checkpoint")]
    UncommittedOutputs {
        /// The first existing output that prevents a fresh run.
        path: PathBuf,
    },
    /// Integration ownership cannot overlap a saved application's manifest and generations.
    #[error("A saved application state at {path} cannot be used as an integration workspace")]
    SavedStateWorkspace {
        /// The competing application-state manifest.
        path: PathBuf,
    },
    /// Run configuration was changed after its first checkpoint was published.
    #[error("Cannot change workspace configuration after a checkpoint was committed")]
    ConfigurationCommitted,
    /// Reference descriptors do not provide one entry for each ordered slot.
    #[error("Workspace reference descriptors must match its {slots} slots")]
    ReferenceSlots {
        /// The required descriptor count.
        slots: usize,
    },
    /// A reference has no resolved center or has an invalid normalization scale.
    #[error(
        "Workspace reference descriptor for slot {slot} must contain a resolved center \
         and finite positive moment scale"
    )]
    ReferenceDescriptor {
        /// The invalid reference's owner.
        slot: SlotIndex,
    },
    /// Requested reference parameters or the observable convention differ from the saved run.
    #[error(
        "Workspace physical/reference configuration or reporting convention \
         does not match the request"
    )]
    ReferenceConfiguration,
    /// The private restart journal cannot authorize an unambiguous output transition.
    #[error("Restart output journal at {path} does not match its workspace")]
    RestartJournal {
        /// The inconsistent journal entry or conflicting output.
        path: PathBuf,
    },
    /// A postcommit failure requires reopening the writer to finish its recovery journal.
    #[error("Reopen integration workspace {path} before writing after a postcommit failure")]
    RecoveryRequired {
        /// The workspace whose last write did not complete cleanup.
        path: PathBuf,
    },
    /// Publication succeeded; resume must advance from this iteration despite cleanup failure.
    #[error(
        "Iteration {iteration} checkpoint committed, \
         but output cleanup or generation pruning failed"
    )]
    PostCommit {
        /// The iteration that is already authoritative.
        iteration: usize,
        /// The original output-transition or pruning failure.
        source: Box<IntegrationCheckpointError>,
    },
}

/// Recovery choices for unreadable checkpoints or incompatible integration requests.
/// See [`IntegrationCheckpointError::recovery_action`] for classification without losing sources.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum IntegrationCheckpointRecovery {
    /// Select a fresh workspace, or explicitly remove the unusable one.
    FreshWorkspace,
    /// Keep application state separate from integration output.
    DedicatedWorkspace,
    /// Restore the saved reference configuration or explicitly restart it.
    RestoreReference,
}

impl IntegrationCheckpointError {
    /// Whether another writer holds the snapshot lock, allowing a retry hint at the CLI boundary.
    ///
    /// ```
    /// # use alphal00p_persistence::SnapshotError;
    /// # use gammalooprs::integrate::IntegrationCheckpointError;
    /// let error = SnapshotError::Busy { path: "run/.snapshot.lock".into() };
    /// assert!(IntegrationCheckpointError::from(error).is_busy());
    /// ```
    #[must_use]
    pub fn is_busy(&self) -> bool {
        matches!(self, Self::Snapshot(source) if source.is_busy())
    }

    /// Classify actionable recovery without replacing the typed error or its source.
    /// I/O and publication failures have no generic destructive recovery suggestion.
    ///
    /// ```
    /// use gammalooprs::integrate::{IntegrationCheckpointError, IntegrationCheckpointRecovery};
    /// let error = IntegrationCheckpointError::UnsupportedVersion { found: 3, expected: 4 };
    /// assert_eq!(error.recovery_action(), Some(IntegrationCheckpointRecovery::FreshWorkspace));
    /// ```
    #[must_use]
    pub fn recovery_action(&self) -> Option<IntegrationCheckpointRecovery> {
        match self {
            Self::LegacyLayout { .. }
            | Self::UnsupportedVersion { .. }
            | Self::ManifestDecode { .. }
            | Self::MissingMember { .. }
            | Self::SettingsDecode { .. }
            | Self::SettingsEncoding { .. }
            | Self::StateDecode { .. }
            | Self::Json { .. }
            | Self::WorkspaceSlots { .. }
            | Self::IncompleteIteration { .. }
            | Self::ObservableMembership { .. }
            | Self::StateDescriptor { .. }
            | Self::TrailingState { .. }
            | Self::MemberDescriptor { .. }
            | Self::UndeclaredObservables { .. }
            | Self::ObservableConfiguration { .. }
            | Self::NotRegularFile { .. }
            | Self::ReferenceSlots { .. }
            | Self::ReferenceDescriptor { .. }
            | Self::UncommittedOutputs { .. } => {
                Some(IntegrationCheckpointRecovery::FreshWorkspace)
            }
            Self::SavedStateWorkspace { .. } => {
                Some(IntegrationCheckpointRecovery::DedicatedWorkspace)
            }
            Self::ReferenceConfiguration => Some(IntegrationCheckpointRecovery::RestoreReference),
            Self::Snapshot(_)
            | Self::Output(_)
            | Self::Io { .. }
            | Self::SettingsEncode { .. }
            | Self::StateEncode { .. }
            | Self::SlotOutOfBounds { .. }
            | Self::SlotTargets
            | Self::CommitSlots
            | Self::IterationOrder { .. }
            | Self::ConfigurationCommitted
            | Self::RestartJournal { .. }
            | Self::RecoveryRequired { .. }
            | Self::PostCommit { .. } => None,
        }
    }
}

type CheckpointResult<T> = Result<T, IntegrationCheckpointError>;

const STATE_FILE: &str = "integration_state.bin";
const SLOT_SETTINGS_FILE: &str = "settings.toml";
const OBSERVABLES_FILE: &str = "observables.json";
const SLOTS_DIRECTORY: &str = "slots";

/// An index in the ordered integration-slot domain, distinct from graph and sampling indices.
/// ```
/// use gammalooprs::integrate::SlotIndex;
/// let slot = SlotIndex::from(2);
/// assert_eq!(slot.index(), 2);
/// ```
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[serde(transparent)]
pub struct SlotIndex(usize);

impl SlotIndex {
    /// The underlying position, as in the [type example](Self).
    #[must_use]
    pub fn index(self) -> usize {
        self.0
    }
}

impl From<usize> for SlotIndex {
    fn from(index: usize) -> Self {
        Self(index)
    }
}

impl Display for SlotIndex {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> FmtResult {
        self.0.fmt(formatter)
    }
}

/// Version 4 commit envelope. The ordered workspace slots identify every numeric slot directory.
/// Unlike a requested workspace configuration, an envelope always names a real generation.
/// Read it through [`IntegrationCheckpoint::open`], as in the [`IntegrationWorkspace`] example.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct IntegrationCheckpointManifest {
    /// Ordered slot identities and the run's validated configuration.
    pub workspace: IntegrationWorkspaceManifest,
    /// The immutable directory selected by this manifest.
    pub generation: GenerationId,
    /// The number of fully accumulated iterations.
    pub completed_iteration: usize,
    /// The number of Monte Carlo points accumulated across those iterations.
    pub num_points: usize,
    /// Sorted, unique slots that own observable payloads in this generation.
    pub observable_slots: Vec<SlotIndex>,
}

/// Decoded payloads from one validated generation, in manifest slot order.
/// See [`IntegrationCheckpoint::load`] and the [`IntegrationWorkspace`] example.
#[derive(Clone, Debug)]
pub struct LoadedIntegrationCheckpoint {
    /// Scalar accumulators, trained grids, and statistics at the completed iteration boundary.
    pub state: IntegrationState,
    /// Pinned per-slot payloads in the manifest's order.
    pub slots: Vec<LoadedIntegrationSlot>,
}

/// Settings and optional observable payload for one ordered slot.
/// See [`LoadedIntegrationCheckpoint`].
#[derive(Clone, Debug, PartialEq)]
pub struct LoadedIntegrationSlot {
    /// The settings that produced the checkpoint's scalar and observable state.
    pub settings: RuntimeSettings,
    /// The declared observable payload, or explicit absence when the slot has no observables.
    pub observables: Option<ObservableSnapshotBundle>,
}

#[derive(Debug, Serialize, Deserialize)]
struct SlotCheckpointMember<Generation, Slot, Payload> {
    generation: Generation,
    completed_iteration: usize,
    num_points: usize,
    slot: Slot,
    payload: Payload,
}

/// A read-only, validated pointer to one immutable completed iteration.
/// Opening never creates directories or lockfiles. See the [`IntegrationWorkspace`] example.
#[derive(Debug)]
pub struct IntegrationCheckpoint {
    manifest: IntegrationCheckpointManifest,
    payload_root: OutputPath,
    root: OutputPath,
}

impl IntegrationCheckpoint {
    /// Pin and validate all referenced members before returning a usable checkpoint.
    ///
    /// # Errors
    /// Missing members, unsupported versions, invalid descriptors, and malformed payloads are
    /// corruption errors; only an absent *manifest* can denote a fresh workspace to a writer.
    /// See the [`IntegrationWorkspace`] example.
    pub fn open(root: &Path) -> CheckpointResult<Self> {
        let root = OutputPath::new(root)?;
        Self::open_checked(&root)
    }

    fn open_checked(root: &OutputPath) -> CheckpointResult<Self> {
        let manifest_path = root.join(INTEGRATION_MANIFEST_FILE)?;
        let bytes = Self::read(&manifest_path)?;
        Self::from_manifest_bytes(root, &bytes)
    }

    fn from_manifest_bytes(root: &OutputPath, bytes: &[u8]) -> CheckpointResult<Self> {
        let checkpoint = Self::decode_manifest(root, bytes)?;
        checkpoint.load()?;
        Ok(checkpoint)
    }

    // Commit ordering needs only the validated envelope, not all grid/observable payloads.
    fn decode_manifest(root: &OutputPath, bytes: &[u8]) -> CheckpointResult<Self> {
        let manifest_path = integrate::workspace_manifest_path(root.as_path());
        let manifest: IntegrationCheckpointManifest =
            serde_json::from_slice(bytes).map_err(|source| {
                // Recognize the old flat layout only to explain refusal; never load or migrate it.
                if serde_json::from_slice::<IntegrationWorkspaceManifest>(bytes).is_ok() {
                    IntegrationCheckpointError::LegacyLayout {
                        path: root.as_path().to_owned(),
                    }
                } else {
                    IntegrationCheckpointError::ManifestDecode {
                        path: manifest_path,
                        source,
                    }
                }
            })?;
        Self::validate_workspace(&manifest.workspace)?;
        let n_slots = manifest.workspace.slots.len();
        if manifest.completed_iteration == 0 || manifest.num_points == 0 {
            return Err(IntegrationCheckpointError::IncompleteIteration {
                iteration: manifest.completed_iteration,
                num_points: manifest.num_points,
            });
        }
        if manifest
            .observable_slots
            .iter()
            .any(|slot| slot.index() >= n_slots)
            || manifest
                .observable_slots
                .windows(2)
                .any(|pair| pair[0] >= pair[1])
        {
            return Err(IntegrationCheckpointError::ObservableMembership { slots: n_slots });
        }
        let payload_root = root
            .join(GENERATIONS_DIRECTORY)?
            .join(manifest.generation.as_str())?;
        Ok(Self {
            payload_root,
            root: root.clone(),
            manifest,
        })
    }

    fn validate_workspace(workspace: &IntegrationWorkspaceManifest) -> CheckpointResult<()> {
        workspace.validate_version()?;
        let n_slots = workspace.slots.len();
        let unique_slots: HashSet<&SlotMeta> = workspace.slots.iter().collect();
        let unique_keys: HashSet<String> = workspace.slots.iter().map(SlotMeta::key).collect();
        for slot in &workspace.slots {
            OutputPath::validate_component(&slot.process_name)?;
            OutputPath::validate_component(&slot.integrand_name)?;
        }
        if n_slots == 0
            || unique_slots.len() != n_slots
            || unique_keys.len() != n_slots
            || workspace.targets.len() != n_slots
            || workspace.effective_model_parameters.len() != n_slots
            || workspace.integrand_fingerprints.len() != n_slots
            || workspace.training_slot.index() >= n_slots
            || workspace.integrator_settings_slot.index() >= n_slots
            || (!workspace.reference_gaussians.is_empty()
                && workspace.reference_gaussians.len() != n_slots)
        {
            return Err(IntegrationCheckpointError::WorkspaceSlots { slots: n_slots });
        }
        Ok(())
    }

    /// The descriptor pinned when this reader opened. See [`IntegrationWorkspace`].
    #[must_use]
    pub fn manifest(&self) -> &IntegrationCheckpointManifest {
        &self.manifest
    }

    /// Whether a committed restart still needs output cleanup under the next writer lock.
    /// Read-only summary callers can warn without changing the checkpoint or journal.
    ///
    /// ```no_run
    /// # use gammalooprs::integrate::IntegrationCheckpoint;
    /// # use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// let checkpoint = IntegrationCheckpoint::open("run".as_ref())?;
    /// let pending = checkpoint.pending_restart_cleanup()?;
    /// # let _ = pending;
    /// # Ok(()) }
    /// ```
    ///
    /// # Errors
    /// Returns an error for unsafe journal paths, malformed journal JSON, or I/O failures.
    pub fn pending_restart_cleanup(&self) -> CheckpointResult<bool> {
        let restart = RestartOutputs::new(&self.root)?;
        Ok(restart.committed_generation()?.as_ref() == Some(&self.manifest.generation))
    }

    /// Scalar moments and trained grids in this generation. See [`IntegrationWorkspace`].
    #[must_use]
    pub fn state_path(&self) -> PathBuf {
        self.payload_root.as_path().join(STATE_FILE)
    }

    /// Settings for the ordered slot index. See [`IntegrationWorkspace`].
    #[must_use]
    pub fn slot_settings_path(&self, slot: SlotIndex) -> PathBuf {
        self.payload_root
            .as_path()
            .join(SLOTS_DIRECTORY)
            .join(slot.to_string())
            .join(SLOT_SETTINGS_FILE)
    }

    /// A path only for declared observable members. See [`IntegrationWorkspace`].
    #[must_use]
    pub fn observable_path(&self, slot: SlotIndex) -> Option<PathBuf> {
        self.manifest
            .observable_slots
            .binary_search(&slot)
            .ok()
            .map(|_| {
                self.payload_root
                    .as_path()
                    .join(SLOTS_DIRECTORY)
                    .join(slot.to_string())
                    .join(OBSERVABLES_FILE)
            })
    }

    /// Decode and check scalar/grid identity against the pinned envelope.
    ///
    /// # Errors
    /// Returns an error for missing, malformed, partial, or inconsistent state. See
    /// [`IntegrationWorkspace`] for an example.
    pub fn read_state(&self) -> CheckpointResult<IntegrationState> {
        let checked = self.payload_root.join(STATE_FILE)?;
        let path = checked.as_path().to_owned();
        let bytes = Self::read_member(&checked)?;
        let (state, consumed): (IntegrationState, usize) =
            bincode::decode_from_slice(&bytes, config::standard()).map_err(|source| {
                IntegrationCheckpointError::StateDecode {
                    path: path.clone(),
                    source,
                }
            })?;
        self.validate_state(&state)?;
        if consumed != bytes.len() {
            return Err(IntegrationCheckpointError::TrailingState {
                path,
                remaining: bytes.len() - consumed,
            });
        }
        Ok(state)
    }

    /// Load a coherent scalar/settings/observable result from this generation.
    ///
    /// # Errors
    /// Missing, malformed, or inconsistent members fail with their path or slot identity.
    /// Explicit cleanup or a successful restart can remove this pinned generation.
    /// See the [`IntegrationWorkspace`] example.
    pub fn load(&self) -> CheckpointResult<LoadedIntegrationCheckpoint> {
        let state = self.read_state()?;
        let slots = (0..self.manifest.workspace.slots.len())
            .map(|index| {
                let slot = SlotIndex::from(index);
                let settings = self.read_slot_settings(slot)?;
                let observables = self.read_observables(slot)?;
                Self::validate_observables(slot, &settings, observables.as_ref())?;
                Ok(LoadedIntegrationSlot {
                    settings,
                    observables,
                })
            })
            .collect::<CheckpointResult<_>>()?;
        Ok(LoadedIntegrationCheckpoint { state, slots })
    }

    fn validate_state(&self, state: &IntegrationState) -> CheckpointResult<()> {
        let workspace = &self.manifest.workspace;
        let n_slots = workspace.slots.len();
        if state.iter != self.manifest.completed_iteration
            || state.num_points != self.manifest.num_points
            || state.slot_metas != workspace.slots
            || state.sampling_correlation_mode != workspace.sampling_correlation_mode
            || state.sampling_states.len()
                != workspace.sampling_correlation_mode.state_count(n_slots)
            || state.all_integrals.len() != n_slots
            || state.slot_stats.len() != n_slots
            || state.slot_re_summaries.len() != n_slots
            || state.slot_im_summaries.len() != n_slots
            || state.slot_absolute_re_summaries.len() != n_slots
            || state.slot_absolute_im_summaries.len() != n_slots
            || state
                .slot_first_non_trivial_discrete_breakdown_metadata
                .len()
                != n_slots
            || state.all_integrals.iter().any(|acc| {
                [&acc.re, &acc.im, &acc.absolute_re, &acc.absolute_im]
                    .iter()
                    .any(|component| {
                        component.processed_samples > state.num_points
                            || component.cur_iter != state.iter
                    })
            })
        {
            return Err(IntegrationCheckpointError::StateDescriptor {
                path: self.state_path(),
            });
        }
        Ok(())
    }

    /// Read one member's settings without resolving a newer manifest.
    ///
    /// # Errors
    /// Invalid indices, missing files, and malformed settings fail. See [`IntegrationWorkspace`].
    pub fn read_slot_settings(&self, slot: SlotIndex) -> CheckpointResult<RuntimeSettings> {
        self.validate_slot(slot)?;
        let checked = self.slot_path(slot)?.join(SLOT_SETTINGS_FILE)?;
        let path = checked.as_path().to_owned();
        let bytes = Self::read_member(&checked)?;
        let text = String::from_utf8(bytes).map_err(|source| {
            IntegrationCheckpointError::SettingsEncoding {
                path: path.clone(),
                source,
            }
        })?;
        let member: SlotCheckpointMember<GenerationId, SlotMeta, RuntimeSettings> =
            toml::from_str(&text).map_err(|source| IntegrationCheckpointError::SettingsDecode {
                path: path.clone(),
                source,
            })?;
        self.validate_member(slot, &path, &member)?;
        Ok(member.payload)
    }

    /// Read exactly the bundle declared for this slot; absence is explicit in the envelope.
    ///
    /// # Errors
    /// Missing declared members, undeclared bundles, or invalid indices fail. See
    /// [`IntegrationWorkspace`] for an example.
    pub fn read_observables(
        &self,
        slot: SlotIndex,
    ) -> CheckpointResult<Option<ObservableSnapshotBundle>> {
        self.validate_slot(slot)?;
        let slot_path = self.slot_path(slot)?;
        if let Some(path) = self.observable_path(slot) {
            let checked = slot_path.join(OBSERVABLES_FILE)?;
            let bytes = Self::read_member(&checked)?;
            let member: SlotCheckpointMember<GenerationId, SlotMeta, ObservableSnapshotBundle> =
                serde_json::from_slice(&bytes).map_err(|source| {
                    IntegrationCheckpointError::Json {
                        path: path.clone(),
                        source,
                    }
                })?;
            self.validate_member(slot, &path, &member)?;
            return Ok(Some(member.payload));
        }
        let path = slot_path.join(OBSERVABLES_FILE)?.as_path().to_owned();
        match fs::symlink_metadata(&path) {
            Err(source) if source.kind() == ErrorKind::NotFound => Ok(None),
            Err(source) => Err(IntegrationCheckpointError::Io { path, source }),
            Ok(_) => Err(IntegrationCheckpointError::UndeclaredObservables { slot, path }),
        }
    }

    fn validate_slot(&self, slot: SlotIndex) -> CheckpointResult<()> {
        let slots = self.manifest.workspace.slots.len();
        if slot.index() >= slots {
            return Err(IntegrationCheckpointError::SlotOutOfBounds { slot, slots });
        }
        Ok(())
    }

    fn slot_path(&self, slot: SlotIndex) -> CheckpointResult<OutputPath> {
        Ok(self
            .payload_root
            .join(SLOTS_DIRECTORY)?
            .join(&slot.to_string())?)
    }

    fn validate_member<T>(
        &self,
        slot: SlotIndex,
        path: &Path,
        member: &SlotCheckpointMember<GenerationId, SlotMeta, T>,
    ) -> CheckpointResult<()> {
        if member.generation != self.manifest.generation
            || member.completed_iteration != self.manifest.completed_iteration
            || member.num_points != self.manifest.num_points
            || member.slot != self.manifest.workspace.slots[slot.index()]
        {
            return Err(IntegrationCheckpointError::MemberDescriptor {
                slot,
                path: path.to_owned(),
            });
        }
        Ok(())
    }

    fn member<'a, T>(
        &'a self,
        slot: SlotIndex,
        payload: &'a T,
    ) -> SlotCheckpointMember<&'a GenerationId, &'a SlotMeta, &'a T> {
        SlotCheckpointMember {
            generation: &self.manifest.generation,
            completed_iteration: self.manifest.completed_iteration,
            num_points: self.manifest.num_points,
            slot: &self.manifest.workspace.slots[slot.index()],
            payload,
        }
    }

    fn validate_observables(
        slot: SlotIndex,
        settings: &RuntimeSettings,
        bundle: Option<&ObservableSnapshotBundle>,
    ) -> CheckpointResult<()> {
        match bundle {
            Some(bundle)
                if !bundle.histograms.is_empty()
                    && bundle.matches_checkpoint_settings(&settings.observables) =>
            {
                Ok(())
            }
            None if settings.observables.is_empty() => Ok(()),
            _ => Err(IntegrationCheckpointError::ObservableConfiguration { slot }),
        }
    }

    // A missing declared payload is corruption. An absent manifest or unrelated I/O failure
    // retains its original classification, so callers do not mistake fresh work for corruption.
    fn read_member(checked: &OutputPath) -> CheckpointResult<Vec<u8>> {
        Self::read(checked).map_err(|error| match error {
            IntegrationCheckpointError::Io { path, source }
                if source.kind() == ErrorKind::NotFound =>
            {
                IntegrationCheckpointError::MissingMember { path, source }
            }
            error => error,
        })
    }

    fn read(checked: &OutputPath) -> CheckpointResult<Vec<u8>> {
        checked.check()?;
        let path = checked.as_path();
        let metadata =
            fs::symlink_metadata(path).map_err(|source| IntegrationCheckpointError::Io {
                path: path.to_owned(),
                source,
            })?;
        if !metadata.is_file() {
            return Err(IntegrationCheckpointError::NotRegularFile {
                path: path.to_owned(),
            });
        }
        fs::read(path).map_err(|source| IntegrationCheckpointError::Io {
            path: path.to_owned(),
            source,
        })
    }
}

/// Owns exclusive writing for a run, including resume, restart preparation, and output assembly.
/// Internal generations are independent of the optional user-facing iteration archives.
///
/// ```no_run
/// # use gammalooprs::integrate::{IntegrationCheckpoint, IntegrationSlot, IntegrationState,
/// #     IntegrationWorkspace, IntegrationWorkspaceManifest, SlotIndex};
/// # use std::error::Error;
/// # use std::path::Path;
/// # fn checkpoint(root: &Path, manifest: IntegrationWorkspaceManifest,
/// #     state: &IntegrationState, slots: &[IntegrationSlot]) -> Result<(), Box<dyn Error>> {
/// let mut workspace = IntegrationWorkspace::acquire(root, manifest)?;
/// // Run workers and merge a complete iteration before publishing it.
/// workspace.commit_completed_iteration(state, slots)?;
/// let pinned = IntegrationCheckpoint::open(workspace.root())?;
/// assert_eq!(pinned.read_state()?.iter, state.iter);
/// let loaded = pinned.load()?;
/// assert_eq!(loaded.state.iter, state.iter);
/// let settings = pinned.read_slot_settings(SlotIndex::from(0))?;
/// let observables = pinned.read_observables(SlotIndex::from(0))?;
/// # Ok(())
/// # }
/// ```
#[derive(Debug)]
pub struct IntegrationWorkspace {
    writer: SnapshotWriter,
    manifest: IntegrationWorkspaceManifest,
    committed: Option<IntegrationCheckpoint>,
    restart_outputs: Option<RestartOutputs>,
    needs_recovery: bool,
}

impl IntegrationWorkspace {
    /// Acquire before reading for resume or preparing restart.
    /// Recovery finishes cleanup for an already committed restart and discards diagnostics
    /// from abandoned attempts. It also removes recognized stale output temporaries.
    ///
    /// # Errors
    /// Writer contention and invalid configuration fail before recovery. Recovery may remove
    /// old outputs and generations before an I/O or journal error is returned; the committed
    /// checkpoint remains authoritative. See the [type example](Self).
    pub fn acquire(root: &Path, manifest: IntegrationWorkspaceManifest) -> CheckpointResult<Self> {
        IntegrationCheckpoint::validate_workspace(&manifest)?;
        let mut workspace = Self {
            writer: SnapshotWriter::acquire(
                root,
                &[ManifestName::new(INTEGRATION_MANIFEST_FILE)?],
            )?,
            manifest,
            committed: None,
            restart_outputs: None,
            needs_recovery: false,
        };
        Self::check_no_saved_state(workspace.writer.checked_root())?;
        workspace.recover_restart_outputs()?;
        workspace.remove_output_temporaries()?;
        let manifest_path = workspace
            .writer
            .checked_root()
            .join(INTEGRATION_MANIFEST_FILE)?;
        match fs::symlink_metadata(manifest_path.as_path()) {
            Err(source) if source.kind() == ErrorKind::NotFound => {
                workspace.check_no_previous_outputs()?;
            }
            Err(source) => {
                return Err(IntegrationCheckpointError::Io {
                    path: manifest_path.as_path().to_owned(),
                    source,
                });
            }
            Ok(_) => (),
        }
        Ok(workspace)
    }

    fn check_no_saved_state(root: &OutputPath) -> CheckpointResult<()> {
        let state_manifest = root.join(STATE_MANIFEST_FILE)?;
        match fs::symlink_metadata(state_manifest.as_path()) {
            Ok(_) => {
                return Err(IntegrationCheckpointError::SavedStateWorkspace {
                    path: state_manifest.as_path().to_owned(),
                });
            }
            Err(source) if source.kind() == ErrorKind::NotFound => (),
            Err(source) => {
                return Err(IntegrationCheckpointError::Io {
                    path: state_manifest.as_path().to_owned(),
                    source,
                });
            }
        }
        Ok(())
    }

    /// Logical workspace root, used for presentation outputs. See the [type example](Self).
    #[must_use]
    pub fn root(&self) -> &Path {
        self.writer.root()
    }

    /// Validated configuration for this run. See the [type example](Self).
    #[must_use]
    pub fn manifest(&self) -> &IntegrationWorkspaceManifest {
        &self.manifest
    }

    /// Preserve saved targets when resuming, before any iteration is committed.
    ///
    /// # Errors
    /// Rejects changes after publication or a target count differing from the slot count.
    /// See the [type example](Self).
    pub fn set_targets(&mut self, targets: Vec<Option<Complex<F<f64>>>>) -> CheckpointResult<()> {
        if self.committed.is_some() {
            return Err(IntegrationCheckpointError::ConfigurationCommitted);
        }
        if targets.len() != self.manifest.slots.len() {
            return Err(IntegrationCheckpointError::SlotTargets);
        }
        self.manifest.targets = targets;
        Ok(())
    }

    /// Store warmed-source fingerprints before the first completed iteration.
    ///
    /// # Errors
    /// Rejects changes after publication or a fingerprint count differing from the slot count.
    /// See the [type example](Self).
    pub fn set_integrand_fingerprints(
        &mut self,
        fingerprints: Vec<String>,
    ) -> CheckpointResult<()> {
        if self.committed.is_some() {
            return Err(IntegrationCheckpointError::ConfigurationCommitted);
        }
        if fingerprints.len() != self.manifest.slots.len() {
            return Err(IntegrationCheckpointError::CommitSlots);
        }
        self.manifest.integrand_fingerprints = fingerprints;
        Ok(())
    }

    /// Destination for this run's diagnostics, isolated until a restart first commits.
    /// The root is resolved once; each writer checks only the descendants it uses.
    ///
    /// # Errors
    /// Returns [`IntegrationCheckpointError::RecoveryRequired`] after an incomplete
    /// postcommit transition. See [`Self::prepare_restart`] and the [type example](Self).
    pub fn output_root(&self) -> CheckpointResult<&OutputPath> {
        if self.needs_recovery {
            return Err(IntegrationCheckpointError::RecoveryRequired {
                path: self.root().to_owned(),
            });
        }
        Ok(match &self.restart_outputs {
            Some(restart) => &restart.outputs,
            None => self.writer.checked_root(),
        })
    }

    /// Stage restart diagnostics without modifying the replaced run's presentation files.
    ///
    /// # Errors
    /// An existing checkpoint must be fully readable before its output ownership can be
    /// replaced. A missing manifest with old outputs is refused. Symlinks and special files
    /// anywhere below the workspace root, and changes after publication, fail.
    /// Call after [`Self::acquire`], before any workers run; see the [type example](Self).
    pub fn prepare_restart(&mut self) -> CheckpointResult<()> {
        if self.committed.is_some() {
            return Err(IntegrationCheckpointError::ConfigurationCommitted);
        }
        if self.restart_outputs.is_some() {
            return Ok(());
        }
        self.writer.checked_root().check_tree()?;
        if self.checkpoint()?.is_none() {
            return Ok(());
        }
        let restart = RestartOutputs::new(self.writer.checked_root())?;
        restart.stage.create_directory()?;
        restart.outputs.create_directory()?;
        self.restart_outputs = Some(restart);
        Ok(())
    }

    /// Resolve once under the held lock; only an absent manifest and outputs mean fresh work.
    ///
    /// # Errors
    /// A present but corrupt checkpoint fails, including missing referenced members. Existing
    /// output files without a manifest are refused. See the [type example](Self).
    pub fn checkpoint(&self) -> CheckpointResult<Option<IntegrationCheckpoint>> {
        let checkpoint = self.checkpoint_metadata()?;
        if let Some(checkpoint) = &checkpoint {
            checkpoint.load()?;
        } else {
            self.check_no_previous_outputs()?;
        }
        Ok(checkpoint)
    }

    fn checkpoint_metadata(&self) -> CheckpointResult<Option<IntegrationCheckpoint>> {
        let root = self.writer.checked_root();
        let checked = root.join(INTEGRATION_MANIFEST_FILE)?;
        let path = checked.as_path();
        match fs::symlink_metadata(path) {
            Ok(_) => {
                let bytes = IntegrationCheckpoint::read(&checked)?;
                IntegrationCheckpoint::decode_manifest(root, &bytes).map(Some)
            }
            Err(source) if source.kind() == ErrorKind::NotFound => Ok(None),
            Err(source) => Err(IntegrationCheckpointError::Io {
                path: path.to_owned(),
                source,
            }),
        }
    }

    fn check_no_previous_outputs(&self) -> CheckpointResult<()> {
        let root = self.writer.checked_root();
        for name in PRESENTATION_ROOTS {
            let checked = root.join(name)?;
            let path = checked.as_path();
            match fs::symlink_metadata(path) {
                Ok(_) => {
                    return Err(IntegrationCheckpointError::UncommittedOutputs {
                        path: path.to_owned(),
                    });
                }
                Err(source) if source.kind() == ErrorKind::NotFound => (),
                Err(source) => {
                    return Err(IntegrationCheckpointError::Io {
                        path: path.to_owned(),
                        source,
                    });
                }
            }
        }
        Ok(())
    }

    /// Remove superseded generations after validating the current checkpoint under its lock.
    /// Returned identities identify removed directories below [`GENERATIONS_DIRECTORY`].
    /// Readers pinned to those generations must finish before this explicit operation.
    ///
    /// ```no_run
    /// # use gammalooprs::integrate::IntegrationWorkspace;
    /// # use std::error::Error;
    /// # fn main() -> Result<(), Box<dyn Error>> {
    /// let removed = IntegrationWorkspace::cleanup_superseded_generations("run".as_ref())?;
    /// assert!(removed.iter().all(|id| id.as_str().starts_with("gen-")));
    /// # Ok(()) }
    /// ```
    ///
    /// # Errors
    /// Refuses an absent or unreadable checkpoint before creating a lock. After acquisition,
    /// revalidates the checkpoint and rejects unsafe generation trees or an active writer.
    /// Filesystem failures retain their affected paths.
    pub fn cleanup_superseded_generations(root: &Path) -> CheckpointResult<Vec<GenerationId>> {
        let checked_root = OutputPath::new(root)?;
        Self::check_no_saved_state(&checked_root)?;
        IntegrationCheckpoint::open_checked(&checked_root)?;
        let mut writer = SnapshotWriter::acquire(
            checked_root.as_path(),
            &[ManifestName::new(INTEGRATION_MANIFEST_FILE)?],
        )?;
        Self::check_no_saved_state(writer.checked_root())?;
        let checkpoint = IntegrationCheckpoint::open_checked(writer.checked_root())?;
        Ok(writer.prune_generations([&checkpoint.manifest.generation])?)
    }

    /// The latest checkpoint published by this run, pinned without rereading a mutable pointer.
    /// A run with no completed work returns `None`; see the [type example](Self).
    #[must_use]
    pub fn committed_checkpoint(&self) -> Option<&IntegrationCheckpoint> {
        self.committed.as_ref()
    }

    /// Stage all members, then replace the manifest exactly once.
    ///
    /// # Errors
    /// Zero or partial iterations and mismatched slots are rejected before publication.
    /// `PostCommit` and `SnapshotError::Committed` identify failures after manifest replacement;
    /// those iterations must not be replayed. See the [type example](Self).
    pub fn commit_completed_iteration(
        &mut self,
        state: &IntegrationState,
        slots: &[IntegrationSlot],
    ) -> CheckpointResult<()> {
        if slots
            .iter()
            .map(|slot| slot.target)
            .ne(self.manifest.targets.iter().copied())
        {
            return Err(IntegrationCheckpointError::SlotTargets);
        }
        let bundles = slots
            .iter()
            .map(|slot| slot.integrand.observable_snapshot_bundle())
            .collect::<Vec<_>>();
        self.commit(
            state,
            slots.iter().map(|slot| (&slot.meta, &slot.settings)),
            &bundles,
            |_| Ok(()),
        )
    }

    fn commit<'a>(
        &mut self,
        state: &IntegrationState,
        slots: impl IntoIterator<Item = (&'a SlotMeta, &'a RuntimeSettings)>,
        bundles: &[Option<ObservableSnapshotBundle>],
        mut before_write: impl FnMut(&Path) -> CheckpointResult<()>,
    ) -> CheckpointResult<()> {
        if self.needs_recovery {
            return Err(IntegrationCheckpointError::RecoveryRequired {
                path: self.root().to_owned(),
            });
        }
        IntegrationCheckpoint::validate_workspace(&self.manifest)?;
        let slots: Vec<_> = slots.into_iter().collect();
        if slots.len() != self.manifest.slots.len()
            || bundles.len() != slots.len()
            || slots
                .iter()
                .map(|(meta, _)| *meta)
                .ne(self.manifest.slots.iter())
        {
            return Err(IntegrationCheckpointError::CommitSlots);
        }
        if state.iter == 0 || state.num_points == 0 {
            return Err(IntegrationCheckpointError::IncompleteIteration {
                iteration: state.iter,
                num_points: state.num_points,
            });
        }
        for (index, ((_, settings), bundle)) in slots.iter().zip(bundles).enumerate() {
            IntegrationCheckpoint::validate_observables(
                SlotIndex::from(index),
                settings,
                bundle.as_ref(),
            )?;
        }
        let root = self.writer.checked_root().clone();
        let previous = self.checkpoint_metadata()?;
        let restarting = self.restart_outputs.is_some();
        if let Some(previous) = previous.as_ref().filter(|_| !restarting)
            && state.iter <= previous.manifest.completed_iteration
        {
            return Err(IntegrationCheckpointError::IterationOrder {
                previous: previous.manifest.completed_iteration,
                iteration: state.iter,
            });
        }
        let manifest_path = integrate::workspace_manifest_path(root.as_path());
        let staged = self.writer.begin()?;
        let manifest = IntegrationCheckpointManifest {
            workspace: self.manifest.clone(),
            generation: staged.id().clone(),
            completed_iteration: state.iter,
            num_points: state.num_points,
            observable_slots: bundles
                .iter()
                .enumerate()
                .filter_map(|(slot, bundle)| bundle.as_ref().map(|_| SlotIndex::from(slot)))
                .collect(),
        };
        let checkpoint = IntegrationCheckpoint {
            manifest,
            payload_root: staged.output_path()?,
            root: root.clone(),
        };
        checkpoint.validate_state(state)?;
        let state_path = checkpoint.state_path();
        let state_bytes = bincode::encode_to_vec(state, config::standard()).map_err(|source| {
            IntegrationCheckpointError::StateEncode {
                path: state_path.clone(),
                source,
            }
        })?;
        let checked_state = checkpoint.payload_root.join(STATE_FILE)?;
        Self::write(&checked_state, &state_bytes, &mut before_write)?;
        for (index, ((_, settings), bundle)) in slots.iter().zip(bundles).enumerate() {
            let slot = SlotIndex::from(index);
            let slot_path = checkpoint.slot_path(slot)?;
            slot_path.ensure_directory()?;
            let checked = slot_path.join(SLOT_SETTINGS_FILE)?;
            let path = checked.as_path().to_owned();
            let text =
                toml::to_string_pretty(&checkpoint.member(slot, *settings)).map_err(|source| {
                    IntegrationCheckpointError::SettingsEncode {
                        path: path.clone(),
                        source,
                    }
                })?;
            Self::write(&checked, text.as_bytes(), &mut before_write)?;
            if let Some(bundle) = bundle {
                let path = checkpoint
                    .observable_path(slot)
                    .expect("bundle membership is declared");
                let bytes = serde_json::to_vec_pretty(&checkpoint.member(slot, bundle)).map_err(
                    |source| IntegrationCheckpointError::Json {
                        path: path.clone(),
                        source,
                    },
                )?;
                let checked = slot_path.join(OBSERVABLES_FILE)?;
                Self::write(&checked, &bytes, &mut before_write)?;
            }
        }
        let bytes = serde_json::to_vec_pretty(&checkpoint.manifest).map_err(|source| {
            IntegrationCheckpointError::Json {
                path: manifest_path.clone(),
                source,
            }
        })?;
        // JSON encoding can succeed while replacing nonfinite numbers with null. Validate the
        // complete serialized candidate through the same owner that will load it after commit.
        // No staged file is modified after this boundary.
        let committed = IntegrationCheckpoint::from_manifest_bytes(&root, &bytes)?;
        before_write(&manifest_path)?;
        if let Some(restart) = &self.restart_outputs {
            restart.prepare_commit(&root, &committed.manifest.generation)?;
        }
        if let Err(error) = staged.publish(&bytes, ManifestName::new(INTEGRATION_MANIFEST_FILE)?) {
            if matches!(&error, SnapshotError::Committed { .. }) {
                self.committed = Some(committed);
                self.needs_recovery = true;
            }
            return Err(error.into());
        }
        self.committed = Some(committed);
        let mut cleanup = || -> CheckpointResult<()> {
            let current = self
                .committed
                .as_ref()
                .expect("publication installed the descriptor");
            if let Some(restart) = &self.restart_outputs {
                before_write(restart.stage.join("cleared")?.as_path())?;
                restart.finish(
                    &mut self.writer,
                    &current.manifest.generation,
                    &mut before_write,
                )?;
                self.restart_outputs = None;
            }
            Ok(())
        };
        match cleanup() {
            Ok(()) => Ok(()),
            Err(source) => {
                self.needs_recovery = true;
                Err(IntegrationCheckpointError::PostCommit {
                    iteration: state.iter,
                    source: Box::new(source),
                })
            }
        }
    }

    fn recover_restart_outputs(&mut self) -> CheckpointResult<()> {
        let root = self.writer.checked_root().clone();
        let retired = root.join(COMPLETED_RESTART_DIRECTORY)?;
        match fs::symlink_metadata(retired.as_path()) {
            Ok(_) => {
                // The active journal is atomically renamed only after output promotion and
                // pruning finish. Its entries may already be partly deleted; never replay it.
                self.checkpoint()?
                    .ok_or_else(|| IntegrationCheckpointError::RestartJournal {
                        path: retired.as_path().to_owned(),
                    })?;
                retired.remove_tree()?;
                root.sync_directory()?;
            }
            Err(source) if source.kind() == ErrorKind::NotFound => (),
            Err(source) => {
                return Err(IntegrationCheckpointError::Io {
                    path: retired.as_path().to_owned(),
                    source,
                });
            }
        }
        let restart = RestartOutputs::new(&root)?;
        let stage = &restart.stage;
        match fs::symlink_metadata(stage.as_path()) {
            Err(source) if source.kind() == ErrorKind::NotFound => return Ok(()),
            Err(source) => {
                return Err(IntegrationCheckpointError::Io {
                    path: stage.as_path().to_owned(),
                    source,
                });
            }
            Ok(_) => stage.check_tree()?,
        }
        let checkpoint =
            self.checkpoint()?
                .ok_or_else(|| IntegrationCheckpointError::RestartJournal {
                    path: stage.as_path().to_owned(),
                })?;
        if restart.committed_generation()?.as_ref() == Some(&checkpoint.manifest.generation) {
            // A previous writer may have stopped after rename but before directory sync.
            root.sync_directory()?;
            restart.finish(
                &mut self.writer,
                &checkpoint.manifest.generation,
                &mut |_| Ok(()),
            )?;
        } else {
            warn!(
                path = %restart.outputs.as_path().display(),
                "Discarding staged diagnostics from an abandoned restart"
            );
            // This attempt never committed. Its private diagnostics cannot replace the old run.
            stage.remove_tree()?;
            root.sync_directory()?;
        }
        Ok(())
    }

    fn remove_output_temporaries(&self) -> CheckpointResult<()> {
        let root = self.writer.checked_root();
        let mut pending = vec![root.clone()];
        while let Some(directory) = pending.pop() {
            for entry in fs::read_dir(directory.as_path()).map_err(|source| {
                IntegrationCheckpointError::Io {
                    path: directory.as_path().to_owned(),
                    source,
                }
            })? {
                let entry = entry.map_err(|source| IntegrationCheckpointError::Io {
                    path: directory.as_path().to_owned(),
                    source,
                })?;
                let name = entry.file_name();
                let Some(name) = name.to_str() else {
                    continue;
                };
                let in_presentation = directory.as_path() != root.as_path();
                if !in_presentation
                    && !PRESENTATION_ROOTS.contains(&name)
                    && !Self::is_output_temporary(name)
                {
                    continue;
                }
                let file_type =
                    entry
                        .file_type()
                        .map_err(|source| IntegrationCheckpointError::Io {
                            path: entry.path(),
                            source,
                        })?;
                // Inspect, but never follow unrelated symlinks or special entries. Actual
                // presentation writes and restart replacement validate the trees they own.
                if file_type.is_dir() {
                    if in_presentation || PRESENTATION_ROOTS.contains(&name) {
                        pending.push(directory.join(name)?);
                    }
                } else if Self::is_output_temporary(name) {
                    directory.join(name)?.remove_tree()?;
                }
            }
        }
        Ok(())
    }

    fn is_output_temporary(name: &str) -> bool {
        let Some(identity) = name
            .strip_prefix(PRESENTATION_TEMP_PREFIX)
            .and_then(|name| name.strip_suffix(PRESENTATION_TEMP_SUFFIX))
        else {
            return false;
        };
        let Some((pid, sequence)) = identity.split_once('-') else {
            return false;
        };
        !pid.is_empty()
            && pid.bytes().all(|byte| byte.is_ascii_digit())
            && !sequence.is_empty()
            && sequence
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    }

    fn write(
        checked: &OutputPath,
        bytes: &[u8],
        before_write: &mut impl FnMut(&Path) -> CheckpointResult<()>,
    ) -> CheckpointResult<()> {
        let path = checked.as_path();
        before_write(path)?;
        checked.check()?;
        let write = || -> Result<(), IoError> {
            let mut file = File::create_new(path)?;
            file.write_all(bytes)?;
            file.flush()
        };
        write().map_err(|source| IntegrationCheckpointError::Io {
            path: path.to_owned(),
            source,
        })
    }
}

impl Drop for IntegrationWorkspace {
    fn drop(&mut self) {
        if self.committed.is_none()
            && let Some(restart) = &self.restart_outputs
        {
            warn!(
                path = %restart.outputs.as_path().display(),
                "Restart did not commit; staged diagnostics remain \
                 until the next writer acquisition"
            );
        }
    }
}

const RESTART_OUTPUT_DIRECTORY: &str = ".restart-output";
const COMPLETED_RESTART_DIRECTORY: &str = ".restart-complete";
// Workspace-root entries written by integration presentation. Failed samples live under
// numerical_stability, and per-slot observables live under integrands/<slot>.
const PRESENTATION_ROOTS: [&str; 5] = [
    INTEGRATION_RESULT_FILE,
    RESULTS_DIRECTORY,
    INTEGRANDS_DIRECTORY,
    NUMERICAL_STABILITY_DIRECTORY,
    INTERRUPTED_DIRECTORY,
];

/// Recovery journal for the transition from a replaced run to a new run's presentation.
/// Its generation is written before publication, but authorizes deletion only when that exact
/// generation has become the validated, durable workspace manifest under the writer lock.
#[derive(Debug)]
struct RestartOutputs {
    stage: OutputPath,
    outputs: OutputPath,
}

impl RestartOutputs {
    fn new(root: &OutputPath) -> CheckpointResult<Self> {
        let stage = root.join(RESTART_OUTPUT_DIRECTORY)?;
        let outputs = stage.join("outputs")?;
        Ok(Self { stage, outputs })
    }

    fn prepare_commit(&self, root: &OutputPath, generation: &GenerationId) -> CheckpointResult<()> {
        let stage = &self.stage;
        // Check the managed tree before either output cleanup or generation pruning can begin.
        root.check_tree()?;
        let path = stage.join("commit.json")?;
        let temporary = stage.join("commit.tmp")?;
        let bytes =
            serde_json::to_vec(generation).map_err(|source| IntegrationCheckpointError::Json {
                path: path.as_path().to_owned(),
                source,
            })?;
        // Only this writer owns the private journal and its temporary.
        temporary.write(&bytes)?;
        // Diagnostics accumulated before the checkpoint must survive promotion too.
        stage.sync_tree()?;
        fs::rename(temporary.as_path(), path.as_path()).map_err(|source| {
            IntegrationCheckpointError::Io {
                path: path.as_path().to_owned(),
                source,
            }
        })?;
        stage.sync_directory()?;
        root.sync_directory()?;
        Ok(())
    }

    fn committed_generation(&self) -> CheckpointResult<Option<GenerationId>> {
        let path = self.stage.join("commit.json")?.as_path().to_owned();
        match fs::read(&path) {
            Ok(bytes) => serde_json::from_slice(&bytes)
                .map(Some)
                .map_err(|source| IntegrationCheckpointError::Json { path, source }),
            Err(source) if source.kind() == ErrorKind::NotFound => Ok(None),
            Err(source) => Err(IntegrationCheckpointError::Io { path, source }),
        }
    }

    fn finish(
        &self,
        writer: &mut SnapshotWriter,
        generation: &GenerationId,
        before_write: &mut impl FnMut(&Path) -> CheckpointResult<()>,
    ) -> CheckpointResult<()> {
        let root = writer.checked_root().clone();
        let stage = &self.stage;
        root.check_tree()?;
        if self.committed_generation()?.as_ref() != Some(generation) {
            return Err(IntegrationCheckpointError::RestartJournal {
                path: stage.as_path().to_owned(),
            });
        }
        // A previous attempt may have renamed the marker but stopped before syncing this
        // directory. Establish its durability before trusting it to skip output deletion.
        stage.sync_directory()?;
        let cleared = stage.join("cleared")?;
        let already_cleared = match fs::read(cleared.as_path()) {
            Ok(bytes) if bytes == generation.as_str().as_bytes() => true,
            Ok(_) => {
                return Err(IntegrationCheckpointError::RestartJournal {
                    path: cleared.as_path().to_owned(),
                });
            }
            Err(source) if source.kind() == ErrorKind::NotFound => false,
            Err(source) => {
                return Err(IntegrationCheckpointError::Io {
                    path: cleared.as_path().to_owned(),
                    source,
                });
            }
        };
        if !already_cleared {
            let old = PRESENTATION_ROOTS
                .into_iter()
                .map(|name| root.join(name))
                .collect::<Result<Vec<_>, _>>()?;
            for path in &old {
                path.check_tree()?;
            }
            for path in &old {
                path.remove_tree()?;
                before_write(path.as_path())?;
            }
            root.sync_directory()?;
            // Persist the barrier before promoting any new output: retries never delete it.
            Self::write_cleared_marker(stage, generation, before_write)?;
        }
        // Keep the journal through pruning. A stopped transition must recover even if the
        // resumed integrator has no remaining iterations and never publishes another checkpoint.
        before_write(root.join(GENERATIONS_DIRECTORY)?.as_path())?;
        writer.prune_generations([generation])?;
        let outputs = &self.outputs;
        for name in PRESENTATION_ROOTS {
            let source_path = outputs.join(name)?;
            let target = root.join(name)?;
            match fs::symlink_metadata(source_path.as_path()) {
                Err(source) if source.kind() == ErrorKind::NotFound => continue,
                Err(source) => {
                    return Err(IntegrationCheckpointError::Io {
                        path: source_path.as_path().to_owned(),
                        source,
                    });
                }
                Ok(_) => (),
            }
            // A completed rename leaves no source; an existing source and target is corruption.
            match fs::symlink_metadata(target.as_path()) {
                Err(source) if source.kind() == ErrorKind::NotFound => (),
                Err(source) => {
                    return Err(IntegrationCheckpointError::Io {
                        path: target.as_path().to_owned(),
                        source,
                    });
                }
                Ok(_) => {
                    return Err(IntegrationCheckpointError::RestartJournal {
                        path: target.as_path().to_owned(),
                    });
                }
            }
            fs::rename(source_path.as_path(), target.as_path()).map_err(|source| {
                IntegrationCheckpointError::Io {
                    path: target.as_path().to_owned(),
                    source,
                }
            })?;
            before_write(target.as_path())?;
            // Both sides must be synchronized before deleting the recovery journal.
            root.sync_directory()?;
            outputs.sync_directory()?;
        }
        // Retire the whole journal atomically before deleting any of its phase markers. A crash
        // during recursive removal must never make completed promotion look like pending cleanup.
        let retired = root.join(COMPLETED_RESTART_DIRECTORY)?;
        fs::rename(stage.as_path(), retired.as_path()).map_err(|source| {
            IntegrationCheckpointError::Io {
                path: retired.as_path().to_owned(),
                source,
            }
        })?;
        root.sync_directory()?;
        before_write(retired.as_path())?;
        retired.remove_tree()?;
        root.sync_directory()?;
        Ok(())
    }

    fn write_cleared_marker(
        stage: &OutputPath,
        generation: &GenerationId,
        before_write: &mut impl FnMut(&Path) -> CheckpointResult<()>,
    ) -> CheckpointResult<()> {
        let temporary = stage.join("cleared.tmp")?;
        let marker = stage.join("cleared")?;
        let bytes = generation.as_str().as_bytes();
        let mut file = temporary.create_file()?;
        file.write_all(&bytes[..bytes.len() / 2])
            .map_err(|source| IntegrationCheckpointError::Io {
                path: temporary.as_path().to_owned(),
                source,
            })?;
        // Test barriers run after a real partial write. Only this unpublished temporary can be
        // incomplete; recovery never mistakes it for permission to promote new presentation.
        before_write(temporary.as_path())?;
        file.write_all(&bytes[bytes.len() / 2..])
            .and_then(|()| file.sync_all())
            .map_err(|source| IntegrationCheckpointError::Io {
                path: temporary.as_path().to_owned(),
                source,
            })?;
        fs::rename(temporary.as_path(), marker.as_path()).map_err(|source| {
            IntegrationCheckpointError::Io {
                path: marker.as_path().to_owned(),
                source,
            }
        })?;
        before_write(marker.as_path())?;
        stage.sync_directory()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrate::{SamplingCorrelationMode, SamplingSlotState};
    use crate::model::SerializableInputParamCard;
    use crate::observables::{
        HistogramAccumulatorState, HistogramSettings, ObservablePhase, ObservableValueTransform,
    };
    use crate::utils::F;
    use spenso::algebra::complex::Complex;
    use std::collections::BTreeMap;
    use std::error::Error as StdError;
    #[cfg(unix)]
    use std::os::unix::fs as unix_fs;
    #[cfg(unix)]
    use std::os::unix::fs::MetadataExt;
    #[cfg(unix)]
    use std::os::unix::net::UnixListener;
    use symbolica::numerical_integration::{ContinuousGrid, Grid, MonteCarloRng, Sample};

    struct Fixture {
        state: IntegrationState,
        settings: Vec<RuntimeSettings>,
        histograms: Vec<HistogramAccumulatorState>,
    }

    impl Fixture {
        fn new() -> Self {
            let metas = (0..2)
                .map(|slot| SlotMeta {
                    process_name: "process".into(),
                    integrand_name: format!("slot{slot}"),
                })
                .collect::<Vec<_>>();
            let settings = (0..2)
                .map(|slot| {
                    let mut settings: RuntimeSettings = toml::from_str(&format!(
                        r#"
                    [observables.histogram]
                    quantity = "energy"
                    title = "slot{slot}"
                    x_min = 0.0
                    x_max = 1.0
                    n_bins = 2
                "#
                    ))
                    .unwrap();
                    settings.integrator.seed = 7331 + slot;
                    settings
                })
                .collect();
            let state = IntegrationState::new_from_settings(
                SamplingCorrelationMode::Correlated,
                vec![SamplingSlotState::new(
                    Grid::Continuous(ContinuousGrid::new(1, 2, 8, None, false).unwrap()),
                    Vec::new(),
                )],
                metas,
                None,
                None,
                None,
                vec![None, None],
            );
            let histograms = (0..2)
                .map(|slot| {
                    HistogramAccumulatorState::continuous(
                        format!("slot{slot}"),
                        "AL".into(),
                        ObservablePhase::Real,
                        ObservableValueTransform::Identity,
                        0.0,
                        1.0,
                        false,
                        true,
                        2,
                    )
                })
                .collect();
            Self {
                state,
                settings,
                histograms,
            }
        }

        fn manifest(&self) -> IntegrationWorkspaceManifest {
            IntegrationWorkspaceManifest {
                version: IntegrationWorkspaceManifest::CURRENT_VERSION,
                slots: self.state.slot_metas.clone(),
                targets: vec![None, None],
                effective_model_parameters: vec![
                    SerializableInputParamCard {
                        data: BTreeMap::new()
                    };
                    2
                ],
                integrand_fingerprints: vec!["fingerprint0".into(), "fingerprint1".into()],
                reference_gaussians: Vec::new(),
                reference_observable_convention: None,
                training_slot: SlotIndex::from(0),
                integrator_settings_slot: SlotIndex::from(0),
                sampling_correlation_mode: SamplingCorrelationMode::Correlated,
            }
        }

        // Exercise the real grid/statistics/observable owners and serializers with deterministic
        // small data. The API differential tests separately certify actual worker/resume wiring.
        fn advance(&mut self) {
            self.advance_with_histogram_weight(None);
        }

        fn advance_with_histogram_weight(&mut self, first_slot1_weight: Option<f64>) {
            let mut grid = self.state.sampling_states[0].grid.clone_without_samples();
            let mut rng = MonteCarloRng::new(7331 + self.state.iter as u64, 0);
            for sample_index in 0..8 {
                let mut sample = Sample::new();
                grid.sample(&mut rng, &mut sample);
                let Sample::Continuous(weight, xs) = &sample else {
                    panic!("continuous fixture")
                };
                let value = xs[0];
                for (slot, (accumulator, histogram)) in self
                    .state
                    .all_integrals
                    .iter_mut()
                    .zip(&mut self.histograms)
                    .enumerate()
                {
                    let result = Complex::new(value * F((slot + 1) as f64), value * value);
                    accumulator.add_sample(result, None, *weight, Some(&sample));
                    let histogram_weight = if slot == 1 && sample_index == 0 {
                        first_slot1_weight.unwrap_or(result.re.0 * weight.0)
                    } else {
                        result.re.0 * weight.0
                    };
                    histogram
                        .fill_continuous_sample(&[
                            (value.0, histogram_weight),
                            (
                                if sample_index % 2 == 0 { -0.25 } else { 1.25 },
                                result.im.0,
                            ),
                        ])
                        .unwrap();
                }
                grid.add_training_sample(&sample, value).unwrap();
            }
            grid.update(F(1.0), F(1.0));
            self.state.sampling_states[0].grid = grid;
            self.state.update_iter(false);
            self.state.iter += 1;
            self.state.num_points += 8;
            for histogram in &mut self.histograms {
                histogram.update_results();
            }
        }

        fn bundles(&self) -> Vec<Option<ObservableSnapshotBundle>> {
            self.histograms
                .iter()
                .map(|histogram| {
                    Some(ObservableSnapshotBundle {
                        histograms: BTreeMap::from([("histogram".into(), histogram.snapshot())]),
                    })
                })
                .collect()
        }

        fn commit(
            &self,
            workspace: &mut IntegrationWorkspace,
            hook: impl FnMut(&Path) -> CheckpointResult<()>,
        ) -> CheckpointResult<()> {
            workspace.commit(
                &self.state,
                self.state.slot_metas.iter().zip(&self.settings),
                &self.bundles(),
                hook,
            )
        }

        fn restore(checkpoint: &IntegrationCheckpoint) -> Self {
            Self {
                state: checkpoint.read_state().unwrap(),
                settings: (0..2)
                    .map(|slot| {
                        checkpoint
                            .read_slot_settings(SlotIndex::from(slot))
                            .unwrap()
                    })
                    .collect(),
                histograms: (0..2)
                    .map(|slot| {
                        checkpoint
                            .read_observables(SlotIndex::from(slot))
                            .unwrap()
                            .unwrap()
                            .histograms
                            .remove("histogram")
                            .unwrap()
                            .into_accumulator_state()
                    })
                    .collect(),
            }
        }

        fn assert_same(&self, other: &Self) {
            assert_eq!(
                bincode::encode_to_vec(&self.state, config::standard()).unwrap(),
                bincode::encode_to_vec(&other.state, config::standard()).unwrap()
            );
            assert_eq!(self.settings, other.settings);
            for (actual, expected) in self.histograms.iter().zip(&other.histograms) {
                assert_eq!(actual.snapshot(), expected.snapshot());
                assert_eq!(actual.sample_count, self.state.num_points);
            }
        }
    }

    #[test]
    fn checkpoint_nonfinite_optional_bin_edge_preserves_two_slot_commit() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.advance();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let prior = Fixture::restore(&IntegrationCheckpoint::open(root.path()).unwrap());
        let manifest = fs::read(root.path().join("manifest.json")).unwrap();

        let mut candidate = Fixture::new();
        let HistogramSettings::Continuous(settings) = &mut candidate.settings[1]
            .observables
            .get_mut("histogram")
            .unwrap()
            .histogram
        else {
            panic!("continuous fixture")
        };
        settings.x_max = 1.0e308;
        candidate.histograms[1].x_max = 1.0e308;
        candidate.advance();
        candidate.advance();
        let snapshot = candidate.histograms[1].snapshot();
        assert_eq!(snapshot.x_min, Some(0.0));
        assert_eq!(snapshot.x_max, Some(1.0e308));
        assert_eq!(snapshot.bins.len(), 2);
        assert!(snapshot.bins[1].x_max.unwrap().is_infinite());
        assert!(
            snapshot
                .bins
                .iter()
                .chain([&snapshot.underflow_bin, &snapshot.overflow_bin])
                .all(|bin| bin.sum_weights.is_finite() && bin.sum_weights_squared.is_finite())
        );

        assert!(candidate.commit(&mut writer, |_| Ok(())).is_err());
        assert_eq!(
            fs::read(root.path().join("manifest.json")).unwrap(),
            manifest
        );
        assert_eq!(
            fs::read_dir(root.path().join("generations"))
                .unwrap()
                .count(),
            1
        );
        drop(writer);
        let checkpoint = IntegrationCheckpoint::open(root.path()).unwrap();
        assert_eq!(
            checkpoint.manifest().observable_slots,
            vec![SlotIndex::from(0), SlotIndex::from(1)]
        );
        Fixture::restore(&checkpoint).assert_same(&prior);
    }

    #[test]
    fn checkpoint_nonfinite_json_preserves_nonzero_two_slot_commit() {
        for invalid_target in [true, false] {
            let root = CheckpointTestDirectory::new();
            let mut fixture = Fixture::new();
            let mut writer =
                IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
            fixture.advance();
            fixture.commit(&mut writer, |_| Ok(())).unwrap();
            let prior = Fixture::restore(&IntegrationCheckpoint::open(root.path()).unwrap());
            let manifest = fs::read(root.path().join("manifest.json")).unwrap();
            if invalid_target {
                fixture.advance();
                writer.manifest.targets[0] = Some(Complex::new(F(f64::INFINITY), F(0.0)));
            } else {
                // A finite actual fill can overflow a raw second moment. The accumulation
                // algorithm is unchanged; refusing its unreadable JSON is a save boundary.
                fixture.advance_with_histogram_weight(Some(1.0e200));
                assert!(
                    fixture.histograms[1]
                        .snapshot()
                        .bins
                        .iter()
                        .any(|bin| bin.sum_weights_squared.is_infinite())
                );
                assert_eq!(fixture.histograms[1].sample_count, 16);
            }
            assert!(fixture.commit(&mut writer, |_| Ok(())).is_err());
            assert_eq!(
                fs::read(root.path().join("manifest.json")).unwrap(),
                manifest
            );
            assert_eq!(
                fs::read_dir(root.path().join("generations"))
                    .unwrap()
                    .count(),
                1
            );
            drop(writer);
            let checkpoint = IntegrationCheckpoint::open(root.path()).unwrap();
            assert_eq!(
                checkpoint.manifest().observable_slots,
                vec![SlotIndex::from(0), SlotIndex::from(1)]
            );
            Fixture::restore(&checkpoint).assert_same(&prior);
        }
    }

    #[test]
    fn checkpoint_payload_fault_matrix_resumes_nonzero_two_slot_histograms() {
        let mut uninterrupted = Fixture::new();
        for _ in 0..3 {
            uninterrupted.advance();
        }
        // scalar/grid, settings0, histogram0, settings1, histogram1, manifest publication
        for fault in 0..6 {
            let root = CheckpointTestDirectory::new();
            let mut fixture = Fixture::new();
            let mut writer =
                IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
            fixture.advance();
            fixture.commit(&mut writer, |_| Ok(())).unwrap();
            let previous_manifest = fs::read(root.path().join("manifest.json")).unwrap();
            let pinned = IntegrationCheckpoint::open(root.path()).unwrap();
            fixture.advance();
            let mut step = 0;
            let result = fixture.commit(&mut writer, |path| {
                let fail = step == fault;
                step += 1;
                if fail {
                    if path.file_name().unwrap() == "manifest.json" {
                        // The shared primitive separately injects its actual manifest write and
                        // rename errors. This owner injects the handoff failure before publish.
                        return Err(IntegrationCheckpointError::Io {
                            path: path.to_owned(),
                            source: IoError::other("injected manifest handoff failure"),
                        });
                    }
                    // Force a real exclusive-file creation error, without permission assumptions.
                    fs::create_dir(path).unwrap();
                }
                Ok(())
            });
            assert!(result.is_err(), "fault {fault} did not fail");
            assert_eq!(
                fs::read(root.path().join("manifest.json")).unwrap(),
                previous_manifest
            );
            assert_eq!(pinned.read_state().unwrap().iter, 1);
            assert_eq!(
                fs::read_dir(root.path().join("generations"))
                    .unwrap()
                    .count(),
                1
            );
            drop(writer);
            let restored = IntegrationCheckpoint::open(root.path()).unwrap();
            let mut resumed = Fixture::restore(&restored);
            let mut writer =
                IntegrationWorkspace::acquire(root.path(), resumed.manifest()).unwrap();
            while resumed.state.iter < 3 {
                resumed.advance();
                resumed.commit(&mut writer, |_| Ok(())).unwrap();
                if resumed.state.iter == 2 {
                    assert_eq!(pinned.read_state().unwrap().iter, 1);
                }
            }
            Fixture::restore(&IntegrationCheckpoint::open(root.path()).unwrap())
                .assert_same(&uninterrupted);
            // A previously pinned reader keeps its old generation after both later publications.
            assert_eq!(pinned.read_state().unwrap().iter, 1);
            assert_eq!(
                fs::read_dir(root.path().join(GENERATIONS_DIRECTORY))
                    .unwrap()
                    .count(),
                3
            );
        }
    }

    #[test]
    fn checkpoint_cleanup_does_not_initialize_missing_or_nonworkspace_targets() {
        let root = CheckpointTestDirectory::new();
        let missing = root.path().join("missing");
        let error = IntegrationWorkspace::cleanup_superseded_generations(&missing).unwrap_err();
        assert!(!error.is_busy());
        assert!(
            matches!(error, IntegrationCheckpointError::Io { source, .. }
            if source.kind() == ErrorKind::NotFound)
        );
        assert!(!missing.exists());

        let unrelated = root.path().join("unrelated");
        fs::create_dir(&unrelated).unwrap();
        fs::write(unrelated.join("keep.txt"), b"unrelated bytes").unwrap();
        assert!(IntegrationWorkspace::cleanup_superseded_generations(&unrelated).is_err());
        assert!(!unrelated.join(".snapshot.lock").exists());
        assert_eq!(
            fs::read(unrelated.join("keep.txt")).unwrap(),
            b"unrelated bytes"
        );

        let mut fixture = Fixture::new();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.advance();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let payload = writer
            .committed_checkpoint()
            .unwrap()
            .payload_root
            .as_path()
            .to_owned();
        assert!(IntegrationWorkspace::cleanup_superseded_generations(&payload).is_err());
        assert!(!payload.join(".snapshot.lock").exists());
        assert_eq!(
            writer
                .checkpoint()
                .unwrap()
                .unwrap()
                .read_state()
                .unwrap()
                .iter,
            1
        );
    }

    #[test]
    fn checkpoint_cleanup_keeps_current_and_ignores_strays_under_writer_lock() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        let mut generations = Vec::new();
        for _ in 0..3 {
            fixture.advance();
            fixture.commit(&mut writer, |_| Ok(())).unwrap();
            generations.push(
                writer
                    .committed_checkpoint()
                    .unwrap()
                    .manifest
                    .generation
                    .clone(),
            );
        }
        let generation_root = root.path().join(GENERATIONS_DIRECTORY);
        let stray = generation_root.join(".DS_Store");
        fs::write(&stray, b"unrelated metadata").unwrap();
        let unrelated_directory = generation_root.join("notes");
        fs::create_dir(&unrelated_directory).unwrap();
        fs::write(unrelated_directory.join("keep.txt"), b"user note").unwrap();
        let pointer = fs::read(root.path().join(INTEGRATION_MANIFEST_FILE)).unwrap();
        let error = IntegrationWorkspace::cleanup_superseded_generations(root.path()).unwrap_err();
        assert!(error.is_busy());
        assert!(matches!(
            error,
            IntegrationCheckpointError::Snapshot(SnapshotError::Busy { .. })
        ));
        for generation in &generations {
            assert!(generation_root.join(generation.as_str()).is_dir());
        }
        drop(writer);
        let removed = IntegrationWorkspace::cleanup_superseded_generations(root.path()).unwrap();
        assert_eq!(removed.len(), 2);
        assert!(removed.contains(&generations[0]));
        assert!(removed.contains(&generations[1]));
        assert!(!removed.contains(&generations[2]));
        for generation in &removed {
            assert!(!generation_root.join(generation.as_str()).exists());
        }
        assert_eq!(fs::read(&stray).unwrap(), b"unrelated metadata");
        assert_eq!(
            fs::read(unrelated_directory.join("keep.txt")).unwrap(),
            b"user note"
        );
        assert_eq!(
            fs::read(root.path().join(INTEGRATION_MANIFEST_FILE)).unwrap(),
            pointer
        );
        assert_eq!(
            IntegrationCheckpoint::open(root.path())
                .unwrap()
                .read_state()
                .unwrap()
                .iter,
            3
        );
        let repeated_cleanup =
            IntegrationWorkspace::cleanup_superseded_generations(root.path()).unwrap();
        assert!(repeated_cleanup.is_empty());
    }

    #[test]
    fn checkpoint_cleanup_validates_current_payload_before_pruning() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.advance();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let first = writer
            .committed_checkpoint()
            .unwrap()
            .manifest
            .generation
            .clone();
        fixture.advance();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let state = writer.committed_checkpoint().unwrap().state_path();
        drop(writer);
        fs::write(&state, []).unwrap();
        let error = IntegrationWorkspace::cleanup_superseded_generations(root.path()).unwrap_err();
        assert!(matches!(
            error,
            IntegrationCheckpointError::StateDecode { path, .. } if path == state
        ));
        assert!(
            root.path()
                .join(GENERATIONS_DIRECTORY)
                .join(first.as_str())
                .is_dir()
        );
    }

    #[test]
    fn checkpoint_restart_refuses_uncommitted_outputs_without_changing_them() {
        for name in PRESENTATION_ROOTS {
            let root = CheckpointTestDirectory::new();
            let fixture = Fixture::new();
            let mut writer =
                IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
            let output = root.path().join(name);
            fs::write(&output, b"old attempt").unwrap();
            let error = writer.prepare_restart().unwrap_err();
            assert!(
                matches!(error, IntegrationCheckpointError::UncommittedOutputs { path }
                if path == output)
            );
            assert_eq!(fs::read(&output).unwrap(), b"old attempt");
            assert!(!root.path().join(RESTART_OUTPUT_DIRECTORY).exists());
            assert!(!root.path().join(INTEGRATION_MANIFEST_FILE).exists());
            drop(writer);
            let error = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap_err();
            assert!(
                matches!(error, IntegrationCheckpointError::UncommittedOutputs { path }
                if path == output)
            );
            assert_eq!(fs::read(&output).unwrap(), b"old attempt");
        }
    }

    #[test]
    fn checkpoint_restart_reports_undecodable_settings_without_staging_outputs() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.advance();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let path = writer
            .committed_checkpoint()
            .unwrap()
            .slot_settings_path(SlotIndex::from(0));
        drop(writer);
        fs::write(&path, [0xff]).unwrap();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        assert!(matches!(
            writer.prepare_restart(),
            Err(IntegrationCheckpointError::SettingsEncoding { path: actual, .. }) if actual == path
        ));
        assert_eq!(fs::read(&path).unwrap(), [0xff]);
        assert!(!root.path().join(RESTART_OUTPUT_DIRECTORY).exists());
    }

    #[test]
    fn checkpoint_commit_reads_previous_envelope_without_decoding_old_payloads() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.advance();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let old_state = writer.committed_checkpoint().unwrap().state_path();
        let bytes = fs::read(&old_state).unwrap();
        // The writer already owns the validated live state; publication need only compare
        // completed-iteration metadata. Restore the old bytes so the retained reader stays valid.
        fs::write(&old_state, []).unwrap();
        fixture.advance();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        fs::write(&old_state, bytes).unwrap();
        assert_eq!(
            IntegrationCheckpoint::open(root.path())
                .unwrap()
                .read_state()
                .unwrap()
                .iter,
            2
        );
    }

    #[test]
    fn checkpoint_acquisition_removes_only_recognized_stale_output_temporaries() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.advance();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        drop(writer);
        let managed = root.path().join("integrands/slot");
        fs::create_dir_all(&managed).unwrap();
        let stale = [
            root.path().join(".gammaloop-output-123-ab.tmp"),
            managed.join(".gammaloop-output-456-0.tmp"),
        ];
        for path in &stale {
            fs::write(path, b"partial presentation").unwrap();
        }
        let preserved = [
            managed.join(".gammaloop-output-not-a-number.tmp"),
            managed.join("observables_final.json"),
        ];
        for path in &preserved {
            fs::write(path, b"preserved").unwrap();
        }
        let writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        for path in &stale {
            assert!(!path.exists(), "stale output survived: {}", path.display());
        }
        for path in &preserved {
            assert_eq!(fs::read(path).unwrap(), b"preserved");
        }
        assert_eq!(
            writer
                .checkpoint()
                .unwrap()
                .unwrap()
                .read_state()
                .unwrap()
                .iter,
            1
        );
    }

    #[test]
    fn checkpoint_missing_member_and_config_mismatch_are_corruption() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.advance();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let checkpoint = IntegrationCheckpoint::open(root.path()).unwrap();
        let manifest_path = root.path().join(INTEGRATION_MANIFEST_FILE);
        let manifest_bytes = fs::read(&manifest_path).unwrap();
        let members = [
            checkpoint.state_path(),
            checkpoint.slot_settings_path(SlotIndex::from(1)),
            checkpoint.observable_path(SlotIndex::from(1)).unwrap(),
        ];
        let contents: Vec<_> = members.iter().map(|path| fs::read(path).unwrap()).collect();
        for (missing, bytes) in members.iter().zip(&contents) {
            fs::remove_file(missing).unwrap();
            for error in [
                writer.checkpoint().unwrap_err(),
                IntegrationCheckpoint::open(root.path()).unwrap_err(),
            ] {
                let IntegrationCheckpointError::MissingMember { path, source } = error else {
                    panic!("unexpected missing-member error: {error}");
                };
                assert_eq!(&path, missing);
                assert_eq!(source.kind(), ErrorKind::NotFound);
            }
            assert!(!missing.exists());
            for (path, expected) in members.iter().zip(&contents) {
                if path != missing {
                    assert_eq!(&fs::read(path).unwrap(), expected);
                }
            }
            assert_eq!(fs::read(&manifest_path).unwrap(), manifest_bytes);
            fs::write(missing, bytes).unwrap();
        }
        let path = checkpoint.observable_path(SlotIndex::from(1)).unwrap();
        let bytes = fs::read(&path).unwrap();
        let mut wrong_layout: SlotCheckpointMember<
            GenerationId,
            SlotMeta,
            ObservableSnapshotBundle,
        > = serde_json::from_slice(&bytes).unwrap();
        wrong_layout
            .payload
            .histograms
            .get_mut("histogram")
            .unwrap()
            .bins
            .pop();
        fs::write(&path, serde_json::to_vec(&wrong_layout).unwrap()).unwrap();
        assert!(
            IntegrationCheckpoint::open(root.path())
                .err()
                .unwrap()
                .to_string()
                .contains("observable configuration")
        );
        fs::write(&path, &bytes).unwrap();
        // Even valid settings/observable payloads from the other slot must not be accepted.
        fs::copy(
            checkpoint.slot_settings_path(SlotIndex::from(0)),
            checkpoint.slot_settings_path(SlotIndex::from(1)),
        )
        .unwrap();
        assert!(
            IntegrationCheckpoint::open(root.path())
                .err()
                .unwrap()
                .to_string()
                .contains("descriptor")
        );
    }

    #[test]
    fn checkpoint_no_observable_slot_is_explicit_and_partial_state_is_refused() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        fixture.settings[1].observables.clear();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        assert!(
            writer
                .commit(
                    &fixture.state,
                    fixture.state.slot_metas.iter().zip(&fixture.settings),
                    &[None, None],
                    |_| Ok(())
                )
                .is_err()
        );
        assert!(!root.path().join("manifest.json").exists());
        fixture.advance();
        let mut bundles = fixture.bundles();
        bundles[1] = None;
        writer
            .commit(
                &fixture.state,
                fixture.state.slot_metas.iter().zip(&fixture.settings),
                &bundles,
                |_| Ok(()),
            )
            .unwrap();
        let checkpoint = IntegrationCheckpoint::open(root.path()).unwrap();
        assert_eq!(
            checkpoint.manifest().observable_slots,
            vec![SlotIndex::from(0)]
        );
        assert!(
            checkpoint
                .read_observables(SlotIndex::from(1))
                .unwrap()
                .is_none()
        );
        let previous = fs::read(root.path().join("manifest.json")).unwrap();
        // An iteration's scalar accumulator has advanced while its descriptor has not.
        fixture.advance();
        let mut bundles = fixture.bundles();
        bundles[1] = None;
        fixture.state.all_integrals[0].re.update_iter(false);
        assert!(matches!(
            writer.commit(
                &fixture.state,
                fixture.state.slot_metas.iter().zip(&fixture.settings),
                &bundles,
                |_| Ok(())
            ),
            Err(IntegrationCheckpointError::StateDescriptor { .. })
        ));
        assert_eq!(
            fs::read(root.path().join("manifest.json")).unwrap(),
            previous
        );
    }

    #[test]
    fn checkpoint_restart_failure_retains_previous_commit_and_lock_inode() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        fixture.advance();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        #[cfg(unix)]
        let lock_inode = fs::metadata(root.path().join(".snapshot.lock"))
            .unwrap()
            .ino();
        let before = fs::read(root.path().join("manifest.json")).unwrap();
        assert!(matches!(
            IntegrationWorkspace::acquire(root.path(), fixture.manifest()),
            Err(IntegrationCheckpointError::Snapshot(
                SnapshotError::Busy { .. }
            ))
        ));
        drop(writer);
        let mut restart = Fixture::new();
        let mut manifest = restart.manifest();
        manifest.integrand_fingerprints[0] = "new-restart-configuration".into();
        let mut writer = IntegrationWorkspace::acquire(root.path(), manifest).unwrap();
        writer.prepare_restart().unwrap();
        assert!(writer.committed_checkpoint().is_none());
        assert_eq!(fs::read(root.path().join("manifest.json")).unwrap(), before);
        restart.advance();
        assert!(
            restart
                .commit(&mut writer, |path| Err(IntegrationCheckpointError::Io {
                    path: path.to_owned(),
                    source: IoError::other("restart serialization failure"),
                }))
                .is_err()
        );
        assert_eq!(fs::read(root.path().join("manifest.json")).unwrap(), before);
        assert!(root.path().join(".snapshot.lock").is_file());
        #[cfg(unix)]
        assert_eq!(
            fs::metadata(root.path().join(".snapshot.lock"))
                .unwrap()
                .ino(),
            lock_inode
        );
        assert_eq!(
            IntegrationCheckpoint::open(root.path())
                .unwrap()
                .read_state()
                .unwrap()
                .iter,
            1
        );
    }

    #[test]
    fn checkpoint_failed_first_publication_is_fresh_and_retryable() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        fixture.advance();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        assert!(
            fixture
                .commit(&mut writer, |path| {
                    if path.file_name().unwrap() == "manifest.json" {
                        Err(IntegrationCheckpointError::Io {
                            path: path.to_owned(),
                            source: IoError::other("first publication failure"),
                        })
                    } else {
                        Ok(())
                    }
                })
                .is_err()
        );
        assert!(!root.path().join("manifest.json").exists());
        assert!(writer.checkpoint().unwrap().is_none());
        assert_eq!(
            fs::read_dir(root.path().join("generations"))
                .unwrap()
                .count(),
            0
        );
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        Fixture::restore(&IntegrationCheckpoint::open(root.path()).unwrap()).assert_same(&fixture);
    }

    #[test]
    fn checkpoint_reader_never_creates_a_lock_or_repairs_presentation() {
        let source = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        fixture.advance();
        let mut writer = IntegrationWorkspace::acquire(source.path(), fixture.manifest()).unwrap();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let pinned = IntegrationCheckpoint::open(source.path()).unwrap();
        let target = CheckpointTestDirectory::new();
        let destination = target
            .path()
            .join("generations")
            .join(pinned.manifest().generation.as_str());
        fs::create_dir_all(destination.join("slots/0")).unwrap();
        fs::create_dir_all(destination.join("slots/1")).unwrap();
        fs::copy(
            source.path().join("manifest.json"),
            target.path().join("manifest.json"),
        )
        .unwrap();
        fs::copy(pinned.state_path(), destination.join(STATE_FILE)).unwrap();
        for slot in 0..2 {
            fs::copy(
                pinned.slot_settings_path(SlotIndex::from(slot)),
                destination.join(format!("slots/{slot}/settings.toml")),
            )
            .unwrap();
            fs::copy(
                pinned.observable_path(SlotIndex::from(slot)).unwrap(),
                destination.join(format!("slots/{slot}/observables.json")),
            )
            .unwrap();
        }
        let read_only = IntegrationCheckpoint::open(target.path()).unwrap();
        Fixture::restore(&read_only).assert_same(&fixture);
        assert!(!target.path().join(".snapshot.lock").exists());
        assert!(!target.path().join("integration_result.json").exists());
    }

    #[test]
    fn checkpoint_summary_and_resume_reject_versions_before_state_decoding() {
        for version in [
            None,
            Some(2),
            Some(IntegrationWorkspaceManifest::CURRENT_VERSION + 1),
        ] {
            let root = CheckpointTestDirectory::new();
            let mut fixture = Fixture::new();
            fixture.advance();
            let mut writer =
                IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
            fixture.commit(&mut writer, |_| Ok(())).unwrap();
            let checkpoint = writer.committed_checkpoint().unwrap();
            fs::write(checkpoint.state_path(), b"deliberately invalid bincode").unwrap();
            let path = root.path().join("manifest.json");
            let mut manifest = serde_json::to_value(checkpoint.manifest()).unwrap();
            if let Some(version) = version {
                manifest["workspace"]["version"] = version.into();
            } else {
                manifest["workspace"]
                    .as_object_mut()
                    .unwrap()
                    .remove("version");
            }
            let bytes = serde_json::to_vec(&manifest).unwrap();
            fs::write(&path, &bytes).unwrap();
            for error in [
                IntegrationCheckpoint::open(root.path()).unwrap_err(),
                writer.checkpoint().unwrap_err(),
            ] {
                let IntegrationCheckpointError::UnsupportedVersion { found, expected } = error
                else {
                    panic!("unexpected version error: {error}");
                };
                assert_eq!(found, version.unwrap_or(0));
                assert_eq!(expected, IntegrationWorkspaceManifest::CURRENT_VERSION);
            }
            assert_eq!(fs::read(&path).unwrap(), bytes);
        }
    }

    #[test]
    fn checkpoint_resume_propagates_state_io_and_decode_errors() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        fixture.advance();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let path = writer.committed_checkpoint().unwrap().state_path();
        let pointer = fs::read(root.path().join("manifest.json")).unwrap();
        fs::remove_file(&path).unwrap();
        for error in [
            IntegrationCheckpoint::open(root.path()).unwrap_err(),
            writer.checkpoint().unwrap_err(),
        ] {
            let IntegrationCheckpointError::MissingMember {
                path: actual,
                source,
            } = error
            else {
                panic!("unexpected missing-member error: {error}");
            };
            assert_eq!(actual, path);
            assert_eq!(source.kind(), ErrorKind::NotFound);
        }
        fs::create_dir(&path).unwrap();
        for error in [
            IntegrationCheckpoint::open(root.path()).unwrap_err(),
            writer.checkpoint().unwrap_err(),
        ] {
            let IntegrationCheckpointError::NotRegularFile { path: actual } = error else {
                panic!("unexpected non-file member error: {error}");
            };
            assert_eq!(actual, path);
        }
        fs::remove_dir(&path).unwrap();
        fs::write(&path, b"invalid current-version scalar state").unwrap();
        for error in [
            IntegrationCheckpoint::open(root.path()).unwrap_err(),
            writer.checkpoint().unwrap_err(),
        ] {
            let IntegrationCheckpointError::StateDecode { path: actual, .. } = error else {
                panic!("unexpected state-decoding error: {error}");
            };
            assert_eq!(actual, path);
        }
        assert_eq!(
            fs::read(root.path().join("manifest.json")).unwrap(),
            pointer
        );
    }

    fn old_presentation(root: &Path) -> Vec<PathBuf> {
        let paths = [
            "integration_result.json",
            "results/integration_result_iter_0027.json",
            "integrands/removed@slot/observables_final.hwu",
            "integrands/removed@slot/observables_final_iter_0027.json",
            "numerical_stability/failed_samples.jsonl",
            "numerical_stability/removed@slot.hwu",
            "interrupted/iteration_0001/numerical_stability/global.json",
        ];
        paths
            .into_iter()
            .map(|name| {
                let path = root.join(name);
                fs::create_dir_all(path.parent().unwrap()).unwrap();
                fs::write(&path, b"previous run bytes").unwrap();
                path
            })
            .collect()
    }

    #[test]
    fn checkpoint_restart_preserves_all_outputs_until_commit_then_removes_replaced_run() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.advance();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        fixture.advance();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let previous = writer.checkpoint().unwrap().unwrap();
        let pointer = fs::read(root.path().join("manifest.json")).unwrap();
        let old_paths = old_presentation(root.path());
        drop(writer);

        let mut restart = Fixture::new();
        restart.advance();
        let mut writer = IntegrationWorkspace::acquire(root.path(), restart.manifest()).unwrap();
        writer.prepare_restart().unwrap();
        let staged_log = writer
            .output_root()
            .unwrap()
            .as_path()
            .join("numerical_stability/failed_samples.jsonl");
        fs::create_dir_all(staged_log.parent().unwrap()).unwrap();
        fs::write(&staged_log, b"new run failure sample\n").unwrap();
        assert!(
            restart
                .commit(&mut writer, |path| Err(IntegrationCheckpointError::Io {
                    path: path.to_owned(),
                    source: IoError::other("before first restart commit")
                }))
                .is_err()
        );
        assert_eq!(
            fs::read(root.path().join("manifest.json")).unwrap(),
            pointer
        );
        assert_eq!(
            fs::read_dir(root.path().join("generations"))
                .unwrap()
                .count(),
            2
        );
        for path in &old_paths {
            assert_eq!(fs::read(path).unwrap(), b"previous run bytes");
        }
        assert_eq!(previous.read_state().unwrap().iter, 2);

        // A stopped attempt's journal is discarded under the next writer lock; the old manifest
        // still selects its previous run, which can be resumed without restarting.
        drop(writer);
        let mut writer = IntegrationWorkspace::acquire(root.path(), restart.manifest()).unwrap();
        assert_eq!(
            writer
                .checkpoint()
                .unwrap()
                .unwrap()
                .read_state()
                .unwrap()
                .iter,
            2
        );
        assert!(!root.path().join(RESTART_OUTPUT_DIRECTORY).exists());
        for path in &old_paths {
            assert_eq!(fs::read(path).unwrap(), b"previous run bytes");
        }
        writer.prepare_restart().unwrap();
        let staged_log = writer
            .output_root()
            .unwrap()
            .as_path()
            .join("numerical_stability/failed_samples.jsonl");
        fs::create_dir_all(staged_log.parent().unwrap()).unwrap();
        fs::write(&staged_log, b"new run failure sample\n").unwrap();

        restart.commit(&mut writer, |_| Ok(())).unwrap();
        assert_eq!(
            writer
                .committed_checkpoint()
                .unwrap()
                .read_state()
                .unwrap()
                .iter,
            1
        );
        assert_eq!(
            fs::read_dir(root.path().join("generations"))
                .unwrap()
                .count(),
            1
        );
        assert!(previous.read_state().is_err());
        for path in &old_paths {
            if path.ends_with("numerical_stability/failed_samples.jsonl") {
                assert_eq!(fs::read(path).unwrap(), b"new run failure sample\n");
            } else {
                assert!(!path.exists(), "stale output survived: {}", path.display());
            }
        }
        assert!(!root.path().join(RESTART_OUTPUT_DIRECTORY).exists());
    }

    #[test]
    fn checkpoint_restart_recovers_committed_cleanup_before_new_writes() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        fixture.advance();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let old_paths = old_presentation(root.path());
        drop(writer);
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        writer.prepare_restart().unwrap();
        let staged = writer
            .output_root()
            .unwrap()
            .as_path()
            .join("numerical_stability/failed_samples.jsonl");
        fs::create_dir_all(staged.parent().unwrap()).unwrap();
        fs::write(&staged, b"new diagnostics").unwrap();
        let error = fixture
            .commit(&mut writer, |path| {
                if path.ends_with("cleared") {
                    Err(IntegrationCheckpointError::Io {
                        path: path.to_owned(),
                        source: IoError::other("interrupted after durable commit"),
                    })
                } else {
                    Ok(())
                }
            })
            .unwrap_err();
        assert!(matches!(
            error,
            IntegrationCheckpointError::PostCommit { iteration: 1, .. }
        ));
        let generation = writer
            .committed_checkpoint()
            .unwrap()
            .manifest
            .generation
            .clone();
        let journal = root
            .path()
            .join(RESTART_OUTPUT_DIRECTORY)
            .join("commit.json");
        let journal_bytes = fs::read(&journal).unwrap();
        let pending = IntegrationCheckpoint::open(root.path()).unwrap();
        assert!(pending.pending_restart_cleanup().unwrap());
        assert_eq!(fs::read(&journal).unwrap(), journal_bytes);
        for path in &old_paths {
            assert_eq!(fs::read(path).unwrap(), b"previous run bytes");
        }
        drop(writer);
        let writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        let recovered = IntegrationCheckpoint::open(root.path()).unwrap();
        assert!(!recovered.pending_restart_cleanup().unwrap());
        assert_eq!(
            writer.checkpoint().unwrap().unwrap().manifest.generation,
            generation
        );
        assert_eq!(
            fs::read(root.path().join("numerical_stability/failed_samples.jsonl")).unwrap(),
            b"new diagnostics"
        );
        assert!(!root.path().join("results").exists());
        assert!(!root.path().join("integrands").exists());
        assert!(!root.path().join("interrupted").exists());
        assert!(!root.path().join(RESTART_OUTPUT_DIRECTORY).exists());
        assert_eq!(
            fs::read_dir(root.path().join("generations"))
                .unwrap()
                .count(),
            1
        );
    }

    #[test]
    fn checkpoint_restart_recovers_marker_rename_before_directory_sync() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        fixture.advance();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let old_generation = writer
            .committed_checkpoint()
            .unwrap()
            .manifest
            .generation
            .clone();
        old_presentation(root.path());
        drop(writer);

        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        writer.prepare_restart().unwrap();
        let stage = root.path().join(RESTART_OUTPUT_DIRECTORY);
        let staged_log = writer
            .output_root()
            .unwrap()
            .as_path()
            .join("numerical_stability/failed_samples.jsonl");
        fs::create_dir_all(staged_log.parent().unwrap()).unwrap();
        fs::write(&staged_log, b"new diagnostics").unwrap();
        let mut marker_calls = 0;
        let error = fixture
            .commit(&mut writer, |path| {
                if path != stage.join("cleared") {
                    return Ok(());
                }
                marker_calls += 1;
                if !path.is_file() {
                    // Preserve the earlier synthetic pre-cleanup hook; this case fails only after
                    // the real rename has exposed the complete marker, before directory sync.
                    return Ok(());
                }
                assert!(!stage.join("cleared.tmp").exists());
                Err(IntegrationCheckpointError::Io {
                    path: path.to_owned(),
                    source: IoError::other(
                        "interrupted after marker rename, before directory sync",
                    ),
                })
            })
            .unwrap_err();
        assert_eq!(marker_calls, 2);
        assert!(matches!(
            error,
            IntegrationCheckpointError::PostCommit { iteration: 1, .. }
        ));
        let generation = writer
            .committed_checkpoint()
            .unwrap()
            .manifest
            .generation
            .clone();
        assert_eq!(
            fs::read(stage.join("cleared")).unwrap(),
            generation.as_str().as_bytes()
        );
        assert!(staged_log.is_file());
        assert!(!root.path().join("numerical_stability").exists());
        assert!(
            root.path()
                .join("generations")
                .join(old_generation.as_str())
                .is_dir()
        );
        drop(writer);

        // Reacquisition must make the visible marker durable before pruning or promotion.
        let writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        assert_eq!(
            writer.checkpoint().unwrap().unwrap().manifest.generation,
            generation
        );
        assert_eq!(
            fs::read(root.path().join("numerical_stability/failed_samples.jsonl")).unwrap(),
            b"new diagnostics"
        );
        assert_eq!(
            fs::read_dir(root.path().join("generations"))
                .unwrap()
                .count(),
            1
        );
        assert!(
            !root
                .path()
                .join("generations")
                .join(old_generation.as_str())
                .exists()
        );
        assert!(!stage.exists());
        assert!(!root.path().join(COMPLETED_RESTART_DIRECTORY).exists());
    }

    #[test]
    fn checkpoint_restart_recovers_partial_marker_prune_failure_and_partial_promotion() {
        let invalid_generation = "gen-00000000-00000000000000000000000000000000-0000000000000000";
        for boundary in [
            "results",
            "cleared.tmp",
            "generations",
            "numerical_stability",
            COMPLETED_RESTART_DIRECTORY,
        ] {
            let root = CheckpointTestDirectory::new();
            let mut fixture = Fixture::new();
            fixture.advance();
            let mut writer =
                IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
            fixture.commit(&mut writer, |_| Ok(())).unwrap();
            let old_generation = writer
                .committed_checkpoint()
                .unwrap()
                .manifest
                .generation
                .clone();
            old_presentation(root.path());
            drop(writer);
            let mut writer =
                IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
            writer.prepare_restart().unwrap();
            let staged_log = writer
                .output_root()
                .unwrap()
                .as_path()
                .join("numerical_stability/failed_samples.jsonl");
            fs::create_dir_all(staged_log.parent().unwrap()).unwrap();
            fs::write(&staged_log, b"new diagnostics").unwrap();
            let error = fixture
                .commit(&mut writer, |path| {
                    if !path.ends_with(boundary) {
                        return Ok(());
                    }
                    if boundary == "numerical_stability"
                        && !root
                            .path()
                            .join(RESTART_OUTPUT_DIRECTORY)
                            .join("cleared")
                            .is_file()
                    {
                        // This case targets promotion, after old-output removal has finished.
                        return Ok(());
                    }
                    if boundary == "generations" {
                        // Stray names are ignored. A valid generation name containing a file
                        // still forces a real prune validation error after the durable marker.
                        fs::write(path.join(invalid_generation), b"not a generation directory")
                            .unwrap();
                        return Ok(());
                    }
                    if boundary == COMPLETED_RESTART_DIRECTORY {
                        // A stopped recursive removal may already have removed a phase marker.
                        // The retired namespace must never authorize another output cleanup.
                        fs::remove_file(path.join("cleared")).unwrap();
                    }
                    Err(IntegrationCheckpointError::Io {
                        path: path.to_owned(),
                        source: IoError::other("interrupted restart transition"),
                    })
                })
                .unwrap_err();
            assert!(matches!(
                error,
                IntegrationCheckpointError::PostCommit { iteration: 1, .. }
            ));
            assert!(matches!(
                writer.output_root(),
                Err(IntegrationCheckpointError::RecoveryRequired { .. })
            ));
            let generation = writer
                .committed_checkpoint()
                .unwrap()
                .manifest
                .generation
                .clone();
            let stage = root.path().join(RESTART_OUTPUT_DIRECTORY);
            let journal = if boundary == COMPLETED_RESTART_DIRECTORY {
                assert!(!stage.exists());
                root.path().join(COMPLETED_RESTART_DIRECTORY)
            } else {
                stage.clone()
            };
            assert!(journal.join("commit.json").is_file());
            if boundary == "cleared.tmp" {
                let partial = fs::read(stage.join("cleared.tmp")).unwrap();
                assert!(!partial.is_empty());
                assert!(partial.len() < generation.as_str().len());
                assert!(!stage.join("cleared").exists());
                assert!(staged_log.is_file());
            } else if boundary != "results" && boundary != COMPLETED_RESTART_DIRECTORY {
                assert_eq!(
                    fs::read(stage.join("cleared")).unwrap(),
                    generation.as_str().as_bytes()
                );
            }
            if boundary == "generations" {
                assert!(
                    root.path()
                        .join("generations")
                        .join(old_generation.as_str())
                        .is_dir()
                );
                assert!(
                    staged_log.is_file(),
                    "new output must not be promoted before pruning"
                );
                fs::remove_file(
                    root.path()
                        .join(GENERATIONS_DIRECTORY)
                        .join(invalid_generation),
                )
                .unwrap();
            }
            if boundary == "numerical_stability" || boundary == COMPLETED_RESTART_DIRECTORY {
                assert!(!staged_log.exists());
                assert_eq!(
                    fs::read(root.path().join("numerical_stability/failed_samples.jsonl")).unwrap(),
                    b"new diagnostics"
                );
            }
            drop(writer);
            // Merely reopening performs recovery; no additional integration or commit is needed.
            let writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
            assert_eq!(
                writer.checkpoint().unwrap().unwrap().manifest.generation,
                generation
            );
            assert_eq!(
                fs::read(root.path().join("numerical_stability/failed_samples.jsonl")).unwrap(),
                b"new diagnostics"
            );
            assert_eq!(
                fs::read_dir(root.path().join("generations"))
                    .unwrap()
                    .count(),
                1
            );
            assert!(
                !root
                    .path()
                    .join("generations")
                    .join(old_generation.as_str())
                    .exists()
            );
            assert!(!root.path().join("results").exists());
            assert!(!root.path().join("integrands").exists());
            assert!(!root.path().join("interrupted").exists());
            assert!(!stage.exists());
            assert!(!root.path().join(COMPLETED_RESTART_DIRECTORY).exists());
        }
    }

    #[test]
    fn checkpoint_identifies_versioned_and_unversioned_flat_layouts() {
        let root = CheckpointTestDirectory::new();
        let fixture = Fixture::new();
        let path = root.path().join(INTEGRATION_MANIFEST_FILE);
        for version in [None, Some(3)] {
            let mut manifest = serde_json::to_value(fixture.manifest()).unwrap();
            if let Some(version) = version {
                manifest["version"] = version.into();
            } else {
                manifest.as_object_mut().unwrap().remove("version");
            }
            let bytes = serde_json::to_vec(&manifest).unwrap();
            fs::write(&path, &bytes).unwrap();
            let error = IntegrationCheckpoint::open(root.path()).unwrap_err();
            assert_eq!(
                error.to_string(),
                format!(
                    "Integration workspace at {} was created by an older GammaLoop version",
                    root.path().display()
                )
            );
            assert!(matches!(
                error,
                IntegrationCheckpointError::LegacyLayout { path: rejected }
                    if rejected == root.path()
            ));
            assert_eq!(fs::read(&path).unwrap(), bytes);
        }
        fs::write(&path, br#"{"version":3,"slots":[]}"#).unwrap();
        assert!(matches!(
            IntegrationCheckpoint::open(root.path()),
            Err(IntegrationCheckpointError::ManifestDecode { path: rejected, .. })
                if rejected == path
        ));
        assert!(!root.path().join(".snapshot.lock").exists());
    }

    #[test]
    fn checkpoint_errors_display_context_and_preserve_io_sources() {
        let path = PathBuf::from("run/generations/gen-example/integration_state.bin");
        let error = IntegrationCheckpointError::Io {
            path: path.clone(),
            source: IoError::other("checkpoint read failed"),
        };
        assert_eq!(
            error.to_string(),
            format!("Integration checkpoint I/O failed at {}", path.display())
        );
        assert_eq!(
            error.source().unwrap().to_string(),
            "checkpoint read failed"
        );
        assert_eq!(error.recovery_action(), None);

        let legacy = IntegrationCheckpointError::LegacyLayout {
            path: PathBuf::from("old-run"),
        };
        assert_eq!(
            legacy.to_string(),
            "Integration workspace at old-run was created by an \
             older GammaLoop version"
        );
        assert_eq!(
            legacy.recovery_action(),
            Some(IntegrationCheckpointRecovery::FreshWorkspace)
        );
        let version = IntegrationCheckpointError::UnsupportedVersion {
            found: 5,
            expected: 4,
        };
        assert_eq!(
            version.to_string(),
            "Integration workspace version 5 is unsupported (expected 4; \
             zero denotes an unversioned workspace)"
        );
    }

    #[cfg(unix)]
    #[test]
    fn checkpoint_accepts_a_symlinked_workspace_root_for_read_write_and_restart() {
        let directory = CheckpointTestDirectory::new();
        let root = directory.path().join("actual");
        let alias = directory.path().join("alias");
        fs::create_dir(&root).unwrap();
        unix_fs::symlink(&root, &alias).unwrap();
        let mut fixture = Fixture::new();
        fixture.advance();
        let mut writer = IntegrationWorkspace::acquire(&alias, fixture.manifest()).unwrap();
        assert_eq!(writer.root(), fs::canonicalize(&root).unwrap());
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        Fixture::restore(&IntegrationCheckpoint::open(&alias).unwrap()).assert_same(&fixture);
        drop(writer);

        let mut writer = IntegrationWorkspace::acquire(&alias, fixture.manifest()).unwrap();
        writer.prepare_restart().unwrap();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        Fixture::restore(&IntegrationCheckpoint::open(&root).unwrap()).assert_same(&fixture);
        assert!(
            fs::symlink_metadata(&alias)
                .unwrap()
                .file_type()
                .is_symlink()
        );
    }

    #[cfg(unix)]
    #[test]
    fn restart_rejects_unsafe_entries_anywhere_below_the_workspace_before_staging() {
        let root = CheckpointTestDirectory::new();
        let outside = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        fixture.advance();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        drop(writer);
        let manifest = root.path().join(INTEGRATION_MANIFEST_FILE);
        let before = fs::read(&manifest).unwrap();
        let unrelated = root.path().join("unrelated");
        fs::create_dir(&unrelated).unwrap();
        let linked = unrelated.join("external");
        unix_fs::symlink(outside.path(), &linked).unwrap();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        let error = writer.prepare_restart().unwrap_err();
        assert!(matches!(
            error,
            IntegrationCheckpointError::Output(OutputPathError::Symlink(path)) if path == linked
        ));
        assert!(!root.path().join(RESTART_OUTPUT_DIRECTORY).exists());
        assert_eq!(fs::read(&manifest).unwrap(), before);

        fs::remove_file(&linked).unwrap();
        let socket = unrelated.join("socket");
        let listener = UnixListener::bind(&socket).unwrap();
        let error = writer.prepare_restart().unwrap_err();
        assert!(matches!(
            error,
            IntegrationCheckpointError::Output(OutputPathError::SpecialEntry(path))
                if path == socket
        ));
        assert!(!root.path().join(RESTART_OUTPUT_DIRECTORY).exists());
        assert_eq!(fs::read(&manifest).unwrap(), before);
        drop(listener);
        fs::remove_file(socket).unwrap();

        writer.prepare_restart().unwrap();
        unix_fs::symlink(outside.path(), &linked).unwrap();
        assert!(matches!(
            fixture.commit(&mut writer, |_| Ok(())),
            Err(IntegrationCheckpointError::Output(OutputPathError::Symlink(path)))
                if path == linked
        ));
        assert_eq!(fs::read(&manifest).unwrap(), before);
    }

    #[cfg(unix)]
    #[test]
    fn checkpoint_rejects_internal_symlinks_before_reading_or_restart_cleanup() {
        let root = CheckpointTestDirectory::new();
        let outside = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        fixture.advance();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let checkpoint = writer.checkpoint().unwrap().unwrap();
        let state = checkpoint.state_path();
        fs::copy(&state, outside.path().join("state.bin")).unwrap();
        fs::remove_file(&state).unwrap();
        unix_fs::symlink(outside.path().join("state.bin"), &state).unwrap();
        assert!(matches!(
            IntegrationCheckpoint::open(root.path()),
            Err(IntegrationCheckpointError::Output(_))
        ));
        fs::remove_file(&state).unwrap();
        fs::copy(outside.path().join("state.bin"), &state).unwrap();
        drop(writer);

        fs::write(outside.path().join("do-not-touch"), b"outside bytes").unwrap();
        let pointer = fs::read(root.path().join("manifest.json")).unwrap();
        unix_fs::symlink(outside.path(), root.path().join("integrands")).unwrap();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        assert!(matches!(
            writer.prepare_restart(),
            Err(IntegrationCheckpointError::Output(OutputPathError::Symlink(path)))
                if path == root.path().join("integrands")
        ));
        assert!(!root.path().join(RESTART_OUTPUT_DIRECTORY).exists());
        assert_eq!(
            fs::read(root.path().join("manifest.json")).unwrap(),
            pointer
        );
        assert_eq!(
            fs::read(outside.path().join("do-not-touch")).unwrap(),
            b"outside bytes"
        );
        assert!(
            fs::symlink_metadata(root.path().join("integrands"))
                .unwrap()
                .file_type()
                .is_symlink()
        );
    }
}
