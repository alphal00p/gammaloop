//! One manifest commits the scalar/grid state, settings, and observables of a completed iteration.
//!
//! Readers pin an immutable generation; writers retain their OS lock throughout a run. Old
//! generations remain available, including during a failed restart. See [`IntegrationWorkspace`].

use super::{IntegrationSlot, IntegrationState, IntegrationWorkspaceManifest, SlotMeta};
use crate::integrate;
use crate::observables::ObservableSnapshotBundle;
use crate::processes::{OutputPath, OutputPathError};
use crate::settings::RuntimeSettings;
use crate::utils::snapshot::{GenerationId, ManifestName, SnapshotError, SnapshotWriter};
use bincode::config;
use bincode::error::{DecodeError, EncodeError};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::io::{Error as IoError, ErrorKind, Write};
use std::path::{Path, PathBuf};
use thiserror::Error;

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
    #[error(transparent)]
    Snapshot(#[from] SnapshotError),
    #[error(transparent)]
    Output(#[from] OutputPathError),
    #[error("Integration checkpoint I/O failed at {path}")]
    Io { path: PathBuf, source: IoError },
    #[error("Invalid integration checkpoint JSON")]
    Json(#[from] serde_json::Error),
    #[error("Invalid integration checkpoint manifest; use --restart for an incompatible workspace")]
    ManifestDecode { source: serde_json::Error },
    #[error("Invalid integration checkpoint settings")]
    SettingsDecode(#[from] toml::de::Error),
    #[error("Could not encode integration checkpoint settings")]
    SettingsEncode(#[from] toml::ser::Error),
    #[error(
        "Could not decode integration checkpoint state; use --restart for an incompatible layout"
    )]
    StateDecode(#[from] DecodeError),
    #[error("Could not encode integration checkpoint state")]
    StateEncode(#[from] EncodeError),
    #[error("Corrupt integration checkpoint: {0}")]
    Inconsistent(String),
}

type CheckpointResult<T> = Result<T, IntegrationCheckpointError>;

/// Version 4 commit envelope. The ordered workspace slots identify every numeric slot directory.
/// Unlike a requested workspace configuration, an envelope always names a real generation.
/// Read it through [`IntegrationCheckpoint::open`], as in the [`IntegrationWorkspace`] example.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct IntegrationCheckpointManifest {
    pub workspace: IntegrationWorkspaceManifest,
    pub generation: GenerationId,
    pub completed_iteration: usize,
    pub num_points: usize,
    pub observable_slots: Vec<usize>,
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
    payload_root: PathBuf,
}

impl IntegrationCheckpoint {
    /// Pin and validate all referenced members before returning a usable checkpoint.
    ///
    /// # Errors
    /// Missing members, unsupported versions, invalid descriptors, and malformed payloads are
    /// corruption errors; only an absent *manifest* can denote a fresh workspace to a writer.
    /// See the [`IntegrationWorkspace`] example.
    pub fn open(root: &Path) -> CheckpointResult<Self> {
        let manifest_path = integrate::workspace_manifest_path(root);
        let bytes = Self::read(&manifest_path)?;
        Self::from_manifest_bytes(root, &bytes)
    }

    fn from_manifest_bytes(root: &Path, bytes: &[u8]) -> CheckpointResult<Self> {
        let manifest: IntegrationCheckpointManifest = serde_json::from_slice(bytes)
            .map_err(|source| IntegrationCheckpointError::ManifestDecode { source })?;
        Self::validate_workspace(&manifest.workspace)?;
        let n_slots = manifest.workspace.slots.len();
        if manifest.completed_iteration == 0
            || manifest.num_points == 0
            || manifest
                .observable_slots
                .iter()
                .any(|&slot| slot >= n_slots)
            || manifest
                .observable_slots
                .windows(2)
                .any(|pair| pair[0] >= pair[1])
        {
            return Err(Self::inconsistent(
                "invalid completed iteration or observable membership",
            ));
        }
        let checkpoint = Self {
            payload_root: root.join("generations").join(manifest.generation.as_str()),
            manifest,
        };
        checkpoint.read_state()?;
        for slot in 0..n_slots {
            let settings = checkpoint.read_slot_settings(slot)?;
            let bundle = checkpoint.read_observables(slot)?;
            Self::validate_observables(&settings, bundle.as_ref())?;
        }
        Ok(checkpoint)
    }

    fn validate_workspace(workspace: &IntegrationWorkspaceManifest) -> CheckpointResult<()> {
        if workspace.version != IntegrationWorkspaceManifest::CURRENT_VERSION {
            return Err(Self::inconsistent(format!(
                "workspace version {} is unsupported (expected {}); use --restart",
                workspace.version,
                IntegrationWorkspaceManifest::CURRENT_VERSION,
            )));
        }
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
            || workspace.training_slot >= n_slots
            || workspace.integrator_settings_slot >= n_slots
            || (!workspace.reference_gaussians.is_empty()
                && workspace.reference_gaussians.len() != n_slots)
        {
            return Err(Self::inconsistent(
                "workspace slot/configuration lengths disagree",
            ));
        }
        Ok(())
    }

    /// The descriptor pinned when this reader opened. See [`IntegrationWorkspace`].
    #[must_use]
    pub fn manifest(&self) -> &IntegrationCheckpointManifest {
        &self.manifest
    }

    /// Scalar moments and trained grids in this generation. See [`IntegrationWorkspace`].
    #[must_use]
    pub fn state_path(&self) -> PathBuf {
        self.payload_root.join("integration_state.bin")
    }

    /// Settings for the ordered slot index. See [`IntegrationWorkspace`].
    #[must_use]
    pub fn slot_settings_path(&self, slot: usize) -> PathBuf {
        self.payload_root
            .join("slots")
            .join(slot.to_string())
            .join("settings.toml")
    }

    /// A path only for declared observable members. See [`IntegrationWorkspace`].
    #[must_use]
    pub fn observable_path(&self, slot: usize) -> Option<PathBuf> {
        self.manifest
            .observable_slots
            .binary_search(&slot)
            .ok()
            .map(|_| {
                self.payload_root
                    .join("slots")
                    .join(slot.to_string())
                    .join("observables.json")
            })
    }

    /// Decode and check scalar/grid identity against the pinned envelope.
    ///
    /// # Errors
    /// Returns an error for missing, malformed, partial, or inconsistent state. See
    /// [`IntegrationWorkspace`] for an example.
    pub fn read_state(&self) -> CheckpointResult<IntegrationState> {
        let bytes = Self::read(&self.state_path())?;
        let (state, consumed): (IntegrationState, usize) =
            bincode::decode_from_slice(&bytes, config::standard())?;
        self.validate_state(&state)?;
        if consumed != bytes.len() {
            return Err(Self::inconsistent("trailing bytes in scalar/grid state"));
        }
        Ok(state)
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
            return Err(Self::inconsistent(
                "scalar/grid state disagrees with its completed-iteration descriptor",
            ));
        }
        Ok(())
    }

    /// Read one member's settings without resolving a newer manifest.
    ///
    /// # Errors
    /// Invalid indices, missing files, and malformed settings fail. See [`IntegrationWorkspace`].
    pub fn read_slot_settings(&self, slot: usize) -> CheckpointResult<RuntimeSettings> {
        if slot >= self.manifest.workspace.slots.len() {
            return Err(Self::inconsistent("settings slot index is out of range"));
        }
        let path = self.slot_settings_path(slot);
        OutputPath::new(&path)?;
        let text = fs::read_to_string(&path)
            .map_err(|source| IntegrationCheckpointError::Io { path, source })?;
        let member: SlotCheckpointMember<GenerationId, SlotMeta, RuntimeSettings> =
            toml::from_str(&text)?;
        self.validate_member(slot, &member)?;
        Ok(member.payload)
    }

    /// Read exactly the bundle declared for this slot; absence is explicit in the envelope.
    ///
    /// # Errors
    /// Missing declared members, undeclared bundles, or invalid indices fail. See
    /// [`IntegrationWorkspace`] for an example.
    pub fn read_observables(
        &self,
        slot: usize,
    ) -> CheckpointResult<Option<ObservableSnapshotBundle>> {
        if slot >= self.manifest.workspace.slots.len() {
            return Err(Self::inconsistent("observable slot index is out of range"));
        }
        if let Some(path) = self.observable_path(slot) {
            let bytes = Self::read(&path)?;
            let member: SlotCheckpointMember<GenerationId, SlotMeta, ObservableSnapshotBundle> =
                serde_json::from_slice(&bytes)?;
            self.validate_member(slot, &member)?;
            return Ok(Some(member.payload));
        }
        let path = self
            .slot_settings_path(slot)
            .with_file_name("observables.json");
        match fs::symlink_metadata(&path) {
            Err(source) if source.kind() == ErrorKind::NotFound => Ok(None),
            Err(source) => Err(IntegrationCheckpointError::Io { path, source }),
            Ok(_) => Err(Self::inconsistent("undeclared observable bundle")),
        }
    }

    fn validate_member<T>(
        &self,
        slot: usize,
        member: &SlotCheckpointMember<GenerationId, SlotMeta, T>,
    ) -> CheckpointResult<()> {
        if member.generation != self.manifest.generation
            || member.completed_iteration != self.manifest.completed_iteration
            || member.num_points != self.manifest.num_points
            || member.slot != self.manifest.workspace.slots[slot]
        {
            return Err(Self::inconsistent(
                "slot member disagrees with its generation/iteration/configuration descriptor",
            ));
        }
        Ok(())
    }

    fn member<'a, T>(
        &'a self,
        slot: usize,
        payload: &'a T,
    ) -> SlotCheckpointMember<&'a GenerationId, &'a SlotMeta, &'a T> {
        SlotCheckpointMember {
            generation: &self.manifest.generation,
            completed_iteration: self.manifest.completed_iteration,
            num_points: self.manifest.num_points,
            slot: &self.manifest.workspace.slots[slot],
            payload,
        }
    }

    fn validate_observables(
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
            _ => Err(Self::inconsistent(
                "observable configuration or membership disagrees with slot settings",
            )),
        }
    }

    fn read(path: &Path) -> CheckpointResult<Vec<u8>> {
        OutputPath::new(path)?;
        let metadata =
            fs::symlink_metadata(path).map_err(|source| IntegrationCheckpointError::Io {
                path: path.to_owned(),
                source,
            })?;
        if !metadata.is_file() {
            return Err(Self::inconsistent(format!(
                "{} is not a regular file",
                path.display()
            )));
        }
        fs::read(path).map_err(|source| IntegrationCheckpointError::Io {
            path: path.to_owned(),
            source,
        })
    }

    fn inconsistent(message: impl Into<String>) -> IntegrationCheckpointError {
        IntegrationCheckpointError::Inconsistent(message.into())
    }
}

/// Owns exclusive writing for a run, including resume, restart preparation, and output assembly.
/// Internal generations are independent of the optional user-facing iteration archives.
///
/// ```no_run
/// # use gammalooprs::integrate::{IntegrationCheckpoint, IntegrationSlot, IntegrationState,
/// #     IntegrationWorkspace, IntegrationWorkspaceManifest};
/// # use std::path::Path;
/// # fn checkpoint(root: &Path, manifest: IntegrationWorkspaceManifest,
/// #     state: &IntegrationState, slots: &[IntegrationSlot]) -> Result<(), Box<dyn std::error::Error>> {
/// let mut workspace = IntegrationWorkspace::acquire(root, manifest)?;
/// // Run workers and merge a complete iteration before publishing it.
/// workspace.commit_completed_iteration(state, slots)?;
/// let pinned = IntegrationCheckpoint::open(workspace.root())?;
/// assert_eq!(pinned.read_state()?.iter, state.iter);
/// let settings = pinned.read_slot_settings(0)?;
/// let observables = pinned.read_observables(0)?;
/// # Ok(())
/// # }
/// ```
#[derive(Debug)]
pub struct IntegrationWorkspace {
    writer: SnapshotWriter,
    pub manifest: IntegrationWorkspaceManifest,
    committed: Option<IntegrationCheckpoint>,
}

impl IntegrationWorkspace {
    /// Acquire before reading for resume or preparing restart. Never removes the old commit.
    ///
    /// # Errors
    /// Writer contention and invalid configuration fail before any payload is written. See the
    /// [type example](Self).
    pub fn acquire(root: &Path, manifest: IntegrationWorkspaceManifest) -> CheckpointResult<Self> {
        IntegrationCheckpoint::validate_workspace(&manifest)?;
        Ok(Self {
            writer: SnapshotWriter::acquire(root)?,
            manifest,
            committed: None,
        })
    }

    /// Logical workspace root, used for presentation outputs. See the [type example](Self).
    #[must_use]
    pub fn root(&self) -> &Path {
        self.writer.root()
    }

    /// Resolve once under the held lock; only an absent manifest means no completed iteration.
    ///
    /// # Errors
    /// A present but corrupt checkpoint fails, including missing referenced members. See the
    /// [type example](Self).
    pub fn checkpoint(&self) -> CheckpointResult<Option<IntegrationCheckpoint>> {
        let path = integrate::workspace_manifest_path(self.root());
        OutputPath::new(&path)?;
        match fs::read(&path) {
            Ok(bytes) => IntegrationCheckpoint::from_manifest_bytes(self.root(), &bytes).map(Some),
            Err(source) if source.kind() == ErrorKind::NotFound => Ok(None),
            Err(source) => Err(IntegrationCheckpointError::Io { path, source }),
        }
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
    /// All failures precede publication and leave the old manifest intact. Zero or partial
    /// iterations and mismatched slots are rejected. See the [type example](Self).
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
            return Err(IntegrationCheckpoint::inconsistent(
                "slot targets disagree with workspace configuration",
            ));
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
        IntegrationCheckpoint::validate_workspace(&self.manifest)?;
        let slots: Vec<_> = slots.into_iter().collect();
        if slots.len() != self.manifest.slots.len()
            || bundles.len() != slots.len()
            || slots
                .iter()
                .map(|(meta, _)| *meta)
                .ne(self.manifest.slots.iter())
            || state.iter == 0
            || state.num_points == 0
        {
            return Err(IntegrationCheckpoint::inconsistent(
                "cannot commit an incomplete iteration or mismatched slots",
            ));
        }
        for ((_, settings), bundle) in slots.iter().zip(bundles) {
            IntegrationCheckpoint::validate_observables(settings, bundle.as_ref())?;
        }
        let root = self.root().to_owned();
        let manifest_path = integrate::workspace_manifest_path(&root);
        let staged = self.writer.begin()?;
        let manifest = IntegrationCheckpointManifest {
            workspace: self.manifest.clone(),
            generation: staged.id().clone(),
            completed_iteration: state.iter,
            num_points: state.num_points,
            observable_slots: bundles
                .iter()
                .enumerate()
                .filter_map(|(slot, bundle)| bundle.as_ref().map(|_| slot))
                .collect(),
        };
        let checkpoint = IntegrationCheckpoint {
            manifest,
            payload_root: staged.path().to_owned(),
        };
        checkpoint.validate_state(state)?;
        Self::write(
            &checkpoint.state_path(),
            &bincode::encode_to_vec(state, config::standard())?,
            &mut before_write,
        )?;
        for (slot, ((_, settings), bundle)) in slots.iter().zip(bundles).enumerate() {
            let path = checkpoint.slot_settings_path(slot);
            fs::create_dir_all(path.parent().expect("slot settings have a parent")).map_err(
                |source| IntegrationCheckpointError::Io {
                    path: path.clone(),
                    source,
                },
            )?;
            Self::write(
                &path,
                toml::to_string_pretty(&checkpoint.member(slot, *settings))?.as_bytes(),
                &mut before_write,
            )?;
            if let Some(bundle) = bundle {
                Self::write(
                    &checkpoint
                        .observable_path(slot)
                        .expect("bundle membership is declared"),
                    &serde_json::to_vec_pretty(&checkpoint.member(slot, bundle))?,
                    &mut before_write,
                )?;
            }
        }
        let bytes = serde_json::to_vec_pretty(&checkpoint.manifest)?;
        // JSON encoding can succeed while replacing nonfinite numbers with null. Validate the
        // complete serialized candidate through the same owner that will load it after commit.
        // No staged file is modified after this boundary.
        let committed = IntegrationCheckpoint::from_manifest_bytes(&root, &bytes)?;
        before_write(&manifest_path)?;
        staged.publish(&bytes, ManifestName::Integration)?;
        self.committed = Some(committed);
        Ok(())
    }

    fn write(
        path: &Path,
        bytes: &[u8],
        before_write: &mut impl FnMut(&Path) -> CheckpointResult<()>,
    ) -> CheckpointResult<()> {
        before_write(path)?;
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
    #[cfg(unix)]
    use std::os::unix::fs::MetadataExt;
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
                training_slot: 0,
                integrator_settings_slot: 0,
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
                    .map(|slot| checkpoint.read_slot_settings(slot).unwrap())
                    .collect(),
                histograms: (0..2)
                    .map(|slot| {
                        checkpoint
                            .read_observables(slot)
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
        assert_eq!(checkpoint.manifest().observable_slots, vec![0, 1]);
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
            assert_eq!(checkpoint.manifest().observable_slots, vec![0, 1]);
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
            }
            Fixture::restore(&IntegrationCheckpoint::open(root.path()).unwrap())
                .assert_same(&uninterrupted);
            // A previously pinned reader keeps its old generation after both later publications.
            assert_eq!(pinned.read_state().unwrap().iter, 1);
        }
    }

    #[test]
    fn checkpoint_missing_member_and_config_mismatch_are_corruption() {
        let root = CheckpointTestDirectory::new();
        let mut fixture = Fixture::new();
        let mut writer = IntegrationWorkspace::acquire(root.path(), fixture.manifest()).unwrap();
        fixture.advance();
        fixture.commit(&mut writer, |_| Ok(())).unwrap();
        let checkpoint = IntegrationCheckpoint::open(root.path()).unwrap();
        let path = checkpoint.observable_path(1).unwrap();
        let bytes = fs::read(&path).unwrap();
        fs::remove_file(&path).unwrap();
        assert!(writer.checkpoint().is_err());
        assert!(IntegrationCheckpoint::open(root.path()).is_err());
        fs::write(&path, &bytes).unwrap();
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
            checkpoint.slot_settings_path(0),
            checkpoint.slot_settings_path(1),
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
        assert_eq!(checkpoint.manifest().observable_slots, vec![0]);
        assert!(checkpoint.read_observables(1).unwrap().is_none());
        let previous = fs::read(root.path().join("manifest.json")).unwrap();
        // An iteration's scalar accumulator has advanced while its descriptor has not.
        fixture.state.all_integrals[0].re.update_iter(false);
        assert!(
            writer
                .commit(
                    &fixture.state,
                    fixture.state.slot_metas.iter().zip(&fixture.settings),
                    &bundles,
                    |_| Ok(())
                )
                .is_err()
        );
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
        fs::copy(
            pinned.state_path(),
            destination.join("integration_state.bin"),
        )
        .unwrap();
        for slot in 0..2 {
            fs::copy(
                pinned.slot_settings_path(slot),
                destination.join(format!("slots/{slot}/settings.toml")),
            )
            .unwrap();
            fs::copy(
                pinned.observable_path(slot).unwrap(),
                destination.join(format!("slots/{slot}/observables.json")),
            )
            .unwrap();
        }
        let read_only = IntegrationCheckpoint::open(target.path()).unwrap();
        Fixture::restore(&read_only).assert_same(&fixture);
        assert!(!target.path().join(".snapshot.lock").exists());
        assert!(!target.path().join("integration_result.json").exists());
    }
}
