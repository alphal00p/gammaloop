//! Attachment of RustRed's cumulative native pack writer. The native owner
//! controls admission, bounded queues, filesystem locking and durable commits.
use std::{path::Path, sync::Arc};

use rustred_app::{
    SavedApplicationCacheBundle, SavedApplicationCacheFile, SavedApplicationCachePack,
    SavedApplicationReadOnlyCache, SavedApplicationRunningCache, SavedApplicationRunningCacheBase,
    SavedApplicationRunningCacheOptions,
};

use super::{SavedFiveLoopCache, five::SavedFiveLoopBacking, operation, resource_checkpoint};
use crate::{RustRedEvaluationError, Vakint, VakintError};

pub use rustred_app::{
    SavedApplicationRunningCacheControl as RustRedRunningCacheControl,
    SavedApplicationRunningCacheReport as RustRedRunningCacheReport,
    SavedApplicationRunningCacheSnapshot as RustRedRunningCacheSnapshot,
};

/// An existing immutable backing may be used with a separate writable delta.
/// Captured records never modify these bases or their frozen native indices.
#[derive(Clone, Debug)]
pub struct RustRedRunningCacheOptions {
    pub native: SavedApplicationRunningCacheOptions,
    pub base_paths: Vec<std::path::PathBuf>,
    pub snapshots: Vec<RustRedRunningCacheSnapshot>,
    /// False captures useful work without any persistent cache lookup. Native
    /// in-process memoization retains its ordinary semantics in either mode.
    pub read_backing: bool,
    /// False reads only the supplied immutable bases, while still capturing
    /// completed records into the separate delta. This keeps a held-out read
    /// snapshot fixed as the same run produces useful records for later use.
    pub read_delta: bool,
}

impl Default for RustRedRunningCacheOptions {
    fn default() -> Self {
        Self {
            native: Default::default(),
            base_paths: Vec::new(),
            snapshots: Vec::new(),
            read_backing: true,
            read_delta: true,
        }
    }
}

impl RustRedRunningCacheOptions {
    fn bases(
        &self,
        binding: &rustred::persistence::application_cache::ApplicationCacheBinding,
    ) -> Result<Vec<SavedApplicationRunningCacheBase>, VakintError> {
        let mut bases = self
            .base_paths
            .iter()
            .map(|path| Self::open_base(path, binding, self.native.storage))
            .collect::<Result<Vec<_>, _>>()?;
        for snapshot in &self.snapshots {
            if &snapshot.binding != binding {
                return Err(super::five::failure("cache snapshot binding mismatch"));
            }
            let reader = SavedApplicationCachePack::open_snapshot(
                &snapshot.directory,
                binding,
                &snapshot.delta,
                self.native.storage,
            )
            .map_err(operation)?;
            bases.push(SavedApplicationRunningCacheBase::Pack(reader));
        }
        Ok(bases)
    }

    pub(super) fn open_base(
        path: &Path,
        binding: &rustred::persistence::application_cache::ApplicationCacheBinding,
        limits: rustred_app::SavedApplicationCachePackLimits,
    ) -> Result<SavedApplicationRunningCacheBase, VakintError> {
        if path.is_dir() {
            if path.join("bundle.rrindex").exists() {
                let bundle = SavedApplicationCacheBundle::open_generated(
                    path,
                    std::slice::from_ref(binding),
                    limits,
                )
                .map_err(operation)?;
                return SavedApplicationRunningCacheBase::from_bundle(&bundle, binding)
                    .map_err(operation)
                    .map_err(Into::into);
            }
            return SavedApplicationCachePack::open_generated(path, binding, limits)
                .map(SavedApplicationRunningCacheBase::Pack)
                .map_err(operation)
                .map_err(Into::into);
        }
        SavedApplicationCacheFile::open(path, binding, super::five::file_limits(limits)?)
            .map(SavedApplicationRunningCacheBase::File)
            .map_err(operation)
            .map_err(Into::into)
    }
}

impl Vakint {
    /// Create or resume one process's cumulative native delta. The returned
    /// control can report, checkpoint or freeze an index during evaluation
    /// without reacquiring the engine's prepared-application mutex.
    /// Control clones deliberately retain the writer lease after detachment;
    /// `control.close()` flushes and closes every shared session. Detachment
    /// flushes without closing an inherited application that is still running.
    pub fn attach_running_cache(
        &self,
        loop_count: usize,
        path: impl AsRef<Path>,
        options: &RustRedRunningCacheOptions,
    ) -> Result<RustRedRunningCacheControl, VakintError> {
        if loop_count != 5 {
            return Err(RustRedEvaluationError::CacheBackendUnavailable { loop_count }.into());
        }
        resource_checkpoint()?;
        let application = self.five_loop_cache_application()?;
        let mut prepared = application.prepared.lock().map_err(super::five::failure)?;
        let binding = prepared.cache_binding(application.cache_digest);
        let mut attached = application.cache.lock().map_err(super::five::failure)?;
        let mut bases = options.bases(&binding)?;
        if bases.is_empty()
            && options.read_backing
            && let Some(previous) = attached.as_ref()
        {
            bases.push(previous.immutable_base()?);
        }
        let context = prepared.family().coefficient_context().clone();
        let read_only = if options.read_backing && !options.read_delta {
            if bases.is_empty() {
                return Err(super::five::failure(
                    "reading only frozen bases requires at least one explicit base or snapshot",
                ));
            }
            Some(
                SavedApplicationReadOnlyCache::new(
                    binding.clone(),
                    context.clone(),
                    bases
                        .iter()
                        .map(SavedApplicationRunningCacheBase::fork_reader)
                        .collect(),
                )
                .map_err(operation)?,
            )
        } else {
            None
        };
        if let Some(previous) = attached.as_ref() {
            previous.flush()?;
        }
        let session = SavedApplicationRunningCache::open(
            binding.clone(),
            context,
            path.as_ref(),
            bases,
            options.native,
        )
        .map_err(operation)?;
        let control = session.control();
        // Everything fallible precedes replacement. A failed opening leaves the
        // previous attachment and its resident results unchanged.
        resource_checkpoint()?;
        prepared.clear_cache().map_err(operation)?;
        *attached = Some(SavedFiveLoopCache {
            binding,
            backing: SavedFiveLoopBacking::Running {
                session: Arc::new(std::sync::Mutex::new(session)),
                read_backing: options.read_backing,
                read_only,
            },
            persistent: Default::default(),
        });
        Ok(control)
    }

    /// Resume an existing cumulative native delta; an absent path is an error.
    pub fn load_running_cache(
        &self,
        loop_count: usize,
        path: impl AsRef<Path>,
        options: &RustRedRunningCacheOptions,
    ) -> Result<RustRedRunningCacheControl, VakintError> {
        if !path.as_ref().is_dir() {
            return Err(super::five::failure(
                "running cache directory does not exist",
            ));
        }
        self.attach_running_cache(loop_count, path, options)
    }

    /// Attach only the specified immutable bases/frozen indices. No writer is
    /// acquired, no records are captured and the snapshot never refreshes.
    pub fn load_cache_snapshot(
        &self,
        loop_count: usize,
        options: &RustRedRunningCacheOptions,
    ) -> Result<(), VakintError> {
        if loop_count != 5 {
            return Err(RustRedEvaluationError::CacheBackendUnavailable { loop_count }.into());
        }
        let application = self.five_loop_cache_application()?;
        let mut prepared = application.prepared.lock().map_err(super::five::failure)?;
        let binding = prepared.cache_binding(application.cache_digest);
        let bases = options.bases(&binding)?;
        if bases.is_empty() || !options.read_backing {
            return Err(super::five::failure(
                "read-only snapshot requires at least one readable base",
            ));
        }
        let reader = SavedApplicationReadOnlyCache::new(
            binding.clone(),
            prepared.family().coefficient_context().clone(),
            bases,
        )
        .map_err(operation)?;
        let mut attached = application.cache.lock().map_err(super::five::failure)?;
        if let Some(previous) = attached.as_ref() {
            previous.flush()?;
        }
        resource_checkpoint()?;
        prepared.clear_cache().map_err(operation)?;
        *attached = Some(SavedFiveLoopCache {
            binding,
            backing: SavedFiveLoopBacking::ReadOnly(reader),
            persistent: Default::default(),
        });
        Ok(())
    }

    /// A durability barrier for this running delta, without copying its bases
    /// or exporting its payload to a legacy single-file cache.
    pub fn save_running_cache(
        &self,
        loop_count: usize,
    ) -> Result<RustRedRunningCacheReport, VakintError> {
        self.running_cache_control(loop_count)?
            .flush()
            .map_err(operation)
            .map_err(Into::into)
    }

    /// Freeze this delta's acknowledged native index. The descriptor does not
    /// copy or freeze its bases: retain their explicit immutable descriptors
    /// alongside it when constructing a reproducible read-only cache snapshot.
    pub fn freeze_running_cache(
        &self,
        loop_count: usize,
        sidecar: impl AsRef<Path>,
    ) -> Result<RustRedRunningCacheSnapshot, VakintError> {
        self.running_cache_control(loop_count)?
            .snapshot(sidecar.as_ref())
            .map_err(operation)
            .map_err(Into::into)
    }

    pub fn running_cache_control(
        &self,
        loop_count: usize,
    ) -> Result<RustRedRunningCacheControl, VakintError> {
        if loop_count != 5 {
            return Err(RustRedEvaluationError::CacheBackendUnavailable { loop_count }.into());
        }
        let application = self.loaded_application()?;
        let attached = application.cache.lock().map_err(super::five::failure)?;
        attached
            .as_ref()
            .and_then(SavedFiveLoopCache::running_control)
            .transpose()?
            .ok_or_else(|| super::five::failure("no running cache is attached"))
    }
}
