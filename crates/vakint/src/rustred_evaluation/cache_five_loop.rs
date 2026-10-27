//! Saved-application cache steering. Native RustRed owns completed nodes,
//! coefficient transport, record validation and the single-file codec.
use std::{
    collections::BTreeMap,
    fs,
    io::{Read, Write},
    path::{Path, PathBuf},
    sync::{
        Arc, Mutex,
        atomic::{AtomicBool, Ordering},
    },
    time::Instant,
};

use rustred::{
    algebra::{Coefficient, CoefficientContext},
    family::IntegralKey,
    persistence::application_cache::{
        ApplicationCacheBinding, ApplicationCacheKey, ApplicationCacheRecord,
    },
    reduction::{
        GuardedCacheBacking, GuardedCacheBackingError, GuardedCacheBackingErrorKind,
        GuardedCaptureStatistics, GuardedIntermediateCapture,
    },
    solver::{
        CandidateRoutedApplicationObserver, CandidateRoutedPersistentCacheStatistics,
        CandidateStatistics, UnresolvedPolicy,
    },
};
use rustred_app::{
    AppError, AppErrorKind, CachePreparationBackend, CachePreparationReport,
    CachePreparationResolvedOptions, CachePreparationStatus, CachePreparationTarget,
    PreparedSavedApplication, SavedApplicationCacheBundle, SavedApplicationCacheFile,
    SavedApplicationCacheFileLimits, SavedApplicationCacheFileStatistics,
    SavedApplicationCachePack, SavedApplicationEntryCoveragePolicy, SavedApplicationOptions,
    SavedApplicationReadOnlyCache, SavedApplicationReduction, SavedApplicationRunningCache,
    SavedApplicationRunningCacheBase, SavedApplicationRunningCacheControl, load_saved_application,
    prepare_application_cache,
};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

use super::{
    RustRedCacheLoadOptions, RustRedCachePreparationOptions, RustRedCacheStatistics,
    RustRedCacheWorkload, operation, resource_checkpoint,
};
use crate::rustred_evaluation::monitor::MonitorContext;
use crate::rustred_evaluation::{data::RustRedData, five_loop::SavedApplication};
use crate::{Vakint, VakintError, VakintSettings};

pub(super) fn failure(error: impl std::fmt::Display) -> VakintError {
    operation(AppError::new(AppErrorKind::Input, error.to_string())).into()
}

fn backing_error(error: impl std::error::Error + 'static) -> GuardedCacheBackingError {
    if let Some(resource) =
        crate::rustred_evaluation::RustRedEvaluationError::resource_cause(&error)
    {
        return GuardedCacheBackingError::resource(resource.clone());
    }
    GuardedCacheBackingError {
        kind: GuardedCacheBackingErrorKind::Io,
        message: error.to_string(),
    }
}

/// An application identity deliberately excludes cache inventory, release names
/// and numerical master coefficients: none changes a raw Route/Apply record.
pub(in crate::rustred_evaluation) fn saved_application_digest(
    data: &RustRedData,
    path: &Path,
) -> Result<[u8; 32], VakintError> {
    let manifest = if path.is_dir() {
        path.join("application.json")
    } else {
        path.to_owned()
    };
    let application = fs::read(manifest).map_err(failure)?;
    let mut digest = Sha256::new();
    digest.update(b"vakint-saved-application-cache-v1\0");
    digest.update((application.len() as u64).to_le_bytes());
    digest.update(application);
    if let Some(pack) = data.loops.get(&5) {
        fn field(digest: &mut Sha256, bytes: &[u8]) {
            digest.update((bytes.len() as u64).to_le_bytes());
            digest.update(bytes);
        }
        field(&mut digest, pack.manifest.dimension.as_bytes());
        field(&mut digest, pack.manifest.normalization.as_bytes());
        let mut families = pack
            .manifest
            .families
            .iter()
            .map(|family| {
                family["fingerprint"]
                    .as_str()
                    .ok_or_else(|| failure("family fingerprint is absent"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        families.sort_unstable();
        digest.update((families.len() as u64).to_le_bytes());
        for family in families {
            field(&mut digest, family.as_bytes());
        }
        let mut masters = pack
            .manifest
            .masters
            .iter()
            .map(|master| {
                let family = master["family_sha256"]
                    .as_str()
                    .ok_or_else(|| failure("master family identity is absent"))?;
                let powers = master["powers"]
                    .as_array()
                    .ok_or_else(|| failure("master powers are absent"))?
                    .iter()
                    .map(|power| {
                        power
                            .as_i64()
                            .ok_or_else(|| failure("invalid master power"))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                Ok((family, powers))
            })
            .collect::<Result<Vec<_>, VakintError>>()?;
        masters.sort_unstable();
        digest.update((masters.len() as u64).to_le_bytes());
        for (family, powers) in masters {
            field(&mut digest, family.as_bytes());
            digest.update((powers.len() as u64).to_le_bytes());
            for power in powers {
                digest.update(power.to_le_bytes());
            }
        }
    }
    Ok(digest.finalize().into())
}

/// Five-loop work counters have different native semantics from the lower-loop
/// guarded-root cache, so both remain named instead of coercing one into another.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RustRedSavedCacheStatistics {
    pub binding: ApplicationCacheBinding,
    pub read_backing: bool,
    pub native: CandidateStatistics,
    /// Cumulative work of successful requests using the attached backing.
    pub persistent: CandidateRoutedPersistentCacheStatistics,
    pub file: Option<SavedApplicationCacheFileStatistics>,
    pub pack: Option<rustred_app::SavedApplicationCachePackStatistics>,
    /// Native report; enqueued records are explicitly distinct from durable
    /// checkpoints. This is a primitive snapshot, not an engine callback.
    pub running: Option<Value>,
}

impl RustRedSavedCacheStatistics {
    pub fn to_document(&self) -> Value {
        let work = self.persistent;
        json!({
            "binding": {
                "artifact_digest": self.binding.artifact_digest,
                "family_fingerprint": self.binding.family_fingerprint,
                "physical_arity": self.binding.physical_arity,
                "layer": self.binding.layer,
            },
            "read_backing": self.read_backing,
            "native": {
                "rule_applications": self.native.rule_applications(),
                "cache_hits": self.native.cache_hits(),
                "coalescing_additions": self.native.coalescing_additions(),
                "cached_integrals": self.native.cached_integrals(),
                "cached_coefficient_terms": self.native.cached_coefficient_terms(),
                "cached_coefficient_bytes": self.native.cached_coefficient_bytes(),
            },
            "persistent": {
                "lookups": work.lookups, "route_hits": work.route_hits, "apply_hits": work.apply_hits,
                "requested_root_hits": work.requested_root_hits, "intermediate_hits": work.intermediate_hits,
                "misses": work.misses, "captured_records": work.captured_records, "skipped_records": work.skipped_records,
            },
            "file": self.file.map(|file| json!({
                "records": file.records, "compressed_bytes": file.compressed_bytes, "decoded_bytes": file.decoded_bytes,
                "index_owned_bytes": file.index_owned_bytes, "decoded_coefficients": file.decoded_coefficients,
                "decoded_coefficient_terms": file.decoded_coefficient_terms, "decoded_coefficient_bytes": file.decoded_coefficient_bytes,
                "hits": file.hits, "misses": file.misses,
            })),
            "pack": self.pack.map(|pack| json!({
                "runs": pack.runs, "chunks": pack.chunks, "records": pack.records,
                "referenced_chunk_bytes": pack.referenced_chunk_bytes,
                "manifest_bytes": pack.manifest_bytes,
                "directory_owned_bytes": pack.directory_owned_bytes,
                "hits": pack.hits, "misses": pack.misses,
            })),
            "running": self.running,
        })
    }
}

pub(in crate::rustred_evaluation) struct SavedFiveLoopCache {
    pub(super) binding: ApplicationCacheBinding,
    pub(super) backing: SavedFiveLoopBacking,
    pub(super) persistent: CandidateRoutedPersistentCacheStatistics,
}

pub(super) enum SavedFiveLoopBacking {
    File(SavedApplicationCacheFile),
    Pack(Box<CheckpointBacking>),
    ReadOnly(SavedApplicationReadOnlyCache),
    Running {
        session: Arc<Mutex<SavedApplicationRunningCache>>,
        read_backing: bool,
        read_only: Option<SavedApplicationReadOnlyCache>,
    },
}

impl SavedFiveLoopCache {
    pub(super) fn running_control(
        &self,
    ) -> Option<Result<SavedApplicationRunningCacheControl, VakintError>> {
        match &self.backing {
            SavedFiveLoopBacking::Running { session, .. } => Some(
                session
                    .lock()
                    .map(|session| session.control())
                    .map_err(failure),
            ),
            _ => None,
        }
    }

    pub(in crate::rustred_evaluation) fn flush(&self) -> Result<(), VakintError> {
        if let Some(control) = self.running_control() {
            let control = control?;
            let report = control.report().map_err(operation)?;
            if report.closed && report.pending_records == 0 && report.writer_error.is_none() {
                return Ok(());
            }
            control.flush().map_err(operation)?;
        }
        Ok(())
    }

    pub(super) fn immutable_base(&self) -> Result<SavedApplicationRunningCacheBase, VakintError> {
        match &self.backing {
            SavedFiveLoopBacking::File(file) => {
                Ok(SavedApplicationRunningCacheBase::File(file.fork_reader()))
            }
            SavedFiveLoopBacking::Pack(pack) => {
                Ok(SavedApplicationRunningCacheBase::Pack(pack.0.fork_reader()))
            }
            _ => Err(failure(
                "detach the active running or composed snapshot attachment before attaching another writer; supply immutable bases explicitly",
            )),
        }
    }

    fn fork(&self) -> Self {
        let backing = match &self.backing {
            SavedFiveLoopBacking::File(file) => SavedFiveLoopBacking::File(file.fork_reader()),
            SavedFiveLoopBacking::Pack(pack) => {
                SavedFiveLoopBacking::Pack(Box::new(CheckpointBacking(pack.0.fork_reader())))
            }
            SavedFiveLoopBacking::ReadOnly(reader) => {
                SavedFiveLoopBacking::ReadOnly(reader.fork_reader())
            }
            SavedFiveLoopBacking::Running {
                session,
                read_backing,
                read_only,
            } => SavedFiveLoopBacking::Running {
                session: Arc::clone(session),
                read_backing: *read_backing,
                read_only: read_only
                    .as_ref()
                    .map(SavedApplicationReadOnlyCache::fork_reader),
            },
        };
        Self {
            binding: self.binding.clone(),
            backing,
            persistent: Default::default(),
        }
    }

    fn add_work(&mut self, work: CandidateRoutedPersistentCacheStatistics) {
        macro_rules! add { ($($field:ident),* $(,)?) => { $(self.persistent.$field = self.persistent.$field.saturating_add(work.$field);)* }; }
        add!(
            lookups,
            route_hits,
            apply_hits,
            requested_root_hits,
            intermediate_hits,
            misses,
            captured_records,
            skipped_records
        );
    }
}

impl SavedApplication {
    fn reduce_routed_terms_with_cache<'storage>(
        prepared: &mut PreparedSavedApplication<15>,
        terms: Vec<(IntegralKey, Coefficient)>,
        binding: &ApplicationCacheBinding,
        backing: Option<&'storage mut (dyn GuardedCacheBacking + 'storage)>,
        capture: Option<&'storage mut (dyn GuardedIntermediateCapture + 'storage)>,
        observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Result<SavedApplicationReduction, AppError> {
        match observation {
            Some(observation) => prepared
                .reduce_routed_linear_combination_with_cache_and_observation(
                    terms,
                    binding,
                    backing,
                    capture,
                    observation,
                    None,
                ),
            None => prepared
                .reduce_routed_linear_combination_with_cache(terms, binding, backing, capture),
        }
    }

    /// Operational settings may rebuild a reducer without changing the exact
    /// program. Keep only a semantically compatible immutable backing; native
    /// resident results and request counters start fresh in the new session.
    pub(in crate::rustred_evaluation) fn inherit_compatible_cache(
        &self,
        previous: &Self,
    ) -> Result<(), VakintError> {
        let prepared = self.prepared.lock().map_err(failure)?;
        let binding = prepared.cache_binding(self.cache_digest);
        let previous = previous.cache.lock().map_err(failure)?;
        if let Some(previous) = previous.as_ref()
            && previous.binding == binding
        {
            *self.cache.lock().map_err(failure)? = Some(previous.fork());
        }
        Ok(())
    }

    pub(in crate::rustred_evaluation) fn reduce_saved(
        &self,
        prepared: &mut PreparedSavedApplication<15>,
        terms: Vec<(IntegralKey, Coefficient)>,
        observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Result<SavedApplicationReduction, VakintError> {
        if log::log_enabled!(target: "vakint::cache_workload", log::Level::Debug) {
            let roots = terms
                .iter()
                .filter(|(_, coefficient)| !coefficient.is_zero())
                .map(|(integral, _)| integral.powers())
                .collect::<Vec<_>>();
            log::debug!(target: "vakint::cache_workload", "{}", json!({"phase":"input_roots", "roots":roots}));
        }
        let mut attached = self.cache.lock().map_err(failure)?;
        let Some(cache) = attached.as_mut() else {
            return match observation {
                Some(observation) => prepared.reduce_routed_linear_combination_with_observation(
                    terms,
                    observation,
                    None,
                ),
                None => prepared.reduce_routed_linear_combination(terms),
            }
            .map_err(operation)
            .map_err(Into::into);
        };
        let binding = cache.binding.clone();
        let result = match &mut cache.backing {
            SavedFiveLoopBacking::File(file) => Self::reduce_routed_terms_with_cache(
                prepared,
                terms,
                &binding,
                Some(file),
                None,
                observation,
            ),
            SavedFiveLoopBacking::Pack(pack) => Self::reduce_routed_terms_with_cache(
                prepared,
                terms,
                &binding,
                Some(pack.as_mut()),
                None,
                observation,
            ),
            SavedFiveLoopBacking::ReadOnly(reader) => Self::reduce_routed_terms_with_cache(
                prepared,
                terms,
                &binding,
                Some(reader),
                None,
                observation,
            ),
            SavedFiveLoopBacking::Running {
                session,
                read_backing,
                read_only,
            } => {
                let mut session = session.lock().map_err(failure)?;
                if *read_backing && read_only.is_none() {
                    session.with_io(|backing, capture| {
                        Self::reduce_routed_terms_with_cache(
                            prepared,
                            terms,
                            &binding,
                            Some(backing),
                            Some(capture),
                            observation,
                        )
                    })
                } else {
                    session.with_capture(|capture| {
                        Self::reduce_routed_terms_with_cache(
                            prepared,
                            terms,
                            &binding,
                            read_only
                                .as_mut()
                                .map(|reader| reader as &mut dyn GuardedCacheBacking),
                            Some(capture),
                            observation,
                        )
                    })
                }
            }
        }
        .map_err(operation)?;
        cache.add_work(result.statistics().persistent_cache);
        Ok(result)
    }
}

impl Vakint {
    pub(super) fn five_loop_cache_application(
        &self,
    ) -> Result<std::sync::Arc<SavedApplication>, VakintError> {
        if self.rustred_application.lock().map_err(failure)?.is_some() {
            self.loaded_application()
        } else {
            self.saved_application(&VakintSettings::default())
        }
    }

    pub(super) fn load_five_loop_cache(
        &self,
        path: &Path,
        options: &RustRedCacheLoadOptions,
    ) -> Result<(), VakintError> {
        if options.lookup_budget_policy.is_some() {
            return Err(failure(
                "lookup_budget_policy applies to lower-loop chunk packs; the five-loop single-file cache has no per-lookup I/O budget",
            ));
        }
        let application = self.five_loop_cache_application()?;
        // Keep the runtime settings/session identity intact. A cache attachment
        // is not permission to silently replace the active arithmetic limits.
        {
            let loaded = self.rustred_application.lock().map_err(failure)?;
            let configured = &loaded
                .as_ref()
                .ok_or_else(|| failure("saved application disappeared"))?
                .options;
            if options
                .native
                .is_some_and(|limits| limits != configured.reduction_limits)
                || options
                    .policy
                    .is_some_and(|policy| policy != configured.cache_policy)
            {
                return Err(failure(
                    "configure five-loop native limits and cache policy through VakintSettings.rustred_application_options before attaching the cache; load_cache preserves the active settings",
                ));
            }
        }
        let mut prepared = application.prepared.lock().map_err(failure)?;
        let binding = prepared.cache_binding(application.cache_digest);
        let backing = if path.is_dir() {
            if options.file_identity.is_some() {
                return Err(failure(
                    "a single-file identity cannot authenticate a native pack directory",
                ));
            }
            match super::running::RustRedRunningCacheOptions::open_base(
                path,
                &binding,
                options.storage,
            )? {
                SavedApplicationRunningCacheBase::File(file) => SavedFiveLoopBacking::File(file),
                SavedApplicationRunningCacheBase::Pack(pack) => {
                    SavedFiveLoopBacking::Pack(Box::new(CheckpointBacking(pack)))
                }
            }
        } else {
            let limits = file_limits(options.storage)?;
            let file = if let Some(identity) = options.file_identity {
                SavedApplicationCacheFile::open_generated(path, &binding, &identity, limits)
            } else {
                SavedApplicationCacheFile::open(path, &binding, limits)
            }
            .map_err(operation)?;
            let mut validation = file.fork_reader();
            for key in file.keys() {
                let record = validation
                    .lookup(key, prepared.family().coefficient_context())
                    .map_err(operation)?
                    .ok_or_else(|| failure("cache index references an absent record"))?;
                prepared
                    .validate_completed_cache_record(&binding, &record)
                    .map_err(operation)?;
            }
            SavedFiveLoopBacking::File(file)
        };
        // No old attachment or resident result is discarded until the complete
        // framing, binding and index have been checked successfully.
        let mut attached = application.cache.lock().map_err(failure)?;
        if let Some(previous) = attached.as_ref() {
            previous.flush()?;
        }
        // All fallible resource work precedes the replacement commit.
        resource_checkpoint()?;
        prepared.clear_cache().map_err(operation)?;
        *attached = Some(SavedFiveLoopCache {
            binding,
            backing,
            persistent: Default::default(),
        });
        Ok(())
    }

    pub(super) fn five_loop_cache_statistics(&self) -> Result<RustRedCacheStatistics, VakintError> {
        let application = self.five_loop_cache_application()?;
        let prepared = application.prepared.lock().map_err(failure)?;
        let attached = application.cache.lock().map_err(failure)?;
        Ok(RustRedCacheStatistics {
            attached: attached.is_some(),
            saved_application: Some(RustRedSavedCacheStatistics {
                binding: prepared.cache_binding(application.cache_digest),
                read_backing: attached.as_ref().is_some_and(|cache| match &cache.backing {
                    SavedFiveLoopBacking::Running { read_backing, .. } => *read_backing,
                    _ => true,
                }),
                native: prepared.cache_statistics(),
                persistent: attached
                    .as_ref()
                    .map_or(Default::default(), |cache| cache.persistent),
                file: match attached.as_ref().map(|cache| &cache.backing) {
                    Some(SavedFiveLoopBacking::File(file)) => {
                        Some(file.statistics().map_err(operation)?)
                    }
                    _ => None,
                },
                pack: match attached.as_ref().map(|cache| &cache.backing) {
                    Some(SavedFiveLoopBacking::Pack(pack)) => Some(pack.0.statistics()),
                    _ => None,
                },
                running: match attached.as_ref().map(|cache| &cache.backing) {
                    Some(SavedFiveLoopBacking::Running {
                        session,
                        read_backing,
                        read_only,
                    }) => {
                        let report = session
                            .lock()
                            .map_err(failure)?
                            .control()
                            .report()
                            .map_err(operation)?;
                        Some(
                            json!({"read_backing": read_backing, "read_delta": read_only.is_none() && *read_backing, "native": report.to_document()}),
                        )
                    }
                    _ => None,
                },
            }),
            ..Default::default()
        })
    }

    pub(super) fn save_five_loop_cache(
        &self,
        output: &Path,
        mode: super::RustRedCacheSaveMode,
        options: &super::RustRedCacheSaveOptions,
    ) -> Result<super::RustRedCacheSaveReport, VakintError> {
        resource_checkpoint()?;
        if options.resume {
            return Err(failure(
                "a single-file five-loop snapshot cannot be appended; use a new output path, or resume prepare_cache through its checkpoint directory",
            ));
        }
        let application = self.five_loop_cache_application()?;
        let prepared = application.prepared.lock().map_err(failure)?;
        let attached = application.cache.lock().map_err(failure)?;
        let binding = prepared.cache_binding(application.cache_digest);
        let limits = options.resolve(5, Default::default()).map_err(operation)?;
        let file_limits = file_limits(limits)?;
        let mut sink = Snapshot::new(binding, file_limits);
        let context = prepared.family().coefficient_context();
        if matches!(mode, super::RustRedCacheSaveMode::FullSnapshot)
            && let Some(cache) = attached.as_ref()
        {
            let mut accept = |record: ApplicationCacheRecord| -> Result<(), AppError> {
                prepared.validate_completed_cache_record(&sink.binding, &record)?;
                sink.store_completed(&record, context).map_err(|error| {
                    if let GuardedCacheBackingErrorKind::Resource(resource) = &error.kind {
                        AppError::resource(resource.clone())
                    } else {
                        AppError::new(AppErrorKind::Execution, error.to_string())
                    }
                })
            };
            match &cache.backing {
                SavedFiveLoopBacking::File(file) => {
                    let mut reader = file.fork_reader();
                    for key in file.keys() {
                        let record = reader
                            .lookup(key, context)
                            .map_err(operation)?
                            .ok_or_else(|| failure("cache index references an absent record"))?;
                        accept(record).map_err(operation)?;
                    }
                }
                SavedFiveLoopBacking::Pack(pack) => {
                    pack.0
                        .fork_reader()
                        .visit_completed_records(context, accept)
                        .map_err(operation)?;
                }
                SavedFiveLoopBacking::ReadOnly(_) | SavedFiveLoopBacking::Running { .. } => {
                    return Err(failure(
                        "use freeze_running_cache for a native running snapshot; composed immutable bases remain explicit and are not copied to a legacy single-file export",
                    ));
                }
            }
        }
        sink.exported_terms = 0;
        sink.exported_bytes = 0;
        let binding = sink.binding.clone();
        let exported = prepared
            .export_completed_with_statistics(&binding, &mut sink)
            .map_err(operation)?;
        let accepted = rustred::reduction::GuardedRootExportStatistics {
            exported_records: exported.exported_records,
            exported_coefficient_terms: sink.exported_terms,
            exported_coefficient_bytes: sink.exported_bytes,
        };
        let records = sink.records.into_values().collect::<Vec<_>>();
        let (encoded, _) =
            SavedApplicationCacheFile::encode_generated(&binding, &records, context, file_limits)
                .map_err(operation)?;
        write_atomic(output, &encoded)?;
        Ok(super::RustRedCacheSaveReport {
            complete: true,
            accepted,
            durable_records: records.len(),
            skipped_partial_records: exported.skipped_partial_records,
            charged_disk_bytes: encoded.len() as u64,
            backing_attached: attached.is_some(),
            mode,
            base_copy: None,
            storage_limits: limits,
            failure: None,
            checkpoint_failure: None,
        })
    }

    pub(super) fn prepare_five_loop_cache(
        &self,
        output: &Path,
        options: &RustRedCachePreparationOptions,
        cancellation: &AtomicBool,
        observer: impl FnMut(&CachePreparationReport) -> Result<(), AppError>,
        observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Result<CachePreparationReport, VakintError> {
        let started = Instant::now();
        resource_checkpoint()?;
        if !matches!(options.workload, RustRedCacheWorkload::Explicit(_)) {
            return Err(failure(
                "five-loop cache preparation requires an explicit selected workload; exhaustive primers are not supported",
            ));
        }
        if output.exists() {
            return Err(failure(
                "five-loop cache output already exists; preserve the completed file and choose a new destination",
            ));
        }
        let application = self.five_loop_cache_application()?;
        let (path, mut native) = {
            let loaded = self.rustred_application.lock().map_err(failure)?;
            let current = loaded
                .as_ref()
                .ok_or_else(|| failure("saved application disappeared"))?;
            (current.path.clone(), current.options)
        };
        let mut resolved = options
            .storage
            .resolve(5, Default::default())
            .map_err(operation)?;
        let output_limits = file_limits(resolved.storage_limits)?;
        let data = self.rustred_data()?;
        if saved_application_digest(&data, &path)? != application.cache_digest {
            return Err(failure(
                "saved application changed since its validated load; reload the data before preparing a cache",
            ));
        }
        native.unresolved_policy = UnresolvedPolicy::Reject;
        native.entry_coverage_policy = SavedApplicationEntryCoveragePolicy::Advisory;
        native.reduction_limits = options.native;
        native.reduction_limits.max_cached_coefficient_bytes = resolved.retained_coefficient_bytes;
        native.cache_policy = resolved.cache_policy;
        let prepared = load_saved_application::<15>(&path, native).map_err(operation)?;
        if saved_application_digest(&data, &path)? != application.cache_digest {
            return Err(failure(
                "saved application changed while preparing its strict cache session",
            ));
        }
        let binding = prepared.cache_binding(application.cache_digest);
        let running = {
            let attached = application.cache.lock().map_err(failure)?;
            match attached.as_ref().map(|cache| &cache.backing) {
                Some(SavedFiveLoopBacking::Running {
                    session,
                    read_backing,
                    read_only,
                }) => Some((
                    Arc::clone(session),
                    *read_backing,
                    read_only
                        .as_ref()
                        .map(SavedApplicationReadOnlyCache::fork_reader),
                )),
                _ => None,
            }
        };
        let mut backend = SavedBackend {
            prepared,
            binding: vec![binding],
            native,
            backing: None,
            coverage: Vec::new(),
            observation,
            running,
        };
        let workload = options.workload.plan(
            &backend,
            &[],
            resolved.diagnostics_enabled.then_some(resolved.max_roots),
        )?;
        backend.coverage = workload.targets().iter().map(|root| {
            Ok(json!({
                "namespace": root.namespace(), "powers": root.integral().powers(),
                "coverage_admitted": backend.prepared.entry_is_admitted(root.integral()).map_err(operation)?,
            }))
        }).collect::<Result<Vec<_>, VakintError>>()?;
        let checkpoint = options.checkpoint_path.clone().unwrap_or_else(|| {
            let mut name = output.as_os_str().to_owned();
            name.push(".checkpoints");
            PathBuf::from(name)
        });
        // Strict application import is part of the caller's preparation
        // allowance. It must not eat the supervisor's finalization reserve.
        let setup_elapsed = started.elapsed();
        resolved.time_limit = resolved.time_limit.saturating_sub(setup_elapsed);
        MonitorContext::native_phase(
            backend.observation.as_deref_mut(),
            "strict_cache_roots",
            None,
        );
        let mut report = prepare_application_cache(
            &workload,
            &resolved,
            &checkpoint,
            cancellation,
            &mut backend,
            observer,
        )
        .map_err(operation)?;
        report.options["strict_application_setup_seconds"] = setup_elapsed.as_secs_f64().into();
        report.elapsed = started.elapsed();
        if report.status == CachePreparationStatus::Complete {
            MonitorContext::native_phase(
                backend.observation.as_deref_mut(),
                "cache_compaction_and_validation",
                None,
            );
            if cancellation.load(Ordering::Relaxed) {
                report.status = CachePreparationStatus::Cancelled;
                return Ok(report);
            }
            let bundle = SavedApplicationCacheBundle::open_generated(
                &checkpoint,
                backend.bindings(),
                resolved.storage_limits,
            )
            .map_err(operation)?;
            if cancellation.load(Ordering::Relaxed) {
                report.status = CachePreparationStatus::Cancelled;
                report.elapsed = started.elapsed();
                return Ok(report);
            }
            let mut reader = bundle.open_reader(&backend.binding[0]).map_err(operation)?;
            let limits = output_limits;
            let mut sink = Snapshot::new(backend.binding[0].clone(), limits);
            reader
                .visit_completed_records(backend.context(0).map_err(operation)?, |record| {
                    backend.validate_completed_record(0, &record)?;
                    sink.store_completed(&record, backend.context(0)?)
                        .map_err(|error| {
                            if let GuardedCacheBackingErrorKind::Resource(resource) = &error.kind {
                                AppError::resource(resource.clone())
                            } else {
                                AppError::new(AppErrorKind::Limit, error.to_string())
                            }
                        })
                })
                .map_err(operation)?;
            let records = sink.records.into_values().collect::<Vec<_>>();
            let (bytes, identity) = SavedApplicationCacheFile::encode_generated(
                &backend.binding[0],
                &records,
                backend.context(0).map_err(operation)?,
                limits,
            )
            .map_err(operation)?;
            if cancellation.load(Ordering::Relaxed) {
                report.status = CachePreparationStatus::Cancelled;
                report.elapsed = started.elapsed();
                return Ok(report);
            }
            // Write only a fully encoded, validated file. The explicit output is
            // installed atomically; interruption leaves reusable checkpoints.
            let decoded_sha256 = decoded_sha256(&bytes)?;
            report.options["single_file"] = json!({
                "format": "rustred.saved-application-cache-file.v1",
                "identity": identity, "decoded_sha256": decoded_sha256,
                "binding": {
                    "artifact_digest": backend.binding[0].artifact_digest,
                    "family_fingerprint": backend.binding[0].family_fingerprint,
                    "physical_arity": backend.binding[0].physical_arity,
                    "layer": backend.binding[0].layer,
                },
                "records": records.len(), "path": output, "checkpoint_path": checkpoint,
            });
            write_atomic(output, &bytes)?;
            report.elapsed = started.elapsed();
        }
        Ok(report)
    }
}

pub(super) fn file_limits(
    storage: rustred_app::SavedApplicationCachePackLimits,
) -> Result<SavedApplicationCacheFileLimits, VakintError> {
    // Pack chunk limits describe each checkpoint chunk, whereas the completed
    // five-loop file contains the whole selection. Their defaults differ.
    // Reject custom chunk bounds instead of silently widening caller limits.
    if storage.chunk != rustred_app::SavedApplicationCachePackLimits::default().chunk {
        return Err(failure(
            "custom per-chunk limits are unsupported for five-loop single-file caches; retain the native bounded file defaults and configure reduction RAM separately",
        ));
    }
    let mut limits = SavedApplicationCacheFileLimits::default();
    limits.max_compressed_bytes = limits
        .max_compressed_bytes
        .min(usize::try_from(storage.max_disk_bytes).unwrap_or(usize::MAX));
    limits.chunk.max_index_owned_bytes = limits
        .chunk
        .max_index_owned_bytes
        .min(storage.max_directory_owned_bytes);
    Ok(limits)
}

fn decoded_sha256(encoded: &[u8]) -> Result<String, VakintError> {
    let mut decoder = flate2::read::GzDecoder::new(encoded);
    let mut digest = Sha256::new();
    let mut block = [0; 65536];
    loop {
        resource_checkpoint()?;
        let read = decoder.read(&mut block).map_err(failure)?;
        if read == 0 {
            break;
        }
        digest.update(&block[..read]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

fn write_atomic(output: &Path, bytes: &[u8]) -> Result<(), VakintError> {
    resource_checkpoint()?;
    let parent = output
        .parent()
        .filter(|path| !path.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    static NEXT: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
    let nonce = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let temporary = parent.join(format!(
        ".vakint-cache-{}-{nonce}.pending",
        std::process::id()
    ));
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .map_err(failure)?;
    let result = (|| {
        file.write_all(bytes).map_err(failure)?;
        file.sync_all().map_err(failure)?;
        resource_checkpoint()?;
        // A same-directory hard link publishes atomically without replacing an
        // existing completed cache, even if another writer raced this one.
        fs::hard_link(&temporary, output).map_err(failure)
    })();
    drop(file);
    let cleanup = fs::remove_file(&temporary).map_err(failure);
    result.and(cleanup)?;
    #[cfg(unix)]
    fs::File::open(parent)
        .and_then(|directory| directory.sync_all())
        .map_err(failure)?;
    Ok(())
}

struct SavedBackend<'observation, 'callback> {
    // A stack-local borrow for this explicit preparation call, never an engine cache owner.
    observation: Option<&'observation mut CandidateRoutedApplicationObserver<'callback>>,
    prepared: PreparedSavedApplication<15>,
    binding: Vec<ApplicationCacheBinding>,
    native: SavedApplicationOptions,
    backing: Option<CheckpointBacking>,
    coverage: Vec<Value>,
    running: Option<(
        Arc<Mutex<SavedApplicationRunningCache>>,
        bool,
        Option<SavedApplicationReadOnlyCache>,
    )>,
}

impl SavedBackend<'_, '_> {
    fn reduce_root_terms_with_cache<'storage>(
        prepared: &mut PreparedSavedApplication<15>,
        terms: [(IntegralKey, Coefficient); 1],
        binding: &ApplicationCacheBinding,
        backing: Option<&'storage mut (dyn GuardedCacheBacking + 'storage)>,
        capture: &'storage mut (dyn GuardedIntermediateCapture + 'storage),
        observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Result<SavedApplicationReduction, AppError> {
        match observation {
            Some(observation) => prepared.reduce_linear_combination_with_cache_and_observation(
                terms,
                binding,
                backing,
                Some(capture),
                observation,
                None,
            ),
            None => prepared.reduce_linear_combination_with_cache(
                terms,
                binding,
                backing,
                Some(capture),
            ),
        }
    }
}

pub(super) struct CheckpointBacking(pub(super) SavedApplicationCachePack);
impl GuardedCacheBacking for CheckpointBacking {
    fn binding(&self) -> &ApplicationCacheBinding {
        self.0.binding()
    }
    fn lookup(
        &mut self,
        key: &ApplicationCacheKey,
        context: &CoefficientContext,
    ) -> Result<Option<ApplicationCacheRecord>, GuardedCacheBackingError> {
        self.0.lookup(key, context).map_err(backing_error)
    }
    fn store_completed(
        &mut self,
        _: &ApplicationCacheRecord,
        _: &CoefficientContext,
    ) -> Result<(), GuardedCacheBackingError> {
        Ok(())
    }
}

/// Two existing immutable backing views; the running delta is searched first.
/// Admission and decoding remain native, and this borrowed view retains no data.
struct BackingOverlay<'first, 'second> {
    first: &'first mut dyn GuardedCacheBacking,
    second: Option<&'second mut dyn GuardedCacheBacking>,
}
impl GuardedCacheBacking for BackingOverlay<'_, '_> {
    fn binding(&self) -> &ApplicationCacheBinding {
        self.first.binding()
    }
    fn lookup(
        &mut self,
        key: &ApplicationCacheKey,
        context: &CoefficientContext,
    ) -> Result<Option<ApplicationCacheRecord>, GuardedCacheBackingError> {
        let result = self.first.lookup(key, context)?;
        if result.is_some() {
            return Ok(result);
        }
        match self.second.as_mut() {
            Some(second) => second.lookup(key, context),
            None => Ok(None),
        }
    }
    fn store_completed(
        &mut self,
        _: &ApplicationCacheRecord,
        _: &CoefficientContext,
    ) -> Result<(), GuardedCacheBackingError> {
        Ok(())
    }
}

/// Forward the same completed native record to preparation's existing sink and
/// the cumulative writer. Native calls wants/store consecutively for one node;
/// retain only that decision, never a second record index or coefficient table.
struct CapturePair<'a> {
    first: &'a mut dyn GuardedIntermediateCapture,
    second: &'a mut dyn GuardedIntermediateCapture,
    wanted: Option<(ApplicationCacheKey, bool, bool)>,
}
impl GuardedIntermediateCapture for CapturePair<'_> {
    fn fail_on_materialization_limit(&self) -> bool {
        self.first.fail_on_materialization_limit() || self.second.fail_on_materialization_limit()
    }

    fn capture_backing_hits(&self) -> bool {
        // Strict preparation may need a root imported from its checkpoint.
        // Preserve that required copy even though the optional running sink
        // can consequently publish a duplicate. These counters are record
        // publications, never a claim of globally unique cache keys.
        true
    }

    fn metadata_limit_bytes(&self) -> Option<usize> {
        match (
            self.first.metadata_limit_bytes(),
            self.second.metadata_limit_bytes(),
        ) {
            // Either sink can admit a completed node. Each still enforces its
            // own queue/storage grants; the optional writer cannot narrow the
            // existing preparation sink's native metadata allowance.
            (Some(a), Some(b)) => Some(a.max(b)),
            _ => None,
        }
    }
    fn report_statistics(&mut self, statistics: GuardedCaptureStatistics) {
        self.first.report_statistics(statistics);
        self.second.report_statistics(statistics);
    }
    fn wants(
        &mut self,
        key: &ApplicationCacheKey,
        context: &CoefficientContext,
    ) -> Result<bool, GuardedCacheBackingError> {
        let first = self.first.wants(key, context)?;
        let second = self.second.wants(key, context)?;
        self.wanted = Some((key.clone(), first, second));
        Ok(first || second)
    }
    fn store_completed(
        &mut self,
        record: &ApplicationCacheRecord,
        context: &CoefficientContext,
    ) -> Result<(), GuardedCacheBackingError> {
        let Some((key, first, second)) = self.wanted.take() else {
            return Err(backing_error(AppError::new(
                AppErrorKind::InternalInvariant,
                "capture store has no preceding admission",
            )));
        };
        if key != record.key {
            return Err(backing_error(AppError::new(
                AppErrorKind::InternalInvariant,
                "capture admission key changed before store",
            )));
        }
        let first = if first {
            self.first.store_completed(record, context)
        } else {
            Ok(())
        };
        let second = if second {
            self.second.store_completed(record, context)
        } else {
            Ok(())
        };
        // Attempt both owners even when preparation's sink fails. Its original
        // failure remains primary; the running control retains its own fault.
        first.and(second)
    }
}

/// Retain the selected native root while forwarding bounded optional captures.
/// A native root is not reconstructed from post-substitution master coefficients.
struct RootCapture<'a> {
    key: ApplicationCacheKey,
    root: Option<ApplicationCacheRecord>,
    inner: &'a mut dyn GuardedIntermediateCapture,
}
impl GuardedIntermediateCapture for RootCapture<'_> {
    fn fail_on_materialization_limit(&self) -> bool {
        self.inner.fail_on_materialization_limit()
    }

    fn capture_backing_hits(&self) -> bool {
        self.inner.capture_backing_hits()
    }

    fn metadata_limit_bytes(&self) -> Option<usize> {
        self.inner.metadata_limit_bytes()
    }
    fn report_statistics(&mut self, statistics: GuardedCaptureStatistics) {
        self.inner.report_statistics(statistics);
    }
    fn wants(
        &mut self,
        key: &ApplicationCacheKey,
        context: &CoefficientContext,
    ) -> Result<bool, GuardedCacheBackingError> {
        if key == &self.key {
            Ok(true)
        } else {
            self.inner.wants(key, context)
        }
    }
    fn store_completed(
        &mut self,
        record: &ApplicationCacheRecord,
        context: &CoefficientContext,
    ) -> Result<(), GuardedCacheBackingError> {
        if record.key == self.key {
            self.root = Some(record.try_clone_with_resources().map_err(backing_error)?);
            Ok(())
        } else {
            self.inner.store_completed(record, context)
        }
    }
}

impl CachePreparationBackend for SavedBackend<'_, '_> {
    fn bindings(&self) -> &[ApplicationCacheBinding] {
        &self.binding
    }
    fn context(&self, namespace: usize) -> Result<&CoefficientContext, AppError> {
        if namespace != 0 {
            return Err(AppError::new(
                AppErrorKind::Input,
                "five-loop saved cache has one common-family namespace",
            ));
        }
        Ok(self.prepared.family().coefficient_context())
    }
    fn root_key(&self, target: &CachePreparationTarget) -> ApplicationCacheKey {
        ApplicationCacheKey::route(target.integral())
    }
    fn configure(&mut self, options: &CachePreparationResolvedOptions) -> Result<(), AppError> {
        if self.native.reduction_limits.max_cached_coefficient_bytes
            != options.retained_coefficient_bytes
            || self.native.cache_policy != options.cache_policy
        {
            return Err(AppError::new(
                AppErrorKind::InternalInvariant,
                "preparation settings differ from the isolated strict application",
            ));
        }
        Ok(())
    }
    fn validate_entry(&self, target: &CachePreparationTarget) -> Result<(), AppError> {
        self.context(target.namespace())?;
        self.prepared
            .validate_entry_for_reduction(target.integral())
    }
    fn validate_completed_record(
        &self,
        namespace: usize,
        record: &ApplicationCacheRecord,
    ) -> Result<(), AppError> {
        self.context(namespace)?;
        self.prepared
            .validate_completed_cache_record(&self.binding[0], record)
    }
    fn attach_completed_cache(
        &mut self,
        namespace: usize,
        reader: SavedApplicationCachePack,
    ) -> Result<(), AppError> {
        self.context(namespace)?;
        self.backing = Some(CheckpointBacking(reader));
        Ok(())
    }
    fn refresh_completed_caches(&mut self) -> Result<(), AppError> {
        if let Some(backing) = self.backing.as_mut() {
            backing.0.refresh()?;
        }
        Ok(())
    }
    fn completed_cache_directory_bytes(&self) -> Result<usize, AppError> {
        Ok(self
            .backing
            .as_ref()
            .map_or(0, |backing| backing.0.statistics().directory_owned_bytes))
    }
    fn reduce_guarded_root(
        &mut self,
        _: &CachePreparationTarget,
    ) -> Result<ApplicationCacheRecord, AppError> {
        Err(AppError::new(
            AppErrorKind::InternalInvariant,
            "saved root preparation requires native completed-node capture",
        ))
    }
    fn reduce_guarded_root_with_capture(
        &mut self,
        target: &CachePreparationTarget,
        capture: &mut dyn GuardedIntermediateCapture,
    ) -> Result<ApplicationCacheRecord, AppError> {
        let one = self.prepared.family().coefficient_context().one();
        let mut capture = RootCapture {
            key: self.root_key(target),
            root: None,
            inner: capture,
        };
        let terms = [(target.integral().clone(), one)];
        let backing = self
            .backing
            .as_mut()
            .map(|backing| backing as &mut dyn GuardedCacheBacking);
        let prepared = &mut self.prepared;
        let binding = &self.binding[0];
        let observation = self.observation.as_deref_mut();
        if let Some((session, read_backing, read_only)) = &mut self.running {
            let mut session = session.lock().map_err(|_| {
                AppError::new(
                    AppErrorKind::InternalInvariant,
                    "running cache session lock poisoned",
                )
            })?;
            if *read_backing && read_only.is_none() {
                session.with_io(|running_backing, running_capture| {
                    let mut backing = BackingOverlay {
                        first: running_backing,
                        second: backing,
                    };
                    let mut captures = CapturePair {
                        first: &mut capture,
                        second: running_capture,
                        wanted: None,
                    };
                    Self::reduce_root_terms_with_cache(
                        prepared,
                        terms,
                        binding,
                        Some(&mut backing),
                        &mut captures,
                        observation,
                    )
                })?;
            } else {
                session.with_capture(|running_capture| {
                    let mut captures = CapturePair {
                        first: &mut capture,
                        second: running_capture,
                        wanted: None,
                    };
                    Self::reduce_root_terms_with_cache(
                        prepared,
                        terms,
                        binding,
                        read_only
                            .as_mut()
                            .map(|reader| reader as &mut dyn GuardedCacheBacking),
                        &mut captures,
                        observation,
                    )
                })?;
            }
        } else {
            Self::reduce_root_terms_with_cache(
                prepared,
                terms,
                binding,
                backing,
                &mut capture,
                observation,
            )?;
        }
        let root = capture.root;
        let root = match root {
            Some(root) => Some(root),
            None => self
                .prepared
                .completed_cache_record(&self.binding[0], &self.root_key(target))?,
        };
        root.ok_or_else(|| {
            AppError::new(
                AppErrorKind::InternalInvariant,
                "native preparation did not capture its completed Route root",
            )
        })
    }
    fn operational_profile(&self) -> Value {
        json!({"kind":"saved-application-route-apply", "strict":true,
            "unresolved_policy":"reject", "entry_coverage":"advisory", "roots":self.coverage,
            "required_root_postcondition":"complete-native-declared-terminals",
            "native":format!("{:?}", self.native)})
    }
    fn statistics(&self) -> Value {
        let native = self.prepared.cache_statistics();
        let work = self.prepared.cache_work_statistics();
        let running = self.running.as_ref().map(|(session, read_backing, read_only)| {
            let report = session
                .lock()
                .map_err(|_| "running cache session lock poisoned".to_owned())
                .and_then(|session| {
                    session
                        .control()
                        .report()
                        .map_err(|error| error.to_string())
                });
            match report {
                Ok(report) => json!({"read_backing":read_backing,"read_delta":read_only.is_none() && *read_backing,"native":report.to_document()}),
                Err(error) => json!({"read_backing":read_backing,"report_error":error}),
            }
        });
        json!({
            "native": {"rule_applications":native.rule_applications(), "cache_hits":native.cache_hits(),
                "coalescing_additions":native.coalescing_additions(), "cached_integrals":native.cached_integrals(),
                "cached_coefficient_terms":native.cached_coefficient_terms(), "cached_coefficient_bytes":native.cached_coefficient_bytes()},
            "cache_work": {"evicted_entries":work.evicted_entries, "evicted_coefficient_terms":work.evicted_coefficient_terms,
                "evicted_coefficient_bytes":work.evicted_coefficient_bytes, "recomputed_nodes":work.recomputed_nodes},
            "running": running,
        })
    }
}

/// Bounded shared-table snapshot input. Native export already excludes partial
/// nodes; this sink owns only unique records and explicit resource accounting.
struct Snapshot {
    binding: ApplicationCacheBinding,
    records: BTreeMap<ApplicationCacheKey, ApplicationCacheRecord>,
    limits: SavedApplicationCacheFileLimits,
    coefficient_terms: usize,
    coefficient_bytes: usize,
    index_bytes: usize,
    exported_terms: usize,
    exported_bytes: usize,
}
impl Snapshot {
    fn new(binding: ApplicationCacheBinding, limits: SavedApplicationCacheFileLimits) -> Self {
        Self {
            binding,
            records: BTreeMap::new(),
            limits,
            coefficient_terms: 0,
            coefficient_bytes: 0,
            index_bytes: 0,
            exported_terms: 0,
            exported_bytes: 0,
        }
    }
}
impl GuardedCacheBacking for Snapshot {
    fn binding(&self) -> &ApplicationCacheBinding {
        &self.binding
    }
    fn lookup(
        &mut self,
        key: &ApplicationCacheKey,
        _: &CoefficientContext,
    ) -> Result<Option<ApplicationCacheRecord>, GuardedCacheBackingError> {
        symbolica::utils::resources::checkpoint().map_err(GuardedCacheBackingError::resource)?;
        self.records
            .get(key)
            .map(ApplicationCacheRecord::try_clone_with_resources)
            .transpose()
            .map_err(backing_error)
    }
    fn store_completed(
        &mut self,
        record: &ApplicationCacheRecord,
        _: &CoefficientContext,
    ) -> Result<(), GuardedCacheBackingError> {
        symbolica::utils::resources::checkpoint().map_err(GuardedCacheBackingError::resource)?;
        let (terms, bytes) = record.coefficient_census().map_err(backing_error)?;
        let (previous_terms, previous_bytes) = self
            .records
            .get(&record.key)
            .map(ApplicationCacheRecord::coefficient_census)
            .transpose()
            .map_err(backing_error)?
            .unwrap_or_default();
        let next_terms = self
            .coefficient_terms
            .saturating_sub(previous_terms)
            .saturating_add(terms);
        let next_bytes = self
            .coefficient_bytes
            .saturating_sub(previous_bytes)
            .saturating_add(bytes);
        let new_key = !self.records.contains_key(&record.key);
        let next_index = self.index_bytes.saturating_add(if new_key {
            record
                .key
                .owned_heap_bytes()
                .saturating_add(std::mem::size_of::<ApplicationCacheRecord>())
        } else {
            0
        });
        if (self.limits.chunk.diagnostics_enabled
            && (self.records.len().saturating_add(usize::from(new_key))
                > self.limits.chunk.max_records
                || next_terms > self.limits.chunk.max_decoded_coefficient_terms))
            || next_bytes > self.limits.chunk.max_decoded_coefficient_bytes
            || next_index > self.limits.chunk.max_index_owned_bytes
        {
            return Err(GuardedCacheBackingError {
                kind: GuardedCacheBackingErrorKind::Limit,
                message: "five-loop cache snapshot exceeds its bounded unique-record allowance"
                    .into(),
            });
        }
        self.coefficient_terms = next_terms;
        self.coefficient_bytes = next_bytes;
        self.index_bytes = next_index;
        self.exported_terms = self.exported_terms.saturating_add(terms);
        self.exported_bytes = self.exported_bytes.saturating_add(bytes);
        let owned = record.try_clone_with_resources().map_err(backing_error)?;
        self.records.insert(owned.key.clone(), owned);
        symbolica::utils::resources::checkpoint().map_err(GuardedCacheBackingError::resource)?;
        Ok(())
    }
}

#[cfg(all(test, unix))]
mod tests {
    use super::*;

    #[test]
    #[ignore = "requires VAKINT_TEST_EXTERNAL_DATA with the authenticated five-loop campaign"]
    fn running_cache_captures_without_reads_and_reopens_as_a_pinned_backing() {
        let source =
            std::env::var_os("VAKINT_TEST_EXTERNAL_DATA").expect("set VAKINT_TEST_EXTERNAL_DATA");
        let scratch = tempfile::tempdir().unwrap();
        let delta = scratch.path().join("delta");
        let mut engine = Vakint::new().unwrap();
        engine.load_external_data(source, 5).unwrap();
        engine.detach_cache(5).unwrap();
        let application = engine.loaded_application().unwrap();
        let target = IntegralKey::try_new([0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 5]).unwrap();
        let reduce = |observed: bool| {
            use rustred::solver::{
                CandidateRoutedApplicationEnvelope, CandidateRoutedApplicationEvent,
                CandidateRoutedApplicationStatus, CandidateRoutedObservationOptions,
                CandidateRoutedObserverControl, CandidateRoutedObserverStatus,
            };

            let mut sequence = 0;
            let mut core_finished = 0;
            let mut request_finished = 0;
            let mut completed_unique = 0;
            let mut node_completions = 0;
            let result = {
                let mut callback = |event: &CandidateRoutedApplicationEnvelope<'_>| {
                    assert!(event.sequence > sequence);
                    sequence = event.sequence;
                    if event.event == CandidateRoutedApplicationEvent::CoreFinished {
                        core_finished += 1;
                        completed_unique = event.metrics.completed_unique_operational_nodes;
                        assert!(!event.request_finished);
                        assert!(event.frontier.is_some());
                        assert_eq!(
                            event.snapshot.unwrap().status,
                            CandidateRoutedApplicationStatus::Succeeded
                        );
                    }
                    if matches!(
                        event.event,
                        CandidateRoutedApplicationEvent::NodeCompleted(_)
                    ) {
                        node_completions += 1;
                    }
                    if event.event == CandidateRoutedApplicationEvent::RequestFinished {
                        request_finished += 1;
                        assert!(event.request_finished);
                        assert_eq!(event.status, CandidateRoutedApplicationStatus::Succeeded);
                        assert!(event.snapshot.is_some());
                    }
                    Ok(CandidateRoutedObserverControl::Continue)
                };
                let mut observation = CandidateRoutedApplicationObserver::new(
                    CandidateRoutedObservationOptions {
                        interval: std::time::Duration::ZERO,
                        ..Default::default()
                    },
                    &mut callback,
                );
                let mut prepared = application.prepared.lock().unwrap();
                let one = prepared.family().coefficient_context().one();
                let result = application
                    .reduce_saved(
                        &mut prepared,
                        vec![(target.clone(), one)],
                        observed.then_some(&mut observation),
                    )
                    .unwrap();
                assert_eq!(observation.status(), CandidateRoutedObserverStatus::Active);
                result
            };
            assert_eq!(core_finished, usize::from(observed));
            assert_eq!(request_finished, usize::from(observed));
            (result, completed_unique, node_completions)
        };
        let reference = reduce(false).0;
        let assert_parity = |result: &SavedApplicationReduction| {
            assert_eq!(reference.terms(), result.terms());
            assert_eq!(reference.nonzero_conditions(), result.nonzero_conditions());
            assert_eq!(reference.unresolved(), result.unresolved());
            let families = |reduction: &SavedApplicationReduction| {
                reduction
                    .terms()
                    .keys()
                    .map(|integral| integral.family_fingerprint().to_owned())
                    .collect::<std::collections::BTreeSet<_>>()
            };
            assert_eq!(families(&reference), families(result));
            assert_eq!(
                reference.terminal_substitution_additions(),
                result.terminal_substitution_additions()
            );
        };
        engine.clear_cache(5).unwrap();
        let options = super::super::RustRedRunningCacheOptions {
            read_backing: false,
            ..Default::default()
        };
        let control = engine.attach_running_cache(5, &delta, &options).unwrap();
        let captured = reduce(false).0;
        assert_parity(&captured);
        assert_eq!(reference.terms(), captured.terms());
        assert_eq!(
            reference.nonzero_conditions(),
            captured.nonzero_conditions()
        );
        let report = control.report().unwrap();
        assert!(report.accepted_records > 0);
        assert_eq!(report.accepted_records, report.committed_records);
        assert_eq!(report.pending_records, 0);
        assert!(report.writer_error.is_none());
        let stats = engine
            .cache_statistics(5)
            .unwrap()
            .saved_application
            .unwrap();
        assert_eq!(stats.persistent.lookups, 0);
        assert!(stats.running.is_some());

        // Reducer reloads must share this exact lease, not reopen the writer.
        {
            let attached = application.cache.lock().unwrap();
            let original = attached.as_ref().unwrap();
            let fork = original.fork();
            let SavedFiveLoopBacking::Running {
                session: original, ..
            } = &original.backing
            else {
                panic!("running attachment expected")
            };
            let SavedFiveLoopBacking::Running { session: fork, .. } = &fork.backing else {
                panic!("running attachment expected")
            };
            assert!(Arc::ptr_eq(original, fork));
        }
        let frozen = engine
            .freeze_running_cache(5, scratch.path().join("frozen.rrindex"))
            .unwrap();
        assert_eq!(frozen.delta.generation, report.generation);
        // A second writer is refused without changing the current attachment.
        assert!(
            engine
                .attach_running_cache(5, &delta, &Default::default())
                .is_err()
        );
        assert_eq!(control.report().unwrap().generation, report.generation);
        engine.detach_cache(5).unwrap();
        assert_eq!(control.close().unwrap().pending_records, 0);

        let readonly = super::super::RustRedRunningCacheOptions {
            snapshots: vec![frozen.clone()],
            ..Default::default()
        };
        engine.load_cache_snapshot(5, &readonly).unwrap();
        assert!(engine.running_cache_control(5).is_err());
        let cached = reduce(true).0;
        assert_parity(&cached);
        assert_eq!(reference.terms(), cached.terms());
        assert_eq!(reference.nonzero_conditions(), cached.nonzero_conditions());
        let work = engine
            .cache_statistics(5)
            .unwrap()
            .saved_application
            .unwrap()
            .persistent;
        assert!(work.route_hits + work.apply_hits > 0);
        engine.detach_cache(5).unwrap();

        // A frozen read set still captures into a distinct active writer. It
        // must keep both observation and capture on the with_capture branch.
        let frozen_reads = super::super::RustRedRunningCacheOptions {
            snapshots: vec![frozen.clone()],
            read_delta: false,
            ..Default::default()
        };
        let frozen_writer = engine
            .attach_running_cache(5, scratch.path().join("frozen-read-delta"), &frozen_reads)
            .unwrap();
        let captured = reduce(true).0;
        assert_parity(&captured);
        let work = engine
            .cache_statistics(5)
            .unwrap()
            .saved_application
            .unwrap()
            .persistent;
        assert!(work.requested_root_hits + work.intermediate_hits > 0);
        let frozen_report = frozen_writer.flush().unwrap();
        assert_eq!(frozen_report.accepted_records, 0);
        assert_eq!(frozen_report.committed_records, 0);
        assert_eq!(frozen_report.pending_records, 0);
        assert!(frozen_report.writer_error.is_none());
        engine.detach_cache(5).unwrap();
        frozen_writer.close().unwrap();

        // Exercise the composed observer + read/write session from a fresh
        // delta, then reopen it after releasing every reducer resident entry.
        let monitored_delta = scratch.path().join("monitored-delta");
        let monitored = engine
            .attach_running_cache(5, &monitored_delta, &Default::default())
            .unwrap();
        let (cold, completed_unique, node_completions) = reduce(true);
        assert_parity(&cold);
        assert!(completed_unique > 0);
        assert!(node_completions > 0);
        let durable = monitored.flush().unwrap();
        assert!(durable.accepted_records > 0);
        assert_eq!(durable.accepted_records, durable.committed_records);
        assert_eq!(durable.pending_records, 0);
        assert!(durable.writer_error.is_none());
        engine.detach_cache(5).unwrap();
        monitored.close().unwrap();
        let monitored = engine
            .load_running_cache(5, &monitored_delta, &Default::default())
            .unwrap();
        let warm = reduce(true).0;
        assert_parity(&warm);
        let work = engine
            .cache_statistics(5)
            .unwrap()
            .saved_application
            .unwrap()
            .persistent;
        assert!(work.requested_root_hits + work.intermediate_hits > 0);
        let reopened = monitored.flush().unwrap();
        assert_eq!(reopened.accepted_records, 0);
        assert_eq!(reopened.pending_records, 0);
        assert_eq!(reopened.generation, durable.generation);
        engine.detach_cache(5).unwrap();
        monitored.close().unwrap();

        let resumed = engine
            .load_running_cache(5, &delta, &Default::default())
            .unwrap();
        assert_eq!(
            resumed.report().unwrap().generation,
            frozen.delta.generation
        );
        let cached = reduce(true).0;
        assert_parity(&cached);
        assert_eq!(reference.terms(), cached.terms());
        assert_eq!(reference.nonzero_conditions(), cached.nonzero_conditions());
        let work = engine
            .cache_statistics(5)
            .unwrap()
            .saved_application
            .unwrap()
            .persistent;
        assert!(work.requested_root_hits + work.intermediate_hits > 0);
        let warm = resumed.flush().unwrap();
        assert_eq!(warm.accepted_records, 0);
        assert_eq!(warm.committed_records, 0);
        assert_eq!(warm.pending_records, 0);
        assert_eq!(warm.generation, frozen.delta.generation);
        let warm_snapshot = resumed
            .snapshot(&scratch.path().join("warm.rrindex"))
            .unwrap();
        assert_eq!(
            warm_snapshot.delta.manifest_digest,
            frozen.delta.manifest_digest
        );
        resumed.close().unwrap();
        engine.detach_cache(5).unwrap();
    }

    #[test]
    fn single_file_limits_preserve_global_bounds_and_reject_chunk_overrides() {
        let storage = rustred_app::SavedApplicationCachePackLimits::default();
        let defaults = file_limits(storage).unwrap();
        assert_eq!(defaults.chunk.max_records, 1_000_000);
        assert_eq!(defaults.max_compressed_bytes, 95 * 1024 * 1024);
        let mut bounded = storage;
        bounded.max_disk_bytes = 1234;
        bounded.max_directory_owned_bytes = 4321;
        let limits = file_limits(bounded).unwrap();
        assert_eq!(limits.max_compressed_bytes, 1234);
        assert_eq!(limits.chunk.max_index_owned_bytes, 4321);
        bounded = storage;
        bounded.chunk.max_records = 0;
        assert!(file_limits(bounded).is_err());
        bounded = storage;
        bounded.chunk.max_decoded_coefficient_bytes = 0;
        assert!(file_limits(bounded).is_err());
    }

    #[test]
    fn atomic_publication_refuses_resources_before_installing_the_cache() {
        use symbolica::utils::resources::{ResourceBudget, ResourceLimitError, ResourceLimits};

        let scratch = tempfile::tempdir().unwrap();
        let output = scratch.path().join("cache.rrcache.gz");
        let directory = scratch.path().to_owned();
        let budget = ResourceBudget::with_memory_monitor(
            ResourceLimits {
                max_process_rss_bytes: Some(1),
                ..Default::default()
            },
            move || {
                let staged = fs::read_dir(&directory)
                    .map_err(|error| error.to_string())?
                    .try_fold(false, |found, entry| {
                        Ok::<_, std::io::Error>(
                            found
                                || entry?
                                    .file_name()
                                    .to_string_lossy()
                                    .starts_with(".vakint-cache-"),
                        )
                    })
                    .map_err(|error| error.to_string())?;
                Ok(if staged { 2 } else { 0 })
            },
        );
        let error = budget
            .run(|| write_atomic(&output, b"complete"))
            .unwrap_err();
        let mut source: Option<&(dyn std::error::Error + 'static)> = Some(&error);
        let mut resource = None;
        while let Some(error) = source {
            if let Some(error) = error.downcast_ref::<ResourceLimitError>() {
                resource = Some(error);
                break;
            }
            source = error.source();
        }
        assert!(matches!(resource, Some(ResourceLimitError::Memory { .. })));
        assert!(!output.exists());
        assert_eq!(fs::read_dir(scratch.path()).unwrap().count(), 0);
        write_atomic(&output, b"complete").unwrap();
        assert_eq!(fs::read(&output).unwrap(), b"complete");
    }

    #[test]
    #[ignore = "requires VAKINT_TEST_EXTERNAL_DATA with the authenticated five-loop campaign"]
    fn five_loop_single_file_cache_preserves_results_and_transactional_attachment() {
        let source =
            std::env::var_os("VAKINT_TEST_EXTERNAL_DATA").expect("set VAKINT_TEST_EXTERNAL_DATA");
        let scratch = tempfile::tempdir().unwrap();
        let output = scratch.path().join("prepared.rrcache.gz");
        let mut engine = Vakint::new().unwrap();
        engine.load_external_data(source, 5).unwrap();
        engine.detach_cache(5).unwrap();
        let application = engine.loaded_application().unwrap();
        let base = IntegralKey::try_new([0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 5]).unwrap();
        let roots = application
            .prepared
            .lock()
            .unwrap()
            .route_entry_to_owner(&base)
            .unwrap()
            .terms()
            .keys()
            .cloned()
            .map(|integral| super::super::RustRedCacheRoot {
                namespace: 0,
                integral,
            })
            .collect::<Vec<_>>();
        assert!(!roots.is_empty());
        let target = roots[0].integral.clone();
        let options = RustRedCachePreparationOptions {
            workload: RustRedCacheWorkload::Explicit(roots),
            checkpoint_path: Some(scratch.path().join("checkpoints")),
            ..Default::default()
        };
        let report = engine
            .prepare_cache(5, &output, &options, &AtomicBool::new(false), |_| Ok(()))
            .unwrap();
        assert_eq!(
            report.status,
            CachePreparationStatus::Complete,
            "{}",
            report.to_document()
        );
        assert_eq!(report.total_roots, report.completed_roots);
        assert!(report.native_profile["strict"].as_bool().unwrap());
        let reference = {
            let mut prepared = application.prepared.lock().unwrap();
            let one = prepared.family().coefficient_context().one();
            application
                .reduce_saved(&mut prepared, vec![(target.clone(), one)], None)
                .unwrap()
        };
        engine.load_cache(5, &output, &Default::default()).unwrap();
        let cached = {
            let mut prepared = application.prepared.lock().unwrap();
            let one = prepared.family().coefficient_context().one();
            application
                .reduce_saved(&mut prepared, vec![(target.clone(), one)], None)
                .unwrap()
        };
        assert_eq!(reference.terms(), cached.terms());
        let guards = |result: &SavedApplicationReduction| {
            result
                .nonzero_conditions()
                .iter()
                .map(Coefficient::to_expression)
                .collect::<std::collections::BTreeSet<_>>()
        };
        assert_eq!(guards(&reference), guards(&cached));
        let statistics = engine.cache_statistics(5).unwrap();
        let saved = statistics.saved_application.unwrap();
        assert!(statistics.attached);
        assert!(saved.persistent.requested_root_hits > 0);
        assert_eq!(cached.statistics().work.rule_applications(), 0);

        // Entry checks pass; the monitor refuses only after the owning commit
        // has acquired Prepared and Cache. Neither clear nor detach may mutate.
        for detach in [false, true] {
            let snapshot = application.clone();
            let budget = symbolica::utils::resources::ResourceBudget::with_memory_monitor(
                symbolica::utils::resources::ResourceLimits {
                    max_process_rss_bytes: Some(1),
                    ..Default::default()
                },
                move || {
                    let committing =
                        snapshot.prepared.try_lock().is_err() && snapshot.cache.try_lock().is_err();
                    Ok(if committing { 2 } else { 0 })
                },
            );
            assert!(
                budget
                    .run(|| if detach {
                        engine.detach_cache(5)
                    } else {
                        engine.clear_cache(5)
                    })
                    .is_err()
            );
            assert_eq!(
                engine
                    .cache_statistics(5)
                    .unwrap()
                    .saved_application
                    .unwrap(),
                saved
            );
        }

        let bad = scratch.path().join("corrupt.rrcache.gz");
        let mut bytes = fs::read(&output).unwrap();
        let index = bytes.len() / 2;
        bytes[index] ^= 1;
        fs::write(&bad, bytes).unwrap();
        assert!(engine.load_cache(5, &bad, &Default::default()).is_err());
        let after = engine
            .cache_statistics(5)
            .unwrap()
            .saved_application
            .unwrap();
        assert_eq!(after, saved);

        let snapshot = scratch.path().join("snapshot.rrcache.gz");
        let report = engine
            .save_cache(
                5,
                &snapshot,
                super::super::RustRedCacheSaveMode::FullSnapshot,
                &Default::default(),
            )
            .unwrap();
        assert!(report.complete, "{}", report.to_document());
        assert!(report.durable_records >= saved.file.unwrap().records);
        engine
            .load_cache(5, &snapshot, &Default::default())
            .unwrap();
        engine.clear_cache(5).unwrap();
        assert!(engine.cache_statistics(5).unwrap().attached);
        engine.detach_cache(5).unwrap();
        assert!(!engine.cache_statistics(5).unwrap().attached);
    }
}
