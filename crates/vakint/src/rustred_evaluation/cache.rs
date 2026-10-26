//! Explicit native cache steering. RustRed owns root admission, exact
//! reduction, persistent records, RAM policy and all filesystem operations.
use std::{
    path::{Path, PathBuf},
    sync::atomic::AtomicBool,
};

use super::monitor::{MonitorContext, RustRedMonitorOptions};
use rustred::solver::CandidateRoutedApplicationObserver;

use rustred::{
    family::IntegralKey,
    reduction::{
        ClosedUnitMassProgram, CompletedCachePolicy, GuardedRootCacheStatistics,
        HomogeneousMasterDecomposition, ReductionLimits,
    },
};
use rustred_app::{
    CacheLookupBudgetPolicy, CachePreparationBackend, CachePreparationOptions,
    CachePreparationPrimerBudget, CachePreparationReport, CachePreparationWorkload,
    CandidateCachePreparationBackend, ClosedApplicationCache, ClosedCachePreparationBackend,
    EntryPowerBudget, SavedApplicationCachePackLimits, SavedApplicationCachePackStatistics,
    prepare_application_cache,
};

#[path = "cache_five_loop.rs"]
mod five;
pub use five::RustRedSavedCacheStatistics;
pub(super) use five::{SavedFiveLoopCache, saved_application_digest};

#[path = "cache_running.rs"]
mod running;
pub use running::{
    RustRedRunningCacheControl, RustRedRunningCacheOptions, RustRedRunningCacheReport,
    RustRedRunningCacheSnapshot,
};

#[path = "cache_save.rs"]
mod save;
pub use save::{RustRedCacheSaveMode, RustRedCacheSaveOptions, RustRedCacheSaveReport};

#[cfg(all(unix, not(target_arch = "wasm32")))]
#[path = "cache_transfer.rs"]
mod transfer;
#[cfg(all(unix, not(target_arch = "wasm32")))]
pub use transfer::{RustRedCacheTransferOptions, RustRedCacheTransferReport};

use super::{
    RustRedEvaluationError,
    artifact::{ArtifactFamily, shipped_family_for_loop_count},
};
use crate::{Vakint, VakintError};

#[derive(Clone, Debug)]
pub struct RustRedCacheRoot {
    /// Namespace in the shipped program bundle. One through three loops each
    /// have one namespace, numbered zero.
    pub namespace: usize,
    pub integral: IntegralKey,
}

#[derive(Clone, Debug, Default)]
pub enum RustRedCacheWorkload {
    /// Native finite operational primer: two extra positive powers above the
    /// surviving unit-power lines and rank at most one. Not physics coverage.
    #[default]
    Primer,
    RelativePrimer {
        max_extra_dots: u32,
        max_numerator_rank: u32,
    },
    Explicit(Vec<RustRedCacheRoot>),
    Finite(EntryPowerBudget),
}

#[derive(Clone, Debug, Default)]
pub struct RustRedCachePreparationOptions {
    pub storage: CachePreparationOptions,
    /// Native arithmetic/work limits; storage selection never raises these.
    pub native: ReductionLimits,
    pub workload: RustRedCacheWorkload,
    /// Five-loop preparation checkpoints remain separate from its final single
    /// gzip file. Omission uses `<output>.checkpoints`; lower loops use output.
    pub checkpoint_path: Option<PathBuf>,
}

#[derive(Clone, Debug, Default)]
pub struct RustRedCacheLoadOptions {
    /// None retains the already configured cache's native limits, or the
    /// ordinary native defaults when there is no existing attachment.
    pub native: Option<ReductionLimits>,
    pub policy: Option<CompletedCachePolicy>,
    /// Omitted reload preserves the active policy. A new optional attachment
    /// falls back only on typed prospective lookup-work exhaustion.
    /// Supported for one through four loops; five loops reject an override.
    pub lookup_budget_policy: Option<CacheLookupBudgetPolicy>,
    /// Five-loop single files support disk and index bounds, with fixed native
    /// file-codec defaults. Custom per-chunk bounds apply only to 1L–4L packs
    /// and are rejected for five loops.
    pub storage: SavedApplicationCachePackLimits,
    /// Authenticated official single-file identity. Explicit trusted-generated
    /// paths may omit it; their embedded digest and semantic binding are checked.
    pub file_identity: Option<rustred_app::SavedApplicationCacheFileIdentity>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RustRedCacheStatistics {
    /// Five-loop Route/Apply counters retain their own native semantics.
    pub saved_application: Option<RustRedSavedCacheStatistics>,
    pub attached: bool,
    /// Retained native coefficient payload, not process RSS.
    pub native: GuardedRootCacheStatistics,
    pub disk: Option<SavedApplicationCachePackStatistics>,
    pub namespaces: Vec<RustRedCacheNamespaceStatistics>,
    pub backing_budget_fallbacks: usize,
}

impl RustRedCacheStatistics {
    /// A capture-only running writer is attached but performs no persistent
    /// lookup. Keep this distinct from both attachment and resident memoization.
    pub fn has_read_backing(&self) -> bool {
        self.attached
            && self
                .saved_application
                .as_ref()
                .is_none_or(|saved| saved.read_backing)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RustRedCacheNamespaceStatistics {
    pub binding: rustred::persistence::application_cache::ApplicationCacheBinding,
    pub native: GuardedRootCacheStatistics,
    pub disk: Option<SavedApplicationCachePackStatistics>,
    pub backing_budget_fallbacks: usize,
}

#[derive(Debug, Default)]
pub(crate) struct LoadedLowerCaches {
    closed: Vec<(ArtifactFamily, ClosedApplicationCache<'static>)>,
    candidate: Option<CandidateCachePreparationBackend<10>>,
}

pub(super) fn operation(error: rustred_app::AppError) -> RustRedEvaluationError {
    RustRedEvaluationError::CacheOperation { source: error }
}

pub(super) fn resource_checkpoint() -> Result<(), RustRedEvaluationError> {
    symbolica::utils::resources::checkpoint()
        .map_err(rustred_app::AppError::resource)
        .map_err(operation)
}
fn program(
    engine: &Vakint,
    loop_count: usize,
) -> Result<(ArtifactFamily, ClosedUnitMassProgram<'static>), RustRedEvaluationError> {
    let family = shipped_family_for_loop_count(loop_count)
        .ok_or(RustRedEvaluationError::CacheBackendUnavailable { loop_count })?;
    let snapshot = engine.rustred_data()?.artifact(family)?;
    Ok((family, snapshot.program()))
}

impl RustRedCacheWorkload {
    fn plan(
        &self,
        backend: &impl CachePreparationBackend,
        parents: &[(usize, Vec<bool>)],
        max_roots: Option<usize>,
    ) -> Result<CachePreparationWorkload, RustRedEvaluationError> {
        let bindings = backend.bindings().to_vec();
        if let Self::Explicit(roots) = self {
            return CachePreparationWorkload::explicit_with_diagnostic_limit(
                bindings,
                roots
                    .iter()
                    .map(|root| (root.namespace, root.integral.clone())),
                max_roots,
            )
            .map_err(operation);
        }
        let budget = match self {
            Self::Primer => CachePreparationPrimerBudget::default(),
            Self::RelativePrimer {
                max_extra_dots,
                max_numerator_rank,
            } => CachePreparationPrimerBudget::Relative {
                max_extra_dots: *max_extra_dots,
                max_numerator_rank: *max_numerator_rank,
            },
            Self::Finite(budget) => CachePreparationPrimerBudget::Absolute(*budget),
            Self::Explicit(_) => unreachable!(),
        };
        CachePreparationWorkload::parent_pinches_with_diagnostic_limit(
            bindings, parents, budget, max_roots,
        )
        .map_err(operation)
    }
}

impl Vakint {
    /// Used by the Python facade to apply a partial operational override
    /// without resetting unrelated active limits. This does not load programs.
    #[cfg(feature = "symbolica_community_module")]
    pub(crate) fn configured_cache_limits(
        &self,
        loop_count: usize,
    ) -> Result<Option<ReductionLimits>, VakintError> {
        if loop_count == 5 {
            let loaded = self
                .rustred_application
                .lock()
                .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
            return Ok(loaded
                .as_ref()
                .map(|current| current.options.reduction_limits));
        }
        let caches = self
            .rustred_lower_caches
            .lock()
            .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
        if loop_count == 4 {
            return Ok(caches.candidate.as_ref().map(|cache| cache.native_limits()));
        }
        let family = shipped_family_for_loop_count(loop_count)
            .ok_or(RustRedEvaluationError::CacheBackendUnavailable { loop_count })?;
        Ok(caches
            .closed
            .iter()
            .find(|(candidate, _)| *candidate == family)
            .map(|(_, cache)| cache.limits()))
    }
    /// Prepare an explicitly requested exact unit-mass cache. No cache is
    /// attached automatically and no numerical master values are computed.
    /// A partial report must be inspected; it is not successful completion.
    pub fn prepare_cache(
        &self,
        loop_count: usize,
        output: impl AsRef<Path>,
        options: &RustRedCachePreparationOptions,
        cancellation: &AtomicBool,
        observer: impl FnMut(&CachePreparationReport) -> Result<(), rustred_app::AppError>,
    ) -> Result<CachePreparationReport, VakintError> {
        self.prepare_cache_observed(
            loop_count,
            output.as_ref(),
            options,
            cancellation,
            observer,
            None,
        )
    }

    /// Prepare the same strict cache workload with optional native per-call diagnostics.
    pub fn prepare_cache_with_monitor(
        &self,
        loop_count: usize,
        output: impl AsRef<Path>,
        options: &RustRedCachePreparationOptions,
        cancellation: &AtomicBool,
        observer: impl FnMut(&CachePreparationReport) -> Result<(), rustred_app::AppError>,
        monitor: &RustRedMonitorOptions,
    ) -> Result<CachePreparationReport, VakintError> {
        MonitorContext::run(monitor, |mut observation| {
            let result = self.prepare_cache_observed(
                loop_count,
                output.as_ref(),
                options,
                cancellation,
                observer,
                observation.as_deref_mut(),
            );
            if let Ok(report) = &result {
                MonitorContext::native_phase(
                    observation,
                    match report.status {
                        rustred_app::CachePreparationStatus::Complete => {
                            "cache_preparation_complete"
                        }
                        rustred_app::CachePreparationStatus::Cancelled => {
                            "cache_preparation_cancelled"
                        }
                        _ => "cache_preparation_incomplete",
                    },
                );
            }
            result
        })
    }

    fn prepare_cache_observed(
        &self,
        loop_count: usize,
        output: &Path,
        options: &RustRedCachePreparationOptions,
        cancellation: &AtomicBool,
        observer: impl FnMut(&CachePreparationReport) -> Result<(), rustred_app::AppError>,
        mut observation: Option<&mut CandidateRoutedApplicationObserver<'_>>,
    ) -> Result<CachePreparationReport, VakintError> {
        MonitorContext::native_phase(observation.as_deref_mut(), "cache_preparation_setup");
        resource_checkpoint()?;
        if loop_count == 5 {
            return self.prepare_five_loop_cache(
                output,
                options,
                cancellation,
                observer,
                observation,
            );
        }
        MonitorContext::native_phase(observation, "lower_loop_cache_preparation_unobserved");
        if options.checkpoint_path.is_some() {
            return Err(operation(rustred_app::AppError::new(
                rustred_app::AppErrorKind::Input,
                "a separate checkpoint path is only used by the five-loop single-file cache",
            ))
            .into());
        }
        let resolved = options
            .storage
            .resolve(loop_count, Default::default())
            .map_err(operation)?;
        if loop_count == 4 {
            let mut backend = CandidateCachePreparationBackend::new(
                super::four_loop::cache::programs(self)?,
                options.native,
            )
            .map_err(operation)?;
            let workload = options.workload.plan(
                &backend,
                &super::four_loop::cache::parent_masks(self)?,
                resolved.diagnostics_enabled.then_some(resolved.max_roots),
            )?;
            return prepare_application_cache(
                &workload,
                &resolved,
                output,
                cancellation,
                &mut backend,
                observer,
            )
            .map_err(operation)
            .map_err(Into::into);
        }
        let (_, program) = program(self, loop_count)?;
        let mut backend =
            ClosedCachePreparationBackend::new(vec![program], options.native).map_err(operation)?;
        // Physical-slot membership is shipped data. Native RustRed owns every
        // pinch, composition, duplicate and bound calculation.
        let parents = backend
            .bindings()
            .iter()
            .enumerate()
            .map(|(n, b)| (n, vec![true; b.physical_arity]))
            .collect::<Vec<_>>();
        let workload = options.workload.plan(
            &backend,
            &parents,
            resolved.diagnostics_enabled.then_some(resolved.max_roots),
        )?;
        prepare_application_cache(
            &workload,
            &resolved,
            output,
            cancellation,
            &mut backend,
            observer,
        )
        .map_err(operation)
        .map_err(Into::into)
    }

    /// Explicitly attach a trusted-generated read-only cache bundle. Existing
    /// configured RAM/native limits are preserved unless options override them.
    /// Misses may fill native RAM; ordinary evaluation never modifies the pack.
    pub fn load_cache(
        &self,
        loop_count: usize,
        path: impl AsRef<Path>,
        options: &RustRedCacheLoadOptions,
    ) -> Result<(), VakintError> {
        resource_checkpoint()?;
        if loop_count == 5 {
            return self.load_five_loop_cache(path.as_ref(), options);
        }
        if options.file_identity.is_some() {
            return Err(operation(rustred_app::AppError::new(
                rustred_app::AppErrorKind::Input,
                "single-file identity is only applicable to a five-loop cache",
            ))
            .into());
        }
        if loop_count == 4 {
            let programs = super::four_loop::cache::programs(self)?;
            let mut caches = self
                .rustred_lower_caches
                .lock()
                .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
            let existing = caches.candidate.as_ref();
            let limits = options
                .native
                .unwrap_or_else(|| existing.map(|c| c.native_limits()).unwrap_or_default());
            let policy = options.policy.unwrap_or_else(|| {
                existing
                    .map(|c| c.cache_policy())
                    .unwrap_or(CompletedCachePolicy::EvictLeastRecentlyUsed)
            });
            let lookup = options.lookup_budget_policy.unwrap_or_else(|| {
                existing
                    .map(|c| c.lookup_budget_policy())
                    .unwrap_or_default()
            });
            let mut loaded =
                CandidateCachePreparationBackend::with_policy(programs, limits, policy)
                    .map_err(operation)?;
            loaded.set_lookup_budget_policy(lookup);
            loaded
                .attach_bundle(path.as_ref(), options.storage)
                .map_err(operation)?;
            resource_checkpoint()?;
            caches.candidate = Some(loaded);
            return Ok(());
        }
        let (family, program) = program(self, loop_count)?;
        let mut caches = self
            .rustred_lower_caches
            .lock()
            .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
        let existing = caches
            .closed
            .iter()
            .position(|(candidate, _)| *candidate == family);
        let limits = options.native.unwrap_or_else(|| {
            existing
                .map(|index| caches.closed[index].1.limits())
                .unwrap_or_default()
        });
        let policy = options.policy.unwrap_or_else(|| {
            existing
                .map(|index| caches.closed[index].1.cache_policy())
                .unwrap_or(CompletedCachePolicy::EvictLeastRecentlyUsed)
        });
        let mut loaded = ClosedApplicationCache::new(program, limits, policy);
        let lookup = options.lookup_budget_policy.unwrap_or_else(|| {
            existing
                .map(|index| caches.closed[index].1.lookup_budget_policy())
                .unwrap_or_default()
        });
        loaded.set_lookup_budget_policy(lookup);
        loaded
            .attach_bundle(path.as_ref(), options.storage)
            .map_err(operation)?;
        resource_checkpoint()?;
        if let Some(index) = existing {
            caches.closed[index].1 = loaded;
        } else {
            caches.closed.push((family, loaded));
        }
        Ok(())
    }

    pub fn cache_statistics(
        &self,
        loop_count: usize,
    ) -> Result<RustRedCacheStatistics, VakintError> {
        if loop_count == 5 {
            return self.five_loop_cache_statistics();
        }
        if loop_count == 4 {
            let caches = self
                .rustred_lower_caches
                .lock()
                .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
            let Some(backend) = caches.candidate.as_ref() else {
                return Ok(Default::default());
            };
            let mut result = RustRedCacheStatistics {
                attached: true,
                ..Default::default()
            };
            for session in backend.sessions() {
                let native = session.cache_statistics();
                result.native.saturating_merge(native);
                let fallbacks = session.backing_budget_fallbacks();
                result.backing_budget_fallbacks =
                    result.backing_budget_fallbacks.saturating_add(fallbacks);
                result.namespaces.push(RustRedCacheNamespaceStatistics {
                    binding: session.program().binding().clone(),
                    native,
                    disk: session.io_statistics(),
                    backing_budget_fallbacks: fallbacks,
                });
            }
            return Ok(result);
        }
        let family = shipped_family_for_loop_count(loop_count)
            .ok_or(RustRedEvaluationError::CacheBackendUnavailable { loop_count })?;
        let caches = self
            .rustred_lower_caches
            .lock()
            .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
        let Some((_, cache)) = caches
            .closed
            .iter()
            .find(|(candidate, _)| *candidate == family)
        else {
            return Ok(RustRedCacheStatistics::default());
        };
        let statistics = cache.cache_statistics();
        let io = cache.io_statistics();
        Ok(RustRedCacheStatistics {
            saved_application: None,
            attached: true,
            native: statistics,
            disk: io,
            namespaces: vec![RustRedCacheNamespaceStatistics {
                binding: cache.program().binding().clone(),
                native: statistics,
                disk: io,
                backing_budget_fallbacks: cache.backing_budget_fallbacks(),
            }],
            backing_budget_fallbacks: cache.backing_budget_fallbacks(),
        })
    }

    /// Clear only attached RAM payloads; preserve the read-only backing and
    /// cumulative work counters. No files are removed.
    pub fn clear_cache(&self, loop_count: usize) -> Result<(), VakintError> {
        resource_checkpoint()?;
        if loop_count == 5 {
            return self.clear_rustred_application_cache();
        }
        if loop_count == 4 {
            let mut caches = self
                .rustred_lower_caches
                .lock()
                .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
            if let Some(backend) = caches.candidate.as_mut() {
                resource_checkpoint()?;
                backend.clear_ram().map_err(operation)?;
            }
            return Ok(());
        }
        let family = shipped_family_for_loop_count(loop_count)
            .ok_or(RustRedEvaluationError::CacheBackendUnavailable { loop_count })?;
        let mut caches = self
            .rustred_lower_caches
            .lock()
            .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
        if let Some((_, cache)) = caches
            .closed
            .iter_mut()
            .find(|(candidate, _)| *candidate == family)
        {
            resource_checkpoint()?;
            cache.clear_ram().map_err(operation)?;
        }
        Ok(())
    }

    /// Remove an attached backing and its resident reduction cache. Immutable
    /// application data remain loaded; no files are modified or removed.
    pub fn detach_cache(&self, loop_count: usize) -> Result<(), VakintError> {
        resource_checkpoint()?;
        if loop_count == 5 {
            return self.loaded_application()?.clear_cache_state(true);
        }
        let mut caches = self
            .rustred_lower_caches
            .lock()
            .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
        resource_checkpoint()?;
        if loop_count == 4 {
            caches.candidate = None;
        } else {
            let family = shipped_family_for_loop_count(loop_count)
                .ok_or(RustRedEvaluationError::CacheBackendUnavailable { loop_count })?;
            caches.closed.retain(|(candidate, _)| *candidate != family);
        }
        Ok(())
    }

    pub(super) fn cached_closed_homogeneous(
        &self,
        family: ArtifactFamily,
        target: &IntegralKey,
    ) -> Result<Option<HomogeneousMasterDecomposition>, RustRedEvaluationError> {
        let mut caches = self
            .rustred_lower_caches
            .lock()
            .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
        let Some((_, cache)) = caches
            .closed
            .iter_mut()
            .find(|(candidate, _)| *candidate == family)
        else {
            return Ok(None);
        };
        cache
            .reduce_with_common_mass_homogeneity(target)
            .map(Some)
            .map_err(operation)
    }

    pub(super) fn cached_candidate_reduction(
        &self,
        namespace: usize,
        target: &IntegralKey,
    ) -> Result<Option<super::experimental::CandidateReduction>, VakintError> {
        let mut caches = self
            .rustred_lower_caches
            .lock()
            .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
        let Some(backend) = caches.candidate.as_mut() else {
            return Ok(None);
        };
        let reduction = backend
            .session_mut(namespace)
            .map_err(operation)?
            .reduce_guarded_root_with_decomposition(target);
        let (result, decomposition) = match reduction {
            Ok(value) => value,
            Err(error)
                if matches!(
                    error.kind(),
                    rustred_app::AppErrorKind::UnresolvedReduction
                        | rustred_app::AppErrorKind::OutsideReductionScope
                ) =>
            {
                return Ok(None);
            }
            Err(error) => return Err(operation(error).into()),
        };
        Ok(Some(super::experimental::CandidateReduction {
            terms: decomposition
                .terms()
                .iter()
                .map(|(key, coefficient)| (key.clone(), coefficient.to_expression()))
                .collect(),
            applied_rules: result.work.rule_applications(),
        }))
    }
}
