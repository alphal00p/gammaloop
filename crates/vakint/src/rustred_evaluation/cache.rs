//! Explicit lower-loop cache steering. RustRed owns root admission, exact
//! reduction, persistent records, RAM policy and all filesystem operations.
use std::{path::Path, sync::atomic::AtomicBool};

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

#[path = "cache_save.rs"]
mod save;
pub use save::{RustRedCacheSaveMode, RustRedCacheSaveOptions, RustRedCacheSaveReport};

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
}

#[derive(Clone, Debug, Default)]
pub struct RustRedCacheLoadOptions {
    /// None retains the already configured cache's native limits, or the
    /// ordinary native defaults when there is no existing attachment.
    pub native: Option<ReductionLimits>,
    pub policy: Option<CompletedCachePolicy>,
    /// Omitted reload preserves the active policy. A new optional attachment
    /// falls back only on typed prospective lookup-work exhaustion.
    pub lookup_budget_policy: Option<CacheLookupBudgetPolicy>,
    pub storage: SavedApplicationCachePackLimits,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RustRedCacheStatistics {
    pub attached: bool,
    /// Retained native coefficient payload, not process RSS.
    pub native: GuardedRootCacheStatistics,
    pub disk: Option<SavedApplicationCachePackStatistics>,
    pub namespaces: Vec<RustRedCacheNamespaceStatistics>,
    pub backing_budget_fallbacks: usize,
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

fn operation(error: rustred_app::AppError) -> RustRedEvaluationError {
    RustRedEvaluationError::CacheOperation { source: error }
}
fn program(
    loop_count: usize,
) -> Result<(ArtifactFamily, ClosedUnitMassProgram<'static>), RustRedEvaluationError> {
    let family = shipped_family_for_loop_count(loop_count)
        .ok_or(RustRedEvaluationError::CacheBackendUnavailable { loop_count })?;
    let snapshot = family.snapshot()?;
    Ok((
        family,
        ClosedUnitMassProgram::new(snapshot.artifact(), snapshot.source_digest()),
    ))
}

impl RustRedCacheWorkload {
    fn plan(
        &self,
        backend: &impl CachePreparationBackend,
        parents: &[(usize, Vec<bool>)],
        max_roots: usize,
    ) -> Result<CachePreparationWorkload, RustRedEvaluationError> {
        let bindings = backend.bindings().to_vec();
        if let Self::Explicit(roots) = self {
            return CachePreparationWorkload::explicit(
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
        CachePreparationWorkload::parent_pinches(bindings, parents, budget, max_roots)
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
        let resolved = options
            .storage
            .resolve(loop_count, Default::default())
            .map_err(operation)?;
        if loop_count == 4 {
            let mut backend = CandidateCachePreparationBackend::new(
                super::four_loop::cache::programs()?,
                options.native,
            )
            .map_err(operation)?;
            let workload = options.workload.plan(
                &backend,
                &super::four_loop::cache::parent_masks()?,
                resolved.max_roots,
            )?;
            return prepare_application_cache(
                &workload,
                &resolved,
                output.as_ref(),
                cancellation,
                &mut backend,
                observer,
            )
            .map_err(operation)
            .map_err(Into::into);
        }
        let (_, program) = program(loop_count)?;
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
        let workload = options
            .workload
            .plan(&backend, &parents, resolved.max_roots)?;
        prepare_application_cache(
            &workload,
            &resolved,
            output.as_ref(),
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
        if loop_count == 4 {
            let programs = super::four_loop::cache::programs()?;
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
            caches.candidate = Some(loaded);
            return Ok(());
        }
        let (family, program) = program(loop_count)?;
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
        if loop_count == 4 {
            let mut caches = self
                .rustred_lower_caches
                .lock()
                .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
            if let Some(backend) = caches.candidate.as_mut() {
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
            cache.clear_ram().map_err(operation)?;
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
        let (result, decomposition) = backend
            .session_mut(namespace)
            .map_err(operation)?
            .reduce_guarded_root_with_decomposition(target)
            .map_err(operation)?;
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
