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
    CachePreparationBackend, CachePreparationOptions, CachePreparationPrimerBudget,
    CachePreparationReport, CachePreparationWorkload, ClosedApplicationCache,
    ClosedCachePreparationBackend, EntryPowerBudget, SavedApplicationCachePackLimits,
    SavedApplicationCachePackStatistics, prepare_application_cache,
};

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
    pub storage: SavedApplicationCachePackLimits,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RustRedCacheStatistics {
    pub attached: bool,
    /// Retained native coefficient payload, not process RSS.
    pub native: GuardedRootCacheStatistics,
    pub disk: Option<SavedApplicationCachePackStatistics>,
}

#[derive(Debug, Default)]
pub(crate) struct LoadedLowerCaches {
    closed: Vec<(ArtifactFamily, ClosedApplicationCache<'static>)>,
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
        backend: &ClosedCachePreparationBackend<'_>,
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
        // The shipped closed K1/K3/K6 parents have only physical denominator
        // slots. RustRed owns every pinch/composition/count calculation.
        let parents = bindings
            .iter()
            .enumerate()
            .map(|(namespace, binding)| (namespace, vec![true; binding.physical_arity]))
            .collect::<Vec<_>>();
        CachePreparationWorkload::parent_pinches(bindings, &parents, budget, max_roots)
            .map_err(operation)
    }
}

impl Vakint {
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
        let (_, program) = program(loop_count)?;
        let resolved = options
            .storage
            .resolve(loop_count, Default::default())
            .map_err(operation)?;
        let mut backend =
            ClosedCachePreparationBackend::new(vec![program], options.native).map_err(operation)?;
        let workload = options.workload.plan(&backend, resolved.max_roots)?;
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
        })
    }

    /// Clear only attached RAM payloads; preserve the read-only backing and
    /// cumulative work counters. No files are removed.
    pub fn clear_cache(&self, loop_count: usize) -> Result<(), VakintError> {
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
}
