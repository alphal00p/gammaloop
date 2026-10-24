//! Explicit resident-root snapshot steering; all I/O remains in RustRed.
use super::*;

pub use rustred_app::{
    CacheSaveMode as RustRedCacheSaveMode, CacheSaveOptions as RustRedCacheSaveOptions,
    CacheSaveReport as RustRedCacheSaveReport,
};

impl Vakint {
    /// Snapshot an explicitly configured native cache without altering its RAM
    /// allowance, eviction policy or arithmetic limits. Inspect the returned
    /// completion/partial-write report. ResidentRoots excludes unloaded bases;
    /// FullSnapshot copies their verified bytes before the resident overlay.
    /// Ordinary evaluation does not implicitly create a cache session.
    pub fn save_cache(
        &self,
        loop_count: usize,
        output: impl AsRef<Path>,
        mode: RustRedCacheSaveMode,
        options: &RustRedCacheSaveOptions,
    ) -> Result<RustRedCacheSaveReport, VakintError> {
        resource_checkpoint()?;
        if loop_count == 5 {
            return self.save_five_loop_cache(output.as_ref(), mode, options);
        }
        let limits = options
            .resolve(loop_count, Default::default())
            .map_err(operation)?;
        let caches = self
            .rustred_lower_caches
            .lock()
            .map_err(|_| RustRedEvaluationError::CacheLockPoisoned)?;
        if loop_count == 4 {
            let Some(backend) = caches.candidate.as_ref() else {
                return Err(operation(rustred_app::AppError::new(
                    rustred_app::AppErrorKind::Input,
                    "no explicitly configured native cache exists for this loop count",
                ))
                .into());
            };
            return backend
                .save_cache(output.as_ref(), mode, limits, options.resume)
                .map_err(operation)
                .map_err(Into::into);
        }
        let family = shipped_family_for_loop_count(loop_count)
            .ok_or(RustRedEvaluationError::CacheBackendUnavailable { loop_count })?;
        let Some((_, cache)) = caches
            .closed
            .iter()
            .find(|(candidate, _)| *candidate == family)
        else {
            return Err(operation(rustred_app::AppError::new(
                rustred_app::AppErrorKind::Input,
                "no explicitly configured native cache exists for this loop count",
            ))
            .into());
        };
        cache
            .save_cache(output.as_ref(), mode, limits, options.resume)
            .map_err(operation)
            .map_err(Into::into)
    }
}
