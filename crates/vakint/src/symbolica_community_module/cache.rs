//! Thin Python steering of the existing native cache APIs, with cooperative
//! Python callbacks only at the runner's ordinary root/checkpoint boundaries.
use std::{
    path::PathBuf,
    sync::atomic::{AtomicBool, Ordering},
    time::{Duration, Instant},
};

use pyo3::{
    prelude::*,
    types::{PyAnyMethods, PyDict, PyDictMethods, PyList, PyListMethods, PyModule},
};
use rustred::{
    family::IntegralKey,
    reduction::{GuardedRootCacheStatistics, ReductionLimits},
};
use rustred_app::{
    AppError, AppErrorKind, CacheLookupBudgetPolicy, CachePreparationReport,
    CachePreparationStatus, SavedApplicationCachePackStatistics,
};

use super::vakint_to_python_error;
use crate::{
    RustRedCachePreparationOptions, RustRedCacheRoot, RustRedCacheSaveMode, RustRedCacheStatistics,
    RustRedCacheWorkload, Vakint,
};

pub(super) fn check_loops(loops: usize) -> PyResult<()> {
    if !(1..=5).contains(&loops) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "native cache preparation/storage supports one through five loops",
        ));
    }
    Ok(())
}

pub(super) fn lookup_policy(value: &str) -> PyResult<CacheLookupBudgetPolicy> {
    match value {
        "strict" => Ok(CacheLookupBudgetPolicy::Strict),
        "cold-on-work-exhaustion" => Ok(CacheLookupBudgetPolicy::ColdOnLookupWorkExhaustion),
        _ => Err(pyo3::exceptions::PyValueError::new_err(
            "lookup_budget_policy must be 'strict' or 'cold-on-work-exhaustion'",
        )),
    }
}
pub(super) fn save_mode(value: &str) -> PyResult<RustRedCacheSaveMode> {
    match value {
        "resident-roots" => Ok(RustRedCacheSaveMode::ResidentRoots),
        "full-snapshot" => Ok(RustRedCacheSaveMode::FullSnapshot),
        _ => Err(pyo3::exceptions::PyValueError::new_err(
            "mode must be 'resident-roots' or 'full-snapshot'",
        )),
    }
}
pub(super) fn apply_native_overrides(
    limits: &mut ReductionLimits,
    ram_bytes: Option<u64>,
    terms: Option<usize>,
    integrals: Option<usize>,
    rules: Option<usize>,
) -> PyResult<()> {
    if let Some(bytes) = ram_bytes {
        limits.max_cached_coefficient_bytes = usize::try_from(bytes).map_err(|_| {
            pyo3::exceptions::PyValueError::new_err("RAM allowance does not fit this platform")
        })?;
    }
    if let Some(value) = terms {
        limits.max_cached_coefficient_terms = value;
    }
    if let Some(value) = integrals {
        limits.max_cached_integrals = value;
    }
    if let Some(value) = rules {
        limits.max_rule_applications = value;
    }
    Ok(())
}

pub(super) fn workload(
    dots: Option<u32>,
    rank: Option<u32>,
    roots: Option<Vec<(usize, Vec<i64>)>>,
) -> PyResult<RustRedCacheWorkload> {
    match roots {
        Some(roots) => {
            if dots.is_some() || rank.is_some() {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "explicit workload cannot be combined with relative primer controls",
                ));
            }
            roots
                .into_iter()
                .map(|(namespace, powers)| {
                    Ok(RustRedCacheRoot {
                        namespace,
                        integral: IntegralKey::try_new(powers)
                            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?,
                    })
                })
                .collect::<PyResult<Vec<_>>>()
                .map(RustRedCacheWorkload::Explicit)
        }
        None if dots.is_none() && rank.is_none() => Ok(RustRedCacheWorkload::Primer),
        None => Ok(RustRedCacheWorkload::RelativePrimer {
            max_extra_dots: dots.unwrap_or(2),
            max_numerator_rank: rank.unwrap_or(1),
        }),
    }
}

pub(super) fn json_report(py: Python<'_>, document: &str) -> PyResult<Py<PyDict>> {
    Ok(PyModule::import(py, "json")?
        .call_method1("loads", (document,))?
        .cast_into::<PyDict>()?
        .unbind())
}

#[cfg(all(unix, not(target_arch = "wasm32")))]
#[allow(clippy::too_many_arguments)]
pub(super) fn transfer_options(
    sha256: &str,
    archive_bytes: Option<u64>,
    extracted_bytes: Option<u64>,
    entries: Option<usize>,
    temporary_bytes: Option<u64>,
    disk_bytes: Option<u64>,
    timeout_seconds: Option<f64>,
) -> PyResult<crate::RustRedCacheTransferOptions> {
    if sha256.len() != 64 || !sha256.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "sha256 must be exactly 64 hexadecimal digits from a trusted source",
        ));
    }
    let mut digest = [0u8; 32];
    for (index, value) in digest.iter_mut().enumerate() {
        *value = u8::from_str_radix(&sha256[index * 2..index * 2 + 2], 16)
            .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
    }
    let mut options = crate::RustRedCacheTransferOptions::new(digest);
    if let Some(value) = archive_bytes {
        options.max_archive_bytes = value;
    }
    if let Some(value) = extracted_bytes {
        options.max_extracted_bytes = value;
    }
    if let Some(value) = entries {
        options.max_entries = value;
    }
    if let Some(value) = temporary_bytes {
        options.max_temporary_bytes = value;
    }
    if let Some(value) = disk_bytes {
        options.storage.max_disk_bytes = value;
    }
    if let Some(value) = timeout_seconds {
        options.timeout = Duration::try_from_secs_f64(value).map_err(|_| {
            pyo3::exceptions::PyValueError::new_err("timeout_seconds must be finite and positive")
        })?;
    }
    Ok(options)
}

#[cfg(all(unix, not(target_arch = "wasm32")))]
pub(super) fn transfer_report(
    py: Python<'_>,
    report: crate::RustRedCacheTransferReport,
) -> PyResult<Py<PyDict>> {
    let result = PyDict::new(py);
    result.set_item("output", report.output.to_string_lossy())?;
    result.set_item(
        "sha256",
        report
            .sha256
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect::<String>(),
    )?;
    result.set_item("archive_bytes", report.archive_bytes)?;
    result.set_item("extracted_bytes", report.extracted_bytes)?;
    result.set_item("entries", report.entries)?;
    result.set_item("namespaces", report.copied.namespaces)?;
    result.set_item("copied_bases", report.copied.copied_bases)?;
    result.set_item("copied_chunks", report.copied.copied_chunks)?;
    result.set_item(
        "preserved_record_references",
        report.copied.preserved_record_references,
    )?;
    result.set_item("compressed_bytes", report.copied.compressed_bytes)?;
    result.set_item("index_bytes", report.copied.index_bytes)?;
    result.set_item("transfer_buffer_bytes", report.copied.transfer_buffer_bytes)?;
    Ok(result.unbind())
}

struct PythonProgress {
    progress: Option<Py<PyAny>>,
    cancel: Option<Py<PyAny>>,
    interval: Duration,
    last: Option<Instant>,
    exception: Option<PyErr>,
}
impl PythonProgress {
    fn observe(
        &mut self,
        report: &CachePreparationReport,
        cancelled: &AtomicBool,
    ) -> Result<(), AppError> {
        // Once a Python error is captured, allow the native final observer to
        // finish checkpointing without invoking Python again or replacing it.
        if self.exception.is_some() {
            return Ok(());
        }
        let result = Python::attach(|py| -> PyResult<bool> {
            py.check_signals()?;
            if report.status == CachePreparationStatus::Running
                && let Some(callback) = self.cancel.as_ref()
                && callback.bind(py).call0()?.extract::<bool>()?
            {
                return Ok(true);
            }
            let due = report.status != CachePreparationStatus::Running
                || self.last.is_none_or(|last| last.elapsed() >= self.interval);
            if due && let Some(callback) = self.progress.as_ref() {
                let value = callback
                    .bind(py)
                    .call1((json_report(py, &report.to_document().to_string())?,))?;
                self.last = Some(Instant::now());
                if !value.is_none()
                    && !value.extract::<bool>()?
                    && report.status == CachePreparationStatus::Running
                {
                    return Ok(true);
                }
            }
            Ok(false)
        });
        match result {
            Ok(false) => Ok(()),
            Ok(true) => {
                cancelled.store(true, Ordering::Relaxed);
                Err(AppError::new(
                    AppErrorKind::Cancelled,
                    "Python callback requested cooperative cancellation",
                ))
            }
            Err(error) => {
                self.exception = Some(error);
                cancelled.store(true, Ordering::Relaxed);
                Err(AppError::new(
                    AppErrorKind::Cancelled,
                    "Python signal or callback interrupted preparation",
                ))
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare(
    py: Python<'_>,
    engine: Vakint,
    loops: usize,
    output: PathBuf,
    options: RustRedCachePreparationOptions,
    progress: Option<Py<PyAny>>,
    cancel: Option<Py<PyAny>>,
    interval: f64,
) -> PyResult<Py<PyDict>> {
    for callback in [&progress, &cancel].into_iter().flatten() {
        if !callback.bind(py).is_callable() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "progress and cancel must be callable",
            ));
        }
    }
    let interval = Duration::try_from_secs_f64(interval).map_err(|_| {
        pyo3::exceptions::PyValueError::new_err(
            "progress_interval_seconds must be finite and nonnegative",
        )
    })?;
    let cancelled = AtomicBool::new(false);
    let mut observer = PythonProgress {
        progress,
        cancel,
        interval,
        last: None,
        exception: None,
    };
    let result = {
        let observer = &mut observer;
        let cancelled = &cancelled;
        py.detach(move || {
            engine.prepare_cache(loops, &output, &options, cancelled, |report| {
                observer.observe(report, cancelled)
            })
        })
    };
    if let Some(error) = observer.exception.take() {
        if let Ok(report) = &result
            && let Ok(document) = json_report(py, &report.to_document().to_string())
        {
            // Preserve the original exception even if a custom exception
            // class refuses an added attribute. The report records any
            // checkpoint failure; only prior successful writes are durable.
            let _ = error.value(py).setattr("cache_report", document);
        }
        return Err(error);
    }
    let report = result.map_err(vakint_to_python_error)?;
    json_report(py, &report.to_document().to_string())
}

fn native_stats(py: Python<'_>, value: GuardedRootCacheStatistics) -> PyResult<Bound<'_, PyDict>> {
    let result = PyDict::new(py);
    for (key, value) in [
        ("ram_hits", value.ram_hits),
        ("backing_hits", value.backing_hits),
        ("descendant_ram_hits", value.descendant_ram_hits),
        ("descendant_backing_hits", value.descendant_backing_hits),
        ("completed_reductions", value.completed_reductions),
        ("payload_validations", value.payload_validations),
        ("evicted_entries", value.evicted_entries),
        ("evicted_coefficient_terms", value.evicted_coefficient_terms),
        ("evicted_coefficient_bytes", value.evicted_coefficient_bytes),
        ("cached_integrals", value.cached_integrals),
        ("cached_coefficient_terms", value.cached_coefficient_terms),
        ("cached_coefficient_bytes", value.cached_coefficient_bytes),
        ("rule_applications", value.work.rule_applications()),
        ("cache_hits", value.work.cache_hits()),
        ("coalescing_additions", value.work.coalescing_additions()),
    ] {
        result.set_item(key, value)?;
    }
    Ok(result)
}
fn disk_stats(
    py: Python<'_>,
    value: SavedApplicationCachePackStatistics,
) -> PyResult<Bound<'_, PyDict>> {
    let result = PyDict::new(py);
    for (key, value) in [
        ("runs", value.runs),
        ("chunks", value.chunks),
        ("records", value.records),
        ("manifest_bytes", value.manifest_bytes),
        ("directory_owned_bytes", value.directory_owned_bytes),
        ("chunks_read", value.chunks_read),
        ("largest_chunk_read", value.largest_chunk_read),
        ("largest_decoded_chunk", value.largest_decoded_chunk),
        ("runs_probed", value.runs_probed),
        ("key_index_skips", value.key_index_skips),
        ("lookup_limit_refusals", value.lookup_limit_refusals),
        ("hits", value.hits),
        ("misses", value.misses),
    ] {
        result.set_item(key, value)?;
    }
    for (key, value) in [
        ("referenced_chunk_bytes", value.referenced_chunk_bytes),
        ("chunk_bytes_read", value.chunk_bytes_read),
        ("decoded_chunk_bytes", value.decoded_chunk_bytes),
        (
            "lookup_encoded_bytes_charged",
            value.lookup_encoded_bytes_charged,
        ),
        (
            "lookup_decoded_bytes_charged",
            value.lookup_decoded_bytes_charged,
        ),
    ] {
        result.set_item(key, value)?;
    }
    Ok(result)
}

pub(super) fn statistics(py: Python<'_>, stats: &RustRedCacheStatistics) -> PyResult<Py<PyDict>> {
    let result = PyDict::new(py);
    result.set_item("attached", stats.attached)?;
    result.set_item("native", native_stats(py, stats.native)?)?;
    result.set_item("backing_budget_fallbacks", stats.backing_budget_fallbacks)?;
    result.set_item("ram_is_process_rss", false)?;
    result.set_item(
        "saved_application",
        stats
            .saved_application
            .as_ref()
            .map(|value| json_report(py, &value.to_document().to_string()))
            .transpose()?,
    )?;
    result.set_item(
        "disk",
        stats.disk.map(|disk| disk_stats(py, disk)).transpose()?,
    )?;
    let namespaces = PyList::empty(py);
    for namespace in &stats.namespaces {
        let item = PyDict::new(py);
        let binding = PyDict::new(py);
        binding.set_item(
            "artifact_digest",
            namespace
                .binding
                .artifact_digest
                .iter()
                .map(|byte| format!("{byte:02x}"))
                .collect::<String>(),
        )?;
        binding.set_item("family_fingerprint", &namespace.binding.family_fingerprint)?;
        binding.set_item("physical_arity", namespace.binding.physical_arity)?;
        binding.set_item("layer", &namespace.binding.layer)?;
        item.set_item("binding", binding)?;
        item.set_item("native", native_stats(py, namespace.native)?)?;
        item.set_item(
            "disk",
            namespace
                .disk
                .map(|disk| disk_stats(py, disk))
                .transpose()?,
        )?;
        item.set_item(
            "backing_budget_fallbacks",
            namespace.backing_budget_fallbacks,
        )?;
        namespaces.append(item)?;
    }
    result.set_item("namespaces", namespaces)?;
    Ok(result.unbind())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn controls_preserve_omitted_limits_and_validate_workload_modes() {
        let mut limits = ReductionLimits {
            max_cached_coefficient_bytes: 42,
            max_rule_applications: 17,
            ..Default::default()
        };
        apply_native_overrides(&mut limits, None, Some(99), None, None).unwrap();
        assert_eq!(limits.max_cached_coefficient_bytes, 42);
        assert_eq!(limits.max_rule_applications, 17);
        assert_eq!(limits.max_cached_coefficient_terms, 99);
        assert!(matches!(
            workload(None, None, None).unwrap(),
            RustRedCacheWorkload::Primer
        ));
        assert!(workload(Some(1), None, Some(vec![(0, vec![1])])).is_err());
        assert!(save_mode("automatic").is_err());
        assert!(lookup_policy("ignore-corruption").is_err());
        for loops in [1, 2, 3, 4, 5] {
            check_loops(loops).unwrap();
        }
        for loops in [0, 6] {
            assert!(check_loops(loops).is_err());
        }
    }
}
