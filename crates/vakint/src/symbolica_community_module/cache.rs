//! Thin Python steering of the existing native cache APIs, with cooperative
//! Python callbacks only at the runner's ordinary root/checkpoint boundaries.
use std::{
    path::PathBuf,
    sync::atomic::{AtomicBool, Ordering},
    time::{Duration, Instant},
};

use pyo3::{
    prelude::*,
    types::{PyAnyMethods, PyBytes, PyDict, PyDictMethods, PyList, PyListMethods, PyModule},
};
use rustred::{
    family::IntegralKey,
    reduction::{GuardedRootCacheStatistics, ReductionLimits},
};
use rustred_app::{
    AppError, AppErrorKind, CacheLookupBudgetPolicy, CachePreparationReport,
    CachePreparationStatus, SavedApplicationCachePackStatistics,
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::*;

use super::vakint_to_python_error;
use crate::{
    RustRedCachePreparationOptions, RustRedCacheRoot, RustRedCacheSaveMode, RustRedCacheStatistics,
    RustRedCacheWorkload, RustRedMonitorDisplay, RustRedMonitorOptions, Vakint,
};

/// Read-only native single-file content identity. Obtain this document from an
/// authenticated manifest; constructing the record does not authenticate its source.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CacheFileIdentity",
    module = "symbolica.community.hepkit.vakint",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct CacheFileIdentityWrapper {
    pub(super) native: rustred_app::SavedApplicationCacheFileIdentity,
}

impl super::ModuleInit for CacheFileIdentityWrapper {}

impl CacheFileIdentityWrapper {
    fn parse(document: serde_json::Value) -> PyResult<Self> {
        serde_json::from_value(document)
            .map(|native| Self { native })
            .map_err(|error| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "invalid cache file identity document: {error}"
                ))
            })
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl CacheFileIdentityWrapper {
    #[new]
    #[pyo3(signature = (*, compressed_bytes, compressed_blake3, decoded_bytes, decoded_blake3))]
    fn new(
        compressed_bytes: u64,
        compressed_blake3: Vec<u8>,
        decoded_bytes: u64,
        decoded_blake3: Vec<u8>,
    ) -> PyResult<Self> {
        Self::parse(serde_json::json!({
            "compressed_bytes": compressed_bytes,
            "compressed_blake3": compressed_blake3,
            "decoded_bytes": decoded_bytes,
            "decoded_blake3": decoded_blake3,
        }))
    }

    /// Parse the exact native identity schema, rejecting missing or extra fields.
    #[staticmethod]
    fn from_document(py: Python<'_>, document: &Bound<'_, PyDict>) -> PyResult<Self> {
        Self::parse(json_document(py, document, "cache file identity")?)
    }

    /// Return a detached document in the native serializer's schema.
    fn to_document(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
        let document = serde_json::to_string(&self.native).map_err(|error| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "cannot serialize cache file identity: {error}"
            ))
        })?;
        json_report(py, &document)
    }

    #[getter]
    fn compressed_bytes(&self) -> u64 {
        self.native.compressed_bytes
    }

    #[getter]
    fn compressed_blake3<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.native.compressed_blake3)
    }

    #[getter]
    fn decoded_bytes(&self) -> u64 {
        self.native.decoded_bytes
    }

    #[getter]
    fn decoded_blake3<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, &self.native.decoded_blake3)
    }
}

/// Explicit native single-file transport and decoded-store grants. Ordinary
/// construction retains the compact distribution defaults.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "CacheFileOptions",
    module = "symbolica.community.hepkit.vakint",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct CacheFileOptionsWrapper {
    pub(super) native: crate::RustRedCacheFileOptions,
}

impl super::ModuleInit for CacheFileOptionsWrapper {}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl CacheFileOptionsWrapper {
    #[new]
    #[pyo3(signature = (*, max_compressed_bytes=None, max_decoded_bytes=None, max_index_bytes=None, decoded_directory=None))]
    fn new(
        max_compressed_bytes: Option<u64>,
        max_decoded_bytes: Option<u64>,
        max_index_bytes: Option<u64>,
        decoded_directory: Option<PathBuf>,
    ) -> PyResult<Self> {
        let mut native = crate::RustRedCacheFileOptions {
            decoded_directory,
            ..Default::default()
        };
        let addressable = |value| {
            usize::try_from(value).map_err(|_| {
                pyo3::exceptions::PyValueError::new_err(
                    "cache byte grant exceeds this platform's address space",
                )
            })
        };
        if let Some(value) = max_compressed_bytes {
            native.limits.max_compressed_bytes = addressable(value)?;
        }
        if let Some(value) = max_decoded_bytes {
            native.limits = native
                .limits
                .with_decoded_bytes(addressable(value)?)
                .map_err(native_error)?;
        }
        if let Some(value) = max_index_bytes {
            native.limits.chunk.max_index_owned_bytes = addressable(value)?;
        }
        Ok(Self { native })
    }

    #[getter]
    fn max_compressed_bytes(&self) -> usize {
        self.native.limits.max_compressed_bytes
    }

    #[getter]
    fn max_decoded_bytes(&self) -> usize {
        self.native.limits.max_decoded_bytes()
    }

    #[getter]
    fn max_index_bytes(&self) -> usize {
        self.native.limits.chunk.max_index_owned_bytes
    }

    #[getter]
    fn decoded_directory(&self) -> Option<PathBuf> {
        self.native.decoded_directory.clone()
    }
}

/// Independent native I/O control; it never borrows a prepared reduction owner.
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "RunningCacheControl",
    module = "symbolica.community.hepkit.vakint",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct RunningCacheControlWrapper {
    pub(super) control: crate::RustRedRunningCacheControl,
    pub(super) resources: crate::ResourceLimits,
}

impl super::ModuleInit for RunningCacheControlWrapper {}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl RunningCacheControlWrapper {
    /// Read native queue and durability counters without locking the reducer.
    pub fn report(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
        super::resources::run(self.resources, || {
            let control = self.control.clone();
            let report = py.detach(move || control.report()).map_err(native_error)?;
            running_report(py, &report)
        })
    }

    /// Acknowledge a native durability barrier, independently of evaluation.
    pub fn flush(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
        super::resources::run_committing(self.resources, || {
            let control = self.control.clone();
            let report = py.detach(move || control.flush()).map_err(native_error)?;
            running_report(py, &report)
        })
    }

    /// Freeze the current delta generation into a native descriptor.
    pub fn snapshot(&self, py: Python<'_>, sidecar_path: PathBuf) -> PyResult<Py<PyDict>> {
        super::resources::run_committing(self.resources, || {
            let control = self.control.clone();
            let snapshot = py
                .detach(move || control.snapshot(&sidecar_path))
                .map_err(native_error)?;
            json_report(py, &snapshot.to_document().to_string())
        })
    }

    /// Flush and release the native writer lock for every cloned control.
    pub fn close(&self, py: Python<'_>) -> PyResult<Py<PyDict>> {
        super::resources::run_committing(self.resources, || {
            let control = self.control.clone();
            let report = py.detach(move || control.close()).map_err(native_error)?;
            running_report(py, &report)
        })
    }
}

fn native_error(source: AppError) -> PyErr {
    vakint_to_python_error(crate::RustRedEvaluationError::CacheOperation { source }.into())
}

pub(super) fn monitor_options(
    monitor: &str,
    progress_snapshot: Option<PathBuf>,
    interval: f64,
) -> PyResult<RustRedMonitorOptions> {
    let display = match monitor {
        "off" => RustRedMonitorDisplay::Off,
        "tui" => RustRedMonitorDisplay::Tui,
        "json" => RustRedMonitorDisplay::Json,
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "monitor must be 'off', 'tui' or 'json'",
            ));
        }
    };
    let interval = Duration::try_from_secs_f64(interval).map_err(|_| {
        pyo3::exceptions::PyValueError::new_err(
            "progress_interval_seconds must be finite and nonnegative",
        )
    })?;
    let options = RustRedMonitorOptions {
        display,
        progress_snapshot,
        interval,
        ..Default::default()
    };
    options
        .validate()
        .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
    Ok(options)
}

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

#[allow(clippy::too_many_arguments)]
pub(super) fn running_options(
    py: Python<'_>,
    read_backing: bool,
    read_delta: bool,
    base_paths: Option<Vec<PathBuf>>,
    snapshots: Option<Vec<Py<PyDict>>>,
    disk_bytes: Option<u64>,
    max_pending_records: Option<usize>,
    max_pending_bytes: Option<usize>,
    batch_records: Option<usize>,
    batch_bytes: Option<usize>,
    checkpoint_interval_seconds: Option<f64>,
) -> PyResult<crate::RustRedRunningCacheOptions> {
    let mut options = crate::RustRedRunningCacheOptions {
        read_backing,
        read_delta,
        base_paths: base_paths.unwrap_or_default(),
        ..Default::default()
    };
    if let Some(snapshots) = snapshots {
        options.snapshots = snapshots
            .into_iter()
            .map(|snapshot| {
                let document = json_document(py, snapshot.bind(py), "running-cache snapshot")?;
                crate::RustRedRunningCacheSnapshot::from_document(&document).map_err(|source| {
                    vakint_to_python_error(
                        crate::RustRedEvaluationError::CacheOperation { source }.into(),
                    )
                })
            })
            .collect::<PyResult<Vec<_>>>()?;
    }
    if let Some(value) = disk_bytes {
        options.native.storage.max_disk_bytes = value;
    }
    if let Some(value) = max_pending_records {
        options.native.max_pending_records = value;
    }
    if let Some(value) = max_pending_bytes {
        options.native.max_pending_bytes = value;
    }
    if let Some(value) = batch_records {
        options.native.batch_records = value;
    }
    if let Some(value) = batch_bytes {
        options.native.batch_bytes = value;
    }
    if let Some(value) = checkpoint_interval_seconds {
        options.native.checkpoint_interval = Duration::try_from_secs_f64(value).map_err(|_| {
            pyo3::exceptions::PyValueError::new_err(
                "checkpoint_interval_seconds must be finite and nonnegative",
            )
        })?;
    }
    Ok(options)
}

pub(super) fn running_report(
    py: Python<'_>,
    report: &crate::RustRedRunningCacheReport,
) -> PyResult<Py<PyDict>> {
    json_report(py, &report.to_document().to_string())
}

pub(super) fn apply_native_overrides(
    limits: &mut ReductionLimits,
    ram_bytes: Option<u64>,
    terms: Option<usize>,
    integrals: Option<usize>,
    rules: Option<usize>,
) -> PyResult<()> {
    if terms.is_some() || integrals.is_some() || rules.is_some() {
        *limits = limits.with_diagnostics();
    }
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

fn json_document(
    py: Python<'_>,
    document: &Bound<'_, PyDict>,
    description: &str,
) -> PyResult<serde_json::Value> {
    let document: String = PyModule::import(py, "json")?
        .call_method1("dumps", (document,))?
        .extract()?;
    serde_json::from_str(&document).map_err(|error| {
        pyo3::exceptions::PyValueError::new_err(format!("invalid {description} document: {error}"))
    })
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
    format: &str,
    manifest_bytes: Option<u64>,
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
    let mut options = match format {
        "archive" => crate::RustRedCacheTransferOptions::new(digest),
        "multipart_manifest" => crate::RustRedCacheTransferOptions::multipart_manifest(digest),
        _ => {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "format must be 'archive' or 'multipart_manifest'",
            ));
        }
    };
    if let Some(value) = manifest_bytes {
        options.max_manifest_bytes = value;
    }
    if let Some(value) = archive_bytes {
        options.max_archive_bytes = value;
    }
    if let Some(value) = extracted_bytes {
        options.max_extracted_bytes = value;
    }
    if let Some(value) = entries {
        options = options.with_diagnostics();
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
    format: &str,
) -> PyResult<Py<PyDict>> {
    let result = PyDict::new(py);
    result.set_item("format", format)?;
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
    monitor: RustRedMonitorOptions,
) -> PyResult<Py<PyDict>> {
    for callback in [&progress, &cancel].into_iter().flatten() {
        if !callback.bind(py).is_callable() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "progress and cancel must be callable",
            ));
        }
    }
    let interval = monitor.interval;
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
            if monitor.enabled() {
                engine.prepare_cache_with_monitor(
                    loops,
                    &output,
                    &options,
                    cancelled,
                    |report| observer.observe(report, cancelled),
                    &monitor,
                )
            } else {
                engine.prepare_cache(loops, &output, &options, cancelled, |report| {
                    observer.observe(report, cancelled)
                })
            }
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
        assert!(limits.diagnostics_enabled);
        assert_eq!(limits.max_cached_coefficient_bytes, 42);
        assert_eq!(limits.max_rule_applications, 17);
        assert_eq!(limits.max_cached_coefficient_terms, 99);
        let mut ordinary = ReductionLimits::default();
        apply_native_overrides(&mut ordinary, Some(1024), None, None, None).unwrap();
        assert!(!ordinary.diagnostics_enabled);
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
