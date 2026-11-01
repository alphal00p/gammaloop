//! Python surface of the persistent Ordered step cache.
//!
//! The depth-first cache operations (preparation, loading, running caches,
//! transport and statistics) were removed together with depth-first reduction.
//! Engine keywords select the Ordered cache settings; the reports below are
//! plain dictionaries.
use pyo3::types::{PyDict, PyDictMethods, PyList, PyListMethods};
use pyo3::{Bound, PyResult, Python};

use crate::{
    OrderedCacheAutosave, OrderedCacheMode, OrderedCacheSaveReport, OrderedCacheSettings,
    OrderedCacheStatistics, OrderedCacheWarmReport,
};

use super::vakint_to_python_error;

/// Validate the `ordered_cache*` constructor keywords.
pub(super) fn settings(
    mode: &str,
    directory: Option<std::path::PathBuf>,
    ram_bytes: Option<usize>,
    autosave: &str,
) -> PyResult<OrderedCacheSettings> {
    let defaults = OrderedCacheSettings::default();
    Ok(OrderedCacheSettings {
        mode: mode
            .parse::<OrderedCacheMode>()
            .map_err(vakint_to_python_error)?,
        directory,
        max_resident_bytes: ram_bytes.or(defaults.max_resident_bytes),
        autosave: autosave
            .parse::<OrderedCacheAutosave>()
            .map_err(vakint_to_python_error)?,
        // Shipped caches seed every mode (Rust default).
        shipped: defaults.shipped,
    })
}

pub(super) fn statistics<'py>(
    py: Python<'py>,
    statistics: &OrderedCacheStatistics,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("records", statistics.records)?;
    result.set_item("hits", statistics.hits)?;
    result.set_item("misses", statistics.misses)?;
    result.set_item("inserts", statistics.inserts)?;
    result.set_item("resident_bytes", statistics.resident_bytes)?;
    result.set_item("file_bytes", statistics.file_bytes)?;
    result.set_item(
        "last_autosave_error",
        statistics.last_autosave_error.clone(),
    )?;
    let programs = PyList::empty(py);
    for program in &statistics.programs {
        let native = &program.native;
        let record = PyDict::new(py);
        record.set_item("loop_count", program.loop_count)?;
        record.set_item("program", &program.program)?;
        record.set_item("data_version", &program.data_version)?;
        record.set_item("rule_tables", &program.rule_tables)?;
        record.set_item("stamp", &program.stamp)?;
        record.set_item("directory", program.directory.clone())?;
        record.set_item("mode", program.mode.as_str())?;
        record.set_item("autosave", program.autosave.as_str())?;
        record.set_item(
            "shipped",
            program
                .shipped
                .as_ref()
                .map(|shipped| shipped.status.as_str()),
        )?;
        record.set_item("records", native.records)?;
        record.set_item("hits", native.hits)?;
        record.set_item("misses", native.misses)?;
        record.set_item("inserts", native.inserts)?;
        record.set_item("skipped_inserts", native.skipped_inserts)?;
        record.set_item("resident_bytes", native.resident_bytes)?;
        record.set_item("file_bytes", native.file_bytes)?;
        record.set_item("load_seconds", native.load_seconds)?;
        record.set_item("shards", native.shards.len())?;
        record.set_item("loaded_shards", native.loaded_shards)?;
        record.set_item("skipped_shards", native.skipped_shards)?;
        record.set_item("rejected_shards", native.rejected_shards)?;
        record.set_item("factor_basis", native.factor_basis)?;
        record.set_item("budget_exhausted", native.budget_exhausted)?;
        programs.append(record)?;
    }
    result.set_item("programs", programs)?;
    Ok(result)
}

pub(super) fn save_report<'py>(
    py: Python<'py>,
    report: &OrderedCacheSaveReport,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("shards_written", report.shards_written)?;
    result.set_item("records_written", report.records_written)?;
    result.set_item("bytes_written", report.bytes_written)?;
    let programs = PyList::empty(py);
    for program in &report.programs {
        let record = PyDict::new(py);
        record.set_item("loop_count", program.loop_count)?;
        record.set_item("program", &program.program)?;
        record.set_item("directory", program.native.directory.clone())?;
        record.set_item("shards_written", program.native.shards_written)?;
        record.set_item("records_written", program.native.records_written)?;
        record.set_item("bytes_written", program.native.bytes_written)?;
        record.set_item("shards", program.native.shards)?;
        record.set_item("records", program.native.records)?;
        record.set_item("bytes", program.native.bytes)?;
        programs.append(record)?;
    }
    result.set_item("programs", programs)?;
    Ok(result)
}

pub(super) fn warm_report<'py>(
    py: Python<'py>,
    report: &OrderedCacheWarmReport,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("integrals", report.integrals)?;
    result.set_item("seconds", report.seconds)?;
    result.set_item("hits", report.hits)?;
    result.set_item("misses", report.misses)?;
    result.set_item("inserts", report.inserts)?;
    result.set_item(
        "save",
        report
            .save
            .as_ref()
            .map(|save| save_report(py, save))
            .transpose()?,
    )?;
    result.set_item("statistics", statistics(py, &report.statistics)?)?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ordered_cache_keywords_are_validated_and_keep_the_default_ceiling() {
        let settings = settings("read-only", Some("cache".into()), None, "never").unwrap();
        assert_eq!(settings.mode, OrderedCacheMode::ReadOnly);
        assert_eq!(settings.autosave, OrderedCacheAutosave::Never);
        assert_eq!(
            settings.max_resident_bytes,
            OrderedCacheSettings::default().max_resident_bytes
        );
        assert_eq!(
            super::settings("read-write", None, Some(1 << 20), "on-drop")
                .unwrap()
                .max_resident_bytes,
            Some(1 << 20)
        );
        assert!(super::settings("sometimes", None, None, "never").is_err());
        assert!(super::settings("off", None, None, "always").is_err());
    }
}
