//! Native Python resource configuration shared by embedding hosts.

use pyo3::types::{PyCFunction, PyModule, PyModuleMethods};
use pyo3::{Bound, PyResult, pyfunction, wrap_pyfunction};
use serde_json::{Value, json};

use crate::{Asset, Data};

#[pyfunction]
fn _set_data_directories(directories: Vec<String>) {
    Data::set_directories(directories.into_iter().map(Into::into).collect());
}

/// Register resource discovery using the host's supported assets.
///
/// ```no_run
/// use pyo3::types::PyModule;
/// use pyo3::{PyErr, Python};
/// use symbolica_data::python;
/// Python::attach(|py| {
///     let core = PyModule::new(py, "core")?;
///     python::register(&core, &[])
/// })?;
/// # Ok::<(), PyErr>(())
/// ```
///
/// # Errors
/// Returns Python errors if functions cannot be registered on the module.
pub fn register(module: &Bound<'_, PyModule>, assets: &'static [Asset]) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(_set_data_directories, module)?)?;
    let manifest = PyCFunction::new_closure(
        module.py(),
        Some(c"_data_manifest"),
        None,
        move |_, _| -> PyResult<String> {
            let records: Vec<Value> = assets
                .iter()
                .map(|asset| {
                    json!({
                        "name": asset.name,
                        "sha256": asset.sha256,
                        "size": asset.size,
                        "description": asset.description,
                        "loaded": asset.is_loaded(),
                    })
                })
                .collect();
            Ok(json!({"schema": 1, "bundle": "data", "assets": records}).to_string())
        },
    )?;
    module.add("_data_manifest", manifest)
}
