use pyo3::prelude::*;
use pyo3::{types::PyModule, Py, PyAny, PyResult, Python};

#[test]

fn python() -> PyResult<()> {
    pyo3::prepare_freethreaded_python();
    let py_api = include_str!("../../python/gammaloop/interface/gammaloop_interface.py");

    let a = Python::with_gil(|py| -> PyResult<()> {
        let interface: Py<PyAny> = PyModule::from_code(py, py_api, "", "")?
            .getattr("Gammaloop")?
            .into();
        let gloop = interface.call0(py)?;
        gloop.call_method(py, "do_import_model", ("help",), None)?;

        Ok(())
    })?;

    Ok(())
}
