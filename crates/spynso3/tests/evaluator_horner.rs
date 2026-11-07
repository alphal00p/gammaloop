use std::ffi::CString;

use pyo3::types::{PyDict, PyDictMethods, PyModule};
use pyo3::{PyResult, Python};
use spynso3::SpensoModule;
use symbolica::api::python;
use symbolica::api::python::SymbolicaCommunityModule;

#[test]
fn tensor_evaluator_forwards_both_horner_scoring_modes() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "horner_core")?;
        python::create_symbolica_module(&core)?;
        let tensor = PyModule::new(py, "horner_tensor")?;
        SpensoModule::register_module(&tensor)?;
        SpensoModule::initialize(py)?;
        let locals = PyDict::new(py);
        locals.set_item("core", core)?;
        locals.set_item("tensor", tensor)?;
        let source = CString::new(include_str!("evaluator_horner.py")).unwrap();
        py.run(&source, Some(&locals), Some(&locals))
    })
    .unwrap();
}
