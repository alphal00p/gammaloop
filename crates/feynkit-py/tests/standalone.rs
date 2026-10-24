use std::ffi::CString;

use feynkit_py::FeynkitModule;
#[cfg(feature = "linnet-python")]
use linnet_py::GraphModule;
#[cfg(not(feature = "linnet-python"))]
use linnet_render_py::GraphModule;
use pyo3::types::{PyDict, PyDictMethods, PyModule};
use pyo3::{PyResult, Python};
use spynso3::SpensoModule;
use symbolica::api::python;
use symbolica::api::python::SymbolicaCommunityModule;

// Keep this host in its own process: installed-host tests register canonical
// Python packages, which would hide a dependency on those import names.
#[test]
fn standalone_host_requires_no_canonical_python_namespaces() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let core = PyModule::new(py, "standalone_core")?;
        python::create_symbolica_module(&core)?;
        let graph = PyModule::new(py, "standalone_graph")?;
        GraphModule::register_module(&graph)?;
        let tensor = PyModule::new(py, "standalone_tensor")?;
        SpensoModule::register_module(&tensor)?;
        let hep = PyModule::new(py, "standalone_hep")?;
        FeynkitModule::register_module(&hep)?;
        GraphModule::initialize(py)?;
        SpensoModule::initialize(py)?;
        FeynkitModule::initialize(py)?;
        let locals = PyDict::new(py);
        locals.set_item("core", core)?;
        locals.set_item("graph", graph)?;
        locals.set_item("tensor", tensor)?;
        locals.set_item("hep", hep)?;
        locals.set_item("linnet_python", cfg!(feature = "linnet-python"))?;
        let source = CString::new(include_str!("standalone_host.py")).unwrap();
        py.run(&source, Some(&locals), Some(&locals))
    })
    .unwrap();
}
