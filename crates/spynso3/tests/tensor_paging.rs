//! Run the public Python paging regressions against this build's native owners.
//! Requires the notebook Python dependencies and configured Symbolica fonts.
//! Run explicitly with `cargo test -p spynso3 --test tensor_paging -- --ignored`.
#![cfg(feature = "native")]

use spynso3::SpensoModule;
use std::ffi::CString;
use std::{env, fs};

use pyo3::types::PyModule;
use pyo3::types::{PyAnyMethods, PyModuleMethods};
use pyo3::{PyResult, Python};
use symbolica::api::python;
use symbolica::api::python::SymbolicaCommunityModule;

#[test]
#[ignore = "requires optional notebook Python dependencies and external Symbolica fonts"]
fn direct_tensor_paging() {
    Python::initialize();
    Python::attach(|py| {
        let result = (|| -> PyResult<()> {
            let modules = py.import("sys")?.getattr("modules")?;
            let symbolica = PyModule::new(py, "symbolica")?;
            symbolica.add("__path__", Vec::<String>::new())?;
            python::create_symbolica_module(&symbolica)?;
            modules.set_item("symbolica", &symbolica)?;
            modules.set_item("symbolica.core", &symbolica)?;
            let community = PyModule::new(py, "symbolica.community")?;
            community.add("__path__", Vec::<String>::new())?;
            symbolica.add("community", &community)?;
            modules.set_item("symbolica.community", &community)?;

            let tensor = PyModule::new(py, "symbolica.community.tensor")?;
            SpensoModule::register_module(&tensor)?;
            SpensoModule::initialize(py)?;
            community.add("tensor", &tensor)?;
            community.add("tensor_native", &tensor)?;
            modules.set_item("symbolica.community.tensor", &tensor)?;
            modules.set_item("symbolica.community.tensor_native", &tensor)?;

            let code = CString::new(include_str!("installed_tensor_paging.py")).unwrap();
            let filename = CString::new(concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/installed_tensor_paging.py"
            ))
            .unwrap();
            let tests = PyModule::from_code(py, &code, &filename, c"installed_tensor_paging")?;
            let unittest = py.import("unittest")?;
            let suite = unittest
                .getattr("defaultTestLoader")?
                .call_method1("loadTestsFromModule", (&tests,))?;
            assert!(suite.call_method0("countTestCases")?.extract::<usize>()? > 0);
            let result = unittest
                .getattr("TextTestRunner")?
                .call0()?
                .call_method1("run", (suite,))?;
            assert!(result.call_method0("wasSuccessful")?.is_truthy()?);
            if let Some(path) = env::var_os("SPENSO_MATHML_BROWSER_FIXTURES") {
                let fixtures = tests.call_method0("browser_fixtures")?;
                let json = py.import("json")?.call_method1("dumps", (fixtures,))?;
                fs::write(path, json.extract::<String>()?)?;
            }
            Ok(())
        })();
        if let Err(error) = result {
            error.print(py);
            panic!("Python tensor paging regressions failed");
        }
    });
}
