use std::ffi::CString;

use feynkit_graph::FeynmanDiagram;
use feynkit_py::{FeynkitModule, PyFeynmanDiagram};
use pyo3::{
    prelude::*,
    types::{PyCFunction, PyDict, PyList, PyModule},
};
use spynso3::SpensoModule;
use symbolica::api::python::{SymbolicaCommunityModule, create_symbolica_module};

const FEYNKIT_WRAPPER: &str = include_str!("../python/symbolica/community/feynkit/__init__.py");
const SPENSO_WRAPPER: &str = "from ..spenso_native import *\n\ninitialize_module()\n";

fn install_package<'py>(py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyModule>> {
    let package = PyModule::new(py, name)?;
    package.setattr("__package__", name)?;
    package.setattr("__path__", PyList::empty(py))?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(name, &package)?;
    Ok(package)
}

fn register_native<'py, M>(core: &Bound<'py, PyModule>) -> PyResult<Bound<'py, PyModule>>
where
    M: SymbolicaCommunityModule + 'static,
{
    let py = core.py();
    let name = M::get_name();
    let native_name = format!("{name}_native");
    let native = PyModule::new(py, &native_name)?;
    let initialize_module =
        PyCFunction::new_closure(py, Some(c"initialize_module"), None, |args, _kwargs| {
            M::initialize(args.py())
        })?;
    native.add("initialize_module", initialize_module)?;
    M::register_module(&native)?;
    core.add_submodule(&native)?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(format!("symbolica.community.{native_name}"), &native)?;
    Ok(native)
}

fn import_wrapper<'py>(
    py: Python<'py>,
    community: &Bound<'py, PyModule>,
    name: &str,
    source: &str,
) -> PyResult<Bound<'py, PyModule>> {
    let full_name = format!("symbolica.community.{name}");
    let wrapper = PyModule::new(py, &full_name)?;
    wrapper.setattr("__package__", &full_name)?;
    wrapper.setattr("__path__", PyList::empty(py))?;
    py.import("sys")?
        .getattr("modules")?
        .set_item(&full_name, &wrapper)?;
    let source = CString::new(source).expect("community wrapper contains a null byte");
    py.run(&source, Some(&wrapper.dict()), Some(&wrapper.dict()))?;
    community.add(name, &wrapper)?;
    PyModule::import(py, &full_name)
}

fn remove_wrapper(py: Python<'_>, community: &Bound<'_, PyModule>, name: &str) -> PyResult<()> {
    py.import("sys")?
        .getattr("modules")?
        .del_item(format!("symbolica.community.{name}"))?;
    community.delattr(name)
}

#[test]
fn feynkit_and_spenso_share_one_symbolica_kernel_in_both_import_orders() {
    Python::initialize();
    Python::attach(|py| -> PyResult<()> {
        let symbolica = install_package(py, "symbolica")?;
        let core = PyModule::new(py, "symbolica.core")?;
        create_symbolica_module(&core)?;
        py.import("sys")?
            .getattr("modules")?
            .set_item("symbolica.core", &core)?;
        symbolica.add("core", &core)?;

        let community = install_package(py, "symbolica.community")?;
        symbolica.add("community", &community)?;
        register_native::<FeynkitModule>(&core)?;
        register_native::<SpensoModule>(&core)?;

        for order in [["feynkit", "spenso"], ["spenso", "feynkit"]] {
            for name in order {
                let source = match name {
                    "feynkit" => FEYNKIT_WRAPPER,
                    "spenso" => SPENSO_WRAPPER,
                    _ => unreachable!(),
                };
                import_wrapper(py, &community, name, source)?;
            }

            let feynkit = PyModule::import(py, "symbolica.community.feynkit")?;
            let spenso = PyModule::import(py, "symbolica.community.spenso")?;
            let diagram = FeynmanDiagram::builder("shared-kernel")
                .overall_factor("x + 1")
                .build()
                .expect("empty test diagram is valid");
            let diagram = Py::new(py, PyFeynmanDiagram::from(diagram))?;
            let locals = PyDict::new(py);
            locals.set_item("core", &core)?;
            locals.set_item("diagram", diagram)?;
            locals.set_item("feynkit", feynkit)?;
            locals.set_item("spenso", spenso)?;
            let assertions = CString::new(
                r#"
feynkit_expression = diagram.overall_factor_expression()
spenso_expression = spenso.TensorName("T").to_expression()

assert diagram.__class__ is feynkit.FeynmanDiagram
assert type(feynkit_expression) is core.Expression
assert type(spenso_expression) is core.Expression
assert type(feynkit_expression + spenso_expression) is core.Expression
"#,
            )
            .unwrap();
            py.run(&assertions, Some(&locals), Some(&locals))?;

            remove_wrapper(py, &community, "feynkit")?;
            remove_wrapper(py, &community, "spenso")?;
        }
        Ok(())
    })
    .unwrap();
}
