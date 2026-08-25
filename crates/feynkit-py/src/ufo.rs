use std::path::PathBuf;

use feynkit_ufo::{LoadedModel, UfoLoadDiagnostics, UfoLoadOptions, UfoLoader};
use pyo3::{prelude::*, types::PyModule};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::{
    error,
    model::{PyModel, PyParameterCard},
};

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "UfoLoadDiagnostics",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyUfoLoadDiagnostics {
    inner: UfoLoadDiagnostics,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyUfoLoadDiagnostics {
    #[getter]
    fn source(&self) -> PathBuf {
        self.inner.source.clone()
    }

    #[getter]
    fn restriction_name(&self) -> Option<String> {
        self.inner.options.restriction_name.clone()
    }

    #[getter]
    fn simplify_model(&self) -> bool {
        self.inner.options.simplify_model
    }

    #[getter]
    fn wrap_indices_in_lorentz_structures(&self) -> bool {
        self.inner.options.wrap_indices_in_lorentz_structures
    }

    #[getter]
    fn order_count(&self) -> usize {
        self.inner.order_count
    }

    #[getter]
    fn model_parameter_count(&self) -> usize {
        self.inner.model_parameter_count
    }

    #[getter]
    fn particle_count(&self) -> usize {
        self.inner.particle_count
    }

    #[getter]
    fn propagator_count(&self) -> usize {
        self.inner.propagator_count
    }

    #[getter]
    fn lorentz_structure_count(&self) -> usize {
        self.inner.lorentz_structure_count
    }

    #[getter]
    fn coupling_count(&self) -> usize {
        self.inner.coupling_count
    }

    #[getter]
    fn vertex_rule_count(&self) -> usize {
        self.inner.vertex_rule_count
    }

    #[getter]
    fn function_count(&self) -> usize {
        self.inner.function_count
    }

    #[getter]
    fn form_factor_count(&self) -> usize {
        self.inner.form_factor_count
    }

    #[getter]
    fn parameter_value_count(&self) -> usize {
        self.inner.parameter_value_count
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "LoadedModel",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyLoadedModel {
    inner: LoadedModel,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyLoadedModel {
    #[getter]
    fn model(&self) -> PyModel {
        self.inner.model.clone().into()
    }

    #[getter]
    fn parameters(&self) -> PyParameterCard {
        self.inner.parameters.clone().into()
    }

    #[getter]
    fn diagnostics(&self) -> PyUfoLoadDiagnostics {
        PyUfoLoadDiagnostics {
            inner: self.inner.diagnostics.clone(),
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "UfoLoader",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyUfoLoader {
    inner: UfoLoader,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyUfoLoader {
    #[new]
    #[pyo3(signature = (*, restriction_name=None, simplify_model=true, wrap_indices_in_lorentz_structures=true))]
    fn new(
        restriction_name: Option<String>,
        simplify_model: bool,
        wrap_indices_in_lorentz_structures: bool,
    ) -> Self {
        Self {
            inner: UfoLoader::with_options(UfoLoadOptions {
                restriction_name,
                simplify_model,
                wrap_indices_in_lorentz_structures,
            }),
        }
    }

    #[getter]
    fn restriction_name(&self) -> Option<String> {
        self.inner.options().restriction_name.clone()
    }

    #[getter]
    fn simplify_model(&self) -> bool {
        self.inner.options().simplify_model
    }

    #[getter]
    fn wrap_indices_in_lorentz_structures(&self) -> bool {
        self.inner.options().wrap_indices_in_lorentz_structures
    }

    /// Load through the caller's attached interpreter without changing Python globals.
    fn load(&self, py: Python<'_>, path: PathBuf) -> PyResult<PyLoadedModel> {
        self.inner
            .load(py, path)
            .map(|inner| PyLoadedModel { inner })
            .map_err(error::ufo)
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyUfoLoadDiagnostics>()?;
    module.add_class::<PyLoadedModel>()?;
    module.add_class::<PyUfoLoader>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use pyo3::types::PyDict;

    use super::*;

    const MODEL_JSON: &str = include_str!("../tests/fixtures/scalars_2p_3p.json");

    #[test]
    fn exposes_frozen_diagnostics_from_the_normalized_loader_output() {
        Python::initialize();
        Python::attach(|py| {
            let module = PyModule::new(py, "symbolica.community.feynkit").unwrap();
            crate::initialize_feynkit(&module).unwrap();
            let locals = PyDict::new(py);
            locals.set_item("fk", &module).unwrap();
            locals.set_item("MODEL_JSON", MODEL_JSON).unwrap();
            let code = CString::new(
                r#"
import sys
import types

class JsonValue:
    def __init__(self, value):
        self.value = value

    def to_json(self):
        return self.value

def load_model(**kwargs):
    return JsonValue(MODEL_JSON), JsonValue('{"mass_scalar_0":[2.5,0.0]}')

package = types.ModuleType("ufo_model_loader")
commands = types.ModuleType("ufo_model_loader.commands")
commands.load_model = load_model
package.commands = commands
sys.modules["ufo_model_loader"] = package
sys.modules["ufo_model_loader.commands"] = commands

loaded = fk.UfoLoader(
    restriction_name="massless",
    simplify_model=False,
    wrap_indices_in_lorentz_structures=False,
).load("/models/scalars")
diagnostics = loaded.diagnostics
assert isinstance(diagnostics, fk.UfoLoadDiagnostics)
assert str(diagnostics.source) == "/models/scalars"
assert diagnostics.restriction_name == "massless"
assert diagnostics.simplify_model is False
assert diagnostics.wrap_indices_in_lorentz_structures is False
assert diagnostics.order_count == 2
assert diagnostics.model_parameter_count == 12
assert diagnostics.particle_count == 3
assert diagnostics.propagator_count == 3
assert diagnostics.lorentz_structure_count == 2
assert diagnostics.coupling_count == 1
assert diagnostics.vertex_rule_count == 16
assert diagnostics.function_count == 0
assert diagnostics.form_factor_count == 0
assert diagnostics.parameter_value_count == 1

try:
    diagnostics.particle_count = 0
except AttributeError:
    pass
else:
    raise AssertionError("UfoLoadDiagnostics is mutable")

del sys.modules["ufo_model_loader.commands"]
del sys.modules["ufo_model_loader"]
"#,
            )
            .unwrap();

            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }
}
