use std::path::PathBuf;

use feynkit_ufo::{LoadedModel, UfoLoadDiagnostics, UfoLoadOptions, UfoLoader};
use pyo3::{
    prelude::*,
    types::{PyAny, PyModule},
    wrap_pyfunction,
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyfunction, gen_stub_pymethods};

use crate::{
    display::escape_html,
    error,
    model::{PyModel, PyParameterCard},
};

/// Provenance, options, and entity counts from loading a UFO model.
///
/// Diagnostics make it easy to verify which restriction was applied and how
/// many particles, interactions, and auxiliary structures were normalized.
///
/// Examples
/// --------
/// >>> diagnostics = loaded.diagnostics
/// >>> print(diagnostics.source, diagnostics.particle_count, diagnostics.vertex_rule_count)
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
    /// Return the UFO source directory.
    #[getter]
    fn source(&self) -> PathBuf {
        self.inner.source.clone()
    }

    /// Return the applied restriction-card name, if any.
    #[getter]
    fn restriction_name(&self) -> Option<String> {
        self.inner.options.restriction_name.clone()
    }

    /// Return whether model expressions were simplified while loading.
    #[getter]
    fn simplify_model(&self) -> bool {
        self.inner.options.simplify_model
    }

    /// Return whether normalized Lorentz indices were wrapped.
    #[getter]
    fn wrap_indices_in_lorentz_structures(&self) -> bool {
        self.inner.options.wrap_indices_in_lorentz_structures
    }

    /// Return the number of coupling orders loaded from the UFO model.
    #[getter]
    fn order_count(&self) -> usize {
        self.inner.order_count
    }

    /// Return the number of model parameters loaded.
    #[getter]
    fn model_parameter_count(&self) -> usize {
        self.inner.model_parameter_count
    }

    /// Return the number of particles loaded.
    #[getter]
    fn particle_count(&self) -> usize {
        self.inner.particle_count
    }

    /// Return the number of propagators loaded.
    #[getter]
    fn propagator_count(&self) -> usize {
        self.inner.propagator_count
    }

    /// Return the number of Lorentz structures loaded.
    #[getter]
    fn lorentz_structure_count(&self) -> usize {
        self.inner.lorentz_structure_count
    }

    /// Return the number of couplings loaded.
    #[getter]
    fn coupling_count(&self) -> usize {
        self.inner.coupling_count
    }

    /// Return the number of vertex rules loaded.
    #[getter]
    fn vertex_rule_count(&self) -> usize {
        self.inner.vertex_rule_count
    }

    /// Return the number of model functions loaded.
    #[getter]
    fn function_count(&self) -> usize {
        self.inner.function_count
    }

    /// Return the number of form factors loaded.
    #[getter]
    fn form_factor_count(&self) -> usize {
        self.inner.form_factor_count
    }

    /// Return the number of parameter values loaded into the parameter card.
    #[getter]
    fn parameter_value_count(&self) -> usize {
        self.inner.parameter_value_count
    }

    /// Return a concise summary of the UFO import diagnostics.
    ///
    /// Examples
    /// --------
    /// >>> print(loaded.diagnostics)
    ///
    fn __repr__(&self) -> String {
        format!(
            "UfoLoadDiagnostics(source={:?}, particles={}, vertices={})",
            self.inner.source, self.inner.particle_count, self.inner.vertex_rule_count,
        )
    }

    /// Render the UFO source and imported object counts as a notebook table.
    ///
    /// Examples
    /// --------
    /// Leave ``loaded.diagnostics`` as the final expression in a notebook cell
    /// to display the import inventory.
    ///
    fn _repr_html_(&self) -> String {
        let restriction = self
            .inner
            .options
            .restriction_name
            .as_deref()
            .map_or_else(|| "none".to_owned(), escape_html);
        format!(
            "<div class=\"feynkit-ufo-diagnostics\" style=\"display:inline-block;max-width:\
             100%;overflow-x:auto\"><strong>UFO import</strong><div style=\"opacity:.75\">\
             {} · restriction: {restriction}</div><table style=\"border-collapse:collapse;\
             margin-top:.3rem\"><thead><tr><th style=\"padding:.2rem .55rem;text-align:right\">\
             particles</th><th style=\"padding:.2rem .55rem;text-align:right\">parameters</th>\
             <th style=\"padding:.2rem .55rem;text-align:right\">couplings</th><th style=\"\
             padding:.2rem .55rem;text-align:right\">vertex rules</th><th style=\"padding:.2rem \
             .55rem;text-align:right\">Lorentz structures</th></tr></thead><tbody><tr><td \
             style=\"padding:.2rem .55rem;text-align:right\">{}</td><td style=\"padding:.2rem \
             .55rem;text-align:right\">{}</td><td style=\"padding:.2rem .55rem;text-align:\
             right\">{}</td><td style=\"padding:.2rem .55rem;text-align:right\">{}</td><td \
             style=\"padding:.2rem .55rem;text-align:right\">{}</td></tr></tbody></table></div>",
            escape_html(&self.inner.source.to_string_lossy()),
            self.inner.particle_count,
            self.inner.model_parameter_count,
            self.inner.coupling_count,
            self.inner.vertex_rule_count,
            self.inner.lorentz_structure_count,
        )
    }

    /// Write the concise diagnostic summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython invokes this method when only a text representation is supported.
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     The IPython pretty-printer object.
    /// cycle : bool
    ///     Whether this object is part of a recursive formatting cycle.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "...".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

/// A normalized model, its parameter card, and its loading diagnostics.
///
/// This is the complete result returned by ``UfoLoader.load``; use ``model``
/// for diagram generation and ``parameters`` for external numerical inputs.
///
/// Examples
/// --------
/// >>> loaded = fk.UfoLoader().load("path/to/MyUFO")
/// >>> generator = fk.Generator(loaded.model)
/// >>> parameters = loaded.parameters
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
    /// Return the normalized FeynKit model.
    ///
    /// Examples
    /// --------
    /// >>> generator = fk.Generator(loaded.model)
    ///
    #[getter]
    fn model(&self) -> PyModel {
        self.inner.model.clone().into()
    }

    /// Return the parameter card loaded alongside the model.
    ///
    /// Examples
    /// --------
    /// >>> parameters = loaded.parameters
    /// >>> shifted = loaded.model.with_parameter_card(parameters)
    ///
    #[getter]
    fn parameters(&self) -> PyParameterCard {
        self.inner.parameters.clone().into()
    }

    /// Return counts and options recorded while loading the UFO model.
    #[getter]
    fn diagnostics(&self) -> PyUfoLoadDiagnostics {
        PyUfoLoadDiagnostics {
            inner: self.inner.diagnostics.clone(),
        }
    }

    /// Return a concise summary of the loaded model and its source.
    ///
    /// Examples
    /// --------
    /// >>> print(loaded)
    ///
    fn __repr__(&self) -> String {
        format!(
            "LoadedModel(name={:?}, source={:?}, particles={})",
            self.inner.model.name(),
            self.inner.diagnostics.source,
            self.inner.model.particles().len(),
        )
    }

    /// Render a compact overview of the normalized UFO model.
    ///
    /// Examples
    /// --------
    /// Leave ``loaded`` as the final expression in a notebook cell to display
    /// the normalized model and import source.
    ///
    fn _repr_html_(&self) -> String {
        format!(
            "<div class=\"feynkit-loaded-model\" style=\"display:inline-block;max-width:100%;\
             overflow-x:auto\"><strong>{}</strong><div style=\"opacity:.75\">loaded from {}\
             </div><table style=\"border-collapse:collapse;margin-top:.3rem\"><thead><tr><th \
             style=\"padding:.2rem .6rem;text-align:right\">particles</th><th style=\"padding:\
             .2rem .6rem;text-align:right\">parameters</th><th style=\"padding:.2rem .6rem;\
             text-align:right\">couplings</th><th style=\"padding:.2rem .6rem;text-align:right\">\
             vertex rules</th></tr></thead><tbody><tr><td style=\"padding:.2rem .6rem;text-align:\
             right\">{}</td><td style=\"padding:.2rem .6rem;text-align:right\">{}</td><td \
             style=\"padding:.2rem .6rem;text-align:right\">{}</td><td style=\"padding:.2rem \
             .6rem;text-align:right\">{}</td></tr></tbody></table></div>",
            escape_html(self.inner.model.name()),
            escape_html(&self.inner.diagnostics.source.to_string_lossy()),
            self.inner.model.particles().len(),
            self.inner.model.parameters().len(),
            self.inner.model.couplings().len(),
            self.inner.model.vertex_rules().len(),
        )
    }

    /// Write the concise loaded-model summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython invokes this method when only a text representation is supported.
    ///
    /// Parameters
    /// ----------
    /// pretty : object
    ///     The IPython pretty-printer object.
    /// cycle : bool
    ///     Whether this object is part of a recursive formatting cycle.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        pretty.call_method1(
            "text",
            (if cycle {
                "...".to_owned()
            } else {
                self.__repr__()
            },),
        )?;
        Ok(())
    }
}

/// Load and normalize a Universal FeynRules Output model.
///
/// The loader delegates UFO parsing to ``ufo_model_loader`` and returns typed
/// FeynKit model entities suitable for process and diagram generation.
///
/// Examples
/// --------
/// >>> import symbolica.community.feynkit as fk
/// >>> loader = fk.UfoLoader(restriction_name="massless")
/// >>> loaded = loader.load("path/to/MyUFO")
///
/// Parameters
/// ----------
/// restriction_name : str, optional
///     Restriction-card name to apply while loading.
/// simplify_model : bool, optional
///     Simplify normalized symbolic model expressions.
/// wrap_indices_in_lorentz_structures : bool, optional
///     Wrap normalized indices in Lorentz structures.
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
    /// Construct a UFO loader with normalization options.
    ///
    /// Examples
    /// --------
    /// >>> loader = UfoLoader(restriction_name="massless")
    ///
    /// Parameters
    /// ----------
    /// restriction_name : str, optional
    ///     Restriction-card name passed to the UFO loader.
    /// simplify_model : bool, optional
    ///     Simplify expressions in the normalized model.
    /// wrap_indices_in_lorentz_structures : bool, optional
    ///     Wrap indices in normalized Lorentz structures.
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

    /// Return the configured restriction-card name.
    #[getter]
    fn restriction_name(&self) -> Option<String> {
        self.inner.options().restriction_name.clone()
    }

    /// Return whether model simplification is enabled.
    #[getter]
    fn simplify_model(&self) -> bool {
        self.inner.options().simplify_model
    }

    /// Return whether normalized Lorentz indices are wrapped.
    #[getter]
    fn wrap_indices_in_lorentz_structures(&self) -> bool {
        self.inner.options().wrap_indices_in_lorentz_structures
    }

    /// Load through the caller's attached interpreter without changing Python globals.
    ///
    /// Examples
    /// --------
    /// >>> loaded = UfoLoader().load("path/to/ufo_model")
    /// >>> model = loaded.model
    ///
    /// Parameters
    /// ----------
    /// path : str or os.PathLike
    ///     Directory containing the UFO Python modules.
    fn load(&self, py: Python<'_>, path: PathBuf) -> PyResult<PyLoadedModel> {
        self.inner
            .load(py, path)
            .map(|inner| PyLoadedModel { inner })
            .map_err(error::ufo)
    }
}

/// Load and normalize a Universal FeynRules Output model.
///
/// The returned native :class:`LoadedModel` bundles the normalized model,
/// parameter card, and diagnostics. Use :meth:`Model.from_path` instead when a
/// model has already been normalized to FeynKit JSON.
///
/// Examples
/// --------
/// Load a UFO directory and inspect its particle content:
///
/// >>> loaded = fk.load_ufo_model("models/sm", restriction="massless")
/// >>> electron = loaded.model.particle("e-")
/// >>> electron.pdg_code
/// 11
///
/// Parameters
/// ----------
/// path : str or os.PathLike
///     Directory containing the UFO Python modules.
/// restriction : str or None, optional
///     Restriction-card name to apply while loading.
/// simplify : bool, optional
///     Simplify normalized symbolic model expressions.
/// wrap_lorentz_indices : bool, optional
///     Wrap normalized indices in Lorentz structures.
#[cfg_attr(
    feature = "python_stubgen",
    gen_stub_pyfunction(module = "symbolica.community.feynkit")
)]
#[pyfunction]
#[pyo3(signature = (path, *, restriction=None, simplify=true, wrap_lorentz_indices=true))]
fn load_ufo_model(
    py: Python<'_>,
    path: PathBuf,
    restriction: Option<String>,
    simplify: bool,
    wrap_lorentz_indices: bool,
) -> PyResult<PyLoadedModel> {
    UfoLoader::with_options(UfoLoadOptions {
        restriction_name: restriction,
        simplify_model: simplify,
        wrap_indices_in_lorentz_structures: wrap_lorentz_indices,
    })
    .load(py, path)
    .map(|inner| PyLoadedModel { inner })
    .map_err(error::ufo)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyUfoLoadDiagnostics>()?;
    module.add_class::<PyLoadedModel>()?;
    module.add_class::<PyUfoLoader>()?;
    let load_ufo_model = wrap_pyfunction!(load_ufo_model, module)?;
    load_ufo_model.setattr("__module__", "symbolica.community.feynkit")?;
    module.add_function(load_ufo_model)?;
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
