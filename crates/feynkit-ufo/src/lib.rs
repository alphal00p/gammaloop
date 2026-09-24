//! A thin adapter from Python's `ufo_model_loader` to [`feynkit_model`].
//!
//! The caller owns the Python interpreter. This crate neither initializes an
//! interpreter nor changes Python logging, standard output, or standard error.

#![forbid(unsafe_code)]

use std::path::{Path, PathBuf};

use feynkit_model::{Model, ModelError, ParameterCard};
use pyo3::{
    Py, PyAny, PyErr, Python,
    prelude::PyAnyMethods,
    types::{PyDict, PyDictMethods, PyModule},
};
use thiserror::Error;

/// Options forwarded to `ufo_model_loader.commands.load_model`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UfoLoadOptions {
    pub restriction_name: Option<String>,
    pub simplify_model: bool,
    pub wrap_indices_in_lorentz_structures: bool,
}

impl Default for UfoLoadOptions {
    fn default() -> Self {
        Self {
            restriction_name: None,
            simplify_model: true,
            wrap_indices_in_lorentz_structures: true,
        }
    }
}

/// Configured adapter for Python's `ufo_model_loader` package.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct UfoLoader {
    options: UfoLoadOptions,
}

impl UfoLoader {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn with_options(options: UfoLoadOptions) -> Self {
        Self { options }
    }

    pub fn options(&self) -> &UfoLoadOptions {
        &self.options
    }

    pub fn restriction_name(mut self, restriction_name: impl Into<String>) -> Self {
        self.options.restriction_name = Some(restriction_name.into());
        self
    }

    pub fn simplify_model(mut self, simplify_model: bool) -> Self {
        self.options.simplify_model = simplify_model;
        self
    }

    pub fn wrap_indices_in_lorentz_structures(mut self, wrap: bool) -> Self {
        self.options.wrap_indices_in_lorentz_structures = wrap;
        self
    }

    /// Load a UFO model using an already-attached Python interpreter.
    pub fn load(
        &self,
        py: Python<'_>,
        path: impl AsRef<Path>,
    ) -> Result<LoadedModel, UfoLoadError> {
        load_model(py, path.as_ref(), &self.options)
    }
}

/// A model and the parameter values returned by the UFO loader.
#[derive(Clone, Debug)]
pub struct LoadedModel {
    pub model: Model,
    pub parameters: ParameterCard,
    pub diagnostics: UfoLoadDiagnostics,
}

/// Requested loader inputs and counts from the normalized output.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UfoLoadDiagnostics {
    pub source: PathBuf,
    pub options: UfoLoadOptions,
    pub order_count: usize,
    pub model_parameter_count: usize,
    pub particle_count: usize,
    pub propagator_count: usize,
    pub lorentz_structure_count: usize,
    pub coupling_count: usize,
    pub vertex_rule_count: usize,
    pub function_count: usize,
    pub form_factor_count: usize,
    pub parameter_value_count: usize,
}

/// Errors produced while loading a UFO model.
#[derive(Debug, Error)]
pub enum UfoLoadError {
    #[error("UFO model path '{}' is not valid UTF-8", path.display())]
    NonUtf8Path { path: PathBuf },
    #[error("could not import 'ufo_model_loader.commands': {source}")]
    Import {
        #[source]
        source: PyErr,
    },
    #[error("could not resolve 'ufo_model_loader.commands.load_model': {source}")]
    ResolveLoader {
        #[source]
        source: PyErr,
    },
    #[error("could not construct arguments for 'ufo_model_loader.commands.load_model': {source}")]
    Arguments {
        #[source]
        source: PyErr,
    },
    #[error("'ufo_model_loader.commands.load_model' failed: {source}")]
    Load {
        #[source]
        source: PyErr,
    },
    #[error("invalid output from 'ufo_model_loader.commands.load_model': {source}")]
    Output {
        #[source]
        source: PyErr,
    },
    #[error("could not deserialize the loaded model: {source}")]
    ModelJson {
        #[source]
        source: ModelError,
    },
    #[error("could not deserialize the loaded parameter card: {source}")]
    ParameterCardJson {
        #[source]
        source: ModelError,
    },
}

fn load_model(
    py: Python<'_>,
    path: &Path,
    options: &UfoLoadOptions,
) -> Result<LoadedModel, UfoLoadError> {
    let input_model_path = path.to_str().ok_or_else(|| UfoLoadError::NonUtf8Path {
        path: path.to_path_buf(),
    })?;
    let commands = PyModule::import(py, "ufo_model_loader.commands")
        .map_err(|source| UfoLoadError::Import { source })?;
    let loader = commands
        .getattr("load_model")
        .map_err(|source| UfoLoadError::ResolveLoader { source })?;

    // kwargs: input_model_path, restriction_name, simplify_model,
    // wrap_indices_in_lorentz_structures
    let kwargs = PyDict::new(py);
    kwargs
        .set_item("input_model_path", input_model_path)
        .map_err(|source| UfoLoadError::Arguments { source })?;
    match options.restriction_name.as_deref() {
        Some(name) => kwargs.set_item("restriction_name", name),
        None => kwargs.set_item("restriction_name", py.None()),
    }
    .and_then(|_| kwargs.set_item("simplify_model", options.simplify_model))
    .and_then(|_| {
        kwargs.set_item(
            "wrap_indices_in_lorentz_structures",
            options.wrap_indices_in_lorentz_structures,
        )
    })
    .map_err(|source| UfoLoadError::Arguments { source })?;

    // Expect a 2-tuple: (model, input_param_card).
    let output = loader
        .call((), Some(&kwargs))
        .map_err(|source| UfoLoadError::Load { source })?;
    let (py_model, py_card): (Py<PyAny>, Py<PyAny>) = output
        .extract()
        .map_err(|source| UfoLoadError::Output { source })?;

    // to_json() -> String
    let mut model_json = to_json(py, &py_model)?;
    let card_json = to_json(py, &py_card)?;

    // Preserve metadata that the upstream loader does not serialize. Read from
    // the imported UFO module, retaining unknown hypercharge as null, not zero.
    let metadata = (|| -> Result<Option<String>, PyErr> {
        let module_name = path
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("invalid UFO module path"))?;
        let modules = py.import("sys")?.getattr("modules")?;
        let modules = modules.cast::<PyDict>()?;
        let Some(module) = modules.get_item(module_name)? else {
            return Ok(None);
        };
        let metadata = PyDict::new(py);
        if module.hasattr("covariant_cut_multiplets")? {
            metadata.set_item(
                "covariant_cut_multiplets",
                module.getattr("covariant_cut_multiplets")?,
            )?;
        }
        let particles = PyDict::new(py);
        if module.hasattr("all_particles")? {
            let fraction = py.import("fractions")?.getattr("Fraction")?;
            for particle in module.getattr("all_particles")?.try_iter()? {
                let particle = particle?;
                let numbers = PyDict::new(py);
                for (field, attribute) in [
                    ("charge", "charge"),
                    ("y_charge", "Y"),
                    ("y_charge_right", "YRight"),
                ] {
                    let attribute = if field == "charge" && particle.hasattr("charge_exact")? {
                        "charge_exact"
                    } else {
                        attribute
                    };
                    if particle.hasattr(attribute)? && !particle.getattr(attribute)?.is_none() {
                        // Fraction strings stay exact. Existing float-valued UFO
                        // metadata retains its decimal value; never guess a charge
                        // from PDG codes or silently round to a small denominator.
                        let value = particle.getattr(attribute)?.str()?;
                        let exact = fraction.call1((value,))?.str()?;
                        numbers.set_item(field, exact)?;
                    } else {
                        numbers.set_item(field, py.None())?;
                    }
                }
                particles.set_item(particle.getattr("name")?, numbers)?;
            }
        }
        metadata.set_item("particles", particles)?;
        py.import("json")?
            .call_method1("dumps", (metadata,))?
            .extract()
            .map(Some)
    })()
    .map_err(|source| UfoLoadError::Output { source })?;
    if let Some(metadata) = metadata {
        let enrich = || -> Result<String, serde_json::Error> {
            let mut serialized: serde_json::Value = serde_json::from_str(&model_json)?;
            let metadata: serde_json::Value = serde_json::from_str(&metadata)?;
            if let Some(declaration) = metadata.get("covariant_cut_multiplets") {
                serialized["covariant_cut_multiplets"] = declaration.clone();
            }
            if let Some(particles) = serialized["particles"].as_array_mut() {
                for particle in particles {
                    if let Some(numbers) =
                        metadata["particles"].get(particle["name"].as_str().unwrap_or_default())
                    {
                        for field in ["charge", "y_charge", "y_charge_right"] {
                            particle[field] = numbers[field].clone();
                        }
                    }
                }
            }
            serde_json::to_string(&serialized)
        };
        model_json = enrich().map_err(|source| UfoLoadError::ModelJson {
            source: source.into(),
        })?;
    }

    // Deserialize into the standalone Rust types.
    let model =
        Model::from_json(&model_json).map_err(|source| UfoLoadError::ModelJson { source })?;
    let parameters = ParameterCard::from_json(&card_json)
        .map_err(|source| UfoLoadError::ParameterCardJson { source })?;

    let diagnostics = UfoLoadDiagnostics {
        source: path.to_path_buf(),
        options: options.clone(),
        order_count: model.orders().len(),
        model_parameter_count: model.parameters().len(),
        particle_count: model.particles().len(),
        propagator_count: model.propagators().len(),
        lorentz_structure_count: model.lorentz_structures().len(),
        coupling_count: model.couplings().len(),
        vertex_rule_count: model.vertex_rules().len(),
        function_count: model.functions().len(),
        form_factor_count: model.form_factors().len(),
        parameter_value_count: parameters.len(),
    };

    Ok(LoadedModel {
        model,
        parameters,
        diagnostics,
    })
}

fn to_json(py: Python<'_>, value: &Py<PyAny>) -> Result<String, UfoLoadError> {
    value
        .call_method0(py, "to_json")
        .and_then(|json| json.extract(py))
        .map_err(|source| UfoLoadError::Output { source })
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use pyo3::types::{PyDict, PyModule};

    use super::*;

    const MODEL_JSON: &str = r#"{
        "name":"scalar","restriction":null,
        "orders":[{"name":"QED","expansion_order":99,"hierarchy":1}],
        "parameters":[
            {"name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal","parameter_type":"real","value":[0.0,0.0],"expression":null},
            {"name":"mass","lhablock":"MASS","lhacode":[1],"nature":"external","parameter_type":"real","value":[1.0,0.0],"expression":null}
        ],
        "particles":[{"pdg_code":1,"name":"s","antiname":"s","spin":1,"color":1,"mass":"mass","width":"ZERO","texname":"s","antitexname":"s","charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0}],
        "propagators":[],"lorentz_structures":[],"couplings":[],"vertex_rules":[]
    }"#;

    #[test]
    fn defaults_match_the_existing_loader_contract() {
        assert_eq!(UfoLoader::new(), UfoLoader::default());
        assert_eq!(
            UfoLoader::with_options(UfoLoadOptions::default()).options(),
            &UfoLoadOptions {
                restriction_name: None,
                simplify_model: true,
                wrap_indices_in_lorentz_structures: true,
            }
        );
    }

    #[test]
    fn forwards_options_and_deserializes_json() {
        Python::initialize();
        Python::attach(|py| {
            let code = format!(
                r#"
import sys
import types

MODEL_JSON = {MODEL_JSON:?}
CARD_JSON = '{{"mass":[2.5,0.0]}}'
seen = {{}}

class JsonValue:
    def __init__(self, value):
        self.value = value

    def to_json(self):
        return self.value

def load_model(**kwargs):
    seen.update(kwargs)
    return JsonValue(MODEL_JSON), JsonValue(CARD_JSON)

package = types.ModuleType("ufo_model_loader")
commands = types.ModuleType("ufo_model_loader.commands")
commands.load_model = load_model
commands.seen = seen
package.commands = commands
sys.modules["ufo_model_loader"] = package
sys.modules["ufo_model_loader.commands"] = commands
"#
            );
            PyModule::from_code(
                py,
                CString::new(code).unwrap().as_c_str(),
                c"feynkit_ufo_test.py",
                c"feynkit_ufo_test",
            )
            .unwrap();

            let loader = UfoLoader::new()
                .restriction_name("massless")
                .simplify_model(false)
                .wrap_indices_in_lorentz_structures(false);
            let loaded = loader.load(py, "/models/scalar").unwrap();
            assert_eq!(loaded.model.name(), "scalar");
            assert_eq!(loaded.parameters["mass"].re, 2.5);
            assert!(loaded.model.particle("s").unwrap().is_propagating());
            assert!(loaded.model.functions().is_empty());
            assert!(loaded.model.form_factors().is_empty());
            assert_eq!(
                loaded.diagnostics,
                UfoLoadDiagnostics {
                    source: PathBuf::from("/models/scalar"),
                    options: loader.options().clone(),
                    order_count: 1,
                    model_parameter_count: 2,
                    particle_count: 1,
                    propagator_count: 0,
                    lorentz_structure_count: 0,
                    coupling_count: 0,
                    vertex_rule_count: 0,
                    function_count: 0,
                    form_factor_count: 0,
                    parameter_value_count: 1,
                }
            );

            let commands = PyModule::import(py, "ufo_model_loader.commands").unwrap();
            let seen = commands
                .getattr("seen")
                .unwrap()
                .cast_into::<PyDict>()
                .unwrap();
            assert_eq!(
                seen.get_item("input_model_path")
                    .unwrap()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "/models/scalar"
            );
            assert_eq!(
                seen.get_item("restriction_name")
                    .unwrap()
                    .unwrap()
                    .extract::<String>()
                    .unwrap(),
                "massless"
            );
            assert!(
                !seen
                    .get_item("simplify_model")
                    .unwrap()
                    .unwrap()
                    .extract::<bool>()
                    .unwrap()
            );
            assert!(
                !seen
                    .get_item("wrap_indices_in_lorentz_structures")
                    .unwrap()
                    .unwrap()
                    .extract::<bool>()
                    .unwrap()
            );

            let modules = PyModule::import(py, "sys")
                .unwrap()
                .getattr("modules")
                .unwrap()
                .cast_into::<PyDict>()
                .unwrap();
            let source = PyModule::from_code(
                py,
                c"from types import SimpleNamespace\nall_particles = [SimpleNamespace(name='s', charge=0.6666666666666666, charge_exact='2/3', Y='1/3', YRight='4/3')]",
                c"scalar.py", c"scalar",
            ).unwrap();
            modules.set_item("scalar", &source).unwrap();
            let exact = loader.load(py, "/models/scalar").unwrap();
            let particle = exact.model.particle("s").unwrap();
            assert_eq!(particle.charge.to_string(), "2/3");
            assert_eq!(particle.y_charge.as_ref().unwrap().to_string(), "1/3");
            assert_eq!(particle.y_charge_right.as_ref().unwrap().to_string(), "4/3");
            let particle = source
                .getattr("all_particles")
                .unwrap()
                .get_item(0)
                .unwrap();
            particle.delattr("Y").unwrap();
            particle.setattr("YRight", py.None()).unwrap();
            let unknown = loader.load(py, "/models/scalar").unwrap();
            assert!(
                unknown
                    .model
                    .particle("s")
                    .unwrap()
                    .weak_isospin()
                    .is_none()
            );
            assert!(
                unknown
                    .model
                    .particle("s")
                    .unwrap()
                    .weak_isospin_right()
                    .is_none()
            );
            modules.del_item("scalar").unwrap();
            modules.del_item("ufo_model_loader.commands").unwrap();
            modules.del_item("ufo_model_loader").unwrap();
        });
    }
}
