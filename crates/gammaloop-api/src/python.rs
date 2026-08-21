use gammaloop_tracing_filter::LogFormat;
use gammalooprs::{
    cff::generation::{generate_cff_expression_from_subgraph, SurfaceCache},
    feyngen::diagram_generator::evaluate_overall_factor,
    graph::{self, FeynmanGraph, Graph, LMBext},
    integrands::evaluation::{
        BatchSampleEvaluationResult, SampleEvaluationResult, SingleSampleEvaluationResult,
    },
    observables::{
        AdditionalWeightKey, DiscreteBinOrdering, Event, EventGroup, GenericAdditionalWeightInfo,
        HistogramAccumulatorState, HistogramSnapshot, HistogramStatisticsSnapshot,
    },
    processes::{DotExportSettings, ProcessCollection},
    settings::{global::OrientationPattern, RuntimeSettings},
    utils::tracing::LogLevel,
};
use idenso::shorthands::{metric::to_dots_impl, schoonschip::Schoonschip};
use linnet::half_edge::{
    involution::{EdgeIndex, Orientation},
    subgraph::{ModifySubSet, SuBitGraph},
};
use numpy::{PyReadonlyArray1, PyReadonlyArray2};
use typed_index_collections::TiVec;

use crate::{
    commands::{
        evaluate_samples::EvaluateSamples, import::model::ImportModel, integrate::Integrate,
        Evaluate,
    },
    integrand_info::{
        IntegrandActiveThresholdCutInfo, IntegrandCutInfo, IntegrandCutThresholdInfo,
        IntegrandGraphGroupInfo, IntegrandGraphInfo, IntegrandInfo, IntegrandLoopMomentumBasisInfo,
        IntegrandOrientationInfo, IntegrandThresholdEsurfaceInfo,
    },
    render_smart_toml,
    session::{display_command, CliSession, CliSessionState},
    settings_tree::{json_type_name, serialize_settings_with_defaults, value_at_path},
    state::{ProcessListExt, ProcessRef, RunHistory, State},
    CLISettings, LoadedState, StateLoadOption,
};
use ahash::{HashMap, HashMapExt};
use clap::ValueEnum;

use color_eyre::Result;
use eyre::eyre;
use serde::Serialize;
use serde_json::Value as JsonValue;

/*
use gammalooprs::feyngen::{
    FeynGenError, FeynGenFilter, FeynGenFilters, GenerationType, SelfEnergyFilterOptions,
    SnailFilterOptions, TadpolesFilterOptions,
};
*/

// use git_version::git_version;
use itertools::{self, Itertools};
use std::{path::PathBuf, str::FromStr};

use symbolica::{atom::AtomCore, parse, printer::PrintOptions};
// const GIT_VERSION: &str = git_version!(fallback = "unavailable");

#[allow(unused)]
use pyo3::{
    exceptions,
    prelude::*,
    pyclass,
    pyclass::CompareOp,
    pyfunction, pymethods, pymodule,
    types::{PyAny, PyComplex, PyComplexMethods, PyDict, PyList, PyModule, PyTuple, PyType},
    wrap_pyfunction, FromPyObject, PyRef, Python,
};

/// Evaluate a graph's symbolic overall factor and return its canonical form.
///
/// Parameters
/// ----------
/// overall_factor : str
///     Symbolica expression used as the graph's overall factor.
///
/// Returns
/// -------
/// str
///     The evaluated factor in canonical Symbolica form.
///
/// Raises
/// ------
/// Exception
///     If the expression cannot be parsed or evaluated.
///
/// Examples
/// --------
/// Evaluate the symmetry and fermion-loop factors attached to a graph::
///
///     evaluate_graph_overall_factor(
///         "AutG(2)^-1*InternalFermionLoopSign(-1)"
///     )
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(name = "evaluate_graph_overall_factor")]
pub(crate) fn evaluate_graph_overall_factor(overall_factor: &str) -> Result<String> {
    let overall_factor = parse!(overall_factor);
    let overall_factor_evaluated = evaluate_overall_factor(overall_factor.as_view());
    Ok(overall_factor_evaluated.to_canonical_string())
}

/// Parse a Symbolica expression and return its canonical string form.
///
/// Parameters
/// ----------
/// atom_str : str
///     Symbolica expression to parse and normalize.
///
/// Returns
/// -------
/// str
///     The parsed expression in canonical string form.
///
/// Raises
/// ------
/// Exception
///     If ``atom_str`` is not a valid Symbolica expression.
///
/// Examples
/// --------
/// Normalize a symbolic expression before comparing or storing it::
///
///     atom_to_canonical_string("x + x")
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(name = "atom_to_canonical_string")]
pub(crate) fn atom_to_canonical_string(atom_str: &str) -> Result<String> {
    Ok(parse!(atom_str).to_canonical_string())
}

/// Rewrite a Symbolica tensor expression into Idenso dot-product notation.
///
/// Parameters
/// ----------
/// atom_str : str
///     Symbolica tensor expression whose repeated indices encode contractions.
///
/// Returns
/// -------
/// str
///     The equivalent expression using Idenso dot-product notation.
///
/// Raises
/// ------
/// Exception
///     If ``atom_str`` cannot be parsed or converted.
///
/// Examples
/// --------
/// Rewrite a contraction before passing it to an Idenso workflow::
///
///     to_dots("p(mu) * q(mu)")
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
#[pyo3(name = "to_dots")]
pub(crate) fn atom_to_dots(atom_str: &str) -> Result<String> {
    let dotted = to_dots_impl(
        parse!(atom_str, default_namespace = "python")
            .to_dots()
            .as_view(),
    );
    Ok(format!(
        "{}",
        dotted.as_view().printer(PrintOptions {
            hide_namespace: Some(std::borrow::Cow::Borrowed("python")),
            ..PrintOptions::file()
        })
    ))
}

#[pymodule(name = "_gammaloop")]
fn python_module(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    gammalooprs::initialisation::initialise().expect("initialization failed");
    gammalooprs::set_interrupt_handler();
    register_python_api(m)
}

fn register_python_api(m: &Bound<PyModule>) -> PyResult<()> {
    m.add_class::<GammaLoopAPI>()?;
    m.add_class::<LogLevel>()?;
    m.add_class::<DotExportSettings>()?;
    m.add_class::<PyEvaluationResult>()?;
    m.add_class::<PyBatchEvaluationResult>()?;
    m.add_class::<PySampleEvaluationResult>()?;
    m.add_class::<PyEventGroup>()?;
    m.add_class::<PyEvent>()?;
    m.add_class::<PyFourMomentum>()?;
    m.add_class::<PyCutInfo>()?;
    m.add_class::<PyIntegrandInfo>()?;
    m.add_class::<PyIntegrandGraphGroupInfo>()?;
    m.add_class::<PyIntegrandGraphInfo>()?;
    m.add_class::<PyIntegrandOrientationInfo>()?;
    m.add_class::<PyIntegrandLoopMomentumBasisInfo>()?;
    m.add_class::<PyIntegrandCutInfo>()?;
    m.add_class::<PyIntegrandCutThresholdInfo>()?;
    m.add_class::<PyIntegrandActiveThresholdCutInfo>()?;
    m.add_class::<PyIntegrandThresholdEsurfaceInfo>()?;
    m.add_class::<PyAdditionalWeight>()?;
    m.add_class::<PyHistogramAccumulator>()?;
    m.add_class::<PyHistogramSnapshot>()?;
    m.add_class::<PyHistogramBinSnapshot>()?;
    m.add_class::<PyHistogramStatisticsSnapshot>()?;
    m.add_class::<PyIntegrationOutput>()?;
    m.add_class::<PyIntegralEstimate>()?;
    m.add_class::<PyIntegrationTableComponentResult>()?;
    m.add_class::<PyIntegrationStatisticsSnapshot>()?;
    m.add_class::<PyMaxWeightInfoEntry>()?;
    m.add_class::<PyDiscreteCoordinate>()?;
    m.add_class::<PyDiscreteBreakdownEntry>()?;
    m.add_class::<PyDiscreteBreakdown>()?;
    m.add_class::<PyComponentDiscreteBreakdown>()?;
    m.add_class::<PySlotIntegrationResult>()?;
    m.add_class::<PyIntegrationResult>()?;
    m.add_class::<PyStabilityResult>()?;
    m.add_class::<PySettingsValue>()?;
    /*
    m.add_class::<PyFeynGenFilters>()?;
    m.add_class::<PySnailFilterOptions>()?;
    m.add_class::<PySewedFilterOptions>()?;
    m.add_class::<PyTadpolesFilterOptions>()?;
    m.add_class::<PySelfEnergyFilterOptions>()?;
    m.add_class::<GlobalSettings>()?;
    m.add_class::<RuntimeSettings>()?;

    // m.add_class::<PyFeynGenOptions>()?;
    m.add_class::<PyNumeratorAwareGroupingOption>()?;
    */
    // m.add("git_version", GIT_VERSION)?;
    m.add_wrapped(wrap_pyfunction!(atom_to_canonical_string))?;
    m.add_wrapped(wrap_pyfunction!(atom_to_dots))?;
    m.add_wrapped(wrap_pyfunction!(evaluate_graph_overall_factor))?;
    Ok(())
}

#[cfg(feature = "python_stubgen")]
#[doc(hidden)]
pub fn register_python_api_for_docs(m: &Bound<PyModule>) -> PyResult<()> {
    register_python_api(m)
}

/// Gather the Python API registered for `gammaloop._gammaloop`.
#[cfg(feature = "python_stubgen")]
pub fn stub_info() -> pyo3_stub_gen::Result<pyo3_stub_gen::StubInfo> {
    pyo3_stub_gen::StubInfo::from_project_root(
        "gammaloop._gammaloop".to_owned(),
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")),
    )
}

/// Read-only, detached view over serialized GammaLoop settings.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "SettingsValue", skip_from_py_object)]
#[derive(Clone)]
pub struct PySettingsValue {
    value: JsonValue,
    path: String,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PySettingsValue {
    /// Dot-separated location of this value within the settings tree.
    #[getter]
    fn path(&self) -> String {
        self.path.clone()
    }

    /// JSON value kind: ``object``, ``array``, ``string``, ``number``, ``boolean``, or ``null``.
    #[getter]
    fn kind(&self) -> &'static str {
        json_type_name(&self.value)
    }

    /// Return whether this value is an object with named children.
    fn is_object(&self) -> bool {
        self.value.is_object()
    }

    /// Return whether this value is an array with indexed children.
    fn is_array(&self) -> bool {
        self.value.is_array()
    }

    /// Return whether this value is a scalar JSON value, including ``null``.
    fn is_scalar(&self) -> bool {
        !self.value.is_object() && !self.value.is_array()
    }

    /// Return this object's child keys in lexical order.
    ///
    /// Raises
    /// ------
    /// TypeError
    ///     If this value is not an object.
    fn keys(&self) -> PyResult<Vec<String>> {
        match &self.value {
            JsonValue::Object(map) => {
                let mut keys = map.keys().cloned().collect::<Vec<_>>();
                keys.sort();
                Ok(keys)
            }
            _ => Err(exceptions::PyTypeError::new_err(format!(
                "Cannot list keys for {} at '{}'",
                json_type_name(&self.value),
                self.path
            ))),
        }
    }

    /// Resolve a dot-separated descendant path.
    ///
    /// An empty path returns this value. Objects become ``SettingsValue`` views,
    /// while scalar values become their ordinary Python equivalents.
    ///
    /// Parameters
    /// ----------
    /// path : str
    ///     Relative path such as ``"sampling.parameterization"``.
    ///
    /// Raises
    /// ------
    /// KeyError
    ///     If the path does not exist below this value.
    fn get<'py>(&self, py: Python<'py>, path: &str) -> PyResult<Bound<'py, PyAny>> {
        let value = if path.is_empty() {
            &self.value
        } else {
            value_at_path(&self.value, path).map_err(|e| {
                exceptions::PyKeyError::new_err(format!(
                    "Could not access '{}': {}",
                    self.extend_path(path),
                    e
                ))
            })?
        };
        py_object_from_settings_value(py, value, &self.extend_path(path))
    }

    /// Convert the complete detached value tree to ordinary Python containers and scalars.
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        py_builtin_from_settings_value(py, &self.value)
    }

    #[gen_stub(skip)]
    fn __getitem__<'py>(
        &self,
        py: Python<'py>,
        key: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Ok(path) = key.extract::<String>() {
            return self.get(py, &path);
        }
        if let Ok(index) = key.extract::<usize>() {
            let items = self.value.as_array().ok_or_else(|| {
                exceptions::PyTypeError::new_err(format!(
                    "Cannot use an integer index on {} at '{}'",
                    json_type_name(&self.value),
                    self.path
                ))
            })?;
            let value = items.get(index).ok_or_else(|| {
                exceptions::PyIndexError::new_err(format!(
                    "Index {} out of bounds for '{}' with len={}",
                    index,
                    self.path,
                    items.len()
                ))
            })?;
            return py_object_from_settings_value(py, value, &child_index_path(&self.path, index));
        }
        Err(exceptions::PyTypeError::new_err(
            "Settings values can only be indexed with a string path or integer index",
        ))
    }

    /// Resolve an object child by attribute name.
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Immediate child key to resolve.
    fn __getattr__<'py>(&self, py: Python<'py>, name: &str) -> PyResult<Bound<'py, PyAny>> {
        match &self.value {
            JsonValue::Object(map) => map.get(name).map_or_else(
                || {
                    Err(exceptions::PyAttributeError::new_err(format!(
                        "No attribute '{}' under '{}'",
                        name, self.path
                    )))
                },
                |value| {
                    py_object_from_settings_value(py, value, &child_field_path(&self.path, name))
                },
            ),
            _ => Err(exceptions::PyAttributeError::new_err(format!(
                "Cannot access attribute '{}' on {} at '{}'",
                name,
                json_type_name(&self.value),
                self.path
            ))),
        }
    }

    /// Return object-child names in lexical order, or an empty list for other kinds.
    fn __dir__(&self) -> Vec<String> {
        match &self.value {
            JsonValue::Object(map) => {
                let mut keys = map.keys().cloned().collect::<Vec<_>>();
                keys.sort();
                keys
            }
            _ => Vec::new(),
        }
    }

    /// Return the number of object or array children.
    fn __len__(&self) -> PyResult<usize> {
        match &self.value {
            JsonValue::Object(map) => Ok(map.len()),
            JsonValue::Array(items) => Ok(items.len()),
            _ => Err(exceptions::PyTypeError::new_err(format!(
                "Cannot take len() of {} at '{}'",
                json_type_name(&self.value),
                self.path
            ))),
        }
    }

    /// Return the detached settings value as pretty-printed JSON.
    fn __str__(&self) -> PyResult<String> {
        serde_json::to_string_pretty(&self.value).map_err(|e| {
            exceptions::PyException::new_err(format!("Could not format settings: {}", e))
        })
    }

    /// Return a compact representation containing this value's path and kind.
    fn __repr__(&self) -> String {
        format!(
            "SettingsValue(path='{}', kind='{}')",
            self.path,
            json_type_name(&self.value)
        )
    }
}

#[cfg(feature = "python_stubgen")]
pyo3_stub_gen::inventory::submit! {
    pyo3_stub_gen::derive::gen_methods_from_python! {
        r#"
        class PySettingsValue:
            def __getitem__(self, key: typing.Any, /) -> typing.Any:
                """Resolve a relative string path or array index.

                Parameters
                ----------
                key : str or int
                    Descendant path for an object, or zero-based index for an array.
                """
                ...
        "#
    }
}

impl PySettingsValue {
    fn from_settings<T: Serialize>(settings: &T, context: &str, path: &str) -> PyResult<Self> {
        let value = serialize_settings_with_defaults(settings, context).map_err(|e| {
            exceptions::PyException::new_err(format!("Could not serialize {}: {}", context, e))
        })?;
        Ok(Self {
            value,
            path: path.to_string(),
        })
    }

    fn extend_path(&self, relative_path: &str) -> String {
        if relative_path.is_empty() {
            self.path.clone()
        } else {
            format!("{}.{}", self.path, relative_path)
        }
    }
}

fn child_field_path(parent: &str, field: &str) -> String {
    format!("{}.{}", parent, field)
}

fn child_index_path(parent: &str, index: usize) -> String {
    format!("{}[{}]", parent, index)
}

fn py_object_from_settings_value<'py>(
    py: Python<'py>,
    value: &JsonValue,
    path: &str,
) -> PyResult<Bound<'py, PyAny>> {
    match value {
        JsonValue::Object(_) | JsonValue::Array(_) => {
            let value = Py::new(
                py,
                PySettingsValue {
                    value: value.clone(),
                    path: path.to_string(),
                },
            )?;
            Ok(value.into_bound(py).into_any())
        }
        _ => py_builtin_from_settings_value(py, value),
    }
}

fn py_builtin_from_settings_value<'py>(
    py: Python<'py>,
    value: &JsonValue,
) -> PyResult<Bound<'py, PyAny>> {
    match value {
        JsonValue::Null => Ok(py.None().into_bound(py)),
        JsonValue::Bool(value) => Ok(pyo3::types::PyBool::new(py, *value).to_owned().into_any()),
        JsonValue::Number(number) => py_builtin_from_json_number(py, number),
        JsonValue::String(value) => Ok(value.clone().into_pyobject(py)?.into_any()),
        JsonValue::Array(items) => {
            let list = PyList::empty(py);
            for item in items {
                list.append(py_builtin_from_settings_value(py, item)?)?;
            }
            Ok(list.into_any())
        }
        JsonValue::Object(map) => {
            let dict = PyDict::new(py);
            for (key, value) in map {
                dict.set_item(key, py_builtin_from_settings_value(py, value)?)?;
            }
            Ok(dict.into_any())
        }
    }
}

fn py_builtin_from_json_number<'py>(
    py: Python<'py>,
    number: &serde_json::Number,
) -> PyResult<Bound<'py, PyAny>> {
    if let Some(value) = number.as_i64() {
        return Ok(value.into_pyobject(py)?.into_any());
    }
    if let Some(value) = number.as_u64() {
        return Ok(value.into_pyobject(py)?.into_any());
    }
    if let Some(value) = number.as_f64() {
        return Ok(value.into_pyobject(py)?.into_any());
    }
    Err(exceptions::PyException::new_err(format!(
        "Unsupported numeric settings value: {}",
        number
    )))
}

#[cfg(test)]
mod settings_wrapper_tests {
    use super::*;
    use pyo3::types::PyAnyMethods;

    #[test]
    fn runtime_settings_wrapper_exposes_nested_values() {
        Python::initialize();

        let mut settings = RuntimeSettings::default();
        settings.general.generate_events = true;
        settings.subtraction.disable_threshold_subtraction = true;

        let wrapped =
            PySettingsValue::from_settings(&settings, "runtime settings", "runtime_settings")
                .unwrap();
        assert!(wrapped.keys().unwrap().contains(&"subtraction".to_string()));

        Python::attach(|py| {
            let subtraction = wrapped.get(py, "subtraction").unwrap();
            assert!(subtraction
                .getattr("disable_threshold_subtraction")
                .unwrap()
                .extract::<bool>()
                .unwrap());

            let settings_dict = wrapped.to_dict(py).unwrap();
            let general = settings_dict.get_item("general").unwrap();
            assert!(general
                .get_item("generate_events")
                .unwrap()
                .extract::<bool>()
                .unwrap());
        });
    }
}

#[pyclass(from_py_object, name = "ComplexValue", get_all)]
#[derive(Clone)]
pub struct PyComplexValue {
    pub re: f64,
    pub im: f64,
}

// `ComplexValue` is returned inside other records but is intentionally not
// registered as a module-level constructor. Keep its generated annotations
// honest by treating that opaque value as `Any` without adding a fake export.
#[cfg(feature = "python_stubgen")]
impl pyo3_stub_gen::PyStubType for PyComplexValue {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        pyo3_stub_gen::TypeInfo::any()
    }
}

#[cfg(feature = "python_stubgen")]
impl pyo3_stub_gen::PyStubType for ProcessRef {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        <i64 as pyo3_stub_gen::PyStubType>::type_output()
            | <String as pyo3_stub_gen::PyStubType>::type_output()
    }

    fn type_input() -> pyo3_stub_gen::TypeInfo {
        <i64 as pyo3_stub_gen::PyStubType>::type_input()
            | <String as pyo3_stub_gen::PyStubType>::type_input()
    }
}

/// One energy-momentum four-vector returned with a generated event.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "FourMomentum", get_all)]
#[derive(Clone)]
pub struct PyFourMomentum {
    /// Energy component.
    pub e: f64,
    /// Momentum along the x axis.
    pub px: f64,
    /// Momentum along the y axis.
    pub py: f64,
    /// Momentum along the z axis.
    pub pz: f64,
}

/// Named auxiliary complex weight attached to a generated event.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "AdditionalWeight", get_all)]
#[derive(Clone)]
pub struct PyAdditionalWeight {
    /// Stable weight name, optionally followed by colon-separated identifiers.
    pub key: String,
    /// Complex auxiliary event weight associated with ``key``.
    pub value: PyComplexValue,
}

/// Process, graph, orientation, and loop-basis identifiers for one cut event.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "CutInfo", get_all)]
#[derive(Clone)]
pub struct PyCutInfo {
    /// PDG identifiers aligned with ``Event.incoming_momenta``.
    pub incoming_pdgs: Vec<isize>,
    /// PDG identifiers aligned with ``Event.outgoing_momenta``.
    pub outgoing_pdgs: Vec<isize>,
    /// Zero-based cut identifier within the selected graph.
    pub cut_id: usize,
    /// Zero-based graph identifier within the integrand.
    pub graph_id: usize,
    /// Graph-group identifier, when group metadata is available.
    pub graph_group_id: Option<usize>,
    /// Causal-flow orientation identifier, when sampled explicitly.
    pub orientation_id: Option<usize>,
    /// Loop-momentum-basis multichannel identifier, when sampled explicitly.
    pub lmb_channel_id: Option<usize>,
    /// Edge identifiers defining the selected loop-momentum basis, when available.
    pub lmb_channel_edge_ids: Option<Vec<usize>>,
}

/// Identity and master-graph status of one graph in an integrand.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegrandGraph", get_all)]
#[derive(Clone)]
pub struct PyIntegrandGraphInfo {
    /// Zero-based graph identifier within the integrand.
    pub graph_id: usize,
    /// Persistent graph name.
    pub name: String,
    /// Whether this graph is the representative graph of its group.
    pub is_master: bool,
}

/// Edge-direction signature for one causal-flow orientation.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegrandOrientation", get_all)]
#[derive(Clone)]
pub struct PyIntegrandOrientationInfo {
    /// Zero-based orientation identifier within the graph group.
    pub orientation_id: usize,
    /// Direction per ``IntegrandGraphGroup.orientation_edge_ids``: ``1``, ``-1``, or ``0``.
    pub signature: Vec<i8>,
}

/// Loop-momentum basis and optional multichannel identifier for a graph group.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegrandLoopMomentumBasis", get_all)]
#[derive(Clone)]
pub struct PyIntegrandLoopMomentumBasisInfo {
    /// Zero-based loop-momentum-basis identifier on the representative graph.
    pub basis_id: usize,
    /// Effective multichannel identifier, or ``None`` when this basis is not sampled.
    pub channel_id: Option<usize>,
    /// Edge identifiers whose momenta form the basis.
    pub edge_ids: Vec<usize>,
    /// Whether this is the loop-momentum basis used during integrand generation.
    pub matches_generation_basis: bool,
}

/// Cut edges, raising power, and threshold associations for one cross-section cut.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegrandCut", get_all)]
#[derive(Clone)]
pub struct PyIntegrandCutInfo {
    /// Zero-based cut identifier within the representative graph.
    pub cut_id: usize,
    /// Edge identifiers crossed by the cut.
    pub edge_ids: Vec<usize>,
    /// Maximum multiplicity of the cut's related threshold group.
    pub raising_power: usize,
    /// Threshold associations on the left side of the cut.
    pub left_thresholds: Vec<PyIntegrandCutThresholdInfo>,
    /// Threshold associations on the right side of the cut.
    pub right_thresholds: Vec<PyIntegrandCutThresholdInfo>,
}

/// Classification and boundary edges of one threshold associated with a cut.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegrandCutThreshold", get_all)]
#[derive(Clone)]
pub struct PyIntegrandCutThresholdInfo {
    /// Threshold E-surface identifier within the graph group.
    pub esurface_id: usize,
    /// Viability classification for this threshold-cut association.
    pub status: String,
    /// Cut edges that bound the associated threshold region.
    pub cut_boundary_edge_ids: Vec<usize>,
    /// Threshold edges that bound the associated threshold region.
    pub threshold_boundary_edge_ids: Vec<usize>,
    /// Whether the invariant-bound criterion applies to this association.
    pub invariant_bound_is_applicable: bool,
}

/// Viable cut associated with a threshold E-surface.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegrandActiveThresholdCut", get_all)]
#[derive(Clone)]
pub struct PyIntegrandActiveThresholdCutInfo {
    /// Zero-based identifier of a cut on which the E-surface remains viable.
    pub cut_id: usize,
    /// Whether this cut can make the E-surface pinched for some kinematics.
    pub can_become_pinched: bool,
}

/// Threshold E-surface metadata and the cuts on which it remains active.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegrandThresholdEsurface", get_all)]
#[derive(Clone)]
pub struct PyIntegrandThresholdEsurfaceInfo {
    /// E-surface identifier within the graph group.
    pub esurface_id: usize,
    /// Graph whose local E-surface supplies ``edge_ids``.
    pub representative_graph_id: usize,
    /// Edge identifiers contributing energies to this E-surface.
    pub edge_ids: Vec<usize>,
    /// Kinematic existence class, or ``None`` when no static classification was computed.
    pub classification: Option<String>,
    /// Cuts on which this E-surface remains viable.
    pub active_cuts: Vec<PyIntegrandActiveThresholdCutInfo>,
}

/// Graphs and shared orientation, loop-basis, threshold, and cut data in one group.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegrandGraphGroup", get_all)]
#[derive(Clone)]
pub struct PyIntegrandGraphGroupInfo {
    /// Zero-based graph-group identifier within the integrand.
    pub group_id: usize,
    /// Graphs in this group, with the representative graph marked as master.
    pub graphs: Vec<PyIntegrandGraphInfo>,
    /// Edge identifiers that establish the ordering of every orientation signature.
    pub orientation_edge_ids: Vec<usize>,
    /// Available causal-flow orientations for the representative graph.
    pub orientations: Vec<PyIntegrandOrientationInfo>,
    /// Available loop-momentum bases for the representative graph.
    pub loop_momentum_bases: Vec<PyIntegrandLoopMomentumBasisInfo>,
    /// Ordered identifiers of threshold E-surfaces retained for the group.
    pub threshold_esurface_ids: Vec<usize>,
    /// Structured metadata for the retained threshold E-surfaces.
    pub threshold_esurfaces: Vec<PyIntegrandThresholdEsurfaceInfo>,
    /// Cross-section cuts; empty for amplitude graph groups.
    pub cuts: Vec<PyIntegrandCutInfo>,
}

/// Structured description of a generated integrand and its evaluation backend.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegrandInfo", get_all)]
#[derive(Clone)]
pub struct PyIntegrandInfo {
    /// Numeric identifier of the process that owns the integrand.
    pub process_id: usize,
    /// Persistent folder name of the owning process.
    pub process_name: String,
    /// Canonical name of the resolved integrand.
    pub integrand_name: String,
    /// ``"amplitude"`` or ``"cross section"``.
    pub kind: String,
    /// Compilation backend frozen into the generated integrand.
    pub generation_backend: String,
    /// Backend-specific compilation options, when configured.
    pub generation_compile_options: Option<String>,
    /// Floating-point evaluator backend currently active for f64 evaluation.
    pub active_f64_backend: String,
    /// Number of amplitude graphs or cross-section supergraphs.
    pub graph_count: usize,
    /// Number of graph groups in the generated integrand.
    pub graph_group_count: usize,
    /// Serialized integrand-record size in bytes, excluding its compiled evaluator.
    pub record_size_bytes: usize,
    /// Per-group graph, orientation, basis, threshold, and cut metadata.
    pub graph_groups: Vec<PyIntegrandGraphGroupInfo>,
}

/// Accepted cut event with momenta, its primary weight, and auxiliary weights.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "Event", get_all)]
#[derive(Clone)]
pub struct PyEvent {
    /// Incoming four-momenta, aligned with ``cut_info.incoming_pdgs``.
    pub incoming_momenta: Vec<PyFourMomentum>,
    /// Outgoing four-momenta, aligned with ``cut_info.outgoing_pdgs``.
    pub outgoing_momenta: Vec<PyFourMomentum>,
    /// Cut, graph, orientation, and loop-basis identifiers for this event.
    pub cut_info: PyCutInfo,
    /// Primary complex event weight.
    pub weight: PyComplexValue,
    /// Named auxiliary complex weights such as threshold-counterterm contributions.
    pub additional_weights: Vec<PyAdditionalWeight>,
}

/// Correlated accepted events produced by one graph group for one sample.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "EventGroup", get_all)]
#[derive(Clone)]
pub struct PyEventGroup {
    /// Correlated accepted events in this group.
    pub events: Vec<PyEvent>,
}

/// Mutable continuous or discrete histogram accumulator with sample-level statistics.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "HistogramAccumulator", skip_from_py_object)]
#[derive(Clone)]
pub struct PyHistogramAccumulator {
    inner: HistogramAccumulatorState,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyHistogramAccumulator {
    /// Create an evenly binned continuous histogram accumulator.
    ///
    /// Parameters
    /// ----------
    /// title : str
    ///     Human-readable title used by snapshots and exported output.
    /// x_min : float
    ///     Lower edge of the histogram range.
    /// x_max : float
    ///     Upper edge of the histogram range.
    /// n_bins : int
    ///     Number of equal-width in-range bins.
    /// type_description : str, default="AL"
    ///     HwU ``TYPE@`` metadata.
    /// phase : {"real", "imag"}, default="real"
    ///     Component of each complex event weight to accumulate.
    /// value_transform : {"identity", "log10"}, default="identity"
    ///     Transformation applied to coordinates before binning.
    /// log_x_axis : bool, default=False
    ///     Request logarithmic horizontal-axis rendering.
    /// log_y_axis : bool, default=True
    ///     Request logarithmic vertical-axis rendering.
    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (title, x_min, x_max, n_bins, type_description="AL".to_string(), phase="real".to_string(), value_transform="identity".to_string(), log_x_axis=false, log_y_axis=true))]
    fn continuous(
        title: String,
        x_min: f64,
        x_max: f64,
        n_bins: usize,
        type_description: String,
        phase: String,
        value_transform: String,
        log_x_axis: bool,
        log_y_axis: bool,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: HistogramAccumulatorState::continuous(
                title,
                type_description,
                py_phase_to_phase(&phase)?,
                py_value_transform_to_transform(&value_transform)?,
                x_min,
                x_max,
                log_x_axis,
                log_y_axis,
                n_bins,
            ),
        })
    }

    /// Create a discrete histogram accumulator over an inclusive integer range.
    ///
    /// Parameters
    /// ----------
    /// title : str
    ///     Human-readable title used by snapshots and exported output.
    /// min_bin_id : int
    ///     Smallest accepted bin identifier.
    /// max_bin_id : int
    ///     Largest accepted bin identifier, inclusive.
    /// ordering : {"ascending_bin_id", "value_descending", "abs_value_descending"}, default="ascending_bin_id"
    ///     Ordering used when bins are returned or exported.
    /// labels : Sequence[str], optional
    ///     One label per bin, ordered by ascending bin identifier.
    /// type_description : str, default="AL"
    ///     HwU ``TYPE@`` metadata.
    /// phase : {"real", "imag"}, default="real"
    ///     Component of each complex event weight to accumulate.
    /// log_y_axis : bool, default=True
    ///     Request logarithmic vertical-axis rendering.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If the range, label count, ordering, or phase is invalid.
    #[staticmethod]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (title, min_bin_id, max_bin_id, ordering="ascending_bin_id".to_string(), labels=None, type_description="AL".to_string(), phase="real".to_string(), log_y_axis=true))]
    fn discrete(
        title: String,
        min_bin_id: isize,
        max_bin_id: isize,
        ordering: String,
        labels: Option<Vec<String>>,
        type_description: String,
        phase: String,
        log_y_axis: bool,
    ) -> PyResult<Self> {
        if max_bin_id < min_bin_id {
            return Err(exceptions::PyValueError::new_err(format!(
                "invalid discrete bin range: max_bin_id ({max_bin_id}) must be >= min_bin_id ({min_bin_id})"
            )));
        }
        let n_bins = max_bin_id
            .checked_sub(min_bin_id)
            .and_then(|delta| delta.checked_add(1))
            .ok_or_else(|| exceptions::PyValueError::new_err("invalid discrete bin range"))?
            as usize;
        let labels = match labels {
            Some(labels) if labels.len() != n_bins => {
                return Err(exceptions::PyValueError::new_err(format!(
                    "expected {} labels for discrete range [{}, {}], got {}",
                    n_bins,
                    min_bin_id,
                    max_bin_id,
                    labels.len()
                )));
            }
            Some(labels) => labels.into_iter().map(Some).collect(),
            None => vec![None; n_bins],
        };
        Ok(Self {
            inner: HistogramAccumulatorState::discrete(
                title,
                type_description,
                py_phase_to_phase(&phase)?,
                min_bin_id,
                max_bin_id,
                py_discrete_ordering(&ordering)?,
                log_y_axis,
                labels,
            )
            .map_err(to_py_value_error)?,
        })
    }

    /// Return an immutable snapshot including committed and pending samples.
    fn snapshot(&self) -> PyHistogramSnapshot {
        py_histogram_snapshot_from_snapshot(self.inner.snapshot())
    }

    /// Move pending samples from another compatible accumulator into this one.
    ///
    /// Parameters
    /// ----------
    /// other : HistogramAccumulator
    ///     Compatible accumulator whose pending samples are consumed.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If histogram definitions differ.
    fn merge_in_place(
        &mut self,
        #[gen_stub(override_type(type_repr = "HistogramAccumulator"))]
        other: &mut PyHistogramAccumulator,
    ) -> PyResult<()> {
        self.inner
            .merge_in_place(&mut other.inner)
            .map_err(to_py_value_error)
    }

    /// Return a copy with each run of adjacent continuous bins combined.
    ///
    /// Discrete histograms are copied unchanged.
    ///
    /// Parameters
    /// ----------
    /// contiguous_bins : int
    ///     Positive number of old bins per new bin; it must divide the bin count.
    fn rebin(&self, contiguous_bins: usize) -> PyResult<Self> {
        self.inner
            .rebin(contiguous_bins)
            .map(|inner| Self { inner })
            .map_err(to_py_value_error)
    }

    /// Multiply all accumulated weights by a constant in place.
    ///
    /// Parameters
    /// ----------
    /// factor : float
    ///     Scale applied to weight sums; squared-weight sums use its square.
    fn rescale(&mut self, factor: f64) {
        self.inner.rescale(factor);
    }

    /// Change the display ordering of a discrete histogram.
    ///
    /// Parameters
    /// ----------
    /// ordering : {"ascending_bin_id", "value_descending", "abs_value_descending"}
    ///     New ordering for snapshots and exported output.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If this is a continuous histogram or the name is invalid.
    fn change_bin_ordering(&mut self, ordering: String) -> PyResult<()> {
        self.inner
            .change_bin_ordering(py_discrete_ordering(&ordering)?)
            .map_err(to_py_value_error)
    }

    /// Commit pending samples to the accumulator's completed-result counters.
    fn update_results(&mut self) {
        self.inner.update_results();
    }

    /// Add one independent continuous sample containing coordinate-weight pairs.
    ///
    /// Parameters
    /// ----------
    /// entries : Sequence[tuple[float, float]]
    ///     ``(coordinate, projected_weight)`` entries for this sample.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If this is not a continuous histogram.
    fn fill_continuous_sample(&mut self, entries: Vec<(f64, f64)>) -> PyResult<()> {
        self.inner
            .fill_continuous_sample(&entries)
            .map_err(to_py_value_error)
    }

    /// Add one independent discrete sample containing bin-weight pairs.
    ///
    /// Parameters
    /// ----------
    /// entries : Sequence[tuple[int, float]]
    ///     ``(bin_id, projected_weight)`` entries for this sample.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     If this is not a discrete histogram.
    fn fill_discrete_sample(&mut self, entries: Vec<(isize, f64)>) -> PyResult<()> {
        self.inner
            .fill_discrete_sample(&entries)
            .map_err(to_py_value_error)
    }
}

/// Raw statistics and optional coordinate metadata for one histogram bin.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "HistogramBin", get_all)]
#[derive(Clone)]
pub struct PyHistogramBinSnapshot {
    /// Lower coordinate edge, or ``None`` for a discrete or overflow bin.
    pub x_min: Option<f64>,
    /// Upper coordinate edge, or ``None`` for a discrete or underflow bin.
    pub x_max: Option<f64>,
    /// Discrete bin identifier, or ``None`` for a continuous bin.
    pub bin_id: Option<isize>,
    /// Optional label of a discrete bin.
    pub label: Option<String>,
    /// Number of event entries accumulated in this bin.
    pub entry_count: usize,
    /// Sum of per-sample projected weights for this bin.
    pub sum_weights: f64,
    /// Sum of squared per-sample projected weights for this bin.
    pub sum_weights_squared: f64,
    /// Number of fills produced by paired boundary-misbinning mitigation.
    pub mitigated_fill_count: usize,
}

/// Aggregate entry, NaN, and misbinning-mitigation counts for a histogram.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "HistogramStats", get_all)]
#[derive(Clone)]
pub struct PyHistogramStatisticsSnapshot {
    /// Number of event entries assigned to in-range bins.
    pub in_range_entry_count: usize,
    /// Number of entries skipped because their coordinate was non-finite.
    pub nan_value_count: usize,
    /// Number of entry pairs processed by boundary-misbinning mitigation.
    pub mitigated_pair_count: usize,
}

/// Immutable, mergeable histogram snapshot containing raw per-bin statistics.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "HistogramSnapshot", get_all)]
#[derive(Clone)]
pub struct PyHistogramSnapshot {
    /// ``"continuous"`` or ``"discrete"``.
    pub kind: String,
    /// Human-readable histogram title.
    pub title: String,
    /// HwU ``TYPE@`` description retained with the histogram.
    pub type_description: String,
    /// Complex-weight component accumulated: ``"real"`` or ``"imag"``.
    pub phase: String,
    /// Coordinate transform: ``"identity"`` or ``"log10"``.
    pub value_transform: String,
    /// Whether paired fills can mitigate bin-boundary instability.
    pub supports_misbinning_mitigation: bool,
    /// Lower bound of a continuous histogram, otherwise ``None``.
    pub x_min: Option<f64>,
    /// Upper bound of a continuous histogram, otherwise ``None``.
    pub x_max: Option<f64>,
    /// Number of statistically independent Monte Carlo samples represented.
    pub sample_count: usize,
    /// Whether compatible renderers should use a logarithmic horizontal axis.
    pub log_x_axis: bool,
    /// Whether compatible renderers should use a logarithmic vertical axis.
    pub log_y_axis: bool,
    /// Smallest discrete bin identifier, otherwise ``None``.
    pub discrete_min_bin_id: Option<isize>,
    /// Discrete display ordering, otherwise ``None``.
    pub discrete_ordering: Option<String>,
    /// Raw statistics for every in-range bin, in display order.
    pub bins: Vec<PyHistogramBinSnapshot>,
    /// Raw statistics for entries below a continuous histogram's range.
    pub underflow_bin: PyHistogramBinSnapshot,
    /// Raw statistics for entries above a continuous histogram's range.
    pub overflow_bin: PyHistogramBinSnapshot,
    /// Histogram-wide entry and mitigation counters.
    pub statistics: PyHistogramStatisticsSnapshot,
}

/// Completed integration result, observable snapshots, and workspace location.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(name = "IntegrationOutput", skip_from_py_object)]
#[derive(Clone)]
pub struct PyIntegrationOutput {
    inner: crate::commands::integrate::IntegrationOutput,
}

/// Complex integral estimate with errors, chi-squared values, and evaluation counts.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegralEstimate", get_all)]
#[derive(Clone)]
pub struct PyIntegralEstimate {
    pub neval: usize,
    pub real_zero: usize,
    pub im_zero: usize,
    pub result: PyComplexValue,
    pub error: PyComplexValue,
    pub real_chisq: f64,
    pub im_chisq: f64,
}

/// One real or imaginary component row from an integration result table.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegrationTableComponentResult", get_all)]
#[derive(Clone)]
pub struct PyIntegrationTableComponentResult {
    pub component: String,
    pub value: f64,
    pub error: f64,
    pub relative_error_percent: Option<f64>,
    pub chi_sq_per_dof: f64,
    pub target_delta_sigma: Option<f64>,
    pub target_delta_percent: Option<f64>,
    pub max_weight_impact: f64,
}

/// Timing, precision, stability, and event-selection statistics for an integration slot.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegrationStatisticsSnapshot", get_all)]
#[derive(Clone)]
pub struct PyIntegrationStatisticsSnapshot {
    pub num_evals: usize,
    pub average_total_time_seconds: f64,
    pub average_parameterization_time_seconds: f64,
    pub average_integrand_time_seconds: f64,
    pub average_evaluator_time_seconds: f64,
    pub average_observable_time_seconds: f64,
    pub average_integrator_time_seconds: f64,
    pub f64_percentage: f64,
    pub f128_percentage: f64,
    pub arb_percentage: f64,
    pub nan_percentage: f64,
    pub nan_or_unstable_percentage: f64,
    pub generated_event_count: usize,
    pub accepted_event_count: usize,
    pub selection_efficiency_percentage: Option<f64>,
}

/// Largest observed contribution for one signed integration component.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "MaxWeightInfoEntry", get_all)]
#[derive(Clone)]
pub struct PyMaxWeightInfoEntry {
    pub component: String,
    pub sign: String,
    pub max_eval: f64,
    pub coordinates: Option<String>,
}

/// One labeled coordinate on a discrete integration axis.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "DiscreteCoordinate", get_all)]
#[derive(Clone)]
pub struct PyDiscreteCoordinate {
    pub axis_label: String,
    pub bin_index: usize,
    pub bin_label: Option<String>,
}

/// Estimate and sampling metadata for one bin of a discrete-axis breakdown.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "DiscreteBreakdownEntry", get_all)]
#[derive(Clone)]
pub struct PyDiscreteBreakdownEntry {
    pub bin_index: usize,
    pub bin_label: Option<String>,
    pub pdf: f64,
    pub value: f64,
    pub error: f64,
    pub chi_sq: f64,
    pub processed_samples: usize,
}

/// Per-bin integration estimates along one discrete axis at fixed other coordinates.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "DiscreteBreakdown", get_all)]
#[derive(Clone)]
pub struct PyDiscreteBreakdown {
    pub axis_label: String,
    pub fixed_coordinates: Vec<PyDiscreteCoordinate>,
    pub entries: Vec<PyDiscreteBreakdownEntry>,
}

/// Optional real and imaginary discrete-axis breakdowns for an integration slot.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "ComponentDiscreteBreakdown", get_all)]
#[derive(Clone)]
pub struct PyComponentDiscreteBreakdown {
    pub re: Option<PyDiscreteBreakdown>,
    pub im: Option<PyDiscreteBreakdown>,
}

/// Result, diagnostics, and discrete breakdown for one named integration slot.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "SlotIntegrationResult", get_all)]
#[derive(Clone)]
pub struct PySlotIntegrationResult {
    pub key: String,
    pub process: String,
    pub integrand: String,
    pub target: Option<PyComplexValue>,
    pub integral: PyIntegralEstimate,
    pub table_results: Vec<PyIntegrationTableComponentResult>,
    pub integration_statistics: PyIntegrationStatisticsSnapshot,
    pub max_weight_info: Vec<PyMaxWeightInfoEntry>,
    pub grid_breakdown: PyComponentDiscreteBreakdown,
}

/// Collection of independently addressable integration-slot results.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "IntegrationResult", get_all)]
#[derive(Clone)]
pub struct PyIntegrationResult {
    pub slots: Vec<PySlotIntegrationResult>,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyHistogramBinSnapshot {
    /// Compute this bin's Monte Carlo mean from its raw weight sum.
    ///
    /// Parameters
    /// ----------
    /// sample_count : int
    ///     Number of statistically independent samples represented by the histogram.
    ///
    /// Returns
    /// -------
    /// float
    ///     ``sum_weights / sample_count``, or zero when ``sample_count`` is zero.
    fn average(&self, sample_count: usize) -> f64 {
        if sample_count == 0 {
            0.0
        } else {
            self.sum_weights / sample_count as f64
        }
    }

    /// Compute the standard error of this bin's Monte Carlo mean.
    ///
    /// Parameters
    /// ----------
    /// sample_count : int
    ///     Number of statistically independent samples represented by the histogram.
    ///
    /// Returns
    /// -------
    /// float
    ///     Sample-mean standard error, or zero when it cannot be estimated reliably.
    fn error(&self, sample_count: usize) -> f64 {
        if sample_count <= 1 {
            return 0.0;
        }
        let n = sample_count as f64;
        let variance_numerator =
            self.sum_weights_squared - (self.sum_weights * self.sum_weights) / n;
        if !variance_numerator.is_finite() || variance_numerator <= 0.0 {
            0.0
        } else {
            (variance_numerator / (n * (n - 1.0))).sqrt()
        }
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyIntegrationOutput {
    #[getter]
    fn workspace_path(&self) -> PathBuf {
        self.inner.workspace_path.clone()
    }

    #[getter]
    fn result(&self) -> PyIntegrationResult {
        py_integration_result_from_result(self.inner.result.clone())
    }

    #[getter]
    fn observables<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let dict = PyDict::new(py);
        for (slot_key, bundle) in &self.inner.observables {
            dict.set_item(slot_key, py_observable_dict_from_bundle(py, bundle)?)?;
        }
        Ok(dict)
    }

    fn __repr__(&self) -> String {
        format!(
            "IntegrationOutput(workspace_path='{}', slots={})",
            self.inner.workspace_path.display(),
            self.inner.result.slots.len()
        )
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyIntegrationResult {
    fn slot(&self, key: &str) -> Option<PySlotIntegrationResult> {
        self.slots.iter().find(|slot| slot.key == key).cloned()
    }

    fn single_slot(&self) -> Option<PySlotIntegrationResult> {
        (self.slots.len() == 1).then(|| self.slots[0].clone())
    }
}

/// Outcome and cost of one numerical-stability precision level.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "StabilityResult", get_all)]
#[derive(Clone)]
pub struct PyStabilityResult {
    /// Numerical precision used for this stability level.
    pub precision: String,
    /// Estimated relative accuracy, or ``None`` when it could not be estimated.
    pub estimated_relative_accuracy: Option<f64>,
    /// Human-readable stability outcome and number of rotated samples.
    pub status: String,
    /// Number of rotated samples evaluated at this precision level.
    pub sample_count: usize,
    /// Total wall-clock time spent at this precision level, in seconds.
    pub total_time_seconds: f64,
}

/// Numerical value, metadata, stability records, and generated events for one sample.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "SampleEvaluationResult")]
#[derive(Clone)]
pub struct PySampleEvaluationResult {
    inner: SampleEvaluationResult,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PySampleEvaluationResult {
    /// Complex integrand value before applying the parameterization Jacobian.
    #[getter]
    fn integrand_result<'py>(&self, py: Python<'py>) -> Bound<'py, PyComplex> {
        PyComplex::from_doubles(
            py,
            self.inner.evaluation.integrand_result.re.0,
            self.inner.evaluation.integrand_result.im.0,
        )
    }

    /// Monte Carlo weight supplied by the integrator, excluding the parameterization Jacobian.
    #[getter]
    fn integrator_weight(&self) -> f64 {
        self.inner.evaluation.integrator_weight.0
    }

    /// Number of candidate events generated, or ``None`` when metadata was omitted.
    #[getter]
    fn generated_event_count(&self) -> Option<usize> {
        self.inner
            .evaluation
            .evaluation_metadata
            .as_ref()
            .map(|metadata| metadata.generated_event_count)
    }

    /// Number of generated events accepted by selectors, or ``None`` without metadata.
    #[getter]
    fn accepted_event_count(&self) -> Option<usize> {
        self.inner
            .evaluation
            .evaluation_metadata
            .as_ref()
            .map(|metadata| metadata.accepted_event_count)
    }

    /// Time spent building and selecting events, in seconds, or ``None`` without metadata.
    #[getter]
    fn event_processing_time_seconds(&self) -> Option<f64> {
        self.inner
            .evaluation
            .evaluation_metadata
            .as_ref()
            .map(|metadata| metadata.event_processing_time.as_secs_f64())
    }

    /// Sample parameterization Jacobian, or ``None`` when no Jacobian was supplied.
    #[getter]
    fn parameterization_jacobian(&self) -> Option<f64> {
        self.inner
            .evaluation
            .parameterization_jacobian
            .map(|jac| jac.0)
    }

    /// Whether evaluation metadata flagged a non-finite result, or ``None`` without metadata.
    #[getter]
    fn is_nan(&self) -> Option<bool> {
        self.inner
            .evaluation
            .evaluation_metadata
            .as_ref()
            .map(|metadata| metadata.is_nan)
    }

    /// Per-precision stability attempts, or ``None`` when metadata was omitted.
    #[getter]
    fn stability_results(&self) -> Option<Vec<PyStabilityResult>> {
        self.inner
            .evaluation
            .evaluation_metadata
            .as_ref()
            .map(|metadata| {
                metadata
                    .stability_results
                    .iter()
                    .map(|result| PyStabilityResult {
                        precision: result.precision.to_string(),
                        estimated_relative_accuracy: result
                            .estimated_relative_accuracy
                            .map(|value| value.0),
                        status: result.status.to_string(),
                        sample_count: result.status.sample_count(),
                        total_time_seconds: result.total_time.as_secs_f64(),
                    })
                    .collect()
            })
    }

    /// Correlated event groups produced while evaluating this sample.
    #[getter]
    fn event_groups(&self) -> Vec<PyEventGroup> {
        self.inner
            .evaluation
            .event_groups
            .iter()
            .map(py_event_group_from_event_group)
            .collect()
    }

    /// Return a human-readable evaluation summary.
    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    /// Return the human-readable evaluation summary.
    fn __repr__(&self) -> String {
        self.__str__()
    }
}

/// Single-sample evaluation paired with the observable snapshot for that sample.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "EvaluationResult")]
#[derive(Clone)]
pub struct PyEvaluationResult {
    inner: SingleSampleEvaluationResult,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyEvaluationResult {
    /// Evaluation record for the requested sample.
    #[getter]
    fn sample(&self) -> PySampleEvaluationResult {
        PySampleEvaluationResult {
            inner: self.inner.sample.clone(),
        }
    }

    /// Complex integrand value before applying the parameterization Jacobian.
    #[getter]
    fn integrand_result<'py>(&self, py: Python<'py>) -> Bound<'py, PyComplex> {
        self.sample().integrand_result(py)
    }

    /// Monte Carlo weight supplied by the integrator, excluding the parameterization Jacobian.
    #[getter]
    fn integrator_weight(&self) -> f64 {
        self.inner.sample.evaluation.integrator_weight.0
    }

    /// Number of candidate events generated, or ``None`` when metadata was omitted.
    #[getter]
    fn generated_event_count(&self) -> Option<usize> {
        self.sample().generated_event_count()
    }

    /// Number of generated events accepted by selectors, or ``None`` without metadata.
    #[getter]
    fn accepted_event_count(&self) -> Option<usize> {
        self.sample().accepted_event_count()
    }

    /// Time spent building and selecting events, in seconds, or ``None`` without metadata.
    #[getter]
    fn event_processing_time_seconds(&self) -> Option<f64> {
        self.sample().event_processing_time_seconds()
    }

    /// Sample parameterization Jacobian, or ``None`` when no Jacobian was supplied.
    #[getter]
    fn parameterization_jacobian(&self) -> Option<f64> {
        self.inner
            .sample
            .evaluation
            .parameterization_jacobian
            .map(|jac| jac.0)
    }

    /// Whether evaluation metadata flagged a non-finite result, or ``None`` without metadata.
    #[getter]
    fn is_nan(&self) -> Option<bool> {
        self.sample().is_nan()
    }

    /// Per-precision stability attempts, or ``None`` when metadata was omitted.
    #[getter]
    fn stability_results(&self) -> Option<Vec<PyStabilityResult>> {
        self.sample().stability_results()
    }

    /// Correlated event groups produced while evaluating this sample.
    #[getter]
    fn event_groups(&self) -> Vec<PyEventGroup> {
        self.sample().event_groups()
    }

    /// Observable snapshots produced from the same Monte Carlo sample.
    #[getter]
    fn observables<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        py_observable_dict_from_bundle(py, &self.inner.observables)
    }

    /// Return a human-readable evaluation and observable summary.
    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    /// Return the human-readable evaluation and observable summary.
    fn __repr__(&self) -> String {
        self.__str__()
    }
}

/// Per-sample evaluations paired with one observable snapshot for the complete batch.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(from_py_object, name = "BatchEvaluationResult")]
#[derive(Clone)]
pub struct PyBatchEvaluationResult {
    inner: BatchSampleEvaluationResult,
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyBatchEvaluationResult {
    /// Evaluation records in the same order as the requested samples.
    #[getter]
    fn samples(&self) -> Vec<PySampleEvaluationResult> {
        self.inner
            .samples
            .iter()
            .cloned()
            .map(|inner| PySampleEvaluationResult { inner })
            .collect()
    }

    /// Observable snapshots accumulated across the complete sample batch.
    #[getter]
    fn observables<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        py_observable_dict_from_bundle(py, &self.inner.observables)
    }

    /// Return a human-readable batch and observable summary.
    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    /// Return the human-readable batch and observable summary.
    fn __repr__(&self) -> String {
        self.__str__()
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyEvent {
    /// Return a human-readable event summary.
    fn __str__(&self) -> String {
        event_from_py_event(self).to_string()
    }

    /// Return the human-readable event summary.
    fn __repr__(&self) -> String {
        self.__str__()
    }
}

#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyEventGroup {
    /// Return a human-readable summary of the grouped events.
    fn __str__(&self) -> String {
        event_group_from_py_event_group(self).to_string()
    }

    /// Return the human-readable grouped-event summary.
    fn __repr__(&self) -> String {
        self.__str__()
    }
}

fn py_complex_from_complex(
    value: spenso::algebra::complex::Complex<gammalooprs::utils::F<f64>>,
) -> PyComplexValue {
    PyComplexValue {
        re: value.re.0,
        im: value.im.0,
    }
}

fn py_four_momentum_from_momentum(
    momentum: &gammalooprs::momentum::FourMomentum<gammalooprs::utils::F<f64>>,
) -> PyFourMomentum {
    PyFourMomentum {
        e: momentum.temporal.value.0,
        px: momentum.spatial.px.0,
        py: momentum.spatial.py.0,
        pz: momentum.spatial.pz.0,
    }
}

fn additional_weight_key_to_string(key: AdditionalWeightKey) -> String {
    match key {
        AdditionalWeightKey::FullMultiplicativeFactor => "full_multiplicative_factor".to_string(),
        AdditionalWeightKey::Original => "original".to_string(),
        AdditionalWeightKey::ThresholdCounterterm { subset_index } => {
            format!("threshold_counterterm:{subset_index}")
        }
        AdditionalWeightKey::AmplitudeThresholdCounterterm {
            esurface_id,
            overlap_group,
        } => format!("threshold_counterterm:{esurface_id}:{overlap_group}"),
    }
}

fn py_event_from_event(event: &Event) -> PyEvent {
    let additional_weights = event
        .additional_weights
        .weights
        .iter()
        .map(|(key, value)| PyAdditionalWeight {
            key: additional_weight_key_to_string(*key),
            value: py_complex_from_complex(*value),
        })
        .collect();

    PyEvent {
        incoming_momenta: event
            .kinematic_configuration
            .0
            .iter()
            .map(py_four_momentum_from_momentum)
            .collect(),
        outgoing_momenta: event
            .kinematic_configuration
            .1
            .iter()
            .map(py_four_momentum_from_momentum)
            .collect(),
        cut_info: PyCutInfo {
            incoming_pdgs: event.cut_info.particle_pdgs.0.iter().copied().collect(),
            outgoing_pdgs: event.cut_info.particle_pdgs.1.iter().copied().collect(),
            cut_id: event.cut_info.cut_id,
            graph_id: event.cut_info.graph_id,
            graph_group_id: event.cut_info.graph_group_id,
            orientation_id: event.cut_info.orientation_id,
            lmb_channel_id: event.cut_info.lmb_channel_id,
            lmb_channel_edge_ids: event
                .cut_info
                .lmb_channel_edge_ids
                .as_ref()
                .map(|edge_ids| edge_ids.iter().copied().collect()),
        },
        weight: py_complex_from_complex(event.weight),
        additional_weights,
    }
}

fn event_from_py_event(event: &PyEvent) -> Event {
    let incoming_momenta = event
        .incoming_momenta
        .iter()
        .map(|momentum| gammalooprs::momentum::FourMomentum {
            temporal: gammalooprs::momentum::Energy {
                value: gammalooprs::utils::F(momentum.e),
            },
            spatial: gammalooprs::momentum::ThreeMomentum {
                px: gammalooprs::utils::F(momentum.px),
                py: gammalooprs::utils::F(momentum.py),
                pz: gammalooprs::utils::F(momentum.pz),
            },
        })
        .collect();
    let outgoing_momenta = event
        .outgoing_momenta
        .iter()
        .map(|momentum| gammalooprs::momentum::FourMomentum {
            temporal: gammalooprs::momentum::Energy {
                value: gammalooprs::utils::F(momentum.e),
            },
            spatial: gammalooprs::momentum::ThreeMomentum {
                px: gammalooprs::utils::F(momentum.px),
                py: gammalooprs::utils::F(momentum.py),
                pz: gammalooprs::utils::F(momentum.pz),
            },
        })
        .collect();
    let additional_weights = event
        .additional_weights
        .iter()
        .map(|weight| {
            let key = match weight.key.as_str() {
                "original" => AdditionalWeightKey::Original,
                "full_multiplicative_factor" => AdditionalWeightKey::FullMultiplicativeFactor,
                _ => match weight.key.strip_prefix("threshold_counterterm:") {
                    Some(indices) => {
                        let mut indices = indices.split(':');
                        let first = indices
                            .next()
                            .unwrap_or_default()
                            .parse()
                            .unwrap_or_default();
                        match indices.next() {
                            Some(second) => AdditionalWeightKey::AmplitudeThresholdCounterterm {
                                esurface_id: first,
                                overlap_group: second.parse().unwrap_or_default(),
                            },
                            None => AdditionalWeightKey::ThresholdCounterterm {
                                subset_index: first,
                            },
                        }
                    }
                    None => AdditionalWeightKey::Original,
                },
            };
            (
                key,
                spenso::algebra::complex::Complex::new(
                    gammalooprs::utils::F(weight.value.re),
                    gammalooprs::utils::F(weight.value.im),
                ),
            )
        })
        .collect();

    Event {
        kinematic_configuration: (incoming_momenta, outgoing_momenta),
        cut_info: gammalooprs::observables::CutInfo {
            particle_pdgs: (
                event.cut_info.incoming_pdgs.iter().copied().collect(),
                event.cut_info.outgoing_pdgs.iter().copied().collect(),
            ),
            cut_id: event.cut_info.cut_id,
            graph_id: event.cut_info.graph_id,
            graph_group_id: event.cut_info.graph_group_id,
            orientation_id: event.cut_info.orientation_id,
            lmb_channel_id: event.cut_info.lmb_channel_id,
            lmb_channel_edge_ids: event
                .cut_info
                .lmb_channel_edge_ids
                .as_ref()
                .map(|edge_ids| edge_ids.iter().copied().collect()),
        },
        weight: spenso::algebra::complex::Complex::new(
            gammalooprs::utils::F(event.weight.re),
            gammalooprs::utils::F(event.weight.im),
        ),
        additional_weights: GenericAdditionalWeightInfo {
            weights: additional_weights,
        },
        derived_observable_data: Default::default(),
    }
}

fn py_event_group_from_event_group(group: &EventGroup) -> PyEventGroup {
    PyEventGroup {
        events: group.iter().map(py_event_from_event).collect(),
    }
}

fn event_group_from_py_event_group(group: &PyEventGroup) -> EventGroup {
    gammalooprs::observables::GenericEventGroup(
        group.events.iter().map(event_from_py_event).collect(),
    )
}

fn py_observable_dict_from_bundle<'py>(
    py: Python<'py>,
    bundle: &gammalooprs::observables::ObservableSnapshotBundle,
) -> PyResult<Bound<'py, PyDict>> {
    let dict = PyDict::new(py);
    for (name, histogram) in &bundle.histograms {
        dict.set_item(name, py_histogram_snapshot_from_snapshot(histogram.clone()))?;
    }
    Ok(dict)
}

fn py_phase_to_phase(phase: &str) -> PyResult<gammalooprs::observables::ObservablePhase> {
    match phase {
        "real" => Ok(gammalooprs::observables::ObservablePhase::Real),
        "imag" => Ok(gammalooprs::observables::ObservablePhase::Imag),
        _ => Err(exceptions::PyValueError::new_err(format!(
            "unknown histogram phase '{}'",
            phase
        ))),
    }
}

fn py_value_transform_to_transform(
    value_transform: &str,
) -> PyResult<gammalooprs::observables::ObservableValueTransform> {
    match value_transform {
        "identity" => Ok(gammalooprs::observables::ObservableValueTransform::Identity),
        "log10" => Ok(gammalooprs::observables::ObservableValueTransform::Log10),
        _ => Err(exceptions::PyValueError::new_err(format!(
            "unknown histogram value_transform '{}'",
            value_transform
        ))),
    }
}

fn py_discrete_ordering(ordering: &str) -> PyResult<DiscreteBinOrdering> {
    match ordering {
        "ascending_bin_id" => Ok(DiscreteBinOrdering::AscendingBinId),
        "value_descending" => Ok(DiscreteBinOrdering::ValueDescending),
        "abs_value_descending" => Ok(DiscreteBinOrdering::AbsValueDescending),
        _ => Err(exceptions::PyValueError::new_err(format!(
            "unknown discrete bin ordering '{}'",
            ordering
        ))),
    }
}

fn to_py_value_error(error: eyre::Report) -> PyErr {
    exceptions::PyValueError::new_err(error.to_string())
}

fn py_histogram_snapshot_from_snapshot(snapshot: HistogramSnapshot) -> PyHistogramSnapshot {
    PyHistogramSnapshot {
        kind: format!("{:?}", snapshot.kind).to_lowercase(),
        title: snapshot.title,
        type_description: snapshot.type_description,
        phase: format!("{:?}", snapshot.phase).to_lowercase(),
        value_transform: format!("{:?}", snapshot.value_transform).to_lowercase(),
        supports_misbinning_mitigation: snapshot.supports_misbinning_mitigation,
        x_min: snapshot.x_min,
        x_max: snapshot.x_max,
        sample_count: snapshot.sample_count,
        log_x_axis: snapshot.log_x_axis,
        log_y_axis: snapshot.log_y_axis,
        discrete_min_bin_id: snapshot.discrete_min_bin_id,
        discrete_ordering: snapshot
            .discrete_ordering
            .map(|ordering| ordering.as_str().to_string()),
        bins: snapshot
            .bins
            .into_iter()
            .map(|bin| PyHistogramBinSnapshot {
                x_min: bin.x_min,
                x_max: bin.x_max,
                bin_id: bin.bin_id,
                label: bin.label,
                entry_count: bin.entry_count,
                sum_weights: bin.sum_weights,
                sum_weights_squared: bin.sum_weights_squared,
                mitigated_fill_count: bin.mitigated_fill_count,
            })
            .collect(),
        underflow_bin: PyHistogramBinSnapshot {
            x_min: snapshot.underflow_bin.x_min,
            x_max: snapshot.underflow_bin.x_max,
            bin_id: snapshot.underflow_bin.bin_id,
            label: snapshot.underflow_bin.label,
            entry_count: snapshot.underflow_bin.entry_count,
            sum_weights: snapshot.underflow_bin.sum_weights,
            sum_weights_squared: snapshot.underflow_bin.sum_weights_squared,
            mitigated_fill_count: snapshot.underflow_bin.mitigated_fill_count,
        },
        overflow_bin: PyHistogramBinSnapshot {
            x_min: snapshot.overflow_bin.x_min,
            x_max: snapshot.overflow_bin.x_max,
            bin_id: snapshot.overflow_bin.bin_id,
            label: snapshot.overflow_bin.label,
            entry_count: snapshot.overflow_bin.entry_count,
            sum_weights: snapshot.overflow_bin.sum_weights,
            sum_weights_squared: snapshot.overflow_bin.sum_weights_squared,
            mitigated_fill_count: snapshot.overflow_bin.mitigated_fill_count,
        },
        statistics: py_histogram_statistics_from_snapshot(snapshot.statistics),
    }
}

fn py_histogram_statistics_from_snapshot(
    stats: HistogramStatisticsSnapshot,
) -> PyHistogramStatisticsSnapshot {
    PyHistogramStatisticsSnapshot {
        in_range_entry_count: stats.in_range_entry_count,
        nan_value_count: stats.nan_value_count,
        mitigated_pair_count: stats.mitigated_pair_count,
    }
}

fn py_integral_estimate_from_estimate(
    estimate: gammalooprs::settings::runtime::IntegralEstimate,
) -> PyIntegralEstimate {
    PyIntegralEstimate {
        neval: estimate.neval,
        real_zero: estimate.real_zero,
        im_zero: estimate.im_zero,
        result: py_complex_from_complex(estimate.result),
        error: py_complex_from_complex(estimate.error),
        real_chisq: estimate.real_chisq.0,
        im_chisq: estimate.im_chisq.0,
    }
}

fn py_table_component_result_from_result(
    result: gammalooprs::settings::runtime::IntegrationTableComponentResult,
) -> PyIntegrationTableComponentResult {
    PyIntegrationTableComponentResult {
        component: result.component,
        value: result.value.0,
        error: result.error.0,
        relative_error_percent: result.relative_error_percent,
        chi_sq_per_dof: result.chi_sq_per_dof,
        target_delta_sigma: result.target_delta_sigma,
        target_delta_percent: result.target_delta_percent,
        max_weight_impact: result.max_weight_impact,
    }
}

fn py_integration_statistics_from_snapshot(
    snapshot: gammalooprs::settings::runtime::IntegrationStatisticsSnapshot,
) -> PyIntegrationStatisticsSnapshot {
    PyIntegrationStatisticsSnapshot {
        num_evals: snapshot.num_evals,
        average_total_time_seconds: snapshot.average_total_time_seconds,
        average_parameterization_time_seconds: snapshot.average_parameterization_time_seconds,
        average_integrand_time_seconds: snapshot.average_integrand_time_seconds,
        average_evaluator_time_seconds: snapshot.average_evaluator_time_seconds,
        average_observable_time_seconds: snapshot.average_observable_time_seconds,
        average_integrator_time_seconds: snapshot.average_integrator_time_seconds,
        f64_percentage: snapshot.f64_percentage,
        f128_percentage: snapshot.f128_percentage,
        arb_percentage: snapshot.arb_percentage,
        nan_percentage: snapshot.nan_percentage,
        nan_or_unstable_percentage: snapshot.nan_or_unstable_percentage,
        generated_event_count: snapshot.generated_event_count,
        accepted_event_count: snapshot.accepted_event_count,
        selection_efficiency_percentage: snapshot.selection_efficiency_percentage,
    }
}

fn py_max_weight_info_entry_from_entry(
    entry: gammalooprs::settings::runtime::MaxWeightInfoEntry,
) -> PyMaxWeightInfoEntry {
    PyMaxWeightInfoEntry {
        component: entry.component,
        sign: entry.sign,
        max_eval: entry.max_eval.0,
        coordinates: entry.coordinates,
    }
}

fn py_discrete_coordinate_from_coordinate(
    coordinate: gammalooprs::settings::runtime::DiscreteCoordinate,
) -> PyDiscreteCoordinate {
    PyDiscreteCoordinate {
        axis_label: coordinate.axis_label,
        bin_index: coordinate.bin_index,
        bin_label: coordinate.bin_label,
    }
}

fn py_discrete_breakdown_entry_from_entry(
    entry: gammalooprs::settings::runtime::DiscreteBreakdownEntry,
) -> PyDiscreteBreakdownEntry {
    PyDiscreteBreakdownEntry {
        bin_index: entry.bin_index,
        bin_label: entry.bin_label,
        pdf: entry.pdf.0,
        value: entry.value.0,
        error: entry.error.0,
        chi_sq: entry.chi_sq.0,
        processed_samples: entry.processed_samples,
    }
}

fn py_discrete_breakdown_from_breakdown(
    breakdown: gammalooprs::settings::runtime::DiscreteBreakdown,
) -> PyDiscreteBreakdown {
    PyDiscreteBreakdown {
        axis_label: breakdown.axis_label,
        fixed_coordinates: breakdown
            .fixed_coordinates
            .into_iter()
            .map(py_discrete_coordinate_from_coordinate)
            .collect(),
        entries: breakdown
            .entries
            .into_iter()
            .map(py_discrete_breakdown_entry_from_entry)
            .collect(),
    }
}

fn py_component_discrete_breakdown_from_breakdown(
    breakdown: gammalooprs::settings::runtime::ComponentDiscreteBreakdown,
) -> PyComponentDiscreteBreakdown {
    PyComponentDiscreteBreakdown {
        re: breakdown.re.map(py_discrete_breakdown_from_breakdown),
        im: breakdown.im.map(py_discrete_breakdown_from_breakdown),
    }
}

fn py_slot_integration_result_from_result(
    result: gammalooprs::settings::runtime::SlotIntegrationResult,
) -> PySlotIntegrationResult {
    PySlotIntegrationResult {
        key: result.key,
        process: result.process,
        integrand: result.integrand,
        target: result.target.map(py_complex_from_complex),
        integral: py_integral_estimate_from_estimate(result.integral),
        table_results: result
            .table_results
            .into_iter()
            .map(py_table_component_result_from_result)
            .collect(),
        integration_statistics: py_integration_statistics_from_snapshot(
            result.integration_statistics,
        ),
        max_weight_info: result
            .max_weight_info
            .into_iter()
            .map(py_max_weight_info_entry_from_entry)
            .collect(),
        grid_breakdown: py_component_discrete_breakdown_from_breakdown(result.grid_breakdown),
    }
}

fn py_integration_result_from_result(
    result: gammalooprs::settings::runtime::IntegrationResult,
) -> PyIntegrationResult {
    PyIntegrationResult {
        slots: result
            .slots
            .into_iter()
            .map(py_slot_integration_result_from_result)
            .collect(),
    }
}

fn py_integrand_graph_info_from_info(graph: IntegrandGraphInfo) -> PyIntegrandGraphInfo {
    PyIntegrandGraphInfo {
        graph_id: graph.graph_id,
        name: graph.name,
        is_master: graph.is_master,
    }
}

fn py_integrand_orientation_info_from_info(
    orientation: IntegrandOrientationInfo,
) -> PyIntegrandOrientationInfo {
    PyIntegrandOrientationInfo {
        orientation_id: orientation.orientation_id,
        signature: orientation.signature,
    }
}

fn py_integrand_loop_momentum_basis_info_from_info(
    basis: IntegrandLoopMomentumBasisInfo,
) -> PyIntegrandLoopMomentumBasisInfo {
    PyIntegrandLoopMomentumBasisInfo {
        basis_id: basis.basis_id,
        channel_id: basis.channel_id,
        edge_ids: basis.edge_ids,
        matches_generation_basis: basis.matches_generation_basis,
    }
}

fn py_integrand_cut_info_from_info(cut: IntegrandCutInfo) -> PyIntegrandCutInfo {
    PyIntegrandCutInfo {
        cut_id: cut.cut_id,
        edge_ids: cut.edge_ids,
        raising_power: cut.raising_power,
        left_thresholds: cut
            .left_thresholds
            .into_iter()
            .map(py_integrand_cut_threshold_info_from_info)
            .collect(),
        right_thresholds: cut
            .right_thresholds
            .into_iter()
            .map(py_integrand_cut_threshold_info_from_info)
            .collect(),
    }
}

fn py_integrand_cut_threshold_info_from_info(
    threshold: IntegrandCutThresholdInfo,
) -> PyIntegrandCutThresholdInfo {
    PyIntegrandCutThresholdInfo {
        esurface_id: threshold.esurface_id,
        status: threshold.status.as_str().to_string(),
        cut_boundary_edge_ids: threshold.cut_boundary_edge_ids,
        threshold_boundary_edge_ids: threshold.threshold_boundary_edge_ids,
        invariant_bound_is_applicable: threshold.invariant_bound_is_applicable,
    }
}

fn py_integrand_active_threshold_cut_info_from_info(
    cut: IntegrandActiveThresholdCutInfo,
) -> PyIntegrandActiveThresholdCutInfo {
    PyIntegrandActiveThresholdCutInfo {
        cut_id: cut.cut_id,
        can_become_pinched: cut.can_become_pinched,
    }
}

fn py_integrand_threshold_esurface_info_from_info(
    threshold: IntegrandThresholdEsurfaceInfo,
) -> PyIntegrandThresholdEsurfaceInfo {
    PyIntegrandThresholdEsurfaceInfo {
        esurface_id: threshold.esurface_id,
        representative_graph_id: threshold.representative_graph_id,
        edge_ids: threshold.edge_ids,
        classification: threshold
            .classification
            .map(|classification| classification.as_str().to_string()),
        active_cuts: threshold
            .active_cuts
            .into_iter()
            .map(py_integrand_active_threshold_cut_info_from_info)
            .collect(),
    }
}

fn py_integrand_graph_group_info_from_info(
    group: IntegrandGraphGroupInfo,
) -> PyIntegrandGraphGroupInfo {
    PyIntegrandGraphGroupInfo {
        group_id: group.group_id,
        graphs: group
            .graphs
            .into_iter()
            .map(py_integrand_graph_info_from_info)
            .collect(),
        orientation_edge_ids: group.orientation_edge_ids,
        orientations: group
            .orientations
            .into_iter()
            .map(py_integrand_orientation_info_from_info)
            .collect(),
        loop_momentum_bases: group
            .loop_momentum_bases
            .into_iter()
            .map(py_integrand_loop_momentum_basis_info_from_info)
            .collect(),
        threshold_esurface_ids: group.threshold_esurface_ids,
        threshold_esurfaces: group
            .threshold_esurfaces
            .into_iter()
            .map(py_integrand_threshold_esurface_info_from_info)
            .collect(),
        cuts: group
            .cuts
            .into_iter()
            .map(py_integrand_cut_info_from_info)
            .collect(),
    }
}

fn py_integrand_info_from_info(info: IntegrandInfo) -> PyIntegrandInfo {
    PyIntegrandInfo {
        process_id: info.process_id,
        process_name: info.process_name,
        integrand_name: info.integrand_name,
        kind: info.kind.to_string(),
        generation_backend: info
            .generation_compilation
            .active_backend_name()
            .to_string(),
        generation_compile_options: info
            .generation_compilation
            .external_options()
            .map(ToString::to_string),
        active_f64_backend: info.active_f64_backend.to_string(),
        graph_count: info.graph_count,
        graph_group_count: info.graph_group_count,
        record_size_bytes: info.record_size_bytes,
        graph_groups: info
            .graph_groups
            .into_iter()
            .map(py_integrand_graph_group_info_from_info)
            .collect(),
    }
}

#[allow(dead_code)] // Retained with the currently disabled Python `integrate` binding.
fn py_process_ref_from_any(process: &Bound<'_, PyAny>) -> PyResult<ProcessRef> {
    if let Ok(process_id) = process.extract::<usize>() {
        return Ok(ProcessRef::Id(process_id));
    }

    let process = process.extract::<String>().map_err(|_| {
        exceptions::PyTypeError::new_err(
            "process selectors must be either an integer process id or a string selector",
        )
    })?;
    ProcessRef::from_str(&process).map_err(exceptions::PyValueError::new_err)
}

#[allow(dead_code)] // Retained with the currently disabled Python `integrate` binding.
fn py_complex_target_from_any(
    target: &Bound<'_, PyAny>,
) -> PyResult<spenso::algebra::complex::Complex<gammalooprs::utils::F<f64>>> {
    if let Ok(complex) = target.cast::<PyComplex>() {
        return Ok(spenso::algebra::complex::Complex::new(
            gammalooprs::utils::F(complex.real()),
            gammalooprs::utils::F(complex.imag()),
        ));
    }

    if let Ok((re, im)) = target.extract::<(f64, f64)>() {
        return Ok(spenso::algebra::complex::Complex::new(
            gammalooprs::utils::F(re),
            gammalooprs::utils::F(im),
        ));
    }

    Err(exceptions::PyTypeError::new_err(
        "targets must be a Python complex number or a two-element numeric tuple/list",
    ))
}

#[allow(dead_code)] // Retained with the currently disabled Python `integrate` binding.
fn resolve_python_slot_key(
    state: &State,
    process: &ProcessRef,
    integrand_name: &str,
) -> PyResult<String> {
    let (process_id, resolved_integrand_name) = state
        .find_integrand_ref(Some(process), Some(&integrand_name.to_string()))
        .map_err(|err| {
            exceptions::PyException::new_err(format!("Could not resolve slot: {err}"))
        })?;
    Ok(format!(
        "{}@{}",
        state.process_list.processes[process_id]
            .definition
            .folder_name,
        resolved_integrand_name
    ))
}

#[allow(dead_code, clippy::too_many_arguments)] // Retained with the disabled Python `integrate` binding.
fn build_python_integrate_command(
    state: &State,
    slots: Option<Vec<(ProcessRef, String)>>,
    process: Option<ProcessRef>,
    integrand_name: Option<String>,
    target: Option<&Bound<'_, PyAny>>,
    n_cores: Option<usize>,
    workspace_path: Option<PathBuf>,
    restart: bool,
    uncorrelated: bool,
    show_max_weight_info: bool,
    no_show_integration_statistics: bool,
    show_phase: String,
    show_top_discrete_grid: bool,
    show_discrete_contributions_sum: bool,
    sort_contributions: String,
    show_max_weight_info_for_discrete_bins: bool,
    show_summary_only: bool,
    no_stream_iterations: bool,
    no_stream_updates: bool,
    batch_size: Option<usize>,
    batch_timing: f64,
    min_time_between_status_updates: f64,
    max_table_width: usize,
    write_results_for_each_iteration: bool,
) -> PyResult<Integrate> {
    let mut integrate = if let Some(slots) = slots {
        if process.is_some() || integrand_name.is_some() {
            return Err(exceptions::PyValueError::new_err(
                "`slots` cannot be combined with `process` or `integrand_name`",
            ));
        }
        Integrate::from_slots(slots)
    } else {
        let mut integrate = Integrate::default();
        if let Some(process) = process {
            integrate.process.push(process);
        }
        if let Some(integrand_name) = integrand_name {
            integrate.integrand_name.push(integrand_name);
        }
        integrate
    };

    integrate.n_cores = n_cores;
    integrate.workspace_path = workspace_path;
    integrate.restart = restart;
    integrate.uncorrelated = uncorrelated;
    integrate.show_max_weight_info = show_max_weight_info;
    integrate.no_show_integration_statistics = no_show_integration_statistics;
    integrate.show_phase = crate::commands::integrate::ShowPhaseOption::from_str(&show_phase, true)
        .map_err(|_| {
            exceptions::PyValueError::new_err(format!("Unknown show_phase value '{show_phase}'"))
        })?;
    integrate.show_top_discrete_grid = show_top_discrete_grid;
    integrate.show_discrete_contributions_sum = show_discrete_contributions_sum;
    integrate.sort_contributions =
        crate::commands::integrate::ContributionSortOption::from_str(&sort_contributions, true)
            .map_err(|_| {
                exceptions::PyValueError::new_err(format!(
                    "Unknown sort_contributions value '{sort_contributions}'"
                ))
            })?;
    integrate.show_max_weight_info_for_discrete_bins = show_max_weight_info_for_discrete_bins;
    integrate.show_summary_only = show_summary_only;
    integrate.no_stream_iterations = no_stream_iterations;
    integrate.no_stream_updates = no_stream_updates;
    integrate.batch_size = batch_size;
    integrate.batch_timing = batch_timing;
    integrate.min_time_between_status_updates = min_time_between_status_updates;
    integrate.max_table_width = max_table_width;
    integrate.write_results_for_each_iteration = write_results_for_each_iteration;
    integrate.renderer = crate::commands::integrate::RendererOption::Tabled;

    if let Some(target) = target {
        if let Ok(target_dict) = target.cast::<PyDict>() {
            let keyed_targets = target_dict
                .iter()
                .map(|(key, value)| {
                    let slot_key = if let Ok(slot_key) = key.extract::<String>() {
                        slot_key
                    } else {
                        let tuple = key.cast::<PyTuple>().map_err(|_| {
                            exceptions::PyTypeError::new_err(
                                "multi-slot target keys must be slot-key strings or `(process, integrand_name)` tuples",
                            )
                        })?;
                        if tuple.len() != 2 {
                            return Err(exceptions::PyValueError::new_err(
                                "multi-slot target tuple keys must contain exactly two items",
                            ));
                        }
                        let process = py_process_ref_from_any(&tuple.get_item(0)?)?;
                        let integrand_name = tuple.get_item(1)?.extract::<String>()?;
                        resolve_python_slot_key(state, &process, &integrand_name)?
                    };
                    Ok((slot_key, py_complex_target_from_any(&value)?))
                })
                .collect::<PyResult<Vec<_>>>()?;
            integrate = integrate.with_keyed_targets(keyed_targets);
        } else {
            integrate = integrate.with_single_target(py_complex_target_from_any(target)?);
        }
    }

    Ok(integrate)
}

/*
pub struct OutputOptions {}

/// There should only be 1 worker per instance of gammaloop.
#[pyclass(from_py_object, name = "Worker", unsendable)]
pub struct PythonWorker {
    pub model: Model,
    pub process_list: ProcessList,
    // pub integrands: HashMap<String, Integrand>,
    pub master_node: Option<MasterNode>,
}

impl Clone for PythonWorker {
    fn clone(&self) -> PythonWorker {
        PythonWorker {
            model: self.model.clone(),
            process_list: self.process_list.clone(),
            // integrands: self.integrands.clone(),
            master_node: self.master_node.clone(),
        }
    }
}

#[pyclass(from_py_object, name = "SnailFilterOptions")]
pub struct PySnailFilterOptions {
    pub filter_options: SnailFilterOptions,
}

#[pymethods]
impl PySnailFilterOptions {
    #[new]
    pub(crate) fn __new__(
        veto_snails_attached_to_massive_lines: Option<bool>,
        veto_snails_attached_to_massless_lines: Option<bool>,
        veto_only_scaleless_snails: Option<bool>,
    ) -> Result<PySnailFilterOptions> {
        Ok(PySnailFilterOptions {
            filter_options: SnailFilterOptions {
                veto_snails_attached_to_massive_lines: veto_snails_attached_to_massive_lines
                    .unwrap_or(false),
                veto_snails_attached_to_massless_lines: veto_snails_attached_to_massless_lines
                    .unwrap_or(true),
                veto_only_scaleless_snails: veto_only_scaleless_snails.unwrap_or(false),
            },
        })
    }

    pub(crate) fn __str__(&self) -> Result<String> {
        Ok(format!("{}", self.filter_options))
    }
}

#[pyclass(from_py_object, name = "SelfEnergyFilterOptions")]
pub struct PySelfEnergyFilterOptions {
    pub filter_options: SelfEnergyFilterOptions,
}

#[pymethods]
impl PySelfEnergyFilterOptions {
    #[new]
    pub(crate) fn __new__(
        veto_self_energy_of_massive_lines: Option<bool>,
        veto_self_energy_of_massless_lines: Option<bool>,
        veto_only_scaleless_self_energy: Option<bool>,
    ) -> Result<PySelfEnergyFilterOptions> {
        Ok(PySelfEnergyFilterOptions {
            filter_options: SelfEnergyFilterOptions {
                veto_self_energy_of_massive_lines: veto_self_energy_of_massive_lines
                    .unwrap_or(true),
                veto_self_energy_of_massless_lines: veto_self_energy_of_massless_lines
                    .unwrap_or(true),
                veto_only_scaleless_self_energy: veto_only_scaleless_self_energy.unwrap_or(false),
            },
        })
    }

    pub(crate) fn __str__(&self) -> Result<String> {
        Ok(format!("{}", self.filter_options))
    }
}

#[pyclass(from_py_object, name = "SewedFilterOptions")]
pub struct PySewedFilterOptions {
    pub filter_options: SewedFilterOptions,
}

#[pymethods]
impl PySewedFilterOptions {
    #[new]
    pub(crate) fn __new__(filter_tadpoles: Option<bool>) -> Result<PySewedFilterOptions> {
        Ok(PySewedFilterOptions {
            filter_options: SewedFilterOptions {
                filter_tadpoles: filter_tadpoles.unwrap_or(false),
            },
        })
    }
}

#[pyclass(from_py_object, name = "TadpolesFilterOptions")]
pub struct PyTadpolesFilterOptions {
    pub filter_options: TadpolesFilterOptions,
}

#[pymethods]
impl PyTadpolesFilterOptions {
    #[new]
    pub(crate) fn __new__(
        veto_tadpoles_attached_to_massive_lines: Option<bool>,
        veto_tadpoles_attached_to_massless_lines: Option<bool>,
        veto_only_scaleless_tadpoles: Option<bool>,
    ) -> Result<PyTadpolesFilterOptions> {
        Ok(PyTadpolesFilterOptions {
            filter_options: TadpolesFilterOptions {
                veto_tadpoles_attached_to_massive_lines: veto_tadpoles_attached_to_massive_lines
                    .unwrap_or(true),
                veto_tadpoles_attached_to_massless_lines: veto_tadpoles_attached_to_massless_lines
                    .unwrap_or(true),
                veto_only_scaleless_tadpoles: veto_only_scaleless_tadpoles.unwrap_or(false),
            },
        })
    }

    pub(crate) fn __str__(&self) -> Result<String> {
        Ok(format!("{}", self.filter_options))
    }
}

#[pyclass(from_py_object, name = "FeynGenFilters")]
pub struct PyFeynGenFilters {
    pub filters: Vec<FeynGenFilter>,
}
impl<'a> FromPyObject<'a> for PyFeynGenFilters {
    fn extract_bound(ob: &Bound<'a, PyAny>) -> PyResult<Self> {
        if let Ok(a) = ob.extract::<PyFeynGenFilters>() {
            Ok(PyFeynGenFilters { filters: a.filters })
        } else {
            Err(exceptions::PyValueError::new_err(
                "Not a valid Feynman generation filter",
            ))
        }
    }
}

#[pymethods]
impl PyFeynGenFilters {
    pub(crate) fn __str__(&self) -> Result<String> {
        Ok(self.filters.iter().map(|f| format!(" > {}", f)).join("\n"))
    }

    #[new]
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn __new__(
        particle_veto: Option<Vec<i64>>,
        max_number_of_bridges: Option<usize>,
        sewed_filter: Option<PyRef<PySewedFilterOptions>>,
        self_energy_filter: Option<PyRef<PySelfEnergyFilterOptions>>,
        tadpoles_filter: Option<PyRef<PyTadpolesFilterOptions>>,
        zero_snails_filter: Option<PyRef<PySnailFilterOptions>>,
        perturbative_orders: Option<HashMap<String, usize>>,
        coupling_orders: Option<HashMap<String, (usize, Option<usize>)>>,
        loop_count_range: Option<(usize, usize)>,
        fermion_loop_count_range: Option<(usize, usize)>,
        factorized_loop_topologies_count_range: Option<(usize, usize)>,
    ) -> Result<PyFeynGenFilters> {
        let mut filters = Vec::new();
        if let Some(sewed_filter) = sewed_filter {
            filters.push(FeynGenFilter::SewedFilter(sewed_filter.filter_options));
        }
        if let Some(self_energy_filter) = self_energy_filter {
            filters.push(FeynGenFilter::SelfEnergyFilter(
                self_energy_filter.filter_options,
            ));
        }
        if let Some(particle_veto) = particle_veto {
            filters.push(FeynGenFilter::ParticleVeto(particle_veto));
        }
        if let Some(max_number_of_bridges) = max_number_of_bridges {
            filters.push(FeynGenFilter::MaxNumberOfBridges(max_number_of_bridges));
        }
        if let Some(tadpoles_filter) = tadpoles_filter {
            filters.push(FeynGenFilter::TadpolesFilter(
                tadpoles_filter.filter_options,
            ));
        }
        if let Some(zero_snails_filter) = zero_snails_filter {
            filters.push(FeynGenFilter::ZeroSnailsFilter(
                zero_snails_filter.filter_options,
            ));
        }
        if let Some(perturbative_orders) = perturbative_orders {
            filters.push(FeynGenFilter::PerturbativeOrders(perturbative_orders));
        }
        if let Some(coupling_orders) = coupling_orders {
            filters.push(FeynGenFilter::CouplingOrders(coupling_orders));
        }
        if let Some(loop_count_range) = loop_count_range {
            filters.push(FeynGenFilter::LoopCountRange(loop_count_range));
        }
        if let Some(fermion_loop_count_range) = fermion_loop_count_range {
            filters.push(FeynGenFilter::FermionLoopCountRange(
                fermion_loop_count_range,
            ));
        }
        if let Some(factorized_loop_topologies_count_range) = factorized_loop_topologies_count_range
        {
            filters.push(FeynGenFilter::FactorizedLoopTopologiesCountRange(
                factorized_loop_topologies_count_range,
            ));
        }
        Ok(PyFeynGenFilters { filters })
    }
}

fn feyngen_to_python_error(error: FeynGenError) -> PyErr {
    exceptions::PyValueError::new_err(format!("Feynman diagram generator error | {error}"))
}

#[pyclass(from_py_object, name = "NumeratorAwareGroupingOption")]
pub struct PyNumeratorAwareGroupingOption {
    pub grouping_options: NumeratorAwareGraphGroupingOption,
}

#[pymethods]
impl PyNumeratorAwareGroupingOption {
    #[new]
    pub(crate) fn __new__(
        numerator_aware_grouping_option: Option<String>,
        compare_canonized_numerator: Option<bool>,
        number_of_samples_for_numerator_comparisons: Option<usize>,
        consider_internal_masses_only_in_numerator_isomorphisms: Option<bool>,
        fully_numerical_substitution_when_comparing_numerators: Option<bool>,
        numerical_samples_seed: Option<u16>,
        symmetric_polarizations: Option<bool>,
    ) -> Result<PyNumeratorAwareGroupingOption> {
        Ok(PyNumeratorAwareGroupingOption {
            grouping_options: NumeratorAwareGraphGroupingOption::new_with_attributes(
                numerator_aware_grouping_option
                    .unwrap_or("group_identical_graphs_up_to_scalar_rescaling".into())
                    .as_str(),
                numerical_samples_seed,
                number_of_samples_for_numerator_comparisons,
                consider_internal_masses_only_in_numerator_isomorphisms,
                fully_numerical_substitution_when_comparing_numerators,
                compare_canonized_numerator,
                symmetric_polarizations,
            )
            .map_err(|e| exceptions::PyException::new_err(e.to_string()))?,
        })
    }

    pub(crate) fn __str__(&self) -> Result<String> {
        Ok(format!("{}", self.grouping_options))
    }
}
*/

/// Load, inspect, and evaluate one mutable GammaLoop session.
///
/// One instance owns a GammaLoop state, run history, CLI settings, default runtime
/// settings, and session state. Calls to ``run`` and the evaluation methods all act
/// on that same in-memory session.
///
/// Notes
/// -----
/// ``read_only_state=True`` prevents writes inside the active state directory; it
/// does not make this Python object immutable. Commands may still change in-memory
/// settings, processes, or run history. Create separate instances when independent
/// sessions are required.
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(unsendable, name = "GammaLoopAPI")
)]
struct GammaLoopAPI {
    gammaloop_state: State,
    cli_settings: CLISettings,
    run_history: RunHistory,
    default_runtime_settings: RuntimeSettings,
    session_state: CliSessionState,
}

// TODO: Improve error broadcasting to Python everywhere so as to show rust backtrace
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl GammaLoopAPI {
    /// Load or create a GammaLoop state and initialize its CLI session.
    ///
    /// Parameters
    /// ----------
    /// state_folder : path-like, optional
    ///     State directory to load or create. The default is ``./gammaloop_state``.
    /// boot_commands_path : path-like, optional
    ///     TOML run card whose commands are applied during startup.
    /// model_file : path-like, optional
    ///     Model file override used while loading or initializing the state.
    /// trace_logs_filename : str, optional
    ///     File receiving native trace records for this session.
    /// level : LogLevel, optional
    ///     Terminal log-level override for this session.
    /// logfile_level : LogLevel, optional
    ///     File log-level override for this session.
    /// logging_prefix : object, optional
    ///     Native logging-prefix configuration.
    /// read_only_state : bool, default=False
    ///     Prevent writes whose target lies inside the active state directory and
    ///     disable file logging there. In-memory session changes remain possible.
    /// settings_global_path : path-like, optional
    ///     TOML file overriding the global settings loaded at startup.
    /// settings_runtime_defaults_path : path-like, optional
    ///     TOML file overriding the default runtime settings loaded at startup.
    /// clean_state : bool, default=False
    ///     Remove the resolved state path before startup. This is destructive and
    ///     cannot be combined with ``read_only_state=True``.
    ///
    /// Returns
    /// -------
    /// GammaLoopAPI
    ///     A stateful API instance sharing one state, history, and settings session.
    ///
    /// Raises
    /// ------
    /// Exception
    ///     If startup options conflict, the state or settings cannot be loaded, a
    ///     boot command fails, or a boot card requests process exit.
    ///
    /// Examples
    /// --------
    /// Open an existing generated state without permitting writes to it::
    ///
    ///     api = GammaLoopAPI(state_folder="./state", read_only_state=True)
    #[new]
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature = (
        state_folder=None,
        boot_commands_path=None,
        model_file=None,
        trace_logs_filename=None,
        level=None,
        logfile_level=None,
        logging_prefix=None,
        read_only_state=false,
        settings_global_path=None,
        settings_runtime_defaults_path=None,
        clean_state=false
    ))]
    pub fn new_python(
        state_folder: Option<PathBuf>,
        boot_commands_path: Option<PathBuf>,
        model_file: Option<PathBuf>,
        trace_logs_filename: Option<String>,
        level: Option<LogLevel>,
        logfile_level: Option<LogLevel>,
        #[gen_stub(override_type(
            type_repr = "builtins.object | None",
            imports = ("builtins")
        ))]
        logging_prefix: Option<LogFormat>,
        read_only_state: bool,
        settings_global_path: Option<PathBuf>,
        settings_runtime_defaults_path: Option<PathBuf>,
        clean_state: bool,
    ) -> Result<Self> {
        let LoadedState {
            state,
            run_history,
            cli_settings,
            default_runtime_settings,
            session_state,
            ..
        } = StateLoadOption {
            clean_state,
            boot_commands_path,
            state_folder,
            model_file,
            trace_logs_filename,
            level,
            logfile_level,
            logging_prefix,
            read_only_state,
            settings_global_path,
            settings_runtime_defaults_path,
        }
        .load()
        .map_err(|e| {
            exceptions::PyException::new_err(format!(
                "Could not load or create GammaLoop state: {}",
                e
            ))
        })?;
        Ok(GammaLoopAPI {
            gammaloop_state: state,
            run_history,
            cli_settings,
            default_runtime_settings,
            session_state,
        })
    }

    /// Evaluate one integration or momentum-space sample.
    ///
    /// Parameters
    /// ----------
    /// point : Sequence[float]
    ///     Coordinates for one sample. In integration space, the length must match
    ///     the selected integrand and ``discrete_dim``. In momentum space, values
    ///     are grouped as ``(px, py, pz)`` with one triplet per independent loop
    ///     momentum. Energy components and external momenta are not accepted.
    /// process_id : int, optional
    ///     Process containing the integrand. Supply this when selection is ambiguous.
    /// integrand_name : str, optional
    ///     Integrand to evaluate. Supply this when selection is ambiguous.
    /// use_arb_prec : bool, default=False
    ///     Legacy compatibility option selecting the configured ``f128`` stability
    ///     level, or arbitrary precision when no ``f128`` level is available.
    ///     Returned numeric fields remain ``float64``.
    /// minimal_output : bool, default=False
    ///     Omit the optional evaluation metadata from the returned sample.
    /// return_events : bool, optional
    ///     Temporarily override event generation for this call. The integrand setting
    ///     is restored afterward.
    /// momentum_space : bool, default=False
    ///     Interpret ``point`` as consecutive spatial loop-momentum ``(px, py, pz)``
    ///     triplets instead of integration-space coordinates.
    /// integrator_weight : float, optional
    ///     Weight associated with this sample. The default is 1.0.
    /// discrete_dim : Sequence[int], optional
    ///     Discrete integration coordinates used to determine the expected dimension.
    /// graph_name : str, optional
    ///     Graph selected for momentum-space evaluation.
    /// orientation : int, optional
    ///     Orientation index for ``graph_name`` in momentum-space evaluation.
    ///
    /// Returns
    /// -------
    /// EvaluationResult
    ///     The sample evaluation and the observable snapshot for its one-row batch.
    ///
    /// Raises
    /// ------
    /// Exception
    ///     If integrand selection is ambiguous or invalid, dimensions do not match,
    ///     graph or orientation selection is invalid, or warm-up/evaluation fails.
    ///
    /// Notes
    /// -----
    /// With ``use_arb_prec=False``, evaluation follows the configured ``f64``,
    /// ``f128``, and arbitrary-precision stability ladder. Despite its historical
    /// name, ``use_arb_prec=True`` selects configured ``f128`` and falls back to
    /// arbitrary precision only when that level is absent. Python-visible numeric
    /// fields use the package's ``float64`` output contract. Evaluation may
    /// warm the integrand and update in-memory caches or observable snapshots even in
    /// a read-only-state session.
    ///
    /// See Also
    /// --------
    /// GammaLoop's sample-evaluation contract in the interface guide and the
    /// maintained events-and-observables example.
    ///
    /// Examples
    /// --------
    /// Evaluate a point whose dimension matches integrand 0::
    ///
    ///     result = api.evaluate_sample(point, process_id=0, minimal_output=True)
    ///     value = result.integrand_result
    #[allow(clippy::too_many_arguments)]
    #[pyo3(
        name = "evaluate_sample",
        signature = (point, process_id=None, integrand_name=None, use_arb_prec=false, minimal_output=false, return_events=None, momentum_space=false, integrator_weight=None, discrete_dim=None, graph_name=None, orientation=None)
    )]
    pub fn evaluate_sample<'py>(
        &mut self,
        _py: Python<'py>,
        point: Vec<f64>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        use_arb_prec: bool,
        minimal_output: bool,
        return_events: Option<bool>,
        momentum_space: bool,
        integrator_weight: Option<f64>,
        discrete_dim: Option<Vec<usize>>,
        graph_name: Option<String>,
        orientation: Option<usize>,
    ) -> Result<PyEvaluationResult> {
        let points =
            ndarray::Array2::from_shape_vec((1, point.len()), point).map_err(|e| eyre!(e))?;
        let integrator_weights = integrator_weight
            .map(|weight| ndarray::Array1::from_shape_vec(1, vec![weight]))
            .transpose()
            .map_err(|e| eyre!(e))?;
        let discrete_dims = discrete_dim.unwrap_or_default();
        let discrete_dims = Some(
            ndarray::Array2::from_shape_vec((1, discrete_dims.len()), discrete_dims)
                .map_err(|e| eyre!(e))?,
        );
        let res = EvaluateSamples {
            process_id,
            integrand_name,
            use_arb_prec,
            minimal_output,
            return_generated_events: return_events,
            momentum_space,
            points: points.view(),
            integrator_weights: integrator_weights.as_ref().map(|weights| weights.view()),
            discrete_dims: discrete_dims.as_ref().map(|dims| dims.view()),
            graph_names: Some(vec![graph_name]),
            orientations: Some(vec![orientation]),
        }
        .run(&mut self.gammaloop_state)?;
        let gammalooprs::integrands::evaluation::BatchSampleEvaluationResult {
            mut samples,
            observables,
        } = res;
        let value = samples.pop().ok_or_else(|| {
            eyre!("evaluate_sample did not return any result for the single input sample")
        })?;
        Ok(PyEvaluationResult {
            inner: SingleSampleEvaluationResult {
                sample: value,
                observables,
            },
        })
    }

    /// Evaluate a batch of integration or momentum-space samples.
    ///
    /// Parameters
    /// ----------
    /// points : numpy.ndarray[numpy.float64]
    ///     Two-dimensional array with one sample per row. Integration-space columns
    ///     must match the selected integrand; momentum-space columns are flattened
    ///     ``(px, py, pz)`` groups with one triplet per independent loop momentum.
    ///     Energy components and external momenta are not accepted.
    /// process_id : int, optional
    ///     Process containing the integrand. Supply this when selection is ambiguous.
    /// integrand_name : str, optional
    ///     Integrand to evaluate. Supply this when selection is ambiguous.
    /// use_arb_prec : bool, default=False
    ///     Legacy compatibility option selecting the configured ``f128`` stability
    ///     level, or arbitrary precision when no ``f128`` level is available.
    ///     Returned numeric fields remain ``float64``.
    /// minimal_output : bool, default=False
    ///     Omit the optional evaluation metadata from every returned sample.
    /// return_events : bool, optional
    ///     Temporarily override event generation for this call. The integrand setting
    ///     is restored afterward.
    /// momentum_space : bool, default=False
    ///     Interpret each row as consecutive spatial loop-momentum ``(px, py, pz)``
    ///     triplets instead of integration-space coordinates.
    /// integrator_weights : numpy.ndarray[numpy.float64], optional
    ///     One weight per row. Defaults to 1.0 for every sample.
    /// discrete_dims : numpy.ndarray[numpy.unsignedinteger], optional
    ///     Two-dimensional array with one row of discrete coordinates per sample.
    /// graph_names : Sequence[str | None], optional
    ///     Momentum-space graph selection for each sample.
    /// orientations : Sequence[int | None], optional
    ///     Momentum-space orientation selection for each sample.
    ///
    /// Returns
    /// -------
    /// BatchEvaluationResult
    ///     Per-sample evaluations and one observable snapshot for the complete batch.
    ///
    /// Raises
    /// ------
    /// Exception
    ///     If integrand selection is ambiguous or invalid, array or option lengths do
    ///     not match, dimensions are invalid, or warm-up/evaluation fails.
    ///
    /// Notes
    /// -----
    /// With ``use_arb_prec=False``, evaluation follows the configured ``f64``,
    /// ``f128``, and arbitrary-precision stability ladder. Despite its historical
    /// name, ``use_arb_prec=True`` selects configured ``f128`` and falls back to
    /// arbitrary precision only when that level is absent. Python-visible numeric
    /// fields use the package's ``float64`` output contract. Evaluation may
    /// update in-memory caches or observable snapshots in a read-only-state session.
    ///
    /// See Also
    /// --------
    /// GammaLoop's sample-evaluation contract in the interface guide and the
    /// maintained events-and-observables example.
    ///
    /// Examples
    /// --------
    /// Evaluate a caller-provided two-dimensional array::
    ///
    ///     points = numpy.asarray(points, dtype=numpy.float64)
    ///     result = api.evaluate_samples(points, process_id=0, minimal_output=True)
    #[allow(clippy::too_many_arguments)]
    #[pyo3(
        name = "evaluate_samples",
        signature = (points, process_id=None, integrand_name=None, use_arb_prec=false, minimal_output=false, return_events=None, momentum_space=false, integrator_weights=None, discrete_dims=None, graph_names=None, orientations=None)
    )]
    pub fn evaluate_samples<'py>(
        &mut self,
        _py: Python<'py>,
        points: PyReadonlyArray2<'py, f64>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        use_arb_prec: bool,
        minimal_output: bool,
        return_events: Option<bool>,
        momentum_space: bool,
        integrator_weights: Option<PyReadonlyArray1<'py, f64>>,
        #[gen_stub(override_type(
            type_repr = "numpy.typing.NDArray[numpy.unsignedinteger] | None",
            imports = ("numpy", "numpy.typing")
        ))]
        discrete_dims: Option<PyReadonlyArray2<'py, usize>>,
        graph_names: Option<Vec<Option<String>>>,
        orientations: Option<Vec<Option<usize>>>,
    ) -> Result<PyBatchEvaluationResult> {
        let points_rust = points.as_array();

        let evaluate_samples = EvaluateSamples {
            process_id,
            integrand_name,
            use_arb_prec,
            minimal_output,
            return_generated_events: return_events,
            momentum_space,
            points: points_rust,
            integrator_weights: integrator_weights.as_ref().map(PyReadonlyArray1::as_array),
            discrete_dims: discrete_dims.as_ref().map(PyReadonlyArray2::as_array),
            graph_names,
            orientations,
        };

        let res = evaluate_samples.run(&mut self.gammaloop_state)?;
        Ok(PyBatchEvaluationResult { inner: res })
    }

    /// Import DOT graphs into a new or existing process/integrand collection.
    ///
    /// Parameters
    /// ----------
    /// graphs : str
    ///     DOT file path when ``format="dot"`` or inline DOT text when
    ///     ``format="string"``.
    /// process_name, process_id : str or int, optional
    ///     Process to create or update. Inline text requires one of these selectors.
    /// integrand_name : str, optional
    ///     Integrand within the selected process.
    /// format : {"dot", "string"}, default="dot"
    ///     Select file-backed or inline input.
    /// overwrite, append : bool, default=False
    ///     Replace an existing collection or append to it; these modes conflict.
    ///
    /// Raises
    /// ------
    /// Exception
    ///     If selectors conflict, the source is missing or malformed, or importing fails.
    #[allow(clippy::too_many_arguments)]
    #[pyo3(name="import_graphs", signature = (graphs, process_name=None, process_id=None, integrand_name=None, format="dot".into(), overwrite=false, append=false))]
    pub(crate) fn import_graphs_python(
        &mut self,
        graphs: String,
        process_name: Option<String>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        format: String,
        overwrite: bool,
        append: bool,
    ) -> Result<()> {
        if overwrite && append {
            return Err(eyre!(
                "Cannot use both overwrite and append options when importing graphs"
            ));
        }
        let (graphs, process_name) = match format.as_str() {
            "dot" => {
                let dot_path = PathBuf::from(&graphs);
                if !dot_path.exists() {
                    return Err(eyre!("Path does not exist: {}", dot_path.display()));
                }

                let process_name = process_name
                    .unwrap_or(dot_path.file_stem().unwrap().to_string_lossy().into_owned());

                let graphs = Graph::from_path(&dot_path, &self.gammaloop_state.model)
                    .map_err(|e| eyre!("Could not parse graphs from path: {}", e.to_string()))?;
                (graphs, Some(process_name))
            }
            "string" => {
                if process_id.is_none() && process_name.is_none() {
                    return Err(eyre!(
                        "When importing graphs from string, either a process id or a process name must be provided"
                    ));
                }

                let graphs = Graph::from_string(&graphs, &self.gammaloop_state.model)
                    .map_err(|e| eyre!("Could not parse graphs from string: {}", e.to_string()))?;
                (graphs, process_name)
            }
            other => {
                return Err(eyre!(
                    "Unknown graph import format: {}. Supported formats are 'yaml' and 'json'",
                    other
                ));
            }
        };

        self.gammaloop_state.import_graphs(
            graphs,
            process_name,
            process_id,
            integrand_name,
            overwrite,
            append,
        )
    }

    /// Generate loop-momentum bases for each supplied graph.
    ///
    /// Parameters
    /// ----------
    /// graphs : str
    ///     DOT file path or inline DOT text, as selected by ``format``.
    /// format : {"dot", "string"}, default="dot"
    ///     Select file-backed or inline input.
    ///
    /// Returns
    /// -------
    /// list
    ///     Per graph, a list of ``(loop_edges, external_edges, edge_signatures)``
    ///     tuples. Each signature contains its loop and external coefficients.
    #[allow(clippy::type_complexity)]
    #[pyo3(name="get_lmbs", signature = (graphs, format="dot".into()))]
    pub(crate) fn get_lmbs(
        &self,
        graphs: String,
        format: String,
    ) -> Result<
        Vec<
            Vec<(
                Vec<usize>, // loop edges
                Vec<usize>, // external edges
                HashMap<usize, (Vec<i8>, Vec<i8>)>,
            )>,
        >,
    > {
        let graphs = match format.as_str() {
            "dot" => {
                let dot_path = PathBuf::from(&graphs);
                if !dot_path.exists() {
                    return Err(eyre!("Path does not exist: {}", dot_path.display()));
                }

                let graphs = Graph::from_path(&dot_path, &self.gammaloop_state.model)
                    .map_err(|e| eyre!("Could not parse graphs from path: {}", e.to_string()))?;
                graphs
            }
            "string" => {
                let graphs = Graph::from_string(&graphs, &self.gammaloop_state.model)
                    .map_err(|e| eyre!("Could not parse graphs from string: {}", e.to_string()))?;
                graphs
            }
            other => {
                return Err(eyre!(
                    "Unknown graph import format: {}. Supported formats are sting and dot",
                    other
                ));
            }
        };

        let lmbs = graphs
            .into_iter()
            .map(|g| g.generate_loop_momentum_bases_of(&g.no_dummy()))
            .collect_vec();

        Ok(lmbs
            .into_iter()
            .map(|lmb_list| {
                lmb_list
                    .into_iter()
                    .map(|lmb| {
                        let loop_edges = lmb.loop_edges.into_iter().map(|e| e.0).collect_vec();
                        let external_edges = lmb.ext_edges.into_iter().map(|e| e.0).collect_vec();
                        let mut edge_signatures = HashMap::new();

                        for (edge_id, signature) in lmb.edge_signatures.into_iter() {
                            let loop_part = signature
                                .internal
                                .into_iter()
                                .map(|e| e as i8)
                                .collect_vec();

                            let ext_part = signature
                                .external
                                .into_iter()
                                .map(|e| e as i8)
                                .collect_vec();

                            edge_signatures.insert(edge_id.0, (loop_part, ext_part));
                        }
                        (loop_edges, external_edges, edge_signatures)
                    })
                    .collect()
            })
            .collect())
    }

    /// Return the causal-flow orientations generated for one graph.
    ///
    /// Each returned dictionary maps an edge id to ``1`` (default), ``-1``
    /// (reversed), or ``0`` (undirected). Supply process and integrand selectors when
    /// the active state does not identify a unique integrand.
    ///
    /// Parameters
    /// ----------
    /// graph_name : str
    ///     Name of the graph within the selected integrand.
    /// process_id : int, optional
    ///     Numeric process identifier; omit when process selection is unambiguous.
    /// integrand_name : str, optional
    ///     Integrand containing the graph; omit when integrand selection is unambiguous.
    ///
    /// Returns
    /// -------
    /// list[dict[int, int]]
    ///     One edge-direction mapping per generated orientation.
    #[pyo3(name="get_orientations", signature = (graph_name, process_id=None, integrand_name=None))]
    pub(crate) fn get_orientations(
        &self,
        graph_name: String,
        process_id: Option<usize>,
        integrand_name: Option<String>,
    ) -> Result<Vec<HashMap<usize, i8>>> {
        let (pid, name) = self
            .gammaloop_state
            .process_list
            .find_integrand(process_id, integrand_name.as_ref())
            .map_err(|e| {
                exceptions::PyException::new_err(format!("Could not find integrand: {}", e))
            })?;

        let orientations = match &self.gammaloop_state.process_list.processes[pid].collection {
            ProcessCollection::Amplitudes(amplitudes) => {
                let cff = amplitudes
                    .get(&name)
                    .unwrap()
                    .graphs
                    .iter()
                    .find(|g| g.graph.name == graph_name)
                    .as_ref()
                    .unwrap()
                    .derived_data
                    .cff_expression
                    .as_ref()
                    .unwrap();

                cff.orientations
                    .iter()
                    .map(|or_data| or_data.data.orientation.clone())
                    .collect_vec()
            }

            ProcessCollection::CrossSections(cross_sections) => {
                let cff = cross_sections
                    .get(&name)
                    .unwrap()
                    .supergraphs
                    .iter()
                    .find(|g| g.graph.name == graph_name)
                    .as_ref()
                    .unwrap()
                    .derived_data
                    .global_cff_expression
                    .as_ref()
                    .unwrap();

                cff.orientations
                    .iter()
                    .map(|or_data| or_data.data.orientation.clone())
                    .collect_vec()
            }
        };

        Ok(orientations
            .into_iter()
            .map(|orientation| {
                let mut result = HashMap::new();
                for (edge_id, direction) in orientation.into_iter() {
                    let direction = match direction {
                        Orientation::Default => 1,
                        Orientation::Reversed => -1,
                        Orientation::Undirected => 0,
                    };
                    result.insert(edge_id.0, direction);
                }
                result
            })
            .collect())
    }

    /// Serialize the active physics model as JSON.
    ///
    /// Returns
    /// -------
    /// str
    ///     JSON representation of the model currently owned by this session.
    #[pyo3(name = "get_model")]
    pub(crate) fn get_model(&self) -> PyResult<String> {
        let serializable_model = self.gammaloop_state.model.to_serializable();
        serde_json::to_string(&serializable_model).map_err(|e| {
            exceptions::PyException::new_err(format!("Could not serialize model: {}", e))
        })
    }

    /// Evaluate one generated graph group as a symbolic or numerical expression.
    ///
    /// Parameters
    /// ----------
    /// process_id : int, optional
    ///     Process containing the requested graph group.
    /// graphs_group_name : str, optional
    ///     Group to evaluate when the current state is ambiguous.
    /// result_path : path-like, optional
    ///     Optional destination for the evaluated expression.
    /// numerical : bool, default=True
    ///     Evaluate numerically instead of retaining a symbolic result.
    /// number_of_terms_in_epsilon_expansion : int, optional
    ///     Truncate the dimensional-regulator expansion to this many terms.
    ///
    /// Returns
    /// -------
    /// str
    ///     Canonical Symbolica representation of the result.
    #[pyo3(name="evaluate", signature = (process_id=None, graphs_group_name=None, result_path=None, numerical=true, number_of_terms_in_epsilon_expansion=None))]
    pub(crate) fn evaluate_python(
        &mut self,
        process_id: Option<usize>,
        graphs_group_name: Option<String>,
        result_path: Option<PathBuf>,
        numerical: bool,
        number_of_terms_in_epsilon_expansion: Option<usize>,
    ) -> PyResult<String> {
        Evaluate {
            process: process_id.map(ProcessRef::Id),
            graphs_group_name,
            result_path,
            numerical,
            number_of_terms_in_epsilon_expansion,
        }
        .run(
            &mut self.gammaloop_state,
            &self.cli_settings,
            &self.default_runtime_settings,
        )
        .map_err(|e| {
            exceptions::PyException::new_err(format!(
                "Could not load or create GammaLoop state: {}",
                e
            ))
        })
        .map(|res| res.to_canonical_string())
    }

    /// Replace the active model from a GammaLoop model file or supported model source.
    ///
    /// ``simplify_model`` applies the standard symbolic simplification pass while
    /// importing. Existing process data may no longer be compatible with a replaced
    /// model, so import the model before generating processes.
    ///
    /// Parameters
    /// ----------
    /// model_specifier : path-like
    ///     Path or model specifier accepted by GammaLoop's model importer.
    /// simplify_model : bool, default=True
    ///     Apply the standard symbolic simplification pass while importing.
    #[pyo3(name="import_model", signature = (model_specifier, simplify_model=true))]
    pub(crate) fn import_model_python(
        &mut self,
        model_specifier: PathBuf,
        simplify_model: bool,
    ) -> PyResult<()> {
        ImportModel {
            path: model_specifier,
            simplify_model,
        }
        .run(&mut self.gammaloop_state)
        .map_err(|e| exceptions::PyException::new_err(format!("Could not import model: {}", e)))
    }

    /// List generated amplitude and cross-section names with their process ids.
    ///
    /// Returns
    /// -------
    /// tuple[dict[str, int], dict[str, int]]
    ///     Amplitude mapping followed by the cross-section mapping.
    #[pyo3(name="list_outputs", signature = ())]
    pub(crate) fn list_outputs(
        &mut self,
    ) -> PyResult<(HashMap<String, usize>, HashMap<String, usize>)> {
        let mut all_amplitudes = HashMap::new();
        let mut all_cross_sections = HashMap::new();
        for (p_id, p) in self
            .gammaloop_state
            .process_list
            .processes
            .iter()
            .enumerate()
        {
            match &p.collection {
                ProcessCollection::Amplitudes(amplitudes) => {
                    amplitudes.keys().cloned().for_each(|amp| {
                        all_amplitudes.insert(amp, p_id);
                    });
                }
                ProcessCollection::CrossSections(cross_sections) => {
                    cross_sections.keys().cloned().for_each(|cs| {
                        all_cross_sections.insert(cs, p_id);
                    });
                }
            }
        }
        Ok((all_amplitudes, all_cross_sections))
    }

    /// Describe the selected generated integrand and its graph structure.
    ///
    /// Parameters
    /// ----------
    /// process_id : int, optional
    ///     Process containing the integrand. Supply this when selection is ambiguous.
    /// integrand_name : str, optional
    ///     Integrand to inspect. Supply this when selection is ambiguous.
    ///
    /// Returns
    /// -------
    /// IntegrandInfo
    ///     Structured process, backend, graph, orientation, cut, and size metadata.
    ///
    /// Raises
    /// ------
    /// Exception
    ///     If no unique generated integrand matches the selection.
    ///
    /// Examples
    /// --------
    /// Inspect the active evaluation backend and graph count::
    ///
    ///     info = api.get_integrand_info(process_id=0)
    ///     print(info.active_f64_backend, info.graph_count)
    #[pyo3(name="get_integrand_info", signature = (process_id=None, integrand_name=None))]
    pub(crate) fn get_integrand_info(
        &self,
        process_id: Option<usize>,
        integrand_name: Option<String>,
    ) -> PyResult<PyIntegrandInfo> {
        let process_ref = process_id.map(ProcessRef::Id);
        let info = self
            .gammaloop_state
            .get_integrand_info(process_ref.as_ref(), integrand_name.as_ref())
            .map_err(|e| {
                exceptions::PyException::new_err(format!("Could not get integrand info: {}", e))
            })?;
        Ok(py_integrand_info_from_info(info))
    }

    /// Return a detached, read-only snapshot of one integrand's settings.
    ///
    /// Parameters
    /// ----------
    /// process_id : int, optional
    ///     Process containing the integrand. Supply this when selection is ambiguous.
    /// integrand_name : str, optional
    ///     Integrand whose settings are required.
    ///
    /// Returns
    /// -------
    /// SettingsValue
    ///     Serialized settings snapshot supporting ``get(path)``, attribute access,
    ///     indexing, and ``to_dict()``. Mutating it does not update the live session.
    ///
    /// Raises
    /// ------
    /// Exception
    ///     If selection fails, the integrand has not been generated, or its settings
    ///     cannot be serialized.
    ///
    /// Examples
    /// --------
    /// Read a nested setting without changing the live integrand::
    ///
    ///     settings = api.get_integrand_settings(process_id=0)
    ///     print(settings.get("general.generate_events"))
    #[pyo3(name="get_integrand_settings", signature = (process_id=None, integrand_name=None))]
    pub(crate) fn get_integrand_settings(
        &mut self,
        process_id: Option<usize>,
        integrand_name: Option<String>,
    ) -> PyResult<PySettingsValue> {
        let (pid, name) = self
            .gammaloop_state
            .process_list
            .find_integrand(process_id, integrand_name.as_ref())
            .map_err(|e| {
                exceptions::PyException::new_err(format!("Could not find integrand: {}", e))
            })?;
        match &self.gammaloop_state.process_list.processes[pid].collection {
            ProcessCollection::Amplitudes(amplitudes) => {
                match &amplitudes.get(&name).unwrap().integrand {
                    Some(integrand) => PySettingsValue::from_settings(
                        integrand.get_settings(),
                        "integrand settings",
                        "integrand_settings",
                    ),
                    None => Err(exceptions::PyException::new_err(
                        "Integrand for selected amplitude not yet generated",
                    )),
                }
            }
            ProcessCollection::CrossSections(cross_sections) => {
                match &cross_sections.get(&name).unwrap().integrand {
                    Some(integrand) => PySettingsValue::from_settings(
                        integrand.get_settings(),
                        "integrand settings",
                        "integrand_settings",
                    ),
                    None => Err(exceptions::PyException::new_err(
                        "Integrand for selected cross-section not yet generated",
                    )),
                }
            }
        }
    }

    /// Render the current in-memory run history as TOML.
    ///
    /// Returns
    /// -------
    /// str
    ///     TOML representation of the history owned by this API instance.
    ///
    /// Raises
    /// ------
    /// Exception
    ///     If the current history cannot be serialized.
    ///
    /// Notes
    /// -----
    /// This reports the live session, including commands run through ``run``; it does
    /// not imply that the history has been persisted to the state directory.
    ///
    /// Examples
    /// --------
    /// Capture a reproducible run card from the current session::
    ///
    ///     run_card_toml = api.get_run_history()
    #[pyo3(name = "get_run_history", signature = ())]
    pub(crate) fn get_run_history(&self) -> PyResult<String> {
        self.run_history
            .to_toml_string(self.cli_settings.try_strings)
            .map_err(|e| {
                exceptions::PyException::new_err(format!(
                    "Could not render the current run history as TOML: {e}"
                ))
            })
    }

    /// Render the current effective CLI and global settings as TOML.
    ///
    /// Returns
    /// -------
    /// str
    ///     TOML representation of the settings used by this API session.
    ///
    /// Raises
    /// ------
    /// Exception
    ///     If the settings cannot be rendered.
    ///
    /// Examples
    /// --------
    /// Record the effective settings after applying startup overrides::
    ///
    ///     settings_toml = api.get_global_settings()
    #[pyo3(name = "get_global_settings", signature = ())]
    pub(crate) fn get_global_settings(&self) -> PyResult<String> {
        render_smart_toml(&self.cli_settings).map_err(|e| {
            exceptions::PyException::new_err(format!(
                "Could not render the current global settings as TOML: {e}"
            ))
        })
    }

    /// Return the named command blocks in the current run history.
    ///
    /// Returns
    /// -------
    /// dict[str, list[str]]
    ///     Mapping from block name to its rendered CLI commands.
    ///
    /// Examples
    /// --------
    /// Inspect the commands currently grouped under each block::
    ///
    ///     for name, commands in api.get_active_command_blocks().items():
    ///         print(name, commands)
    #[pyo3(name = "get_active_command_blocks", signature = ())]
    pub(crate) fn get_active_command_blocks(
        &self,
    ) -> PyResult<std::collections::BTreeMap<String, Vec<String>>> {
        Ok(self
            .run_history
            .command_blocks
            .iter()
            .map(|block| {
                (
                    block.name.clone(),
                    block.commands.iter().map(display_command).collect(),
                )
            })
            .collect())
    }

    /// Return a detached, read-only snapshot of the default runtime settings.
    ///
    /// Returns
    /// -------
    /// SettingsValue
    ///     Serialized settings including defaults. Use ``get(path)`` or ``to_dict()``
    ///     to inspect it; changes to derived Python values do not affect the session.
    ///
    /// Raises
    /// ------
    /// Exception
    ///     If the settings cannot be serialized.
    ///
    /// Examples
    /// --------
    /// Inspect a runtime default by its documented settings path::
    ///
    ///     runtime = api.get_default_runtime_settings()
    ///     print(runtime.get("integrator.n_start"))
    #[pyo3(name = "get_default_runtime_settings", signature = ())]
    pub(crate) fn get_default_runtime_settings(&self) -> PyResult<PySettingsValue> {
        PySettingsValue::from_settings(
            &self.default_runtime_settings,
            "default runtime settings",
            "default_runtime_settings",
        )
    }

    /// Render a selected amplitude or cross section as Graphviz DOT text.
    ///
    /// Parameters
    /// ----------
    /// process : int or str, optional
    ///     Process id or name; omit only when selection is unambiguous.
    /// integrand_name : str, optional
    ///     Integrand to render.
    /// settings : DotExportSettings, optional
    ///     Controls diagram combination, UV terms, algebra, and generated fields.
    ///
    /// Returns
    /// -------
    /// str
    ///     DOT source suitable for Graphviz or GammaLoop's drawing pipeline.
    #[pyo3(name="get_dot_files", signature = (process=None, integrand_name=None, settings=DotExportSettings::default()))]
    pub(crate) fn get_dot_files(
        &mut self,
        process: Option<ProcessRef>,
        integrand_name: Option<String>,
        settings: DotExportSettings,
    ) -> PyResult<String> {
        let (pid, name) = self
            .gammaloop_state
            .process_list
            .find_integrand_ref(process.as_ref(), integrand_name.as_ref())
            .map_err(|e| {
                exceptions::PyException::new_err(format!("Could not find integrand: {}", e))
            })?;
        let mut dot_output = String::new();
        match &self.gammaloop_state.process_list.processes[pid].collection {
            ProcessCollection::Amplitudes(amplitudes) => match &amplitudes.get(&name) {
                Some(amplitude) => {
                    amplitude
                        .write_dot_fmt(&mut dot_output, &settings)
                        .map_err(|e| {
                            exceptions::PyException::new_err(format!(
                                "Could not write DOT format for amplitude {}: {}",
                                name, e
                            ))
                        })
                }
                None => Err(exceptions::PyException::new_err(format!(
                    "Could not find amplitude named {}",
                    name
                ))),
            },
            ProcessCollection::CrossSections(cross_sections) => match &cross_sections.get(&name) {
                Some(cross_section) => cross_section
                    .write_dot_fmt(&mut dot_output, &settings)
                    .map_err(|e| {
                        exceptions::PyException::new_err(format!(
                            "Could not write DOT format for amplitude {}: {}",
                            name, e
                        ))
                    }),
                None => Err(exceptions::PyException::new_err(format!(
                    "Could not find cross-section named {}",
                    name
                ))),
            },
        }?;
        Ok(dot_output)
    }

    /// Parse and execute a CLI command in this API instance's session.
    ///
    /// Parameters
    /// ----------
    /// command : str
    ///     GammaLoop CLI command text.
    ///
    /// Returns
    /// -------
    /// None
    ///     The command's normal output is emitted through the configured CLI/logging
    ///     sinks; inspect structured state through the corresponding getter methods.
    ///
    /// Raises
    /// ------
    /// Exception
    ///     If the command cannot be parsed or execution fails.
    ///
    /// Notes
    /// -----
    /// Commands share and may mutate the instance's in-memory state, settings, and run
    /// history. ``read_only_state=True`` only blocks writes inside the active state
    /// directory. The API does not automatically persist the session after this call;
    /// persistence depends on the executed command's explicit output behavior.
    ///
    /// Examples
    /// --------
    /// Display the processes loaded in the current state::
    ///
    ///     api.run("display processes")
    #[pyo3(name="run", signature = (command,))]
    pub(crate) fn run_command(&mut self, command: String) -> PyResult<()> {
        let command_history =
            crate::state::CommandHistory::from_raw_string(&command).map_err(|e| {
                exceptions::PyException::new_err(format!(
                    "Failed to parse command '{}': {}",
                    command, e
                ))
            })?;
        let mut session = CliSession::new(
            &mut self.gammaloop_state,
            &mut self.run_history,
            &mut self.cli_settings,
            &mut self.default_runtime_settings,
            &mut self.session_state,
        );

        session
            .execute_command(command_history)
            .map_err(|e| {
                exceptions::PyException::new_err(format!(
                    "Execution of command '{}' failed: {}",
                    command, e
                ))
            })
            .map(|_| ())
    }

    /// Build a causal-flow expression from an inline DOT graph or one of its subgraphs.
    ///
    /// Parameters
    /// ----------
    /// dot_string : str
    ///     Inline DOT graph using particles from the active model.
    /// subgraph_nodes : Sequence[str]
    ///     Vertex names retained in the subgraph; an empty sequence selects all nodes.
    /// reverse_dangling : Sequence[int]
    ///     Dangling edge ids whose orientation is reversed.
    /// orientation_pattern : str, optional
    ///     Pattern restricting returned causal-flow orientations.
    ///
    /// Returns
    /// -------
    /// list[tuple[dict[int, int], str]]
    ///     Edge-direction maps paired with their energy-denominator expressions.
    #[pyo3(name = "generate_cff", signature = (dot_string, subgraph_nodes, reverse_dangling,orientation_pattern=None))]
    pub(crate) fn generate_cff(
        &self,
        dot_string: String,
        subgraph_nodes: Vec<String>,
        reverse_dangling: Vec<usize>,
        orientation_pattern: Option<String>,
    ) -> PyResult<Vec<(HashMap<usize, i32>, String)>> {
        let graph = Graph::from_string(dot_string, &self.gammaloop_state.model)
            .unwrap()
            .pop()
            .unwrap();

        let reverse_dangling = reverse_dangling
            .into_iter()
            .map(EdgeIndex::from)
            .collect_vec();

        let subgraph: SuBitGraph = if subgraph_nodes.is_empty() {
            graph.full_filter()
        } else {
            let mut result: SuBitGraph = graph.empty_subgraph();
            for (_node_id, neighbors, vertex) in graph.iter_nodes() {
                if subgraph_nodes.contains(&vertex.name.to_string()) {
                    neighbors.for_each(|hedge| result.add(hedge));
                }
            }
            result
        };

        let mut surface_cache = SurfaceCache {
            esurface_cache: TiVec::new(),
            hsurface_cache: TiVec::new(),
        };

        let cff = generate_cff_expression_from_subgraph(
            &graph.underlying,
            &subgraph,
            &None,
            &reverse_dangling,
            &graph.get_edges_in_initial_state_cut(),
            &mut surface_cache,
        )
        .map_err(|e| {
            exceptions::PyException::new_err(format!("Could not generate CFF expression: {}", e))
        })?;

        let or_pattern = orientation_pattern
            .as_deref()
            .map(OrientationPattern::from_user_pattern)
            .transpose()
            .map_err(|error| exceptions::PyException::new_err(error.to_string()))?
            .unwrap_or_default();

        let atoms = cff.get_orientation_atoms_with_data(or_pattern);
        let inverse_energies = graph::get_cff_inverse_energy_product_impl(&graph, &subgraph, &[]);

        let result = atoms
            .into_iter()
            .map(|(atom, orientation_data)| {
                let energy_sub = cff.surfaces.substitute_energies(&atom, &[]) * &inverse_energies;
                let string_atom = energy_sub.to_string();

                let mut orientation_as_hashmap = HashMap::new();
                for (edge_id, direction) in orientation_data.orientation.into_iter() {
                    let direction = match direction {
                        Orientation::Default => 1,
                        Orientation::Reversed => -1,
                        Orientation::Undirected => 0,
                    };
                    orientation_as_hashmap.insert(edge_id.0, direction);
                }
                (orientation_as_hashmap, string_atom)
            })
            .collect_vec();

        Ok(result)
    }

    /// Serialize a causal-flow expression and its surfaces as JSON.
    ///
    /// This accepts the same graph, subgraph, and dangling-edge inputs as
    /// ``generate_cff``. The current JSON representation is intended for GammaLoop
    /// tooling and may contain internal structural details.
    ///
    /// Parameters
    /// ----------
    /// dot_string : str
    ///     Inline DOT graph using particles from the active model.
    /// subgraph_nodes : Sequence[str]
    ///     Vertex names retained in the subgraph; an empty sequence selects all nodes.
    /// reverse_dangling : Sequence[int]
    ///     Dangling edge ids whose orientation is reversed.
    /// orientation_pattern : str, optional
    ///     Pattern restricting returned causal-flow orientations.
    ///
    /// Returns
    /// -------
    /// str
    ///     JSON representation of the causal-flow expression and E-surfaces.
    #[pyo3(
        name = "generate_cff_as_json_string",
        signature = (dot_string, subgraph_nodes, reverse_dangling, orientation_pattern = None)
    )]
    pub(crate) fn generate_cff_as_json_string(
        &self,
        dot_string: String,
        subgraph_nodes: Vec<String>,
        reverse_dangling: Vec<usize>,
        orientation_pattern: Option<String>,
    ) -> PyResult<String> {
        let _ = orientation_pattern;
        let graph = Graph::from_string(dot_string, &self.gammaloop_state.model)
            .unwrap()
            .pop()
            .unwrap();

        let reverse_dangling = reverse_dangling
            .into_iter()
            .map(EdgeIndex::from)
            .collect_vec();

        let subgraph: SuBitGraph = if subgraph_nodes.is_empty() {
            graph.full_filter()
        } else {
            let mut result: SuBitGraph = graph.empty_subgraph();
            for (_node_id, neighbors, vertex) in graph.iter_nodes() {
                if subgraph_nodes.contains(&vertex.name.to_string()) {
                    neighbors.for_each(|hedge| result.add(hedge));
                }
            }
            result
        };

        let mut surface_cache = SurfaceCache {
            esurface_cache: TiVec::new(),
            hsurface_cache: TiVec::new(),
        };

        let cff = generate_cff_expression_from_subgraph(
            &graph.underlying,
            &subgraph,
            &None,
            &reverse_dangling,
            &graph.get_edges_in_initial_state_cut(),
            &mut surface_cache,
        )
        .map_err(|e| {
            exceptions::PyException::new_err(format!("Could not generate CFF expression: {}", e))
        })?;

        let json_string = serde_json::to_string(&cff).map_err(|e| {
            exceptions::PyException::new_err(format!("Could not serialize cff to json: {}", e))
        })?;

        Ok(json_string)
    }

    /*

    #[allow(clippy::too_many_arguments)]
    #[pyo3(name="generate", signature = (process_name, numerator_aware_isomorphism_grouping=None, filter_self_loop=None, graph_prefix=None, selected_graphs=None, vetoed_graphs=None, loop_momentum_bases=None, global_prefactor_color=None, global_prefactor_colorless=None, num_threads=None))]
    pub(crate) fn generate_diagrams(
        &mut self,
        //generation_options: PyRef<PyFeynGenOptions>,
        process_name: String,
        numerator_aware_isomorphism_grouping: Option<PyRef<PyNumeratorAwareGroupingOption>>,
        filter_self_loop: Option<bool>,
        graph_prefix: Option<String>,
        selected_graphs: Option<Vec<String>>,
        vetoed_graphs: Option<Vec<String>>,
        loop_momentum_bases: Option<HashMap<String, Vec<String>>>,
        global_prefactor_color: Option<String>,
        global_prefactor_colorless: Option<String>,
        num_threads: Option<usize>,
    ) -> Result<Vec<String>> {
        todo!();
        //if self.model.is_empty() {
        //    return Err(eyre!(
        //        "A physics model must be loaded before generating diagrams"
        //    ));
        //}
        //let feyngen_options = generation_options.options.clone();

        //// clone some of the options that will be used in the process definition
        //let initial_pdgs = feyngen_options.initial_pdgs.clone();
        //let final_pdgs_lists = feyngen_options.final_pdgs_lists.clone();
        //let generation_type = feyngen_options.generation_type.clone();
        //let amplitude_filters = feyngen_options.amplitude_filters.clone();
        //let cross_section_filters = feyngen_options.cross_section_filters.clone();

        //let diagram_generator = FeynGen::new(feyngen_options);

        //let mut global_prefactor = GlobalPrefactor::default();
        //if let Some(global_prefactor_color) = global_prefactor_color {
        //    global_prefactor.num = parse!(&global_prefactor_color);
        //}
        //if let Some(global_prefactor_colorless) = global_prefactor_colorless {
        //    global_prefactor.projector = parse!(&global_prefactor_colorless);
        //}

        //let _diagrams = diagram_generator
        //    .generate(
        //        &self.model,
        //        &numerator_aware_isomorphism_grouping
        //            .map(|o| o.grouping_options.clone())
        //            .unwrap_or(NumeratorAwareGraphGroupingOption::NoGrouping),
        //        filter_self_loop.unwrap_or(false),
        //        graph_prefix.unwrap_or("GL".to_string()),
        //        selected_graphs,
        //        vetoed_graphs,
        //        loop_momentum_bases,
        //        global_prefactor,
        //        num_threads,
        //    )
        //    .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        // let res = Ok(diagrams
        //     .iter()
        //     .map(|d| serde_yaml::to_string(&SerializableGraph::from_graph(d)).unwrap())
        //     .collect());

        // load everything into processlist
        //let (n_unresolved, unresolved_cut_content) =
        //    diagram_generator.unresolved_cut_content(&self.model);

        //let process_definition = ProcessDefinition {
        //    initial_pdgs,
        //    final_pdgs_lists,
        //    n_unresolved,
        //    unresolved_cut_content: unresolved_cut_content.into_iter().collect(),
        //    amplitude_filters,
        //    cross_section_filters,
        //};

        //let process = Process::from_graph_list(
        //    proccess_name.into(),
        //    vec![],
        //    generation_type,
        //    process_definition,
        //    None,
        //)
        //.map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

        //self.process_list.add_process(process);

        //Ok(Vec::new())
    }

    // pub(crate) fn export_expressions(
    //     &mut self,
    //     export_root: &str,
    //     amplitued_list: Vec<String>,
    //     format: &str,
    //     export_yaml_str: &str,
    // ) -> Result<String> {
    //     let export_settings: ProcessSettings = serde_yaml::from_str(export_yaml_str)
    //         .map_err(|e| exceptions::PyException::new_err(e.to_string()))
    //         .unwrap();

    //     for amplitude in amplitued_list.into_iter() {
    //         match Amplitude::from_yaml_str(&self.model, amplitude) {
    //             Ok(amp) => {
    //                 amp.map(|a| a.map(|ag| ag.forget_type()))
    //                     .export_expressions(
    //                         export_root,
    //                         Self::printer_options(format),
    //                         &export_settings,
    //                     )
    //                     .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;
    //             }
    //             Err(e) => return Err(exceptions::PyException::new_err(e.to_string())),
    //         }
    //     }
    //     Ok("Exported expressions".to_string())
    // }

    // pub(crate) fn export_coupling_replacement_rules(
    //     &self,
    //     export_root: &str,
    //     format: &str,
    // ) -> Result<String> {
    //     self.model
    //         .export_coupling_replacement_rules(export_root, Self::printer_options(format))
    //         .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;
    //     Ok("Exported coupling substitutions".to_string())
    // }

    #[allow(clippy::too_many_arguments)]
    #[pyo3(
        name = "integrate",
        signature = (
            slots = None,
            process = None,
            integrand_name = None,
            n_cores = None,
            workspace_path = None,
            target = None,
            restart = false,
            uncorrelated = false,
            show_max_weight_info = true,
            no_show_integration_statistics = false,
            show_phase = "both".to_string(),
            show_top_discrete_grid = false,
            show_discrete_contributions_sum = false,
            sort_contributions = "error".to_string(),
            show_max_weight_info_for_discrete_bins = false,
            show_summary_only = false,
            no_stream_iterations = false,
            no_stream_updates = false,
            batch_size = None,
            batch_timing = 5.0,
            min_time_between_status_updates = 0.0,
            max_table_width = 250,
            write_results_for_each_iteration = false
        )
    )]
    pub fn integrate(
        &mut self,
        py: Python<'_>,
        slots: Option<Vec<(ProcessRef, String)>>,
        process: Option<ProcessRef>,
        integrand_name: Option<String>,
        n_cores: Option<usize>,
        workspace_path: Option<PathBuf>,
        target: Option<PyObject>,
        restart: bool,
        uncorrelated: bool,
        show_max_weight_info: bool,
        no_show_integration_statistics: bool,
        show_phase: String,
        show_top_discrete_grid: bool,
        show_discrete_contributions_sum: bool,
        sort_contributions: String,
        show_max_weight_info_for_discrete_bins: bool,
        show_summary_only: bool,
        no_stream_iterations: bool,
        no_stream_updates: bool,
        batch_size: Option<usize>,
        batch_timing: f64,
        min_time_between_status_updates: f64,
        max_table_width: usize,
        write_results_for_each_iteration: bool,
    ) -> Result<PyIntegrationOutput> {
        let integrate = build_python_integrate_command(
            &self.gammaloop_state,
            slots,
            process,
            integrand_name,
            target.as_ref().map(|target| target.bind(py)),
            n_cores,
            workspace_path,
            restart,
            uncorrelated,
            show_max_weight_info,
            no_show_integration_statistics,
            show_phase,
            show_top_discrete_grid,
            show_discrete_contributions_sum,
            sort_contributions,
            show_max_weight_info_for_discrete_bins,
            show_summary_only,
            no_stream_iterations,
            no_stream_updates,
            batch_size,
            batch_timing,
            min_time_between_status_updates,
            max_table_width,
            write_results_for_each_iteration,
        )?;

        Ok(PyIntegrationOutput {
            inner: integrate.run(&mut self.gammaloop_state, &self.cli_settings)?,
        })
    }

    // pub(crate) fn inspect_lmw_integrand(
    //     &mut self,
    //     integrand_name: &str,
    //     workspace_path: &str,
    //     use_f128: bool,
    // ) -> Result<()> {
    //     let integrand = self.process_list.get_integrand_mut(0, integrand_name)?;

    //     let settings = integrand.get_settings().clone();
    //     let workspace_path = PathBuf::from(workspace_path);
    //     let path_to_state = workspace_path.join("integration_state");

    //     // match fs::read(path_to_state) {
    //     //     Ok(state_bytes) => {
    //     //         let integration_state: IntegrationState =
    //     //             bincode::decode_from_slice::<IntegrationState, _>(
    //     //                 &state_bytes,
    //     //                 bincode::config::standard(),
    //     //             )
    //     //             .expect("failed to obtain state")
    //     //             .0;
    //     //         let path_to_workspace_settings = workspace_path.join("settings.yaml");
    //     //         let workspace_settings_string = fs::read_to_string(path_to_workspace_settings)
    //     //             .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

    //     //         let mut workspace_settings: Settings =
    //     //             serde_yaml::from_str(&workspace_settings_string)
    //     //                 .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

    //     //         workspace_settings.general.debug = new_settings.general.debug;

    //     //         let max_weight_samples = vec![
    //     //             integration_state.integral.re.max_eval_positive_xs,
    //     //             integration_state.integral.re.max_eval_negative_xs,
    //     //             integration_state.integral.im.max_eval_positive_xs,
    //     //             integration_state.integral.im.max_eval_negative_xs,
    //     //         ];

    //     //         for max_weight_sample in max_weight_samples
    //     //             .into_iter()
    //     //             .filter_map(std::convert::identity)
    //     //         {
    //     //             // bypass inspect function as it does not take a symbolica sample as input
    //     //             let eval_result = integrand_struct.evaluate_sample(
    //     //                 &max_weight_sample,
    //     //                 F(0.0),
    //     //                 1,
    //     //                 use_f128,
    //     //                 Complex::new_zero(),
    //     //             );

    //     //             let eval = eval_result.integrand_result;

    //     //             info!(
    //     //                     "\nFor input point xs: \n\n{}\n\nThe evaluation of integrand '{}' is:\n\n{}\n",
    //     //                     format!(
    //     //                         "( {:?} )",
    //     //                         max_weight_sample,
    //     //                     )
    //     //                     .blue(),
    //     //                     integrand_name,
    //     //                     format!("( {:+.16e}, {:+.16e} i)", eval.re, eval.im).blue(),
    //     //                 );
    //     //         }

    //     //         Ok(())
    //     //     }
    //     // }
    //     //
    //     Ok(())
    // }

    // pub(crate) fn load_master_node(&mut self, integrand: &str) -> Result<String> {
    //     let selected_integrand = self.process_list.get_integrand_mut(0, integrand)?;

    //     let grid = selected_integrand.create_grid();
    //     let integrator_settings = selected_integrand.get_integrator_settings();

    //     let master_node = MasterNode::new(grid, integrator_settings);
    //     self.master_node = Some(master_node);

    //     Ok(format!("Initialized master grid for {}", integrand))
    // }

    // pub(crate) fn write_batch_input(
    //     &mut self,
    //     num_cores: usize,
    //     num_samples: usize,
    //     export_grid: bool,
    //     output_accumulator: bool,
    //     workspace_path: &str,
    //     job_id: usize,
    // ) -> Result<String> {
    //     let master_node = self
    //         .master_node
    //         .as_mut()
    //         .expect("Could not get master node");

    //     // extract the integrated phase in a hacky way
    //     match master_node
    //         .write_batch_input(
    //             num_cores,
    //             num_samples,
    //             export_grid,
    //             output_accumulator,
    //             workspace_path,
    //             job_id,
    //         )
    //         .map_err(|e| exceptions::PyException::new_err(e.to_string()))
    //     {
    //         Ok(_) => Ok(format!("Wrote batch input for job {}", job_id)),
    //         Err(e) => Err(e),
    //     }
    // }

    // pub(crate) fn process_batch_output(
    //     &mut self,
    //     workspace_path: &str,
    //     job_id: usize,
    // ) -> Result<String> {
    //     let master_node = self
    //         .master_node
    //         .as_mut()
    //         .expect("could not get master node");

    //     let job_out_name = format!("job_{}_out", job_id);
    //     let job_out_path = Path::new(workspace_path).join(job_out_name);

    //     let output_file = std::fs::read(job_out_path)?;
    //     let batch_result: BatchResult =
    //         bincode::decode_from_slice(&output_file, bincode::config::standard())
    //             .expect("Could not deserialize batch")
    //             .0;

    //     master_node
    //         .process_batch_output(batch_result)
    //         .map_err(|e| exceptions::PyException::new_err(e.to_string()))?;

    //     Ok(format!("Processed job {}", job_id))
    // }

    // pub(crate) fn display_master_node_status(&self) {
    //     if let Some(master_node) = &self.master_node {
    //         master_node.display_status();
    //     }
    // }

    // pub(crate) fn update_iter(&mut self) {
    //     if let Some(master_node) = &mut self.master_node {
    //         master_node.update_iter();
    //     }
    // }


    */
}
