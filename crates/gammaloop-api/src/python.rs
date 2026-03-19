use gammalooprs::{
    cff::generation::{generate_cff_expression_from_subgraph, SurfaceCache},
    feyngen::diagram_generator::evaluate_overall_factor,
    graph::{self, FeynmanGraph, Graph, LMBext},
    integrands::evaluation::{
        BatchSampleEvaluationResult, SampleEvaluationResult, SingleSampleEvaluationResult,
    },
    observables::{
        AdditionalWeightKey, Event, EventGroup, GenericAdditionalWeightInfo, HistogramSnapshot,
        HistogramStatisticsSnapshot,
    },
    processes::{DotExportSettings, ProcessCollection},
    settings::{global::OrientationPattern, RuntimeSettings},
    utils::tracing::{LogFormat, LogLevel},
};
use linnet::half_edge::{
    involution::{EdgeIndex, Orientation},
    subgraph::{ModifySubSet, SuBitGraph},
};
use numpy::PyReadonlyArray2;
use typed_index_collections::TiVec;

use crate::{
    commands::{
        evaluate_samples::EvaluateSamples, import::model::ImportModel, integrate::Integrate,
        Evaluate,
    },
    session::{CliSession, CliSessionState},
    state::{ProcessRef, RunHistory, State},
    CLISettings, LoadedState, StateLoadOption,
};
use ahash::{HashMap, HashMapExt};
use clap::ValueEnum;

use color_eyre::Result;
use eyre::eyre;

/*
use gammalooprs::feyngen::{
    FeynGenError, FeynGenFilter, FeynGenFilters, GenerationType, SelfEnergyFilterOptions,
    SnailFilterOptions, TadpolesFilterOptions,
};
*/

use git_version::git_version;
use itertools::{self, Itertools};
use std::{path::PathBuf, str::FromStr};

use symbolica::{atom::AtomCore, parse};
const GIT_VERSION: &str = git_version!(fallback = "unavailable");

#[allow(unused)]
use pyo3::{
    exceptions,
    prelude::*,
    pyclass,
    pyclass::CompareOp,
    pyfunction, pymethods, pymodule,
    types::{PyComplex, PyComplexMethods, PyDict, PyModule, PyTuple, PyType},
    wrap_pyfunction, FromPyObject, PyRef, Python,
};

#[pyfunction]
#[pyo3(name = "evaluate_graph_overall_factor")]
pub(crate) fn evaluate_graph_overall_factor(overall_factor: &str) -> Result<String> {
    let overall_factor = parse!(overall_factor);
    let overall_factor_evaluated = evaluate_overall_factor(overall_factor.as_view());
    Ok(overall_factor_evaluated.to_canonical_string())
}

#[pyfunction]
#[pyo3(name = "atom_to_canonical_string")]
pub(crate) fn atom_to_canonical_string(atom_str: &str) -> Result<String> {
    Ok(parse!(atom_str).to_canonical_string())
}

#[pymodule(name = "_gammaloop")]
fn python_module(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    gammalooprs::initialisation::initialise().expect("initialization failed");
    gammalooprs::set_interrupt_handler();
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
    m.add_class::<PyAdditionalWeight>()?;
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
    m.add("git_version", GIT_VERSION)?;
    m.add_wrapped(wrap_pyfunction!(atom_to_canonical_string))?;
    m.add_wrapped(wrap_pyfunction!(evaluate_graph_overall_factor))?;
    Ok(())
}

#[pyclass(from_py_object, name = "ComplexValue", get_all)]
#[derive(Clone)]
pub struct PyComplexValue {
    pub re: f64,
    pub im: f64,
}

#[pyclass(from_py_object, name = "FourMomentum", get_all)]
#[derive(Clone)]
pub struct PyFourMomentum {
    pub e: f64,
    pub px: f64,
    pub py: f64,
    pub pz: f64,
}

#[pyclass(from_py_object, name = "AdditionalWeight", get_all)]
#[derive(Clone)]
pub struct PyAdditionalWeight {
    pub key: String,
    pub value: PyComplexValue,
}

#[pyclass(from_py_object, name = "CutInfo", get_all)]
#[derive(Clone)]
pub struct PyCutInfo {
    pub incoming_pdgs: Vec<isize>,
    pub outgoing_pdgs: Vec<isize>,
    pub cut_id: usize,
    pub graph_id: usize,
}

#[pyclass(from_py_object, name = "Event", get_all)]
#[derive(Clone)]
pub struct PyEvent {
    pub incoming_momenta: Vec<PyFourMomentum>,
    pub outgoing_momenta: Vec<PyFourMomentum>,
    pub cut_info: PyCutInfo,
    pub weight: PyComplexValue,
    pub additional_weights: Vec<PyAdditionalWeight>,
}

#[pyclass(from_py_object, name = "EventGroup", get_all)]
#[derive(Clone)]
pub struct PyEventGroup {
    pub events: Vec<PyEvent>,
}

#[pyclass(from_py_object, name = "HistogramBin", get_all)]
#[derive(Clone)]
pub struct PyHistogramBinSnapshot {
    pub x_min: Option<f64>,
    pub x_max: Option<f64>,
    pub entry_count: usize,
    pub sum_weights: f64,
    pub sum_weights_squared: f64,
    pub mitigated_fill_count: usize,
}

#[pyclass(from_py_object, name = "HistogramStats", get_all)]
#[derive(Clone)]
pub struct PyHistogramStatisticsSnapshot {
    pub in_range_entry_count: usize,
    pub nan_value_count: usize,
    pub mitigated_pair_count: usize,
}

#[pyclass(from_py_object, name = "HistogramSnapshot", get_all)]
#[derive(Clone)]
pub struct PyHistogramSnapshot {
    pub title: String,
    pub phase: String,
    pub value_transform: String,
    pub supports_misbinning_mitigation: bool,
    pub x_min: f64,
    pub x_max: f64,
    pub sample_count: usize,
    pub log_x_axis: bool,
    pub log_y_axis: bool,
    pub bins: Vec<PyHistogramBinSnapshot>,
    pub underflow_bin: PyHistogramBinSnapshot,
    pub overflow_bin: PyHistogramBinSnapshot,
    pub statistics: PyHistogramStatisticsSnapshot,
}

#[pyclass(name = "IntegrationOutput", skip_from_py_object)]
#[derive(Clone)]
pub struct PyIntegrationOutput {
    inner: crate::commands::integrate::IntegrationOutput,
}

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

#[pyclass(from_py_object, name = "MaxWeightInfoEntry", get_all)]
#[derive(Clone)]
pub struct PyMaxWeightInfoEntry {
    pub component: String,
    pub sign: String,
    pub max_eval: f64,
    pub coordinates: Option<String>,
}

#[pyclass(from_py_object, name = "DiscreteCoordinate", get_all)]
#[derive(Clone)]
pub struct PyDiscreteCoordinate {
    pub axis_label: String,
    pub bin_index: usize,
    pub bin_label: Option<String>,
}

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

#[pyclass(from_py_object, name = "DiscreteBreakdown", get_all)]
#[derive(Clone)]
pub struct PyDiscreteBreakdown {
    pub axis_label: String,
    pub fixed_coordinates: Vec<PyDiscreteCoordinate>,
    pub entries: Vec<PyDiscreteBreakdownEntry>,
}

#[pyclass(from_py_object, name = "ComponentDiscreteBreakdown", get_all)]
#[derive(Clone)]
pub struct PyComponentDiscreteBreakdown {
    pub re: Option<PyDiscreteBreakdown>,
    pub im: Option<PyDiscreteBreakdown>,
}

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

#[pyclass(from_py_object, name = "IntegrationResult", get_all)]
#[derive(Clone)]
pub struct PyIntegrationResult {
    pub slots: Vec<PySlotIntegrationResult>,
}

#[pymethods]
impl PyHistogramBinSnapshot {
    fn average(&self, sample_count: usize) -> f64 {
        if sample_count == 0 {
            0.0
        } else {
            self.sum_weights / sample_count as f64
        }
    }

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

#[pymethods]
impl PyIntegrationResult {
    fn slot(&self, key: &str) -> Option<PySlotIntegrationResult> {
        self.slots.iter().find(|slot| slot.key == key).cloned()
    }

    fn single_slot(&self) -> Option<PySlotIntegrationResult> {
        (self.slots.len() == 1).then(|| self.slots[0].clone())
    }
}

#[pyclass(from_py_object, name = "StabilityResult", get_all)]
#[derive(Clone)]
pub struct PyStabilityResult {
    pub precision: String,
    pub estimated_relative_accuracy: Option<f64>,
    pub accepted_as_stable: bool,
    pub total_time_seconds: f64,
}

#[pyclass(from_py_object, name = "SampleEvaluationResult")]
#[derive(Clone)]
pub struct PySampleEvaluationResult {
    inner: SampleEvaluationResult,
}

#[pymethods]
impl PySampleEvaluationResult {
    #[getter]
    fn integrand_result<'py>(&self, py: Python<'py>) -> Bound<'py, PyComplex> {
        PyComplex::from_doubles(
            py,
            self.inner.evaluation.integrand_result.re.0,
            self.inner.evaluation.integrand_result.im.0,
        )
    }

    #[getter]
    fn integrator_weight(&self) -> f64 {
        self.inner.evaluation.integrator_weight.0
    }

    #[getter]
    fn generated_event_count(&self) -> Option<usize> {
        self.inner
            .evaluation
            .evaluation_metadata
            .as_ref()
            .map(|metadata| metadata.generated_event_count)
    }

    #[getter]
    fn accepted_event_count(&self) -> Option<usize> {
        self.inner
            .evaluation
            .evaluation_metadata
            .as_ref()
            .map(|metadata| metadata.accepted_event_count)
    }

    #[getter]
    fn event_processing_time_seconds(&self) -> Option<f64> {
        self.inner
            .evaluation
            .evaluation_metadata
            .as_ref()
            .map(|metadata| metadata.event_processing_time.as_secs_f64())
    }

    #[getter]
    fn parameterization_jacobian(&self) -> Option<f64> {
        self.inner
            .evaluation
            .parameterization_jacobian
            .map(|jac| jac.0)
    }

    #[getter]
    fn is_nan(&self) -> Option<bool> {
        self.inner
            .evaluation
            .evaluation_metadata
            .as_ref()
            .map(|metadata| metadata.is_nan)
    }

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
                        accepted_as_stable: result.accepted_as_stable,
                        total_time_seconds: result.total_time.as_secs_f64(),
                    })
                    .collect()
            })
    }

    #[getter]
    fn event_groups(&self) -> Vec<PyEventGroup> {
        self.inner
            .evaluation
            .event_groups
            .iter()
            .map(py_event_group_from_event_group)
            .collect()
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}

#[pyclass(from_py_object, name = "EvaluationResult")]
#[derive(Clone)]
pub struct PyEvaluationResult {
    inner: SingleSampleEvaluationResult,
}

#[pymethods]
impl PyEvaluationResult {
    #[getter]
    fn sample(&self) -> PySampleEvaluationResult {
        PySampleEvaluationResult {
            inner: self.inner.sample.clone(),
        }
    }

    #[getter]
    fn integrand_result<'py>(&self, py: Python<'py>) -> Bound<'py, PyComplex> {
        self.sample().integrand_result(py)
    }

    #[getter]
    fn integrator_weight(&self) -> f64 {
        self.inner.sample.evaluation.integrator_weight.0
    }

    #[getter]
    fn generated_event_count(&self) -> Option<usize> {
        self.sample().generated_event_count()
    }

    #[getter]
    fn accepted_event_count(&self) -> Option<usize> {
        self.sample().accepted_event_count()
    }

    #[getter]
    fn event_processing_time_seconds(&self) -> Option<f64> {
        self.sample().event_processing_time_seconds()
    }

    #[getter]
    fn parameterization_jacobian(&self) -> Option<f64> {
        self.inner
            .sample
            .evaluation
            .parameterization_jacobian
            .map(|jac| jac.0)
    }

    #[getter]
    fn is_nan(&self) -> Option<bool> {
        self.sample().is_nan()
    }

    #[getter]
    fn stability_results(&self) -> Option<Vec<PyStabilityResult>> {
        self.sample().stability_results()
    }

    #[getter]
    fn event_groups(&self) -> Vec<PyEventGroup> {
        self.sample().event_groups()
    }

    #[getter]
    fn observables<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        py_observable_dict_from_bundle(py, &self.inner.observables)
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}

#[pyclass(from_py_object, name = "BatchEvaluationResult")]
#[derive(Clone)]
pub struct PyBatchEvaluationResult {
    inner: BatchSampleEvaluationResult,
}

#[pymethods]
impl PyBatchEvaluationResult {
    #[getter]
    fn samples(&self) -> Vec<PySampleEvaluationResult> {
        self.inner
            .samples
            .iter()
            .cloned()
            .map(|inner| PySampleEvaluationResult { inner })
            .collect()
    }

    #[getter]
    fn observables<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        py_observable_dict_from_bundle(py, &self.inner.observables)
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}

#[pymethods]
impl PyEvent {
    fn __str__(&self) -> String {
        event_from_py_event(self).to_string()
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}

#[pymethods]
impl PyEventGroup {
    fn __str__(&self) -> String {
        event_group_from_py_event_group(self).to_string()
    }

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
                    Some(subset_index) => AdditionalWeightKey::ThresholdCounterterm {
                        subset_index: subset_index.parse().unwrap_or_default(),
                    },
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

fn py_histogram_snapshot_from_snapshot(snapshot: HistogramSnapshot) -> PyHistogramSnapshot {
    PyHistogramSnapshot {
        title: snapshot.title,
        phase: format!("{:?}", snapshot.phase).to_lowercase(),
        value_transform: format!("{:?}", snapshot.value_transform).to_lowercase(),
        supports_misbinning_mitigation: snapshot.supports_misbinning_mitigation,
        x_min: snapshot.x_min,
        x_max: snapshot.x_max,
        sample_count: snapshot.sample_count,
        log_x_axis: snapshot.log_x_axis,
        log_y_axis: snapshot.log_y_axis,
        bins: snapshot
            .bins
            .into_iter()
            .map(|bin| PyHistogramBinSnapshot {
                x_min: bin.x_min,
                x_max: bin.x_max,
                entry_count: bin.entry_count,
                sum_weights: bin.sum_weights,
                sum_weights_squared: bin.sum_weights_squared,
                mitigated_fill_count: bin.mitigated_fill_count,
            })
            .collect(),
        underflow_bin: PyHistogramBinSnapshot {
            x_min: snapshot.underflow_bin.x_min,
            x_max: snapshot.underflow_bin.x_max,
            entry_count: snapshot.underflow_bin.entry_count,
            sum_weights: snapshot.underflow_bin.sum_weights,
            sum_weights_squared: snapshot.underflow_bin.sum_weights_squared,
            mitigated_fill_count: snapshot.underflow_bin.mitigated_fill_count,
        },
        overflow_bin: PyHistogramBinSnapshot {
            x_min: snapshot.overflow_bin.x_min,
            x_max: snapshot.overflow_bin.x_max,
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

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(unsendable, name = "GammaLoopAPI")
)]
struct GammaLoopAPI {
    gammaloop_state: State,
    cli_settings: CLISettings,
    #[allow(unused)]
    run_history: RunHistory,
    default_runtime_settings: RuntimeSettings,
    session_state: CliSessionState,
}

// TODO: Improve error broadcasting to Python everywhere so as to show rust backtrace
#[pymethods]
impl GammaLoopAPI {
    #[new]
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
        fresh_state=false
    ))]
    pub fn new_python(
        state_folder: Option<PathBuf>,
        boot_commands_path: Option<PathBuf>,
        model_file: Option<PathBuf>,
        trace_logs_filename: Option<String>,
        level: Option<LogLevel>,
        logfile_level: Option<LogLevel>,
        logging_prefix: Option<LogFormat>,
        read_only_state: bool,
        settings_global_path: Option<PathBuf>,
        settings_runtime_defaults_path: Option<PathBuf>,
        fresh_state: bool,
    ) -> Result<Self> {
        let LoadedState {
            state,
            run_history,
            cli_settings,
            default_runtime_settings,
            session_state,
            ..
        } = StateLoadOption {
            fresh_state,
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

    #[pyo3(
        name = "evaluate_sample",
        signature = (point, process_id=None, integrand_name=None, use_arb_prec=false, force_radius=false, minimal_output=false, momentum_space=false, discrete_dim=None, graph_name=None, orientation=None)
    )]
    pub fn evaluate_sample<'py>(
        &mut self,
        _py: Python<'py>,
        point: Vec<f64>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        use_arb_prec: bool,
        force_radius: bool,
        minimal_output: bool,
        momentum_space: bool,
        discrete_dim: Option<Vec<usize>>,
        graph_name: Option<String>,
        orientation: Option<usize>,
    ) -> Result<PyEvaluationResult> {
        let points =
            ndarray::Array2::from_shape_vec((1, point.len()), point).map_err(|e| eyre!(e))?;
        let discrete_dims = discrete_dim.unwrap_or_default();
        let discrete_dims = if momentum_space {
            None
        } else {
            Some(
                ndarray::Array2::from_shape_vec((1, discrete_dims.len()), discrete_dims)
                    .map_err(|e| eyre!(e))?,
            )
        };
        let res = EvaluateSamples {
            process_id,
            integrand_name,
            use_arb_prec,
            force_radius,
            minimal_output,
            momentum_space,
            points: points.view(),
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

    #[pyo3(
        name = "evaluate_samples",
        signature = (points, process_id=None, integrand_name=None, use_arb_prec=false, minimal_output=false, momentum_space=false, discrete_dims=None, force_radius=false, graph_names=None, orientations=None)
    )]
    pub fn evaluate_samples<'py>(
        &mut self,
        _py: Python<'py>,
        points: PyReadonlyArray2<'py, f64>,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        use_arb_prec: bool,
        minimal_output: bool,
        momentum_space: bool,
        discrete_dims: Option<PyReadonlyArray2<'py, usize>>,
        force_radius: bool,
        graph_names: Option<Vec<Option<String>>>,
        orientations: Option<Vec<Option<usize>>>,
    ) -> Result<PyBatchEvaluationResult> {
        let points_rust = points.as_array();

        let evaluate_samples = EvaluateSamples {
            process_id,
            integrand_name,
            use_arb_prec,
            force_radius,
            minimal_output,
            momentum_space,
            points: points_rust,
            discrete_dims: discrete_dims.as_ref().map(PyReadonlyArray2::as_array),
            graph_names,
            orientations,
        };

        let res = evaluate_samples.run(&mut self.gammaloop_state)?;
        Ok(PyBatchEvaluationResult { inner: res })
    }

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
                    .cff_expression
                    .as_ref()
                    .unwrap();

                cff.orientation_data
                    .iter()
                    .map(|or_data| or_data.orientation.clone())
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

    #[pyo3(name = "get_model")]
    pub(crate) fn get_model(&self) -> PyResult<String> {
        let serializable_model = self.gammaloop_state.model.to_serializable();
        serde_json::to_string(&serializable_model).map_err(|e| {
            exceptions::PyException::new_err(format!("Could not serialize model: {}", e))
        })
    }

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

    #[pyo3(name="get_integrand_settings", signature = (process_id=None, integrand_name=None))]
    pub(crate) fn get_integrand_settings(
        &mut self,
        process_id: Option<usize>,
        integrand_name: Option<String>,
    ) -> PyResult<RuntimeSettings> {
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
                    Some(integrand) => Ok::<_, PyErr>(integrand.get_settings().clone()),
                    None => Err(exceptions::PyException::new_err(
                        "Integrand for selected amplitude not yet generated",
                    )),
                }
            }
            ProcessCollection::CrossSections(cross_sections) => {
                match &cross_sections.get(&name).unwrap().integrand {
                    Some(integrand) => Ok::<_, PyErr>(integrand.get_settings().clone()),
                    None => Err(exceptions::PyException::new_err(
                        "Integrand for selected cross-section not yet generated",
                    )),
                }
            }
        }
    }

    #[pyo3(name="get_dot_files", signature = (process_id=None, integrand_name=None,settings=DotExportSettings::default()))]
    pub(crate) fn get_dot_files(
        &mut self,
        process_id: Option<usize>,
        integrand_name: Option<String>,
        settings: DotExportSettings,
    ) -> PyResult<String> {
        let (pid, name) = self
            .gammaloop_state
            .process_list
            .find_integrand(process_id, integrand_name.as_ref())
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

        let or_pattern: OrientationPattern = match orientation_pattern {
            Some(pattern) => {
                serde_json::from_str(format!("{{\"pat\": \"{}\"}}", pattern).as_str()).unwrap()
            }
            None => OrientationPattern { pat: None },
        };

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
