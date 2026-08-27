use std::{collections::BTreeMap, path::PathBuf, sync::Arc};

use feynkit_model::{
    ComplexValue, Coupling, EvaluatedValues, EvaluationRequest, LorentzStructure, Model,
    ModelEvaluator, ModelExpression, ModelFormFactor, ModelFunction, Parameter, ParameterCard,
    ParameterNature, ParameterType, Particle, Propagator, VertexRule,
};
use pyo3::{
    prelude::*,
    types::{PyAny, PyModule},
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

use crate::{display::escape_html, error};

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Particle",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyParticle {
    inner: Particle,
}

impl From<Particle> for PyParticle {
    fn from(inner: Particle) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyParticle {
    /// Return the particle name used by the model.
    ///
    /// Examples
    /// --------
    /// >>> particle.name
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Return the name of the corresponding antiparticle.
    ///
    /// Examples
    /// --------
    /// >>> particle.antiname
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn antiname(&self) -> &str {
        &self.inner.antiname
    }

    /// Return the particle's PDG code.
    ///
    /// Examples
    /// --------
    /// >>> particle.pdg_code
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn pdg_code(&self) -> i64 {
        self.inner.pdg_code
    }

    /// Return the UFO spin representation code.
    ///
    /// Examples
    /// --------
    /// >>> particle.spin
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn spin(&self) -> i64 {
        self.inner.spin
    }

    /// Return the UFO color representation code.
    ///
    /// Examples
    /// --------
    /// >>> particle.color
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn color(&self) -> i64 {
        self.inner.color
    }

    /// Return the name of the particle's mass parameter.
    ///
    /// Examples
    /// --------
    /// >>> particle.mass_parameter
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn mass_parameter(&self) -> &str {
        &self.inner.mass
    }

    /// Return the name of the particle's width parameter.
    ///
    /// Examples
    /// --------
    /// >>> particle.width_parameter
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn width_parameter(&self) -> &str {
        &self.inner.width
    }

    /// Return the particle's electric charge.
    ///
    /// Examples
    /// --------
    /// >>> particle.charge
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn charge(&self) -> f64 {
        self.inner.charge
    }

    /// Report whether this object represents an antiparticle.
    ///
    /// Examples
    /// --------
    /// >>> particle.is_antiparticle
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn is_antiparticle(&self) -> bool {
        self.inner.is_antiparticle()
    }

    /// Report whether the particle is its own antiparticle.
    ///
    /// Examples
    /// --------
    /// >>> particle.is_self_antiparticle
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn is_self_antiparticle(&self) -> bool {
        self.inner.is_self_antiparticle()
    }

    /// Report whether the particle has fermionic spin.
    ///
    /// Examples
    /// --------
    /// >>> particle.is_fermion
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn is_fermion(&self) -> bool {
        self.inner.is_fermion()
    }

    /// Report whether the particle's mass parameter is zero.
    ///
    /// Examples
    /// --------
    /// >>> particle.is_massless
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn is_massless(&self) -> bool {
        self.inner.is_massless()
    }

    /// Return a concise representation containing the name and PDG code.
    ///
    /// Examples
    /// --------
    /// >>> repr(particle)
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        format!(
            "Particle(name='{}', pdg_code={})",
            self.inner.name, self.inner.pdg_code
        )
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    name = "ParameterNature",
    module = "symbolica.community.feynkit",
    rename_all = "SCREAMING_SNAKE_CASE",
    eq,
    eq_int,
    from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PyParameterNature {
    External,
    Internal,
}

impl From<ParameterNature> for PyParameterNature {
    fn from(value: ParameterNature) -> Self {
        match value {
            ParameterNature::External => Self::External,
            ParameterNature::Internal => Self::Internal,
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass_enum)]
#[pyclass(
    name = "ParameterType",
    module = "symbolica.community.feynkit",
    rename_all = "SCREAMING_SNAKE_CASE",
    eq,
    eq_int,
    from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PyParameterType {
    Real,
    Complex,
}

impl From<ParameterType> for PyParameterType {
    fn from(value: ParameterType) -> Self {
        match value {
            ParameterType::Real => Self::Real,
            ParameterType::Complex => Self::Complex,
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Parameter",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyParameter {
    inner: Parameter,
}

impl From<Parameter> for PyParameter {
    fn from(inner: Parameter) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyParameter {
    /// Return the parameter name.
    ///
    /// Examples
    /// --------
    /// >>> parameter.name
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Return the Les Houches block name, when defined.
    ///
    /// Examples
    /// --------
    /// >>> parameter.lhablock
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn lhablock(&self) -> Option<String> {
        self.inner.lhablock.clone()
    }

    /// Return the Les Houches entry indices, when defined.
    ///
    /// Examples
    /// --------
    /// >>> parameter.lhacode
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn lhacode(&self) -> Option<Vec<usize>> {
        self.inner.lhacode.clone()
    }

    /// Return whether the parameter is external or internal.
    ///
    /// Examples
    /// --------
    /// >>> parameter.nature
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn nature(&self) -> PyParameterNature {
        self.inner.nature.clone().into()
    }

    /// Return whether the parameter is real or complex.
    ///
    /// Examples
    /// --------
    /// >>> parameter.parameter_type
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn parameter_type(&self) -> PyParameterType {
        self.inner.parameter_type.clone().into()
    }

    /// Return the evaluated complex value, when available.
    ///
    /// Examples
    /// --------
    /// >>> parameter.value
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn value(&self) -> Option<(f64, f64)> {
        self.inner.value.map(Into::into)
    }

    /// Return the defining expression for an internal parameter.
    ///
    /// Examples
    /// --------
    /// >>> parameter.expression
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn expression(&self) -> Option<String> {
        self.inner.expression.clone()
    }

    /// Return a concise representation containing the parameter name.
    ///
    /// Examples
    /// --------
    /// >>> repr(parameter)
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        format!("Parameter(name='{}')", self.inner.name)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Coupling",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyCoupling {
    inner: Coupling,
}

impl From<Coupling> for PyCoupling {
    fn from(inner: Coupling) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyCoupling {
    /// Return the coupling name.
    ///
    /// Examples
    /// --------
    /// >>> coupling.name
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Return the expression defining the coupling.
    ///
    /// Examples
    /// --------
    /// >>> coupling.expression
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn expression(&self) -> &str {
        &self.inner.expression
    }

    /// Return the coupling-order powers keyed by order name.
    ///
    /// Examples
    /// --------
    /// >>> coupling.orders
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn orders(&self) -> BTreeMap<String, usize> {
        self.inner.orders.clone()
    }

    /// Return the evaluated complex coupling value, when available.
    ///
    /// Examples
    /// --------
    /// >>> coupling.value
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn value(&self) -> Option<(f64, f64)> {
        self.inner.value.map(Into::into)
    }

    /// Return a concise representation containing the coupling name.
    ///
    /// Examples
    /// --------
    /// >>> repr(coupling)
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        format!("Coupling(name='{}')", self.inner.name)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "VertexRule",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyVertexRule {
    inner: VertexRule,
}

impl From<VertexRule> for PyVertexRule {
    fn from(inner: VertexRule) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyVertexRule {
    /// Return the vertex-rule name.
    ///
    /// Examples
    /// --------
    /// >>> vertex.name
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Return the ordered particle names attached to the vertex.
    ///
    /// Examples
    /// --------
    /// >>> vertex.particles
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn particles(&self) -> Vec<String> {
        self.inner.particles.clone()
    }

    /// Return the color structures used by the vertex.
    ///
    /// Examples
    /// --------
    /// >>> vertex.color_structures
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn color_structures(&self) -> Vec<String> {
        self.inner.color_structures.clone()
    }

    /// Return the Lorentz-structure names used by the vertex.
    ///
    /// Examples
    /// --------
    /// >>> vertex.lorentz_structures
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn lorentz_structures(&self) -> Vec<String> {
        self.inner.lorentz_structures.clone()
    }

    /// Return the coupling-name matrix indexed by color and Lorentz structure.
    ///
    /// Examples
    /// --------
    /// >>> vertex.couplings
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn couplings(&self) -> Vec<Vec<Option<String>>> {
        self.inner.couplings.clone()
    }

    /// Combine the coupling-order powers referenced by this vertex.
    ///
    /// Examples
    /// --------
    /// >>> vertex.coupling_orders(model)
    ///
    /// Parameters
    /// ----------
    /// model : Model
    ///     Model containing the referenced couplings.
    fn coupling_orders(&self, model: &PyModel) -> PyResult<BTreeMap<String, usize>> {
        self.inner
            .coupling_orders(&model.inner)
            .map_err(error::model)
    }

    /// Return a concise representation containing the vertex-rule name.
    ///
    /// Examples
    /// --------
    /// >>> repr(vertex)
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        format!("VertexRule(name='{}')", self.inner.name)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "LorentzStructure",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyLorentzStructure {
    inner: LorentzStructure,
}

impl From<LorentzStructure> for PyLorentzStructure {
    fn from(inner: LorentzStructure) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyLorentzStructure {
    /// Return the Lorentz-structure name.
    ///
    /// Examples
    /// --------
    /// >>> lorentz.name
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Return the spin codes of the structure's particle slots.
    ///
    /// Examples
    /// --------
    /// >>> lorentz.spins
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn spins(&self) -> Vec<i64> {
        self.inner.spins.clone()
    }

    /// Return the symbolic Lorentz expression.
    ///
    /// Examples
    /// --------
    /// >>> lorentz.structure
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn structure(&self) -> &str {
        &self.inner.structure
    }

    /// Return a concise representation containing the structure name.
    ///
    /// Examples
    /// --------
    /// >>> repr(lorentz)
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        format!("LorentzStructure(name='{}')", self.inner.name)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Propagator",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyPropagator {
    inner: Propagator,
}

impl From<Propagator> for PyPropagator {
    fn from(inner: Propagator) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyPropagator {
    /// Return the propagator name.
    ///
    /// Examples
    /// --------
    /// >>> propagator.name
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Return the name of the particle using this propagator.
    ///
    /// Examples
    /// --------
    /// >>> propagator.particle
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn particle(&self) -> &str {
        &self.inner.particle
    }

    /// Return the symbolic propagator numerator.
    ///
    /// Examples
    /// --------
    /// >>> propagator.numerator
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn numerator(&self) -> &str {
        &self.inner.numerator
    }

    /// Return the symbolic propagator denominator.
    ///
    /// Examples
    /// --------
    /// >>> propagator.denominator
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn denominator(&self) -> &str {
        &self.inner.denominator
    }

    /// Return a concise representation containing the propagator name.
    ///
    /// Examples
    /// --------
    /// >>> repr(propagator)
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        format!("Propagator(name='{}')", self.inner.name)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ModelFunction",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyModelFunction {
    inner: ModelFunction,
}

impl From<ModelFunction> for PyModelFunction {
    fn from(inner: ModelFunction) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyModelFunction {
    /// Return the function name.
    ///
    /// Examples
    /// --------
    /// >>> function.name
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Return the ordered argument names.
    ///
    /// Examples
    /// --------
    /// >>> function.arguments
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn arguments(&self) -> Vec<String> {
        self.inner.arguments.clone()
    }

    /// Return the function body, when one is defined by the model.
    ///
    /// Examples
    /// --------
    /// >>> function.expression
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn expression(&self) -> Option<String> {
        self.inner.expression.clone()
    }

    /// Return a concise representation containing the function name.
    ///
    /// Examples
    /// --------
    /// >>> repr(function)
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        format!("ModelFunction(name='{}')", self.inner.name)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "FormFactor",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyFormFactor {
    inner: ModelFormFactor,
}

impl From<ModelFormFactor> for PyFormFactor {
    fn from(inner: ModelFormFactor) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyFormFactor {
    /// Return the form-factor name.
    ///
    /// Examples
    /// --------
    /// >>> form_factor.name
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Return the model-defined form-factor type, when present.
    ///
    /// Examples
    /// --------
    /// >>> form_factor.type_name
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn type_name(&self) -> Option<String> {
        self.inner.type_name.clone()
    }

    /// Return the symbolic form-factor value, when present.
    ///
    /// Examples
    /// --------
    /// >>> form_factor.value
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn value(&self) -> Option<String> {
        self.inner.value.clone()
    }

    /// Return a concise representation containing the form-factor name.
    ///
    /// Examples
    /// --------
    /// >>> repr(form_factor)
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        format!("FormFactor(name='{}')", self.inner.name)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ModelExpression",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyModelExpression {
    inner: ModelExpression,
}

impl From<ModelExpression> for PyModelExpression {
    fn from(inner: ModelExpression) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyModelExpression {
    /// Return the name assigned to the expression.
    ///
    /// Examples
    /// --------
    /// >>> expression.name
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Return the symbolic expression text.
    ///
    /// Examples
    /// --------
    /// >>> expression.expression
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn expression(&self) -> &str {
        &self.inner.expression
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "EvaluationRequest",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyEvaluationRequest {
    inner: EvaluationRequest,
}

impl From<EvaluationRequest> for PyEvaluationRequest {
    fn from(inner: EvaluationRequest) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyEvaluationRequest {
    /// Return the known parameter values keyed by name.
    ///
    /// Examples
    /// --------
    /// >>> request.known_parameters
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn known_parameters(&self) -> BTreeMap<String, (f64, f64)> {
        self.inner
            .known_parameters
            .iter()
            .map(|(name, value)| (name.clone(), (*value).into()))
            .collect()
    }

    /// Return the internal expressions that must be evaluated.
    ///
    /// Examples
    /// --------
    /// >>> request.internal_parameters
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn internal_parameters(&self) -> Vec<PyModelExpression> {
        self.inner
            .internal_parameters
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Return the coupling expressions that must be evaluated.
    ///
    /// Examples
    /// --------
    /// >>> request.couplings
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn couplings(&self) -> Vec<PyModelExpression> {
        self.inner
            .couplings
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Return the function definitions available during evaluation.
    ///
    /// Examples
    /// --------
    /// >>> request.functions
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn functions(&self) -> Vec<PyModelFunction> {
        self.inner
            .functions
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Return the form-factor definitions available during evaluation.
    ///
    /// Examples
    /// --------
    /// >>> request.form_factors
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn form_factors(&self) -> Vec<PyFormFactor> {
        self.inner
            .form_factors
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "EvaluatedValues",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyEvaluatedValues {
    inner: EvaluatedValues,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyEvaluatedValues {
    /// Create values returned by a custom model evaluator.
    ///
    /// Examples
    /// --------
    /// >>> values = EvaluatedValues(couplings={"GC_1": (1.0, 0.0)})
    ///
    /// Parameters
    /// ----------
    /// internal_parameters : dict[str, tuple[float, float]] or None
    ///     Evaluated internal parameters keyed by name.
    /// couplings : dict[str, tuple[float, float]] or None
    ///     Evaluated couplings keyed by name.
    #[new]
    #[pyo3(signature = (*, internal_parameters=None, couplings=None))]
    fn new(
        internal_parameters: Option<BTreeMap<String, (f64, f64)>>,
        couplings: Option<BTreeMap<String, (f64, f64)>>,
    ) -> Self {
        Self {
            inner: EvaluatedValues {
                internal_parameters: internal_parameters
                    .unwrap_or_default()
                    .into_iter()
                    .map(|(name, value)| (name, value.into()))
                    .collect(),
                couplings: couplings
                    .unwrap_or_default()
                    .into_iter()
                    .map(|(name, value)| (name, value.into()))
                    .collect(),
            },
        }
    }

    /// Return the evaluated internal parameters keyed by name.
    ///
    /// Examples
    /// --------
    /// >>> values.internal_parameters
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn internal_parameters(&self) -> BTreeMap<String, (f64, f64)> {
        self.inner
            .internal_parameters
            .iter()
            .map(|(name, value)| (name.clone(), (*value).into()))
            .collect()
    }

    /// Return the evaluated couplings keyed by name.
    ///
    /// Examples
    /// --------
    /// >>> values.couplings
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn couplings(&self) -> BTreeMap<String, (f64, f64)> {
        self.inner
            .couplings
            .iter()
            .map(|(name, value)| (name.clone(), (*value).into()))
            .collect()
    }
}

struct PythonModelEvaluator<'a, 'py> {
    callable: &'a Bound<'py, PyAny>,
}

impl ModelEvaluator for PythonModelEvaluator<'_, '_> {
    type Error = PyErr;

    fn evaluate(&mut self, request: EvaluationRequest) -> Result<EvaluatedValues, Self::Error> {
        Ok(self
            .callable
            .call1((PyEvaluationRequest::from(request),))?
            .extract::<PyEvaluatedValues>()
            .map(|values| values.inner)?)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "ParameterCard",
    module = "symbolica.community.feynkit",
    from_py_object
)]
#[derive(Clone, Default)]
pub struct PyParameterCard {
    pub(crate) inner: ParameterCard,
}

impl From<ParameterCard> for PyParameterCard {
    fn from(inner: ParameterCard) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyParameterCard {
    /// Create an empty parameter card.
    ///
    /// Examples
    /// --------
    /// >>> card = ParameterCard()
    ///
    /// Parameters
    /// ----------
    /// None
    #[new]
    fn new() -> Self {
        Self::default()
    }

    /// Parse a parameter card from its JSON representation.
    ///
    /// Examples
    /// --------
    /// >>> card = ParameterCard.from_json('{"mass": [1.0, 0.0]}')
    ///
    /// Parameters
    /// ----------
    /// json : str
    ///     Serialized parameter-card object.
    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        ParameterCard::from_json(json)
            .map(Self::from)
            .map_err(error::model)
    }

    /// Load a parameter card from a JSON file.
    ///
    /// Examples
    /// --------
    /// >>> card = ParameterCard.from_path("parameters.json")
    ///
    /// Parameters
    /// ----------
    /// path : str or os.PathLike
    ///     Path to the JSON parameter card.
    #[staticmethod]
    fn from_path(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        py.detach(move || ParameterCard::from_path(path))
            .map(Self::from)
            .map_err(error::model)
    }

    /// Return a parameter's complex value, or ``None`` when it is absent.
    ///
    /// Examples
    /// --------
    /// >>> card.get("mass")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Parameter name to look up.
    fn get(&self, name: &str) -> Option<(f64, f64)> {
        self.inner.get(name).copied().map(Into::into)
    }

    /// Store a complex parameter value.
    ///
    /// Examples
    /// --------
    /// >>> card.set("mass", 1.0)
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Parameter name to store.
    /// real : float
    ///     Real component of the value.
    /// imaginary : float
    ///     Imaginary component of the value.
    #[pyo3(signature = (name, real, imaginary=0.0))]
    fn set(&mut self, name: String, real: f64, imaginary: f64) {
        self.inner.insert(name, ComplexValue::new(real, imaginary));
    }

    /// Remove and return a parameter value, if present.
    ///
    /// Examples
    /// --------
    /// >>> card.remove("mass")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Parameter name to remove.
    fn remove(&mut self, name: &str) -> Option<(f64, f64)> {
        self.inner.remove(name).map(Into::into)
    }

    /// Return all parameter-card entries.
    ///
    /// Examples
    /// --------
    /// >>> card.items()
    ///
    /// Parameters
    /// ----------
    /// None
    fn items(&self) -> Vec<(String, (f64, f64))> {
        self.inner
            .iter()
            .map(|(name, value)| (name.clone(), (*value).into()))
            .collect()
    }

    /// Serialize the parameter card as JSON.
    ///
    /// Examples
    /// --------
    /// >>> card.to_json(pretty=False)
    ///
    /// Parameters
    /// ----------
    /// pretty : bool
    ///     Indent the output when true.
    #[pyo3(signature = (pretty=true))]
    fn to_json(&self, pretty: bool) -> PyResult<String> {
        if pretty {
            self.inner.to_json_pretty()
        } else {
            self.inner.to_json()
        }
        .map_err(error::model)
    }

    /// Write the parameter card to a JSON file.
    ///
    /// Examples
    /// --------
    /// >>> card.write_json("parameters.json")
    ///
    /// Parameters
    /// ----------
    /// path : str or os.PathLike
    ///     Destination for the JSON parameter card.
    fn write_json(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        let card = self.inner.clone();
        py.detach(move || card.write_json(path))
            .map_err(error::model)
    }

    /// Return the number of entries in the parameter card.
    ///
    /// Examples
    /// --------
    /// >>> len(card)
    ///
    /// Parameters
    /// ----------
    /// None
    fn __len__(&self) -> usize {
        self.inner.len()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "Model",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyModel {
    pub(crate) inner: Arc<Model>,
}

impl From<Model> for PyModel {
    fn from(inner: Model) -> Self {
        Self {
            inner: Arc::new(inner),
        }
    }
}

impl From<Arc<Model>> for PyModel {
    fn from(inner: Arc<Model>) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyModel {
    /// Parse a model from its JSON representation.
    ///
    /// Examples
    /// --------
    /// >>> model = Model.from_json(model_json)
    ///
    /// Parameters
    /// ----------
    /// json : str
    ///     Serialized model object.
    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        Model::from_json(json).map(Self::from).map_err(error::model)
    }

    /// Load a model from a JSON file.
    ///
    /// Examples
    /// --------
    /// >>> model = Model.from_path("model.json")
    ///
    /// Parameters
    /// ----------
    /// path : str or os.PathLike
    ///     Path to the JSON model.
    #[staticmethod]
    fn from_path(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        py.detach(move || Model::from_path(path))
            .map(Self::from)
            .map_err(error::model)
    }

    /// Return the model name.
    ///
    /// Examples
    /// --------
    /// >>> model.name
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    /// Return the applied restriction name, when present.
    ///
    /// Examples
    /// --------
    /// >>> model.restriction
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn restriction(&self) -> Option<String> {
        self.inner.restriction().map(str::to_owned)
    }

    /// Return all particles in model order.
    ///
    /// Examples
    /// --------
    /// >>> model.particles
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn particles(&self) -> Vec<PyParticle> {
        self.inner
            .particles()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Look up a particle by name.
    ///
    /// Examples
    /// --------
    /// >>> particle = model.particle("electron")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Particle name or antiname.
    fn particle(&self, name: &str) -> PyResult<PyParticle> {
        self.inner
            .particle(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    /// Look up a particle by PDG code.
    ///
    /// Examples
    /// --------
    /// >>> particle = model.particle_by_pdg(11)
    ///
    /// Parameters
    /// ----------
    /// pdg : int
    ///     Signed PDG particle code.
    fn particle_by_pdg(&self, pdg: i64) -> PyResult<PyParticle> {
        self.inner
            .particle_by_pdg(pdg)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    /// Return all parameters in model order.
    ///
    /// Examples
    /// --------
    /// >>> model.parameters
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn parameters(&self) -> Vec<PyParameter> {
        self.inner
            .parameters()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Look up a parameter by name.
    ///
    /// Examples
    /// --------
    /// >>> parameter = model.parameter("mass")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Parameter name.
    fn parameter(&self, name: &str) -> PyResult<PyParameter> {
        self.inner
            .parameter(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    /// Return all couplings in model order.
    ///
    /// Examples
    /// --------
    /// >>> model.couplings
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn couplings(&self) -> Vec<PyCoupling> {
        self.inner
            .couplings()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Look up a coupling by name.
    ///
    /// Examples
    /// --------
    /// >>> coupling = model.coupling("GC_1")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Coupling name.
    fn coupling(&self, name: &str) -> PyResult<PyCoupling> {
        self.inner
            .coupling(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    /// Return all vertex rules in model order.
    ///
    /// Examples
    /// --------
    /// >>> model.vertex_rules
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn vertex_rules(&self) -> Vec<PyVertexRule> {
        self.inner
            .vertex_rules()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Look up a vertex rule by name.
    ///
    /// Examples
    /// --------
    /// >>> vertex = model.vertex_rule("V_1")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Vertex-rule name.
    fn vertex_rule(&self, name: &str) -> PyResult<PyVertexRule> {
        self.inner
            .vertex_rule(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    /// Return all Lorentz structures in model order.
    ///
    /// Examples
    /// --------
    /// >>> model.lorentz_structures
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn lorentz_structures(&self) -> Vec<PyLorentzStructure> {
        self.inner
            .lorentz_structures()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Look up a Lorentz structure by name.
    ///
    /// Examples
    /// --------
    /// >>> lorentz = model.lorentz_structure("FFV1")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Lorentz-structure name.
    fn lorentz_structure(&self, name: &str) -> PyResult<PyLorentzStructure> {
        self.inner
            .lorentz_structure(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    /// Return all propagators in model order.
    ///
    /// Examples
    /// --------
    /// >>> model.propagators
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn propagators(&self) -> Vec<PyPropagator> {
        self.inner
            .propagators()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Look up a propagator by name.
    ///
    /// Examples
    /// --------
    /// >>> propagator = model.propagator("electron")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Propagator name.
    fn propagator(&self, name: &str) -> PyResult<PyPropagator> {
        self.inner
            .propagator(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    /// Return all model functions in model order.
    ///
    /// Examples
    /// --------
    /// >>> model.functions
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn functions(&self) -> Vec<PyModelFunction> {
        self.inner
            .functions()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Look up a model function by name.
    ///
    /// Examples
    /// --------
    /// >>> function = model.function("complexconjugate")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Function name.
    fn function(&self, name: &str) -> PyResult<PyModelFunction> {
        self.inner
            .function(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    /// Return all form factors in model order.
    ///
    /// Examples
    /// --------
    /// >>> model.form_factors
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn form_factors(&self) -> Vec<PyFormFactor> {
        self.inner
            .form_factors()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    /// Look up a form factor by name.
    ///
    /// Examples
    /// --------
    /// >>> form_factor = model.form_factor("FF_1")
    ///
    /// Parameters
    /// ----------
    /// name : str
    ///     Form-factor name.
    fn form_factor(&self, name: &str) -> PyResult<PyFormFactor> {
        self.inner
            .form_factor(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    /// Build a parameter card from the model's current external values.
    ///
    /// Examples
    /// --------
    /// >>> card = model.default_parameter_card()
    ///
    /// Parameters
    /// ----------
    /// None
    fn default_parameter_card(&self) -> PyResult<PyParameterCard> {
        self.inner
            .default_parameter_card()
            .map(Into::into)
            .map_err(error::model)
    }

    /// Return a copy with a parameter card applied atomically.
    ///
    /// Examples
    /// --------
    /// >>> updated = model.with_parameter_card(card)
    ///
    /// Parameters
    /// ----------
    /// card : ParameterCard
    ///     External parameter values to apply.
    /// evaluator : Callable[[EvaluationRequest], EvaluatedValues] or None
    ///     Optional evaluator used to recompute dependent values.
    #[pyo3(signature = (card, evaluator=None))]
    fn with_parameter_card(
        &self,
        card: &PyParameterCard,
        #[gen_stub(override_type(
            type_repr = "collections.abc.Callable[[EvaluationRequest], EvaluatedValues] | None",
            imports = ("collections.abc")
        ))]
        evaluator: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Self> {
        let mut model = self.inner.as_ref().clone();
        if let Some(callable) = evaluator {
            let mut evaluator = PythonModelEvaluator { callable };
            model
                .apply_parameter_card_with(&card.inner, &mut evaluator)
                .map_err(error::recompute)?;
        } else {
            model
                .apply_parameter_card(&card.inner)
                .map_err(error::model)?;
        }
        Ok(model.into())
    }

    /// Return a copy with all dependent values recomputed by a callback.
    ///
    /// Examples
    /// --------
    /// >>> recomputed = model.recompute_with(evaluate)
    ///
    /// Parameters
    /// ----------
    /// evaluator : Callable[[EvaluationRequest], EvaluatedValues]
    ///     Callback that evaluates every expression in a request.
    fn recompute_with(
        &self,
        #[gen_stub(override_type(
            type_repr = "collections.abc.Callable[[EvaluationRequest], EvaluatedValues]",
            imports = ("collections.abc")
        ))]
        evaluator: &Bound<'_, PyAny>,
    ) -> PyResult<Self> {
        let mut model = self.inner.as_ref().clone();
        let mut evaluator = PythonModelEvaluator {
            callable: evaluator,
        };
        model
            .recompute_with(&mut evaluator)
            .map_err(error::recompute)?;
        Ok(model.into())
    }

    /// Serialize the model as JSON.
    ///
    /// Examples
    /// --------
    /// >>> model.to_json(pretty=False)
    ///
    /// Parameters
    /// ----------
    /// pretty : bool
    ///     Indent the output when true.
    #[pyo3(signature = (pretty=true))]
    fn to_json(&self, pretty: bool) -> PyResult<String> {
        if pretty {
            self.inner.to_json_pretty()
        } else {
            self.inner.to_json()
        }
        .map_err(error::model)
    }

    /// Write the model to a JSON file.
    ///
    /// Examples
    /// --------
    /// >>> model.write_json("model.json")
    ///
    /// Parameters
    /// ----------
    /// path : str or os.PathLike
    ///     Destination for the JSON model.
    fn write_json(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        let model = self.inner.clone();
        py.detach(move || model.write_json(path))
            .map_err(error::model)
    }

    /// Return a concise representation containing the name and particle count.
    ///
    /// Examples
    /// --------
    /// >>> repr(model)
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        format!(
            "Model(name='{}', particles={})",
            self.inner.name(),
            self.inner.particles().len()
        )
    }

    /// Render a compact inventory of the model in notebook frontends.
    ///
    /// Examples
    /// --------
    /// Leave ``model`` as the final expression in a notebook cell to display
    /// its particle-content and interaction counts.
    ///
    /// Parameters
    /// ----------
    /// None
    fn _repr_html_(&self) -> String {
        let restriction = self
            .inner
            .restriction()
            .map_or_else(|| "none".to_owned(), escape_html);
        format!(
            "<div class=\"feynkit-model\" style=\"display:inline-block;max-width:100%;\
             overflow-x:auto\"><strong>{}</strong><span style=\"margin-left:.45rem;\
             opacity:.7\">restriction: {restriction}</span><table style=\"border-collapse:\
             collapse;margin-top:.35rem\"><thead><tr><th style=\"padding:.2rem .65rem;\
             text-align:right\">particles</th><th style=\"padding:.2rem .65rem;text-align:\
             right\">parameters</th><th style=\"padding:.2rem .65rem;text-align:right\">\
             couplings</th><th style=\"padding:.2rem .65rem;text-align:right\">vertex rules\
             </th><th style=\"padding:.2rem .65rem;text-align:right\">Lorentz structures\
             </th></tr></thead><tbody><tr><td style=\"padding:.2rem .65rem;text-align:right\">\
             {}</td><td style=\"padding:.2rem .65rem;text-align:right\">{}</td><td style=\"\
             padding:.2rem .65rem;text-align:right\">{}</td><td style=\"padding:.2rem .65rem;\
             text-align:right\">{}</td><td style=\"padding:.2rem .65rem;text-align:right\">\
             {}</td></tr></tbody></table></div>",
            escape_html(self.inner.name()),
            self.inner.particles().len(),
            self.inner.parameters().len(),
            self.inner.couplings().len(),
            self.inner.vertex_rules().len(),
            self.inner.lorentz_structures().len(),
        )
    }

    /// Write the concise model summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython invokes this method when only a text representation is supported.
    ///
    /// Parameters
    /// ----------
    /// pretty
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

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyParticle>()?;
    module.add_class::<PyParameterNature>()?;
    module.add_class::<PyParameterType>()?;
    module.add_class::<PyParameter>()?;
    module.add_class::<PyCoupling>()?;
    module.add_class::<PyVertexRule>()?;
    module.add_class::<PyLorentzStructure>()?;
    module.add_class::<PyPropagator>()?;
    module.add_class::<PyModelFunction>()?;
    module.add_class::<PyFormFactor>()?;
    module.add_class::<PyModelExpression>()?;
    module.add_class::<PyEvaluationRequest>()?;
    module.add_class::<PyEvaluatedValues>()?;
    module.add_class::<PyParameterCard>()?;
    module.add_class::<PyModel>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::ffi::CString;

    use pyo3::types::PyDict;

    use super::*;

    const MODEL_JSON: &str = r#"{
        "name": "scalar",
        "restriction": null,
        "orders": [{"name":"QED","expansion_order":99,"hierarchy":1}],
        "parameters": [{
            "name":"ZERO","lhablock":null,"lhacode":null,"nature":"internal",
            "parameter_type":"real","value":[0.0,0.0],"expression":null
        },{
            "name":"mass","lhablock":"MASS","lhacode":[1],"nature":"external",
            "parameter_type":"real","value":[1.0,0.0],"expression":null
        },{
            "name":"double_mass","lhablock":null,"lhacode":null,"nature":"internal",
            "parameter_type":"complex","value":[2.0,0.0],"expression":"2*mass"
        }],
        "particles": [{
            "pdg_code":1,"name":"s","antiname":"s","spin":1,"color":1,
            "mass":"mass","width":"ZERO","texname":"s","antitexname":"s",
            "charge":0.0,"ghost_number":0,"lepton_number":0,"y_charge":0,
            "propagator":"s_prop"
        }],
        "propagators": [{
            "name":"s_prop","particle":"s","numerator":"1","denominator":"P^2-mass^2"
        }],
        "lorentz_structures": [{"name":"L1","spins":[1,1,1],"structure":"1"}],
        "couplings": [{
            "name":"GC1","expression":"double_mass","orders":[["QED",1]],"value":[2.0,0.0]
        }],
        "vertex_rules": [{
            "name":"V1","particles":["s","s","s"],"color_structures":["1"],
            "lorentz_structures":["L1"],"couplings":[["GC1"]]
        }],
        "functions":[{"name":"twice","arguments":["x"],"expression":"2*x"}],
        "form_factors":[{"name":"FF1","type":"complex","value":"twice(P(1)^2)"}]
    }"#;

    fn registered_module<'py>(py: Python<'py>) -> Bound<'py, PyModule> {
        let module = PyModule::new(py, "symbolica.community.feynkit").unwrap();
        crate::initialize_feynkit(&module).unwrap();
        module
    }

    #[test]
    fn exposes_typed_model_entities_and_atomic_python_evaluation() {
        Python::initialize();
        Python::attach(|py| {
            let module = registered_module(py);
            let locals = PyDict::new(py);
            locals.set_item("fk", &module).unwrap();
            locals.set_item("MODEL_JSON", MODEL_JSON).unwrap();
            let code = CString::new(
                r#"
model = fk.Model.from_json(MODEL_JSON)

assert len(model.parameters) == 3
mass = model.parameter("mass")
assert isinstance(mass, fk.Parameter)
assert mass.name == "mass"
assert mass.lhablock == "MASS"
assert mass.lhacode == [1]
assert mass.nature == fk.ParameterNature.EXTERNAL
assert mass.parameter_type == fk.ParameterType.REAL
assert mass.value == (1.0, 0.0)
try:
    mass.name = "changed"
except AttributeError:
    pass
else:
    raise AssertionError("model entity wrappers must be immutable")

internal = model.parameter("double_mass")
assert internal.nature == fk.ParameterNature.INTERNAL
assert internal.parameter_type == fk.ParameterType.COMPLEX
assert internal.expression == "2*mass"

coupling = model.coupling("GC1")
assert isinstance(coupling, fk.Coupling)
assert coupling.orders == {"QED": 1}
assert coupling.value == (2.0, 0.0)

vertex = model.vertex_rule("V1")
assert isinstance(vertex, fk.VertexRule)
assert vertex.particles == ["s", "s", "s"]
assert vertex.color_structures == ["1"]
assert vertex.lorentz_structures == ["L1"]
assert vertex.couplings == [["GC1"]]
assert vertex.coupling_orders(model) == {"QED": 1}

lorentz = model.lorentz_structure("L1")
assert isinstance(lorentz, fk.LorentzStructure)
assert lorentz.spins == [1, 1, 1]
assert lorentz.structure == "1"

propagator = model.propagator("s_prop")
assert isinstance(propagator, fk.Propagator)
assert propagator.particle == "s"
assert propagator.numerator == "1"
assert propagator.denominator == "P^2-mass^2"

function = model.function("twice")
assert isinstance(function, fk.ModelFunction)
assert function.arguments == ["x"]
assert function.expression == "2*x"

form_factor = model.form_factor("FF1")
assert isinstance(form_factor, fk.FormFactor)
assert form_factor.type_name == "complex"
assert form_factor.value == "twice(P(1)^2)"

assert all(isinstance(value, fk.Parameter) for value in model.parameters)
assert all(isinstance(value, fk.Coupling) for value in model.couplings)
assert all(isinstance(value, fk.VertexRule) for value in model.vertex_rules)
assert all(isinstance(value, fk.LorentzStructure) for value in model.lorentz_structures)
assert all(isinstance(value, fk.Propagator) for value in model.propagators)
assert all(isinstance(value, fk.ModelFunction) for value in model.functions)
assert all(isinstance(value, fk.FormFactor) for value in model.form_factors)

requests = []
def evaluate(request):
    assert isinstance(request, fk.EvaluationRequest)
    assert all(isinstance(value, fk.ModelExpression) for value in request.internal_parameters)
    assert all(isinstance(value, fk.ModelExpression) for value in request.couplings)
    assert request.internal_parameters[0].name == "double_mass"
    assert request.internal_parameters[0].expression == "2*mass"
    assert request.couplings[0].name == "GC1"
    assert request.functions[0].name == "twice"
    assert request.form_factors[0].name == "FF1"
    requests.append(request)
    real, imaginary = request.known_parameters["mass"]
    value = (2.0 * real, 2.0 * imaginary)
    return fk.EvaluatedValues(
        internal_parameters={"double_mass": value},
        couplings={"GC1": value},
    )

recomputed = model.recompute_with(evaluate)
assert isinstance(recomputed, fk.Model)
assert recomputed.parameter("double_mass").value == (2.0, 0.0)
assert recomputed.coupling("GC1").value == (2.0, 0.0)
assert model.parameter("double_mass").value == (2.0, 0.0)

card = fk.ParameterCard()
card.set("mass", 3.0, 0.5)
updated = model.with_parameter_card(card, evaluator=evaluate)
assert updated.parameter("mass").value == (3.0, 0.5)
assert updated.parameter("double_mass").value == (6.0, 1.0)
assert updated.coupling("GC1").value == (6.0, 1.0)
assert requests[-1].known_parameters["mass"] == (3.0, 0.5)

invalidated = model.with_parameter_card(card)
assert invalidated.parameter("mass").value == (3.0, 0.5)
assert invalidated.parameter("double_mass").value is None
assert invalidated.coupling("GC1").value is None

before = model.to_json(pretty=False)
def fail(_request):
    raise ValueError("evaluator failed")

try:
    model.with_parameter_card(card, evaluator=fail)
except ValueError as error:
    assert str(error) == "evaluator failed"
else:
    raise AssertionError("the evaluator exception must be preserved")
assert model.to_json(pretty=False) == before

try:
    model.recompute_with(lambda _request: fk.EvaluatedValues())
except fk.ModelError as error:
    assert "did not return a value" in str(error)
else:
    raise AssertionError("incomplete evaluated values must be rejected")
assert model.to_json(pretty=False) == before

try:
    model.recompute_with(lambda _request: {})
except TypeError:
    pass
else:
    raise AssertionError("the evaluator must return EvaluatedValues")
assert model.to_json(pretty=False) == before
"#,
            )
            .unwrap();
            py.run(&code, Some(&locals), Some(&locals)).unwrap();
        });
    }
}
