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

use crate::error;

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
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[getter]
    fn antiname(&self) -> &str {
        &self.inner.antiname
    }

    #[getter]
    fn pdg_code(&self) -> i64 {
        self.inner.pdg_code
    }

    #[getter]
    fn spin(&self) -> i64 {
        self.inner.spin
    }

    #[getter]
    fn color(&self) -> i64 {
        self.inner.color
    }

    #[getter]
    fn mass_parameter(&self) -> &str {
        &self.inner.mass
    }

    #[getter]
    fn width_parameter(&self) -> &str {
        &self.inner.width
    }

    #[getter]
    fn charge(&self) -> f64 {
        self.inner.charge
    }

    #[getter]
    fn is_antiparticle(&self) -> bool {
        self.inner.is_antiparticle()
    }

    #[getter]
    fn is_self_antiparticle(&self) -> bool {
        self.inner.is_self_antiparticle()
    }

    #[getter]
    fn is_fermion(&self) -> bool {
        self.inner.is_fermion()
    }

    #[getter]
    fn is_massless(&self) -> bool {
        self.inner.is_massless()
    }

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
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[getter]
    fn lhablock(&self) -> Option<String> {
        self.inner.lhablock.clone()
    }

    #[getter]
    fn lhacode(&self) -> Option<Vec<usize>> {
        self.inner.lhacode.clone()
    }

    #[getter]
    fn nature(&self) -> PyParameterNature {
        self.inner.nature.clone().into()
    }

    #[getter]
    fn parameter_type(&self) -> PyParameterType {
        self.inner.parameter_type.clone().into()
    }

    #[getter]
    fn value(&self) -> Option<(f64, f64)> {
        self.inner.value.map(Into::into)
    }

    #[getter]
    fn expression(&self) -> Option<String> {
        self.inner.expression.clone()
    }

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
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[getter]
    fn expression(&self) -> &str {
        &self.inner.expression
    }

    #[getter]
    fn orders(&self) -> BTreeMap<String, usize> {
        self.inner.orders.clone()
    }

    #[getter]
    fn value(&self) -> Option<(f64, f64)> {
        self.inner.value.map(Into::into)
    }

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
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[getter]
    fn particles(&self) -> Vec<String> {
        self.inner.particles.clone()
    }

    #[getter]
    fn color_structures(&self) -> Vec<String> {
        self.inner.color_structures.clone()
    }

    #[getter]
    fn lorentz_structures(&self) -> Vec<String> {
        self.inner.lorentz_structures.clone()
    }

    #[getter]
    fn couplings(&self) -> Vec<Vec<Option<String>>> {
        self.inner.couplings.clone()
    }

    fn coupling_orders(&self, model: &PyModel) -> PyResult<BTreeMap<String, usize>> {
        self.inner
            .coupling_orders(&model.inner)
            .map_err(error::model)
    }

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
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[getter]
    fn spins(&self) -> Vec<i64> {
        self.inner.spins.clone()
    }

    #[getter]
    fn structure(&self) -> &str {
        &self.inner.structure
    }

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
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[getter]
    fn particle(&self) -> &str {
        &self.inner.particle
    }

    #[getter]
    fn numerator(&self) -> &str {
        &self.inner.numerator
    }

    #[getter]
    fn denominator(&self) -> &str {
        &self.inner.denominator
    }

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
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[getter]
    fn arguments(&self) -> Vec<String> {
        self.inner.arguments.clone()
    }

    #[getter]
    fn expression(&self) -> Option<String> {
        self.inner.expression.clone()
    }

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
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[getter]
    fn type_name(&self) -> Option<String> {
        self.inner.type_name.clone()
    }

    #[getter]
    fn value(&self) -> Option<String> {
        self.inner.value.clone()
    }

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
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

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
    #[getter]
    fn known_parameters(&self) -> BTreeMap<String, (f64, f64)> {
        self.inner
            .known_parameters
            .iter()
            .map(|(name, value)| (name.clone(), (*value).into()))
            .collect()
    }

    #[getter]
    fn internal_parameters(&self) -> Vec<PyModelExpression> {
        self.inner
            .internal_parameters
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    #[getter]
    fn couplings(&self) -> Vec<PyModelExpression> {
        self.inner
            .couplings
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    #[getter]
    fn functions(&self) -> Vec<PyModelFunction> {
        self.inner
            .functions
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

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

    #[getter]
    fn internal_parameters(&self) -> BTreeMap<String, (f64, f64)> {
        self.inner
            .internal_parameters
            .iter()
            .map(|(name, value)| (name.clone(), (*value).into()))
            .collect()
    }

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
    #[new]
    fn new() -> Self {
        Self::default()
    }

    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        ParameterCard::from_json(json)
            .map(Self::from)
            .map_err(error::model)
    }

    #[staticmethod]
    fn from_path(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        py.detach(move || ParameterCard::from_path(path))
            .map(Self::from)
            .map_err(error::model)
    }

    fn get(&self, name: &str) -> Option<(f64, f64)> {
        self.inner.get(name).copied().map(Into::into)
    }

    #[pyo3(signature = (name, real, imaginary=0.0))]
    fn set(&mut self, name: String, real: f64, imaginary: f64) {
        self.inner.insert(name, ComplexValue::new(real, imaginary));
    }

    fn remove(&mut self, name: &str) -> Option<(f64, f64)> {
        self.inner.remove(name).map(Into::into)
    }

    fn items(&self) -> Vec<(String, (f64, f64))> {
        self.inner
            .iter()
            .map(|(name, value)| (name.clone(), (*value).into()))
            .collect()
    }

    #[pyo3(signature = (pretty=true))]
    fn to_json(&self, pretty: bool) -> PyResult<String> {
        if pretty {
            self.inner.to_json_pretty()
        } else {
            self.inner.to_json()
        }
        .map_err(error::model)
    }

    fn write_json(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        let card = self.inner.clone();
        py.detach(move || card.write_json(path))
            .map_err(error::model)
    }

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
    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        Model::from_json(json).map(Self::from).map_err(error::model)
    }

    #[staticmethod]
    fn from_path(py: Python<'_>, path: PathBuf) -> PyResult<Self> {
        py.detach(move || Model::from_path(path))
            .map(Self::from)
            .map_err(error::model)
    }

    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    #[getter]
    fn restriction(&self) -> Option<String> {
        self.inner.restriction().map(str::to_owned)
    }

    #[getter]
    fn particles(&self) -> Vec<PyParticle> {
        self.inner
            .particles()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    fn particle(&self, name: &str) -> PyResult<PyParticle> {
        self.inner
            .particle(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    fn particle_by_pdg(&self, pdg: i64) -> PyResult<PyParticle> {
        self.inner
            .particle_by_pdg(pdg)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    #[getter]
    fn parameters(&self) -> Vec<PyParameter> {
        self.inner
            .parameters()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    fn parameter(&self, name: &str) -> PyResult<PyParameter> {
        self.inner
            .parameter(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    #[getter]
    fn couplings(&self) -> Vec<PyCoupling> {
        self.inner
            .couplings()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    fn coupling(&self, name: &str) -> PyResult<PyCoupling> {
        self.inner
            .coupling(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    #[getter]
    fn vertex_rules(&self) -> Vec<PyVertexRule> {
        self.inner
            .vertex_rules()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    fn vertex_rule(&self, name: &str) -> PyResult<PyVertexRule> {
        self.inner
            .vertex_rule(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    #[getter]
    fn lorentz_structures(&self) -> Vec<PyLorentzStructure> {
        self.inner
            .lorentz_structures()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    fn lorentz_structure(&self, name: &str) -> PyResult<PyLorentzStructure> {
        self.inner
            .lorentz_structure(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    #[getter]
    fn propagators(&self) -> Vec<PyPropagator> {
        self.inner
            .propagators()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    fn propagator(&self, name: &str) -> PyResult<PyPropagator> {
        self.inner
            .propagator(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    #[getter]
    fn functions(&self) -> Vec<PyModelFunction> {
        self.inner
            .functions()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    fn function(&self, name: &str) -> PyResult<PyModelFunction> {
        self.inner
            .function(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    #[getter]
    fn form_factors(&self) -> Vec<PyFormFactor> {
        self.inner
            .form_factors()
            .iter()
            .cloned()
            .map(Into::into)
            .collect()
    }

    fn form_factor(&self, name: &str) -> PyResult<PyFormFactor> {
        self.inner
            .form_factor(name)
            .cloned()
            .map(Into::into)
            .map_err(error::model)
    }

    fn default_parameter_card(&self) -> PyResult<PyParameterCard> {
        self.inner
            .default_parameter_card()
            .map(Into::into)
            .map_err(error::model)
    }

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

    #[pyo3(signature = (pretty=true))]
    fn to_json(&self, pretty: bool) -> PyResult<String> {
        if pretty {
            self.inner.to_json_pretty()
        } else {
            self.inner.to_json()
        }
        .map_err(error::model)
    }

    fn write_json(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        let model = self.inner.clone();
        py.detach(move || model.write_json(path))
            .map_err(error::model)
    }

    fn __repr__(&self) -> String {
        format!(
            "Model(name='{}', particles={})",
            self.inner.name(),
            self.inner.particles().len()
        )
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
