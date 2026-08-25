use feynkit_graph::{
    DiagramEdge, DiagramVertex, FeynmanDiagram, LoopMomentumBasis, MomentumSignature,
};
use pyo3::{prelude::*, types::PyModule};
use symbolica::{api::python::PythonExpression, atom::Atom, parser::ParseSettings, wrap_input};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::{error, model::PyModel};

pub(crate) fn parse_symbolic_annotation(value: &str) -> Result<PythonExpression, String> {
    Atom::parse_with_default_namespace(wrap_input!(value), ParseSettings::default())
        .map(|expr| PythonExpression { expr })
        .map_err(|parse_error| parse_error.to_string())
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "DiagramVertex",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyDiagramVertex {
    #[pyo3(get)]
    id: usize,
    inner: DiagramVertex,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyDiagramVertex {
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    #[getter]
    fn interaction(&self) -> Option<String> {
        self.inner.interaction.clone()
    }

    #[getter]
    fn numerator(&self) -> Option<String> {
        self.inner.numerator.clone()
    }

    fn numerator_expression(&self) -> PyResult<Option<PythonExpression>> {
        self.inner
            .numerator
            .as_deref()
            .map(parse_symbolic_annotation)
            .transpose()
            .map_err(error::DiagramError::new_err)
    }

    #[getter]
    fn external_index(&self) -> Option<usize> {
        self.inner.external.as_ref().map(|leg| leg.index)
    }

    #[getter]
    fn external_state(&self) -> Option<&'static str> {
        self.inner.external.as_ref().map(|leg| match leg.state {
            feynkit_graph::ExternalState::Incoming => "incoming",
            feynkit_graph::ExternalState::Outgoing => "outgoing",
        })
    }

    #[getter]
    fn is_external(&self) -> bool {
        self.inner.is_external()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "DiagramEdge",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyDiagramEdge {
    #[pyo3(get)]
    id: usize,
    #[pyo3(get)]
    source: usize,
    #[pyo3(get)]
    target: usize,
    inner: DiagramEdge,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyDiagramEdge {
    #[getter]
    fn particle_name(&self) -> &str {
        &self.inner.particle.name
    }

    #[getter]
    fn particle_pdg(&self) -> i64 {
        self.inner.particle.pdg
    }

    #[getter]
    fn directed(&self) -> bool {
        self.inner.directed
    }

    #[getter]
    fn numerator(&self) -> Option<String> {
        self.inner.numerator.clone()
    }

    fn numerator_expression(&self) -> PyResult<Option<PythonExpression>> {
        self.inner
            .numerator
            .as_deref()
            .map(parse_symbolic_annotation)
            .transpose()
            .map_err(error::DiagramError::new_err)
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "MomentumSignature",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyMomentumSignature {
    inner: MomentumSignature,
}

impl From<MomentumSignature> for PyMomentumSignature {
    fn from(inner: MomentumSignature) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyMomentumSignature {
    #[getter]
    fn loops(&self) -> Vec<isize> {
        self.inner.loops.integer_coefficients()
    }

    #[getter]
    fn external(&self) -> Vec<isize> {
        self.inner.external.integer_coefficients()
    }

    fn integer_coefficients(&self) -> (Vec<isize>, Vec<isize>) {
        self.inner.integer_coefficients()
    }

    fn format_momentum(&self) -> String {
        self.inner.format_momentum()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "LoopMomentumBasis",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyLoopMomentumBasis {
    inner: LoopMomentumBasis,
}

impl From<LoopMomentumBasis> for PyLoopMomentumBasis {
    fn from(inner: LoopMomentumBasis) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyLoopMomentumBasis {
    #[getter]
    fn tree_edges(&self) -> Vec<usize> {
        self.inner.tree_edges.iter().map(|id| id.0).collect()
    }

    #[getter]
    fn loop_edges(&self) -> Vec<usize> {
        self.inner.loop_edges.iter().map(|id| id.0).collect()
    }

    #[getter]
    fn external_edges(&self) -> Vec<usize> {
        self.inner.external_edges.iter().map(|id| id.0).collect()
    }

    #[getter]
    fn dependent_externals(&self) -> Vec<usize> {
        self.inner
            .dependent_externals
            .iter()
            .map(|id| id.0)
            .collect()
    }

    #[getter]
    fn edge_signatures(&self) -> Vec<(usize, PyMomentumSignature)> {
        self.inner
            .edge_signatures
            .iter()
            .map(|(edge, signature)| (edge.0, signature.clone().into()))
            .collect()
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "FeynmanDiagram",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyFeynmanDiagram {
    pub(crate) inner: FeynmanDiagram,
}

impl From<FeynmanDiagram> for PyFeynmanDiagram {
    fn from(inner: FeynmanDiagram) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyFeynmanDiagram {
    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        FeynmanDiagram::from_json(json)
            .map(Into::into)
            .map_err(error::diagram)
    }

    #[staticmethod]
    fn from_dot(dot: &str) -> PyResult<Self> {
        FeynmanDiagram::from_dot(dot)
            .map(Into::into)
            .map_err(error::diagram)
    }

    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    #[getter]
    fn symmetry_factor(&self) -> u64 {
        self.inner.symmetry_factor()
    }

    #[getter]
    fn overall_factor(&self) -> &str {
        self.inner.overall_factor()
    }

    #[getter]
    fn numerator(&self) -> Option<String> {
        self.inner.numerator().map(str::to_owned)
    }

    fn numerator_expression(&self) -> PyResult<Option<PythonExpression>> {
        self.inner
            .numerator()
            .map(parse_symbolic_annotation)
            .transpose()
            .map_err(error::DiagramError::new_err)
    }

    fn overall_factor_expression(&self) -> PyResult<PythonExpression> {
        parse_symbolic_annotation(self.inner.overall_factor()).map_err(error::DiagramError::new_err)
    }

    #[getter]
    fn loop_count(&self) -> usize {
        self.inner.loop_count()
    }

    #[getter]
    fn vertices(&self) -> Vec<PyDiagramVertex> {
        self.inner
            .vertices()
            .map(|(id, inner)| PyDiagramVertex {
                id: id.0,
                inner: inner.clone(),
            })
            .collect()
    }

    #[getter]
    fn edges(&self) -> Vec<PyDiagramEdge> {
        self.inner
            .edges()
            .map(|(id, endpoints, inner)| PyDiagramEdge {
                id: id.0,
                source: endpoints.source.0,
                target: endpoints.target.0,
                inner: inner.clone(),
            })
            .collect()
    }

    fn validate(&self, model: &PyModel) -> PyResult<()> {
        self.inner.validate(&model.inner).map_err(error::diagram)
    }

    fn to_json(&self) -> PyResult<String> {
        self.inner.to_json().map_err(error::diagram)
    }

    fn to_dot(&self) -> PyResult<String> {
        self.inner.to_dot().map_err(error::diagram)
    }

    #[pyo3(signature = (limit=None))]
    fn loop_momentum_bases(
        &self,
        py: Python<'_>,
        limit: Option<usize>,
    ) -> PyResult<Vec<PyLoopMomentumBasis>> {
        let diagram = self.inner.clone();
        py.detach(move || match limit {
            Some(limit) => diagram.loop_momentum_bases_with_limit(limit),
            None => diagram.loop_momentum_bases(),
        })
        .map(|bases| bases.into_iter().map(Into::into).collect())
        .map_err(error::diagram)
    }

    fn __repr__(&self) -> String {
        format!(
            "FeynmanDiagram(name='{}', loops={})",
            self.inner.name(),
            self.inner.loop_count()
        )
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyDiagramVertex>()?;
    module.add_class::<PyDiagramEdge>()?;
    module.add_class::<PyMomentumSignature>()?;
    module.add_class::<PyLoopMomentumBasis>()?;
    module.add_class::<PyFeynmanDiagram>()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::parse_symbolic_annotation;

    #[test]
    fn symbolic_annotations_are_parsed_fallibly() {
        assert!(parse_symbolic_annotation("g^2*(x+y)").is_ok());
        assert!(parse_symbolic_annotation("g+").is_err());
    }
}
