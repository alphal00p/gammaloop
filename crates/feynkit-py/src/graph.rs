use feynkit_graph::{
    DiagramEdge, DiagramVertex, FeynmanDiagram, LoopMomentumBasis, MomentumSignature,
};
use pyo3::{
    prelude::*,
    types::{PyAny, PyModule},
};
use symbolica::{api::python::PythonExpression, atom::Atom, parser::ParseSettings, wrap_input};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

use crate::{display::escape_html, error, model::PyModel};

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
    /// Return this vertex's integer ID within the diagram.
    ///
    /// Examples
    /// --------
    /// >>> vertex = diagram.vertices[0]
    /// >>> isinstance(vertex.id, int)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[pyo3(get)]
    id: usize,
    inner: DiagramVertex,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyDiagramVertex {
    /// Return the vertex name stored in the diagram.
    ///
    /// Examples
    /// --------
    /// Obtain a vertex from a diagram and inspect its name:
    ///
    /// >>> vertex = diagram.vertices[0]
    /// >>> isinstance(vertex.name, str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Return the interaction name for an internal vertex, if present.
    ///
    /// Examples
    /// --------
    /// >>> vertex = next(vertex for vertex in diagram.vertices if not vertex.is_external)
    /// >>> vertex.interaction is None or isinstance(vertex.interaction, str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn interaction(&self) -> Option<String> {
        self.inner.interaction.clone()
    }

    /// Return the symbolic numerator annotation as source text, if present.
    ///
    /// Examples
    /// --------
    /// >>> vertex = diagram.vertices[0]
    /// >>> vertex.numerator is None or isinstance(vertex.numerator, str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn numerator(&self) -> Option<String> {
        self.inner.numerator.clone()
    }

    /// Parse the optional numerator annotation as a Symbolica expression.
    ///
    /// Examples
    /// --------
    /// >>> vertex = diagram.vertices[0]
    /// >>> numerator = vertex.numerator_expression()
    /// >>> numerator is None or isinstance(str(numerator), str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn numerator_expression(&self) -> PyResult<Option<PythonExpression>> {
        self.inner
            .numerator
            .as_deref()
            .map(parse_symbolic_annotation)
            .transpose()
            .map_err(error::DiagramError::new_err)
    }

    /// Return the external-leg index, or ``None`` for an internal vertex.
    ///
    /// Examples
    /// --------
    /// >>> vertex = diagram.vertices[0]
    /// >>> vertex.external_index is None or vertex.external_index >= 0
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn external_index(&self) -> Option<usize> {
        self.inner.external.as_ref().map(|leg| leg.index)
    }

    /// Return ``"incoming"`` or ``"outgoing"`` for an external vertex.
    ///
    /// Examples
    /// --------
    /// >>> vertex = next(vertex for vertex in diagram.vertices if vertex.is_external)
    /// >>> vertex.external_state in {"incoming", "outgoing"}
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn external_state(&self) -> Option<&'static str> {
        self.inner.external.as_ref().map(|leg| match leg.state {
            feynkit_graph::ExternalState::Incoming => "incoming",
            feynkit_graph::ExternalState::Outgoing => "outgoing",
        })
    }

    /// Report whether this vertex represents an external particle state.
    ///
    /// Examples
    /// --------
    /// >>> vertex = diagram.vertices[0]
    /// >>> isinstance(vertex.is_external, bool)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn is_external(&self) -> bool {
        self.inner.is_external()
    }

    /// Return a concise, unambiguous description of the vertex.
    ///
    /// Examples
    /// --------
    /// >>> vertex = diagram.vertices[0]
    /// >>> repr(vertex).startswith("DiagramVertex(")
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        if let Some(external) = &self.inner.external {
            format!(
                "DiagramVertex(id={}, name={:?}, external_state={:?}, external_index={})",
                self.id,
                self.inner.name,
                match external.state {
                    feynkit_graph::ExternalState::Incoming => "incoming",
                    feynkit_graph::ExternalState::Outgoing => "outgoing",
                },
                external.index,
            )
        } else {
            format!(
                "DiagramVertex(id={}, name={:?}, interaction={:?})",
                self.id, self.inner.name, self.inner.interaction,
            )
        }
    }

    /// Write the vertex summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython calls this method when formatting a vertex for text display.
    ///
    /// Parameters
    /// ----------
    /// pretty:
    ///     The IPython pretty-printer object.
    /// cycle:
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

#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "DiagramEdge",
    module = "symbolica.community.feynkit",
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyDiagramEdge {
    /// Return this edge's integer ID within the diagram.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> isinstance(edge.id, int)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[pyo3(get)]
    id: usize,
    /// Return the source vertex ID.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> isinstance(edge.source, int)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[pyo3(get)]
    source: usize,
    /// Return the target vertex ID.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> isinstance(edge.target, int)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[pyo3(get)]
    target: usize,
    inner: DiagramEdge,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[pymethods]
impl PyDiagramEdge {
    /// Return the particle name associated with this propagator edge.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> isinstance(edge.particle_name, str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn particle_name(&self) -> &str {
        &self.inner.particle.name
    }

    /// Return the PDG code of the particle carried by this edge.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> isinstance(edge.particle_pdg, int)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn particle_pdg(&self) -> i64 {
        self.inner.particle.pdg
    }

    /// Report whether the edge has a particle-flow direction.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> isinstance(edge.directed, bool)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn directed(&self) -> bool {
        self.inner.directed
    }

    /// Return the symbolic numerator annotation as source text, if present.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> edge.numerator is None or isinstance(edge.numerator, str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn numerator(&self) -> Option<String> {
        self.inner.numerator.clone()
    }

    /// Parse the optional numerator annotation as a Symbolica expression.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> numerator = edge.numerator_expression()
    /// >>> numerator is None or isinstance(str(numerator), str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn numerator_expression(&self) -> PyResult<Option<PythonExpression>> {
        self.inner
            .numerator
            .as_deref()
            .map(parse_symbolic_annotation)
            .transpose()
            .map_err(error::DiagramError::new_err)
    }

    /// Return a concise description of the edge and its endpoints.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> repr(edge).startswith("DiagramEdge(")
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        let connector = if self.inner.directed { "->" } else { "--" };
        format!(
            "DiagramEdge(id={}, {}{}{}, particle={:?}, pdg={})",
            self.id,
            self.source,
            connector,
            self.target,
            self.inner.particle.name,
            self.inner.particle.pdg,
        )
    }

    /// Write the edge summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython calls this method when formatting an edge for text display.
    ///
    /// Parameters
    /// ----------
    /// pretty:
    ///     The IPython pretty-printer object.
    /// cycle:
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
    /// Return the integer coefficients of the independent loop momenta.
    ///
    /// Examples
    /// --------
    /// Obtain a signature from a loop-momentum basis:
    ///
    /// >>> _, signature = basis.edge_signatures[0]
    /// >>> all(isinstance(coefficient, int) for coefficient in signature.loops)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn loops(&self) -> Vec<isize> {
        self.inner.loops.integer_coefficients()
    }

    /// Return the integer coefficients of the external momenta.
    ///
    /// Examples
    /// --------
    /// >>> _, signature = basis.edge_signatures[0]
    /// >>> all(isinstance(coefficient, int) for coefficient in signature.external)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn external(&self) -> Vec<isize> {
        self.inner.external.integer_coefficients()
    }

    /// Return loop and external coefficients as a pair of integer lists.
    ///
    /// Examples
    /// --------
    /// >>> _, signature = basis.edge_signatures[0]
    /// >>> loops, external = signature.integer_coefficients()
    /// >>> (loops, external) == (signature.loops, signature.external)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn integer_coefficients(&self) -> (Vec<isize>, Vec<isize>) {
        self.inner.integer_coefficients()
    }

    /// Format the signature as a conventional momentum sum.
    ///
    /// Examples
    /// --------
    /// >>> _, signature = basis.edge_signatures[0]
    /// >>> isinstance(signature.format_momentum(), str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn format_momentum(&self) -> String {
        self.inner.format_momentum()
    }

    /// Return the conventional momentum sum for ``str(signature)``.
    ///
    /// Examples
    /// --------
    /// >>> _, signature = basis.edge_signatures[0]
    /// >>> str(signature) == signature.format_momentum()
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn __str__(&self) -> String {
        self.inner.format_momentum()
    }

    /// Return an unambiguous momentum-signature summary.
    ///
    /// Examples
    /// --------
    /// >>> _, signature = basis.edge_signatures[0]
    /// >>> repr(signature).startswith("MomentumSignature(")
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        format!("MomentumSignature({})", self.inner.format_momentum())
    }

    /// Write the momentum signature to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython calls this method when formatting a signature for text display.
    ///
    /// Parameters
    /// ----------
    /// pretty:
    ///     The IPython pretty-printer object.
    /// cycle:
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
    /// Return the edge identifiers belonging to the spanning tree.
    ///
    /// Examples
    /// --------
    /// >>> basis = diagram.loop_momentum_bases(limit=1)[0]
    /// >>> all(isinstance(edge_id, int) for edge_id in basis.tree_edges)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn tree_edges(&self) -> Vec<usize> {
        self.inner.tree_edges.iter().map(|id| id.0).collect()
    }

    /// Return the edge identifiers chosen as independent loop momenta.
    ///
    /// Examples
    /// --------
    /// >>> basis = diagram.loop_momentum_bases(limit=1)[0]
    /// >>> len(basis.loop_edges) == diagram.loop_count
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn loop_edges(&self) -> Vec<usize> {
        self.inner.loop_edges.iter().map(|id| id.0).collect()
    }

    /// Return the identifiers of edges attached to external states.
    ///
    /// Examples
    /// --------
    /// >>> basis = diagram.loop_momentum_bases(limit=1)[0]
    /// >>> all(isinstance(edge_id, int) for edge_id in basis.external_edges)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn external_edges(&self) -> Vec<usize> {
        self.inner.external_edges.iter().map(|id| id.0).collect()
    }

    /// Return external-edge identifiers fixed by momentum conservation.
    ///
    /// Examples
    /// --------
    /// >>> basis = diagram.loop_momentum_bases(limit=1)[0]
    /// >>> set(basis.dependent_externals) <= set(basis.external_edges)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn dependent_externals(&self) -> Vec<usize> {
        self.inner
            .dependent_externals
            .iter()
            .map(|id| id.0)
            .collect()
    }

    /// Return each edge identifier together with its momentum signature.
    ///
    /// Examples
    /// --------
    /// >>> basis = diagram.loop_momentum_bases(limit=1)[0]
    /// >>> edge_id, signature = basis.edge_signatures[0]
    /// >>> isinstance(edge_id, int) and isinstance(signature.format_momentum(), str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn edge_signatures(&self) -> Vec<(usize, PyMomentumSignature)> {
        self.inner
            .edge_signatures
            .iter()
            .map(|(edge, signature)| (edge.0, signature.clone().into()))
            .collect()
    }

    /// Return a concise description of the selected momentum basis.
    ///
    /// Examples
    /// --------
    /// >>> basis = diagram.loop_momentum_bases(limit=1)[0]
    /// >>> repr(basis).startswith("LoopMomentumBasis(")
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn __repr__(&self) -> String {
        format!(
            "LoopMomentumBasis(loop_edges={:?}, tree_edges={:?}, external_edges={:?})",
            self.inner
                .loop_edges
                .iter()
                .map(|id| id.0)
                .collect::<Vec<_>>(),
            self.inner
                .tree_edges
                .iter()
                .map(|id| id.0)
                .collect::<Vec<_>>(),
            self.inner
                .external_edges
                .iter()
                .map(|id| id.0)
                .collect::<Vec<_>>(),
        )
    }

    /// Return an HTML table of edge momentum assignments for notebooks.
    ///
    /// Examples
    /// --------
    /// Leave ``basis`` as the final expression in a Jupyter or Marimo cell, or
    /// obtain the HTML explicitly:
    ///
    /// >>> "<table" in basis._repr_html_()
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn _repr_html_(&self) -> String {
        let rows = self
            .inner
            .edge_signatures
            .iter()
            .map(|(edge, signature)| {
                format!(
                    "<tr><td style=\"padding:.2rem .6rem\">{}</td>\
                     <td style=\"padding:.2rem .6rem\"><code>{}</code></td></tr>",
                    edge.0,
                    escape_html(&signature.format_momentum()),
                )
            })
            .collect::<String>();
        format!(
            "<div class=\"feynkit-loop-momentum-basis\" style=\"display:inline-block;\
             max-width:100%;overflow-x:auto\"><strong>Loop-momentum basis</strong>\
             <table style=\"border-collapse:collapse;margin-top:.25rem\"><thead><tr>\
             <th style=\"padding:.2rem .6rem;text-align:left\">edge</th>\
             <th style=\"padding:.2rem .6rem;text-align:left\">momentum</th>\
             </tr></thead><tbody>{rows}</tbody></table></div>"
        )
    }

    /// Write the basis summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython calls this method when formatting a basis for text display.
    ///
    /// Parameters
    /// ----------
    /// pretty:
    ///     The IPython pretty-printer object.
    /// cycle:
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
    /// Deserialize a Feynman diagram from its JSON representation.
    ///
    /// Examples
    /// --------
    /// >>> encoded = diagram.to_json()
    /// >>> restored = FeynmanDiagram.from_json(encoded)
    /// >>> restored.name == diagram.name
    /// True
    ///
    /// Parameters
    /// ----------
    /// json:
    ///     JSON text produced by :meth:`FeynmanDiagram.to_json` or another
    ///     schema-compatible producer.
    #[staticmethod]
    fn from_json(json: &str) -> PyResult<Self> {
        FeynmanDiagram::from_json(json)
            .map(Into::into)
            .map_err(error::diagram)
    }

    /// Parse a Feynman diagram from Graphviz DOT text.
    ///
    /// Examples
    /// --------
    /// >>> dot = diagram.to_dot()
    /// >>> restored = FeynmanDiagram.from_dot(dot)
    /// >>> restored.loop_count == diagram.loop_count
    /// True
    ///
    /// Parameters
    /// ----------
    /// dot:
    ///     DOT text containing the diagram topology and FeynKit annotations.
    #[staticmethod]
    fn from_dot(dot: &str) -> PyResult<Self> {
        FeynmanDiagram::from_dot(dot)
            .map(Into::into)
            .map_err(error::diagram)
    }

    /// Return the model-assigned diagram name.
    ///
    /// Examples
    /// --------
    /// >>> isinstance(diagram.name, str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    /// Return the diagram symmetry factor as a positive integer denominator.
    ///
    /// Examples
    /// --------
    /// >>> diagram.symmetry_factor >= 1
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn symmetry_factor(&self) -> u64 {
        self.inner.symmetry_factor()
    }

    /// Return the diagram-wide multiplicative factor as Symbolica source text.
    ///
    /// Examples
    /// --------
    /// >>> isinstance(diagram.overall_factor, str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn overall_factor(&self) -> &str {
        self.inner.overall_factor()
    }

    /// Return the diagram numerator annotation as source text, if present.
    ///
    /// Examples
    /// --------
    /// >>> diagram.numerator is None or isinstance(diagram.numerator, str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn numerator(&self) -> Option<String> {
        self.inner.numerator().map(str::to_owned)
    }

    /// Parse the optional diagram numerator as a Symbolica expression.
    ///
    /// Examples
    /// --------
    /// >>> numerator = diagram.numerator_expression()
    /// >>> numerator is None or isinstance(str(numerator), str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn numerator_expression(&self) -> PyResult<Option<PythonExpression>> {
        self.inner
            .numerator()
            .map(parse_symbolic_annotation)
            .transpose()
            .map_err(error::DiagramError::new_err)
    }

    /// Parse the diagram-wide factor as a Symbolica expression.
    ///
    /// Examples
    /// --------
    /// >>> factor = diagram.overall_factor_expression()
    /// >>> isinstance(str(factor), str)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn overall_factor_expression(&self) -> PyResult<PythonExpression> {
        parse_symbolic_annotation(self.inner.overall_factor()).map_err(error::DiagramError::new_err)
    }

    /// Return the number of independent loops in the diagram topology.
    ///
    /// Examples
    /// --------
    /// >>> diagram.loop_count >= 0
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    #[getter]
    fn loop_count(&self) -> usize {
        self.inner.loop_count()
    }

    /// Return the diagram vertices with stable integer identifiers.
    ///
    /// Examples
    /// --------
    /// >>> vertices = diagram.vertices
    /// >>> all(vertex.id == index for index, vertex in enumerate(vertices))
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
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

    /// Return the diagram edges with their endpoint identifiers.
    ///
    /// Examples
    /// --------
    /// >>> edges = diagram.edges
    /// >>> all(edge.source >= 0 and edge.target >= 0 for edge in edges)
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
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

    /// Validate particle and interaction references against a physics model.
    ///
    /// Examples
    /// --------
    /// Validation returns ``None`` and raises ``DiagramError`` for invalid
    /// references:
    ///
    /// >>> diagram.validate(model) is None
    /// True
    ///
    /// Parameters
    /// ----------
    /// model:
    ///     The :class:`Model` whose particle and interaction definitions should
    ///     be used for validation.
    fn validate(&self, model: &PyModel) -> PyResult<()> {
        self.inner.validate(&model.inner).map_err(error::diagram)
    }

    /// Serialize the complete diagram to JSON text.
    ///
    /// Examples
    /// --------
    /// >>> encoded = diagram.to_json()
    /// >>> encoded.lstrip().startswith("{")
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn to_json(&self) -> PyResult<String> {
        self.inner.to_json().map_err(error::diagram)
    }

    /// Serialize the diagram topology and annotations to Graphviz DOT text.
    ///
    /// Examples
    /// --------
    /// >>> dot = diagram.to_dot()
    /// >>> "graph" in dot.lower()
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn to_dot(&self) -> PyResult<String> {
        self.inner.to_dot().map_err(error::diagram)
    }

    /// Render the diagram as a self-contained, responsive SVG.
    ///
    /// Examples
    /// --------
    /// >>> svg = diagram.to_svg()
    /// >>> svg.startswith("<svg")
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn to_svg(&self) -> String {
        self.inner.to_svg()
    }

    /// Render the diagram as a self-contained HTML figure.
    ///
    /// Examples
    /// --------
    /// Embed the returned markup in a web page or notebook component:
    ///
    /// >>> html = diagram.to_html()
    /// >>> "<figure" in html and "<svg" in html
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn to_html(&self) -> String {
        self.inner.to_html()
    }

    /// Render the diagram as HTML in Marimo, Jupyter, and IPython.
    ///
    /// Examples
    /// --------
    /// Leave `diagram` as the final expression in a notebook cell to render it.
    ///
    /// Parameters
    /// ----------
    /// None
    fn _repr_html_(&self) -> String {
        self.inner.to_html()
    }

    /// Return the raw SVG representation used by rich notebook frontends.
    ///
    /// Examples
    /// --------
    /// >>> "<svg" in diagram._repr_svg_()
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
    fn _repr_svg_(&self) -> String {
        self.inner.to_svg()
    }

    /// Write a concise summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython invokes this method when only a text representation is supported.
    ///
    /// Parameters
    /// ----------
    /// pretty:
    ///     The IPython pretty-printer object.
    /// cycle:
    ///     Whether this object is part of a recursive formatting cycle.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        let text = if cycle {
            "...".to_owned()
        } else {
            format!(
                "FeynmanDiagram(name='{}', loops={}, vertices={}, edges={})",
                self.inner.name(),
                self.inner.loop_count(),
                self.inner.vertices().count(),
                self.inner.edges().count(),
            )
        };
        pretty.call_method1("text", (text,))?;
        Ok(())
    }

    /// Enumerate valid loop-momentum bases for this diagram.
    ///
    /// Examples
    /// --------
    /// Limit exploratory calculations to the first basis:
    ///
    /// >>> bases = diagram.loop_momentum_bases(limit=1)
    /// >>> len(bases) <= 1
    /// True
    ///
    /// Parameters
    /// ----------
    /// limit:
    ///     Maximum number of bases to return. Pass ``None`` to enumerate every
    ///     valid basis.
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

    /// Return a concise description of the diagram and its loop count.
    ///
    /// Examples
    /// --------
    /// >>> repr(diagram).startswith("FeynmanDiagram(")
    /// True
    ///
    /// Parameters
    /// ----------
    /// None
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
