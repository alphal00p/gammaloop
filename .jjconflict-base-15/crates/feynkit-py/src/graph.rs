use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use feynkit_graph::expressions::GraphExpressions;
use feynkit_graph::{
    DiagramCut, DiagramCutSide, DiagramEdge, DiagramEndpoint, DiagramHalfEdge,
    DiagramThresholdCandidate, DiagramVertex, FeynmanDiagram, LoopMomentumBasis, MomentumSignature,
};
use feynkit_model::Model;
use feynkit_tensor::FeynmanDiagramTensorExt;
use linnet::half_edge::{
    NodeIndex,
    subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
};
use pyo3::{
    PyTraverseError, PyVisit,
    exceptions::{PyTypeError, PyValueError},
    prelude::*,
    types::{PyAny, PyBytes, PyDict, PyModule, PyTuple},
};
use spynso3::expression::TensorExpression;
use symbolica::{
    api::python::{ConvertibleToExpression, PythonExpression},
    atom::Atom,
    parser::ParseSettings,
    wrap_input,
};

#[cfg(feature = "python_stubgen")]
use pyo3_stub_gen::{
    PyStubType, TypeInfo,
    derive::{gen_methods_from_python, gen_stub_pyclass, gen_stub_pymethods},
    inventory::submit,
};

use crate::{
    cff::{PyCffResult, PyCutPropagator, build_cff_for_diagram},
    display::{escape_html, render_diagram_html, render_diagram_svg},
    error,
    graph_interop::LinnetCache,
    integrals::PyIntegralFamily,
    kinematics::{PyFourMomentum, PyKinematics, PyThreeMomentum},
    model::{PyModel, PyParticle},
    tensor::PyTensorReducer,
};

pub(crate) fn parse_symbolic_annotation(value: &str) -> Result<PythonExpression, String> {
    Atom::parse_with_default_namespace(wrap_input!(value), ParseSettings::default())
        .map(|expr| PythonExpression { expr })
        .map_err(|parse_error| parse_error.to_string())
}

/// An interaction point in a Feynman diagram.
///
/// External states belong to edges; every public vertex references an interaction.
///
/// Examples
/// --------
/// >>> vertex = next(iter(diagram.vertices))
/// >>> interaction = model.vertex_rule(vertex.interaction)
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
    #[pyo3(get)]
    id: usize,
    inner: DiagramVertex,
    model: Arc<Model>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyDiagramVertex {
    /// Return the stable vertex label stored in the diagram.
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }

    /// Return the interaction name.
    ///
    /// Raises :class:`DiagramError` for an incomplete imported diagram.
    ///
    /// Examples
    /// --------
    /// >>> vertex = diagram.vertices[0]
    /// >>> rule = model.vertex_rule(vertex.interaction)
    ///
    #[getter]
    fn interaction(&self) -> PyResult<String> {
        let interaction = self.inner.interaction.ok_or_else(|| {
            error::DiagramError::new_err(format!(
                "vertex {} has no interaction annotation",
                self.id
            ))
        })?;
        self.model
            .vertex_rule_by_id(interaction)
            .map(|rule| rule.name.clone())
            .map_err(error::model)
    }

    /// Return the symbolic numerator annotation as source text.
    ///
    /// Raises :class:`DiagramError` when the diagram does not carry an
    /// instantiated Feynman-rule numerator for this vertex.
    #[getter]
    fn numerator(&self) -> PyResult<String> {
        Ok(self.inner.numerator.to_plain_string())
    }

    /// Return the numerator annotation as a Spenso TensorExpression.
    ///
    /// Examples
    /// --------
    /// >>> vertex = diagram.vertices[0]
    /// >>> vertex_factor = vertex.numerator_expression()
    /// >>> weighted_vertex_factor = diagram.overall_factor_expression() * vertex_factor
    ///
    fn numerator_expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface(py, self.inner.numerator.clone(), None)
    }

    /// Return a concise, unambiguous description of the vertex.
    ///
    /// Examples
    /// --------
    /// >>> vertex = diagram.vertices[0]
    /// >>> print(vertex)  # includes its external state or interaction rule
    ///
    fn __repr__(&self) -> String {
        format!(
            "DiagramVertex(id={}, name={:?}, interaction={:?})",
            self.id, self.inner.name, self.inner.interaction
        )
    }

    /// Write the vertex summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython calls this method when formatting a vertex for text display.
    ///
    /// Parameters
    /// ----------
    /// pretty : Any
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

/// A particle propagator joining two vertices in a Feynman diagram.
///
/// Edges retain their particle identity, endpoints, flow direction, and an
/// symbolic numerator contribution when Feynman rules have been instantiated.
///
/// Examples
/// --------
/// >>> edge = next(iter(diagram.edges))
/// >>> particle = model.particle_by_pdg(edge.particle_pdg)
/// >>> source, target = diagram.vertices[edge.source], diagram.vertices[edge.target]
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
    #[pyo3(get)]
    id: usize,
    /// Return the source interaction ID, or None at a dangling sink.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> source_vertex = diagram.vertices[edge.source]
    ///
    #[pyo3(get)]
    source: Option<usize>,
    /// Return the target interaction ID, or None at a dangling source.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> target_vertex = diagram.vertices[edge.target]
    ///
    #[pyo3(get)]
    target: Option<usize>,
    inner: DiagramEdge,
    model: Arc<Model>,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyDiagramEdge {
    /// Return the external-leg index.
    ///
    /// Raises :class:`DiagramError` for an internal edge. Use
    /// :attr:`is_external` before accessing external-state metadata.
    ///
    /// Examples
    /// --------
    /// >>> external = [v for v in diagram.vertices if v.is_external]
    /// >>> sorted(v.external_index for v in external) == list(range(len(external)))
    /// True
    ///
    #[getter]
    fn external_index(&self) -> PyResult<usize> {
        self.inner
            .external
            .as_ref()
            .map(|leg| leg.index)
            .ok_or_else(|| {
                error::DiagramError::new_err(format!("edge {} is not external", self.id))
            })
    }

    /// Return ``"incoming"`` or ``"outgoing"`` for an external edge.
    ///
    /// Raises :class:`DiagramError` for an internal edge.
    #[getter]
    fn external_state(&self) -> PyResult<&'static str> {
        self.inner
            .external
            .as_ref()
            .map(|leg| match leg.state {
                feynkit_graph::ExternalState::Incoming => "incoming",
                feynkit_graph::ExternalState::Outgoing => "outgoing",
            })
            .ok_or_else(|| {
                error::DiagramError::new_err(format!("edge {} is not external", self.id))
            })
    }

    /// Report whether this edge represents an external particle state.
    ///
    /// Examples
    /// --------
    /// >>> external_edges = [edge for edge in diagram.edges if edge.is_external]
    ///
    #[getter]
    fn is_external(&self) -> bool {
        self.inner.external.is_some()
    }

    /// Return the shared sewing identity of an external momentum carrier.
    #[getter]
    fn external_connection(&self) -> Option<usize> {
        self.inner.external.as_ref().map(|leg| leg.connection)
    }

    /// Return the external-state label retained during finalization.
    #[getter]
    fn external_name(&self) -> Option<String> {
        self.inner.external.as_ref().map(|leg| leg.name.clone())
    }

    /// Report whether this line is a dummy graph attachment.
    #[getter]
    fn is_dummy(&self) -> bool {
        self.inner.is_dummy
    }

    /// Report whether this line has only one incident interaction vertex.
    #[getter]
    fn is_dangling(&self) -> bool {
        self.source.is_none() || self.target.is_none()
    }

    /// Return the particle name associated with this propagator edge.
    #[getter]
    fn particle_name(&self) -> PyResult<String> {
        self.model
            .particle_by_id(self.inner.particle)
            .map(|particle| particle.name.clone())
            .map_err(error::model)
    }

    /// Return the PDG code of the particle carried by this edge.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> particle = model.particle_by_pdg(edge.particle_pdg)
    ///
    #[getter]
    fn particle_pdg(&self) -> PyResult<i64> {
        self.model
            .particle_by_id(self.inner.particle)
            .map(|particle| particle.pdg_code)
            .map_err(error::model)
    }

    /// Report whether the edge carries an oriented particle-flow arrow.
    ///
    /// Examples
    /// --------
    /// >>> fermion_edges = [edge for edge in diagram.edges if edge.directed]
    ///
    #[getter]
    fn directed(&self) -> bool {
        self.inner.directed
    }

    /// Return the symbolic numerator annotation as source text.
    ///
    /// Raises :class:`DiagramError` for an incomplete imported diagram that has
    /// no instantiated propagator numerator.
    #[getter]
    fn numerator(&self) -> PyResult<String> {
        Ok(self.inner.numerator.to_plain_string())
    }

    /// Return the numerator annotation as a Spenso TensorExpression.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> propagator_factor = edge.numerator_expression()
    /// >>> weighted_propagator = diagram.overall_factor_expression() * propagator_factor
    ///
    fn numerator_expression(&self, py: Python<'_>) -> PyResult<Py<TensorExpression>> {
        TensorExpression::from_atom_interface(py, self.inner.numerator.clone(), None)
    }

    /// Return a concise description of the edge and its endpoints.
    ///
    /// Examples
    /// --------
    /// >>> edge = diagram.edges[0]
    /// >>> print(edge)  # shows particle identity, endpoints, and flow direction
    ///
    fn __repr__(&self) -> String {
        let connector = if self.inner.directed { "->" } else { "--" };
        let particle = self
            .model
            .particle_by_id(self.inner.particle)
            .expect("diagram particle IDs are model-validated");
        format!(
            "DiagramEdge(id={}, {:?}{}{:?}, particle={:?}, pdg={})",
            self.id, self.source, connector, self.target, particle.name, particle.pdg_code,
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
    /// pretty : Any
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

/// One amplitude region of a physical cut, with generation metadata.
///
/// Examples
/// --------
/// >>> side = diagram.cuts[0].left
/// >>> side_numerator = side.subgraph.numerator_expression()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "DiagramCutSide",
    module = "symbolica.community.feynkit",
    frozen
)]
pub struct PyDiagramCutSide {
    diagram: Py<PyFeynmanDiagram>,
    inner: DiagramCutSide,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyDiagramCutSide {
    /// Return a physics view of this amplitude side, retaining its original diagram.
    #[getter]
    fn subgraph(&self, py: Python<'_>) -> PyResult<Py<PySubgraph>> {
        self.diagram
            .borrow(py)
            .half_edge_selection(py, &self.inner.half_edges)
    }
    /// Return coupling powers of this amplitude side.
    #[getter]
    fn coupling_orders(&self) -> BTreeMap<String, usize> {
        self.inner.coupling_orders.clone()
    }
    /// Return the number of loops in this amplitude side.
    #[getter]
    fn loop_count(&self) -> usize {
        self.inner.loop_count
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.diagram)
    }
}

/// A generated physical final-state cut; each side retains its physics metadata.
///
/// Examples
/// --------
/// >>> cut = diagram.cuts[0]
/// >>> factors = cut.propagators()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "DiagramCut", module = "symbolica.community.feynkit", frozen)]
pub struct PyDiagramCut {
    diagram: Py<PyFeynmanDiagram>,
    inner: DiagramCut,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyDiagramCut {
    /// Return the left amplitude and its generation metadata.
    #[getter]
    fn left(&self, py: Python<'_>) -> PyDiagramCutSide {
        PyDiagramCutSide {
            diagram: self.diagram.clone_ref(py),
            inner: self.inner.left.clone(),
        }
    }
    /// Return the right amplitude and its generation metadata.
    #[getter]
    fn right(&self, py: Python<'_>) -> PyDiagramCutSide {
        PyDiagramCutSide {
            diagram: self.diagram.clone_ref(py),
            inner: self.inner.right.clone(),
        }
    }
    /// Return the oriented crossing half-edges as a reusable selection.
    #[getter]
    fn subgraph(&self, py: Python<'_>) -> PyResult<Py<PySubgraph>> {
        self.diagram
            .borrow(py)
            .half_edge_selection(py, &self.inner.cut)
    }
    /// Return the crossing particle lines in their stored cut order.
    #[getter]
    fn edges(&self, py: Python<'_>) -> Vec<PyDiagramEdge> {
        let edges = self.diagram.borrow(py).whole().edges();
        self.inner
            .cut
            .iter()
            .map(|half| edges[half.edge.0].clone())
            .collect()
    }
    /// Return physical final-state particles in the stored cut-edge order.
    ///
    /// The finalized left side is the conjugate amplitude, so source-oriented
    /// lines correspond to antiparticles of their stored edge species.
    /// Momentum orientations retain the native positive-energy routing.
    #[getter]
    fn particles(&self, py: Python<'_>) -> PyResult<Vec<PyParticle>> {
        let diagram = self.diagram.borrow(py);
        let model = diagram.inner.model_arc();
        Ok(diagram
            .inner
            .cut_particles(&self.inner)
            .map_err(error::diagram)?
            .into_iter()
            .map(|particle| PyParticle::new(particle, Arc::clone(&model)))
            .collect())
    }

    /// Orient crossing lines from the left side to the right side.
    #[getter]
    fn orientations(&self) -> BTreeMap<usize, i32> {
        self.inner
            .cut
            .iter()
            .map(|half| {
                (
                    half.edge.0,
                    match half.endpoint {
                        DiagramEndpoint::Source => 1,
                        DiagramEndpoint::Target => -1,
                    },
                )
            })
            .collect()
    }
    /// Return stored edge momenta in the diagram's selected routing.
    /// Multiply by :attr:`orientations` for momenta directed from left to right.
    #[getter]
    fn momentum_signatures(&self, py: Python<'_>) -> BTreeMap<usize, PyMomentumSignature> {
        let diagram = self.diagram.borrow(py);
        self.inner
            .cut
            .iter()
            .map(|half| {
                (
                    half.edge.0,
                    diagram.inner.loop_momentum_basis().edge_signatures[&half.edge]
                        .clone()
                        .into(),
                )
            })
            .collect()
    }
    /// Construct the oriented cut-propagator distributions, separately from
    /// the remaining numerator, uncut propagators, and global factors.
    ///
    /// Examples
    /// --------
    /// >>> factors = diagram.cuts[0].propagators(edge_powers={2: 2})
    /// >>> distributions = [factor.to_expression() for factor in factors]
    ///
    /// Parameters
    /// ----------
    /// edge_powers : mapping[int, int] or None, optional
    ///     Positive propagator powers by diagram edge ID; omitted edges have power one.
    /// prescription : int, optional
    ///     Sign of the imaginary prescription, either +1 or -1.
    /// normalization : Expression or None, optional
    ///     Explicit replacement normalization; defaults to -2*pi*i.
    #[pyo3(signature = (*, edge_powers=None, prescription=1, normalization=None))]
    fn propagators(
        &self,
        edge_powers: Option<BTreeMap<usize, usize>>,
        prescription: i8,
        normalization: Option<ConvertibleToExpression>,
    ) -> PyResult<Vec<PyCutPropagator>> {
        let powers = edge_powers.unwrap_or_default();
        let normalization = normalization.map_or_else(
            || -2 * Atom::var(symbolica::atom::Symbol::PI) * Atom::i(),
            |value| value.to_expression().expr,
        );
        self.inner
            .cut
            .iter()
            .map(|half| {
                let edge = feynkit_cff::EdgeId::new(half.edge.0);
                let inner = feynkit_cff::CutPropagator {
                    energy: feynkit_cff::symbols::energy_atom(edge),
                    on_shell_energy: feynkit_cff::symbols::on_shell_atom(edge),
                    power: powers.get(&half.edge.0).copied().unwrap_or(1),
                    orientation: match half.endpoint {
                        DiagramEndpoint::Source => 1,
                        DiagramEndpoint::Target => -1,
                    },
                    prescription,
                    normalization: normalization.clone(),
                };
                inner.validate().map_err(error::cff)?;
                Ok(PyCutPropagator { inner })
            })
            .collect()
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.diagram)
    }
}

/// A topology threshold partition, independent of the requested physical final state.
///
/// Examples
/// --------
/// >>> threshold = diagram.topology_threshold_candidates[0]
/// >>> crossing_lines = threshold.edges
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "DiagramThresholdCandidate",
    module = "symbolica.community.feynkit",
    frozen
)]
pub struct PyDiagramThresholdCandidate {
    diagram: Py<PyFeynmanDiagram>,
    inner: DiagramThresholdCandidate,
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyDiagramThresholdCandidate {
    /// Return the left topology selection.
    #[getter]
    fn left(&self, py: Python<'_>) -> PyResult<Py<PySubgraph>> {
        self.diagram
            .borrow(py)
            .half_edge_selection(py, &self.inner.left)
    }
    /// Return the right topology selection.
    #[getter]
    fn right(&self, py: Python<'_>) -> PyResult<Py<PySubgraph>> {
        self.diagram
            .borrow(py)
            .half_edge_selection(py, &self.inner.right)
    }
    /// Return the oriented crossing half-edge selection.
    #[getter]
    fn subgraph(&self, py: Python<'_>) -> PyResult<Py<PySubgraph>> {
        self.diagram
            .borrow(py)
            .half_edge_selection(py, &self.inner.cut)
    }
    /// Return threshold-crossing lines in their stored order.
    #[getter]
    fn edges(&self, py: Python<'_>) -> Vec<PyDiagramEdge> {
        let edges = self.diagram.borrow(py).whole().edges();
        self.inner
            .cut
            .iter()
            .map(|half| edges[half.edge.0].clone())
            .collect()
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.diagram)
    }
}

/// Integer coefficients expressing one edge momentum in a chosen basis.
///
/// A signature separates coefficients of independent loop momenta from those
/// of external momenta, for example ``-l1 + p2``.
///
/// Examples
/// --------
/// >>> basis = next(iter(diagram.loop_momentum_bases(limit=1)))
/// >>> edge_id, signature = next(iter(basis.edge_signatures.items()))
/// >>> print(f"q_{edge_id} = {signature.format_momentum()}")
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

#[derive(FromPyObject)]
enum MomentumInput {
    Three(PyThreeMomentum),
    Four(PyFourMomentum),
}

#[derive(IntoPyObject)]
enum RoutedMomentum {
    Three(PyThreeMomentum),
    Four(PyFourMomentum),
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for MomentumInput {
    fn type_input() -> TypeInfo {
        PyThreeMomentum::type_input() | PyFourMomentum::type_input()
    }

    fn type_output() -> TypeInfo {
        PyThreeMomentum::type_output() | PyFourMomentum::type_output()
    }
}

#[cfg(feature = "python_stubgen")]
impl PyStubType for RoutedMomentum {
    fn type_input() -> TypeInfo {
        PyThreeMomentum::type_input() | PyFourMomentum::type_input()
    }

    fn type_output() -> TypeInfo {
        PyThreeMomentum::type_output() | PyFourMomentum::type_output()
    }
}

fn components_close<const N: usize>(left: [f64; N], right: [f64; N]) -> bool {
    left.into_iter().zip(right).all(|(left, right)| {
        let scale = 1.0 + left.abs().max(right.abs());
        (left - right).abs() <= 1.0e-10 * scale
    })
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyMomentumSignature {
    /// Return the integer coefficients of the independent loop momenta.
    #[getter]
    fn loops(&self) -> Vec<isize> {
        self.inner.loops.integer_coefficients()
    }

    /// Return the integer coefficients of the external momenta.
    #[getter]
    fn external(&self) -> Vec<isize> {
        self.inner.external.integer_coefficients()
    }

    /// Return loop and external coefficients as a pair of integer lists.
    ///
    /// Examples
    /// --------
    /// >>> signature = next(iter(basis.edge_signatures.values()))
    /// >>> loops, external = signature.integer_coefficients()
    /// >>> print("loop coefficients:", loops, "external coefficients:", external)
    ///
    fn integer_coefficients(&self) -> (Vec<isize>, Vec<isize>) {
        self.inner.integer_coefficients()
    }

    /// Format the signature as a conventional momentum sum.
    ///
    /// Examples
    /// --------
    /// Print the momentum routing assigned to every propagator:
    ///
    /// >>> for edge_id, signature in basis.edge_signatures.items():
    /// ...     print(edge_id, signature.format_momentum())
    ///
    fn format_momentum(&self) -> String {
        self.inner.format_momentum()
    }

    /// Return the conventional momentum sum for ``str(signature)``.
    ///
    /// Examples
    /// --------
    /// >>> signature = next(iter(basis.edge_signatures.values()))
    /// >>> print(f"propagator momentum: {signature}")
    ///
    fn __str__(&self) -> String {
        self.inner.format_momentum()
    }

    /// Return an unambiguous momentum-signature summary.
    ///
    /// Examples
    /// --------
    /// >>> edge_id, signature = next(iter(basis.edge_signatures.items()))
    /// >>> print(f"edge {edge_id}: {signature!r}")
    ///
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
    /// pretty : Any
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

/// A consistent routing of independent loop and external momenta.
///
/// Chords of the spanning tree define the loop momenta; every diagram edge is
/// then assigned an integer ``MomentumSignature`` by momentum conservation.
///
/// Examples
/// --------
/// >>> basis = next(iter(diagram.loop_momentum_bases(limit=1)))
/// >>> len(basis.loop_edges) == diagram.loop_count
/// True
/// >>> assignments = basis.edge_signatures
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
    owner: Arc<()>,
}

impl PyLoopMomentumBasis {
    fn from_diagram(inner: LoopMomentumBasis, diagram: &PyFeynmanDiagram) -> Self {
        Self {
            inner,
            owner: Arc::clone(&diagram.owner),
        }
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyLoopMomentumBasis {
    /// Return the edge identifiers belonging to the spanning tree.
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
    #[getter]
    fn loop_edges(&self) -> Vec<usize> {
        self.inner.loop_edges.iter().map(|id| id.0).collect()
    }

    /// Return the identifiers of edges attached to external states.
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
    #[getter]
    fn dependent_externals(&self) -> Vec<usize> {
        self.inner
            .dependent_externals
            .iter()
            .map(|id| id.0)
            .collect()
    }

    /// Return momentum signatures keyed by stable diagram edge ID.
    ///
    /// Examples
    /// --------
    /// >>> basis = diagram.loop_momentum_bases(limit=1)[0]
    /// >>> momentum_by_edge = {
    /// ...     edge_id: signature.format_momentum()
    /// ...     for edge_id, signature in basis.edge_signatures.items()
    /// ... }
    /// >>> for edge in diagram.edges:
    /// ...     print(edge.particle_name, momentum_by_edge[edge.id])
    ///
    #[getter]
    fn edge_signatures(&self) -> BTreeMap<usize, PyMomentumSignature> {
        self.inner
            .edge_signatures
            .iter()
            .map(|(edge, signature)| (edge.0, signature.clone().into()))
            .collect()
    }

    /// Return Symbolica Replacement objects using the same routing as numerical momenta.
    /// Tensor index arguments are retained on every routed vector.
    ///
    /// Examples
    /// --------
    /// >>> rules = diagram.loop_momentum_basis.momentum_replacements()
    /// >>> routed = diagram.numerator_expression().replace_multiple(rules)
    ///
    fn momentum_replacements(&self, py: Python<'_>) -> PyResult<Vec<Py<PyAny>>> {
        let constructor = py.import("symbolica.core")?.getattr("Replacement")?;
        self.inner
            .momentum_replacements()
            .into_iter()
            .map(|rule| {
                let pattern = rule.pat.to_atom().map_err(error::DiagramError::new_err)?;
                let symbolica::id::ReplaceWith::Pattern(rhs) = rule.rhs else {
                    unreachable!("momentum substitutions are symbolic patterns")
                };
                let rhs = rhs.to_atom().map_err(error::DiagramError::new_err)?;
                Ok(constructor
                    .call1((
                        PythonExpression { expr: pattern },
                        PythonExpression { expr: rhs },
                    ))?
                    .unbind())
            })
            .collect()
    }

    /// Express edge momenta in this basis while retaining any tensor index arguments.
    ///
    /// Examples
    /// --------
    /// >>> routed = diagram.loop_momentum_basis.route_expression(diagram.numerator_expression())
    ///
    /// Parameters
    /// ----------
    /// expression : Expression or TensorExpression
    ///     Expression with canonical indexed edge momenta to route.
    fn route_expression(&self, expression: ConvertibleToExpression) -> PythonExpression {
        PythonExpression {
            expr: self
                .inner
                .route_expression(&expression.to_expression().expr),
        }
    }

    /// Route loop and external momenta through every diagram edge.
    ///
    /// Inputs must be uniformly :class:`ThreeMomentum` or uniformly
    /// :class:`FourMomentum`. Loop momenta follow :attr:`loop_edges`; external
    /// momenta follow :attr:`external_edges`. The result is keyed by stable edge
    /// ID. FeynKit checks the supplied external momenta against momentum
    /// conservation, including each dependent external leg.
    ///
    /// Examples
    /// --------
    /// Route a one-loop spatial momentum through the complete graph:
    ///
    /// >>> routed = basis.route([loop_momentum], external_spatial_momenta)
    /// >>> internal_momentum = routed[basis.loop_edges[0]]
    ///
    /// Parameters
    /// ----------
    /// loop_momenta : sequence[ThreeMomentum] or sequence[FourMomentum]
    ///     Independent loop momenta in :attr:`loop_edges` order.
    /// external_momenta : sequence[ThreeMomentum] or sequence[FourMomentum]
    ///     External momenta in :attr:`external_edges` order.
    #[gen_stub(skip)]
    fn route(
        &self,
        loop_momenta: Vec<MomentumInput>,
        external_momenta: Vec<MomentumInput>,
    ) -> PyResult<BTreeMap<usize, RoutedMomentum>> {
        let use_four_momenta = loop_momenta
            .iter()
            .chain(&external_momenta)
            .next()
            .map(|momentum| matches!(momentum, MomentumInput::Four(_)))
            .ok_or_else(|| {
                PyTypeError::new_err(
                    "cannot infer momentum dimension when both basis inputs are empty",
                )
            })?;

        if use_four_momenta {
            let loops = loop_momenta
                .into_iter()
                .map(|momentum| match momentum {
                    MomentumInput::Four(momentum) => Ok(momentum.inner),
                    MomentumInput::Three(_) => Err(PyTypeError::new_err(
                        "loop_momenta and external_momenta must use one momentum dimension",
                    )),
                })
                .collect::<PyResult<Vec<_>>>()?;
            let external = external_momenta
                .into_iter()
                .map(|momentum| match momentum {
                    MomentumInput::Four(momentum) => Ok(momentum.inner),
                    MomentumInput::Three(_) => Err(PyTypeError::new_err(
                        "loop_momenta and external_momenta must use one momentum dimension",
                    )),
                })
                .collect::<PyResult<Vec<_>>>()?;
            let mut routed = BTreeMap::new();
            for (edge, signature) in &self.inner.edge_signatures {
                let momentum = signature
                    .apply(&loops, &external)
                    .map_err(|error| error::DiagramError::new_err(error.to_string()))?
                    .ok_or_else(|| {
                        error::DiagramError::new_err(format!(
                            "edge {} has an all-zero momentum signature",
                            edge.0
                        ))
                    })?;
                if let Some(index) = self
                    .inner
                    .external_edges
                    .iter()
                    .position(|external_edge| external_edge == edge)
                    && !components_close(<[f64; 4]>::from(momentum), external[index].into())
                {
                    return Err(error::DiagramError::new_err(format!(
                        "external momentum at edge {} violates the basis momentum-conservation relation",
                        edge.0
                    )));
                }
                routed.insert(edge.0, RoutedMomentum::Four(momentum.into()));
            }
            Ok(routed)
        } else {
            let loops = loop_momenta
                .into_iter()
                .map(|momentum| match momentum {
                    MomentumInput::Three(momentum) => Ok(momentum.inner),
                    MomentumInput::Four(_) => Err(PyTypeError::new_err(
                        "loop_momenta and external_momenta must use one momentum dimension",
                    )),
                })
                .collect::<PyResult<Vec<_>>>()?;
            let external = external_momenta
                .into_iter()
                .map(|momentum| match momentum {
                    MomentumInput::Three(momentum) => Ok(momentum.inner),
                    MomentumInput::Four(_) => Err(PyTypeError::new_err(
                        "loop_momenta and external_momenta must use one momentum dimension",
                    )),
                })
                .collect::<PyResult<Vec<_>>>()?;
            let mut routed = BTreeMap::new();
            for (edge, signature) in &self.inner.edge_signatures {
                let momentum = signature
                    .apply(&loops, &external)
                    .map_err(|error| error::DiagramError::new_err(error.to_string()))?
                    .ok_or_else(|| {
                        error::DiagramError::new_err(format!(
                            "edge {} has an all-zero momentum signature",
                            edge.0
                        ))
                    })?;
                if let Some(index) = self
                    .inner
                    .external_edges
                    .iter()
                    .position(|external_edge| external_edge == edge)
                    && !components_close(<[f64; 3]>::from(momentum), external[index].into())
                {
                    return Err(error::DiagramError::new_err(format!(
                        "external momentum at edge {} violates the basis momentum-conservation relation",
                        edge.0
                    )));
                }
                routed.insert(edge.0, RoutedMomentum::Three(momentum.into()));
            }
            Ok(routed)
        }
    }

    /// Return a concise description of the selected momentum basis.
    ///
    /// Examples
    /// --------
    /// >>> basis = diagram.loop_momentum_bases(limit=1)[0]
    /// >>> print(basis)
    ///
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
    /// Leave ``basis`` as the final expression in a Jupyter or Marimo cell to
    /// inspect every propagator's loop- and external-momentum assignment.
    ///
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
    /// pretty : Any
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

#[cfg(feature = "python_stubgen")]
submit! {
    gen_methods_from_python! {
        r#"
        import typing

        class PyLoopMomentumBasis:
            @typing.overload
            def route(
                self,
                loop_momenta: typing.Sequence[ThreeMomentum],
                external_momenta: typing.Sequence[ThreeMomentum],
            ) -> dict[int, ThreeMomentum]:
                """Route three-momenta through every diagram edge.

                Examples
                --------
                >>> routed = basis.route([loop_momentum], external_spatial_momenta)
                >>> internal_momentum = routed[basis.loop_edges[0]]

                Parameters
                ----------
                loop_momenta : sequence[ThreeMomentum]
                    Independent loop momenta in ``basis.loop_edges`` order.
                external_momenta : sequence[ThreeMomentum]
                    External momenta in ``basis.external_edges`` order.
                """

            @typing.overload
            def route(
                self,
                loop_momenta: typing.Sequence[FourMomentum],
                external_momenta: typing.Sequence[FourMomentum],
            ) -> dict[int, FourMomentum]:
                """Route four-momenta through every diagram edge.

                Examples
                --------
                >>> routed = basis.route([loop_momentum], external_four_momenta)
                >>> internal_momentum = routed[basis.loop_edges[0]]

                Parameters
                ----------
                loop_momenta : sequence[FourMomentum]
                    Independent loop momenta in ``basis.loop_edges`` order.
                external_momenta : sequence[FourMomentum]
                    External momenta in ``basis.external_edges`` order.
                """
        "#
    }
}

/// A typed Feynman graph with model and symbolic physics annotations.
///
/// Diagrams expose vertices, propagator edges, loop-momentum routings, symmetry
/// factors, Symbolica expressions, and Linnest/Typst notebook rendering.
///
/// Examples
/// --------
/// >>> diagram = next(iter(result.diagrams))
/// >>> diagram.validate()
/// >>> diagram  # renders as a Linnest graph in Jupyter or Marimo
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(
    name = "FeynmanDiagram",
    module = "symbolica.community.feynkit",
    subclass,
    frozen,
    from_py_object
)]
#[derive(Clone)]
pub struct PyFeynmanDiagram {
    pub(crate) inner: Arc<FeynmanDiagram>,
    linnet: LinnetCache,
    owner: Arc<()>,
    selected_region: Option<DiagramSelection>,
}

impl From<FeynmanDiagram> for PyFeynmanDiagram {
    fn from(inner: FeynmanDiagram) -> Self {
        Self {
            inner: Arc::new(inner),
            linnet: LinnetCache::default(),
            owner: Arc::new(()),
            selected_region: None,
        }
    }
}

#[derive(Clone)]
struct DiagramSelection {
    hedges: SuBitGraph,
    isolated: BTreeSet<usize>,
}

impl PyFeynmanDiagram {
    pub(crate) fn whole(&self) -> Self {
        let mut whole = self.clone();
        whole.selected_region = None;
        whole
    }

    pub(crate) fn selection(&self) -> SuBitGraph {
        self.selected_region.as_ref().map_or_else(
            || self.inner.underlying().full_filter(),
            |region| region.hedges.clone(),
        )
    }

    fn isolated_nodes(&self) -> BTreeSet<usize> {
        self.selected_region.as_ref().map_or_else(
            || {
                self.inner
                    .underlying()
                    .iter_nodes()
                    .filter_map(|(node, mut crown, _)| crown.next().is_none().then_some(node.0))
                    .collect()
            },
            |region| region.isolated.clone(),
        )
    }

    fn region_argument(
        &self,
        py: Python<'_>,
        value: &Bound<'_, PyAny>,
    ) -> PyResult<DiagramSelection> {
        if let Ok(diagram) = value.extract::<PyRef<'_, PyFeynmanDiagram>>() {
            if !Arc::ptr_eq(&self.inner, &diagram.inner) {
                return Err(PyValueError::new_err(
                    "subgraph belongs to a different original diagram",
                ));
            }
            return Ok(DiagramSelection {
                hedges: diagram.selection(),
                isolated: diagram.isolated_nodes(),
            });
        }
        Ok(DiagramSelection {
            hedges: self.linnet.selection(py, self, Some(value))?,
            isolated: value
                .call_method0("isolated_node_indices")?
                .extract::<Vec<usize>>()?
                .into_iter()
                .collect(),
        })
    }

    pub(crate) fn selection_argument(
        &self,
        py: Python<'_>,
        value: &Bound<'_, PyAny>,
    ) -> PyResult<SuBitGraph> {
        Ok(self.region_argument(py, value)?.hedges)
    }

    fn view(&self, py: Python<'_>, region: DiagramSelection) -> PyResult<Py<PySubgraph>> {
        // Share the canonical export holder with the original diagram.
        let mut base = self.clone();
        base.selected_region = Some(region);
        Py::new(py, (PySubgraph {}, base))
    }

    fn wrap_selection(
        &self,
        py: Python<'_>,
        selection: &Bound<'_, PyAny>,
    ) -> PyResult<Py<PySubgraph>> {
        self.view(py, self.region_argument(py, selection)?)
    }

    fn restrict_selection(
        &self,
        py: Python<'_>,
        selection: &Bound<'_, PyAny>,
    ) -> PyResult<Py<PySubgraph>> {
        let region = self.region_argument(py, selection)?;
        self.view(
            py,
            DiagramSelection {
                hedges: region.hedges.intersection(&self.selection()),
                isolated: region
                    .isolated
                    .intersection(&self.isolated_nodes())
                    .copied()
                    .collect(),
            },
        )
    }

    fn half_edge_selection(
        &self,
        py: Python<'_>,
        halves: &[DiagramHalfEdge],
    ) -> PyResult<Py<PySubgraph>> {
        let graph = self.inner.underlying();
        let mut selected = SuBitGraph::empty(graph.n_hedges());
        for half in halves {
            selected.add(self.inner.half_edge_id(*half).ok_or_else(|| {
                error::DiagramError::new_err("cut references a missing half-edge")
            })?);
        }
        self.view(
            py,
            DiagramSelection {
                hedges: selected,
                isolated: BTreeSet::new(),
            },
        )
    }

    fn local_numerator(&self, without: Option<&DiagramSelection>) -> Atom {
        let graph = self.inner.underlying();
        let excluded = without.map_or_else(
            || SuBitGraph::empty(graph.n_hedges()),
            |region| region.hedges.clone(),
        );
        let isolated = self.isolated_nodes();
        graph.numerator_of(
            &self.selection(),
            &excluded,
            |node, vertex| {
                if graph.iter_crown(node).next().is_none()
                    && (!isolated.contains(&node.0)
                        || without.is_some_and(|region| region.isolated.contains(&node.0)))
                {
                    Atom::one()
                } else {
                    vertex.numerator.clone()
                }
            },
            |edge| edge.numerator.clone(),
            |edge| edge.is_dummy,
        )
    }

    fn is_complete(&self) -> bool {
        self.selected_region.as_ref().is_none_or(|region| {
            region.hedges == self.inner.underlying().full_filter()
                && self
                    .inner
                    .underlying()
                    .iter_nodes()
                    .all(|(node, mut crown, _)| {
                        crown.next().is_some() || region.isolated.contains(&node.0)
                    })
        })
    }

    fn require_complete(&self) -> PyResult<()> {
        if !self.is_complete() {
            return Err(error::DiagramError::new_err(
                "this operation needs an independent diagram; call subgraph.excise() first",
            ));
        }
        Ok(())
    }

    fn includes_vertex(&self, vertex: usize) -> bool {
        self.selected_region.as_ref().is_none_or(|region| {
            region.isolated.contains(&vertex)
                || self
                    .inner
                    .underlying()
                    .iter_crown(NodeIndex(vertex))
                    .any(|hedge| region.hedges.includes(&hedge))
        })
    }

    fn includes_edge(&self, edge: feynkit_graph::EdgeId) -> bool {
        self.selected_region.as_ref().is_none_or(|region| {
            [DiagramEndpoint::Source, DiagramEndpoint::Target]
                .into_iter()
                .any(|endpoint| {
                    self.inner
                        .half_edge_id(DiagramHalfEdge { edge, endpoint })
                        .is_some_and(|hedge| region.hedges.includes(&hedge))
                })
        })
    }
}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyFeynmanDiagram {
    /// Return the canonical installed Linnet graph with physics objects as payloads.
    ///
    /// Node, edge, and half-edge identities are mapped explicitly. Structural
    /// edits affect this analysis graph only; the next export starts a fresh
    /// graph, and selections from the modified topology cannot be used here.
    ///
    /// Examples
    /// --------
    /// >>> graph = diagram.to_linnet()
    /// >>> gluons = graph.filter(edge=lambda edge: edge.data.particle_name == "g")
    ///
    #[gen_stub(override_return_type(type_repr="linnet.Graph", imports=("linnet")))]
    fn to_linnet(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.linnet.graph(py, self)
    }

    /// Return the canonical Linnet selection representing this physics region.
    #[getter]
    #[gen_stub(override_return_type(type_repr="linnet.Subgraph", imports=("linnet")))]
    fn linnet_selection(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        self.linnet
            .export_selection(py, self, &self.selection(), &self.isolated_nodes())
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        self.linnet.traverse(visit)
    }

    #[gen_stub(skip)]
    fn __clear__(&self) {
        self.linnet.clear();
    }

    /// Select graph elements using canonical Linnet IDs, or import a graph-bound selection.
    /// Nested selections intersect this region and retain its immutable original diagram.
    ///
    /// Examples
    /// --------
    /// >>> region = diagram.subgraph(edges=[0, 1])
    /// >>> numerator = region.numerator_expression()
    ///
    /// Parameters
    /// ----------
    /// selection : Subgraph or linnet.Subgraph or None, optional
    ///     Existing selection from the same original diagram; exclusive with element IDs.
    /// nodes : list[int] or None, optional
    ///     Canonical Linnet nodes IDs to include.
    /// edges : list[int] or None, optional
    ///     Canonical Linnet edges IDs to include.
    /// half_edges : list[int] or None, optional
    ///     Canonical Linnet half-edges IDs to include.
    #[pyo3(signature = (selection=None, *, nodes=None, edges=None, half_edges=None))]
    fn subgraph(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="Subgraph | linnet.Subgraph | None", imports=("linnet")))]
        selection: Option<&Bound<'_, PyAny>>,
        nodes: Option<Vec<usize>>,
        edges: Option<Vec<usize>>,
        half_edges: Option<Vec<usize>>,
    ) -> PyResult<Py<PySubgraph>> {
        if let Some(selection) = selection {
            if nodes.is_some() || edges.is_some() || half_edges.is_some() {
                return Err(PyTypeError::new_err(
                    "pass either a selection or element IDs, not both",
                ));
            }
            return self.restrict_selection(py, selection);
        }
        let graph = self.to_linnet(py)?;
        let kwargs = PyDict::new(py);
        for (name, indices) in [
            ("nodes", nodes),
            ("edges", edges),
            ("half_edges", half_edges),
        ] {
            if let Some(indices) = indices {
                kwargs.set_item(name, indices)?;
            }
        }
        let selected = graph.bind(py).call_method("subgraph", (), Some(&kwargs))?;
        self.restrict_selection(py, &selected)
    }

    /// Select by predicates on Linnet views; their ``data`` is a physics object.
    ///
    /// Examples
    /// --------
    /// >>> region = diagram.filter(edge=lambda edge: edge.data.particle_name == "g")
    ///
    /// Parameters
    /// ----------
    /// node : callable or None, optional
    ///     Predicate on canonical Linnet node views.
    /// edge : callable or None, optional
    ///     Predicate on canonical Linnet edge views.
    /// half_edge : callable or None, optional
    ///     Predicate on canonical Linnet half-edge views.
    #[pyo3(signature = (*, node=None, edge=None, half_edge=None))]
    fn filter(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Callable[[linnet.Node], bool] | None", imports=("typing", "linnet")))]
        node: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[linnet.Edge], bool] | None", imports=("typing", "linnet")))]
        edge: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[linnet.HalfEdge], bool] | None", imports=("typing", "linnet")))]
        half_edge: Option<Py<PyAny>>,
    ) -> PyResult<Py<PySubgraph>> {
        let graph = self.to_linnet(py)?;
        let kwargs = PyDict::new(py);
        kwargs.set_item("node", node)?;
        kwargs.set_item("edge", edge)?;
        kwargs.set_item("half_edge", half_edge)?;
        let selected = graph.bind(py).call_method("filter", (), Some(&kwargs))?;
        self.restrict_selection(py, &selected)
    }

    /// Return boundaries around the selected interaction region.
    ///
    /// The current region determines which interaction boundaries are returned.
    ///
    /// Examples
    /// --------
    /// >>> region = diagram.subgraph(nodes=[0])
    /// >>> boundary = region.boundary()
    fn boundary(&self, py: Python<'_>) -> PyResult<Py<PySubgraph>> {
        let selected = self.linnet_selection(py)?;
        let result = self
            .to_linnet(py)?
            .bind(py)
            .call_method1("boundary", (selected,))?;
        self.wrap_selection(py, &result)
    }

    /// Return connected interaction regions as reusable selections.
    ///
    /// Examples
    /// --------
    /// >>> components = diagram.connected_components()
    fn connected_components(&self, py: Python<'_>) -> PyResult<Vec<Py<PySubgraph>>> {
        let selected = self.linnet_selection(py)?;
        self.to_linnet(py)?
            .bind(py)
            .call_method1("connected_components", (selected,))?
            .try_iter()?
            .map(|value| self.wrap_selection(py, &value?))
            .collect()
    }

    /// Test connectivity of the current diagram or selected region.
    ///
    /// Examples
    /// --------
    /// >>> connected = diagram.is_connected()
    fn is_connected(&self, py: Python<'_>) -> PyResult<bool> {
        let selected = self.linnet_selection(py)?;
        self.to_linnet(py)?
            .bind(py)
            .call_method1("is_connected", (selected,))?
            .extract()
    }

    /// Return lines whose removal disconnects the selection.
    ///
    /// Examples
    /// --------
    /// >>> bridges = diagram.bridges()
    fn bridges(&self, py: Python<'_>) -> PyResult<Py<PySubgraph>> {
        let selected = self.linnet_selection(py)?;
        let result = self
            .to_linnet(py)?
            .bind(py)
            .call_method1("bridges", (selected,))?;
        self.wrap_selection(py, &result)
    }

    /// Return a cycle basis and its covered half-edges.
    ///
    /// Examples
    /// --------
    /// >>> cycles, covered = diagram.cycle_basis()
    #[gen_stub(override_return_type(type_repr="tuple[list[linnet.Cycle], Subgraph]", imports=("linnet")))]
    fn cycle_basis(&self, py: Python<'_>) -> PyResult<(Py<PyAny>, Py<PySubgraph>)> {
        let selected = self.linnet_selection(py)?;
        let result = self
            .to_linnet(py)?
            .bind(py)
            .call_method1("cycle_basis", (selected,))?;
        let result = result.cast::<PyTuple>()?;
        Ok((
            result.get_item(0)?.unbind(),
            self.wrap_selection(py, &result.get_item(1)?)?,
        ))
    }

    /// Enumerate spanning forests within the selected topology.
    ///
    /// Examples
    /// --------
    /// >>> forests = diagram.all_spanning_forests()
    fn all_spanning_forests(&self, py: Python<'_>) -> PyResult<Vec<Py<PySubgraph>>> {
        let selected = self.linnet_selection(py)?;
        self.to_linnet(py)?
            .bind(py)
            .call_method1("all_spanning_forests", (selected,))?
            .try_iter()?
            .map(|value| self.wrap_selection(py, &value?))
            .collect()
    }

    /// Enumerate minimal cutsets, independently of physical final-state cuts.
    ///
    /// Examples
    /// --------
    /// >>> bonds = diagram.all_bonds(min_size=2, max_size=3)
    ///
    /// Parameters
    /// ----------
    /// min_size : int or None, optional
    ///     Minimum number of crossing edges.
    /// max_size : int or None, optional
    ///     Maximum number of crossing edges.
    #[pyo3(signature = (*, min_size=None, max_size=None))]
    fn all_bonds(
        &self,
        py: Python<'_>,
        min_size: Option<usize>,
        max_size: Option<usize>,
    ) -> PyResult<Vec<Py<PySubgraph>>> {
        let kwargs = PyDict::new(py);
        kwargs.set_item("subgraph", self.linnet_selection(py)?)?;
        kwargs.set_item("min_size", min_size)?;
        kwargs.set_item("max_size", max_size)?;
        self.to_linnet(py)?
            .bind(py)
            .call_method("all_bonds", (), Some(&kwargs))?
            .try_iter()?
            .map(|value| self.wrap_selection(py, &value?))
            .collect()
    }

    /// Enumerate separating partitions between disjoint interaction vertex groups.
    ///
    /// Examples
    /// --------
    /// >>> partitions = diagram.all_cuts([0], [1])
    ///
    /// Parameters
    /// ----------
    /// source : list[int]
    ///     Canonical Linnet vertices required on the first side.
    /// target : list[int]
    ///     Canonical Linnet vertices required on the opposite side.
    #[gen_stub(override_return_type(type_repr="list[linnet.CutPartition]", imports=("linnet")))]
    fn all_cuts(
        &self,
        py: Python<'_>,
        source: Vec<usize>,
        target: Vec<usize>,
    ) -> PyResult<Py<PyAny>> {
        self.require_complete()?;
        Ok(self
            .to_linnet(py)?
            .bind(py)
            .call_method1("all_cuts", (source, target))?
            .unbind())
    }

    /// Traverse a selected interaction region in depth-first order.
    ///
    /// Examples
    /// --------
    /// >>> tree = diagram.depth_first_traverse(0)
    ///
    /// Parameters
    /// ----------
    /// root : int
    ///     Canonical Linnet vertex ID at which traversal starts.
    /// include : int or None, optional
    ///     Canonical Linnet half-edge ID to prioritize at the root.
    #[pyo3(signature = (root, *, include=None))]
    #[gen_stub(override_return_type(type_repr="linnet.TraversalTree", imports=("linnet")))]
    fn depth_first_traverse(
        &self,
        py: Python<'_>,
        root: usize,
        include: Option<usize>,
    ) -> PyResult<Py<PyAny>> {
        let kwargs = PyDict::new(py);
        kwargs.set_item("subgraph", self.linnet_selection(py)?)?;
        kwargs.set_item("include", include)?;
        Ok(self
            .to_linnet(py)?
            .bind(py)
            .call_method("depth_first_traverse", (root,), Some(&kwargs))?
            .unbind())
    }

    /// Traverse a selected interaction region in breadth-first order.
    ///
    /// Examples
    /// --------
    /// >>> tree = diagram.breadth_first_traverse(0)
    ///
    /// Parameters
    /// ----------
    /// root : int
    ///     Canonical Linnet vertex ID at which traversal starts.
    /// include : int or None, optional
    ///     Canonical Linnet half-edge ID to prioritize at the root.
    #[pyo3(signature = (root, *, include=None))]
    #[gen_stub(override_return_type(type_repr="linnet.TraversalTree", imports=("linnet")))]
    fn breadth_first_traverse(
        &self,
        py: Python<'_>,
        root: usize,
        include: Option<usize>,
    ) -> PyResult<Py<PyAny>> {
        let kwargs = PyDict::new(py);
        kwargs.set_item("subgraph", self.linnet_selection(py)?)?;
        kwargs.set_item("include", include)?;
        Ok(self
            .to_linnet(py)?
            .bind(py)
            .call_method("breadth_first_traverse", (root,), Some(&kwargs))?
            .unbind())
    }

    /// Deserialize a Feynman diagram from its JSON representation.
    ///
    /// Examples
    /// --------
    /// >>> encoded = diagram.to_json()
    /// >>> restored = FeynmanDiagram.from_json(model, encoded)
    /// >>> restored.validate()
    /// >>> restored  # render the recovered graph in a notebook
    ///
    /// Parameters
    /// ----------
    /// model : Model
    ///     Model whose stable IDs and fingerprint are referenced by the diagram.
    /// json : str
    ///     JSON text produced by :meth:`FeynmanDiagram.to_json` or another
    ///     schema-compatible producer.
    #[staticmethod]
    fn from_json(model: &PyModel, json: &str) -> PyResult<Self> {
        FeynmanDiagram::from_json(Arc::clone(&model.inner), json)
            .map(Into::into)
            .map_err(error::diagram)
    }

    /// Parse a Feynman diagram from Graphviz DOT text.
    ///
    /// Examples
    /// --------
    /// >>> dot = diagram.to_dot()
    /// >>> restored = FeynmanDiagram.from_dot(model, dot)
    /// >>> restored.validate()
    /// >>> restored  # preserve topology through a Graphviz workflow
    ///
    /// Parameters
    /// ----------
    /// model : Model
    ///     Model used to resolve particles and interactions or validate annotated IDs.
    /// dot : str
    ///     Compact physics DOT or annotated FeynKit DOT. Compact cross-sections pair
    ///     initial-state legs with is_cut and specify comma-separated final_state particles.
    #[staticmethod]
    fn from_dot(model: &PyModel, dot: &str) -> PyResult<Self> {
        FeynmanDiagram::from_dot(Arc::clone(&model.inner), dot)
            .map(Into::into)
            .map_err(error::diagram)
    }

    /// Return the deterministic name assigned during diagram generation.
    #[getter]
    fn name(&self) -> &str {
        self.inner.name()
    }

    /// Return the stable content-derived hexadecimal diagram ID.
    #[getter]
    fn id(&self) -> String {
        self.inner.id().to_string()
    }

    /// Return the positive integer denominator of the graph symmetry factor.
    ///
    /// Examples
    /// --------
    /// >>> graph_weight = 1 / diagram.symmetry_factor
    ///
    #[getter]
    fn symmetry_factor(&self) -> u64 {
        self.inner.symmetry_factor()
    }

    /// Return the diagram-wide multiplicative factor as Symbolica source text.
    ///
    /// This string form is useful for serialization. For algebra or notebook
    /// display, prefer ``overall_factor_expression()`` so Symbolica's native
    /// rich formatting is preserved.
    ///
    /// Examples
    /// --------
    /// >>> source = diagram.overall_factor
    /// >>> factor = diagram.overall_factor_expression()
    /// >>> factor  # rich Symbolica output in a notebook
    ///
    #[getter]
    fn overall_factor(&self) -> String {
        self.inner.overall_factor().to_plain_string()
    }

    /// Return the diagram-wide multiplicative factor as a Symbolica expression.
    ///
    /// ``evaluate=True`` evaluates the generator's sign, multiplicity and
    /// symmetry annotations using the shared graph-factor evaluator. Other
    /// symbolic factors remain unchanged.
    ///
    /// Examples
    /// --------
    /// >>> factor = diagram.overall_factor_expression()
    /// >>> weighted_numerator = factor * diagram.numerator_expression()
    /// >>> weighted_numerator
    ///
    /// Parameters
    /// ----------
    /// evaluate : bool
    ///     Evaluate known graph-factor annotations while preserving other symbols.
    #[pyo3(signature = (*, evaluate=false))]
    fn overall_factor_expression(&self, evaluate: bool) -> PythonExpression {
        PythonExpression {
            expr: if evaluate {
                feynkit_graph::expressions::evaluate_overall_factor(
                    self.inner.overall_factor().as_view(),
                )
            } else {
                self.inner.overall_factor().clone()
            },
        }
    }

    /// Return the diagram numerator annotation as source text.
    ///
    #[getter]
    fn numerator(&self) -> String {
        self.local_numerator(None).to_plain_string()
    }

    /// Return the product of internal propagator denominators as a scalar TensorExpression.
    ///
    /// Each factor retains GammaLoop's four-argument ``denom`` annotation around
    /// q_e² - m_e², with matching momentum labels and symbolic masses. Dimension
    /// defaults to GammaLoop's symbolic dimension; powers may be signed.
    /// The default region excludes external carriers and dummy edges. An explicit
    /// region follows GammaLoop and includes every complete selected edge. Widths,
    /// an imaginary prescription, and custom UFO denominator formulas are excluded.
    /// A selection without internal edges returns one.
    ///
    /// Examples
    /// --------
    /// >>> denominator = diagram.denominator_expression(in_lmb=True)
    /// >>> integrand = diagram.numerator_expression(in_lmb=True) / denominator
    /// >>> basis = diagram.loop_momentum_bases()[0]
    /// >>> denominator = diagram.denominator_expression(lmb=basis)
    ///
    /// Parameters
    /// ----------
    /// A complete diagram defaults to internal propagators; a Subgraph uses its region.
    /// edge_powers : mapping[int, int] or None, optional
    ///     Signed propagator powers by diagram edge ID; omitted edges have power one.
    /// dimension : Expression or int or None, optional
    ///     Lorentz dimension; defaults to the shared symbolic dimension.
    /// in_lmb : bool, optional
    ///     Express edge momenta in the diagram's stored loop-momentum basis.
    /// lmb : LoopMomentumBasis or None, optional
    ///     Basis from this diagram instance. Supplying it enables routing and
    ///     takes precedence over ``in_lmb``, including for a selected region.
    #[pyo3(signature = (*, edge_powers=None, dimension=None, in_lmb=false, lmb=None))]
    fn denominator_expression(
        &self,
        py: Python<'_>,
        edge_powers: Option<BTreeMap<usize, isize>>,
        dimension: Option<ConvertibleToExpression>,
        in_lmb: bool,
        lmb: Option<&PyLoopMomentumBasis>,
    ) -> PyResult<Py<TensorExpression>> {
        let selected = self.selected_region.as_ref().map_or_else(
            || self.inner.internal_subgraph(),
            |region| region.hedges.clone(),
        );
        let powers = edge_powers
            .unwrap_or_default()
            .into_iter()
            .map(|(edge, power)| (feynkit_graph::EdgeId(edge), power))
            .collect();
        let denominator = match dimension {
            Some(dimension) => {
                let dimension = dimension.to_expression();
                let dimension = dimension.expr.as_view().try_into().map_err(|_| {
                    PyTypeError::new_err("dimension must be a positive integer or a symbol")
                })?;
                self.inner
                    .denominator_of_in_dimension(&selected, &powers, dimension)
            }
            None => self.inner.denominator_of(&selected, &powers),
        }
        .map_err(error::diagram)?;
        let denominator = match lmb {
            Some(basis) => {
                if !Arc::ptr_eq(&basis.owner, &self.owner) {
                    return Err(error::DiagramError::new_err(
                        "momentum basis belongs to a different diagram",
                    ));
                }
                basis.inner.route_expression(&denominator)
            }
            None if in_lmb => self
                .inner
                .loop_momentum_basis()
                .route_expression(&denominator),
            None => denominator,
        };
        TensorExpression::from_atom_interface(py, denominator, None)
    }

    /// Build a scalar integral family using the diagram's stored momentum routing.
    ///
    /// Reuses the shared propagator builder and model masses. Denominators follow
    /// ascending internal edge IDs, as in ``internal_edges``, retaining bridges
    /// and repeated propagators. Dependent external coordinates are eliminated;
    /// external carriers and dummy edges are excluded. Widths, prescriptions and
    /// custom UFO denominator formulas are not inferred. Tree diagrams raise
    /// DiagramError because they contain no loop integral.
    ///
    /// Examples
    /// --------
    /// >>> family = diagram.integral_family()
    /// >>> p = family.external_momenta[0]
    /// >>> kin = family.kinematics.with_scalar_product(p, p, s)
    /// >>> family = diagram.integral_family(kinematics=kin)
    /// >>> U, F = family.symanzik(parameters)
    ///
    /// Parameters
    /// ----------
    /// kinematics : Kinematics or None, optional
    ///     Assumptions on routed momentum names and the Lorentz dimension.
    ///     None uses the shared symbolic dimension with no on-shell assumptions.
    #[pyo3(signature = (*, kinematics=None))]
    fn integral_family(&self, kinematics: Option<&PyKinematics>) -> PyResult<PyIntegralFamily> {
        self.require_complete()?;
        let default;
        let kinematics = match kinematics {
            Some(context) => &context.inner,
            None => {
                default = feynkit_kinematics::Kinematics::in_dimension(&Atom::var(
                    feynkit_graph::symbols::dimension(),
                ))
                .expect("the shared dimension is a Symbolica symbol");
                &default
            }
        };
        Ok(PyIntegralFamily {
            inner: self
                .inner
                .integral_family(kinematics)
                .map_err(error::diagram)?,
        })
    }

    /// Return the diagram numerator as a Spenso TensorExpression.
    ///
    /// Examples
    /// --------
    /// >>> numerator = diagram.numerator_expression(in_lmb=True)
    /// >>> basis = diagram.loop_momentum_bases()[0]
    /// >>> numerator = diagram.numerator_expression(lmb=basis)
    /// >>> integrand_numerator = diagram.overall_factor_expression() * numerator
    /// >>> integrand_numerator  # native Symbolica algebra and rich display
    ///
    /// Parameters
    /// ----------
    /// without : Subgraph or linnet.Subgraph or None, optional
    ///     Ignored region, using GammaLoop boundary and local-factor selection semantics.
    /// in_lmb : bool, optional
    ///     Express edge momenta in the diagram's stored loop-momentum basis.
    /// lmb : LoopMomentumBasis or None, optional
    ///     Basis from this diagram instance. Supplying it enables routing and
    ///     takes precedence over ``in_lmb``, including for a selected region.
    #[pyo3(signature = (*, without=None, in_lmb=false, lmb=None))]
    fn numerator_expression(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="Subgraph | linnet.Subgraph | None", imports=("linnet")))]
        without: Option<&Bound<'_, PyAny>>,
        in_lmb: bool,
        lmb: Option<&PyLoopMomentumBasis>,
    ) -> PyResult<Py<TensorExpression>> {
        let without = without
            .map(|region| self.region_argument(py, region))
            .transpose()?;
        let numerator = self.local_numerator(without.as_ref());
        let numerator = match lmb {
            Some(basis) => {
                if !Arc::ptr_eq(&basis.owner, &self.owner) {
                    return Err(error::DiagramError::new_err(
                        "momentum basis belongs to a different diagram",
                    ));
                }
                basis.inner.route_expression(&numerator)
            }
            None if in_lmb => self
                .inner
                .loop_momentum_basis()
                .route_expression(&numerator),
            None => numerator,
        };
        TensorExpression::from_atom_interface(py, numerator, None)
    }

    /// Expand the local integrand through its UV degree of divergence.
    ///
    /// The selected region's loop momenta are scaled together. Propagators are
    /// expanded about the auxiliary mass ``uv_mass``, retaining every power
    /// through logarithmic divergence in ``dimension`` spacetime dimensions.
    /// The result uses edge momenta and the same tagged ``denom`` convention as
    /// :meth:`denominator_expression`; the diagram is unchanged.
    ///
    /// ``numerator`` optionally replaces the local numerator (in edge momenta),
    /// for example after contracting a projector. Overall factors, numerator
    /// prefactors and projectors remain separate unless supplied in it.
    /// Empty, tree and UV-convergent regions return zero. This performs one UV
    /// limit; it does not enumerate forests or integrate the counterterm.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> mass = S("mUV", is_scalar=True)
    /// >>> expansion = diagram.uv_expansion(mass)
    /// >>> region = diagram.filter(edge=lambda e: e.data.id in selected_edge_ids)
    /// >>> local_ct = region.uv_counterterm(mass)
    ///
    /// Parameters
    /// ----------
    /// uv_mass : Expression or int
    ///     Auxiliary mass used for the propagator expansion.
    /// dimension : int, optional
    ///     Positive spacetime dimension for UV power counting; defaults to four.
    /// numerator : Expression or TensorExpression or None, optional
    ///     Prepared numerator in edge momenta; None uses the selected local numerator.
    #[pyo3(signature = (uv_mass, *, dimension=4, numerator=None))]
    fn uv_expansion(
        &self,
        py: Python<'_>,
        uv_mass: ConvertibleToExpression,
        dimension: i32,
        numerator: Option<ConvertibleToExpression>,
    ) -> PyResult<Py<TensorExpression>> {
        let selected = self.selection();
        let numerator = numerator.map_or_else(
            || self.local_numerator(None),
            |value| value.to_expression().expr,
        );
        let expanded = self
            .inner
            .uv_expansion_of(
                &selected,
                &uv_mass.to_expression().expr,
                dimension,
                Some(&numerator),
            )
            .map_err(error::diagram)?;
        TensorExpression::from_atom_interface(py, expanded, None)
    }

    /// Return the additive local UV counterterm, the negative of ``uv_expansion``.
    ///
    /// Arguments and selection semantics are those of :meth:`uv_expansion`.
    /// Add this unintegrated expression to the selected integrand to subtract
    /// its simultaneous UV limit. Subdivergences require separate forest terms.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import S
    /// >>> mass = S("mUV", is_scalar=True)
    /// >>> counterterm = diagram.uv_counterterm(mass)
    ///
    /// Parameters
    /// ----------
    /// uv_mass : Expression or int
    ///     Auxiliary mass used for the propagator expansion.
    /// dimension : int, optional
    ///     Positive spacetime dimension for UV power counting; defaults to four.
    /// numerator : Expression or TensorExpression or None, optional
    ///     Prepared numerator in edge momenta; None uses the selected local numerator.
    #[pyo3(signature = (uv_mass, *, dimension=4, numerator=None))]
    fn uv_counterterm(
        &self,
        py: Python<'_>,
        uv_mass: ConvertibleToExpression,
        dimension: i32,
        numerator: Option<ConvertibleToExpression>,
    ) -> PyResult<Py<TensorExpression>> {
        let expanded = self.uv_expansion(py, uv_mass, dimension, numerator)?;
        let counterterm = -&expanded.borrow(py).as_super().expr;
        TensorExpression::from_atom_interface(py, counterterm, None)
    }

    /// Return the request-wide numerator multiplier as a Symbolica expression.
    ///
    /// Examples
    /// --------
    /// >>> prefactor = diagram.numerator_prefactor_expression()
    /// >>> weighted_numerator = prefactor * diagram.numerator_expression()
    /// >>> weighted_numerator
    ///
    fn numerator_prefactor_expression(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.numerator_prefactor().clone(),
        }
    }

    /// Return the external-state projector as a Symbolica expression.
    ///
    /// Examples
    /// --------
    /// >>> projector = diagram.projector_expression()
    /// >>> projected_numerator = projector * diagram.numerator_expression()
    /// >>> projected_numerator
    ///
    fn projector_expression(&self) -> PythonExpression {
        PythonExpression {
            expr: self.inner.projector().clone(),
        }
    }

    /// Reduce the numerator and projector using the diagram's internal edge momenta.
    ///
    /// Four-dimensional Lorentz slots are promoted to `D` before reduction.
    /// External momenta remain projector vectors, and the scalar numerator
    /// prefactor remains separate. Requires at least one internal edge.
    /// A supplied expression replaces the numerator and projector, for example
    /// after a UV expansion; it must use the diagram's edge momentum names.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import E
    /// >>> reduced = diagram.tensor_reduce(E("D"))
    ///
    /// Parameters
    /// ----------
    /// dimension : Expression or int
    ///     Lorentz dimension used for the reduction and four-dimensional input slots.
    /// expression : Expression or TensorExpression, optional
    ///     Complete prepared input to reduce instead of the numerator and projector.
    ///     Cannot be combined with an explicit projector.
    /// projector : Expression or None, optional
    ///     Explicit external tensor projector, required for partial regions unless expression is supplied.
    #[pyo3(signature = (dimension, *, expression=None, projector=None))]
    fn tensor_reduce(
        &self,
        py: Python<'_>,
        dimension: ConvertibleToExpression,
        expression: Option<ConvertibleToExpression>,
        projector: Option<ConvertibleToExpression>,
    ) -> PyResult<Py<TensorExpression>> {
        if expression.is_some() && projector.is_some() {
            return Err(error::DiagramError::new_err(
                "tensor reduction accepts either an expression or a projector, not both",
            ));
        }
        let selected = self.selection();
        let expression = expression.map(|expression| expression.to_expression().expr);
        let projector = match projector {
            Some(projector) => projector.to_expression().expr,
            None if expression.is_some() => Atom::one(),
            None if self.is_complete() => self.inner.projector().clone(),
            None => {
                return Err(error::DiagramError::new_err(
                    "partial-subgraph tensor reduction requires an explicit projector",
                ));
            }
        };
        let expression = expression.unwrap_or_else(|| self.local_numerator(None) * &projector);
        let reduction = self
            .inner
            .tensor_reduce_of(
                &selected,
                dimension.to_expression().expr,
                &projector,
                Some(&expression),
            )
            .map_err(error::tensor)?;
        TensorExpression::from_atom_interface(py, reduction.into_expression(), None)
    }

    /// Reduce the finalized numerator and projector with a tensor reducer.
    ///
    /// The result is a native Symbolica expression using ``spenso::dot`` and
    /// ``spenso::g``. For a fully projected vacuum numerator, all Lorentz
    /// indices disappear and only scalar invariants remain. Residual free
    /// projector indices are preserved as metric tensors; use
    /// :meth:`reduce_tensor_graphs` when scalar graph contributions are
    /// required. The scalar numerator prefactor remains separate and is not
    /// included in this expression.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import E
    /// >>> reducer = fk.TensorReducer.feynkit(E("4"))
    /// >>> scalar_numerator = vacuum_diagram.reduce_tensor_numerator(reducer)
    /// >>> scalar_numerator
    ///
    /// Parameters
    /// ----------
    /// reducer : TensorReducer
    ///     Tensor projector and integrated-momentum selection to apply.
    fn reduce_tensor_numerator(&self, reducer: &PyTensorReducer) -> PyResult<PythonExpression> {
        self.require_complete()?;
        self.inner
            .reduce_tensor_numerator(&reducer.inner)
            .map(|reduction| PythonExpression {
                expr: reduction.into_expression(),
            })
            .map_err(error::tensor)
    }

    /// Split the tensor numerator into scalar contributions on this topology.
    ///
    /// Each returned diagram has one compact reduction term as its numerator,
    /// including that term's exact projector coefficient. The graph topology,
    /// model, denominator data, and topology ID are preserved; deterministic
    /// ``.tensor[index]`` name suffixes distinguish the derived contributions.
    /// The external-state projector is consumed and reset to one, while the
    /// scalar numerator prefactor is retained. The projection must be fully
    /// contracted: any residual indexed Minkowski slot raises
    /// ``TensorReductionError``. Use
    /// :meth:`reduce_tensor_numerator` when residual free Lorentz indices are
    /// intentional.
    ///
    /// Examples
    /// --------
    /// >>> from symbolica import E
    /// >>> reducer = fk.TensorReducer.feynkit(E("4"))
    /// >>> scalar_graphs = vacuum_diagram.reduce_tensor_graphs(reducer)
    /// >>> for contribution in scalar_graphs:
    /// ...     print(contribution.name, contribution.numerator_expression())
    ///
    /// Parameters
    /// ----------
    /// reducer : TensorReducer
    ///     Tensor projector and integrated-momentum selection to apply.
    fn reduce_tensor_graphs(&self, reducer: &PyTensorReducer) -> PyResult<Vec<PyFeynmanDiagram>> {
        self.require_complete()?;
        self.inner
            .reduce_tensor_graphs(&reducer.inner)
            .map(|diagrams| diagrams.into_iter().map(Into::into).collect())
            .map_err(error::tensor)
    }

    /// Return physical final-state cuts selected during generation.
    #[getter]
    fn cuts(slf: Py<Self>, py: Python<'_>) -> Vec<PyDiagramCut> {
        let diagram = slf.borrow(py);
        let selected = diagram.selection();
        diagram
            .inner
            .cuts()
            .iter()
            .filter(|cut| {
                cut.left
                    .half_edges
                    .iter()
                    .chain(&cut.right.half_edges)
                    .chain(&cut.cut)
                    .all(|half| {
                        diagram
                            .inner
                            .half_edge_id(*half)
                            .is_some_and(|hedge| selected.includes(&hedge))
                    })
            })
            .map(|inner| PyDiagramCut {
                diagram: slf.clone_ref(py),
                inner: inner.clone(),
            })
            .collect()
    }

    /// Return topology threshold candidates separately from physical cuts.
    #[getter]
    fn topology_threshold_candidates(
        slf: Py<Self>,
        py: Python<'_>,
    ) -> Vec<PyDiagramThresholdCandidate> {
        let diagram = slf.borrow(py);
        let selected = diagram.selection();
        diagram
            .inner
            .topology_threshold_candidates()
            .iter()
            .filter(|cut| {
                cut.left
                    .iter()
                    .chain(&cut.right)
                    .chain(&cut.cut)
                    .all(|half| {
                        diagram
                            .inner
                            .half_edge_id(*half)
                            .is_some_and(|hedge| selected.includes(&hedge))
                    })
            })
            .map(|inner| PyDiagramThresholdCandidate {
                diagram: slf.clone_ref(py),
                inner: inner.clone(),
            })
            .collect()
    }

    /// Return a diagram using the specified loop-edge coordinates, in the supplied order.
    ///
    /// Examples
    /// --------
    /// >>> rerouted = diagram.with_loop_momentum_edges(diagram.loop_momentum_basis.loop_edges)
    ///
    /// Parameters
    /// ----------
    /// edges : list[int]
    ///     Diagram edge IDs specifying the requested routing coordinates.
    fn with_loop_momentum_edges(&self, edges: Vec<usize>) -> PyResult<Self> {
        self.require_complete()?;
        self.inner
            .as_ref()
            .clone()
            .with_loop_momentum_edges(
                &edges
                    .into_iter()
                    .map(feynkit_graph::EdgeId)
                    .collect::<Vec<_>>(),
            )
            .map(Into::into)
            .map_err(error::diagram)
    }

    /// Return a diagram whose momentum routing uses the specified spanning forest.
    ///
    /// Examples
    /// --------
    /// >>> rerouted = diagram.with_loop_momentum_tree_edges(diagram.loop_momentum_basis.tree_edges)
    ///
    /// Parameters
    /// ----------
    /// edges : list[int]
    ///     Diagram edge IDs specifying the requested routing coordinates.
    fn with_loop_momentum_tree_edges(&self, edges: Vec<usize>) -> PyResult<Self> {
        self.require_complete()?;
        self.inner
            .as_ref()
            .clone()
            .with_loop_momentum_tree_edges(
                &edges
                    .into_iter()
                    .map(feynkit_graph::EdgeId)
                    .collect::<Vec<_>>(),
            )
            .map(Into::into)
            .map_err(error::diagram)
    }

    /// Construct the canonical routing of the selected interaction region.
    ///
    /// Examples
    /// --------
    /// >>> basis = diagram.momentum_basis()
    /// >>> routed = basis.route_expression(diagram.numerator_expression())
    #[pyo3(signature = (*))]
    fn momentum_basis(&self) -> PyResult<PyLoopMomentumBasis> {
        self.inner
            .momentum_basis_of(&self.selection())
            .map(|basis| PyLoopMomentumBasis::from_diagram(basis, self))
            .map_err(error::diagram)
    }

    /// Reuse a parent basis's loop coordinates wherever the selected topology permits it.
    ///
    /// Examples
    /// --------
    /// >>> basis = diagram.compatible_momentum_basis(diagram.loop_momentum_basis)
    ///
    /// Parameters
    /// ----------
    /// parent : LoopMomentumBasis
    ///     Parent coordinates belonging to this same diagram instance.
    #[pyo3(signature = (parent))]
    fn compatible_momentum_basis(
        &self,
        parent: &PyLoopMomentumBasis,
    ) -> PyResult<PyLoopMomentumBasis> {
        if !Arc::ptr_eq(&parent.owner, &self.owner) {
            return Err(error::DiagramError::new_err(
                "parent momentum basis belongs to a different diagram",
            ));
        }
        self.inner
            .compatible_momentum_basis_of(&self.selection(), &parent.inner)
            .map(|basis| PyLoopMomentumBasis::from_diagram(basis, self))
            .map_err(error::diagram)
    }

    /// Route the selected region after contracting complete internal edges.
    ///
    /// Examples
    /// --------
    /// >>> contracted = diagram.filter(edge=lambda edge: edge.data.is_dummy)
    /// >>> basis = diagram.contracted_momentum_basis(contracted)
    ///
    /// Parameters
    /// ----------
    /// contracted : Subgraph or linnet.Subgraph
    ///     Complete internal edges to contract, selected from this diagram.
    #[pyo3(signature = (contracted))]
    fn contracted_momentum_basis(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="Subgraph | linnet.Subgraph", imports=("linnet")))]
        contracted: &Bound<'_, PyAny>,
    ) -> PyResult<PyLoopMomentumBasis> {
        let selected = self.selection();
        let contracted = self.selection_argument(py, contracted)?;
        let contracted = linnet::half_edge::subgraph::InternalSubGraph::try_new(
            contracted,
            self.inner.underlying(),
        )
        .ok_or_else(|| {
            error::DiagramError::new_err("contracted region must contain complete internal edges")
        })?;
        self.inner
            .contracted_momentum_basis_of(&selected, &contracted)
            .map(|basis| PyLoopMomentumBasis::from_diagram(basis, self))
            .map_err(error::diagram)
    }

    /// Return the loop-momentum routing selected during generation.
    #[getter]
    fn loop_momentum_basis(&self) -> PyResult<PyLoopMomentumBasis> {
        let basis = if self.selected_region.is_some() {
            self.inner
                .compatible_momentum_basis_of(&self.selection(), self.inner.loop_momentum_basis())
                .map_err(error::diagram)?
        } else {
            self.inner.loop_momentum_basis().clone()
        };
        Ok(PyLoopMomentumBasis::from_diagram(basis, self))
    }

    /// Return the number of independent loops in the diagram topology.
    ///
    /// Examples
    /// --------
    /// >>> basis = diagram.loop_momentum_bases(limit=1)[0]
    /// >>> len(basis.loop_edges) == diagram.loop_count
    /// True
    ///
    #[getter]
    fn loop_count(&self) -> PyResult<usize> {
        if self.selected_region.is_none() {
            return Ok(self.inner.loop_count());
        }
        self.inner
            .momentum_basis_of(&self.selection())
            .map(|basis| basis.loop_edges.len())
            .map_err(error::diagram)
    }

    /// Return the local superficial UV degree of divergence.
    ///
    /// Counts ``dimension * loops`` plus vertex momentum powers and internal
    /// propagator numerator powers minus two per internal propagator. Uses the
    /// stored local numerators; excludes external legs, projectors, and global
    /// prefactors. Vertex momenta scale together, before tensor cancellations.
    /// Zero is logarithmic, positive is power divergent, and negative is
    /// superficially convergent. Subdivergences are not tested.
    ///
    /// Examples
    /// --------
    /// >>> degree = diagram.superficial_degree_of_divergence()
    /// >>> degree_in_six_dimensions = diagram.superficial_degree_of_divergence(dimension=6)
    ///
    /// Parameters
    /// ----------
    /// dimension : int, optional
    ///     Spacetime dimension for each loop integration measure; defaults to four.
    #[pyo3(signature = (*, dimension=4))]
    fn superficial_degree_of_divergence(&self, dimension: i32) -> PyResult<i32> {
        self.inner
            .superficial_degree_of_divergence_of(&self.selection(), dimension)
            .map_err(error::diagram)
    }

    /// Return the diagram's interaction vertices with stable integer identifiers.
    #[getter]
    pub(crate) fn vertices(&self) -> Vec<PyDiagramVertex> {
        let model = self.inner.model_arc();
        self.inner
            .vertices()
            .filter(|(id, _)| self.includes_vertex(id.0))
            .map(|(id, inner)| PyDiagramVertex {
                id: id.0,
                inner: inner.clone(),
                model: Arc::clone(&model),
            })
            .collect()
    }

    /// Return every particle line, including dangling and sewn external carriers.
    #[getter]
    pub(crate) fn edges(&self) -> Vec<PyDiagramEdge> {
        let model = self.inner.model_arc();
        self.inner
            .edges()
            .filter(|(id, _, _)| self.includes_edge(*id))
            .map(|(id, endpoints, inner)| PyDiagramEdge {
                id: id.0,
                source: endpoints.source.map(|vertex| vertex.0),
                target: endpoints.target.map(|vertex| vertex.0),
                inner: inner.clone(),
                model: Arc::clone(&model),
            })
            .collect()
    }

    /// Return native Linnet half-edge views; ``data`` records their native diagram IDs.
    #[getter]
    #[gen_stub(override_return_type(type_repr="list[linnet.HalfEdge]", imports=("linnet")))]
    fn half_edges(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        Ok(self
            .linnet_selection(py)?
            .bind(py)
            .call_method0("to_half_edges")?
            .unbind())
    }

    /// Return propagators excluding external momentum carriers and dummy edges.
    ///
    /// Examples
    /// --------
    /// >>> propagators = diagram.internal_edges
    /// >>> on_shell = {edge.id: energy[edge.id] for edge in propagators}
    #[getter]
    fn internal_edges(&self) -> Vec<PyDiagramEdge> {
        let complete: BTreeSet<_> = self
            .inner
            .underlying()
            .iter_edges_of(&self.selection())
            .filter_map(|(pair, id, _)| pair.is_paired().then_some(id.0))
            .collect();
        self.edges()
            .into_iter()
            .filter(|edge| {
                !edge.is_external() && !edge.inner.is_dummy && complete.contains(&edge.id)
            })
            .collect()
    }

    /// Return incoming or outgoing external-state momentum carriers.
    ///
    /// Examples
    /// --------
    /// >>> external_particles = [edge.particle_name for edge in diagram.external_edges]
    #[getter]
    fn external_edges(&self) -> Vec<PyDiagramEdge> {
        self.edges()
            .into_iter()
            .filter(PyDiagramEdge::is_external)
            .collect()
    }

    /// Validate particle and interaction references against a physics model.
    ///
    /// Examples
    /// --------
    /// Validation raises ``DiagramError`` for invalid particle or interaction
    /// references:
    ///
    /// >>> diagram.validate()
    fn validate(&self) -> PyResult<()> {
        self.inner.validate().map_err(error::diagram)
    }

    /// Build the diagram's Cross-Free Family representation.
    ///
    /// Edge constraints use the stable integer IDs exposed by
    /// ``diagram.edges``. ``False`` fixes an edge in its stored direction and
    /// ``True`` reverses it.
    ///
    /// Examples
    /// --------
    /// Construct and display the causal denominators of a one-loop diagram:
    ///
    /// >>> cff = diagram.build_cff(max_orientations=10_000)
    /// >>> cff.to_expression()  # native Symbolica display in a notebook
    ///
    /// Parameters
    /// ----------
    /// max_orientations : int or None, optional
    ///     Maximum number of candidate orientations to inspect.
    /// fixed_orientations : mapping[int, bool] or None, optional
    ///     Edge IDs mapped to stored (false) or reversed (true) directions.
    /// contracted_edges : iterable[int], optional
    ///     Edge IDs to contract before constructing denominator surfaces.
    /// initial_state_edges : iterable[int], optional
    ///     Edge IDs to classify as incoming external lines.
    #[pyo3(signature = (*, max_orientations=None, fixed_orientations=None, contracted_edges=None, initial_state_edges=None))]
    fn build_cff(
        &self,
        py: Python<'_>,
        max_orientations: Option<usize>,
        fixed_orientations: Option<BTreeMap<usize, bool>>,
        contracted_edges: Option<Vec<usize>>,
        initial_state_edges: Option<Vec<usize>>,
    ) -> PyResult<PyCffResult> {
        build_cff_for_diagram(
            py,
            self,
            max_orientations,
            fixed_orientations,
            contracted_edges,
            initial_state_edges,
        )
    }

    /// Serialize the complete diagram to JSON text.
    ///
    /// Examples
    /// --------
    /// >>> encoded = diagram.to_json()
    /// >>> restored = FeynmanDiagram.from_json(model, encoded)
    /// >>> restored.validate()
    ///
    fn to_json(&self) -> PyResult<String> {
        self.require_complete()?;
        self.inner.to_json().map_err(error::diagram)
    }

    /// Serialize the diagram topology and annotations to Graphviz DOT text.
    ///
    /// Examples
    /// --------
    /// >>> dot = diagram.to_dot()
    /// >>> restored = FeynmanDiagram.from_dot(model, dot)
    /// >>> restored.validate()
    ///
    fn to_dot(&self) -> PyResult<String> {
        self.require_complete()?;
        self.inner.to_dot().map_err(error::diagram)
    }

    /// Emit the exact Typst source used by ``render`` without compiling it.
    ///
    /// Uses the same ``config``, ``momenta``, ``lmb`` and ``highlight`` settings
    /// as ``render``. The shared Linnest/Kurvst and physics assets must be available
    /// beneath the Typst project root when compiling this source separately.
    ///
    /// Examples
    /// --------
    /// >>> from pathlib import Path
    /// >>> Path("diagram.typ").write_text(diagram.to_linnest(momenta=True))
    ///
    /// Parameters
    /// ----------
    /// config : linnet.RenderConfig or None, optional
    ///     Layout, drawing, style and physics settings, as in ``render``.
    /// momenta : bool, optional
    ///     Draw momentum arrows and labels in the stored basis.
    /// lmb : LoopMomentumBasis or None, optional
    ///     Routing from this diagram; also enables momentum display.
    /// highlight : Subgraph or linnet.Subgraph or None, optional
    ///     Region to highlight in the complete diagram.
    #[pyo3(signature = (*, config=None, momenta=false, lmb=None, highlight=None))]
    fn to_linnest(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="linnet.RenderConfig | None", imports=("linnet")))]
        config: Option<&Bound<'_, PyAny>>,
        momenta: bool,
        lmb: Option<&PyLoopMomentumBasis>,
        #[gen_stub(override_type(type_repr="Subgraph | linnet.Subgraph | None", imports=("linnet")))]
        highlight: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<String> {
        self.prepare_render(py, config, momenta, lmb, highlight)?
            .getattr("typst_source")?
            .extract()
    }

    /// Render an interactive, transparent SVG using the shared physics renderer.
    ///
    /// Examples
    /// --------
    /// >>> import linnet as ln
    /// >>> svg = diagram.render(momenta=True, config=ln.RenderConfig(
    /// ...     layouts=ln.LayoutOptions(external_label_length_scale=0.7),
    /// ...     template_options={"show-particle": False},
    /// ... ))
    /// >>> svg = diagram.render(lmb=next(iter(diagram.loop_momentum_bases())))
    ///
    /// Parameters
    /// ----------
    /// config : linnet.RenderConfig or None, optional
    ///     Typed ``layouts``, ``drawing`` and ``style`` groups. Physics controls
    ///     use ``template_options`` with the same names as ``just draw --input``:
    ///     ``show-particle``, ``show-edge-index``, ``show-node-index``, ``debug``,
    ///     ``momentum-arrows`` and the ``momentum-arrow-*``/``momentum-label-*`` options.
    ///     Cross sections open their initial-state connections by default; set
    ///     ``split-initial-state`` to ``False`` to draw the sewn graph.
    /// momenta : bool, optional
    ///     Show momentum arrows and labels routed in the diagram's stored basis.
    ///     Explicit physics settings in ``config`` override these display defaults.
    /// lmb : LoopMomentumBasis or None, optional
    ///     Explicit routing from this diagram; also enables momentum display.
    ///     Rendering never changes the diagram's stored loop-momentum basis.
    /// highlight : Subgraph or linnet.Subgraph or None, optional
    ///     Highlight a region while preserving the full diagram as muted context.
    ///     A Subgraph highlights its own region by default.
    ///
    #[pyo3(signature = (*, config=None, momenta=false, lmb=None, highlight=None))]
    fn render(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="linnet.RenderConfig | None", imports=("linnet")))]
        config: Option<&Bound<'_, PyAny>>,
        momenta: bool,
        lmb: Option<&PyLoopMomentumBasis>,
        #[gen_stub(override_type(type_repr="Subgraph | linnet.Subgraph | None", imports=("linnet")))]
        highlight: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<String> {
        let prepared = self.prepare_render(py, config, momenta, lmb, highlight)?;
        render_diagram_svg(py, &prepared)
    }

    /// Render an HTML figure with the same options and hover information as ``render``.
    ///
    /// Examples
    /// --------
    /// >>> import marimo as mo
    /// >>> mo.iframe(diagram.to_html(momenta=True))
    ///
    /// Parameters
    /// ----------
    /// config : linnet.RenderConfig or None, optional
    ///     Layout, drawing, style and physics settings, as in ``render``.
    /// momenta : bool, optional
    ///     Draw momentum arrows and labels in the stored basis.
    /// lmb : LoopMomentumBasis or None, optional
    ///     Routing from this diagram; also enables momentum display.
    /// highlight : Subgraph or linnet.Subgraph or None, optional
    ///     Region to highlight in the complete diagram.
    #[pyo3(signature = (*, config=None, momenta=false, lmb=None, highlight=None))]
    fn to_html(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="linnet.RenderConfig | None", imports=("linnet")))]
        config: Option<&Bound<'_, PyAny>>,
        momenta: bool,
        lmb: Option<&PyLoopMomentumBasis>,
        #[gen_stub(override_type(type_repr="Subgraph | linnet.Subgraph | None", imports=("linnet")))]
        highlight: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<String> {
        let svg = self.render(py, config, momenta, lmb, highlight)?;
        Ok(render_diagram_html(&self.inner, &svg))
    }

    /// Render the diagram as HTML in Marimo, Jupyter, and IPython.
    ///
    /// Examples
    /// --------
    /// Leave `diagram` as the final expression in a notebook cell to render it.
    pub(crate) fn _repr_html_(&self, py: Python<'_>) -> PyResult<String> {
        self.to_html(py, None, false, None, None)
    }

    /// Return the raw SVG representation used by rich notebook frontends.
    ///
    /// Examples
    /// --------
    /// >>> from IPython.display import SVG
    /// >>> SVG(diagram._repr_svg_())
    fn _repr_svg_(&self, py: Python<'_>) -> PyResult<String> {
        self.render(py, None, false, None, None)
    }

    /// Write a concise summary to an IPython pretty printer.
    ///
    /// Examples
    /// --------
    /// IPython invokes this method when only a text representation is supported.
    ///
    /// Parameters
    /// ----------
    /// pretty : Any
    ///     The IPython pretty-printer object.
    /// cycle : bool
    ///     Whether this object is part of a recursive formatting cycle.
    fn _repr_pretty_(&self, pretty: &Bound<'_, PyAny>, cycle: bool) -> PyResult<()> {
        let text = if cycle {
            "...".to_owned()
        } else {
            self.__repr__()?
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
    /// >>> basis = diagram.loop_momentum_bases(limit=1)[0]
    /// >>> routing = {
    /// ...     edge_id: signature.format_momentum()
    /// ...     for edge_id, signature in basis.edge_signatures.items()
    /// ... }
    ///
    /// Parameters
    /// ----------
    /// limit : int or None
    ///     Maximum number of bases to return. Pass ``None`` to enumerate every
    ///     valid basis.
    #[pyo3(signature = (limit=None))]
    fn loop_momentum_bases(
        &self,
        py: Python<'_>,
        limit: Option<usize>,
    ) -> PyResult<Vec<PyLoopMomentumBasis>> {
        let selected = self.selection();
        let diagram = self.inner.clone();
        py.detach(move || diagram.loop_momentum_bases_of(&selected, limit.unwrap_or(usize::MAX)))
            .map(|bases| {
                bases
                    .into_iter()
                    .map(|basis| PyLoopMomentumBasis::from_diagram(basis, self))
                    .collect()
            })
            .map_err(error::diagram)
    }

    /// Return a concise description of the diagram and its loop count.
    ///
    /// Examples
    /// --------
    /// >>> print(diagram)
    ///
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "{}(name='{}', loops={}, vertices={}, edges={})",
            if self.selected_region.is_some() {
                "Subgraph"
            } else {
                "FeynmanDiagram"
            },
            self.inner.name(),
            self.loop_count()?,
            self.vertices().len(),
            self.edges().len()
        ))
    }
}

/// A physics region retaining an Arc to its immutable original Feynman diagram.
///
/// Construct with ``diagram.subgraph(...)`` or ``diagram.filter(...)``. Inherited
/// physics operations use this selection, with original edge IDs and momentum
/// coordinates. Display keeps the full original graph as muted context.
/// ``excise()`` creates independent topology and routing when those are needed.
///
/// Examples
/// --------
/// >>> region = diagram.filter(edge=lambda edge: not edge.data.is_external)
/// >>> numerator = region.numerator_expression()
/// >>> independent = region.excise()
#[cfg_attr(feature = "python_stubgen", gen_stub_pyclass)]
#[pyclass(name = "Subgraph", module = "symbolica.community.feynkit", extends=PyFeynmanDiagram, frozen)]
pub struct PySubgraph {}

#[cfg_attr(feature = "python_stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PySubgraph {
    /// Return a full view of the immutable original diagram, sharing its Arc.
    #[getter]
    fn original(slf: PyRef<'_, Self>) -> PyFeynmanDiagram {
        slf.into_super().whole()
    }

    /// Materialize an independent diagram with dangling legs at the region boundary.
    ///
    /// Complete interaction crowns are retained. Only selected internal pairs stay
    /// paired; remaining incident halves become boundary legs. Proper regions get
    /// independent momentum coordinates and unit global weights and projector.
    /// Excising the complete original preserves its factors, cuts and routing.
    ///
    /// Examples
    /// --------
    /// >>> region = diagram.subgraph(nodes=[0])
    /// >>> independent = region.excise()
    /// >>> independent.validate()
    fn excise(slf: PyRef<'_, Self>) -> PyResult<PyFeynmanDiagram> {
        let diagram = slf.into_super();
        diagram
            .inner
            .excise(
                &diagram.selection(),
                &diagram
                    .isolated_nodes()
                    .into_iter()
                    .map(feynkit_graph::VertexId)
                    .collect::<Vec<_>>(),
            )
            .map(Into::into)
            .map_err(error::diagram)
    }

    /// Return canonical Linnet half-edge IDs of this region.
    ///
    /// Examples
    /// --------
    /// >>> selected_half_edges = region.half_edge_indices()
    /// >>> canonical = region.to_linnet().subgraph(half_edges=selected_half_edges)
    fn half_edge_indices(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Vec<usize>> {
        slf.into_super()
            .linnet_selection(py)?
            .bind(py)
            .call_method0("half_edge_indices")?
            .extract()
    }

    /// Return explicitly selected vertices with no incident half-edges.
    ///
    /// Examples
    /// --------
    /// >>> isolated_vertices = region.isolated_node_indices()
    fn isolated_node_indices(slf: PyRef<'_, Self>) -> Vec<usize> {
        slf.into_super().isolated_nodes().into_iter().collect()
    }

    /// Return the number of selected half-edges, excluding isolated vertices.
    #[getter]
    fn n_half_edges(slf: PyRef<'_, Self>) -> usize {
        slf.into_super().selection().n_included()
    }

    /// Count selected half-edges and explicitly selected isolated vertices.
    ///
    /// Examples
    /// --------
    /// >>> size = len(region)
    fn __len__(slf: PyRef<'_, Self>) -> usize {
        let diagram = slf.into_super();
        diagram.selection().n_included() + diagram.isolated_nodes().len()
    }

    /// Test whether any half-edge or isolated vertex is selected.
    ///
    /// Examples
    /// --------
    /// >>> nonempty = bool(region)
    fn __bool__(slf: PyRef<'_, Self>) -> bool {
        Self::__len__(slf) != 0
    }

    /// Intersect two regions belonging to the same original diagram.
    ///
    /// Examples
    /// --------
    /// >>> common = left & right
    ///
    /// Parameters
    /// ----------
    /// other : FeynmanDiagram
    ///     Region or full diagram sharing this region's immutable original.
    fn __and__(
        slf: PyRef<'_, Self>,
        py: Python<'_>,
        other: &PyFeynmanDiagram,
    ) -> PyResult<Py<PySubgraph>> {
        let diagram = slf.into_super();
        if !Arc::ptr_eq(&diagram.inner, &other.inner) {
            return Err(PyValueError::new_err(
                "subgraphs belong to different original diagrams",
            ));
        }
        diagram.view(
            py,
            DiagramSelection {
                hedges: diagram.selection().intersection(&other.selection()),
                isolated: diagram
                    .isolated_nodes()
                    .intersection(&other.isolated_nodes())
                    .copied()
                    .collect(),
            },
        )
    }

    /// Unite two regions belonging to the same original diagram.
    ///
    /// Examples
    /// --------
    /// >>> combined = left | right
    ///
    /// Parameters
    /// ----------
    /// other : FeynmanDiagram
    ///     Region or full diagram sharing this region's immutable original.
    fn __or__(
        slf: PyRef<'_, Self>,
        py: Python<'_>,
        other: &PyFeynmanDiagram,
    ) -> PyResult<Py<PySubgraph>> {
        let diagram = slf.into_super();
        if !Arc::ptr_eq(&diagram.inner, &other.inner) {
            return Err(PyValueError::new_err(
                "subgraphs belong to different original diagrams",
            ));
        }
        diagram.view(
            py,
            DiagramSelection {
                hedges: diagram.selection().union(&other.selection()),
                isolated: diagram
                    .isolated_nodes()
                    .union(&other.isolated_nodes())
                    .copied()
                    .collect(),
            },
        )
    }

    /// Select elements present in exactly one of two regions.
    ///
    /// Examples
    /// --------
    /// >>> exclusive = left ^ right
    ///
    /// Parameters
    /// ----------
    /// other : FeynmanDiagram
    ///     Region or full diagram sharing this region's immutable original.
    fn __xor__(
        slf: PyRef<'_, Self>,
        py: Python<'_>,
        other: &PyFeynmanDiagram,
    ) -> PyResult<Py<PySubgraph>> {
        let diagram = slf.into_super();
        if !Arc::ptr_eq(&diagram.inner, &other.inner) {
            return Err(PyValueError::new_err(
                "subgraphs belong to different original diagrams",
            ));
        }
        diagram.view(
            py,
            DiagramSelection {
                hedges: diagram.selection().sym_diff(&other.selection()),
                isolated: diagram
                    .isolated_nodes()
                    .symmetric_difference(&other.isolated_nodes())
                    .copied()
                    .collect(),
            },
        )
    }

    /// Remove the other region's elements from this selection.
    ///
    /// Examples
    /// --------
    /// >>> remaining = left - right
    ///
    /// Parameters
    /// ----------
    /// other : FeynmanDiagram
    ///     Region or full diagram sharing this region's immutable original.
    fn __sub__(
        slf: PyRef<'_, Self>,
        py: Python<'_>,
        other: &PyFeynmanDiagram,
    ) -> PyResult<Py<PySubgraph>> {
        let diagram = slf.into_super();
        if !Arc::ptr_eq(&diagram.inner, &other.inner) {
            return Err(PyValueError::new_err(
                "subgraphs belong to different original diagrams",
            ));
        }
        diagram.view(
            py,
            DiagramSelection {
                hedges: diagram.selection().subtract(&other.selection()),
                isolated: diagram
                    .isolated_nodes()
                    .difference(&other.isolated_nodes())
                    .copied()
                    .collect(),
            },
        )
    }

    /// Complement the region relative to its immutable original diagram.
    ///
    /// Examples
    /// --------
    /// >>> outside = ~region
    fn __invert__(slf: PyRef<'_, Self>, py: Python<'_>) -> PyResult<Py<PySubgraph>> {
        let diagram = slf.into_super();
        diagram.view(
            py,
            DiagramSelection {
                hedges: diagram
                    .inner
                    .underlying()
                    .full_filter()
                    .subtract(&diagram.selection()),
                isolated: diagram
                    .whole()
                    .isolated_nodes()
                    .difference(&diagram.isolated_nodes())
                    .copied()
                    .collect(),
            },
        )
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyDiagramVertex>()?;
    module.add_class::<PyDiagramEdge>()?;
    module.add_class::<PyMomentumSignature>()?;
    module.add_class::<PyLoopMomentumBasis>()?;
    module.add_class::<PyFeynmanDiagram>()?;
    module.add_class::<PySubgraph>()?;
    module.add_class::<PyDiagramCut>()?;
    module.add_class::<PyDiagramCutSide>()?;
    module.add_class::<PyDiagramThresholdCandidate>()?;
    Ok(())
}

impl PyFeynmanDiagram {
    fn prepare_render<'py>(
        &self,
        py: Python<'py>,
        config: Option<&Bound<'_, PyAny>>,
        momenta: bool,
        lmb: Option<&PyLoopMomentumBasis>,
        highlight: Option<&Bound<'_, PyAny>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        if let Some(basis) = lmb
            && !Arc::ptr_eq(&basis.owner, &self.owner)
        {
            return Err(error::DiagramError::new_err(
                "momentum basis belongs to a different diagram",
            ));
        }
        let selected = highlight
            .map(|value| self.region_argument(py, value))
            .transpose()?
            .or_else(|| self.selected_region.clone());
        let isolated = selected
            .as_ref()
            .map_or_else(BTreeSet::new, |region| region.isolated.clone());
        let selected = selected.as_ref().map(|region| &region.hedges);
        let source = self
            .inner
            .to_linnest(
                selected,
                &isolated,
                lmb.map(|basis| &basis.inner),
                "_linnet_config",
            )
            .map_err(error::diagram)?;
        let linnet = py.import("linnet").map_err(|error| {
            if error.is_instance_of::<pyo3::exceptions::PyImportError>(py) {
                pyo3::exceptions::PyImportError::new_err(format!(
                    "diagram rendering requires linnet and typst-py: {error}"
                ))
            } else {
                error
            }
        })?;
        let options = PyDict::new(py);
        if momenta || lmb.is_some() {
            options.set_item("momentum-arrows", true)?;
            options.set_item("show-momentum", true)?;
        }
        let kwargs = PyDict::new(py);
        kwargs.set_item("template_options", options)?;
        let mut effective = linnet.getattr("RenderConfig")?.call((), Some(&kwargs))?;
        if let Some(config) = config {
            effective = effective.call_method1("overlay", (config,))?;
        }
        let sources = PyDict::new(py);
        sources.set_item("main.typ", PyBytes::new(py, source.as_bytes()))?;
        for (path, source) in [
            (
                "assets/embedded/drawing/templates/layout-core.typ",
                include_bytes!("../../../assets/embedded/drawing/templates/layout-core.typ")
                    .as_slice(),
            ),
            (
                "assets/embedded/drawing/templates/physics-edge-style.typ",
                include_bytes!("../../../assets/embedded/drawing/templates/physics-edge-style.typ")
                    .as_slice(),
            ),
            (
                "assets/embedded/drawing/templates/impl/physics-edge-style.typ",
                include_bytes!(
                    "../../../assets/embedded/drawing/templates/impl/physics-edge-style.typ"
                )
                .as_slice(),
            ),
        ] {
            sources.set_item(path, PyBytes::new(py, source))?;
        }
        let kwargs = PyDict::new(py);
        kwargs.set_item("config", effective)?;
        linnet
            .getattr("PreparedRender")?
            .call_method("from_sources", (sources,), Some(&kwargs))
    }
}

#[cfg(test)]
mod tests {
    use super::parse_symbolic_annotation;

    #[test]
    fn canonical_denominator_has_a_scalar_interface() {
        use symbolica::atom::{Atom, FunctionBuilder};

        let denominator = feynkit_graph::symbols::denominator();
        assert!(denominator.is_scalar());
        let momentum = FunctionBuilder::new(feynkit_graph::momentum_symbol())
            .add_arg(0)
            .finish();
        let atom = FunctionBuilder::new(denominator)
            .add_arg(0)
            .add_arg(momentum)
            .add_arg(0)
            .add_arg(Atom::one())
            .finish();
        let symbolica::atom::AtomView::Fun(function) = atom.as_view() else {
            unreachable!()
        };
        assert!(function.get_symbol().is_scalar());
        pyo3::Python::initialize();
        pyo3::Python::attach(|py| {
            spynso3::expression::TensorExpression::from_atom_interface(py, atom, None).unwrap();
        });
    }

    #[test]
    fn symbolic_annotations_are_parsed_fallibly() {
        assert!(parse_symbolic_annotation("g^2*(x+y)").is_ok());
        assert!(parse_symbolic_annotation("g+").is_err());
    }
}
