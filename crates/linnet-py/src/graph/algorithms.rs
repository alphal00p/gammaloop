//! Graphica algorithms over a lossless, colored incidence graph.
use super::*;
use std::sync::{Arc, Mutex};
use symbolica::{
    api::python::{ConvertibleToExpression, PythonExpression},
    atom::Atom,
    graph::{GenerationSettings, Graph, HalfEdge},
};

/// A labeled port used in graph generation signatures.
/// Direction is outgoing (True), incoming (False), or undirected (None).
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pyclass)]
#[pyclass(
    module = "symbolica.community.graph",
    name = "EdgeSignature",
    from_py_object,
    frozen
)]
#[derive(Clone)]
pub struct PyEdgeSignature {
    edge: HalfEdge<Atom>,
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[pymethods]
impl PyEdgeSignature {
    /// A labeled outgoing (True), incoming (False), or undirected (None) port.
    #[new]
    #[pyo3(signature=(data, direction=None))]
    fn new(data: ConvertibleToExpression, direction: Option<bool>) -> Self {
        let data = data.to_expression().expr;
        Self {
            edge: match direction {
                Some(true) => HalfEdge::outgoing(data),
                Some(false) => HalfEdge::incoming(data),
                None => HalfEdge::undirected(data),
            },
        }
    }
    fn flip(&self) -> Self {
        Self {
            edge: self.edge.flip(),
        }
    }
    #[getter]
    fn data(&self) -> PythonExpression {
        self.edge.data.clone().into()
    }
    #[getter]
    fn direction(&self) -> Option<bool> {
        self.edge.direction
    }
}

type Incidence = Graph<(u8, Atom), ()>;
type CanonicalGraph = (PyGraph, Vec<usize>, Py<PyAny>, Vec<usize>);
struct Encoding {
    graph: Incidence,
    nodes: Vec<usize>,
    edges: Vec<usize>,
    hedges: Vec<usize>,
}
impl Encoding {
    fn key(py: Python<'_>, data: &Py<PyAny>, callback: Option<&Py<PyAny>>) -> PyResult<Atom> {
        let value = match callback {
            Some(f) => f.call1(py, (data,))?,
            None => data.clone_ref(py),
        };
        if value.is_none(py) {
            return Ok(Atom::num(0));
        }
        value.extract::<ConvertibleToExpression>(py).map(|x|x.to_expression().expr)
            .map_err(|_|PyTypeError::new_err("graph labels must be expression-compatible; supply an element key callback for custom payloads"))
    }
    fn new(
        py: Python<'_>,
        state: &GraphState,
        node_key: Option<&Py<PyAny>>,
        edge_key: Option<&Py<PyAny>>,
        half_edge_key: Option<&Py<PyAny>>,
    ) -> PyResult<Self> {
        let mut graph = Graph::new();
        let mut nodes = vec![0; state.graph.n_nodes()];
        let mut edges = vec![0; state.graph.n_edges()];
        let mut hedges = vec![0; state.graph.n_hedges()];
        for (index, _, node) in state.graph.iter_nodes() {
            nodes[index.0] = graph.add_node((0, Self::key(py, &node.data, node_key)?));
        }
        for (pair, index, edge) in state.graph.iter_edges() {
            let center = graph.add_node((1, Self::key(py, &edge.data.data, edge_key)?));
            edges[index.0] = center;
            let ends = match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    vec![(source, Flow::Source), (sink, Flow::Sink)]
                }
                HedgePair::Unpaired { hedge, flow } => vec![(hedge, flow)],
            };
            for (hedge, flow) in ends {
                // Undirected source/sink storage order is not a graph invariant.
                let role = match edge.orientation {
                    Orientation::Undirected => 2,
                    Orientation::Default => {
                        if flow == Flow::Source {
                            3
                        } else {
                            4
                        }
                    }
                    Orientation::Reversed => {
                        if flow == Flow::Source {
                            4
                        } else {
                            3
                        }
                    }
                };
                let port = graph.add_node((
                    role,
                    Self::key(py, &state.graph[hedge].data, half_edge_key)?,
                ));
                hedges[hedge.0] = port;
                graph
                    .add_edge(center, port, false, ())
                    .map_err(PyValueError::new_err)?;
                graph
                    .add_edge(port, nodes[state.graph.node_id(hedge).0], false, ())
                    .map_err(PyValueError::new_err)?;
            }
        }
        Ok(Self {
            graph,
            nodes,
            edges,
            hedges,
        })
    }
    fn order(indices: &[usize], map: &[usize]) -> Vec<usize> {
        let mut order = (0..indices.len()).collect::<Vec<_>>();
        order.sort_by_key(|&i| map[indices[i]]);
        order
    }
}
fn python_integer(py: Python<'_>, value: impl std::fmt::Display) -> PyResult<Py<PyAny>> {
    Ok(py
        .import("builtins")?
        .getattr("int")?
        .call1((value.to_string(),))?
        .unbind())
}
impl PyGraph {
    /// Construct the common Python graph directly from a native Graphica graph.
    pub fn from_symbolica(py: Python<'_>, graph: &Graph<Atom, Atom>) -> PyResult<Self> {
        let mut builder = HedgeGraphBuilder::new();
        let nodes = graph
            .nodes()
            .iter()
            .map(|node| {
                Ok(builder.add_node(NodeRecord {
                    name: None,
                    data: Py::new(py, PythonExpression::from(node.data.clone()))?.into_any(),
                    drawing: drawing_dict(py, NODE_FIELDS, None)?,
                }))
            })
            .collect::<PyResult<Vec<_>>>()?;
        for edge in graph.edges() {
            let endpoint = |index| -> PyResult<_> {
                Ok(HedgeData {
                    is_in_subgraph: false,
                    node: nodes[index],
                    data: HalfEdgeRecord {
                        data: py.None(),
                        drawing: drawing_dict(py, HEDGE_FIELDS, None)?,
                    },
                })
            };
            builder.add_edge(
                endpoint(edge.vertices.0)?,
                endpoint(edge.vertices.1)?,
                EdgeRecord {
                    name: None,
                    data: Py::new(py, PythonExpression::from(edge.data.clone()))?.into_any(),
                    drawing: drawing_dict(py, EDGE_FIELDS, None)?,
                },
                if edge.directed {
                    Orientation::Default
                } else {
                    Orientation::Undirected
                },
            );
        }
        Ok(Self::from_state(GraphState {
            graph: PyHedgeGraph::build(builder, PyNodeStore::Vec),
            name: None,
            global_data: PyGlobalData::default(),
            codec: None,
            render_config: crate::typst::default_render_config(py)?,
            revision: 0,
        }))
    }
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyGraph {
    fn __copy__(&self, py: Python<'_>) -> PyResult<Self> {
        Ok(Self::from_state(
            self.state
                .borrow()
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
                .copy(py)?,
        ))
    }
    /// Canonical graph, old-to-new vertex map, automorphism-group size, and vertex orbit.
    /// Keys receive element payloads; names and drawing metadata do not affect equivalence.
    #[pyo3(signature=(*, node_key=None, edge_key=None, half_edge_key=None))]
    #[gen_stub(override_return_type(type_repr = "tuple[Graph, list[int], int, list[int]]"))]
    fn canonize(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] node_key: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] edge_key: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] half_edge_key: Option<Py<PyAny>>,
    ) -> PyResult<CanonicalGraph> {
        // Release the original state before invoking user key callbacks.
        let state = self
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy(py)?;
        let encoded = Encoding::new(
            py,
            &state,
            node_key.as_ref(),
            edge_key.as_ref(),
            half_edge_key.as_ref(),
        )?;
        let canonical = py.detach(|| encoded.graph.canonize());
        let node_order = Encoding::order(&encoded.nodes, &canonical.vertex_map);
        let mut vertex_map = vec![0; node_order.len()];
        for (new, &old) in node_order.iter().enumerate() {
            vertex_map[old] = new;
        }
        let canonical_nodes = node_order
            .iter()
            .map(|&i| canonical.vertex_map[encoded.nodes[i]])
            .collect::<Vec<_>>();
        let orbit = canonical_nodes
            .iter()
            .map(|&v| {
                canonical_nodes
                    .iter()
                    .position(|&n| n == canonical.orbit[v])
                    .expect("vertex orbit preserves the node color")
            })
            .collect();
        let result = Self::from_state(state);
        {
            let mut borrow = result.state.borrow_mut();
            let graph = &mut borrow.as_mut().unwrap().graph;
            let flips = graph
                .iter_edges()
                .filter_map(|(pair, index, edge)| {
                    let flip = match pair {
                        HedgePair::Paired { source, sink } => match edge.orientation {
                            Orientation::Reversed => Some(source),
                            Orientation::Undirected
                                if canonical.vertex_map[encoded.hedges[source.0]]
                                    > canonical.vertex_map[encoded.hedges[sink.0]] =>
                            {
                                Some(source)
                            }
                            _ => None,
                        },
                        HedgePair::Unpaired { hedge, .. }
                            if edge.orientation == Orientation::Reversed =>
                        {
                            Some(hedge)
                        }
                        _ => None,
                    };
                    flip.map(|h| (index, h, edge.orientation))
                })
                .collect::<Vec<_>>();
            for (index, hedge, orientation) in flips {
                graph.set_flow(hedge, -graph.flow(hedge));
                graph.set_orientation(
                    index,
                    if orientation == Orientation::Reversed {
                        Orientation::Default
                    } else {
                        orientation
                    },
                );
            }
        }
        result.reorder_nodes(node_order)?;
        result.reorder_edges(Encoding::order(&encoded.edges, &canonical.vertex_map))?;
        Ok((
            result,
            vertex_map,
            python_integer(py, canonical.automorphism_group_size)?,
            orbit,
        ))
    }
    #[pyo3(signature=(other, *, node_key=None, edge_key=None, half_edge_key=None))]
    fn is_isomorphic(
        &self,
        py: Python<'_>,
        other: &Self,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] node_key: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] edge_key: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] half_edge_key: Option<Py<PyAny>>,
    ) -> PyResult<bool> {
        let copy = |g: &Self| {
            g.state
                .borrow()
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
                .copy(py)
        };
        let a = Encoding::new(
            py,
            &copy(self)?,
            node_key.as_ref(),
            edge_key.as_ref(),
            half_edge_key.as_ref(),
        )?;
        let b = Encoding::new(
            py,
            &copy(other)?,
            node_key.as_ref(),
            edge_key.as_ref(),
            half_edge_key.as_ref(),
        )?;
        Ok(py.detach(|| a.graph.is_isomorphic(&b.graph)))
    }
    /// Sort edges canonically while fixing every vertex.
    #[pyo3(signature=(*, edge_key=None, half_edge_key=None))]
    fn canonize_edges(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] edge_key: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] half_edge_key: Option<Py<PyAny>>,
    ) -> PyResult<()> {
        let mut state = self
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy(py)?;
        let node_ids = state
            .graph
            .iter_nodes()
            .map(|(index, _, _)| index)
            .collect::<Vec<_>>();
        for index in node_ids {
            state.graph[index].data = py.None();
        }
        let mut encoded =
            Encoding::new(py, &state, None, edge_key.as_ref(), half_edge_key.as_ref())?;
        for (index, &v) in encoded.nodes.iter().enumerate() {
            let _ = encoded.graph.set_node_data(v, (0, Atom::num(index as i64)));
        }
        let canonical = py.detach(|| encoded.graph.canonize());
        self.check_revision(state.revision)?;
        self.reorder_edges(Encoding::order(&encoded.edges, &canonical.vertex_map))
    }
    /// Generate connected labeled graphs and their automorphism-group sizes.
    #[staticmethod]
    #[pyo3(signature=(external_edges,vertex_signatures,*,max_vertices=None,max_loops=None,max_bridges=None,allow_self_loops=None,allow_zero_flow_edges=None,filter_fn=None,progress_fn=None))]
    #[allow(clippy::too_many_arguments)]
    #[gen_stub(override_return_type(type_repr = "list[tuple[Graph, int]]"))]
    fn generate(
        py: Python<'_>,
        external_edges: Vec<(ConvertibleToExpression, PyEdgeSignature)>,
        vertex_signatures: Vec<Vec<PyEdgeSignature>>,
        max_vertices: Option<usize>,
        max_loops: Option<usize>,
        max_bridges: Option<usize>,
        allow_self_loops: Option<bool>,
        allow_zero_flow_edges: Option<bool>,
        #[gen_stub(override_type(type_repr="typing.Callable[[Graph, int], bool] | None", imports=("typing")))]
        filter_fn: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[Graph], bool] | None", imports=("typing")))]
        progress_fn: Option<Py<PyAny>>,
    ) -> PyResult<Vec<(Self, Py<PyAny>)>> {
        if max_vertices.is_none() && max_loops.is_none() {
            return Err(PyValueError::new_err("set max_vertices or max_loops"));
        }
        let external = external_edges
            .into_iter()
            .map(|(v, e)| (v.to_expression().expr, e.edge))
            .collect::<Vec<_>>();
        let signatures = vertex_signatures
            .into_iter()
            .map(|s| s.into_iter().map(|e| e.edge).collect())
            .collect::<Vec<_>>();
        let mut settings = GenerationSettings::new();
        if let Some(v) = max_vertices {
            settings = settings.max_vertices(v);
        }
        if let Some(v) = max_loops {
            settings = settings.max_loops(v);
        }
        if let Some(v) = max_bridges {
            settings = settings.max_bridges(v);
        }
        if let Some(v) = allow_self_loops {
            settings = settings.allow_self_loops(v);
        }
        if let Some(v) = allow_zero_flow_edges {
            settings = settings.allow_zero_flow_edges(v);
        }
        let error: Arc<Mutex<Option<PyErr>>> = Arc::new(Mutex::new(None));
        if let Some(callback) = filter_fn {
            let error = error.clone();
            settings = settings.filter_fn(Box::new(move |g, n| {
                Python::attach(|py| {
                    let result = Self::from_symbolica(py, g)
                        .and_then(|g| callback.call1(py, (g, n))?.extract::<bool>(py));
                    match result {
                        Ok(v) => v,
                        Err(e) => {
                            *error.lock().unwrap() = Some(e);
                            false
                        }
                    }
                })
            }));
        }
        if let Some(callback) = progress_fn {
            let error = error.clone();
            settings = settings.progress_fn(Box::new(move |g| {
                Python::attach(|py| {
                    let result = Self::from_symbolica(py, g)
                        .and_then(|g| callback.call1(py, (g,))?.extract::<bool>(py));
                    match result {
                        Ok(v) => v,
                        Err(e) => {
                            *error.lock().unwrap() = Some(e);
                            false
                        }
                    }
                })
            }));
        }
        let abort = error.clone();
        settings = settings.abort_check(Box::new(move || {
            if abort.lock().unwrap().is_some() {
                return true;
            }
            Python::attach(|py| match py.check_signals() {
                Ok(()) => false,
                Err(e) => {
                    *abort.lock().unwrap() = Some(e);
                    true
                }
            })
        }));
        let generated = py.detach(|| {
            Graph::generate(&external, &signatures, settings).unwrap_or_else(|partial| partial)
        });
        if let Some(error) = error.lock().unwrap().take() {
            if !error.is_instance_of::<pyo3::exceptions::PyKeyboardInterrupt>(py) {
                return Err(error);
            }
        }
        let mut generated = generated.into_iter().collect::<Vec<_>>();
        generated.sort_by(|a, b| a.0.cmp(&b.0));
        generated
            .into_iter()
            .map(|(graph, size)| Ok((Self::from_symbolica(py, &graph)?, python_integer(py, size)?)))
            .collect()
    }
    /// Export a Mermaid flowchart, with escaped display labels and dangling terminals.
    fn to_mermaid(&self, py: Python<'_>) -> PyResult<String> {
        let state = self
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy(py)?;
        let mut lines = vec!["flowchart LR".to_owned()];
        let label = |data: &Py<PyAny>| -> PyResult<String> {
            Ok(data
                .bind(py)
                .str()?
                .to_string()
                .replace('&', "&amp;")
                .replace('"', "&quot;")
                .replace('<', "&lt;")
                .replace('>', "&gt;")
                .replace('\n', " "))
        };
        for (i, _, node) in state.graph.iter_nodes() {
            lines.push(format!("  n{}[\"{}\"]", i.0, label(&node.data)?));
        }
        for (pair, index, edge) in state.graph.iter_edges() {
            let (mut a, mut b) = match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => (
                    format!("n{}", state.graph.node_id(source).0),
                    format!("n{}", state.graph.node_id(sink).0),
                ),
                HedgePair::Unpaired { hedge, flow } => {
                    let n = format!("n{}", state.graph.node_id(hedge).0);
                    let t = format!("e{}", index.0);
                    lines.push(format!("  {t}(( ))"));
                    if flow == Flow::Source {
                        (n, t)
                    } else {
                        (t, n)
                    }
                }
            };
            if edge.orientation == Orientation::Reversed {
                std::mem::swap(&mut a, &mut b);
            }
            let arrow = if edge.orientation == Orientation::Undirected {
                "---"
            } else {
                "-->"
            };
            lines.push(format!(
                "  {a} {arrow}|\"{}\"| {b}",
                label(&edge.data.data)?
            ));
        }
        Ok(lines.join("\n"))
    }
}
