//! Graphica algorithms over a colored graph on the Linnet nodes.
use super::*;
use std::sync::{Arc, Mutex};
use symbolica::{
    api::python::{ConvertibleToExpression, PythonExpression},
    atom::Atom,
    domains::integer::Integer,
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

/// `(edge key, tail half-edge key, head half-edge key, direction)` as dense ranks. For a
/// paired edge the direction is 0 when undirected and 1 when oriented from tail to head; for a
/// dangling half-edge (head key 0) it is 0 when undirected, 1 out of and 2 into its node.
type EdgeKey = (u32, u32, u32, u8);
/// Graphica graph on the Linnet nodes. Dangling half-edges fold into the node colour, and
/// the half-edge keys of a paired edge into its edge colour.
type Folded = Graph<(u32, Vec<EdgeKey>), EdgeKey>;
type CanonicalGraph = (PyGraph, Vec<usize>, Py<PyAny>, Vec<usize>);

/// A paired edge or a dangling half-edge with its ends ordered independently of the stored
/// source/sink roles.
struct Normal {
    tail: Hedge,
    head: Option<Hedge>,
    key: EdgeKey,
}

/// Element keys indexed by node, edge and half-edge.
struct Keys<T> {
    nodes: Vec<T>,
    edges: Vec<T>,
    hedges: Vec<T>,
}
impl Keys<Atom> {
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
    fn read<'a>(
        py: Python<'_>,
        data: impl Iterator<Item = &'a Py<PyAny>>,
        callback: Option<&Py<PyAny>>,
    ) -> PyResult<Vec<Atom>> {
        data.map(|data| Self::key(py, data, callback)).collect()
    }
    /// Edge and half-edge keys only; node payloads are not read.
    fn edges(
        py: Python<'_>,
        graph: &PyHedgeGraph,
        edge_key: Option<&Py<PyAny>>,
        half_edge_key: Option<&Py<PyAny>>,
    ) -> PyResult<Self> {
        Ok(Self {
            nodes: vec![],
            edges: Self::read(
                py,
                graph.iter_edges().map(|(_, _, e)| &e.data.data),
                edge_key,
            )?,
            hedges: Self::read(py, graph.iter_hedges().map(|(_, h)| &h.data), half_edge_key)?,
        })
    }
    fn new(
        py: Python<'_>,
        graph: &PyHedgeGraph,
        node_key: Option<&Py<PyAny>>,
        edge_key: Option<&Py<PyAny>>,
        half_edge_key: Option<&Py<PyAny>>,
    ) -> PyResult<Self> {
        Ok(Self {
            nodes: Self::read(py, graph.iter_nodes().map(|(_, _, n)| &n.data), node_key)?,
            ..Self::edges(py, graph, edge_key, half_edge_key)?
        })
    }
    /// Dense ranks in the intrinsic key order, shared by all graphs, so that the encoding
    /// never depends on the input order.
    fn rank<const N: usize>(keys: [Self; N]) -> [Keys<u32>; N] {
        let mut sorted = keys
            .iter()
            .flat_map(|k| k.nodes.iter().chain(&k.edges).chain(&k.hedges))
            .collect::<Vec<_>>();
        sorted.sort_unstable();
        sorted.dedup();
        let rank = |keys: &[Atom]| {
            keys.iter()
                .map(|k| sorted.binary_search(&k).expect("key was collected") as u32)
                .collect()
        };
        keys.each_ref().map(|k| Keys {
            nodes: rank(&k.nodes),
            edges: rank(&k.edges),
            hedges: rank(&k.hedges),
        })
    }
}
impl Keys<u32> {
    /// Order the ends of an edge by its physical direction, else by half-edge key, else by
    /// node `position`.
    fn normal(
        &self,
        graph: &PyHedgeGraph,
        pair: HedgePair,
        edge: EdgeIndex,
        orientation: Orientation,
        position: impl Fn(NodeIndex) -> usize,
    ) -> Normal {
        let data = self.edges[edge.0];
        match pair {
            HedgePair::Unpaired { hedge, flow } => Normal {
                tail: hedge,
                head: None,
                key: (
                    data,
                    self.hedges[hedge.0],
                    0,
                    match (orientation, flow) {
                        (Orientation::Undirected, _) => 0,
                        (Orientation::Default, Flow::Source)
                        | (Orientation::Reversed, Flow::Sink) => 1,
                        _ => 2,
                    },
                ),
            },
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                let rank = |h: Hedge| (self.hedges[h.0], position(graph.node_id(h)));
                let (tail, head) = match orientation {
                    Orientation::Default => (source, sink),
                    Orientation::Reversed => (sink, source),
                    Orientation::Undirected if rank(source) <= rank(sink) => (source, sink),
                    Orientation::Undirected => (sink, source),
                };
                Normal {
                    tail,
                    head: Some(head),
                    key: (
                        data,
                        self.hedges[tail.0],
                        self.hedges[head.0],
                        (orientation != Orientation::Undirected).into(),
                    ),
                }
            }
        }
    }
    /// The folded encoding, and the permutations of identical dangling half-edges at a node
    /// that it hides from the automorphism group.
    fn fold(&self, graph: &PyHedgeGraph) -> (Folded, Integer) {
        let mut legs = vec![vec![]; graph.n_nodes()];
        let mut links = vec![];
        for (pair, index, edge) in graph.iter_edges() {
            let normal = self.normal(graph, pair, index, edge.orientation, |n| n.0);
            let tail = graph.node_id(normal.tail).0;
            match normal.head {
                None => legs[tail].push(normal.key),
                Some(head) => links.push((tail, graph.node_id(head).0, normal.key)),
            }
        }
        let mut folded = Folded::new();
        let mut hidden = Integer::one();
        for (&node, mut legs) in self.nodes.iter().zip(legs) {
            legs.sort_unstable();
            for run in legs.chunk_by(|a, b| a == b) {
                hidden *= Integer::factorial(run.len() as u32);
            }
            folded.add_node((node, legs));
        }
        for (tail, head, key) in links {
            // Distinct end keys orient an undirected edge, so swapping its ends is no symmetry.
            folded
                .add_edge(tail, head, key.3 != 0 || key.1 != key.2, key)
                .expect("folded nodes exist");
        }
        (folded, hidden)
    }
    /// Edge order by normalized endpoint positions and keys.
    fn edge_order(
        &self,
        graph: &PyHedgeGraph,
        position: impl Fn(NodeIndex) -> usize + Copy,
    ) -> Vec<usize> {
        let mut order = graph
            .iter_edges()
            .map(|(pair, index, edge)| {
                let normal = self.normal(graph, pair, index, edge.orientation, position);
                let end = |h: Hedge| position(graph.node_id(h));
                (
                    (end(normal.tail), normal.head.map(end), normal.key),
                    index.0,
                )
            })
            .collect::<Vec<_>>();
        order.sort_unstable();
        order.into_iter().map(|(_, index)| index).collect()
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
    /// The group counts permutations of identical dangling half-edges at a node, and swaps
    /// the ends of an undirected edge only when they have equal half-edge keys.
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
        let mut state = self
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy(py)?;
        let [keys] = Keys::rank([Keys::new(
            py,
            &state.graph,
            node_key.as_ref(),
            edge_key.as_ref(),
            half_edge_key.as_ref(),
        )?]);
        let (folded, hidden) = keys.fold(&state.graph);
        let canonical = py.detach(|| folded.canonize());
        let position = |n: NodeIndex| canonical.vertex_map[n.0];
        let graph = &mut state.graph;
        let edge_order = keys.edge_order(graph, position);
        // Store every edge in its normal direction, so the result only depends on its class.
        let flips = graph
            .iter_edges()
            .filter_map(|(pair, index, edge)| {
                let normal = keys.normal(graph, pair, index, edge.orientation, position);
                let (hedge, flow) = match pair {
                    HedgePair::Unpaired { hedge, flow } => (hedge, flow),
                    HedgePair::Paired { source, .. } | HedgePair::Split { source, .. } => {
                        (source, Flow::Source)
                    }
                };
                let wanted = if normal.tail == hedge && normal.key.3 != 2 {
                    Flow::Source
                } else {
                    Flow::Sink
                };
                let flip = wanted != flow || edge.orientation == Orientation::Reversed;
                flip.then_some((index, hedge, wanted, edge.orientation))
            })
            .collect::<Vec<_>>();
        for (index, hedge, flow, orientation) in flips {
            graph.set_flow(hedge, flow);
            if orientation == Orientation::Reversed {
                graph.set_orientation(index, Orientation::Default);
            }
        }
        let mut node_order = vec![0; canonical.vertex_map.len()];
        for (old, &new) in canonical.vertex_map.iter().enumerate() {
            node_order[new] = old;
        }
        let result = Self::from_state(state);
        result.reorder_nodes(node_order)?;
        result.reorder_edges(edge_order)?;
        Ok((
            result,
            canonical.vertex_map,
            python_integer(py, canonical.automorphism_group_size * hidden)?,
            canonical.orbit,
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
        let (a, b) = (copy(self)?, copy(other)?);
        let read = |state: &GraphState| {
            Keys::new(
                py,
                &state.graph,
                node_key.as_ref(),
                edge_key.as_ref(),
                half_edge_key.as_ref(),
            )
        };
        let [ka, kb] = Keys::rank([read(&a)?, read(&b)?]);
        let (fa, fb) = (ka.fold(&a.graph).0, kb.fold(&b.graph).0);
        Ok(py.detach(|| fa.is_isomorphic(&fb)))
    }
    /// Sort edges by their endpoints and keys while fixing every vertex.
    #[pyo3(signature=(*, edge_key=None, half_edge_key=None))]
    fn canonize_edges(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] edge_key: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] half_edge_key: Option<Py<PyAny>>,
    ) -> PyResult<()> {
        let state = self
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy(py)?;
        let graph = &state.graph;
        let [keys] = Keys::rank([Keys::edges(
            py,
            graph,
            edge_key.as_ref(),
            half_edge_key.as_ref(),
        )?]);
        self.check_revision(state.revision)?;
        self.reorder_edges(keys.edge_order(graph, |n| n.0))
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
