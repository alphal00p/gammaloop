//! Canonization with Linnet's native canonizer.
use super::keys::{ElementKeys, PayloadKeys};
use super::*;
use pyo3::exceptions::PyOverflowError;

/// `(edge key, tail half-edge key, head half-edge key, direction)` as dense ranks. For a
/// paired edge the direction is 0 when undirected and 1 when oriented from tail to head; for a
/// dangling half-edge (head key 0) it is 0 when undirected, 1 out of and 2 into its node.
type EdgeKey = (u32, u32, u32, u8);
type CanonicalGraph = (PyGraph, Vec<usize>, u128, Vec<usize>);

/// An automorphism-group order, which Linnet counts up to `u128::MAX`.
pub(super) fn automorphism_order(count: Option<u128>) -> PyResult<u128> {
    count.ok_or_else(|| PyOverflowError::new_err("the automorphism group order exceeds 2**128"))
}

/// A paired edge or a dangling half-edge with its ends ordered independently of the stored
/// source/sink roles.
struct Normal {
    tail: Hedge,
    head: Option<Hedge>,
    key: EdgeKey,
}

impl ElementKeys<u32> {
    /// The topology with the key ids or ranks as payloads.
    fn ranked(&self, graph: &PyHedgeGraph) -> HedgeGraph<u32, u32, u32> {
        let mut nodes = self.nodes.iter().copied();
        graph.map_data_ref(
            |_, _, _| nodes.next().expect("one key per node"),
            |_, edge, _, data| data.map(|_| self.edges[edge.0]),
            |hedge, _| self.hedges[hedge.0],
        )
    }
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
    /// Store every edge in its normal direction, so that the stored roles only depend on its
    /// class: Reversed becomes Default, and the tail of an undirected edge its source.
    fn normalize_directions(
        &self,
        graph: &mut PyHedgeGraph,
        position: impl Fn(NodeIndex) -> usize,
    ) {
        let flips = graph
            .iter_edges()
            .filter_map(|(pair, index, edge)| {
                let normal = self.normal(graph, pair, index, edge.orientation, &position);
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
    }
}
impl PyGraph {
    /// A copy of the state, so that key callbacks run while no borrow is held.
    fn state_copy(&self, py: Python<'_>) -> PyResult<GraphState> {
        self.state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy(py)
    }
}
#[cfg_attr(feature = "python_stubgen", pyo3_stub_gen::derive::gen_stub_pymethods)]
#[cfg_attr(not(feature = "python_stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl PyGraph {
    fn __copy__(&self, py: Python<'_>) -> PyResult<Self> {
        Ok(Self::from_state(self.state_copy(py)?))
    }
    /// Canonical graph, old-to-new vertex map, automorphism-group size, and vertex orbits
    /// indexed by canonical vertex, each mapped to the smallest canonical vertex of its orbit.
    /// Keys are the element payloads, or what the key callbacks return for them; names and
    /// drawing metadata do not affect equivalence. Keys must be hashable; keys of the same kind
    /// compare by Python equality, Symbolica expressions structurally, and keys of different
    /// kinds never compare equal. Distinct keys sort as None, numbers (`numbers.Real`), strings,
    /// bytes, tuples and expressions, then other values by type name and repr, which must be
    /// the same for equal keys and differ for distinct ones; the default repr is rejected.
    /// The group counts permutations of identical dangling half-edges at a node, and swaps
    /// the ends of an undirected edge only when they have equal half-edge keys.
    #[pyo3(signature=(*, node_key=None, edge_key=None, half_edge_key=None))]
    #[gen_stub(override_return_type(type_repr = "tuple[Graph, list[int], int, list[int]]"))]
    fn canonize(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] node_key: Option<Bound<'_, PyAny>>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] edge_key: Option<Bound<'_, PyAny>>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] half_edge_key: Option<Bound<'_, PyAny>>,
    ) -> PyResult<CanonicalGraph> {
        let mut state = self.state_copy(py)?;
        let mut keys = PayloadKeys::new(
            py,
            node_key.as_ref(),
            edge_key.as_ref(),
            half_edge_key.as_ref(),
        );
        let read = keys.read(&state.graph)?;
        let keys = keys.rank(read)?;
        let ranked = keys.ranked(&state.graph);
        let canonical = py.detach(|| ranked.canonize_by(|&v| v, |&e| e, |&h| h));
        let size = automorphism_order(canonical.automorphism_count())?;
        let position = |n: NodeIndex| canonical.canonical_node(n).0;
        let n = state.graph.n_nodes();
        let vertex_map = (0..n).map(|v| position(NodeIndex(v))).collect::<Vec<_>>();
        let mut node_order = vec![0; n];
        for (old, &new) in vertex_map.iter().enumerate() {
            node_order[new] = old;
        }
        let orbits = canonical.node_orbits();
        let mut smallest = vec![usize::MAX; n];
        for (old, &new) in vertex_map.iter().enumerate() {
            let first = &mut smallest[orbits[NodeIndex(old)].0];
            *first = (*first).min(new);
        }
        let orbit = node_order
            .iter()
            .map(|&old| smallest[orbits[NodeIndex(old)].0])
            .collect();
        let edge_order = keys.edge_order(&state.graph, position);
        keys.normalize_directions(&mut state.graph, position);
        state.reorder_nodes(node_order)?;
        state.reorder_edges(edge_order)?;
        let result = GraphState {
            revision: 0,
            element_epoch: 0,
            ..state
        };
        Ok((Self::from_state(result), vertex_map, size, orbit))
    }
    #[pyo3(signature=(other, *, node_key=None, edge_key=None, half_edge_key=None))]
    fn is_isomorphic(
        &self,
        py: Python<'_>,
        other: &Self,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] node_key: Option<Bound<'_, PyAny>>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] edge_key: Option<Bound<'_, PyAny>>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] half_edge_key: Option<Bound<'_, PyAny>>,
    ) -> PyResult<bool> {
        let (a, b) = (self.state_copy(py)?, other.state_copy(py)?);
        let mut keys = PayloadKeys::new(
            py,
            node_key.as_ref(),
            edge_key.as_ref(),
            half_edge_key.as_ref(),
        );
        // Equal ids for equal keys suffice: isomorphism does not depend on the key order.
        let (ka, kb) = (keys.read(&a.graph)?, keys.read(&b.graph)?);
        let (a, b) = (ka.ranked(&a.graph), kb.ranked(&b.graph));
        Ok(py.detach(|| a.is_isomorphic_by(&b, |&v| v, |&e| e, |&h| h)))
    }
    /// Sort edges by their endpoints and keys, ordered as in `canonize`, while fixing every
    /// vertex.
    #[pyo3(signature=(*, edge_key=None, half_edge_key=None))]
    fn canonize_edges(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] edge_key: Option<Bound<'_, PyAny>>,
        #[gen_stub(override_type(type_repr = "_LabelKey"))] half_edge_key: Option<Bound<'_, PyAny>>,
    ) -> PyResult<()> {
        let state = self.state_copy(py)?;
        let mut keys = PayloadKeys::new(py, None, edge_key.as_ref(), half_edge_key.as_ref());
        let read = keys.read_edges(&state.graph)?;
        let keys = keys.rank(read)?;
        self.check_revision(state.revision)?;
        self.reorder_edges(keys.edge_order(&state.graph, |n| n.0))
    }
    /// Export a Mermaid flowchart, with escaped display labels and dangling terminals.
    fn to_mermaid(&self, py: Python<'_>) -> PyResult<String> {
        let state = self.state_copy(py)?;
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
