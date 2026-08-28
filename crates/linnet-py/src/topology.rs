//! Python views of Linnet's native half-edge topology and algorithms.

use std::collections::BTreeSet;

use linnet::half_edge::algorithms::DirectionBasis;
use linnet::half_edge::involution::{EdgeIndex, Hedge};
use linnet::half_edge::subgraph::{
    cut::OrientedCut, cycle::Cycle, Inclusion, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps,
};
use linnet::half_edge::tree::SimpleTraversalTree;
use linnet::half_edge::NodeIndex;
use pyo3::class::gc::{PyTraverseError, PyVisit};
use pyo3::exceptions::{PyIndexError, PyReferenceError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyAny;
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};

use crate::graph::{PyEdge, PyGraph, PyHalfEdge, PyNode, PyOrientation};
use crate::native_graph::PyHedgeGraph;

/// Selects which notion of edge direction a directed algorithm follows.
#[gen_stub_pyclass_enum]
#[pyclass(from_py_object, eq, eq_int, name = "DirectionBasis")]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PyDirectionBasis {
    /// Follow the source/sink roles stored by the half-edge involution.
    #[default]
    Underlying,
    /// Apply each edge's superficial orientation and ignore undirected edges.
    Superficial,
}

impl From<PyDirectionBasis> for DirectionBasis {
    fn from(value: PyDirectionBasis) -> Self {
        match value {
            PyDirectionBasis::Underlying => Self::Underlying,
            PyDirectionBasis::Superficial => Self::Superficial,
        }
    }
}

/// A graph-bound structural selection.
#[gen_stub_pyclass]
#[pyclass(unsendable, name = "Subgraph")]
pub struct PySubgraph {
    graph: Option<Py<PyGraph>>,
    revision: u64,
    pub(crate) subgraph: SuBitGraph,
    pub(crate) isolated_nodes: BTreeSet<usize>,
}

impl PySubgraph {
    pub(crate) fn new(graph: Py<PyGraph>, revision: u64, subgraph: SuBitGraph) -> Self {
        Self::with_isolated_nodes(graph, revision, subgraph, BTreeSet::new())
    }

    pub(crate) fn with_isolated_nodes(
        graph: Py<PyGraph>,
        revision: u64,
        subgraph: SuBitGraph,
        isolated_nodes: BTreeSet<usize>,
    ) -> Self {
        Self {
            graph: Some(graph),
            revision,
            subgraph,
            isolated_nodes,
        }
    }

    fn owner<'py>(&self, py: Python<'py>) -> PyResult<&Bound<'py, PyGraph>> {
        let graph = self
            .graph
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("subgraph has been cleared"))?;
        graph.borrow(py).check_revision(self.revision)?;
        Ok(graph.bind(py))
    }

    fn ensure_owner(&self, py: Python<'_>, graph: &Py<PyGraph>, revision: u64) -> PyResult<()> {
        self.owner(py)?;
        if self.revision != revision
            || !self
                .graph
                .as_ref()
                .expect("checked")
                .bind(py)
                .is(graph.bind(py))
        {
            return Err(PyValueError::new_err(
                "subgraph belongs to a different graph revision",
            ));
        }
        Ok(())
    }

    pub(crate) fn selection_for<'a>(
        &'a self,
        py: Python<'_>,
        graph: &Py<PyGraph>,
        revision: u64,
    ) -> PyResult<(&'a SuBitGraph, &'a BTreeSet<usize>)> {
        self.ensure_owner(py, graph, revision)?;
        Ok((&self.subgraph, &self.isolated_nodes))
    }

    fn ensure_compatible(&self, py: Python<'_>, other: &Self) -> PyResult<()> {
        let graph = self.owner(py)?;
        other.ensure_owner(py, graph.as_unbound(), self.revision)
    }

    fn with_selection(
        &self,
        py: Python<'_>,
        subgraph: SuBitGraph,
        isolated_nodes: BTreeSet<usize>,
    ) -> PyResult<Self> {
        self.owner(py)?;
        Ok(Self::with_isolated_nodes(
            self.graph.as_ref().expect("checked").clone_ref(py),
            self.revision,
            subgraph,
            isolated_nodes,
        ))
    }
}

fn isolated_node_indices(graph: &PyHedgeGraph) -> BTreeSet<usize> {
    graph
        .iter_nodes()
        .filter_map(|(node, mut crown, _)| crown.next().is_none().then_some(node.0))
        .collect()
}

#[gen_stub_pymethods]
#[pymethods]
impl PySubgraph {
    #[getter]
    fn graph(&self, py: Python<'_>) -> PyResult<Py<PyGraph>> {
        self.owner(py)?;
        Ok(self.graph.as_ref().expect("checked").clone_ref(py))
    }

    #[getter]
    fn revision(&self, py: Python<'_>) -> PyResult<u64> {
        self.owner(py)?;
        Ok(self.revision)
    }

    fn __len__(&self, py: Python<'_>) -> PyResult<usize> {
        self.owner(py)?;
        Ok(self.subgraph.n_included() + self.isolated_nodes.len())
    }

    fn __bool__(&self, py: Python<'_>) -> PyResult<bool> {
        self.owner(py)?;
        Ok(!self.subgraph.is_empty() || !self.isolated_nodes.is_empty())
    }

    #[getter]
    fn size(&self, py: Python<'_>) -> PyResult<usize> {
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        Ok(self.subgraph.size()
            + isolated_node_indices(&state.as_ref().expect("checked").graph).len())
    }

    #[getter]
    fn n_half_edges(&self, py: Python<'_>) -> PyResult<usize> {
        self.owner(py)?;
        Ok(self.subgraph.n_included())
    }

    #[getter]
    fn n_isolated_nodes(&self, py: Python<'_>) -> PyResult<usize> {
        self.owner(py)?;
        Ok(self.isolated_nodes.len())
    }

    fn includes(&self, py: Python<'_>, half_edge: usize) -> PyResult<bool> {
        self.owner(py)?;
        if half_edge >= self.subgraph.size() {
            return Err(PyIndexError::new_err("half-edge index out of range"));
        }
        Ok(self.subgraph.includes(&Hedge(half_edge)))
    }

    /// Return stable indices suitable for constructing another selection on this revision.
    fn half_edge_indices(&self, py: Python<'_>) -> PyResult<Vec<usize>> {
        self.owner(py)?;
        Ok(self.subgraph.included_iter().map(|hedge| hedge.0).collect())
    }

    /// Return the explicitly selected zero-crown node indices.
    fn isolated_node_indices(&self, py: Python<'_>) -> PyResult<Vec<usize>> {
        self.owner(py)?;
        Ok(self.isolated_nodes.iter().copied().collect())
    }

    fn includes_node(&self, py: Python<'_>, node: usize) -> PyResult<bool> {
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        let graph = &state.as_ref().expect("checked").graph;
        if node >= graph.n_nodes() {
            return Err(PyIndexError::new_err("node index out of range"));
        }
        Ok(self.isolated_nodes.contains(&node)
            || graph
                .iter_crown(NodeIndex(node))
                .any(|hedge| self.subgraph.includes(&hedge)))
    }

    fn to_half_edges(&self, py: Python<'_>) -> PyResult<Vec<Py<PyHalfEdge>>> {
        let graph = self.owner(py)?;
        self.subgraph
            .included_iter()
            .map(|hedge| {
                Py::new(
                    py,
                    PyHalfEdge::new(graph.as_unbound().clone_ref(py), hedge.0, self.revision),
                )
            })
            .collect()
    }

    fn is_subset(&self, py: Python<'_>, other: &Self) -> PyResult<bool> {
        self.ensure_compatible(py, other)?;
        Ok(other.subgraph.includes(&self.subgraph)
            && self.isolated_nodes.is_subset(&other.isolated_nodes))
    }

    fn is_superset(&self, py: Python<'_>, other: &Self) -> PyResult<bool> {
        self.ensure_compatible(py, other)?;
        Ok(self.subgraph.includes(&other.subgraph)
            && other.isolated_nodes.is_subset(&self.isolated_nodes))
    }

    fn is_disjoint(&self, py: Python<'_>, other: &Self) -> PyResult<bool> {
        self.ensure_compatible(py, other)?;
        Ok(!self.subgraph.intersects(&other.subgraph)
            && self.isolated_nodes.is_disjoint(&other.isolated_nodes))
    }

    fn union(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.ensure_compatible(py, other)?;
        self.with_selection(
            py,
            self.subgraph.union(&other.subgraph),
            self.isolated_nodes
                .union(&other.isolated_nodes)
                .copied()
                .collect(),
        )
    }

    fn intersection(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.ensure_compatible(py, other)?;
        self.with_selection(
            py,
            self.subgraph.intersection(&other.subgraph),
            self.isolated_nodes
                .intersection(&other.isolated_nodes)
                .copied()
                .collect(),
        )
    }

    fn symmetric_difference(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.ensure_compatible(py, other)?;
        self.with_selection(
            py,
            self.subgraph.sym_diff(&other.subgraph),
            self.isolated_nodes
                .symmetric_difference(&other.isolated_nodes)
                .copied()
                .collect(),
        )
    }

    fn difference(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.ensure_compatible(py, other)?;
        self.with_selection(
            py,
            self.subgraph.subtract(&other.subgraph),
            self.isolated_nodes
                .difference(&other.isolated_nodes)
                .copied()
                .collect(),
        )
    }

    fn complement(&self, py: Python<'_>) -> PyResult<Self> {
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        let filter = state
            .as_ref()
            .expect("checked")
            .graph
            .complement(&self.subgraph);
        let isolated_nodes = isolated_node_indices(&state.as_ref().expect("checked").graph)
            .difference(&self.isolated_nodes)
            .copied()
            .collect();
        drop(state);
        drop(graph);
        self.with_selection(py, filter, isolated_nodes)
    }

    fn __contains__(&self, py: Python<'_>, half_edge: usize) -> PyResult<bool> {
        self.includes(py, half_edge)
    }

    fn __or__(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.union(py, other)
    }

    fn __and__(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.intersection(py, other)
    }

    fn __xor__(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.symmetric_difference(py, other)
    }

    fn __sub__(&self, py: Python<'_>, other: &Self) -> PyResult<Self> {
        self.difference(py, other)
    }

    fn __invert__(&self, py: Python<'_>) -> PyResult<Self> {
        self.complement(py)
    }

    fn __eq__(&self, py: Python<'_>, other: &Self) -> PyResult<bool> {
        let graph = self.owner(py)?;
        let other_graph = other.owner(py)?;
        Ok(graph.is(other_graph)
            && self.subgraph == other.subgraph
            && self.isolated_nodes == other.isolated_nodes)
    }

    fn __le__(&self, py: Python<'_>, other: &Self) -> PyResult<bool> {
        self.is_subset(py, other)
    }

    fn __lt__(&self, py: Python<'_>, other: &Self) -> PyResult<bool> {
        Ok(self.is_subset(py, other)?
            && (self.subgraph != other.subgraph || self.isolated_nodes != other.isolated_nodes))
    }

    fn __ge__(&self, py: Python<'_>, other: &Self) -> PyResult<bool> {
        self.is_superset(py, other)
    }

    fn __gt__(&self, py: Python<'_>, other: &Self) -> PyResult<bool> {
        Ok(self.is_superset(py, other)?
            && (self.subgraph != other.subgraph || self.isolated_nodes != other.isolated_nodes))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.owner(py)?;
        let size = self.size(py)?;
        Ok(format!(
            "Subgraph(half_edges={}, isolated_nodes={}, size={})",
            self.subgraph.n_included(),
            self.isolated_nodes.len(),
            size,
        ))
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.graph)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.graph = None;
    }
}

/// A cycle in a particular graph revision.
#[gen_stub_pyclass]
#[pyclass(unsendable, name = "Cycle")]
pub struct PyCycle {
    graph: Option<Py<PyGraph>>,
    revision: u64,
    cycle: Cycle,
}

impl PyCycle {
    fn owner<'py>(&self, py: Python<'py>) -> PyResult<&Bound<'py, PyGraph>> {
        let graph = self
            .graph
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("cycle has been cleared"))?;
        graph.borrow(py).check_revision(self.revision)?;
        Ok(graph.bind(py))
    }

    fn ensure_owner(&self, py: Python<'_>, graph: &Py<PyGraph>, revision: u64) -> PyResult<()> {
        let owner = self.owner(py)?;
        if self.revision != revision || !owner.is(graph.bind(py)) {
            return Err(PyValueError::new_err(
                "cycle belongs to a different graph revision",
            ));
        }
        Ok(())
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl PyCycle {
    #[getter]
    fn filter(&self, py: Python<'_>) -> PyResult<PySubgraph> {
        self.owner(py)?;
        Ok(PySubgraph::new(
            self.graph.as_ref().expect("checked").clone_ref(py),
            self.revision,
            self.cycle.filter.clone(),
        ))
    }

    #[getter]
    fn loop_count(&self, py: Python<'_>) -> PyResult<Option<usize>> {
        self.filter(py)?;
        Ok(self.cycle.loop_count)
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.filter(py)?;
        Ok(format!(
            "Cycle(included={}, loop_count={:?})",
            self.cycle.filter.n_included(),
            self.cycle.loop_count
        ))
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.graph)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.graph = None;
    }
}

/// An oriented cut in a particular graph revision.
#[gen_stub_pyclass]
#[pyclass(unsendable, name = "OrientedCut")]
pub struct PyOrientedCut {
    graph: Option<Py<PyGraph>>,
    revision: u64,
    cut: OrientedCut,
}

impl PyOrientedCut {
    fn owner<'py>(&self, py: Python<'py>) -> PyResult<&Bound<'py, PyGraph>> {
        let graph = self
            .graph
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("cut has been cleared"))?;
        graph.borrow(py).check_revision(self.revision)?;
        Ok(graph.bind(py))
    }

    fn edge_ids(&self, py: Python<'_>) -> PyResult<Vec<usize>> {
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        let mut edges = state
            .as_ref()
            .expect("checked")
            .graph
            .iter_edges_of(&self.cut.left)
            .map(|(_, edge, _)| edge.0)
            .collect::<Vec<_>>();
        edges.sort_unstable();
        Ok(edges)
    }
}

/// One source-side, oriented-boundary, target-side cut partition.
#[gen_stub_pyclass]
#[pyclass(unsendable, name = "CutPartition")]
pub struct PyCutPartition {
    graph: Option<Py<PyGraph>>,
    revision: u64,
    source_side: SuBitGraph,
    boundary: OrientedCut,
    target_side: SuBitGraph,
}

impl PyCutPartition {
    fn owner<'py>(&self, py: Python<'py>) -> PyResult<&Bound<'py, PyGraph>> {
        let graph = self
            .graph
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("cut partition has been cleared"))?;
        graph.borrow(py).check_revision(self.revision)?;
        Ok(graph.bind(py))
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl PyCutPartition {
    #[getter]
    fn source_side(&self, py: Python<'_>) -> PyResult<PySubgraph> {
        self.owner(py)?;
        Ok(PySubgraph::new(
            self.graph.as_ref().expect("checked").clone_ref(py),
            self.revision,
            self.source_side.clone(),
        ))
    }

    #[getter]
    fn boundary(&self, py: Python<'_>) -> PyResult<PyOrientedCut> {
        self.owner(py)?;
        Ok(PyOrientedCut {
            graph: Some(self.graph.as_ref().expect("checked").clone_ref(py)),
            revision: self.revision,
            cut: self.boundary.clone(),
        })
    }

    #[getter]
    fn target_side(&self, py: Python<'_>) -> PyResult<PySubgraph> {
        self.owner(py)?;
        Ok(PySubgraph::new(
            self.graph.as_ref().expect("checked").clone_ref(py),
            self.revision,
            self.target_side.clone(),
        ))
    }

    #[getter]
    fn edges(&self, py: Python<'_>) -> PyResult<Vec<Py<PyEdge>>> {
        let graph = self.owner(py)?;
        let ids = {
            let graph = graph.borrow();
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .iter_edges_of(&self.boundary.left)
                .map(|(_, edge, _)| edge.0)
                .collect::<Vec<_>>()
        };
        ids.into_iter()
            .map(|edge| {
                Py::new(
                    py,
                    PyEdge::new(
                        self.graph.as_ref().expect("checked").clone_ref(py),
                        edge,
                        self.revision,
                    ),
                )
            })
            .collect()
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.owner(py)?;
        Ok(format!(
            "CutPartition(source_half_edges={}, boundary_edges={}, target_half_edges={})",
            self.source_side.n_included(),
            self.boundary.left.n_included(),
            self.target_side.n_included(),
        ))
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.graph)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.graph = None;
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl PyOrientedCut {
    fn side(&self, py: Python<'_>, left: bool) -> PyResult<PySubgraph> {
        self.owner(py)?;
        Ok(PySubgraph::new(
            self.graph.as_ref().expect("checked").clone_ref(py),
            self.revision,
            if left {
                self.cut.left.clone()
            } else {
                self.cut.right.clone()
            },
        ))
    }

    #[getter]
    fn left(&self, py: Python<'_>) -> PyResult<PySubgraph> {
        self.side(py, true)
    }

    #[getter]
    fn right(&self, py: Python<'_>) -> PyResult<PySubgraph> {
        self.side(py, false)
    }

    /// Return the cut edges in stable graph order.
    #[getter]
    fn edges(&self, py: Python<'_>) -> PyResult<Vec<Py<PyEdge>>> {
        self.edge_ids(py)?
            .into_iter()
            .map(|edge| {
                Py::new(
                    py,
                    PyEdge::new(
                        self.graph.as_ref().expect("checked").clone_ref(py),
                        edge,
                        self.revision,
                    ),
                )
            })
            .collect()
    }

    /// Return one edge's orientation relative to this cut.
    #[pyo3(signature = (edge))]
    fn orientation(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="Edge | builtins.int | builtins.str", imports=("builtins")))]
        edge: &Bound<'_, PyAny>,
    ) -> PyResult<PyOrientation> {
        self.owner(py)?;
        let graph_object = self
            .graph
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("cut has been cleared"))?;
        let edge = resolve_edge_target(graph_object, py, self.revision, edge)?;
        if !self.edge_ids(py)?.contains(&edge) {
            return Err(PyValueError::new_err("edge is not part of this cut"));
        }
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        let pair = state.as_ref().expect("checked").graph[&EdgeIndex(edge)].1;
        Ok(self.cut.get_from_pair(pair).into())
    }

    /// Return the signed intersection count with a cycle from the same graph revision.
    fn winding_number(&self, py: Python<'_>, cycle: &PyCycle) -> PyResult<i32> {
        let graph = self
            .graph
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("cut has been cleared"))?;
        self.owner(py)?;
        cycle.ensure_owner(py, graph, self.revision)?;
        Ok(self.cut.winding_number(&cycle.cycle))
    }

    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        self.left(py)?;
        Ok(format!("OrientedCut({})", self.cut))
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.graph)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.graph = None;
    }
}

/// A graph-bound DFS or BFS traversal tree.
#[gen_stub_pyclass]
#[pyclass(unsendable, name = "TraversalTree")]
pub struct PyTraversalTree {
    graph: Option<Py<PyGraph>>,
    revision: u64,
    tree: SimpleTraversalTree,
    isolated_root: Option<usize>,
}

impl PyTraversalTree {
    fn owner<'py>(&self, py: Python<'py>) -> PyResult<&Bound<'py, PyGraph>> {
        let graph = self
            .graph
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("traversal tree has been cleared"))?;
        graph.borrow(py).check_revision(self.revision)?;
        Ok(graph.bind(py))
    }

    fn resolve_node(&self, py: Python<'_>, key: &Bound<'_, PyAny>) -> PyResult<usize> {
        self.owner(py)?;
        let graph = self
            .graph
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("traversal tree has been cleared"))?;
        let node = resolve_node_target(graph, py, self.revision, key)?;
        if self.isolated_root == Some(node)
            || (self.isolated_root.is_none() && self.tree.node_data(NodeIndex(node)).includes())
        {
            Ok(node)
        } else {
            Err(PyValueError::new_err(
                "node is not part of this traversal tree",
            ))
        }
    }

    fn node_view(&self, py: Python<'_>, node: usize) -> PyResult<Py<PyNode>> {
        Py::new(
            py,
            PyNode::new(
                self.graph.as_ref().expect("checked").clone_ref(py),
                node,
                self.revision,
            ),
        )
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl PyTraversalTree {
    #[getter]
    fn subgraph(&self, py: Python<'_>) -> PyResult<PySubgraph> {
        self.owner(py)?;
        Ok(PySubgraph::with_isolated_nodes(
            self.graph.as_ref().expect("checked").clone_ref(py),
            self.revision,
            self.tree.tree_subgraph.clone(),
            self.isolated_root.into_iter().collect(),
        ))
    }

    #[getter]
    fn nodes(&self, py: Python<'_>) -> PyResult<Vec<Py<PyNode>>> {
        let graph = self.owner(py)?;
        if let Some(node) = self.isolated_root {
            return Ok(vec![Py::new(
                py,
                PyNode::new(graph.as_unbound().clone_ref(py), node, self.revision),
            )?]);
        }
        self.tree
            .node_order()
            .into_iter()
            .map(|node| {
                Py::new(
                    py,
                    PyNode::new(graph.as_unbound().clone_ref(py), node.0, self.revision),
                )
            })
            .collect()
    }

    fn covers(&self, py: Python<'_>, subgraph: &PySubgraph) -> PyResult<PySubgraph> {
        let graph = self.owner(py)?;
        subgraph.ensure_owner(py, graph.as_unbound(), self.revision)?;
        Ok(PySubgraph::with_isolated_nodes(
            self.graph.as_ref().expect("checked").clone_ref(py),
            self.revision,
            self.tree.covers(&subgraph.subgraph),
            self.isolated_root
                .filter(|node| subgraph.isolated_nodes.contains(node))
                .into_iter()
                .collect(),
        ))
    }

    #[getter]
    fn half_edges(&self, py: Python<'_>) -> PyResult<Vec<Py<PyHalfEdge>>> {
        let graph = self.owner(py)?;
        self.tree
            .tree_subgraph
            .included_iter()
            .map(|hedge| {
                Py::new(
                    py,
                    PyHalfEdge::new(graph.as_unbound().clone_ref(py), hedge.0, self.revision),
                )
            })
            .collect()
    }

    /// Return the immediate parent, or `None` for the traversal root.
    fn parent(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="Node | builtins.int | builtins.str", imports=("builtins")))]
        node: &Bound<'_, PyAny>,
    ) -> PyResult<Option<Py<PyNode>>> {
        let node = self.resolve_node(py, node)?;
        if self.isolated_root.is_some() {
            return Ok(None);
        }
        let graph = self.owner(py)?.borrow();
        let state = graph.state.borrow();
        state
            .as_ref()
            .expect("checked")
            .graph
            .traversal_parent(&self.tree, NodeIndex(node))
            .map(|parent| self.node_view(py, parent.0))
            .transpose()
    }

    /// Return the immediate children in traversal discovery order.
    fn children(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="Node | builtins.int | builtins.str", imports=("builtins")))]
        node: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Py<PyNode>>> {
        let node = self.resolve_node(py, node)?;
        if self.isolated_root.is_some() {
            return Ok(Vec::new());
        }
        let children = {
            let graph = self.owner(py)?.borrow();
            let state = graph.state.borrow();
            let graph = &state.as_ref().expect("checked").graph;
            self.tree
                .node_order()
                .into_iter()
                .filter(|candidate| {
                    graph.traversal_parent(&self.tree, *candidate) == Some(NodeIndex(node))
                })
                .collect::<Vec<_>>()
        };
        children
            .into_iter()
            .map(|child| self.node_view(py, child.0))
            .collect()
    }

    /// Return strict ancestors from the immediate parent through the root.
    fn ancestors(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="Node | builtins.int | builtins.str", imports=("builtins")))]
        node: &Bound<'_, PyAny>,
    ) -> PyResult<Vec<Py<PyNode>>> {
        let node = self.resolve_node(py, node)?;
        if self.isolated_root.is_some() {
            return Ok(Vec::new());
        }
        let ancestors = {
            let graph = self.owner(py)?.borrow();
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .traversal_ancestors(&self.tree, NodeIndex(node))
        };
        ancestors
            .into_iter()
            .map(|ancestor| self.node_view(py, ancestor.0))
            .collect()
    }

    /// Return the fundamental cycle closed by an internal edge, or `None` for a tree edge.
    fn fundamental_cycle(
        &self,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr = "_HalfEdgeTarget"))] half_edge: &Bound<'_, PyAny>,
    ) -> PyResult<Option<PyCycle>> {
        self.owner(py)?;
        let graph_object = self
            .graph
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("traversal tree has been cleared"))?;
        let hedge = Hedge(graph_object.borrow(py).resolve_half_edge(
            py,
            graph_object,
            half_edge,
        )?);
        let cycle = {
            let graph = self.owner(py)?.borrow();
            let state = graph.state.borrow();
            let graph = &state.as_ref().expect("checked").graph;
            let inverse = graph.inv(hedge);
            if inverse == hedge {
                return Err(PyValueError::new_err(
                    "an external half-edge cannot close a fundamental cycle",
                ));
            }
            if !self.tree.node_data(graph.node_id(hedge)).includes()
                || !self.tree.node_data(graph.node_id(inverse)).includes()
            {
                return Err(PyValueError::new_err(
                    "both edge endpoints must be part of this traversal tree",
                ));
            }
            graph.traversal_cycle(&self.tree, hedge)
        };
        Ok(cycle.map(|cycle| PyCycle {
            graph: Some(graph_object.clone_ref(py)),
            revision: self.revision,
            cycle,
        }))
    }

    #[gen_stub(skip)]
    fn __traverse__(&self, visit: PyVisit<'_>) -> Result<(), PyTraverseError> {
        visit.call(&self.graph)
    }

    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.graph = None;
    }
}

fn selected_filter(
    graph: &Py<PyGraph>,
    py: Python<'_>,
    revision: u64,
    subgraph: Option<&PySubgraph>,
) -> PyResult<(SuBitGraph, BTreeSet<usize>)> {
    if let Some(subgraph) = subgraph {
        let (filter, isolated_nodes) = subgraph.selection_for(py, graph, revision)?;
        return Ok((filter.clone(), isolated_nodes.clone()));
    }
    let graph = graph.borrow(py);
    let state = graph.state.borrow();
    let graph = &state.as_ref().expect("checked").graph;
    Ok((graph.full_filter(), isolated_node_indices(graph)))
}

fn resolve_node_target(
    graph: &Py<PyGraph>,
    py: Python<'_>,
    revision: u64,
    key: &Bound<'_, PyAny>,
) -> PyResult<usize> {
    if let Ok(node) = key.extract::<PyRef<'_, PyNode>>() {
        node.index_for_graph(py, graph, revision)
    } else {
        graph.borrow(py).resolve_node(key)
    }
}

fn resolve_edge_target(
    graph: &Py<PyGraph>,
    py: Python<'_>,
    revision: u64,
    key: &Bound<'_, PyAny>,
) -> PyResult<usize> {
    if let Ok(edge) = key.extract::<PyRef<'_, PyEdge>>() {
        edge.index_for_graph(py, graph, revision)
    } else {
        graph.borrow(py).resolve_edge(key)
    }
}

fn selection_includes_node(
    graph: &PyHedgeGraph,
    filter: &SuBitGraph,
    isolated_nodes: &BTreeSet<usize>,
    node: usize,
) -> bool {
    isolated_nodes.contains(&node)
        || graph
            .iter_crown(NodeIndex(node))
            .any(|hedge| filter.includes(&hedge))
}

fn check_hedge(index: Hedge, n_hedges: usize) -> PyResult<Hedge> {
    (index.0 < n_hedges)
        .then_some(index)
        .ok_or_else(|| PyIndexError::new_err("half-edge index out of range"))
}

#[gen_stub_pymethods]
#[pymethods]
impl PyGraph {
    #[pyo3(signature = ())]
    fn full_subgraph(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<PySubgraph> {
        let revision = slf.borrow(py).revision()?;
        let (filter, isolated_nodes) = selected_filter(&slf, py, revision, None)?;
        Ok(PySubgraph::with_isolated_nodes(
            slf.clone_ref(py),
            revision,
            filter,
            isolated_nodes,
        ))
    }

    #[pyo3(signature = ())]
    fn empty_subgraph(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<PySubgraph> {
        let revision = slf.borrow(py).revision()?;
        let size = slf
            .borrow(py)
            .state
            .borrow()
            .as_ref()
            .expect("checked")
            .graph
            .n_hedges();
        Ok(PySubgraph::new(
            slf.clone_ref(py),
            revision,
            SuBitGraph::empty(size),
        ))
    }

    #[pyo3(signature = (*, nodes=Vec::new(), edges=Vec::new(), half_edges=Vec::new()))]
    fn subgraph(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Sequence[builtins.int | builtins.str]", imports=("builtins", "typing")))]
        nodes: Vec<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Sequence[builtins.int | builtins.str]", imports=("builtins", "typing")))]
        edges: Vec<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Sequence[builtins.int]", imports=("builtins", "typing")))]
        half_edges: Vec<usize>,
    ) -> PyResult<PySubgraph> {
        let revision = slf.borrow(py).revision()?;
        let node_ids = nodes
            .iter()
            .map(|key| slf.borrow(py).resolve_node(key.bind(py)))
            .collect::<PyResult<Vec<_>>>()?;
        let edge_ids = edges
            .iter()
            .map(|key| slf.borrow(py).resolve_edge(key.bind(py)))
            .collect::<PyResult<Vec<_>>>()?;
        let graph_ref = slf.borrow(py);
        let state = graph_ref.state.borrow();
        let hedge_graph = &state.as_ref().expect("checked").graph;
        let n_hedges = hedge_graph.n_hedges();
        let mut filter = SuBitGraph::empty(n_hedges);
        let mut isolated_nodes = BTreeSet::new();
        for hedge in half_edges {
            filter.add(check_hedge(Hedge(hedge), n_hedges)?);
        }
        for node_id in node_ids {
            let crown = hedge_graph
                .iter_crown(NodeIndex(node_id))
                .collect::<Vec<_>>();
            if crown.is_empty() {
                isolated_nodes.insert(node_id);
            } else {
                filter.union_with_iter(crown.into_iter());
            }
        }
        for edge_id in edge_ids {
            let pair = hedge_graph
                .iter_edges()
                .find_map(|(pair, edge, _)| (edge.0 == edge_id).then_some(pair))
                .expect("resolved edge");
            filter.add(pair);
        }
        drop(state);
        drop(graph_ref);
        Ok(PySubgraph::with_isolated_nodes(
            slf.clone_ref(py),
            revision,
            filter,
            isolated_nodes,
        ))
    }

    #[pyo3(signature = (*, node=None, edge=None, half_edge=None))]
    fn filter(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Callable[[Node], builtins.bool] | None", imports=("builtins", "typing")))]
        node: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[Edge], builtins.bool] | None", imports=("builtins", "typing")))]
        edge: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdge], builtins.bool] | None", imports=("builtins", "typing")))]
        half_edge: Option<Py<PyAny>>,
    ) -> PyResult<PySubgraph> {
        let revision = slf.borrow(py).revision()?;
        let (node_crowns, edge_pairs, n_hedges) = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            let graph = &state.as_ref().expect("checked").graph;
            (
                graph
                    .iter_node_ids()
                    .map(|node| (node.0, graph.iter_crown(node).collect::<Vec<_>>()))
                    .collect::<Vec<_>>(),
                graph
                    .iter_edges()
                    .map(|(pair, edge, _)| (edge.0, pair))
                    .collect::<Vec<_>>(),
                graph.n_hedges(),
            )
        };
        let mut selected = SuBitGraph::empty(n_hedges);
        let mut isolated_nodes = BTreeSet::new();

        if let Some(predicate) = node {
            for (index, crown) in node_crowns {
                let view = Py::new(py, PyNode::new(slf.clone_ref(py), index, revision))?;
                let include = predicate.call1(py, (view,))?.bind(py).is_truthy()?;
                slf.borrow(py).check_revision(revision)?;
                if include {
                    if crown.is_empty() {
                        isolated_nodes.insert(index);
                    } else {
                        selected.union_with_iter(crown.into_iter());
                    }
                }
            }
        }
        if let Some(predicate) = edge {
            for (index, pair) in edge_pairs {
                let view = Py::new(py, PyEdge::new(slf.clone_ref(py), index, revision))?;
                let include = predicate.call1(py, (view,))?.bind(py).is_truthy()?;
                slf.borrow(py).check_revision(revision)?;
                if include {
                    selected.add(pair);
                }
            }
        }
        if let Some(predicate) = half_edge {
            for index in 0..n_hedges {
                let view = Py::new(py, PyHalfEdge::new(slf.clone_ref(py), index, revision))?;
                let include = predicate.call1(py, (view,))?.bind(py).is_truthy()?;
                slf.borrow(py).check_revision(revision)?;
                if include {
                    selected.add(Hedge(index));
                }
            }
        }
        slf.borrow(py).check_revision(revision)?;
        Ok(PySubgraph::with_isolated_nodes(
            slf.clone_ref(py),
            revision,
            selected,
            isolated_nodes,
        ))
    }

    #[pyo3(signature = (subgraph))]
    fn nodes_of(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<Vec<Py<PyNode>>> {
        let revision = slf.borrow(py).revision()?;
        subgraph.ensure_owner(py, &slf, revision)?;
        let mut ids = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .iter_nodes_of(&subgraph.subgraph)
                .map(|(node, _, _)| node.0)
                .collect::<BTreeSet<_>>()
        };
        ids.extend(subgraph.isolated_nodes.iter().copied());
        ids.into_iter()
            .map(|node| Py::new(py, PyNode::new(slf.clone_ref(py), node, revision)))
            .collect()
    }

    #[pyo3(signature = (subgraph))]
    fn edges_of(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<Vec<Py<PyEdge>>> {
        let revision = slf.borrow(py).revision()?;
        subgraph.ensure_owner(py, &slf, revision)?;
        let ids = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .iter_edges_of(&subgraph.subgraph)
                .map(|(_, edge, _)| edge)
                .collect::<Vec<_>>()
        };
        ids.into_iter()
            .map(|edge| Py::new(py, PyEdge::new(slf.clone_ref(py), edge.0, revision)))
            .collect()
    }

    #[pyo3(signature = (subgraph))]
    fn half_edges_of(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<Vec<Py<PyHalfEdge>>> {
        let revision = slf.borrow(py).revision()?;
        subgraph.ensure_owner(py, &slf, revision)?;
        subgraph.to_half_edges(py)
    }

    /// Return selected half-edges that lie on the selection's internal boundary.
    #[pyo3(signature = (subgraph))]
    fn internal_boundary(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<PySubgraph> {
        let revision = slf.borrow(py).revision()?;
        subgraph.ensure_owner(py, &slf, revision)?;
        let boundary = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .internal_boundary(&subgraph.subgraph)
        };
        Ok(PySubgraph::new(slf.clone_ref(py), revision, boundary))
    }

    /// Return all boundary half-edges incident to nodes touched by the selection.
    #[pyo3(signature = (subgraph))]
    fn boundary(slf: Py<PyGraph>, py: Python<'_>, subgraph: &PySubgraph) -> PyResult<PySubgraph> {
        let revision = slf.borrow(py).revision()?;
        subgraph.ensure_owner(py, &slf, revision)?;
        let boundary = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .boundary(&subgraph.subgraph)
        };
        Ok(PySubgraph::new(slf.clone_ref(py), revision, boundary))
    }

    /// Return dangling/external half-edges as a composable structural selection.
    #[pyo3(signature = ())]
    fn external_half_edges(slf: Py<PyGraph>, py: Python<'_>) -> PyResult<PySubgraph> {
        let revision = slf.borrow(py).revision()?;
        let external = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state.as_ref().expect("checked").graph.external_filter()
        };
        Ok(PySubgraph::new(slf.clone_ref(py), revision, external))
    }

    #[pyo3(signature = (subgraph=None))]
    fn connected_components(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<Vec<PySubgraph>> {
        let revision = slf.borrow(py).revision()?;
        let (filter, isolated_nodes) = selected_filter(&slf, py, revision, subgraph)?;
        let components = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .connected_components(&filter)
        };
        let mut components = components
            .into_iter()
            .map(|filter| PySubgraph::new(slf.clone_ref(py), revision, filter))
            .collect::<Vec<_>>();
        components.extend(isolated_nodes.into_iter().map(|node| {
            PySubgraph::with_isolated_nodes(
                slf.clone_ref(py),
                revision,
                SuBitGraph::empty(filter.size()),
                [node].into_iter().collect(),
            )
        }));
        Ok(components)
    }

    #[pyo3(signature = (subgraph=None))]
    fn count_connected_components(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<usize> {
        let revision = slf.borrow(py).revision()?;
        let (filter, isolated_nodes) = selected_filter(&slf, py, revision, subgraph)?;
        let graph = slf.borrow(py);
        let state = graph.state.borrow();
        Ok(state
            .as_ref()
            .expect("checked")
            .graph
            .count_connected_components(&filter)
            + isolated_nodes.len())
    }

    #[pyo3(signature = (subgraph=None))]
    fn is_connected(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<bool> {
        let revision = slf.borrow(py).revision()?;
        let (filter, isolated_nodes) = selected_filter(&slf, py, revision, subgraph)?;
        let graph = slf.borrow(py);
        let state = graph.state.borrow();
        let graph = &state.as_ref().expect("checked").graph;
        Ok(graph.count_connected_components(&filter) + isolated_nodes.len() <= 1)
    }

    #[pyo3(signature = (subgraph=None))]
    fn cyclomatic_number(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<usize> {
        let revision = slf.borrow(py).revision()?;
        let (filter, _) = selected_filter(&slf, py, revision, subgraph)?;
        let graph = slf.borrow(py);
        let state = graph.state.borrow();
        Ok(state
            .as_ref()
            .expect("checked")
            .graph
            .cyclotomatic_number(&filter))
    }

    /// Test directed reachability within an optional structural selection.
    ///
    /// Both endpoints must belong to the selected subgraph; otherwise this raises `ValueError`.
    #[pyo3(signature = (source, target, *, subgraph=None, direction=PyDirectionBasis::Underlying))]
    fn is_reachable(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="Node | builtins.int | builtins.str", imports=("builtins")))]
        source: &Bound<'_, PyAny>,
        #[gen_stub(override_type(type_repr="Node | builtins.int | builtins.str", imports=("builtins")))]
        target: &Bound<'_, PyAny>,
        subgraph: Option<&PySubgraph>,
        direction: PyDirectionBasis,
    ) -> PyResult<bool> {
        let revision = slf.borrow(py).revision()?;
        let (filter, isolated_nodes) = selected_filter(&slf, py, revision, subgraph)?;
        let source = resolve_node_target(&slf, py, revision, source)?;
        let target = resolve_node_target(&slf, py, revision, target)?;
        let graph = slf.borrow(py);
        let state = graph.state.borrow();
        let graph = &state.as_ref().expect("checked").graph;
        if !selection_includes_node(graph, &filter, &isolated_nodes, source) {
            return Err(PyValueError::new_err(
                "source node is not part of the selected subgraph",
            ));
        }
        if !selection_includes_node(graph, &filter, &isolated_nodes, target) {
            return Err(PyValueError::new_err(
                "target node is not part of the selected subgraph",
            ));
        }
        if source == target {
            return Ok(true);
        }
        if isolated_nodes.contains(&source) || isolated_nodes.contains(&target) {
            return Ok(false);
        }
        Ok(graph.is_reachable_of(
            &filter,
            NodeIndex(source),
            NodeIndex(target),
            direction.into(),
        ))
    }

    /// Return a deterministic topological order for the selected directed graph.
    ///
    /// Explicitly selected zero-crown nodes appear first in ascending graph order.
    #[pyo3(signature = (*, subgraph=None, direction=PyDirectionBasis::Underlying))]
    fn topological_order(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
        direction: PyDirectionBasis,
    ) -> PyResult<Vec<Py<PyNode>>> {
        let revision = slf.borrow(py).revision()?;
        let (filter, isolated_nodes) = selected_filter(&slf, py, revision, subgraph)?;
        let mut order = isolated_nodes
            .into_iter()
            .map(NodeIndex)
            .collect::<Vec<_>>();
        let directed = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .topological_order_of(&filter, direction.into())
                .map_err(PyValueError::new_err)?
        };
        order.extend(directed);
        order
            .into_iter()
            .map(|node| Py::new(py, PyNode::new(slf.clone_ref(py), node.0, revision)))
            .collect()
    }

    /// Return a new graph with directionally redundant DAG edges removed.
    #[pyo3(signature = (*, direction=PyDirectionBasis::Underlying))]
    fn transitive_reduction(
        slf: Py<PyGraph>,
        py: Python<'_>,
        direction: PyDirectionBasis,
    ) -> PyResult<PyGraph> {
        let candidate = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            let state = state
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
            let reduced = state
                .copy_graph(py)?
                .transitive_reduction(direction.into())
                .map_err(PyValueError::new_err)?;
            state.derived(py, reduced)?
        };
        Ok(PyGraph::from_state(candidate))
    }

    #[pyo3(signature = (root, *, subgraph=None, include=None))]
    fn depth_first_traverse(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="builtins.int | builtins.str", imports=("builtins")))]
        root: &Bound<'_, PyAny>,
        subgraph: Option<&PySubgraph>,
        include: Option<usize>,
    ) -> PyResult<PyTraversalTree> {
        Self::traverse(slf, py, root, subgraph, include, true)
    }

    #[pyo3(signature = (root, *, subgraph=None, include=None))]
    fn breadth_first_traverse(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="builtins.int | builtins.str", imports=("builtins")))]
        root: &Bound<'_, PyAny>,
        subgraph: Option<&PySubgraph>,
        include: Option<usize>,
    ) -> PyResult<PyTraversalTree> {
        Self::traverse(slf, py, root, subgraph, include, false)
    }

    #[pyo3(signature = (subgraph=None))]
    fn bridges(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<PySubgraph> {
        let revision = slf.borrow(py).revision()?;
        let (filter, _) = selected_filter(&slf, py, revision, subgraph)?;
        let bridges = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state.as_ref().expect("checked").graph.bridges_of(&filter)
        };
        Ok(PySubgraph::new(slf.clone_ref(py), revision, bridges))
    }

    #[pyo3(signature = (subgraph=None))]
    fn cycle_basis(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<(Vec<PyCycle>, PySubgraph)> {
        let revision = slf.borrow(py).revision()?;
        let (filter, _) = selected_filter(&slf, py, revision, subgraph)?;
        let (cycles, covered) = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .cycle_basis_of(&filter)
        };
        Ok((
            cycles
                .into_iter()
                .map(|cycle| PyCycle {
                    graph: Some(slf.clone_ref(py)),
                    revision,
                    cycle,
                })
                .collect(),
            PySubgraph::new(slf.clone_ref(py), revision, covered),
        ))
    }

    #[pyo3(signature = (subgraph=None))]
    fn all_spanning_forests(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<Vec<PySubgraph>> {
        let revision = slf.borrow(py).revision()?;
        let (filter, isolated_nodes) = selected_filter(&slf, py, revision, subgraph)?;
        let mut forests = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .all_spanning_forests_of(&filter)
        };
        if forests.is_empty() && !isolated_nodes.is_empty() {
            forests.push(SuBitGraph::empty(filter.size()));
        }
        Ok(forests
            .into_iter()
            .map(|forest| {
                PySubgraph::with_isolated_nodes(
                    slf.clone_ref(py),
                    revision,
                    forest,
                    isolated_nodes.clone(),
                )
            })
            .collect())
    }

    /// Enumerate inclusion-minimal cutsets with an inclusive boundary-size range.
    ///
    /// Bounds default to `[1, unbounded]`; `min_size=0` is rejected because a bond is non-empty.
    #[pyo3(signature = (*, subgraph=None, min_size=None, max_size=None))]
    fn all_bonds(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
        min_size: Option<usize>,
        max_size: Option<usize>,
    ) -> PyResult<Vec<PySubgraph>> {
        let revision = slf.borrow(py).revision()?;
        let (filter, isolated_nodes) = selected_filter(&slf, py, revision, subgraph)?;
        let minimum = min_size.unwrap_or(1);
        if minimum == 0 {
            return Err(PyValueError::new_err("min_size must be at least 1"));
        }
        let maximum = max_size.unwrap_or(usize::MAX);
        if minimum > maximum {
            return Err(PyValueError::new_err(
                "min_size must be less than or equal to max_size",
            ));
        }
        let mut bonds = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .all_bonds_of(&filter, minimum, maximum)
        };
        // Isolated nodes have no edge boundary and therefore contribute no bond.
        drop(isolated_nodes);
        bonds.sort_by_key(|bond| {
            bond.included_iter()
                .map(|half_edge| half_edge.0)
                .collect::<Vec<_>>()
        });
        Ok(bonds
            .into_iter()
            .map(|bond| PySubgraph::new(slf.clone_ref(py), revision, bond))
            .collect())
    }

    /// Find one native inclusion-minimal cutset within an inclusive size range.
    #[pyo3(signature = (*, subgraph=None, min_size=None, max_size=None))]
    fn find_bond(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: Option<&PySubgraph>,
        min_size: Option<usize>,
        max_size: Option<usize>,
    ) -> PyResult<Option<PySubgraph>> {
        let revision = slf.borrow(py).revision()?;
        let (filter, _) = selected_filter(&slf, py, revision, subgraph)?;
        let minimum = min_size.unwrap_or(1);
        if minimum == 0 {
            return Err(PyValueError::new_err("min_size must be at least 1"));
        }
        let maximum = max_size.unwrap_or(usize::MAX);
        if minimum > maximum {
            return Err(PyValueError::new_err(
                "min_size must be less than or equal to max_size",
            ));
        }
        let bond = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .a_bond_of(&filter, minimum, maximum)
        };
        Ok(bond.map(|bond| PySubgraph::new(slf.clone_ref(py), revision, bond)))
    }

    /// Enumerate circuits only when the complete cycle-space candidate count fits the bound.
    ///
    /// `max_results` bounds all non-empty combinations of a cycle basis, not only
    /// the combinations that ultimately form circuits.
    #[pyo3(signature = (*, max_results, subgraph=None))]
    fn all_cycles(
        slf: Py<PyGraph>,
        py: Python<'_>,
        max_results: usize,
        subgraph: Option<&PySubgraph>,
    ) -> PyResult<Vec<PyCycle>> {
        let revision = slf.borrow(py).revision()?;
        let (filter, _) = selected_filter(&slf, py, revision, subgraph)?;
        let rank = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .cyclotomatic_number(&filter)
        };
        let Some(candidates) = (rank < usize::BITS as usize).then(|| (1usize << rank) - 1) else {
            return Err(PyValueError::new_err(format!(
                "cycle-space powerset at rank {rank} is not representable or safely enumerable"
            )));
        };
        if candidates > max_results {
            return Err(PyValueError::new_err(format!(
                "cycle space has {candidates} candidate combinations, exceeding max_results={max_results}"
            )));
        }
        let cycles = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .all_cycles_of(&filter, rank)
                .map_err(PyValueError::new_err)?
        };
        Ok(cycles
            .into_iter()
            .map(|cycle| PyCycle {
                graph: Some(slf.clone_ref(py)),
                revision,
                cycle,
            })
            .collect())
    }

    /// Enumerate tadpole components after identifying a non-empty terminal node set.
    #[pyo3(signature = (externals))]
    fn tadpoles(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Sequence[Node | builtins.int | builtins.str]", imports=("builtins", "typing")))]
        externals: Vec<Py<PyAny>>,
    ) -> PyResult<Vec<PySubgraph>> {
        let revision = slf.borrow(py).revision()?;
        let externals = externals
            .iter()
            .map(|key| resolve_node_target(&slf, py, revision, key.bind(py)).map(NodeIndex))
            .collect::<PyResult<Vec<_>>>()?;
        let mut tadpoles = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            state
                .as_ref()
                .expect("checked")
                .graph
                .tadpoles(&externals)
                .map_err(PyValueError::new_err)?
        };
        tadpoles.sort_by_key(|tadpole| {
            tadpole
                .included_iter()
                .map(|half_edge| half_edge.0)
                .collect::<Vec<_>>()
        });
        Ok(tadpoles
            .into_iter()
            .map(|tadpole| PySubgraph::new(slf.clone_ref(py), revision, tadpole))
            .collect())
    }

    /// Enumerate all native separating partitions between two disjoint, non-empty node groups.
    #[pyo3(signature = (source, target))]
    fn all_cuts(
        slf: Py<PyGraph>,
        py: Python<'_>,
        #[gen_stub(override_type(type_repr="typing.Sequence[builtins.int | builtins.str]", imports=("builtins", "typing")))]
        source: Vec<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Sequence[builtins.int | builtins.str]", imports=("builtins", "typing")))]
        target: Vec<Py<PyAny>>,
    ) -> PyResult<Vec<PyCutPartition>> {
        if source.is_empty() || target.is_empty() {
            return Err(PyValueError::new_err(
                "source and target must each contain at least one node",
            ));
        }
        let revision = slf.borrow(py).revision()?;
        let source = source
            .iter()
            .map(|key| slf.borrow(py).resolve_node(key.bind(py)).map(NodeIndex))
            .collect::<PyResult<Vec<_>>>()?;
        let target = target
            .iter()
            .map(|key| slf.borrow(py).resolve_node(key.bind(py)).map(NodeIndex))
            .collect::<PyResult<Vec<_>>>()?;
        if source.iter().any(|node| target.contains(node)) {
            return Err(PyValueError::new_err(
                "source and target node groups must be disjoint",
            ));
        }
        let mut cuts = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            let graph = &state.as_ref().expect("checked").graph;
            if source
                .iter()
                .chain(&target)
                .any(|node| graph.iter_crown(*node).next().is_none())
            {
                return Err(PyValueError::new_err(
                    "every source and target node must have an incident half-edge",
                ));
            }
            let source = graph.combine_to_single_hedgenode(&source);
            let target = graph.combine_to_single_hedgenode(&target);
            if source.hairs.is_empty() || target.hairs.is_empty() {
                return Err(PyValueError::new_err(
                    "source and target node groups must each have a boundary half-edge",
                ));
            }
            graph.all_cuts(source, target)
        };
        cuts.sort_by_key(|(source, boundary, target)| {
            (
                source
                    .included_iter()
                    .map(|half_edge| half_edge.0)
                    .collect::<Vec<_>>(),
                boundary
                    .left
                    .included_iter()
                    .map(|half_edge| half_edge.0)
                    .collect::<Vec<_>>(),
                target
                    .included_iter()
                    .map(|half_edge| half_edge.0)
                    .collect::<Vec<_>>(),
            )
        });
        Ok(cuts
            .into_iter()
            .map(|(source_side, boundary, target_side)| PyCutPartition {
                graph: Some(slf.clone_ref(py)),
                revision,
                source_side,
                boundary,
                target_side,
            })
            .collect())
    }
}

impl PyGraph {
    fn traverse(
        slf: Py<PyGraph>,
        py: Python<'_>,
        root: &Bound<'_, PyAny>,
        subgraph: Option<&PySubgraph>,
        include: Option<usize>,
        depth_first: bool,
    ) -> PyResult<PyTraversalTree> {
        let revision = slf.borrow(py).revision()?;
        let (filter, isolated_nodes) = selected_filter(&slf, py, revision, subgraph)?;
        let root = NodeIndex(slf.borrow(py).resolve_node(root)?);
        let (tree, isolated_root) = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            let graph = &state.as_ref().expect("checked").graph;
            let include = include
                .map(|hedge| check_hedge(Hedge(hedge), graph.n_hedges()))
                .transpose()?;
            if isolated_nodes.contains(&root.0) {
                if include.is_some() {
                    return Err(PyValueError::new_err(
                        "an isolated traversal root has no includable half-edge",
                    ));
                }
                (graph.empty_traversal_tree(), Some(root.0))
            } else if depth_first {
                graph
                    .depth_first_traverse(&filter, &root, include)
                    .map(|tree| (tree, None))
                    .map_err(|error| PyValueError::new_err(error.to_string()))?
            } else {
                graph
                    .breadth_first_traverse(&filter, &root, include)
                    .map(|tree| (tree, None))
                    .map_err(|error| PyValueError::new_err(error.to_string()))?
            }
        };
        Ok(PyTraversalTree {
            graph: Some(slf.clone_ref(py)),
            revision,
            tree,
            isolated_root,
        })
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PyDirectionBasis>()?;
    module.add_class::<PySubgraph>()?;
    module.add_class::<PyCycle>()?;
    module.add_class::<PyOrientedCut>()?;
    module.add_class::<PyCutPartition>()?;
    module.add_class::<PyTraversalTree>()?;
    Ok(())
}
