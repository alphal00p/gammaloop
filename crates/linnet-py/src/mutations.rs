//! Python graph transformations backed by Linnet's owning topology operations.

use std::{cell::RefCell, collections::BTreeSet};

use linnet::half_edge::builder::HedgeGraphBuilder;
use linnet::half_edge::involution::{EdgeData, Flow};
use linnet::half_edge::subgraph::{Inclusion, SuBitGraph, SubSetLike};
use linnet::half_edge::NodeIndex;
use pyo3::exceptions::{PyReferenceError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use pyo3_stub_gen::derive::gen_stub_pymethods;

use crate::dot::{PyEdgeValue, PyNodeValue};
use crate::drawing::copy_drawing;
use crate::graph::{
    EdgeRecord, GraphState, NodeRecord, PyFlow, PyGraph, PyHalfEdge, PyOrientation,
};
use crate::native_graph::PyHedgeGraph;
use crate::topology::PySubgraph;

impl PyGraph {
    fn selected(
        graph: &Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<(u64, SuBitGraph, BTreeSet<usize>)> {
        let revision = graph.borrow(py).revision()?;
        let (filter, isolated_nodes) = subgraph.selection_for(py, graph, revision)?;
        Ok((revision, filter.clone(), isolated_nodes.clone()))
    }

    fn split_isolated_nodes(
        graph: &mut PyHedgeGraph,
        isolated_nodes: &BTreeSet<usize>,
    ) -> PyHedgeGraph {
        if isolated_nodes.is_empty() {
            return PyHedgeGraph::empty(graph.node_store());
        }
        // Augmented subgraphs only store zero-crown nodes, so neither edge
        // mapping callback can be reached during this native node extraction.
        graph.extract_nodes(
            isolated_nodes.iter().copied().map(NodeIndex),
            |_| unreachable!("an isolated node cannot split an edge"),
            |_| unreachable!("an isolated node cannot own an edge"),
        )
    }

    fn split_off(
        py: Python<'_>,
        graph: &mut PyHedgeGraph,
        subgraph: &SuBitGraph,
        isolated_nodes: &BTreeSet<usize>,
    ) -> PyResult<PyHedgeGraph> {
        let node_store = graph.node_store();
        let isolated = Self::split_isolated_nodes(graph, isolated_nodes);

        let error = RefCell::new(None);
        let mut extracted = if subgraph.is_empty() {
            PyHedgeGraph::empty(node_store)
        } else {
            graph.extract(
                subgraph,
                |edge| EdgeData {
                    orientation: edge.orientation,
                    data: edge.data.copy(py).unwrap_or_else(|copy_error| {
                        if error.borrow().is_none() {
                            *error.borrow_mut() = Some(copy_error);
                        }
                        EdgeRecord {
                            name: edge.data.name.clone(),
                            data: edge.data.data.clone_ref(py),
                            drawing: pyo3::types::PyDict::new(py).unbind(),
                        }
                    }),
                },
                |edge| edge,
                |node| {
                    node.copy(py).unwrap_or_else(|copy_error| {
                        if error.borrow().is_none() {
                            *error.borrow_mut() = Some(copy_error);
                        }
                        NodeRecord {
                            name: node.name.clone(),
                            data: node.data.clone_ref(py),
                            drawing: pyo3::types::PyDict::new(py).unbind(),
                        }
                    })
                },
                |node| node,
            )
        };
        match error.into_inner() {
            Some(error) => Err(error),
            None => {
                if !isolated_nodes.is_empty() {
                    extracted
                        .append_disconnected_mut(isolated)
                        .map_err(|error| PyValueError::new_err(error.to_string()))?;
                }
                Ok(extracted)
            }
        }
    }

    fn validate_names(graph: &PyHedgeGraph) -> PyResult<()> {
        let mut names = std::collections::BTreeSet::new();
        for (_, _, node) in graph.iter_nodes() {
            if let Some(name) = &node.name {
                if !names.insert(name) {
                    return Err(PyValueError::new_err(format!(
                        "duplicate node name {name:?} after graph composition"
                    )));
                }
            }
        }
        names.clear();
        for (_, _, edge) in graph.iter_edges() {
            if let Some(name) = &edge.data.name {
                if !names.insert(name) {
                    return Err(PyValueError::new_err(format!(
                        "duplicate edge name {name:?} after graph composition"
                    )));
                }
            }
        }
        Ok(())
    }

    fn appended_state(
        left: &Py<PyGraph>,
        right: &Py<PyGraph>,
        py: Python<'_>,
    ) -> PyResult<(u64, GraphState)> {
        let revision = left.borrow(py).revision()?;
        let mut candidate = left
            .borrow(py)
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy(py)?;
        let other = right
            .borrow(py)
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy_graph(py)?;
        candidate
            .graph
            .append_disconnected_mut(other)
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        Self::validate_names(&candidate.graph)?;
        candidate.revision = 0;
        Ok((revision, candidate))
    }

    fn joined_state(
        left: &Py<PyGraph>,
        right: &Py<PyGraph>,
        py: Python<'_>,
        matching: &Py<PyAny>,
        merge: &Py<PyAny>,
    ) -> PyResult<(u64, GraphState)> {
        let left_revision = left.borrow(py).revision()?;
        let right_revision = right.borrow(py).revision()?;
        let candidate = left
            .borrow(py)
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy(py)?;
        let right_graph = right
            .borrow(py)
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy_graph(py)?;
        let GraphState {
            graph: left_graph,
            name,
            global_data,
            codec,
            render_config,
            ..
        } = candidate;
        let left_hedges = left_graph.n_hedges();
        let error = RefCell::<Option<PyErr>>::new(None);
        let matched = RefCell::<Option<(usize, usize)>>::new(None);
        let joined = left_graph.join_with_hedge_data(
            right_graph,
            |left_hedge, _, _, _, right_hedge, _, _, _| {
                if error.borrow().is_some() {
                    return false;
                }
                let Some(right_index) = right_hedge.0.checked_sub(left_hedges) else {
                    *error.borrow_mut() = Some(PyValueError::new_err(
                        "Linnet returned an invalid joined half-edge index",
                    ));
                    return false;
                };
                let left_view = Py::new(
                    py,
                    PyHalfEdge::new(left.clone_ref(py), left_hedge.0, left_revision),
                );
                let right_view = Py::new(
                    py,
                    PyHalfEdge::new(right.clone_ref(py), right_index, right_revision),
                );
                let result = left_view.and_then(|left_view| {
                    right_view.and_then(|right_view| {
                        matching.call1(py, (left_view, right_view))?.extract(py)
                    })
                });
                match result {
                    Ok(true) => {
                        *matched.borrow_mut() = Some((left_hedge.0, right_index));
                        true
                    }
                    Ok(false) => false,
                    Err(callback_error) => {
                        *error.borrow_mut() = Some(callback_error);
                        false
                    }
                }
            },
            |left_flow, left_edge, _, _right_edge| {
                if error.borrow().is_some() {
                    return (left_flow, left_edge);
                }
                let Some((left_index, right_index)) = matched.borrow_mut().take() else {
                    *error.borrow_mut() = Some(PyValueError::new_err(
                        "Linnet merged dangling edges without a matched pair",
                    ));
                    return (left_flow, left_edge);
                };
                let left_view = Py::new(
                    py,
                    PyHalfEdge::new(left.clone_ref(py), left_index, left_revision),
                );
                let right_view = Py::new(
                    py,
                    PyHalfEdge::new(right.clone_ref(py), right_index, right_revision),
                );
                let result = left_view.and_then(|left_view| {
                    right_view.and_then(|right_view| {
                        let value = merge.call1(py, (left_view, right_view))?;
                        Self::merge_result(py, value.bind(py))
                    })
                });
                match result {
                    Ok(merged) => merged,
                    Err(callback_error) => {
                        *error.borrow_mut() = Some(callback_error);
                        (left_flow, left_edge)
                    }
                }
            },
        );

        if let Some(error) = error.into_inner() {
            return Err(error);
        }
        left.borrow(py).check_revision(left_revision)?;
        right.borrow(py).check_revision(right_revision)?;
        let graph = joined.map_err(|error| PyValueError::new_err(error.to_string()))?;
        Self::validate_names(&graph)?;
        Ok((
            left_revision,
            GraphState {
                graph,
                name,
                global_data,
                codec,
                render_config,
                revision: 0,
            },
        ))
    }

    fn snapshot_edge_value(
        py: Python<'_>,
        value: &Bound<'_, PyAny>,
        name: Option<String>,
    ) -> PyResult<EdgeRecord> {
        let value = value.extract::<PyRef<'_, PyEdgeValue>>().map_err(|_| {
            PyTypeError::new_err("join merge callback must return EdgeValue in item 4")
        })?;
        Ok(EdgeRecord {
            name,
            data: value
                .data
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("EdgeValue has been cleared"))?
                .clone_ref(py),
            drawing: copy_drawing(
                py,
                value
                    .drawing
                    .as_ref()
                    .ok_or_else(|| PyReferenceError::new_err("EdgeValue has been cleared"))?,
            )?,
        })
    }

    fn merge_result(
        py: Python<'_>,
        value: &Bound<'_, PyAny>,
    ) -> PyResult<(Flow, EdgeData<EdgeRecord>)> {
        let value = value.cast::<PyTuple>().map_err(|_| {
            PyTypeError::new_err(
                "join merge callback must return (Flow, Orientation, name, EdgeValue)",
            )
        })?;
        if value.len() != 4 {
            return Err(PyTypeError::new_err(
                "join merge callback must return four items",
            ));
        }
        let flow = value.get_item(0)?.extract::<PyFlow>()?;
        let orientation = value.get_item(1)?.extract::<PyOrientation>()?;
        let name = value.get_item(2)?.extract::<Option<String>>()?;
        let edge = Self::snapshot_edge_value(py, &value.get_item(3)?, name)?;
        Ok((flow.into(), EdgeData::new(edge, orientation.into())))
    }
}

#[gen_stub_pymethods]
#[pymethods]
impl PyGraph {
    /// Copy a structural subgraph into an independent graph.
    #[pyo3(signature = (subgraph))]
    fn concretize(slf: Py<PyGraph>, py: Python<'_>, subgraph: &PySubgraph) -> PyResult<PyGraph> {
        let (_, filter, isolated_nodes) = Self::selected(&slf, py, subgraph)?;
        let state = slf.borrow(py);
        let state = state.state.borrow();
        let state = state
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        let mut graph = state.copy_graph(py)?;
        let extracted = Self::split_off(py, &mut graph, &filter, &isolated_nodes)?;
        Ok(PyGraph::from_state(state.derived(py, extracted)?))
    }

    /// Remove a structural subgraph and return it as an independent graph.
    #[pyo3(signature = (subgraph))]
    fn extract(slf: Py<PyGraph>, py: Python<'_>, subgraph: &PySubgraph) -> PyResult<PyGraph> {
        let (revision, filter, isolated_nodes) = Self::selected(&slf, py, subgraph)?;
        if filter.is_empty() && isolated_nodes.is_empty() {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            let state = state
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
            return Ok(PyGraph::from_state(
                state.derived(py, PyHedgeGraph::empty(state.graph.node_store()))?,
            ));
        }
        let mut candidate = slf
            .borrow(py)
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy(py)?;
        let extracted = Self::split_off(py, &mut candidate.graph, &filter, &isolated_nodes)?;
        let result = candidate.derived(py, extracted)?;
        slf.borrow(py).check_revision(revision)?;
        candidate.revision = revision + 1;
        *slf.borrow(py).state.borrow_mut() = Some(candidate);
        Ok(PyGraph::from_state(result))
    }

    /// Delete a structural subgraph.
    #[pyo3(signature = (subgraph))]
    fn delete(slf: Py<PyGraph>, py: Python<'_>, subgraph: &PySubgraph) -> PyResult<()> {
        let (revision, filter, isolated_nodes) = Self::selected(&slf, py, subgraph)?;
        if filter.is_empty() && isolated_nodes.is_empty() {
            return Ok(());
        }
        let mut candidate = slf
            .borrow(py)
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy(py)?;
        Self::split_isolated_nodes(&mut candidate.graph, &isolated_nodes);
        candidate.graph.delete_hedges(&filter);
        slf.borrow(py).check_revision(revision)?;
        candidate.revision = revision + 1;
        *slf.borrow(py).state.borrow_mut() = Some(candidate);
        Ok(())
    }

    /// Contract a non-empty subgraph into one replacement node.
    #[pyo3(signature = (subgraph, replacement, *, name=None))]
    fn contract(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
        replacement: PyRef<'_, PyNodeValue>,
        name: Option<String>,
    ) -> PyResult<()> {
        let (revision, filter, isolated_nodes) = Self::selected(&slf, py, subgraph)?;
        if filter.is_empty() && isolated_nodes.is_empty() {
            return Err(PyValueError::new_err("cannot contract an empty subgraph"));
        }
        let replacement = NodeRecord {
            name,
            data: replacement
                .data
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("NodeValue has been cleared"))?
                .clone_ref(py),
            drawing: copy_drawing(
                py,
                replacement
                    .drawing
                    .as_ref()
                    .ok_or_else(|| PyReferenceError::new_err("NodeValue has been cleared"))?,
            )?,
        };
        let mut candidate = slf
            .borrow(py)
            .state
            .borrow()
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
            .copy(py)?;
        Self::split_isolated_nodes(&mut candidate.graph, &isolated_nodes);
        if filter.is_empty() {
            let mut builder = HedgeGraphBuilder::new();
            builder.add_node(replacement);
            candidate
                .graph
                .append_disconnected_mut(PyHedgeGraph::build(builder, candidate.graph.node_store()))
                .map_err(|error| PyValueError::new_err(error.to_string()))?;
            Self::validate_names(&candidate.graph)?;
            slf.borrow(py).check_revision(revision)?;
            candidate.revision = revision + 1;
            *slf.borrow(py).state.borrow_mut() = Some(candidate);
            return Ok(());
        }
        // Move unrelated empty-crown nodes aside because node stores differ in
        // how contraction compaction treats them. Restore the replacement when
        // the contracted island has no boundary hairs.
        let unrelated_isolated = candidate
            .graph
            .iter_nodes()
            .filter_map(|(node, mut crown, _)| crown.next().is_none().then_some(node.0))
            .collect::<BTreeSet<_>>();
        let mut preserved = Self::split_isolated_nodes(&mut candidate.graph, &unrelated_isolated);
        let involved = candidate
            .graph
            .iter_nodes_of(&filter)
            .map(|(node, _, _)| node)
            .collect::<Vec<NodeIndex>>();
        let has_boundary = involved.iter().any(|node| {
            candidate
                .graph
                .iter_crown(*node)
                .any(|hedge| !filter.includes(&hedge))
        });
        let replacement_backup = (!has_boundary).then(|| replacement.copy(py)).transpose()?;
        candidate.graph.contract_subgraph(&filter, replacement);
        if let Some(node) = replacement_backup {
            let mut builder = HedgeGraphBuilder::new();
            builder.add_node(node);
            preserved
                .append_disconnected_mut(PyHedgeGraph::build(builder, candidate.graph.node_store()))
                .map_err(|error| PyValueError::new_err(error.to_string()))?;
        }
        if preserved.n_nodes() != 0 {
            candidate
                .graph
                .append_disconnected_mut(preserved)
                .map_err(|error| PyValueError::new_err(error.to_string()))?;
        }
        Self::validate_names(&candidate.graph)?;
        slf.borrow(py).check_revision(revision)?;
        candidate.revision = revision + 1;
        *slf.borrow(py).state.borrow_mut() = Some(candidate);
        Ok(())
    }

    /// Append another graph without matching dangling half-edges.
    #[pyo3(signature = (other))]
    fn append(slf: Py<PyGraph>, py: Python<'_>, other: Py<PyGraph>) -> PyResult<PyGraph> {
        let (_, candidate) = Self::appended_state(&slf, &other, py)?;
        Ok(PyGraph::from_state(candidate))
    }

    /// Append another graph in place without matching dangling half-edges.
    #[pyo3(signature = (other))]
    fn append_mut(slf: Py<PyGraph>, py: Python<'_>, other: Py<PyGraph>) -> PyResult<()> {
        let (revision, mut candidate) = Self::appended_state(&slf, &other, py)?;
        slf.borrow(py).check_revision(revision)?;
        candidate.revision = revision + 1;
        *slf.borrow(py).state.borrow_mut() = Some(candidate);
        Ok(())
    }

    /// Join dangling half-edges selected by Python callbacks.
    #[pyo3(signature = (other, *, matching, merge))]
    fn join(
        slf: Py<PyGraph>,
        py: Python<'_>,
        other: Py<PyGraph>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdge, HalfEdge], builtins.bool]", imports=("builtins", "typing")))]
        matching: Py<PyAny>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdge, HalfEdge], tuple[Flow, Orientation, builtins.str | None, EdgeValue]]", imports=("builtins", "typing")))]
        merge: Py<PyAny>,
    ) -> PyResult<PyGraph> {
        let (_, candidate) = Self::joined_state(&slf, &other, py, &matching, &merge)?;
        Ok(PyGraph::from_state(candidate))
    }

    /// Join another graph in place while leaving the other graph unchanged.
    #[pyo3(signature = (other, *, matching, merge))]
    fn join_mut(
        slf: Py<PyGraph>,
        py: Python<'_>,
        other: Py<PyGraph>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdge, HalfEdge], builtins.bool]", imports=("builtins", "typing")))]
        matching: Py<PyAny>,
        #[gen_stub(override_type(type_repr="typing.Callable[[HalfEdge, HalfEdge], tuple[Flow, Orientation, builtins.str | None, EdgeValue]]", imports=("builtins", "typing")))]
        merge: Py<PyAny>,
    ) -> PyResult<()> {
        let (revision, mut candidate) = Self::joined_state(&slf, &other, py, &matching, &merge)?;
        slf.borrow(py).check_revision(revision)?;
        candidate.revision = revision + 1;
        *slf.borrow(py).state.borrow_mut() = Some(candidate);
        Ok(())
    }
}
