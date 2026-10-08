//! Python graph transformations backed by Linnet's owning topology operations.

use std::{cell::RefCell, collections::BTreeSet};

use linnet::half_edge::builder::HedgeGraphBuilder;
use linnet::half_edge::involution::{EdgeData, Flow, Hedge, Orientation};
use linnet::half_edge::subgraph::{Inclusion, SuBitGraph, SubSetLike};
use linnet::half_edge::NodeIndex;
use pyo3::exceptions::{PyReferenceError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use crate::dot::{PyEdgeValue, PyNodeValue};
use crate::drawing::copy_drawing;
use crate::graph::{
    EdgeRecord, GraphState, IncrementalEdgeRecords, NodeRecord, PyEdge, PyEdgeSpec, PyFlow,
    PyGraph, PyHalfEdge, PyHedgeGraph, PyNode, PyNodeSpec, PyOrientation, SETTLED,
};
use crate::topology::PySubgraph;

impl PyGraph {
    pub(crate) fn selected(
        graph: &Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<(u64, SuBitGraph, BTreeSet<usize>)> {
        let revision = graph.borrow(py).revision()?;
        let (filter, isolated_nodes) = subgraph.selection_for(py, graph, revision)?;
        Ok((revision, filter.clone(), isolated_nodes.clone()))
    }

    pub(crate) fn split_isolated_nodes(
        graph: &mut PyHedgeGraph,
        isolated_nodes: &BTreeSet<usize>,
    ) -> PyHedgeGraph {
        if isolated_nodes.is_empty() {
            return HedgeGraphBuilder::new().build();
        }
        // Augmented subgraphs only store zero-crown nodes, so neither edge
        // mapping callback can be reached during this native node extraction.
        graph
            .extract_nodes(
                isolated_nodes.iter().copied().map(NodeIndex),
                |_| unreachable!("an isolated node cannot split an edge"),
                |_| unreachable!("an isolated node cannot own an edge"),
            )
            .expect(SETTLED)
            .0
    }

    pub(crate) fn split_off(
        py: Python<'_>,
        graph: &mut PyHedgeGraph,
        subgraph: &SuBitGraph,
        isolated_nodes: &BTreeSet<usize>,
    ) -> PyResult<PyHedgeGraph> {
        let isolated = Self::split_isolated_nodes(graph, isolated_nodes);

        let error = RefCell::new(None);
        let mut extracted = if subgraph.is_empty() {
            HedgeGraphBuilder::new().build()
        } else {
            graph
                .extract(
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
                .expect(SETTLED)
                .0
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

    pub(crate) fn validate_names(graph: &PyHedgeGraph) -> PyResult<()> {
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

    pub(crate) fn appended_state(
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

    pub(crate) fn joined_state(
        left: &Py<PyGraph>,
        right: &Py<PyGraph>,
        py: Python<'_>,
        matching: &Py<PyAny>,
        merge: &Py<PyAny>,
    ) -> PyResult<(u64, GraphState)> {
        let left_revision = left.borrow(py).revision()?;
        let right_revision = right.borrow(py).revision()?;
        let left_epoch = left.borrow(py).element_epoch()?;
        let right_epoch = right.borrow(py).element_epoch()?;
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
                    PyHalfEdge::new(left.clone_ref(py), left_hedge.0, left_epoch),
                );
                let right_view = Py::new(
                    py,
                    PyHalfEdge::new(right.clone_ref(py), right_index, right_epoch),
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
                    PyHalfEdge::new(left.clone_ref(py), left_index, left_epoch),
                );
                let right_view = Py::new(
                    py,
                    PyHalfEdge::new(right.clone_ref(py), right_index, right_epoch),
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
                element_epoch: left_epoch,
            },
        ))
    }

    pub(crate) fn snapshot_edge_value(
        py: Python<'_>,
        value: &Bound<'_, PyAny>,
        name: Option<String>,
    ) -> PyResult<EdgeRecord> {
        let value = value.extract::<PyRef<'_, PyEdgeValue>>().map_err(|_| {
            PyTypeError::new_err("join merge callback must return EdgeValue in item 4")
        })?;
        Self::snapshot_edge_record(py, &value, name)
    }

    pub(crate) fn snapshot_edge_record(
        py: Python<'_>,
        value: &PyEdgeValue,
        name: Option<String>,
    ) -> PyResult<EdgeRecord> {
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

    pub(crate) fn merge_result(
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

impl PyGraph {
    /// Append one declarative node without copying the existing graph or its metadata.
    pub(crate) fn add_node(
        slf: Py<PyGraph>,
        py: Python<'_>,
        spec: PyRef<'_, PyNodeSpec>,
    ) -> PyResult<Py<PyNode>> {
        let revision = slf.borrow(py).revision()?;
        let epoch = slf.borrow(py).element_epoch()?;
        let node = spec.snapshot(py)?;
        let index = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            let graph = &state.as_ref().expect("checked").graph;
            if let Some(name) = &node.name {
                if graph
                    .iter_node_ids()
                    .any(|index| graph[index].name.as_ref() == Some(name))
                {
                    return Err(PyValueError::new_err(format!(
                        "duplicate node name {name:?}"
                    )));
                }
            }
            graph.node_bound()
        };
        // Complete Python allocation and validation before changing live storage.
        let result = Py::new(py, PyNode::new(slf.clone_ref(py), index, epoch))?;
        slf.borrow(py).check_revision(revision)?;
        let graph = slf.borrow(py);
        let mut state = graph.state.borrow_mut();
        let state = state.as_mut().expect("checked");
        let added = state.graph.add_node(node);
        debug_assert_eq!(added.0, index);
        state.revision = revision + 1;
        Ok(result)
    }

    /// Append one edge in place after resolving and snapshotting its endpoints.
    pub(crate) fn add_edge(
        slf: Py<PyGraph>,
        py: Python<'_>,
        spec: PyRef<'_, PyEdgeSpec>,
    ) -> PyResult<Py<PyEdge>> {
        let revision = slf.borrow(py).revision()?;
        let epoch = slf.borrow(py).element_epoch()?;
        let IncrementalEdgeRecords {
            source,
            sink,
            edge,
            orientation,
        } = spec.incremental_records(py, &slf, epoch)?;
        let index = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            let graph = &state.as_ref().expect("checked").graph;
            if let Some(name) = &edge.name {
                if graph
                    .iter_edges()
                    .any(|(_, _, edge)| edge.data.name.as_ref() == Some(name))
                {
                    return Err(PyValueError::new_err(format!(
                        "duplicate edge name {name:?}"
                    )));
                }
            }
            graph.n_edges()
        };
        let result = Py::new(py, PyEdge::new(slf.clone_ref(py), index, epoch))?;
        slf.borrow(py).check_revision(revision)?;
        let graph = slf.borrow(py);
        let mut state = graph.state.borrow_mut();
        let state = state.as_mut().expect("checked");
        match (source, sink) {
            (Some(source), Some(sink)) => state
                .graph
                .add_pair_mut(source, sink, edge, Orientation::from(orientation))
                .map(|_| ()),
            (Some(source), None) => state
                .graph
                .add_dangling_edge_mut(source, edge, Flow::Source, Orientation::from(orientation))
                .map(|_| ()),
            (None, Some(sink)) => state
                .graph
                .add_dangling_edge_mut(sink, edge, Flow::Sink, Orientation::from(orientation))
                .map(|_| ()),
            (None, None) => unreachable!("validated edge has an endpoint"),
        }
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
        state.revision = revision + 1;
        Ok(result)
    }

    /// Split a paired edge at one half-edge, preserving the opposite edge value.
    pub(crate) fn split_edge(
        slf: Py<PyGraph>,
        py: Python<'_>,
        at: &Bound<'_, PyAny>,
        replacement: PyRef<'_, PyEdgeValue>,
        name: Option<String>,
        orientation: PyOrientation,
    ) -> PyResult<(Py<PyEdge>, Py<PyEdge>)> {
        let revision = slf.borrow(py).revision()?;
        let at = slf.borrow(py).resolve_half_edge(py, &slf, at)?;
        let opposite = {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            let graph = &state.as_ref().expect("checked").graph;
            let opposite = graph.inv(Hedge(at));
            if opposite.0 == at {
                return Err(PyValueError::new_err("cannot split a dangling edge"));
            }
            opposite
        };
        let edge = Self::snapshot_edge_record(py, &replacement, name)?;
        let mut candidate = slf
            .borrow(py)
            .state
            .borrow()
            .as_ref()
            .expect("checked")
            .copy(py)?;
        candidate
            .graph
            .split_edge(Hedge(at), EdgeData::new(edge, orientation.into()))
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        let selected = candidate.graph[&Hedge(at)].0;
        let opposite = candidate.graph[&opposite].0;
        Self::validate_names(&candidate.graph)?;
        let selected = Py::new(
            py,
            PyEdge::new(slf.clone_ref(py), selected, candidate.element_epoch + 1),
        )?;
        let opposite = Py::new(
            py,
            PyEdge::new(slf.clone_ref(py), opposite, candidate.element_epoch + 1),
        )?;
        slf.borrow(py).check_revision(revision)?;
        candidate.revision = revision + 1;
        candidate.element_epoch += 1;
        *slf.borrow(py).state.borrow_mut() = Some(candidate);
        Ok((selected, opposite))
    }

    /// Connect two dangling half-edges; the first becomes the source endpoint.
    pub(crate) fn connect(
        slf: Py<PyGraph>,
        py: Python<'_>,
        source: &Bound<'_, PyAny>,
        sink: &Bound<'_, PyAny>,
        replacement: PyRef<'_, PyEdgeValue>,
        name: Option<String>,
        orientation: PyOrientation,
    ) -> PyResult<Py<PyEdge>> {
        let revision = slf.borrow(py).revision()?;
        let source = slf.borrow(py).resolve_half_edge(py, &slf, source)?;
        let sink = slf.borrow(py).resolve_half_edge(py, &slf, sink)?;
        if source == sink {
            return Err(PyValueError::new_err(
                "cannot connect a dangling half-edge to itself",
            ));
        }
        {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            let graph = &state.as_ref().expect("checked").graph;
            if graph.inv(Hedge(source)).0 != source || graph.inv(Hedge(sink)).0 != sink {
                return Err(PyValueError::new_err(
                    "connect() requires two dangling half-edges",
                ));
            }
        }
        let edge = Self::snapshot_edge_record(py, &replacement, name)?;
        let mut candidate = slf
            .borrow(py)
            .state
            .borrow()
            .as_ref()
            .expect("checked")
            .copy(py)?;
        let edge = RefCell::new(Some(EdgeData::new(edge, orientation.into())));
        candidate
            .graph
            .connect_identities(Hedge(source), Hedge(sink), |_, _, _, _| {
                (Flow::Source, edge.take().expect("merge is called once"))
            });
        let index = candidate.graph[&Hedge(source)].0;
        Self::validate_names(&candidate.graph)?;
        let result = Py::new(
            py,
            PyEdge::new(slf.clone_ref(py), index, candidate.element_epoch + 1),
        )?;
        slf.borrow(py).check_revision(revision)?;
        candidate.revision = revision + 1;
        candidate.element_epoch += 1;
        *slf.borrow(py).state.borrow_mut() = Some(candidate);
        Ok(result)
    }

    /// Copy a structural subgraph into an independent graph.
    pub(crate) fn concretize(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<PyGraph> {
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
    pub(crate) fn extract(
        slf: Py<PyGraph>,
        py: Python<'_>,
        subgraph: &PySubgraph,
    ) -> PyResult<PyGraph> {
        let (revision, filter, isolated_nodes) = Self::selected(&slf, py, subgraph)?;
        if filter.is_empty() && isolated_nodes.is_empty() {
            let graph = slf.borrow(py);
            let state = graph.state.borrow();
            let state = state
                .as_ref()
                .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
            return Ok(PyGraph::from_state(
                state.derived(py, HedgeGraphBuilder::new().build())?,
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
        candidate.element_epoch += 1;
        *slf.borrow(py).state.borrow_mut() = Some(candidate);
        Ok(PyGraph::from_state(result))
    }

    /// Delete a structural subgraph.
    pub(crate) fn delete(slf: Py<PyGraph>, py: Python<'_>, subgraph: &PySubgraph) -> PyResult<()> {
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
        candidate.graph.delete_hedges(&filter).expect(SETTLED);
        slf.borrow(py).check_revision(revision)?;
        candidate.revision = revision + 1;
        candidate.element_epoch += 1;
        *slf.borrow(py).state.borrow_mut() = Some(candidate);
        Ok(())
    }

    /// Contract a non-empty subgraph into one replacement node.
    pub(crate) fn contract(
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
        // Contraction keeps the replacement only on a half-edge that leaves the
        // island; a closed island is deleted and the replacement added instead.
        let graph = &mut candidate.graph;
        let closed = graph
            .iter_nodes_of(&filter)
            .all(|(_, mut crown, _)| crown.all(|hedge| filter.includes(&hedge)));
        if closed {
            graph.delete_hedges(&filter).expect(SETTLED);
            graph.add_node(replacement);
        } else {
            graph.contract_subgraph(&filter, replacement);
        }
        Self::validate_names(&candidate.graph)?;
        slf.borrow(py).check_revision(revision)?;
        candidate.revision = revision + 1;
        candidate.element_epoch += 1;
        *slf.borrow(py).state.borrow_mut() = Some(candidate);
        Ok(())
    }

    /// Append another graph without matching dangling half-edges.
    pub(crate) fn append(
        slf: Py<PyGraph>,
        py: Python<'_>,
        other: Py<PyGraph>,
    ) -> PyResult<PyGraph> {
        let (_, mut candidate) = Self::appended_state(&slf, &other, py)?;
        candidate.element_epoch = 0;
        Ok(PyGraph::from_state(candidate))
    }

    /// Append another graph in place without matching dangling half-edges.
    pub(crate) fn append_mut(slf: Py<PyGraph>, py: Python<'_>, other: Py<PyGraph>) -> PyResult<()> {
        let (revision, mut candidate) = Self::appended_state(&slf, &other, py)?;
        slf.borrow(py).check_revision(revision)?;
        candidate.revision = revision + 1;
        candidate.element_epoch += 1;
        *slf.borrow(py).state.borrow_mut() = Some(candidate);
        Ok(())
    }

    /// Join dangling half-edges selected by Python callbacks.
    pub(crate) fn join(
        slf: Py<PyGraph>,
        py: Python<'_>,
        other: Py<PyGraph>,
        matching: Py<PyAny>,
        merge: Py<PyAny>,
    ) -> PyResult<PyGraph> {
        let (_, mut candidate) = Self::joined_state(&slf, &other, py, &matching, &merge)?;
        candidate.element_epoch = 0;
        Ok(PyGraph::from_state(candidate))
    }

    /// Join another graph in place while leaving the other graph unchanged.
    pub(crate) fn join_mut(
        slf: Py<PyGraph>,
        py: Python<'_>,
        other: Py<PyGraph>,
        matching: Py<PyAny>,
        merge: Py<PyAny>,
    ) -> PyResult<()> {
        let (revision, mut candidate) = Self::joined_state(&slf, &other, py, &matching, &merge)?;
        slf.borrow(py).check_revision(revision)?;
        candidate.revision = revision + 1;
        candidate.element_epoch += 1;
        *slf.borrow(py).state.borrow_mut() = Some(candidate);
        Ok(())
    }
}
