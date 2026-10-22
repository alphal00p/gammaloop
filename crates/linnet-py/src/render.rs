use std::collections::BTreeMap;

use linnest::{
    encode_graph_spec_bytes, TypstEdgeSpec, TypstEndpointSpec, TypstGraphSpec, TypstNodeSpec,
};
use linnet::half_edge::builder::{HedgeData, HedgeGraphBuilder};
use linnet::half_edge::involution::{Flow, Hedge, HedgePair, Orientation};
use linnet::half_edge::subgraph::{Inclusion, SubSetLike};
use pyo3::exceptions::{PyReferenceError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyAny, PyDict, PyDictMethods, PyList, PyListMethods};

use crate::dot::PyGlobalData;
use crate::drawing::{element_drawing, DrawingKind};
use crate::graph::{
    EdgeRecord, GraphState, HalfEdgeRecord, NodeRecord, PyEdge, PyGraph, PyHalfEdge, PyNode,
};
use crate::native_graph::{PyHedgeGraph, PyNodeStore};
use crate::topology::PySubgraph;
use crate::typst::{
    evaluate_selector, render_config_transport, PyRenderSettings, RenderSettingsTransport,
    SelectorCallbacks,
};

fn base_elements(py: Python<'_>, graph: &PyGraph) -> PyResult<Py<PyDict>> {
    let state = graph.state.borrow();
    let state = state
        .as_ref()
        .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
    let nodes = PyList::empty(py);
    for (_, _, node) in state.graph.iter_nodes() {
        nodes.append(element_drawing(
            py,
            node.drawing.bind(py),
            DrawingKind::Node,
        )?)?;
    }
    let edges = PyList::empty(py);
    let half_edges = PyList::empty(py);
    for (_, _, edge) in state.graph.iter_edges() {
        edges.append(element_drawing(
            py,
            edge.data.drawing.bind(py),
            DrawingKind::Edge,
        )?)?;
    }
    for index in 0..state.graph.n_hedges() {
        half_edges.append(element_drawing(
            py,
            state.graph[Hedge(index)].drawing.bind(py),
            DrawingKind::HalfEdge,
        )?)?;
    }
    let elements = PyDict::new(py);
    elements.set_item("graph", PyDict::new(py))?;
    elements.set_item("nodes", nodes)?;
    elements.set_item("edges", edges)?;
    elements.set_item("hedges", half_edges)?;
    Ok(elements.unbind())
}

fn apply_selector(
    py: Python<'_>,
    selector: Option<&Py<PyAny>>,
    elements: &Bound<'_, PyList>,
    views: impl IntoIterator<Item = PyResult<Py<PyAny>>>,
    kind: DrawingKind,
    override_drawing: bool,
) -> PyResult<()> {
    let Some(selector) = selector else {
        return Ok(());
    };
    for (index, view) in views.into_iter().enumerate() {
        let Some(value) = evaluate_selector(py, selector, view?.bind(py), kind)? else {
            continue;
        };
        let patch = element_drawing(py, value.bind(py), kind)?;
        let element = elements.get_item(index)?.cast_into::<PyDict>()?;
        for (key, value) in patch.iter() {
            if override_drawing
                || element
                    .get_item(&key)?
                    .as_ref()
                    .is_none_or(|value| crate::typst::is_inherit(value))
            {
                element.set_item(key, value)?;
            }
        }
    }
    Ok(())
}

fn apply_selectors(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    elements: &Bound<'_, PyDict>,
    selectors: &SelectorCallbacks,
    override_drawing: bool,
) -> PyResult<()> {
    let graph_ref = graph.borrow(py);
    let (revision, node_count, edge_count, half_edge_count, half_edge_roles) = {
        let state = graph_ref.state.borrow();
        let state = state
            .as_ref()
            .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
        (
            state.revision,
            state.graph.n_nodes(),
            state.graph.n_edges(),
            state.graph.n_hedges(),
            (0..state.graph.n_hedges())
                .map(|index| state.graph.flow(Hedge(index)) == Flow::Source)
                .collect::<Vec<_>>(),
        )
    };
    drop(graph_ref);

    let nodes = elements
        .get_item("nodes")?
        .expect("base elements")
        .cast_into::<PyList>()?;
    apply_selector(
        py,
        selectors.node.as_ref(),
        &nodes,
        (0..node_count).map(|index| {
            Ok(Py::new(py, PyNode::new(graph.clone_ref(py), index, revision))?.into_any())
        }),
        DrawingKind::Node,
        override_drawing,
    )?;

    let edges = elements
        .get_item("edges")?
        .expect("base elements")
        .cast_into::<PyList>()?;
    apply_selector(
        py,
        selectors.edge.as_ref(),
        &edges,
        (0..edge_count).map(|index| {
            Ok(Py::new(py, PyEdge::new(graph.clone_ref(py), index, revision))?.into_any())
        }),
        DrawingKind::Edge,
        override_drawing,
    )?;

    let half_edges = elements
        .get_item("hedges")?
        .expect("base elements")
        .cast_into::<PyList>()?;
    for (index, is_source) in half_edge_roles
        .into_iter()
        .enumerate()
        .take(half_edge_count)
    {
        let view = Py::new(py, PyHalfEdge::new(graph.clone_ref(py), index, revision))?.into_any();
        let selector = if is_source {
            selectors.source.as_ref()
        } else {
            selectors.sink.as_ref()
        };
        if let Some(selector) = selector {
            let Some(value) =
                evaluate_selector(py, selector, view.bind(py), DrawingKind::HalfEdge)?
            else {
                continue;
            };
            let patch = element_drawing(py, value.bind(py), DrawingKind::HalfEdge)?;
            let element = half_edges.get_item(index)?.cast_into::<PyDict>()?;
            for (key, value) in patch.iter() {
                if override_drawing
                    || element
                        .get_item(&key)?
                        .as_ref()
                        .is_none_or(|value| crate::typst::is_inherit(value))
                {
                    element.set_item(key, value)?;
                }
            }
        }
    }
    graph.borrow(py).check_revision(revision).map_err(|_| {
        PyReferenceError::new_err("drawing selectors must not mutate graph topology")
    })?;
    Ok(())
}

fn request(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    overlay: Option<&Bound<'_, PyAny>>,
    override_drawing: bool,
) -> PyResult<(Vec<u8>, RenderSettingsTransport)> {
    let graph_ref = graph.borrow(py);
    let elements = base_elements(py, &graph_ref)?;
    let base = graph_ref
        .state
        .borrow()
        .as_ref()
        .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?
        .render_config
        .clone_ref(py);
    let initial = render_config_transport(py, &base, overlay, elements.bind(py))?;
    drop(graph_ref);
    apply_selectors(
        py,
        graph,
        elements.bind(py),
        &initial.selectors,
        override_drawing,
    )?;
    let transport = render_config_transport(py, &base, overlay, elements.bind(py))?;
    let topology = topology_spec(&graph.borrow(py))?;
    Ok((topology, transport))
}

fn topology_endpoint(graph: &PyHedgeGraph, hedge: Hedge) -> TypstEndpointSpec {
    TypstEndpointSpec {
        node: graph.node_id(hedge).0,
        statement: None,
        id: Some(hedge.0),
        data: None,
        port_label: None,
        route_points: Vec::new(),
        compass: None,
        in_subgraph: false,
    }
}

fn topology_spec(graph: &PyGraph) -> PyResult<Vec<u8>> {
    let state = graph.state.borrow();
    let state = state
        .as_ref()
        .ok_or_else(|| PyReferenceError::new_err("graph has been cleared"))?;
    let nodes = state
        .graph
        .iter_nodes()
        .map(|(index, _, node)| TypstNodeSpec {
            name: node.name.clone(),
            index: Some(index.0),
            data: None,
            pos: None,
            statements: BTreeMap::new(),
        })
        .collect();
    let mut edges = Vec::with_capacity(state.graph.n_edges());
    for (pair, index, edge) in state.graph.iter_edges() {
        let (source, sink, flow) = match pair {
            HedgePair::Paired { source, sink } => (
                Some(topology_endpoint(&state.graph, source)),
                Some(topology_endpoint(&state.graph, sink)),
                None,
            ),
            HedgePair::Unpaired { hedge, flow } => match flow {
                Flow::Source => (
                    Some(topology_endpoint(&state.graph, hedge)),
                    None,
                    Some("source".to_owned()),
                ),
                Flow::Sink => (
                    None,
                    Some(topology_endpoint(&state.graph, hedge)),
                    Some("sink".to_owned()),
                ),
            },
            HedgePair::Split { .. } => {
                return Err(PyValueError::new_err(
                    "a full graph cannot render a split edge",
                ));
            }
        };
        let orientation = match edge.orientation {
            Orientation::Default => "default",
            Orientation::Reversed => "reversed",
            Orientation::Undirected => "undirected",
        };
        edges.push(TypstEdgeSpec {
            name: edge.data.name.clone(),
            source,
            sink,
            data: None,
            orientation: Some(orientation.to_owned()),
            flow,
            id: Some(index.0),
            pos: None,
            statements: BTreeMap::new(),
        });
    }
    // GlobalData is a DOT-codec concern. Rendering transports only topology;
    // typed drawing state travels separately in the V1 configuration.
    let spec = TypstGraphSpec {
        name: Some(state.name.clone().unwrap_or_else(|| "linnet".to_owned())),
        data: None,
        statements: BTreeMap::new(),
        default_edge_statements: BTreeMap::new(),
        default_node_statements: BTreeMap::new(),
        nodes,
        edges,
    };
    encode_graph_spec_bytes(&spec).map_err(PyRuntimeError::new_err)
}

pub(crate) fn render_graph(
    py: Python<'_>,
    graph: &Py<PyGraph>,
    config: Option<&Bound<'_, PyAny>>,
    subgraph: Option<&PySubgraph>,
) -> PyResult<crate::PyDiagramView> {
    // Validate before callbacks and again after preparation of the graph records.
    if let Some(subgraph) = subgraph {
        subgraph.selection_for(py, graph, graph.borrow(py).revision()?)?;
    }
    let (topology, transport) = request(py, graph, config, false)?;
    let selection = if let Some(subgraph) = subgraph {
        let owner = graph.borrow(py);
        let (selected, isolated) = subgraph.selection_for(py, graph, owner.revision()?)?;
        let state = owner.state.borrow();
        let state = state.as_ref().expect("checked selection owner");
        let hedges = (0..selected.size())
            .map(|index| selected.includes(&Hedge(index)))
            .collect();
        let nodes = state
            .graph
            .iter_nodes()
            .filter_map(|(node, mut crown, _)| {
                (isolated.contains(&node.0) || crown.any(|hedge| selected.includes(&hedge)))
                    .then_some(node.0)
            })
            .collect();
        Some((hedges, nodes))
    } else {
        None
    };
    let envelope: linnest::TypstGraphSpecEnvelope = ciborium::from_reader(topology.as_slice())
        .map_err(|error| PyRuntimeError::new_err(error.to_string()))?;
    let scene = transport.native_scene(envelope.graph, selection.as_ref(), None)?;
    crate::PyDiagramView::from_scene(py, scene, None)
}

/// Supply live generic graph views only when the optional Linnet bindings are registered.
pub(crate) fn scene_selectors(
    py: Python<'_>,
    spec: &TypstGraphSpec,
    elements: &Bound<'_, PyDict>,
    settings: &PyRenderSettings,
) -> PyResult<RenderSettingsTransport> {
    let node_drawings = elements
        .get_item("nodes")?
        .expect("scene nodes")
        .cast_into::<PyList>()?;
    let edge_drawings = elements
        .get_item("edges")?
        .expect("scene edges")
        .cast_into::<PyList>()?;
    let details = |drawing: &Bound<'_, PyDict>| -> PyResult<Py<PyAny>> {
        drawing
            .get_item("extensions")?
            .map(|extensions| extensions.cast_into::<PyDict>()?.get_item("inspection"))
            .transpose()
            .map(|value| {
                value
                    .flatten()
                    .map(Bound::unbind)
                    .unwrap_or_else(|| py.None())
            })
    };
    let mut builder = HedgeGraphBuilder::new();
    let mut nodes = BTreeMap::new();
    for (index, node) in spec.nodes.iter().enumerate() {
        let drawing = node_drawings.get_item(index)?.cast_into::<PyDict>()?;
        nodes.insert(
            node.index.unwrap_or(index),
            builder.add_node(NodeRecord {
                name: node.name.clone(),
                data: details(&drawing)?,
                drawing: drawing.unbind(),
            }),
        );
    }
    for (index, edge) in spec.edges.iter().enumerate() {
        let drawing = edge_drawings.get_item(index)?.cast_into::<PyDict>()?;
        let record = EdgeRecord {
            name: edge.name.clone(),
            data: details(&drawing)?,
            drawing: drawing.unbind(),
        };
        let endpoint = |end: &TypstEndpointSpec| HedgeData {
            node: nodes[&end.node],
            data: HalfEdgeRecord {
                data: py.None(),
                drawing: PyDict::new(py).unbind(),
            },
            is_in_subgraph: end.in_subgraph,
        };
        let orientation = match edge.orientation.as_deref() {
            Some("undirected") => Orientation::Undirected,
            Some("reversed") => Orientation::Reversed,
            _ => Orientation::Default,
        };
        match (&edge.source, &edge.sink) {
            (Some(source), Some(sink)) => {
                builder.add_edge(endpoint(source), endpoint(sink), record, orientation);
            }
            (Some(source), None) => {
                builder.add_external_edge(endpoint(source), record, orientation, Flow::Source);
            }
            (None, Some(sink)) => {
                builder.add_external_edge(endpoint(sink), record, orientation, Flow::Sink);
            }
            _ => return Err(PyRuntimeError::new_err("scene edge has no endpoints")),
        }
    }
    let graph = Py::new(
        py,
        PyGraph::from_state(GraphState {
            graph: PyHedgeGraph::build(builder, PyNodeStore::Vec),
            name: spec.name.clone(),
            global_data: PyGlobalData::default(),
            codec: None,
            render_config: Py::new(py, settings.clone())?.into_any(),
            revision: 0,
        }),
    )?;
    request(py, &graph, None, true).map(|(_, transport)| transport)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn measured_graph_inspection_preserves_native_drawing() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let modules = py.import("sys")?.getattr("modules")?;
            let symbolica = pyo3::types::PyModule::new(py, "symbolica")?;
            let community = pyo3::types::PyModule::new(py, "symbolica.community")?;
            symbolica.add("community", &community)?;
            modules.set_item("symbolica", &symbolica)?;
            modules.set_item("symbolica.community", &community)?;
            let linnet = pyo3::types::PyModule::new(py, "symbolica.community.graph")?;
            community.add("graph", &linnet)?;
            crate::register(&linnet)?;
            py.import("sys")?
                .getattr("modules")?
                .set_item("symbolica.community.graph", &linnet)?;
            let metrics_source =
                std::ffi::CString::new(include_str!("../tests/rendering_metrics.py")).unwrap();
            let metrics = PyModule::from_code(
                py,
                &metrics_source,
                c"rendering_metrics.py",
                c"rendering_metrics",
            )?;
            modules.set_item("rendering_metrics", &metrics)?;
            let source =
                std::ffi::CString::new(include_str!("../tests/test_svg_interaction.py")).unwrap();
            let fixture = PyModule::from_code(
                py,
                &source,
                c"test_svg_interaction.py",
                c"test_svg_interaction",
            )?;
            let case = fixture.getattr("SvgInteractionTests")?.call0()?;
            case.call_method0("setUp")?;
            let graph = case.getattr("graph")?.extract::<Py<PyGraph>>()?;
            let view = render_graph(py, &graph, None, None)?;
            let drawing = metrics
                .call_method1("external_metrics_snapshot", (Py::new(py, view)?,))?
                .extract::<PyRef<'_, crate::PyDiagramRender>>()?
                .clone();
            let pages = drawing.to_svg_pages()?;
            assert_eq!(pages.len(), 1);
            case.call_method1(
                "assert_native_drawing_unchanged",
                (&pages[0], drawing.to_svg()?),
            )?;
            Ok(())
        })
        .unwrap();
    }
}
