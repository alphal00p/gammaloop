use std::collections::{BTreeMap, BTreeSet};

use feynkit_graph::{ExternalState, FeynmanDiagram};
use serde::{Deserialize, Deserializer, Serialize, de};

use crate::surface::MAX_VERTEX_COUNT;
use crate::{CffError, EdgeId, VertexId, VertexSet};

/// Direction of an edge at the vertex to which it is attached.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EdgeFlow {
    Incoming,
    Outgoing,
}

/// Topological role of an edge in CFF generation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EdgeKind {
    /// An energy-carrying edge whose orientation is generated.
    Internal { source: VertexId, sink: VertexId },
    /// A physical external momentum entering or leaving one vertex.
    External { vertex: VertexId, flow: EdgeFlow },
    /// A severed internal edge at a subgraph boundary.
    ///
    /// Boundary edges retain an on-shell energy in the generated surface, but
    /// unlike internal edges their orientation is fixed by `flow`.
    Boundary { vertex: VertexId, flow: EdgeFlow },
    /// A sewn initial-state cut edge attached to both endpoints.
    ///
    /// It contributes external-energy shifts at both vertices but is excluded
    /// from cycle and orientation generation.
    InitialState { source: VertexId, sink: VertexId },
}

/// An edge with a stable caller-provided identifier.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CffEdge {
    pub id: EdgeId,
    pub kind: EdgeKind,
}

impl CffEdge {
    pub const fn internal(id: EdgeId, source: VertexId, sink: VertexId) -> Self {
        Self {
            id,
            kind: EdgeKind::Internal { source, sink },
        }
    }

    pub const fn external(id: EdgeId, vertex: VertexId, flow: EdgeFlow) -> Self {
        Self {
            id,
            kind: EdgeKind::External { vertex, flow },
        }
    }

    pub const fn boundary(id: EdgeId, vertex: VertexId, flow: EdgeFlow) -> Self {
        Self {
            id,
            kind: EdgeKind::Boundary { vertex, flow },
        }
    }

    pub const fn initial_state(id: EdgeId, source: VertexId, sink: VertexId) -> Self {
        Self {
            id,
            kind: EdgeKind::InitialState { source, sink },
        }
    }

    pub const fn is_internal(self) -> bool {
        matches!(self.kind, EdgeKind::Internal { .. })
    }
}

/// Minimal graph input needed by the CFF contraction algorithm.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct CffGraph {
    vertex_sets: Vec<VertexSet>,
    edges: Vec<CffEdge>,
}

#[derive(Deserialize)]
struct CffGraphSerde {
    vertex_sets: Vec<VertexSet>,
    edges: Vec<CffEdge>,
}

impl<'de> Deserialize<'de> for CffGraph {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let graph = CffGraphSerde::deserialize(deserializer)?;
        Self::with_vertex_sets(graph.vertex_sets, graph.edges).map_err(de::Error::custom)
    }
}

impl CffGraph {
    pub fn new(
        vertex_count: usize,
        edges: impl IntoIterator<Item = CffEdge>,
    ) -> Result<Self, CffError> {
        Self::with_vertex_sets(
            (0..vertex_count).map(|vertex| VertexSet::singleton(VertexId::new(vertex))),
            edges,
        )
    }

    /// Build a graph whose dense working vertices retain caller-provided identities.
    ///
    /// This lets subgraph generation intern surfaces directly into a shared arena:
    /// generated [`VertexSet`] metadata already refers to the containing graph and
    /// does not require a post-generation adapter.
    pub fn with_vertex_sets(
        vertex_sets: impl IntoIterator<Item = VertexSet>,
        edges: impl IntoIterator<Item = CffEdge>,
    ) -> Result<Self, CffError> {
        let graph = Self {
            vertex_sets: vertex_sets.into_iter().collect(),
            edges: edges.into_iter().collect(),
        };
        graph.validate()?;
        Ok(graph)
    }

    pub const fn vertex_count(&self) -> usize {
        self.vertex_sets.len()
    }

    pub fn vertex_sets(&self) -> &[VertexSet] {
        &self.vertex_sets
    }

    pub fn edges(&self) -> &[CffEdge] {
        &self.edges
    }

    pub fn edge(&self, id: EdgeId) -> Option<&CffEdge> {
        self.edges.iter().find(|edge| edge.id == id)
    }

    pub fn internal_edge_ids(&self) -> impl Iterator<Item = EdgeId> + '_ {
        self.edges
            .iter()
            .filter(|edge| edge.is_internal())
            .map(|edge| edge.id)
    }

    pub fn edge_count(&self) -> usize {
        self.edges
            .iter()
            .map(|edge| edge.id.index() + 1)
            .max()
            .unwrap_or(0)
    }

    fn validate(&self) -> Result<(), CffError> {
        if self.vertex_sets.is_empty() {
            return Err(CffError::EmptyGraph);
        }
        if self.vertex_sets.len() > MAX_VERTEX_COUNT {
            return Err(CffError::TooManyVertices {
                actual: self.vertex_sets.len(),
                maximum: MAX_VERTEX_COUNT,
            });
        }

        let mut edge_ids = BTreeSet::new();
        for edge in &self.edges {
            if !edge_ids.insert(edge.id) {
                return Err(CffError::DuplicateEdge(edge.id));
            }
            match edge.kind {
                EdgeKind::Internal { source, sink } => {
                    self.validate_vertex(edge.id, source)?;
                    self.validate_vertex(edge.id, sink)?;
                }
                EdgeKind::External { vertex, .. } | EdgeKind::Boundary { vertex, .. } => {
                    self.validate_vertex(edge.id, vertex)?;
                }
                EdgeKind::InitialState { source, sink } => {
                    self.validate_vertex(edge.id, source)?;
                    self.validate_vertex(edge.id, sink)?;
                }
            }
        }

        if self.vertex_sets.len() > 1 && !self.internal_topology_is_connected() {
            return Err(CffError::DisconnectedGraph);
        }
        Ok(())
    }

    fn validate_vertex(&self, edge: EdgeId, vertex: VertexId) -> Result<(), CffError> {
        if vertex.index() < self.vertex_sets.len() {
            Ok(())
        } else {
            Err(CffError::UnknownVertex {
                edge,
                vertex,
                vertex_count: self.vertex_sets.len(),
            })
        }
    }

    fn internal_topology_is_connected(&self) -> bool {
        let mut adjacency = BTreeMap::<VertexId, Vec<VertexId>>::new();
        for edge in &self.edges {
            if let EdgeKind::Internal { source, sink } = edge.kind {
                adjacency.entry(source).or_default().push(sink);
                adjacency.entry(sink).or_default().push(source);
            }
        }

        let seed = VertexId::new(0);
        let mut visited = BTreeSet::from([seed]);
        let mut pending = vec![seed];
        while let Some(vertex) = pending.pop() {
            for neighbour in adjacency.get(&vertex).into_iter().flatten() {
                if visited.insert(*neighbour) {
                    pending.push(*neighbour);
                }
            }
        }
        visited.len() == self.vertex_sets.len()
    }
}

/// Convert a complete Feynman diagram while preserving its edge identifiers.
///
/// Interaction vertices are densely remapped and external-leg metadata defines
/// incoming/outgoing flow. Subgraph boundary edges cannot be inferred from a
/// complete diagram and must still be supplied explicitly as [`CffEdge::boundary`].
impl TryFrom<&FeynmanDiagram> for CffGraph {
    type Error = CffError;

    fn try_from(diagram: &FeynmanDiagram) -> Result<Self, Self::Error> {
        diagram
            .validate()
            .map_err(|error| CffError::Invariant(error.to_string()))?;
        let vertex_ids = diagram
            .vertices()
            .enumerate()
            .map(|(dense, (id, _))| (id, VertexId::new(dense)))
            .collect::<BTreeMap<_, _>>();
        let mut edges = Vec::new();
        for (id, endpoints, data) in diagram.edges() {
            if data.is_dummy {
                continue;
            }
            let edge = EdgeId::new(id.0);
            let source = endpoints.source.and_then(|id| vertex_ids.get(&id)).copied();
            let sink = endpoints.target.and_then(|id| vertex_ids.get(&id)).copied();
            edges.push(match (source, sink, &data.external) {
                (Some(source), Some(sink), Some(_)) => CffEdge::initial_state(edge, source, sink),
                (Some(source), Some(sink), None) => CffEdge::internal(edge, source, sink),
                (Some(vertex), None, Some(external)) | (None, Some(vertex), Some(external)) => {
                    CffEdge::external(edge, vertex, EdgeFlow::from(external.state))
                }
                (Some(_), None, None) | (None, Some(_), None) => {
                    return Err(CffError::MissingExternalMetadata(edge));
                }
                (None, None, _) => return Err(CffError::ExternalToExternalEdge(edge)),
            });
        }
        Self::with_vertex_sets(
            vertex_ids
                .keys()
                .map(|id| VertexSet::singleton(VertexId::new(id.0))),
            edges,
        )
    }
}

impl From<ExternalState> for EdgeFlow {
    fn from(state: ExternalState) -> Self {
        match state {
            ExternalState::Incoming => Self::Incoming,
            ExternalState::Outgoing => Self::Outgoing,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use super::*;
    use feynkit_graph::{DiagramEdge, DiagramError, DiagramVertex, ExternalLeg};
    use feynkit_model::Model;

    #[test]
    fn rejects_duplicate_edge_ids() {
        let duplicate = EdgeId::new(4);
        let error = CffGraph::new(
            2,
            [
                CffEdge::internal(duplicate, VertexId::new(0), VertexId::new(1)),
                CffEdge::external(duplicate, VertexId::new(0), EdgeFlow::Incoming),
            ],
        )
        .unwrap_err();

        assert_eq!(error, CffError::DuplicateEdge(duplicate));
    }

    #[test]
    fn rejects_disconnected_internal_topology() {
        let error = CffGraph::new(
            3,
            [CffEdge::internal(
                EdgeId::new(0),
                VertexId::new(0),
                VertexId::new(1),
            )],
        )
        .unwrap_err();

        assert_eq!(error, CffError::DisconnectedGraph);
    }

    #[test]
    fn deserialization_revalidates_graph_invariants() {
        let graph = CffGraph::new(
            2,
            [CffEdge::internal(
                EdgeId::new(0),
                VertexId::new(0),
                VertexId::new(1),
            )],
        )
        .unwrap();
        let mut json = serde_json::to_value(graph).unwrap();
        json["vertex_sets"] = serde_json::json!([VertexSet::singleton(VertexId::new(0))]);

        assert!(serde_json::from_value::<CffGraph>(json).is_err());
    }

    fn scalar_model() -> Arc<Model> {
        Arc::new(
            Model::from_json(include_str!(
                "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
            ))
            .unwrap(),
        )
    }

    fn scalar_edge(model: &Model) -> DiagramEdge {
        DiagramEdge::new(model.particle_id("scalar_0").unwrap(), false)
    }

    fn external_scalar(model: &Model, index: usize, state: ExternalState) -> DiagramEdge {
        let mut edge = scalar_edge(model);
        edge.external = Some(ExternalLeg {
            name: format!("p{index}"),
            index,
            state,
            connection: index,
        });
        edge
    }

    fn bubble_diagram() -> FeynmanDiagram {
        let model = scalar_model();
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "bubble");
        let rule = model.vertex_rule_id("V_3_SCALAR_000").unwrap();
        let left = builder.add_vertex(DiagramVertex::interaction("left", rule));
        let right = builder.add_vertex(DiagramVertex::interaction("right", rule));
        builder
            .add_edge(
                None,
                left,
                external_scalar(&model, 0, ExternalState::Incoming),
            )
            .unwrap();
        builder
            .add_edge(
                right,
                None,
                external_scalar(&model, 1, ExternalState::Outgoing),
            )
            .unwrap();
        builder.add_edge(left, right, scalar_edge(&model)).unwrap();
        builder.add_edge(left, right, scalar_edge(&model)).unwrap();
        builder.build().unwrap()
    }

    #[test]
    fn converts_diagram_with_dense_internal_vertices_and_external_flow() {
        let graph = CffGraph::try_from(&bubble_diagram()).unwrap();

        assert_eq!(graph.vertex_count(), 2);
        assert_eq!(
            graph.edges(),
            &[
                CffEdge::external(EdgeId::new(0), VertexId::new(0), EdgeFlow::Incoming),
                CffEdge::external(EdgeId::new(1), VertexId::new(1), EdgeFlow::Outgoing),
                CffEdge::internal(EdgeId::new(2), VertexId::new(0), VertexId::new(1)),
                CffEdge::internal(EdgeId::new(3), VertexId::new(0), VertexId::new(1)),
            ]
        );
    }

    #[test]
    fn native_external_flow_and_metadata_agree() {
        let model = scalar_model();
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "external-flow");
        let internal = builder.add_vertex(DiagramVertex::interaction(
            "v",
            model.vertex_rule_id("V_3_SCALAR_000").unwrap(),
        ));
        builder
            .add_edge(
                None,
                internal,
                external_scalar(&model, 0, ExternalState::Incoming),
            )
            .unwrap();
        builder
            .add_edge(
                internal,
                None,
                external_scalar(&model, 1, ExternalState::Outgoing),
            )
            .unwrap();
        builder
            .add_edge(
                internal,
                None,
                external_scalar(&model, 2, ExternalState::Outgoing),
            )
            .unwrap();
        let diagram = builder.build().unwrap();

        let graph = CffGraph::try_from(&diagram).unwrap();
        assert_eq!(
            graph.edges().first(),
            Some(&CffEdge::external(
                EdgeId::new(0),
                VertexId::new(0),
                EdgeFlow::Incoming,
            ))
        );
    }

    #[test]
    fn diagram_builder_rejects_an_edge_without_any_attached_endpoint() {
        let model = scalar_model();
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "external-only");
        assert!(matches!(
            builder
                .add_edge(None, None, scalar_edge(&model))
                .unwrap_err(),
            DiagramError::Invariant {
                operation: "adding an edge",
                ..
            }
        ));
    }
}
