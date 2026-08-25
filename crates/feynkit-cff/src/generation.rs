use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

use crate::{
    CffError, CffExpression, CffGraph, EdgeFlow, EdgeId, EdgeKind, EdgeOrientation, EnergySurface,
    ExpressionTree, ExternalShift, GraphOrientation, HSurface, OrientationExpression,
    OrientationId, Surface, SurfaceCache, SurfaceId, VertexId, VertexSet,
};

/// Replace one dependent external momentum by a canonical linear combination.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ShiftRewrite {
    pub dependent_momentum: EdgeId,
    pub replacement: ExternalShift,
}

/// Explicit controls for topology-only CFF generation.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CffOptions {
    fixed_orientations: BTreeMap<EdgeId, EdgeOrientation>,
    contracted_edges: BTreeSet<EdgeId>,
    initial_state_edges: BTreeSet<EdgeId>,
    shift_rewrite: Option<ShiftRewrite>,
    max_orientations: usize,
}

impl Default for CffOptions {
    fn default() -> Self {
        Self {
            fixed_orientations: BTreeMap::new(),
            contracted_edges: BTreeSet::new(),
            initial_state_edges: BTreeSet::new(),
            shift_rewrite: None,
            max_orientations: usize::MAX,
        }
    }
}

impl CffOptions {
    pub fn with_fixed_orientation(mut self, edge: EdgeId, orientation: EdgeOrientation) -> Self {
        self.fixed_orientations.insert(edge, orientation);
        self
    }

    pub fn with_contracted_edge(mut self, edge: EdgeId) -> Self {
        self.contracted_edges.insert(edge);
        self
    }

    /// Treat an on-shell energy as an external shift after each surface is generated.
    pub fn with_initial_state_edge(mut self, edge: EdgeId) -> Self {
        self.initial_state_edges.insert(edge);
        self
    }

    pub fn with_shift_rewrite(mut self, rewrite: ShiftRewrite) -> Self {
        self.shift_rewrite = Some(rewrite);
        self
    }

    pub fn with_max_orientations(mut self, maximum: usize) -> Self {
        self.max_orientations = maximum;
        self
    }

    pub fn fixed_orientations(&self) -> &BTreeMap<EdgeId, EdgeOrientation> {
        &self.fixed_orientations
    }

    pub fn contracted_edges(&self) -> &BTreeSet<EdgeId> {
        &self.contracted_edges
    }

    pub fn initial_state_edges(&self) -> &BTreeSet<EdgeId> {
        &self.initial_state_edges
    }

    pub fn shift_rewrite(&self) -> Option<&ShiftRewrite> {
        self.shift_rewrite.as_ref()
    }

    pub const fn max_orientations(&self) -> usize {
        self.max_orientations
    }
}

/// Summary of the finite combinatorial search performed by a generator call.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct CffReport {
    pub candidate_orientations: usize,
    pub acyclic_orientations: usize,
    pub unfolded_terms: usize,
    pub interned_surfaces: usize,
}

/// Typed CFF expression together with its generation-local surfaces and report.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct CffResult {
    pub expression: CffExpression,
    pub surfaces: SurfaceCache,
    pub report: CffReport,
}

/// Reusable CFF generator configuration.
///
/// Surface caches are intentionally not stored here. Every call creates a new
/// cache, so identifiers in a result cannot depend on earlier generations.
#[derive(Clone, Debug, Default)]
pub struct CffGenerator {
    options: CffOptions,
}

impl CffGenerator {
    pub const fn new(options: CffOptions) -> Self {
        Self { options }
    }

    pub const fn options(&self) -> &CffOptions {
        &self.options
    }

    pub fn generate(&self, graph: &CffGraph) -> Result<CffResult, CffError> {
        self.validate_options(graph)?;

        let mut seed_orientation = graph
            .edges()
            .iter()
            .map(|edge| (edge.id, EdgeOrientation::Default))
            .collect::<BTreeMap<_, _>>();
        let mut seed =
            WorkingGraph::from_graph(graph, GraphOrientation::new(seed_orientation.clone()));
        for edge in &self.options.contracted_edges {
            seed.contract_edge(*edge)?;
            seed_orientation.insert(*edge, EdgeOrientation::Undirected);
        }

        let remaining_internal = seed.internal_edge_ids();
        let free_edges = remaining_internal
            .iter()
            .filter(|edge| !self.options.fixed_orientations.contains_key(edge))
            .copied()
            .collect::<Vec<_>>();
        let free_edge_count =
            u32::try_from(free_edges.len()).map_err(|_| CffError::OrientationCountOverflow)?;
        let candidate_count = 1_usize
            .checked_shl(free_edge_count)
            .ok_or(CffError::OrientationCountOverflow)?;
        if candidate_count > self.options.max_orientations {
            return Err(CffError::TooManyOrientations {
                candidate_count,
                maximum: self.options.max_orientations,
            });
        }

        let mut cache = SurfaceCache::default();
        let mut orientations = Vec::new();
        for identifier in 0..candidate_count {
            let mut orientation = seed_orientation.clone();
            for (edge, fixed) in &self.options.fixed_orientations {
                orientation.insert(*edge, *fixed);
            }
            for (bit, edge) in free_edges.iter().enumerate() {
                let direction = if identifier & (1_usize << bit) == 0 {
                    EdgeOrientation::Default
                } else {
                    EdgeOrientation::Reversed
                };
                orientation.insert(*edge, direction);
            }

            let orientation = GraphOrientation::new(orientation);
            let mut oriented_graph = seed.clone();
            oriented_graph.apply_orientation(orientation.clone());
            if oriented_graph.has_directed_cycle() {
                continue;
            }

            let id = OrientationId(orientations.len());
            let expression = self.generate_tree(oriented_graph, &mut cache, identifier)?;
            orientations.push(OrientationExpression {
                id,
                orientation,
                expression,
            });
        }

        let expression = CffExpression::new(orientations);
        let report = CffReport {
            candidate_orientations: candidate_count,
            acyclic_orientations: expression.orientations().len(),
            unfolded_terms: expression.unfolded_term_count(),
            interned_surfaces: cache.len(),
        };
        Ok(CffResult {
            expression,
            surfaces: cache,
            report,
        })
    }

    fn validate_options(&self, graph: &CffGraph) -> Result<(), CffError> {
        for edge in self
            .options
            .fixed_orientations
            .keys()
            .chain(&self.options.contracted_edges)
        {
            let graph_edge = graph.edge(*edge).ok_or(CffError::UnknownEdge(*edge))?;
            if !graph_edge.is_internal() {
                return Err(CffError::NotInternalEdge(*edge));
            }
        }
        for edge in &self.options.contracted_edges {
            if self.options.fixed_orientations.contains_key(edge) {
                return Err(CffError::OrientedContractedEdge(*edge));
            }
        }
        for (edge, orientation) in &self.options.fixed_orientations {
            if *orientation == EdgeOrientation::Undirected {
                return Err(CffError::UndirectedFixedEdge(*edge));
            }
        }
        for edge in &self.options.initial_state_edges {
            if graph.edge(*edge).is_none() {
                return Err(CffError::UnknownEdge(*edge));
            }
        }
        if let Some(rewrite) = &self.options.shift_rewrite {
            if graph.edge(rewrite.dependent_momentum).is_none() {
                return Err(CffError::UnknownEdge(rewrite.dependent_momentum));
            }
            for (edge, _) in rewrite.replacement.iter() {
                if graph.edge(edge).is_none() {
                    return Err(CffError::UnknownEdge(edge));
                }
            }
        }
        Ok(())
    }

    fn generate_tree(
        &self,
        graph: WorkingGraph,
        cache: &mut SurfaceCache,
        orientation: usize,
    ) -> Result<ExpressionTree<SurfaceId>, CffError> {
        let mut tree = ExpressionTree::from_root(GenerationNode {
            graph,
            surface: None,
        });

        loop {
            let leaves = tree.leaf_ids();
            let generated = leaves
                .iter()
                .map(|leaf| {
                    let graph = &tree
                        .node(*leaf)
                        .ok_or_else(|| {
                            CffError::Invariant(format!(
                                "leaf node {} is missing from its expression tree",
                                leaf.index()
                            ))
                        })?
                        .data
                        .graph;
                    let (children, surface) = graph.generate_children(orientation)?;
                    let surface = self.prepare_surface(surface)?;
                    Ok((children, cache.intern(surface)))
                })
                .collect::<Result<Vec<_>, CffError>>()?;

            for (leaf, (_, surface)) in leaves.iter().zip(&generated) {
                tree.data_mut(*leaf)
                    .ok_or_else(|| {
                        CffError::Invariant(format!(
                            "leaf node {} disappeared from its expression tree",
                            leaf.index()
                        ))
                    })?
                    .surface = Some(*surface);
            }

            let finished = generated.iter().all(|(children, _)| children.is_none());
            let continuing = generated.iter().all(|(children, _)| children.is_some());
            if !finished && !continuing {
                return Err(CffError::UnevenBranchDepth { orientation });
            }
            if finished {
                break;
            }

            for (leaf, (children, _)) in leaves.into_iter().zip(generated) {
                let children = children.ok_or_else(|| {
                    CffError::Invariant(format!(
                        "leaf node {} ended while sibling branches continued",
                        leaf.index()
                    ))
                })?;
                for child in children {
                    tree.insert(
                        leaf,
                        GenerationNode {
                            graph: child,
                            surface: None,
                        },
                    )
                    .ok_or_else(|| {
                        CffError::Invariant(format!(
                            "cannot insert a child below missing leaf node {}",
                            leaf.index()
                        ))
                    })?;
                }
            }
        }

        tree.try_map(|node| {
            node.surface.ok_or_else(|| {
                CffError::Invariant("completed generation node has no surface".to_owned())
            })
        })
    }

    fn prepare_surface(&self, surface: Surface) -> Result<Surface, CffError> {
        let surface = move_initial_state_energies(surface, &self.options.initial_state_edges)?;
        match (surface, &self.options.shift_rewrite) {
            (Surface::Energy(mut surface), Some(rewrite)) => {
                surface
                    .external_shift
                    .rewrite(rewrite.dependent_momentum, &rewrite.replacement)?;
                Ok(Surface::Energy(surface))
            }
            (surface, _) => Ok(surface),
        }
    }
}

#[derive(Clone, Debug)]
struct GenerationNode {
    graph: WorkingGraph,
    surface: Option<SurfaceId>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum WorkingEdgeKind {
    External,
    Internal,
    Boundary,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct WorkingEdge {
    id: EdgeId,
    kind: WorkingEdgeKind,
}

#[derive(Clone, Debug)]
struct WorkingVertex {
    vertices: VertexSet,
    incoming: Vec<WorkingEdge>,
    outgoing: Vec<WorkingEdge>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum VertexType {
    Source,
    Sink,
    Neither,
}

impl WorkingVertex {
    fn vertex_type(&self) -> VertexType {
        let is_sink = self
            .outgoing
            .iter()
            .all(|edge| edge.kind != WorkingEdgeKind::Internal);
        let is_source = self
            .incoming
            .iter()
            .all(|edge| edge.kind != WorkingEdgeKind::Internal);
        if is_sink {
            VertexType::Sink
        } else if is_source {
            VertexType::Source
        } else {
            VertexType::Neither
        }
    }

    fn generates_energy_surface(&self) -> bool {
        match self.vertex_type() {
            VertexType::Sink => self
                .outgoing
                .iter()
                .all(|edge| edge.kind == WorkingEdgeKind::External),
            VertexType::Source => self
                .incoming
                .iter()
                .all(|edge| edge.kind == WorkingEdgeKind::External),
            VertexType::Neither => false,
        }
    }

    fn edges(&self) -> impl Iterator<Item = WorkingEdge> + '_ {
        self.incoming.iter().chain(&self.outgoing).copied()
    }
}

#[derive(Clone, Debug)]
struct WorkingGraph {
    vertices: Vec<WorkingVertex>,
    orientation: GraphOrientation,
}

impl WorkingGraph {
    fn from_graph(graph: &CffGraph, orientation: GraphOrientation) -> Self {
        let mut vertices = (0..graph.vertex_count())
            .map(|vertex| WorkingVertex {
                vertices: VertexSet::singleton(VertexId::new(vertex)),
                incoming: Vec::new(),
                outgoing: Vec::new(),
            })
            .collect::<Vec<_>>();

        for edge in graph.edges() {
            match edge.kind {
                EdgeKind::Internal { source, sink } => {
                    let working = WorkingEdge {
                        id: edge.id,
                        kind: WorkingEdgeKind::Internal,
                    };
                    vertices[source.index()].outgoing.push(working);
                    vertices[sink.index()].incoming.push(working);
                }
                EdgeKind::External { vertex, flow } => {
                    Self::attach(
                        &mut vertices[vertex.index()],
                        WorkingEdge {
                            id: edge.id,
                            kind: WorkingEdgeKind::External,
                        },
                        flow,
                    );
                }
                EdgeKind::Boundary { vertex, flow } => {
                    Self::attach(
                        &mut vertices[vertex.index()],
                        WorkingEdge {
                            id: edge.id,
                            kind: WorkingEdgeKind::Boundary,
                        },
                        flow,
                    );
                }
            }
        }
        Self {
            vertices,
            orientation,
        }
    }

    fn attach(vertex: &mut WorkingVertex, edge: WorkingEdge, flow: EdgeFlow) {
        match flow {
            EdgeFlow::Incoming => vertex.incoming.push(edge),
            EdgeFlow::Outgoing => vertex.outgoing.push(edge),
        }
    }

    fn apply_orientation(&mut self, orientation: GraphOrientation) {
        for (edge, direction) in orientation.iter() {
            if direction != EdgeOrientation::Reversed {
                continue;
            }
            for vertex in &mut self.vertices {
                if let Some(position) = vertex.outgoing.iter().position(|item| item.id == edge) {
                    let edge = vertex.outgoing.remove(position);
                    vertex.incoming.push(edge);
                } else if let Some(position) =
                    vertex.incoming.iter().position(|item| item.id == edge)
                {
                    let edge = vertex.incoming.remove(position);
                    vertex.outgoing.push(edge);
                }
            }
        }
        self.orientation = orientation;
    }

    fn internal_edge_ids(&self) -> Vec<EdgeId> {
        self.vertices
            .iter()
            .flat_map(WorkingVertex::edges)
            .filter(|edge| edge.kind == WorkingEdgeKind::Internal)
            .map(|edge| edge.id)
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect()
    }

    fn vertex(&self, vertices: VertexSet) -> Option<&WorkingVertex> {
        self.vertices
            .iter()
            .find(|vertex| vertex.vertices == vertices)
    }

    fn directed_neighbours(&self, vertices: VertexSet) -> Vec<VertexSet> {
        let Some(vertex) = self.vertex(vertices) else {
            return Vec::new();
        };
        vertex
            .outgoing
            .iter()
            .filter(|edge| edge.kind == WorkingEdgeKind::Internal)
            .filter_map(|edge| {
                self.vertices
                    .iter()
                    .find(|candidate| candidate.incoming.contains(edge))
                    .map(|candidate| candidate.vertices)
            })
            .collect()
    }

    fn undirected_neighbours(&self, vertices: VertexSet) -> Vec<VertexSet> {
        self.vertices
            .iter()
            .filter(|other| other.vertices != vertices)
            .filter(|other| {
                self.are_directed_adjacent(vertices, other.vertices)
                    || self.are_directed_adjacent(other.vertices, vertices)
            })
            .map(|other| other.vertices)
            .collect()
    }

    fn are_directed_adjacent(&self, left: VertexSet, right: VertexSet) -> bool {
        let (Some(left), Some(right)) = (self.vertex(left), self.vertex(right)) else {
            return false;
        };
        left.outgoing
            .iter()
            .filter(|edge| edge.kind == WorkingEdgeKind::Internal)
            .any(|edge| right.incoming.contains(edge))
    }

    fn has_directed_cycle(&self) -> bool {
        self.vertices.iter().any(|vertex| {
            self.depth_first_cycle(vertex.vertices, &mut BTreeSet::new(), &mut Vec::new())
        })
    }

    fn depth_first_cycle(
        &self,
        vertex: VertexSet,
        visited: &mut BTreeSet<VertexSet>,
        stack: &mut Vec<VertexSet>,
    ) -> bool {
        if visited.contains(&vertex) {
            return stack.contains(&vertex);
        }
        visited.insert(vertex);
        stack.push(vertex);
        if self
            .directed_neighbours(vertex)
            .into_iter()
            .any(|neighbour| self.depth_first_cycle(neighbour, visited, stack))
        {
            return true;
        }
        stack.pop();
        false
    }

    fn has_connected_complement(&self, removed: VertexSet) -> bool {
        if self.vertices.len() == 1 {
            return true;
        }
        let Some(seed) = self
            .vertices
            .iter()
            .find(|vertex| vertex.vertices != removed)
            .map(|vertex| vertex.vertices)
        else {
            return true;
        };
        let mut visited = BTreeSet::from([seed]);
        let mut pending = vec![seed];
        while let Some(vertex) = pending.pop() {
            for neighbour in self.undirected_neighbours(vertex) {
                if neighbour != removed && visited.insert(neighbour) {
                    pending.push(neighbour);
                }
            }
        }
        visited.len() == self.vertices.len() - 1
    }

    fn valid_source_or_sink(&self, vertex: &WorkingVertex) -> bool {
        vertex.vertex_type() != VertexType::Neither
            && self.has_connected_complement(vertex.vertices)
    }

    fn contract_edge(&mut self, edge: EdgeId) -> Result<(), CffError> {
        let source = self
            .vertices
            .iter()
            .find(|vertex| vertex.outgoing.iter().any(|item| item.id == edge))
            .map(|vertex| vertex.vertices)
            .ok_or(CffError::NotInternalEdge(edge))?;
        let sink = self
            .vertices
            .iter()
            .find(|vertex| vertex.incoming.iter().any(|item| item.id == edge))
            .map(|vertex| vertex.vertices)
            .ok_or(CffError::NotInternalEdge(edge))?;
        if source == sink {
            self.remove_edge(edge);
        } else {
            *self = self.contract_vertices(source, sink, Some(edge))?;
        }
        Ok(())
    }

    fn remove_edge(&mut self, edge: EdgeId) {
        for vertex in &mut self.vertices {
            vertex.incoming.retain(|item| item.id != edge);
            vertex.outgoing.retain(|item| item.id != edge);
        }
    }

    fn contract_vertices(
        &self,
        left: VertexSet,
        right: VertexSet,
        remove_single: Option<EdgeId>,
    ) -> Result<Self, CffError> {
        let left_vertex = self
            .vertex(left)
            .ok_or_else(|| CffError::Invariant("left contraction vertex is missing".to_owned()))?;
        let right_vertex = self
            .vertex(right)
            .ok_or_else(|| CffError::Invariant("right contraction vertex is missing".to_owned()))?;
        let shared_edges = left_vertex
            .edges()
            .filter(|left_edge| {
                right_vertex
                    .edges()
                    .any(|right_edge| left_edge == &right_edge)
            })
            .map(|edge| edge.id)
            .collect::<BTreeSet<_>>();
        let retain = |edge: &WorkingEdge| match remove_single {
            Some(remove) => edge.id != remove,
            None => !shared_edges.contains(&edge.id),
        };
        let mut incoming = left_vertex
            .incoming
            .iter()
            .chain(&right_vertex.incoming)
            .filter(|edge| retain(edge))
            .copied()
            .collect::<Vec<_>>();
        let mut outgoing = left_vertex
            .outgoing
            .iter()
            .chain(&right_vertex.outgoing)
            .filter(|edge| retain(edge))
            .copied()
            .collect::<Vec<_>>();
        incoming.sort_by_key(|edge| edge.id);
        outgoing.sort_by_key(|edge| edge.id);

        let mut vertices = self
            .vertices
            .iter()
            .filter(|vertex| vertex.vertices != left && vertex.vertices != right)
            .cloned()
            .collect::<Vec<_>>();
        vertices.push(WorkingVertex {
            vertices: left.union(right),
            incoming,
            outgoing,
        });
        vertices.sort_by_key(|vertex| vertex.vertices);
        Ok(Self {
            vertices,
            orientation: self.orientation.clone(),
        })
    }

    fn generate_children(
        &self,
        orientation: usize,
    ) -> Result<(Option<Vec<Self>>, Surface), CffError> {
        if self.vertices.len() < 2 {
            return Ok((None, Surface::Unit));
        }
        let vertex = self
            .vertices
            .iter()
            .find(|vertex| self.valid_source_or_sink(vertex))
            .ok_or(CffError::NoSourceOrSink { orientation })?;
        let surface = self.surface_for(vertex)?;

        if self.vertices.len() == 2 {
            return Ok((None, surface));
        }
        let children = self
            .undirected_neighbours(vertex.vertices)
            .into_iter()
            .map(|adjacent| self.contract_vertices(vertex.vertices, adjacent, None))
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .filter(|graph| !graph.has_directed_cycle())
            .collect::<Vec<_>>();
        if children.is_empty() {
            Err(CffError::NoContractionBranch { orientation })
        } else {
            Ok((Some(children), surface))
        }
    }

    fn surface_for(&self, vertex: &WorkingVertex) -> Result<Surface, CffError> {
        let vertex_type = vertex.vertex_type();
        if vertex_type == VertexType::Neither {
            return Err(CffError::Invariant(
                "only a source or sink may generate a surface".to_owned(),
            ));
        }
        let mut external_shift = ExternalShift::default();
        for edge in vertex
            .incoming
            .iter()
            .filter(|edge| edge.kind == WorkingEdgeKind::External)
        {
            external_shift.add(
                edge.id,
                match vertex_type {
                    VertexType::Source => -1,
                    VertexType::Sink => 1,
                    VertexType::Neither => 0,
                },
            )?;
        }
        for edge in vertex
            .outgoing
            .iter()
            .filter(|edge| edge.kind == WorkingEdgeKind::External)
        {
            external_shift.add(
                edge.id,
                match vertex_type {
                    VertexType::Source => 1,
                    VertexType::Sink => -1,
                    VertexType::Neither => 0,
                },
            )?;
        }

        let mut positive_energies = vertex
            .edges()
            .filter(|edge| edge.kind == WorkingEdgeKind::Internal)
            .map(|edge| edge.id)
            .collect::<Vec<_>>();
        if vertex.generates_energy_surface() {
            positive_energies.extend(
                vertex
                    .edges()
                    .filter(|edge| edge.kind == WorkingEdgeKind::Boundary)
                    .map(|edge| edge.id),
            );
            positive_energies.sort();
            return Ok(Surface::Energy(EnergySurface {
                energies: positive_energies,
                external_shift,
                vertices: vertex.vertices,
            }));
        }

        let (positive_boundary, negative_energies): (Vec<_>, Vec<_>) = match vertex_type {
            VertexType::Sink => (
                vertex
                    .incoming
                    .iter()
                    .filter(|edge| edge.kind == WorkingEdgeKind::Boundary)
                    .map(|edge| edge.id)
                    .collect(),
                vertex
                    .outgoing
                    .iter()
                    .filter(|edge| edge.kind == WorkingEdgeKind::Boundary)
                    .map(|edge| edge.id)
                    .collect(),
            ),
            VertexType::Source => (
                vertex
                    .outgoing
                    .iter()
                    .filter(|edge| edge.kind == WorkingEdgeKind::Boundary)
                    .map(|edge| edge.id)
                    .collect(),
                vertex
                    .incoming
                    .iter()
                    .filter(|edge| edge.kind == WorkingEdgeKind::Boundary)
                    .map(|edge| edge.id)
                    .collect(),
            ),
            VertexType::Neither => {
                return Err(CffError::Invariant(
                    "only a source or sink may generate a surface".to_owned(),
                ));
            }
        };
        positive_energies.extend(positive_boundary);
        positive_energies.sort();
        let mut negative_energies = negative_energies;
        negative_energies.sort();
        Ok(Surface::H(HSurface {
            positive_energies,
            negative_energies,
            external_shift,
            vertices: vertex.vertices,
        }))
    }
}

fn move_initial_state_energies(
    surface: Surface,
    initial: &BTreeSet<EdgeId>,
) -> Result<Surface, CffError> {
    match surface {
        Surface::Energy(mut surface) => {
            let moved = surface
                .energies
                .iter()
                .filter(|edge| initial.contains(edge))
                .copied()
                .collect::<Vec<_>>();
            surface.energies.retain(|edge| !initial.contains(edge));
            for edge in moved {
                surface.external_shift.add(edge, 1)?;
            }
            Ok(Surface::Energy(surface))
        }
        Surface::H(mut surface) => {
            let moved_positive = surface
                .positive_energies
                .iter()
                .filter(|edge| initial.contains(edge))
                .copied()
                .collect::<Vec<_>>();
            let moved_negative = surface
                .negative_energies
                .iter()
                .filter(|edge| initial.contains(edge))
                .copied()
                .collect::<Vec<_>>();
            surface
                .positive_energies
                .retain(|edge| !initial.contains(edge));
            surface
                .negative_energies
                .retain(|edge| !initial.contains(edge));
            for edge in moved_positive {
                surface.external_shift.add(edge, 1)?;
            }
            for edge in moved_negative {
                surface.external_shift.add(edge, -1)?;
            }
            if surface.negative_energies.is_empty() {
                Ok(Surface::Energy(EnergySurface {
                    energies: surface.positive_energies,
                    external_shift: surface.external_shift,
                    vertices: surface.vertices,
                }))
            } else {
                Ok(Surface::H(surface))
            }
        }
        surface => Ok(surface),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::CffEdge;
    use feynkit_graph::{
        DiagramEdge, DiagramVertex, ExternalState, FeynmanDiagram, ParticleReference,
    };

    fn edge(index: usize) -> EdgeId {
        EdgeId::new(index)
    }

    fn vertex(index: usize) -> VertexId {
        VertexId::new(index)
    }

    fn bubble() -> CffGraph {
        CffGraph::new(
            2,
            [
                CffEdge::external(edge(0), vertex(0), EdgeFlow::Incoming),
                CffEdge::external(edge(1), vertex(1), EdgeFlow::Outgoing),
                CffEdge::internal(edge(2), vertex(0), vertex(1)),
                CffEdge::internal(edge(3), vertex(0), vertex(1)),
            ],
        )
        .unwrap()
    }

    #[test]
    fn generates_only_acyclic_bubble_orientations() {
        let result = CffGenerator::default().generate(&bubble()).unwrap();

        assert_eq!(result.report.candidate_orientations, 4);
        assert_eq!(result.report.acyclic_orientations, 2);
        assert_eq!(result.report.unfolded_terms, 2);
        assert_eq!(result.surfaces.energy_surfaces().len(), 2);
        assert!(result.surfaces.h_surfaces().is_empty());
        assert!(
            result
                .expression
                .orientations()
                .iter()
                .all(|orientation| orientation.denominator_products().len() == 1)
        );
    }

    #[test]
    fn triangle_preserves_the_source_sink_contraction_tree() {
        let graph = CffGraph::new(
            3,
            [
                CffEdge::internal(edge(0), vertex(0), vertex(1)),
                CffEdge::internal(edge(1), vertex(1), vertex(2)),
                CffEdge::internal(edge(2), vertex(0), vertex(2)),
            ],
        )
        .unwrap();
        let result = CffGenerator::default().generate(&graph).unwrap();

        assert_eq!(result.report.candidate_orientations, 8);
        assert_eq!(result.report.acyclic_orientations, 6);
        // Of the two apparent contractions from a triangle source/sink, one
        // would orient the remaining parallel edges into a directed cycle.
        assert_eq!(result.report.unfolded_terms, 6);
        assert!(
            result
                .expression
                .orientations()
                .iter()
                .all(|orientation| orientation.denominator_products().len() == 1)
        );
    }

    #[test]
    fn fixed_orientation_prunes_the_search_before_generation() {
        let options = CffOptions::default()
            .with_fixed_orientation(edge(2), EdgeOrientation::Default)
            .with_fixed_orientation(edge(3), EdgeOrientation::Default);
        let result = CffGenerator::new(options).generate(&bubble()).unwrap();

        assert_eq!(result.report.candidate_orientations, 1);
        assert_eq!(result.report.acyclic_orientations, 1);
    }

    #[test]
    fn cache_is_local_to_each_generation() {
        let generator = CffGenerator::default();
        let first = generator.generate(&bubble()).unwrap();
        let second = generator.generate(&bubble()).unwrap();

        assert_eq!(first.surfaces, second.surfaces);
        assert_eq!(first.expression, second.expression);
    }

    #[test]
    fn boundary_edges_generate_h_surfaces() {
        let graph = CffGraph::new(
            2,
            [
                CffEdge::boundary(edge(0), vertex(0), EdgeFlow::Incoming),
                CffEdge::external(edge(1), vertex(1), EdgeFlow::Outgoing),
                CffEdge::internal(edge(2), vertex(0), vertex(1)),
            ],
        )
        .unwrap();
        let result = CffGenerator::default().generate(&graph).unwrap();

        assert!(!result.surfaces.h_surfaces().is_empty());
        assert!(result.surfaces.h_surfaces().iter().any(|surface| {
            surface.negative_energies == vec![edge(0)]
                || surface.positive_energies.contains(&edge(0))
        }));
    }

    #[test]
    fn initial_state_energy_and_dependent_shift_are_canonicalized() {
        let graph = CffGraph::new(
            2,
            [
                CffEdge::external(edge(0), vertex(0), EdgeFlow::Incoming),
                CffEdge::external(edge(1), vertex(1), EdgeFlow::Outgoing),
                CffEdge::internal(edge(2), vertex(0), vertex(1)),
            ],
        )
        .unwrap();
        let options = CffOptions::default()
            .with_initial_state_edge(edge(2))
            .with_shift_rewrite(ShiftRewrite {
                dependent_momentum: edge(2),
                replacement: ExternalShift::new([(edge(0), -1), (edge(1), 1)]).unwrap(),
            });
        let result = CffGenerator::new(options).generate(&graph).unwrap();

        assert!(result.surfaces.energy_surfaces().iter().all(|surface| {
            !surface.energies.contains(&edge(2)) && surface.external_shift.get(edge(2)) == 0
        }));
    }

    #[test]
    fn reports_orientation_budget_before_searching() {
        let error = CffGenerator::new(CffOptions::default().with_max_orientations(3))
            .generate(&bubble())
            .unwrap_err();

        assert_eq!(
            error,
            CffError::TooManyOrientations {
                candidate_count: 4,
                maximum: 3,
            }
        );
    }

    #[test]
    fn rejects_undirected_fixed_orientation() {
        let error = CffGenerator::new(
            CffOptions::default().with_fixed_orientation(edge(2), EdgeOrientation::Undirected),
        )
        .generate(&bubble())
        .unwrap_err();

        assert_eq!(error, CffError::UndirectedFixedEdge(edge(2)));
    }

    #[test]
    fn result_round_trips_through_json() {
        let result = CffGenerator::default().generate(&bubble()).unwrap();
        let json = serde_json::to_string(&result).unwrap();
        let decoded: CffResult = serde_json::from_str(&json).unwrap();

        assert_eq!(decoded, result);
    }

    #[test]
    fn generates_from_a_feynman_diagram_conversion() {
        let mut builder = FeynmanDiagram::builder("bubble");
        let incoming =
            builder.add_vertex(DiagramVertex::external("p1", 0, ExternalState::Incoming));
        let outgoing =
            builder.add_vertex(DiagramVertex::external("p2", 1, ExternalState::Outgoing));
        let left = builder.add_vertex(DiagramVertex::interaction("left", "V_1"));
        let right = builder.add_vertex(DiagramVertex::interaction("right", "V_1"));
        let scalar = || DiagramEdge::new(ParticleReference::new("phi", 25), false);
        builder.add_edge(incoming, left, scalar()).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(right, outgoing, scalar()).unwrap();
        let diagram = builder.build().unwrap();

        let graph = CffGraph::try_from(&diagram).unwrap();
        let result = CffGenerator::default().generate(&graph).unwrap();

        assert_eq!(result.report.candidate_orientations, 4);
        assert_eq!(result.report.acyclic_orientations, 2);
        assert_eq!(result.report.unfolded_terms, 2);
    }
}
