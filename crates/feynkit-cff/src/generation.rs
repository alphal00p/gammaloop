use std::collections::{BTreeMap, BTreeSet};

use feynkit_graph::FeynmanDiagram;
use linnet::half_edge::{
    HedgeGraph, NodeIndex,
    involution::{EdgeIndex, Flow, HedgePair, Orientation},
    nodestore::NodeStorageOps,
    subgraph::SubGraphLike,
};
use serde::{Deserialize, Serialize};

use crate::{
    CffEdge, CffError, CffExpression, CffGraph, EdgeFlow, EdgeId, EdgeKind, EdgeOrientation,
    EnergySurface, ExpressionTree, ExternalShift, GraphOrientation, HSurface, OrientationData,
    OrientationExpression, OrientationId, Surface, SurfaceCache, SurfaceId, VertexId, VertexSet,
};

/// Replace one dependent external momentum by a canonical linear combination.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, bincode::Encode, bincode::Decode)]
pub struct ShiftRewrite {
    pub dependent_momentum: EdgeId,
    pub replacement: ExternalShift,
}

/// Explicit controls for topology-only CFF generation.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, bincode::Encode, bincode::Decode)]
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
#[derive(
    Clone,
    Copy,
    Debug,
    Default,
    PartialEq,
    Eq,
    Serialize,
    Deserialize,
    bincode::Encode,
    bincode::Decode,
)]
pub struct CffReport {
    pub candidate_orientations: usize,
    pub acyclic_orientations: usize,
    pub unfolded_terms: usize,
    pub interned_surfaces: usize,
}

/// Typed CFF expression together with its generation-local surfaces and report.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, bincode::Encode, bincode::Decode)]
pub struct CffResult {
    pub expression: CffExpression,
    pub surfaces: SurfaceCache,
    pub report: CffReport,
}

impl std::ops::Deref for CffResult {
    type Target = CffExpression;

    fn deref(&self) -> &Self::Target {
        &self.expression
    }
}

impl std::ops::DerefMut for CffResult {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.expression
    }
}

impl CffResult {
    /// Select every non-zero residue order while retaining this result's arena.
    pub fn select_energy_surface_residue(
        &self,
        group: &crate::RaisedEnergySurfaceGroup,
    ) -> Vec<Self> {
        self.expression
            .select_energy_surface_residue(group)
            .into_iter()
            .map(|expression| Self {
                expression,
                surfaces: self.surfaces.clone(),
                report: self.report,
            })
            .collect()
    }

    pub fn normalize_raised_surfaces(&mut self, raised: &crate::RaisedEnergySurfaceData) {
        self.expression.normalize_raised_surfaces(raised);
    }
}

/// A generated expression whose IDs refer to a caller-owned [`SurfaceCache`].
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, bincode::Encode, bincode::Decode)]
pub struct CffGeneration {
    pub expression: CffExpression,
    pub report: CffReport,
}

/// Reusable CFF generator configuration.
///
/// Surface caches are intentionally not stored here. [`Self::generate`] creates
/// a fresh arena, while [`Self::generate_into`] makes reuse of a caller-owned
/// arena explicit.
#[derive(Clone, Debug, Default)]
pub struct CffGenerator {
    options: CffOptions,
}

/// CFF construction directly on FeynKit diagrams without wrapping the graph type.
pub trait FeynmanDiagramCffExt {
    fn build_cff(&self, options: CffOptions) -> Result<CffResult, CffError>;
}

impl FeynmanDiagramCffExt for FeynmanDiagram {
    fn build_cff(&self, options: CffOptions) -> Result<CffResult, CffError> {
        let graph = CffGraph::try_from(self)?;
        CffGenerator::new(options).generate(&graph)
    }
}

/// CFF construction directly on a selected Linnet topology.
///
/// The selected vertices retain their original Linnet identities in every
/// generated surface, and orientations remain indexed by the original edge
/// IDs. `reversed_dangling` changes the flow of boundary edges produced by the
/// selection. Initial-state edges and shift rewrites are supplied through
/// [`CffOptions`]. The caller owns `cache`; the returned result contains a
/// snapshot of that arena after generation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HedgeEdgeRole {
    /// Include the edge with the role implied by its Linnet topology.
    #[default]
    Standard,
    /// Leave the edge out of CFF topology and return it as undirected.
    Omitted,
    /// Include an unpaired external edge in CFF surfaces but leave its output
    /// orientation undirected.
    UnorientedExternal,
    /// Treat a paired edge as a sewn initial-state connection.
    InitialState,
}

pub trait HedgeGraphCffExt {
    fn build_cff_from_subgraph<S: SubGraphLike>(
        &self,
        subgraph: &S,
        options: CffOptions,
        reversed_dangling: &[EdgeId],
        cache: &mut SurfaceCache,
    ) -> Result<CffResult, CffError>;

    /// Build CFF while assigning application-specific roles to selected edges.
    ///
    /// The callback classifies original Linnet edge IDs without constructing an
    /// intermediate graph wrapper. This is useful for finalized runtime graphs
    /// that retain bookkeeping edges which must not participate in CFF, or sewn
    /// initial-state edges which contribute shifts but no generated orientation.
    fn build_cff_from_subgraph_with_edge_roles<S, F>(
        &self,
        subgraph: &S,
        options: CffOptions,
        reversed_dangling: &[EdgeId],
        edge_role: F,
        cache: &mut SurfaceCache,
    ) -> Result<CffResult, CffError>
    where
        S: SubGraphLike,
        F: FnMut(EdgeId) -> HedgeEdgeRole;
}

impl<E, V, H, N> HedgeGraphCffExt for HedgeGraph<E, V, H, N>
where
    N: NodeStorageOps<NodeData = V>,
{
    fn build_cff_from_subgraph<S: SubGraphLike>(
        &self,
        subgraph: &S,
        options: CffOptions,
        reversed_dangling: &[EdgeId],
        cache: &mut SurfaceCache,
    ) -> Result<CffResult, CffError> {
        self.build_cff_from_subgraph_with_edge_roles(
            subgraph,
            options,
            reversed_dangling,
            |_| HedgeEdgeRole::Standard,
            cache,
        )
    }

    fn build_cff_from_subgraph_with_edge_roles<S, F>(
        &self,
        subgraph: &S,
        options: CffOptions,
        reversed_dangling: &[EdgeId],
        edge_role: F,
        cache: &mut SurfaceCache,
    ) -> Result<CffResult, CffError>
    where
        S: SubGraphLike,
        F: FnMut(EdgeId) -> HedgeEdgeRole,
    {
        let input = HedgeSubgraphCffInput::new(self, subgraph, reversed_dangling, edge_role)?;
        let contracted_edges = options.contracted_edges().clone();
        let generation = CffGenerator::new(options).generate_into(&input.graph, cache)?;
        let mut expression = generation.expression;

        for orientation in expression.orientations_mut() {
            for edge in &input.internal_edges {
                if contracted_edges.contains(edge) {
                    continue;
                }
                match orientation.data.orientation.get(*edge) {
                    Some(EdgeOrientation::Default | EdgeOrientation::Reversed) => {}
                    Some(EdgeOrientation::Undirected) | None => {
                        return Err(CffError::Invariant(format!(
                            "generated orientation does not direct represented edge {edge}"
                        )));
                    }
                }
            }

            orientation.data.orientation = self.new_edgevec(|_, edge, _| {
                if !input.represented_edges.contains(&edge) {
                    return Orientation::Undirected;
                }
                if contracted_edges.contains(&edge)
                    || input.unoriented_external_edges.contains(&edge)
                {
                    return Orientation::Undirected;
                }
                if input.reversed_boundary_edges.contains(&edge) {
                    return Orientation::Reversed;
                }
                match orientation.data.orientation.get(edge) {
                    Some(EdgeOrientation::Reversed) if input.internal_edges.contains(&edge) => {
                        Orientation::Reversed
                    }
                    _ => Orientation::Default,
                }
            });
        }

        Ok(CffResult {
            expression,
            surfaces: cache.clone(),
            report: generation.report,
        })
    }
}

struct HedgeSubgraphCffInput {
    graph: CffGraph,
    represented_edges: BTreeSet<EdgeIndex>,
    internal_edges: BTreeSet<EdgeIndex>,
    reversed_boundary_edges: BTreeSet<EdgeIndex>,
    unoriented_external_edges: BTreeSet<EdgeIndex>,
}

impl HedgeSubgraphCffInput {
    fn new<E, V, H, N, S, F>(
        graph: &HedgeGraph<E, V, H, N>,
        subgraph: &S,
        reversed_dangling: &[EdgeId],
        mut edge_role: F,
    ) -> Result<Self, CffError>
    where
        N: NodeStorageOps<NodeData = V>,
        S: SubGraphLike,
        F: FnMut(EdgeId) -> HedgeEdgeRole,
    {
        let mut original_vertices = graph
            .iter_nodes_of(subgraph)
            .map(|(vertex, _, _)| vertex)
            .collect::<Vec<_>>();
        original_vertices.sort_unstable();
        if let Some(vertex) = original_vertices
            .iter()
            .find(|vertex| vertex.0 >= u64::BITS as usize)
        {
            return Err(CffError::VertexIdentityOutOfRange {
                vertex: *vertex,
                maximum: u64::BITS as usize - 1,
            });
        }

        let dense_vertices = original_vertices
            .iter()
            .enumerate()
            .map(|(dense, original)| (*original, VertexId::new(dense)))
            .collect::<BTreeMap<_, _>>();
        let dense_vertex = |vertex: NodeIndex, edge: EdgeIndex| {
            dense_vertices
                .get(&vertex)
                .copied()
                .ok_or(CffError::SubgraphEndpointOutsideSelection { edge, vertex })
        };

        let reversed_dangling = reversed_dangling.iter().copied().collect::<BTreeSet<_>>();
        let mut edges = Vec::new();
        let mut represented_edges = BTreeSet::new();
        let mut internal_edges = BTreeSet::new();
        let mut reversed_boundary_edges = BTreeSet::new();
        let mut unoriented_external_edges = BTreeSet::new();
        for (pair, edge, _) in graph.iter_edges_of(subgraph) {
            let role = edge_role(edge);
            if role == HedgeEdgeRole::Omitted {
                continue;
            }
            represented_edges.insert(edge);
            let cff_edge = match (role, pair) {
                (HedgeEdgeRole::Standard, HedgePair::Unpaired { hedge, flow }) => {
                    CffEdge::external(
                        edge,
                        dense_vertex(graph.node_id(hedge), edge)?,
                        match flow {
                            Flow::Sink => EdgeFlow::Incoming,
                            Flow::Source => EdgeFlow::Outgoing,
                        },
                    )
                }
                (HedgeEdgeRole::UnorientedExternal, HedgePair::Unpaired { hedge, flow }) => {
                    unoriented_external_edges.insert(edge);
                    CffEdge::external(
                        edge,
                        dense_vertex(graph.node_id(hedge), edge)?,
                        match flow {
                            Flow::Sink => EdgeFlow::Incoming,
                            Flow::Source => EdgeFlow::Outgoing,
                        },
                    )
                }
                (HedgeEdgeRole::Standard, HedgePair::Paired { source, sink }) => {
                    internal_edges.insert(edge);
                    CffEdge::internal(
                        edge,
                        dense_vertex(graph.node_id(source), edge)?,
                        dense_vertex(graph.node_id(sink), edge)?,
                    )
                }
                (
                    HedgeEdgeRole::Standard,
                    HedgePair::Split {
                        source,
                        sink,
                        split,
                    },
                ) => {
                    let reversed = reversed_dangling.contains(&edge);
                    if reversed {
                        reversed_boundary_edges.insert(edge);
                    }
                    let (vertex, flow) = match (split, reversed) {
                        (Flow::Source, false) => (graph.node_id(source), EdgeFlow::Outgoing),
                        (Flow::Source, true) => (graph.node_id(source), EdgeFlow::Incoming),
                        (Flow::Sink, false) => (graph.node_id(sink), EdgeFlow::Incoming),
                        (Flow::Sink, true) => (graph.node_id(sink), EdgeFlow::Outgoing),
                    };
                    CffEdge::boundary(edge, dense_vertex(vertex, edge)?, flow)
                }
                (HedgeEdgeRole::InitialState, HedgePair::Paired { source, sink }) => {
                    CffEdge::initial_state(
                        edge,
                        dense_vertex(graph.node_id(source), edge)?,
                        dense_vertex(graph.node_id(sink), edge)?,
                    )
                }
                (HedgeEdgeRole::InitialState, _) => {
                    return Err(CffError::Invariant(format!(
                        "initial-state edge {edge} is not paired"
                    )));
                }
                (HedgeEdgeRole::UnorientedExternal, _) => {
                    return Err(CffError::Invariant(format!(
                        "unoriented external edge {edge} is not unpaired"
                    )));
                }
                (HedgeEdgeRole::Omitted, _) => unreachable!("omitted edges are skipped"),
            };
            edges.push(cff_edge);
        }

        let vertex_sets = original_vertices
            .iter()
            .map(|vertex| VertexSet::from_usize(vertex.0));
        Ok(Self {
            graph: CffGraph::with_vertex_sets(vertex_sets, edges)?,
            represented_edges,
            internal_edges,
            reversed_boundary_edges,
            unoriented_external_edges,
        })
    }
}

impl CffGenerator {
    pub const fn new(options: CffOptions) -> Self {
        Self { options }
    }

    pub const fn options(&self) -> &CffOptions {
        &self.options
    }

    pub fn generate(&self, graph: &CffGraph) -> Result<CffResult, CffError> {
        let mut surfaces = SurfaceCache::default();
        let generation = self.generate_into(graph, &mut surfaces)?;
        Ok(CffResult {
            expression: generation.expression,
            surfaces,
            report: generation.report,
        })
    }

    /// Generate into a caller-owned arena so related expressions can share surface IDs.
    pub fn generate_into(
        &self,
        graph: &CffGraph,
        cache: &mut SurfaceCache,
    ) -> Result<CffGeneration, CffError> {
        self.validate_options(graph)?;
        let initial_surface_count = cache.len();

        let mut seed_orientation = GraphOrientation::from_iter(std::iter::repeat_n(
            EdgeOrientation::Undirected,
            graph.edge_count(),
        ));
        for edge in graph.edges() {
            seed_orientation[edge.id] = EdgeOrientation::Default;
        }
        let mut seed = WorkingGraph::from_graph(graph, seed_orientation.clone());
        for edge in &self.options.contracted_edges {
            seed.contract_edge(*edge)?;
            seed_orientation[*edge] = EdgeOrientation::Undirected;
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

        let mut orientations = Vec::new();
        for identifier in 0..candidate_count {
            let mut orientation = seed_orientation.clone();
            for (edge, fixed) in &self.options.fixed_orientations {
                orientation[*edge] = *fixed;
            }
            for (bit, edge) in free_edges.iter().enumerate() {
                let direction = if identifier & (1_usize << bit) == 0 {
                    EdgeOrientation::Default
                } else {
                    EdgeOrientation::Reversed
                };
                orientation[*edge] = direction;
            }

            let mut oriented_graph = seed.clone();
            oriented_graph.apply_orientation(orientation.clone());
            if oriented_graph.has_directed_cycle() {
                continue;
            }

            let id = OrientationId(orientations.len());
            let expression = self.generate_tree(oriented_graph, cache, identifier)?;
            orientations.push(OrientationExpression {
                id,
                data: OrientationData { orientation },
                expression,
            });
        }

        let expression = CffExpression::new(orientations);
        let report = CffReport {
            candidate_orientations: candidate_count,
            acyclic_orientations: expression.orientations().len(),
            unfolded_terms: expression.unfolded_term_count(),
            interned_surfaces: cache.len() - initial_surface_count,
        };
        Ok(CffGeneration { expression, report })
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
                if graph.edge(*edge).is_none() {
                    return Err(CffError::UnknownEdge(*edge));
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
        let mut vertices = graph
            .vertex_sets()
            .iter()
            .copied()
            .map(|vertices| WorkingVertex {
                vertices,
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
                EdgeKind::InitialState { source, sink } => {
                    let edge = WorkingEdge {
                        id: edge.id,
                        kind: WorkingEdgeKind::External,
                    };
                    Self::attach(&mut vertices[source.index()], edge, EdgeFlow::Outgoing);
                    Self::attach(&mut vertices[sink.index()], edge, EdgeFlow::Incoming);
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
            if *direction != EdgeOrientation::Reversed {
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
                vertex_set: vertex.vertices,
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
            vertex_set: vertex.vertices,
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
                    vertex_set: surface.vertex_set,
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
    use std::sync::Arc;

    use super::*;
    use crate::{CffEdge, VertexId};
    use feynkit_graph::{DiagramEdge, DiagramVertex, ExternalState, FeynmanDiagram};
    use feynkit_model::Model;
    use linnet::{
        half_edge::subgraph::{ModifySubSet, SuBitGraph},
        parser::DotGraph,
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
    fn shared_arena_preserves_existing_expression_ids_and_raised_data() {
        let generator = CffGenerator::default();
        let mut surfaces = SurfaceCache::default();
        let first = generator.generate_into(&bubble(), &mut surfaces).unwrap();
        let first_expression = first.expression.clone();
        let first_ids = first_expression
            .orientations()
            .iter()
            .flat_map(|orientation| orientation.expression.nodes())
            .map(|node| node.data)
            .collect::<BTreeSet<_>>();
        let first_resolved = first_ids
            .iter()
            .map(|id| (*id, surfaces.get(*id).unwrap()))
            .collect::<Vec<_>>();
        let representative = |edge| {
            if edge == EdgeId::new(3) {
                EdgeId::new(2)
            } else {
                edge
            }
        };
        let raised_before = first_expression
            .raised_energy_surfaces(&surfaces, representative)
            .unwrap();
        let first_surface_count = surfaces.len();

        let triangle = CffGraph::new(
            3,
            [
                CffEdge::internal(edge(4), vertex(0), vertex(1)),
                CffEdge::internal(edge(5), vertex(1), vertex(2)),
                CffEdge::internal(edge(6), vertex(0), vertex(2)),
            ],
        )
        .unwrap();
        generator.generate_into(&triangle, &mut surfaces).unwrap();

        assert!(surfaces.len() > first_surface_count);
        assert_eq!(first.expression, first_expression);
        assert_eq!(
            first_expression
                .raised_energy_surfaces(&surfaces, representative)
                .unwrap(),
            raised_before
        );
        for (id, surface) in first_resolved {
            assert_eq!(surfaces.get(id), Some(surface));
        }
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
        let model = Arc::new(
            Model::from_json(include_str!(
                "../../feynkit-model/tests/fixtures/scalars_2p_3p.json"
            ))
            .unwrap(),
        );
        let mut builder = FeynmanDiagram::builder(Arc::clone(&model), "bubble");
        let incoming =
            builder.add_vertex(DiagramVertex::external("p1", 0, ExternalState::Incoming));
        let outgoing =
            builder.add_vertex(DiagramVertex::external("p2", 1, ExternalState::Outgoing));
        let rule = model.vertex_rule_id("V_3_SCALAR_000").unwrap();
        let left = builder.add_vertex(DiagramVertex::interaction("left", rule));
        let right = builder.add_vertex(DiagramVertex::interaction("right", rule));
        let scalar = || DiagramEdge::new(model.particle_id("scalar_0").unwrap(), false);
        builder.add_edge(incoming, left, scalar()).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(left, right, scalar()).unwrap();
        builder.add_edge(right, outgoing, scalar()).unwrap();
        let diagram = builder.build().unwrap();

        let result = diagram.build_cff(CffOptions::default()).unwrap();

        assert_eq!(result.report.candidate_orientations, 4);
        assert_eq!(result.report.acyclic_orientations, 2);
        assert_eq!(result.report.unfolded_terms, 2);
    }

    #[test]
    fn linnet_subgraph_keeps_original_ids_boundary_flow_and_caller_cache() {
        let graph: DotGraph = DotGraph::from_string(
            r#"
            digraph triangle {
                ext_in [style=invis]
                ext_out [style=invis]
                A [id=0]
                B [id=1]
                C [id=2]
                A -> C [id=0]
                B -> C [id=1]
                A -> B [id=2]
                C -> ext_out [id=3]
                ext_in -> A [id=4]
            }
            "#,
        )
        .unwrap();
        let mut selected: SuBitGraph = graph.empty_subgraph();
        for (_, neighbours, vertex) in graph.iter_nodes() {
            if matches!(vertex.name(), Some("B" | "C")) {
                neighbours.for_each(|hedge| selected.add(hedge));
            }
        }

        let mut cache = SurfaceCache::default();
        cache.intern(Surface::Energy(EnergySurface {
            energies: vec![edge(99)],
            external_shift: ExternalShift::default(),
            vertex_set: VertexSet::dummy(),
        }));
        let result = graph
            .build_cff_from_subgraph(
                &selected,
                CffOptions::default().with_initial_state_edge(edge(1)),
                &[edge(2)],
                &mut cache,
            )
            .unwrap();

        assert_eq!(result.surfaces, cache);
        assert!(cache.len() > 1);
        for orientation in result.orientations() {
            assert_eq!(orientation.data.orientation[edge(2)], Orientation::Reversed);
            assert_eq!(
                orientation.data.orientation[edge(4)],
                Orientation::Undirected
            );
            assert_ne!(
                orientation.data.orientation[edge(1)],
                Orientation::Undirected
            );
        }
        assert!(cache.energy_surfaces().iter().skip(1).all(|surface| {
            !surface.energies.contains(&edge(1))
                && surface
                    .vertex_set
                    .iter()
                    .all(|vertex| [1, 2].contains(&vertex.0))
        }));
        assert!(cache.h_surfaces().iter().all(|surface| {
            !surface.positive_energies.contains(&edge(1))
                && !surface.negative_energies.contains(&edge(1))
                && surface
                    .vertex_set
                    .iter()
                    .all(|vertex| [1, 2].contains(&vertex.0))
        }));
        assert!(
            cache
                .energy_surfaces()
                .iter()
                .skip(1)
                .any(|surface| surface.external_shift.get(edge(1)) != 0)
                || cache
                    .h_surfaces()
                    .iter()
                    .any(|surface| surface.external_shift.get(edge(1)) != 0)
        );
    }

    #[test]
    fn linnet_edge_roles_cover_runtime_bookkeeping_without_an_adapter_graph() {
        let graph: DotGraph = DotGraph::from_string(
            r#"
            digraph bubble {
                ext_in [style=invis]
                ext_out [style=invis]
                ext_in -> A [id=2]
                B -> ext_out [id=3]
                A -> B [id=0]
                A -> B [id=1]
                A -> B [id=4]
            }
            "#,
        )
        .unwrap();
        let mut cache = SurfaceCache::default();
        let result = graph
            .build_cff_from_subgraph_with_edge_roles(
                &graph.full_filter(),
                CffOptions::default().with_contracted_edge(edge(0)),
                &[],
                |id| match id.index() {
                    1 => HedgeEdgeRole::InitialState,
                    2 | 3 => HedgeEdgeRole::UnorientedExternal,
                    4 => HedgeEdgeRole::Omitted,
                    _ => HedgeEdgeRole::Standard,
                },
                &mut cache,
            )
            .unwrap();

        assert_eq!(result.report.candidate_orientations, 1);
        assert_eq!(result.orientations().len(), 1);
        let orientation = &result.orientations()[0].data.orientation;
        assert_eq!(orientation[edge(0)], Orientation::Undirected);
        assert_eq!(orientation[edge(1)], Orientation::Default);
        assert_eq!(orientation[edge(2)], Orientation::Undirected);
        assert_eq!(orientation[edge(3)], Orientation::Undirected);
        assert_eq!(orientation[edge(4)], Orientation::Undirected);
    }
}
