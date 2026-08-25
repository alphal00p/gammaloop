use std::collections::{BTreeMap, BTreeSet};

use color_eyre::Result;
use eyre::eyre;
use feynkit_cff::{
    CffEdge, CffGenerator, CffGraph, CffOptions, EdgeFlow, EdgeId,
    EdgeOrientation as FeynkitEdgeOrientation, ExpressionTree, ExternalShift as FeynkitShift,
    GraphOrientation as FeynkitGraphOrientation, Surface as FeynkitSurface,
    SurfaceCache as FeynkitSurfaceCache, SurfaceId as FeynkitSurfaceId, VertexId,
};
use linnet::half_edge::{
    HedgeGraph, NodeIndex,
    involution::{EdgeIndex, EdgeVec, Flow, HedgePair, Orientation},
    subgraph::SubGraphLike,
};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cff_graph::VertexSet,
        esurface::{Esurface, EsurfaceID, add_external_shifts},
        expression::{
            CFFExpression, OrientationData,
            OrientationExpression as GammaLoopOrientationExpression, OrientationID,
        },
        generation::{ShiftRewrite, SurfaceCache},
        hsurface::{Hsurface, HsurfaceID},
        surface::HybridSurfaceID,
        tree::{NodeId, Tree},
    },
    graph::Graph,
    settings::global::OrientationPattern,
};

struct FeynkitCffInput {
    graph: CffGraph,
    options: CffOptions,
    dummy_edges: BTreeSet<EdgeIndex>,
    initial_edges: BTreeSet<EdgeIndex>,
    synthetic_initial_edges: BTreeMap<EdgeId, EdgeIndex>,
}

struct FeynkitSubgraphInput {
    graph: CffGraph,
    options: CffOptions,
    represented_edges: BTreeSet<EdgeIndex>,
    internal_edges: BTreeSet<EdgeIndex>,
    reversed_boundary_edges: BTreeSet<EdgeIndex>,
    original_vertices: Vec<NodeIndex>,
}

pub(super) fn generate_cff_expression_from_subgraph_with_feynkit<E, V, H, S: SubGraphLike>(
    graph: &HedgeGraph<E, V, H>,
    subgraph: &S,
    canonize_esurface: &Option<ShiftRewrite>,
    reversed_dangling: &[EdgeIndex],
    edges_in_initial_state_cut: &[EdgeIndex],
    surface_cache: &mut SurfaceCache,
) -> Result<CFFExpression<OrientationID>> {
    let FeynkitSubgraphInput {
        graph: cff_graph,
        options,
        represented_edges,
        internal_edges,
        reversed_boundary_edges,
        original_vertices,
    } = feynkit_subgraph_input(
        graph,
        subgraph,
        reversed_dangling,
        edges_in_initial_state_cut,
    )?;
    let result = CffGenerator::new(options).generate(&cff_graph)?;
    let synthetic_initial_edges = BTreeMap::new();
    let mut orientations = TiVec::new();

    for orientation in result.expression.orientations() {
        for edge in &internal_edges {
            match orientation.orientation.get(EdgeId::new(edge.0)) {
                Some(FeynkitEdgeOrientation::Default | FeynkitEdgeOrientation::Reversed) => {}
                Some(FeynkitEdgeOrientation::Undirected) | None => {
                    return Err(eyre!(
                        "FeynKit returned no directed orientation for represented edge {edge}"
                    ));
                }
            }
        }

        let global_orientation = graph.new_edgevec(|_, edge, _| {
            if !represented_edges.contains(&edge) {
                return Orientation::Undirected;
            }
            if reversed_boundary_edges.contains(&edge) {
                return Orientation::Reversed;
            }
            match orientation.orientation.get(EdgeId::new(edge.0)) {
                Some(FeynkitEdgeOrientation::Reversed) if internal_edges.contains(&edge) => {
                    Orientation::Reversed
                }
                _ => Orientation::Default,
            }
        });
        let expression = Graph::gammaloop_expression_tree(
            &orientation.expression,
            &result.surfaces,
            surface_cache,
            &synthetic_initial_edges,
            &original_vertices,
            canonize_esurface,
        )?;
        orientations.push(GammaLoopOrientationExpression {
            data: OrientationData {
                orientation: global_orientation,
            },
            expression,
        });
    }

    Ok(CFFExpression {
        orientations,
        surfaces: surface_cache.clone(),
    })
}

fn feynkit_subgraph_input<E, V, H, S: SubGraphLike>(
    graph: &HedgeGraph<E, V, H>,
    subgraph: &S,
    reversed_dangling: &[EdgeIndex],
    edges_in_initial_state_cut: &[EdgeIndex],
) -> Result<FeynkitSubgraphInput> {
    let mut original_vertices = graph
        .iter_nodes_of(subgraph)
        .map(|(vertex, _, _)| vertex)
        .collect::<Vec<_>>();
    original_vertices.sort_unstable();
    if let Some(vertex) = original_vertices
        .iter()
        .find(|vertex| vertex.0 >= u64::BITS as usize)
    {
        return Err(eyre!(
            "GammaLoop CFF surfaces support at most {} vertices, received vertex {}",
            u64::BITS,
            vertex.0
        ));
    }
    let dense_vertices = original_vertices
        .iter()
        .enumerate()
        .map(|(dense, original)| (*original, VertexId::new(dense)))
        .collect::<BTreeMap<_, _>>();
    let dense_vertex = |vertex: NodeIndex, edge: EdgeIndex| {
        dense_vertices.get(&vertex).copied().ok_or_else(|| {
            eyre!("edge {edge} refers to vertex {vertex:?} outside the selected subgraph")
        })
    };

    let reversed_dangling = reversed_dangling.iter().copied().collect::<BTreeSet<_>>();
    let mut edges = Vec::new();
    let mut represented_edges = BTreeSet::new();
    let mut internal_edges = BTreeSet::new();
    let mut reversed_boundary_edges = BTreeSet::new();
    for (pair, edge, _) in graph.iter_edges_of(subgraph) {
        let id = EdgeId::new(edge.0);
        represented_edges.insert(edge);
        let cff_edge = match pair {
            HedgePair::Unpaired { hedge, flow } => CffEdge::external(
                id,
                dense_vertex(graph.node_id(hedge), edge)?,
                match flow {
                    Flow::Sink => EdgeFlow::Incoming,
                    Flow::Source => EdgeFlow::Outgoing,
                },
            ),
            HedgePair::Paired { source, sink } => {
                internal_edges.insert(edge);
                CffEdge::internal(
                    id,
                    dense_vertex(graph.node_id(source), edge)?,
                    dense_vertex(graph.node_id(sink), edge)?,
                )
            }
            HedgePair::Split {
                source,
                sink,
                split,
            } => {
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
                CffEdge::boundary(id, dense_vertex(vertex, edge)?, flow)
            }
        };
        edges.push(cff_edge);
    }

    let cff_graph = CffGraph::new(original_vertices.len(), edges)?;
    let options = edges_in_initial_state_cut
        .iter()
        .filter(|edge| represented_edges.contains(edge))
        .fold(CffOptions::default(), |options, edge| {
            options.with_initial_state_edge(EdgeId::new(edge.0))
        });

    Ok(FeynkitSubgraphInput {
        graph: cff_graph,
        options,
        represented_edges,
        internal_edges,
        reversed_boundary_edges,
        original_vertices,
    })
}

impl Graph {
    pub(crate) fn generate_cff_with_feynkit(
        &mut self,
        contract_edges: &[EdgeIndex],
        canonize_esurface: &Option<ShiftRewrite>,
        orientation_pattern: &OrientationPattern,
    ) -> Result<CFFExpression<OrientationID>> {
        let FeynkitCffInput {
            graph,
            options,
            dummy_edges,
            initial_edges,
            synthetic_initial_edges,
        } = self.feynkit_cff_input(contract_edges)?;
        let result = CffGenerator::new(options).generate(&graph)?;
        let original_vertices = (0..self.n_nodes()).map(NodeIndex::from).collect::<Vec<_>>();
        let mut orientations = TiVec::new();

        for orientation in result.expression.orientations() {
            let global_orientation =
                self.gammaloop_orientation(&orientation.orientation, &dummy_edges, &initial_edges)?;
            if !orientation_pattern.filter(&global_orientation) {
                continue;
            }

            let expression = Self::gammaloop_expression_tree(
                &orientation.expression,
                &result.surfaces,
                &mut self.surface_cache,
                &synthetic_initial_edges,
                &original_vertices,
                canonize_esurface,
            )?;
            orientations.push(GammaLoopOrientationExpression {
                data: OrientationData {
                    orientation: global_orientation,
                },
                expression,
            });
        }

        Ok(CFFExpression {
            orientations,
            surfaces: self.surface_cache.clone(),
        })
    }

    fn feynkit_cff_input(&self, contract_edges: &[EdgeIndex]) -> Result<FeynkitCffInput> {
        let initial_edges = self
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge, _)| edge)
            .collect::<BTreeSet<_>>();
        let contracted_edges = contract_edges.iter().copied().collect::<BTreeSet<_>>();
        if let Some(edge) = contracted_edges.intersection(&initial_edges).next() {
            return Err(eyre!(
                "initial-state edge {edge} cannot be contracted during root CFF generation"
            ));
        }

        let mut edges = Vec::new();
        let mut dummy_edges = BTreeSet::new();
        let mut synthetic_initial_edges = BTreeMap::new();
        let mut next_synthetic_edge = match self.iter_edges().map(|(_, edge, _)| edge.0).max() {
            Some(maximum) => maximum
                .checked_add(1)
                .ok_or_else(|| eyre!("cannot allocate a synthetic initial-state edge ID"))?,
            None => 0,
        };
        for (pair, edge_id, edge_data) in self.iter_edges() {
            if edge_data.data.is_dummy {
                // Dummy edges are absent from the standalone CFF graph.
                dummy_edges.insert(edge_id);
                continue;
            }

            let edge = match pair {
                HedgePair::Unpaired { hedge, flow } => CffEdge::external(
                    EdgeId::new(edge_id.0),
                    VertexId::new(self.node_id(hedge).0),
                    match flow {
                        Flow::Sink => EdgeFlow::Incoming,
                        Flow::Source => EdgeFlow::Outgoing,
                    },
                ),
                HedgePair::Paired { source, sink } if initial_edges.contains(&edge_id) => {
                    // GammaLoop's root CFF graph treats a sewn initial-state edge as one
                    // external edge attached to both endpoints. FeynKit external edges have
                    // one endpoint, so represent both half-edges explicitly and map their
                    // shifts back to the original edge below. This also keeps the cut edge out
                    // of cycle detection, matching the historical CFF semantics.
                    let outgoing = EdgeId::new(next_synthetic_edge);
                    next_synthetic_edge = next_synthetic_edge.checked_add(1).ok_or_else(|| {
                        eyre!("cannot allocate a synthetic initial-state edge ID")
                    })?;
                    let incoming = EdgeId::new(next_synthetic_edge);
                    next_synthetic_edge = next_synthetic_edge.checked_add(1).ok_or_else(|| {
                        eyre!("cannot allocate a synthetic initial-state edge ID")
                    })?;
                    synthetic_initial_edges.insert(outgoing, edge_id);
                    synthetic_initial_edges.insert(incoming, edge_id);
                    edges.push(CffEdge::external(
                        outgoing,
                        VertexId::new(self.node_id(source).0),
                        EdgeFlow::Outgoing,
                    ));
                    CffEdge::external(
                        incoming,
                        VertexId::new(self.node_id(sink).0),
                        EdgeFlow::Incoming,
                    )
                }
                HedgePair::Paired { source, sink } => CffEdge::internal(
                    EdgeId::new(edge_id.0),
                    VertexId::new(self.node_id(source).0),
                    VertexId::new(self.node_id(sink).0),
                ),
                HedgePair::Split { .. } => {
                    return Err(eyre!(
                        "split edge {edge_id} cannot be represented by a root CFF graph"
                    ));
                }
            };
            edges.push(edge);
        }

        let graph = CffGraph::new(self.n_nodes(), edges)?;
        let mut options = CffOptions::default();
        for edge in contract_edges {
            options = options.with_contracted_edge(EdgeId::new(edge.0));
        }

        Ok(FeynkitCffInput {
            graph,
            options,
            dummy_edges,
            initial_edges,
            synthetic_initial_edges,
        })
    }

    fn gammaloop_orientation(
        &self,
        orientation: &FeynkitGraphOrientation,
        dummy_edges: &BTreeSet<EdgeIndex>,
        initial_edges: &BTreeSet<EdgeIndex>,
    ) -> Result<EdgeVec<Orientation>> {
        // Build a fresh orientation instead of mutating graph edges; applying a mutating
        // orientation twice could otherwise flip the same edge twice.
        for (pair, edge, _) in self.iter_edges() {
            if matches!(pair, HedgePair::Paired { .. })
                && !dummy_edges.contains(&edge)
                && !initial_edges.contains(&edge)
                && orientation.get(EdgeId::new(edge.0)).is_none()
            {
                return Err(eyre!(
                    "FeynKit returned no orientation for represented edge {edge}"
                ));
            }
        }

        Ok(self.new_edgevec(|_, edge, pair| {
            if pair.is_unpaired() || dummy_edges.contains(&edge) {
                return Orientation::Undirected;
            }
            if initial_edges.contains(&edge) {
                return Orientation::Default;
            }
            match orientation.get(EdgeId::new(edge.0)) {
                Some(FeynkitEdgeOrientation::Default) => Orientation::Default,
                Some(FeynkitEdgeOrientation::Reversed) => Orientation::Reversed,
                Some(FeynkitEdgeOrientation::Undirected) | None => Orientation::Undirected,
            }
        }))
    }

    fn gammaloop_expression_tree(
        expression: &ExpressionTree<FeynkitSurfaceId>,
        source_cache: &FeynkitSurfaceCache,
        target_cache: &mut SurfaceCache,
        synthetic_initial_edges: &BTreeMap<EdgeId, EdgeIndex>,
        original_vertices: &[NodeIndex],
        shift_rewrite: &Option<ShiftRewrite>,
    ) -> Result<Tree<HybridSurfaceID>> {
        let mut nodes = expression.nodes().iter();
        let root = nodes
            .next()
            .ok_or_else(|| eyre!("FeynKit returned an empty CFF expression tree"))?;
        if root.id.index() != 0 || root.parent.is_some() {
            return Err(eyre!("FeynKit returned a malformed CFF expression root"));
        }
        let mut tree = Tree::from_root(Self::intern_feynkit_surface(
            root.data,
            source_cache,
            target_cache,
            synthetic_initial_edges,
            original_vertices,
            shift_rewrite,
        )?);

        for (expected, node) in nodes.enumerate() {
            let expected = expected + 1;
            if node.id.index() != expected {
                return Err(eyre!(
                    "FeynKit CFF tree node {} appeared at position {expected}",
                    node.id.index()
                ));
            }
            let parent = node
                .parent
                .ok_or_else(|| eyre!("FeynKit CFF tree node {expected} has no parent"))?;
            if parent.index() >= expected {
                return Err(eyre!(
                    "FeynKit CFF tree node {expected} has non-ancestral parent {}",
                    parent.index()
                ));
            }
            tree.insert_node(
                NodeId::from(parent.index()),
                Self::intern_feynkit_surface(
                    node.data,
                    source_cache,
                    target_cache,
                    synthetic_initial_edges,
                    original_vertices,
                    shift_rewrite,
                )?,
            );
        }
        Ok(tree)
    }

    fn intern_feynkit_surface(
        id: FeynkitSurfaceId,
        source_cache: &FeynkitSurfaceCache,
        target_cache: &mut SurfaceCache,
        synthetic_initial_edges: &BTreeMap<EdgeId, EdgeIndex>,
        original_vertices: &[NodeIndex],
        shift_rewrite: &Option<ShiftRewrite>,
    ) -> Result<HybridSurfaceID> {
        match source_cache
            .get(id)
            .ok_or_else(|| eyre!("FeynKit CFF expression references missing surface {id:?}"))?
        {
            FeynkitSurface::Energy(surface) => {
                let mut surface = Esurface {
                    energies: surface
                        .energies
                        .into_iter()
                        .map(|edge| EdgeIndex(edge.index()))
                        .collect(),
                    external_shift: Self::gammaloop_shift(
                        &surface.external_shift,
                        synthetic_initial_edges,
                    ),
                    vertex_set: Self::gammaloop_vertex_set(surface.vertices, original_vertices)?,
                };
                if let Some(rewrite) = shift_rewrite {
                    surface.canonicalize_shift(rewrite);
                }
                let id = target_cache
                    .esurface_cache
                    .position(|cached| cached == &surface)
                    .unwrap_or_else(|| {
                        let id = EsurfaceID(target_cache.esurface_cache.len());
                        target_cache.esurface_cache.push(surface);
                        id
                    });
                Ok(HybridSurfaceID::Esurface(id))
            }
            FeynkitSurface::H(surface) => {
                let surface = Hsurface {
                    positive_energies: surface
                        .positive_energies
                        .into_iter()
                        .map(|edge| EdgeIndex(edge.index()))
                        .collect(),
                    negative_energies: surface
                        .negative_energies
                        .into_iter()
                        .map(|edge| EdgeIndex(edge.index()))
                        .collect(),
                    external_shift: Self::gammaloop_shift(
                        &surface.external_shift,
                        synthetic_initial_edges,
                    ),
                    vertex_set: Self::gammaloop_vertex_set(surface.vertices, original_vertices)?,
                };
                let id = target_cache
                    .hsurface_cache
                    .position(|cached| cached == &surface)
                    .unwrap_or_else(|| {
                        let id = HsurfaceID::from(target_cache.hsurface_cache.len());
                        target_cache.hsurface_cache.push(surface);
                        id
                    });
                Ok(HybridSurfaceID::Hsurface(id))
            }
            FeynkitSurface::Unit => Ok(HybridSurfaceID::Unit),
            FeynkitSurface::Infinite => Ok(HybridSurfaceID::Infinite),
        }
    }

    fn gammaloop_shift(
        shift: &FeynkitShift,
        synthetic_initial_edges: &BTreeMap<EdgeId, EdgeIndex>,
    ) -> Vec<(EdgeIndex, i64)> {
        let mapped = shift
            .iter()
            .map(|(edge, coefficient)| {
                (
                    synthetic_initial_edges
                        .get(&edge)
                        .copied()
                        .unwrap_or(EdgeIndex(edge.index())),
                    coefficient,
                )
            })
            .collect();
        add_external_shifts(&Vec::new(), &mapped)
    }

    fn gammaloop_vertex_set(
        vertices: feynkit_cff::VertexSet,
        original_vertices: &[NodeIndex],
    ) -> Result<VertexSet> {
        vertices
            .iter()
            .try_fold(VertexSet::dummy(), |vertices, vertex| {
                let original = original_vertices.get(vertex.index()).ok_or_else(|| {
                    eyre!(
                        "FeynKit CFF surface references missing dense vertex {}",
                        vertex.index()
                    )
                })?;
                if original.0 >= u64::BITS as usize {
                    return Err(eyre!(
                        "GammaLoop CFF surfaces support at most {} vertices, received vertex {}",
                        u64::BITS,
                        original.0
                    ));
                }
                Ok(vertices.join(&VertexSet::from_usize(original.0)))
            })
    }
}

#[cfg(test)]
mod tests {
    use feynkit_cff::CffError;
    use linnet::half_edge::{
        involution::EdgeIndex,
        subgraph::{ModifySubSet, SuBitGraph},
    };

    use crate::{
        cff::{
            cff_graph::VertexSet,
            esurface::Esurface,
            generation::{ShiftRewrite, generate_cff_expression_from_subgraph},
            hsurface::Hsurface,
        },
        dot,
        graph::{Graph, parse::IntoGraph},
        initialisation::test_initialise,
        settings::global::OrientationPattern,
    };

    fn bubble() -> Graph {
        dot!(
            digraph bubble {
                edge [particle=scalar_1]
                node [num=1]
                ext [style=invis]
                ext -> A:0 [id=3]
                B:1 -> ext [id=2]
                A -> B [id=1]
                A -> B [id=0]
            },
            "scalars"
        )
        .unwrap()
    }

    fn cross_section_bubble() -> Graph {
        dot!(
            digraph cross_section_bubble {
                edge [particle=scalar_1]
                node [num=1]
                ext [style=invis is_cut=0]
                ext -> A:0 [id=3]
                B:1 -> ext [id=2]
                A -> B [id=1]
                A -> B [id=0]
            },
            "scalars"
        )
        .unwrap()
    }

    fn triangle() -> Graph {
        dot!(
            digraph triangle {
                edge [particle=scalar_1]
                node [num=1]
                ext [style=invis]
                ext -> A:0 [id=4]
                C:1 -> ext [id=3]
                A -> B [id=2]
                B -> C [id=1]
                A -> C [id=0]
            },
            "scalars"
        )
        .unwrap()
    }

    fn triangle_left_subgraph(graph: &Graph) -> SuBitGraph {
        let mut subgraph: SuBitGraph = graph.empty_subgraph();
        for (_, neighbours, vertex) in graph.iter_nodes() {
            if vertex.name.as_str() == "A" || vertex.name.as_str() == "B" {
                neighbours.for_each(|hedge| subgraph.add(hedge));
            }
        }
        subgraph
    }

    #[test]
    fn root_amplitude_default_snapshot() {
        test_initialise().unwrap();
        let mut graph = bubble();
        let expression = graph
            .generate_cff(&[], &None, &OrientationPattern::default())
            .unwrap();

        insta::assert_json_snapshot!(
            "root_amplitude_default",
            serde_json::json!({
                "expression": expression,
                "surface_cache": graph.surface_cache,
            })
        );
    }

    #[test]
    fn root_amplitude_selected_orientation_snapshot() {
        test_initialise().unwrap();
        // FeynKit enumerates internal-edge orientations from a bitset: unset bits retain
        // the graph direction and set bits reverse it. Fixed external and boundary directions
        // pad the GammaLoop orientation, and cyclic orientations are filtered before conversion.
        let selected = OrientationPattern::from_user_pattern("(1,1,0,0)").unwrap();
        let mut graph = bubble();
        let expression = graph.generate_cff(&[], &None, &selected).unwrap();
        assert_eq!(expression.orientations.len(), 1);

        insta::assert_json_snapshot!(
            "root_amplitude_selected_orientation",
            serde_json::json!({
                "expression": expression,
                "surface_cache": graph.surface_cache,
            })
        );
    }

    #[test]
    fn root_tree_snapshot() {
        test_initialise().unwrap();
        let mut graph = triangle();
        let expression = graph
            .generate_cff(&[], &None, &OrientationPattern::default())
            .unwrap();

        insta::assert_json_snapshot!(
            "root_tree",
            serde_json::json!({
                "expression": expression,
                "surface_cache": graph.surface_cache,
            })
        );
    }

    #[test]
    fn root_contracted_edge_snapshot() {
        test_initialise().unwrap();
        let mut graph = triangle();
        let expression = graph
            .generate_cff(&[EdgeIndex(1)], &None, &OrientationPattern::default())
            .unwrap();

        insta::assert_json_snapshot!(
            "root_contracted_edge",
            serde_json::json!({
                "expression": expression,
                "surface_cache": graph.surface_cache,
            })
        );
    }

    #[test]
    fn root_initial_state_shift_and_cache_snapshot() {
        test_initialise().unwrap();
        let mut graph = cross_section_bubble();
        graph.surface_cache.esurface_cache.push(Esurface {
            energies: vec![EdgeIndex(99)],
            external_shift: vec![],
            vertex_set: VertexSet::dummy(),
        });
        graph.surface_cache.hsurface_cache.push(Hsurface {
            positive_energies: vec![EdgeIndex(98)],
            negative_energies: vec![EdgeIndex(99)],
            external_shift: vec![],
            vertex_set: VertexSet::dummy(),
        });
        let initial_edge = graph
            .iter_edges_of(&graph.initial_state_cut)
            .next()
            .unwrap()
            .1;
        let shift_rewrite = Some(ShiftRewrite {
            dependent_momentum: initial_edge,
            dependent_momentum_expr: vec![],
        });

        let expression = graph
            .generate_cff(&[], &shift_rewrite, &OrientationPattern::default())
            .unwrap();
        insta::assert_json_snapshot!(
            "root_initial_state_shift_and_cache",
            serde_json::json!({
                "expression": expression,
                "surface_cache": graph.surface_cache,
            })
        );
    }

    #[test]
    fn root_generation_rejects_contracted_initial_state_edges() {
        test_initialise().unwrap();
        let mut graph = cross_section_bubble();
        let initial_edge = graph
            .iter_edges_of(&graph.initial_state_cut)
            .next()
            .unwrap()
            .1;

        let error = graph
            .generate_cff(&[initial_edge], &None, &OrientationPattern::default())
            .unwrap_err();
        assert_eq!(
            error.to_string(),
            format!(
                "initial-state edge {initial_edge} cannot be contracted during root CFF generation"
            )
        );
    }

    #[test]
    fn subgraph_reversed_boundary_snapshot() {
        test_initialise().unwrap();
        let graph = triangle();
        let subgraph = triangle_left_subgraph(&graph);
        let reversed_dangling = [EdgeIndex(0)];
        let initial_edges = graph.get_edges_in_initial_state_cut();
        let mut surface_cache = graph.surface_cache.clone();
        let expression = generate_cff_expression_from_subgraph(
            &graph.underlying,
            &subgraph,
            &None,
            &reversed_dangling,
            &initial_edges,
            &mut surface_cache,
        )
        .unwrap();

        insta::assert_json_snapshot!(
            "subgraph_reversed_boundary",
            serde_json::json!({
                "expression": expression,
                "surface_cache": surface_cache,
            })
        );
    }

    #[test]
    fn empty_subgraph_reports_the_typed_feynkit_error() {
        test_initialise().unwrap();
        let graph = triangle();
        let empty: SuBitGraph = graph.empty_subgraph();
        let mut cache = graph.surface_cache.clone();

        let error = generate_cff_expression_from_subgraph(
            &graph.underlying,
            &empty,
            &None,
            &[],
            &[],
            &mut cache,
        )
        .unwrap_err();

        assert_eq!(
            error.downcast_ref::<CffError>(),
            Some(&CffError::EmptyGraph)
        );
    }
}
