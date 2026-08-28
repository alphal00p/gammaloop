use std::{
    collections::{BTreeMap, BTreeSet},
    ops::Deref,
    path::{Path, PathBuf},
};

use crate::{
    graph::{
        FinalizedCut, GraphGroup, GroupId, LoopMomentumBasis,
        attribute_warnings::warn_about_unknown_attributes, edge::EdgeExtraData,
    },
    integrands::process::ParamBuilder,
    model::{Model, ParticleId, ParticleIdGammaLoopExt},
    momentum::{
        SignOrZero,
        sample::LoopIndex,
        signature::{LoopExtSignature, SignatureLike},
    },
    numerator::{
        GlobalPrefactor,
        aind::{Aind, NewAind},
    },
    processes::DotExportSettings,
    utils::symbolica_ext::DOD,
    uv::UltravioletGraph,
};
use feynkit_cff::SurfaceCache;
use feynkit_graph::{
    DiagramEndpoint as FeynkitDiagramEndpoint, DiagramHalfEdge as FeynkitDiagramHalfEdge,
    EdgeId as FeynkitEdgeId, ExternalState as FeynkitExternalState, FeynmanDiagram,
    VertexId as FeynkitVertexId,
};
use idenso::{
    color::{ColorSimplifier, ColorSimplifySettings},
    dirac::PS,
    tensor::SymbolicNetParse,
};
use spenso::shadowing::symbolica_utils::LogPrint;

use color_eyre::{Report, Result, Section};

use eyre::{Context, Ok, eyre};
use itertools::Itertools;
use linnet::{
    half_edge::{
        HedgeGraph, NodeIndex,
        builder::HedgeGraphBuilder,
        involution::{EdgeData, EdgeIndex, Flow, Hedge, HedgePair, Orientation},
        nodestore::NodeStorageVec,
        subgraph::{ModifySubSet, OrientedCut, SuBitGraph, SubSetOps},
        swap::Swap,
    },
    parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GraphSet, HedgeParseError},
    permutation::Permutation,
};
use spenso::{network::parsing::ParseSettings, structure::slot::IsAbstractSlot};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder},
    symbol,
};
use tracing::instrument;
use tracing::{debug, warn};
use typed_index_collections::TiVec;

use super::{
    Autogen, Edge, Graph, HedgeData, LMBext, Vertex,
    edge::{EdgeMass, ParseEdge},
    global::ParseData,
    hedge_data::{NumIndices, ParseHedgeData},
    vertex::ParseVertex,
};

/// Extract oriented particles from hedges, filtering out dummy edges
pub fn extract_oriented_particles_from_vertex_hedges<I, V>(
    graph: &HedgeGraph<ParseEdge, V, ParseHedgeData>,
    hedges: I,
    model: &Model,
) -> Vec<ParticleId>
where
    I: Iterator<Item = Hedge>,
{
    hedges
        .filter_map(|h| {
            let eid = graph[&h];
            if graph[eid].is_dummy {
                return None;
            }
            let particle = graph[eid].particle.particle()?;
            Some(if graph.flow(h) != Flow::Sink {
                particle.antiparticle(model)
            } else {
                particle
            })
        })
        .collect()
}

// Type aliases for cleaner code
type NumGraph = HedgeGraph<ParseEdge, ParseVertex, HedgeData>;
type UnderlyingGraph = HedgeGraph<Edge, Vertex, HedgeData>;

pub mod string_utils;
pub use string_utils::{FromStripedStr, StripParse, ToQuoted};

#[derive(Clone, Debug)]
pub struct ParseGraph {
    pub global_data: ParseData,
    pub graph: HedgeGraph<ParseEdge, ParseVertex, ParseHedgeData>,
}

impl Deref for ParseGraph {
    type Target = HedgeGraph<ParseEdge, ParseVertex, ParseHedgeData>;

    fn deref(&self) -> &Self::Target {
        &self.graph
    }
}

impl ParseGraph {
    pub fn n_fermion_loops(&self, model: &Model) -> usize {
        let fermions: SuBitGraph = self
            .graph
            .from_filter(|edge| edge.particle.is_fermion(model));

        self.graph.cyclotomatic_number(&fermions)
    }
    pub fn n_external_fermion_loops(&mut self, model: &Model) -> Result<usize> {
        let internal = self.n_fermion_loops(model);
        self.graph
            .sew(
                |_, ae, _, be| {
                    if let (Some(a), Some(b)) = (ae.data.is_cut, be.data.is_cut) {
                        a == b
                    } else {
                        false
                    }
                },
                |af, ae, bf, be| match (af, bf) {
                    (Flow::Sink, Flow::Source) => (Flow::Sink, ae),
                    (Flow::Source, Flow::Sink) => (Flow::Source, be),
                    _ => panic!("Cannot sew hedges with flow {:?} and {:?}", af, bf),
                },
            )
            .map_err(|e| eyre::eyre!("Graph sewing failed: {:?}", e))?;

        Ok(self.n_fermion_loops(model) - internal)
    }

    pub fn debug_dot(&self) -> String {
        DotGraph::from(self).debug_dot()
    }
    /// Return the explicit UFO slot recorded for every half-edge.
    ///
    /// A finalized runtime artifact must retain the rule-leg assignment made
    /// by FeynKit. Inferring it from particles here would reintroduce a second
    /// physics-generation path and is therefore rejected.
    pub(crate) fn hedge_order(&self) -> Result<Vec<u8>> {
        (0..self.n_hedges())
            .map(|index| {
                self.graph[Hedge(index)].ufo_order.ok_or_else(|| {
                    eyre!(
                        "finalized runtime DOT graph '{}' is missing ufo_order for half-edge {index}",
                        self.global_data.name
                    )
                })
            })
            .collect()
    }

    /// Build the mechanical GammaLoop half-edge representation of a finalized
    /// FeynKit diagram.
    ///
    /// This conversion never resolves a vertex rule or regenerates a numerator.
    /// FeynKit has already made those physics decisions. The temporary parse
    /// graph is used only to reuse the mature cut sewing and runtime-cache
    /// construction below.
    fn from_feynkit_diagram(
        diagram: &FeynmanDiagram,
        model: &Model,
        external_connections: &[(Option<usize>, Option<usize>)],
    ) -> Result<Self> {
        #[derive(Clone, Copy)]
        struct ExternalEdge {
            edge: FeynkitEdgeId,
            internal: FeynkitVertexId,
            internal_flow: Flow,
        }

        fn particle_orientation(
            model: &Model,
            particle: feynkit_model::ParticleId,
            directed: bool,
        ) -> Result<Orientation> {
            let particle = model.particle_by_id(particle)?;
            Ok(if !directed {
                Orientation::Undirected
            } else if particle.is_antiparticle() {
                Orientation::Reversed
            } else {
                Orientation::Default
            })
        }

        fn external_index(tag: usize) -> Result<usize> {
            tag.checked_sub(1)
                .ok_or_else(|| eyre!("external connection tags are one-based, found zero"))
        }

        let mut builder = HedgeGraphBuilder::new();
        let mut vertex_map = BTreeMap::new();
        for (vertex_id, vertex) in diagram.vertices() {
            if vertex.external.is_some() {
                continue;
            }
            let runtime = builder.add_node(ParseVertex {
                name: Some(vertex.name.clone()),
                feynkit_id: Some(vertex_id),
                vertex_rule: vertex.interaction,
                num: Some(vertex.numerator.clone()),
                dod: None,
            });
            vertex_map.insert(vertex_id, runtime);
        }

        let mut external_edges = BTreeMap::<usize, ExternalEdge>::new();
        let mut edges = BTreeMap::new();
        for (edge_id, endpoints, edge) in diagram.edges() {
            edges.insert(edge_id, (endpoints, edge));
            let source_external = diagram
                .vertex(endpoints.source)
                .and_then(|vertex| vertex.external.as_ref());
            let target_external = diagram
                .vertex(endpoints.target)
                .and_then(|vertex| vertex.external.as_ref());
            let external = match (source_external, target_external) {
                (Some(external), None) => Some((
                    external,
                    ExternalEdge {
                        edge: edge_id,
                        internal: endpoints.target,
                        internal_flow: Flow::Sink,
                    },
                )),
                (None, Some(external)) => Some((
                    external,
                    ExternalEdge {
                        edge: edge_id,
                        internal: endpoints.source,
                        internal_flow: Flow::Source,
                    },
                )),
                (None, None) => None,
                (Some(_), Some(_)) => {
                    return Err(eyre!(
                        "FeynKit diagram '{}' has an edge between two external vertices",
                        diagram.name()
                    ));
                }
            };
            if let Some((external, edge)) = external
                && external_edges.insert(external.index, edge).is_some()
            {
                return Err(eyre!(
                    "FeynKit diagram '{}' has duplicate external index {}",
                    diagram.name(),
                    external.index
                ));
            }
        }

        let mut seen_edges = BTreeSet::new();
        let add_external =
            |external: ExternalEdge,
             flow: Flow,
             cut: Option<Hedge>,
             seen_edges: &mut BTreeSet<FeynkitEdgeId>,
             builder: &mut HedgeGraphBuilder<ParseEdge, ParseVertex, ParseHedgeData>|
             -> Result<()> {
                if !seen_edges.insert(external.edge) {
                    return Err(eyre!(
                        "external FeynKit edge {} is connected more than once",
                        external.edge.0
                    ));
                }
                let (_, edge) = edges
                    .get(&external.edge)
                    .ok_or_else(|| eyre!("missing FeynKit edge {}", external.edge.0))?;
                let mut orientation = particle_orientation(model, edge.particle, edge.directed)?;
                let mut data = ParseEdge::new(edge.particle)
                    .with_label(format!("feynkit_edge_{}", external.edge.0))
                    .with_num(edge.numerator.clone());
                data.feynkit_id = Some(external.edge);
                data.is_cut = cut;

                // A dangling runtime edge carries the process flow. If the
                // FeynKit internal half-edge has the opposite flow, reverse only
                // the particle representation; the finalized symbolic numerator
                // is translated independently below.
                if external.internal_flow != flow {
                    data.particle = data.particle.reverse(model);
                    orientation = orientation.reverse();
                }
                let node = *vertex_map.get(&external.internal).ok_or_else(|| {
                    eyre!(
                        "external FeynKit edge {} is not attached to an internal vertex",
                        external.edge.0
                    )
                })?;
                builder.add_external_edge(node, data, orientation, flow);
                Ok(())
            };

        // Keep the same cut ordering expected by the runtime: incoming cut
        // hedges first, ordinary outgoing amplitude legs next, internal edges,
        // and finally the outgoing partners of initial-state cuts.
        for (connection, (incoming, outgoing)) in external_connections.iter().enumerate() {
            if let Some(incoming) = incoming {
                let index = external_index(*incoming)?;
                let edge = *external_edges
                    .get(&index)
                    .ok_or_else(|| eyre!("missing incoming FeynKit external leg {index}"))?;
                let state = diagram
                    .vertex(if edge.internal_flow == Flow::Sink {
                        edges[&edge.edge].0.source
                    } else {
                        edges[&edge.edge].0.target
                    })
                    .and_then(|vertex| vertex.external.as_ref())
                    .map(|external| external.state);
                if state != Some(FeynkitExternalState::Incoming) {
                    return Err(eyre!("FeynKit external leg {index} is not incoming"));
                }
                add_external(
                    edge,
                    Flow::Sink,
                    outgoing.map(|_| Hedge(connection)),
                    &mut seen_edges,
                    &mut builder,
                )?;
            }
        }
        for (incoming, outgoing) in external_connections {
            if incoming.is_none()
                && let Some(outgoing) = outgoing
            {
                let index = external_index(*outgoing)?;
                let edge = *external_edges
                    .get(&index)
                    .ok_or_else(|| eyre!("missing outgoing FeynKit external leg {index}"))?;
                add_external(edge, Flow::Source, None, &mut seen_edges, &mut builder)?;
            }
        }
        for (edge_id, (endpoints, edge)) in &edges {
            if seen_edges.contains(edge_id) {
                continue;
            }
            if !vertex_map.contains_key(&endpoints.source)
                || !vertex_map.contains_key(&endpoints.target)
            {
                continue;
            }
            let orientation = particle_orientation(model, edge.particle, edge.directed)?;
            let mut data = ParseEdge::new(edge.particle)
                .with_label(format!("feynkit_edge_{}", edge_id.0))
                .with_num(edge.numerator.clone());
            data.feynkit_id = Some(*edge_id);
            builder.add_edge(
                vertex_map[&endpoints.source],
                vertex_map[&endpoints.target],
                data,
                orientation,
            );
            seen_edges.insert(*edge_id);
        }
        for (connection, (incoming, outgoing)) in external_connections.iter().enumerate() {
            if incoming.is_some()
                && let Some(outgoing) = outgoing
            {
                let index = external_index(*outgoing)?;
                let edge = *external_edges
                    .get(&index)
                    .ok_or_else(|| eyre!("missing outgoing FeynKit external leg {index}"))?;
                add_external(
                    edge,
                    Flow::Source,
                    Some(Hedge(connection)),
                    &mut seen_edges,
                    &mut builder,
                )?;
            }
        }
        if seen_edges.len() != edges.len() {
            let missing = edges
                .keys()
                .filter(|edge| !seen_edges.contains(edge))
                .map(|edge| edge.0)
                .collect_vec();
            return Err(eyre!(
                "FeynKit diagram '{}' contains edges outside the requested external connections: {missing:?}",
                diagram.name()
            ));
        }

        let mut graph: HedgeGraph<ParseEdge, ParseVertex, ParseHedgeData> = builder.into();
        let mut hedge_orders = Vec::new();
        for (pair, runtime_edge, edge) in graph.iter_edges() {
            let feynkit_id = edge.data.feynkit_id.ok_or_else(|| {
                eyre!("runtime edge {runtime_edge} lost its transient FeynKit identity")
            })?;
            let (endpoints, source_edge) = edges
                .get(&feynkit_id)
                .ok_or_else(|| eyre!("unknown FeynKit edge {}", feynkit_id.0))?;
            let source_order = u8::try_from(source_edge.source_slot().0)
                .map_err(|_| eyre!("FeynKit source slot does not fit in u8"))?;
            let target_order = u8::try_from(source_edge.target_slot().0)
                .map_err(|_| eyre!("FeynKit target slot does not fit in u8"))?;
            match pair {
                HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                    hedge_orders.push((source, source_order));
                    hedge_orders.push((sink, target_order));
                }
                HedgePair::Unpaired { hedge, .. } => {
                    let order = if diagram
                        .vertex(endpoints.source)
                        .is_some_and(|vertex| !vertex.is_external())
                    {
                        source_order
                    } else {
                        target_order
                    };
                    hedge_orders.push((hedge, order));
                }
            }
        }
        for (hedge, order) in hedge_orders {
            graph[hedge].ufo_order = Some(order);
        }

        Ok(Self {
            global_data: ParseData {
                name: diagram.name().to_owned(),
                overall_factor: diagram.overall_factor().clone(),
                projectors: Some(diagram.projector().clone()),
                num: diagram.numerator_prefactor().clone(),
                ..Default::default()
            },
            graph,
        })
    }
}

impl ParseGraph {
    pub(crate) fn from_parsed(graph: DotGraph, model: &Model) -> Result<Self> {
        warn_about_unknown_attributes(&graph);
        if graph
            .global_data
            .statements
            .contains_key("canonical_cuts_required")
        {
            return Err(eyre!(
                "cross-section runtime DOT does not carry canonical physical cuts; import the canonical FeynmanDiagram DOT artifact through FeynKit instead"
            ));
        }
        let global_data = graph.global_data.into();
        let graph = graph
            .graph
            .map_data_ref_result(
                |_, _, v| Ok(v),
                ParseEdge::parse(model),
                ParseHedgeData::parse(),
            )?
            .map_data_ref_result(
                ParseVertex::parse(model, &global_data),
                |_, _, _, e| Ok(e.map(Clone::clone)),
                |(_, h)| Ok(h.clone()),
            )?;

        Ok(Self { graph, global_data })
    }
}

/// Helper struct to hold initial data extracted from ParseGraph
struct InitialGraphData {
    overall_factor: Atom,
    global_prefactor: GlobalPrefactor,
    additional_params: Vec<Atom>,
    add_polarizations: bool,
    group_id: Option<GroupId>,
    is_group_master: bool,
    name: String,
}

/// Result of processing cut edges
struct CutProcessingResult {
    lmb_ids: BTreeMap<LoopIndex, EdgeIndex>,
    xs_ext_id: BTreeMap<Hedge, (EdgeIndex, Hedge)>,
    initial_hedges: SuBitGraph,
    full_cut: SuBitGraph,
}

impl CutProcessingResult {
    fn permute(&mut self, graph: &mut NumGraph) -> Result<()> {
        let (h_perm, edge_perm): (Vec<_>, Vec<_>) = self
            .xs_ext_id
            .iter()
            .enumerate()
            .map(|(target_pos, (_, (edge_idx, h_id)))| {
                ((h_id.0, target_pos), (edge_idx.0, target_pos))
            })
            .unzip();

        let per = Permutation::from_mappings(edge_perm, graph.n_edges()).unwrap();
        let perh = Permutation::from_mappings(h_perm, graph.n_hedges()).unwrap();

        debug!("Before: {}", graph.dot(&self.initial_hedges));
        <HedgeGraph<_, _, _> as Swap<Hedge>>::permute(graph, &perh);
        let trans = perh.transpositions();

        for (i, j) in trans.into_iter().rev() {
            self.full_cut.swap(Hedge(i), Hedge(j));
            // self.initial_hedges.swap(i, j);// initial hedges is already assuming permuted hedges
        }

        debug!("Before after: {}", graph.dot(&self.initial_hedges));
        <HedgeGraph<_, _, _> as Swap<EdgeIndex>>::permute(graph, &per);

        debug!(" after: {}", graph.dot(&self.initial_hedges));
        Ok(())
    }
}

/// Mechanical correspondence between the finalized FeynKit graph and the
/// post-sewing GammaLoop half-edge storage. Physics expressions are translated
/// through this table exactly once; no rule lookup or numerator regeneration is
/// permitted on this path.
struct FeynkitRuntimeMap {
    edges: BTreeMap<FeynkitEdgeId, EdgeIndex>,
    vertices: BTreeMap<FeynkitVertexId, NodeIndex>,
    half_edges: BTreeMap<(FeynkitEdgeId, Flow), Hedge>,
    momentum_signs: BTreeMap<FeynkitEdgeId, i8>,
}

fn feynkit_external_edges(diagram: &FeynmanDiagram) -> BTreeMap<usize, FeynkitEdgeId> {
    let mut result = BTreeMap::new();
    for (edge_id, endpoints, _) in diagram.edges() {
        let external = diagram
            .vertex(endpoints.source)
            .and_then(|vertex| vertex.external.as_ref())
            .or_else(|| {
                diagram
                    .vertex(endpoints.target)
                    .and_then(|vertex| vertex.external.as_ref())
            });
        if let Some(external) = external {
            result.insert(external.index, edge_id);
        }
    }
    result
}

fn feynkit_external_connections(
    diagram: &FeynmanDiagram,
) -> Result<Vec<(Option<usize>, Option<usize>)>> {
    let mut connections = BTreeMap::<usize, (Option<usize>, Option<usize>)>::new();
    for (_, vertex) in diagram.vertices() {
        let Some(external) = &vertex.external else {
            continue;
        };
        let tag = external
            .index
            .checked_add(1)
            .ok_or_else(|| eyre!("external index {} cannot be incremented", external.index))?;
        let connection = connections.entry(external.connection).or_default();
        match external.state {
            FeynkitExternalState::Incoming => connection.0 = Some(tag),
            FeynkitExternalState::Outgoing => connection.1 = Some(tag),
        }
    }
    Ok(connections.into_values().collect())
}

fn feynkit_tag_index(tag: usize) -> Result<usize> {
    tag.checked_sub(1)
        .ok_or_else(|| eyre!("external connection tags are one-based, found zero"))
}

fn feynkit_runtime_map(
    diagram: &FeynmanDiagram,
    graph: &NumGraph,
    external_connections: &[(Option<usize>, Option<usize>)],
) -> Result<FeynkitRuntimeMap> {
    let mut vertices = BTreeMap::new();
    for (node, _, vertex) in graph.iter_nodes() {
        if let Some(id) = vertex.feynkit_id {
            vertices.insert(id, node);
        }
    }

    let mut edges = BTreeMap::new();
    for (_, edge_id, edge) in graph.iter_edges() {
        if let Some(id) = edge.data.feynkit_id {
            edges.insert(id, edge_id);
        }
    }

    // Sewing retains one of the two cut-edge payloads. Restore both canonical
    // FeynKit IDs as aliases of the resulting runtime edge.
    let external_edges = feynkit_external_edges(diagram);
    for (connection, (incoming, outgoing)) in external_connections.iter().enumerate() {
        if let (Some(incoming), Some(outgoing)) = (incoming, outgoing) {
            let incoming = external_edges
                .get(&feynkit_tag_index(*incoming)?)
                .ok_or_else(|| eyre!("missing incoming FeynKit external edge"))?;
            let outgoing = external_edges
                .get(&feynkit_tag_index(*outgoing)?)
                .ok_or_else(|| eyre!("missing outgoing FeynKit external edge"))?;
            let runtime = EdgeIndex::from(connection);
            edges.insert(*incoming, runtime);
            edges.insert(*outgoing, runtime);
        }
    }

    let choose_hedge = |vertex: FeynkitVertexId, flow: Flow, pair: HedgePair| -> Result<Hedge> {
        let candidates = match pair {
            HedgePair::Paired { source, sink } | HedgePair::Split { source, sink, .. } => {
                vec![source, sink]
            }
            HedgePair::Unpaired { hedge, .. } => vec![hedge],
        };
        if let Some(runtime_node) = vertices.get(&vertex) {
            let attached = candidates
                .iter()
                .copied()
                .filter(|hedge| graph.node_id(*hedge) == *runtime_node)
                .collect_vec();
            if let Some(hedge) = attached
                .iter()
                .copied()
                .find(|hedge| graph.flow(*hedge) == flow)
                .or_else(|| attached.first().copied())
            {
                return Ok(hedge);
            }
        }
        candidates
            .iter()
            .copied()
            .find(|hedge| graph.flow(*hedge) == flow)
            .or_else(|| candidates.first().copied())
            .ok_or_else(|| eyre!("runtime edge has no half-edge"))
    };

    let mut half_edges = BTreeMap::new();
    let mut momentum_signs = BTreeMap::new();
    for (edge_id, endpoints, _) in diagram.edges() {
        let runtime = *edges.get(&edge_id).ok_or_else(|| {
            eyre!(
                "finalized FeynKit edge {} was lost during mechanical sewing",
                edge_id.0
            )
        })?;
        let pair = graph[&runtime].1;
        let source = choose_hedge(endpoints.source, Flow::Source, pair)?;
        let sink = choose_hedge(endpoints.target, Flow::Sink, pair)?;
        half_edges.insert((edge_id, Flow::Source), source);
        half_edges.insert((edge_id, Flow::Sink), sink);

        let (internal_flow, internal_hedge) = if diagram
            .vertex(endpoints.source)
            .is_some_and(|vertex| !vertex.is_external())
        {
            (Flow::Source, source)
        } else {
            (Flow::Sink, sink)
        };
        momentum_signs.insert(
            edge_id,
            if graph.flow(internal_hedge) == internal_flow {
                1
            } else {
                -1
            },
        );
    }

    Ok(FeynkitRuntimeMap {
        edges,
        vertices,
        half_edges,
        momentum_signs,
    })
}

fn feynkit_runtime_half_edge(
    half_edge: FeynkitDiagramHalfEdge,
    mapping: &FeynkitRuntimeMap,
) -> Result<Hedge> {
    let flow = match half_edge.endpoint {
        FeynkitDiagramEndpoint::Source => Flow::Source,
        FeynkitDiagramEndpoint::Target => Flow::Sink,
    };
    mapping
        .half_edges
        .get(&(half_edge.edge, flow))
        .copied()
        .ok_or_else(|| {
            eyre!(
                "finalized FeynKit cut references missing {:?} endpoint of edge {}",
                half_edge.endpoint,
                half_edge.edge.0
            )
        })
}

fn feynkit_runtime_cuts(
    diagram: &FeynmanDiagram,
    graph: &NumGraph,
    mapping: &FeynkitRuntimeMap,
) -> Result<Vec<FinalizedCut>> {
    diagram
        .cuts()
        .iter()
        .map(|cut| {
            let mut left = graph.empty_subgraph::<SuBitGraph>();
            for half_edge in &cut.left.half_edges {
                left.add(feynkit_runtime_half_edge(*half_edge, mapping)?);
            }
            let mut right = graph.empty_subgraph::<SuBitGraph>();
            for half_edge in &cut.right.half_edges {
                right.add(feynkit_runtime_half_edge(*half_edge, mapping)?);
            }
            let mut oriented_left = graph.empty_subgraph::<SuBitGraph>();
            for half_edge in &cut.cut {
                oriented_left.add(feynkit_runtime_half_edge(*half_edge, mapping)?);
            }
            let cut = OrientedCut::from_underlying_strict(oriented_left, graph)
                .map_err(|error| eyre!("invalid finalized FeynKit cut orientation: {error}"))?;
            Ok(FinalizedCut { cut, left, right })
        })
        .collect()
}

fn translate_feynkit_atom(atom: &Atom, mapping: &FeynkitRuntimeMap) -> Result<Atom> {
    fn natural(argument: AtomView<'_>, owner: &str) -> Result<usize> {
        let value = i64::try_from(argument)
            .map_err(|_| eyre!("{owner} expects a non-negative integer, found {argument}"))?;
        usize::try_from(value)
            .map_err(|_| eyre!("{owner} expects a non-negative integer, found {value}"))
    }

    let source_index = symbol!("FeynKit::SourceIndex");
    let sink_index = symbol!("FeynKit::SinkIndex");
    let edge_dummy = symbol!("FeynKit::EdgeDummy");
    let vertex_dummy = symbol!("FeynKit::VertexDummy");
    let momentum = symbol!("FeynKit::Momentum");

    let mut error = None;
    let indexed = atom.replace_map(|term, _, out| {
        if error.is_some() {
            return;
        }
        let AtomView::Fun(function) = term else {
            return;
        };
        let symbol = function.get_symbol();
        if ![source_index, sink_index, edge_dummy, vertex_dummy].contains(&symbol) {
            return;
        }
        let arguments = function.iter().collect_vec();
        let translated = (|| -> Result<Atom> {
            let [owner, local] = arguments.as_slice() else {
                return Err(eyre!("{} expects two arguments", term));
            };
            let owner = natural(*owner, "a FeynKit index owner")?;
            let local = u16::try_from(natural(*local, "a FeynKit local index")?)
                .map_err(|_| eyre!("FeynKit local index does not fit in u16"))?;
            if symbol == vertex_dummy {
                let node = mapping
                    .vertices
                    .get(&FeynkitVertexId(owner))
                    .ok_or_else(|| eyre!("unknown FeynKit vertex {owner}"))?;
                return Ok(node.aind(local).into());
            }
            let edge = FeynkitEdgeId(owner);
            if symbol == edge_dummy {
                let edge = mapping
                    .edges
                    .get(&edge)
                    .ok_or_else(|| eyre!("unknown FeynKit edge {owner}"))?;
                return Ok(edge.aind(local).into());
            }
            let flow = if symbol == source_index {
                Flow::Source
            } else {
                Flow::Sink
            };
            let hedge = mapping
                .half_edges
                .get(&(edge, flow))
                .ok_or_else(|| eyre!("unknown {flow:?} half-edge for FeynKit edge {owner}"))?;
            Ok(hedge.aind(local).into())
        })();
        match translated {
            std::result::Result::Ok(replacement) => **out = replacement,
            Err(reason) => error = Some(reason),
        }
    });
    if let Some(error) = error {
        return Err(error);
    }

    error = None;
    let translated = indexed.replace_map(|term, _, out| {
        if error.is_some() {
            return;
        }
        let AtomView::Fun(function) = term else {
            return;
        };
        let function_symbol = function.get_symbol();
        if [PS.u, PS.ubar, PS.v, PS.vbar, PS.eps, PS.ebar].contains(&function_symbol) {
            let arguments = function.iter().collect_vec();
            let replacement = (|| -> Result<Atom> {
                let [edge, index] = arguments.as_slice() else {
                    return Err(eyre!(
                        "a FeynKit external wavefunction expects an edge and an index"
                    ));
                };
                let edge = FeynkitEdgeId(natural(*edge, "a FeynKit external-state edge")?);
                let runtime = mapping
                    .edges
                    .get(&edge)
                    .ok_or_else(|| eyre!("unknown FeynKit external-state edge {}", edge.0))?;
                let head = match function_symbol {
                    symbol if symbol == PS.u => crate::utils::GS.u,
                    symbol if symbol == PS.ubar => crate::utils::GS.ubar,
                    symbol if symbol == PS.v => crate::utils::GS.v,
                    symbol if symbol == PS.vbar => crate::utils::GS.vbar,
                    symbol if symbol == PS.eps => crate::utils::GS.epsilon,
                    symbol if symbol == PS.ebar => crate::utils::GS.epsilonbar,
                    _ => unreachable!("the polarization symbol was checked above"),
                };
                Ok(FunctionBuilder::new(head)
                    .add_arg(
                        i64::try_from(runtime.0).map_err(|_| eyre!("runtime edge is too large"))?,
                    )
                    .add_arg(index.to_owned())
                    .finish())
            })();
            match replacement {
                std::result::Result::Ok(replacement) => **out = replacement,
                Err(reason) => error = Some(reason),
            }
            return;
        }
        if function_symbol != momentum {
            return;
        }
        let arguments = function.iter().collect_vec();
        let replacement = (|| -> Result<Atom> {
            let [edge, index] = arguments.as_slice() else {
                return Err(eyre!("FeynKit::Momentum expects an edge and an index"));
            };
            let edge = FeynkitEdgeId(natural(*edge, "FeynKit::Momentum edge")?);
            let runtime = mapping
                .edges
                .get(&edge)
                .ok_or_else(|| eyre!("unknown FeynKit momentum edge {}", edge.0))?;
            let momentum = FunctionBuilder::new(crate::utils::GS.emr_mom)
                .add_arg(i64::try_from(runtime.0).map_err(|_| eyre!("runtime edge is too large"))?)
                .add_arg(index.to_owned())
                .finish();
            Ok(if mapping.momentum_signs.get(&edge) == Some(&-1) {
                -momentum
            } else {
                momentum
            })
        })();
        match replacement {
            std::result::Result::Ok(replacement) => **out = replacement,
            Err(reason) => error = Some(reason),
        }
    });
    if let Some(error) = error {
        return Err(error);
    }
    Ok(translated)
}

fn add_feynkit_sign(left: SignOrZero, right: SignOrZero) -> Result<SignOrZero> {
    match (left, right) {
        (SignOrZero::Zero, value) | (value, SignOrZero::Zero) => Ok(value),
        (SignOrZero::Plus, SignOrZero::Minus) | (SignOrZero::Minus, SignOrZero::Plus) => {
            Ok(SignOrZero::Zero)
        }
        _ => Err(eyre!(
            "identifying cross-section external momenta produced a coefficient outside -1, 0, 1"
        )),
    }
}

fn feynkit_runtime_lmb(
    diagram: &FeynmanDiagram,
    graph: &NumGraph,
    mapping: &FeynkitRuntimeMap,
    external_connections: &[(Option<usize>, Option<usize>)],
) -> Result<LoopMomentumBasis> {
    let source = diagram.loop_momentum_basis();
    let external_edges = feynkit_external_edges(diagram);

    let mut tree: SuBitGraph = graph.empty_subgraph();
    let mut seen_tree_edges = BTreeSet::new();
    for edge in &source.tree_edges {
        let runtime = *mapping
            .edges
            .get(edge)
            .ok_or_else(|| eyre!("unknown FeynKit tree edge {}", edge.0))?;
        if seen_tree_edges.insert(runtime) {
            tree.add(graph[&runtime].1);
        }
    }

    let loop_edges = source
        .loop_edges
        .iter()
        .map(|edge| {
            mapping
                .edges
                .get(edge)
                .copied()
                .ok_or_else(|| eyre!("unknown FeynKit loop edge {}", edge.0))
        })
        .collect::<Result<Vec<_>>>()?;

    // Runtime cross sections identify each incoming/outgoing pair as one cut
    // edge. Amplitudes retain the one canonical external edge per connection.
    let mut external_groups = Vec::<Vec<FeynkitEdgeId>>::new();
    let mut ext_edges = Vec::new();
    for (connection, (incoming, outgoing)) in external_connections.iter().enumerate() {
        let mut group = Vec::new();
        for tag in [incoming, outgoing].into_iter().flatten() {
            let edge = *external_edges
                .get(&feynkit_tag_index(*tag)?)
                .ok_or_else(|| eyre!("missing FeynKit external edge for tag {tag}"))?;
            if !group.contains(&edge) {
                group.push(edge);
            }
        }
        let runtime = if incoming.is_some() && outgoing.is_some() {
            EdgeIndex::from(connection)
        } else {
            *mapping
                .edges
                .get(group.first().ok_or_else(|| {
                    eyre!("external connection {connection} contains no FeynKit leg")
                })?)
                .ok_or_else(|| eyre!("missing runtime external edge for connection {connection}"))?
        };
        external_groups.push(group);
        ext_edges.push(runtime);
    }

    let source_external_positions = source
        .external_edges
        .iter()
        .enumerate()
        .map(|(index, edge)| (*edge, index))
        .collect::<BTreeMap<_, _>>();
    let mut origins = BTreeMap::<EdgeIndex, FeynkitEdgeId>::new();
    for (source_edge, runtime_edge) in &mapping.edges {
        origins.entry(*runtime_edge).or_insert(*source_edge);
    }
    for (connection, group) in external_groups.iter().enumerate() {
        if let Some(incoming) = group.first() {
            origins.insert(ext_edges[connection], *incoming);
        }
    }

    let mut edge_signatures = Vec::with_capacity(graph.n_edges());
    for runtime_index in 0..graph.n_edges() {
        let runtime = EdgeIndex(runtime_index);
        let origin = origins
            .get(&runtime)
            .ok_or_else(|| eyre!("runtime edge {runtime} has no FeynKit origin"))?;
        let signature = source.edge_signatures.get(origin).ok_or_else(|| {
            eyre!(
                "FeynKit loop-momentum basis has no signature for edge {}",
                origin.0
            )
        })?;
        let negate = mapping.momentum_signs.get(origin) == Some(&-1);
        let internal = SignatureLike::from_iter(
            signature
                .loops
                .iter()
                .map(|sign| if negate { -sign } else { sign }),
        );
        let mut external = Vec::with_capacity(external_groups.len());
        for group in &external_groups {
            let mut coefficient = SignOrZero::Zero;
            for edge in group {
                let Some(position) = source_external_positions.get(edge) else {
                    continue;
                };
                let sign = signature.external.get(*position).ok_or_else(|| {
                    eyre!("FeynKit external signature is shorter than its external-edge basis")
                })?;
                coefficient = add_feynkit_sign(coefficient, sign)?;
            }
            external.push(if negate { -coefficient } else { coefficient });
        }
        edge_signatures.push(LoopExtSignature {
            internal,
            external: SignatureLike::from_iter(external),
        });
    }

    Ok(LoopMomentumBasis {
        tree,
        loop_edges: loop_edges.into(),
        ext_edges: ext_edges.into(),
        edge_signatures: edge_signatures.into(),
    })
}

fn display_graph_source_path(path: &Path) -> PathBuf {
    if path.is_absolute() {
        return path.canonicalize().unwrap_or_else(|_| path.to_path_buf());
    }

    let joined = std::env::current_dir()
        .map(|cwd| cwd.join(path))
        .unwrap_or_else(|_| path.to_path_buf());
    joined.canonicalize().unwrap_or(joined)
}

impl Graph {
    pub fn dot_serialize(&self, settings: &DotExportSettings) -> String {
        let mut out = String::new();
        self.dot_serialize_fmt(&mut out, settings).unwrap();
        out
    }

    pub(crate) fn dot_serialize_io(
        &self,
        writer: &mut impl std::io::Write,
        settings: &DotExportSettings,
    ) -> Result<(), std::io::Error> {
        let g = self.to_dot_graph_with_settings(settings);
        g.write_io(writer)
    }

    #[allow(dead_code)]
    pub(crate) fn dot_split_serialize_io(
        &self,
        writer: &mut impl std::io::Write,
    ) -> Result<(), std::io::Error> {
        let g = self.to_split_dotgraph();
        g.write_io(writer)
    }

    pub fn dot_serialize_fmt(
        &self,
        writer: &mut impl std::fmt::Write,
        settings: &DotExportSettings,
    ) -> Result<(), std::fmt::Error> {
        let g = self.to_dot_graph_with_settings(settings);
        g.write_fmt(writer)
    }

    pub(crate) fn from_parsed_with_validation(graph: ParseGraph, model: &Model) -> Result<Self> {
        let res = Self::from_parsed(graph, model)?;
        res.validate_full_numerator_tensor_network()
            .with_context(|| {
                format!(
                    "Failed to validate full numerator tensor network for graph {}",
                    res.name
                )
            })?;
        Ok(res)
    }

    /// Enrich a finalized FeynKit diagram with GammaLoop runtime caches.
    ///
    /// This is deliberately a mechanical downward conversion: it copies the
    /// canonical topology, symbolic fragments, factors, projector, and routing.
    /// It never canonicalizes again, resolves a vertex rule, generates a
    /// numerator, or chooses a second loop-momentum basis.
    pub(crate) fn from_feynkit(
        diagram: &FeynmanDiagram,
        group_id: Option<GroupId>,
        is_group_master: bool,
    ) -> Result<Self> {
        diagram.validate().with_context(|| {
            format!(
                "FeynKit diagram '{}' is not finalized consistently",
                diagram.name()
            )
        })?;

        let model = diagram.model();
        let external_connections = feynkit_external_connections(diagram)?;

        let mut parsed = ParseGraph::from_feynkit_diagram(diagram, model, &external_connections)?;
        parsed.global_data.group_id = group_id;
        parsed.global_data.is_group_master = is_group_master;
        let (mut initial_data, mut graph) = Self::extract_initial_data(&parsed, model)?;

        graph
            .sew(
                |_, left, _, right| {
                    matches!((left.data.is_cut, right.data.is_cut), (Some(a), Some(b)) if a == b)
                },
                |left_flow, left, right_flow, right| match (left_flow, right_flow) {
                    (Flow::Sink, Flow::Source) => (Flow::Sink, left),
                    (Flow::Source, Flow::Sink) => (Flow::Source, right),
                    _ => panic!(
                        "cannot sew FeynKit cut hedges with flows {left_flow:?} and {right_flow:?}"
                    ),
                },
            )
            .map_err(|error| eyre!("FeynKit graph sewing failed: {error:?}"))?;
        let mut cut_result = Self::process_cut_edges(&graph)?;
        cut_result.permute(&mut graph)?;

        let mapping = feynkit_runtime_map(diagram, &graph, &external_connections)?;
        for (_, _, vertex) in graph.iter_nodes_mut() {
            let numerator = vertex.num.as_ref().ok_or_else(|| {
                eyre!("FeynKit runtime vertex is missing its finalized numerator")
            })?;
            vertex.num = Some(translate_feynkit_atom(numerator, &mapping)?);
        }
        let edge_ids = (0..graph.n_edges()).map(EdgeIndex::from).collect_vec();
        for edge in edge_ids {
            let numerator = graph[edge]
                .num
                .as_ref()
                .ok_or_else(|| eyre!("FeynKit runtime edge {edge} is missing its numerator"))?;
            graph[edge].num = Some(translate_feynkit_atom(numerator, &mapping)?);
        }
        initial_data.overall_factor =
            translate_feynkit_atom(&initial_data.overall_factor, &mapping)?;
        initial_data.global_prefactor.num =
            translate_feynkit_atom(&initial_data.global_prefactor.num, &mapping)?;
        initial_data.global_prefactor.projector =
            translate_feynkit_atom(&initial_data.global_prefactor.projector, &mapping)?;
        let translated_numerator = translate_feynkit_atom(diagram.numerator(), &mapping)?;
        let initial_state_cut =
            OrientedCut::from_underlying_strict(cut_result.initial_hedges, &graph)?;
        let mut finalized_cuts = feynkit_runtime_cuts(diagram, &graph, &mapping)?;
        for cut in &mut finalized_cuts {
            // Before sewing, FeynKit cut sides contain ordinary external legs.
            // Sewing identifies each incoming/outgoing pair and would make
            // those two legs look like an extra loop on one amplitude side.
            // GammaLoop represents the same momenta as the initial-state
            // external basis, so discard the sewn partner half-edges from the
            // side subgraphs while retaining the physical oriented cut.
            cut.left.subtract_with(&initial_state_cut.right);
            cut.right.subtract_with(&initial_state_cut.right);
        }

        let loop_momentum_basis =
            feynkit_runtime_lmb(diagram, &graph, &mapping, &external_connections)?;
        let global_prefactor = initial_data.global_prefactor;
        let polarizations = global_prefactor.polarizations();
        let param_builder = ParamBuilder::new(
            &(&polarizations, &graph),
            model,
            &loop_momentum_basis,
            initial_data.additional_params.clone(),
        );

        let underlying =
            Self::build_underlying_graph(graph, model, &param_builder).with_context(|| {
                format!(
                    "failed to build GammaLoop runtime storage for finalized FeynKit diagram {}",
                    initial_data.name
                )
            })?;

        let mut full_without_initials = underlying.full_filter();
        full_without_initials.subtract_with(&initial_state_cut.left);
        let mut tree_edges = underlying.bridges_of(&full_without_initials);
        tree_edges.union_with(&initial_state_cut.left);
        let mut result = Graph {
            overall_factor: initial_data.overall_factor,
            polarizations,
            global_prefactor,
            tree_edges,
            name: initial_data.name,
            loop_momentum_basis,
            initial_state_cut,
            underlying,
            surface_cache: SurfaceCache::default(),
            group_id: initial_data.group_id,
            is_group_master: initial_data.is_group_master,
            param_builder,
            finalized_cuts,
        };
        result.param_builder = ParamBuilder::new(
            &result,
            model,
            &result.loop_momentum_basis,
            initial_data.additional_params,
        );
        let runtime_numerator = result
            .numerator(&result.full_filter(), &result.empty_subgraph())
            .get_single_atom()?;
        if runtime_numerator.expand() != translated_numerator.expand() {
            return Err(eyre!(
                "GammaLoop runtime conversion of FeynKit diagram '{}' did not preserve its finalized numerator",
                diagram.name()
            ));
        }
        Ok(result)
    }

    #[instrument(skip_all, fields(graph= %graph.debug_dot(),name = %graph.global_data.name.as_str()))]
    pub(crate) fn from_parsed(graph: ParseGraph, model: &Model) -> Result<Self> {
        if graph.global_data.projectors.is_none() {
            return Err(eyre!(
                "finalized DOT graph '{}' must provide an explicit projector (use `1` when no projector is required)",
                graph.global_data.name
            ));
        }
        for (vertex, _, data) in graph.graph.iter_nodes() {
            if data.num.is_none() {
                return Err(eyre!(
                    "finalized DOT graph '{}' is missing the numerator for vertex {vertex}",
                    graph.global_data.name
                ));
            }
        }
        for (_, edge, data) in graph.graph.iter_edges() {
            if data.data.num.is_none() {
                return Err(eyre!(
                    "finalized DOT graph '{}' is missing the numerator for edge {edge}",
                    graph.global_data.name
                ));
            }
            if data.data.is_cut.is_some() {
                return Err(eyre!(
                    "finalized runtime DOT graph '{}' contains cross-section sewing metadata but no canonical physical cuts; import the canonical FeynmanDiagram DOT artifact through FeynKit instead",
                    graph.global_data.name
                ));
            }
        }

        let (initial_data, mut graph) = Self::extract_initial_data(&graph, model)?;

        // Sew the graph based on cut edges
        graph
            .sew(
                |_, ae, _, be| {
                    if let (Some(a), Some(b)) = (ae.data.is_cut, be.data.is_cut) {
                        a == b
                    } else {
                        false
                    }
                },
                |af, ae, bf, be| match (af, bf) {
                    (Flow::Sink, Flow::Source) => (Flow::Sink, ae),
                    (Flow::Source, Flow::Sink) => (Flow::Source, be),
                    _ => panic!("Cannot sew hedges with flow {:?} and {:?}", af, bf),
                },
            )
            .map_err(|e| eyre::eyre!("Graph sewing failed: {:?}", e))?;

        let mut cut_result = Self::process_cut_edges(&graph)?;

        cut_result.permute(&mut graph)?;

        let initial_state_cut =
            OrientedCut::from_underlying_strict(cut_result.initial_hedges, &graph)?;

        debug!("Initial state cut: {}", graph.dot(&initial_state_cut.left));
        debug_assert!(!initial_data.add_polarizations);
        let global_prefactor = initial_data.global_prefactor;
        let polarizations = global_prefactor.polarizations();
        let loop_momentum_basis = Self::materialize_explicit_loop_momentum_basis(
            &graph,
            &cut_result.full_cut,
            &cut_result.lmb_ids,
            &cut_result.xs_ext_id,
        )
        .with_context(|| format!("Failed to build lmb for graph {}", initial_data.name))?;
        let param_builder = ParamBuilder::new(
            &(&polarizations, &graph),
            model,
            &loop_momentum_basis,
            initial_data.additional_params.clone(),
        );

        let underlying = Self::build_underlying_graph(graph, model, &param_builder)
            .with_context(|| format!("Failed to build underlying graph {}", initial_data.name))?;

        let mut full_without_initials = underlying.full_filter();
        full_without_initials.subtract_with(&initial_state_cut.left);
        let mut tree_edges = underlying.bridges_of(&full_without_initials);
        tree_edges.union_with(&initial_state_cut.left);

        let mut g = Graph {
            overall_factor: initial_data.overall_factor,
            polarizations: global_prefactor.polarizations(),
            global_prefactor,
            tree_edges,
            name: initial_data.name,
            loop_momentum_basis,
            initial_state_cut,
            underlying,
            surface_cache: SurfaceCache::default(),
            group_id: initial_data.group_id,
            is_group_master: initial_data.is_group_master,
            param_builder,
            finalized_cuts: Vec::new(),
        };

        let external_momentum_edge_order = g.external_momentum_edge_order();
        g.loop_momentum_basis
            .canonicalize_external_order(&external_momentum_edge_order);

        let updated_param_builder_with_lmb = ParamBuilder::new(
            &g,
            model,
            &g.loop_momentum_basis,
            initial_data.additional_params,
        );

        debug!(
            "Updated param builder with LMB: {}\n{}",
            g.loop_momentum_basis,
            updated_param_builder_with_lmb.table(),
        );

        g.param_builder = updated_param_builder_with_lmb;

        debug!("{}", g.debug_dot());

        Ok(g)
    }

    fn validate_full_numerator_tensor_network(&self) -> Result<()> {
        let full_num = self
            .numerator(&self.full_filter(), &self.empty_subgraph())
            .get_single_atom()
            .unwrap()
            * &self.global_prefactor.num
            * &self.global_prefactor.projector
            * &self.overall_factor;
        let color_simplified = full_num
            .as_view()
            .simplify_color_with(ColorSimplifySettings::default().with_cof_dimension_invariants());
        if !full_num.is_zero() && color_simplified.is_zero() {
            warn!(
                "Full numerator for graph '{}' becomes zero after color algebra. The graph/projector color structure likely annihilates the amplitude.",
                self.name
            );
        }
        let net = full_num
            .parse_to_symbolic_net::<Aind>(&ParseSettings::default())
            .map_err(Report::from)?;
        let dangling = net.graph.dangling_indices();
        if !dangling.is_empty() {
            return Err(eyre!(
                "Full numerator still has dangling tensor indices: \n{}",
                dangling
                    .iter()
                    .map(|slot| format!(
                        "{}:{}",
                        slot.to_atom().log_print(None),
                        slot.to_atom().to_plain_string()
                    ))
                    .join(",\n")
            ));
        }

        Ok(())
    }

    fn extract_initial_data(
        parse_graph: &ParseGraph,
        model: &Model,
    ) -> Result<(InitialGraphData, NumGraph)> {
        let hedge_order = parse_graph.hedge_order()?;
        let global_data = &parse_graph.global_data;

        let initial_data = InitialGraphData {
            additional_params: global_data.parameters.clone(),
            overall_factor: global_data.overall_factor.clone(),
            global_prefactor: GlobalPrefactor {
                num: global_data.num.clone(),
                projector: global_data.projectors.clone().unwrap_or(Atom::one()),
            },
            add_polarizations: global_data.projectors.is_none(),
            group_id: global_data.group_id,
            is_group_master: global_data.is_group_master,
            name: global_data.name.clone(),
        };

        let num_graph = parse_graph.graph.map_data_ref(
            |_, _, v| v.clone(),
            |_, _, _, e| e.map(|e| e.clone()),
            |h, hd| HedgeData {
                num_indices: NumIndices::parse(parse_graph, model)(h, hd),
                ufo_order: Autogen::explicit(hedge_order[h.0]),
            },
        );

        Ok((initial_data, num_graph))
    }

    fn process_cut_edges(graph: &NumGraph) -> Result<CutProcessingResult> {
        let mut lmb_ids: BTreeMap<LoopIndex, EdgeIndex> = BTreeMap::new();
        let mut xs_ext_id: BTreeMap<Hedge, (EdgeIndex, Hedge)> = BTreeMap::new();
        let mut full_cut: SuBitGraph = graph.full_filter();

        for (p, eid, e) in graph.iter_edges() {
            let HedgePair::Paired { sink, .. } = p else {
                if e.data.is_cut.is_some() {
                    //As we have already sewn the graph, all cut edges must be paired, failure to do so would indicate a bug
                    return Err(eyre!("Cut edge must be paired"));
                } else {
                    continue;
                }
            };

            if let Some(lmb_id) = e.data.lmb_id {
                if let Some(old_value) = lmb_ids.insert(lmb_id, eid) {
                    return Err(eyre!(
                        "lmb_id {lmb_id:?} already exists with value {old_value:?}",
                    ));
                }
                debug!("Cutting {eid} for lmb_id{lmb_id}");
                full_cut.sub(p);
            } else if let Some(h) = e.data.is_cut {
                if let Some(old_value) = xs_ext_id.insert(h, (eid, sink)) {
                    return Err(eyre!("h {h:?} already exists with value {old_value:?}",));
                }
                full_cut.sub(p);
            }
        }

        // debug!("Graph now:{}", graph.dot(full_cut));

        let mut initial_hedges: SuBitGraph = graph.empty_subgraph();
        for (target_pos, _) in xs_ext_id.iter().enumerate() {
            initial_hedges.add(Hedge(target_pos));
        }

        Ok(CutProcessingResult {
            full_cut,
            lmb_ids,
            xs_ext_id,
            initial_hedges,
        })
    }

    fn build_underlying_graph(
        graph: NumGraph,
        model: &Model,
        param_builder: &ParamBuilder,
    ) -> Result<UnderlyingGraph> {
        let intermediate: UnderlyingGraph = graph.map_result(
            |_, i, v| {
                let num = Autogen::explicit(v.num.ok_or_else(|| {
                    eyre!("finalized graph is missing the numerator for vertex {i}")
                })?);

                let dod = match v.dod {
                    Some(dod) => Autogen::explicit(dod),
                    None => Autogen::generated(num.all_dod()),
                };

                Ok(Vertex {
                    name: Autogen::from_option_or_generate(v.name, || i.to_string()),
                    num,
                    dod,
                    vertex_rule: v.vertex_rule,
                })
            },
            |_, _, _, eid, ed| {
                let e = ed.data;
                if e.particle.is_fermion(model)
                    && !e.particle.is_self_antiparticle(model)
                    && e.particle.orientation(model) != ed.orientation
                {
                    return Err(eyre!(
                        "Edge orientation {:?} does not match particle orientation {:?} for edge {},{}",
                        ed.orientation,
                        e.particle.orientation(model),
                        eid,
                        e
                    ));
                }

                let mass = EdgeMass::from_atom(e.particle.mass_atom(model), model, param_builder)?;

                let num = Autogen::explicit(e.num.ok_or_else(|| {
                    eyre!("finalized graph is missing the numerator for edge {eid}")
                })?);

                let dod = match e.dod {
                    Some(dod) => Autogen::explicit(dod),
                    None => Autogen::generated(num.edge_dod(eid) - 2),
                };

                Ok(EdgeData::new(
                    Edge {
                        mass,
                        is_dummy: e.is_dummy,
                        name: Autogen::from_option_or_generate(e.name, || eid.to_string()),
                        particle: e.particle,
                        num,
                        dod,
                        extra_data: EdgeExtraData {
                            momtrop_edge_power: e.momtrop_edge_power,
                            vakint_edge_power: e.vakint_edge_power,
                        }
                    },
                    ed.orientation,
                ))
            },
            |_, h| Ok(h),
        )?;

        Ok(intermediate)
    }

    /// Materialize signatures from the exact loop edges selected by `lmb_id`.
    ///
    /// The spanning-forest routine only propagates momenta through that fixed
    /// complement. Missing, extra, or substituted loop edges are rejected, so
    /// the DOT runtime import cannot silently choose a second basis.
    fn materialize_explicit_loop_momentum_basis(
        graph: &NumGraph,
        full_cut: &SuBitGraph,
        lmb_ids: &BTreeMap<LoopIndex, EdgeIndex>,
        xs_ext_id: &BTreeMap<Hedge, (EdgeIndex, Hedge)>,
    ) -> Result<LoopMomentumBasis> {
        debug!("{}", graph.dot(full_cut));

        let mut full = graph.full_filter();
        for (pair, _, edge) in graph.iter_edges() {
            if edge.data.is_dummy {
                full.sub(pair);
            }
        }
        let total_loops = graph.cyclotomatic_number(&full);
        let explicit_loops = total_loops.checked_sub(xs_ext_id.len()).ok_or_else(|| {
            eyre!(
                "graph has {total_loops} loops but {} cut edges were marked as external",
                xs_ext_id.len()
            )
        })?;
        let expected_ids = (0..explicit_loops).map(LoopIndex).collect_vec();
        let actual_ids = lmb_ids.keys().copied().collect_vec();
        if actual_ids != expected_ids {
            return Err(eyre!(
                "finalized DOT graph must label exactly {explicit_loops} loop edges with contiguous lmb_id values 0..{explicit_loops}; found {actual_ids:?}"
            ));
        }

        let external = graph.internal_crown(&full);
        let mut loop_momentum_basis = graph.lmb_impl(&full, full_cut, external)?;

        for e in 0..xs_ext_id.len() {
            let (l, _) = loop_momentum_basis
                .loop_edges
                .iter()
                .find_position(|a| *a == &EdgeIndex(e))
                .ok_or_else(|| {
                    eyre!("cut edge {e} is not a loop edge in the explicit momentum basis")
                })?;

            loop_momentum_basis.put_loop_to_ext(LoopIndex(l));
        }

        let materialized_edges = loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .collect::<BTreeSet<_>>();
        let selected_edges = lmb_ids.values().copied().collect::<BTreeSet<_>>();
        if materialized_edges != selected_edges {
            return Err(eyre!(
                "lmb_id edges {selected_edges:?} do not form a loop-momentum basis; materialized edges were {materialized_edges:?}"
            ));
        }

        // Put the explicitly labelled edges in their requested order.
        for (target, edge) in lmb_ids {
            let current = loop_momentum_basis
                .loop_edges
                .iter()
                .position(|candidate| candidate == edge)
                .ok_or_else(|| {
                    eyre!(
                        "explicit loop edge {edge} disappeared while ordering the finalized basis"
                    )
                })?;
            if current != target.0 {
                loop_momentum_basis.swap_loops(LoopIndex(current), *target);
            }
        }

        Ok(loop_momentum_basis)
    }

    /// Import a fully finalized amplitude runtime artifact.
    ///
    /// The artifact must contain explicit numerator fragments, projector, UFO
    /// half-edge slots, and loop-momentum-basis IDs. Cross-section DOT must be
    /// imported as a canonical [`FeynmanDiagram`] so its typed cuts survive.
    pub fn from_finalized_runtime_dot(graph: DotGraph, model: &Model) -> Result<Self> {
        Self::from_parsed(ParseGraph::from_parsed(graph, model)?, model)
    }

    /// Import finalized amplitude runtime artifacts from one DOT file.
    pub fn from_finalized_runtime_file<P>(p: P, model: &Model) -> Result<Vec<Self>>
    where
        P: AsRef<Path>,
    {
        Self::from_finalized_runtime_path(p, model)
    }

    /// Import finalized amplitude runtime artifacts from a file or directory.
    pub fn from_finalized_runtime_path<P>(p: P, model: &Model) -> Result<Vec<Self>>
    where
        P: AsRef<Path>,
    {
        let path = p.as_ref();

        if path.is_dir() {
            // Load all .dot files from directory
            let mut all_graphs = Vec::new();
            let entries = std::fs::read_dir(path)
                .with_context(|| format!("Failed to read directory: {}", path.display()))?;

            let mut dot_files = Vec::new();
            for entry in entries {
                let entry = entry?;
                let file_path = entry.path();
                if file_path.is_file() && file_path.extension().is_some_and(|ext| ext == "dot") {
                    dot_files.push(file_path);
                }
            }

            // Sort files for consistent ordering
            dot_files.sort();

            for dot_file in dot_files {
                let graphs = Self::from_single_finalized_runtime_file(&dot_file, model)?;
                all_graphs.extend(graphs);
            }

            if all_graphs.is_empty() {
                return Err(eyre!(
                    "No .dot files found in directory: {}",
                    path.display()
                ));
            }

            Ok(all_graphs)
        } else {
            // Load single file
            Self::from_single_finalized_runtime_file(path, model)
        }
    }

    fn from_single_finalized_runtime_file<P>(p: P, model: &Model) -> Result<Vec<Self>>
    where
        P: AsRef<Path>,
    {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_file(p.as_ref()).map_err(|a| match a {
            HedgeParseError::GraphFromFile(e) => match e.as_ref() {
                dot_parser::ast::GraphFromFileError::FileError(e) => eyre!(e.to_string())
                    .with_note(|| {
                        format!(
                            "Tried to access the file at: {}",
                            display_graph_source_path(p.as_ref()).display()
                        )
                    }),
                dot_parser::ast::GraphFromFileError::ParseError(e) => {
                    eyre!("Dot parsing error: {}", e)
                }
                dot_parser::ast::GraphFromFileError::PestParseError(e) => {
                    eyre!(e.to_string())
                }
            },
            HedgeParseError::ParseError(i) => color_eyre::Report::from(i),
            _ => {
                eyre!("Hedge parse error")
            }
        })?;
        Self::from_finalized_runtime_graph_set(hedge_graph_set, model)
    }

    /// Import finalized amplitude runtime artifacts from a DOT string.
    pub fn from_finalized_runtime_string<Str: AsRef<str>>(
        s: Str,
        model: &Model,
    ) -> Result<Vec<Self>> {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_string(s).map_err(|a| match a {
            HedgeParseError::GraphFromFile(e) => match e.as_ref() {
                dot_parser::ast::GraphFromFileError::FileError(e) => {
                    eyre!(e.to_string())
                }
                dot_parser::ast::GraphFromFileError::ParseError(e) => {
                    eyre!("Dot parsing error: {}", e)
                }
                dot_parser::ast::GraphFromFileError::PestParseError(e) => {
                    eyre!(e.to_string())
                }
            },
            HedgeParseError::ParseError(i) => color_eyre::Report::from(i),
            _ => {
                eyre!("Hedge parse error")
            }
        })?;

        Self::from_finalized_runtime_graph_set(hedge_graph_set, model)
    }

    fn from_finalized_runtime_graph_set(
        set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        >,
        model: &Model,
    ) -> Result<Vec<Self>> {
        let mut graphs = Vec::new();

        for (graph, global_data) in set.set.into_iter().zip(set.global_data) {
            let graph = DotGraph { global_data, graph };
            debug!("Parsing: \n{}", graph.debug_dot());
            graphs.push(Graph::from_parsed(
                ParseGraph::from_parsed(graph, model)?,
                model,
            )?);
        }
        Ok(graphs)
    }
}

pub mod serialization;

/// completes and extract the user defined group structure on a lis of graphs
pub(crate) fn complete_group_parsing(graphs: &mut [Graph]) -> Result<TiVec<GroupId, GraphGroup>> {
    // validate the input
    let defined_group_ids = graphs
        .iter()
        .filter_map(|graph| graph.group_id)
        .sorted()
        .dedup()
        .collect_vec();

    let expected_group_ids = (0..defined_group_ids.len()).map(GroupId).collect_vec();

    if defined_group_ids != expected_group_ids {
        return Err(eyre!(
            "invalid group ids, group ids must start at 0 and contain no gaps"
        ));
    }
    // now set the remaining group ids
    let mut current_group_id = defined_group_ids.len();
    for graph in graphs.iter_mut() {
        if graph.group_id.is_none() {
            graph.group_id = Some(GroupId(current_group_id));
            graph.is_group_master = true;
            current_group_id += 1;
        }
    }

    let num_groups = current_group_id;

    // build the groups
    (0..num_groups)
        .map(|group_id| {
            let group_id = GroupId(group_id);
            let graphs_in_group = graphs
                .iter()
                .enumerate()
                .filter_map(|(i, g)| {
                    if g.group_id == Some(group_id) {
                        Some(i)
                    } else {
                        None
                    }
                })
                .collect_vec();

            // the special case of a single graph in the group is easy
            if graphs_in_group.len() == 1 {
                graphs[graphs_in_group[0]].is_group_master = true;
                Ok(GraphGroup {
                    master: graphs_in_group[0],
                    remaining: vec![],
                })
            } else {
                // see if a master is defined
                let master = graphs_in_group
                    .iter()
                    .find(|&&i| graphs[i].is_group_master)
                    .copied();

                if let Some(master) = master {
                    // find the remaining graphs and make sure no other master is defined
                    let remaining = graphs_in_group
                        .into_iter()
                        .filter(|&i| i != master)
                        .collect_vec();

                    let duplicate_master = remaining.iter().any(|&i| graphs[i].is_group_master);

                    if duplicate_master {
                        return Err(eyre!(
                            "Multiple group masters defined for group {group_id:?}"
                        ));
                    }
                    Ok(GraphGroup { master, remaining })
                } else {
                    // no master defined, take the first graph as master
                    let master = graphs_in_group[0];
                    graphs[master].is_group_master = true;
                    Ok(GraphGroup {
                        master,
                        remaining: graphs_in_group[1..].to_vec(),
                    })
                }
            }
        })
        .collect::<Result<TiVec<GroupId, GraphGroup>>>()
}

pub mod from_dot;
pub use from_dot::*;
#[cfg(test)]
pub mod tests;
