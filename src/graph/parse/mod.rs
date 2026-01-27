use std::{collections::BTreeMap, ops::Deref, path::Path};

use crate::{
    feyngen::{
        GenerationType,
        diagram_generator::{EdgeColor, NodeColorWithVertexRule},
    },
    gammaloop_integrand::ParamBuilder,
    graph::{GraphGroup, GroupId, LoopMomentumBasis, edge::EdgeExtraData},
    model::Model,
    momentum_sample::LoopIndex,
    numerator::{
        GlobalPrefactor,
        aind::{Aind, NewAind},
        graph::GeneratePolarizations,
        ufo::UFO,
    },
    processes::DotExportSettings,
    utils::symbolica_ext::DOD,
};
use ahash::{AHashMap, AHashSet};

use color_eyre::{Result, Section};

use eyre::{Context, Ok, eyre};
use itertools::Itertools;
// use symbolica::{atom::Atom, graph::Graph as SymbolicaGraph};

use linnet::{
    half_edge::{
        HedgeGraph, NodeIndex,
        builder::HedgeGraphBuilder,
        involution::{EdgeData, EdgeIndex, EdgeVec, Flow, Hedge, HedgePair},
        nodestore::NodeStorageVec,
        subgraph::{Inclusion, ModifySubSet, OrientedCut, SuBitGraph, SubSetLike, SubSetOps},
        swap::Swap,
    },
    parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GraphSet, HedgeParseError},
    permutation::Permutation,
};
use spenso::{
    contraction::Contract,
    network::library::TensorLibraryData,
    structure::{HasStructure, OrderedStructure, representation::Euclidean},
    tensors::{data::StorageTensor, parametric::ParamTensor},
};
use symbolica::{
    atom::{Atom, AtomOrView},
    graph::Graph as SymbolicaGraph,
};
use tracing::instrument;
use tracing::{debug, info};
use typed_index_collections::TiVec;

use super::{
    Edge, Graph, LMBext, NumHedgeData, Vertex,
    edge::{EdgeMass, ParseEdge},
    global::ParseData,
    hedge_data::{NumIndices, ParseHedge},
    vertex::ParseVertex,
};

/// Extract oriented particles from hedges, filtering out dummy edges
pub fn extract_oriented_particles_from_vertex_hedges<I, V>(
    graph: &HedgeGraph<ParseEdge, V, ParseHedge>,
    hedges: I,
    model: &Model,
) -> Vec<crate::model::ArcParticle>
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
                particle.get_anti_particle(model)
            } else {
                particle
            })
        })
        .collect()
}

// Type aliases for cleaner code
type NumGraph = HedgeGraph<ParseEdge, ParseVertex, NumHedgeData>;
type UnderlyingGraph = HedgeGraph<Edge, Vertex, NumHedgeData>;

pub mod string_utils;
pub use string_utils::{FromStripedStr, StripParse, ToQuoted};

#[derive(Clone, Debug)]
pub struct ParseGraph {
    pub global_data: ParseData,
    pub graph: HedgeGraph<ParseEdge, ParseVertex, ParseHedge>,
}

impl Deref for ParseGraph {
    type Target = HedgeGraph<ParseEdge, ParseVertex, ParseHedge>;

    fn deref(&self) -> &Self::Target {
        &self.graph
    }
}

impl ParseGraph {
    pub fn n_fermion_loops(&self) -> usize {
        let fermions: SuBitGraph = self.graph.from_filter(|a| a.particle.is_fermion());

        self.graph.cyclotomatic_number(&fermions)
    }
    pub fn n_external_fermion_loops(&mut self) -> Result<usize> {
        let internal = self.n_fermion_loops();
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

        Ok(self.n_fermion_loops() - internal)
    }

    pub fn debug_dot(&self) -> String {
        DotGraph::from(self).debug_dot()
    }
    pub(crate) fn hedge_order(&self, model: &Model) -> Result<Vec<u8>> {
        let mut hedges = vec![None; self.n_hedges()];

        for (_, neighs, v) in self.iter_nodes() {
            self.process_vertex_hedges(&mut hedges, neighs, v, model)?;
        }

        hedges
            .into_iter()
            .collect::<Option<Vec<u8>>>()
            .ok_or_else(|| eyre!("Nodes do not cover hedges"))
    }

    fn process_vertex_hedges(
        &self,
        hedges: &mut [Option<u8>],
        neighs: impl Iterator<Item = Hedge>,
        vertex: &ParseVertex,
        model: &Model,
    ) -> Result<()> {
        let (mut particles, vertex_name) = Self::extract_vertex_particles(vertex);
        let mut other_order = particles.len();

        let hedge_vec: Vec<_> = neighs.collect();
        let oriented_particles =
            extract_oriented_particles_from_vertex_hedges(self, hedge_vec.iter().copied(), model);

        // Create iterator that pairs each hedge with its oriented particle (if any)
        let mut particle_iter = oriented_particles.into_iter();

        for h in hedge_vec {
            let eid = self[&h];

            let order = if self[eid].is_dummy || self[eid].particle.particle().is_none() {
                other_order += 1;
                (other_order - 1) as u8
            } else {
                // This hedge has a particle, so get the next oriented particle
                let oriented_particle = particle_iter
                    .next()
                    .expect("Mismatch between hedges and oriented particles");
                debug!("Oriented particle: {h} : {}", oriented_particle.name);
                // Try to match with vertex rule particles
                if let Some(name) = &vertex_name {
                    if let Some((pos, matched_particle)) = particles
                        .iter_mut()
                        .find_position(|p| **p == Some(oriented_particle.clone()))
                    {
                        *matched_particle = None;
                        pos as u8
                    } else {
                        return Err(eyre!(
                            "Particle {} not in vertex rule {}",
                            oriented_particle.name.to_string(),
                            name
                        ));
                    }
                } else {
                    other_order += 1;
                    (other_order - 1) as u8
                }
            };

            hedges[h.0] = Some(order);
        }

        // Verify all particles in vertex rule were matched
        if particles.iter().any(|p| p.is_some()) {
            return Err(eyre!(
                "Particles to vertex rules no match for set: {:?}",
                particles
            ));
        }

        Ok(())
    }

    fn extract_vertex_particles(
        vertex: &ParseVertex,
    ) -> (Vec<Option<crate::model::ArcParticle>>, Option<String>) {
        if let Some(vertex_rule) = &vertex.vertex_rule {
            let particles = vertex_rule
                .particles
                .iter()
                .map(|p| Some(p.clone()))
                .collect();
            (particles, Some(vertex_rule.name.to_string()))
        } else {
            (Vec::new(), None)
        }
    }

    #[instrument(skip_all, fields(graph= %graph.to_dot(),name = %graph_name.as_ref(),external_connections = ?external_connections))]
    pub(crate) fn from_symbolica_graph(
        model: &Model,
        graph_name: impl AsRef<str>,
        graph: &SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
        symmetry_factor: Atom,
        external_connections: &[(Option<usize>, Option<usize>)],
    ) -> Result<Self> {
        fn mark_edge_as_seen(seen_edges: &mut AHashSet<usize>, edge_idx: usize) -> Result<()> {
            if !seen_edges.insert(edge_idx) {
                return Err(eyre!(
                    "External connections must be unique: edge {} already used",
                    edge_idx
                ));
            }
            Ok(())
        }

        fn validate_edge_compatibility(
            out_edge: &symbolica::graph::Edge<EdgeColor>,
            in_edge: &symbolica::graph::Edge<EdgeColor>,
            in_id: usize,
            out_id: usize,
        ) -> Result<()> {
            if out_edge.directed != in_edge.directed {
                return Err(eyre!(
                    "External edges must have the same directedness, for edge ids {} and {} found {:?} and {:?}",
                    in_id,
                    out_id,
                    in_edge,
                    out_edge
                ));
            }

            if in_edge.data.pdg.abs() != out_edge.data.pdg.abs() {
                return Err(eyre!(
                    "External edges must have the same pdg in abs, for edge ids {} and {} found {:?} and {:?}",
                    in_id,
                    out_id,
                    in_edge,
                    out_edge
                ));
            }

            Ok(())
        }

        /// Add dangling hedges based on external connections
        /// external connections is provided in  the process definition order
        /// it maps the external_tag s of the external degree 1 nodes together
        #[allow(clippy::too_many_arguments)]
        fn process_single_connection_internal(
            edge_idx: usize,
            flow: Flow,
            hedge: Option<Hedge>,
            vertex_map: &AHashMap<usize, NodeIndex>,
            seen_edges: &mut AHashSet<usize>,
            graph: &SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
            model: &Model,
            builder: &mut HedgeGraphBuilder<ParseEdge, ParseVertex, ParseHedge>,
        ) -> Result<()> {
            // Only mark as seen if not already processed (for bidirectional case)
            if !seen_edges.contains(&edge_idx) {
                mark_edge_as_seen(seen_edges, edge_idx)?;
            }

            let edge = &graph.edges()[edge_idx];
            let mut data = ParseEdge::from_symbolica_edge(model, &edge.data, hedge);

            let (out_vertex, in_vertex) = edge.vertices;
            let mut orientation = data.particle.orientation();

            // Determine which vertex to connect to and adjust particle/orientation if needed
            let (node_idx, final_orientation, final_flow) = match flow {
                Flow::Source => {
                    if let Some(&sink_node) = vertex_map.get(&in_vertex) {
                        data.particle = data.particle.reverse(model);
                        orientation = orientation.reverse();
                        (sink_node, orientation, flow)
                    } else if let Some(&source_node) = vertex_map.get(&out_vertex) {
                        (source_node, orientation, flow)
                    } else {
                        return Err(eyre!(
                            "Outgoing external edges must be attached to an external node (degree 1)"
                        ));
                    }
                }
                Flow::Sink => {
                    if let Some(&sink_node) = vertex_map.get(&in_vertex) {
                        (sink_node, orientation, flow)
                    } else if let Some(&source_node) = vertex_map.get(&out_vertex) {
                        data.particle = data.particle.reverse(model);
                        orientation = orientation.reverse();
                        (source_node, orientation, flow)
                    } else {
                        return Err(eyre!(
                            "Incoming external edges must be attached to an external node (degree 1)"
                        ));
                    }
                }
            };

            // debug!(node =  %node_idx, edge_data = ?data ,"adding_external_edge");
            builder.add_external_edge(node_idx, data, final_orientation, final_flow);
            Ok(())
        }

        // debug!("Input:{}", graph.to_dot());
        let mut builder = HedgeGraphBuilder::new();

        let mut tags_to_edge_id = BTreeMap::new();
        let mut vertex_map = AHashMap::new();
        for (i, n) in graph.nodes().iter().enumerate() {
            if n.edges.len() == 1 {
                tags_to_edge_id.insert(n.data.external_tag, n.edges[0]);
            } else {
                vertex_map.insert(i, builder.add_node(ParseVertex::from(&n.data)));
            }
        }

        let mut seen_edges = AHashSet::new();
        let mut generation_type: Option<GenerationType> = None;

        // first add incoming edges
        for (i, &(in_tag, out_tag)) in external_connections.iter().enumerate() {
            if let Some(in_tag) = in_tag {
                let in_edge_idx = tags_to_edge_id[&(in_tag as i32)];
                mark_edge_as_seen(&mut seen_edges, in_edge_idx)?;
                let in_edge = &graph.edges()[in_edge_idx];

                let is_cut_hedge = if let Some(out_id) = out_tag {
                    let out_edge_idx = tags_to_edge_id[&(out_id as i32)];
                    mark_edge_as_seen(&mut seen_edges, out_edge_idx)?;

                    let out_edge = &graph.edges()[out_edge_idx];

                    validate_edge_compatibility(out_edge, in_edge, in_tag, out_id)?;
                    if let Some(existing_type) = &generation_type {
                        if *existing_type != GenerationType::CrossSection {
                            return Err(eyre!(
                                "Cannot have both incoming and outgoing external connections for amplitudes"
                            ));
                        }
                    } else {
                        generation_type = Some(GenerationType::CrossSection);
                    }
                    Some(Hedge(i))
                } else {
                    None
                };

                process_single_connection_internal(
                    in_edge_idx,
                    Flow::Sink,
                    is_cut_hedge,
                    &vertex_map,
                    &mut seen_edges,
                    graph,
                    model,
                    &mut builder,
                )?;
            } else if out_tag.is_some() {
                if let Some(existing_type) = &generation_type {
                    if *existing_type != GenerationType::Amplitude {
                        return Err(eyre!(
                            "Cannot mix single directional connections with bidirectional ones for cross sections"
                        ));
                    }
                } else {
                    generation_type = Some(GenerationType::Amplitude);
                }
            }
        }

        // then add outgoing amplitude edges
        for (in_id, out_id) in external_connections.iter() {
            if let Some(out_id) = out_id
                && in_id.is_none()
            {
                let out_edge_idx = tags_to_edge_id[&(*out_id as i32)];
                mark_edge_as_seen(&mut seen_edges, out_edge_idx)?;

                process_single_connection_internal(
                    out_edge_idx,
                    Flow::Source,
                    None,
                    &vertex_map,
                    &mut seen_edges,
                    graph,
                    model,
                    &mut builder,
                )?;
            }
        }

        // Add internal edges
        for (i, edge) in graph.edges().iter().enumerate() {
            if seen_edges.contains(&i) {
                continue;
            }
            let (source_v, sink_v) = edge.vertices;

            let source = vertex_map[&source_v];
            let sink = vertex_map[&sink_v];
            let data = ParseEdge::from_symbolica_edge(model, &edge.data, None);
            let orientation = data.particle.orientation();
            builder.add_edge(source, sink, data, orientation);
        }

        // then add outgoing cross_section edges
        for (i, &(in_id, out_id)) in external_connections.iter().enumerate() {
            if let Some(out_id) = out_id
                && in_id.is_some()
            {
                let out_edge_idx = tags_to_edge_id[&(out_id as i32)];

                process_single_connection_internal(
                    out_edge_idx,
                    Flow::Source,
                    Some(Hedge(i)),
                    &vertex_map,
                    &mut seen_edges,
                    graph,
                    model,
                    &mut builder,
                )?;
            }
        }

        let mut parsed = ParseGraph {
            global_data: ParseData {
                name: graph_name.as_ref().into(),
                overall_factor: symmetry_factor,
                ..Default::default()
            },
            graph: builder.into(),
        };

        debug!("Parsing {}", parsed.debug_dot());
        parsed.fix_cp_vertex_rules(model)?;
        debug!("Parsing fixed{}", parsed.debug_dot());
        Ok(parsed)
    }

    fn fix_cp_vertex_rules(&mut self, model: &Model) -> Result<()> {
        let mut new_nodes = AHashMap::new();
        for (node_id, neighs, v) in self.graph.iter_nodes() {
            let (particles, vertex_name) = Self::extract_vertex_particles(v);
            let Some(vertex) = vertex_name.map(|a| model.get_vertex_rule(a)) else {
                continue;
            };
            let hedge_vec: Vec<_> = neighs.collect();
            let particles = particles.into_iter().flatten().sorted().collect_vec();
            let mut oriented_particles = extract_oriented_particles_from_vertex_hedges(
                self,
                hedge_vec.iter().copied(),
                model,
            );
            oriented_particles.sort();

            let couplings = vertex.coupling_orders(model);
            // let particles_n = particles.iter().map(|p| p.name.as_str()).collect_vec();
            // let particles_vn = oriented_particles
            //     .iter()
            //     .map(|p| p.name.as_str())
            //     .collect_vec();

            // debug!(
            //     "Comparing  vertex rules particles {:?} with incoming particles {:?}",
            //     particles_n, particles_vn
            // );

            if particles != oriented_particles {
                debug!("Need to change");
                let cp_particles: Vec<_> = oriented_particles
                    .iter()
                    .map(|a| a.get_anti_particle(model))
                    .sorted()
                    .collect();
                if cp_particles == particles {
                    let res = model
                        .particle_set_to_vertex_rules_map
                        .get(&oriented_particles);

                    if let Some(res) = res {
                        let possible: Vec<_> = res
                            .iter()
                            .filter(|a| a.coupling_orders(model) == couplings)
                            .collect();

                        if possible.len() == 1 {
                            new_nodes.insert(node_id, possible[0].clone());
                        } else {
                            let particles = particles.iter().map(|p| p.name.as_str()).collect_vec();

                            return Err(eyre!(
                                "Multiple compatible  vertex rules for {:?}",
                                particles
                            ));
                        }
                    } else {
                        return Err(eyre!(
                            "Failed to find CP vertex rule for particles: {:?} for node {node_id} in graph {}",
                            particles,
                            self.global_data.name,
                        ));
                    }
                } else {
                    let particles = particles.iter().map(|p| p.name.as_str()).collect_vec();
                    return Err(eyre!(
                        "Failed to find CP vertex rule for particles: {:?} for node {node_id} in graph {}",
                        particles,
                        self.global_data.name,
                    ));
                }
            }
        }

        for (node_id, vr) in new_nodes {
            debug!("New vr for {node_id}:{}", vr.name);
            self.graph[node_id].vertex_rule = Some(vr);
        }
        Ok(())
    }
}

impl ParseGraph {
    pub(crate) fn from_parsed(graph: DotGraph, model: &Model) -> Result<Self> {
        let global_data = graph.global_data.into();
        let graph = graph
            .graph
            .map_data_ref_result(
                |_, _, v| Ok(v),
                ParseEdge::parse(model),
                ParseHedge::parse(),
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

/// Edge and vertex numerators
struct NumeratorData {
    color_edge: EdgeVec<Atom>,
    spin_edge: EdgeVec<Atom>,
    color_vertex: Vec<Option<ParamTensor<OrderedStructure<Euclidean, Aind>>>>,
    spin_vertex: Vec<Option<ParamTensor<OrderedStructure<Euclidean, Aind>>>>,
}

impl Graph {
    pub fn dot_serialize(&self) -> String {
        let mut out = String::new();
        self.dot_serialize_fmt(&mut out).unwrap();
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
    ) -> Result<(), std::fmt::Error> {
        let g = DotGraph::from(self);
        g.write_fmt(writer)
    }

    pub(crate) fn from_symbolica_graph(
        model: &Model,
        graph_name: impl AsRef<str>,
        graph: &SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
        symmetry_factor: Atom,
        external_connections: &[(Option<usize>, Option<usize>)],
    ) -> Result<Self> {
        let mut parse_graph = ParseGraph::from_symbolica_graph(
            model,
            graph_name,
            graph,
            symmetry_factor,
            external_connections,
        )?;

        parse_graph
            .graph
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
        Graph::from_parsed(parse_graph, model)
    }

    #[instrument(skip_all, fields(graph= %graph.debug_dot(),name = %graph.global_data.name.as_str()))]
    pub(crate) fn from_parsed(graph: ParseGraph, model: &Model) -> Result<Self> {
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

        let numerators = Self::generate_numerators(&graph, model);

        let initial_state_cut =
            OrientedCut::from_underlying_strict(cut_result.initial_hedges, &graph)?;

        let (global_prefactor, param_builder) = Self::setup_global_prefactor_and_params(
            initial_data.global_prefactor,
            initial_data.add_polarizations,
            initial_data.additional_params.clone(),
            &initial_state_cut,
            &graph,
            model,
        )?;

        let underlying = Self::build_underlying_graph(
            graph,
            &initial_state_cut,
            &numerators,
            model,
            &param_builder,
        )?;

        let loop_momentum_basis = Self::setup_loop_momentum_basis(
            &underlying,
            &cut_result.full_cut,
            &cut_result.lmb_ids,
            &cut_result.xs_ext_id,
        )?;

        let mut g = Graph {
            overall_factor: initial_data.overall_factor,
            polarizations: global_prefactor.polarizations(),
            global_prefactor,
            name: initial_data.name,
            loop_momentum_basis,
            initial_state_cut,
            underlying,
            group_id: initial_data.group_id,
            is_group_master: initial_data.is_group_master,
            param_builder,
        };

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

    fn extract_initial_data(
        parse_graph: &ParseGraph,
        model: &Model,
    ) -> Result<(InitialGraphData, NumGraph)> {
        let hedge_order = parse_graph.hedge_order(model)?;
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
            |h, hd| NumHedgeData {
                num_indices: NumIndices::parse(parse_graph)(h, hd),
                node_order: hedge_order[h.0],
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

    fn generate_numerators(graph: &NumGraph, model: &Model) -> NumeratorData {
        let mut color_edge: EdgeVec<_> = vec![Atom::num(1); graph.n_edges()].into();
        let mut spin_edge: EdgeVec<_> = vec![Atom::i(); graph.n_edges()].into();

        for (p, eid, e) in graph.iter_edges() {
            if let HedgePair::Paired { source, sink } = p {
                let prop = e
                    .data
                    .particle
                    .particle()
                    .map(|p| model.get_propagator_for_particle(&p.name).numerator.clone())
                    .unwrap_or(Atom::num(1));

                color_edge[eid] = graph[source].color_kronekers(&graph[sink]);

                let spin_slots = [
                    &graph[source].num_indices.spin_indices.edge_indices,
                    &graph[sink].num_indices.spin_indices.edge_indices,
                ];

                let momenta = [(Flow::Source, graph[&source]), (Flow::Sink, graph[&sink])];
                let spin_nume = UFO.reindex_spin(&spin_slots, &momenta, prop, |i| {
                    Aind::Edge(usize::from(eid) as u16, i as u16)
                });

                spin_edge[eid] = spin_nume;
            }
        }

        let (color_vertex, spin_vertex) = Self::generate_vertex_numerators(graph);

        NumeratorData {
            color_edge,
            spin_edge,
            color_vertex,
            spin_vertex,
        }
    }

    fn setup_global_prefactor_and_params<'a, A: Into<AtomOrView<'a>>, P: IntoIterator<Item = A>>(
        mut global_prefactor: GlobalPrefactor,
        add_polarizations: bool,
        params: P,
        initial_state_cut: &OrientedCut,
        graph: &NumGraph,
        model: &Model,
    ) -> Result<(GlobalPrefactor, ParamBuilder)> {
        debug!("Initial state cut: {}", graph.dot(&initial_state_cut.left));

        if add_polarizations {
            let external_edges = initial_state_cut
                .left
                .union(&initial_state_cut.right)
                .union(&graph.external_filter());
            let polarizations = graph.generate_polarizations_of(&external_edges);
            global_prefactor.projector *= polarizations;
        }

        let polarizations = global_prefactor.polarizations();
        let param_builder =
            ParamBuilder::new(&(&polarizations, graph), model, &graph.lmb(), params);

        Ok((global_prefactor, param_builder))
    }

    #[allow(clippy::type_complexity)]
    fn generate_vertex_numerators(
        graph: &NumGraph,
    ) -> (
        Vec<Option<ParamTensor<OrderedStructure<Euclidean, Aind>>>>,
        Vec<Option<ParamTensor<OrderedStructure<Euclidean, Aind>>>>,
    ) {
        let mut color_vertex: Vec<Option<ParamTensor<OrderedStructure<Euclidean, Aind>>>> =
            vec![None; graph.n_nodes()];
        let mut spin_vertex = color_vertex.clone();

        for (ni, c, v) in graph.iter_nodes() {
            let mut color_slots = vec![];
            let mut spin_slots = vec![];
            let mut order = vec![];
            let mut momenta = vec![];
            for h in c {
                color_slots.push(&graph[h].num_indices.color_indices.vertex_indices);
                spin_slots.push(&graph[h].num_indices.spin_indices.vertex_indices);
                order.push(graph[h].node_order);
                momenta.push((graph.flow(h), graph[&h]));
            }

            let perm = Permutation::sort(&order);
            perm.apply_slice_in_place(&mut color_slots);
            perm.apply_slice_in_place(&mut spin_slots);
            perm.apply_slice_in_place(&mut momenta);

            let Some(vertex_rule) = &v.vertex_rule else {
                continue;
            };

            let [mut color_structure, couplings, mut spin_structure] =
                vertex_rule.tensors(ni.aind(1), ni.aind(0));

            spin_structure.map_data_mut(|a| {
                *a = UFO.reindex_spin(&spin_slots, &momenta, (*a).clone(), |u| ni.aind(u as u16));
            });

            // couplings.map_data_mut(|a| *a = UFO.normalize_complex((*a).clone()));

            color_structure.map_data_mut(|a| {
                *a = UFO.reindex_color(&color_slots, (*a).clone(), |u| ni.aind(u as u16));
            });

            spin_vertex[ni.0] = Some(spin_structure.contract(&couplings).unwrap());
            color_vertex[ni.0] = Some(color_structure);
        }

        (color_vertex, spin_vertex)
    }

    fn build_underlying_graph(
        graph: NumGraph,
        initial_state_cut: &OrientedCut,
        numerators: &NumeratorData,
        model: &Model,
        param_builder: &ParamBuilder,
    ) -> Result<UnderlyingGraph> {
        let mut vertex_color_nums = numerators.color_vertex.clone();
        let mut vertex_spin_nums = numerators.spin_vertex.clone();
        graph.map_result(
            |_, i, v| {
                let num = if let Some(num) = v.num {
                    num
                } else {
                    vertex_spin_nums[i.0]
                        .take()
                        .unwrap()
                        .contract(&vertex_color_nums[i.0].take().unwrap())
                        .unwrap()
                        .scalar()
                        .unwrap()
                };

                let dod = if let Some(d) = v.dod { d } else { DOD::dod(&num) };

                Ok(Vertex {
                    name: v.name.unwrap_or(i.to_string()),
                    num,
                    dod,
                    vertex_rule: v.vertex_rule,
                })
            },
            |_, _, p, eid, ed| {
                let e = ed.data;
                if e.particle.is_fermion() && !e.particle.is_self_antiparticle()&&  e.particle.orientation() != ed.orientation {
                    return Err(eyre!(
                        "Edge orientation {:?} does not match particle orientation {:?} for edge {},{}",
                        ed.orientation,
                        e.particle.orientation(),
                        eid,
                        e
                    ));
                    }


                let num = e.num.unwrap_or(if initial_state_cut.left.intersects(&p) {
                    numerators.color_edge[eid].clone()
                } else {
                    &numerators.color_edge[eid] * &numerators.spin_edge[eid]
                });

                let dod = if let Some(d) = e.dod {
                    d
                } else {
                    DOD::dod(&num) - 2
                };

                Ok(EdgeData::new(
                    Edge {
                        mass: EdgeMass::from_atom(e.particle.mass_atom(), model, param_builder)?,
                        is_dummy: e.is_dummy,
                        name: e.name.unwrap_or(eid.to_string()),
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
        )
    }

    fn setup_loop_momentum_basis(
        underlying: &UnderlyingGraph,
        full_cut: &SuBitGraph,
        lmb_ids: &BTreeMap<LoopIndex, EdgeIndex>,
        xs_ext_id: &BTreeMap<Hedge, (EdgeIndex, Hedge)>,
    ) -> Result<LoopMomentumBasis> {
        debug!("{}", underlying.dot(full_cut));

        let mut loop_momentum_basis = if let Some(_) = full_cut.included_iter().next() {
            // let tree = SimpleTraversalTree::depth_first_traverse(
            //     underlying,
            //     full_cut,
            //     &underlying.node_id(i),
            //     None,
            // )
            // .unwrap();

            let mut full = underlying.full_filter();

            for (p, _, i) in underlying.iter_edges() {
                if i.data.is_dummy {
                    full.sub(p);
                }
            }

            // let covers = tree.covers(&full);
            // assert_eq!(
            //     full,
            //     covers,
            //     "Tree does not cover all: {}, lmb specification must be wrong",
            //     underlying.dot(&covers)
            // );
            let external = underlying.internal_crown(&full);
            underlying.lmb_impl(&full, &full_cut, external)
        } else {
            return Err(eyre!(
                "No included edges found in full_cut for loop momentum basis setup"
            ));
        };

        let inv_lmb_ids: BTreeMap<_, _> = lmb_ids
            .iter()
            .map(|(k, v)| {
                // debug!("v{v}k{k}");
                (*v, *k)
            })
            .collect();

        for e in 0..xs_ext_id.len() {
            let (l, _) = loop_momentum_basis
                .loop_edges
                .iter()
                .find_position(|a| *a == &EdgeIndex(e))
                .unwrap();

            loop_momentum_basis.put_loop_to_ext(LoopIndex(l));
        }

        // Process swaps until no more changes needed
        let mut swapped = true;
        while swapped {
            swapped = false;
            for i in 0..loop_momentum_basis.loop_edges.len() {
                if let Some(&target_pos) =
                    inv_lmb_ids.get(&loop_momentum_basis.loop_edges[LoopIndex(i)])
                    && target_pos.0 < loop_momentum_basis.loop_edges.len()
                    && target_pos.0 != i
                {
                    loop_momentum_basis.swap_loops(LoopIndex(i), target_pos);
                    swapped = true;
                    break;
                }
            }
        }

        Ok(loop_momentum_basis)
    }

    pub fn from_dot(graph: DotGraph, model: &Model) -> Result<Self> {
        Self::from_parsed(ParseGraph::from_parsed(graph, model)?, model)
    }
    pub fn from_file<P>(p: P, model: &Model) -> Result<Vec<Self>>
    where
        P: AsRef<Path>,
    {
        Self::from_path(p, model)
    }

    pub fn from_path<P>(p: P, model: &Model) -> Result<Vec<Self>>
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
                let graphs = Self::from_single_file(&dot_file, model)?;
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
            Self::from_single_file(path, model)
        }
    }

    fn from_single_file<P>(p: P, model: &Model) -> Result<Vec<Self>>
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
            HedgeParseError::GraphFromFile(e) => match e {
                dot_parser::ast::GraphFromFileError::FileError(e) => color_eyre::Report::from(e)
                    .with_note(|| format!("Tried to access the file at:{}", p.as_ref().display())),
                dot_parser::ast::GraphFromFileError::ParseError(e) => {
                    eyre!("Dot parsing error: {}", e)
                }
                dot_parser::ast::GraphFromFileError::PestParseError(e) => {
                    color_eyre::Report::from(e)
                }
            },
            HedgeParseError::ParseError(i) => color_eyre::Report::from(i),
            _ => {
                eyre!("Hedge parse error")
            }
        })?;
        Self::from_hedge_graph_set(hedge_graph_set, model)
    }

    pub fn from_string<Str: AsRef<str>>(s: Str, model: &Model) -> Result<Vec<Self>> {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_string(s).map_err(|a| match a {
            HedgeParseError::GraphFromFile(e) => match e {
                dot_parser::ast::GraphFromFileError::FileError(e) => color_eyre::Report::from(e),
                dot_parser::ast::GraphFromFileError::ParseError(e) => {
                    eyre!("Dot parsing error: {}", e)
                }
                dot_parser::ast::GraphFromFileError::PestParseError(e) => {
                    color_eyre::Report::from(e)
                }
            },
            HedgeParseError::ParseError(i) => color_eyre::Report::from(i),
            _ => {
                eyre!("Hedge parse error")
            }
        })?;

        Self::from_hedge_graph_set(hedge_graph_set, model)
    }

    fn from_hedge_graph_set(
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

        for (graph, global_data) in set.set.into_iter().zip(set.global_data.into_iter()) {
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
