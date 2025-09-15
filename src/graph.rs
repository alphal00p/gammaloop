use std::ops::Index;

use ahash::{AHashMap, AHashSet, HashSet};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use global::ParseData;
use itertools::Itertools;
use linnet::{
    half_edge::{
        builder::HedgeGraphBuilder,
        involution::{EdgeData, EdgeIndex, Flow, Hedge, HedgePair, Orientation},
        subgraph::{OrientedCut, SubGraph},
        HedgeGraph,
    },
    parser::DotGraph,
};
use log::debug;

use parse::ParseGraph;
use rand::{rngs::SmallRng, Rng, SeedableRng};
// use petgraph::Direction::Outgoing;
use symbolica::{atom::Atom, graph::Graph as SymbolicaGraph};
use typed_index_collections::TiVec;

use crate::{
    define_index,
    feyngen::{
        diagram_generator::{EdgeColor, NodeColorWithVertexRule},
        GenerationType,
    },
    gammaloop_integrand::{LmbMultiChannelingSetup, ParamBuilder},
    graph::{edge::ParseEdge, vertex::ParseVertex},
    model::Model,
    momentum::{Dep, ExternalMomenta, PolDef},
    momentum_sample::ExternalIndex,
    numerator::{symbolica_ext::AtomCoreExt, GlobalPrefactor, ParsingNet},
    settings::runtime::kinematic::Externals,
    utils::{ose_atom_from_index, Length, F, GS},
};

pub mod global;

#[derive(Clone, Copy, bincode_trait_derive::Encode, bincode_trait_derive::Decode, Default)]
pub struct VertexOrder(pub u8);

define_index! {pub struct GroupId;}

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub struct Graph {
    pub overall_factor: Atom,
    pub name: String,
    pub group_id: Option<GroupId>,
    pub is_group_master: bool,
    pub underlying: HedgeGraph<Edge, Vertex, NumHedgeData>,
    pub loop_momentum_basis: LoopMomentumBasis,
    pub param_builder: ParamBuilder,
    pub global_prefactor: GlobalPrefactor,
    /// The cross section initial state cut
    /// Only relevant for cross sections, but stored here for the parsing
    pub initial_state_cut: OrientedCut,
    pub polarizations: Vec<(PolDef, Atom)>,
}

// impl Deref for Graph {
//     type Target = HedgeGraph<Edge, Vertex>;
// }

pub mod feynman_graph;
pub use feynman_graph::FeynmanGraph;
pub mod ext;
impl Graph {
    pub fn debug_dot(&self) -> String {
        DotGraph::from(self).debug_dot()
    }

    /// With wrapped color, so that it doesn't enter the network as a tensor. Can unwrap using `unwrap_function`
    /// Contains the parametric sign on the OSE
    pub(crate) fn global_network(&self) -> ParsingNet {
        let net =
            (&self.global_prefactor.num * &self.global_prefactor.projector * &self.overall_factor)
                .wrap_color(GS.color_wrap)
                .parse_into_net()
                .unwrap();

        let mut reps = Vec::new();
        for (p, eid, _) in self.iter_edges() {
            if p.is_paired() {
                reps.push(GS.add_parametric_sign(eid));
            }
        }
        net.replace_multiple(&reps)
    }

    pub(crate) fn random_externals(&self, seed: u64) -> Externals {
        let mut rng = SmallRng::seed_from_u64(seed);
        let mom_range = -10.0..10.0;

        let mut momenta = vec![ExternalMomenta::Dependent(Dep::Dep)];
        let mut helicities = vec![];
        let ext = self.external_filter();

        for (_, _, d) in self.iter_edges_of(&ext) {
            let hel = d.data.random_helicity(seed);
            helicities.push(hel);
            if helicities.len() == 2 {
                continue;
            }
            let mom = ExternalMomenta::Independent([
                F(rng.random_range(mom_range.clone())),
                F(rng.random_range(mom_range.clone())),
                F(rng.random_range(mom_range.clone())),
                F(rng.random_range(mom_range.clone())),
            ]);

            momenta.push(mom);
        }
        Externals::Constant {
            momenta,
            helicities,
        }
    }

    // pub(crate) fn new(
    //     name: SmartString<LazyCompact>,
    //     multiplicity: Atom,
    //     underlying: HedgeGraph<Edge, Vertex, NumHedgeData>,
    // ) -> Result<Self> {
    //     Ok(Self {
    //         name: name.to_string(),
    //         overall_factor: multiplicity,
    //         loop_momentum_basis: underlying.lmb(&underlying.full_filter()),
    //         group_id: None,
    //         is_group_master: false,
    //         underlying,
    //         global_prefactor: GlobalPrefactor::default(),
    //         polarizations: vec![],
    //     })
    // }

    pub(crate) fn build_multi_channeling_channels(
        &self,
        lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    ) -> LmbMultiChannelingSetup {
        let mut channels: Vec<LmbIndex> = Vec::new();

        // Filter out channels that are non-singular, or have the same singularities as another channel already included
        for (lmb_index, lmb) in lmbs.iter_enumerated() {
            let massless_edges = lmb
                .loop_edges
                .iter()
                .filter(|&edge_id| self.underlying[*edge_id].particle.is_massless())
                .collect_vec();

            if massless_edges.is_empty() {
                continue;
            }

            if channels.iter().any(|channel| {
                let basis_of_included_channel = &lmbs[*channel].loop_edges;
                massless_edges
                    .iter()
                    .all(|edge_id| basis_of_included_channel.contains(edge_id))
            }) {
                continue;
            }

            // only for 1 to n for now, assuming center of mass
            if channels.iter().any(|channel| {
                let massless_edges_of_included_channel = lmbs[*channel]
                    .loop_edges
                    .iter()
                    .filter(|&edge_id| self.underlying[*edge_id].particle.is_massless())
                    .collect_vec();

                let loop_signatures_of_massless_edges_of_included_channel =
                    massless_edges_of_included_channel
                        .iter()
                        .map(|edge_index| {
                            self.loop_momentum_basis.edge_signatures[**edge_index]
                                .internal
                                .first_abs()
                        })
                        .collect::<HashSet<_>>();

                let loop_signatures_of_massless_edges_of_potential_channel = massless_edges
                    .iter()
                    .map(|edge_index| {
                        self.loop_momentum_basis.edge_signatures[**edge_index]
                            .internal
                            .first_abs()
                    })
                    .collect::<HashSet<_>>();

                loop_signatures_of_massless_edges_of_included_channel
                    == loop_signatures_of_massless_edges_of_potential_channel
            }) {
                continue;
            }

            channels.push(lmb_index);
        }

        let channels: TiVec<_, _> = if !channels.is_empty() {
            channels.into_iter().sorted().collect()
        } else {
            let current_lmb_index = lmbs
                .iter_enumerated()
                .find(|(_lmb_index, lmb)| lmb.loop_edges == self.loop_momentum_basis.loop_edges)
                .map(|(lmb_index, _)| lmb_index)
                .unwrap();

            vec![current_lmb_index].into()
        };

        debug!(
            "number of lmbs: {}, number of channels: {}",
            lmbs.len(),
            channels.len()
        );

        LmbMultiChannelingSetup { channels }
    }

    pub(crate) fn iter_loop_edges(
        &self,
    ) -> impl Iterator<Item = (HedgePair, EdgeIndex, EdgeData<&Edge>)> {
        self.underlying.iter_edges().filter(|(_, edge_index, _)| {
            self.loop_momentum_basis.edge_signatures[*edge_index]
                .internal
                .iter()
                .any(|sign| sign.is_sign())
        })
    }

    pub(crate) fn iter_non_loop_edges(
        &self,
    ) -> impl Iterator<Item = (HedgePair, EdgeIndex, EdgeData<&Edge>)> {
        self.underlying.iter_edges().filter(|(_, edge_index, _)| {
            self.loop_momentum_basis.edge_signatures[*edge_index]
                .internal
                .iter()
                .all(|sign| sign.is_zero())
        })
    }

    pub(crate) fn from_symbolica_graph(
        model: &Model,
        graph_name: impl AsRef<str>,
        graph: &SymbolicaGraph<NodeColorWithVertexRule, EdgeColor>,
        symmetry_factor: Atom,
        external_connections: &[(Option<usize>, Option<usize>)],
    ) -> Result<Self> {
        // println!("Input:{}", graph.to_dot());

        //let builder = HedgeGraphBuilder::new();

        // let input_graph: HedgeGraph<_, _, ()> = graph.clone().into();
        //

        let mut builder = HedgeGraphBuilder::new();
        let mut map = AHashMap::new();

        let mut external_tags = Vec::new();

        for (i, node) in graph.nodes().iter().enumerate() {
            if node.edges.len() == 1 {
                external_tags.push((i, node.edges[0]));
                continue;
            }
            map.insert(i, builder.add_node(ParseVertex::from(&node.data)));
        }

        external_tags.sort_by_key(|&(i, _)| i);

        let mut seen = AHashSet::new();
        let mut gen_type = None;

        for (i, (in_id, out_id)) in external_connections.iter().enumerate() {
            // println!(
            //     "External connection {}: in_id={:?}, out_id={:?}",
            //     i + 1,
            //     in_id,
            //     out_id
            // );
            match (in_id, out_id) {
                (Some(in_id), Some(out_id)) => {
                    let out_ind = external_tags[*out_id - 1].1;
                    let in_ind = external_tags[*in_id - 1].1;
                    if !seen.insert(in_ind) || !seen.insert(out_ind) {
                        return Err(eyre::eyre!("External connections must be unique {in_ind} or {out_ind} have been seen before"));
                    }
                    if let Some(GenerationType::Amplitude) = gen_type {
                        return Err(eyre::eyre!("Cannot have both incoming and outgoing external connections for amplitudes"));
                    }
                    gen_type = Some(GenerationType::CrossSection);

                    let out_edge = &graph.edges()[out_ind];
                    let (out_out_v, out_in_v) = out_edge.vertices;
                    let orientation = Orientation::from(out_edge.directed);

                    if out_edge.data.pdg < 0 {
                        orientation.reverse();
                    }

                    let edata =
                        ParseEdge::from_symbolica_edge(model, &out_edge.data, Some(Hedge(i)));

                    let in_edge = &graph.edges()[in_ind];
                    let (in_out_v, in_in_v) = in_edge.vertices;
                    assert_eq!(in_edge.directed, out_edge.directed, "External edges must have the same directedness, for edge ids {in_id} and {out_id} found {in_edge:?} and {out_edge:?}");
                    assert_eq!(in_edge.data.pdg.abs(), out_edge.data.pdg.abs(), "External edges must have the same pdg in abs, for edge ids {in_id} and {out_id} found {in_edge:?} and {out_edge:?}");

                    if let Some(sink) = map.get(&out_in_v) {
                        builder.add_external_edge(
                            *sink,
                            edata.clone(),
                            orientation.reverse(),
                            Flow::Source,
                        );
                    } else if let Some(source) = map.get(&out_out_v) {
                        builder.add_external_edge(
                            *source,
                            edata.clone(),
                            orientation,
                            Flow::Source,
                        );
                    } else {
                        return Err(eyre::eyre!("Outgoing external edges must be attached to an external node (degree 1)"));
                    }

                    if let Some(sink) = map.get(&in_in_v) {
                        builder.add_external_edge(*sink, edata.clone(), orientation, Flow::Sink);
                    } else if let Some(source) = map.get(&in_out_v) {
                        builder.add_external_edge(
                            *source,
                            edata.clone(),
                            orientation.reverse(),
                            Flow::Sink,
                        );
                    } else {
                        return Err(eyre::eyre!("Incoming external edges must be attached to an external node (degree 1)"));
                    }
                }
                (None, None) => {
                    return Err(eyre::eyre!("External connections must have at least one of incoming or outgoing defined"));
                }
                //Incoming external edge, needs to have incoming momentum
                (Some(in_id), None) => {
                    let in_ind = external_tags[*in_id - 1].1;
                    if !seen.insert(in_ind) {
                        return Err(eyre::eyre!(
                            "External connections must be unique {in_ind} has been seen before"
                        ));
                    }
                    if let Some(GenerationType::CrossSection) = gen_type {
                        return Err(eyre::eyre!("Cannot single incoming for crossection"));
                    }
                    gen_type = Some(GenerationType::Amplitude);

                    let edge = &graph.edges()[in_ind];
                    let edata = ParseEdge::from_symbolica_edge(model, &edge.data, None);
                    let (out_v, in_v) = edge.vertices;
                    let orientation = Orientation::from(edge.directed);

                    if let Some(sink) = map.get(&in_v) {
                        builder.add_external_edge(*sink, edata, orientation, Flow::Sink);
                    } else if let Some(source) = map.get(&out_v) {
                        builder.add_external_edge(
                            *source,
                            edata,
                            orientation.reverse(),
                            Flow::Sink,
                        );
                    } else {
                        return Err(eyre::eyre!("Incoming external edges must be attached to an external node (degree 1)"));
                    }
                }
                //Outgoing  external edge, needs to have outgoing momentum
                (None, Some(out_id)) => {
                    let out_ind = external_tags[*out_id - 1].1;
                    if !seen.insert(out_ind) {
                        return Err(eyre::eyre!(
                            "External connections must be unique {out_ind} has been seen before"
                        ));
                    }
                    if let Some(GenerationType::CrossSection) = gen_type {
                        return Err(eyre::eyre!("Cannot single outgoing for crossection"));
                    }
                    gen_type = Some(GenerationType::Amplitude);

                    let edge = &graph.edges()[out_ind];
                    let edata = ParseEdge::from_symbolica_edge(model, &edge.data, None);
                    let (out_v, in_v) = edge.vertices;
                    let orientation = Orientation::from(edge.directed);

                    if let Some(sink) = map.get(&in_v) {
                        builder.add_external_edge(
                            *sink,
                            edata,
                            orientation.reverse(),
                            Flow::Source,
                        );
                    } else if let Some(source) = map.get(&out_v) {
                        builder.add_external_edge(*source, edata, orientation, Flow::Source);
                    } else {
                        return Err(eyre::eyre!("Incoming external edges must be attached to an external node (degree 1)"));
                    }
                }
            }
        }

        for (i, edge) in graph.edges().iter().enumerate() {
            if seen.contains(&i) {
                continue;
            }
            let vertices = edge.vertices;
            let source = map[&vertices.0];
            let sink = map[&vertices.1];
            builder.add_edge(
                source,
                sink,
                ParseEdge::from_symbolica_edge(model, &edge.data, None),
                edge.directed,
            );
        }

        Graph::from_parsed(
            ParseGraph {
                global_data: ParseData {
                    name: graph_name.as_ref().into(),
                    overall_factor: symmetry_factor,
                    ..Default::default()
                },
                graph: builder.into(),
            },
            model,
        )
    }
}

impl ParseEdge {
    pub fn from_symbolica_edge(
        model: &Model,
        edge_color: &EdgeColor,
        is_cut: Option<Hedge>,
    ) -> Self {
        let particle = model.get_particle_from_pdg(edge_color.pdg);
        let mut e = ParseEdge::new(particle);
        e.is_cut = is_cut;
        e
    }
}

impl From<&NodeColorWithVertexRule> for ParseVertex {
    fn from(value: &NodeColorWithVertexRule) -> Self {
        value.vertex_rule.clone().into()
    }
}

impl Graph {
    // pub(crate) fn apply_vertex_rule(&self, node_id: NodeIndex) -> Option<[DataTensor<Atom>; 3]> {
    //     // self.underlying[node_id].vertex_info.apply_vertex_rule(
    //     //     &self.underlying.add_signs_to_edges(node_id),
    //     //     Into::<usize>::into(node_id),
    //     //     &self.vertex_slots[node_id],
    //     // )
    // }
}

pub mod edge;
pub mod parse;
pub use edge::Edge;
pub mod hedge_data;
pub use hedge_data::NumHedgeData;
pub mod vertex;
pub use vertex::Vertex;

pub mod lmb;
pub use lmb::{LMBext, LmbIndex, LoopMomentumBasis};

#[derive(Debug, Clone, Encode, Decode, PartialEq, Eq)]
pub struct GraphGroup {
    master: usize,
    remaining: Vec<usize>,
}

impl GraphGroup {
    pub(crate) fn master(&self) -> usize {
        self.master
    }

    pub(crate) fn iter_enumerated(&self) -> impl Iterator<Item = (GraphGroupPosition, usize)> + '_ {
        self.into_iter()
            .enumerate()
            .map(|(i, graph_id)| (GraphGroupPosition(i), graph_id))
    }

    pub(crate) fn find_position(&self, graph_id: usize) -> Option<GraphGroupPosition> {
        self.into_iter()
            .position(|id| id == graph_id)
            .map(GraphGroupPosition)
    }
}

impl<'a> IntoIterator for &'a GraphGroup {
    type Item = usize;
    type IntoIter = std::iter::Chain<
        std::array::IntoIter<usize, 1>,
        std::iter::Copied<std::slice::Iter<'a, usize>>,
    >;

    fn into_iter(self) -> Self::IntoIter {
        [self.master]
            .into_iter()
            .chain(self.remaining.iter().copied())
    }
}

impl Length for GraphGroup {
    fn len(&self) -> usize {
        1 + self.remaining.len()
    }
}

impl Index<GraphGroupPosition> for GraphGroup {
    type Output = usize;

    fn index(&self, index: GraphGroupPosition) -> &Self::Output {
        if index.0 == 0 {
            &self.master
        } else {
            &self.remaining[index.0 - 1]
        }
    }
}

define_index! {pub struct GraphGroupPosition;}

// impl From<BareGraph> for HedgeGraph<Edge, Vertex, VertexOrder> {
//     fn from(value: BareGraph) -> Self {
//         let mut node_map = AHashMap::default();
//         let mut external_nodes_to_be_dropped = AHashSet::default();
//         let mut builder = HedgeGraphBuilder::new();

//         for (bare_node_id, bare_node) in value.vertices.iter().enumerate() {
//             let is_external = bare_node.edges.len() == 1;
//             if is_external {
//                 external_nodes_to_be_dropped.insert(bare_node_id);
//                 continue;
//             }

//             let node_id = builder.add_node(bare_node.clone().into());
//             node_map.insert(bare_node_id, node_id);
//         }

//         for (bare_edge_id, bare_edge) in value.edges.into_iter().enumerate() {
//             match bare_edge.edge_type {
//                 EdgeType::Incoming => {
//                     assert!(external_nodes_to_be_dropped.contains(&bare_edge.vertices[0]));
//                     builder.add_external_edge(
//                         node_map[&bare_edge.vertices[1]].add_data(VertexOrder(0)),
//                         bare_edge.into(),
//                         Orientation::Default,
//                         Flow::Sink,
//                     );
//                 }
//                 EdgeType::Outgoing => {
//                     assert!(external_nodes_to_be_dropped.contains(&bare_edge.vertices[1]));
//                     builder.add_external_edge(
//                         node_map[&bare_edge.vertices[0]].add_data(VertexOrder(0)),
//                         bare_edge.into(),
//                         Orientation::Default,
//                         Flow::Source,
//                     );
//                 }
//                 EdgeType::Virtual => {
//                     let source = &bare_edge.vertices[0];
//                     let (source_pos, _) = value.vertices[*source]
//                         .edges
//                         .iter()
//                         .find_position(|i| &bare_edge_id == *i)
//                         .unwrap();
//                     let sink = &bare_edge.vertices[1];
//                     let (sink_pos, _) = value.vertices[*sink]
//                         .edges
//                         .iter()
//                         .find_position(|i| &bare_edge_id == *i)
//                         .unwrap();

//                     builder.add_edge(
//                         node_map[&bare_edge.vertices[0]].add_data(VertexOrder(source_pos as u8)),
//                         node_map[&bare_edge.vertices[1]].add_data(VertexOrder(sink_pos as u8)),
//                         bare_edge.into(),
//                         Orientation::Default,
//                     );
//                 }
//             }
//         }

//         builder.build()
//     }
// }

#[derive(
    Debug, Copy, Clone, PartialEq, Eq, bincode_trait_derive::Encode, bincode_trait_derive::Decode,
)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct ExternalConnection {
    pub incoming_index: ExternalIndex,
    pub outgoing_index: ExternalIndex,
}

pub(crate) fn get_cff_inverse_energy_product_impl<E, V, H, S: SubGraph>(
    graph: &HedgeGraph<E, V, H>,
    subgraph: &S,
    contract_edges: &[EdgeIndex],
) -> Atom {
    Atom::num(1)
        / graph
            .iter_edges_of(subgraph)
            .filter_map(|(pair, edge_index, _)| match pair {
                HedgePair::Paired { .. } => {
                    if contract_edges.contains(&edge_index) {
                        None
                    } else {
                        Some(-Atom::num(2) * ose_atom_from_index(edge_index))
                    }
                }
                _ => None,
            })
            .reduce(|acc, x| acc * x)
            .unwrap_or_else(|| Atom::num(1))
}
