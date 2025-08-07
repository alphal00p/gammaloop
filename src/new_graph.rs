use std::borrow::Borrow;

use ahash::HashSet;
use bitvec::vec::BitVec;
use color_eyre::Result;
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeData, EdgeIndex, EdgeVec, Flow, HedgePair},
    subgraph::{InternalSubGraph, ModifySubgraph, SubGraph, SubGraphOps},
    HedgeGraph, NodeIndex,
};
use log::debug;
use momtrop::float::MomTropFloat;

use rand::{rngs::SmallRng, Rng, SeedableRng};
use spenso::algebra::{algebraic_traits::IsZero, complex::Complex};
// use petgraph::Direction::Outgoing;
use idenso::metric::MS;
use smartstring::{LazyCompact, SmartString};
use spenso::structure::representation::Minkowski;
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    id::Replacement,
    parse,
};
use typed_index_collections::TiVec;

use crate::{
    cff::generation::ShiftRewrite,
    graph::BareGraph,
    model::{ArcParticle, EdgeSlots, VertexSlots},
    momentum::{Dep, ExternalMomenta, Helicity, SignOrZero},
    momentum_sample::{ExternalFourMomenta, ExternalIndex, LoopMomenta},
    new_gammaloop_integrand::LmbMultiChannelingSetup,
    numerator::GlobalPrefactor,
    signature::{ExternalSignature, SignatureLike},
    utils::{external_energy_atom_from_index, ose_atom_from_index, FloatLike, F, GS},
    Externals,
};

pub mod global;

#[derive(Clone, Copy, bincode_trait_derive::Encode, bincode_trait_derive::Decode, Default)]
pub struct VertexOrder(pub u8);

#[derive(Clone, bincode_trait_derive::Encode, bincode_trait_derive::Decode)]
#[trait_decode(trait = crate::GammaLoopContext)]
pub struct Graph {
    pub overall_factor: Atom,
    pub name: String,
    pub underlying: HedgeGraph<Edge, Vertex, NumHedgeData>,
    pub loop_momentum_basis: LoopMomentumBasis,
    pub global_prefactor: GlobalPrefactor,
}

// impl Deref for Graph {
//     type Target = HedgeGraph<Edge, Vertex>;
// }

impl From<BareGraph> for Graph {
    fn from(value: BareGraph) -> Self {
        todo!()
    }
}

pub mod feynman_graph;
pub use feynman_graph::FeynmanGraph;

impl Graph {
    pub(crate) fn random_externals(&self, seed: u64) -> Externals {
        let mut rng = SmallRng::seed_from_u64(seed);
        let mom_range = -10.0..10.0;

        let mut momenta = vec![ExternalMomenta::Dependent(Dep::Dep)];
        let mut helicities = vec![];
        let ext = self.external_filter();

        for (p, i, d) in self.iter_edges_of(&ext) {
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

    pub(crate) fn new(
        name: SmartString<LazyCompact>,
        multiplicity: Atom,
        underlying: HedgeGraph<Edge, Vertex, NumHedgeData>,
    ) -> Result<Self> {
        Ok(Self {
            name: name.to_string(),
            overall_factor: multiplicity,
            loop_momentum_basis: underlying.lmb(&underlying.full_filter()),
            underlying,
            global_prefactor: GlobalPrefactor::default(),
        })
    }

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
