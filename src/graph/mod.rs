use std::ops::Index;

use ahash::AHashSet;
use bincode_trait_derive::{Decode, Encode};
use idenso::color::ColorSimplifier;
use itertools::Itertools;
use linnet::{
    half_edge::{
        HedgeGraph,
        involution::{EdgeData, EdgeIndex, Hedge, HedgePair},
        subgraph::{
            HedgeNode, Inclusion, ModifySubSet, OrientedCut, SuBitGraph, SubSetLike, SubSetOps,
            subset::SubSet,
        },
    },
    parser::DotGraph,
};
use tracing::debug;

use rand::{Rng, SeedableRng, rngs::SmallRng};
// use petgraph::Direction::Outgoing;
use symbolica::atom::Atom;
use tracing::warn;
use typed_index_collections::TiVec;

use crate::{
    define_index,
    feyngen::diagram_generator::evaluate_overall_factor,
    integrands::process::{LmbMultiChannelingSetup, ParamBuilder},
    settings::runtime::kinematic::improvement::PhaseSpaceImprovementSettings,
    momentum::{Dep, ExternalMomenta, PolDef},
    momentum::sample::ExternalIndex,
    numerator::{GlobalPrefactor, ParsingNet, symbolica_ext::AtomCoreExt},
    settings::runtime::kinematic::Externals,
    utils::{F, GS, Length, ose_atom_from_index},
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
    pub tree_edges: SuBitGraph,
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
        let net = (&self.global_prefactor.num
            * &self.global_prefactor.projector
            * evaluate_overall_factor(self.overall_factor.as_view()))
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
        let ext: SuBitGraph = self.external_filter();

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
            improvement_settings: PhaseSpaceImprovementSettings::default(),
            f_64_cache: None,
            f_128_cache: None,
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
            //if channels.iter().any(|channel| {
            //    let massless_edges_of_included_channel = lmbs[*channel]
            //        .loop_edges
            //        .iter()
            //        .filter(|&edge_id| self.underlying[*edge_id].particle.is_massless())
            //        .collect_vec();

            //    let loop_signatures_of_massless_edges_of_included_channel =
            //        massless_edges_of_included_channel
            //            .iter()
            //            .map(|edge_index| {
            //                self.loop_momentum_basis.edge_signatures[**edge_index]
            //                    .internal
            //                    .first_abs()
            //            })
            //            .collect::<HashSet<_>>();

            //    let loop_signatures_of_massless_edges_of_potential_channel = massless_edges
            //        .iter()
            //        .map(|edge_index| {
            //            self.loop_momentum_basis.edge_signatures[**edge_index]
            //                .internal
            //                .first_abs()
            //        })
            //        .collect::<HashSet<_>>();

            //    loop_signatures_of_massless_edges_of_included_channel
            //        == loop_signatures_of_massless_edges_of_potential_channel
            //}) {
            //    continue;
            //}

            channels.push(lmb_index);
        }

        let channels: TiVec<_, _> = if !channels.is_empty() {
            channels.into_iter().sorted().collect()
        } else {
            let current_lmb_index = lmbs
                .iter_enumerated()
                .find(|(_lmb_index, lmb)| lmb.loop_edges == self.loop_momentum_basis.loop_edges)
                .map(|(lmb_index, _)| lmb_index);

            if let Some(current_lmb_index) = current_lmb_index {
                vec![current_lmb_index].into()
            } else {
                warn!("lmb not in the list of all lmb, please check if lmbs are sorted!");
                vec![].into()
            }
        };

        debug!(
            "number of lmbs: {}, number of channels: {}",
            lmbs.len(),
            channels.len()
        );

        LmbMultiChannelingSetup {
            channels,
            graph: self.clone(),
            all_bases: lmbs.clone(),
        }
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

    pub(crate) fn get_source_and_target(&self) -> (HedgeNode, HedgeNode) {
        let mut source_nodes = AHashSet::new();
        let mut target_nodes = AHashSet::new();

        for (hedge_pair, _, _) in self.underlying.iter_edges_of(&self.initial_state_cut) {
            match hedge_pair {
                HedgePair::Split { source, sink, .. } => {
                    let source_node = self.underlying.node_id(sink);
                    let sink_node = self.underlying.node_id(source);

                    source_nodes.insert(source_node);
                    target_nodes.insert(sink_node);
                }
                _ => {
                    unreachable!();
                }
            }
        }

        // They don't need to be the same!
        // assert_eq!(
        //     source_nodes.len(),
        //     target_nodes.len(),
        //     "The number of source and target nodes should be the same{}",
        //     self.debug_dot()
        // );

        let source_node_vec = source_nodes.into_iter().collect_vec();
        let target_node_vec = target_nodes.into_iter().collect_vec();

        //println!("source nodes: {:?}", source_node_vec);
        //println!("target nodes: {:?}", target_node_vec);
        //panic!("stop");

        let source_node = self
            .underlying
            .combine_to_single_hedgenode(&source_node_vec);

        let target_node = self
            .underlying
            .combine_to_single_hedgenode(&target_node_vec);

        (source_node, target_node)
    }

    pub(crate) fn edge_name_to_index(&self, name: &str) -> Option<EdgeIndex> {
        for (_, edge_index, edge_data) in self.underlying.iter_edges() {
            if edge_data.data.name == name {
                return Some(edge_index);
            }
        }

        None
    }

    pub(crate) fn get_initial_state_tree(&self) -> (SuBitGraph, Vec<EdgeIndex>) {
        let mut tree_like_edges = Vec::new();
        let full_graph = self.underlying.full_filter();
        let full_is_cut = self
            .initial_state_cut
            .left
            .union(&self.initial_state_cut.right);

        let full_graph_without_initial_state_cut = full_graph.subtract(&full_is_cut);

        let mut result: SubSet<Hedge> = self.underlying.empty_subgraph();

        for (pair, edge_id, _) in self
            .underlying
            .iter_edges_of(&full_graph_without_initial_state_cut)
        {
            if let HedgePair::Paired { source, sink } = pair {
                let loop_signature = &self.loop_momentum_basis.edge_signatures[edge_id];
                let is_tree_like = loop_signature.internal.iter().all(|sign| sign.is_zero());
                if is_tree_like {
                    tree_like_edges.push(edge_id);
                    let source_node = self.underlying.node_id(source);
                    let sink_node = self.underlying.node_id(sink);

                    let source_connects_initial_state =
                        self.underlying.iter_crown(source_node).any(|hedge| {
                            let mut single_hedge_subgraph: SubSet<Hedge> =
                                self.underlying.empty_subgraph();

                            single_hedge_subgraph.add(hedge);

                            single_hedge_subgraph.intersects(&full_is_cut)
                        });

                    if source_connects_initial_state {
                        for hedge in self.underlying.iter_crown(source_node) {
                            result.add(hedge);
                        }
                    } else {
                        for hedge in self.underlying.iter_crown(sink_node) {
                            result.add(hedge);
                        }
                    }
                }
            }
        }

        (result, tree_like_edges)
    }

    pub(crate) fn get_raised_edge_groups(&self) -> Vec<Vec<EdgeIndex>> {
        let mut result = Vec::<Vec<EdgeIndex>>::new();

        for (_, edge_index, _) in self.iter_loop_edges() {
            let group_position = result.iter().position(|group| {
                group.iter().all(|e| {
                    self.loop_momentum_basis.edges_are_raised(*e, edge_index)
                        && self[edge_index].mass == self[*e].mass
                })
            });

            if let Some(pos) = group_position {
                result[pos].push(edge_index);
            } else {
                result.push(vec![edge_index]);
            }
        }

        result.iter_mut().for_each(|group| group.sort());

        result
    }

    pub fn get_edges_in_initial_state_cut(&self) -> Vec<EdgeIndex> {
        self.underlying
            .iter_edges_of(&self.initial_state_cut)
            .map(|(_, edge_index, _)| edge_index)
            .collect_vec()
    }
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

#[derive(
    Debug, Copy, Clone, PartialEq, Eq, bincode_trait_derive::Encode, bincode_trait_derive::Decode,
)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct ExternalConnection {
    pub incoming_index: ExternalIndex,
    pub outgoing_index: ExternalIndex,
}

pub fn get_cff_inverse_energy_product_impl<E, V, H, S: SubSetLike>(
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
