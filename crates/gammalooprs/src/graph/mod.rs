use std::ops::Index;

use ahash::AHashSet;
use bincode_trait_derive::{Decode, Encode};
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
use color_eyre::Result;
use eyre::eyre;
use spenso::algebra::complex::Complex;
use symbolica::atom::{Atom, AtomCore};
use three_dimensional_reps::ThreeDGraphSource;
use tracing::warn;
use typed_index_collections::TiVec;

use crate::{
    cff::surface::SurfaceCache,
    define_index,
    feyngen::diagram_generator::evaluate_overall_factor,
    graph::edge::EdgeMass,
    integrands::process::{LmbMultiChannelingSetup, ParamBuilder},
    model::Model,
    momentum::{Dep, ExternalMomenta, PolDef, sample::ExternalIndex},
    numerator::GlobalPrefactor,
    processes::DotExportSettings,
    settings::runtime::kinematic::{Externals, improvement::PhaseSpaceImprovementSettings},
    utils::{F, Length, ose_atom_from_index, symbolica_ext::LogPrint},
};

pub(crate) mod attribute_warnings;
pub mod autogen;
pub mod cuts;
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
    pub underlying: HedgeGraph<Edge, Vertex, HedgeData>,
    pub loop_momentum_basis: LoopMomentumBasis,
    pub param_builder: ParamBuilder,
    pub global_prefactor: GlobalPrefactor,
    pub surface_cache: SurfaceCache,
    /// The cross section initial state cut
    /// Only relevant for cross sections, but stored here for the parsing
    pub initial_state_cut: OrientedCut,
    pub polarizations: Vec<(PolDef, Atom)>,
}

#[derive(Clone, Debug)]
pub struct ThreeDRepMassShift {
    pub group_index: usize,
    pub split_index: usize,
    pub local_edge_id: usize,
    pub edge_id: EdgeIndex,
    pub base_mass: f64,
    pub shifted_mass: f64,
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

    pub fn debug_dot_with_settings(&self, settings: &DotExportSettings) -> String {
        self.to_dot_graph_with_settings(settings).debug_dot()
    }

    pub fn rebuild_param_builder(&mut self, model: &Model) {
        let additional_params = self.param_builder.pairs.additional_params.params.clone();
        let loop_momentum_basis = self.loop_momentum_basis.clone();
        self.param_builder =
            ParamBuilder::new(self, model, &loop_momentum_basis, additional_params);
    }

    pub fn split_repeated_masses_for_three_drep(
        &self,
        model: &Model,
        epsilon: f64,
    ) -> Result<(Self, Vec<ThreeDRepMassShift>)> {
        let parsed = self
            .to_three_d_parsed_graph()
            .map_err(|error| eyre!("could not inspect repeated 3Drep propagators: {error}"))?;
        let repeated_groups = three_dimensional_reps::repeated_groups(&parsed);
        let source_edges = self
            .underlying
            .iter_edges()
            .filter(|(pair, _, _)| matches!(pair, HedgePair::Paired { .. }))
            .map(|(_, edge, _)| edge)
            .sorted()
            .collect::<Vec<_>>();

        let mut shifted = self.clone();
        let mut records = Vec::new();
        for (group_index, group) in repeated_groups.iter().enumerate() {
            let center = (group.edge_ids.len().saturating_sub(1)) as f64 / 2.0;
            for (split_index, local_edge_id) in group.edge_ids.iter().copied().enumerate() {
                let Some(edge_id) = source_edges.get(local_edge_id).copied() else {
                    return Err(eyre!(
                        "repeated 3Drep edge id {local_edge_id} does not map to a GammaLoop edge"
                    ));
                };
                let base_mass = self.underlying[edge_id]
                    .mass
                    .value::<f64>(model, &self.param_builder)
                    .map(|value| value.re.0)
                    .ok_or_else(|| {
                        eyre!(
                            "could not evaluate mass for repeated 3Drep edge {} (group {group_index}, split {split_index}, local edge id {local_edge_id})",
                            edge_id.0
                        )
                    })?;
                let shifted_mass = base_mass + (split_index as f64 - center) * epsilon;
                shifted.underlying[edge_id].particle = shifted.underlying[edge_id]
                    .particle
                    .clone()
                    .override_mass(Some(Atom::num(shifted_mass)));
                shifted.underlying[edge_id].mass =
                    EdgeMass::Value(Complex::new_re(F(shifted_mass)));
                records.push(ThreeDRepMassShift {
                    group_index,
                    split_index,
                    local_edge_id,
                    edge_id,
                    base_mass,
                    shifted_mass,
                });
            }
        }
        shifted.rebuild_param_builder(model);
        Ok((shifted, records))
    }

    pub fn pretty_dot(&self) -> String {
        self.dot_impl(
            &self.full_filter(),
            "",
            &|_a| None,
            &|e| Some(format!("label=\"{}\"", e.num.log_print(None))),
            &|v| Some(format!("label=\"{}\"", v.num.log_print(None))),
        )
    }

    pub(crate) fn global_atom(&self) -> Atom {
        &self.global_prefactor.num
            * &self.global_prefactor.projector
            * evaluate_overall_factor(self.overall_factor.as_view())
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
        override_lmb_heuristics: bool,
    ) -> LmbMultiChannelingSetup {
        let mut channels: Vec<LmbIndex> = Vec::new();

        if override_lmb_heuristics {
            channels.extend(lmbs.iter_enumerated().map(|(lmb_index, _)| lmb_index));

            return LmbMultiChannelingSetup {
                channels: channels.into_iter().sorted().collect(),
                graph: self.clone(),
                all_bases: lmbs.clone(),
            };
        }

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

    pub(crate) fn get_edge_subgraph(&self, edge: EdgeIndex) -> SuBitGraph {
        let mut subgraph: SuBitGraph = self.underlying.empty_subgraph();

        match self[&edge].1 {
            HedgePair::Paired { source, sink } => {
                subgraph.add(source);
                subgraph.add(sink);
            }
            HedgePair::Unpaired { hedge, .. } => {
                subgraph.add(hedge);
            }
            HedgePair::Split { .. } => unreachable!(),
        }

        subgraph
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
            if edge_data.data.name.value == name {
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
    pub(crate) fn is_always_pinch(
        &self,
        sandwich: &SuBitGraph,
        cut_a: &OrientedCut,
        cut_b: &OrientedCut,
    ) -> bool {
        let cut_a_all_edges = cut_a.left.union(&cut_a.right);
        let cut_b_all_edges = cut_b.left.union(&cut_b.right);

        let cut_a_sandwich = cut_a_all_edges.intersection(sandwich);
        let cut_b_sandwich = cut_b_all_edges.intersection(sandwich);

        let n_edges_cut_a_sandwich = cut_a_sandwich.n_included();
        let n_edges_cut_b_sandwich = cut_b_sandwich.n_included();

        if n_edges_cut_a_sandwich > 1 && n_edges_cut_b_sandwich > 1 {
            return false;
        }

        let mut cut_a_mass_sum = Atom::new();
        for (_, _, edge_data) in self.iter_edges_of(&cut_a_sandwich) {
            match edge_data.data.mass {
                EdgeMass::Zero => {}
                EdgeMass::ModelVar(m) => {
                    cut_a_mass_sum += m;
                }
                _ => {
                    // If there is a non-zero, non-model-var mass, we can't make any statements about pinches without runtime information, so we return false
                    return false;
                }
            }
        }

        let mut cut_b_mass_sum = Atom::new();
        for (_, _, edge_data) in self.iter_edges_of(&cut_b_sandwich) {
            match edge_data.data.mass {
                EdgeMass::Zero => {}
                EdgeMass::ModelVar(m) => {
                    cut_b_mass_sum += m;
                }
                _ => {
                    // If there is a non-zero, non-model-var mass, we can't make any statements about pinches without runtime information, so we return false
                    return false;
                }
            }
        }

        (cut_a_mass_sum - cut_b_mass_sum).expand().is_zero()
    }
}

pub mod edge;
pub mod parse;
pub use autogen::Autogen;
pub use edge::Edge;
pub mod hedge_data;
pub use hedge_data::HedgeData;
pub(crate) mod three_d_source;
pub(crate) use three_d_source::GraphThreeDSource;
pub mod vertex;
pub use vertex::Vertex;

pub mod lmb;
pub use lmb::{LMBext, LmbError, LmbIndex, LoopMomentumBasis};

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
