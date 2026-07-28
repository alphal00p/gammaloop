use std::{
    collections::{BTreeMap, BTreeSet},
    ops::Index,
};

use ahash::AHashSet;
use bincode_trait_derive::{Decode, Encode};
use gammaloop_tracing_filter::LogMessage;
use itertools::Itertools;
use linnet::{
    half_edge::{
        HedgeGraph,
        involution::{EdgeData, EdgeIndex, Hedge, HedgePair},
        subgraph::{
            HedgeNode, Inclusion, ModifySubSet, OrientedCut, SuBitGraph, SubGraphLike, SubSetLike,
            SubSetOps, subset::SubSet,
        },
    },
    parser::DotGraph,
};
use tracing::debug;

use rand::{Rng, SeedableRng, rngs::SmallRng};
// use petgraph::Direction::Outgoing;
use symbolica::atom::{Atom, AtomCore};
use tracing::warn;
use typed_index_collections::TiVec;

use crate::{
    cff::generation::SurfaceCache,
    define_index,
    feyngen::diagram_generator::evaluate_overall_factor,
    integrands::process::{ChannelIndex, LmbMultiChannelingSetup, ParamBuilder},
    momentum::{Dep, ExternalMomenta, PolDef, sample::ExternalIndex},
    numerator::GlobalPrefactor,
    processes::DotExportSettings,
    settings::runtime::kinematic::{Externals, improvement::PhaseSpaceImprovementSettings},
    utils::{F, Length, ose_atom_from_index, symbolica_ext::LogPrint},
    uv::uv_graph::UVE,
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

impl LogMessage for Graph {
    fn log_display(&self) -> String {
        self.name.to_string()
    }
}

// impl Deref for Graph {
//     type Target = HedgeGraph<Edge, Vertex>;
// }

pub mod feynman_graph;
pub use feynman_graph::FeynmanGraph;
pub mod ext;

#[derive(Clone, Copy)]
pub(crate) enum LmbChannelFallback {
    CurrentGraphBasis,
    FirstBasis,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum ThresholdPinchStatus {
    Always,
    CanBecome,
    NotProven,
}

impl Graph {
    pub fn debug_dot(&self) -> String {
        DotGraph::from(self).debug_dot()
    }

    pub fn debug_dot_with_settings(&self, settings: &DotExportSettings) -> String {
        self.to_dot_graph_with_settings(settings).debug_dot()
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

    pub(crate) fn external_momentum_edge_order(&self) -> Vec<EdgeIndex> {
        if self.initial_state_cut.nedges(&self.underlying) == 0 {
            let external_filter: SuBitGraph = self.external_filter();
            external_filter
                .included_iter()
                .sorted()
                .map(|hedge| self.underlying[&hedge])
                .collect_vec()
        } else {
            self.initial_state_cut
                .iter_left_hedges()
                .sorted()
                .map(|hedge| self.underlying[&hedge])
                .collect_vec()
        }
    }

    pub(crate) fn canonicalize_lmb_external_order(&self, lmb: &mut LoopMomentumBasis) {
        lmb.canonicalize_external_order(&self.external_momentum_edge_order());
    }

    pub(crate) fn dummy_stripped_external_flows_of<S: SubGraphLike>(&self, subgraph: &S) -> S::Base
    where
        S::Base: ModifySubSet<HedgePair> + ModifySubSet<Hedge>,
    {
        let mut externals = self.underlying.full_crown(subgraph);
        for (pair, _, edge) in self.underlying.iter_edges() {
            if edge.data.is_dummy {
                externals.sub(pair);
            }
        }
        externals
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
        let channels = self.select_amplitude_lmb_channel_indices(
            lmbs,
            override_lmb_heuristics,
            LmbChannelFallback::CurrentGraphBasis,
        );

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

    pub(crate) fn select_amplitude_lmb_channel_indices(
        &self,
        lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        override_lmb_heuristics: bool,
        fallback: LmbChannelFallback,
    ) -> TiVec<ChannelIndex, LmbIndex> {
        if override_lmb_heuristics {
            return lmbs
                .iter_enumerated()
                .map(|(lmb_index, _)| lmb_index)
                .collect();
        }

        let Some((_, first_lmb)) = lmbs.iter_enumerated().next() else {
            return TiVec::new();
        };
        let num_loops = first_lmb.loop_edges.len();
        if num_loops == 0 {
            return self.fallback_lmb_channel_indices(lmbs, fallback);
        }

        let mut universe = BTreeMap::<Vec<EdgeIndex>, usize>::new();
        let mut candidate_covers = Vec::<(LmbIndex, BTreeSet<usize>)>::new();
        for (lmb_index, lmb) in lmbs.iter_enumerated() {
            let mut cover = BTreeSet::<usize>::new();
            let massless_loop_edges = lmb
                .loop_edges
                .iter()
                .copied()
                .filter(|edge_id| self.underlying[*edge_id].particle.is_massless())
                .collect_vec();

            for mut combination in massless_loop_edges.into_iter().combinations(num_loops) {
                combination.sort();
                let next_universe_id = universe.len();
                let universe_id = *universe.entry(combination).or_insert(next_universe_id);
                cover.insert(universe_id);
            }

            if !cover.is_empty() {
                candidate_covers.push((lmb_index, cover));
            }
        }

        if universe.is_empty() {
            return self.fallback_lmb_channel_indices(lmbs, fallback);
        }

        let channels = Self::minimum_lmb_set_cover(universe.len(), &candidate_covers);
        if channels.is_empty() {
            self.fallback_lmb_channel_indices(lmbs, fallback)
        } else {
            channels.into_iter().sorted().collect()
        }
    }

    fn fallback_lmb_channel_indices(
        &self,
        lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
        fallback: LmbChannelFallback,
    ) -> TiVec<ChannelIndex, LmbIndex> {
        let fallback_index = match fallback {
            LmbChannelFallback::CurrentGraphBasis => lmbs
                .iter_enumerated()
                .find(|(_, lmb)| lmb.loop_edges == self.loop_momentum_basis.loop_edges)
                .map(|(lmb_index, _)| lmb_index)
                .or_else(|| {
                    warn!(
                        "lmb not in the list of all lmb, falling back to the first generated basis"
                    );
                    lmbs.iter_enumerated()
                        .next()
                        .map(|(lmb_index, _)| lmb_index)
                }),
            LmbChannelFallback::FirstBasis => lmbs
                .iter_enumerated()
                .next()
                .map(|(lmb_index, _)| lmb_index),
        };

        fallback_index.into_iter().collect()
    }

    fn minimum_lmb_set_cover(
        universe_len: usize,
        candidate_covers: &[(LmbIndex, BTreeSet<usize>)],
    ) -> Vec<LmbIndex> {
        let mut candidates_by_element = vec![Vec::<usize>::new(); universe_len];
        for (candidate_pos, (_, cover)) in candidate_covers.iter().enumerate() {
            for element in cover {
                candidates_by_element[*element].push(candidate_pos);
            }
        }

        let remaining = (0..universe_len).collect::<BTreeSet<_>>();
        let mut chosen = Vec::<usize>::new();
        let mut best = None;
        Self::search_minimum_lmb_set_cover(
            remaining,
            &mut chosen,
            &mut best,
            candidate_covers,
            &candidates_by_element,
        );
        best.unwrap_or_default()
    }

    fn search_minimum_lmb_set_cover(
        remaining: BTreeSet<usize>,
        chosen: &mut Vec<usize>,
        best: &mut Option<Vec<LmbIndex>>,
        candidate_covers: &[(LmbIndex, BTreeSet<usize>)],
        candidates_by_element: &[Vec<usize>],
    ) {
        if remaining.is_empty() {
            let mut candidate = chosen
                .iter()
                .map(|candidate_pos| candidate_covers[*candidate_pos].0)
                .collect_vec();
            candidate.sort();
            let better = best.as_ref().is_none_or(|current_best| {
                candidate.len() < current_best.len()
                    || (candidate.len() == current_best.len() && candidate < *current_best)
            });
            if better {
                *best = Some(candidate);
            }
            return;
        }

        if let Some(current_best) = best.as_ref()
            && chosen.len() >= current_best.len()
        {
            return;
        }

        let Some(element) = remaining
            .iter()
            .copied()
            .min_by_key(|element| candidates_by_element[*element].len())
        else {
            return;
        };
        if candidates_by_element[element].is_empty() {
            return;
        }

        for candidate_pos in &candidates_by_element[element] {
            if chosen.contains(candidate_pos) {
                continue;
            }
            let cover = &candidate_covers[*candidate_pos].1;
            if cover.is_disjoint(&remaining) {
                continue;
            }

            let mut next_remaining = remaining.clone();
            for covered_element in cover {
                next_remaining.remove(covered_element);
            }
            chosen.push(*candidate_pos);
            Self::search_minimum_lmb_set_cover(
                next_remaining,
                chosen,
                best,
                candidate_covers,
                candidates_by_element,
            );
            chosen.pop();
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
        self.initial_state_cut
            .iter_left_hedges()
            .map(|hedge| self.underlying[&hedge])
            .collect_vec()
    }
    pub(crate) fn classify_threshold_pinch(
        &self,
        cut_boundary_edges: &[EdgeIndex],
        threshold_boundary_edges: &[EdgeIndex],
    ) -> ThresholdPinchStatus {
        // The caller has already intersected both cut boundaries with the connected sandwich.
        // These edge sets therefore carry the same graph-relative information that the former
        // `is_always_pinch` implementation derived from the two cuts and sandwich internally.
        // Keep every mass representation symbolic here. A non-zero difference is not a
        // model-independent conclusion; resolved values handle that case during generation.
        let boundary_mass_sum = |edges: &[EdgeIndex]| {
            edges
                .iter()
                .map(|edge_id| self[*edge_id].mass_atom())
                .fold(Atom::new(), |sum, mass| sum + mass)
        };

        let mass_sums_are_identical = (boundary_mass_sum(cut_boundary_edges)
            - boundary_mass_sum(threshold_boundary_edges))
        .expand()
        .is_zero();

        if !mass_sums_are_identical {
            return ThresholdPinchStatus::NotProven;
        }

        if cut_boundary_edges.len() > 1 && threshold_boundary_edges.len() > 1 {
            ThresholdPinchStatus::CanBecome
        } else {
            ThresholdPinchStatus::Always
        }
    }
}

pub mod edge;
pub mod parse;
pub use autogen::Autogen;
pub use edge::Edge;
pub mod hedge_data;
pub use hedge_data::HedgeData;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        dot, graph::parse::from_dot::IntoGraph, initialisation::test_initialise,
        momentum::signature::LoopExtSignature,
    };
    use std::sync::OnceLock;

    fn test_lmb(graph: &Graph, loop_edges: &[usize]) -> LoopMomentumBasis {
        LoopMomentumBasis {
            tree: graph.underlying.empty_subgraph(),
            loop_edges: loop_edges.iter().copied().map(EdgeIndex::from).collect(),
            ext_edges: Vec::new().into(),
            edge_signatures: graph.underlying.new_edgevec(|_, _, _| {
                LoopExtSignature::from((Vec::<isize>::new(), Vec::<isize>::new()))
            }),
        }
    }

    fn selector_test_graph() -> Graph {
        static GRAPH: OnceLock<Graph> = OnceLock::new();
        GRAPH
            .get_or_init(|| {
                test_initialise().unwrap();
                dot!(
                    digraph lmb_selector {
                        edge [num=1 mass=0]
                        node [num=1]
                        A -> B [id=0]
                        A -> B [id=1]
                        A -> B [id=2]
                        A -> B [id=3]
                        A -> B [id=4 mass=1]
                    }
                )
                .unwrap()
            })
            .clone()
    }

    #[test]
    fn exact_lmb_set_cover_prefers_smallest_deterministic_solution() {
        let candidates = [
            (LmbIndex::from(0), BTreeSet::from([0])),
            (LmbIndex::from(1), BTreeSet::from([1])),
            (LmbIndex::from(2), BTreeSet::from([0, 1])),
            (LmbIndex::from(3), BTreeSet::from([0, 1])),
        ];

        assert_eq!(
            Graph::minimum_lmb_set_cover(2, &candidates),
            vec![LmbIndex::from(2)]
        );
    }

    #[test]
    fn amplitude_lmb_selector_covers_massless_combinations_and_skips_duplicates() {
        let mut graph = selector_test_graph();
        graph.loop_momentum_basis = test_lmb(&graph, &[0, 1]);
        let lmbs: TiVec<LmbIndex, LoopMomentumBasis> = vec![
            test_lmb(&graph, &[0, 1]),
            test_lmb(&graph, &[2, 3]),
            test_lmb(&graph, &[0, 2]),
            test_lmb(&graph, &[0, 1]),
            test_lmb(&graph, &[0, 4]),
        ]
        .into();

        let selected = graph.select_amplitude_lmb_channel_indices(
            &lmbs,
            false,
            LmbChannelFallback::CurrentGraphBasis,
        );

        assert_eq!(
            selected.into_iter().collect_vec(),
            vec![LmbIndex::from(0), LmbIndex::from(1), LmbIndex::from(2)]
        );
    }

    #[test]
    fn amplitude_lmb_selector_override_keeps_all_bases() {
        let graph = selector_test_graph();
        let lmbs: TiVec<LmbIndex, LoopMomentumBasis> = vec![
            test_lmb(&graph, &[0, 1]),
            test_lmb(&graph, &[2, 3]),
            test_lmb(&graph, &[0, 4]),
        ]
        .into();

        let selected = graph.select_amplitude_lmb_channel_indices(
            &lmbs,
            true,
            LmbChannelFallback::CurrentGraphBasis,
        );

        assert_eq!(
            selected.into_iter().collect_vec(),
            vec![LmbIndex::from(0), LmbIndex::from(1), LmbIndex::from(2)]
        );
    }

    #[test]
    fn amplitude_lmb_selector_falls_back_when_no_massless_combination_exists() {
        let mut graph = selector_test_graph();
        graph.loop_momentum_basis = test_lmb(&graph, &[0, 4]);
        let lmbs: TiVec<LmbIndex, LoopMomentumBasis> =
            vec![test_lmb(&graph, &[1, 4]), test_lmb(&graph, &[0, 4])].into();

        let selected_current = graph.select_amplitude_lmb_channel_indices(
            &lmbs,
            false,
            LmbChannelFallback::CurrentGraphBasis,
        );
        let selected_first = graph.select_amplitude_lmb_channel_indices(
            &lmbs,
            false,
            LmbChannelFallback::FirstBasis,
        );

        assert_eq!(
            selected_current.into_iter().collect_vec(),
            vec![LmbIndex::from(1)]
        );
        assert_eq!(
            selected_first.into_iter().collect_vec(),
            vec![LmbIndex::from(0)]
        );
    }

    #[test]
    fn threshold_pinch_classification_distinguishes_fixed_and_multiparticle_boundaries() {
        let graph = selector_test_graph();

        assert_eq!(
            graph.classify_threshold_pinch(&[EdgeIndex::from(0)], &[EdgeIndex::from(1)],),
            ThresholdPinchStatus::Always,
        );
        assert_eq!(
            graph.classify_threshold_pinch(
                &[EdgeIndex::from(0), EdgeIndex::from(1)],
                &[EdgeIndex::from(2), EdgeIndex::from(3)],
            ),
            ThresholdPinchStatus::CanBecome,
        );
        assert_eq!(
            graph.classify_threshold_pinch(&[EdgeIndex::from(0)], &[EdgeIndex::from(4)],),
            ThresholdPinchStatus::NotProven,
        );
    }
}
