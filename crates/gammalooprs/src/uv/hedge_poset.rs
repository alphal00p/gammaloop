use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Display,
    sync::Arc,
};

#[cfg(test)]
use crate::uv::Integrands;
#[cfg(test)]
use crate::{
    cff::expression::OrientationID,
    uv::approx::direct_3d::{DirectCoordinateFrame, DirectResidueBranches, DirectSector},
};
#[cfg(test)]
use std::cmp::Reverse;

use std::collections::btree_map::Entry;

use ahash::AHashMap;
use eyre::{WrapErr, ensure, eyre};
use gammaloop_tracing_filter::LogMessage;
use idenso::{color::ColorSimplifier, shorthands::schoonschip::Schoonschip};
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph, NoData, NodeIndex,
    algorithms::{
        DirectionBasis,
        trace_unfold::{HiddenData, Independence, TraceKey, TraceUnfold, UnfoldedTraceGraph},
    },
    involution::{EdgeIndex, Flow, HedgePair},
    nodestore::{NodeStorageOps, NodeStorageVec},
    subgraph::{Inclusion, InternalSubGraph, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    function,
};
use tracing::debug;
use vakint::Vakint;

use crate::uv::spinney::WoodRole;
use crate::uv::thermal_zero_mode::{ThermalZeroModeRegion, ZeroModeHostError, ZeroModeHostKind};
use crate::{
    debug_tags,
    graph::{
        Graph, LMBext, LoopMomentumBasis,
        cuts::CutSet,
        parse::string_utils::{ToOrderedSimple, dot_attr_value},
    },
    settings::global::GenerationSettings,
    utils::{GS, W_},
    uv::{
        ApproximationType, RenormalizationPart, Spinney, UVgenerationSettings, UltravioletGraph,
        approx::{
            CutStructure, ForestNodeLike, OrientationProjection, Rooted, UVCtx,
            direct_3d::{Direct3dApproximation, Direct3dCts},
            final_integrand::{FinalIntegrandBuilder, FinalIntegrands},
            integrated::{Integrated, IntegratedCts},
            local_3d::{Local3DCts, Localizer},
            local_4d::{self, Full4dCts, Local4dCts},
            projected_4d::{Local4dProjectionContext, Projected4dApproximation},
        },
        export::{
            UVForestLocalProvenance, UVForestNodeExpression, UVForestNodeProvenance,
            UVForestProvenance,
        },
        forest::ParametricIntegrands,
        marker::UvMarker,
    },
};
use color_eyre::Result;
use spenso::shadowing::symbolica_utils::{LogPrint, SpensoPrintSettings};

pub struct Wood {
    pub graph: HedgeGraph<SuBitGraph, Spinney>,
    pub root: NodeIndex,
    pub vakint_settings: vakint::VakintSettings,
    cuts: CutStructure,
}

impl Independence<HiddenData<SuBitGraph, EdgeIndex>> for Wood {
    fn independent(
        &self,
        a: &HiddenData<SuBitGraph, EdgeIndex>,
        b: &HiddenData<SuBitGraph, EdgeIndex>,
    ) -> bool {
        !a.order.intersects(&b.order)
    }
}

impl TraceUnfold<SuBitGraph> for Wood {
    type EdgeData = SuBitGraph;
    type HedgeData = NoData;
    type NodeData = Spinney;
    type NodeStorage = NodeStorageVec<Spinney>;

    fn graph(
        &self,
    ) -> &HedgeGraph<Self::EdgeData, Self::NodeData, Self::HedgeData, Self::NodeStorage> {
        &self.graph
    }

    fn key(&self, e: EdgeIndex) -> SuBitGraph {
        self.graph[e].clone()
    }

    /// Treats a wood node as a factorized join when its incoming sink edges are exactly the
    /// disjoint connected components of the target spinney.
    ///
    /// The returned `SuBitGraph`s are the required branch factors for the generic trace unfold:
    /// their union must be the target filter, they must be pairwise disjoint, and there must be
    /// one factor per connected component.
    fn join_factors(&self, target: NodeIndex) -> Option<BTreeSet<SuBitGraph>> {
        if self.graph[target].n_components() < 2 {
            return None;
        }

        let factors = self
            .graph
            .iter_crown(target)
            .filter(|hedge| self.graph.flow(*hedge) == Flow::Sink)
            .map(|hedge| self.graph[self.graph[&hedge]].clone())
            .collect::<BTreeSet<_>>();
        if factors.len() != self.graph[target].n_components() {
            return None;
        }
        if factors.iter().any(|factor| {
            !self
                .graph
                .iter_nodes()
                .any(|(_, _, spinney)| spinney.n_components() == 1 && spinney.filter() == factor)
        }) {
            return None;
        }

        let mut cover: Option<SuBitGraph> = None;
        for factor in &factors {
            if let Some(acc) = &mut cover {
                if acc.intersects(factor) {
                    return None;
                }
                acc.union_with(factor);
            } else {
                cover = Some(factor.clone());
            }
        }

        (cover.as_ref() == Some(self.graph[target].filter())).then_some(factors)
    }
}

impl Wood {
    pub fn current_given_pair<'a>(
        &'a self,
        edge_id: EdgeIndex,
        order: usize,
    ) -> (ForestNode<'a>, ForestNode<'a>) {
        let HedgePair::Paired { source, sink } = self.graph[&edge_id].1 else {
            panic!("edge in self is not paired");
        };

        // get this hedge's forest node from the self. This is the node that has already been computed (as it is a parent to this edge)
        let given = ForestNode {
            spinney: &self.graph[self.graph.node_id(source)],
            topo_order: order,
        };

        // this is the current node, which should be the same for all union edges (since they all have the same sink)
        let current_for_h = self.graph.node_id(sink);

        // this is the current node, which we want to compute with
        let current = ForestNode {
            spinney: &self.graph[current_for_h],
            topo_order: order,
        };

        (current, given)
    }

    pub(crate) fn new(
        cuts: CutStructure,
        graph: &Graph,
        settings: &GenerationSettings,
    ) -> Result<Self> {
        let mut spinneys = Vec::new();
        for cut_spinneys in graph.classified_cut_spinneys(&cuts, &settings.uv, &settings.medium) {
            spinneys.extend(cut_spinneys?);
        }

        Self::from_spinneys(spinneys, graph, cuts, &settings.uv).wrap_err_with(|| {
            format!(
                "graph '{}': failed to build the hedge-poset wood",
                graph.name
            )
        })
    }

    /// Builds the wood of `s`, one node per subgraph and [`WoodRole`]: the
    /// first spinney of each subgraph and role is kept. A zero-mode union host
    /// is a second node on the subgraph of its plain factorized join, above
    /// it, and one collective operation.
    ///
    /// # Errors
    ///
    /// Duplicates of one node with different zero-mode regions, and a union
    /// host whose scope equals that of a collective plain union.
    pub(crate) fn from_spinneys<I: IntoIterator<Item = Spinney>>(
        s: I,
        graph: &Graph,
        cuts: CutStructure,
        settings: &UVgenerationSettings,
    ) -> Result<Self> {
        let spinneys = Self::nodes(s, graph)?;
        let max_loops = spinneys
            .values()
            .map(|spinney| graph.n_loops(spinney.filter()))
            .max()
            .unwrap_or(0);
        let mut vakint_settings = settings.vakint.true_settings();
        vakint_settings.project_onto_tensor_integrals =
            settings.project_integrated_uv_cts_onto_tensor_integrals;
        // Retain enough positive epsilon powers for finite terms formed when disconnected
        // integrated counterterms are multiplied.
        vakint_settings.number_of_terms_in_epsilon_expansion = max_loops as i64 + 1;

        let mut unions = BTreeSet::new();
        let g: HedgeGraph<_, _> = HedgeGraph::poset(spinneys.into_values());
        let mut poset = g.map(
            |_, _, v| v,
            |_, n, pair, _, e| {
                let HedgePair::Paired { source, sink } = pair else {
                    return e.map(|_| graph.as_ref().empty_subgraph());
                };
                let nsource = n.node_id_ref(source);
                let nsink = n.node_id_ref(sink);

                let source_subgraph = &n.get_node_data(nsource).subgraph;
                let sink_subgraph = &n.get_node_data(nsink).subgraph;
                let reduced_subgraph = sink_subgraph.subtract(source_subgraph).filter;
                let hairy_source = graph.as_ref().full_crown(source_subgraph);

                if graph.as_ref().bridges_of(&reduced_subgraph).is_empty()
                    && !hairy_source.intersects(&reduced_subgraph)
                {
                    // if the reduced graph is bridgless,and has no node overlap with the source, then the source subgraph is cycle independent of reduced subgraph
                    // if the source is not empty then this is a disjoint union
                    if !source_subgraph.is_empty() {
                        unions.insert(nsink);
                    }
                    e.map(|_| reduced_subgraph)
                } else {
                    e.map(|_| sink_subgraph.filter.clone())
                }
            },
            |_, d| d,
        );

        let mut to_remove: SuBitGraph = poset.empty_subgraph();
        let classified_filters = poset
            .iter_nodes()
            .map(|(_, _, spinney)| spinney.filter().clone())
            .collect::<BTreeSet<_>>();
        let mut factorized_scopes = BTreeSet::new();

        // Not quite transitive closure. A disconnected target factorizes only when every
        // connected component is both an independently classified UV spinney and supplied by
        // an incoming edge. In that case, keep only those component edges so trace unfolding
        // builds the canonical product. If either condition fails, the target is a genuinely
        // collective UV region whose atomic Taylor term and nested forest terms are both
        // required. Label every retained transition by that full Taylor scope: the reduced
        // complement remains available from the edge's source/sink nodes, while using it as
        // the trace order would incorrectly make a collective operation commute with its
        // nested prefix and could coalesce distinct forests into a false factorized join.
        // A zero-mode union host is always collective: its operation S0 acts on its whole
        // subgraph. The transitions from its factors and its plain twin make it a union, so
        // every transition into it, from those, the nested prefixes of its factors and the
        // root, gets its scope here; the transition from the root already carries it.
        for &u in &unions {
            let collective_scope = poset[u].subgraph.filter.clone();
            let comps: BTreeSet<_> = graph
                .as_ref()
                .connected_components(&collective_scope)
                .into_iter()
                .collect();
            let incoming = poset
                .iter_crown(u)
                .filter(|c| poset.flow(*c) == Flow::Sink)
                .collect::<Vec<_>>();
            let supplied = incoming
                .iter()
                .map(|c| &poset[poset[c]])
                .filter(|factor| comps.contains(*factor))
                .cloned()
                .collect::<BTreeSet<_>>();
            let independently_classified = comps
                .iter()
                .all(|component| classified_filters.contains(component));
            let factorized = poset[u].wood_role() == WoodRole::Plain
                && independently_classified
                && supplied == comps;
            if factorized {
                for c in incoming {
                    if comps.contains(&poset[poset[&c]]) {
                        continue;
                    }
                    to_remove.add(c);
                    to_remove.add(poset.inv(c));
                }
                factorized_scopes.insert(collective_scope);
            } else {
                for c in incoming {
                    let edge = poset[&c];
                    poset[edge] = collective_scope.clone();
                }
            }
        }

        Self::check_union_host_scopes(&poset, &unions, &factorized_scopes, graph)?;

        poset.delete_hedges(&to_remove);
        let root = poset
            .iter_nodes()
            .find(|(_, _, s)| s.subgraph.is_empty())
            .map(|(n, _, _)| n)
            .ok_or_else(|| eyre!("internal error: the hedge-poset wood has no empty spinney"))?;

        Ok(Wood {
            graph: poset,
            root,
            cuts,
            vakint_settings,
        })
    }

    /// The nodes of the wood of `spinneys`, with the empty spinney, keyed by
    /// subgraph and [`WoodRole`]: the first spinney of each key is kept.
    ///
    /// # Errors
    ///
    /// Two spinneys of one key with different zero-mode regions, an internal
    /// error.
    fn nodes(
        spinneys: impl IntoIterator<Item = Spinney>,
        graph: &Graph,
    ) -> Result<BTreeMap<(SuBitGraph, WoodRole), Spinney>> {
        let mut nodes = BTreeMap::new();
        for spinney in spinneys {
            match nodes.entry((spinney.filter().clone(), spinney.wood_role())) {
                Entry::Vacant(entry) => {
                    entry.insert(spinney);
                }
                Entry::Occupied(entry) => ensure!(
                    entry.get().thermal_zero_mode() == spinney.thermal_zero_mode(),
                    "internal error: the spinney {} occurs with different zero-mode regions",
                    spinney.filter().string_label()
                ),
            }
        }
        let empty = Spinney::empty(graph);
        nodes
            .entry((empty.filter().clone(), empty.wood_role()))
            .or_insert(empty);
        Ok(nodes)
    }

    /// Trace operations are told apart by their scope alone, so a union host
    /// and a collective plain union of one subgraph would merge. Every union
    /// host of `poset` must be one of its `unions`, and its scope one of the
    /// `factorized_scopes` whenever a plain node shares it.
    ///
    /// # Errors
    ///
    /// A union host without a transition from its plain twin or a factor, an
    /// internal error, and a union host whose scope equals that of a
    /// collective plain union.
    fn check_union_host_scopes(
        poset: &HedgeGraph<SuBitGraph, Spinney>,
        unions: &BTreeSet<NodeIndex>,
        factorized_scopes: &BTreeSet<SuBitGraph>,
        graph: &Graph,
    ) -> Result<()> {
        for (host, _, spinney) in poset.iter_nodes() {
            if spinney.wood_role() != WoodRole::UnionHost {
                continue;
            }
            let scope = spinney.filter();
            ensure!(
                unions.contains(&host),
                "internal error: the zero-mode union host {} has no transition from its plain \
                 twin or a factor",
                scope.string_label()
            );
            let collective_twin = poset.iter_nodes().any(|(_, _, spinney)| {
                spinney.wood_role() == WoodRole::Plain && spinney.filter() == scope
            }) && !factorized_scopes.contains(scope);
            if collective_twin {
                return Err(ZeroModeHostError::CollectiveUnionScope {
                    host: graph.paired_edges(scope).into_iter().collect(),
                }
                .into());
            }
        }
        Ok(())
    }

    fn unfold_with_cached_node_label_atoms(self, cache_node_label_atoms: bool) -> Forests {
        let unfolded = self.trace_unfold::<NodeStorageVec<_>>(self.root);
        let graph = unfolded.map(|_, _, key| OperationNode { key });

        let mut cuts: Vec<(SuBitGraph, CutSet)> = Vec::new();
        for c in &self.cuts.cuts {
            let mut compatible: SuBitGraph = graph.empty_subgraph();
            for (_, crown, s) in graph.iter_nodes() {
                if s.is_compatible_with(c) {
                    for h in crown {
                        compatible.add(h);
                    }
                }
            }
            cuts.push((compatible, c.clone()));
        }
        let root = graph
            .iter_nodes()
            .find(|(_, _, operation)| operation.key.is_empty())
            .map(|(node, _, _)| node)
            .expect("no empty trace key found in unfolded hedge-poset forest");

        let forests = Forests {
            graph,
            cuts,
            root,
            cached_node_label_atoms: cache_node_label_atoms.then(AHashMap::new),
            compute_store: ComputeStore::default(),
            wood: self,
        };

        let forests = if cache_node_label_atoms {
            forests.with_cached_node_label_atoms()
        } else {
            forests
        };
        debug_tags!(#uv, #generation, #wood;
            trace_nodes = forests.graph.n_nodes(),
            file.wood = %forests.wood,
            file.forest = %forests,
            "Unfolded the hedge-poset wood"
        );
        forests
    }

    pub fn unfold(self) -> Forests {
        self.unfold_with_cached_node_label_atoms(true)
    }

    pub fn unfold_uncached(self) -> Forests {
        self.unfold_with_cached_node_label_atoms(false)
    }

    /// Whether some node of this wood hosts a thermal zero-mode region.
    pub(crate) fn hosts_zero_mode(&self) -> bool {
        self.graph
            .iter_nodes()
            .any(|(_, _, spinney)| spinney.thermal_zero_mode().is_some())
    }
}

impl Display for Wood {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.graph.dot_impl_fmt(
            f,
            &self.graph.full_filter(),
            "start=2;\n",
            &|_| None,
            &|a| Some(format!("label=\"{}\"", a.string_label())),
            &|v| {
                let role = match v.wood_role() {
                    WoodRole::Plain => "",
                    WoodRole::UnionHost => " (union host)",
                };
                Some(format!("label=\"{}{role}\"", v.subgraph.string_label()))
            },
        )
    }
}
// pub type SpinneyGraph = HedgeGraph<SuBitGraph, Spinney>;

// impl Key<&SuBitGraph> for SpinneyGraph {
//     fn key(&self, e: linnet::half_edge::involution::EdgeIndex) -> K {
//         &self[e]
//     }
// }

#[derive(Default)]
pub struct ComputeStore {
    entries: AHashMap<OperationNode, ComputeNode>,
}

impl ComputeStore {
    fn get(&self, key: &OperationNode) -> Option<&ComputeNode> {
        self.entries.get(key)
    }

    fn require(&self, key: &OperationNode) -> Result<&ComputeNode> {
        // Union operations, like individual Taylor steps, must be computed
        // before their typed values are consumed by another forest node.
        self.get(key)
            .ok_or_else(|| eyre!("{key} not yet added to compute store"))
    }

    fn entry(
        &mut self,
        key: OperationNode,
    ) -> std::collections::hash_map::Entry<'_, OperationNode, ComputeNode> {
        self.entries.entry(key)
    }
}

pub struct Forests {
    pub graph: UnfoldedTraceGraph<Wood, OperationNode, NodeStorageVec<OperationNode>>,
    pub root: NodeIndex,
    /// Wood subgraph that has compatible
    cuts: Vec<(SuBitGraph, CutSet)>,
    cached_node_label_atoms: Option<AHashMap<NodeIndex, Atom>>,
    pub compute_store: ComputeStore,
    wood: Wood,
}

#[derive(Clone, Debug, Hash, PartialEq, Eq)]
pub struct OperationNode {
    pub key: TraceKey<SuBitGraph, EdgeIndex>,
}

pub struct ForestNode<'a> {
    pub spinney: &'a Spinney,
    pub topo_order: usize,
}

pub struct OwnedForestNode {
    pub spinney: Spinney,
    pub topo_order: usize,
}

#[cfg(test)]
struct LocalLeafOperation {
    op: HiddenData<SuBitGraph, EdgeIndex>,
    frontier: NodeIndex,
}

struct UnionReplayState {
    integrated: NodeIndex,
    local_edges: Vec<EdgeIndex>,
}

#[cfg(test)]
impl LocalLeafOperation {
    fn new(op: &HiddenData<SuBitGraph, EdgeIndex>, frontier: NodeIndex) -> Self {
        Self {
            op: op.clone(),
            frontier,
        }
    }
}

impl OperationNode {
    pub fn is_compatible_with(&self, cut: &CutSet) -> bool {
        self.covers().is_none_or(|c| !c.intersects(&cut.union))
    }

    pub fn current<'a>(&'a self, wood: &'a Wood, topo_order: usize) -> Option<Vec<ForestNode<'a>>> {
        if self.key.is_empty() {
            return None;
        }

        Some(
            self.key
                .iter_leaf_ops()
                .map(|op| {
                    let HedgePair::Paired { sink, .. } = wood.graph[&op.data].1 else {
                        panic!("edge in trace key is not paired");
                    };
                    let spinney = &wood.graph[wood.graph.node_id(sink)];
                    ForestNode {
                        spinney,
                        topo_order,
                    }
                })
                .collect::<Vec<_>>(),
        )
    }
}

impl AsRef<TraceKey<SuBitGraph, EdgeIndex>> for OperationNode {
    fn as_ref(&self) -> &TraceKey<SuBitGraph, EdgeIndex> {
        &self.key
    }
}

impl LogMessage for ForestNode<'_> {
    fn log_display(&self) -> String {
        format!(
            "subgraph={}, topo_order={}, dod={}",
            self.spinney.filter().string_label(),
            self.topo_order,
            self.spinney.dod
        )
    }
}

impl ForestNodeLike for ForestNode<'_> {
    fn dod(&self) -> i32 {
        self.spinney.dod
    }
    fn renormalization_scheme(&self) -> crate::uv::ApproximationType {
        self.spinney.renormalization_scheme
    }
    fn lmb(&self) -> &LoopMomentumBasis {
        &self.spinney.lmb
    }
    fn reduced_subgraph(&self, given: &Self) -> SuBitGraph {
        self.spinney
            .subgraph
            .subtract(&given.spinney.subgraph)
            .filter
    }
    fn subgraph(&self) -> &SuBitGraph {
        self.spinney.filter()
    }
    fn topo_order(&self) -> usize {
        self.topo_order
    }
    fn thermal_zero_mode(&self) -> Option<&ThermalZeroModeRegion> {
        self.spinney.thermal_zero_mode()
    }
}

impl LogMessage for OwnedForestNode {
    fn log_display(&self) -> String {
        format!(
            "subgraph={}, topo_order={}, dod={}",
            self.spinney.filter().string_label(),
            self.topo_order,
            self.spinney.dod
        )
    }
}

impl ForestNodeLike for OwnedForestNode {
    fn dod(&self) -> i32 {
        self.spinney.dod
    }

    fn renormalization_scheme(&self) -> crate::uv::ApproximationType {
        self.spinney.renormalization_scheme
    }

    fn lmb(&self) -> &LoopMomentumBasis {
        &self.spinney.lmb
    }

    fn reduced_subgraph(&self, given: &Self) -> SuBitGraph {
        self.spinney
            .subgraph
            .subtract(&given.spinney.subgraph)
            .filter
    }

    fn subgraph(&self) -> &SuBitGraph {
        self.spinney.filter()
    }

    fn topo_order(&self) -> usize {
        self.topo_order
    }

    fn thermal_zero_mode(&self) -> Option<&ThermalZeroModeRegion> {
        self.spinney.thermal_zero_mode()
    }
}

impl Display for OperationNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.key.is_empty() {
            write!(f, "∅")
        } else {
            self.key.write_foata_like(f, |op| op.string_label())
        }
    }
}

impl OperationNode {
    fn foata_level_labels(&self) -> String {
        self.key
            .iter_levels_top_down()
            .map(|level| {
                level
                    .iter_leaf_ops()
                    .map(|op| op.order.string_label())
                    .join(",")
            })
            .join(";")
    }

    pub fn covers(&self) -> Option<SuBitGraph> {
        let mut acc: Option<SuBitGraph> = None;

        for level in self.key.iter_levels_top_down() {
            for op in level.iter_leaf_ops() {
                if let Some(a) = &mut acc {
                    a.union_with(&op.order);
                } else {
                    acc = Some(op.order.clone());
                }
            }
        }

        acc
    }

    /// The node of this operation's cover, for its integrated value. It is a
    /// plain spinney, so `source`, the wood node the operation ends on, must
    /// not host a zero-mode region that it would drop.
    fn forest_node(
        &self,
        graph: &Graph,
        source: &Spinney,
        topo_order: usize,
    ) -> Result<OwnedForestNode> {
        ensure!(
            source.thermal_zero_mode().is_none(),
            "graph '{}': operation {self} ends on the zero-mode host {}, whose region a forest \
             node of its cover would drop",
            graph.name,
            source.filter().string_label()
        );
        let spinney = match self.covers() {
            Some(cover) if !cover.is_empty() => {
                let subgraph = InternalSubGraph::cleaned_filter_optimist(cover, graph.as_ref());
                Spinney::new(subgraph, graph, &graph.loop_momentum_basis)
                    .wrap_err_with(|| {
                        format!(
                            "graph '{}': failed to construct forest node for operation {self}",
                            graph.name
                        )
                    })?
                    .ok_or_else(|| {
                        eyre!(
                            "graph '{}': operation {self} cover has negative degree of divergence",
                            graph.name
                        )
                    })?
            }
            _ => Spinney::empty(graph),
        };

        Ok(OwnedForestNode {
            spinney,
            topo_order,
        })
    }

    pub fn to_atom(&self) -> Atom {
        let mut acc = Atom::one();

        let approx = FunctionBuilder::new(GS.uv_approx);
        let mut levels = self.key.iter_levels_top_down();
        let Some(first_level) = levels.next() else {
            return acc;
        };

        let Some(first_op) = first_level.iter_leaf_ops().next() else {
            return acc;
        };

        let mut last = SuBitGraph::empty(first_op.order.size());
        for l in std::iter::once(first_level).chain(levels) {
            let last_sym = if last.is_empty() {
                Atom::Zero
            } else {
                last.symbol().to_atom()
            };

            let mut mul = Atom::one();

            for op in l.iter_leaf_ops() {
                let new = function!(op.order.symbol(), usize::from(op.data));
                mul *= approx.clone().add_arg((new - &last_sym) * &acc).finish();
                last.union_with(&op.order);
            }

            acc = mul
        }

        acc
    }

    // Four-dimensional and per-cut local terms are composed by `Forests` from typed
    // dependency-frontier values. Empty frontiers start from the typed roots, and
    // `Direct3dApproximation::run` applies the direct-3D subtraction signs.
}

#[derive(Default)]
pub struct ComputeNode {
    local_4d: Option<Local4dCts>,
    integrated: Option<IntegratedCts>,
    cuts: AHashMap<CutSet, CutComputation>,
}

pub struct CutComputation {
    local_3d: Local3DCts,
    final_integrands: FinalIntegrands,
}

impl ComputeNode {
    fn local_4d(&self, operation: &OperationNode) -> Result<&Local4dCts> {
        self.local_4d
            .as_ref()
            .ok_or_else(|| eyre!("{operation} has no computed local 4D counterterm"))
    }

    fn integrated(&self, operation: &OperationNode) -> Result<&IntegratedCts> {
        self.integrated
            .as_ref()
            .ok_or_else(|| eyre!("{operation} has no computed integrated counterterm"))
    }

    fn cut(&self, operation: &OperationNode, cutset: &CutSet) -> Result<&CutComputation> {
        self.cuts.get(cutset).ok_or_else(|| {
            eyre!("{operation} has no computed local counterterms for cut {cutset:?}")
        })
    }
}

impl Forests {
    fn source_spinney(&self, node: NodeIndex) -> &Spinney {
        &self.wood.graph[self.graph.source_node(node)]
    }

    fn recursion_input_4d(&self, node: NodeIndex) -> Result<Full4dCts> {
        // The typed 4D value is the sequential accumulator. It already encodes that an empty
        // dependency frontier starts at the typed root, that other prefixes enter through their
        // reduced branch, and that both local subtraction signs are applied without a raw
        // Foata-level product. Disconnected values likewise retain their factorized local and
        // integrated-prefix cross terms, so no separate active/frozen root-path replay is needed.
        let operation = &self.graph[node];
        let computed = self.compute_store.require(operation)?;
        if self.graph.is_disjoint_union(node) {
            return Ok(Full4dCts::from_factorized_local(
                computed.local_4d(operation)?,
            ));
        }
        Full4dCts::recursion_input(
            computed.local_4d(operation)?,
            computed.integrated(operation)?,
            self.source_spinney(node).renormalization_scheme,
            node == self.root,
            &self.source_spinney(node).lmb,
        )
    }

    fn cached_node_label_atom(&self, node: NodeIndex) -> Option<Atom> {
        self.cached_node_label_atoms
            .as_ref()
            .and_then(|labels| labels.get(&node))
            .cloned()
    }

    fn node_label_atom_factor(
        &self,
        frontier: NodeIndex,
        op: &HiddenData<SuBitGraph, EdgeIndex>,
    ) -> Atom {
        let approx = FunctionBuilder::new(GS.uv_approx);
        let frontier_atom = self.node_label_atom(frontier);

        let current = function!(op.order.symbol(), usize::from(op.data));
        let argument = if self.graph[frontier].covers().is_none() {
            current
        } else {
            let previous = self.graph[frontier]
                .covers()
                .expect("non-empty frontier cover must exist")
                .symbol();
            (current - previous) * frontier_atom
        };
        approx.add_arg(argument).finish()
    }

    fn node_label_atom(&self, node: NodeIndex) -> Atom {
        if self.graph[node].key.is_empty() {
            return Atom::one();
        }

        if let Some(cached) = self.cached_node_label_atom(node) {
            return cached;
        }

        self.graph
            .leaf_op_dependency_frontiers(node, &self.wood)
            .fold(Atom::one(), |acc, (op, frontier)| {
                acc * self.node_label_atom_factor(frontier, op)
            })
    }

    fn with_cached_node_label_atoms(mut self) -> Self {
        self.cache_node_label_atoms();
        self
    }

    fn cache_node_label_atoms(&mut self) {
        self.cached_node_label_atoms = Some(AHashMap::new());
        for nidx in self
            .graph
            .topo_sort_kahn(DirectionBasis::Underlying)
            .unwrap()
        {
            let atom = self.node_label_atom(nidx);
            self.cached_node_label_atoms
                .as_mut()
                .expect("node-label cache was initialized")
                .insert(nidx, atom);
        }
    }

    fn node_label(&self, node: NodeIndex) -> String {
        let key = &self.graph[node];
        if key.key.is_empty() {
            return "∅".to_string();
        }

        let mut foata = String::new();
        key.key
            .write_foata_like(&mut foata, |op| op.string_label())
            .expect("writing a trace key into a string must succeed");
        let atom = self.node_label_atom(node);
        format!("{foata}: {}", atom.to_ordered_simple())
    }

    fn dot_serialize_expr_atom() -> Atom {
        Atom::var(GS.expr)
    }

    fn dot_serialize_node_atom_factor(
        &self,
        frontier: NodeIndex,
        op: &HiddenData<SuBitGraph, EdgeIndex>,
    ) -> Atom {
        let frontier_depth = self.graph[frontier].key.op_count();
        let (current, given) = self.wood.current_given_pair(op.data, frontier_depth);
        // Structured forest-DOT approximation records use the same operation and subgraph
        // markers as the computed UV expressions.
        function!(
            GS.uv_approx,
            UvMarker::subgraph(current.subgraph(), given.subgraph())
                * self.dot_serialize_node_atom(frontier)
        )
    }

    fn dot_serialize_node_atom(&self, node: NodeIndex) -> Atom {
        if self.graph[node].key.is_empty() {
            return Self::dot_serialize_expr_atom();
        }

        self.graph
            .leaf_op_dependency_frontiers(node, &self.wood)
            .fold(Atom::one(), |acc, (op, frontier)| {
                acc * self.dot_serialize_node_atom_factor(frontier, op)
            })
    }

    fn dot_serialize_node_attrs(&self, node: NodeIndex) -> String {
        let key = &self.graph[node];
        let label = self
            .dot_serialize_node_atom(node)
            .printer(SpensoPrintSettings::typst_options())
            .to_string();
        let cover = key
            .covers()
            .unwrap_or_else(|| self.graph.empty_subgraph())
            .string_label();

        format!(
            "label={} foata={} cover={}",
            dot_attr_value(&label),
            dot_attr_value(&key.foata_level_labels()),
            dot_attr_value(&cover),
        )
    }

    pub fn dot_serialize(&self) -> String {
        let mut output = String::new();
        self.dot_serialize_fmt(&mut output)
            .expect("writing hedge-poset forest DOT into a string must succeed");
        output
    }

    pub fn dot_serialize_fmt(&self, writer: &mut impl std::fmt::Write) -> std::fmt::Result {
        let attrs: AHashMap<_, _> = self
            .graph
            .iter_nodes()
            .map(|(node, _, key)| (key.clone(), self.dot_serialize_node_attrs(node)))
            .collect();

        self.graph.dot_impl_fmt(
            writer,
            &self.graph.full_filter(),
            "start=2;\n",
            &|_| None,
            &|_| None,
            &|v| Some(attrs[v].clone()),
        )
    }

    fn compatible_topological_order(&self, subset: &SuBitGraph) -> Result<Vec<NodeIndex>> {
        let mut order = self
            .graph
            .topo_sort_kahn_of(subset, DirectionBasis::Underlying)?;

        if let Some(root_position) = order.iter().position(|node| *node == self.root) {
            if root_position != 0 {
                let root = order.remove(root_position);
                order.insert(0, root);
            }
        } else {
            order.insert(0, self.root);
        }

        Ok(order)
    }

    #[cfg(test)]
    fn local_leaf_operations(&self, node: NodeIndex) -> Vec<LocalLeafOperation> {
        let mut leaves = self
            .graph
            .leaf_op_dependency_frontiers(node, &self.wood)
            .map(|(op, frontier)| LocalLeafOperation::new(op, frontier))
            .collect::<Vec<_>>();

        leaves.sort_by_key(|leaf| {
            (
                Reverse(self.graph[leaf.frontier].key.op_count()),
                leaf.op.order.clone(),
                usize::from(leaf.op.data),
            )
        });
        leaves
    }

    fn disconnected_component_nodes(&self, node: NodeIndex) -> Result<Vec<NodeIndex>> {
        let operation = &self.graph[node];
        let factors = self
            .wood
            .join_factors(self.graph.source_node(node))
            .ok_or_else(|| eyre!("{operation} has no disconnected component factors"))?;
        let target_ops = operation
            .key
            .iter_levels_top_down()
            .flat_map(|level| level.iter_leaf_ops())
            .cloned()
            .collect::<Vec<_>>();

        for op in &target_ops {
            let owners = factors
                .iter()
                .filter(|factor| factor.includes(&op.order))
                .count();
            if owners != 1 {
                return Err(eyre!(
                    "operation {} in {operation} belongs to {owners} disconnected factors",
                    op.order.string_label()
                ));
            }
        }

        factors
            .into_iter()
            .map(|factor| {
                let key = target_ops
                    .iter()
                    .filter(|op| factor.includes(&op.order))
                    .fold(TraceKey::empty(), |key, op| {
                        key.push(&self.wood, op.clone())
                    });
                self.graph
                    .iter_nodes()
                    .find_map(|(component, _, candidate)| {
                        (candidate.key == key && self.source_spinney(component).filter() == &factor)
                            .then_some(component)
                    })
                    .ok_or_else(|| {
                        eyre!(
                            "component {} of {operation} is not in the unfolded forest",
                            factor.string_label()
                        )
                    })
            })
            .collect()
    }

    fn operation_node_index(&self, operation: &OperationNode) -> Result<NodeIndex> {
        let cover = operation.covers();
        self.graph
            .iter_nodes()
            .find_map(|(node, _, candidate)| {
                let source_matches = cover.as_ref().map_or_else(
                    || self.source_spinney(node).filter().is_empty(),
                    |cover| self.source_spinney(node).filter() == cover,
                );
                (candidate == operation && source_matches).then_some(node)
            })
            .ok_or_else(|| eyre!("{operation} is not in the unfolded forest"))
    }

    /// Orders the factors `components` of a factorized join for replay. A
    /// factor holding a connected zero-mode host whose soft set meets another
    /// factor replays first: the host's zero mode keeps the static lines of
    /// that factor, whose own U then acts on it (`U_B S0_A`), while a U of
    /// the factor before the host leaves it no zero mode (`S0_A U_B = 0`).
    /// Unconstrained factors keep their order.
    ///
    /// # Errors
    ///
    /// Factors that must each replay before another one.
    fn pinned_replay_order(
        &self,
        graph: &Graph,
        mut pending: Vec<NodeIndex>,
    ) -> Result<Vec<NodeIndex>, ZeroModeHostError> {
        let hosts = self
            .wood
            .graph
            .iter_nodes()
            .filter_map(|(_, _, spinney)| {
                spinney
                    .thermal_zero_mode()
                    .filter(|region| region.kind() == ZeroModeHostKind::Connected)
                    .map(|region| (spinney.filter(), region.sigma()))
            })
            .collect::<Vec<_>>();
        let precedes = |first: NodeIndex, second: NodeIndex| {
            let factor = self.source_spinney(first).filter();
            let other = graph.paired_edges(self.source_spinney(second).filter());
            hosts.iter().any(|(host, sigma)| {
                factor.includes(*host) && other.iter().any(|edge| sigma.contains(edge))
            })
        };

        let mut ordered = Vec::with_capacity(pending.len());
        while !pending.is_empty() {
            let Some(next) = pending.iter().position(|candidate| {
                !pending
                    .iter()
                    .any(|other| other != candidate && precedes(*other, *candidate))
            }) else {
                return Err(ZeroModeHostError::MutuallyCoveringJoin {
                    factors: pending
                        .iter()
                        .map(|component| {
                            graph
                                .paired_edges(self.source_spinney(*component).filter())
                                .into_iter()
                                .collect()
                        })
                        .collect(),
                });
            };
            ordered.push(pending.remove(next));
        }
        Ok(ordered)
    }

    /// The replay states of `node`: each pairs an integrated prefix with the
    /// local operations replayed on it, in order. The factors of a
    /// factorized join replay in their pinned order
    /// ([`Self::pinned_replay_order`]).
    fn union_replay_states(&self, graph: &Graph, node: NodeIndex) -> Result<Vec<UnionReplayState>> {
        let operation = &self.graph[node];
        if operation.key.is_empty() {
            return Ok(vec![UnionReplayState {
                integrated: node,
                local_edges: Vec::new(),
            }]);
        }

        if !self.graph.is_disjoint_union(node) {
            let (parent, edge) = self
                .graph
                .unique_parent(node)
                .ok_or_else(|| eyre!("{operation} has no unique replay parent"))?;
            let mut states = self.union_replay_states(graph, parent)?;
            for state in &mut states {
                state.local_edges.push(edge);
            }
            states.push(UnionReplayState {
                integrated: node,
                local_edges: Vec::new(),
            });
            return Ok(states);
        }

        let components = self.disconnected_component_nodes(node)?;
        let component_states = self
            .pinned_replay_order(graph, components)
            .wrap_err_with(|| format!("while ordering the factors of {operation} for replay"))?
            .into_iter()
            .map(|component| self.union_replay_states(graph, component))
            .collect::<Result<Vec<_>>>()?;
        let mut states = Vec::new();
        for components in component_states
            .iter()
            .map(|states| states.iter())
            .multi_cartesian_product()
        {
            let key = TraceKey::try_foata_join(
                components
                    .iter()
                    .map(|component| &self.graph[component.integrated].key),
                &self.wood,
            )
            .ok_or_else(|| eyre!("cannot join integrated prefixes for {operation}"))?;
            let mut local_edges = Vec::new();
            for component in components {
                local_edges.extend_from_slice(&component.local_edges);
            }
            let integrated = self.operation_node_index(&OperationNode { key })?;
            states.push(UnionReplayState {
                integrated,
                local_edges,
            });
        }

        Ok(states)
    }

    /// The plain U ([`Spinney::to_plain_u`]) that a replay applies in place
    /// of the connected zero-mode host `current` to a sector whose integrated
    /// prefix `prefix` meets the host's soft set, if it does.
    ///
    /// The factor of a join that holds a connected host replays before the
    /// factors in the host's soft set ([`Self::pinned_replay_order`]), so the
    /// U of such a factor acts on the host's zero mode on static rows, as a
    /// static U without an integrated counterpart. The integrated prefix of
    /// that factor therefore meets the host's U alone, also where it
    /// contracted only part of a frequency class of the host, which would
    /// leave the zero-mode projection incomplete.
    ///
    /// # Errors
    ///
    /// A union host whose soft set meets `prefix`, an internal error
    /// ([`ZeroModeHostError::PrefixInUnionSoftSet`]): a factorized join does
    /// not produce it, and the host has no U to keep.
    fn replayed_plain_u(
        &self,
        graph: &Graph,
        current: &ForestNode<'_>,
        prefix: &SuBitGraph,
    ) -> Result<Option<Spinney>> {
        let Some(region) = current.thermal_zero_mode() else {
            return Ok(None);
        };
        let prefix = graph.paired_edges(prefix);
        if !prefix.iter().any(|edge| region.sigma().contains(edge)) {
            return Ok(None);
        }
        match region.kind() {
            ZeroModeHostKind::Connected => Ok(Some(current.spinney.to_plain_u())),
            ZeroModeHostKind::Union => Err(ZeroModeHostError::PrefixInUnionSoftSet {
                host: graph.paired_edges(current.subgraph()).into_iter().collect(),
                prefix: prefix.into_iter().collect(),
            }
            .into()),
        }
    }

    /// Whether `node` has a counterpart in the zero-temperature integrated
    /// lane ([`Self::compute_4d_for_node`]). Only UV counterterms with the
    /// vacuum mass have integrated counterparts. A zero-mode operation and
    /// the static U's above it are local only, so a union host and every
    /// node above one have none; a connected host enters by its U. Under
    /// `vacuum_subtraction` an operation covering every loop has none
    /// either: its add-back would sit on a loopless cograph, whose
    /// vacuum-subtracted weight vanishes, as its local U does.
    fn has_zero_temperature_counterpart(
        &self,
        node: NodeIndex,
        graph: &Graph,
        vacuum_subtraction: bool,
    ) -> bool {
        let operation = &self.graph[node];
        let above_union_host = operation
            .key
            .iter_levels_top_down()
            .flat_map(|level| level.iter_leaf_ops())
            .any(|op| {
                let HedgePair::Paired { sink, .. } = self.wood.graph[&op.data].1 else {
                    panic!("edge in trace key is not paired");
                };
                self.wood.graph[self.wood.graph.node_id(sink)].wood_role() == WoodRole::UnionHost
            });
        let covers_every_loop = operation
            .covers()
            .is_some_and(|cover| graph.n_loops(&cover) == graph.n_loops(&graph.full_filter()));
        !above_union_host && !(vacuum_subtraction && covers_every_loop)
    }

    /// The zero-temperature local counterterm of `node` and its integrated
    /// value, from those of the nodes below it. A node without a counterpart
    /// there ([`Self::has_zero_temperature_counterpart`]) gets zero, without
    /// an integration. A connected zero-mode host is its plain U
    /// ([`Spinney::to_plain_u`]): the replay of its prefix dresses that
    /// add-back with the host's own zero mode
    /// ([`Direct3dApproximation::localize_prefix`]).
    fn compute_4d_for_node(
        &self,
        node: NodeIndex,
        graph: &Graph,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
        vacuum_subtraction: bool,
    ) -> Result<(Local4dCts, IntegratedCts)> {
        let operation = &self.graph[node];
        if operation.key.is_empty() {
            return Ok((Local4dCts::root(), IntegratedCts::root()));
        }
        if !self.has_zero_temperature_counterpart(node, graph, vacuum_subtraction) {
            return Ok((Local4dCts::zero(), IntegratedCts::root()));
        }
        let source = self.source_spinney(node);
        let cover = operation
            .covers()
            .expect("a non-root operation has a cover");

        if self.graph.is_disjoint_union(node) {
            let components = self.disconnected_component_nodes(node)?;
            let mut full_components = Vec::with_capacity(components.len());
            let mut integrated_components = Vec::with_capacity(components.len());
            for component in components {
                let component_operation = &self.graph[component];
                full_components.push(self.recursion_input_4d(component).wrap_err_with(|| {
                    format!("while loading 4D component {component_operation} for {operation}")
                })?);
                integrated_components.push(
                    self.compute_store
                        .require(component_operation)?
                        .integrated(component_operation)?,
                );
            }

            let depth = graph.n_loops(&cover) + 1;
            return Ok((
                Local4dCts::from_full_product(full_components),
                IntegratedCts::factorized_product(integrated_components, depth)?,
            ));
        }

        let (parent, edge) = self
            .graph
            .unique_parent(node)
            .ok_or_else(|| eyre!("{operation} has no unique 4D parent"))?;
        let parent_operation = &self.graph[parent];
        let full = self
            .recursion_input_4d(parent)
            .wrap_err_with(|| format!("while loading the 4D parent for {operation}"))?;
        let step_order = parent_operation.key.op_count();
        let (current, given) = self.wood.current_given_pair(edge, step_order);
        let plain_u;
        let current = if source.thermal_zero_mode().is_some() {
            plain_u = source.to_plain_u();
            ForestNode {
                spinney: &plain_u,
                topo_order: current.topo_order,
            }
        } else {
            current
        };
        let ctx = UVCtx::new(graph, settings);
        let integrated_approximation = Integrated::new(vakint, &self.wood.vakint_settings);
        let local = local_4d::uv_limit(&full, &ctx, &current, &given, &current, &given)?;
        let integrated = if settings.generate_integrated {
            integrated_approximation.run(&local, &ctx, &current, &given, &current, &given)?
        } else {
            IntegratedCts::root()
        };

        Ok((local, integrated))
    }

    fn local_3d_for_node(
        &self,
        node: NodeIndex,
        graph: &mut Graph,
        cutset: &CutSet,
        localizer: Localizer<'_>,
        settings: &UVgenerationSettings,
        projection_context: &mut Local4dProjectionContext,
    ) -> Result<CutComputation> {
        let operation = &self.graph[node];
        let forest_node = ForestNode {
            spinney: self.source_spinney(node),
            topo_order: operation.key.op_count(),
        };
        let local_3d = if operation.key.is_empty() {
            Local3DCts::Direct(Direct3dCts::root(graph, localizer)?)
        } else if settings.local_uv_cts_from_expanded_4d_integrands {
            // This is the only route that may produce `Projected4d`: the 4D
            // Taylor coefficient is kept factorized and final assembly later
            // attaches its untouched outer CFF. It never enters any of the
            // `Direct3dApproximation` replay operations below.
            let local_4d = self.compute_store.require(operation)?.local_4d(operation)?;
            Local3DCts::Projected4d(
                Projected4dApproximation::new(localizer, graph, settings)
                    .project_local_4d(local_4d, projection_context)?,
            )
        } else if self.graph.is_disjoint_union(node) {
            // Both direct variants replay the Taylor operators on the complete
            // post-energy-integration CFF. `explicit_orientation_sum_only`
            // changes selector materialization only, not this construction.
            let mut active_sectors = Vec::new();
            for state in self
                .union_replay_states(graph, node)?
                .into_iter()
                .filter(|state| !state.local_edges.is_empty())
            {
                let integrated_operation = &self.graph[state.integrated];
                let prefix_cover = integrated_operation
                    .covers()
                    .unwrap_or_else(|| graph.empty_subgraph());
                // An empty integrated prefix starts from the per-cut root integrand. Every
                // other prefix enters through the reduced branch of its first local operation.
                let mut sector = if integrated_operation.key.is_empty() {
                    Some(Direct3dCts::root(graph, localizer)?)
                } else {
                    None
                };
                for (offset, edge) in state.local_edges.into_iter().enumerate() {
                    let step_order = integrated_operation.key.op_count() + offset;
                    let (current, given) = self.wood.current_given_pair(edge, step_order);
                    let plain_u = self.replayed_plain_u(graph, &current, &prefix_cover)?;
                    let current = match &plain_u {
                        Some(spinney) => ForestNode {
                            spinney,
                            topo_order: step_order,
                        },
                        None => current,
                    };
                    sector = Some(match sector {
                        Some(sector) => Direct3dApproximation::new(localizer, graph, settings)
                            .run_local(&sector, &current, &given, &current, &given)?,
                        None => {
                            // Retain each connected coefficient's owner before the
                            // product enters a local suffix. Equal physical masses in
                            // a disjoint sibling must not acquire the parent's grading.
                            // Keep disconnected vacuum masses owned as well; the
                            // complete union still determines the localization kernel.
                            let components = if self.graph.is_disjoint_union(state.integrated) {
                                self.disconnected_component_nodes(state.integrated)?
                            } else {
                                vec![state.integrated]
                            };
                            let finite_components = components
                                .into_iter()
                                .map(|component| {
                                    let operation = &self.graph[component];
                                    let integrated = self
                                        .compute_store
                                        .require(operation)?
                                        .integrated(operation)?;
                                    Ok((
                                        integrated.physical_finite_counterterm_atom(),
                                        ForestNode {
                                            spinney: self.source_spinney(component),
                                            topo_order: operation.key.op_count(),
                                        },
                                    ))
                                })
                                .collect::<Result<Vec<_>>>()?;
                            let finite_counterterms = finite_components
                                .iter()
                                .map(|(coefficient, owner)| (coefficient, owner))
                                .collect::<Vec<_>>();
                            // A disabled (or algebraically vanishing) integrated prefix
                            // still owns its typed component-path frame. Replay that
                            // zero coefficient through the local suffix; it remains
                            // algebraically zero without deleting the union sector.
                            let prefix_node = ForestNode {
                                spinney: self.source_spinney(state.integrated),
                                topo_order: integrated_operation.key.op_count(),
                            };
                            Direct3dApproximation::new(localizer, graph, settings).run_integrated(
                                &finite_counterterms,
                                &prefix_node,
                                &current,
                                &given,
                                &current,
                                &given,
                            )?
                        }
                    });
                }
                let sector = sector.expect("a proper integrated prefix has a local suffix");

                // Keep each active/frozen split after its root-path replay so any connected
                // descendants rescale only the loop variables still active in that sector.
                active_sectors.extend(sector.sectors()?.iter().cloned());
            }

            Local3DCts::Direct(
                Direct3dCts::from_sectors(active_sectors)
                    .wrap_err_with(|| format!("{operation} has no proper integrated prefixes"))?,
            )
        } else {
            let (parent, edge) = self
                .graph
                .unique_parent(node)
                .ok_or_else(|| eyre!("{operation} has no unique local parent"))?;
            let parent_operation = &self.graph[parent];
            // An empty dependency frontier starts from the per-cut root integrand;
            // otherwise its typed local result remains the sequential accumulator.
            let root;
            let parent_local = if parent_operation.key.is_empty() {
                root = Direct3dCts::root(graph, localizer)?;
                &root
            } else {
                self.compute_store
                    .require(parent_operation)?
                    .cut(parent_operation, cutset)?
                    .local_3d
                    .direct()?
            };
            let parent_integrated = self
                .compute_store
                .require(parent_operation)?
                .integrated(parent_operation)
                .wrap_err_with(|| {
                    format!("while loading integrated parent {parent_operation} for {operation}")
                })?;
            let step_order = parent_operation.key.op_count();
            let (current, given) = self.wood.current_given_pair(edge, step_order);
            // `run` applies the next Taylor operator to the complete CFF and
            // applies both subtraction signs; no projected 4D coefficient,
            // external sign, or raw Foata-level product enters this route.
            Local3DCts::Direct(Direct3dApproximation::new(localizer, graph, settings).run(
                parent_local,
                parent_integrated,
                &current,
                &given,
                &current,
                &given,
            )?)
        };

        let integrated = self
            .compute_store
            .require(operation)?
            .integrated(operation)?;
        let final_builder = FinalIntegrandBuilder::new(localizer, settings);
        let final_integrands = match &local_3d {
            Local3DCts::Direct(direct) => {
                final_builder.build_direct(graph, &forest_node, direct, integrated)?
            }
            Local3DCts::Projected4d(projected) => {
                final_builder.build_projected(graph, &forest_node, projected, integrated)?
            }
        };

        Ok(CutComputation {
            local_3d,
            final_integrands,
        })
    }

    /// Compute the zero-temperature local counterterm and the integrated
    /// value of every node ([`Self::compute_4d_for_node`]).
    /// `vacuum_subtraction` is that of the local 3D lane localizing these
    /// values, if any.
    pub fn integrate(
        &mut self,
        graph: &Graph,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
        vacuum_subtraction: bool,
    ) -> Result<()> {
        for (order, nidx) in self
            .graph
            .topo_sort_kahn(DirectionBasis::Underlying)?
            .into_iter()
            .enumerate()
        {
            debug!(order, nidx=%nidx, key=%self.graph[nidx], "Computing hedge-poset 4D term");
            let operation = self.graph[nidx].clone();
            let (local_4d, integrated) =
                self.compute_4d_for_node(nidx, graph, vakint, settings, vacuum_subtraction)?;
            let cover = operation
                .covers()
                .unwrap_or_else(|| self.graph.empty_subgraph())
                .string_label();
            let source = self.source_spinney(nidx).filter().string_label();
            debug_tags!(#generation, #uv, #integrated, #graph, #term, #inspect;
                stage = "hedge_poset_4d_node_done",
                order,
                node_index = %nidx,
                forest_term = %operation,
                cover = %cover,
                source = %source,
                is_union = self.graph.is_disjoint_union(nidx),
                log.local_4d = local_4d.atom(),
                log.integrated_pole = integrated.physical_pole_atom(),
                log.integrated_finite = integrated.physical_finite_counterterm_atom(),
                "Computed hedge-poset 4D node"
            );
            let computed = self.compute_store.entry(operation).or_default();
            computed.local_4d = Some(local_4d);
            computed.integrated = Some(integrated);
        }

        Ok(())
    }

    /// Return the measured representation-independent part of this computation.
    pub(crate) fn compute(
        &mut self,
        graph: &mut Graph,
        vakint: &Vakint,
        orientation: OrientationProjection<'_>,
        settings: &UVgenerationSettings,
        compute_four_d: bool,
        projection_context: &mut Local4dProjectionContext,
    ) -> Result<std::time::Duration> {
        let shared_before = projection_context.shared_preparation_time;
        if compute_four_d {
            let started = std::time::Instant::now();
            if settings.generate_integrated || settings.local_uv_cts_from_expanded_4d_integrands {
                // Local-4D projection also needs these coefficients when analytic addbacks are off.
                self.integrate(graph, vakint, settings, orientation.vacuum_subtraction)?;
            } else {
                // Direct local-3D consumes only the zero integrated lane. Seed every
                // operation because a disconnected replay can visit a prefix first.
                let zero = IntegratedCts::root();
                for (_, _, operation) in self.graph.iter_nodes() {
                    self.compute_store
                        .entry(operation.clone())
                        .or_default()
                        .integrated = Some(zero.clone());
                }
            }
            projection_context.shared_preparation_time += started.elapsed();
        }

        for (compatible_subset, cutset) in self.cuts.clone() {
            let localizer = Localizer::new(&cutset, orientation);
            for (order, nidx) in self
                .compatible_topological_order(&compatible_subset)?
                .into_iter()
                .enumerate()
            {
                debug!(order, nidx=%nidx, key=%self.graph[nidx], "Computing hedge-poset per-cut term");
                let operation = self.graph[nidx].clone();
                // Direct local-3D nodes Taylor-expand the complete post-energy-integration CFF.
                // Expanded-4D nodes instead project their typed local coefficients and attach
                // the outer CFF only during final assembly.
                let cut_computation = self.local_3d_for_node(
                    nidx,
                    graph,
                    &cutset,
                    localizer,
                    settings,
                    projection_context,
                )?;
                self.compute_store
                    .entry(operation)
                    .or_default()
                    .cuts
                    .insert(cutset.clone(), cut_computation);
            }
        }

        Ok(projection_context.shared_preparation_time - shared_before)
    }

    pub(crate) fn orientation_parametric_exprs(
        &self,
        graph: &Graph,
        _settings: &UVgenerationSettings,
    ) -> Result<Vec<ParametricIntegrands>> {
        let split_momentum_replacements = graph
            .iter_edges_of(
                &graph
                    .full_filter()
                    .subtract(&graph.initial_state_cut)
                    .subtract(&graph.tree_edges),
            )
            .filter_map(|(pair, edge_index, _)| {
                (!matches!(pair, HedgePair::Unpaired { .. }))
                    .then(|| GS.split_mom_pattern_simple(edge_index))
            })
            .collect::<Vec<_>>();
        let mut expressions = Vec::with_capacity(self.cuts.len());

        for (compatible_subset, cutset) in &self.cuts {
            let mut sum: Option<FinalIntegrands> = None;
            for nidx in self.compatible_topological_order(compatible_subset)? {
                let operation = &self.graph[nidx];
                let terms = self
                    .compute_store
                    .require(operation)?
                    .cut(operation, cutset)?
                    .final_integrands
                    .map_expressions(|integrand| Ok(integrand.clone().collect_color()))?;
                sum = Some(match sum {
                    Some(sum) => sum.zip_add([terms]).wrap_err_with(|| {
                        format!("while aggregating hedge-poset term {operation} for cut {cutset:?}")
                    })?,
                    None => terms,
                });
            }
            // A per-node final integrand retains inspect/export provenance but
            // is not independently evaluator-ready.  In particular, temporal
            // components in localized integrated counterterms can cancel only
            // in the completed signed forest.  Preserve that exact cancellation
            // before the final orientation-local momentum split.
            let integrands = sum
                .ok_or_else(|| eyre!("No terms in hedge-poset forest for cut {cutset:?}"))?
                .map_expressions(|integrand| {
                    Ok(integrand
                        .replace_multiple(&split_momentum_replacements)
                        .replace(function!(GS.den, W_.a_, W_.b_, W_.c_, W_.d_))
                        .with(W_.d_))
                })?;
            expressions.push(ParametricIntegrands::from_final(integrands, cutset.clone()));
        }

        Ok(expressions)
    }

    pub(crate) fn export_node_expressions(
        &self,
        forest_index: usize,
        post_process: &mut impl FnMut(Atom) -> Atom,
    ) -> Result<Vec<UVForestNodeExpression>> {
        let (cut_compatible_forest_subset, cutset) = self
            .cuts
            .first()
            .ok_or_else(|| eyre!("No cuts in hedge-poset forest export"))?;
        let mut terms = Vec::new();
        for (node_index, nidx) in self
            .compatible_topological_order(cut_compatible_forest_subset)?
            .into_iter()
            .enumerate()
        {
            let operation = &self.graph[nidx];
            let computed = self.compute_store.require(operation)?;
            let cut_computation = computed.cut(operation, cutset)?;
            let final_integrands = cut_computation
                .final_integrands
                .clone()
                .into_integrands()
                .resolved()?
                .map(|numerator| post_process(numerator.clone()));
            let node_key = operation.to_string();
            let provenance = self.node_export_3d_provenance(nidx, &cut_computation.local_3d);
            for (term_index, (&residue_index, numerator)) in final_integrands.iter().enumerate() {
                terms.push(UVForestNodeExpression {
                    forest_index,
                    node_index,
                    node_key: node_key.clone(),
                    term_index,
                    residue_index,
                    numerator: numerator.clone(),
                    forest_provenance: Some(provenance.clone()),
                });
            }
        }

        Ok(terms)
    }

    pub(crate) fn export_4d_node_expressions(
        &self,
        forest_index: usize,
    ) -> Result<Vec<UVForestNodeExpression>> {
        let (cut_compatible_forest_subset, _) = self
            .cuts
            .first()
            .ok_or_else(|| eyre!("No cuts in hedge-poset forest export"))?;
        self.compatible_topological_order(cut_compatible_forest_subset)?
            .into_iter()
            .enumerate()
            .map(|(node_index, nidx)| {
                let operation = &self.graph[nidx];
                let local_4d = self.compute_store.require(operation)?.local_4d(operation)?;
                Ok(UVForestNodeExpression {
                    forest_index,
                    node_index,
                    node_key: operation.to_string(),
                    term_index: 0,
                    residue_index: crate::cff::CutCFFIndex::new_all_none(),
                    numerator: local_4d.atom().clone(),
                    forest_provenance: Some(self.node_export_4d_provenance(nidx, local_4d)),
                })
            })
            .collect()
    }

    fn node_parent_keys(&self, node: NodeIndex) -> Vec<String> {
        let mut parent_keys = self
            .graph
            .parents(node)
            .map(|(parent, _)| self.graph[parent].to_string())
            .collect::<Vec<_>>();
        parent_keys.sort();
        parent_keys.dedup();
        parent_keys
    }

    fn node_export_4d_provenance(
        &self,
        node: NodeIndex,
        local_4d: &Local4dCts,
    ) -> Arc<UVForestProvenance> {
        Arc::new(UVForestProvenance {
            node: UVForestNodeProvenance {
                parent_keys: self.node_parent_keys(node),
                local: UVForestLocalProvenance::FourD {
                    components: local_4d.provenance_components().to_vec(),
                    branches: local_4d.provenance_branches().to_vec(),
                },
            },
        })
    }

    fn node_export_3d_provenance(
        &self,
        node: NodeIndex,
        local_3d: &Local3DCts,
    ) -> Arc<UVForestProvenance> {
        Arc::new(UVForestProvenance {
            node: UVForestNodeProvenance {
                parent_keys: self.node_parent_keys(node),
                local: UVForestLocalProvenance::ThreeD {
                    projection_paths: local_3d.projection_paths(),
                },
            },
        })
    }

    pub(crate) fn renormalization_part_of_ends(
        &self,
        graph: &Graph,
        settings: &UVgenerationSettings,
    ) -> Result<RenormalizationPart> {
        let mut sum = Atom::Zero;
        let marker = UvMarker::new(settings);

        let wild = Atom::var(W_.x___);

        let replacements =
            graph.integrand_replacement(&graph.full_filter(), &graph.loop_momentum_basis, &[wild]);
        for (node, mut crown, key) in self.graph.iter_nodes() {
            if crown.any(|r| self.graph.flow(r).is_source()) {
                continue;
            }

            let forest_node =
                key.forest_node(graph, self.source_spinney(node), key.key.op_count())?;
            // A disconnected terminal can combine components with different schemes. Select
            // each component's projection before multiplying; the aggregate integrated value
            // only represents the homogeneous all-finite and all-pole projections.
            let components = if self.graph.is_disjoint_union(node) {
                self.disconnected_component_nodes(node)?
            } else {
                vec![node]
            };
            let physical = components.into_iter().try_fold(
                Atom::one(),
                |product, component| -> Result<Atom> {
                    let component_key = &self.graph[component];
                    let integrated = self
                        .compute_store
                        .require(component_key)?
                        .integrated(component_key)?;
                    let projection = match self.source_spinney(component).renormalization_scheme {
                        ApproximationType::MUV | ApproximationType::IR => {
                            integrated.physical_finite_counterterm_atom()
                        }
                        ApproximationType::PolePart => integrated.physical_pole_atom(),
                        scheme => {
                            return Err(eyre!("No terminal counterterm projection for {scheme}"));
                        }
                    };
                    Ok(product * projection)
                },
            )?;
            let atom = marker.prefix(&graph.full_filter(), forest_node.subgraph(), &physical);
            debug!(
                key=%key,
               expr = % atom.expand_num().log_print(None),"Term before simplification"
            );
            let atom = (&atom
                * &graph.global_prefactor.projector
                * &graph.global_prefactor.num
                * &graph.overall_factor)
                .simplify_color()
                .expand_num()
                .to_dots();

            debug!(
                key=%key,
               expr = % atom.log_print(None),"Term"
            );
            sum += atom;
        }

        Ok(RenormalizationPart::new(
            sum.replace_multiple(&replacements)
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum),
            self.graph.n_nodes(),
        ))
    }

    pub fn debug_walk(&self) {
        let mut cover_groups: BTreeMap<SuBitGraph, Vec<NodeIndex>> = BTreeMap::new();

        self.graph
            .topo_sort_kahn(DirectionBasis::Underlying)
            .unwrap()
            .iter()
            .for_each(|nidx| {
                let trace_key = &self.graph[*nidx];
                cover_groups
                    .entry(trace_key.covers().unwrap_or(self.graph.empty_subgraph()))
                    .and_modify(|e| e.push(*nidx))
                    .or_insert_with(|| vec![*nidx]);

                println!("Node {}:{}", nidx, self.node_label(*nidx));
            });

        println!("edge [constraint=true style=invis];");
        for (a, b) in cover_groups.values().tuple_windows() {
            println!("{}->{}", a.first().unwrap(), b.first().unwrap())
        }
        println!("edge [style=solid];");

        for (s, g) in cover_groups.iter() {
            println!("subgraph group_{} {{rank=same; ", s.string_label());
            for n in g {
                println!("{};", n.0);
            }
            println!("}}");
        }
    }
}

impl Display for Forests {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let labels: AHashMap<_, _> = self
            .graph
            .iter_nodes()
            .map(|(node, _, key)| (key.clone(), self.node_label(node)))
            .collect();
        self.graph.dot_impl_fmt(
            f,
            &self.graph.full_filter(),
            "start=2;\n",
            &|_| None,
            &|_| None,
            &|v| Some(format!("label=\"{}\"", labels[v])),
        )
    }
}

#[cfg(test)]
mod tests {
    use crate::initialisation;
    use crate::numerator::aind::Aind;
    use crate::utils;
    use crate::uv::uv_graph::UVE;
    use crate::{
        dot,
        graph::{FeynmanGraph, Graph, parse::IntoGraph},
        initialisation::test_initialise,
        processes::DotExportSettings,
        uv::{
            UltravioletGraph, Wood as OldWood, approx::local_4d::Local4dBranch,
            settings::RenormalizationPrescriptionSettings,
        },
    };
    use idenso::IndexTooling;
    use idenso::color::ColorSimplifySettings;
    use idenso::shorthands::metric::MetricSimplifier;

    use super::*;
    use crate::uv::approx::local_3d::{Local3DConceptualBranch, Local3DLoopRescaling};
    use crate::uv::test_fixtures::{Ring, edge_set, ring_uv_settings, thermal_medium};
    use crate::uv::thermal_zero_mode::tests::{connected_region, zero_mode_policy};
    use color_eyre::Result;
    use std::fmt::{Formatter, Result as FmtResult};
    use symbolica::atom::Symbol;
    use symbolica::id::Replacement;
    use symbolica::symbol;
    use three_dimensional_reps::RepresentationMode;

    /// The unfolded wood of a thermal ring under the zero-mode host policy:
    /// the IR prescription on the direct local-3D route at finite temperature,
    /// with vacuum subtraction.
    fn thermal_ring_forests(graph: &Graph) -> Result<Forests> {
        let settings = GenerationSettings {
            uv: ring_uv_settings(ApproximationType::IR),
            medium: thermal_medium(true),
            ..Default::default()
        };
        Ok(Wood::new(CutStructure::empty(graph), graph, &settings)?.unfold())
    }

    /// The Foata levels of the trace at `node`, in order of application, each
    /// as the set of its operation scopes.
    fn trace_levels(forests: &Forests, node: NodeIndex) -> Vec<BTreeSet<SuBitGraph>> {
        forests.graph[node]
            .key
            .iter_levels_top_down()
            .map(|level| level.iter_leaf_ops().map(|op| op.order.clone()).collect())
            .collect()
    }

    /// The factorized join of the plain union `union` of a ring's bubbles.
    fn ring_join(forests: &Forests, union: &SuBitGraph) -> NodeIndex {
        forests
            .graph
            .iter_nodes()
            .find_map(|(node, _, _)| {
                (forests.graph.is_disjoint_union(node)
                    && forests.source_spinney(node).filter() == union)
                    .then_some(node)
            })
            .expect("the ring has the factorized join of its two bubbles")
    }

    impl Forests {
        /// Seeds the integrated counterterm of every node with the symbolic
        /// integral `(c<node> M^dod + e<node> mUVexp^2 M^(dod - 2)) (1 + ε c'<node>)`
        /// ([`seed_symbol`], [`expansion_seed_symbol`]) in place of its analytic
        /// integration, as the zero-temperature lane decides
        /// ([`Self::compute_4d_for_node`]): zero without a counterpart there,
        /// the product of the components' for a factorized join. Like vakint's
        /// counterterms, the vertex holds the vacuum mass both as `mUV` and as
        /// the expansion mass `mUVexp` that the integrated U leaves in its
        /// numerators.
        fn seed_symbolic_integrated(
            &mut self,
            graph: &Graph,
            vacuum_subtraction: bool,
        ) -> Result<()> {
            for node in self.graph.topo_sort_kahn(DirectionBasis::Underlying)? {
                let operation = self.graph[node].clone();
                let integrated = if operation.key.is_empty()
                    || !self.has_zero_temperature_counterpart(node, graph, vacuum_subtraction)
                {
                    IntegratedCts::root()
                } else {
                    let depth = graph.n_loops(&operation.covers().expect("a proper cover")) + 1;
                    if self.graph.is_disjoint_union(node) {
                        let factors = self
                            .disconnected_component_nodes(node)?
                            .into_iter()
                            .map(|component| {
                                let operation = &self.graph[component];
                                self.compute_store.require(operation)?.integrated(operation)
                            })
                            .collect::<Result<Vec<_>>>()?;
                        IntegratedCts::factorized_product(factors, depth)?
                    } else {
                        // The integrated vertex has the mass dimension of its
                        // subgraph, which the soft charts read off.
                        let dod = i64::from(self.source_spinney(node).dod);
                        let mass = Atom::var(GS.m_uv_vacuum);
                        let finite = Atom::var(seed_symbol(node)) * mass.pow(dod)
                            + Atom::var(expansion_seed_symbol(node))
                                * Atom::var(GS.m_uv_expansion).pow(2)
                                * mass.pow(dod - 2);
                        let positive =
                            Atom::var(symbol!(&format!("hedge_poset_test_eps_c{}", node.0)));
                        IntegratedCts::from_integral(
                            &(&finite * (Atom::one() + Atom::var(GS.dim_epsilon) * positive)),
                            depth,
                        )?
                    }
                };
                self.compute_store.entry(operation).or_default().integrated = Some(integrated);
            }
            Ok(())
        }

        /// Names and subgraphs for the string labels of the wood's subgraphs
        /// and their differences: `A` = {e2, e3}, `B` the second loop, `AB`
        /// their union, `G` the whole graph, `D<edges>` a two-loop spinney and
        /// `∅` the empty subgraph; a difference is named after the spinney it
        /// equals, else as `X-Y`.
        fn subgraph_names(&self, graph: &Graph) -> BTreeMap<String, (String, SuBitGraph)> {
            let filters = self
                .wood
                .graph
                .iter_nodes()
                .map(|(_, _, spinney)| spinney.filter().clone())
                .collect::<BTreeSet<_>>();
            let edges = |subgraph: &SuBitGraph| {
                graph
                    .paired_edges(subgraph)
                    .into_iter()
                    .map(usize::from)
                    .collect::<Vec<_>>()
            };
            let all = (0..graph.underlying.n_edges()).collect::<Vec<_>>();
            let name = |subgraph: &SuBitGraph| match edges(subgraph).as_slice() {
                [] => "∅".to_string(),
                [2, 3] => "A".to_string(),
                [4, 5] | [4] => "B".to_string(),
                [2, 3, 4, 5] | [2, 3, 4] => "AB".to_string(),
                edges if edges == all.as_slice() => "G".to_string(),
                edges => format!("D{}", edges.iter().join("")),
            };
            let mut names = BTreeMap::new();
            for outer in &filters {
                for inner in &filters {
                    let difference = outer.subtract(inner);
                    let named = if filters.contains(&difference) || difference.is_empty() {
                        name(&difference)
                    } else if outer.includes(inner) {
                        format!("{}-{}", name(outer), name(inner))
                    } else {
                        continue;
                    };
                    names
                        .entry(difference.string_label())
                        .or_insert((named, difference));
                }
            }
            names
        }

        /// The rescaled support of an operation on `component` as the
        /// component less its integrated part: `X` for the whole of `X`,
        /// `X-Y` once the part `Y` is integrated, `∅` once all of it is.
        fn support_name(
            names: &BTreeMap<String, (String, SuBitGraph)>,
            component: &str,
            support: &str,
        ) -> String {
            let (component_name, component) = &names[component];
            let (_, support) = &names[support];
            if support == component {
                component_name.clone()
            } else if support.is_empty() {
                "∅".to_string()
            } else {
                format!(
                    "{component_name}-{}",
                    names[&component.subtract(support).string_label()].0
                )
            }
        }

        /// The trace of `node` as its Foata levels, each the sorted names
        /// of its operation scopes.
        fn trace_name(
            &self,
            node: NodeIndex,
            names: &BTreeMap<String, (String, SuBitGraph)>,
        ) -> String {
            trace_levels(self, node)
                .iter()
                .map(|level| {
                    level
                        .iter()
                        .map(|scope| names[&scope.string_label()].0.clone())
                        .sorted()
                        .join(",")
                })
                .join(";")
        }

        /// The shapes of the sectors of the local 3D counterterm of `node`
        /// ([`SectorShape`]), vanishing ones included, sorted.
        fn sector_shapes(
            &self,
            node: NodeIndex,
            names: &BTreeMap<String, (String, SuBitGraph)>,
        ) -> Result<Vec<SectorShape>> {
            let operation = &self.graph[node];
            if operation.key.is_empty() {
                return Ok(Vec::new());
            }
            let cutset = self.cuts.first().expect("a cut").1.clone();
            let local = self
                .compute_store
                .require(operation)?
                .cut(operation, &cutset)?
                .local_3d
                .direct()?;
            let seeded = self
                .graph
                .iter_nodes()
                .map(|(node, _, _)| {
                    (
                        seed_symbol(node),
                        expansion_seed_symbol(node),
                        self.trace_name(node, names),
                    )
                })
                .collect::<Vec<_>>();
            // The expansion masses of the seeded vertices, which no local U
            // introduced.
            let without_seeded_expansion_masses = seeded
                .iter()
                .map(|(_, expansion, _)| {
                    Replacement::new(Atom::var(*expansion).to_pattern(), Atom::Zero.to_pattern())
                })
                .collect::<Vec<_>>();
            let mut shapes = Vec::new();
            for sector in local.sectors()? {
                let atom = sector.active.factorized_sum();
                let [path] = sector.projection_paths.as_slice() else {
                    panic!("a ring sector has one projection path");
                };
                let steps = path
                    .steps
                    .iter()
                    .map(|step| {
                        let branches = step
                            .conceptual_branches
                            .iter()
                            .map(|branch| branch.branch)
                            .collect::<Vec<_>>();
                        let kind = match branches.as_slice() {
                            [Local3DConceptualBranch::U] => StepKind::U,
                            [Local3DConceptualBranch::S] => StepKind::S0,
                            _ => StepKind::H,
                        };
                        (
                            kind,
                            names[&step.current_component].0.clone(),
                            Self::support_name(
                                names,
                                &step.current_component,
                                &step.rescaled_subgraph,
                            ),
                        )
                    })
                    .collect();
                let holds = |terminal: Symbol, expansion: Symbol| {
                    atom.contains_symbol(terminal) || atom.contains_symbol(expansion)
                };
                shapes.push(SectorShape {
                    local: path
                        .steps
                        .first()
                        .is_some_and(|step| step.rescaling == Local3DLoopRescaling::FullSubgraph),
                    vanishes: atom.is_zero(),
                    prefix: seeded
                        .iter()
                        .filter(|(symbol, _, _)| atom.contains_symbol(*symbol))
                        .map(|(_, _, trace)| trace.clone())
                        .collect(),
                    keeps_expansion_masses: seeded
                        .iter()
                        .filter(|(symbol, _, _)| atom.contains_symbol(*symbol))
                        .all(|(_, expansion, _)| atom.contains_symbol(*expansion)),
                    steps,
                    static_mass: holds(GS.m_ir_vacuum, GS.m_ir_expansion),
                    vacuum_mass: holds(GS.m_uv_vacuum, GS.m_uv_expansion),
                    vacuum_expansion_mass: atom
                        .replace_multiple(&without_seeded_expansion_masses)
                        .contains_symbol(GS.m_uv_expansion),
                });
            }
            shapes.sort();
            Ok(shapes)
        }

        fn normalized_node_label(&self, node: NodeIndex) -> String {
            let label = self.node_label(node);
            let mut normalized = String::with_capacity(label.len());
            let mut chars = label.chars().peekable();

            while let Some(ch) = chars.next() {
                normalized.push(ch);
                if ch != '(' {
                    continue;
                }

                let mut digits = String::new();
                while chars.peek().is_some_and(|next| next.is_ascii_digit()) {
                    digits.push(chars.next().expect("peeked digit must exist"));
                }

                if digits.is_empty() || !matches!(chars.peek(), Some(')')) {
                    normalized.push_str(&digits);
                } else {
                    normalized.push('_');
                }
            }

            normalized
        }

        fn normalized_node_labels_with_cover(&self, cover_label: &str) -> Vec<String> {
            let mut labels = self
                .graph
                .iter_nodes()
                .filter(|(_, _, operation)| {
                    operation
                        .covers()
                        .is_some_and(|cover| cover.string_label() == cover_label)
                })
                .map(|(node, _, _)| self.normalized_node_label(node))
                .collect::<Vec<_>>();
            labels.sort();
            labels
        }
    }

    #[test]
    fn three_d_compute_skips_unused_four_d_atoms_when_integration_is_disabled() -> Result<()> {
        test_initialise().unwrap();
        for project_local_4d in [false, true] {
            let mut graph: Graph = dot!(
                digraph scalar_self_energy {
                    edge [particle="scalar_1"];
                    node [num=1];
                    ext [style=invis];
                    ext -> a;
                    b -> ext;
                    a -> b;
                    a -> b;
                },
                "scalars"
            )?;
            let settings = GenerationSettings {
                uv: UVgenerationSettings {
                    generate_integrated: false,
                    local_uv_cts_from_expanded_4d_integrands: project_local_4d,
                    ..Default::default()
                },
                ..Default::default()
            };
            let orientation_pattern = crate::settings::global::OrientationPattern::default();
            let options = graph.denominator_only_cff_3d_expression_options();
            let production = graph.generate_3d_expression_for_integrand(
                &[],
                &graph.get_esurface_canonization(&graph.loop_momentum_basis),
                &options,
                Some(&Atom::one()),
            )?;
            let mut forests = Wood::new(CutStructure::empty(&graph), &graph, &settings)?.unfold();

            let mut context = Local4dProjectionContext::default();
            forests.compute(
                &mut graph,
                crate::utils::vakint()?,
                OrientationProjection::exact_expression(
                    &production,
                    &options,
                    &orientation_pattern,
                    true,
                ),
                &settings.uv,
                true,
                &mut context,
            )?;

            if project_local_4d {
                let mut legacy =
                    crate::uv::wood::CutWoods::new(CutStructure::empty(&graph), &graph, &settings)?
                        .unfold(&graph);
                let mut legacy_context = Local4dProjectionContext::default();
                legacy.compute(
                    &mut graph,
                    crate::utils::vakint()?,
                    OrientationProjection::exact_expression(
                        &production,
                        &options,
                        &orientation_pattern,
                        true,
                    ),
                    &settings.uv,
                    true,
                    &mut legacy_context,
                )?;
                let canonical_builds = [
                    context.canonical_preparation_builds,
                    legacy_context.canonical_preparation_builds,
                ];
                assert!(canonical_builds.iter().all(|count| *count > 0));
                for representation in [
                    three_dimensional_reps::RepresentationMode::Ltd,
                    three_dimensional_reps::RepresentationMode::Cff,
                ] {
                    let options = three_dimensional_reps::Generate3DExpressionOptions {
                        representation,
                        ..options.clone()
                    };
                    let production = graph.generate_3d_expression_for_integrand(
                        &[],
                        &graph.get_esurface_canonization(&graph.loop_momentum_basis),
                        &options,
                        Some(&Atom::one()),
                    )?;
                    assert_eq!(production.representation, representation);
                    let orientation = OrientationProjection::exact_expression(
                        &production,
                        &options,
                        &orientation_pattern,
                        true,
                    );
                    let before = [&context, &legacy_context].map(|context| {
                        (
                            context.shared_preparation_time,
                            context.preparations.hits,
                            context.source_preparation_builds,
                            context.generation_cache.native_generations,
                        )
                    });
                    let shared = [
                        forests.compute(
                            &mut graph,
                            crate::utils::vakint()?,
                            orientation,
                            &settings.uv,
                            false,
                            &mut context,
                        )?,
                        legacy.compute(
                            &mut graph,
                            crate::utils::vakint()?,
                            orientation,
                            &settings.uv,
                            false,
                            &mut legacy_context,
                        )?,
                    ];
                    for (index, context) in [&context, &legacy_context].into_iter().enumerate() {
                        assert_eq!(
                            context.canonical_preparation_builds,
                            canonical_builds[index]
                        );
                        assert!(context.preparations.hits > before[index].1);
                        assert_eq!(
                            shared[index],
                            context.shared_preparation_time - before[index].0
                        );
                        if representation == three_dimensional_reps::RepresentationMode::Ltd {
                            assert!(context.source_preparation_builds > before[index].2);
                            assert!(context.generation_cache.native_generations > before[index].3);
                        } else {
                            assert_eq!(context.source_preparation_builds, before[index].2);
                            assert_eq!(
                                context.generation_cache.native_generations,
                                before[index].3
                            );
                        }
                    }
                }
            }

            assert_eq!(forests.compute_store.entries.len(), forests.graph.n_nodes());
            for (_, _, operation) in forests.graph.iter_nodes() {
                let computed = forests.compute_store.require(operation)?;
                assert_eq!(
                    computed.local_4d.is_some(),
                    project_local_4d,
                    "4D coefficients must be constructed exactly when the projected route needs them for {operation}"
                );
                let integrated = computed.integrated(operation)?;
                assert!(integrated.physical_pole_atom().is_zero());
                assert!(integrated.physical_finite_counterterm_atom().is_zero());
                assert!(
                    !computed.cuts.is_empty(),
                    "both routes must compute a cut for {operation}"
                );
                assert!(
                    computed
                        .cuts
                        .values()
                        .all(|cut| cut.final_integrands.iter().next().is_some()),
                    "every computed cut must retain a 3D final integrand for {operation}"
                );
            }
        }

        Ok(())
    }

    #[test]
    fn vacuum_subtraction_prunes_only_full_observable_spinneys() -> Result<()> {
        use crate::settings::global::MediumMode;

        test_initialise()?;
        let graphs: Vec<Graph> = dot!(
            digraph sunrise {
                edge [particle="scalar_1"];
                v1 -> v2; v1 -> v2; v1 -> v2;
            }
            digraph dumbbell {
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v2;
                v1 -> v1; v1 -> v1;
            }
            digraph sunrise_with_bridge {
                edge [particle="scalar_1"];
                ext [style=invis];
                v1 -> v2; v1 -> v2; v1 -> v2;
                v2 -> v3;
                ext -> v3; v3 -> ext;
            }, "scalars"
        )?;
        for graph in graphs {
            // The counterterm covers the loop part, leaving any bridge untouched.
            let full = InternalSubGraph::cleaned_filter_pessimist(
                graph.full_filter().subtract(&graph.tree_edges),
                graph.as_ref(),
            );
            for subtract_uv in [false, true] {
                let mut settings = GenerationSettings::default();
                settings.uv.subtract_uv = subtract_uv;
                let ordinary = Wood::new(CutStructure::empty(&graph), &graph, &settings)?;
                let all = ordinary
                    .graph
                    .iter_nodes()
                    .map(|(_, _, spinney)| spinney.filter().clone())
                    .collect::<BTreeSet<_>>();
                let proper = all
                    .iter()
                    .filter(|filter| **filter != full.filter)
                    .cloned()
                    .collect::<BTreeSet<_>>();
                assert!(proper.iter().any(|filter| filter.is_empty()));
                if subtract_uv {
                    assert!(all.contains(&full.filter), "{}", graph.name);
                    assert!(proper.len() > 1, "{}", graph.name);
                }
                for mode in [
                    MediumMode::Vacuum,
                    MediumMode::ThermodynamicEquilibrium,
                    MediumMode::ZeroTemperatureEquilibrium,
                ] {
                    settings.medium.mode = mode;
                    for vacuum_subtraction in [false, true] {
                        if mode == MediumMode::Vacuum && vacuum_subtraction {
                            continue; // Generation validation rejects this combination.
                        }
                        settings.medium.vacuum_subtraction = vacuum_subtraction;
                        let expected = if vacuum_subtraction { &proper } else { &all };
                        let wood = Wood::new(CutStructure::empty(&graph), &graph, &settings)?;
                        assert_eq!(
                            &wood
                                .graph
                                .iter_nodes()
                                .map(|(_, _, spinney)| spinney.filter().clone())
                                .collect::<BTreeSet<_>>(),
                            expected,
                            "{}: {mode:?}, subtract_uv={subtract_uv}, vacuum_subtraction={vacuum_subtraction}",
                            graph.name
                        );
                        // The shared filter also supplies the legacy executor.
                        let mut spinneys = graph.classified_spinneys(
                            &graph.full_filter(),
                            &settings.uv,
                            &graph.loop_momentum_basis,
                        )?;
                        graph.restrict_full_observable_spinney(&mut spinneys, &settings.medium)?;
                        assert_eq!(
                            &spinneys
                                .iter()
                                .map(|spinney| spinney.filter().clone())
                                .collect::<BTreeSet<_>>(),
                            expected,
                            "legacy {}: {mode:?}, subtract_uv={subtract_uv}, vacuum_subtraction={vacuum_subtraction}",
                            graph.name
                        );
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn vacuum_subtraction_commutes_with_local_and_integrated_forests() -> Result<()> {
        use crate::{
            integrands::process::param_builder::ParamBuilderGraph,
            numerator::{
                soft_energy_certificate::{EnergyZeroOutcome, SoftEnergyAlgebra},
                symbolica_ext::NumeratorAtomExt,
            },
            settings::global::MediumMode,
        };
        use symbolica::atom::AtomView;
        use three_dimensional_reps::thermal::OCCUPATION_FORM;
        test_initialise()?;
        let graphs: Vec<Graph> = dot!(
            digraph nested {
                edge [particle="scalar_1"]; node [num=1];
                v1 -> v2; v1 -> v2; v1 -> v2;
                v1 -> v1;
            }
            digraph disjoint {
                edge [particle="scalar_1"]; node [num=1];
                v1 -> v2; v1 -> v2;
                v3 -> v4; v3 -> v4;
                v2 -> v3; v1 -> v4;
            }, "scalars"
        )?;
        for graph in graphs {
            for medium in [
                MediumMode::ThermodynamicEquilibrium,
                MediumMode::ZeroTemperatureEquilibrium,
            ] {
                for generate_integrated in [false, true] {
                    let mut settings = GenerationSettings::default();
                    settings.medium.mode = medium;
                    settings.medium.vacuum_subtraction = true;
                    settings.uv.softct = false;
                    settings.uv.generate_integrated = generate_integrated;
                    let mut graph = graph.clone();
                    let options = graph.production_3d_expression_options(
                        &settings,
                        three_dimensional_reps::RepresentationMode::Cff,
                    )?;
                    let contract =
                        graph.paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut));
                    let production = graph
                        .generate_3d_expression_for_integrand(&contract, &None, &options, None)?;
                    let mut results = Vec::new();
                    for early_subtraction in [false, true] {
                        let projection = OrientationProjection::exact_expression(
                            &production,
                            &options,
                            &settings.orientation_pattern,
                            true,
                        )
                        .with_vacuum_subtraction(early_subtraction);
                        // Both routes use precisely the same proper forests. Only the
                        // boundary where the full vacuum projection acts is different.
                        let mut forests =
                            Wood::new(CutStructure::empty(&graph), &graph, &settings)?.unfold();
                        forests.compute(
                            &mut graph,
                            crate::utils::vakint()?,
                            projection,
                            &settings.uv,
                            true,
                            &mut Local4dProjectionContext::default(),
                        )?;
                        results.push(forests.orientation_parametric_exprs(&graph, &settings.uv)?);
                    }
                    assert_eq!(results[0].len(), results[1].len());
                    for (late, early) in results[0].iter().zip(&results[1]) {
                        let late = late.integrands.resolved()?;
                        let early = early.integrands.resolved()?;
                        assert_eq!(late.iter().count(), early.iter().count());
                        for ((late_key, raw), (early_key, actual)) in late.iter().zip(early.iter())
                        {
                            assert_eq!(late_key, early_key);
                            let vacuum = graph.make_thermal_distributions_explicit(
                                raw,
                                MediumMode::Vacuum,
                                graph.iter_edges().map(|(_, edge, _)| edge),
                            )?;
                            assert!(
                                raw.is_zero() || raw.contains_symbol(GS.thermal_weight_wrapper)
                            );
                            assert!(
                                actual.is_zero()
                                    || actual.contains_symbol(GS.thermal_weight_wrapper)
                            );
                            // Resolve only this scalar comparison. Final energy factors can
                            // group as E^-5 or E^-1/P^2, with E^2 = P, after the two routes.
                            let difference = (actual - (raw - vacuum))
                                .replace(function!(GS.thermal_weight_wrapper, W_.a_))
                                .with(W_.a_)
                                .replace_multiple(&*OCCUPATION_FORM)
                                .expand();
                            // Certify each occupation coefficient separately so unrelated
                            // energy denominators do not multiply the proof's term budget.
                            let mut occupations = BTreeSet::new();
                            difference.visitor(&mut |view| {
                                if let AtomView::Fun(call) = view
                                    && call.get_symbol() == GS.thermal_distribution
                                {
                                    occupations.insert(view.to_owned());
                                    return false;
                                }
                                true
                            });
                            let coefficients = difference
                                .coefficient_list_exact(&occupations.into_iter().collect_vec())
                                .expect("scalar difference must be polynomial in occupations");
                            for (occupation, coefficient) in coefficients {
                                let certificate =
                                    SoftEnergyAlgebra::default().scalar_zero(&coefficient)?;
                                assert!(
                                    matches!(certificate, EnergyZeroOutcome::Zero(_)),
                                    "{}: {medium:?}, integrated={generate_integrated}, residue={late_key:?}, occupation={occupation}: {certificate:?}: {coefficient}",
                                    graph.name
                                );
                            }
                        }
                    }
                }
            }
        }
        Ok(())
    }

    #[test]
    fn union_host_sits_above_the_factorized_join() -> Result<()> {
        test_initialise()?;
        // The mixed ring is the design case; the others share the union host.
        // The tadpole ring has three two-loop spinneys instead of four.
        for (ring, trace_nodes) in [
            (Ring::Mixed, 32),
            (Ring::Quark, 32),
            (Ring::Ghost, 32),
            (Ring::Tadpole, 28),
        ] {
            let graph = ring.physical_graph()?;
            let forests = thermal_ring_forests(&graph)?;
            let second = ring.second_loop_edges();
            let a = graph.edge_subgraph(&[2, 3]);
            let b = graph.edge_subgraph(second);
            let h = graph.edge_subgraph(&[[2, 3].as_slice(), second].concat());
            let g = graph.edge_subgraph(&(0..graph.underlying.n_edges()).collect::<Vec<_>>());
            let level = |scopes: &[&SuBitGraph]| {
                scopes
                    .iter()
                    .map(|scope| (*scope).clone())
                    .collect::<BTreeSet<_>>()
            };

            // The host sub-forest, each chain also under the whole-graph U:
            // the union host acts after the empty prefix, A, B and their
            // product.
            let traces = forests
                .graph
                .iter_nodes()
                .map(|(node, _, _)| trace_levels(&forests, node))
                .collect::<BTreeSet<_>>();
            let chains = [
                vec![level(&[&a])],
                vec![level(&[&b])],
                vec![level(&[&a, &b])],
                vec![level(&[&h])],
                vec![level(&[&a]), level(&[&h])],
                vec![level(&[&b]), level(&[&h])],
                vec![level(&[&a, &b]), level(&[&h])],
            ];
            for chain in chains {
                let mut enclosed = chain.clone();
                enclosed.push(level(&[&g]));
                assert!(traces.contains(&chain), "{ring:?}: missing {chain:?}");
                assert!(traces.contains(&enclosed), "{ring:?}: missing {enclosed:?}");
            }
            // The ordinary forests and the 8 chains that end on, or pass, the
            // host.
            assert_eq!(forests.graph.n_nodes(), trace_nodes, "{ring:?}: {forests}");
            let host_steps = forests
                .graph
                .iter_nodes()
                .filter_map(|(node, _, operation)| {
                    trace_levels(&forests, node)
                        .iter()
                        .any(|level| level.contains(&h))
                        .then_some((node, operation))
                })
                .collect::<Vec<_>>();
            assert_eq!(host_steps.len(), 8, "{ring:?}");
            for (node, operation) in host_steps {
                let source = forests.source_spinney(node);
                if source.filter() == &h {
                    assert_eq!(source.wood_role(), WoodRole::UnionHost, "{operation}");
                } else {
                    assert_eq!(source.filter(), &g, "{operation}");
                }
            }

            // The host is one collective operation: no join of its factors,
            // and every transition into it carries its whole scope.
            let wood = &forests.wood;
            let node_of = |filter: &SuBitGraph, role: WoodRole| {
                wood.graph
                    .iter_nodes()
                    .find_map(|(node, _, spinney)| {
                        (spinney.filter() == filter && spinney.wood_role() == role).then_some(node)
                    })
                    .expect("the wood has the node")
            };
            let host = node_of(&h, WoodRole::UnionHost);
            assert!(wood.join_factors(host).is_none(), "{ring:?}");
            assert_eq!(
                wood.join_factors(node_of(&h, WoodRole::Plain)),
                Some(BTreeSet::from([a.clone(), b.clone()])),
                "{ring:?}"
            );
            let incoming = wood
                .graph
                .iter_crown(host)
                .filter(|hedge| wood.graph.flow(*hedge) == Flow::Sink)
                .map(|hedge| {
                    assert_eq!(wood.graph[wood.graph[&hedge]], h, "{ring:?}");
                    let source = &wood.graph[wood
                        .graph
                        .involved_node_id(hedge)
                        .expect("a wood edge is paired")];
                    (source.filter().clone(), source.wood_role())
                })
                .collect::<BTreeSet<_>>();
            assert_eq!(
                incoming,
                BTreeSet::from([
                    (graph.empty_subgraph::<SuBitGraph>(), WoodRole::Plain),
                    (a, WoodRole::Plain),
                    (b, WoodRole::Plain),
                    (h, WoodRole::Plain),
                ]),
                "{ring:?}"
            );
            assert_eq!(
                forests.graph.validate_invariant::<SuBitGraph>(),
                Ok(()),
                "{ring:?}"
            );
        }
        Ok(())
    }

    #[test]
    fn factorized_join_replays_the_covering_host_first() -> Result<()> {
        test_initialise()?;
        for ring in [Ring::Mixed, Ring::Ghost] {
            let graph = ring.physical_graph()?;
            let forests = thermal_ring_forests(&graph)?;
            let a = graph.edge_subgraph(&[2, 3]);
            let b = graph.edge_subgraph(&[4, 5]);
            let join = ring_join(&forests, &graph.edge_subgraph(&[2, 3, 4, 5]));
            // The unmodified fixtures list B first, so the order is pinned
            // rather than inherited.
            let factors = forests
                .disconnected_component_nodes(join)?
                .into_iter()
                .map(|component| forests.source_spinney(component).filter().clone())
                .collect::<Vec<_>>();
            assert_eq!(factors, [b.clone(), a.clone()], "{ring:?}");

            // The soft set of the connected host A meets B: H_A replays before
            // U_B, which then acts on the zero mode of A.
            let scope = |edge: EdgeIndex| {
                forests
                    .wood
                    .current_given_pair(edge, 0)
                    .0
                    .subgraph()
                    .clone()
            };
            let replays = forests
                .union_replay_states(&graph, join)?
                .into_iter()
                .filter(|state| !state.local_edges.is_empty())
                .map(|state| {
                    (
                        trace_levels(&forests, state.integrated),
                        state
                            .local_edges
                            .iter()
                            .map(|edge| scope(*edge))
                            .collect::<Vec<_>>(),
                    )
                })
                .collect::<BTreeSet<_>>();
            assert_eq!(
                replays,
                BTreeSet::from([
                    (Vec::new(), vec![a.clone(), b.clone()]),
                    (vec![BTreeSet::from([b.clone()])], vec![a.clone()]),
                    (vec![BTreeSet::from([a.clone()])], vec![b.clone()]),
                ]),
                "{ring:?}"
            );
        }
        Ok(())
    }

    /// The kind of a local operation in a sector's projection path: an
    /// ordinary U, a zero-mode operation alone, or the operation `S0 + U - U S0`
    /// of a connected host.
    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
    enum StepKind {
        U,
        S0,
        H,
    }

    /// One sector of a node's local 3D counterterm with the
    /// zero-temperature lane seeded by symbols
    /// ([`Forests::seed_symbolic_integrated`]).
    #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
    struct SectorShape {
        /// Whether the sector is local, its first operation acting on the
        /// complete integrand rather than on an integrated prefix.
        local: bool,
        /// Whether the sector vanishes.
        vanishes: bool,
        /// The traces of the nodes whose seeds the sector holds: its
        /// integrated prefix, one trace per connected component.
        prefix: Vec<String>,
        /// Whether the sector holds the expansion-mass seed of every node of
        /// its prefix, the part of the integrated vertex with `mUVexp`.
        keeps_expansion_masses: bool,
        /// The local operations in order, each with its component and its
        /// rescaled support.
        steps: Vec<(StepKind, String, String)>,
        /// Whether the sector holds the static-theory mass `mIR`.
        static_mass: bool,
        /// Whether the sector holds the vacuum mass `mUV`.
        vacuum_mass: bool,
        /// Whether the sector holds the expansion mass `mUVexp` that a U
        /// introduces with the vacuum mass, besides that of the seeded
        /// integrated vertices.
        vacuum_expansion_mass: bool,
    }

    impl SectorShape {
        /// Whether a U follows a zero-mode operation in this sector's word,
        /// the only way it can hold the static-theory mass: a U after an
        /// `S0`, or after the operation of a connected host, whose rows of
        /// `S0` it expands as a static U.
        fn has_u_after_zero_mode(&self) -> bool {
            self.steps
                .iter()
                .skip_while(|(kind, _, _)| *kind == StepKind::U)
                .skip(1)
                .any(|(kind, _, _)| *kind != StepKind::S0)
        }

        /// Whether this sector's word starts with a zero-mode operation alone,
        /// after which every row holds a completed zero-mode branch.
        fn starts_on_zero_mode(&self) -> bool {
            self.steps
                .first()
                .is_some_and(|(kind, _, _)| *kind == StepKind::S0)
        }

        /// The components of the operations of the integrated prefix.
        fn prefix_components(&self) -> BTreeSet<&str> {
            self.prefix
                .iter()
                .flat_map(|trace| trace.split([';', ',']))
                .collect()
        }

        /// Whether a local sector `partner` is the word of this integrated
        /// sector with the U's of the prefix in place: without the
        /// operations of the prefix's components, its steps are this sector's
        /// suffix, the self-suffix `S0` of a connected host on the empty
        /// support standing for the host's own U, which the partner keeps.
        /// A connected host whose plain U alone replays after a prefix in its
        /// soft set matches the host's operation.
        fn is_integrated_form_of(&self, partner: &Self) -> bool {
            let components = self.prefix_components();
            let suffix = self
                .steps
                .iter()
                .filter(|(kind, component, support)| {
                    !(*kind == StepKind::S0
                        && support == "∅"
                        && components.contains(component.as_str()))
                })
                .collect::<Vec<_>>();
            let local = partner
                .steps
                .iter()
                .filter(|(_, component, _)| !components.contains(component.as_str()))
                .collect::<Vec<_>>();
            partner.local
                && suffix.len() == local.len()
                && suffix.iter().zip(&local).all(
                    |((kind, component, _), (partner_kind, partner_component, _))| {
                        component == partner_component
                            && (kind == partner_kind
                                || (*kind == StepKind::U && *partner_kind == StepKind::H))
                    },
                )
        }
    }

    impl Display for SectorShape {
        fn fmt(&self, f: &mut Formatter<'_>) -> FmtResult {
            let prefix = self
                .prefix
                .iter()
                .map(|trace| format!("∫[{trace}]"))
                .join(" ");
            let steps = self
                .steps
                .iter()
                .map(|(kind, component, support)| format!("{kind:?}_{component}|{support}"))
                .join(", ");
            let masses = [(self.static_mass, "mIR"), (self.vacuum_mass, "mUV")]
                .into_iter()
                .filter_map(|(held, name)| held.then_some(name))
                .join(",");
            write!(f, "({prefix}; {steps}) [{masses}]")
        }
    }

    /// The symbol `c<node>` that stands for the finite counterterm of `node`
    /// in [`Forests::seed_symbolic_integrated`].
    fn seed_symbol(node: NodeIndex) -> Symbol {
        symbol!(&format!("hedge_poset_test_c{}", node.0))
    }

    /// The symbol `e<node>` of the part with the expansion mass `mUVexp` of
    /// the finite counterterm of `node` in
    /// [`Forests::seed_symbolic_integrated`].
    fn expansion_seed_symbol(node: NodeIndex) -> Symbol {
        symbol!(&format!("hedge_poset_test_e{}", node.0))
    }

    /// The integrated-route settings of the zero-mode pins on the ring
    /// topologies.
    fn ring_topology_settings(vacuum_subtraction: bool) -> GenerationSettings {
        GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: true,
                ..ring_uv_settings(ApproximationType::MUV)
            },
            medium: thermal_medium(vacuum_subtraction),
            ..Default::default()
        }
    }

    /// The spinneys of the host policy on the physical `ring`, with its
    /// schemes and degrees, transplanted onto `graph`, the topology of `ring`
    /// with unit numerators, whose expansions and integrals are cheap. The
    /// hosts carry their regions on the topology ([`Ring::expected_regions`]).
    /// Returns the forest of plain U's on the same subgraphs and the forest
    /// of the host policy.
    fn ring_topology_spinneys(
        ring: Ring,
        graph: &Graph,
        settings: &GenerationSettings,
    ) -> Result<(Vec<Spinney>, Vec<Spinney>)> {
        let physical = ring.physical_graph()?;
        let policy = zero_mode_policy(
            &physical,
            &UVgenerationSettings {
                generate_integrated: true,
                ..ring_uv_settings(ApproximationType::IR)
            },
            &settings.medium,
        )?;
        let regions = ring.expected_regions(graph, settings.medium.vacuum_subtraction);
        let transplant = |spinney: &Spinney, scheme| -> Result<Spinney> {
            let edges = physical
                .paired_edges(spinney.filter())
                .into_iter()
                .map(usize::from)
                .collect::<Vec<_>>();
            if edges.is_empty() {
                return Ok(Spinney::empty(graph));
            }
            let subgraph = InternalSubGraph::try_new(graph.edge_subgraph(&edges), graph.as_ref())
                .ok_or_else(|| eyre!("{ring:?}: {edges:?} is no internal subgraph"))?;
            Spinney::with_scheme(
                subgraph,
                graph,
                &graph.loop_momentum_basis,
                scheme,
                spinney.dod,
            )?
            .ok_or_else(|| {
                eyre!(
                    "{ring:?}: {edges:?} is no spinney of degree {}",
                    spinney.dod
                )
            })
        };
        let mut plain = Vec::new();
        let mut hosted = Vec::new();
        for spinney in &policy {
            let transplanted = transplant(spinney, spinney.renormalization_scheme)?;
            if let Some(region) = spinney.thermal_zero_mode() {
                let edges = physical
                    .paired_edges(spinney.filter())
                    .into_iter()
                    .collect::<BTreeSet<_>>();
                assert_eq!(regions[&edges].kind(), region.kind(), "{ring:?}");
                // The plain join of a union host is a spinney of its own.
                if region.kind() == ZeroModeHostKind::Connected {
                    plain.push(transplant(spinney, ApproximationType::MUV)?);
                }
                hosted.push(transplanted.with_thermal_zero_mode(regions[&edges].clone()));
            } else {
                plain.push(transplant(spinney, ApproximationType::MUV)?);
                hosted.push(transplanted);
            }
        }
        assert_eq!(hosted.len(), plain.len() + 1, "{ring:?}: one union host");
        Ok((plain, hosted))
    }

    /// The sector shapes of every node of the physical `ring` on the
    /// integrated route under the host policy, with the zero-temperature lane
    /// seeded by symbols ([`Forests::seed_symbolic_integrated`]), keyed by
    /// the node's trace ([`Forests::trace_name`]).
    fn symbolic_integrated_inventory(
        ring: Ring,
        vacuum_subtraction: bool,
    ) -> Result<BTreeMap<String, Vec<SectorShape>>> {
        let mut graph = ring.physical_graph()?;
        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: true,
                ..ring_uv_settings(ApproximationType::IR)
            },
            medium: thermal_medium(vacuum_subtraction),
            ..Default::default()
        };
        let options = graph.production_3d_expression_options(&settings, RepresentationMode::Cff)?;
        let contract = graph.paired_edges(&graph.tree_edges.subtract(&graph.initial_state_cut));
        let production =
            graph.generate_3d_expression_for_integrand(&contract, &None, &options, None)?;
        let projection = OrientationProjection::exact_expression(
            &production,
            &options,
            &settings.orientation_pattern,
            true,
        )
        .with_vacuum_subtraction(vacuum_subtraction);
        let mut forests = Wood::new(CutStructure::empty(&graph), &graph, &settings)?.unfold();
        forests.seed_symbolic_integrated(&graph, vacuum_subtraction)?;
        forests.compute(
            &mut graph,
            utils::vakint()?,
            projection,
            &settings.uv,
            false,
            &mut Local4dProjectionContext::default(),
        )?;
        let names = forests.subgraph_names(&graph);
        forests
            .graph
            .iter_nodes()
            .map(|(node, _, _)| {
                Ok((
                    forests.trace_name(node, &names),
                    forests.sector_shapes(node, &names)?,
                ))
            })
            .collect()
    }

    /// The zero-temperature integrated lane counts a connected zero-mode host
    /// by its plain U and gives nothing to a union host and to every node
    /// above one: their zero-mode operations are local only. Under vacuum
    /// subtraction it gives nothing to a node covering every loop either,
    /// whose add-back would sit on a loopless cograph. On the ring topologies
    /// with unit numerators, whose integrals are cheap, with the hosts of the
    /// host policy ([`ring_topology_spinneys`]), every other node has the
    /// local and integrated counterterms of the same trace in the forest of
    /// plain U's, the factorized join the product of its factors'.
    #[test]
    fn zero_mode_hosts_integrate_their_plain_u() -> Result<()> {
        test_initialise()?;
        for ring in Ring::ALL {
            let graph = ring.topology()?;
            let settings = ring_topology_settings(true);
            let (plain, hosted) = ring_topology_spinneys(ring, &graph, &settings)?;
            let vakint = utils::vakint()?;
            let mut forests =
                Wood::from_spinneys(hosted, &graph, CutStructure::empty(&graph), &settings.uv)?
                    .unfold();
            let mut reference =
                Wood::from_spinneys(plain, &graph, CutStructure::empty(&graph), &settings.uv)?
                    .unfold();
            forests.integrate(&graph, vakint, &settings.uv, true)?;
            reference.integrate(&graph, vakint, &settings.uv, true)?;

            let union =
                graph.edge_subgraph(&[[2, 3].as_slice(), ring.second_loop_edges()].concat());
            let a = graph.edge_subgraph(&[2, 3]);
            let mut compared = BTreeMap::new();
            for (node, _, operation) in forests.graph.iter_nodes() {
                if operation.key.is_empty() {
                    continue;
                }
                let levels = trace_levels(&forests, node);
                let context = format!("{ring:?}: {}", forests.node_label(node));
                let computed = forests.compute_store.require(operation)?;
                let local = computed.local_4d(operation)?.atom().clone();
                let integrated = computed.integrated(operation)?;
                let above_union_host = levels.iter().any(|level| level.contains(&union));
                let covers_every_loop = graph.n_loops(&operation.covers().expect("a proper cover"))
                    == graph.n_loops(&graph.full_filter());
                assert_eq!(
                    forests.has_zero_temperature_counterpart(node, &graph, false),
                    !above_union_host,
                    "{context}: without vacuum subtraction"
                );
                if above_union_host || covers_every_loop {
                    assert!(local.is_zero(), "{context}: {local}");
                    assert!(
                        integrated.physical_finite_counterterm_atom().is_zero()
                            && integrated.physical_pole_atom().is_zero(),
                        "{context}"
                    );
                    continue;
                }
                let (twin, _, twin_operation) = reference
                    .graph
                    .iter_nodes()
                    .find(|(twin, _, _)| trace_levels(&reference, *twin) == levels)
                    .ok_or_else(|| eyre!("{context}: no plain twin"))?;
                let twin_computed = reference.compute_store.require(twin_operation)?;
                let twin_integrated = twin_computed.integrated(twin_operation)?;
                assert_eq!(
                    &local,
                    twin_computed.local_4d(twin_operation)?.atom(),
                    "{context}: local, against {}",
                    reference.node_label(twin)
                );
                assert_eq!(
                    integrated.physical_finite_counterterm_atom(),
                    twin_integrated.physical_finite_counterterm_atom(),
                    "{context}: finite counterterm"
                );
                assert_eq!(
                    integrated.physical_pole_atom(),
                    twin_integrated.physical_pole_atom(),
                    "{context}: pole"
                );
                compared.insert(levels, integrated.physical_finite_counterterm_atom());
            }
            // The host itself, its bubble partner and their join integrate to
            // nonzero counterterms.
            let b = graph.edge_subgraph(ring.second_loop_edges());
            for levels in [
                vec![BTreeSet::from([a.clone()])],
                vec![BTreeSet::from([b.clone()])],
                vec![BTreeSet::from([a.clone(), b.clone()])],
            ] {
                assert!(
                    compared
                        .get(&levels)
                        .is_some_and(|finite| !finite.is_zero()),
                    "{ring:?}: {levels:?} is not compared or vanishes"
                );
            }
        }
        Ok(())
    }

    /// Checks the sectors of the physical `ring` on the integrated route,
    /// with or without `vacuum_subtraction`, against the integrated sectors
    /// of the forest words in `expected`: one line
    /// `[trace] (∫[prefix trace] ...; operations) [masses]` per nonzero
    /// sector with an integrated prefix
    /// ([`symbolic_integrated_inventory`]).
    ///
    /// It also checks that the add-back of every U with the vacuum mass
    /// cancels its mass per counterterm family: every such sector has its
    /// local partner at its node, the same word with the prefix's U's in
    /// place ([`SectorShape::is_integrated_form_of`]). No U after a zero-mode
    /// operation introduces the vacuum mass: only a sector with a U after one
    /// holds `mIR`, and a sector whose word starts with a zero-mode operation
    /// alone, all of whose rows hold a completed branch, holds no `mUVexp`.
    fn check_integrated_sectors(
        ring: Ring,
        vacuum_subtraction: bool,
        expected: &str,
    ) -> Result<()> {
        test_initialise()?;
        let mut integrated = Vec::new();
        for (trace, shapes) in symbolic_integrated_inventory(ring, vacuum_subtraction)? {
            for shape in shapes.iter().filter(|shape| !shape.vanishes) {
                let context = format!(
                    "{ring:?}, vacuum subtraction {vacuum_subtraction}, [{trace}]: {shape}"
                );
                if !shape.local {
                    assert!(!shape.prefix.is_empty(), "{context}: no integrated prefix");
                    assert!(
                        shapes
                            .iter()
                            .any(|partner| shape.is_integrated_form_of(partner)),
                        "{context}: no local partner"
                    );
                    integrated.push(format!("[{trace}] {shape}"));
                }
                assert!(
                    !shape.static_mass || shape.has_u_after_zero_mode(),
                    "{context}: mIR without a U after a zero-mode operation"
                );
                assert!(
                    !shape.starts_on_zero_mode() || !shape.vacuum_expansion_mass,
                    "{context}: a U after a zero-mode operation introduces mUV"
                );
                assert!(
                    !shape.starts_on_zero_mode() || shape.keeps_expansion_masses,
                    "{context}: the zero-mode image drops the expansion mass of its integrated \
                     vertex, and with it part of the vertex's static mode"
                );
            }
        }
        integrated.sort();
        let expected = expected
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .collect::<Vec<_>>();
        assert_eq!(
            integrated,
            expected,
            "{ring:?}, vacuum subtraction {vacuum_subtraction}: the integrated sectors are\n{}",
            integrated.join("\n")
        );
        Ok(())
    }

    #[test]
    fn quark_ring_integrated_sectors_follow_the_forest_words() -> Result<()> {
        check_integrated_sectors(
            Ring::Quark,
            true,
            r"
            [A,B;AB;G] (∫[A]; U_B|B, S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A,B;AB;G] (∫[B] ∫[A]; S0_AB|∅, U_G|G-AB) [mIR,mUV]
            [A,B;AB;G] (∫[B]; U_A|A, S0_AB|AB-B, U_G|G-B) [mIR,mUV]
            [A,B;AB] (∫[A]; U_B|B, S0_AB|AB-A) [mUV]
            [A,B;AB] (∫[B] ∫[A]; S0_AB|∅) [mUV]
            [A,B;AB] (∫[B]; U_A|A, S0_AB|AB-B) [mUV]
            [A,B] (∫[A]; U_B|B) [mUV]
            [A,B] (∫[B]; U_A|A) [mUV]
            [A;AB;G] (∫[A]; S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A;AB] (∫[A]; S0_AB|AB-A) [mUV]
            [A;D01234] (∫[A]; U_D01234|D01234-A) [mUV]
            [A;D01235] (∫[A]; U_D01235|D01235-A) [mUV]
            [B;AB;G] (∫[B]; S0_AB|AB-B, U_G|G-B) [mIR,mUV]
            [B;AB] (∫[B]; S0_AB|AB-B) [mUV]
            [B;D01245] (∫[B]; U_D01245|D01245-B) [mUV]
            [B;D01345] (∫[B]; U_D01345|D01345-B) [mUV]
            ",
        )
    }

    #[test]
    fn quark_ring_integrated_sectors_without_vacuum_subtraction() -> Result<()> {
        check_integrated_sectors(
            Ring::Quark,
            false,
            r"
            [A,B;AB;G] (∫[A]; U_B|B, S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A,B;AB;G] (∫[B]; U_A|A, S0_AB|AB-B, U_G|G-B) [mIR,mUV]
            [A,B;AB] (∫[A]; U_B|B, S0_AB|AB-A) [mUV]
            [A,B;AB] (∫[B] ∫[A]; S0_AB|∅) [mUV]
            [A,B;AB] (∫[B]; U_A|A, S0_AB|AB-B) [mUV]
            [A,B;G] (∫[A]; U_B|B, U_G|G-A) [mUV]
            [A,B;G] (∫[B] ∫[A]; U_G|G-AB) [mUV]
            [A,B;G] (∫[B]; U_A|A, U_G|G-B) [mUV]
            [A,B] (∫[A]; U_B|B) [mUV]
            [A,B] (∫[B]; U_A|A) [mUV]
            [A;AB;G] (∫[A]; S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A;AB] (∫[A]; S0_AB|AB-A) [mUV]
            [A;D01234;G] (∫[A;D01234]; U_G|G-D01234) [mUV]
            [A;D01234;G] (∫[A]; U_D01234|D01234-A, U_G|G-A) [mUV]
            [A;D01234] (∫[A]; U_D01234|D01234-A) [mUV]
            [A;D01235;G] (∫[A;D01235]; U_G|G-D01235) [mUV]
            [A;D01235;G] (∫[A]; U_D01235|D01235-A, U_G|G-A) [mUV]
            [A;D01235] (∫[A]; U_D01235|D01235-A) [mUV]
            [A;G] (∫[A]; U_G|G-A) [mUV]
            [B;AB;G] (∫[B]; S0_AB|AB-B, U_G|G-B) [mIR,mUV]
            [B;AB] (∫[B]; S0_AB|AB-B) [mUV]
            [B;D01245;G] (∫[B;D01245]; U_G|G-D01245) [mUV]
            [B;D01245;G] (∫[B]; U_D01245|D01245-B, U_G|G-B) [mUV]
            [B;D01245] (∫[B]; U_D01245|D01245-B) [mUV]
            [B;D01345;G] (∫[B;D01345]; U_G|G-D01345) [mUV]
            [B;D01345;G] (∫[B]; U_D01345|D01345-B, U_G|G-B) [mUV]
            [B;D01345] (∫[B]; U_D01345|D01345-B) [mUV]
            [B;G] (∫[B]; U_G|G-B) [mUV]
            [D01234;G] (∫[D01234]; U_G|G-D01234) [mUV]
            [D01235;G] (∫[D01235]; U_G|G-D01235) [mUV]
            [D01245;G] (∫[D01245]; U_G|G-D01245) [mUV]
            [D01345;G] (∫[D01345]; U_G|G-D01345) [mUV]
            ",
        )
    }

    #[test]
    fn tadpole_ring_integrated_sectors_follow_the_forest_words() -> Result<()> {
        check_integrated_sectors(
            Ring::Tadpole,
            true,
            r"
            [A,B;AB;G] (∫[A]; U_B|B, S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A,B;AB;G] (∫[B] ∫[A]; S0_AB|∅, U_G|G-AB) [mIR,mUV]
            [A,B;AB;G] (∫[B]; U_A|A, S0_AB|AB-B, U_G|G-B) [mIR,mUV]
            [A,B;AB] (∫[A]; U_B|B, S0_AB|AB-A) [mUV]
            [A,B;AB] (∫[B] ∫[A]; S0_AB|∅) [mUV]
            [A,B;AB] (∫[B]; U_A|A, S0_AB|AB-B) [mUV]
            [A,B] (∫[A]; U_B|B) [mUV]
            [A,B] (∫[B]; U_A|A) [mUV]
            [A;AB;G] (∫[A]; S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A;AB] (∫[A]; S0_AB|AB-A) [mUV]
            [A;D0123] (∫[A]; U_D0123|D0123-A) [mUV]
            [B;AB;G] (∫[B]; S0_AB|AB-B, U_G|G-B) [mIR,mUV]
            [B;AB] (∫[B]; S0_AB|AB-B) [mUV]
            [B;D0124] (∫[B]; U_D0124|D0124-B) [mUV]
            [B;D0134] (∫[B]; U_D0134|D0134-B) [mUV]
            ",
        )
    }

    #[test]
    fn tadpole_ring_integrated_sectors_without_vacuum_subtraction() -> Result<()> {
        check_integrated_sectors(
            Ring::Tadpole,
            false,
            r"
            [A,B;AB;G] (∫[B]; U_A|A, S0_AB|AB-B, U_G|G-B) [mIR,mUV]
            [A,B;AB] (∫[A]; U_B|B, S0_AB|AB-A) [mUV]
            [A,B;AB] (∫[B] ∫[A]; S0_AB|∅) [mUV]
            [A,B;AB] (∫[B]; U_A|A, S0_AB|AB-B) [mUV]
            [A,B;G] (∫[A]; U_B|B, U_G|G-A) [mUV]
            [A,B;G] (∫[B] ∫[A]; U_G|G-AB) [mUV]
            [A,B;G] (∫[B]; U_A|A, U_G|G-B) [mUV]
            [A,B] (∫[A]; U_B|B) [mUV]
            [A,B] (∫[B]; U_A|A) [mUV]
            [A;AB] (∫[A]; S0_AB|AB-A) [mUV]
            [A;D0123;G] (∫[A;D0123]; U_G|G-D0123) [mUV]
            [A;D0123;G] (∫[A]; U_D0123|D0123-A, U_G|G-A) [mUV]
            [A;D0123] (∫[A]; U_D0123|D0123-A) [mUV]
            [A;G] (∫[A]; U_G|G-A) [mUV]
            [B;AB;G] (∫[B]; S0_AB|AB-B, U_G|G-B) [mIR,mUV]
            [B;AB] (∫[B]; S0_AB|AB-B) [mUV]
            [B;D0124;G] (∫[B;D0124]; U_G|G-D0124) [mUV]
            [B;D0124;G] (∫[B]; U_D0124|D0124-B, U_G|G-B) [mUV]
            [B;D0124] (∫[B]; U_D0124|D0124-B) [mUV]
            [B;D0134;G] (∫[B;D0134]; U_G|G-D0134) [mUV]
            [B;D0134;G] (∫[B]; U_D0134|D0134-B, U_G|G-B) [mUV]
            [B;D0134] (∫[B]; U_D0134|D0134-B) [mUV]
            [B;G] (∫[B]; U_G|G-B) [mUV]
            [D0123;G] (∫[D0123]; U_G|G-D0123) [mUV]
            [D0124;G] (∫[D0124]; U_G|G-D0124) [mUV]
            [D0134;G] (∫[D0134]; U_G|G-D0134) [mUV]
            ",
        )
    }

    /// The gluon ring has the forest of the ghost ring, but its local
    /// counterterms take longer to build than a unit test may run.
    #[test]
    fn ghost_ring_integrated_sectors_follow_the_forest_words() -> Result<()> {
        check_integrated_sectors(
            Ring::Ghost,
            true,
            r"
            [A,B;AB;G] (∫[A]; S0_A|∅, U_B|B, S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A,B;AB;G] (∫[A]; U_B|B, S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A,B;AB;G] (∫[B] ∫[A]; S0_AB|∅, U_G|G-AB) [mIR,mUV]
            [A,B;AB;G] (∫[B]; U_A|A, S0_AB|AB-B, U_G|G-B) [mIR,mUV]
            [A,B;AB] (∫[A]; S0_A|∅, U_B|B, S0_AB|AB-A) [mIR,mUV]
            [A,B;AB] (∫[A]; U_B|B, S0_AB|AB-A) [mUV]
            [A,B;AB] (∫[B] ∫[A]; S0_AB|∅) [mUV]
            [A,B;AB] (∫[B]; U_A|A, S0_AB|AB-B) [mUV]
            [A,B;G] (∫[A]; S0_A|∅, U_B|B, U_G|G-A) [mIR,mUV]
            [A,B] (∫[A]; S0_A|∅, U_B|B) [mIR,mUV]
            [A,B] (∫[A]; U_B|B) [mUV]
            [A,B] (∫[B]; U_A|A) [mUV]
            [A;AB;G] (∫[A]; S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A;AB;G] (∫[A]; S0_A|∅, S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A;AB] (∫[A]; S0_AB|AB-A) [mUV]
            [A;AB] (∫[A]; S0_A|∅, S0_AB|AB-A) [mUV]
            [A;D01234] (∫[A]; U_D01234|D01234-A) [mUV]
            [A;D01235] (∫[A]; U_D01235|D01235-A) [mUV]
            [A;G] (∫[A]; S0_A|∅, U_G|G-A) [mIR,mUV]
            [B;AB;G] (∫[B]; S0_AB|AB-B, U_G|G-B) [mIR,mUV]
            [B;AB] (∫[B]; S0_AB|AB-B) [mUV]
            [B;D01245] (∫[B]; U_D01245|D01245-B) [mUV]
            [B;D01345] (∫[B]; U_D01345|D01345-B) [mUV]
            ",
        )
    }

    #[test]
    fn ghost_ring_integrated_sectors_without_vacuum_subtraction() -> Result<()> {
        check_integrated_sectors(
            Ring::Ghost,
            false,
            r"
            [A,B;AB;G] (∫[A]; S0_A|∅, U_B|B, S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A,B;AB;G] (∫[A]; U_B|B, S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A,B;AB;G] (∫[B]; U_A|A, S0_AB|AB-B, U_G|G-B) [mIR,mUV]
            [A,B;AB] (∫[A]; S0_A|∅, U_B|B, S0_AB|AB-A) [mIR,mUV]
            [A,B;AB] (∫[A]; U_B|B, S0_AB|AB-A) [mUV]
            [A,B;AB] (∫[B] ∫[A]; S0_AB|∅) [mUV]
            [A,B;AB] (∫[B]; U_A|A, S0_AB|AB-B) [mUV]
            [A,B;G] (∫[A]; S0_A|∅, U_B|B, U_G|G-A) [mIR,mUV]
            [A,B;G] (∫[A]; U_B|B, U_G|G-A) [mUV]
            [A,B;G] (∫[B] ∫[A]; U_G|G-AB) [mUV]
            [A,B;G] (∫[B]; U_A|A, U_G|G-B) [mUV]
            [A,B] (∫[A]; S0_A|∅, U_B|B) [mIR,mUV]
            [A,B] (∫[A]; U_B|B) [mUV]
            [A,B] (∫[B]; U_A|A) [mUV]
            [A;AB;G] (∫[A]; S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A;AB;G] (∫[A]; S0_A|∅, S0_AB|AB-A, U_G|G-A) [mIR,mUV]
            [A;AB] (∫[A]; S0_AB|AB-A) [mUV]
            [A;AB] (∫[A]; S0_A|∅, S0_AB|AB-A) [mUV]
            [A;D01234;G] (∫[A;D01234]; U_G|G-D01234) [mUV]
            [A;D01234;G] (∫[A]; U_D01234|D01234-A, U_G|G-A) [mUV]
            [A;D01234] (∫[A]; U_D01234|D01234-A) [mUV]
            [A;D01235;G] (∫[A;D01235]; U_G|G-D01235) [mUV]
            [A;D01235;G] (∫[A]; U_D01235|D01235-A, U_G|G-A) [mUV]
            [A;D01235] (∫[A]; U_D01235|D01235-A) [mUV]
            [A;G] (∫[A]; S0_A|∅, U_G|G-A) [mIR,mUV]
            [A;G] (∫[A]; U_G|G-A) [mUV]
            [B;AB;G] (∫[B]; S0_AB|AB-B, U_G|G-B) [mIR,mUV]
            [B;AB] (∫[B]; S0_AB|AB-B) [mUV]
            [B;D01245;G] (∫[B;D01245]; U_G|G-D01245) [mUV]
            [B;D01245;G] (∫[B]; U_D01245|D01245-B, U_G|G-B) [mUV]
            [B;D01245] (∫[B]; U_D01245|D01245-B) [mUV]
            [B;D01345;G] (∫[B;D01345]; U_G|G-D01345) [mUV]
            [B;D01345;G] (∫[B]; U_D01345|D01345-B, U_G|G-B) [mUV]
            [B;D01345] (∫[B]; U_D01345|D01345-B) [mUV]
            [B;G] (∫[B]; U_G|G-B) [mUV]
            [D01234;G] (∫[D01234]; U_G|G-D01234) [mUV]
            [D01235;G] (∫[D01235]; U_G|G-D01235) [mUV]
            [D01245;G] (∫[D01245]; U_G|G-D01245) [mUV]
            [D01345;G] (∫[D01345]; U_G|G-D01345) [mUV]
            ",
        )
    }

    #[test]
    fn mutually_covering_hosts_in_a_join_are_rejected() -> Result<()> {
        test_initialise()?;
        let graph = Ring::Mixed.physical_graph()?;
        let settings = ring_uv_settings(ApproximationType::IR);
        let b = graph.edge_subgraph(&[4, 5]);
        // A synthetic connected host on B whose soft set meets A, as the soft
        // set of the host A meets B: neither factor can replay first.
        let spinneys = zero_mode_policy(&graph, &settings, &thermal_medium(true))?
            .into_iter()
            .map(|spinney| {
                if spinney.filter() == &b {
                    spinney.with_thermal_zero_mode(connected_region(&[0, 1, 2, 3]))
                } else {
                    spinney
                }
            })
            .collect::<Vec<_>>();
        let forests =
            Wood::from_spinneys(spinneys, &graph, CutStructure::empty(&graph), &settings)?.unfold();
        let join = ring_join(&forests, &graph.edge_subgraph(&[2, 3, 4, 5]));
        let Err(error) = forests.union_replay_states(&graph, join) else {
            panic!("mutually covering hosts must not replay");
        };
        let Some(ZeroModeHostError::MutuallyCoveringJoin { factors }) =
            error.root_cause().downcast_ref::<ZeroModeHostError>()
        else {
            panic!("unexpected error: {error:#}");
        };
        assert_eq!(
            factors.iter().cloned().collect::<BTreeSet<_>>(),
            BTreeSet::from([edge_set(&[2, 3]), edge_set(&[4, 5])])
        );
        assert!(
            format!("{error:#}").contains("incompatible zero-mode hosts in a factorized join"),
            "{error:#}"
        );
        Ok(())
    }

    #[test]
    fn collective_host_scope_cannot_equal_a_collective_union_scope() -> Result<()> {
        test_initialise()?;
        let graph = Ring::Mixed.physical_graph()?;
        let settings = ring_uv_settings(ApproximationType::IR);
        let b = graph.edge_subgraph(&[4, 5]);
        // Without B, the plain union {A, B} is a collective region whose
        // trace operations carry the scope of the union host.
        let spinneys = zero_mode_policy(&graph, &settings, &thermal_medium(true))?
            .into_iter()
            .filter(|spinney| spinney.filter() != &b)
            .collect::<Vec<_>>();
        let Err(error) =
            Wood::from_spinneys(spinneys, &graph, CutStructure::empty(&graph), &settings)
        else {
            panic!("a union host must not share the scope of a collective union");
        };
        assert!(
            matches!(
                error.root_cause().downcast_ref::<ZeroModeHostError>(),
                Some(ZeroModeHostError::CollectiveUnionScope { host })
                    if *host == edge_set(&[2, 3, 4, 5])
            ),
            "{error:#}"
        );
        Ok(())
    }

    #[test]
    fn local_leaf_operations_follow_dependency_frontiers() -> Result<()> {
        test_initialise().unwrap();
        let dumbell: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v2;
                v1 -> v1;v1 -> v1;
            },"scalars"
        )?;

        let forests = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &GenerationSettings::default(),
        )?
        .unfold();

        let frontier_label = |operation: &OperationNode| {
            operation
                .covers()
                .map_or_else(|| "∅".to_string(), |cover| cover.string_label())
        };
        let leaf_labels = |node| {
            forests
                .local_leaf_operations(node)
                .iter()
                .map(|leaf| {
                    (
                        leaf.op.order.string_label(),
                        frontier_label(&forests.graph[leaf.frontier]),
                    )
                })
                .collect::<Vec<_>>()
        };

        let root_disconnected = forests
            .graph
            .iter_nodes()
            .find_map(|(node, _, _)| {
                forests
                    .normalized_node_label(node)
                    .starts_with("{36,F}:")
                    .then_some(node)
            })
            .expect("lopsided dumbbell should contain a root disconnected frontier");
        assert_eq!(
            leaf_labels(root_disconnected),
            vec![
                ("36".to_string(), "∅".to_string()),
                ("F".to_string(), "∅".to_string())
            ]
        );
        let structured_dot = forests.dot_serialize();
        let mut join_markers = Vec::new();
        for (node, _, _) in forests.graph.iter_nodes() {
            for leaf in forests.local_leaf_operations(node) {
                let frontier_depth = forests.graph[leaf.frontier].key.op_count();
                let (current, given) = forests
                    .wood
                    .current_given_pair(leaf.op.data, frontier_depth);
                if current.subgraph() == &leaf.op.order {
                    continue;
                }
                join_markers.push(
                    UvMarker::subgraph(current.subgraph(), given.subgraph())
                        .printer(SpensoPrintSettings::typst_options())
                        .to_string(),
                );
            }
        }
        assert!(
            !join_markers.is_empty(),
            "lopsided dumbbell must contain a disconnected join"
        );
        for marker in join_markers {
            assert!(
                structured_dot.contains(&marker),
                "structured DOT must use the runtime current/given marker {marker}"
            );
        }

        let dependent_disconnected = forests
            .graph
            .iter_nodes()
            .find_map(|(node, _, _)| {
                forests
                    .normalized_node_label(node)
                    .starts_with("{C} · {36,F}:")
                    .then_some(node)
            })
            .expect("lopsided dumbbell should contain a C-dependent disconnected frontier");
        assert_eq!(
            leaf_labels(dependent_disconnected),
            vec![
                ("F".to_string(), "C".to_string()),
                ("36".to_string(), "∅".to_string())
            ]
        );
        let replay_states = forests.union_replay_states(&dumbell, dependent_disconnected)?;
        assert_eq!(replay_states.len(), 6);
        assert_eq!(
            replay_states
                .iter()
                .filter(|state| !state.local_edges.is_empty())
                .count(),
            5
        );
        Ok(())
    }

    #[test]
    fn union_terms_project_factorized_typed_4d_values() -> Result<()> {
        test_initialise().unwrap();
        let graph: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v2;
                v1 -> v1;v1 -> v1;
            },"scalars"
        )?;
        // Disable integrated terms so each union equality isolates factorized
        // local 4D composition.
        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: false,
                local_uv_cts_from_expanded_4d_integrands: true,
                ..Default::default()
            },
            ..Default::default()
        };
        let cut_structure = CutStructure::empty(&graph);
        let mut forests = Wood::new(cut_structure, &graph, &settings)?.unfold();
        forests.integrate(
            &graph,
            crate::utils::vakint()?,
            &settings.uv,
            settings.medium.vacuum_subtraction,
        )?;

        let unions = forests
            .graph
            .iter_nodes()
            .filter_map(|(node, _, operation)| {
                (forests.graph.is_disjoint_union(node) && !operation.key.is_empty()).then_some(node)
            })
            .collect::<Vec<_>>();
        assert!(!unions.is_empty());

        for union in unions {
            let components = forests.disconnected_component_nodes(union)?;
            let expected = Local4dCts::from_full_product(
                components
                    .into_iter()
                    .map(|component| forests.recursion_input_4d(component))
                    .collect::<Result<Vec<_>>>()?,
            );
            let operation = &forests.graph[union];
            let local = forests
                .compute_store
                .require(operation)?
                .local_4d(operation)?;

            // Compare the complete typed union expression; no per-cut parent cache
            // or component-path replay participates.
            assert!((local.atom().collect_factors() - expected.atom().collect_factors()).is_zero());
        }
        Ok(())
    }

    #[test]
    fn union_terms_replay_component_paths_from_typed_roots() -> Result<()> {
        use crate::graph::feynman_graph::FeynmanGraph;

        test_initialise().unwrap();
        let mut graph: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v2;
                v1 -> v1;v1 -> v1;
            },"scalars"
        )?;
        let settings = GenerationSettings::default();
        let cut_structure = CutStructure::empty(&graph);
        let cutset = cut_structure
            .cuts
            .first()
            .expect("empty cut structure has one cut")
            .clone();
        let mut forests = Wood::new(cut_structure, &graph, &settings)?.unfold();

        let root_disconnected = forests
            .graph
            .iter_nodes()
            .find_map(|(node, _, _)| {
                forests
                    .normalized_node_label(node)
                    .starts_with("{36,F}:")
                    .then_some(node)
            })
            .expect("lopsided dumbbell should contain a root disconnected frontier");
        let dependent_disconnected = forests
            .graph
            .iter_nodes()
            .find_map(|(node, _, _)| {
                forests
                    .normalized_node_label(node)
                    .starts_with("{C} · {36,F}:")
                    .then_some(node)
            })
            .expect("lopsided dumbbell should contain a C-dependent disconnected frontier");
        let non_root_frontier = forests
            .local_leaf_operations(dependent_disconnected)
            .first()
            .expect("dependent disconnected operation has a leaf")
            .frontier;
        assert_eq!(
            forests.graph[non_root_frontier]
                .covers()
                .expect("non-root frontier has a cover")
                .string_label(),
            "C"
        );

        let operations = forests
            .graph
            .iter_nodes()
            .map(|(_, _, operation)| operation.clone())
            .collect::<Vec<_>>();
        // Zero integrated terms isolate component-path replay in the local accumulator.
        for operation in operations {
            forests
                .compute_store
                .entry(operation)
                .or_default()
                .integrated = Some(IntegratedCts::root());
        }

        let orientation_pattern = crate::settings::global::OrientationPattern::default();
        let options = graph.denominator_only_cff_3d_expression_options();
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        let contract_edges = graph
            .iter_edges_of(&graph.tree_edges)
            .map(|(_, edge, _)| edge)
            .collect_vec();
        let production = graph.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonization,
            &options,
            Some(&Atom::one()),
        )?;
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact_expression(
                &production,
                &options,
                &orientation_pattern,
                false,
            ),
        );
        let seed = forests.local_3d_for_node(
            forests.root,
            &mut graph,
            &cutset,
            localizer,
            &settings.uv,
            &mut Local4dProjectionContext::default(),
        )?;
        let root_store_marker = symbolica::symbol!("root_store_marker");
        let frontier_store_marker = symbolica::symbol!("frontier_store_marker");
        // Mark cached local terms so the union result proves that replay starts from typed roots.
        let marked_computation = |marker| -> Result<CutComputation> {
            Ok(CutComputation {
                local_3d: seed.local_3d.map(|_| Ok(Atom::var(marker)))?,
                final_integrands: seed.final_integrands.clone(),
            })
        };

        let root_operation = forests.graph[forests.root].clone();
        forests
            .compute_store
            .entry(root_operation)
            .or_default()
            .cuts
            .insert(cutset.clone(), marked_computation(root_store_marker)?);
        forests
            .compute_store
            .entry(forests.graph[non_root_frontier].clone())
            .or_default()
            .cuts
            .insert(cutset.clone(), marked_computation(frontier_store_marker)?);

        let root_result = forests.local_3d_for_node(
            root_disconnected,
            &mut graph,
            &cutset,
            localizer,
            &settings.uv,
            &mut Local4dProjectionContext::default(),
        )?;
        assert!(
            root_result
                .local_3d
                .direct()?
                .branches()?
                .iter_keys()
                .all(|(_, integrands)| integrands
                    .iter()
                    .all(|(_, term)| !term.contains_symbol(root_store_marker))),
            "an empty dependency frontier must start from the typed root"
        );

        let frontier_result = forests.local_3d_for_node(
            dependent_disconnected,
            &mut graph,
            &cutset,
            localizer,
            &settings.uv,
            &mut Local4dProjectionContext::default(),
        )?;
        let replay_states = forests.union_replay_states(&graph, dependent_disconnected)?;
        let state = replay_states
            .iter()
            .find(|state| {
                state.local_edges.len() == 2
                    && forests
                        .normalized_node_label(state.integrated)
                        .starts_with("{36}:")
            })
            .expect("the integrated 36 prefix has the two-step C-to-F suffix");
        let expected_active = state
            .local_edges
            .iter()
            .enumerate()
            .map(|(offset, edge)| {
                let step_order = forests.graph[state.integrated].key.op_count() + offset;
                let (current, given) = forests.wood.current_given_pair(*edge, step_order);
                current.reduced_subgraph(&given)
            })
            .reduce(|active, reduced| active.union(&reduced))
            .expect("the selected replay state has a local suffix");
        assert!(
            expected_active.empty_intersection(forests.source_spinney(state.integrated).filter())
        );
        let sector = frontier_result
            .local_3d
            .direct()?
            .sectors()?
            .iter()
            .find(|sector| sector.active_subgraph == expected_active)
            .expect("the replay frontier contains the expected active sector");
        assert!(
            !sector.projection_paths.is_empty()
                && sector.projection_paths.iter().all(|path| {
                    path.steps.len() == state.local_edges.len()
                        && path
                            .steps
                            .windows(2)
                            .all(|pair| pair[0].topo_order <= pair[1].topo_order)
                }),
            "each disconnected sector must retain its own ordered local replay path"
        );
        assert!(
            frontier_result
                .local_3d
                .direct()?
                .branches()?
                .iter_keys()
                .all(|(_, integrands)| integrands.iter().all(|(_, term)| {
                    !term.contains_symbol(root_store_marker)
                        && !term.contains_symbol(frontier_store_marker)
                })),
            "a union must replay every component path from its typed root"
        );
        Ok(())
    }

    #[test]
    fn disconnected_soft_union_is_the_product_of_completed_component_atoms() -> Result<()> {
        test_initialise().unwrap();
        let mut graph: Graph =
            include_str!("../../../../tests/resources/graphs/local_ir_dod2_scalar_spectacles.dot")
                .into_graph(&crate::utils::load_generic_model("scalars"))?;
        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: false,
                renormalization_prescription: RenormalizationPrescriptionSettings {
                    log_divergent: ApproximationType::IR,
                    massive_power_divergent: ApproximationType::IR,
                    massless_power_divergent: ApproximationType::IR,
                    ..Default::default()
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let mut forests = Wood::new(CutStructure::empty(&graph), &graph, &settings)?.unfold();
        forests.integrate(
            &graph,
            crate::utils::vakint()?,
            &settings.uv,
            settings.medium.vacuum_subtraction,
        )?;

        let union = forests
            .graph
            .iter_nodes()
            .find_map(|(node, _, _)| forests.graph.is_disjoint_union(node).then_some(node))
            .expect("the spectacles fixture must contain a disconnected UV union");
        assert_eq!(
            forests.source_spinney(union).renormalization_scheme,
            ApproximationType::MUV,
            "a disconnected union carries only the neutral aggregate scheme"
        );
        let components = forests.disconnected_component_nodes(union)?;
        assert_eq!(components.len(), 2);
        assert!(components.iter().all(|component| {
            forests.source_spinney(*component).renormalization_scheme == ApproximationType::IR
        }));

        let expected =
            components
                .iter()
                .try_fold(Atom::one(), |product, component| -> Result<Atom> {
                    let completed = forests.recursion_input_4d(*component)?;
                    assert!(
                        !completed.atom().is_zero(),
                        "each disconnected soft component must be completed before the join"
                    );
                    Ok(product * completed.atom())
                })?;
        let operation = &forests.graph[union];
        let actual = forests
            .compute_store
            .require(operation)?
            .local_4d(operation)?;
        assert_eq!(actual.atom(), &expected);
        let provenance = forests.node_export_4d_provenance(union, actual);
        assert_eq!(provenance.node.parent_keys.len(), 2);
        let UVForestLocalProvenance::FourD {
            components: provenance_components,
            branches,
        } = &provenance.node.local
        else {
            panic!("a computed 4D union must carry 4D provenance")
        };
        assert_eq!(provenance_components.len(), 2);
        assert!(provenance_components.iter().all(|component| {
            component.scheme == ApproximationType::IR
                && !component.route_loop_edges.is_empty()
                && !component.route_signatures.is_empty()
        }));
        assert_eq!(branches.len(), 1);
        assert_eq!(branches[0].branch, Local4dBranch::Combined);
        assert!(branches[0].byte_size > 0);

        // Exercise the same factorization directly in one parent-CFF
        // orientation. With integrated generation disabled, only the replay
        // sector containing both local component operations may survive at
        // the disconnected-union node.
        let cutset = CutSet::empty(graph.n_hedges());
        let orientation_pattern = crate::settings::global::OrientationPattern::default();
        let options = graph.production_3d_expression_options(
            &Default::default(),
            three_dimensional_reps::RepresentationMode::Cff,
        )?;
        let canonization = graph.get_esurface_canonization(&graph.loop_momentum_basis);
        // Match the production root: its structural bridge remains a 4D
        // denominator outside the CFF and its loop-energy contour.
        let contract_edges = graph
            .iter_edges_of(
                &graph
                    .tree_edges
                    .subtract(&graph.initial_state_cut)
                    .subtract(&graph.external_filter::<SuBitGraph>()),
            )
            .map(|(_, edge, _)| edge)
            .collect::<Vec<_>>();
        let mut production = graph.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonization,
            &options,
            None,
        )?;
        assert!(
            !production.expression.orientations.is_empty(),
            "the spectacles parent CFF has at least one orientation"
        );
        production.expression.orientations.truncate(1);
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact_expression(
                &production,
                &options,
                &orientation_pattern,
                false,
            ),
        );
        let union_3d = forests.local_3d_for_node(
            union,
            &mut graph,
            &cutset,
            localizer,
            &settings.uv,
            &mut Local4dProjectionContext::default(),
        )?;
        let sectors = union_3d
            .local_3d
            .direct()?
            .sectors()
            .expect("a disconnected local 3D union retains its replay sectors");
        assert_eq!(sectors.len(), 3);
        let combined_sectors = sectors
            .iter()
            .map(DirectSector::combine)
            .collect::<Result<Vec<_>>>()?;
        let nonzero_sectors = combined_sectors
            .iter()
            .filter(|branches| {
                branches
                    .iter_keys()
                    .any(|(_, integrands)| integrands.iter().any(|(_, atom)| !atom.is_zero()))
            })
            .collect::<Vec<_>>();
        assert_eq!(
            nonzero_sectors.len(),
            1,
            "zero integrated prefixes must leave only the two-local-component cross sector"
        );
        let cross = nonzero_sectors[0];

        let assert_same_branches = |left: &DirectResidueBranches,
                                    right: &DirectResidueBranches,
                                    label: &str|
         -> Result<()> {
            let left = left.iter_keys().collect::<Vec<_>>();
            let right = right.iter_keys().collect::<Vec<_>>();
            assert_eq!(
                left.len(),
                right.len(),
                "{label} differs in residue-map count"
            );
            for ((left_key, left), (right_key, right)) in left.into_iter().zip(right) {
                assert_eq!(
                    left_key, right_key,
                    "{label} changed the complete residue map"
                );
                // Independently replayed forests have different numerator-family
                // scopes. Compare their values using the retained definitions.
                let left = left.resolved()?;
                let right = right.resolved()?;
                left.checked_zip(&right, |residue, left, right| {
                    // Compare the exact factorized forest values before selectors
                    // are materialized or any graph numerator is distributed.
                    let difference = (left - right).together().cancel();
                    assert!(
                        difference.is_zero(),
                        "{label} differs in residue {residue:?}: {difference}"
                    );
                    Ok(Atom::Zero)
                })?;
            }
            Ok(())
        };
        assert_same_branches(
            &union_3d.local_3d.direct()?.branches()?,
            cross,
            "the disconnected aggregate and its cross sector",
        )?;

        let component_edges = components
            .iter()
            .map(|component| {
                let (parent, edge) = forests
                    .graph
                    .unique_parent(*component)
                    .expect("a primitive disconnected component has one root edge");
                assert!(
                    forests.graph[parent].key.is_empty(),
                    "the spectacles components must both be primitive"
                );
                edge
            })
            .collect::<Vec<_>>();
        let mut replay = |edges: &[EdgeIndex]| -> Result<Direct3dCts> {
            let mut local = Direct3dCts::root(&graph, localizer)?;
            for (order, edge) in edges.iter().enumerate() {
                let (current, given) = forests.wood.current_given_pair(*edge, order);
                local = Direct3dApproximation::new(localizer, &mut graph, &settings.uv)
                    .run_local(&local, &current, &given, &current, &given)?;
            }
            assert!(local.branches()?.iter_keys().all(|(_, integrands)| {
                integrands
                    .iter()
                    .all(|(_, atom)| !atom.contains_symbol(GS.rescale))
            }));
            Ok(local)
        };
        let forward = replay(&component_edges)?;
        let reverse = replay(&component_edges.iter().rev().copied().collect::<Vec<_>>())?;
        assert_same_branches(
            &forward.branches()?,
            cross,
            "forward component-local replay and the produced cross sector",
        )?;
        assert_same_branches(
            &reverse.branches()?,
            cross,
            "reverse component-local replay and the produced cross sector",
        )?;
        Ok(())
    }

    #[test]
    fn nested_soft_top_self_energy_preserves_completed_parent_poles() -> Result<()> {
        initialisation::test_initialise()?;
        let model = utils::load_generic_model("sm");
        let graph: Graph = include_str!(
            "../../../../tests/resources/graphs/local_nested_soft_top_self_energy.dot"
        )
        .into_graph(&model)?;
        let external: SuBitGraph = graph.external_filter();
        let parent = graph.full_filter().subtract(&external);
        assert_eq!(graph.n_loops(&parent), 2);
        let physical_mass = graph[EdgeIndex::from(2)].mass_atom();
        assert!(!physical_mass.is_zero());
        assert_ne!(physical_mass, Atom::var(GS.m_uv_vacuum));

        let mut completed_poles = Vec::with_capacity(2);
        for scheme in [ApproximationType::MUV, ApproximationType::IR] {
            let settings = GenerationSettings {
                uv: UVgenerationSettings {
                    generate_integrated: true,
                    renormalization_prescription: RenormalizationPrescriptionSettings {
                        log_divergent: scheme,
                        massive_power_divergent: scheme,
                        massless_power_divergent: scheme,
                        ..Default::default()
                    },
                    ..Default::default()
                },
                ..Default::default()
            };
            let mut forests = Wood::new(CutStructure::empty(&graph), &graph, &settings)?.unfold();
            forests.integrate(
                &graph,
                utils::vakint()?,
                &settings.uv,
                settings.medium.vacuum_subtraction,
            )?;
            let mut pole = Atom::Zero;
            let mut bare_parent_count = 0;
            let mut nested_parent_count = 0;
            let mut finite_prefix_retained = false;
            for (_, _, operation) in forests.graph.iter_nodes() {
                if operation.covers().as_ref() != Some(&parent) {
                    continue;
                }
                let computed = forests.compute_store.require(operation)?;
                pole += computed.integrated(operation)?.physical_pole_atom();
                if operation.key.op_count() == 1 {
                    bare_parent_count += 1;
                } else {
                    nested_parent_count += 1;
                    finite_prefix_retained |=
                        computed.local_4d(operation)?.active_sectors().len() > 1;
                }
            }
            assert_eq!(bare_parent_count, 1);
            assert!(nested_parent_count > 0);
            assert!(
                finite_prefix_retained,
                "the nested {scheme} parent must include the integrated child's finite completion"
            );
            // Compare the complete parent forest, keeping its Feynman-rule numerator and
            // finite-prefix epsilon jets. Resolve the model's color representation
            // invariants before comparing; individual nested terms need not agree.
            let pole = pole
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum)
                .replace(GS.dim)
                .with(4)
                .simplify_metrics()
                .to_dots()
                .simplify_color_with(
                    ColorSimplifySettings::default().with_cof_dimension_invariants(),
                )
                .canonize(Aind::Dummy)?
                .collect_factors();
            assert!(
                !pole.is_zero(),
                "the completed {scheme} parent pole must be nonzero"
            );
            completed_poles.push(pole);
        }
        let difference = (&completed_poles[0] - &completed_poles[1])
            .together()
            .cancel();
        assert!(
            difference.is_zero(),
            "MUV and H must give the same completed parent poles: {difference}"
        );
        Ok(())
    }

    #[test]
    fn heterogeneous_terminal_projects_components_before_multiplying() -> Result<()> {
        test_initialise().unwrap();
        let graph: Graph = dot!(
            digraph G {
                num = "1";
                projector = "1";
                overall_factor = "1";
                edge [particle = "scalar_1"];
                node [num = "1"];
                v1 -> v2;
                v1 -> v1;
                v2 -> v3;
                v2 -> v3;
            },
            "scalars"
        )?;
        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                softct: false,
                renormalization_prescription: RenormalizationPrescriptionSettings {
                    log_divergent: ApproximationType::MUV,
                    massive_power_divergent: ApproximationType::PolePart,
                    massless_power_divergent: ApproximationType::PolePart,
                    ..Default::default()
                },
                ..Default::default()
            },
            ..Default::default()
        };
        let mut forests = Wood::new(CutStructure::empty(&graph), &graph, &settings)?.unfold();
        forests.integrate(
            &graph,
            crate::utils::vakint()?,
            &settings.uv,
            settings.medium.vacuum_subtraction,
        )?;

        let terminals = forests
            .graph
            .iter_nodes()
            .filter_map(|(node, mut crown, _)| {
                (!crown.any(|hedge| forests.graph.flow(hedge).is_source())).then_some(node)
            })
            .collect::<Vec<_>>();
        assert_eq!(terminals.len(), 1);
        let terminal = terminals[0];
        assert!(forests.graph.is_disjoint_union(terminal));

        let components = forests.disconnected_component_nodes(terminal)?;
        assert_eq!(components.len(), 2);
        let schemes = components
            .iter()
            .map(|component| forests.source_spinney(*component).renormalization_scheme)
            .collect::<Vec<_>>();
        assert!(schemes.contains(&ApproximationType::MUV));
        assert!(schemes.contains(&ApproximationType::PolePart));

        let expected =
            components
                .iter()
                .try_fold(Atom::one(), |product, &component| -> Result<Atom> {
                    let component_key = &forests.graph[component];
                    let integrated = forests
                        .compute_store
                        .require(component_key)?
                        .integrated(component_key)?;
                    let projection = match forests.source_spinney(component).renormalization_scheme
                    {
                        ApproximationType::MUV => integrated.physical_finite_counterterm_atom(),
                        ApproximationType::PolePart => integrated.physical_pole_atom(),
                        scheme => panic!("unexpected component scheme {scheme}"),
                    };
                    Ok(product * projection)
                })?;
        let terminal_key = &forests.graph[terminal];
        let aggregate = forests
            .compute_store
            .require(terminal_key)?
            .integrated(terminal_key)?;
        assert_ne!(
            expected.expand(),
            aggregate.physical_finite_counterterm_atom().expand()
        );
        assert_ne!(expected.expand(), aggregate.physical_pole_atom().expand());

        let wild = Atom::var(W_.x___);
        let replacements =
            graph.integrand_replacement(&graph.full_filter(), &graph.loop_momentum_basis, &[wild]);
        let expected = expected
            .simplify_color()
            .expand_num()
            .to_dots()
            .replace_multiple(&replacements)
            .replace(GS.m_uv_expansion)
            .with(GS.m_uv_vacuum);
        let actual = forests
            .renormalization_part_of_ends(&graph, &settings.uv)?
            .expression;
        assert_eq!(actual.expand(), expected.expand());

        Ok(())
    }

    #[test]
    fn triple_tadpole() -> Result<()> {
        test_initialise().unwrap();
        let dumbell: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v3;
                v3 -> v3;
                v2 -> v2;
                v1 -> v1;
            },"scalars"
        )?;

        let spinneys: Vec<_> = dumbell
            .spinneys(&dumbell.full_filter())
            .into_iter()
            .filter_map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis).transpose())
            .collect::<Result<_>>()?;
        let f = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &GenerationSettings::default(),
        )?;

        println!("{}", f);

        insta::assert_snapshot!(
            f.graph.n_nodes(),
            @"8",
            // format!("Wood does not have correct number of spinneys: \n{}",f)
        );

        let three_component_target = f
            .graph
            .iter_nodes()
            .find_map(|(node, _, spinney)| (spinney.n_components() == 3).then_some(node))
            .expect("triple tadpole should contain its fully factorized target");
        let target_filter = f.graph[three_component_target].filter().clone();
        let components = dumbell
            .as_ref()
            .connected_components(&target_filter)
            .into_iter()
            .collect::<BTreeSet<_>>();
        assert_eq!(components.len(), 3);
        let incoming = f
            .graph
            .iter_crown(three_component_target)
            .filter(|hedge| f.graph.flow(*hedge) == Flow::Sink)
            .map(|hedge| f.graph[f.graph[&hedge]].clone())
            .collect::<BTreeSet<_>>();
        assert_eq!(incoming, components);
        assert!(
            !incoming.contains(&target_filter),
            "a fully supplied product must not retain a duplicate atomic full-target transition"
        );

        for (_, _, d) in f.graph.iter_nodes() {
            println!(
                "//Node {}: \n{}",
                d.subgraph.string_label(),
                dumbell.dot(&d.subgraph)
            );
        }
        let _ff = OldWood::from_spinneys(spinneys, &dumbell); //.unfold(&g, &g.loop_momentum_basis);

        // println!("{}", ff.dot(&dumbell));

        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
            f.graph.n_nodes(),
            @"8");

        let three_component_union = f
            .graph
            .iter_nodes()
            .find_map(|(node, _, operation)| {
                (operation.key.op_count() == 3
                    && f.wood
                        .join_factors(f.graph.source_node(node))
                        .is_some_and(|factors| factors.len() == 3))
                .then_some(node)
            })
            .expect("triple tadpole should contain a three-component union");
        assert_eq!(
            f.disconnected_component_nodes(three_component_union)?.len(),
            3
        );
        let replay_states = f.union_replay_states(&dumbell, three_component_union)?;
        assert_eq!(replay_states.len(), 8);
        assert_eq!(
            replay_states
                .iter()
                .filter(|state| !state.local_edges.is_empty())
                .count(),
            7
        );

        Ok(())
    }

    #[test]
    fn saclay() -> Result<()> {
        test_initialise().unwrap();
        let dt: Graph = dot!(digraph GL16{

        num = "spenso::g(spenso::coad(8,gammalooprs::hedge(8)),spenso::coad(8,gammalooprs::hedge(11)))"
                ext	 [style=invis];
                ext	-> 0  [dir=none id=0 particle="g"];
                2	-> ext  [dir=none id=1 particle="g"];
                0	-> 1 -> 2->3->0  [ particle="d"];
                      1 ->3 [particle = "g"]
                    })?;

        let f = Wood::new(
            CutStructure::empty(&dt),
            &dt,
            &GenerationSettings::default(),
        )?;

        println!("{}", dt.dot_serialize(&DotExportSettings::default()));

        insta::assert_snapshot!(
            f.graph.n_nodes(),
            @"5",
            // format!("Wood does not have correct number of spinneys: \n{}",f)
        );

        let f = f.unfold();
        let structured_dot = f.dot_serialize();
        assert!(structured_dot.contains(r#"label="K[#expr S_44⊛0]""#));
        assert!(!structured_dot.contains("#T("));
        println!("{structured_dot}");
        insta::assert_snapshot!(
            f.graph.n_nodes(),
            @"8");

        Ok(())
    }

    #[test]
    fn dumbells() -> Result<()> {
        test_initialise().unwrap();
        let dumbell: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v2;
                v1 -> v1;
            },"scalars"
        )?;

        let spinneys: Vec<_> = dumbell
            .spinneys(&dumbell.full_filter())
            .into_iter()
            .filter_map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis).transpose())
            .collect::<Result<_>>()?;
        let f = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &GenerationSettings::default(),
        )?;

        println!("{}", f);

        insta::assert_snapshot!(
            f.graph.n_nodes(),
            @"4",
            // format!("Wood does not have correct number of spinneys: \n{}",f)
        );

        for (_, _, d) in f.graph.iter_nodes() {
            println!(
                "//Node {}: \n{}",
                d.subgraph.string_label(),
                dumbell.dot(&d.subgraph)
            );
        }
        let _ff = OldWood::from_spinneys(spinneys, &dumbell); //.unfold(&g, &g.loop_momentum_basis);

        // println!("{}", ff.dot(&dumbell));

        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
            f.graph.n_nodes(),
            @"4");

        Ok(())
    }

    #[test]
    fn bugblatter() -> Result<()> {
        test_initialise().unwrap();

        match dot!(
            digraph G{
                A1 -> A2 [particle="t"];
                A2 -> A3 [particle="t"];
                A3 -> A1 [particle="t"];
                B1 -> B2 [particle="t"];
                B2 -> B3 [particle="t"];
                B3 -> B1 [particle="t"];
                A1 -> B1 [particle="a"];
                A2 -> B2 [particle="a"];
                A3 -> B3 [particle="a"];
            },"sm"
        ) {
            Ok(g) => {
                let g: Graph = g;
                let spinneys: Vec<_> = g
                    .spinneys(&g.full_filter())
                    .into_iter()
                    .filter_map(|a| Spinney::new(a, &g, &g.loop_momentum_basis).transpose())
                    .collect::<Result<_>>()?;
                let f = Wood::new(CutStructure::empty(&g), &g, &GenerationSettings::default())?;

                println!("{}", f);

                assert_eq!(
                    20,
                    f.graph.n_nodes(),
                    "Wood does not have correct number of spinneys: \n{}",
                    f
                );

                for (_, _, d) in f.graph.iter_nodes() {
                    println!(
                        "//Node {}: \n{}",
                        d.subgraph.string_label(),
                        g.dot_lmb_of(&d.subgraph, &d.lmb)
                    );
                }
                let _ff = OldWood::from_spinneys(spinneys, &g); //.unfold(&g, &g.loop_momentum_basis);

                // println!("{}", ff.dot(&g));

                let f = f.unfold();
                f.debug_walk();
                println!("{}", f);
                assert_eq!(
                    152,
                    f.graph.n_nodes(),
                    "Forest unfolds into the wrong number of terms :\n{}",
                    f
                );

                // println!("{}", f)
            }
            Err(e) => {
                eprintln!("{}", e);
            }
        }

        // let f = SpinneyWood::from_spinneys(
        //     g.spinneys(&g.full_filter()).into_iter().map(|a| a.filter),
        //     &g,
        // )
        // .unfold();

        // println!("{}", f.graph.base_dot());
        Ok(())
    }

    #[test]
    fn mercedes() -> Result<()> {
        test_initialise().unwrap();

        let mercedes: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v3;
                v3 -> v1;
                v1 -> v4;
                v2 -> v4;
                v3 -> v4;
            },"scalars"
        )?;

        let f = Wood::new(
            CutStructure::empty(&mercedes),
            &mercedes,
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"2",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"2",
         );

        Ok(())
    }

    #[test]
    fn sunrise() -> Result<()> {
        test_initialise().unwrap();

        let sunrise: Graph = dot!( digraph sunrise{
            node [num = "1"]
            edge [particle=scalar_1]
            e        [style=invis]
            e -> A:0   [ id=3 ]
            B:1 -> e   [ id=4 ]
            A -> B    [ id=0 ]
            A -> B    [ id=1 ]
            A -> B    [ id=2 ]
        },"scalars")?;
        // let spinneys = spectacles.spinneys(&spectacles.full_filter());
        let f = Wood::new(
            CutStructure::empty(&sunrise),
            &sunrise,
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"5",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"8",
         );

        Ok(())
    }
    #[test]
    fn dotted_sunrise() -> Result<()> {
        test_initialise().unwrap();

        let sunrise: Graph = dot!( digraph sunrise{
            edge [particle=scalar_1]
            e        [style=invis]
            e -> A:0   [ id=3]
            B:1 -> e   [ id=4]

            A -> C    [ id=0]
            C -> e
            C -> B
            A -> B    [ id=1]
            A -> B    [ id=2]
        },"scalars")?;
        // let spinneys = spectacles.spinneys(&spectacles.full_filter());
        let f = Wood::new(
            CutStructure::empty(&sunrise),
            &sunrise,
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"3",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"4",
         );

        Ok(())
    }

    #[test]
    fn dotted() -> Result<()> {
        test_initialise().unwrap();

        let sunrise: Graph = dot!( digraph sunrise{
            edge [particle=scalar_1]
            e        [style=invis]
            e -> A:0   [ id=4]
            B:1 -> e   [ id=5]
            C:2 -> e   [ id=6]

            A -> B    [ id=0]
            B -> C     [ id=1]
            C -> A   [ id=2]
            B -> C    [ id=3]

        },"scalars")?;
        // let spinneys = spectacles.spinneys(&spectacles.full_filter());
        let f = Wood::new(
            CutStructure::empty(&sunrise),
            &sunrise,
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"3",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"4",
         );

        Ok(())
    }

    #[test]
    fn spectacles() -> Result<()> {
        use crate::graph::feynman_graph::FeynmanGraph;

        test_initialise().unwrap();

        let mut spectacles: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v1 -> v2;

                v3 -> v4;
                v3 -> v4;

                v2 -> v3;
                v1 -> v4;
            },"scalars"
        )?;

        // let spinneys = spectacles.spinneys(&spectacles.full_filter());
        let settings = GenerationSettings::default();
        let cut_structure = CutStructure::empty(&spectacles);
        let cutset = cut_structure
            .cuts
            .first()
            .expect("empty cut structure has one cut")
            .clone();
        let f = Wood::new(cut_structure, &spectacles, &settings)?;
        println!("{}", f);
        let mut f = f.unfold();
        println!("{}", f);

        let (union, edge) = f
            .graph
            .iter_nodes()
            .find_map(|(child, _, _)| {
                let (parent, edge) = f.graph.unique_parent(child)?;
                (!f.graph.is_disjoint_union(child) && f.graph.is_disjoint_union(parent))
                    .then_some((parent, edge))
            })
            .expect("spectacles has a connected child above its disconnected union");
        let operations = f
            .graph
            .iter_nodes()
            .map(|(_, _, operation)| operation.clone())
            .collect::<Vec<_>>();
        // Zero integrated terms isolate component-path replay in the local accumulator.
        for operation in operations {
            f.compute_store.entry(operation).or_default().integrated = Some(IntegratedCts::root());
        }
        let orientation_pattern = crate::settings::global::OrientationPattern::default();
        // The replay fixture keeps its synthetic unit branch, while the
        // existing exact source supplies a valid energy map for that key.
        let options = spectacles.denominator_only_cff_3d_expression_options();
        let canonization = spectacles.get_esurface_canonization(&spectacles.loop_momentum_basis);
        let contract_edges = spectacles.paired_edges(
            &spectacles
                .tree_edges
                .subtract(&spectacles.initial_state_cut),
        );
        let production = spectacles.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonization,
            &options,
            None,
        )?;
        assert!(!production.expression.orientations.is_empty());
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::exact_expression(
                &production,
                &options,
                &orientation_pattern,
                false,
            ),
        );
        let replay_states = f
            .union_replay_states(&spectacles, union)?
            .into_iter()
            .filter(|state| !state.local_edges.is_empty())
            .collect::<Vec<_>>();
        let union_active = replay_states
            .iter()
            .map(|state| {
                state
                    .local_edges
                    .iter()
                    .enumerate()
                    .map(|(offset, edge)| {
                        let step_order = f.graph[state.integrated].key.op_count() + offset;
                        let (current, given) = f.wood.current_given_pair(*edge, step_order);
                        current.reduced_subgraph(&given)
                    })
                    .reduce(|active, reduced| active.union(&reduced))
                    .expect("a proper union replay state has a local suffix")
            })
            .collect::<Vec<_>>();
        let computed = f.local_3d_for_node(
            union,
            &mut spectacles,
            &cutset,
            localizer,
            &settings.uv,
            &mut Local4dProjectionContext::default(),
        )?;
        let sectors = computed.local_3d.direct()?.sectors()?;
        assert_eq!(
            sectors.len(),
            replay_states.len(),
            "every proper integrated prefix must produce one disconnected replay sector"
        );
        assert_eq!(
            sectors
                .iter()
                .map(|sector| sector.active_subgraph.clone())
                .collect::<Vec<_>>(),
            union_active,
            "the completed union must retain the active component support of every replay path"
        );
        assert!(
            computed.final_integrands.iter().next().is_some(),
            "the successfully replayed disconnected union must reach final-integrand construction"
        );
        let union_local = Direct3dCts::from_sectors(
            union_active
                .iter()
                .cloned()
                .map(|active_subgraph| {
                    let coordinate_frames = spectacles
                        .as_ref()
                        .connected_components(&active_subgraph)
                        .into_iter()
                        .filter(|component| spectacles.n_loops(component) != 0)
                        .map(|component| DirectCoordinateFrame {
                            lmb: spectacles.lmb_of(&component),
                            active_subgraph: component,
                        })
                        .collect();
                    Ok(DirectSector {
                        active_subgraph,
                        coordinate_frames,
                        active: DirectResidueBranches::production(
                            OrientationID(0),
                            Integrands::root(),
                        )?,
                        frozen_integrands: Integrands::root(),
                        projection_paths: Vec::new(),
                    })
                })
                .collect::<Result<_>>()?,
        )?;

        let step_order = f.graph[union].key.op_count();
        let (current, given) = f.wood.current_given_pair(edge, step_order);
        let reduced = current.reduced_subgraph(&given);
        let numerator = spectacles
            .numerator(&reduced, given.subgraph())
            .get_single_atom()
            .expect("the scalar reduced numerator is available");
        assert!(!numerator.contains_symbol(GS.emr_mom));
        // A Taylor projection leaves these momentum-independent source values
        // unchanged and supplies the subtraction sign. Compare the complete
        // result so retained and quotient coordinates can use any compatible
        // basis without fixing the number or storage of coordinate frames.
        let expected = union_local
            .branches()?
            .materialize(false)?
            .resolved()?
            .map(|atom| -(atom * &numerator));
        let given_internal =
            InternalSubGraph::try_new(given.subgraph().clone(), spectacles.as_ref())
                .expect("the connected child has an internal given prefix");
        let quotient_rank = spectacles
            .shrunken_sub_lmb(
                current.subgraph(),
                &given_internal,
                spectacles.dummy_stripped_external_flows_of(current.subgraph()),
                None,
            )?
            .loop_edges
            .len();
        let expected_frames = union_local
            .sectors()?
            .iter()
            .map(|sector| {
                let retained_region = sector.active_subgraph.intersection(given.subgraph());
                let mut retained = sector
                    .coordinate_frames
                    .iter()
                    .filter(|frame| retained_region.includes(&frame.active_subgraph))
                    .flat_map(|frame| frame.lmb.loop_edges.iter().copied())
                    .collect::<Vec<_>>();
                retained.sort();
                retained.dedup();
                let expected_rank = retained.len() + quotient_rank;
                (retained, expected_rank)
            })
            .collect::<Vec<_>>();
        let external: SuBitGraph = spectacles.external_filter();
        let full_internal = spectacles.full_filter().subtract(&external);
        let canonical_route = if current.subgraph() == &full_internal
            && spectacles.n_loops(&full_internal) == spectacles.loop_momentum_basis.loop_edges.len()
        {
            spectacles.loop_momentum_basis.clone()
        } else {
            spectacles.try_compatible_sub_lmb(
                current.subgraph(),
                spectacles.dummy_less_full_crown(current.subgraph()),
                &spectacles.loop_momentum_basis,
            )?
        };
        let expected_active = union_active
            .iter()
            .map(|active| active.union(&reduced))
            .collect::<Vec<_>>();
        let child_local = Direct3dApproximation::new(localizer, &mut spectacles, &settings.uv)
            .run(
                &union_local,
                &IntegratedCts::root(),
                &current,
                &given,
                &current,
                &given,
            )?;
        let actual = child_local.branches()?.materialize(false)?.resolved()?;
        assert_eq!(
            actual.map(|atom| atom.collect_factors()),
            expected.map(|atom| atom.collect_factors()),
            "the connected child must retain every complete source value across its enclosing Taylor operation",
        );
        let child_sectors = child_local.sectors()?;
        let child_active = child_sectors
            .iter()
            .map(|sector| sector.active_subgraph.clone())
            .collect::<Vec<_>>();
        assert_eq!(child_active.len(), 3);
        assert_eq!(child_active, expected_active);
        for (sector, (retained, expected_rank)) in child_sectors.iter().zip(expected_frames.iter())
        {
            assert_eq!(
                sector.coordinate_frames.len(),
                1,
                "one connected enclosing Taylor operation must merge every contained child frame"
            );
            let child_carriers = sector
                .coordinate_frames
                .iter()
                .flat_map(|frame| frame.lmb.loop_edges.iter())
                .copied()
                .collect::<BTreeSet<_>>();
            assert!(
                retained.iter().all(|edge| child_carriers.contains(edge)),
                "the connected child must retain every carrier from the contained union frames"
            );
            assert!(
                child_carriers.iter().all(|edge| {
                    spectacles.loop_momentum_basis.edge_signatures[*edge]
                        .external
                        .iter()
                        .all(|sign| sign.is_zero())
                }),
                "the connected child route must be homogeneous in graph-external momenta"
            );
            assert_eq!(
                sector
                    .coordinate_frames
                    .iter()
                    .map(|frame| frame.lmb.loop_edges.len())
                    .sum::<usize>(),
                *expected_rank,
                "the connected child rank is the retained prefix rank plus its quotient rank"
            );
            let canonical_carriers = canonical_route
                .loop_edges
                .iter()
                .copied()
                .collect::<BTreeSet<_>>();
            let canonical_is_admissible = retained
                .iter()
                .all(|edge| canonical_carriers.contains(edge))
                && canonical_carriers.len() == *expected_rank
                && canonical_carriers.iter().all(|edge| {
                    spectacles.loop_momentum_basis.edge_signatures[*edge]
                        .external
                        .iter()
                        .all(|sign| sign.is_zero())
                });
            if canonical_is_admissible {
                assert_eq!(
                    child_carriers, canonical_carriers,
                    "an admissible graph-canonical route must win before equivalent fallbacks"
                );
            }
        }

        Ok(())
    }

    #[test]
    fn collective_two_component_region_requires_independently_divergent_factors() -> Result<()> {
        test_initialise().unwrap();

        // The two-line component has DOD two, while the scalar triangle has DOD
        // minus two. Their disconnected union has logarithmic simultaneous
        // power counting, but the triangle is not a UV region on its own.
        // Componentwise classification therefore retains only the bubble.
        let graph: Graph = dot!(
            digraph G {
                node [num = "1"];
                edge [particle = "scalar_1"];
                a0 -> a1 [id = 0, dod = "-1", num = "gammalooprs::Q(0,spenso::mink(4,0))"];
                a0 -> a1 [id = 1, dod = "-1", num = "gammalooprs::Q(1,spenso::mink(4,0))"];
                a1 -> b0 [id = 2];
                b0 -> b1 [id = 3];
                b1 -> b2 [id = 4];
                b2 -> b0 [id = 5];
            },
            "scalars"
        )?;

        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: false,
                ..Default::default()
            },
            ..Default::default()
        };
        let wood = Wood::new(CutStructure::empty(&graph), &graph, &settings)?;
        let collective = graph
            .spinneys(&graph.full_filter())
            .into_iter()
            .find(|spinney| graph.as_ref().connected_components(spinney).len() == 2)
            .expect("the raw cycle inventory contains the collective region");
        let collective_filter = collective.filter;
        let components = graph
            .as_ref()
            .connected_components(&collective_filter)
            .into_iter()
            .collect::<BTreeSet<_>>();
        assert_eq!(components.len(), 2);
        assert_eq!(graph.compute_dod(&collective_filter), 0);
        assert_eq!(
            components
                .iter()
                .map(|component| graph.compute_dod(component))
                .sorted()
                .collect::<Vec<_>>(),
            vec![-2, 2]
        );

        let root = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let divergent_filter = components
            .iter()
            .find(|component| graph.compute_dod(*component) >= 0)
            .expect("the bubble is the independently divergent component")
            .clone();
        let divergent = OwnedForestNode {
            spinney: Spinney::new(
                InternalSubGraph::cleaned_filter_optimist(divergent_filter, graph.as_ref()),
                &graph,
                &graph.loop_momentum_basis,
            )?
            .expect("the divergent component has a compatible source chart"),
            topo_order: 0,
        };
        let collective_source = OwnedForestNode {
            spinney: Spinney::new(
                InternalSubGraph::cleaned_filter_optimist(
                    collective_filter.clone(),
                    graph.as_ref(),
                ),
                &graph,
                &graph.loop_momentum_basis,
            )?
            .expect("the collective region has a compatible source chart"),
            topo_order: 1,
        };
        let ctx = UVCtx::new(&graph, &settings.uv);
        let input = Full4dCts::from_factorized_local(&Local4dCts::root());
        let prefix = local_4d::uv_limit(&input, &ctx, &divergent, &root, &divergent, &root)?;
        // The two physical forests containing U are {U} and {A,U}:
        // (-T_U) I + (-T_U)(-T_A) I. Their full Taylor scope is U;
        // its convergent complement does not supply a separate counterterm.
        let atomic = local_4d::uv_limit(
            &input,
            &ctx,
            &collective_source,
            &root,
            &collective_source,
            &root,
        )?;
        let nested = local_4d::uv_limit(
            &Full4dCts::from_factorized_local(&prefix),
            &ctx,
            &collective_source,
            &divergent,
            &collective_source,
            &divergent,
        )?;
        let expected = atomic.atom() + nested.atom();
        let mut forests = wood.unfold();
        forests.integrate(
            &graph,
            crate::utils::vakint()?,
            &settings.uv,
            settings.medium.vacuum_subtraction,
        )?;
        let actual = forests
            .graph
            .iter_nodes()
            .filter(|(_, _, operation)| operation.covers().as_ref() == Some(&collective_filter))
            .try_fold(Atom::Zero, |sum, (_, _, operation)| -> Result<_> {
                Ok(sum
                    + forests
                        .compute_store
                        .require(operation)?
                        .local_4d(operation)?
                        .atom())
            })?;
        // These vacuum components have no nonzero crown shifts. Put both
        // complete coefficients into the same formal K coordinates, retaining
        // factorization and every denominator multiplicity before comparison.
        let replacements = graph.integrand_replacement(
            &graph.full_filter(),
            &graph.loop_momentum_basis,
            &[Atom::var(W_.x___)],
        );
        let [actual, expected] = [actual, expected].map(|expression| {
            GS.erase_uv_momentum_provenance(&expression)
                .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
                .with(W_.prop_)
                .replace_multiple(&replacements)
                .normalize_dots()
                .collect_factors()
        });
        assert!(
            !atomic.atom().is_zero(),
            "the atomic Taylor term must be nonzero"
        );
        assert!((actual - expected).collect_factors().is_zero());

        assert!(
            forests
                .wood
                .graph
                .iter_nodes()
                .all(|(_, _, spinney)| { spinney.filter() != &collective_filter }),
            "a convergent connected factor cannot enter a disconnected UV spinney"
        );
        let convergent = components
            .iter()
            .find(|component| graph.compute_dod(*component) < 0)
            .expect("the collective region has one convergent component");
        let divergent = components
            .iter()
            .find(|component| graph.compute_dod(*component) >= 0)
            .expect("the collective region has one independently divergent component");
        assert!(
            forests
                .wood
                .graph
                .iter_nodes()
                .all(|(_, _, spinney)| spinney.filter() != convergent)
        );
        assert!(
            forests
                .wood
                .graph
                .iter_nodes()
                .any(|(_, _, spinney)| spinney.filter() == divergent)
        );

        assert!(
            forests.graph.iter_nodes().all(|(_, _, operation)| {
                operation.covers().as_ref() != Some(&collective_filter)
            }),
            "neither atomic nor nested collective counterterms may bypass factor selection"
        );
        let mut cut = CutSet::empty(graph.n_hedges());
        cut.union = convergent.clone();
        assert!(
            forests.graph.iter_nodes().any(|(_, _, operation)| {
                operation.covers().as_ref() == Some(divergent) && operation.is_compatible_with(&cut)
            }),
            "a disjoint independently divergent prefix remains cut-compatible"
        );
        cut.union = divergent.clone();
        assert!(
            forests
                .graph
                .iter_nodes()
                .filter(|(_, _, operation)| { operation.covers().as_ref() == Some(divergent) })
                .all(|(_, _, operation)| !operation.is_compatible_with(&cut)),
            "cutting the retained divergent factor must exclude its counterterms"
        );

        Ok(())
    }

    #[test]
    fn collective_three_component_region_keeps_only_its_divergent_join() -> Result<()> {
        test_initialise().unwrap();

        // The two bubbles have DOD two and the scalar triangle has DOD minus two.
        // Every component nevertheless has a divergent complement because the
        // other two components have non-negative combined DOD. The triangle itself
        // supplies no independent counterterm: this aggregate edge-label cover
        // cannot replace the collective Taylor operation by a factorized product.
        // Only the two independently divergent bubbles form a factorized join;
        // their complete Taylor cancellation permits omitting the collective
        // transition rather than classifying the triangle as a UV factor.
        let graph: Graph = dot!(
            digraph G {
                node [num = "1"];
                edge [particle = "scalar_1"];
                a0 -> a1 [id = 0, dod = "-1", num = "gammalooprs::Q(0,spenso::mink(4,0))"];
                a0 -> a1 [id = 1, dod = "-1", num = "gammalooprs::Q(1,spenso::mink(4,0))"];
                a1 -> b0 [id = 2];
                b0 -> b1 [id = 3, dod = "-1", num = "gammalooprs::Q(3,spenso::mink(4,1))"];
                b0 -> b1 [id = 4, dod = "-1", num = "gammalooprs::Q(4,spenso::mink(4,1))"];
                b1 -> c0 [id = 5];
                c0 -> c1 [id = 6];
                c1 -> c2 [id = 7];
                c2 -> c0 [id = 8];
            },
            "scalars"
        )?;

        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: false,
                ..Default::default()
            },
            ..Default::default()
        };
        let wood = Wood::new(CutStructure::empty(&graph), &graph, &settings)?;
        let collective = graph
            .spinneys(&graph.full_filter())
            .into_iter()
            .find(|spinney| graph.as_ref().connected_components(spinney).len() == 3)
            .expect("the raw cycle inventory contains the collective region");
        let collective_filter = collective.filter;
        let components = graph
            .as_ref()
            .connected_components(&collective_filter)
            .into_iter()
            .collect::<BTreeSet<_>>();
        assert_eq!(components.len(), 3);
        assert_eq!(graph.compute_dod(&collective_filter), 2);
        assert_eq!(
            components
                .iter()
                .map(|component| graph.compute_dod(component))
                .sorted()
                .collect::<Vec<_>>(),
            vec![-2, 2, 2]
        );
        assert!(
            wood.graph.iter_nodes().all(|(_, _, spinney)| {
                graph
                    .as_ref()
                    .connected_components(spinney.filter())
                    .iter()
                    .all(|component| graph.compute_dod(component) >= 0)
            }),
            "a convergent connected factor cannot enter a disconnected UV spinney"
        );
        assert!(
            components.iter().all(|component| {
                graph.compute_dod(&collective_filter.subtract(component)) >= 0
            })
        );
        let retained_components = components
            .iter()
            .filter(|component| {
                wood.graph
                    .iter_nodes()
                    .any(|(_, _, spinney)| spinney.filter() == *component)
            })
            .cloned()
            .collect::<BTreeSet<_>>();
        assert_eq!(
            retained_components.len(),
            2,
            "complete complement coverage must not stand in for independent component counterterms"
        );
        assert!(
            retained_components
                .iter()
                .all(|component| graph.compute_dod(component) >= 0)
        );
        let retained_union = retained_components
            .iter()
            .cloned()
            .reduce(|left, right| left.union(&right))
            .unwrap();
        let join = wood
            .graph
            .iter_nodes()
            .find_map(|(node, _, spinney)| (spinney.filter() == &retained_union).then_some(node))
            .expect("the two independently divergent bubbles retain their union");
        assert_eq!(wood.join_factors(join), Some(retained_components));

        let root = OwnedForestNode {
            spinney: Spinney::empty(&graph),
            topo_order: 0,
        };
        let root_local = Local4dCts::root();
        let input = Full4dCts::from_factorized_local(&root_local);
        let ctx = UVCtx::new(&graph, &settings.uv);
        let convergent = components
            .iter()
            .find(|component| graph.compute_dod(*component) < 0)
            .expect("the triangle is the convergent component");
        // Use one common physical chart only because these vacuum components
        // have no nonzero crown shifts. This is not a general replacement for
        // the source reconstruction certificates used by projected CFF.
        let replacements = graph.integrand_replacement(
            &graph.full_filter(),
            &graph.loop_momentum_basis,
            &[Atom::var(W_.x___)],
        );
        let mut independent = Vec::<(OwnedForestNode, Local4dCts)>::new();
        for component in components
            .iter()
            .filter(|component| graph.compute_dod(*component) >= 0)
        {
            let source = OwnedForestNode {
                spinney: Spinney::new(
                    InternalSubGraph::cleaned_filter_optimist(component.clone(), graph.as_ref()),
                    &graph,
                    &graph.loop_momentum_basis,
                )?
                .expect("each divergent bubble has a compatible source chart"),
                topo_order: 0,
            };
            let local = local_4d::uv_limit(&input, &ctx, &source, &root, &source, &root)?;
            let mixed = OwnedForestNode {
                spinney: Spinney::new(
                    InternalSubGraph::cleaned_filter_optimist(
                        component.union(convergent),
                        graph.as_ref(),
                    ),
                    &graph,
                    &graph.loop_momentum_basis,
                )?
                .expect("the bubble-triangle collective region is logarithmic"),
                topo_order: 1,
            };
            let atomic = local_4d::uv_limit(&input, &ctx, &mixed, &root, &mixed, &root)?;
            let nested = local_4d::uv_limit(
                &Full4dCts::from_factorized_local(&local),
                &ctx,
                &mixed,
                &source,
                &mixed,
                &source,
            )?;
            assert!(
                !atomic.atom().is_zero(),
                "the mixed atomic Taylor term must be nonzero"
            );
            // T_AC(1-T_A)I and T_BC(1-T_B)I vanish: the degree-two
            // bubble subtraction leaves no divergent scaling with the triangle.
            // Their complete contributions therefore cancel even if a forest
            // representation includes or removes zero-quotient outer operations.
            let mixed_remainder = GS
                .erase_uv_momentum_provenance(&(atomic.atom() + nested.atom()))
                .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
                .with(W_.prop_)
                .replace_multiple(&replacements)
                .normalize_dots()
                .collect_factors();
            assert!(
                mixed_remainder.is_zero(),
                "mixed collective Taylor terms must cancel"
            );
            independent.push((source, local));
        }
        let joint_scope = independent.iter().fold(
            graph.empty_subgraph::<SuBitGraph>(),
            |scope, (source, _)| scope.union(source.subgraph()),
        );
        let joint = OwnedForestNode {
            spinney: Spinney::new(
                InternalSubGraph::cleaned_filter_optimist(joint_scope, graph.as_ref()),
                &graph,
                &graph.loop_momentum_basis,
            )?
            .expect("the disjoint divergent bubbles have a compatible joint source chart"),
            topo_order: 1,
        };
        let product = Local4dCts::from_full_product(
            independent
                .iter()
                .map(|(_, local)| Full4dCts::from_factorized_local(local)),
        );
        let collective_source = OwnedForestNode {
            spinney: Spinney::new(
                InternalSubGraph::cleaned_filter_optimist(
                    collective_filter.clone(),
                    graph.as_ref(),
                ),
                &graph,
                &graph.loop_momentum_basis,
            )?
            .expect("the collective region has a compatible source chart"),
            topo_order: 2,
        };
        // The surviving complete value is -T_U(1-T_A)(1-T_B)I.
        // AB supplies a product of independent counterterms, while the
        // convergent triangle supplies none; no incoming-edge or join catalogue
        // determines these four physical prefix values.
        let mut expected = Atom::Zero;
        for (given, local) in std::iter::once((&root, &root_local))
            .chain(independent.iter().map(|(source, local)| (source, local)))
            .chain(std::iter::once((&joint, &product)))
        {
            let term = local_4d::uv_limit(
                &Full4dCts::from_factorized_local(local),
                &ctx,
                &collective_source,
                given,
                &collective_source,
                given,
            )?;
            if given.subgraph().is_empty() {
                assert!(
                    !term.atom().is_zero(),
                    "the full atomic Taylor term must be nonzero"
                );
            }
            expected += term.atom();
        }
        let mut forests = wood.unfold();
        forests.integrate(
            &graph,
            crate::utils::vakint()?,
            &settings.uv,
            settings.medium.vacuum_subtraction,
        )?;
        let actual = forests
            .graph
            .iter_nodes()
            .filter(|(_, _, operation)| operation.covers().as_ref() == Some(&collective_filter))
            .try_fold(Atom::Zero, |sum, (_, _, operation)| -> Result<_> {
                Ok(sum
                    + forests
                        .compute_store
                        .require(operation)?
                        .local_4d(operation)?
                        .atom())
            })?;
        let [actual, expected] = [actual, expected].map(|expression| {
            GS.erase_uv_momentum_provenance(&expression)
                .replace(GS.den(W_.a_, W_.mom_, W_.mass_, W_.prop_))
                .with(W_.prop_)
                .replace_multiple(&replacements)
                .normalize_dots()
        });
        assert!(
            (actual - expected).collect_factors().cancel().is_zero(),
            "the complete collective Taylor value must equal -T_U(1-T_A)(1-T_B)I"
        );

        assert!(
            forests.graph.iter_nodes().all(|(_, _, operation)| {
                operation.covers().is_none_or(|scope| {
                    graph
                        .as_ref()
                        .connected_components(&scope)
                        .iter()
                        .all(|component| graph.compute_dod(component) >= 0)
                })
            }),
            "no atomic or nested collective operation may include the convergent triangle"
        );
        assert!(
            forests.graph.iter_nodes().any(|(node, _, operation)| {
                operation.covers().as_ref() == Some(&retained_union)
                    && forests.graph.is_disjoint_union(node)
            }),
            "the complete independently divergent cover remains a factorized join"
        );

        Ok(())
    }

    #[test]
    fn spectacles_typed_4d_local_construction_matches_uv_limit() -> Result<()> {
        test_initialise().unwrap();

        let spectacles: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v1 -> v2;

                v3 -> v4;
                v3 -> v4;

                v2 -> v3;
                v1 -> v4;
            },"scalars"
        )?;

        // let spinneys = spectacles.spinneys(&spectacles.full_filter());
        let settings = GenerationSettings {
            uv: UVgenerationSettings {
                generate_integrated: false,
                local_uv_cts_from_expanded_4d_integrands: true,
                ..Default::default()
            },
            ..Default::default()
        };
        let cut_structure = CutStructure::empty(&spectacles);
        let f = Wood::new(cut_structure, &spectacles, &settings)?;
        println!("{}", f);
        let mut f = f.unfold();
        println!("{}", f);

        let (union, child, edge) = f
            .graph
            .iter_nodes()
            .find_map(|(child, _, _)| {
                let (parent, edge) = f.graph.unique_parent(child)?;
                (!f.graph.is_disjoint_union(child) && f.graph.is_disjoint_union(parent))
                    .then_some((parent, child, edge))
            })
            .expect("spectacles has a connected child above its disconnected union");
        f.integrate(
            &spectacles,
            crate::utils::vakint()?,
            &settings.uv,
            settings.medium.vacuum_subtraction,
        )?;

        let step_order = f.graph[union].key.op_count();
        let (current, given) = f.wood.current_given_pair(edge, step_order);
        let expected = local_4d::uv_limit(
            &f.recursion_input_4d(union)?,
            &UVCtx::new(&spectacles, &settings.uv),
            &current,
            &given,
            &current,
            &given,
        )?;
        let operation = &f.graph[child];
        let local = f.compute_store.require(operation)?.local_4d(operation)?;
        assert!((local.atom().collect_factors() - expected.atom().collect_factors()).is_zero());

        Ok(())
    }

    #[test]
    fn basketball() -> Result<()> {
        test_initialise().unwrap();

        let basketball: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v1 -> v2;
                v1 -> v2;
                v1 -> v2;
            },"scalars"
        )?;

        let f = Wood::new(
            CutStructure::empty(&basketball),
            &basketball,
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"12",
              );
        let f = f.unfold();
        println!("{}", f.dot_serialize());
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"46");
        Ok(())
    }

    #[test]
    fn fourloop_b() -> Result<()> {
        test_initialise().unwrap();

        let fourloop_b: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v1 -> v2;

                v3 -> v4;
                v3 -> v4;

                v2 -> v3;
                v1 -> v3;
                v1 -> v4;
            },"scalars"
        )?;

        // let spinneys = fourloop_b.spinneys(&fourloop_b.full_filter());
        let f = Wood::new(
            CutStructure::empty(&fourloop_b),
            &fourloop_b,
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"14",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"80",
         );

        Ok(())
    }

    #[test]
    fn four_loop_a() -> Result<()> {
        test_initialise().unwrap();

        let four_loop_a: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v1 -> v2;
                v2 -> v3;
                v3 -> v1;
                v1 -> v4;
                v2 -> v4;
                v3 -> v4;
            },"scalars"
        )?;

        let f = Wood::new(
            CutStructure::empty(&four_loop_a),
            &four_loop_a,
            &GenerationSettings::default(),
        )?;
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"12",
        );
        let f = f.unfold();
        println!("{}", f);
        insta::assert_snapshot!(
        f.graph.n_nodes(),
        @"60",
         );

        Ok(())
    }

    #[test]
    fn triple_double_tadpole() -> Result<()> {
        test_initialise().unwrap();
        let dumbell: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v3;
                v3 -> v3;v3 -> v3;
                v2 -> v2; v2 -> v2;
                v1 -> v1;v1 -> v1;
            },"scalars"
        )?;

        let f = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &GenerationSettings::default(),
        )?;

        insta::assert_snapshot!(
            f.graph.n_nodes(),
            @"64");

        let f = f.unfold_uncached();
        insta::assert_snapshot!(
            f.graph.n_nodes(),
            @"307");

        Ok(())
    }

    mod failing {
        use super::*;

        #[test]
        fn lobsided_double_dumbell() -> Result<()> {
            test_initialise().unwrap();
            let dumbell: Graph = dot!(
                digraph G{
                    edge [particle="scalar_1"];
                    v1 -> v2;
                    v2 -> v2; //v2 -> v2;
                    v1 -> v1;v1 -> v1;
                },"scalars"
            )?;

            let spinneys: Vec<_> = dumbell
                .spinneys(&dumbell.full_filter())
                .into_iter()
                .filter_map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis).transpose())
                .collect::<Result<_>>()?;
            let f = Wood::new(
                CutStructure::empty(&dumbell),
                &dumbell,
                &GenerationSettings::default(),
            )?;

            println!("{}", f);

            insta::assert_snapshot!(
                f.graph.n_nodes(),
                @"8",
                // format!("Wood does not have correct number of spinneys: \n{}",f)
            );

            for (_, _, d) in f.graph.iter_nodes() {
                println!(
                    "//Node {}: \n{}",
                    d.subgraph.string_label(),
                    dumbell.dot(&d.subgraph)
                );
            }
            let _ff = OldWood::from_spinneys(spinneys, &dumbell); //.unfold(&g, &g.loop_momentum_basis);

            // println!("{}", ff.dot(&dumbell));

            let f = f.unfold();
            f.debug_walk();
            println!("{}", f);
            insta::assert_snapshot!(
                f.normalized_node_labels_with_cover("3L").join("\n"),
                @r###"
{36,F}: Approx(S_36(_))*Approx(S_F(_))
{3} · {36,F}: Approx((-1*S_3+S_F(_))*Approx(S_3(_)))*Approx(S_36(_))
{C} · {36,F}: Approx((-1*S_C+S_F(_))*Approx(S_C(_)))*Approx(S_36(_))
"###
            );
            insta::assert_snapshot!(
                f.graph.n_nodes(),
                @"12");

            let f = Wood::new(
                CutStructure::empty(&dumbell),
                &dumbell,
                &GenerationSettings::default(),
            )?
            .unfold_uncached();
            assert!(f.compute_store.entries.is_empty());
            insta::assert_snapshot!(
                f.normalized_node_labels_with_cover("3L").join("\n"),
                @r###"
{36,F}: Approx(S_36(_))*Approx(S_F(_))
{3} · {36,F}: Approx((-1*S_3+S_F(_))*Approx(S_3(_)))*Approx(S_36(_))
{C} · {36,F}: Approx((-1*S_C+S_F(_))*Approx(S_C(_)))*Approx(S_36(_))
"###
            );

            Ok(())
        }

        #[test]
        fn double_double_dumbell() -> Result<()> {
            test_initialise().unwrap();
            let dumbell: Graph = dot!(
                digraph G{
                    edge [particle="scalar_1"];
                    v1 -> v2;
                    v2 -> v2;v2 -> v2; //v2 -> v2;
                    v1 -> v1;v1 -> v1;
                },"scalars"
            )?;

            let _spinneys: Vec<_> = dumbell
                .spinneys(&dumbell.full_filter())
                .into_iter()
                .filter_map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis).transpose())
                .collect::<Result<_>>()?;
            let f = Wood::new(
                CutStructure::empty(&dumbell),
                &dumbell,
                &GenerationSettings::default(),
            )?;

            println!("{}", f);

            insta::assert_snapshot!(
                f.graph.n_nodes(),
                @"16",
                // format!("Wood does not have correct number of spinneys: \n{}",f)
            );

            let f = f.unfold();

            insta::assert_snapshot!(
                f.graph.n_nodes(),
                @"36");

            let f = Wood::new(
                CutStructure::empty(&dumbell),
                &dumbell,
                &GenerationSettings::default(),
            )?
            .unfold_uncached();

            let foata_labels = f
                .graph
                .iter_nodes()
                .map(|(_, _, key)| key.foata_level_labels())
                .collect::<Vec<_>>();
            assert!(
                foata_labels
                    .iter()
                    .any(|label| label == "3,36;FU,F" || label == "36,3;FU,F"),
                "expected 3;F and 36;FU branch histories to combine, got:\n{}",
                foata_labels.join("\n")
            );

            Ok(())
        }
    }
}

#[cfg(test)]
#[path = "overlap_poles.rs"]
pub(crate) mod overlap_poles;
