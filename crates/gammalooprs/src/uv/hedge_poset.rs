use std::{
    cmp::Reverse,
    collections::{BTreeMap, BTreeSet},
    fmt::Display,
};

use ahash::AHashMap;
use eyre::{WrapErr, eyre};
use gammaloop_tracing_filter::LogMessage;
use idenso::{color::ColorSimplifier, shorthands::schoonschip::Schoonschip};
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph, NoData, NodeIndex,
    algorithms::trace_unfold::{
        HiddenData, Independence, TraceKey, TraceUnfold, UnfoldedTraceGraph,
    },
    involution::{EdgeIndex, Flow, HedgePair},
    nodestore::{NodeStorageOps, NodeStorageVec},
    subgraph::{Inclusion, InternalSubGraph, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    function,
    printer::PrintOptions,
    symbol,
};
use tracing::debug;
use vakint::Vakint;

use crate::{
    cff::ResidueSelectedTerms,
    graph::{Graph, LMBext, LoopMomentumBasis, cuts::CutSet, parse::string_utils::ToOrderedSimple},
    momentum::Sign,
    utils::{GS, W_},
    uv::{
        RenormalizationPart, Spinney, UVgenerationSettings, UltravioletGraph,
        approx::{
            ApproximationKernel, CutStructure, ForestNodeLike, OrientationProjection,
            ResidueProjection, UVCtx,
            final_integrand::{FinalIntegrand, IntegratedCtTerms},
            integrated::Integrated,
            local_3d::Local3DApproximation,
        },
        export::UVForestNodeExpression,
        forest::ParametricIntegrands,
        settings::VakintSettings,
    },
};
use color_eyre::Result;
use spenso::shadowing::symbolica_utils::LogPrint;

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

    pub(crate) fn new(cuts: CutStructure, graph: &Graph, settings: &UVgenerationSettings) -> Self {
        let mut subgraph = graph.full_filter();
        subgraph.subtract_with(&graph.initial_state_cut.left);
        let mut spinneys = Vec::new();

        for cut in cuts.cuts.iter() {
            let cut_sub = subgraph.subtract(&cut.union);
            spinneys.extend(graph.classified_spinneys(
                &cut_sub,
                settings,
                &graph.loop_momentum_basis,
            ));
        }

        Self::from_spinneys(spinneys, graph, cuts, &settings.vakint)
    }

    pub(crate) fn from_spinneys<I: IntoIterator<Item = Spinney>>(
        s: I,
        graph: &Graph,
        cuts: CutStructure,
        vakint_settings: &VakintSettings,
    ) -> Self {
        let mut max_loops = 0;
        let mut spinneys = BTreeMap::new();
        for spinney in s {
            max_loops = max_loops.max(spinney.max_comp_loop_count());
            spinneys.entry(spinney.filter().clone()).or_insert(spinney);
        }
        let empty = Spinney::empty(graph);
        spinneys.entry(empty.filter().clone()).or_insert(empty);
        let mut vakint_settings = vakint_settings.true_settings();
        // Set the number of terms in epsilon expansion to max number of loops across all components + 1
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

        // Not quite transitive closure. For disjoint unions, only keep edges that add a
        // single connected component of the sink; those are the only ones that can be
        // composed canonically by trace unfolding.
        for u in unions {
            // println!("//{u}:{}", poset[u].subgraph.string_label());
            let mut comps: BTreeSet<_> = graph
                .as_ref()
                .connected_components(&poset[u].subgraph)
                .into_iter()
                .collect();
            for c in poset.iter_crown(u) {
                let Flow::Sink = poset.flow(c) else {
                    continue;
                };
                let edge_id = poset[&c];
                if comps.contains(&poset[edge_id]) {
                    comps.remove(&poset[edge_id]);
                } else {
                    to_remove.add(c);
                    to_remove.add(poset.inv(c));
                }
            }
        }

        poset.delete_hedges(&to_remove);
        let root = poset
            .iter_nodes()
            .find(|(_, _, s)| s.subgraph.is_empty())
            .map(|(n, _, _)| n);

        Wood {
            graph: poset,
            root: root.expect("no empty spinney found"),
            cuts,
            vakint_settings,
        }
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
            cached_node_label_atoms: cache_node_label_atoms,
            compute_store: ComputeStore::default(),
            wood: self,
        };

        if cache_node_label_atoms {
            forests.with_cached_node_label_atoms()
        } else {
            forests
        }
    }

    pub fn unfold(self) -> Forests {
        self.unfold_with_cached_node_label_atoms(true)
    }

    pub fn unfold_uncached(self) -> Forests {
        self.unfold_with_cached_node_label_atoms(false)
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
            &|v| Some(format!("label=\"{}\"", v.subgraph.string_label())),
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
    pub kernel_hits: usize,
}

impl ComputeStore {
    fn get(&self, key: &OperationNode) -> Option<&ComputeNode> {
        self.entries.get(key)
    }

    fn require(&self, key: &OperationNode) -> Result<&ComputeNode> {
        self.get(key)
            .ok_or_else(|| eyre!("{} not yet added to store", key))
    }

    fn local_terms(&self, key: &OperationNode, cutset: &CutSet) -> Result<&ResidueSelectedTerms> {
        let computed = self.require(key)?;
        let integrands = computed
            .local_3d
            .get(cutset)
            .ok_or_else(|| eyre!("{} local_3d not computed for current cutset", key))?;
        let Integrands::Multiple(terms) = integrands else {
            return Err(eyre!("{} local_3d not computed yet", key));
        };
        Ok(terms)
    }

    fn entry(
        &mut self,
        key: OperationNode,
    ) -> std::collections::hash_map::Entry<'_, OperationNode, ComputeNode> {
        self.entries.entry(key)
    }

    fn record_kernel_hit(&mut self) {
        self.kernel_hits += 1;
    }
}

pub struct Forests {
    pub graph: UnfoldedTraceGraph<Wood, OperationNode, NodeStorageVec<OperationNode>>,
    pub root: NodeIndex,
    /// Wood subgraph that has compatible
    cuts: Vec<(SuBitGraph, CutSet)>,
    cached_node_label_atoms: bool,
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

struct LocalLeafOperation {
    op: HiddenData<SuBitGraph, EdgeIndex>,
    frontier: OperationNode,
}

impl LocalLeafOperation {
    fn new(op: &HiddenData<SuBitGraph, EdgeIndex>, frontier: &OperationNode) -> Self {
        Self {
            op: op.clone(),
            frontier: frontier.clone(),
        }
    }

    fn frontier_depth(&self) -> usize {
        self.frontier.key.op_count()
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

    fn forest_node(&self, graph: &Graph, topo_order: usize) -> OwnedForestNode {
        let spinney = match self.covers() {
            Some(cover) if !cover.is_empty() => {
                let subgraph = InternalSubGraph::cleaned_filter_optimist(cover, graph.as_ref());
                Spinney::new(subgraph, graph, &graph.loop_momentum_basis)
                    .expect("operation cover should define a valid spinney")
            }
            _ => Spinney::empty(graph),
        };

        OwnedForestNode {
            spinney,
            topo_order,
        }
    }

    pub fn to_atom(&self) -> Atom {
        let mut acc = Atom::one();

        let approx = FunctionBuilder::new(GS.t_op);
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

    fn integrated_4d_from_store(
        &self,
        compute_store: &ComputeStore,
        settings: &UVgenerationSettings,
    ) -> Result<Option<(ResidueSelectedTerms, Sign)>> {
        if !settings.generate_integrated || self.key.is_empty() {
            return Ok(None);
        }

        let computed = compute_store.require(self)?;
        let Integrand::Single(atom) = &computed.integrated_4d else {
            return Err(eyre!("{} integrated_4d not computed yet", self));
        };

        Ok(Some((
            ResidueSelectedTerms::all_none(atom.clone()),
            Sign::Positive,
        )))
    }

    pub fn integrated(
        &self,
        graph: &Graph,
        compute_store: &mut ComputeStore,
        wood: &Wood,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
    ) -> Result<Atom> {
        let mut acc = Atom::one();
        let integrated_orchestrator = Integrated::new(vakint, &wood.vakint_settings);
        let uvctx = UVCtx { graph, settings };

        let mut order = 0;
        let mut levels = self.key.view();
        if settings.cached
            && let Some((prefix, leaf_level)) = self.key.split_last_level()
        {
            let prefix_key = OperationNode {
                key: prefix.to_owned(),
            };
            if let Some(computed) = compute_store.get(&prefix_key) {
                let Integrand::Single(cached) = &computed.integrated_4d else {
                    return Err(eyre!("{} integrated_4d not computed yet", prefix_key));
                };
                acc = cached.clone();
                order = prefix.op_count();
                levels = leaf_level;
            }
        }

        for l in levels.iter_levels_top_down() {
            let mut mul = Atom::one();

            for op in l.iter_leaf_ops() {
                let (current, given) = wood.current_given_pair(op.data, order);
                order += 1;
                compute_store.record_kernel_hit();
                let raw_integrated =
                    integrated_orchestrator.kernel(&uvctx, &current, &given, &acc)?;
                let integrated = if given.subgraph().is_empty() {
                    raw_integrated
                } else {
                    -raw_integrated
                };
                mul *= integrated;
            }

            acc = mul
        }

        Ok(acc)
    }
}

pub enum Integrand {
    NotComputed,
    Single(Atom),
}
pub enum Integrands {
    NotComputed,
    Multiple(ResidueSelectedTerms),
}

pub struct ComputeNode {
    pub local_3d: AHashMap<CutSet, Integrands>,
    pub final_integrand: Integrands,
    pub integrated_4d: Integrand, //4d
    pub simple: Integrand,
    pub node_label_atom: Option<Atom>,
}

impl Default for ComputeNode {
    fn default() -> Self {
        ComputeNode {
            local_3d: AHashMap::new(),
            final_integrand: Integrands::NotComputed,
            integrated_4d: Integrand::NotComputed,
            simple: Integrand::NotComputed,
            node_label_atom: None,
        }
    }
}

impl ResidueProjection<'_> {
    fn initial_accumulator(
        &self,
        compute_store: &ComputeStore,
        frontier: &OperationNode,
        operation: &OperationNode,
        root_terms: &ResidueSelectedTerms,
    ) -> Result<ResidueSelectedTerms> {
        if frontier.key.is_empty() {
            return Ok(root_terms.clone());
        }

        compute_store
            .local_terms(frontier, self.cutset)
            .cloned()
            .wrap_err_with(|| {
                format!(
                    "while loading local frontier {} for {}",
                    frontier, operation
                )
            })
    }

    fn apply_local_step(
        &self,
        graph: &mut Graph,
        settings: &UVgenerationSettings,
        terms: &ResidueSelectedTerms,
        current: &ForestNode<'_>,
        given: &ForestNode<'_>,
        integrated_finite: Option<&Atom>,
    ) -> Result<ResidueSelectedTerms> {
        let integrated_terms = integrated_finite
            .map(|finite| {
                FinalIntegrand::new(*self, None, false)
                    .localized_finite_integrated_ct(graph, given, finite)
            })
            .transpose()?;
        let uvctx = UVCtx { graph, settings };
        let local = Local3DApproximation::full();

        terms.try_map(|index, atom| {
            let mut stepped = -local.kernel(&uvctx, current, given, atom)?;
            if let Some(integrated_terms) = &integrated_terms {
                let active_term = integrated_terms
                    .active
                    .get(index)
                    .wrap_err("while applying hedge-poset local subtraction")?;
                let frozen = integrated_terms
                    .frozen_integrands
                    .get(index)
                    .wrap_err("while applying hedge-poset local subtraction")?;
                stepped -=
                    Local3DApproximation::reduced().kernel(&uvctx, current, given, active_term)?
                        * frozen;
            }
            Ok(stepped)
        })
    }
}

impl Forests {
    fn compatible_topological_order(&self, subset: &SuBitGraph) -> Result<Vec<NodeIndex>> {
        let mut order = self.graph.topo_sort_kahn_of(subset)?;

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

    fn local_leaf_operations(&self, node: NodeIndex) -> Vec<LocalLeafOperation> {
        let mut leaves = self
            .graph
            .leaf_op_dependency_frontiers(node, &self.wood)
            .map(|(op, frontier)| LocalLeafOperation::new(op, &self.graph[frontier]))
            .collect::<Vec<_>>();

        leaves.sort_by_key(|leaf| {
            (
                Reverse(leaf.frontier_depth()),
                leaf.op.order.clone(),
                usize::from(leaf.op.data),
            )
        });
        leaves
    }

    fn dependent_integrated_finite(
        &self,
        frontier: &OperationNode,
        settings: &UVgenerationSettings,
    ) -> Result<Option<Atom>> {
        frontier
            .integrated_4d_from_store(&self.compute_store, settings)?
            .map(|integrated_4d| {
                let integrated_ct = IntegratedCtTerms::from(Some(integrated_4d));
                FinalIntegrand::finite_integrated_ct(integrated_ct)
            })
            .transpose()
    }

    fn local_3d_for_node(
        &mut self,
        node: NodeIndex,
        graph: &mut Graph,
        projection: ResidueProjection<'_>,
        root_terms: &ResidueSelectedTerms,
        settings: &UVgenerationSettings,
    ) -> Result<ResidueSelectedTerms> {
        let operation = self.graph[node].clone();
        let leaves = self.local_leaf_operations(node);
        let Some(first_leaf) = leaves.first() else {
            if operation.key.is_empty() {
                return Ok(root_terms.clone());
            }
            return Err(eyre!("{} has no local leaf operations", operation));
        };

        let mut acc = projection.initial_accumulator(
            &self.compute_store,
            &first_leaf.frontier,
            &operation,
            root_terms,
        )?;

        for (order, leaf) in leaves.iter().enumerate() {
            let integrated_finite = self.dependent_integrated_finite(&leaf.frontier, settings)?;
            let (current, given) = self.wood.current_given_pair(leaf.op.data, order);
            self.compute_store.record_kernel_hit();
            acc = projection.apply_local_step(
                graph,
                settings,
                &acc,
                &current,
                &given,
                integrated_finite.as_ref(),
            )?;
        }

        Ok(acc)
    }

    fn cached_node_label_atom(&self, node: NodeIndex) -> Option<Atom> {
        self.compute_store
            .get(&self.graph[node])
            .and_then(|computed| computed.node_label_atom.as_ref())
            .cloned()
    }

    fn node_label_atom_factor(
        &self,
        frontier: NodeIndex,
        op: &HiddenData<SuBitGraph, EdgeIndex>,
    ) -> Atom {
        let approx = FunctionBuilder::new(GS.t_op);
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

        if self.cached_node_label_atoms
            && let Some(cached) = self.cached_node_label_atom(node)
        {
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
        for nidx in self.graph.topo_sort_kahn().unwrap() {
            let atom = self.node_label_atom(nidx);

            self.compute_store
                .entry(self.graph[nidx].clone())
                .or_default()
                .node_label_atom = Some(atom);
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
        Atom::var(symbol!("expr"))
    }

    fn dot_serialize_node_atom_factor(
        &self,
        frontier: NodeIndex,
        op: &HiddenData<SuBitGraph, EdgeIndex>,
    ) -> Atom {
        function!(
            GS.t_op,
            function!(op.order.symbol(), usize::from(op.data)),
            self.graph[frontier]
                .covers()
                .map_or(Atom::Zero, |cover| cover.symbol().to_atom()),
            self.dot_serialize_node_atom(frontier)
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
            .printer(PrintOptions::typst())
            .to_string();
        let cover = key
            .covers()
            .unwrap_or_else(|| self.graph.empty_subgraph())
            .string_label();

        format!(
            "label={} foata={} cover={}",
            Self::dot_attr_value(&label),
            Self::dot_attr_value(&key.foata_level_labels()),
            Self::dot_attr_value(&cover),
        )
    }

    fn dot_attr_value(value: &str) -> String {
        let mut escaped = String::with_capacity(value.len() + 2);
        escaped.push('"');
        for c in value.chars() {
            match c {
                '\\' => escaped.push_str("\\\\"),
                '"' => escaped.push_str("\\\""),
                '\n' => escaped.push_str("\\n"),
                '\r' => escaped.push_str("\\r"),
                '\t' => escaped.push_str("\\t"),
                _ => escaped.push(c),
            }
        }
        escaped.push('"');
        escaped
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

    pub fn integrate(
        &mut self,
        graph: &Graph,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        for (order, nidx) in self.graph.topo_sort_kahn().unwrap().iter().enumerate() {
            debug!(order=%order,cache=%settings.cached,nidx=%nidx,key=%self.graph[*nidx],"One integrated step");
            let integrand = self.graph[*nidx].integrated(
                graph,
                &mut self.compute_store,
                &self.wood,
                vakint,
                settings,
            )?;
            self.compute_store
                .entry(self.graph[*nidx].clone())
                .or_default()
                .integrated_4d = Integrand::Single(integrand);
        }

        Ok(())
    }

    pub(crate) fn compute(
        &mut self,
        graph: &mut Graph,
        vakint: &Vakint,
        orientation: OrientationProjection<'_>,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        if settings.generate_integrated {
            self.integrate(graph, vakint, settings)?;
        }

        for (cut_compatible_forest_subset, cuts) in self.cuts.clone() {
            let projection = orientation.for_cut(&cuts);
            let root = Local3DApproximation::root(graph, projection)?;

            for (order, nidx) in self
                .compatible_topological_order(&cut_compatible_forest_subset)?
                .into_iter()
                .enumerate()
            {
                debug!(order=%order,cache=%settings.cached,nidx=%nidx,key=%self.graph[nidx],"One integrated step");
                let operation = self.graph[nidx].clone();
                let integrands =
                    self.local_3d_for_node(nidx, graph, projection, &root, settings)?;
                let forest_node = operation.forest_node(graph, order);
                let integrated_ct = IntegratedCtTerms::from(
                    operation.integrated_4d_from_store(&self.compute_store, settings)?,
                );
                let uv_marker =
                    (settings.add_sigma && settings.keep_sigma).then(|| operation.to_atom());
                let final_integrand = FinalIntegrand::new(projection, uv_marker.as_ref(), true)
                    .build(
                        graph,
                        &forest_node,
                        &integrands,
                        Sign::Positive,
                        integrated_ct,
                    )?;

                self.compute_store
                    .entry(operation)
                    .or_default()
                    .local_3d
                    .insert(cuts.clone(), Integrands::Multiple(integrands));
                self.compute_store
                    .entry(self.graph[nidx].clone())
                    .or_default()
                    .final_integrand = Integrands::Multiple(final_integrand);
            }
        }

        Ok(())
    }

    pub(crate) fn orientation_parametric_exprs(
        &self,
        graph: &Graph,
        settings: &UVgenerationSettings,
    ) -> Result<Vec<ParametricIntegrands>> {
        let mut exprs = Vec::new();

        for (cut_compatible_forest_subset, cuts) in &self.cuts {
            let mut sum: Option<ResidueSelectedTerms> = None;

            for nidx in self.compatible_topological_order(cut_compatible_forest_subset)? {
                let operation = &self.graph[nidx];
                let computed = self.compute_store.require(operation)?;
                let Integrands::Multiple(final_integrand) = &computed.final_integrand else {
                    return Err(eyre!("{} final_integrand not computed yet", operation));
                };

                let sum = sum.get_or_insert_with(ResidueSelectedTerms::empty);
                for (cut_index, integrand) in final_integrand {
                    sum.add_assign_term(*cut_index, integrand.clone());
                }
            }

            let mut integrands = sum.ok_or_else(|| eyre!("No terms in hedge-poset forest"))?;

            for (pair, edge_index, _) in graph.iter_edges_of(
                &graph
                    .full_filter()
                    .subtract(&graph.initial_state_cut)
                    .subtract(&graph.tree_edges),
            ) {
                if matches!(pair, HedgePair::Unpaired { .. }) {
                    continue;
                }

                for expr in integrands.values_mut() {
                    *expr = expr.replace_multiple(&[GS.split_mom_pattern_simple(edge_index)]);
                }
            }

            for expr in integrands.values_mut() {
                *expr = expr
                    .replace(function!(GS.den, W_.a_, W_.b_, W_.c_, W_.d_))
                    .with(W_.d_);
                if !settings.keep_sigma {
                    *expr = expr
                        .replace(function!(GS.uv_local, W_.a___))
                        .with(Atom::num(1))
                        .replace(function!(GS.uv_integrated, W_.a___))
                        .with(Atom::num(1));
                }
            }

            exprs.push(ParametricIntegrands {
                integrands,
                cuts: cuts.clone(),
            });
        }

        Ok(exprs)
    }

    pub(crate) fn export_node_expressions(
        &self,
        forest_index: usize,
        post_process: &mut impl FnMut(Atom) -> Atom,
    ) -> Result<Vec<UVForestNodeExpression>> {
        let (cut_compatible_forest_subset, _) = self
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
            let Integrands::Multiple(final_integrand) = &computed.final_integrand else {
                return Err(eyre!("{} final_integrand not computed yet", operation));
            };
            let node_key = operation.to_string();
            for (term_index, (&residue_index, numerator)) in final_integrand.iter().enumerate() {
                terms.push(UVForestNodeExpression {
                    forest_index,
                    node_index,
                    node_key: node_key.clone(),
                    term_index,
                    residue_index,
                    numerator: post_process(numerator.clone()),
                });
            }
        }

        Ok(terms)
    }

    // pub fn compute(
    //     &mut self,
    //     graph: &mut Graph,
    //     wood: &Wood,
    //     settings: &UVgenerationSettings,
    // ) -> Result<()> {
    //     let local_orchestrator = Local3DApproximation::full();
    //     for (cut_compatible_forest_subset, c) in &self.cuts {
    //         let mut integrands = Some(Local3DApproximation::root(graph, c)?);

    //         let uvctx = UVCtx { graph, settings };
    //         for (order, nidx) in self
    //             .graph
    //             .topo_sort_kahn_of(cut_compatible_forest_subset)
    //             .unwrap()
    //             .iter()
    //             .enumerate()
    //         {
    //             for h in self.iter_parents(*nidx, order, wood) {
    //                 let (computed, current, given, parent_key, is_union) = h?;

    //                 let Integrands::Multiple(a) = &computed.local_3d else {
    //                     return Err(eyre!("{} integrated_4d not computed yet", parent_key));
    //                 };

    //                 if is_union {
    //                     // integrand *= a;
    //                 } else {
    //                     integrands = Some(
    //                         a.iter()
    //                             .map(|a| local_orchestrator.kernel(&uvctx, &current, &given, a))
    //                             .collect::<Result<_>>()?,
    //                     );
    //                 }
    //             }

    //             self.compute_store
    //                 .entry(self.graph[*nidx].clone())
    //                 .or_default()
    //                 .local_3d = Integrands::Multiple(
    //                 integrands
    //                     .take()
    //                     .ok_or(eyre!("Taken integrand not filled in"))?,
    //             );
    //         }
    //     }
    //     Ok(())
    // }

    pub(crate) fn pole_part_of_ends(&self, graph: &Graph) -> Result<RenormalizationPart> {
        let mut sum = Atom::Zero;

        let wild = Atom::var(W_.x___);

        let replacements =
            graph.integrand_replacement(&graph.full_filter(), &graph.loop_momentum_basis, &[wild]);
        for (_, mut crown, key) in self.graph.iter_nodes() {
            if crown.any(|r| self.graph.flow(r).is_source()) {
                continue;
            }

            let computed = self
                .compute_store
                .get(key)
                .ok_or(eyre!("{} not yet added to store", key))?;

            let Integrand::Single(atom) = &computed.integrated_4d else {
                return Err(eyre!("{} integrated_4d not computed yet", key));
            };
            debug!(
                key=%key,
               expr = % atom.expand_num().log_print(None),"Term before simplification"
            );
            let atom = (atom
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

        let n_loops = graph.n_loops(&graph.full_filter());
        let pole_stripped = sum
            .series(GS.dim_epsilon, Atom::Zero, n_loops as i64 + 1)
            .unwrap();

        sum = Atom::Zero;

        for (power, p) in pole_stripped.terms() {
            if power < 0 {
                sum += p * Atom::var(GS.dim_epsilon).pow(power);
            }
        }

        Ok(RenormalizationPart::new(
            sum.replace_multiple(&replacements)
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum),
            self.compute_store.kernel_hits,
            self.graph.n_nodes(),
        ))
    }

    pub fn debug_walk(&self) {
        let mut cover_groups: BTreeMap<SuBitGraph, Vec<NodeIndex>> = BTreeMap::new();

        self.graph
            .topo_sort_kahn()
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
    use crate::{
        cff::CutCFFIndex,
        dot,
        graph::{Graph, parse::IntoGraph},
        initialisation::test_initialise,
        processes::DotExportSettings,
        settings::global::OrientationPattern,
        uv::{UltravioletGraph, Wood as OldWood},
    };

    use super::*;
    use color_eyre::Result;
    use linnet::half_edge::involution::Hedge;
    use symbolica::symbol;

    fn test_subgraph(n_edges: usize, edges: &[usize]) -> SuBitGraph {
        let mut subgraph = SuBitGraph::empty(2 * n_edges);
        for edge in edges {
            subgraph.add(Hedge(2 * *edge));
            subgraph.add(Hedge(2 * *edge + 1));
        }
        subgraph
    }

    fn test_op(order: SuBitGraph, data: usize) -> HiddenData<SuBitGraph, EdgeIndex> {
        HiddenData {
            order,
            data: EdgeIndex(data),
        }
    }

    impl Forests {
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
            &UVgenerationSettings::default(),
        )
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
                .map(|leaf| (leaf.op.order.string_label(), frontier_label(&leaf.frontier)))
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
        Ok(())
    }

    #[test]
    fn local_frontier_terms_read_non_root_compute_store() -> Result<()> {
        test_initialise().unwrap();
        let graph: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v2;
                v1 -> v1;
            },"scalars"
        )?;
        let cut_structure = CutStructure::empty(&graph);
        let cutset = cut_structure
            .cuts
            .first()
            .expect("empty cut structure has one cut")
            .clone();
        let root = ResidueSelectedTerms::all_none(Atom::var(symbol!("root")));
        let cached = ResidueSelectedTerms::all_none(Atom::var(symbol!("cached")));
        let frontier = OperationNode {
            key: TraceKey::from_levels(vec![vec![test_op(test_subgraph(4, &[1]), 0)]]),
        };
        let operation = OperationNode {
            key: TraceKey::from_levels(vec![
                vec![test_op(test_subgraph(4, &[1]), 0)],
                vec![test_op(test_subgraph(4, &[0, 1]), 1)],
            ]),
        };
        let mut compute_store = ComputeStore::default();
        compute_store
            .entry(frontier.clone())
            .or_default()
            .local_3d
            .insert(cutset.clone(), Integrands::Multiple(cached.clone()));
        let orientation_pattern = OrientationPattern::default();
        let orientation = OrientationProjection::new(&[], &orientation_pattern);
        let projection = orientation.for_cut(&cutset);

        assert_eq!(
            projection
                .initial_accumulator(&compute_store, &frontier, &operation, &root)?
                .get(CutCFFIndex::new_all_none())?
                .to_canonical_string(),
            cached
                .get(CutCFFIndex::new_all_none())?
                .to_canonical_string()
        );
        assert_eq!(
            projection
                .initial_accumulator(
                    &compute_store,
                    &OperationNode {
                        key: TraceKey::empty(),
                    },
                    &operation,
                    &root,
                )?
                .get(CutCFFIndex::new_all_none())?
                .to_canonical_string(),
            root.get(CutCFFIndex::new_all_none())?.to_canonical_string()
        );
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
            .filter_map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis))
            .collect();
        let f = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &UVgenerationSettings::default(),
        );

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
        println!("{}", f);
        insta::assert_snapshot!(
            f.graph.n_nodes(),
            @"8");

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
            &UVgenerationSettings::default(),
        );

        println!("{}", dt.dot_serialize(&DotExportSettings::default()));

        insta::assert_snapshot!(
            f.graph.n_nodes(),
            @"5",
            // format!("Wood does not have correct number of spinneys: \n{}",f)
        );

        let f = f.unfold();
        println!("{}", f.dot_serialize());
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
            .filter_map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis))
            .collect();
        let f = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &UVgenerationSettings::default(),
        );

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
                    .filter_map(|a| Spinney::new(a, &g, &g.loop_momentum_basis))
                    .collect();
                let f = Wood::new(
                    CutStructure::empty(&g),
                    &g,
                    &UVgenerationSettings::default(),
                );

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
            &UVgenerationSettings::default(),
        );
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
            &UVgenerationSettings::default(),
        );
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
            &UVgenerationSettings::default(),
        );
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
            &UVgenerationSettings::default(),
        );
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
        let f = Wood::new(
            CutStructure::empty(&spectacles),
            &spectacles,
            &UVgenerationSettings::default(),
        );
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
    fn basketball() -> Result<()> {
        test_initialise().unwrap();

        let basketball: Graph = dot!(
            digraph basketball{
                node [num = "1" ]
                edge [particle=scalar_1]
                e    [style=invis]

                e -> A:1   [ id=5]
                B:0 -> e   [ id=4]
                A -> B    [ id=0 lmb_id=0]
                A -> B    [ id=1 lmb_id=1]
                A -> B    [ id=2 lmb_id=2]
                A -> B    [ id=3 ]
            },"scalars"
        )?;

        let f = Wood::new(
            CutStructure::empty(&basketball),
            &basketball,
            &UVgenerationSettings::default(),
        );
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
            &UVgenerationSettings::default(),
        );
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
            &UVgenerationSettings::default(),
        );
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
            &UVgenerationSettings::default(),
        );

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
        use crate::processes::DotExportSettings;

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
                .filter_map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis))
                .collect();
            let f = Wood::new(
                CutStructure::empty(&dumbell),
                &dumbell,
                &UVgenerationSettings::default(),
            );

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
{36,F}: T(S_36(_))*T(S_F(_))
{3} · {36,F}: T((-1*S_3+S_F(_))*T(S_3(_)))*T(S_36(_))
{C} · {36,F}: T((-1*S_C+S_F(_))*T(S_C(_)))*T(S_36(_))
"###
            );
            insta::assert_snapshot!(
                f.graph.n_nodes(),
                @"12");

            let f = Wood::new(
                CutStructure::empty(&dumbell),
                &dumbell,
                &UVgenerationSettings::default(),
            )
            .unfold_uncached();
            assert!(f.compute_store.entries.is_empty());
            insta::assert_snapshot!(
                f.normalized_node_labels_with_cover("3L").join("\n"),
                @r###"
{36,F}: T(S_36(_))*T(S_F(_))
{3} · {36,F}: T((-1*S_3+S_F(_))*T(S_3(_)))*T(S_36(_))
{C} · {36,F}: T((-1*S_C+S_F(_))*T(S_C(_)))*T(S_36(_))
"###
            );

            println!(
                "{}",
                dumbell.dot_serialize(&DotExportSettings {
                    ..Default::default()
                })
            );
            println!("{}", f.dot_serialize());

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
                .map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis))
                .collect();
            let f = Wood::new(
                CutStructure::empty(&dumbell),
                &dumbell,
                &UVgenerationSettings::default(),
            );

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
                &UVgenerationSettings::default(),
            )
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

            println!(
                "{}",
                dumbell.dot_serialize(&DotExportSettings {
                    ..Default::default()
                })
            );
            println!("{}", f.dot_serialize());

            Ok(())
        }
    }
}
