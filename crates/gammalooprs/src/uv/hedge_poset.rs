use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Display,
};

#[cfg(test)]
use std::cmp::Reverse;

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
};
use tracing::debug;
use vakint::Vakint;

use crate::{
    graph::{Graph, LMBext, LoopMomentumBasis, cuts::CutSet, parse::string_utils::ToOrderedSimple},
    utils::{GS, W_},
    uv::{
        ApproximationType, Integrands, RenormalizationPart, Spinney, UVgenerationSettings,
        UltravioletGraph,
        approx::{
            CutStructure, ForestNodeLike, OrientationProjection, Rooted, UVCtx,
            final_integrand::{FinalIntegrandBuilder, FinalIntegrands},
            integrated::{Integrated, IntegratedCts},
            local_3d::{Local3DApproximation, Local3DCts, Localizer},
            local_4d::{self, Full4DCts, Local4dCts},
        },
        export::UVForestNodeExpression,
        forest::ParametricIntegrands,
        marker::UvMarker,
        settings::VakintSettings,
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
            max_loops = max_loops.max(graph.n_loops(spinney.filter()));
            spinneys.entry(spinney.filter().clone()).or_insert(spinney);
        }
        let empty = Spinney::empty(graph);
        spinneys.entry(empty.filter().clone()).or_insert(empty);
        let mut vakint_settings = vakint_settings.true_settings();
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
            cached_node_label_atoms: cache_node_label_atoms.then(AHashMap::new),
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
}

impl ComputeStore {
    fn get(&self, key: &OperationNode) -> Option<&ComputeNode> {
        self.entries.get(key)
    }

    fn require(&self, key: &OperationNode) -> Result<&ComputeNode> {
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
    // `Local3DApproximation::run` applies the local subtraction signs directly.
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

    fn recursion_input_4d(&self, node: NodeIndex) -> Result<Full4DCts> {
        let operation = &self.graph[node];
        let computed = self.compute_store.require(operation)?;
        if self.graph.is_disjoint_union(node) {
            return Ok(Full4DCts::from_factorized_local(
                computed.local_4d(operation)?,
            ));
        }
        Full4DCts::recursion_input(
            computed.local_4d(operation)?,
            computed.integrated(operation)?,
            self.source_spinney(node).renormalization_scheme,
            node == self.root,
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
        for nidx in self.graph.topo_sort_kahn().unwrap() {
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
        let mut options = PrintOptions::typst();
        options.custom_print_mode = SpensoPrintSettings::typst().into();
        let label = self
            .dot_serialize_node_atom(node)
            .printer(options)
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

    fn union_replay_states(&self, node: NodeIndex) -> Result<Vec<UnionReplayState>> {
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
            let mut states = self.union_replay_states(parent)?;
            for state in &mut states {
                state.local_edges.push(edge);
            }
            states.push(UnionReplayState {
                integrated: node,
                local_edges: Vec::new(),
            });
            return Ok(states);
        }

        let component_states = self
            .disconnected_component_nodes(node)?
            .into_iter()
            .map(|component| self.union_replay_states(component))
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

    fn compute_4d_for_node(
        &self,
        node: NodeIndex,
        graph: &Graph,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
    ) -> Result<(Local4dCts, IntegratedCts)> {
        let operation = &self.graph[node];
        if operation.key.is_empty() {
            return Ok((Local4dCts::root(), IntegratedCts::root()));
        }

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

            let depth = graph.n_loops(
                &operation
                    .covers()
                    .expect("a non-root operation has a cover"),
            ) + 1;
            return Ok((
                Local4dCts::from_full_product(full_components),
                IntegratedCts::factorized_product(integrated_components, depth)?,
            ));
        }

        let (parent, edge) = self
            .graph
            .unique_parent(node)
            .ok_or_else(|| eyre!("{operation} has no unique 4D parent"))?;
        let full = self
            .recursion_input_4d(parent)
            .wrap_err_with(|| format!("while loading the 4D parent for {operation}"))?;
        let step_order = self.graph[parent].key.op_count();
        let (current, given) = self.wood.current_given_pair(edge, step_order);
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
    ) -> Result<CutComputation> {
        let operation = &self.graph[node];
        let local_3d = if operation.key.is_empty() {
            Local3DCts::root(graph, localizer)?
        } else if self.graph.is_disjoint_union(node) {
            let mut active_sectors = Vec::new();
            for state in self
                .union_replay_states(node)?
                .into_iter()
                .filter(|state| !state.local_edges.is_empty())
            {
                let integrated_operation = &self.graph[state.integrated];
                let mut edges = state.local_edges.into_iter().enumerate();
                let (offset, first_edge) = edges
                    .next()
                    .expect("a proper integrated prefix has a local suffix");
                let step_order = integrated_operation.key.op_count() + offset;
                let (current, given) = self.wood.current_given_pair(first_edge, step_order);

                // An empty integrated prefix starts from the per-cut root integrand. Every
                // other prefix enters through the reduced branch of its first local operation.
                let mut sector = if integrated_operation.key.is_empty() {
                    let root = Local3DCts::root(graph, localizer)?;
                    Local3DApproximation::new(localizer, graph, settings)
                        .run_local(&root, &current, &given, &current, &given)?
                } else {
                    let integrated = self
                        .compute_store
                        .require(integrated_operation)?
                        .integrated(integrated_operation)?;
                    let prefix_node = ForestNode {
                        spinney: self.source_spinney(state.integrated),
                        topo_order: integrated_operation.key.op_count(),
                    };
                    Local3DApproximation::new(localizer, graph, settings).run_integrated(
                        integrated,
                        &prefix_node,
                        &current,
                        &given,
                        &current,
                        &given,
                    )?
                };

                for (offset, edge) in edges {
                    let step_order = integrated_operation.key.op_count() + offset;
                    let (current, given) = self.wood.current_given_pair(edge, step_order);
                    sector = Local3DApproximation::new(localizer, graph, settings)
                        .run_local(&sector, &current, &given, &current, &given)?;
                }

                // Keep each active/frozen split after its root-path replay so any connected
                // descendants rescale only the loop variables still active in that sector.
                active_sectors.extend(
                    sector
                        .active_sectors()
                        .expect("a replayed union sector retains its active subgraph")
                        .iter()
                        .cloned(),
                );
            }

            Local3DCts::from_active_sectors(active_sectors)
                .wrap_err_with(|| format!("{operation} has no proper integrated prefixes"))?
        } else {
            let (parent, edge) = self
                .graph
                .unique_parent(node)
                .ok_or_else(|| eyre!("{operation} has no unique local parent"))?;
            let parent_operation = &self.graph[parent];
            // An empty dependency frontier starts from the per-cut root integrand;
            // otherwise its typed local result remains the sequential accumulator.
            let parent_local = if parent_operation.key.is_empty() {
                Local3DCts::root(graph, localizer)?
            } else {
                self.compute_store
                    .require(parent_operation)?
                    .cut(parent_operation, cutset)?
                    .local_3d
                    .clone()
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
            // `run` applies both subtraction signs; no external sign or raw
            // Foata-level product is introduced for an ordinary single-parent node.
            Local3DApproximation::new(localizer, graph, settings).run(
                &parent_local,
                parent_integrated,
                &current,
                &given,
                &current,
                &given,
            )?
        };

        let integrated = self
            .compute_store
            .require(operation)?
            .integrated(operation)?;
        let forest_node = ForestNode {
            spinney: self.source_spinney(node),
            topo_order: operation.key.op_count(),
        };
        let final_integrands = FinalIntegrandBuilder::new(localizer, settings).build_3d(
            graph,
            &forest_node,
            &local_3d,
            integrated,
        )?;

        Ok(CutComputation {
            local_3d,
            final_integrands,
        })
    }

    pub fn integrate(
        &mut self,
        graph: &Graph,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        for (order, nidx) in self.graph.topo_sort_kahn()?.into_iter().enumerate() {
            debug!(order, nidx=%nidx, key=%self.graph[nidx], "Computing hedge-poset 4D term");
            let operation = self.graph[nidx].clone();
            let (local_4d, integrated) = self.compute_4d_for_node(nidx, graph, vakint, settings)?;
            let computed = self.compute_store.entry(operation).or_default();
            computed.local_4d = Some(local_4d);
            computed.integrated = Some(integrated);
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
        self.integrate(graph, vakint, settings)?;

        for (compatible_subset, cutset) in self.cuts.clone() {
            let localizer = Localizer::new(&cutset, orientation);
            for (order, nidx) in self
                .compatible_topological_order(&compatible_subset)?
                .into_iter()
                .enumerate()
            {
                debug!(order, nidx=%nidx, key=%self.graph[nidx], "Computing hedge-poset per-cut term");
                let operation = self.graph[nidx].clone();
                let cut_computation =
                    self.local_3d_for_node(nidx, graph, &cutset, localizer, settings)?;
                self.compute_store
                    .entry(operation)
                    .or_default()
                    .cuts
                    .insert(cutset.clone(), cut_computation);
            }
        }

        Ok(())
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
            let mut sum: Option<Integrands> = None;
            for nidx in self.compatible_topological_order(compatible_subset)? {
                let operation = &self.graph[nidx];
                let terms: Integrands = self
                    .compute_store
                    .require(operation)?
                    .cut(operation, cutset)?
                    .final_integrands
                    .iter()
                    .map(|(index, integrand)| (*index, integrand.clone().collect_color()))
                    .collect();
                sum = Some(match sum {
                    Some(sum) => sum.zip_add(&terms).wrap_err_with(|| {
                        format!("while aggregating hedge-poset term {operation} for cut {cutset:?}")
                    })?,
                    None => terms,
                });
            }

            let integrands = sum
                .ok_or_else(|| eyre!("No terms in hedge-poset forest for cut {cutset:?}"))?
                .map(|integrand| {
                    integrand
                        .replace_multiple(&split_momentum_replacements)
                        .replace(function!(GS.den, W_.a_, W_.b_, W_.c_, W_.d_))
                        .with(W_.d_)
                });
            expressions.push(ParametricIntegrands {
                integrands,
                cuts: cutset.clone(),
            });
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
            let final_integrands = &computed.cut(operation, cutset)?.final_integrands;
            let node_key = operation.to_string();
            for (term_index, (&residue_index, numerator)) in final_integrands.iter().enumerate() {
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
    //     let local_orchestrator = Local3DApproximation {};
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

            let forest_node = key.forest_node(graph, key.key.op_count());
            let computed = self.compute_store.require(key)?;
            let integrated = computed.integrated(key)?;
            let physical = match self.source_spinney(node).renormalization_scheme {
                ApproximationType::MUV => integrated.physical_finite_counterterm_atom(),
                ApproximationType::PolePart => integrated.physical_pole_atom(),
                scheme => return Err(eyre!("No terminal counterterm projection for {scheme}")),
            };
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
        dot,
        graph::{Graph, parse::IntoGraph},
        initialisation::test_initialise,
        processes::DotExportSettings,
        uv::{UltravioletGraph, Wood as OldWood},
    };

    use super::*;
    use color_eyre::Result;

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
                let mut options = PrintOptions::typst();
                options.custom_print_mode = SpensoPrintSettings::typst().into();
                join_markers.push(
                    UvMarker::subgraph(current.subgraph(), given.subgraph())
                        .printer(options)
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
        let replay_states = forests.union_replay_states(dependent_disconnected)?;
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
    fn union_terms_replay_component_paths_from_typed_roots() -> Result<()> {
        test_initialise().unwrap();
        let mut graph: Graph = dot!(
            digraph G{
                edge [particle="scalar_1"];
                v1 -> v2;
                v2 -> v2;
                v1 -> v1;v1 -> v1;
            },"scalars"
        )?;
        let settings = UVgenerationSettings::default();
        let cut_structure = CutStructure::empty(&graph);
        let cutset = cut_structure
            .cuts
            .first()
            .expect("empty cut structure has one cut")
            .clone();
        let mut forests = Wood::new(cut_structure, &graph, &settings).unfold();

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
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::new(&[], &orientation_pattern),
        );
        let seed =
            forests.local_3d_for_node(forests.root, &mut graph, &cutset, localizer, &settings)?;
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
            &settings,
        )?;
        assert_eq!(
            root_result
                .local_3d
                .active_sectors()
                .expect("a root union keeps its active sectors")
                .len(),
            3
        );
        assert!(
            root_result
                .local_3d
                .integrands()
                .iter()
                .all(|(_, term)| !term.contains_symbol(root_store_marker)),
            "an empty dependency frontier must start from the typed root"
        );

        let frontier_result = forests.local_3d_for_node(
            dependent_disconnected,
            &mut graph,
            &cutset,
            localizer,
            &settings,
        )?;
        assert_eq!(
            frontier_result
                .local_3d
                .active_sectors()
                .expect("a dependent union keeps its active sectors")
                .len(),
            5
        );
        let replay_states = forests.union_replay_states(dependent_disconnected)?;
        let (state, (active_subgraph, _)) = replay_states
            .iter()
            .filter(|state| !state.local_edges.is_empty())
            .zip(
                frontier_result
                    .local_3d
                    .active_sectors()
                    .expect("a dependent union keeps its active sectors"),
            )
            .find(|(state, _)| {
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
        assert_eq!(active_subgraph, &expected_active);
        assert!(
            active_subgraph.empty_intersection(forests.source_spinney(state.integrated).filter())
        );
        assert!(
            frontier_result
                .local_3d
                .integrands()
                .iter()
                .all(|(_, term)| !term.contains_symbol(root_store_marker)
                    && !term.contains_symbol(frontier_store_marker)),
            "a union must replay every component path from its typed root"
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
        let replay_states = f.union_replay_states(three_component_union)?;
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
            &UVgenerationSettings::default(),
        );

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
        let settings = UVgenerationSettings::default();
        let cut_structure = CutStructure::empty(&spectacles);
        let cutset = cut_structure
            .cuts
            .first()
            .expect("empty cut structure has one cut")
            .clone();
        let f = Wood::new(cut_structure, &spectacles, &settings);
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

        let (union, edge) = f
            .graph
            .iter_nodes()
            .find_map(|(child, _, _)| {
                let (parent, edge) = f.graph.unique_parent(child)?;
                (!f.graph.is_disjoint_union(child) && f.graph.is_disjoint_union(parent))
                    .then_some((parent, edge))
            })
            .expect("spectacles has a connected child above its disconnected union");
        let orientation_pattern = crate::settings::global::OrientationPattern::default();
        let localizer = Localizer::new(
            &cutset,
            OrientationProjection::new(&[], &orientation_pattern),
        );
        let union_active = f
            .union_replay_states(union)?
            .into_iter()
            .filter(|state| !state.local_edges.is_empty())
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
        assert_eq!(union_active.len(), 3);
        let union_local = Local3DCts::from_active_sectors(
            union_active
                .iter()
                .cloned()
                .map(|active| (active, Integrands::root()))
                .collect(),
        )?;

        let step_order = f.graph[union].key.op_count();
        let (current, given) = f.wood.current_given_pair(edge, step_order);
        let reduced = current.reduced_subgraph(&given);
        let expected_active = union_active
            .iter()
            .map(|active| active.union(&reduced))
            .chain(std::iter::once(reduced.clone()))
            .collect::<Vec<_>>();
        let child_local = Local3DApproximation::new(localizer, &mut spectacles, &settings).run(
            &union_local,
            &IntegratedCts::root(),
            &current,
            &given,
            &current,
            &given,
        )?;
        let child_active = child_local
            .active_sectors()
            .expect("a connected child keeps its parent's active sectors")
            .iter()
            .map(|(active, _)| active.clone())
            .collect::<Vec<_>>();
        assert_eq!(child_active.len(), 4);
        assert_eq!(child_active, expected_active);

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

            Ok(())
        }
    }
}
