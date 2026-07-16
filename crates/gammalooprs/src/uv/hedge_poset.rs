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
    function, symbol,
};
use tracing::debug;
use vakint::Vakint;

use crate::{
    graph::{Graph, LMBext, LoopMomentumBasis, cuts::CutSet, parse::string_utils::ToOrderedSimple},
    utils::{GS, W_, symbolica_ext::LogPrint},
    uv::{
        Integrands, RenormalizationPart, Spinney, UVgenerationSettings, UltravioletGraph,
        approx::{
            CutStructure, ForestNodeLike, OrientationProjection, Rooted, UVCtx,
            final_integrand::{FinalIntegrandBuilder, FinalIntegrands},
            integrated::{Integrated, IntegratedCts},
            local_3d::{Local3DApproximation, Local3DCts, Localizer},
            local_4d::{self, Full4DCts, Local4dCts},
        },
        forest::ParametricIntegrands,
        marker::UvMarker,
        settings::VakintSettings,
    },
};
use color_eyre::Result;

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
        let mut spinneyset: BTreeSet<_> = s
            .into_iter()
            .inspect(|a| {
                max_loops = max_loops.max(a.max_comp_loop_count());
            })
            .collect();

        spinneyset.insert(Spinney::empty(&graph));
        let mut vakint_settings = vakint_settings.true_settings();
        // Set the number of terms in epsilon expansion to max number of loops across all components + 1
        vakint_settings.number_of_terms_in_epsilon_expansion = max_loops as i64 + 1;

        let mut unions = BTreeSet::new();
        let g: HedgeGraph<_, _> = HedgeGraph::poset(spinneyset);
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

        let approx = FunctionBuilder::new(symbol!("T"));
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
                Atom::var(symbol!(format!("S_{}", last.string_label())))
            };

            let mut mul = Atom::one();

            for op in l.iter_leaf_ops() {
                let new = function!(
                    symbol!(format!("S_{}", op.order.string_label())),
                    usize::from(op.data)
                );
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

    fn recursion_input_4d(&self, operation: &OperationNode, pole_part: bool) -> Result<Full4DCts> {
        Ok(Full4DCts::recursion_input(
            self.local_4d(operation)?,
            self.integrated(operation)?,
            pole_part,
            operation.key.is_empty(),
        ))
    }

    fn cut(&self, operation: &OperationNode, cutset: &CutSet) -> Result<&CutComputation> {
        self.cuts.get(cutset).ok_or_else(|| {
            eyre!("{operation} has no computed local counterterms for cut {cutset:?}")
        })
    }
}

impl Forests {
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
        let approx = FunctionBuilder::new(symbol!("T"));
        let frontier_atom = self.node_label_atom(frontier);

        let current = function!(
            symbol!(format!("S_{}", op.order.string_label())),
            usize::from(op.data)
        );
        let argument = if self.graph[frontier].covers().is_none() {
            current
        } else {
            let previous = Atom::var(symbol!(format!(
                "S_{}",
                self.graph[frontier]
                    .covers()
                    .expect("non-empty frontier cover must exist")
                    .string_label()
            )));
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

        let leaves = self.local_leaf_operations(node);
        let first_leaf = leaves
            .first()
            .ok_or_else(|| eyre!("{operation} has no leaf operations"))?;
        let mut full = self
            .compute_store
            .require(&first_leaf.frontier)?
            .recursion_input_4d(&first_leaf.frontier, settings.pole_part)
            .wrap_err_with(|| {
                format!(
                    "while loading 4D dependency frontier {} for {operation}",
                    first_leaf.frontier
                )
            })?;
        let mut local_4d = None;
        let mut integrated = None;
        let ctx = UVCtx::new(graph, settings);
        let integrated_approximation = Integrated::new(vakint, &self.wood.vakint_settings);

        for (step_order, leaf) in (first_leaf.frontier_depth()..).zip(leaves) {
            let (current, given) = self.wood.current_given_pair(leaf.op.data, step_order);
            let next_local = local_4d::uv_limit(&full, &ctx, &current, &given)?;
            let next_integrated = if settings.generate_integrated {
                integrated_approximation.run(&next_local, &ctx, &current, &given)?
            } else {
                IntegratedCts::root()
            };
            full = Full4DCts::recursion_input(
                &next_local,
                &next_integrated,
                settings.pole_part,
                false,
            );
            local_4d = Some(next_local);
            integrated = Some(next_integrated);
        }

        Ok((
            local_4d.expect("a non-root operation has at least one leaf"),
            integrated.expect("a non-root operation has at least one leaf"),
        ))
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
        } else {
            let leaves = self.local_leaf_operations(node);
            let first_leaf = leaves
                .first()
                .ok_or_else(|| eyre!("{operation} has no local leaf operations"))?;
            // An empty dependency frontier starts from the per-cut root integrand;
            // otherwise its typed local result is the sequential accumulator.
            let mut accumulator = if first_leaf.frontier.key.is_empty() {
                Local3DCts::root(graph, localizer)?
            } else {
                self.compute_store
                    .require(&first_leaf.frontier)?
                    .cut(&first_leaf.frontier, cutset)?
                    .local_3d
                    .clone()
            };
            for (step_order, leaf) in (first_leaf.frontier_depth()..).zip(leaves) {
                let frontier_integrated = self
                    .compute_store
                    .require(&leaf.frontier)?
                    .integrated(&leaf.frontier)
                    .wrap_err_with(|| {
                        format!(
                            "while loading integrated dependency frontier {} for {operation}",
                            leaf.frontier
                        )
                    })?;
                let (current, given) = self.wood.current_given_pair(leaf.op.data, step_order);
                // `run` applies the local and integrated subtraction signs; no external
                // sign or raw Foata-level product is introduced here.
                accumulator = Local3DApproximation::new(localizer, graph, settings).run(
                    &accumulator,
                    frontier_integrated,
                    &current,
                    &given,
                )?;
            }
            accumulator
        };

        let integrated = self
            .compute_store
            .require(operation)?
            .integrated(operation)?;
        let forest_node = operation.forest_node(graph, operation.key.op_count());
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

    // pub fn local_subtract(
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
        for (_, mut crown, key) in self.graph.iter_nodes() {
            if crown.any(|r| self.graph.flow(r).is_source()) {
                continue;
            }

            let forest_node = key.forest_node(graph, key.key.op_count());
            let atom = marker.prefix(
                &graph.full_filter(),
                forest_node.subgraph(),
                &self
                    .compute_store
                    .require(key)?
                    .integrated(key)?
                    .physical_atom(),
            );
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
        uv::{UltravioletGraph, Wood as OldWood},
    };

    use super::*;
    use color_eyre::Result;

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
        println!("{}", f);
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

    mod failing {
        use super::*;

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
                @"64",
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
                f.graph.n_nodes(),
                @"160");

            Ok(())
        }

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
                f.node_label(NodeIndex(10)),
                @"{C} · {36,F}: T((-1*S_C+S_F(11))*T(S_C(1)))*T(S_36(2))"
            );
            insta::assert_snapshot!(
                f.node_label(NodeIndex(11)),
                @"{3} · {36,F}: T((-1*S_3+S_F(4))*T(S_3(3)))*T(S_36(2))"
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
                f.node_label(NodeIndex(10)),
                @"{C} · {36,F}: T((-1*S_C+S_F(11))*T(S_C(1)))*T(S_36(2))"
            );
            insta::assert_snapshot!(
                f.node_label(NodeIndex(11)),
                @"{3} · {36,F}: T((-1*S_3+S_F(4))*T(S_3(3)))*T(S_36(2))"
            );

            Ok(())
        }
    }
}
