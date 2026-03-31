use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Display,
};

use ahash::AHashMap;
use eyre::eyre;
use idenso::{color::ColorSimplifier, metric::MetricSimplifier};
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph, NoData, NodeIndex,
    algorithms::trace_unfold::{
        HiddenData, Independence, TraceKey, TraceUnfold, UnfoldedTraceGraph,
    },
    involution::{EdgeIndex, Flow, HedgePair},
    nodestore::{NodeStorageOps, NodeStorageVec},
    subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
};
use spenso::network::library::TensorLibraryData;
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    function, symbol,
};
use tracing::debug;
use vakint::Vakint;

use crate::{
    graph::{Graph, LMBext, LoopMomentumBasis, cuts::CutSet, parse::string_utils::ToOrderedSimple},
    utils::{W_, symbolica_ext::LogPrint},
    uv::{
        RenormalizationPart, Spinney, UVgenerationSettings, UltravioletGraph,
        approx::{
            ApproximationKernel, CutStructure, ForestNodeLike, UVCtx, integrated::Integrated,
            local_3d::Local3DApproximation,
        },
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

    pub(crate) fn new(cuts: CutStructure, graph: &Graph, vakint_settings: &VakintSettings) -> Self {
        let mut subgraph = graph.full_filter();
        subgraph.subtract_with(&graph.initial_state_cut.left);
        let mut spinneys = Vec::new();

        for cut in cuts.cuts.iter() {
            let cut_sub = subgraph.subtract(&cut.union);
            spinneys.extend(graph.spinneys(&cut_sub));
        }
        spinneys.sort_by(|a, b| a.filter.cmp(&b.filter));
        spinneys.dedup_by(|a, b| a.filter == b.filter);

        Self::from_spinneys(
            spinneys
                .into_iter()
                .map(|a| Spinney::new(a, graph, &graph.loop_momentum_basis)),
            graph,
            cuts,
            vakint_settings,
        )
    }

    pub(crate) fn new_with_settings(
        cuts: CutStructure,
        graph: &Graph,
        settings: &UVgenerationSettings,
    ) -> Self {
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
        let forests = Forests {
            graph,
            cuts,
            root: self.root,
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

impl ForestNodeLike for ForestNode<'_> {
    fn dod(&self) -> i32 {
        self.spinney.dod
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
        if settings.cached_integrated
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
                mul *= -integrated_orchestrator.kernel(&uvctx, &current, &given, &acc)?;
            }

            acc = mul
        }

        Ok(acc)
    }

    pub fn local(
        &self,
        graph: &Graph,
        cutset: &CutSet,
        compute_store: &mut ComputeStore,
        wood: &Wood,
        settings: &UVgenerationSettings,
    ) -> Result<Vec<Atom>> {
        let mut acc = None;
        let local = Local3DApproximation {};
        let uvctx = UVCtx { graph, settings };

        let mut order = 0;
        let mut levels = self.key.view();
        if settings.cached_integrated
            && let Some((prefix, leaf_level)) = self.key.split_last_level()
        {
            let prefix_key = OperationNode {
                key: prefix.to_owned(),
            };
            if let Some(computed) = compute_store.get(&prefix_key) {
                let Integrands::Multiple(cached) = &computed.local_3d[cutset] else {
                    return Err(eyre!("{} local_3d not computed yet", prefix_key));
                };
                acc = Some(cached.clone());
                order = prefix.op_count();
                levels = leaf_level;
            }
        }
        // let acc = acc.unwrap_or(Local3DApproximation::root(graph, cutset)?);

        for l in levels.iter_levels_top_down() {
            let mut mul = Atom::one();

            for op in l.iter_leaf_ops() {
                let (current, given) = wood.current_given_pair(op.data, order);
                order += 1;
                compute_store.record_kernel_hit();
                // mul *= -local.kernel(&uvctx, &current, &given, &acc)?;
            }

            // acc = mul
        }

        Ok(vec![])
        // Ok(acc)
    }
}

pub enum Integrand {
    NotComputed,
    Single(Atom),
}
pub enum Integrands {
    NotComputed,
    Multiple(Vec<Atom>),
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

impl Forests {
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

    pub fn integrate(
        &mut self,
        graph: &Graph,
        vakint: &Vakint,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        for (order, nidx) in self.graph.topo_sort_kahn().unwrap().iter().enumerate() {
            debug!(order=%order,cache=%settings.cached_integrated,nidx=%nidx,key=%self.graph[*nidx],"One integrated step");
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

    pub fn local_subtract(
        &mut self,
        graph: &mut Graph,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        for (cut_compatible_forest_subset, cuts) in &self.cuts {
            let mut first = true;

            for (order, nidx) in self
                .graph
                .topo_sort_kahn_of(cut_compatible_forest_subset)
                .unwrap()
                .iter()
                .enumerate()
            {
                debug!(order=%order,cache=%settings.cached_integrated,nidx=%nidx,key=%self.graph[*nidx],"One integrated step");
                let integrands = if first {
                    first = false;
                    Local3DApproximation::root(graph, cuts)
                } else {
                    self.graph[*nidx].local(
                        graph,
                        cuts,
                        &mut self.compute_store,
                        &self.wood,
                        settings,
                    )
                }?;

                self.compute_store
                    .entry(self.graph[*nidx].clone())
                    .or_default()
                    .local_3d
                    .insert(cuts.clone(), Integrands::Multiple(integrands));
            }
        }

        Ok(())
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

        Ok(RenormalizationPart::new(
            sum.replace_multiple(&replacements),
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
            .map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis))
            .collect();
        let f = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &VakintSettings::default(),
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
            .map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis))
            .collect();
        let f = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &VakintSettings::default(),
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
            .map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis))
            .collect();
        let f = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &VakintSettings::default(),
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
            &VakintSettings::default(),
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
            .map(|a| Spinney::new(a, &dumbell, &dumbell.loop_momentum_basis))
            .collect();
        let f = Wood::new(
            CutStructure::empty(&dumbell),
            &dumbell,
            &VakintSettings::default(),
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
                    .map(|a| Spinney::new(a, &g, &g.loop_momentum_basis))
                    .collect();
                let f = Wood::new(CutStructure::empty(&g), &g, &VakintSettings::default());

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
            &VakintSettings::default(),
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
            &VakintSettings::default(),
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
            &VakintSettings::default(),
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
            &VakintSettings::default(),
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
            &VakintSettings::default(),
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
            &VakintSettings::default(),
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
            &VakintSettings::default(),
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
            &VakintSettings::default(),
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
}
