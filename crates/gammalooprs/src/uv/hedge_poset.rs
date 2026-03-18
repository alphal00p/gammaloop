use std::{collections::BTreeMap, fmt::Display};

use ahash::{AHashMap, AHashSet};
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::{
    HedgeGraph, NoData, NodeIndex,
    algorithms::trace_unfold::{HiddenData, Independence, TraceKey, TraceUnfold},
    involution::{EdgeIndex, Flow, HedgePair},
    nodestore::{NodeStorageOps, NodeStorageVec},
    subgraph::{Inclusion, ModifySubSet, SuBitGraph, SubSetLike, SubSetOps},
};
use spenso::network::library::TensorLibraryData;
use symbolica::{
    atom::{Atom, FunctionBuilder},
    function, symbol,
};
use vakint::Vakint;

use crate::{
    graph::{Graph, LMBext, LoopMomentumBasis, cuts::CutSet},
    uv::{
        UVgenerationSettings, UltravioletGraph,
        approx::{ApproximationKernel, ForestNodeLike, UVCtx, integrated::Integrated},
    },
};
use color_eyre::Result;

#[allow(clippy::derived_hash_with_manual_eq)]
#[derive(Clone, Debug, Hash, Eq)]
pub struct Spinney {
    pub subgraph: SuBitGraph,
    components: Vec<SuBitGraph>,
    pub dod: i32,
    pub lmb: LoopMomentumBasis,
    // pub topo_order: Option<usize>,
}

impl Spinney {
    pub fn compatible_with(&self, cut: &CutSet) -> bool {
        !self.subgraph.intersects(&cut.union)
    }

    pub fn empty<E, V, H, G: AsRef<HedgeGraph<E, V, H>> + LMBext>(g: &G) -> Self {
        Spinney {
            components: vec![],
            dod: 0,

            subgraph: g.as_ref().empty_subgraph(),
            lmb: g.empty_lmb(),
        }
    }

    pub fn forest_node<'a>(&'a self, topo_order: usize) -> ForestNode<'a> {
        ForestNode {
            spinney: self,
            topo_order,
        }
    }
    pub fn new<E, V, H, G: UltravioletGraph + AsRef<HedgeGraph<E, V, H>>>(
        subgraph: SuBitGraph,
        g: &G,
        lmb: &LoopMomentumBasis,
    ) -> Self {
        Spinney {
            components: g.as_ref().connected_components(&subgraph),
            dod: g.dod(&subgraph),
            lmb: g.compatible_sub_lmb(&subgraph, g.dummy_less_full_crown(&subgraph), lmb),
            subgraph,
        }
    }
}

impl PartialEq for Spinney {
    fn eq(&self, other: &Self) -> bool {
        self.subgraph == other.subgraph
    }
}

impl PartialOrd for Spinney {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self.subgraph.includes(&other.subgraph) {
            Some(std::cmp::Ordering::Greater)
        } else if other.subgraph.includes(&self.subgraph) {
            Some(std::cmp::Ordering::Less)
        } else {
            None
        }
    }
}

pub struct Wood {
    pub graph: HedgeGraph<SuBitGraph, Spinney>,
    pub root: NodeIndex,
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
    pub(crate) fn from_spinneys<E, V, H, I: IntoIterator<Item = Spinney>>(
        s: I,
        graph: impl AsRef<HedgeGraph<E, V, H>> + LMBext,
    ) -> Self {
        let mut spinneyset: AHashSet<_> = s.into_iter().collect();
        spinneyset.insert(Spinney::empty(&graph));

        let mut unions = AHashSet::new();
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
                let reduced_subgraph = sink_subgraph.subtract(source_subgraph);
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
                    e.map(|_| sink_subgraph.clone())
                }
            },
            |_, d| d,
        );

        let mut to_remove: SuBitGraph = poset.empty_subgraph();

        //Not quite transitive closure. We need to remove all edges that point to unions that are not single components
        for u in unions {
            // println!("//{u}:{}", poset[u].subgraph.string_label());
            let mut comps: AHashSet<_> = graph
                .as_ref()
                .connected_components(&poset[u].subgraph)
                .into_iter()
                .collect();
            for c in poset.iter_crown(u) {
                let Flow::Sink = poset.flow(c) else {
                    continue;
                };

                let Some(nid) = poset.involved_node_id(c) else {
                    continue;
                };
                if comps.contains(&poset[nid].subgraph) {
                    comps.remove(&poset[nid].subgraph);
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
        }
    }

    fn compatible_with(&self, cut: &CutSet) -> SuBitGraph {
        let mut compatible: SuBitGraph = self.graph.empty_subgraph();
        for (_, crown, s) in self.graph.iter_nodes() {
            if s.compatible_with(cut) {
                for h in crown {
                    compatible.add(h);
                }
            }
        }
        compatible
    }

    pub fn unfold(&self) -> Forests {
        Forests {
            graph: self.trace_unfold::<NodeStorageVec<_>>(self.root).map(
                |_, _, key| OperationNode { key },
                |_, _, _, _, e| e,
                |_, h| h,
            ),
            root: self.root,
            compute_store: AHashMap::new(),
        }
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

pub struct Forests {
    pub graph: HedgeGraph<EdgeIndex, OperationNode>,
    pub root: NodeIndex,
    pub compute_store: AHashMap<OperationNode, ComputeNode>,
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
    pub fn current<'a>(&'a self, wood: &'a Wood, topo_order: usize) -> Option<Vec<ForestNode<'a>>> {
        Some(
            self.key
                .levels
                .last()?
                .iter()
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
                .collect(),
        )
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
        self.spinney.subgraph.subtract(&given.spinney.subgraph)
    }
    fn subgraph(&self) -> &SuBitGraph {
        &self.spinney.subgraph
    }
    fn topo_order(&self) -> usize {
        self.topo_order
    }
}

impl Display for OperationNode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.key.levels.is_empty() {
            write!(f, "∅")
        } else {
            let mut acc: Option<SuBitGraph> = None;

            if f.alternate() {
                self.key.write_foata_like(f, |op| {
                    if let Some(a) = &mut acc {
                        a.union_with(op);
                    } else {
                        acc = Some(op.clone());
                    }
                    op.string_label()
                })?;
            }
            write!(f, "{}", self.to_atom())?;

            Ok(())
        }
    }
}

impl OperationNode {
    pub fn covers(&self) -> Option<SuBitGraph> {
        let mut acc: Option<SuBitGraph> = None;

        if self.key.levels.is_empty() {
            return acc;
        }

        for level in self.key.levels.iter() {
            for op in level.iter() {
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
        if self.key.levels.is_empty() {
            return acc;
        }

        let mut last = SuBitGraph::empty(self.key.levels[0][0].order.size());
        for l in &self.key.levels {
            let last_sym = if last.is_empty() {
                Atom::Zero
            } else {
                Atom::var(symbol!(format!("S_{}", last.string_label())))
            };

            let mut mul = Atom::one();

            for op in l {
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

    // fn simple_op()
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
    pub local_3d: Integrands, //3d denoms
    pub final_integrand: Integrands,
    pub integrated_4d: Integrand, //4d
    pub simple: Integrand,
}

impl Default for ComputeNode {
    fn default() -> Self {
        ComputeNode {
            local_3d: Integrands::NotComputed,
            final_integrand: Integrands::NotComputed,
            integrated_4d: Integrand::NotComputed,
            simple: Integrand::NotComputed,
        }
    }
}

impl Forests {
    pub fn iter_parents<'a>(
        &'a self,
        node: NodeIndex,
        order: usize,
        wood: &'a Wood,
    ) -> impl Iterator<
        Item = Result<(
            &'a ComputeNode,
            ForestNode<'a>,
            ForestNode<'a>,
            &'a OperationNode,
            bool,
        )>,
    > + 'a {
        let mut current = None;
        let mut is_union = false;
        self.graph
            .iter_crown(node)
            .filter(|h| self.graph.flow(*h).is_sink())
            .map(move |h| {
                // iterate over the sink-half-edges of the forest, i.e. the incoming half-edges to the current node
                // most of the time this will be a single half-edge, but in the case of a union, there may be multiple

                let wood_eid = self.graph[self.graph[&self.graph[&h]].0];

                let HedgePair::Paired { source, sink } = wood.graph[&wood_eid].1 else {
                    panic!("edge in wood is not paired");
                };

                // get this hedge's forest node from the wood. This is the node that has already been computed (as it is a parent to this edge)
                let given = wood.graph[wood.graph.node_id(source)].forest_node(order);

                // this is the current node, which should be the same for all union edges (since they all have the same sink)
                let current_for_h = wood.graph.node_id(sink);
                if let Some(current) = &current {
                    if current != &current_for_h {
                        return Err(eyre!("Mismatched current nodes"));
                    } else {
                        is_union = true;
                    }
                } else {
                    current = Some(current_for_h);
                }

                // this is the current node, which we want to compute with
                let current = wood.graph[current_for_h].forest_node(order);

                // this is the parent node, in the forest, which has already been computed and we want to get the computed value
                let parent_node = self.graph.node_id(h);
                let parent_key = &self.graph[parent_node];
                let computed = self
                    .compute_store
                    .get(parent_key)
                    .ok_or(eyre!("{} not yet added to store", parent_key))?;

                Ok((computed, current, given, parent_key, is_union))
            })
    }

    pub fn integrate(
        &mut self,
        graph: &mut Graph,
        wood: &Wood,
        vakint: (&Vakint, &vakint::VakintSettings),
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        let integrated_orchestrator = Integrated::new(vakint.0, vakint.1);
        let uvctx = UVCtx {
            graph: &*graph,
            settings,
        };
        for (order, nidx) in self.graph.topo_sort_kahn().unwrap().iter().enumerate() {
            let mut integrand = Atom::num(1);

            for h in self.iter_parents(*nidx, order, wood) {
                let (computed, current, given, parent_key, is_union) = h?;

                let Integrand::Single(a) = &computed.integrated_4d else {
                    return Err(eyre!("{} integrated_4d not computed yet", parent_key));
                };
                if is_union {
                    integrand *= a;
                } else {
                    integrand = integrated_orchestrator.kernel(&uvctx, &current, &given, a)?;
                }
            }

            self.compute_store
                .entry(self.graph[*nidx].clone())
                .or_insert(ComputeNode::default())
                .integrated_4d = Integrand::Single(integrand);
        }

        Ok(())
    }

    pub fn local_subtract(
        &mut self,
        graph: &mut Graph,
        wood: &Wood,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        Ok(())
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

                println!("Node {}: {}:{:#}", nidx, trace_key, trace_key);
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
        self.graph.dot_impl_fmt(
            f,
            &self.graph.full_filter(),
            "start=2;\n",
            &|_| None,
            &|_| None,
            &|v| Some(format!("label=\"{}\"", v)),
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

        let spinneys = dumbell.spinneys(&dumbell.full_filter());
        let f = Wood::from_spinneys(
            dumbell
                .spinneys(&dumbell.full_filter())
                .into_iter()
                .map(|a| Spinney::new(a.filter, &dumbell, &dumbell.loop_momentum_basis)),
            &dumbell,
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
                let spinneys = g.spinneys(&g.full_filter());
                let f = Wood::from_spinneys(
                    g.spinneys(&g.full_filter())
                        .into_iter()
                        .map(|a| Spinney::new(a.filter, &g, &g.loop_momentum_basis)),
                    &g,
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

        let f = Wood::from_spinneys(
            mercedes
                .spinneys(&mercedes.full_filter())
                .into_iter()
                .map(|a| Spinney::new(a.filter, &mercedes, &mercedes.loop_momentum_basis)),
            &mercedes,
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
        let f = Wood::from_spinneys(
            sunrise
                .spinneys(&sunrise.full_filter())
                .into_iter()
                .map(|a| Spinney::new(a.filter, &sunrise, &sunrise.loop_momentum_basis)),
            &sunrise,
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
        let f = Wood::from_spinneys(
            sunrise
                .spinneys(&sunrise.full_filter())
                .into_iter()
                .map(|a| Spinney::new(a.filter, &sunrise, &sunrise.loop_momentum_basis)),
            &sunrise,
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
        let f = Wood::from_spinneys(
            sunrise
                .spinneys(&sunrise.full_filter())
                .into_iter()
                .map(|a| Spinney::new(a.filter, &sunrise, &sunrise.loop_momentum_basis)),
            &sunrise,
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
        let f = Wood::from_spinneys(
            spectacles
                .spinneys(&spectacles.full_filter())
                .into_iter()
                .map(|a| Spinney::new(a.filter, &spectacles, &spectacles.loop_momentum_basis)),
            &spectacles,
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

        let f = Wood::from_spinneys(
            basketball
                .spinneys(&basketball.full_filter())
                .into_iter()
                .map(|a| Spinney::new(a.filter, &basketball, &basketball.loop_momentum_basis)),
            &basketball,
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
        let f = Wood::from_spinneys(
            fourloop_b
                .spinneys(&fourloop_b.full_filter())
                .into_iter()
                .map(|a| Spinney::new(a.filter, &fourloop_b, &fourloop_b.loop_momentum_basis)),
            &fourloop_b,
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

        let f = Wood::from_spinneys(
            four_loop_a
                .spinneys(&four_loop_a.full_filter())
                .into_iter()
                .map(|a| Spinney::new(a.filter, &four_loop_a, &four_loop_a.loop_momentum_basis)),
            &four_loop_a,
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
