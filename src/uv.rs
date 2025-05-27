use std::{
    cmp::Ordering,
    collections::VecDeque,
    hash::Hash,
    ops::{Deref, Index},
    sync::LazyLock,
};

use crate::{
    cff::{
        expression::OrientationData,
        generation::{generate_uv_cff, ShiftRewrite},
    },
    cross_section::SuperGraph,
    graph::{Graph, VertexInfo},
    model::ArcParticle,
    momentum::Sign,
    new_graph::{self, no_filter, Edge, LMBext, LoopMomentumBasis, Vertex},
    utils::{GS, W_},
};
use ahash::AHashSet;
use bitvec::vec::BitVec;
use eyre::eyre;
use itertools::Itertools;
use pathfinding::prelude::BfsReachable;
use serde::{Deserialize, Serialize};
use spenso::{
    network::parsing::ShadowedStructure,
    shadowing::symbolica_utils::SerializableAtom,
    structure::{
        abstract_index::AbstractIndex,
        dimension::Dimension,
        representation::{Minkowski, RepName},
        HasStructure, NamedStructure, OrderedStructure, ToSymbolic,
    },
    tensors::parametric::{atomcore::PatternReplacement, ParamTensor},
};
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder, Symbol},
    domains::integer::Integer,
    function,
    id::Replacement,
    parse,
    printer::PrintOptions,
    symbol,
};
use symbolica_community::physics::algebraic_simplification::metric::{MetricSimplifier, MS};

use linnet::half_edge::subgraph::{InternalSubGraph, SubGraph, SubGraphOps};
use linnet::half_edge::{
    involution::{EdgeIndex, Hedge, HedgePair, SignOrZero},
    subgraph::{Cycle, ModifySubgraph},
    HedgeGraph,
};
use typed_index_collections::TiVec;
// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

use crate::{
    graph::{BareEdge, BareGraph, BareVertex},
    model::normalise_complex,
};

// pub static VAKINT: LazyLock<Vakint> = LazyLock::new(|| {
//     Vakint::new(Some(VakintSettings {
//         evaluation_order: EvaluationOrder::alphaloop_only(),
//         integral_normalization_factor: LoopNormalizationFactor::MSbar,
//         run_time_decimal_precision: 16,
//         ..VakintSettings::default()
//     }))
//     .unwrap()
// });

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct UVEdge {
    og_edge: usize,
    dod: i32,
    particle: ArcParticle,
    // prop:ArcPropagator.
    num: Atom,
    den: Atom,
}

impl UVEdge {
    pub fn from_edge(edge: &BareEdge, id: usize, bare_graph: &BareGraph) -> Self {
        let [colorless, _] = edge.color_separated_numerator(bare_graph, id);
        UVEdge {
            particle: edge.particle.clone(),
            og_edge: id,
            dod: edge.dod() as i32,
            num: normalise_complex(&colorless).into(),
            den: Atom::num(1).into(), // edge.full_den(bare_graph, index).into(),
        }
    }

    pub fn numerator_rank(&self) -> isize {
        0
    }
}

#[derive(Clone, Debug)]
pub struct UVNode {
    dod: i32,
    num: Atom,
    color: Option<ParamTensor>,
}

impl UVNode {
    pub fn from_vertex(vertex: &BareVertex, graph: &BareGraph) -> Self {
        if let Some((colorless, color)) = vertex.contracted_colorless_vertex_rule(graph) {
            UVNode {
                dod: vertex.dod() as i32,
                num: normalise_complex(&colorless).into(),
                color: Some(ParamTensor::param(
                    color.map_structure(OrderedStructure::from),
                )),
            }
        } else {
            UVNode {
                dod: vertex.dod() as i32,
                num: Atom::num(1).into(),
                color: None,
            }
        }
    }

    pub fn numerator_rank(&self) -> isize {
        0
    }
}

#[derive(Clone, Debug)]
pub struct UVGraph {
    hedge_graph: HedgeGraph<UVEdge, UVNode>,
    lmb: LoopMomentumBasis,
    lmb_replacement: Vec<Replacement>,
}

impl Deref for UVGraph {
    type Target = HedgeGraph<UVEdge, UVNode>;
    fn deref(&self) -> &Self::Target {
        &self.hedge_graph
    }
}

impl AsRef<HedgeGraph<UVEdge, UVNode>> for UVGraph {
    fn as_ref(&self) -> &HedgeGraph<UVEdge, UVNode> {
        &self.hedge_graph
    }
}
pub fn spenso_lor(
    tag: i32,
    ind: impl Into<AbstractIndex>,
    dim: impl Into<Dimension>,
) -> ShadowedStructure {
    let mink = Minkowski {}.new_slot(dim, ind);
    NamedStructure::from_iter([mink], GS.emr_mom, Some(vec![Atom::num(tag)])).structure
}

pub fn spenso_lor_atom(tag: i32, ind: impl Into<AbstractIndex>, dim: impl Into<Dimension>) -> Atom {
    spenso_lor(tag, ind, dim).to_symbolic().unwrap()
}

#[allow(dead_code)]
impl UVGraph {
    pub fn from_supergraph(sg: &new_graph::Graph) -> Self {
        let hedge_graph = sg.underlying.map_data_ref(
            |_, _, n| UVNode {
                dod: n.dod,
                num: n.num.clone(),
                color: None,
            },
            |_, eid, _, n| {
                n.map(|d| {
                    let m2 = parse!(d.particle.mass.name).npow(2);
                    UVEdge {
                        og_edge: 1,
                        dod: d.dod,
                        particle: d.particle.clone(),
                        num: d.num.clone(),
                        den: spenso_lor_atom(usize::from(eid) as i32, 1, GS.dim)
                            .npow(2)
                            //.to_dots()
                            - m2,
                    }
                })
            },
        );

        let reps = hedge_graph.normal_emr_replacement(
            &hedge_graph.full_filter(),
            &sg.loop_momentum_basis,
            &[W_.x___],
            no_filter,
        );
        UVGraph {
            hedge_graph,
            lmb: sg.loop_momentum_basis.clone(),
            lmb_replacement: reps,
        }
    }

    pub fn from_underlying(hedge_graph: &HedgeGraph<Edge, Vertex>) -> Self {
        Self::from_hedge(hedge_graph.map_data_ref(
            |_, _, n| UVNode {
                dod: n.dod,
                num: n.num.clone(),
                color: None,
            },
            |_, eid, _, n| {
                n.map(|d| {
                    let m2 = parse!(d.particle.mass.name).npow(2);
                    UVEdge {
                        og_edge: 1,
                        dod: d.dod,
                        particle: d.particle.clone(),
                        num: d.num.clone(),
                        den: spenso_lor_atom(usize::from(eid) as i32, 1, GS.dim)
                            .npow(2)
                            //.to_dots()
                            - m2,
                    }
                })
            },
        ))
    }

    pub fn from_hedge(hedge_graph: HedgeGraph<UVEdge, UVNode>) -> Self {
        let full = hedge_graph.full_filter();

        let lmb = hedge_graph.lmb(&full);

        // println!("//lmb{lmb:?} for \n{}", hedge_graph.dot_lmb(&full, &lmb));
        let reps = hedge_graph.normal_emr_replacement(&full, &lmb, &[W_.x___], no_filter);
        UVGraph {
            hedge_graph,
            lmb,
            lmb_replacement: reps,
        }
    }

    pub fn from_graph(graph: &BareGraph) -> Self {
        let mut excised: BitVec = graph.hedge_representation.empty_subgraph();

        for (_, n, d) in graph.hedge_representation.iter_nodes() {
            if matches!(
                graph.vertices[*d].vertex_info,
                VertexInfo::ExternalVertexInfo(_)
            ) {
                excised.union_with(&n.into())
            }
        }

        excised = excised.complement(&graph.hedge_representation);

        let excised = graph.hedge_representation.concretize(&excised).map(
            |_, _, d| UVNode::from_vertex(&graph.vertices[*d], graph),
            |_, _, _, e| e.map(|d| UVEdge::from_edge(&graph.edges[*d], *d, graph)),
        );

        UVGraph::from_hedge(excised)
    }

    fn n_loops<S: SubGraph>(&self, subgraph: &S) -> usize {
        // if let Some(loop_count) = subgraph.loopcount {
        //     // println!("found loop_count nloops: {}", loop_count);
        //     loop_count
        // } else {
        self.cyclotomatic_number(subgraph)
        // }
    }

    fn numerator<S: SubGraph>(&self, subgraph: &S) -> Atom {
        let mut num = Atom::num(1);

        for (_, _, n) in self.iter_nodes_of(subgraph) {
            num = num * &n.num;
        }

        for (pair, _eid, d) in self.iter_edges_of(subgraph) {
            if matches!(pair, HedgePair::Paired { .. }) {
                num = num * &d.data.num;
            }
        }

        num.into()
    }

    fn oriented_numerator<S: SubGraph>(
        &self,
        subgraph: &S,
        // orientationdata: &OrientationData,
    ) -> Atom {
        let mut num = Atom::num(1);

        for (_, _, n) in self.iter_nodes_of(subgraph) {
            num = num * &n.num;
        }

        for (pair, _eid, d) in self.iter_edges_of(subgraph) {
            if matches!(pair, HedgePair::Paired { .. }) {
                num = num * &d.data.num;
            }
        }

        num
    }

    fn denominator<S: SubGraph>(&self, subgraph: &S) -> Atom {
        let mut den = Atom::num(1);

        for (pair, eid, d) in self.iter_edges_of(subgraph) {
            if matches!(pair, HedgePair::Paired { .. }) {
                let m2 = parse!(d.data.particle.mass.name).npow(2);
                den = den
                    * function!(
                        GS.den,
                        usize::from(eid) as i64,
                        function!(GS.emr_mom, usize::from(eid) as i64),
                        m2,
                        spenso_lor_atom(usize::from(eid) as i32, 1, GS.dim)
                            .npow(2)
                            //.to_dots()
                            - m2
                    );
            }
        }

        den.into()
    }

    fn dot<S: SubGraph>(&self, subgraph: &S) -> String {
        self.dot_impl(
            subgraph,
            format!("dod_of_subgraph={}", self.dod(subgraph)),
            &|e| Some(format!("label=\"{}\"", e.dod)),
            &|n| Some(format!("label=\"{}\"", n.dod)),
        )
    }

    fn t_op<I: Iterator<Item = InternalSubGraph>>(&self, mut subgraph_iter: I) -> Atom {
        if let Some(subgraph) = subgraph_iter.next() {
            let t = self.t_op(subgraph_iter);
            FunctionBuilder::new(symbol!("Top"))
                .add_arg(&Atom::num(self.dod(&subgraph)))
                .add_arg(&t)
                .finish()
        } else {
            Atom::num(1)
        }
    }
}

pub struct Wood {
    poset: Poset<InternalSubGraph, ()>,
    additional_unions: SecondaryMap<PosetNode, Vec<PosetNode>>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct IntegrandExpr {
    integrand: Atom,
    // add_arg: Option<Atom>,
}

impl IntegrandExpr {
    pub fn from_subgraph<S: SubGraph>(subgraph: &S, graph: &UVGraph) -> Self {
        let num = graph.numerator(subgraph);

        let den = graph.denominator(subgraph);

        IntegrandExpr {
            integrand: num / den,
        }
    }
}

// pub fn limit(&)

pub trait UltravioletGraph: LMBext {
    fn all_cycle_unions<E, V, S: SubGraph<Base = BitVec>>(
        &self,
        subgraph: &S,
    ) -> AHashSet<InternalSubGraph>
    where
        Self: AsRef<HedgeGraph<E, V>>,
    {
        let ref_graph = self.as_ref();
        let init_node = ref_graph.iter_nodes_of(subgraph).next().unwrap().0;
        let all_subcycles: Vec<_> = Cycle::all_sum_powerset_filter_map(
            &ref_graph
                .paton_cycle_basis(subgraph, &init_node, None)
                .unwrap()
                .0,
            &Some,
        )
        .map(|a| a.into_iter().map(|c| c.internal_graph(ref_graph)).collect())
        .unwrap();

        // println!("{}", self.base_dot());
        let spinneys: AHashSet<_> = InternalSubGraph::all_unions_iterative(&all_subcycles);

        spinneys
    }
    fn limit<E, V, S: SubGraph>(&self, subgraph: &S, expr: &Atom, expansion: Symbol) -> Atom
    where
        Self: AsRef<HedgeGraph<E, V>>,
    {
        let mut expr = expr.clone();
        for (p, eid, _) in self.as_ref().iter_edges_of(subgraph) {
            if matches!(p, HedgePair::Paired { .. }) {
                expr = expr
                    .replace(function!(GS.emr_mom, usize::from(eid) as i32, W_.x___))
                    .with(function!(GS.emr_mom, usize::from(eid) as i32, W_.x___) * expansion);
            }
        }
        expr
    }

    fn wood<E, V, S: SubGraph<Base = BitVec>>(&self, subgraph: &S) -> Wood
    where
        Self: AsRef<HedgeGraph<E, V>>,
    {
        Wood::from_spinneys(self.spinneys(subgraph), self)
    }

    fn dod<S: SubGraph>(&self, subgraph: &S) -> i32;

    fn spinneys<E, V, S: SubGraph<Base = BitVec>>(&self, subgraph: &S) -> AHashSet<InternalSubGraph>
    where
        Self: AsRef<HedgeGraph<E, V>>,
    {
        let ref_graph = self.as_ref();
        let init_node = ref_graph.iter_nodes_of(subgraph).next().unwrap().0;
        let all_subcycles: Vec<_> = Cycle::all_sum_powerset_filter_map(
            &ref_graph
                .paton_cycle_basis(subgraph, &init_node, None)
                .unwrap()
                .0,
            &Some,
        )
        .map(|a| a.into_iter().map(|c| c.internal_graph(ref_graph)).collect())
        .unwrap();

        // println!("{}", self.base_dot());
        let mut spinneys: AHashSet<_> = InternalSubGraph::all_ops_iterative_filter_map(
            &all_subcycles,
            &|a, b| a.union(b),
            &|union| {
                if self.dod(&union) >= 0 {
                    Some(union)
                } else {
                    None
                }
            },
        );

        spinneys.insert(ref_graph.empty_subgraph());

        spinneys
    }
}

impl LMBext for UVGraph {
    fn empty_lmb(&self) -> LoopMomentumBasis {
        self.as_ref().empty_lmb()
    }

    fn dot_lmb<S: SubGraph>(&self, subgraph: &S, lmb: &LoopMomentumBasis) -> String {
        self.as_ref().dot_lmb(subgraph, lmb)
    }
    fn lmb<S: SubGraph<Base = BitVec>>(&self, subgraph: &S) -> LoopMomentumBasis {
        self.as_ref().lmb(subgraph)
    }

    fn lmb_impl<S: SubGraph + SubGraphOps + ModifySubgraph<HedgePair> + ModifySubgraph<Hedge>>(
        &self,
        subgraph: &S,
        tree: &S,
        externals: &S,
    ) -> LoopMomentumBasis {
        self.as_ref().lmb_impl(subgraph, tree, externals)
    }

    fn compatible_sub_lmb<S: SubGraph>(
        &self,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
    ) -> LoopMomentumBasis
    where
        S::Base: SubGraph<Base = S::Base>
            + SubGraphOps
            + Clone
            + ModifySubgraph<HedgePair>
            + ModifySubgraph<Hedge>,
    {
        self.as_ref().compatible_sub_lmb(subgraph, lmb)
    }

    fn replacement_impl<'a, S: SubGraph, I>(
        &self,
        rep: impl Fn(EdgeIndex, Atom, Atom) -> Replacement,
        subgraph: &S,
        lmb: &LoopMomentumBasis,
        loop_symbol: Symbol,
        ext_symbol: Symbol,
        loop_args: &'a [I],
        ext_args: &'a [I],
        filter_pair: fn(&HedgePair) -> bool,
        emr_id: bool,
    ) -> Vec<Replacement>
    where
        &'a I: Into<symbolica::atom::AtomOrView<'a>>,
    {
        self.as_ref().replacement_impl(
            rep,
            subgraph,
            lmb,
            loop_symbol,
            ext_symbol,
            loop_args,
            ext_args,
            filter_pair,
            emr_id,
        )
    }

    fn generate_loop_momentum_bases<S: SubGraph>(
        &self,
        subgraph: &S,
    ) -> TiVec<crate::new_graph::LmbIndex, LoopMomentumBasis>
    where
        S::Base: SubGraph<Base = S::Base>
            + SubGraphOps
            + Clone
            + ModifySubgraph<HedgePair>
            + ModifySubgraph<Hedge>,
    {
        self.as_ref().generate_loop_momentum_bases(subgraph)
    }
}

impl UltravioletGraph for UVGraph {
    fn dod<S: SubGraph>(&self, subgraph: &S) -> i32 {
        let mut dod: i32 = 4 * self.n_loops(subgraph) as i32;
        // println!("nloops: {}", dod / 4);

        for (p, _, e) in self.iter_edges_of(subgraph) {
            if p.is_paired() {
                dod += e.data.dod;
            }
        }

        for (_, _, n) in self.iter_nodes_of(subgraph) {
            dod += n.dod;
        }

        dod
    }
}

impl Wood {
    pub fn n_spinneys(&self) -> usize {
        self.poset.n_nodes()
    }

    pub fn from_spinneys<E, V, I: IntoIterator<Item = InternalSubGraph>>(
        s: I,
        graph: impl AsRef<HedgeGraph<E, V>>,
    ) -> Self {
        let mut poset = Poset::from_iter(s.into_iter().map(|s| (s, ())));
        let ref_graph = graph.as_ref();

        poset.invert();
        poset.compute_topological_order();

        let mut unions = SecondaryMap::new();

        for (i, sg) in poset.nodes.iter() {
            let cs = ref_graph.connected_components(&sg.data);

            if cs.len() > 1 {
                // sg is a disjoint union of spinneys (at the level of half-edges) (strongly disjoint)
                let mut union = vec![];

                for &c in sg.parents.iter() {
                    let mut is_in = 0;
                    for comp in &cs {
                        let comp =
                            InternalSubGraph::cleaned_filter_optimist(comp.clone(), ref_graph);
                        if comp == poset.nodes[c].data {
                            // find the components in the wood that this union is made of
                            union.push(c);
                            is_in += 1;
                        }
                    }

                    if is_in > 1 {
                        panic!("is in too many components")
                    }
                }

                unions.insert(i, union.clone());
                // sg.children = union;
            }
        }

        // let coverset = poset.to_cover_set();
        Wood {
            poset,
            additional_unions: unions,
        }
    }

    fn unfold_bfs<E, V, G>(
        &self,
        graph: &G,
        lmb: &LoopMomentumBasis,
        dag: &mut DAG<Approximation, DagNode, ()>,
        unions: &mut SecondaryMap<PosetNode, Option<Vec<(PosetNode, Option<DagNode>)>>>,
        root: PosetNode,
    ) -> DagNode
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V>>,
    {
        // let graph = graph.as_ref();
        let mut search_front = VecDeque::new();

        let tree_root = dag.add_node(Approximation::new(
            self.poset.data(root).clone(),
            graph,
            lmb,
        ));
        search_front.push_front((root, tree_root));

        while let Some((node, parent)) = search_front.pop_front() {
            for c in &self.poset.nodes[node].children {
                if let Some(tagged_union) = unions.get_mut(*c) {
                    // Is this node a disjoint union of spinneys
                    if let Some(mut union) = tagged_union.take() {
                        let mut all_supplied = true;
                        for (p, d) in &mut union {
                            if self.poset.nodes[*p].data == self.poset.nodes[node].data {
                                *d = Some(parent);
                            }
                            if d.is_none() {
                                all_supplied = false;
                            }
                        }
                        if all_supplied {
                            let child = dag.add_node(Approximation::new(
                                self.poset.data(*c).clone(),
                                graph,
                                lmb,
                            ));
                            for (_, d) in union {
                                dag.add_edge(d.unwrap(), child);
                            }
                            search_front.push_front((*c, child));
                        } else {
                            *tagged_union = Some(union);
                        }
                    }
                } else {
                    let child =
                        dag.add_node(Approximation::new(self.poset.data(*c).clone(), graph, lmb));
                    dag.add_edge(parent, child);
                    search_front.push_front((*c, child));
                }
            }
        }
        tree_root
    }

    pub fn unfold<E, V, G>(&self, graph: &G, lmb: &LoopMomentumBasis) -> Forest
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V>>,
    {
        let mut dag: DAG<Approximation, DagNode, ()> = DAG::new();

        let root = self.poset.minimum().unwrap();

        let mut unions = SecondaryMap::new();

        for (p, u) in self.additional_unions.iter() {
            let union: Vec<(PosetNode, Option<DagNode>)> = u.iter().map(|i| (*i, None)).collect();
            unions.insert(p, Some(union));
        }

        let _ = self.unfold_bfs(graph, lmb, &mut dag, &mut unions, root);

        Forest { dag }
    }

    pub fn dot(&self, graph: &UVGraph) -> String {
        let shift = self.poset.shift();
        self.poset.to_dot_impl(&|n| {
            format!(
                "label={}, dod={}, n_edges = {},topo_order = {}",
                n.dot_id(shift),
                graph.dod(&n.data),
                graph.count_internal_edges(&n.data),
                n.order.unwrap()
            )
        })
    }

    pub fn dot_spinneys(&self, graph: &UVGraph) {
        for s in self.poset.node_values() {
            println!(
                "found {} loop spinney with dod {}:{} ",
                graph.n_loops(s),
                graph.dod(s),
                graph.dot(&s.to_hairy_subgraph(&graph))
            );
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum ApproxOp {
    NotComputed,
    Union {
        sign: Sign,
        t_args: Vec<IntegrandExpr>,
        subgraphs: Vec<InternalSubGraph>,
    },
    Dependent {
        sign: Sign,
        t_arg: IntegrandExpr,
        subgraph: InternalSubGraph,
    },
    Root,
}

pub struct SimpleApprox {
    t_args: Vec<Atom>,
    pub sign: Sign,
    graph: InternalSubGraph,
}

impl SimpleApprox {
    fn subgraph_shadow(graph: &BitVec, subgraph: &InternalSubGraph) -> Symbol {
        symbol!(&format!(
            "S_{}⊛{}",
            graph.string_label(),
            subgraph.string_label()
        ))
    }

    pub fn expr(&self, bigger_graph: &BitVec) -> Atom {
        let reduced = Atom::var(Self::subgraph_shadow(bigger_graph, &self.graph));
        let mut mul = Atom::num(1);
        for i in &self.t_args {
            mul = mul * i
        }
        reduced * mul
    }

    pub fn t_op(&self, bigger_graph: &BitVec) -> Atom {
        function!(GS.top, self.expr(bigger_graph))
    }

    pub fn root(subgraph: InternalSubGraph) -> Self {
        if !subgraph.is_empty() {
            panic!(
                "Root approximation must be empty {} {:?}",
                subgraph.string_label(),
                subgraph
            )
        }
        SimpleApprox {
            t_args: vec![],
            sign: Sign::Positive,
            graph: subgraph,
        }
    }

    pub fn dependent(&self, bigger_graph: InternalSubGraph) -> Self {
        Self {
            t_args: vec![self.t_op(&bigger_graph.filter)],
            sign: -self.sign,
            graph: bigger_graph,
        }
    }

    pub fn union<'a>(
        subgraph: InternalSubGraph,
        union: impl IntoIterator<Item = &'a Self>,
    ) -> Self {
        let mut t_args = vec![];
        let mut sign = Sign::Positive;
        for u in union {
            if u.t_args.len() != 1 {
                panic!("Union can only be applied to dependent approximations");
            }
            t_args.push(u.t_args[0].clone());
            sign = sign * u.sign;
        }
        SimpleApprox {
            t_args,
            sign,
            graph: subgraph,
        }
    }
}

pub struct Approximation {
    // The union of all spinneys, remaining graph is full graph minus subgraph
    subgraph: InternalSubGraph,
    dod: i32,
    lmb: LoopMomentumBasis,
    pub cff_expr: ApproxOp, //3d denoms
    pub t_op: ApproxOp,     //4d
    pub simple_approx: Option<SimpleApprox>,
}

impl Approximation {
    pub fn simplify_notation(expr: &Atom) -> Atom {
        let replacements = [(function!(GS.den, W_.a_, W_.x_), Atom::var(W_.x_))];

        let reps: Vec<_> = replacements
            .into_iter()
            .map(|(a, b)| Replacement::new(a.to_pattern(), b))
            .collect();

        expr.replace_multiple_repeat(&reps)
    }

    pub fn new<G, E, V>(
        spinney: InternalSubGraph,
        graph: &G,
        lmb: &LoopMomentumBasis,
    ) -> Approximation
    where
        G: UltravioletGraph + AsRef<HedgeGraph<E, V>>,
    {
        let lmb = graph.compatible_sub_lmb(&spinney, lmb);
        println!("//lmb for spinney \n{}", graph.dot_lmb(&spinney, &lmb));
        Approximation {
            dod: graph.dod(&spinney),
            subgraph: spinney,
            lmb,
            simple_approx: None,
            cff_expr: ApproxOp::NotComputed,
            t_op: ApproxOp::NotComputed,
        }
    }

    pub fn reduced_graph(&self, subgraph: &InternalSubGraph) -> InternalSubGraph {
        self.subgraph.subtract(subgraph)
    }

    pub fn structure_approximate(&self, dependent: &Self) -> SimpleApprox {
        dependent
            .simple_approx
            .as_ref()
            .unwrap()
            .dependent(self.subgraph.clone())
    }

    pub fn approximate(&self, dependent: &Self, graph: &UVGraph) -> ApproxOp {
        let reduced = self.subgraph.subtract(&dependent.subgraph);

        if let Some((inner_t, sign)) = dependent.t_op.expr() {
            let t_arg = IntegrandExpr::from_subgraph(&reduced, graph);

            let mut atomarg = t_arg.integrand * inner_t;

            // only apply replacements for edges in the reduced graph
            let mom_reps = graph.uv_wrapped_replacement(&reduced, &self.lmb, &[W_.x___]);

            println!("Reps:");
            for r in &mom_reps {
                println!("{r}");
            }

            println!(
                "Expand-prerep {} with dod={} in {:?}",
                atomarg, self.dod, self.lmb.ext_edges
            );

            // rewrite the inner_t as well
            atomarg = atomarg.replace_multiple(&mom_reps);

            println!(
                "Expand {} with dod={} in {:?}",
                atomarg, self.dod, self.lmb.ext_edges
            );
            for e in &self.lmb.ext_edges {
                atomarg = atomarg
                    .replace(function!(GS.emr_mom, usize::from(*e) as i64, W_.x___))
                    .with(function!(GS.emr_mom, usize::from(*e) as i64, W_.x___) * GS.rescale);
            }

            let soft_ct = graph.full_crown(&self.subgraph).count_ones() == 2 && self.dod > 0;

            let mut masses = AHashSet::new();
            masses.insert(Atom::var(GS.m_uv));
            // scale all masses, including UV masses from subgraphs

            for (p, _, e) in graph.iter_edges_of(&self.subgraph) {
                if p.is_paired() {
                    let e_mass = parse!(&e.data.particle.mass.name);
                    masses.insert(e_mass);
                }
            }

            if !soft_ct {
                for m in &masses {
                    let rescaled = m.clone() * GS.rescale;
                    atomarg = atomarg.replace(m.clone()).with(rescaled);
                }

                // expand the propagator around a propagator with a UV mass
                atomarg = atomarg
                    .replace(parse!("den(n_,q_,mass_,prop_)"))
                    .with(parse!(
                        "den(n_,q_,mass_ + mUV^2 - t^2*mUV^2, prop_- mUV^2 + t^2*mUV^2)"
                    ));
            }

            atomarg = atomarg
                .replace(function!(MS.dot, GS.rescale * W_.x_, W_.y_))
                .repeat()
                .with(function!(MS.dot, W_.x_, W_.y_) * GS.rescale);

            // println!("atomarg:{}", atomarg);

            // den(..) tags a propagator, its first derivative is 1 and the rest is 0
            let mut a = atomarg
                .series(GS.rescale, Atom::Zero, self.dod.into(), true)
                .unwrap()
                .to_atom()
                .replace(parse!("der(0,0,0,1, den(y__))"))
                .with(Atom::num(1))
                .replace(parse!("der(x__, den(y__))"))
                .with(Atom::num(0));

            if soft_ct {
                let coeffs = a.coefficient_list::<u8>(&[Atom::var(GS.rescale)]);
                let mut b = Atom::Zero;
                let dod_pow = Atom::var(GS.rescale).npow(self.dod);
                for (pow, mut i) in coeffs {
                    if pow == dod_pow {
                        // set the masses in the t=dod term to 0
                        // UV rearrange the denominators
                        for m in &masses {
                            i = i.replace(m.clone()).with(Atom::Zero);
                        }

                        i = i
                            .replace(parse!("den(n_,q_,mass_,prop_)"))
                            .with(parse!("den(n_,q_,mUV^2,prop_-mUV^2)"));
                    }

                    b += i;
                }

                a = b;
            } else {
                a = a.replace(GS.rescale).with(Atom::num(1));
            }

            // strip the momentum wrapper from the denominator
            a = a
                .replace(function!(
                    GS.den,
                    W_.prop_,
                    function!(GS.emr_mom, W_.prop_, W_.mom_),
                    W_.x__
                ))
                .with(function!(GS.den, W_.prop_, W_.mom_, W_.x__));

            println!("Expanded: {:>}", a.expand());

            /*let mut integrand_vakint = a.expand();
            let mut propagator_id = 1;
            for (pair, index, _data) in graph.iter_edges_of(&reduced) {
                // FIXME: not a good way to check for internal edges?
                if let HedgePair::Paired { source, sink } = pair {
                    integrand_vakint *= vakint_parse!(format!(
                        "prop({},edge({},{}),{},mUV,0)",
                        propagator_id,
                        usize::from(graph.node_id(source)),
                        usize::from(graph.node_id(sink)),
                        function!(GS.emr_mom, usize::from(index) as i64)
                            .replace_multiple(mom_reps)
                            .to_plain_string()
                    ))
                    .unwrap();
                    propagator_id += 1;
                }
            }

            integrand_vakint = integrand_vakint
                .replace(function!(GS.emr_mom, GS.x__))
                .with(function!(vakint_symbol!("k"), GS.x__))
                .replace(function!(GS.loop_mom, GS.x__))
                .with(function!(vakint_symbol!("k"), GS.x__))
                .expand()
                .replace(parse!("vk::prop(x_,y___,p_)*den(x_,z___)^n_").unwrap())
                .repeat()
                .with(parse!("vk::prop(x_,y___,p_-n_)").unwrap())
                .replace(parse!("den(x___)").unwrap())
                .with(Atom::num(1));

            let mut dummy_index = 1;
            while let Some(r) = integrand_vakint
                .replace(parse!("symbolica_community::dot(vk::k(x_),vk::k(y_))").unwrap())
                .iter(parse!(format!("vk::k(x_,{0})*vk::k(x_,{0})", dummy_index)).unwrap())
                .next()
            {
                integrand_vakint = r;
                dummy_index += 1;
            }

            integrand_vakint = integrand_vakint
                .replace(parse!("vk::prop(x__)").unwrap())
                .with(parse!("vk::topo(vk::prop(x__))").unwrap())
                .replace(parse!("vk::topo(x_)*vk::topo(y_)").unwrap())
                .with(parse!("vk::topo(x_ * y_)").unwrap());

            // map loop momenta and external momenta
            for (ei, e) in external_edges.iter().enumerate() {
                integrand_vakint = integrand_vakint
                    .replace(parse!(format!("vk::k({},x__)", usize::from(*e))).unwrap())
                    .with(parse!(format!("vk::p({},x__)", ei)).unwrap());
            }
            for (ei, e) in lmb.iter().enumerate() {
                // TODO: check if e is in the reduced graph
                // we need to rewrite loop momentum combinations to a new loop momentum
                integrand_vakint = integrand_vakint
                    .replace(parse!(format!("vk::k({},x__)", usize::from(*e))).unwrap())
                    .with(parse!(format!("vk::k({},x__)", ei)).unwrap());
            }

            println!("Integrand vakint: {}", integrand_vakint);

            let mut vakint_expr = VakintExpression::try_from(integrand_vakint.clone()).unwrap();
            println!("\nVakint expression:\n{}", vakint_expr);*/

            //vakint_expr.evaluate_integral(&VAKINT).unwrap();

            // Convert the numerator of the first integral to a dot notation
            /*vakint_expr.0[0].numerator =
                Vakint::convert_to_dot_notation(vakint_expr.0[0].numerator.as_view());
            println!("\nInput integral in dot notation:\n{}\n", vakint_expr);

            //let integral = VAKINT.evaluate(integral.as_view()).unwrap();
            //println!("Evaluated integral:\n{}\n", integral.clone());

            // Set some value for the mass parameters
            let params = VAKINT.params_from_f64(
                &[("muvsq".into(), 3.0), ("mursq".into(), 5.0)]
                    .iter()
                    .cloned()
                    .collect(),
            );

            // And for the external momenta part of the numerator
            let externals = VAKINT.externals_from_f64(
                &(1..=2)
                    .map(|i| {
                        (
                            i,
                            (
                                0.17 * ((i + 1) as f64),
                                0.4 * ((i + 2) as f64),
                                0.3 * ((i + 3) as f64),
                                0.12 * ((i + 4) as f64),
                            ),
                        )
                    })
                    .collect(),
            );

            let (eval, error) = VAKINT
                .numerical_evaluation(integral.as_view(), &params, Some(&externals))
                .unwrap();
            println!("Numerical evaluation:\n{}\n", eval);
            let eval_atom = eval.to_atom(vakint_symbol!(VAKINT.settings.epsilon_symbol.clone()));
            println!("Numerical evaluation as atom:\n{}\n", eval_atom);
            #[rustfmt::skip]
            let target_eval =  NumericalEvaluationResult::from_vec(
            vec![
                    (-3, ( "0.0".into(),  "-10202.59860843888064555902993586".into()),),
                    (-2, ( "0.0".into(),  "62122.38565651740465978420334366".into()),),
                    (-1, ( "0.0".into(),  "-188670.2193437045050954664088623".into()),),
                    ( 0, ( "0.0".into(),  "148095.4883501202267659938351786".into()),),
                ],
                &VAKINT.settings);
            let (matches, match_msg) = target_eval.does_approx_match(
                &eval,
                error.as_ref(),
                10.0_f64.powi(-((VAKINT.settings.run_time_decimal_precision - 4) as i32)),
                1.0,
            );
            if matches {
                println!("Numerical evaluation matches target result.");
            } else {
                println!(
                    "Numerical evaluation does not match target result:\n{}",
                    match_msg
                );
            }*/

            ApproxOp::Dependent {
                t_arg: IntegrandExpr { integrand: a },
                sign: -sign,
                subgraph: reduced,
            }
        } else {
            ApproxOp::NotComputed
        }
    }

    pub fn compute_cff(
        &self,
        dependent: &Self,
        orientation: &OrientationData,
        canonize_esurface: &Option<ShiftRewrite>,
        graph: &UVGraph,
    ) -> ApproxOp {
        ApproxOp::cff_expr(
            dependent,
            &self.subgraph,
            orientation,
            canonize_esurface,
            graph,
        )
    }

    pub fn union(&self, dependents: &[&Self], _graph: &UVGraph) -> (ApproxOp, SimpleApprox) {
        let simple_approx = SimpleApprox::union(
            self.subgraph.clone(),
            dependents.iter().map(|s| s.simple_approx.as_ref().unwrap()),
        );
        (ApproxOp::union(dependents).unwrap(), simple_approx)
    }

    /// Get final expression in the forest sum
    pub fn unwrapped_expr(
        &self,
        graph: &UVGraph,
        amplitude: &InternalSubGraph,
    ) -> Option<SerializableAtom> {
        let (t, s) = self.t_op.expr()?;

        let contracted =
            s * IntegrandExpr::from_subgraph(&amplitude.subtract(&self.subgraph), graph).integrand;

        // FIXME: we are likely mapping in too many momenta, only replace the contracted graph momenta
        Some(
            Self::simplify_notation(&(t * contracted))
                .replace_multiple(&graph.lmb_replacement)
                .into(),
        )
    }

    pub fn final_expr<S: SubGraph<Base = BitVec>>(
        &self,
        graph: &UVGraph,
        amplitude: &S,
    ) -> Option<Atom> {
        let (t, s) = self.t_op.expr()?;

        let reduced = amplitude.included().subtract(&self.subgraph.included());

        let contracted = s * IntegrandExpr::from_subgraph(&reduced, graph).integrand;

        let mut res = t * contracted;

        // set the momenta flowing through the reduced graph edges to the identity wrt the supergraph
        for (e, ei, _) in graph.iter_edges_of(&reduced) {
            let edge_id = usize::from(ei) as i64;
            if e.is_paired() {
                res = res
                    .replace(function!(GS.emr_mom, edge_id, W_.y_))
                    .with(function!(
                        GS.emr_mom,
                        edge_id,
                        function!(GS.emr_mom, edge_id),
                        W_.y_
                    ));
            }
        }

        Some(res)
    }

    pub fn final_cff<S: SubGraph>(
        &self,
        graph: &UVGraph,
        canonize_esurface: &Option<ShiftRewrite>,
        amplitude: &S,
        orientation: &OrientationData,
    ) -> Option<Atom> {
        let (t, s) = self.cff_expr.expr()?;

        let contracted_edges = graph
            .iter_edges_of(&self.subgraph)
            .filter_map(|(p, e, _)| {
                if matches!(p, HedgePair::Paired { .. }) {
                    Some(e)
                } else {
                    None
                }
            })
            .collect_vec();

        let contracted = s * generate_uv_cff(
            graph,
            amplitude,
            canonize_esurface,
            &contracted_edges,
            &orientation.orientation,
        )
        .unwrap();

        Some(t * contracted)
    }

    pub fn simple_expr(
        &self,
        graph: &UVGraph,
        amplitude: &InternalSubGraph,
    ) -> Option<SerializableAtom> {
        let simple_approx = self.simple_approx.as_ref()?;

        Some((simple_approx.sign * simple_approx.expr(&amplitude.filter)).into())
    }
}

impl ApproxOp {
    pub fn sign(&self) -> Option<Sign> {
        match self {
            ApproxOp::NotComputed => None,
            ApproxOp::Union { sign, .. } => Some(*sign),
            ApproxOp::Dependent { sign, .. } => Some(*sign),
            ApproxOp::Root => Some(Sign::Positive),
        }
    }

    pub fn expr(&self) -> Option<(Atom, Sign)> {
        match self {
            ApproxOp::NotComputed => None,
            ApproxOp::Union { t_args, sign, .. } => {
                let mut mul = Atom::num(1);
                for t in t_args {
                    mul = mul * &t.integrand;
                }
                Some((mul, *sign))
            }
            ApproxOp::Dependent { t_arg, sign, .. } => Some((t_arg.integrand.clone(), *sign)),
            ApproxOp::Root => Some((Atom::num(1).into(), Sign::Positive)),
        }
    }

    pub fn cff_expr(
        dependent: &Approximation,
        subgraph: &InternalSubGraph,
        orientation: &OrientationData,
        canonize_esurface: &Option<ShiftRewrite>,
        graph: &UVGraph,
    ) -> Self {
        let reduced = subgraph.subtract(&dependent.subgraph);

        if let Some((inner_t, sign)) = dependent.cff_expr.expr() {
            let contracted_edges = graph
                .iter_edges_of(&dependent.subgraph)
                .filter_map(|(p, e, _)| {
                    if matches!(p, HedgePair::Paired { .. }) {
                        Some(e)
                    } else {
                        None
                    }
                })
                .collect_vec();

            let cff = generate_uv_cff(
                graph,
                subgraph,
                canonize_esurface,
                &contracted_edges,
                &orientation.orientation,
            )
            .unwrap()
                * inner_t;

            Self::Dependent {
                t_arg: IntegrandExpr { integrand: cff },
                sign: -sign,
                subgraph: reduced,
            }
        } else {
            Self::NotComputed
        }
    }

    pub fn union(dependent: &[&Approximation]) -> Option<Self> {
        let mut t_args = vec![];
        let mut subgraphs = vec![];

        let mut final_sign = Sign::Positive;
        for d in dependent {
            match &d.t_op {
                ApproxOp::Dependent {
                    t_arg,
                    sign,
                    subgraph,
                } => {
                    t_args.push(t_arg.clone());
                    final_sign = final_sign * *sign;
                    subgraphs.push(subgraph.clone())
                }
                _ => return None,
            }
        }

        Some(Self::Union {
            t_args,
            sign: final_sign,
            subgraphs,
        })
    }

    pub fn cff_union(dependent: &[&Approximation]) -> Option<Self> {
        let mut t_args = vec![];
        let mut subgraphs = vec![];

        let mut final_sign = Sign::Positive;
        for d in dependent {
            match &d.cff_expr {
                ApproxOp::Dependent {
                    t_arg,
                    sign,
                    subgraph,
                } => {
                    t_args.push(t_arg.clone());
                    final_sign = final_sign * *sign;
                    subgraphs.push(subgraph.clone())
                }
                _ => return None,
            }
        }

        Some(Self::Union {
            t_args,
            sign: final_sign,
            subgraphs,
        })
    }

    pub fn is_computed(&self) -> bool {
        match self {
            ApproxOp::NotComputed => false,
            _ => true,
        }
    }
}

pub struct Forest {
    dag: DAG<Approximation, DagNode, ()>,
}

impl Forest {
    pub fn n_terms(&self) -> usize {
        self.dag.nodes.len()
    }

    pub fn compute(&mut self, graph: &UVGraph) {
        let order = self.dag.compute_topological_order();

        for n in order {
            let current = &self.dag.nodes[n];

            let parents = current.parents.clone();

            let (approx, simple_approx) = if parents.is_empty() {
                let subgraph: InternalSubGraph = graph.empty_subgraph();
                (ApproxOp::Root, SimpleApprox::root(subgraph))
            } else if parents.len() == 1 {
                let approx = current
                    .data
                    .approximate(&self.dag.nodes[parents[0]].data, graph);
                (
                    approx,
                    current
                        .data
                        .structure_approximate(&self.dag.nodes[parents[0]].data),
                )
            } else {
                let mut dependents = vec![];
                for p in parents {
                    dependents.push(&self.dag.nodes[p].data);
                }
                current.data.union(&dependents, graph)
            };
            // if let Some(cff) = cff {
            //     self.dag.nodes[n].data.cff = cff;
            // }
            self.dag.nodes[n].data.t_op = approx;
            self.dag.nodes[n].data.simple_approx = Some(simple_approx);
        }
    }

    pub fn compute_cff(
        &mut self,
        graph: &UVGraph,
        orientation: &OrientationData,
        canonize_esurface: &Option<ShiftRewrite>,
    ) {
        let order = self.dag.compute_topological_order();

        for n in order {
            let current = &self.dag.nodes[n];

            let parents = current.parents.clone();

            let approx = if parents.is_empty() {
                ApproxOp::Root
            } else if parents.len() == 1 {
                current.data.compute_cff(
                    &self.dag.nodes[parents[0]].data,
                    orientation,
                    canonize_esurface,
                    graph,
                )
            } else {
                let mut dependents = vec![];
                for p in parents {
                    dependents.push(&self.dag.nodes[p].data);
                }
                ApproxOp::cff_union(&dependents).unwrap()
            };
            self.dag.nodes[n].data.cff_expr = approx;
        }
    }

    ///4d expr with unwrapped dens
    pub fn expr(&self, graph: &UVGraph, amplitude: &InternalSubGraph) -> Option<SerializableAtom> {
        let mut sum = Atom::num(0).into();
        for (_, n) in &self.dag.nodes {
            sum = sum + n.data.unwrapped_expr(graph, amplitude)?;
        }

        Some(sum)
    }

    pub fn local_expr<S: SubGraph<Base = BitVec>>(
        &self,
        graph: &UVGraph,
        amplitude: &S,
        canonize_esurface: &Option<ShiftRewrite>,
        orientation: &OrientationData,
    ) -> Atom {
        let mut sum = Atom::new();

        for (_, n) in &self.dag.nodes {
            sum += n.data.final_expr(graph, amplitude).unwrap()
                * n.data
                    .final_cff(graph, canonize_esurface, amplitude, orientation)
                    .unwrap();
        }

        println!("SUM {:>}", sum.expand());

        sum.expand().map_terms_single_core(|t| {
            let mut expr = t.to_owned();
            let mut data = Vec::new();
            for m in t.pattern_match(
                &function!(GS.den, W_.edgeid_, W_.mom_, W_.mass_, W_.x_)
                    .pow(Atom::var(W_.c_))
                    .to_pattern(),
                None,
                None,
            ) {
                let eid: i64 = m.get(&W_.edgeid_).unwrap().try_into().unwrap();
                let pow = -i64::try_from(m.get(&W_.c_).unwrap()).unwrap();
                let mom = m.get(&W_.mom_).unwrap().clone();
                let mass = m.get(&W_.mass_).unwrap().clone();

                expr = expr
                    .replace(function!(GS.den, eid, W_.y__))
                    .with(Atom::num(1));

                data.push((eid, pow, mom, mass));
            }

            let orientation = orientation.clone();

            // split momenta into energies and spatial part
            // modified OSE represented as OSE(edge_id, momentum, mass_sq, index)
            // in the next step, a mass will be added
            // take derivative of raised propagators
            for (edge_id, pow, momentum, new_mass_sq) in data {
                let orientation = orientation.orientation.clone();
                let new_mass_sqm = new_mass_sq.clone();
                let momentumm = momentum.clone();
                expr = expr
                    .replace(function!(GS.emr_mom, edge_id, W_.mom_, W_.x___))
                    .with_map(move |m| {
                        let momentum2 = m.get(W_.mom_).unwrap().to_atom();
                        assert_eq!(
                            momentumm, momentum2,
                            "{} vs {} for edge {}",
                            momentumm, momentum2, edge_id
                        );
                        let index = m.get(W_.x___).unwrap().to_atom();

                        let sign = SignOrZero::from(
                            orientation[EdgeIndex::from(edge_id as usize)].clone(),
                        ) * 1;

                        function!(GS.ose, edge_id, momentumm, new_mass_sqm, index) * sign
                            + function!(GS.emr_vec, momentumm, index)
                    });

                for _ in 1..pow {
                    expr = expr
                        .replace(function!(GS.ose, edge_id, W_.x___))
                        .with(function!(GS.ose, edge_id, W_.x___) * GS.rescale) // TODO: check if derivative is in the correct quantity
                        .derivative(GS.rescale)
                        .replace(GS.rescale)
                        .with(Atom::num(1));
                }

                expr = expr / Integer::factorial(pow as u32 - 1);

                // set OSE from CFF with the proper momentum and mass
                expr = expr.replace(function!(GS.ose, edge_id)).with(function!(
                    GS.ose,
                    edge_id,
                    momentum,
                    new_mass_sq
                ));
            }

            expr
        })
    }

    pub fn simple_expr(
        &self,
        graph: &UVGraph,
        amplitude: &InternalSubGraph,
    ) -> Option<SerializableAtom> {
        let mut sum = Atom::num(0).into();
        for (_, n) in &self.dag.nodes {
            sum = sum + n.data.simple_expr(graph, amplitude)?;
        }

        Some(sum)
    }

    pub fn structure_and_res(&self, graph: &UVGraph, amplitude: &InternalSubGraph) -> String {
        let mut out = String::new();

        for (_, n) in &self.dag.nodes {
            out.push_str(&format!(
                "{}:{}\n",
                n.data.simple_expr(graph, amplitude).unwrap(),
                n.data.final_expr(graph, amplitude).unwrap()
            ));
        }
        out
    }

    pub fn graphs(&self) -> String {
        self.dag
            .to_dot_impl(&|n| format!("label=S_{}", n.data.subgraph.string_label()))
    }

    pub fn show_structure(&self, graph: &UVGraph, amplitude: &InternalSubGraph) -> Option<String> {
        let mut out = String::new();

        out.push_str(
            &self
                .simple_expr(graph, amplitude)?
                .0
                .printer(PrintOptions {
                    terms_on_new_line: true,
                    ..Default::default()
                })
                .to_string(),
        );
        Some(out)
    }
}

impl Wood {
    pub fn show_graphs(&self, graph: &UVGraph) -> String {
        let mut out = String::new();
        out.push_str("Poset structure:\n");
        out.push_str(&self.poset.dot_structure());

        out.push_str("Graphs:\n");
        for (k, n) in self.poset.nodes.iter() {
            out.push_str(&graph.dot_impl(
                &n.data.to_hairy_subgraph(&graph),
                format!(
                    "dod={};nodeid ={};\n",
                    graph.dod(&n.data),
                    self.poset.dot_id(k)
                ),
                &|e| Some(format!("dod={}", e.dod)),
                &|n| Some(format!("dod={}", n.dod)),
            ));
            out.push('\n');
        }
        out
    }
}

use slotmap::{new_key_type, Key, SecondaryMap, SlotMap};
use std::collections::HashSet;

// Define a new key type for the Poset
new_key_type! {
    pub struct PosetNode;
    pub struct DagNode;
    pub struct UnfoldedWoodNode;
    pub struct CoverSetNode;
}

/// Trait to define DOT attributes for node data.
pub trait DotAttrs {
    fn dot_attrs(&self) -> String;
}

/// A node in the poset, storing generic data, edges to child nodes, and references to parent odes.
#[derive(Debug, Eq)]
pub struct SlotNode<T, R: Key> {
    data: T,
    order: Option<u64>,
    id: R,
    parents: Vec<R>,  // References to parent nodes by key
    children: Vec<R>, // Edges to child nodes by key
}

impl<T: Hash, R: Key> Hash for SlotNode<T, R> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.data.hash(state);
        self.id.hash(state);
    }
}

impl<T: PartialEq, R: Key> PartialEq for SlotNode<T, R> {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data && self.id == other.id
    }
}

impl<T: PartialEq, R: Key> PartialOrd for SlotNode<T, R> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.order.unwrap().partial_cmp(&other.order.unwrap())
    }
}

impl<T: PartialEq + Eq, R: Key> Ord for SlotNode<T, R> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.order.unwrap().cmp(&other.order.unwrap())
    }
}

impl<T, R: Key> SlotNode<T, R> {
    pub fn dot_id(&self, shift: u64) -> String {
        base62::encode(self.id() - shift)
    }

    pub fn to_topo_ordered(&self) -> Result<TopoOrdered<T>>
    where
        T: Clone,
    {
        Ok(TopoOrdered::new(
            self.data.clone(),
            self.order
                .ok_or_else(|| eyre!("Node has no topological order"))?,
        ))
    }

    pub fn in_degree(&self) -> usize {
        self.parents.len()
    }

    pub fn out_degree(&self) -> usize {
        self.children.len()
    }

    pub fn id(&self) -> u64 {
        self.id.data().as_ffi()
    }

    pub fn add_parent(&mut self, parent: R) {
        self.parents.push(parent);
    }

    pub fn add_child(&mut self, child: R) {
        self.children.push(child);
    }

    pub fn remove_child(&mut self, child: R) {
        self.children.retain(|&c| c != child);
    }

    pub fn remove_parent(&mut self, parent: R) {
        self.parents.retain(|&c| c != parent);
    }

    pub fn is_parent_of(&self, child: R) -> bool {
        self.children.contains(&child)
    }

    pub fn new(data: T, id: R) -> Self {
        SlotNode {
            data,
            id,
            order: None,
            parents: Vec::new(),
            children: Vec::new(),
        }
    }

    pub fn compare_data(&self, other: &Self) -> Option<std::cmp::Ordering>
    where
        T: PartialOrd,
    {
        self.data.partial_cmp(&other.data)
    }
}

/// A partially ordered set (poset) that can be built from an iterator and a slotmap.
pub struct DAG<T, R: Key, D = ()> {
    pub nodes: SlotMap<R, SlotNode<T, R>>,
    associated_data: SecondaryMap<R, D>,
}

impl<T, R: Key, D> DAG<T, R, D> {
    pub fn n_nodes(&self) -> usize {
        self.nodes.len()
    }
}

pub type Poset<T, D> = DAG<T, PosetNode, D>;
pub type HasseDiagram<T, D> = DAG<T, CoverSetNode, D>;

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct TopoOrdered<T> {
    pub data: T,
    order: u64,
}

impl<T> TopoOrdered<T> {
    pub fn new(data: T, order: u64) -> Self {
        TopoOrdered { data, order }
    }
}

#[allow(clippy::non_canonical_partial_ord_impl)]
impl<T> PartialOrd for TopoOrdered<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.order.partial_cmp(&other.order)
    }
}

impl<T> PartialEq for TopoOrdered<T> {
    fn eq(&self, other: &Self) -> bool {
        self.order == other.order
    }
}

impl<T> Eq for TopoOrdered<T> {}

impl<T> Ord for TopoOrdered<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.order.cmp(&other.order)
    }
}

use color_eyre::Result;

impl<T, R: Key, D> DAG<T, R, D> {
    pub fn new() -> Self {
        DAG {
            nodes: SlotMap::with_key(),
            associated_data: SecondaryMap::new(),
        }
    }

    pub fn children(&self, key: R) -> impl Iterator<Item = R> + '_ {
        self.nodes.get(key).unwrap().children.iter().copied()
    }

    pub fn dot_id(&self, key: R) -> String {
        self.nodes.get(key).unwrap().dot_id(self.shift())
    }

    pub fn node_values(&self) -> impl Iterator<Item = &T> {
        self.nodes.values().map(|node| &node.data)
    }

    pub fn add_edge(&mut self, from: R, to: R) {
        let from_node = self.nodes.get_mut(from).unwrap();
        from_node.add_child(to);
        let to_node = self.nodes.get_mut(to).unwrap();
        to_node.add_parent(from);
    }

    pub fn add_edge_if_new(&mut self, from: R, to: R) {
        if !self.nodes.get(from).unwrap().is_parent_of(to) {
            self.add_edge(from, to);
        }
    }

    pub fn remove_edge(&mut self, from: R, to: R) {
        let from_node = self.nodes.get_mut(from).unwrap();
        from_node.remove_child(to);
        let to_node = self.nodes.get_mut(to).unwrap();
        to_node.remove_parent(from);
    }

    pub fn bfs_reach<'a>(
        &'a self,
        start: &'a R,
    ) -> BfsReachable<&'a R, impl FnMut(&'a &'a R) -> &'a [R]> {
        pathfinding::directed::bfs::bfs_reach(start, |&s| self.succesors(*s))
    }

    pub fn maximum(&self) -> Option<R>
    where
        T: Eq,
    {
        Some(self.nodes.values().max()?.id)
    }

    pub fn minimum(&self) -> Option<R>
    where
        T: Eq,
    {
        Some(self.nodes.values().min()?.id)
    }

    pub fn succesors(&self, node_key: R) -> &[R] {
        &self.nodes.get(node_key).unwrap().children
    }

    /// Returns an iterator over all paths starting from the root node, traversed in BFS order.
    pub fn bfs_paths(&self) -> BfsPaths<T, R>
    where
        T: Eq,
    {
        let mut queue = VecDeque::new();
        queue.push_back(vec![self.minimum().unwrap()]);
        BfsPaths {
            queue,
            visited: HashSet::new(),
            nodes: &self.nodes,
        }
    }

    pub fn invert(&mut self) {
        self.nodes.iter_mut().for_each(|(_, node)| {
            std::mem::swap(&mut node.children, &mut node.parents);
        });
    }

    pub fn bfs_paths_inv(&self) -> BfsPaths<T, R> {
        let mut queue = VecDeque::new();
        let mut maximal_elements = Vec::new();
        let mut has_incoming = HashSet::new();

        for node in self.nodes.values() {
            for &parent in node.parents.iter() {
                has_incoming.insert(parent);
            }
        }

        for (key, _node) in self.nodes.iter() {
            if !has_incoming.contains(&key) {
                maximal_elements.push(key);
            }
        }

        for &max in &maximal_elements {
            queue.push_back(vec![max]);
        }

        BfsPaths {
            queue,
            visited: HashSet::new(),
            nodes: &self.nodes,
        }
    }

    pub fn data(&self, key: R) -> &T {
        &self.nodes.get(key).unwrap().data
    }

    fn shift(&self) -> u64 {
        self.nodes.iter().next().unwrap().1.id()
    }

    pub fn to_dot(&self, label: &impl Fn(&T) -> String) -> String {
        self.to_dot_impl(&|node| label(&node.data))
    }

    fn to_dot_impl(&self, label: &impl Fn(&SlotNode<T, R>) -> String) -> String {
        let mut dot = String::new();
        dot.push_str("digraph Poset {\n");
        dot.push_str("    node [shape=circle];\n");

        let shift = self.shift();

        for node in self.nodes.values() {
            let node_id = node.dot_id(shift);
            dot.push_str(&format!("n{} [{}];\n", node_id, label(node)));
            for &child in node.children.iter() {
                dot.push_str(&format!(
                    "n{} -> n{};\n",
                    node_id,
                    self.nodes.get(child).unwrap().dot_id(shift)
                ));
            }
        }

        dot.push_str("}\n");
        dot
    }

    pub fn dot_structure(&self) -> String {
        let shift = self.shift();
        self.to_dot_impl(&|n| format!("label={}", n.dot_id(shift)))
    }

    pub fn add_node(&mut self, data: T) -> R {
        self.nodes.insert_with_key(|key| SlotNode::new(data, key))
    }

    /// Returns all descendants of the node, used for propagating transitive relations.
    fn get_all_descendants(&self, node_key: R) -> HashSet<R> {
        let mut descendants = HashSet::new();
        let mut stack = vec![node_key];
        while let Some(current_key) = stack.pop() {
            let current_node = self.nodes.get(current_key).unwrap();
            for &child_key in current_node.children.iter() {
                if descendants.insert(child_key) {
                    stack.push(child_key);
                }
            }
        }
        descendants
    }

    /// Returns all ancestors of the node, used for propagating transitive relations.
    fn get_all_ancestors(&self, node_key: R) -> HashSet<R> {
        let mut ancestors = HashSet::new();
        let mut stack = vec![node_key];
        while let Some(current_key) = stack.pop() {
            let current_node = self.nodes.get(current_key).unwrap();
            for &parent_key in current_node.parents.iter() {
                if ancestors.insert(parent_key) {
                    stack.push(parent_key);
                }
            }
        }
        ancestors
    }

    pub fn transitive_edges(&self, a: R) -> Vec<(R, R)> {
        let mut edges = Vec::new();

        let children: AHashSet<_> = self
            .nodes
            .get(a)
            .unwrap()
            .children
            .iter()
            .cloned()
            .collect();

        for &child in &children {
            for descendant in self.get_all_descendants(child) {
                if children.contains(&descendant) {
                    edges.push((a, descendant));
                }
            }
        }
        edges
    }

    pub fn in_degree(&self, key: R) -> usize {
        self.nodes.get(key).unwrap().in_degree()
    }

    pub fn compute_topological_order(&mut self) -> Vec<R> {
        // Initialize queue with nodes having in-degree zero
        let mut queue = VecDeque::new(); //S in the wikipedia article

        let mut indegrees = SecondaryMap::new();

        for (key, node) in self.nodes.iter() {
            indegrees.insert(key, node.in_degree());
            if node.in_degree() == 0 {
                queue.push_back(key);
            }
        }

        let mut order = vec![];
        while let Some(node_key) = queue.pop_front() {
            // Assign the order number to the node
            if let Some(node) = self.nodes.get_mut(node_key) {
                node.order = Some(order.len() as u64);
                order.push(node_key);
            }

            // For each child, decrement its in-degree
            for &child_key in &self.nodes.get(node_key).unwrap().children {
                indegrees[child_key] -= 1;
                if indegrees[child_key] == 0 {
                    queue.push_back(child_key);
                }
            }
        }

        // Optional: Check if graph has cycles
        if order.len() != self.nodes.len() {
            panic!("The graph contains a cycle!");
        }
        order
    }
}

impl<T, R: Key, D> Default for DAG<T, R, D> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T, D> Poset<T, D> {
    pub fn poset_family(&self, data: &T) -> [Vec<PosetNode>; 2]
    where
        T: PartialOrd,
    {
        let mut parents = Vec::new();
        let mut children = Vec::new();
        for (key, node) in self.nodes.iter() {
            match data.partial_cmp(&node.data) {
                Some(std::cmp::Ordering::Greater) => {
                    children.push(key);
                }
                Some(std::cmp::Ordering::Less) => {
                    parents.push(key);
                }
                _ => {}
            }
        }
        [parents, children]
    }
    pub fn poset_push(&mut self, data: T, associated_data: D, flip: bool) -> PosetNode
    where
        T: PartialOrd,
    {
        let id = self.nodes.insert_with_key(|key| SlotNode::new(data, key));

        self.associated_data.insert(id, associated_data);
        let new_node = self.nodes.get(id).unwrap();

        let [mut parents, mut children] = self.poset_family(&new_node.data);

        if flip {
            std::mem::swap(&mut parents, &mut children);
        }

        for &parent_key in &parents {
            self.add_edge(parent_key, id);
        }

        for &child_key in &children {
            self.add_edge(id, child_key);
        }

        self.update_transitive_closure(id);

        id
    }

    /// Updates the transitive closure of the poset by propagating the relationships.
    fn update_transitive_closure(&mut self, new_node_key: PosetNode) {
        // Propagate relationships for all descendants of new_node
        let descendants = self.get_all_descendants(new_node_key);
        for &descendant_key in &descendants {
            self.add_edge_if_new(new_node_key, descendant_key);
        }

        // Propagate relationships for all ancestors of new_node
        let ancestors = self.get_all_ancestors(new_node_key);
        for &ancestor_key in &ancestors {
            self.add_edge_if_new(ancestor_key, new_node_key);
        }
    }

    pub fn remove_transitive_edges(mut self) -> HasseDiagram<T, D> {
        let edges_to_remove: Vec<_> = self
            .nodes
            .keys()
            .flat_map(|node_key| self.transitive_edges(node_key))
            .collect();

        // let shift = self.shift();
        for (a, b) in edges_to_remove {
            // println!(
            //     "removing edge from {} to {}",
            //     base62::encode(a.0.as_ffi() - shift),
            //     base62::encode(b.0.as_ffi() - shift)
            // );
            self.remove_edge(a, b);
        }

        let mut hasse = HasseDiagram::new();

        let mut new_map = SecondaryMap::new();

        for (key, node) in self.nodes.into_iter() {
            new_map.insert(key, hasse.add_node(node.data));
            if let Some(d) = self.associated_data.remove(key) {
                hasse.associated_data.insert(new_map[key], d);
            }
            for child in node.children {
                hasse.add_edge(new_map[key], new_map[child]);
            }
        }

        hasse
    }
}

impl<T: PartialOrd, D> FromIterator<(T, D)> for Poset<T, D> {
    fn from_iter<I: IntoIterator<Item = (T, D)>>(iter: I) -> Self {
        let mut poset = Poset::new();
        for (data, assoc) in iter {
            poset.poset_push(data, assoc, false);
        }
        poset
    }
}

/// An iterator over paths in the poset, traversed in BFS order.
pub struct BfsPaths<'a, T, R: Key> {
    queue: VecDeque<Vec<R>>,
    visited: HashSet<Vec<R>>,
    nodes: &'a SlotMap<R, SlotNode<T, R>>,
}

impl<'a, T, R: Key> Iterator for BfsPaths<'a, T, R> {
    type Item = Vec<R>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(path) = self.queue.pop_front() {
            if !self.visited.insert(path.clone()) {
                continue;
            }

            let last_node_key = path.last().unwrap();
            let last_node = self.nodes.get(*last_node_key).unwrap();

            for &child_key in last_node.children.iter() {
                if !path.contains(&child_key) {
                    let mut new_path = path.clone();
                    new_path.push(child_key);
                    self.queue.push_back(new_path);
                }
            }

            return Some(path);
        }
        None
    }
}
#[cfg(test)]
mod tests;
