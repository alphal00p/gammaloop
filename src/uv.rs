use std::{
    cmp::Ordering,
    collections::{HashMap, VecDeque},
    fmt::{self, Display, Formatter},
    hash::Hash,
    ops::{Deref, Index},
};

use crate::{
    cff::generation::generate_cff_expression_from_subgraph_to_ose_atom, graph::VertexInfo,
    model::ArcParticle, momentum::Sign, new_graph::LoopMomentumBasis, utils::GS,
};
use ahash::{AHashMap, AHashSet};
use bitvec::vec::BitVec;
use color_eyre::Report;
use eyre::eyre;
use indexmap::IndexMap;
use itertools::Itertools;
use pathfinding::prelude::BfsReachable;
use rand::{rngs::SmallRng, Rng, SeedableRng};
use ref_ops::RefNeg;
use serde::{Deserialize, Serialize};
use spenso::{
    data::StorageTensor,
    parametric::{atomcore::PatternReplacement, ParamTensor},
    structure::{HasStructure, VecStructure},
    symbolica_utils::SerializableAtom,
};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, Symbol},
    function,
    id::{AtomMatchIterator, Pattern, Replacement},
    parse,
    printer::PrintOptions,
    state::State,
    symbol,
};
use symbolica_community::physics::algebraic_simplification::metric::MetricSimplifier;
use trie_rs::{try_collect::TryFromIterator, Trie, TrieBuilder};

use linnet::half_edge::{
    builder::HedgeGraphBuilder,
    involution::{EdgeIndex, Hedge, HedgePair, SignOrZero},
    subgraph::{cycle::SignedCycle, ModifySubgraph},
    tree::SimpleTraversalTree,
    EdgeAccessors, HedgeGraph, NodeIndex, PowersetIterator,
};
use linnet::half_edge::{
    involution::Flow,
    subgraph::{HedgeNode, InternalSubGraph, SubGraph, SubGraphOps},
};

use crate::{
    graph::{BareEdge, BareGraph, BareVertex, EdgeType},
    model::normalise_complex,
};

use bitvec::prelude::*;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct UVEdge {
    og_edge: usize,
    dod: i32,
    particle: ArcParticle,
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
            den: Atom::new_num(1).into(), // edge.full_den(bare_graph, index).into(),
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
                color: Some(ParamTensor::param(color.map_structure(VecStructure::from))),
            }
        } else {
            UVNode {
                dod: vertex.dod() as i32,
                num: Atom::new_num(1).into(),
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
    cut_edges: BitVec,
    lmb_replacement: Vec<Replacement>,
}

impl Deref for UVGraph {
    type Target = HedgeGraph<UVEdge, UVNode>;
    fn deref(&self) -> &Self::Target {
        &self.hedge_graph
    }
}

#[allow(dead_code)]
impl UVGraph {
    pub fn from_hedge(hedge_graph: HedgeGraph<UVEdge, UVNode>) -> Self {
        let cut_edges = hedge_graph.cycle_basis().1.tree_subgraph;

        let mut uv_graph = UVGraph {
            hedge_graph,
            cut_edges,
            lmb_replacement: vec![],
        };

        let reps = uv_graph
            .lmb_reps_for_subgraph(&uv_graph.full_graph())
            .0
            .iter()
            .map(|(edge, mom)| {
                let r = Replacement::new(
                    function!(GS.emr_mom, usize::from(*edge) as i64).to_pattern(),
                    mom.to_pattern(),
                );
                r
            })
            .collect();

        uv_graph.lmb_replacement = reps;
        uv_graph
    }

    /// Returns momentum assignment, external edges and lmb edges
    pub fn lmb_reps_for_subgraph(
        &self,
        subgraph: &InternalSubGraph,
        // lmb: &[EdgeIndex],
    ) -> (HashMap<EdgeIndex, Atom>, Vec<EdgeIndex>, Vec<EdgeIndex>) {
        let lmb = self.select_lmb(subgraph);
        let mut edge_rep: HashMap<_, Atom> = HashMap::default();

        let Some(i) = subgraph.filter.first_one() else {
            return (edge_rep, vec![], lmb);
        };

        let mut tree: BitVec = self.full_filter();

        let externals = self.nesting_node_from_subgraph(subgraph.clone()).hairs;

        let root_node = self.node_id(Hedge(externals.first_one().unwrap_or(i)));

        for e in &lmb {
            let (_, p) = self[e];
            tree.sub(p);
        }

        let trav_tree =
            SimpleTraversalTree::depth_first_traverse(&self, &tree, &root_node, None).unwrap();

        for loop_edge in &lmb {
            let (_, p) = self[loop_edge];
            let loop_id: usize = (*loop_edge).into();
            let loop_mom = function!(GS.loop_mom, loop_id as i64);

            if let Some(c) = trav_tree.get_cycle(p.any_hedge(), &self.hedge_graph) {
                if let HedgePair::Paired { source, .. } = p {
                    if let Some(signed) = SignedCycle::from_cycle(c, source, self) {
                        for h in signed.filter.included_iter() {
                            let eid = self[&h];
                            let flow = self.flow(h);

                            match flow {
                                Flow::Source => {
                                    *edge_rep.entry(eid).or_insert(Atom::Zero) += &loop_mom;
                                    //Validated the sign on 13.05.2025
                                }
                                Flow::Sink => {
                                    *edge_rep.entry(eid).or_insert(Atom::Zero) -= &loop_mom;
                                }
                            }
                        }
                    }
                }
            }
        }

        let mut externals_ids = vec![];
        for ext in externals.included_iter().skip(1) {
            let ext_sign: SignOrZero = self.flow(ext).into();
            externals_ids.push(self[&ext]);
            let ext_id: usize = self[&ext].into();
            let ext_mom = match ext_sign {
                SignOrZero::Plus => {
                    function!(GS.emr_mom, ext_id as i64)
                }
                SignOrZero::Minus => -function!(GS.emr_mom, ext_id as i64),
                SignOrZero::Zero => {
                    panic!("Missing external momentum sign")
                }
            };

            for h in trav_tree
                .ancestor_iter_hedge(ext, &self.hedge_graph.as_ref())
                .step_by(2)
            {
                let eid = self[&h];
                let (_, p) = &self[&eid];

                if let HedgePair::Paired { source, sink } = p {
                    if h == *source {
                        *edge_rep.entry(eid).or_insert(Atom::Zero) -= &ext_mom;
                    } else if h == *sink {
                        *edge_rep.entry(eid).or_insert(Atom::Zero) += &ext_mom;
                    } else {
                        panic!("Should be in HedgePair");
                    }
                }
            }
        }

        (edge_rep, externals_ids, lmb)
    }

    pub fn select_lmb(&self, subgraph: &InternalSubGraph) -> Vec<EdgeIndex> {
        let n_loops = self.n_loops(subgraph);

        if n_loops == 0 {
            return vec![];
        }

        // the subgraph may have disconnected components in case the of disjoint graphs in a spinney
        let components = self.count_connected_components(subgraph);
        let cut_edges_in_subgraph = subgraph.filter.intersection(&self.cut_edges);

        for v in self
            .iter_edges(&cut_edges_in_subgraph)
            .combinations(n_loops)
        {
            let mut cut_subgraph = subgraph.filter.clone();
            let mut lmb = vec![];

            for (p, e, _) in v {
                cut_subgraph.sub(p);
                lmb.push(e);
            }

            if self.count_connected_components(&cut_subgraph) == components {
                return lmb;
            }
        }

        panic!(
            "No lmb found for {} and cut edges {}",
            self.dot(subgraph),
            self.dot(&self.cut_edges)
        )
    }

    pub fn from_graph(graph: &BareGraph) -> Self {
        let mut excised: BitVec = graph.hedge_representation.empty_subgraph();

        for (n, _, d) in graph.hedge_representation.iter_nodes() {
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

        let cut_edges = SimpleTraversalTree::depth_first_traverse(
            &excised,
            &excised.full_filter(),
            &NodeIndex(0),
            None,
        )
        .unwrap()
        .tree_subgraph;

        let mut uv_graph = UVGraph {
            hedge_graph: excised,
            cut_edges,
            lmb_replacement: vec![],
        };

        let reps = uv_graph
            .lmb_reps_for_subgraph(&uv_graph.full_graph())
            .0
            .iter()
            .map(|(edge, mom)| {
                let r = Replacement::new(
                    function!(GS.emr_mom, usize::from(*edge) as i64).to_pattern(),
                    mom.to_pattern(),
                );
                // println!("Rep:{r}");
                r
            })
            .collect();

        uv_graph.lmb_replacement = reps;

        uv_graph
    }

    fn spinneys(&self) -> AHashSet<InternalSubGraph> {
        // println!("{}", self.base_dot());
        let mut spinneys: AHashSet<_> = InternalSubGraph::all_ops_iterative_filter_map(
            &self.all_cycle_sym_diffs().unwrap(),
            &|a, b| a.union(b),
            &|union| {
                if self.dod(&union) >= 0 {
                    Some(union)
                } else {
                    None
                }
            },
        );

        spinneys.insert(self.empty_subgraph());

        spinneys
    }

    fn wood(&self) -> Wood {
        Wood::from_spinneys(self.spinneys(), self)
    }

    fn n_loops<S: SubGraph>(&self, subgraph: &S) -> usize {
        // if let Some(loop_count) = subgraph.loopcount {
        //     // println!("found loop_count nloops: {}", loop_count);
        //     loop_count
        // } else {
        self.cyclotomatic_number(subgraph)
        // }
    }

    fn dod<S: SubGraph>(&self, subgraph: &S) -> i32 {
        let mut dod: i32 = 4 * self.n_loops(subgraph) as i32;
        // println!("nloops: {}", dod / 4);

        for e in self.iter_internal_edge_data(subgraph) {
            dod += e.data.dod;
        }

        for (_, _, n) in self.iter_node_data(subgraph) {
            dod += n.dod;
        }

        dod
    }

    fn numerator<S: SubGraph>(&self, subgraph: &S) -> Atom {
        let mut num = Atom::new_num(1);

        for (_, _, n) in self.iter_node_data(subgraph) {
            num = num * &n.num;
        }

        for e in self.iter_internal_edge_data(subgraph) {
            num = num * &e.data.num;
        }

        num.into()
    }

    fn denominator(&self, subgraph: &InternalSubGraph) -> Atom {
        let mut den = Atom::new_num(1);

        for (pair, eid, d) in self.iter_edges(subgraph) {
            if matches!(pair, HedgePair::Paired { .. }) {
                den = den * function!(GS.den, usize::from(eid) as i64, &d.data.den);
            }
        }

        den.into()
    }

    fn dot<S: SubGraph>(&self, subgraph: &S) -> String {
        self.dot_impl(
            subgraph,
            format!("{}", self.dod(subgraph)),
            &|e| Some(format!("label=\"{}\"", e.dod)),
            &|n| Some(format!("label=\"{}\"", n.dod)),
        )
    }

    fn t_op<I: Iterator<Item = InternalSubGraph>>(&self, mut subgraph_iter: I) -> Atom {
        if let Some(subgraph) = subgraph_iter.next() {
            let t = self.t_op(subgraph_iter);
            FunctionBuilder::new(symbol!("Top"))
                .add_arg(&Atom::new_num(self.dod(&subgraph)))
                .add_arg(&t)
                .finish()
        } else {
            Atom::new_num(1)
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
    pub fn from_subgraph(subgraph: &InternalSubGraph, graph: &UVGraph) -> Self {
        let num = graph.numerator(subgraph).to_dots();

        let den = graph.denominator(subgraph);

        // println!("{}", graph.denominator(subgraph).0);
        IntegrandExpr {
            integrand: num / den,
        }
    }
}

impl Wood {
    pub fn n_spinneys(&self) -> usize {
        self.poset.n_nodes()
    }
    pub fn from_spinneys<I: IntoIterator<Item = InternalSubGraph>>(s: I, graph: &UVGraph) -> Self {
        let mut poset = Poset::from_iter(s.into_iter().map(|s| (s, ())));

        poset.invert();
        poset.compute_topological_order();

        let mut unions = SecondaryMap::new();

        for (i, sg) in poset.nodes.iter() {
            let cs = graph.connected_components(&sg.data);

            if cs.len() > 1 {
                // sg is a disjoint union of spinneys (at the level of half-edges) (strongly disjoint)
                let mut union = vec![];

                for &c in sg.parents.iter() {
                    let mut is_in = 0;
                    for comp in &cs {
                        let comp = InternalSubGraph::cleaned_filter_optimist(comp.clone(), &graph);
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

    fn unfold_bfs(
        &self,
        graph: &UVGraph,
        dag: &mut DAG<Approximation, DagNode, ()>,
        unions: &mut SecondaryMap<PosetNode, Option<Vec<(PosetNode, Option<DagNode>)>>>,
        root: PosetNode,
    ) -> DagNode {
        let mut search_front = VecDeque::new();

        let tree_root = dag.add_node(Approximation::new(self.poset.data(root).clone(), graph));
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
                            let child = dag
                                .add_node(Approximation::new(self.poset.data(*c).clone(), graph));
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
                        dag.add_node(Approximation::new(self.poset.data(*c).clone(), graph));
                    dag.add_edge(parent, child);
                    search_front.push_front((*c, child));
                }
            }
        }
        tree_root
    }

    pub fn unfold(&self, graph: &UVGraph) -> Forest {
        let mut dag: DAG<Approximation, DagNode, ()> = DAG::new();

        let root = self.poset.minimum().unwrap();

        let mut unions = SecondaryMap::new();

        for (p, u) in self.additional_unions.iter() {
            let union: Vec<(PosetNode, Option<DagNode>)> = u.iter().map(|i| (*i, None)).collect();
            unions.insert(p, Some(union));
        }

        let _ = self.unfold_bfs(graph, &mut dag, &mut unions, root);

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
        let reduced = Atom::new_var(Self::subgraph_shadow(bigger_graph, &self.graph));
        let mut mul = Atom::new_num(1);
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
    lmb: Vec<EdgeIndex>,
    externals: Vec<EdgeIndex>,
    pub approx_op: ApproxOp,
    momentum_assignment: HashMap<EdgeIndex, Atom>,
    pub mom_rep: Vec<Replacement>,
    pub simple_approx: Option<SimpleApprox>,
}

impl Approximation {
    pub fn simplify_notation(expr: &SerializableAtom) -> SerializableAtom {
        let replacements = [(function!(GS.den, GS.a_, GS.x_), Atom::new_var(GS.x_))];

        let reps: Vec<_> = replacements
            .into_iter()
            .map(|(a, b)| Replacement::new(a.to_pattern(), b))
            .collect();

        expr.replace_multiple_repeat(&reps)
    }

    pub fn new(spinney: InternalSubGraph, graph: &UVGraph) -> Approximation {
        // let lmb = graph.select_lmb(&spinney);

        let (momentum_assignment, externals, lmb) = graph.lmb_reps_for_subgraph(&spinney);

        Approximation {
            dod: graph.dod(&spinney),

            subgraph: spinney,
            lmb,
            mom_rep: momentum_assignment
                .iter()
                .map(|(edge, mom)| {
                    let r = Replacement::new(
                        function!(GS.emr_mom, usize::from(*edge) as i64).to_pattern(),
                        mom.to_pattern(),
                    );
                    // println!("Rep:{r}");
                    r
                })
                .collect(),
            momentum_assignment,
            externals,
            simple_approx: None,
            approx_op: ApproxOp::NotComputed,
        }
    }

    pub fn reduced_graph(&self, subgraph: &InternalSubGraph) -> InternalSubGraph {
        self.subgraph.subtract(subgraph)
    }

    pub fn dependent(&self, dependent: &Self, graph: &UVGraph) -> (ApproxOp, SimpleApprox) {
        let simple_approx = dependent
            .simple_approx
            .as_ref()
            .unwrap()
            .dependent(self.subgraph.clone());

        (
            ApproxOp::dependent(
                dependent,
                &self.subgraph,
                &self.externals,
                &self.mom_rep,
                self.dod,
                graph,
            ),
            simple_approx,
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
    pub fn expr(&self, graph: &UVGraph) -> Option<SerializableAtom> {
        let (t, s) = self.approx_op.expr()?;

        let contracted = s * IntegrandExpr::from_subgraph(
            &graph.full_node().internal_graph.subtract(&self.subgraph),
            graph,
        )
        .integrand;

        Some(Self::simplify_notation(
            &(t * contracted)
                .replace_multiple(&graph.lmb_replacement)
                .into(),
        ))
    }

    pub fn simple_expr(&self, graph: &UVGraph) -> Option<SerializableAtom> {
        let simple_approx = self.simple_approx.as_ref()?;

        Some((simple_approx.sign * simple_approx.expr(&graph.full_filter())).into())
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
                let mut mul = Atom::new_num(1);
                for t in t_args {
                    mul = mul * &t.integrand;
                }
                Some((mul, *sign))
            }
            ApproxOp::Dependent { t_arg, sign, .. } => Some((t_arg.integrand.clone(), *sign)),
            ApproxOp::Root => Some((Atom::new_num(1).into(), Sign::Positive)),
        }
    }

    pub fn dependent(
        dependent: &Approximation,   //is smaller than subgraph
        subgraph: &InternalSubGraph, //is not necessarily full_graph
        external_edges: &[EdgeIndex],
        mom_reps: &[Replacement],
        dod: i32,
        graph: &UVGraph,
    ) -> Self {
        let reduced = subgraph.subtract(&dependent.subgraph);

        if let Some((inner_t, sign)) = dependent.approx_op.expr() {
            let t_arg = IntegrandExpr::from_subgraph(&reduced, graph);

            let mut atomarg = t_arg.integrand * inner_t;

            // rewrite the inner_t as well
            atomarg = atomarg.replace_multiple(mom_reps);

            // println!(
            //     "Expand {} with dod={} in {:?}",
            //     atomarg, dod, external_edges
            // );
            for e in external_edges {
                atomarg = atomarg
                    .replace(function!(GS.emr_mom, usize::from(*e) as i64))
                    .with(function!(GS.emr_mom, usize::from(*e) as i64) * GS.rescale);
            }

            let soft_ct = graph
                .nesting_node_from_subgraph(subgraph.clone())
                .hairs
                .count_ones()
                == 2
                && dod > 0;

            let mut masses = AHashSet::new();
            masses.insert(Atom::new_var(GS.m_uv));
            // scale all masses, including UV masses from subgraphs
            for e in graph.iter_internal_edge_data(subgraph) {
                let e_mass = parse!(&e.data.particle.mass.name).unwrap();
                masses.insert(e_mass);
            }

            if !soft_ct {
                for m in &masses {
                    let rescaled = m.clone() * GS.rescale;
                    atomarg = atomarg.replace(m.clone()).with(rescaled);
                }

                // expand the propagator around a propagator with a UV mass
                atomarg = atomarg
                    .replace(parse!("den(n_,x_)").unwrap())
                    .with(parse!("den(n_,x_- mUV^2 + t^2*mUV^2)").unwrap());
            }

            atomarg = atomarg
                .replace(parse!("symbolica_community::dot(t*x__,y_)").unwrap())
                .repeat()
                .with(parse!("t*symbolica_community::dot(x__,y_)").unwrap());

            // println!("atomarg:{}", atomarg);

            // den(..) tags a propagator, its first derivative is 1 and the rest is 0
            let mut a = atomarg
                .series(GS.rescale, Atom::Zero, dod.into(), true)
                .unwrap()
                .to_atom()
                .replace(parse!("der(0,1, den(y__))").unwrap())
                .with(Atom::new_num(1))
                .replace(parse!("der(x__, den(y__))").unwrap())
                .with(Atom::new_num(0));

            if soft_ct {
                let coeffs = a.coefficient_list::<u8>(&[Atom::new_var(GS.rescale)]);
                let mut b = Atom::Zero;
                let dod_pow = Atom::new_var(GS.rescale).npow(dod);
                for (pow, mut i) in coeffs {
                    if pow == dod_pow {
                        // set the masses in the t=dod term to 0
                        // UV rearrange the denominators
                        for m in &masses {
                            i = i.replace(m.clone()).with(Atom::Zero);
                        }

                        i = i
                            .replace(parse!("den(n_,x_)").unwrap())
                            .with(parse!("den(n_,x_ - mUV^2)").unwrap());
                    }

                    b += i;
                }

                a = b;
            } else {
                a = a.replace(GS.rescale).with(Atom::new_num(1));
            }

            // println!("Expanded: {:>}", a.expand());

            Self::Dependent {
                t_arg: IntegrandExpr { integrand: a },
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
            match &d.approx_op {
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
                current
                    .data
                    .dependent(&self.dag.nodes[parents[0]].data, graph)
            } else {
                let mut dependents = vec![];
                for p in parents {
                    dependents.push(&self.dag.nodes[p].data);
                }
                current.data.union(&dependents, graph)
            };

            self.dag.nodes[n].data.approx_op = approx;
            self.dag.nodes[n].data.simple_approx = Some(simple_approx);
        }
    }

    pub fn expr(&self, graph: &UVGraph) -> Option<SerializableAtom> {
        let mut sum = Atom::new_num(0).into();
        for (_, n) in &self.dag.nodes {
            sum = sum + n.data.expr(graph)?;
        }

        Some(sum)
    }

    pub fn simple_expr(&self, graph: &UVGraph) -> Option<SerializableAtom> {
        let mut sum = Atom::new_num(0).into();
        for (_, n) in &self.dag.nodes {
            sum = sum + n.data.simple_expr(graph)?;
        }

        Some(sum)
    }

    pub fn structure_and_res(&self, graph: &UVGraph) -> String {
        let mut out = String::new();

        for (_, n) in &self.dag.nodes {
            out.push_str(&format!(
                "{}:{}\n",
                n.data.simple_expr(graph).unwrap(),
                n.data.expr(graph).unwrap()
            ));
        }
        out
    }

    pub fn graphs(&self) -> String {
        self.dag
            .to_dot_impl(&|n| format!("label=S_{}", n.data.subgraph.string_label()))
    }

    pub fn show_structure(&self, graph: &UVGraph) -> Option<String> {
        let mut out = String::new();

        out.push_str(
            &self
                .simple_expr(graph)?
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
