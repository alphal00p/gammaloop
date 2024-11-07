use std::{
    cmp::Ordering,
    collections::VecDeque,
    fmt::{self, Display, Formatter},
    hash::Hash,
    ops::{Deref, Index},
};

use crate::{momentum::Sign, utils::GS};
use ahash::{AHashMap, AHashSet};
use bitvec::{slice::IterOnes, vec::BitVec};
use color_eyre::Report;
use eyre::eyre;
use indexmap::IndexMap;
use itertools::Itertools;
use pathfinding::prelude::BfsReachable;
use rand::{rngs::SmallRng, Rng, SeedableRng};
use ref_ops::RefNeg;
use serde::{Deserialize, Serialize};
use spenso::{parametric::ParamTensor, structure::VecStructure, symbolica_utils::SerializableAtom};
use symbolica::{
    atom::{Atom, AtomView, FunctionBuilder, Symbol},
    fun,
    id::{Pattern, Replacement},
    printer::PrintOptions,
    state::State,
    symb,
};
use trie_rs::{try_collect::TryFromIterator, Trie, TrieBuilder};

use crate::{
    graph::{
        half_edge::{
            subgraph::{HedgeNode, InternalSubGraph, SubGraph},
            HedgeGraph, HedgeGraphBuilder, NodeIndex, PowersetIterator,
        },
        BareGraph, Edge, EdgeType, LoopMomentumBasis, Vertex,
    },
    model::normalise_complex,
};

use bitvec::prelude::*;

pub struct PoSet<N>
where
    N: PartialOrd,
{
    greater_than: Vec<Vec<usize>>,
    nodes: Vec<N>,
}

pub struct BfsIterator<N: Index<usize>> {
    order: Vec<usize>,
    nodes: N,
}

impl<N: Index<usize, Output: Sized + Clone>> Iterator for BfsIterator<N> {
    type Item = N::Output;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.nodes[self.order.pop()?].clone())
    }
}

pub struct CoverSet<N>
where
    N: PartialOrd,
{
    covers: Vec<Vec<usize>>,
    nodes: Vec<N>,
}

impl<N: PartialOrd> CoverSet<N> {
    pub fn build_predecessors(&self) -> Vec<Vec<usize>> {
        let mut predecessors = vec![Vec::new(); self.nodes.len()];

        for (node_index, successors) in self.covers.iter().enumerate() {
            for &succ_index in successors {
                predecessors[succ_index].push(node_index);
            }
        }

        predecessors
    }

    pub fn covers(&self, a: usize, b: usize) -> bool {
        self.covers[a].contains(&b)
    }

    pub fn dot(&self) -> String {
        let mut out = "digraph {\n ".to_string();
        out.push_str("  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\";\n ");

        for (i, c) in self.covers.iter().enumerate() {
            for j in c {
                out.push_str(&format!("{} -> {};\n", i, j));
            }
        }

        out += "}";
        out
    }

    pub fn bfs(&self, start: usize) -> Vec<usize> {
        let mut visited = bitvec![usize, Lsb0; 0; self.nodes.len()];
        let mut queue = VecDeque::new();
        let mut order = Vec::new();

        queue.push_back(start);

        while let Some(node) = queue.pop_front() {
            if visited[node] {
                continue;
            }

            visited.set(node, true);
            order.push(node);

            for &n in &self.covers[node] {
                queue.push_back(n);
            }
        }

        order
    }

    pub fn bfs_iter(&self, start: usize) -> BfsIterator<Vec<N>>
    where
        N: Clone,
    {
        BfsIterator {
            order: self.bfs(start),
            nodes: self.nodes.clone(),
        }
    }
}

impl<N> PoSet<N>
where
    N: Clone + PartialOrd,
{
    /// Returns true if a covers b
    /// a covers b if a > b and there is no c such that a > c > b
    ///
    pub fn covers(&self, a: usize, b: usize) -> bool {
        let mut covers = false;
        let mut has_c = false;
        for c in &self.greater_than[b] {
            if *c == a {
                covers = true;
            }

            if self.nodes[a] > self.nodes[*c] {
                has_c = true;
            }
        }
        covers && !has_c
    }
    pub fn to_cover_set(&self) -> CoverSet<N> {
        let mut covers: Vec<Vec<usize>> = Vec::new();
        let nodes = self.nodes.clone();

        for (i, _) in nodes.iter().enumerate() {
            let mut coversacents = Vec::new();

            for j in (i + 1)..nodes.len() {
                if self.covers(j, i) {
                    coversacents.push(j);
                }
            }
            covers.push(coversacents);
        }

        CoverSet { covers, nodes }
    }

    pub fn dot(&self) -> String {
        let mut out = "digraph {\n ".to_string();
        out.push_str("  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\";\n ");

        for (i, c) in self.greater_than.iter().enumerate() {
            for j in c {
                out.push_str(&format!("{} -> {};\n", i, j));
            }
        }

        for (i, _n) in self.nodes.iter().enumerate() {
            out.push_str(&format!("{};\n", i));
        }

        out += "}";
        out
    }

    pub fn greater_than(&self, a: usize) -> Vec<usize> {
        self.greater_than[a].clone()
    }

    pub fn build_predecessors(&self) -> Vec<Vec<usize>> {
        let mut predecessors = vec![Vec::new(); self.nodes.len()];

        for (node_index, successors) in self.greater_than.iter().enumerate() {
            for &succ_index in successors {
                predecessors[succ_index].push(node_index);
            }
        }

        predecessors
    }
}

impl<N> FromIterator<N> for PoSet<N>
where
    N: PartialOrd,
{
    fn from_iter<T: IntoIterator<Item = N>>(iter: T) -> Self {
        let mut nodes = Vec::new();
        let mut greater_than: Vec<Vec<usize>> = Vec::new();

        // Insert topologically
        // we want that if n < m then n comes before m in the nodes list
        for n in iter {
            let mut pos = nodes.len();
            for (i, m) in nodes.iter().enumerate().rev() {
                if n < *m {
                    pos = i;
                }
            }
            nodes.insert(pos, n);
        }

        for (i, n) in nodes.iter().enumerate() {
            let mut greater_thanacents = Vec::new();

            for (j, m) in nodes.iter().enumerate().skip(i) {
                if *m > *n {
                    greater_thanacents.push(j);
                }
            }
            greater_than.push(greater_thanacents);
        }

        PoSet {
            greater_than,
            nodes,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct BitFilter {
    data: u128,
}

impl PartialOrd for BitFilter {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        if self.eq(other) {
            Some(Ordering::Equal)
        } else if self.data & other.data == self.data {
            Some(Ordering::Greater)
        } else if self.data & other.data == other.data {
            Some(Ordering::Less)
        } else {
            None
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
struct UVEdge {
    og_edge: usize,
    dod: i32,
    num: SerializableAtom,
    den: SerializableAtom,
}

impl UVEdge {
    pub fn from_edge(edge: &Edge, id: usize, bare_graph: &BareGraph) -> Self {
        let index = (id) as i32;
        let [colorless, _] = edge.color_separated_numerator(bare_graph);
        UVEdge {
            og_edge: id,
            dod: edge.dod() as i32,
            num: normalise_complex(&colorless).into(),
            den: edge.full_den(bare_graph, index).into(),
        }
    }

    pub fn numerator_rank(&self) -> isize {
        0
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct UVNode {
    dod: i32,
    num: SerializableAtom,
    color: Option<ParamTensor>,
}

impl UVNode {
    pub fn from_vertex(vertex: &Vertex, graph: &BareGraph) -> Self {
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

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct UVGraph(HedgeGraph<UVEdge, UVNode>);

#[allow(dead_code)]
impl UVGraph {
    pub fn from_graph(graph: &BareGraph) -> Self {
        let mut uv_graph = HedgeGraphBuilder::new();

        for n in &graph.vertices {
            uv_graph.add_node(UVNode::from_vertex(n, graph));
        }

        for (i, edge) in graph.edges.iter().enumerate() {
            let ves = edge.vertices;

            let sink = NodeIndex(ves[1]);
            let source = NodeIndex(ves[0]);

            match edge.edge_type {
                EdgeType::Virtual => {
                    uv_graph.add_edge(source, sink, UVEdge::from_edge(edge, i, graph));
                }
                EdgeType::Outgoing => {
                    uv_graph.add_external_edge(source, UVEdge::from_edge(edge, i, graph));
                }
                EdgeType::Incoming => {
                    uv_graph.add_external_edge(sink, UVEdge::from_edge(edge, i, graph));
                }
            }
        }

        UVGraph(uv_graph.into())
    }

    pub fn half_edge_id(&self, g: &BareGraph, id: usize) -> usize {
        self.0
            .involution
            .find_from_data(&UVEdge::from_edge(&g.edges[id], id, g))
            .unwrap()
    }

    pub fn node_id(&self, g: &BareGraph, id: usize) -> HedgeNode {
        let e_id = g.vertices[id].edges.first().unwrap();
        self.0
            .involution
            .get_node_id(self.half_edge_id(g, *e_id))
            .clone()
    }

    pub fn cycle_basis_from_lmb(&self, lmb: &LoopMomentumBasis) -> Vec<HedgeNode> {
        let mut cycle_basis = Vec::new();
        let cut_edges = self.0.paired_filter_from_pos(&lmb.basis);

        for v in lmb.basis.iter() {
            let loop_edge = self.0.paired_filter_from_pos(&[*v]);

            let mut hairy_loop =
                self.0
                    .nesting_node_from_subgraph(InternalSubGraph::cleaned_filter_optimist(
                        !cut_edges.clone() | &loop_edge,
                        &self.0,
                    ));

            self.0.cut_branches(&mut hairy_loop);
            cycle_basis.push(hairy_loop);
        }
        cycle_basis
    }

    fn spanning_forest_from_lmb(&self, lmb: LoopMomentumBasis) -> HedgeNode {
        let cutting_edges = self.0.paired_filter_from_pos(&lmb.basis);

        self.0
            .nesting_node_from_subgraph(InternalSubGraph::cleaned_filter_optimist(
                !cutting_edges,
                &self.0,
            ))
    }

    fn connected_spinneys(&self) -> AHashSet<InternalSubGraph> {
        let cycles = self.0.all_cycles();

        let mut spinneys = AHashSet::new();

        for (ci, cj) in cycles.iter().tuple_combinations().map(|(a, b)| {
            (
                self.0.nesting_node_from_subgraph(a.clone()),
                self.0.nesting_node_from_subgraph(b.clone()),
            )
        }) {
            if !ci.strongly_disjoint(&cj) {
                let union = ci.internal_graph.union(&cj.internal_graph); //pairwise unions uffices up to 5 loops, then you need to do three unions.
                if self.dod(&union) >= 0 {
                    spinneys.insert(union);
                } else {
                    // println!("not dod >=0 spinney :{}", self.dod(&union));
                }
            }
        }

        for c in cycles {
            if self.dod(&c) >= 0 {
                spinneys.insert(c);
            }
        }

        spinneys.insert(self.0.empty_subgraph());
        spinneys
    }

    fn spinneys(&self) -> AHashSet<InternalSubGraph> {
        let cycles = self.0.all_cycles();

        let mut spinneys = AHashSet::new();

        for (ci, cj) in cycles.iter().tuple_combinations() {
            let union = ci.union(cj);

            if self.dod(&union) >= 0 {
                spinneys.insert(union);
            } else {
                // println!("not dod >=0 spinney :{}", self.dod(&union));
            }
        }

        for c in cycles {
            if self.dod(&c) >= 0 {
                spinneys.insert(c);
            }
        }

        spinneys.insert(self.0.empty_subgraph());

        println!("n_spinneys:{}", spinneys.len());
        spinneys
    }

    fn wood(&self) -> Wood {
        Wood::from_spinneys(self.spinneys(), self)
    }

    fn n_loops(&self, subgraph: &InternalSubGraph) -> usize {
        if let Some(loop_count) = subgraph.loopcount {
            // println!("found loop_count nloops: {}", loop_count);
            loop_count
        } else {
            self.0.cyclotomatic_number(subgraph)
        }
    }

    fn dod(&self, subgraph: &InternalSubGraph) -> i32 {
        let mut dod: i32 = 4 * self.n_loops(subgraph) as i32;
        // println!("nloops: {}", dod / 4);

        for e in self.0.iter_internal_edge_data(subgraph) {
            dod += e.dod;
        }

        for (_, n) in self.0.iter_node_data(subgraph) {
            dod += n.dod;
        }

        dod
    }

    fn numerator(&self, subgraph: &InternalSubGraph) -> SerializableAtom {
        let mut num = Atom::new_num(1);
        for (_, n) in self.0.iter_node_data(subgraph) {
            num = num * &n.num.0;
        }

        for e in self.0.iter_internal_edge_data(subgraph) {
            num = num * &e.num.0;
        }

        FunctionBuilder::new(State::get_symbol("num"))
            .add_arg(&num)
            .finish()
            .into()
    }

    fn denominator(&self, subgraph: &InternalSubGraph) -> SerializableAtom {
        let mut den = Atom::new_num(1);

        for e in self.0.iter_internal_edge_data(subgraph) {
            den = den * &e.den.0;
        }

        FunctionBuilder::new(State::get_symbol("den"))
            .add_arg(&den)
            .finish()
            .into()
    }

    fn dot(&self, subgraph: &InternalSubGraph) -> String {
        self.0.dot_impl(
            &subgraph.to_hairy_subgraph(&self.0),
            format!("{}", self.dod(subgraph)),
            &|e| Some(format!("label=\"{}\"", e.dod)),
            &|n| Some(format!("label=\"{}\"", n.dod)),
        )
    }

    fn t_op<I: Iterator<Item = InternalSubGraph>>(&self, mut subgraph_iter: I) -> Atom {
        if let Some(subgraph) = subgraph_iter.next() {
            let t = self.t_op(subgraph_iter);
            FunctionBuilder::new(State::get_symbol("Top"))
                .add_arg(&Atom::new_num(self.dod(&subgraph)))
                .add_arg(&t)
                .finish()
        } else {
            Atom::new_num(1)
        }
    }
}

pub struct Wood {
    poset: Poset<InternalSubGraph, Option<Top>>,
    additional_unions: SecondaryMap<PosetNode, Vec<PosetNode>>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TrieWood {
    trie: Trie<Top>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct UnfoldedWoodEl {
    expr: SerializableAtom,
    // graphs: Vec<PosetNode>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct UnfoldedWood {
    elements: Vec<UnfoldedWoodEl>,
}

impl UnfoldedWood {
    pub fn n_elements(&self) -> usize {
        self.elements.len()
    }

    pub fn show_structure(&self, wood: &Wood, graph: &UVGraph) -> String {
        let mut out = wood.show_graphs(graph);
        for e in &self.elements {
            out.push('\n');
            // for g in &e.graphs {
            //     out.push_str(&format!("{} ->", wood.poset.dot_id(*g)));
            // }
            out.push('\n');
            out.push_str(&format!("+{} \n", e.expr));
        }
        out
    }
}

impl Display for UnfoldedWood {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        for e in &self.elements {
            writeln!(f, "{}", e.expr)?;
            writeln!(f, " + ")?;
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct IntegrandExpr {
    dod: i32,
    num: SerializableAtom,
    den: SerializableAtom,
    add_arg: Option<SerializableAtom>,
}

impl IntegrandExpr {
    pub fn from_subgraph(
        subgraph: &InternalSubGraph,
        graph: &UVGraph,
        dod: i32,
        // add_arg: Option<IntegrandExpr>,
    ) -> Self {
        IntegrandExpr {
            num: graph.numerator(subgraph),
            den: graph.denominator(subgraph),
            add_arg: None,
            dod,
        }
    }

    pub fn from_subgraph_with_arg(
        subgraph: &InternalSubGraph,
        graph: &UVGraph,
        add_arg: Option<SerializableAtom>,
        dod: i32,
    ) -> Self {
        IntegrandExpr {
            num: graph.numerator(subgraph),
            den: graph.denominator(subgraph),
            add_arg,
            dod,
        }
    }

    pub fn top(&self) -> SerializableAtom {
        if let Some(add_arg) = &self.add_arg {
            FunctionBuilder::new(GS.top)
                .add_arg(&Atom::new_num(self.dod))
                .add_arg(&self.num.0)
                .add_arg(&self.den.0)
                .add_arg(&add_arg.0)
                .finish()
                // .ref_neg()
                .into()
        } else {
            FunctionBuilder::new(GS.top)
                .add_arg(&Atom::new_num(self.dod))
                .add_arg(&self.num.0)
                .add_arg(&self.den.0)
                .finish()
                // .ref_neg()
                .into()
        }
    }

    pub fn bare(&self) -> SerializableAtom {
        (&self.num.0 * &self.den.0).into()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Top {
    dod: i32,
    // graph_ref: PosetNode,
    pub t_arg: Option<IntegrandExpr>,
    pub contracted: IntegrandExpr,
    graph: TopoOrdered<InternalSubGraph>,
}

impl Top {
    pub fn root(graph: &UVGraph) -> Self {
        Top {
            dod: 0,
            // graph_ref: root_ref,
            t_arg: None,
            contracted: IntegrandExpr::from_subgraph(&graph.0.full_graph(), graph, 0),
            graph: TopoOrdered::new(graph.0.empty_subgraph(), 0),
        }
    }

    // pub fn new(&self,subgraph:SubGraph,graph:&UVGraph)->Self{

    // }

    pub fn from_subgraph(
        &self,
        subgraph: TopoOrdered<InternalSubGraph>,
        graph: &UVGraph,
        // graph_ref: PosetNode,
    ) -> Self {
        let reduced = self
            .graph
            .data
            .complement(&graph.0)
            .intersection(&subgraph.data);

        let dod = graph.dod(&subgraph.data);
        Top {
            // graph_ref,
            dod,
            t_arg: Some(IntegrandExpr::from_subgraph(&reduced, graph, dod)),
            contracted: IntegrandExpr::from_subgraph(
                &subgraph.data.complement(&graph.0),
                graph,
                dod,
            ),
            graph: subgraph,
        }
    }

    pub fn to_atom(&self, other: Option<AtomView>) -> SerializableAtom {
        if let Some(o) = other {
            self.to_fn_builder().add_arg(o)
        } else {
            self.to_fn_builder()
        }
        .finish()
        .into()
    }

    pub fn to_fn_builder(&self) -> FunctionBuilder {
        if let Some(integrand) = &self.t_arg {
            let num = &integrand.num.0;
            let den = &integrand.den.0;
            let dod = self.dod;
            FunctionBuilder::new(State::get_symbol("Top"))
                .add_arg(&Atom::new_num(dod))
                .add_arg(num)
                .add_arg(den)
        } else {
            FunctionBuilder::new(State::get_symbol("Top"))
        }
    }
}

#[allow(clippy::non_canonical_partial_ord_impl)]
impl PartialOrd for Top {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.graph.partial_cmp(&other.graph)
    }
}

impl Ord for Top {
    fn cmp(&self, other: &Self) -> Ordering {
        self.graph.cmp(&other.graph)
    }
}

impl TryFromIterator<Top, ()> for UnfoldedWoodEl {
    type Error = Report;

    fn try_from_iter<T>(iter: T) -> Result<Self, Self::Error>
    where
        Self: Sized,
        T: IntoIterator<Item = Top>,
    {
        let mut iter = iter.into_iter();
        let first_top = iter.next().ok_or(eyre!("empty iter.."))?;

        let mut expr = first_top.to_atom(None);
        let mut sign = 1;

        // let mut graphs = vec![first_top.graph_ref];

        let mut last = first_top.contracted;

        for t in iter {
            expr = t.to_atom(Some(expr.0.as_view()));
            // graphs.push(t.graph_ref);
            sign *= -1;
            last = t.contracted;
        }

        expr = (&expr.0 * &last.bare().0 * sign).into();

        Ok(UnfoldedWoodEl { expr })
    }
}

impl Wood {
    pub fn from_spinneys<I: IntoIterator<Item = InternalSubGraph>>(s: I, graph: &UVGraph) -> Self {
        let mut poset = Poset::from_iter(s.into_iter().map(|s| (s, None)));

        poset.invert();
        poset.compute_topological_order();

        let mut unions = SecondaryMap::new();

        for (i, sg) in poset.nodes.iter() {
            let cs = graph.0.connected_components(&sg.data);

            if cs.len() > 1 {
                let mut union = vec![];

                for &c in sg.parents.iter() {
                    let mut is_in = 0;
                    for comp in &cs {
                        if *comp == poset.nodes[c].data {
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
        let tree_root = dag.add_node(Approximation::new(
            self.poset.data(root).clone(),
            graph,
            root,
        ));
        search_front.push_front((root, tree_root));

        while let Some((node, parent)) = search_front.pop_front() {
            for c in &self.poset.nodes[node].children {
                if let Some(tagged_union) = unions.get_mut(*c) {
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
                                *c,
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
                        dag.add_node(Approximation::new(self.poset.data(*c).clone(), graph, *c));
                    dag.add_edge(parent, child);
                    search_front.push_front((*c, child));
                }
            }
        }
        tree_root
    }

    pub fn unfold_impl(&self, graph: &UVGraph) -> Forest {
        let mut dag: DAG<Approximation, DagNode, ()> = DAG::new();

        let root = self.poset.minimum().unwrap();

        let mut unions = SecondaryMap::new();

        for (p, u) in self.additional_unions.iter() {
            let union: Vec<(PosetNode, Option<DagNode>)> = u.iter().map(|i| (*i, None)).collect();
            unions.insert(p, Some(union));
        }

        let root = self.unfold_bfs(graph, &mut dag, &mut unions, root);

        Forest {
            dag,
            orig_shift: self.poset.shift(),
            all_graphs: AHashMap::new(),
        }
    }

    pub fn dot(&self, graph: &UVGraph) -> String {
        let shift = self.poset.shift();
        self.poset.to_dot_impl(&|n| {
            format!(
                "label={}, dod={}, n_edges = {},topo_order = {}",
                n.dot_id(shift),
                graph.dod(&n.data),
                graph.0.count_internal_edges(&n.data),
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
                graph.0.dot(&s.to_hairy_subgraph(&graph.0))
            );
        }
    }

    // pub fn unfold(&self, graph: &UVGraph) -> Atom {
    //     for (i, _) in self.coverset.nodes.iter().enumerate() {

    //         let mut t = graph.t_op(self.coverset.bfs_iter(v));
    //         t = FunctionBuilder::new(State::get_symbol("Top"))
    //             .add_arg(&Atom::new_num(v.numerator_rank()))
    //             .add_arg(&t)
    //             .finish();
    //     }
    //     // for s in self.coverset.bfs_iter(0){

    //     // }
    // }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
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
        symb!(&format!(
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
        fun!(GS.top, self.expr(bigger_graph))
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
    subgraph: InternalSubGraph,
    orig: PosetNode,
    dod: i32,
    pub approx_op: ApproxOp,
    pub simple_approx: Option<SimpleApprox>,
}

impl Approximation {
    pub fn simplify_notation(expr: &SerializableAtom) -> SerializableAtom {
        let replacements = [
            (
                Pattern::parse("ZERO").unwrap(),
                Atom::new_num(0).into_pattern().into(),
            ),
            (
                fun!(GS.num, Atom::new_num(1)).into_pattern(),
                Atom::new_num(1).into_pattern().into(),
            ),
            (
                fun!(GS.den, Atom::new_num(1)).into_pattern(),
                Atom::new_num(1).into_pattern().into(),
            ),
        ];

        let mut expr = expr.clone();

        let reps: Vec<_> = replacements
            .iter()
            .map(|(a, b)| Replacement::new(a, b))
            .collect();

        expr.replace_repeat_multiple(&reps);
        expr
    }

    pub fn new(spinney: InternalSubGraph, graph: &UVGraph, orig: PosetNode) -> Approximation {
        Approximation {
            dod: graph.dod(&spinney),
            orig,
            subgraph: spinney,
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
            ApproxOp::dependent(dependent, &self.subgraph, self.dod, graph),
            simple_approx,
        )
    }

    pub fn union(&self, dependents: &[&Self], graph: &UVGraph) -> (ApproxOp, SimpleApprox) {
        let simple_approx = SimpleApprox::union(
            self.subgraph.clone(),
            dependents.iter().map(|s| s.simple_approx.as_ref().unwrap()),
        );
        (ApproxOp::union(dependents).unwrap(), simple_approx)
    }

    pub fn expr(&self, graph: &UVGraph) -> Option<SerializableAtom> {
        let (t, s) = self.approx_op.expr()?;

        let contracted = s * IntegrandExpr::from_subgraph(
            &self.reduced_graph(&graph.0.full_graph()),
            graph,
            self.dod,
        )
        .bare()
        .0;

        Some(Self::simplify_notation(&(t.0 * contracted).into()))
    }

    pub fn simple_expr(&self, graph: &UVGraph) -> Option<SerializableAtom> {
        let simple_approx = self.simple_approx.as_ref()?;

        Some((simple_approx.sign * simple_approx.expr(&graph.0.full_filter())).into())
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

    pub fn expr(&self) -> Option<(SerializableAtom, Sign)> {
        match self {
            ApproxOp::NotComputed => None,
            ApproxOp::Union { t_args, sign, .. } => {
                let mut mul = Atom::new_num(1);
                for t in t_args {
                    mul = mul * t.top().0;
                }
                Some((mul.into(), *sign))
            }
            ApproxOp::Dependent { t_arg, sign, .. } => Some((t_arg.top().0.into(), *sign)),
            ApproxOp::Root => Some((Atom::new_num(1).into(), Sign::Positive)),
        }
    }

    pub fn dependent(
        dependent: &Approximation,
        subgraph: &InternalSubGraph,
        dod: i32,
        graph: &UVGraph,
    ) -> Self {
        let reduced = dependent.reduced_graph(subgraph);

        if let Some((expr, sign)) = dependent.approx_op.expr() {
            Self::Dependent {
                t_arg: IntegrandExpr::from_subgraph_with_arg(&reduced, graph, Some(expr), dod),
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
    all_graphs: AHashMap<InternalSubGraph, ()>,
    // root: DagNode,
    orig_shift: u64,
}

impl Forest {
    pub fn compute(&mut self, graph: &UVGraph) {
        let order = self.dag.compute_topological_order();

        for n in order {
            let current = &self.dag.nodes[n];

            let parents = current.parents.clone();

            let (approx, simple_approx) = if parents.is_empty() {
                let subgraph: InternalSubGraph = graph.0.empty_subgraph();
                self.all_graphs.insert(subgraph.clone(), ());
                (ApproxOp::Root, SimpleApprox::root(subgraph))
            } else if parents.len() == 1 {
                self.all_graphs.insert(
                    self.dag.nodes[parents[0]]
                        .data
                        .reduced_graph(&current.data.subgraph),
                    (),
                );
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
            out.push_str(&graph.0.dot_impl(
                &n.data.to_hairy_subgraph(&graph.0),
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
    pub fn unfold(&self, graph: &UVGraph) -> UnfoldedWood {
        let mut trie_builder = TrieBuilder::new();
        let _shift = self.poset.shift();

        for p in self.poset.bfs_paths().map(|s| {
            s.into_iter()
                .map(|k| {
                    (
                        self.poset.nodes.get(k).unwrap().to_topo_ordered().unwrap(),
                        k,
                    )
                })
                .collect::<Vec<_>>()
        }) {
            let root = Top::root(graph);

            let mut top_chain = vec![];

            for (s, k) in p {
                let new_top = top_chain.last().unwrap_or(&root).from_subgraph(s, graph);
                top_chain.push(new_top);
            }
            trie_builder.push(top_chain);
        }
        let triewood = TrieWood {
            trie: trie_builder.build(),
        };

        UnfoldedWood {
            elements: triewood.trie.iter::<UnfoldedWoodEl, ()>().collect(),
        }
    }

    // fn compute_t_operator(
    //     &self,
    //     subgraph: &SubGraph,
    //     graph: &UVGraph,
    //     t_results: &AHashMap<SubGraph, Atom>,
    //     predecessors: &[Vec<usize>],
    //     node_index: usize,
    // ) -> Atom {
    //     let t = FunctionBuilder::new(State::get_symbol("Top"))
    //         .add_arg(&Atom::new_num(graph.dod(subgraph)));

    //     let mut result = graph.numerator(subgraph).0 * graph.denominator(subgraph).0;
    //     for &pred_index in &predecessors[node_index] {
    //         let predecessor = &self.poset.nodes[pred_index];
    //         let pred_t_value = t_results.get(predecessor).unwrap();

    //         let complement = subgraph.complement().intersection(subgraph);

    //         let comp_num = graph.numerator(&complement).0;
    //         let comp_den = graph.denominator(&complement).0;

    //         result = result + comp_den * comp_num * pred_t_value;
    //     }

    //     -t.add_arg(&result).finish()
    // }
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
            dot.push_str(&format!("{} [{}];\n", node_id, label(node)));
            for &child in node.children.iter() {
                dot.push_str(&format!(
                    "{} -> {};\n",
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
