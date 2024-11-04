use std::{
    cmp::Ordering,
    collections::VecDeque,
    fmt::{self, Display, Formatter},
    hash::Hash,
    ops::{Deref, Index},
};

use ahash::{AHashMap, AHashSet};
use bitvec::{slice::IterOnes, vec::BitVec};
use color_eyre::Report;
use eyre::eyre;
use indexmap::IndexMap;
use itertools::Itertools;
use pathfinding::prelude::BfsReachable;
use rand::{rngs::SmallRng, Rng, SeedableRng};
use serde::{Deserialize, Serialize};
use spenso::symbolica_utils::SerializableAtom;
use symbolica::{
    atom::{Atom, AtomView, FunctionBuilder},
    state::State,
};
use trie_rs::{try_collect::TryFromIterator, Trie, TrieBuilder};

use crate::{
    graph::{
        half_edge::{
            subgraph::HedgeNode, subgraph::InternalSubGraph, HedgeGraph, HedgeGraphBuilder,
            NodeIndex, PowersetIterator,
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
        println!("name: {}", edge.name);
        print!("dod: {}", edge.dod());
        let index = (id * 100) as i32;
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
}

impl UVNode {
    pub fn from_vertex(vertex: &Vertex, graph: &BareGraph) -> Self {
        println!("name: {}", vertex.name);
        print!("dod: {}", vertex.dod());
        if let Some(colorless) = vertex.contracted_vertex_rule(graph) {
            UVNode {
                dod: vertex.dod() as i32,
                num: normalise_complex(&colorless).into(),
            }
        } else {
            UVNode {
                dod: vertex.dod() as i32,
                num: Atom::new_num(1).into(),
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

    pub fn node_id(&self, g: &BareGraph, id: usize) -> HalfEdgeNode {
        let e_id = g.vertices[id].edges.first().unwrap();
        self.0
            .involution
            .get_node_id(self.half_edge_id(g, *e_id))
            .clone()
    }

    pub fn cycle_basis_from_lmb(&self, lmb: &LoopMomentumBasis) -> Vec<HalfEdgeNode> {
        let mut cycle_basis = Vec::new();
        let cut_edges = self.0.paired_filter_from_pos(&lmb.basis);

        for v in lmb.basis.iter() {
            let loop_edge = self.0.paired_filter_from_pos(&[*v]);

            let mut hairy_loop = self
                .0
                .nesting_node_from_subgraph(SubGraph::from(!cut_edges.clone() | &loop_edge));

            self.0.cut_branches(&mut hairy_loop);
            cycle_basis.push(hairy_loop);
        }
        cycle_basis
    }

    fn spanning_forest_from_lmb(&self, lmb: LoopMomentumBasis) -> HalfEdgeNode {
        let cutting_edges = self.0.paired_filter_from_pos(&lmb.basis);

        self.0
            .nesting_node_from_subgraph(SubGraph::from(!cutting_edges))
    }

    fn numerator_rank(&self, subgraph: &SubGraph) -> isize {
        let mut rank = 0;

        for (id, n) in self.0.iter_node_data(subgraph) {
            rank += n.numerator_rank();
            for e in id.externalhedges.iter_ones() {
                rank += self.0.get_edge_data(e).numerator_rank();
            }
        }

        rank
    }

    fn connected_spinneys(&self) -> AHashSet<SubGraph> {
        let cycles = self.0.all_cycles();

        let mut spinneys = AHashSet::new();

        for (ci, cj) in cycles.iter().tuple_combinations() {
            if !ci.strongly_disjoint(cj) {
                let union = ci.internal_graph.union(&cj.internal_graph); //pairwise unions suffices up to 5 loops, then you need to do three unions.
                if self.dod(&union) >= 0 {
                    spinneys.insert(union);
                } else {
                    // println!("not dod >=0 spinney :{}", self.dod(&union));
                }
            }
        }

        for c in cycles {
            if self.dod(&c.internal_graph) >= 0 {
                spinneys.insert(c.internal_graph);
            }
        }

        spinneys.insert(self.0.empty_filter().into());
        spinneys
    }

    fn spinneys(&self) -> AHashSet<SubGraph> {
        let cycles = self.0.all_cycles();

        let mut spinneys = AHashSet::new();

        for (ci, cj) in cycles.iter().tuple_combinations() {
            let union = ci.internal_graph.union(&cj.internal_graph);
            if self.dod(&union) >= 0 {
                spinneys.insert(union);
            } else {
                // println!("not dod >=0 spinney :{}", self.dod(&union));
            }
        }

        for c in cycles {
            if self.dod(&c.internal_graph) >= 0 {
                spinneys.insert(c.internal_graph);
            }
        }

        spinneys.insert(self.0.empty_filter().into());
        spinneys
    }

    fn wood(&self) -> Wood {
        Wood::from_spinneys(self.spinneys(), self)
    }

    fn n_loops(&self, subgraph: &SubGraph) -> usize {
        if let Some(loop_count) = subgraph.loopcount {
            // println!("found loop_count nloops: {}", loop_count);
            loop_count
        } else {
            self.0.cyclotomatic_number(subgraph)
        }
    }

    fn dod(&self, subgraph: &SubGraph) -> i32 {
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

    fn numerator(&self, subgraph: &SubGraph) -> SerializableAtom {
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

    fn denominator(&self, subgraph: &SubGraph) -> SerializableAtom {
        let mut den = Atom::new_num(1);

        for e in self.0.iter_internal_edge_data(subgraph) {
            den = den * &e.den.0;
        }

        FunctionBuilder::new(State::get_symbol("den"))
            .add_arg(&den)
            .finish()
            .into()
    }

    fn dot(&self, subgraph: &SubGraph) -> String {
        self.0.dot_impl(
            &subgraph.to_nesting_node(&self.0),
            format!("{}", self.dod(subgraph)),
            &|e| format!("label=\"{}\"", e.dod),
            &|n| format!("label=\"{}\"", n.dod),
        )
    }

    fn t_op<I: Iterator<Item = SubGraph>>(&self, mut subgraph_iter: I) -> Atom {
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
    poset: Poset<SubGraph, Option<Top>>,
    additional_unions: SecondaryMap<NodeRef, Vec<NodeRef>>,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TrieWood {
    trie: Trie<Top>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct UnfoldedWoodEl {
    expr: SerializableAtom,
    graphs: Vec<NodeRef>,
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
            for g in &e.graphs {
                out.push_str(&format!("{} ->", wood.poset.dot_id(*g)));
            }
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
    num: SerializableAtom,
    den: SerializableAtom,
}

impl IntegrandExpr {
    pub fn from_subgraph(
        subgraph: &SubGraph,
        graph: &UVGraph,
        // add_arg: Option<IntegrandExpr>,
    ) -> Self {
        IntegrandExpr {
            num: graph.numerator(subgraph),
            den: graph.denominator(subgraph),
        }
    }

    pub fn top(&self, dod: i32) -> SerializableAtom {
        FunctionBuilder::new(State::get_symbol("Top"))
            .add_arg(&Atom::new_num(dod))
            .add_arg(&self.num.0)
            .add_arg(&self.den.0)
            .finish()
            .into()
    }

    pub fn bare(&self) -> SerializableAtom {
        (&self.num.0 * &self.den.0).into()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Top {
    dod: i32,
    graph_ref: NodeRef,
    pub t_arg: Option<IntegrandExpr>,
    pub contracted: IntegrandExpr,
    graph: TopoOrdered<SubGraph>,
}

impl Top {
    pub fn root(graph: &UVGraph, root_ref: NodeRef) -> Self {
        Top {
            dod: 0,
            graph_ref: root_ref,
            t_arg: None,
            contracted: IntegrandExpr::from_subgraph(&graph.0.full_graph(), graph),
            graph: TopoOrdered::new(graph.0.full_graph().complement(), 0),
        }
    }

    pub fn from_subgraph(
        &self,
        subgraph: TopoOrdered<SubGraph>,
        graph: &UVGraph,
        graph_ref: NodeRef,
    ) -> Self {
        let reduced = self.graph.data.complement().intersection(&subgraph.data);

        Top {
            graph_ref,
            dod: graph.dod(&subgraph.data),
            t_arg: Some(IntegrandExpr::from_subgraph(&reduced, graph)),
            contracted: IntegrandExpr::from_subgraph(&subgraph.data.complement(), graph),
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

        let mut graphs = vec![first_top.graph_ref];

        let mut last = first_top.contracted;

        for t in iter {
            expr = t.to_atom(Some(expr.0.as_view()));
            graphs.push(t.graph_ref);
            sign *= -1;
            last = t.contracted;
        }

        expr = (&expr.0 * &last.bare().0 * sign).into();

        Ok(UnfoldedWoodEl { expr, graphs })
    }
}

impl Wood {
    pub fn from_spinneys<I: IntoIterator<Item = SubGraph>>(s: I, graph: &UVGraph) -> Self {
        let mut poset = Poset::from_iter(s.into_iter().map(|s| (s, None)));

        poset.invert();
        poset.compute_topological_order();

        let mut unions = SecondaryMap::new();

        for (i, sg) in poset.nodes.iter() {
            let cs = graph.0.connected_components(&sg.data);

            if cs.len() > 1 {
                let mut union = vec![];

                for &c in sg.children.iter() {
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
                graph.0.dot(&s.to_nesting_node(&graph.0))
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

impl Wood {
    pub fn show_graphs(&self, graph: &UVGraph) -> String {
        let mut out = String::new();
        out.push_str("Poset structure:\n");
        out.push_str(&self.poset.dot_structure());

        out.push_str("Graphs:\n");
        for (k, n) in self.poset.nodes.iter() {
            out.push_str(&graph.0.dot_impl(
                &n.data.to_nesting_node(&graph.0),
                format!(
                    "dod={};nodeid ={};\n",
                    graph.dod(&n.data),
                    self.poset.dot_id(k)
                ),
                &|e| format!("dod={}", e.dod),
                &|n| format!("dod={}", n.dod),
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
            let root = Top::root(graph, self.poset.minimum().unwrap());

            let mut top_chain = vec![];

            for (s, k) in p {
                let new_top = top_chain.last().unwrap_or(&root).from_subgraph(s, graph, k);
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

use slotmap::{new_key_type, SecondaryMap, SlotMap};
use std::collections::HashSet;

// Define a new key type for the Poset
new_key_type! {
    pub struct NodeRef;
}

/// Trait to define DOT attributes for node data.
pub trait DotAttrs {
    fn dot_attrs(&self) -> String;
}

/// A node in the poset, storing generic data, edges to child nodes, and references to parent nodes.
#[derive(Debug, Eq)]
pub struct SlotNode<T> {
    data: T,
    order: Option<u64>,
    id: NodeRef,
    parents: Vec<NodeRef>,  // References to parent nodes by key
    children: Vec<NodeRef>, // Edges to child nodes by key
}

impl<T: Hash> Hash for SlotNode<T> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.data.hash(state);
        self.id.hash(state);
    }
}

impl<T: PartialEq> PartialEq for SlotNode<T> {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data && self.id == other.id
    }
}

impl<T: PartialEq> PartialOrd for SlotNode<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.order.unwrap().partial_cmp(&other.order.unwrap())
    }
}

impl<T: PartialEq + Eq> Ord for SlotNode<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.order.unwrap().cmp(&other.order.unwrap())
    }
}

impl<T> SlotNode<T> {
    pub fn dot_id(&self, shift: u64) -> String {
        format!("node{}", base62::encode(self.id() - shift))
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
        self.id.0.as_ffi()
    }

    pub fn add_parent(&mut self, parent: NodeRef) {
        self.parents.push(parent);
    }

    pub fn add_child(&mut self, child: NodeRef) {
        self.children.push(child);
    }

    pub fn remove_child(&mut self, child: NodeRef) {
        self.children.retain(|&c| c != child);
    }

    pub fn remove_parent(&mut self, parent: NodeRef) {
        self.parents.retain(|&c| c != parent);
    }

    pub fn is_parent_of(&self, child: NodeRef) -> bool {
        self.children.contains(&child)
    }

    pub fn new(data: T, id: NodeRef) -> Self {
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
pub struct DAG<T, D = ()> {
    pub nodes: SlotMap<NodeRef, SlotNode<T>>,
    associated_data: SecondaryMap<NodeRef, D>,
}

#[allow(dead_code)]
pub struct TreeNode<T> {
    data: T,
    order: Option<u64>,
    id: NodeRef,
    parent: NodeRef,
    children: Vec<NodeRef>, // Edges to child nodes by key
}

#[allow(dead_code)]
pub struct Tree<T, D = ()> {
    nodes: SlotMap<NodeRef, TreeNode<T>>,
    associated_data: SecondaryMap<NodeRef, D>,
}

pub type Poset<T, D> = DAG<T, D>;
pub type HasseDiagram<T, D> = DAG<T, D>;

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

impl<T, D> Poset<T, D> {
    pub fn new() -> Self {
        Poset {
            nodes: SlotMap::with_key(),
            associated_data: SecondaryMap::new(),
        }
    }

    pub fn dot_id(&self, key: NodeRef) -> String {
        self.nodes.get(key).unwrap().dot_id(self.shift())
    }

    pub fn node_values(&self) -> impl Iterator<Item = &T> {
        self.nodes.values().map(|node| &node.data)
    }

    pub fn add_edge(&mut self, from: NodeRef, to: NodeRef) {
        let from_node = self.nodes.get_mut(from).unwrap();
        from_node.add_child(to);
        let to_node = self.nodes.get_mut(to).unwrap();
        to_node.add_parent(from);
    }

    pub fn add_edge_if_new(&mut self, from: NodeRef, to: NodeRef) {
        if !self.nodes.get(from).unwrap().is_parent_of(to) {
            self.add_edge(from, to);
        }
    }

    pub fn remove_edge(&mut self, from: NodeRef, to: NodeRef) {
        let from_node = self.nodes.get_mut(from).unwrap();
        from_node.remove_child(to);
        let to_node = self.nodes.get_mut(to).unwrap();
        to_node.remove_parent(from);
    }

    pub fn bfs_reach<'a>(
        &'a self,
        start: &'a NodeRef,
    ) -> BfsReachable<&'a NodeRef, impl FnMut(&'a &'a NodeRef) -> &'a [NodeRef]> {
        pathfinding::directed::bfs::bfs_reach(start, |&s| self.succesors(*s))
    }

    pub fn maximum(&self) -> Option<NodeRef>
    where
        T: Eq,
    {
        Some(self.nodes.values().max()?.id)
    }

    pub fn minimum(&self) -> Option<NodeRef>
    where
        T: Eq,
    {
        Some(self.nodes.values().min()?.id)
    }

    pub fn succesors(&self, node_key: NodeRef) -> &[NodeRef] {
        &self.nodes.get(node_key).unwrap().children
    }

    /// Returns an iterator over all paths starting from the root node, traversed in BFS order.
    pub fn bfs_paths(&self) -> BfsPaths<T>
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

    pub fn bfs_paths_inv(&self) -> BfsPaths<T> {
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

    pub fn data(&self, key: NodeRef) -> &T {
        &self.nodes.get(key).unwrap().data
    }

    fn shift(&self) -> u64 {
        self.nodes.iter().next().unwrap().1.id()
    }

    pub fn to_dot(&self, label: &impl Fn(&T) -> String) -> String {
        self.to_dot_impl(&|node| label(&node.data))
    }

    fn to_dot_impl(&self, label: &impl Fn(&SlotNode<T>) -> String) -> String {
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

    pub fn poset_family(&self, data: &T) -> [Vec<NodeRef>; 2]
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

    pub fn poset_push(&mut self, data: T, associated_data: D, flip: bool) -> NodeRef
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
    fn update_transitive_closure(&mut self, new_node_key: NodeRef) {
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

    /// Returns all descendants of the node, used for propagating transitive relations.
    fn get_all_descendants(&self, node_key: NodeRef) -> HashSet<NodeRef> {
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
    fn get_all_ancestors(&self, node_key: NodeRef) -> HashSet<NodeRef> {
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

        HasseDiagram {
            nodes: self.nodes,
            associated_data: self.associated_data,
        }
    }

    pub fn transitive_edges(&self, a: NodeRef) -> Vec<(NodeRef, NodeRef)> {
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

    pub fn in_degree(&self, key: NodeRef) -> usize {
        self.nodes.get(key).unwrap().in_degree()
    }

    pub fn compute_topological_order(&mut self) {
        // Initialize queue with nodes having in-degree zero
        let mut queue = VecDeque::new(); //S in the wikipedia article

        let mut indegrees = SecondaryMap::new();

        for (key, node) in self.nodes.iter() {
            indegrees.insert(key, node.in_degree());
            if node.in_degree() == 0 {
                queue.push_back(key);
            }
        }

        let mut order = 0;
        while let Some(node_key) = queue.pop_front() {
            // Assign the order number to the node
            if let Some(node) = self.nodes.get_mut(node_key) {
                node.order = Some(order);
                order += 1;
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
        if order as usize != self.nodes.len() {
            panic!("The graph contains a cycle!");
        }
    }
}

impl<T, D> Default for Poset<T, D> {
    fn default() -> Self {
        Self::new()
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
pub struct BfsPaths<'a, T> {
    queue: VecDeque<Vec<NodeRef>>,
    visited: HashSet<Vec<NodeRef>>,
    nodes: &'a SlotMap<NodeRef, SlotNode<T>>,
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct SubGraphChainLink {
    subgraph: SubGraph,
}

impl Deref for SubGraphChainLink {
    type Target = SubGraph;

    fn deref(&self) -> &Self::Target {
        &self.subgraph
    }
}

impl From<SubGraph> for SubGraphChainLink {
    fn from(subgraph: SubGraph) -> Self {
        SubGraphChainLink { subgraph }
    }
}

#[allow(clippy::non_canonical_partial_ord_impl)]
impl PartialOrd for SubGraphChainLink {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.subgraph.partial_cmp(&other.subgraph)
    }
}

impl Ord for SubGraphChainLink {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl<'a, T> Iterator for BfsPaths<'a, T> {
    type Item = Vec<NodeRef>;

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
