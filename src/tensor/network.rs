use ahash::AHashMap;

use serde::{Deserialize, Serialize};
use slotmap::{new_key_type, DenseSlotMap, Key, SecondaryMap};
use symbolica::{
    domains::float::Complex,
    representations::{Atom, AtomView, Symbol},
    state::State,
};

use self::parametric::{MixedTensor, MixedTensors};
use self::structure::HistoryStructure;

use super::{
    arithmetic::ScalarMul, parametric, structure, Contract, DataTensor, HasName, Shadowable, Slot,
    TensorStructure, TracksCount,
};
use smartstring::alias::String;
use std::fmt::{Debug, Display};

new_key_type! {
    pub struct NodeId;
    pub struct HedgeId;
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HalfEdgeGraph<N, E> {
    pub edges: DenseSlotMap<HedgeId, E>,
    pub involution: SecondaryMap<HedgeId, HedgeId>,
    pub neighbors: SecondaryMap<HedgeId, HedgeId>,
    pub nodes: DenseSlotMap<NodeId, N>,
    pub nodemap: SecondaryMap<HedgeId, NodeId>,
    pub reverse_nodemap: SecondaryMap<NodeId, HedgeId>,
}

struct IncidentIterator<'a> {
    neighbors: &'a SecondaryMap<HedgeId, HedgeId>,
    current: Option<HedgeId>,
    start: HedgeId,
}

impl<'a> Iterator for IncidentIterator<'a> {
    type Item = HedgeId;
    fn next(&mut self) -> Option<HedgeId> {
        let current = self.current?;

        self.current = Some(self.neighbors[current]);

        if self.current == Some(self.start) {
            self.current = None;
        }

        Some(current)
    }
}

impl<'a> IncidentIterator<'a> {
    fn new<N, E>(graph: &'a HalfEdgeGraph<N, E>, initial: HedgeId) -> Self {
        IncidentIterator {
            neighbors: &graph.neighbors,
            current: Some(initial),
            start: initial,
        }
    }
}
#[allow(dead_code)]
impl<N, E> HalfEdgeGraph<N, E> {
    fn new() -> Self {
        HalfEdgeGraph {
            involution: SecondaryMap::new(),
            nodemap: SecondaryMap::new(),
            neighbors: SecondaryMap::new(),
            reverse_nodemap: SecondaryMap::new(),
            nodes: DenseSlotMap::with_key(),
            edges: DenseSlotMap::with_key(),
        }
    }

    fn dot(&self) -> std::string::String {
        let mut out = "graph {".to_string();
        out.push_str("  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\";");

        // for (i, n) in &self.nodes {
        //     out.push_str(&format!("\n {}", i.data().as_ffi()));
        // }
        for (i, _) in &self.neighbors {
            match i.cmp(&self.involution[i]) {
                std::cmp::Ordering::Greater => {
                    out.push_str(&format!(
                        "\n {} -- {}",
                        self.nodemap[i].data().as_ffi(),
                        self.nodemap[self.involution[i]].data().as_ffi()
                    ));
                }
                std::cmp::Ordering::Equal => {
                    out.push_str(&format!(
                        "ext{} [shape=none, label=\"\"];",
                        i.data().as_ffi()
                    ));
                    out.push_str(&format!(
                        "\n {} -- ext{};",
                        self.nodemap[i].data().as_ffi(),
                        i.data().as_ffi()
                    ));
                }
                _ => {}
            }
        }

        out += "}";
        out
    }

    fn add_node(&mut self, data: N) -> NodeId {
        self.nodes.insert(data)
    }

    fn node_indices(&self) -> slotmap::dense::Keys<'_, NodeId, N> {
        self.nodes.keys()
    }

    /// Add a node with a list of edget with associated data. Matches edges by equality.
    fn add_node_with_edges(&mut self, data: N, edges: &[E]) -> NodeId
    where
        E: Eq + Clone,
    {
        let idx = self.add_node(data);
        for e in edges {
            let mut found_match = false;
            for (i, other_e) in &self.edges {
                if *e == *other_e && self.involution[i] == i {
                    found_match = true;
                    let eid = self.edges.insert(e.clone());
                    self.involution.insert(eid, i);
                    self.involution.insert(i, eid);
                    self.nodemap.insert(eid, idx);
                    if let Some(prev_eid) = self.reverse_nodemap.insert(idx, eid) {
                        let next_eid = self.neighbors.insert(prev_eid, eid).unwrap();
                        self.neighbors.insert(eid, next_eid);
                    } else {
                        self.neighbors.insert(eid, eid);
                    }
                    break;
                }
            }
            if !found_match {
                let eid = self.edges.insert(e.clone());
                self.involution.insert(eid, eid);
                self.nodemap.insert(eid, idx);
                if let Some(prev_eid) = self.reverse_nodemap.insert(idx, eid) {
                    let next_eid = self.neighbors.insert(prev_eid, eid).unwrap();
                    self.neighbors.insert(eid, next_eid);
                } else {
                    self.neighbors.insert(eid, eid);
                }
            }
        }

        idx
    }

    pub fn validate_neighbors(&self) -> bool {
        for (i, n) in &self.reverse_nodemap {
            for j in IncidentIterator::new(self, *n) {
                if self.nodemap[j] != i {
                    return false;
                }
            }
        }
        true
    }

    fn node_labels(&self) -> String
    where
        N: Display,
    {
        let mut out = String::new();
        for (i, n) in &self.nodes {
            out.push_str(&format!("{}[label= \"{}\"]\n", i.data().as_ffi(), n));
        }
        out
    }

    #[allow(clippy::too_many_lines)]
    fn merge_nodes(&mut self, a: NodeId, b: NodeId, data: N) -> NodeId {
        let c = self.nodes.insert(data);

        // New initial edge for reverse_nodemap, that does not link to b
        // if none is found, all incident edges are link to b and must be removed from the neighbors list
        let mut new_initial_a = self
            .edges_incident(a)
            .find(|x| self.nodemap[self.involution[*x]] != b && self.involution[*x] != *x);

        if new_initial_a.is_none() {
            new_initial_a = self
                .edges_incident(a)
                .find(|x| self.nodemap[self.involution[*x]] != b);
        }

        let mut a_edge = new_initial_a;

        if a_edge.is_none() {
            // all edges link to b, and must be removed
            let initial = self.reverse_nodemap[a];
            let mut current = Some(initial);
            loop {
                if current.is_none() {
                    break;
                }
                let next = self.neighbors.remove(current.unwrap());

                if next == Some(initial) {
                    current = None;
                } else {
                    current = next;
                }
            }
        } else {
            loop {
                let mut next = self.neighbors[a_edge.unwrap()];

                while self.nodemap[self.involution[next]] == b {
                    next = self.neighbors.remove(next).unwrap();
                }

                self.nodemap.insert(a_edge.unwrap(), c);
                self.neighbors.insert(a_edge.unwrap(), next);

                if new_initial_a == Some(next) {
                    break;
                }

                a_edge = Some(next);
            }
        }

        let mut new_initial_b = self
            .edges_incident(b)
            .find(|x| self.nodemap[self.involution[*x]] != a && self.involution[*x] != *x);

        if new_initial_b.is_none() {
            new_initial_b = self
                .edges_incident(b)
                .find(|x| self.nodemap[self.involution[*x]] != a);
        }
        let mut b_edge = new_initial_b;

        if b_edge.is_none() {
            let initial = self.reverse_nodemap[b];
            let mut current = Some(initial);
            loop {
                if current.is_none() {
                    break;
                }
                let next = self.neighbors.remove(current.unwrap());

                if next == Some(initial) {
                    current = None;
                } else {
                    current = next;
                }
            }
        } else {
            loop {
                let mut next = self.neighbors[b_edge.unwrap()];

                while self.nodemap[self.involution[next]] == a {
                    next = self.neighbors.remove(next).unwrap();
                }

                self.nodemap.insert(b_edge.unwrap(), c);
                self.neighbors.insert(b_edge.unwrap(), next);

                if new_initial_b == Some(next) {
                    break;
                }

                b_edge = Some(next);
            }
        }

        match (new_initial_a, new_initial_b) {
            (Some(new_edge_a), Some(new_edge_b)) => {
                self.reverse_nodemap.insert(c, new_edge_a);
                self.reverse_nodemap.remove(a);
                self.reverse_nodemap.remove(b);
                let old_neig = self.neighbors.insert(new_edge_a, new_edge_b).unwrap();
                self.neighbors.insert(b_edge.unwrap(), old_neig);
            }
            (Some(new_edge_a), None) => {
                self.reverse_nodemap.insert(c, new_edge_a);
                self.reverse_nodemap.remove(a);
                self.reverse_nodemap.remove(b);
            }
            (None, Some(new_edge_b)) => {
                self.reverse_nodemap.insert(c, new_edge_b);
                self.reverse_nodemap.remove(a);
                self.reverse_nodemap.remove(b);
            }
            (None, None) => {
                self.reverse_nodemap.remove(b);
                self.reverse_nodemap.remove(a);
            }
        }
        self.nodes.remove(a);
        self.nodes.remove(b);
        c
    }

    /// Add an internal edge between two nodes.
    fn add_edge(&mut self, a: NodeId, b: NodeId, data: E) -> HedgeId
    where
        E: Clone,
    {
        let hedge_id_a = self.edges.insert(data.clone());
        let hedge_id_b = self.edges.insert(data);
        self.involution.insert(hedge_id_a, hedge_id_b);
        self.involution.insert(hedge_id_b, hedge_id_a);
        self.nodemap.insert(hedge_id_a, a);
        if let Some(prev_eid) = self.reverse_nodemap.insert(a, hedge_id_a) {
            let next_eid = self.neighbors.insert(prev_eid, hedge_id_a).unwrap();
            self.neighbors.insert(hedge_id_a, next_eid).unwrap();
        } else {
            self.neighbors.insert(hedge_id_a, hedge_id_a);
        }
        self.nodemap.insert(hedge_id_b, b);
        if let Some(prev_eid) = self.reverse_nodemap.insert(b, hedge_id_b) {
            let next_eid = self.neighbors.insert(prev_eid, hedge_id_b).unwrap();
            self.neighbors.insert(hedge_id_b, next_eid).unwrap();
        } else {
            self.neighbors.insert(hedge_id_b, hedge_id_b);
        }
        hedge_id_a
    }

    /// Add external, as a fixed point involution half edge.
    fn add_external(&mut self, a: NodeId, data: E) -> HedgeId {
        let id = self.edges.insert(data);
        self.involution.insert(id, id);
        self.nodemap.insert(id, a);
        if let Some(prev_eid) = self.reverse_nodemap.insert(a, id) {
            let next_eid = self.neighbors.insert(prev_eid, id).unwrap();
            self.neighbors.insert(id, next_eid).unwrap();
        } else {
            self.neighbors.insert(id, id);
        }
        id
    }

    fn edges_incident(&self, node: NodeId) -> impl Iterator<Item = HedgeId> + '_ {
        IncidentIterator::new(self, self.reverse_nodemap[node])
    }

    fn edges_between(&self, a: NodeId, b: NodeId) -> impl Iterator<Item = HedgeId> + '_ {
        self.edges_incident(a)
            .filter(move |&i| self.nodemap[self.involution[i]] == b)
    }

    fn internal_edges_incident(&self, node: NodeId) -> impl Iterator<Item = HedgeId> + '_ {
        self.edges_incident(node)
            .filter(move |&i| self.nodemap[self.involution[i]] != node)
    }

    fn external_edges_incident(&self, node: NodeId) -> impl Iterator<Item = HedgeId> + '_ {
        self.edges_incident(node)
            .filter(move |&i| self.nodemap[self.involution[i]] == node)
    }

    fn degree(&self, node: NodeId) -> usize {
        self.edges_incident(node).collect::<Vec<_>>().len()
    }

    fn neighbors(&self, node: NodeId) -> impl Iterator<Item = NodeId> + '_ {
        self.edges_incident(node)
            .map(move |i| self.nodemap[self.involution[i]])
    }

    // fn map_nodes<F, U>(&self, f: F) -> HalfEdgeGraph<U, E>
    // where
    //     F: Fn(&N) -> U,
    //     E: Clone,
    // {
    //     let edges = self.edges.clone();
    //     let involution = self.involution.clone();

    //     let mut nodes = DenseSlotMap::with_key();
    //     let mut nodemap = SecondaryMap::new();

    //     for n in &self.nodes {
    //         let nid = nodes.insert(f(n.1));
    //         for e in self.edges_incident(n.0) {
    //             nodemap.insert(e, nid);
    //         }
    //     }

    //     HalfEdgeGraph {
    //         edges,
    //         involution,
    //         nodes,
    //         nodemap,
    //     }
    // }
}

#[test]
fn merge() {
    let mut graph = HalfEdgeGraph::new();
    let a = graph.add_node_with_edges(1, &[1, 2, 3, 4, 5]);
    let b = graph.add_node_with_edges(2, &[1, 2, 6, 7, 8]);
    let c = graph.add_node_with_edges(4, &[4, 6, 9, 10, 11]);

    println!("{}", graph.dot());
    println!("{}", graph.degree(a));
    println!("{}", graph.degree(b));

    for (i, n) in &graph.neighbors {
        println!("{} {}", graph.edges[i], graph.edges[*n]);
    }

    let d = graph.merge_nodes(a, b, 3);

    // for (i, n) in &graph.neighbors {
    //     println!("{} {}", graph.edges[i], graph.edges[*n]);
    // }

    println!("{}", graph.dot());
    println!("{}", graph.degree(c));
    println!("{}", graph.neighbors.len());

    let e = graph.merge_nodes(c, d, 5);

    println!("{}", graph.dot());
    println!("{}", graph.degree(e));
    println!("{}", graph.neighbors.len());

    let mut graph = HalfEdgeGraph::new();
    let a = graph.add_node_with_edges("a", &[10, 2, 3]);
    let b = graph.add_node_with_edges("b", &[20, 3, 4]);
    let c = graph.add_node_with_edges("c", &[30, 4, 2]);
    let d = graph.add_node_with_edges("d", &[20]);
    let e = graph.add_node_with_edges("e", &[30]);

    println!("Test {}", graph.dot());
    println!("{}", graph.degree(a));
    println!("{}", graph.degree(b));

    for (i, n) in &graph.neighbors {
        println!("{} {}", graph.edges[i], graph.edges[*n]);
    }

    let d = graph.merge_nodes(d, b, "bd");

    // for (i, n) in &graph.neighbors {
    //     println!("{} {}", graph.edges[i], graph.edges[*n]);
    // }

    println!("{}", graph.degree(c));
    println!("{}", graph.neighbors.len());

    println!("{}", graph.dot());

    let e = graph.merge_nodes(c, e, "ce");

    if graph.validate_neighbors() {
        println!("valid");
    } else {
        println!("invalid");
    }

    println!("{}", graph.dot());
    let f = graph.merge_nodes(d, e, "de");

    if graph.validate_neighbors() {
        println!("valid");
    } else {
        println!("invalid");
    }

    println!("{}", graph.dot());
    println!("{}", graph.node_labels());
    println!("{}", graph.degree(a));
    println!("{}", graph.neighbors.len());

    let g = graph.merge_nodes(a, f, "af");

    if graph.validate_neighbors() {
        println!("valid");
    } else {
        println!("invalid");
    }

    println!("{}", graph.dot());
    println!("{}", graph.neighbors.len());
    println!("{}", graph.degree(g));

    // println!("{}", graph.degree(b));
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TensorNetwork<T> {
    pub graph: HalfEdgeGraph<T, Slot>,
    pub params: Vec<T>,
}

impl<T> TensorNetwork<T> {
    fn edge_to_min_degree_node(&self) -> Option<HedgeId> {
        let mut neighs = self.graph.reverse_nodemap.clone();
        if neighs.is_empty() {
            return None;
        }

        loop {
            let mut all_ext = true;
            for (node, initial) in &mut neighs {
                *initial = self.graph.neighbors[*initial];
                let start = self.graph.reverse_nodemap[node];

                if self.graph.involution[start] != start {
                    all_ext = false;
                    if *initial == start {
                        return Some(start);
                    }
                }
            }
            if all_ext {
                return None;
            }
        }
    }

    pub fn to_vec(&self) -> Vec<&T> {
        self.graph.nodes.values().collect()
    }
}

impl<N> TensorNetwork<MixedTensor<N>>
where
    N: Debug + TensorStructure,
{
    pub fn to_symbolic_tensor_vec(mut self) -> Vec<DataTensor<Atom, N>> {
        self.graph
            .nodes
            .drain()
            .filter(|(_, n)| n.is_symbolic())
            .map(|(_, n)| n.try_into_symbolic().unwrap())
            .collect()
    }

    // pub fn evaluate<'a, D>(
    //     &'a self,
    //     const_map: &AHashMap<AtomView<'a>, D>,
    // ) -> TensorNetwork<DataTensor<D, N>>
    // where
    //     D: Clone,
    //     N: Clone + TensorStructure,
    // {
    //     // let mut evaluated_net = TensorNetwork::new();
    //     // for (id,t) in &self.graph.nodes {

    //     //     let evaluated_tensor = match t{
    //     //         MixedTensor::Complex(t)=> t.into(),
    //     //         MixedTensor::Symbolic(t)=> t.evaluate(const_map),
    //     //         MixedTensor::Float(t)=> t.into(),
    //     //     }
    //     //     evaluated_net.push(evaluated_tensor);
    //     // }

    //     // evaluated_net
    // }

    pub fn evaluate_float<'a>(&'a mut self, const_map: &AHashMap<AtomView<'a>, f64>)
    where
        N: Clone,
    {
        for (_, n) in &mut self.graph.nodes {
            n.evaluate_float(const_map);
        }
    }

    pub fn evaluate_complex<'a>(&'a mut self, const_map: &AHashMap<AtomView<'a>, Complex<f64>>)
    where
        N: Clone,
    {
        for (_, n) in &mut self.graph.nodes {
            n.evaluate_complex(const_map);
        }
    }
}

impl<T> From<Vec<T>> for TensorNetwork<T>
where
    T: TensorStructure,
{
    fn from(tensors: Vec<T>) -> Self {
        TensorNetwork {
            graph: Self::generate_network_graph(tensors),
            params: Vec::new(),
        }
    }
}

impl<T> Default for TensorNetwork<T>
where
    T: TensorStructure,
{
    fn default() -> Self {
        Self::new()
    }
}

impl<T> TensorNetwork<T>
where
    T: TensorStructure,
{
    pub fn new() -> Self {
        TensorNetwork {
            graph: HalfEdgeGraph::new(),
            params: Vec::new(),
        }
    }

    pub fn scalar_mul<U>(&mut self, scalar: U)
    where
        for<'a> &'a T: ScalarMul<U, Output = T>,
    {
        if let Some((id, tensor)) = self.graph.nodes.iter().next() {
            self.graph.nodes[id] = tensor.scalar_mul(scalar).unwrap();
        }
    }

    pub fn push(&mut self, tensor: T) -> NodeId {
        let slots = tensor.external_structure().to_vec();
        self.graph.add_node_with_edges(tensor, &slots)
    }

    fn generate_network_graph(tensors: Vec<T>) -> HalfEdgeGraph<T, Slot> {
        let mut graph = HalfEdgeGraph::<T, Slot>::new();

        for tensor in tensors {
            let slots = tensor.external_structure().to_vec();
            graph.add_node_with_edges(tensor, &slots);
        }

        graph
    }

    pub fn edge_to_min_degree_node_with_depth(&self, depth: usize) -> Option<HedgeId>
    where
        T: TracksCount,
    {
        let mut neighs: SecondaryMap<NodeId, HedgeId> = self
            .graph
            .reverse_nodemap
            .clone()
            .into_iter()
            .filter(|(n, _e)| self.graph.nodes[*n].contractions_num() < depth)
            .collect();
        if neighs.is_empty() {
            return None;
        }

        loop {
            let mut all_ext = true;
            for (node, initial) in &mut neighs {
                *initial = self.graph.neighbors[*initial];
                let start = self.graph.reverse_nodemap[node];

                if self.graph.involution[start] != start
                    && self.graph.nodes[self.graph.nodemap[self.graph.involution[start]]]
                        .contractions_num()
                        < depth
                {
                    all_ext = false;
                    if *initial == start {
                        return Some(start);
                    }
                }
            }
            if all_ext {
                return None;
            }
        }
    }
}
impl<T> TensorNetwork<T>
where
    T: Clone,
{
    pub fn result(&self) -> T {
        self.graph.nodes.iter().next().unwrap().1.clone()
    }
}

impl<T> TensorNetwork<T> {
    pub fn dot(&self) -> std::string::String {
        self.graph.dot()
    }
}

impl<I> TensorNetwork<MixedTensor<I>>
where
    I: TensorStructure + Clone,
{
    pub fn generate_params(&mut self) {
        for (_i, n) in &self.graph.nodes {
            if n.is_symbolic() {
                self.params.push(n.clone());
            }
        }
    }
}

impl<T> TensorNetwork<T>
where
    T: TensorStructure<Structure = HistoryStructure<Symbol>> + Clone,
{
    pub fn symbolic_shadow(&mut self, name: &str) -> TensorNetwork<MixedTensors> {
        {
            for (i, n) in &mut self.graph.nodes {
                n.mut_structure().set_name(&State::get_symbol(format!(
                    "{}{}",
                    name,
                    i.data().as_ffi()
                )));
            }
        }

        let edges = self.graph.edges.clone();
        let involution = self.graph.involution.clone();
        let neighbors = self.graph.neighbors.clone();

        let mut nodes = DenseSlotMap::with_key();
        let mut nodemap = SecondaryMap::new();
        let mut reverse_nodemap = SecondaryMap::new();
        let mut params = Vec::new();

        for (i, n) in &self.graph.nodes {
            let node = n.structure().clone().shadow().unwrap();

            let nid = nodes.insert(MixedTensor::<HistoryStructure<Symbol>>::from(node.clone()));

            params.push(node.into());
            let mut first = true;
            for e in self.graph.edges_incident(i) {
                if first {
                    reverse_nodemap.insert(nid, e);
                    first = false;
                }
                nodemap.insert(e, nid);
            }
        }

        let g = HalfEdgeGraph {
            edges,
            involution,
            reverse_nodemap,
            neighbors,
            nodes,
            nodemap,
        };

        TensorNetwork { graph: g, params }
    }
}

impl<T> TensorNetwork<T>
where
    T: HasName,
{
    pub fn name(&mut self, name: T::Name)
    where
        T::Name: From<std::string::String> + Display,
    {
        for (id, n) in &mut self.graph.nodes {
            n.set_name(&format!("{}{}", name, id.data().as_ffi()).into());
        }
    }
}

impl<T> TensorNetwork<T>
where
    T: HasName<Name = Symbol>,
{
    pub fn namesym(&mut self, name: &str) {
        for (id, n) in &mut self.graph.nodes {
            n.set_name(&State::get_symbol(format!(
                "{}{}",
                name,
                id.data().as_ffi()
            )));
        }
    }
}

impl<T> TensorNetwork<T>
where
    T: Contract<T, LCM = T> + TensorStructure,
{
    pub fn contract_algo(&mut self, edge_choice: fn(&TensorNetwork<T>) -> Option<HedgeId>) {
        if let Some(e) = edge_choice(self) {
            self.contract_edge(e);
            // println!("{}", self.dot());
            self.contract_algo(edge_choice);
        }
    }
    fn contract_edge(&mut self, edge_idx: HedgeId) {
        let a = self.graph.nodemap[edge_idx];
        let b = self.graph.nodemap[self.graph.involution[edge_idx]];

        let ai = self.graph.nodes.get(a).unwrap();
        let bi = self.graph.nodes.get(b).unwrap();

        let f = ai.contract(bi).unwrap();

        self.graph.merge_nodes(a, b, f);
    }

    pub fn contract(&mut self) {
        self.contract_algo(Self::edge_to_min_degree_node)
    }
}
