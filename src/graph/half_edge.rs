use std::{collections::VecDeque, fmt::Display, hash::Hash, num::TryFromIntError};

use ahash::{AHashMap, AHashSet};
use bitvec::{slice::IterOnes, vec::BitVec};
use indexmap::IndexMap;
use itertools::Itertools;
use pathfinding::undirected::connected_components;
use rand::{rngs::SmallRng, Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use bitvec::prelude::*;
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct NodeIndex(pub usize);

impl std::fmt::Display for NodeIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

pub struct PowersetIterator {
    size: usize,
    current: usize,
}

impl PowersetIterator {
    pub fn new(n_elements: u8) -> Self {
        PowersetIterator {
            size: 1 << n_elements,
            current: 0,
        }
    }
}

impl Iterator for PowersetIterator {
    type Item = BitVec;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.size {
            let out = BitVec::<_, Lsb0>::from_element(self.current);

            self.current += 1;
            Some(out)
        } else {
            None
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum InvolutiveMapping<E> {
    Identity(Option<E>),
    Source((Option<E>, usize)),
    Sink(usize),
    // Undirected(usize),
    // UndirectedData((Option<E>, usize)),
}

impl<E> InvolutiveMapping<E> {
    pub fn is_identity(&self) -> bool {
        matches!(self, InvolutiveMapping::Identity(_))
    }

    pub fn map_data_ref<E2, F: FnMut(&E) -> E2>(&self, f: F) -> InvolutiveMapping<E2> {
        match self {
            InvolutiveMapping::Identity(d) => InvolutiveMapping::Identity(d.as_ref().map(f)),
            InvolutiveMapping::Source((d, s)) => InvolutiveMapping::Source((d.as_ref().map(f), *s)),
            InvolutiveMapping::Sink(s) => InvolutiveMapping::Sink(*s),
        }
    }

    pub fn map_data_none<E2>(&self) -> InvolutiveMapping<E2> {
        match self {
            InvolutiveMapping::Identity(_) => InvolutiveMapping::Identity(None),
            InvolutiveMapping::Source((_, s)) => InvolutiveMapping::Source((None, *s)),
            InvolutiveMapping::Sink(s) => InvolutiveMapping::Sink(*s),
        }
    }

    pub fn is_internal(&self) -> bool {
        !matches!(self, InvolutiveMapping::Identity(_))
    }

    pub fn make_source(&mut self, sink: usize) -> Option<InvolutiveMapping<E>> {
        let data = match self {
            InvolutiveMapping::Identity(d) => d.take(),
            _ => None,
        };
        Some(InvolutiveMapping::Source((data, sink)))
    }

    pub fn new_identity(data: E) -> Self {
        InvolutiveMapping::Identity(Some(data))
    }

    pub fn identity_dot(edge_id: usize, source: usize, attr: Option<&GVEdgeAttrs>) -> String {
        let mut out = "".to_string();
        out.push_str(&format!("ext{} [shape=none, label=\"\"];\n ", edge_id));
        if let Some(attr) = attr {
            out.push_str(&format!("{} -- ext{} {};\n ", source, edge_id, attr));
        } else {
            out.push_str(&format!("{} -- ext{} ;\n ", source, edge_id));
        }
        out
    }

    pub fn pair_dot(source: usize, sink: usize, attr: Option<&GVEdgeAttrs>) -> String {
        let mut out = "".to_string();
        if let Some(attr) = attr {
            out.push_str(&format!("  {} -- {} {};\n", source, sink, attr));
        } else {
            out.push_str(&format!(
                "  {} -- {} [color=\"red:blue;0.5 \" ];\n",
                source, sink
            ));
        }
        out
    }

    pub fn default_dot(
        &self,
        edge_id: usize,
        source: Option<usize>,
        sink: Option<usize>,
        attr: Option<&GVEdgeAttrs>,
    ) -> String {
        let mut out = "".to_string();
        match self {
            InvolutiveMapping::Identity(_) => {
                out.push_str(&Self::identity_dot(edge_id, source.unwrap(), attr));
            }
            InvolutiveMapping::Source((_, _)) => {
                out.push_str(&Self::pair_dot(source.unwrap(), sink.unwrap(), attr));
            }
            InvolutiveMapping::Sink(_) => {}
        }
        out
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Involution<N, E> {
    pub inv: Vec<(N, InvolutiveMapping<E>)>,
}

impl<N, E> Display for Involution<N, E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut out = "".to_string();
        for (i, (_, e)) in self.inv.iter().enumerate() {
            match e {
                InvolutiveMapping::Identity(_) => {
                    out.push_str(&format!("{}\n", i));
                }
                InvolutiveMapping::Source((_, s)) => {
                    out.push_str(&format!("{}->{}\n", i, s));
                }
                InvolutiveMapping::Sink(s) => {
                    out.push_str(&format!("{}<-{}\n", i, s));
                }
            }
        }
        write!(f, "{}", out)
    }
}

#[allow(dead_code)]
impl<N, E> Involution<N, E> {
    fn new() -> Self {
        Involution { inv: Vec::new() }
    }

    pub fn flip(&mut self, hedge: usize) {
        let pair = self.inv(hedge);

        self.inv.swap(hedge, pair);

        match &mut self.inv[hedge].1 {
            InvolutiveMapping::Identity(_) => {}
            InvolutiveMapping::Source((_, s)) => {
                *s = pair;
            }
            InvolutiveMapping::Sink(s) => {
                *s = pair;
            }
        }

        match &mut self.inv[pair].1 {
            InvolutiveMapping::Identity(_) => {}
            InvolutiveMapping::Source((_, s)) => {
                *s = hedge;
            }
            InvolutiveMapping::Sink(s) => {
                *s = hedge;
            }
        }
    }

    fn print<S: SubGraph>(
        &self,
        subgraph: &S,
        n_label: &impl Fn(&N) -> Option<String>,
        h_label: &impl Fn(&E) -> Option<String>,
    ) -> String {
        let mut out = "".to_string();
        for (i, (n, e)) in self.inv.iter().enumerate() {
            if !subgraph.is_included(i) {
                continue;
            }
            if let Some(l) = n_label(n) {
                out.push_str(&format!("{}:", l));
            }
            match e {
                InvolutiveMapping::Identity(_) => {
                    out.push_str(&format!("{}\n", i));
                }
                InvolutiveMapping::Source((d, s)) => {
                    if let Some(d) = d {
                        if let Some(l) = h_label(d) {
                            out.push_str(&format!("{}-{}->{}\n", i, s, l));
                        } else {
                            out.push_str(&format!("{}->{}\n", i, s));
                        }
                    } else {
                        out.push_str(&format!("{}->{}\n", i, s));
                    }
                }
                InvolutiveMapping::Sink(s) => {
                    out.push_str(&format!("{}<-{}\n", i, s));
                }
            }
        }
        out
    }

    fn map_data_ref<F, G, N2, E2>(&self, mut f: F, g: &G) -> Involution<N2, E2>
    where
        F: FnMut(&N) -> N2,
        G: FnMut(&E) -> E2 + Clone,
    {
        let inv = self
            .inv
            .iter()
            .map(|(n, e)| (f(n), e.map_data_ref(g.clone())))
            .collect_vec();

        Involution { inv }
    }

    fn forgetful_map_node_data_ref<N2, E2>(
        &self,
        mut f: impl FnMut(&N) -> N2,
    ) -> Involution<N2, E2> {
        let inv = self
            .inv
            .iter()
            .map(|(n, e)| (f(n), e.map_data_none()))
            .collect_vec();

        Involution { inv }
    }

    fn n_internals<S: SubGraph>(&self, subgraph: &S) -> usize {
        subgraph
            .included()
            .filter(|i| self.is_internal(*i, subgraph))
            .count()
            / 2
    }

    fn first_internal(&self) -> Option<usize> {
        self.inv
            .iter()
            .position(|(_, e)| !matches!(e, InvolutiveMapping::Identity(_)))
    }

    fn random(len: usize, seed: u64) -> Involution<Option<N>, ()> {
        let mut rng = SmallRng::seed_from_u64(seed);

        let mut inv = Involution::new();

        for _ in 0..len {
            let r = rng.gen_bool(0.1);
            if r {
                inv.add_identity((), None);
            } else {
                inv.add_pair((), None, None);
            }
        }

        inv
    }

    fn add_identity(&mut self, data: E, id: N) -> usize {
        let index = self.inv.len();
        self.inv.push((id, InvolutiveMapping::new_identity(data)));
        index
    }

    fn set_id(&mut self, index: usize, id: N) {
        self.inv[index].0 = id;
    }

    fn get_data(&self, index: usize) -> &E {
        match &self.inv[index].1 {
            InvolutiveMapping::Identity(data) => data.as_ref().unwrap(),
            InvolutiveMapping::Source((data, _)) => data.as_ref().unwrap(),
            InvolutiveMapping::Sink(i) => self.get_data(*i),
        }
    }

    fn set_data(&mut self, index: usize, data: E) -> Option<E> {
        if let InvolutiveMapping::Sink(i) = &self.inv[index].1 {
            return self.set_data(*i, data);
        }

        match &mut self.inv[index].1 {
            InvolutiveMapping::Identity(old_data) => {
                let old = old_data.take();
                *old_data = Some(data);
                old
            }
            InvolutiveMapping::Source((old_data, _)) => {
                let old = old_data.take();
                *old_data = Some(data);
                old
            }
            InvolutiveMapping::Sink(_) => unreachable!(),
        }
    }

    fn is_internal<S: SubGraph>(&self, index: usize, subgraph: &S) -> bool {
        if !subgraph.is_included(index) {
            return false;
        }
        match &self.inv[index].1 {
            InvolutiveMapping::Identity(_) => false,
            InvolutiveMapping::Source((_, i)) => subgraph.is_included(*i),
            InvolutiveMapping::Sink(i) => subgraph.is_included(*i),
        }
    }

    fn get_smart_data<S: SubGraph>(&self, index: usize, subgraph: &S) -> Option<&E> {
        //should be a better enum (none,internal,external)
        if subgraph.is_included(index) {
            match &self.inv[index].1 {
                InvolutiveMapping::Identity(data) => Some(data.as_ref().unwrap()),
                InvolutiveMapping::Source((data, _)) => Some(data.as_ref().unwrap()),
                InvolutiveMapping::Sink(i) => {
                    if subgraph.is_included(*i) {
                        None
                    } else {
                        Some(self.get_data(*i))
                    }
                }
            }
        } else {
            None
        }
    }

    pub fn find_from_data(&self, data: &E) -> Option<usize>
    where
        E: PartialEq,
    {
        (0..self.inv.len()).find(|&i| self.get_data(i) == data)
    }

    pub fn get_node_id(&self, index: usize) -> &N {
        &self.inv[index].0
    }

    pub fn get_connected_node_id(&self, index: usize) -> Option<&N> {
        match &self.inv[index].1 {
            InvolutiveMapping::Source((_, i)) => Some(&self.inv[*i].0),
            InvolutiveMapping::Sink(i) => Some(&self.inv[*i].0),
            InvolutiveMapping::Identity(_) => None,
        }
    }
    fn connect_to_identity(&mut self, source: usize, id: N) -> usize {
        let sink = self.inv.len();
        self.inv.push((id, InvolutiveMapping::Sink(source)));
        self.inv[source].1 = self.inv[source]
            .1
            .make_source(sink)
            .unwrap_or_else(|| panic!("Source must be an identity mapping"));
        sink
    }

    fn connect(&mut self, source: usize, sink: usize)
    where
        E: Clone,
    {
        self.inv[source].1 = self.inv[source].1.make_source(sink).unwrap();

        if let InvolutiveMapping::Identity(_) = &self.inv[sink].1 {
            self.inv[sink].1 = InvolutiveMapping::Sink(source);
        } else {
            panic!("Sink must be an identity mapping")
        }
    }

    fn add_pair(&mut self, data: E, source_id: N, sink_id: N) -> (usize, usize) {
        let source = self.add_identity(data, source_id);
        let sink = self.connect_to_identity(source, sink_id);
        (source, sink)
    }

    fn len(&self) -> usize {
        self.inv.len()
    }

    pub fn inv(&self, index: usize) -> usize {
        match &self.inv[index].1 {
            InvolutiveMapping::Identity(_) => index,
            InvolutiveMapping::Source((_, i)) => *i,
            InvolutiveMapping::Sink(i) => *i,
        }
    }

    pub fn is_identity(&self, index: usize) -> bool {
        matches!(&self.inv[index].1, InvolutiveMapping::Identity(_))
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Default)]
pub struct GVEdgeAttrs {
    pub label: Option<String>,
    pub color: Option<String>,
    pub other: Option<String>,
}

impl std::fmt::Display for GVEdgeAttrs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let out = format!(
            "[{}]",
            [
                ("label=", self.label.as_ref()),
                ("color=", self.color.as_ref()),
                ("", self.other.as_ref())
            ]
            .iter()
            .filter_map(|(prefix, x)| x.map(|s| format!("{}{}", prefix, s)))
            .join(",")
        );
        write!(f, "{}", out)
    }
}
pub mod subgraph;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HedgeGraph<E, V> {
    nodes: IndexMap<HedgeNode, V>, // Forest of nodes, that contain half-edges, with data E
    base_nodes: usize,
    pub involution: Involution<HedgeNode, E>, // Involution of half-edges
}

#[derive(Clone, Debug)]
pub struct HedgeNodeBuilder<V> {
    data: V,
    hedges: Vec<usize>,
}

#[derive(Clone, Debug)]
pub struct HedgeGraphBuilder<E, V> {
    nodes: Vec<HedgeNodeBuilder<V>>,
    involution: Involution<NodeIndex, E>,
}

pub struct NodeIterator<'a, E, V, I = IterOnes<'a, usize, Lsb0>> {
    graph: &'a HedgeGraph<E, V>,
    edges: I,
    seen: BitVec,
}

impl<'a, E, V, I: Iterator<Item = usize>> Iterator for NodeIterator<'a, E, V, I> {
    type Item = (&'a HedgeNode, &'a V);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(next) = self.edges.next() {
            let node = self.graph.get_incident_node_id(next);
            let node_pos = self.graph.get_node_pos(node);

            if self.seen[node_pos] {
                self.next()
            } else {
                self.seen.set(node_pos, true);
                Some((node, &self.graph.nodes[node_pos]))
            }
        } else {
            None
        }
    }
}

impl<E, V> HedgeGraphBuilder<E, V> {
    pub fn new() -> Self {
        HedgeGraphBuilder {
            nodes: Vec::new(),
            involution: Involution::new(),
        }
    }

    pub fn build(self) -> HedgeGraph<E, V> {
        self.into()
    }

    pub fn add_node(&mut self, data: V) -> NodeIndex {
        let index = self.nodes.len();
        self.nodes.push(HedgeNodeBuilder {
            data,
            hedges: Vec::new(),
        });
        NodeIndex(index)
    }

    pub fn add_edge(&mut self, source: NodeIndex, sink: NodeIndex, data: E) {
        let id = self.involution.add_pair(data, source, sink);
        self.nodes[source.0].hedges.push(id.0);

        self.nodes[sink.0].hedges.push(id.1);
    }

    pub fn add_external_edge(&mut self, source: NodeIndex, data: E) {
        let id = self.involution.add_identity(data, source);
        self.nodes[source.0].hedges.push(id);
    }
}

impl<E, V> Default for HedgeGraphBuilder<E, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<E, V> From<HedgeGraphBuilder<E, V>> for HedgeGraph<E, V> {
    fn from(builder: HedgeGraphBuilder<E, V>) -> Self {
        let len = builder.involution.len();
        let involution = Involution {
            inv: builder
                .involution
                .inv
                .into_iter()
                .map(|(n, i)| (HedgeNode::from_builder(&builder.nodes[n.0], len), i))
                .collect(),
        };
        let nodes: IndexMap<HedgeNode, V> = builder
            .nodes
            .into_iter()
            .map(|x| (HedgeNode::from_builder(&x, len), x.data))
            .collect();
        HedgeGraph {
            base_nodes: nodes.len(),
            nodes,
            involution,
        }
    }
}

impl<N: Clone, E: Clone> From<symbolica::graph::Graph<N, E>> for HedgeGraph<E, N> {
    fn from(graph: symbolica::graph::Graph<N, E>) -> Self {
        let mut builder = HedgeGraphBuilder::new();
        let mut map = AHashMap::new();

        for (i, node) in graph.nodes().iter().enumerate() {
            map.insert(i, builder.add_node(node.data.clone()));
        }

        for edge in graph.edges() {
            let vertices = edge.vertices;
            let source = map[&vertices.0];
            let sink = map[&vertices.1];
            builder.add_edge(source, sink, edge.data.clone());
        }

        builder.into()
    }
}

impl<'a> From<&'a BareGraph> for HedgeGraph<usize, usize> {
    fn from(value: &BareGraph) -> Self {
        let mut builder = HedgeGraphBuilder::new();
        let mut map = AHashMap::new();

        for (i, _) in value.vertices.iter().enumerate() {
            map.insert(i, builder.add_node(i));
        }

        for (i, edge) in value.edges.iter().enumerate() {
            let source = map[&edge.vertices[0]];
            let sink = map[&edge.vertices[1]];
            builder.add_edge(source, sink, i);
        }

        builder.into()
    }
}

use subgraph::{HedgeNode, InternalSubGraph, SubGraph};
use thiserror::Error;

use super::BareGraph;

#[derive(Error, Debug)]

pub enum HedgeError {
    #[error("Invalid start node")]
    InvalidStart,
}

impl<E, V> HedgeGraph<E, V> {
    // fn node_out_edges(&self, node_index: NodeIndex) -> Vec<usize> {
    //     let leaves = AHashSet::from_iter(
    //         self.nodes
    //             .get_leaves_from_node(node_index)
    //             .into_iter()
    //             .map(|x| *x),
    //     );

    //     let inv_leaves = leaves
    //         .clone()
    //         .into_iter()
    //         .map(|x| self.involutions.inv(x))
    //         .collect::<AHashSet<_>>();

    //     let externals: Vec<usize> = leaves.difference(&inv_leaves).map(|x| *x).collect();

    //     externals
    // }

    // fn out_degree(&self, node_index: NodeIndex) -> usize {
    //     self.node_out_edges(node_index).len()
    // }

    pub fn iter_egde_data<'a>(
        &'a self,
        subgraph: &'a InternalSubGraph,
    ) -> impl Iterator<Item = &E> + '_ {
        subgraph
            .filter
            .iter_ones()
            .map(|i| self.involution.get_data(i))
    }

    pub fn iter_internal_edge_data<'a>(
        &'a self,
        subgraph: &'a InternalSubGraph,
    ) -> impl Iterator<Item = &E> + '_ {
        subgraph
            .filter
            .iter_ones()
            .flat_map(|i| self.involution.get_smart_data(i, subgraph))
    }

    // pub fn strongly_disjoint(&self, graphs: impl IntoIterator<Item = SubGraph>) -> bool {
    //     let mut union: SubGraph = self.empty_filter().into();

    //     let mut nodes: SubGraph = self.empty_filter().into();

    //     for subgraph in graphs {
    //         if union.intersection(&subgraph).is_empty() {
    //             return false;
    //         }
    //         union.union_with(&subgraph);

    //         let nodes = subgraph.to_nesting_node(graph)
    //     }
    //     let internals = self.internal_graph.filter.clone() & &other.internal_graph.filter;

    //     let externals_in_self = self.internal_graph.filter.clone() & &other.;
    //     let externals_in_other = self.externalhedges.clone() & &other.internal_graph.filter;

    //     internals.count_ones() == 0
    //         && externals_in_self.count_ones() == 0
    //         && externals_in_other.count_ones() == 0
    // }

    pub fn is_connected(&self, subgraph: &InternalSubGraph) -> bool {
        let n_hedges = subgraph.filter.count_ones();
        let start = subgraph.filter.first_one();

        if let Some(start) = start {
            self.dfs_reach(subgraph, start).filter.count_ones() == n_hedges
        } else {
            true
        }
    }

    pub fn count_connected_components<S: SubGraph>(&self, subgraph: &S) -> usize {
        self.connected_components(subgraph).len()
    }

    pub fn connected_components<S: SubGraph>(&self, subgraph: &S) -> Vec<BitVec> {
        let mut visited_edges = self.empty_filter();
        let mut components = vec![];
        // Iterate over all edges in the subgraph
        for hedge_index in subgraph.included() {
            if !visited_edges[hedge_index] {
                // Perform DFS to find all reachable edges from this edge
                //
                let root_node = self.involution.get_node_id(hedge_index);
                let reachable_edges = TraversalTree::dfs(self, subgraph, root_node).covers();
                visited_edges.union_with(&reachable_edges);
                // let component = self.clean_subgraph(reachable_edges);
                components.push(reachable_edges);
            }
        }
        components
    }

    pub fn cyclotomatic_number(&self, subgraph: &InternalSubGraph) -> usize {
        let n_hedges = self.count_internal_edges(subgraph);
        // println!("n_hedges: {}", n_hedges);
        let n_nodes = self.number_of_nodes_in_subgraph(subgraph);
        // println!("n_nodes: {}", n_nodes);
        let n_components = self.count_connected_components(subgraph);
        // println!("n_components: {}", n_components);

        n_hedges - n_nodes + n_components
    }

    pub fn count_internal_edges(&self, subgraph: &InternalSubGraph) -> usize {
        let mut internal_edge_count = 0;

        // Iterate over all half-edges in the subgraph
        for hedge_index in subgraph.filter.iter_ones() {
            let inv_hedge_index = self.involution.inv(hedge_index);

            // Check if the involuted half-edge is also in the subgraph
            if subgraph.filter[inv_hedge_index] {
                // To avoid double-counting, only count when hedge_index < inv_hedge_index
                if hedge_index < inv_hedge_index {
                    internal_edge_count += 1;
                }
            }
        }

        internal_edge_count
    }

    pub fn dfs_reach(&self, subgraph: &InternalSubGraph, start: usize) -> InternalSubGraph {
        let mut dfs_reach = self.empty_filter();
        for i in pathfinding::directed::dfs::dfs_reach(start, |i| {
            self.succesors(subgraph, *i).into_iter()
        }) {
            dfs_reach.set(i, true);
        }

        InternalSubGraph::try_new(dfs_reach, self).unwrap()
    }

    pub fn succesors(&self, subgraph: &InternalSubGraph, start: usize) -> Vec<usize> {
        let mut succesors = Vec::new();
        if let Some(connected) = self.involution.get_connected_node_id(start) {
            for i in connected.hairs.iter_ones().filter(|x| subgraph.filter[*x]) {
                {
                    succesors.push(i);
                }
            }
        }
        succesors
    }

    fn combinations<const K: usize>(n: usize) -> Vec<[usize; K]> {
        let mut result = Vec::new();
        let mut current = [0; K]; // Fixed-size array to store combinations
        Self::generate_combinations::<K>(0, 0, n, &mut current, &mut result);
        result
    }

    fn generate_combinations<const K: usize>(
        depth: usize,
        start: usize,
        n: usize,
        current: &mut [usize; K],
        result: &mut Vec<[usize; K]>,
    ) {
        if depth == K {
            result.push(*current); // Push a copy of the current array
            return;
        }

        for i in start..n {
            current[depth] = i;
            Self::generate_combinations::<K>(depth + 1, i + 1, n, current, result);
        }
    }

    fn non_cut_edges_impl(
        &self,
        connected_components: usize,
        cyclotomatic_number: usize,
        start: usize,
        current: &mut BitVec,
        set: &mut AHashSet<BitVec>,
    ) {
        if current.count_ones() >= 2 * cyclotomatic_number {
            return;
        }
        if self.count_connected_components(&current.complement(self)) == connected_components {
            set.insert(current.clone());
        }

        for i in start..self.involution.len() {
            let j = self.involution.inv(i);

            if i != j {
                current.set(i, true);
                current.set(j, true);
            }
            self.non_cut_edges_impl(
                connected_components,
                cyclotomatic_number,
                i + 1,
                current,
                set,
            );

            current.set(i, false);
            current.set(j, false);
        }
    }

    pub fn non_cut_edges(&self) -> AHashSet<BitVec> {
        let connected_components = self.count_connected_components(&self.full_filter());

        let cyclotomatic_number = self.cyclotomatic_number(&self.full_node().internal_graph);

        let mut current = self.empty_filter();
        let mut set = AHashSet::new();

        self.non_cut_edges_impl(
            connected_components,
            cyclotomatic_number,
            0,
            &mut current,
            &mut set,
        );

        set
    }

    pub fn all_connectivity_spanning_subgraphs<S: SubGraph>(&self, subgraph: &S) -> AHashSet<S> {
        let n = self.count_connected_components(subgraph);
        let mut spanning_subgraphs = AHashSet::new();

        for n in self.iter_node_data(subgraph) {}
        spanning_subgraphs
    }

    pub fn neighbors<S: SubGraph>(&self, subgraph: &S, pos: usize) -> BitVec {
        subgraph.hairs(self.involution.get_node_id(pos))
    }

    pub fn connected_neighbors<S: SubGraph>(&self, subgraph: &S, pos: usize) -> Option<BitVec> {
        Some(subgraph.hairs(self.involution.get_connected_node_id(pos)?))
    }

    pub fn iter_egde_node<'a, S: SubGraph>(
        &'a self,
        subgraph: &'a S,
    ) -> impl Iterator<Item = &HedgeNode> + '_ {
        subgraph.included().map(|i| self.involution.get_node_id(i))
    }

    pub fn iter_node_data<'a, S: SubGraph>(
        &'a self,
        subgraph: &'a S,
    ) -> impl Iterator<Item = (&'a HedgeNode, &'a V)> {
        NodeIterator {
            graph: self,
            edges: subgraph.included(),
            seen: bitvec![usize, Lsb0; 0; self.base_nodes],
        }
    }

    pub fn get_node_pos(&self, node: &HedgeNode) -> usize {
        self.nodes.get_index_of(node).unwrap()
    }

    pub fn get_incident_node_id(&self, edge: usize) -> &HedgeNode {
        self.involution.get_node_id(edge)
    }

    pub fn n_hedges(&self) -> usize {
        self.involution.len()
    }

    pub fn n_nodes(&self) -> usize {
        self.base_nodes
    }

    pub fn n_externals(&self) -> usize {
        self.involution
            .inv
            .iter()
            .filter(|(_, e)| e.is_identity())
            .count()
    }

    pub fn n_internals(&self) -> usize {
        self.involution
            .inv
            .iter()
            .filter(|(_, e)| e.is_internal())
            .count()
            / 2
    }

    pub fn n_base_nodes(&self) -> usize {
        self.nodes.iter().filter(|(n, _)| n.is_node()).count()
    }

    pub fn random(nodes: usize, edges: usize, seed: u64) -> HedgeGraph<(), ()> {
        let mut inv: Involution<Option<HedgeNode>, ()> =
            Involution::<HedgeNode, ()>::random(edges, seed);

        let mut rng = SmallRng::seed_from_u64(seed);

        let mut externals = Vec::new();
        let mut sources = Vec::new();
        let mut sinks = Vec::new();

        for (i, e) in inv.inv.iter().enumerate() {
            let nodeid = HedgeNode::node_from_pos(&[i], inv.inv.len());
            match e.1 {
                InvolutiveMapping::Identity(_) => externals.push(nodeid),
                InvolutiveMapping::Source((_, _)) => sources.push(nodeid),
                InvolutiveMapping::Sink(_) => sinks.push(nodeid),
            }
        }

        while !externals.is_empty() {
            if rng.gen_bool(0.5) {
                let source_i = rng.gen_range(0..sources.len());

                sources[source_i].union_with(&externals.pop().unwrap());
            } else {
                let sink_i = rng.gen_range(0..sinks.len());

                sinks[sink_i].union_with(&externals.pop().unwrap());
            }
        }

        let mut lengthone = false;

        while sources.len() + sinks.len() > nodes {
            if rng.gen_bool(0.5) {
                if sources.len() <= 1 {
                    if lengthone {
                        break;
                    }
                    lengthone = true;
                    continue;
                }

                let idx1 = rng.gen_range(0..sources.len());
                let idx2 = rng.gen_range(0..sources.len() - 1);

                let n_i = sources.swap_remove(idx1);
                sources[idx2].union_with(&n_i);
            } else {
                if sinks.len() <= 1 {
                    if lengthone {
                        break;
                    }
                    lengthone = true;
                    continue;
                }
                let idx1 = rng.gen_range(0..sinks.len());
                let idx2 = rng.gen_range(0..sinks.len() - 1);
                let n_i = sinks.swap_remove(idx1);
                sinks[idx2].union_with(&n_i);
            }
        }

        for sink in &sinks {
            for i in sink.hairs.iter_ones() {
                inv.set_id(i, Some(sink.clone()));
            }
        }

        for source in &sources {
            for i in source.hairs.iter_ones() {
                inv.set_id(i, Some(source.clone()));
            }
        }

        let new_inv = Involution {
            inv: inv.inv.into_iter().map(|(n, e)| (n.unwrap(), e)).collect(),
        };

        let mut nodes = IndexMap::new();

        for n in sources {
            nodes.insert(n.clone(), ());
        }

        for n in sinks {
            nodes.insert(n.clone(), ());
        }

        HedgeGraph {
            base_nodes: nodes.len(),
            nodes,
            involution: new_inv,
        }
    }

    pub fn base_dot(&self) -> String {
        let mut out = "graph {\n ".to_string();
        out.push_str("  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\";\n ");
        for (i, (n, e)) in self.involution.inv.iter().enumerate() {
            out.push_str(
                &e.default_dot(
                    i,
                    self.nodes.get_index_of(n),
                    self.involution
                        .get_connected_node_id(i)
                        .map(|x| self.nodes.get_index_of(x).unwrap()),
                    None,
                ),
            );
        }
        out += "}";
        out
    }

    pub fn nesting_node_from_subgraph(&self, internal_graph: InternalSubGraph) -> HedgeNode {
        let mut hairs = bitvec![usize, Lsb0; 0; self.involution.len()];

        if !internal_graph.valid::<E, V>(self) {
            panic!("Invalid subgraph")
        }

        for i in internal_graph.included() {
            hairs |= &self.involution.inv[i].0.hairs;
        }

        HedgeNode {
            hairs: !(!hairs | &internal_graph.filter),
            internal_graph,
        }
    }

    fn nesting_node_fix(&self, node: &mut HedgeNode) {
        let mut externalhedges = bitvec![usize, Lsb0; 0; self.involution.len()];

        for i in node.internal_graph.filter.iter_ones() {
            externalhedges |= &self.involution.inv[i].0.hairs;
        }

        node.hairs = !(!externalhedges | &node.internal_graph.filter);
    }

    pub fn paired_filter_from_pos(&self, pos: &[usize]) -> BitVec {
        let mut filter = bitvec![usize, Lsb0; 0; self.involution.len()];

        for &i in pos {
            filter.set(i, true);
            filter.set(self.involution.inv(i), true);
        }

        filter
    }

    // pub fn filter_from_pos(&self, pos: &[usize]) -> BitVec {
    //     Nested<T>::filter_from_pos(pos, self.involution.len())
    // }

    // pub fn nesting_node_from_pos(&self, pos: &[usize]) -> Nested<T> {
    //     self.nesting_node_from_subgraph(SubGraph::from(self.filter_from_pos(pos)))
    // }

    fn remove_externals(&self, subgraph: &mut HedgeNode) {
        let externals = self.external_filter();

        subgraph.internal_graph.filter &= !externals;
    }

    fn external_filter(&self) -> BitVec {
        let mut filter = bitvec![usize, Lsb0; 0; self.involution.len()];

        for (i, (_, edge)) in self.involution.inv.iter().enumerate() {
            if let InvolutiveMapping::Identity(_) = edge {
                filter.set(i, true);
            }
        }

        filter
    }

    pub fn full_filter(&self) -> BitVec {
        bitvec![usize, Lsb0; 1; self.involution.len()]
    }

    pub fn empty_filter(&self) -> BitVec {
        bitvec![usize, Lsb0; 0; self.involution.len()]
    }

    pub fn clean_subgraph(&self, filter: BitVec) -> InternalSubGraph {
        InternalSubGraph::cleaned_filter_optimist(filter, self)
    }

    pub fn full_node(&self) -> HedgeNode {
        self.nesting_node_from_subgraph(self.full_graph())
    }

    pub fn cycle_basis(&self) -> (Vec<InternalSubGraph>, TraversalTree) {
        self.paton_cycle_basis(&self.full_graph(), self.involution.get_node_id(0))
            .unwrap()
    }

    ///Read, R.C. and Tarjan, R.E. (1975), Bounds on Backtrack Algorithms for Listing Cycles, Paths, and Spanning Trees. Networks, 5: 237-252. https://doi.org/10.1002/net.1975.5.3.237
    // pub fn read_tarjan(&self) -> Vec<HedgeNode> {
    //     todo!("Implement")
    // }

    pub fn backtrack(
        &self,
        options: &mut InternalSubGraph,
        tree: &mut InternalSubGraph,
        current: usize,
    ) -> Option<usize> {
        tree.filter.set(self.involution.inv(current), false);
        tree.filter.set(current, false);

        let current_node = &self.involution.get_node_id(current).hairs;

        let current = current_node.iter_ones().find(|&i| options.filter[i]);

        if let Some(current) = current {
            tree.filter.set(current, true);
            options.filter.set(current, false);
            Some(current)
        } else {
            let current = current_node.iter_ones().find_map(|i| {
                if tree.filter[i] {
                    Some(self.involution.inv(i))
                } else {
                    None
                }
            });
            self.backtrack(options, tree, current?)
        }
    }

    pub fn order_basis(&self, basis: &[HedgeNode]) -> Vec<Vec<InternalSubGraph>> {
        let mut seen = vec![basis[0].internal_graph.clone()];
        let mut partitions = vec![seen.clone()];

        for cycle in basis.iter() {
            if seen
                .iter()
                .any(|p| !p.empty_intersection(&cycle.internal_graph))
            {
                partitions
                    .last_mut()
                    .unwrap()
                    .push(cycle.internal_graph.clone());
            } else {
                for p in partitions.last().unwrap() {
                    seen.push(p.clone());
                }
                partitions.push(vec![cycle.internal_graph.clone()]);
            }
        }

        partitions
    }

    pub fn all_cycles(&self) -> Vec<InternalSubGraph> {
        InternalSubGraph::all_sym_diff_powerset(&self.cycle_basis().0, &|mut c| {
            if self.is_cycle(&c) {
                c.loopcount = Some(1);
                Some(c)
            } else {
                None
            }
        })
        .unwrap()
        .into_iter()
        .collect()
    }

    pub fn all_cycle_sym_diffs(&self) -> Result<Vec<InternalSubGraph>, TryFromIntError> {
        InternalSubGraph::all_sym_diff_powerset(&self.cycle_basis().0, &Some)
            .map(|a| a.into_iter().collect())
    }

    pub fn all_cycle_unions(&self) -> AHashSet<InternalSubGraph> {
        InternalSubGraph::all_unions_iterative(&self.all_cycle_sym_diffs().unwrap())
    }

    fn is_cycle(&self, subgraph: &InternalSubGraph) -> bool {
        let mut is = true;
        for e in self.iter_egde_node(subgraph) {
            let adgacent = subgraph.filter.clone() & &e.hairs;
            if adgacent.count_ones() != 2 {
                is = false;
                break;
            }
        }
        is
    }

    pub fn full_graph(&self) -> InternalSubGraph {
        InternalSubGraph::cleaned_filter_optimist(self.full_filter(), self)
    }

    pub fn empty_subgraph<S: SubGraph>(&self) -> S {
        S::empty(self.n_hedges())
    }

    pub fn paton_count_loops(
        &self,
        subgraph: &InternalSubGraph,
        start: &HedgeNode,
    ) -> Result<usize, HedgeError> {
        let tree = TraversalTree::dfs(self, subgraph, start);

        let cuts = subgraph.subtract(&tree.tree);
        Ok(self.involution.n_internals(&cuts))
    }

    pub fn number_of_nodes_in_subgraph(&self, subgraph: &InternalSubGraph) -> usize {
        self.iter_node_data(subgraph).count()
    }

    pub fn node_degrees_in_subgraph(&self, subgraph: &InternalSubGraph) -> AHashMap<usize, usize> {
        let mut degrees = AHashMap::new();

        for (node, _) in self.iter_node_data(subgraph) {
            let node_pos = self.get_node_pos(node);

            // Count the number of edges in the subgraph incident to this node
            let incident_edges = node.hairs.clone() & &subgraph.filter;
            let degree = incident_edges.count_ones();

            degrees.insert(node_pos, degree);
        }

        degrees
    }

    pub fn hairy_from_filter(&self, filter: BitVec) -> HedgeNode {
        self.nesting_node_from_subgraph(InternalSubGraph::cleaned_filter_pessimist(filter, self))
    }

    fn paton_cycle_basis(
        &self,
        subgraph: &InternalSubGraph,
        start: &HedgeNode,
    ) -> Result<(Vec<InternalSubGraph>, TraversalTree), HedgeError> {
        if subgraph.filter.empty_intersection(&start.hairs) {
            return Err(HedgeError::InvalidStart);
        }

        let tree = TraversalTree::dfs(self, subgraph, start);

        let cuts = subgraph.subtract(&tree.tree);

        let mut cycle_basis = Vec::new();

        for c in cuts.included() {
            if c > self.involution.inv(c) {
                cycle_basis.push(InternalSubGraph::try_new(tree.cycle(c).unwrap(), self).unwrap());

                cycle_basis.last_mut().unwrap().loopcount = Some(1);
            }
        }

        Ok((cycle_basis, tree))
    }

    pub fn dot_impl<S: SubGraph>(
        &self,
        node_as_graph: &S,
        graph_info: String,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> String {
        node_as_graph.dot(self, graph_info, edge_attr, node_attr)
    }
    pub fn dot<S: SubGraph>(&self, node_as_graph: &S) -> String {
        self.dot_impl(node_as_graph, "".to_string(), &|_| None, &|_| None)
    }

    pub fn cut_branches(&self, subgraph: &mut HedgeNode) {
        let nodes = AHashSet::<&HedgeNode>::from_iter(
            subgraph
                .internal_graph
                .filter
                .iter_ones()
                .map(|i| self.involution.get_node_id(i)),
        );
        self.remove_externals(subgraph);

        let mut has_branch = true;
        while has_branch {
            has_branch = false;

            for n in &nodes {
                let int = n.hairs.clone() & &subgraph.internal_graph.filter;

                if int.count_ones() == 1 {
                    subgraph
                        .internal_graph
                        .filter
                        .set(int.first_one().unwrap(), false);
                    subgraph
                        .internal_graph
                        .filter
                        .set(self.involution.inv(int.first_one().unwrap()), false);
                    has_branch = true;
                }
            }
        }

        self.nesting_node_fix(subgraph);
    }

    pub fn get_edge_data(&self, edge: usize) -> &E {
        self.involution.get_data(edge)
    }

    pub fn all_spinneys_with_basis(&self, basis: &[&InternalSubGraph]) -> AHashSet<HedgeNode> {
        let mut spinneys = AHashSet::new();
        let mut base_cycle: InternalSubGraph = self.empty_subgraph();

        for cycle in basis {
            base_cycle.sym_diff_with(cycle);
        }

        spinneys.insert(self.nesting_node_from_subgraph(base_cycle.clone()));

        if basis.len() == 1 {
            return spinneys;
        }

        for i in 0..basis.len() {
            for s in self.all_spinneys_with_basis(
                &basis
                    .iter()
                    .enumerate()
                    .filter_map(|(j, s)| if j != i { Some(*s) } else { None })
                    .collect_vec(),
            ) {
                spinneys
                    .insert(self.nesting_node_from_subgraph(s.internal_graph.union(&base_cycle)));
                spinneys.insert(s);
            }
        }

        spinneys
    }

    pub fn all_spinneys_rec(&self, spinneys: &mut AHashSet<HedgeNode>, cycle_sums: Vec<HedgeNode>) {
        let _len = spinneys.len();

        let mut pset = PowersetIterator::new(cycle_sums.len() as u8);

        pset.next(); //Skip empty set

        for (ci, cj) in cycle_sums.iter().tuple_combinations() {
            let _union = ci.internal_graph.union(&cj.internal_graph);

            // spinneys.insert(union);
        }
    }

    pub fn all_spinneys(
        &self,
    ) -> AHashMap<InternalSubGraph, Vec<(InternalSubGraph, Option<InternalSubGraph>)>> {
        let mut cycles = self.cycle_basis().0;

        let mut all_combinations = PowersetIterator::new(cycles.len() as u8);
        all_combinations.next(); //Skip empty set

        let mut spinneys: AHashMap<
            InternalSubGraph,
            Vec<(InternalSubGraph, Option<InternalSubGraph>)>,
        > = AHashMap::new();

        for p in all_combinations {
            let mut base_cycle: InternalSubGraph = self.empty_subgraph();

            for i in p.iter_ones() {
                base_cycle.sym_diff_with(&cycles[i]);
            }

            cycles.push(base_cycle);
        }

        for (ci, cj) in cycles.iter().tuple_combinations() {
            let union = ci.union(cj);

            if let Some(v) = spinneys.get_mut(&union) {
                v.push((ci.clone(), Some(cj.clone())));
            } else {
                spinneys.insert(union, vec![(ci.clone(), Some(cj.clone()))]);
            }
        }

        for c in cycles {
            spinneys.insert(c.clone(), vec![(c.clone(), None)]);
        }
        spinneys
    }

    pub fn all_spinneys_alt(&self) -> AHashSet<InternalSubGraph> {
        let mut spinneys = AHashSet::new();
        let cycles = self.all_cycles();

        let mut pset = PowersetIterator::new(cycles.len() as u8);
        pset.next(); //Skip empty set

        for p in pset {
            let mut union: InternalSubGraph = self.empty_subgraph();

            for i in p.iter_ones() {
                union.union_with(&cycles[i]);
            }

            spinneys.insert(union);
        }

        for c in cycles {
            spinneys.insert(c);
        }
        spinneys
    }

    pub fn all_s_t_cuts(
        &self,
        s: &HedgeNode,
        t: &HedgeNode,
        regions: &mut AHashSet<InternalSubGraph>,
    ) {
        let mut new_internals = vec![];
        for h in s.hairs.iter_ones() {
            let invh = self.involution.inv(h);

            if h > invh && s.hairs[self.involution.inv(h)] {
                new_internals.push(h);
            }
        }

        let mut new_node = s.clone();

        for h in new_internals {
            new_node.hairs.set(h, false);
            new_node.hairs.set(self.involution.inv(h), false);
            new_node.internal_graph.filter.set(h, true);
            new_node
                .internal_graph
                .filter
                .set(self.involution.inv(h), true);
        }

        let complement = new_node.complement(self).hairs;

        let count = self.count_connected_components(&complement);

        if count == 1 && !regions.insert(new_node.internal_graph) {
            return;
        }

        for h in new_node.hairs.iter_ones() {
            let invh = self.involution.inv(h);

            if invh != h && !t.hairs[invh] {
                let mut new_node = s.clone();
                new_node
                    .hairs
                    .union_with(&self.get_incident_node_id(invh).hairs);

                new_node.hairs.set(h, false);
                new_node.hairs.set(invh, false);
                new_node.internal_graph.filter.set(h, true);
                new_node.internal_graph.filter.set(invh, true);
                self.all_s_t_cuts(&new_node, t, regions);
            }
        }
    }
}

pub struct TraversalTree {
    pub traversal: Vec<usize>,
    inv: Involution<Parent, ()>,
    pub tree: InternalSubGraph,
}

pub enum Parent {
    Unset,
    Root,
    Hedge(usize, usize),
}

impl TraversalTree {
    fn covers(&self) -> BitVec {
        let mut covers = <BitVec as SubGraph>::empty(self.inv.inv.len());
        for (i, (m, _)) in self.inv.inv.iter().enumerate() {
            match m {
                Parent::Unset => {}
                _ => {
                    covers.set(i, true);
                }
            }
        }
        covers
    }

    fn path_to_root(&self, start: usize) -> BitVec {
        let mut path = <BitVec as SubGraph>::empty(self.inv.inv.len());
        let mut current = start;
        path.set(current, true);

        while let Parent::Hedge(p, _) = self.inv.get_node_id(current) {
            path.set(*p, true);
            current = self.inv.inv(*p);
            path.set(current, true);
        }
        path
    }

    pub fn cycle(&self, cut: usize) -> Option<BitVec> {
        match self.inv.get_node_id(cut) {
            Parent::Hedge(p, _) => {
                if *p == cut {
                    return None;
                }
            }
            Parent::Root => {}
            _ => return None,
        }

        let cut_pair = self.inv.inv(cut);
        match self.inv.get_node_id(cut_pair) {
            Parent::Hedge(p, _) => {
                if *p == cut_pair {
                    return None;
                }
            }
            Parent::Root => {}
            _ => return None,
        }

        let mut cycle = self.path_to_root(cut);
        cycle.sym_diff_with(&self.path_to_root(cut_pair));
        Some(cycle)
    }

    pub fn bfs<E, V, S: SubGraph>(
        graph: &HedgeGraph<E, V>,
        subgraph: &S,
        root_node: &HedgeNode,
        // target: Option<&HedgeNode>,
    ) -> Self {
        let mut queue = VecDeque::new();
        let mut seen = subgraph.hairs(root_node);

        let mut traversal: Vec<usize> = Vec::new();
        let mut involution: Involution<Parent, ()> = graph
            .involution
            .forgetful_map_node_data_ref(|_| Parent::Unset);

        // add all hedges from root node that are not self loops
        // to the queue
        // They are all potential branches
        for i in seen.iter_ones() {
            involution.set_id(i, Parent::Root);
            if !seen[graph.involution.inv(i)] {
                // if not self loop
                queue.push_back(i)
            }
        }
        while let Some(hedge) = queue.pop_front() {
            // if the hedge is not external get the neighbors of the paired hedge
            if let Some(cn) = graph.connected_neighbors(subgraph, hedge) {
                let connected = involution.inv(hedge);

                if !seen[connected] {
                    // if this new hedge hasn't been seen before, it means the node it belongs to
                    // is a new node in the traversal
                    traversal.push(connected);
                } else {
                    continue;
                }
                // mark the new node as seen
                seen.union_with(&cn);

                // for all hedges in this new node, they have a parent, the initial hedge
                for i in cn.iter_ones() {
                    if let Parent::Unset = involution.inv[i].0 {
                        involution.set_id(i, Parent::Hedge(connected, traversal.len()));
                    }
                    // if they lead to a new node, they are potential branches, add them to the queue
                    if !seen[involution.inv(i)] {
                        queue.push_back(i);
                    }
                }
            }
        }

        TraversalTree::new(graph, traversal, involution)
    }

    pub fn new<E, V>(
        graph: &HedgeGraph<E, V>,
        traversal: Vec<usize>,
        inv: Involution<Parent, ()>,
    ) -> Self {
        let mut tree = graph.empty_filter();

        for (i, j) in traversal.iter().map(|x| (*x, inv.inv(*x))) {
            tree.set(i, true);
            tree.set(j, true);
        }

        TraversalTree {
            traversal,
            inv,
            tree: InternalSubGraph::cleaned_filter_optimist(tree, graph),
        }
    }

    pub fn dfs<E, V, S: SubGraph>(
        graph: &HedgeGraph<E, V>,
        subgraph: &S,
        root_node: &HedgeNode,
        // target: Option<&HedgeNode>,
    ) -> Self {
        let mut stack = Vec::new();
        let mut seen = subgraph.hairs(root_node);

        let mut traversal: Vec<usize> = Vec::new();
        let mut involution: Involution<Parent, ()> = graph
            .involution
            .forgetful_map_node_data_ref(|_| Parent::Unset);

        // add all hedges from root node that are not self loops
        // to the stack
        // They are all potential branches
        for i in seen.iter_ones() {
            involution.set_id(i, Parent::Root);
            if !seen[graph.involution.inv(i)] {
                // if not self loop
                stack.push(i)
            }
        }
        while let Some(hedge) = stack.pop() {
            // if the hedge is not external get the neighbors of the paired hedge
            if let Some(cn) = graph.connected_neighbors(subgraph, hedge) {
                let connected = involution.inv(hedge);

                if !seen[connected] {
                    // if this new hedge hasn't been seen before, it means the node it belongs to
                    // is a new node in the traversal
                    traversal.push(connected);
                } else {
                    continue;
                }

                // mark the new node as seen
                seen.union_with(&cn);

                for i in cn.iter_ones() {
                    if let Parent::Unset = involution.inv[i].0 {
                        involution.set_id(i, Parent::Hedge(connected, traversal.len()));
                    }

                    if !seen[involution.inv(i)] {
                        stack.push(i);
                    }
                }
            }
        }

        TraversalTree::new(graph, traversal, involution)
    }
}

#[cfg(test)]
mod tests;
