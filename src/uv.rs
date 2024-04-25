use std::{clone, cmp::Ordering, collections::VecDeque, hash::Hash, iter::Map, ops::Sub};

use ahash::AHashSet;
use bitvec::{
    slice::IterOnes,
    vec::{self, BitVec},
};
use indexmap::IndexMap;
use itertools::Itertools;
use petgraph::{algo::Cycle, graph};
use rand::{Rng, SeedableRng};
use rand_xoshiro::Xoroshiro64Star;
use serde::{Deserialize, Serialize};

use crate::{
    cross_section::SuperGraph,
    graph::{Edge, EdgeType, Graph, LoopMomentumBasis, Vertex},
};

use bitvec::prelude::*;
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct NodeIndex(usize);

impl std::fmt::Display for NodeIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

struct PowersetIterator {
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
enum InvolutiveMapping<E> {
    Identity(Option<E>),
    Source((Option<E>, usize)),
    Sink(usize),
}

impl<E> InvolutiveMapping<E> {
    pub fn is_identity(&self) -> bool {
        matches!(self, InvolutiveMapping::Identity(_))
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
struct Involution<N, E> {
    inv: Vec<(N, InvolutiveMapping<E>)>,
}

#[allow(dead_code)]
impl<N, E> Involution<N, E> {
    fn new() -> Self {
        Involution { inv: Vec::new() }
    }

    fn first_internal(&self) -> Option<usize> {
        self.inv
            .iter()
            .position(|(_, e)| !matches!(e, InvolutiveMapping::Identity(_)))
    }

    fn random(len: usize, seed: u64) -> Involution<Option<N>, ()> {
        let mut rng = Xoroshiro64Star::seed_from_u64(seed);

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

#[derive(Clone, Debug, Serialize, Deserialize)]
struct GVEdgeAttrs {
    label: Option<String>,
    color: Option<String>,
}

impl Default for GVEdgeAttrs {
    fn default() -> Self {
        GVEdgeAttrs {
            label: None,
            color: None,
        }
    }
}

impl std::fmt::Display for GVEdgeAttrs {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut out = "[".to_string();
        if let Some(label) = &self.label {
            out.push_str(&format!("label=\"{}\",", label));
        }
        if let Some(color) = &self.color {
            out.push_str(&format!("color=\"{}\",", color));
        }
        out.push(']');
        write!(f, "{}", out)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct NestingGraph<E, V> {
    nodes: IndexMap<NestingNode, V>, // Forest of nodes, that contain half-edges, with data E
    base_nodes: usize,
    involution: Involution<NestingNode, E>, // Involution of half-edges
}

#[derive(Clone, Debug, Serialize, Deserialize, Eq)]
pub struct SubGraph {
    filter: BitVec,
    loopcount: Option<usize>,
}

impl Hash for SubGraph {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.filter.hash(state);
    }
}

impl PartialEq for SubGraph {
    fn eq(&self, other: &Self) -> bool {
        self.filter == other.filter
    }
}

impl From<BitVec> for SubGraph {
    fn from(filter: BitVec) -> Self {
        SubGraph {
            filter,
            loopcount: None,
        }
    }
}

impl PartialOrd for SubGraph {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self == other {
            Some(std::cmp::Ordering::Equal)
        } else if self.filter.clone() | &other.filter == self.filter {
            Some(std::cmp::Ordering::Greater)
        } else if self.filter.clone() | &other.filter == other.filter {
            Some(std::cmp::Ordering::Less)
        } else {
            None
        }
    }
}

impl SubGraph {
    pub fn is_empty(&self) -> bool {
        self.filter.count_ones() == 0
    }

    pub fn set_loopcount<E, V>(&mut self, graph: &NestingGraph<E, V>) {
        self.loopcount = Some(
            graph
                .paton_count_loops(self, self.filter.first_one().unwrap())
                .unwrap(),
        );
    }

    pub fn intersect_with(&mut self, other: &SubGraph) {
        self.filter &= &other.filter;
    }

    pub fn union_with(&mut self, other: &SubGraph) {
        self.filter |= &other.filter;
    }

    pub fn sym_diff_with(&mut self, other: &SubGraph) {
        self.filter ^= &other.filter;
    }

    pub fn intersection(&self, other: &SubGraph) -> SubGraph {
        let mut new = self.clone();
        new.intersect_with(other);
        new
    }

    pub fn union(&self, other: &SubGraph) -> SubGraph {
        let mut new = self.clone();
        new.union_with(other);
        new
    }

    pub fn sym_diff(&self, other: &SubGraph) -> SubGraph {
        let mut new = self.clone();
        new.sym_diff_with(other);
        new
    }

    pub fn empty_intersection(&self, other: &SubGraph) -> bool {
        (self.filter.clone() & &other.filter).count_ones() == 0
    }

    pub fn empty_union(&self, other: &SubGraph) -> bool {
        (self.filter.clone() | &other.filter).count_ones() == 0
    }

    pub fn complement(&self) -> SubGraph {
        SubGraph {
            filter: !self.filter.clone(),
            loopcount: None,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct NestingNode {
    internal_graph: SubGraph,
    externalhedges: BitVec, // essentially all hedges that are not in the internal graph, but are connected to it
}

impl PartialOrd for NestingNode {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        if self == other {
            return Some(std::cmp::Ordering::Equal);
        }

        if self.internal_graph.filter.len() != other.internal_graph.filter.len() {
            return None;
        }

        if self.externalhedges.len() != other.externalhedges.len() {
            return None;
        }

        let self_all = self.all_edges();
        let other_all = other.all_edges();

        let all_order = if self_all == other_all {
            Some(std::cmp::Ordering::Equal)
        } else if self_all.clone() | &other_all == self_all {
            Some(std::cmp::Ordering::Greater)
        } else if self_all.clone() | &other_all == other_all {
            Some(std::cmp::Ordering::Less)
        } else {
            None
        };

        let internal_order = self.internal_graph.partial_cmp(&other.internal_graph);

        match (all_order, internal_order) {
            (Some(std::cmp::Ordering::Equal), Some(std::cmp::Ordering::Equal)) => {
                Some(std::cmp::Ordering::Equal)
            }
            (Some(all), Some(internal)) => {
                if all == internal {
                    Some(all)
                } else {
                    None
                }
            }
            _ => None,
        }
    }
}

impl NestingNode {
    pub fn new(len: usize) -> Self {
        NestingNode {
            internal_graph: bitvec![usize, Lsb0; 0; len].into(),
            externalhedges: bitvec![usize, Lsb0; 0; len],
        }
    }

    fn all_edges(&self) -> BitVec {
        self.internal_graph.filter.clone() | &self.externalhedges
    }

    pub fn union_with(&mut self, other: &NestingNode) {
        self.internal_graph.filter |= &other.internal_graph.filter;
        self.externalhedges |= &other.externalhedges;
        self.externalhedges &= !self.internal_graph.filter.clone();
    }

    pub fn intersect_with(&mut self, other: &NestingNode) {
        self.internal_graph.filter &= &other.internal_graph.filter;
        self.externalhedges &= &other.externalhedges;
    }

    #[must_use]
    pub fn intersection(&self, other: &NestingNode) -> NestingNode {
        let mut new = self.clone();
        new.intersect_with(other);
        new
    }

    #[must_use]
    pub fn union(&self, other: &NestingNode) -> NestingNode {
        let mut new = self.clone();
        new.union_with(other);
        new
    }

    pub fn empty_intersection(&self, other: &NestingNode) -> bool {
        let internals = self.internal_graph.filter.clone() & &other.internal_graph.filter;

        let externals = self.externalhedges.clone() & &other.externalhedges;

        internals.count_ones() == 0 && externals.count_ones() == 0
    }

    pub fn weakly_disjoint(&self, other: &NestingNode) -> bool {
        let internals = self.internal_graph.filter.clone() & &other.internal_graph.filter;

        internals.count_ones() == 0
    }

    pub fn strongly_disjoint(&self, other: &NestingNode) -> bool {
        let internals = self.internal_graph.filter.clone() & &other.internal_graph.filter;

        let externals = self.externalhedges.clone() & &other.externalhedges;

        internals.count_ones() == 0 && externals.count_ones() == 0
    }

    pub fn node_from_pos(pos: &[usize], len: usize) -> NestingNode {
        NestingNode {
            externalhedges: NestingNode::filter_from_pos(pos, len),
            internal_graph: bitvec![usize, Lsb0; 0; len].into(),
        }
    }

    pub fn filter_from_pos(pos: &[usize], len: usize) -> BitVec {
        let mut filter = bitvec![usize, Lsb0; 0; len];

        for &i in pos {
            filter.set(i, true);
        }

        filter
    }

    pub fn internal_graph_union(&self, other: &NestingNode) -> SubGraph {
        SubGraph {
            filter: self.internal_graph.filter.clone() | &other.internal_graph.filter,
            loopcount: None,
        }
    }

    pub fn from_builder<V>(builder: &NestingNodeBuilder<V>, len: usize) -> Self {
        let internal_graph = bitvec![usize, Lsb0; 0; len].into();
        let mut externalhedges = bitvec![usize, Lsb0; 0; len];

        for hedge in &builder.hedges {
            let mut bit = externalhedges.get_mut(*hedge).unwrap();
            *bit = true;
        }

        NestingNode {
            internal_graph,
            externalhedges,
        }
    }

    pub fn is_node(&self) -> bool {
        self.internal_graph.is_empty()
    }

    pub fn is_subgraph(&self) -> bool {
        !self.is_node()
    }
}

#[derive(Clone, Debug)]
pub struct NestingNodeBuilder<V> {
    data: V,
    hedges: Vec<usize>,
}

pub struct NestingGraphBuilder<E, V> {
    nodes: Vec<NestingNodeBuilder<V>>,
    involution: Involution<NodeIndex, E>,
}

pub struct NodeIterator<'a, E, V> {
    graph: &'a NestingGraph<E, V>,
    edges: IterOnes<'a, usize, Lsb0>,
    seen: BitVec,
}

impl<'a, E, V> Iterator for NodeIterator<'a, E, V> {
    type Item = (&'a NestingNode, &'a V);

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

impl<E, V> NestingGraphBuilder<E, V> {
    pub fn new() -> Self {
        NestingGraphBuilder {
            nodes: Vec::new(),
            involution: Involution::new(),
        }
    }

    pub fn build(self) -> NestingGraph<E, V> {
        self.into()
    }

    pub fn add_node(&mut self, data: V) -> NodeIndex {
        let index = self.nodes.len();
        self.nodes.push(NestingNodeBuilder {
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

impl<E, V> Default for NestingGraphBuilder<E, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<E, V> From<NestingGraphBuilder<E, V>> for NestingGraph<E, V> {
    fn from(builder: NestingGraphBuilder<E, V>) -> Self {
        let len = builder.involution.len();
        let involution = Involution {
            inv: builder
                .involution
                .inv
                .into_iter()
                .map(|(n, i)| (NestingNode::from_builder(&builder.nodes[n.0], len), i))
                .collect(),
        };
        let nodes: IndexMap<NestingNode, V> = builder
            .nodes
            .into_iter()
            .map(|x| (NestingNode::from_builder(&x, len), x.data))
            .collect();
        NestingGraph {
            base_nodes: nodes.len(),
            nodes,
            involution,
        }
    }
}

use thiserror::Error;

#[derive(Error, Debug)]

enum NestingNodeError {
    #[error("Invalid start node")]
    InvalidStart,
}

impl<E, V> NestingGraph<E, V> {
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

    pub fn iter_egde_data<'a>(&'a self, subgraph: &'a SubGraph) -> impl Iterator<Item = &E> + '_ {
        subgraph
            .filter
            .iter_ones()
            .map(|i| self.involution.get_data(i))
    }

    pub fn is_connected(&self, subgraph: &SubGraph) -> bool {
        let n_hedges = subgraph.filter.count_ones();
        let start = subgraph.filter.first_one();

        if let Some(start) = start {
            self.dfs_reach(subgraph, start).len() == n_hedges
        } else {
            true
        }
    }

    pub fn dfs_reach(&self, subgraph: &SubGraph, start: usize) -> Vec<usize> {
        pathfinding::directed::dfs::dfs_reach(start, |i| self.succesors(subgraph, *i).into_iter())
            .collect()
    }

    pub fn succesors(&self, subgraph: &SubGraph, start: usize) -> Vec<usize> {
        let mut succesors = Vec::new();
        if let Some(connected) = self.involution.get_connected_node_id(start) {
            for i in connected
                .externalhedges
                .iter_ones()
                .filter(|x| subgraph.filter[*x])
            {
                {
                    succesors.push(i);
                }
            }
        }
        succesors
    }

    pub fn iter_egde_node<'a>(
        &'a self,
        subgraph: &'a SubGraph,
    ) -> impl Iterator<Item = &NestingNode> + '_ {
        subgraph
            .filter
            .iter_ones()
            .map(|i| self.involution.get_node_id(i))
    }

    pub fn iter_node_data<'a>(&'a self, subgraph: &'a SubGraph) -> NodeIterator<'a, E, V> {
        NodeIterator {
            graph: self,
            edges: subgraph.filter.iter_ones(),
            seen: bitvec![usize, Lsb0; 0; self.base_nodes],
        }
    }

    pub fn get_node_pos(&self, node: &NestingNode) -> usize {
        self.nodes.get_index_of(node).unwrap()
    }

    pub fn get_incident_node_id(&self, edge: usize) -> &NestingNode {
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

    pub fn random(nodes: usize, edges: usize, seed: u64) -> NestingGraph<(), ()> {
        let mut inv: Involution<Option<NestingNode>, ()> =
            Involution::<NestingNode, ()>::random(edges, seed);

        let mut rng = Xoroshiro64Star::seed_from_u64(seed);

        let mut externals = Vec::new();
        let mut sources = Vec::new();
        let mut sinks = Vec::new();

        for (i, e) in inv.inv.iter().enumerate() {
            let nodeid = NestingNode::node_from_pos(&[i], inv.inv.len());
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
            for i in sink.externalhedges.iter_ones() {
                inv.set_id(i, Some(sink.clone()));
            }
        }

        for source in &sources {
            for i in source.externalhedges.iter_ones() {
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

        NestingGraph {
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

    fn nesting_node_from_subgraph(&self, internal_graph: SubGraph) -> NestingNode {
        let mut externalhedges = bitvec![usize, Lsb0; 0; self.involution.len()];

        for i in internal_graph.filter.iter_ones() {
            externalhedges |= &self.involution.inv[i].0.externalhedges;
        }

        NestingNode {
            externalhedges: !(!externalhedges | &internal_graph.filter),
            internal_graph,
        }
    }

    fn nesting_node_fix(&self, node: &mut NestingNode) {
        let mut externalhedges = bitvec![usize, Lsb0; 0; self.involution.len()];

        for i in node.internal_graph.filter.iter_ones() {
            externalhedges |= &self.involution.inv[i].0.externalhedges;
        }

        node.externalhedges = !(!externalhedges | &node.internal_graph.filter);
    }

    fn paired_filter_from_pos(&self, pos: &[usize]) -> BitVec {
        let mut filter = bitvec![usize, Lsb0; 0; self.involution.len()];

        for &i in pos {
            filter.set(i, true);
            filter.set(self.involution.inv(i), true);
        }

        filter
    }

    fn filter_from_pos(&self, pos: &[usize]) -> BitVec {
        NestingNode::filter_from_pos(pos, self.involution.len())
    }

    fn nesting_node_from_pos(&self, pos: &[usize]) -> NestingNode {
        self.nesting_node_from_subgraph(SubGraph::from(self.filter_from_pos(pos)))
    }

    fn remove_externals(&self, subgraph: &mut NestingNode) {
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

    fn full_filter(&self) -> BitVec {
        bitvec![usize, Lsb0; 1; self.involution.len()]
    }

    fn empty_filter(&self) -> BitVec {
        bitvec![usize, Lsb0; 0; self.involution.len()]
    }

    fn full_node(&self) -> NestingNode {
        NestingNode {
            internal_graph: (self.full_filter() & !self.external_filter()).into(),
            externalhedges: self.external_filter(),
        }
    }

    fn cycle_basis(&self) -> Vec<NestingNode> {
        let i = self.involution.first_internal().unwrap();

        self.paton_cycle_basis(&self.full_filter().into(), i)
            .unwrap()
    }

    ///Read, R.C. and Tarjan, R.E. (1975), Bounds on Backtrack Algorithms for Listing Cycles, Paths, and Spanning Trees. Networks, 5: 237-252. https://doi.org/10.1002/net.1975.5.3.237

    fn read_tarjan(&self) -> Vec<NestingNode> {
        let mut cycles = Vec::new();
        for (i, (_, e)) in self.involution.inv.iter().enumerate() {
            if matches!(e, InvolutiveMapping::Source(_)) {
                cycles.extend(self.increasing_cycles_from(i));
            }
        }
        cycles
    }

    fn backtrack(
        &self,
        options: &mut SubGraph,
        tree: &mut SubGraph,
        current: usize,
    ) -> Option<usize> {
        tree.filter.set(self.involution.inv(current), false);
        tree.filter.set(current, false);

        let current_node = &self.involution.get_node_id(current).externalhedges;

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

    fn increasing_cycles_from(&self, start: usize) -> Vec<NestingNode> {
        let mut cycles = Vec::new();

        let mut tree: SubGraph = self.empty_filter().into();
        let mut options: SubGraph = SubGraph::from(self.external_filter()).complement();

        for i in 0..=start {
            options.filter.set(i, false);
        }

        tree.filter.set(start, true);
        options.filter.set(start, false);

        let mut current = start;

        loop {
            let next_possible = SubGraph::from(
                self.involution
                    .get_connected_node_id(current)
                    .unwrap()
                    .externalhedges
                    .clone(),
            );

            if next_possible.empty_intersection(&tree) {
                tree.filter.set(self.involution.inv(current), true);
                options.filter.set(self.involution.inv(current), false);

                if let Some(i) = next_possible
                    .filter
                    .iter_ones()
                    .find(|&i| i != self.involution.inv(start) && options.filter[i])
                {
                    tree.filter.set(i, true);
                    options.filter.set(i, false);

                    current = i;
                } else if current == start {
                    //No edges greater than start-> abort
                    break;
                } else {
                    current = self
                        .backtrack(&mut options, &mut tree, current)
                        .unwrap_or(start);

                    if current == start {
                        break;
                    }
                }
            } else {
                let int = next_possible.intersection(&tree);
                if int.filter[start] {
                    tree.filter.set(self.involution.inv(current), true);
                    options.filter.set(self.involution.inv(current), false);

                    let mut cycle = self.nesting_node_from_subgraph(tree.clone());

                    self.cut_branches(&mut cycle);

                    cycle.internal_graph.loopcount = Some(1);

                    cycles.push(cycle);
                }
                current = self
                    .backtrack(&mut options, &mut tree, current)
                    .unwrap_or(start);

                if current == start {
                    break;
                }
            }
        }

        cycles
    }

    fn order_basis(&self, basis: &[NestingNode]) -> Vec<Vec<SubGraph>> {
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

    fn welchs_violating(partitions: &[Vec<SubGraph>]) -> Vec<SubGraph> {
        let mut violations = Vec::new();

        for (b1, b2) in partitions.iter().tuple_windows() {
            for (_, _, ck) in
                b1.iter()
                    .chain(b2.iter())
                    .tuple_combinations()
                    .filter(|(ci, cj, ck)| {
                        ci.empty_intersection(cj)
                            && !ci.empty_intersection(ck)
                            && !cj.empty_intersection(ck)
                    })
            {
                violations.push(ck.clone());
            }
        }

        violations
    }

    fn all_cycles(&self) -> Vec<NestingNode> {
        self.all_composite_cycles_with_basis(&self.cycle_basis())
    }

    fn all_cycle_unions(&self) -> AHashSet<SubGraph> {
        let cycles = self.read_tarjan();
        let mut spinneys = AHashSet::new();

        let pset = PowersetIterator::new(cycles.len() as u8);

        for p in pset {
            // let union = p
            //     .iter_ones()
            //     .map(|i| cycles[i].clone())
            //     .reduce(|acc, n| acc.union(&n));

            // if let Some(union) = union {
            //     spinneys.push(union);
            // }
            let mut union: Option<SubGraph> = None;

            for i in p.iter_ones() {
                if let Some(union) = &mut union {
                    union.union_with(&cycles[i].internal_graph);
                } else {
                    union = Some(cycles[i].internal_graph.clone());
                }
            }

            if let Some(union) = union {
                spinneys.insert(union);
            }
        }

        spinneys
    }

    fn all_composite_cycles_with_basis(&self, basis: &[NestingNode]) -> Vec<NestingNode> {
        let mut cycles = Vec::new();

        // println!("loops: {}", basis.len());

        let mut pset = PowersetIterator::new(basis.len() as u8);

        pset.next(); //Skip empty set

        for p in pset {
            let mut cycle: SubGraph = self.empty_filter().into();
            for c in p.iter_ones().map(|i| &basis[i]) {
                cycle.sym_diff_with(&c.internal_graph);
            }

            if self.is_cycle(&cycle) {
                let cycle = self.nesting_node_from_subgraph(cycle);
                cycles.push(cycle);
            }
        }

        cycles
    }

    fn is_cycle(&self, subgraph: &SubGraph) -> bool {
        let mut is = true;
        for e in self.iter_egde_node(subgraph) {
            let adgacent = subgraph.filter.clone() & &e.externalhedges;
            if adgacent.count_ones() != 2 {
                is = false;
                break;
            }
        }
        is
    }

    fn paton_count_loops(
        &self,
        subgraph: &SubGraph,
        start: usize,
    ) -> Result<usize, NestingNodeError> {
        if !subgraph.filter[start] {
            return Err(NestingNodeError::InvalidStart);
        }
        let mut loopcount = 0;

        let n = self.involution.get_node_id(start);

        let mut tree: SubGraph = n.clone().externalhedges.into();

        let mut left: SubGraph = (self.full_filter() & !self.external_filter()).into();

        loop {
            let z = tree
                .filter
                .iter()
                .enumerate()
                .find(|(i, x)| *x.as_ref() & left.filter[*i])
                .map(|(i, _)| i);

            if let Some(z) = z {
                left.filter.set(z, false); // has been visited
                left.filter.set(self.involution.inv(z), false);

                let w: SubGraph = self
                    .involution
                    .get_connected_node_id(z)
                    .unwrap()
                    .externalhedges
                    .clone()
                    .into();

                if w.empty_intersection(&tree) {
                    // w is not yet in the tree
                    tree.union_with(&w);
                } else {
                    loopcount += 1;
                    tree.filter.set(self.involution.inv(z), false);
                    tree.filter.set(z, false);
                }
            } else {
                break;
            }
        }

        Ok(loopcount)
    }

    fn paton_cycle_basis(
        &self,
        subgraph: &SubGraph,
        start: usize,
    ) -> Result<Vec<NestingNode>, NestingNodeError> {
        if !subgraph.filter[start] {
            return Err(NestingNodeError::InvalidStart);
        }
        let mut cycle_basis = Vec::new();

        let n = self.involution.get_node_id(start);

        let mut tree: SubGraph = n.clone().externalhedges.into();

        let mut left: SubGraph = (self.full_filter() & !self.external_filter()).into();

        loop {
            let z = tree
                .filter
                .iter()
                .enumerate()
                .find(|(i, x)| *x.as_ref() & left.filter[*i])
                .map(|(i, _)| i);

            if let Some(z) = z {
                left.filter.set(z, false); // has been visited
                left.filter.set(self.involution.inv(z), false);

                let w: SubGraph = self
                    .involution
                    .get_connected_node_id(z)
                    .unwrap()
                    .externalhedges
                    .clone()
                    .into();

                if w.empty_intersection(&tree) {
                    // w is not yet in the tree
                    tree.union_with(&w);
                } else {
                    let mut cycle = left.complement();
                    cycle.intersect_with(&tree);
                    let mut cycle = self.nesting_node_from_subgraph(cycle);
                    cycle.internal_graph.filter.set(z, true);
                    cycle
                        .internal_graph
                        .filter
                        .set(self.involution.inv(z), true);

                    cycle.internal_graph.loopcount = Some(1);

                    self.cut_branches(&mut cycle);

                    cycle_basis.push(cycle);

                    tree.filter.set(self.involution.inv(z), false);
                    tree.filter.set(z, false);
                }
            } else {
                break;
            }
        }

        Ok(cycle_basis)
    }

    pub fn dot(&self, node_as_graph: &NestingNode) -> String {
        let mut out = "graph {\n ".to_string();
        out.push_str("  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\";\n ");

        for (hedge_id, (incident_node, edge)) in self.involution.inv.iter().enumerate() {
            match &edge {
                //Internal graphs never have unpaired edges
                InvolutiveMapping::Identity(_) => {
                    let attr = if *node_as_graph.internal_graph.filter.get(hedge_id).unwrap() {
                        None
                    } else if *node_as_graph.externalhedges.get(hedge_id).unwrap() {
                        Some(GVEdgeAttrs {
                            color: Some("gray50".to_string()),
                            label: None,
                        })
                    } else {
                        Some(GVEdgeAttrs {
                            color: Some("gray75".to_string()),
                            label: None,
                        })
                    };
                    out.push_str(&InvolutiveMapping::<()>::identity_dot(
                        hedge_id,
                        self.nodes.get_index_of(incident_node).unwrap(),
                        attr.as_ref(),
                    ));
                }
                InvolutiveMapping::Source((_, _)) => {
                    let attr = if *node_as_graph.internal_graph.filter.get(hedge_id).unwrap() {
                        None
                    } else if *node_as_graph.externalhedges.get(hedge_id).unwrap()
                        && !*node_as_graph
                            .externalhedges
                            .get(self.involution.inv(hedge_id))
                            .unwrap()
                    {
                        Some(GVEdgeAttrs {
                            color: Some("gray50:gray75;0.5".to_string()),
                            label: None,
                        })
                    } else if *node_as_graph
                        .externalhedges
                        .get(self.involution.inv(hedge_id))
                        .unwrap()
                        && !*node_as_graph.externalhedges.get(hedge_id).unwrap()
                    {
                        Some(GVEdgeAttrs {
                            color: Some("gray75:gray50;0.5".to_string()),
                            label: None,
                        })
                    } else if *node_as_graph
                        .externalhedges
                        .get(self.involution.inv(hedge_id))
                        .unwrap()
                        && *node_as_graph.externalhedges.get(hedge_id).unwrap()
                    {
                        Some(GVEdgeAttrs {
                            color: Some("gray50".to_string()),
                            label: None,
                        })
                    } else {
                        Some(GVEdgeAttrs {
                            color: Some("gray75".to_string()),
                            label: None,
                        })
                    };
                    out.push_str(&InvolutiveMapping::<()>::pair_dot(
                        self.nodes.get_index_of(incident_node).unwrap(),
                        self.nodes
                            .get_index_of(self.involution.get_connected_node_id(hedge_id).unwrap())
                            .unwrap(),
                        attr.as_ref(),
                    ));
                }
                InvolutiveMapping::Sink(i) => {
                    if *node_as_graph.internal_graph.filter.get(hedge_id).unwrap()
                        && !*node_as_graph.internal_graph.filter.get(*i).unwrap()
                    {
                        panic!("Internal graph has a dangling sink edge")
                    }
                }
            }
        }

        out += "}";
        out
    }

    fn cut_branches(&self, subgraph: &mut NestingNode) {
        let nodes = AHashSet::<&NestingNode>::from_iter(
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
                let int = n.externalhedges.clone() & &subgraph.internal_graph.filter;

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

    pub fn nesting_poset(&self) -> PoSet<NestingNode> {
        let iter = self.nodes.keys().cloned();

        PoSet::from_iter(iter)
    }

    fn all_spinneys_with_basis(&self, basis: &[&SubGraph]) -> AHashSet<NestingNode> {
        let mut spinneys = AHashSet::new();
        let mut base_cycle: SubGraph = self.empty_filter().into();

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

    fn all_spinneys(&self) -> AHashSet<SubGraph> {
        let mut cycles = self.cycle_basis();

        let mut all_combinations = PowersetIterator::new(cycles.len() as u8);
        all_combinations.next(); //Skip empty set

        let mut spinneys = AHashSet::new();

        for p in all_combinations {
            let mut base_cycle: SubGraph = self.empty_filter().into();

            for i in p.iter_ones() {
                base_cycle.sym_diff_with(&cycles[i].internal_graph);
            }

            cycles.push(self.nesting_node_from_subgraph(base_cycle.clone()));
        }

        for (ci, cj) in cycles.iter().tuple_combinations() {
            spinneys.insert(ci.internal_graph.union(&cj.internal_graph));
        }

        for c in cycles {
            spinneys.insert(c.internal_graph);
        }
        spinneys
    }

    fn all_spinneys_alt(&self) -> AHashSet<SubGraph> {
        let mut spinneys = AHashSet::new();
        let mut cycles = self.all_cycles();

        let mut pset = PowersetIterator::new(cycles.len() as u8);
        pset.next(); //Skip empty set

        for p in pset {
            let mut union: SubGraph = self.empty_filter().into();

            for i in p.iter_ones() {
                union.union_with(&cycles[i].internal_graph);
            }

            spinneys.insert(union);
        }

        for c in cycles {
            spinneys.insert(c.internal_graph);
        }
        spinneys
    }
}

pub struct PoSet<N>
where
    N: PartialOrd,
{
    greater_than: Vec<Vec<usize>>,
    nodes: Vec<N>,
}

pub struct CoverSet<N>
where
    N: PartialOrd,
{
    covers: Vec<Vec<usize>>,
    nodes: Vec<N>,
}

impl<N: PartialOrd> CoverSet<N> {
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

        out += "}";
        out
    }

    pub fn greater_than(&self, a: usize) -> Vec<usize> {
        self.greater_than[a].clone()
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
}

impl UVEdge {
    pub fn from_edge(edge: &Edge, id: usize) -> Self {
        UVEdge { og_edge: id }
    }

    pub fn numerator_rank(&self) -> isize {
        0
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
struct UVNode {}

impl UVNode {
    pub fn from_vertex(vertex: &Vertex) -> Self {
        UVNode {}
    }

    pub fn numerator_rank(&self) -> isize {
        0
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct UVGraph(NestingGraph<UVEdge, UVNode>);

#[allow(dead_code)]
impl UVGraph {
    pub fn from_graph(graph: &Graph) -> Self {
        let mut uv_graph = NestingGraphBuilder::new();

        for n in &graph.vertices {
            uv_graph.add_node(UVNode::from_vertex(n));
        }

        for (i, edge) in graph.edges.iter().enumerate() {
            let ves = edge.vertices;

            let sink = NodeIndex(ves[1]);
            let source = NodeIndex(ves[0]);

            match edge.edge_type {
                EdgeType::Virtual => {
                    uv_graph.add_edge(source, sink, UVEdge::from_edge(edge, i));
                }
                EdgeType::Outgoing => {
                    uv_graph.add_external_edge(source, UVEdge::from_edge(edge, i));
                }
                EdgeType::Incoming => {
                    uv_graph.add_external_edge(sink, UVEdge::from_edge(edge, i));
                }
            }
        }

        UVGraph(uv_graph.into())
    }

    pub fn half_edge_id(&self, g: &Graph, id: usize) -> usize {
        self.0
            .involution
            .find_from_data(&UVEdge::from_edge(&g.edges[id], id))
            .unwrap()
    }

    pub fn node_id(&self, g: &Graph, id: usize) -> NestingNode {
        let e_id = g.vertices[id].edges.first().unwrap();
        self.0
            .involution
            .get_node_id(self.half_edge_id(g, *e_id))
            .clone()
    }

    pub fn cycle_basis_from_lmb(&self, lmb: &LoopMomentumBasis) -> Vec<NestingNode> {
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

    fn spanning_forest_from_lmb(&self, lmb: LoopMomentumBasis) -> NestingNode {
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

    fn wood(&self) -> Vec<SubGraph> {
        let cycles = self.0.read_tarjan();
        let mut spinneys = Vec::with_capacity(1 << cycles.len());

        let pset = PowersetIterator::new(cycles.len() as u8);

        for p in pset {
            let mut union: Option<SubGraph> = None;

            for i in p.iter_ones() {
                if let Some(union) = &mut union {
                    union.union_with(&cycles[i].internal_graph);
                } else {
                    union = Some(cycles[i].internal_graph.clone());
                }
            }

            if let Some(union) = union {
                if self.dod(&union) >= 0 {
                    spinneys.push(union);
                }
            }
        }

        spinneys
    }

    fn dod(&self, subgraph: &SubGraph) -> isize {
        let loop_count = if let Some(loop_count) = subgraph.loopcount {
            loop_count
        } else {
            self.0
                .paton_count_loops(subgraph, subgraph.filter.first_one().unwrap())
                .unwrap()
        };

        isize::try_from(4 * loop_count).unwrap() + isize::try_from(self.0.n_externals()).unwrap()
            - isize::try_from(subgraph.filter.count_ones()).unwrap()
            + self.numerator_rank(subgraph)
    }
}

#[cfg(test)]
mod tests;
