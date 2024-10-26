use std::hash::Hash;

use ahash::{AHashMap, AHashSet};
use bitvec::{slice::IterOnes, vec::BitVec};
use indexmap::IndexMap;
use itertools::Itertools;
use rand::{rngs::SmallRng, Rng, SeedableRng};
use serde::{Deserialize, Serialize};

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
    // Undirected(usize),
    // UndirectedData((Option<E>, usize)),
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

    fn n_internals(&self, subgraph: &SubGraph) -> usize {
        subgraph
            .filter
            .iter_ones()
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

    fn is_internal(&self, index: usize, subgraph: &SubGraph) -> bool {
        if !subgraph.filter[index] {
            return false;
        }
        match &self.inv[index].1 {
            InvolutiveMapping::Identity(_) => false,
            InvolutiveMapping::Source((_, i)) => subgraph.filter[*i],
            InvolutiveMapping::Sink(i) => subgraph.filter[*i],
        }
    }

    fn get_smart_data(&self, index: usize, subgraph: &SubGraph) -> Option<&E> {
        //should be a better enum (none,internal,external)
        if subgraph.filter[index] {
            match &self.inv[index].1 {
                InvolutiveMapping::Identity(data) => Some(data.as_ref().unwrap()),
                InvolutiveMapping::Source((data, _)) => Some(data.as_ref().unwrap()),
                InvolutiveMapping::Sink(i) => {
                    if subgraph.filter[*i] {
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
struct GVEdgeAttrs {
    label: Option<String>,
    color: Option<String>,
    other: Option<String>,
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
        if let Some(other) = &self.other {
            out.push_str(&format!("{},", other));
        }
        out.push(']');
        write!(f, "{}", out)
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct HedgeGraph<E, V> {
    nodes: IndexMap<HalfEdgeNode, V>, // Forest of nodes, that contain half-edges, with data E
    base_nodes: usize,
    involution: Involution<HalfEdgeNode, E>, // Involution of half-edges
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
    pub fn to_nesting_node<E, V>(&self, graph: &HedgeGraph<E, V>) -> HalfEdgeNode {
        graph.nesting_node_from_subgraph(self.clone())
    }
    pub fn is_empty(&self) -> bool {
        self.filter.count_ones() == 0
    }

    pub fn set_loopcount<E, V>(&mut self, graph: &HedgeGraph<E, V>) {
        self.loopcount = Some(graph.cyclotomatic_number(self));
    }

    pub fn cycle_basis<E, V>(&self, graph: &HedgeGraph<E, V>) -> Vec<HalfEdgeNode> {
        graph
            .paton_cycle_basis(self, self.filter.first_one().unwrap())
            .unwrap()
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
pub struct HalfEdgeNode {
    internal_graph: SubGraph,
    externalhedges: BitVec, // essentially all hedges that are not in the internal graph, but are connected to it
}

impl PartialOrd for HalfEdgeNode {
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

impl HalfEdgeNode {
    pub fn new(len: usize) -> Self {
        HalfEdgeNode {
            internal_graph: bitvec![usize, Lsb0; 0; len].into(),
            externalhedges: bitvec![usize, Lsb0; 0; len],
        }
    }

    fn all_edges(&self) -> BitVec {
        self.internal_graph.filter.clone() | &self.externalhedges
    }

    pub fn union_with(&mut self, other: &HalfEdgeNode) {
        self.internal_graph.filter |= &other.internal_graph.filter;
        self.externalhedges |= &other.externalhedges;
        self.externalhedges &= !self.internal_graph.filter.clone();
    }

    pub fn intersect_with(&mut self, other: &HalfEdgeNode) {
        self.internal_graph.filter &= &other.internal_graph.filter;
        self.externalhedges &= &other.externalhedges;
    }

    #[must_use]
    pub fn intersection(&self, other: &HalfEdgeNode) -> HalfEdgeNode {
        let mut new = self.clone();
        new.intersect_with(other);
        new
    }

    #[must_use]
    pub fn union(&self, other: &HalfEdgeNode) -> HalfEdgeNode {
        let mut new = self.clone();
        new.union_with(other);
        new
    }

    pub fn empty_intersection(&self, other: &HalfEdgeNode) -> bool {
        let internals = self.internal_graph.filter.clone() & &other.internal_graph.filter;

        let externals = self.externalhedges.clone() & &other.externalhedges;

        internals.count_ones() == 0 && externals.count_ones() == 0
    }

    pub fn weakly_disjoint(&self, other: &HalfEdgeNode) -> bool {
        let internals = self.internal_graph.filter.clone() & &other.internal_graph.filter;

        internals.count_ones() == 0
    }

    pub fn strongly_disjoint(&self, other: &HalfEdgeNode) -> bool {
        let internals = self.internal_graph.filter.clone() & &other.internal_graph.filter;

        let externals = self.externalhedges.clone() & &other.externalhedges;

        internals.count_ones() == 0 && externals.count_ones() == 0
    }

    pub fn node_from_pos(pos: &[usize], len: usize) -> HalfEdgeNode {
        HalfEdgeNode {
            externalhedges: HalfEdgeNode::filter_from_pos(pos, len),
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

    pub fn internal_graph_union(&self, other: &HalfEdgeNode) -> SubGraph {
        SubGraph {
            filter: self.internal_graph.filter.clone() | &other.internal_graph.filter,
            loopcount: None,
        }
    }

    pub fn from_builder<V>(builder: &HedgeNodeBuilder<V>, len: usize) -> Self {
        let internal_graph = bitvec![usize, Lsb0; 0; len].into();
        let mut externalhedges = bitvec![usize, Lsb0; 0; len];

        for hedge in &builder.hedges {
            let mut bit = externalhedges.get_mut(*hedge).unwrap();
            *bit = true;
        }

        HalfEdgeNode {
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
pub struct HedgeNodeBuilder<V> {
    data: V,
    hedges: Vec<usize>,
}

pub struct HedgeGraphBuidler<E, V> {
    nodes: Vec<HedgeNodeBuilder<V>>,
    involution: Involution<NodeIndex, E>,
}

pub struct NodeIterator<'a, E, V> {
    graph: &'a HedgeGraph<E, V>,
    edges: IterOnes<'a, usize, Lsb0>,
    seen: BitVec,
}

impl<'a, E, V> Iterator for NodeIterator<'a, E, V> {
    type Item = (&'a HalfEdgeNode, &'a V);

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

impl<E, V> HedgeGraphBuidler<E, V> {
    pub fn new() -> Self {
        HedgeGraphBuidler {
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

impl<E, V> Default for HedgeGraphBuidler<E, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<E, V> From<HedgeGraphBuidler<E, V>> for HedgeGraph<E, V> {
    fn from(builder: HedgeGraphBuidler<E, V>) -> Self {
        let len = builder.involution.len();
        let involution = Involution {
            inv: builder
                .involution
                .inv
                .into_iter()
                .map(|(n, i)| (HalfEdgeNode::from_builder(&builder.nodes[n.0], len), i))
                .collect(),
        };
        let nodes: IndexMap<HalfEdgeNode, V> = builder
            .nodes
            .into_iter()
            .map(|x| (HalfEdgeNode::from_builder(&x, len), x.data))
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
        let mut builder = HedgeGraphBuidler::new();
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

impl From<BareGraph> for HedgeGraph<usize, usize> {
    fn from(value: BareGraph) -> Self {
        let mut builder = HedgeGraphBuidler::new();
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

use thiserror::Error;

use super::BareGraph;

#[derive(Error, Debug)]

pub enum NestingNodeError {
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

    pub fn iter_egde_data<'a>(&'a self, subgraph: &'a SubGraph) -> impl Iterator<Item = &E> + '_ {
        subgraph
            .filter
            .iter_ones()
            .map(|i| self.involution.get_data(i))
    }

    pub fn iter_internal_edge_data<'a>(
        &'a self,
        subgraph: &'a SubGraph,
    ) -> impl Iterator<Item = &E> + '_ {
        subgraph
            .filter
            .iter_ones()
            .flat_map(|i| self.involution.get_smart_data(i, subgraph))
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

    pub fn count_connected_components(&self, subgraph: &SubGraph) -> usize {
        let mut visited_edges = bitvec![usize, Lsb0; 0; self.n_hedges()];
        let mut component_count = 0;

        // Iterate over all edges in the subgraph
        for hedge_index in subgraph.filter.iter_ones() {
            if !visited_edges[hedge_index] {
                // Start a new component
                component_count += 1;

                // Perform DFS to find all reachable edges from this edge
                let reachable_edges = self.dfs_reach(subgraph, hedge_index);

                // Mark all reachable edges as visited
                for edge in reachable_edges {
                    visited_edges.set(edge, true);
                }
            }
        }

        component_count
    }

    pub fn cyclotomatic_number(&self, subgraph: &SubGraph) -> usize {
        let n_hedges = self.count_internal_edges(subgraph);
        let n_nodes = self.number_of_nodes_in_subgraph(subgraph);
        let n_components = self.count_connected_components(subgraph);

        n_hedges - n_nodes + n_components
    }

    pub fn count_internal_edges(&self, subgraph: &SubGraph) -> usize {
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
    ) -> impl Iterator<Item = &HalfEdgeNode> + '_ {
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

    pub fn get_node_pos(&self, node: &HalfEdgeNode) -> usize {
        self.nodes.get_index_of(node).unwrap()
    }

    pub fn get_incident_node_id(&self, edge: usize) -> &HalfEdgeNode {
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
        let mut inv: Involution<Option<HalfEdgeNode>, ()> =
            Involution::<HalfEdgeNode, ()>::random(edges, seed);

        let mut rng = SmallRng::seed_from_u64(seed);

        let mut externals = Vec::new();
        let mut sources = Vec::new();
        let mut sinks = Vec::new();

        for (i, e) in inv.inv.iter().enumerate() {
            let nodeid = HalfEdgeNode::node_from_pos(&[i], inv.inv.len());
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

    fn nesting_node_from_subgraph(&self, internal_graph: SubGraph) -> HalfEdgeNode {
        let mut externalhedges = bitvec![usize, Lsb0; 0; self.involution.len()];

        for i in internal_graph.filter.iter_ones() {
            externalhedges |= &self.involution.inv[i].0.externalhedges;
        }

        HalfEdgeNode {
            externalhedges: !(!externalhedges | &internal_graph.filter),
            internal_graph,
        }
    }

    fn nesting_node_fix(&self, node: &mut HalfEdgeNode) {
        let mut externalhedges = bitvec![usize, Lsb0; 0; self.involution.len()];

        for i in node.internal_graph.filter.iter_ones() {
            externalhedges |= &self.involution.inv[i].0.externalhedges;
        }

        node.externalhedges = !(!externalhedges | &node.internal_graph.filter);
    }

    pub fn paired_filter_from_pos(&self, pos: &[usize]) -> BitVec {
        let mut filter = bitvec![usize, Lsb0; 0; self.involution.len()];

        for &i in pos {
            filter.set(i, true);
            filter.set(self.involution.inv(i), true);
        }

        filter
    }

    pub fn filter_from_pos(&self, pos: &[usize]) -> BitVec {
        HalfEdgeNode::filter_from_pos(pos, self.involution.len())
    }

    pub fn nesting_node_from_pos(&self, pos: &[usize]) -> HalfEdgeNode {
        self.nesting_node_from_subgraph(SubGraph::from(self.filter_from_pos(pos)))
    }

    fn remove_externals(&self, subgraph: &mut HalfEdgeNode) {
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

    fn full_node(&self) -> HalfEdgeNode {
        HalfEdgeNode {
            internal_graph: (self.full_filter() & !self.external_filter()).into(),
            externalhedges: self.external_filter(),
        }
    }

    fn cycle_basis(&self) -> Vec<HalfEdgeNode> {
        let i = self.involution.first_internal().unwrap();

        self.paton_cycle_basis(&self.full_filter().into(), i)
            .unwrap()
    }

    ///Read, R.C. and Tarjan, R.E. (1975), Bounds on Backtrack Algorithms for Listing Cycles, Paths, and Spanning Trees. Networks, 5: 237-252. https://doi.org/10.1002/net.1975.5.3.237

    pub fn read_tarjan(&self) -> Vec<HalfEdgeNode> {
        let mut cycles = Vec::new();
        for (i, (_, e)) in self.involution.inv.iter().enumerate() {
            if matches!(e, InvolutiveMapping::Source(_)) {
                cycles.extend(self.increasing_cycles_from(i));
            }
        }
        cycles
    }

    pub fn backtrack(
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

    pub fn increasing_cycles_from(&self, start: usize) -> Vec<HalfEdgeNode> {
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

    pub fn order_basis(&self, basis: &[HalfEdgeNode]) -> Vec<Vec<SubGraph>> {
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

    pub fn welchs_violating(partitions: &[Vec<SubGraph>]) -> Vec<SubGraph> {
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

    pub fn all_cycles(&self) -> Vec<HalfEdgeNode> {
        self.all_composite_cycles_with_basis(&self.cycle_basis())
    }

    pub fn all_cycle_unions(&self) -> AHashSet<SubGraph> {
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

    pub fn all_basis_sym_diffs(&self, basis: &[HalfEdgeNode]) -> Vec<HalfEdgeNode> {
        let mut diffs = Vec::new();

        let mut pset = PowersetIterator::new(basis.len() as u8);

        pset.next(); //Skip empty set

        for p in pset {
            let mut cycle: SubGraph = self.empty_filter().into();
            for c in p.iter_ones().map(|i| &basis[i]) {
                cycle.sym_diff_with(&c.internal_graph);
            }

            let cycle = self.nesting_node_from_subgraph(cycle);
            diffs.push(cycle);
        }

        diffs
    }

    pub fn all_composite_cycles_with_basis(&self, basis: &[HalfEdgeNode]) -> Vec<HalfEdgeNode> {
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

    pub fn full_graph(&self) -> SubGraph {
        self.full_node().internal_graph
    }

    pub fn paton_count_loops(
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

    pub fn number_of_nodes_in_subgraph(&self, subgraph: &SubGraph) -> usize {
        self.iter_node_data(subgraph).count()
    }

    pub fn node_degrees_in_subgraph(&self, subgraph: &SubGraph) -> AHashMap<usize, usize> {
        let mut degrees = AHashMap::new();

        for (node, _) in self.iter_node_data(subgraph) {
            let node_pos = self.get_node_pos(node);

            // Count the number of edges in the subgraph incident to this node
            let incident_edges = node.externalhedges.clone() & &subgraph.filter;
            let degree = incident_edges.count_ones();

            degrees.insert(node_pos, degree);
        }

        degrees
    }

    fn paton_cycle_basis(
        &self,
        subgraph: &SubGraph,
        start: usize,
    ) -> Result<Vec<HalfEdgeNode>, NestingNodeError> {
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

    pub fn dot_impl(
        &self,
        node_as_graph: &HalfEdgeNode,
        graph_info: String,
        edge_attr: &impl Fn(&E) -> String,
        node_attr: &impl Fn(&V) -> String,
    ) -> String {
        let mut out = "graph {\n ".to_string();
        out.push_str("  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\";\n ");

        out.push_str(graph_info.as_str());

        for (n, v) in self.iter_node_data(&node_as_graph.internal_graph) {
            out.push_str(
                format!(
                    "  {} [{}];\n",
                    self.nodes.get_index_of(n).unwrap(),
                    node_attr(v)
                )
                .as_str(),
            );
        }

        for (hedge_id, (incident_node, edge)) in self.involution.inv.iter().enumerate() {
            match &edge {
                //Internal graphs never have unpaired edges
                InvolutiveMapping::Identity(data) => {
                    let attr = if *node_as_graph.internal_graph.filter.get(hedge_id).unwrap() {
                        None
                    } else if *node_as_graph.externalhedges.get(hedge_id).unwrap() {
                        Some(GVEdgeAttrs {
                            color: Some("gray50".to_string()),
                            label: None,
                            other: data.as_ref().map(edge_attr),
                        })
                    } else {
                        Some(GVEdgeAttrs {
                            color: Some("gray75".to_string()),
                            label: None,
                            other: data.as_ref().map(edge_attr),
                        })
                    };
                    out.push_str(&InvolutiveMapping::<()>::identity_dot(
                        hedge_id,
                        self.nodes.get_index_of(incident_node).unwrap(),
                        attr.as_ref(),
                    ));
                }
                InvolutiveMapping::Source((data, _)) => {
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
                            other: data.as_ref().map(edge_attr),
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
                            other: data.as_ref().map(edge_attr),
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
                            other: data.as_ref().map(edge_attr),
                        })
                    } else {
                        Some(GVEdgeAttrs {
                            color: Some("gray75".to_string()),
                            label: None,
                            other: data.as_ref().map(edge_attr),
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
    pub fn dot(&self, node_as_graph: &HalfEdgeNode) -> String {
        self.dot_impl(node_as_graph, "".to_string(), &|_| "".to_string(), &|_| {
            "".to_string()
        })
    }

    fn cut_branches(&self, subgraph: &mut HalfEdgeNode) {
        let nodes = AHashSet::<&HalfEdgeNode>::from_iter(
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

    pub fn all_spinneys_with_basis(&self, basis: &[&SubGraph]) -> AHashSet<HalfEdgeNode> {
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

    pub fn all_spinneys_rec(
        &self,
        spinneys: &mut AHashSet<HalfEdgeNode>,
        cycle_sums: Vec<HalfEdgeNode>,
    ) {
        let _len = spinneys.len();

        let mut pset = PowersetIterator::new(cycle_sums.len() as u8);

        pset.next(); //Skip empty set

        for (ci, cj) in cycle_sums.iter().tuple_combinations() {
            let _union = ci.internal_graph.union(&cj.internal_graph);

            // spinneys.insert(union);
        }
    }

    pub fn all_spinneys(&self) -> AHashMap<SubGraph, Vec<(HalfEdgeNode, Option<HalfEdgeNode>)>> {
        let mut cycles = self.cycle_basis();

        let mut all_combinations = PowersetIterator::new(cycles.len() as u8);
        all_combinations.next(); //Skip empty set

        let mut spinneys: AHashMap<SubGraph, Vec<(HalfEdgeNode, Option<HalfEdgeNode>)>> =
            AHashMap::new();

        for p in all_combinations {
            let mut base_cycle: SubGraph = self.empty_filter().into();

            for i in p.iter_ones() {
                base_cycle.sym_diff_with(&cycles[i].internal_graph);
            }

            cycles.push(self.nesting_node_from_subgraph(base_cycle.clone()));
        }

        for (ci, cj) in cycles.iter().tuple_combinations() {
            let union = ci.internal_graph.union(&cj.internal_graph);

            if let Some(v) = spinneys.get_mut(&union) {
                v.push((ci.clone(), Some(cj.clone())));
            } else {
                spinneys.insert(union, vec![(ci.clone(), Some(cj.clone()))]);
            }
        }

        for c in cycles {
            spinneys.insert(c.internal_graph.clone(), vec![(c.clone(), None)]);
        }
        spinneys
    }

    pub fn all_spinneys_alt(&self) -> AHashSet<SubGraph> {
        let mut spinneys = AHashSet::new();
        let cycles = self.all_cycles();

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

#[cfg(test)]
mod test {

    #[test]
    fn threeloop() {
        let mut builder: HedgeGraphBuidler<(), ()> = HedgeGraphBuidler::new();
        let a = builder.add_node(());
        let b = builder.add_node(());
        let c = builder.add_node(());
        let d = builder.add_node(());

        builder.add_edge(a, b, ());
        builder.add_edge(b, a, ());
        builder.add_edge(a, b, ());

        builder.add_edge(b, c, ());
        builder.add_edge(c, d, ());
        builder.add_edge(d, a, ());

        let graph = builder.build();

        insta::assert_snapshot!("three_loop_dot", graph.base_dot());
        insta::assert_ron_snapshot!("three_loop", graph);

        for i in 0..graph.n_hedges() {
            assert_eq!(
                3,
                graph
                    .paton_cycle_basis(&graph.full_filter().into(), i)
                    .unwrap()
                    .len()
            );
        }

        let cycles = graph.cycle_basis();

        assert_eq!(3, cycles.len());

        let all_cycles = graph.read_tarjan();

        assert_eq!(6, all_cycles.len());

        insta::assert_ron_snapshot!("three_loop_cycles", cycles);
    }

    #[test]
    fn hairythreeloop() {
        let mut builder: HedgeGraphBuidler<(), ()> = HedgeGraphBuidler::new();
        let a = builder.add_node(());
        let b = builder.add_node(());
        let c = builder.add_node(());
        let d = builder.add_node(());

        builder.add_edge(a, b, ());
        builder.add_edge(b, a, ());
        builder.add_edge(a, b, ());
        builder.add_external_edge(a, ());
        builder.add_external_edge(b, ());
        builder.add_external_edge(b, ());

        builder.add_edge(b, c, ());
        builder.add_edge(c, d, ());
        builder.add_edge(d, a, ());

        assert_eq!(builder.involution.len(), 15);
        let graph = builder.build();

        insta::assert_snapshot!("hairy_three_loop_dot", graph.base_dot());
        insta::assert_ron_snapshot!("hairy_three_loop", graph);
        insta::assert_snapshot!(
            "hairy_three_loop_dot_internal",
            graph.dot(&graph.full_node())
        );

        for i in graph.full_node().internal_graph.filter.iter_ones() {
            assert_eq!(
                3,
                graph
                    .paton_cycle_basis(&graph.full_filter().into(), i)
                    .unwrap()
                    .len()
            );
        }

        let cycles = graph.cycle_basis();

        insta::assert_ron_snapshot!("hairy_three_loop_cycles", cycles);
    }

    #[test]
    fn cube() {
        let mut builder: HedgeGraphBuidler<(), ()> = HedgeGraphBuidler::new();
        let a = builder.add_node(());
        let b = builder.add_node(());
        let c = builder.add_node(());
        let d = builder.add_node(());
        let e = builder.add_node(());
        let f = builder.add_node(());
        let g = builder.add_node(());
        let h = builder.add_node(());

        builder.add_edge(a, b, ());
        builder.add_edge(b, c, ());
        builder.add_edge(c, d, ());
        builder.add_edge(d, a, ());

        builder.add_edge(e, f, ());
        builder.add_edge(f, g, ());
        builder.add_edge(g, h, ());
        builder.add_edge(h, e, ());

        builder.add_edge(a, e, ());
        builder.add_edge(b, f, ());
        builder.add_edge(c, g, ());
        builder.add_edge(d, h, ());

        let graph = builder.build();

        insta::assert_snapshot!("cube_dot", graph.base_dot());

        // let mut all_spinneys = graph.all_spinneys().into_iter().collect_vec();
        // all_spinneys.sort_by(|a, b| a.filter.cmp(&b.filter));

        // assert_eq!(162, all_spinneys.len());
    }

    #[test]
    #[ignore]
    fn alt_vs_pair() {
        for s in 0..100 {
            let rand_graph = HedgeGraph::<(), ()>::random(10, 14, s);

            let before = Instant::now();
            let all_spinneys = rand_graph.all_spinneys();
            let after = before.elapsed();
            let before = Instant::now();
            let all_spinneys_alt = rand_graph.all_spinneys_alt();
            let after_alt = before.elapsed();
            println!("{s} {after:?} {after_alt:?}");

            assert_eq!(
                all_spinneys.len(),
                all_spinneys_alt.len(),
                "{}",
                rand_graph.base_dot()
            );
        }
        // let rand_graph = NestingGraph::<(), ()>::random(6, 9, 8);

        // println!("{}", rand_graph.base_dot());

        // println!("loops {}", rand_graph.cycle_basis().len());

        // let all_spinneys_other = rand_graph.all_spinneys();

        // // println!("all spinneys read tarjan {}", all_spinneys.len());

        // println!("all spinneys {}", all_spinneys_other.len());
        // println!("all spinneys alt{}", rand_graph.all_spinneys_alt().len());
    }

    #[test]
    #[should_panic]
    fn read_tarjan_vs_cycle_space() {
        for s in 0..100 {
            let rand_graph = HedgeGraph::<(), ()>::random(6, 9, s);

            let all_cycles = rand_graph.read_tarjan();
            let all_cycles_alt = rand_graph.all_cycles();

            assert_eq!(
                all_cycles.len(),
                all_cycles_alt.len(),
                "{} with seed {s}",
                rand_graph.base_dot()
            );
        }
    }

    // #[test]
    // fn random_graph() {
    //     let rand_graph = NestingGraph::<(), ()>::random(6, 9, 3);

    //     println!(
    //         "{} loop graph: \n {}",
    //         rand_graph.cycle_basis().len(),
    //         rand_graph.base_dot()
    //     );

    //     for c in rand_graph.all_cycles() {
    //         println!(" {}", rand_graph.dot(&c));
    //     }

    //     for c in rand_graph.read_tarjan() {
    //         println!("{}", rand_graph.dot(&c));
    //     }
    // }
    #[allow(non_snake_case)]
    #[test]
    fn K33() {
        let mut builder: HedgeGraphBuidler<(), ()> = HedgeGraphBuidler::new();
        let a = builder.add_node(());
        let b = builder.add_node(());
        let c = builder.add_node(());
        let d = builder.add_node(());
        let e = builder.add_node(());
        let f = builder.add_node(());

        builder.add_edge(a, d, ());
        builder.add_edge(a, e, ());
        builder.add_edge(a, f, ());

        builder.add_edge(b, d, ());
        builder.add_edge(b, e, ());
        builder.add_edge(b, f, ());

        builder.add_edge(c, d, ());
        builder.add_edge(c, e, ());
        builder.add_edge(c, f, ());

        let graph = builder.build();

        println!("{}", graph.dot(&graph.full_node()));

        for (s, _) in graph.all_spinneys() {
            println!("cyclotomatic_number: {}", graph.cyclotomatic_number(&s));
            println!(
                "paton_count_loops {}",
                graph
                    .paton_count_loops(&s, s.filter.iter_ones().next().unwrap())
                    .unwrap()
            );

            println!("paton_cycle_basislen {}", s.cycle_basis(&graph).len());
            println!("{}", graph.dot(&s.to_nesting_node(&graph)));
        }

        assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());
    }

    #[test]
    fn petersen() {
        let mut builder: HedgeGraphBuidler<(), ()> = HedgeGraphBuidler::new();
        let a = builder.add_node(());
        let b = builder.add_node(());
        let c = builder.add_node(());
        let d = builder.add_node(());
        let e = builder.add_node(());
        let f = builder.add_node(());
        let g = builder.add_node(());
        let h = builder.add_node(());
        let i = builder.add_node(());
        let j = builder.add_node(());

        builder.add_edge(a, b, ());
        builder.add_edge(a, f, ());
        builder.add_edge(a, e, ());

        builder.add_edge(b, c, ());
        builder.add_edge(b, g, ());

        builder.add_edge(c, d, ());
        builder.add_edge(c, h, ());

        builder.add_edge(d, e, ());
        builder.add_edge(d, i, ());

        builder.add_edge(e, j, ());

        builder.add_edge(f, h, ());
        builder.add_edge(f, i, ());

        builder.add_edge(g, i, ());
        builder.add_edge(g, j, ());

        builder.add_edge(h, j, ());

        let graph = builder.build();

        println!("{}", graph.base_dot());

        // assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());

        println!("loop count {}", graph.cycle_basis().len());
        println!("cycle count {}", graph.all_cycles().len());
        if let Some((s, v)) = graph
            .all_spinneys()
            .iter()
            .find(|(s, _)| graph.full_filter() == s.filter)
        {
            println!(
                "{}",
                graph.dot(&graph.nesting_node_from_subgraph(s.clone()))
            );
            for (ci, cj) in v {
                println!("{}", graph.dot(ci));
                println!("{}", graph.dot(cj.as_ref().unwrap()));
            }
        } else {
            println!("not found");
        }
    }

    #[test]
    fn wagner_graph() {
        let mut builder: HedgeGraphBuidler<(), ()> = HedgeGraphBuidler::new();
        let n1 = builder.add_node(());
        let n2 = builder.add_node(());
        let n3 = builder.add_node(());
        let n4 = builder.add_node(());
        let n5 = builder.add_node(());
        let n6 = builder.add_node(());
        let n7 = builder.add_node(());
        let n8 = builder.add_node(());

        builder.add_edge(n1, n2, ());
        builder.add_edge(n1, n5, ());

        builder.add_edge(n2, n3, ());
        builder.add_edge(n2, n6, ());

        builder.add_edge(n3, n4, ());
        builder.add_edge(n3, n7, ());

        builder.add_edge(n4, n5, ());
        builder.add_edge(n4, n8, ());

        builder.add_edge(n5, n6, ());

        builder.add_edge(n6, n7, ());

        builder.add_edge(n7, n8, ());

        builder.add_edge(n8, n1, ());

        let graph = builder.build();

        println!("{}", graph.base_dot());

        // assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());

        println!("loop count {}", graph.cycle_basis().len());
        println!("cycle count {}", graph.all_cycles().len());
        if let Some((s, v)) = graph
            .all_spinneys()
            .iter()
            .find(|(s, _)| graph.full_filter() == s.filter)
        {
            println!(
                "{}",
                graph.dot(&graph.nesting_node_from_subgraph(s.clone()))
            );
            for (ci, cj) in v {
                println!("{}", graph.dot(ci));
                println!("{}", graph.dot(cj.as_ref().unwrap()));
            }
        } else {
            println!("not found");
        }
    }

    #[test]
    fn flower_snark() {
        let mut builder: HedgeGraphBuidler<(), ()> = HedgeGraphBuidler::new();
        let n1 = builder.add_node(());
        let n2 = builder.add_node(());
        let n3 = builder.add_node(());
        let n4 = builder.add_node(());
        let n5 = builder.add_node(());
        let n6 = builder.add_node(());
        let n7 = builder.add_node(());
        let n8 = builder.add_node(());
        let n9 = builder.add_node(());
        let n10 = builder.add_node(());
        let n11 = builder.add_node(());
        let n12 = builder.add_node(());
        let n13 = builder.add_node(());
        let n14 = builder.add_node(());
        let n15 = builder.add_node(());
        let n16 = builder.add_node(());
        let n17 = builder.add_node(());
        let n18 = builder.add_node(());
        let n19 = builder.add_node(());
        let n20 = builder.add_node(());

        builder.add_edge(n1, n2, ());
        builder.add_edge(n2, n3, ());
        builder.add_edge(n3, n4, ());
        builder.add_edge(n4, n5, ());
        builder.add_edge(n5, n1, ());

        builder.add_edge(n6, n1, ()); // center
        builder.add_edge(n6, n7, ()); //next

        builder.add_edge(n7, n17, ()); //+10
        builder.add_edge(n7, n8, ()); //next

        builder.add_edge(n8, n13, ()); //+5
        builder.add_edge(n8, n9, ()); //next

        builder.add_edge(n9, n2, ()); //center
        builder.add_edge(n9, n10, ()); //next

        builder.add_edge(n10, n20, ()); //+10
        builder.add_edge(n10, n11, ()); //next

        builder.add_edge(n11, n16, ()); //+5
        builder.add_edge(n11, n12, ()); //next

        builder.add_edge(n12, n3, ()); //center
        builder.add_edge(n12, n13, ()); //next

        builder.add_edge(n13, n14, ()); //next

        builder.add_edge(n14, n19, ()); //+5
        builder.add_edge(n14, n15, ()); //next

        builder.add_edge(n15, n4, ()); //center
        builder.add_edge(n15, n16, ()); //next

        builder.add_edge(n16, n17, ()); //next

        builder.add_edge(n17, n18, ()); //next

        builder.add_edge(n18, n5, ()); //center
        builder.add_edge(n18, n19, ()); //next

        builder.add_edge(n19, n20, ()); //next

        builder.add_edge(n20, n6, ()); //next

        let graph = builder.build();

        println!("{}", graph.base_dot());

        // assert_eq!(graph.all_spinneys().len(), graph.all_spinneys_alt().len());

        println!("loop count {}", graph.cycle_basis().len());
        println!("cycle count {}", graph.all_cycles().len());
        println!(
            "loop count {}",
            graph.paton_count_loops(&graph.full_graph(), 0).unwrap()
        );
        if let Some((s, v)) = graph
            .all_spinneys()
            .iter()
            .find(|(s, _)| graph.full_filter() == s.filter)
        {
            println!(
                "{}",
                graph.dot(&graph.nesting_node_from_subgraph(s.clone()))
            );
            for (ci, cj) in v {
                println!("{}", graph.dot(ci));
                println!("{}", graph.dot(cj.as_ref().unwrap()));
            }
        } else {
            println!("not found");
        }
    }

    use std::time::Instant;

    use super::*;
}
