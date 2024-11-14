use std::{
    hash::Hash,
    ops::{BitAndAssign, BitOrAssign, BitXorAssign, Index},
};

use bitvec::{bitvec, order::Lsb0, vec::BitVec};
use serde::{Deserialize, Serialize};

use super::{GVEdgeAttrs, HedgeGraph, HedgeNodeBuilder, InvolutiveMapping, TraversalTree};

const BASE62_ALPHABET: &[u8; 62] =
    b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

pub trait SubGraph: Clone {
    fn string_label(&self) -> String;
    fn included(&self) -> impl Iterator<Item = usize>;
    fn is_included(&self, i: usize) -> bool;

    fn dot<E, V>(
        &self,
        graph: &HedgeGraph<E, V>,
        graph_info: String,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> String;

    fn hairs(&self, node: &HedgeNode) -> BitVec;
    fn empty(size: usize) -> Self;
    fn is_empty(&self) -> bool;
    fn intersect_with(&mut self, other: &Self);
    fn union_with(&mut self, other: &Self);
    fn sym_diff_with(&mut self, other: &Self);
    fn empty_intersection(&self, other: &Self) -> bool;
    fn empty_union(&self, other: &Self) -> bool;
    fn intersection(&self, other: &Self) -> Self {
        let mut new = self.clone();
        new.intersect_with(other);
        new
    }
    fn union(&self, other: &Self) -> Self {
        let mut new = self.clone();
        new.union_with(other);

        new
    }
    fn sym_diff(&self, other: &Self) -> Self {
        let mut new = self.clone();
        new.sym_diff_with(other);
        new
    }

    fn subtract_with(&mut self, other: &Self);
    fn subtract(&self, other: &Self) -> Self {
        let mut new = self.clone();
        new.subtract_with(other);
        new
    }

    fn complement<E, V>(&self, graph: &HedgeGraph<E, V>) -> Self;
}

impl SubGraph for BitVec {
    fn is_included(&self, i: usize) -> bool {
        self[i]
    }

    fn dot<E, V>(
        &self,
        graph: &HedgeGraph<E, V>,
        graph_info: String,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> String {
        let mut out = "graph {\n ".to_string();
        out.push_str(
            "  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\"; layout=\"neato\";\n ",
        );

        out.push_str(graph_info.as_str());

        for (n, v) in graph.iter_node_data(self) {
            out.push_str(
                format!(
                    "  {} [{}];\n",
                    graph.nodes.get_index_of(n).unwrap(),
                    node_attr(v).map_or("".into(), |x| x).as_str()
                )
                .as_str(),
            );
        }

        for (hedge_id, (incident_node, edge)) in graph.involution.inv.iter().enumerate() {
            match &edge {
                //Internal graphs never have unpaired edges
                InvolutiveMapping::Identity(data) => {
                    let attr = if *self.get(hedge_id).unwrap() {
                        Some(GVEdgeAttrs {
                            color: Some("\"green\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    };
                    out.push_str(&InvolutiveMapping::<()>::identity_dot(
                        hedge_id,
                        graph.nodes.get_index_of(incident_node).unwrap(),
                        attr.as_ref(),
                    ));
                }
                InvolutiveMapping::Source((data, _)) => {
                    let attr = if *self.get(hedge_id).unwrap()
                        && !*self.get(graph.involution.inv(hedge_id)).unwrap()
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75:blue;0.5\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else if *self.get(graph.involution.inv(hedge_id)).unwrap()
                        && !*self.get(hedge_id).unwrap()
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"red:gray75;0.5\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else if *self.get(graph.involution.inv(hedge_id)).unwrap()
                        && *self.get(hedge_id).unwrap()
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else {
                        None
                    };
                    out.push_str(&InvolutiveMapping::<()>::pair_dot(
                        graph.nodes.get_index_of(incident_node).unwrap(),
                        graph
                            .nodes
                            .get_index_of(graph.involution.get_connected_node_id(hedge_id).unwrap())
                            .unwrap(),
                        attr.as_ref(),
                    ));
                }
                InvolutiveMapping::Sink(_) => {}
            }
        }

        out += "}";
        out
    }

    fn hairs(&self, node: &HedgeNode) -> BitVec {
        node.hairs.intersection(self)
    }

    fn included(&self) -> impl Iterator<Item = usize> {
        self.iter_ones()
    }
    fn empty(size: usize) -> Self {
        bitvec![usize, Lsb0; 0; size]
    }

    fn string_label(&self) -> String {
        if self.is_empty() {
            return "0".to_string();
        }

        let mut digits = vec![0u8]; // Initialize with a single zero digit

        // Iterate over the bits from MSB to LSB
        for bit in self.iter().by_vals().rev() {
            let mut carry = 0u8;

            // Multiply existing digits by 2 (shift left)
            for digit in &mut digits {
                let temp = (*digit as u16) * 2 + carry as u16;
                *digit = (temp % 62) as u8;
                carry = (temp / 62) as u8;
            }

            if carry > 0 {
                digits.push(carry);
            }

            // Add the current bit (if it's 1)
            if bit {
                let mut carry = 1u8;
                for digit in &mut digits {
                    let temp = *digit as u16 + carry as u16;
                    *digit = (temp % 62) as u8;
                    carry = (temp / 62) as u8;

                    if carry == 0 {
                        break;
                    }
                }
                if carry > 0 {
                    digits.push(carry);
                }
            }
        }

        // Map digits to base62 characters and reverse the result
        let base62_string: String = digits
            .iter()
            .rev()
            .map(|&d| BASE62_ALPHABET[d as usize] as char)
            .collect();

        base62_string
    }

    fn is_empty(&self) -> bool {
        self.count_ones() == 0
    }

    fn intersect_with(&mut self, other: &Self) {
        self.bitand_assign(other)
    }

    fn union_with(&mut self, other: &Self) {
        self.bitor_assign(other)
    }

    fn sym_diff_with(&mut self, other: &Self) {
        self.bitxor_assign(other)
    }

    fn empty_union(&self, other: &Self) -> bool {
        self.union(other).count_ones() == 0
    }

    fn empty_intersection(&self, other: &Self) -> bool {
        self.intersection(other).count_ones() == 0
    }

    fn complement<E, V>(&self, _graph: &HedgeGraph<E, V>) -> Self {
        !self.clone()
    }

    fn subtract_with(&mut self, other: &Self) {
        self.bitand_assign(!other.clone());
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Eq)]
pub struct InternalSubGraph {
    // cannot be hairy. I.e. it must always have paired hedges.
    // To represent a hairy subgraph, use a ContractedSubGraph
    pub filter: BitVec,
    pub loopcount: Option<usize>,
}

impl Hash for InternalSubGraph {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.filter.hash(state);
    }
}

impl PartialEq for InternalSubGraph {
    fn eq(&self, other: &Self) -> bool {
        self.filter == other.filter
    }
}

impl Index<usize> for InternalSubGraph {
    type Output = bool;

    fn index(&self, index: usize) -> &Self::Output {
        self.filter.index(index)
    }
}

impl PartialOrd for InternalSubGraph {
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

impl SubGraph for InternalSubGraph {
    fn empty(size: usize) -> Self {
        InternalSubGraph {
            filter: BitVec::empty(size),
            loopcount: Some(0),
        }
    }
    fn dot<E, V>(
        &self,
        graph: &HedgeGraph<E, V>,
        graph_info: String,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> String {
        let mut out = "graph {\n ".to_string();
        out.push_str(
            "  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\"; layout=\"neato\";\n ",
        );

        out.push_str(graph_info.as_str());

        for (n, v) in graph.iter_node_data(self) {
            out.push_str(
                format!(
                    "  {} [{}];\n",
                    graph.nodes.get_index_of(n).unwrap(),
                    node_attr(v).map_or("".into(), |x| x).as_str()
                )
                .as_str(),
            );
        }

        for (hedge_id, (incident_node, edge)) in graph.involution.inv.iter().enumerate() {
            match &edge {
                //Internal graphs never have unpaired edges
                InvolutiveMapping::Identity(data) => {
                    let attr = if *self.filter.get(hedge_id).unwrap() {
                        panic!("Internal subgraphs should never have unpaired edges")
                    } else {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    };
                    out.push_str(&InvolutiveMapping::<()>::identity_dot(
                        hedge_id,
                        graph.nodes.get_index_of(incident_node).unwrap(),
                        attr.as_ref(),
                    ));
                }
                InvolutiveMapping::Source((data, _)) => {
                    let attr = if *self.filter.get(hedge_id).unwrap()
                        && !*self.filter.get(graph.involution.inv(hedge_id)).unwrap()
                        || *self.filter.get(graph.involution.inv(hedge_id)).unwrap()
                            && !*self.filter.get(hedge_id).unwrap()
                    {
                        panic!("Internal subgraphs should never have unpaired edges")
                    } else if !*self.filter.get(graph.involution.inv(hedge_id)).unwrap()
                        && !*self.filter.get(hedge_id).unwrap()
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else {
                        None
                    };
                    out.push_str(&InvolutiveMapping::<()>::pair_dot(
                        graph.nodes.get_index_of(incident_node).unwrap(),
                        graph
                            .nodes
                            .get_index_of(graph.involution.get_connected_node_id(hedge_id).unwrap())
                            .unwrap(),
                        attr.as_ref(),
                    ));
                }
                InvolutiveMapping::Sink(_) => {}
            }
        }

        out += "}";
        out
    }

    fn hairs(&self, node: &HedgeNode) -> BitVec {
        node.hairs.intersection(&self.filter)
    }

    fn is_included(&self, i: usize) -> bool {
        self.filter[i]
    }
    fn included(&self) -> impl Iterator<Item = usize> {
        self.filter.iter_ones()
    }

    fn string_label(&self) -> String {
        self.filter.string_label()
    }
    fn is_empty(&self) -> bool {
        self.filter.count_ones() == 0
    }

    fn intersect_with(&mut self, other: &InternalSubGraph) {
        self.filter &= &other.filter;
        self.loopcount = None;
    }

    fn union_with(&mut self, other: &InternalSubGraph) {
        self.filter |= &other.filter;
        self.loopcount = None;
    }

    fn sym_diff_with(&mut self, other: &InternalSubGraph) {
        self.filter ^= &other.filter;
        self.loopcount = None;
    }

    fn empty_intersection(&self, other: &InternalSubGraph) -> bool {
        (self.filter.clone() & &other.filter).count_ones() == 0
    }

    fn empty_union(&self, other: &InternalSubGraph) -> bool {
        (self.filter.clone() | &other.filter).count_ones() == 0
    }

    fn complement<E, V>(&self, graph: &HedgeGraph<E, V>) -> Self {
        InternalSubGraph {
            filter: !self.filter.clone() & !graph.external_filter(),
            loopcount: None,
        }
    }

    fn subtract_with(&mut self, other: &Self) {
        self.filter = !other.filter.clone() & &self.filter;
        self.loopcount = None;
    }
}

impl InternalSubGraph {
    fn valid_filter<E, V>(filter: &BitVec, graph: &HedgeGraph<E, V>) -> bool {
        for i in filter.iter_ones() {
            if !filter[graph.involution.inv(i)] {
                return false;
            }
        }
        true
    }

    pub fn add_edge<E, V>(&mut self, hedge: usize, graph: &HedgeGraph<E, V>) {
        if !graph.involution.is_identity(hedge) {
            self.filter.set(hedge, true);
            self.filter.set(graph.involution.inv(hedge), true);
        }
    }

    pub fn remove_edge<E, V>(&mut self, hedge: usize, graph: &HedgeGraph<E, V>) {
        if !graph.involution.is_identity(hedge) {
            self.filter.set(hedge, false);
            self.filter.set(graph.involution.inv(hedge), false);
        }
    }

    pub fn try_new<E, V>(filter: BitVec, graph: &HedgeGraph<E, V>) -> Option<Self> {
        if filter.len() != graph.involution.len() {
            return None;
        }
        if !Self::valid_filter(&filter, graph) {
            return None;
        }

        Some(InternalSubGraph {
            filter,
            loopcount: None,
        })
    }

    pub fn cleaned_filter_optimist<E, V>(mut filter: BitVec, graph: &HedgeGraph<E, V>) -> Self {
        for (i, (_, m)) in graph.involution.inv.iter().enumerate() {
            match m {
                InvolutiveMapping::Identity(_) => filter.set(i, false),
                InvolutiveMapping::Sink(j) => {
                    if filter[i] {
                        filter.set(*j, true);
                    }
                }
                InvolutiveMapping::Source((_, j)) => {
                    if filter[i] {
                        filter.set(*j, true);
                    }
                }
            }
        }
        InternalSubGraph {
            filter,
            loopcount: None,
        }
    }

    pub fn cleaned_filter_pessimist<E, V>(mut filter: BitVec, graph: &HedgeGraph<E, V>) -> Self {
        for (i, (_, m)) in graph.involution.inv.iter().enumerate() {
            match m {
                InvolutiveMapping::Identity(_) => filter.set(i, false),
                InvolutiveMapping::Sink(j) => {
                    if !filter[i] {
                        filter.set(*j, false);
                    }
                }
                InvolutiveMapping::Source((_, j)) => {
                    if !filter[i] {
                        filter.set(*j, false);
                    }
                }
            }
        }

        InternalSubGraph {
            filter,
            loopcount: None,
        }
    }

    pub fn valid<E, V>(&self, graph: &HedgeGraph<E, V>) -> bool {
        Self::valid_filter(&self.filter, graph)
    }

    pub fn to_hairy_subgraph<E, V>(&self, graph: &HedgeGraph<E, V>) -> HedgeNode {
        graph.nesting_node_from_subgraph(self.clone())
    }

    pub fn set_loopcount<E, V>(&mut self, graph: &HedgeGraph<E, V>) {
        self.loopcount = Some(graph.cyclotomatic_number(self));
    }

    pub fn cycle_basis<E, V>(
        &self,
        graph: &HedgeGraph<E, V>,
    ) -> (Vec<InternalSubGraph>, TraversalTree) {
        let node = graph
            .nodes
            .get_index(self.filter.first_one().unwrap())
            .unwrap()
            .0;
        graph.paton_cycle_basis(self, node).unwrap()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct ContractedSubGraph {
    pub internal_graph: InternalSubGraph, // cannot have any external hedges (i.e. unpaired hedges)
    pub allhedges: BitVec,                // all hedges , including that are in the internal graph.
}

impl SubGraph for ContractedSubGraph {
    fn is_included(&self, i: usize) -> bool {
        self.allhedges[i] && !self.internal_graph.filter[i]
    }

    fn dot<E, V>(
        &self,
        graph: &HedgeGraph<E, V>,
        graph_info: String,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> String {
        let mut out = "graph {\n ".to_string();
        out.push_str(
            "  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\"; layout=\"neato\";\n ",
        );

        out.push_str(graph_info.as_str());

        for (n, v) in graph.iter_node_data(self) {
            out.push_str(
                format!(
                    "  {} [{}];\n",
                    graph.nodes.get_index_of(n).unwrap(),
                    node_attr(v).map_or("".into(), |x| x).as_str()
                )
                .as_str(),
            );
        }

        for (hedge_id, (incident_node, edge)) in graph.involution.inv.iter().enumerate() {
            match &edge {
                //External edges cannot be in the contracted subgraph
                InvolutiveMapping::Identity(data) => {
                    let attr = if *self.allhedges.get(hedge_id).unwrap() {
                        Some(GVEdgeAttrs {
                            color: Some("\"green\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    };
                    out.push_str(&InvolutiveMapping::<()>::identity_dot(
                        hedge_id,
                        graph.nodes.get_index_of(incident_node).unwrap(),
                        attr.as_ref(),
                    ));
                }
                InvolutiveMapping::Source((data, sink)) => {
                    let attr = if self.internal_graph.is_included(hedge_id) {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray5\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else if self.allhedges.is_included(hedge_id)
                        && !self.allhedges.is_included(*sink)
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"red:gray75;0.5\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else if !self.allhedges.is_included(hedge_id)
                        && self.allhedges.is_included(*sink)
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75:blue;0.5\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else if !self.allhedges.is_included(hedge_id)
                        && !self.allhedges.is_included(*sink)
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else {
                        None
                    };
                    out.push_str(&InvolutiveMapping::<()>::pair_dot(
                        graph.nodes.get_index_of(incident_node).unwrap(),
                        graph
                            .nodes
                            .get_index_of(graph.involution.get_connected_node_id(hedge_id).unwrap())
                            .unwrap(),
                        attr.as_ref(),
                    ));
                }
                InvolutiveMapping::Sink(_) => {}
            }
        }

        out += "}";
        out
    }

    fn hairs(&self, node: &HedgeNode) -> BitVec {
        let mut hairs = self.allhedges.intersection(&node.hairs);
        hairs.subtract_with(&self.internal_graph.filter);
        hairs
    }

    fn included(&self) -> impl Iterator<Item = usize> {
        self.allhedges
            .iter_ones()
            .filter(move |&i| self.is_included(i))
    }

    fn string_label(&self) -> String {
        self.allhedges.string_label() + "⊛" + self.internal_graph.string_label().as_str()
    }

    fn is_empty(&self) -> bool {
        self.allhedges.is_empty()
    }

    fn complement<E, V>(&self, graph: &HedgeGraph<E, V>) -> Self {
        let externalhedges = !self.allhedges.clone() & !self.internal_graph.filter.clone();

        Self {
            internal_graph: InternalSubGraph::empty(graph.n_hedges()),
            allhedges: externalhedges,
        }
    }

    fn empty(size: usize) -> Self {
        Self {
            internal_graph: InternalSubGraph::empty(size),
            allhedges: BitVec::empty(size),
        }
    }

    fn union_with(&mut self, other: &Self) {
        // union is the intersection of the internal graphs, and the union of the external graph.
        self.internal_graph.intersect_with(&other.internal_graph);
        self.allhedges.union_with(&other.allhedges);
    }

    fn intersect_with(&mut self, other: &Self) {
        // intersection is the union of the internal graphs, and the intersection of the external graph.
        self.internal_graph.union_with(&other.internal_graph);
        self.allhedges.intersect_with(&other.allhedges);
    }

    fn sym_diff_with(&mut self, other: &Self) {
        // external hedges that are only present in one of the two graphs.
        // contracted parts unioned
        self.internal_graph.union_with(&other.internal_graph);
        self.allhedges.sym_diff_with(&other.allhedges);
    }

    fn empty_intersection(&self, other: &Self) -> bool {
        self.allhedges.empty_intersection(&other.allhedges)
    }

    fn empty_union(&self, other: &Self) -> bool {
        self.allhedges.empty_union(&other.allhedges)
    }

    fn subtract_with(&mut self, other: &Self) {
        self.internal_graph.union_with(&other.internal_graph);
        self.allhedges.subtract_with(&other.allhedges);
    }
}

impl ContractedSubGraph {
    pub fn all_edges(&self) -> BitVec {
        self.internal_graph.filter.clone() | &self.allhedges
    }

    pub fn weakly_disjoint(&self, other: &ContractedSubGraph) -> bool {
        let internals = self.internal_graph.filter.clone() & &other.internal_graph.filter;

        internals.count_ones() == 0
    }

    pub fn strongly_disjoint(&self, other: &ContractedSubGraph) -> bool {
        let internals = self.internal_graph.filter.clone() & &other.internal_graph.filter;

        let externals_in_self = self.internal_graph.filter.clone() & &other.allhedges;
        let externals_in_other = self.allhedges.clone() & &other.internal_graph.filter;

        internals.count_ones() == 0
            && externals_in_self.count_ones() == 0
            && externals_in_other.count_ones() == 0
    }

    pub fn node_from_pos(pos: &[usize], len: usize) -> ContractedSubGraph {
        ContractedSubGraph {
            allhedges: ContractedSubGraph::filter_from_pos(pos, len),
            internal_graph: InternalSubGraph::empty(len),
        }
    }

    pub fn filter_from_pos(pos: &[usize], len: usize) -> BitVec {
        let mut filter = bitvec![usize, Lsb0; 0; len];

        for &i in pos {
            filter.set(i, true);
        }

        filter
    }

    pub fn internal_graph_union(&self, other: &ContractedSubGraph) -> InternalSubGraph {
        InternalSubGraph {
            filter: self.internal_graph.filter.clone() | &other.internal_graph.filter,
            loopcount: None,
        }
    }

    pub fn from_builder<V>(builder: &HedgeNodeBuilder<V>, len: usize) -> Self {
        let internal_graph = InternalSubGraph::empty(len);
        let mut externalhedges = bitvec![usize, Lsb0; 0; len];

        for hedge in &builder.hedges {
            let mut bit = externalhedges.get_mut(*hedge).unwrap();
            *bit = true;
        }

        ContractedSubGraph {
            internal_graph,
            allhedges: externalhedges,
        }
    }

    pub fn is_node(&self) -> bool {
        self.internal_graph.is_empty()
    }

    pub fn is_subgraph(&self) -> bool {
        !self.is_node()
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct HedgeNode {
    pub internal_graph: InternalSubGraph,
    pub hairs: BitVec,
}

impl SubGraph for HedgeNode {
    fn is_included(&self, i: usize) -> bool {
        self.hairs[i] || !self.internal_graph.filter[i]
    }

    fn dot<E, V>(
        &self,
        graph: &HedgeGraph<E, V>,
        graph_info: String,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> String {
        let mut out = "graph {\n ".to_string();
        out.push_str(
            "  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\"; layout=\"neato\";\n ",
        );

        out.push_str(graph_info.as_str());

        for (n, v) in graph.iter_node_data(&self.internal_graph) {
            out.push_str(
                format!(
                    "  {} [{}];\n",
                    graph.nodes.get_index_of(n).unwrap(),
                    node_attr(v).map_or("".into(), |x| x).as_str()
                )
                .as_str(),
            );
        }

        for (hedge_id, (incident_node, edge)) in graph.involution.inv.iter().enumerate() {
            match &edge {
                //Internal graphs never have unpaired edges
                InvolutiveMapping::Identity(data) => {
                    let attr = if *self.internal_graph.filter.get(hedge_id).unwrap() {
                        None
                    } else if *self.hairs.get(hedge_id).unwrap() {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray50\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    };
                    out.push_str(&InvolutiveMapping::<()>::identity_dot(
                        hedge_id,
                        graph.nodes.get_index_of(incident_node).unwrap(),
                        attr.as_ref(),
                    ));
                }
                InvolutiveMapping::Source((data, _)) => {
                    let attr = if *self.internal_graph.filter.get(hedge_id).unwrap() {
                        None
                    } else if *self.hairs.get(hedge_id).unwrap()
                        && !*self.hairs.get(graph.involution.inv(hedge_id)).unwrap()
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray50:gray75;0.5\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else if *self.hairs.get(graph.involution.inv(hedge_id)).unwrap()
                        && !*self.hairs.get(hedge_id).unwrap()
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75:gray50;0.5\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else if *self.hairs.get(graph.involution.inv(hedge_id)).unwrap()
                        && *self.hairs.get(hedge_id).unwrap()
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray50\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    } else {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75\"".to_string()),
                            label: None,
                            other: data.as_ref().and_then(edge_attr),
                        })
                    };
                    out.push_str(&InvolutiveMapping::<()>::pair_dot(
                        graph.nodes.get_index_of(incident_node).unwrap(),
                        graph
                            .nodes
                            .get_index_of(graph.involution.get_connected_node_id(hedge_id).unwrap())
                            .unwrap(),
                        attr.as_ref(),
                    ));
                }
                InvolutiveMapping::Sink(i) => {
                    if *self.internal_graph.filter.get(hedge_id).unwrap()
                        && !*self.internal_graph.filter.get(*i).unwrap()
                    {
                        panic!("Internal graph has a dangling sink edge")
                    }
                }
            }
        }

        out += "}";
        out
    }

    fn hairs(&self, node: &HedgeNode) -> BitVec {
        let mut hairs = node.all_edges();
        hairs.intersect_with(&node.hairs);
        hairs
    }

    fn included(&self) -> impl Iterator<Item = usize> {
        self.hairs.iter_ones().chain(self.internal_graph.included())
    }

    fn string_label(&self) -> String {
        (self.hairs.string_label() + "⦻") + self.internal_graph.string_label().as_str()
    }

    fn is_empty(&self) -> bool {
        self.hairs.is_empty() && self.internal_graph.is_empty()
    }

    fn complement<E, V>(&self, graph: &HedgeGraph<E, V>) -> Self {
        let externalhedges = !self.hairs.clone() & !self.internal_graph.filter.clone();

        Self {
            internal_graph: InternalSubGraph::empty(graph.n_hedges()),
            hairs: externalhedges,
        }
    }

    fn empty(size: usize) -> Self {
        Self {
            internal_graph: InternalSubGraph::empty(size),
            hairs: BitVec::empty(size),
        }
    }

    fn union_with(&mut self, other: &Self) {
        // union is the intersection of the internal graphs, and the union of the external graph.
        self.internal_graph.intersect_with(&other.internal_graph);
        self.hairs.union_with(&other.hairs);
    }

    fn intersect_with(&mut self, other: &Self) {
        // intersection is the union of the internal graphs, and the intersection of the external graph.
        self.internal_graph.union_with(&other.internal_graph);
        self.hairs.intersect_with(&other.hairs);
    }

    fn sym_diff_with(&mut self, other: &Self) {
        // external hedges that are only present in one of the two graphs.
        // contracted parts unioned
        self.internal_graph.union_with(&other.internal_graph);
        self.hairs.sym_diff_with(&other.hairs);
    }

    fn empty_intersection(&self, other: &Self) -> bool {
        self.hairs.empty_intersection(&other.hairs)
    }

    fn empty_union(&self, other: &Self) -> bool {
        self.hairs.empty_union(&other.hairs)
    }

    fn subtract_with(&mut self, other: &Self) {
        self.internal_graph.union_with(&other.internal_graph);
        self.hairs.subtract_with(&other.hairs);
    }
}
impl HedgeNode {
    fn all_edges(&self) -> BitVec {
        self.internal_graph.filter.clone() | &self.hairs
    }

    pub fn weakly_disjoint(&self, other: &HedgeNode) -> bool {
        let internals = self.internal_graph.filter.clone() & &other.internal_graph.filter;

        internals.count_ones() == 0
    }

    pub fn strongly_disjoint(&self, other: &HedgeNode) -> bool {
        let internals = self.internal_graph.filter.clone() & &other.internal_graph.filter;

        let externals_in_self = self.internal_graph.filter.clone() & &other.hairs;
        let externals_in_other = self.hairs.clone() & &other.internal_graph.filter;

        internals.count_ones() == 0
            && externals_in_self.count_ones() == 0
            && externals_in_other.count_ones() == 0
    }

    pub fn node_from_pos(pos: &[usize], len: usize) -> HedgeNode {
        HedgeNode {
            hairs: HedgeNode::filter_from_pos(pos, len),
            internal_graph: InternalSubGraph::empty(len),
        }
    }

    pub fn filter_from_pos(pos: &[usize], len: usize) -> BitVec {
        let mut filter = bitvec![usize, Lsb0; 0; len];

        for &i in pos {
            filter.set(i, true);
        }

        filter
    }

    pub fn internal_graph_union(&self, other: &HedgeNode) -> InternalSubGraph {
        InternalSubGraph {
            filter: self.internal_graph.filter.clone() | &other.internal_graph.filter,
            loopcount: None,
        }
    }

    pub fn from_builder<V>(builder: &HedgeNodeBuilder<V>, len: usize) -> Self {
        let internal_graph = InternalSubGraph::empty(len);
        let mut externalhedges = bitvec![usize, Lsb0; 0; len];

        for hedge in &builder.hedges {
            let mut bit = externalhedges.get_mut(*hedge).unwrap();
            *bit = true;
        }

        HedgeNode {
            internal_graph,
            hairs: externalhedges,
        }
    }

    pub fn is_node(&self) -> bool {
        self.internal_graph.is_empty()
    }

    pub fn valid<E, V>(&self, _graph: &HedgeGraph<E, V>) -> bool {
        true
    }

    pub fn is_subgraph(&self) -> bool {
        !self.is_node()
    }
}

impl From<HedgeNode> for ContractedSubGraph {
    fn from(value: HedgeNode) -> Self {
        ContractedSubGraph {
            internal_graph: value.internal_graph,
            allhedges: value.hairs,
        }
    }
}
