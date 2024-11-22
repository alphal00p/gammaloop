use bitvec::vec::BitVec;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::hash::Hash;
use std::ops::Index;

use crate::graph::half_edge::{GVEdgeAttrs, HedgeGraph, InvolutiveMapping, TraversalTree};

use super::{node::HedgeNode, SubGraph, SubGraphOps};

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

        for (hedge_id, (incident_node, edge)) in graph
            .involution
            .hedge_data
            .iter()
            .zip_eq(graph.involution.inv.iter())
            .enumerate()
        {
            match &edge {
                //Internal graphs never have unpaired edges
                InvolutiveMapping::Identity(data) => {
                    let attr = if *self.filter.get(hedge_id).unwrap() {
                        panic!("Internal subgraphs should never have unpaired edges")
                    } else {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75\"".to_string()),
                            label: None,
                            other: data.as_ref().data.and_then(edge_attr),
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
                            other: data.as_ref().data.and_then(edge_attr),
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

    fn includes(&self, i: usize) -> bool {
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
}

impl SubGraphOps for InternalSubGraph {
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
        for (i, m) in graph.involution.inv.iter().enumerate() {
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
        for (i, m) in graph.involution.inv.iter().enumerate() {
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
