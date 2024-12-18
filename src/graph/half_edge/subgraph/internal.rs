use bitvec::slice::BitSlice;
use bitvec::vec::BitVec;
use serde::{Deserialize, Serialize};
use std::hash::Hash;
use std::ops::Index;

use crate::graph::half_edge::{GVEdgeAttrs, Hedge, HedgeGraph, InvolutiveMapping, TraversalTree};

use super::{node::HedgeNode, Cycle, Inclusion, SubGraph, SubGraphOps};

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

impl Inclusion<Hedge> for InternalSubGraph {
    fn includes(&self, hedge_id: &Hedge) -> bool {
        self.filter[hedge_id.0]
    }

    fn intersects(&self, other: &Hedge) -> bool {
        self.includes(other)
    }
}

impl Inclusion<InternalSubGraph> for InternalSubGraph {
    fn includes(&self, other: &InternalSubGraph) -> bool {
        self.filter.intersection(&other.filter) == other.filter
    }

    fn intersects(&self, other: &InternalSubGraph) -> bool {
        self.filter.intersection(&other.filter).count_ones() > 0
    }
}

impl Inclusion<BitVec> for InternalSubGraph {
    fn includes(&self, other: &BitVec) -> bool {
        &self.filter.intersection(other) == other
    }

    fn intersects(&self, other: &BitVec) -> bool {
        self.filter.intersection(other).count_ones() > 0
    }
}

impl SubGraph for InternalSubGraph {
    fn nedges<E, V>(&self, _graph: &HedgeGraph<E, V>) -> usize {
        self.nhedges() / 2
    }

    fn nhedges(&self) -> usize {
        self.filter.count_ones()
    }

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
        let mut out = "digraph {\n ".to_string();
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

        for (hedge_id, (edge, incident_node)) in graph.involution.iter() {
            match &edge {
                //Internal graphs never have unpaired edges
                InvolutiveMapping::Identity { data, underlying } => {
                    let attr = if self.filter.includes(&hedge_id) {
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
                        data.orientation,
                        *underlying,
                    ));
                }
                InvolutiveMapping::Source { data, .. } => {
                    let attr = if self.filter.includes(&hedge_id)
                        && !self.filter.includes(&graph.involution.inv(hedge_id))
                        || self.filter.includes(&graph.involution.inv(hedge_id))
                            && !self.filter.includes(&hedge_id)
                    {
                        panic!("Internal subgraphs should never have unpaired edges")
                    } else if !self.filter.includes(&graph.involution.inv(hedge_id))
                        && !self.filter.includes(&hedge_id)
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
                            .get_index_of(graph.involved_node_id(hedge_id).unwrap())
                            .unwrap(),
                        attr.as_ref(),
                        data.orientation,
                    ));
                }
                InvolutiveMapping::Sink { .. } => {}
            }
        }

        out += "}";
        out
    }

    fn hairs(&self, node: &HedgeNode) -> BitVec {
        node.hairs.intersection(&self.filter)
    }

    fn included(&self) -> &BitSlice {
        self.filter.included()
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
        for i in filter.included_iter() {
            if !filter.includes(&graph.involution.inv(i)) {
                return false;
            }
        }
        true
    }

    pub fn add_edge<E, V>(&mut self, hedge: Hedge, graph: &HedgeGraph<E, V>) {
        if !graph.involution.is_identity(hedge) {
            self.filter.set(hedge.0, true);
            self.filter.set(graph.involution.inv(hedge).0, true);
        }
    }

    pub fn remove_edge<E, V>(&mut self, hedge: Hedge, graph: &HedgeGraph<E, V>) {
        if !graph.involution.is_identity(hedge) {
            self.filter.set(hedge.0, false);
            self.filter.set(graph.involution.inv(hedge).0, false);
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
                InvolutiveMapping::Identity { .. } => filter.set(i, false),
                InvolutiveMapping::Sink { source_idx } => {
                    if filter.includes(&Hedge(i)) {
                        filter.set(source_idx.0, true);
                    }
                }
                InvolutiveMapping::Source { sink_idx, .. } => {
                    if filter.includes(&Hedge(i)) {
                        filter.set(sink_idx.0, true);
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
                InvolutiveMapping::Identity { .. } => filter.set(i, false),
                InvolutiveMapping::Sink { source_idx } => {
                    if !filter.includes(&Hedge(i)) {
                        filter.set(source_idx.0, false);
                    }
                }
                InvolutiveMapping::Source { sink_idx, .. } => {
                    if !filter.includes(&Hedge(i)) {
                        filter.set(sink_idx.0, false);
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

    pub fn cycle_basis<E, V>(&self, graph: &HedgeGraph<E, V>) -> (Vec<Cycle>, TraversalTree) {
        let node = graph
            .nodes
            .get_index(self.filter.first_one().unwrap())
            .unwrap()
            .0;
        graph.paton_cycle_basis(self, node, None).unwrap()
    }
}
