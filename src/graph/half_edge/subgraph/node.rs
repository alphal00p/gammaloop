use bitvec::slice::BitSlice;
use bitvec::vec::BitVec;
use bitvec::{bitvec, order::Lsb0};
use serde::{Deserialize, Serialize};

use crate::graph::half_edge::{
    GVEdgeAttrs, Hedge, HedgeGraph, HedgeNodeBuilder, InvolutiveMapping,
};

use super::contracted::ContractedSubGraph;
use super::Inclusion;
use super::{internal::InternalSubGraph, SubGraph, SubGraphOps};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct HedgeNode {
    pub internal_graph: InternalSubGraph,
    pub hairs: BitVec,
}

impl Inclusion<Hedge> for HedgeNode {
    fn includes(&self, hedge_id: &Hedge) -> bool {
        self.internal_graph.includes(hedge_id) || self.hairs.includes(hedge_id)
    }

    fn intersects(&self, other: &Hedge) -> bool {
        self.includes(other)
    }
}

impl Inclusion<HedgeNode> for HedgeNode {
    fn includes(&self, other: &HedgeNode) -> bool {
        self.internal_graph.includes(&other.internal_graph)
    }

    fn intersects(&self, other: &HedgeNode) -> bool {
        self.hairs.intersects(&other.hairs)
    }
}

impl Inclusion<BitVec> for HedgeNode {
    fn includes(&self, other: &BitVec) -> bool {
        self.internal_graph.includes(other) || self.hairs.includes(other)
    }

    fn intersects(&self, other: &BitVec) -> bool {
        self.hairs.intersects(other)
    }
}

impl SubGraph for HedgeNode {
    fn nedges<E, V>(&self, graph: &HedgeGraph<E, V>) -> usize {
        self.internal_graph.nedges(graph)
    }

    fn nhedges(&self) -> usize {
        self.hairs.nhedges()
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

        for (n, v) in graph.iter_node_data(&self.internal_graph) {
            if let Some(a) = node_attr(v) {
                out.push_str(
                    format!("  {} [{}];\n", graph.id_from_hairs(n).unwrap().0, a).as_str(),
                );
            }
        }

        for (hedge_id, (edge, incident_node)) in graph.involution.iter() {
            match &edge {
                //Internal graphs never have unpaired edges
                InvolutiveMapping::Identity { data, underlying } => {
                    let attr = if self.internal_graph.filter.includes(&hedge_id) {
                        None
                    } else if self.hairs.includes(&hedge_id) {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray50\"".to_string()),
                            label: None,
                            other: data.as_ref().data.and_then(edge_attr),
                        })
                    } else {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75\"".to_string()),
                            label: None,
                            other: data.as_ref().data.and_then(edge_attr),
                        })
                    };
                    out.push_str(&InvolutiveMapping::<()>::identity_dot(
                        hedge_id,
                        incident_node.0,
                        attr.as_ref(),
                        data.orientation,
                        *underlying,
                    ));
                }
                InvolutiveMapping::Source { data, .. } => {
                    let attr = if self.internal_graph.filter.includes(&hedge_id) {
                        None
                    } else if self.hairs.includes(&hedge_id)
                        && !self.hairs.includes(&graph.involution.inv(hedge_id))
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray50:gray75;0.5\"".to_string()),
                            label: None,
                            other: data.as_ref().data.and_then(edge_attr),
                        })
                    } else if self.hairs.includes(&graph.involution.inv(hedge_id))
                        && !self.hairs.includes(&hedge_id)
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75:gray50;0.5\"".to_string()),
                            label: None,
                            other: data.as_ref().data.and_then(edge_attr),
                        })
                    } else if self.hairs.includes(&graph.involution.inv(hedge_id))
                        && self.hairs.includes(&hedge_id)
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray50\"".to_string()),
                            label: None,
                            other: data.as_ref().data.and_then(edge_attr),
                        })
                    } else {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75\"".to_string()),
                            label: None,
                            other: data.as_ref().data.and_then(edge_attr),
                        })
                    };
                    out.push_str(&InvolutiveMapping::<()>::pair_dot(
                        incident_node.0,
                        graph.involved_node_id(hedge_id).unwrap().0,
                        attr.as_ref(),
                        data.orientation,
                    ));
                }
                InvolutiveMapping::Sink { source_idx } => {
                    if self.internal_graph.filter.includes(&hedge_id)
                        && !self.internal_graph.filter.includes(source_idx)
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

    fn included(&self) -> &BitSlice {
        self.hairs.included()
    }

    fn string_label(&self) -> String {
        (self.hairs.string_label() + "⦻") + self.internal_graph.string_label().as_str()
    }

    fn is_empty(&self) -> bool {
        self.hairs.is_empty() && self.internal_graph.is_empty()
    }

    fn empty(size: usize) -> Self {
        Self {
            internal_graph: InternalSubGraph::empty(size),
            hairs: BitVec::empty(size),
        }
    }
}

impl SubGraphOps for HedgeNode {
    fn complement<E, V>(&self, graph: &HedgeGraph<E, V>) -> Self {
        Self::from_internal_graph(self.internal_graph.complement(graph), graph)
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

    pub fn from_internal_graph<E, V>(subgraph: InternalSubGraph, graph: &HedgeGraph<E, V>) -> Self {
        graph.nesting_node_from_subgraph(subgraph)
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
            let mut bit = externalhedges.get_mut(hedge.0).unwrap();
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
