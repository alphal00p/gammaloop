use bitvec::vec::BitVec;
use bitvec::{bitvec, order::Lsb0};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::graph::half_edge::{GVEdgeAttrs, HedgeGraph, HedgeNodeBuilder, InvolutiveMapping};

use super::contracted::ContractedSubGraph;
use super::{internal::InternalSubGraph, SubGraph, SubGraphOps};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct HedgeNode {
    pub internal_graph: InternalSubGraph,
    pub hairs: BitVec,
}

impl SubGraph for HedgeNode {
    fn includes(&self, i: usize) -> bool {
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
                    let attr = if *self.internal_graph.filter.get(hedge_id).unwrap() {
                        None
                    } else if *self.hairs.get(hedge_id).unwrap() {
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
                            other: data.as_ref().data.and_then(edge_attr),
                        })
                    } else if *self.hairs.get(graph.involution.inv(hedge_id)).unwrap()
                        && !*self.hairs.get(hedge_id).unwrap()
                    {
                        Some(GVEdgeAttrs {
                            color: Some("\"gray75:gray50;0.5\"".to_string()),
                            label: None,
                            other: data.as_ref().data.and_then(edge_attr),
                        })
                    } else if *self.hairs.get(graph.involution.inv(hedge_id)).unwrap()
                        && *self.hairs.get(hedge_id).unwrap()
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

    fn empty(size: usize) -> Self {
        Self {
            internal_graph: InternalSubGraph::empty(size),
            hairs: BitVec::empty(size),
        }
    }
}

impl SubGraphOps for HedgeNode {
    fn complement<E, V>(&self, graph: &HedgeGraph<E, V>) -> Self {
        let externalhedges = !self.hairs.clone() & !self.internal_graph.filter.clone();

        Self {
            internal_graph: InternalSubGraph::empty(graph.n_hedges()),
            hairs: externalhedges,
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
