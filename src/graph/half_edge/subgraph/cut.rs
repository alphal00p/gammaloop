use ahash::AHashMap;
use bitvec::vec::BitVec;
use by_address::ByAddress;
use indexmap::IndexMap;
use std::hash::Hash;

use crate::graph::half_edge::{
    layout::{LayoutEdge, LayoutParams, LayoutSettings, LayoutVertex},
    EdgeData, HedgeGraph, HedgeGraphBuilder, InvolutiveMapping, Orientation,
};

use super::{SubGraph, SubGraphOps};

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct OrientedCut {
    reference: BitVec,
    cut: BitVec,
}

impl OrientedCut {
    pub fn iter_edges<'a, E, V>(
        &'a self,
        graph: &'a HedgeGraph<E, V>,
    ) -> impl Iterator<Item = (Orientation, EdgeData<&'a E>)> {
        self.cut
            .iter_ones()
            .map(|i| (self.orientation(i, graph), graph.involution.get_data(i)))
    }

    pub fn orientation<E, V>(&self, i: usize, graph: &HedgeGraph<E, V>) -> Orientation {
        let orientation = match graph.involution.get_data(i).orientation {
            Orientation::Undirected => Orientation::Default,
            o => o,
        };

        match (self.cut[i], self.reference[i]) {
            (true, true) => orientation,
            (true, false) => orientation.reverse(),
            _ => Orientation::Undirected,
        }
    }
}

impl SubGraph for OrientedCut {
    fn empty(size: usize) -> Self {
        OrientedCut {
            reference: BitVec::empty(size),
            cut: BitVec::empty(size),
        }
    }

    fn dot<E, V>(
        &self,
        graph: &crate::graph::half_edge::HedgeGraph<E, V>,
        graph_info: String,
        edge_attr: &impl Fn(&E) -> Option<String>,
        node_attr: &impl Fn(&V) -> Option<String>,
    ) -> String {
        let mut out = String::new();

        out
    }

    fn hairs(&self, node: &super::HedgeNode) -> BitVec {
        self.cut.hairs(node)
    }

    fn included(&self) -> impl Iterator<Item = usize> {
        self.cut.included()
    }

    fn includes(&self, i: usize) -> bool {
        self.cut.includes(i)
    }

    fn is_empty(&self) -> bool {
        self.cut.count_ones() == 0
    }

    fn string_label(&self) -> String {
        self.cut.string_label()
    }
}

impl OrientedCut {
    pub fn layout<E, V>(
        self,
        graph: &HedgeGraph<E, V>,
        params: LayoutParams,
        seed: u64,
        temperature: f64,
        edge: f64,
    ) -> HedgeGraph<LayoutEdge<&E>, LayoutVertex<&V>> {
        let mut left = vec![];
        let mut leftright_map = IndexMap::new();
        let mut right = vec![];

        for (j, i) in graph.involution.inv.iter().enumerate() {
            if let InvolutiveMapping::Identity(e) = i {
                let data = ByAddress(e.as_ref().data.unwrap());
                match e.orientation {
                    Orientation::Default => {
                        leftright_map.entry(data).or_insert_with(|| [Some(j), None])[0] = Some(j)
                    }
                    Orientation::Reversed => {
                        leftright_map.entry(data).or_insert_with(|| [None, Some(j)])[1] = Some(j)
                    }
                    _ => {}
                }
            }
        }

        for (_, [i, j]) in leftright_map {
            if let Some(i) = i {
                left.push(i);
            }
            if let Some(j) = j {
                right.push(j);
            }
        }

        let settings =
            LayoutSettings::left_right_square(graph, params, seed, temperature, edge, left, right);
        self.to_owned_graph(graph).layout(settings)
    }

    fn to_owned_graph<E, V>(self, graph: &HedgeGraph<E, V>) -> HedgeGraph<&E, &V> {
        let mut builder = HedgeGraphBuilder::new();

        let mut nodeidmap = AHashMap::new();
        for (n, v) in &graph.nodes {
            nodeidmap.insert(n, builder.add_node(v));
        }

        let complement = self.cut.complement(graph);
        for i in complement.included() {
            let source = graph.get_incident_node_id(i);
            match &graph.involution.inv[i] {
                InvolutiveMapping::Identity(e) => builder.add_external_edge(
                    nodeidmap[source],
                    e.as_ref().data.unwrap(),
                    e.orientation,
                ),
                InvolutiveMapping::Source((e, h)) => {
                    let sink = graph.get_incident_node_id(*h);
                    if complement.includes(*h) {
                    } else {
                        let orientation = match e.orientation {
                            Orientation::Undirected => Orientation::Default,
                            o => o,
                        };
                        builder.add_external_edge(
                            nodeidmap[source],
                            e.as_ref().data.unwrap(),
                            orientation,
                        );
                        builder.add_external_edge(
                            nodeidmap[sink],
                            e.as_ref().data.unwrap(),
                            orientation.reverse(),
                        );
                    }
                }
                InvolutiveMapping::Sink(h) => {
                    let sink = graph.get_incident_node_id(*h);
                    let e = graph.involution.get_data(*h);
                    if !complement.includes(*h) {
                        let orientation = match e.orientation {
                            Orientation::Undirected => Orientation::Reversed,
                            o => o.reverse(),
                        };
                        builder.add_external_edge(
                            nodeidmap[source],
                            e.as_ref().data.unwrap(),
                            orientation,
                        );
                        builder.add_external_edge(
                            nodeidmap[sink],
                            e.as_ref().data.unwrap(),
                            orientation.reverse(),
                        );
                    }
                }
            }
        }

        builder.build()
    }
}
