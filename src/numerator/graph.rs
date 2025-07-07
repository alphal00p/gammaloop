use linnet::half_edge::{
    involution::{Hedge, Orientation},
    EdgeAccessors, HedgeGraph,
};
use spenso::structure::{
    abstract_index::AbstractIndex, OrderedStructure, PermutedStructure, TensorStructure,
};


use crate::new_graph::{
        parse::{ParseGraph, ParseHedge},
        Edge, Vertex,
    };

pub struct NumGraph {
    graph: HedgeGraph<Edge, Vertex, HedgeData>,
}

#[derive(Clone)]
pub struct HedgeIndices {
    vertex_indices: OrderedStructure,
    edge_indices: OrderedStructure,
}

#[derive(Clone)]
pub struct NumIndices {
    pub color_indices: HedgeIndices,
    pub spin_indices: HedgeIndices,
}

impl NumIndices {
    pub fn parse<'a>(graph: &'a ParseGraph) -> impl FnMut(Hedge, &'a ParseHedge) -> Self {
        |h, _| {
            let eid = graph[&h];
            let flow = graph.flow(h);
            let orientation = graph.orientation(h);

            let flow = match orientation {
                Orientation::Default => flow,
                Orientation::Reversed => -flow,
                Orientation::Undirected => flow,
            };

            let creps = graph[eid].particle.color_reps(flow);
            let sreps = graph[eid].particle.spin_reps();

            let mut init = 0;
            let mut last = None;
            let color_structure: PermutedStructure<_> = creps
                .external_reps_iter()
                .map(|r| {
                    if let Some(l) = last {
                        if l != r {
                            last = Some(r);
                            init = 0;
                        } else {
                            init += 1;
                        }
                    } else {
                        last = Some(r);
                        init = 0;
                    }
                    r.slot(AbstractIndex::Double(h.0 as u16, init))
                })
                .collect();

            let spin_structure: PermutedStructure<_> = sreps
                .external_reps_iter()
                .map(|r| {
                    if let Some(l) = last {
                        if l != r {
                            last = Some(r);
                            init = 0;
                        } else {
                            init += 1;
                        }
                    } else {
                        last = Some(r);
                        init = 0;
                    }
                    r.slot(AbstractIndex::Double(h.0 as u16, init))
                })
                .collect();

            NumIndices {
                color_indices: HedgeIndices::new(color_structure.structure),
                spin_indices: HedgeIndices::new(spin_structure.structure),
            }
        }
    }
}

pub struct HedgeData {
    pub num_indices: NumIndices,
    pub node_order: u8,
}

impl HedgeIndices {
    pub fn new(edge_indices: OrderedStructure) -> Self {
        Self {
            vertex_indices: edge_indices.clone().dual(),
            edge_indices,
        }
    }
}
