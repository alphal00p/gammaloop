use ahash::AHashMap;
use itertools::Itertools;
use linnet::half_edge::{
    builder::HedgeGraphBuilder,
    involution::{HedgePair, Orientation},
    nodestore::NodeStorageOps,
    HedgeGraph, HedgeGraphError, NoData, NodeIndex,
};
use symbolica::graph::Graph as SymbolicaGraph;

pub trait HedgeGraphExt<N, E> {
    type Error;
    fn from_sym(graph: SymbolicaGraph<N, E>) -> Self;

    fn to_sym(graph: &Self) -> Result<SymbolicaGraph<&N, &E>, Self::Error>;
}

impl<N: Clone, E: Clone, S: NodeStorageOps<NodeData = N>> HedgeGraphExt<N, E>
    for HedgeGraph<E, N, NoData, S>
{
    fn from_sym(graph: SymbolicaGraph<N, E>) -> Self {
        let mut builder = HedgeGraphBuilder::new();
        let mut map = AHashMap::new();

        for (i, node) in graph.nodes().iter().enumerate() {
            map.insert(i, builder.add_node(node.data.clone()));
        }

        // let mut edges = graph.edges().to_vec();
        // edges.sort_by(|a, b| a.vertices.cmp(&b.vertices));

        for edge in graph
            .edges()
            .iter()
            .sorted_by(|a, b| a.vertices.cmp(&b.vertices))
        {
            let vertices = edge.vertices;
            let source = map[&vertices.0];
            let sink = map[&vertices.1];
            builder.add_edge(source, sink, edge.data.clone(), edge.directed);
        }

        builder.into()
    }

    type Error = HedgeGraphError;

    fn to_sym(value: &HedgeGraph<E, N, NoData, S>) -> Result<SymbolicaGraph<&N, &E>, Self::Error> {
        let mut graph = SymbolicaGraph::new();
        let mut map = AHashMap::new();

        for (n, (_, _, node)) in value.iter_nodes().enumerate() {
            map.insert(NodeIndex(n), graph.add_node(node));
        }

        for (i, _, d) in value.iter_edges() {
            if let HedgePair::Paired { source, sink } = i {
                let source = map[&value.node_id(source)];
                let sink = map[&value.node_id(sink)];

                let data = d.data;
                let orientation = d.orientation;

                match orientation {
                    Orientation::Default => {
                        graph
                            .add_edge(source, sink, true, data)
                            .map_err(HedgeGraphError::SymbolicaError)?;
                    }
                    Orientation::Reversed => {
                        graph
                            .add_edge(sink, source, true, data)
                            .map_err(HedgeGraphError::SymbolicaError)?;
                    }
                    Orientation::Undirected => {
                        graph
                            .add_edge(source, sink, false, data)
                            .map_err(HedgeGraphError::SymbolicaError)?;
                    }
                }
            } else {
                return Err(HedgeGraphError::HasIdentityHedge);
            }
        }

        Ok(graph)
    }
}
