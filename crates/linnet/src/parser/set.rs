use std::path::Path;

#[cfg(feature = "serde")]
use figment;
use figment::Figment;

use crate::half_edge::{
    nodestore::{DefaultNodeStore, NodeStorage, NodeStorageOps},
    HedgeGraph,
};

use super::{
    subgraph_free::SubGraphFreeGraph, DotEdgeData, DotGraph, DotHedgeData, DotVertexData,
    GlobalData, HedgeParseError,
};

pub struct GraphSet<E, V, H, G, S: NodeStorage<NodeData = V>> {
    pub global_data: Vec<G>,
    pub set: Vec<HedgeGraph<E, V, H, S>>,
}

pub type DotGraphSet = GraphSet<
    DotEdgeData,
    DotVertexData,
    DotHedgeData,
    GlobalData,
    DefaultNodeStore<DotVertexData>,
>;

impl<S: NodeStorageOps<NodeData = DotVertexData>>
    GraphSet<DotEdgeData, DotVertexData, DotHedgeData, GlobalData, S>
{
    #[allow(clippy::result_large_err, clippy::type_complexity)]
    pub fn from_file<'a, P>(p: P) -> Result<Self, HedgeParseError<'a, (), (), (), ()>>
    where
        P: AsRef<Path>,
    {
        let ast_graphs = dot_parser::ast::Graphs::from_file(p)?;
        let mut set = Vec::with_capacity(ast_graphs.graphs.len());
        let mut global_data = Vec::new();
        for g in ast_graphs.graphs {
            let can_graph = SubGraphFreeGraph::from(g);
            let graph = DotGraph::from((can_graph, Figment::new()));

            set.push(graph.graph);
            global_data.push(graph.global_data);
        }
        Ok(GraphSet { set, global_data })
    }

    #[allow(clippy::result_large_err, clippy::type_complexity)]
    pub fn from_string<'a, Str: AsRef<str>>(
        s: Str,
    ) -> Result<Self, HedgeParseError<'a, (), (), (), ()>> {
        let ast_graphs = dot_parser::ast::Graphs::try_from(s.as_ref())?;

        let mut global_data = Vec::new();
        let mut set = Vec::with_capacity(ast_graphs.graphs.len());
        for g in ast_graphs.graphs {
            let can_graph =
                SubGraphFreeGraph::from(g.filter_map(&|a| Some((a.0.into(), a.1.into()))));

            let graph = DotGraph::from((can_graph, Figment::new()));
            global_data.push(graph.global_data);
            set.push(graph.graph);
        }

        Ok(GraphSet { set, global_data })
    }

    #[cfg(feature = "serde")]
    #[allow(clippy::result_large_err, clippy::type_complexity)]
    pub fn from_string_with_figment<'a, Str: AsRef<str>>(
        s: Str,
        figment: figment::Figment,
    ) -> Result<Self, HedgeParseError<'a, (), (), (), ()>> {
        let ast_graphs = dot_parser::ast::Graphs::try_from(s.as_ref())?;

        let mut global_data = Vec::new();
        let mut set = Vec::with_capacity(ast_graphs.graphs.len());
        for g in ast_graphs.graphs {
            let can_graph =
                SubGraphFreeGraph::from(g.filter_map(&|a| Some((a.0.into(), a.1.into()))));

            let graph = DotGraph::from((can_graph, figment.clone()));
            global_data.push(graph.global_data);
            set.push(graph.graph);
        }

        Ok(GraphSet { set, global_data })
    }
}

impl<S: NodeStorageOps<NodeData = DotVertexData>> IntoIterator
    for GraphSet<DotEdgeData, DotVertexData, DotHedgeData, GlobalData, S>
{
    type Item = DotGraph<S>;
    type IntoIter = std::iter::Map<
        std::iter::Zip<
            std::vec::IntoIter<HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData, S>>,
            std::vec::IntoIter<GlobalData>,
        >,
        fn(
            (
                HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData, S>,
                GlobalData,
            ),
        ) -> DotGraph<S>,
    >;

    fn into_iter(self) -> Self::IntoIter {
        self.set
            .into_iter()
            .zip(self.global_data)
            .map(|(graph, global_data)| DotGraph { global_data, graph })
    }
}
