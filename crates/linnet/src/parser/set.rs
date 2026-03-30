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

#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Archive, rkyv::Serialize, rkyv::Deserialize)
)]
pub struct GraphSet<E, V, H, G, S: NodeStorage<NodeData = V>> {
    pub global_data: Vec<G>,
    pub set: Vec<HedgeGraph<E, V, H, S>>,
}

pub type DotGraphSet =
    GraphSet<DotEdgeData, DotVertexData, DotHedgeData, GlobalData, DefaultNodeStore<DotVertexData>>;

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

impl<E, V, H, G, S: NodeStorage<NodeData = V>> GraphSet<E, V, H, G, S> {
    #[cfg(feature = "rkyv")]
    pub fn to_rkyv_bytes<const N: usize>(&self) -> Result<rkyv::AlignedVec, String>
    where
        Self: rkyv::Serialize<rkyv::ser::serializers::AllocSerializer<N>>,
    {
        rkyv::to_bytes::<_, N>(self).map_err(|err| err.to_string())
    }

    #[cfg(feature = "rkyv")]
    pub unsafe fn archived_from_bytes<'a>(bytes: &'a [u8]) -> &'a <Self as rkyv::Archive>::Archived
    where
        Self: rkyv::Archive,
    {
        unsafe { rkyv::archived_root::<Self>(bytes) }
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

#[cfg(all(test, feature = "rkyv"))]
mod tests {
    use super::DotGraphSet;

    #[test]
    fn dot_graph_set_archives_without_reparse() {
        let graphs =
            DotGraphSet::from_string(r#"digraph first { a -> b } digraph second { x -> y }"#)
                .unwrap();

        let bytes = graphs.to_rkyv_bytes::<1024>().unwrap();
        let archived = unsafe { DotGraphSet::archived_from_bytes(&bytes) };

        assert_eq!(archived.global_data.len(), 2);
        let names = archived.global_data.as_slice();
        assert_eq!(names[0].name.as_str(), "first");
        assert_eq!(names[1].name.as_str(), "second");
    }
}
