use linnet::{
    half_edge::nodestore::NodeStorageVec,
    parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GraphSet},
};

use color_eyre::Result;

use crate::{graph::Graph, model::Model, processes::AmplitudeGraph};

use super::ParseGraph;

pub trait IntoGraph<T> {
    fn into_graph(self, model: &Model) -> Result<T>;
}

impl IntoGraph<AmplitudeGraph> for String {
    fn into_graph(self, model: &Model) -> Result<AmplitudeGraph> {
        let g: Graph = self.into_graph(model)?;
        Ok(AmplitudeGraph::new(g))
    }
}
impl IntoGraph<AmplitudeGraph> for &str {
    fn into_graph(self, model: &Model) -> Result<AmplitudeGraph> {
        let g: Graph = self.into_graph(model)?;
        Ok(AmplitudeGraph::new(g))
    }
}

impl IntoGraph<Vec<Graph>> for String {
    fn into_graph(self, model: &Model) -> Result<Vec<Graph>> {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_string(self).unwrap();

        Graph::from_hedge_graph_set(hedge_graph_set, model)
    }
}

impl IntoGraph<Graph> for String {
    fn into_graph(self, model: &Model) -> Result<Graph> {
        let graph: DotGraph<NodeStorageVec<DotVertexData>> = DotGraph::from_string(self).unwrap();

        Graph::from_parsed(ParseGraph::from_parsed(graph, model)?, model)
    }
}

impl IntoGraph<Vec<Graph>> for &str {
    fn into_graph(self, model: &Model) -> Result<Vec<Graph>> {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_string(self).unwrap();

        Graph::from_hedge_graph_set(hedge_graph_set, model)
    }
}

impl IntoGraph<Graph> for &str {
    fn into_graph(self, model: &Model) -> Result<Graph> {
        let graph: DotGraph<NodeStorageVec<DotVertexData>> = DotGraph::from_string(self).unwrap();

        Graph::from_parsed(ParseGraph::from_parsed(graph, model)?, model)
    }
}

#[macro_export]
macro_rules! dot {
    // ------------------ Internal Rules (Do not call directly) ------------------

    (@internal [$($code:tt)*], $model:literal) => {

        stringify!($($code)*).into_graph(&$crate::utils::test_utils::load_generic_model($model))
    };

    // Internal rule: End of parsing, with an optional argument.
    // This is matched when the accumulator has collected the code block and we hit a comma.
    (@internal [$($code:tt)*], $model:expr) => {
       stringify!($($code)*).into_graph($model)
    };

    // Internal rule: End of parsing, no optional argument.
    // This is matched when the accumulator has run out of tokens to process.
    (@internal [$($code:tt)*]) => {
        stringify!($($code)*).into_graph(&$crate::utils::test_utils::load_generic_model("sm"))
    };

    // Internal rule: The "accumulator".
    // It takes the next token ($next), adds it to the $code accumulator,
    // and recursively calls the macro with the rest of the tokens ($($rest)*).
    (@internal [$($code:tt)*] $next:tt $($rest:tt)*) => {
        dot!(@internal [$($code)* $next] $($rest)*)
    };

    // ------------------ Public Entry Point ------------------

    // This is the only rule users should call. It kicks off the process by
    // invoking the internal rules with an empty accumulator `[]`.
    ($($all_tokens:tt)+) => {
        dot!(@internal [] $($all_tokens)+)
    };
}
