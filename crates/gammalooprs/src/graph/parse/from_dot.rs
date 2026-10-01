use linnet::{
    half_edge::nodestore::NodeStorageVec,
    parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData, GraphSet},
};

use color_eyre::Result;

use crate::{graph::Graph, model::Model, processes::AmplitudeGraph};

use super::ParseGraph;

/// Import a fully finalized amplitude runtime DOT artifact.
pub trait IntoFinalizedRuntimeGraph<T> {
    fn into_finalized_runtime_graph(self, model: &Model) -> Result<T>;
}

impl IntoFinalizedRuntimeGraph<AmplitudeGraph> for String {
    fn into_finalized_runtime_graph(self, model: &Model) -> Result<AmplitudeGraph> {
        let g: Graph = self.into_finalized_runtime_graph(model)?;
        Ok(AmplitudeGraph::new(g))
    }
}
impl IntoFinalizedRuntimeGraph<AmplitudeGraph> for &str {
    fn into_finalized_runtime_graph(self, model: &Model) -> Result<AmplitudeGraph> {
        let g: Graph = self.into_finalized_runtime_graph(model)?;
        Ok(AmplitudeGraph::new(g))
    }
}

impl IntoFinalizedRuntimeGraph<Vec<Graph>> for String {
    fn into_finalized_runtime_graph(self, model: &Model) -> Result<Vec<Graph>> {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_string(self)
            .map_err(|error| color_eyre::eyre::eyre!("DOT parsing error: {error}"))?;

        Graph::from_finalized_runtime_graph_set(hedge_graph_set, model)
    }
}

impl IntoFinalizedRuntimeGraph<Graph> for String {
    fn into_finalized_runtime_graph(self, model: &Model) -> Result<Graph> {
        let graph: DotGraph<NodeStorageVec<DotVertexData>> = DotGraph::from_string(self)
            .map_err(|error| color_eyre::eyre::eyre!("DOT parsing error: {error}"))?;

        Graph::from_parsed_with_validation(ParseGraph::from_parsed(graph, model)?, model)
    }
}

impl IntoFinalizedRuntimeGraph<Vec<Graph>> for &str {
    fn into_finalized_runtime_graph(self, model: &Model) -> Result<Vec<Graph>> {
        let hedge_graph_set: GraphSet<
            DotEdgeData,
            DotVertexData,
            DotHedgeData,
            linnet::parser::GlobalData,
            NodeStorageVec<DotVertexData>,
        > = GraphSet::from_string(self)
            .map_err(|error| color_eyre::eyre::eyre!("DOT parsing error: {error}"))?;

        Graph::from_finalized_runtime_graph_set(hedge_graph_set, model)
    }
}

impl IntoFinalizedRuntimeGraph<Graph> for &str {
    fn into_finalized_runtime_graph(self, model: &Model) -> Result<Graph> {
        let graph: DotGraph<NodeStorageVec<DotVertexData>> = DotGraph::from_string(self)
            .map_err(|error| color_eyre::eyre::eyre!("DOT parsing error: {error}"))?;

        Graph::from_parsed_with_validation(ParseGraph::from_parsed(graph, model)?, model)
    }
}

#[macro_export]
macro_rules! finalized_runtime_dot {
    // ------------------ Internal Rules (Do not call directly) ------------------

    (@internal [$($code:tt)*], $model:literal) => {

        stringify!($($code)*).into_finalized_runtime_graph(&$crate::utils::load_generic_model($model))
    };

    // Internal rule: End of parsing, with an optional argument.
    // This is matched when the accumulator has collected the code block and we hit a comma.
    (@internal [$($code:tt)*], $model:expr) => {
       stringify!($($code)*).into_finalized_runtime_graph($model)
    };

    // Internal rule: End of parsing, no optional argument.
    // This is matched when the accumulator has run out of tokens to process.
    (@internal [$($code:tt)*]) => {
        stringify!($($code)*).into_finalized_runtime_graph(&$crate::utils::load_generic_model("sm"))
    };

    // Internal rule: The "accumulator".
    // It takes the next token ($next), adds it to the $code accumulator,
    // and recursively calls the macro with the rest of the tokens ($($rest)*).
    (@internal [$($code:tt)*] $next:tt $($rest:tt)*) => {
        finalized_runtime_dot!(@internal [$($code)* $next] $($rest)*)
    };

    // ------------------ Public Entry Point ------------------

    // This is the only rule users should call. It kicks off the process by
    // invoking the internal rules with an empty accumulator `[]`.
    ($($all_tokens:tt)+) => {
        finalized_runtime_dot!(@internal [] $($all_tokens)+)
    };
}

#[cfg(test)]
mod tests {
    use crate::{graph::Graph, utils::load_generic_model};

    use super::IntoFinalizedRuntimeGraph;

    #[test]
    fn finalized_runtime_dot_requires_explicit_ufo_slots() {
        let model = load_generic_model("scalars");
        let error = IntoFinalizedRuntimeGraph::<Graph>::into_finalized_runtime_graph(
            r#"digraph strict_amplitude {
                graph [projector="1"];
                node [num="1"];
                edge [pdg=1000, num="1", dir=none];
                a -> b;
            }"#,
            &model,
        )
        .err()
        .expect("missing UFO slots must reject a finalized runtime artifact");

        assert!(error.to_string().contains("missing ufo_order"));
    }

    #[test]
    fn finalized_runtime_dot_rejects_cross_section_sewing_without_canonical_cuts() {
        let model = load_generic_model("scalars");
        let error = IntoFinalizedRuntimeGraph::<Graph>::into_finalized_runtime_graph(
            r#"digraph strict_cross_section {
                graph [projector="1"];
                node [num="1"];
                edge [pdg=1000, num="1", dir=none];
                a -> b [is_cut=0];
            }"#,
            &model,
        )
        .err()
        .expect("cross-section runtime DOT without canonical cuts must be rejected");

        assert!(error.to_string().contains("no canonical physical cuts"));
        assert!(error.to_string().contains("FeynmanDiagram DOT"));
    }
}
