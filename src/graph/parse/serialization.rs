use linnet::{
    half_edge::{subgraph::SubGraph, HedgeGraph},
    parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData},
};

use crate::{
    graph::Graph,
    utils::{GS, W_},
};

use super::{string_utils::ToQuoted, ParseGraph};

impl From<&Graph> for DotGraph {
    fn from(value: &Graph) -> Self {
        let global_data = value.global_data();
        let mut graph: HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData> =
            value.underlying.map_data_ref(
                |_, _, v| v.into(),
                |_, _, _, e| e.map(|e| e.into()),
                |_, d| d.into(),
            );

        // value.normal_emr_replacement(subgraph, lmb, rep_args, filter_pair)
        for (_, i, _) in value.iter_edges() {
            let loop_expr = value
                .loop_momentum_basis
                .loop_atom(i, GS.loop_mom, &[W_.a___], false);
            let external_expr =
                value
                    .loop_momentum_basis
                    .ext_atom(i, GS.external_mom, &[W_.a___], false);

            graph[i].add_statement("lmb_rep", (loop_expr + external_expr).to_quoted());
        }

        for i in value.initial_state_cut.left.included_iter() {
            let eid = graph[&i];
            graph[eid].add_statement("is_cut", i);
        }

        for (l, i) in value.loop_momentum_basis.loop_edges.iter_enumerated() {
            graph[i].0.add_statement("lmb_id", l);
        }

        DotGraph { global_data, graph }
    }
}

impl From<&ParseGraph> for DotGraph {
    fn from(value: &ParseGraph) -> Self {
        let global_data = value.global_data();
        let graph: HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData> = value.graph.map_data_ref(
            |_, _, v| v.into(),
            |_, _, _, e| e.map(|e| e.into()),
            |_, d| d.into(),
        );

        DotGraph { global_data, graph }
    }
}
