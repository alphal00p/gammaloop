use linnet::{
    half_edge::{
        HedgeGraph,
        involution::{Flow, HedgePair},
        subgraph::SubSetLike,
    },
    parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData},
};

use crate::{
    graph::Graph,
    utils::{GS, W_},
};

use super::{ParseGraph, string_utils::ToQuoted};

impl From<&Graph> for DotGraph {
    fn from(value: &Graph) -> Self {
        let global_data = value.global_data();

        let mut graph: HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData> =
            if value.initial_state_cut.is_empty() {
                value.underlying.map_data_ref(
                    |_, _, v| v.into(),
                    |_, _, p, e| {
                        if let HedgePair::Unpaired { flow, .. } = p {
                            e.map(|e| {
                                let mut dot: DotEdgeData = e.into();
                                match flow {
                                    Flow::Sink => {
                                        dot.add_statement("pin", format!("\"x:@right\""));
                                    }
                                    Flow::Source => {
                                        dot.add_statement("pin", format!("\"x:@left\""));
                                    }
                                }
                                dot
                            })
                        } else {
                            e.map(|e| e.into())
                        }
                    },
                    |_, d| d.into(),
                )
            } else {
                value
                    .initial_state_cut
                    .clone()
                    .to_owned_graph_ref(&value.underlying)
                    .map(
                        |_, _, v| v.into(),
                        |_, _, _, _, e| {
                            e.map(|e| {
                                let mut dot: DotEdgeData = (*e.edge_data()).into();
                                match e.flow() {
                                    Some(Flow::Sink) => {
                                        dot.add_statement(
                                            "pin",
                                            format!("\"x:@right,y:@edge{}\"", e.index),
                                        );
                                    }
                                    Some(Flow::Source) => {
                                        dot.add_statement(
                                            "pin",
                                            format!("\"x:@left,y:@edge{}\"", e.index),
                                        );
                                    }
                                    None => {}
                                }
                                dot
                            })
                        },
                        |_, d| d.into(),
                    )
            };

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
