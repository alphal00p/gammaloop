use idenso::{color::ColorSimplifier, gamma::GammaSimplifier};
use linnet::{
    half_edge::{
        HedgeGraph,
        involution::{Flow, HedgePair},
        subgraph::{Inclusion, SubSetLike},
    },
    parser::{DotEdgeData, DotGraph, DotHedgeData, DotVertexData},
};

use crate::{
    graph::Graph,
    processes::DotExportSettings,
    utils::{GS, W_},
    uv::UltravioletGraph,
};

use super::{ParseGraph, string_utils::ToQuoted};

impl Graph {
    pub fn to_dot_graph_with_settings(&self, settings: &DotExportSettings) -> DotGraph {
        let mut dotgraph = if settings.split_xs_by_initial_states {
            self.to_split_dotgraph()
        } else {
            DotGraph::from(self)
        };

        if settings.output_full_numerator {
            let mut num = self
                .numerator(&self.full_filter())
                .get_single_atom()
                .unwrap();

            if settings.do_color_algebra {
                num = num.simplify_color();
            }

            if settings.do_gamma_algebra {
                num = num.simplify_gamma();
            }

            dotgraph
                .global_data
                .statements
                .insert("full_num".into(), num.to_quoted());
        }

        dotgraph
    }

    pub fn to_split_dotgraph(&self) -> DotGraph {
        let global_data = self.global_data();

        let mut graph: HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData> = if self
            .initial_state_cut
            .is_empty()
        {
            self.underlying.map_data_ref(
                |_, _, v| v.into(),
                |_, _, p, e| {
                    if let HedgePair::Unpaired { flow, .. } = p {
                        e.map(|e| {
                            let mut dot: DotEdgeData = e.into();
                            match flow {
                                Flow::Sink => {
                                    dot.add_statement("pin", format!("\"x:@+right\""));
                                }
                                Flow::Source => {
                                    dot.add_statement("pin", format!("\"x:@-left\""));
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
            self.initial_state_cut
                .clone()
                .to_owned_graph_ref(&self.underlying)
                .map(
                    |_, _, v| v.into(),
                    |_, _, p, _, e| {
                        e.map(|e| {
                            let mut dot: DotEdgeData = (*e.edge_data()).into();

                            match e.flow() {
                                Some(Flow::Sink) => {
                                    assert!(
                                        self.initial_state_cut
                                            .left
                                            .includes(&self.inv(p.any_hedge()))
                                    );
                                    dot.add_statement("is_cut", self.inv(p.any_hedge()));
                                    dot.add_statement(
                                        "pin",
                                        format!("\"x:@+right,y:@edge{}\"", e.index),
                                    );
                                }
                                Some(Flow::Source) => {
                                    assert!(self.initial_state_cut.left.includes(&p.any_hedge()));
                                    dot.add_statement("is_cut", p.any_hedge());
                                    dot.add_statement(
                                        "pin",
                                        format!("\"x:@-left,y:@edge{}\"", e.index),
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

        // self.normal_emr_replacement(subgraph, lmb, rep_args, filter_pair)
        for (_, i, _) in self.iter_edges() {
            let loop_expr = self
                .loop_momentum_basis
                .loop_atom(i, GS.loop_mom, &[W_.a___], false);
            let external_expr =
                self.loop_momentum_basis
                    .ext_atom(i, GS.external_mom, &[W_.a___], false);

            graph[i].add_statement("lmb_rep", (loop_expr + external_expr).to_quoted());
        }

        for (l, i) in self.loop_momentum_basis.loop_edges.iter_enumerated() {
            graph[i].0.add_statement("lmb_id", l);
        }

        DotGraph { global_data, graph }
    }
}

impl From<&Graph> for DotGraph {
    fn from(value: &Graph) -> Self {
        let global_data = value.global_data();

        let mut graph: HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData> =
            value.underlying.map_data_ref(
                |_, _, v| v.into(),
                |_, _, p, e| {
                    if let HedgePair::Unpaired { .. } = p {
                        e.map(|e| e.into())
                    } else {
                        e.map(|e| e.into())
                    }
                },
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

impl ParseGraph {
    fn to_simple_dot(&self) -> DotGraph {
        let global_data = self.global_data();
        let graph: HedgeGraph<DotEdgeData, DotVertexData, DotHedgeData> = self.graph.map_data_ref(
            |_, _, v| v.into(),
            |_, _, _, e| e.map(|e| e.into()),
            |_, d| d.into(),
        );

        DotGraph { global_data, graph }
    }
}
