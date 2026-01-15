use crate::{
    cff::{
        expression::{GraphOrientation, OrientationData, OrientationID},
        generation::{ConstraintData, ShiftRewrite},
    },
    graph::{Edge, Graph, Vertex},
    momentum::SignOrZero,
    numerator::{ParsingNet, symbolica_ext::AtomCoreExt},
    utils::{GS, ose_atom_from_index, symbolica_ext::CallSymbol},
    uv::approx::CFFapprox,
};
use spenso::{
    network::{Network, store::TensorScalarStoreMapping},
    structure::abstract_index::AIND_SYMBOLS,
    tensors::{data::StorageTensor, parametric::ParamOrConcrete},
};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
    id::Replacement,
};

use linnet::half_edge::{
    HedgeGraph,
    involution::EdgeIndex,
    subgraph::{SuBitGraph, SubGraphLike, SubSetLike},
};
use std::fmt::Write;
use tracing::{debug, instrument};

use typed_index_collections::TiVec;
use vakint::Vakint;

use super::{
    UVgenerationSettings, UltravioletGraph,
    approx::Approximation,
    poset::{DAG, DagNode},
};

pub struct Forest {
    pub dag: DAG<Approximation, DagNode, ()>,
}

impl Forest {
    pub(crate) fn _n_terms(&self) -> usize {
        self.dag.nodes.len()
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn compute<
        H,
        G: UltravioletGraph + AsRef<HedgeGraph<Edge, Vertex, H>>,
        S: SubGraphLike<Base = SuBitGraph>,
        OID: OrientationID,
        O: GraphOrientation,
    >(
        &mut self,
        graph: &G,
        amplitude_subgraph: &S,
        vakint: &Vakint,
        orientations: &TiVec<OID, O>,
        canonize_esurface: &Option<ShiftRewrite>,
        cut_edges: &[EdgeIndex],
        constraint_data: Option<ConstraintData>,
        settings: &UVgenerationSettings,
        _conjugate: bool,
    ) {
        let order = self.dag.compute_topological_order();

        for n in order {
            match self.dag.nodes[n].parents.len() {
                0 => {
                    self.dag.nodes[n].data.root(
                        graph,
                        amplitude_subgraph,
                        canonize_esurface,
                        orientations,
                        cut_edges,
                        constraint_data,
                    );
                }
                1 => {
                    let [current, parent] = &mut self
                        .dag
                        .nodes
                        .get_disjoint_mut([n, self.dag.nodes[n].parents[0]])
                        .unwrap();

                    assert!(matches!(parent.data.local_3d, CFFapprox::Dependent { .. }));
                    if settings.generate_integrated {
                        current.data.compute_integrated(
                            graph,
                            vakint,
                            amplitude_subgraph,
                            &parent.data,
                            settings,
                        );
                    }
                    current.data.compute(
                        graph,
                        amplitude_subgraph,
                        canonize_esurface,
                        orientations,
                        cut_edges,
                        constraint_data,
                        &parent.data,
                        settings,
                    );
                }
                _ => {
                    unimplemented!("Union not implemented");
                    // let (local_4d, simple, local_3d) = {
                    //     let current = &self.dag.nodes[n];
                    //     let mut dependents =
                    //         current.parents.iter().map(|p| &self.dag.nodes[*p].data);

                    //     let mut dep_vec: Vec<&Approximation> = Vec::new();
                    //     dep_vec.push(dependents.next().unwrap());

                    //     let Some((mut approx, _)) = dep_vec[0].local_3d.expr() else {
                    //         panic!("local 3d not computed for parents")
                    //     };

                    //     let mut iterative_approx = dep_vec[0].clone();

                    //     for p in dependents {
                    //         dep_vec.push(p);

                    //         approx = p.local_3d(&iterative_approx, graph, approx);

                    //         iterative_approx = p.clone();
                    //     }
                    //     (
                    //         ApproxOp::union(&dep_vec).unwrap(),
                    //         SimpleApprox::union(
                    //             self.dag.nodes[n].data.subgraph.clone(),
                    //             dep_vec.iter().map(|s| s.simple_approx.as_ref().unwrap()),
                    //         ),
                    //         iterative_approx.local_3d,
                    //     )
                    // };

                    // self.dag.nodes[n].data.local_4d = local_4d;
                    // self.dag.nodes[n].data.simple_approx = Some(simple);
                    // self.dag.nodes[n].data.local_3d = local_3d;
                }
            }
        }
    }

    #[instrument(skip_all)]
    pub(crate) fn orientation_parametric_expr(
        &self,
        cut_edges: Option<&SuBitGraph>,
        graph: &Graph,
        add_sigma: bool,
    ) -> ParsingNet {
        let mut sum: Option<ParsingNet> = None;

        for (_, n) in &self.dag.nodes {
            let mut net = n.data.final_integrand.clone().unwrap();

            if add_sigma {
                debug!(
                    "{}",
                    n.data
                        .simple_approx
                        .as_ref()
                        .unwrap()
                        .expr(&graph.full_filter())
                );
                net = net
                    * Network::from_scalar(function!(
                        GS.if_sigma,
                        n.data
                            .simple_approx
                            .as_ref()
                            .unwrap()
                            .expr(&graph.full_filter())
                    ))
            };

            let Some(s) = &mut sum else {
                sum = Some(net);
                continue;
            };

            *s += net;
        }

        let Some(mut s) = sum else {
            return ParsingNet::zero();
        };
        if let Some(cut) = cut_edges {
            // add Feynman rules of cut edges
            for (_p, edge_index, d) in graph.iter_edges_of(cut) {
                s *= (&d.data.num / (-Atom::num(2) * ose_atom_from_index(edge_index)))
                    .wrap_color(GS.color_wrap)
                    .parse_into_net()
                    .unwrap();

                s = s.replace_multiple(&[GS.add_parametric_sign(edge_index)]);
            }
        }
        s
    }

    pub(crate) fn _local_expr(
        &self,
        orientation: &OrientationData,
        cut_edges: Option<&SuBitGraph>,
        graph: &Graph,
    ) -> ParsingNet {
        let mut sum: Option<ParsingNet> = None;

        for (_, n) in &self.dag.nodes {
            let net = n.data.final_integrand.as_ref().unwrap().map_ref(
                |a| orientation.select(a),
                |a| match a {
                    ParamOrConcrete::Param(a) => {
                        ParamOrConcrete::Param(a.map_data_ref_self(|a| orientation.select(a)))
                    }
                    a => a.clone(),
                },
            );
            let Some(s) = &mut sum else {
                sum = Some(net);
                continue;
            };
            *s += net;
        }

        let Some(mut s) = sum else {
            return ParsingNet::zero();
        };
        if let Some(cut) = cut_edges {
            // add Feynman rules of cut edges
            for (_p, edge_index, d) in graph.iter_edges_of(cut) {
                let edge_id = usize::from(edge_index) as i64;

                // do not set the cut momenta generated in the amplitude to their OSE values
                // yet in order to do 4d scaling tests
                let temp_orientation = orientation.clone();

                s *= d
                    .data
                    .num
                    .wrap_color(GS.color_wrap)
                    .parse_into_net()
                    .unwrap();
                s = s.replace_multiple(&[
                    Replacement::new(
                        function!(GS.emr_mom, edge_id, AIND_SYMBOLS.cind.f(&[Atom::Zero]))
                            .to_pattern(),
                        SignOrZero::from(temp_orientation.orientation[edge_index])
                            * function!(GS.ose, edge_id),
                    ),
                    // Replacement::new(
                    //     function!(GS.color_wrap, W_.a___).to_pattern(),
                    //     Atom::var(W_.a___),
                    // ),
                ]);
            }
        }
        s
    }

    // pub(crate) fn simple_expr(
    //     &self,
    //     graph: &UVGraph,
    //     amplitude: &InternalSubGraph,
    // ) -> Option<SerializableAtom> {
    //     let mut sum = Atom::num(0).into();
    //     for (_, n) in &self.dag.nodes {
    //         sum = sum + n.data.simple_expr(graph, amplitude)?;
    //     }

    //     Some(sum)
    // }

    pub(crate) fn graphs(&self, graph: &Graph) -> String {
        let mut out = String::new();
        self.dag.nodes.iter().for_each(|(_, a)| {
            writeln!(
                out,
                "//S_{}:\n{}",
                a.data.subgraph.string_label(),
                a.data.dot(graph)
            )
            .unwrap()
        });
        writeln!(
            out,
            "{}",
            self.dag
                .to_dot_impl(&|n| format!("label=S_{}", n.data.subgraph.string_label()))
        )
        .unwrap();

        out
    }

    // pub(crate) fn show_structure(&self, graph: &UVGraph, amplitude: &InternalSubGraph) -> Option<String> {
    //     let mut out = String::new();

    //     out.push_str(
    //         &self
    //             .simple_expr(graph, amplitude)?
    //             .0
    //             .printer(PrintOptions {
    //                 terms_on_new_line: true,
    //                 ..Default::default()
    //             })
    //             .to_string(),
    //     );
    //     Some(out)
    // }
}
