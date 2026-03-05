use crate::{
    cff::{
        expression::{GraphOrientation, OrientationData, OrientationID},
        generation::{PostProcessingSetup, ShiftRewrite},
    },
    graph::{Edge, Graph, LMBext, Vertex},
    momentum::SignOrZero,
    numerator::{ParsingNet, symbolica_ext::AtomCoreExt},
    utils::{
        GS, W_, ose_atom_from_index,
        symbolica_ext::{CallSymbol, LogPrint},
    },
    uv::approx::CFFapprox,
};
use color_eyre::Result;
use eyre::eyre;
use idenso::{color::ColorSimplifier, metric::MetricSimplifier};
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
    subgraph::{SuBitGraph, SubGraphLike, SubSetLike, SubSetOps},
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
        S: SubGraphLike<Base = SuBitGraph> + SubSetOps,
        OID: OrientationID,
        O: GraphOrientation,
    >(
        &mut self,
        graph: &G,
        tree_edges: &S,
        amplitude_subgraph: &S,
        vakint: (&Vakint, &vakint::VakintSettings),
        orientations: &TiVec<OID, O>,
        canonize_esurface: &Option<ShiftRewrite>,
        cut_edges: &[EdgeIndex],
        edges_in_initial_state_cut: &[EdgeIndex],
        post_processing: PostProcessingSetup<'_>,
        settings: &UVgenerationSettings,
        _conjugate: bool,
    ) -> Result<()> {
        let order = self.dag.compute_topological_order();

        for (i, n) in order.into_iter().enumerate() {
            match self.dag.nodes[n].parents.len() {
                0 => {
                    self.dag.nodes[n].data.root(
                        graph,
                        tree_edges,
                        amplitude_subgraph,
                        canonize_esurface,
                        orientations,
                        cut_edges,
                        edges_in_initial_state_cut,
                        post_processing,
                    );
                }
                1 => {
                    // debug!("")
                    let [current, parent] = &mut self
                        .dag
                        .nodes
                        .get_disjoint_mut([n, self.dag.nodes[n].parents[0]])
                        .unwrap();

                    let Some(a) = &parent.data.simple_approx else {
                        panic!("Should have computed the simple approx");
                    };
                    current.data.simple_approx = Some(a.dependent(current.data.subgraph.clone()));

                    if settings.generate_integrated {
                        current.data.compute_integrated(
                            graph,
                            vakint,
                            amplitude_subgraph,
                            &parent.data,
                            i,
                            settings,
                        )?;
                    }
                    if settings.only_integrated {
                        continue;
                    }
                    assert!(matches!(parent.data.local_3d, CFFapprox::Dependent { .. }));
                    current.data.compute(
                        graph,
                        tree_edges,
                        amplitude_subgraph,
                        canonize_esurface,
                        orientations,
                        cut_edges,
                        edges_in_initial_state_cut,
                        post_processing,
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
        Ok(())
    }

    pub(crate) fn pole_part_of_ends(&self, graph: &Graph) -> Result<Atom> {
        let mut sum = Atom::Zero;

        let wild = Atom::var(W_.x___);

        let replacements =
            graph.integrand_replacement(&graph.full_filter(), &graph.loop_momentum_basis, &[wild]);
        for (_, n) in &self.dag.nodes {
            if !n.children.is_empty() {
                continue;
            }

            let (mut atom, sign) = n.data.integrated_pole_part.expr().ok_or(eyre!(
                "Integrated pole part not computed for {} of graph {}",
                n.data
                    .simple_approx
                    .as_ref()
                    .unwrap()
                    .expr(&graph.full_filter()),
                graph.name
            ))?;

            atom = sign * atom;

            // debug!(
            //     forest_term=%
            //     n.data
            //         .simple_approx
            //         .as_ref()
            //         .unwrap()
            //         .expr(&graph.full_filter()),
            //    expr = % atom,"Term before simplification"
            // );
            //
            debug!(
                forest_term=%
                n.data
                    .simple_approx
                    .as_ref()
                    .unwrap()
                    .expr(&graph.full_filter()),
               expr = % atom.expand_num().log_print(),"Term before simplification"
            );

            let atom = (atom
                * &graph.global_prefactor.projector
                * &graph.global_prefactor.num
                * &graph.overall_factor)
                .simplify_color()
                .expand_num()
                .to_dots();
            // .replace(GS.dim)
            // .max_level(0)
            // .with(4); //.with(Atom::var(GS.dim_epsilon) * (-2) + 4);

            debug!(
                forest_term=%
                n.data
                    .simple_approx
                    .as_ref()
                    .unwrap()
                    .expr(&graph.full_filter()),
               expr = % atom.log_print(),"Term"
            );
            sum += atom;
        }
        // let n_loops = graph.n_loops(&graph.full_filter());
        // let pole_stripped = sum
        //     .series(
        //         GS.dim_epsilon,
        //         Atom::Zero,
        //         (n_loops as i64 + 1).into(),
        //         true,
        //     )
        //     .unwrap();

        // let mut sum = Atom::Zero;

        // for (power, p) in pole_stripped.terms() {
        //     if power < 0 {
        //         sum += p * Atom::var(GS.dim_epsilon).pow(power);
        //     }
        // }
        Ok(sum.replace_multiple(&replacements))
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
            }
        }
        for (_, edge_index, _) in graph.iter_edges() {
            s = s.replace_multiple(&[GS.add_parametric_sign(edge_index)]);
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
