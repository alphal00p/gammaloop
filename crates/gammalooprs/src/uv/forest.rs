use crate::{
    GammaLoopContext,
    graph::{Graph, LMBext, cuts::CutSet},
    utils::{GS, W_, symbolica_ext::LogPrint},
    uv::approx::{CFFapprox, CutStructure, ForestNodeLike},
};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::eyre;
use idenso::{color::ColorSimplifier, metric::MetricSimplifier};
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeVec, Orientation};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};

use linnet::half_edge::{involution::HedgePair, subgraph::SubSetOps};
use tracing::{debug, instrument};

use vakint::Vakint;

use super::{
    RenormalizationPart, UVgenerationSettings,
    approx::Approximation,
    poset::{DAG, DagNode},
};

pub struct CutForests {
    pub cuts: CutStructure,
    pub forests: Vec<Forest>,
    pub settings: Vec<vakint::VakintSettings>,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ParametricIntegrands {
    pub integrands: Vec<Atom>,
    pub cuts: CutSet,
}

impl ParametricIntegrands {
    pub fn map<F: FnMut(Atom) -> Atom>(self, map: F) -> Self {
        Self {
            integrands: self.integrands.into_iter().map(map).collect(),
            cuts: self.cuts,
        }
    }

    pub fn zero_like(&self) -> Self {
        Self {
            integrands: vec![Atom::Zero; self.integrands.len()],
            cuts: self.cuts.clone(),
        }
    }
}

impl CutForests {
    pub(crate) fn compute(
        &mut self,
        graph: &mut Graph,
        vakint: &Vakint,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        for ((forest, cuts), vakint_settings) in &mut self
            .forests
            .iter_mut()
            .zip(self.cuts.cuts.iter())
            .zip(self.settings.iter())
        {
            forest.compute(
                graph,
                (vakint, vakint_settings),
                cuts,
                valid_orientations,
                settings,
            )?;
        }
        Ok(())
    }
    #[instrument(skip_all)]
    pub(crate) fn orientation_parametric_exprs(
        self,
        graph: &Graph,
        add_sigma: bool,
    ) -> Result<Vec<ParametricIntegrands>> {
        let mut exprs = vec![];

        for (i, (forest, cuts)) in self
            .forests
            .iter()
            .zip(self.cuts.cuts.into_iter())
            .enumerate()
        {
            let integrands = forest.orientation_parametric_expr(graph, add_sigma)?;

            debug!(integrands=%integrands.iter().map(|s| s.to_canonical_string()).join("\n\n"),
                "Orientation Parametric integrand {i},with {} terms for \n{}\n{}",
                forest.n_terms(),
                graph.dot(&cuts.union),
                integrands.iter().map(|s| s.log_print(Some(100))).join("\n"),

            );
            exprs.push(ParametricIntegrands { integrands, cuts });
        }
        Ok(exprs)
    }
}

pub struct Forest {
    pub dag: DAG<Approximation, DagNode, ()>,
}

impl Forest {
    pub(crate) fn n_terms(&self) -> usize {
        self.dag.nodes.len()
    }

    #[allow(clippy::too_many_arguments)]
    #[instrument(skip_all)]
    pub(crate) fn compute(
        &mut self,
        graph: &mut Graph,
        vakint: (&Vakint, &vakint::VakintSettings),
        cut_data: &CutSet,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        let order = self.dag.compute_topological_order();

        for (i, n) in order.into_iter().enumerate() {
            match self.dag.nodes[n].parents.len() {
                0 => {
                    self.dag.nodes[n].data.topo_order = i;
                    self.dag.nodes[n]
                        .data
                        .root(graph, cut_data, valid_orientations, settings)?;
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
                    current.data.simple_approx =
                        Some(a.dependent(current.data.spinney.subgraph.clone()));
                    current.data.update_filtered_integrated_uv_chain_state(
                        graph,
                        &parent.data,
                        settings,
                    );

                    current.data.topo_order = i;
                    if settings.generate_integrated {
                        current
                            .data
                            .compute_integrated(graph, vakint, &parent.data, settings)?;
                    }
                    if settings.only_integrated {
                        continue;
                    }
                    assert!(matches!(parent.data.local_3d, CFFapprox::Dependent { .. }));
                    current.data.compute(
                        graph,
                        cut_data,
                        &parent.data,
                        valid_orientations,
                        settings,
                    )?;
                }
                _ => {
                    unimplemented!("Union not implemented");
                }
            }
        }
        Ok(())
    }

    pub(crate) fn pole_part_of_ends(&self, graph: &Graph) -> Result<RenormalizationPart> {
        let mut sum = Atom::Zero;

        let wild = Atom::var(W_.x___);

        let replacements =
            graph.integrand_replacement(&graph.full_filter(), &graph.loop_momentum_basis, &[wild]);
        for (_, n) in &self.dag.nodes {
            if !n.children.is_empty() {
                continue;
            }

            let (atom, sign) = n.data.integrated_4d.expr().ok_or(eyre!(
                "Integrated pole part not computed for {} of graph {}",
                n.data
                    .simple_approx
                    .as_ref()
                    .unwrap()
                    .expr(&graph.full_filter()),
                graph.name
            ))?;

            let mut atom = atom[0].clone();

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
               expr = % atom.expand_num().log_print(None),"Term before simplification"
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
               expr = % atom.log_print(None),"Term"
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
        Ok(RenormalizationPart::legacy(
            sum.replace_multiple(&replacements)
                .replace(GS.m_uv_int)
                .with(GS.m_uv),
        ))
    }

    #[instrument(skip_all)]
    pub(crate) fn orientation_parametric_expr(
        &self,
        graph: &Graph,
        add_sigma: bool,
    ) -> Result<Vec<Atom>> {
        let mut sum = None;

        for (_, n) in &self.dag.nodes {
            debug!(
                dod = %n.data.dod(),
                simple = %
                n.data
                    .simple_approx
                    .as_ref()
                    .unwrap()
                    .expr(&graph.full_filter()),"Terms"
            );
            let mut first = false;

            if sum.is_none() {
                first = true;
                sum = Some(vec![]);
            }

            let sum = sum.as_mut().unwrap();

            for (i, integrand) in n
                .data
                .final_integrand
                .as_ref()
                .ok_or(eyre!("Final integrand not computed"))?
                .iter()
                .enumerate()
            {
                let a = if add_sigma {
                    debug!(
                        "{}",
                        n.data
                            .simple_approx
                            .as_ref()
                            .unwrap()
                            .expr(&graph.full_filter())
                    );
                    integrand
                        * function!(
                            GS.if_sigma,
                            n.data
                                .simple_approx
                                .as_ref()
                                .unwrap()
                                .expr(&graph.full_filter())
                        )
                } else {
                    integrand.clone()
                };
                if first {
                    sum.push(a);
                } else {
                    sum[i] += a;
                }
            }
        }

        let mut sum = sum.ok_or(eyre!("No terms in forest"))?;

        for (pair, edge_index, _) in graph.iter_edges_of(
            &graph
                .full_filter()
                .subtract(&graph.initial_state_cut)
                .subtract(&graph.tree_edges),
        ) {
            if matches!(pair, HedgePair::Unpaired { .. }) {
                continue;
            }

            for s in &mut sum {
                *s = s.replace_multiple(&[GS.split_mom_pattern_simple(edge_index)]);
            }
        }

        for s in &mut sum {
            *s = s
                .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                .with(W_.d_)
                .collect_factors();
        }
        Ok(sum)
    }

    // pub(crate) fn graphs(&self, graph: &Graph) -> String {
    //     let mut out = String::new();
    //     self.dag.nodes.iter().for_each(|(_, a)| {
    //         writeln!(
    //             out,
    //             "//S_{}:\n{}",
    //             a.data.subgraph.string_label(),
    //             a.data.dot(graph)
    //         )
    //         .unwrap()
    //     });
    //     writeln!(
    //         out,
    //         "{}",
    //         self.dag
    //             .to_dot_impl(&|n| format!("label=S_{}", n.data.subgraph.string_label()))
    //     )
    //     .unwrap();

    //     out
    // }
}

#[cfg(test)]
mod tests {
    use super::ParametricIntegrands;
    use crate::graph::cuts::CutSet;
    use symbolica::{atom::Atom, symbol};

    #[test]
    fn zero_like_preserves_shape_and_cuts() {
        let integrands = ParametricIntegrands {
            integrands: vec![Atom::var(symbol!("x")), Atom::num(7)],
            cuts: CutSet::empty(3),
        };

        let zeroed = integrands.zero_like();

        assert_eq!(zeroed.integrands, vec![Atom::Zero, Atom::Zero]);
        assert_eq!(zeroed.cuts, CutSet::empty(3));
    }
}
