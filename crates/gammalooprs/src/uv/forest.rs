use crate::{
    GammaLoopContext, debug_tags,
    graph::{Graph, LMBext, cuts::CutSet},
    utils::{GS, W_, symbolica_ext::LogPrint},
    uv::{
        Integrands,
        approx::{CutStructure, ForestNodeLike, OrientationProjection, local_3d::Localizer},
        settings::FinalIntegrandDimension,
    },
};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::{WrapErr, eyre};
use gammaloop_tracing_filter::{LogMessage, debug_instrument};
use idenso::{color::ColorSimplifier, shorthands::schoonschip::Schoonschip};

use symbolica::atom::{Atom, AtomCore};

use linnet::half_edge::{involution::HedgePair, subgraph::SubSetOps};
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
    pub integrands: Integrands,
    pub cuts: CutSet,
}

impl ParametricIntegrands {
    pub fn map<F: FnMut(Atom) -> Atom>(self, mut map: F) -> Self {
        Self {
            integrands: self.integrands.map(|atom| map(atom.clone())),
            cuts: self.cuts,
        }
    }

    pub fn zero_like(&self) -> Self {
        Self {
            integrands: self.integrands.map(|_| Atom::Zero),
            cuts: self.cuts.clone(),
        }
    }
}

impl CutForests {
    #[debug_instrument(graph = %graph.log_display())]
    pub(crate) fn compute(
        &mut self,
        graph: &mut Graph,
        vakint: &Vakint,
        orientation: OrientationProjection<'_>,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        for ((forest, cuts), vakint_settings) in &mut self
            .forests
            .iter_mut()
            .zip(self.cuts.cuts.iter())
            .zip(self.settings.iter())
        {
            let localizer = Localizer::new(cuts, orientation);
            debug_tags!(#forest,#uv;
                n_terms = %forest.n_terms(),
                "Computing cut forest");
            forest.compute(graph, (vakint, vakint_settings), localizer, settings)?;
        }
        Ok(())
    }
    #[debug_instrument(graph = %graph.log_display())]
    pub(crate) fn orientation_parametric_exprs(
        self,
        graph: &Graph,
        _settings: &UVgenerationSettings,
    ) -> Result<Vec<ParametricIntegrands>> {
        let started = std::time::Instant::now();
        let CutForests {
            cuts,
            forests,
            settings: _vakint_settings,
        } = self;
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "orientation_parametric_exprs_start",
            forest_count = forests.len(),
            "Building orientation parametric integrands"
        );
        let mut exprs = vec![];

        for (forest, cuts) in forests.iter().zip(cuts.cuts.into_iter()) {
            exprs.push(ParametricIntegrands {
                integrands: forest.orientation_parametric_expr(graph)?,
                cuts,
            });
        }
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "orientation_parametric_exprs_done",
            result_count = exprs.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Built orientation parametric integrands"
        );

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
    #[debug_instrument(graph = %graph.log_display(), forest_terms = self.n_terms())]
    pub(crate) fn compute(
        &mut self,
        graph: &mut Graph,
        vakint: (&Vakint, &vakint::VakintSettings),
        localizer: Localizer<'_>,
        settings: &UVgenerationSettings,
    ) -> Result<()> {
        let started = std::time::Instant::now();
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "forest_compute_start",
            generate_integrated = settings.generate_integrated,
            final_integrand = %settings.final_integrand,
            "Computing UV forest"
        );
        let order = self.dag.compute_topological_order();

        for (i, n) in order.into_iter().enumerate() {
            let node_started = std::time::Instant::now();
            let parent_count = self.dag.nodes[n].parents.len();
            let dod = self.dag.nodes[n].data.spinney.dod;
            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                stage = "forest_node_start",
                node = ?n,
                topo_index = i,
                parent_count,
                dod,
                "Computing UV forest node"
            );
            match self.dag.nodes[n].parents.len() {
                0 => {
                    self.dag.nodes[n].data.topo_order = i;
                    let root_started = std::time::Instant::now();
                    self.dag.nodes[n].data.root(graph, localizer, settings)?;
                    debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                        stage = "forest_node_root_done",
                        elapsed_ms = root_started.elapsed().as_secs_f64() * 1000.0,
                        "Computed root UV forest node"
                    );
                }
                1 => {
                    // debug!("")
                    let parent_id = self.dag.nodes[n].parents[0];
                    let [current, parent] =
                        &mut self.dag.nodes.get_disjoint_mut([n, parent_id]).unwrap();

                    let Some(a) = &parent.data.simple_approx else {
                        panic!("Should have computed the simple approx");
                    };
                    current.data.simple_approx =
                        Some(a.dependent(current.data.spinney.subgraph.clone()));

                    debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                        stage = "computing forest node",
                        simple = %current.data.simple_approx.as_ref().map(|a| a.expr(&graph.full_filter()).to_string()).unwrap(),
                        log.current = current.data.spinney,
                        log.given = parent.data.spinney,
                        "Computing UV forest node single"
                    );

                    current.data.topo_order = i;
                    current
                        .data
                        .compute_4d(graph, vakint, &parent.data, settings)?;

                    match settings.final_integrand {
                        FinalIntegrandDimension::FourD => {
                            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                                "Skipping local UV forest node"
                            );
                            continue;
                        }
                        FinalIntegrandDimension::ThreeD => {
                            current
                                .data
                                .compute_3d(&parent.data, graph, localizer, settings)?;
                        }
                    }
                }
                _ => {
                    unimplemented!("Union not implemented");
                }
            }
            debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                stage = "forest_node_done",
                elapsed_ms = node_started.elapsed().as_secs_f64() * 1000.0,
                "Computed UV forest node"
            );
        }
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "forest_compute_done",
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Computed UV forest"
        );
        Ok(())
    }

    pub(crate) fn renormalization_part_of_ends(
        &self,
        graph: &Graph,
    ) -> Result<RenormalizationPart> {
        let mut sum = Atom::Zero;

        let wild = Atom::var(W_.x___);

        let replacements =
            graph.integrand_replacement(&graph.full_filter(), &graph.loop_momentum_basis, &[wild]);
        for (_, n) in &self.dag.nodes {
            if !n.children.is_empty() {
                continue;
            }

            let atom = n.data.integrated(graph)?.physical_atom();

            let expanded_atom = atom.expand_num();
            debug_tags!(#generation, #uv, #graph, #term;
                forest_term = %n.data.simple_display(graph),
                log.expr = expanded_atom,
                "Term before simplification"
            );

            let atom = &atom
                * &graph.global_prefactor.projector
                * &graph.global_prefactor.num
                * &graph.overall_factor;
            debug_tags!(#generation, #uv, #inspect, #dump;
                log.expression = atom,
                dod = n.data.spinney.dod,
                "Dumped pole part color simplification input"
            );
            let atom = atom.simplify_color().expand_num().to_dots();
            // .replace(GS.dim)
            // .max_level(0)
            // .with(4); //.with(Atom::var(GS.dim_epsilon) * (-2) + 4);

            debug_tags!(#generation, #uv, #graph, #term;
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

        Ok(RenormalizationPart::legacy(
            sum.replace_multiple(&replacements)
                .replace(GS.m_uv_expansion)
                .with(GS.m_uv_vacuum),
        ))
    }

    #[debug_instrument(graph = %graph.log_display())]
    pub(crate) fn orientation_parametric_expr(&self, graph: &Graph) -> Result<Integrands> {
        let mut sum: Option<Integrands> = None;

        for (_, n) in &self.dag.nodes {
            debug_tags!(#generation, #uv, #graph, #term;
                dod = %n.data.dod(),
                log.graph = %graph.dot_lmb_of(&n.data.spinney.subgraph,&n.data.spinney.lmb),
                simple = %n.data.simple_display(graph),"Terms"
            );
            let terms = n
                .data
                .final_integrand(graph)?
                .iter()
                .map(|(cut_index, integrand)| (*cut_index, integrand.clone().collect_color()))
                .collect();
            sum = Some(match sum {
                Some(sum) => sum.zip_add(&terms).wrap_err_with(|| {
                    format!(
                        "while aggregating legacy UV forest term {}",
                        n.data.simple_display(graph)
                    )
                })?,
                None => terms,
            });
        }

        let sum = sum.ok_or(eyre!("No terms in forest"))?;
        let split_momentum_replacements = graph
            .iter_edges_of(
                &graph
                    .full_filter()
                    .subtract(&graph.initial_state_cut)
                    .subtract(&graph.tree_edges),
            )
            .filter_map(|(pair, edge_index, _)| {
                (!matches!(pair, HedgePair::Unpaired { .. }))
                    .then(|| GS.split_mom_pattern_simple(edge_index))
            })
            .collect::<Vec<_>>();

        Ok(sum.map(|integrand| {
            integrand
                .replace_multiple(&split_momentum_replacements)
                .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                .with(W_.d_)
            // .collect_factors(); Really bad ! Turns
        }))
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
    use crate::{cff::CutCFFIndex, graph::cuts::CutSet, uv::Integrands};
    use symbolica::{atom::Atom, symbol};

    #[test]
    fn zero_like_preserves_shape_and_cuts() {
        let integrands = ParametricIntegrands {
            integrands: Integrands::from_iter([
                (CutCFFIndex::new_all_none(), Atom::var(symbol!("x"))),
                (
                    CutCFFIndex {
                        left_threshold_order: Some(1),
                        right_threshold_order: None,
                        lu_cut_order: None,
                    },
                    Atom::num(7),
                ),
            ]),
            cuts: CutSet::empty(3),
        };

        let zeroed = integrands.zero_like();

        assert_eq!(
            zeroed.integrands,
            Integrands::from_iter([
                (CutCFFIndex::new_all_none(), Atom::Zero),
                (
                    CutCFFIndex {
                        left_threshold_order: Some(1),
                        right_threshold_order: None,
                        lu_cut_order: None,
                    },
                    Atom::Zero,
                ),
            ]),
        );
        assert_eq!(zeroed.cuts, CutSet::empty(3));
    }
}
