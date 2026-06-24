use std::collections::BTreeMap;

use crate::{
    GammaLoopContext,
    cff::CutCFFIndex,
    debug_tags,
    graph::{Graph, LMBext, cuts::CutSet},
    settings::global::OrientationPattern,
    utils::{GS, W_, symbolica_ext::LogPrint},
    uv::approx::{CFFapprox, CutStructure, ForestNodeLike},
};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::eyre;
use gammaloop_tracing_filter::{LogMessage, debug_instrument};
use idenso::{color::ColorSimplifier, shorthands::schoonschip::Schoonschip};
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeVec, Orientation};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};

use linnet::half_edge::{involution::HedgePair, subgraph::SubSetOps};
use vakint::Vakint;

use super::{
    RenormalizationPart, UVgenerationSettings, UltravioletGraph,
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
    pub integrands: BTreeMap<CutCFFIndex, Atom>,
    pub cuts: CutSet,
}

impl ParametricIntegrands {
    pub fn map<F: FnMut(Atom) -> Atom>(self, mut map: F) -> Self {
        Self {
            integrands: self
                .integrands
                .into_iter()
                .map(|(index, atom)| (index, map(atom)))
                .collect(),
            cuts: self.cuts,
        }
    }

    pub fn zero_like(&self) -> Self {
        Self {
            integrands: self
                .integrands
                .keys()
                .map(|index| (*index, Atom::Zero))
                .collect(),
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
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &UVgenerationSettings,
        orientation_pattern: &OrientationPattern,
    ) -> Result<()> {
        for ((forest, cuts), vakint_settings) in &mut self
            .forests
            .iter_mut()
            .zip(self.cuts.cuts.iter())
            .zip(self.settings.iter())
        {
            debug_tags!(#forest,#uv;
                n_terms = %forest.n_terms(),
                "Computing cut forest");
            forest.compute(
                graph,
                (vakint, vakint_settings),
                cuts,
                valid_orientations,
                settings,
                orientation_pattern,
            )?;
        }
        Ok(())
    }
    #[debug_instrument(graph = %graph.log_display())]
    pub(crate) fn orientation_parametric_exprs(
        self,
        graph: &Graph,
        settings: &UVgenerationSettings,
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

        for (i, (forest, cuts)) in forests.iter().zip(cuts.cuts.into_iter()).enumerate() {
            let forest_started = std::time::Instant::now();
            debug_tags!(#generation, #profile, #uv, #graph, #cut, #summary;
                stage = "orientation_parametric_expr_start",
                forest_index = i,
                forest_terms = forest.n_terms(),
                "Building orientation parametric forest"
            );
            let mut integrands = forest.orientation_parametric_expr(graph)?;
            debug_tags!(#generation, #profile, #uv, #graph, #cut, #summary;
                stage = "orientation_parametric_expr_sum_done",
                forest_index = i,
                forest_terms = forest.n_terms(),
                integrand_count = integrands.len(),
                elapsed_ms = forest_started.elapsed().as_secs_f64() * 1000.0,
                "Built orientation parametric forest sum"
            );

            debug_tags!(#generation, #uv, #graph, #dump;
                n_terms =%forest.n_terms(),
                log.graph = %graph.dot(&cuts.union),
                integrands=%integrands.iter().enumerate().map(|(i, s)| format!("{}: {}", i, s.1.log_print(Some(100)))).join("\n"),
                file.integrands = %integrands.iter().map(|s| s.1.to_canonical_string()).join(";"),
                "Orientation Parametric integrand {i}",
            );
            if !settings.keep_sigma {
                integrands.values_mut().for_each(|s| {
                    *s = s
                        .replace(function!(GS.uv_local, W_.a___))
                        .with(Atom::num(1))
                        .replace(function!(GS.uv_integrated, W_.a___))
                        .with(Atom::num(1))
                });
            }
            debug_tags!(#generation, #profile, #uv, #graph, #cut, #summary;
                stage = "orientation_parametric_expr_done",
                forest_index = i,
                forest_terms = forest.n_terms(),
                integrand_count = integrands.len(),
                elapsed_ms = forest_started.elapsed().as_secs_f64() * 1000.0,
                total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                "Built orientation parametric forest"
            );
            exprs.push(ParametricIntegrands { integrands, cuts });
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
        cut_data: &CutSet,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &UVgenerationSettings,
        orientation_pattern: &OrientationPattern,
    ) -> Result<()> {
        let started = std::time::Instant::now();
        debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "forest_compute_start",
            generate_integrated = settings.generate_integrated,
            only_integrated = settings.only_integrated,
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
                    self.dag.nodes[n].data.root(
                        graph,
                        cut_data,
                        valid_orientations,
                        settings,
                        orientation_pattern,
                    )?;
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
                    current.data.set_filter_state(graph, &parent.data, settings);

                    debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                        stage = "computing forest node",
                        simple = %current.data.simple_approx.as_ref().map(|a| a.expr(&graph.full_filter()).to_string()).unwrap(),
                        log.current = current.data.spinney,
                        log.given = parent.data.spinney,
                        "Computing UV forest node single"
                    );

                    current.data.topo_order = i;
                    if settings.generate_integrated {
                        let integrated_started = std::time::Instant::now();
                        debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                            "Computing integrated UV forest node"
                        );
                        current
                            .data
                            .compute_integrated(graph, vakint, &parent.data, settings)?;
                        debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                            elapsed_ms = integrated_started.elapsed().as_secs_f64() * 1000.0,
                            "Computed integrated UV forest node"
                        );
                    }
                    if settings.only_integrated {
                        debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                            "Skipping local UV forest node"
                        );
                        continue;
                    }
                    assert!(matches!(parent.data.local_3d, CFFapprox::Dependent { .. }));
                    let local_started = std::time::Instant::now();
                    debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                        stage = "forest_node_compute_local_3d_start",
                        "Computing local 3D UV forest node"
                    );
                    current.data.compute(
                        graph,
                        cut_data,
                        &parent.data,
                        valid_orientations,
                        settings,
                        orientation_pattern,
                    )?;
                    debug_tags!(#generation, #profile, #uv, #graph, #term, #summary;
                        elapsed_ms = local_started.elapsed().as_secs_f64() * 1000.0,
                        "Computed local 3D UV forest node"
                    );
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

    pub(crate) fn pole_part_of_ends(
        &self,
        graph: &Graph,
        pole_part: bool,
    ) -> Result<RenormalizationPart> {
        let mut sum = Atom::Zero;

        let wild = Atom::var(W_.x___);

        let replacements =
            graph.integrand_replacement(&graph.full_filter(), &graph.loop_momentum_basis, &[wild]);
        for (_, n) in &self.dag.nodes {
            if !n.children.is_empty() {
                continue;
            }

            let atom = n.data.integrated_4d.expr().ok_or(eyre!(
                "Integrated pole part not computed for {} of graph {}",
                n.data
                    .simple_approx
                    .as_ref()
                    .unwrap()
                    .expr(&graph.full_filter()),
                graph.name
            ))?;

            let atom = atom[&CutCFFIndex::new_all_none()].clone();

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

            let forest_term = n
                .data
                .simple_approx
                .as_ref()
                .unwrap()
                .expr(&graph.full_filter());
            let expanded_atom = atom.expand_num();
            debug_tags!(#generation, #uv, #graph, #term;
                log.forest_term = forest_term,
                log.expr = expanded_atom,
                "Term before simplification"
            );

            let atom = atom
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
        if pole_part {
            let n_loops = graph.n_loops(&graph.full_filter());
            let pole_stripped = sum
                .series(GS.dim_epsilon, Atom::Zero, n_loops as i64 + 1)
                .unwrap();

            sum = Atom::Zero;

            for (power, p) in pole_stripped.terms() {
                if power < 0 {
                    sum += p * Atom::var(GS.dim_epsilon).pow(power);
                }
            }
        }
        Ok(RenormalizationPart::legacy(
            sum.replace_multiple(&replacements)
                .replace(GS.m_uv_int)
                .with(GS.m_uv),
        ))
    }

    #[debug_instrument(graph = %graph.log_display())]
    pub(crate) fn orientation_parametric_expr(
        &self,
        graph: &Graph,
    ) -> Result<BTreeMap<CutCFFIndex, Atom>> {
        let mut sum = None;

        for (_, n) in &self.dag.nodes {
            debug_tags!(#generation, #uv, #graph, #term;

                dod = %n.data.dod(),
                log.graph = %graph.dot_lmb_of(&n.data.spinney.subgraph,&n.data.spinney.lmb),
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
                sum = Some(BTreeMap::new());
            }

            let sum = sum.as_mut().unwrap();

            for (cut_index, integrand) in n
                .data
                .final_integrand
                .as_ref()
                .ok_or(eyre!("Final integrand not computed"))?
                .iter()
            {
                let a = integrand.clone().collect_color();

                if first {
                    sum.insert(*cut_index, a);
                } else {
                    *sum.get_mut(cut_index).unwrap() += a;
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

            for s in &mut sum.values_mut() {
                *s = s.replace_multiple(&[GS.split_mom_pattern_simple(edge_index)]);
            }
        }

        for s in &mut sum.values_mut() {
            *s = s.replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_)).with(W_.d_);
            // .collect_factors(); Really bad ! Turns
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
    use std::collections::BTreeMap;

    use super::ParametricIntegrands;
    use crate::{cff::CutCFFIndex, graph::cuts::CutSet};
    use symbolica::{atom::Atom, symbol};

    #[test]
    fn zero_like_preserves_shape_and_cuts() {
        let integrands = ParametricIntegrands {
            integrands: BTreeMap::from([
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
            BTreeMap::from([
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
