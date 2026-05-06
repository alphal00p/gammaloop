use crate::{
    GammaLoopContext,
    cff::expression::{CFFExpression, GammaLoopGraphOrientation, OrientationID},
    graph::{Graph, cuts::CutSet},
    settings::global::{GenerationSettings, ThreeDRepresentation},
    utils::{GS, W_, symbolica_ext::LogPrint},
    uv::approx::{CFFapprox, CutStructure, ForestNodeLike},
};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeVec, Orientation};
use std::collections::BTreeSet;
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};

use linnet::half_edge::{
    involution::HedgePair,
    subgraph::{SubSetLike, SubSetOps},
};
use tracing::{debug, instrument};

use vakint::Vakint;

use super::{
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
    pub explicitly_summed_orientations: bool,
}

impl ParametricIntegrands {
    pub fn map<F: FnMut(Atom) -> Atom>(self, map: F) -> Self {
        Self {
            integrands: self.integrands.into_iter().map(map).collect(),
            cuts: self.cuts,
            explicitly_summed_orientations: self.explicitly_summed_orientations,
        }
    }

    pub fn zero_like(&self) -> Self {
        Self {
            integrands: vec![Atom::Zero; self.integrands.len()],
            cuts: self.cuts.clone(),
            explicitly_summed_orientations: self.explicitly_summed_orientations,
        }
    }

    pub fn sum_orientations_explicitly(self, orientations: &[EdgeVec<Orientation>]) -> Self {
        let mut summed = self.map(|atom| explicit_orientation_sum_atom(&atom, orientations));
        summed.explicitly_summed_orientations = true;
        summed
    }
}

pub(crate) fn explicit_orientation_sum_atom(
    atom: &Atom,
    orientations: &[EdgeVec<Orientation>],
) -> Atom {
    orientations
        .iter()
        .map(|orientation| orientation.select_gs(atom.as_atom_view()))
        .fold(Atom::Zero, |acc, term| acc + term)
        .collect_factors()
}

fn uv_reduced_subgraph_has_duplicate_leading_denominators(
    graph: &Graph,
    current: &Approximation,
    parent: &Approximation,
) -> bool {
    // Local UV rescaling drops external shifts and replaces all masses in the
    // reduced subgraph by mUV. Equal loop signatures up to an overall sign
    // therefore become confluent denominators even when the original graph did
    // not have a repeated edge group.
    let reduced = current.reduced_subgraph(parent);
    let mut seen = BTreeSet::new();

    for (pair, edge_id, _) in graph.iter_edges_of(&reduced) {
        if !matches!(pair, HedgePair::Paired { .. }) {
            continue;
        }

        let mut leading_signature = current.lmb().edge_signatures[edge_id]
            .internal
            .to_momtrop_format();
        if leading_signature.iter().all(|sign| *sign == 0) {
            continue;
        }
        if leading_signature
            .iter()
            .find(|sign| **sign != 0)
            .is_some_and(|first_non_zero| *first_non_zero < 0)
        {
            for sign in &mut leading_signature {
                *sign = -*sign;
            }
        }

        if !seen.insert(leading_signature) {
            return true;
        }
    }

    false
}

/// Decide whether a forest term must be built from the 4D-expanded integrand
/// before projecting to CFF/LTD. This includes the explicit user setting, LTD
/// generation, and explicit-orientation-sum CFF cases where raised/repeated
/// residues make a local expansion of the already-projected 3D expression
/// representation-dependent.
pub(crate) fn uses_expanded_4d_forest_path(
    graph: &Graph,
    cut_data: &CutSet,
    settings: &GenerationSettings,
    representation: ThreeDRepresentation,
) -> bool {
    let explicit_sum_cff_residue_needs_expanded_4d_local_uv = settings
        .explicit_orientation_sum_only
        && representation == ThreeDRepresentation::Cff
        && (cut_data
            .residue_selector
            .lu_cut
            .as_ref()
            .is_some_and(|lu_cut| lu_cut.max_occurence > 1)
            || graph
                .get_raised_edge_groups()
                .iter()
                .any(|group| group.len() > 1));

    (settings.uv.subtract_uv
        && (settings.uv.local_uv_cts_from_expanded_4d_integrands
            || explicit_sum_cff_residue_needs_expanded_4d_local_uv))
        || representation == ThreeDRepresentation::Ltd
}

pub(crate) fn cff_explicit_sum_needs_outer_orientation_projection(
    graph: &Graph,
    cut_data: &CutSet,
    settings: &GenerationSettings,
    representation: ThreeDRepresentation,
) -> bool {
    settings.explicit_orientation_sum_only
        && representation == ThreeDRepresentation::Cff
        && !uses_expanded_4d_forest_path(graph, cut_data, settings, representation)
}

impl CutForests {
    pub(crate) fn compute(
        &mut self,
        graph: &mut Graph,
        vakint: &Vakint,
        valid_orientations: &[EdgeVec<Orientation>],
        settings: &GenerationSettings,
        root_expression: Option<&CFFExpression<OrientationID>>,
        representation: ThreeDRepresentation,
        explicit_orientation_sum_only: bool,
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
                root_expression,
                representation,
                explicit_orientation_sum_only,
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
            exprs.push(ParametricIntegrands {
                integrands,
                cuts,
                explicitly_summed_orientations: forest.explicitly_summed_orientations,
            });
        }
        Ok(exprs)
    }
}

pub struct Forest {
    pub dag: DAG<Approximation, DagNode, ()>,
    pub explicitly_summed_orientations: bool,
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
        settings: &GenerationSettings,
        root_expression: Option<&CFFExpression<OrientationID>>,
        representation: ThreeDRepresentation,
        explicit_orientation_sum_only: bool,
    ) -> Result<()> {
        // Raised LU residues and repeated loop denominators differentiate the
        // residue-selected 3D expression. The local-UV expansion of that
        // already-differentiated CFF object does not commute with the 4D local
        // forest in general. Use the expanded-4D bridge for these cases so CFF
        // and LTD project the same local counterterm before
        // representation-specific residue evaluation.
        let use_expanded_4d_local_uv =
            uses_expanded_4d_forest_path(graph, cut_data, settings, representation)
                || self.cff_explicit_sum_has_uv_leading_duplicate_denominators(
                    graph,
                    settings,
                    representation,
                );
        self.explicitly_summed_orientations = settings.explicit_orientation_sum_only
            && representation == ThreeDRepresentation::Cff
            && use_expanded_4d_local_uv;
        let order = self.dag.compute_topological_order();

        for (i, n) in order.into_iter().enumerate() {
            match self.dag.nodes[n].parents.len() {
                0 => {
                    self.dag.nodes[n].data.topo_order = i;
                    if use_expanded_4d_local_uv {
                        self.dag.nodes[n].data.root_expanded_4d(
                            graph,
                            cut_data,
                            valid_orientations,
                            settings,
                            representation,
                            root_expression,
                        )?;
                    } else {
                        self.dag.nodes[n].data.root(
                            graph,
                            cut_data,
                            valid_orientations,
                            settings,
                            root_expression,
                        )?;
                    }
                }
                1 => {
                    // debug!("")
                    let parent_id = self.dag.nodes[n].parents[0];
                    let [current, parent] = &mut self
                        .dag
                        .nodes
                        .get_disjoint_mut([n, parent_id])
                        .ok_or_else(|| {
                            eyre!(
                                "UV forest node {n:?} and parent {parent_id:?} could not be borrowed independently"
                            )
                        })?;

                    let Some(a) = &parent.data.simple_approx else {
                        return Err(eyre!(
                            "UV forest parent {parent_id:?} simple approximation was not computed before node {n:?}"
                        ));
                    };
                    current.data.simple_approx =
                        Some(a.dependent(current.data.spinney.subgraph.clone()));
                    current.data.update_filtered_integrated_uv_chain_state(
                        graph,
                        &parent.data,
                        &settings.uv,
                    );

                    current.data.topo_order = i;
                    if settings.uv.generate_integrated {
                        current.data.compute_integrated(
                            graph,
                            vakint,
                            &parent.data,
                            &settings.uv,
                        )?;
                    }
                    if settings.uv.only_integrated {
                        continue;
                    }
                    if use_expanded_4d_local_uv {
                        current.data.compute_expanded_4d(
                            graph,
                            cut_data,
                            &parent.data,
                            valid_orientations,
                            settings,
                            representation,
                        )?;
                    } else {
                        if !matches!(parent.data.local_3d, CFFapprox::Dependent { .. }) {
                            return Err(eyre!(
                                "UV forest parent {parent_id:?} local 3D approximation was not computed before node {n:?}"
                            ));
                        }
                        current.data.compute(
                            graph,
                            cut_data,
                            &parent.data,
                            valid_orientations,
                            &settings.uv,
                            explicit_orientation_sum_only,
                        )?;
                    }
                }
                _ => {
                    return Err(eyre!(
                        "UV forest union nodes are not supported in local counterterm generation"
                    ));
                }
            }
        }
        Ok(())
    }

    fn cff_explicit_sum_has_uv_leading_duplicate_denominators(
        &self,
        graph: &Graph,
        settings: &GenerationSettings,
        representation: ThreeDRepresentation,
    ) -> bool {
        if !settings.uv.subtract_uv
            || !settings.explicit_orientation_sum_only
            || representation != ThreeDRepresentation::Cff
        {
            return false;
        }

        self.dag.nodes.values().any(|node| {
            if node.parents.len() != 1 {
                return false;
            }
            let Some(parent) = self.dag.nodes.get(node.parents[0]) else {
                return false;
            };
            uv_reduced_subgraph_has_duplicate_leading_denominators(graph, &node.data, &parent.data)
        })
    }

    #[instrument(skip_all)]
    pub(crate) fn orientation_parametric_expr(
        &self,
        graph: &Graph,
        add_sigma: bool,
    ) -> Result<Vec<Atom>> {
        let mut sum = None;
        let mut node_dump_terms = std::env::var("GAMMALOOP_DUMP_UV_FOREST_NODE_TERMS")
            .ok()
            .map(|dump_dir| (dump_dir, Vec::<(String, Vec<Atom>)>::new()));
        let node_filter = std::env::var("GAMMALOOP_UV_FOREST_NODE_FILTER").ok();

        for (node_id, n) in &self.dag.nodes {
            let simple_approx = n.data.simple_approx.as_ref().ok_or_else(|| {
                eyre!("UV forest node {node_id:?} simple approximation was not computed")
            })?;
            let simple_expr = simple_approx.expr(&graph.full_filter());
            let node_header = format!(
                "node {node_id:?}; topo_order {}; dod {}; subgraph {}; simple {}",
                n.data.topo_order,
                n.data.dod(),
                n.data.spinney.subgraph.string_label(),
                simple_expr.to_canonical_string()
            );
            if node_filter
                .as_ref()
                .is_some_and(|filter| !node_header.contains(filter))
            {
                continue;
            }
            debug!(dod = %n.data.dod(), simple = %simple_expr, "Terms");
            let mut node_terms = Vec::new();

            let first = sum.is_none();
            let sum = sum.get_or_insert_with(Vec::new);
            let final_integrands = n.data.final_integrand.as_ref().ok_or_else(|| {
                eyre!("UV forest node {node_id:?} final integrand was not computed")
            })?;
            if !first && final_integrands.len() != sum.len() {
                return Err(eyre!(
                    "UV forest node {node_id:?} produced {} residue integrands, but previous nodes produced {}",
                    final_integrands.len(),
                    sum.len()
                ));
            }

            for (i, integrand) in final_integrands.iter().enumerate() {
                let a = if add_sigma {
                    debug!("{}", simple_expr);
                    integrand * function!(GS.if_sigma, simple_expr.clone())
                } else {
                    integrand.clone()
                };
                if first {
                    sum.push(a.clone());
                } else {
                    sum[i] += a.clone();
                }
                node_terms.push(a);
            }

            if let Some((_, dump_terms)) = node_dump_terms.as_mut() {
                dump_terms.push((node_header, node_terms));
            }
        }

        let mut sum = if let Some(sum) = sum {
            sum
        } else if node_filter.is_some() {
            vec![Atom::Zero]
        } else {
            return Err(eyre!("No terms in forest"));
        };

        Self::finalize_parametric_terms(graph, &mut sum);

        if let Some((dump_dir, mut dump_terms)) = node_dump_terms {
            static FOREST_DUMP_COUNTER: std::sync::atomic::AtomicUsize =
                std::sync::atomic::AtomicUsize::new(0);
            let dump_id = FOREST_DUMP_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let mode = std::env::var("GAMMALOOP_DUMP_UV_FOREST_NODE_MODE")
                .unwrap_or_else(|_| "unknown".to_string());
            let path = std::path::Path::new(&dump_dir)
                .join(format!("{}_{}_forest_{dump_id:03}.txt", graph.name, mode));
            let mut dump = String::new();
            use std::fmt::Write;
            writeln!(dump, "graph {}", graph.name).ok();
            writeln!(dump, "mode {mode}").ok();
            writeln!(dump, "forest {dump_id}").ok();
            for (header, terms) in &mut dump_terms {
                Self::finalize_parametric_terms(graph, terms);
                writeln!(dump, "{header}").ok();
                for (term_id, term) in terms.iter().enumerate() {
                    writeln!(dump, "term {term_id}: {}", term.to_canonical_string()).ok();
                }
            }
            writeln!(dump, "sum").ok();
            for (term_id, term) in sum.iter().enumerate() {
                writeln!(dump, "term {term_id}: {}", term.to_canonical_string()).ok();
            }
            std::fs::create_dir_all(&dump_dir).ok();
            std::fs::write(path, dump).ok();
        }

        Ok(sum)
    }

    fn finalize_parametric_terms(graph: &Graph, terms: &mut [Atom]) {
        for (pair, edge_index, _) in graph.iter_edges_of(
            &graph
                .full_filter()
                .subtract(&graph.initial_state_cut)
                .subtract(&graph.tree_edges),
        ) {
            if matches!(pair, HedgePair::Unpaired { .. }) {
                continue;
            }

            for s in terms.iter_mut() {
                *s = s.replace_multiple(&[GS.split_mom_pattern_simple(edge_index)]);
            }
        }

        for s in terms.iter_mut() {
            *s = s
                .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
                .with(W_.d_)
                .collect_factors();
        }
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
            explicitly_summed_orientations: false,
        };

        let zeroed = integrands.zero_like();

        assert_eq!(zeroed.integrands, vec![Atom::Zero, Atom::Zero]);
        assert_eq!(zeroed.cuts, CutSet::empty(3));
    }
}
