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
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};

use linnet::half_edge::{involution::HedgePair, subgraph::SubSetOps};
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

    pub fn sum_orientations_explicitly(self, orientations: &[EdgeVec<Orientation>]) -> Self {
        self.map(|atom| explicit_orientation_sum_atom(&atom, orientations))
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
            uses_expanded_4d_forest_path(graph, cut_data, settings, representation);
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

    #[instrument(skip_all)]
    pub(crate) fn orientation_parametric_expr(
        &self,
        graph: &Graph,
        add_sigma: bool,
    ) -> Result<Vec<Atom>> {
        let mut sum = None;

        for (node_id, n) in &self.dag.nodes {
            let simple_approx = n.data.simple_approx.as_ref().ok_or_else(|| {
                eyre!("UV forest node {node_id:?} simple approximation was not computed")
            })?;
            let simple_expr = simple_approx.expr(&graph.full_filter());
            debug!(dod = %n.data.dod(), simple = %simple_expr, "Terms");

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
