use std::{collections::BTreeMap, fmt::Display};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
    subgraph::{SubGraphLike, SubSetLike, SubSetOps},
};
use serde::{Deserialize, Serialize};
use symbolica::atom::{Atom, AtomCore};
use three_dimensional_reps::RepresentationMode;

use crate::{
    cff::{
        expression::{
            GammaLoopGraphOrientation, GammaLoopOrientationExpression, GammaLoopThreeDExpression,
            OrientationID, OrientationSelector, ResidualDenominator, ThreeDExpression,
        },
        surface::GammaLoopSurfaceCache,
    },
    graph::{FeynmanGraph, Graph, cuts::CutSet, get_cff_inverse_energy_product_impl},
    settings::global::OrientationPattern,
    utils::GS,
    uv::UltravioletGraph,
};
use color_eyre::Result;

//pub mod cut_expression;
pub mod esurface;
pub mod expression;
pub mod generation;
pub mod hsurface;
pub mod surface;
pub mod tree;
mod vertex_set;
pub(crate) use vertex_set::VertexSet;

pub struct CFFTerm {
    // One per orientation
    pub expression: Vec<Atom>,
    pub orientations: Vec<EdgeVec<Orientation>>,
}

impl CFFTerm {
    pub fn expression_with_selectors(&self) -> Atom {
        let mut result = Atom::Zero;
        for (expr, orient) in self.expression.iter().zip(self.orientations.iter()) {
            result += expr.clone() * orient.orientation_thetas_gs();
        }
        result
    }
}

#[derive(
    Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Encode, Decode, Serialize, Deserialize,
)]
// This describes the combinations of residues that are selected.
pub struct CutCFFIndex {
    pub left_threshold_order: Option<usize>,
    pub right_threshold_order: Option<usize>,
    pub lu_cut_order: Option<usize>,
}

impl CutCFFIndex {
    pub fn new_all_none() -> Self {
        Self {
            left_threshold_order: None,
            right_threshold_order: None,
            lu_cut_order: None,
        }
    }
}

impl Display for CutCFFIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut parts = vec![];
        if let Some(order) = self.lu_cut_order {
            parts.push(format!("lu_cut_{}", order));
        }

        if let Some(order) = self.left_threshold_order {
            parts.push(format!("left_th_{}", order));
        }

        if let Some(order) = self.right_threshold_order {
            parts.push(format!("right_th_{}", order));
        }

        if parts.is_empty() {
            write!(f, "")
        } else {
            write!(f, "{}", parts.join("_"))
        }
    }
}

pub struct CutCFF {
    pub terms: BTreeMap<CutCFFIndex, CFFTerm>,
}

impl CutCFF {
    pub fn expression_with_selectors(&self) -> BTreeMap<CutCFFIndex, Atom> {
        self.terms
            .iter()
            .map(|(index, term)| (*index, term.expression_with_selectors()))
            .collect()
    }
}

#[derive(Clone, Copy)]
enum CutCffResidueAxis {
    RightThreshold,
    LeftThreshold,
    LuCut,
}

impl CutCffResidueAxis {
    fn set_order(self, index: &mut CutCFFIndex, order: usize) {
        match self {
            Self::RightThreshold => index.right_threshold_order = Some(order),
            Self::LeftThreshold => index.left_threshold_order = Some(order),
            Self::LuCut => index.lu_cut_order = Some(order),
        }
    }
}

fn apply_indexed_residue_selection<F>(
    residues: Vec<(CutCFFIndex, ThreeDExpression<OrientationID>)>,
    axis: CutCffResidueAxis,
    mut select: F,
) -> Vec<(CutCFFIndex, ThreeDExpression<OrientationID>)>
where
    F: FnMut(ThreeDExpression<OrientationID>) -> Vec<ThreeDExpression<OrientationID>>,
{
    residues
        .into_iter()
        .flat_map(|(index, expression)| {
            select(expression)
                .into_iter()
                .enumerate()
                .map(move |(i, residue)| {
                    let mut new_index = index;
                    axis.set_order(&mut new_index, i + 1);
                    (new_index, residue)
                })
        })
        .collect()
}

impl Graph {
    fn residual_denominator_factor_gs(
        &self,
        residual_denominators: &[ResidualDenominator],
        wrap_tree_denoms: bool,
    ) -> Atom {
        if residual_denominators.is_empty() {
            return Atom::num(1);
        }

        let denominator = residual_denominators
            .iter()
            .map(|residual| {
                let edge_subgraph = self.get_edge_subgraph(residual.edge_id);
                self.denominator(&edge_subgraph, |_| -(residual.power as isize))
            })
            .reduce(|acc, atom| acc * atom)
            .unwrap_or_else(|| Atom::num(1));

        if wrap_tree_denoms {
            GS.wrap_tree_denoms(denominator)
        } else {
            denominator
        }
    }

    pub(crate) fn three_d_expression_parametric_atom_with_numerator_gs(
        &self,
        expression: &ThreeDExpression<OrientationID>,
        numerator: &Atom,
        representation: RepresentationMode,
        explicit_orientation_sum_only: bool,
        pattern: &OrientationPattern,
    ) -> Atom {
        let atom = if explicit_orientation_sum_only {
            expression.diagnostic_parametric_atom_with_numerator_gs(
                self,
                numerator,
                &OrientationPattern::default(),
            )
        } else {
            expression.parametric_atom_with_numerator_gs(self, numerator, pattern)
        };
        let atom =
            atom * self.residual_denominator_factor_gs(&expression.residual_denominators, true);

        if representation == RepresentationMode::Cff
            && !cff_expression_uses_local_half_edges(expression)
        {
            let graph_without_is_cut = self
                .underlying
                .full_filter()
                .subtract(&self.initial_state_cut.left)
                .subtract(&self.initial_state_cut.right);
            let mut excluded_edges =
                self.preserved_4d_denominator_edges_for_3d_expression(RepresentationMode::Cff);
            excluded_edges.sort_unstable();
            excluded_edges.dedup();
            atom * get_cff_inverse_energy_product_impl(self, &graph_without_is_cut, &excluded_edges)
        } else {
            atom
        }
    }

    pub fn cff<S: SubGraphLike + SubSetLike>(
        &mut self,
        contract_subgraph: &S,
        cutset: &CutSet,
        orientation_pattern: &OrientationPattern,
    ) -> Result<CutCFF> {
        let mut contract_edges = vec![];

        for (p, eid, _) in self.iter_edges_of(contract_subgraph) {
            if p.is_paired() {
                contract_edges.push(eid);
            }
        }
        contract_edges.sort_unstable();
        contract_edges.dedup();

        if self.get_loop_number() == 0 {
            let graph_without_is_cut = self
                .underlying
                .full_filter()
                .subtract(&self.initial_state_cut.left)
                .subtract(&self.initial_state_cut.right);
            let residual_denominators = self
                .underlying
                .iter_edges_of(&graph_without_is_cut)
                .filter_map(|(pair, edge_id, _)| {
                    (pair.is_paired() && !contract_edges.contains(&edge_id))
                        .then(|| ResidualDenominator::new(edge_id, None))
                })
                .collect::<Vec<_>>();
            let orientation = self
                .underlying
                .new_edgevec(|_, _, _| Orientation::Undirected);
            let (expression, orientations) = if orientation_pattern.filter_orientation(&orientation)
            {
                (
                    vec![self.residual_denominator_factor_gs(&residual_denominators, true)],
                    vec![orientation],
                )
            } else {
                (Vec::new(), Vec::new())
            };
            return Ok(CutCFF {
                terms: BTreeMap::from([(
                    CutCFFIndex::new_all_none(),
                    CFFTerm {
                        expression,
                        orientations,
                    },
                )]),
            });
        }

        let canonize_esurface = self.get_esurface_canonization(&self.loop_momentum_basis);

        let mut inverse_energy_excluded_edges = contract_edges.clone();
        inverse_energy_excluded_edges
            .extend(self.preserved_4d_denominator_edges_for_3d_expression(RepresentationMode::Cff));
        inverse_energy_excluded_edges.sort_unstable();
        inverse_energy_excluded_edges.dedup();

        let cff_options = self.denominator_only_cff_3d_expression_options();
        let repeated_active_cff_source =
            self.cff_source_has_repeated_active_denominators(&contract_edges, &cff_options)?;
        let pure_repeated_threshold_residue =
            cutset.residue_selector.is_pure_raised_threshold_residue()
                && repeated_active_cff_source;
        let use_confluent_cff = (cutset.residue_selector.right_th_cut.is_none()
            && cutset.residue_selector.left_th_cut.is_none()
            && cutset.residue_selector.lu_cut().is_none())
            || pure_repeated_threshold_residue;
        let cff = self.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonize_esurface,
            &cff_options,
            use_confluent_cff,
        )?;

        let mut residues = vec![(CutCFFIndex::new_all_none(), cff)];

        if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
            residues = apply_indexed_residue_selection(
                residues,
                CutCffResidueAxis::RightThreshold,
                |expression| {
                    if pure_repeated_threshold_residue
                        && cutset.residue_selector.is_pure_right_threshold_residue()
                    {
                        // The confluent source is already expressed in the
                        // generated repeated-channel coordinate for this
                        // single raised threshold. Re-reading the same
                        // repeated pole in the canonical family would put the
                        // Laurent coefficient back on the unresolved
                        // threshold surface.
                        expression.select_esurface_residue_in_generated_basis(right_threshold)
                    } else {
                        expression.select_esurface_residue(right_threshold)
                    }
                },
            );
        }

        if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
            residues = apply_indexed_residue_selection(
                residues,
                CutCffResidueAxis::LeftThreshold,
                |expression| {
                    if pure_repeated_threshold_residue
                        && cutset.residue_selector.is_pure_left_threshold_residue()
                    {
                        // See the right-threshold branch above. This is the
                        // amplitude-like raised-threshold case after the
                        // CutCFFIndex migration made the selected axis
                        // explicit instead of storing it as a synthetic LU
                        // residue.
                        expression.select_esurface_residue_in_generated_basis(left_threshold)
                    } else {
                        expression.select_esurface_residue(left_threshold)
                    }
                },
            );
        }

        if let Some(lu_cut) = cutset.residue_selector.lu_cut() {
            residues =
                apply_indexed_residue_selection(residues, CutCffResidueAxis::LuCut, |expression| {
                    if cutset.residue_selector.is_threshold_esurface_residue() {
                        if use_confluent_cff {
                            // The confluent source is already expressed in the
                            // generated repeated-channel coordinate. Applying
                            // the canonical selected-denominator sign again
                            // would over-rotate odd repeated threshold poles.
                            expression.select_esurface_residue_in_generated_basis(lu_cut)
                        } else {
                            expression.select_esurface_residue(lu_cut)
                        }
                    } else {
                        expression.select_esurface_residue_with_cut_edges(
                            lu_cut,
                            cutset.residue_selector.lu_cut_edge_sets(),
                        )
                    }
                });
        }

        // println!("residue orders: {}", residue.len());

        let graph_without_is_cut = self
            .underlying
            .full_filter()
            .subtract(&self.initial_state_cut.left)
            .subtract(&self.initial_state_cut.right);

        // The CFF carries the measure normalization for the loop variables that remain after
        // contracting a UV subgraph. Fully contracted integrated CTs therefore get no extra CFF
        // measure factor, while ordinary root terms get the full graph-loop factor.
        let cff_loop_number = self
            .get_loop_number()
            .saturating_sub(self.cyclotomatic_number(contract_subgraph));
        let cff_phase = (-Atom::i()).pow(cff_loop_number as i64);
        let cff_normalization = cff_phase / (Atom::var(GS.pi) * 2).pow(3 * cff_loop_number as i64);
        crate::debug_tags!(#cff, #trace;
            stage = "graph_cff_normalization",
            graph = %self.name,
            cff_loop_number = cff_loop_number,
            log.cff_normalization = cff_normalization,
            "Graph CFF normalization"
        );

        let mut terms = BTreeMap::new();

        let replacement_rules = if cutset.canonicalize_external_shifts {
            self.surface_cache
                .get_all_replacements_gs_in_lmb(&[], &self.loop_momentum_basis)
        } else {
            self.surface_cache.get_all_replacements_gs(&[])
        };

        for (cut_cff_index, expr) in residues {
            let mut cff_term = CFFTerm {
                expression: vec![],
                orientations: vec![],
            };
            for orientation in expr.orientations.iter().filter(|orientation| {
                orientation_pattern.filter_orientation(&orientation.data.orientation)
            }) {
                let eta_expr = orientation.to_atom_gs();
                let mut ose_expr = eta_expr.replace_multiple(&replacement_rules);
                ose_expr *= self.residual_denominator_factor_gs(&expr.residual_denominators, true);
                ose_expr *= get_cff_inverse_energy_product_impl(
                    self,
                    &graph_without_is_cut,
                    &inverse_energy_excluded_edges,
                );

                ose_expr *= cff_normalization.clone();

                crate::debug_tags!(#cff, #trace;
                    stage = "graph_cff_term_expr",
                    graph = %self.name,
                    cut_index = ?cut_cff_index,
                    log.expr = ose_expr,
                    "Graph CFF term expression"
                );
                // println!("ose expr :{}", ose_expr);
                cff_term.expression.push(ose_expr);
                cff_term
                    .orientations
                    .push(orientation.data.orientation.clone());
            }
            terms.insert(cut_cff_index, cff_term);
        }

        let cut_cff = CutCFF { terms };
        Ok(cut_cff)
    }
}

pub(crate) fn cff_expression_uses_local_half_edges(
    expression: &ThreeDExpression<OrientationID>,
) -> bool {
    expression.orientations.iter().any(|orientation| {
        orientation
            .variants
            .iter()
            .any(|variant| !variant.half_edges.is_empty())
    })
}
