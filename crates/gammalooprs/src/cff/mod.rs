use std::{collections::BTreeMap, fmt::Display};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::subgraph::{SubGraphLike, SubSetLike, SubSetOps};
use serde::{Deserialize, Serialize};
use symbolica::atom::{Atom, AtomCore};

use crate::{
    cff::{
        expression::{
            GammaLoopOrientationExpression, OrientationExpression, OrientationID,
            OrientationSelector, ThreeDExpression,
        },
        orientations::GraphOrientation,
        surface::GammaLoopSurfaceCache,
    },
    graph::{FeynmanGraph, Graph, cuts::CutSet, get_cff_inverse_energy_product_impl},
    settings::global::OrientationPattern,
    utils::GS,
    uv::Integrands,
};
use color_eyre::Result;
use three_dimensional_reps::Generate3DExpressionOptions;

pub mod orientations;
//pub mod cut_expression;
pub mod esurface;
pub mod expression;
pub mod generation;
pub mod hsurface;
pub mod surface;
pub mod tree;
mod vertex_set;
pub(crate) use vertex_set::VertexSet;

pub(crate) struct CFFOrientationTerm {
    pub(crate) expression: Atom,
    pub(crate) orientation: OrientationExpression,
}

pub struct CFFTerm {
    // One denominator expression per reduced-graph energy map. Keep the map
    // until the UV localizer can assign its full-graph OrientationID.
    pub(crate) orientations: Vec<CFFOrientationTerm>,
}

impl CFFTerm {
    pub fn expression_with_selectors(&self) -> Atom {
        self.orientations
            .iter()
            .map(|term| {
                term.expression.clone() * term.orientation.data.orientation.orientation_thetas()
            })
            .reduce(|left, right| left + right)
            .unwrap_or(Atom::Zero)
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
    pub fn expression_with_selectors(&self) -> Integrands {
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
    pub fn cff<S: SubGraphLike + SubSetLike>(
        &mut self,
        contract_subgraph: &S,
        cutset: &CutSet,
        orientation_pattern: &OrientationPattern,
        options: &Generate3DExpressionOptions,
        analysis_numerator: Option<&Atom>,
    ) -> Result<CutCFF> {
        let mut contract_edges = vec![];

        for (p, eid, _) in self.iter_edges_of(contract_subgraph) {
            if p.is_paired() {
                contract_edges.push(eid);
            }
        }
        contract_edges.sort_unstable();
        contract_edges.dedup();

        let canonize_esurface = self.get_esurface_canonization(&self.loop_momentum_basis);

        // Reduced UV graphs use the same CFF dispatcher and representation family as the
        // production graph. The old
        // repeated-channel branch switched only this side to a generated residue coordinate, so
        // its surviving edge-energy maps could not be exact restrictions of production maps.
        // Repeated thresholds are selected below in the canonical E-surface family; a distinct
        // confluent/LTD representation remains deferred rather than being selected implicitly.
        let cff = self.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonize_esurface,
            options,
            analysis_numerator,
            false,
        )?;

        let mut residues = vec![(CutCFFIndex::new_all_none(), cff)];

        if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
            residues = apply_indexed_residue_selection(
                residues,
                CutCffResidueAxis::RightThreshold,
                |expression| expression.select_esurface_residue(right_threshold),
            );
        }

        if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
            residues = apply_indexed_residue_selection(
                residues,
                CutCffResidueAxis::LeftThreshold,
                |expression| expression.select_esurface_residue(left_threshold),
            );
        }

        if let Some(lu_cut) = cutset.residue_selector.lu_cut.as_ref() {
            residues =
                apply_indexed_residue_selection(residues, CutCffResidueAxis::LuCut, |expression| {
                    // The old branch distinguished threshold E-surface residues here so a
                    // confluent source could remain in its generated repeated-channel coordinate
                    // without reapplying the canonical selected-denominator sign. It also filtered
                    // lower-sector variants against each individual Cutkosky cut. Main's selector
                    // exposes neither that classification nor those edge-set alternatives, and its
                    // union subgraph is not equivalent to them. Retain main's canonical LU
                    // residue semantics until neutral metadata has an agreed owner: residues are
                    // assembled in GammaLoop's positive-energy Cutkosky convention, so no
                    // generated E-surface selection sign is applied here.
                    expression.select_esurface_residue(lu_cut)
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

        for (cut_cff_index, expr) in residues {
            let replacement_rules = if cutset.canonicalize_external_shifts {
                expr.surfaces
                    .get_all_replacements_gs_in_lmb(&[], &self.loop_momentum_basis)
            } else {
                expr.surfaces.get_all_replacements_gs(&[])
            };
            let mut cff_term = CFFTerm {
                orientations: vec![],
            };
            for orientation in expr.orientations.iter().filter(|orientation| {
                orientation_pattern.filter_orientation(&orientation.data.orientation)
            }) {
                let eta_expr = orientation.to_atom_gs();
                let mut ose_expr = eta_expr.replace_multiple(&replacement_rules);
                ose_expr *= get_cff_inverse_energy_product_impl(
                    self,
                    &graph_without_is_cut,
                    &contract_edges,
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
                cff_term.orientations.push(CFFOrientationTerm {
                    expression: ose_expr,
                    orientation: orientation.clone(),
                });
            }
            terms.insert(cut_cff_index, cff_term);
        }

        let cut_cff = CutCFF { terms };
        Ok(cut_cff)
    }
}
