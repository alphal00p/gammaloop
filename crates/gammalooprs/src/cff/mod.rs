use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
    subgraph::{SubSetLike, SubSetOps},
};
use symbolica::atom::{Atom, AtomCore};
use three_dimensional_reps::RepresentationMode;

use crate::{
    cff::{
        expression::{
            GammaLoopGraphOrientation, GammaLoopOrientationExpression, ResidualDenominator,
        },
        surface::GammaLoopSurfaceCache,
    },
    graph::{FeynmanGraph, Graph, cuts::CutSet, get_cff_inverse_energy_product_impl},
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

pub struct CutCFF {
    pub terms: Vec<CFFTerm>, //index is the power of the esurface +1
}

impl CutCFF {
    pub fn expression_with_selectors(&self) -> Vec<Atom> {
        self.terms
            .iter()
            .map(|t| t.expression_with_selectors())
            .collect()
    }
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

    pub fn cff<S: SubSetLike>(&mut self, contract_subgraph: &S, cutset: &CutSet) -> Result<CutCFF> {
        let canonize_esurface = self.get_esurface_canonization(&self.loop_momentum_basis);
        let mut contract_edges = vec![];

        for (p, eid, _) in self.iter_edges_of(contract_subgraph) {
            if p.is_paired() {
                contract_edges.push(eid);
            }
        }
        contract_edges.sort_unstable();
        contract_edges.dedup();

        let mut inverse_energy_excluded_edges = contract_edges.clone();
        inverse_energy_excluded_edges
            .extend(self.preserved_4d_denominator_edges_for_3d_expression(RepresentationMode::Cff));
        inverse_energy_excluded_edges.sort_unstable();
        inverse_energy_excluded_edges.dedup();

        let cff_options = self.denominator_only_cff_3d_expression_options();
        let cff = self.generate_3d_expression_for_integrand(
            &contract_edges,
            &canonize_esurface,
            &cff_options,
        )?;
        let residue = if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
            cff.select_esurface_residue(right_threshold).pop().unwrap()
        } else {
            cff
        };

        let residue = if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
            residue
                .select_esurface_residue(left_threshold)
                .pop()
                .unwrap()
        } else {
            residue
        };

        let residue = if let Some(lu_cut) = cutset.residue_selector.lu_cut.as_ref() {
            residue.select_esurface_residue(lu_cut)
        } else {
            vec![residue]
        };

        // println!("residue orders: {}", residue.len());

        let graph_without_is_cut = self
            .underlying
            .full_filter()
            .subtract(&self.initial_state_cut.left)
            .subtract(&self.initial_state_cut.right);
        let mut terms = vec![];

        let replacement_rules = self.surface_cache.get_all_replacements_gs(&[]);

        for expr in residue.into_iter() {
            let mut cff_term = CFFTerm {
                expression: vec![],
                orientations: vec![],
            };
            for orientation in expr.orientations.iter() {
                let eta_expr = orientation.to_atom_gs();
                let mut ose_expr = eta_expr.replace_multiple(&replacement_rules);
                ose_expr *= self.residual_denominator_factor_gs(&expr.residual_denominators, true);

                let inverse_energies = get_cff_inverse_energy_product_impl(
                    self,
                    &graph_without_is_cut,
                    &inverse_energy_excluded_edges,
                );

                ose_expr *= inverse_energies;

                // println!("ose expr :{}", ose_expr);
                cff_term.expression.push(ose_expr);
                cff_term
                    .orientations
                    .push(orientation.data.orientation.clone());
            }
            terms.push(cff_term);
        }

        let cut_cff = CutCFF { terms };
        Ok(cut_cff)
    }
}
