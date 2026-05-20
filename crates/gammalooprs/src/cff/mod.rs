use std::{collections::BTreeMap, fmt::Display};

use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
    subgraph::{SubSetLike, SubSetOps},
};
use symbolica::atom::{Atom, AtomCore};

use crate::{
    cff::expression::GraphOrientation,
    graph::{FeynmanGraph, Graph, cuts::CutSet, get_cff_inverse_energy_product_impl},
    settings::global::OrientationPattern,
};
use color_eyre::Result;

pub mod cff_graph;
//pub mod cut_expression;
pub mod esurface;
pub mod expression;
pub mod generation;
pub mod hsurface;
pub mod surface;
pub mod tree;

pub struct CFFTerm {
    // One per orientation
    pub expression: Vec<Atom>,
    pub orientations: Vec<EdgeVec<Orientation>>,
}

impl CFFTerm {
    pub fn expression_with_selectors(&self) -> Atom {
        let mut result = Atom::Zero;
        for (expr, orient) in self.expression.iter().zip(self.orientations.iter()) {
            result += expr.clone() * orient.orientation_thetas();
        }
        result
    }
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug, Hash, Encode, Decode)]
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

impl Graph {
    pub fn cff<S: SubSetLike>(
        &mut self,
        contract_subgraph: &S,
        cutset: &CutSet,
        orientation_pattern: &OrientationPattern,
    ) -> Result<CutCFF> {
        let canonize_esurface = self.get_esurface_canonization(&self.loop_momentum_basis);
        let mut contract_edges = vec![];

        for (p, eid, _) in self.iter_edges_of(contract_subgraph) {
            if p.is_paired() {
                contract_edges.push(eid);
            }
        }

        let cff = [(
            CutCFFIndex::new_all_none(),
            self.generate_cff(&contract_edges, &canonize_esurface, orientation_pattern)?,
        )];

        let mut residues = BTreeMap::new();

        cff.into_iter()
            .flat_map(|(index, cff_expression)| {
                if let Some(right_threshold) = cutset.residue_selector.right_th_cut.as_ref() {
                    cff_expression
                        .select_esurface_residue(right_threshold)
                        .into_iter()
                        .enumerate()
                        .map(|(i, residue)| {
                            let mut new_index = index;
                            new_index.right_threshold_order = Some(i + 1);
                            (new_index, residue)
                        })
                        .collect()
                } else {
                    vec![(index, cff_expression)]
                }
            })
            .flat_map(|(index, cff_expression)| {
                if let Some(left_threshold) = cutset.residue_selector.left_th_cut.as_ref() {
                    cff_expression
                        .select_esurface_residue(left_threshold)
                        .into_iter()
                        .enumerate()
                        .map(|(i, residue)| {
                            let mut new_index = index;
                            new_index.left_threshold_order = Some(i + 1);
                            (new_index, residue)
                        })
                        .collect()
                } else {
                    vec![(index, cff_expression)]
                }
            })
            .flat_map(|(index, cff_expression)| {
                if let Some(lu_cut) = cutset.residue_selector.lu_cut.as_ref() {
                    cff_expression
                        .select_esurface_residue(lu_cut)
                        .into_iter()
                        .enumerate()
                        .map(|(i, residue)| {
                            let mut new_index = index;
                            new_index.lu_cut_order = Some(i + 1);
                            (new_index, residue)
                        })
                        .collect()
                } else {
                    vec![(index, cff_expression)]
                }
            })
            .for_each(|(index, residue)| {
                residues.insert(index, residue);
            });

        // println!("residue orders: {}", residue.len());

        let graph_without_is_cut = self
            .underlying
            .full_filter()
            .subtract(&self.initial_state_cut.left)
            .subtract(&self.initial_state_cut.right);

        let mut terms = BTreeMap::new();

        let replacement_rules = self.surface_cache.get_all_replacements(&[]);

        for (cut_cff_index, expr) in residues.into_iter() {
            let mut cff_term = CFFTerm {
                expression: vec![],
                orientations: vec![],
            };
            for orientation in expr.orientations.iter() {
                let eta_expr = orientation.expression.to_atom_inv();
                let mut ose_expr = eta_expr.replace_multiple(&replacement_rules);

                let inverse_energies = get_cff_inverse_energy_product_impl(
                    self,
                    &graph_without_is_cut,
                    &contract_edges,
                );

                ose_expr *= inverse_energies;

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
