use std::{
    collections::{BTreeMap, btree_map},
    fmt::Display,
};

use bincode_trait_derive::{Decode, Encode};
use eyre::eyre;
use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
    subgraph::{SubGraphLike, SubSetLike, SubSetOps},
};
use symbolica::atom::{Atom, AtomCore};

use crate::{
    GammaLoopContext,
    cff::expression::GraphOrientation,
    graph::{FeynmanGraph, Graph, cuts::CutSet, get_cff_inverse_energy_product_impl},
    settings::global::OrientationPattern,
    utils::GS,
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

#[derive(Clone, Debug, PartialEq, Eq, Hash, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ResidueSelectedTerms {
    terms: BTreeMap<CutCFFIndex, Atom>,
}

impl ResidueSelectedTerms {
    pub fn new(terms: BTreeMap<CutCFFIndex, Atom>) -> Self {
        Self { terms }
    }

    pub fn empty() -> Self {
        Self::new(BTreeMap::new())
    }

    pub fn all_none(atom: Atom) -> Self {
        Self::new(BTreeMap::from([(CutCFFIndex::new_all_none(), atom)]))
    }

    pub fn zero_for_keys(allowed_keys: &[CutCFFIndex]) -> Self {
        allowed_keys.iter().map(|&key| (key, Atom::Zero)).collect()
    }

    pub fn zero_like(&self) -> Self {
        self.terms.keys().map(|&key| (key, Atom::Zero)).collect()
    }

    pub fn map<F: FnMut(Atom) -> Atom>(self, mut map: F) -> Self {
        self.terms
            .into_iter()
            .map(|(index, atom)| (index, map(atom)))
            .collect()
    }

    pub fn try_map<F>(&self, mut map: F) -> Result<Self>
    where
        F: FnMut(CutCFFIndex, &Atom) -> Result<Atom>,
    {
        self.terms
            .iter()
            .map(|(&index, atom)| Ok((index, map(index, atom)?)))
            .collect()
    }

    pub fn try_zip_with<F>(&self, other: &Self, mut combine: F) -> Result<Self>
    where
        F: FnMut(CutCFFIndex, &Atom, &Atom) -> Result<Atom>,
    {
        if self.len() != other.len() {
            return Err(eyre!(
                "Mismatched residue-selected term counts: {} vs {}",
                self.len(),
                other.len()
            ));
        }

        self.terms
            .iter()
            .map(|(&index, left)| {
                let right = other.get(index)?;
                Ok((index, combine(index, left, right)?))
            })
            .collect()
    }

    pub fn multiply_assign_aligned(&mut self, other: Self) -> Result<()> {
        for (index, atom) in other {
            let entry = self.get_mut(index)?;
            *entry *= atom;
        }
        Ok(())
    }

    pub fn add_assign_term(&mut self, index: CutCFFIndex, term: Atom) {
        self.terms
            .entry(index)
            .and_modify(|current| *current += &term)
            .or_insert(term);
    }

    pub fn len(&self) -> usize {
        self.terms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.terms.is_empty()
    }

    pub fn iter(&self) -> btree_map::Iter<'_, CutCFFIndex, Atom> {
        self.terms.iter()
    }

    pub fn iter_mut(&mut self) -> btree_map::IterMut<'_, CutCFFIndex, Atom> {
        self.terms.iter_mut()
    }

    pub fn keys(&self) -> btree_map::Keys<'_, CutCFFIndex, Atom> {
        self.terms.keys()
    }

    pub fn values_mut(&mut self) -> btree_map::ValuesMut<'_, CutCFFIndex, Atom> {
        self.terms.values_mut()
    }

    pub fn values(&self) -> btree_map::Values<'_, CutCFFIndex, Atom> {
        self.terms.values()
    }

    pub fn pop_first(&mut self) -> Option<(CutCFFIndex, Atom)> {
        self.terms.pop_first()
    }

    pub fn get(&self, index: CutCFFIndex) -> Result<&Atom> {
        self.terms
            .get(&index)
            .ok_or_else(|| eyre!("Missing residue-selected term key {:?}", index))
    }

    pub fn get_mut(&mut self, index: CutCFFIndex) -> Result<&mut Atom> {
        self.terms
            .get_mut(&index)
            .ok_or_else(|| eyre!("Missing residue-selected term key {:?}", index))
    }

    pub fn into_only_all_none(mut self) -> Result<Atom> {
        if self.len() != 1 {
            return Err(eyre!(
                "Expected exactly one residue-selected term, got {}",
                self.len()
            ));
        }
        let scalar_key = CutCFFIndex::new_all_none();
        self.terms
            .remove(&scalar_key)
            .ok_or_else(|| eyre!("Missing residue-selected term key {:?}", scalar_key))
    }

    pub fn into_inner(self) -> BTreeMap<CutCFFIndex, Atom> {
        self.terms
    }
}

impl FromIterator<(CutCFFIndex, Atom)> for ResidueSelectedTerms {
    fn from_iter<T: IntoIterator<Item = (CutCFFIndex, Atom)>>(iter: T) -> Self {
        Self::new(iter.into_iter().collect())
    }
}

impl IntoIterator for ResidueSelectedTerms {
    type Item = (CutCFFIndex, Atom);
    type IntoIter = btree_map::IntoIter<CutCFFIndex, Atom>;

    fn into_iter(self) -> Self::IntoIter {
        self.terms.into_iter()
    }
}

impl<'a> IntoIterator for &'a ResidueSelectedTerms {
    type Item = (&'a CutCFFIndex, &'a Atom);
    type IntoIter = btree_map::Iter<'a, CutCFFIndex, Atom>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a> IntoIterator for &'a mut ResidueSelectedTerms {
    type Item = (&'a CutCFFIndex, &'a mut Atom);
    type IntoIter = btree_map::IterMut<'a, CutCFFIndex, Atom>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

#[cfg(test)]
mod tests {
    use super::{CutCFFIndex, ResidueSelectedTerms};
    use eyre::WrapErr;
    use std::collections::BTreeMap;
    use symbolica::{atom::Atom, symbol};

    fn left_threshold_key() -> CutCFFIndex {
        CutCFFIndex {
            left_threshold_order: Some(1),
            right_threshold_order: None,
            lu_cut_order: None,
        }
    }

    #[test]
    fn residue_selected_terms_zero_like_preserves_keys() {
        let terms = ResidueSelectedTerms::new(BTreeMap::from([
            (CutCFFIndex::new_all_none(), Atom::var(symbol!("x"))),
            (left_threshold_key(), Atom::num(3)),
        ]));

        assert_eq!(
            terms.zero_like(),
            ResidueSelectedTerms::new(BTreeMap::from([
                (CutCFFIndex::new_all_none(), Atom::Zero),
                (left_threshold_key(), Atom::Zero),
            ]))
        );
    }

    #[test]
    fn residue_selected_terms_zip_combines_matching_keys() {
        let left = ResidueSelectedTerms::all_none(Atom::var(symbol!("x")));
        let right = ResidueSelectedTerms::all_none(Atom::num(2));

        let combined = left
            .try_zip_with(&right, |_, lhs, rhs| Ok(lhs.clone() + rhs.clone()))
            .unwrap();

        assert_eq!(
            combined.get(CutCFFIndex::new_all_none()).unwrap(),
            &(Atom::var(symbol!("x")) + Atom::num(2))
        );
    }

    #[test]
    fn residue_selected_terms_zip_reports_missing_key_with_context() {
        let left = ResidueSelectedTerms::all_none(Atom::var(symbol!("x")));
        let right =
            ResidueSelectedTerms::new(BTreeMap::from([(left_threshold_key(), Atom::num(2))]));

        let err = left
            .try_zip_with(&right, |_, lhs, rhs| Ok(lhs.clone() + rhs.clone()))
            .wrap_err("while aligning residue-selected terms for test zip")
            .unwrap_err();
        let report = format!("{err:?}");

        assert!(report.contains("while aligning residue-selected terms for test zip"));
        assert!(report.contains("Missing residue-selected term key"));
    }
}

pub struct CutCFF {
    pub terms: BTreeMap<CutCFFIndex, CFFTerm>,
}

impl CutCFF {
    pub fn expression_with_selectors(&self) -> ResidueSelectedTerms {
        self.terms
            .iter()
            .map(|(index, term)| (*index, term.expression_with_selectors()))
            .collect()
    }
}

impl Graph {
    pub fn cff<S: SubGraphLike + SubSetLike>(
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
                .get_all_replacements_in_lmb(&[], &self.loop_momentum_basis)
        } else {
            self.surface_cache.get_all_replacements(&[])
        };

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
