use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{SuBitGraph, SubSetLike},
};

use crate::cff::{esurface::RaisedEsurfaceGroup, surface::EsurfaceID};

#[derive(Debug, Clone, Encode, Decode, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub struct CutSet {
    pub residue_selector: ResidueSelector,
    pub union: SuBitGraph,
}

#[derive(Debug, Clone, Encode, Decode, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub struct ResidueSelector {
    pub lu_cut: Option<RaisedEsurfaceGroup>,
    /// Original denominator-edge support for each Cutkosky cut represented by
    /// `lu_cut`. A 3D-expression lower-sector/contact variant contributes to a
    /// Cutkosky residue only if it still contains all denominator edges of at
    /// least one of these alternatives.
    pub lu_cut_edge_sets: Vec<Vec<EdgeIndex>>,
    /// LTD-only selected E-surface variables keyed by generated Cutkosky
    /// E-surface id.
    ///
    /// For simple LU cuts this contains the full bridge from the branch-local
    /// LTD selected variable to GammaLoop's positive-energy simultaneous
    /// Cutkosky convention. For repeated/confluent LU cuts the generated
    /// denominator variable is kept canonical and the remaining bridge is
    /// carried by `ltd_lu_cut_residue_prefactor_sign`.
    pub ltd_lu_cut_esurface_signs: Vec<(EsurfaceID, i64)>,
    /// LTD-only residue-basis bridge that cannot be represented as a sign of
    /// one selected denominator variable. CFF residues are already assembled in
    /// GammaLoop's Cutkosky convention and ignore this factor.
    pub ltd_lu_cut_residue_prefactor_sign: i64,
    pub left_th_cut: Option<RaisedEsurfaceGroup>,
    pub right_th_cut: Option<RaisedEsurfaceGroup>,
}

impl CutSet {
    pub fn empty(size: usize) -> Self {
        CutSet {
            residue_selector: ResidueSelector {
                lu_cut: None,
                lu_cut_edge_sets: Vec::new(),
                ltd_lu_cut_esurface_signs: Vec::new(),
                ltd_lu_cut_residue_prefactor_sign: 1,
                left_th_cut: None,
                right_th_cut: None,
            },
            union: SuBitGraph::empty(size),
        }
    }
}

impl ResidueSelector {
    pub(crate) fn ltd_threshold_residue_prefactor_sign(&self) -> i64 {
        // The threshold counterterm radial variable is oriented oppositely to
        // the canonical selected E-surface variable used by the LTD expression.
        // Each selected threshold residue therefore contributes one parity
        // conversion. LU residues use their own Cutkosky-oriented bridge above.
        let threshold_residue_count =
            usize::from(self.left_th_cut.is_some()) + usize::from(self.right_th_cut.is_some());
        if threshold_residue_count.is_multiple_of(2) {
            1
        } else {
            -1
        }
    }

    pub(crate) fn ltd_residue_prefactor_sign(&self) -> i64 {
        self.ltd_lu_cut_residue_prefactor_sign * self.ltd_threshold_residue_prefactor_sign()
    }
}
