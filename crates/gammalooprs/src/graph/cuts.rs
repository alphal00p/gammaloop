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
    /// LTD-only effective selected-variable signs keyed by generated Cutkosky
    /// E-surface id. They encode the sign of the canonical E-surface variable
    /// relative to the positive-energy Cutkosky direction of the corresponding
    /// cut edges. For simple LU cuts, the simultaneous-to-dual residue parity
    /// and resolved LU-basis orientation can be pushed onto this single
    /// selected denominator variable and are included here.
    pub ltd_lu_cut_esurface_signs: Vec<(EsurfaceID, i64)>,
    /// LTD-only residue prefactor mapping branch-local LTD residues to
    /// GammaLoop's resolved Cutkosky-basis residue convention when the bridge
    /// cannot be folded into the selected variable signs. For simple LU cuts
    /// this is one because the bridge is carried by `ltd_lu_cut_esurface_signs`.
    /// For confluent LU cuts the repeated-pole LTD construction already carries
    /// the Cauchy derivative convention, so this prefactor carries the
    /// remaining resolved LU-basis and variable-orientation bridge. CFF
    /// residues are already assembled in the GammaLoop convention and ignore
    /// this factor.
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
