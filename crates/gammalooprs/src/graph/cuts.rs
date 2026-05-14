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
    /// LTD-only residue bridge keyed by generated Cutkosky E-surface id.
    /// GammaLoop cross-section LU cuts are simultaneous Cutkosky residues,
    /// while LTD selects the same support from branch-local dual residues.
    /// For a cut with n on-shell propagators the bridge contains the generic
    /// simultaneous-residue parity (-1)^(n-1), multiplied by the product of
    /// cut-edge orientation signs that converts the graph source-side edge
    /// flow to the left-to-right positive-energy Cutkosky convention. Keeping
    /// this on the E-surface lets raised Cutkosky groups containing several cut
    /// alternatives retain branch-local signs before the E-surfaces are
    /// normalized to their representative. CFF residues are already assembled
    /// in the GammaLoop Cutkosky convention and ignore this.
    pub ltd_lu_cut_esurface_signs: Vec<(EsurfaceID, i64)>,
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
                left_th_cut: None,
                right_th_cut: None,
            },
            union: SuBitGraph::empty(size),
        }
    }
}
