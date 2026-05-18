use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{SuBitGraph, SubSetLike},
};

use crate::cff::{esurface::RaisedEsurfaceGroup, surface::EsurfaceID};

#[derive(Debug, Clone, Copy, Encode, Decode, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum LuResidueSelectionBasis {
    /// Residues are assembled in GammaLoop's positive-energy Cutkosky
    /// convention, so selected generated E-surface signs are not applied during
    /// combinatorial residue selection.
    PositiveEnergyCutkosky,
    /// Residues are selected in the generated E-surface variables. The
    /// selected-denominator signs stored in the residue selector orient those
    /// variables against the positive-energy Cutkosky convention.
    GeneratedEsurface,
}

#[derive(Debug, Clone, Encode, Decode, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub struct CutSet {
    pub residue_selector: ResidueSelector,
    pub union: SuBitGraph,
}

#[derive(Debug, Clone, Copy, Encode, Decode, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub struct LuLocalSeriesCoordinate {
    pub esurface_id: EsurfaceID,
    pub external_edge: EdgeIndex,
}

#[derive(Debug, Clone, Encode, Decode, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub struct LuResiduePlan {
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
    /// LTD-only selected E-surface variables for direct Laurent extraction.
    ///
    /// These signs orient the generated E-surface parameter used in
    /// `ltd_lu_residue_parameter` against GammaLoop's positive-energy
    /// Cutkosky coordinate. They are distinct from
    /// `ltd_lu_cut_esurface_signs`, which are the signs needed by the
    /// combinatorial residue selector in the generated E-surface basis.
    pub ltd_lu_cut_local_series_esurface_signs: Vec<(EsurfaceID, i64)>,
    /// LTD-only residue-basis bridge that cannot be represented as a sign of
    /// one selected denominator variable. CFF residues are already assembled in
    /// GammaLoop's Cutkosky convention and ignore this factor.
    pub ltd_lu_cut_residue_prefactor_sign: i64,
    /// LTD-only bridge for the direct local-series extraction path. For simple
    /// LU cuts this is paired with
    /// `ltd_lu_cut_local_series_esurface_signs`: the per-surface signs orient
    /// the local Laurent parameter, while this prefactor carries the full
    /// surface-family and projection bridge. Repeated/confluent cuts carry the
    /// orientation of the canonical Laurent parameter here; the repeated-pole
    /// derivative bridge belongs only to the combinatorial residue-selection
    /// path above.
    pub ltd_lu_cut_local_series_prefactor_sign: i64,
    /// LTD-only bridge for the direct original/root integrand extracted with a
    /// Laurent parameter. Simple LU cuts use the generated selected
    /// E-surface orientation; the auxiliary local-series denominator
    /// orientation is not part of the final original-integrand residue.
    /// Repeated/confluent cuts use the ordinary residue bridge.
    pub ltd_lu_cut_direct_original_prefactor_sign: i64,
    /// LTD-only direct local-series coordinates. Each entry identifies the
    /// generated E-surface variable and the external-energy localization
    /// variable solved for when introducing the Laurent parameter. The separate
    /// local-series prefactor bridges that coordinate to GammaLoop's
    /// positive-energy Cutkosky variable.
    pub ltd_lu_cut_local_series_coordinates: Vec<LuLocalSeriesCoordinate>,
}

#[derive(Debug, Clone, Encode, Decode, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub struct ResidueSelector {
    pub lu_plan: Option<LuResiduePlan>,
    pub left_th_cut: Option<RaisedEsurfaceGroup>,
    pub right_th_cut: Option<RaisedEsurfaceGroup>,
}

impl CutSet {
    pub fn empty(size: usize) -> Self {
        CutSet {
            residue_selector: ResidueSelector {
                lu_plan: None,
                left_th_cut: None,
                right_th_cut: None,
            },
            union: SuBitGraph::empty(size),
        }
    }
}

impl LuResiduePlan {
    pub(crate) fn new_lu_cut(
        lu_cut: RaisedEsurfaceGroup,
        lu_cut_edge_sets: Vec<Vec<EdgeIndex>>,
        ltd_lu_cut_esurface_signs: Vec<(EsurfaceID, i64)>,
        ltd_lu_cut_local_series_esurface_signs: Vec<(EsurfaceID, i64)>,
        ltd_lu_cut_residue_prefactor_sign: i64,
        ltd_lu_cut_local_series_prefactor_sign: i64,
        ltd_lu_cut_direct_original_prefactor_sign: i64,
        ltd_lu_cut_local_series_coordinates: Vec<LuLocalSeriesCoordinate>,
    ) -> Self {
        Self {
            lu_cut: Some(lu_cut),
            lu_cut_edge_sets,
            ltd_lu_cut_esurface_signs,
            ltd_lu_cut_local_series_esurface_signs,
            ltd_lu_cut_residue_prefactor_sign,
            ltd_lu_cut_local_series_prefactor_sign,
            ltd_lu_cut_direct_original_prefactor_sign,
            ltd_lu_cut_local_series_coordinates,
        }
    }

    pub(crate) fn new_threshold_esurface(lu_cut: RaisedEsurfaceGroup) -> Self {
        Self {
            lu_cut: Some(lu_cut),
            lu_cut_edge_sets: Vec::new(),
            ltd_lu_cut_esurface_signs: Vec::new(),
            ltd_lu_cut_local_series_esurface_signs: Vec::new(),
            ltd_lu_cut_residue_prefactor_sign: 1,
            ltd_lu_cut_local_series_prefactor_sign: 1,
            ltd_lu_cut_direct_original_prefactor_sign: 1,
            ltd_lu_cut_local_series_coordinates: Vec::new(),
        }
    }

    pub(crate) fn lu_cut(&self) -> Option<&RaisedEsurfaceGroup> {
        self.lu_cut.as_ref()
    }
}

impl ResidueSelector {
    pub(crate) fn new_lu_cut(
        lu_cut: RaisedEsurfaceGroup,
        lu_cut_edge_sets: Vec<Vec<EdgeIndex>>,
        ltd_lu_cut_esurface_signs: Vec<(EsurfaceID, i64)>,
        ltd_lu_cut_local_series_esurface_signs: Vec<(EsurfaceID, i64)>,
        ltd_lu_cut_residue_prefactor_sign: i64,
        ltd_lu_cut_local_series_prefactor_sign: i64,
        ltd_lu_cut_direct_original_prefactor_sign: i64,
        ltd_lu_cut_local_series_coordinates: Vec<LuLocalSeriesCoordinate>,
        left_th_cut: Option<RaisedEsurfaceGroup>,
        right_th_cut: Option<RaisedEsurfaceGroup>,
    ) -> Self {
        Self {
            lu_plan: Some(LuResiduePlan::new_lu_cut(
                lu_cut,
                lu_cut_edge_sets,
                ltd_lu_cut_esurface_signs,
                ltd_lu_cut_local_series_esurface_signs,
                ltd_lu_cut_residue_prefactor_sign,
                ltd_lu_cut_local_series_prefactor_sign,
                ltd_lu_cut_direct_original_prefactor_sign,
                ltd_lu_cut_local_series_coordinates,
            )),
            left_th_cut,
            right_th_cut,
        }
    }

    pub(crate) fn new_threshold_esurface(lu_cut: RaisedEsurfaceGroup) -> Self {
        Self {
            lu_plan: Some(LuResiduePlan::new_threshold_esurface(lu_cut)),
            left_th_cut: None,
            right_th_cut: None,
        }
    }

    pub(crate) fn lu_plan_mut(&mut self) -> Option<&mut LuResiduePlan> {
        self.lu_plan.as_mut()
    }

    pub(crate) fn lu_cut(&self) -> Option<&RaisedEsurfaceGroup> {
        self.lu_plan.as_ref().and_then(LuResiduePlan::lu_cut)
    }

    pub(crate) fn lu_cut_edge_sets(&self) -> &[Vec<EdgeIndex>] {
        self.lu_plan
            .as_ref()
            .map(|plan| plan.lu_cut_edge_sets.as_slice())
            .unwrap_or(&[])
    }

    pub(crate) fn ltd_lu_cut_esurface_signs(&self) -> &[(EsurfaceID, i64)] {
        self.lu_plan
            .as_ref()
            .map(|plan| plan.ltd_lu_cut_esurface_signs.as_slice())
            .unwrap_or(&[])
    }

    pub(crate) fn ltd_lu_cut_local_series_esurface_signs(&self) -> &[(EsurfaceID, i64)] {
        self.lu_plan
            .as_ref()
            .map(|plan| plan.ltd_lu_cut_local_series_esurface_signs.as_slice())
            .unwrap_or(&[])
    }

    pub(crate) fn ltd_lu_cut_local_series_coordinates(&self) -> &[LuLocalSeriesCoordinate] {
        self.lu_plan
            .as_ref()
            .map(|plan| plan.ltd_lu_cut_local_series_coordinates.as_slice())
            .unwrap_or(&[])
    }

    pub(crate) fn is_threshold_esurface_residue(&self) -> bool {
        self.lu_cut().is_some()
            && self.lu_cut_edge_sets().is_empty()
            && self.left_th_cut.is_none()
            && self.right_th_cut.is_none()
    }

    pub(crate) fn has_lu_cut_residue(&self) -> bool {
        self.lu_cut().is_some() && !self.lu_cut_edge_sets().is_empty()
    }

    pub(crate) fn ltd_threshold_residue_prefactor_sign(&self) -> i64 {
        if self.left_th_cut.is_some() || self.right_th_cut.is_some() {
            // Cross-section threshold counterterms are first localized on the
            // threshold E-surface and then on the LU surface. The threshold
            // Cauchy orientation therefore bridges to the same branch-local
            // coordinate as the direct LU Laurent series, canceling that
            // coordinate orientation from the threshold term without
            // reintroducing a graph-dependent rule.
            self.lu_plan
                .as_ref()
                .map(|plan| plan.ltd_lu_cut_local_series_prefactor_sign)
                .unwrap_or(1)
        } else if self.is_threshold_esurface_residue() {
            // Amplitude threshold counterterms store the selected threshold
            // E-surface in `lu_cut` to reuse the raised-surface residue
            // machinery, but they are not Cutkosky/LU cuts. After the LTD
            // energy residue, the remaining threshold denominator has the
            // opposite Cauchy orientation from the CFF threshold-counterterm
            // convention used by the radial helper. This is a representation
            // bridge for threshold E-surface residues; true LU cuts without a
            // threshold residue carry their own Cutkosky orientation data.
            -1
        } else {
            1
        }
    }

    pub(crate) fn ltd_residue_prefactor_sign(&self) -> i64 {
        self.lu_plan
            .as_ref()
            .map(|plan| plan.ltd_lu_cut_residue_prefactor_sign)
            .unwrap_or(1)
            * self.ltd_threshold_residue_prefactor_sign()
    }

    pub(crate) fn generated_esurface_lu_residue_prefactor_sign(&self) -> i64 {
        if self.has_lu_cut_residue() {
            self.lu_plan
                .as_ref()
                .map(|plan| plan.ltd_lu_cut_local_series_prefactor_sign)
                .unwrap_or(1)
        } else {
            1
        }
    }

    pub(crate) fn ltd_local_series_residue_prefactor_sign(&self) -> i64 {
        self.lu_plan
            .as_ref()
            .map(|plan| plan.ltd_lu_cut_local_series_prefactor_sign)
            .unwrap_or(1)
            * self.ltd_threshold_residue_prefactor_sign()
    }

    pub(crate) fn ltd_direct_original_residue_prefactor_sign(&self) -> i64 {
        self.lu_plan
            .as_ref()
            .map(|plan| plan.ltd_lu_cut_direct_original_prefactor_sign)
            .unwrap_or(1)
            * self.ltd_threshold_residue_prefactor_sign()
    }
}
