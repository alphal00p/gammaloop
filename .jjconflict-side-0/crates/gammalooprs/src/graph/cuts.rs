use bincode_trait_derive::{Decode, Encode};
use gammaloop_tracing_filter::LogMessage;
use itertools::Itertools;
use linnet::half_edge::{
    involution::EdgeIndex,
    subgraph::{SuBitGraph, SubSetLike},
};

use crate::cff::{CutCFFIndex, esurface::RaisedEsurfaceGroup};

#[derive(Debug, Clone, Encode, Decode, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub struct CutSet {
    pub residue_selector: ResidueSelector,
    pub union: SuBitGraph,
    pub canonicalize_external_shifts: bool,
}

impl LogMessage for CutSet {
    fn log_display(&self) -> String {
        format!("CutSet: union={}", self.union.string_label())
    }
}

#[derive(Debug, Clone, Encode, Decode, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub struct LuCutSelection {
    pub raised_group: RaisedEsurfaceGroup,
    /// Normalized physical edge support for each actual Cutkosky cut in this
    /// raised group. Alternative order remains the physical-cut order so a
    /// future LTD policy can attach representation-specific data without
    /// reconstructing alternatives from the lossy union subgraph.
    pub cut_edge_alternatives: Vec<Vec<EdgeIndex>>,
}

#[derive(Debug, Clone, Encode, Decode, PartialEq, Hash, Eq, PartialOrd, Ord)]
pub struct ResidueSelector {
    pub lu: Option<LuCutSelection>,
    pub left_th_cut: Option<RaisedEsurfaceGroup>,
    pub right_th_cut: Option<RaisedEsurfaceGroup>,
}

impl ResidueSelector {
    pub fn generate_allowed_keys(&self) -> Vec<CutCFFIndex> {
        let allowed_keys = vec![CutCFFIndex::new_all_none()];
        allowed_keys
            .into_iter()
            .flat_map(|index| {
                if let Some(lu) = &self.lu {
                    (1..=lu.raised_group.max_occurence)
                        .map(|lu_cut_index| {
                            let mut new_index = index;
                            new_index.lu_cut_order = Some(lu_cut_index);
                            new_index
                        })
                        .collect_vec()
                } else {
                    vec![index]
                }
            })
            .flat_map(|index| {
                if let Some(left_th_cut) = &self.left_th_cut {
                    (1..=left_th_cut.max_occurence)
                        .map(|left_th_cut_index| {
                            let mut new_index = index;
                            new_index.left_threshold_order = Some(left_th_cut_index);
                            new_index
                        })
                        .collect_vec()
                } else {
                    vec![index]
                }
            })
            .flat_map(|index| {
                if let Some(right_th_cut) = &self.right_th_cut {
                    (1..=right_th_cut.max_occurence)
                        .map(|right_th_cut_index| {
                            let mut new_index = index;
                            new_index.right_threshold_order = Some(right_th_cut_index);
                            new_index
                        })
                        .collect_vec()
                } else {
                    vec![index]
                }
            })
            .collect()
    }
}

impl CutSet {
    pub fn empty(size: usize) -> Self {
        CutSet {
            residue_selector: ResidueSelector {
                lu: None,
                left_th_cut: None,
                right_th_cut: None,
            },
            union: SuBitGraph::empty(size),
            canonicalize_external_shifts: false,
        }
    }
}
