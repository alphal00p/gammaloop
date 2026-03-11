use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};

use crate::cff::esurface::{EsurfaceID, RaisedEsurfaceGroup};

#[derive(Debug, Clone)]
pub struct CutSet {
    pub esurfaces: RaisedEsurfaceGroup,
    pub union: SuBitGraph,
}

impl CutSet {
    pub fn empty(size: usize) -> Self {
        CutSet {
            esurfaces: RaisedEsurfaceGroup {
                esurface_ids: vec![],
                max_occurence: 0,
            },
            union: SuBitGraph::empty(size),
        }
    }
}
