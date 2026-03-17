use bincode_trait_derive::{Decode, Encode};
use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};

use crate::cff::esurface::RaisedEsurfaceGroup;

#[derive(Debug, Clone, Encode, Decode)]
pub struct CutSet {
    pub esurfaces: Option<RaisedEsurfaceGroup>,
    pub union: SuBitGraph,
}

impl CutSet {
    pub fn empty(size: usize) -> Self {
        CutSet {
            esurfaces: None,
            union: SuBitGraph::empty(size),
        }
    }
}
