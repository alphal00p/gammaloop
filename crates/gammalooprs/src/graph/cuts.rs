use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};

use crate::cff::esurface::EsurfaceID;

#[derive(Debug, Clone)]
pub struct CutSet {
    pub esurfaces: Vec<EsurfaceID>,
    pub union: SuBitGraph,
}

impl CutSet {
    pub fn empty(size: usize) -> Self {
        CutSet {
            esurfaces: vec![],
            union: SuBitGraph::empty(size),
        }
    }
}
