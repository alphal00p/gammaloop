use linnet::half_edge::subgraph::SuBitGraph;

use crate::cff::esurface::EsurfaceID;

pub struct CutSet {
    pub esurfaces: Vec<EsurfaceID>,
    pub union: SuBitGraph,
}
