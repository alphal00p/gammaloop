use derive_more::{From, Into};
use linnet::half_edge::{hedgevec::HedgeVec, involution::Orientation};
use serde::{Deserialize, Serialize};
use typed_index_collections::TiVec;

use crate::new_cs::CutId;

use super::{generation::SurfaceCache, surface::HybridSurfaceID, tree::Tree};

#[derive(Debug, Clone, Serialize, Deserialize, From, Into, Hash, PartialEq, Eq, Copy)]
pub struct OrientationID(usize);

pub struct OrientationData {
    pub orientation: HedgeVec<Orientation>,
    pub cuts: Vec<CutId>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CutOrientationExpression {
    pub left: Tree<HybridSurfaceID>,
    pub right: Tree<HybridSurfaceID>,
}

pub struct OrientationExpression {
    pub data: OrientationData,
    pub expressions: Vec<CutOrientationExpression>,
}

pub struct CFFCutExpression {
    pub orientations: TiVec<OrientationID, OrientationExpression>,
    pub surfaces: SurfaceCache,
}
