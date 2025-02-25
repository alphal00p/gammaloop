use derive_more::{From, Into};
use linnet::half_edge::{hedgevec::HedgeVec, involution::Orientation};
use serde::{Deserialize, Serialize};

use crate::new_cs::CutId;

use super::{surface::HybridSurfaceID, tree::Tree};

#[derive(Debug, Clone, Serialize, Deserialize, From, Into, Hash, PartialEq, Eq)]
pub struct OrientationID(usize);

pub struct OrientationData {
    pub orientation: HedgeVec<Orientation>,
    pub cuts: Vec<CutId>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CutOrientationExpression {
    pub orientation: OrientationID,
    pub left: Tree<HybridSurfaceID>,
    pub right: Tree<HybridSurfaceID>,
}
