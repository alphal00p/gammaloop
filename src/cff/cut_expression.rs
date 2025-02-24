use linnet::half_edge::{hedgevec::HedgeVec, involution::Orientation};
use serde::{Deserialize, Serialize};

use super::{surface::HybridSurfaceID, tree::Tree};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CutOrientationExpression {
    pub orientation: HedgeVec<Orientation>,
    pub left: Tree<HybridSurfaceID>,
    pub right: Tree<HybridSurfaceID>,
}
