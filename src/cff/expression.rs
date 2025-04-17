use bincode::{Decode, Encode};
use linnet::half_edge::{hedgevec::HedgeVec, involution::Orientation};
use serde::{Deserialize, Serialize};
use symbolica::atom::Atom;
use typed_index_collections::TiVec;

use super::{
    cut_expression::OrientationID, generation::SurfaceCache, surface::HybridSurfaceID, tree::Tree,
};

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct OrientationData {
    #[bincode(with_serde)]
    pub orientation: HedgeVec<Orientation>,
}

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct OrientationExpression {
    pub data: OrientationData,
    #[bincode(with_serde)]
    pub expression: Tree<HybridSurfaceID>,
}
#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct CFFExpression {
    #[bincode(with_serde)]
    pub orientations: TiVec<OrientationID, OrientationExpression>,
    pub surfaces: SurfaceCache,
}

impl CFFExpression {
    pub fn new_empty() -> Self {
        Self {
            orientations: TiVec::new(),
            surfaces: SurfaceCache {
                esurface_cache: TiVec::new(),
                hsurface_cache: TiVec::new(),
            },
        }
    }

    pub fn to_atom(&self) -> Atom {
        self.orientations
            .iter()
            .map(|orientation| orientation.expression.to_atom_inv())
            .reduce(|a, b| a + b)
            .unwrap_or_default()
    }

    pub fn get_orientation_atoms(&self) -> TiVec<OrientationID, Atom> {
        self.orientations
            .iter()
            .map(|orientation| orientation.expression.to_atom_inv())
            .collect()
    }

    pub fn num_unfolded_terms(&self) -> usize {
        self.orientations
            .iter()
            .map(|o| o.expression.get_bottom_layer().len())
            .sum()
    }
}
