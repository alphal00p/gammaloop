use std::borrow::Borrow;

use bincode::{Decode, Encode};
use linnet::half_edge::{hedgevec::HedgeVec, involution::Orientation};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::Atom,
    id::{Pattern, Replacement},
    parse,
};
use typed_index_collections::TiVec;

use crate::utils::ose_atom_from_index;

use super::{
    cut_expression::OrientationID, generation::SurfaceCache, surface::HybridSurfaceID, tree::Tree,
};

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct OrientationData {
    #[bincode(with_serde)]
    pub orientation: HedgeVec<Orientation>,
}

impl OrientationData {
    pub fn get_ose_replacements(&self) -> Vec<Replacement> {
        self.orientation
            .borrow()
            .into_iter()
            .filter_map(|(edge_index, orientation)| {
                if matches!(orientation, Orientation::Reversed) {
                    let energy_atom = ose_atom_from_index(edge_index);
                    let neg_energy_atom = -&energy_atom;

                    Some(Replacement::new(
                        Pattern::from(energy_atom),
                        Pattern::from(neg_energy_atom),
                    ))
                } else {
                    None
                }
            })
            .collect()
    }
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

    pub fn get_orientation_atoms_with_data(&self) -> TiVec<OrientationID, (Atom, OrientationData)> {
        self.orientations
            .iter()
            .map(|orientation| {
                let atom = orientation.expression.to_atom_inv();
                let data = orientation.data.clone();
                (atom, data)
            })
            .collect()
    }

    pub fn num_unfolded_terms(&self) -> usize {
        self.orientations
            .iter()
            .map(|o| o.expression.get_bottom_layer().len())
            .sum()
    }
}
