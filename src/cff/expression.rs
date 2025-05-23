use std::borrow::Borrow;

use crate::utils::ose_atom_from_index;
use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use linnet::half_edge::{
    hedgevec::HedgeVec,
    involution::Orientation,
    nodestore::{NodeStorage, NodeStorageOps},
    GVEdgeAttrs, HedgeGraph,
};
use serde::{Deserialize, Serialize};
use spenso::structure::concrete_index::FlatIndex;
use std::fmt::Write;
use symbolica::{
    atom::Atom,
    function,
    id::{Pattern, Replacement},
    parse,
};
use typed_index_collections::TiVec;

use super::{
    cut_expression::SuperGraphOrientationID, generation::SurfaceCache, surface::HybridSurfaceID,
    tree::Tree,
};

#[derive(
    Debug, Clone, Serialize, Deserialize, From, Into, Hash, PartialEq, Eq, Copy, Encode, Decode,
)]
pub struct AmplitudeOrientationID(pub usize);

#[derive(
    Debug, Clone, Serialize, Deserialize, From, Into, Hash, PartialEq, Eq, Copy, Encode, Decode,
)]
pub struct SubgraphOrientationID(pub usize);

pub trait OrientationID: From<usize> + Into<usize> {}

impl OrientationID for AmplitudeOrientationID {}
impl OrientationID for SubgraphOrientationID {}
impl OrientationID for SuperGraphOrientationID {}

#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct OrientationData {
    pub orientation: HedgeVec<Orientation>,
}

impl OrientationData {
    pub fn dot<E, V, N: NodeStorageOps<NodeData = V>>(&self, graph: &HedgeGraph<E, V, N>) {
        let mut writer = String::new();
        writer.push_str("digraph {{");

        writer.push_str(&format!(
            "  node [shape=circle,height=0.1,label=\"\"];  overlap=\"scale\"; layout=\"neato\";",
        ));

        for (hedge_pair, id, _) in graph.iter_all_edges() {
            let attr = GVEdgeAttrs {
                color: None,
                label: None,
                other: None,
            };
            writer.push_str("  ");

            let attr = hedge_pair.fill_color(attr);
            hedge_pair
                .dot_fmt(&mut writer, graph, self.orientation[id], attr)
                .unwrap();
        }
        writer.push_str("}}");
    }

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
    pub expression: Tree<HybridSurfaceID>,
}
#[derive(Clone, Debug, Serialize, Deserialize, Encode, Decode)]
pub struct CFFExpression<O: OrientationID> {
    pub orientations: TiVec<O, OrientationExpression>,
    pub surfaces: SurfaceCache,
}

impl<O: OrientationID> CFFExpression<O>
where
    usize: From<O>,
{
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

    pub fn get_orientation_atoms(&self) -> TiVec<AmplitudeOrientationID, Atom> {
        self.orientations
            .iter()
            .map(|orientation| orientation.expression.to_atom_inv())
            .collect()
    }

    pub fn get_orientation_atoms_with_data(
        &self,
    ) -> TiVec<AmplitudeOrientationID, (Atom, OrientationData)> {
        self.orientations
            .iter()
            .map(|orientation| {
                let atom = orientation.expression.to_atom_inv();
                let data = orientation.data.clone();
                (atom, data)
            })
            .collect()
    }

    pub fn get_orientation_atom(&self, orientation_id: O) -> Atom {
        self.orientations[orientation_id].expression.to_atom_inv()
    }

    pub fn num_unfolded_terms(&self) -> usize {
        self.orientations
            .iter()
            .map(|o| o.expression.get_bottom_layer().len())
            .sum()
    }
}

impl From<CFFExpression<AmplitudeOrientationID>> for CFFExpression<SubgraphOrientationID> {
    fn from(value: CFFExpression<AmplitudeOrientationID>) -> Self {
        Self {
            orientations: value.orientations.raw.into(),
            surfaces: value.surfaces,
        }
    }
}

impl From<CFFExpression<SubgraphOrientationID>> for CFFExpression<AmplitudeOrientationID> {
    fn from(value: CFFExpression<SubgraphOrientationID>) -> Self {
        Self {
            orientations: value.orientations.raw.into(),
            surfaces: value.surfaces,
        }
    }
}
