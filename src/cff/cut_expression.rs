use ahash::HashMap;
use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use linnet::half_edge::{hedgevec::HedgeVec, involution::Orientation};
use serde::{Deserialize, Serialize};
use std::ops::Index;
use symbolica::atom::Atom;
use typed_index_collections::TiVec;

use crate::new_cs::CutId;

use super::{
    expression::{AmplitudeOrientationID, CFFExpression},
    generation::SurfaceCache,
    surface::HybridSurfaceID,
    tree::Tree,
};

#[derive(
    Debug, Clone, Serialize, Deserialize, From, Into, Hash, PartialEq, Eq, Copy, Encode, Decode,
)]
pub struct SuperGraphOrientationID(pub usize);

#[derive(Debug, Clone, Serialize, Deserialize, Encode)]
pub struct CutOrientationData {
    pub orientation: HedgeVec<Orientation>,
    pub cuts: Vec<CutId>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode)]
pub struct SingleCutOrientationExpression {
    pub left: Tree<HybridSurfaceID>,
    pub right: Tree<HybridSurfaceID>,
}

impl From<&SingleCutOrientationExpression> for Atom {
    fn from(value: &SingleCutOrientationExpression) -> Self {
        let left_atom = value.left.to_atom_inv();
        let right_atom = value.right.to_atom_inv();
        left_atom * right_atom
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode)]
pub struct OrientationMap {
    map: HashMap<SuperGraphOrientationID, (AmplitudeOrientationID, AmplitudeOrientationID)>,
    revesed_map: HashMap<(AmplitudeOrientationID, AmplitudeOrientationID), SuperGraphOrientationID>,
}

impl Index<SuperGraphOrientationID> for OrientationMap {
    type Output = (AmplitudeOrientationID, AmplitudeOrientationID);

    fn index(&self, index: SuperGraphOrientationID) -> &Self::Output {
        &self.map[&index]
    }
}

impl Index<(AmplitudeOrientationID, AmplitudeOrientationID)> for OrientationMap {
    type Output = SuperGraphOrientationID;

    fn index(&self, index: (AmplitudeOrientationID, AmplitudeOrientationID)) -> &Self::Output {
        &self.revesed_map[&index]
    }
}

impl OrientationMap {
    pub fn new() -> Self {
        Self {
            map: HashMap::default(),
            revesed_map: HashMap::default(),
        }
    }

    pub fn insert(
        &mut self,
        orientation_id: SuperGraphOrientationID,
        left: AmplitudeOrientationID,
        right: AmplitudeOrientationID,
    ) {
        self.map.insert(orientation_id, (left, right));
        self.revesed_map.insert((left, right), orientation_id);
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode)]
pub struct SingleCutExpression {
    pub left_amplitude: CFFExpression,
    pub right_amplitude: CFFExpression,
    pub orientation_map: OrientationMap,
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode)]
pub struct CFFCutsExpression {
    pub cut_expressions: TiVec<CutId, SingleCutExpression>,
    pub surfaces: SurfaceCache,
    pub orientation_data: TiVec<SuperGraphOrientationID, CutOrientationData>,
}

impl CFFCutsExpression {
    pub fn new_empty() -> Self {
        Self {
            cut_expressions: TiVec::new(),
            surfaces: SurfaceCache {
                esurface_cache: TiVec::new(),
                hsurface_cache: TiVec::new(),
            },
            orientation_data: TiVec::new(),
        }
    }

    pub fn to_atom_for_cut(&self, cut: CutId) -> Atom {
        let left_amplitude_atom = self.cut_expressions[cut].left_amplitude.to_atom();

        let right_amplitude_atom = self.cut_expressions[cut].right_amplitude.to_atom();

        left_amplitude_atom * right_amplitude_atom
    }

    pub fn get_orientation_atom(&self, orientation_id: SuperGraphOrientationID) -> Vec<Atom> {
        self.orientation_data[orientation_id]
            .cuts
            .iter()
            .map(|cut| {
                let cut_expression = &self.cut_expressions[*cut];
                let (left_orientation_id, right_orientation_id) =
                    cut_expression.orientation_map[orientation_id];

                let left_amplitude_atom = cut_expression
                    .left_amplitude
                    .get_orientation_atom(left_orientation_id);
                let right_amplitude_atom = cut_expression
                    .right_amplitude
                    .get_orientation_atom(right_orientation_id);

                left_amplitude_atom * right_amplitude_atom
            })
            .collect()
    }

    pub fn get_orientation_atoms(&self) -> TiVec<SuperGraphOrientationID, Vec<Atom>> {
        self.orientation_data
            .iter_enumerated()
            .map(|(orientation_id, _)| self.get_orientation_atom(orientation_id))
            .collect()
    }
}

// merge orientations, If one of the entries is undirected, the orientation of the other entry is set. If both are undirected, return None,
// if two entries are different, return None
pub fn amplitude_orientations_to_sg_orientaion(
    left: &HedgeVec<Orientation>,
    right: &HedgeVec<Orientation>,
) -> Option<HedgeVec<Orientation>> {
    let mut result = Vec::with_capacity(left.len());

    for ((_, left_entry), (_, right_entry)) in left.into_iter().zip(right.into_iter()) {
        match (left_entry, right_entry) {
            (Orientation::Undirected, Orientation::Undirected) => {
                return None;
            }
            (Orientation::Default, Orientation::Reversed) => {
                return None;
            }
            (Orientation::Reversed, Orientation::Default) => {
                return None;
            }
            (Orientation::Default, Orientation::Default) => {
                result.push(Orientation::Default);
            }
            (Orientation::Reversed, Orientation::Reversed) => {
                result.push(Orientation::Reversed);
            }
            (Orientation::Undirected, _) => {
                result.push(*right_entry);
            }
            (_, Orientation::Undirected) => {
                result.push(*left_entry);
            }
        }
    }

    Some(HedgeVec::from_raw(result))
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use linnet::half_edge::{
        builder::HedgeGraphBuilder,
        involution::{Flow, Orientation},
        nodestore::NodeStorageVec,
    };
    use typed_index_collections::TiVec;

    use crate::{
        cff,
        new_cs::{CrossSectionCut, CutId},
    };

    #[test]
    fn test_double_triangle() {
        let mut hedge_graph_builder = HedgeGraphBuilder::new();

        let nodes = (0..4)
            .map(|_| hedge_graph_builder.add_node(()))
            .collect_vec();

        hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[0], nodes[2], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[1], nodes[2], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[1], nodes[3], (), Orientation::Undirected);
        hedge_graph_builder.add_edge(nodes[2], nodes[3], (), Orientation::Undirected);

        hedge_graph_builder.add_external_edge(nodes[0], (), Orientation::Undirected, Flow::Sink);
        hedge_graph_builder.add_external_edge(nodes[3], (), Orientation::Undirected, Flow::Source);

        let hedge_graph = hedge_graph_builder.build::<NodeStorageVec<_>>();
        let node_0 = hedge_graph.iter_crown(nodes[0]).into();
        let node_3 = hedge_graph.iter_crown(nodes[3]).into();

        let cuts: TiVec<CutId, CrossSectionCut> = hedge_graph
            .all_cuts(node_0, node_3)
            .into_iter()
            .map(|(left, cut, right)| CrossSectionCut { left, cut, right })
            .collect();

        for cut in cuts.iter() {
            let edges_in_cut = hedge_graph
                .iter_edges(&cut.cut)
                .map(|(_, id, _)| id)
                .collect_vec();

            println!("{:?}", edges_in_cut);
            let left_dot = hedge_graph.dot(&cut.left);
            let right_dot = hedge_graph.dot(&cut.right);

            println!("Left: {}", left_dot);
            println!("Right: {}", right_dot);
        }

        let cut_expression =
            cff::generation::generate_cff_with_cuts(&hedge_graph, &None, &cuts).unwrap();
        let atom_cut_1 = cut_expression.to_atom_for_cut(CutId::from(2));
        let atom_with_energies = cut_expression.surfaces.substitute_energies(&atom_cut_1);

        println!("{}", atom_cut_1);
        println!("{}", atom_with_energies);
    }
}
