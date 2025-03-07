use derive_more::{From, Into};
use linnet::half_edge::{hedgevec::HedgeVec, involution::Orientation};
use serde::{Deserialize, Serialize};
use symbolica::atom::Atom;
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

impl From<&CutOrientationExpression> for Atom {
    fn from(value: &CutOrientationExpression) -> Self {
        let left_atom = value.left.to_atom_inv();
        let right_atom = value.right.to_atom_inv();
        left_atom * right_atom
    }
}

pub struct OrientationExpression {
    pub data: OrientationData,
    pub expressions: Vec<CutOrientationExpression>,
}

pub struct CFFCutExpression {
    pub orientations: TiVec<OrientationID, OrientationExpression>,
    pub surfaces: SurfaceCache,
}

impl CFFCutExpression {
    pub fn to_atom_for_cut(&self, cut: CutId) -> Atom {
        self.orientations
            .iter()
            .filter_map(|orientation| {
                orientation
                    .data
                    .cuts
                    .iter()
                    .position(|c| *c == cut)
                    .map(|cut_orientation_index| {
                        Atom::from(&orientation.expressions[cut_orientation_index])
                    })
            })
            .fold(Atom::new(), |acc, x| acc + x)
    }
}

#[cfg(test)]
mod tests {
    use itertools::Itertools;
    use linnet::half_edge::{builder::HedgeGraphBuilder, involution::Orientation};

    #[test]
    fn test_double_triangle() {
        let mut hedge_graph_builder = HedgeGraphBuilder::new();

        let nodes = (0..4)
            .map(|_| hedge_graph_builder.add_node(()))
            .collect_vec();

        hedge_graph_builder.add_edge(nodes[0], nodes[1], (), Orientation::Default);
        todo!()
    }
}
