use derive_more::{From, Into};
use linnet::half_edge::{hedgevec::HedgeVec, involution::Orientation};
use serde::{Deserialize, Serialize};
use symbolica::atom::Atom;
use typed_index_collections::TiVec;

use crate::new_cs::CutId;

use super::{generation::SurfaceCache, surface::HybridSurfaceID, tree::Tree};

#[derive(Debug, Clone, Serialize, Deserialize, From, Into, Hash, PartialEq, Eq, Copy)]
pub struct OrientationID(usize);

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CutOrientationData {
    pub orientation: HedgeVec<Orientation>,
    pub cuts: Vec<CutId>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CutOrientationExpression {
    pub data: CutOrientationData,
    pub expressions: Vec<SingleCutOrientationExpression>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CFFCutExpression {
    pub orientations: TiVec<OrientationID, CutOrientationExpression>,
    pub surfaces: SurfaceCache,
}

impl CFFCutExpression {
    pub fn new_empty() -> Self {
        Self {
            orientations: TiVec::new(),
            surfaces: SurfaceCache {
                esurface_cache: TiVec::new(),
                hsurface_cache: TiVec::new(),
            },
        }
    }

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

    pub fn get_orientation_atom(&self, orientation_id: OrientationID) -> Vec<Atom> {
        self.orientations[orientation_id]
            .expressions
            .iter()
            .map(|expr| Atom::from(expr))
            .collect()
    }

    pub fn get_orientation_atoms(&self) -> TiVec<OrientationID, Vec<Atom>> {
        self.orientations
            .iter()
            .map(|orientation| {
                let atoms = orientation
                    .expressions
                    .iter()
                    .map(|expr| Atom::from(expr))
                    .collect();
                atoms
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use core::panic;

    use itertools::Itertools;
    use linnet::half_edge::{
        builder::HedgeGraphBuilder,
        involution::{Flow, Orientation},
        nodestorage::NodeStorageVec,
    };
    use typed_index_collections::TiVec;

    use crate::{
        cff, disable,
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

        let cuts: TiVec<CutId, CrossSectionCut> = hedge_graph
            .all_cuts(
                hedge_graph[&nodes[0]].clone(),
                hedge_graph[&nodes[3]].clone(),
            )
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

        let cut_expression = cff::generation::generate_cff_with_cuts(&hedge_graph, &None, &cuts);
        let atom_cut_1 = cut_expression.to_atom_for_cut(CutId::from(2));
        let atom_with_energies = cut_expression.surfaces.substitute_energies(&atom_cut_1);

        println!("{}", atom_cut_1);
        println!("{}", atom_with_energies);
    }
}
