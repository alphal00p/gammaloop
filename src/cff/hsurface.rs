use crate::cff::cff_graph::VertexSet;
use crate::cff::esurface::add_external_shifts;
use crate::utils::{
    cut_energy, external_energy_atom_from_index, ose_atom_from_index, FloatLike, F,
};
use bincode_trait_derive::{Decode, Encode};
use bitvec::vec::BitVec;
use derive_more::{From, Into};
use itertools::Itertools;
use linnet::half_edge::involution::{EdgeIndex, EdgeVec};
use linnet::half_edge::subgraph::InternalSubGraph;
use serde::{Deserialize, Serialize};
use symbolica::atom::Atom;
use symbolica::parse;
use typed_index_collections::TiVec;

use super::esurface::ExternalShift;
use super::{esurface::Esurface, surface};

#[derive(From, Into, Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize, Encode, Decode)]
pub struct HsurfaceID(usize);
pub type HsurfaceCollection = TiVec<HsurfaceID, Hsurface>;
pub type HsurfaceCache<T> = TiVec<HsurfaceID, T>;

#[derive(Serialize, Deserialize, Debug, Clone)]
/// H-surface of the supergraph, is most likely E-surface of the amplitude, kind of badly named.
pub struct Hsurface {
    pub positive_energies: Vec<EdgeIndex>,
    pub negative_energies: Vec<EdgeIndex>,
    pub external_shift: ExternalShift,
}

impl PartialEq for Hsurface {
    fn eq(&self, other: &Self) -> bool {
        self.positive_energies == other.positive_energies
            && self.negative_energies == other.negative_energies
    }
}

impl Hsurface {
    pub(crate) fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        let (symbolic_positive_energies, symbolic_negative_energies) =
            [&self.positive_energies, &self.negative_energies]
                .iter()
                .map(|energies| {
                    energies
                        .iter()
                        .map(|i| {
                            if cut_edges.contains(i) {
                                cut_energy(*i)
                            } else {
                                ose_atom_from_index(*i)
                            }
                        })
                        .collect_vec()
                })
                .collect_tuple()
                .unwrap_or_else(|| unreachable!());

        let symbolic_shift = self
            .external_shift
            .iter()
            .fold(Atom::new(), |sum, (i, sign)| {
                Atom::num(*sign) * external_energy_atom_from_index(*i) + &sum
            });

        let symbolic_sum_positive_energies = symbolic_positive_energies
            .iter()
            .fold(Atom::new(), |sum, e| sum + e);

        let symbolic_sum_negative_energies = symbolic_negative_energies
            .iter()
            .fold(Atom::new(), |sum, e| sum + e);

        symbolic_sum_positive_energies - &symbolic_sum_negative_energies + &symbolic_shift
    }
}

impl From<HsurfaceID> for Atom {
    fn from(value: HsurfaceID) -> Self {
        parse!(&format!("H({})", Into::<usize>::into(value)))
    }
}

#[cfg(test)]
mod tests {
    use bitvec::vec::BitVec;
    use linnet::half_edge::{
        builder::HedgeGraphBuilder,
        involution::{EdgeIndex, Orientation},
        subgraph::InternalSubGraph,
        HedgeGraph, NodeIndex,
    };
    use symbolica::atom::{Atom, AtomCore};

    use symbolica::parse;

    use crate::{
        cff::{cff_graph::VertexSet, esurface::Esurface},
        utils::F,
    };

    use super::Hsurface;

    fn dummy_hedge_graph(num_edges: usize) -> HedgeGraph<(), ()> {
        let mut graph = HedgeGraphBuilder::new();
        graph.add_node(());

        for _ in 0..num_edges {
            graph.add_edge(NodeIndex(0), NodeIndex(0), (), Orientation::Default);
        }

        graph.build()
    }

    #[test]
    fn test_to_atom() {
        let external_shift = vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(5), 1)];
        let h_surface = Hsurface {
            positive_energies: vec![EdgeIndex::from(0), EdgeIndex::from(1)],
            negative_energies: vec![EdgeIndex::from(2), EdgeIndex::from(3)],
            external_shift,
        };

        let h_surface_atom = h_surface.to_atom(&[]);
        let expected_atom = parse!("Q(0, cind(0)) + Q(1, cind(0)) - Q(2, cind(0)) - Q(3, cind(0)) - P(4, cind(0)) + P(5, cind(0))");
        let diff = h_surface_atom - &expected_atom;
        let diff = diff.expand();
        assert_eq!(diff, Atom::new());
    }
}
