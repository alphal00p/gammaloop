use crate::cff::cff_graph::VertexSet;
use crate::utils::{cut_energy, external_energy_atom_from_index, ose_atom_from_index};
use bincode_trait_derive::{Decode, Encode};

use derive_more::{From, Into};
use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use serde::{Deserialize, Serialize};
use symbolica::atom::Atom;
use symbolica::parse;
use tracing::warn;
use typed_index_collections::TiVec;

use crate::graph::LoopMomentumBasis;
use crate::momentum::SignOrZero;

use super::esurface::Esurface;
use super::esurface::ExternalShift;

#[derive(
    From, Into, Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize, Encode, Decode, Hash,
)]
pub struct HsurfaceID(usize);
pub type HsurfaceCollection = TiVec<HsurfaceID, Hsurface>;
pub type HsurfaceCache<T> = TiVec<HsurfaceID, T>;

#[derive(Serialize, Deserialize, Debug, Clone)]
/// H-surface of the supergraph, is most likely E-surface of the amplitude, kind of badly named.
pub struct Hsurface {
    pub positive_energies: Vec<EdgeIndex>,
    pub negative_energies: Vec<EdgeIndex>,
    pub external_shift: ExternalShift,
    pub vertex_set: VertexSet,
}

impl PartialEq for Hsurface {
    fn eq(&self, other: &Self) -> bool {
        self.positive_energies == other.positive_energies
            && self.negative_energies == other.negative_energies
    }
}

impl Hsurface {
    pub(crate) fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        self.to_atom_impl(cut_edges, external_energy_atom_from_index)
    }

    pub(crate) fn to_atom_in_lmb(&self, cut_edges: &[EdgeIndex], lmb: &LoopMomentumBasis) -> Atom {
        self.to_atom_impl(cut_edges, |edge| {
            lmb.edge_signatures[edge].external.iter_enumerated().fold(
                Atom::Zero,
                |sum, (external_index, sign)| {
                    let atom = external_energy_atom_from_index(EdgeIndex::from(usize::from(
                        external_index,
                    )));
                    match sign {
                        SignOrZero::Zero => sum,
                        SignOrZero::Plus => sum + atom,
                        SignOrZero::Minus => sum - atom,
                    }
                },
            )
        })
    }

    fn to_atom_impl(
        &self,
        cut_edges: &[EdgeIndex],
        external_shift_atom: impl Fn(EdgeIndex) -> Atom,
    ) -> Atom {
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
                Atom::num(*sign) * external_shift_atom(*i) + &sum
            });

        let symbolic_sum_positive_energies = symbolic_positive_energies
            .iter()
            .fold(Atom::new(), |sum, e| sum + e);

        let symbolic_sum_negative_energies = symbolic_negative_energies
            .iter()
            .fold(Atom::new(), |sum, e| sum + e);

        symbolic_sum_positive_energies - &symbolic_sum_negative_energies + &symbolic_shift
    }

    #[allow(dead_code)]
    pub(crate) fn equality_under_energy_conservation(
        &self,
        other: &Esurface,
        constraints: &[&Esurface],
    ) -> Option<bool> {
        if !self.external_shift.is_empty() {
            warn!("this is not handled yet");
            return None;
        }
        constraints
            .iter()
            .find(|esurface| {
                self.negative_energies
                    .iter()
                    .all(|index| esurface.energies.contains(index))
            })
            .map(|constraint| {
                let energies_to_be_added = constraint
                    .energies
                    .iter()
                    .filter(|index| !self.negative_energies.contains(index));

                let mut new_positive_energies = self.positive_energies.clone();
                new_positive_energies.extend(energies_to_be_added);
                new_positive_energies.sort();

                let new_esurface = Esurface {
                    energies: new_positive_energies,
                    external_shift: constraint.external_shift.clone(),
                    vertex_set: VertexSet::dummy(),
                };

                other == &new_esurface
            })
    }

    pub fn equality_by_try_convert(&self, other: &Esurface) -> bool {
        let negative_as_external_shift = self
            .negative_energies
            .iter()
            .map(|e| (*e, -1))
            .collect_vec();

        let self_as_esurface = Esurface {
            energies: self.positive_energies.clone(),
            external_shift: negative_as_external_shift,
            vertex_set: VertexSet::dummy(),
        };

        self_as_esurface == *other
    }
}

impl From<HsurfaceID> for Atom {
    fn from(value: HsurfaceID) -> Self {
        parse!(&format!("H({})", Into::<usize>::into(value)))
    }
}

#[cfg(test)]
mod tests {

    use linnet::half_edge::involution::EdgeIndex;
    use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};
    use symbolica::atom::{Atom, AtomCore};

    use symbolica::parse;

    use crate::cff::{cff_graph::VertexSet, esurface::Esurface};
    use crate::graph::LoopMomentumBasis;
    use crate::momentum::signature::LoopExtSignature;
    use crate::utils::{external_energy_atom_from_index, test_utils::dummy_hedge_graph};

    use super::Hsurface;

    #[test]
    fn test_equality_under_energy_conservation() {
        let constraint = Esurface {
            energies: vec![EdgeIndex::from(1), EdgeIndex::from(2)],
            external_shift: vec![(EdgeIndex::from(0), -1)],
            vertex_set: VertexSet::dummy(),
        };

        let hsurface = Hsurface {
            positive_energies: vec![EdgeIndex::from(3), EdgeIndex::from(5)],
            negative_energies: vec![EdgeIndex::from(2)],
            external_shift: vec![],
            vertex_set: VertexSet::dummy(),
        };

        let other = Esurface {
            energies: vec![EdgeIndex::from(1), EdgeIndex::from(3), EdgeIndex::from(5)],
            external_shift: vec![(EdgeIndex::from(0), -1)],
            vertex_set: VertexSet::dummy(),
        };

        assert!(
            hsurface
                .equality_under_energy_conservation(&other, &[&constraint])
                .unwrap()
        );
    }

    #[test]
    fn to_atom_in_lmb_uses_external_slots_not_carrier_edges() {
        let dummy_graph = dummy_hedge_graph(7);
        let mut edge_signatures = dummy_graph
            .new_edgevec_from_iter(
                (0..7).map(|_| LoopExtSignature::from((Vec::<isize>::new(), vec![0]))),
            )
            .unwrap();
        edge_signatures[EdgeIndex::from(6)] =
            LoopExtSignature::from((Vec::<isize>::new(), vec![1]));
        let lmb = LoopMomentumBasis {
            tree: SuBitGraph::empty(0),
            loop_edges: vec![].into(),
            ext_edges: vec![EdgeIndex::from(6)].into(),
            edge_signatures,
        };
        let hsurface = Hsurface {
            positive_energies: vec![],
            negative_energies: vec![],
            external_shift: vec![(EdgeIndex::from(6), -1)],
            vertex_set: VertexSet::dummy(),
        };

        let atom = hsurface.to_atom_in_lmb(&[], &lmb).expand();
        let expected =
            (Atom::num(-1) * external_energy_atom_from_index(EdgeIndex::from(0))).expand();

        assert_eq!(atom.to_canonical_string(), expected.to_canonical_string());
    }

    mod failing {
        use super::*;

        #[test]
        fn test_to_atom() {
            let external_shift = vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(5), 1)];
            let h_surface = Hsurface {
                positive_energies: vec![EdgeIndex::from(0), EdgeIndex::from(1)],
                negative_energies: vec![EdgeIndex::from(2), EdgeIndex::from(3)],
                external_shift,
                vertex_set: VertexSet::dummy(),
            };

            let h_surface_atom = h_surface.to_atom(&[]);
            let expected_atom = parse!(
                "Q(0, cind(0)) + Q(1, cind(0)) - Q(2, cind(0)) - Q(3, cind(0)) - P(4, cind(0)) + P(5, cind(0))"
            );
            let diff = h_surface_atom - &expected_atom;
            let diff = diff.expand();
            assert_eq!(diff, Atom::new());
        }
    }
}
