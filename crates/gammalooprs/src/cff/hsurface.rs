use crate::utils::{cut_energy, external_energy_atom_from_index, ose_atom_from_index};
use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use symbolica::atom::Atom;

use crate::graph::LoopMomentumBasis;
use crate::momentum::SignOrZero;

use feynkit_cff::HSurface;

/// GammaLoop numerical and symbolic operations on a FeynKit H-surface.
pub trait HSurfaceExt {
    fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom;
    fn to_atom_in_lmb(&self, cut_edges: &[EdgeIndex], lmb: &LoopMomentumBasis) -> Atom;
    fn to_atom_impl(
        &self,
        cut_edges: &[EdgeIndex],
        external_shift_atom: impl Fn(EdgeIndex) -> Atom,
    ) -> Atom;
}

impl HSurfaceExt for HSurface {
    fn to_atom(&self, cut_edges: &[EdgeIndex]) -> Atom {
        self.to_atom_impl(cut_edges, external_energy_atom_from_index)
    }

    fn to_atom_in_lmb(&self, cut_edges: &[EdgeIndex], lmb: &LoopMomentumBasis) -> Atom {
        self.to_atom_impl(cut_edges, |edge| {
            lmb.edge_signatures[edge].external.iter_enumerated().fold(
                Atom::Zero,
                |sum, (external_index, sign)| {
                    let atom = external_energy_atom_from_index(lmb.ext_edges[external_index]);
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
}

#[cfg(test)]
mod tests {

    use linnet::half_edge::involution::EdgeIndex;
    use linnet::half_edge::subgraph::{SuBitGraph, SubSetLike};
    use symbolica::atom::{Atom, AtomCore};

    use symbolica::parse;

    use crate::graph::LoopMomentumBasis;
    use crate::momentum::signature::LoopExtSignature;
    use crate::utils::{external_energy_atom_from_index, test_utils::dummy_hedge_graph};
    use feynkit_cff::VertexSet;

    use super::{HSurface, HSurfaceExt};

    #[test]
    fn to_atom_in_lmb_uses_canonical_external_edges_not_carrier_edges() {
        let dummy_graph = dummy_hedge_graph(9);
        let mut edge_signatures = dummy_graph
            .new_edgevec_from_iter(
                (0..9).map(|_| LoopExtSignature::from((Vec::<isize>::new(), vec![0, 0]))),
            )
            .unwrap();
        edge_signatures[EdgeIndex::from(8)] =
            LoopExtSignature::from((Vec::<isize>::new(), vec![0, 1]));
        let lmb = LoopMomentumBasis {
            tree: SuBitGraph::empty(0),
            loop_edges: vec![].into(),
            ext_edges: vec![EdgeIndex::from(2), EdgeIndex::from(6)].into(),
            edge_signatures,
        };
        let hsurface = HSurface {
            positive_energies: vec![],
            negative_energies: vec![],
            external_shift: vec![(EdgeIndex::from(8), -1)].into(),
            vertex_set: VertexSet::dummy(),
        };

        let atom = hsurface.to_atom_in_lmb(&[], &lmb).expand();
        let expected =
            (Atom::num(-1) * external_energy_atom_from_index(EdgeIndex::from(6))).expand();

        assert_eq!(atom.to_canonical_string(), expected.to_canonical_string());
    }

    mod failing {
        use super::*;

        #[test]
        fn test_to_atom() {
            let external_shift = vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(5), 1)].into();
            let h_surface = HSurface {
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
