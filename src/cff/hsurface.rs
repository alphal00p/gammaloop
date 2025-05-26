use crate::cff::esurface::add_external_shifts;
use crate::cff::surface::Surface;
use crate::utils::{external_energy_atom_from_index, ose_atom_from_index, FloatLike, F};
use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use itertools::Itertools;
use linnet::half_edge::hedgevec::HedgeVec;
use linnet::half_edge::involution::EdgeIndex;
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
pub struct Hsurface {
    pub positive_energies: Vec<EdgeIndex>,
    pub negative_energies: Vec<EdgeIndex>,
    pub external_shift: ExternalShift,
}

impl Surface for Hsurface {
    fn get_positive_energies(&self) -> impl Iterator<Item = &EdgeIndex> {
        self.positive_energies.iter()
    }

    fn get_negative_energies(&self) -> impl Iterator<Item = &EdgeIndex> {
        self.negative_energies.iter()
    }

    fn get_external_shift(&self) -> impl Iterator<Item = &(EdgeIndex, i64)> {
        self.external_shift.iter()
    }
}

impl PartialEq for Hsurface {
    fn eq(&self, other: &Self) -> bool {
        self.positive_energies == other.positive_energies
            && self.negative_energies == other.negative_energies
    }
}

impl Hsurface {
    #[inline]
    pub fn compute_value<T: FloatLike>(&self, energy_cache: &HedgeVec<F<T>>) -> F<T> {
        surface::compute_value(self, energy_cache)
    }

    #[inline]
    pub fn compute_shift_part<T: FloatLike>(&self, energy_cache: &HedgeVec<F<T>>) -> F<T> {
        surface::compute_shift_part(self, energy_cache)
    }

    pub fn to_atom(&self) -> Atom {
        let (symbolic_positive_energies, symbolic_negative_energies) =
            [&self.positive_energies, &self.negative_energies]
                .iter()
                .map(|energies| {
                    energies
                        .iter()
                        .map(|i| ose_atom_from_index(*i))
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

    pub fn to_atom_with_rewrite(&self, esurface: &Esurface) -> Option<Atom> {
        let possible = self
            .negative_energies
            .iter()
            .all(|energy| esurface.energies.contains(energy));

        if !possible {
            return None;
        }

        let additional_positive_energies = esurface
            .energies
            .iter()
            .filter(|energy| !self.negative_energies.contains(energy))
            .copied();

        let energies = self
            .positive_energies
            .iter()
            .copied()
            .chain(additional_positive_energies)
            .sorted()
            .collect_vec();

        let external_shift = add_external_shifts(&self.external_shift, &esurface.external_shift);

        let dummy_esurface = Esurface {
            energies,
            external_shift,
        };

        Some(dummy_esurface.to_atom())
    }
}

pub fn compute_hsurface_cache<T: FloatLike>(
    hsurfaces: &HsurfaceCollection,
    energy_cache: &HedgeVec<F<T>>,
) -> HsurfaceCache<F<T>> {
    hsurfaces
        .iter()
        .map(|hsurface| hsurface.compute_value(energy_cache))
        .collect_vec()
        .into()
}

impl From<HsurfaceID> for Atom {
    fn from(value: HsurfaceID) -> Self {
        parse!(&format!("H({})", Into::<usize>::into(value)))
    }
}

#[cfg(test)]
mod tests {
    use linnet::half_edge::{
        builder::HedgeGraphBuilder,
        involution::{EdgeIndex, Orientation},
        HedgeGraph, NodeIndex,
    };
    use symbolica::atom::{Atom, AtomCore};

    use symbolica::parse;

    use crate::{cff::esurface::Esurface, utils::F};

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
    fn test_compute_shift_part() {
        let dummy_graph = dummy_hedge_graph(4);
        let external_shift = vec![(EdgeIndex::from(2), 1), (EdgeIndex::from(3), -1)];

        let h_surface = Hsurface {
            positive_energies: vec![EdgeIndex::from(0)],
            negative_energies: vec![EdgeIndex::from(1)],
            external_shift,
        };

        let energy_cache = dummy_graph
            .new_hedgevec_from_iter(vec![F(1.0), F(2.0), F(3.0), F(4.0)])
            .unwrap();
        let shift_part = h_surface.compute_shift_part(&energy_cache);
        assert_eq!(shift_part.0, -1.0);
    }

    #[test]
    fn test_compute_value() {
        let dummy_graph = dummy_hedge_graph(6);
        let external_shift = vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(5), 1)];

        let h_surface = Hsurface {
            positive_energies: vec![EdgeIndex::from(0), EdgeIndex::from(1)],
            negative_energies: vec![EdgeIndex::from(2), EdgeIndex::from(3)],
            external_shift,
        };

        let energy_cache = dummy_graph
            .new_hedgevec_from_iter(vec![F(1.0), F(2.0), F(3.0), F(4.0), F(5.0), F(6.0)])
            .unwrap();
        let value = h_surface.compute_value(&energy_cache);
        assert_eq!(value.0, 1.0 + 2.0 - 3.0 - 4.0 - 5.0 + 6.0);
    }

    #[test]
    fn test_to_atom() {
        let external_shift = vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(5), 1)];
        let h_surface = Hsurface {
            positive_energies: vec![EdgeIndex::from(0), EdgeIndex::from(1)],
            negative_energies: vec![EdgeIndex::from(2), EdgeIndex::from(3)],
            external_shift,
        };

        let h_surface_atom = h_surface.to_atom();
        let expected_atom = parse!("Q(0, cind(0)) + Q(1, cind(0)) - Q(2, cind(0)) - Q(3, cind(0)) - P(4, cind(0)) + P(5, cind(0))");
        let diff = h_surface_atom - &expected_atom;
        let diff = diff.expand();
        assert_eq!(diff, Atom::new());
    }

    #[test]
    fn test_to_atom_with_rewrite() {
        let external_shift = vec![(EdgeIndex::from(4), -1), (EdgeIndex::from(5), 1)];

        let h_surface = Hsurface {
            positive_energies: vec![EdgeIndex::from(0), EdgeIndex::from(1)],
            negative_energies: vec![EdgeIndex::from(2), EdgeIndex::from(3)],
            external_shift: external_shift.clone(),
        };

        let e_surface = Esurface {
            energies: vec![EdgeIndex::from(2), EdgeIndex::from(3), EdgeIndex::from(6)],
            external_shift: vec![(EdgeIndex::from(4), 1)],
        };

        let h_surface_atom = h_surface.to_atom().expand();
        let e_surface_atom = e_surface.to_atom();

        println!("rewriting {}", h_surface_atom);
        println!("with {}", e_surface_atom);

        let rewritten = h_surface.to_atom_with_rewrite(&e_surface).unwrap();
        let target = parse!("Q(0, cind(0)) + Q(1, cind(0)) + Q(6, cind(0)) + P(5, cind(0))");
        let diff = rewritten - &target;
        let diff = diff.expand();

        assert_eq!(diff, Atom::new(), "diff: {}", diff);

        let h_surface_2 = Hsurface {
            positive_energies: vec![EdgeIndex::from(1), EdgeIndex::from(7)],
            negative_energies: vec![EdgeIndex::from(3), EdgeIndex::from(6)],
            external_shift,
        };

        let rewritten = h_surface_2.to_atom_with_rewrite(&e_surface).unwrap();
        let target = parse!("Q(1, cind(0)) + Q(7, cind(0)) + Q(2, cind(0)) + P(5, cind(0))");
        let diff = rewritten - &target;
        let diff = diff.expand();

        assert_eq!(diff, Atom::new(), "diff: {}", diff);
    }
}
