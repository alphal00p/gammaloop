use crate::utils::thermal_distribution_atom_from_index;
use bincode_trait_derive::{Decode, Encode};
use derive_more::{From, Into};
use linnet::{half_edge::involution::EdgeIndex, num_traits::Sign};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use symbolica::{atom::Atom, parse};
use typed_index_collections::TiVec;

#[derive(From, Into, Copy, Clone, Debug, Eq, PartialEq, Serialize, Deserialize, Encode, Decode)]
pub struct ThermalNumeratorID(usize);
pub type ThermalNumeratorCollection = TiVec<ThermalNumeratorID, ThermalNumerator>;
pub type ThermalNumeratorCache<T> = TiVec<ThermalNumeratorID, T>;

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct ThermalNumerator {
    pub positive_energies: Vec<EdgeIndex>,
    pub negative_energies: Vec<EdgeIndex>,
}

impl ThermalNumerator {
    pub(crate) fn from_edge_lists_canonicalized(
        mut positive_energies: Vec<EdgeIndex>,
        mut negative_energies: Vec<EdgeIndex>,
    ) -> (Self, Sign) {
        positive_energies.sort();
        negative_energies.sort();

        let should_swap = match negative_energies.len().cmp(&positive_energies.len()) {
            Ordering::Greater => true,
            Ordering::Less => false,
            Ordering::Equal => positive_energies > negative_energies,
        };

        let (positive_energies, negative_energies, sign) = if should_swap {
            (negative_energies, positive_energies, Sign::Negative)
        } else {
            (positive_energies, negative_energies, Sign::Positive)
        };

        (
            Self {
                positive_energies,
                negative_energies,
            },
            sign,
        )
    }

    pub(crate) fn is_trivial(&self) -> bool {
        self.positive_energies.len() == 1 && self.negative_energies.is_empty()
    }

    pub(crate) fn to_atom(&self, _cut_edges: &[EdgeIndex]) -> Atom {
        let product = |positive_part_sign: Sign, negative_part_sign: Sign| {
            let positive_part = self
                .positive_energies
                .iter()
                .fold(Atom::num(1), |acc, &edge| {
                    acc * thermal_distribution_atom_from_index(edge, positive_part_sign)
                });

            let negative_part = self
                .negative_energies
                .iter()
                .fold(Atom::num(1), |acc, &edge| {
                    acc * thermal_distribution_atom_from_index(edge, negative_part_sign)
                });

            positive_part * negative_part
        };

        product(Sign::Positive, Sign::Negative) - product(Sign::Negative, Sign::Positive)
    }
}

impl From<ThermalNumeratorID> for Atom {
    fn from(id: ThermalNumeratorID) -> Self {
        parse!(&format!("Tnum({})", Into::<usize>::into(id.0)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use linnet::num_traits::Sign;
    use symbolica::atom::AtomCore;

    #[test]
    fn thermal_numerator_to_atom_has_expected_sign_structure() {
        let numerator = ThermalNumerator {
            positive_energies: vec![EdgeIndex::from(1), EdgeIndex::from(3)],
            negative_energies: vec![EdgeIndex::from(2)],
        };

        let expected = thermal_distribution_atom_from_index(EdgeIndex::from(1), Sign::Positive)
            * thermal_distribution_atom_from_index(EdgeIndex::from(3), Sign::Positive)
            * thermal_distribution_atom_from_index(EdgeIndex::from(2), Sign::Negative)
            - thermal_distribution_atom_from_index(EdgeIndex::from(1), Sign::Negative)
                * thermal_distribution_atom_from_index(EdgeIndex::from(3), Sign::Negative)
                * thermal_distribution_atom_from_index(EdgeIndex::from(2), Sign::Positive);

        assert_eq!(
            numerator.to_atom(&[]).to_canonical_string(),
            expected.to_canonical_string()
        );
    }

    #[test]
    fn canonicalization_swaps_and_tracks_minus_sign() {
        let (canonicalized, sign) = ThermalNumerator::from_edge_lists_canonicalized(
            vec![EdgeIndex::from(1)],
            vec![EdgeIndex::from(4), EdgeIndex::from(2)],
        );

        assert_eq!(sign, Sign::Negative);
        assert_eq!(
            canonicalized.positive_energies,
            vec![EdgeIndex::from(2), EdgeIndex::from(4)]
        );
        assert_eq!(canonicalized.negative_energies, vec![EdgeIndex::from(1)]);
    }

    #[test]
    fn canonicalization_keeps_positive_side_when_already_canonical() {
        let (canonicalized, sign) = ThermalNumerator::from_edge_lists_canonicalized(
            vec![EdgeIndex::from(1), EdgeIndex::from(4)],
            vec![EdgeIndex::from(2), EdgeIndex::from(5)],
        );

        assert_eq!(sign, Sign::Positive);
        assert_eq!(
            canonicalized,
            ThermalNumerator {
                positive_energies: vec![EdgeIndex::from(1), EdgeIndex::from(4)],
                negative_energies: vec![EdgeIndex::from(2), EdgeIndex::from(5)],
            }
        );
    }

    #[test]
    fn trivial_numerator_detection_runs_after_sign_fixup() {
        let (canonicalized, sign) =
            ThermalNumerator::from_edge_lists_canonicalized(vec![], vec![EdgeIndex::from(7)]);

        assert_eq!(sign, Sign::Negative);
        assert!(canonicalized.is_trivial());
    }
}
