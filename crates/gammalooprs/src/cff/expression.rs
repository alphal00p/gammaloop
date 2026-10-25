//! GammaLoop extensions for FeynKit-owned CFF expressions and orientations.

use std::borrow::Borrow;

use feynkit_cff::{CffResult, OrientationData, OrientationExpression, OrientationId};
use linnet::half_edge::involution::{EdgeVec, Orientation};
use symbolica::{
    atom::Atom,
    id::{Pattern, Replacement},
};
use typed_index_collections::TiVec;

use crate::{
    cff::orientations::GraphOrientation, settings::global::OrientationPattern,
    utils::ose_atom_from_index,
};

impl GraphOrientation for OrientationData {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        &self.orientation
    }
}

impl GraphOrientation for feynkit_cff::GraphOrientation {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        self
    }
}

impl GraphOrientation for OrientationExpression {
    fn orientation(&self) -> &EdgeVec<Orientation> {
        &self.data.orientation
    }
}

/// GammaLoop diagnostics and symbolic rewrites for a FeynKit orientation.
pub trait OrientationDataExt {
    fn get_ose_replacements(&self) -> Vec<Replacement>;
}

impl OrientationDataExt for OrientationData {
    fn get_ose_replacements(&self) -> Vec<Replacement> {
        self.orientation
            .borrow()
            .into_iter()
            .filter(|(_, orientation)| **orientation == Orientation::Reversed)
            .map(|(edge, _)| {
                let energy = ose_atom_from_index(edge);
                Replacement::new(Pattern::from(energy.clone()), Pattern::from(-energy))
            })
            .collect()
    }
}

/// GammaLoop's symbolic lowering and residue API over a canonical FeynKit result.
pub trait CffExpressionExt {
    fn to_atom(&self, pattern: OrientationPattern) -> Atom;
    fn get_orientation_atoms(&self, pattern: OrientationPattern) -> TiVec<OrientationId, Atom>;
    fn get_orientation_atoms_with_data(
        &self,
        pattern: OrientationPattern,
    ) -> TiVec<OrientationId, (Atom, OrientationData)>;
}

impl CffExpressionExt for CffResult {
    fn to_atom(&self, pattern: OrientationPattern) -> Atom {
        self.orientations
            .iter()
            .filter(|orientation| pattern.filter(orientation.orientation()))
            .map(|orientation| orientation.to_atom())
            .reduce(|sum, term| sum + term)
            .unwrap_or_default()
    }

    fn get_orientation_atoms(&self, pattern: OrientationPattern) -> TiVec<OrientationId, Atom> {
        self.orientations
            .iter()
            .map(|orientation| {
                if pattern.filter(orientation.orientation()) {
                    orientation.to_atom()
                } else {
                    Atom::Zero
                }
            })
            .collect()
    }

    fn get_orientation_atoms_with_data(
        &self,
        pattern: OrientationPattern,
    ) -> TiVec<OrientationId, (Atom, OrientationData)> {
        self.orientations
            .iter()
            .map(|orientation| {
                let atom = if pattern.filter(orientation.orientation()) {
                    orientation.to_atom()
                } else {
                    Atom::Zero
                };
                (atom, orientation.data.clone())
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use feynkit_cff::{CffEdge, CffGenerator, CffGraph, EdgeId, VertexId};

    use super::*;
    use crate::initialisation::test_initialise;

    #[test]
    fn all_root_lowerings_preserve_free_tadpole_weights_and_selection() {
        test_initialise().unwrap();
        assert!(symbolica::license::LicenseManager::is_licensed());

        for count in [1, 2] {
            let graph = CffGraph::new(
                1,
                (0..count).map(|edge| {
                    CffEdge::internal(EdgeId::new(edge), VertexId::new(0), VertexId::new(0))
                }),
            )
            .unwrap();
            let result = CffGenerator::default().generate(&graph).unwrap();
            assert_eq!(result.orientations.len(), 1 << count);
            let weight = Atom::num((1_i64, 1_i64 << count));
            let selected =
                OrientationPattern::from_orientation(&result.orientations[OrientationId(0)]);
            let rejected = OrientationPattern::from(Atom::Zero);

            for (pattern, expected_total) in [
                (OrientationPattern::default(), Atom::num(1)),
                (selected, weight.clone()),
                (rejected, Atom::Zero),
            ] {
                // A free self-loop contributes one half per contour direction.
                // This sum is authored independently of the owner's lowering API.
                let explicit_sum = result
                    .orientations
                    .iter()
                    .filter(|orientation| pattern.filter(orientation.orientation()))
                    .fold(Atom::Zero, |sum, _| sum + weight.clone());
                assert_eq!(explicit_sum, expected_total);
                assert_eq!(
                    CffExpressionExt::to_atom(&result, pattern.clone()),
                    expected_total
                );

                let atoms = result.get_orientation_atoms(pattern.clone());
                let with_data = result.get_orientation_atoms_with_data(pattern.clone());
                assert_eq!(atoms.len(), result.orientations.len());
                assert_eq!(with_data.len(), result.orientations.len());
                for (index, orientation) in result.orientations.iter().enumerate() {
                    let expected = if pattern.filter(orientation.orientation()) {
                        weight.clone()
                    } else {
                        Atom::Zero
                    };
                    let id = OrientationId(index);
                    assert_eq!(atoms[id], expected);
                    assert_eq!(with_data[id].0, expected);
                    assert_eq!(with_data[id].1, orientation.data);
                }
                assert_eq!(
                    atoms.into_iter().fold(Atom::Zero, |sum, term| sum + term),
                    expected_total
                );
            }
        }
    }
}
