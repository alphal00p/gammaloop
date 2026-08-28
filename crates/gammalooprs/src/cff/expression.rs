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
            .map(|orientation| orientation.expression.to_atom_inverse())
            .reduce(|sum, term| sum + term)
            .unwrap_or_default()
    }

    fn get_orientation_atoms(&self, pattern: OrientationPattern) -> TiVec<OrientationId, Atom> {
        self.orientations
            .iter()
            .map(|orientation| {
                if pattern.filter(orientation.orientation()) {
                    orientation.expression.to_atom_inverse()
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
                    orientation.expression.to_atom_inverse()
                } else {
                    Atom::Zero
                };
                (atom, orientation.data.clone())
            })
            .collect()
    }
}
