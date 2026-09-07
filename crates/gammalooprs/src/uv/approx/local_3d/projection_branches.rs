use std::ops::Neg;

use color_eyre::Result;
use symbolica::atom::Atom;

use crate::{
    cff::{CutCFFIndex, expression::OrientationID, surface::LinearEnergyExpr},
    uv::{Integrands, approx::Rooted},
};

use super::{FrozenActiveCt, OrientationIntegrandBranch, OrientationIntegrands};

impl FrozenActiveCt {
    pub(crate) fn combine(&self) -> Result<OrientationIntegrands> {
        self.active.zip_mul_unmapped(&self.frozen_integrands)
    }
}

impl From<OrientationIntegrands> for FrozenActiveCt {
    fn from(active: OrientationIntegrands) -> Self {
        let frozen_integrands = active
            .0
            .first()
            .map(|branch| {
                branch
                    .integrands
                    .iter()
                    .map(|(index, _)| (*index, Atom::one()))
                    .collect()
            })
            .unwrap_or_else(Integrands::root);
        Self {
            active,
            frozen_integrands,
        }
    }
}

impl From<Integrands> for FrozenActiveCt {
    fn from(integrands: Integrands) -> Self {
        OrientationIntegrands::from(integrands).into()
    }
}

impl Neg for FrozenActiveCt {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            active: -self.active,
            frozen_integrands: self.frozen_integrands,
        }
    }
}
impl From<Integrands> for OrientationIntegrands {
    fn from(integrands: Integrands) -> Self {
        Self(vec![OrientationIntegrandBranch {
            selector_id: OrientationID(0),
            source_edge_energy_map: None,
            integrands,
        }])
    }
}

impl OrientationIntegrands {
    pub(super) fn from_ids_and_indices(
        ids: impl IntoIterator<Item = OrientationID>,
        indices: &[CutCFFIndex],
    ) -> Self {
        Self(
            ids.into_iter()
                .map(|selector_id| OrientationIntegrandBranch {
                    selector_id,
                    source_edge_energy_map: None,
                    integrands: indices.iter().map(|index| (*index, Atom::Zero)).collect(),
                })
                .collect(),
        )
    }

    pub(crate) fn zip_mul_unmapped(&self, other: &Integrands) -> Result<Self> {
        self.0
            .iter()
            .map(|branch| {
                Ok(OrientationIntegrandBranch {
                    selector_id: branch.selector_id,
                    source_edge_energy_map: branch.source_edge_energy_map.clone(),
                    integrands: branch.integrands.zip_mul(other)?,
                })
            })
            .collect::<Result<Vec<_>>>()
            .map(Self)
    }

    pub(crate) fn multiply_mapped(
        &self,
        mut map: impl FnMut(OrientationID, Option<&[LinearEnergyExpr]>) -> Result<Atom>,
    ) -> Result<Self> {
        self.0
            .iter()
            .map(|branch| {
                let mapped = map(branch.selector_id, branch.source_edge_energy_map.as_deref())?;
                Ok(OrientationIntegrandBranch {
                    selector_id: branch.selector_id,
                    source_edge_energy_map: branch.source_edge_energy_map.clone(),
                    integrands: branch.integrands.map(|atom| atom * &mapped),
                })
            })
            .collect::<Result<Vec<_>>>()
            .map(Self)
    }

    pub(crate) fn map(&self, mut f: impl FnMut(&Atom) -> Atom) -> Self {
        Self(
            self.0
                .iter()
                .map(|branch| OrientationIntegrandBranch {
                    selector_id: branch.selector_id,
                    source_edge_energy_map: branch.source_edge_energy_map.clone(),
                    integrands: branch.integrands.map(&mut f),
                })
                .collect(),
        )
    }

    #[cfg(test)]
    pub(crate) fn iter(&self) -> impl Iterator<Item = (&CutCFFIndex, &Atom)> {
        self.0.iter().flat_map(|branch| branch.integrands.iter())
    }

    pub(crate) fn iter_orientations(
        &self,
    ) -> impl Iterator<Item = (OrientationID, Option<&[LinearEnergyExpr]>, &Integrands)> {
        self.0.iter().map(|branch| {
            (
                branch.selector_id,
                branch.source_edge_energy_map.as_deref(),
                &branch.integrands,
            )
        })
    }
}

impl FromIterator<(OrientationID, Integrands)> for OrientationIntegrands {
    fn from_iter<T: IntoIterator<Item = (OrientationID, Integrands)>>(iter: T) -> Self {
        Self(
            iter.into_iter()
                .map(|(selector_id, integrands)| OrientationIntegrandBranch {
                    selector_id,
                    source_edge_energy_map: None,
                    integrands,
                })
                .collect(),
        )
    }
}

impl Neg for OrientationIntegrands {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(
            self.0
                .into_iter()
                .map(|branch| OrientationIntegrandBranch {
                    selector_id: branch.selector_id,
                    source_edge_energy_map: branch.source_edge_energy_map,
                    integrands: -branch.integrands,
                })
                .collect(),
        )
    }
}
