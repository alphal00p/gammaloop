use std::ops::Neg;

use color_eyre::Result;
use symbolica::{atom::Atom, symbol};

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
            active_lmb: None,
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
            active_lmb: self.active_lmb,
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

    pub(crate) fn zip_add(&self, other: &Self) -> Result<Self> {
        // Deferred integrands expose nonzero evaluator calls through `iter`, so
        // this predicate only prunes compact branches proven to be zero.
        let is_zero = |branch: &OrientationIntegrandBranch| {
            branch.integrands.iter().all(|(_, atom)| atom.is_zero())
        };
        let fallback_zero = self
            .0
            .iter()
            .chain(&other.0)
            .find(|branch| is_zero(branch))
            .cloned();
        let mut branches = self
            .0
            .iter()
            .filter(|branch| !is_zero(branch))
            .cloned()
            .collect::<Vec<_>>();
        for right in other.0.iter().filter(|branch| !is_zero(branch)) {
            if let Some(left) = branches.iter_mut().find(|left| {
                left.selector_id == right.selector_id
                    && left.source_edge_energy_map == right.source_edge_energy_map
            }) {
                left.integrands = left.integrands.clone().zip_add(right.integrands.clone())?;
            } else {
                branches.push(right.clone());
            }
        }
        if branches.is_empty() {
            branches.extend(fallback_zero);
        }
        Ok(Self(branches))
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

    /// Multiply branches hosted by the same production selector. A branch
    /// absent on either side contributes zero, so independently projected
    /// factors retain only their selector intersection. The mapper sees the
    /// complete active factor only after its host has been selected.
    pub(crate) fn zip_mul_mapped_factor(
        &self,
        factor: &Self,
        mut map: impl FnMut(OrientationID, Option<&[LinearEnergyExpr]>, &Atom) -> Result<Atom>,
    ) -> Result<Self> {
        self.0
            .iter()
            .filter_map(|outer| {
                let matching = factor
                    .0
                    .iter()
                    .filter(|inner| inner.selector_id == outer.selector_id)
                    .collect::<Vec<_>>();
                (!matching.is_empty()).then(|| {
                    let mut product = outer.integrands.zero_like();
                    for inner in matching {
                        let mapped = inner.integrands.fallible_map(|atom| {
                            map(
                                outer.selector_id,
                                outer.source_edge_energy_map.as_deref(),
                                atom,
                            )
                        })?;
                        product = product.zip_add(outer.integrands.zip_mul(&mapped)?)?;
                    }
                    Ok(OrientationIntegrandBranch {
                        selector_id: outer.selector_id,
                        source_edge_energy_map: outer.source_edge_energy_map.clone(),
                        integrands: product,
                    })
                })
            })
            .collect::<Result<Vec<_>>>()
            .map(Self)
    }

    /// A factorized additive projection for coefficient diagnostics. Repeated
    /// selector hosts of the same active expression are included only once.
    #[cfg(test)]
    pub(crate) fn factorized_sum(&self) -> Atom {
        let mut distinct = Vec::new();
        for atom in self
            .0
            .iter()
            .flat_map(|branch| branch.integrands.iter().map(|(_, atom)| atom))
        {
            if !distinct.contains(atom) {
                distinct.push(atom.clone());
            }
        }
        distinct
            .into_iter()
            .fold(Atom::Zero, |sum, atom| sum + atom)
    }

    /// Keep independently evaluated production branches algebraically
    /// independent while deriving one conservative outer-CFF capacity.
    /// Branch tags are analysis-only scalar coefficients: they neither expand
    /// the factorized atoms nor enter the mapped production numerator.
    pub(crate) fn factorized_capacity_envelope(&self) -> Atom {
        self.0
            .iter()
            .flat_map(|branch| branch.integrands.iter().map(|(_, atom)| atom))
            .filter(|atom| !atom.is_zero())
            .enumerate()
            .fold(Atom::Zero, |sum, (branch, atom)| {
                let tag = Atom::var(symbol!(format!(
                    "__gammaloop_outer_cff_capacity_branch_{branch}"
                )));
                sum + tag * atom
            })
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

    pub(in crate::uv::approx) fn fallible_map(
        &self,
        mut f: impl FnMut(OrientationID, Option<&[LinearEnergyExpr]>, &Atom) -> Result<Atom>,
    ) -> Result<Self> {
        self.0
            .iter()
            .map(|branch| {
                Ok(OrientationIntegrandBranch {
                    selector_id: branch.selector_id,
                    source_edge_energy_map: branch.source_edge_energy_map.clone(),
                    integrands: branch.integrands.fallible_map(|atom| {
                        f(
                            branch.selector_id,
                            branch.source_edge_energy_map.as_deref(),
                            atom,
                        )
                    })?,
                })
            })
            .collect::<Result<Vec<_>>>()
            .map(Self)
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
