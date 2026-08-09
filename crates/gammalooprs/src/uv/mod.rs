use std::{
    collections::BTreeMap,
    hash::Hash,
    ops::{Mul, Neg},
};

use crate::{
    GammaLoopContext, cff::CutCFFIndex, numerator::aind::Aind, utils::GS, uv::approx::Rooted,
};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::eyre;
use itertools::{EitherOrBoth, Itertools};
use spenso::{
    network::parsing::ShadowedStructure,
    structure::{
        NamedStructure, ToSymbolic,
        dimension::Dimension,
        representation::{Minkowski, RepName},
    },
};
use symbolica::{atom::Atom, function};

use linnet::half_edge::involution::HedgePair;

// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

pub(crate) fn spenso_lor(
    tag: i32,
    ind: impl Into<Aind>,
    dim: impl Into<Dimension>,
) -> ShadowedStructure<Aind> {
    let mink = Minkowski {}.new_slot(dim, ind);
    NamedStructure::from_iter([mink], GS.emr_mom, Some(vec![Atom::num(tag)])).structure
}

pub(crate) fn spenso_lor_atom(tag: i32, ind: impl Into<Aind>, dim: impl Into<Dimension>) -> Atom {
    spenso_lor(tag, ind, dim).to_symbolic(None).unwrap()
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct IntegrandExpr {
    integrands: BTreeMap<CutCFFIndex, Atom>,
    // add_arg: Option<Atom>,
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Integrands {
    compact: BTreeMap<CutCFFIndex, Atom>,
    deferred: Option<DeferredIntegrands>,
}

impl Integrands {
    pub fn map<F: FnMut(&Atom) -> Atom>(&self, mut f: F) -> Self {
        match &self.deferred {
            Some(deferred) => Self::from_deferred(deferred.map(f)),
            None => Self {
                compact: self
                    .compact
                    .iter()
                    .map(|(key, atom)| (*key, f(atom)))
                    .collect(),
                deferred: None,
            },
        }
    }

    pub fn fallible_map<F: FnMut(&Atom) -> Result<Atom>>(&self, mut f: F) -> Result<Self> {
        match &self.deferred {
            Some(deferred) => deferred.fallible_map(f).map(Self::from_deferred),
            None => Ok(Self {
                compact: self
                    .compact
                    .iter()
                    .map(|(key, atom)| Ok((*key, f(atom)?)))
                    .collect::<Result<_>>()?,
                deferred: None,
            }),
        }
    }

    /// Iterate over the compact expression passed to evaluators. Deferred
    /// projected-CFF bodies remain available through [`Self::deferred_terms`].
    pub fn iter(&self) -> impl Iterator<Item = (&CutCFFIndex, &Atom)> {
        self.compact.iter()
    }

    pub(crate) fn from_deferred(deferred: DeferredIntegrands) -> Self {
        let compact = deferred.compact_calls();
        Self {
            compact,
            deferred: Some(deferred),
        }
    }

    pub(crate) fn deferred_terms(&self, index: &CutCFFIndex) -> Option<&[Atom]> {
        self.deferred.as_ref()?.terms(index)
    }

    pub(crate) fn materialize(&self) -> Self {
        self.deferred
            .as_ref()
            .map_or_else(|| self.clone(), DeferredIntegrands::materialize)
    }

    pub(crate) fn zero_like(&self) -> Self {
        Self {
            compact: self
                .compact
                .keys()
                .map(|index| (*index, Atom::Zero))
                .collect(),
            deferred: None,
        }
    }

    pub fn checked_zip(
        &self,
        other: &Integrands,
        mut map: impl FnMut(&CutCFFIndex, &Atom, &Atom) -> Result<Atom>,
    ) -> Result<Integrands> {
        if self.deferred.is_some() || other.deferred.is_some() {
            return Err(eyre!(
                "checked integrand zips require materialized deferred projected-CFF sums"
            ));
        }
        self.compact
            .iter()
            .merge_join_by(&other.compact, |(left_key, _), (right_key, _)| {
                left_key.cmp(right_key)
            })
            .map(|pair| match pair {
                EitherOrBoth::Both((key, left), (_, right)) => Ok((*key, map(key, left, right)?)),
                EitherOrBoth::Left((key, _)) => {
                    Err(eyre!("right integrands are missing key {key:?}"))
                }
                EitherOrBoth::Right((key, _)) => {
                    Err(eyre!("left integrands are missing key {key:?}"))
                }
            })
            .collect()
    }

    pub fn zip_mul(&self, other: &Integrands) -> Result<Integrands> {
        self.checked_zip(other, |_, v1, v2| Ok(v1 * v2))
    }
    pub fn zip_add(mut self, mut other: Integrands) -> Result<Integrands> {
        match (self.deferred.take(), other.deferred.take()) {
            (None, None) => self.checked_zip(&other, |_, left, right| Ok(left + right)),
            (Some(left), Some(right)) => left.zip_add(right).map(Self::from_deferred),
            (Some(left), None) => left
                .zip_add(DeferredIntegrands::from_compact(other.compact))
                .map(Self::from_deferred),
            (None, Some(right)) => DeferredIntegrands::from_compact(self.compact)
                .zip_add(right)
                .map(Self::from_deferred),
        }
    }
}

impl FromIterator<(CutCFFIndex, Atom)> for Integrands {
    fn from_iter<I: IntoIterator<Item = (CutCFFIndex, Atom)>>(iter: I) -> Self {
        Self {
            compact: BTreeMap::from_iter(iter),
            deferred: None,
        }
    }
}

impl Mul<Atom> for Integrands {
    type Output = Self;

    fn mul(self, rhs: Atom) -> Self::Output {
        self.map(|a| a * &rhs)
    }
}

impl Neg for Integrands {
    type Output = Self;

    fn neg(self) -> Self::Output {
        self.map(|a| a.neg())
    }
}

impl Mul<&Atom> for Integrands {
    type Output = Self;

    fn mul(self, rhs: &Atom) -> Self::Output {
        self.map(|a| a * rhs)
    }
}

impl Rooted for Integrands {
    fn root() -> Self {
        Self {
            compact: BTreeMap::from([(CutCFFIndex::new_all_none(), Atom::num(1))]),
            deferred: None,
        }
    }
}

/// An additive integrand whose branches stay separate until evaluator
/// construction. This keeps post-4D CFF orientation sums from being expanded
/// while UV markers and final tensor replacements still need to see each body.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct DeferredIntegrands(BTreeMap<CutCFFIndex, Vec<Atom>>);

impl DeferredIntegrands {
    pub(crate) fn from_indices(indices: impl IntoIterator<Item = CutCFFIndex>) -> Self {
        Self(
            indices
                .into_iter()
                .map(|index| (index, Vec::new()))
                .collect(),
        )
    }

    fn from_compact(integrands: BTreeMap<CutCFFIndex, Atom>) -> Self {
        Self(
            integrands
                .into_iter()
                .map(|(index, atom)| (index, vec![atom]))
                .collect(),
        )
    }

    pub(crate) fn push(&mut self, index: CutCFFIndex, atom: Atom) -> Result<()> {
        self.0
            .get_mut(&index)
            .ok_or_else(|| eyre!("deferred integrands are missing key {index:?}"))?
            .push(atom);
        Ok(())
    }

    pub(crate) fn map(&self, mut f: impl FnMut(&Atom) -> Atom) -> Self {
        Self(
            self.0
                .iter()
                .map(|(index, atoms)| (*index, atoms.iter().map(&mut f).collect::<Vec<_>>()))
                .collect(),
        )
    }

    pub(crate) fn fallible_map(&self, mut f: impl FnMut(&Atom) -> Result<Atom>) -> Result<Self> {
        self.0
            .iter()
            .map(|(index, atoms)| {
                Ok((
                    *index,
                    atoms.iter().map(&mut f).collect::<Result<Vec<_>>>()?,
                ))
            })
            .collect::<Result<BTreeMap<_, _>>>()
            .map(Self)
    }

    pub(crate) fn zip_add(self, other: Self) -> Result<Self> {
        self.0
            .into_iter()
            .merge_join_by(other.0, |(left, _), (right, _)| left.cmp(right))
            .map(|pair| match pair {
                EitherOrBoth::Both((index, mut left), (_, right)) => {
                    left.extend(right);
                    Ok((index, left))
                }
                EitherOrBoth::Left((index, _)) => {
                    Err(eyre!("right deferred integrands are missing key {index:?}"))
                }
                EitherOrBoth::Right((index, _)) => {
                    Err(eyre!("left deferred integrands are missing key {index:?}"))
                }
            })
            .collect::<Result<BTreeMap<_, _>>>()
            .map(Self)
    }

    pub(crate) fn materialize(&self) -> Integrands {
        self.0
            .iter()
            .map(|(index, atoms)| {
                (
                    *index,
                    atoms
                        .iter()
                        .cloned()
                        .fold(Atom::Zero, |sum, atom| sum + atom),
                )
            })
            .collect()
    }

    fn compact_calls(&self) -> BTreeMap<CutCFFIndex, Atom> {
        self.0
            .iter()
            .map(|(index, atoms)| {
                let calls = atoms
                    .iter()
                    .enumerate()
                    .map(|(tag, _)| function!(GS.projected_cff_sum, tag))
                    .fold(Atom::Zero, |sum, call| sum + call);
                (*index, calls)
            })
            .collect()
    }

    pub(crate) fn terms(&self, index: &CutCFFIndex) -> Option<&[Atom]> {
        self.0.get(index).map(Vec::as_slice)
    }
}
#[allow(dead_code)]
pub(crate) fn is_not_paired(pair: &HedgePair) -> bool {
    !pair.is_paired()
}

pub mod hedge_poset;
mod marker;
mod orchestrator;
pub mod renormalization;
pub use renormalization::{RenormalizationPart, RenormalizationStats};
pub mod settings;
pub use settings::{
    ApproximationType, CTIdentifier, CTRenormalizationRule, RenormalizationPrescriptionSettings,
    UVOrchestrator, UVgenerationSettings,
};
pub mod uv_graph;
pub use uv_graph::UltravioletGraph;

pub mod spinney;
pub use spinney::Spinney;

pub mod poset;
pub use poset::Poset;

pub mod wood;
pub use wood::Wood;

pub mod approx;
pub use approx::ApproxOp;

pub mod export;

pub mod forest;
pub use forest::Forest;

pub mod profile;
pub use profile::{UVProfile, UVProfileAnalysis, UVProfilePassFail};

#[cfg(test)]
mod deferred_integrands_tests {
    use super::*;

    #[test]
    fn deferred_zip_add_preserves_keys_and_branch_order() {
        let root = CutCFFIndex::new_all_none();
        let second = CutCFFIndex {
            left_threshold_order: Some(1),
            right_threshold_order: None,
            lu_cut_order: None,
        };
        let mut left = DeferredIntegrands::from_indices([root, second]);
        left.push(root, Atom::num(1)).unwrap();
        left.push(root, Atom::num(2)).unwrap();
        left.push(second, Atom::num(10)).unwrap();
        let mut right = DeferredIntegrands::from_indices([root, second]);
        right.push(root, Atom::num(3)).unwrap();
        right.push(second, Atom::num(20)).unwrap();
        right.push(second, Atom::num(30)).unwrap();

        let sum = Integrands::from_deferred(left)
            .zip_add(Integrands::from_deferred(right))
            .unwrap();

        assert_eq!(
            sum.iter().map(|(index, _)| *index).collect::<Vec<_>>(),
            vec![root, second]
        );
        assert_eq!(
            sum.deferred_terms(&root).unwrap(),
            [Atom::num(1), Atom::num(2), Atom::num(3)]
        );
        assert_eq!(
            sum.deferred_terms(&second).unwrap(),
            [Atom::num(10), Atom::num(20), Atom::num(30)]
        );
    }

    #[test]
    fn deferred_materialize_and_zero_like_return_direct_integrands() {
        let index = CutCFFIndex::new_all_none();
        let mut deferred = DeferredIntegrands::from_indices([index]);
        deferred.push(index, Atom::num(2)).unwrap();
        deferred.push(index, Atom::num(3)).unwrap();
        let integrands = Integrands::from_deferred(deferred);

        assert_eq!(
            integrands.materialize(),
            [(index, Atom::num(5))].into_iter().collect()
        );
        assert_eq!(
            integrands.zero_like(),
            [(index, Atom::Zero)].into_iter().collect()
        );
    }
}

#[cfg(test)]
mod tests;
