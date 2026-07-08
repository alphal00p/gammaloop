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
use symbolica::atom::Atom;

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
pub struct Integrands(BTreeMap<CutCFFIndex, Atom>);

impl Integrands {
    pub fn map<F: FnMut(&Atom) -> Atom>(&self, mut f: F) -> Self {
        Integrands(self.0.iter().map(|(k, v)| (*k, f(v))).collect())
    }

    pub fn fallible_map<F: FnMut(&Atom) -> Result<Atom>>(&self, mut f: F) -> Result<Self> {
        self.0.iter().map(|(k, v)| Ok((*k, f(v)?))).collect()
    }

    pub fn iter(&self) -> impl Iterator<Item = (&CutCFFIndex, &Atom)> {
        self.0.iter()
    }

    pub fn checked_zip(
        &self,
        other: &Integrands,
        mut map: impl FnMut(&CutCFFIndex, &Atom, &Atom) -> Result<Atom>,
    ) -> Result<Integrands> {
        self.0
            .iter()
            .merge_join_by(&other.0, |(left_key, _), (right_key, _)| {
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
    pub fn zip_add(&self, other: &Integrands) -> Result<Integrands> {
        self.checked_zip(other, |_, v1, v2| Ok(v1 + v2))
    }
}

impl FromIterator<(CutCFFIndex, Atom)> for Integrands {
    fn from_iter<I: IntoIterator<Item = (CutCFFIndex, Atom)>>(iter: I) -> Self {
        Integrands(BTreeMap::from_iter(iter))
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
        Integrands(BTreeMap::from([(
            CutCFFIndex::new_all_none(),
            Atom::num(1),
        )]))
    }
}
#[allow(dead_code)]
pub(crate) fn is_not_paired(pair: &HedgePair) -> bool {
    !pair.is_paired()
}

pub mod hedge_poset;
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

pub mod forest;
pub use forest::Forest;

pub mod profile;
pub use profile::{UVProfile, UVProfileAnalysis, UVProfilePassFail};

#[cfg(test)]
mod tests;
