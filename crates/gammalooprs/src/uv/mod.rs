use std::{
    collections::BTreeMap,
    hash::Hash,
    ops::{Mul, Neg},
    sync::Arc,
};

use crate::{
    GammaLoopContext, cff::CutCFFIndex, integrands::process::param_builder::FnMapEntry,
    numerator::aind::Aind, utils::GS, uv::approx::Rooted,
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
use symbolica::{
    atom::{Atom, AtomCore},
    symbol,
};

use linnet::half_edge::involution::HedgePair;

// use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

pub(crate) fn spenso_lor(
    tag: i32,
    ind: impl Into<Aind>,
    dim: impl Into<Dimension>,
) -> ShadowedStructure<Aind> {
    let mink = Minkowski {}.new_slot(dim, ind);
    NamedStructure::from_iter([mink], GS.emr_mom, Some(vec![Atom::num(tag)])).into_canonical()
}

pub(crate) fn spenso_lor_atom(tag: i32, ind: impl Into<Aind>, dim: impl Into<Dimension>) -> Atom {
    spenso_lor(tag, ind, dim).to_symbolic(None).unwrap()
}

/// Cut-indexed factorized integrands. UV markers and final tensor replacements
/// act on these expressions before evaluator construction. Shared tensor-family
/// bodies are retained separately; ordinary root maps never multiply or mark them.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Integrands {
    atoms: BTreeMap<CutCFFIndex, Atom>,
    numerators: Vec<Arc<FnMapEntry>>,
}

impl Integrands {
    pub fn map<F: FnMut(&Atom) -> Atom>(&self, mut f: F) -> Self {
        Self {
            atoms: self.iter().map(|(key, atom)| (*key, f(atom))).collect(),
            numerators: self.numerators.clone(),
        }
    }

    pub fn fallible_map<F: FnMut(&Atom) -> Result<Atom>>(&self, mut f: F) -> Result<Self> {
        Ok(Self {
            atoms: self
                .iter()
                .map(|(key, atom)| Ok((*key, f(atom)?)))
                .collect::<Result<_>>()?,
            numerators: self.numerators.clone(),
        })
    }

    /// Iterate over the factorized expressions passed to evaluators.
    pub fn iter(&self) -> impl Iterator<Item = (&CutCFFIndex, &Atom)> {
        self.atoms.iter()
    }

    pub(crate) fn numerators(&self) -> &[Arc<FnMapEntry>] {
        &self.numerators
    }

    /// Replace the complete flat definition store. Equal definitions share one
    /// entry; a reused call with a different body or binding is an error.
    pub(crate) fn with_numerators(
        mut self,
        numerators: impl IntoIterator<Item = Arc<FnMapEntry>>,
    ) -> Result<Self> {
        let family = symbol!("gammalooprs::uv::numerator_family");
        let mut definitions: BTreeMap<Vec<Atom>, Arc<FnMapEntry>> = BTreeMap::new();
        for numerator in numerators {
            let parameters = numerator
                .args
                .iter()
                .cloned()
                .map(Atom::from)
                .collect::<Vec<_>>();
            if numerator.tags.is_empty()
                || numerator.lhs
                    != family.call_args(
                        numerator
                            .tags
                            .iter()
                            .cloned()
                            .chain(parameters.iter().cloned()),
                    )
                // Sequential formal substitution must not capture a fixed tag
                // or alter another formal key before that key is substituted.
                || parameters.iter().enumerate().any(|(index, parameter)| {
                    numerator.tags.iter().any(|tag| tag.contains(parameter))
                        || parameters[..index].iter().any(|other| {
                            other.contains(parameter) || parameter.contains(other)
                        })
                })
            {
                return Err(eyre!(
                    "invalid retained numerator binding {}",
                    numerator.lhs
                ));
            }
            if let Some(existing) = definitions.get(&numerator.tags) {
                if existing != &numerator {
                    return Err(eyre!(
                        "conflicting retained numerator definition for {}",
                        numerator.lhs
                    ));
                }
                continue;
            }
            if numerator.rhs.contains_symbol(family) {
                return Err(eyre!(
                    "retained numerator definitions must be flat: {} contains a family call",
                    numerator.lhs
                ));
            }
            definitions.insert(numerator.tags.clone(), numerator);
        }
        self.numerators = definitions.into_values().collect();
        Ok(self)
    }

    /// Transform shared bodies explicitly, preserving their call signatures and
    /// the cut-indexed roots. Unchanged bodies retain their existing Arc.
    pub(crate) fn map_numerators(
        &self,
        mut map: impl FnMut(&Atom) -> Result<Atom>,
    ) -> Result<Self> {
        let numerators = self
            .numerators
            .iter()
            .map(|entry| {
                let rhs = map(&entry.rhs)?;
                Ok(if rhs == entry.rhs {
                    Arc::clone(entry)
                } else {
                    Arc::new(FnMapEntry {
                        rhs,
                        ..entry.as_ref().clone()
                    })
                })
            })
            .collect::<Result<Vec<_>>>()?;
        self.clone().with_numerators(numerators)
    }

    /// Materialize a semantic view with the existing argument-aware replacement
    /// rules. A single substitution pass suffices for this flat store; any
    /// remaining family call is unresolved or cyclic and must not escape.
    pub(crate) fn resolved(&self) -> Result<Self> {
        let validated = self
            .clone()
            .with_numerators(self.numerators.iter().cloned())?;
        let replacements = validated
            .numerators
            .iter()
            .map(|entry| entry.replacement())
            .collect::<Vec<_>>();
        let mut resolved = validated.map(|atom| atom.replace_multiple(&replacements));
        let family = symbol!("gammalooprs::uv::numerator_family");
        if resolved
            .iter()
            .any(|(_, atom)| atom.contains_symbol(family))
        {
            return Err(eyre!("unresolved or cyclic retained numerator family call"));
        }
        resolved.numerators.clear();
        Ok(resolved)
    }

    pub(crate) fn zero_like(&self) -> Self {
        self.map(|_| Atom::Zero)
    }

    pub fn checked_zip(
        &self,
        other: &Integrands,
        mut map: impl FnMut(&CutCFFIndex, &Atom, &Atom) -> Result<Atom>,
    ) -> Result<Integrands> {
        let atoms = self
            .iter()
            .merge_join_by(other.iter(), |(left_key, _), (right_key, _)| {
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
            .collect::<Result<_>>()?;
        Self {
            atoms,
            numerators: Vec::new(),
        }
        .with_numerators(self.numerators.iter().chain(&other.numerators).cloned())
    }

    pub fn zip_mul(&self, other: &Integrands) -> Result<Integrands> {
        self.checked_zip(other, |_, v1, v2| Ok(v1 * v2))
    }

    /// Validate every cut-key shape, then merge each symbolic sum once. Pairwise
    /// accumulation repeatedly copies the already assembled residue numerator.
    pub fn zip_add(self, others: impl IntoIterator<Item = Self>) -> Result<Self> {
        let mut numerators = self.numerators;
        let mut terms = self
            .atoms
            .into_iter()
            .map(|(key, atom)| (key, vec![atom]))
            .collect::<BTreeMap<_, _>>();
        for other in others {
            numerators.extend(other.numerators);
            for pair in terms
                .iter_mut()
                .merge_join_by(other.atoms, |(left, _), (right, _)| (*left).cmp(right))
            {
                match pair {
                    EitherOrBoth::Both((_, terms), (_, atom)) => terms.push(atom),
                    EitherOrBoth::Left((key, _)) => {
                        return Err(eyre!("right integrands are missing key {key:?}"));
                    }
                    EitherOrBoth::Right((key, _)) => {
                        return Err(eyre!("left integrands are missing key {key:?}"));
                    }
                }
            }
        }
        Self {
            atoms: terms
                .into_iter()
                .map(|(key, terms)| (key, Atom::add_many(terms)))
                .collect(),
            numerators: Vec::new(),
        }
        .with_numerators(numerators)
    }
}

impl FromIterator<(CutCFFIndex, Atom)> for Integrands {
    fn from_iter<I: IntoIterator<Item = (CutCFFIndex, Atom)>>(iter: I) -> Self {
        Self {
            atoms: BTreeMap::from_iter(iter),
            numerators: Vec::new(),
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
        [(CutCFFIndex::new_all_none(), Atom::num(1))]
            .into_iter()
            .collect()
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

pub mod export;

pub mod forest;
pub use forest::Forest;

pub mod profile;
pub use profile::{UVProfile, UVProfileAnalysis, UVProfilePassFail};

#[cfg(test)]
mod tests;
