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
    atom::{Atom, AtomCore, AtomView, Symbol},
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
/// act on these expressions before evaluator construction. Shared tensor and
/// scalar bodies are retained separately; root maps never multiply or mark them.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Integrands {
    atoms: BTreeMap<CutCFFIndex, Atom>,
    numerators: Vec<Arc<FnMapEntry>>,
    scalars: Vec<Arc<FnMapEntry>>,
}

impl Integrands {
    pub fn map<F: FnMut(&Atom) -> Atom>(&self, mut f: F) -> Self {
        Self {
            atoms: self.iter().map(|(key, atom)| (*key, f(atom))).collect(),
            numerators: self.numerators.clone(),
            scalars: self.scalars.clone(),
        }
    }

    pub fn fallible_map<F: FnMut(&Atom) -> Result<Atom>>(&self, mut f: F) -> Result<Self> {
        Ok(Self {
            atoms: self
                .iter()
                .map(|(key, atom)| Ok((*key, f(atom)?)))
                .collect::<Result<_>>()?,
            numerators: self.numerators.clone(),
            scalars: self.scalars.clone(),
        })
    }

    /// Iterate over the factorized expressions passed to evaluators.
    pub fn iter(&self) -> impl Iterator<Item = (&CutCFFIndex, &Atom)> {
        self.atoms.iter()
    }

    pub(crate) fn numerators(&self) -> &[Arc<FnMapEntry>] {
        &self.numerators
    }

    pub(crate) fn scalar_symbol() -> Symbol {
        symbol!("gammalooprs::uv::scalar_coefficient"; Scalar)
    }

    pub(crate) fn scalar_definitions(&self) -> &[Arc<FnMapEntry>] {
        &self.scalars
    }

    /// Scalar bodies contain only arithmetic in formal variables. All physical
    /// functions, momenta, masses, and energy owners remain in the call arguments,
    /// where enclosing forest substitutions can still reach them.
    pub(crate) fn with_scalar_definitions(
        mut self,
        scalars: impl IntoIterator<Item = Arc<FnMapEntry>>,
    ) -> Result<Self> {
        let head = Self::scalar_symbol();
        let mut definitions: BTreeMap<Vec<Atom>, Arc<FnMapEntry>> = BTreeMap::new();
        for entry in scalars {
            if let Some(existing) = definitions.get(&entry.tags) {
                if Arc::ptr_eq(existing, &entry) || existing == &entry {
                    continue;
                }
                return Err(eyre!(
                    "conflicting retained scalar definition for {}",
                    entry.lhs
                ));
            }
            let parameters = entry
                .args
                .iter()
                .cloned()
                .map(Atom::from)
                .collect::<Vec<_>>();
            if entry.tags.is_empty()
                || entry.is_alias
                || entry.tags.iter().any(|tag| tag.contains_symbol(head))
                || entry.lhs != head.call_args(entry.tags.iter().chain(&parameters))
                || parameters.iter().enumerate().any(|(index, parameter)| {
                    !matches!(parameter.as_view(), AtomView::Var(_))
                        || entry.tags.iter().any(|tag| tag.contains(parameter))
                        || parameters[..index].contains(parameter)
                })
            {
                return Err(eyre!("invalid retained scalar binding {}", entry.lhs));
            }
            let mut invalid = false;
            entry.rhs.visitor(&mut |view| {
                let valid = match view {
                    AtomView::Num(_) | AtomView::Add(_) | AtomView::Mul(_) | AtomView::Pow(_) => {
                        true
                    }
                    AtomView::Var(_) => parameters
                        .iter()
                        .any(|parameter| parameter.as_view() == view),
                    _ => false,
                };
                invalid |= !valid;
                valid
            });
            if invalid {
                return Err(eyre!(
                    "retained scalar bodies must contain only arithmetic in their formals: {}",
                    entry.lhs
                ));
            }
            definitions.insert(entry.tags.clone(), entry);
        }
        self.scalars = definitions.into_values().collect();
        Ok(self)
    }

    /// Prepare bindings once for a complete residue scope. The returned resolver
    /// shares instantiated scalar bodies across every root and numerator in it.
    pub(crate) fn scalar_resolver(
        scalars: impl IntoIterator<Item = Arc<FnMapEntry>>,
    ) -> Result<impl FnMut(&Atom) -> Result<Atom>> {
        let validated = Self::from_iter([]).with_scalar_definitions(scalars)?;
        let mut replacements = BTreeMap::new();
        for entry in &validated.scalars {
            replacements
                .entry(entry.tags.len())
                .or_insert_with(BTreeMap::new)
                .insert(
                    entry.tags.clone(),
                    (entry.tags.len() + entry.args.len(), entry.replacement()),
                );
        }
        let head = Self::scalar_symbol();
        let mut cache = BTreeMap::<Atom, Atom>::new();
        Ok(move |atom: &Atom| -> Result<Atom> {
            if !atom.contains_symbol(head) {
                return Ok(atom.clone());
            }
            let mut error = None;
            let result = atom.replace_map_bottom_up(|view, _, out| {
                let AtomView::Fun(call) = view else {
                    return;
                };
                if call.get_symbol() != head || error.is_some() {
                    return;
                }
                let original = view.to_owned();
                if let Some(cached) = cache.get(&original) {
                    **out = cached.clone();
                    return;
                }
                let mut matched = None;
                for (tag_count, definitions) in &replacements {
                    let tags = call
                        .iter()
                        .take(*tag_count)
                        .map(|arg| arg.to_owned())
                        .collect::<Vec<_>>();
                    if let Some((arity, replacement)) = definitions.get(&tags)
                        && *arity == call.get_nargs()
                    {
                        if matched.is_some() {
                            error = Some(eyre!("ambiguous retained scalar call {view}"));
                            return;
                        }
                        matched = Some(replacement);
                    }
                }
                let Some(replacement) = matched else {
                    error = Some(eyre!("unresolved retained scalar call {view}"));
                    return;
                };
                let resolved = view.replace_multiple([replacement]);
                if resolved.contains_symbol(head) {
                    error = Some(eyre!("unresolved or cyclic retained scalar call {view}"));
                    return;
                }
                cache.insert(original, resolved.clone());
                **out = resolved;
            });
            if let Some(error) = error {
                Err(error)
            } else {
                Ok(result)
            }
        })
    }

    #[cfg(test)]
    pub(crate) fn resolve_scalar(&self, atom: &Atom) -> Result<Atom> {
        Self::scalar_resolver(self.scalars.iter().cloned())?(atom)
    }

    /// Resolve innermost scalar calls first, without distributing products or
    /// materializing shared graph numerators. Flat arithmetic definitions cannot
    /// create new calls; nested calls may occur only in actual arguments.
    pub(crate) fn resolved_scalars(&self) -> Result<Self> {
        let mut resolve = Self::scalar_resolver(self.scalars.iter().cloned())?;
        let mut resolved = self.fallible_map(&mut resolve)?.map_numerators(resolve)?;
        resolved.scalars.clear();
        Ok(resolved)
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

    /// Resolve scalar arguments and bodies before materializing tensor families.
    /// A single numerator substitution pass then suffices for its flat store;
    /// unresolved or cyclic calls must not escape into a semantic view.
    pub(crate) fn resolved(&self) -> Result<Self> {
        let validated = self.resolved_scalars()?;
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
            scalars: Vec::new(),
        }
        .with_numerators(self.numerators.iter().chain(&other.numerators).cloned())?
        .with_scalar_definitions(self.scalars.iter().chain(&other.scalars).cloned())
    }

    pub fn zip_mul(&self, other: &Integrands) -> Result<Integrands> {
        self.checked_zip(other, |_, v1, v2| Ok(v1 * v2))
    }

    /// Validate every cut-key shape, then merge each symbolic sum once. Pairwise
    /// accumulation repeatedly copies the already assembled residue numerator.
    pub fn zip_add(self, others: impl IntoIterator<Item = Self>) -> Result<Self> {
        let mut numerators = self.numerators;
        let mut scalars = self.scalars;
        let mut terms = self
            .atoms
            .into_iter()
            .map(|(key, atom)| (key, vec![atom]))
            .collect::<BTreeMap<_, _>>();
        for other in others {
            numerators.extend(other.numerators);
            scalars.extend(other.scalars);
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
            scalars: Vec::new(),
        }
        .with_numerators(numerators)?
        .with_scalar_definitions(scalars)
    }
}

impl FromIterator<(CutCFFIndex, Atom)> for Integrands {
    fn from_iter<I: IntoIterator<Item = (CutCFFIndex, Atom)>>(iter: I) -> Self {
        Self {
            atoms: BTreeMap::from_iter(iter),
            numerators: Vec::new(),
            scalars: Vec::new(),
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
#[cfg(test)]
pub(crate) mod overlap_control;
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
