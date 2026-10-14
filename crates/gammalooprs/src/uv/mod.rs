use std::{
    collections::BTreeMap,
    hash::Hash,
    ops::{Mul, Neg},
    sync::Arc,
};

use crate::{
    GammaLoopContext,
    cff::CutCFFIndex,
    integrands::process::param_builder::FnMapEntry,
    numerator::{aind::Aind, symbolica_ext::NumeratorAtomExt},
    utils::GS,
    uv::approx::Rooted,
};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use eyre::eyre;
use idenso::IndexTooling;
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
    atom::{Atom, AtomCore, AtomView, Indeterminate, Symbol},
    domains::rational::Rational,
    id::Replacement,
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

/// One complete scalar-coefficient row. Its native source keys remain available
/// until all physical mapping and cut selection have finished.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct ParametricIntegrandRow {
    pub(crate) coefficients: Vec<Rational>,
    pub(crate) carrier: Atom,
    pub(crate) source_keys: Vec<approx::direct_3d::DirectResidueKey>,
}

/// One additive numerator product and its Cartesian residue rows. UV operators
/// act on this ordinary generic body, never on its concrete row substitutions.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub(crate) struct ParametricIntegrandTerm {
    pub(crate) numerator: Arc<FnMapEntry>,
    pub(crate) rows: Vec<ParametricIntegrandRow>,
}

impl ParametricIntegrandTerm {
    pub(crate) fn new(
        numerator: Arc<FnMapEntry>,
        rows: Vec<ParametricIntegrandRow>,
    ) -> Result<Self> {
        eyre::ensure!(
            rows.iter()
                .all(|row| row.coefficients.len() == numerator.args.len()),
            "parametric numerator row does not match its formal coefficients"
        );
        Ok(Self { numerator, rows })
    }

    fn product(&self, other: &Self) -> Result<Self> {
        let (scope, tag) = approx::direct_3d::DirectResidueBranches::numerator_scope();
        let parameters = (0..self.numerator.args.len() + other.numerator.args.len())
            .map(|index| {
                symbol!("gammalooprs::uv::numerator_coefficient"; Scalar)
                    .call_args([Atom::num(scope), Atom::num(index)])
            })
            .collect::<Vec<_>>();
        let bind = |entry: &FnMapEntry, parameters: &[Atom]| -> Result<Atom> {
            let interface = entry.rhs.list_dangling::<Aind>()?.into_iter().collect();
            // Each independent factor owns its private contractions. Only its
            // free tensor slots may join the other factor in the product.
            let body = entry.rhs.freshen_private_indices(&interface, || {
                approx::direct_3d::DirectResidueBranches::numerator_scope().1
            });
            Ok(body.replace_multiple(
                entry
                    .args
                    .iter()
                    .cloned()
                    .map(Atom::from)
                    .zip(parameters)
                    .map(|(formal, parameter)| {
                        Replacement::new(formal.to_pattern(), parameter.clone())
                    }),
            ))
        };
        let split = self.numerator.args.len();
        let numerator = Arc::new(FnMapEntry {
            is_alias: false,
            lhs: symbol!("gammalooprs::uv::numerator_family")
                .call_args(std::iter::once(tag.clone()).chain(parameters.iter().cloned())),
            rhs: bind(&self.numerator, &parameters[..split])?
                * bind(&other.numerator, &parameters[split..])?,
            args: parameters
                .into_iter()
                .map(Indeterminate::try_from)
                .collect::<std::result::Result<_, _>>()
                .map_err(|error| eyre!(error))?,
            tags: vec![tag],
            inlining: Default::default(),
        });
        let rows = self
            .rows
            .iter()
            .flat_map(|left| {
                other.rows.iter().map(move |right| ParametricIntegrandRow {
                    coefficients: left
                        .coefficients
                        .iter()
                        .chain(&right.coefficients)
                        .cloned()
                        .collect(),
                    carrier: &left.carrier * &right.carrier,
                    source_keys: left
                        .source_keys
                        .iter()
                        .chain(&right.source_keys)
                        .cloned()
                        .collect(),
                })
            })
            .collect();
        Self::new(numerator, rows)
    }
}

/// Cut-indexed additive integrands. Scalar contributions and parametric terms
/// are disjoint: no concrete copy of a retained row is stored in `atoms`.
#[derive(Clone, Debug, PartialEq, Eq, Hash, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Integrands {
    atoms: BTreeMap<CutCFFIndex, Atom>,
    numerators: Vec<Arc<FnMapEntry>>,
    scalars: Vec<Arc<FnMapEntry>>,
    parametric_terms: BTreeMap<CutCFFIndex, Vec<ParametricIntegrandTerm>>,
}

impl Integrands {
    pub fn map<F: FnMut(&Atom) -> Atom>(&self, mut f: F) -> Self {
        self.fallible_map(|atom| Ok(f(atom)))
            .expect("infallible integrand map")
    }

    pub fn fallible_map<F: FnMut(&Atom) -> Result<Atom>>(&self, mut f: F) -> Result<Self> {
        Ok(Self {
            atoms: self
                .atoms
                .iter()
                .map(|(key, atom)| Ok((*key, f(atom)?)))
                .collect::<Result<_>>()?,
            numerators: self.numerators.clone(),
            scalars: self.scalars.clone(),
            parametric_terms: self
                .parametric_terms
                .iter()
                .map(|(key, terms)| {
                    Ok((
                        *key,
                        terms
                            .iter()
                            .map(|term| {
                                Ok(ParametricIntegrandTerm {
                                    numerator: Arc::clone(&term.numerator),
                                    rows: term
                                        .rows
                                        .iter()
                                        .map(|row| {
                                            Ok(ParametricIntegrandRow {
                                                carrier: f(&row.carrier)?,
                                                ..row.clone()
                                            })
                                        })
                                        .collect::<Result<_>>()?,
                                })
                            })
                            .collect::<Result<_>>()?,
                    ))
                })
                .collect::<Result<_>>()?,
        })
    }

    /// Scalar additive contributions only; consumers must also handle the rows.
    pub fn iter(&self) -> impl Iterator<Item = (&CutCFFIndex, &Atom)> {
        self.atoms.iter()
    }

    /// Inspect stored expressions without specializing residue rows or expanding
    /// shared definitions. Bodies occur once per definition, followed by scalar
    /// roots and row carriers; these are structural pieces, not additive terms.
    pub fn iter_expressions(&self) -> impl Iterator<Item = &Atom> {
        self.numerators
            .iter()
            .chain(&self.scalars)
            .map(|entry| &entry.rhs)
            .chain(self.atoms.values())
            .chain(
                self.parametric_terms
                    .values()
                    .flatten()
                    .flat_map(|term| term.rows.iter().map(|row| &row.carrier)),
            )
    }

    pub(crate) fn atom(&self, index: &CutCFFIndex) -> Option<&Atom> {
        self.atoms.get(index)
    }

    pub(crate) fn parametric_terms(&self, index: &CutCFFIndex) -> &[ParametricIntegrandTerm] {
        self.parametric_terms.get(index).map_or(&[], Vec::as_slice)
    }

    pub(crate) fn cut_indices(&self) -> impl Iterator<Item = &CutCFFIndex> {
        self.atoms.keys()
    }

    pub(crate) fn is_zero(&self) -> bool {
        self.atoms.values().all(Atom::is_zero)
            && self.parametric_terms.values().flatten().all(|term| {
                term.numerator.rhs.is_zero() || term.rows.iter().all(|row| row.carrier.is_zero())
            })
    }

    pub(crate) fn with_parametric_terms(
        mut self,
        terms: impl IntoIterator<Item = (CutCFFIndex, ParametricIntegrandTerm)>,
    ) -> Result<Self> {
        for (index, term) in terms {
            eyre::ensure!(
                term.rows
                    .iter()
                    .all(|row| row.coefficients.len() == term.numerator.args.len()),
                "parametric numerator row does not match its formal coefficients"
            );
            self.atoms.entry(index).or_insert(Atom::Zero);
            self.parametric_terms.entry(index).or_default().push(term);
        }
        let definitions = self
            .numerators
            .iter()
            .cloned()
            .chain(
                self.parametric_terms
                    .values()
                    .flatten()
                    .map(|term| Arc::clone(&term.numerator)),
            )
            .collect::<Vec<_>>();
        self.with_numerators(definitions)
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
                || numerator.lhs != family.call_args(numerator.tags.iter().cloned().chain(parameters.iter().cloned()))
                // Sequential substitution must not capture tags or another key.
                || parameters.iter().enumerate().any(|(index, parameter)| {
                    numerator.tags.iter().any(|tag| tag.contains(parameter))
                        || parameters[..index].iter().any(|other| other.contains(parameter) || parameter.contains(other))
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
        for term in self.parametric_terms.values_mut().flatten() {
            term.numerator =
                Arc::clone(definitions.get(&term.numerator.tags).ok_or_else(|| {
                    eyre!(
                        "missing parametric numerator definition {}",
                        term.numerator.lhs
                    )
                })?);
        }
        for terms in self.parametric_terms.values_mut() {
            let mut families = BTreeMap::<Vec<Atom>, ParametricIntegrandTerm>::new();
            for term in std::mem::take(terms) {
                if let Some(existing) = families.get_mut(&term.numerator.tags) {
                    existing.rows.extend(term.rows);
                } else {
                    families.insert(term.numerator.tags.clone(), term);
                }
            }
            *terms = families.into_values().collect();
        }
        self.numerators = definitions.into_values().collect();
        Ok(self)
    }

    /// Transform each body once, retaining its call signature and all row data.
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

    /// Diagnostic semantic materialization; production preprocessing keeps rows.
    pub(crate) fn resolved(&self) -> Result<Self> {
        let mut resolved = self.resolved_scalars()?;
        for (index, terms) in &resolved.parametric_terms {
            let root = resolved
                .atoms
                .get_mut(index)
                .expect("row cut is registered");
            for term in terms {
                let parameters = term
                    .numerator
                    .args
                    .iter()
                    .cloned()
                    .map(Atom::from)
                    .collect::<Vec<_>>();
                for row in &term.rows {
                    let body = term.numerator.rhs.replace_multiple(
                        parameters
                            .iter()
                            .zip(&row.coefficients)
                            .map(|(parameter, value)| {
                                Replacement::new(parameter.to_pattern(), Atom::num(value.clone()))
                            }),
                    );
                    *root += body * &row.carrier;
                }
            }
        }
        resolved.parametric_terms.clear();
        let replacements = resolved
            .numerators
            .iter()
            .map(|entry| entry.replacement())
            .collect::<Vec<_>>();
        resolved = resolved.map(|atom| atom.replace_multiple(&replacements));
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
        Self {
            atoms: self
                .cut_indices()
                .map(|index| (*index, Atom::Zero))
                .collect(),
            numerators: self.numerators.clone(),
            scalars: self.scalars.clone(),
            parametric_terms: BTreeMap::new(),
        }
    }

    pub fn checked_zip(
        &self,
        other: &Integrands,
        mut map: impl FnMut(&CutCFFIndex, &Atom, &Atom) -> Result<Atom>,
    ) -> Result<Integrands> {
        eyre::ensure!(
            self.parametric_terms.is_empty() && other.parametric_terms.is_empty(),
            "arbitrary scalar zip cannot discard parametric residue rows"
        );
        let atoms = self
            .iter()
            .merge_join_by(other.iter(), |(left, _), (right, _)| left.cmp(right))
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
            parametric_terms: BTreeMap::new(),
        }
        .with_numerators(self.numerators.iter().chain(&other.numerators).cloned())?
        .with_scalar_definitions(self.scalars.iter().chain(&other.scalars).cloned())
    }

    pub fn zip_mul(&self, other: &Integrands) -> Result<Integrands> {
        if let Some(index) = self
            .cut_indices()
            .find(|index| !other.atoms.contains_key(index))
        {
            return Err(eyre!("right integrands are missing key {index:?}"));
        }
        if let Some(index) = other
            .cut_indices()
            .find(|index| !self.atoms.contains_key(index))
        {
            return Err(eyre!("left integrands are missing key {index:?}"));
        }
        let output = Self::from_iter(
            self.atoms
                .iter()
                .map(|(index, left)| (*index, left * &other.atoms[index])),
        );
        let mut terms = Vec::new();
        for index in self.cut_indices() {
            for (source, scalar) in [(self, &other.atoms[index]), (other, &self.atoms[index])] {
                if !scalar.is_zero() {
                    terms.extend(source.parametric_terms(index).iter().map(|term| {
                        (
                            *index,
                            ParametricIntegrandTerm {
                                numerator: Arc::clone(&term.numerator),
                                rows: term
                                    .rows
                                    .iter()
                                    .map(|row| ParametricIntegrandRow {
                                        carrier: &row.carrier * scalar,
                                        ..row.clone()
                                    })
                                    .collect(),
                            },
                        )
                    }));
                }
            }
            for left in self.parametric_terms(index) {
                for right in other.parametric_terms(index) {
                    terms.push((*index, left.product(right)?));
                }
            }
        }
        output
            .with_numerators(self.numerators.iter().chain(&other.numerators).cloned())?
            .with_parametric_terms(terms)?
            .with_scalar_definitions(self.scalars.iter().chain(&other.scalars).cloned())
    }

    /// Validate cut shape, then concatenate additive terms. Products alone form
    /// Cartesian rows; adding forests never takes a Cartesian product.
    pub fn zip_add(self, others: impl IntoIterator<Item = Self>) -> Result<Self> {
        let mut numerators = self.numerators;
        let mut scalars = self.scalars;
        let mut parametric_terms = self.parametric_terms;
        let mut terms = self
            .atoms
            .into_iter()
            .map(|(key, atom)| (key, vec![atom]))
            .collect::<BTreeMap<_, _>>();
        for other in others {
            numerators.extend(other.numerators);
            scalars.extend(other.scalars);
            for (index, rows) in other.parametric_terms {
                parametric_terms.entry(index).or_default().extend(rows);
            }
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
            parametric_terms,
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
            parametric_terms: BTreeMap::new(),
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
