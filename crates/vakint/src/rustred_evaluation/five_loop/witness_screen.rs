//! Exact key support of a verified witness's weighted lowering, without
//! native collection over symbolic coefficients.
//!
//! Weighted lowering is linear over the scalar coefficients. A numerator that
//! is a sum of terms `q * s * P`, with `q` rational, `s` a polynomial or rational
//! function of plain variables free of loop momenta, and `P` a homogeneous
//! polynomial in loop scalar products, is grouped as `sum_(s, r) s * P_(s, r)`,
//! where `P_(s, r)` has rational coefficients and scalar-product degree `r`.
//! Each group then follows the unchanged matcher route (momentum replacement,
//! parent routing) and lowers through RustRed's rational path at unit mass.
//!
//! The family coordinate change is affine with constant terms proportional to
//! the common mass, so a degree-`r` group contributes `m^(r - |b|)` times its
//! unit-mass weight to the key with denominator shift `b`. If the distinct `s`
//! are linearly independent over Q (proven by a nonsingular exact evaluation
//! matrix) and the mass is a nonzero rational or a non-constant function of
//! variables absent from every `s`, a lowered key vanishes exactly if and only
//! if every rational group weight of that key vanishes. The result is the exact
//! key support of the full weighted lowering: its admission and its smallest
//! key decide most witness comparisons. Inputs outside this grammar return
//! `None`, and the caller uses the full lowering instead.
//!
//! Tensor-reduced numerators carry loop-free tensor structures in their
//! coefficients: external vectors with open indices `p(i, hedge(..))`, metrics
//! `g(hedge(..), hedge(..))` and scalar products of external vectors
//! `dot(p(i), p(j))`. Such an atom is never lowered and never looked inside: the
//! full key's zero test of a lowered weight (RustRed's common-factor and
//! rational collection, `collect_factors` and `together`) treats each distinct
//! function atom as an independent indeterminate. The screen therefore treats
//! them as opaque coefficient atoms: linear independence of the distinct `s` is
//! proven over Q with every opaque atom evaluated as its own variable, which is
//! exactly linear independence in the ring of that zero test. Other function
//! atoms (`log`, ...) stay outside the grammar, as does any atom with a loop
//! momentum inside. Tensor-reduced numerators come as products of sums
//! (a sum of coefficient * scalar-product terms times a scalar product and a
//! metric, say); in this tensor grammar such a product is distributed over the
//! scalar products (`expand_in(dot)`, coefficients kept whole) before the
//! grouping, which changes the form of the numerator, not its value. A
//! numerator without opaque atoms keeps the plain grammar unchanged.
//! `VAKINT_WITNESS_SCREEN_OPAQUE=0` (or `off`) restores the plain grammar for
//! every numerator (debug switch for the A/B of the selection).

use std::collections::{BTreeMap, BTreeSet};

use rustred::family::IntegralKey;
use rustred::scalar_numerator::FamilyScalarNumeratorService;
use symbolica::atom::{Atom, AtomCore, AtomView, Symbol};
use symbolica::domains::rational::Rational;
use symbolica::id::Replacement;

use super::super::numerator::{MatchedScalarInput, ScalarNumeratorForm};
use crate::symbols::S;
use crate::{ReplacementRules, VakintError, VakintSettings, VakintTerm};

/// One `P_(s, r)`: rational scalar-product polynomial of degree `degree`
/// multiplying the distinct coefficient `coefficient`.
struct Group {
    coefficient: usize,
    degree: u32,
    part: Atom,
}

/// Witness-independent grouping of one scalar numerator.
pub(super) struct NumeratorGroups {
    coefficients: Vec<Atom>,
    /// Plain variables of all coefficients, for the mass independence check.
    symbols: BTreeSet<Symbol>,
    /// Distinct opaque coefficient atoms, each an independent variable of the
    /// independence proof.
    opaque: Vec<Atom>,
    groups: Vec<Group>,
}

/// Exact weights of one lowered key, by independent slot: distinct
/// coefficient, and group degree when the mass is symbolic.
type KeyWeights = BTreeMap<(usize, u32), Atom>;

/// Why a numerator is outside the exactly decidable grammar (the selection
/// then lowers every witness in full). Diagnostics only.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) enum Refusal {
    /// Screening is off by the settings (identification checks, diagnostics).
    DisabledBySettings,
    /// A factor with loop momenta that is not a homogeneous rational
    /// polynomial in loop scalar products.
    LoopFactor(String),
    /// A number that is not rational.
    NonRationalNumber(String),
    /// A loop-free factor outside the plain-coefficient grammar.
    Coefficient(String),
    /// The distinct coefficients were not proven independent over Q.
    DependentCoefficients(usize),
}

impl std::fmt::Display for Refusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::DisabledBySettings => write!(f, "disabled by settings"),
            Self::LoopFactor(factor) => write!(f, "loop factor outside the grammar: {factor}"),
            Self::NonRationalNumber(number) => write!(f, "non-rational number: {number}"),
            Self::Coefficient(factor) => {
                write!(f, "coefficient factor outside the grammar: {factor}")
            }
            Self::DependentCoefficients(count) => {
                write!(f, "{count} coefficients not proven independent")
            }
        }
    }
}

/// A short excerpt of an atom for a refusal.
fn excerpt(view: AtomView) -> String {
    const LIMIT: usize = 160;
    let text = view.to_string();
    if text.chars().count() <= LIMIT {
        text
    } else {
        text.chars().take(LIMIT).chain("...".chars()).collect()
    }
}

/// Exact admission decision for one witness.
pub(super) enum Screen {
    /// A certainly nonzero lowered key lies outside the saved entry union.
    Rejected,
    /// The complete sorted key support; every key is admitted.
    Admitted(Vec<IntegralKey>),
}

fn rational(view: AtomView) -> Option<Rational> {
    Rational::try_from(view).ok()
}

/// Whether the opaque coefficient atoms are on (debug switch
/// `VAKINT_WITNESS_SCREEN_OPAQUE=0` or `off` turns them off).
pub(super) fn opaque_coefficients_enabled() -> bool {
    !matches!(
        std::env::var("VAKINT_WITNESS_SCREEN_OPAQUE").as_deref(),
        Ok("0") | Ok("off") | Ok("false")
    )
}

/// A factor that belongs to the scalar-product part: with `opaque`, one that
/// contains a loop momentum; otherwise also any scalar product.
fn contains_loop_momentum(view: AtomView, opaque: bool) -> bool {
    let mut found = false;
    view.visitor(&mut |part| {
        if let AtomView::Fun(function) = part
            && (function.get_symbol() == S.k || !opaque && function.get_symbol() == S.dot)
        {
            found = true;
        }
        !found
    });
    found
}

/// A loop-free tensor structure of a tensor-reduced coefficient: an external
/// vector `p(..)`, a metric `g(..)` or a scalar product `dot(..)`.
fn is_opaque_atom(view: AtomView) -> bool {
    matches!(view, AtomView::Fun(function)
        if [S.p, S.g, S.dot].contains(&function.get_symbol()))
        && !contains_loop_momentum(view, true)
}

/// [`is_plain`], with opaque coefficient atoms as further variables.
fn is_plain_or_opaque(view: AtomView) -> bool {
    match view {
        AtomView::Num(_) => rational(view).is_some(),
        AtomView::Var(_) => true,
        AtomView::Fun(_) => is_opaque_atom(view),
        AtomView::Add(sum) => sum.iter().all(is_plain_or_opaque),
        AtomView::Mul(product) => product.iter().all(is_plain_or_opaque),
        AtomView::Pow(power) => {
            rational(power.get_exp()).is_some_and(|exponent| exponent.is_integer())
                && is_plain_or_opaque(power.get_base())
        }
    }
}

/// Whether `view` contains an opaque coefficient atom (a tensor-reduced
/// numerator).
fn has_opaque_atom(view: AtomView) -> bool {
    let mut found = false;
    view.visitor(&mut |part| {
        if is_opaque_atom(part) {
            found = true;
        }
        !found
    });
    found
}

/// Distinct function atoms of a coefficient in the [`is_plain_or_opaque`]
/// grammar (its opaque atoms), in first-appearance order.
fn collect_opaque(view: AtomView, opaque: &mut Vec<Atom>) {
    match view {
        AtomView::Fun(_) => {
            if !opaque.iter().any(|known| known.as_view() == view) {
                opaque.push(view.to_owned());
            }
        }
        AtomView::Add(sum) => sum.iter().for_each(|term| collect_opaque(term, opaque)),
        AtomView::Mul(product) => product
            .iter()
            .for_each(|factor| collect_opaque(factor, opaque)),
        AtomView::Pow(power) => collect_opaque(power.get_base(), opaque),
        AtomView::Num(_) | AtomView::Var(_) => {}
    }
}

/// Replace every opaque atom of a coefficient by its value.
fn substitute_opaque(view: AtomView, opaque: &[Atom], values: &[Atom]) -> Atom {
    match view {
        AtomView::Fun(_) => opaque
            .iter()
            .position(|known| known.as_view() == view)
            .map_or_else(|| view.to_owned(), |index| values[index].clone()),
        AtomView::Add(sum) => Atom::add_many(
            sum.iter()
                .map(|term| substitute_opaque(term, opaque, values))
                .collect::<Vec<_>>(),
        ),
        AtomView::Mul(product) => product.iter().fold(Atom::one(), |value, factor| {
            value * substitute_opaque(factor, opaque, values)
        }),
        AtomView::Pow(power) => {
            substitute_opaque(power.get_base(), opaque, values).pow(power.get_exp().to_owned())
        }
        AtomView::Num(_) | AtomView::Var(_) => view.to_owned(),
    }
}

/// Polynomial or rational function of plain variables with rational numbers
/// and integer exponents. Function atoms and non-rational numbers are refused.
fn is_plain(view: AtomView) -> bool {
    match view {
        AtomView::Num(_) => rational(view).is_some(),
        AtomView::Var(_) => true,
        AtomView::Fun(_) => false,
        AtomView::Add(sum) => sum.iter().all(is_plain),
        AtomView::Mul(product) => product.iter().all(is_plain),
        AtomView::Pow(power) => {
            rational(power.get_exp()).is_some_and(|exponent| exponent.is_integer())
                && is_plain(power.get_base())
        }
    }
}

/// A loop momentum combination: rational multiples and sums of `k(i)`.
fn is_loop_vector(view: AtomView) -> bool {
    match view {
        AtomView::Num(_) => rational(view).is_some(),
        AtomView::Fun(function) => {
            function.get_symbol() == S.k
                && function.get_nargs() == 1
                && function.iter().all(|axis| rational(axis).is_some())
        }
        AtomView::Add(sum) => sum.iter().all(is_loop_vector),
        AtomView::Mul(product) => product.iter().all(is_loop_vector),
        AtomView::Var(_) | AtomView::Pow(_) => false,
    }
}

/// Degree in loop scalar products of a homogeneous factor with rational
/// coefficients, or `None` outside that grammar.
fn scalar_product_degree(view: AtomView) -> Option<u32> {
    match view {
        AtomView::Num(_) => rational(view).map(|_| 0),
        AtomView::Fun(function) => (function.get_symbol() == S.dot
            && function.get_nargs() == 2
            && function.iter().all(is_loop_vector))
        .then_some(1),
        AtomView::Var(_) => None,
        AtomView::Add(sum) => {
            let mut degree = None;
            for term in sum.iter() {
                let term = scalar_product_degree(term)?;
                if degree.is_some_and(|degree| degree != term) {
                    return None;
                }
                degree = Some(term);
            }
            degree
        }
        AtomView::Mul(product) => product.iter().try_fold(0u32, |sum, factor| {
            sum.checked_add(scalar_product_degree(factor)?)
        }),
        AtomView::Pow(power) => {
            let exponent = rational(power.get_exp())?;
            if !exponent.is_integer() || exponent.is_negative() {
                return None;
            }
            let exponent = u32::try_from(exponent.numerator().to_i64()?).ok()?;
            scalar_product_degree(power.get_base())?.checked_mul(exponent)
        }
    }
}

/// Evaluate a coefficient at a rational point (its opaque atoms at `values`),
/// or `None` at a pole.
fn evaluate(
    coefficient: &Atom,
    point: &[Replacement],
    opaque: &[Atom],
    values: &[Atom],
) -> Option<Atom> {
    let value = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        if opaque.is_empty() {
            coefficient.replace_multiple(point)
        } else {
            substitute_opaque(coefficient.as_view(), opaque, values).replace_multiple(point)
        }
    }))
    .ok()?;
    rational(value.as_view()).map(|_| value)
}

/// Exact rank test over Q: Gaussian elimination on normalized numeric atoms.
fn full_row_rank(mut rows: Vec<Vec<Atom>>) -> bool {
    let n = rows.len();
    let columns = rows.first().map_or(0, Vec::len);
    let mut rank = 0;
    for column in 0..columns {
        if rank == n {
            break;
        }
        let Some(pivot) = (rank..n).find(|row| !rows[*row][column].is_zero()) else {
            continue;
        };
        rows.swap(rank, pivot);
        let pivot = rows[rank].clone();
        for row in rows.iter_mut().skip(rank + 1) {
            if row[column].is_zero() {
                continue;
            }
            let factor = &row[column] / &pivot[column];
            for (entry, value) in row.iter_mut().zip(&pivot).skip(column) {
                *entry = &*entry - &factor * value;
            }
        }
        rank += 1;
    }
    rank == n
}

impl NumeratorGroups {
    /// Group a scalar-product-form numerator, or `None` if it is outside the
    /// exactly decidable grammar described in the module documentation.
    #[cfg(test)]
    pub(super) fn new(numerator: AtomView) -> Option<Self> {
        Self::try_new(numerator).ok()
    }

    /// [`Self::new`], with the reason of a refusal; opaque coefficient atoms
    /// as set by the debug switch ([`opaque_coefficients_enabled`]).
    pub(super) fn try_new(numerator: AtomView) -> Result<Self, Refusal> {
        Self::try_new_with(numerator, opaque_coefficients_enabled())
    }

    /// [`Self::try_new`]; with `opaque`, a tensor-reduced numerator (one
    /// with an opaque atom) uses the tensor grammar of the module
    /// documentation: opaque coefficient atoms, and products of sums
    /// distributed over the scalar products. Every other numerator, and every
    /// numerator without `opaque`, uses the plain grammar unchanged.
    pub(super) fn try_new_with(numerator: AtomView, opaque: bool) -> Result<Self, Refusal> {
        let tensor = opaque && has_opaque_atom(numerator);
        match Self::group(numerator, tensor) {
            // A product of sums (a sum of coefficient * scalar-product terms
            // times further factors) is distributed over the scalar products;
            // the coefficients stay whole.
            Err(Refusal::LoopFactor(_)) if tensor => {
                let distributed = numerator.to_owned().expand_in(Atom::var(S.dot));
                Self::group(distributed.as_view(), true)
            }
            grouped => grouped,
        }
    }

    /// The grouping proper, in the tensor grammar when `opaque`.
    fn group(numerator: AtomView, opaque: bool) -> Result<Self, Refusal> {
        let mut coefficients = Vec::<Atom>::new();
        let mut parts = BTreeMap::<(usize, u32), Vec<Atom>>::new();
        let terms = match numerator {
            AtomView::Add(sum) => sum.iter().collect::<Vec<_>>(),
            term => vec![term],
        };
        for term in terms {
            let factors = match term {
                AtomView::Mul(product) => product.iter().collect::<Vec<_>>(),
                factor => vec![factor],
            };
            let mut number = Atom::one();
            let mut coefficient = Atom::one();
            let mut scalar_products = Atom::one();
            let mut degree = 0u32;
            for factor in factors {
                if contains_loop_momentum(factor, opaque) {
                    degree = scalar_product_degree(factor)
                        .and_then(|factor| degree.checked_add(factor))
                        .ok_or_else(|| Refusal::LoopFactor(excerpt(factor)))?;
                    scalar_products *= factor.to_owned();
                } else if let AtomView::Num(_) = factor {
                    rational(factor).ok_or_else(|| Refusal::NonRationalNumber(excerpt(factor)))?;
                    number *= factor.to_owned();
                } else if is_plain(factor) || opaque && is_plain_or_opaque(factor) {
                    coefficient *= factor.to_owned();
                } else {
                    return Err(Refusal::Coefficient(excerpt(factor)));
                }
            }
            if coefficient.is_zero() || number.is_zero() {
                continue;
            }
            let index = match coefficients.iter().position(|known| *known == coefficient) {
                Some(index) => index,
                None => {
                    coefficients.push(coefficient);
                    coefficients.len() - 1
                }
            };
            parts
                .entry((index, degree))
                .or_default()
                .push(number * scalar_products);
        }
        let symbols = coefficients
            .iter()
            .flat_map(|coefficient| coefficient.get_all_symbols(true))
            .collect::<BTreeSet<_>>();
        let mut opaque_atoms = Vec::new();
        for coefficient in &coefficients {
            collect_opaque(coefficient.as_view(), &mut opaque_atoms);
        }
        let groups = parts
            .into_iter()
            .map(|((coefficient, degree), terms)| Group {
                coefficient,
                degree,
                part: Atom::add_many(terms),
            })
            .collect();
        let decomposition = Self {
            coefficients,
            symbols,
            opaque: opaque_atoms,
            groups,
        };
        if decomposition.coefficients_are_independent() {
            Ok(decomposition)
        } else {
            Err(Refusal::DependentCoefficients(
                decomposition.coefficients.len(),
            ))
        }
    }

    /// Distinct coefficients and groups, for diagnostics.
    pub(super) fn shape(&self) -> (usize, usize) {
        (self.coefficients.len(), self.groups.len())
    }

    /// Number of distinct opaque coefficient atoms, for diagnostics.
    pub(super) fn opaque_atoms(&self) -> usize {
        self.opaque.len()
    }

    /// Prove linear independence of the distinct coefficients over Q with a
    /// nonsingular matrix of exact values at deterministic rational points.
    fn coefficients_are_independent(&self) -> bool {
        let count = self.coefficients.len();
        if count <= 1 {
            return true;
        }
        let symbols = self.symbols.iter().copied().collect::<Vec<_>>();
        // Two extra points tolerate an unlucky singular square sample.
        let points = (0..count + 2)
            .map(|point| {
                symbols
                    .iter()
                    .enumerate()
                    .map(|(index, symbol)| {
                        // A large prime denominator avoids the small rational
                        // poles of physical coefficients. A pole or an unlucky
                        // singular sample only fails the proof (full lowering).
                        let numerator = 7919 * (index as i64 + 1)
                            + 104_729 * (point as i64 + 1)
                            + (index as i64 + 3) * (point as i64 + 5) * (point as i64 + 7);
                        Replacement::new(
                            Atom::var(*symbol),
                            Atom::num(numerator) / Atom::num(1_000_003),
                        )
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        // Each opaque atom is a further independent variable.
        let opaque_values = (0..count + 2)
            .map(|point| {
                (0..self.opaque.len())
                    .map(|atom| {
                        let index = (symbols.len() + atom) as i64;
                        let numerator = 7919 * (index + 1)
                            + 104_729 * (point as i64 + 1)
                            + (index + 3) * (point as i64 + 5) * (point as i64 + 7);
                        Atom::num(numerator) / Atom::num(1_000_003)
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut rows = Vec::with_capacity(count);
        for coefficient in &self.coefficients {
            let mut row = Vec::with_capacity(points.len());
            for (point, values) in points.iter().zip(&opaque_values) {
                let Some(value) = evaluate(coefficient, point, &self.opaque, values) else {
                    return false;
                };
                row.push(value);
            }
            rows.push(row);
        }
        full_row_rank(rows)
    }

    /// Decide the admission of `witness` from the exact key support of the
    /// full weighted lowering of the numerator through it: `Rejected` as soon
    /// as a key that is certainly nonzero lies outside the saved entry union,
    /// otherwise `Admitted` with the complete sorted support. `None` if this
    /// route cannot decide it (unsupported mass, or any refusal or error of the
    /// unit-mass lowering); the caller then runs the unchanged full lowering,
    /// which reproduces any genuine error. Admission errors propagate.
    ///
    /// Groups are ordered by coefficient and then degree, so each independent
    /// weight slot (a group for a symbolic mass, a coefficient for a numeric
    /// mass) is complete before the next one starts. A key with a nonzero
    /// complete slot weight is certainly in the support, which permits the
    /// early rejection.
    pub(super) fn screen(
        &self,
        settings: &VakintSettings,
        witness: &ReplacementRules,
        scalar: &FamilyScalarNumeratorService<'_>,
        admitted: &(dyn Fn(&IntegralKey) -> Result<bool, VakintError> + Sync),
    ) -> Result<Option<Screen>, VakintError> {
        let mut matched = witness.clone();
        if matched.apply_replacement_rules().is_err() {
            return Ok(None);
        }
        let Ok(descriptor) = super::select(&matched.canonical_topology) else {
            return Ok(None);
        };
        // Numeric masses scale each group weight; a symbolic mass separates
        // group degrees instead and must be independent of all coefficients.
        let mut numeric_mass = None;
        let mut weights = BTreeMap::<IntegralKey, KeyWeights>::new();
        let mut checked = BTreeSet::<IntegralKey>::new();
        for (index, group) in self.groups.iter().enumerate() {
            let Some(lowered) = self.lower_group(
                settings,
                &matched,
                descriptor,
                group,
                scalar,
                &mut numeric_mass,
            ) else {
                return Ok(None);
            };
            let numeric = numeric_mass.as_ref().expect("mass mode is set");
            for (key, weight) in lowered {
                let slot = match numeric {
                    Some(_) => (group.coefficient, 0),
                    None => (group.coefficient, group.degree),
                };
                let total = weights
                    .entry(key)
                    .or_default()
                    .entry(slot)
                    .or_insert(Atom::Zero);
                *total = &*total + &weight;
            }
            let slot_complete = numeric.is_none()
                || self
                    .groups
                    .get(index + 1)
                    .is_none_or(|next| next.coefficient != group.coefficient);
            if !slot_complete {
                continue;
            }
            let slot = match numeric {
                Some(_) => (group.coefficient, 0),
                None => (group.coefficient, group.degree),
            };
            for (key, weights) in &weights {
                if !checked.contains(key)
                    && weights.get(&slot).is_some_and(|weight| !weight.is_zero())
                {
                    if !admitted(key)? {
                        return Ok(Some(Screen::Rejected));
                    }
                    checked.insert(key.clone());
                }
            }
        }
        // Every key with a nonzero slot was checked when its slot completed.
        Ok(Some(Screen::Admitted(
            weights
                .into_iter()
                .filter(|(_, weights)| weights.values().any(|weight| !weight.is_zero()))
                .map(|(key, _)| key)
                .collect(),
        )))
    }

    /// Lower one group through the unchanged matcher route at unit mass and
    /// restore its exact rational weights (times the mass power for a numeric
    /// mass). Sets the mass mode on the first group.
    fn lower_group(
        &self,
        settings: &VakintSettings,
        matched: &ReplacementRules,
        descriptor: &super::ParentDescriptor,
        group: &Group,
        scalar: &FamilyScalarNumeratorService<'_>,
        numeric_mass: &mut Option<Option<Atom>>,
    ) -> Option<Vec<(IntegralKey, Atom)>> {
        let mut term = VakintTerm {
            integral: Atom::one(),
            numerator: group.part.clone(),
            vectors: Vec::new(),
        };
        term.apply_numerator_replacement_rules(matched, settings)
            .ok()?;
        let input = MatchedScalarInput::prepare(
            term.numerator.as_view(),
            matched,
            &descriptor.momenta,
            descriptor.physical_slots,
            15,
            settings,
            ScalarNumeratorForm::PreserveProducts,
        )
        .ok()?;
        if numeric_mass.is_none() {
            *numeric_mass = Some(self.mass_mode(&input.mass_squared)?);
        }
        let mass = numeric_mass.as_ref().expect("mass mode is set");
        let lowered = scalar
            .lower_weighted(&input.numerator, &input.target, &Atom::one())
            .ok()?;
        lowered
            .terms()
            .iter()
            .map(|lowered| {
                let weight = lowered.coefficient().to_expression() * lowered.scalar_spectator();
                rational(weight.as_view())?;
                let shift = input
                    .target
                    .powers()
                    .iter()
                    .zip(lowered.integral().powers())
                    .try_fold(0i64, |sum, (target, key)| {
                        sum.checked_add(target.checked_sub(*key)?)
                    })?;
                let mass_power = u32::try_from(i64::from(group.degree).checked_sub(shift)?).ok()?;
                let weight = match mass {
                    Some(mass) => weight * mass.pow(Atom::num(i64::from(mass_power))),
                    None => weight,
                };
                Some((lowered.integral().clone(), weight))
            })
            .collect()
    }

    /// The exact sorted key support, without admission (for tests).
    #[cfg(test)]
    fn support(
        &self,
        settings: &VakintSettings,
        witness: &ReplacementRules,
        scalar: &FamilyScalarNumeratorService<'_>,
    ) -> Option<Vec<IntegralKey>> {
        match self.screen(settings, witness, scalar, &|_| Ok(true)) {
            Ok(Some(Screen::Admitted(support))) => Some(support),
            _ => None,
        }
    }

    /// `Some(Some(m))` for a nonzero rational mass, `Some(None)` for a plain
    /// non-constant mass in variables absent from every coefficient.
    fn mass_mode(&self, mass: &Atom) -> Option<Option<Atom>> {
        if let AtomView::Num(_) = mass.as_view() {
            return (rational(mass.as_view()).is_some() && !mass.is_zero())
                .then(|| Some(mass.clone()));
        }
        let symbols = mass.get_all_symbols(true);
        (is_plain(mass.as_view())
            && !symbols.is_empty()
            && symbols.iter().all(|symbol| !self.symbols.contains(symbol)))
        .then_some(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Vakint;
    use crate::utils::vakint_macros::vk_parse;

    #[test]
    fn groups_split_coefficients_and_scalar_product_degrees() {
        Vakint::initialize_vakint_symbols();
        let numerator = vk_parse!(
            "2/3*gs^10*(4-2*ε)*dot(k(1),k(2))^2 + gs^10*(4-2*ε)*(dot(k(1),k(1))-3*dot(k(2),k(3)))*dot(k(4),k(4)) + 5*gs^10*dot(k(5),k(5))/(3-2*ε) - dot(k(1),k(1))"
        )
        .unwrap();
        let groups = NumeratorGroups::new(numerator.as_view()).unwrap();
        // gs^10*(4-2ε), gs^10/(3-2ε) and 1; degrees 2, 1 and 1.
        assert_eq!(groups.shape(), (3, 3));
        let rebuilt = groups.groups.iter().fold(Atom::Zero, |sum, group| {
            sum + &groups.coefficients[group.coefficient] * &group.part
        });
        assert!((rebuilt - numerator).expand().is_zero());
    }

    #[test]
    fn refuses_inhomogeneous_function_and_dependent_coefficients() {
        Vakint::initialize_vakint_symbols();
        for refused in [
            // Mixed scalar-product degree inside one factor.
            "(dot(k(1),k(1))+1)*dot(k(2),k(2))",
            // Function coefficient: no exact independence proof by evaluation.
            "log(x)*dot(k(1),k(1))",
            // x*(1+x) = x + x^2 is a rational combination of x and x^2.
            "x*dot(k(1),k(1)) + x^2*dot(k(2),k(2)) + x*(1+x)*dot(k(3),k(3))",
            // Non-rational number.
            "1.5*dot(k(1),k(1))",
        ] {
            let numerator = vk_parse!(refused).unwrap();
            assert!(
                NumeratorGroups::new(numerator.as_view()).is_none(),
                "{refused} must use the full lowering"
            );
        }
        let independent =
            vk_parse!("x*dot(k(1),k(1)) + x^2*dot(k(2),k(2)) + (1+x)^-1*dot(k(3),k(3))").unwrap();
        assert_eq!(
            NumeratorGroups::new(independent.as_view()).unwrap().shape(),
            (3, 3)
        );
    }

    /// Tensor-reduced coefficients: external vectors with open indices,
    /// metrics and external scalar products are opaque coefficient atoms;
    /// without them (the debug switch off) the numerator is refused as before.
    #[test]
    fn opaque_tensor_coefficients_group_and_rebuild() {
        Vakint::initialize_vakint_symbols();
        let numerator = vk_parse!(
            "gs^9*(4-2*ε)^-1*g(hedge(10,1),hedge(23,1))*dot(k(1),k(3))^2
            + 2*gs^9*p(6,hedge(10,1))*p(6,hedge(23,1))*dot(k(1),k(1))*dot(k(3),k(3))
            - 3*gs^9*dot(p(6),p(6))*g(hedge(10,1),hedge(23,1))*dot(k(2),k(4))
            + gs^9*(5-2*ε)*dot(p(6),p(6))*g(hedge(10,1),hedge(23,1))*dot(k(2),k(2))"
        )
        .unwrap();
        let groups = NumeratorGroups::try_new_with(numerator.as_view(), true).unwrap();
        // gs^9*g/(4-2ε), gs^9*p*p, gs^9*dot(p,p)*g and gs^9*(5-2ε)*dot(p,p)*g;
        // opaque atoms g(..), p(6,hedge(10,1)), p(6,hedge(23,1)), dot(p(6),p(6)).
        assert_eq!(groups.shape(), (4, 4));
        assert_eq!(groups.opaque_atoms(), 4);
        let rebuilt = groups.groups.iter().fold(Atom::Zero, |sum, group| {
            sum + &groups.coefficients[group.coefficient] * &group.part
        });
        assert!((rebuilt - &numerator).expand().is_zero());
        assert!(matches!(
            NumeratorGroups::try_new_with(numerator.as_view(), false),
            Err(Refusal::Coefficient(_) | Refusal::LoopFactor(_))
        ));
    }

    /// Opaque atoms are only loop-free p, g and dot atoms; the grammar's
    /// other refusals are unchanged.
    #[test]
    fn opaque_coefficients_keep_loop_momenta_and_other_functions_out() {
        Vakint::initialize_vakint_symbols();
        for (refused, why) in [
            (
                "log(x)*g(hedge(1,1),hedge(2,1))*dot(k(1),k(1))",
                "a function coefficient other than p, g and dot",
            ),
            (
                "p(6,k(1))*dot(k(2),k(2))",
                "a loop momentum in an external slot",
            ),
            ("dot(p(6),k(1))*dot(k(2),k(2))", "a mixed scalar product"),
            (
                "g(hedge(1,1),hedge(2,1))*(dot(k(1),k(1))+log(x))",
                "a function coefficient other than p, g and dot after distribution",
            ),
            (
                "x*g(hedge(1,1),hedge(2,1))*dot(k(1),k(1)) + x^2*g(hedge(1,1),hedge(2,1))*dot(k(2),k(2)) + x*(1+x)*g(hedge(1,1),hedge(2,1))*dot(k(3),k(3))",
                "dependent coefficients",
            ),
        ] {
            let numerator = vk_parse!(refused).unwrap();
            assert!(
                NumeratorGroups::try_new_with(numerator.as_view(), true).is_err(),
                "{refused} ({why}) must use the full lowering"
            );
        }
        for refused in [
            "(dot(k(1),k(1))+1)*dot(k(2),k(2))",
            "log(x)*dot(k(1),k(1))",
            "x*dot(k(1),k(1)) + x^2*dot(k(2),k(2)) + x*(1+x)*dot(k(3),k(3))",
            "1.5*dot(k(1),k(1))",
        ] {
            let numerator = vk_parse!(refused).unwrap();
            assert!(NumeratorGroups::try_new_with(numerator.as_view(), true).is_err());
        }
        // Distinct opaque monomials are independent; equal ones share a slot.
        let accepted = vk_parse!(
            "g(hedge(1,1),hedge(2,1))*dot(k(1),k(1)) + p(6,hedge(1,1))*p(6,hedge(2,1))*dot(k(2),k(2)) + dot(p(6),p(6))*g(hedge(1,1),hedge(2,1))*dot(k(3),k(3)) + 3*g(hedge(2,1),hedge(1,1))*dot(k(4),k(4))"
        )
        .unwrap();
        let groups = NumeratorGroups::try_new_with(accepted.as_view(), true).unwrap();
        assert_eq!(groups.shape(), (3, 3));
    }

    /// A tensor-reduced product of sums is distributed over the scalar
    /// products; the same product without tensor atoms keeps the plain
    /// grammar's refusal.
    #[test]
    fn tensor_products_of_sums_are_distributed() {
        Vakint::initialize_vakint_symbols();
        let tensor = vk_parse!(
            "(gs*ε*dot(k(1),k(1)) + 3*gs*dot(k(2),k(2)) + gs*ε*dot(k(2),k(4)))*(4-2*ε)^-1*dot(k(1),k(3))*g(hedge(10,1),hedge(23,1))"
        )
        .unwrap();
        let groups = NumeratorGroups::try_new_with(tensor.as_view(), true).unwrap();
        // gs*ε*g/(4-2ε) and gs*g/(4-2ε), both of degree 2.
        assert_eq!(groups.shape(), (2, 2));
        let rebuilt = groups.groups.iter().fold(Atom::Zero, |sum, group| {
            sum + &groups.coefficients[group.coefficient] * &group.part
        });
        assert!((rebuilt - &tensor).expand().is_zero());
        let plain = vk_parse!(
            "(gs*ε*dot(k(1),k(1)) + 3*gs*dot(k(2),k(2)) + gs*ε*dot(k(2),k(4)))*(4-2*ε)^-1*dot(k(1),k(3))"
        )
        .unwrap();
        assert!(matches!(
            NumeratorGroups::try_new_with(plain.as_view(), true),
            Err(Refusal::LoopFactor(_))
        ));
        // Mixed degrees and an external scalar product inside one loop factor
        // become separate homogeneous groups.
        let mixed =
            vk_parse!("g(hedge(1,1),hedge(2,1))*(dot(k(1),k(1))+dot(p(6),p(6)))*dot(k(2),k(2))")
                .unwrap();
        assert_eq!(
            NumeratorGroups::try_new_with(mixed.as_view(), true)
                .unwrap()
                .shape(),
            (2, 2)
        );
    }

    /// The grouped support of every verified witness of the ladder equals the
    /// key set of the full weighted lowering, for coefficients with tensor
    /// structures (and cancellations inside each coefficient class).
    #[test]
    fn grouped_support_equals_full_weighted_lowering_support_with_tensor_coefficients() {
        Vakint::initialize_vakint_symbols();
        assert_grouped_support_equals_full_weighted_lowering_support(
            "256/81*gs^10*g(hedge(1,1),hedge(2,1))*(4-2*ε)*(dot(k(2),k(3))-dot(k(1),k(2)))*(dot(k(3),k(4))-dot(k(1),k(4)))*dot(k(4),k(4))*dot(k(5),k(5))
            -32/81*gs^10*p(6,hedge(1,1))*p(6,hedge(2,1))*(4-2*ε)^2/((3-2*ε)*(6-2*ε))*(dot(k(1),k(1))-2*dot(k(1),k(4))+dot(k(4),k(4)))^2*dot(k(2),k(5))
            +32/81*gs^10*p(6,hedge(1,1))*p(6,hedge(2,1))*(4-2*ε)^2/((3-2*ε)*(6-2*ε))*dot(k(1),k(1))^2*dot(k(2),k(5))
            +5*gs^10*dot(p(6),p(6))*g(hedge(1,1),hedge(2,1))*(dot(k(1),k(5))-dot(k(3),k(5)))*dot(k(2),k(4))*dot(k(4),k(4))*dot(k(1),k(2))
            -5*gs^10*dot(p(6),p(6))*g(hedge(1,1),hedge(2,1))*dot(k(1),k(5))*dot(k(2),k(4))*dot(k(4),k(4))*dot(k(1),k(2))
            +7/3*gs^10*g(hedge(1,1),hedge(2,1))*(4-2*ε)*dot(k(1),k(2))^2*dot(k(3),k(5))^2",
            Some((3, 3)),
        );
        // A product of sums, as tensor reduction writes it (distributed), with
        // a cancellation that appears only after the distribution.
        assert_grouped_support_equals_full_weighted_lowering_support(
            "(gs^10*ε*dot(k(1),k(1))*dot(k(4),k(4)) - 3*gs^10*(5-2*ε)*dot(k(2),k(5))*dot(k(3),k(3)) + 2*gs^10*dot(k(1),k(2))^2 + gs^10*ε^2*dot(k(3),k(5))*dot(k(1),k(4)))*(4-2*ε)^-1*dot(k(1),k(3))*g(hedge(10,1),hedge(23,1))
            + (7*dot(k(2),k(2))*dot(k(4),k(5)) + ε*dot(k(1),k(5))^2)*gs^10*dot(p(6),p(6))*p(6,hedge(10,1))*p(6,hedge(23,1))*dot(k(2),k(3))
            - 2*gs^10*dot(k(1),k(2))^2*(4-2*ε)^-1*dot(k(1),k(3))*g(hedge(10,1),hedge(23,1))",
            None,
        );
    }

    /// The ladder of [`grouped_support_equals_full_weighted_lowering_support`]
    /// with `numerator`: grouped (opaque atoms on, expected `shape` if given) and full
    /// supports agree on every verified witness.
    fn assert_grouped_support_equals_full_weighted_lowering_support(
        numerator: &str,
        shape: Option<(usize, usize)>,
    ) {
        let settings = VakintSettings::default();
        let topologies = crate::Topologies::generate_topologies().unwrap();
        topologies
            .install_five_loop_catalogue(crate::topologies::FiveLoopCatalogue::prepare().unwrap());
        let family = super::super::common_basis_family("tide_five_loop_common_basis").unwrap();
        let scalar = FamilyScalarNumeratorService::try_new(
            &family,
            S.dot,
            (1..=5).map(|axis| S.k.call(axis)).collect(),
            settings.rustred_scalar_numerator_limits,
        )
        .unwrap();
        let input = vk_parse!(
            "topo(prop(0,edge(1,7),k(1),mUV^2,2)*prop(1,edge(3,5),k(2),mUV^2,2)*prop(10,edge(3,5),-k(2)-k(3)+k(4)+k(5),mUV^2,1)*prop(11,edge(1,7),-k(1)+k(4),mUV^2,1)*prop(2,edge(0,1),k(3),mUV^2,1)*prop(3,edge(0,3),-k(3)+k(4),mUV^2,2)*prop(4,edge(7,0),k(4),mUV^2,2)*prop(7,edge(1,3),k(5),mUV^2,1)*prop(8,edge(5,1),-k(3)+k(4)+k(5),mUV^2,2))"
        )
        .unwrap();
        let numerator = vk_parse!(numerator).unwrap();
        let scalar_form = VakintTerm::scalar_product_form(&settings, numerator.as_view()).unwrap();
        let groups = NumeratorGroups::try_new_with(scalar_form.as_view(), true).unwrap();
        if let Some(shape) = shape {
            assert_eq!(groups.shape(), shape);
        }
        let mut compared = 0;
        topologies
            .match_topologies_to_user_input_with_selection(
                input.as_view(),
                false,
                Default::default(),
                |witnesses| {
                    for witness in witnesses {
                        let grouped = groups.support(&settings, witness, &scalar).unwrap();
                        let mut matched = witness.clone();
                        matched.apply_replacement_rules().unwrap();
                        let mut term = VakintTerm {
                            integral: Atom::one(),
                            numerator: numerator.clone(),
                            vectors: Vec::new(),
                        };
                        term.apply_scalar_numerator_replacement_rules(&matched, &settings)
                            .unwrap();
                        let descriptor = super::super::select(&matched.canonical_topology).unwrap();
                        let full = MatchedScalarInput::prepare(
                            term.numerator.as_view(),
                            &matched,
                            &descriptor.momenta,
                            descriptor.physical_slots,
                            15,
                            &settings,
                            ScalarNumeratorForm::PreserveProducts,
                        )
                        .unwrap();
                        let lowered = scalar
                            .lower_weighted(&full.numerator, &full.target, &full.mass_squared)
                            .unwrap();
                        let expected = lowered
                            .terms()
                            .iter()
                            .map(|term| term.integral().clone())
                            .collect::<Vec<_>>();
                        assert!(!expected.is_empty());
                        assert_eq!(grouped, expected);
                        compared += 1;
                    }
                    Ok(None)
                },
            )
            .unwrap()
            .unwrap();
        assert!(compared > 1, "the ladder has several verified witnesses");
    }

    /// Beta's property on synthetic tensor-reduced numerators with the
    /// official package: the screened selection equals the first minimal full
    /// key (terminal-only first), also for a master input with a tensor
    /// coefficient (its support must stay that master).
    #[test]
    #[ignore = "loads the official five-loop saved application (minutes)"]
    fn screened_selection_matches_first_minimal_full_key_with_tensor_coefficients() {
        let directory = super::super::super::test_data::external_path_through(5).join("5L/ibps");
        let mut settings = VakintSettings::default();
        let engine = super::super::super::test_data::external_engine()
            .with_rustred_application(&mut settings, &directory)
            .unwrap();
        let input = vk_parse!(
            "topo(prop(0,edge(1,7),k(1),mUV^2,2)*prop(1,edge(3,5),k(2),mUV^2,2)*prop(10,edge(3,5),-k(2)-k(3)+k(4)+k(5),mUV^2,1)*prop(11,edge(1,7),-k(1)+k(4),mUV^2,1)*prop(2,edge(0,1),k(3),mUV^2,1)*prop(3,edge(0,3),-k(3)+k(4),mUV^2,2)*prop(4,edge(7,0),k(4),mUV^2,2)*prop(7,edge(1,3),k(5),mUV^2,1)*prop(8,edge(5,1),-k(3)+k(4)+k(5),mUV^2,2))"
        )
        .unwrap();
        for numerator in [
            "g(hedge(1,1),hedge(2,1))*gs^10*(4-2*ε)*dot(k(1),k(1))",
            "p(6,hedge(1,1))*p(6,hedge(2,1))*gs^10*dot(k(2),k(2))-dot(p(6),p(6))*g(hedge(1,1),hedge(2,1))*gs^10/(3-2*ε)*dot(k(1),k(2))",
            "256/81*gs^10*g(hedge(1,1),hedge(2,1))*(4-2*ε)*(dot(k(2),k(3))-dot(k(1),k(2)))*(dot(k(3),k(4))-dot(k(1),k(4)))*dot(k(4),k(4))*dot(k(5),k(5))+7/3*gs^10*p(6,hedge(1,1))*p(6,hedge(2,1))/(3-2*ε)*dot(k(1),k(2))^2*dot(k(3),k(5))^2",
        ] {
            let numerator = vk_parse!(numerator).unwrap();
            let scalar_form =
                VakintTerm::scalar_product_form(&settings, numerator.as_view()).unwrap();
            assert!(NumeratorGroups::try_new_with(scalar_form.as_view(), true).is_ok());
            assert!(screened_selection_matches_reference(
                &engine,
                &settings,
                input.as_view(),
                numerator.as_view(),
                None,
            ));
        }
        let powers = [1, 1, 0, 1, 2, -1, 1, 1, 1, 1, 0, 1, 0, 0, 1];
        let key = IntegralKey::try_new(powers.to_vec()).unwrap();
        let master = super::super::admission_census_tests::input_for_key(&key).unwrap();
        let numerator =
            &master.numerator * &vk_parse!("g(hedge(1,1),hedge(2,1))*p(6,hedge(3,1))").unwrap();
        assert!(screened_selection_matches_reference(
            &engine,
            &settings,
            master.integral.as_view(),
            numerator.as_view(),
            Some(&key),
        ));
    }

    /// Every verified witness of a five-loop ladder input: the grouped
    /// support equals the key set of the unchanged full weighted lowering.
    #[test]
    fn grouped_support_equals_full_weighted_lowering_support() {
        Vakint::initialize_vakint_symbols();
        let settings = VakintSettings::default();
        let topologies = crate::Topologies::generate_topologies().unwrap();
        topologies
            .install_five_loop_catalogue(crate::topologies::FiveLoopCatalogue::prepare().unwrap());
        let family = super::super::common_basis_family("tide_five_loop_common_basis").unwrap();
        let scalar = FamilyScalarNumeratorService::try_new(
            &family,
            S.dot,
            (1..=5).map(|axis| S.k.call(axis)).collect(),
            settings.rustred_scalar_numerator_limits,
        )
        .unwrap();
        let input = vk_parse!(
            "topo(prop(0,edge(1,7),k(1),mUV^2,2)*prop(1,edge(3,5),k(2),mUV^2,2)*prop(10,edge(3,5),-k(2)-k(3)+k(4)+k(5),mUV^2,1)*prop(11,edge(1,7),-k(1)+k(4),mUV^2,1)*prop(2,edge(0,1),k(3),mUV^2,1)*prop(3,edge(0,3),-k(3)+k(4),mUV^2,2)*prop(4,edge(7,0),k(4),mUV^2,2)*prop(7,edge(1,3),k(5),mUV^2,1)*prop(8,edge(5,1),-k(3)+k(4)+k(5),mUV^2,2))"
        )
        .unwrap();
        // Coefficients and degrees as in the projected ladder numerators,
        // with deliberate cancellations between coefficient classes.
        let numerator = vk_parse!(
            "256/81*gs^10*(4-2*ε)*(dot(k(2),k(3))-dot(k(1),k(2)))*(dot(k(3),k(4))-dot(k(1),k(4)))*dot(k(4),k(4))*dot(k(5),k(5))
            -32/81*gs^10*(4-2*ε)^2/((3-2*ε)*(6-2*ε))*(dot(k(1),k(1))-2*dot(k(1),k(4))+dot(k(4),k(4)))^2*dot(k(2),k(5))
            +32/81*gs^10*(4-2*ε)^2/((3-2*ε)*(6-2*ε))*dot(k(1),k(1))^2*dot(k(2),k(5))
            +5*gs^10*(dot(k(1),k(5))-dot(k(3),k(5)))*dot(k(2),k(4))*dot(k(4),k(4))*dot(k(1),k(2))
            -5*gs^10*dot(k(1),k(5))*dot(k(2),k(4))*dot(k(4),k(4))*dot(k(1),k(2))
            +7/3*gs^10*(4-2*ε)*dot(k(1),k(2))^2*dot(k(3),k(5))^2"
        )
        .unwrap();
        let scalar_form = VakintTerm::scalar_product_form(&settings, numerator.as_view()).unwrap();
        let groups = NumeratorGroups::new(scalar_form.as_view()).unwrap();
        assert_eq!(groups.shape(), (3, 3));
        let mut compared = 0;
        topologies
            .match_topologies_to_user_input_with_selection(
                input.as_view(),
                false,
                Default::default(),
                |witnesses| {
                    for witness in witnesses {
                        let grouped = groups.support(&settings, witness, &scalar).unwrap();
                        let mut matched = witness.clone();
                        matched.apply_replacement_rules().unwrap();
                        let mut term = VakintTerm {
                            integral: Atom::one(),
                            numerator: numerator.clone(),
                            vectors: Vec::new(),
                        };
                        term.apply_scalar_numerator_replacement_rules(&matched, &settings)
                            .unwrap();
                        let descriptor = super::super::select(&matched.canonical_topology).unwrap();
                        let full = MatchedScalarInput::prepare(
                            term.numerator.as_view(),
                            &matched,
                            &descriptor.momenta,
                            descriptor.physical_slots,
                            15,
                            &settings,
                            ScalarNumeratorForm::PreserveProducts,
                        )
                        .unwrap();
                        let lowered = scalar
                            .lower_weighted(&full.numerator, &full.target, &full.mass_squared)
                            .unwrap();
                        let expected = lowered
                            .terms()
                            .iter()
                            .map(|term| term.integral().clone())
                            .collect::<Vec<_>>();
                        assert!(!expected.is_empty());
                        assert_eq!(grouped, expected);
                        compared += 1;
                    }
                    Ok(None)
                },
            )
            .unwrap()
            .unwrap();
        assert!(compared > 1, "the ladder has several verified witnesses");
    }

    /// The exact selection equals the first minimal full key over all verified
    /// witnesses, with the official five-loop package's real admission, among
    /// the witnesses whose whole support is remaining terminals when there are
    /// any (terminal-first preference), otherwise among all of them.
    #[test]
    #[ignore = "loads the official five-loop saved application (minutes)"]
    fn screened_selection_matches_first_minimal_full_key_with_the_official_package() {
        let directory = super::super::super::test_data::external_path_through(5).join("5L/ibps");
        let mut settings = VakintSettings::default();
        let engine = super::super::super::test_data::external_engine()
            .with_rustred_application(&mut settings, &directory)
            .unwrap();
        let input = vk_parse!(
            "topo(prop(0,edge(1,7),k(1),mUV^2,2)*prop(1,edge(3,5),k(2),mUV^2,2)*prop(10,edge(3,5),-k(2)-k(3)+k(4)+k(5),mUV^2,1)*prop(11,edge(1,7),-k(1)+k(4),mUV^2,1)*prop(2,edge(0,1),k(3),mUV^2,1)*prop(3,edge(0,3),-k(3)+k(4),mUV^2,2)*prop(4,edge(7,0),k(4),mUV^2,2)*prop(7,edge(1,3),k(5),mUV^2,1)*prop(8,edge(5,1),-k(3)+k(4)+k(5),mUV^2,2))"
        )
        .unwrap();
        for numerator in [
            "1",
            "gs^10*(4-2*ε)*dot(k(1),k(1))",
            "256/81*gs^10*(4-2*ε)*(dot(k(2),k(3))-dot(k(1),k(2)))*(dot(k(3),k(4))-dot(k(1),k(4)))*dot(k(4),k(4))*dot(k(5),k(5))+7/3*gs^10/(3-2*ε)*dot(k(1),k(2))^2*dot(k(3),k(5))^2",
        ] {
            let numerator = vk_parse!(numerator).unwrap();
            assert!(screened_selection_matches_reference(
                &engine,
                &settings,
                input.as_view(),
                numerator.as_view(),
                None,
            ));
        }
        // Master inputs: the integrand of a remaining terminal must be sent as
        // that terminal itself, also when a symmetric image has a smaller key
        // (BASIS-v2 P387 = 640-pack P387, [1,1,0,1,2,-1,1,1,1,1,0,1,0,0,1]; its
        // image [-1,1,0,1,1,1,1,1,1,1,0,2,0,0,1] costs a 1.57 M-key walk).
        for powers in [
            [1, 1, 0, 1, 2, -1, 1, 1, 1, 1, 0, 1, 0, 0, 1],
            [1, 1, 0, 1, 1, -1, 1, 1, 1, 1, 0, 1, 0, 0, 1],
        ] {
            let key = IntegralKey::try_new(powers.to_vec()).unwrap();
            let master = super::super::admission_census_tests::input_for_key(&key).unwrap();
            assert!(
                screened_selection_matches_reference(
                    &engine,
                    &settings,
                    master.integral.as_view(),
                    master.numerator.as_view(),
                    Some(&key),
                ),
                "master {powers:?}"
            );
        }
    }

    /// Compares the screened selection with the brute-force reference (full
    /// keys of every witness, terminal-only witnesses first). With `master`,
    /// the input is that remaining terminal's own integrand and the selected
    /// support must be exactly that key. Returns whether a comparison ran.
    fn screened_selection_matches_reference(
        engine: &Vakint,
        settings: &VakintSettings,
        input: AtomView,
        numerator: AtomView,
        master: Option<&IntegralKey>,
    ) -> bool {
        // The same terminal set as the implementation: remaining terminals of
        // the package's own family (family fingerprint filter).
        let terminals = {
            let application = engine.saved_application(settings).unwrap();
            let prepared = application.prepared.lock().unwrap();
            let family = prepared.family().fingerprint();
            prepared
                .remaining_terminals()
                .iter()
                .filter(|terminal| terminal.family_fingerprint() == family)
                .map(|terminal| terminal.integral().clone())
                .collect::<BTreeSet<IntegralKey>>()
        };
        if let Some(master) = master {
            assert!(
                terminals.contains(master),
                "{master:?} is not a remaining terminal"
            );
        }
        let mut compared = false;
        {
            engine
                .topologies
                .match_topologies_to_user_input_with_selection(
                    input,
                    false,
                    Default::default(),
                    |witnesses| {
                        // Terminal-only witnesses sort first (false < true).
                        let expected =
                            crate::topologies::first_minimal_acceptance(witnesses, |witness| {
                                Ok(engine
                                    .rustred_matched_input_key(settings, numerator, witness)?
                                    .map(|key| {
                                        let outside =
                                            !key.iter().all(|(key, _)| terminals.contains(key));
                                        (outside, key)
                                    }))
                            })?;
                        let (selected, support) =
                            engine.select_rustred_witness(settings, numerator, witnesses)?;
                        assert_eq!(selected, expected, "numerator {numerator}");
                        if let Some(index) = selected {
                            let key = engine
                                .rustred_matched_input_key(settings, numerator, &witnesses[index])?
                                .unwrap();
                            let support = support.unwrap();
                            assert!(key.iter().map(|(key, _)| key).eq(support.iter()));
                            if let Some(master) = master {
                                assert_eq!(support.iter().collect::<Vec<_>>(), vec![master]);
                            }
                        } else {
                            assert!(support.is_none());
                            assert!(master.is_none(), "a master input must be admitted");
                        }
                        compared = true;
                        Ok(selected)
                    },
                )
                .unwrap();
        }
        compared
    }

    #[test]
    fn symbolic_masses_must_be_independent_of_the_coefficients() {
        Vakint::initialize_vakint_symbols();
        let numerator = vk_parse!("x*dot(k(1),k(1)) + dot(k(2),k(2))").unwrap();
        let groups = NumeratorGroups::new(numerator.as_view()).unwrap();
        assert_eq!(
            groups.mass_mode(&vk_parse!("2").unwrap()),
            Some(Some(vk_parse!("2").unwrap()))
        );
        assert_eq!(groups.mass_mode(&vk_parse!("mUV^2").unwrap()), Some(None));
        assert_eq!(groups.mass_mode(&vk_parse!("x*mUV").unwrap()), None);
        assert_eq!(groups.mass_mode(&vk_parse!("0").unwrap()), None);
        assert_eq!(groups.mass_mode(&vk_parse!("2^(1/2)").unwrap()), None);
    }

    /// Recorded selection inputs: the `event: "input"` lines of a
    /// `VAKINT_WITNESS_LOG` (integral and numerator as the selector saw them,
    /// e.g. from the six-loop computations 640, 642, 213 and 216), read from
    /// the file `VAKINT_WS_RECORDED_INPUTS`; only the input ordinals listed in
    /// `VAKINT_WS_RECORDED_LINES` (comma separated) when that is set.
    fn recorded_inputs() -> Vec<(usize, serde_json::Value)> {
        let path = std::env::var_os("VAKINT_WS_RECORDED_INPUTS")
            .expect("VAKINT_WS_RECORDED_INPUTS names a witness-selection record");
        let wanted = std::env::var("VAKINT_WS_RECORDED_LINES").ok().map(|lines| {
            lines
                .split(',')
                .map(|line| line.trim().parse::<usize>().unwrap())
                .collect::<BTreeSet<_>>()
        });
        std::fs::read_to_string(path)
            .unwrap()
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| serde_json::from_str::<serde_json::Value>(line).unwrap())
            .filter(|record| record["event"] == "input")
            .enumerate()
            .filter(|(ordinal, _)| {
                wanted
                    .as_ref()
                    .is_none_or(|wanted| wanted.contains(ordinal))
            })
            .collect()
    }

    /// The official five-loop saved application with the recorded regulator
    /// symbol and `VAKINT_WS_THREADS` frontend threads (default 1).
    fn recorded_engine(epsilon: &str) -> (Vakint, VakintSettings) {
        let directory = super::super::super::test_data::external_path_through(5).join("5L/ibps");
        let threads = std::env::var("VAKINT_WS_THREADS")
            .ok()
            .map_or(1, |threads| threads.parse().unwrap());
        let mut settings = VakintSettings {
            epsilon_symbol: epsilon.into(),
            n_cores: threads,
            ..Default::default()
        };
        let engine = super::super::super::test_data::external_engine()
            .with_rustred_application(&mut settings, &directory)
            .unwrap();
        (engine, settings)
    }

    fn recorded_atom(record: &serde_json::Value, field: &str) -> Atom {
        symbolica::try_parse!(record[field].as_str().unwrap()).unwrap()
    }

    /// Beta's property on real inputs: on every recorded input the screened
    /// selection is the first minimal full key over all verified witnesses
    /// (terminal-only witnesses first), with the same exact support.
    #[test]
    #[ignore = "needs VAKINT_WS_RECORDED_INPUTS and the official five-loop package (minutes per input)"]
    fn screened_selection_matches_first_minimal_full_key_on_recorded_inputs() {
        Vakint::initialize_vakint_symbols();
        let inputs = recorded_inputs();
        assert!(!inputs.is_empty(), "no recorded input selected");
        let epsilon = inputs[0].1["epsilon_symbol"].as_str().unwrap().to_owned();
        let (engine, settings) = recorded_engine(&epsilon);
        for (ordinal, record) in &inputs {
            assert_eq!(record["epsilon_symbol"].as_str(), Some(epsilon.as_str()));
            let integral = recorded_atom(record, "integral");
            let numerator = recorded_atom(record, "numerator");
            let scalar_form =
                VakintTerm::scalar_product_form(&settings, numerator.as_view()).unwrap();
            let grouping = NumeratorGroups::try_new(scalar_form.as_view());
            let started = std::time::Instant::now();
            assert!(
                screened_selection_matches_reference(
                    &engine,
                    &settings,
                    integral.as_view(),
                    numerator.as_view(),
                    None,
                ),
                "input {ordinal}: no witness selection ran"
            );
            println!(
                "{}",
                serde_json::json!({
                    "input": ordinal,
                    "numerator_digest": record["numerator_digest"],
                    "groups": grouping.as_ref().ok().map(NumeratorGroups::shape),
                    "refusal": grouping.as_ref().err().map(ToString::to_string),
                    "equal_to_reference": true,
                    "seconds": started.elapsed().as_secs_f64(),
                })
            );
        }
    }

    /// The production selection call on recorded inputs, timed (one JSON
    /// line per input; with `VAKINT_WITNESS_LOG` also the selection record).
    #[test]
    #[ignore = "timing driver: needs VAKINT_WS_RECORDED_INPUTS and the official five-loop package"]
    fn witness_selection_timing_on_recorded_inputs() {
        Vakint::initialize_vakint_symbols();
        let inputs = recorded_inputs();
        assert!(!inputs.is_empty(), "no recorded input selected");
        let epsilon = inputs[0].1["epsilon_symbol"].as_str().unwrap().to_owned();
        let (engine, settings) = recorded_engine(&epsilon);
        for (ordinal, record) in &inputs {
            let integral = recorded_atom(record, "integral");
            let numerator = recorded_atom(record, "numerator");
            let started = std::time::Instant::now();
            let selected = engine
                .match_rustred_input_selected(&settings, integral.as_view(), numerator.as_view())
                .unwrap();
            let seconds = started.elapsed().as_secs_f64();
            let (witness, support) = match &selected {
                Some((witness, support)) => (
                    Some(super::super::selection_record::witness_text(witness)),
                    support.clone().flatten().map(|support| {
                        support
                            .iter()
                            .map(|key| key.powers().to_vec())
                            .collect::<Vec<_>>()
                    }),
                ),
                None => (None, None),
            };
            println!(
                "{}",
                serde_json::json!({
                    "input": ordinal,
                    "numerator_digest": record["numerator_digest"],
                    "threads": settings.n_cores,
                    "seconds": seconds,
                    "selected_witness_digest": witness.as_deref().map(super::super::selection_record::digest),
                    "selected_witness": witness,
                    "support": support,
                })
            );
        }
    }

    /// Diagnostics of the grouping on recorded inputs (one JSON line each):
    /// the scalar form's top level, the refusal or shape without and with the
    /// opaque atoms, and the same after expansion in the scalar products.
    #[test]
    #[ignore = "diagnostics: needs VAKINT_WS_RECORDED_INPUTS"]
    fn grouping_shape_of_recorded_inputs() {
        Vakint::initialize_vakint_symbols();
        let describe = |grouping: Result<NumeratorGroups, Refusal>| match grouping {
            Ok(groups) => serde_json::json!({
                "shape": groups.shape(),
                "opaque_atoms": groups.opaque_atoms(),
            }),
            Err(refusal) => serde_json::json!({ "refusal": refusal.to_string() }),
        };
        let top = |view: AtomView| match view {
            AtomView::Add(sum) => format!("add({})", sum.iter().count()),
            AtomView::Mul(product) => format!("mul({})", product.iter().count()),
            _ => "other".to_owned(),
        };
        for (ordinal, record) in recorded_inputs() {
            let settings = VakintSettings {
                epsilon_symbol: record["epsilon_symbol"].as_str().unwrap().into(),
                ..Default::default()
            };
            let numerator = recorded_atom(&record, "numerator");
            let started = std::time::Instant::now();
            let scalar_form =
                VakintTerm::scalar_product_form(&settings, numerator.as_view()).unwrap();
            let scalar_seconds = started.elapsed().as_secs_f64();
            let started = std::time::Instant::now();
            let expanded = scalar_form.expand_in(Atom::var(S.dot));
            let expand_seconds = started.elapsed().as_secs_f64();
            let started = std::time::Instant::now();
            let expanded_opaque = describe(NumeratorGroups::try_new_with(expanded.as_view(), true));
            let group_seconds = started.elapsed().as_secs_f64();
            println!(
                "{}",
                serde_json::json!({
                    "input": ordinal,
                    "numerator_digest": record["numerator_digest"],
                    "scalar_form_top": top(scalar_form.as_view()),
                    "scalar_form_chars": scalar_form.to_string().len(),
                    "plain": describe(NumeratorGroups::try_new_with(scalar_form.as_view(), false)),
                    "opaque": describe(NumeratorGroups::try_new_with(scalar_form.as_view(), true)),
                    "expanded_top": top(expanded.as_view()),
                    "expanded_opaque": expanded_opaque,
                    "scalar_seconds": scalar_seconds,
                    "expand_seconds": expand_seconds,
                    "group_seconds": group_seconds,
                })
            );
        }
    }
}
