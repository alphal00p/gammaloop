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
    groups: Vec<Group>,
}

/// Exact weights of one lowered key, by independent slot: distinct
/// coefficient, and group degree when the mass is symbolic.
type KeyWeights = BTreeMap<(usize, u32), Atom>;

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

fn contains_loop_momentum(view: AtomView) -> bool {
    let mut found = false;
    view.visitor(&mut |part| {
        if let AtomView::Fun(function) = part
            && (function.get_symbol() == S.dot || function.get_symbol() == S.k)
        {
            found = true;
        }
        !found
    });
    found
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

/// Evaluate a plain coefficient at a rational point, or `None` at a pole.
fn evaluate(coefficient: &Atom, point: &[Replacement]) -> Option<Atom> {
    let value = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        coefficient.replace_multiple(point)
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
    pub(super) fn new(numerator: AtomView) -> Option<Self> {
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
                if contains_loop_momentum(factor) {
                    degree = degree.checked_add(scalar_product_degree(factor)?)?;
                    scalar_products *= factor.to_owned();
                } else if let AtomView::Num(_) = factor {
                    rational(factor)?;
                    number *= factor.to_owned();
                } else if is_plain(factor) {
                    coefficient *= factor.to_owned();
                } else {
                    return None;
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
            groups,
        };
        decomposition
            .coefficients_are_independent()
            .then_some(decomposition)
    }

    /// Distinct coefficients and groups, for diagnostics.
    pub(super) fn shape(&self) -> (usize, usize) {
        (self.coefficients.len(), self.groups.len())
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
        let mut rows = Vec::with_capacity(count);
        for coefficient in &self.coefficients {
            let mut row = Vec::with_capacity(points.len());
            for point in &points {
                let Some(value) = evaluate(coefficient, point) else {
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
    /// witnesses, with the official five-loop package's real admission.
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
            let mut compared = false;
            engine
                .topologies
                .match_topologies_to_user_input_with_selection(
                    input.as_view(),
                    false,
                    Default::default(),
                    |witnesses| {
                        let expected =
                            crate::topologies::first_minimal_acceptance(witnesses, |witness| {
                                engine.rustred_matched_input_key(
                                    &settings,
                                    numerator.as_view(),
                                    witness,
                                )
                            })?;
                        let (selected, support) = engine.select_rustred_witness(
                            &settings,
                            numerator.as_view(),
                            witnesses,
                        )?;
                        assert_eq!(selected, expected, "numerator {numerator}");
                        if let Some(index) = selected {
                            let key = engine
                                .rustred_matched_input_key(
                                    &settings,
                                    numerator.as_view(),
                                    &witnesses[index],
                                )?
                                .unwrap();
                            assert!(key.iter().map(|(key, _)| key).eq(support.unwrap().iter()));
                        } else {
                            assert!(support.is_none());
                        }
                        compared = true;
                        Ok(selected)
                    },
                )
                .unwrap();
            assert!(compared);
        }
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
}
