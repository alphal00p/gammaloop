//! Exact energy identities in a final, test-only scalar projection view.
//!
//! Callers must retain graph numerators outside this algebra. Only scalar
//! denominator coefficients may enter the bounded polynomial conversion.

use std::{
    collections::{BTreeMap, BTreeSet},
    sync::Arc,
};

use color_eyre::eyre::{Result, ensure, eyre};
use spenso::network::parsing::{AtomStructureExt, StrictTensorFilter};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::ConvertToRing,
    domains::{Ring, algebraic::AlgebraicExtension, rational::Q},
    id::Replacement,
    symbol,
};

use super::{
    pole_certificate::{Polynomial, ScalarPoleCertificate, ScalarPoleFailure},
    symbolica_ext::NumeratorAtomExt,
};
use crate::{cff::expression::OrientationID, utils::GS};

const MAX_INVARIANT_BYTES: usize = 8 * 1024;
const MAX_SCALAR_BYTES: usize = 1024 * 1024;
const MAX_TERMS: usize = 4096;
const MAX_LEAVES: usize = 256;
const MAX_ENERGIES: usize = 128;

#[derive(Debug)]
pub(crate) struct EnergyZeroWitness {
    pub original_numerator: Atom,
    pub numerator_remainder: Atom,
    pub denominator_remainder: Atom,
    /// Each entry witnesses N_previous = quotient * (E² - P) + N_next.
    pub relations: Vec<(Atom, Atom)>,
}

#[derive(Debug)]
pub(crate) enum EnergyZeroOutcome {
    Zero(EnergyZeroWitness),
    /// A nonzero remainder means this sufficient certificate did not prove
    /// zero. Further physical relations can still make the expression vanish.
    Nonzero(EnergyZeroWitness),
    Unproven(String),
}

#[derive(Default)]
pub(crate) struct SoftEnergyAlgebra {
    invariants: BTreeMap<Atom, Atom>,
    parameters: BTreeSet<Symbol>,
}

impl SoftEnergyAlgebra {
    /// This is only valid after the final forest operation: owner provenance
    /// is discarded and energies with the same physical invariant are merged.
    /// The soft parameter approaches zero from positive values.
    pub(crate) fn normalize_energies(&mut self, atom: &Atom, lambda: Symbol) -> Result<Atom> {
        self.parameters.insert(lambda);
        let mut calls = Vec::new();
        atom.visitor(&mut |view| {
            if let AtomView::Fun(call) = view
                && call.get_symbol() == GS.energy_surface
            {
                if !calls.iter().any(|old: &Atom| old.as_view() == view) {
                    calls.push(view.to_owned());
                }
                return false;
            }
            true
        });
        let parameter = Atom::var(lambda);
        let mut replacements = Vec::with_capacity(calls.len());
        for call in calls {
            let AtomView::Fun(energy) = call.as_view() else {
                unreachable!()
            };
            ensure!(
                energy.get_nargs() == 2,
                ScalarPoleFailure::Unsupported("energy arity")
            );
            let invariant = Self::invariant(energy.get(1))?;
            let terms = invariant
                .coefficient_list_exact(std::slice::from_ref(&parameter))
                .ok_or_else(|| eyre!(ScalarPoleFailure::Budget("energy invariant coefficients")))?;
            let mut powers = Vec::new();
            for (monomial, coefficient) in terms {
                if coefficient.is_zero() {
                    continue;
                }
                let power = if monomial.is_one() {
                    0
                } else if monomial == parameter {
                    1
                } else if let AtomView::Pow(power) = monomial.as_view()
                    && power.get_base() == parameter.as_view()
                {
                    i64::try_from(power.get_exp()).map_err(|_| {
                        eyre!(ScalarPoleFailure::Unsupported("noninteger invariant power"))
                    })?
                } else {
                    return Err(eyre!(ScalarPoleFailure::Unsupported(
                        "energy Laurent monomial"
                    )));
                };
                ensure!(
                    power >= 0 && !coefficient.contains_symbol(lambda),
                    ScalarPoleFailure::Unsupported("non-polynomial energy invariant")
                );
                powers.push((power, coefficient));
            }
            let valuation = powers
                .iter()
                .map(|(power, _)| *power)
                .min()
                .ok_or_else(|| eyre!(ScalarPoleFailure::Unsupported("identically zero energy")))?;
            ensure!(
                valuation % 2 == 0,
                ScalarPoleFailure::Unsupported("half-integer soft energy order")
            );
            let reduced = powers.iter().fold(Atom::Zero, |sum, (power, coefficient)| {
                sum + coefficient * parameter.pow(*power - valuation)
            });
            let base = Self::invariant(reduced.replace(lambda).with(0).as_view())?;
            let base_energy = GS.energy_surface.call_args([Atom::Zero, base.clone()]);
            ensure!(
                self.invariants.contains_key(&base_energy) || self.invariants.len() < MAX_ENERGIES,
                ScalarPoleFailure::Budget("physical energies")
            );
            self.invariants.insert(base_energy, base);
            let normalized =
                parameter.pow(valuation / 2) * GS.energy_surface.call_args([Atom::Zero, reduced]);
            replacements.push(Replacement::new(call.to_pattern(), normalized));
        }
        Ok(atom.replace_multiple(&replacements))
    }

    /// A zero result is exact ideal membership over Q(i), valid away from the
    /// original denominators. Nonzero is an unresolved remainder, not a claim
    /// that all possible physical energy identities have been exhausted.
    pub(crate) fn scalar_zero(&self, scalar: &Atom) -> Result<EnergyZeroOutcome> {
        match self.certificate(scalar) {
            Ok(witness) if witness.numerator_remainder.is_zero() => {
                Ok(EnergyZeroOutcome::Zero(witness))
            }
            Ok(witness) => Ok(EnergyZeroOutcome::Nonzero(witness)),
            Err(error)
                if matches!(
                    error.downcast_ref::<ScalarPoleFailure>(),
                    Some(ScalarPoleFailure::Unsupported(_) | ScalarPoleFailure::Budget(_))
                ) =>
            {
                Ok(EnergyZeroOutcome::Unproven(error.to_string()))
            }
            Err(error) => Err(error),
        }
    }

    /// Sufficient generic nonzero certificate after the caller has put every
    /// physical scalar leaf into one independent coordinate/parameter chart.
    /// Multiplying by all energy conjugates gives a rational polynomial norm;
    /// a nonzero norm rules out zero on every square-root branch. A zero norm
    /// or exceeded bound remains inconclusive, including dependent radicals.
    pub(crate) fn scalar_nonzero(&self, scalar: &Atom) -> Result<bool> {
        let witness = match self.certificate(scalar) {
            Ok(witness) => witness,
            Err(error)
                if matches!(
                    error.downcast_ref::<ScalarPoleFailure>(),
                    Some(ScalarPoleFailure::Unsupported(_) | ScalarPoleFailure::Budget(_))
                ) =>
            {
                return Ok(false);
            }
            Err(error) => return Err(error),
        };
        for mut polynomial in [witness.numerator_remainder, witness.denominator_remainder] {
            let mut energies = Vec::new();
            polynomial.visitor(&mut |part| {
                if matches!(part, AtomView::Fun(call) if call.get_symbol() == GS.energy_surface) {
                    if !energies.iter().any(|old: &Atom| old.as_view() == part) {
                        energies.push(part.to_owned());
                    }
                    return false;
                }
                true
            });
            for energy in energies {
                let conjugate = polynomial.replace(energy.to_pattern()).with(-&energy);
                let product = &polynomial * conjugate;
                polynomial = match self.certificate(&product) {
                    Ok(witness) => {
                        ensure!(
                            witness.denominator_remainder.is_one(),
                            ScalarPoleFailure::Failed("energy norm acquired a denominator")
                        );
                        witness.numerator_remainder
                    }
                    Err(error)
                        if matches!(
                            error.downcast_ref::<ScalarPoleFailure>(),
                            Some(ScalarPoleFailure::Unsupported(_) | ScalarPoleFailure::Budget(_))
                        ) =>
                    {
                        return Ok(false);
                    }
                    Err(error) => return Err(error),
                };
                ensure!(
                    !polynomial.contains(&energy),
                    ScalarPoleFailure::Failed("energy conjugation did not eliminate its root")
                );
            }
            if polynomial.is_zero() {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn certificate(&self, scalar: &Atom) -> Result<EnergyZeroWitness> {
        ensure!(
            scalar.as_view().get_byte_size() <= MAX_SCALAR_BYTES,
            ScalarPoleFailure::Budget("scalar bytes")
        );
        let field = AlgebraicExtension::complex(Q);
        let mut scalar_leaves = Vec::new();
        Self::leaves(scalar.as_view(), &field, &mut scalar_leaves)?;
        let energies = scalar_leaves.iter().filter_map(|leaf| {
            if matches!(leaf.as_view(), AtomView::Fun(call) if call.get_symbol() == GS.energy_surface) {
                Some(leaf.clone())
            } else {
                None
            }
        }).collect::<Vec<_>>();
        let mut leaves = energies.clone();
        let mut invariants = Vec::new();
        for energy in &energies {
            let invariant = if let Some(invariant) = self.invariants.get(energy) {
                invariant.clone()
            } else {
                let AtomView::Fun(call) = energy.as_view() else {
                    unreachable!()
                };
                ensure!(
                    call.get_nargs() == 2 && call.get(0).is_zero(),
                    ScalarPoleFailure::Unsupported("energy owner not finalized")
                );
                Self::invariant(call.get(1))?
            };
            ensure!(
                self.parameters
                    .iter()
                    .all(|parameter| !invariant.contains_symbol(*parameter)),
                ScalarPoleFailure::Unsupported("Taylor-dependent energy in scalar coefficient")
            );
            Self::leaves(invariant.as_view(), &field, &mut scalar_leaves)?;
            invariants.push(invariant);
        }
        for leaf in scalar_leaves {
            if !leaves.contains(&leaf) {
                leaves.push(leaf);
            }
        }
        ensure!(
            leaves.len() <= MAX_LEAVES,
            ScalarPoleFailure::Budget("scalar leaves")
        );
        let variables = leaves
            .iter()
            .cloned()
            .map(TryInto::try_into)
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(|_| eyre!(ScalarPoleFailure::Unsupported("non-indeterminate leaf")))?;
        let template = Polynomial::new(&field, None, Arc::new(variables));
        let (mut numerator, mut denominator) =
            ScalarPoleCertificate::fraction(scalar.as_view(), &leaves, &template)?;
        let original_numerator = numerator.to_expression();
        let mut relations = Vec::new();
        // With energies ordered first, e_i² is each relation's leading
        // monomial. Invariants contain no energies, so these monomials are
        // pairwise coprime and reducing one never reintroduces another.
        for (index, invariant) in invariants.iter().enumerate() {
            let (radicand, radicand_denominator) =
                ScalarPoleCertificate::fraction(invariant.as_view(), &leaves, &template)?;
            ensure!(
                radicand_denominator.is_one(),
                ScalarPoleFailure::Unsupported("rational energy invariant")
            );
            let mut exponents = vec![0; leaves.len()];
            exponents[index] = 2;
            let relation = template.monomial(field.one(), exponents) - radicand.clone();
            for polynomial in [&numerator, &denominator] {
                let bound = polynomial.exponents_iter().fold(0usize, |sum, powers| {
                    sum.saturating_add(radicand.nterms().saturating_pow(powers[index] / 2))
                });
                ensure!(
                    bound <= MAX_TERMS,
                    ScalarPoleFailure::Budget("energy reduction")
                );
            }
            let (quotient, remainder) = numerator.quot_rem(&relation, false);
            ensure!(
                ScalarPoleCertificate::sum(
                    &ScalarPoleCertificate::multiply(&quotient, &relation)?,
                    &remainder
                )? == numerator,
                ScalarPoleFailure::Failed("energy division reconstruction")
            );
            numerator = remainder;
            denominator = denominator.rem(&relation);
            if !quotient.is_zero() {
                relations.push((relation.to_expression(), quotient.to_expression()));
            }
        }
        ensure!(
            !denominator.is_zero(),
            ScalarPoleFailure::Unsupported("denominator vanishes under energy identities")
        );
        Ok(EnergyZeroWitness {
            original_numerator,
            numerator_remainder: numerator.to_expression(),
            denominator_remainder: denominator.to_expression(),
            relations,
        })
    }

    fn invariant(atom: AtomView<'_>) -> Result<Atom> {
        ensure!(
            atom.get_byte_size() <= MAX_INVARIANT_BYTES,
            ScalarPoleFailure::Budget("energy invariant bytes")
        );
        let field = AlgebraicExtension::complex(Q);
        Self::polynomial_size(atom, &field)?;
        // Only the squared scalar energy invariant is expanded; surrounding
        // numerator factors and their tensor/vertex definitions stay untouched.
        Ok(atom.expand())
    }

    fn polynomial_size(atom: AtomView<'_>, field: &AlgebraicExtension<Q>) -> Result<(usize, u32)> {
        let size = match atom {
            AtomView::Num(number) => {
                field
                    .try_element_from_coefficient_view(number.get_coeff_view())
                    .map_err(|_| {
                        eyre!(ScalarPoleFailure::Unsupported(
                            "inexact invariant coefficient"
                        ))
                    })?;
                (1, 0)
            }
            AtomView::Var(_) | AtomView::Fun(_) => {
                let mut leaves = Vec::new();
                Self::leaves(atom, field, &mut leaves)?;
                ensure!(
                    !atom.contains_symbol(GS.energy_surface),
                    ScalarPoleFailure::Unsupported("nested energy invariant")
                );
                (1, 1)
            }
            AtomView::Pow(power) => {
                let exponent = u32::try_from(i64::try_from(power.get_exp()).map_err(|_| {
                    eyre!(ScalarPoleFailure::Unsupported("noninteger invariant power"))
                })?)
                .map_err(|_| eyre!(ScalarPoleFailure::Unsupported("negative invariant power")))?;
                ensure!(exponent <= 16, ScalarPoleFailure::Budget("invariant power"));
                let (terms, degree) = Self::polynomial_size(power.get_base(), field)?;
                (
                    terms.saturating_pow(exponent),
                    degree.saturating_mul(exponent),
                )
            }
            AtomView::Add(sum) => {
                let mut result = (0usize, 0u32);
                for term in sum {
                    let size = Self::polynomial_size(term, field)?;
                    result = (result.0.saturating_add(size.0), result.1.max(size.1));
                }
                result
            }
            AtomView::Mul(product) => {
                let mut result = (1usize, 0u32);
                for term in product {
                    let size = Self::polynomial_size(term, field)?;
                    result = (
                        result.0.saturating_mul(size.0),
                        result.1.saturating_add(size.1),
                    );
                }
                result
            }
        };
        ensure!(
            size.0 <= MAX_TERMS && size.1 <= 16,
            ScalarPoleFailure::Budget("invariant polynomial")
        );
        Ok(size)
    }

    fn leaves(
        atom: AtomView<'_>,
        field: &AlgebraicExtension<Q>,
        leaves: &mut Vec<Atom>,
    ) -> Result<()> {
        match atom {
            AtomView::Num(number) => {
                field
                    .try_element_from_coefficient_view(number.get_coeff_view())
                    .map_err(|_| {
                        eyre!(ScalarPoleFailure::Unsupported(
                            "coefficient outside exact Q(i)"
                        ))
                    })?;
            }
            AtomView::Var(_) | AtomView::Fun(_) => {
                let explicit_component = matches!(atom, AtomView::Fun(call)
                    if [GS.emr_vec, GS.emr_mom].contains(&call.get_symbol())
                    && call.get_nargs() == 2
                    && (0..=3).any(|index| call.get(1) == GS.cind(index).as_view()));
                let energy =
                    matches!(atom, AtomView::Fun(call) if call.get_symbol() == GS.energy_surface);
                ensure!(
                    energy
                        || explicit_component
                        || !atom.is_tensorial(StrictTensorFilter::ContainsReps),
                    ScalarPoleFailure::Unsupported("tensor in scalar energy coefficient")
                );
                for head in [
                    Symbol::IF,
                    OrientationID::symbol(),
                    GS.theta,
                    GS.orientation_delta,
                    symbol!("gammalooprs::uv::numerator_family"),
                ] {
                    ensure!(
                        !atom.contains_symbol(head),
                        ScalarPoleFailure::Unsupported("guard or retained numerator")
                    );
                }
                if !leaves.iter().any(|leaf| leaf.as_view() == atom) {
                    ensure!(
                        leaves.len() < MAX_LEAVES,
                        ScalarPoleFailure::Budget("scalar leaves")
                    );
                    leaves.push(atom.to_owned());
                }
            }
            AtomView::Pow(power) => Self::leaves(power.get_base(), field, leaves)?,
            AtomView::Add(sum) => {
                for term in sum {
                    Self::leaves(term, field, leaves)?;
                }
            }
            AtomView::Mul(product) => {
                for factor in product {
                    Self::leaves(factor, field, leaves)?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::{function, parse_lit};

    #[test]
    fn energy_zero_proves_quadratic_relation_and_keeps_nonzero_remainder() -> Result<()> {
        crate::initialisation::test_initialise()?;
        let lambda = symbol!("soft_energy_certificate::lambda");
        let p = parse_lit!(x ^ 2 + y ^ 2 + m ^ 2);
        let mut algebra = SoftEnergyAlgebra::default();
        let energy = algebra.normalize_energies(&function!(GS.energy_surface, 7, &p), lambda)?;
        let EnergyZeroOutcome::Zero(witness) =
            algebra.scalar_zero(&((&energy * &energy - &p) / (&energy + 1)))?
        else {
            panic!("the exact energy relation must vanish")
        };
        assert!(!witness.relations.is_empty());
        assert!(matches!(
            algebra.scalar_zero(&(energy.pow(2) - p + 1))?,
            EnergyZeroOutcome::Nonzero(_)
        ));
        Ok(())
    }

    #[test]
    fn energy_normalization_extracts_positive_soft_scale_and_merges_final_owners() -> Result<()> {
        crate::initialisation::test_initialise()?;
        let lambda = symbol!("soft_energy_certificate::lambda");
        let t = Atom::var(lambda);
        let p = parse_lit!(x ^ 2 + y ^ 2);
        let mut algebra = SoftEnergyAlgebra::default();
        let left =
            algebra.normalize_energies(&function!(GS.energy_surface, 7, t.pow(2) * &p), lambda)?;
        let right = algebra.normalize_energies(&function!(GS.energy_surface, 9, &p), lambda)?;
        assert_eq!(left, &t * right);
        let hard = algebra.normalize_energies(
            &function!(
                GS.energy_surface,
                3,
                (parse_lit!(x) + &t * parse_lit!(y)).pow(2) + 1
            ),
            lambda,
        )?;
        let coefficients =
            super::super::exact_soft_jet::ExactSoftJet::new(lambda, 2).coefficients(&hard)?;
        assert_eq!(coefficients.keys().copied().collect::<Vec<_>>(), [0, 1, 2]);
        let base = function!(GS.energy_surface, 0, parse_lit!(x ^ 2 + 1));
        let expected = [
            base.clone(),
            parse_lit!(x * y) / &base,
            parse_lit!(y ^ 2) / (base.pow(3) * 2),
        ];
        for (degree, expected) in expected.iter().enumerate() {
            let difference = &coefficients[&(degree as i64)] - expected;
            assert!(
                matches!(
                    algebra.scalar_zero(&difference)?,
                    EnergyZeroOutcome::Zero(_)
                ),
                "incorrect finite energy coefficient at degree {degree}: {difference}"
            );
        }
        // Independently convolve the three exact coefficients. This retains
        // the complete second-order identity, including its nonzero y² term.
        let invariant = [
            parse_lit!(x ^ 2 + 1),
            parse_lit!(2 * x * y),
            parse_lit!(y ^ 2),
        ];
        for (degree, expected) in invariant.iter().enumerate() {
            let squared = Atom::add_many((0..=degree).map(|left| {
                &coefficients[&(left as i64)] * &coefficients[&((degree - left) as i64)]
            }));
            let difference = squared - expected;
            assert!(
                matches!(
                    algebra.scalar_zero(&difference)?,
                    EnergyZeroOutcome::Zero(_)
                ),
                "energy-square identity fails at degree {degree}: {difference}"
            );
        }
        Ok(())
    }

    #[test]
    fn energy_certificate_accepts_explicit_components_and_autoregisters_final_energies()
    -> Result<()> {
        crate::initialisation::test_initialise()?;
        let q = GS.emr_vec(linnet::half_edge::involution::EdgeIndex(6), GS.cind(1));
        let invariant = q.pow(2) + 1;
        let energy = function!(GS.energy_surface, 0, &invariant);
        let algebra = SoftEnergyAlgebra::default();
        assert!(matches!(
            algebra.scalar_zero(&(energy.pow(2) - &invariant))?,
            EnergyZeroOutcome::Zero(_)
        ));
        assert!(matches!(
            algebra.scalar_zero(&(Atom::one() / (energy.pow(2) - invariant)))?,
            EnergyZeroOutcome::Unproven(_)
        ));
        Ok(())
    }

    #[test]
    fn energy_norm_certifies_generic_nonzero_but_not_a_dependent_root() -> Result<()> {
        crate::initialisation::test_initialise()?;
        let algebra = SoftEnergyAlgebra::default();
        let first = function!(GS.energy_surface, 0, parse_lit!(x ^ 2 + y ^ 2));
        let second = function!(GS.energy_surface, 0, parse_lit!(x ^ 2 + z ^ 2));
        assert!(algebra.scalar_nonzero(&(&first - &second))?);
        assert!(algebra.scalar_nonzero(&(Atom::one() / (&first + &second)))?);
        let dependent = function!(GS.energy_surface, 0, parse_lit!(x ^ 2));
        assert!(!algebra.scalar_nonzero(&(dependent - parse_lit!(x)))?);
        Ok(())
    }

    #[test]
    fn energy_certificate_reports_unsupported_and_budget_limits() -> Result<()> {
        crate::initialisation::test_initialise()?;
        let lambda = symbol!("soft_energy_certificate::lambda");
        let mut algebra = SoftEnergyAlgebra::default();
        assert!(
            algebra
                .normalize_energies(&function!(GS.energy_surface, 0, Atom::var(lambda)), lambda)
                .is_err()
        );
        assert!(
            algebra
                .normalize_energies(
                    &function!(GS.energy_surface, 0, parse_lit!((x + y) ^ 100)),
                    lambda
                )
                .is_err()
        );
        for scalar in [
            parse_lit!(x ^ (1 / 2)),
            Atom::num(0.125f64),
            function!(symbol!("gammalooprs::uv::numerator_family"), 0),
            function!(GS.energy_surface, 9, parse_lit!(x)),
        ] {
            assert!(matches!(
                algebra.scalar_zero(&scalar)?,
                EnergyZeroOutcome::Unproven(_)
            ));
        }
        Ok(())
    }
}
