//! Offline certificates for scalar source pairs, never part of generation.
//!
//! Keep unreduced polynomial pairs: a normal rational conversion would cancel
//! the very denominator factors whose disappearance we want to establish.

use std::sync::Arc;

use color_eyre::eyre::{Result, bail, ensure, eyre};
use spenso::network::parsing::{AtomStructureExt, StrictTensorFilter};
use symbolica::{
    atom::{Atom, AtomCore, AtomView, Symbol},
    coefficient::ConvertToRing,
    domains::{Ring, algebraic::AlgebraicExtension, rational::Q},
    poly::polynomial::MultivariatePolynomial,
    symbol,
};

use crate::{cff::expression::OrientationID, utils::GS};

pub(super) type Polynomial = MultivariatePolynomial<AlgebraicExtension<Q>, u32>;
type Fraction = (Polynomial, Polynomial);

const MAX_BYTES: usize = 16 * 1024;
const MAX_LEAVES: usize = 128;
const MAX_TERMS: usize = 4096;
const MAX_DEGREE: u32 = 64;

/// Callers may skip unsupported/budget-limited pairs, but a failed proof must
/// fail the audit. The variants remain available through eyre's downcast_ref.
#[derive(Debug, thiserror::Error)]
pub(crate) enum ScalarPoleFailure {
    #[error("unsupported scalar certificate: {0}")]
    Unsupported(&'static str),
    #[error("unproven scalar certificate: budget exceeded ({0})")]
    Budget(&'static str),
    #[error("failed scalar certificate: {0}")]
    Failed(&'static str),
}

#[derive(Debug)]
pub(crate) struct ScalarPoleCertificate {
    pub identity_verified: bool,
    /// An identically zero coefficient has no distinguished removed pole.
    pub scalar_zero: bool,
    /// The monic product of lost denominator factors, with multiplicities.
    /// It is one for no loss and for the separately classified zero case.
    pub removed_factor: Atom,
    pub removed_factor_degree: u32,
    pub numerator_divisibility_verified: bool,
    pub before_num_terms: usize,
    pub before_den_terms: usize,
    pub after_num_terms: usize,
    pub after_den_terms: usize,
    pub before_den_degree: u32,
    pub after_den_degree: u32,
    pub cross_product_terms: [usize; 2],
    /// Complete functions retain their original owner and argument atoms.
    pub opaque_leaves: Vec<Atom>,
}

impl ScalarPoleCertificate {
    pub(crate) fn check(before: &Atom, after: &Atom) -> Result<Self> {
        let field = AlgebraicExtension::complex(Q);
        let mut leaves = Vec::new();
        for expression in [before, after] {
            ensure!(
                expression.as_view().get_byte_size() <= MAX_BYTES,
                ScalarPoleFailure::Budget("bytes")
            );
            for head in [
                Symbol::IF,
                OrientationID::symbol(),
                GS.theta,
                GS.orientation_delta,
                symbol!("gammalooprs::uv::numerator_family"),
            ] {
                ensure!(
                    !expression.contains_symbol(head),
                    ScalarPoleFailure::Unsupported("guard or retained numerator")
                );
            }
            Self::collect_leaves(expression.as_view(), &field, &mut leaves)?;
        }
        let variables = leaves
            .iter()
            .cloned()
            .map(TryInto::try_into)
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(|_| eyre!(ScalarPoleFailure::Unsupported("non-indeterminate leaf")))?;
        let template = Polynomial::new(&field, None, Arc::new(variables));
        let (nb, db) = Self::fraction(before.as_view(), &leaves, &template)?;
        let (na, da) = Self::fraction(after.as_view(), &leaves, &template)?;
        let left = Self::multiply(&nb, &da)?;
        let right = Self::multiply(&na, &db)?;
        let cross_product_terms = [left.nterms(), right.nterms()];
        ensure!(
            Self::sum(&left, &(-right))?.is_zero(),
            ScalarPoleFailure::Failed("nonzero exact identity cross-product")
        );

        let scalar_zero = nb.is_zero();
        let removed = if scalar_zero {
            // Do not assign arbitrary denominator factors to the identity 0=0.
            ensure!(na.is_zero(), ScalarPoleFailure::Failed("zero mismatch"));
            template.one()
        } else {
            let common = db.gcd(&da);
            let removed = Self::divide(&db, &common)?.make_monic();
            let quotient = Self::divide(&nb, &removed)?;
            ensure!(
                Self::multiply(&removed, &quotient)? == nb,
                ScalarPoleFailure::Failed("numerator divisibility witness")
            );
            removed
        };
        Ok(Self {
            identity_verified: true,
            scalar_zero,
            removed_factor: removed.to_expression(),
            removed_factor_degree: Self::degree(&removed),
            numerator_divisibility_verified: !scalar_zero,
            before_num_terms: nb.nterms(),
            before_den_terms: db.nterms(),
            after_num_terms: na.nterms(),
            after_den_terms: da.nterms(),
            before_den_degree: Self::degree(&db),
            after_den_degree: Self::degree(&da),
            cross_product_terms,
            opaque_leaves: leaves,
        })
    }

    fn collect_leaves(
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
                let owned_energy =
                    matches!(atom, AtomView::Fun(call) if call.get_symbol() == GS.energy_surface);
                ensure!(
                    owned_energy || !atom.is_tensorial(StrictTensorFilter::ContainsReps),
                    ScalarPoleFailure::Unsupported("tensor in scalar coefficient")
                );
                if !leaves.iter().any(|leaf| leaf.as_view() == atom) {
                    ensure!(
                        leaves.len() < MAX_LEAVES,
                        ScalarPoleFailure::Budget("opaque leaves")
                    );
                    leaves.push(atom.to_owned());
                }
                // In particular, never inspect or rewrite an energy's arguments.
            }
            AtomView::Pow(power) => {
                Self::exponent(power.get_exp())?;
                Self::collect_leaves(power.get_base(), field, leaves)?;
            }
            AtomView::Add(sum) => {
                for term in sum {
                    Self::collect_leaves(term, field, leaves)?;
                }
            }
            AtomView::Mul(product) => {
                for factor in product {
                    Self::collect_leaves(factor, field, leaves)?;
                }
            }
        }
        Ok(())
    }

    fn exponent(atom: AtomView<'_>) -> Result<i64> {
        let exponent = i64::try_from(atom)
            .map_err(|_| eyre!(ScalarPoleFailure::Unsupported("noninteger power")))?;
        ensure!(
            exponent.unsigned_abs() <= u64::from(MAX_DEGREE),
            ScalarPoleFailure::Budget("exponent")
        );
        Ok(exponent)
    }

    pub(super) fn fraction(
        atom: AtomView<'_>,
        leaves: &[Atom],
        template: &Polynomial,
    ) -> Result<Fraction> {
        match atom {
            AtomView::Num(number) => Ok((
                template.constant(
                    template
                        .ring()
                        .try_element_from_coefficient_view(number.get_coeff_view())
                        .map_err(|_| {
                            eyre!(ScalarPoleFailure::Unsupported(
                                "coefficient outside exact Q(i)"
                            ))
                        })?,
                ),
                template.one(),
            )),
            AtomView::Var(_) | AtomView::Fun(_) => {
                let index = leaves
                    .iter()
                    .position(|leaf| leaf.as_view() == atom)
                    .ok_or_else(|| eyre!(ScalarPoleFailure::Failed("missing opaque leaf")))?;
                let mut exponents = vec![0; leaves.len()];
                exponents[index] = 1;
                Ok((
                    template.monomial(template.ring().one(), exponents),
                    template.one(),
                ))
            }
            AtomView::Pow(power) => {
                let exponent = Self::exponent(power.get_exp())?;
                let (mut num, mut den) = Self::fraction(power.get_base(), leaves, template)?;
                if exponent < 0 {
                    std::mem::swap(&mut num, &mut den);
                }
                ensure!(
                    !den.is_zero(),
                    ScalarPoleFailure::Failed("identically zero denominator")
                );
                Ok((
                    Self::power(num, exponent.unsigned_abs() as u32)?,
                    Self::power(den, exponent.unsigned_abs() as u32)?,
                ))
            }
            AtomView::Mul(product) => {
                let mut result = (template.one(), template.one());
                for factor in product {
                    let (num, den) = Self::fraction(factor, leaves, template)?;
                    result = (
                        Self::multiply(&result.0, &num)?,
                        Self::multiply(&result.1, &den)?,
                    );
                }
                Ok(result)
            }
            AtomView::Add(sum) => {
                let mut result = (template.zero(), template.one());
                for term in sum {
                    let (num, den) = Self::fraction(term, leaves, template)?;
                    // LCM avoids manufacturing a removable x from 1/x+1/x.
                    let gcd = result.1.gcd(&den);
                    let left_scale = Self::divide(&den, &gcd)?;
                    let right_scale = Self::divide(&result.1, &gcd)?;
                    result = (
                        Self::sum(
                            &Self::multiply(&result.0, &left_scale)?,
                            &Self::multiply(&num, &right_scale)?,
                        )?,
                        Self::multiply(&result.1, &left_scale)?,
                    );
                }
                Ok(result)
            }
        }
    }

    fn degree(poly: &Polynomial) -> u32 {
        if poly.nvars() == 0 {
            return 0;
        }
        poly.exponents_iter()
            .map(|powers| powers.iter().copied().sum())
            .max()
            .unwrap_or(0)
    }

    pub(super) fn multiply(left: &Polynomial, right: &Polynomial) -> Result<Polynomial> {
        ensure!(
            left.nterms().saturating_mul(right.nterms()) <= MAX_TERMS
                && Self::degree(left).saturating_add(Self::degree(right)) <= MAX_DEGREE,
            ScalarPoleFailure::Budget("polynomial product")
        );
        Ok(left * right)
    }

    pub(super) fn sum(left: &Polynomial, right: &Polynomial) -> Result<Polynomial> {
        ensure!(
            left.nterms().saturating_add(right.nterms()) <= MAX_TERMS,
            ScalarPoleFailure::Budget("polynomial sum")
        );
        Ok(left + right)
    }

    fn power(poly: Polynomial, exponent: u32) -> Result<Polynomial> {
        ensure!(
            poly.nterms().saturating_pow(exponent) <= MAX_TERMS
                && Self::degree(&poly).saturating_mul(exponent) <= MAX_DEGREE,
            ScalarPoleFailure::Budget("polynomial power")
        );
        Ok(poly.pow(exponent as usize))
    }

    fn divide(num: &Polynomial, den: &Polynomial) -> Result<Polynomial> {
        let Some(quotient) = num.try_div(den) else {
            bail!(ScalarPoleFailure::Failed("exact polynomial division"));
        };
        ensure!(
            Self::multiply(&quotient, den)? == *num,
            ScalarPoleFailure::Failed("exact division reconstruction")
        );
        Ok(quotient)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::{coefficient::Coefficient, function, parse_lit};

    #[test]
    fn certifies_real_removable_pole_and_complex_coefficients() {
        let before = parse_lit!((1 / m - 1 / (m + s)) / s);
        let after = parse_lit!(1 / (m * (m + s)));
        for phase in [Atom::num(1), parse_lit!(1i)] {
            let proof =
                ScalarPoleCertificate::check(&(&phase * &before), &(&phase * &after)).unwrap();
            assert!(proof.identity_verified && proof.numerator_divisibility_verified);
            assert!(!proof.scalar_zero);
            assert_eq!(proof.removed_factor, parse_lit!(s));
            assert_eq!(proof.removed_factor_degree, 1);
            assert_eq!(proof.before_den_degree, proof.after_den_degree + 1);
        }
    }

    #[test]
    fn common_denominators_and_formatting_do_not_fake_pole_removal() {
        for (before, after) in [
            (parse_lit!(1 / x + 1 / x), parse_lit!(2 / x)),
            (parse_lit!(1 / x + 1 / (2 * x)), parse_lit!(3 / (2 * x))),
            (parse_lit!(1 / x + 1 / y), parse_lit!((x + y) / (x * y))),
            (
                parse_lit!(1 / (x * (x + 1)) + 1 / (x * (x + 2))),
                parse_lit!((2 * x + 3) / (x * (x + 1) * (x + 2))),
            ),
        ] {
            let proof = ScalarPoleCertificate::check(&before, &after).unwrap();
            assert_eq!(proof.removed_factor, Atom::num(1));
            assert_eq!(proof.before_den_degree, proof.after_den_degree);
        }
    }

    #[test]
    fn tiny_nonzero_coefficients_and_unequal_pairs_remain_distinct() {
        let before = parse_lit!(1 / (x + 1) - 1 / (x + 1 + 1 / 10 ^ 30));
        let after = parse_lit!(1 / (10 ^ 30 * (x + 1) * (x + 1 + 1 / 10 ^ 30)));
        let proof = ScalarPoleCertificate::check(&before, &after).unwrap();
        assert!(!proof.scalar_zero);
        let failure = ScalarPoleCertificate::check(&before, &Atom::Zero).unwrap_err();
        assert!(matches!(
            failure.downcast_ref::<ScalarPoleFailure>(),
            Some(ScalarPoleFailure::Failed(_))
        ));
        assert!(ScalarPoleCertificate::check(&parse_lit!(1 / x), &parse_lit!(2 / x)).is_err());
        let zero = parse_lit!((x + y) ^ 2 - x ^ 2 - 2 * x * y - y ^ 2);
        let proof = ScalarPoleCertificate::check(&(zero / parse_lit!(s)), &Atom::Zero).unwrap();
        assert!(proof.scalar_zero);
        assert_eq!(proof.removed_factor, Atom::num(1));
        assert!(!proof.numerator_divisibility_verified);
    }

    #[test]
    fn floating_and_undefined_coefficients_are_unsupported() {
        for expression in [
            Atom::num(0.125f64),
            Atom::num(Coefficient::Indeterminate),
            Atom::num(Coefficient::positive_infinity()),
        ] {
            let failure = ScalarPoleCertificate::check(&expression, &expression).unwrap_err();
            assert!(matches!(
                failure.downcast_ref::<ScalarPoleFailure>(),
                Some(ScalarPoleFailure::Unsupported(_))
            ));
        }
    }

    #[test]
    fn owned_energies_keep_their_full_arguments() {
        crate::initialisation::test_initialise().unwrap();
        let argument = parse_lit!(m ^ 2 + (q + t) ^ 2);
        let energy = function!(GS.energy_surface, 3, &argument);
        let s = parse_lit!(s);
        let before = (energy.pow(-1) - (&energy + &s).pow(-1)) / &s;
        let after = (&energy * (&energy + &s)).pow(-1);
        let proof = ScalarPoleCertificate::check(&before, &after).unwrap();
        assert_eq!(proof.removed_factor, s);
        assert!(proof.opaque_leaves.contains(&energy));
        assert!(!proof.opaque_leaves.contains(&argument));
        assert!(ScalarPoleCertificate::check(&energy.pow(2), &argument).is_err());
        assert!(
            ScalarPoleCertificate::check(&energy, &function!(GS.energy_surface, 7, &argument))
                .is_err()
        );
    }

    #[test]
    fn unsupported_inputs_and_budgets_are_unproven() {
        for expression in [
            parse_lit!(x ^ (1 / 2)),
            parse_lit!(x ^ 100000),
            parse_lit!((x + y) ^ 64),
            function!(
                symbol!("gammalooprs::uv::numerator_family"),
                0,
                parse_lit!(q)
            ),
            Symbol::IF.call_args([parse_lit!(q), parse_lit!(1 / q), Atom::Zero]),
        ] {
            assert!(ScalarPoleCertificate::check(&expression, &expression).is_err());
        }
    }
}
