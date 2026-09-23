use std::collections::BTreeMap;

use symbolica::{
    atom::{Atom, AtomCore},
    tensors::matrix::Matrix,
};

use super::{IntegralFamily, IntegralFamilyError};

/// An affine relation `D_pivot - sum(c_i D_i) = constant`, with `i < pivot`.
#[derive(Clone)]
struct Dependence {
    pivot: usize,
    coefficients: Vec<(usize, Atom)>,
    constant: Atom,
}

impl IntegralFamily {
    /// Decompose a product into terms with independent positive-power propagators.
    ///
    /// Each returned pair contains a coefficient and signed propagator powers
    /// in the original family order. Negative powers remain numerator factors.
    /// The identity is algebraic: no loop shifts or scaleless-term removal are
    /// performed. Generic nonzero external invariants may appear in coefficients;
    /// specialize kinematics before constructing the family for degenerate cases.
    ///
    /// `max_states` bounds processed exponent vectors, including intermediate
    /// states. Exceeding the budget returns an error rather than a partial result.
    pub fn partial_fraction(
        &self,
        powers: &[i32],
        max_states: usize,
    ) -> Result<Vec<(Atom, Vec<i32>)>, IntegralFamilyError> {
        if powers.len() != self.denominators.len() {
            return Err(IntegralFamilyError::InvalidPowers);
        }
        if self
            .denominators
            .iter()
            .zip(powers)
            .any(|(d, p)| d.is_zero() && *p > 0)
        {
            return Err(IntegralFamilyError::InvalidBasis(
                "a positive-power denominator is identically zero".into(),
            ));
        }
        let mut pending = BTreeMap::from([(powers.to_vec(), Atom::one())]);
        let mut relations = BTreeMap::new();
        let mut result = Vec::new();
        let mut processed = 0;
        while let Some((powers, coefficient)) = pending.pop_last() {
            processed += 1;
            if processed > max_states {
                return Err(IntegralFamilyError::PartialFractionLimit(max_states));
            }
            let coefficient = coefficient.together();
            if coefficient.is_zero() {
                continue;
            }
            let support = powers
                .iter()
                .enumerate()
                .filter_map(|(i, p)| (*p > 0).then_some(i))
                .collect::<Vec<_>>();
            if !relations.contains_key(&support) {
                relations.insert(support.clone(), self.first_dependence(&support)?);
            }
            let Some(relation) = &relations[&support] else {
                result.push((coefficient, powers));
                continue;
            };
            let branches = if relation.constant.is_zero() {
                if relation.coefficients.is_empty() {
                    return Err(IntegralFamilyError::InvalidBasis(
                        "a positive-power denominator is identically zero".into(),
                    ));
                }
                // D_j = sum(c_i D_i). Dividing by D_j transfers one power
                // from an earlier denominator to j. Total degree is fixed,
                // while the exponent vector decreases lexicographically.
                relation
                    .coefficients
                    .iter()
                    .map(|(i, c)| (*i, c.clone()))
                    .collect::<Vec<_>>()
            } else {
                // 1 = (D_j - sum(c_i D_i))/constant reduces total degree.
                std::iter::once((relation.pivot, Atom::one() / &relation.constant))
                    .chain(
                        relation
                            .coefficients
                            .iter()
                            .map(|(i, c)| (*i, -c / &relation.constant)),
                    )
                    .collect()
            };
            for (removed, factor) in branches {
                if factor.is_zero() {
                    continue;
                }
                let mut next = powers.clone();
                next[removed] -= 1;
                if relation.constant.is_zero() {
                    next[relation.pivot] = next[relation.pivot]
                        .checked_add(1)
                        .ok_or(IntegralFamilyError::PowerOverflow)?;
                }
                // Processing largest vectors first coalesces all paths to a
                // state before visiting it, and avoids recursive stack growth.
                debug_assert!(next < powers);
                *pending.entry(next).or_insert_with(Atom::new) += &coefficient * factor;
                if processed.saturating_add(pending.len()) > max_states {
                    return Err(IntegralFamilyError::PartialFractionLimit(max_states));
                }
            }
        }
        Ok(result)
    }

    fn first_dependence(
        &self,
        support: &[usize],
    ) -> Result<Option<Dependence>, IntegralFamilyError> {
        if support.is_empty() {
            return Ok(None);
        }
        let denominators = support
            .iter()
            .map(|i| &self.denominators[*i])
            .collect::<Vec<_>>();
        let (matrix, _) = Self::affine_system(&denominators, &self.scalar_products)
            .map_err(|e| IntegralFamilyError::InvalidBasis(e.to_string()))?;
        let mut rows = Vec::new();
        for (position, row) in matrix.row_iter().enumerate() {
            rows.push(row.to_vec());
            let candidate = Matrix::from_nested_vec(rows.clone(), matrix.field().clone())
                .map_err(IntegralFamilyError::InvalidBasis)?;
            if candidate.rank() == rows.len() {
                continue;
            }
            rows.pop();
            let coefficients = if rows.is_empty() {
                Vec::new()
            } else {
                let basis = Matrix::from_nested_vec(rows, matrix.field().clone())
                    .map_err(IntegralFamilyError::InvalidBasis)?
                    .transpose();
                let target = Matrix::new_vec(row.to_vec(), matrix.field().clone());
                basis
                    .solve(&target)
                    .map_err(|e| IntegralFamilyError::InvalidBasis(e.to_string()))?
                    .into_vec()
                    .into_iter()
                    .enumerate()
                    .map(|(i, c)| (support[i], c.to_expression()))
                    .filter(|(_, c)| !c.is_zero())
                    .collect::<Vec<_>>()
            };
            let pivot = support[position];
            let constant = (&self.denominators[pivot]
                - coefficients
                    .iter()
                    .map(|(i, c)| c * &self.denominators[*i])
                    .sum::<Atom>())
            .together();
            return Ok(Some(Dependence {
                pivot,
                coefficients,
                constant,
            }));
        }
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use feynkit_kinematics::Kinematics;
    use symbolica::parse;

    fn verify(
        denominators: Vec<Atom>,
        powers: Vec<i32>,
        loops: Vec<Atom>,
        external: Vec<Atom>,
        kin: &Kinematics,
    ) {
        let family = IntegralFamily::new(loops, external, denominators.clone(), kin).unwrap();
        let terms = family.partial_fraction(&powers, 10_000).unwrap();
        let expression = |powers: &[i32]| {
            denominators
                .iter()
                .zip(powers)
                .map(|(d, p)| d.clone().pow(-i64::from(*p)))
                .product::<Atom>()
        };
        let mut reconstructed = Atom::Zero;
        for (coefficient, term_powers) in terms {
            let active = denominators
                .iter()
                .zip(&term_powers)
                .filter_map(|(d, p)| (*p > 0).then_some(d.clone()))
                .collect::<Vec<_>>();
            assert_eq!(family.rank_of(&active).unwrap(), active.len());
            reconstructed += coefficient * expression(&term_powers);
        }
        assert!((reconstructed - expression(&powers)).together().is_zero());
    }

    #[test]
    fn mass_shifts_and_repeated_propagators_preserve_the_rational_function() {
        let k = parse!("apart_mass::k");
        let kin = Kinematics::new();
        let x = kin.scalar_product(&k, &k).unwrap();
        verify(
            vec![x.clone(), &x - parse!("m2")],
            vec![2, 3],
            vec![k.clone()],
            vec![],
            &kin,
        );
        verify(
            vec![x.clone(), Atom::num(2) * &x, &x - parse!("m2")],
            vec![1, 2, 1],
            vec![k.clone()],
            vec![],
            &kin,
        );
        verify(
            vec![Atom::num(7), x.clone(), x],
            vec![2, 1, 1],
            vec![k],
            vec![],
            &kin,
        );
    }

    #[test]
    fn homogeneous_and_inhomogeneous_multivariate_relations() {
        let k = parse!("apart_multi::k");
        let p = parse!("apart_multi::p");
        let kin = Kinematics::new();
        let x = kin.scalar_product(&k, &k).unwrap();
        let y = kin.scalar_product(&k, &p).unwrap();
        for constant in [Atom::Zero, parse!("s")] {
            verify(
                vec![x.clone(), y.clone(), &x + &y + constant, &x - &y],
                vec![2, 1, 2, 1],
                vec![k.clone()],
                vec![p.clone()],
                &kin,
            );
        }
        verify(
            vec![x.clone(), y.clone(), &x + &y],
            vec![-2, 1, 1],
            vec![k],
            vec![p],
            &kin,
        );
    }

    #[test]
    fn budgets_and_invalid_powers_are_explicit() {
        let k = parse!("apart_limit::k");
        let kin = Kinematics::new();
        let x = kin.scalar_product(&k, &k).unwrap();
        let family =
            IntegralFamily::new(vec![k.clone()], vec![], vec![x.clone(), &x - 1], &kin).unwrap();
        assert!(matches!(
            family.partial_fraction(&[1], 10),
            Err(IntegralFamilyError::InvalidPowers)
        ));
        assert!(matches!(
            family.partial_fraction(&[2, 2], 1),
            Err(IntegralFamilyError::PartialFractionLimit(1))
        ));
        let zero = IntegralFamily::new(vec![k], vec![], vec![Atom::Zero], &kin).unwrap();
        assert!(zero.partial_fraction(&[1], 10).is_err());
        assert_eq!(
            zero.partial_fraction(&[0], 10).unwrap(),
            vec![(Atom::one(), vec![0])]
        );
    }

    #[test]
    fn fractional_propagator_coefficients_reconstruct_exactly() {
        let k = parse!("apart_fraction::k");
        let kin = Kinematics::new();
        let square = kin.scalar_product(&k, &k).unwrap();
        verify(
            vec![&square / 2 - parse!("m2"), &square / 3 - parse!("M2")],
            vec![2, 1],
            vec![k],
            vec![],
            &kin,
        );
    }
}
