//! Emit the shared trace recurrence as an expanded scalar polynomial.
//! FactoredTrace stores products and sums until one coefficient-list emission.

use super::{TraceAlgebra, factored::FactoredTrace, is_four_dimension};
use std::collections::HashMap;
use symbolica::{
    atom::{Atom, AtomCore, AtomView},
    domains::integer::{Integer, Z},
    poly::{PolyVariable, polynomial::MultivariatePolynomial},
};

// Bound the dense exponent rows and in-Vec coefficient slots. This excludes
// Integer heap payloads, the recipe, normalized variables and sorting scratch.
// The ordinary free 14-factor trace (135135 rows, 91 variables) still fits.
const MAX_ROW_BYTES: usize = 16 * 1024 * 1024;

pub(super) struct SparseOutput {
    variables: Vec<PolyVariable>,
    positions: HashMap<Atom, usize>,
    recipe: FactoredTrace,
    trace_unit: Option<Integer>,
}
impl Default for SparseOutput {
    fn default() -> Self {
        Self {
            variables: Vec::new(),
            positions: HashMap::new(),
            recipe: FactoredTrace::unit_recipe(),
            trace_unit: None,
        }
    }
}
impl SparseOutput {
    fn variable(&mut self, atom: AtomView<'_>) -> usize {
        if let Some(&position) = self.positions.get(&atom.to_owned()) {
            return position;
        }
        let position = self.variables.len();
        let atom = atom.to_owned();
        self.variables.push(atom.clone().try_into().unwrap());
        self.positions.insert(atom, position);
        position
    }

    fn exponent_capacity(count: usize, variables: usize, row_budget: usize) -> Option<usize> {
        let exponents = count.checked_mul(variables)?;
        let coefficients = count.checked_mul(std::mem::size_of::<Integer>())?;
        (exponents.checked_add(coefficients)? <= row_budget).then_some(exponents)
    }

    fn emit(self, value: usize, row_budget: usize) -> Atom {
        let Some(unit) = self.trace_unit else {
            return Atom::Zero;
        };
        let allocation = self.recipe.recipe_leaf_count(value).and_then(|count| {
            Self::exponent_capacity(count, self.variables.len(), row_budget)
                .map(|exponents| (count, exponents))
        });
        let Some((count, exponents)) = allocation else {
            // Reuse the completed recurrence and its already normalized atoms.
            // A refused allocation must not construct the user's metrics twice.
            let mut factors = vec![Atom::Zero; self.variables.len()];
            for (atom, position) in self.positions {
                factors[position] = atom;
            }
            return self
                .recipe
                .evaluate(value, &factors, Atom::num(unit).as_view())
                .expand();
        };
        if count == 0 {
            return Atom::Zero;
        }
        let mut coefficients = Vec::with_capacity(count);
        let mut exponents = Vec::with_capacity(exponents);
        let mut monomial = vec![0u8; self.variables.len()];
        self.recipe
            .visit_recipe_leaves(value, &unit, &mut |factors, coefficient| {
                for &factor in factors {
                    monomial[factor] += 1;
                }
                coefficients.push(coefficient.clone());
                exponents.extend_from_slice(&monomial);
                for &factor in factors {
                    monomial[factor] -= 1;
                }
            });
        MultivariatePolynomial::<_, u8>::from_coefficient_list(
            coefficients,
            exponents,
            self.variables.into(),
            &Z,
        )
        .to_expression()
    }
}
impl TraceAlgebra for SparseOutput {
    type Value = usize;
    fn unit(&mut self, trace_unit: AtomView<'_>) -> usize {
        self.trace_unit
            .get_or_insert_with(|| Integer::from(i64::try_from(trace_unit).unwrap()));
        0
    }
    fn scale(&mut self, coefficient: i32, value: usize) -> usize {
        self.recipe.push_product(Vec::new(), coefficient, value)
    }
    fn metric_product(&mut self, metric: AtomView<'_>, coefficient: i32, value: usize) -> usize {
        let variable = self.variable(metric);
        self.recipe.push_product(vec![variable], coefficient, value)
    }
    fn dimension_product(
        &mut self,
        dimension: AtomView<'_>,
        constant: i32,
        linear: i32,
        value: usize,
    ) -> usize {
        if is_four_dimension(dimension) {
            return self.scale(constant + 4 * linear, value);
        }
        let variable = self.variable(dimension);
        let mut terms = Vec::with_capacity(2);
        if constant != 0 {
            terms.push(self.recipe.push_product(Vec::new(), constant, value));
        }
        if linear != 0 {
            terms.push(self.recipe.push_product(vec![variable], linear, value));
        }
        self.recipe.push_sum(terms)
    }
    fn add(&mut self, left: usize, right: usize) -> usize {
        self.recipe.push_sum(vec![left, right])
    }
    fn subtract(&mut self, left: usize, right: usize) -> usize {
        let right = self.scale(-1, right);
        self.add(left, right)
    }
    fn sum(&mut self, values: Vec<usize>) -> usize {
        self.recipe.push_sum(values)
    }
    fn finish(self, value: usize) -> Atom {
        self.emit(value, MAX_ROW_BYTES)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn row_allocation_checks_overflow_and_preserves_the_free_fourteen_bound() {
        assert_eq!(
            SparseOutput::exponent_capacity(135135, 91, MAX_ROW_BYTES),
            Some(135135 * 91)
        );
        let bytes = 3 * (2 + std::mem::size_of::<Integer>());
        assert_eq!(SparseOutput::exponent_capacity(3, 2, bytes), Some(6));
        assert_eq!(SparseOutput::exponent_capacity(3, 2, bytes - 1), None);
        assert_eq!(
            SparseOutput::exponent_capacity(2, usize::MAX, usize::MAX),
            None
        );
        assert_eq!(
            SparseOutput::exponent_capacity(usize::MAX, 0, usize::MAX),
            None
        );
    }

    #[test]
    fn refused_rows_reuse_the_requested_recipe_root_and_trace_unit() {
        let x = Atom::var(symbolica::symbol!("sparse_budget_x"));
        let d = Atom::var(symbolica::symbol!("sparse_budget_d"));
        let expected = (Atom::num(4)
            * (Atom::num(2) * &x * &x * (&d - Atom::num(4)) + Atom::num(3) * &d))
            .expand();
        for budget in [0, MAX_ROW_BYTES] {
            let mut output = SparseOutput::default();
            let unit = output.unit(Atom::num(4).as_view());
            let x2 = output.metric_product(x.as_view(), 2, unit);
            let x2 = output.metric_product(x.as_view(), 1, x2);
            let left = output.dimension_product(d.as_view(), -4, 1, x2);
            let right = output.metric_product(d.as_view(), 3, unit);
            let root = output.add(left, right);
            // Runtime roots need not equal the stored template root or last node.
            output.scale(7, unit);
            assert_eq!(output.emit(root, budget), expected);
        }
        let mut output = SparseOutput::default();
        let unit = output.unit(Atom::num(4).as_view());
        let zero = output.scale(0, unit);
        assert_eq!(output.emit(zero, 0), Atom::Zero);
    }

    #[test]
    fn overflowing_leaf_count_emits_the_exact_small_polynomial() {
        let mut output = SparseOutput::default();
        let mut root = output.unit(Atom::num(4).as_view());
        for _ in 0..usize::BITS {
            root = output.add(root, root);
        }
        assert_eq!(output.recipe.recipe_leaf_count(root), None);
        let expected = Atom::num(4) * Atom::num(2).as_view().pow(i64::from(usize::BITS));
        assert_eq!(output.finish(root), expected);
    }
}
