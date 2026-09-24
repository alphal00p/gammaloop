//! Emit the shared trace recurrence as an expanded scalar polynomial.
//! FactoredTrace stores products and sums until one coefficient-list emission.

use super::{TraceAlgebra, factored::FactoredTrace, is_four_dimension};
use std::collections::HashMap;
use symbolica::{
    atom::{Atom, AtomView},
    domains::integer::{Integer, Z},
    poly::{PolyVariable, polynomial::MultivariatePolynomial},
};

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
        let Some(unit) = self.trace_unit else {
            return Atom::Zero;
        };
        let count = self.recipe.recipe_leaf_count(value);
        let mut coefficients = Vec::with_capacity(count);
        let mut exponents = Vec::with_capacity(count * self.variables.len());
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
