//! Explicit polynomial materialization of the existing alias registry.

use std::{borrow::Cow, sync::Arc};

use symbolica::{
    domains::rational::Q,
    poly::{PolyVariable, polynomial::MultivariatePolynomial},
};

use super::*;

type Polynomial = MultivariatePolynomial<Q, u32>;

impl SymbolicTensor<AliasInterfaces, AliasedAtom> {
    pub(super) fn polynomial_materialization(&self) -> Result<Option<Atom>> {
        let order = self.dependency_order()?;
        let definitions = self.expression.get_aliases();
        let mut leaves = HashSet::new();
        for source in definitions.values().chain([self.expression.get_root()]) {
            if !InterfaceInference::normalization_is_intrinsic(source.as_view())
                || !Self::polynomial_leaves(source.as_view(), definitions, &mut leaves)
            {
                return Ok(None);
            }
        }
        let mut values = HashMap::new();
        for leaf in leaves {
            let Ok(value) = leaf.try_to_polynomial::<_, u32>(&Q, None) else {
                return Ok(None);
            };
            values.insert(leaf, value);
        }
        for handle in order {
            let Some(value) = Self::polynomial_definition(definitions[handle].as_view(), &values)
            else {
                return Ok(None);
            };
            let value = value.into_owned();
            values.insert(handle.clone(), value);
        }
        let root = self.expression.get_root().as_view();
        if let AtomView::Add(sum) = root {
            // The final result is an Atom. Each summand is already expanded;
            // merging them here needs no polynomial with their combined variables.
            let terms = sum
                .iter()
                .map(|term| {
                    Self::polynomial_definition(term, &values).map(|value| value.to_expression())
                })
                .collect::<Option<Vec<_>>>();
            Ok(terms.map(Atom::add_many))
        } else {
            Ok(Self::polynomial_definition(root, &values).map(|value| value.to_expression()))
        }
    }

    fn polynomial_leaves(
        source: AtomView<'_>,
        definitions: &ahash::HashMap<Atom, Atom>,
        leaves: &mut HashSet<Atom>,
    ) -> bool {
        if definitions.contains_key(source.get_data()) {
            return true;
        }
        match source {
            AtomView::Add(sum) => sum
                .iter()
                .all(|term| Self::polynomial_leaves(term, definitions, leaves)),
            AtomView::Mul(product) => product
                .iter()
                .all(|factor| Self::polynomial_leaves(factor, definitions, leaves)),
            AtomView::Pow(power)
                if i64::try_from(power.get_exp()).is_ok_and(|power| power >= 0) =>
            {
                Self::polynomial_leaves(power.get_base(), definitions, leaves)
            }
            _ => {
                // An alias inside a function or a non-polynomial power needs
                // Symbolica's ordinary literal resolution before expansion.
                let mut nested = false;
                source.visitor(&mut |node| {
                    nested |= definitions.contains_key(node.get_data());
                    !nested
                });
                if !nested {
                    leaves.insert(source.to_owned());
                }
                !nested
            }
        }
    }

    /// Establish one variable table per definition before doing its arithmetic.
    /// Cached definitions retain their smaller tables outside this invocation.
    fn polynomial_definition<'a>(
        source: AtomView<'_>,
        values: &'a HashMap<Atom, Polynomial>,
    ) -> Option<Cow<'a, Polynomial>> {
        if let Some(value) = values.get(source.get_data()) {
            return Some(Cow::Borrowed(value));
        }
        let mut needed = HashMap::new();
        source.visitor(&mut |node| {
            if let Some(value) = values.get(node.get_data()) {
                needed.insert(node.to_owned(), value);
                false
            } else {
                true
            }
        });
        let mut variables = needed
            .values()
            .flat_map(|value| value.get_vars_ref().iter().cloned())
            .collect::<HashSet<_>>()
            .into_iter()
            .collect::<Vec<PolyVariable>>();
        variables.sort_unstable();
        let variables = Arc::new(variables);
        let positions = variables
            .iter()
            .enumerate()
            .map(|(index, variable)| (variable, index))
            .collect::<HashMap<_, _>>();
        let mut local = HashMap::new();
        for (leaf, value) in needed {
            let polynomial = if value.get_vars_ref() == variables.as_slice() {
                Cow::Borrowed(value)
            } else {
                let map = value
                    .get_vars_ref()
                    .iter()
                    .map(|variable| positions[variable])
                    .collect::<Vec<_>>();
                let expanded = if map.windows(2).all(|positions| positions[0] < positions[1]) {
                    // Inserting zero columns preserves monomial order.
                    let mut expanded =
                        Polynomial::new(&Q, Some(value.nterms()), Arc::clone(&variables));
                    let mut exponents = vec![0; variables.len()];
                    for term in value {
                        for (&position, &power) in map.iter().zip(term.exponents) {
                            exponents[position] = power;
                        }
                        expanded.append_monomial_back(term.coefficient.clone(), &exponents);
                    }
                    expanded
                } else {
                    value.rearrange_with_growth(&variables).ok()?
                };
                Cow::Owned(expanded)
            };
            local.insert(leaf, polynomial);
        }
        Self::polynomial_value(source, &local).map(|value| Cow::Owned(value.into_owned()))
    }

    /// Merge the complete sum once using the definition's established variables.
    fn merge_polynomial_terms(terms: &[Cow<'_, Polynomial>]) -> Polynomial {
        if let [left, right] = terms {
            // Both inputs are sorted already, so a binary sum needs one merge.
            return left.as_ref() + right.as_ref();
        }
        let variables = terms[0].get_vars();
        let count = terms.iter().map(|term| term.nterms()).sum();
        let mut coefficients = Vec::with_capacity(count);
        let mut exponents = Vec::with_capacity(count * variables.len());
        for polynomial in terms {
            debug_assert_eq!(polynomial.get_vars_ref(), variables.as_slice());
            for term in polynomial.as_ref() {
                coefficients.push(term.coefficient.clone());
                exponents.extend_from_slice(term.exponents);
            }
        }
        if coefficients.is_empty() {
            Polynomial::new(&Q, None, variables)
        } else {
            Polynomial::from_coefficient_list(coefficients, exponents, variables, &Q)
        }
    }

    fn polynomial_value<'a>(
        source: AtomView<'_>,
        values: &'a HashMap<Atom, Cow<'_, Polynomial>>,
    ) -> Option<Cow<'a, Polynomial>> {
        if let Some(value) = values.get(source.get_data()) {
            return Some(Cow::Borrowed(value.as_ref()));
        }
        match source {
            AtomView::Add(sum) => {
                let terms = sum
                    .iter()
                    .map(|term| Self::polynomial_value(term, values))
                    .collect::<Option<Vec<_>>>()?;
                Some(Cow::Owned(Self::merge_polynomial_terms(&terms)))
            }
            AtomView::Mul(product) => {
                let mut factors = product
                    .iter()
                    .map(|factor| Self::polynomial_value(factor, values))
                    .collect::<Option<Vec<_>>>()?;
                factors.sort_unstable_by_key(|factor| factor.nterms());
                // Combine small factors first without growing one left-folded
                // variable table for every original factor.
                while factors.len() > 1 {
                    let mut next = Vec::with_capacity(factors.len().div_ceil(2));
                    let mut pairs = factors.into_iter();
                    while let Some(left) = pairs.next() {
                        next.push(if let Some(right) = pairs.next() {
                            Cow::Owned(left.as_ref() * right.as_ref())
                        } else {
                            left
                        });
                    }
                    factors = next;
                }
                factors.pop()
            }
            AtomView::Pow(power) => Some(Cow::Owned(
                Self::polynomial_value(power.get_base(), values)?
                    .pow(usize::try_from(i64::try_from(power.get_exp()).ok()?).ok()?),
            )),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::structure::partial::PartialStructureExt;
    use symbolica::{function, symbol};

    fn scalar(expression: Atom) -> SymbolicTensor<PartialStructure> {
        SymbolicTensor::checked_parts(expression, PartialStructure::from_logical_slots([])).unwrap()
    }

    #[test]
    fn polynomial_forward_pass_matches_nested_literal_resolution() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbol!("alias_materialize::x"));
        let y = Atom::var(symbol!("alias_materialize::y"));
        let first = scalar(&x + &y);
        let a = first.alias_handle().unwrap();
        let second = scalar(&a.expression * &a.expression + &x);
        let b = second.alias_handle().unwrap();
        let root = scalar((&b.expression + &a.expression) * (&b.expression - &y));
        let value = root.with_aliases([(a, first), (b, second)]).unwrap();
        let polynomial = value.polynomial_materialization().unwrap().unwrap();
        assert_eq!(polynomial, value.resolved().unwrap().expression.expand());
        assert_eq!(value.expanded().unwrap().expression, polynomial);
    }

    #[test]
    fn polynomial_materialization_defers_aliases_inside_opaque_functions() {
        crate::test_support::test_initialize();
        let body = scalar(Atom::var(symbol!("alias_materialize::z")) + 1);
        let handle = body.alias_handle().unwrap();
        let opaque = symbol!("alias_materialize::opaque"; Scalar);
        let root = scalar(function!(opaque, &handle.expression));
        let value = root.with_aliases([(handle, body)]).unwrap();
        assert!(value.polynomial_materialization().unwrap().is_none());
        assert_eq!(
            value.expanded().unwrap().expression,
            value.resolved().unwrap().expression.expand()
        );
    }

    #[test]
    fn polynomial_materialization_leaves_callback_order_to_symbolica() {
        crate::test_support::test_initialize();
        let callback = symbol!("alias_materialize::callback", norm = |_, _| {});
        let body = scalar(function!(
            callback,
            Atom::var(symbol!("alias_materialize::w"))
        ));
        let handle = body.alias_handle().unwrap();
        let value = handle.clone().with_aliases([(handle, body)]).unwrap();
        assert!(value.polynomial_materialization().unwrap().is_none());
        assert_eq!(value.expanded().unwrap(), value.resolved().unwrap());
    }
    #[test]
    fn materialization_merges_local_variables_and_cancelling_definitions() {
        crate::test_support::test_initialize();
        let x = Atom::var(symbol!("alias_materialize_local::x"));
        let y = Atom::var(symbol!("alias_materialize_local::y"));
        let z = Atom::var(symbol!("alias_materialize_local::z"));
        let first = scalar(&x + &z);
        let a = first.alias_handle().unwrap();
        let second = scalar(&y + &z);
        let b = second.alias_handle().unwrap();
        let duplicate = scalar(&x + &z);
        let c = scalar(Atom::var(symbol!(
            "alias_materialize_local::duplicate_handle"
        )));
        for root in [
            scalar(&a.expression - &c.expression),
            scalar(&a.expression * &b.expression + &c.expression),
            scalar((&a.expression - &c.expression) * &b.expression),
        ] {
            let value = root
                .with_aliases([
                    (a.clone(), first.clone()),
                    (b.clone(), second.clone()),
                    (c.clone(), duplicate.clone()),
                ])
                .unwrap();
            assert_eq!(
                value.expanded().unwrap().expression,
                value.resolved().unwrap().expression.expand()
            );
        }
    }
}
