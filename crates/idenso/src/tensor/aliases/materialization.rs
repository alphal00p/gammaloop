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
        let mut variables = HashSet::new();
        let mut leaves = HashSet::new();
        for source in definitions.values().chain([self.expression.get_root()]) {
            if !InterfaceInference::normalization_is_intrinsic(source.as_view())
                || !Self::polynomial_leaves(source.as_view(), definitions, &mut leaves)
            {
                return Ok(None);
            }
        }
        for leaf in &leaves {
            if !matches!(leaf.as_view(), AtomView::Num(_)) {
                let Ok(variable) = PolyVariable::try_from(leaf.clone()) else {
                    return Ok(None);
                };
                variables.insert(variable);
            }
        }
        let mut variables = variables.into_iter().collect::<Vec<_>>();
        variables.sort_unstable();
        let variables = Arc::new(variables);
        let mut values = HashMap::new();
        for leaf in leaves {
            let Ok(value) = leaf.try_to_polynomial::<_, u32>(&Q, Some(variables.clone())) else {
                return Ok(None);
            };
            values.insert(leaf, value);
        }
        for handle in order {
            let Some(value) = Self::polynomial_value(definitions[handle].as_view(), &values) else {
                return Ok(None);
            };
            let value = value.into_owned();
            values.insert(handle.clone(), value);
        }
        Ok(
            Self::polynomial_value(self.expression.get_root().as_view(), &values)
                .map(|value| value.to_expression()),
        )
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

    fn polynomial_value<'a>(
        source: AtomView<'_>,
        values: &'a HashMap<Atom, Polynomial>,
    ) -> Option<Cow<'a, Polynomial>> {
        if let Some(value) = values.get(source.get_data()) {
            return Some(Cow::Borrowed(value));
        }
        match source {
            AtomView::Add(sum) => {
                let mut terms = sum.iter();
                let mut value = Self::polynomial_value(terms.next()?, values)?;
                for term in terms {
                    value =
                        Cow::Owned(value.as_ref() + Self::polynomial_value(term, values)?.as_ref());
                }
                Some(value)
            }
            AtomView::Mul(product) => {
                let mut factors = product.iter();
                let mut value = Self::polynomial_value(factors.next()?, values)?;
                for factor in factors {
                    value = Cow::Owned(
                        value.as_ref() * Self::polynomial_value(factor, values)?.as_ref(),
                    );
                }
                Some(value)
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
}
