//! Exact scalar specialization shared by native model consumers.

use std::collections::BTreeMap;

use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    domains::rational::Rational,
    id::{Pattern, Replacement},
    symbol,
};

use crate::{Model, ModelError, ParameterCard, ParameterNature};

/// Close a set of named scalar definitions using literal Symbolica replacements.
///
/// Keys are literal symbols, including names ending in underscores. Values may
/// retain constants or symbols outside this map; this operation resolves only
/// dependencies among its keys. It neither numerically evaluates nor expands
/// the expressions. Cyclic or unresolved dependencies between keys are errors.
pub fn resolve_scalar_bindings(
    mut values: BTreeMap<Symbol, Atom>,
) -> Result<BTreeMap<Symbol, Atom>, ModelError> {
    for _ in 0..=values.len() {
        let replacements = values
            .iter()
            .map(|(symbol, value)| {
                Replacement::new(
                    Pattern::Literal(Atom::var(*symbol)),
                    Pattern::Literal(value.clone()),
                )
            })
            .collect::<Vec<_>>();
        let next = values
            .iter()
            .map(|(symbol, value)| (*symbol, value.replace_multiple(&replacements)))
            .collect::<BTreeMap<_, _>>();
        if next == values {
            break;
        }
        values = next;
    }
    if values.values().any(|value| {
        values
            .keys()
            .any(|symbol| value.contains(Atom::var(*symbol).as_view()))
    }) {
        return Err(ModelError::UnresolvedScalarBindings);
    }
    Ok(values)
}

impl Model {
    /// Return exact scalar bindings for this model without changing it.
    ///
    /// The optional card is applied to a copy. External parameters and explicit
    /// card values become exact complex rationals representing their binary64
    /// values. Internal parameters retain their analytic definitions unless
    /// explicitly supplied by the card; cached dependent values are never used
    /// instead of those definitions. Couplings remain analytic until the named
    /// dependencies are substituted, without floating-point recomputation.
    ///
    /// Exact `overrides` take precedence over both model and card definitions
    /// before closure. Additional scalar keys are allowed. The returned map has
    /// no dependencies on its own keys, but may retain other symbols/functions;
    /// it is not an evaluator for arbitrary UFO functions or form factors.
    ///
    /// Pass the restriction card explicitly even for a model to which it was
    /// previously applied: the model does not retain internal-override history.
    pub fn scalar_bindings(
        &self,
        card: Option<&ParameterCard>,
        overrides: &BTreeMap<Symbol, Atom>,
    ) -> Result<BTreeMap<Symbol, Atom>, ModelError> {
        let mut updated;
        let model = if let Some(card) = card {
            updated = self.clone();
            updated.apply_parameter_card(card)?;
            &updated
        } else {
            self
        };
        let mut values = BTreeMap::new();
        for parameter in model.parameters() {
            let analytic = parameter.expression.as_ref().filter(|_| {
                parameter.nature == ParameterNature::Internal
                    && !card.is_some_and(|card| card.contains_key(&parameter.name))
            });
            let value = if let Some(expression) = analytic {
                expression.clone()
            } else if let Some(value) = parameter.value {
                let exact = |number| {
                    Rational::try_from(number).map_err(|_| ModelError::NonFiniteScalarBinding {
                        name: parameter.name.clone(),
                    })
                };
                Atom::num(exact(value.re)?) + Atom::num(exact(value.im)?) * Atom::i()
            } else if let Some(expression) = &parameter.expression {
                expression.clone()
            } else {
                continue;
            };
            values.insert(symbol!(&format!("UFO::{}", parameter.name)), value);
        }
        for coupling in model.couplings() {
            values.insert(
                symbol!(&format!("UFO::{}", coupling.name)),
                coupling.expression.clone(),
            );
        }
        values.extend(
            overrides
                .iter()
                .map(|(symbol, value)| (*symbol, value.clone())),
        );
        resolve_scalar_bindings(values)
    }
}

#[cfg(test)]
mod tests;
