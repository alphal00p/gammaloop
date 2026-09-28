//! Shared scope and explicit-index observation. This is a syntax-only pass:
//! no normalization callback, materialization, or dummy allocation is permitted.

use super::*;
use spenso::network::parsing::ChainNestingError;

#[derive(Clone, Copy, Default, Eq, Hash, PartialEq)]
pub(in crate::tensor) struct ObservationScope {
    pub(super) owner: Option<Symbol>,
    pub(super) inside_factor: bool,
    pub(super) validate: bool,
    pub(super) count_indices: bool,
}

type Cache = HashMap<ObservationScope, HashMap<Vec<u8>, ExplicitIndexOccurrences>>;

impl InterfaceInference {
    pub(in crate::tensor) fn index_occurrences(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        cache: &mut Cache,
    ) -> ExplicitIndexOccurrences {
        Self::observe_indices(
            value,
            slots,
            cache,
            ObservationScope {
                count_indices: true,
                ..ObservationScope::default()
            },
        )
        .expect("index-only observation cannot reject a placeholder scope")
    }

    pub(super) fn scoped_occurrences(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
    ) -> InferenceResult<ExplicitIndexOccurrences> {
        #[cfg(test)]
        super::tests::SCOPE_VALIDATIONS.with(|count| count.set(count.get() + 1));
        Self::observe_indices(
            value,
            slots,
            &mut HashMap::new(),
            ObservationScope {
                validate: true,
                count_indices: true,
                ..ObservationScope::default()
            },
        )
    }

    pub(super) fn observe_indices(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        cache: &mut Cache,
        scope: ObservationScope,
    ) -> InferenceResult<ExplicitIndexOccurrences> {
        let reusable = matches!(value, AtomView::Fun(_))
            && !value.needs_normalization()
            && value.get_byte_size() <= Self::CACHE_KEY_BYTES;
        if reusable
            && let Some(summary) = cache
                .get(&scope)
                .and_then(|values| values.get(value.get_data()))
        {
            return Ok(summary.clone());
        }
        let result = Self::observe_indices_uncached(value, slots, cache, scope)?;
        if reusable && cache.values().map(HashMap::len).sum::<usize>() < Self::CACHE_ENTRIES {
            // Scope is part of the key: a placeholder proven inside a chain is
            // not legal when the same literal subtree occurs outside it.
            cache
                .entry(scope)
                .or_default()
                .insert(value.get_data().to_vec(), result.clone());
        }
        Ok(result)
    }

    fn observe_indices_uncached(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        cache: &mut Cache,
        mut scope: ObservationScope,
    ) -> InferenceResult<ExplicitIndexOccurrences> {
        #[cfg(test)]
        if scope.count_indices && matches!(value, AtomView::Fun(_)) {
            super::tests::OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(count.get() + 1));
        }
        let mut result = ExplicitIndexOccurrences::default();
        if scope.count_indices
            && let Ok(slot) = slots.parse::<LibraryRep, AbstractIndex>(value)
        {
            result.add(slot.rep().base().slot(slot.aind()), 1);
            // Index payloads are opaque to multiplicity, but the scope contract
            // still checks placeholders and nested chains inside those payloads.
            scope.count_indices = false;
        }
        if !scope.count_indices && !scope.validate {
            return Ok(result);
        }
        match value {
            AtomView::Var(_) if scope.validate && Self::is_chain_placeholder(value) => {
                return Err(TensorInferenceError::invalid(
                    "Spenso chain placeholders are only valid inside chain or trace tensor factors",
                ));
            }
            AtomView::Add(sum) => {
                for term in sum.iter() {
                    let branch = Self::observe_indices(term, slots, cache, scope)?;
                    for (slot, count) in branch.0 {
                        let maximum = result.0.entry(slot).or_default();
                        *maximum = (*maximum).max(count);
                    }
                }
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    result.append(Self::observe_indices(factor, slots, cache, scope)?);
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                // Preserve the occurrence owner's branch bound. Tensor power
                // multiplicity remains the separate interface inference check.
                for term in [base, exponent] {
                    let branch = Self::observe_indices(term, slots, cache, scope)?;
                    for (slot, count) in branch.0 {
                        let maximum = result.0.entry(slot).or_default();
                        *maximum = (*maximum).max(count);
                    }
                }
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                scope.validate &= !symbol.is_scalar();
                if scope.validate {
                    scope.owner = ChainNestingError::enter(scope.owner, symbol)
                        .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
                }
                let tensor_leaf = symbol.has_tag(&SPENSO_TAG.tensor);
                let count_arguments = scope.count_indices
                    && (tensor_leaf || composition::is_composite_head(symbol, &SPENSO_TAG));
                let mut inputs = 0;
                let mut outputs = 0;
                for (position, argument) in function.iter().enumerate() {
                    let mut child = scope;
                    if symbol == SPENSO_TAG.chain {
                        child.inside_factor = position >= 2;
                    } else if symbol == SPENSO_TAG.trace {
                        child.inside_factor = position >= 1;
                    } else if scope.validate && Self::is_chain_placeholder(argument) {
                        if !(scope.inside_factor && tensor_leaf) {
                            return Err(TensorInferenceError::invalid(
                                "Spenso chain placeholders are only valid as tensor ports inside chain or trace factors",
                            ));
                        }
                        let AtomView::Var(variable) = argument else {
                            unreachable!()
                        };
                        if variable.get_symbol() == SPENSO_TAG.chain_in {
                            inputs += 1;
                        } else {
                            outputs += 1;
                        }
                        continue;
                    }
                    child.count_indices = count_arguments
                        && (slots.parse::<LibraryRep, AbstractIndex>(argument).is_ok()
                            || Representation::<LibraryRep>::try_from(argument).is_ok()
                            || argument.is_tensorial(StrictTensorFilter::Tagged)
                            || matches!(argument, AtomView::Fun(nested)
                                if nested.get_symbol() == *shadowing::SYM
                                || nested.get_symbol() == *shadowing::ANTISYM
                                || nested.get_symbol() == *shadowing::CYCLIC));
                    if child.validate || child.count_indices {
                        result.append(Self::observe_indices(argument, slots, cache, child)?);
                    }
                }
                if inputs + outputs > 0 && (inputs != 1 || outputs != 1) {
                    return Err(TensorInferenceError::invalid(format!(
                        "tensor factor `{symbol}` requires exactly one direct `in` and one direct `out` placeholder",
                    )));
                }
            }
            _ => {}
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use spenso::structure::representation::ExtendibleReps;

    #[test]
    fn one_observation_checks_scopes_and_keeps_branch_multiplicity() {
        crate::test_support::test_initialize();
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let a = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99401));
        let leaf = |name| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(a.to_atom())
                .finish()
        };
        let alternatives = leaf("observation_A") + leaf("observation_B");
        let valid = &alternatives * leaf("observation_C");
        let invalid = valid.clone() * leaf("observation_D");
        for (value, expected) in [(alternatives, 1), (valid, 2), (invalid, 3)] {
            let occurrences = InterfaceInference::scoped_occurrences(
                value.as_view(),
                &mut SlotMatcher::default(),
            )
            .unwrap();
            assert_eq!(occurrences.0[&a], expected);
            assert_eq!(occurrences.validate().is_ok(), expected <= 2);
        }
    }

    #[test]
    fn cached_scope_proof_is_local_and_scalar_metadata_is_opaque() {
        crate::test_support::test_initialize();
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let input = Atom::var(SPENSO_TAG.chain_in);
        let output = Atom::var(SPENSO_TAG.chain_out);
        let leaf = FunctionBuilder::new(SPENSO_TAG.tensor_symbol("observation_channel"))
            .add_arg(&input)
            .add_arg(&output)
            .finish();
        let mut slots = SlotMatcher::default();
        let mut cache = Cache::new();
        let in_factor = ObservationScope {
            validate: true,
            inside_factor: true,
            count_indices: true,
            ..ObservationScope::default()
        };
        assert!(
            InterfaceInference::observe_indices(leaf.as_view(), &mut slots, &mut cache, in_factor)
                .is_ok()
        );
        assert!(
            InterfaceInference::observe_indices(
                leaf.as_view(),
                &mut slots,
                &mut cache,
                ObservationScope {
                    inside_factor: false,
                    ..in_factor
                }
            )
            .is_err()
        );
        let chain = FunctionBuilder::new(SPENSO_TAG.chain)
            .add_arg(rep.slot::<AbstractIndex, _>(99403).to_atom())
            .add_arg(rep.slot::<AbstractIndex, _>(99405).to_atom())
            .add_arg(leaf)
            .finish();
        assert!(InterfaceInference::scoped_occurrences(chain.as_view(), &mut slots).is_ok());
        let nested = FunctionBuilder::new(SPENSO_TAG.trace)
            .add_arg(rep.to_symbolic([]))
            .add_arg(chain)
            .finish();
        assert!(InterfaceInference::scoped_occurrences(nested.as_view(), &mut slots).is_err());
        let scalar = FunctionBuilder::new(symbolica::symbol!("observation_metadata"; Scalar))
            .add_arg(nested)
            .add_arg(input)
            .finish();
        assert!(
            InterfaceInference::scoped_occurrences(scalar.as_view(), &mut slots)
                .unwrap()
                .0
                .is_empty()
        );
    }
}
