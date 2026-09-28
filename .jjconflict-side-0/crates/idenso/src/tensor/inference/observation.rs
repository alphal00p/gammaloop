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
    pub(super) pair_representation: Option<LibraryRep>,
}

type Cache = HashMap<ObservationScope, HashMap<Vec<u8>, ExplicitIndexOccurrences>>;

impl InterfaceInference {
    /// Detect explicit index pairs in any additive branch without distributing it.
    ///
    /// Exclusions suppress ports of the named tensor leaves, without changing
    /// their payloads or recursing into scalar metadata. Closed powers and
    /// inverses keep their internal dummy scope; malformed multiplicities and
    /// placeholder scopes remain errors.
    pub fn has_explicit_index_pairs(
        value: AtomView<'_>,
        representation: LibraryRep,
        excluded_tensor_heads: &[Symbol],
    ) -> InferenceResult<bool> {
        let occurrences = Self::observe_indices(
            value,
            &mut SlotMatcher::default(),
            &mut HashMap::new(),
            ObservationScope {
                validate: true,
                count_indices: true,
                pair_representation: Some(representation),
                ..ObservationScope::default()
            },
            excluded_tensor_heads,
        )?;
        occurrences.validate()?;
        Ok(occurrences.selected_pair)
    }

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
            &[],
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
            &[],
        )
    }

    /// Keep only the boundary ports of one independently scoped expression.
    /// Observing its interface cannot normalize callbacks or allocate indices.
    fn scope_explicit_occurrences(
        value: AtomView<'_>,
        occurrences: &mut ExplicitIndexOccurrences,
        slots: &mut SlotMatcher,
    ) -> InferenceResult<PartialStructure> {
        let mut inference = InterfaceInference {
            slots: std::mem::take(slots),
            leaf_inference: LeafInference::Observe,
            ..InterfaceInference::default()
        };
        let interface = inference
            .infer_view(value)
            .and_then(|interface| Self::merge_explicit_interface_sequence(&[interface]));
        *slots = inference.slots;
        let interface = interface?;
        let boundary = interface.logical_slots();
        occurrences.counts.retain(|slot, _| {
            boundary.iter().any(|port| {
                port.rep().base() == slot.rep().base()
                    && port.aind == PartialIndex::Explicit(slot.aind())
            })
        });
        Ok(interface)
    }

    pub(super) fn observe_indices(
        value: AtomView<'_>,
        slots: &mut SlotMatcher,
        cache: &mut Cache,
        scope: ObservationScope,
        excluded_tensor_heads: &[Symbol],
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
        let result =
            Self::observe_indices_uncached(value, slots, cache, scope, excluded_tensor_heads)?;
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
        excluded_tensor_heads: &[Symbol],
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
                    let mut branch =
                        Self::observe_indices(term, slots, cache, scope, excluded_tensor_heads)?;
                    // Each alternative owns its internal contractions. Only
                    // surviving ports can meet factors outside this sum. Keep
                    // malformed local counts visible to the existing validator.
                    // Pair queries allow unlike branch interfaces: scope only
                    // proven interfaces there, without adding typed admission.
                    if !branch.counts.is_empty()
                        && branch.counts.values().all(|count| *count <= 2)
                        && let Err(error) =
                            Self::scope_explicit_occurrences(term, &mut branch, slots)
                        && scope.validate
                        && scope.pair_representation.is_none()
                    {
                        return Err(error);
                    }
                    result.selected_pair |= branch.selected_pair;
                    for (slot, count) in branch.counts {
                        let maximum = result.counts.entry(slot).or_default();
                        *maximum = (*maximum).max(count);
                    }
                }
            }
            AtomView::Mul(product) => {
                for factor in product.iter() {
                    result.append(Self::observe_indices(
                        factor,
                        slots,
                        cache,
                        scope,
                        excluded_tensor_heads,
                    )?);
                }
            }
            AtomView::Pow(power) => {
                let (base, exponent) = power.get_base_exp();
                // A closed base owns its dummy indices for every power. An
                // open base repeats its explicit ports, so those counts remain
                // visible to the surrounding product. Observe never normalizes
                // leaves or allocates dummy indices.
                for (position, term) in [base, exponent].into_iter().enumerate() {
                    let mut branch =
                        Self::observe_indices(term, slots, cache, scope, excluded_tensor_heads)?;
                    if position == 0
                        && !branch.counts.is_empty()
                        && branch.counts.values().all(|count| *count <= 2)
                    {
                        result.selected_pair |= branch.selected_pair;
                        let interface =
                            match Self::scope_explicit_occurrences(base, &mut branch, slots) {
                                Ok(interface) => Some(interface),
                                Err(error) if scope.validate => return Err(error),
                                Err(_) => None,
                            };
                        if interface
                            .as_ref()
                            .is_some_and(|interface| interface.canonical().is_scalar())
                        {
                            continue;
                        }
                        if scope.validate
                            && let Some(interface) = &interface
                        {
                            Self::power_interface(interface, exponent, &value.to_string())?;
                        }
                        if let Ok(exponent) = Rational::try_from(exponent)
                            && exponent.denominator() == 1
                        {
                            let repetitions = exponent.numerator().abs();
                            let repetitions = usize::try_from(repetitions).unwrap_or(3).min(3);
                            for count in branch.counts.values_mut() {
                                *count = count.saturating_mul(repetitions).min(3);
                            }
                        }
                    }
                    result.selected_pair |= branch.selected_pair;
                    for (slot, count) in branch.counts {
                        let maximum = result.counts.entry(slot).or_default();
                        *maximum = (*maximum).max(count);
                    }
                }
            }
            AtomView::Fun(function) => {
                let symbol = function.get_symbol();
                scope.validate &= !symbol.is_scalar();
                if scope.validate {
                    scope.owner = ChainNestingError::enter(scope.owner, symbol, &SPENSO_TAG)
                        .map_err(|error| TensorInferenceError::invalid(error.to_string()))?;
                }
                let tensor_leaf = symbol.has_tag(&SPENSO_TAG.tensor);
                let count_arguments = scope.count_indices
                    && (tensor_leaf || composition::is_composite_head(symbol, &SPENSO_TAG))
                    && !(composition::is_tensor_leaf_head(symbol)
                        && excluded_tensor_heads.contains(&symbol));
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
                        result.append(Self::observe_indices(
                            argument,
                            slots,
                            cache,
                            child,
                            excluded_tensor_heads,
                        )?);
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
        if let Some(representation) = scope.pair_representation {
            result.selected_pair |= result
                .counts
                .iter()
                .any(|(slot, count)| slot.rep().base().rep == representation.base() && *count >= 2);
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
            assert_eq!(occurrences.counts[&a], expected);
            assert_eq!(occurrences.validate().is_ok(), expected <= 2);
        }
    }

    #[test]
    fn additive_branches_keep_internal_dummy_indices_local() {
        crate::test_support::test_initialize();
        let mu = spenso::mink!(4, 99451);
        let nu = spenso::mink!(4, 99453);
        let p = FunctionBuilder::new(spenso::vector_symbol!("add_scope_p"))
            .add_arg(&mu)
            .finish();
        let q = FunctionBuilder::new(spenso::vector_symbol!("add_scope_q"))
            .add_arg(&mu)
            .finish();
        let r = FunctionBuilder::new(spenso::vector_symbol!("add_scope_r"))
            .add_arg(&nu)
            .finish();
        let metric = FunctionBuilder::new(ETS.metric)
            .add_args([&mu, &nu])
            .finish();
        let x = Atom::var(symbolica::symbol!("add_scope_x"));
        for expression in [
            (&x + p.pow(2)) * &p,
            (&x + &p * &q) * &p,
            (&metric * &p + &r) * &p,
        ] {
            SymbolicTensor::<PartialStructure>::validate_atom(&expression).unwrap();
            let value = SymbolicTensor::infer(expression.clone()).unwrap();
            assert_eq!(value.expression(), &expression);
            let occurrences = InterfaceInference::scoped_occurrences(
                expression.as_view(),
                &mut SlotMatcher::default(),
            )
            .unwrap();
            occurrences.validate().unwrap();
            assert!(
                InterfaceInference::has_explicit_index_pairs(
                    expression.as_view(),
                    Minkowski {}.into(),
                    &[],
                )
                .unwrap()
            );
        }
        // An explicitly supplied cubic term has a different index scope.
        // Construct it directly; do not expand a factored tensor expression.
        let cubic = &x * &p + p.pow(3);
        assert!(SymbolicTensor::infer(cubic).is_err());
        let invalid_branch = (&x + p.pow(3)) * &q;
        assert!(SymbolicTensor::<PartialStructure>::validate_atom(&invalid_branch).is_err());
        assert!(
            InterfaceInference::has_explicit_index_pairs(
                invalid_branch.as_view(),
                Minkowski {}.into(),
                &[],
            )
            .is_err()
        );
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
            InterfaceInference::observe_indices(
                leaf.as_view(),
                &mut slots,
                &mut cache,
                in_factor,
                &[]
            )
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
                },
                &[],
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
                .counts
                .is_empty()
        );
    }
    #[test]
    fn scalar_inverse_scopes_keep_invalid_local_indices_and_placeholders_visible() {
        crate::test_support::test_initialize();
        let rep = ExtendibleReps::MINKOWSKI.new_rep(Dimension::Concrete(4));
        let slot = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99417));
        let leaf = |name| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(slot.to_atom())
                .finish()
        };
        let a = leaf("inverse_scope_A");
        let b = leaf("inverse_scope_B");
        let c = leaf("inverse_scope_C");
        for denominator in [&a * &b, a.pow(2), &a * &b + Atom::one()] {
            let inverse = spenso::bracket!(denominator).pow(-1);
            let value = &a * inverse;
            let occurrences = InterfaceInference::scoped_occurrences(
                value.as_view(),
                &mut SlotMatcher::default(),
            )
            .unwrap();
            assert_eq!(occurrences.counts[&slot], 1);
            occurrences.validate().unwrap();
        }
        let invalid = spenso::bracket!(&a * &b * &c).pow(-1);
        let occurrences =
            InterfaceInference::scoped_occurrences(invalid.as_view(), &mut SlotMatcher::default())
                .unwrap();
        assert!(occurrences.validate().is_err());
        let malformed = spenso::bracket!(&a * &b * Atom::var(SPENSO_TAG.chain_in)).pow(-1);
        assert!(
            InterfaceInference::scoped_occurrences(
                malformed.as_view(),
                &mut SlotMatcher::default()
            )
            .is_err()
        );
    }

    #[test]
    fn powers_share_closed_dummy_scopes_and_open_port_multiplicity() {
        crate::test_support::test_initialize();
        let slot = spenso::mink!(4, 99431);
        let leaf = |name| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(&slot)
                .finish()
        };
        let a = leaf("power_scope_A");
        let b = leaf("power_scope_B");
        let outer = leaf("power_scope_outer");
        let closed = spenso::bracket!(&a * &b);
        for exponent in [-2, -1, 2, 7] {
            let expression = closed.pow(exponent) * &outer;
            SymbolicTensor::<PartialStructure>::validate_atom(&expression).unwrap();
            let value = SymbolicTensor::infer(expression.clone()).unwrap();
            assert_eq!(value.expression(), &expression);
            assert_eq!(value.structure().canonical().external_structure().len(), 1);
            let counts = InterfaceInference::index_occurrences(
                expression.as_view(),
                &mut SlotMatcher::default(),
                &mut Cache::default(),
            );
            assert_eq!(counts.counts.values().copied().collect::<Vec<_>>(), vec![1]);
            assert!(
                InterfaceInference::has_explicit_index_pairs(
                    expression.as_view(),
                    Minkowski {}.into(),
                    &[]
                )
                .unwrap()
            );
        }
        let invalid = a.pow(2) * &outer;
        assert!(SymbolicTensor::<PartialStructure>::validate_atom(&invalid).is_err());
        assert!(SymbolicTensor::infer(invalid.clone()).is_err());
        assert!(
            InterfaceInference::has_explicit_index_pairs(
                invalid.as_view(),
                Minkowski {}.into(),
                &[]
            )
            .is_err()
        );
        assert!(
            InterfaceInference::index_occurrences(
                invalid.as_view(),
                &mut SlotMatcher::default(),
                &mut Cache::default()
            )
            .validate()
            .is_err()
        );
    }

    #[test]
    fn powered_composites_repeat_only_boundary_ports() {
        crate::test_support::test_initialize();
        let mu = spenso::mink!(4, 99441);
        let nu = spenso::mink!(4, 99443);
        let p = FunctionBuilder::new(spenso::vector_symbol!("power_boundary_p"))
            .add_arg(&mu)
            .finish();
        let outer = FunctionBuilder::new(spenso::vector_symbol!("power_boundary_outer"))
            .add_arg(&mu)
            .finish();
        let metric = FunctionBuilder::new(ETS.metric)
            .add_arg(&mu)
            .add_arg(&nu)
            .finish();
        let base = spenso::bracket!(metric * &p);
        let x = Atom::var(symbolica::symbol!("power_boundary_x"));
        let unbracketed = FunctionBuilder::new(ETS.metric)
            .add_arg(&mu)
            .add_arg(&nu)
            .finish()
            * &p;
        for (base, normalized_base) in [
            (base.clone(), base.clone()),
            (&base * &x + &base, &unbracketed * &x + &unbracketed),
        ] {
            let expression = base.pow(2) * &outer;
            SymbolicTensor::<PartialStructure>::validate_atom(&expression).unwrap();
            let value = SymbolicTensor::infer(expression.clone()).unwrap();
            assert_eq!(value.expression(), &(normalized_base.pow(2) * &outer));
            assert_eq!(value.structure().canonical().external_structure().len(), 1);
            assert!(
                InterfaceInference::has_explicit_index_pairs(
                    expression.as_view(),
                    Minkowski {}.into(),
                    &[]
                )
                .unwrap()
            );
            let contracted = value
                .contract(Default::default())
                .unwrap()
                .resolved()
                .unwrap()
                .to_dots()
                .unwrap();
            SymbolicTensor::<PartialStructure>::validate_atom(contracted.expression()).unwrap();
            assert_eq!(contracted.structure(), value.structure());
        }
    }

    #[test]
    fn explicit_pair_query_preserves_factored_branches_and_exclusions() {
        crate::test_support::test_initialize();
        let mu = spenso::mink!(4, 99501);
        let nu = spenso::mink!(4, 99503);
        let leaf = |name, slot: &Atom| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(slot)
                .finish()
        };
        let a = leaf("incidence_A", &mu);
        let b = leaf("incidence_B", &mu);
        let c = leaf("incidence_C", &mu);
        let d = leaf("incidence_D", &nu);
        let x = Atom::var(symbolica::symbol!("incidence_x"));
        let y = Atom::var(symbolica::symbol!("incidence_y"));
        let pairs = |value: &Atom| {
            InterfaceInference::has_explicit_index_pairs(value.as_view(), Minkowski {}.into(), &[])
        };
        assert!(!pairs(&(&a + &b)).unwrap());
        assert!(pairs(&((&a + &b) * &c)).unwrap());
        assert!(!pairs(&((&a + Atom::one()) * &d)).unwrap());
        assert!(pairs(&((&a + Atom::one()) * &b)).unwrap());
        assert!(pairs(&((&a + &b) * &b * &c)).is_err());
        let weighted = (&x - &y) * &a * &b * (&x + &y).pow(30);
        let original = weighted.clone();
        assert!(pairs(&weighted).unwrap());
        assert_eq!(weighted, original);

        let momentum = SPENSO_TAG.tensor_symbol("incidence_momentum");
        let p = FunctionBuilder::new(momentum).add_arg(&mu).finish();
        assert!(
            !InterfaceInference::has_explicit_index_pairs(
                (&a * &p).as_view(),
                Minkowski {}.into(),
                &[momentum],
            )
            .unwrap()
        );
        // Equalized momentum leaves must not cancel two distinct branches.
        let q = FunctionBuilder::new(momentum)
            .add_arg(1)
            .add_arg(&mu)
            .finish();
        let weighted = (&p - &q) * &a * &b;
        assert!(
            InterfaceInference::has_explicit_index_pairs(
                weighted.as_view(),
                Minkowski {}.into(),
                &[momentum],
            )
            .unwrap()
        );
        // The diagnostic query also accepts nested alternatives whose
        // interfaces differ. Failed scope inference must leave their original
        // occurrences available, rather than introduce typed admission here.
        let nested = ((&a + Atom::one()) * &d + &d) * &b;
        assert!(pairs(&nested).unwrap());
        let nested_filtered = ((&a * &p + Atom::one()) * &d + &d) * &b;
        assert!(
            !InterfaceInference::has_explicit_index_pairs(
                nested_filtered.as_view(),
                Minkowski {}.into(),
                &[momentum],
            )
            .unwrap()
        );
        let scalar = FunctionBuilder::new(symbolica::symbol!("incidence_metadata"; Scalar))
            .add_arg(&b * &c)
            .finish();
        assert!(!pairs(&(&a * scalar)).unwrap());
        assert!(
            InterfaceInference::has_explicit_index_pairs(
                spenso::bracket!(&a * &b).as_view(),
                Minkowski {}.into(),
                &[SPENSO_TAG.bracket],
            )
            .unwrap()
        );
    }

    #[test]
    fn explicit_pair_query_keeps_power_scopes_and_rejects_malformed_indices() {
        crate::test_support::test_initialize();
        let mu = spenso::mink!(4, 99511);
        let leaf = |name| {
            FunctionBuilder::new(SPENSO_TAG.tensor_symbol(name))
                .add_arg(&mu)
                .finish()
        };
        let a = leaf("incidence_power_A");
        let b = leaf("incidence_power_B");
        let c = leaf("incidence_power_C");
        let pairs = |value: &Atom| {
            InterfaceInference::has_explicit_index_pairs(value.as_view(), Minkowski {}.into(), &[])
        };
        assert!(pairs(&a.pow(2)).unwrap());
        assert!(pairs(&a.pow(3)).is_err());
        let scalar = spenso::bracket!(&a * &b);
        for power in [-2, -1, 2, 7] {
            assert!(pairs(&(scalar.pow(power) * &c)).unwrap());
        }
        assert!(pairs(&spenso::bracket!(&a * &b * &c).pow(-1)).is_err());
        assert!(pairs(&Atom::var(SPENSO_TAG.chain_in)).is_err());
        let momentum = SPENSO_TAG.tensor_symbol("incidence_power_momentum");
        let p = FunctionBuilder::new(momentum).add_arg(&mu).finish();
        for power in [-1, 2, 7] {
            let value = spenso::bracket!(&a * &p).pow(power) * &b;
            assert!(
                !InterfaceInference::has_explicit_index_pairs(
                    value.as_view(),
                    Minkowski {}.into(),
                    &[momentum],
                )
                .unwrap()
            );
        }
    }
}
