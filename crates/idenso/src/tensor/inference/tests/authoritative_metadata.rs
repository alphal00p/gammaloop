use super::*;
use crate::gamma;
use crate::tensor::{ContractSettings, SymbolicNetExt, SymbolicNetParse};
use crate::test_support;
use spenso::network::parsing::ParseSettings;
use spenso::shadowing::TensorCollectFilter;
use spenso::{euc, mink, tensor_symbol, trace, vector_symbol};
use symbolica::{function, symbol};

#[test]
fn materializers_reuse_exact_publication_and_release_replaced_roots() {
    test_support::test_initialize();
    let leaf = function!(
        tensor_symbol!("authoritative_metadata::prepared_leaf"),
        mink!(4, 99721)
    );
    let coefficient = Atom::var(symbol!("authoritative_metadata::old_coefficient"));
    let source = SymbolicTensor::infer(&leaf * &coefficient).unwrap();
    let previous = Arc::downgrade(source.proofs.admitted_interfaces.get().unwrap());
    let unchanged = source
        .borrowed_parts()
        .with_generated_expression("unchanged materializer", source.expression.clone(), &[])
        .unwrap();
    assert!(Arc::ptr_eq(
        unchanged.proofs.admitted_interfaces.get().unwrap(),
        source.proofs.admitted_interfaces.get().unwrap(),
    ));
    let expression = Atom::num(2) * &leaf;
    let mut proofs = source.proofs.clone();
    proofs.publish_interfaces(expression.as_view(), source.structure.clone(), []);
    let prepared = Arc::clone(proofs.admitted_interfaces.get().unwrap());
    let leaves = proofs.parser_leaf_interfaces();
    let result = source
        .with_published_materialization("changed materializer", expression, proofs)
        .unwrap();
    assert!(Arc::ptr_eq(
        result.proofs.admitted_interfaces.get().unwrap(),
        &prepared
    ));
    assert!(Arc::ptr_eq(
        result.proofs.leaf_interfaces.get().unwrap(),
        &leaves
    ));
    assert!(
        result
            .proofs
            .retained_interface(coefficient.as_view())
            .is_none()
    );
    drop(unchanged);
    drop(source);
    assert!(previous.upgrade().is_none());
}

#[test]
fn unchanged_materializer_keeps_the_declared_order_over_a_prepared_storage_order() {
    test_support::test_initialize();
    let expression = function!(
        tensor_symbol!("authoritative_metadata::ordered_materializer"),
        mink!(4, 99727),
        mink!(4, 99729)
    );
    let source = SymbolicTensor::infer(expression.clone()).unwrap();
    let reversed =
        PartialStructure::from_logical_slots(source.structure.logical_slots().into_iter().rev());
    assert_ne!(reversed, source.structure);
    let mut proofs = source.proofs.clone();
    proofs.publish_interfaces(expression.as_view(), reversed, []);
    let result = source
        .with_published_materialization("prepared order", expression, proofs)
        .unwrap();
    assert_eq!(result.structure, source.structure);
    assert_eq!(
        result
            .proofs
            .retained_interface(result.expression.as_view()),
        Some(result.structure.clone()),
    );
}

#[test]
fn a_prepared_materializer_cannot_certify_its_own_callback_output() {
    test_support::test_initialize();
    let target = tensor_symbol!("authoritative_metadata::callback_result");
    let head = tensor_symbol!(
        "authoritative_metadata::materializer_callback",
        norm = move |node, out| {
            let AtomView::Fun(function) = node else {
                return;
            };
            let mut args = function.iter();
            let Some(selector) = args.next() else { return };
            if selector == Atom::one().as_view() {
                **out = Atom::one();
            } else if selector == Atom::num(2).as_view()
                && let (Some(first), Some(second)) = (args.next(), args.next())
            {
                **out = function!(target, second, first);
            }
        }
    );
    let first = mink!(4, 99739);
    let second = mink!(4, 99741);
    let source = SymbolicTensor::infer(function!(head, 0, &first, &second)).unwrap();
    for selector in [1, 2] {
        let expression = function!(head, selector, &first, &second);
        let actual = SymbolicTensor::infer(expression.clone()).unwrap();
        let mut prepared = source.proofs.clone();
        prepared.publish_interfaces(
            expression.as_view(),
            source.structure.clone(),
            [&actual.proofs],
        );
        let result = source.with_published_materialization("callback", expression, prepared);
        if selector == 1 {
            assert!(actual.is_scalar());
            assert!(
                result
                    .unwrap_err()
                    .to_string()
                    .contains("compatible tensor interface")
            );
        } else {
            assert_ne!(source.structure, actual.structure);
            let value = result.unwrap();
            assert_eq!(value.structure, source.structure);
            assert_eq!(
                value.proofs.retained_interface(value.expression.as_view()),
                Some(source.structure.clone())
            );
            assert_eq!(
                value.proofs.parser_leaf_interfaces().lock().unwrap()[&value.expression].1,
                [second.clone(), first.clone()]
            );
        }
    }
}

#[test]
fn chainify_rejects_a_callback_erasing_a_spectator_from_its_generated_composite() {
    test_support::test_initialize();
    let calls = Arc::new(AtomicUsize::new(0));
    let count = Arc::clone(&calls);
    let template = tensor_symbol!("authoritative_metadata::changed_chain_template");
    let head = tensor_symbol!(
        "authoritative_metadata::chainify_callback",
        norm = move |node, out| {
            if let AtomView::Fun(function) = node
                && function.iter().any(|arg| matches!(arg, AtomView::Var(var) if var.get_symbol() == SPENSO_TAG.chain_in))
            {
                count.fetch_add(1, Ordering::Relaxed);
                **out = function!(template, SPENSO_TAG.chain_in, SPENSO_TAG.chain_out);
            }
        }
    );
    let spectator = LibraryRep::from(ColorFundamental {})
        .new_rep(3)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99747))
        .to_atom();
    let source =
        SymbolicTensor::infer(function!(head, mink!(4, 99743), mink!(4, 99745), spectator))
            .unwrap();
    let initial_leaves = source
        .proofs
        .parser_leaf_interfaces()
        .lock()
        .unwrap()
        .clone();
    assert_eq!(source.rank(), 3);
    assert!(source.chainified(LibraryRep::from(Minkowski {})).is_err());
    assert!(calls.load(Ordering::Relaxed) > 0);
    assert_eq!(
        *source.proofs.parser_leaf_interfaces().lock().unwrap(),
        initial_leaves
    );
}

#[test]
fn independent_validation_checks_the_same_ports_without_collecting_admission_facts() {
    test_support::test_initialize();
    let calls = Arc::new(AtomicUsize::new(0));
    let counter = Arc::clone(&calls);
    let head = tensor_symbol!(
        "authoritative_metadata::validation_only",
        norm = move |_, _| {
            counter.fetch_add(1, Ordering::Relaxed);
        }
    );
    let leaf = function!(head, mink!(4, 99719));
    let expression = Atom::num(3) * &leaf;
    calls.store(0, Ordering::Relaxed);
    let observe = |purpose| {
        SymbolicTensor::observed_interface(&expression, LeafInference::Observe, purpose).unwrap()
    };
    let checked = observe(ObservationPurpose::Validate);
    let admitted = observe(ObservationPurpose::Publish);
    assert_eq!(checked.interface, admitted.interface);
    assert_eq!(
        checked
            .occurrences
            .unwrap()
            .oriented_entries()
            .collect::<Vec<_>>(),
        admitted
            .occurrences
            .unwrap()
            .oriented_entries()
            .collect::<Vec<_>>(),
    );
    assert!(checked.leaf_interfaces.lock().unwrap().is_empty());
    assert!(admitted.leaf_interfaces.lock().unwrap().contains_key(&leaf));
    let retained = |interfaces| TensorProofs {
        admitted_interfaces: OnceLock::from(Arc::new(interfaces)),
        ..TensorProofs::default()
    };
    assert!(!retained(checked.interfaces).owns_expression(expression.as_view()));
    assert_eq!(
        retained(admitted.interfaces).retained_interface(leaf.as_view()),
        Some(admitted.interface),
    );
    assert_eq!(calls.load(Ordering::Relaxed), 0);
}

#[test]
fn publishing_an_admitted_parser_leaf_keeps_callback_ports_scoped() {
    test_support::test_initialize();
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&calls);
    let head = tensor_symbol!(
        "authoritative_metadata::parser_callback_leaf",
        norm = move |_, _| {
            observed.fetch_add(1, Ordering::Relaxed);
        }
    );
    let ports = [mink!(4, 99417), mink!(4, 99415)];
    let leaf = function!(head, &ports[0], &ports[1]);
    assert!(calls.load(Ordering::Relaxed) > 0);
    let sibling = function!(
        tensor_symbol!("authoritative_metadata::parser_callback_sibling"),
        mink!(4, 99419)
    );
    let source = SymbolicTensor::infer(&leaf * &sibling).unwrap();
    calls.store(0, Ordering::Relaxed);
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));

    let network = source
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    let parsed = network
        .store
        .tensors
        .iter()
        .find(|tensor| tensor.expression == leaf)
        .unwrap();
    assert!(parsed.proofs.owns_expression(leaf.as_view()));
    assert_eq!(
        parsed
            .proofs
            .retained_interface(leaf.as_view())
            .unwrap()
            .logical_slots()
            .len(),
        ports.len()
    );
    assert!(
        parsed
            .proofs
            .retained_interface(sibling.as_view())
            .is_none()
    );
    let leaves = parsed.proofs.leaf_interfaces.get().unwrap().lock().unwrap();
    assert_eq!(leaves[&leaf].1, ports);
    assert!(!leaves.contains_key(&sibling));
    assert_eq!(calls.load(Ordering::Relaxed), 0);
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
    assert!(source.with_checked_expression(Atom::one()).is_err());
}

#[test]
fn renamed_generic_leaves_survive_collection_and_network_execution_without_inference() {
    test_support::test_initialize();
    let rep = LibraryRep::from(Minkowski {}).new_rep(4);
    let old = AbstractIndex::Normal(99631);
    let new = AbstractIndex::Normal(99633);
    let left = tensor_symbol!("authoritative_metadata::pipeline_left");
    let right = tensor_symbol!("authoritative_metadata::pipeline_right");
    let original = function!(left, rep.slot::<AbstractIndex, _>(old).to_atom())
        * function!(right, rep.slot::<AbstractIndex, _>(old).to_atom());
    let expected = function!(left, rep.slot::<AbstractIndex, _>(new).to_atom())
        * function!(right, rep.slot::<AbstractIndex, _>(new).to_atom());
    let source = SymbolicTensor::infer(original).unwrap();
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));

    let renamed = source
        .renamed_indices(&HashMap::from([(old, new)]))
        .unwrap();
    assert_eq!(renamed.expression(), &expected);
    let rows = renamed
        .coefficient_list(TensorCollectFilter::Reps([rep.rep]))
        .unwrap();
    assert!(!rows.is_empty());
    let mut reconstructed = Atom::Zero;
    for (tensor, coefficient) in rows {
        let value = tensor.multiply(&coefficient).unwrap();
        let network = value
            .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
            .unwrap();
        reconstructed += network.simple_execute::<()>().unwrap();
    }
    assert_eq!(reconstructed, expected);
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

#[test]
fn expansion_keeps_coalesced_oriented_tensor_copies() {
    test_support::test_initialize();
    let slot = LibraryRep::from(ColorFundamental {})
        .new_rep(3)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99489));
    let vector = SymbolicTensor::infer(function!(
        tensor_symbol!("authoritative_metadata::expanded_copy"),
        slot.to_atom()
    ))
    .unwrap();
    let x = Atom::var(symbol!("authoritative_metadata::expanded_x"));
    let y = Atom::var(symbol!("authoritative_metadata::expanded_y"));
    let coefficient = SymbolicTensor::infer(&x + &y).unwrap();
    let weighted = vector
        .multiply(&vector)
        .unwrap()
        .multiply(&coefficient)
        .unwrap();
    assert_eq!(weighted.rank(), 2);
    let expected = weighted.expression.expand();
    OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
    SCOPE_VALIDATIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));
    let expanded = weighted.expanded(None, false).unwrap();
    assert_eq!(expanded.expression, expected);
    assert_eq!(expanded.structure, weighted.structure);
    assert_eq!(
        expanded.explicit_index_occurrences().unwrap().count(slot),
        2
    );
    assert!(matches!(
        expanded.multiply(&vector),
        Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
    ));
    assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
    assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

#[test]
fn zero_algebra_maps_hidden_pair_before_later_explicit_collisions() {
    test_support::test_initialize();
    let index = AbstractIndex::Normal(99471);
    let minkowski = LibraryRep::from(Minkowski {});
    let physical = minkowski.new_rep(4).slot::<AbstractIndex, _>(index);
    let extended = minkowski.new_rep(6).slot::<AbstractIndex, _>(index);
    // Ordinary tensor heads keep the admitted pair in product notation.
    let source = SymbolicTensor::infer(
        function!(
            tensor_symbol!("authoritative_metadata::hidden_a"),
            physical.to_atom()
        ) * function!(
            tensor_symbol!("authoritative_metadata::hidden_b"),
            physical.to_atom()
        ),
    )
    .unwrap();
    let physical_leaf = SymbolicTensor::infer(function!(
        tensor_symbol!("authoritative_metadata::physical_leaf"),
        physical.to_atom()
    ))
    .unwrap();
    let extended_leaf = SymbolicTensor::infer(function!(
        tensor_symbol!("authoritative_metadata::extended_leaf"),
        extended.to_atom()
    ))
    .unwrap();
    assert!(source.is_scalar());
    assert_eq!(
        source.explicit_index_occurrences().unwrap().count(physical),
        2
    );
    physical_leaf.explicit_index_occurrences().unwrap();
    extended_leaf.explicit_index_occurrences().unwrap();
    OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
    SCOPE_VALIDATIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));

    let zero = source.with_algebra_result(Atom::Zero).unwrap();
    assert!(zero.is_zero());
    assert!(zero.is_scalar());
    assert_eq!(
        zero.explicit_index_occurrences().unwrap().count(physical),
        2
    );
    assert!(matches!(
        zero.multiply(&physical_leaf),
        Err(TensorCompositionError::InvalidExplicitMultiplicity { index: found, occurrences: 3 })
            if found == index
    ));

    let mapped = zero.with_lorentz_dimension(Dimension::Concrete(6)).unwrap();
    assert!(mapped.is_zero());
    assert!(mapped.is_scalar());
    assert_eq!(
        mapped.explicit_index_occurrences().unwrap().count(extended),
        2
    );
    assert_eq!(
        mapped.explicit_index_occurrences().unwrap().count(physical),
        0
    );
    assert!(matches!(
        mapped.multiply(&extended_leaf),
        Err(TensorCompositionError::InvalidExplicitMultiplicity { index: found, occurrences: 3 })
            if found == index
    ));
    let distinct_dimension = mapped.multiply(&physical_leaf).unwrap();
    assert!(distinct_dimension.is_zero());
    assert_eq!(
        distinct_dimension.structure.logical_slots(),
        [physical.rep().slot(PartialIndex::Explicit(index))]
    );
    let occurrences = distinct_dimension.explicit_index_occurrences().unwrap();
    assert_eq!(occurrences.count(extended), 2);
    assert_eq!(occurrences.count(physical), 1);
    assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
    assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

#[test]
fn dual_free_dimension_mapping_preserves_orientation_and_later_collision_history() {
    test_support::test_initialize();
    let index = AbstractIndex::Normal(99473);
    let dual_three = LibraryRep::from(ColorAntiFundamental {}).new_rep(3);
    let dual_four = LibraryRep::from(ColorAntiFundamental {}).new_rep(4);
    let source_slot = dual_three.slot::<AbstractIndex, _>(index);
    let mapped_slot = dual_four.slot::<AbstractIndex, _>(index);
    let partner_slot = dual_four.dual().slot::<AbstractIndex, _>(index);
    let head = tensor_symbol!("authoritative_metadata::dual_free");
    let source = SymbolicTensor::infer(function!(head, source_slot.to_atom())).unwrap();
    let partner = SymbolicTensor::infer(function!(
        tensor_symbol!("authoritative_metadata::base_partner"),
        partner_slot.to_atom()
    ))
    .unwrap();
    assert_eq!(
        source
            .explicit_index_occurrences()
            .unwrap()
            .count(source_slot),
        1
    );
    partner.explicit_index_occurrences().unwrap();
    OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
    SCOPE_VALIDATIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));

    let mut calls = Vec::new();
    let mapped = source
        .map_representation_dimensions(|representation, dimension| {
            calls.push((representation, dimension.to_owned()));
            if representation == dual_three.rep && dimension == Atom::num(3).as_view() {
                Atom::num(4)
            } else {
                dimension.to_owned()
            }
        })
        .unwrap();
    // The base orientation is absent from this source. Mapping its singleton
    // must neither introduce that callback nor retain dimension three.
    assert_eq!(calls, [(dual_three.rep, Atom::num(3))]);
    assert_eq!(mapped.expression(), &function!(head, mapped_slot.to_atom()));
    assert_eq!(
        mapped.structure.logical_slots(),
        [dual_four.slot(PartialIndex::Explicit(index))]
    );
    let occurrences = mapped.explicit_index_occurrences().unwrap();
    assert_eq!(occurrences.count(source_slot), 0);
    assert_eq!(occurrences.count(mapped_slot), 1);
    assert_eq!(occurrences.count(partner_slot), 1);
    assert_eq!(
        occurrences.entries().collect::<Vec<_>>(),
        [(mapped_slot, 1)]
    );

    let pair = mapped.multiply(&partner).unwrap();
    assert!(pair.is_scalar());
    assert_eq!(
        pair.explicit_index_occurrences()
            .unwrap()
            .count(mapped_slot),
        2
    );
    assert!(matches!(
        pair.multiply(&mapped),
        Err(TensorCompositionError::InvalidExplicitMultiplicity { index: found, occurrences: 3 })
            if found == index
    ));
    assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
    assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

#[test]
fn admitted_arithmetic_owns_canonical_named_ports_and_composes_without_observation() {
    test_support::test_initialize();
    let minkowski = LibraryRep::from(Minkowski {}).new_rep(4);
    let indices = [99479, 99477].map(AbstractIndex::Normal);
    let slots = indices.map(|index| minkowski.slot(index));
    let ports = indices.map(|index| minkowski.slot(PartialIndex::Explicit(index)));
    // Include a callback-sensitive leaf and a logical order different from the
    // other operand's order. Admission owns both before the counters are reset.
    let a = SymbolicTensor::from_declared_parts(
        FunctionBuilder::new(tensor_symbol!(
            "authoritative_metadata::arithmetic_a",
            norm = |_, _| {}
        ))
        .add_args(ports.map(composition::port_atom))
        .finish(),
        PartialStructure::from_logical_slots(ports),
    )
    .unwrap();
    let b = SymbolicTensor::from_declared_parts(
        FunctionBuilder::new(tensor_symbol!("authoritative_metadata::arithmetic_b"))
            .add_args(ports.into_iter().rev().map(composition::port_atom))
            .finish(),
        PartialStructure::from_logical_slots(ports.into_iter().rev()),
    )
    .unwrap();
    let scalar =
        SymbolicTensor::infer(Atom::var(symbol!("authoritative_metadata::weight"))).unwrap();
    let scalar_zero = SymbolicTensor::infer(Atom::Zero).unwrap();
    let probe = SymbolicTensor::infer(function!(
        tensor_symbol!("authoritative_metadata::arithmetic_probe"),
        slots[0].to_atom()
    ))
    .unwrap();
    for source in [&a, &b, &scalar, &scalar_zero, &probe] {
        source.explicit_index_occurrences().unwrap();
    }
    OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
    SCOPE_VALIDATIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));

    let sum = a.try_add(&b).unwrap();
    let difference = a.try_sub(&b).unwrap();
    let negated = -sum.clone();
    let weighted = sum.multiply(&scalar).unwrap();
    let cancelled = a.try_sub(&a).unwrap();
    let scaled_zero = sum.multiply(&scalar_zero).unwrap();
    assert!(cancelled.is_zero());
    assert!(scaled_zero.is_zero());
    for result in [sum, difference, negated, weighted, cancelled, scaled_zero] {
        assert_eq!(result.structure.logical_slots(), [ports[1], ports[0]]);
        let occurrences = result.explicit_index_occurrences().unwrap();
        assert_eq!(occurrences.entries().count(), 2);
        for slot in slots {
            assert_eq!(occurrences.count(slot), 1);
        }

        // The retained public port must contract once, including for a zero,
        // and that pair must reject a third use in a later composition.
        let contracted = result.multiply(&probe).unwrap();
        assert_eq!(contracted.structure.logical_slots(), [ports[1]]);
        assert_eq!(
            contracted
                .explicit_index_occurrences()
                .unwrap()
                .count(slots[0]),
            2
        );
        assert!(matches!(
            contracted.multiply(&probe),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { index: found, occurrences: 3 })
                if found == indices[0]
        ));
    }
    assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
    assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

#[test]
fn additive_branches_keep_maximum_public_multiplicity_in_both_operand_orders() {
    test_support::test_initialize();
    let index = AbstractIndex::Normal(99481);
    let slot = LibraryRep::from(Minkowski {})
        .new_rep(4)
        .slot::<AbstractIndex, _>(index);
    // Explicit declared intake retains the slot-use history; tensor-free inference is opaque.
    let scalar = SymbolicTensor::from_declared_parts(
        slot.to_atom(),
        PartialStructure::from_logical_slots([]),
    )
    .unwrap();
    let vector = SymbolicTensor::infer(function!(
        tensor_symbol!("authoritative_metadata::branch_vector"),
        slot.to_atom()
    ))
    .unwrap();
    let extra = SymbolicTensor::infer(function!(
        tensor_symbol!("authoritative_metadata::branch_extra"),
        slot.to_atom()
    ))
    .unwrap();
    assert!(scalar.is_scalar());
    assert_eq!(scalar.explicit_index_occurrences().unwrap().count(slot), 1);
    assert_eq!(vector.explicit_index_occurrences().unwrap().count(slot), 1);
    extra.explicit_index_occurrences().unwrap();
    OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
    SCOPE_VALIDATIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));

    let branch = scalar.multiply(&vector).unwrap();
    let interface = [slot.rep().slot(PartialIndex::Explicit(index))];
    assert_eq!(branch.structure.logical_slots(), interface);
    assert_eq!(branch.explicit_index_occurrences().unwrap().count(slot), 2);
    for result in [
        vector.try_add(&branch).unwrap(),
        branch.try_add(&vector).unwrap(),
        vector.try_sub(&branch).unwrap(),
        branch.try_sub(&vector).unwrap(),
    ] {
        assert_eq!(result.structure.logical_slots(), interface);
        assert_eq!(
            result
                .explicit_index_occurrences()
                .unwrap()
                .entries()
                .collect::<Vec<_>>(),
            [(slot, 2)]
        );
        assert!(matches!(
            result.multiply(&extra),
            Err(TensorCompositionError::InvalidExplicitMultiplicity { index: found, occurrences: 3 })
                if found == index
        ));
    }
    assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
    assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

#[test]
fn mapping_hidden_base_use_separates_it_from_the_free_dual_port() {
    test_support::test_initialize();
    let index = AbstractIndex::Normal(99483);
    let base_three = LibraryRep::from(ColorFundamental {}).new_rep(3);
    let base_four = LibraryRep::from(ColorFundamental {}).new_rep(4);
    let base_slot = base_three.slot::<AbstractIndex, _>(index);
    let mapped_base_slot = base_four.slot::<AbstractIndex, _>(index);
    let dual_slot = base_three.dual().slot::<AbstractIndex, _>(index);
    let scalar_base = SymbolicTensor::from_declared_parts(
        base_slot.to_atom(),
        PartialStructure::from_logical_slots([]),
    )
    .unwrap();
    let head = tensor_symbol!("authoritative_metadata::hidden_base_dual");
    let free_dual = SymbolicTensor::infer(function!(head, dual_slot.to_atom())).unwrap();
    let compatible = SymbolicTensor::infer(function!(
        tensor_symbol!("authoritative_metadata::hidden_base_compatible"),
        base_slot.to_atom()
    ))
    .unwrap();
    // Declared scalar intake retains the hidden base use; opaque inference does not.
    assert!(scalar_base.is_scalar());
    assert_eq!(
        scalar_base
            .explicit_index_occurrences()
            .unwrap()
            .count(base_slot),
        1
    );
    free_dual.explicit_index_occurrences().unwrap();
    compatible.explicit_index_occurrences().unwrap();
    OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
    SCOPE_VALIDATIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));

    let source = free_dual.multiply(&scalar_base).unwrap();
    let interface = [dual_slot.rep().slot(PartialIndex::Explicit(index))];
    assert_eq!(source.structure.logical_slots(), interface);
    assert_eq!(
        source
            .explicit_index_occurrences()
            .unwrap()
            .oriented_entries()
            .collect::<Vec<_>>(),
        [(dual_slot, 2, 1, 1)]
    );
    let mapped = source
        .map_representation_dimensions(|representation, dimension| {
            if representation == base_three.rep && dimension == Atom::num(3).as_view() {
                Atom::num(4)
            } else {
                dimension.to_owned()
            }
        })
        .unwrap();
    assert_eq!(mapped.structure.logical_slots(), interface);
    assert_eq!(
        mapped.expression(),
        &(function!(head, dual_slot.to_atom()) * mapped_base_slot.to_atom())
    );
    let occurrences = mapped.explicit_index_occurrences().unwrap();
    assert_eq!(occurrences.count(base_slot), 1);
    assert_eq!(occurrences.count(mapped_base_slot), 1);
    assert_eq!(
        occurrences.oriented_entries().collect::<Vec<_>>(),
        [(dual_slot, 1, 0, 1), (mapped_base_slot, 1, 1, 0)]
    );

    // The old dimension now has only the free dual use, so a new base use
    // contracts it instead of encountering a spurious third occurrence.
    let contracted = mapped.multiply(&compatible).unwrap();
    assert!(contracted.is_scalar());
    let occurrences = contracted.explicit_index_occurrences().unwrap();
    assert_eq!(occurrences.count(base_slot), 2);
    assert_eq!(occurrences.count(mapped_base_slot), 1);
    assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
    assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);

    // Protect the exact external variance, rather than hiding the opposite
    // orientation's use merely because it shares this collision key.
    let scope = symbol!("authoritative_metadata::hidden_base_scope");
    let scoped_slot = base_three.slot::<AbstractIndex, _>(index.scoped(scope));
    let scoped = source.wrap_indices(scope, true).unwrap();
    assert_eq!(scoped.structure.logical_slots(), interface);
    assert_eq!(
        scoped
            .explicit_index_occurrences()
            .unwrap()
            .count(base_slot),
        1
    );
    assert_eq!(
        scoped
            .explicit_index_occurrences()
            .unwrap()
            .count(scoped_slot),
        1
    );
    assert_eq!(
        scoped.expression(),
        &(function!(head, dual_slot.to_atom()) * scoped_slot.to_atom())
    );
    let probe = SymbolicTensor::infer(function!(
        tensor_symbol!("authoritative_metadata::scoped_base_probe"),
        scoped_slot.to_atom()
    ))
    .unwrap();
    let twice = scoped.multiply(&probe).unwrap();
    assert!(matches!(twice.multiply(&probe),
        Err(TensorCompositionError::InvalidExplicitMultiplicity { index: found, occurrences: 3 })
            if found == index.scoped(scope)
    ));
}

#[test]
fn dummy_scoping_moves_a_pair_closed_inside_a_sum_together() {
    test_support::test_initialize();
    let scope = symbol!("authoritative_metadata::closed_sum_scope");
    let rep = LibraryRep::from(ColorFundamental {}).new_rep(3);
    let index = AbstractIndex::Normal(99493);
    let scoped_index = index.scoped(scope);
    let base = rep.slot::<AbstractIndex, _>(index);
    let dual = rep.dual().slot::<AbstractIndex, _>(index);
    let spectator = rep.slot::<AbstractIndex, _>(AbstractIndex::Normal(99494));
    let leaf = |name: &str, slot: Slot<LibraryRep, AbstractIndex>| {
        function!(
            SPENSO_TAG.tensor_symbol(&format!("authoritative_metadata::closed_sum_{name}")),
            slot.to_atom()
        )
    };
    let x = Atom::var(symbol!("authoritative_metadata::closed_sum_x"));
    let y = Atom::var(symbol!("authoritative_metadata::closed_sum_y"));
    // `a(ī)·b(i)` contracts inside its sum, whose free neighbour `v` reuses the name.
    let closed = |index: AbstractIndex, free: Slot<LibraryRep, AbstractIndex>| {
        &y + leaf("a", free.rep().dual().slot(index)) * leaf("b", free.rep().slot(index))
    };
    let cases = [
        (
            (&x + closed(index, base)) * leaf("v", base),
            (&x + closed(scoped_index, base)) * leaf("v", base),
            vec![base],
        ),
        (
            (&x + closed(index, dual)) * leaf("v", dual),
            (&x + closed(scoped_index, dual)) * leaf("v", dual),
            vec![dual],
        ),
        (
            (&x * leaf("c", spectator) + closed(index, base) * leaf("d", spectator))
                * leaf("v", base),
            (&x * leaf("c", spectator) + closed(scoped_index, base) * leaf("d", spectator))
                * leaf("v", base),
            vec![base, spectator],
        ),
    ];
    for (source, expected, ports) in cases {
        let source = SymbolicTensor::infer(source).unwrap();
        assert!(!source.is_scalar());
        assert_eq!(source.rank(), ports.len());
        let scoped = source.wrap_indices(scope, true).unwrap();
        // Both members of the closed pair are scoped; the free ports are not.
        assert_eq!(scoped.expression(), &expected);
        assert_eq!(scoped.structure, source.structure);
        assert!(!scoped.is_scalar());
        assert_eq!(scoped.rank(), ports.len());
        let interface = scoped
            .structure
            .logical_slots()
            .into_iter()
            .map(|slot| match slot.aind {
                PartialIndex::Explicit(index) => slot.rep().slot::<AbstractIndex, _>(index),
                PartialIndex::Open(_) => panic!("closed-sum fixture has only explicit ports"),
            })
            .collect::<HashSet<_>>();
        assert_eq!(interface, ports.into_iter().collect::<HashSet<_>>());
        assert_eq!(
            SymbolicTensor::infer(expected).unwrap().structure,
            scoped.structure
        );
    }
    // The scoped sum on its own remains a scalar.
    assert!(
        SymbolicTensor::infer(closed(scoped_index, base))
            .unwrap()
            .is_scalar()
    );
}

#[test]
fn dummy_scoping_moves_a_pair_closed_inside_a_scalar_power_together() {
    test_support::test_initialize();
    let scope = symbol!("authoritative_metadata::closed_power_scope");
    let rep = LibraryRep::from(ColorFundamental {}).new_rep(3);
    let index = AbstractIndex::Normal(99495);
    let scoped_index = index.scoped(scope);
    let leaf = |name: &str, slot: Slot<LibraryRep, AbstractIndex>| {
        SymbolicTensor::infer(function!(
            SPENSO_TAG.tensor_symbol(&format!("authoritative_metadata::closed_power_{name}")),
            slot.to_atom()
        ))
        .unwrap()
    };
    let closed_scope = |value: &SymbolicTensor<PartialStructure>| {
        SymbolicTensor::<PartialStructure>::is_closed_scalar_power(
            value.expression().as_view(),
            &mut SlotMatcher::default(),
        )
    };
    // `bracket(a(ī)·b(i))^n` closes its pair, while `v` carries a free port of that name.
    let power = |index: AbstractIndex, free: Representation<LibraryRep>, exponent: i64| {
        let pair = leaf("a", free.dual().slot(index))
            .multiply(&leaf("b", free.slot(index)))
            .unwrap();
        assert!(pair.is_scalar());
        pair.pow(&SymbolicTensor::infer(Atom::num(exponent)).unwrap())
            .unwrap()
    };
    for free in [rep, rep.dual()] {
        let v = leaf("v", free.slot(index));
        let interface =
            [free.slot::<PartialIndex<AbstractIndex>, _>(PartialIndex::Explicit(index))];
        for exponent in [2, 3] {
            let closed = power(index, free, exponent);
            let scoped_closed = power(scoped_index, free, exponent);
            for scalar in [&closed, &scoped_closed] {
                assert!(scalar.is_scalar());
                assert_eq!(scalar.rank(), 0);
                let AtomView::Pow(power) = scalar.expression().as_view() else {
                    panic!("closed pair power is not a power: {}", scalar.expression());
                };
                assert!(matches!(
                    power.get_base_exp().0,
                    AtomView::Fun(base) if base.get_symbol() == SPENSO_TAG.bracket
                ));
                assert!(closed_scope(scalar));
            }
            let source = closed.multiply(&v).unwrap();
            let expected = scoped_closed.multiply(&v).unwrap();
            assert!(!source.is_scalar());
            assert_eq!(source.rank(), 1);
            assert_eq!(source.structure.logical_slots(), interface);
            let scoped = source.wrap_indices(scope, true).unwrap();
            // Both members of the closed pair are scoped; the free port is not.
            assert_eq!(scoped.expression(), expected.expression());
            assert_eq!(scoped.structure, source.structure);
            assert!(!scoped.is_scalar());
            assert_eq!(scoped.rank(), 1);
            assert_eq!(scoped.structure.logical_slots(), interface);
            assert_eq!(
                SymbolicTensor::infer(expected.expression().clone())
                    .unwrap()
                    .structure,
                scoped.structure
            );
        }
    }

    // A tensor-valued power base repeats its ports into the surrounding
    // product, so it is no closed scope and keeps the exact-variance rule.
    let mu = LibraryRep::from(Minkowski {})
        .new_rep(4)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99496));
    let tensor_power = leaf("t", mu)
        .pow(&SymbolicTensor::infer(Atom::num(2)).unwrap())
        .unwrap();
    assert!(matches!(
        tensor_power.expression().as_view(),
        AtomView::Pow(_)
    ));
    assert!(!closed_scope(&tensor_power));
    let free = rep.slot::<AbstractIndex, _>(index);
    let v = leaf("v", free);
    let source = tensor_power.multiply(&v).unwrap();
    assert!(!source.is_scalar());
    assert_eq!(source.rank(), 1);
    let scoped = source.wrap_indices(scope, true).unwrap();
    let exact_variance =
        AbstractIndex::wrap_expression(source.expression().as_view(), scope, |slot| *slot != free);
    assert_eq!(scoped.expression(), &exact_variance);
    assert_eq!(
        scoped.expression(),
        &(leaf("t", mu.rep().slot(mu.aind().scoped(scope)))
            .expression()
            .pow(Atom::num(2))
            * v.expression())
    );
    assert_eq!(scoped.structure, source.structure);
    assert_eq!(scoped.rank(), 1);
    assert_eq!(
        scoped.structure.logical_slots(),
        [free
            .rep()
            .slot::<PartialIndex<AbstractIndex>, _>(PartialIndex::Explicit(index))]
    );
}

#[test]
fn intake_installs_hidden_history_before_zero_product_finishing() {
    test_support::test_initialize();
    let slot = LibraryRep::from(Minkowski {})
        .new_rep(4)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99591));
    let left = function!(tensor_symbol!("publication::left"), slot.to_atom());
    let right = function!(tensor_symbol!("publication::right"), slot.to_atom());
    let expression = FunctionBuilder::new(SPENSO_TAG.bracket)
        .add_arg(Atom::Zero)
        .add_arg(left)
        .add_arg(right)
        .finish();
    let value = SymbolicTensor::new(expression, PartialStructure::from_logical_slots([])).unwrap();
    assert!(value.is_zero());
    assert_eq!(value.explicit_index_occurrences().unwrap().count(slot), 2);
    let probe = SymbolicTensor::infer(function!(
        tensor_symbol!("publication::probe"),
        slot.to_atom()
    ))
    .unwrap();
    assert!(matches!(
        value.multiply(&probe),
        Err(TensorCompositionError::InvalidExplicitMultiplicity { occurrences: 3, .. })
    ));
}

#[test]
fn publication_borrows_supplied_history_without_observing_the_payload() {
    test_support::test_initialize();
    let slot = LibraryRep::from(Minkowski {})
        .new_rep(4)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99593));
    let mut history = ExplicitIndexOccurrences::default();
    history.add(slot, 2);
    let history = Arc::new(history);
    let value = SymbolicTensor::from_validated_parts(
        Atom::Zero,
        PartialStructure::from_logical_slots([]),
        Arc::clone(&history),
    );
    let scalar = SymbolicTensor::infer(Atom::num(3)).unwrap();
    OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
    SCOPE_VALIDATIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));
    assert_eq!(value.explicit_index_occurrences().unwrap().count(slot), 2);
    assert!(Arc::ptr_eq(
        value.proofs.occurrences.get().unwrap(),
        &history
    ));
    let result = SymbolicTensor::multiply_many([&value, &scalar]).unwrap();
    assert_eq!(result.explicit_index_occurrences().unwrap().count(slot), 2);
    let mut unpublished = value;
    unpublished.proofs = TensorProofs::default();
    assert!(unpublished.explicit_index_occurrences().is_err());
    assert_eq!(OCCURRENCE_FUNCTION_VISITS.with(|count| count.get()), 0);
    assert_eq!(SCOPE_VALIDATIONS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

#[test]
fn checked_constructors_preserve_two_declared_zero_ports_with_the_same_index() {
    test_support::test_initialize();
    let slot = LibraryRep::from(Minkowski {})
        .new_rep(4)
        .slot(PartialIndex::Explicit(AbstractIndex::Normal(99597)));
    let structure = PartialStructure::from_logical_slots([slot, slot]);
    for value in [
        SymbolicTensor::new(Atom::Zero, structure.clone()).unwrap(),
        SymbolicTensor::from_produced_parts(Atom::Zero, structure.clone()).unwrap(),
    ] {
        assert_eq!(value.rank(), 2);
        assert_eq!(value.structure.logical_slots(), [slot, slot]);
        assert_eq!(
            value
                .explicit_index_occurrences()
                .unwrap()
                .count(slot.rep().slot(AbstractIndex::Normal(99597))),
            2
        );
        let product = SymbolicTensor::multiply_many([&value]).unwrap();
        assert_eq!(product.rank(), 2);
        assert_eq!(product.structure.logical_slots(), [slot, slot]);
    }
}

#[test]
fn dummy_scoping_without_nested_scopes_adds_no_callback_replay() {
    test_support::test_initialize();
    let calls = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&calls);
    let head = spenso::vector_symbol!(
        "scope_wrap_callback",
        norm = move |_, _| {
            observed.fetch_add(1, Ordering::Relaxed);
        }
    );
    let source = SymbolicTensor::infer(function!(head, spenso::mink!(4, 99965))).unwrap();
    calls.store(0, Ordering::Relaxed);
    let wrapped = source
        .wrap_indices(symbol!("scope_wrap_callback_scope"), true)
        .unwrap();
    assert_eq!(wrapped.expression, source.expression);
    assert_eq!(wrapped.structure, source.structure);
    // The original flat wrapper calls this normalizer once. Scope discovery
    // must not add another call when no nested binding scope exists.
    assert!(
        calls.load(Ordering::Relaxed) <= 1,
        "calls={}",
        calls.load(Ordering::Relaxed)
    );
}

/// Scope separation in expansion freshens closed dummies, also inside a
/// self-loop leaf. Its syntax cannot reproduce its ports, so the exact
/// renaming carries its raw ports, and the result parses without inference.
#[test]
fn expansion_carries_raw_ports_of_freshened_self_loops() {
    test_support::test_initialize();
    let mink = LibraryRep::from(Minkowski {}).new_rep(4);
    let [i, k] = [99641, 99643].map(|index| mink.slot::<AbstractIndex, _>(index));
    let colour = ColorFundamental {}
        .new_rep(3)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99645));
    let trace = tensor_symbol!("authoritative_metadata::self_loop");
    let left = tensor_symbol!("authoritative_metadata::loop_left");
    let right = tensor_symbol!("authoritative_metadata::loop_right");
    let open = tensor_symbol!("authoritative_metadata::loop_open");
    let other = tensor_symbol!("authoritative_metadata::loop_other");
    let x = Atom::var(symbol!("authoritative_metadata::loop_x"));
    for (self_loop, pair) in [
        (
            function!(trace, i.to_atom(), i.to_atom()),
            function!(left, i.to_atom()) * function!(right, i.to_atom()),
        ),
        (
            function!(trace, colour.to_atom(), colour.dual().to_atom()),
            function!(left, colour.to_atom()) * function!(right, colour.dual().to_atom()),
        ),
    ] {
        let mut renamed = false;
        for (case, source) in [
            // The raised loop is one bracketed contraction, `bracket(T(i,i))^2`.
            (&self_loop + &x).pow(2),
            // The neighbouring pair reuses the loop's closed name.
            (&self_loop + &x) * &pair,
            // The square of an open sum freshens each copy's closed loop.
            (&self_loop * function!(open, k.to_atom()) + function!(other, k.to_atom())).pow(2),
        ]
        .into_iter()
        .enumerate()
        {
            let typed = SymbolicTensor::infer(source.clone()).unwrap();
            for via_poly in [false, true] {
                let result = typed.expanded(None, via_poly).unwrap();
                assert_eq!(result.structure, typed.structure, "{source}");
                if case == 0 {
                    let power = function!(SPENSO_TAG.bracket, &self_loop).pow(2);
                    assert!(
                        result.expression.contains(power.as_view()),
                        "{}",
                        result.expression
                    );
                }
                let loops = leaves_with_head(&result.expression, trace);
                assert!(!loops.is_empty(), "{source}");
                for leaf in &loops {
                    renamed |= *leaf != self_loop;
                    assert!(
                        result.proofs.leaf_interfaces.fact(leaf.as_view()).is_some(),
                        "{source} -> {leaf}"
                    );
                }
                parses_without_inference(&result);
                // The carried facts agree with an independent admission.
                let admitted = SymbolicTensor::infer(result.expression.clone()).unwrap();
                assert_eq!(
                    result.canonize().unwrap().expression,
                    admitted.canonize().unwrap().expression,
                    "{source}"
                );
            }
        }
        assert!(renamed, "{self_loop}");
    }
}

/// Each leaf of `value` with head `head`.
fn leaves_with_head(value: &Atom, head: Symbol) -> Vec<Atom> {
    let mut leaves = Vec::new();
    value.as_view().visitor(&mut |node| {
        let found = node.get_symbol() == Some(head);
        if found {
            leaves.push(node.to_owned());
        }
        !found
    });
    leaves
}

/// Parse `value` from its retained facts alone, without any inference.
fn parses_without_inference(value: &SymbolicTensor<PartialStructure>) {
    INFERENCE_CALLS.with(|count| count.set(0));
    INFERENCE_VISITS.with(|count| count.set(0));
    value
        .parse_to_symbolic_net::<AbstractIndex>(&ParseSettings::default())
        .unwrap();
    assert_eq!(INFERENCE_CALLS.with(|count| count.get()), 0);
    assert_eq!(INFERENCE_VISITS.with(|count| count.get()), 0);
}

/// Scope separation carries raw ports at each leaf occurrence that it renames.
/// Copies of a freshened scope rename its leaves again. A neighboring
/// self-loop outside that scope is never rebuilt.
#[test]
fn expansion_carries_raw_ports_at_each_renamed_leaf_occurrence() {
    test_support::test_initialize();
    let mink = LibraryRep::from(Minkowski {}).new_rep(4);
    let i = mink.slot::<AbstractIndex, _>(99671);
    let x = Atom::var(symbol!("authoritative_metadata::occurrence_x"));
    let trace = tensor_symbol!("authoritative_metadata::occurrence_loop");
    let open = tensor_symbol!("authoritative_metadata::occurrence_open");
    let other = tensor_symbol!("authoritative_metadata::occurrence_other");
    let self_loop = function!(trace, i.to_atom(), i.to_atom());

    // The inner product freshens the loop once; the square of the outer sum
    // then freshens that loop again in each copy.
    let nested =
        ((&self_loop + &x) * function!(open, i.to_atom()) + function!(other, i.to_atom())).pow(2);
    let typed = SymbolicTensor::infer(nested.clone()).unwrap();
    for via_poly in [false, true] {
        let result = typed.expanded(None, via_poly).unwrap();
        let loops = leaves_with_head(&result.expression, trace);
        let renamed = loops
            .iter()
            .filter(|leaf| **leaf != self_loop)
            .collect::<HashSet<_>>();
        assert_eq!(renamed.len(), 2, "{nested} -> {}", result.expression);
        for leaf in &loops {
            assert!(
                result.proofs.leaf_interfaces.fact(leaf.as_view()).is_some(),
                "{nested} -> {leaf}"
            );
        }
        let leaf_facts = result.proofs.parser_leaf_interfaces();
        for leaf in leaf_facts.lock().unwrap().keys() {
            if leaf.as_view().get_symbol() == Some(trace) {
                assert!(loops.contains(leaf), "intermediate loop {leaf}");
            }
        }
        parses_without_inference(&result);
    }

    // The neighbouring self-loop keeps its name, so it is never rebuilt and
    // no fact is published for a renamed copy of it.
    let rebuilt = Arc::new(AtomicUsize::new(0));
    let observed = Arc::clone(&rebuilt);
    let original = i.to_atom();
    let counted = tensor_symbol!(
        "authoritative_metadata::occurrence_counted",
        norm = move |node, _| {
            if let AtomView::Fun(function) = node
                && !function
                    .iter()
                    .all(|argument| argument == original.as_view())
            {
                observed.fetch_add(1, Ordering::Relaxed);
            }
        }
    );
    let neighbour = function!(counted, i.to_atom(), i.to_atom());
    let source = (&self_loop + &x) * &neighbour;
    let typed = SymbolicTensor::infer(source.clone()).unwrap();
    for via_poly in [false, true] {
        rebuilt.store(0, Ordering::Relaxed);
        let result = typed.expanded(None, via_poly).unwrap();
        assert_eq!(rebuilt.load(Ordering::Relaxed), 0, "{source}");
        let neighbours = leaves_with_head(&result.expression, counted);
        assert!(!neighbours.is_empty());
        assert!(
            neighbours.iter().all(|leaf| *leaf == neighbour),
            "{}",
            result.expression
        );
        assert!(
            leaves_with_head(&result.expression, trace)
                .iter()
                .any(|leaf| *leaf != self_loop)
        );
        let facts = result.proofs.parser_leaf_interfaces();
        assert!(
            facts
                .lock()
                .unwrap()
                .keys()
                .all(|leaf| { leaf.as_view().get_symbol() != Some(counted) || *leaf == neighbour })
        );
        parses_without_inference(&result);
    }
}

#[test]
fn expansion_carries_raw_ports_through_an_antisymmetric_sign_flip() {
    test_support::test_initialize();
    let colour = ColorFundamental {}.new_rep(3);
    let [inner, spectator] = [10, 50].map(|index| colour.slot::<AbstractIndex, _>(index));
    let [up, fixed, down] = [inner.to_atom(), spectator.to_atom(), inner.dual().to_atom()];
    let antisymmetric = tensor_symbol!(
        "authoritative_metadata::occurrence_antisymmetric"; Antisymmetric
    );
    let [vector, left, right] = [
        tensor_symbol!("authoritative_metadata::occurrence_vector"),
        tensor_symbol!("authoritative_metadata::occurrence_left"),
        tensor_symbol!("authoritative_metadata::occurrence_right"),
    ];
    // Keep coefficient +1 initially. Both up slots share a head, while
    // the down wrapper sorts by its own head and stays in place.
    let mut args = [&up, &fixed, &down];
    args.sort_unstable_by(|a, b| a.as_view().cmp(&b.as_view()));
    let leaf = function!(antisymmetric, args[0], args[1], args[2]);
    assert!(matches!(leaf.as_view(), AtomView::Fun(_)));
    assert!(up.as_view().cmp(&fixed.as_view()).is_lt());

    let x = Atom::var(symbol!("authoritative_metadata::occurrence_sign_x"));
    let spectator = function!(vector, &fixed);
    let left = function!(left, &up);
    let right = function!(right, &down);
    let source = (&leaf + &x * &spectator) * &left * &right;
    let typed = SymbolicTensor::infer(source).unwrap();
    for via_poly in [false, true] {
        let result = typed.expanded(None, via_poly).unwrap();
        let leaves = leaves_with_head(&result.expression, antisymmetric);
        let [renamed] = leaves.as_slice() else {
            panic!("expected one antisymmetric leaf in {}", result.expression);
        };
        let AtomView::Fun(function) = renamed.as_view() else {
            unreachable!("a tensor leaf is a function")
        };
        let ports = function
            .iter()
            .map(|port| port.to_owned())
            .collect::<Vec<_>>();
        let renamed_up = ports
            .iter()
            .find(|port| {
                **port != fixed
                    && Slot::<LibraryRep, AbstractIndex>::try_from(port.as_view())
                        .is_ok_and(|slot| !slot.rep_name().is_dual())
            })
            .unwrap();
        assert!(fixed.as_view().cmp(&renamed_up.as_view()).is_lt());
        assert_eq!(
            result.expression,
            &x * &spectator * &left * &right - renamed * &left * &right,
            "the dummy crosses the spectator in one transposition"
        );
        assert_eq!(
            result.proofs.leaf_interfaces.fact(renamed.as_view()),
            Some((true, ports))
        );
        parses_without_inference(&result);
    }
}

/// Expansion renames the dummies of direct user leaves, which keep no raw-port
/// fact. Regional reads still compose their written interfaces, also through
/// a raised closed product, so contraction finds every region it opens.
#[test]
fn expansion_renamed_direct_leaves_keep_regional_interfaces() {
    test_support::test_initialize();
    let mink = LibraryRep::from(Minkowski {}).new_rep(4);
    let i = mink.slot::<AbstractIndex, _>(99681);
    let [p, q, s] = [
        tensor_symbol!("authoritative_metadata::regional_p"),
        tensor_symbol!("authoritative_metadata::regional_q"),
        tensor_symbol!("authoritative_metadata::regional_s"),
    ];
    let x = Atom::var(symbol!("authoritative_metadata::regional_x"));
    let pair = function!(p, i.to_atom()) * function!(q, i.to_atom());
    let spectator = function!(s, i.to_atom());
    let open =
        PartialStructure::from_logical_slots([i.rep().slot(PartialIndex::Explicit(i.aind()))]);
    // Typed powers keep the closed product atomic, `bracket(P(i) Q(i))^2`.
    let [x, spectator] = [x, spectator].map(|value| SymbolicTensor::infer(value).unwrap());
    let pair = SymbolicTensor::infer(pair).unwrap();
    let two = SymbolicTensor::infer(Atom::num(2)).unwrap();
    for typed in [
        pair.try_add(&x).unwrap().multiply(&spectator).unwrap(),
        pair.pow(&two)
            .unwrap()
            .try_add(&x)
            .unwrap()
            .multiply(&spectator)
            .unwrap(),
    ] {
        let source = &typed.expression;
        for via_poly in [false, true] {
            let result = typed.expanded(None, via_poly).unwrap();
            let leaves = leaves_with_head(&result.expression, p);
            assert!(
                leaves.iter().all(|leaf| *leaf != function!(p, i.to_atom())),
                "{source}"
            );
            for leaf in &leaves {
                let AtomView::Fun(function) = leaf.as_view() else {
                    unreachable!("a tensor leaf is a function")
                };
                let slot = Slot::<LibraryRep, AbstractIndex>::try_from(function.get(0)).unwrap();
                let written = PartialStructure::from_logical_slots([slot
                    .rep()
                    .slot(PartialIndex::Explicit(slot.aind()))]);
                assert!(result.proofs.leaf_interfaces.fact(leaf.as_view()).is_none());
                assert_eq!(
                    result.proofs.inspect_interface(leaf.as_view()),
                    Some(written.clone())
                );
                assert_eq!(
                    result.proofs.admitted_interface(leaf.as_view()),
                    Some(written)
                );
            }
            // Each enclosing term, including a raised closed product, keeps
            // the spectator's open port.
            let AtomView::Add(sum) = result.expression.as_view() else {
                panic!("expansion keeps two terms: {}", result.expression)
            };
            for term in sum.iter() {
                assert_eq!(
                    result.proofs.inspect_interface(term),
                    Some(open.clone()),
                    "{term}"
                );
                assert_eq!(
                    result.proofs.admitted_interface(term),
                    Some(open.clone()),
                    "{term}"
                );
            }
            let contracted = result.contract(ContractSettings::default()).unwrap();
            assert_eq!(contracted.structure, typed.structure);
            parses_without_inference(&result);
        }
    }
}

/// A symmetric leaf whose declared port order differs from its written one
/// keeps that order when canonical relabelling moves its dummy across its
/// fixed arguments: the template binds modulo the argument permutation.
#[test]
fn canonical_relabelling_keeps_the_declared_order_of_symmetric_leaves() {
    test_support::test_initialize();
    let mink = LibraryRep::from(Minkowski {}).new_rep(4);
    let [i, j, k, h] = [10, 50, 70, 90].map(|index| mink.slot::<AbstractIndex, _>(index));
    let explicit = |slot: Slot<LibraryRep, AbstractIndex>| {
        slot.rep().slot(PartialIndex::Explicit(slot.aind()))
    };
    let symmetric = tensor_symbol!("authoritative_metadata::declared_symmetric"; Symmetric);
    let [p, u, q] = [
        tensor_symbol!("authoritative_metadata::declared_p"),
        tensor_symbol!("authoritative_metadata::declared_u"),
        tensor_symbol!("authoritative_metadata::declared_q"),
    ];
    let leaf = function!(symmetric, j.to_atom(), k.to_atom(), i.to_atom());
    let AtomView::Fun(written) = leaf.as_view() else {
        unreachable!("a tensor leaf is a function")
    };
    let initial_position = written
        .iter()
        .position(|argument| argument == i.to_atom().as_view());
    let declared = [j, i, k];
    assert!(
        written
            .iter()
            .ne(declared.iter().map(|slot| slot.to_atom()))
    );
    let declared_ports = declared.map(|slot| slot.to_atom()).to_vec();
    let tensor = SymbolicTensor::from_declared_parts(
        leaf.clone(),
        PartialStructure::from_logical_slots(declared.map(explicit)),
    )
    .unwrap();
    tensor
        .proofs
        .retain_leaf_ports(leaf.as_view(), declared_ports.clone())
        .unwrap();
    let first = tensor
        .multiply(&SymbolicTensor::infer(function!(p, i.to_atom())).unwrap())
        .unwrap();
    let second = SymbolicTensor::infer(
        function!(u, j.to_atom(), k.to_atom(), h.to_atom()) * function!(q, h.to_atom()),
    )
    .unwrap();
    let sum = first.try_add(&second).unwrap();
    assert_eq!(
        sum.proofs.leaf_interfaces.fact(leaf.as_view()),
        Some((true, declared_ports))
    );
    let canonical = sum.canonize().unwrap();
    let leaves = leaves_with_head(&canonical.expression, symmetric);
    let [relabelled] = leaves.as_slice() else {
        panic!("one symmetric leaf: {}", canonical.expression)
    };
    let AtomView::Fun(function) = relabelled.as_view() else {
        unreachable!("a tensor leaf is a function")
    };
    let dummy = function
        .iter()
        .find(|argument| {
            ![j.to_atom(), k.to_atom()]
                .iter()
                .any(|fixed| fixed.as_view() == *argument)
        })
        .unwrap()
        .to_owned();
    // Relabelling moves the dummy in the written argument list.
    assert_ne!(
        function
            .iter()
            .position(|argument| argument == dummy.as_view()),
        initial_position
    );
    assert_eq!(
        canonical.proofs.leaf_interfaces.fact(relabelled.as_view()),
        Some((true, vec![j.to_atom(), dummy, k.to_atom()]))
    );
    parses_without_inference(&canonical);
}

#[test]
fn canonical_relabelling_rotates_cyclic_self_loops() {
    test_support::test_initialize();
    let head = tensor_symbol!("authoritative_metadata::cyclic_loop"; Cyclesymmetric);
    let [i, j] = [mink!(4, 1), mink!(4, 2)];
    let source = SymbolicTensor::infer(function!(head, &i, &j, &i)).unwrap();
    let canonical = source.canonize().unwrap();
    let leaves = leaves_with_head(&canonical.expression, head);
    let [leaf] = leaves.as_slice() else {
        panic!("one cyclic leaf: {}", canonical.expression);
    };
    let AtomView::Fun(function) = leaf.as_view() else {
        unreachable!("a tensor leaf is a function");
    };
    let args = function.iter().collect::<Vec<_>>();
    assert_eq!(args[0], j.as_view());
    assert_eq!(args[1], args[2]);
    assert_ne!(args[1], i.as_view());
    let (tensorial, ports) = canonical
        .proofs
        .leaf_interfaces
        .fact(leaf.as_view())
        .unwrap();
    assert!(tensorial);
    assert_eq!(ports.iter().map(Atom::as_view).collect::<Vec<_>>(), args);
    assert_eq!(canonical.structure, source.structure);
    parses_without_inference(&canonical);
}

#[test]
fn canonical_relabelling_retains_written_ports_of_symmetric_and_cyclic_self_loops() {
    test_support::test_initialize();
    let symmetric = tensor_symbol!("authoritative_metadata::symmetric_loops"; Symmetric);
    let cyclic = tensor_symbol!("authoritative_metadata::cyclic_loops"; Cyclesymmetric);
    let p = tensor_symbol!("authoritative_metadata::loop_p");
    let q = tensor_symbol!("authoritative_metadata::loop_q");
    let [i, k, l] = [mink!(4, 1), mink!(4, 2), mink!(4, 3)];

    for (head, expression) in [
        (symmetric, function!(symmetric, &i, &i, &k, &k)),
        (
            symmetric,
            function!(symmetric, &i, &i, &k, &l) * function!(p, &k) * function!(q, &l),
        ),
        (cyclic, function!(cyclic, &i, &i, &k, &k)),
        (cyclic, function!(cyclic, &i, &k, &i, &k)),
    ] {
        let source = SymbolicTensor::infer(expression).unwrap();
        let canonical = source
            .canonize()
            .unwrap_or_else(|error| panic!("{}: {error}", source.expression));
        let leaves = leaves_with_head(&canonical.expression, head);
        let [leaf] = leaves.as_slice() else {
            panic!("one symmetric or cyclic leaf: {}", canonical.expression);
        };
        let AtomView::Fun(function) = leaf.as_view() else {
            unreachable!("a tensor leaf is a function");
        };
        let written = function.iter().map(|port| port.to_owned()).collect();
        assert_eq!(
            canonical.proofs.leaf_interfaces.fact(leaf.as_view()),
            Some((true, written))
        );
        assert_eq!(canonical.structure, source.structure);
        parses_without_inference(&canonical);
    }
}

#[test]
fn canonical_relabelling_rejects_ambiguous_declared_symmetric_ports() {
    test_support::test_initialize();
    let mink = LibraryRep::from(Minkowski {}).new_rep(4);
    let symmetric = tensor_symbol!("authoritative_metadata::ambiguous_symmetric"; Symmetric);
    let p = tensor_symbol!("authoritative_metadata::ambiguous_p");
    let q = tensor_symbol!("authoritative_metadata::ambiguous_q");

    for indices in [
        [AbstractIndex::Normal(1), AbstractIndex::Normal(2)],
        [AbstractIndex::Dummy(0), AbstractIndex::Dummy(1)],
    ] {
        let [i, j] = indices.map(|index| mink.slot::<AbstractIndex, _>(index).to_atom());
        let leaf = function!(symmetric, &i, &j);
        let AtomView::Fun(function) = leaf.as_view() else {
            unreachable!("a tensor leaf is a function");
        };
        let mut ports = function
            .iter()
            .map(|port| port.to_owned())
            .collect::<Vec<_>>();
        ports.reverse();
        let symmetric_tensor = SymbolicTensor::infer(leaf.clone()).unwrap();
        symmetric_tensor
            .proofs
            .retain_leaf_ports(leaf.as_view(), ports)
            .unwrap();
        let p = SymbolicTensor::infer(function!(p, &i)).unwrap();
        let q = SymbolicTensor::infer(function!(q, &j)).unwrap();
        let source = symmetric_tensor.multiply(&p).unwrap().multiply(&q).unwrap();

        let error = source.canonize().unwrap_err();
        assert!(
            error
                .to_string()
                .contains("ambiguous across canonicalized dummies"),
            "{error:?}"
        );
    }
}

#[test]
fn producer_routing_carries_gamma_ports_when_a_momentum_sign_moves_outside() {
    let representations = test_support::test_initialize();
    let spin = representations.bis4.to_symbolic([]);
    let mink = representations.mink4.to_symbolic([]);
    let q = vector_symbol!("authoritative_metadata::edge_momentum");
    let p = vector_symbol!("authoritative_metadata::routed_momentum");
    let momentum = function!(q, 1, &mink);
    let routed = function!(p, 0, &mink);
    let source = SymbolicTensor::infer(trace!(&spin, gamma!(&momentum))).unwrap();
    let mapped = source
        .map_produced_expression(&HashMap::new(), |value| {
            Ok::<_, TensorInferenceError>((value == momentum.as_view()).then(|| -&routed))
        })
        .unwrap();

    let gamma = gamma!(&routed);
    assert!(mapped.expression.contains(gamma.as_view()));
    assert_eq!(
        mapped.proofs.leaf_interfaces.fact(gamma.as_view()),
        Some((
            true,
            vec![
                Atom::var(SPENSO_TAG.chain_in),
                Atom::var(SPENSO_TAG.chain_out),
            ],
        ))
    );
    assert!(
        mapped
            .proofs
            .inspect_interface((-gamma).as_view())
            .is_some_and(|interface| interface.canonical().is_scalar())
    );
}

#[test]
fn canonical_relabelling_rejects_ambiguous_ports_with_another_matching_template() {
    test_support::test_initialize();
    let symmetric = tensor_symbol!("authoritative_metadata::shared_symmetric_template"; Symmetric);
    let p = tensor_symbol!("authoritative_metadata::template_p");
    let q = tensor_symbol!("authoritative_metadata::template_q");
    let r = tensor_symbol!("authoritative_metadata::template_r");
    let t = tensor_symbol!("authoritative_metadata::template_t");
    let [i, j, k, l] = [mink!(4, 1), mink!(4, 2), mink!(4, 3), mink!(4, 4)];
    let reversed = function!(symmetric, &i, &j);
    let written = function!(symmetric, &k, &l);
    let expression = &reversed
        * function!(p, &i)
        * function!(q, &j)
        * &written
        * function!(r, &k)
        * function!(t, &l);
    let source = SymbolicTensor::infer(expression).unwrap();
    let Some((true, mut ports)) = source.proofs.leaf_interfaces.fact(reversed.as_view()) else {
        panic!("missing source port fact for {reversed}");
    };
    ports.reverse();
    source
        .proofs
        .retain_leaf_ports(reversed.as_view(), ports.clone())
        .unwrap();
    assert_eq!(
        source.proofs.leaf_interfaces.fact(reversed.as_view()),
        Some((true, ports))
    );
    assert!(
        source
            .proofs
            .leaf_interfaces
            .fact(written.as_view())
            .is_some()
    );

    let error = source.canonize().unwrap_err();
    assert!(
        error
            .to_string()
            .contains("ambiguous across canonicalized dummies"),
        "{error:?}"
    );
}

#[test]
fn selected_expansion_brackets_closed_tensor_functions_with_selected_metadata() {
    test_support::test_initialize();
    let c = Atom::var(symbol!("authoritative_metadata::selected_coefficient"));
    let x = Atom::var(symbol!("authoritative_metadata::selected_metadata"));
    let head = tensor_symbol!("authoritative_metadata::selected_loop");
    let i = mink!(4, 1);
    let closed = function!(head, &x, &i, &i);
    let source = SymbolicTensor::infer((&c * &closed + &x).pow(2)).unwrap();
    let squared =
        SymbolicTensor::<PartialStructure>::scalar_power_expression(closed, Some(true)).pow(2);
    for via_poly in [false, true] {
        let result = source.expanded(Some(x.as_view()), via_poly).unwrap();
        let mut found = false;
        result.expression.as_view().visitor(&mut |value| {
            found |= value == squared.as_view();
            !found
        });
        assert!(found, "poly={via_poly}: {}", result.expression);
        SymbolicTensor::infer(result.expression.clone())
            .unwrap()
            .canonize()
            .unwrap();
        result.canonize().unwrap();
    }
}

#[test]
fn canonical_relabelling_checks_exact_byte_hits_against_other_summand_templates() {
    test_support::test_initialize();
    let mink = LibraryRep::from(Minkowski {}).new_rep(4);
    // Use canonical dummy indices so collapsing both summands reproduces one
    // input leaf's exact bytes in the output's retained-interface table.
    let [d0, d1] = [0, 1].map(|index| {
        mink.slot::<AbstractIndex, _>(AbstractIndex::Dummy(index))
            .to_atom()
    });
    let external = euc!(3, 99681);
    let p = tensor_symbol!("authoritative_metadata::exact_byte_p");
    let symmetric = tensor_symbol!("authoritative_metadata::exact_byte_sym"; Symmetric);
    let ordinary = tensor_symbol!("authoritative_metadata::exact_byte_plain");

    for head in [symmetric, ordinary] {
        let first = function!(head, &d0, &external);
        let second = function!(head, &d1, &external);
        let left = SymbolicTensor::infer(&first * function!(p, &d0)).unwrap();
        let right = SymbolicTensor::infer(&second * function!(p, &d1)).unwrap();

        for reverse_first in [false, true] {
            let conflicting = left.try_add(&right).unwrap();
            let flipped = if reverse_first { &first } else { &second };
            let mut ports = conflicting
                .proofs
                .leaf_interfaces
                .fact(flipped.as_view())
                .unwrap()
                .1;
            ports.reverse();
            conflicting
                .proofs
                .retain_leaf_ports(flipped.as_view(), ports)
                .unwrap();
            let error = conflicting.canonize().unwrap_err();
            assert!(
                error
                    .to_string()
                    .contains("ambiguous across canonicalized dummies"),
                "{head}: {error:?}"
            );
        }

        let consistent = left.try_add(&right).unwrap();
        let mut expected = consistent
            .proofs
            .leaf_interfaces
            .fact(first.as_view())
            .unwrap()
            .1;
        expected.reverse();
        for leaf in [&first, &second] {
            let mut ports = consistent
                .proofs
                .leaf_interfaces
                .fact(leaf.as_view())
                .unwrap()
                .1;
            ports.reverse();
            consistent
                .proofs
                .retain_leaf_ports(leaf.as_view(), ports)
                .unwrap();
        }
        let canonical = consistent.canonize().unwrap();
        assert_eq!(
            canonical.expression,
            2 * &first * function!(p, &d0),
            "{head}"
        );
        assert_eq!(
            canonical.proofs.leaf_interfaces.fact(first.as_view()),
            Some((true, expected)),
            "{head}"
        );
    }
}
