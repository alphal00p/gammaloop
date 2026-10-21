use super::*;
use symbolica::function;

#[test]
fn expansion_keeps_coalesced_oriented_tensor_copies() {
    crate::test_support::test_initialize();
    let slot = LibraryRep::from(ColorFundamental {})
        .new_rep(3)
        .slot::<AbstractIndex, _>(AbstractIndex::Normal(99489));
    let vector = SymbolicTensor::infer(function!(
        spenso::tensor_symbol!("authoritative_metadata::expanded_copy"),
        slot.to_atom()
    ))
    .unwrap();
    let x = Atom::var(symbolica::symbol!("authoritative_metadata::expanded_x"));
    let y = Atom::var(symbolica::symbol!("authoritative_metadata::expanded_y"));
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
}

#[test]
fn zero_algebra_maps_hidden_pair_before_later_explicit_collisions() {
    crate::test_support::test_initialize();
    let index = AbstractIndex::Normal(99471);
    let minkowski = LibraryRep::from(Minkowski {});
    let physical = minkowski.new_rep(4).slot::<AbstractIndex, _>(index);
    let extended = minkowski.new_rep(6).slot::<AbstractIndex, _>(index);
    // Ordinary tensor heads keep the admitted pair in product notation.
    let source = SymbolicTensor::infer(
        function!(
            spenso::tensor_symbol!("authoritative_metadata::hidden_a"),
            physical.to_atom()
        ) * function!(
            spenso::tensor_symbol!("authoritative_metadata::hidden_b"),
            physical.to_atom()
        ),
    )
    .unwrap();
    let physical_leaf = SymbolicTensor::infer(function!(
        spenso::tensor_symbol!("authoritative_metadata::physical_leaf"),
        physical.to_atom()
    ))
    .unwrap();
    let extended_leaf = SymbolicTensor::infer(function!(
        spenso::tensor_symbol!("authoritative_metadata::extended_leaf"),
        extended.to_atom()
    ))
    .unwrap();
    assert!(source.is_scalar());
    assert_eq!(
        source.explicit_index_occurrences().unwrap().count(physical),
        2
    );
    physical_leaf.ensure_validated().unwrap();
    extended_leaf.ensure_validated().unwrap();
    OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
    SCOPE_VALIDATIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));

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
}

#[test]
fn dual_free_dimension_mapping_preserves_orientation_and_later_collision_history() {
    crate::test_support::test_initialize();
    let index = AbstractIndex::Normal(99473);
    let dual_three = LibraryRep::from(ColorAntiFundamental {}).new_rep(3);
    let dual_four = LibraryRep::from(ColorAntiFundamental {}).new_rep(4);
    let source_slot = dual_three.slot::<AbstractIndex, _>(index);
    let mapped_slot = dual_four.slot::<AbstractIndex, _>(index);
    let partner_slot = dual_four.dual().slot::<AbstractIndex, _>(index);
    let head = spenso::tensor_symbol!("authoritative_metadata::dual_free");
    let source = SymbolicTensor::infer(function!(head, source_slot.to_atom())).unwrap();
    let partner = SymbolicTensor::infer(function!(
        spenso::tensor_symbol!("authoritative_metadata::base_partner"),
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
    partner.ensure_validated().unwrap();
    OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
    SCOPE_VALIDATIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));

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
}

#[test]
fn admitted_arithmetic_keeps_logical_ports_and_composes_without_observation() {
    crate::test_support::test_initialize();
    let minkowski = LibraryRep::from(Minkowski {}).new_rep(4);
    let indices = [99479, 99477].map(AbstractIndex::Normal);
    let slots = indices.map(|index| minkowski.slot(index));
    let ports = indices.map(|index| minkowski.slot(PartialIndex::Explicit(index)));
    // Include a callback-sensitive leaf and a logical order different from the
    // other operand's order. Admission owns both before the counters are reset.
    let a = SymbolicTensor::from_declared_parts(
        FunctionBuilder::new(spenso::tensor_symbol!(
            "authoritative_metadata::arithmetic_a",
            norm = |_, _| {}
        ))
        .add_args(ports.map(composition::port_atom))
        .finish(),
        PartialStructure::from_logical_slots(ports),
    )
    .unwrap();
    let b = SymbolicTensor::from_declared_parts(
        FunctionBuilder::new(spenso::tensor_symbol!(
            "authoritative_metadata::arithmetic_b"
        ))
        .add_args(ports.into_iter().rev().map(composition::port_atom))
        .finish(),
        PartialStructure::from_logical_slots(ports.into_iter().rev()),
    )
    .unwrap();
    let scalar = SymbolicTensor::infer(Atom::var(symbolica::symbol!(
        "authoritative_metadata::weight"
    )))
    .unwrap();
    let scalar_zero = SymbolicTensor::infer(Atom::Zero).unwrap();
    let probe = SymbolicTensor::infer(function!(
        spenso::tensor_symbol!("authoritative_metadata::arithmetic_probe"),
        slots[0].to_atom()
    ))
    .unwrap();
    for source in [&a, &b, &scalar, &scalar_zero, &probe] {
        source.ensure_validated().unwrap();
    }
    OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
    SCOPE_VALIDATIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));

    let sum = a.try_add(&b).unwrap();
    let difference = a.try_sub(&b).unwrap();
    let negated = -sum.clone();
    let weighted = sum.multiply(&scalar).unwrap();
    let cancelled = a.try_sub(&a).unwrap();
    let scaled_zero = sum.multiply(&scalar_zero).unwrap();
    assert!(cancelled.is_zero());
    assert!(scaled_zero.is_zero());
    for result in [sum, difference, negated, weighted, cancelled, scaled_zero] {
        assert_eq!(result.structure.logical_slots(), ports);
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
}

#[test]
fn additive_branches_keep_maximum_public_multiplicity_in_both_operand_orders() {
    crate::test_support::test_initialize();
    let index = AbstractIndex::Normal(99481);
    let slot = LibraryRep::from(Minkowski {})
        .new_rep(4)
        .slot::<AbstractIndex, _>(index);
    let scalar = SymbolicTensor::infer(slot.to_atom()).unwrap();
    let vector = SymbolicTensor::infer(function!(
        spenso::tensor_symbol!("authoritative_metadata::branch_vector"),
        slot.to_atom()
    ))
    .unwrap();
    let extra = SymbolicTensor::infer(function!(
        spenso::tensor_symbol!("authoritative_metadata::branch_extra"),
        slot.to_atom()
    ))
    .unwrap();
    assert!(scalar.is_scalar());
    assert_eq!(scalar.explicit_index_occurrences().unwrap().count(slot), 1);
    assert_eq!(vector.explicit_index_occurrences().unwrap().count(slot), 1);
    extra.ensure_validated().unwrap();
    OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
    SCOPE_VALIDATIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));

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
}

#[test]
fn mapping_hidden_base_use_separates_it_from_the_free_dual_port() {
    crate::test_support::test_initialize();
    let index = AbstractIndex::Normal(99483);
    let base_three = LibraryRep::from(ColorFundamental {}).new_rep(3);
    let base_four = LibraryRep::from(ColorFundamental {}).new_rep(4);
    let base_slot = base_three.slot::<AbstractIndex, _>(index);
    let mapped_base_slot = base_four.slot::<AbstractIndex, _>(index);
    let dual_slot = base_three.dual().slot::<AbstractIndex, _>(index);
    let scalar_base = SymbolicTensor::infer(base_slot.to_atom()).unwrap();
    let head = spenso::tensor_symbol!("authoritative_metadata::hidden_base_dual");
    let free_dual = SymbolicTensor::infer(function!(head, dual_slot.to_atom())).unwrap();
    let compatible = SymbolicTensor::infer(function!(
        spenso::tensor_symbol!("authoritative_metadata::hidden_base_compatible"),
        base_slot.to_atom()
    ))
    .unwrap();
    // A bare base slot is scalar metadata; its dual is not the same fixture.
    assert!(scalar_base.is_scalar());
    assert_eq!(
        scalar_base
            .explicit_index_occurrences()
            .unwrap()
            .count(base_slot),
        1
    );
    free_dual.ensure_validated().unwrap();
    compatible.ensure_validated().unwrap();
    OCCURRENCE_FUNCTION_VISITS.with(|count| count.set(0));
    SCOPE_VALIDATIONS.with(|count| count.set(0));
    INFERENCE_CALLS.with(|count| count.set(0));

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

    // Protect the exact external variance, rather than hiding the opposite
    // orientation's use merely because it shares this collision key.
    let scope = symbolica::symbol!("authoritative_metadata::hidden_base_scope");
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
        spenso::tensor_symbol!("authoritative_metadata::scoped_base_probe"),
        scoped_slot.to_atom()
    ))
    .unwrap();
    let twice = scoped.multiply(&probe).unwrap();
    assert!(matches!(twice.multiply(&probe),
        Err(TensorCompositionError::InvalidExplicitMultiplicity { index: found, occurrences: 3 })
            if found == index.scoped(scope)
    ));
}
