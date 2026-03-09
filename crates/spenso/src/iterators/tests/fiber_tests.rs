//! Tests for fiber abstractions
//!
//! Contains tests for fiber behavior and operations.

use crate::iterators::{AbstractFiber, AbstractFiberIndex, Fiber, FiberClass, FiberIndex};
use crate::structure::representation::RepName;
use crate::structure::TensorStructure;
use crate::structure::{representation::Euclidean, OrderedStructure};

#[test]
fn test_fiber_creation() {
    let rep = Euclidean {};
    let structure: OrderedStructure<Euclidean> = OrderedStructure::new(vec![
        rep.new_slot(3, 0),
        rep.new_slot(4, 0),
        rep.new_slot(5, 0),
    ])
    .structure;

    // Test creation from boolean filter
    let filter = [true, false, true];
    let fiber = Fiber::from_filter(&filter, &structure);

    // Check that the pattern of fixed/free indices matches the filter
    assert_eq!(fiber.bitvec().len(), filter.len());
    for (i, &is_free) in filter.iter().enumerate() {
        assert_eq!(fiber[i].is_free(), is_free);
    }

    // Test creation from flat index
    let flat_idx = 5.into(); // some arbitrary flat index
    let fiber_from_flat = Fiber::from_flat(flat_idx, &structure);

    // All indices should be fixed when created from flat index
    for i in 0..structure.order() {
        assert!(fiber_from_flat[i].is_fixed());
    }

    // Test zeros constructor
    let zeros_fiber = Fiber::zeros(&structure);
    for i in 0..structure.order() {
        assert!(zeros_fiber[i].is_fixed());
        if let FiberIndex::Fixed(val) = zeros_fiber[i] {
            assert_eq!(val, 0);
        }
    }
}

#[test]
fn test_fiber_modification() {
    let rep = Euclidean {};
    let structure: OrderedStructure<Euclidean> = OrderedStructure::new(vec![
        rep.new_slot(3, 0),
        rep.new_slot(4, 0),
        rep.new_slot(5, 0),
    ])
    .structure;

    // Start with all zeros
    let mut fiber = Fiber::zeros(&structure);

    // Modify one index to be free
    fiber.free(1);
    assert!(fiber[1].is_free());

    // Fix an index to a specific value
    fiber.fix(0, 2);
    if let FiberIndex::Fixed(val) = fiber[0] {
        assert_eq!(val, 2);
    } else {
        panic!("Expected index 0 to be fixed");
    }

    // Test is_single() for detecting exactly one free index
    let single_idx = fiber.is_single();
    if let FiberIndex::Fixed(val) = single_idx {
        assert_eq!(val, 1); // Index 1 should be the only free one
    } else {
        panic!("Expected is_single() to identify index 1 as the only free index");
    }

    // Now make another index free
    fiber.free(2);
    assert!(fiber[2].is_free());

    // is_single() should now return Free since there are multiple free indices
    // assert!(matches!(fiber.is_single(), FiberIndex::Free));
}

#[test]
fn test_fiber_class_conversion() {
    let rep = Euclidean {};
    let structure: OrderedStructure<Euclidean> = OrderedStructure::new(vec![
        rep.new_slot(3, 0),
        rep.new_slot(4, 0),
        rep.new_slot(5, 0),
    ])
    .structure;

    // Create a fiber with mixed fixed/free indices
    let filter = [true, false, true];
    let fiber = Fiber::from_filter(&filter, &structure);

    // Convert to fiber class
    let fiber_class: FiberClass<_> = fiber.clone().into();

    // The bitvector pattern should be inverted for fiber classes
    let fiber_bits = fiber.bitvec();
    let class_bits = fiber_class.bitvec();

    assert_eq!(fiber_bits.len(), class_bits.len());
    for i in 0..fiber_bits.len() {
        assert_eq!(fiber_bits[i], !class_bits[i]);
    }

    // Convert back to fiber
    let fiber_back: Fiber<_> = fiber_class.into();

    // The pattern should be preserved
    let back_bits = fiber_back.bitvec();
    assert_eq!(fiber_bits, back_bits);
}

#[test]
fn test_fiber_display() {
    let rep = Euclidean {};
    let structure: OrderedStructure<Euclidean> =
        OrderedStructure::new(vec![rep.new_slot(3, 0), rep.new_slot(4, 0)]).structure;

    // Create a fiber with specific pattern
    let mut fiber = Fiber::zeros(&structure);
    fiber.fix(0, 2);
    fiber.free(1);

    // Display should show "2 :"
    let display = format!("{}", fiber);
    assert_eq!(display.trim(), "2 :");
}

#[test]
fn test_single_fiber_detection() {
    let rep = Euclidean {};
    let structure: OrderedStructure<Euclidean> = OrderedStructure::new(vec![
        rep.new_slot(3, 0),
        rep.new_slot(4, 0),
        rep.new_slot(5, 0),
    ])
    .structure;

    // Create fibers with different numbers of free indices
    let mut zero_free = Fiber::zeros(&structure);
    let mut one_free = Fiber::zeros(&structure);
    let mut two_free = Fiber::zeros(&structure);

    one_free.free(1);
    two_free.free(0);
    two_free.free(2);

    // Test single detection
    assert!(matches!(zero_free.is_single(), FiberIndex::Free));

    if let FiberIndex::Fixed(val) = one_free.is_single() {
        assert_eq!(val, 1);
    } else {
        panic!("Expected is_single() to identify index 1");
    }

    assert!(matches!(two_free.is_single(), FiberIndex::Free));
}
