//! Tests for core iterator implementations
//!
//! Contains tests for the low-level iterator types that traverse tensor dimensions.

use crate::iterators::{
    CoreExpandedFiberIterator, CoreFlatFiberIterator, Fiber, IteratesAlongFibers,
    IteratesAlongPermutedFibers, MetricFiberIterator, MetricItem, ResetableIterator,
    ShiftableIterator,
};
use crate::structure::representation::{Minkowski, RepName};
use crate::structure::PermutedStructure;
use crate::structure::{concrete_index::FlatIndex, representation::Euclidean, OrderedStructure};
use linnet::permutation::Permutation;

#[test]
fn test_core_flat_iterator() {
    let rep = Euclidean {};
    let structure: OrderedStructure<Euclidean> =
        OrderedStructure::new(vec![rep.new_slot(2, 0), rep.new_slot(3, 0)]).structure;

    // Create a fiber with one free index
    let mut fiber = Fiber::zeros(&structure);
    fiber.free(0); // Free the first dimension (of size 2)

    // Create a flat iterator
    let mut iter = CoreFlatFiberIterator::new(&fiber, false);

    // Collect all indices
    let indices: Vec<FlatIndex> = std::iter::from_fn(|| iter.next()).collect();

    // Should have 2 indices (since dimension 0 has size 2)
    assert_eq!(indices.len(), 2);
    assert_eq!(indices[0], 0.into());
    assert_eq!(indices[1], 3.into());

    // Test iteration with conjugate
    let mut iter_conj = CoreFlatFiberIterator::new(&fiber, true);
    let indices_conj: Vec<FlatIndex> = std::iter::from_fn(|| iter_conj.next()).collect();

    // Should still have the right number of indices but different pattern
    assert_eq!(indices_conj.len(), 3); // Size of dimension 1
}

#[test]
fn test_core_expanded_iterator() {
    let rep = Euclidean {};
    let structure: OrderedStructure =
        PermutedStructure::from_iter([rep.new_slot(3, 0), rep.new_slot(5, 0), rep.new_slot(2, 1)])
            .structure;

    let mut fiber = Fiber::zeros(&structure);
    fiber.free(1);
    fiber.fix(2, 1);
    fiber.fix(0, 2);

    // Create an expanded iterator
    let mut iter = CoreExpandedFiberIterator::new(&fiber, false);

    // Collect all indices
    let indices: Vec<usize> = std::iter::from_fn(|| iter.next().map(|a| a.into())).collect();

    // Should iterate through all elements (2*3 = 6)
    assert_eq!(indices, vec![0, 5, 10]);

    // Test reset functionality
    iter.reset();
    let first_idx = iter.next().unwrap();
    assert_eq!(first_idx, 0.into());

    // Test shifting
    iter.reset();
    iter.shift(10);
    let shifted_idx = iter.next().unwrap();
    assert_eq!(shifted_idx, 10.into());
}

#[test]
fn test_metric_iterator() {
    let rep = Minkowski {};
    let structure: OrderedStructure<Minkowski> =
        OrderedStructure::new(vec![rep.new_slot(2, 0), rep.new_slot(2, 1)]).structure;

    // Create a fiber with all indices free
    let mut fiber = Fiber::zeros(&structure);
    fiber.free(0);
    fiber.free(1);

    // Create a metric iterator
    let iter = MetricFiberIterator::new(&fiber, false);

    // Collect all indices and their signs
    let mut results = Vec::new();
    for MetricItem { neg, item } in iter {
        results.push((item, neg));
    }

    // Should iterate through all elements (2*2 = 4)
    assert_eq!(results.len(), 4);

    // Some elements should have negative sign due to the pseudo-Euclidean dimension
    assert!(results.iter().any(|(_, neg)| *neg));
    assert!(results.iter().any(|(_, neg)| !*neg));
}

#[test]
fn test_permuted_iterator() {
    let rep = Euclidean {};
    let structure: OrderedStructure<Euclidean> =
        OrderedStructure::new(vec![rep.new_slot(2, 0), rep.new_slot(3, 0)]).structure;

    // Create a fiber with all indices free
    let mut fiber = Fiber::zeros(&structure);
    fiber.free(0);
    fiber.free(1);

    // Create a permutation that swaps dimensions 0 and 1
    let perm = Permutation::from_map(vec![1, 0]);

    // Create an expanded iterator with permutation
    let mut iter = CoreExpandedFiberIterator::new_permuted(&fiber, false, perm);

    // Collect all indices
    let indices: Vec<FlatIndex> = std::iter::from_fn(|| iter.next()).collect();

    // Should still have all elements (2*3 = 6)
    assert_eq!(indices.len(), 6);

    // The order should be different due to permutation
    // But without direct access to the expansion, it's hard to test the exact order
    // This is just testing that iteration completes
}

#[test]
fn test_paired_conjugates() {
    let rep = Euclidean {};
    let structure: OrderedStructure<Euclidean> =
        OrderedStructure::new(vec![rep.new_slot(2, 0), rep.new_slot(3, 0)]).structure;

    // Create a fiber with specific pattern
    let mut fiber = Fiber::zeros(&structure);
    fiber.free(0);

    // Get paired conjugate iterators
    let (mut iter_conj, mut iter) = CoreFlatFiberIterator::new_paired_conjugates(&fiber);

    // The regular iterator should iterate over dimension 0 (size 2)
    let indices: Vec<FlatIndex> = std::iter::from_fn(|| iter.next()).collect();
    assert_eq!(indices.len(), 2);

    // The conjugate iterator should iterate over dimension 1 (size 3)
    let indices_conj: Vec<FlatIndex> = std::iter::from_fn(|| iter_conj.next()).collect();
    assert_eq!(indices_conj.len(), 3);
}
