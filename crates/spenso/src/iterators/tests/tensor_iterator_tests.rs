//! Tests for tensor-specific iterators
//!
//! Contains tests for iterators designed for different tensor types.

use crate::{
    iterators::{
        DenseTensorIterator, DenseTensorLinearIterator, FiberData, IteratableTensor,
        SparseTensorIterator, SparseTensorLinearIterator, TensorStructureIndexIterator,
    },
    structure::{
        OrderedStructure, PermutedStructure,
        concrete_index::{ExpandedIndex, FlatIndex},
        representation::{Euclidean, RepName},
    },
    tensors::data::{DenseTensor, SetTensorData, SparseTensor},
};
use std::collections::HashSet;

#[test]
fn test_dense_tensor_iterators() {
    let rep = Euclidean {};
    let structure: OrderedStructure<Euclidean> =
        OrderedStructure::new(vec![rep.new_slot(2, 0), rep.new_slot(3, 0)]).structure;

    // Create a dense tensor with values 0...5
    let mut tensor = DenseTensor::<i32, _>::zero(structure);
    for i in 0..6 {
        tensor.data[i] = i as i32;
    }

    // Test expanded iterator
    let expanded_iter = DenseTensorIterator::new(&tensor);
    let expanded_results: Vec<(ExpandedIndex, &i32)> = expanded_iter.collect();

    assert_eq!(expanded_results.len(), 6); // 2x3 tensor

    // Test linear iterator
    let linear_iter = DenseTensorLinearIterator::new(&tensor);
    let linear_results: Vec<(FlatIndex, &i32)> = linear_iter.collect();

    assert_eq!(linear_results.len(), 6);
    for (i, &(_, &val)) in linear_results.iter().enumerate() {
        assert_eq!(val, i as i32);
    }

    // Test IntoIterator
    let into_iter_results: Vec<(ExpandedIndex, &i32)> = (&tensor).into_iter().collect();
    assert_eq!(into_iter_results.len(), 6);

    // Test IteratableTensor trait
    let fiber_data = [true, false].as_slice().into(); // Iterate along first dimension
    let fiber = tensor.fiber(fiber_data);
    let fiber_iter = fiber.iter();
    let fiber_results: Vec<(&i32, ())> = fiber_iter.collect();

    assert_eq!(fiber_results.len(), 2); // Size of dimension 0
}

#[test]
fn test_sparse_tensor_iterators() {
    let rep = Euclidean {};
    let structure: OrderedStructure<Euclidean> =
        OrderedStructure::new(vec![rep.new_slot(2, 0), rep.new_slot(3, 0)]).structure;

    // Create a sparse tensor with a few values
    let mut tensor = SparseTensor::<i32, _>::empty(structure, 0);

    // Set a few values
    tensor.set(&[0, 0], 1).unwrap();
    tensor.set(&[1, 1], 2).unwrap();
    tensor.set(&[1, 2], 3).unwrap();

    // Test linear iterator
    let linear_iter = SparseTensorLinearIterator::new(&tensor);
    let linear_results: Vec<(FlatIndex, &i32)> = linear_iter.collect();

    assert_eq!(linear_results.len(), 3); // Only 3 values set

    // Test expanded iterator
    let expanded_iter = SparseTensorIterator::new(&tensor);
    let expanded_results: Vec<(ExpandedIndex, &i32)> = expanded_iter.collect();

    assert_eq!(expanded_results.len(), 3);

    // Verify each index has the expected value
    for &(ref indices, &val) in &expanded_results {
        match indices.indices.as_slice() {
            [0, 0] => assert_eq!(val, 1),
            [1, 1] => assert_eq!(val, 2),
            [1, 2] => assert_eq!(val, 3),
            _ => panic!("Unexpected index"),
        }
    }
}

#[test]
fn test_fiber_iterators() {
    let rep = Euclidean {};
    let structure: PermutedStructure<OrderedStructure<Euclidean>> =
        OrderedStructure::new(vec![rep.new_slot(2, 0), rep.new_slot(3, 0)]);

    let mut tensor = DenseTensor::<i32, _>::zero((structure.structure).clone());
    for i in 0..6 {
        tensor.data[i] = i as i32;
    }

    // Create a fiber that iterates over the first dimension
    let fiber = tensor.fiber(FiberData::BoolFilter(&[true, false]));

    // Use FiberIterator to traverse the tensor
    let fiber_iter = fiber.iter();
    let results: Vec<(&i32, ())> = fiber_iter.collect();

    assert_eq!(results.len(), 2);
    assert_eq!(*results[0].0, 0); // First element in dimension 0
    assert_eq!(*results[1].0, 3); // Second element in dimension 0, fixing dimension 1 to 0
}

#[test]
fn test_fiber_class_iterators() {
    let rep = Euclidean {};
    let structure: PermutedStructure<OrderedStructure<Euclidean>> =
        OrderedStructure::new(vec![rep.new_slot(2, 0), rep.new_slot(3, 0)]);

    let mut tensor = DenseTensor::<i32, _>::zero((structure.structure).clone());
    for i in 0..6 {
        tensor.data[i] = i as i32;
    }

    // Create a fiber class for the second dimension
    let fiber_class = tensor.fiber_class(FiberData::BoolFilter(&[false, true]));

    // Create class iterator
    let mut class_iter = fiber_class.iter();

    // The class iterator should give us 2 fibers (one for each value in dimension 0)
    let fiber1_opt = class_iter.next();
    assert!(fiber1_opt.is_some());

    let fiber2_opt = class_iter.next();
    assert!(fiber2_opt.is_some());

    let fiber3_opt = class_iter.next();
    assert!(fiber3_opt.is_none()); // No more fibers

    // Now test that we can iterate through each fiber
    let fiber1 = fiber1_opt.unwrap();
    let results1: Vec<(&i32, ())> = fiber1.collect();
    assert_eq!(results1.len(), 3); // Dimension 1 has size 3
}

#[test]
fn test_tensor_structure_index_iterator() {
    let rep = Euclidean {};
    let structure: OrderedStructure<Euclidean> =
        OrderedStructure::new(vec![rep.new_slot(2, 0), rep.new_slot(3, 0)]).structure;

    // Create an index iterator
    let index_iter = TensorStructureIndexIterator::new(&structure);
    let indices: Vec<ExpandedIndex> = index_iter.collect();

    // Should have one index for each element in the tensor
    assert_eq!(indices.len(), 6);

    // Check that all indices are unique
    let unique_indices: HashSet<_> = indices.iter().collect();
    assert_eq!(unique_indices.len(), 6);

    // Verify the bounds of the indices
    for idx in indices {
        assert!(idx[0] < 2);
        assert!(idx[1] < 3);
    }
}
