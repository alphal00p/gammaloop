//! Traits defining iterator behaviors for tensors
//!
//! This module provides the fundamental traits that define how iterators
//! operate on tensor data structures.

use super::fiber::{Fiber, FiberClass, FiberClassMut, FiberMut};
use super::indices::{AbstractFiberIndex, FiberData};
use crate::structure::{
    concrete_index::{ExpandedIndex, FlatIndex},
    dimension::Dimension,
    representation::{RepName, Representation},
    HasStructure, TensorStructure,
};
use bitvec::vec::BitVec;
use linnet::permutation::Permutation;

/// Trait for items yielded by fiber iterators
///
/// Defines the interface for items returned by fiber iterators, including
/// their flat index representation and any additional data they carry.
pub trait FiberIteratorItem {
    /// Additional data carried by the iterator item
    type OtherData;

    /// Returns the flat index of the item
    fn flat_idx(&self) -> FlatIndex;

    /// Extracts additional data from the item
    fn other_data(self) -> Self::OtherData;
}

impl FiberIteratorItem for FlatIndex {
    type OtherData = ();
    fn flat_idx(&self) -> FlatIndex {
        *self
    }

    fn other_data(self) -> Self::OtherData {}
}

pub trait ResetableIterator {
    /// Reset the iterator to its initial position
    fn reset(&mut self);
}

pub trait ShiftableIterator {
    /// Shift the iterator by the specified amount
    ///
    /// # Arguments
    ///
    /// * `shift` - The amount to shift by
    fn shift(&mut self, shift: usize);
}

/// Trait for iterators that traverse along fibers
///
/// Defines the core operations for iterators that move along tensor fibers,
/// with support for resetting and creating paired conjugate iterators.
pub trait IteratesAlongFibers<R: RepName>:
    Iterator + ShiftableIterator + ResetableIterator
{
    /// Create a new iterator from a fiber
    ///
    /// # Arguments
    ///
    /// * `fiber` - The fiber to iterate along
    /// * `conj` - Whether to use conjugate iteration
    fn new<I, J>(fiber: &I, conj: bool) -> Self
    where
        I: AbstractFiber<J, Repr = R>,
        J: AbstractFiberIndex,
        Self: Sized;

    /// Create a pair of conjugate iterators from a fiber
    ///
    /// # Arguments
    ///
    /// * `fiber` - The fiber to iterate along
    fn new_paired_conjugates<I, J>(fiber: &I) -> (Self, Self)
    where
        I: AbstractFiber<J, Repr = R>,
        J: AbstractFiberIndex,
        Self: Sized;
}

/// Trait for iterators that traverse permuted fibers
///
/// Extends `IteratesAlongFibers` with support for permutations.
pub trait IteratesAlongPermutedFibers<R: RepName>: IteratesAlongFibers<R> {
    /// Create a new iterator from a fiber with a permutation
    ///
    /// # Arguments
    ///
    /// * `fiber` - The fiber to iterate along
    /// * `conj` - Whether to use conjugate iteration
    /// * `permutation` - The permutation to apply
    fn new_permuted<I, J>(fiber: &I, conj: bool, permutation: Permutation) -> Self
    where
        I: AbstractFiber<J, Repr = R>,
        J: AbstractFiberIndex;
}

/// Abstract trait defining fiber properties
///
/// This trait defines the properties of a fiber, which is a slice through a tensor
/// along specific dimensions.
pub trait AbstractFiber<Out: AbstractFiberIndex>: std::ops::Index<usize, Output = Out> {
    /// Representation type for the fiber
    type Repr: RepName;

    /// Returns the strides of the fiber
    fn strides(&self) -> Vec<usize>;

    /// Returns the shape of the fiber
    fn shape(&self) -> Vec<Dimension>;

    /// Returns the representations of the fiber dimensions
    fn reps(&self) -> Vec<Representation<Self::Repr>>;

    /// Returns the order (number of dimensions) of the fiber
    fn order(&self) -> usize;

    /// Returns the single varying dimension if there is one
    fn single(&self) -> Option<usize>;

    /// Returns a bitvector indicating free (true) and fixed (false) dimensions
    fn bitvec(&self) -> BitVec;
}

/// Trait for tensor types that can be iterated
///
/// Defines the interface for tensor types that support various forms of iteration.
pub trait IteratableTensor: HasStructure + Sized + TensorStructure {
    /// Type of data yielded by iterators over this tensor
    type Data<'a>
    where
        Self: 'a;

    /// Returns an iterator over expanded indices and data
    fn iter_expanded(&self) -> impl Iterator<Item = (ExpandedIndex, Self::Data<'_>)>;

    /// Returns an iterator over flat indices and data
    fn iter_flat(&self) -> impl Iterator<Item = (FlatIndex, Self::Data<'_>)>;

    /// Creates a fiber from the tensor
    ///
    /// # Arguments
    ///
    /// * `fiber_data` - Data specifying which dimensions are fixed/free
    fn fiber<'a>(&'a self, fiber_data: FiberData<'_>) -> Fiber<'a, Self> {
        Fiber::from(fiber_data, self)
    }

    /// Creates a mutable fiber from the tensor
    ///
    /// # Arguments
    ///
    /// * `fiber_data` - Data specifying which dimensions are fixed/free
    fn fiber_mut<'a>(&'a mut self, fiber_data: FiberData<'_>) -> FiberMut<'a, Self> {
        FiberMut::from(fiber_data, self)
    }

    /// Creates a fiber class from the tensor
    ///
    /// # Arguments
    ///
    /// * `fiber_data` - Data specifying which dimensions are fixed/free
    fn fiber_class<'a>(&'a self, fiber_data: FiberData<'_>) -> FiberClass<'a, Self> {
        Fiber::from(fiber_data, self).into()
    }

    /// Creates a mutable fiber class from the tensor
    ///
    /// # Arguments
    ///
    /// * `fiber_data` - Data specifying which dimensions are fixed/free
    fn fiber_class_mut<'a>(&'a mut self, fiber_data: FiberData<'_>) -> FiberClassMut<'a, Self> {
        FiberMut::from(fiber_data, self).into()
    }
}
