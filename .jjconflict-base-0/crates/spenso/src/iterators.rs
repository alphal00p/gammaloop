//! Iterators for tensors
//!
//! This module provides various iterators for tensor operations, organized as:
//!
//! - `indices`: Types and traits related to tensor indices and fibers
//! - `fiber`: Fiber abstractions for traversing tensor data
//! - `core_iterators`: Low-level iterators implementing core iteration logic
//! - `fiber_iterators`: Iterators for fibers (fixed index or class-based)
//! - `tensor_iterators`: Iterators specific to different tensor types
//! - `traits`: Common traits defining iterator behaviors
//!
//! Iterators for tensors are used to iterate over the elements of a tensor.
//! More specialized iterators are provided that fix a certain subset of indices, and iterate over the remaining indices.
//! At each iteration, the iterator returns a vector of references to the elements of the tensor along the fixed indices (so called fibers).
//!
//! The iterators are built using the basic index iterators provided by the `TensorStructure+ScalarStructureIterator`s.

mod core_iterators;
mod fiber;
mod fiber_iterators;
mod indices;
mod tensor_iterators;
mod traits;

#[cfg(test)]
mod tests;

// Re-export common types and traits
pub use core_iterators::{
    CoreExpandedFiberIterator, CoreFlatFiberIterator, MetricFiberIterator, MetricItem,
};
pub use fiber::{Fiber, FiberClass, FiberClassMut, FiberMut};
pub use fiber_iterators::{FiberClassIterator, FiberIterator, MutFiberIterator};
pub use indices::{AbstractFiberIndex, FiberClassIndex, FiberData, FiberIndex, IteratorEnum};
pub use tensor_iterators::{
    DataTensorExpandedIterator, DataTensorLinearIterator, DenseTensorIntoIterator,
    DenseTensorIterator, DenseTensorLinearIterator, DenseTensorTraceIterator, SparseTensorIterator,
    SparseTensorLinearIterator, SparseTensorTraceIterator, TensorStructureIndexIterator,
};
pub use traits::{
    AbstractFiber, FiberIteratorItem, IteratableTensor, IteratesAlongFibers,
    IteratesAlongPermutedFibers, ResetableIterator, ShiftableIterator,
};

// Re-export helper types
pub use core_iterators::{MultiStrideShift, SingleStrideShift, SkippingItem, StrideShift};
