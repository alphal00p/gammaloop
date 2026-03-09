//! Tensor index and fiber abstractions
//!
//! This module provides the core abstractions for tensor indices and fibers,
//! which are used to navigate tensor data structures efficiently.

use crate::structure::concrete_index::FlatIndex;
use bitvec::vec::BitVec;
use serde::{Deserialize, Serialize};
use std::fmt::{Debug, Display};

/// Abstract trait for fiber indices
///
/// This trait defines methods for determining whether an index is free or fixed,
/// which is fundamental to how fibers are defined and traversed.
pub trait AbstractFiberIndex {
    /// Returns true if the index is free (varying during iteration)
    fn is_free(&self) -> bool;

    /// Returns true if the index is fixed (constant during iteration)
    fn is_fixed(&self) -> bool {
        !self.is_free()
    }
}

/// Represents a classification of fiber indices (free or fixed)
///
/// Used for fiber class iterators to categorize indices.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum FiberClassIndex {
    /// A free index that varies during iteration
    Free,
    /// A fixed index that remains constant during iteration
    Fixed,
}

impl AbstractFiberIndex for FiberClassIndex {
    fn is_free(&self) -> bool {
        matches!(self, FiberClassIndex::Free)
    }
}

/// Represents a specific fiber index, either free or fixed to a particular value
///
/// Used for precise control of which dimensions are fixed and which vary during iteration.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum FiberIndex {
    /// A free index that varies during iteration
    Free,
    /// A fixed index with a specific value
    Fixed(usize),
}

impl AbstractFiberIndex for FiberIndex {
    fn is_free(&self) -> bool {
        matches!(self, FiberIndex::Free)
    }
}

impl From<usize> for FiberIndex {
    fn from(value: usize) -> Self {
        FiberIndex::Fixed(value)
    }
}

impl Display for FiberIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            FiberIndex::Fixed(val) => {
                write!(f, "{val}")
            }
            FiberIndex::Free => {
                write!(f, ":")
            }
        }
    }
}

/// Enumeration of different ways to specify fiber data
///
/// This is used as a flexible input mechanism to create fibers from various data sources.
pub enum FiberData<'a> {
    /// A single index position
    Single(usize),
    /// A flat index directly
    Flat(FlatIndex),
    /// A boolean filter where true indicates free indices
    BoolFilter(&'a [bool]),
    /// Position data where negative values mean free indices
    Pos(&'a [isize]),
    /// Integer filter where non-zero values indicate free indices
    IntFilter(&'a [u8]),
    /// A bitvec filter where true indicates free indices
    BitVec(&'a BitVec),
}

impl From<usize> for FiberData<'_> {
    fn from(value: usize) -> Self {
        Self::Single(value)
    }
}

impl<'a> From<&'a [bool]> for FiberData<'a> {
    fn from(value: &'a [bool]) -> Self {
        Self::BoolFilter(value)
    }
}

impl<'a> From<&'a BitVec> for FiberData<'a> {
    fn from(value: &'a BitVec) -> Self {
        Self::BitVec(value)
    }
}

impl<'a> From<&'a [u8]> for FiberData<'a> {
    fn from(value: &'a [u8]) -> Self {
        Self::IntFilter(value)
    }
}

impl<'a> From<&'a [isize]> for FiberData<'a> {
    fn from(value: &'a [isize]) -> Self {
        Self::Pos(value)
    }
}

impl From<FlatIndex> for FiberData<'_> {
    fn from(value: FlatIndex) -> Self {
        Self::Flat(value)
    }
}

/// Enumeration for switching between iterator implementations
pub enum IteratorEnum<A, B> {
    /// First iterator type
    A(A),
    /// Second iterator type
    B(B),
}

impl<A, B> Iterator for IteratorEnum<A, B>
where
    A: Iterator,
    B: Iterator<Item = A::Item>,
{
    type Item = A::Item;
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            IteratorEnum::A(a) => a.next(),
            IteratorEnum::B(b) => b.next(),
        }
    }
}
