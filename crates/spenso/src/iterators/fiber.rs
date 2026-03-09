//! Fiber abstractions for tensor iteration
//!
//! This module contains implementations of fiber types that represent fixed and free
//! tensor dimensions, enabling efficient iteration through tensor elements.

use bitvec::vec::BitVec;
use std::{
    fmt::{Debug, Display},
    ops::Index,
};

use crate::structure::{
    TensorStructure, concrete_index::FlatIndex, dimension::Dimension,
    representation::Representation, slot::IsAbstractSlot,
};

use super::indices::{AbstractFiberIndex, FiberClassIndex, FiberData, FiberIndex};
use super::traits::AbstractFiber;

/// Represents the bare internal structure of a fiber
///
/// This is a low-level implementation detail that tracks which indices
/// are fixed and which are free. It's used by the higher-level fiber types.
#[derive(Debug, Clone)]
struct BareFiber {
    indices: Vec<FiberIndex>,
    is_single: FiberIndex,
}

impl BareFiber {
    /// Creates a bare fiber from fiber data and a tensor structure
    pub fn from<I: TensorStructure>(data: FiberData, structure: &I) -> Self {
        match data {
            FiberData::Flat(i) => Self::from_flat(i, structure),
            FiberData::BoolFilter(b) => Self::from_filter(b),
            FiberData::BitVec(b) => Self::from_bitvec(b),
            FiberData::Single(i) => {
                let mut out = Self::zeros(structure);
                out.free(i);
                out
            }
            FiberData::IntFilter(i) => {
                let mut out = Self::zeros(structure);
                for (pos, val) in i.iter().enumerate() {
                    if *val > 0 {
                        out.free(pos);
                    }
                }
                out
            }
            FiberData::Pos(i) => {
                let mut out = Self::zeros(structure);
                for (pos, val) in i.iter().enumerate() {
                    if *val < 0 {
                        out.free(pos);
                    } else {
                        out.fix(pos, *val as usize);
                    }
                }
                out
            }
        }
    }

    /// Returns a bitvector where true represents free indices
    pub fn bitvec(&self) -> BitVec {
        self.indices.iter().map(|x| x.is_free()).collect()
    }

    /// Returns a bitvector where true represents fixed indices
    pub fn bitvecinv(&self) -> BitVec {
        self.indices.iter().map(|x| x.is_fixed()).collect()
    }

    /// Creates a bare fiber from a flat index and tensor structure
    pub fn from_flat<I>(flat: FlatIndex, structure: &I) -> BareFiber
    where
        I: TensorStructure,
    {
        let expanded = structure.expanded_index(flat).unwrap();

        BareFiber {
            indices: expanded.into_iter().map(FiberIndex::from).collect(),
            is_single: FiberIndex::Free,
        }
    }

    /// Creates a bare fiber from a boolean filter
    ///
    /// # Arguments
    ///
    /// * `filter` - A slice of booleans where true indicates a free index
    pub fn from_filter(filter: &[bool]) -> BareFiber {
        let mut f = BareFiber {
            indices: filter
                .iter()
                .map(|i| {
                    if *i {
                        FiberIndex::Free
                    } else {
                        FiberIndex::Fixed(0)
                    }
                })
                .collect(),
            is_single: FiberIndex::Free,
        };
        f.is_single();
        f
    }

    /// Creates a bare fiber from a boolean filter
    ///
    /// # Arguments
    ///
    /// * `filter` - A slice of booleans where true indicates a free index
    pub fn from_bitvec(filter: &BitVec) -> BareFiber {
        let mut f = BareFiber {
            indices: filter
                .iter()
                .map(|i| {
                    if *i {
                        FiberIndex::Free
                    } else {
                        FiberIndex::Fixed(0)
                    }
                })
                .collect(),
            is_single: FiberIndex::Free,
        };
        f.is_single();
        f
    }

    /// Creates a bare fiber with all indices fixed to zero
    pub fn zeros<I: TensorStructure>(structure: &I) -> BareFiber {
        BareFiber {
            indices: vec![FiberIndex::Fixed(0); structure.order()],
            is_single: FiberIndex::Free,
        }
    }

    /// Fixes an index to a specific value
    ///
    /// # Arguments
    ///
    /// * `pos` - The position to fix
    /// * `val` - The value to fix it to
    pub fn fix(&mut self, pos: usize, val: usize) {
        if let FiberIndex::Fixed(single_pos) = self.is_single
            && single_pos == pos
        {
            self.is_single = FiberIndex::Free;
        }

        self.indices[pos] = val.into();
    }

    /// Determines if the fiber has exactly one free index and returns its position
    pub fn is_single(&mut self) -> FiberIndex {
        if let FiberIndex::Fixed(pos) = self.is_single {
            FiberIndex::Fixed(pos)
        } else {
            let mut has_one = false;
            let mut has_two = false;
            let mut pos = 0;
            for (posi, index) in self.indices.iter().enumerate() {
                if let FiberIndex::Free = index {
                    if !has_one {
                        has_one = true;
                        pos = posi;
                    } else {
                        has_two = true;
                    }
                }
            }
            if has_one && !has_two {
                self.is_single = FiberIndex::Fixed(pos);
                return FiberIndex::Fixed(pos);
            }
            self.is_single = FiberIndex::Free;
            FiberIndex::Free
        }
    }

    /// Sets an index to be free (variable during iteration)
    ///
    /// # Arguments
    ///
    /// * `pos` - The position to make free
    pub fn free(&mut self, pos: usize) {
        self.indices[pos] = FiberIndex::Free;
    }
}

impl Index<usize> for BareFiber {
    type Output = FiberIndex;

    fn index(&self, index: usize) -> &Self::Output {
        &(self.indices[index])
    }
}

/// Represents a fiber through a tensor (immutable)
///
/// A fiber is a slice through a tensor where some indices are fixed and others vary.
/// It facilitates iterating over specific patterns within tensors.
#[derive(Debug)]
pub struct Fiber<'a, I: TensorStructure> {
    pub(crate) structure: &'a I,
    bare_fiber: BareFiber,
}

impl<I: TensorStructure> Clone for Fiber<'_, I> {
    fn clone(&self) -> Self {
        Fiber {
            structure: self.structure,
            bare_fiber: self.bare_fiber.clone(),
        }
    }
}

impl<I> Index<usize> for Fiber<'_, I>
where
    I: TensorStructure,
{
    type Output = FiberIndex;

    fn index(&self, index: usize) -> &Self::Output {
        &(self.bare_fiber[index])
    }
}

impl<I> AbstractFiber<FiberIndex> for Fiber<'_, I>
where
    I: TensorStructure,
{
    type Repr = <I::Slot as IsAbstractSlot>::R;
    fn strides(&self) -> Vec<usize> {
        self.structure.strides().unwrap()
    }

    fn reps(&self) -> Vec<Representation<Self::Repr>> {
        self.structure.reps()
    }

    fn shape(&self) -> Vec<Dimension> {
        self.structure.shape()
    }

    fn order(&self) -> usize {
        self.structure.order()
    }

    fn single(&self) -> Option<usize> {
        if let FiberIndex::Fixed(pos) = self.bare_fiber.is_single {
            Some(pos)
        } else {
            None
        }
    }

    fn bitvec(&self) -> BitVec {
        self.bare_fiber.bitvec()
    }
}

impl<'a, S> Fiber<'a, S>
where
    S: TensorStructure,
{
    /// Creates a conjugate fiber
    pub fn conj(self) -> Self {
        self
    }

    /// Creates an iterator over the fiber
    pub fn iter(
        self,
    ) -> super::fiber_iterators::FiberIterator<'a, S, super::core_iterators::CoreFlatFiberIterator>
    {
        super::fiber_iterators::FiberIterator::new(self, false)
    }

    /// Creates a conjugate iterator over the fiber
    pub fn iter_conj(
        self,
    ) -> super::fiber_iterators::FiberIterator<'a, S, super::core_iterators::CoreFlatFiberIterator>
    {
        super::fiber_iterators::FiberIterator::new(self, true)
    }

    /// Creates a permuted iterator over the fiber
    pub fn iter_perm(
        self,
        permutation: linnet::permutation::Permutation,
    ) -> super::fiber_iterators::FiberIterator<
        'a,
        S,
        super::core_iterators::CoreExpandedFiberIterator<<S::Slot as IsAbstractSlot>::R>,
    > {
        super::fiber_iterators::FiberIterator::new_permuted(self, permutation, false)
    }

    /// Creates a metric iterator over the fiber
    pub fn iter_metric(
        self,
    ) -> super::fiber_iterators::FiberIterator<
        'a,
        S,
        super::core_iterators::MetricFiberIterator<<S::Slot as IsAbstractSlot>::R>,
    > {
        super::fiber_iterators::FiberIterator::new(self, false)
    }

    /// Creates a permuted metric iterator over the fiber
    pub fn iter_perm_metric(
        self,
        permutation: linnet::permutation::Permutation,
    ) -> super::fiber_iterators::FiberIterator<
        'a,
        S,
        super::core_iterators::MetricFiberIterator<<S::Slot as IsAbstractSlot>::R>,
    > {
        super::fiber_iterators::FiberIterator::new_permuted(self, permutation, false)
    }

    /// Creates a fiber from data and a tensor structure
    pub fn from<'b>(data: impl Into<FiberData<'b>>, structure: &'a S) -> Self {
        Fiber {
            bare_fiber: BareFiber::from(data.into(), structure),
            structure,
        }
    }

    /// Returns a bitvector where true represents free indices
    pub fn bitvec(&self) -> BitVec {
        self.bare_fiber.bitvec()
    }

    /// Returns a bitvector where true represents fixed indices
    pub fn bitvecinv(&self) -> BitVec {
        self.bare_fiber.bitvecinv()
    }

    /// Creates a fiber from a flat index
    pub fn from_flat(flat: FlatIndex, structure: &'a S) -> Fiber<'a, S> {
        Fiber {
            bare_fiber: BareFiber::from_flat(flat, structure),
            structure,
        }
    }

    /// Creates a fiber from a boolean filter
    pub fn from_filter(filter: &[bool], structure: &'a S) -> Fiber<'a, S> {
        Fiber {
            bare_fiber: BareFiber::from_filter(filter),
            structure,
        }
    }

    /// Creates a fiber with all indices fixed to zero
    pub fn zeros(structure: &'a S) -> Fiber<'a, S> {
        Fiber {
            bare_fiber: BareFiber::zeros(structure),
            structure,
        }
    }

    /// Fixes an index to a specific value
    pub fn fix(&mut self, pos: usize, val: usize) {
        self.bare_fiber.fix(pos, val);
    }

    /// Determines if the fiber has exactly one free index and returns its position
    pub fn is_single(&mut self) -> FiberIndex {
        self.bare_fiber.is_single()
    }

    /// Sets an index to be free (variable during iteration)
    pub fn free(&mut self, pos: usize) {
        self.bare_fiber.free(pos);
    }
}

impl<I: TensorStructure> Display for Fiber<'_, I> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for index in self.bare_fiber.indices.iter() {
            write!(f, "{} ", index)?
        }
        Ok(())
    }
}

/// Represents a mutable fiber through a tensor
///
/// Similar to `Fiber` but allows mutation of the underlying tensor.
#[derive(Debug)]
pub struct FiberMut<'a, I: TensorStructure> {
    pub structure: &'a mut I,
    bare_fiber: BareFiber,
}

impl<I> Index<usize> for FiberMut<'_, I>
where
    I: TensorStructure,
{
    type Output = FiberIndex;

    fn index(&self, index: usize) -> &Self::Output {
        &(self.bare_fiber[index])
    }
}

impl<I> AbstractFiber<FiberIndex> for FiberMut<'_, I>
where
    I: TensorStructure,
{
    type Repr = <I::Slot as IsAbstractSlot>::R;
    fn strides(&self) -> Vec<usize> {
        self.structure.strides().unwrap()
    }

    fn reps(&self) -> Vec<Representation<Self::Repr>> {
        self.structure.reps()
    }

    fn shape(&self) -> Vec<Dimension> {
        self.structure.shape()
    }

    fn order(&self) -> usize {
        self.structure.order()
    }

    fn single(&self) -> Option<usize> {
        if let FiberIndex::Fixed(pos) = self.bare_fiber.is_single {
            Some(pos)
        } else {
            None
        }
    }

    fn bitvec(&self) -> BitVec {
        self.bare_fiber.bitvec()
    }
}

impl<'a, I> FiberMut<'a, I>
where
    I: TensorStructure,
{
    /// Creates a fiber from data and a mutable tensor structure
    pub fn from<'b>(data: FiberData<'b>, structure: &'a mut I) -> Self {
        FiberMut {
            bare_fiber: BareFiber::from(data, &*structure),
            structure,
        }
    }

    /// Creates a conjugate fiber
    pub fn conj(self) -> Self {
        self
    }

    /// Returns a bitvector where true represents free indices
    pub fn bitvec(&self) -> BitVec {
        self.bare_fiber.bitvec()
    }

    /// Returns a bitvector where true represents fixed indices
    pub fn bitvecinv(&self) -> BitVec {
        self.bare_fiber.bitvecinv()
    }

    /// Creates a fiber with all indices fixed to zero
    pub fn zeros(structure: &'a I) -> Fiber<'a, I> {
        Fiber {
            bare_fiber: BareFiber::zeros(structure),
            structure,
        }
    }

    /// Fixes an index to a specific value
    pub fn fix(&mut self, pos: usize, val: usize) {
        self.bare_fiber.fix(pos, val);
    }

    /// Determines if the fiber has exactly one free index and returns its position
    pub fn is_single(&mut self) -> FiberIndex {
        self.bare_fiber.is_single()
    }

    /// Sets an index to be free (variable during iteration)
    pub fn free(&mut self, pos: usize) {
        self.bare_fiber.free(pos);
    }
}

impl<I: TensorStructure> Display for FiberMut<'_, I> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for index in self.bare_fiber.indices.iter() {
            write!(f, "{} ", index)?
        }
        Ok(())
    }
}

impl<'a, I: TensorStructure> FiberMut<'a, I> {
    /// Creates an iterator over the mutable fiber
    pub fn iter(
        self,
    ) -> super::fiber_iterators::MutFiberIterator<'a, I, super::core_iterators::CoreFlatFiberIterator>
    {
        super::fiber_iterators::MutFiberIterator::new(self, false)
    }
}

/// Represents a class of fibers with similar structure
///
/// Fiber classes group fibers that share a pattern of fixed and free indices,
/// allowing iteration over all fibers of that class.
pub struct FiberClass<'a, I: TensorStructure> {
    structure: &'a I,
    bare_fiber: BareFiber, // A representant of the class
}

impl<I: TensorStructure> Clone for FiberClass<'_, I> {
    fn clone(&self) -> Self {
        FiberClass {
            bare_fiber: self.bare_fiber.clone(),
            structure: self.structure,
        }
    }
}

impl<I> Index<usize> for FiberClass<'_, I>
where
    I: TensorStructure,
{
    type Output = FiberClassIndex;

    fn index(&self, index: usize) -> &Self::Output {
        if self.bare_fiber[index].is_fixed() {
            &FiberClassIndex::Free
        } else {
            &FiberClassIndex::Fixed
        }
    }
}

impl<'a, I: TensorStructure> From<Fiber<'a, I>> for FiberClass<'a, I> {
    fn from(fiber: Fiber<'a, I>) -> Self {
        FiberClass {
            bare_fiber: fiber.bare_fiber,
            structure: fiber.structure,
        }
    }
}

impl<'a, I: TensorStructure> From<FiberClass<'a, I>> for Fiber<'a, I> {
    fn from(fiber: FiberClass<'a, I>) -> Self {
        Fiber {
            bare_fiber: fiber.bare_fiber,
            structure: fiber.structure,
        }
    }
}

impl<I: TensorStructure> AbstractFiber<FiberClassIndex> for FiberClass<'_, I> {
    type Repr = <I::Slot as IsAbstractSlot>::R;
    fn strides(&self) -> Vec<usize> {
        self.structure.strides().unwrap()
    }

    fn shape(&self) -> Vec<Dimension> {
        self.structure.shape()
    }

    fn reps(&self) -> Vec<Representation<Self::Repr>> {
        self.structure.reps()
    }

    fn order(&self) -> usize {
        self.structure.order()
    }

    fn single(&self) -> Option<usize> {
        match self.bare_fiber.is_single {
            FiberIndex::Fixed(i) => Some(i),
            _ => None,
        }
    }

    fn bitvec(&self) -> BitVec {
        !self.bare_fiber.bitvec()
    }
}

impl<'a, S: TensorStructure> FiberClass<'a, S> {
    /// Creates an iterator over the fiber class
    pub fn iter(self) -> super::fiber_iterators::FiberClassIterator<'a, S> {
        super::fiber_iterators::FiberClassIterator::new(self)
    }

    /// Creates a permuted iterator over the fiber class
    pub fn iter_perm(
        self,
        permutation: linnet::permutation::Permutation,
    ) -> super::fiber_iterators::FiberClassIterator<
        'a,
        S,
        super::core_iterators::CoreExpandedFiberIterator<<S::Slot as IsAbstractSlot>::R>,
    > {
        super::fiber_iterators::FiberClassIterator::new_permuted(self, permutation)
    }

    /// Creates a permuted metric iterator over the fiber class
    pub fn iter_perm_metric(
        self,
        permutation: linnet::permutation::Permutation,
    ) -> super::fiber_iterators::FiberClassIterator<
        'a,
        S,
        super::core_iterators::MetricFiberIterator<<S::Slot as IsAbstractSlot>::R>,
    > {
        super::fiber_iterators::FiberClassIterator::new_permuted(self, permutation)
    }
}

/// Represents a mutable class of fibers with similar structure
///
/// Similar to `FiberClass` but allows mutation of the underlying tensor.
pub struct FiberClassMut<'a, I: TensorStructure> {
    structure: &'a mut I,
    bare_fiber: BareFiber, // A representant of the class
}

impl<I> Index<usize> for FiberClassMut<'_, I>
where
    I: TensorStructure,
{
    type Output = FiberClassIndex;

    fn index(&self, index: usize) -> &Self::Output {
        if self.bare_fiber[index].is_fixed() {
            &FiberClassIndex::Free
        } else {
            &FiberClassIndex::Fixed
        }
    }
}

impl<'a, I: TensorStructure> From<FiberMut<'a, I>> for FiberClassMut<'a, I> {
    fn from(fiber: FiberMut<'a, I>) -> Self {
        FiberClassMut {
            bare_fiber: fiber.bare_fiber,
            structure: fiber.structure,
        }
    }
}

impl<'a, I: TensorStructure> From<FiberClassMut<'a, I>> for FiberMut<'a, I> {
    fn from(fiber: FiberClassMut<'a, I>) -> Self {
        FiberMut {
            bare_fiber: fiber.bare_fiber,
            structure: fiber.structure,
        }
    }
}

impl<I: TensorStructure> AbstractFiber<FiberClassIndex> for FiberClassMut<'_, I> {
    type Repr = <I::Slot as IsAbstractSlot>::R;
    fn strides(&self) -> Vec<usize> {
        self.structure.strides().unwrap()
    }

    fn shape(&self) -> Vec<Dimension> {
        self.structure.shape()
    }

    fn reps(&self) -> Vec<Representation<Self::Repr>> {
        self.structure.reps()
    }

    fn order(&self) -> usize {
        self.structure.order()
    }

    fn single(&self) -> Option<usize> {
        match self.bare_fiber.is_single {
            FiberIndex::Fixed(i) => Some(i),
            _ => None,
        }
    }

    fn bitvec(&self) -> BitVec {
        !self.bare_fiber.bitvec()
    }
}
