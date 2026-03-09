//! High-level fiber iterators for tensor traversal
//!
//! This module contains iterators that build on the core iterators to provide
//! high-level fiber iteration capabilities for various tensor types.

use gat_lending_iterator::LendingIterator;
use linnet::permutation::Permutation;
use std::fmt::Debug;

use crate::{
    structure::{representation::LibraryRep, slot::IsAbstractSlot, TensorStructure},
    tensors::data::{DenseTensor, GetTensorData, SparseTensor},
};

use super::{
    core_iterators::CoreFlatFiberIterator,
    fiber::{Fiber, FiberClass, FiberMut},
    traits::ResetableIterator,
    FiberIteratorItem, IteratesAlongFibers, IteratesAlongPermutedFibers,
};

/// Iterator for traversing tensor fibers
///
/// This high-level iterator uses a core iterator to traverse fibers in tensors,
/// returning references to tensor elements at each position.
#[derive(Debug)]
pub struct FiberIterator<
    'a,
    S: TensorStructure,
    I: IteratesAlongFibers<<S::Slot as IsAbstractSlot>::R>,
> {
    /// The fiber being iterated
    pub fiber: Fiber<'a, S>,
    /// The underlying core iterator
    pub iter: I,
    /// Number of indices skipped
    pub skipped: usize,
}

impl<S: TensorStructure, I: IteratesAlongFibers<<S::Slot as IsAbstractSlot>::R> + Clone> Clone
    for FiberIterator<'_, S, I>
{
    fn clone(&self) -> Self {
        FiberIterator {
            fiber: self.fiber.clone(),
            iter: self.iter.clone(),
            skipped: self.skipped,
        }
    }
}

impl<'a, S: TensorStructure, I: IteratesAlongFibers<<S::Slot as IsAbstractSlot>::R>>
    FiberIterator<'a, S, I>
{
    /// Creates a new fiber iterator
    ///
    /// # Arguments
    ///
    /// * `fiber` - The fiber to iterate over
    /// * `conj` - Whether to use conjugate iteration
    pub fn new(fiber: Fiber<'a, S>, conj: bool) -> Self {
        FiberIterator {
            iter: I::new(&fiber, conj),
            fiber,
            skipped: 0,
        }
    }

    /// Resets the iterator to its initial state
    pub fn reset(&mut self) {
        self.iter.reset();
        self.skipped = 0;
    }

    /// Shifts the iterator by the given amount
    ///
    /// # Arguments
    ///
    /// * `shift` - The amount to shift by
    pub fn shift(&mut self, shift: usize) {
        self.reset();
        self.iter.shift(shift);
    }
}

impl<'a, S: TensorStructure, I: IteratesAlongPermutedFibers<<S::Slot as IsAbstractSlot>::R>>
    FiberIterator<'a, S, I>
{
    /// Creates a new fiber iterator with a permutation
    ///
    /// # Arguments
    ///
    /// * `fiber` - The fiber to iterate over
    /// * `permutation` - The permutation to apply
    /// * `conj` - Whether to use conjugate iteration
    pub fn new_permuted(fiber: Fiber<'a, S>, permutation: Permutation, conj: bool) -> Self {
        FiberIterator {
            iter: I::new_permuted(&fiber, conj, permutation),
            fiber,
            skipped: 0,
        }
    }
}

impl<I: IteratesAlongFibers<LibraryRep>> Iterator
    for FiberIterator<'_, crate::structure::OrderedStructure, I>
{
    type Item = I::Item;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next()
    }
}

impl<
        'a,
        I: IteratesAlongFibers<<S::Slot as IsAbstractSlot>::R, Item = It>,
        S: TensorStructure,
        T,
        It,
    > Iterator for FiberIterator<'a, DenseTensor<T, S>, I>
where
    It: FiberIteratorItem,
{
    type Item = (&'a T, It::OtherData);
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|x| {
            // println!(
            //     "DenseTensor: flat_idx: {}, size: {:?}",
            //     x.flat_idx(),
            //     self.fiber.structure.size()
            // );
            if let Some(t) = self.fiber.structure.get_ref_linear(x.flat_idx()) {
                (t, x.other_data())
            } else {
                panic!(
                    "DenseTensor: Out of bounds {} {}",
                    x.flat_idx(),
                    self.fiber.structure.size().unwrap()
                )
            }
        })
    }
}

impl<
        'a,
        I: IteratesAlongFibers<<S::Slot as IsAbstractSlot>::R, Item = It>,
        S: TensorStructure,
        T,
        It,
    > Iterator for FiberIterator<'a, SparseTensor<T, S>, I>
where
    It: FiberIteratorItem,
{
    type Item = (&'a T, usize, It::OtherData);
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(i) = self.iter.next() {
            if let Some(t) = self.fiber.structure.get_ref_linear(i.flat_idx()) {
                let skipped = self.skipped;
                self.skipped = 0;
                return Some((t, skipped, i.other_data()));
            } else {
                self.skipped += 1;
                return self.next();
            }
        }
        None
    }
}

/// Mutable iterator for traversing tensor fibers
///
/// Similar to `FiberIterator` but returns mutable references to tensor elements.
pub struct MutFiberIterator<
    'a,
    S: TensorStructure,
    I: IteratesAlongFibers<<S::Slot as IsAbstractSlot>::R>,
> {
    /// The underlying core iterator
    iter: I,
    /// The fiber being iterated
    fiber: FiberMut<'a, S>,
    /// Number of indices skipped
    skipped: usize,
}

impl<
        I: IteratesAlongFibers<<S::Slot as IsAbstractSlot>::R, Item = It>,
        S: TensorStructure,
        T,
        It,
    > LendingIterator for MutFiberIterator<'_, SparseTensor<T, S>, I>
where
    It: FiberIteratorItem,
{
    type Item<'r>
        = (&'r mut T, usize, It::OtherData)
    where
        Self: 'r;
    fn next(&mut self) -> Option<Self::Item<'_>> {
        let flat = self.iter.next()?;
        if self.fiber.structure.is_empty_at_flat(flat.flat_idx()) {
            let skipped = self.skipped;
            self.skipped = 0;
            Some((
                self.fiber
                    .structure
                    .get_mut_linear(flat.flat_idx())
                    .unwrap(),
                skipped,
                flat.other_data(),
            ))
        } else {
            self.skipped += 1;
            self.next()
        }
    }
}

impl<
        I: IteratesAlongFibers<<S::Slot as IsAbstractSlot>::R, Item = It>,
        S: TensorStructure,
        T,
        It,
    > LendingIterator for MutFiberIterator<'_, DenseTensor<T, S>, I>
where
    It: FiberIteratorItem,
{
    type Item<'r>
        = (&'r mut T, It::OtherData)
    where
        Self: 'r;
    fn next(&mut self) -> Option<Self::Item<'_>> {
        self.iter.next().map(|x| {
            (
                self.fiber.structure.get_mut_linear(x.flat_idx()).unwrap(),
                x.other_data(),
            )
        })
    }
}

impl<'a, S: TensorStructure, I: IteratesAlongFibers<<S::Slot as IsAbstractSlot>::R>>
    MutFiberIterator<'a, S, I>
{
    /// Creates a new mutable fiber iterator
    ///
    /// # Arguments
    ///
    /// * `fiber` - The fiber to iterate over
    /// * `conj` - Whether to use conjugate iteration
    pub fn new(fiber: FiberMut<'a, S>, conj: bool) -> Self {
        MutFiberIterator {
            iter: I::new(&fiber, conj),
            fiber,
            skipped: 0,
        }
    }

    /// Resets the iterator to its initial state
    pub fn reset(&mut self) {
        self.iter.reset();
        self.skipped = 0;
    }

    /// Shifts the iterator by the given amount
    ///
    /// # Arguments
    ///
    /// * `shift` - The amount to shift by
    pub fn shift(&mut self, shift: usize) {
        self.iter.shift(shift);
    }
}

impl<'a, S: TensorStructure, I: IteratesAlongPermutedFibers<<S::Slot as IsAbstractSlot>::R>>
    MutFiberIterator<'a, S, I>
{
    /// Creates a new mutable fiber iterator with a permutation
    ///
    /// # Arguments
    ///
    /// * `fiber` - The fiber to iterate over
    /// * `permutation` - The permutation to apply
    /// * `conj` - Whether to use conjugate iteration
    pub fn new_permuted(fiber: FiberMut<'a, S>, permutation: Permutation, conj: bool) -> Self {
        MutFiberIterator {
            iter: I::new_permuted(&fiber, conj, permutation),
            fiber,
            skipped: 0,
        }
    }
}

/// Iterator for traversing fiber classes
///
/// Iterates over all fibers in a fiber class, returning an iterator for each fiber.
pub struct FiberClassIterator<
    'b,
    S: TensorStructure,
    I: IteratesAlongFibers<<S::Slot as IsAbstractSlot>::R> = CoreFlatFiberIterator,
> {
    /// Iterator over fibers within the class
    pub fiber_iter: FiberIterator<'b, S, I>,
    /// Iterator over indices of fibers in the class
    pub class_iter: CoreFlatFiberIterator,
}

impl<'b, N: TensorStructure> FiberClassIterator<'b, N, CoreFlatFiberIterator> {
    /// Creates a new fiber class iterator
    ///
    /// # Arguments
    ///
    /// * `class` - The fiber class to iterate over
    pub fn new(class: FiberClass<'b, N>) -> Self {
        let (iter, iter_conj) = CoreFlatFiberIterator::new_paired_conjugates(&class);

        let fiber = FiberIterator {
            fiber: class.into(),
            iter,
            skipped: 0,
        };

        FiberClassIterator {
            fiber_iter: fiber,
            class_iter: iter_conj,
        }
    }
}

impl<N: TensorStructure, I: IteratesAlongFibers<<N::Slot as IsAbstractSlot>::R>>
    FiberClassIterator<'_, N, I>
{
    /// Resets the iterator to its initial state
    pub fn reset(&mut self) {
        self.class_iter.reset();
        self.fiber_iter.reset();
        self.fiber_iter.shift(0);
    }
}

impl<'b, N: TensorStructure, I: IteratesAlongPermutedFibers<<N::Slot as IsAbstractSlot>::R>>
    FiberClassIterator<'b, N, I>
{
    /// Creates a new fiber class iterator with a permutation
    ///
    /// # Arguments
    ///
    /// * `class` - The fiber class to iterate over
    /// * `permutation` - The permutation to apply
    pub fn new_permuted(class: FiberClass<'b, N>, permutation: Permutation) -> Self {
        let iter = CoreFlatFiberIterator::new(&class, false);

        let fiber = FiberIterator::new_permuted(class.into(), permutation, false);

        FiberClassIterator {
            fiber_iter: fiber,
            class_iter: iter,
        }
    }
}

impl<
        'a,
        S: TensorStructure + 'a,
        I: IteratesAlongFibers<<S::Slot as IsAbstractSlot>::R> + Clone + Debug,
    > Iterator for FiberClassIterator<'a, S, I>
{
    type Item = FiberIterator<'a, S, I>;

    fn next(&mut self) -> Option<Self::Item> {
        let shift = self.class_iter.next()?;
        self.fiber_iter.reset();
        self.fiber_iter.shift(shift.into());
        Some(self.fiber_iter.clone())
    }
}

impl<'a, S: TensorStructure + 'a, I: IteratesAlongFibers<<S::Slot as IsAbstractSlot>::R>>
    LendingIterator for FiberClassIterator<'a, S, I>
{
    type Item<'r>
        = &'r mut FiberIterator<'a, S, I>
    where
        Self: 'r;

    fn next(&mut self) -> Option<Self::Item<'_>> {
        let shift = self.class_iter.next()?;
        self.fiber_iter.reset();
        self.fiber_iter.shift(shift.into());
        Some(&mut self.fiber_iter)
    }
}

#[cfg(test)]
mod tests {

    use crate::structure::{
        representation::{Euclidean, RepName},
        OrderedStructure, PermutedStructure,
    };

    use super::*;

    #[test]
    fn weaved_iterator() {
        let strct: DenseTensor<u32, OrderedStructure<Euclidean>> = DenseTensor::zero(
            PermutedStructure::from_iter([
                Euclidean {}.new_slot(4, 1),
                Euclidean {}.new_slot(4, 2),
                Euclidean {}.new_slot(4, 3),
                Euclidean {}.new_slot(4, 4),
            ])
            .structure,
        );

        let fiber_spec = [true, false, true, false];
        let self_fiber_class = Fiber::from(fiber_spec.as_slice(), &strct.structure); //We use the partition as a filter here, for indices that belong to self, vs those that belong to other
        let (self_fiber_class_iter, mut _other_fiber_class_iter) =
            CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class); // these are iterators over the open indices of self and other, except expressed in the flat indices of the resulting structure

        for i in self_fiber_class_iter {
            println!("{}-> {:?}", i, strct.expanded_index(i))
        }
    }
}
