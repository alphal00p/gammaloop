//! Core iterator implementations for tensor iteration
//!
//! This module provides the fundamental iterator types that implement the core logic
//! for iterating through tensor dimensions efficiently.

use std::fmt::Debug;

use super::indices::AbstractFiberIndex;
use super::traits::{
    AbstractFiber, FiberIteratorItem, IteratesAlongFibers, IteratesAlongPermutedFibers,
    ResetableIterator, ShiftableIterator,
};
use crate::structure::{
    concrete_index::ConcreteIndex, concrete_index::FlatIndex, dimension::Dimension,
    representation::RepName, representation::Representation,
};
use bitvec::vec::BitVec;
use linnet::permutation::Permutation;

/// Represents a single stride and shift for fiber iteration
///
/// Used by [CoreFlatFiberIterator] to efficiently iterate through fibers.
#[derive(Debug, Clone, Copy)]
pub struct SingleStrideShift {
    /// The stride to use for iteration
    pub stride: usize,
    /// The shift to apply after each stride
    pub shift: usize,
}

/// Represents multiple strides and shifts for complex fiber iteration
///
/// Used for more complex iteration patterns where a single stride/shift is insufficient.
#[derive(Debug, Clone)]
pub struct MultiStrideShift {
    /// The strides to use for iteration
    pub strides: Vec<usize>,
    /// The shifts to apply after each stride
    pub shifts: Vec<usize>,
}

/// Enumeration of stride/shift patterns for fiber iteration
///
/// Represents either a simple stride/shift pattern or a more complex multi-dimensional one.
#[derive(Debug, Clone)]
pub enum StrideShift {
    /// A single stride/shift pair (or none)
    Single(Option<SingleStrideShift>),
    /// Multiple stride/shift pairs
    Multi(MultiStrideShift),
}

impl StrideShift {
    /// Creates a new single stride/shift pattern
    ///
    /// # Arguments
    ///
    /// * `stride` - Optional stride value
    /// * `shift` - Optional shift value
    pub fn new_single(stride: Option<usize>, shift: Option<usize>) -> Self {
        let stride_shift = stride.zip(shift);
        if let Some(stride_shift) = stride_shift {
            StrideShift::Single(Some(SingleStrideShift {
                stride: stride_shift.0,
                shift: stride_shift.1,
            }))
        } else {
            StrideShift::Single(None)
        }
    }

    /// Creates a new empty single stride/shift pattern
    pub fn new_single_none() -> Self {
        StrideShift::Single(None)
    }

    /// Creates a new multi stride/shift pattern
    ///
    /// # Arguments
    ///
    /// * `strides` - The strides to use
    /// * `shifts` - The shifts to apply
    pub fn new_multi(strides: Vec<usize>, shifts: Vec<usize>) -> Self {
        StrideShift::Multi(MultiStrideShift { strides, shifts })
    }
}

/// Wrapper for items that can be skipped during iteration
///
/// Used by sparse tensor iterators to keep track of skipped positions.
pub struct SkippingItem<I: FiberIteratorItem> {
    /// Number of indices skipped
    pub skips: usize,
    /// The underlying item
    pub item: I,
}

impl<I: FiberIteratorItem> FiberIteratorItem for SkippingItem<I> {
    type OtherData = (usize, I::OtherData);
    fn flat_idx(&self) -> FlatIndex {
        self.item.flat_idx()
    }
    fn other_data(self) -> Self::OtherData {
        (self.skips, self.item.other_data())
    }
}

/// Wrapper for items that include sign/metric information
///
/// Used by metric iterators to track sign information.
pub struct MetricItem<I: FiberIteratorItem> {
    /// Whether the metric is negative
    pub neg: bool,
    /// The underlying item
    pub item: I,
}

impl<I: FiberIteratorItem> FiberIteratorItem for MetricItem<I> {
    type OtherData = (bool, I::OtherData);
    fn flat_idx(&self) -> FlatIndex {
        self.item.flat_idx()
    }
    fn other_data(self) -> Self::OtherData {
        (self.neg, self.item.other_data())
    }
}

/// Core flat fiber iterator for efficient traversal of tensor fibers
///
/// Optimized for iterating along a single dimension using flat indices.
#[derive(Debug, Clone)]
pub struct CoreFlatFiberIterator {
    /// The current fiber index
    pub varying_fiber_index: FlatIndex,
    /// The increment to apply at each step
    pub increment: FlatIndex,
    /// The stride and shift pattern
    pub stride_shift: StrideShift,
    /// The maximum index value
    pub max: FlatIndex,
    /// The starting index offset
    pub zero_index: FlatIndex,
}

impl CoreFlatFiberIterator {
    /// Initializes a multi-fiber iterator
    ///
    /// Used for iterating over multiple fibers simultaneously.
    fn init_multi_fiber_iter<I, J>(
        strides: Vec<usize>,
        dims: Vec<Dimension>,
        order: usize,
        fiber: &I,
        conj: bool,
    ) -> (FlatIndex, Vec<usize>, Vec<usize>, FlatIndex)
    where
        I: std::ops::Index<usize, Output = J>,
        J: AbstractFiberIndex,
    {
        let mut max = 0;

        let mut increment = 1;

        let mut fixed_strides = vec![];
        let mut shifts = vec![];

        let mut before = true;
        let mut has_seen_stride = false;
        let mut first = true;

        for pos in (0..order).rev() {
            let fi = &fiber[pos];

            if fi.is_fixed() ^ conj && !before && !first {
                has_seen_stride = true;
                fixed_strides.push(strides[pos]);
            }
            if fi.is_free() ^ conj && before && has_seen_stride {
                shifts.push(strides[pos]);
            }

            if fi.is_free() ^ conj {
                let dimminus1: usize = match usize::try_from(dims[pos]).unwrap() {
                    0 => 0,
                    _ => usize::try_from(dims[pos]).unwrap() - 1,
                };
                max += dimminus1 * strides[pos];
                if first {
                    increment = strides[pos];
                    first = false;
                }
            }

            before = fi.is_fixed() ^ conj;
        }

        if fixed_strides.len() > shifts.len() {
            fixed_strides.pop();
        }
        (increment.into(), fixed_strides, shifts, max.into())
    }

    /// Initializes a single fiber iterator
    ///
    /// Optimized for iterating along a single dimension.
    fn init_single_fiber_iter(
        strides: Vec<usize>,
        fiber_position: usize,
        dims: Vec<Dimension>,
        conj: bool,
    ) -> (FlatIndex, Option<usize>, Option<usize>, FlatIndex) {
        let fiber_stride = strides[fiber_position];
        let dim: usize = dims[fiber_position].try_into().unwrap();
        let size = dims
            .iter()
            .map(|x| usize::try_from(*x).unwrap())
            .product::<usize>();
        let mut stride = None;
        let mut shift = None;

        if conj {
            let max = size - fiber_stride * (dim - 1) - 1;

            let mut increment = 1;

            if fiber_position == dims.len() - 1 {
                increment = *strides.get(dims.len().wrapping_sub(2)).unwrap_or(&1);
            } else if fiber_position != 0 {
                shift = Some(strides[fiber_position - 1]);
                stride = Some(strides[fiber_position]);
            }

            (increment.into(), stride, shift, max.into())
        } else {
            let increment = fiber_stride;
            let max = fiber_stride * (dim - 1);

            (increment.into(), stride, shift, max.into())
        }
    }
}

impl Iterator for CoreFlatFiberIterator {
    type Item = FlatIndex;
    fn next(&mut self) -> Option<Self::Item> {
        if self.varying_fiber_index > self.max {
            return None;
        }
        let index = self.varying_fiber_index + self.zero_index;

        self.varying_fiber_index += self.increment;

        match self.stride_shift {
            StrideShift::Multi(ref ss) => {
                for (i, s) in ss.strides.iter().enumerate() {
                    if self.varying_fiber_index % s == 0.into() {
                        self.varying_fiber_index += (ss.shifts[i] - s).into();
                    } else {
                        break;
                    }
                }
            }
            StrideShift::Single(Some(ss)) => {
                if self.varying_fiber_index % ss.stride == 0.into() {
                    self.varying_fiber_index += (ss.shift - ss.stride).into();
                }
            }
            _ => {}
        }
        Some(index)
    }
}

impl ResetableIterator for CoreFlatFiberIterator {
    fn reset(&mut self) {
        self.varying_fiber_index = 0.into();
    }
}

impl ShiftableIterator for CoreFlatFiberIterator {
    fn shift(&mut self, shift: usize) {
        self.zero_index = shift.into();
    }
}

impl<R: RepName> IteratesAlongFibers<R> for CoreFlatFiberIterator {
    fn new<I, J>(fiber: &I, conj: bool) -> Self
    where
        I: AbstractFiber<J, Repr = R>,
        J: AbstractFiberIndex,
        Self: Sized,
    {
        if let Some(single) = fiber.single() {
            let (increment, fixed_strides, shifts, max) =
                Self::init_single_fiber_iter(fiber.strides(), single, fiber.shape(), conj);

            CoreFlatFiberIterator {
                increment,
                stride_shift: StrideShift::new_single(fixed_strides, shifts),
                max,
                zero_index: 0.into(),
                varying_fiber_index: 0.into(),
            }
        } else {
            let (increment, fixed_strides, shifts, max) = Self::init_multi_fiber_iter(
                fiber.strides(),
                fiber.shape(),
                fiber.order(),
                fiber,
                conj,
            );

            CoreFlatFiberIterator {
                increment,
                stride_shift: StrideShift::new_multi(fixed_strides, shifts),
                max,
                zero_index: 0.into(),
                varying_fiber_index: 0.into(),
            }
        }
    }

    fn new_paired_conjugates<I, J>(fiber: &I) -> (Self, Self)
    where
        I: AbstractFiber<J, Repr = R>,
        J: AbstractFiberIndex,
    {
        let strides = fiber.strides();
        let dims = fiber.shape();
        let order = fiber.order();
        let mut max = 0;

        let mut increment = 1;

        let mut fixed_strides = vec![];
        let mut fixed_strides_conj = vec![];
        let mut shifts = vec![];
        let mut shifts_conj = vec![];

        let mut before = true;
        let mut has_seen_stride = false;
        let mut has_seen_stride_conj = false;
        let mut first = true;
        let mut first_conj = true;
        let mut increment_conj = 1;

        let mut max_conj = 0;

        for pos in (0..order).rev() {
            let fi = &fiber[pos];

            if fi.is_fixed() && !before {
                if !first {
                    has_seen_stride = true;
                    fixed_strides.push(strides[pos]);
                }

                if has_seen_stride_conj {
                    shifts_conj.push(strides[pos]);
                }
            }
            if fi.is_free() && before {
                if has_seen_stride {
                    shifts.push(strides[pos]);
                }

                if !first_conj {
                    fixed_strides_conj.push(strides[pos]);
                    has_seen_stride_conj = true;
                }
            }

            if fi.is_fixed() {
                max_conj += (usize::try_from(dims[pos]).unwrap() - 1) * strides[pos];
                if first_conj {
                    increment_conj = strides[pos];
                    first_conj = false;
                }
            } else {
                max += (usize::try_from(dims[pos]).unwrap() - 1) * strides[pos];
                if first {
                    increment = strides[pos];
                    first = false;
                }
            }

            before = fi.is_fixed();
        }

        if fixed_strides.len() > shifts.len() {
            fixed_strides.pop();
        }

        if fixed_strides_conj.len() > shifts_conj.len() {
            fixed_strides_conj.pop();
        }

        (
            CoreFlatFiberIterator {
                varying_fiber_index: 0.into(),
                stride_shift: StrideShift::new_multi(fixed_strides_conj, shifts_conj),
                increment: increment_conj.into(),
                max: max_conj.into(),
                zero_index: 0.into(),
            },
            CoreFlatFiberIterator {
                varying_fiber_index: 0.into(),
                increment: increment.into(),
                stride_shift: StrideShift::new_multi(fixed_strides, shifts),
                max: max.into(),
                zero_index: 0.into(),
            },
        )
    }
}

/// Core expanded fiber iterator for traversing tensor fibers using expanded indices
///
/// More flexible than `CoreFlatFiberIterator` as it can handle permutations and
/// dimension-specific logic.
#[derive(Debug, Clone)]
pub struct CoreExpandedFiberIterator<R: RepName> {
    /// The current expanded indices
    pub varying_fiber_index: Vec<ConcreteIndex>,
    /// The dimensions of the fiber
    pub dims: Vec<Representation<R>>,
    /// The strides for each dimension
    pub strides: Vec<usize>,
    /// The starting index offset
    pub zero_index: FlatIndex,
    /// The current flat index
    pub flat: FlatIndex,
    /// Whether iteration is complete
    exhausted: bool,
}

impl<R: RepName> CoreExpandedFiberIterator<R> {
    /// Initializes a new expanded fiber iterator
    /// The fixed indices of the fiber are not taken into account
    fn init_iter<I, J>(fiber: &I, conj: bool, permutation: Option<Permutation>) -> Self
    where
        I: AbstractFiber<J, Repr = R>,
        J: AbstractFiberIndex,
    {
        let varying_indices = fiber.bitvec();
        let mut dims = Self::filter(&varying_indices, &fiber.reps(), conj);

        let mut strides = Self::filter(&varying_indices, &fiber.strides(), conj);
        let varying_fiber_index = vec![0; dims.len()];

        if let Some(perm) = permutation {
            perm.apply_slice_in_place(&mut dims);
            perm.apply_slice_in_place(&mut strides);
        }

        CoreExpandedFiberIterator {
            varying_fiber_index,
            zero_index: 0.into(),
            dims,
            strides,
            flat: 0.into(),
            exhausted: false,
        }
    }

    /// Filters elements from a vector based on a bitvector filter
    fn filter<T: Clone>(filter: &BitVec, vec: &[T], conj: bool) -> Vec<T> {
        let mut res = vec![];
        for (i, x) in filter.iter().enumerate() {
            if conj ^ *x {
                res.push(vec[i].clone());
            }
        }
        res
    }
}

impl<R: RepName> Iterator for CoreExpandedFiberIterator<R> {
    type Item = FlatIndex;
    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            return None;
        }

        let current_flat = self.flat + self.zero_index; // Store the current flat value before modifications

        let mut carry = true;
        for ((pos, dim), stride) in self
            .varying_fiber_index
            .iter_mut()
            .zip(self.dims.iter())
            .zip(self.strides.iter())
            .rev()
        {
            if carry {
                *pos += 1;
                if *pos >= usize::try_from(*dim).unwrap() {
                    *pos = 0;
                    self.flat -= (stride * (usize::try_from(*dim).unwrap() - 1)).into();
                } else {
                    self.flat += (*stride).into();
                    carry = false;
                }
            }
        }

        if carry {
            self.exhausted = true; // Set the flag to prevent further iterations after this one
        }

        Some(current_flat)
    }
}

impl<R: RepName> ShiftableIterator for CoreExpandedFiberIterator<R> {
    fn shift(&mut self, shift: usize) {
        self.zero_index = shift.into();
    }
}

impl<R: RepName> ResetableIterator for CoreExpandedFiberIterator<R> {
    fn reset(&mut self) {
        self.flat = 0.into();
        self.exhausted = false;
        self.varying_fiber_index = vec![0; self.dims.len()];
    }
}

impl<R: RepName> IteratesAlongFibers<R> for CoreExpandedFiberIterator<R> {
    fn new<I, J>(fiber: &I, conj: bool) -> Self
    where
        I: AbstractFiber<J, Repr = R>,
        J: AbstractFiberIndex,
    {
        Self::init_iter(fiber, conj, None)
    }

    fn new_paired_conjugates<I, J>(fiber: &I) -> (Self, Self)
    where
        I: AbstractFiber<J, Repr = R>,
        J: AbstractFiberIndex,
        Self: Sized,
    {
        (
            Self::init_iter(fiber, true, None),
            Self::init_iter(fiber, false, None),
        )
    }
}

impl<R: RepName> IteratesAlongPermutedFibers<R> for CoreExpandedFiberIterator<R> {
    fn new_permuted<I, J>(fiber: &I, conj: bool, permutation: Permutation) -> Self
    where
        I: AbstractFiber<J, Repr = R>,
        J: AbstractFiberIndex,
    {
        Self::init_iter(fiber, conj, Some(permutation))
    }
}

/// Metric fiber iterator for traversing tensor fibers with metric information
///
/// Extends `CoreExpandedFiberIterator` to track sign/metric information during iteration.
#[derive(Debug, Clone)]
pub struct MetricFiberIterator<R: RepName> {
    /// The underlying expanded iterator
    pub iter: CoreExpandedFiberIterator<R>,
    /// Whether the current item has negative sign
    neg: bool,
}

impl<R: RepName> ResetableIterator for MetricFiberIterator<R> {
    fn reset(&mut self) {
        self.iter.reset();
    }
}

impl<R: RepName> ShiftableIterator for MetricFiberIterator<R> {
    fn shift(&mut self, shift: usize) {
        self.iter.shift(shift);
    }
}

impl<R: RepName> IteratesAlongFibers<R> for MetricFiberIterator<R> {
    fn new<I, J>(fiber: &I, conj: bool) -> Self
    where
        I: AbstractFiber<J, Repr = R>,
        J: AbstractFiberIndex,
    {
        MetricFiberIterator {
            iter: CoreExpandedFiberIterator::new(fiber, conj),
            neg: false,
        }
    }

    fn new_paired_conjugates<I, J>(fiber: &I) -> (Self, Self)
    where
        I: AbstractFiber<J, Repr = R>,
        J: AbstractFiberIndex,
        Self: Sized,
    {
        (
            MetricFiberIterator {
                iter: CoreExpandedFiberIterator::new(fiber, true),
                neg: false,
            },
            MetricFiberIterator {
                iter: CoreExpandedFiberIterator::new(fiber, false),
                neg: false,
            },
        )
    }
}

impl<R: RepName> IteratesAlongPermutedFibers<R> for MetricFiberIterator<R> {
    fn new_permuted<I, J>(fiber: &I, conj: bool, permutation: Permutation) -> Self
    where
        I: AbstractFiber<J, Repr = R>,
        J: AbstractFiberIndex,
    {
        MetricFiberIterator {
            iter: CoreExpandedFiberIterator::new_permuted(fiber, conj, permutation),
            neg: false,
        }
    }
}

impl<R: RepName> Iterator for MetricFiberIterator<R> {
    type Item = MetricItem<FlatIndex>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.iter.exhausted {
            return None;
        }

        let current_flat = self.iter.flat + self.iter.zero_index; // Store the current flat value before modifications

        let mut carry = true;
        self.neg = false;
        for ((pos, dim), stride) in self
            .iter
            .varying_fiber_index
            .iter_mut()
            .zip(self.iter.dims.iter())
            .zip(self.iter.strides.iter())
            .rev()
        {
            self.neg ^= dim.is_neg(*pos);
            if carry {
                *pos += 1;
                if *pos >= usize::try_from(*dim).unwrap() {
                    *pos = 0;
                    self.iter.flat -= (stride * (usize::try_from(*dim).unwrap() - 1)).into();
                } else {
                    self.iter.flat += (*stride).into();
                    carry = false;
                }
            }
        }

        if carry {
            self.iter.exhausted = true; // Set the flag to prevent further iterations after this one
        }

        Some(MetricItem {
            neg: self.neg,
            item: current_flat,
        })
    }
}
