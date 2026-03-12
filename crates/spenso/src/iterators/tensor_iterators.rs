//! Iterators for specific tensor types
//!
//! This module provides iterators specific to different tensor implementations,
//! including dense, sparse, and data tensors.

use std::collections::hash_map;

use crate::{
    algebra::algebraic_traits::RefZero,
    algebra::upgrading_arithmetic::{FallibleAddAssign, FallibleSubAssign},
    contraction::ContractableWith,
    structure::{
        concrete_index::{ConcreteIndex, ExpandedIndex, FlatIndex},
        TensorStructure,
    },
    tensors::data::{DataTensor, DenseTensor, GetTensorData, SparseTensor},
};

use super::traits::IteratableTensor;

/// Iterator over all the elements of a sparse tensor
///
/// Returns the expanded indices and the element at that index
pub struct SparseTensorIterator<'a, T, N> {
    iter: hash_map::Iter<'a, FlatIndex, T>,
    structure: &'a N,
}

impl<'a, T, N> SparseTensorIterator<'a, T, N> {
    /// Creates a new sparse tensor iterator
    ///
    /// # Arguments
    ///
    /// * `tensor` - The tensor to iterate over
    pub fn new(tensor: &'a SparseTensor<T, N>) -> Self {
        SparseTensorIterator {
            iter: tensor.elements.iter(),
            structure: &tensor.structure,
        }
    }
}

impl<'a, T, N> Iterator for SparseTensorIterator<'a, T, N>
where
    N: TensorStructure,
{
    type Item = (ExpandedIndex, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some((k, v)) = self.iter.next() {
            let indices = self.structure.expanded_index(*k).unwrap();
            Some((indices, v))
        } else {
            None
        }
    }
}

/// Iterator over all the elements of a sparse tensor
///
/// Returns the flat index and the element at that index
pub struct SparseTensorLinearIterator<'a, T> {
    iter: hash_map::Iter<'a, FlatIndex, T>,
}

impl<'a, T> SparseTensorLinearIterator<'a, T> {
    /// Creates a new sparse tensor linear iterator
    ///
    /// # Arguments
    ///
    /// * `tensor` - The tensor to iterate over
    pub fn new<N>(tensor: &'a SparseTensor<T, N>) -> Self {
        SparseTensorLinearIterator {
            iter: tensor.elements.iter(),
        }
    }
}

impl<'a, T> Iterator for SparseTensorLinearIterator<'a, T> {
    type Item = (FlatIndex, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|(k, v)| (*k, v))
    }
}

/// Iterator over all but two indices of a sparse tensor, where the two indices are traced
///
/// The iterator next returns the value of the trace at the current indices, and the current indices
pub struct SparseTensorTraceIterator<'a, T, I> {
    tensor: &'a SparseTensor<T, I>,
    trace_indices: [usize; 2],
    current_indices: Vec<ConcreteIndex>,
    done: bool,
}

impl<'a, T, I> SparseTensorTraceIterator<'a, T, I>
where
    I: TensorStructure,
{
    /// Create a new trace iterator
    ///
    /// # Arguments
    ///
    /// * `tensor` - A reference to the tensor
    /// * `trace_indices` - The indices to be traced
    ///
    /// # Panics
    ///
    /// Panics if the trace indices do not point to the same dimension
    pub fn new(tensor: &'a SparseTensor<T, I>, trace_indices: [usize; 2]) -> Self {
        //trace positions must point to the same dimension
        assert!(
            trace_indices
                .iter()
                .map(|&pos| tensor.get_rep(pos).unwrap())
                .collect::<Vec<_>>()
                .iter()
                .all(|&sig| sig == tensor.get_rep(trace_indices[0]).unwrap()),
            "Trace indices must point to the same dimension"
        );
        SparseTensorTraceIterator {
            tensor,
            trace_indices,
            current_indices: vec![0; tensor.order()],
            done: false,
        }
    }

    /// Increments the current indices, skipping trace indices
    ///
    /// Returns true if incrementation succeeded, false if we've reached the end
    fn increment_indices(&mut self) -> bool {
        for (i, index) in self
            .current_indices
            .iter_mut()
            .enumerate()
            .rev()
            .filter(|(pos, _)| !self.trace_indices.contains(pos))
        // Filter out the trace indices
        {
            *index += 1;
            // If the index goes beyond the shape boundary, wrap around to 0
            if index >= &mut usize::try_from(self.tensor.shape()[i]).unwrap() {
                *index = 0;
                continue; // carry over to the next dimension
            }
            return true; // We've successfully found the next combination
        }
        false // No more combinations left
    }
}

impl<T, I> Iterator for SparseTensorTraceIterator<'_, T, I>
where
    T: ContractableWith<T> + FallibleAddAssign<T> + FallibleSubAssign<T> + Clone + RefZero,
    I: TensorStructure + Clone,
{
    type Item = (Vec<ConcreteIndex>, T);
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let trace_dimension = self.tensor.get_rep(self.trace_indices[0]).unwrap();
        let trace_sign = trace_dimension.negative().unwrap();
        let mut iter = trace_sign.iter().enumerate();
        let mut indices = self.current_indices.clone();
        let (i, mut sign) = iter.next().unwrap(); //First element (to eliminate the need for default)

        indices[self.trace_indices[0]] = i;
        indices[self.trace_indices[1]] = i;

        // Data might not exist at that concrete index usize, we advance it till it does, and if not we skip
        while self.tensor.is_empty_at_expanded(&indices) {
            let Some((i, signint)) = iter.next() else {
                self.done = !self.increment_indices();
                return self.next(); // skip
            };
            indices[self.trace_indices[0]] = i;
            indices[self.trace_indices[1]] = i;
            sign = signint;
        }

        let value = self.tensor.get_ref(&indices).unwrap(); //Should now be safe to unwrap
        let zero = value.ref_zero();

        let mut trace = if *sign {
            let mut zero = zero.clone();
            zero.sub_assign_fallible(value);
            zero
        } else {
            value.clone()
        };

        for (i, sign) in iter {
            indices[self.trace_indices[0]] = i;
            indices[self.trace_indices[1]] = i;
            if let Ok(value) = self.tensor.get_ref(&indices) {
                if *sign {
                    trace.sub_assign_fallible(value);
                } else {
                    trace.add_assign_fallible(value);
                }
            }
        }

        //make a vector without the trace indices
        let trace_indices: Vec<ConcreteIndex> = self
            .current_indices
            .clone()
            .into_iter()
            .enumerate()
            .filter(|&(i, _)| !self.trace_indices.contains(&i))
            .map(|(_, x)| x)
            .collect();

        self.done = !self.increment_indices();

        Some((trace_indices, trace))
    }
}

impl<T, S: TensorStructure> SparseTensor<T, S> {
    /// Creates an iterator that traces over two indices
    ///
    /// # Arguments
    ///
    /// * `trace_indices` - The indices to trace over
    pub fn iter_trace(&self, trace_indices: [usize; 2]) -> SparseTensorTraceIterator<'_, T, S> {
        SparseTensorTraceIterator::new(self, trace_indices)
    }
}

/// Iterator over all the elements of a dense tensor
///
/// Returns the expanded indices and the element at that index
pub struct DenseTensorIterator<'a, T, I> {
    tensor: &'a DenseTensor<T, I>,
    current_flat_index: FlatIndex,
}

impl<'a, T, I> DenseTensorIterator<'a, T, I> {
    /// Create a new dense tensor iterator
    ///
    /// # Arguments
    ///
    /// * `tensor` - A reference to the tensor
    pub fn new(tensor: &'a DenseTensor<T, I>) -> Self {
        DenseTensorIterator {
            tensor,
            current_flat_index: 0.into(),
        }
    }
}

impl<'a, T, I> Iterator for DenseTensorIterator<'a, T, I>
where
    I: TensorStructure,
{
    type Item = (ExpandedIndex, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        if let Ok(indices) = self.tensor.expanded_index(self.current_flat_index) {
            let value = self.tensor.get_ref_linear(self.current_flat_index).unwrap();

            self.current_flat_index += 1.into();

            Some((indices, value))
        } else {
            None
        }
    }
}

/// Iterator over all the elements of a dense tensor
///
/// Returns the flat index and the element at that index
pub struct DenseTensorLinearIterator<'a, T, I> {
    tensor: &'a DenseTensor<T, I>,
    current_flat_index: FlatIndex,
}

impl<'a, T, I> DenseTensorLinearIterator<'a, T, I> {
    /// Creates a new dense tensor linear iterator
    ///
    /// # Arguments
    ///
    /// * `tensor` - The tensor to iterate over
    pub fn new(tensor: &'a DenseTensor<T, I>) -> Self {
        DenseTensorLinearIterator {
            tensor,
            current_flat_index: 0.into(),
        }
    }
}

impl<'a, T, I> Iterator for DenseTensorLinearIterator<'a, T, I>
where
    I: TensorStructure,
{
    type Item = (FlatIndex, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        let value = self.tensor.get_ref_linear(self.current_flat_index)?;
        let index = self.current_flat_index;
        self.current_flat_index += 1.into();
        Some((index, value))
    }
}

impl<'a, T, I> IntoIterator for &'a DenseTensor<T, I>
where
    I: TensorStructure,
{
    type Item = (ExpandedIndex, &'a T);
    type IntoIter = DenseTensorIterator<'a, T, I>;

    fn into_iter(self) -> Self::IntoIter {
        DenseTensorIterator::new(self)
    }
}

/// An consuming iterator over the elements of a dense tensor
///
/// Returns the expanded indices and the element at that index
pub struct DenseTensorIntoIterator<T, I> {
    tensor: DenseTensor<T, I>,
    current_flat_index: FlatIndex,
}

impl<T, I> DenseTensorIntoIterator<T, I> {
    /// Creates a new consuming iterator over a dense tensor
    ///
    /// # Arguments
    ///
    /// * `tensor` - The tensor to consume
    pub fn new(tensor: DenseTensor<T, I>) -> Self {
        DenseTensorIntoIterator {
            tensor,
            current_flat_index: 0.into(),
        }
    }
}

impl<T, I> Iterator for DenseTensorIntoIterator<T, I>
where
    I: TensorStructure,
{
    type Item = (ExpandedIndex, T);

    fn next(&mut self) -> Option<Self::Item> {
        if let Ok(indices) = self.tensor.expanded_index(self.current_flat_index) {
            let indices = indices.clone();
            let value = self.tensor.data.remove(self.current_flat_index.into());

            self.current_flat_index += 1.into();

            Some((indices, value))
        } else {
            None
        }
    }
}

impl<T, I> IntoIterator for DenseTensor<T, I>
where
    I: TensorStructure,
{
    type Item = (ExpandedIndex, T);
    type IntoIter = DenseTensorIntoIterator<T, I>;

    fn into_iter(self) -> Self::IntoIter {
        DenseTensorIntoIterator::new(self)
    }
}

/// Iterator over all indices of a dense tensor, keeping two indices fixed and tracing over them
///
/// The next method returns the value of the trace at the current indices, and the current indices
pub struct DenseTensorTraceIterator<'a, T, I> {
    tensor: &'a DenseTensor<T, I>,
    trace_indices: [usize; 2],
    current_indices: Vec<ConcreteIndex>,
    done: bool,
}

impl<'a, T, I> DenseTensorTraceIterator<'a, T, I>
where
    I: TensorStructure,
{
    /// Create a new trace iterator
    ///
    /// # Arguments
    ///
    /// * `tensor` - A reference to the tensor
    /// * `trace_indices` - The indices to be traced
    ///
    pub fn new(tensor: &'a DenseTensor<T, I>, trace_indices: [usize; 2]) -> Self {
        assert!(trace_indices.len() >= 2, "Invalid trace indices");
        //trace positions must point to the same dimension
        assert!(
            trace_indices
                .iter()
                .map(|&pos| tensor.get_rep(pos))
                .collect::<Vec<_>>()
                .iter()
                .all(|&sig| sig == tensor.get_rep(trace_indices[0])),
            "Trace indices must point to the same dimension"
        );
        DenseTensorTraceIterator {
            tensor,
            trace_indices,
            current_indices: vec![0; tensor.order()],
            done: false,
        }
    }

    /// Increments the current indices, skipping trace indices
    ///
    /// Returns true if incrementation succeeded, false if we've reached the end
    fn increment_indices(&mut self) -> bool {
        for (i, index) in self
            .current_indices
            .iter_mut()
            .enumerate()
            .rev()
            .filter(|(pos, _)| !self.trace_indices.contains(pos))
        {
            *index += 1;
            // If the index goes beyond the shape boundary, wrap around to 0
            if index >= &mut usize::try_from(self.tensor.shape()[i]).unwrap() {
                *index = 0;
                continue; // carry over to the next dimension
            }
            return true; // We've successfully found the next combination
        }
        false // No more combinations left
    }
}

impl<T, I> Iterator for DenseTensorTraceIterator<'_, T, I>
where
    T: ContractableWith<T, Out = T> + FallibleAddAssign<T> + FallibleSubAssign<T> + Clone + RefZero,
    I: TensorStructure,
{
    type Item = (Vec<ConcreteIndex>, T);
    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }

        let trace_dimension = self.tensor.get_rep(self.trace_indices[0]).unwrap();
        let trace_sign = trace_dimension.negative().unwrap();

        let mut iter = trace_sign.iter().enumerate();
        let mut indices = self.current_indices.clone();
        let (_, sign) = iter.next().unwrap(); //First sign

        for pos in self.trace_indices {
            indices[pos] = 0;
        }

        let value = self.tensor.get_ref(&indices).unwrap();
        let zero = value.ref_zero();

        let mut trace = if *sign {
            let mut zero = zero.clone();
            zero.sub_assign_fallible(value);
            zero
        } else {
            value.clone()
        };

        for (i, sign) in iter {
            for pos in self.trace_indices {
                indices[pos] = i;
            }

            if let Ok(value) = self.tensor.get_ref(&indices) {
                if *sign {
                    trace.sub_assign_fallible(value);
                } else {
                    trace.add_assign_fallible(value);
                }
            }
        }

        //make a vector without the trace indices
        let trace_indices: Vec<ConcreteIndex> = self
            .current_indices
            .clone()
            .into_iter()
            .enumerate()
            .filter(|&(i, _)| !self.trace_indices.contains(&i))
            .map(|(_, x)| x)
            .collect();

        self.done = !self.increment_indices();

        Some((trace_indices, trace))
    }
}

impl<T, I> DenseTensor<T, I>
where
    I: TensorStructure,
{
    /// Creates an iterator that traces over two indices
    ///
    /// # Arguments
    ///
    /// * `trace_indices` - The indices to trace over
    pub fn iter_trace(&self, trace_indices: [usize; 2]) -> DenseTensorTraceIterator<'_, T, I> {
        DenseTensorTraceIterator::new(self, trace_indices)
    }
}

/// Linear iterator that works with both dense and sparse tensors
///
/// Returns flat indices and elements.
pub enum DataTensorLinearIterator<'a, T, S> {
    /// Iterator for dense tensors
    Dense(DenseTensorLinearIterator<'a, T, S>),
    /// Iterator for sparse tensors
    Sparse(SparseTensorLinearIterator<'a, T>),
}

impl<'a, T, S> From<DenseTensorLinearIterator<'a, T, S>> for DataTensorLinearIterator<'a, T, S> {
    fn from(value: DenseTensorLinearIterator<'a, T, S>) -> Self {
        DataTensorLinearIterator::Dense(value)
    }
}

impl<'a, T, S> From<SparseTensorLinearIterator<'a, T>> for DataTensorLinearIterator<'a, T, S> {
    fn from(value: SparseTensorLinearIterator<'a, T>) -> Self {
        DataTensorLinearIterator::Sparse(value)
    }
}

impl<'a, T, S: TensorStructure> Iterator for DataTensorLinearIterator<'a, T, S> {
    type Item = (FlatIndex, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            DataTensorLinearIterator::Dense(iter) => iter.next(),
            DataTensorLinearIterator::Sparse(iter) => iter.next(),
        }
    }
}

/// Expanded iterator that works with both dense and sparse tensors
///
/// Returns expanded indices and elements.
pub enum DataTensorExpandedIterator<'a, T, S> {
    /// Iterator for dense tensors
    Dense(DenseTensorIterator<'a, T, S>),
    /// Iterator for sparse tensors
    Sparse(SparseTensorIterator<'a, T, S>),
}

impl<'a, T, S> From<DenseTensorIterator<'a, T, S>> for DataTensorExpandedIterator<'a, T, S> {
    fn from(value: DenseTensorIterator<'a, T, S>) -> Self {
        DataTensorExpandedIterator::Dense(value)
    }
}

impl<'a, T, S> From<SparseTensorIterator<'a, T, S>> for DataTensorExpandedIterator<'a, T, S> {
    fn from(value: SparseTensorIterator<'a, T, S>) -> Self {
        DataTensorExpandedIterator::Sparse(value)
    }
}

impl<'a, T, S: TensorStructure> Iterator for DataTensorExpandedIterator<'a, T, S> {
    type Item = (ExpandedIndex, &'a T);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            DataTensorExpandedIterator::Dense(iter) => iter.next(),
            DataTensorExpandedIterator::Sparse(iter) => iter.next(),
        }
    }
}

impl<T: Clone, S: TensorStructure> IteratableTensor for DataTensor<T, S> {
    type Data<'a>
        = &'a T
    where
        Self: 'a;

    fn iter_expanded(&self) -> impl Iterator<Item = (ExpandedIndex, Self::Data<'_>)> {
        match self {
            DataTensor::Dense(tensor) => {
                DataTensorExpandedIterator::Dense(DenseTensorIterator::new(tensor))
            }
            DataTensor::Sparse(tensor) => {
                DataTensorExpandedIterator::Sparse(SparseTensorIterator::new(tensor))
            }
        }
    }

    fn iter_flat(&self) -> impl Iterator<Item = (FlatIndex, Self::Data<'_>)> {
        match self {
            DataTensor::Dense(tensor) => {
                DataTensorLinearIterator::Dense(DenseTensorLinearIterator::new(tensor))
            }
            DataTensor::Sparse(tensor) => {
                DataTensorLinearIterator::Sparse(SparseTensorLinearIterator::new(tensor))
            }
        }
    }
}

impl<T, S: TensorStructure> IteratableTensor for SparseTensor<T, S> {
    type Data<'a>
        = &'a T
    where
        Self: 'a;

    fn iter_expanded(&self) -> impl Iterator<Item = (ExpandedIndex, Self::Data<'_>)> {
        SparseTensorIterator::new(self)
    }

    fn iter_flat(&self) -> impl Iterator<Item = (FlatIndex, Self::Data<'_>)> {
        SparseTensorLinearIterator::new(self)
    }
}

impl<T, S> IteratableTensor for DenseTensor<T, S>
where
    S: TensorStructure,
{
    type Data<'a>
        = &'a T
    where
        Self: 'a;

    fn iter_expanded(&self) -> impl Iterator<Item = (ExpandedIndex, &T)> {
        DenseTensorIterator::new(self)
    }

    fn iter_flat(&self) -> impl Iterator<Item = (FlatIndex, &T)> {
        DenseTensorLinearIterator::new(self)
    }
}

/// An iterator over all indices of a tensor structure
///
/// Returns expanded indices
pub struct TensorStructureIndexIterator<'a, I: TensorStructure> {
    structure: &'a I,
    current_flat_index: FlatIndex,
}

impl<I: TensorStructure> Iterator for TensorStructureIndexIterator<'_, I> {
    type Item = ExpandedIndex;
    fn next(&mut self) -> Option<Self::Item> {
        if let Ok(indices) = self.structure.expanded_index(self.current_flat_index) {
            self.current_flat_index += 1.into();

            Some(indices)
        } else {
            None
        }
    }
}

impl<'a, I: TensorStructure> TensorStructureIndexIterator<'a, I> {
    /// Creates a new tensor structure index iterator
    ///
    /// # Arguments
    ///
    /// * `structure` - The tensor structure to iterate over
    #[must_use]
    pub fn new(structure: &'a I) -> Self {
        TensorStructureIndexIterator {
            structure,
            current_flat_index: 0.into(),
        }
    }
}
