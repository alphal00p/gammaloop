use std::{collections::hash_map::Entry, ops::AddAssign};

use crate::{
    algebra::{algebraic_traits::RefZero, complex::Complex},
    iterators::IteratableTensor,
    structure::{ScalarStructure, TensorStructure},
    tensors::{
        complex::RealOrComplexTensor,
        data::{DataTensor, DenseTensor, SparseOrDense, SparseTensor},
    },
};

impl<T, U, I> AddAssign<DenseTensor<T, I>> for DenseTensor<U, I>
where
    U: for<'a> AddAssign<&'a T>,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: DenseTensor<T, I>) {
        *self += &rhs;
    }
}

impl<T, U, I> AddAssign<&DenseTensor<T, I>> for DenseTensor<U, I>
where
    U: for<'a> AddAssign<&'a T>,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: &DenseTensor<T, I>) {
        for (u, t) in self.data.iter_mut().zip(rhs.data.iter()) {
            *u += t;
        }
    }
}

impl<T, U, I> AddAssign<DenseTensor<T, I>> for SparseTensor<U, I>
where
    U: for<'a> AddAssign<&'a T> + Clone,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: DenseTensor<T, I>) {
        *self += &rhs;
    }
}

impl<T, U, I> AddAssign<&DenseTensor<T, I>> for SparseTensor<U, I>
where
    U: for<'a> AddAssign<&'a T> + Clone,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: &DenseTensor<T, I>) {
        for (i, u) in rhs.iter_flat() {
            *self.elements.entry(i).or_insert(self.zero.clone()) += u;
        }
    }
}

impl<T, U, I> AddAssign<SparseTensor<T, I>> for SparseTensor<U, I>
where
    U: for<'a> AddAssign<&'a T> + Clone,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: SparseTensor<T, I>) {
        *self += &rhs;
    }
}

impl<T, U, I> AddAssign<&SparseTensor<T, I>> for SparseTensor<U, I>
where
    U: for<'a> AddAssign<&'a T> + Clone,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: &SparseTensor<T, I>) {
        for (flat_index, rhs_value) in &rhs.elements {
            match self.elements.entry(*flat_index) {
                Entry::Occupied(mut entry) => {
                    *entry.get_mut() += rhs_value;
                }
                Entry::Vacant(entry) => {
                    let mut new = self.zero.clone();
                    new += rhs_value;
                    entry.insert(new);
                }
            }
        }
    }
}

impl<T, U, I> AddAssign<SparseTensor<T, I>> for DenseTensor<U, I>
where
    U: for<'a> AddAssign<&'a T>,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: SparseTensor<T, I>) {
        *self += &rhs;
    }
}

impl<T, U, I> AddAssign<&SparseTensor<T, I>> for DenseTensor<U, I>
where
    U: for<'a> AddAssign<&'a T>,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: &SparseTensor<T, I>) {
        for (i, u) in rhs.elements.iter() {
            self[*i] += u;
        }
    }
}

impl<T, U, I> AddAssign<DataTensor<T, I>> for DataTensor<U, I>
where
    U: for<'a> AddAssign<&'a T> + Clone + Default + PartialEq,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: DataTensor<T, I>) {
        if self.is_sparse() && rhs.is_dense() {
            self.to_dense_mut();
        }
        match (self, rhs) {
            (DataTensor::Dense(a), DataTensor::Dense(b)) => {
                *a += b;
            }
            (DataTensor::Sparse(a), DataTensor::Sparse(b)) => {
                *a += b;
            }
            (DataTensor::Dense(a), DataTensor::Sparse(b)) => {
                *a += b;
            }

            _ => {
                unreachable!("Sparse dense add assign should be turned into dense dense addassign")
            }
        }
    }
}

impl<T, U, I> AddAssign<&DataTensor<T, I>> for DataTensor<U, I>
where
    U: for<'a> AddAssign<&'a T> + Clone + Default + PartialEq,
    I: TensorStructure + Clone,
{
    fn add_assign(&mut self, rhs: &DataTensor<T, I>) {
        if self.is_sparse() && rhs.is_dense() {
            self.to_dense_mut();
        }

        match (self, rhs) {
            (DataTensor::Dense(a), DataTensor::Dense(b)) => {
                *a += b;
            }
            (DataTensor::Sparse(a), DataTensor::Sparse(b)) => {
                *a += b;
            }
            (DataTensor::Dense(a), DataTensor::Sparse(b)) => {
                *a += b;
            }
            _ => {
                unreachable!("Sparse dense add assign should be turned into dense dense addassign")
            }
        }
    }
}

impl<T, U, I> AddAssign<RealOrComplexTensor<T, I>> for RealOrComplexTensor<U, I>
where
    U: for<'a> AddAssign<&'a T> + RefZero + Clone + Default + PartialEq,
    Complex<U>: for<'a> AddAssign<&'a Complex<T>> + for<'a> AddAssign<&'a T>,
    I: TensorStructure + Clone + ScalarStructure,
{
    fn add_assign(&mut self, rhs: RealOrComplexTensor<T, I>) {
        match (self, rhs) {
            (RealOrComplexTensor::Complex(a), RealOrComplexTensor::Complex(b)) => {
                *a += b;
            }
            (RealOrComplexTensor::Real(a), RealOrComplexTensor::Real(b)) => {
                *a += b;
            }
            (RealOrComplexTensor::Complex(a), RealOrComplexTensor::Real(b)) => {
                *a += b;
            }
            (a, b) => {
                a.to_complex();
                *a += b;
            }
        }
    }
}

impl<T, U, I> AddAssign<&RealOrComplexTensor<T, I>> for RealOrComplexTensor<U, I>
where
    U: for<'a> AddAssign<&'a T> + RefZero + Clone + Default + PartialEq,
    Complex<U>: for<'a> AddAssign<&'a Complex<T>> + for<'a> AddAssign<&'a T>,
    I: TensorStructure + Clone + ScalarStructure,
{
    fn add_assign(&mut self, rhs: &RealOrComplexTensor<T, I>) {
        match (self, rhs) {
            (RealOrComplexTensor::Complex(a), RealOrComplexTensor::Complex(b)) => {
                *a += b;
            }
            (RealOrComplexTensor::Real(a), RealOrComplexTensor::Real(b)) => {
                *a += b;
            }
            (RealOrComplexTensor::Complex(a), RealOrComplexTensor::Real(b)) => {
                *a += b;
            }
            (a, b) => {
                a.to_complex();
                *a += b;
            }
        }
    }
}
