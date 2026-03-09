use std::{iter::Sum, ops::AddAssign};

use gat_lending_iterator::LendingIterator;

use crate::{structure::TensorStructure, tensors::data::DataTensor};

pub mod add;
pub mod add_assign;
pub mod algebraic_traits;
pub mod complex;
pub mod fallible_add;
pub mod fallible_sub;
pub mod neg;
pub mod ref_zero;
pub mod scalar_mul;
pub mod sub;
pub mod sub_assign;
pub mod upgrading_arithmetic;

impl<T, S: TensorStructure + Clone> Sum for DataTensor<T, S>
where
    T: for<'a> AddAssign<&'a T> + Clone + Default + PartialEq,
{
    fn sum<I: Iterator<Item = Self>>(mut iter: I) -> Self {
        if let Some(mut i) = iter.next() {
            for j in iter {
                i += j;
            }
            i
        } else {
            panic!("Empty iterator in sum");
        }
    }
}

// pub trait LendingSum: LendingIterator
// where
//     for<'p> Self::Item<'p>: Clone + for<'a> AddAssign<Self::Item>,
// {
//     fn sum_ref<I>(mut iter: I) -> Self {
//         if let Some(mut i) = iter.next() {
//             let sum = i.clone();

//             while let Some(j) = iter.next() {
//                 sum += j;
//             }
//             sum
//         } else {
//             panic!("Empty iterator in sum");
//         }
//     }
// }

impl<T, S: TensorStructure + Clone> DataTensor<T, S>
where
    T: for<'a> AddAssign<&'a T> + Clone + Default + PartialEq,
{
    pub fn sum_ref<I>(mut iter: I) -> Self
    where
        for<'a> I: LendingIterator<Item<'a> = &'a DataTensor<T, S>>,
    {
        if let Some(i) = iter.next() {
            let mut sum = i.clone();

            while let Some(j) = iter.next() {
                sum += j;
            }
            sum
        } else {
            panic!("Empty iterator in sum");
        }
    }

    pub fn sum<I>(mut iter: I) -> Self
    where
        for<'a> I: LendingIterator<Item<'a> = DataTensor<T, S>>,
    {
        if let Some(i) = iter.next() {
            let mut sum = i;

            while let Some(j) = iter.next() {
                sum += j;
            }
            sum
        } else {
            panic!("Empty iterator in sum");
        }
    }
}

pub trait ScalarMul<T> {
    type Output;
    fn scalar_mul(&self, rhs: &T) -> Option<Self::Output>;
}

pub trait ScalarMulMut<T> {
    type Output;
    fn scalar_mul_mut(&self, rhs: &T) -> Option<Self::Output>;
}
