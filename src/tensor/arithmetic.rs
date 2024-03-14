use crate::tensor::{structure, ConcreteIndex, GetTensorData, SetTensorData};

use super::{DataTensor, DenseTensor, FallibleAdd, FallibleSub, SparseTensor, TensorStructure};
use ahash::AHashMap;
use num::traits::Num;
use std::ops::{Add, Mul, Sub};

impl<'a, T, U, I, Out> FallibleAdd<&DenseTensor<T, I>> for &'a DenseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleAdd<&'c T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn add_fallible(self, rhs: &DenseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let permutation = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = rhs.structure().clone();

        let data: Option<Vec<Out>> = rhs
            .data
            .iter()
            .enumerate()
            .map(|(i, u)| {
                let indices = structure.expanded_index(i).unwrap();
                let permuted_indices: Vec<ConcreteIndex> =
                    permutation.iter().map(|&index| indices[index]).collect();
                let self_value = self.get(&permuted_indices).unwrap();
                self_value.add_fallible(u)
            })
            .collect();

        if let Some(data) = data {
            Some(DenseTensor { structure, data })
        } else {
            None
        }
    }
}

impl<'a, T, U, I, Out> FallibleAdd<&SparseTensor<T, I>> for &'a DenseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleAdd<&'c T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn add_fallible(self, rhs: &SparseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let permutation = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = rhs.structure().clone();

        let data: Option<Vec<Out>> = rhs
            .iter()
            .map(|(indices, value)| {
                let permuted_indices: Vec<ConcreteIndex> =
                    permutation.iter().map(|&index| indices[index]).collect();
                let self_value = self.get(&permuted_indices).unwrap();
                self_value.add_fallible(value)
            })
            .collect();

        if let Some(data) = data {
            Some(DenseTensor { structure, data })
        } else {
            None
        }
    }
}

impl<'a, T, U, I, Out> FallibleAdd<&DenseTensor<T, I>> for &'a SparseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleAdd<&'c T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn add_fallible(self, rhs: &DenseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let permutation = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = rhs.structure().clone();

        let data: Option<Vec<Out>> = rhs
            .iter()
            .map(|(indices, value)| {
                let permuted_indices: Vec<ConcreteIndex> =
                    permutation.iter().map(|&index| indices[index]).collect();
                let self_value = self.get(&permuted_indices).unwrap();
                self_value.add_fallible(value)
            })
            .collect();

        if let Some(data) = data {
            Some(DenseTensor { structure, data })
        } else {
            None
        }
    }
}

impl<'a, T, U, I, Out> FallibleAdd<&SparseTensor<T, I>> for &'a SparseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleAdd<&'c T, Output = Out>,
    I: TensorStructure + Clone,
    U: Default + Clone,
{
    type Output = SparseTensor<Out, I>;
    fn add_fallible(self, rhs: &SparseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let permutation = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = rhs.structure().clone();
        let mut data = SparseTensor::empty(structure);
        for (indices, value) in rhs.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self.smart_get(&permuted_indices).unwrap();
            data.set(&indices, self_value.add_fallible(value)?).unwrap();
        }

        Some(data)
    }
}

impl<'a, T, U, I, Out> FallibleSub<&DenseTensor<T, I>> for &'a DenseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleSub<&'c T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn sub_fallible(self, rhs: &DenseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let permutation = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = rhs.structure().clone();

        let data: Option<Vec<Out>> = rhs
            .data
            .iter()
            .enumerate()
            .map(|(i, u)| {
                let indices = structure.expanded_index(i).unwrap();
                let permuted_indices: Vec<ConcreteIndex> =
                    permutation.iter().map(|&index| indices[index]).collect();
                let self_value = self.get(&permuted_indices).unwrap();
                self_value.sub_fallible(u)
            })
            .collect();

        if let Some(data) = data {
            Some(DenseTensor { structure, data })
        } else {
            None
        }
    }
}

impl<'a, T, U, I, Out> FallibleSub<&DenseTensor<T, I>> for &'a SparseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleSub<&'c T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn sub_fallible(self, rhs: &DenseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let permutation = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = rhs.structure().clone();

        let data: Option<Vec<Out>> = rhs
            .data
            .iter()
            .enumerate()
            .map(|(i, u)| {
                let indices = structure.expanded_index(i).unwrap();
                let permuted_indices: Vec<ConcreteIndex> =
                    permutation.iter().map(|&index| indices[index]).collect();
                let self_value = self.get(&permuted_indices).unwrap();
                self_value.sub_fallible(u)
            })
            .collect();

        if let Some(data) = data {
            Some(DenseTensor { structure, data })
        } else {
            None
        }
    }
}

impl<'a, T, U, I, Out> FallibleSub<&SparseTensor<T, I>> for &'a DenseTensor<U, I>
where
    for<'b, 'c> &'b T: FallibleSub<&'c U, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn sub_fallible(self, rhs: &SparseTensor<T, I>) -> Option<Self::Output> {
        rhs.sub_fallible(self)
    }
}

impl<'a, T, U, I, Out> FallibleSub<&SparseTensor<T, I>> for &'a SparseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleSub<&'c T, Output = Out>,
    I: TensorStructure + Clone,
    U: Default + Clone,
{
    type Output = SparseTensor<Out, I>;
    fn sub_fallible(self, rhs: &SparseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let permutation = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = rhs.structure().clone();
        let mut data = SparseTensor::empty(structure);
        for (indices, value) in rhs.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self.smart_get(&permuted_indices).unwrap();
            data.set(&indices, self_value.sub_fallible(value)?).unwrap();
        }

        Some(data)
    }
}

// impl<'a, T, I> Sub<&DenseTensor<T, I>> for &'a DenseTensor<T, I>
// where
//     for<'b, 'c> &'b T: Sub<&'c T, Output = T>,
//     T: Clone,
//     I: TensorStructure + Clone,
// {
//     type Output = DenseTensor<T, I>;
//     fn sub(self, other: &DenseTensor<T, I>) -> DenseTensor<T, I> {
//         assert!(self.structure().same_external(other.structure()));
//         let permutation = self
//             .structure()
//             .find_permutation(other.structure())
//             .unwrap();

//         let mut result = self.clone();

//         for (indices, value) in other.iter() {
//             let permuted_indices: Vec<ConcreteIndex> =
//                 permutation.iter().map(|&index| indices[index]).collect();
//             let self_value = self.get(&permuted_indices).unwrap();
//             result.set(&indices, self_value - value).unwrap();
//         }
//         result
//     }
// }

// impl<'a, T, I> Sub<&SparseTensor<T, I>> for &'a SparseTensor<T, I>
// where
//     for<'b, 'c> &'b T: Sub<&'c T, Output = T>,
//     T: Clone + Default,
//     I: TensorStructure + Clone,
// {
//     type Output = SparseTensor<T, I>;
//     fn sub(self, other: &SparseTensor<T, I>) -> SparseTensor<T, I> {
//         assert!(self.structure().same_external(other.structure()));
//         let permutation = self
//             .structure()
//             .find_permutation(other.structure())
//             .unwrap();

//         let mut result = self.clone();

//         for (indices, value) in other.iter() {
//             let permuted_indices: Vec<ConcreteIndex> =
//                 permutation.iter().map(|&index| indices[index]).collect();
//             let self_value = self.smart_get(&permuted_indices).unwrap();
//             result.set(&indices, self_value.as_ref() - value).unwrap();
//         }

//         result
//     }
// }

// impl<'a, T, I> Sub<&SparseTensor<T, I>> for &'a DenseTensor<T, I>
// where
//     for<'b, 'c> &'b T: Sub<&'c T, Output = T>,
//     T: Clone,
//     I: TensorStructure + Clone,
// {
//     type Output = DenseTensor<T, I>;
//     fn sub(self, other: &SparseTensor<T, I>) -> DenseTensor<T, I> {
//         assert!(self.structure().same_external(other.structure()));
//         let permutation = self
//             .structure()
//             .find_permutation(other.structure())
//             .unwrap();

//         let mut result = self.clone();

//         for (indices, value) in other.iter() {
//             let permuted_indices: Vec<ConcreteIndex> =
//                 permutation.iter().map(|&index| indices[index]).collect();
//             let self_value = self.get(&permuted_indices).unwrap();
//             result.set(&indices, self_value - value).unwrap();
//         }
//         result
//     }
// }

// impl<'a, T, I> Sub<&DenseTensor<T, I>> for &'a SparseTensor<T, I>
// where
//     for<'b, 'c> &'b T: Sub<&'c T, Output = T>,
//     T: Clone,
//     I: TensorStructure + Clone,
// {
//     type Output = DenseTensor<T, I>;
//     fn sub(self, other: &DenseTensor<T, I>) -> DenseTensor<T, I> {
//         other - self
//     }
// }

// impl<'a, T, I> Mul<&T> for &'a DenseTensor<T, I>
// where
//     T: Num + Clone + Copy,
//     I: TensorStructure + Clone,
// {
//     type Output = DenseTensor<T, I>;
//     fn mul(self, other: &T) -> DenseTensor<T, I> {
//         let mut result = self.clone();

//         for (indices, value) in self.iter() {
//             result.set(&indices, *other * *value).unwrap();
//         }
//         result
//     }
// }

// impl<'a, T, I> Mul<&T> for &'a SparseTensor<T, I>
// where
//     T: Num + Copy + Default,
//     I: TensorStructure + Clone,
// {
//     type Output = SparseTensor<T, I>;
//     fn mul(self, other: &T) -> Self::Output {
//         let mut result = self.clone();

//         for (indices, value) in self.iter() {
//             result.set(&indices, *other * *value).unwrap();
//         }
//         result
//     }
// }

impl<'a, T, U, Out, I> FallibleAdd<&DataTensor<T, I>> for &'a DataTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleAdd<&'c T, Output = Out>,

    U: Default + Clone,
    I: TensorStructure + Clone,
{
    type Output = DataTensor<Out, I>;
    fn add_fallible(self, rhs: &DataTensor<T, I>) -> Option<Self::Output> {
        match (self, rhs) {
            (DataTensor::Dense(a), DataTensor::Dense(b)) => {
                Some(DataTensor::Dense(a.add_fallible(b)?))
            }
            (DataTensor::Sparse(a), DataTensor::Sparse(b)) => {
                Some(DataTensor::Sparse(a.add_fallible(b)?))
            }
            (DataTensor::Dense(a), DataTensor::Sparse(b)) => {
                Some(DataTensor::Dense(a.add_fallible(b)?))
            }
            (DataTensor::Sparse(a), DataTensor::Dense(b)) => {
                Some(DataTensor::Dense(a.add_fallible(b)?))
            }
        }
    }
}
