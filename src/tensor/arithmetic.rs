use crate::tensor::{ConcreteIndex, GetTensorData, SetTensorData};

use super::{DenseTensor, SparseTensor, TensorStructure};
use num::traits::Num;
use std::ops::{Add, Mul, Sub};

impl<T, I> Add<DenseTensor<T, I>> for DenseTensor<T, I>
where
    T: Num + Clone + Copy,

    I: TensorStructure + Clone,
{
    type Output = DenseTensor<T, I>;
    fn add(self, other: DenseTensor<T, I>) -> DenseTensor<T, I> {
        assert!(self.structure().same_external(other.structure()));

        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self.get(&permuted_indices).unwrap();
            let _ = result.set(&indices, *self_value + *value);
        }
        result
    }
}

impl<T, I> Add<SparseTensor<T, I>> for SparseTensor<T, I>
where
    T: Num + Clone + Default + Copy,
    I: TensorStructure + Clone,
{
    type Output = SparseTensor<T, I>;
    fn add(self, other: SparseTensor<T, I>) -> SparseTensor<T, I> {
        assert!(self.structure().same_external(other.structure()));
        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self.smart_get(&permuted_indices).unwrap().into_owned();
            result.set(&indices, self_value + *value).unwrap();
        }

        result
    }
}

impl<T, I> Add<SparseTensor<T, I>> for DenseTensor<T, I>
where
    T: Num + Clone + Default + Copy,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<T, I>;
    fn add(self, other: SparseTensor<T, I>) -> DenseTensor<T, I> {
        assert!(self.structure().same_external(other.structure()));

        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self.get(&permuted_indices).unwrap();
            result.set(&indices, *self_value + *value).unwrap();
        }
        result
    }
}

impl<T, I> Add<DenseTensor<T, I>> for SparseTensor<T, I>
where
    T: Num + Clone + Default + Copy,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<T, I>;
    fn add(self, other: DenseTensor<T, I>) -> DenseTensor<T, I> {
        other + self
    }
}

impl<T, I> Sub<DenseTensor<T, I>> for DenseTensor<T, I>
where
    T: Num + Clone + Copy,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<T, I>;
    fn sub(self, other: DenseTensor<T, I>) -> DenseTensor<T, I> {
        assert!(self.structure().same_external(other.structure()));
        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self.get(&permuted_indices).unwrap();
            result.set(&indices, *self_value - *value).unwrap();
        }
        result
    }
}

impl<T, I> Sub<SparseTensor<T, I>> for SparseTensor<T, I>
where
    T: Num + Default + Copy,
    I: TensorStructure + Clone,
{
    type Output = SparseTensor<T, I>;
    fn sub(self, other: SparseTensor<T, I>) -> SparseTensor<T, I> {
        assert!(self.structure().same_external(other.structure()));
        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self.smart_get(&permuted_indices).unwrap().into_owned();
            result.set(&indices, self_value - *value).unwrap();
        }

        result
    }
}

impl<T, I> Sub<SparseTensor<T, I>> for DenseTensor<T, I>
where
    T: Num + Clone + Default + Copy,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<T, I>;
    fn sub(self, other: SparseTensor<T, I>) -> DenseTensor<T, I> {
        assert!(self.structure().same_external(other.structure()));
        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self.get(&permuted_indices).unwrap();
            result.set(&indices, *self_value - *value).unwrap();
        }
        result
    }
}

impl<T, I> Sub<DenseTensor<T, I>> for SparseTensor<T, I>
where
    T: Num + Clone + Default + Copy,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<T, I>;
    fn sub(self, other: DenseTensor<T, I>) -> DenseTensor<T, I> {
        other - self
    }
}

impl<T, I> Mul<T> for DenseTensor<T, I>
where
    T: Num + Clone + Copy,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<T, I>;
    fn mul(self, other: T) -> DenseTensor<T, I> {
        let mut result = self.clone();

        for (indices, value) in self.iter() {
            result.set(&indices, other * *value).unwrap();
        }
        result
    }
}

impl<T, I> Mul<T> for SparseTensor<T, I>
where
    T: Num + Copy + Default,
    I: TensorStructure + Clone,
{
    type Output = SparseTensor<T, I>;
    fn mul(self, other: T) -> Self::Output {
        let mut result = self.clone();

        for (indices, value) in self.iter() {
            result.set(&indices, other * *value).unwrap();
        }
        result
    }
}
