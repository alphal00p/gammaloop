use crate::tensor::ConcreteIndex;

use super::{DenseTensor, HasTensorStructure, SparseTensor, VecSlotExtension};
use num::traits::Num;
use std::ops::{Add, Mul, Sub};

impl<T> Add<DenseTensor<T>> for DenseTensor<T>
where
    T: Num + Clone + Copy,
{
    type Output = DenseTensor<T>;
    fn add(self, other: DenseTensor<T>) -> DenseTensor<T> {
        assert!(self.structure().same_content(other.structure()));

        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self.get(&permuted_indices).unwrap();
            result.set(&indices, *self_value + *value);
        }
        result
    }
}

impl<T> Add<SparseTensor<T>> for SparseTensor<T>
where
    T: Num + Clone + Default + Copy,
{
    type Output = SparseTensor<T>;
    fn add(self, other: SparseTensor<T>) -> SparseTensor<T> {
        assert!(self.structure().same_content(other.structure()));
        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self
                .get_with_defaults(&permuted_indices)
                .unwrap()
                .into_owned();
            result.set(indices, self_value + *value).unwrap();
        }

        result
    }
}

impl<T> Add<SparseTensor<T>> for DenseTensor<T>
where
    T: Num + Clone + Default + Copy,
{
    type Output = DenseTensor<T>;
    fn add(self, other: SparseTensor<T>) -> DenseTensor<T> {
        assert!(self.structure().same_content(other.structure()));

        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self.get(&permuted_indices).unwrap();
            result.set(indices, *self_value + *value);
        }
        result
    }
}

impl<T> Add<DenseTensor<T>> for SparseTensor<T>
where
    T: Num + Clone + Default + Copy,
{
    type Output = DenseTensor<T>;
    fn add(self, other: DenseTensor<T>) -> DenseTensor<T> {
        other + self
    }
}

impl<T> Sub<DenseTensor<T>> for DenseTensor<T>
where
    T: Num + Clone + Copy,
{
    type Output = DenseTensor<T>;
    fn sub(self, other: DenseTensor<T>) -> DenseTensor<T> {
        assert!(self.structure().same_content(other.structure()));
        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self.get(&permuted_indices).unwrap();
            result.set(&indices, *self_value - *value);
        }
        result
    }
}

impl<T> Sub<SparseTensor<T>> for SparseTensor<T>
where
    T: Num + Default + Copy,
{
    type Output = SparseTensor<T>;
    fn sub(self, other: SparseTensor<T>) -> SparseTensor<T> {
        assert!(self.structure().same_content(other.structure()));
        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self
                .get_with_defaults(&permuted_indices)
                .unwrap()
                .into_owned();
            result.set(indices, self_value - *value).unwrap();
        }

        result
    }
}

impl<T> Sub<SparseTensor<T>> for DenseTensor<T>
where
    T: Num + Clone + Default + Copy,
{
    type Output = DenseTensor<T>;
    fn sub(self, other: SparseTensor<T>) -> DenseTensor<T> {
        assert!(self.structure().same_content(other.structure()));
        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let self_value = self.get(&permuted_indices).unwrap();
            result.set(indices, *self_value - *value);
        }
        result
    }
}

impl<T> Sub<DenseTensor<T>> for SparseTensor<T>
where
    T: Num + Clone + Default + Copy,
{
    type Output = DenseTensor<T>;
    fn sub(self, other: DenseTensor<T>) -> DenseTensor<T> {
        other - self
    }
}

impl<T> Mul<T> for DenseTensor<T>
where
    T: Num + Clone + Copy,
{
    type Output = DenseTensor<T>;
    fn mul(self, other: T) -> DenseTensor<T> {
        let mut result = self.clone();

        for (indices, value) in self.iter() {
            result.set(&indices, other * *value);
        }
        result
    }
}

impl<T> Mul<T> for SparseTensor<T>
where
    T: Num + Copy + Default,
{
    type Output = SparseTensor<T>;
    fn mul(self, other: T) -> Self::Output {
        let mut result = self.clone();

        for (indices, value) in self.iter() {
            result.set(indices, other * *value).unwrap();
        }
        result
    }
}
