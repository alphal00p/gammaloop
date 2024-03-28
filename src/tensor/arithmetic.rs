use std::fmt::Debug;

use crate::tensor::{ConcreteIndex, GetTensorData, SetTensorData};

use super::{
    DataTensor, DenseTensor, FallibleAdd, FallibleMul, FallibleSub, MixedTensor, SparseTensor,
    TensorStructure, TryFromUpgrade, TryIntoUpgrade,
};

use symbolica::domains::float::Complex;
use symbolica::representations::Atom;

impl<'a, T, U, I, Out> FallibleAdd<&DenseTensor<T, I>> for &'a DenseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleAdd<&'c T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn add_fallible(self, rhs: &DenseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        // Makes rhs into self ,when applied.
        let permutation = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = self.structure().clone();

        let data: Option<Vec<Out>> = self
            .iter()
            .map(|(indices, u)| {
                let permuted_indices: Vec<ConcreteIndex> =
                    permutation.iter().map(|&index| indices[index]).collect();
                let t = rhs.get(&permuted_indices).unwrap();
                u.add_fallible(t)
            })
            .collect();

        data.map(|data| DenseTensor { structure, data })
    }
}

impl<'a, T, U, I, Out> FallibleAdd<&SparseTensor<T, I>> for &'a DenseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleAdd<&'c T, Output = Out>,
    I: TensorStructure + Clone,
    T: Default + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn add_fallible(self, rhs: &SparseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let permutation = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = self.structure().clone();

        let data: Option<Vec<Out>> = self
            .iter()
            .map(|(indices, u)| {
                let permuted_indices: Vec<ConcreteIndex> =
                    permutation.iter().map(|&index| indices[index]).collect();
                let t = rhs.smart_get(&permuted_indices).unwrap();
                u.add_fallible(&t)
            })
            .collect();

        data.map(|data| DenseTensor { structure, data })
    }
}

impl<'a, T, U, I, Out> FallibleAdd<&DenseTensor<T, I>> for &'a SparseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleAdd<&'c T, Output = Out>,
    I: TensorStructure + Clone,
    U: Default + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn add_fallible(self, rhs: &DenseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let permutation = rhs.structure().find_permutation(self.structure()).unwrap();
        let structure = rhs.structure().clone();

        let data: Option<Vec<Out>> = rhs
            .iter()
            .map(|(indices, t)| {
                let permuted_indices: Vec<ConcreteIndex> =
                    permutation.iter().map(|&index| indices[index]).collect();
                let u = self.smart_get(&permuted_indices).unwrap();
                u.add_fallible(t)
            })
            .collect();

        data.map(|data| DenseTensor { structure, data })
    }
}

impl<'a, T, U, I, Out> FallibleAdd<&SparseTensor<T, I>> for &'a SparseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleAdd<&'c T, Output = Out>,
    I: TensorStructure + Clone,
    T: Default + Clone + Debug,
    Out: Default + PartialEq + TryFromUpgrade<T>,
{
    type Output = SparseTensor<Out, I>;
    fn add_fallible(self, rhs: &SparseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let permutation = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = self.structure().clone();
        let mut data = SparseTensor::empty(structure);
        for (indices, u) in self.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let t = rhs.smart_get(&permuted_indices).unwrap();
            println!("{:?}", t);
            data.smart_set(&indices, u.add_fallible(&t)?).unwrap();
        }

        let permutation: Vec<usize> = rhs.structure().find_permutation(self.structure()).unwrap();
        for (i, t) in rhs.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| i[index]).collect();

            if self.get(&permuted_indices).is_err() {
                data.smart_set(&i, t.clone().try_into_upgrade()?).unwrap();
            }
        }

        Some(data)
    }
}

impl<'a, T, U, Out, I> FallibleAdd<&DataTensor<T, I>> for &'a DataTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleAdd<&'c T, Output = Out>,
    U: Default + Clone,
    T: Default + Clone + Debug,
    Out: Default + PartialEq + TryFromUpgrade<T>,
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

impl<'a, I> FallibleAdd<&MixedTensor<I>> for &'a MixedTensor<I>
where
    I: TensorStructure + Clone,
{
    type Output = MixedTensor<I>;
    fn add_fallible(self, rhs: &MixedTensor<I>) -> Option<Self::Output> {
        match (self, rhs) {
            (MixedTensor::Float(a), MixedTensor::Float(b)) => {
                Some(MixedTensor::Float(a.add_fallible(b)?))
            }
            (MixedTensor::Complex(a), MixedTensor::Complex(b)) => {
                Some(MixedTensor::Complex(a.add_fallible(b)?))
            }
            (MixedTensor::Float(a), MixedTensor::Complex(b)) => {
                Some(MixedTensor::Complex(a.add_fallible(b)?))
            }
            (MixedTensor::Complex(a), MixedTensor::Float(b)) => {
                Some(MixedTensor::Complex(a.add_fallible(b)?))
            }
            (MixedTensor::Symbolic(a), MixedTensor::Symbolic(b)) => {
                Some(MixedTensor::Symbolic(a.add_fallible(b)?))
            }
            (MixedTensor::Symbolic(a), MixedTensor::Float(b)) => {
                Some(MixedTensor::Symbolic(a.add_fallible(b)?))
            }
            (MixedTensor::Symbolic(a), MixedTensor::Complex(b)) => {
                Some(MixedTensor::Symbolic(a.add_fallible(b)?))
            }
            (MixedTensor::Float(a), MixedTensor::Symbolic(b)) => {
                Some(MixedTensor::Symbolic(a.add_fallible(b)?))
            }
            (MixedTensor::Complex(a), MixedTensor::Symbolic(b)) => {
                Some(MixedTensor::Symbolic(a.add_fallible(b)?))
            }
        }
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
        let structure = self.structure().clone();

        let data: Option<Vec<Out>> = self
            .iter()
            .map(|(indices, u)| {
                let permuted_indices: Vec<ConcreteIndex> =
                    permutation.iter().map(|&index| indices[index]).collect();
                let t = rhs.get(&permuted_indices).unwrap();
                u.sub_fallible(t)
            })
            .collect();

        data.map(|data| DenseTensor { structure, data })
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
        let permutation = rhs.structure().find_permutation(self.structure()).unwrap();
        let structure = rhs.structure().clone();

        let data: Option<Vec<Out>> = rhs
            .iter()
            .map(|(indices, t)| {
                let permuted_indices: Vec<ConcreteIndex> =
                    permutation.iter().map(|&index| indices[index]).collect();
                let u = self.get(&permuted_indices).unwrap();
                u.sub_fallible(t)
            })
            .collect();

        data.map(|data| DenseTensor { structure, data })
    }
}

impl<'a, T, U, I, Out> FallibleSub<&SparseTensor<T, I>> for &'a DenseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleSub<&'c T, Output = Out>,
    I: TensorStructure + Clone,
    T: Default + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn sub_fallible(self, rhs: &SparseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let permutation = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = self.structure().clone();

        let data: Option<Vec<Out>> = self
            .iter()
            .map(|(indices, u)| {
                let permuted_indices: Vec<ConcreteIndex> =
                    permutation.iter().map(|&index| indices[index]).collect();
                let t = rhs.smart_get(&permuted_indices).unwrap();
                u.sub_fallible(&t)
            })
            .collect();

        data.map(|data| DenseTensor { structure, data })
    }
}

impl<'a, T, U, I, Out> FallibleSub<&SparseTensor<T, I>> for &'a SparseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleSub<&'c T, Output = Out>,
    I: TensorStructure + Clone,
    T: Default + Clone,
    U: Default,
    Out: Default + PartialEq,
{
    type Output = SparseTensor<Out, I>;
    fn sub_fallible(self, rhs: &SparseTensor<T, I>) -> Option<Self::Output> {
        assert!(self.structure().same_external(rhs.structure()));
        let permutation = self.structure().find_permutation(rhs.structure()).unwrap();
        let structure = self.structure().clone();
        let mut data = SparseTensor::empty(structure);
        for (indices, u) in self.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| indices[index]).collect();
            let t = rhs.smart_get(&permuted_indices).unwrap();
            data.smart_set(&indices, u.sub_fallible(&t)?).unwrap();
        }
        let permutation: Vec<usize> = rhs.structure().find_permutation(self.structure()).unwrap();
        for (i, t) in rhs.iter() {
            let permuted_indices: Vec<ConcreteIndex> =
                permutation.iter().map(|&index| i[index]).collect();

            if self.get(&permuted_indices).is_err() {
                data.smart_set(&i, U::default().sub_fallible(t)?).unwrap();
            }
        }

        Some(data)
    }
}

impl<'a, T, U, Out, I> FallibleSub<&DataTensor<T, I>> for &'a DataTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleSub<&'c T, Output = Out>,
    U: Default + Clone,
    T: Default + Clone,
    Out: Default + PartialEq,
    I: TensorStructure + Clone,
{
    type Output = DataTensor<Out, I>;
    fn sub_fallible(self, rhs: &DataTensor<T, I>) -> Option<Self::Output> {
        match (self, rhs) {
            (DataTensor::Dense(a), DataTensor::Dense(b)) => {
                Some(DataTensor::Dense(a.sub_fallible(b)?))
            }
            (DataTensor::Sparse(a), DataTensor::Sparse(b)) => {
                Some(DataTensor::Sparse(a.sub_fallible(b)?))
            }
            (DataTensor::Dense(a), DataTensor::Sparse(b)) => {
                Some(DataTensor::Dense(a.sub_fallible(b)?))
            }
            (DataTensor::Sparse(a), DataTensor::Dense(b)) => {
                Some(DataTensor::Dense(a.sub_fallible(b)?))
            }
        }
    }
}

impl<'a, I> FallibleSub<&MixedTensor<I>> for &'a MixedTensor<I>
where
    I: TensorStructure + Clone,
{
    type Output = MixedTensor<I>;
    fn sub_fallible(self, rhs: &MixedTensor<I>) -> Option<Self::Output> {
        match (self, rhs) {
            (MixedTensor::Float(a), MixedTensor::Float(b)) => {
                Some(MixedTensor::Float(a.sub_fallible(b)?))
            }
            (MixedTensor::Complex(a), MixedTensor::Complex(b)) => {
                Some(MixedTensor::Complex(a.sub_fallible(b)?))
            }
            (MixedTensor::Float(a), MixedTensor::Complex(b)) => {
                Some(MixedTensor::Complex(a.sub_fallible(b)?))
            }
            (MixedTensor::Complex(a), MixedTensor::Float(b)) => {
                Some(MixedTensor::Complex(a.sub_fallible(b)?))
            }
            (MixedTensor::Symbolic(a), MixedTensor::Symbolic(b)) => {
                Some(MixedTensor::Symbolic(a.sub_fallible(b)?))
            }
            (MixedTensor::Symbolic(a), MixedTensor::Float(b)) => {
                Some(MixedTensor::Symbolic(a.sub_fallible(b)?))
            }
            (MixedTensor::Symbolic(a), MixedTensor::Complex(b)) => {
                Some(MixedTensor::Symbolic(a.sub_fallible(b)?))
            }
            (MixedTensor::Float(a), MixedTensor::Symbolic(b)) => {
                Some(MixedTensor::Symbolic(a.sub_fallible(b)?))
            }
            (MixedTensor::Complex(a), MixedTensor::Symbolic(b)) => {
                Some(MixedTensor::Symbolic(a.sub_fallible(b)?))
            }
        }
    }
}

pub trait ScalarMul<T> {
    type Output;
    fn scalar_mul(self, rhs: T) -> Option<Self::Output>;
}

impl<'a, T, U, I, Out> ScalarMul<&T> for &'a DenseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleMul<&'c T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn scalar_mul(self, rhs: &T) -> Option<Self::Output> {
        let data: Option<Vec<Out>> = self.iter_flat().map(|(_, u)| u.mul_fallible(rhs)).collect();

        data.map(|data| DenseTensor {
            structure: self.structure().clone(),
            data,
        })
    }
}

impl<'a, T, U, I, Out> ScalarMul<&T> for &'a SparseTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleMul<&'c T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = SparseTensor<Out, I>;
    fn scalar_mul(self, rhs: &T) -> Option<Self::Output> {
        let mut data = SparseTensor::empty(self.structure().clone());
        for (indices, u) in self.iter_flat() {
            data.set_flat(indices, u.mul_fallible(rhs)?).unwrap();
        }
        Some(data)
    }
}

impl<'a, T, U, I, Out> ScalarMul<&T> for &'a DataTensor<U, I>
where
    for<'b, 'c> &'b U: FallibleMul<&'c T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DataTensor<Out, I>;
    fn scalar_mul(self, rhs: &T) -> Option<Self::Output> {
        match self {
            DataTensor::Dense(a) => Some(DataTensor::Dense(a.scalar_mul(rhs)?)),
            DataTensor::Sparse(a) => Some(DataTensor::Sparse(a.scalar_mul(rhs)?)),
        }
    }
}

impl<'a, I> ScalarMul<&f64> for &'a MixedTensor<I>
where
    I: TensorStructure + Clone,
{
    type Output = MixedTensor<I>;
    fn scalar_mul(self, rhs: &f64) -> Option<Self::Output> {
        match self {
            MixedTensor::Float(a) => Some(MixedTensor::Float(a.scalar_mul(rhs)?)),
            MixedTensor::Complex(a) => Some(MixedTensor::Complex(a.scalar_mul(rhs)?)),
            MixedTensor::Symbolic(a) => Some(MixedTensor::Symbolic(a.scalar_mul(rhs)?)),
        }
    }
}

impl<'a, I> ScalarMul<&Complex<f64>> for &'a MixedTensor<I>
where
    I: TensorStructure + Clone,
{
    type Output = MixedTensor<I>;
    fn scalar_mul(self, rhs: &Complex<f64>) -> Option<Self::Output> {
        match self {
            MixedTensor::Float(a) => Some(MixedTensor::Complex(a.scalar_mul(rhs)?)),
            MixedTensor::Complex(a) => Some(MixedTensor::Complex(a.scalar_mul(rhs)?)),
            MixedTensor::Symbolic(a) => Some(MixedTensor::Symbolic(a.scalar_mul(rhs)?)),
        }
    }
}

impl<'a, I> ScalarMul<&Atom> for &'a MixedTensor<I>
where
    I: TensorStructure + Clone,
{
    type Output = MixedTensor<I>;
    fn scalar_mul(self, rhs: &Atom) -> Option<Self::Output> {
        match self {
            MixedTensor::Float(a) => Some(MixedTensor::Symbolic(a.scalar_mul(rhs)?)),
            MixedTensor::Complex(a) => Some(MixedTensor::Symbolic(a.scalar_mul(rhs)?)),
            MixedTensor::Symbolic(a) => Some(MixedTensor::Symbolic(a.scalar_mul(rhs)?)),
        }
    }
}
