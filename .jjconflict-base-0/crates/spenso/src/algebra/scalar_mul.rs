use super::ScalarMul;
use crate::{
    algebra::complex::{Complex, RealOrComplex},
    algebra::upgrading_arithmetic::FallibleMul,
    iterators::IteratableTensor,
    structure::{HasStructure, TensorStructure},
    tensors::complex::RealOrComplexTensor,
    tensors::data::{DataTensor, DenseTensor, SetTensorData, SparseTensor},
};

impl<T, U, I, Out> ScalarMul<T> for DenseTensor<U, I>
where
    U: FallibleMul<T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = DenseTensor<Out, I>;
    fn scalar_mul(&self, rhs: &T) -> Option<Self::Output> {
        let data: Option<Vec<Out>> = self.iter_flat().map(|(_, u)| u.mul_fallible(rhs)).collect();

        data.map(|data| DenseTensor {
            structure: self.structure().clone(),
            data,
        })
    }
}

impl<T, U, I, Out> ScalarMul<T> for SparseTensor<U, I>
where
    U: FallibleMul<T, Output = Out>,
    I: TensorStructure + Clone,
{
    type Output = SparseTensor<Out, I>;
    fn scalar_mul(&self, rhs: &T) -> Option<Self::Output> {
        let mut data = SparseTensor::empty(self.structure().clone(), self.zero.mul_fallible(rhs)?);
        for (indices, u) in self.iter_flat() {
            data.set_flat(indices, u.mul_fallible(rhs)?).unwrap();
        }
        Some(data)
    }
}

impl<T, U, I, Out> ScalarMul<T> for DataTensor<U, I>
where
    DenseTensor<U, I>: ScalarMul<T, Output = DenseTensor<Out, I>>,
    SparseTensor<U, I>: ScalarMul<T, Output = SparseTensor<Out, I>>,
    I: TensorStructure + Clone,
{
    type Output = DataTensor<Out, I>;
    fn scalar_mul(&self, rhs: &T) -> Option<Self::Output> {
        match self {
            DataTensor::Dense(a) => Some(DataTensor::Dense(a.scalar_mul(rhs)?)),
            DataTensor::Sparse(a) => Some(DataTensor::Sparse(a.scalar_mul(rhs)?)),
        }
    }
}

impl<T, U, S, Out> ScalarMul<RealOrComplex<T>> for RealOrComplexTensor<U, S>
where
    DataTensor<U, S>: ScalarMul<T, Output = DataTensor<Out, S>>
        + ScalarMul<Complex<T>, Output = DataTensor<Complex<Out>, S>>,
    DataTensor<Complex<U>, S>: ScalarMul<T, Output = DataTensor<Complex<Out>, S>>
        + ScalarMul<Complex<T>, Output = DataTensor<Complex<Out>, S>>,
    S: TensorStructure + Clone,
{
    type Output = RealOrComplexTensor<Out, S>;

    fn scalar_mul(&self, rhs: &RealOrComplex<T>) -> Option<Self::Output> {
        Some(match (self, rhs) {
            (RealOrComplexTensor::Real(rt), RealOrComplex::Real(rs)) => {
                RealOrComplexTensor::Real(rt.scalar_mul(rs)?)
            }
            (RealOrComplexTensor::Complex(rt), RealOrComplex::Real(rs)) => {
                RealOrComplexTensor::Complex(rt.scalar_mul(rs)?)
            }
            (RealOrComplexTensor::Real(rt), RealOrComplex::Complex(rs)) => {
                RealOrComplexTensor::Complex(rt.scalar_mul(rs)?)
            }
            (RealOrComplexTensor::Complex(rt), RealOrComplex::Complex(rs)) => {
                RealOrComplexTensor::Complex(rt.scalar_mul(rs)?)
            }
        })
    }
}

impl<U, S> ScalarMul<U> for RealOrComplexTensor<U, S>
where
    DataTensor<U, S>: ScalarMul<U, Output = DataTensor<U, S>>,
    DataTensor<Complex<U>, S>: ScalarMul<U, Output = DataTensor<Complex<U>, S>>,
    S: TensorStructure + Clone,
{
    type Output = RealOrComplexTensor<U, S>;

    fn scalar_mul(&self, rhs: &U) -> Option<Self::Output> {
        Some(match self {
            RealOrComplexTensor::Real(rt) => RealOrComplexTensor::Real(rt.scalar_mul(rhs)?),
            RealOrComplexTensor::Complex(rt) => RealOrComplexTensor::Complex(rt.scalar_mul(rhs)?),
        })
    }
}

impl<U, S> ScalarMul<Complex<U>> for RealOrComplexTensor<U, S>
where
    DataTensor<U, S>: ScalarMul<Complex<U>, Output = DataTensor<Complex<U>, S>>,
    DataTensor<Complex<U>, S>: ScalarMul<Complex<U>, Output = DataTensor<Complex<U>, S>>,
    S: TensorStructure + Clone,
{
    type Output = RealOrComplexTensor<U, S>;

    fn scalar_mul(&self, rhs: &Complex<U>) -> Option<Self::Output> {
        Some(match self {
            RealOrComplexTensor::Real(rt) => RealOrComplexTensor::Complex(rt.scalar_mul(rhs)?),
            RealOrComplexTensor::Complex(rt) => RealOrComplexTensor::Complex(rt.scalar_mul(rhs)?),
        })
    }
}
