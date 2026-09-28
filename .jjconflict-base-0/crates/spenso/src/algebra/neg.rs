use crate::{
    algebra::complex::Complex,
    structure::TensorStructure,
    tensors::complex::RealOrComplexTensor,
    tensors::data::{DataTensor, DenseTensor, SparseTensor},
};

impl<T, U, S> std::ops::Neg for DenseTensor<T, S>
where
    T: std::ops::Neg<Output = U>,
    S: TensorStructure + Clone,
{
    type Output = DenseTensor<U, S>;
    fn neg(self) -> Self::Output {
        DenseTensor {
            structure: self.structure.clone(),
            data: self.data.into_iter().map(|x| x.neg()).collect(),
        }
    }
}

impl<T, U, S> std::ops::Neg for SparseTensor<T, S>
where
    T: std::ops::Neg<Output = U>,
    S: TensorStructure + Clone,
{
    type Output = SparseTensor<U, S>;
    fn neg(self) -> Self::Output {
        SparseTensor {
            zero: self.zero.neg(),
            structure: self.structure.clone(),
            elements: self
                .elements
                .into_iter()
                .map(|x| (x.0, x.1.neg()))
                .collect(),
        }
    }
}

impl<T, U, S> std::ops::Neg for DataTensor<T, S>
where
    T: std::ops::Neg<Output = U>,
    S: TensorStructure + Clone,
{
    type Output = DataTensor<U, S>;
    fn neg(self) -> Self::Output {
        match self {
            DataTensor::Dense(d) => DataTensor::Dense(-d),
            DataTensor::Sparse(s) => DataTensor::Sparse(-s),
        }
    }
}

impl<T, U, S> std::ops::Neg for RealOrComplexTensor<T, S>
where
    T: std::ops::Neg<Output = U>,
    Complex<T>: std::ops::Neg<Output = Complex<U>>,
    S: TensorStructure + Clone,
{
    type Output = RealOrComplexTensor<U, S>;
    fn neg(self) -> Self::Output {
        match self {
            RealOrComplexTensor::Real(d) => RealOrComplexTensor::Real(-d),
            RealOrComplexTensor::Complex(s) => RealOrComplexTensor::Complex(-s),
        }
    }
}
