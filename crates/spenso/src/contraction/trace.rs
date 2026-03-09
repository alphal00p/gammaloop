use crate::{
    algebra::algebraic_traits::{IsZero, RefZero},
    algebra::upgrading_arithmetic::{FallibleAddAssign, FallibleSubAssign},
    structure::{StructureContract, TensorStructure},
    tensors::data::{DataTensor, DenseTensor, SetTensorData, SparseTensor},
};

use std::iter::Iterator;

use super::{ContractableWith, Trace};

impl<T, I> Trace for DenseTensor<T, I>
where
    T: ContractableWith<T, Out = T> + Clone + RefZero + FallibleAddAssign<T> + FallibleSubAssign<T>,
    I: TensorStructure + Clone + StructureContract,
{
    /// Contract the tensor with itself, i.e. trace over all matching indices.
    fn internal_contract(&self) -> Self {
        let mut result: DenseTensor<T, I> = self.clone();
        for trace in self.traces() {
            let mut new_structure = self.structure.clone();
            new_structure.trace(trace[0], trace[1]);

            let mut new_result = DenseTensor::from_data_coerced(&self.data, new_structure)
                .unwrap_or_else(|_| unreachable!());
            for (idx, t) in result.iter_trace(trace) {
                new_result.set(&idx, t).unwrap_or_else(|_| unreachable!());
            }
            result = new_result;
        }
        result
    }
}

impl<T, I> Trace for SparseTensor<T, I>
where
    T: ContractableWith<T, Out = T>
        + Clone
        + RefZero
        + IsZero
        + FallibleAddAssign<T>
        + FallibleSubAssign<T>,
    I: TensorStructure + Clone + StructureContract,
{
    /// Contract the tensor with itself, i.e. trace over all matching indices.
    fn internal_contract(&self) -> Self {
        let trace = if let Some(e) = self.traces().first() {
            *e
        } else {
            return self.clone();
        };

        // println!("trace {:?}", trace);
        let mut new_structure = self.structure.clone();
        // println!("{}", new_structure);
        new_structure.trace(trace[0], trace[1]);

        let mut new_result = SparseTensor::empty(new_structure, self.zero.clone());
        for (idx, t) in self.iter_trace(trace).filter(|(_, t)| !t.is_zero()) {
            new_result.set(&idx, t).unwrap();
        }

        if new_result.traces().is_empty() {
            new_result
        } else {
            new_result.internal_contract()
        }
    }
}

impl<T, I> Trace for DataTensor<T, I>
where
    T: ContractableWith<T, Out = T>,
    T: FallibleAddAssign<T> + FallibleSubAssign<T> + Clone + RefZero + IsZero,
    I: TensorStructure + Clone + StructureContract,
{
    fn internal_contract(&self) -> Self {
        match self {
            DataTensor::Dense(d) => DataTensor::Dense(d.internal_contract()),
            DataTensor::Sparse(s) => DataTensor::Sparse(s.internal_contract()),
        }
    }
}
