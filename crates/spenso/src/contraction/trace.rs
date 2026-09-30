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
        // Trace positions refer to the structure being traced, so each pair is
        // located only after the previous trace has been applied.
        while let Some(&trace) = result.traces().first() {
            let mut new_structure = result.structure.clone();
            new_structure.trace(trace[0], trace[1]);

            let mut new_result =
                DenseTensor::from_storage_data_coerced(&result.data, new_structure)
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

#[cfg(test)]
mod tests {
    use crate::{
        contraction::Trace,
        structure::{
            Canonicalized, HasStructure, OrderedStructure,
            representation::{Euclidean, LibraryRep, Lorentz, RepName},
            slot::{DualSlotTo, IsAbstractSlot},
        },
        tensors::data::DenseTensor,
    };

    #[test]
    fn dense_trace_contracts_every_pair() {
        // Pairs of different dimensions make each trace depend on the structure
        // left behind by the previous one.
        let euclidean = Euclidean {}.new_slot(2, 0).to_lib();
        let lorentz = Lorentz {}.new_slot(3, 1).to_lib();
        let structure: OrderedStructure<LibraryRep> =
            Canonicalized::from_iter([euclidean, euclidean, lorentz, lorentz.dual()])
                .into_canonical();
        let tensor =
            DenseTensor::from_storage_data((1..=36).map(f64::from).collect(), structure).unwrap();

        // The row-major diagonal components T[e, e, l, l] are 1 + 27 e + 4 l.
        let expected = (0..2)
            .flat_map(|e| (0..3).map(move |l| f64::from(1 + 27 * e + 4 * l)))
            .sum::<f64>();
        assert_eq!(tensor.internal_contract().scalar(), Some(expected));
        assert_eq!(
            tensor.to_sparse().internal_contract().scalar(),
            Some(expected)
        );
    }
}
