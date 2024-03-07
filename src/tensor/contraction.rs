use ahash::AHashMap;

use num::Zero;

use serde::{Deserialize, Serialize};
use slotmap::{new_key_type, DenseSlotMap, Key, SecondaryMap};
use symbolica::state::{State, Workspace};

use self::parametric::{MixedTensor, MixedTensors, SymbolicContract};
use self::structure::HistoryStructure;

use super::{
    parametric, structure, DataIterator, DataTensor, DenseTensor, HasName, HasTensorData,
    NumTensor, Representation, SetTensorData, Shadowable, Slot, SparseTensor, StructureContract,
    TensorStructure, TracksCount,
};
use smartstring::alias::String;
use std::{
    fmt::{Debug, Display},
    // intrinsics::needs_drop,
    ops::Neg,
};

trait LeastCommonStorage<Other: HasTensorData + SetTensorData>: HasTensorData + SetTensorData {
    type OutStorage<LCMData>: SetTensorData<SetData = LCMData>;
    fn least_common_storage<LCMData>(&self, other: &Other) -> Self::OutStorage<LCMData>
    where
        for<'a, 'b> &'a Self::Data: std::ops::Mul<&'b Other::Data, Output = LCMData>,
        LCMData: num::Zero + Clone;

    fn empty<LCMData>(structure: Self::Structure) -> Self::OutStorage<LCMData>
    where
        for<'a, 'b> &'a Self::Data: std::ops::Mul<&'b Other::Data, Output = LCMData>,
        LCMData: Zero + Clone;
}

impl<T, U, I> LeastCommonStorage<DenseTensor<T, I>> for DenseTensor<U, I>
where
    T: Clone,
    U: Clone,
    I: TensorStructure + StructureContract + Clone,
{
    type OutStorage<LCMData> = DenseTensor<LCMData, I>;
    fn least_common_storage<LCMData>(&self, other: &DenseTensor<T, I>) -> Self::OutStorage<LCMData>
    where
        for<'a, 'b> &'a Self::Data: std::ops::Mul<&'b T, Output = LCMData>,
        LCMData: Zero + Clone,
    {
        let mut final_structure = self.structure().clone();
        final_structure.merge(other.structure());
        DenseTensor::zero(final_structure)
    }

    fn empty<LCMData>(structure: Self::Structure) -> Self::OutStorage<LCMData>
    where
        for<'a, 'b> &'a Self::Data:
            std::ops::Mul<&'b <DenseTensor<T, I> as HasTensorData>::Data, Output = LCMData>,
        LCMData: Zero + Clone,
    {
        DenseTensor::zero(structure)
    }
}

impl<T, U, I> LeastCommonStorage<DenseTensor<T, I>> for SparseTensor<U, I>
where
    T: Clone,
    U: Clone,
    I: TensorStructure + StructureContract + Clone,
{
    type OutStorage<LCMData> = DenseTensor<LCMData, I>;
    fn least_common_storage<LCMData>(&self, other: &DenseTensor<T, I>) -> Self::OutStorage<LCMData>
    where
        for<'a, 'b> &'a Self::Data: std::ops::Mul<&'b T, Output = LCMData>,
        LCMData: Zero + Clone,
    {
        let mut final_structure = self.structure().clone();
        final_structure.merge(other.structure());
        DenseTensor::zero(final_structure)
    }

    fn empty<LCMData>(structure: Self::Structure) -> Self::OutStorage<LCMData>
    where
        for<'a, 'b> &'a Self::Data:
            std::ops::Mul<&'b <DenseTensor<T, I> as HasTensorData>::Data, Output = LCMData>,
        LCMData: Zero + Clone,
    {
        DenseTensor::zero(structure)
    }
}

impl<T, U, I> LeastCommonStorage<SparseTensor<T, I>> for DenseTensor<U, I>
where
    T: Clone,
    U: Clone,
    I: TensorStructure + StructureContract + Clone,
{
    type OutStorage<LCMData> = DenseTensor<LCMData, I>;
    fn least_common_storage<LCMData>(&self, other: &SparseTensor<T, I>) -> Self::OutStorage<LCMData>
    where
        for<'a, 'b> &'a Self::Data: std::ops::Mul<&'b T, Output = LCMData>,
        LCMData: Zero + Clone,
    {
        let mut final_structure = self.structure().clone();
        final_structure.merge(other.structure());
        DenseTensor::zero(final_structure)
    }

    fn empty<LCMData>(structure: Self::Structure) -> Self::OutStorage<LCMData>
    where
        for<'a, 'b> &'a Self::Data:
            std::ops::Mul<&'b <SparseTensor<T, I> as HasTensorData>::Data, Output = LCMData>,
        LCMData: Zero + Clone,
    {
        DenseTensor::zero(structure)
    }
}

impl<T, U, I> LeastCommonStorage<SparseTensor<T, I>> for SparseTensor<U, I>
where
    T: Clone,
    U: Clone,
    I: TensorStructure + StructureContract + Clone,
{
    type OutStorage<LCMData> = SparseTensor<LCMData, I>;
    fn least_common_storage<LCMData>(&self, other: &SparseTensor<T, I>) -> Self::OutStorage<LCMData>
    where
        for<'a, 'b> &'a Self::Data: std::ops::Mul<&'b T, Output = LCMData>,
        LCMData: Zero + Clone,
    {
        let mut final_structure = self.structure().clone();
        final_structure.merge(other.structure());
        SparseTensor::empty(final_structure)
    }

    fn empty<LCMData>(structure: Self::Structure) -> Self::OutStorage<LCMData>
    where
        for<'a, 'b> &'a Self::Data:
            std::ops::Mul<&'b <SparseTensor<T, I> as HasTensorData>::Data, Output = LCMData>,
        LCMData: Zero + Clone,
    {
        SparseTensor::empty(structure)
    }
}

trait ExteriorProduct<T> {
    type Out;
    fn exterior_product(&self, other: &T) -> Self::Out;
}

impl<T, U, LCMData> ExteriorProduct<T> for U
where
    U: LeastCommonStorage<T> + DataIterator<U::Data>,
    for<'a, 'b> &'a U::Data: std::ops::Mul<&'b T::Data, Output = LCMData>,
    T: LeastCommonStorage<U, OutStorage<LCMData> = U::OutStorage<LCMData>> + DataIterator<T::Data>,
    LCMData: Zero + Clone,
{
    type Out = U::OutStorage<LCMData>;

    fn exterior_product(&self, other: &T) -> Self::Out {
        let mut out = self.least_common_storage::<LCMData>(other);

        let stride = other.size();

        for (i, u) in self.flat_iter() {
            for (j, t) in other.flat_iter() {
                let _ = out.set_flat(i * stride + j, u * t);
            }
        }

        out
    }
}

impl<T, I> DenseTensor<T, I>
where
    T: for<'a> std::ops::AddAssign<&'a T>
        + for<'b> std::ops::SubAssign<&'b T>
        + Neg<Output = T>
        + Clone
        + std::fmt::Debug,
    I: TensorStructure + Clone + StructureContract,
{
    #[must_use]

    /// Contract the tensor with itself, i.e. trace over all matching indices.
    pub fn internal_contract(&self) -> Self {
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

impl<T, I> SparseTensor<T, I>
where
    T: for<'a> std::ops::AddAssign<&'a T>
        + for<'b> std::ops::SubAssign<&'b T>
        + std::ops::Neg<Output = T>
        + Clone
        + Default
        + PartialEq,
    I: TensorStructure + Clone + StructureContract,
{
    #[must_use]
    /// Contract the tensor with itself, i.e. trace over all matching indices.
    pub fn internal_contract(&self) -> Self {
        let trace = self.traces()[0];

        // println!("trace {:?}", trace);
        let mut new_structure = self.structure.clone();
        new_structure.trace(trace[0], trace[1]);

        let mut new_result = SparseTensor::empty(new_structure);
        for (idx, t) in self.iter_trace(trace).filter(|(_, t)| *t != T::default()) {
            new_result.set(&idx, t).unwrap_or_else(|_| unreachable!());
        }

        if new_result.traces().is_empty() {
            new_result
        } else {
            new_result.internal_contract()
        }
    }
}
pub trait Contract<T> {
    type LCM;
    fn contract(&self, other: &T) -> Option<Self::LCM>;
}

pub trait SingleContract<T> {
    type LCM;
    fn single_contract(&self, other: &T, i: usize, j: usize) -> Option<Self::LCM>;
}

trait MultiContract<T> {
    type LCM;
    fn multi_contract(&self, other: &T) -> Option<Self::LCM>;
}
pub trait DotProduct<T> {
    type Out: std::ops::AddAssign<Self::Out>
        + std::ops::SubAssign<Self::Out>
        + Neg<Output = Self::Out>
        + Clone
        + num::traits::Zero
        + for<'a> std::ops::AddAssign<&'a Self::Out>
        + for<'b> std::ops::SubAssign<&'b Self::Out>;
}

impl<T, U, O> DotProduct<T> for U
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = O>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = O>,
    O: std::ops::AddAssign<O>
        + std::ops::SubAssign<O>
        + Neg<Output = O>
        + Clone
        + num::traits::Zero
        + for<'a> std::ops::AddAssign<&'a O>
        + for<'b> std::ops::SubAssign<&'b O>,
{
    type Out = O;
}

impl<T, U, I> SingleContract<DenseTensor<T, I>> for DenseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = U::Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = U::Out>,
    U: DotProduct<T>,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;

    fn single_contract(&self, other: &DenseTensor<T, I>, i: usize, j: usize) -> Option<Self::LCM> {
        let final_structure = self.structure.merge_at(&other.structure, (i, j));
        let mut result_data = vec![U::Out::zero(); final_structure.size()];
        let mut result_index = 0;

        let mut self_iter = self.iter_fiber(i);
        let mut other_iter = other.iter_fiber(j);

        let fiber_representation: Representation = self.reps()[i];

        for fiber_a in self_iter.by_ref() {
            for fiber_b in other_iter.by_ref() {
                for k in 0..usize::from(fiber_representation) {
                    if fiber_representation.is_neg(k) {
                        result_data[result_index] -= fiber_a[k] * fiber_b[k];
                    } else {
                        result_data[result_index] += fiber_a[k] * fiber_b[k];
                    }
                }
                result_index += 1;
            }
            let _ = other_iter.reset();
        }
        let result = DenseTensor {
            data: result_data,
            structure: final_structure,
        };

        Some(result)
    }
}

impl<T, U, I> MultiContract<DenseTensor<T, I>> for DenseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = U::Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = U::Out>,
    U: DotProduct<T>,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;
    fn multi_contract(&self, other: &DenseTensor<T, I>) -> Option<Self::LCM> {
        let (permutation, self_matches, other_matches) =
            self.structure().match_indices(other.structure()).unwrap();

        let mut final_structure = self.structure.clone();
        final_structure.merge(&other.structure);

        // Initialize result tensor with default values
        let mut result_data = vec![U::Out::zero(); final_structure.size()];
        let mut result_index = 0;

        let mut selfiter = self.iter_multi_fibers_metric(&self_matches, permutation);

        let mut other_iter = other.iter_multi_fibers(&other_matches);
        while let Some(fiber_a) = selfiter.next() {
            for fiber_b in other_iter.by_ref() {
                for k in 0..fiber_a.len() {
                    if fiber_a[k].1 {
                        result_data[result_index] -= fiber_b[selfiter.map[k]] * fiber_a[k].0;
                    } else {
                        result_data[result_index] += fiber_b[selfiter.map[k]] * fiber_a[k].0;
                    }
                }
                result_index += 1;
            }
            let _ = other_iter.reset();
        }
        let result: DenseTensor<U::Out, I> = DenseTensor {
            data: result_data,
            structure: final_structure,
        };

        Some(result)
    }
}

impl<T, U> Contract<T> for U
where
    U: SingleContract<T>
        + MultiContract<T, LCM = <U as SingleContract<T>>::LCM>
        + ExteriorProduct<T, Out = <U as SingleContract<T>>::LCM>
        + TensorStructure,
    U::Structure: TensorStructure,
    T: SingleContract<U, LCM = <U as SingleContract<T>>::LCM>
        + MultiContract<U, LCM = <U as SingleContract<T>>::LCM>
        + ExteriorProduct<U, Out = <U as SingleContract<T>>::LCM>
        + TensorStructure<Structure = U::Structure>,
{
    type LCM = <U as SingleContract<T>>::LCM;
    fn contract(&self, other: &T) -> Option<Self::LCM> {
        if let Some((single, i, j)) = self.structure().match_index(other.structure()) {
            if i >= j {
                if single {
                    // println!("single");
                    return self.single_contract(other, i, j);
                }
                // println!("multi");
                return self.multi_contract(other);
            }
            // println!("flip");
            return other.contract(self);
        }
        // println!("exterior");
        let result = self.exterior_product(other);
        Some(result)
    }
}

impl<T, U, I> SingleContract<DenseTensor<T, I>> for SparseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = U::Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = U::Out>,
    U: DotProduct<T>,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;

    fn single_contract(&self, other: &DenseTensor<T, I>, i: usize, j: usize) -> Option<Self::LCM> {
        let final_structure = self.structure.merge_at(&other.structure, (i, j));
        let mut result_data = vec![U::Out::zero(); final_structure.size()];
        let mut result_index = 0;

        let mut self_iter = self.iter_fiber(i);
        let mut other_iter = other.iter_fiber(j);

        let fiber_representation: Representation = self.reps()[i];

        let pos = self.order();
        let stride = *final_structure
            .strides()
            .get(pos.wrapping_sub(2))
            .unwrap_or(&1);

        for (skipped, nonzeros, fiber_a) in self_iter.by_ref() {
            result_index += skipped * stride;

            for fiber_b in other_iter.by_ref() {
                for (i, k) in nonzeros.iter().enumerate() {
                    if fiber_representation.is_neg(*k) {
                        result_data[result_index] -= fiber_a[i] * fiber_b[*k];
                    } else {
                        result_data[result_index] += fiber_a[i] * fiber_b[*k];
                    }
                }
                result_index += 1;
            }

            result_index += other_iter.reset();
        }

        let result = DenseTensor {
            data: result_data,
            structure: final_structure,
        };

        Some(result)
    }
}

impl<T, U, I> SingleContract<SparseTensor<T, I>> for DenseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = U::Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = U::Out>,
    U: DotProduct<T>,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;

    fn single_contract(&self, other: &SparseTensor<T, I>, i: usize, j: usize) -> Option<Self::LCM> {
        let final_structure = self.structure.merge_at(&other.structure, (i, j));
        let mut result_data = vec![U::Out::zero(); final_structure.size()];
        let mut result_index = 0;

        let mut self_iter = self.iter_fiber(i);
        let mut other_iter = other.iter_fiber(j);

        let fiber_representation: Representation = self.reps()[i];

        for fiber_a in self_iter.by_ref() {
            for (skipped, nonzeros, fiber_b) in other_iter.by_ref() {
                result_index += skipped;

                for (i, k) in nonzeros.iter().enumerate() {
                    if fiber_representation.is_neg(*k) {
                        result_data[result_index] -= fiber_a[*k] * fiber_b[i];
                    } else {
                        result_data[result_index] += fiber_a[*k] * fiber_b[i];
                    }
                }
                result_index += 1;
            }
            result_index += other_iter.reset();
        }

        let result = DenseTensor {
            data: result_data,
            structure: final_structure,
        };

        Some(result)
    }
}

impl<T, U, I> MultiContract<DenseTensor<T, I>> for SparseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = U::Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = U::Out>,
    U: DotProduct<T>,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;
    fn multi_contract(&self, other: &DenseTensor<T, I>) -> Option<Self::LCM> {
        let (permutation, self_matches, other_matches) =
            self.structure().match_indices(other.structure()).unwrap();

        let mut final_structure = self.structure.clone();
        let mergepoint = final_structure.merge(&other.structure);

        // println!("mergepoint {:?}", mergepoint);
        let stride = if let Some(pos) = mergepoint {
            *final_structure
                .strides()
                .get(pos.saturating_sub(1))
                .unwrap_or(&1)
        } else {
            1
        };

        let mut result_data = vec![U::Out::zero(); final_structure.size()];
        let mut result_index = 0;

        let mut selfiter = self.iter_multi_fibers_metric(&self_matches, permutation);
        let mut other_iter = other.iter_multi_fibers(&other_matches);

        let mut max_skip = 0;
        while let Some((skipped, nonzeros, fiber_a)) = selfiter.next() {
            result_index += skipped * stride;
            if skipped > max_skip {
                max_skip = skipped;
            }
            for fiber_b in other_iter.by_ref() {
                for (i, k) in nonzeros.iter().enumerate() {
                    if fiber_a[i].1 {
                        result_data[result_index] -= fiber_a[i].0 * fiber_b[selfiter.map[*k]];
                    } else {
                        result_data[result_index] += fiber_a[i].0 * fiber_b[selfiter.map[*k]];
                    }
                }
                result_index += 1;
            }
            let _ = other_iter.reset();
        }
        if max_skip > 0 {
            // println!("skippedDS {max_skip}");
        }
        let result = DenseTensor {
            data: result_data,
            structure: final_structure,
        };

        Some(result)
    }
}

impl<T, U, I> MultiContract<SparseTensor<T, I>> for DenseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = U::Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = U::Out>,
    U: DotProduct<T>,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;

    fn multi_contract(&self, other: &SparseTensor<T, I>) -> Option<Self::LCM> {
        let (permutation, self_matches, other_matches) =
            self.structure().match_indices(other.structure()).unwrap();

        let mut final_structure = self.structure.clone();
        final_structure.merge(&other.structure);

        let mut result_data = vec![U::Out::zero(); final_structure.size()];
        let mut result_index = 0;

        let selfiter = self.iter_multi_fibers_metric(&self_matches, permutation.clone());
        let mut other_iter = other.iter_multi_fibers_metric(
            &other_matches,
            permutation.clone().inverse().normalize(true),
        );

        let mut max_skip = 0;
        for fiber_a in selfiter {
            while let Some((skipped, nonzeros, fiber_b)) = other_iter.next() {
                result_index += skipped;
                if skipped > max_skip {
                    max_skip = skipped;
                }

                for (i, k) in nonzeros.iter().enumerate() {
                    if fiber_a[other_iter.map[*k]].1 {
                        result_data[result_index] -= fiber_a[other_iter.map[*k]].0 * fiber_b[i].0;
                    } else {
                        result_data[result_index] += fiber_a[other_iter.map[*k]].0 * fiber_b[i].0;
                    }
                }
                result_index += 1;
            }
            result_index += other_iter.reset();
        }

        if max_skip > 0 {
            // println!("skippedSD {max_skip}");
        }

        let result = DenseTensor {
            data: result_data,
            structure: final_structure,
        };

        Some(result)
    }
}

impl<T, U, I> SingleContract<SparseTensor<T, I>> for SparseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = U::Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = U::Out>,
    U: DotProduct<T>,
    I: TensorStructure + Clone + StructureContract,
    U::Out: PartialEq,
{
    type LCM = SparseTensor<U::Out, I>;

    fn single_contract(&self, other: &SparseTensor<T, I>, i: usize, j: usize) -> Option<Self::LCM> {
        let final_structure = self.structure.merge_at(&other.structure, (i, j));
        let mut result_data = AHashMap::default();
        let mut result_index = 0;
        let pos = self.order();

        let mut self_iter = self.iter_fiber(i);
        let mut other_iter = other.iter_fiber(j);

        let metric = self.external_structure()[i].representation.negative();

        let stride = *final_structure
            .strides()
            .get(pos.wrapping_sub(2))
            .unwrap_or(&1);

        for (skipped_a, nonzeros_a, fiber_a) in self_iter.by_ref() {
            result_index += skipped_a * stride;
            for (skipped_b, nonzeros_b, fiber_b) in other_iter.by_ref() {
                result_index += skipped_b;
                let mut value = U::Out::zero();
                let mut nonzero = false;

                for (i, j, x) in nonzeros_a
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &x)| nonzeros_b.binary_search(&x).ok().map(|j| (i, j, x)))
                {
                    if metric[x] {
                        value -= fiber_a[i] * fiber_b[j];
                    } else {
                        value += fiber_a[i] * fiber_b[j];
                    }
                    nonzero = true;
                }

                if nonzero && value != U::Out::zero() {
                    result_data.insert(result_index, value);
                }
                result_index += 1;
            }
            result_index += other_iter.reset();
        }

        let result = SparseTensor {
            elements: result_data,
            structure: final_structure,
        };

        Some(result)
    }
}

impl<T, U, I> MultiContract<SparseTensor<T, I>> for SparseTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = U::Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = U::Out>,
    U: DotProduct<T>,
    I: TensorStructure + Clone + StructureContract,
    U::Out: PartialEq,
{
    type LCM = SparseTensor<U::Out, I>;

    fn multi_contract(&self, other: &SparseTensor<T, I>) -> Option<Self::LCM> {
        let (permutation, self_matches, other_matches) =
            self.structure().match_indices(other.structure()).unwrap();

        let mut final_structure = self.structure.clone();
        let mergepoint = final_structure.merge(&other.structure);
        let mut result_data = AHashMap::default();

        // println!("mergepoint {:?}", mergepoint);
        let stride = if let Some(pos) = mergepoint {
            *final_structure
                .strides()
                .get(pos.saturating_sub(1))
                .unwrap_or(&1)
        } else {
            1
        };

        let mut result_index = 0;

        let mut self_iter = self.iter_multi_fibers_metric(&self_matches, permutation);
        let mut other_iter = other.iter_multi_fibers(&other_matches);
        while let Some((skipped_a, nonzeros_a, fiber_a)) = self_iter.next() {
            result_index += skipped_a * stride;

            for (skipped_b, nonzeros_b, fiber_b) in other_iter.by_ref() {
                result_index += skipped_b;
                let mut value = U::Out::zero();
                let mut nonzero = false;

                for (i, j, _x) in nonzeros_a.iter().enumerate().filter_map(|(i, &x)| {
                    nonzeros_b
                        .binary_search(&self_iter.map[x])
                        .ok()
                        .map(|j| (i, j, x))
                }) {
                    if fiber_a[i].1 {
                        value -= fiber_a[i].0 * fiber_b[j];
                    } else {
                        value += fiber_a[i].0 * fiber_b[j];
                    }
                    nonzero = true;
                }

                if nonzero && value != U::Out::zero() {
                    result_data.insert(result_index, value);
                }
                result_index += 1;
            }
            result_index += other_iter.reset();
        }

        let result = SparseTensor {
            elements: result_data,
            structure: final_structure,
        };

        Some(result)
    }
}

impl<T, U, I, Out> Contract<DataTensor<T, I>> for DataTensor<U, I>
where
    for<'a, 'b> &'a U: std::ops::Mul<&'b T, Output = Out>,
    for<'a, 'b> &'a T: std::ops::Mul<&'b U, Output = Out>,
    Out: std::ops::AddAssign<Out>
        + std::ops::SubAssign<Out>
        + Neg<Output = Out>
        + Clone
        + Default
        + num::traits::Zero
        + for<'a> std::ops::AddAssign<&'a Out>
        + for<'b> std::ops::SubAssign<&'b Out>
        + std::fmt::Debug
        + std::cmp::PartialEq,
    I: Clone + Debug + TensorStructure + StructureContract,
    T: Clone + Debug,
    U: Clone + Debug,
{
    type LCM = DataTensor<Out, I>;
    fn contract(&self, other: &DataTensor<T, I>) -> Option<DataTensor<Out, I>> {
        match (self, other) {
            (DataTensor::Dense(s), DataTensor::Dense(o)) => Some(DataTensor::Dense(s.contract(o)?)),
            (DataTensor::Dense(s), DataTensor::Sparse(o)) => {
                Some(DataTensor::Dense(s.contract(o)?))
            }
            (DataTensor::Sparse(s), DataTensor::Dense(o)) => {
                Some(DataTensor::Dense(s.contract(o)?))
            }
            (DataTensor::Sparse(s), DataTensor::Sparse(o)) => {
                Some(DataTensor::Sparse(s.contract(o)?))
            }
        }
    }
}

impl<I> Contract<NumTensor<I>> for NumTensor<I>
where
    I: Clone + Debug + TensorStructure + StructureContract,
{
    type LCM = NumTensor<I>;
    fn contract(&self, other: &NumTensor<I>) -> Option<Self::LCM> {
        match (self, other) {
            (NumTensor::Float(a), NumTensor::Float(b)) => Some(NumTensor::Float(a.contract(b)?)),
            (NumTensor::Float(a), NumTensor::Complex(b)) => {
                Some(NumTensor::Complex(a.contract(b)?))
            }
            (NumTensor::Complex(a), NumTensor::Float(b)) => {
                Some(NumTensor::Complex(a.contract(b)?))
            }
            (NumTensor::Complex(a), NumTensor::Complex(b)) => {
                Some(NumTensor::Complex(a.contract(b)?))
            }
        }
    }
}
