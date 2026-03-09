use std::collections::HashMap;

use crate::{
    algebra::{
        algebraic_traits::{IsZero, RefZero},
        upgrading_arithmetic::{FallibleAddAssign, FallibleSubAssign, TrySmallestUpgrade},
    },
    iterators::{
        CoreFlatFiberIterator, Fiber, IteratableTensor, IteratesAlongFibers, ResetableIterator,
    },
    structure::{HasStructure, StructureContract, TensorStructure, concrete_index::ExpandedIndex},
    tensors::data::{DataIterator, DenseTensor, GetTensorData, SetTensorData, SparseTensor},
};

use std::iter::Iterator;

use super::{ContractableWith, ContractionError, ExteriorProduct, ExteriorProductInterleaved};

impl<T, U, I> ExteriorProduct<SparseTensor<T, I>> for SparseTensor<U, I>
where
    U: ContractableWith<
            T,
            Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
        >,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = SparseTensor<U::Out, I>;

    fn exterior_product(
        &self,
        other: &SparseTensor<T, I>,
        final_structure: I,
    ) -> Result<Self::LCM, ContractionError> {
        let mut out = SparseTensor::empty(
            final_structure,
            self.zero.mul_fallible(&other.zero).unwrap(),
        );
        let stride = other.size()?;

        for (i, u) in self.flat_iter() {
            for (j, t) in other.flat_iter() {
                let _ = out.set_flat(i * stride + j, u.mul_fallible(t).unwrap());
            }
        }

        Ok(out)
    }
}

impl<T, U, I> ExteriorProduct<DenseTensor<T, I>> for DenseTensor<U, I>
where
    U: ContractableWith<
            T,
            Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
        >,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;

    fn exterior_product(
        &self,
        other: &DenseTensor<T, I>,
        final_structure: I,
    ) -> Result<Self::LCM, ContractionError> {
        let zero = self.data[0].try_upgrade().unwrap().into_owned().ref_zero();
        let mut out = DenseTensor {
            data: vec![zero.clone(); final_structure.size()?],
            structure: final_structure,
        };

        let stride = other.size()?;

        for (i, u) in self.flat_iter() {
            for (j, t) in other.flat_iter() {
                let _ = out.set_flat(i * stride + j, u.mul_fallible(t).unwrap());
            }
        }

        Ok(out)
    }
}

impl<T, U, I> ExteriorProduct<SparseTensor<T, I>> for DenseTensor<U, I>
where
    U: ContractableWith<
            T,
            Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
        >,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;

    fn exterior_product(
        &self,
        other: &SparseTensor<T, I>,
        final_structure: I,
    ) -> Result<Self::LCM, ContractionError> {
        let zero = self.data[0].try_upgrade().unwrap().into_owned().ref_zero();
        let mut out = DenseTensor {
            data: vec![zero.clone(); final_structure.size()?],
            structure: final_structure,
        };
        let stride = other.size()?;

        for (i, u) in self.flat_iter() {
            for (j, t) in other.flat_iter() {
                let _ = out.set_flat(i * stride + j, u.mul_fallible(t).unwrap());
            }
        }

        Ok(out)
    }
}

impl<T, U, I> ExteriorProduct<DenseTensor<T, I>> for SparseTensor<U, I>
where
    U: ContractableWith<
            T,
            Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
        >,
    T: TrySmallestUpgrade<U, LCM = U::Out>,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;

    fn exterior_product(
        &self,
        other: &DenseTensor<T, I>,
        final_structure: I,
    ) -> Result<Self::LCM, ContractionError> {
        let zero = if let Some((_, s)) = self.flat_iter().next() {
            s.try_upgrade().unwrap().as_ref().ref_zero()
        } else if let Some((_, o)) = other.iter_flat().next() {
            o.try_upgrade().unwrap().as_ref().ref_zero()
        } else {
            return Err(ContractionError::EmptySparse);
        };
        let mut out = DenseTensor {
            data: vec![zero.clone(); final_structure.size()?],
            structure: final_structure,
        };

        let stride = other.size()?;

        for (i, u) in self.flat_iter() {
            for (j, t) in other.flat_iter() {
                let _ = out.set_flat(i * stride + j, u.mul_fallible(t).unwrap());
            }
        }

        Ok(out)
    }
}

impl<T, U, I> ExteriorProductInterleaved<SparseTensor<T, I>> for SparseTensor<U, I>
where
    U: ContractableWith<
            T,
            Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
        >,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = SparseTensor<U::Out, I>;

    fn exterior_product_interleaved(
        &self,
        other: &SparseTensor<T, I>,
        resulting_structure: <Self::LCM as crate::structure::HasStructure>::Structure,
        partition: bitvec::prelude::BitVec,
    ) -> Result<Self::LCM, ContractionError> {
        let mut result_data = HashMap::default();
        if let Some((_, _)) = self.flat_iter().next() {
            let self_fiber_class = Fiber::from(&partition, &resulting_structure); //We use the partition as a filter here, for indices that belong to self, vs those that belong to other

            let (mut self_fiber_class_iter, mut other_fiber_class_iter) =
                CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class); // these are iterators over the open indices of self and other, except expressed in the flat indices of the resulting structure

            //We first iterate over the free indices (self_fiber_class)
            for fiber_class_a_id in self_fiber_class_iter.by_ref() {
                for fiber_class_b_id in other_fiber_class_iter.by_ref() {
                    // This is the index in the resulting structure for these two class indices
                    let result_index = fiber_class_a_id + fiber_class_b_id;

                    //To obtain the corresponding flat indices for the self and other we partition the expanded index
                    let ((_, expa), (_, expb)): ((Vec<_>, ExpandedIndex), (Vec<_>, ExpandedIndex)) =
                        resulting_structure
                            .expanded_index(result_index)?
                            .into_iter()
                            .enumerate()
                            .partition(|(i, _)| partition[*i]);

                    // And now we flatten
                    let i = self.structure().flat_index(expa).unwrap();
                    let j = other.structure().flat_index(expb).unwrap();

                    if let Some(u) = self.get_ref_linear(i)
                        && let Some(v) = other.get_ref_linear(j)
                    {
                        result_data.insert(result_index, u.mul_fallible(v).unwrap());
                    }
                }
                other_fiber_class_iter.reset();
            }
        }
        let result = SparseTensor {
            zero: self.zero.mul_fallible(&other.zero).unwrap(),
            elements: result_data,
            structure: resulting_structure,
        };

        Ok(result)
    }
}

impl<T, U, I> ExteriorProductInterleaved<DenseTensor<T, I>> for DenseTensor<U, I>
where
    U: ContractableWith<
            T,
            Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
        >,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;

    fn exterior_product_interleaved(
        &self,
        other: &DenseTensor<T, I>,
        resulting_structure: <Self::LCM as crate::structure::HasStructure>::Structure,
        partition: bitvec::prelude::BitVec,
    ) -> Result<Self::LCM, ContractionError> {
        let zero = self.data[0].try_upgrade().unwrap().into_owned().ref_zero();
        let mut result_data = vec![zero.clone(); resulting_structure.size()?];

        let self_fiber_class = Fiber::from(&partition, &resulting_structure); //We use the partition as a filter here, for indices that belong to self, vs those that belong to other
        let (mut self_fiber_class_iter, mut other_fiber_class_iter) =
            CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class); // these are iterators over the open indices of self and other, except expressed in the flat indices of the resulting structure

        //We first iterate over the free indices (self_fiber_class)
        for fiber_class_a_id in self_fiber_class_iter.by_ref() {
            for fiber_class_b_id in other_fiber_class_iter.by_ref() {
                // This is the index in the resulting structure for these two class indices
                let result_index = fiber_class_a_id + fiber_class_b_id;

                //To obtain the corresponding flat indices for the self and other we partition the expanded index
                let ((_, expa), (_, expb)): ((Vec<_>, ExpandedIndex), (Vec<_>, ExpandedIndex)) =
                    resulting_structure
                        .expanded_index(result_index)?
                        .into_iter()
                        .enumerate()
                        .partition(|(i, _)| partition[*i]);

                // And now we flatten
                let i = self.structure().flat_index(expa).unwrap();
                let j = other.structure().flat_index(expb).unwrap();

                result_data[usize::from(result_index)] = self
                    .get_ref_linear(i)
                    .unwrap()
                    .mul_fallible(other.get_ref_linear(j).unwrap())
                    .unwrap();
            }
            other_fiber_class_iter.reset();
        }
        let result = DenseTensor {
            data: result_data,
            structure: resulting_structure,
        };

        Ok(result)
    }
}

impl<T, U, I> ExteriorProductInterleaved<SparseTensor<T, I>> for DenseTensor<U, I>
where
    U: ContractableWith<
            T,
            Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
        >,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;

    fn exterior_product_interleaved(
        &self,
        other: &SparseTensor<T, I>,
        resulting_structure: <Self::LCM as crate::structure::HasStructure>::Structure,
        partition: bitvec::prelude::BitVec,
    ) -> Result<Self::LCM, ContractionError> {
        let zero = self.data[0].try_upgrade().unwrap().into_owned().ref_zero();
        let mut result_data = vec![zero.clone(); resulting_structure.size()?];

        let self_fiber_class = Fiber::from(&partition, &resulting_structure); //We use the partition as a filter here, for indices that belong to self, vs those that belong to other
        let (mut self_fiber_class_iter, mut other_fiber_class_iter) =
            CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class); // these are iterators over the open indices of self and other, except expressed in the flat indices of the resulting structure

        //We first iterate over the free indices (self_fiber_class)
        for fiber_class_a_id in self_fiber_class_iter.by_ref() {
            for fiber_class_b_id in other_fiber_class_iter.by_ref() {
                // This is the index in the resulting structure for these two class indices
                let result_index = fiber_class_a_id + fiber_class_b_id;

                //To obtain the corresponding flat indices for the self and other we partition the expanded index
                let ((_, expa), (_, expb)): ((Vec<_>, ExpandedIndex), (Vec<_>, ExpandedIndex)) =
                    resulting_structure
                        .expanded_index(result_index)?
                        .into_iter()
                        .enumerate()
                        .partition(|(i, _)| partition[*i]);

                // And now we flatten
                let i = self.structure().flat_index(expa).unwrap();
                let j = other.structure().flat_index(expb).unwrap();

                if let Some(t) = other.get_ref_linear(j) {
                    result_data[usize::from(result_index)] =
                        self.get_ref_linear(i).unwrap().mul_fallible(t).unwrap();
                }
            }
            other_fiber_class_iter.reset();
        }
        let result = DenseTensor {
            data: result_data,
            structure: resulting_structure,
        };

        Ok(result)
    }
}

impl<T, U, I> ExteriorProductInterleaved<DenseTensor<T, I>> for SparseTensor<U, I>
where
    U: ContractableWith<
            T,
            Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
        >,
    T: TrySmallestUpgrade<U, LCM = U::Out>,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;

    fn exterior_product_interleaved(
        &self,
        other: &DenseTensor<T, I>,
        resulting_structure: <Self::LCM as crate::structure::HasStructure>::Structure,
        partition: bitvec::prelude::BitVec,
    ) -> Result<Self::LCM, ContractionError> {
        let zero = if let Some((_, s)) = self.flat_iter().next() {
            s.try_upgrade().unwrap().as_ref().ref_zero()
        } else if let Some((_, o)) = other.iter_flat().next() {
            o.try_upgrade().unwrap().as_ref().ref_zero()
        } else {
            return Err(ContractionError::EmptySparse);
        };
        let mut result_data = vec![zero.clone(); resulting_structure.size()?];

        let self_fiber_class = Fiber::from(&partition, &resulting_structure); //We use the partition as a filter here, for indices that belong to self, vs those that belong to other
        let (mut self_fiber_class_iter, mut other_fiber_class_iter) =
            CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class); // these are iterators over the open indices of self and other, except expressed in the flat indices of the resulting structure

        //We first iterate over the free indices (self_fiber_class)
        for fiber_class_a_id in self_fiber_class_iter.by_ref() {
            for fiber_class_b_id in other_fiber_class_iter.by_ref() {
                // This is the index in the resulting structure for these two class indices
                let result_index = fiber_class_a_id + fiber_class_b_id;

                //To obtain the corresponding flat indices for the self and other we partition the expanded index
                let ((_, expa), (_, expb)): ((Vec<_>, ExpandedIndex), (Vec<_>, ExpandedIndex)) =
                    resulting_structure
                        .expanded_index(result_index)?
                        .into_iter()
                        .enumerate()
                        .partition(|(i, _)| partition[*i]);

                // And now we flatten
                let i = self.structure().flat_index(expa).unwrap();
                let j = other.structure().flat_index(expb).unwrap();

                if let Some(u) = self.get_ref_linear(i) {
                    result_data[usize::from(result_index)] =
                        u.mul_fallible(other.get_ref_linear(j).unwrap()).unwrap();
                }
            }
            other_fiber_class_iter.reset();
        }
        let result = DenseTensor {
            data: result_data,
            structure: resulting_structure,
        };

        Ok(result)
    }
}
