use bitvec::vec::BitVec;
use log::trace;
use std::{cmp::Ordering, collections::HashMap};
// use num::Zero;
use crate::{
    algebra::{
        algebraic_traits::{IsZero, RefZero},
        upgrading_arithmetic::{FallibleAddAssign, FallibleSubAssign, TrySmallestUpgrade},
    },
    iterators::{
        CoreFlatFiberIterator, Fiber, FiberData, IteratableTensor, IteratesAlongFibers,
        ResetableIterator,
    },
    structure::{
        concrete_index::{ExpandedIndex, FlatIndex},
        HasStructure, StructureContract, TensorStructure,
    },
    tensors::data::{DataIterator, DenseTensor, SparseTensor},
};

use std::iter::Iterator;

use super::{ContractableWith, ContractionError, MultiContract, MultiContractInterleaved};

impl<T, U, I> MultiContract<DenseTensor<T, I>> for DenseTensor<U, I>
where
    U: ContractableWith<
        T,
        Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
    >,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;
    fn multi_contract(
        &self,
        other: &DenseTensor<T, I>,
        final_structure: I,
    ) -> Result<Self::LCM, ContractionError> {
        trace!("multi contract dense dense");
        let zero = self.data[0].try_upgrade().unwrap().into_owned().ref_zero();
        let (permutation, self_matches, other_matches) =
            self.structure().match_indices(other.structure()).unwrap();

        // Initialize result tensor with default values
        let mut result_data = vec![zero.clone(); final_structure.size()?];
        let mut result_index = 0;

        let selfiter = self
            .fiber_class(self_matches.as_slice().into())
            .iter_perm_metric(permutation);
        let mut other_iter = other.fiber_class(other_matches.as_slice().into()).iter();

        for mut fiber_a in selfiter {
            for mut fiber_b in other_iter.by_ref() {
                for (a, (neg, _)) in fiber_a.by_ref() {
                    if let Some((b, _)) = fiber_b.next() {
                        if neg {
                            result_data[result_index]
                                .sub_assign_fallible(&a.mul_fallible(b).unwrap());
                        } else {
                            result_data[result_index]
                                .add_assign_fallible(&a.mul_fallible(b).unwrap());
                        }
                    }
                }
                result_index += 1;
                fiber_a.reset();
            }
            other_iter.reset();
        }
        let result: DenseTensor<U::Out, I> = DenseTensor {
            data: result_data,
            structure: final_structure,
        };

        Ok(result)
    }
}

impl<T, U, I> MultiContract<DenseTensor<T, I>> for SparseTensor<U, I>
where
    U: ContractableWith<
        T,
        Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
    >,
    T: TrySmallestUpgrade<U, LCM = U::Out>,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;
    fn multi_contract(
        &self,
        other: &DenseTensor<T, I>,
        final_structure: I,
    ) -> Result<Self::LCM, ContractionError> {
        trace!("multi contract sparse dense");
        let zero = if let Some((_, s)) = self.flat_iter().next() {
            s.try_upgrade().unwrap().as_ref().ref_zero()
        } else if let Some((_, o)) = other.iter_flat().next() {
            o.try_upgrade().unwrap().as_ref().ref_zero()
        } else {
            return Err(ContractionError::EmptySparse);
        };
        // let zero = other.data[0].try_upgrade().unwrap().as_ref().ref_zero();
        let (permutation, self_matches, other_matches) =
            self.structure().match_indices(other.structure()).unwrap();

        // let mut final_structure = self.structure.clone();
        // let _ = final_structure.merge(&other.structure);

        let mut result_data = vec![zero.clone(); final_structure.size()?];
        let mut result_index = 0;

        let selfiter = self
            .fiber_class(self_matches.as_slice().into())
            .iter_perm_metric(permutation);
        let mut other_iter = other.fiber_class(other_matches.as_slice().into()).iter();

        for mut fiber_a in selfiter {
            for mut fiber_b in other_iter.by_ref() {
                for (a, skip, (neg, _)) in fiber_a.by_ref() {
                    if let Some((b, _)) = fiber_b.by_ref().nth(skip) {
                        if neg {
                            result_data[result_index]
                                .sub_assign_fallible(&a.mul_fallible(b).unwrap());
                        } else {
                            result_data[result_index]
                                .add_assign_fallible(&a.mul_fallible(b).unwrap());
                        }
                    }
                }
                result_index += 1;
                fiber_a.reset();
            }
            other_iter.reset();
        }

        let result = DenseTensor {
            data: result_data,
            structure: final_structure,
        };

        Ok(result)
    }
}

impl<T, U, I> MultiContract<SparseTensor<T, I>> for DenseTensor<U, I>
where
    U: ContractableWith<
        T,
        Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
    >,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;

    fn multi_contract(
        &self,
        other: &SparseTensor<T, I>,
        final_structure: I,
    ) -> Result<Self::LCM, ContractionError> {
        trace!("multi contract dense sparse");
        let zero = self.data[0].try_upgrade().unwrap().as_ref().ref_zero();
        let (permutation, self_matches, other_matches) =
            self.structure().match_indices(other.structure()).unwrap();

        let mut result_data = vec![zero.clone(); final_structure.size()?];
        let mut result_index = 0;

        let selfiter = self
            .fiber_class(self_matches.as_slice().into())
            .iter_perm_metric(permutation);
        let mut other_iter = other.fiber_class(other_matches.as_slice().into()).iter();

        for mut fiber_a in selfiter {
            for mut fiber_b in other_iter.by_ref() {
                for (b, skip, _) in fiber_b.by_ref() {
                    if let Some((a, (neg, _))) = fiber_a.by_ref().nth(skip) {
                        if neg {
                            result_data[result_index]
                                .sub_assign_fallible(&a.mul_fallible(b).unwrap());
                        } else {
                            result_data[result_index]
                                .add_assign_fallible(&a.mul_fallible(b).unwrap());
                        }
                    }
                }
                result_index += 1;
                fiber_a.reset();
            }
            other_iter.reset();
        }

        let result = DenseTensor {
            data: result_data,
            structure: final_structure,
        };

        Ok(result)
    }
}

impl<T, U, I> MultiContract<SparseTensor<T, I>> for SparseTensor<U, I>
where
    U: ContractableWith<
        T,
        Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
    >,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = SparseTensor<U::Out, I>;
    #[allow(clippy::comparison_chain)]
    fn multi_contract(
        &self,
        other: &SparseTensor<T, I>,
        final_structure: I,
    ) -> Result<Self::LCM, ContractionError> {
        trace!("multi contract sparse sparse");
        let (permutation, self_matches, other_matches) =
            self.structure().match_indices(other.structure()).unwrap();

        // let mut final_structure = self.structure.clone();
        // let _ = final_structure.merge(&other.structure);
        let mut result_data = HashMap::default();

        let zero = self.zero.try_upgrade().unwrap().as_ref().ref_zero();
        if self.flat_iter().next().is_some() {
            let mut result_index = 0;

            let self_iter = self
                .fiber_class(self_matches.as_slice().into())
                .iter_perm_metric(permutation);
            let mut other_iter = other.fiber_class(other_matches.as_slice().into()).iter();

            for mut fiber_a in self_iter {
                for mut fiber_b in other_iter.by_ref() {
                    let mut items = fiber_a
                        .next()
                        .map(|(a, skip, (neg, _))| (a, skip, neg))
                        .zip(fiber_b.next().map(|(b, skip, _)| (b, skip)));

                    let mut value = zero.clone();
                    let mut nonzero = false;

                    while let Some(((a, skip_a, neg), (b, skip_b))) = items {
                        if skip_a > skip_b {
                            let b = fiber_b
                                .by_ref()
                                .next()
                                .map(|(b, skip, _)| (b, skip + skip_b + 1));
                            items = Some((a, skip_a, neg)).zip(b);
                        } else if skip_b > skip_a {
                            let a = fiber_a
                                .by_ref()
                                .next()
                                .map(|(a, skip, (neg, _))| (a, skip + skip_a + 1, neg));
                            items = a.zip(Some((b, skip_b)));
                        } else {
                            // println!("v{:?}", value);
                            if neg {
                                value.sub_assign_fallible(&a.mul_fallible(b).unwrap());
                            } else {
                                value.add_assign_fallible(&a.mul_fallible(b).unwrap());
                            }
                            let b = fiber_b
                                .by_ref()
                                .next()
                                .map(|(b, skip, _)| (b, skip + skip_b + 1));
                            let a = fiber_a
                                .by_ref()
                                .next()
                                .map(|(a, skip, (neg, _))| (a, skip + skip_a + 1, neg));
                            items = a.zip(b);
                            nonzero = true;
                        }
                    }
                    if nonzero && value.is_non_zero() {
                        result_data.insert(result_index.into(), value);
                    }
                    result_index += 1;
                    fiber_a.reset();
                }
                other_iter.reset();
            }
        }
        let result = SparseTensor {
            zero,
            elements: result_data,
            structure: final_structure,
        };

        Ok(result)
    }
}

impl<T, U, I> MultiContractInterleaved<DenseTensor<T, I>> for DenseTensor<U, I>
where
    U: ContractableWith<
        T,
        Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
    >,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;
    fn multi_contract_interleaved(
        &self,
        other: &DenseTensor<T, I>,
        pos_self: BitVec,
        pos_other: BitVec,
        resulting_structure: <Self::LCM as HasStructure>::Structure,
        resulting_partition: bitvec::prelude::BitVec,
    ) -> Result<Self::LCM, ContractionError> {
        // Assume we're contracting the first positions for now - this needs to be updated
        let zero = self.data[0].try_upgrade().unwrap().into_owned().ref_zero();
        // println!("Interleaving dense dense multi"); // Initialize result tensor with default values
        let mut result_data = vec![zero.clone(); resulting_structure.size()?];

        let self_fiber_class = Fiber::from(&resulting_partition, &resulting_structure); //We use the partition as a filter here, for indices that belong to self, vs those that belong to other
        let (mut self_fiber_class_iter, mut other_fiber_class_iter) =
            CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class); // these are iterators over the open indices of self and other, except expressed in the flat indices of the resulting structure

        let mut iter_self = self.fiber(FiberData::from(&pos_self)).iter_metric(); //The summed over index comes from the actual self structure (and is a single index)
        let mut iter_other = other.fiber(FiberData::from(&pos_other)).iter(); // same for other

        let mut exp_self = self.expanded_index(FlatIndex::from(0)).unwrap();
        let mut exp_other = other.expanded_index(FlatIndex::from(0)).unwrap();

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
                        .partition(|(i, _)| resulting_partition[*i]);

                for (i, v) in pos_self.iter_zeros().zip(expa) {
                    exp_self.indices[i] = v;
                }

                for (i, v) in pos_other.iter_zeros().zip(expb) {
                    exp_other.indices[i] = v;
                }
                // And now we flatten
                let shift_a = self.structure().flat_index(&exp_self).unwrap();
                let shift_b = other.structure().flat_index(&exp_other).unwrap();

                // println!("shift_a: {:?}, shift_b: {:?}", shift_a, shift_b);

                // we shift the fiber start by the flattened indices. This also resets the iterator
                iter_self.shift(shift_a.into());
                iter_other.shift(shift_b.into());

                for ((a, (neg, _)), (b, _)) in (iter_self.by_ref()).zip(iter_other.by_ref()) {
                    if neg {
                        result_data[usize::from(result_index)]
                            .sub_assign_fallible(&(a.mul_fallible(b).unwrap()));
                    } else {
                        result_data[usize::from(result_index)]
                            .add_assign_fallible(&a.mul_fallible(b).unwrap());
                    }
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

impl<T, U, I> MultiContractInterleaved<DenseTensor<T, I>> for SparseTensor<U, I>
where
    U: ContractableWith<
        T,
        Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
    >,
    T: TrySmallestUpgrade<U, LCM = U::Out>,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;
    fn multi_contract_interleaved(
        &self,
        other: &DenseTensor<T, I>,
        pos_self: BitVec,
        pos_other: BitVec,
        resulting_structure: <Self::LCM as HasStructure>::Structure,
        resulting_partition: bitvec::prelude::BitVec,
    ) -> Result<Self::LCM, ContractionError> {
        let zero = if let Some((_, s)) = self.flat_iter().next() {
            s.try_upgrade().unwrap().as_ref().ref_zero()
        } else if let Some((_, o)) = other.iter_flat().next() {
            o.try_upgrade().unwrap().as_ref().ref_zero()
        } else {
            return Err(ContractionError::EmptySparse);
        };

        let mut result_data = vec![zero.clone(); resulting_structure.size()?];

        let self_fiber_class = Fiber::from(&resulting_partition, &resulting_structure); //We use the partition as a filter here, for indices that belong to self, vs those that belong to other
        let (mut self_fiber_class_iter, mut other_fiber_class_iter) =
            CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class); // these are iterators over the open indices of self and other, except expressed in the flat indices of the resulting structure
        let mut iter_self = self.fiber(FiberData::from(&pos_self)).iter_metric(); //The summed over index comes from the actual self structure (and is a single index)
        let mut iter_other = other.fiber(FiberData::from(&pos_other)).iter(); // same for other

        let mut exp_self = self.expanded_index(FlatIndex::from(0)).unwrap();
        let mut exp_other = other.expanded_index(FlatIndex::from(0)).unwrap();

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
                        .partition(|(i, _)| resulting_partition[*i]);

                // println!("expa: {:?}", expa);
                for (i, v) in pos_self.iter_zeros().zip(expa) {
                    exp_self.indices[i] = v;
                }

                for (i, v) in pos_other.iter_zeros().zip(expb) {
                    exp_other.indices[i] = v;
                }

                // println!("exp_self: {:?}", exp_self);
                // println!("exp_other: {:?}", exp_other);
                // println!("")
                // And now we flatten
                let shift_a = self.structure().flat_index(&exp_self).unwrap();
                let shift_b = other.structure().flat_index(&exp_other).unwrap();

                // we shift the fiber start by the flattened indices. This also resets the iterator
                iter_self.shift(shift_a.into());
                iter_other.shift(shift_b.into());

                for (a, skip, (neg, _)) in iter_self.by_ref() {
                    if let Some((b, _)) = iter_other.by_ref().nth(skip) {
                        if neg {
                            result_data[usize::from(result_index)]
                                .sub_assign_fallible(&a.mul_fallible(b).unwrap());
                        } else {
                            result_data[usize::from(result_index)]
                                .add_assign_fallible(&a.mul_fallible(b).unwrap());
                        }
                    }
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

impl<T, U, I> MultiContractInterleaved<SparseTensor<T, I>> for DenseTensor<U, I>
where
    U: ContractableWith<
        T,
        Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
    >,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = DenseTensor<U::Out, I>;

    fn multi_contract_interleaved(
        &self,
        other: &SparseTensor<T, I>,
        pos_self: BitVec,
        pos_other: BitVec,
        resulting_structure: <Self::LCM as HasStructure>::Structure,
        resulting_partition: bitvec::prelude::BitVec,
    ) -> Result<Self::LCM, ContractionError> {
        let zero = self.data[0].try_upgrade().unwrap().into_owned().ref_zero();
        let mut result_data = vec![zero.clone(); resulting_structure.size()?];

        let self_fiber_class = Fiber::from(&resulting_partition, &resulting_structure); //We use the partition as a filter here, for indices that belong to self, vs those that belong to other
        let (mut self_fiber_class_iter, mut other_fiber_class_iter) =
            CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class); // these are iterators over the open indices of self and other, except expressed in the flat indices of the resulting structure

        let mut iter_self = self.fiber(FiberData::from(&pos_self)).iter_metric(); //The summed over index comes from the actual self structure (and is a single index)
        let mut iter_other = other.fiber(FiberData::from(&pos_other)).iter(); // same for other

        let mut exp_self = self.expanded_index(FlatIndex::from(0)).unwrap();
        let mut exp_other = other.expanded_index(FlatIndex::from(0)).unwrap();

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
                        .partition(|(i, _)| resulting_partition[*i]);

                // println!("expa: {:?}", expa);
                for (i, v) in pos_self.iter_zeros().zip(expa) {
                    exp_self.indices[i] = v;
                }

                for (i, v) in pos_other.iter_zeros().zip(expb) {
                    exp_other.indices[i] = v;
                }

                // println!("exp_self: {:?}", exp_self);
                // println!("exp_other: {:?}", exp_other);
                // println!("")
                // And now we flatten
                let shift_a = self.structure().flat_index(&exp_self).unwrap();
                let shift_b = other.structure().flat_index(&exp_other).unwrap();

                // we shift the fiber start by the flattened indices. This also resets the iterator
                iter_self.shift(shift_a.into());
                iter_other.shift(shift_b.into());

                for (b, skip, _) in iter_other.by_ref() {
                    if let Some((a, (neg, _))) = iter_self.by_ref().nth(skip) {
                        if neg {
                            result_data[usize::from(result_index)]
                                .sub_assign_fallible(&a.mul_fallible(b).unwrap());
                        } else {
                            result_data[usize::from(result_index)]
                                .add_assign_fallible(&a.mul_fallible(b).unwrap());
                        }
                    }
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

impl<T, U, I> MultiContractInterleaved<SparseTensor<T, I>> for SparseTensor<U, I>
where
    U: ContractableWith<
        T,
        Out: FallibleAddAssign<U::Out> + FallibleSubAssign<U::Out> + Clone + RefZero + IsZero,
    >,
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = SparseTensor<U::Out, I>;

    fn multi_contract_interleaved(
        &self,
        other: &SparseTensor<T, I>,
        pos_self: BitVec,
        pos_other: BitVec,
        resulting_structure: <Self::LCM as HasStructure>::Structure,
        resulting_partition: bitvec::prelude::BitVec,
    ) -> Result<Self::LCM, ContractionError> {
        let mut result_data = HashMap::default();
        let zero = self.zero.try_upgrade().unwrap().as_ref().ref_zero();
        if self.flat_iter().next().is_some() {
            let self_fiber_class = Fiber::from(&resulting_partition, &resulting_structure); //We use the partition as a filter here, for indices that belong to self, vs those that belong to other
            let (mut self_fiber_class_iter, mut other_fiber_class_iter) =
                CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class); // these are iterators over the open indices of self and other, except expressed in the flat indices of the resulting structure

            let mut iter_self = self.fiber(FiberData::from(&pos_self)).iter_metric(); //The summed over index comes from the actual self structure (and is a single index)
            let mut iter_other = other.fiber(FiberData::from(&pos_other)).iter(); // same for other

            let mut exp_self = self.expanded_index(FlatIndex::from(0)).unwrap();
            let mut exp_other = other.expanded_index(FlatIndex::from(0)).unwrap();

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
                            .partition(|(i, _)| resulting_partition[*i]);

                    // println!("expa: {:?}", expa);
                    for (i, v) in pos_self.iter_zeros().zip(expa) {
                        exp_self.indices[i] = v;
                    }

                    for (i, v) in pos_other.iter_zeros().zip(expb) {
                        exp_other.indices[i] = v;
                    }

                    // println!("exp_self: {:?}", exp_self);
                    // println!("exp_other: {:?}", exp_other);
                    // println!("")
                    // And now we flatten
                    let shift_a = self.structure().flat_index(&exp_self).unwrap();
                    let shift_b = other.structure().flat_index(&exp_other).unwrap();

                    // we shift the fiber start by the flattened indices. This also resets the iterator
                    iter_self.shift(shift_a.into());
                    iter_other.shift(shift_b.into());

                    let mut items = iter_self
                        .next()
                        .map(|(a, skip, (neg, _))| (a, skip, neg))
                        .zip(iter_other.next().map(|(b, skip, _)| (b, skip)));

                    let mut value = zero.clone();
                    let mut nonzero = false;

                    while let Some(((a, skip_a, neg), (b, skip_b))) = items {
                        match skip_a.cmp(&skip_b) {
                            std::cmp::Ordering::Greater => {
                                let b = iter_other
                                    .by_ref()
                                    .next()
                                    .map(|(b, skip, _)| (b, skip + skip_b + 1));
                                items = Some((a, skip_a, neg)).zip(b);
                            }
                            Ordering::Less => {
                                //Advance iter_self
                                let a = iter_self
                                    .by_ref()
                                    .next()
                                    .map(|(a, skip, (neg, _))| (a, skip + skip_a + 1, neg));
                                items = a.zip(Some((b, skip_b)));
                            }
                            _ => {
                                // skip_a == skip_b
                                if neg {
                                    value.sub_assign_fallible(&a.mul_fallible(b).unwrap());
                                } else {
                                    value.add_assign_fallible(&a.mul_fallible(b).unwrap());
                                }
                                let b = iter_other
                                    .by_ref()
                                    .next()
                                    .map(|(b, skip, _)| (b, skip + skip_b + 1));
                                let a = iter_self
                                    .by_ref()
                                    .next()
                                    .map(|(a, skip, (neg, _))| (a, skip + skip_a + 1, neg));
                                items = a.zip(b);
                                nonzero = true;
                            }
                        }
                    }

                    if nonzero && !value.is_zero() {
                        result_data.insert(result_index, value);
                    }
                }
                other_fiber_class_iter.reset();
            }
        }
        let result = SparseTensor {
            zero,
            elements: result_data,
            structure: resulting_structure,
        };

        Ok(result)
    }
}
