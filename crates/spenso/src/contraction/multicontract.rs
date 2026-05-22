use log::trace;
use std::collections::HashMap;
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
        HasStructure, SlotIndex, StructureContract, TensorStructure,
        concrete_index::{ExpandedIndex, FlatIndex},
        slot::{DualSlotTo, IsAbstractSlot},
    },
    tensors::data::{DataIterator, DenseTensor, SparseTensor},
};

use linnet::half_edge::subgraph::{SubSetLike, subset::SubSet};

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
        pos_self: SubSet<SlotIndex>,
        pos_other: SubSet<SlotIndex>,
        resulting_structure: <Self::LCM as HasStructure>::Structure,
        resulting_partition: SubSet<SlotIndex>,
    ) -> Result<Self::LCM, ContractionError> {
        // Assume we're contracting the first positions for now - this needs to be updated
        let zero = self.data[0].try_upgrade().unwrap().into_owned().ref_zero();

        let pos_self_inv = !&pos_self;
        let pos_other_inv = !&pos_other;
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
                        .partition(|(i, _)| resulting_partition[SlotIndex(*i)]);

                for (i, v) in pos_self_inv.included_iter().zip(expa) {
                    exp_self.indices[i.0] = v;
                }

                for (i, v) in pos_other_inv.included_iter().zip(expb) {
                    exp_other.indices[i.0] = v;
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
        pos_self: SubSet<SlotIndex>,
        pos_other: SubSet<SlotIndex>,
        resulting_structure: <Self::LCM as HasStructure>::Structure,
        resulting_partition: SubSet<SlotIndex>,
    ) -> Result<Self::LCM, ContractionError> {
        let zero = if let Some((_, s)) = self.flat_iter().next() {
            s.try_upgrade().unwrap().as_ref().ref_zero()
        } else if let Some((_, o)) = other.iter_flat().next() {
            o.try_upgrade().unwrap().as_ref().ref_zero()
        } else {
            return Err(ContractionError::EmptySparse);
        };
        let pos_self_inv = !&pos_self;
        let pos_other_inv = !&pos_other;

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
                        .partition(|(i, _)| resulting_partition[SlotIndex(*i)]);

                // println!("expa: {:?}", expa);
                for (i, v) in pos_self_inv.included_iter().zip(expa) {
                    exp_self.indices[i.0] = v;
                }

                for (i, v) in pos_other_inv.included_iter().zip(expb) {
                    exp_other.indices[i.0] = v;
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
        pos_self: SubSet<SlotIndex>,
        pos_other: SubSet<SlotIndex>,
        resulting_structure: <Self::LCM as HasStructure>::Structure,
        resulting_partition: SubSet<SlotIndex>,
    ) -> Result<Self::LCM, ContractionError> {
        let zero = self.data[0].try_upgrade().unwrap().into_owned().ref_zero();
        let mut result_data = vec![zero.clone(); resulting_structure.size()?];
        let pos_self_inv = !&pos_self;
        let pos_other_inv = !&pos_other;
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
                        .partition(|(i, _)| resulting_partition[SlotIndex(*i)]);

                // println!("expa: {:?}", expa);
                for (i, v) in pos_self_inv.included_iter().zip(expa) {
                    exp_self.indices[i.0] = v;
                }

                for (i, v) in pos_other_inv.included_iter().zip(expb) {
                    exp_other.indices[i.0] = v;
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
        pos_self: SubSet<SlotIndex>,
        pos_other: SubSet<SlotIndex>,
        resulting_structure: <Self::LCM as HasStructure>::Structure,
        resulting_partition: SubSet<SlotIndex>,
    ) -> Result<Self::LCM, ContractionError> {
        let mut result_data = HashMap::default();
        let zero = self.zero.try_upgrade().unwrap().as_ref().ref_zero();
        let self_slots = self
            .structure()
            .external_structure_iter()
            .collect::<Vec<_>>();
        let other_slots = other
            .structure()
            .external_structure_iter()
            .collect::<Vec<_>>();
        let matched_slots = pos_self
            .included_iter()
            .map(|self_position| {
                let self_slot = self_slots[self_position.0];
                let other_position = pos_other
                    .included_iter()
                    .find(|other_position| other_slots[other_position.0].dual() == self_slot)
                    .ok_or_else(|| {
                        ContractionError::Other(eyre::eyre!(
                            "missing matching interleaved sparse contraction slot for {self_slot}"
                        ))
                    })?;
                Ok((self_position, other_position))
            })
            .collect::<Result<Vec<_>, ContractionError>>()?;
        let self_free = (!&pos_self).included_iter().collect::<Vec<_>>();
        let other_free = (!&pos_other).included_iter().collect::<Vec<_>>();

        for (self_flat, self_value) in self.flat_iter() {
            let self_expanded = self.expanded_index(self_flat)?;
            for (other_flat, other_value) in other.flat_iter() {
                let other_expanded = other.expanded_index(other_flat)?;
                let mut neg = false;
                let mut matches = true;
                for (self_position, other_position) in &matched_slots {
                    let self_index = self_expanded.indices[self_position.0];
                    if self_index != other_expanded.indices[other_position.0] {
                        matches = false;
                        break;
                    }
                    if self_slots[self_position.0].rep().negative()?[self_index] {
                        neg = !neg;
                    }
                }
                if !matches {
                    continue;
                }

                let mut self_free_iter = self_free.iter();
                let mut other_free_iter = other_free.iter();
                let result_expanded = (0..resulting_structure.order())
                    .map(|position| {
                        if resulting_partition[SlotIndex(position)] {
                            let position = self_free_iter.next().ok_or_else(|| {
                                ContractionError::Other(eyre::eyre!(
                                    "missing self free slot while building interleaved sparse contraction result"
                                ))
                            })?;
                            Ok(self_expanded.indices[position.0])
                        } else {
                            let position = other_free_iter.next().ok_or_else(|| {
                                ContractionError::Other(eyre::eyre!(
                                    "missing other free slot while building interleaved sparse contraction result"
                                ))
                            })?;
                            Ok(other_expanded.indices[position.0])
                        }
                    })
                    .collect::<Result<ExpandedIndex, ContractionError>>()?;
                let result_index = resulting_structure.flat_index(&result_expanded)?;
                let value = result_data
                    .entry(result_index)
                    .or_insert_with(|| zero.clone());
                if neg {
                    value.sub_assign_fallible(&self_value.mul_fallible(other_value).unwrap());
                } else {
                    value.add_assign_fallible(&self_value.mul_fallible(other_value).unwrap());
                }
            }
        }
        result_data.retain(|_, value| !value.is_zero());
        let result = SparseTensor {
            zero,
            elements: result_data,
            structure: resulting_structure,
        };

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        contraction::Contract,
        structure::{
            OrderedStructure, PermutedStructure, StructureContract,
            representation::{LibraryRep, Lorentz, RepName},
            slot::{DualSlotTo, IsAbstractSlot},
        },
        tensors::data::{DenseTensor, GetTensorData, SetTensorData},
    };

    #[test]
    fn dense_interleaved_multi_contract_uses_other_free_slots() {
        let rep = Lorentz {};
        let self_free_left = rep.new_slot(2, 0).to_lib();
        let other_free = rep.new_slot(2, 2).to_lib();
        let self_free_right = rep.new_slot(2, 4).to_lib();
        let contracted_left = rep.new_slot(2, 10).to_lib();
        let contracted_right = rep.new_slot(2, 12).to_lib();

        let self_structure: OrderedStructure<LibraryRep> = PermutedStructure::from_iter([
            self_free_left,
            self_free_right,
            contracted_left,
            contracted_right,
        ])
        .structure;
        let other_structure: OrderedStructure<LibraryRep> = PermutedStructure::from_iter([
            other_free,
            contracted_left.dual(),
            contracted_right.dual(),
        ])
        .structure;

        let (resulting_structure, pos_self, pos_other, mergeinfo) =
            self_structure.merge(&other_structure).unwrap();
        println!("self structure: {self_structure}");
        println!("other structure: {other_structure}");
        println!("result structure: {resulting_structure}");
        println!("pos_self: {pos_self:?}");
        println!("pos_other: {pos_other:?}");
        println!("mergeinfo: {mergeinfo:?}");

        let mut lhs = DenseTensor::<i32, _>::zero(self_structure);
        let mut rhs = DenseTensor::<i32, _>::zero(other_structure);

        lhs.set(&[0, 0, 1, 1], 2).unwrap();
        rhs.set(&[1, 1, 1], 5).unwrap();

        let result = lhs.contract(&rhs).unwrap();
        println!("result data: {:?}", result.data);
        println!("result[0, 1, 0]: {}", result.get_ref(&[0, 1, 0]).unwrap());

        assert_eq!(*result.get_ref(&[0, 1, 0]).unwrap(), 10);
    }
}
