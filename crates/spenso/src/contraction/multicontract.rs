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
        slot::DualSlotTo,
    },
    tensors::data::{DataIterator, DenseTensor, SparseTensor},
};

use linnet::{
    half_edge::subgraph::{SubSetLike, subset::SubSet},
    permutation::Permutation,
};

use std::iter::Iterator;

use super::{ContractableWith, ContractionError, MultiContract, MultiContractInterleaved};

/// Matched slots of `left` and `right`, with the permutation that orders the
/// matched axes of a `left` fiber like the matched axes of a `right` fiber.
///
/// [`TensorStructure::match_indices`] returns the permutation that aligns the
/// right-hand matched axes with the left-hand order, which permuted left
/// fibers must apply inversely. The two coincide only for involutions, so a
/// direct use is wrong once three or more matched axes form a longer cycle.
pub(crate) fn matched_left_fibers<I: TensorStructure>(
    left: &I,
    right: &I,
) -> Option<(Permutation, Vec<bool>, Vec<bool>)> {
    left.match_indices(right)
        .map(|(permutation, left_matches, right_matches)| {
            (permutation.inverse(), left_matches, right_matches)
        })
}

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
            matched_left_fibers(self.structure(), other.structure()).unwrap();

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
            matched_left_fibers(self.structure(), other.structure()).unwrap();

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
            matched_left_fibers(self.structure(), other.structure()).unwrap();

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
            matched_left_fibers(self.structure(), other.structure()).unwrap();

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
        let contracted_slots_aligned = pos_self.included_iter().count()
            == pos_other.included_iter().count()
            && pos_self.included_iter().zip(pos_other.included_iter()).all(
                |(self_position, other_position)| {
                    let self_slot = self.structure().get_slot(self_position).unwrap();
                    let other_slot = other.structure().get_slot(other_position).unwrap();
                    self_slot == other_slot.dual()
                },
            );

        let pos_self_inv = !&pos_self;
        let pos_other_inv = !&pos_other;
        // println!("Interleaving dense dense multi"); // Initialize result tensor with default values
        let mut result_data = vec![zero.clone(); resulting_structure.size()?];

        let self_fiber_class = Fiber::from(&resulting_partition, &resulting_structure); //We use the partition as a filter here, for indices that belong to self, vs those that belong to other
        let (mut self_fiber_class_iter, mut other_fiber_class_iter) =
            CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class); // these are iterators over the open indices of self and other, except expressed in the flat indices of the resulting structure

        let mut exp_self = self.expanded_index(FlatIndex::from(0)).unwrap();
        let mut exp_other = other.expanded_index(FlatIndex::from(0)).unwrap();

        let mut iter_self = if contracted_slots_aligned {
            self.fiber(FiberData::from(&pos_self)).iter_metric()
        } else {
            let contraction_permutation = matched_left_fibers(self.structure(), other.structure())
                .map(|(permutation, _, _)| permutation)
                .ok_or_else(|| {
                    ContractionError::Other(eyre::eyre!(
                        "missing matching interleaved dense contraction slots"
                    ))
                })?;
            self.fiber(FiberData::from(&pos_self))
                .iter_perm_metric(contraction_permutation)
        };
        let mut iter_other = other.fiber(FiberData::from(&pos_other)).iter(); // same for other

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
        let contracted_slots_aligned = pos_self.included_iter().count()
            == pos_other.included_iter().count()
            && pos_self.included_iter().zip(pos_other.included_iter()).all(
                |(self_position, other_position)| {
                    let self_slot = self.structure().get_slot(self_position).unwrap();
                    let other_slot = other.structure().get_slot(other_position).unwrap();
                    self_slot == other_slot.dual()
                },
            );
        let pos_self_inv = !&pos_self;
        let pos_other_inv = !&pos_other;

        let mut result_data = vec![zero.clone(); resulting_structure.size()?];

        let self_fiber_class = Fiber::from(&resulting_partition, &resulting_structure); //We use the partition as a filter here, for indices that belong to self, vs those that belong to other
        let (mut self_fiber_class_iter, mut other_fiber_class_iter) =
            CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class); // these are iterators over the open indices of self and other, except expressed in the flat indices of the resulting structure
        let mut iter_self = if contracted_slots_aligned {
            self.fiber(FiberData::from(&pos_self)).iter_metric()
        } else {
            let contraction_permutation = matched_left_fibers(self.structure(), other.structure())
                .map(|(permutation, _, _)| permutation)
                .ok_or_else(|| {
                    ContractionError::Other(eyre::eyre!(
                        "missing matching interleaved sparse-dense contraction slots"
                    ))
                })?;
            self.fiber(FiberData::from(&pos_self))
                .iter_perm_metric(contraction_permutation)
        };
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
        let contracted_slots_aligned = pos_self.included_iter().count()
            == pos_other.included_iter().count()
            && pos_self.included_iter().zip(pos_other.included_iter()).all(
                |(self_position, other_position)| {
                    let self_slot = self.structure().get_slot(self_position).unwrap();
                    let other_slot = other.structure().get_slot(other_position).unwrap();
                    self_slot == other_slot.dual()
                },
            );
        let mut result_data = vec![zero.clone(); resulting_structure.size()?];
        let pos_self_inv = !&pos_self;
        let pos_other_inv = !&pos_other;
        let self_fiber_class = Fiber::from(&resulting_partition, &resulting_structure); //We use the partition as a filter here, for indices that belong to self, vs those that belong to other
        let (mut self_fiber_class_iter, mut other_fiber_class_iter) =
            CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class); // these are iterators over the open indices of self and other, except expressed in the flat indices of the resulting structure

        let mut iter_self = if contracted_slots_aligned {
            self.fiber(FiberData::from(&pos_self)).iter_metric()
        } else {
            let contraction_permutation = matched_left_fibers(self.structure(), other.structure())
                .map(|(permutation, _, _)| permutation)
                .ok_or_else(|| {
                    ContractionError::Other(eyre::eyre!(
                        "missing matching interleaved dense-sparse contraction slots"
                    ))
                })?;
            self.fiber(FiberData::from(&pos_self))
                .iter_perm_metric(contraction_permutation)
        };
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

        if self.flat_iter().next().is_some() {
            let contracted_slots_aligned = pos_self.included_iter().count()
                == pos_other.included_iter().count()
                && pos_self.included_iter().zip(pos_other.included_iter()).all(
                    |(self_position, other_position)| {
                        let self_slot = self.structure().get_slot(self_position).unwrap();
                        let other_slot = other.structure().get_slot(other_position).unwrap();
                        self_slot == other_slot.dual()
                    },
                );
            let pos_self_inv = !&pos_self;
            let pos_other_inv = !&pos_other;

            let self_fiber_class = Fiber::from(&resulting_partition, &resulting_structure);
            let (mut self_fiber_class_iter, mut other_fiber_class_iter) =
                CoreFlatFiberIterator::new_paired_conjugates(&self_fiber_class);

            let mut iter_self = if contracted_slots_aligned {
                self.fiber(FiberData::from(&pos_self)).iter_metric()
            } else {
                let contraction_permutation =
                    matched_left_fibers(self.structure(), other.structure())
                        .map(|(permutation, _, _)| permutation)
                        .ok_or_else(|| {
                            ContractionError::Other(eyre::eyre!(
                                "missing matching interleaved sparse contraction slots"
                            ))
                        })?;
                self.fiber(FiberData::from(&pos_self))
                    .iter_perm_metric(contraction_permutation)
            };
            let mut iter_other = other.fiber(FiberData::from(&pos_other)).iter();

            let mut exp_self = self.expanded_index(FlatIndex::from(0)).unwrap();
            let mut exp_other = other.expanded_index(FlatIndex::from(0)).unwrap();

            for fiber_class_a_id in self_fiber_class_iter.by_ref() {
                for fiber_class_b_id in other_fiber_class_iter.by_ref() {
                    let result_index = fiber_class_a_id + fiber_class_b_id;

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

                    let shift_a = self.structure().flat_index(&exp_self).unwrap();
                    let shift_b = other.structure().flat_index(&exp_other).unwrap();

                    iter_self.shift(shift_a.into());
                    iter_other.shift(shift_b.into());

                    let mut items = iter_self
                        .next()
                        .map(|(a, skip, (neg, _))| (a, skip, neg))
                        .zip(iter_other.next().map(|(b, skip, _)| (b, skip)));

                    let mut value = zero.clone();
                    let mut nonzero = false;

                    while let Some(((a, skip_a, neg), (b, skip_b))) = items {
                        if skip_a > skip_b {
                            let b = iter_other
                                .by_ref()
                                .next()
                                .map(|(b, skip, _)| (b, skip + skip_b + 1));
                            items = Some((a, skip_a, neg)).zip(b);
                        } else if skip_b > skip_a {
                            let a = iter_self
                                .by_ref()
                                .next()
                                .map(|(a, skip, (neg, _))| (a, skip + skip_a + 1, neg));
                            items = a.zip(Some((b, skip_b)));
                        } else {
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
                    if nonzero && value.is_non_zero() {
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

#[cfg(test)]
mod tests {
    use crate::{
        contraction::Contract,
        structure::{
            Canonicalized, OrderedStructure, StructureContract,
            representation::{LibraryRep, Lorentz, RepName},
            slot::{DualSlotTo, IsAbstractSlot},
        },
        tensors::data::{DenseTensor, GetTensorData, SetTensorData, SparseTensor},
    };

    #[test]
    fn dense_interleaved_multi_contract_uses_other_free_slots() {
        let rep = Lorentz {};
        let self_free_left = rep.new_slot(2, 0).to_lib();
        let other_free = rep.new_slot(2, 2).to_lib();
        let self_free_right = rep.new_slot(2, 4).to_lib();
        let contracted_left = rep.new_slot(2, 10).to_lib();
        let contracted_right = rep.new_slot(2, 12).to_lib();

        let self_structure: OrderedStructure<LibraryRep> = Canonicalized::from_iter([
            self_free_left,
            self_free_right,
            contracted_left,
            contracted_right,
        ])
        .into_canonical();
        let other_structure: OrderedStructure<LibraryRep> =
            Canonicalized::from_iter([other_free, contracted_left.dual(), contracted_right.dual()])
                .into_canonical();

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
        println!("result[0, 1, 0]: {}", result.get_ref([0, 1, 0]).unwrap());

        assert_eq!(*result.get_ref([0, 1, 0]).unwrap(), 10);
    }

    #[test]
    fn dense_interleaved_multi_contract_matches_dual_slots_across_base_dual_order() {
        let rep = Lorentz {};
        let self_free_left = rep.new_slot(2, 0).to_lib();
        let other_free = rep.new_slot(2, 2).to_lib();
        let self_free_right = rep.new_slot(2, 4).to_lib();
        let contracted_base = rep.new_slot(2, 10).to_lib();
        let contracted_other_base = rep.new_slot(2, 12).to_lib();

        let self_structure: OrderedStructure<LibraryRep> = Canonicalized::from_iter([
            self_free_left,
            self_free_right,
            contracted_base,
            contracted_other_base.dual(),
        ])
        .into_canonical();
        let other_structure: OrderedStructure<LibraryRep> =
            Canonicalized::from_iter([other_free, contracted_other_base, contracted_base.dual()])
                .into_canonical();

        let (resulting_structure, pos_self, pos_other, mergeinfo) =
            self_structure.merge(&other_structure).unwrap();
        println!("mixed-orientation dense self structure: {self_structure}");
        println!("mixed-orientation dense other structure: {other_structure}");
        println!("mixed-orientation dense result structure: {resulting_structure}");
        println!("mixed-orientation dense pos_self: {pos_self:?}");
        println!("mixed-orientation dense pos_other: {pos_other:?}");
        println!("mixed-orientation dense mergeinfo: {mergeinfo:?}");

        let mut lhs = DenseTensor::<i32, _>::zero(self_structure);
        let mut rhs = DenseTensor::<i32, _>::zero(other_structure);

        lhs.set(&[0, 0, 0, 1], 2).unwrap();
        rhs.set(&[1, 1, 0], 5).unwrap();

        let result = lhs.contract(&rhs).unwrap();
        println!("mixed-orientation dense result data: {:?}", result.data);
        println!(
            "mixed-orientation dense result[0, 1, 0]: {}",
            result.get_ref([0, 1, 0]).unwrap()
        );

        assert_eq!(*result.get_ref([0, 1, 0]).unwrap(), 10);
    }

    #[test]
    fn sparse_interleaved_multi_contract_matches_dual_slots_across_base_dual_order() {
        let rep = Lorentz {};
        let self_free_left = rep.new_slot(2, 0).to_lib();
        let other_free = rep.new_slot(2, 2).to_lib();
        let self_free_right = rep.new_slot(2, 4).to_lib();
        let contracted_base = rep.new_slot(2, 10).to_lib();
        let contracted_other_base = rep.new_slot(2, 12).to_lib();

        let self_structure: OrderedStructure<LibraryRep> = Canonicalized::from_iter([
            self_free_left,
            self_free_right,
            contracted_base,
            contracted_other_base.dual(),
        ])
        .into_canonical();
        let other_structure: OrderedStructure<LibraryRep> =
            Canonicalized::from_iter([other_free, contracted_other_base, contracted_base.dual()])
                .into_canonical();

        let (resulting_structure, pos_self, pos_other, mergeinfo) =
            self_structure.merge(&other_structure).unwrap();
        println!("mixed-orientation self structure: {self_structure}");
        println!("mixed-orientation other structure: {other_structure}");
        println!("mixed-orientation result structure: {resulting_structure}");
        println!("mixed-orientation pos_self: {pos_self:?}");
        println!("mixed-orientation pos_other: {pos_other:?}");
        println!("mixed-orientation mergeinfo: {mergeinfo:?}");

        let mut lhs = SparseTensor::<i32, _>::empty(self_structure, 0);
        let mut rhs = SparseTensor::<i32, _>::empty(other_structure, 0);

        lhs.set(&[0, 0, 0, 1], 2).unwrap();
        rhs.set(&[1, 1, 0], 5).unwrap();

        let result = lhs.contract(&rhs).unwrap();
        println!("mixed-orientation result elements: {:?}", result.elements);
        println!(
            "mixed-orientation result[0, 1, 0]: {}",
            result.get_ref([0, 1, 0]).unwrap()
        );

        assert_eq!(*result.get_ref([0, 1, 0]).unwrap(), 10);
    }

    mod matched_slot_cycles {
        use std::collections::BTreeMap;

        use crate::{
            contraction::Contract,
            structure::{
                Canonicalized, HasStructure, OrderedStructure, TensorStructure,
                abstract_index::AbstractIndex,
                representation::{Euclidean, LibraryRep, LibrarySlot, Lorentz, Minkowski, RepName},
                slot::IsAbstractSlot,
            },
            tensors::data::{DataTensor, DenseTensor, GetTensorData, SetTensorData},
        };

        type Slot = LibrarySlot<AbstractIndex>;
        type Tensor = DenseTensor<i32, OrderedStructure<LibraryRep>>;

        fn up(label: usize) -> Slot {
            Lorentz {}.new_slot(2, label).to_lib()
        }
        fn down(label: usize) -> Slot {
            Lorentz {}.dual().new_slot(2, label).to_lib()
        }
        fn mink(label: usize) -> Slot {
            Minkowski {}.new_slot(2, label).to_lib()
        }
        fn euc(label: usize) -> Slot {
            Euclidean {}.new_slot(2, label).to_lib()
        }

        /// A component determined by the index labels alone, so the reference
        /// is independent of canonical storage order. It is not affine in the
        /// indices, so a wrong pairing of matched axes cannot cancel out.
        fn value(seed: i32, labels: &BTreeMap<usize, usize>) -> i32 {
            labels.iter().fold(seed, |acc, (label, index)| {
                acc * (2 + *index as i32) + (*label as i32 % 5 + seed) * (*index as i32 + 1)
            })
        }

        fn labels(slots: &[Slot], indices: &[usize]) -> BTreeMap<usize, usize> {
            slots
                .iter()
                .map(|slot| usize::from(slot.aind()))
                .zip(indices.iter().copied())
                .collect()
        }

        fn assignments(rank: usize) -> impl Iterator<Item = Vec<usize>> {
            (0..1_usize << rank).map(move |bits| (0..rank).map(|axis| (bits >> axis) & 1).collect())
        }

        fn tensor(slots: &[Slot], seed: i32) -> Tensor {
            let structure: OrderedStructure<LibraryRep> =
                Canonicalized::from_iter(slots.iter().copied()).into_canonical();
            let stored = structure.external_structure();
            let mut tensor = Tensor::zero(structure);
            for indices in assignments(stored.len()) {
                tensor
                    .set(&indices, value(seed, &labels(&stored, &indices)))
                    .unwrap();
            }
            tensor
        }

        /// `sum` over the shared labels of the two label-defined factors, with
        /// the Minkowski metric diag(+, -) on shared `mink` labels. Each factor
        /// is given with the seed its components were built from.
        fn reference(
            (left, left_seed): (&[Slot], i32),
            (right, right_seed): (&[Slot], i32),
            free: &BTreeMap<usize, usize>,
        ) -> i32 {
            let shared = left
                .iter()
                .filter(|slot| right.iter().any(|other| other.aind() == slot.aind()))
                .copied()
                .collect::<Vec<_>>();
            assignments(shared.len())
                .map(|indices| {
                    let mut all = free.clone();
                    let mut sign = 1;
                    for (slot, index) in shared.iter().zip(&indices) {
                        all.insert(usize::from(slot.aind()), *index);
                        if *slot == mink(usize::from(slot.aind())) && *index > 0 {
                            sign = -sign;
                        }
                    }
                    let restrict = |slots: &[Slot]| {
                        slots
                            .iter()
                            .map(|slot| {
                                let label = usize::from(slot.aind());
                                (label, all[&label])
                            })
                            .collect::<BTreeMap<_, _>>()
                    };
                    sign * value(left_seed, &restrict(left)) * value(right_seed, &restrict(right))
                })
                .sum()
        }

        fn check(left: (&[Slot], i32), right: (&[Slot], i32), result: &Tensor) {
            let stored = result.external_structure();
            for indices in assignments(stored.len()) {
                let free = labels(&stored, &indices);
                assert_eq!(
                    *result.get_ref(&indices).unwrap(),
                    reference(left, right, &free),
                    "{left:?} . {right:?} at {free:?}"
                );
            }
        }

        fn contract_all_storages(left: &[Slot], right: &[Slot]) {
            let (a, b) = (tensor(left, 1), tensor(right, 2));
            // The matched dual slots must form a cycle longer than two, where
            // a permutation and its inverse differ.
            let (permutation, _, _) = a.structure().match_indices(b.structure()).unwrap();
            assert_ne!(permutation, permutation.inverse());
            let storages = |tensor: &Tensor| {
                [
                    DataTensor::Dense(tensor.clone()),
                    DataTensor::Sparse(tensor.to_sparse()),
                ]
            };
            for lhs in storages(&a) {
                for rhs in storages(&b) {
                    for (result, l, r) in [
                        (lhs.contract(&rhs).unwrap(), (left, 1), (right, 2)),
                        (rhs.contract(&lhs).unwrap(), (right, 2), (left, 1)),
                    ] {
                        let dense = match result {
                            DataTensor::Dense(dense) => dense,
                            DataTensor::Sparse(sparse) => sparse.to_dense(),
                        };
                        check(l, r, &dense);
                    }
                }
            }
        }

        #[test]
        fn three_cycle_of_dual_pairs_beside_a_free_index() {
            contract_all_storages(&[up(1), up(2), down(3), euc(9)], &[up(3), down(1), down(2)]);
        }

        #[test]
        fn three_cycle_of_dual_pairs_with_interleaved_free_indices() {
            contract_all_storages(
                &[up(1), up(2), down(3), euc(5), euc(9)],
                &[up(3), down(1), down(2), euc(7)],
            );
        }

        /// Canonical order sorts by representation first, so a self-dual
        /// metric pair is a fixed point of the cycle; its sign must still
        /// follow its own axis.
        #[test]
        fn three_cycle_beside_a_metric_pair() {
            contract_all_storages(
                &[up(1), up(2), down(3), mink(4), euc(9)],
                &[up(3), down(1), down(2), mink(4)],
            );
        }

        #[test]
        fn fully_contracted_three_cycle() {
            contract_all_storages(&[up(1), up(2), down(3)], &[up(3), down(1), down(2)]);
        }
    }
}
