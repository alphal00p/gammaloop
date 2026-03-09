use ahash::AHashMap;
use eyre::eyre;
use eyre::{Result};
use bitvec::vec::BitVec;
use concrete_index::ConcreteIndex;
use concrete_index::DualConciousExpandedIndex;
use concrete_index::DualConciousIndex;
use concrete_index::ExpandedIndex;
use concrete_index::FlatIndex;
use delegate::delegate;
use dimension::Dimension;

use thiserror::Error;

#[cfg(feature = "shadowing")]
use crate::structure::slot::ParseableAind;
#[cfg(feature = "shadowing")]
use crate::{
    shadowing::symbolica_utils::{IntoArgs, IntoSymbol},
    tensors::data::DenseTensor,
    tensors::parametric::{ExpandedCoefficent, FlatCoefficent, TensorCoefficient},
};
use linnet::permutation::Permutation;
use representation::{RepName, Representation};
use serde::Deserialize;
use serde::Serialize;
use slot::DualSlotTo;
use slot::IsAbstractSlot;
use slot::SlotError;
use std::fmt::Debug;
#[cfg(feature = "shadowing")]
use symbolica::atom::{Atom, FunctionBuilder, Symbol};

use crate::iterators::TensorStructureIndexIterator;
use std::collections::HashMap;
use std::collections::HashSet;

pub mod abstract_index;
pub mod concrete_index;
pub mod dimension;
pub mod indexless;
pub use indexless::{IndexLess, IndexlessNamedStructure};
pub mod named;
pub use named::NamedStructure;
pub mod permuted;
pub use permuted::PermutedStructure;
pub mod ordered;
pub mod representation;
pub mod slot;
pub use ordered::OrderedStructure;

pub mod smart_shadow;
pub use smart_shadow::SmartShadowStructure;

pub trait ScalarTensor: HasStructure<Structure: ScalarStructure> {
    fn new_scalar(scalar: Self::Scalar) -> Self;
}
pub trait HasStructure {
    type Store<S>: HasStructure<Structure = S>
    where
        S: TensorStructure;

    type Structure: TensorStructure;
    type Scalar;
    type ScalarRef<'a>
    where
        Self: 'a;

    fn map_structure<O: TensorStructure>(self, f: impl Fn(Self::Structure) -> O) -> Self::Store<O>;

    fn map_structure_result<O: TensorStructure, Er>(
        self,
        f: impl Fn(Self::Structure) -> Result<O, Er>,
    ) -> std::result::Result<Self::Store<O>, Er>;
    fn structure(&self) -> &Self::Structure;
    fn mut_structure(&mut self) -> &mut Self::Structure;
    fn scalar(self) -> Option<Self::Scalar>;
    fn scalar_ref(&self) -> Option<Self::ScalarRef<'_>>;
    fn map_same_structure(self, f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self;

    fn set_structure_name<N>(&mut self, name: N)
    where
        Self::Structure: HasName<Name = N>,
    {
        self.mut_structure().set_name(name);
    }
    fn structure_name(&self) -> Option<<Self::Structure as HasName>::Name>
    where
        Self::Structure: HasName,
    {
        self.structure().name()
    }
    fn structure_id(&self) -> Option<<Self::Structure as HasName>::Args>
    where
        Self::Structure: HasName,
    {
        self.structure().args()
    }
    // fn cast_structure<O, S>(self) -> O
    // where
    //     O: HasStructure<Structure = S, Scalar = Self::Scalar>,
    //     S: TensorStructure + From<Self::Structure>;
}

pub trait CastStructure<O: HasStructure<Structure: From<Self::Structure>>>: HasStructure {
    fn cast_structure(self) -> O;
}

#[cfg(feature = "shadowing")]
impl<T: HasName> ToSymbolic for PermutedStructure<T>
where
    T: TensorStructure,
    T::Name: IntoSymbol,
    <T::Slot as IsAbstractSlot>::Aind: ParseableAind,
    T::Args: IntoArgs,
{
    fn concrete_atom(&self, id: FlatIndex) -> ExpandedCoefficent<()> {
        ExpandedCoefficent {
            name: self.structure.name().map(|n| n.ref_into_symbol()),
            index: self.structure.co_expanded_index(id).unwrap(),
            args: None,
        }
    }

    fn flat_atom(&self, id: FlatIndex) -> FlatCoefficent<()> {
        FlatCoefficent {
            name: self.structure.name().map(|n| n.ref_into_symbol()),
            index: id,
            args: None,
        }
    }

    fn to_dense_labeled<R>(
        self,
        index_to_atom: impl Fn(&Self, FlatIndex) -> R,
    ) -> Result<DenseTensor<Atom, Self>>
    where
        Self: Sized,
        R: TensorCoefficient,
    {
        let mut data = vec![];
        for index in 0..self.structure.size()? {
            data.push(index_to_atom(&self, index.into()).to_atom().unwrap());
        }

        Ok(DenseTensor {
            data,
            structure: self,
        })
    }

    fn to_dense_labeled_complex<R>(
        self,
        index_to_atom: impl Fn(&Self, FlatIndex) -> R,
    ) -> Result<DenseTensor<Atom, Self>>
    where
        Self: Sized,
        R: TensorCoefficient,
    {
        let mut data = vec![];
        for index in 0..self.structure.size()? {
            let re = index_to_atom(&self, index.into()).to_atom_re().unwrap();
            let im = index_to_atom(&self, index.into()).to_atom_im().unwrap();
            let i = Atom::i();
            data.push(&re + i * &im);
        }

        Ok(DenseTensor {
            data,
            structure: self,
        })
    }

    fn to_symbolic_with(
        &self,
        name: Symbol,
        args: &[Atom],
        permutation: Option<Permutation>,
    ) -> Atom {
        let mut slots = self
            .structure
            .external_structure_iter()
            .map(|slot| slot.to_atom())
            .collect::<Vec<_>>();
        if let Some(perm) = permutation {
            perm.apply_slice_in_place(&mut slots);
        }
        FunctionBuilder::new(name.ref_into_symbol())
            .add_args(args)
            .add_args(&slots)
            .finish()
    }
}

#[cfg(feature = "shadowing")]
impl<T: HasName> ToSymbolic for T
where
    T: TensorStructure,
    <T::Slot as IsAbstractSlot>::Aind: ParseableAind,
    T::Name: IntoSymbol,
    T::Args: IntoArgs,
{
    fn concrete_atom(&self, id: FlatIndex) -> ExpandedCoefficent<()> {
        ExpandedCoefficent {
            name: self.name().map(|n| n.ref_into_symbol()),
            index: self.co_expanded_index(id).unwrap(),
            args: None,
        }
    }

    fn flat_atom(&self, id: FlatIndex) -> FlatCoefficent<()> {
        FlatCoefficent {
            name: self.name().map(|n| n.ref_into_symbol()),
            index: id,
            args: None,
        }
    }

    fn to_dense_labeled<R>(
        self,
        index_to_atom: impl Fn(&Self, FlatIndex) -> R,
    ) -> Result<DenseTensor<Atom, Self>>
    where
        Self: Sized,
        R: TensorCoefficient,
    {
        let mut data = vec![];
        for index in 0..self.size()? {
            data.push(index_to_atom(&self, index.into()).to_atom().unwrap());
        }

        Ok(DenseTensor {
            data,
            structure: self,
        })
    }

    fn to_dense_labeled_complex<R>(
        self,
        index_to_atom: impl Fn(&Self, FlatIndex) -> R,
    ) -> Result<DenseTensor<Atom, Self>>
    where
        Self: Sized,
        R: TensorCoefficient,
    {
        let mut data = vec![];
        for index in 0..self.size()? {
            let re = index_to_atom(&self, index.into()).to_atom_re().unwrap();
            let im = index_to_atom(&self, index.into()).to_atom_im().unwrap();
            let i = Atom::i();
            data.push(&re + i * &im);
        }

        Ok(DenseTensor {
            data,
            structure: self,
        })
    }

    fn to_symbolic_with(&self, name: Symbol, args: &[Atom], perm: Option<Permutation>) -> Atom {
        let mut slots = self
            .external_structure_iter()
            .map(|slot| slot.to_atom())
            .collect::<Vec<_>>();
        if let Some(perm) = perm {
            perm.apply_slice_in_place(&mut slots);
        }

        FunctionBuilder::new(name.ref_into_symbol())
            .add_args(args)
            .add_args(&slots)
            .finish()
    }
}

#[cfg(feature = "shadowing")]
pub trait ToSymbolic {
    fn concrete_atom(&self, id: FlatIndex) -> ExpandedCoefficent<()>;
    // {
    //     ExpandedCoefficent {
    //         name: None,
    //         index: self.co_expanded_index(id).unwrap(),
    //         args: None,
    //     }
    // }

    fn flat_atom(&self, id: FlatIndex) -> FlatCoefficent<()> {
        FlatCoefficent {
            name: None,
            index: id,
            args: None,
        }
    }

    fn to_dense_expanded_labels(self) -> Result<DenseTensor<Atom, Self>>
    where
        Self: std::marker::Sized + Clone,
    {
        self.to_dense_labeled(Self::concrete_atom)
    }

    fn to_dense_labeled<T>(
        self,
        index_to_atom: impl Fn(&Self, FlatIndex) -> T,
    ) -> Result<DenseTensor<Atom, Self>>
    where
        Self: Sized,
        T: TensorCoefficient;

    fn to_dense_labeled_complex<T>(
        self,
        index_to_atom: impl Fn(&Self, FlatIndex) -> T,
    ) -> Result<DenseTensor<Atom, Self>>
    where
        Self: Sized,
        T: TensorCoefficient;

    fn to_dense_flat_labels(self) -> Result<DenseTensor<Atom, Self>>
    where
        Self: std::marker::Sized + Clone,
    {
        self.to_dense_labeled(Self::flat_atom)
    }

    fn to_symbolic(&self, perm: Option<Permutation>) -> Option<Atom>
    where
        Self: HasName<Name: IntoSymbol, Args: IntoArgs>,
    {
        let args = self.args().map(|s| s.args()).unwrap_or_default();

        Some(self.to_symbolic_with(self.name()?.ref_into_symbol(), &args, perm))
    }

    fn to_symbolic_with(&self, name: Symbol, args: &[Atom], perm: Option<Permutation>) -> Atom;
}

pub trait ScalarStructure {
    fn scalar_structure() -> Self;
}
pub trait TensorStructure {
    type Slot: IsAbstractSlot + DualSlotTo<Dual = Self::Slot>;
    type Indexed: TensorStructure<Indexed = Self::Indexed, Slot = Self::Slot>;
    // type R: Rep;
    //

    fn reindex(
        self,
        indices: &[<Self::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::Indexed>, StructureError>;
    fn dual(self) -> Self;

    fn string_rep(&self) -> String {
        self.external_structure_iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>()
            .join("")
    }

    fn external_structure_iter(&self) -> impl Iterator<Item = Self::Slot>;
    fn external_dims_iter(&self) -> impl Iterator<Item = Dimension>;
    fn external_reps_iter(
        &self,
    ) -> impl Iterator<Item = Representation<<Self::Slot as IsAbstractSlot>::R>>;

    fn external_indices_iter(&self) -> impl Iterator<Item = <Self::Slot as IsAbstractSlot>::Aind>;
    fn get_aind(&self, i: usize) -> Option<<Self::Slot as IsAbstractSlot>::Aind>;
    fn get_rep(&self, i: usize) -> Option<Representation<<Self::Slot as IsAbstractSlot>::R>>;
    fn get_dim(&self, i: usize) -> Option<Dimension>;
    fn get_slot(&self, i: usize) -> Option<Self::Slot>;
    fn order(&self) -> usize;
    /// returns the list of slots that are the external indices of the tensor
    fn external_structure(&self) -> Vec<Self::Slot> {
        self.external_structure_iter().collect()
    }

    fn to_shell(self) -> TensorShell<Self>
    where
        Self: Sized,
    {
        TensorShell::new(self)
    }

    fn contains_matching(&self, slot: &Self::Slot) -> bool {
        self.external_structure_iter().any(|s| s.matches(slot))
    }

    fn external_reps(&self) -> Vec<Representation<<Self::Slot as IsAbstractSlot>::R>> {
        self.external_reps_iter().collect()
    }

    fn external_indices(&self) -> Vec<<Self::Slot as IsAbstractSlot>::Aind> {
        self.external_indices_iter().collect()
    }

    // fn iter_index_along_fiber(&self,fiber_position: &[bool]  )-> TensorStructureMultiFiberIterator where Self: Sized{
    //     TensorStructureMultiFiberIterator::new(self, fiber_position)
    // }

    // fn single_fiber_at(&self,fiber_pos:usize)->Fiber{
    //     let mut  f =Fiber::zeros(self.external_structure().len());
    //     f.free(fiber_pos);
    //     f.is_single();
    //     f
    // }

    /// checks if the tensor has the same exact structure as another tensor
    fn same_content(&self, other: &Self) -> bool {
        self.same_external(other)
    }

    /// Given two [`TensorStructure`]s, returns the index of the first matching slot in each external index list, along with a boolean indicating if there is a single match
    fn match_index(&self, other: &Self) -> Option<(bool, usize, usize)> {
        let posmap = self
            .external_structure_iter()
            .enumerate()
            .map(|(i, slot)| (slot, i))
            .collect::<AHashMap<_, _>>();

        let mut first_pair: Option<(usize, usize)> = None;

        for (j, slot) in other.external_structure_iter().enumerate() {
            if let Some(&i) = posmap.get(&slot.dual()) {
                if let Some((i, j)) = first_pair {
                    // Found a second match, return early with false indicating non-unique match
                    return Some((false, i, j));
                }
                first_pair = Some((i, j));
            }
        }

        first_pair.map(|(i, j)| (true, i, j)) // Maps the found pair to Some with true indicating a unique match, or None if no match was found
    }

    /// Given two [`TensorStructure`]s, returns the index of the first matching slot in each external index list
    fn match_indices(&self, other: &Self) -> Option<(Permutation, Vec<bool>, Vec<bool>)> {
        let mut self_matches = vec![false; self.order()];
        let mut perm = Vec::new();
        let mut other_matches = vec![false; other.order()];

        let posmap = self
            .external_structure_iter()
            .enumerate()
            .map(|(i, slot)| (slot, i))
            .collect::<AHashMap<_, _>>();

        for (j, slot_other) in other.external_structure_iter().enumerate() {
            if let Some(&i) = posmap.get(&slot_other.dual()) {
                self_matches[i] = true;
                other_matches[j] = true;
                perm.push(i);
            }
        }

        if perm.is_empty() {
            None
        } else {
            let p: Permutation = Permutation::sort(&perm);
            Some((p, self_matches, other_matches))
        }
    }
    /// Identify the repeated slots in the external index list
    fn traces(&self) -> Vec<[usize; 2]> {
        let mut positions: HashMap<<Self as TensorStructure>::Slot, Vec<usize>> = HashMap::new();

        // Track the positions of each element
        for (index, key) in self.external_structure_iter().enumerate() {
            if let Some(v) = positions.get_mut(&key.dual()) {
                v.push(index);
            } else {
                positions.insert(key, vec![index]);
            }
        }

        // Collect only the positions of repeated elements
        positions
            .into_iter()
            .filter_map(|(_, indices)| {
                if indices.len() == 2 {
                    Some([indices[0], indices[1]])
                } else {
                    None
                }
            })
            .collect()
    }

    /// yields the (outwards facing) shape of the tensor as a list of dimensions
    fn shape(&self) -> Vec<Dimension> {
        self.external_dims_iter().collect()
    }

    fn reps(&self) -> Vec<Representation<<Self::Slot as IsAbstractSlot>::R>> {
        self.external_reps_iter().collect()
    }

    /// checks if externally, the two tensors are the same
    fn same_external(&self, other: &Self) -> bool {
        let set1: HashSet<_> = self.external_structure_iter().collect();
        let set2: HashSet<_> = other.external_structure_iter().collect();
        set1 == set2
    }

    /// find the permutation of the external indices that would make the two tensors the same. Applying the permutation to other should make it the same as self
    fn find_permutation(&self, other: &Self) -> Result<Permutation> {
        if self.order() != other.order() {
            return Err(eyre!(
                "Mismatched order: {} vs {}",
                self.order(),
                other.order()
            ));
        }
        let other_structure = other.external_structure();
        let self_structure = self.external_structure();

        let other_sort = Permutation::sort(&other_structure);
        let self_sort = Permutation::sort(&self_structure);

        if other_sort.apply_slice(&other_structure) == self_sort.apply_slice(&self_structure) {
            Ok(other_sort.compose(&self_sort.inverse()))
        } else {
            Err(eyre!("Mismatched structure"))
        }

        // let mut index_map = HashMap::new();
        // for (i, item) in other.external_structure_iter().enumerate() {
        //     index_map.entry(item).or_insert_with(Vec::new).push(i);
        // }

        // let mut permutation = Vec::with_capacity(self.order());
        // let mut used_indices = HashSet::new();
        // for item in self.external_structure_iter() {
        //     if let Some(indices) = index_map.get_mut(&item) {
        //         // Find an index that hasn't been used yet
        //         if let Some(&index) = indices.iter().find(|&&i| !used_indices.contains(&i)) {
        //             permutation.push(index);
        //             used_indices.insert(index);
        //         } else {
        //             // No available index for this item
        //             return Err(eyre!("No available index for {:?}", item));
        //         }
        //     } else {
        //         // Item not found in other
        //         return Err(eyre!("Item {:?} not found in other", item));
        //     }
        // }

        // Ok(permutation)
    }

    /// yields the strides of the tensor in column major order
    fn strides_column_major(&self) -> Result<Vec<usize>> {
        let mut strides: Vec<usize> = vec![1; self.order()];

        if self.order() == 0 {
            return Ok(strides);
        }

        for i in 0..self.order() - 1 {
            strides[i + 1] = strides[i] * usize::try_from(self.shape()[i])?;
        }

        Ok(strides)
    }

    /// yields the strides of the tensor in row major order
    fn strides_row_major(&self) -> Result<Vec<usize>> {
        let mut strides = vec![1; self.order()];
        if self.order() == 0 {
            return Ok(strides);
        }

        for i in (0..self.order() - 1).rev() {
            strides[i] = strides[i + 1] * usize::try_from(self.shape()[i + 1])?;
        }

        Ok(strides)
    }

    /// By default, the strides are row major
    fn strides(&self) -> Result<Vec<usize>> {
        self.strides_row_major()
    }

    /// Verifies that the list of indices provided are valid for the tensor
    ///
    /// # Errors
    ///
    /// `Mismatched order` = if the length of the indices is different from the order of the tensor,
    ///
    /// `Index out of bounds` = if the index is out of bounds for the dimension of that index
    ///
    fn verify_indices<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<()> {
        if indices.as_ref().len() != self.order() {
            return Err(eyre!(
                "Mismatched order: {} indices, vs order {}",
                indices.as_ref().len(),
                self.order()
            ));
        }

        for (i, dim_len) in self
            .external_structure_iter()
            .map(|slot| slot.dim())
            .enumerate()
        {
            if indices.as_ref()[i] >= usize::try_from(dim_len)? {
                return Err(eyre!(
                    "Index {} out of bounds for dimension {} of size {}",
                    indices.as_ref()[i],
                    i,
                    usize::try_from(dim_len)?
                ));
            }
        }
        Ok(())
    }

    /// yields the flat index of the tensor given a list of indices
    ///
    /// # Errors
    ///
    /// Same as [`Self::verify_indices`]
    fn flat_index<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<FlatIndex> {
        let strides = self.strides()?;
        self.verify_indices(&indices)?;

        let mut idx = 0;
        for (i, &index) in indices.as_ref().iter().enumerate() {
            idx += index * strides[i];
        }
        Ok(idx.into())
    }

    /// yields the expanded index of the tensor given a flat index
    ///
    /// # Errors
    ///
    /// `Index out of bounds` = if the flat index is out of bounds for the tensor
    fn expanded_index(&self, flat_index: FlatIndex) -> Result<ExpandedIndex> {
        let mut indices = vec![];
        let mut index: usize = flat_index.into();
        for &stride in &self.strides()? {
            indices.push(index / stride);
            index %= stride;
        }
        if usize::from(flat_index) < self.size()? {
            Ok(indices.into())
        } else {
            Err(eyre!("Index {flat_index} out of bounds"))
        }
    }

    fn co_expanded_index(&self, flat_index: FlatIndex) -> Result<DualConciousExpandedIndex> {
        let mut indices = vec![];

        for (r, i) in self
            .external_reps_iter()
            .zip(self.expanded_index(flat_index)?.iter())
        {
            if r.rep.is_base() && r.rep.is_dual() {
                indices.push(DualConciousIndex::SelfDual(*i));
            } else if r.rep.is_base() {
                indices.push(DualConciousIndex::Up(*i));
            } else {
                indices.push(DualConciousIndex::Down(*i));
            }
        }
        Ok(indices.into())
    }

    /// yields an iterator over the indices of the tensor
    fn index_iter(&self) -> TensorStructureIndexIterator<'_, Self>
    where
        Self: Sized,
    {
        TensorStructureIndexIterator::new(self)
    }

    /// if the tensor has no (external) indices, it is a scalar
    fn is_scalar(&self) -> bool {
        self.order() == 0
    }

    fn is_fully_self_dual(&self) -> bool;

    // /// get the metric along the i-th index
    // fn get_ith_metric(&self, i: usize) -> Result<Vec<bool>> {
    //     self.get_rep(i)
    //         .ok_or(eyre!("out of bounds access"))?
    //         .negative()
    // }

    /// yields the size of the tensor, i.e. the product of the dimensions. This is the length of the vector of the data in a dense tensor
    fn size(&self) -> Result<usize> {
        if self.order() == 0 {
            return Ok(1);
        }
        let mut size = 1;
        for dim in self.shape() {
            size *= usize::try_from(dim)?;
        }
        Ok(size)
    }
}

pub trait TracksCount {
    fn contractions_num(&self) -> usize;

    fn is_composite(&self) -> bool {
        self.contractions_num() > 0
    }
}

/// An enum to represent information about a clean split in merged structures
#[derive(Debug, PartialEq, Eq)]
pub enum MergeInfo {
    /// Clean merge where all elements of first structure come before second structure
    FirstBeforeSecond,
    /// Clean merge where all elements of second structure come before first structure
    SecondBeforeFirst,
    /// No clean merge, elements are interleaved with the following partition (first is true)
    Interleaved(BitVec),
}

impl From<BitVec> for MergeInfo {
    fn from(bitvec: BitVec) -> Self {
        let mut n_transitions = 0;
        if bitvec.is_empty() {
            return MergeInfo::FirstBeforeSecond;
        }
        for i in 0..(bitvec.len() - 1) {
            if bitvec[i] != bitvec[i + 1] {
                n_transitions += 1;
            }

            if n_transitions > 1 {
                return MergeInfo::Interleaved(bitvec);
            }
        }

        if let Some(first) = bitvec.first() {
            if *first {
                MergeInfo::FirstBeforeSecond
            } else {
                MergeInfo::SecondBeforeFirst
            }
        } else {
            MergeInfo::FirstBeforeSecond
        }
    }
}

/// A trait for a structure that can be traced and merged, during a contraction.
pub trait StructureContract: Sized {
    fn trace(&mut self, i: usize, j: usize);

    fn trace_out(&mut self);

    /// Merges the two structures, yielding a new ordered structure.
    ///
    /// Assumes both are ordered structures
    ///
    /// Errors if self or other contained non traced out indices
    ///
    /// Otherwise returns a new structure without any common indices, the positions of the common indices self, and the positions in other, along with a MergeInfo specifying whether the merger is clean or not.
    ///
    fn merge(&self, other: &Self) -> Result<(Self, BitVec, BitVec, MergeInfo), StructureError>;

    // fn concat(&mut self, other: Self);

    // #[must_use]
    // fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self;
}

#[derive(Error, Debug)]
pub enum StructureError {
    #[error("SlotError: {0}")]
    SlotError(#[from] SlotError),
    #[error("empty structure {0}")]
    EmptyStructure(SlotError),
    #[error("wrong number of arguments {0}, expected {1}")]
    WrongNumberOfArguments(usize, usize),
    // #[error("Non traced out indices before merger {0}")]
    // NonTracedOut(#[from] DuplicateItemError),
    #[error("Parsing error: expected function view found {0}")]
    ParsingError(String),
}

pub trait HasName {
    type Name: Clone;
    type Args: Clone;
    fn name(&self) -> Option<Self::Name>;
    fn args(&self) -> Option<Self::Args>;
    fn set_name(&mut self, name: Self::Name);

    #[cfg(feature = "shadowing")]
    fn expanded_coef(&self, id: FlatIndex) -> ExpandedCoefficent<Self::Args>
    where
        Self: TensorStructure,
        Self::Name: IntoSymbol,
        Self::Args: IntoArgs,
    {
        ExpandedCoefficent {
            name: self.name().map(|n| n.ref_into_symbol()),
            index: self.co_expanded_index(id).unwrap(),
            args: self.args(),
        }
    }
    #[cfg(feature = "shadowing")]
    fn expanded_coef_perm(
        &self,
        id: FlatIndex,
        permutation: &Permutation,
    ) -> ExpandedCoefficent<Self::Args>
    where
        Self: TensorStructure,

        Self::Name: IntoSymbol,
        Self::Args: IntoArgs,
    {
        let mut index = self.co_expanded_index(id).unwrap();
        index.permute(permutation);
        ExpandedCoefficent {
            name: self.name().map(|n| n.ref_into_symbol()),
            index,
            args: self.args(),
        }
    }

    #[cfg(feature = "shadowing")]
    fn flat_coef(&self, id: FlatIndex) -> FlatCoefficent<Self::Args>
    where
        Self: TensorStructure,
        Self::Name: IntoSymbol,
        Self::Args: IntoArgs,
    {
        FlatCoefficent {
            name: self.name().map(|n| n.ref_into_symbol()),
            index: id,
            args: self.args(),
        }
    }
}

#[derive(
    Clone,
    PartialEq,
    Eq,
    Debug,
    Serialize,
    Deserialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
)]
pub struct TensorShell<S> {
    pub(crate) structure: S,
}

impl<S: TensorStructure + ScalarStructure> ScalarTensor for TensorShell<S> {
    fn new_scalar(_scalar: Self::Scalar) -> Self {
        TensorShell {
            structure: S::scalar_structure(),
        }
    }
}

impl<T> TensorStructure for TensorShell<T>
where
    T: TensorStructure,
{
    // type R = <T::Structure as TensorStructure>::R;
    type Indexed = TensorShell<T::Indexed>;
    type Slot = T::Slot;

    fn reindex(
        self,
        indices: &[<Self::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::Indexed>, StructureError> {
        let res = self.structure.reindex(indices)?;
        Ok(PermutedStructure {
            structure: TensorShell {
                structure: res.structure,
            },
            index_permutation: res.index_permutation,
            rep_permutation: res.rep_permutation,
        })
    }

    fn dual(self) -> Self {
        self.map_same_structure(|s| s.dual())
    }

    delegate! {
        to self.structure() {
            fn is_fully_self_dual(&self)-> bool;
            fn external_reps_iter(&self)-> impl Iterator<Item = Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn external_indices_iter(&self)-> impl Iterator<Item = <Self::Slot as IsAbstractSlot>::Aind>;
            fn external_dims_iter(&self)-> impl Iterator<Item = Dimension>;
            fn external_structure_iter(&self)-> impl Iterator<Item = Self::Slot>;
            fn get_slot(&self, i: usize)-> Option<Self::Slot>;
            fn get_rep(&self, i: usize)-> Option<Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn get_dim(&self, i: usize)-> Option<Dimension>;
            fn get_aind(&self, i: usize)-> Option<<Self::Slot as IsAbstractSlot>::Aind>;
            fn order(&self)-> usize;
        }
    }
}

impl<S: TensorStructure> HasStructure for TensorShell<S> {
    type Structure = S;
    type Scalar = ();
    type ScalarRef<'a>
        = &'a ()
    where
        Self: 'a;
    type Store<U>
        = TensorShell<U>
    where
        U: TensorStructure;

    fn map_structure<O: TensorStructure>(
        self,
        f: impl FnOnce(Self::Structure) -> O,
    ) -> Self::Store<O> {
        TensorShell {
            structure: f(self.structure),
        }
    }

    fn map_structure_result<O: TensorStructure, Er>(
        self,
        f: impl FnOnce(Self::Structure) -> Result<O, Er>,
    ) -> std::result::Result<Self::Store<O>, Er> {
        Ok(TensorShell {
            structure: f(self.structure)?,
        })
    }

    fn structure(&self) -> &S {
        &self.structure
    }

    fn map_same_structure(mut self, f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self {
        self.structure = f(self.structure);
        self
    }

    fn scalar(self) -> Option<Self::Scalar> {
        if self.structure.is_scalar() {
            Some(())
        } else {
            None
        }
    }

    fn scalar_ref(&self) -> Option<&Self::Scalar> {
        None
    }
    fn mut_structure(&mut self) -> &mut S {
        &mut self.structure
    }
}

impl<S: TensorStructure> HasName for TensorShell<S>
where
    S: HasName,
{
    type Args = S::Args;
    type Name = S::Name;

    fn args(&self) -> Option<Self::Args> {
        self.structure.args()
    }

    fn name(&self) -> Option<Self::Name> {
        self.structure.name()
    }

    fn set_name(&mut self, name: Self::Name) {
        self.structure.set_name(name);
    }
}

// impl<I> HasName for I
// where
//     I: HasStructure,
//     I::Structure: HasName,
// {
//     type Name = <I::Structure as HasName>::Name;
//     fn name(&self) -> Option<Cow<Self::Name>> {
//         self.structure().name()
//     }
//     fn set_name(&mut self, name: &Self::Name) {
//         self.mut_structure().set_name(name);
//     }
// }

impl<S> TensorShell<S> {
    pub fn new(structure: S) -> Self {
        Self { structure }
    }
}

impl<S: TensorStructure, O: From<S> + TensorStructure> CastStructure<TensorShell<O>>
    for TensorShell<S>
{
    fn cast_structure(self) -> TensorShell<O> {
        TensorShell {
            structure: self.structure.into(),
        }
    }
}

impl<S: TensorStructure> From<S> for TensorShell<S> {
    fn from(structure: S) -> Self {
        Self::new(structure)
    }
}

#[cfg(test)]
#[cfg(feature = "shadowing")]
mod shadowing_tests {}
