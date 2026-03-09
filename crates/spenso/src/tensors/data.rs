use crate::{
    algebra::{
        algebraic_traits::RefZero,
        upgrading_arithmetic::{TryFromUpgrade, TrySmallestUpgrade},
    },
    iterators::{DenseTensorLinearIterator, SparseTensorLinearIterator},
    structure::{
        CastStructure, HasName, HasStructure, IndexLess, OrderedStructure, PermutedStructure,
        ScalarStructure, ScalarTensor, StructureContract, TensorStructure, TracksCount,
        concrete_index::{ConcreteIndex, ExpandedIndex, FlatIndex},
        permuted::PermuteTensor,
        representation::RepName,
        slot::{AbsInd, Slot},
    },
};

use crate::structure::StructureError;
use crate::structure::dimension::Dimension;
use crate::structure::representation::Representation;
use crate::structure::slot::IsAbstractSlot;
use delegate::delegate;

#[cfg(feature = "shadowing")]
use crate::{
    shadowing::{
        ShadowMapping, Shadowable,
        symbolica_utils::{IntoArgs, IntoSymbol},
    },
    structure::slot::ParseableAind,
    tensors::parametric::TensorCoefficient,
};

use eyre::{Result, eyre};

use bincode::{Decode, Encode};
use derive_more::From;
use enum_try_as_inner::EnumTryAsInner;
use indexmap::IndexMap;
use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use std::{borrow::Cow, fmt::Display, hash::Hash};

#[cfg(feature = "shadowing")]
use symbolica::{atom::Atom, atom::Symbol};
pub trait DataIterator<T> {
    type FlatIter<'a>: Iterator<Item = (FlatIndex, &'a T)>
    where
        Self: 'a,
        T: 'a;

    fn flat_iter(&self) -> Self::FlatIter<'_>;
}

impl<T, I> DataIterator<T> for SparseTensor<T, I> {
    type FlatIter<'a>
        = SparseTensorLinearIterator<'a, T>
    where
        I: 'a,
        T: 'a;

    fn flat_iter(&self) -> Self::FlatIter<'_> {
        SparseTensorLinearIterator::new(self)
    }
}

impl<T, I: TensorStructure> DataIterator<T> for DenseTensor<T, I> {
    type FlatIter<'a>
        = DenseTensorLinearIterator<'a, T, I>
    where
        I: 'a,
        T: 'a;

    fn flat_iter(&self) -> Self::FlatIter<'_> {
        DenseTensorLinearIterator::new(self)
    }
}

#[allow(dead_code)]
trait Settable {
    type SetData;
    fn set(&mut self, index: usize, data: Self::SetData);
}

impl<T> Settable for Vec<T> {
    type SetData = T;
    fn set(&mut self, index: usize, data: T) {
        self[index] = data;
    }
}

impl<T> Settable for HashMap<usize, T> {
    type SetData = T;
    fn set(&mut self, index: usize, data: T) {
        self.insert(index, data);
    }
}

/// Trait for getting the data of a tensor
pub trait HasTensorData: HasStructure {
    type Data: Clone;
    // type Storage: Settable<SetData = Self::Data>;
    /// Returns all the data in the tensor, withouth any structure
    fn data(&self) -> Vec<Self::Data>;

    /// Returns all the indices of the tensor, the order of the indices is the same as the order of the data
    fn indices(&self) -> Vec<ExpandedIndex>;

    /// Returns a hashmap of the data, with the (expanded) indices as keys
    fn hashmap(&self) -> IndexMap<ExpandedIndex, Self::Data>;

    /// Returns a hashmap of the data, with the the shadowed indices as keys
    #[cfg(feature = "shadowing")]
    fn symhashmap(&self, name: Symbol, args: &[Atom]) -> HashMap<Atom, Self::Data>;

    // fn map(&self, f: impl Fn(&Self::Data) -> Self::Data) -> Self;
}

/// Trait for setting the data of a tensor
pub trait SetTensorData {
    type SetData;
    /// Set the data at the given indices, returns an error if the indices are out of bounds
    ///
    /// # Errors
    ///
    /// Forwards the error from [`TensorStructure::verify_indices`]
    ///
    fn set(&mut self, indices: &[ConcreteIndex], value: Self::SetData) -> Result<()>;

    fn set_flat(&mut self, index: FlatIndex, value: Self::SetData) -> Result<()>;
}

/// Trait for getting the data of a tensor
pub trait GetTensorData {
    type GetDataRef<'a>
    where
        Self: 'a;

    type GetDataRefMut<'a>
    where
        Self: 'a;

    type GetDataOwned;

    fn get_ref<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<Self::GetDataRef<'_>>;

    fn get_ref_linear(&self, index: FlatIndex) -> Option<Self::GetDataRef<'_>>;

    fn get_mut_linear(&mut self, index: FlatIndex) -> Option<Self::GetDataRefMut<'_>>;

    fn get_owned<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<Self::GetDataOwned>
    where
        Self::GetDataOwned: Clone;

    fn get_owned_linear(&self, index: FlatIndex) -> Option<Self::GetDataOwned>
    where
        Self::GetDataOwned: Clone;
}

pub mod dense;
pub use dense::DenseTensor;
pub mod sparse;
pub use sparse::SparseTensor;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct Tensor<Store, Structure> {
    pub elements: Store,
    pub structure: Structure,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct SparseStore<T> {
    pub elements: HashMap<FlatIndex, T>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct DenseStore<T> {
    pub elements: Vec<T>,
}

pub trait CastData<O: HasTensorData>: HasStructure<Structure = O::Structure> {
    fn cast_data(self) -> O;
}

/// Enum for storing either a dense or a sparse tensor, with the same structure
#[derive(
    Debug, Clone, EnumTryAsInner, Serialize, Deserialize, From, Hash, PartialEq, Eq, Encode, Decode,
)]
#[derive_err(Debug)]
pub enum DataTensor<T, I = OrderedStructure> {
    Dense(DenseTensor<T, I>),
    Sparse(SparseTensor<T, I>),
}

impl<T, S> crate::network::Ref for DataTensor<T, S> {
    type Ref<'a>
        = &'a DataTensor<T, S>
    where
        Self: 'a;

    fn refer(&self) -> Self::Ref<'_> {
        self
    }
}

impl<T: RefZero, S: TensorStructure> DataTensor<T, S> {
    pub fn ref_zero(&self) -> T {
        match self {
            DataTensor::Dense(d) => d.ref_zero(),
            DataTensor::Sparse(s) => s.ref_zero(),
        }
    }
}

pub trait SparseOrDense {
    fn to_sparse(self) -> Self;

    fn to_sparse_mut(&mut self);
    fn to_dense_mut(&mut self);
    fn to_dense(self) -> Self;
}

impl<T: Clone, U: From<T> + Clone, I: TensorStructure + Clone> CastData<DataTensor<U, I>>
    for DataTensor<T, I>
{
    fn cast_data(self) -> DataTensor<U, I> {
        match self {
            Self::Dense(d) => DataTensor::Dense(d.cast_data()),

            Self::Sparse(d) => DataTensor::Sparse(d.cast_data()),
        }
    }
}

impl<T: Clone, S: TensorStructure, O: From<S> + TensorStructure> CastStructure<DataTensor<T, O>>
    for DataTensor<T, S>
{
    fn cast_structure(self) -> DataTensor<T, O> {
        match self {
            DataTensor::Dense(d) => DataTensor::Dense(d.cast_structure()),
            DataTensor::Sparse(d) => DataTensor::Sparse(d.cast_structure()),
        }
    }
}

#[cfg(feature = "shadowing")]
impl<T: Clone, S: TensorStructure> Shadowable for DataTensor<T, S>
where
    S: HasName + Clone,
    S::Name: IntoSymbol,
    S::Args: IntoArgs,
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
}

#[cfg(feature = "shadowing")]
impl<T: Clone, S: TensorStructure, R> ShadowMapping<R> for DataTensor<T, S>
where
    S: HasName + Clone,
    S::Name: IntoSymbol,
    S::Args: IntoArgs,
    R: From<T>,
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    // fn shadow_with_map<'a, U>(
    //     &'a self,
    //     fn_map: &mut symbolica::evaluate::FunctionMap<'a, R>,
    //     index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> U,
    // ) -> Option<ParamTensor<Self::Structure>>
    // where
    //     U: TensorCoefficient,
    // {
    //     match self {
    //         DataTensor::Dense(d) => d.shadow_with_map(fn_map, index_to_atom),
    //         DataTensor::Sparse(s) => s.shadow_with_map(fn_map, index_to_atom),
    //     }
    // }

    fn append_map<U>(
        &self,
        fn_map: &mut symbolica::evaluate::FunctionMap<R>,
        index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> U,
    ) where
        U: TensorCoefficient,
    {
        match self {
            DataTensor::Dense(d) => d.append_map(fn_map, index_to_atom),
            DataTensor::Sparse(s) => s.append_map(fn_map, index_to_atom),
        }
    }
}

impl<T, I> DataTensor<T, I>
where
    I: TensorStructure + Clone,
{
    pub fn actual_size(&self) -> usize {
        match self {
            DataTensor::Dense(d) => d.data.len(),
            DataTensor::Sparse(s) => s.elements.len(),
        }
    }
    pub fn to_bare_sparse(self) -> SparseTensor<T, I>
    where
        T: Clone + Default + PartialEq,
    {
        match self {
            DataTensor::Dense(d) => d.to_sparse(),
            DataTensor::Sparse(s) => s,
        }
    }

    pub fn to_bare_dense(self) -> DenseTensor<T, I>
    where
        T: Clone + PartialEq,
    {
        match self {
            DataTensor::Dense(d) => d,
            DataTensor::Sparse(s) => s.to_dense(),
        }
    }
}

impl<T, I> SparseOrDense for DataTensor<T, I>
where
    I: TensorStructure + Clone,
    T: Clone + PartialEq + Default,
{
    fn to_dense(self) -> Self {
        DataTensor::Dense(self.to_bare_dense())
    }

    fn to_dense_mut(&mut self) {
        let s = std::mem::replace(
            self,
            DataTensor::Dense(DenseTensor {
                structure: self.structure().clone(),
                data: vec![],
            }),
        );

        *self = s.to_dense();
    }

    fn to_sparse_mut(&mut self) {
        let s = std::mem::replace(
            self,
            DataTensor::Dense(DenseTensor {
                structure: self.structure().clone(),
                data: vec![],
            }),
        );

        *self = s.to_sparse();
    }

    fn to_sparse(self) -> Self {
        DataTensor::Sparse(self.to_bare_sparse())
    }
}

impl<T, I> HasTensorData for DataTensor<T, I>
where
    I: TensorStructure + Clone,
    T: Clone,
{
    type Data = T;
    fn data(&self) -> Vec<T> {
        match self {
            DataTensor::Dense(d) => d.data(),
            DataTensor::Sparse(s) => s.data(),
        }
    }

    fn indices(&self) -> Vec<ExpandedIndex> {
        match self {
            DataTensor::Dense(d) => d.indices(),
            DataTensor::Sparse(s) => s.indices(),
        }
    }

    fn hashmap(&self) -> IndexMap<ExpandedIndex, T> {
        match self {
            DataTensor::Dense(d) => d.hashmap(),
            DataTensor::Sparse(s) => s.hashmap(),
        }
    }
    #[cfg(feature = "shadowing")]
    fn symhashmap(&self, name: Symbol, args: &[Atom]) -> HashMap<Atom, T> {
        match self {
            DataTensor::Dense(d) => d.symhashmap(name, args),
            DataTensor::Sparse(s) => s.symhashmap(name, args),
        }
    }
}

impl<T: Clone, Aind: AbsInd, S: Clone + Into<IndexLess<R, Aind>>, R: RepName<Dual = R>>
    PermuteTensor for DataTensor<T, S>
where
    S: TensorStructure<Slot = Slot<R, Aind>> + PermuteTensor<IdSlot = Slot<R, Aind>, Id = S>,
{
    type Id = DataTensor<T, S>;
    type IdSlot = (T, Slot<R, Aind>);
    type Permuted = DataTensor<T, S>;

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
        DataTensor::Sparse(SparseTensor::id(i, j))
    }

    fn permute_inds(self, permutation: &linnet::permutation::Permutation) -> Self::Permuted {
        match self {
            DataTensor::Dense(d) => DataTensor::Dense(d.permute_inds(permutation)),
            DataTensor::Sparse(s) => DataTensor::Sparse(s.permute_inds(permutation)),
        }
    }

    fn permute_reps(self, rep_perm: &linnet::permutation::Permutation) -> Self::Permuted {
        match self {
            DataTensor::Dense(d) => DataTensor::Dense(d.permute_reps(rep_perm)),
            DataTensor::Sparse(s) => DataTensor::Sparse(s.permute_reps(rep_perm)),
        }
    }
}

impl<T, S> TensorStructure for DataTensor<T, S>
where
    S: TensorStructure,
{
    // type R = <T::Structure as TensorStructure>::R;
    type Indexed = DataTensor<T, S::Indexed>;
    type Slot = S::Slot;

    fn reindex(
        self,
        indices: &[<Self::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::Indexed>, StructureError> {
        Ok(match self {
            DataTensor::Dense(d) => {
                let res = d.reindex(indices)?;
                PermutedStructure {
                    structure: DataTensor::Dense(res.structure),
                    index_permutation: res.index_permutation,
                    rep_permutation: res.rep_permutation,
                }
            }
            DataTensor::Sparse(d) => {
                let res = d.reindex(indices)?;
                PermutedStructure {
                    structure: DataTensor::Sparse(res.structure),
                    index_permutation: res.index_permutation,
                    rep_permutation: res.rep_permutation,
                }
            }
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

impl<T, I> HasStructure for DataTensor<T, I>
where
    I: TensorStructure,
{
    type Scalar = T;
    type ScalarRef<'a>
        = &'a T
    where
        Self: 'a;
    type Structure = I;

    type Store<S>
        = DataTensor<T, S>
    where
        S: TensorStructure;

    fn map_structure<O: TensorStructure>(self, f: impl Fn(Self::Structure) -> O) -> Self::Store<O> {
        match self {
            DataTensor::Dense(d) => DataTensor::Dense(d.map_structure(f)),
            DataTensor::Sparse(s) => DataTensor::Sparse(s.map_structure(f)),
        }
    }

    fn map_structure_result<O: TensorStructure, Er>(
        self,
        f: impl Fn(Self::Structure) -> Result<O, Er>,
    ) -> std::result::Result<Self::Store<O>, Er> {
        Ok(match self {
            DataTensor::Dense(d) => match d.map_structure_result(f) {
                Ok(d) => DataTensor::Dense(d),
                Err(e) => return Err(e),
            },
            DataTensor::Sparse(s) => match s.map_structure_result(f) {
                Ok(s) => DataTensor::Sparse(s),
                Err(e) => return Err(e),
            },
        })
    }

    fn structure(&self) -> &Self::Structure {
        match self {
            DataTensor::Dense(d) => d.structure(),
            DataTensor::Sparse(s) => s.structure(),
        }
    }
    fn mut_structure(&mut self) -> &mut Self::Structure {
        match self {
            DataTensor::Dense(d) => d.mut_structure(),
            DataTensor::Sparse(s) => s.mut_structure(),
        }
    }

    fn scalar(self) -> Option<Self::Scalar> {
        match self {
            DataTensor::Dense(d) => d.scalar(),
            DataTensor::Sparse(s) => s.scalar(),
        }
    }

    fn scalar_ref(&self) -> Option<Self::ScalarRef<'_>> {
        match self {
            DataTensor::Dense(d) => d.scalar_ref(),
            DataTensor::Sparse(s) => s.scalar_ref(),
        }
    }

    fn map_same_structure(self, f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self {
        match self {
            DataTensor::Dense(d) => DataTensor::Dense(d.map_same_structure(f)),
            DataTensor::Sparse(s) => DataTensor::Sparse(s.map_same_structure(f)),
        }
    }
}

impl<T, I> ScalarTensor for DataTensor<T, I>
where
    I: TensorStructure + ScalarStructure,
{
    fn new_scalar(scalar: Self::Scalar) -> Self {
        DataTensor::Dense(DenseTensor::new_scalar(scalar))
    }
}

impl<T: Display, S: TensorStructure> Display for DataTensor<T, S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DataTensor::Dense(d) => write!(f, "{}", d),
            DataTensor::Sparse(s) => write!(f, "{}", s),
        }
    }
}

impl<T, I> TracksCount for DataTensor<T, I>
where
    I: TracksCount,
    T: Clone,
    I: TensorStructure,
{
    fn contractions_num(&self) -> usize {
        match self {
            DataTensor::Dense(d) => d.contractions_num(),
            DataTensor::Sparse(s) => s.contractions_num(),
        }
    }
}

impl<T, S> HasName for DataTensor<T, S>
where
    S: HasName + TensorStructure,
{
    type Args = S::Args;
    type Name = S::Name;
    fn name(&self) -> Option<Self::Name> {
        match self {
            DataTensor::Dense(d) => d.name(),
            DataTensor::Sparse(s) => s.name(),
        }
    }

    fn args(&self) -> Option<Self::Args> {
        match self {
            DataTensor::Dense(d) => d.args(),
            DataTensor::Sparse(s) => s.args(),
        }
    }

    fn set_name(&mut self, name: Self::Name) {
        match self {
            DataTensor::Dense(d) => d.set_name(name),
            DataTensor::Sparse(s) => s.set_name(name),
        }
    }
}
impl<U, I> DataTensor<U, I>
where
    I: TensorStructure + Clone,
{
    pub fn try_upgrade<T>(&self) -> Option<Cow<'_, DataTensor<U::LCM, I>>>
    where
        U: TrySmallestUpgrade<T>,
        U::LCM: Clone,
    {
        match self {
            DataTensor::Dense(d) => d
                .try_upgrade()
                .map(|x| Cow::Owned(DataTensor::Dense(x.into_owned()))),
            DataTensor::Sparse(s) => s
                .try_upgrade()
                .map(|x| Cow::Owned(DataTensor::Sparse(x.into_owned()))),
        }
    }
}

impl<T, U, S> TryFromUpgrade<DataTensor<U, S>> for DataTensor<T, S>
where
    U: TrySmallestUpgrade<T, LCM = T>,
    S: TensorStructure + Clone,
    T: Clone,
{
    fn try_from_upgrade(data: &DataTensor<U, S>) -> Option<DataTensor<T, S>> {
        data.try_upgrade().map(Cow::into_owned)
    }
}

impl<T, S> SetTensorData for DataTensor<T, S>
where
    S: TensorStructure,
{
    type SetData = T;

    fn set(&mut self, indices: &[ConcreteIndex], value: Self::SetData) -> Result<()> {
        match self {
            DataTensor::Dense(d) => d.set(indices, value),
            DataTensor::Sparse(s) => s.set(indices, value),
        }
    }

    fn set_flat(&mut self, index: FlatIndex, value: Self::SetData) -> Result<()> {
        match self {
            DataTensor::Dense(d) => d.set_flat(index, value),
            DataTensor::Sparse(s) => s.set_flat(index, value),
        }
    }
}

impl<T, S> GetTensorData for DataTensor<T, S>
where
    S: TensorStructure,
{
    type GetDataRef<'a>
        = &'a T
    where
        Self: 'a;
    type GetDataRefMut<'a>
        = &'a mut T
    where
        Self: 'a;

    type GetDataOwned = T;
    fn get_ref<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<&T> {
        match self {
            DataTensor::Dense(d) => d.get_ref(indices),
            DataTensor::Sparse(s) => s.get_ref(indices),
        }
    }

    fn get_ref_linear(&self, index: FlatIndex) -> Option<&T> {
        match self {
            DataTensor::Dense(d) => d.get_ref_linear(index),
            DataTensor::Sparse(s) => s.get_ref_linear(index),
        }
    }

    fn get_mut_linear(&mut self, index: FlatIndex) -> Option<&mut T> {
        match self {
            DataTensor::Dense(d) => d.get_mut_linear(index),
            DataTensor::Sparse(s) => s.get_mut_linear(index),
        }
    }

    fn get_owned<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<Self::GetDataOwned>
    where
        Self::GetDataOwned: Clone,
    {
        self.get_ref(indices).cloned()
    }

    fn get_owned_linear(&self, index: FlatIndex) -> Option<Self::GetDataOwned>
    where
        Self::GetDataOwned: Clone,
    {
        self.get_ref_linear(index).cloned()
    }
}

impl<T, S: TensorStructure + StructureContract + Clone> DenseTensor<DenseTensor<T, S>, S> {
    pub fn flatten(self) -> Result<DenseTensor<T, S>> {
        // Check that the data is not empty
        if self.data.is_empty() {
            return Err(eyre!("Cannot flatten an empty tensor"));
        }

        // Verify that all inner tensors have the same structure
        let first_inner_structure = &self.data[0].structure;
        // for tensor in &self.data {
        //     if tensor.structure != *first_inner_structure {
        //         return Err(eyre!("Inner tensors have different structures"));
        //     }
        // }

        // Concatenate the outer and inner structures
        let combined_structure = self.structure.merge(first_inner_structure)?.0;

        // Flatten the data by concatenating inner tensors' data
        let data = self
            .data
            .into_iter()
            .flat_map(|tensor| tensor.data.into_iter())
            .collect();

        // Create the new flattened tensor
        Ok(DenseTensor {
            data,
            structure: combined_structure,
        })
    }
}

impl<T: Clone, S: TensorStructure + StructureContract + Clone> DataTensor<DataTensor<T, S>, S> {
    pub fn flatten(self, fill: &T) -> Result<DataTensor<T, S>> {
        let densified = self.map_data(|a| match a {
            DataTensor::Dense(d) => d,
            DataTensor::Sparse(s) => s.to_dense_with(fill),
        });
        match densified {
            DataTensor::Dense(d) => d.flatten().map(DataTensor::Dense),
            DataTensor::Sparse(s) => {
                let dense_fill = DenseTensor::fill(s.structure().clone(), fill.clone());
                s.to_dense_with(&dense_fill)
                    .flatten()
                    .map(DataTensor::Dense)
            }
        }
    }
}

pub trait StorageTensor: Sized + HasStructure<Structure: Clone> {
    // type ContainerStructure<S: TensorStructure>: HasStructure<Structure = S>;
    type ContainerData<Data>: HasStructure<Structure = Self::Structure>;
    type Data;

    // fn map_structure<S>(self, f: impl Fn(Self::Structure) -> S) -> Self::ContainerStructure<S>
    // where
    //     S: TensorStructure;

    fn map_data_ref<U>(&self, f: impl Fn(&Self::Data) -> U) -> Self::ContainerData<U>;

    fn map_data_ref_self(&self, f: impl Fn(&Self::Data) -> Self::Data) -> Self;

    fn map_data_ref_result<U, E>(
        &self,
        f: impl Fn(&Self::Data) -> Result<U, E>,
    ) -> Result<Self::ContainerData<U>, E>;

    fn map_data_ref_result_self<E>(
        &self,
        f: impl Fn(&Self::Data) -> Result<Self::Data, E>,
    ) -> Result<Self, E>;

    fn map_data_ref_mut<U>(
        &mut self,
        f: impl FnMut(&mut Self::Data) -> U,
    ) -> Self::ContainerData<U>;

    fn map_data_ref_mut_result<U, E>(
        &mut self,
        f: impl FnMut(&mut Self::Data) -> Result<U, E>,
    ) -> Result<Self::ContainerData<U>, E>;

    fn map_data_ref_mut_self(&mut self, f: impl FnMut(&mut Self::Data) -> Self::Data) -> Self;

    fn map_data_mut(&mut self, f: impl FnMut(&mut Self::Data));

    fn map_data<U>(self, f: impl Fn(Self::Data) -> U) -> Self::ContainerData<U>;

    fn map_data_self(self, f: impl Fn(Self::Data) -> Self::Data) -> Self;
}

impl<S: TensorStructure + Clone, T> StorageTensor for DataTensor<T, S> {
    type Data = T;
    type ContainerData<Data> = DataTensor<Data, S>;

    fn map_data_self(self, f: impl Fn(Self::Data) -> Self::Data) -> Self {
        self.map_data(f)
    }

    fn map_data_ref_mut_result<U, E>(
        &mut self,
        f: impl FnMut(&mut Self::Data) -> Result<U, E>,
    ) -> Result<Self::ContainerData<U>, E> {
        match self {
            DataTensor::Dense(d) => Ok(DataTensor::Dense(d.map_data_ref_mut_result(f)?)),
            DataTensor::Sparse(s) => Ok(DataTensor::Sparse(s.map_data_ref_mut_result(f)?)),
        }
    }

    fn map_data_ref_self(&self, f: impl Fn(&Self::Data) -> Self::Data) -> Self {
        self.map_data_ref(f)
    }

    fn map_data_ref_mut_self(&mut self, f: impl FnMut(&mut Self::Data) -> Self::Data) -> Self {
        self.map_data_ref_mut(f)
    }

    fn map_data_ref_result_self<E>(
        &self,
        f: impl Fn(&Self::Data) -> Result<Self::Data, E>,
    ) -> Result<Self, E> {
        self.map_data_ref_result(f)
    }

    // fn map_structure<S2: TensorStructure>(self, f: impl Fn(S) -> S2) -> DataTensor<T, S2> {
    //     match self {
    //         DataTensor::Dense(d) => DataTensor::Dense(d.map_structure(f)),
    //         DataTensor::Sparse(s) => DataTensor::Sparse(s.map_structure(f)),
    //     }
    // }

    fn map_data_ref_result<U, E>(
        &self,
        f: impl Fn(&T) -> Result<U, E>,
    ) -> Result<DataTensor<U, S>, E> {
        match self {
            DataTensor::Dense(d) => Ok(DataTensor::Dense(d.map_data_ref_result(f)?)),
            DataTensor::Sparse(s) => Ok(DataTensor::Sparse(s.map_data_ref_result(f)?)),
        }
    }

    fn map_data_ref<U>(&self, f: impl Fn(&T) -> U) -> DataTensor<U, S> {
        match self {
            DataTensor::Dense(d) => DataTensor::Dense(d.map_data_ref(f)),
            DataTensor::Sparse(s) => DataTensor::Sparse(s.map_data_ref(f)),
        }
    }

    fn map_data<U>(self, f: impl Fn(T) -> U) -> DataTensor<U, S> {
        match self {
            DataTensor::Dense(d) => DataTensor::Dense(d.map_data(f)),
            DataTensor::Sparse(s) => DataTensor::Sparse(s.map_data(f)),
        }
    }

    fn map_data_mut(&mut self, f: impl FnMut(&mut T)) {
        match self {
            DataTensor::Dense(d) => d.map_data_mut(f),
            DataTensor::Sparse(s) => s.map_data_mut(f),
        }
    }

    fn map_data_ref_mut<U>(&mut self, f: impl FnMut(&mut T) -> U) -> DataTensor<U, S> {
        match self {
            DataTensor::Dense(d) => DataTensor::Dense(d.map_data_ref_mut(f)),
            DataTensor::Sparse(s) => DataTensor::Sparse(s.map_data_ref_mut(f)),
        }
    }
}
