use crate::{
    algebra::{
        algebraic_traits::RefZero,
        upgrading_arithmetic::{TryFromUpgrade, TrySmallestUpgrade},
    },
    iterators::IteratableTensor,
    structure::{
        concrete_index::{ConcreteIndex, ExpandedIndex, FlatIndex},
        permuted::PermuteTensor,
        representation::RepName,
        slot::{AbsInd, Slot},
        CastStructure, HasName, HasStructure, IndexLess, OrderedStructure, PermutedStructure,
        ScalarStructure, ScalarTensor, TensorStructure, TracksCount,
    },
};
use crate::structure::dimension::Dimension;
use crate::structure::representation::Representation;
use crate::structure::slot::IsAbstractSlot;
use crate::structure::StructureError;
use delegate::delegate;

#[cfg(feature = "shadowing")]
use crate::{
    shadowing::{
        symbolica_utils::{atomic_expanded_label_id, IntoArgs, IntoSymbol},
        ShadowMapping, Shadowable,
    },
    structure::slot::ParseableAind,
    tensors::{data::DataIterator, parametric::TensorCoefficient},
};
use eyre::{eyre, Result};
use bincode::{Decode, Encode};
use indexmap::IndexMap;
use num::Zero;

use super::{CastData, GetTensorData, HasTensorData, SetTensorData, SparseTensor, StorageTensor};

use serde::{Deserialize, Serialize};
use std::{
    borrow::Cow,
    fmt::{Display, LowerExp},
    hash::Hash,
    ops::{Index, IndexMut},
};

#[cfg(feature = "shadowing")]
use symbolica::{atom::Atom, atom::Symbol};

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize, Hash, Eq, Encode, Decode)]
pub struct DenseTensor<T, S = OrderedStructure> {
    pub data: Vec<T>,
    pub structure: S,
}

impl<T, S> crate::network::Ref for DenseTensor<T, S> {
    type Ref<'a>
        = &'a DenseTensor<T, S>
    where
        Self: 'a;

    fn refer(&self) -> Self::Ref<'_> {
        self
    }
}

impl<T: RefZero, S: TensorStructure> DenseTensor<T, S> {
    pub fn ref_zero(&self) -> T {
        self.data.first().unwrap().ref_zero()
    }
}

impl<Aind: AbsInd, T: Clone, S: Clone + Into<IndexLess<R, Aind>>, R: RepName<Dual = R>>
    PermuteTensor for DenseTensor<T, S>
where
    S: TensorStructure<Slot = Slot<R, Aind>> + PermuteTensor<IdSlot = Slot<R, Aind>, Id = S>,
{
    type Id = DenseTensor<T, S>;
    type IdSlot = (T, Slot<R, Aind>);
    type Permuted = DenseTensor<T, S>;

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
        let (zero, i) = i;
        let (one, j) = j;
        let s = S::id(i, j);
        let mut data = vec![zero; s.size().unwrap()];
        for i in 0..usize::try_from(i.dim()).unwrap() {
            data[usize::from(s.flat_index([i, i]).unwrap())] = one.clone();
        }
        DenseTensor { data, structure: s }
    }

    fn permute_inds(self, permutation: &linnet::permutation::Permutation) -> Self::Permuted {
        let mut permuteds: IndexLess<R, Aind> = self.structure.clone().into();
        permutation.apply_slice_in_place(&mut permuteds.structure);

        let mut permuted = self.clone();
        for (i, d) in self.iter_expanded() {
            permuted
                .set_flat(
                    permuteds
                        .flat_index(i.apply_permutation(permutation))
                        .unwrap(),
                    d.clone(),
                )
                .unwrap();
        }
        permuted
    }

    fn permute_reps(self, _rep_perm: &linnet::permutation::Permutation) -> Self::Permuted {
        todo!()
    }
}

impl<T, S> TensorStructure for DenseTensor<T, S>
where
    S: TensorStructure,
{
    // type R = <T::Structure as TensorStructure>::R;
    type Indexed = DenseTensor<T, S::Indexed>;
    type Slot = S::Slot;

    fn reindex(
        self,
        indices: &[<Self::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::Indexed>, StructureError> {
        let res = self.structure.reindex(indices)?;

        Ok(PermutedStructure {
            structure: DenseTensor {
                structure: res.structure,
                data: self.data,
            },
            rep_permutation: res.rep_permutation,
            index_permutation: res.index_permutation,
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

impl<T, U: From<T> + Clone, I: TensorStructure + Clone> CastData<DenseTensor<U, I>>
    for DenseTensor<T, I>
{
    fn cast_data(self) -> DenseTensor<U, I> {
        DenseTensor {
            data: self.data.into_iter().map(|v| v.into()).collect(),
            structure: self.structure,
        }
    }
}

impl<T, S: TensorStructure, O: From<S> + TensorStructure> CastStructure<DenseTensor<T, O>>
    for DenseTensor<T, S>
{
    fn cast_structure(self) -> DenseTensor<T, O> {
        DenseTensor {
            data: self.data,
            structure: self.structure.into(),
        }
    }
}

#[cfg(feature = "shadowing")]
impl<T: Clone, S: TensorStructure> Shadowable for DenseTensor<T, S>
where
    S: HasName + Clone,
    S::Name: IntoSymbol,
    S::Args: IntoArgs,
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
}

#[cfg(feature = "shadowing")]
impl<T: Clone, S: TensorStructure, R> ShadowMapping<R> for DenseTensor<T, S>
where
    S: HasName + Clone,
    S::Name: IntoSymbol,
    S::Args: IntoArgs,
    R: From<T>,
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    // fn shadow_with_map<'a, U>(
    //     &self,
    //     fn_map: &mut symbolica::evaluate::FunctionMap<'a, R>,
    //     index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> U,
    // ) -> Option<ParamTensor<Self::Structure>>
    // where
    //     U: TensorCoefficient,
    //     R: From<T>,
    // {
    //     let mut data = vec![];
    //     for (i, d) in self.flat_iter() {
    //         let labeled_coef = index_to_atom(self.structure(), i).to_atom().unwrap();
    //         fn_map.add_constant(labeled_coef.clone().into(), d.clone().into());
    //         data.push(labeled_coef);
    //     }

    //     let param = DenseTensor {
    //         data,
    //         structure: self.structure.clone(),
    //     };

    //     Some(ParamTensor::Param(param.into()))
    // }

    fn append_map<U>(
        &self,
        fn_map: &mut symbolica::evaluate::FunctionMap<R>,
        index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> U,
    ) where
        U: TensorCoefficient,
    {
        for (i, d) in self.flat_iter() {
            let labeled_coef = index_to_atom(self.structure(), i).to_atom().unwrap();
            fn_map.add_constant(labeled_coef.clone(), d.clone().into());
        }
    }
}

impl<T: Display, I: TensorStructure> std::fmt::Display for DenseTensor<T, I> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        for (i, v) in self.iter_expanded() {
            s.push_str(&format!("{}: {}\n", i, v));
        }
        write!(f, "{}", s)
    }
}

impl<T: LowerExp, I: TensorStructure> std::fmt::LowerExp for DenseTensor<T, I> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        for (i, v) in self.iter_expanded() {
            s.push_str(&format!("{}: {:+e}\n", i, v));
        }
        write!(f, "{}", s)
    }
}

impl<T, I> Index<FlatIndex> for DenseTensor<T, I> {
    type Output = T;

    fn index(&self, index: FlatIndex) -> &Self::Output {
        let i: usize = index.into();
        &self.data[i]
    }
}

impl<T, I> IndexMut<FlatIndex> for DenseTensor<T, I> {
    fn index_mut(&mut self, index: FlatIndex) -> &mut Self::Output {
        let i: usize = index.into();
        &mut self.data[i]
    }
}

impl<T, I> HasStructure for DenseTensor<T, I>
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
        = DenseTensor<T, S>
    where
        S: TensorStructure;

    fn map_structure<O: TensorStructure>(
        self,
        f: impl FnOnce(Self::Structure) -> O,
    ) -> Self::Store<O> {
        DenseTensor {
            structure: f(self.structure),
            data: self.data,
        }
    }

    fn map_structure_result<O: TensorStructure, Er>(
        self,
        f: impl FnOnce(Self::Structure) -> Result<O, Er>,
    ) -> std::result::Result<Self::Store<O>, Er> {
        Ok(DenseTensor {
            structure: f(self.structure)?,
            data: self.data,
        })
    }

    fn structure(&self) -> &Self::Structure {
        &self.structure
    }

    fn map_same_structure(self, f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self {
        DenseTensor {
            data: self.data,
            structure: f(self.structure),
        }
    }
    fn mut_structure(&mut self) -> &mut Self::Structure {
        &mut self.structure
    }
    fn scalar(mut self) -> Option<Self::Scalar> {
        if self.is_scalar() {
            self.data.drain(0..).next()
        } else {
            None
        }
    }

    fn scalar_ref(&self) -> Option<Self::ScalarRef<'_>> {
        if self.is_scalar() {
            self.data.first()
        } else {
            None
        }
    }
}

impl<T, I> ScalarTensor for DenseTensor<T, I>
where
    I: TensorStructure + ScalarStructure,
{
    fn new_scalar(scalar: Self::Scalar) -> Self {
        DenseTensor {
            data: vec![scalar],
            structure: I::scalar_structure(),
        }
    }
}

impl<T, I> TracksCount for DenseTensor<T, I>
where
    I: TracksCount + TensorStructure,
{
    fn contractions_num(&self) -> usize {
        self.structure.contractions_num()
    }
}

impl<T, S> HasName for DenseTensor<T, S>
where
    S: HasName + TensorStructure,
{
    type Args = S::Args;
    type Name = S::Name;
    fn name(&self) -> Option<Self::Name> {
        self.structure.name()
    }

    fn args(&self) -> Option<Self::Args> {
        self.structure.args()
    }

    fn set_name(&mut self, name: Self::Name) {
        self.structure.set_name(name);
    }
}

impl<T, I> DenseTensor<T, I>
where
    I: TensorStructure,
{
    pub fn default(structure: I) -> Self
    where
        T: Default + Clone,
    {
        Self::fill(structure, T::default())
    }

    pub fn fill(structure: I, fill: T) -> Self
    where
        T: Clone,
    {
        let length = if structure.is_scalar() {
            1
        } else {
            structure.size().unwrap()
        };
        DenseTensor {
            data: vec![fill; length],
            structure,
        }
    }
}

impl<T, S: TensorStructure> DenseTensor<T, S> {
    pub fn repeat(structure: S, r: T) -> Self
    where
        T: Clone,
    {
        let length = if structure.is_scalar() {
            1
        } else {
            structure.size().unwrap()
        };
        DenseTensor {
            data: vec![r; length],
            structure,
        }
    }
}

impl<T: Zero + Clone, I> DenseTensor<T, I>
where
    I: TensorStructure,
{
    pub fn zero(structure: I) -> Self {
        let length = if structure.is_scalar() {
            1
        } else {
            structure.size().unwrap()
        };
        DenseTensor {
            data: vec![T::zero(); length],
            structure,
        }
    }
}

// impl<T,S:TensorStructure> DenseTensor<T,S>{}

impl<U, I> DenseTensor<U, I>
where
    I: TensorStructure + Clone,
{
    pub fn try_upgrade<T>(&self) -> Option<Cow<'_, DenseTensor<U::LCM, I>>>
    where
        U: TrySmallestUpgrade<T>,
        U::LCM: Clone,
    {
        let structure = self.structure.clone();
        let data: Option<Vec<U::LCM>> = self
            .data
            .iter()
            .map(|v| match v.try_upgrade() {
                Some(Cow::Owned(u)) => Some(u),
                Some(Cow::Borrowed(u)) => Some(u.clone()),
                None => None,
            })
            .collect();
        Some(Cow::Owned(DenseTensor {
            data: data?,
            structure,
        }))
    }
}

impl<T: Clone, I> DenseTensor<T, I>
where
    I: TensorStructure,
{
    /// Generates a new dense tensor from the given data and structure
    pub fn from_data(data: Vec<T>, structure: I) -> Result<Self> {
        if data.len() != structure.size()? && !(data.len() == 1 && structure.is_scalar()) {
            return Err(eyre!("Data length does not match shape"));
        }
        Ok(DenseTensor { data, structure })
    }

    pub fn cast<U>(&self) -> DenseTensor<U, I>
    where
        U: Clone + From<T>,
        I: Clone,
    {
        let data = self.data.iter().map(|x| x.clone().into()).collect();
        DenseTensor {
            data,
            structure: self.structure.clone(),
        }
    }

    /// Generates a new dense tensor from the given data and structure, truncating the data if it is too long with respect to the structure
    pub fn from_data_coerced(data: &[T], structure: I) -> Result<Self> {
        if data.len() < structure.size()? {
            return Err(eyre!("Data length is too small"));
        }
        let mut data = data.to_vec();
        if structure.is_scalar() {
            data.truncate(1);
        } else {
            data.truncate(structure.size()?);
        }
        Ok(DenseTensor { data, structure })
    }
}

impl<T, U, S> TryFromUpgrade<DenseTensor<U, S>> for DenseTensor<T, S>
where
    U: TrySmallestUpgrade<T, LCM = T>,
    S: TensorStructure + Clone,
    T: Clone,
{
    fn try_from_upgrade(data: &DenseTensor<U, S>) -> Option<DenseTensor<T, S>> {
        data.try_upgrade().map(Cow::into_owned)
    }
}

impl<T, I> DenseTensor<T, I>
where
    I: TensorStructure + Clone,
{
    /// converts the dense tensor to a sparse tensor, with the same structure
    pub fn to_sparse(&self) -> SparseTensor<T, I>
    where
        T: Clone + Default + PartialEq,
    {
        let mut sparse = SparseTensor::empty(self.structure.clone(), T::default());
        for (i, value) in self.iter_expanded() {
            if *value != T::default() {
                let _ = sparse.set(&i, value.clone());
            }
        }
        sparse
    }
}
// why no specialization? :(
// impl<T, U> From<DenseTensor<U>> for DenseTensor<T>
// where
//     U: Into<T>,
// {
//     fn from(other: DenseTensor<U>) -> Self {
//         let data = other.data.into_iter().map(|x| x.into()).collect();
//         DenseTensor {
//             data,
//             structure: other.structure,
//         }
//     }
// }

impl<T, I> DenseTensor<T, I>
where
    I: Clone + TensorStructure,
{
    pub fn convert_to<U>(&self) -> DenseTensor<U, I>
    where
        U: for<'a> From<&'a T>,
    {
        let data = self.data.iter().map(|x| x.into()).collect();
        DenseTensor {
            data,
            structure: self.structure.clone(),
        }
    }
}

impl<T, I> HasTensorData for DenseTensor<T, I>
where
    T: Clone,
    I: TensorStructure + Clone,
{
    type Data = T;
    fn data(&self) -> Vec<T> {
        self.data.clone()
    }

    fn indices(&self) -> Vec<ExpandedIndex> {
        let mut indices = Vec::new();
        for i in 0..self.size().unwrap() {
            indices.push(self.expanded_index(i.into()).unwrap());
        }
        indices
    }

    fn hashmap(&self) -> IndexMap<ExpandedIndex, T> {
        let mut hashmap = IndexMap::new();
        for (k, v) in self.iter_expanded() {
            hashmap.insert(k.clone(), v.clone());
        }
        hashmap
    }
    #[cfg(feature = "shadowing")]
    fn symhashmap(&self, name: Symbol, args: &[Atom]) -> std::collections::HashMap<Atom, T> {
        let mut hashmap = std::collections::HashMap::new();

        for (k, v) in self.iter_expanded() {
            hashmap.insert(atomic_expanded_label_id(&k, name, args), v.clone());
        }
        hashmap
    }
}

impl<T, I> SetTensorData for DenseTensor<T, I>
where
    I: TensorStructure,
{
    type SetData = T;
    fn set(&mut self, indices: &[ConcreteIndex], value: T) -> Result<()> {
        self.verify_indices(indices)?;
        let idx = self.flat_index(indices);
        if let Ok(i) = idx {
            self[i] = value;
        }
        Ok(())
    }

    fn set_flat(&mut self, index: FlatIndex, value: T) -> Result<()> {
        if index < self.size()?.into() {
            self[index] = value;
        } else {
            return Err(eyre!("Index out of bounds"));
        }
        Ok(())
    }
}

impl<T, I> GetTensorData for DenseTensor<T, I>
where
    I: TensorStructure,
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
    fn get_ref_linear(&self, index: FlatIndex) -> Option<&T> {
        let i: usize = index.into();
        self.data.get(i)
    }

    fn get_ref<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<&T> {
        if let Ok(idx) = self.flat_index(&indices) {
            Ok(&self[idx])
        } else if self.structure.is_scalar() && indices.as_ref().is_empty() {
            Ok(&self.data[0])
        } else {
            Err(eyre!("Index out of bounds"))
        }
    }

    fn get_mut_linear(&mut self, index: FlatIndex) -> Option<&mut T> {
        let i: usize = index.into();
        self.data.get_mut(i)
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

impl<S: TensorStructure + Clone, D> StorageTensor for DenseTensor<D, S> {
    type ContainerData<Data> = DenseTensor<Data, S>;

    type Data = D;

    fn map_data_self(self, f: impl Fn(Self::Data) -> Self::Data) -> Self {
        self.map_data(f)
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

    fn map_data_ref<U>(&self, f: impl Fn(&D) -> U) -> DenseTensor<U, S> {
        let data = self.data.iter().map(f).collect();
        DenseTensor {
            data,
            structure: self.structure.clone(),
        }
    }

    fn map_data_ref_result<U, E>(
        &self,
        f: impl Fn(&D) -> Result<U, E>,
    ) -> Result<DenseTensor<U, S>, E> {
        let data: Result<Vec<U>, E> = self.data.iter().map(f).collect();
        Ok(DenseTensor {
            data: data?,
            structure: self.structure.clone(),
        })
    }

    fn map_data_ref_mut<U>(&mut self, f: impl FnMut(&mut D) -> U) -> DenseTensor<U, S> {
        let data = self.data.iter_mut().map(f).collect();
        DenseTensor {
            data,
            structure: self.structure.clone(),
        }
    }

    fn map_data_ref_mut_result<U, E>(
        &mut self,
        f: impl FnMut(&mut Self::Data) -> Result<U, E>,
    ) -> Result<Self::ContainerData<U>, E> {
        let data: Result<Vec<U>, E> = self.data.iter_mut().map(f).collect();
        Ok(DenseTensor {
            data: data?,
            structure: self.structure.clone(),
        })
    }

    fn map_data_mut(&mut self, f: impl FnMut(&mut D)) {
        self.data.iter_mut().for_each(f);
    }

    fn map_data<U>(self, f: impl Fn(D) -> U) -> DenseTensor<U, S> {
        let data = self.data.into_iter().map(f).collect();
        DenseTensor {
            data,
            structure: self.structure,
        }
    }
}
