use crate::algebra::algebraic_traits::RefZero;
use crate::structure::dimension::Dimension;
use crate::structure::permuted::PermuteTensor;
use crate::structure::representation::{RepName, Representation};
#[cfg(feature = "shadowing")]
use crate::structure::slot::ParseableAind;
use crate::structure::slot::{AbsInd, IsAbstractSlot, Slot};
use crate::structure::{IndexLess, PermutedStructure, StructureError};
use crate::{
    algebra::algebraic_traits::IsZero,
    algebra::upgrading_arithmetic::{TryFromUpgrade, TrySmallestUpgrade},
    iterators::IteratableTensor,
    structure::{
        concrete_index::{ConcreteIndex, ExpandedIndex, FlatIndex},
        CastStructure, HasName, HasStructure, OrderedStructure, ScalarStructure, ScalarTensor,
        TensorStructure, TracksCount,
    },
};
use eyre::{eyre, Result};
use delegate::delegate;

#[cfg(feature = "shadowing")]
use crate::{
    shadowing::symbolica_utils::{atomic_expanded_label_id, IntoArgs, IntoSymbol},
    shadowing::{ShadowMapping, Shadowable},
    tensors::parametric::{ExpandedCoefficent, FlatCoefficent, TensorCoefficient},
};

use bincode::{Decode, Encode};
use indexmap::IndexMap;
use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use std::{borrow::Cow, fmt::Display, hash::Hash};

#[cfg(feature = "shadowing")]
use symbolica::{atom::Atom, atom::Symbol};

use super::{
    CastData, DataIterator, DenseTensor, GetTensorData, HasTensorData, SetTensorData, StorageTensor,
};

/// Sparse data tensor, generic on storage type `T`, and structure type `I`.
///
/// Stores data in a hashmap of usize, using ahash's hashmap.
/// The usize key is the flattened index of the corresponding position in the dense tensor
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Encode, Decode)]
pub struct SparseTensor<T, I = OrderedStructure> {
    // #[bincode(with_serde)]
    pub elements: std::collections::HashMap<FlatIndex, T>,
    pub zero: T,
    pub structure: I,
}

impl<T, S> crate::network::Ref for SparseTensor<T, S> {
    type Ref<'a>
        = &'a SparseTensor<T, S>
    where
        Self: 'a;

    fn refer(&self) -> Self::Ref<'_> {
        self
    }
}

impl<T: RefZero, S: TensorStructure> SparseTensor<T, S> {
    pub fn ref_zero(&self) -> T {
        self.zero.ref_zero()
    }
}

impl<Aind: AbsInd, T: Clone, S: Clone + Into<IndexLess<R, Aind>>, R: RepName<Dual = R>>
    PermuteTensor for SparseTensor<T, S>
where
    S: TensorStructure<Slot = Slot<R, Aind>> + PermuteTensor<IdSlot = Slot<R, Aind>, Id = S>,
{
    type Id = SparseTensor<T, S>;
    type IdSlot = (T, Slot<R, Aind>);
    type Permuted = SparseTensor<T, S>;

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
        let (zero, i) = i;
        let (one, j) = j;
        let s = S::id(i, j);
        let mut elements = std::collections::HashMap::new();
        for i in 0..usize::try_from(i.dim()).unwrap() {
            elements.insert(s.flat_index([i, i]).unwrap(), one.clone());
        }
        SparseTensor {
            elements,
            zero,
            structure: s,
        }
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

impl<T, S> TensorStructure for SparseTensor<T, S>
where
    S: TensorStructure,
{
    // type R = <T::Structure as TensorStructure>::R;
    type Indexed = SparseTensor<T, S::Indexed>;
    type Slot = S::Slot;

    fn reindex(
        self,
        indices: &[<Self::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::Indexed>, StructureError> {
        let res = self.structure.reindex(indices)?;

        Ok(PermutedStructure {
            structure: SparseTensor {
                zero: self.zero,
                structure: res.structure,
                elements: self.elements,
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

impl<T: Hash, I: Hash> Hash for SparseTensor<T, I> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let mut vecel: Vec<_> = self.elements.iter().collect();

        vecel.sort_by(|(i, _), (j, _)| i.cmp(j));

        vecel.hash(state);
        self.structure.hash(state);
    }
}

impl<T, U: From<T> + Clone, I: TensorStructure + Clone> CastData<SparseTensor<U, I>>
    for SparseTensor<T, I>
{
    fn cast_data(self) -> SparseTensor<U, I> {
        SparseTensor {
            zero: self.zero.into(),
            elements: self
                .elements
                .into_iter()
                .map(|(k, v)| (k, v.into()))
                .collect(),
            structure: self.structure,
        }
    }
}

impl<T, S: TensorStructure, O: From<S> + TensorStructure> CastStructure<SparseTensor<T, O>>
    for SparseTensor<T, S>
{
    fn cast_structure(self) -> SparseTensor<T, O> {
        SparseTensor {
            zero: self.zero,
            elements: self.elements,
            structure: self.structure.into(),
        }
    }
}

#[cfg(feature = "shadowing")]
impl<T: Clone, S: TensorStructure> Shadowable for SparseTensor<T, S>
where
    S: HasName + Clone,
    S::Name: IntoSymbol,
    S::Args: IntoArgs,
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
}

#[cfg(feature = "shadowing")]
impl<T: Clone, S: TensorStructure, R> ShadowMapping<R> for SparseTensor<T, S>
where
    S: HasName + Clone,
    S::Name: IntoSymbol,
    S::Args: IntoArgs,
    R: From<T>,
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    // fn shadow_with_map<'a, C>(
    //     &self,
    //     fn_map: &mut symbolica::evaluate::FunctionMap<'a, R>,
    //     index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> C,
    // ) -> Option<ParamTensor<Self::Structure>>
    // where
    //     C: TensorCoefficient,
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

impl<T, S> HasName for SparseTensor<T, S>
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

    #[cfg(feature = "shadowing")]
    fn expanded_coef(&self, id: FlatIndex) -> ExpandedCoefficent<Self::Args>
    where
        Self: TensorStructure,
        Self::Name: IntoSymbol,
        Self::Args: IntoArgs,
    {
        self.structure.expanded_coef(id)
    }

    #[cfg(feature = "shadowing")]
    fn flat_coef(&self, id: FlatIndex) -> FlatCoefficent<Self::Args>
    where
        Self: TensorStructure,
        Self::Name: IntoSymbol,
        Self::Args: IntoArgs,
    {
        self.structure.flat_coef(id)
    }
}

impl<T, I> HasTensorData for SparseTensor<T, I>
where
    T: Clone,
    I: TensorStructure + Clone,
{
    type Data = T;
    // type Storage = AHashMap<usize, T>;

    fn data(&self) -> Vec<T> {
        let mut d: Vec<(FlatIndex, T)> = self.iter_flat().map(|(i, v)| (i, v.clone())).collect();
        d.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        d.into_iter().map(|(_, v)| v).collect()
    }

    fn indices(&self) -> Vec<ExpandedIndex> {
        self.elements
            .keys()
            .map(|k| self.expanded_index(*k).unwrap())
            .collect()
    }

    fn hashmap(&self) -> IndexMap<ExpandedIndex, T> {
        let mut hashmap = IndexMap::new();
        for (k, v) in self.iter_expanded() {
            hashmap.insert(k.clone(), v.clone());
        }
        hashmap
    }
    #[cfg(feature = "shadowing")]
    fn symhashmap(&self, name: Symbol, args: &[Atom]) -> HashMap<Atom, T> {
        let mut hashmap = HashMap::new();

        for (k, v) in &self.elements {
            hashmap.insert(
                atomic_expanded_label_id(&self.expanded_index(*k).unwrap(), name, args),
                v.clone(),
            );
        }
        hashmap
    }
}

impl<T, S: TensorStructure> SparseTensor<T, S> {
    pub fn map_structure<S2>(self, f: impl Fn(S) -> S2) -> SparseTensor<T, S2>
    where
        S2: TensorStructure,
    {
        SparseTensor {
            zero: self.zero,
            elements: self.elements,
            structure: f(self.structure),
        }
    }

    pub fn map_structure_fallible<S2, E>(
        self,
        f: impl Fn(S) -> Result<S2, E>,
    ) -> Result<SparseTensor<T, S2>, E>
    where
        S2: TensorStructure,
    {
        Ok(SparseTensor {
            zero: self.zero,
            elements: self.elements,
            structure: f(self.structure)?,
        })
    }

    pub fn map_data_ref<U>(&self, f: impl Fn(&T) -> U) -> SparseTensor<U, S>
    where
        // T: Clone,
        // U: Clone,
        S: Clone,
    {
        let elements = self.flat_iter().map(|(k, v)| (k, f(v))).collect();
        SparseTensor {
            zero: f(&self.zero),
            elements,
            structure: self.structure.clone(),
        }
    }

    pub fn map_data_ref_result<U, E>(
        &self,
        f: impl Fn(&T) -> Result<U, E>,
    ) -> Result<SparseTensor<U, S>, E>
    where
        // T: Clone,
        // U: Clone,
        S: Clone,
    {
        let elements: Result<HashMap<FlatIndex, _>, E> = self
            .flat_iter()
            .map(|(k, v)| f(v).map(|v| (k, v)))
            .collect();
        Ok(SparseTensor {
            zero: f(&self.zero)?,
            elements: elements?,
            structure: self.structure.clone(),
        })
    }

    pub fn map_data<U>(self, f: impl Fn(T) -> U) -> SparseTensor<U, S> {
        let elements = self.elements.into_iter().map(|(k, v)| (k, f(v))).collect();
        SparseTensor {
            zero: f(self.zero),
            elements,
            structure: self.structure,
        }
    }

    pub fn map_data_ref_mut<U>(&mut self, mut f: impl FnMut(&mut T) -> U) -> SparseTensor<U, S>
    where
        // T: Clone,
        // U: Clone,
        S: Clone,
    {
        let elements = self.elements.iter_mut().map(|(k, v)| (*k, f(v))).collect();
        SparseTensor {
            zero: f(&mut self.zero),
            elements,
            structure: self.structure.clone(),
        }
    }

    pub fn map_data_mut(&mut self, f: impl FnMut(&mut T))
    where
        // T: Clone,
        // U: Clone,
        S: Clone,
    {
        self.elements.values_mut().for_each(f);
    }
}

impl<T, S> Display for SparseTensor<T, S>
where
    T: Display,
    S: TensorStructure,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut s = String::new();
        for (i, v) in self.iter_expanded() {
            s.push_str(&format!("{}: {}\n", i, v));
        }
        write!(f, "{}", s)
    }
}

impl<T, U, S> TryFromUpgrade<SparseTensor<U, S>> for SparseTensor<T, S>
where
    U: TrySmallestUpgrade<T, LCM = T>,
    S: TensorStructure + Clone,
    T: Clone,
{
    fn try_from_upgrade(data: &SparseTensor<U, S>) -> Option<SparseTensor<T, S>> {
        data.try_upgrade().map(Cow::into_owned)
    }
}

// #[derive(Error, Debug)]
// pub enum DataTensorError {
//     #[error("Data length does not match shape")]
//     DataLengthMismatch,
// }

impl<T, I> SetTensorData for SparseTensor<T, I>
where
    I: TensorStructure,
{
    type SetData = T;
    /// falible set method, returns an error if the indices are out of bounds.
    /// Does not check if the inserted value is zero.
    fn set(&mut self, indices: &[ConcreteIndex], value: T) -> Result<()> {
        self.verify_indices(indices)?;
        self.elements
            .insert(self.flat_index(indices).unwrap(), value);
        Ok(())
    }

    /// falible set given a flat index, returns an error if the indices are out of bounds.
    fn set_flat(&mut self, index: FlatIndex, value: T) -> Result<()> {
        if index >= self.size()?.into() {
            return Err(eyre!("Index out of bounds"));
        }
        self.elements.insert(index, value);
        Ok(())
    }
}
impl<T, I> GetTensorData for SparseTensor<T, I>
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
    fn get_ref<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<&T> {
        if let Ok(idx) = self.flat_index(&indices) {
            self.elements
                .get(&idx)
                .ok_or(eyre!("No elements at that spot"))
        } else if self.structure.is_scalar() && indices.as_ref().is_empty() {
            self.elements
                .iter()
                .next()
                .map(|(_, v)| v)
                .ok_or(eyre!("err"))
        } else {
            Err(eyre!("Index out of bounds"))
        }
    }

    fn get_ref_linear(&self, index: FlatIndex) -> Option<&T> {
        self.elements.get(&index)
    }

    fn get_mut_linear(&mut self, index: FlatIndex) -> Option<&mut T> {
        self.elements.get_mut(&index)
    }

    fn get_owned<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<Self::GetDataOwned>
    where
        T: Clone,
    {
        self.get_ref(indices).cloned()
    }

    fn get_owned_linear(&self, index: FlatIndex) -> Option<Self::GetDataOwned>
    where
        T: Clone,
    {
        self.elements.get(&index).cloned()
    }
}

impl<T: RefZero, I> ScalarTensor for SparseTensor<T, I>
where
    I: TensorStructure + ScalarStructure,
{
    fn new_scalar(scalar: Self::Scalar) -> Self {
        let mut elements = HashMap::new();
        let zero = scalar.ref_zero();
        elements.insert(0.into(), scalar);
        SparseTensor {
            zero,
            elements,
            structure: I::scalar_structure(),
        }
    }
}

impl<T, I> HasStructure for SparseTensor<T, I>
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
        = SparseTensor<T, S>
    where
        S: TensorStructure;

    fn map_structure<O: TensorStructure>(
        self,
        f: impl FnOnce(Self::Structure) -> O,
    ) -> Self::Store<O> {
        SparseTensor {
            zero: self.zero,
            structure: f(self.structure),
            elements: self.elements,
        }
    }

    fn map_structure_result<O: TensorStructure, Er>(
        self,
        f: impl FnOnce(Self::Structure) -> Result<O, Er>,
    ) -> std::result::Result<Self::Store<O>, Er> {
        Ok(SparseTensor {
            zero: self.zero,
            structure: f(self.structure)?,
            elements: self.elements,
        })
    }
    fn structure(&self) -> &Self::Structure {
        &self.structure
    }

    fn map_same_structure(self, f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self {
        SparseTensor {
            zero: self.zero,
            elements: self.elements,
            structure: f(self.structure),
        }
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        &mut self.structure
    }

    fn scalar(mut self) -> Option<Self::Scalar> {
        if self.structure.is_scalar() {
            self.elements.drain().next().map(|(_, v)| v)
        } else {
            None
        }
    }

    fn scalar_ref(&self) -> Option<Self::ScalarRef<'_>> {
        if self.structure.is_scalar() {
            self.elements.values().next()
        } else {
            None
        }
    }
}

impl<T, I> TracksCount for SparseTensor<T, I>
where
    I: TracksCount + TensorStructure,
{
    fn contractions_num(&self) -> usize {
        self.structure.contractions_num()
    }
}

impl<U, I> SparseTensor<U, I>
where
    I: TensorStructure + Clone,
{
    pub fn try_upgrade<T>(&self) -> Option<Cow<'_, SparseTensor<U::LCM, I>>>
    where
        U: TrySmallestUpgrade<T>,
        U::LCM: Clone,
    {
        let structure = self.structure.clone();
        let elements: Option<HashMap<FlatIndex, U::LCM>> = self
            .elements
            .iter()
            .map(|(k, v)| match v.try_upgrade() {
                Some(Cow::Owned(u)) => Some((*k, u)),
                Some(Cow::Borrowed(u)) => Some((*k, u.clone())),
                None => None,
            })
            .collect();
        Some(Cow::Owned(SparseTensor {
            zero: self.zero.try_upgrade()?.into_owned(),
            elements: elements?,
            structure,
        }))
    }
}

impl<T, I> SparseTensor<T, I>
where
    I: TensorStructure,
{
    /// Create a new empty sparse tensor with the given structure
    pub fn empty(structure: I, zero: T) -> Self {
        SparseTensor {
            zero,
            elements: HashMap::default(),
            structure,
        }
    }

    /// Checks if there is a value at the given indices
    pub fn is_empty_at_expanded(&self, indices: &[ConcreteIndex]) -> bool {
        !self
            .elements
            .contains_key(&self.flat_index(indices).unwrap())
    }

    pub fn is_empty_at_flat(&self, index: FlatIndex) -> bool {
        !self.elements.contains_key(&index)
    }
    /// Calulates how dense the tensor is, i.e. the ratio of non-zero elements to total elements
    pub fn density(&self) -> f64 {
        f64::from(self.elements.len() as u32) / f64::from(self.size().unwrap() as u32)
    }

    /// Converts the sparse tensor to a dense tensor, with the same structure
    pub fn to_dense(&self) -> DenseTensor<T, I>
    where
        T: Clone,
        I: Clone,
    {
        self.to_dense_with(&self.zero)
    }

    pub fn to_dense_with(&self, zero: &T) -> DenseTensor<T, I>
    where
        T: Clone,
        I: Clone,
    {
        let mut dense = DenseTensor::fill(self.structure.clone(), zero.clone());
        for (indices, value) in self.elements.iter() {
            let _ = dense.set_flat(*indices, value.clone());
        }
        dense
    }

    /// fallible smart set method, returns an error if the indices are out of bounds.
    /// If the value is zero, it removes the element at the given indices.
    pub fn smart_set(&mut self, indices: &[ConcreteIndex], value: T) -> Result<()>
    where
        T: IsZero,
    {
        self.verify_indices(indices)?;
        if value.is_zero() {
            _ = self.elements.remove(&self.flat_index(indices).unwrap());
            return Ok(());
        }
        self.elements
            .insert(self.flat_index(indices).unwrap(), value);
        Ok(())
    }

    /// Generates a new sparse tensor from the given data and structure
    pub fn from_data(
        data: impl IntoIterator<Item = (Vec<ConcreteIndex>, T)>,
        structure: I,
        zero: T,
    ) -> Result<Self> {
        let mut elements = HashMap::default();
        for (index, value) in data {
            if index.len() != structure.order() {
                return Err(eyre!("Mismatched order"));
            }
            elements.insert(structure.flat_index(&index).unwrap(), value);
        }

        Ok(SparseTensor {
            zero,
            elements,
            structure,
        })
    }

    /// fallible smart get method, returns an error if the indices are out of bounds.
    /// If the index is in the bTree return the value, else return zero.
    pub fn smart_get(&self, indices: &[ConcreteIndex]) -> Result<Cow<'_, T>>
    where
        T: Default + Clone,
    {
        self.verify_indices(indices)?;
        // if the index is in the bTree return the value, else return default, lazily allocating the default
        Ok(
            match self.elements.get(&self.flat_index(indices).unwrap()) {
                Some(value) => Cow::Borrowed(value),
                None => Cow::Owned(T::default()),
            },
        )
    }
}

impl<T, I> SparseTensor<T, I>
where
    I: Clone + TensorStructure,
{
    pub fn convert_to<U>(&self) -> SparseTensor<U, I>
    where
        U: for<'a> From<&'a T>,
    {
        let elements = self.elements.iter().map(|(k, v)| (*k, v.into())).collect();
        SparseTensor {
            zero: (&self.zero).into(),
            elements,
            structure: self.structure.clone(),
        }
    }
}

impl<S: TensorStructure + Clone, T> StorageTensor for SparseTensor<T, S> {
    type Data = T;
    type ContainerData<Data> = SparseTensor<Data, S>;

    fn map_data_self(self, f: impl Fn(Self::Data) -> Self::Data) -> Self {
        self.map_data(f)
    }

    fn map_data_ref_mut_result<U, E>(
        &mut self,
        mut f: impl FnMut(&mut Self::Data) -> Result<U, E>,
    ) -> Result<Self::ContainerData<U>, E> {
        let elements: Result<HashMap<FlatIndex, _>, E> = self
            .elements
            .iter_mut()
            .map(|(k, v)| f(v).map(|v| (*k, v)))
            .collect();
        Ok(SparseTensor {
            zero: f(&mut self.zero)?,
            elements: elements?,
            structure: self.structure.clone(),
        })
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

    fn map_data_ref<U>(&self, f: impl Fn(&T) -> U) -> SparseTensor<U, S> {
        let elements = self.flat_iter().map(|(k, v)| (k, f(v))).collect();
        SparseTensor {
            zero: f(&self.zero),
            elements,
            structure: self.structure.clone(),
        }
    }

    fn map_data_ref_result<U, E>(
        &self,
        f: impl Fn(&T) -> Result<U, E>,
    ) -> Result<SparseTensor<U, S>, E> {
        let elements: Result<HashMap<FlatIndex, _>, E> = self
            .flat_iter()
            .map(|(k, v)| f(v).map(|v| (k, v)))
            .collect();
        Ok(SparseTensor {
            zero: f(&self.zero)?,
            elements: elements?,
            structure: self.structure.clone(),
        })
    }

    fn map_data<U>(self, f: impl Fn(T) -> U) -> SparseTensor<U, S> {
        let elements = self.elements.into_iter().map(|(k, v)| (k, f(v))).collect();
        SparseTensor {
            zero: f(self.zero),
            elements,
            structure: self.structure,
        }
    }

    fn map_data_ref_mut<U>(&mut self, mut f: impl FnMut(&mut T) -> U) -> SparseTensor<U, S> {
        let elements = self.elements.iter_mut().map(|(k, v)| (*k, f(v))).collect();
        SparseTensor {
            zero: f(&mut self.zero),
            elements,
            structure: self.structure.clone(),
        }
    }

    fn map_data_mut(&mut self, f: impl FnMut(&mut T)) {
        self.elements.values_mut().for_each(f);
    }
}
