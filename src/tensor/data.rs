use super::{
    atomic_expanded_label_id, ConcreteIndex, DenseTensorLinearIterator, HasName, Slot,
    SparseTensorLinearIterator, TensorStructure, TracksCount, TrySmallestUpgrade, VecStructure,
};
use ahash::AHashMap;
use derive_more::From;
use enum_try_as_inner::EnumTryAsInner;
use indexmap::IndexMap;
use num::Zero;

use serde::{Deserialize, Serialize};
use smartstring::alias::String;
use std::{borrow::Cow, collections::HashMap};
use symbolica::domains::float::Complex;
use symbolica::{representations::Atom, representations::Symbol};

pub trait DataIterator<T> {
    type FlatIter<'a>: Iterator<Item = (usize, &'a T)>
    where
        Self: 'a,
        T: 'a;

    fn flat_iter(&self) -> Self::FlatIter<'_>;
}

impl<T, I> DataIterator<T> for SparseTensor<T, I> {
    type FlatIter<'a> = SparseTensorLinearIterator<'a, T> where  I:'a,T: 'a;

    fn flat_iter(&self) -> Self::FlatIter<'_> {
        SparseTensorLinearIterator::new(self)
    }
}

impl<T, I: TensorStructure> DataIterator<T> for DenseTensor<T, I> {
    type FlatIter<'a> = DenseTensorLinearIterator<'a,T,I> where  I:'a,T: 'a;

    fn flat_iter(&self) -> Self::FlatIter<'_> {
        DenseTensorLinearIterator::new(self)
    }
}

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

impl<T> Settable for AHashMap<usize, T> {
    type SetData = T;
    fn set(&mut self, index: usize, data: T) {
        self.insert(index, data);
    }
}

/// Trait for getting the data of a tensor
pub trait HasTensorData: TensorStructure {
    type Data: Clone;
    // type Storage: Settable<SetData = Self::Data>;
    /// Returns all the data in the tensor, withouth any structure
    fn data(&self) -> Vec<Self::Data>;

    /// Returns all the indices of the tensor, the order of the indices is the same as the order of the data
    fn indices(&self) -> Vec<Vec<ConcreteIndex>>;

    /// Returns a hashmap of the data, with the (expanded) indices as keys
    fn hashmap(&self) -> IndexMap<Vec<ConcreteIndex>, Self::Data>;

    /// Returns a hashmap of the data, with the the shadowed indices as keys
    fn symhashmap(&self, id: Symbol) -> HashMap<Atom, Self::Data>;
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
    fn set(&mut self, indices: &[ConcreteIndex], value: Self::SetData) -> Result<(), String>;

    fn set_flat(&mut self, index: usize, value: Self::SetData) -> Result<(), String>;
}

/// Trait for getting the data of a tensor
pub trait GetTensorData {
    type GetData;

    fn get(&self, indices: &[ConcreteIndex]) -> Result<&Self::GetData, String>;

    fn get_linear(&self, index: usize) -> Option<&Self::GetData>;
}

/// Sparse data tensor, generic on storage type `T`, and structure type `I`.  
///
/// Stores data in a hashmap of usize, using ahash's hashmap.
/// The usize key is the flattened index of the corresponding position in the dense tensor
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SparseTensor<T, I = VecStructure> {
    pub elements: AHashMap<usize, T>,
    pub structure: I,
}

impl<T, I> HasTensorData for SparseTensor<T, I>
where
    T: Clone,
    I: TensorStructure,
{
    type Data = T;
    // type Storage = AHashMap<usize, T>;

    fn data(&self) -> Vec<T> {
        self.elements.values().cloned().collect()
    }

    fn indices(&self) -> Vec<Vec<ConcreteIndex>> {
        self.elements
            .keys()
            .map(|k| self.expanded_index(*k).unwrap())
            .collect()
    }

    fn hashmap(&self) -> IndexMap<Vec<ConcreteIndex>, T> {
        let mut hashmap = IndexMap::new();
        for (k, v) in self.iter() {
            hashmap.insert(k.clone(), v.clone());
        }
        hashmap
    }

    fn symhashmap(&self, id: Symbol) -> HashMap<Atom, T> {
        let mut hashmap = HashMap::new();

        for (k, v) in &self.elements {
            hashmap.insert(
                atomic_expanded_label_id(&self.expanded_index(*k).unwrap(), id),
                v.clone(),
            );
        }
        hashmap
    }
}

impl<T, I> SetTensorData for SparseTensor<T, I>
where
    I: TensorStructure,
{
    type SetData = T;
    /// falible set method, returns an error if the indices are out of bounds.
    /// Does not check if the inserted value is zero.
    fn set(&mut self, indices: &[ConcreteIndex], value: T) -> Result<(), String> {
        self.verify_indices(indices)?;
        self.elements
            .insert(self.flat_index(indices).unwrap(), value);
        Ok(())
    }

    /// falible set given a flat index, returns an error if the indices are out of bounds.
    fn set_flat(&mut self, index: usize, value: T) -> Result<(), String> {
        if index >= self.size() {
            return Err("Index out of bounds".into());
        }
        self.elements.insert(index, value);
        Ok(())
    }
}
impl<T, I> GetTensorData for SparseTensor<T, I>
where
    I: TensorStructure,
{
    type GetData = T;
    fn get(&self, indices: &[ConcreteIndex]) -> Result<&T, String> {
        if let Ok(idx) = self.flat_index(indices) {
            self.elements
                .get(&idx)
                .ok_or("No elements at that spot".into())
        } else if self.structure.is_scalar() && indices.is_empty() {
            self.elements
                .iter()
                .next()
                .map(|(_, v)| v)
                .ok_or("err".into())
        } else {
            Err("Index out of bounds".into())
        }
    }

    fn get_linear(&self, index: usize) -> Option<&T> {
        self.elements.get(&index)
    }
}

impl<T, I> TensorStructure for SparseTensor<T, I>
where
    I: TensorStructure,
{
    type Structure = I;
    fn structure(&self) -> &Self::Structure {
        &self.structure
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        &mut self.structure
    }
    fn external_structure(&self) -> &[Slot] {
        self.structure.external_structure()
    }
}

impl<T, I> HasName for SparseTensor<T, I>
where
    I: HasName,
{
    type Name = I::Name;

    fn name(&self) -> Option<Cow<'_, <I as HasName>::Name>> {
        self.structure.name()
    }
    fn set_name(&mut self, name: &Self::Name) {
        self.structure.set_name(name);
    }
}

impl<T, I> TracksCount for SparseTensor<T, I>
where
    I: TracksCount,
{
    fn contractions_num(&self) -> usize {
        self.structure.contractions_num()
    }
}

impl<T, U, I> TrySmallestUpgrade<SparseTensor<T, I>> for SparseTensor<U, I>
where
    U: TrySmallestUpgrade<T>,
    U::LCM: Clone,
    I: TensorStructure + Clone,
{
    type LCM = SparseTensor<U::LCM, I>;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
    where
        Self::LCM: Clone,
    {
        let structure = self.structure.clone();
        let elements: Option<AHashMap<usize, U::LCM>> = self
            .elements
            .iter()
            .map(|(k, v)| match v.try_upgrade() {
                Some(Cow::Owned(u)) => Some((*k, u)),
                Some(Cow::Borrowed(u)) => Some((*k, u.clone())),
                None => None,
            })
            .collect();
        Some(Cow::Owned(SparseTensor {
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
    pub fn empty(structure: I) -> Self {
        SparseTensor {
            elements: AHashMap::default(),
            structure,
        }
    }

    /// Checks if there is a value at the given indices
    pub fn is_empty_at(&self, indices: &[ConcreteIndex]) -> bool {
        !self
            .elements
            .contains_key(&self.flat_index(indices).unwrap())
    }

    /// Calulates how dense the tensor is, i.e. the ratio of non-zero elements to total elements
    pub fn density(&self) -> f64 {
        f64::from(self.elements.len() as u32) / f64::from(self.size() as u32)
    }

    /// Converts the sparse tensor to a dense tensor, with the same structure
    pub fn to_dense(&self) -> DenseTensor<T, I>
    where
        T: Clone + Default,
        I: Clone,
    {
        let mut dense = DenseTensor::default(self.structure.clone());
        for (indices, value) in self.elements.iter() {
            let _ = dense.set_flat(*indices, value.clone());
        }
        dense
    }

    /// fallible smart set method, returns an error if the indices are out of bounds.
    /// If the value is zero, it removes the element at the given indices.
    pub fn smart_set(&mut self, indices: &[ConcreteIndex], value: T) -> Result<(), String>
    where
        T: Default + PartialEq,
    {
        self.verify_indices(indices)?;
        if value == T::default() {
            _ = self.elements.remove(&self.flat_index(indices).unwrap());
            return Ok(());
        }
        self.elements
            .insert(self.flat_index(indices).unwrap(), value);
        Ok(())
    }

    /// Generates a new sparse tensor from the given data and structure
    pub fn from_data(data: &[(Vec<ConcreteIndex>, T)], structure: I) -> Result<Self, String>
    where
        T: Clone,
    {
        let mut dimensions = vec![0; structure.order()];
        for (index, _) in data {
            if index.len() != structure.order() {
                return Err("Mismatched order".into());
            }
            for (i, &idx) in index.iter().enumerate() {
                if idx >= dimensions[i] {
                    dimensions[i] = idx + 1;
                }
            }
        }
        let mut elements = AHashMap::default();
        for (index, value) in data {
            elements.insert(structure.flat_index(index).unwrap(), value.clone());
        }

        Ok(SparseTensor {
            elements,
            structure,
        })
    }

    /// fallible smart get method, returns an error if the indices are out of bounds.
    /// If the index is in the bTree return the value, else return zero.
    pub fn smart_get(&self, indices: &[ConcreteIndex]) -> Result<Cow<T>, String>
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

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct DenseTensor<T, I = VecStructure> {
    pub data: Vec<T>,
    pub structure: I,
}

impl<T, I> TensorStructure for DenseTensor<T, I>
where
    I: TensorStructure,
{
    type Structure = I;
    fn structure(&self) -> &Self::Structure {
        &self.structure
    }
    fn mut_structure(&mut self) -> &mut Self::Structure {
        &mut self.structure
    }
    fn external_structure(&self) -> &[Slot] {
        self.structure.external_structure()
    }
}

impl<T, I> HasName for DenseTensor<T, I>
where
    I: HasName,
{
    type Name = I::Name;

    fn name(&self) -> Option<Cow<'_, <I as HasName>::Name>> {
        self.structure.name()
    }
    fn set_name(&mut self, name: &Self::Name) {
        self.structure.set_name(name);
    }
}

impl<T, I> TracksCount for DenseTensor<T, I>
where
    I: TracksCount,
{
    fn contractions_num(&self) -> usize {
        self.structure.contractions_num()
    }
}

impl<T: Default + Clone, I> DenseTensor<T, I>
where
    I: TensorStructure,
{
    pub fn default(structure: I) -> Self {
        let length = if structure.is_scalar() {
            1
        } else {
            structure.size()
        };
        DenseTensor {
            data: vec![T::default(); length],
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
            structure.size()
        };
        DenseTensor {
            data: vec![T::zero(); length],
            structure,
        }
    }
}

impl<T, U, I> TrySmallestUpgrade<DenseTensor<T, I>> for DenseTensor<U, I>
where
    U: TrySmallestUpgrade<T>,
    U::LCM: Clone,
    I: TensorStructure + Clone,
{
    type LCM = DenseTensor<U::LCM, I>;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
    where
        Self::LCM: Clone,
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
    pub fn from_data(data: &[T], structure: I) -> Result<Self, String> {
        if data.len() != structure.size() && !(data.len() == 1 && structure.is_scalar()) {
            return Err("Data length does not match shape".into());
        }
        Ok(DenseTensor {
            data: data.to_vec(),
            structure,
        })
    }

    /// Generates a new dense tensor from the given data and structure, truncating the data if it is too long with respect to the structure
    pub fn from_data_coerced(data: &[T], structure: I) -> Result<Self, String> {
        if data.len() < structure.size() {
            return Err("Data length is too small".into());
        }
        let mut data = data.to_vec();
        if structure.is_scalar() {
            data.truncate(1);
        } else {
            data.truncate(structure.size());
        }
        Ok(DenseTensor { data, structure })
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
        let mut sparse = SparseTensor::empty(self.structure.clone());
        for (i, value) in self.iter() {
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
    I: Clone,
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

impl<T, I> SparseTensor<T, I>
where
    I: Clone,
{
    pub fn convert_to<U>(&self) -> SparseTensor<U, I>
    where
        U: for<'a> From<&'a T>,
    {
        let elements = self.elements.iter().map(|(k, v)| (*k, v.into())).collect();
        SparseTensor {
            elements,
            structure: self.structure.clone(),
        }
    }
}

impl<T, I> HasTensorData for DenseTensor<T, I>
where
    T: Clone,
    I: TensorStructure,
{
    type Data = T;
    fn data(&self) -> Vec<T> {
        self.data.clone()
    }

    fn indices(&self) -> Vec<Vec<ConcreteIndex>> {
        let mut indices = Vec::new();
        for i in 0..self.size() {
            indices.push(self.expanded_index(i).unwrap());
        }
        indices
    }

    fn hashmap(&self) -> IndexMap<Vec<ConcreteIndex>, T> {
        let mut hashmap = IndexMap::new();
        for (k, v) in self.iter() {
            hashmap.insert(k.clone(), v.clone());
        }
        hashmap
    }

    fn symhashmap(&self, id: Symbol) -> HashMap<Atom, T> {
        let mut hashmap = HashMap::new();

        for (k, v) in self.iter() {
            hashmap.insert(atomic_expanded_label_id(&k, id), v.clone());
        }
        hashmap
    }
}

impl<T, I> SetTensorData for DenseTensor<T, I>
where
    I: TensorStructure,
{
    type SetData = T;
    fn set(&mut self, indices: &[ConcreteIndex], value: T) -> Result<(), String> {
        self.verify_indices(indices)?;
        let idx = self.flat_index(indices);
        if let Ok(i) = idx {
            self.data[i] = value;
        }
        Ok(())
    }

    fn set_flat(&mut self, index: usize, value: T) -> Result<(), String> {
        if index < self.size() {
            self.data[index] = value;
        } else {
            return Err("Index out of bounds".into());
        }
        Ok(())
    }
}

impl<T, I> GetTensorData for DenseTensor<T, I>
where
    I: TensorStructure,
{
    type GetData = T;
    fn get_linear(&self, index: usize) -> Option<&T> {
        self.data.get(index)
    }

    fn get(&self, indices: &[ConcreteIndex]) -> Result<&T, String> {
        if let Ok(idx) = self.flat_index(indices) {
            Ok(&self.data[idx])
        } else if self.structure.is_scalar() && indices.is_empty() {
            Ok(&self.data[0])
        } else {
            Err("Index out of bounds".into())
        }
    }
}

/// Enum for storing either a dense or a sparse tensor, with the same structure
#[derive(Debug, Clone, EnumTryAsInner, Serialize, Deserialize, From)]
#[derive_err(Debug)]
pub enum DataTensor<T, I: TensorStructure = VecStructure> {
    Dense(DenseTensor<T, I>),
    Sparse(SparseTensor<T, I>),
}

impl<T, I> DataTensor<T, I>
where
    I: TensorStructure + Clone,
{
    pub fn to_sparse(self) -> SparseTensor<T, I>
    where
        T: Clone + Default + PartialEq,
    {
        match self {
            DataTensor::Dense(d) => d.to_sparse(),
            DataTensor::Sparse(s) => s,
        }
    }

    pub fn to_dense(self) -> DenseTensor<T, I>
    where
        T: Clone + Default + PartialEq,
    {
        match self {
            DataTensor::Dense(d) => d,
            DataTensor::Sparse(s) => s.to_dense(),
        }
    }
}

impl<T, I> HasTensorData for DataTensor<T, I>
where
    I: TensorStructure,
    T: Clone,
{
    type Data = T;
    fn data(&self) -> Vec<T> {
        match self {
            DataTensor::Dense(d) => d.data(),
            DataTensor::Sparse(s) => s.data(),
        }
    }

    fn indices(&self) -> Vec<Vec<ConcreteIndex>> {
        match self {
            DataTensor::Dense(d) => d.indices(),
            DataTensor::Sparse(s) => s.indices(),
        }
    }

    fn hashmap(&self) -> IndexMap<Vec<ConcreteIndex>, T> {
        match self {
            DataTensor::Dense(d) => d.hashmap(),
            DataTensor::Sparse(s) => s.hashmap(),
        }
    }

    fn symhashmap(&self, id: Symbol) -> HashMap<Atom, T> {
        match self {
            DataTensor::Dense(d) => d.symhashmap(id),
            DataTensor::Sparse(s) => s.symhashmap(id),
        }
    }
}

impl<T, I> TensorStructure for DataTensor<T, I>
where
    I: TensorStructure,
    T: Clone,
{
    type Structure = I;
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
    fn external_structure(&self) -> &[Slot] {
        match self {
            DataTensor::Dense(d) => d.external_structure(),

            DataTensor::Sparse(s) => s.external_structure(),
        }
    }
}

impl<T, I> HasName for DataTensor<T, I>
where
    I: HasName,
    T: Clone,
    I: TensorStructure,
{
    type Name = I::Name;

    fn name(&self) -> Option<Cow<'_, <I as HasName>::Name>> {
        match self {
            DataTensor::Dense(d) => d.name(),
            DataTensor::Sparse(s) => s.name(),
        }
    }
    fn set_name(&mut self, name: &Self::Name) {
        match self {
            DataTensor::Dense(d) => d.set_name(name),
            DataTensor::Sparse(s) => s.set_name(name),
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

impl<T, U> TrySmallestUpgrade<DataTensor<T>> for DataTensor<U>
where
    U: TrySmallestUpgrade<T>,
    U::LCM: Clone,
{
    type LCM = DataTensor<U::LCM>;
    fn try_upgrade(&self) -> Option<Cow<Self::LCM>>
    where
        Self::LCM: Clone,
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

impl<T, S> SetTensorData for DataTensor<T, S>
where
    S: TensorStructure,
{
    type SetData = T;

    fn set(&mut self, indices: &[ConcreteIndex], value: Self::SetData) -> Result<(), String> {
        match self {
            DataTensor::Dense(d) => d.set(indices, value),
            DataTensor::Sparse(s) => s.set(indices, value),
        }
    }

    fn set_flat(&mut self, index: usize, value: Self::SetData) -> Result<(), String> {
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
    type GetData = T;

    fn get(&self, indices: &[ConcreteIndex]) -> Result<&Self::GetData, String> {
        match self {
            DataTensor::Dense(d) => d.get(indices),
            DataTensor::Sparse(s) => s.get(indices),
        }
    }

    fn get_linear(&self, index: usize) -> Option<&Self::GetData> {
        match self {
            DataTensor::Dense(d) => d.get_linear(index),
            DataTensor::Sparse(s) => s.get_linear(index),
        }
    }
}

/// Enum for a datatensor with specific numeric data type, generic on the structure type `I`
#[derive(Debug, Clone, EnumTryAsInner, Serialize, Deserialize)]
#[derive_err(Debug)]
pub enum NumTensor<T: TensorStructure = VecStructure> {
    Float(DataTensor<f64, T>),
    Complex(DataTensor<Complex<f64>, T>),
}

impl<T> TensorStructure for NumTensor<T>
where
    T: TensorStructure,
{
    type Structure = T;
    fn structure(&self) -> &Self::Structure {
        match self {
            NumTensor::Float(f) => f.structure(),
            NumTensor::Complex(c) => c.structure(),
        }
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        match self {
            NumTensor::Float(f) => f.mut_structure(),
            NumTensor::Complex(c) => c.mut_structure(),
        }
    }

    fn external_structure(&self) -> &[Slot] {
        match self {
            NumTensor::Float(f) => f.external_structure(),
            NumTensor::Complex(c) => c.external_structure(),
        }
    }
}

impl<T> HasName for NumTensor<T>
where
    T: HasName + TensorStructure,
{
    type Name = T::Name;

    fn name(&self) -> Option<Cow<'_, <T as HasName>::Name>> {
        match self {
            NumTensor::Float(f) => f.name(),
            NumTensor::Complex(c) => c.name(),
        }
    }
    fn set_name(&mut self, name: &Self::Name) {
        match self {
            NumTensor::Float(f) => f.set_name(name),
            NumTensor::Complex(c) => c.set_name(name),
        }
    }
}

impl<T> TracksCount for NumTensor<T>
where
    T: TracksCount + TensorStructure,
{
    fn contractions_num(&self) -> usize {
        match self {
            NumTensor::Float(f) => f.contractions_num(),
            NumTensor::Complex(c) => c.contractions_num(),
        }
    }
}

impl<T> From<DenseTensor<f64, T>> for NumTensor<T>
where
    T: TensorStructure,
{
    fn from(other: DenseTensor<f64, T>) -> Self {
        NumTensor::Float(DataTensor::Dense(other))
    }
}

impl<T> From<SparseTensor<f64, T>> for NumTensor<T>
where
    T: TensorStructure,
{
    fn from(other: SparseTensor<f64, T>) -> Self {
        NumTensor::Float(DataTensor::Sparse(other))
    }
}

impl<T> From<DenseTensor<Complex<f64>, T>> for NumTensor<T>
where
    T: TensorStructure,
{
    fn from(other: DenseTensor<Complex<f64>, T>) -> Self {
        NumTensor::Complex(DataTensor::Dense(other))
    }
}

impl<T> From<SparseTensor<Complex<f64>, T>> for NumTensor<T>
where
    T: TensorStructure,
{
    fn from(other: SparseTensor<Complex<f64>, T>) -> Self {
        NumTensor::Complex(DataTensor::Sparse(other))
    }
}
