use super::{AbstractIndex, ConcreteIndex, Dimension, HasTensorStructure, TensorSkeleton};
use crate::tensor::IntoId;
use ahash::AHashMap;
use enum_dispatch::enum_dispatch;
use enum_try_as_inner::EnumTryAsInner;
use num::Complex;
use std::{borrow::Cow, collections::HashMap};
use symbolica::{
    representations::Atom,
    representations::{Identifier, Num},
    state::{State, Workspace},
};

#[derive(Debug, Clone)]
pub struct SparseTensor<T, I = String> {
    pub elements: AHashMap<Vec<ConcreteIndex>, T>,
    pub structure: TensorSkeleton<I>,
}

impl<T, I> HasTensorStructure for SparseTensor<T, I> {
    type Name = I;
    fn structure(&self) -> &TensorSkeleton<I> {
        &self.structure
    }

    fn mut_structure(&mut self) -> &mut TensorSkeleton<I> {
        &mut self.structure
    }
}

impl<T, I> SparseTensor<T, I> {
    pub fn empty(structure: TensorSkeleton<I>) -> Self {
        SparseTensor {
            elements: AHashMap::new(),
            structure,
        }
    }

    pub fn empty_from_integers(slots: &[(AbstractIndex, Dimension)], name: I) -> Self
    where
        I: Clone,
    {
        let structure = TensorSkeleton::<I>::from_integers(slots, name);
        SparseTensor {
            elements: AHashMap::new(),
            structure,
        }
    }

    pub fn is_empty_at(&self, indices: &[ConcreteIndex]) -> bool {
        !self.elements.contains_key(indices)
    }

    pub fn set(&mut self, indices: &[ConcreteIndex], value: T) -> Result<(), String> {
        self.verify_indices(indices)?;
        self.elements.insert(indices.to_vec(), value);
        Ok(())
    }

    pub fn density(&self) -> f64 {
        f64::from(self.elements.len() as u32) / f64::from(self.size() as u32)
    }
}

impl<T, I> SparseTensor<T, I> {
    pub fn to_dense(&self) -> DenseTensor<T, I>
    where
        T: Clone + Default,
        I: Clone,
    {
        let mut dense = DenseTensor::default(self.structure.clone());
        for (indices, value) in &self.elements {
            dense.set(indices, value.clone());
        }
        dense
    }
}

impl<T, I> SparseTensor<T, I>
where
    T: Default + std::cmp::PartialEq,
{
    pub fn smart_set(&mut self, indices: &[ConcreteIndex], value: T) -> Result<(), String> {
        self.verify_indices(indices)?;
        if value == T::default() {
            _ = self.elements.remove(indices);
            return Ok(());
        }
        self.elements.insert(indices.to_vec(), value);
        Ok(())
    }
}

impl<T, I> SparseTensor<T, I>
where
    T: Clone,
{
    pub fn from_data(
        data: &[(Vec<ConcreteIndex>, T)],
        structure: TensorSkeleton<I>,
    ) -> Result<Self, String> {
        let mut dimensions = vec![0; structure.order()];
        for (index, _) in data {
            if index.len() != structure.order() {
                return Err("Mismatched order".to_string());
            }
            for (i, &idx) in index.iter().enumerate() {
                if idx >= dimensions[i] {
                    dimensions[i] = idx + 1;
                }
            }
        }
        Ok(SparseTensor {
            elements: AHashMap::from_iter(data.iter().cloned()),
            structure,
        })
    }
}

impl<T, I> SparseTensor<T, I>
where
    T: Clone + Default,
    I: Clone,
{
    pub fn get(&self, indices: &[ConcreteIndex]) -> Result<&T, String> {
        self.verify_indices(indices)?;
        self.elements
            .get(indices)
            .ok_or("No elements at that spot".to_string())
    }
    pub fn get_with_defaults(&self, indices: &[ConcreteIndex]) -> Result<Cow<T>, String> {
        self.verify_indices(indices)?;
        // if the index is in the bTree return the value, else return default, lazily allocating the default
        Ok(match self.elements.get(indices) {
            Some(value) => Cow::Borrowed(value),
            None => Cow::Owned(T::default()),
        })
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct DenseTensor<T, I = String> {
    pub data: Vec<T>,
    pub structure: TensorSkeleton<I>,
}

impl<T, I> HasTensorStructure for DenseTensor<T, I> {
    type Name = I;
    fn structure(&self) -> &TensorSkeleton<I> {
        &self.structure
    }
    fn mut_structure(&mut self) -> &mut TensorSkeleton<I> {
        &mut self.structure
    }
}

impl<T: Default + Clone, I> DenseTensor<T, I> {
    pub fn default(structure: TensorSkeleton<I>) -> Self {
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

    pub fn default_from_integers(slots: &[(AbstractIndex, Dimension)], name: I) -> Self
    where
        I: Clone,
    {
        let structure = TensorSkeleton::<I>::from_integers(slots, name);
        DenseTensor::default(structure)
    }
}

impl<T: Clone, I> DenseTensor<T, I> {
    pub fn from_data(data: &[T], structure: TensorSkeleton<I>) -> Result<Self, String> {
        if data.len() != structure.size() && !(data.len() == 1 && structure.is_scalar()) {
            return Err("Data length does not match shape".to_string());
        }
        Ok(DenseTensor {
            data: data.to_vec(),
            structure,
        })
    }

    pub fn from_data_coerced(data: &[T], structure: TensorSkeleton<I>) -> Result<Self, String> {
        if data.len() < structure.size() {
            return Err("Data length is too small".to_string());
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

impl<T, I> DenseTensor<T, I> {
    pub fn set(&mut self, indices: &[ConcreteIndex], value: T) {
        let idx = self.flat_index(indices);
        if let Ok(i) = idx {
            self.data[i] = value;
        }
    }

    pub fn get_linear(&self, index: usize) -> Option<&T> {
        self.data.get(index)
    }

    pub fn get(&self, indices: &[ConcreteIndex]) -> Option<&T> {
        if let Ok(idx) = self.flat_index(indices) {
            Some(&self.data[idx])
        } else if self.structure.is_scalar() && indices.is_empty() {
            Some(&self.data[0])
        } else {
            None
        }
    }
}

impl<T, I> DenseTensor<T, I>
where
    I: Clone,
{
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
    pub fn conver_to<U>(&self) -> DenseTensor<U, I>
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
        let elements = self
            .elements
            .iter()
            .map(|(k, v)| (k.clone(), v.into()))
            .collect();
        SparseTensor {
            elements,
            structure: self.structure.clone(),
        }
    }
}
#[enum_dispatch]
pub trait HasTensorData<T> {
    fn data(&self) -> Vec<T>;

    fn indices(&self) -> Vec<Vec<ConcreteIndex>>;

    fn hashmap(&self) -> AHashMap<Vec<ConcreteIndex>, T>;

    fn symhashmap(&self, id: Identifier, state: &mut State, ws: &Workspace) -> HashMap<Atom, T>;
}

impl<T, I> HasTensorData<T> for DenseTensor<T, I>
where
    T: Clone,
{
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

    fn hashmap(&self) -> AHashMap<Vec<ConcreteIndex>, T> {
        let mut hashmap = AHashMap::new();
        for (k, v) in self.iter() {
            hashmap.insert(k.clone(), v.clone());
        }
        hashmap
    }

    fn symhashmap(&self, id: Identifier, state: &mut State, ws: &Workspace) -> HashMap<Atom, T> {
        let mut hashmap = HashMap::new();

        for (k, v) in self.iter() {
            hashmap.insert(self.atomic_expanded_label_id(&k, id, state, ws), v.clone());
        }
        hashmap
    }
}

impl<T, I> HasTensorData<T> for SparseTensor<T, I>
where
    T: Clone,
{
    fn data(&self) -> Vec<T> {
        self.elements.values().cloned().collect()
    }

    fn indices(&self) -> Vec<Vec<ConcreteIndex>> {
        self.elements.keys().cloned().collect()
    }

    fn hashmap(&self) -> AHashMap<Vec<ConcreteIndex>, T> {
        self.elements.clone()
    }

    fn symhashmap(&self, id: Identifier, state: &mut State, ws: &Workspace) -> HashMap<Atom, T> {
        let mut hashmap = HashMap::new();

        for (k, v) in self.elements.iter() {
            hashmap.insert(self.atomic_expanded_label_id(&k, id, state, ws), v.clone());
        }
        hashmap
    }
}

#[derive(Debug, Clone, EnumTryAsInner)]
pub enum NumTensor<U, T> {
    Dense(DenseTensor<U, T>),
    Sparse(SparseTensor<U, T>),
}

impl<T, I> HasTensorStructure for NumTensor<T, I>
where
    I: Clone,
{
    type Name = I;
    fn structure(&self) -> &TensorSkeleton<I> {
        match self {
            NumTensor::Dense(d) => d.structure(),
            NumTensor::Sparse(s) => s.structure(),
        }
    }

    fn mut_structure(&mut self) -> &mut TensorSkeleton<I> {
        match self {
            NumTensor::Dense(d) => d.mut_structure(),
            NumTensor::Sparse(s) => s.mut_structure(),
        }
    }
}

impl<T, I> HasTensorData<T> for NumTensor<T, I>
where
    T: Clone,
{
    fn data(&self) -> Vec<T> {
        match self {
            NumTensor::Dense(d) => d.data(),
            NumTensor::Sparse(s) => s.data(),
        }
    }

    fn indices(&self) -> Vec<Vec<ConcreteIndex>> {
        match self {
            NumTensor::Dense(d) => d.indices(),
            NumTensor::Sparse(s) => s.indices(),
        }
    }

    fn hashmap(&self) -> AHashMap<Vec<ConcreteIndex>, T> {
        match self {
            NumTensor::Dense(d) => d.hashmap(),
            NumTensor::Sparse(s) => s.hashmap(),
        }
    }

    fn symhashmap(&self, name: Identifier, state: &mut State, ws: &Workspace) -> HashMap<Atom, T> {
        match self {
            NumTensor::Dense(d) => d.symhashmap(name, state, ws),
            NumTensor::Sparse(s) => s.symhashmap(name, state, ws),
        }
    }
}

#[derive(Debug, Clone, EnumTryAsInner)]
#[derive_err(Debug)]
pub enum NumTensors<T = String> {
    Float(NumTensor<f64, T>),
    Complex(NumTensor<Complex<f64>, T>),
}

impl<T> HasTensorStructure for NumTensors<T>
where
    T: Clone,
{
    type Name = T;
    fn structure(&self) -> &TensorSkeleton<T> {
        match self {
            NumTensors::Float(f) => f.structure(),
            NumTensors::Complex(c) => c.structure(),
        }
    }

    fn mut_structure(&mut self) -> &mut TensorSkeleton<T> {
        match self {
            NumTensors::Float(f) => f.mut_structure(),
            NumTensors::Complex(c) => c.mut_structure(),
        }
    }
}

impl<T> From<DenseTensor<f64, T>> for NumTensors<T> {
    fn from(other: DenseTensor<f64, T>) -> Self {
        NumTensors::Float(NumTensor::Dense(other))
    }
}

impl<T> From<SparseTensor<f64, T>> for NumTensors<T> {
    fn from(other: SparseTensor<f64, T>) -> Self {
        NumTensors::Float(NumTensor::Sparse(other))
    }
}

impl<T> From<DenseTensor<Complex<f64>, T>> for NumTensors<T> {
    fn from(other: DenseTensor<Complex<f64>, T>) -> Self {
        NumTensors::Complex(NumTensor::Dense(other))
    }
}

impl<T> From<SparseTensor<Complex<f64>, T>> for NumTensors<T> {
    fn from(other: SparseTensor<Complex<f64>, T>) -> Self {
        NumTensors::Complex(NumTensor::Sparse(other))
    }
}

// trait UpcastableTo<T> {}

// impl UpcastableTo<Complex<f64>> for f64 {}

// impl<T> NumTensor<T> {
//     pub fn contract<U>(&self, other: &NumTensor<U>) -> Option<Self>
//     where
//         U: UpcastableTo<T>,
//     {
//         match (self, other) {
//             (NumTensor::Dense(s), NumTensor::Dense(o)) => {
//                 s.contract_with_dense(o).map(NumTensor::Dense)
//             }
//             (NumTensor::Dense(s), NumTensor::Sparse(o)) => {
//                 s.contract_with_sparse(o).map(NumTensor::Dense)
//             }
//             (NumTensor::Sparse(s), NumTensor::Dense(o)) => {
//                 s.contract_with_dense(o).map(NumTensor::Dense)
//             }
//             (NumTensor::Sparse(s), NumTensor::Sparse(o)) => {
//                 s.contract_with_sparse(o).map(NumTensor::Sparse)
//             }
//         }
//     }
// }
