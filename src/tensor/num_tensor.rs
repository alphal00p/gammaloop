use enum_dispatch::enum_dispatch;
use std::{borrow::Cow, collections::BTreeMap};

use super::{
    AbstractIndex, ConcreteIndex, Dimension, HasTensorStructure, Position, Slot, TensorStructure,
    VecSlotExtension,
};

#[derive(Debug, Clone)]
pub struct SparseTensor<T> {
    pub elements: BTreeMap<Vec<usize>, T>,
    pub structure: Vec<Slot>,
}

impl<T> HasTensorStructure for SparseTensor<T> {
    fn structure(&self) -> &Vec<Slot> {
        &self.structure
    }
}

impl<T> SparseTensor<T> {
    pub fn empty(structure: TensorStructure) -> Self {
        SparseTensor {
            elements: BTreeMap::new(),
            structure,
        }
    }

    pub fn empty_from_integers(indices: &[AbstractIndex], dims: &[Dimension]) -> Self {
        let structure = TensorStructure::from_integers(indices, dims);
        SparseTensor {
            elements: BTreeMap::new(),
            structure,
        }
    }

    pub fn is_empty_at(&self, indices: &[ConcreteIndex]) -> bool {
        !self.elements.contains_key(indices)
    }

    pub fn set(&mut self, indices: &[usize], value: T) -> Result<(), String> {
        self.verify_indices(indices)?;
        self.elements.insert(indices.to_vec(), value);
        Ok(())
    }
}

impl<T> SparseTensor<T>
where
    T: Default + std::cmp::PartialEq,
{
    pub fn smart_set(&mut self, indices: &[usize], value: T) -> Result<(), String> {
        self.verify_indices(indices)?;
        if value == T::default() {
            _ = self.elements.remove(indices);
            return Ok(());
        }
        self.elements.insert(indices.to_vec(), value);
        Ok(())
    }
}

impl<T> SparseTensor<T>
where
    T: Clone,
{
    pub fn from_data(
        data: &[(Vec<ConcreteIndex>, T)],
        indices: &[AbstractIndex],
    ) -> Result<Self, String> {
        let mut dimensions = vec![0; indices.len()];
        for (index, _) in data {
            if index.len() != indices.len() {
                return Err("Mismatched order".to_string());
            }
            for (i, &idx) in index.iter().enumerate() {
                if idx >= dimensions[i] {
                    dimensions[i] = idx + 1;
                }
            }
        }
        Ok(SparseTensor {
            elements: BTreeMap::from_iter(data.iter().cloned()),
            structure: TensorStructure::from_integers(indices, &dimensions),
        })
    }
}

impl<T> SparseTensor<T>
where
    T: Clone + Default,
{
    pub fn get(&self, indices: &[usize]) -> Result<&T, String> {
        self.verify_indices(indices)?;
        self.elements
            .get(indices)
            .ok_or("No elements at that spot".to_string())
    }
    pub fn get_with_defaults(&self, indices: &[usize]) -> Result<Cow<T>, String> {
        self.verify_indices(indices)?;
        // if the index is in the bTree return the value, else return default, lazily allocating the default
        Ok(match self.elements.get(indices) {
            Some(value) => Cow::Borrowed(value),
            None => Cow::Owned(T::default()),
        })
    }
}

#[derive(Debug, Clone)]
pub struct DenseTensor<T> {
    pub data: Vec<T>,
    pub structure: TensorStructure,
}

impl<T> HasTensorStructure for DenseTensor<T> {
    fn structure(&self) -> &TensorStructure {
        &self.structure
    }
}

impl<T: Default + Clone> DenseTensor<T> {
    pub fn default(structure: TensorStructure) -> Self {
        DenseTensor {
            data: vec![T::default(); structure.size()],
            structure,
        }
    }

    pub fn default_from_integers(indices: &[AbstractIndex], dims: &[usize]) -> Self {
        let structure = TensorStructure::from_integers(indices, dims);
        DenseTensor {
            data: vec![T::default(); structure.size()],
            structure,
        }
    }
}

impl<T: Clone> DenseTensor<T> {
    pub fn from_data(data: &[T], structure: TensorStructure) -> Result<Self, String> {
        if data.len() != structure.size() {
            return Err("Data length does not match shape".to_string());
        }
        Ok(DenseTensor {
            data: data.to_vec(),
            structure,
        })
    }

    pub fn from_data_coerced(data: &[T], structure: TensorStructure) -> Result<Self, String> {
        if data.len() < structure.size() {
            return Err("Data length is too small".to_string());
        }
        let mut data = data.to_vec();
        data.truncate(structure.size());
        Ok(DenseTensor { data, structure })
    }
}

impl<T> DenseTensor<T> {
    pub fn set(&mut self, indices: &[usize], value: T) {
        let idx = self.flat_index(indices);
        if let Ok(i) = idx {
            self.data[i] = value;
        }
    }

    pub fn get_linear(&self, index: usize) -> Option<&T> {
        self.data.get(index)
    }

    pub fn get(&self, indices: &[usize]) -> Option<&T> {
        if let Ok(idx) = self.flat_index(indices) {
            Some(&self.data[idx])
        } else {
            None
        }
    }
}

trait NumericTensor<T> {}

impl<T> NumericTensor<T> for DenseTensor<T> {}

impl<T> NumericTensor<T> for SparseTensor<T> {}

#[enum_dispatch(HasTensorStructure)]
pub enum NumTensor<T> {
    Dense(DenseTensor<T>),
    Sparse(SparseTensor<T>),
}

// impl<T> NumTensor<T>
// where
//     T: Default
//         + Clone
//         + std::ops::AddAssign
//         + std::ops::SubAssign
//         + std::ops::Mul<Output = T>
//         + Copy
//         + std::cmp::PartialEq
//         + std::fmt::Debug,
// {
//     pub fn contract(&self, other: &Self) -> Option<Self> {
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