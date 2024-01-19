use enum_dispatch::enum_dispatch;
use num::Complex;
use std::{borrow::Cow, collections::BTreeMap};

use super::{
    AbstractIndex, ConcreteIndex, Dimension, HasTensorStructure, Slot, TensorStructure,
    VecSlotExtension,
};

#[derive(Debug, Clone)]
pub struct SparseTensor<T> {
    pub elements: BTreeMap<Vec<ConcreteIndex>, T>,
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

    pub fn set(&mut self, indices: &[ConcreteIndex], value: T) -> Result<(), String> {
        self.verify_indices(indices)?;
        self.elements.insert(indices.to_vec(), value);
        Ok(())
    }
}

impl<T> SparseTensor<T>
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

impl<T> SparseTensor<T>
where
    T: Clone,
{
    pub fn from_data(
        data: &[(Vec<ConcreteIndex>, T)],
        structure: TensorStructure,
    ) -> Result<Self, String> {
        let mut dimensions = vec![0; structure.len()];
        for (index, _) in data {
            if index.len() != structure.len() {
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
            structure,
        })
    }
}

impl<T> SparseTensor<T>
where
    T: Clone + Default,
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

    pub fn default_from_integers(indices: &[AbstractIndex], dims: &[usize]) -> Self {
        let structure = TensorStructure::from_integers(indices, dims);
        DenseTensor::default(structure)
    }
}

impl<T: Clone> DenseTensor<T> {
    pub fn from_data(data: &[T], structure: TensorStructure) -> Result<Self, String> {
        if data.len() != structure.size() && !(data.len() == 1 && structure.is_scalar()) {
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
        if structure.is_scalar() {
            data.truncate(1);
        } else {
            data.truncate(structure.size());
        }
        Ok(DenseTensor { data, structure })
    }
}

impl<T> DenseTensor<T> {
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

impl<T> DenseTensor<T> {
    pub fn conver_to<U>(&self) -> DenseTensor<U>
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

impl<T> SparseTensor<T> {
    pub fn convert_to<U>(&self) -> SparseTensor<U>
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

#[enum_dispatch(HasTensorStructure)]
pub enum NumTensor<T> {
    Dense(DenseTensor<T>),
    Sparse(SparseTensor<T>),
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
