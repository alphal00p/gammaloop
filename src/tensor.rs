use enum_dispatch::enum_dispatch;
use num::traits::Num;
use std::{
    borrow::Cow,
    collections::BTreeMap,
    iter::FromIterator,
    ops::{Add, Bound::Included, Mul, Sub},
};
pub mod tensor_structure;
pub use tensor_structure::*;

#[derive(Debug, Clone)]
pub struct SparseTensor<T> {
    elements: BTreeMap<Vec<usize>, T>,
    structure: Vec<Slot>,
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

    pub fn set(&mut self, indices: &[usize], value: T) -> Result<(), String> {
        self.verify_indices(indices)?;
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
    pub fn get(&self, indices: &[usize]) -> Result<Cow<T>, String> {
        self.verify_indices(indices)?;
        // if the index is in the bTree return the value, else return default, lazily allocating the default
        Ok(match self.elements.get(indices) {
            Some(value) => Cow::Borrowed(value),
            None => Cow::Owned(T::default()),
        })
    }
}

pub mod tensor_iterator;

#[derive(Debug, Clone)]
pub struct DenseTensor<T> {
    data: Vec<T>,
    structure: TensorStructure,
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

impl<T> Add<DenseTensor<T>> for DenseTensor<T>
where
    T: Num + Clone + Copy,
{
    type Output = DenseTensor<T>;
    fn add(self, other: DenseTensor<T>) -> DenseTensor<T> {
        assert!(self.structure().same_content(other.structure()));

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let current_value = self.get(&indices).unwrap();
            result.set(&indices, *current_value + *value);
        }
        result
    }
}

impl<T> Add<SparseTensor<T>> for SparseTensor<T>
where
    T: Num + Clone + Default + Copy,
{
    type Output = SparseTensor<T>;
    fn add(self, other: SparseTensor<T>) -> SparseTensor<T> {
        assert!(self.structure().same_content(other.structure()));

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let current_value = self.get(indices).unwrap().into_owned();
            result.set(indices, current_value + *value).unwrap();
        }

        for (indices, value) in self.iter() {
            result.elements.entry(indices.clone()).or_insert(*value);
        }
        result
    }
}

impl<T> Add<SparseTensor<T>> for DenseTensor<T>
where
    T: Num + Clone + Default + Copy,
{
    type Output = DenseTensor<T>;
    fn add(self, other: SparseTensor<T>) -> DenseTensor<T> {
        assert!(self.structure().same_content(other.structure()));

        let mut result = self.clone();

        for (indices, value) in self.iter() {
            let current_value = other.get(&indices).unwrap().into_owned();
            result.set(&indices, current_value + *value);
        }
        result
    }
}

impl<T> Add<DenseTensor<T>> for SparseTensor<T>
where
    T: Num + Clone + Default + Copy,
{
    type Output = DenseTensor<T>;
    fn add(self, other: DenseTensor<T>) -> DenseTensor<T> {
        other + self
    }
}

impl<T> Sub<DenseTensor<T>> for DenseTensor<T>
where
    T: Num + Clone + Copy,
{
    type Output = DenseTensor<T>;
    fn sub(self, other: DenseTensor<T>) -> DenseTensor<T> {
        assert!(self.structure().same_content(other.structure()));

        let mut result = self.clone();

        for (indices, value) in other.iter() {
            let current_value = self.get(&indices).unwrap();
            result.set(&indices, *current_value - *value);
        }
        result
    }
}

impl<T> Sub<SparseTensor<T>> for SparseTensor<T>
where
    T: Num + Default + Copy,
{
    type Output = SparseTensor<T>;
    fn sub(self, other: SparseTensor<T>) -> SparseTensor<T> {
        assert!(self.structure().same_content(other.structure()));
        let permutation = self
            .structure()
            .find_permutation(other.structure())
            .unwrap();

        let mut result = self.clone();

        for (indices, value) in self.iter() {
            result.elements.entry(indices.clone()).or_insert(*value);
        }
        for (indices, value) in self.iter() {
            let permuted_indices: Vec<usize> =
                permutation.iter().map(|&index| indices[index]).collect();
            let other_value = other.get(&permuted_indices).unwrap().into_owned();
            result.set(indices, other_value - *value).unwrap();
            if other_value - *value == T::default() {
                result.elements.remove(indices);
            }
        }

        result
    }
}

impl<T> Sub<SparseTensor<T>> for DenseTensor<T>
where
    T: Num + Clone + Default + Copy,
{
    type Output = DenseTensor<T>;
    fn sub(self, other: SparseTensor<T>) -> DenseTensor<T> {
        assert!(self.structure().same_content(other.structure()));

        let mut result = self.clone();

        for (indices, value) in self.iter() {
            let current_value = other.get(&indices).unwrap().into_owned();
            result.set(&indices, current_value - *value);
        }
        result
    }
}

impl<T> Sub<DenseTensor<T>> for SparseTensor<T>
where
    T: Num + Clone + Default + Copy,
{
    type Output = DenseTensor<T>;
    fn sub(self, other: DenseTensor<T>) -> DenseTensor<T> {
        other - self
    }
}

impl<T> Mul<T> for DenseTensor<T>
where
    T: Num + Clone + Copy,
{
    type Output = DenseTensor<T>;
    fn mul(self, other: T) -> DenseTensor<T> {
        let mut result = self.clone();

        for (indices, value) in self.iter() {
            result.set(&indices, other * *value);
        }
        result
    }
}

impl<T> Mul<T> for SparseTensor<T>
where
    T: Num + Copy,
{
    type Output = SparseTensor<T>;
    fn mul(self, other: T) -> Self::Output {
        let mut result = self.clone();

        for (indices, value) in self.iter() {
            result.set(indices, other * *value).unwrap();
        }
        result
    }
}

#[allow(dead_code)]
struct SymbolicTensor {
    structure: TensorStructure,
    expression: String,
}

impl HasTensorStructure for SymbolicTensor {
    fn structure(&self) -> &TensorStructure {
        &self.structure
    }
}

#[allow(dead_code)]
#[enum_dispatch(HasTensorStructure)]
enum Tensor<T> {
    Num(NumTensor<T>),
    Symbolic(SymbolicTensor),
}
pub mod tensor_contraction;
pub use tensor_contraction::*;
pub mod ufo_spin_tensors;

#[cfg(test)]
mod tests;
