use enum_dispatch::enum_dispatch;

use std::collections::HashSet;
use std::{cmp::Ordering, collections::HashMap};

pub type AbstractIndex = usize;
pub type Dimension = usize;
pub type ConcreteIndex = usize;
pub type Position = usize;

#[derive(PartialEq, Eq, Clone, Copy, Debug, Hash)]
pub enum Representation {
    Euclidean(Dimension),
    Lorentz(Dimension),
}

impl PartialOrd for Representation {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Representation {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Representation::Euclidean(dim1), Representation::Euclidean(dim2))
            | (Representation::Lorentz(dim1), Representation::Lorentz(dim2)) => dim1.cmp(dim2),
            (Representation::Euclidean(_), Representation::Lorentz(_)) => Ordering::Less,
            (Representation::Lorentz(_), Representation::Euclidean(_)) => Ordering::Greater,
        }
    }
}

impl Representation {
    pub fn negative(&self) -> Vec<bool> {
        match self {
            Representation::Euclidean(value) => vec![false; *value],
            Representation::Lorentz(value) => std::iter::once(false)
                .chain(std::iter::repeat(true).take(*value - 1))
                .collect::<Vec<_>>(),
        }
    }
}

impl From<Dimension> for Representation {
    fn from(value: Dimension) -> Self {
        Representation::Euclidean(value)
    }
}

impl<'a> std::iter::FromIterator<&'a Representation> for Vec<Dimension> {
    fn from_iter<T: IntoIterator<Item = &'a Representation>>(iter: T) -> Self {
        iter.into_iter()
            .map(|&rep| -> Dimension { (&rep).into() })
            .collect()
    }
}

impl From<&Representation> for Dimension {
    fn from(rep: &Representation) -> Self {
        match rep {
            Representation::Euclidean(value) => *value,
            Representation::Lorentz(value) => *value,
        }
    }
}

impl From<Representation> for Dimension {
    fn from(rep: Representation) -> Self {
        match rep {
            Representation::Euclidean(value) => value,
            Representation::Lorentz(value) => value,
        }
    }
}

#[derive(Debug, Copy, Clone, Hash, PartialEq, Eq)]
pub struct Slot {
    index: AbstractIndex,
    pub representation: Representation,
}

impl PartialOrd for Slot {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Slot {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.representation.cmp(&other.representation) {
            Ordering::Equal => self.index.cmp(&other.index),
            other => other,
        }
    }
}

impl From<(AbstractIndex, Representation)> for Slot {
    fn from(value: (AbstractIndex, Representation)) -> Self {
        Slot {
            index: value.0,
            representation: value.1,
        }
    }
}

pub type TensorStructure = Vec<Slot>;

pub trait VecSlotExtension {
    fn from_idxsing(indices: &[AbstractIndex], dims: &[Representation]) -> Self;
    fn from_integers(indices: &[AbstractIndex], dims: &[usize]) -> Self;
    fn match_index(&self, other: &Self) -> Option<(Position, Position)>;
    fn traces(&self) -> Vec<[Position; 2]>;
    fn merge_at(&self, other: &Self, positions: (Position, Position)) -> Self;
    fn shape(&self) -> Vec<Dimension>;
    fn order(&self) -> usize;
    fn same_content(&self, other: &Self) -> bool;
    fn find_permutation(&self, other: &Self) -> Option<Vec<usize>>;
    fn strides_row_major(&self) -> Vec<usize>;
    fn strides_column_major(&self) -> Vec<usize>;
    fn strides(&self) -> Vec<usize>;
    fn verify_indices(&self, indices: &[ConcreteIndex]) -> Result<(), String>;
    fn flat_index(&self, indices: &[ConcreteIndex]) -> Result<usize, String>;
    fn expanded_index(&self, flat_index: usize) -> Result<Vec<ConcreteIndex>, String>;
    fn size(&self) -> usize;
}

impl VecSlotExtension for TensorStructure {
    fn from_idxsing(indices: &[AbstractIndex], signatures: &[Representation]) -> Self {
        indices
            .iter()
            .zip(signatures.iter())
            .map(|(&index, &dim)| Slot::from((index, dim)))
            .collect()
    }

    fn same_content(&self, other: &Self) -> bool {
        let set1: HashSet<_> = self.iter().collect();
        let set2: HashSet<_> = other.iter().collect();
        set1 == set2
    }

    fn find_permutation(&self, other: &Self) -> Option<Vec<usize>> {
        if self.len() != other.len() {
            return None;
        }

        let mut index_map = HashMap::new();
        for (i, item) in other.iter().enumerate() {
            index_map.entry(item).or_insert_with(Vec::new).push(i);
        }

        let mut permutation = Vec::with_capacity(self.len());
        let mut used_indices = HashSet::new();
        for item in self {
            if let Some(indices) = index_map.get_mut(item) {
                // Find an index that hasn't been used yet
                if let Some(&index) = indices.iter().find(|&&i| !used_indices.contains(&i)) {
                    permutation.push(index);
                    used_indices.insert(index);
                } else {
                    // No available index for this item
                    return None;
                }
            } else {
                // Item not found in other
                return None;
            }
        }

        Some(permutation)
    }

    fn from_integers(indices: &[AbstractIndex], dims: &[usize]) -> Self {
        indices
            .iter()
            .zip(dims.iter())
            .map(|(&index, &dim)| Slot::from((index, Representation::Euclidean(dim))))
            .collect()
    }

    fn match_index(&self, other: &Self) -> Option<(Position, Position)> {
        for (i, slot_a) in self.iter().enumerate().rev() {
            for (j, slot_b) in other.iter().enumerate() {
                if slot_a == slot_b {
                    return Some((i, j));
                }
            }
        }
        None
    }

    fn traces(&self) -> Vec<[Position; 2]> {
        let mut positions = HashMap::new();

        // Track the positions of each element
        for (index, &value) in self.iter().enumerate() {
            positions.entry(value).or_insert_with(Vec::new).push(index);
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

    fn merge_at(&self, other: &Self, positions: (Position, Position)) -> Self {
        let mut slots_b = other.clone();
        let mut slots_a = self.clone();

        slots_a.remove(positions.0);
        slots_b.remove(positions.1);

        slots_a.append(&mut slots_b);
        slots_a
    }

    fn shape(&self) -> Vec<Dimension> {
        self.iter().map(|slot| &slot.representation).collect()
    }

    fn order(&self) -> usize {
        //total valence (or misnamed : rank)
        self.len()
    }

    fn strides_column_major(&self) -> Vec<usize> {
        let mut strides: Vec<usize> = vec![1; self.order()];

        if self.order() == 0 {
            return strides;
        }

        for i in 0..self.order() - 1 {
            strides[i + 1] = strides[i] * usize::from(self[i].representation);
        }

        strides
    }

    fn strides_row_major(&self) -> Vec<usize> {
        let mut strides = vec![1; self.order()];
        if self.order() == 0 {
            return strides;
        }

        for i in (0..self.order() - 1).rev() {
            strides[i] = strides[i + 1] * usize::from(self[i + 1].representation);
        }

        strides
    }

    fn strides(&self) -> Vec<usize> {
        self.strides_row_major()
    }

    fn verify_indices(&self, indices: &[ConcreteIndex]) -> Result<(), String> {
        if indices.len() != self.order() {
            return Err("Mismatched order".to_string());
        }

        for (i, &dim_len) in self.iter().map(|slot| &slot.representation).enumerate() {
            if indices[i] >= usize::from(dim_len) {
                return Err(format!(
                    "Index {} out of bounds for dimension {} of size {}",
                    indices[i],
                    i,
                    usize::from(dim_len)
                ));
            }
        }
        Ok(())
    }

    fn flat_index(&self, indices: &[usize]) -> Result<usize, String> {
        let strides = self.strides();
        self.verify_indices(indices)?;

        let mut idx = 0;
        for (i, &index) in indices.iter().enumerate() {
            idx += index * strides[i];
        }
        Ok(idx)
    }

    fn expanded_index(&self, flat_index: usize) -> Result<Vec<usize>, String> {
        let mut indices = vec![];
        let mut index = flat_index;
        for &stride in self.strides().iter() {
            indices.push(index / stride);
            index %= stride;
        }
        if index == 0 {
            Ok(indices)
        } else {
            Err(format!("Index {} out of bounds", flat_index))
        }
    }

    fn size(&self) -> usize {
        self.shape().iter().product()
    }
}

#[enum_dispatch]
pub trait HasTensorStructure {
    fn structure(&self) -> &Vec<Slot>;
    // inline
    fn order(&self) -> usize {
        self.structure().order()
    }

    fn shape(&self) -> Vec<usize> {
        self.structure().shape()
    }

    fn size(&self) -> usize {
        self.structure().size()
    }

    fn verify_indices(&self, indices: &[usize]) -> Result<(), String> {
        self.structure().verify_indices(indices)
    }

    fn strides(&self) -> Vec<usize> {
        self.structure().strides()
    }

    fn flat_index(&self, indices: &[usize]) -> Result<usize, String> {
        self.structure().flat_index(indices)
    }

    fn expanded_index(&self, flat_index: usize) -> Result<Vec<usize>, String> {
        self.structure().expanded_index(flat_index)
    }

    fn match_index(&self, other: &dyn HasTensorStructure) -> Option<(usize, usize)> {
        self.structure().match_index(other.structure())
    }

    fn traces(&self) -> Vec<[Position; 2]> {
        self.structure().traces()
    }
}
