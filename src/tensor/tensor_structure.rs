use enum_dispatch::enum_dispatch;
use symbolica::representations::{AsAtomView, Atom, FunctionBuilder, Identifier};
use symbolica::state::{State, Workspace};

use std::collections::HashSet;
use std::{cmp::Ordering, collections::HashMap};

use super::TensorStructureIndexIterator;
/// usize is used as label/id for index of tensor
pub type AbstractIndex = usize;
/// usize is used as a Dimension
pub type Dimension = usize;
/// usize is used as a concrete index, i.e. the concrete usize/index of the corresponding abstract index
pub type ConcreteIndex = usize;

#[derive(PartialEq, Eq, Clone, Copy, Debug, Hash, PartialOrd, Ord)]
pub enum Representation {
    Euclidean(Dimension),
    Lorentz(Dimension),
}

impl Representation {
    #[inline]
    // this could be implemented directly in the fiberiterator.
    pub fn negative(&self) -> Vec<bool> {
        match self {
            Representation::Euclidean(value) => vec![false; *value],
            Representation::Lorentz(value) => std::iter::once(false)
                .chain(std::iter::repeat(true).take(*value - 1))
                .collect::<Vec<_>>(),
        }
    }

    pub fn to_fnbuilder<'a, 'b: 'a>(
        &'a self,
        ws: &'b Workspace,
        state: &'b mut State,
    ) -> FunctionBuilder<'a> {
        let (value, id) = match self {
            Representation::Euclidean(value) => (*value, state.get_or_insert_fn("euc", None)),
            Representation::Lorentz(value) => (*value, state.get_or_insert_fn("lor", None)),
        };

        let mut value_builder = FunctionBuilder::new(id.unwrap(), state, ws);

        value_builder = value_builder.add_arg(Atom::new_num(value as i64).as_atom_view());

        value_builder
    }

    pub fn to_symbolic(&self, ws: &Workspace, state: &mut State) -> Atom {
        self.to_fnbuilder(ws, state).finish().into_atom()
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

impl Slot {
    pub fn to_symbolic(&self, ws: &Workspace, state: &mut State) -> Atom {
        let mut value_builder = self.representation.to_fnbuilder(ws, state);
        value_builder = value_builder.add_arg(Atom::new_num(self.index as i64).as_atom_view());
        value_builder.finish().into_atom()
    }
}

pub type TensorStructure = Vec<Slot>;

pub trait VecSlotExtension {
    fn from_idxsing(indices: &[AbstractIndex], dims: &[Representation]) -> Self;
    fn from_integers(indices: &[AbstractIndex], dims: &[Dimension]) -> Self;
    fn match_index(&self, other: &Self) -> Option<(usize, usize)>;
    fn traces(&self) -> Vec<[usize; 2]>;
    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self;
    fn shape(&self) -> Vec<Dimension>;
    fn order(&self) -> usize;
    fn same_content(&self, other: &Self) -> bool;
    fn find_permutation(&self, other: &Self) -> Option<Vec<ConcreteIndex>>;
    fn strides_row_major(&self) -> Vec<usize>;
    fn strides_column_major(&self) -> Vec<usize>;
    fn strides(&self) -> Vec<usize>;
    fn verify_indices(&self, indices: &[ConcreteIndex]) -> Result<(), String>;
    fn flat_index(&self, indices: &[ConcreteIndex]) -> Result<usize, String>;
    fn expanded_index(&self, flat_index: usize) -> Result<Vec<ConcreteIndex>, String>;
    fn size(&self) -> usize;
    fn index_iter(&self) -> TensorStructureIndexIterator;
    fn to_symbolic(&self, label: Identifier, ws: &Workspace, state: &mut State) -> Atom;
    fn is_scalar(&self) -> bool {
        self.order() == 0
    }
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

    fn find_permutation(&self, other: &Self) -> Option<Vec<ConcreteIndex>> {
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

    fn from_integers(indices: &[AbstractIndex], dims: &[Dimension]) -> Self {
        indices
            .iter()
            .zip(dims.iter())
            .map(|(&index, &dim)| Slot::from((index, Representation::Euclidean(dim))))
            .collect()
    }

    fn match_index(&self, other: &Self) -> Option<(usize, usize)> {
        for (i, slot_a) in self.iter().enumerate().rev() {
            for (j, slot_b) in other.iter().enumerate() {
                if slot_a == slot_b {
                    return Some((i, j));
                }
            }
        }
        None
    }

    fn traces(&self) -> Vec<[usize; 2]> {
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

    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self {
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

    fn flat_index(&self, indices: &[ConcreteIndex]) -> Result<usize, String> {
        let strides = self.strides();
        self.verify_indices(indices)?;

        let mut idx = 0;
        for (i, &index) in indices.iter().enumerate() {
            idx += index * strides[i];
        }
        Ok(idx)
    }

    fn expanded_index(&self, flat_index: usize) -> Result<Vec<ConcreteIndex>, String> {
        let mut indices = vec![];
        let mut index = flat_index;
        for &stride in self.strides().iter() {
            indices.push(index / stride);
            index %= stride;
        }
        if flat_index < self.size() {
            Ok(indices)
        } else {
            Err(format!("Index {} out of bounds", flat_index))
        }
    }

    fn size(&self) -> usize {
        self.shape().iter().product()
    }

    fn index_iter(&self) -> TensorStructureIndexIterator {
        TensorStructureIndexIterator::new(self)
    }

    fn to_symbolic(&self, label: Identifier, ws: &Workspace, state: &mut State) -> Atom {
        let atoms = self
            .iter()
            .map(|slot| slot.to_symbolic(ws, state))
            .collect::<Vec<_>>();

        let mut value_builder = FunctionBuilder::new(label, state, ws);
        for atom in atoms {
            value_builder = value_builder.add_arg(atom.as_atom_view());
        }
        value_builder.finish().into_atom()
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

    fn verify_indices(&self, indices: &[ConcreteIndex]) -> Result<(), String> {
        self.structure().verify_indices(indices)
    }

    fn strides(&self) -> Vec<usize> {
        self.structure().strides()
    }

    fn flat_index(&self, indices: &[ConcreteIndex]) -> Result<usize, String> {
        self.structure().flat_index(indices)
    }

    fn expanded_index(&self, flat_index: usize) -> Result<Vec<ConcreteIndex>, String> {
        self.structure().expanded_index(flat_index)
    }

    fn match_index(&self, other: &dyn HasTensorStructure) -> Option<(usize, usize)> {
        self.structure().match_index(other.structure())
    }

    fn traces(&self) -> Vec<[usize; 2]> {
        self.structure().traces()
    }

    fn is_scalar(&self) -> bool {
        self.structure().is_scalar()
    }
}
