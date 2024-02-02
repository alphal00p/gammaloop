use enum_dispatch::enum_dispatch;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::fmt::Debug;
use std::ops::Index;
use std::ops::Range;

use symbolica::representations::{AsAtomView, Atom, FunctionBuilder, Identifier};
use symbolica::state::{State, Workspace};

use std::collections::HashSet;
use std::{cmp::Ordering, collections::HashMap};

use super::DenseTensor;
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

impl std::fmt::Display for Representation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Representation::Euclidean(value) => write!(f, "e{}", value),
            Representation::Lorentz(value) => write!(f, "l{}", value),
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

impl std::fmt::Display for Slot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.index, self.representation)
    }
}

pub type TensorStructure = Vec<Slot>;
#[derive(Clone, PartialEq, Debug)]
pub struct TensorSkeleton<N> {
    internal: TensorStructure,
    external: TensorStructure,
    pub names: HashMap<Range<usize>, N>, //ideally this is a named partion.. maybe a btreemap<usize, N>, and the range is from previous to next
    pub global_name: Option<N>,
}

impl<N> TensorSkeleton<N> {
    /// Since each time we contract, we merge the name maps, the amount of contractions, is the size of the name map
    /// This function returns the number of contractions thus computed
    pub fn contractions_num(&self) -> usize {
        self.names.len()
    }

    /// Given two TensorSkeletons, returns the index of the first matching slot in each external index list
    pub fn match_index(&self, other: &Self) -> Option<(usize, usize)> {
        for (i, slot_a) in self.external.iter().enumerate().rev() {
            for (j, slot_b) in other.external.iter().enumerate() {
                if slot_a == slot_b {
                    return Some((i, j));
                }
            }
        }
        None
    }

    /// Identify the repeated slots in the external index list
    fn traces(&self) -> Vec<[usize; 2]> {
        let mut positions = HashMap::new();

        // Track the positions of each element
        for (index, &value) in self.external.iter().enumerate() {
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

    /// yields the (outwards facing) shape of the tensor as a list of dimensions
    pub fn shape(&self) -> Vec<Dimension> {
        self.external
            .iter()
            .map(|slot| &slot.representation)
            .collect()
    }

    /// yields the order/total valence of the tensor, i.e. the number of indices
    /// (or misnamed : rank)
    pub fn order(&self) -> usize {
        //total valence (or misnamed : rank)
        self.external.len()
    }

    /// checks if externally, the two tensors are the same
    pub fn same_external(&self, other: &Self) -> bool {
        let set1: HashSet<_> = self.external.iter().collect();
        let set2: HashSet<_> = other.external.iter().collect();
        set1 == set2
    }

    /// checks if internally, the two tensors are the same. This implies that the external indices are the same
    pub fn same_content(&self, other: &Self) -> bool {
        let set1: HashSet<_> = self.internal.iter().collect();
        let set2: HashSet<_> = other.internal.iter().collect();
        set1 == set2
        // TODO: check names
    }

    /// find the permutation of the external indices that would make the two tensors the same
    pub fn find_permutation(&self, other: &Self) -> Option<Vec<ConcreteIndex>> {
        if self.external.len() != other.external.len() {
            return None;
        }

        let mut index_map = HashMap::new();
        for (i, item) in other.external.iter().enumerate() {
            index_map.entry(item).or_insert_with(Vec::new).push(i);
        }

        let mut permutation = Vec::with_capacity(self.external.len());
        let mut used_indices = HashSet::new();
        for item in self.external.iter() {
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

    /// yields the strides of the tensor in column major order
    pub fn strides_column_major(&self) -> Vec<usize> {
        let mut strides: Vec<usize> = vec![1; self.order()];

        if self.order() == 0 {
            return strides;
        }

        for i in 0..self.order() - 1 {
            strides[i + 1] = strides[i] * usize::from(self.external[i].representation);
        }

        strides
    }

    /// yields the strides of the tensor in row major order
    pub fn strides_row_major(&self) -> Vec<usize> {
        let mut strides = vec![1; self.order()];
        if self.order() == 0 {
            return strides;
        }

        for i in (0..self.order() - 1).rev() {
            strides[i] = strides[i + 1] * usize::from(self.external[i + 1].representation);
        }

        strides
    }

    /// By default, the strides are row major
    pub fn strides(&self) -> Vec<usize> {
        self.strides_row_major()
    }

    /// Verifies that the list of indices provided are valid for the tensor
    pub fn verify_indices(&self, indices: &[ConcreteIndex]) -> Result<(), String> {
        if indices.len() != self.order() {
            return Err("Mismatched order".to_string());
        }

        for (i, &dim_len) in self
            .external
            .iter()
            .map(|slot| &slot.representation)
            .enumerate()
        {
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

    /// yields the flat index of the tensor given a list of indices
    pub fn flat_index(&self, indices: &[ConcreteIndex]) -> Result<usize, String> {
        let strides = self.strides();
        self.verify_indices(indices)?;

        let mut idx = 0;
        for (i, &index) in indices.iter().enumerate() {
            idx += index * strides[i];
        }
        Ok(idx)
    }

    /// yields the expanded index of the tensor given a flat index
    pub fn expanded_index(&self, flat_index: usize) -> Result<Vec<ConcreteIndex>, String> {
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

    /// yields the size of the tensor, i.e. the product of the dimensions. This is the length of the vector of the data in a dense tensor
    pub fn size(&self) -> usize {
        self.shape().iter().product()
    }

    /// yields an iterator over the indices of the tensor
    pub fn index_iter(&self) -> TensorStructureIndexIterator {
        TensorStructureIndexIterator::new(&self.external)
    }

    /// if the tensor has no (external) indices, it is a scalar
    pub fn is_scalar(&self) -> bool {
        self.order() == 0
    }

    /// get the metric along the i-th index
    pub fn get_ith_metric(&self, i: usize) -> Option<Vec<bool>> {
        Some(self.external.get(i)?.representation.negative())
    }

    /// make the indices in the internal index list of self independent from the indices in the internal index list of other
    /// This is done by shifting the indices in the internal index list of self by the the maximum index present.
    pub fn independentize_internal(&mut self, other: &Self) {
        let internal_set: HashSet<Slot> = self
            .internal
            .clone()
            .into_iter()
            .filter(|s| self.external.contains(s))
            .collect();

        let other_set: HashSet<Slot> = other.internal.clone().into_iter().collect();

        let mut replacement_value = internal_set
            .union(&other_set)
            .map(|s| s.index)
            .max()
            .unwrap_or(0)
            + 1;

        for item in self.internal.iter_mut() {
            if other_set.contains(item) {
                item.index = replacement_value;
                replacement_value += 1;
            }
        }
    }

    /// remove the repeated indices in the external index list
    fn trace_out(&mut self) {
        self.external = self.external.clone().into_iter().unique().collect();
    }

    /// remove the given indices from the external index list
    pub fn trace(&mut self, i: usize, j: usize) {
        if i < j {
            self.trace(j, i);
            return;
        }
        let a = self.external.remove(i);
        let b = self.external.remove(j);
        assert_eq!(a, b);
    }

    /// If the tensor has internal structure, some contraction or trace has been performed. It is then a composite tensor, and the outwards facing indices are not the same as the internal indices
    pub fn is_composite(&self) -> bool {
        self.internal.len() != self.external.len()
    }

    /// get the global name of the tensor
    pub fn global_name(&self) -> Option<&N> {
        match self.global_name {
            Some(ref name) => Some(name),
            None => None,
        }
    }

    /// set the global name of the tensor
    pub fn set_global_name(&mut self, name: N) {
        self.global_name = Some(name);
    }

    pub fn atomic_expanded_label<I: IntoId>(
        &self,
        indices: &[ConcreteIndex],
        name: I,
        state: &mut State,
        ws: &Workspace,
    ) -> Atom {
        let id = name.into_id(state);
        self.atomic_expanded_label_id(indices, id, state, ws)
    }

    pub fn atomic_flat_label<I: IntoId>(
        &self,
        index: usize,
        name: I,
        state: &mut State,
        ws: &Workspace,
    ) -> Atom {
        let id = name.into_id(state);
        self.atomic_flat_label_id(index, id, state, ws)
    }

    pub fn atomic_flat_label_id(
        &self,
        index: usize,
        id: Identifier,
        state: &mut State,
        ws: &Workspace,
    ) -> Atom {
        let mut value_builder = FunctionBuilder::new(id, state, ws);
        value_builder = value_builder.add_arg(Atom::new_num(index as i64).as_atom_view());
        value_builder.finish().into_atom()
    }

    pub fn atomic_expanded_label_id(
        &self,
        indices: &[ConcreteIndex],
        id: Identifier,
        state: &mut State,
        ws: &Workspace,
    ) -> Atom {
        let mut value_builder = FunctionBuilder::new(id, state, ws);
        for &index in indices {
            value_builder = value_builder.add_arg(Atom::new_num(index as i64).as_atom_view());
        }
        value_builder.finish().into_atom()
    }
}

impl<N> TensorSkeleton<N>
where
    N: Clone,
{
    /// Constructs a new TensorSkeleton from a list of tuples of indices and dimension (assumes they are all euclidean), along with a name
    pub fn from_integers(slots: &[(AbstractIndex, Dimension)], name: N) -> Self {
        let slots: Vec<(AbstractIndex, Representation)> = slots
            .iter()
            .map(|(index, dim)| (*index, Representation::Euclidean(*dim)))
            .collect();
        Self::from_idxsing(&slots, name)
    }
    /// Constructs a new TensorSkeleton from a list of tuples of indices and representations, along with a name
    pub fn from_idxsing(slots: &[(AbstractIndex, Representation)], name: N) -> Self {
        let structure: Vec<Slot> = slots
            .iter()
            .map(|(index, representation)| Slot::from((*index, *representation)))
            .collect();

        let name_map = HashMap::from([(0..structure.len(), name.clone())]);

        TensorSkeleton {
            internal: structure.clone(),
            external: structure,
            names: name_map,
            global_name: Some(name),
        }
    }
    /// essentially contract.
    pub fn merge(&mut self, other: &Self) {
        let shift = self.internal.len();
        for (range, name) in other.names.iter() {
            self.names
                .insert((range.start + shift)..(range.end + shift), name.clone());
        }
        self.external.append(&mut other.external.clone());
        self.trace_out();
        self.independentize_internal(other);
        self.internal.append(&mut other.internal.clone());
    }

    /// Merge two TensorSkeletons at the given positions of the external index list. Ideally the internal index list should be independentized before merging
    /// This is essentially a contraction of only one index. The name maps are merged, and shifted accordingly. The global name is lost, since the resulting tensor is composite
    /// The global name can be set again with the set_global_name function
    pub fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self {
        let mut slots_b = other.external.clone();
        let mut slots_a = self.external.clone();

        slots_a.remove(positions.0);
        slots_b.remove(positions.1);

        let mut slots_a_int = self.internal.clone();
        let mut slots_b_int = other.internal.clone();
        slots_a_int.append(&mut slots_b_int);

        let mut names = self.names.clone();
        let shift = self.internal.len();
        for (range, name) in other.names.iter() {
            names.insert((range.start + shift)..(range.end + shift), name.clone());
        }
        slots_a.append(&mut slots_b);
        TensorSkeleton {
            internal: slots_a_int,
            external: slots_a,
            names,
            global_name: None,
        }
    }
}

pub trait IntoId {
    fn into_id(self, state: &mut State) -> Identifier;
}

impl IntoId for String {
    fn into_id(self, state: &mut State) -> Identifier {
        state.get_or_insert_fn(self, None).unwrap()
    }
}

impl IntoId for Identifier {
    fn into_id(self, _state: &mut State) -> Identifier {
        self
    }
}

impl IntoId for &str {
    fn into_id(self, state: &mut State) -> Identifier {
        state.get_or_insert_fn(self, None).unwrap()
    }
}

impl<N> TensorSkeleton<N>
where
    N: IntoId + Clone,
{
    /// turn the structure into a dense tensor made up of atoms, with the global name as the name of labels
    pub fn to_dense(self, state: &mut State, ws: &Workspace) -> Option<DenseTensor<Atom, N>> {
        let name = self.global_name()?.clone();
        Some(self.to_dense_shadow(state, ws, name))
    }

    pub fn dense_shadow(
        &mut self,
        state: &mut State,
        ws: &Workspace,
        name: N,
    ) -> DenseTensor<Atom, N> {
        self.set_global_name(name.clone());
        self.clone().to_dense_shadow(state, ws, name)
    }

    pub fn to_dense_shadow(
        mut self,
        state: &mut State,
        ws: &Workspace,
        name: N,
    ) -> DenseTensor<Atom, N> {
        let f_id = name.clone().into_id(state);
        let mut data = vec![];
        for index in self.index_iter() {
            data.push(self.atomic_expanded_label_id(&index, f_id, state, ws));
        }

        self.internal = self.external.clone();
        self.names = HashMap::new();
        self.names.insert(0..self.internal.len(), name);
        DenseTensor {
            data,
            structure: self,
        }
    }
}
impl TensorSkeleton<String> {
    pub fn to_symbolic(&self, label: Identifier, ws: &Workspace, state: &mut State) -> Atom {
        let atoms = self
            .external
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

impl std::fmt::Display for TensorSkeleton<String> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let mut string = String::new();
        if let Some(global_name) = self.global_name() {
            string.push_str(&format!("{}:", global_name));
        }
        for (range, name) in self
            .names
            .iter()
            .filter(|(r, _)| *r != &(0..self.internal.len()) || !self.is_composite())
        {
            string.push_str(&format!("{}(", name));
            for slot in self.internal[range.clone()].iter() {
                string.push_str(&format!("{},", slot));
            }
            string.pop();
            string.push(')');
        }
        write!(f, "{}", string)
    }
}

impl TensorSkeleton<Identifier> {
    pub fn to_string(&self, state: &State) -> String {
        let mut string = String::new();
        if let Some(global_name) = self.global_name() {
            string.push_str(&format!("{}:", state.get_name(*global_name)));
        }
        for (range, name) in self
            .names
            .iter()
            .filter(|(r, _)| *r != &(0..self.internal.len()) || !self.is_composite())
        {
            string.push_str(&format!("{}(", state.get_name(*name)));
            for slot in self.internal[range.clone()].iter() {
                string.push_str(&format!("{},", slot));
            }
            string.pop();
            string.push(')');
        }
        string
    }
}

impl<T> Index<usize> for TensorSkeleton<T> {
    type Output = Slot;

    fn index(&self, index: usize) -> &Self::Output {
        &self.external[index]
    }
}

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
    type Name;
    fn mut_structure(&mut self) -> &mut TensorSkeleton<Self::Name>;
    fn structure(&self) -> &TensorSkeleton<Self::Name>;
    // inline
    fn order(&self) -> usize {
        self.structure().order()
    }

    fn contractions_num(&self) -> usize {
        self.structure().contractions_num()
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

    fn traces(&self) -> Vec<[usize; 2]> {
        self.structure().traces()
    }

    fn is_scalar(&self) -> bool {
        self.structure().is_scalar()
    }

    fn get_ith_metric(&self, i: usize) -> Option<Vec<bool>> {
        self.structure().get_ith_metric(i)
    }

    fn is_composite(&self) -> bool {
        self.structure().is_composite()
    }

    fn global_name(&self) -> Option<&Self::Name> {
        self.structure().global_name()
    }

    fn atomic_expanded_label<I: IntoId>(
        &self,
        indices: &[ConcreteIndex],
        name: I,
        state: &mut State,
        ws: &Workspace,
    ) -> Atom {
        self.structure()
            .atomic_expanded_label(indices, name, state, ws)
    }

    fn atomic_flat_label<I: IntoId>(
        &self,
        index: usize,
        name: I,
        state: &mut State,
        ws: &Workspace,
    ) -> Atom {
        self.structure().atomic_flat_label(index, name, state, ws)
    }

    fn atomic_flat_label_id(
        &self,
        index: usize,
        id: Identifier,
        state: &mut State,
        ws: &Workspace,
    ) -> Atom {
        self.structure().atomic_flat_label_id(index, id, state, ws)
    }

    fn atomic_expanded_label_id(
        &self,
        indices: &[ConcreteIndex],
        id: Identifier,
        state: &mut State,
        ws: &Workspace,
    ) -> Atom {
        self.structure()
            .atomic_expanded_label_id(indices, id, state, ws)
    }
}
