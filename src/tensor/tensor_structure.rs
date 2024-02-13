use ahash::AHashMap;
use duplicate::duplicate;

use indexmap::IndexMap;
use smartstring::LazyCompact;
use smartstring::SmartString;
use std::fmt::Debug;
use std::ops::Index;
use std::ops::Range;

use permutation::Permutation;

use symbolica::representations::{AsAtomView, Atom, FunctionBuilder, Identifier};
use symbolica::state::{State, Workspace};

use std::collections::HashSet;
use std::{cmp::Ordering, collections::HashMap};

use super::DenseTensor;
use super::TensorStructureIndexIterator;
use smartstring::alias::String;
/// usize is used as label/id for index of tensor
pub type AbstractIndex = usize;
/// usize is used as a Dimension
pub type Dimension = usize;
/// usize is used as a concrete index, i.e. the concrete usize/index of the corresponding abstract index
pub type ConcreteIndex = usize;

/// Enum for the Representation/Dimension of the index.
#[derive(PartialEq, Eq, Clone, Copy, Debug, Hash, PartialOrd, Ord)]
pub enum Representation {
    /// Represents a Euclidean space of the given dimension, with metric diag(1,1,1,1,...)
    Euclidean(Dimension),
    /// Represents a Minkowski space of the given dimension, with metric diag(1,-1,-1,-1,...)
    Lorentz(Dimension),
    /// Represents a Spinor space of the given dimension
    Spin(Dimension),
    /// Represents a Color Fundamental space of the given dimension
    ColorFundamental(Dimension),
    /// Represents a Color Anti-Fundamental space of the given dimension
    ColorAntiFundamental(Dimension),
    /// Represents a Color Adjoint space of the given dimension
    ColorAdjoint(Dimension),
    /// Represents a Color Sextet space of the given dimension
    ColorSextet(Dimension),
    /// Represents a Color Anti-Sextet space of the given dimension
    ColorAntiSextet(Dimension),
}

impl Representation {
    #[inline]
    // this could be implemented directly in the fiberiterator.
    /// gives the vector of booleans, saying which concrete index along a Dimension/Abstract Index should have a minus sign during contraction.
    ///
    /// # Example
    /// ```
    /// # use _gammaloop::tensor::Representation;
    /// let spin = Representation::Spin(5);
    ///
    /// let metric_diag = spin.negative();
    ///
    /// let mut agree= true;
    ///
    ///     for (i,r) in metric_diag.iter().enumerate(){
    ///         if (r ^ spin.is_neg(i)) {
    ///             agree = false;
    ///         }
    ///     }
    ///
    /// assert!(agree);
    /// ```
    pub fn negative(&self) -> Vec<bool> {
        match self {
            Representation::Lorentz(value) => std::iter::once(false)
                .chain(std::iter::repeat(true).take(*value - 1))
                .collect::<Vec<_>>(),
            Representation::Euclidean(value)
            | Representation::Spin(value)
            | Representation::ColorAdjoint(value)
            | Representation::ColorFundamental(value)
            | Representation::ColorAntiFundamental(value)
            | Representation::ColorSextet(value)
            | Representation::ColorAntiSextet(value) => vec![false; *value],
        }
    }

    /// for the given concrete index, says whether it should have a minus sign during contraction
    ///
    /// for example see [`Self::negative`]
    #[inline]
    pub fn is_neg(&self, i: usize) -> bool {
        match self {
            Representation::Lorentz(_) => i > 0,
            _ => false,
        }
    }

    /// yields a function builder for the representation, adding a first variable: the dimension.
    ///
    /// for example see [Slot::to_symbolic]
    pub fn to_fnbuilder<'a, 'b: 'a>(
        &'a self,
        state: &'b mut State,
        ws: &'b Workspace,
    ) -> FunctionBuilder<'a> {
        let (value, id) = match self {
            Representation::Euclidean(value) => (*value, state.get_or_insert_fn("euc", None)),
            Representation::Lorentz(value) => (*value, state.get_or_insert_fn("lor", None)),
            Representation::Spin(value) => (*value, state.get_or_insert_fn("spin", None)),
            Representation::ColorAdjoint(value) => (*value, state.get_or_insert_fn("CAdj", None)),
            Representation::ColorFundamental(value) => (*value, state.get_or_insert_fn("CF", None)),
            Representation::ColorAntiFundamental(value) => {
                (*value, state.get_or_insert_fn("CAF", None))
            }
            Representation::ColorSextet(value) => (*value, state.get_or_insert_fn("CS", None)),
            Representation::ColorAntiSextet(value) => (*value, state.get_or_insert_fn("CAS", None)),
        };

        let mut value_builder = FunctionBuilder::new(id.unwrap(), state, ws);

        value_builder = value_builder.add_arg(Atom::new_num(value as i64).as_atom_view());

        value_builder
    }

    /// Finishes the function builder into an Atom
    ///
    /// # Example
    ///
    /// ```
    /// # use symbolica::state::{State, Workspace};
    /// # use _gammaloop::tensor::Representation;
    /// # let mut state = State::new();
    /// # let ws = Workspace::new();
    /// let mink = Representation::Lorentz(4);
    ///
    /// assert_eq!("lor(4)",format!("{}",mink.to_symbolic(&mut state,&ws).printer(&state)));
    /// assert_eq!("l4",format!("{}",mink));
    /// ```
    pub fn to_symbolic(&self, state: &mut State, ws: &Workspace) -> Atom {
        self.to_fnbuilder(state, ws).finish().into_atom()
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
            Representation::Euclidean(value)
            | Representation::Lorentz(value)
            | Representation::Spin(value)
            | Representation::ColorAdjoint(value)
            | Representation::ColorFundamental(value)
            | Representation::ColorAntiFundamental(value)
            | Representation::ColorSextet(value)
            | Representation::ColorAntiSextet(value) => *value,
        }
    }
}

impl From<Representation> for Dimension {
    fn from(rep: Representation) -> Self {
        match rep {
            Representation::Euclidean(value)
            | Representation::Lorentz(value)
            | Representation::Spin(value)
            | Representation::ColorAdjoint(value)
            | Representation::ColorFundamental(value)
            | Representation::ColorAntiFundamental(value)
            | Representation::ColorSextet(value)
            | Representation::ColorAntiSextet(value) => value,
        }
    }
}

impl std::fmt::Display for Representation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Representation::Euclidean(value) => write!(f, "e{}", value),
            Representation::Lorentz(value) => write!(f, "l{}", value),
            Representation::Spin(value) => write!(f, "s{}", value),
            Representation::ColorAdjoint(value) => write!(f, "cad{}", value),
            Representation::ColorFundamental(value) => write!(f, "cf{}", value),
            Representation::ColorAntiFundamental(value) => write!(f, "caf{}", value),
            Representation::ColorSextet(value) => write!(f, "cs{}", value),
            Representation::ColorAntiSextet(value) => write!(f, "cas{}", value),
        }
    }
}

/// A [`Slot`] is an index, identified by a `usize` and a [`Representation`].
///
/// A vector of slots thus identifies the shape and type of the tensor.
/// Two indices are considered matching if *both* the `usize` and the [`Representation`] matches.
///
/// # Example
///
/// It can be built from a tuple of `usize` and `Representation`
/// ```
/// # use _gammaloop::tensor::{Representation,Slot};
/// let mink = Representation::Lorentz(4);
/// let mu = Slot::from((0,mink));
///
/// assert_eq!("0l4",format!("{}",mu));
/// ```
///
/// It can also be built from a tuple of `usize` and `usize`, where we default to `Representation::Euclidean`
/// ```
/// # use _gammaloop::tensor::{Representation,Slot};
/// let mu = Slot::from((0,4));
/// assert_eq!("0e4",format!("{}",mu));
/// ```
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

impl From<(usize, usize)> for Slot {
    fn from(value: (usize, usize)) -> Self {
        Slot {
            index: value.0,
            representation: value.1.into(),
        }
    }
}

impl Slot {
    /// using the function builder of the representation add the abstract index as an argument, and finish it to an Atom.
    /// # Example
    ///
    /// ```
    /// # use symbolica::state::{State, Workspace};
    /// # use _gammaloop::tensor::{Representation,Slot};
    /// # let mut state = State::new();
    /// # let ws = Workspace::new();
    /// let mink = Representation::Lorentz(4);
    /// let mu = Slot::from((0,mink));
    ///
    /// assert_eq!("lor(4,0)",format!("{}",mu.to_symbolic(&mut state,&ws).printer(&state)));
    /// assert_eq!("0l4",format!("{}",mu));
    /// ```
    pub fn to_symbolic(&self, state: &mut State, ws: &Workspace) -> Atom {
        let mut value_builder = self.representation.to_fnbuilder(state, ws);
        value_builder = value_builder.add_arg(Atom::new_num(self.index as i64).as_atom_view());
        value_builder.finish().into_atom()
    }
}

impl std::fmt::Display for Slot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.index, self.representation)
    }
}

/// A trait for a any struct that functions as a tensor structure.
/// Only three methods are required to be implemented, the rest are default implementations.
///
/// The associated type `Structure` is the type of the structure. This is usefull for containers of structures, like a datatensor.
/// The two methods `structure` and `mut_structure` are used to get a reference to the structure, and a mutable reference to the structure.
///
pub trait TensorStructure {
    type Structure;
    /// returns the list of slots that are the external indices of the tensor
    fn external_structure(&self) -> &[Slot];

    fn structure(&self) -> &Self::Structure;

    fn mut_structure(&mut self) -> &mut Self::Structure;

    /// checks if the tensor has the same exact structure as another tensor
    fn same_content(&self, other: &Self) -> bool {
        self.same_external(other)
    }

    /// Given two TensorStructures, returns the index of the first matching slot in each external index list, along with a boolean indicating if there is a single match
    fn match_index(&self, other: &Self) -> Option<(bool, usize, usize)> {
        let posmap = AHashMap::from_iter(
            self.external_structure()
                .iter()
                .enumerate()
                .map(|(i, slot)| (slot, i)),
        );

        let mut first_pair: Option<(usize, usize)> = None;

        for (j, slot) in other.external_structure().iter().enumerate() {
            if let Some(&i) = posmap.get(slot) {
                if let Some((i, j)) = first_pair {
                    // Found a second match, return early with false indicating non-unique match
                    return Some((false, i, j));
                }
                first_pair = Some((i, j));
            }
        }

        first_pair.map(|(i, j)| (true, i, j)) // Maps the found pair to Some with true indicating a unique match, or None if no match was found
    }

    /// Given two TensorStructures, returns the index of the first matching slot in each external index list
    fn match_indices(&self, other: &Self) -> Option<(Permutation, Vec<bool>, Vec<bool>)> {
        let mut self_matches = vec![false; self.order()];
        let mut perm = Vec::new();
        let mut other_matches = vec![false; other.order()];

        let posmap = AHashMap::from_iter(
            self.external_structure()
                .iter()
                .enumerate()
                .map(|(i, slot)| (slot, i)),
        );

        for (j, slot_other) in other.external_structure().iter().enumerate() {
            if let Some(&i) = posmap.get(slot_other) {
                self_matches[i] = true;
                other_matches[j] = true;
                perm.push(i);
            }
        }

        if !perm.is_empty() {
            let p: Permutation = permutation::sort(&mut perm);
            Some((p, self_matches, other_matches))
        } else {
            None
        }
    }
    /// Identify the repeated slots in the external index list
    fn traces(&self) -> Vec<[usize; 2]> {
        let mut positions = HashMap::new();

        // Track the positions of each element
        for (index, &value) in self.external_structure().iter().enumerate() {
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
    fn shape(&self) -> Vec<Dimension> {
        self.external_structure()
            .iter()
            .map(|slot| &slot.representation)
            .collect()
    }

    fn reps(&self) -> Vec<Representation> {
        self.external_structure()
            .iter()
            .map(|slot| slot.representation)
            .collect()
    }

    /// yields the order/total valence of the tensor, i.e. the number of indices
    /// (or misnamed : rank)
    fn order(&self) -> usize {
        //total valence (or misnamed : rank)
        self.external_structure().len()
    }

    /// checks if externally, the two tensors are the same
    fn same_external(&self, other: &Self) -> bool {
        let set1: HashSet<_> = self.external_structure().iter().collect();
        let set2: HashSet<_> = other.external_structure().iter().collect();
        set1 == set2
    }

    /// find the permutation of the external indices that would make the two tensors the same
    fn find_permutation(&self, other: &Self) -> Option<Vec<ConcreteIndex>> {
        if self.external_structure().len() != other.external_structure().len() {
            return None;
        }

        let mut index_map = HashMap::new();
        for (i, item) in other.external_structure().iter().enumerate() {
            index_map.entry(item).or_insert_with(Vec::new).push(i);
        }

        let mut permutation = Vec::with_capacity(self.external_structure().len());
        let mut used_indices = HashSet::new();
        for item in self.external_structure().iter() {
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
    fn strides_column_major(&self) -> Vec<usize> {
        let mut strides: Vec<usize> = vec![1; self.order()];

        if self.order() == 0 {
            return strides;
        }

        for i in 0..self.order() - 1 {
            strides[i + 1] = strides[i] * usize::from(self.external_structure()[i].representation);
        }

        strides
    }

    /// yields the strides of the tensor in row major order
    fn strides_row_major(&self) -> Vec<usize> {
        let mut strides = vec![1; self.order()];
        if self.order() == 0 {
            return strides;
        }

        for i in (0..self.order() - 1).rev() {
            strides[i] =
                strides[i + 1] * usize::from(self.external_structure()[i + 1].representation);
        }

        strides
    }

    /// By default, the strides are row major
    fn strides(&self) -> Vec<usize> {
        self.strides_row_major()
    }

    /// Verifies that the list of indices provided are valid for the tensor
    fn verify_indices(&self, indices: &[ConcreteIndex]) -> Result<(), String> {
        if indices.len() != self.order() {
            return Err("Mismatched order".into());
        }

        for (i, &dim_len) in self
            .external_structure()
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
                )
                .into());
            }
        }
        Ok(())
    }

    /// yields the flat index of the tensor given a list of indices
    fn flat_index(&self, indices: &[ConcreteIndex]) -> Result<usize, String> {
        let strides = self.strides();
        self.verify_indices(indices)?;

        let mut idx = 0;
        for (i, &index) in indices.iter().enumerate() {
            idx += index * strides[i];
        }
        Ok(idx)
    }

    /// yields the expanded index of the tensor given a flat index
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
            Err(format!("Index {} out of bounds", flat_index).into())
        }
    }

    /// yields an iterator over the indices of the tensor
    fn index_iter(&self) -> TensorStructureIndexIterator {
        TensorStructureIndexIterator::new(self.external_structure())
    }

    /// if the tensor has no (external) indices, it is a scalar
    fn is_scalar(&self) -> bool {
        self.order() == 0
    }

    /// get the metric along the i-th index
    fn get_ith_metric(&self, i: usize) -> Option<Vec<bool>> {
        Some(self.external_structure().get(i)?.representation.negative())
    }

    /// yields the size of the tensor, i.e. the product of the dimensions. This is the length of the vector of the data in a dense tensor
    fn size(&self) -> usize {
        self.shape().iter().product()
    }
}

impl<'a> TensorStructure for &'a [Slot] {
    type Structure = &'a [Slot];

    fn external_structure(&self) -> &[Slot] {
        self
    }

    fn structure(&self) -> &Self::Structure {
        self
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        self
    }
}

impl TensorStructure for Vec<Slot> {
    type Structure = Vec<Slot>;

    fn structure(&self) -> &Self::Structure {
        self
    }
    fn mut_structure(&mut self) -> &mut Self::Structure {
        self
    }
    fn external_structure(&self) -> &[Slot] {
        self
    }
}

/// A trait for a structure that can be traced and merged, during a contraction.
pub trait StructureContract {
    fn trace(&mut self, i: usize, j: usize);

    fn trace_out(&mut self);

    fn merge(&mut self, other: &Self);

    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self;
}

impl StructureContract for Vec<Slot> {
    fn trace(&mut self, i: usize, j: usize) {
        if i < j {
            self.trace(j, i);
            return;
        }
        let a = self.remove(i);
        let b = self.remove(j);
        assert_eq!(a, b);
    }

    fn trace_out(&mut self) {
        let mut positions = IndexMap::new();

        // Track the positions of each element
        for (index, &value) in self.iter().enumerate() {
            positions.entry(value).or_insert_with(Vec::new).push(index);
        }
        // Collect only the positions of non- repeated elements

        *self = positions
            .into_iter()
            .filter_map(|(value, indices)| {
                if indices.len() == 1 {
                    Some(value)
                } else {
                    None
                }
            })
            .collect();
    }

    fn merge(&mut self, other: &Self) {
        self.append(&mut other.clone());
        self.trace_out();
    }

    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self {
        let mut slots_b = other.clone();
        let mut slots_a = self.clone();

        slots_a.remove(positions.0);
        slots_b.remove(positions.1);

        slots_a.append(&mut slots_b);
        slots_a
    }
}
/// A named structure is a structure with a global name, and a list of slots
///
/// It is useful when you want to shadow tensors, to nest tensor network contraction operations.
#[derive(Clone, PartialEq, Debug)]
pub struct NamedStructure {
    pub structure: Vec<Slot>,
    pub global_name: Option<SmartString<LazyCompact>>,
}

impl NamedStructure {
    /// Constructs a new TensorSkeleton from a list of tuples of indices and dimension (assumes they are all euclidean), along with a name
    pub fn from_integers(slots: &[(AbstractIndex, Dimension)], name: &str) -> Self {
        let slots: Vec<(AbstractIndex, Representation)> = slots
            .iter()
            .map(|(index, dim)| (*index, Representation::Euclidean(*dim)))
            .collect();
        Self::new(&slots, name)
    }
    /// Constructs a new TensorSkeleton from a list of tuples of indices and representations, along with a name
    pub fn new(slots: &[(AbstractIndex, Representation)], name: &str) -> Self {
        let structure: Vec<Slot> = slots
            .iter()
            .map(|(index, representation)| Slot::from((*index, *representation)))
            .collect();

        NamedStructure {
            structure,
            global_name: Some(name.into()),
        }
    }
}

/// A trait for a structure that has a name
pub trait HasName {
    type Name;
    fn name<'a>(&self) -> Option<&Self::Name>;
    fn set_name(&mut self, name: &Self::Name);
}

impl HasName for NamedStructure {
    type Name = SmartString<LazyCompact>;
    fn name(&self) -> Option<&Self::Name> {
        self.global_name.as_ref()
    }
    fn set_name(&mut self, name: &Self::Name) {
        self.global_name = Some(name.clone());
    }
}

impl TensorStructure for NamedStructure {
    type Structure = NamedStructure;
    fn structure(&self) -> &Self::Structure {
        self
    }
    fn mut_structure(&mut self) -> &mut Self::Structure {
        self
    }
    fn external_structure(&self) -> &[Slot] {
        &self.structure
    }
}

impl StructureContract for NamedStructure {
    fn merge(&mut self, other: &Self) {
        self.structure.merge(&other.structure);
    }

    fn trace_out(&mut self) {
        self.structure.trace_out();
    }

    /// when merging two named structures, the global name is lost
    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self {
        NamedStructure {
            structure: self.structure.merge_at(&other.structure, positions),
            global_name: None,
        }
    }

    fn trace(&mut self, i: usize, j: usize) {
        self.structure.trace(i, j);
    }
}

/// A contraction count structure
///
/// Useful for tensor network contraction algorithm.
#[derive(Clone, PartialEq, Debug)]
pub struct ContractionCountStructure {
    pub structure: Vec<Slot>,
    pub contractions: usize,
}

impl ContractionCountStructure {
    /// Constructs a new TensorSkeleton from a list of tuples of indices and dimension (assumes they are all euclidean), along with a name
    pub fn from_integers(slots: &[(AbstractIndex, Dimension)]) -> Self {
        let slots: Vec<(AbstractIndex, Representation)> = slots
            .iter()
            .map(|(index, dim)| (*index, Representation::Euclidean(*dim)))
            .collect();
        Self::new(&slots)
    }
    /// Constructs a new TensorSkeleton from a list of tuples of indices and representations, along with a name
    pub fn new(slots: &[(AbstractIndex, Representation)]) -> Self {
        let structure: Vec<Slot> = slots
            .iter()
            .map(|(index, representation)| Slot::from((*index, *representation)))
            .collect();

        ContractionCountStructure {
            structure,
            contractions: 0,
        }
    }
}

pub trait TracksCount {
    fn contractions_num(&self) -> usize;

    fn is_composite(&self) -> bool {
        self.contractions_num() > 0
    }
}

impl TracksCount for ContractionCountStructure {
    fn contractions_num(&self) -> usize {
        self.contractions
    }
}

impl TensorStructure for ContractionCountStructure {
    type Structure = ContractionCountStructure;
    fn structure(&self) -> &Self::Structure {
        self
    }
    fn mut_structure(&mut self) -> &mut Self::Structure {
        self
    }
    fn external_structure(&self) -> &[Slot] {
        &self.structure
    }
}

impl StructureContract for ContractionCountStructure {
    fn merge(&mut self, other: &Self) {
        self.structure.merge(&other.structure);
        self.contractions += other.contractions;
    }

    fn trace_out(&mut self) {
        self.structure.trace_out();
    }

    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self {
        ContractionCountStructure {
            structure: self.structure.merge_at(&other.structure, positions),
            contractions: self.contractions + other.contractions,
        }
    }

    fn trace(&mut self, i: usize, j: usize) {
        self.structure.trace(i, j);
    }
}

/// A structure to enable smart shadowing of tensors in a tensor network contraction algorithm.
#[derive(Clone, PartialEq, Debug)]
pub struct SmartShadowStructure {
    pub structure: Vec<Slot>,
    pub contractions: usize,
    pub global_name: Option<SmartString<LazyCompact>>,
}

impl SmartShadowStructure {
    /// Constructs a new TensorSkeleton from a list of tuples of indices and dimension (assumes they are all euclidean), along with a name
    pub fn from_integers(slots: &[(AbstractIndex, Dimension)], name: &str) -> Self {
        let slots: Vec<(AbstractIndex, Representation)> = slots
            .iter()
            .map(|(index, dim)| (*index, Representation::Euclidean(*dim)))
            .collect();
        Self::new(&slots, name)
    }
    /// Constructs a new TensorSkeleton from a list of tuples of indices and representations, along with a name
    pub fn new(slots: &[(AbstractIndex, Representation)], name: &str) -> Self {
        let structure: Vec<Slot> = slots
            .iter()
            .map(|(index, representation)| Slot::from((*index, *representation)))
            .collect();

        SmartShadowStructure {
            structure,
            contractions: 0,
            global_name: Some(name.into()),
        }
    }
}

impl HasName for SmartShadowStructure {
    type Name = SmartString<LazyCompact>;
    fn name(&self) -> Option<&SmartString<LazyCompact>> {
        self.global_name.as_ref()
    }
    fn set_name(&mut self, name: &SmartString<LazyCompact>) {
        self.global_name = Some(name.clone());
    }
}

impl TracksCount for SmartShadowStructure {
    fn contractions_num(&self) -> usize {
        self.contractions
    }
}

impl TensorStructure for SmartShadowStructure {
    type Structure = SmartShadowStructure;
    fn structure(&self) -> &Self::Structure {
        self
    }
    fn mut_structure(&mut self) -> &mut Self::Structure {
        self
    }
    fn external_structure(&self) -> &[Slot] {
        &self.structure
    }
}

impl StructureContract for SmartShadowStructure {
    fn merge(&mut self, other: &Self) {
        self.structure.merge(&other.structure);
        self.contractions += other.contractions;
    }

    fn trace_out(&mut self) {
        self.structure.trace_out();
    }

    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self {
        SmartShadowStructure {
            structure: self.structure.merge_at(&other.structure, positions),
            contractions: self.contractions + other.contractions,
            global_name: None,
        }
    }

    fn trace(&mut self, i: usize, j: usize) {
        self.structure.trace(i, j);
    }
}

/// A tracking structure
///
/// It contains two vecs of [`Slot`]s, one for the internal structure, simply extended during each contraction, and one external, coresponding to all the free indices
///
/// It enables keeping track of the contraction history of the tensor, mostly for debugging and display purposes.
/// A [`SymbolicTensor`] can also be used in this way, however it needs a symbolica state and workspace during contraction.
#[derive(Clone, PartialEq, Debug)]
pub struct HistoryStructure<N> {
    internal: Vec<Slot>,
    pub external: Vec<Slot>,
    pub names: AHashMap<Range<usize>, N>, //ideally this is a named partion.. maybe a btreemap<usize, N>, and the range is from previous to next
    pub global_name: Option<N>,
}

impl<N> HistoryStructure<N> {
    /// Constructs a new TensorSkeleton from a list of tuples of indices and dimension (assumes they are all euclidean), along with a name
    pub fn from_integers(slots: &[(AbstractIndex, Dimension)], name: N) -> Self
    where
        N: Clone,
    {
        let slots: Vec<(AbstractIndex, Representation)> = slots
            .iter()
            .map(|(index, dim)| (*index, Representation::Euclidean(*dim)))
            .collect();
        Self::new(&slots, name)
    }
    /// Constructs a new TensorSkeleton from a list of tuples of indices and representations, along with a name
    pub fn new(slots: &[(AbstractIndex, Representation)], name: N) -> Self
    where
        N: Clone,
    {
        let structure: Vec<Slot> = slots
            .iter()
            .map(|(index, representation)| Slot::from((*index, *representation)))
            .collect();

        let name_map = AHashMap::from([(0..structure.len(), name.clone())]);

        HistoryStructure {
            internal: structure.clone(),
            external: structure,
            names: name_map,
            global_name: Some(name),
        }
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
}

impl<N> HasName for HistoryStructure<N>
where
    N: Clone,
{
    type Name = N;
    fn name(&self) -> Option<&N> {
        self.global_name.as_ref()
    }
    fn set_name(&mut self, name: &N) {
        self.global_name = Some(name.clone());
    }
}

impl<N> TracksCount for HistoryStructure<N> {
    /// Since each time we contract, we merge the name maps, the amount of contractions, is the size of the name map
    /// This function returns the number of contractions thus computed
    fn contractions_num(&self) -> usize {
        self.names.len()
    }
}

impl<N> TensorStructure for HistoryStructure<N> {
    type Structure = HistoryStructure<N>;

    fn structure(&self) -> &Self::Structure {
        self
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        self
    }
    fn external_structure(&self) -> &[Slot] {
        &self.external
    }
    /// checks if internally, the two tensors are the same. This implies that the external indices are the same
    fn same_content(&self, other: &Self) -> bool {
        let set1: HashSet<_> = self.internal.iter().collect();
        let set2: HashSet<_> = other.internal.iter().collect();
        set1 == set2
        // TODO: check names
    }
}

// impl TensorStructure for [Slot] {
//     type Structure = [Slot];

//     fn external_structure(&self) -> &[Slot] {
//         self
//     }
// }

impl<N> StructureContract for HistoryStructure<N>
where
    N: Clone,
{
    /// remove the repeated indices in the external index list
    fn trace_out(&mut self) {
        let mut positions = IndexMap::new();

        // Track the positions of each element
        for (index, &value) in self.external.iter().enumerate() {
            positions.entry(value).or_insert_with(Vec::new).push(index);
        }
        // Collect only the positions of non- repeated elements

        self.external = positions
            .into_iter()
            .filter_map(|(value, indices)| {
                if indices.len() == 1 {
                    Some(value)
                } else {
                    None
                }
            })
            .collect();
    }

    /// remove the given indices from the external index list
    fn trace(&mut self, i: usize, j: usize) {
        if i < j {
            self.trace(j, i);
            return;
        }
        let a = self.external.remove(i);
        let b = self.external.remove(j);
        assert_eq!(a, b);
    }

    /// essentially contract.
    fn merge(&mut self, other: &Self) {
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
    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self {
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
        HistoryStructure {
            internal: slots_a_int,
            external: slots_a,
            names,
            global_name: None,
        }
    }
}

pub fn atomic_expanded_label<I: IntoId>(
    indices: &[ConcreteIndex],
    name: I,
    state: &mut State,
    ws: &Workspace,
) -> Atom {
    let id = name.into_id(state);
    atomic_expanded_label_id(indices, id, state, ws)
}

pub fn atomic_flat_label<I: IntoId>(
    index: usize,
    name: I,
    state: &mut State,
    ws: &Workspace,
) -> Atom {
    let id = name.into_id(state);
    atomic_flat_label_id(index, id, state, ws)
}

pub fn atomic_flat_label_id(
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
pub trait IntoId {
    fn into_id(self, state: &mut State) -> Identifier;
}

impl IntoId for SmartString<LazyCompact> {
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

impl IntoId for std::string::String {
    fn into_id(self, state: &mut State) -> Identifier {
        state.get_or_insert_fn(self, None).unwrap()
    }
}

/// Trait that enables shadowing of a tensor
///
/// This creates a dense tensor of atoms, where the atoms are the expanded indices of the tensor, with the global name as the name of the labels.
pub trait Shadowable: TensorStructure {
    type Name: IntoId + Clone;
    fn shadow(self, state: &mut State, ws: &Workspace) -> Option<DenseTensor<Atom, Self::Structure>>
    where
        Self: std::marker::Sized + HasName<Name = <Self as Shadowable>::Name>,
        Self::Structure: Clone,
    {
        let name = self.name()?.clone();
        Some(self.shadow_with(name, state, ws))
    }

    fn shadow_with(
        self,
        name: Self::Name,
        state: &mut State,
        ws: &Workspace,
    ) -> DenseTensor<Atom, Self::Structure>
    where
        Self: std::marker::Sized,
        Self::Structure: Clone,
    {
        let f_id = name.clone().into_id(state);
        let mut data = vec![];
        for index in self.index_iter() {
            data.push(atomic_expanded_label_id(&index, f_id, state, ws));
        }

        DenseTensor {
            data,
            structure: self.structure().clone(),
        }
    }

    fn to_symbolic(&self, state: &mut State, ws: &Workspace) -> Option<Atom>
    where
        Self: HasName<Name = <Self as Shadowable>::Name>,
    {
        Some(self.to_symbolic_with(self.name()?.clone(), state, ws))
    }

    fn to_symbolic_with(&self, name: Self::Name, state: &mut State, ws: &Workspace) -> Atom {
        let atoms = self
            .external_structure()
            .iter()
            .map(|slot| slot.to_symbolic(state, ws))
            .collect::<Vec<_>>();

        let mut value_builder = FunctionBuilder::new(name.into_id(state), state, ws);
        for atom in atoms {
            value_builder = value_builder.add_arg(atom.as_atom_view());
        }
        value_builder.finish().into_atom()
    }
}

impl<N> Shadowable for N
where
    N: TensorStructure + HasName,
    N::Name: IntoId + Clone,
{
    type Name = N::Name;
}

duplicate! {[
  N;
[HistoryStructure<std::string::String>];
[HistoryStructure<SmartString<LazyCompact>>];
]
impl std::fmt::Display for N
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let mut string = String::new();
        if let Some(global_name) = self.name() {
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
}
impl HistoryStructure<Identifier> {
    pub fn to_string(&self, state: &State) -> String {
        let mut string = String::new();
        if let Some(global_name) = self.name() {
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
