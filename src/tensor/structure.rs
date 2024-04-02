use ahash::AHashMap;
use derive_more::Add;
use derive_more::AddAssign;
use derive_more::Display;
use derive_more::From;
use derive_more::Into;
use duplicate::duplicate;
use indexmap::IndexMap;
use serde::Deserialize;
use serde::Serialize;

use smartstring::LazyCompact;
use smartstring::SmartString;
use std::borrow::Cow;
use std::fmt::Debug;
use symbolica::representations::ListIterator;

use std::i64;

use std::ops::Range;

use symbolica::coefficient::CoefficientView;

use symbolica::representations::AtomView;

use permutation::Permutation;

use symbolica::representations::{AsAtomView, Atom, FunctionBuilder, Symbol};
use symbolica::state::{State, Workspace};

use std::collections::HashSet;
use std::{cmp::Ordering, collections::HashMap};

use super::ufo;
use super::DenseTensor;
use super::MixedTensor;
use super::TensorStructureIndexIterator;
use smartstring::alias::String;
/// A type that represents the name of an index in a tensor.
#[derive(
    Debug,
    Copy,
    Clone,
    Ord,
    PartialOrd,
    Eq,
    PartialEq,
    Hash,
    Serialize,
    Deserialize,
    From,
    Into,
    Display,
    Add,
    AddAssign,
)]
#[display(fmt = "id{}", _0)]
pub struct AbstractIndex(pub usize);

impl TryFrom<AtomView<'_>> for AbstractIndex {
    type Error = String;

    fn try_from(view: AtomView<'_>) -> Result<Self, Self::Error> {
        if let AtomView::Var(v) = view {
            Ok(AbstractIndex(v.get_symbol().get_id() as usize))
        } else {
            Err("Not a var".to_string().into())
        }
    }
}

impl TryFrom<std::string::String> for AbstractIndex {
    type Error = String;

    fn try_from(value: std::string::String) -> Result<Self, Self::Error> {
        let atom = Atom::parse(&value)?;
        Self::try_from(atom.as_view())
    }
}

impl TryFrom<&'_ str> for AbstractIndex {
    type Error = String;

    fn try_from(value: &'_ str) -> Result<Self, Self::Error> {
        let atom = Atom::parse(value)?;
        Self::try_from(atom.as_view())
    }
}

/// A Dimension
#[derive(
    Debug,
    Copy,
    Clone,
    Ord,
    PartialOrd,
    Eq,
    PartialEq,
    Hash,
    Serialize,
    Deserialize,
    From,
    Into,
    Add,
    Display,
)]
#[into(owned, ref, ref_mut)]
#[display(fmt = "{}", _0)]
pub struct Dimension(pub usize);

impl PartialEq<usize> for Dimension {
    fn eq(&self, other: &usize) -> bool {
        self.0 == *other
    }
}

impl PartialEq<Dimension> for usize {
    fn eq(&self, other: &Dimension) -> bool {
        *self == other.0
    }
}

impl PartialOrd<usize> for Dimension {
    fn partial_cmp(&self, other: &usize) -> Option<Ordering> {
        self.0.partial_cmp(other)
    }
}

impl PartialOrd<Dimension> for usize {
    fn partial_cmp(&self, other: &Dimension) -> Option<Ordering> {
        self.partial_cmp(&other.0)
    }
}

/// A  concrete index, i.e. the concrete usize/index of the corresponding abstract index

pub type ConcreteIndex = usize;

pub const EUCLIDEAN: &str = "euc";
pub const LORENTZ: &str = "lor";
pub const BISPINOR: &str = "bis";
pub const SPINFUND: &str = "spin";
pub const SPINANTIFUND: &str = "spina";
pub const COLORADJ: &str = "coad";
pub const COLORFUND: &str = "cof";
pub const COLORANTIFUND: &str = "coaf";
pub const COLORSEXT: &str = "cos";
pub const COLORANTISEXT: &str = "coas";

/// A Representation/Dimension of the index.
#[derive(PartialEq, Eq, Clone, Copy, Debug, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Representation {
    /// Represents a Euclidean space of the given dimension, with metric diag(1,1,1,1,...)
    Euclidean(Dimension),
    /// Represents a Minkowski space of the given dimension, with metric diag(1,-1,-1,-1,...)
    Lorentz(Dimension),
    Bispinor(Dimension),
    /// Represents a Spinor Fundamental space of the given dimension
    SpinFundamental(Dimension),
    /// Represents a Spinor Adjoint space of the given dimension
    SpinAntiFundamental(Dimension),
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
    /// # use _gammaloop::tensor::Dimension;
    /// let spin = Representation::Bispinor(Dimension(5));
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
    #[must_use]
    pub fn negative(&self) -> Vec<bool> {
        match *self {
            Self::Lorentz(value) => std::iter::once(false)
                .chain(std::iter::repeat(true).take(value.0 - 1))
                .collect::<Vec<_>>(),
            Self::Euclidean(value)
            | Self::Bispinor(value)
            | Self::SpinFundamental(value)
            | Self::SpinAntiFundamental(value) => {
                vec![false; value.into()]
            }
            Self::ColorAdjoint(value)
            | Self::ColorFundamental(value)
            | Self::ColorAntiFundamental(value)
            | Self::ColorSextet(value)
            | Self::ColorAntiSextet(value) => vec![false; value.into()],
        }
    }

    /// for the given concrete index, says whether it should have a minus sign during contraction
    ///
    /// for example see [`Self::negative`]
    #[inline]
    #[must_use]
    pub const fn is_neg(&self, i: usize) -> bool {
        match self {
            Self::Lorentz(_) => i > 0,
            _ => false,
        }
    }

    /// yields a function builder for the representation, adding a first variable: the dimension.
    ///
    /// for example see [`Slot::to_symbolic`]
    #[allow(clippy::cast_possible_wrap)]
    pub fn to_fnbuilder<'a, 'b: 'a>(&'a self) -> FunctionBuilder {
        let (value, id) = match *self {
            Self::Euclidean(value) => (value, State::get_symbol(EUCLIDEAN)),
            Self::Lorentz(value) => (value, State::get_symbol(LORENTZ)),
            Self::Bispinor(value) => (value, State::get_symbol(BISPINOR)),
            Self::SpinFundamental(value) => (value, State::get_symbol(SPINFUND)),
            Self::SpinAntiFundamental(value) => (value, State::get_symbol(SPINANTIFUND)),
            Self::ColorAdjoint(value) => (value, State::get_symbol(COLORADJ)),
            Self::ColorFundamental(value) => (value, State::get_symbol(COLORFUND)),
            Self::ColorAntiFundamental(value) => (value, State::get_symbol(COLORANTIFUND)),
            Self::ColorSextet(value) => (value, State::get_symbol(COLORSEXT)),
            Self::ColorAntiSextet(value) => (value, State::get_symbol(COLORANTISEXT)),
        };

        let mut value_builder = FunctionBuilder::new(id);

        value_builder =
            value_builder.add_arg(Atom::new_num(usize::from(value) as i64).as_atom_view());

        value_builder
    }

    /// Finishes the function builder into an Atom
    ///
    /// # Example
    ///
    /// ```
    /// # use symbolica::state::{State, Workspace};
    /// # use _gammaloop::tensor::Representation;
    /// # use _gammaloop::tensor::Dimension;
    ///
    /// let mink = Representation::Lorentz(Dimension(4));
    ///
    /// assert_eq!("lor(4)",format!("{}",mink.to_symbolic()));
    /// assert_eq!("lor4",format!("{}",mink));
    /// ```
    pub fn to_symbolic(&self) -> Atom {
        self.to_fnbuilder().finish()
    }
}

impl From<Dimension> for Representation {
    fn from(value: Dimension) -> Self {
        Self::Euclidean(value)
    }
}

impl From<usize> for Representation {
    fn from(value: usize) -> Self {
        Self::Euclidean(value.into())
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
            | Representation::Bispinor(value)
            | Representation::SpinFundamental(value)
            | Representation::SpinAntiFundamental(value) => *value,
            Representation::ColorAdjoint(value) => *value, //Dimension(8),
            Representation::ColorFundamental(value)
            | Representation::ColorAntiFundamental(value) => {
                *value // Dimension(3)
            }
            Representation::ColorSextet(value) | Representation::ColorAntiSextet(value) => *value,
        }
    }
}

impl From<&Representation> for usize {
    fn from(value: &Representation) -> Self {
        usize::from(Dimension::from(value))
    }
}

impl From<Representation> for Dimension {
    fn from(rep: Representation) -> Self {
        match rep {
            Representation::Euclidean(value)
            | Representation::Lorentz(value)
            | Representation::Bispinor(value)
            | Representation::SpinFundamental(value)
            | Representation::SpinAntiFundamental(value) => value,
            Representation::ColorAdjoint(value) => value,
            Representation::ColorFundamental(value)
            | Representation::ColorAntiFundamental(value) => {
                value // Dimension(3)
            }
            Representation::ColorSextet(value) | Representation::ColorAntiSextet(value) => value, //Dimension(6),
        }
    }
}

impl From<Representation> for usize {
    fn from(value: Representation) -> Self {
        usize::from(Dimension::from(value))
    }
}

impl std::fmt::Display for Representation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Euclidean(value) => write!(f, "{EUCLIDEAN}{value}"),
            Self::Lorentz(value) => write!(f, "{LORENTZ}{value}"),
            Self::Bispinor(value) => write!(f, "{BISPINOR}{value}"),
            Self::SpinFundamental(value) => write!(f, "{SPINFUND}{value}"),
            Self::SpinAntiFundamental(value) => write!(f, "{SPINANTIFUND}{value}"),
            Self::ColorAdjoint(value) => write!(f, "{COLORADJ}{value}"),
            Self::ColorFundamental(value) => write!(f, "{COLORFUND}{value}"),
            Self::ColorAntiFundamental(value) => write!(f, "{COLORANTIFUND}{value}"),
            Self::ColorSextet(value) => write!(f, "{COLORSEXT}{value}"),
            Self::ColorAntiSextet(value) => write!(f, "{COLORANTISEXT}{value}"),
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
/// # use _gammaloop::tensor::{Representation,Slot,Dimension,AbstractIndex};
/// let mink = Representation::Lorentz(Dimension(4));
/// let mu = Slot::from((AbstractIndex(0),mink));
///
/// assert_eq!("id0lor4",format!("{}",mu));
/// ```
///
/// It can also be built from a tuple of `usize` and `usize`, where we default to `Representation::Euclidean`
/// ```
/// # use _gammaloop::tensor::{Representation,Slot};
/// let mu = Slot::from((0,4));
/// assert_eq!("id0euc4",format!("{}",mu));
/// ```
#[derive(Debug, Copy, Clone, Hash, PartialEq, Eq, Serialize, Deserialize)]
pub struct Slot {
    pub index: AbstractIndex,
    pub representation: Representation,
}

/// Can possibly constuct a Slot from an `AtomView`, if it is of the form: <representation>(<dimension>,<index>)
///
/// # Example
///
/// ```
///
/// # use _gammaloop::tensor::{Representation,Slot,Dimension,AbstractIndex};
/// # use symbolica::representations::AtomView;

///    let mink = Representation::Lorentz(Dimension(4));
///    let mu = Slot::from((AbstractIndex(0), mink));
///    let atom = mu.to_symbolic();
///    let slot = Slot::try_from(atom.as_view()).unwrap();
///    assert_eq!(slot, mu);
/// ```
impl TryFrom<AtomView<'_>> for Slot {
    type Error = &'static str;

    fn try_from(value: AtomView<'_>) -> Result<Self, Self::Error> {
        fn extract_num(iter: &mut ListIterator) -> Result<i64, &'static str> {
            if let Some(a) = iter.next() {
                if let AtomView::Num(n) = a {
                    if let CoefficientView::Natural(n, 1) = n.get_coeff_view() {
                        return Ok(n);
                    }
                    return Err("Argument is not a natural number");
                }
                Err("Argument is not a number")
            } else {
                Err("No more arguments")
            }
        }

        let mut iter = if let AtomView::Fun(f) = value {
            f.iter()
        } else {
            return Err("Not a slot, is composite");
        };

        let dim: Dimension = usize::try_from(extract_num(&mut iter)?)
            .or(Err("Dimension too large"))?
            .into();
        let index: AbstractIndex = usize::try_from(extract_num(&mut iter)?)
            .or(Err("Dimension too large"))?
            .into();

        if extract_num(&mut iter).is_ok() {
            return Err("Too many arguments");
        }

        let euc = State::get_symbol(EUCLIDEAN);

        let lor = State::get_symbol(LORENTZ);
        let bis = State::get_symbol(BISPINOR);
        let spin = State::get_symbol(SPINFUND);
        let spina = State::get_symbol(SPINANTIFUND);
        let coad = State::get_symbol(COLORADJ);
        let cof = State::get_symbol(COLORFUND);
        let coaf = State::get_symbol(COLORANTIFUND);
        let cos = State::get_symbol(COLORSEXT);
        let coas = State::get_symbol(COLORANTISEXT);

        let representation = if let AtomView::Fun(f) = value {
            let sym = f.get_symbol();
            match sym {
                _ if sym == euc => Representation::Euclidean(dim),
                _ if sym == lor => Representation::Lorentz(dim),
                _ if sym == bis => Representation::Bispinor(dim),
                _ if sym == spin => Representation::SpinFundamental(dim),
                _ if sym == spina => Representation::SpinAntiFundamental(dim),
                _ if sym == coad => Representation::ColorAdjoint(dim),
                _ if sym == cof => Representation::ColorFundamental(dim),
                _ if sym == coaf => Representation::ColorAntiFundamental(dim),
                _ if sym == cos => Representation::ColorSextet(dim),
                _ if sym == coas => Representation::ColorAntiSextet(dim),
                _ => return Err("Not a slot, isn't a representation"),
            }
        } else {
            return Err("Not a slot, is composite");
        };

        Ok(Slot {
            index,
            representation,
        })
    }
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
        Self {
            index: value.0,
            representation: value.1,
        }
    }
}

impl From<(usize, usize)> for Slot {
    fn from(value: (usize, usize)) -> Self {
        Self {
            index: value.0.into(),
            representation: value.1.into(),
        }
    }
}

#[allow(clippy::cast_possible_wrap)]
impl Slot {
    /// using the function builder of the representation add the abstract index as an argument, and finish it to an Atom.
    /// # Example
    ///
    /// ```
    /// # use symbolica::state::{State, Workspace};
    /// # use _gammaloop::tensor::{Representation,Slot,Dimension,AbstractIndex};
    /// let mink = Representation::Lorentz(Dimension(4));
    /// let mu = Slot::from((AbstractIndex(0),mink));
    ///
    /// assert_eq!("lor(4,0)",format!("{}",mu.to_symbolic()));
    /// assert_eq!("id0lor4",format!("{}",mu));
    /// ```
    pub fn to_symbolic(&self) -> Atom {
        let mut value_builder = self.representation.to_fnbuilder();
        value_builder =
            value_builder.add_arg(Atom::new_num(usize::from(self.index) as i64).as_atom_view());
        value_builder.finish()
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

    /// Given two [`TensorStructure`]s, returns the index of the first matching slot in each external index list, along with a boolean indicating if there is a single match
    fn match_index(&self, other: &Self) -> Option<(bool, usize, usize)> {
        let posmap = self
            .external_structure()
            .iter()
            .enumerate()
            .map(|(i, slot)| (slot, i))
            .collect::<AHashMap<_, _>>();

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

    /// Given two [`TensorStructure`]s, returns the index of the first matching slot in each external index list
    fn match_indices(&self, other: &Self) -> Option<(Permutation, Vec<bool>, Vec<bool>)> {
        let mut self_matches = vec![false; self.order()];
        let mut perm = Vec::new();
        let mut other_matches = vec![false; other.order()];

        let posmap = self
            .external_structure()
            .iter()
            .enumerate()
            .map(|(i, slot)| (slot, i))
            .collect::<AHashMap<_, _>>();

        for (j, slot_other) in other.external_structure().iter().enumerate() {
            if let Some(&i) = posmap.get(slot_other) {
                self_matches[i] = true;
                other_matches[j] = true;
                perm.push(j);
            }
        }

        if perm.is_empty() {
            None
        } else {
            let p: Permutation = permutation::sort(&mut perm);
            Some((p, self_matches, other_matches))
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

    /// find the permutation of the external indices that would make the two tensors the same. Applying the permutation to other should make it the same as self
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
        for item in self.external_structure() {
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
    ///
    /// # Errors
    ///
    /// `Mismatched order` = if the length of the indices is different from the order of the tensor,
    ///
    /// `Index out of bounds` = if the index is out of bounds for the dimension of that index   
    ///
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
    ///
    /// # Errors
    ///
    /// Same as [`Self::verify_indices`]
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
    ///
    /// # Errors
    ///
    /// `Index out of bounds` = if the flat index is out of bounds for the tensor
    fn expanded_index(&self, flat_index: usize) -> Result<Vec<ConcreteIndex>, String> {
        let mut indices = vec![];
        let mut index = flat_index;
        for &stride in &self.strides() {
            indices.push(index / stride);
            index %= stride;
        }
        if flat_index < self.size() {
            Ok(indices)
        } else {
            Err(format!("Index {flat_index} out of bounds").into())
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
        self.shape().iter().map(|x| usize::from(*x)).product()
    }

    fn shadow_with(self, f_id: Symbol) -> DenseTensor<Atom, Self>
    where
        Self: std::marker::Sized + Clone,
    {
        let mut data = vec![];
        for index in self.index_iter() {
            data.push(atomic_expanded_label_id(&index, f_id));
        }

        DenseTensor {
            data,
            structure: self,
        }
    }

    fn to_explicit_rep(self, f_id: Symbol) -> MixedTensor<Self>
    where
        Self: std::marker::Sized + Clone + TensorStructure,
    {
        let id = State::get_symbol("id");
        let gamma = State::get_symbol("γ");
        let gamma5 = State::get_symbol("γ5");
        let proj_m = State::get_symbol("ProjM");
        let proj_p = State::get_symbol("ProjP");
        let sigma = State::get_symbol("σ");

        match f_id {
            _ if f_id == id => ufo::identity_data::<f64, Self>(self).into(),

            _ if f_id == gamma => ufo::gamma_data(self).into(),
            _ if f_id == gamma5 => ufo::gamma5_data(self).into(),
            _ if f_id == proj_m => ufo::proj_m_data(self).into(),
            _ if f_id == proj_p => ufo::proj_p_data(self).into(),
            _ if f_id == sigma => ufo::sigma_data(self).into(),
            name => self.shadow_with(name).into(),
        }
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
    type Structure = Self;

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

    fn merge(&mut self, other: &Self) -> Option<usize>;

    #[must_use]
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

    fn merge(&mut self, other: &Self) -> Option<usize> {
        let mut positions = IndexMap::new();
        let mut i = 0;

        self.retain(|x| {
            let e = positions.get(x);
            if e.is_some() {
                return false;
            }
            positions.insert(*x, (Some(i), None));
            i += 1;
            true
        });

        let mut first = true;
        let mut first_other = 0;

        for (index, &value) in self.iter().enumerate() {
            positions.entry(value).or_insert((Some(index), None));
        }

        for (index, &value) in other.iter().enumerate() {
            let e = positions.get(&value);
            if let Some((Some(selfi), None)) = e {
                positions.insert(value, (Some(*selfi), Some(index)));
            } else {
                positions.insert(value, (None, Some(index)));
                self.push(value);
            }
        }

        let mut i = 0;

        self.retain(|x| {
            let pos = positions.get(x).unwrap();
            if pos.1.is_none() {
                i += 1;
                return true;
            }
            if pos.0.is_none() {
                if first {
                    first = false;
                    first_other = i;
                }
                return true;
            }
            false
        });

        if first {
            None
        } else {
            Some(first_other)
        }
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

/// A trait for a structure that can be traced and merged, during a contraction, maybe using symbolic state and workspace.

pub trait SymbolicStructureContract {
    fn trace_sym(&mut self, i: usize, j: usize, state: &State, ws: &Workspace);

    fn trace_out_sym(&mut self, state: &State, ws: &Workspace);

    fn merge_sym(&mut self, other: &Self, state: &State, ws: &Workspace);

    #[must_use]
    fn merge_at_sym(
        &self,
        other: &Self,
        positions: (usize, usize),
        state: &State,
        ws: &Workspace,
    ) -> Self;
}

impl<T> SymbolicStructureContract for T
where
    T: StructureContract,
{
    fn trace_sym(&mut self, i: usize, j: usize, _state: &State, _ws: &Workspace) {
        self.trace(i, j);
    }

    fn trace_out_sym(&mut self, _state: &State, _ws: &Workspace) {
        self.trace_out();
    }

    fn merge_sym(&mut self, other: &Self, _state: &State, _ws: &Workspace) {
        self.merge(other);
    }

    fn merge_at_sym(
        &self,
        other: &Self,
        positions: (usize, usize),
        _state: &State,
        _ws: &Workspace,
    ) -> Self {
        self.merge_at(other, positions)
    }
}

#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct VecStructure {
    pub structure: Vec<Slot>,
}

impl TryFrom<AtomView<'_>> for VecStructure {
    type Error = &'static str;
    fn try_from(value: AtomView) -> Result<Self, Self::Error> {
        let mut structure: Vec<Slot> = vec![];
        if let AtomView::Fun(f) = value {
            for arg in f.iter() {
                structure.push(arg.try_into()?);
            }
        } else {
            return Err("Not a valid expression");
        }
        Ok(structure.into())
    }
}

impl FromIterator<Slot> for VecStructure {
    fn from_iter<T: IntoIterator<Item = Slot>>(iter: T) -> Self {
        Self {
            structure: iter.into_iter().collect(),
        }
    }
}

impl VecStructure {
    pub fn new(structure: Vec<Slot>) -> Self {
        Self { structure }
    }

    pub fn to_named(self, name: &str) -> NamedStructure {
        NamedStructure::from_slots(self.structure, name)
    }
}

impl From<ContractionCountStructure> for VecStructure {
    fn from(structure: ContractionCountStructure) -> Self {
        Self {
            structure: structure.structure,
        }
    }
}

impl From<VecStructure> for ContractionCountStructure {
    fn from(structure: VecStructure) -> Self {
        Self {
            structure: structure.structure,
            contractions: 0,
        }
    }
}

impl From<Vec<Slot>> for VecStructure {
    fn from(structure: Vec<Slot>) -> Self {
        Self { structure }
    }
}

impl From<VecStructure> for Vec<Slot> {
    fn from(structure: VecStructure) -> Self {
        structure.structure
    }
}

// const IDPRINTER: Lazy<BlockId<char>> = Lazy::new(|| BlockId::new(Alphabet::alphanumeric(), 1, 1));

impl std::fmt::Display for VecStructure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (index, item) in self.structure.iter().enumerate() {
            if index != 0 {
                // To avoid a newline at the start
                writeln!(f)?;
            }
            write!(
                f,
                "{:<3} ({})",
                usize::from(item.index),
                // IDPRINTER
                //     .encode_string(usize::from(item.index) as u64)
                //     .unwrap(),
                item.representation
            )?;
        }
        Ok(())
    }
}

impl TensorStructure for VecStructure {
    type Structure = VecStructure;
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

impl StructureContract for VecStructure {
    fn merge(&mut self, other: &Self) -> Option<usize> {
        self.structure.merge(&other.structure)
    }

    fn trace_out(&mut self) {
        self.structure.trace_out();
    }

    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self {
        Self {
            structure: self.structure.merge_at(&other.structure, positions),
        }
    }

    fn trace(&mut self, i: usize, j: usize) {
        self.structure.trace(i, j);
    }
}

/// A named structure is a structure with a global name, and a list of slots
///
/// It is useful when you want to shadow tensors, to nest tensor network contraction operations.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct NamedStructure {
    pub structure: Vec<Slot>,
    pub global_name: Option<SmartString<LazyCompact>>,
}

impl NamedStructure {
    /// Constructs a new [`NamedStructure`] from a list of tuples of indices and dimension (assumes they are all euclidean), along with a name
    #[must_use]
    pub fn from_integers(slots: &[(AbstractIndex, Dimension)], name: &str) -> Self {
        let slots: Vec<(AbstractIndex, Representation)> = slots
            .iter()
            .map(|(index, dim)| (*index, Representation::Euclidean(*dim)))
            .collect();
        Self::new(&slots, name)
    }
    /// Constructs a new [`NamedStructure`] from a list of tuples of indices and representations, along with a name
    #[must_use]
    pub fn new(slots: &[(AbstractIndex, Representation)], name: &str) -> Self {
        let structure: Vec<Slot> = slots
            .iter()
            .map(|(index, representation)| Slot::from((*index, *representation)))
            .collect();

        Self {
            structure,
            global_name: Some(name.into()),
        }
    }

    pub fn from_slots(slots: Vec<Slot>, name: &str) -> Self {
        Self {
            structure: slots,
            global_name: Some(name.into()),
        }
    }
}

/// A trait for a structure that has a name
pub trait HasName {
    type Name: Clone;
    fn name(&self) -> Option<Cow<Self::Name>>;
    fn set_name(&mut self, name: &Self::Name);
}

impl HasName for NamedStructure {
    type Name = SmartString<LazyCompact>;
    fn name(&self) -> Option<Cow<Self::Name>> {
        self.global_name.as_ref().map(Cow::Borrowed)
    }
    fn set_name(&mut self, name: &Self::Name) {
        self.global_name = Some(name.clone());
    }
}

impl TensorStructure for NamedStructure {
    type Structure = Self;
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
    fn merge(&mut self, other: &Self) -> Option<usize> {
        self.structure.merge(&other.structure)
    }

    fn trace_out(&mut self) {
        self.structure.trace_out();
    }

    /// when merging two named structures, the global name is lost
    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self {
        Self {
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
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct ContractionCountStructure {
    pub structure: Vec<Slot>,
    pub contractions: usize,
}

impl FromIterator<Slot> for ContractionCountStructure {
    fn from_iter<T: IntoIterator<Item = Slot>>(iter: T) -> Self {
        Self {
            structure: iter.into_iter().collect(),
            contractions: 0,
        }
    }
}

impl ContractionCountStructure {
    /// Constructs a new [`ContractionCountStructure`] from a list of tuples of indices and dimension (assumes they are all euclidean), along with a name
    #[must_use]
    pub fn from_integers(slots: &[(AbstractIndex, Dimension)]) -> Self {
        let slots: Vec<(AbstractIndex, Representation)> = slots
            .iter()
            .map(|(index, dim)| (*index, Representation::Euclidean(*dim)))
            .collect();
        Self::new(&slots)
    }
    /// Constructs a new [`ContractionCountStructure`] from a list of tuples of indices and representations, along with a name
    #[must_use]
    pub fn new(slots: &[(AbstractIndex, Representation)]) -> Self {
        let structure: Vec<Slot> = slots
            .iter()
            .map(|(index, representation)| Slot::from((*index, *representation)))
            .collect();

        Self {
            structure,
            contractions: 0,
        }
    }

    pub fn from_slots(slots: Vec<Slot>) -> Self {
        Self {
            structure: slots,
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
    fn merge(&mut self, other: &Self) -> Option<usize> {
        self.contractions += other.contractions + 1;
        self.structure.merge(&other.structure)
    }

    fn trace_out(&mut self) {
        self.structure.trace_out();
    }

    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self {
        Self {
            structure: self.structure.merge_at(&other.structure, positions),
            contractions: self.contractions + other.contractions + 1,
        }
    }

    fn trace(&mut self, i: usize, j: usize) {
        self.structure.trace(i, j);
    }
}

/// A structure to enable smart shadowing of tensors in a tensor network contraction algorithm.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub struct SmartShadowStructure {
    pub structure: Vec<Slot>,
    pub contractions: usize,
    pub global_name: Option<SmartString<LazyCompact>>,
}

impl SmartShadowStructure {
    /// Constructs a new [`SmartShadow`] from a list of tuples of indices and dimension (assumes they are all euclidean), along with a name
    #[must_use]
    pub fn from_integers(slots: &[(AbstractIndex, Dimension)], name: &str) -> Self {
        let slots: Vec<(AbstractIndex, Representation)> = slots
            .iter()
            .map(|(index, dim)| (*index, Representation::Euclidean(*dim)))
            .collect();
        Self::new(&slots, name)
    }
    /// Constructs a new [`SmartShadow`] from a list of tuples of indices and representations, along with a name
    #[must_use]
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
    fn name(&self) -> Option<Cow<SmartString<LazyCompact>>> {
        self.global_name.as_ref().map(Cow::Borrowed)
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
    fn merge(&mut self, other: &Self) -> Option<usize> {
        self.contractions += other.contractions;
        self.structure.merge(&other.structure)
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
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct HistoryStructure<N> {
    internal: Vec<Slot>,
    pub external: Vec<Slot>,
    pub names: AHashMap<Range<usize>, N>, //ideally this is a named partion.. maybe a btreemap<usize, N>, and the range is from previous to next
    pub global_name: Option<N>,
}

impl<N> HistoryStructure<N> {
    /// Constructs a new [`HistoryStructure`] from a list of tuples of indices and dimension (assumes they are all euclidean), along with a name
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
    /// Constructs a new [`HistoryStructure`] from a list of tuples of indices and representations, along with a name
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
            .unwrap_or(0.into())
            + 1.into();

        for item in &mut self.internal {
            if other_set.contains(item) {
                item.index = replacement_value;
                replacement_value += 1.into();
            }
        }
    }
}

impl<N> HasName for HistoryStructure<N>
where
    N: Clone,
{
    type Name = N;
    fn name(&self) -> Option<Cow<N>> {
        self.global_name.as_ref().map(|name| Cow::Borrowed(name))
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
    fn merge(&mut self, other: &Self) -> Option<usize> {
        let shift = self.internal.len();
        for (range, name) in &other.names {
            self.names
                .insert((range.start + shift)..(range.end + shift), name.clone());
        }
        self.trace_out();
        self.independentize_internal(other);
        self.internal.append(&mut other.internal.clone());
        self.external.merge(&other.external)
    }

    /// Merge two [`HistoryStructure`] at the given positions of the external index list. Ideally the internal index list should be independentized before merging
    /// This is essentially a contraction of only one index. The name maps are merged, and shifted accordingly. The global name is lost, since the resulting tensor is composite
    /// The global name can be set again with the [`Self::set_global_name`] function
    fn merge_at(&self, other: &Self, positions: (usize, usize)) -> Self {
        let mut slots_other = other.external.clone();
        let mut slots_self: Vec<Slot> = self.external.clone();

        slots_self.remove(positions.0);
        slots_other.remove(positions.1);

        let mut slots_self_int = self.internal.clone();
        let mut slots_other_int = other.internal.clone();
        slots_self_int.append(&mut slots_other_int);

        let mut names = self.names.clone();
        let shift = self.internal.len();
        for (range, name) in &other.names {
            names.insert((range.start + shift)..(range.end + shift), name.clone());
        }
        slots_self.append(&mut slots_other);
        HistoryStructure {
            internal: slots_self_int,
            external: slots_self,
            names,
            global_name: None,
        }
    }
}

pub fn atomic_expanded_label<I: IntoId>(
    indices: &[ConcreteIndex],
    name: I,
    _state: &mut State,
    _ws: &Workspace,
) -> Atom {
    let id = name.into_id();
    atomic_expanded_label_id(indices, id)
}

pub fn atomic_flat_label<I: IntoId>(index: usize, name: I) -> Atom {
    let id = name.into_id();
    atomic_flat_label_id(index, id)
}

#[allow(clippy::cast_possible_wrap)]
pub fn atomic_flat_label_id(index: usize, id: Symbol) -> Atom {
    let mut value_builder = FunctionBuilder::new(id);
    value_builder = value_builder.add_arg(Atom::new_num(index as i64).as_atom_view());
    value_builder.finish()
}

#[allow(clippy::cast_possible_wrap)]
pub fn atomic_expanded_label_id(indices: &[ConcreteIndex], id: Symbol) -> Atom {
    let mut value_builder = FunctionBuilder::new(id);
    for &index in indices {
        value_builder = value_builder.add_arg(Atom::new_num(index as i64).as_atom_view());
    }
    value_builder.finish()
}
pub trait IntoId {
    fn into_id(self) -> Symbol;
}

impl IntoId for SmartString<LazyCompact> {
    fn into_id(self) -> Symbol {
        State::get_symbol(self)
    }
}

impl IntoId for Symbol {
    fn into_id(self) -> Symbol {
        self
    }
}

impl IntoId for &str {
    fn into_id(self) -> Symbol {
        State::get_symbol(self)
    }
}

impl IntoId for std::string::String {
    fn into_id(self) -> Symbol {
        State::get_symbol(self)
    }
}

/// Trait that enables shadowing of a tensor
///
/// This creates a dense tensor of atoms, where the atoms are the expanded indices of the tensor, with the global name as the name of the labels.
pub trait Shadowable: TensorStructure {
    type Name: IntoId + Clone;
    fn shadow(self) -> Option<DenseTensor<Atom, Self::Structure>>
    where
        Self: std::marker::Sized + HasName<Name = <Self as Shadowable>::Name>,
        Self::Structure: Clone + TensorStructure,
    {
        let name = self.name()?.into_owned();

        Some(self.structure().clone().shadow_with(name.into_id()))
    }

    fn smart_shadow(self) -> Option<MixedTensor<Self::Structure>>
    where
        Self: std::marker::Sized + HasName<Name = <Self as Shadowable>::Name>,
        Self::Structure: Clone + TensorStructure,
    {
        let name = self.name()?.into_owned();
        Some(self.structure().clone().to_explicit_rep(name.into_id()))
    }

    fn to_symbolic(&self) -> Option<Atom>
    where
        Self: HasName<Name = <Self as Shadowable>::Name>,
    {
        Some(self.to_symbolic_with(self.name()?.into_owned()))
    }

    fn to_symbolic_with(&self, name: Self::Name) -> Atom {
        let atoms = self
            .external_structure()
            .iter()
            .map(|slot| slot.to_symbolic())
            .collect::<Vec<_>>();

        let mut value_builder = FunctionBuilder::new(name.into_id());
        for atom in atoms {
            value_builder = value_builder.add_arg(atom.as_atom_view());
        }
        value_builder.finish()
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
            string.push_str(&format!("{global_name}:"));
        }
        for (range, name) in self
            .names
            .iter()
            .filter(|(r, _)| *r != &(0..self.internal.len()) || !self.is_composite())
        {
            string.push_str(&format!("{name}("));
            for slot in &self.internal[range.clone()] {
                string.push_str(&format!("{slot},"));
            }
            string.pop();
            string.push(')');
        }
        write!(f, "{string}")
    }
}
}
impl HistoryStructure<Symbol> {
    #[must_use]
    pub fn to_string(&self, _state: &State) -> String {
        let mut string = String::new();
        if let Some(global_name) = self.name() {
            string.push_str(&format!("{:?}:", global_name));
        }
        for (range, name) in self
            .names
            .iter()
            .filter(|(r, _)| *r != &(0..self.internal.len()) || !self.is_composite())
        {
            string.push_str(&format!("{:?}(", name));
            for slot in &self.internal[range.clone()] {
                string.push_str(&format!("{slot},"));
            }
            string.pop();
            string.push(')');
        }
        string
    }
}
