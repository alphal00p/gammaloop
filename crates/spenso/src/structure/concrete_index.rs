use std::ops::Deref;

use bincode::Decode;
use bincode::Encode;
use derive_more::Add;
use derive_more::AddAssign;
use derive_more::Display;
use derive_more::From;
use derive_more::Index;
use derive_more::Into;
use derive_more::IntoIterator;
use derive_more::Mul;
use derive_more::MulAssign;
use derive_more::Rem;
use derive_more::RemAssign;
use derive_more::Sub;
use derive_more::SubAssign;
#[cfg(feature = "shadowing")]
use symbolica::{
    atom::{Atom, AtomCore, FunctionBuilder},
    {function, symbol},
};

use serde::{Deserialize, Serialize};

use linnet::permutation::Permutation;

pub const CONCRETEIND: &str = "cind";
pub const FLATIND: &str = "find";
pub const UP: &str = "u";
pub const DOWN: &str = "d";

/// A  concrete index, i.e. the concrete usize/index of the corresponding abstract index
pub type ConcreteIndex = usize;

#[derive(
    Debug, Copy, Clone, Ord, PartialOrd, Eq, PartialEq, Hash, Serialize, Deserialize, Display,
)]
pub enum DualConciousIndex {
    Up(ConcreteIndex),
    Down(ConcreteIndex),
    SelfDual(ConcreteIndex),
}

#[cfg(feature = "shadowing")]
impl From<DualConciousIndex> for Atom {
    fn from(value: DualConciousIndex) -> Self {
        match value {
            DualConciousIndex::Up(s) => Atom::num(s as i64),
            DualConciousIndex::Down(s) => function!(symbol!(DOWN), Atom::num(s as i64)),
            DualConciousIndex::SelfDual(s) => Atom::num(s as i64),
        }
    }
}

#[derive(
    Debug,
    Clone,
    Ord,
    PartialOrd,
    Eq,
    PartialEq,
    Hash,
    Index,
    Serialize,
    Deserialize,
    From,
    Into,
    Display,
    IntoIterator,
)]
#[display(fmt = "{:?}", indices)]
pub struct DualConciousExpandedIndex {
    indices: Vec<DualConciousIndex>,
}

impl Deref for DualConciousExpandedIndex {
    type Target = [DualConciousIndex];

    fn deref(&self) -> &Self::Target {
        &self.indices
    }
}

impl DualConciousExpandedIndex {
    pub fn permute(&mut self, perm: &Permutation) {
        perm.apply_slice_in_place(&mut self.indices);
    }
}

#[cfg(feature = "shadowing")]
impl From<DualConciousExpandedIndex> for Atom {
    fn from(value: DualConciousExpandedIndex) -> Self {
        let mut cind = FunctionBuilder::new(symbol!(CONCRETEIND));
        for i in value.iter() {
            cind = cind.add_arg(Atom::from(*i).as_atom_view());
        }
        cind.finish()
    }
}

#[derive(
    Debug,
    Clone,
    Ord,
    PartialOrd,
    Eq,
    PartialEq,
    Hash,
    Index,
    Serialize,
    Deserialize,
    From,
    Into,
    Display,
    IntoIterator,
    Default,
    Encode,
    Decode,
)]
#[display(fmt = "{:?}", indices)]
pub struct ExpandedIndex {
    pub indices: Vec<ConcreteIndex>,
}

impl Extend<ConcreteIndex> for ExpandedIndex {
    fn extend<T: IntoIterator<Item = ConcreteIndex>>(&mut self, iter: T) {
        self.indices.extend(iter);
    }
}
impl ExpandedIndex {
    pub fn apply_permutation(&self, permutation: &Permutation) -> Self {
        ExpandedIndex {
            indices: permutation.apply_slice(&self.indices),
        }
    }

    pub fn apply_inverse_permutation(&self, permutation: &Permutation) -> Self {
        ExpandedIndex {
            indices: permutation.apply_slice_inv(&self.indices),
        }
    }
}

impl AsRef<[ConcreteIndex]> for ExpandedIndex {
    fn as_ref(&self) -> &[ConcreteIndex] {
        &self.indices
    }
}

#[cfg(feature = "shadowing")]
impl From<ExpandedIndex> for Atom {
    fn from(value: ExpandedIndex) -> Self {
        let mut cind = FunctionBuilder::new(super::abstract_index::AIND_SYMBOLS.cind);
        for i in value.iter() {
            cind = cind.add_arg(Atom::num(*i as i64).as_atom_view());
        }
        cind.finish()
    }
}

impl Deref for ExpandedIndex {
    type Target = [ConcreteIndex];

    fn deref(&self) -> &Self::Target {
        &self.indices
    }
}

impl FromIterator<ConcreteIndex> for ExpandedIndex {
    fn from_iter<T: IntoIterator<Item = ConcreteIndex>>(iter: T) -> Self {
        ExpandedIndex {
            indices: iter.into_iter().collect(),
        }
    }
}

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
    Rem,
    RemAssign,
    Into,
    Display,
    Add,
    Mul,
    MulAssign,
    AddAssign,
    Sub,
    SubAssign,
    Encode,
    Decode,
)]
#[display(fmt = "{}", index)]
pub struct FlatIndex {
    index: usize,
}

#[cfg(feature = "shadowing")]
impl From<FlatIndex> for Atom {
    fn from(value: FlatIndex) -> Self {
        let mut cind = FunctionBuilder::new(symbol!(FLATIND));
        cind = cind.add_arg(Atom::num(value.index as i64).as_atom_view());
        cind.finish()
    }
}
