use bincode::{Decode, Encode};
use derive_more::{From, Into};
use typed_index_collections::TiVec;

use crate::momentum::{FourMomentum, ThreeMomentum};

#[derive(From, Into, Copy, Clone, Hash, Eq, Ord, PartialEq, PartialOrd, Encode, Decode)]
pub struct LoopIndex(pub usize);
#[derive(From, Into, Copy, Clone, Hash, Eq, Ord, PartialEq, PartialOrd, Encode, Decode)]
pub struct ExternalIndex(pub usize);

pub type LoopMomenta<T> = TiVec<LoopIndex, ThreeMomentum<T>>;
pub type ExternalThreeMomenta<T> = TiVec<ExternalIndex, ThreeMomentum<T>>;
pub type ExternalFourMomenta<T> = TiVec<ExternalIndex, FourMomentum<T>>;
