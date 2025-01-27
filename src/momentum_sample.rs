use typed_index_collections::TiVec;

use crate::momentum::{FourMomentum, ThreeMomentum};

pub struct LoopIndex(pub usize);
pub struct ExternalIndex(pub usize);

pub type LoopMomenta<T> = TiVec<LoopIndex, ThreeMomentum<T>>;
pub type ExternalThreeMomenta<T> = TiVec<ExternalIndex, ThreeMomentum<T>>;
pub type ExternalFourMomenta<T> = TiVec<ExternalIndex, FourMomentum<T>>;
