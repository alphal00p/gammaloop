use crate::momentum::SignOrZero;
use crate::momentum_sample::{ExternalIndex, LoopIndex};
use bincode::{Decode, Encode};
use serde::{Deserialize, Serialize};
use std::fmt::Display;

#[derive(
    Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, PartialOrd, Eq, Ord, Hash,
)]
pub struct LoopSignature(Vec<SignOrZero>);
#[derive(
    Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, PartialOrd, Eq, Ord, Hash,
)]
pub struct ExternalSignature(Vec<SignOrZero>);

#[derive(
    Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, PartialOrd, Eq, Ord, Hash,
)]
pub struct NewLoopExtSignature {
    pub internal: LoopSignature,
    pub external: ExternalSignature,
}
