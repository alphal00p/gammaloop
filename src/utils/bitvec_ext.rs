use std::ops::{Deref, DerefMut};

use bincode::Encode;
use bitvec::vec::BitVec;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct BinVec(pub BitVec);

impl Deref for BinVec {
    type Target = BitVec;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for BinVec {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Encode for BinVec {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        bincode::Encode::encode(&self.0.as_raw_slice(), encoder)
    }
}

impl<Context> bincode::Decode<Context> for BinVec {
    fn decode<D: bincode::de::Decoder<Context = Context>>(
        decoder: &mut D,
    ) -> core::result::Result<Self, bincode::error::DecodeError> {
        let vec: Vec<usize> = bincode::Decode::decode(decoder)?;
        Ok(BinVec(BitVec::from_vec(vec)))
    }
}
