//! Interned immutable scopes keep abstract indices cheap to copy. Serialization
//! writes the scope and original index, never a process-local interner identity.

use std::{
    collections::HashSet,
    sync::{LazyLock, RwLock},
};

use serde::{Deserialize, Serialize};
use symbolica::{atom::Symbol, state::HasStateMap};
use symbolica_utils::SerializableSymbol;

use super::AbstractIndex;

#[derive(
    Debug,
    Copy,
    Clone,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Hash,
    Serialize,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    Deserialize,
)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
struct ScopedIndexData {
    scope: SerializableSymbol,
    index: AbstractIndex,
}

#[derive(
    Debug,
    Copy,
    Clone,
    Eq,
    PartialEq,
    Ord,
    PartialOrd,
    Hash,
    Serialize,
    bincode_trait_derive::Encode,
)]
pub struct ScopedIndex(&'static ScopedIndexData);

static SCOPES: LazyLock<RwLock<HashSet<&'static ScopedIndexData>>> =
    LazyLock::new(|| RwLock::new(HashSet::new()));

impl ScopedIndex {
    pub(super) fn new(scope: Symbol, index: AbstractIndex) -> Self {
        let value = ScopedIndexData {
            scope: scope.into(),
            index,
        };
        if let Some(value) = SCOPES.read().unwrap().get(&value) {
            return Self(value);
        }
        let mut scopes = SCOPES.write().unwrap();
        if let Some(value) = scopes.get(&value) {
            return Self(value);
        }
        // Like Symbolica symbols, immutable index identities live for the process.
        let value = Box::leak(Box::new(value));
        scopes.insert(value);
        Self(value)
    }

    pub fn scope(self) -> Symbol {
        self.0.scope.into()
    }
    pub fn index(self) -> AbstractIndex {
        self.0.index
    }
}

impl<'de> Deserialize<'de> for ScopedIndex {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let data = ScopedIndexData::deserialize(deserializer)?;
        Ok(Self::new(data.scope.into(), data.index))
    }
}

impl<C: HasStateMap> bincode::Decode<C> for ScopedIndex {
    fn decode<D: bincode::de::Decoder<Context = C>>(
        decoder: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        let data = ScopedIndexData::decode(decoder)?;
        Ok(Self::new(data.scope.into(), data.index))
    }
}
