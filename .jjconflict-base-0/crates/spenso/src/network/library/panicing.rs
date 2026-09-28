use std::{fmt::Display, marker::PhantomData};

use crate::network::library::{DummyKey, FunctionLibrary, FunctionLibraryError};

pub struct ErroringLibrary<K = DummyKey> {
    key: PhantomData<K>,
}

impl<K> Default for ErroringLibrary<K> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K> ErroringLibrary<K> {
    pub fn new() -> Self {
        Self { key: PhantomData }
    }
}

impl<T, S, K: Display + Clone> FunctionLibrary<T, S> for ErroringLibrary<K> {
    type Key = K;

    fn apply(&self, key: &Self::Key, _tensor: T) -> Result<T, FunctionLibraryError<K>> {
        Err(FunctionLibraryError::NotFound(key.clone()))
    }

    fn apply_scalar(
        &self,
        key: &Self::Key,
        _scalar: S,
    ) -> eyre::Result<S, FunctionLibraryError<Self::Key>> {
        Err(FunctionLibraryError::NotFound(key.clone()))
    }
}
