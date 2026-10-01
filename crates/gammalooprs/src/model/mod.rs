//! GammaLoop runtime behavior for the canonical [`feynkit_model::Model`].
//!
//! Serializable model records, validation, and stable entity identifiers live
//! in `feynkit-model`. This module only adds GammaLoop-specific symbolic and
//! tensor operations through extension traits.

mod canonical;

pub use canonical::*;
