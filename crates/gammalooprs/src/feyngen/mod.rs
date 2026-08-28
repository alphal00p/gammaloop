//! GammaLoop's downstream runtime extension for canonical FeynKit generation.

pub mod feynkit;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum FeynGenError {
    #[error("generation interrupted by user")]
    Interrupted,
    #[error(transparent)]
    Runtime(#[from] feynkit::FeynkitRuntimeError),
    #[error(transparent)]
    Eyre(#[from] color_eyre::Report),
}
