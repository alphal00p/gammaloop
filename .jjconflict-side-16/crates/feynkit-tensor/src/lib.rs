//! Symmetry-aware covariant tensor reduction for loop integrals.
//!
//! A configured external basis separates longitudinal loop components from
//! their orthogonal complement. Symbolica inverts the external Gram matrix;
//! the same vacuum projector then reduces the transverse tensors in `D-n`
//! dimensions. Scalar products remain explicit for later integral reduction.
//!
//! The universal coefficient engine uses orthogonal Weingarten functions.
//! Their value depends only on the coset type of two metric pairings, so a
//! rank-20 table requires 42 classes rather than a `19!!`-by-`19!!` Gram
//! matrix. Fully contracted reductions additionally project onto the repeated-
//! vector orbit spaces on both sides and solve the smaller one. For example,
//! rank-20 multiplicities `[2, 18]` on each side require a 2-by-2 exact solve.

#![forbid(unsafe_code)]

mod orbit_projector;
mod reduction;
mod weingarten;

pub use reduction::{
    ContractionOrbit, FeynmanDiagramTensorExt, TensorReducer, TensorReduction,
    TensorReductionError, TensorReductionTerm,
};
pub use weingarten::{CosetType, OrthogonalWeingarten, WeingartenError};
