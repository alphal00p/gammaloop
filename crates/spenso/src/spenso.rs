/*!

Contains all the tooling for working with arbitrary rank tensors, symbolically, numerically, and parametrically.

It includes special support for a minkowski metric, and a way to add any custom diagonal (sign based) metric.

All tensor types make use of a tensor structure type, either  the minimum `Vec` of [`Slot`]s or a more complete (but slightly more computationally heavy) [`HistoryStructure`].
Data is then added, to make parametric, or fully numeric tensors.
If no data is added, some [`TensorStructure`]s behave like symbolic tensors: namely [`HistoryStructure`] and [`SymbolicTensor`]

There are two main types of data tensors, [`DenseTensor`] and [`SparseTensor`].
They each implement a different type of storage for data.

All types of tensors can be contracted together using the [`Contract`] trait.
This can be done manually, or using a [`TensorNetwork`] and specifying a contraction algorithm.

Several Enums are defined to be able to store heterogenous tensors.
Namely
- [`DataTensor`]
- [`NumTensor`]
- [`MixedTensor`]

*/
extern crate self as spenso;
/// All tooling for tensor structures, indices and representations
// pub use structure::*;
pub mod algebra;
//
/// Tensor contraction
pub mod contraction;
pub mod iterators;
pub mod network;
pub mod structure;
pub mod tensors;

#[cfg(feature = "shadowing")]
pub mod shadowing;

pub mod utils;
