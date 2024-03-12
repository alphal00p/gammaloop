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

/// All tooling for tensor structures, indices and representations
pub mod structure;
pub use structure::*;

/// More ergonomic, and smart arithmatic with symbolic types
pub mod upgrading_arithmetic;
pub use upgrading_arithmetic::*;

/// Tensors with data
pub mod data;
pub use data::*;

/// Parametric tensor contraction
pub mod parametric;
pub use parametric::*;

/// Symbolic tensors
pub mod symbolic;
pub use symbolic::*;
/// Iterators on fibers of tensors
pub mod iterators;
pub use iterators::*;

/// Tensor contraction
pub mod contraction;
pub use contraction::*;

/// Tensor networks
pub mod network;
pub use network::*;
/// Adding, subtracting, scalar multiplication of tensors
pub mod arithmetic;

/// Tensors as defined in the UFO format
pub mod ufo;

#[cfg(test)]
mod tests;
