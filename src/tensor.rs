/*!

 Contains all the tooling for working with arbitrary rank tensors, symbolically, numerically, and parametrically.

 It includes special support for a minkowski metric, and a way to add any custom diagonal (sign based) metric.

 All tensor types make use of a tensor structure type, either  the minimum `Vec` of [`Slot`]s or a more complete (but slightly more computationally heavy) [`TensorSkeleton`].
 Data is then added, to make parametric, or fully numeric tensors.
 If no data is added, a TensorSkeleton behaves like a symbolic tensor, contractions amounting to appending slots, and keeping track of labels

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


use symbolica::{
    representations::{Atom, AtomBuilder},
    state::BufferHandle,
};
pub type Expr<'a> = AtomBuilder<'a, BufferHandle<'a, Atom>>;

/// All tooling for tensor structures, indices and representations
pub mod tensor_structure;
pub use tensor_structure::*;

/// More ergonomic, and smart arithmatic with symbolic types
pub mod upgrading_arithmetic;
pub use upgrading_arithmetic::*;

/// Tensors with data
pub mod data_tensor;
pub use data_tensor::*;

/// Parametric tensor contraction
pub mod mixed_tensor;
pub use mixed_tensor::*;

/// Symbolic tensors
pub mod symbolic_tensor;
/// Iterators on fibers of tensors
pub mod tensor_iterator;
pub use tensor_iterator::*;

/// Tensor contraction
pub mod tensor_contraction;
pub use tensor_contraction::*;
/// Adding, subtracting, scalar multiplication of tensors
pub mod tensor_arithmetic;

/// Tensors as defined in the UFO format
pub mod ufo_spin_tensors;

#[cfg(test)]
mod tests;
