use std::{collections::BTreeMap, ops::Bound::Included};
pub mod tensor_structure;
use symbolica::{
    representations::{Atom, AtomBuilder},
    state::BufferHandle,
};
pub type Expr<'a> = AtomBuilder<'a, BufferHandle<'a, Atom>>;

pub mod upgrading_arithmetic;
pub use upgrading_arithmetic::*;

pub mod fourthroots;
pub use fourthroots::*;

pub use tensor_structure::*;
pub mod num_tensor;
pub use num_tensor::*;
pub mod tensor_iterator;
pub use tensor_iterator::*;

// #[allow(dead_code)]
// #[enum_dispatch(HasTensorStructure)]
// enum Tensor<T> {
//     Num(NumTensor<T>),
//     Symbolic(SymbolicTensor),
// }

pub mod mixed_tensor;

pub mod symbolic_tensor;
pub mod tensor_arithmetic;
pub mod tensor_contraction;
pub use tensor_contraction::*;
pub mod ufo_spin_tensors;

#[cfg(test)]
mod tests;
