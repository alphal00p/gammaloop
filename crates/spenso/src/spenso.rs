/*!

Contains all the tooling for working with arbitrary rank tensors, symbolically, numerically, and parametrically.

It includes special support for a minkowski metric, and a way to add any custom diagonal (sign based) metric.

All tensor types use a [`TensorStructure`](structure::TensorStructure), from a
minimal `Vec` of [`Slot`](structure::slot::Slot)s to named
[`OrderedStructure`](structure::OrderedStructure) and
[`NamedStructure`](structure::NamedStructure) values. Data is then added to make
parametric or fully numeric tensors.

There are two main types of data tensors, [`DenseTensor`](tensors::data::DenseTensor)
and [`SparseTensor`](tensors::data::SparseTensor).
They each implement a different type of storage for data.

Tensor types can be contracted together using the [`Contract`](contraction::Contract)
trait. This can be done manually or through a [`Network`](network::Network) with
a selected contraction strategy.

[`DataTensor`](tensors::data::DataTensor) stores either dense or sparse data.
With the `shadowing` feature, [`MixedTensor`](tensors::parametric::MixedTensor)
stores heterogeneous numeric and parametric tensors.

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
#[cfg(feature = "shadowing")]
pub mod symbolic_parallelism;
#[cfg(feature = "shadowing")]
pub mod symbolica_init;

pub mod utils;
