use std::{borrow::Cow, fmt::Display, marker::PhantomData, ops::Neg};

use crate::{
    algebra::complex::{Complex, RealOrComplex},
    network::StructureLessDisplay,
    structure::{
        HasStructure, PermutedStructure, StructureError, TensorStructure,
        concrete_index::ConcreteIndex, dimension::Dimension, representation::Representation,
        slot::IsAbstractSlot,
    },
    tensors::{
        complex::RealOrComplexTensor,
        data::{DataTensor, DenseTensor, SetTensorData, SparseTensor},
    },
};

use eyre::Result;

#[cfg(feature = "shadowing")]
pub mod function_lib;
pub mod panicing;
use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Debug, Clone, Error)]
pub enum LibraryError<Key: Display> {
    #[error("Not found {0}")]
    NotFound(Key),
    #[error("Multiple keys with the same name:{0}")]
    MultipleKeys(String),

    #[error("Invalid key")]
    InvalidKey,
}

#[derive(Debug, Error)]
pub enum FunctionLibraryError<Key: Display> {
    #[error("Not found {0}")]
    NotFound(Key),

    #[error("Invalid key")]
    InvalidKey,
    #[error(transparent)]
    Other(#[from] eyre::Error),
}

pub trait FunctionLibrary<T, S> {
    type Key: Display;

    fn apply(&self, key: &Self::Key, tensor: T) -> Result<T, FunctionLibraryError<Self::Key>>;

    fn apply_scalar(
        &self,
        key: &Self::Key,
        scalar: S,
    ) -> Result<S, FunctionLibraryError<Self::Key>>;
}

pub trait Library<S> {
    type Key: Display;
    type Value: Clone;
    // type Structure: TensorStructure;

    fn key_for_structure(
        &self,
        structure: &PermutedStructure<S>,
    ) -> Result<Self::Key, LibraryError<Self::Key>>
    where
        S: TensorStructure;
    fn get<'a>(&'a self, key: &Self::Key) -> Result<Cow<'a, Self::Value>, LibraryError<Self::Key>>;
}

#[derive(Serialize, Deserialize)]
pub struct DummyLibrary<V, K = DummyKey> {
    key: PhantomData<K>,
    value: PhantomData<V>,
}

impl StructureLessDisplay for DummyKey {}

impl<V, K> Default for DummyLibrary<V, K> {
    fn default() -> Self {
        DummyLibrary {
            key: PhantomData,
            value: PhantomData,
        }
    }
}

#[allow(clippy::non_canonical_clone_impl)]
impl<V, K> Clone for DummyLibrary<V, K> {
    fn clone(&self) -> Self {
        DummyLibrary {
            key: PhantomData,
            value: PhantomData,
        }
    }
}

impl<V, K> Copy for DummyLibrary<V, K> {}

impl<V, K> DummyLibrary<V, K> {
    pub fn new() -> Self {
        DummyLibrary {
            key: PhantomData,
            value: PhantomData,
        }
    }
}

#[derive(Default, Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct DummyKey {}

impl Display for DummyKey {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "DUMMY")
    }
}

impl<K: Display + Clone, V: Clone, S> Library<S> for DummyLibrary<V, K> {
    type Key = K;
    type Value = PermutedStructure<DummyLibraryTensor<V>>;
    fn get<'a>(&'a self, key: &Self::Key) -> Result<Cow<'a, Self::Value>, LibraryError<Self::Key>> {
        Err(LibraryError::NotFound(key.clone()))
    }

    fn key_for_structure(
        &self,
        _structure: &PermutedStructure<S>,
    ) -> Result<Self::Key, LibraryError<Self::Key>>
    where
        S: TensorStructure,
    {
        Err(LibraryError::InvalidKey)
    }
}

#[derive(Default, Clone, Copy, Debug)]
pub struct DummyLibraryTensor<T> {
    with_indices: PhantomData<T>,
}

pub enum DummyIter<T> {
    None,
    Some(T),
}

impl<T> Iterator for DummyIter<T> {
    type Item = T;
    fn next(&mut self) -> Option<Self::Item> {
        None
    }
}

impl<T: TensorStructure> TensorStructure for DummyLibraryTensor<T> {
    type Slot = T::Slot;
    type Indexed = T::Indexed;

    fn is_fully_self_dual(&self) -> bool {
        false
    }
    fn reindex(
        self,
        _indices: &[<Self::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::Indexed>, StructureError> {
        unimplemented!()
    }
    fn dual(self) -> Self {
        unimplemented!()
    }

    fn external_structure_iter(&self) -> impl Iterator<Item = Self::Slot> {
        DummyIter::<T::Slot>::None
    }
    fn external_dims_iter(&self) -> impl Iterator<Item = Dimension> {
        DummyIter::<Dimension>::None
    }
    fn external_reps_iter(
        &self,
    ) -> impl Iterator<Item = Representation<<Self::Slot as IsAbstractSlot>::R>> {
        DummyIter::<Representation<<T::Slot as IsAbstractSlot>::R>>::None
    }

    fn external_indices_iter(&self) -> impl Iterator<Item = <Self::Slot as IsAbstractSlot>::Aind> {
        DummyIter::<<Self::Slot as IsAbstractSlot>::Aind>::None
    }
    fn get_aind(&self, _i: usize) -> Option<<Self::Slot as IsAbstractSlot>::Aind> {
        unimplemented!()
    }
    fn get_rep(&self, _i: usize) -> Option<Representation<<Self::Slot as IsAbstractSlot>::R>> {
        unimplemented!()
    }
    fn get_dim(&self, _i: usize) -> Option<Dimension> {
        unimplemented!()
    }
    fn get_slot(&self, _i: usize) -> Option<Self::Slot> {
        unimplemented!()
    }
    fn order(&self) -> usize {
        unimplemented!()
    }
}

impl<T: HasStructure + TensorStructure> HasStructure for DummyLibraryTensor<T> {
    type Store<S>
        = T::Store<S>
    where
        S: TensorStructure;

    type Scalar = T::Scalar;
    type ScalarRef<'a>
        = T::ScalarRef<'a>
    where
        Self: 'a;
    type Structure = T::Structure;

    fn map_structure<O: TensorStructure>(
        self,
        _f: impl Fn(Self::Structure) -> O,
    ) -> Self::Store<O> {
        unimplemented!()
    }

    fn map_structure_result<O: TensorStructure, Er>(
        self,
        _f: impl Fn(Self::Structure) -> Result<O, Er>,
    ) -> std::result::Result<Self::Store<O>, Er> {
        unimplemented!()
    }

    fn structure(&self) -> &Self::Structure {
        unimplemented!()
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        unimplemented!()
    }

    fn scalar(self) -> Option<Self::Scalar> {
        unimplemented!()
    }

    fn scalar_ref(&self) -> Option<Self::ScalarRef<'_>> {
        unimplemented!()
    }

    fn map_same_structure(self, _f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self {
        unimplemented!()
    }
}

impl<T: SetTensorData> SetTensorData for DummyLibraryTensor<T> {
    type SetData = T::SetData;

    fn set(&mut self, _indices: &[ConcreteIndex], _value: Self::SetData) -> Result<()> {
        unimplemented!()
    }

    fn set_flat(
        &mut self,
        _index: crate::structure::concrete_index::FlatIndex,
        _value: Self::SetData,
    ) -> Result<()> {
        unimplemented!()
    }
}

impl<T: HasStructure + TensorStructure> LibraryTensor for DummyLibraryTensor<T> {
    type WithIndices = T;
    type Data = ();

    fn empty(_key: Self::Structure, _zero: ()) -> Self {
        unimplemented!()
    }

    fn from_dense(_key: Self::Structure, _data: Vec<Self::Data>) -> Result<Self> {
        unimplemented!()
    }

    fn from_sparse(
        _key: Self::Structure,
        _data: impl IntoIterator<Item = (Vec<ConcreteIndex>, Self::Data)>,
        _zero: (),
    ) -> Result<Self> {
        unimplemented!()
    }

    fn with_indices(
        &self,
        _indices: &[<<<Self::WithIndices as HasStructure>::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::WithIndices>, StructureError> {
        unimplemented!()
    }
}

pub trait LibraryTensor: HasStructure + Sized + TensorStructure {
    type WithIndices: HasStructure;
    type Data;
    fn empty(key: Self::Structure, zero: Self::Data) -> Self;

    fn from_dense(key: Self::Structure, data: Vec<Self::Data>) -> Result<Self>;

    fn from_sparse(
        key: Self::Structure,
        data: impl IntoIterator<Item = (Vec<ConcreteIndex>, Self::Data)>,
        zero: Self::Data,
    ) -> Result<Self>;

    fn with_indices(
        &self,
        indices: &[<<<Self::WithIndices as HasStructure>::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::WithIndices>, StructureError>;
}

impl<D: Clone, S: TensorStructure + Clone> LibraryTensor for DataTensor<D, S> {
    type WithIndices = DataTensor<D, S::Indexed>;
    type Data = D;

    fn empty(key: S, zero: D) -> Self {
        DataTensor::Sparse(SparseTensor::empty(key, zero))
    }

    fn from_dense(key: S, data: Vec<Self::Data>) -> Result<Self> {
        Ok(DataTensor::Dense(DenseTensor::from_data(data, key)?))
    }

    fn from_sparse(
        key: S,
        data: impl IntoIterator<Item = (Vec<ConcreteIndex>, Self::Data)>,
        zero: Self::Data,
    ) -> Result<Self> {
        Ok(DataTensor::Sparse(SparseTensor::from_data(
            data, key, zero,
        )?))
    }

    fn with_indices(
        &self,
        indices: &[<<<Self::WithIndices as HasStructure>::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::WithIndices>, StructureError> {
        self.clone().reindex(indices)
        // let new_structure = self.structure().clone().reindex(indices)?;

        // Ok(match self {
        //     DataTensor::Dense(d) => DataTensor::Dense(DenseTensor {
        //         data: d.data.clone(),
        //         structure: new_structure,
        //     }),
        //     DataTensor::Sparse(s) => DataTensor::Sparse(SparseTensor {
        //         elements: s.elements.clone(),
        //         structure: new_structure,
        //     }),
        // })
    }
}

impl<D: Clone + Default, S: TensorStructure + Clone> LibraryTensor for RealOrComplexTensor<D, S> {
    type WithIndices = RealOrComplexTensor<D, S::Indexed>;
    type Data = RealOrComplex<D>;

    fn empty(key: S, zero: Self::Data) -> Self {
        RealOrComplexTensor::Real(DataTensor::Sparse(SparseTensor::empty(
            key,
            zero.to_complex().re,
        )))
    }

    fn from_dense(key: S, data: Vec<Self::Data>) -> Result<Self> {
        let complex_data = data.into_iter().map(|a| a.to_complex()).collect();
        let complex_tensor =
            <DataTensor<Complex<D>, S> as LibraryTensor>::from_dense(key, complex_data)?;

        Ok(RealOrComplexTensor::Complex(complex_tensor))
    }

    fn from_sparse(
        key: S,
        data: impl IntoIterator<Item = (Vec<ConcreteIndex>, Self::Data)>,
        zero: Self::Data,
    ) -> Result<Self> {
        let complex_tensor = <DataTensor<Complex<D>, S> as LibraryTensor>::from_sparse(
            key,
            data.into_iter()
                .map(|(indices, value)| (indices, value.to_complex())),
            zero.to_complex(),
        )?;

        Ok(RealOrComplexTensor::Complex(complex_tensor))
    }

    fn with_indices(
        &self,
        indices: &[<<<Self::WithIndices as HasStructure>::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::WithIndices>, StructureError> {
        match self {
            RealOrComplexTensor::Real(real_tensor) => {
                let new_real_tensor =
                    <DataTensor<D, S> as LibraryTensor>::with_indices(real_tensor, indices)?;
                Ok(PermutedStructure {
                    structure: RealOrComplexTensor::Real(new_real_tensor.structure),
                    rep_permutation: new_real_tensor.rep_permutation,
                    index_permutation: new_real_tensor.index_permutation,
                })
            }
            RealOrComplexTensor::Complex(complex_tensor) => {
                let new_complex_tensor =
                    <DataTensor<Complex<D>, S> as LibraryTensor>::with_indices(
                        complex_tensor,
                        indices,
                    )?;
                Ok(PermutedStructure {
                    structure: RealOrComplexTensor::Complex(new_complex_tensor.structure),
                    index_permutation: new_complex_tensor.index_permutation,
                    rep_permutation: new_complex_tensor.rep_permutation,
                })
            }
        }
    }
}

#[cfg(feature = "shadowing")]
pub mod symbolic;

pub trait TensorLibraryData: Neg<Output = Self> {
    fn one() -> Self;
    fn minus_one() -> Self;
    fn zero() -> Self;
}
macro_rules! impl_tensor_library_data {
    ( $target_type:ty, $zero_val:expr, $one_val:expr ) => {
        // Assuming TensorLibraryData is defined in spenso::tensor_library
        // If the macro is defined *inside* tensor_library.rs, `TensorLibraryData` might suffice.
        // Using `$crate::spenso::tensor_library::TensorLibraryData` is more robust
        // if the macro can be called from other modules within the spenso crate.
        // Adjust path as needed based on your crate structure.
        impl TensorLibraryData for $target_type {
            /// Returns the additive identity element (zero) of this type.
            #[inline]
            fn zero() -> Self {
                $zero_val
            }

            /// Returns the additive inverse element (minus one) of this type.
            #[inline]
            fn minus_one() -> Self {
                -Self::one()
            }

            /// Returns the multiplicative identity element (one) of this type.
            #[inline]
            fn one() -> Self {
                $one_val
            }
        }
    };
}

impl_tensor_library_data!(f32, 0.0, 1.0);
impl_tensor_library_data!(f64, 0.0, 1.0);
impl_tensor_library_data!(i32, 0, 1);
impl_tensor_library_data!(i64, 0, 1);
#[cfg(feature = "shadowing")]
impl_tensor_library_data!(
    symbolica::atom::Atom,
    symbolica::atom::Atom::Zero,
    symbolica::atom::Atom::num(1)
);

impl<T: TensorLibraryData> TensorLibraryData for RealOrComplex<T> {
    fn zero() -> Self {
        RealOrComplex::Real(T::zero())
    }

    fn minus_one() -> Self {
        RealOrComplex::Real(T::minus_one())
    }

    fn one() -> Self {
        RealOrComplex::Real(T::one())
    }
}

#[cfg(test)]
mod test {}
