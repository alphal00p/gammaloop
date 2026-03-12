use std::fmt::{Debug, Display};

use crate::{
    iterators::IteratorEnum,
    structure::{
        IndexLess, PermutedStructure, SlotIndex,
        concrete_index::ConcreteIndex,
        permuted::PermuteTensor,
        representation::RepName,
        slot::{AbsInd, Slot},
    },
    tensors::data::{SparseTensor, StorageTensor},
};
use enum_try_as_inner::EnumTryAsInner;
use eyre::{Result, eyre};
#[cfg(feature = "shadowing")]
use symbolica::{atom::Atom, evaluate::FunctionMap};

use crate::structure::StructureError;
use crate::structure::dimension::Dimension;
use crate::structure::representation::Representation;
use crate::structure::slot::IsAbstractSlot;
use delegate::delegate;

#[cfg(feature = "shadowing")]
use crate::{
    shadowing::{
        ShadowMapping, Shadowable,
        symbolica_utils::{IntoArgs, IntoSymbol},
    },
    structure::{ToSymbolic, slot::ParseableAind},
    tensors::{
        data::{DataIterator, DenseTensor},
        parametric::TensorCoefficient,
    },
};

use crate::{
    algebra::algebraic_traits::{IsZero, RefZero},
    algebra::complex::{Complex, RealOrComplex, RealOrComplexMut, RealOrComplexRef},
    algebra::upgrading_arithmetic::{FallibleAddAssign, FallibleMul, FallibleSubAssign},
    contraction::{Contract, ContractableWith, ContractionError, Trace},
    iterators::IteratableTensor,
    structure::{
        CastStructure, HasName, HasStructure, ScalarStructure, ScalarTensor, StructureContract,
        TensorStructure, TracksCount,
        concrete_index::{ExpandedIndex, FlatIndex},
    },
    tensors::data::{DataTensor, GetTensorData, HasTensorData, SetTensorData, SparseOrDense},
};

#[derive(
    Clone,
    Debug,
    EnumTryAsInner,
    bincode_trait_derive::Encode,
    bincode_trait_derive::Decode,
    PartialEq,
    Eq,
)]
#[cfg_attr(
feature = "shadowing",
trait_decode(trait = symbolica::state::HasStateMap),
)]
#[derive_err(Debug)]
pub enum RealOrComplexTensor<T, S: TensorStructure> {
    Real(DataTensor<T, S>),
    Complex(DataTensor<Complex<T>, S>),
}

impl<T: RefZero, S: TensorStructure> RealOrComplexTensor<T, S> {
    pub fn ref_zero(&self) -> T {
        match self {
            RealOrComplexTensor::Real(a) => a.ref_zero(),
            RealOrComplexTensor::Complex(a) => a.ref_zero().re,
        }
    }
}

impl<T, S: TensorStructure> crate::network::Ref for RealOrComplexTensor<T, S> {
    type Ref<'a>
        = &'a RealOrComplexTensor<T, S>
    where
        Self: 'a;

    fn refer(&self) -> Self::Ref<'_> {
        self
    }
}

impl<T: RefZero, S: TensorStructure + ScalarStructure + Clone> RealOrComplexTensor<T, S> {
    pub fn to_complex(&mut self) {
        if self.is_real() {
            let zero = self.ref_zero();

            let old = std::mem::replace(
                self,
                RealOrComplexTensor::Real(DataTensor::Sparse(SparseTensor::empty(
                    S::scalar_structure(),
                    zero,
                ))),
            );

            if let RealOrComplexTensor::Real(r) = old {
                *self = RealOrComplexTensor::Complex(r.map_data(|a| Complex::new_re(a)));
            }
        }
    }
}

impl<T: Clone, S: TensorStructure> SetTensorData for RealOrComplexTensor<T, S> {
    type SetData = RealOrComplex<T>;

    fn set(
        &mut self,
        indices: &[crate::structure::concrete_index::ConcreteIndex],
        value: Self::SetData,
    ) -> Result<()> {
        match self {
            RealOrComplexTensor::Real(d) => d.set(
                indices,
                value.try_into_real().map_err(|r| eyre!(r.to_string()))?,
            )?,
            RealOrComplexTensor::Complex(d) => d.set(
                indices,
                value.try_into_complex().map_err(|r| eyre!(r.to_string()))?,
            )?,
        }
        Ok(())
    }

    fn set_flat(&mut self, index: FlatIndex, value: Self::SetData) -> Result<()> {
        match self {
            RealOrComplexTensor::Real(d) => d.set_flat(
                index,
                value.try_into_real().map_err(|r| eyre!(r.to_string()))?,
            )?,
            RealOrComplexTensor::Complex(d) => d.set_flat(
                index,
                value.try_into_complex().map_err(|r| eyre!(r.to_string()))?,
            )?,
        }
        Ok(())
    }
}

impl<T: Clone, S: TensorStructure> GetTensorData for RealOrComplexTensor<T, S> {
    type GetDataRef<'a>
        = RealOrComplexRef<'a, T>
    where
        Self: 'a;

    type GetDataRefMut<'a>
        = RealOrComplexMut<'a, T>
    where
        Self: 'a;

    type GetDataOwned = RealOrComplex<T>;

    fn get_ref<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<Self::GetDataRef<'_>> {
        match self {
            RealOrComplexTensor::Real(d) => Ok(RealOrComplexRef::Real(d.get_ref(indices)?)),
            RealOrComplexTensor::Complex(d) => Ok(RealOrComplexRef::Complex(d.get_ref(indices)?)),
        }
    }

    fn get_ref_linear(&self, index: FlatIndex) -> Option<Self::GetDataRef<'_>> {
        match self {
            RealOrComplexTensor::Real(d) => d.get_ref_linear(index).map(RealOrComplexRef::Real),
            RealOrComplexTensor::Complex(d) => {
                d.get_ref_linear(index).map(RealOrComplexRef::Complex)
            }
        }
    }

    fn get_mut_linear(&mut self, index: FlatIndex) -> Option<Self::GetDataRefMut<'_>> {
        match self {
            RealOrComplexTensor::Real(d) => d.get_mut_linear(index).map(RealOrComplexMut::Real),
            RealOrComplexTensor::Complex(d) => {
                d.get_mut_linear(index).map(RealOrComplexMut::Complex)
            }
        }
    }

    fn get_owned<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<Self::GetDataOwned>
    where
        Self::GetDataOwned: Clone,
    {
        match self {
            RealOrComplexTensor::Real(d) => Ok(RealOrComplex::Real(d.get_owned(indices)?)),
            RealOrComplexTensor::Complex(d) => Ok(RealOrComplex::Complex(d.get_owned(indices)?)),
        }
    }

    fn get_owned_linear(&self, index: FlatIndex) -> Option<Self::GetDataOwned>
    where
        Self::GetDataOwned: Clone,
    {
        match self {
            RealOrComplexTensor::Real(d) => Some(RealOrComplex::Real(d.get_owned_linear(index)?)),
            RealOrComplexTensor::Complex(d) => {
                Some(RealOrComplex::Complex(d.get_owned_linear(index)?))
            }
        }
    }
}

impl<T: Clone, S: TensorStructure + Clone> HasTensorData for RealOrComplexTensor<T, S> {
    type Data = RealOrComplex<T>;

    fn data(&self) -> Vec<Self::Data> {
        match self {
            RealOrComplexTensor::Real(d) => d.data().into_iter().map(RealOrComplex::Real).collect(),
            RealOrComplexTensor::Complex(d) => {
                d.data().into_iter().map(RealOrComplex::Complex).collect()
            }
        }
    }

    fn hashmap(&self) -> indexmap::IndexMap<ExpandedIndex, Self::Data> {
        match self {
            RealOrComplexTensor::Real(d) => d
                .hashmap()
                .into_iter()
                .map(|(k, v)| (k, RealOrComplex::Real(v)))
                .collect(),
            RealOrComplexTensor::Complex(d) => d
                .hashmap()
                .into_iter()
                .map(|(k, v)| (k, RealOrComplex::Complex(v)))
                .collect(),
        }
    }

    fn indices(&self) -> Vec<ExpandedIndex> {
        match self {
            RealOrComplexTensor::Real(d) => d.indices(),
            RealOrComplexTensor::Complex(d) => d.indices(),
        }
    }

    #[cfg(feature = "shadowing")]
    fn symhashmap(
        &self,
        name: symbolica::atom::Symbol,
        args: &[Atom],
    ) -> std::collections::HashMap<Atom, Self::Data> {
        match self {
            RealOrComplexTensor::Real(d) => d
                .symhashmap(name, args)
                .into_iter()
                .map(|(k, v)| (k, RealOrComplex::Real(v)))
                .collect(),
            RealOrComplexTensor::Complex(d) => d
                .symhashmap(name, args)
                .into_iter()
                .map(|(k, v)| (k, RealOrComplex::Complex(v)))
                .collect(),
        }
    }
}

// use crate::data::StorageTensor;

impl<T: Display, S: TensorStructure> Display for RealOrComplexTensor<T, S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RealOrComplexTensor::Real(d) => d.fmt(f),
            RealOrComplexTensor::Complex(d) => d.fmt(f),
        }
    }
}

impl<T: Clone, S: TensorStructure, O: From<S> + TensorStructure>
    CastStructure<RealOrComplexTensor<T, O>> for RealOrComplexTensor<T, S>
{
    fn cast_structure(self) -> RealOrComplexTensor<T, O> {
        match self {
            RealOrComplexTensor::Real(d) => RealOrComplexTensor::Real(d.cast_structure()),
            RealOrComplexTensor::Complex(d) => RealOrComplexTensor::Complex(d.cast_structure()),
        }
    }
}

#[cfg(feature = "shadowing")]
impl<T: Clone, S: TensorStructure> Shadowable for RealOrComplexTensor<T, S>
where
    S: HasName + Clone,
    S::Name: IntoSymbol,
    S::Args: IntoArgs,
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn shadow<C>(
        &self,
        index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> C,
    ) -> Result<DenseTensor<Atom, Self::Structure>>
    where
        C: TensorCoefficient,
    {
        match self {
            RealOrComplexTensor::Real(r) => r.shadow(index_to_atom),
            RealOrComplexTensor::Complex(r) => Ok(r
                .structure()
                .clone()
                .to_dense_labeled_complex(index_to_atom)?),
        }
        // Some(self.structure().clone().to_dense_labeled(index_to_atom))
    }
}

impl<T: Default + Clone + PartialEq, S: TensorStructure + Clone> SparseOrDense
    for RealOrComplexTensor<T, S>
{
    fn to_sparse(self) -> Self {
        match self {
            RealOrComplexTensor::Real(d) => RealOrComplexTensor::Real(d.to_sparse()),
            RealOrComplexTensor::Complex(d) => RealOrComplexTensor::Complex(d.to_sparse()),
        }
    }

    fn to_dense_mut(&mut self) {
        match self {
            RealOrComplexTensor::Real(d) => d.to_dense_mut(),
            RealOrComplexTensor::Complex(d) => d.to_dense_mut(),
        }
    }
    fn to_sparse_mut(&mut self) {
        match self {
            RealOrComplexTensor::Real(d) => d.to_sparse_mut(),
            RealOrComplexTensor::Complex(d) => d.to_sparse_mut(),
        }
    }

    fn to_dense(self) -> Self {
        match self {
            RealOrComplexTensor::Real(d) => RealOrComplexTensor::Real(d.to_dense()),
            RealOrComplexTensor::Complex(d) => RealOrComplexTensor::Complex(d.to_dense()),
        }
    }
}

#[cfg(feature = "shadowing")]
impl<T: Clone + RefZero, S: TensorStructure, R> ShadowMapping<R> for RealOrComplexTensor<T, S>
where
    S: HasName + Clone,
    S::Name: IntoSymbol,
    S::Args: IntoArgs,
    R: From<T>,
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn append_map<C>(
        &self,
        fn_map: &mut FunctionMap<R>,
        index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> C,
    ) where
        C: TensorCoefficient,
    {
        match self {
            RealOrComplexTensor::Real(c) => c.append_map(fn_map, index_to_atom),
            RealOrComplexTensor::Complex(p) => match p {
                DataTensor::Dense(d) => {
                    for (i, c) in d.flat_iter() {
                        let labeled_coef_re =
                            index_to_atom(self.structure(), i).to_atom_re().unwrap();
                        let labeled_coef_im =
                            index_to_atom(self.structure(), i).to_atom_im().unwrap();
                        fn_map.add_constant(labeled_coef_re.clone(), c.re.clone().into());
                        fn_map.add_constant(labeled_coef_im.clone(), c.re.clone().into());
                    }
                }
                DataTensor::Sparse(d) => {
                    for (i, c) in d.flat_iter() {
                        let labeled_coef_re =
                            index_to_atom(self.structure(), i).to_atom_re().unwrap();
                        let labeled_coef_im =
                            index_to_atom(self.structure(), i).to_atom_im().unwrap();
                        fn_map.add_constant(labeled_coef_re.clone(), c.re.clone().into());
                        fn_map.add_constant(labeled_coef_im.clone(), c.re.clone().into());
                    }
                }
            }, // p.append_map(fn_map, index_to_atom),
        }
    }
}

impl<T: Clone, Aind: AbsInd, S: Clone + Into<IndexLess<R, Aind>>, R: RepName<Dual = R>>
    PermuteTensor for RealOrComplexTensor<T, S>
where
    S: TensorStructure<Slot = Slot<R, Aind>> + PermuteTensor<IdSlot = Slot<R, Aind>, Id = S>,
{
    type Id = RealOrComplexTensor<T, S>;
    type IdSlot = (T, Slot<R, Aind>);
    type Permuted = RealOrComplexTensor<T, S>;

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
        RealOrComplexTensor::Real(DataTensor::id(i, j))
    }

    fn permute_inds(self, permutation: &linnet::permutation::Permutation) -> Self::Permuted {
        match self {
            RealOrComplexTensor::Real(d) => RealOrComplexTensor::Real(d.permute_inds(permutation)),
            RealOrComplexTensor::Complex(s) => {
                RealOrComplexTensor::Complex(s.permute_inds(permutation))
            }
        }
    }

    fn permute_reps(self, rep_perm: &linnet::permutation::Permutation) -> Self::Permuted {
        match self {
            RealOrComplexTensor::Real(d) => RealOrComplexTensor::Real(d.permute_reps(rep_perm)),
            RealOrComplexTensor::Complex(s) => {
                RealOrComplexTensor::Complex(s.permute_reps(rep_perm))
            }
        }
    }
}

impl<T: Clone, S> TensorStructure for RealOrComplexTensor<T, S>
where
    S: TensorStructure,
{
    // type R = <T::Structure as TensorStructure>::R;
    type Indexed = RealOrComplexTensor<T, S::Indexed>;
    type Slot = S::Slot;

    fn reindex(
        self,
        indices: &[<Self::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::Indexed>, StructureError> {
        Ok(match self {
            RealOrComplexTensor::Complex(d) => {
                let res = d.reindex(indices)?;
                PermutedStructure {
                    structure: RealOrComplexTensor::Complex(res.structure),
                    rep_permutation: res.rep_permutation,
                    index_permutation: res.index_permutation,
                }
            }
            RealOrComplexTensor::Real(d) => {
                let res = d.reindex(indices)?;
                PermutedStructure {
                    structure: RealOrComplexTensor::Real(res.structure),
                    rep_permutation: res.rep_permutation,
                    index_permutation: res.index_permutation,
                }
            }
        })
    }

    fn dual(self) -> Self {
        self.map_same_structure(|s| s.dual())
    }

    delegate! {
        to self.structure() {
            fn is_fully_self_dual(&self)-> bool;
            fn external_reps_iter(&self)-> impl Iterator<Item = Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn external_indices_iter(&self)-> impl Iterator<Item = <Self::Slot as IsAbstractSlot>::Aind>;
            fn external_dims_iter(&self)-> impl Iterator<Item = Dimension>;
            fn external_structure_iter(&self)-> impl Iterator<Item = Self::Slot>;
            fn get_slot(&self, i: impl Into<SlotIndex>)-> Option<Self::Slot>;
            fn get_rep(&self, i: impl Into<SlotIndex>)-> Option<Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn get_dim(&self, i: impl Into<SlotIndex>)-> Option<Dimension>;
            fn get_aind(&self, i: impl Into<SlotIndex>)-> Option<<Self::Slot as IsAbstractSlot>::Aind>;
            fn order(&self)-> usize;
        }
    }
}

impl<T: Clone, S: TensorStructure> HasStructure for RealOrComplexTensor<T, S> {
    type Scalar = RealOrComplex<T>;
    type ScalarRef<'a>
        = RealOrComplexRef<'a, T>
    where
        Self: 'a;
    type Structure = S;
    type Store<U>
        = RealOrComplexTensor<T, U>
    where
        U: TensorStructure;

    fn map_structure<S2: TensorStructure>(self, f: impl Fn(S) -> S2) -> RealOrComplexTensor<T, S2> {
        match self {
            RealOrComplexTensor::Real(d) => RealOrComplexTensor::Real(d.map_structure(f)),
            RealOrComplexTensor::Complex(d) => RealOrComplexTensor::Complex(d.map_structure(f)),
        }
    }

    fn map_structure_result<S2: TensorStructure, E>(
        self,
        f: impl Fn(S) -> Result<S2, E>,
    ) -> Result<RealOrComplexTensor<T, S2>, E> {
        Ok(match self {
            RealOrComplexTensor::Real(d) => RealOrComplexTensor::Real(d.map_structure_result(f)?),
            RealOrComplexTensor::Complex(d) => {
                RealOrComplexTensor::Complex(d.map_structure_result(f)?)
            }
        })
    }

    fn structure(&self) -> &Self::Structure {
        match self {
            RealOrComplexTensor::Real(r) => r.structure(),
            RealOrComplexTensor::Complex(r) => r.structure(),
        }
    }

    fn map_same_structure(self, f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self {
        match self {
            RealOrComplexTensor::Real(d) => RealOrComplexTensor::Real(d.map_same_structure(f)),
            RealOrComplexTensor::Complex(d) => {
                RealOrComplexTensor::Complex(d.map_same_structure(f))
            }
        }
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        match self {
            RealOrComplexTensor::Real(r) => r.mut_structure(),
            RealOrComplexTensor::Complex(r) => r.mut_structure(),
        }
    }

    fn scalar(self) -> Option<Self::Scalar> {
        match self {
            RealOrComplexTensor::Real(r) => r.scalar().map(|x| RealOrComplex::Real(x)),
            RealOrComplexTensor::Complex(r) => r.scalar().map(|x| RealOrComplex::Complex(x)),
        }
    }

    fn scalar_ref(&self) -> Option<Self::ScalarRef<'_>> {
        match self {
            RealOrComplexTensor::Real(r) => r.scalar_ref().map(|x| RealOrComplexRef::Real(x)),
            RealOrComplexTensor::Complex(r) => r.scalar_ref().map(|x| RealOrComplexRef::Complex(x)),
        }
    }
}

impl<T: Clone, S: TensorStructure + ScalarStructure> ScalarTensor for RealOrComplexTensor<T, S> {
    fn new_scalar(scalar: Self::Scalar) -> Self {
        match scalar {
            RealOrComplex::Real(r) => RealOrComplexTensor::Real(DataTensor::new_scalar(r)),
            RealOrComplex::Complex(r) => RealOrComplexTensor::Complex(DataTensor::new_scalar(r)),
        }
    }
}

impl<T, S> TracksCount for RealOrComplexTensor<T, S>
where
    S: TensorStructure + TracksCount,
    T: Clone,
{
    fn contractions_num(&self) -> usize {
        match self {
            RealOrComplexTensor::Real(r) => r.contractions_num(),
            RealOrComplexTensor::Complex(r) => r.contractions_num(),
        }
    }
}

impl<T, S> HasName for RealOrComplexTensor<T, S>
where
    S: TensorStructure + HasName,
    T: Clone,
{
    type Args = S::Args;
    type Name = S::Name;

    fn name(&self) -> Option<S::Name> {
        match self {
            RealOrComplexTensor::Real(r) => r.name(),
            RealOrComplexTensor::Complex(r) => r.name(),
        }
    }

    fn set_name(&mut self, name: Self::Name) {
        match self {
            RealOrComplexTensor::Real(r) => r.set_name(name),
            RealOrComplexTensor::Complex(r) => r.set_name(name),
        }
    }

    fn args(&self) -> Option<S::Args> {
        match self {
            RealOrComplexTensor::Real(r) => r.args(),
            RealOrComplexTensor::Complex(r) => r.args(),
        }
    }
}

impl<T: Clone, S: TensorStructure> IteratableTensor for RealOrComplexTensor<T, S> {
    type Data<'a>
        = RealOrComplexRef<'a, T>
    where
        Self: 'a;

    fn iter_expanded(&self) -> impl Iterator<Item = (ExpandedIndex, Self::Data<'_>)> {
        match self {
            RealOrComplexTensor::Real(x) => IteratorEnum::A(
                x.iter_expanded()
                    .map(|(i, x)| (i, RealOrComplexRef::Real(x))),
            ),
            RealOrComplexTensor::Complex(x) => IteratorEnum::B(
                x.iter_expanded()
                    .map(|(i, x)| (i, RealOrComplexRef::Complex(x))),
            ),
        }
    }

    fn iter_flat(&self) -> impl Iterator<Item = (FlatIndex, Self::Data<'_>)> {
        match self {
            RealOrComplexTensor::Real(x) => {
                IteratorEnum::A(x.iter_flat().map(|(i, x)| (i, RealOrComplexRef::Real(x))))
            }
            RealOrComplexTensor::Complex(x) => IteratorEnum::B(
                x.iter_flat()
                    .map(|(i, x)| (i, RealOrComplexRef::Complex(x))),
            ),
        }
    }
}

impl<S, T> Trace for RealOrComplexTensor<T, S>
where
    S: TensorStructure + Clone + StructureContract,
    T: ContractableWith<T, Out = T>
        + Clone
        + FallibleMul<Output = T>
        + FallibleAddAssign<T>
        + FallibleSubAssign<T>
        + RefZero
        + IsZero,
    Complex<T>: ContractableWith<Complex<T>, Out = Complex<T>>
        + Clone
        + FallibleMul<Output = Complex<T>>
        + FallibleAddAssign<Complex<T>>
        + FallibleSubAssign<Complex<T>>
        + RefZero
        + IsZero,
{
    fn internal_contract(&self) -> Self {
        match self {
            RealOrComplexTensor::Real(x) => RealOrComplexTensor::Real(x.internal_contract()),
            RealOrComplexTensor::Complex(x) => RealOrComplexTensor::Complex(x.internal_contract()),
        }
    }
}

impl<S, T> Contract<RealOrComplexTensor<T, S>> for RealOrComplexTensor<T, S>
where
    S: TensorStructure + Clone + StructureContract,
    T: ContractableWith<T, Out = T>
        + ContractableWith<Complex<T>, Out = Complex<T>>
        + Clone
        + FallibleMul<Output = T>
        + FallibleAddAssign<T>
        + FallibleSubAssign<T>
        + RefZero
        + IsZero,
    Complex<T>: ContractableWith<T, Out = Complex<T>>
        + ContractableWith<Complex<T>, Out = Complex<T>>
        + Clone
        + FallibleMul<Output = Complex<T>>
        + FallibleAddAssign<Complex<T>>
        + FallibleSubAssign<Complex<T>>
        + RefZero
        + IsZero,
{
    type LCM = RealOrComplexTensor<T, S>;
    fn contract(&self, other: &RealOrComplexTensor<T, S>) -> Result<Self::LCM, ContractionError> {
        match (self, other) {
            (RealOrComplexTensor::Real(s), RealOrComplexTensor::Real(o)) => {
                Ok(RealOrComplexTensor::Real(s.contract(o)?))
            }
            (RealOrComplexTensor::Real(s), RealOrComplexTensor::Complex(o)) => {
                Ok(RealOrComplexTensor::Complex(s.contract(o)?))
            }
            (RealOrComplexTensor::Complex(s), RealOrComplexTensor::Real(o)) => {
                Ok(RealOrComplexTensor::Complex(s.contract(o)?))
            }
            (RealOrComplexTensor::Complex(s), RealOrComplexTensor::Complex(o)) => {
                Ok(RealOrComplexTensor::Complex(s.contract(o)?))
            }
        }
    }
}
