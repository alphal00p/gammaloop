extern crate derive_more;

use std::{
    fmt::{Debug, Display},
    io::Cursor,
    path::Path,
};

use crate::structure::{IndexLess, SlotIndex, slot::IsAbstractSlot};
use crate::structure::{StructureError, slot::AbsInd};
use crate::structure::{
    permuted::PermuteTensor, representation::Representation, slot::ParseableAind,
};
use crate::{
    algebra::algebraic_traits::RefOne,
    structure::{dimension::Dimension, representation::RepName},
};
use crate::{
    algebra::complex::symbolica_traits::CompiledComplexEvaluatorSpenso,
    structure::PermutedStructure,
};
use ahash::HashMap;
use delegate::delegate;

use eyre::Result;
use eyre::eyre;

use atomcore::{ReplaceBuilderGeneric, TensorAtomMaps};

use enum_try_as_inner::EnumTryAsInner;
use log::trace;
use serde::{Deserialize, Serialize};
use to_param::ToAtom;

use crate::{
    algebra::algebraic_traits::{IsZero, RefZero},
    algebra::complex::{Complex, RealOrComplex},
    algebra::upgrading_arithmetic::{
        FallibleAddAssign, FallibleMul, FallibleSubAssign, TrySmallestUpgrade,
    },
    contraction::{Contract, ContractableWith, ContractionError, Trace},
    iterators::{IteratableTensor, IteratorEnum},
    shadowing::symbolica_utils::{IntoArgs, IntoSymbol},
    shadowing::{ShadowMapping, Shadowable},
    structure::{
        CastStructure, HasName, HasStructure, NamedStructure, OrderedStructure, ScalarStructure,
        ScalarTensor, StructureContract, TensorStructure, TracksCount,
        concrete_index::{ConcreteIndex, DualConciousExpandedIndex, ExpandedIndex, FlatIndex},
        slot::Slot,
    },
    tensors::complex::RealOrComplexTensor,
    tensors::data::{
        DataIterator, DataTensor, DenseTensor, GetTensorData, HasTensorData, SetTensorData,
        SparseOrDense, SparseTensor, StorageTensor,
    },
};
use bincode::{Decode, Encode};

use symbolica::{
    atom::{Atom, AtomCore, AtomView, FunctionBuilder, KeyLookup, Symbol},
    coefficient::Coefficient,
    domains::{
        InternalOrdering,
        float::{FloatLike, Real, SingleFloat},
        rational::Rational,
    },
    evaluate::{
        CompileOptions, CompiledCode, CompiledComplexEvaluator, CompiledNumber, EvalTree,
        EvaluationFn, ExportNumber, ExportSettings, ExportedCode, Expression, ExpressionEvaluator,
        FunctionMap, OptimizationSettings,
    },
    id::Pattern,
    state::State,
    symbol,
    utils::BorrowedOrOwned,
};

use std::hash::Hash;

use symbolica::domains::float::Complex as SymComplex;

// impl RefZero for Atom {
//     fn ref_zero(&self) -> Self {
//         Atom::num(0)
//     }
// }

pub trait TensorCoefficient: Display {
    fn cooked_name(&self) -> Option<String>;
    fn name(&self) -> Option<Symbol>;
    fn tags(&self) -> Vec<Atom>;
    fn to_atom(&self) -> Option<Atom>;
    fn to_atom_re(&self) -> Option<Atom>;
    fn to_atom_im(&self) -> Option<Atom>;
    fn add_tagged_function<T>(
        &self,
        fn_map: &mut FunctionMap<T>,
        body: Atom,
    ) -> Result<(), String> {
        let (name, cooked_name) = self
            .name()
            .zip(self.cooked_name())
            .ok_or(format!("unnamed {}", self))?;

        fn_map.add_tagged_function::<Symbol>(name, self.tags(), cooked_name, vec![], body)
    }
}

#[derive(Debug)]
pub struct FlatCoefficent<Args: IntoArgs> {
    pub name: Option<Symbol>,
    pub index: FlatIndex,
    pub args: Option<Args>,
}

impl<Arg: IntoArgs> Display for FlatCoefficent<Arg> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(name) = self.name {
            write!(f, "{}", name)?
        }
        write!(f, "(")?;
        if let Some(ref args) = self.args {
            let args: Vec<String> = args.ref_into_args().map(|s| s.to_string()).collect();
            write!(f, "{},", args.join(","))?
        }

        write!(f, "{})", self.index)?;
        Result::Ok(())
    }
}

impl<Args: IntoArgs> TensorCoefficient for FlatCoefficent<Args> {
    fn name(&self) -> Option<Symbol> {
        self.name
    }

    fn cooked_name(&self) -> Option<String> {
        let mut name = self.name?.to_string();
        if let Some(ref args) = self.args {
            name += args.cooked_name().as_str();
        }
        Some(name)
    }

    fn tags(&self) -> Vec<Atom> {
        let mut tags: Vec<Atom> = if let Some(ref args) = self.args {
            args.args()
        } else {
            vec![]
        };
        tags.push(Atom::from(self.index));
        tags
    }

    fn to_atom(&self) -> Option<Atom> {
        let mut fn_builder = FunctionBuilder::new(self.name?);
        if let Some(ref args) = self.args {
            for arg in args.ref_into_args() {
                fn_builder = fn_builder.add_arg(arg.as_view());
            }
        }
        fn_builder = fn_builder.add_arg(Atom::from(self.index).as_view());
        Some(fn_builder.finish())
    }

    fn to_atom_re(&self) -> Option<Atom> {
        let name = symbol!(self.name?.to_string() + "_re");

        let mut fn_builder = FunctionBuilder::new(name);
        if let Some(ref args) = self.args {
            for arg in args.ref_into_args() {
                fn_builder = fn_builder.add_arg(arg.as_view());
            }
        }
        fn_builder = fn_builder.add_arg(Atom::from(self.index).as_view());
        Some(fn_builder.finish())
    }

    fn to_atom_im(&self) -> Option<Atom> {
        let name = symbol!(self.name?.to_string() + "_im");

        let mut fn_builder = FunctionBuilder::new(name);
        if let Some(ref args) = self.args {
            for arg in args.ref_into_args() {
                fn_builder = fn_builder.add_arg(arg.as_view());
            }
        }
        fn_builder = fn_builder.add_arg(Atom::from(self.index).as_view());
        Some(fn_builder.finish())
    }
}

pub struct ExpandedCoefficent<Args: IntoArgs> {
    pub name: Option<Symbol>,
    pub index: DualConciousExpandedIndex,
    pub args: Option<Args>,
}

impl<Arg: IntoArgs> Display for ExpandedCoefficent<Arg> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(name) = self.name {
            write!(f, "{}", name)?
        }
        write!(f, "(")?;
        if let Some(ref args) = self.args {
            let args: Vec<String> = args.ref_into_args().map(|s| s.to_string()).collect();
            write!(f, "{},", args.join(","))?
        }
        write!(f, "{})", self.index)?;
        Result::Ok(())
    }
}

impl<Args: IntoArgs> TensorCoefficient for ExpandedCoefficent<Args> {
    fn name(&self) -> Option<Symbol> {
        self.name
    }
    fn cooked_name(&self) -> Option<String> {
        let mut name = self.name?.to_string();
        if let Some(ref args) = self.args {
            name += args.cooked_name().as_str();
        }
        Some(name)
    }

    fn tags(&self) -> Vec<Atom> {
        let mut tags: Vec<Atom> = if let Some(ref args) = self.args {
            args.args()
        } else {
            vec![]
        };
        tags.push(Atom::from(self.index.clone()));
        tags
    }

    fn to_atom(&self) -> Option<Atom> {
        let mut fn_builder = FunctionBuilder::new(self.name?);
        if let Some(ref args) = self.args {
            for arg in args.ref_into_args() {
                fn_builder = fn_builder.add_arg(arg.as_view());
            }
        }
        fn_builder = fn_builder.add_arg(Atom::from(self.index.clone()).as_view());
        Some(fn_builder.finish())
    }
    fn to_atom_re(&self) -> Option<Atom> {
        let name = symbol!(self.name?.to_string() + "_re");

        let mut fn_builder = FunctionBuilder::new(name);
        if let Some(ref args) = self.args {
            for arg in args.ref_into_args() {
                fn_builder = fn_builder.add_arg(arg.as_view());
            }
        }
        fn_builder = fn_builder.add_arg(Atom::from(self.index.clone()).as_view());
        Some(fn_builder.finish())
    }

    fn to_atom_im(&self) -> Option<Atom> {
        let name = symbol!(self.name?.to_string() + "_im");

        let mut fn_builder = FunctionBuilder::new(name);
        if let Some(ref args) = self.args {
            for arg in args.ref_into_args() {
                fn_builder = fn_builder.add_arg(arg.as_view());
            }
        }
        fn_builder = fn_builder.add_arg(Atom::from(self.index.clone()).as_view());
        Some(fn_builder.finish())
    }
}

// impl<'a> TryFrom<FunView<'a>> for DenseTensor<Atom, NamedStructure<Symbol, Vec<Atom>>> {
//     type Error = Error;

//     fn try_from(f: FunView<'a>) -> Result<Self> {
//         let mut structure: Vec<Slot<LibraryRep>> = vec![];
//         let f_id = f.get_symbol();
//         let mut args = vec![];

//         for arg in f.iter() {
//             if let Ok(arg) = arg.try_into() {
//                 structure.push(arg);
//             } else {
//                 args.push(arg.to_owned());
//             }
//         }
//         let s = NamedStructure::from_iter(structure, f_id, Some(args));
//         s.to_dense_expanded_labels()
//     }
// }

#[derive(
    Clone, Debug, PartialEq, Eq, Hash, bincode_trait_derive::Encode, bincode_trait_derive::Decode,
)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub struct ParamTensor<S = OrderedStructure> {
    pub tensor: DataTensor<Atom, S>,
    pub param_type: ParamOrComposite,
    // Param(DataTensor<Atom, S>),
    // // Concrete(DataTensor<T, S>),
    // Composite(DataTensor<Atom, S>),
}

pub mod add_assign;
pub mod div;
pub mod mul_assign;
pub mod neg;
pub mod scalar_mul;

impl<Aind: AbsInd, S: Clone + Into<IndexLess<R, Aind>>, R: RepName<Dual = R>> PermuteTensor
    for ParamTensor<S>
where
    S: TensorStructure<Slot = Slot<R, Aind>> + PermuteTensor<IdSlot = Slot<R, Aind>, Id = S>,
{
    type Id = ParamTensor<S>;
    type IdSlot = Slot<R, Aind>;
    type Permuted = ParamTensor<S>;

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
        ParamTensor {
            tensor: DataTensor::id((Atom::Zero, i), (Atom::num(1), j)),
            param_type: ParamOrComposite::Composite,
        }
    }

    fn permute_inds(self, permutation: &linnet::permutation::Permutation) -> Self::Permuted {
        ParamTensor {
            tensor: self.tensor.permute_inds(permutation),
            param_type: self.param_type,
        }
    }

    fn permute_reps(self, rep_perm: &linnet::permutation::Permutation) -> Self::Permuted {
        ParamTensor {
            tensor: self.tensor.permute_reps(rep_perm),
            param_type: self.param_type,
        }
    }
}

impl<S> TensorStructure for ParamTensor<S>
where
    S: TensorStructure,
{
    // type R = <T::Structure as TensorStructure>::R;
    type Indexed = ParamTensor<S::Indexed>;
    type Slot = S::Slot;

    fn reindex(
        self,
        indices: &[<Self::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::Indexed>, StructureError> {
        let res = self.tensor.reindex(indices)?;
        Ok(PermutedStructure {
            structure: ParamTensor {
                tensor: res.structure,
                param_type: self.param_type,
            },
            rep_permutation: res.rep_permutation,
            index_permutation: res.index_permutation,
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

impl<S: TensorStructure + Clone> StorageTensor for ParamTensor<S> {
    type Data = Atom;
    type ContainerData<Data> = DataTensor<Data, S>;

    fn map_data<U>(self, f: impl Fn(Self::Data) -> U) -> Self::ContainerData<U> {
        self.tensor.map_data(f)
    }

    fn map_data_ref_mut_result<U, E>(
        &mut self,
        f: impl FnMut(&mut Self::Data) -> Result<U, E>,
    ) -> Result<Self::ContainerData<U>, E> {
        self.tensor.map_data_ref_mut_result(f)
    }

    fn map_data_ref_result_self<E>(
        &self,
        f: impl Fn(&Self::Data) -> Result<Self::Data, E>,
    ) -> Result<Self, E> {
        Ok(ParamTensor {
            tensor: self.tensor.map_data_ref_result_self(f)?,
            param_type: self.param_type,
        })
    }

    fn map_data_mut(&mut self, f: impl FnMut(&mut Self::Data)) {
        self.tensor.map_data_mut(f)
    }

    fn map_data_ref<U>(&self, f: impl Fn(&Self::Data) -> U) -> Self::ContainerData<U> {
        self.tensor.map_data_ref(f)
    }

    fn map_data_ref_mut<U>(
        &mut self,
        f: impl FnMut(&mut Self::Data) -> U,
    ) -> Self::ContainerData<U> {
        self.tensor.map_data_ref_mut(f)
    }

    fn map_data_ref_self(&self, f: impl Fn(&Self::Data) -> Self::Data) -> Self {
        ParamTensor {
            param_type: self.param_type,
            tensor: self.tensor.map_data_ref(f),
        }
    }

    fn map_data_ref_result<U, E>(
        &self,
        f: impl Fn(&Self::Data) -> Result<U, E>,
    ) -> Result<Self::ContainerData<U>, E> {
        self.tensor.map_data_ref_result(f)
    }

    fn map_data_ref_mut_self(&mut self, f: impl FnMut(&mut Self::Data) -> Self::Data) -> Self {
        ParamTensor {
            param_type: self.param_type,
            tensor: self.tensor.map_data_ref_mut_self(f),
        }
    }

    fn map_data_self(self, f: impl Fn(Self::Data) -> Self::Data) -> Self {
        ParamTensor {
            param_type: self.param_type,
            tensor: self.tensor.map_data(f),
        }
    }
}

impl<S> From<DataTensor<Atom, S>> for ParamTensor<S>
where
    S: TensorStructure + Clone,
{
    fn from(tensor: DataTensor<Atom, S>) -> Self {
        ParamTensor {
            tensor,
            param_type: ParamOrComposite::Composite,
        }
    }
}

impl<S> From<SparseTensor<Atom, S>> for ParamTensor<S>
where
    S: TensorStructure + Clone,
{
    fn from(tensor: SparseTensor<Atom, S>) -> Self {
        ParamTensor {
            tensor: tensor.into(),
            param_type: ParamOrComposite::Composite,
        }
    }
}

impl<S> From<DenseTensor<Atom, S>> for ParamTensor<S>
where
    S: TensorStructure + Clone,
{
    fn from(tensor: DenseTensor<Atom, S>) -> Self {
        ParamTensor {
            tensor: tensor.into(),
            param_type: ParamOrComposite::Composite,
        }
    }
}

// impl<C: HasStructure<Structure = S> + Clone, S: TensorStructure + Clone> ScalarMul<Atom>
//     for ParamOrConcrete<C, S> where C
// {
//     type Output = ParamTensor<S>;
//     fn scalar_mul(&self, rhs: &Atom) -> Option<Self::Output> {
//         match self{

//         }
//     }
// }

// impl<Structure: TensorStructure + Serialize + Clone> Serialize for ParamTensor<Structure> {
//     fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
//     where
//         S: serde::Serializer,
//     {
//         let mut state = serializer.serialize_struct("ParamTensor", 3)?;

//         state.serialize_field("param_type", &self.param_type)?;

//         let serialized_tensor = self.tensor.map_data_ref(|a| {
//             let mut v = Vec::new();
//             a.as_view().write(&mut v).unwrap();
//             v
//         });
//         state.serialize_field("tensor", &serialized_tensor)?;

//         let mut symbolica_state = Vec::new();

//         State::export(&mut symbolica_state).unwrap();

//         state.serialize_field("state", &symbolica_state)?;
//         state.end()
//     }
// }

impl<'a, Structure: TensorStructure + Deserialize<'a> + Clone> Deserialize<'a>
    for ParamTensor<Structure>
{
    fn deserialize<D>(deserializer: D) -> std::result::Result<Self, D::Error>
    where
        D: serde::Deserializer<'a>,
    {
        #[derive(Deserialize)]
        struct ParamTensorHelper<Structure: TensorStructure> {
            param_type: ParamOrComposite,
            tensor: DataTensor<Vec<u8>, Structure>,
            state: Vec<u8>,
        }

        let helper = ParamTensorHelper::deserialize(deserializer)?;

        let state = helper.state;

        let mut export = vec![];
        State::export(&mut export).unwrap();

        let map = State::import(&mut Cursor::new(&state), None).unwrap();

        Ok(ParamTensor {
            tensor: helper
                .tensor
                .map_data_ref_result(|a| Atom::import_with_map(&mut a.as_slice(), &map))
                .map_err(serde::de::Error::custom)?,
            param_type: helper.param_type,
        })
    }
}

impl<S: TensorStructure + Clone> HasTensorData for ParamTensor<S> {
    type Data = Atom;
    fn data(&self) -> Vec<Self::Data> {
        self.tensor.data()
    }

    fn hashmap(&self) -> indexmap::IndexMap<ExpandedIndex, Self::Data> {
        self.tensor.hashmap()
    }

    fn indices(&self) -> Vec<ExpandedIndex> {
        self.tensor.indices()
    }

    fn symhashmap(
        &self,
        name: Symbol,
        args: &[Atom],
    ) -> std::collections::HashMap<Atom, Self::Data> {
        self.tensor.symhashmap(name, args)
    }
}
pub enum TensorSet<T: HasStructure> {
    Tensors(Vec<T>),
    Scalars(Vec<T::Scalar>),
}

impl<T: HasStructure> TensorSet<T> {
    pub fn len(&self) -> usize {
        match self {
            TensorSet::Tensors(t) => t.len(),
            TensorSet::Scalars(t) => t.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        match self {
            TensorSet::Tensors(t) => t.is_empty(),
            TensorSet::Scalars(t) => t.is_empty(),
        }
    }

    pub fn push(&mut self, tensor: T) {
        match self {
            TensorSet::Scalars(t) => t.push(tensor.scalar().unwrap()),
            TensorSet::Tensors(t) => t.push(tensor),
        }
    }
}

impl<T: HasStructure + Clone> FromIterator<T> for TensorSet<T> {
    fn from_iter<I: IntoIterator<Item = T>>(iter: I) -> Self {
        let (scalars, set): (Vec<Option<T::Scalar>>, Vec<T>) =
            iter.into_iter().map(|t| (t.clone().scalar(), t)).collect();
        let scalars: Option<Vec<_>> = scalars.into_iter().collect();
        if let Some(s) = scalars {
            TensorSet::Scalars(s)
        } else {
            TensorSet::Tensors(set)
        }
    }
}

pub struct ParamTensorSet<S: TensorStructure + Clone> {
    pub tensors: TensorSet<ParamTensor<S>>,
    size: usize,
}

impl<S: TensorStructure + Clone> ParamTensorSet<S> {
    pub fn new(tensors: Vec<ParamTensor<S>>) -> Self {
        let size = tensors
            .iter()
            .map(|t| t.tensor.actual_size())
            .reduce(|acc, a| acc + a)
            .unwrap();

        ParamTensorSet {
            tensors: tensors.into_iter().collect(),
            size,
        }
    }

    pub fn empty() -> Self {
        ParamTensorSet {
            tensors: [].into_iter().collect(),
            size: 0,
        }
    }

    // pub fn push(&mut self, tensor: ParamTensor<S>) {
    //     self.size += tensor.tensor.actual_size();
    //     self.tensors.push(tensor);
    // }
}

impl<S: TensorStructure> ParamTensor<S> {
    pub fn param(tensor: DataTensor<Atom, S>) -> Self {
        ParamTensor {
            tensor,
            param_type: ParamOrComposite::Param,
        }
    }

    pub fn composite(tensor: DataTensor<Atom, S>) -> Self {
        ParamTensor {
            tensor,
            param_type: ParamOrComposite::Composite,
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, Copy, PartialEq, Eq, Hash, Encode, Decode)]
pub enum ParamOrComposite {
    Param,
    Composite,
}

impl<S: TensorStructure, O: From<S> + TensorStructure> CastStructure<ParamTensor<O>>
    for ParamTensor<S>
{
    fn cast_structure(self) -> ParamTensor<O> {
        ParamTensor {
            tensor: self.tensor.cast_structure(),
            param_type: self.param_type,
        }
    }
}

impl<S: TensorStructure> Shadowable for ParamTensor<S>
where
    S: HasName + Clone,
    S::Name: IntoSymbol,
    S::Args: IntoArgs,
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
}

impl<S: TensorStructure, Const> ShadowMapping<Const> for ParamTensor<S>
where
    S: HasName + Clone,
    S::Name: IntoSymbol,
    S::Args: IntoArgs,
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    fn append_map<T>(
        &self,
        fn_map: &mut FunctionMap<Const>,
        index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> T,
    ) where
        T: TensorCoefficient,
    {
        match self.param_type {
            ParamOrComposite::Param => {}
            ParamOrComposite::Composite => match &self.tensor {
                DataTensor::Dense(d) => {
                    for (i, a) in d.flat_iter() {
                        let labeled_coef = index_to_atom(self.structure(), i);

                        labeled_coef.add_tagged_function(fn_map, a.clone()).unwrap();
                    }
                }
                DataTensor::Sparse(d) => {
                    for (i, a) in d.flat_iter() {
                        let labeled_coef = index_to_atom(self.structure(), i);

                        labeled_coef.add_tagged_function(fn_map, a.clone()).unwrap();
                    }
                }
            },
        }
    }
}

pub mod atomcore;

use symbolica::id::BorrowReplacement;

impl<S: TensorStructure + Clone> ParamTensorSet<S> {
    pub fn eval_tree(
        &self,
        fn_map: &FunctionMap,
        params: &[Atom],
    ) -> Result<EvalTreeTensorSet<SymComplex<Rational>, S>> {
        match &self.tensors {
            TensorSet::Scalars(s) => {
                trace!("turning {} scalars into eval tree", s.len());
                let exprs = s.iter().map(|a| a.as_view()).collect::<Vec<_>>();
                Ok(EvalTreeTensorSet {
                    tensors: TensorsOrScalars::Scalars,
                    eval: (
                        AtomView::to_eval_tree_multiple(&exprs, fn_map, params)
                            .map_err(|s| eyre!(s))?,
                        None,
                    ),
                    size: self.size,
                })
            }
            TensorSet::Tensors(in_tensors) => {
                trace!("turning {} tensors into eval tree", in_tensors.len());
                let mut tensors = vec![];

                let mut atoms = vec![];
                let mut id = 0;
                for t in in_tensors.iter() {
                    let structure = t.structure().clone();
                    let usize_tensor = match &t.tensor {
                        DataTensor::Dense(d) => {
                            let oldid = id;
                            id += d.size().unwrap();
                            for (_, a) in d.flat_iter() {
                                atoms.push(a.as_view());
                            }
                            DataTensor::Dense(DenseTensor::from_data(
                                Vec::from_iter(oldid..id),
                                structure,
                            )?)
                        }
                        DataTensor::Sparse(s) => {
                            let mut t = SparseTensor::empty(structure, 0);
                            for (i, a) in s.flat_iter() {
                                t.set_flat(i, id)?;
                                atoms.push(a.as_view());
                                id += 1;
                            }
                            DataTensor::Sparse(t)
                        }
                    };
                    tensors.push(usize_tensor);
                }

                Ok(EvalTreeTensorSet {
                    tensors: TensorsOrScalars::Tensors(tensors),
                    eval: (
                        AtomView::to_eval_tree_multiple(&atoms, fn_map, params)
                            .map_err(|s| eyre!(s))?,
                        None,
                    ),
                    size: self.size,
                })
            }
        }
    }
}

impl<S: TensorStructure> IteratableTensor for ParamTensor<S> {
    type Data<'a>
        = AtomView<'a>
    where
        Self: 'a;

    fn iter_expanded(&self) -> impl Iterator<Item = (ExpandedIndex, Self::Data<'_>)> {
        self.tensor.iter_expanded().map(|(i, x)| (i, x.as_view()))
    }

    fn iter_flat(&self) -> impl Iterator<Item = (FlatIndex, Self::Data<'_>)> {
        self.tensor.iter_flat().map(|(i, x)| (i, x.as_view()))
    }
}

impl<S> Display for ParamTensor<S>
where
    S: TensorStructure,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.tensor)
    }
}

impl<S: TensorStructure> HasName for ParamTensor<S>
where
    S: HasName,
{
    type Args = S::Args;
    type Name = S::Name;

    fn args(&self) -> Option<Self::Args> {
        self.tensor.args()
    }

    fn name(&self) -> Option<Self::Name> {
        self.tensor.name()
    }

    fn set_name(&mut self, name: Self::Name) {
        if let ParamOrComposite::Composite = self.param_type {
            self.tensor.set_name(name);
        } // never set the name of a param tensor, it is always set by construction
    }
}

#[derive(
    Clone, Debug, PartialEq, Eq, Hash, bincode_trait_derive::Encode, bincode_trait_derive::Decode,
)]
#[trait_decode(trait = symbolica::state::HasStateMap)]
pub enum ParamOrConcrete<C, S> {
    Concrete(C),

    Param(ParamTensor<S>),
}

impl<C, S> crate::network::Ref for ParamOrConcrete<C, S> {
    type Ref<'a>
        = &'a ParamOrConcrete<C, S>
    where
        Self: 'a;

    fn refer(&self) -> Self::Ref<'_> {
        self
    }
}

impl<C: HasStructure<Structure = S> + Clone, S: TensorStructure + Clone> From<ParamTensor<S>>
    for ParamOrConcrete<C, S>
{
    fn from(tensor: ParamTensor<S>) -> Self {
        ParamOrConcrete::Param(tensor)
    }
}

impl<C: HasStructure<Structure = S> + Clone, S: TensorStructure + Clone> ParamOrConcrete<C, S> {
    pub fn replace<'b, P: Into<BorrowedOrOwned<'b, Pattern>>>(
        &self,
        pattern: P,
    ) -> ReplaceBuilderGeneric<'b, &'_ Self, Self> {
        ReplaceBuilderGeneric::new(self, pattern)
    }

    pub fn replace_multiple<T: BorrowReplacement>(&self, replacements: &[T]) -> Self {
        match self {
            ParamOrConcrete::Param(p) => ParamOrConcrete::Param(p.replace_multiple(replacements)),
            _ => self.clone(),
        }
    }
}

impl<C: Display + HasStructure<Structure = S> + Clone, S: TensorStructure + Clone> Display
    for ParamOrConcrete<C, S>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParamOrConcrete::Concrete(c) => c.fmt(f),
            ParamOrConcrete::Param(p) => p.fmt(f),
        }
    }
}

impl<
    U: HasStructure<Structure = O> + Clone,
    C: CastStructure<U> + HasStructure<Structure = S> + Clone,
    S: TensorStructure + Clone,
    O: From<S> + TensorStructure + Clone,
> CastStructure<ParamOrConcrete<U, O>> for ParamOrConcrete<C, S>
{
    fn cast_structure(self) -> ParamOrConcrete<U, O> {
        match self {
            ParamOrConcrete::Concrete(c) => ParamOrConcrete::Concrete(c.cast_structure()),
            ParamOrConcrete::Param(p) => ParamOrConcrete::Param(p.cast_structure()),
        }
    }
}

impl<
    C: HasStructure<Structure = S> + Clone + Shadowable,
    S: TensorStructure + Clone + HasName<Args: IntoArgs, Name: IntoSymbol>,
> Shadowable for ParamOrConcrete<C, S>
where
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
}

impl<
    U,
    C: HasStructure<Structure = S> + Clone + ShadowMapping<U>,
    S: TensorStructure + Clone + HasName<Args: IntoArgs, Name: IntoSymbol>,
> ShadowMapping<U> for ParamOrConcrete<C, S>
where
    <<Self::Structure as TensorStructure>::Slot as IsAbstractSlot>::Aind: ParseableAind,
{
    // fn shadow_with_map<'a, T>(
    //     &'a self,
    //     fn_map: &mut FunctionMap<'a, U>,
    //     index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> T,
    // ) -> Option<ParamTensor<Self::Structure>>
    // where
    //     T: TensorCoefficient,
    // {
    //     match self {
    //         ParamOrConcrete::Concrete(c) => c.shadow_with_map(fn_map, index_to_atom),
    //         ParamOrConcrete::Param(p) => p.shadow_with_map(fn_map, index_to_atom),
    //     }
    // }

    fn append_map<T>(
        &self,
        fn_map: &mut FunctionMap<U>,
        index_to_atom: impl Fn(&Self::Structure, FlatIndex) -> T,
    ) where
        T: TensorCoefficient,
    {
        match self {
            ParamOrConcrete::Concrete(c) => c.append_map(fn_map, index_to_atom),
            ParamOrConcrete::Param(p) => p.append_map(fn_map, index_to_atom),
        }
    }
}

pub enum AtomViewOrConcrete<'a, T> {
    Atom(AtomView<'a>),
    Concrete(T),
}

pub enum AtomOrConcrete<T> {
    Atom(Atom),
    Concrete(T),
}

pub trait Concrete {}

impl<T> From<Atom> for AtomOrConcrete<T> {
    fn from(value: Atom) -> Self {
        AtomOrConcrete::Atom(value)
    }
}

impl<T: Concrete> From<T> for AtomOrConcrete<T> {
    fn from(value: T) -> Self {
        AtomOrConcrete::Concrete(value)
    }
}

impl<C: HasStructure<Structure = S> + Clone, S: TensorStructure + Clone> ParamOrConcrete<C, S> {
    pub fn is_parametric(&self) -> bool {
        matches!(self, ParamOrConcrete::Param(_))
    }

    pub fn try_into_parametric(self) -> Result<ParamTensor<S>, Self> {
        match self {
            ParamOrConcrete::Param(x) => Ok(x),
            _ => Err(self),
        }
    }

    pub fn try_into_concrete(self) -> Result<C, Self> {
        match self {
            ParamOrConcrete::Concrete(x) => Ok(x),
            _ => Err(self),
        }
    }
}

#[derive(Clone, Debug, EnumTryAsInner)]
#[derive_err(Debug)]
pub enum ConcreteOrParam<C> {
    Concrete(C),
    Param(Atom),
}

impl<C> RefOne for ConcreteOrParam<C> {
    fn ref_one(&self) -> Self {
        ConcreteOrParam::Param(Atom::num(1))
    }
}

impl<T: Clone> From<&ConcreteOrParam<T>> for ConcreteOrParam<T> {
    fn from(value: &ConcreteOrParam<T>) -> Self {
        value.clone()
    }
}

impl<C> From<&Atom> for ConcreteOrParam<C> {
    fn from(value: &Atom) -> Self {
        ConcreteOrParam::Param(value.clone())
    }
}

impl<C> crate::algebra::algebraic_traits::One for ConcreteOrParam<C> {
    fn one() -> Self {
        ConcreteOrParam::Param(Atom::num(1))
    }
}

impl<C> crate::algebra::algebraic_traits::Zero for ConcreteOrParam<C> {
    fn zero() -> Self {
        ConcreteOrParam::Param(Atom::Zero)
    }
}

pub mod to_param;

impl<C> ConcreteOrParam<C> {
    pub fn to_param(&mut self)
    where
        C: ToAtom,
    {
        if self.is_concrete() {
            let old = std::mem::replace(self, ConcreteOrParam::Param(Atom::Zero));

            if let ConcreteOrParam::Concrete(r) = old {
                *self = ConcreteOrParam::Param(r.to_atom());
            }
        }
    }
}

impl<C> From<AtomView<'_>> for ConcreteOrParam<C> {
    fn from(value: AtomView<'_>) -> Self {
        ConcreteOrParam::Param(value.into())
    }
}

#[derive(Clone, Debug, EnumTryAsInner)]
#[derive_err(Debug)]
pub enum ConcreteOrParamRef<'a, C> {
    Concrete(C),
    Param(AtomView<'a>),
}

impl crate::network::Ref for Atom {
    type Ref<'a>
        = AtomView<'a>
    where
        Self: 'a;

    fn refer(&self) -> Self::Ref<'_> {
        self.as_view()
    }
}

impl<C: crate::network::Ref> crate::network::Ref for ConcreteOrParam<C> {
    type Ref<'a>
        = ConcreteOrParamRef<'a, C::Ref<'a>>
    where
        Self: 'a;

    fn refer(&self) -> Self::Ref<'_> {
        match self {
            ConcreteOrParam::Concrete(c) => ConcreteOrParamRef::Concrete(c.refer()),
            ConcreteOrParam::Param(p) => ConcreteOrParamRef::Param(p.refer()),
        }
    }
}

impl<C: Default> Default for ConcreteOrParam<C> {
    fn default() -> Self {
        ConcreteOrParam::Concrete(C::default())
    }
}

impl<C: std::ops::Neg<Output = C>> std::ops::Neg for ConcreteOrParam<C> {
    type Output = ConcreteOrParam<C>;

    fn neg(self) -> Self {
        match self {
            ConcreteOrParam::Concrete(c) => ConcreteOrParam::Concrete(-c),
            ConcreteOrParam::Param(p) => ConcreteOrParam::Param(-p),
        }
    }
}

impl<C: std::ops::Neg<Output = C>> std::ops::Neg for ConcreteOrParamRef<'_, C> {
    type Output = ConcreteOrParam<C>;

    fn neg(self) -> Self::Output {
        match self {
            ConcreteOrParamRef::Concrete(c) => ConcreteOrParam::Concrete(-c),
            ConcreteOrParamRef::Param(p) => ConcreteOrParam::Param(-p),
        }
    }
}

#[derive(Clone, Debug, EnumTryAsInner)]
#[derive_err(Debug)]
pub enum ConcreteOrParamView<'a, C> {
    Concrete(C),
    Param(AtomView<'a>),
}

#[derive(Debug, EnumTryAsInner)]
#[derive_err(Debug)]
pub enum ConcreteOrParamViewMut<'a, C> {
    Concrete(C),
    Param(&'a mut Atom),
}

impl<D: Display> Display for ConcreteOrParam<D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ConcreteOrParam::Concrete(c) => c.fmt(f),
            ConcreteOrParam::Param(p) => write!(f, "{}", p),
        }
    }
}

impl<S: TensorStructure + Clone> SparseOrDense for ParamTensor<S> {
    fn to_dense(self) -> Self {
        ParamTensor {
            tensor: self.tensor.to_dense(),
            param_type: self.param_type,
        }
    }

    fn to_dense_mut(&mut self) {
        self.tensor.to_dense_mut();
    }
    fn to_sparse_mut(&mut self) {
        self.tensor.to_sparse_mut();
    }

    fn to_sparse(self) -> Self {
        ParamTensor {
            tensor: self.tensor.to_sparse(),
            param_type: self.param_type,
        }
    }
}

impl<C, S> SparseOrDense for ParamOrConcrete<C, S>
where
    C: SparseOrDense + Clone + HasStructure<Structure = S>,
    S: TensorStructure + Clone,
{
    fn to_dense(self) -> Self {
        match self {
            ParamOrConcrete::Concrete(x) => ParamOrConcrete::Concrete(x.to_dense()),
            ParamOrConcrete::Param(x) => ParamOrConcrete::Param(x.to_dense()),
        }
    }

    fn to_dense_mut(&mut self) {
        match self {
            ParamOrConcrete::Concrete(x) => x.to_dense_mut(),
            ParamOrConcrete::Param(x) => x.to_dense_mut(),
        }
    }
    fn to_sparse_mut(&mut self) {
        match self {
            ParamOrConcrete::Concrete(x) => x.to_dense_mut(),
            ParamOrConcrete::Param(x) => x.to_dense_mut(),
        }
    }

    fn to_sparse(self) -> Self {
        match self {
            ParamOrConcrete::Concrete(x) => ParamOrConcrete::Concrete(x.to_sparse()),
            ParamOrConcrete::Param(x) => ParamOrConcrete::Param(x.to_sparse()),
        }
    }
}

impl<
    Aind: AbsInd,
    C: PermuteTensor<Id = C, Permuted = C>,
    S: Clone + Into<IndexLess<R, Aind>>,
    R: RepName<Dual = R>,
> PermuteTensor for ParamOrConcrete<C, S>
where
    S: TensorStructure<Slot = Slot<R, Aind>> + PermuteTensor<IdSlot = Slot<R, Aind>, Id = S>,
    C: HasStructure<Structure = S> + TensorStructure,
{
    type Id = ParamOrConcrete<C, S>;
    type IdSlot = C::IdSlot;
    type Permuted = ParamOrConcrete<C, S>;

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
        ParamOrConcrete::Concrete(C::id(i, j))
    }

    fn permute_inds(self, permutation: &linnet::permutation::Permutation) -> Self::Permuted {
        match self {
            ParamOrConcrete::Param(d) => ParamOrConcrete::Param(d.permute_inds(permutation)),
            ParamOrConcrete::Concrete(s) => ParamOrConcrete::Concrete(s.permute_inds(permutation)),
        }
    }

    fn permute_reps(self, rep_perm: &linnet::permutation::Permutation) -> Self::Permuted {
        match self {
            ParamOrConcrete::Param(d) => ParamOrConcrete::Param(d.permute_reps(rep_perm)),
            ParamOrConcrete::Concrete(s) => ParamOrConcrete::Concrete(s.permute_reps(rep_perm)),
        }
    }
}

impl<C, S> TensorStructure for ParamOrConcrete<C, S>
where
    C: HasStructure<Structure = S> + TensorStructure,
    C::Indexed: HasStructure<Structure = S::Indexed>
        + TensorStructure<Slot = <S::Indexed as TensorStructure>::Slot>,
    S: TensorStructure<Slot = C::Slot>,
{
    // type R = <T::Structure as TensorStructure>::R;
    type Indexed = ParamOrConcrete<C::Indexed, S::Indexed>;
    type Slot = S::Slot;

    fn reindex(
        self,
        indices: &[<Self::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::Indexed>, StructureError> {
        Ok(match self {
            ParamOrConcrete::Concrete(c) => {
                let res = c.reindex(indices)?;
                PermutedStructure {
                    rep_permutation: res.rep_permutation,
                    index_permutation: res.index_permutation,
                    structure: ParamOrConcrete::Concrete(res.structure),
                }
            }
            ParamOrConcrete::Param(p) => {
                let res = p.reindex(indices)?;
                PermutedStructure {
                    rep_permutation: res.rep_permutation,
                    index_permutation: res.index_permutation,
                    structure: ParamOrConcrete::Param(res.structure),
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

impl<C, S> HasStructure for ParamOrConcrete<C, S>
where
    C: HasStructure<Structure = S>,
    S: TensorStructure,
{
    type Scalar = ConcreteOrParam<C::Scalar>;
    type ScalarRef<'a>
        = ConcreteOrParamRef<'a, C::ScalarRef<'a>>
    where
        Self: 'a;
    type Structure = S;
    type Store<U>
        = ParamOrConcrete<C::Store<U>, U>
    where
        U: TensorStructure;

    fn map_structure<O: TensorStructure>(self, f: impl Fn(Self::Structure) -> O) -> Self::Store<O> {
        match self {
            ParamOrConcrete::Concrete(c) => ParamOrConcrete::Concrete(c.map_structure(f)),
            ParamOrConcrete::Param(p) => ParamOrConcrete::Param(p.map_structure(f)),
        }
    }

    fn map_structure_result<O: TensorStructure, Er>(
        self,
        _f: impl Fn(Self::Structure) -> Result<O, Er>,
    ) -> std::result::Result<Self::Store<O>, Er> {
        match self {
            ParamOrConcrete::Concrete(c) => {
                Ok(ParamOrConcrete::Concrete(c.map_structure_result(_f)?))
            }
            ParamOrConcrete::Param(p) => Ok(ParamOrConcrete::Param(p.map_structure_result(_f)?)),
        }
    }

    fn structure(&self) -> &Self::Structure {
        match self {
            ParamOrConcrete::Concrete(x) => x.structure(),
            ParamOrConcrete::Param(x) => x.structure(),
        }
    }

    fn map_same_structure(self, f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self {
        match self {
            ParamOrConcrete::Concrete(x) => ParamOrConcrete::Concrete(x.map_same_structure(f)),
            ParamOrConcrete::Param(x) => ParamOrConcrete::Param(x.map_same_structure(f)),
        }
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        match self {
            ParamOrConcrete::Concrete(x) => x.mut_structure(),
            ParamOrConcrete::Param(x) => x.mut_structure(),
        }
    }

    fn scalar(self) -> Option<Self::Scalar> {
        match self {
            ParamOrConcrete::Concrete(x) => x.scalar().map(ConcreteOrParam::Concrete),
            ParamOrConcrete::Param(x) => x.scalar().map(ConcreteOrParam::Param),
        }
    }

    fn scalar_ref(&self) -> Option<Self::ScalarRef<'_>> {
        match self {
            ParamOrConcrete::Concrete(x) => x.scalar_ref().map(ConcreteOrParamRef::Concrete),
            ParamOrConcrete::Param(x) => x.scalar_ref().map(ConcreteOrParamRef::Param),
        }
    }
}

impl<C, S> SetTensorData for ParamOrConcrete<C, S>
where
    C: HasStructure<Structure = S> + Clone + SetTensorData,
    S: TensorStructure + Clone,
{
    type SetData = ConcreteOrParam<C::SetData>;

    fn set(
        &mut self,
        indices: &[crate::structure::concrete_index::ConcreteIndex],
        value: Self::SetData,
    ) -> Result<()> {
        match self {
            ParamOrConcrete::Concrete(x) => x.set(
                indices,
                value
                    .try_into_concrete()
                    .map_err(|r| eyre!(r.to_string()))?,
            ),
            ParamOrConcrete::Param(x) => x.set(
                indices,
                value.try_into_param().map_err(|r| eyre!(r.to_string()))?,
            ),
        }
    }

    fn set_flat(&mut self, index: FlatIndex, value: Self::SetData) -> Result<()> {
        match self {
            ParamOrConcrete::Concrete(x) => x.set_flat(
                index,
                value
                    .try_into_concrete()
                    .map_err(|r| eyre!(r.to_string()))?,
            ),
            ParamOrConcrete::Param(x) => x.set_flat(
                index,
                value.try_into_param().map_err(|r| eyre!(r.to_string()))?,
            ),
        }
    }
}

impl<C, S> GetTensorData for ParamOrConcrete<C, S>
where
    C: HasStructure<Structure = S> + Clone + GetTensorData<GetDataOwned: Clone>,
    S: TensorStructure + Clone,
{
    type GetDataRef<'a>
        = ConcreteOrParamView<'a, C::GetDataRef<'a>>
    where
        Self: 'a;
    type GetDataRefMut<'a>
        = ConcreteOrParamViewMut<'a, C::GetDataRefMut<'a>>
    where
        Self: 'a;
    type GetDataOwned = ConcreteOrParam<C::GetDataOwned>;
    fn get_ref<D: AsRef<[ConcreteIndex]>>(&self, indices: D) -> Result<Self::GetDataRef<'_>> {
        match self {
            ParamOrConcrete::Concrete(x) => x.get_ref(indices).map(ConcreteOrParamView::Concrete),
            ParamOrConcrete::Param(x) => x.get_ref(indices).map(ConcreteOrParamView::Param),
        }
    }

    fn get_ref_linear(&self, index: FlatIndex) -> Option<Self::GetDataRef<'_>> {
        match self {
            ParamOrConcrete::Concrete(x) => {
                x.get_ref_linear(index).map(ConcreteOrParamView::Concrete)
            }
            ParamOrConcrete::Param(x) => x.get_ref_linear(index).map(ConcreteOrParamView::Param),
        }
    }

    fn get_mut_linear(&mut self, index: FlatIndex) -> Option<Self::GetDataRefMut<'_>> {
        match self {
            ParamOrConcrete::Concrete(x) => x
                .get_mut_linear(index)
                .map(ConcreteOrParamViewMut::Concrete),
            ParamOrConcrete::Param(x) => x.get_mut_linear(index).map(ConcreteOrParamViewMut::Param),
        }
    }

    fn get_owned<D: AsRef<[ConcreteIndex]>>(&self, indices: D) -> Result<Self::GetDataOwned>
    where
        Self::GetDataOwned: Clone,
    {
        match self {
            ParamOrConcrete::Concrete(x) => x.get_owned(indices).map(ConcreteOrParam::Concrete),
            ParamOrConcrete::Param(x) => x.get_owned(indices).map(ConcreteOrParam::Param),
        }
    }

    fn get_owned_linear(&self, index: FlatIndex) -> Option<Self::GetDataOwned>
    where
        Self::GetDataOwned: Clone,
    {
        match self {
            ParamOrConcrete::Concrete(x) => {
                x.get_owned_linear(index).map(ConcreteOrParam::Concrete)
            }
            ParamOrConcrete::Param(x) => x.get_owned_linear(index).map(ConcreteOrParam::Param),
        }
    }
}

impl<T> From<Atom> for ConcreteOrParam<T> {
    fn from(value: Atom) -> Self {
        ConcreteOrParam::Param(value)
    }
}

impl<T: Into<Atom>> From<ConcreteOrParam<T>> for Atom {
    fn from(value: ConcreteOrParam<T>) -> Self {
        match value {
            ConcreteOrParam::Concrete(x) => x.into(),
            ConcreteOrParam::Param(x) => x,
        }
    }
}

impl<T: Into<Coefficient>> From<RealOrComplex<T>> for Atom {
    fn from(value: RealOrComplex<T>) -> Self {
        match value {
            RealOrComplex::Real(x) => Atom::num(x),
            RealOrComplex::Complex(x) => {
                let (re, im) = (Atom::num(x.re), Atom::num(x.im));
                let i = Atom::i();
                re + im * i
            }
        }
    }
}

impl<C, S> ScalarTensor for ParamOrConcrete<C, S>
where
    C: HasStructure<Structure = S> + Clone + ScalarTensor + TensorStructure,
    C::Indexed: HasStructure + TensorStructure,
    S: TensorStructure + ScalarStructure + Clone,
{
    fn new_scalar(scalar: Self::Scalar) -> Self {
        match scalar {
            ConcreteOrParam::Concrete(x) => ParamOrConcrete::Concrete(C::new_scalar(x)),
            ConcreteOrParam::Param(x) => ParamOrConcrete::Param(ParamTensor::new_scalar(x)),
        }
    }
}

impl<C, S> TracksCount for ParamOrConcrete<C, S>
where
    C: TracksCount + HasStructure<Structure = S> + Clone,
    S: TensorStructure + TracksCount + Clone,
{
    fn contractions_num(&self) -> usize {
        match self {
            ParamOrConcrete::Concrete(x) => x.contractions_num(),
            ParamOrConcrete::Param(x) => x.contractions_num(),
        }
    }
}

impl<C, S> HasName for ParamOrConcrete<C, S>
where
    C: HasName + HasStructure<Structure = S> + Clone,
    S: TensorStructure + HasName<Name = C::Name, Args = C::Args> + Clone,
{
    type Args = C::Args;
    type Name = C::Name;

    fn args(&self) -> Option<Self::Args> {
        match self {
            ParamOrConcrete::Concrete(x) => x.args(),
            ParamOrConcrete::Param(x) => x.args(),
        }
    }

    fn name(&self) -> Option<Self::Name> {
        match self {
            ParamOrConcrete::Concrete(x) => x.name(),
            ParamOrConcrete::Param(x) => x.name(),
        }
    }

    fn set_name(&mut self, name: Self::Name) {
        match self {
            ParamOrConcrete::Concrete(x) => x.set_name(name),
            ParamOrConcrete::Param(x) => x.set_name(name),
        }
    }
}

impl<C: IteratableTensor + Clone, S: TensorStructure + Clone> IteratableTensor
    for ParamOrConcrete<C, S>
where
    C: HasStructure<Structure = S> + TensorStructure,
    C::Indexed: HasStructure<Structure = S::Indexed>
        + TensorStructure<Slot = <S::Indexed as TensorStructure>::Slot>,
    S: TensorStructure<Slot = C::Slot>,
{
    type Data<'a>
        = AtomViewOrConcrete<'a, C::Data<'a>>
    where
        Self: 'a;

    fn iter_flat(&self) -> impl Iterator<Item = (FlatIndex, Self::Data<'_>)> {
        match self {
            ParamOrConcrete::Concrete(x) => IteratorEnum::A(
                x.iter_flat()
                    .map(|(i, x)| (i, AtomViewOrConcrete::Concrete(x))),
            ),
            ParamOrConcrete::Param(x) => {
                IteratorEnum::B(x.iter_flat().map(|(i, x)| (i, AtomViewOrConcrete::Atom(x))))
            }
        }
    }

    fn iter_expanded(&self) -> impl Iterator<Item = (ExpandedIndex, Self::Data<'_>)> {
        match self {
            ParamOrConcrete::Concrete(x) => IteratorEnum::A(
                x.iter_expanded()
                    .map(|(i, x)| (i, AtomViewOrConcrete::Concrete(x))),
            ),
            ParamOrConcrete::Param(x) => IteratorEnum::B(
                x.iter_expanded()
                    .map(|(i, x)| (i, AtomViewOrConcrete::Atom(x))),
            ),
        }
    }
}

pub type MixedTensor<T = f64, S = NamedStructure<Symbol, Vec<Atom>>> =
    ParamOrConcrete<RealOrComplexTensor<T, S>, S>;

impl<I: TensorStructure + Clone, T: Clone> MixedTensor<T, I> {
    pub fn evaluate_real<A: AtomCore + KeyLookup, F: Fn(&Rational) -> T + Copy>(
        &mut self,
        coeff_map: F,
        const_map: &HashMap<A, T>,
        function_map: &HashMap<Symbol, EvaluationFn<A, T>>,
    ) where
        T: Real + for<'c> From<&'c Rational>,
    {
        let content = match self {
            MixedTensor::Param(x) => Some(x),
            _ => None,
        };

        if let Some(x) = content {
            *self = MixedTensor::Concrete(RealOrComplexTensor::Real(
                x.evaluate(coeff_map, const_map, function_map).unwrap(),
            ));
        }
    }

    pub fn evaluate_complex<A: AtomCore + KeyLookup, F: Fn(&Rational) -> SymComplex<T> + Copy>(
        &mut self,
        coeff_map: F,
        const_map: &HashMap<A, SymComplex<T>>,
        function_map: &HashMap<Symbol, EvaluationFn<A, SymComplex<T>>>,
    ) where
        T: Real + for<'c> From<&'c Rational>,
        SymComplex<T>: Real + for<'c> From<&'c Rational>,
    {
        let content = match self {
            MixedTensor::Param(x) => Some(x),
            _ => None,
        };

        if let Some(x) = content {
            *self = MixedTensor::Concrete(RealOrComplexTensor::Complex(
                x.evaluate(coeff_map, const_map, function_map)
                    .unwrap()
                    .map_data(|c| c.into()),
            ));
        }
    }
}

impl<I> DenseTensor<Atom, I>
where
    I: Clone + TensorStructure,
{
    pub fn append_const_map<'a, 'b, T, U>(
        &'a self,
        data: &DenseTensor<T, I>,
        const_map: &mut HashMap<AtomView<'b>, U>,
    ) where
        I: TensorStructure,
        T: Copy,
        U: From<T>,
        'a: 'b,
    {
        for ((i, a), (j, v)) in self.flat_iter().zip(data.flat_iter()) {
            assert_eq!(i, j);
            const_map.insert(a.as_view(), (*v).into());
        }
    }
}

impl<S> HasStructure for ParamTensor<S>
where
    S: TensorStructure,
{
    type Structure = S;
    type Scalar = Atom;
    type ScalarRef<'a>
        = AtomView<'a>
    where
        Self: 'a;
    type Store<U>
        = ParamTensor<U>
    where
        U: TensorStructure;

    fn map_structure<Ss>(self, f: impl Fn(Self::Structure) -> Ss) -> Self::Store<Ss>
    where
        Ss: TensorStructure,
    {
        ParamTensor {
            param_type: self.param_type,
            tensor: self.tensor.map_structure(f),
        }
    }

    fn map_structure_result<O: TensorStructure, Er>(
        self,
        f: impl Fn(Self::Structure) -> Result<O, Er>,
    ) -> std::result::Result<Self::Store<O>, Er> {
        Ok(ParamTensor {
            param_type: self.param_type,
            tensor: self.tensor.map_structure_result(f)?,
        })
    }

    fn structure(&self) -> &Self::Structure {
        self.tensor.structure()
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        self.tensor.mut_structure()
    }

    fn map_same_structure(self, f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self {
        ParamTensor {
            tensor: self.tensor.map_same_structure(f),
            ..self
        }
    }

    fn scalar(self) -> Option<Self::Scalar> {
        self.tensor.scalar()
    }

    fn scalar_ref(&self) -> Option<Self::ScalarRef<'_>> {
        self.tensor.scalar_ref().map(|a| a.as_view())
    }
}

impl<S: TensorStructure> GetTensorData for ParamTensor<S> {
    type GetDataRef<'a>
        = AtomView<'a>
    where
        Self: 'a;
    type GetDataRefMut<'a>
        = &'a mut Atom
    where
        Self: 'a;
    type GetDataOwned = Atom;
    fn get_ref<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<Self::GetDataRef<'_>> {
        self.tensor.get_ref(indices).map(|x| x.as_view())
    }

    fn get_ref_linear(&self, index: FlatIndex) -> Option<Self::GetDataRef<'_>> {
        self.tensor.get_ref_linear(index).map(|x| x.as_view())
    }

    fn get_mut_linear(&mut self, index: FlatIndex) -> Option<&mut Atom> {
        self.tensor.get_mut_linear(index)
    }

    fn get_owned<C: AsRef<[ConcreteIndex]>>(&self, indices: C) -> Result<Self::GetDataOwned>
    where
        Self::GetDataOwned: Clone,
    {
        self.tensor.get_owned(indices)
    }

    fn get_owned_linear(&self, index: FlatIndex) -> Option<Self::GetDataOwned>
    where
        Self::GetDataOwned: Clone,
    {
        self.tensor.get_owned_linear(index)
    }
}

impl<S: TensorStructure> SetTensorData for ParamTensor<S> {
    type SetData = Atom;

    fn set(
        &mut self,
        indices: &[crate::structure::concrete_index::ConcreteIndex],
        value: Self::SetData,
    ) -> Result<()> {
        self.tensor.set(indices, value)
    }

    fn set_flat(&mut self, index: FlatIndex, value: Self::SetData) -> Result<()> {
        self.tensor.set_flat(index, value)
    }
}

impl<S> ScalarTensor for ParamTensor<S>
where
    S: TensorStructure + ScalarStructure,
{
    fn new_scalar(scalar: Self::Scalar) -> Self {
        ParamTensor {
            tensor: DataTensor::new_scalar(scalar),
            param_type: ParamOrComposite::Composite,
        }
    }
}

impl<S> crate::network::Ref for ParamTensor<S> {
    type Ref<'a>
        = &'a ParamTensor<S>
    where
        Self: 'a;

    fn refer(&self) -> Self::Ref<'_> {
        self
    }
}

impl<S> TracksCount for ParamTensor<S>
where
    S: TensorStructure + TracksCount,
{
    fn contractions_num(&self) -> usize {
        self.tensor.contractions_num()
    }
}

// pub type MixedTensors = MixedTensor<HistoryStructure<Symbol>>;

impl<I, T: Clone> From<DenseTensor<T, I>> for MixedTensor<T, I>
where
    I: TensorStructure + Clone,
{
    fn from(other: DenseTensor<T, I>) -> Self {
        MixedTensor::Concrete(RealOrComplexTensor::Real(DataTensor::Dense(other)))
    }
}

impl<I, T: Clone> From<SparseTensor<T, I>> for MixedTensor<T, I>
where
    I: TensorStructure + Clone,
{
    fn from(other: SparseTensor<T, I>) -> Self {
        MixedTensor::Concrete(RealOrComplexTensor::Real(DataTensor::Sparse(other)))
    }
}

impl<I, T: Clone> From<DenseTensor<Complex<T>, I>> for MixedTensor<T, I>
where
    I: TensorStructure + Clone,
{
    fn from(other: DenseTensor<Complex<T>, I>) -> Self {
        MixedTensor::Concrete(RealOrComplexTensor::Complex(DataTensor::Dense(other)))
    }
}

impl<I, T: Clone> From<SparseTensor<Complex<T>, I>> for MixedTensor<T, I>
where
    I: TensorStructure + Clone,
{
    fn from(other: SparseTensor<Complex<T>, I>) -> Self {
        MixedTensor::Concrete(RealOrComplexTensor::Complex(DataTensor::Sparse(other)))
    }
}

impl<I, T: Clone> MixedTensor<T, I>
where
    I: TensorStructure + Clone,
{
    pub fn param(other: DataTensor<Atom, I>) -> Self {
        MixedTensor::Param(ParamTensor::param(other))
    }

    pub fn composite(other: DataTensor<Atom, I>) -> Self {
        MixedTensor::Param(ParamTensor::composite(other))
    }
}

impl<I> Trace for ParamTensor<I>
where
    I: TensorStructure + Clone + StructureContract,
{
    fn internal_contract(&self) -> Self {
        ParamTensor {
            tensor: self.tensor.internal_contract(),
            param_type: self.param_type,
        }
    }
}

impl<I> Contract<ParamTensor<I>> for ParamTensor<I>
where
    I: TensorStructure + Clone + StructureContract,
{
    type LCM = ParamTensor<I>;
    fn contract(&self, other: &ParamTensor<I>) -> Result<Self::LCM, ContractionError> {
        let s = self.tensor.contract(&other.tensor)?;

        match (self.param_type, other.param_type) {
            (ParamOrComposite::Param, ParamOrComposite::Param) => Ok(ParamTensor::param(s)),
            (ParamOrComposite::Composite, ParamOrComposite::Composite) => {
                Ok(ParamTensor::composite(s))
            }
            (ParamOrComposite::Param, ParamOrComposite::Composite) => Ok(ParamTensor::composite(s)),
            (ParamOrComposite::Composite, ParamOrComposite::Param) => Ok(ParamTensor::composite(s)),
        }
    }
}

impl<I, T> Trace for ParamOrConcrete<DataTensor<T, I>, I>
where
    I: TensorStructure + Clone + StructureContract,
    T: ContractableWith<Atom, Out = Atom>
        + ContractableWith<T, Out = T>
        + Clone
        + FallibleMul<Output = T>
        + FallibleAddAssign<T>
        + FallibleSubAssign<T>
        + RefZero
        + IsZero,
{
    fn internal_contract(&self) -> Self {
        match self {
            ParamOrConcrete::Param(p) => ParamOrConcrete::Param(p.internal_contract()),
            ParamOrConcrete::Concrete(c) => ParamOrConcrete::Concrete(c.internal_contract()),
        }
    }
}

impl<I, T> Contract<ParamOrConcrete<DataTensor<T, I>, I>> for ParamOrConcrete<DataTensor<T, I>, I>
where
    I: TensorStructure + Clone + StructureContract,
    T: ContractableWith<Atom, Out = Atom>
        + ContractableWith<T, Out = T>
        + Clone
        + FallibleMul<Output = T>
        + FallibleAddAssign<T>
        + FallibleSubAssign<T>
        + RefZero
        + IsZero,
    Atom: TrySmallestUpgrade<T, LCM = Atom>, // Atom: ContractableWith<T, Out = Atom> + ContractableWith<Atom, Out = Atom>,
{
    type LCM = ParamOrConcrete<DataTensor<T, I>, I>;
    fn contract(
        &self,
        other: &ParamOrConcrete<DataTensor<T, I>, I>,
    ) -> Result<Self::LCM, ContractionError> {
        match (self, other) {
            (ParamOrConcrete::Param(s), ParamOrConcrete::Param(o)) => {
                Ok(ParamOrConcrete::Param(s.contract(o)?))
            }
            (ParamOrConcrete::Param(s), ParamOrConcrete::Concrete(o)) => match s.param_type {
                ParamOrComposite::Composite => Ok(ParamOrConcrete::Param(ParamTensor::composite(
                    s.tensor.contract(o)?,
                ))),
                ParamOrComposite::Param => Ok(ParamOrConcrete::Param(ParamTensor::composite(
                    s.tensor.contract(o)?,
                ))),
            },
            (ParamOrConcrete::Concrete(s), ParamOrConcrete::Param(o)) => match o.param_type {
                ParamOrComposite::Composite => Ok(ParamOrConcrete::Param(ParamTensor::composite(
                    s.contract(&o.tensor)?,
                ))),
                ParamOrComposite::Param => Ok(ParamOrConcrete::Param(ParamTensor::composite(
                    s.contract(&o.tensor)?,
                ))),
            },
            (ParamOrConcrete::Concrete(s), ParamOrConcrete::Concrete(o)) => {
                Ok(ParamOrConcrete::Concrete(s.contract(o)?))
            }
        }
    }
}
impl<I, T> Trace for ParamOrConcrete<RealOrComplexTensor<T, I>, I>
where
    I: TensorStructure + Clone + StructureContract,
    T: ContractableWith<Atom, Out = Atom>
        + ContractableWith<T, Out = T>
        + ContractableWith<Complex<T>, Out = Complex<T>>
        + Clone
        + FallibleMul<Output = T>
        + FallibleAddAssign<T>
        + FallibleSubAssign<T>
        + RefZero
        + IsZero,
    Complex<T>: ContractableWith<Atom, Out = Atom>
        + ContractableWith<T, Out = Complex<T>>
        + ContractableWith<Complex<T>, Out = Complex<T>>
        + Clone
        + FallibleMul<Output = Complex<T>>
        + FallibleAddAssign<Complex<T>>
        + FallibleSubAssign<Complex<T>>
        + RefZero
        + IsZero,
    Atom: TrySmallestUpgrade<T, LCM = Atom> + TrySmallestUpgrade<Complex<T>, LCM = Atom>,
{
    fn internal_contract(&self) -> Self {
        match self {
            ParamOrConcrete::Param(p) => ParamOrConcrete::Param(p.internal_contract()),
            ParamOrConcrete::Concrete(c) => ParamOrConcrete::Concrete(c.internal_contract()),
        }
    }
}

impl<I, T> Contract<ParamOrConcrete<RealOrComplexTensor<T, I>, I>>
    for ParamOrConcrete<RealOrComplexTensor<T, I>, I>
where
    I: TensorStructure + Clone + StructureContract,
    T: ContractableWith<Atom, Out = Atom>
        + ContractableWith<T, Out = T>
        + ContractableWith<Complex<T>, Out = Complex<T>>
        + Clone
        + FallibleMul<Output = T>
        + FallibleAddAssign<T>
        + FallibleSubAssign<T>
        + RefZero
        + IsZero,
    Complex<T>: ContractableWith<Atom, Out = Atom>
        + ContractableWith<T, Out = Complex<T>>
        + ContractableWith<Complex<T>, Out = Complex<T>>
        + Clone
        + FallibleMul<Output = Complex<T>>
        + FallibleAddAssign<Complex<T>>
        + FallibleSubAssign<Complex<T>>
        + RefZero
        + IsZero,
    Atom: TrySmallestUpgrade<T, LCM = Atom> + TrySmallestUpgrade<Complex<T>, LCM = Atom>,
{
    type LCM = ParamOrConcrete<RealOrComplexTensor<T, I>, I>;
    fn contract(
        &self,
        other: &ParamOrConcrete<RealOrComplexTensor<T, I>, I>,
    ) -> Result<Self::LCM, ContractionError> {
        match (self, other) {
            (ParamOrConcrete::Param(s), ParamOrConcrete::Param(o)) => {
                Ok(ParamOrConcrete::Param(s.contract(o)?))
            }
            (ParamOrConcrete::Param(s), ParamOrConcrete::Concrete(o)) => match (s.param_type, o) {
                (ParamOrComposite::Composite, RealOrComplexTensor::Real(o)) => Ok(
                    ParamOrConcrete::Param(ParamTensor::composite(s.tensor.contract(o)?)),
                ),
                (ParamOrComposite::Composite, RealOrComplexTensor::Complex(o)) => Ok(
                    ParamOrConcrete::Param(ParamTensor::composite(s.tensor.contract(o)?)),
                ),
                (ParamOrComposite::Param, RealOrComplexTensor::Real(o)) => Ok(
                    ParamOrConcrete::Param(ParamTensor::composite(s.tensor.contract(o)?)),
                ),
                (ParamOrComposite::Param, RealOrComplexTensor::Complex(o)) => Ok(
                    ParamOrConcrete::Param(ParamTensor::composite(s.tensor.contract(o)?)),
                ),
            },
            (ParamOrConcrete::Concrete(s), ParamOrConcrete::Param(o)) => match (o.param_type, s) {
                (ParamOrComposite::Composite, RealOrComplexTensor::Real(s)) => Ok(
                    ParamOrConcrete::Param(ParamTensor::composite(o.tensor.contract(s)?)),
                ),
                (ParamOrComposite::Composite, RealOrComplexTensor::Complex(s)) => Ok(
                    ParamOrConcrete::Param(ParamTensor::composite(o.tensor.contract(s)?)),
                ),
                (ParamOrComposite::Param, RealOrComplexTensor::Real(s)) => Ok(
                    ParamOrConcrete::Param(ParamTensor::composite(o.tensor.contract(s)?)),
                ),
                (ParamOrComposite::Param, RealOrComplexTensor::Complex(s)) => Ok(
                    ParamOrConcrete::Param(ParamTensor::composite(o.tensor.contract(s)?)),
                ),
            },
            (ParamOrConcrete::Concrete(s), ParamOrConcrete::Concrete(o)) => {
                Ok(ParamOrConcrete::Concrete(s.contract(o)?))
            }
        }
    }
}

pub type EvalTreeTensor<T, S> = EvalTensor<EvalTree<T>, S>;

pub type EvalTreeTensorSet<T, S> = EvalTensorSet<(EvalTree<T>, Option<Vec<Expression<T>>>), S>;

impl<S: Clone + TensorStructure> EvalTreeTensorSet<SymComplex<Rational>, S> {
    pub fn horner_scheme(&mut self) {
        self.eval.0.horner_scheme()
    }
    pub fn optimize_horner_scheme(&mut self, settings: &OptimizationSettings) {
        let _scheme = self.eval.1.take();
        self.eval.1 = Some(self.eval.0.optimize_horner_scheme(settings));
    }
}

impl<S: TensorStructure> EvalTreeTensorSet<SymComplex<Rational>, S> {
    #[allow(clippy::type_complexity)]
    pub fn linearize(
        mut self,
        settings: &OptimizationSettings,
    ) -> EvalTensorSet<
        (
            ExpressionEvaluator<SymComplex<Rational>>,
            Option<Vec<Expression<SymComplex<Rational>>>>,
        ),
        S,
    > {
        EvalTensorSet {
            eval: (self.eval.0.linearize(settings), self.eval.1),
            tensors: self.tensors,
            size: self.size,
        }
    }
}

impl<T, S: TensorStructure> EvalTreeTensorSet<T, S> {
    pub fn map_coeff<T2, F: Fn(&T) -> T2>(&self, f: &F) -> EvalTreeTensorSet<T2, S>
    where
        T: Clone + PartialEq,
        S: Clone,
    {
        EvalTreeTensorSet {
            eval: (self.eval.0.map_coeff(f), None),
            tensors: self.tensors.clone(),
            size: self.size,
        }
        // self.map_data_ref(|x| x.map_coeff(f))
    }

    pub fn common_subexpression_elimination(&mut self)
    where
        T: Debug + Hash + Eq + Clone + Default + InternalOrdering,
    {
        self.eval.0.common_subexpression_elimination()
    }

    pub fn evaluate(&mut self, params: &[T]) -> TensorSet<DataTensor<T, S>>
    where
        T: Real,
        S: TensorStructure + Clone,
    {
        let zero = params[0].zero();

        let mut elements = vec![zero; self.size];
        self.eval.0.evaluate(params, &mut elements);

        match &self.tensors {
            TensorsOrScalars::Scalars => TensorSet::Scalars(elements),
            TensorsOrScalars::Tensors(t) => {
                let mut out_tensors = Vec::with_capacity(t.len());
                for t in t.iter() {
                    out_tensors.push(t.map_data_ref(|&i| elements[i].clone()));
                }
                TensorSet::Tensors(out_tensors)
            }
        }
    }
}

impl<S: Clone> EvalTreeTensor<SymComplex<Rational>, S> {
    pub fn horner_scheme(&mut self) {
        self.eval.horner_scheme()
    }

    pub fn optimize_horner_scheme(
        &mut self,
        settings: &OptimizationSettings,
    ) -> Vec<Expression<SymComplex<Rational>>> {
        self.eval.optimize_horner_scheme(settings)
    }

    pub fn optimize(
        &mut self,
        settings: &OptimizationSettings,
    ) -> EvalTensor<ExpressionEvaluator<SymComplex<Rational>>, S> {
        let _ = self.optimize_horner_scheme(settings);
        self.common_subexpression_elimination();
        self.clone().linearize(settings)
    }
}

impl<S: Clone> EvalTreeTensor<SymComplex<Rational>, S> {
    pub fn from_dense(
        dense: &DenseTensor<Atom, S>,
        fn_map: &FunctionMap,
        params: &[Atom],
    ) -> Result<Self, String> {
        let atomviews: Vec<AtomView> = dense.data.iter().map(|a| a.as_view()).collect();
        let eval = AtomView::to_eval_tree_multiple(&atomviews, fn_map, params)?;

        Ok(EvalTreeTensor {
            eval,
            indexmap: None,
            structure: dense.structure.clone(),
        })
    }

    pub fn from_sparse(
        dense: &SparseTensor<Atom, S>,
        fn_map: &FunctionMap,
        params: &[Atom],
    ) -> Result<Self, String> {
        let atomviews: (Vec<FlatIndex>, Vec<AtomView>) = dense
            .elements
            .iter()
            .map(|(k, a)| (*k, a.as_view()))
            .unzip();
        let eval = AtomView::to_eval_tree_multiple(&atomviews.1, fn_map, params)?;

        Ok(EvalTreeTensor {
            eval,
            indexmap: Some(atomviews.0),
            structure: dense.structure.clone(),
        })
    }

    pub fn from_data(
        data: &DataTensor<Atom, S>,
        fn_map: &FunctionMap,
        params: &[Atom],
    ) -> Result<Self, String>
    where
        S: TensorStructure,
    {
        match data {
            DataTensor::Dense(d) => Self::from_dense(d, fn_map, params),
            DataTensor::Sparse(s) => Self::from_sparse(s, fn_map, params),
        }
    }
}

impl<S: Clone> EvalTreeTensor<SymComplex<Rational>, S> {
    pub fn linearize(
        mut self,
        settings: &OptimizationSettings,
    ) -> EvalTensor<ExpressionEvaluator<SymComplex<Rational>>, S> {
        EvalTensor {
            eval: self.eval.linearize(settings),
            structure: self.structure,
            indexmap: self.indexmap,
        }
    }
}

impl<S: Clone, T> EvalTreeTensor<T, S> {
    pub fn map_coeff<T2, F: Fn(&T) -> T2>(&self, f: &F) -> EvalTreeTensor<T2, S>
    where
        T: Clone + PartialEq,
    {
        EvalTreeTensor {
            eval: self.eval.map_coeff(f),
            indexmap: self.indexmap.clone(),
            structure: self.structure.clone(),
        }
        // self.map_data_ref(|x| x.map_coeff(f))
    }

    pub fn common_subexpression_elimination(&mut self)
    where
        T: Debug + Hash + Eq + InternalOrdering + Clone + Default,
    {
        self.eval.common_subexpression_elimination()
    }

    pub fn evaluate(&mut self, params: &[T]) -> DataTensor<T, S>
    where
        T: Real,
        S: TensorStructure,
    {
        let zero = params[0].zero();
        if let Some(ref indexmap) = self.indexmap {
            let mut elements = vec![zero.clone(); indexmap.len()];
            self.eval.evaluate(params, &mut elements);
            let s = SparseTensor {
                zero: zero.clone(),
                elements: indexmap.iter().cloned().zip(elements.drain(0..)).collect(),
                structure: self.structure.clone(),
            };
            DataTensor::Sparse(s)
        } else {
            let mut out_data = DenseTensor::repeat(self.structure.clone(), zero.clone());
            self.eval.evaluate(params, &mut out_data.data);
            DataTensor::Dense(out_data)
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct EvalTensor<T, S> {
    eval: T,
    indexmap: Option<Vec<FlatIndex>>,
    structure: S,
}

impl<T, S: TensorStructure + Clone> EvalTensor<T, S> {
    pub fn usize_tensor(&self, shift: usize) -> DataTensor<usize, S> {
        if let Some(ref indexmap) = self.indexmap {
            let mut sparse_tensor = SparseTensor::empty(self.structure.clone(), 0);
            for (i, idx) in indexmap.iter().enumerate() {
                sparse_tensor.elements.insert(*idx, shift + i);
            }
            DataTensor::Sparse(sparse_tensor)
        } else {
            let data: Vec<usize> = (shift..shift + self.structure.size().unwrap()).collect();
            let dense_tensor = DenseTensor::from_data(data, self.structure.clone()).unwrap();
            DataTensor::Dense(dense_tensor)
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub enum TensorsOrScalars<T, S: TensorStructure> {
    Tensors(Vec<DataTensor<T, S>>),
    Scalars,
}

impl<T, S: TensorStructure> TensorsOrScalars<T, S> {
    pub fn push(&mut self, tensor: DataTensor<T, S>) {
        if tensor.is_scalar() {
            if let TensorsOrScalars::Tensors(_) = self {
                panic!("Trying to push a scalar to a list of tensors")
            }
        } else {
            match self {
                TensorsOrScalars::Tensors(t) => t.push(tensor),
                TensorsOrScalars::Scalars => {
                    panic!("Trying to push a tensor to a list of scalars")
                }
            }
        }
    }

    pub fn len(&self) -> usize {
        match self {
            TensorsOrScalars::Tensors(t) => t.len(),
            TensorsOrScalars::Scalars => 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        match self {
            TensorsOrScalars::Tensors(t) => t.is_empty(),
            TensorsOrScalars::Scalars => true,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct EvalTensorSet<T, S: TensorStructure> {
    pub tensors: TensorsOrScalars<usize, S>,
    eval: T,
    size: usize, //
}

impl<S: TensorStructure, T> EvalTensorSet<T, S> {
    pub fn len(&self) -> usize {
        match &self.tensors {
            TensorsOrScalars::Tensors(t) => t.len(),
            TensorsOrScalars::Scalars => self.size,
        }
    }

    pub fn is_empty(&self) -> bool {
        match &self.tensors {
            TensorsOrScalars::Tensors(t) => t.is_empty(),
            TensorsOrScalars::Scalars => self.size == 0,
        }
    }
}

impl<T, S> TensorStructure for EvalTensor<T, S>
where
    S: TensorStructure,
{
    // type R = <T::Structure as TensorStructure>::R;
    type Indexed = EvalTensor<T, S::Indexed>;
    type Slot = S::Slot;

    fn reindex(
        self,
        indices: &[<Self::Slot as IsAbstractSlot>::Aind],
    ) -> Result<PermutedStructure<Self::Indexed>, StructureError> {
        let res = self.structure.reindex(indices)?;
        Ok(PermutedStructure {
            structure: EvalTensor {
                eval: self.eval,
                indexmap: self.indexmap,
                structure: res.structure,
            },
            rep_permutation: res.rep_permutation,
            index_permutation: res.index_permutation,
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

impl<T, S: TensorStructure> HasStructure for EvalTensor<T, S> {
    type Scalar = T;
    type ScalarRef<'a>
        = &'a T
    where
        Self: 'a;
    type Structure = S;
    type Store<U>
        = EvalTensor<T, U>
    where
        U: TensorStructure;

    fn map_structure<O: TensorStructure>(self, f: impl Fn(Self::Structure) -> O) -> Self::Store<O> {
        EvalTensor {
            structure: f(self.structure),
            indexmap: self.indexmap,
            eval: self.eval,
        }
    }

    fn map_structure_result<O: TensorStructure, Er>(
        self,
        f: impl Fn(Self::Structure) -> Result<O, Er>,
    ) -> std::result::Result<Self::Store<O>, Er> {
        Ok(EvalTensor {
            structure: f(self.structure)?,
            indexmap: self.indexmap,
            eval: self.eval,
        })
    }

    fn structure(&self) -> &Self::Structure {
        &self.structure
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        &mut self.structure
    }

    fn map_same_structure(self, f: impl FnOnce(Self::Structure) -> Self::Structure) -> Self {
        Self {
            eval: self.eval,
            indexmap: self.indexmap,
            structure: f(self.structure),
        }
    }

    fn scalar(self) -> Option<Self::Scalar> {
        if self.is_scalar() {
            Some(self.eval)
        } else {
            None
        }
    }

    fn scalar_ref(&self) -> Option<Self::ScalarRef<'_>> {
        if self.is_scalar() {
            Some(&self.eval)
        } else {
            None
        }
    }
}

impl<S: TensorStructure + Clone>
    EvalTensorSet<
        (
            ExpressionEvaluator<SymComplex<Rational>>,
            Option<Vec<Expression<SymComplex<Rational>>>>,
        ),
        S,
    >
{
    pub fn push_optimize(
        &mut self,
        mut tensor: EvalTreeTensor<SymComplex<Rational>, S>,
        settings: &OptimizationSettings,
    ) {
        let usize_tensor = tensor.usize_tensor(self.size);
        trace!("adding a tensor to the list of {} tensors", self.len());
        self.size += usize_tensor.actual_size();
        self.tensors.push(usize_tensor);
        self.eval.1 = Some(tensor.optimize_horner_scheme(settings));
        tensor.common_subexpression_elimination();

        self.eval
            .0
            .merge(tensor.linearize(settings).eval, None)
            .unwrap();
    }
}

pub type LinearizedEvalTensor<T, S> = EvalTensor<ExpressionEvaluator<T>, S>;

impl<T, S> EvalTensor<ExpressionEvaluator<T>, S> {
    pub fn map_coeff<T2, F: Fn(&T) -> T2>(self, f: &F) -> LinearizedEvalTensor<T2, S>
    where
        T: Clone + PartialEq + Default,
        S: Clone,
    {
        LinearizedEvalTensor {
            eval: self.eval.map_coeff(f),
            indexmap: self.indexmap,
            structure: self.structure,
        }
        // self.map_data_ref(|x| x.map_coeff(f))
    }
    pub fn export_cpp<F: CompiledNumber>(
        &self,
        path: impl AsRef<Path>,
        function_name: &str,
        settings: ExportSettings,
    ) -> Result<EvalTensor<ExportedCode<F>, S>, std::io::Error>
    where
        T: ExportNumber + SingleFloat,
        S: Clone,
    {
        Ok(EvalTensor {
            eval: self.eval.export_cpp(path, function_name, settings)?,
            indexmap: self.indexmap.clone(),
            structure: self.structure.clone(),
        })
    }

    pub fn evaluate(&mut self, params: &[T]) -> DataTensor<T, S>
    where
        T: Real,
        S: TensorStructure + Clone,
    {
        let zero = params[0].zero();
        if let Some(ref indexmap) = self.indexmap {
            let mut elements = vec![zero.clone(); indexmap.len()];
            self.eval.evaluate(params, &mut elements);
            let s = SparseTensor {
                zero: zero.clone(),
                elements: indexmap.iter().cloned().zip(elements.drain(0..)).collect(),
                structure: self.structure.clone(),
            };
            DataTensor::Sparse(s)
        } else {
            let mut out_data = DenseTensor::repeat(self.structure.clone(), zero.clone());
            self.eval.evaluate(params, &mut out_data.data);
            DataTensor::Dense(out_data)
        }
    }
}

impl<S: TensorStructure, F: CompiledNumber> EvalTensor<ExportedCode<F>, S> {
    pub fn compile(
        &self,
        out: &str,
        options: CompileOptions,
    ) -> Result<EvalTensor<CompiledCode<F>, S>, std::io::Error>
    where
        S: Clone,
    {
        Ok(EvalTensor {
            eval: self.eval.compile(out, options)?,
            indexmap: self.indexmap.clone(),
            structure: self.structure.clone(),
        })
    }
}

impl<S: TensorStructure, F: CompiledNumber> EvalTensor<CompiledCode<F>, S> {
    pub fn load(&self) -> Result<EvalTensor<F::Evaluator, S>, String>
    where
        S: Clone,
    {
        Ok(EvalTensor {
            eval: self.eval.load()?,
            indexmap: self.indexmap.clone(),
            structure: self.structure.clone(),
        })
    }
}

pub type LinearizedEvalTensorSet<T, S> =
    EvalTensorSet<(ExpressionEvaluator<T>, Option<Vec<Expression<T>>>), S>;

impl<T, S: TensorStructure> LinearizedEvalTensorSet<T, S> {
    pub fn map_coeff<T2, F: Fn(&T) -> T2>(self, f: &F) -> LinearizedEvalTensorSet<T2, S>
    where
        T: Clone + PartialEq + Default,
        S: Clone,
    {
        LinearizedEvalTensorSet {
            eval: (self.eval.0.map_coeff(f), None),
            tensors: self.tensors,
            size: self.size,
        }
        // self.map_data_ref(|x| x.map_coeff(f))
    }

    pub fn export_cpp<F: CompiledNumber>(
        &self,
        path: impl AsRef<Path>,
        function_name: &str,
        settings: ExportSettings,
    ) -> Result<EvalTensorSet<ExportedCode<F>, S>, std::io::Error>
    where
        T: ExportNumber + SingleFloat,
        S: Clone,
    {
        Ok(EvalTensorSet {
            eval: self.eval.0.export_cpp(path, function_name, settings)?,
            tensors: self.tensors.clone(),
            size: self.size,
        })
    }

    pub fn evaluate(&mut self, params: &[T]) -> TensorSet<DataTensor<T, S>>
    where
        T: Real,
        S: TensorStructure + Clone,
    {
        let zero = params[0].zero();

        let mut elements = vec![zero; self.size];
        self.eval.0.evaluate(params, &mut elements);

        match &self.tensors {
            TensorsOrScalars::Scalars => TensorSet::Scalars(elements),
            TensorsOrScalars::Tensors(t) => {
                let mut out_tensors = Vec::with_capacity(t.len());
                trace!("Evaluating {} tensors", t.len());
                for t in t.iter() {
                    out_tensors.push(t.map_data_ref(|&i| elements[i].clone()));
                }
                TensorSet::Tensors(out_tensors)
            }
        }
    }
}

impl<S: TensorStructure, F: CompiledNumber> EvalTensorSet<ExportedCode<F>, S> {
    pub fn compile(
        &self,
        out: impl AsRef<Path>,
        options: CompileOptions,
    ) -> Result<EvalTensorSet<CompiledCode<F>, S>, std::io::Error>
    where
        S: Clone,
    {
        Ok(EvalTensorSet {
            eval: self.eval.compile(out, options)?,
            tensors: self.tensors.clone(),
            size: self.size,
        })
    }
}

impl<S: TensorStructure, F: CompiledNumber> EvalTensorSet<CompiledCode<F>, S> {
    pub fn load(&self) -> Result<EvalTensorSet<F::Evaluator, S>, String>
    where
        S: Clone,
    {
        Ok(EvalTensorSet {
            eval: self.eval.load()?,
            tensors: self.tensors.clone(),
            size: self.size,
        })
    }
}

pub type CompiledEvalTensor<S> = EvalTensor<CompiledComplexEvaluatorSpenso, S>;

impl<S> EvalTensor<CompiledComplexEvaluatorSpenso, S> {
    pub fn evaluate(&mut self, params: &[Complex<f64>]) -> DataTensor<Complex<f64>, S>
    where
        S: TensorStructure + Clone,
    {
        if let Some(ref indexmap) = self.indexmap {
            let mut elements: Vec<Complex<f64>> = vec![Complex::default(); indexmap.len()];
            self.eval.evaluate(params, &mut elements);
            let s = SparseTensor {
                zero: Complex::new_zero(),
                elements: indexmap.iter().cloned().zip(elements.drain(0..)).collect(),
                structure: self.structure.clone(),
            };
            DataTensor::Sparse(s)
        } else {
            let mut out_data = DenseTensor::repeat(self.structure.clone(), Complex::default());
            self.eval.evaluate(params, &mut out_data.data);
            DataTensor::Dense(out_data)
        }
    }
}

impl<S> EvalTensor<CompiledComplexEvaluator, S> {
    pub fn evaluate(&mut self, params: &[SymComplex<f64>]) -> DataTensor<SymComplex<f64>, S>
    where
        S: TensorStructure + Clone,
    {
        if let Some(ref indexmap) = self.indexmap {
            let mut elements: Vec<SymComplex<f64>> = vec![SymComplex::default(); indexmap.len()];
            self.eval.evaluate(params, &mut elements);
            let s = SparseTensor {
                zero: SymComplex::new_zero(),
                elements: indexmap.iter().cloned().zip(elements.drain(0..)).collect(),
                structure: self.structure.clone(),
            };
            DataTensor::Sparse(s)
        } else {
            let mut out_data = DenseTensor::repeat(self.structure.clone(), SymComplex::default());
            self.eval.evaluate(params, &mut out_data.data);
            DataTensor::Dense(out_data)
        }
    }
}

pub type CompiledEvalTensorSet<S> = EvalTensorSet<CompiledComplexEvaluatorSpenso, S>;

impl<S: TensorStructure> CompiledEvalTensorSet<S> {
    pub fn evaluate(&mut self, params: &[Complex<f64>]) -> TensorSet<DataTensor<Complex<f64>, S>>
    where
        S: TensorStructure + Clone,
    {
        let zero = params[0].zero();

        let mut elements = vec![zero; self.size];
        self.eval.evaluate(params, &mut elements);

        match &self.tensors {
            TensorsOrScalars::Scalars => TensorSet::Scalars(elements),
            TensorsOrScalars::Tensors(t) => {
                let mut out_tensors = Vec::with_capacity(t.len());
                for t in t.iter() {
                    out_tensors.push(t.map_data_ref(|&i| elements[i]));
                }
                TensorSet::Tensors(out_tensors)
            }
        }
    }
}

// impl FallibleMul for Atom {
//     type Output = Atom;
//     fn mul_fallible(&self, rhs: &Self) -> Option<Self::Output> {
//         Some(self * rhs)
//     }
// }

// impl FallibleAdd<Atom> for Atom {
//     type Output = Atom;
//     fn add_fallible(&self, rhs: &Self) -> Option<Self::Output> {
//         Some(self + rhs)
//     }
// }

#[cfg(test)]
pub mod test {
    use symbolica::atom::Atom;

    use crate::{
        structure::{
            OrderedStructure, PermutedStructure, TensorStructure,
            representation::{Minkowski, RepName},
        },
        tensors::data::{DataTensor, SparseTensor},
    };

    use super::MixedTensor;

    #[test]
    fn tensor_structure() {
        let a: MixedTensor<f64, OrderedStructure> =
            MixedTensor::param(DataTensor::Sparse(SparseTensor::empty(
                PermutedStructure::from_iter([Minkowski {}.new_slot(2, 1)]).structure,
                Atom::Zero,
            )));

        assert_eq!(a.size().unwrap(), 2);
    }
}
