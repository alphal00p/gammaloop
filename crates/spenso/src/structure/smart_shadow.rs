use linnet::permutation::Permutation;

use delegate::delegate;

use super::{
    HasName, MergeInfo, NamedStructure, OrderedStructure, PermutedStructure, ScalarStructure,
    StructureContract, StructureError, TensorStructure, TracksCount,
    abstract_index::AbstractIndex,
    dimension::Dimension,
    named::IdentityName,
    permuted::PermuteTensor,
    representation::{LibraryRep, RepName, Representation},
    slot::{AbsInd, DummyAind, IsAbstractSlot, Slot},
};
use bitvec::vec::BitVec;

use eyre::Result;

#[cfg(feature = "shadowing")]
use crate::shadowing::symbolica_utils::{IntoArgs, IntoSymbol};

#[cfg(not(feature = "shadowing"))]
use serde::{Deserialize, Serialize};

/// A structure to enable smart shadowing of tensors in a tensor network contraction algorithm.
#[derive(
    Clone, PartialEq, Eq, Debug, Hash, bincode_trait_derive::Encode, bincode_trait_derive::Decode,
)]
#[cfg_attr(not(feature = "shadowing"), derive(Serialize, Deserialize))]
#[cfg_attr(
feature = "shadowing",
trait_decode(trait = symbolica::state::HasStateMap),
)]
pub struct SmartShadowStructure<
    Name = String,
    Args = usize,
    R: RepName = LibraryRep,
    Aind: AbsInd = AbstractIndex,
> {
    pub structure: OrderedStructure<R, Aind>,
    pub contractions: usize,
    pub global_name: Option<Name>,
    additional_args: Option<Args>,
}

impl<Name, Args, R: RepName, Aind: AbsInd> SmartShadowStructure<Name, Args, R, Aind> {
    pub fn map_underlying_structure(
        self,
        mut f: impl FnMut(OrderedStructure<R, Aind>) -> OrderedStructure<R, Aind>,
    ) -> Self {
        Self {
            structure: f(self.structure),
            contractions: self.contractions,
            global_name: self.global_name,
            additional_args: self.additional_args,
        }
    }

    /// Constructs a new [`SmartShadow`] from a list of tuples of indices and dimension (assumes they are all euclidean), along with a name
    #[must_use]
    pub fn from_iter<I, T>(
        iter: T,
        name: Option<Name>,
        args: Option<Args>,
    ) -> PermutedStructure<Self>
    where
        I: Into<Slot<R, Aind>>,
        T: IntoIterator<Item = I>,
    {
        let res = iter
            .into_iter()
            .map(I::into)
            .collect::<PermutedStructure<_>>();
        PermutedStructure {
            structure: Self {
                structure: res.structure,
                global_name: name,
                additional_args: args,
                contractions: 0,
            },
            rep_permutation: res.rep_permutation,
            index_permutation: res.index_permutation,
        }
    }
}

impl<N, A, R: RepName, Aind: AbsInd> HasName for SmartShadowStructure<N, A, R, Aind>
where
    N: Clone,
    A: Clone,
{
    type Name = N;
    type Args = A;

    fn name(&self) -> Option<Self::Name> {
        self.global_name.clone()
    }
    fn set_name(&mut self, name: Self::Name) {
        self.global_name = Some(name);
    }
    fn args(&self) -> Option<Self::Args> {
        self.additional_args.clone()
    }
}

#[cfg(feature = "shadowing")]
impl<N: IntoSymbol, A: IntoArgs, R: RepName, Aind: AbsInd> std::fmt::Display
    for SmartShadowStructure<N, A, R, Aind>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(ref name) = self.global_name {
            write!(f, "{}", name.ref_into_symbol())?
        }
        write!(f, "(")?;
        if let Some(ref args) = self.additional_args {
            let args: Vec<std::string::String> =
                args.ref_into_args().map(|s| s.to_string()).collect();
            write!(f, "{},", args.join(","))?
        }

        write!(f, "{})", self.structure)?;
        Result::Ok(())
    }
}

impl<N, A, R: RepName, Aind: AbsInd> ScalarStructure for SmartShadowStructure<N, A, R, Aind> {
    fn scalar_structure() -> Self {
        SmartShadowStructure {
            structure: OrderedStructure::default(),
            contractions: 0,
            global_name: None,
            additional_args: None,
        }
    }
}

impl<N: IdentityName, A, R: RepName<Dual = R>, Aind: AbsInd + DummyAind> PermuteTensor
    for SmartShadowStructure<N, A, R, Aind>
{
    type Id = Self;
    type IdSlot = Slot<R, Aind>;
    type Permuted = (
        SmartShadowStructure<N, A, LibraryRep, Aind>,
        Vec<SmartShadowStructure<N, A, LibraryRep, Aind>>,
    );

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
        Self {
            contractions: 0,
            structure: OrderedStructure::id(i, j),
            global_name: Some(N::id()),
            additional_args: None,
        }
    }

    fn permute_inds(self, permutation: &Permutation) -> Self::Permuted {
        let mut dummy_structure = Vec::new();
        let mut ids = Vec::new();

        for s in permutation.iter_slice_inv(&self.structure.structure) {
            let d = s.to_dummy_rep();
            let ogs = s.to_lib();
            dummy_structure.push(d);
            ids.push(SmartShadowStructure::id(d, ogs));
        }
        let strct = OrderedStructure::new(dummy_structure);
        if !strct.index_permutation.is_identity() {
            panic!("should be identity")
        }

        (
            SmartShadowStructure {
                contractions: self.contractions,
                global_name: self.global_name,
                additional_args: self.additional_args,
                structure: strct.structure,
            },
            ids,
        )
    }

    fn permute_reps(self, rep_perm: &Permutation) -> Self::Permuted {
        let mut dummy_structure = Vec::new();
        let mut og_reps = Vec::new();
        let mut ids = Vec::new();

        if rep_perm.is_identity() {
            return (
                SmartShadowStructure {
                    contractions: self.contractions,
                    global_name: self.global_name,
                    additional_args: self.additional_args,
                    structure: PermutedStructure::from_iter(
                        self.structure.into_iter().map(|s| s.to_lib()),
                    )
                    .structure,
                },
                vec![],
            );
        }
        for s in rep_perm.iter_slice(&self.structure.structure) {
            og_reps.push(s.rep.to_lib());
            let d = s.to_dummy_rep();
            dummy_structure.push(d);
        }

        for (i, s) in rep_perm.iter_slice(&self.structure.structure).enumerate() {
            let d = dummy_structure[i];
            let new_slot = og_reps[i].slot(s.aind);

            ids.push(SmartShadowStructure::id(d, new_slot));
        }
        let strct = OrderedStructure::new(dummy_structure);
        if !strct.index_permutation.is_identity() {
            panic!("should be identity")
        }
        (
            SmartShadowStructure {
                contractions: self.contractions,
                global_name: self.global_name,
                additional_args: self.additional_args,
                structure: strct.structure,
            },
            ids,
        )
    }
}

impl<N, A, R: RepName<Dual = R>> TensorStructure for SmartShadowStructure<N, A, R> {
    type Slot = Slot<R>;
    type Indexed = Self;

    fn reindex(
        self,
        indices: &[AbstractIndex],
    ) -> Result<PermutedStructure<Self::Indexed>, StructureError> {
        let res = self.structure.reindex(indices)?;

        Ok(PermutedStructure {
            structure: Self {
                contractions: self.contractions,
                global_name: self.global_name,
                additional_args: self.additional_args,
                structure: res.structure,
            },
            rep_permutation: res.rep_permutation,
            index_permutation: res.index_permutation,
        })
    }
    // type R = PhysicalReps;
    //
    fn dual(self) -> Self {
        SmartShadowStructure {
            structure: self.structure.dual(),
            contractions: self.contractions,
            global_name: self.global_name,
            additional_args: self.additional_args,
        }
    }

    delegate! {
        to self.structure{
            fn is_fully_self_dual(&self) -> bool;
           fn external_reps_iter(&self) -> impl Iterator<Item = Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn external_indices_iter(&self) -> impl Iterator<Item = AbstractIndex>;
            fn external_dims_iter(&self)->impl Iterator<Item=Dimension>;
            fn external_structure_iter(&self) -> impl Iterator<Item = Self::Slot>;
            fn order(&self) -> usize;
            fn get_slot(&self, i: usize) -> Option<Self::Slot>;
            fn get_rep(&self, i: usize) -> Option<Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn get_aind(&self,i:usize)->Option<AbstractIndex>;
            fn get_dim(&self, i: usize) -> Option<Dimension>;
        }
    }
}

impl<N, A, R: RepName> TracksCount for SmartShadowStructure<N, A, R> {
    fn contractions_num(&self) -> usize {
        self.contractions
    }
}

impl<N, A, R: RepName<Dual = R>> StructureContract for SmartShadowStructure<N, A, R> {
    fn merge(&self, other: &Self) -> Result<(Self, BitVec, BitVec, MergeInfo), StructureError> {
        let contractions = self.contractions + other.contractions;
        let (structure, pos_self, pos_other, mergeinfo) = self.structure.merge(&other.structure)?;
        Ok((
            Self {
                contractions,
                structure,
                global_name: None,
                additional_args: None,
            },
            pos_self,
            pos_other,
            mergeinfo,
        ))
    }

    delegate! {
        to self.structure{
            fn trace_out(&mut self);
            fn trace(&mut self, i: usize, j: usize);
        }
    }
}

impl<N, A, R: RepName<Dual = R>> From<NamedStructure<N, A, R>> for SmartShadowStructure<N, A, R> {
    fn from(value: NamedStructure<N, A, R>) -> Self {
        Self {
            structure: value.structure,
            contractions: 0,
            global_name: value.global_name,
            additional_args: value.additional_args,
        }
    }
}
