use linnet::permutation::Permutation;
use tabled::{builder::Builder, settings::Style};

use super::{
    abstract_index::AbstractIndex,
    dimension::Dimension,
    permuted::PermuteTensor,
    representation::{LibraryRep, RepName, Representation},
    slot::{AbsInd, DummyAind, IsAbstractSlot, Slot},
    HasName, MergeInfo, OrderedStructure, PermutedStructure, ScalarStructure, StructureContract,
    StructureError, TensorStructure,
};
use bitvec::vec::BitVec;

use eyre::Result;
use delegate::delegate;

#[cfg(feature = "shadowing")]
use crate::shadowing::symbolica_utils::{SerializableAtom, SerializableSymbol};

#[cfg(not(feature = "shadowing"))]
use serde::{Deserialize, Serialize};
/// A named structure is a structure with a global name, and a list of slots
///
/// It is useful when you want to shadow tensors, to nest tensor network contraction operations.
#[derive(
    Clone, PartialEq, Eq, Default, Hash, bincode_trait_derive::Encode, bincode_trait_derive::Decode,
)]
#[cfg_attr(not(feature = "shadowing"), derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "shadowing",
    trait_decode(trait = symbolica::state::HasStateMap),
)]
pub struct NamedStructure<
    Name = String,
    Args = usize,
    R: RepName = LibraryRep,
    Aind = AbstractIndex,
> {
    pub structure: OrderedStructure<R, Aind>,
    pub global_name: Option<Name>,
    pub additional_args: Option<Args>,
}

#[cfg(feature = "shadowing")]
pub type AtomStructure<R, A> = NamedStructure<SerializableSymbol, Vec<SerializableAtom>, R, A>;

impl<Name, Args, R: RepName, Aind: AbsInd> NamedStructure<Name, Args, R, Aind> {
    #[must_use]
    pub fn from_iter<I, T>(iter: T, name: Name, args: Option<Args>) -> PermutedStructure<Self>
    where
        R: From<I>,
        I: RepName,
        T: IntoIterator<Item = Slot<I, Aind>>,
    {
        iter.into_iter()
            .map(|a| a.cast())
            .collect::<PermutedStructure<_>>()
            .map_structure(move |structure| Self {
                structure,
                global_name: Some(name),
                additional_args: args,
            })
    }
}
impl<N, A, R: RepName, Aind: AbsInd> From<OrderedStructure<R, Aind>>
    for NamedStructure<N, A, R, Aind>
{
    fn from(value: OrderedStructure<R, Aind>) -> Self {
        Self {
            structure: value,
            global_name: None,
            additional_args: None,
        }
    }
}

/// A trait for a structure that has a name
impl<N, A, R: RepName, Aind: AbsInd> HasName for NamedStructure<N, A, R, Aind>
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

impl<N, A, R: RepName, Aind: AbsInd> ScalarStructure for NamedStructure<N, A, R, Aind>
where
    N: Clone,
    A: Clone,
{
    fn scalar_structure() -> Self {
        NamedStructure {
            structure: OrderedStructure::default(),
            global_name: None,
            additional_args: None,
        }
    }
}

pub trait IdentityName {
    fn id() -> Self;
}
pub const METRIC_NAME: &str = "g";
impl IdentityName for String {
    fn id() -> Self {
        METRIC_NAME.to_string()
    }
}

impl<N: IdentityName, A, R: RepName<Dual = R>, Aind: AbsInd + DummyAind> PermuteTensor
    for NamedStructure<N, A, R, Aind>
{
    type Id = Self;
    type IdSlot = Slot<R, Aind>;
    type Permuted = (
        NamedStructure<N, A, LibraryRep, Aind>,
        Vec<NamedStructure<N, A, LibraryRep, Aind>>,
    );

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
        Self {
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
            ids.push(NamedStructure::id(d, ogs));
        }
        let strct = OrderedStructure::new(dummy_structure);
        if !strct.index_permutation.is_identity() {
            panic!("should be identity")
        }

        (
            NamedStructure {
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
                NamedStructure {
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

            ids.push(NamedStructure::id(d, new_slot));
        }
        let strct = OrderedStructure::new(dummy_structure);
        if !strct.index_permutation.is_identity() {
            panic!("should be identity")
        }
        (
            NamedStructure {
                global_name: self.global_name,
                additional_args: self.additional_args,
                structure: strct.structure,
            },
            ids,
        )
    }
}

impl<N, A, R: RepName<Dual = R>, Aind: AbsInd> TensorStructure for NamedStructure<N, A, R, Aind> {
    type Slot = Slot<R, Aind>;
    // type R = PhysicalReps;
    type Indexed = Self;

    fn is_fully_self_dual(&self) -> bool {
        self.structure.is_fully_self_dual()
    }

    // fn id(i: Self::Slot, j: Self::Slot) -> Self {
    //     Self {
    //         structure: OrderedStructure::id(i, j),
    //         global_name: Some(N::id()),
    //         additional_args: None,
    //     }
    // }

    fn reindex(self, indices: &[Aind]) -> Result<PermutedStructure<Self::Indexed>, StructureError> {
        let res = self.structure.reindex(indices)?;

        Ok(PermutedStructure {
            rep_permutation: res.rep_permutation,
            structure: Self {
                global_name: self.global_name,
                additional_args: self.additional_args,
                structure: res.structure,
            },
            index_permutation: res.index_permutation,
        })
    }

    fn dual(self) -> Self {
        NamedStructure {
            structure: self.structure.dual(),
            global_name: self.global_name,
            additional_args: self.additional_args,
        }
    }
    delegate! {
        to self.structure{
            fn external_reps_iter(&self) -> impl Iterator<Item = Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn external_indices_iter(&self) -> impl Iterator<Item = Aind>;
            fn external_dims_iter(&self)->impl Iterator<Item=Dimension>;
            fn external_structure_iter(&self) -> impl Iterator<Item = Self::Slot>;
            fn order(&self) -> usize;
            fn get_slot(&self, i: usize) -> Option<Self::Slot>;
            fn get_rep(&self, i: usize) -> Option<Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn get_aind(&self,i:usize)->Option<Aind>;
            fn get_dim(&self, i: usize) -> Option<Dimension>;
        }
    }
}

pub trait ArgDisplay {
    fn arg_display(&self) -> String;

    fn arg_debug(&self) -> String;
}

pub trait ArgDisplayMarker {}

impl<T: std::fmt::Display + ArgDisplayMarker + std::fmt::Debug> ArgDisplay for T {
    fn arg_display(&self) -> String {
        self.to_string()
    }

    fn arg_debug(&self) -> String {
        format!("{:?}", self)
    }
}

impl ArgDisplayMarker for () {}
impl ArgDisplayMarker for usize {}
impl ArgDisplayMarker for isize {}
impl ArgDisplayMarker for f64 {}
impl ArgDisplayMarker for f32 {}
impl ArgDisplayMarker for i8 {}
impl ArgDisplayMarker for i16 {}
impl ArgDisplayMarker for i32 {}
impl ArgDisplayMarker for i64 {}

impl ArgDisplayMarker for u8 {}
impl ArgDisplayMarker for u16 {}
impl ArgDisplayMarker for u32 {}
impl ArgDisplayMarker for u64 {}

#[cfg(feature = "shadowing")]
impl ArgDisplayMarker for symbolica::atom::Atom {}
#[cfg(feature = "shadowing")]
impl ArgDisplayMarker for symbolica::atom::Symbol {}

#[cfg(feature = "shadowing")]
impl ArgDisplay for Vec<symbolica::atom::Atom> {
    fn arg_display(&self) -> String {
        self.iter()
            .map(|a| a.arg_display())
            .collect::<Vec<String>>()
            .join(", ")
    }

    fn arg_debug(&self) -> String {
        self.iter()
            .map(|a| a.arg_debug())
            .collect::<Vec<String>>()
            .join(", ")
    }
}

impl<N: std::fmt::Display, A: ArgDisplay, R: RepName, Aind: AbsInd> std::fmt::Display
    for NamedStructure<N, A, R, Aind>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table = Builder::new();

        table.push_record(&[
            self.global_name
                .as_ref()
                .map(|a| format!("{a}"))
                .unwrap_or("NO NAME".to_string()),
            self.additional_args
                .as_ref()
                .map(|a| a.arg_display())
                .unwrap_or("".to_string()),
        ]);
        for item in self.structure.structure.iter() {
            if item.rep.rep.is_self_dual() {
                table.push_record(&[item.rep.to_string(), format!("{}", item.aind)]);
            } else if item.rep.rep.is_base() {
                table.push_record(&[item.rep.to_string(), format!("{:+}", item.aind)]);
            } else {
                table.push_record(&[item.rep.to_string(), format!("{:-}", item.aind)]);
            }
        }
        writeln!(f)?;
        table.build().with(Style::rounded()).fmt(f)
    }
}
impl<N: std::fmt::Debug, A: ArgDisplay, R: RepName, Aind: AbsInd> std::fmt::Debug
    for NamedStructure<N, A, R, Aind>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table = Builder::new();

        table.push_record(&[
            self.global_name
                .as_ref()
                .map(|a| format!("{:?}", a))
                .unwrap_or("NO NAME".to_string()),
            self.additional_args
                .as_ref()
                .map(|a| a.arg_debug())
                .unwrap_or("".to_string()),
        ]);
        for (index, item) in self.structure.structure.iter().enumerate() {
            table.push_record(&[
                index.to_string(),
                format!("{:?}", item.rep.rep),
                format!("{:?}", item.rep.dim),
                format!("{:?}", item.aind),
            ]);
        }
        writeln!(f)?;
        write!(f, "{}", table.build().with(Style::rounded()))
    }
}
impl<N, A, R: RepName<Dual = R>, Aind: AbsInd> StructureContract for NamedStructure<N, A, R, Aind> {
    delegate! {
        to self.structure{
            fn trace_out(&mut self);
            fn trace(&mut self, i: usize, j: usize);

        }
    }

    fn merge(&self, other: &Self) -> Result<(Self, BitVec, BitVec, MergeInfo), StructureError> {
        let (structure, pos_self, pos_other, mergeinfo) = self.structure.merge(&other.structure)?;
        Ok((
            Self {
                structure,
                global_name: None,
                additional_args: None,
            },
            pos_self,
            pos_other,
            mergeinfo,
        ))
    }
}
