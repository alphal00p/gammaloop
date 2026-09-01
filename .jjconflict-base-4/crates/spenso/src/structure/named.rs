use linnet::half_edge::subgraph::subset::SubSet;
use tabled::{builder::Builder, settings::Style};

use super::{
    ApplyPendingIndexPermutation, Canonicalized, HasName, MergeInfo, OrderedStructure,
    PendingIndexPermutation, Reindexed, ScalarStructure, StructureContract, StructureError,
    TensorIdentity, TensorStructure,
    abstract_index::AbstractIndex,
    dimension::Dimension,
    representation::{LibraryRep, RepName, Representation},
    slot::{AbsInd, IsAbstractSlot, Slot},
};

use delegate::delegate;

use crate::structure::SlotIndex;
#[cfg(feature = "shadowing")]
use symbolica::atom::Atom;
#[cfg(feature = "shadowing")]
use symbolica_utils::SerializableSymbol;

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
pub type AtomStructure<R, A> = NamedStructure<SerializableSymbol, Vec<Atom>, R, A>;

impl<Name, Args, R: RepName, Aind: AbsInd> NamedStructure<Name, Args, R, Aind> {
    #[must_use]
    pub fn from_iter<I, T>(iter: T, name: Name, args: Option<Args>) -> Canonicalized<Self>
    where
        R: From<I>,
        I: RepName,
        T: IntoIterator<Item = Slot<I, Aind>>,
    {
        iter.into_iter()
            .map(|a| a.cast())
            .collect::<Canonicalized<_>>()
            .map_canonical(move |structure| Self {
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

impl<N: IdentityName, A, R: RepName<Dual = R>, Aind: AbsInd + super::slot::DummyAind> TensorIdentity
    for NamedStructure<N, A, R, Aind>
{
    type Id = Self;
    type IdSlot = Slot<R, Aind>;

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
        Self {
            structure: OrderedStructure::id(i, j),
            global_name: Some(N::id()),
            additional_args: None,
        }
    }
}

impl<N, A, R: RepName, Aind> ApplyPendingIndexPermutation for NamedStructure<N, A, R, Aind> {
    type Output = Self;

    fn apply_pending_index_permutation(self, _pending: &PendingIndexPermutation) -> Self::Output {
        self
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

    fn reindex_storage(self, indices: &[Aind]) -> Result<Reindexed<Self::Indexed>, StructureError> {
        self.structure.reindex_storage(indices).map(|reindexed| {
            reindexed.map_target(|structure| Self {
                global_name: self.global_name,
                additional_args: self.additional_args,
                structure,
            })
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
            fn get_slot(&self,i: impl Into<SlotIndex>) -> Option<Self::Slot>;
            fn get_rep(&self,i: impl Into<SlotIndex>) -> Option<Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn get_aind(&self,i: impl Into<SlotIndex>)->Option<Aind>;
            fn get_dim(&self,i: impl Into<SlotIndex>) -> Option<Dimension>;
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

    fn merge(
        &self,
        other: &Self,
    ) -> Result<(Self, SubSet<SlotIndex>, SubSet<SlotIndex>, MergeInfo), StructureError> {
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
