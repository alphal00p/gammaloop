use std::marker::PhantomData;

use linnet::permutation::Permutation;
use tabled::{builder::Builder, settings::Style};

use super::{
    HasName, NamedStructure, OrderedStructure, PermutedStructure, ScalarStructure,
    SmartShadowStructure, StructureError, TensorStructure,
    abstract_index::AbstractIndex,
    dimension::Dimension,
    named::{ArgDisplay, IdentityName},
    permuted::PermuteTensor,
    representation::{LibraryRep, RepName, Representation},
    slot::{AbsInd, IsAbstractSlot, Slot},
};

use delegate::delegate;
use eyre::{Result, eyre};

use crate::network::StructureLessDisplay;
#[cfg(feature = "shadowing")]
use crate::{
    shadowing::symbolica_utils::IntoSymbol,
    structure::ToSymbolic,
    tensors::parametric::{ExpandedCoefficent, TensorCoefficient},
};

#[cfg(feature = "shadowing")]
use crate::{structure::FlatIndex, tensors::data::DenseTensor};

#[cfg(not(feature = "shadowing"))]
use serde::{Deserialize, Serialize};
#[cfg(feature = "shadowing")]
use symbolica::atom::{Atom, FunctionBuilder, Symbol};

#[derive(
    Clone, PartialEq, Eq, Default, Hash, bincode_trait_derive::Encode, bincode_trait_derive::Decode,
)]
#[cfg_attr(not(feature = "shadowing"), derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "shadowing",
    trait_decode(trait = symbolica::state::HasStateMap),
)]
// #[cfg_attr(not(feature = "shadowing"), derive(bincode::Decode))]
pub struct IndexLess<T: RepName = LibraryRep, Aind = AbstractIndex> {
    pub structure: Vec<Representation<T>>,
    ind: PhantomData<Aind>,
}

impl<T: RepName, Aind> StructureLessDisplay for IndexLess<T, Aind> {
    fn display(&self) -> String {
        String::new()
    }
}

impl<R: RepName<Dual = R>, Aind: AbsInd> PermuteTensor for IndexLess<R, Aind> {
    type Id = Self;
    type Permuted = (IndexLess<LibraryRep>, Vec<IndexLess<LibraryRep, Aind>>);
    type IdSlot = Slot<R, Aind>;

    fn permute_inds(self, _permutation: &Permutation) -> Self::Permuted {
        todo!()
    }

    fn permute_reps(self, _rep_perm: &Permutation) -> Self::Permuted {
        todo!()
    }

    fn id(i: Slot<R, Aind>, j: Slot<R, Aind>) -> Self::Id {
        if i.dim() == j.dim() {
            IndexLess::new(vec![i.rep, j.rep])
        } else {
            panic!("Not same dimension for ID")
        }
    }
}

impl<R: RepName, Aind> std::fmt::Display for IndexLess<R, Aind> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table = Builder::new();

        table.push_record(&["".to_string()]);
        for item in self.structure.iter() {
            table.push_record(&[item.rep.to_string(), item.dim.to_string()]);
        }
        writeln!(f)?;
        table.build().with(Style::rounded()).fmt(f)
    }
}

impl<R: RepName, Aind> std::fmt::Debug for IndexLess<R, Aind> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table = Builder::new();

        table.push_record(&["IndexLess".to_string()]);
        for (index, item) in self.structure.iter().enumerate() {
            table.push_record(&[
                index.to_string(),
                format!("{:?}", item.rep),
                format!("{:?}", item.dim),
            ]);
        }
        writeln!(f)?;
        write!(f, "{}", table.build().with(Style::rounded()))
    }
}

impl<R: RepName, Aind> FromIterator<Representation<R>> for PermutedStructure<IndexLess<R, Aind>> {
    fn from_iter<I: IntoIterator<Item = Representation<R>>>(iter: I) -> Self {
        let structure: Vec<_> = iter.into_iter().collect();
        structure.into()
    }
}

impl<R: RepName, Aind> From<Vec<Representation<R>>> for PermutedStructure<IndexLess<R, Aind>> {
    fn from(mut structure: Vec<Representation<R>>) -> Self {
        let permutation = Permutation::sort(&structure);
        permutation.apply_slice_in_place(&mut structure);

        PermutedStructure {
            rep_permutation: permutation,
            index_permutation: Permutation::id(structure.len()),
            structure: IndexLess {
                structure,
                ind: PhantomData,
            },
        }
    }
}

impl<R: RepName, Aind> From<OrderedStructure<R, Aind>> for IndexLess<R, Aind> {
    fn from(structure: OrderedStructure<R, Aind>) -> Self {
        IndexLess {
            structure: structure.into_iter().map(|a| a.rep).collect(),
            ind: PhantomData,
        }
    }
}

impl<N, A, R: RepName, Aind> From<NamedStructure<N, A, R, Aind>> for IndexLess<R, Aind> {
    fn from(structure: NamedStructure<N, A, R, Aind>) -> Self {
        structure.structure.into()
    }
}

impl<N, A, R: RepName, Aind> From<IndexlessNamedStructure<N, A, R, Aind>> for IndexLess<R, Aind> {
    fn from(structure: IndexlessNamedStructure<N, A, R, Aind>) -> Self {
        structure.structure
    }
}

impl<N, A, R: RepName, Aind: AbsInd> From<SmartShadowStructure<N, A, R, Aind>>
    for IndexLess<R, Aind>
{
    fn from(structure: SmartShadowStructure<N, A, R, Aind>) -> Self {
        structure.structure.into()
    }
}

impl<T: RepName, Aind: AbsInd> IndexLess<T, Aind> {
    pub fn new(structure: Vec<Representation<T>>) -> Self {
        Self {
            structure,
            ind: PhantomData,
        }
    }

    pub fn empty() -> Self {
        Self {
            structure: vec![],
            ind: PhantomData,
        }
    }

    pub fn to_indexed(self, indices: &[Aind]) -> Result<Vec<Slot<T, Aind>>, StructureError> {
        if self.structure.len() != indices.len() {
            return Err(StructureError::WrongNumberOfArguments(
                indices.len(),
                self.structure.len(),
            ));
        }

        Ok(indices
            .iter()
            .cloned()
            .zip(self.structure.iter().cloned())
            .map(|(i, r)| Representation::slot(&r, i))
            .collect())
    }
}

impl<T: RepName<Dual = T>, Aind: AbsInd> ScalarStructure for IndexLess<T, Aind> {
    fn scalar_structure() -> Self {
        Self::empty()
    }
}

impl<T: RepName<Dual = T>, Aind: AbsInd> TensorStructure for IndexLess<T, Aind> {
    type Slot = Slot<T, Aind>;
    // type R = T;
    type Indexed = OrderedStructure<T, Aind>;

    fn is_fully_self_dual(&self) -> bool {
        self.structure.iter().all(|r| r.rep.is_self_dual())
    }

    fn reindex(
        self,
        indices: &[Aind],
    ) -> Result<PermutedStructure<OrderedStructure<T, Aind>>, StructureError> {
        if self.structure.len() != indices.len() {
            return Err(StructureError::WrongNumberOfArguments(
                indices.len(),
                self.structure.len(),
            ));
        }

        Ok(indices
            .iter()
            .cloned()
            .zip(self.structure.iter().cloned())
            .map(|(i, r)| Representation::slot(&r, i))
            .collect())
    }
    fn dual(self) -> Self {
        self.structure
            .into_iter()
            .map(|r| r.dual())
            .collect::<PermutedStructure<_>>()
            .structure
    }

    fn external_reps_iter(
        &self,
    ) -> impl Iterator<Item = Representation<<Self::Slot as IsAbstractSlot>::R>> {
        self.structure.iter().copied()
    }

    fn external_dims_iter(&self) -> impl Iterator<Item = Dimension> {
        self.structure.iter().map(|s| s.dim)
    }

    fn get_aind(
        &self,
        _i: impl Into<super::SlotIndex>,
    ) -> Option<<Self::Slot as IsAbstractSlot>::Aind> {
        None
    }

    fn external_indices_iter(&self) -> impl Iterator<Item = Aind> {
        [].iter().cloned()
    }

    fn external_structure_iter(&self) -> impl Iterator<Item = Self::Slot> {
        [].iter().cloned()
    }

    fn order(&self) -> usize {
        self.structure.len()
    }

    fn get_slot(&self, _i: impl Into<super::SlotIndex>) -> Option<Self::Slot> {
        None
    }

    fn find_permutation(&self, other: &Self) -> Result<Permutation> {
        if self.order() != other.order() {
            return Err(eyre!(
                "Mismatched order: {} vs {}",
                self.order(),
                other.order()
            ));
        }
        let other_structure = &other.structure;
        let self_structure = &self.structure;

        let other_sort = Permutation::sort(other_structure);
        let self_sort = Permutation::sort(self_structure);

        if other_sort.apply_slice(other_structure) == self_sort.apply_slice(self_structure) {
            Ok(other_sort.compose(&self_sort.inverse()))
        } else {
            Err(eyre!("Mismatched structure"))
        }

        // let mut index_map = HashMap::new();
        // for (i, item) in other.structure.iter().enumerate() {
        //     index_map.entry(item).or_insert_with(Vec::new).push(i);
        // }

        // let mut permutation = Vec::with_capacity(self.order());
        // let mut used_indices = HashSet::new();
        // for item in self.structure.iter() {
        //     if let Some(indices) = index_map.get_mut(&item) {
        //         // Find an index that hasn't been used yet
        //         if let Some(&index) = indices.iter().find(|&&i| !used_indices.contains(&i)) {
        //             permutation.push(index);
        //             used_indices.insert(index);
        //         } else {
        //             // No available index for this item
        //             return Err(eyre!("No available index for {:?}", item));
        //         }
        //     } else {
        //         // Item not found in other
        //         return Err(eyre!("Item {:?} not found in other", item));
        //     }
        // }

        // Ok(permutation)
    }

    fn get_rep(
        &self,
        i: impl Into<super::SlotIndex>,
    ) -> Option<Representation<<Self::Slot as IsAbstractSlot>::R>> {
        self.structure.get(i.into().0).copied()
    }

    fn get_dim(&self, i: impl Into<super::SlotIndex>) -> Option<Dimension> {
        self.structure.get(i.into().0).map(|&r| r.into())
    }
}

#[cfg(feature = "shadowing")]
impl<T: RepName<Dual = T>, Aind: AbsInd> ToSymbolic for IndexLess<T, Aind> {
    fn concrete_atom(&self, id: FlatIndex) -> ExpandedCoefficent<()> {
        ExpandedCoefficent {
            name: None,
            index: self.co_expanded_index(id).unwrap(),
            args: None,
        }
    }

    fn to_dense_labeled<R>(
        self,
        index_to_atom: impl Fn(&Self, FlatIndex) -> R,
    ) -> Result<DenseTensor<Atom, Self>>
    where
        Self: Sized,
        R: TensorCoefficient,
    {
        let mut data = vec![];
        for index in 0..self.size()? {
            data.push(index_to_atom(&self, index.into()).to_atom().unwrap());
        }

        Ok(DenseTensor {
            data,
            structure: self,
        })
    }

    fn to_dense_labeled_complex<R>(
        self,
        index_to_atom: impl Fn(&Self, FlatIndex) -> R,
    ) -> Result<DenseTensor<Atom, Self>>
    where
        Self: Sized,
        R: TensorCoefficient,
    {
        let mut data = vec![];
        for index in 0..self.size()? {
            let re = index_to_atom(&self, index.into()).to_atom_re().unwrap();
            let im = index_to_atom(&self, index.into()).to_atom_im().unwrap();
            let i = Atom::i();
            data.push(&re + i * &im);
        }

        Ok(DenseTensor {
            data,
            structure: self,
        })
    }

    fn to_symbolic_with(&self, name: Symbol, args: &[Atom], perm: Option<Permutation>) -> Atom {
        let mut slots = self
            .external_reps_iter()
            .map(|slot| slot.to_symbolic([]))
            .collect::<Vec<_>>();
        if let Some(p) = perm {
            p.apply_slice_in_place(&mut slots);
        }
        FunctionBuilder::new(name.ref_into_symbol())
            .add_args(args)
            .add_args(&slots)
            .finish()
    }
}

#[derive(
    Clone, PartialEq, Eq, Default, Hash, bincode_trait_derive::Encode, bincode_trait_derive::Decode,
)]
#[cfg_attr(not(feature = "shadowing"), derive(Serialize, Deserialize))]
#[cfg_attr(
    feature = "shadowing",
    trait_decode(trait = symbolica::state::HasStateMap),
)]
pub struct IndexlessNamedStructure<
    Name = String,
    Args = usize,
    R: RepName = LibraryRep,
    Aind = AbstractIndex,
> {
    pub structure: IndexLess<R, Aind>,
    pub global_name: Option<Name>,
    pub additional_args: Option<Args>,
}

impl<N: IdentityName, A, R: RepName<Dual = R>, Aind: AbsInd> PermuteTensor
    for IndexlessNamedStructure<N, A, R, Aind>
{
    type Id = Self;
    type IdSlot = Slot<R, Aind>;
    type Permuted = (
        IndexlessNamedStructure<N, A, LibraryRep, Aind>,
        Vec<IndexlessNamedStructure<N, A, LibraryRep, Aind>>,
    );

    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
        Self {
            structure: IndexLess::id(i, j),
            global_name: Some(N::id()),
            additional_args: None,
        }
    }

    fn permute_inds(self, _permutation: &Permutation) -> Self::Permuted {
        todo!()
    }

    fn permute_reps(self, _rep_perm: &Permutation) -> Self::Permuted {
        todo!()
    }
}
// impl<N: IdentityName, A, R: RepName<Dual = R>> PermuteTensor for IndexlessNamedStructure<N, A, R> {
//     type Id = NamedStructure<N, A, R>;
//     type IdSlot = Slot<R>;
//     type Permuted = (
//         NamedStructure<N, A, LibraryRep>,
//         Vec<NamedStructure<N, A, LibraryRep>>,
//     );

//     fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id {
//         NamedStructure {
//             structure: OrderedStructure::id(i, j),
//             global_name: Some(N::id()),
//             additional_args: None,
//         }
//     }

//     fn permute(self, permutation: &Permutation) -> Self::Permuted {
//         let mut dummy_structure = Vec::new();
//         let mut ids = Vec::new();

//         for s in permutation.iter_slice_inv(&self.structure.structure) {
//             let dind = AbstractIndex::new_dummy();
//             let d = s.to_dummy().to_lib().slot(dind);
//             let ogs = s.to_lib().slot(dind);
//             dummy_structure.push(d);
//             ids.push(NamedStructure::id(d, ogs));
//         }
//         let strct = OrderedStructure::new(dummy_structure);
//         if !strct.permutation.is_identity() {
//             panic!("should be identity")
//         }

//         (
//             NamedStructure {
//                 global_name: self.global_name,
//                 additional_args: self.additional_args,
//                 structure: strct.structure,
//             },
//             ids,
//         )
//     }
// }

impl<Name, Args, R: RepName<Dual = R>, Aind: AbsInd> TensorStructure
    for IndexlessNamedStructure<Name, Args, R, Aind>
{
    type Slot = Slot<R, Aind>;
    type Indexed = NamedStructure<Name, Args, R, Aind>;

    fn reindex(
        self,
        indices: &[Aind],
    ) -> Result<PermutedStructure<NamedStructure<Name, Args, R, Aind>>, StructureError> {
        let res = self.structure.reindex(indices)?;

        Ok(PermutedStructure {
            rep_permutation: Permutation::id(res.structure.order()),
            structure: NamedStructure {
                global_name: self.global_name,
                additional_args: self.additional_args,
                structure: res.structure,
            },
            index_permutation: res.index_permutation,
        })
    }

    fn dual(self) -> Self {
        Self {
            structure: self.structure.dual(),
            global_name: self.global_name,
            additional_args: self.additional_args,
        }
    }

    delegate! {
        to self.structure{
            fn is_fully_self_dual(&self) -> bool;
            fn external_reps_iter(&self) -> impl Iterator<Item = Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn external_indices_iter(&self) -> impl Iterator<Item = Aind>;
            fn external_dims_iter(&self)->impl Iterator<Item=Dimension>;
            fn external_structure_iter(&self) -> impl Iterator<Item = Self::Slot>;
            fn order(&self) -> usize;
            fn get_slot(&self, i: impl Into<super::SlotIndex>) -> Option<Self::Slot>;
            fn get_rep(&self, i: impl Into<super::SlotIndex>) -> Option<Representation<<Self::Slot as IsAbstractSlot>::R>>;
            fn get_aind(&self,i: impl Into<super::SlotIndex>)->Option<Aind>;
            fn get_dim(&self, i: impl Into<super::SlotIndex>) -> Option<Dimension>;
        }
    }
}

impl<N, A, T: RepName<Dual = T>, Aind: AbsInd> ScalarStructure
    for IndexlessNamedStructure<N, A, T, Aind>
{
    fn scalar_structure() -> Self {
        IndexlessNamedStructure {
            structure: IndexLess::scalar_structure(),
            global_name: None,
            additional_args: None,
        }
    }
}

impl<Name, Args, R: RepName, Aind> IndexlessNamedStructure<Name, Args, R, Aind> {
    #[must_use]
    pub fn from_iter<I, T>(iter: T, name: Name, args: Option<Args>) -> PermutedStructure<Self>
    where
        I: RepName,
        R: From<I>,
        T: IntoIterator<Item = Representation<I>>,
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

impl<N, A, R: RepName, Aind> HasName for IndexlessNamedStructure<N, A, R, Aind>
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

impl<N, A, R: RepName, Aind> From<IndexLess<R, Aind>> for IndexlessNamedStructure<N, A, R, Aind> {
    fn from(value: IndexLess<R, Aind>) -> Self {
        IndexlessNamedStructure {
            structure: value,
            global_name: None,
            additional_args: None,
        }
    }
}

impl<N, A, R: RepName, Aind> From<NamedStructure<N, A, R, Aind>>
    for IndexlessNamedStructure<N, A, R, Aind>
{
    fn from(value: NamedStructure<N, A, R, Aind>) -> Self {
        IndexlessNamedStructure {
            structure: value.structure.into(),
            global_name: value.global_name,
            additional_args: value.additional_args,
        }
    }
}

impl<N: std::fmt::Display, A: ArgDisplay, R: RepName, Aind> std::fmt::Display
    for IndexlessNamedStructure<N, A, R, Aind>
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table = Builder::new();

        table.push_record(&[
            self.global_name
                .as_ref()
                .map(|a| format!("{}", a))
                .unwrap_or("NO NAME".to_string()),
            self.additional_args
                .as_ref()
                .map(|a| a.arg_display())
                .unwrap_or("".to_string()),
        ]);
        for item in self.structure.structure.iter() {
            table.push_record(&[item.rep.to_string(), item.dim.to_string()]);
        }
        writeln!(f)?;
        table.build().with(Style::rounded()).fmt(f)
    }
}
impl<N: std::fmt::Debug, A: ArgDisplay, R: RepName, Aind> std::fmt::Debug
    for IndexlessNamedStructure<N, A, R, Aind>
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
                format!("{:?}", item.rep),
                format!("{:?}", item.dim),
            ]);
        }
        writeln!(f)?;
        write!(f, "{}", table.build().with(Style::rounded()))
    }
}
