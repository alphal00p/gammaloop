use std::fmt::Display;

use bincode_trait_derive::{Decode, Encode};
use linnet::permutation::Permutation;
use serde::{Deserialize, Serialize};

use super::{StructureError, TensorStructure, slot::IsAbstractSlot};

#[derive(Debug, Clone, PartialEq, Eq, Hash, Deserialize, Serialize, Encode, Decode)]
pub struct PermutedStructure<S> {
    pub structure: S,
    pub rep_permutation: Permutation,
    pub index_permutation: Permutation,
}

impl<S: Display> Display for PermutedStructure<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}{}",
            self.rep_permutation, self.index_permutation, self.structure
        )
    }
}

impl<S> PermutedStructure<S> {
    pub fn identity(structure: S) -> Self
    where
        S: TensorStructure,
    {
        let size = structure.order();
        PermutedStructure {
            structure,
            rep_permutation: Permutation::id(size),
            index_permutation: Permutation::id(size),
        }
    }
    pub fn map_structure<U>(self, f: impl FnOnce(S) -> U) -> PermutedStructure<U> {
        PermutedStructure {
            structure: f(self.structure),
            rep_permutation: self.rep_permutation,
            index_permutation: self.index_permutation,
        }
    }

    pub fn reindex<I: IntoIterator<Item: Into<<S::Slot as IsAbstractSlot>::Aind>>>(
        self,
        indices: I,
    ) -> Result<PermutedStructure<S::Indexed>, StructureError>
    where
        S: TensorStructure,
    {
        let mut indices = indices.into_iter().map(|i| i.into()).collect::<Vec<_>>();
        if self.structure.order() != indices.len() {
            return Err(StructureError::WrongNumberOfArguments(
                self.structure.order(),
                indices.len(),
            ));
        }

        self.rep_permutation.apply_slice_in_place(&mut indices);
        let mut structure = self.structure.reindex(&indices)?;
        // println!("Rep:{}", self.rep_permutation);
        // println!("Ind:{}", structure.index_permutation);
        structure.rep_permutation = self.rep_permutation;
        Ok(structure)
    }
}
pub trait Perm: Sized {
    type Permuted;
    type Wrapped<P>;
    fn permute_inds(self) -> Self::Permuted;
    fn permute_reps(self) -> Self::Permuted;
    fn permute_inds_wrapped(self) -> Self::Wrapped<Self::Permuted>;
    fn permute_reps_wrapped(self) -> Self::Wrapped<Self::Permuted>;
}

impl<S> Perm for PermutedStructure<S>
where
    S: PermuteTensor,
{
    type Wrapped<P> = PermutedStructure<P>;
    type Permuted = S::Permuted;
    fn permute_inds(self) -> Self::Permuted {
        self.structure.permute_inds(&self.index_permutation)
    }
    fn permute_inds_wrapped(self) -> Self::Wrapped<Self::Permuted> {
        PermutedStructure {
            structure: self.structure.permute_inds(&self.index_permutation),
            index_permutation: self.index_permutation,
            rep_permutation: self.rep_permutation,
        }
    }

    fn permute_reps_wrapped(self) -> Self::Wrapped<Self::Permuted> {
        PermutedStructure {
            structure: self.structure.permute_reps(&self.rep_permutation),
            index_permutation: self.index_permutation,
            rep_permutation: self.rep_permutation,
        }
    }

    fn permute_reps(self) -> Self::Permuted {
        self.structure.permute_reps(&self.rep_permutation)
    }
}

pub trait PermuteTensor {
    type Permuted;
    type Id;
    type IdSlot;

    fn permute_inds(self, permutation: &Permutation) -> Self::Permuted;

    fn permute_reps(self, rep_perm: &Permutation) -> Self::Permuted;
    fn id(i: Self::IdSlot, j: Self::IdSlot) -> Self::Id;
}
