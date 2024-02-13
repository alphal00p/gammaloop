use super::{
    Expr, HasName, HistoryStructure, IntoId, Shadowable, Slot, StructureContract, SymbolicContract,
    TensorStructure,
};
use smartstring::alias::String;
use symbolica::{
    representations::{default::Linear, AsAtomView, Atom, Identifier},
    state::{ResettableBuffer, State, Workspace},
};

/// A fully symbolic tensor, with no concrete values.
///
/// This tensor is used to represent the structure of a tensor, and is used to perform symbolic contraction.
/// Currently contraction is just a multiplication of the atoms, but in the future this will ensure that internal indices are independent accross the contraction.
///
/// Additionally, this can also be used as a tensor structure, that tracks the history, much like [`HistoryStructure`].
#[derive(Debug)]
pub struct SymbolicTensor {
    structure: Vec<Slot>,
    expression: symbolica::representations::Atom,
}

impl TensorStructure for SymbolicTensor {
    type Structure = Vec<Slot>;

    fn structure(&self) -> &Self::Structure {
        &self.structure
    }

    fn mut_structure(&mut self) -> &mut Self::Structure {
        &mut self.structure
    }
    fn external_structure(&self) -> &[Slot] {
        &self.structure
    }
}

impl SymbolicTensor {
    pub fn from_named<N>(structure: &N, state: &mut State, ws: &Workspace) -> Option<Self>
    where
        N: TensorStructure + HasName,
        N::Name: IntoId + Clone,
    {
        Some(SymbolicTensor {
            expression: structure.to_symbolic(state, ws)?,
            structure: structure.external_structure().to_vec(),
        })
    }

    pub fn get_atom(&self) -> &Atom {
        &self.expression
    }
}
/// Symbolic contraction of two symbolic tensors is just a multiplication of the atoms.
///
/// This is more general than the other implementations, e.g. for Densetensor, as it also does external products (i.e. no contraction).
impl SymbolicContract<SymbolicTensor> for SymbolicTensor {
    type LCM = SymbolicTensor;
    fn contract_sym(
        &self,
        other: &SymbolicTensor,
        state: &State,
        ws: &Workspace,
    ) -> Option<Self::LCM> {
        // if self
        //     .external_structure()
        //     .match_indices(&other.external_structure())
        //     .is_some()
        // {
        let mut out: Atom<Linear> = Atom::new();
        other.expression.mul(state, ws, &self.expression, &mut out);

        let mut new_structure = self.structure.clone();
        new_structure.merge(&other.structure);
        Some(SymbolicTensor {
            expression: out,
            structure: new_structure,
        })
        // }
        // None
    }
}
