use symbolica::{
    representations::{AsAtomView, Atom, Identifier},
    state::{State, Workspace},
};

use super::{Expr, HasTensorStructure, TensorStructure, VecSlotExtension};

#[derive(Debug)]
pub struct SymbolicTensor {
    structure: TensorStructure,
    expression: Atom,
}
#[derive(Debug)]
pub struct SymbolicTensorBuilder<'a> {
    structure: TensorStructure,
    expression: Expr<'a>,
}

impl HasTensorStructure for SymbolicTensor {
    fn structure(&self) -> &TensorStructure {
        &self.structure
    }
}

impl<'a> HasTensorStructure for SymbolicTensorBuilder<'a> {
    fn structure(&self) -> &TensorStructure {
        &self.structure
    }
}

impl<'a> SymbolicTensorBuilder<'a> {
    pub fn new(
        structure: TensorStructure,
        label: Identifier,
        ws: &'a Workspace,
        state: &'a mut State,
    ) -> Self {
        SymbolicTensorBuilder {
            expression: structure
                .clone()
                .to_symbolic(label, ws, state)
                .builder(state, ws),
            structure,
        }
    }

    pub fn from_symbolic_tensor(
        tensor: SymbolicTensor,
        state: &'a State,
        ws: &'a Workspace,
    ) -> Self {
        SymbolicTensorBuilder {
            expression: tensor.expression.builder(state, ws),
            structure: tensor.structure,
        }
    }

    pub fn finish(self) -> SymbolicTensor {
        SymbolicTensor {
            expression: self.expression.into_atom(),
            structure: self.structure,
        }
    }

    // pub fn get_builder<'b: 'a>(mut self) -> &'b mut Expr<'a> {
    //     self.expression
    // }

    pub fn contract(mut self, other: &SymbolicTensor) -> SymbolicTensorBuilder<'a> {
        if let Some((i, j)) = self.match_index(other) {
            self.structure = self.structure().merge_at(other.structure(), (i, j));
        } else {
            self.structure.append(&mut other.structure().clone());
        }
        self.expression = self.expression * other.get_atom().as_view();
        self
    }
}

impl SymbolicTensor {
    pub fn new(
        structure: TensorStructure,
        label: Identifier,
        ws: &Workspace,
        state: &mut State,
    ) -> Self {
        SymbolicTensor {
            expression: structure.clone().to_symbolic(label, ws, state),
            structure,
        }
    }

    pub fn get_atom(&self) -> &Atom {
        &self.expression
    }

    pub fn builder<'a>(self, state: &'a State, ws: &'a Workspace) -> SymbolicTensorBuilder<'a> {
        SymbolicTensorBuilder::from_symbolic_tensor(self, state, ws)
    }
}
