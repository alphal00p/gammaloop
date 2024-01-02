use symbolica::{
    representations::{AsAtomView, Atom},
    state::{State, Workspace},
};

use super::{ContractableWithDense, Expr, HasTensorStructure, TensorStructure, VecSlotExtension};

#[allow(dead_code)]
struct SymbolicTensor {
    structure: TensorStructure,
    expression: Atom,
}

struct SymbolicTensorBuilder<'a> {
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
    pub fn new(structure: TensorStructure, expression: Expr<'a>) -> Self {
        SymbolicTensorBuilder {
            structure,
            expression,
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
    pub fn new(structure: TensorStructure, expression: Atom) -> Self {
        SymbolicTensor {
            structure,
            expression,
        }
    }

    pub fn get_atom(&self) -> &Atom {
        &self.expression
    }

    pub fn builder<'a>(self, state: &'a State, ws: &'a Workspace) -> SymbolicTensorBuilder<'a> {
        SymbolicTensorBuilder::new(self.structure, self.expression.builder(state, ws))
    }
}
