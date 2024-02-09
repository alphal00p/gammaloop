use super::{
    Expr, HasTensorSkeleton, HasTensorStructure, MutTensorStructure, TensorSkeleton,
    TensorStructure,
};
use smartstring::alias::String;
use symbolica::{
    representations::{AsAtomView, Atom, Identifier},
    state::{State, Workspace},
};

#[derive(Debug)]
pub struct SymbolicTensor {
    pub structure: TensorSkeleton<String>,
    expression: Atom,
}
#[derive(Debug)]
pub struct SymbolicTensorBuilder<'a> {
    pub structure: TensorSkeleton<String>,
    expression: Expr<'a>,
}

impl HasTensorSkeleton for SymbolicTensor {
    type Name = String;
    fn skeleton(&self) -> &TensorSkeleton<Self::Name> {
        &self.structure
    }
    fn mut_skeleton(&mut self) -> &mut TensorSkeleton<Self::Name> {
        &mut self.structure
    }
}

impl<'a> HasTensorSkeleton for SymbolicTensorBuilder<'a> {
    type Name = String;
    fn skeleton(&self) -> &TensorSkeleton<Self::Name> {
        &self.structure
    }
    fn mut_skeleton(&mut self) -> &mut TensorSkeleton<Self::Name> {
        &mut self.structure
    }
}

impl<'a> SymbolicTensorBuilder<'a> {
    pub fn new(
        structure: TensorSkeleton<String>,
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
        if let Some((_, i, j)) = self.structure.match_index(&other.structure) {
            self.structure = self.structure.merge_at(&other.structure, (i, j));
        } else {
            self.structure.merge(&other.structure.clone());
        }
        self.expression = self.expression * other.get_atom().as_view();
        self
    }

    #[allow(dead_code)]
    fn internal_contract(&mut self) {}
}

impl SymbolicTensor {
    pub fn new(
        structure: TensorSkeleton<String>,
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
