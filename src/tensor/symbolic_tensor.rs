use super::{Expr, HasTensorStructure, TensorStructure};

#[allow(dead_code)]
struct SymbolicTensor<'a> {
    structure: TensorStructure,
    expression: Expr<'a>,
}

impl<'a> HasTensorStructure for SymbolicTensor<'a> {
    fn structure(&self) -> &TensorStructure {
        &self.structure
    }
}
