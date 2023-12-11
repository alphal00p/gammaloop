use super::{HasTensorStructure, TensorStructure};

#[allow(dead_code)]
struct SymbolicTensor {
    structure: TensorStructure,
    expression: String,
}

impl HasTensorStructure for SymbolicTensor {
    fn structure(&self) -> &TensorStructure {
        &self.structure
    }
}
