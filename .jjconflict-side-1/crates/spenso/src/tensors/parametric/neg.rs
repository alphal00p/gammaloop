use crate::structure::TensorStructure;

use super::{ParamOrComposite, ParamOrConcrete, ParamTensor};

impl<T, U, S> std::ops::Neg for ParamOrConcrete<T, S>
where
    T: std::ops::Neg<Output = U>,
    S: TensorStructure + Clone,
{
    type Output = ParamOrConcrete<U, S>;
    fn neg(self) -> Self::Output {
        match self {
            ParamOrConcrete::Concrete(d) => ParamOrConcrete::Concrete(-d),
            ParamOrConcrete::Param(s) => ParamOrConcrete::Param(-s),
        }
    }
}

impl<S> std::ops::Neg for ParamTensor<S>
where
    S: TensorStructure + Clone,
{
    type Output = ParamTensor<S>;
    fn neg(self) -> Self::Output {
        ParamTensor {
            param_type: ParamOrComposite::Composite,
            tensor: -self.tensor,
        }
    }
}
