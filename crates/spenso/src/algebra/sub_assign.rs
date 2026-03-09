use std::ops::SubAssign;

use crate::{structure::TensorStructure, tensors::data::DenseTensor};

impl<T, U, I> SubAssign<DenseTensor<T, I>> for DenseTensor<U, I>
where
    U: for<'a> SubAssign<&'a T>,
    I: TensorStructure + Clone,
{
    fn sub_assign(&mut self, rhs: DenseTensor<T, I>) {
        for (u, t) in self.data.iter_mut().zip(rhs.data.iter()) {
            *u -= t;
        }
    }
}
