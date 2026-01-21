use spenso::algebra::{algebraic_traits::RefZero, complex::Complex};
use symbolica::domains::dual::HyperDual;

pub(crate) fn new_constant<T: Clone + RefZero>(shape: &HyperDual<T>, value: &T) -> HyperDual<T> {
    let mut new = shape.clone();
    new.values[0] = value.clone();
    let new_values_iter = new.values.iter_mut().skip(1);
    for v in new_values_iter {
        *v = value.ref_zero();
    }
    new
}
