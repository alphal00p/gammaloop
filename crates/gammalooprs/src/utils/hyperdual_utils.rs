use std::fmt::{Display, LowerExp};
use std::ops::AddAssign;

use crate::utils::{F, FloatLike, PrecisionUpgradable};
use spenso::algebra::{algebraic_traits::RefZero, complex::Complex};
use symbolica::domains::dual::{DualNumberStructure, HyperDual};
use symbolica::domains::float::FloatLike as SymbolicaFloatLike;

pub(crate) fn new_constant<T: Clone + RefZero>(shape: &HyperDual<T>, value: &T) -> HyperDual<T> {
    let mut new = shape.clone();
    new.values[0] = value.clone();
    let new_values_iter = new.values.iter_mut().skip(1);
    for v in new_values_iter {
        *v = value.ref_zero();
    }
    new
}

pub(crate) fn new_from_values<T: Clone>(shape: &HyperDual<T>, values: &[T]) -> HyperDual<T> {
    let mut new = shape.clone();
    new.values = values.to_vec();
    new
}

pub(crate) fn shape_for_t_derivatives(num_derivatives: usize) -> Vec<Vec<usize>> {
    (0..=num_derivatives).map(|order| vec![order]).collect()
}

impl<T> PrecisionUpgradable for HyperDual<T>
where
    T: PrecisionUpgradable,
    T::Higher: Default,
    T::Lower: Default,
{
    type Higher = HyperDual<T::Higher>;
    type Lower = HyperDual<T::Lower>;

    fn higher(&self) -> Self::Higher {
        let shape = self.get_shape();
        let owned_shape = shape
            .into_iter()
            .map(|slice| slice.to_vec())
            .collect::<Vec<_>>();

        let mut new_hyperdual = HyperDual::<T::Higher>::new(owned_shape);

        new_hyperdual
            .values
            .iter_mut()
            .zip(self.values.iter())
            .for_each(|(new_v, old_v)| {
                *new_v = old_v.higher();
            });

        new_hyperdual
    }

    fn lower(&self) -> Self::Lower {
        let shape = self.get_shape();
        let owned_shape = shape
            .into_iter()
            .map(|slice| slice.to_vec())
            .collect::<Vec<_>>();

        let mut new_hyperdual = HyperDual::<T::Lower>::new(owned_shape);

        new_hyperdual
            .values
            .iter_mut()
            .zip(self.values.iter())
            .for_each(|(new_v, old_v)| {
                *new_v = old_v.lower();
            });

        new_hyperdual
    }
}

pub enum DualOrNot<T> {
    Dual(HyperDual<T>),
    NonDual(T),
}

impl<T: LowerExp> Display for DualOrNot<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DualOrNot::Dual(dual) => write!(
                f,
                "[{}]",
                dual.values
                    .iter()
                    .map(|v| format!("{:+16e}", v))
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
            DualOrNot::NonDual(value) => write!(f, "{:+16e}", value),
        }
    }
}

impl<T: symbolica::domains::float::FloatLike> AddAssign for DualOrNot<T> {
    fn add_assign(&mut self, rhs: Self) {
        match (self, rhs) {
            (DualOrNot::Dual(x), DualOrNot::Dual(y)) => {
                *x += y;
            }
            (DualOrNot::NonDual(x), DualOrNot::NonDual(y)) => {
                *x += y;
            }
            _ => panic!("Cannot add DualOrNot of different types"),
        }
    }
}

impl<T: Clone> DualOrNot<T> {
    pub(crate) fn new_from_slice(shape: &Option<HyperDual<T>>, values: &[T]) -> Self {
        match shape {
            Some(dual_shape) => {
                let hyperdual = new_from_values(dual_shape, values);
                DualOrNot::Dual(hyperdual)
            }
            None => {
                assert_eq!(values.len(), 1);
                DualOrNot::NonDual(values[0].clone())
            }
        }
    }
}

// this function assumes that the HyperDual has the correct shape for t-derivatives
pub(crate) fn extract_t_derivatives_complex<T: FloatLike>(
    dual: HyperDual<Complex<F<T>>>,
) -> Vec<Complex<F<T>>> {
    let mut result = Vec::with_capacity(dual.values.len());
    let mut n_factorial = Complex::new_re(dual.values[0].re.one());

    for (order, value) in dual.values.iter().enumerate() {
        if order > 0 {
            n_factorial = &n_factorial * Complex::new_re(value.re.from_usize(order));
        }
        result.push(value.clone() * &n_factorial);
    }

    result
}

// this function assumes that the HyperDual has the correct shape for t-derivatives
pub(crate) fn extract_t_derivatives<T: FloatLike>(dual: HyperDual<F<T>>) -> Vec<F<T>> {
    let mut result = Vec::with_capacity(dual.values.len());
    let mut n_factorial = dual.values[0].one();

    for (order, value) in dual.values.iter().enumerate() {
        if order > 0 {
            n_factorial = &n_factorial * value.from_usize(order);
        }
        result.push(value.clone() * &n_factorial);
    }

    result
}
