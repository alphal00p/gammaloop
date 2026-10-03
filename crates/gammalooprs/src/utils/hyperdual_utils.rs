use std::fmt::{Display, LowerExp};
use std::ops::AddAssign;

use crate::cff::CutCFFIndex;
use crate::utils::{F, FloatLike, PrecisionUpgradable};
use itertools::Itertools;
use spenso::algebra::{algebraic_traits::RefZero, complex::Complex};
use symbolica::domains::dual::{DualNumberStructure, HyperDual};

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

pub(crate) fn simple_n_deriv_shape(num_derivatives: usize) -> Vec<Vec<usize>> {
    mixed_derivative_shape(&[num_derivatives])
}

/// Retain every mixed Taylor coefficient up to each axis's requested order.
/// Linear coefficients follow axis order; higher degrees use lexicographic order,
/// preserving the existing LU and threshold evaluator layouts. Zero-order axes
/// remain present, so callers must locate active variables by their multi-index.
pub(crate) fn mixed_derivative_shape(max_orders: &[usize]) -> Vec<Vec<usize>> {
    let mut shape = max_orders
        .iter()
        .map(|order| 0..=*order)
        .multi_cartesian_product()
        .collect_vec();
    shape.sort_by(|left, right| {
        let left_degree: usize = left.iter().sum();
        let right_degree: usize = right.iter().sum();
        left_degree.cmp(&right_degree).then_with(|| {
            if left_degree == 1 {
                right.cmp(left)
            } else {
                left.cmp(right)
            }
        })
    });
    shape
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct CutCFFVariableIndices {
    pub lu_cut: Option<usize>,
    pub left_threshold: Option<usize>,
    pub right_threshold: Option<usize>,
}

pub(crate) fn variable_indices_from_cut_cff_index(
    cut_cff_index: &CutCFFIndex,
) -> CutCFFVariableIndices {
    let mut next_index = 0;

    let lu_cut = cut_cff_index
        .lu_cut_order
        .filter(|order| *order > 1)
        .map(|_| {
            let index = next_index;
            next_index += 1;
            index
        });
    let left_threshold = cut_cff_index
        .left_threshold_order
        .filter(|order| *order > 1)
        .map(|_| {
            let index = next_index;
            next_index += 1;
            index
        });
    let right_threshold = cut_cff_index
        .right_threshold_order
        .filter(|order| *order > 1)
        .map(|_| {
            let index = next_index;
            next_index += 1;
            index
        });

    CutCFFVariableIndices {
        lu_cut,
        left_threshold,
        right_threshold,
    }
}

pub(crate) fn shape_from_cut_cff_index(cut_cff_index: &CutCFFIndex) -> Option<Vec<Vec<usize>>> {
    let max_derivative_shape = [
        cut_cff_index.lu_cut_order,
        cut_cff_index.left_threshold_order,
        cut_cff_index.right_threshold_order,
    ]
    .into_iter()
    .flatten()
    .filter(|order| *order > 1)
    .map(|order| order - 1)
    .collect_vec();

    (!max_derivative_shape.is_empty()).then(|| mixed_derivative_shape(&max_derivative_shape))
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

#[derive(Clone, Debug)]
pub enum DualOrNot<T> {
    Dual(HyperDual<T>),
    NonDual(T),
}

impl<T> DualOrNot<T> {
    pub fn unwrap_real(self) -> T {
        match self {
            DualOrNot::Dual(_dual) => panic!("Cannot unwrap real value from Dual variant"),
            DualOrNot::NonDual(value) => value,
        }
    }
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

impl<T> AddAssign for DualOrNot<T>
where
    T: AddAssign<T>,
{
    fn add_assign(&mut self, rhs: Self) {
        match (self, rhs) {
            (DualOrNot::Dual(x), DualOrNot::Dual(y)) => {
                assert_eq!(x.values.len(), y.values.len());
                for (lhs, rhs) in x.values.iter_mut().zip(y.values) {
                    *lhs += rhs;
                }
            }
            (DualOrNot::NonDual(x), DualOrNot::NonDual(y)) => {
                *x += y;
            }
            _ => panic!("Cannot add DualOrNot of different types"),
        }
    }
}

impl<T: Clone> DualOrNot<T> {
    pub fn new_from_slice(shape: &Option<HyperDual<T>>, values: &[T]) -> Self {
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

pub(crate) fn dualize_dual_t_to_dual_r_t<T: FloatLike>(
    t_dual: HyperDual<F<T>>,
    target_shape: HyperDual<F<T>>,
    variable: usize,
) -> HyperDual<F<T>> {
    let mut new_dual = new_constant(&target_shape, &t_dual.values[0]);
    let n_variables_of_target_shape = target_shape.get_shape()[0].len();

    debug_assert!(variable < n_variables_of_target_shape);

    for (i, value) in t_dual.values.iter().enumerate().skip(1) {
        let dual_shape_to_find = {
            let mut shape = vec![0; n_variables_of_target_shape];
            shape[variable] = i;
            shape
        };

        #[allow(clippy::expect_fun_call)]
        let index_of_derivative = target_shape
            .get_shape()
            .iter()
            .position(|shape| shape == &dual_shape_to_find)
            .expect(&format!(
                "Could not find derivative shape: {:?} in {:?}",
                dual_shape_to_find,
                target_shape.get_shape()
            ));

        new_dual.values[index_of_derivative] = value.clone();
    }

    new_dual
}

pub(crate) fn extract_coefficient_t_duals<T: Clone + RefZero + Default>(
    dual: &HyperDual<T>,
    t_variable: usize,
) -> (Vec<Vec<usize>>, Vec<HyperDual<T>>) {
    let mixed_shape = dual
        .get_shape()
        .into_iter()
        .map(|orders| orders.to_vec())
        .collect_vec();
    let n_variables = mixed_shape.first().map_or(0, Vec::len);

    debug_assert!(t_variable < n_variables);

    let max_t_order = mixed_shape
        .iter()
        .map(|orders| orders[t_variable])
        .max()
        .unwrap_or(0);
    let t_dual_shape = HyperDual::new(simple_n_deriv_shape(max_t_order));

    let mut coefficient_orders = Vec::<Vec<usize>>::new();
    let mut coefficient_duals = Vec::<HyperDual<T>>::new();

    for (orders, value) in mixed_shape.iter().zip(dual.values.iter()) {
        let t_order = orders[t_variable];
        let projected_orders = orders
            .iter()
            .enumerate()
            .filter_map(|(index, order)| (index != t_variable).then_some(*order))
            .collect_vec();

        let coefficient_index = if let Some(existing_index) = coefficient_orders
            .iter()
            .position(|existing_orders| existing_orders == &projected_orders)
        {
            existing_index
        } else {
            coefficient_orders.push(projected_orders);
            coefficient_duals.push(new_constant(&t_dual_shape, &value.ref_zero()));
            coefficient_duals.len() - 1
        };

        coefficient_duals[coefficient_index].values[t_order] = value.clone();
    }

    (coefficient_orders, coefficient_duals)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mixed_derivative_shape_supports_constant_axes() {
        assert_eq!(mixed_derivative_shape(&[]), vec![Vec::<usize>::new()]);
        assert_eq!(mixed_derivative_shape(&[0, 0]), vec![vec![0, 0]]);
        assert_eq!(simple_n_deriv_shape(0), vec![vec![0]]);
        assert_eq!(
            mixed_derivative_shape(&[0, 2, 0, 1]),
            vec![
                vec![0, 0, 0, 0],
                vec![0, 1, 0, 0],
                vec![0, 0, 0, 1],
                vec![0, 1, 0, 1],
                vec![0, 2, 0, 0],
                vec![0, 2, 0, 1],
            ]
        );

        for orders in [vec![], vec![0, 0], vec![0, 2, 0, 1]] {
            let shape = mixed_derivative_shape(&orders);
            let mut values = vec![F(0.0_f64); shape.len()];
            values[0] = F(3.0_f64);
            let constant = HyperDual::from_values(shape, values);
            let square = constant.clone() * constant;
            assert_eq!(square.values[0], F(9.0_f64));
            assert!(square.values[1..].iter().all(|value| *value == F(0.0_f64)));
        }
    }

    #[test]
    fn mixed_derivative_shape_preserves_four_axis_polynomial_coefficients() {
        let shape = mixed_derivative_shape(&[2, 1, 1, 2]);
        let variables = [2.0_f64, 3.0_f64, 5.0_f64, 7.0_f64]
            .into_iter()
            .enumerate()
            .map(|(axis, base)| {
                let mut values = vec![F(0.0_f64); shape.len()];
                values[0] = F(base);
                let mut linear_order = vec![0; 4];
                linear_order[axis] = 1;
                let linear_index = shape
                    .iter()
                    .position(|order| *order == linear_order)
                    .unwrap();
                values[linear_index] = F(1.0_f64);
                HyperDual::from_values(shape.clone(), values)
            })
            .collect_vec();
        let polynomial = variables[0].clone()
            * &variables[0]
            * &variables[1]
            * &variables[2]
            * &variables[3]
            * &variables[3];

        // Coefficients, rather than ordinary derivatives: the highest mixed
        // derivative is 2! * 1! * 1! * 2! times its stored coefficient.
        for (orders, expected) in [
            ([0, 0, 0, 0], 2940.0_f64),
            ([1, 0, 1, 1], 168.0_f64),
            ([2, 1, 1, 2], 1.0_f64),
        ] {
            let index = shape
                .iter()
                .position(|order| order.as_slice() == orders)
                .unwrap();
            assert_eq!(polynomial.values[index], F(expected));
        }
    }

    #[test]
    fn dualize_dual_t_to_dual_r_t_uses_requested_target_variable() {
        let t_dual = HyperDual::new(simple_n_deriv_shape(1));
        let t_dual = new_from_values(&t_dual, &[F(3.0_f64), F(5.0_f64)]);

        let target_shape = HyperDual::new(
            shape_from_cut_cff_index(&CutCFFIndex {
                left_threshold_order: Some(2),
                right_threshold_order: None,
                lu_cut_order: Some(2),
            })
            .unwrap(),
        );

        let dualized = dualize_dual_t_to_dual_r_t(t_dual, target_shape, 1);

        assert_eq!(
            dualized.values,
            vec![F(3.0_f64), F(0.0_f64), F(5.0_f64), F(0.0_f64)]
        );
    }

    #[test]
    fn dualize_dual_t_to_dual_r_t_keeps_constant_inputs_constant() {
        let t_dual = HyperDual::new(simple_n_deriv_shape(0));
        let t_dual = new_from_values(&t_dual, &[F(7.0_f64)]);

        let target_shape = HyperDual::new(
            shape_from_cut_cff_index(&CutCFFIndex {
                left_threshold_order: Some(2),
                right_threshold_order: None,
                lu_cut_order: Some(2),
            })
            .unwrap(),
        );

        let dualized = dualize_dual_t_to_dual_r_t(t_dual, target_shape, 1);

        assert_eq!(
            dualized.values,
            vec![F(7.0_f64), F(0.0_f64), F(0.0_f64), F(0.0_f64)]
        );
    }

    #[test]
    fn extract_coefficient_t_duals_factorizes_single_threshold_mixed_dual() {
        let mixed_shape = HyperDual::new(
            shape_from_cut_cff_index(&CutCFFIndex {
                left_threshold_order: Some(2),
                right_threshold_order: None,
                lu_cut_order: Some(2),
            })
            .unwrap(),
        );
        let mixed_dual = new_from_values(
            &mixed_shape,
            &[F(2.0_f64), F(3.0_f64), F(5.0_f64), F(7.0_f64)],
        );

        let (coefficient_orders, coefficient_duals) = extract_coefficient_t_duals(&mixed_dual, 0);

        assert_eq!(coefficient_orders, vec![vec![0], vec![1]]);
        assert_eq!(coefficient_duals[0].values, vec![F(2.0_f64), F(3.0_f64)]);
        assert_eq!(coefficient_duals[1].values, vec![F(5.0_f64), F(7.0_f64)]);
    }

    #[test]
    fn extract_coefficient_t_duals_factorizes_iterated_mixed_dual() {
        let mixed_shape = HyperDual::new(
            shape_from_cut_cff_index(&CutCFFIndex {
                left_threshold_order: Some(2),
                right_threshold_order: Some(2),
                lu_cut_order: Some(2),
            })
            .unwrap(),
        );
        let mixed_dual = new_from_values(
            &mixed_shape,
            &[
                F(1.0_f64),
                F(2.0_f64),
                F(3.0_f64),
                F(4.0_f64),
                F(5.0_f64),
                F(6.0_f64),
                F(7.0_f64),
                F(8.0_f64),
            ],
        );

        let (coefficient_orders, coefficient_duals) = extract_coefficient_t_duals(&mixed_dual, 0);

        assert_eq!(
            coefficient_orders,
            vec![vec![0, 0], vec![1, 0], vec![0, 1], vec![1, 1]]
        );
        assert_eq!(coefficient_duals[0].values, vec![F(1.0_f64), F(2.0_f64)]);
        assert_eq!(coefficient_duals[1].values, vec![F(3.0_f64), F(7.0_f64)]);
        assert_eq!(coefficient_duals[2].values, vec![F(4.0_f64), F(6.0_f64)]);
        assert_eq!(coefficient_duals[3].values, vec![F(5.0_f64), F(8.0_f64)]);
    }

    #[test]
    fn shape_from_cut_cff_index_preserves_canonical_t_left_right_order() {
        let three_axis_shape = shape_from_cut_cff_index(&CutCFFIndex {
            left_threshold_order: Some(3),
            right_threshold_order: Some(2),
            lu_cut_order: Some(2),
        })
        .unwrap();

        assert_eq!(
            three_axis_shape,
            vec![
                vec![0, 0, 0],
                vec![1, 0, 0],
                vec![0, 1, 0],
                vec![0, 0, 1],
                vec![0, 1, 1],
                vec![0, 2, 0],
                vec![1, 0, 1],
                vec![1, 1, 0],
                vec![0, 2, 1],
                vec![1, 1, 1],
                vec![1, 2, 0],
                vec![1, 2, 1],
            ]
        );

        assert_eq!(
            shape_from_cut_cff_index(&CutCFFIndex {
                left_threshold_order: Some(3),
                right_threshold_order: None,
                lu_cut_order: Some(3),
            }),
            Some(vec![
                vec![0, 0],
                vec![1, 0],
                vec![0, 1],
                vec![0, 2],
                vec![1, 1],
                vec![2, 0],
                vec![1, 2],
                vec![2, 1],
                vec![2, 2],
            ])
        );

        assert_eq!(
            shape_from_cut_cff_index(&CutCFFIndex {
                left_threshold_order: Some(1),
                right_threshold_order: Some(2),
                lu_cut_order: Some(2),
            }),
            Some(vec![vec![0, 0], vec![1, 0], vec![0, 1], vec![1, 1]])
        );

        assert_eq!(
            shape_from_cut_cff_index(&CutCFFIndex {
                left_threshold_order: Some(2),
                right_threshold_order: Some(1),
                lu_cut_order: Some(1),
            }),
            Some(vec![vec![0], vec![1]])
        );
    }

    #[test]
    fn variable_indices_from_cut_cff_index_tracks_active_axes_without_reordering() {
        assert_eq!(
            variable_indices_from_cut_cff_index(&CutCFFIndex {
                left_threshold_order: Some(2),
                right_threshold_order: Some(2),
                lu_cut_order: Some(2),
            }),
            CutCFFVariableIndices {
                lu_cut: Some(0),
                left_threshold: Some(1),
                right_threshold: Some(2),
            }
        );

        assert_eq!(
            variable_indices_from_cut_cff_index(&CutCFFIndex {
                left_threshold_order: Some(2),
                right_threshold_order: Some(2),
                lu_cut_order: Some(1),
            }),
            CutCFFVariableIndices {
                lu_cut: None,
                left_threshold: Some(0),
                right_threshold: Some(1),
            }
        );

        assert_eq!(
            variable_indices_from_cut_cff_index(&CutCFFIndex {
                left_threshold_order: Some(1),
                right_threshold_order: Some(2),
                lu_cut_order: Some(2),
            }),
            CutCFFVariableIndices {
                lu_cut: Some(0),
                left_threshold: None,
                right_threshold: Some(1),
            }
        );
    }
}
