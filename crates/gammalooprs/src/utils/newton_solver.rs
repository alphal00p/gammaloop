use serde::Serialize;
use std::fmt::{Display, Formatter};
use symbolica::domains::dual::HyperDual;

use crate::utils::{F, FloatLike};

#[derive(Clone, Debug)]
pub(crate) enum SafeguardedNewtonError<T: FloatLike> {
    InvalidInside {
        radius: F<T>,
        value: F<T>,
    },
    InvalidOutside {
        radius: F<T>,
        value: F<T>,
        bracket_expansions: usize,
    },
    InvalidDerivative {
        result: NewtonIterationResult<T>,
    },
    DidNotConverge {
        result: NewtonIterationResult<T>,
        lower_bound: F<T>,
        upper_bound: F<T>,
    },
}

impl<T: FloatLike> Display for SafeguardedNewtonError<T> {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidInside { radius, value } => write!(
                formatter,
                "inside radius {radius} has non-finite or non-negative value {value}"
            ),
            Self::InvalidOutside {
                radius,
                value,
                bracket_expansions,
            } => write!(
                formatter,
                "outside radius {radius} has invalid value {value} after {bracket_expansions} bracket expansions"
            ),
            Self::InvalidDerivative { result } => write!(
                formatter,
                "root candidate {} has invalid derivative {} with residual {} after {} iterations",
                result.solution,
                result.derivative_at_solution,
                result.error_of_function,
                result.num_iterations_used,
            ),
            Self::DidNotConverge {
                result,
                lower_bound,
                upper_bound,
            } => write!(
                formatter,
                "did not converge in bracket [{lower_bound}, {upper_bound}]: candidate {}, derivative {}, residual {} after {} iterations",
                result.solution,
                result.derivative_at_solution,
                result.error_of_function,
                result.num_iterations_used,
            ),
        }
    }
}

fn is_finite<T: FloatLike>(value: &F<T>) -> bool {
    !value.is_nan() && !value.is_infinite()
}

/// Find a radial root while preserving a finite negative/positive bracket.
///
/// Newton steps are used only when their derivative and proposed radius are
/// finite and the proposed radius stays strictly inside the bracket. Otherwise
/// the next iterate is the bracket midpoint.
pub(crate) fn safeguarded_newton_iteration_and_derivative<T: FloatLike>(
    inside_radius: &F<T>,
    outside_radius_guess: &F<T>,
    f_x_and_df_x: impl Fn(&F<T>) -> (F<T>, F<T>),
    tolerance: &F<T>,
    max_iterations: usize,
    max_bracket_expansions: usize,
    e_cm: &F<T>,
) -> Result<NewtonIterationResult<T>, SafeguardedNewtonError<T>> {
    let zero = inside_radius.zero();
    let two = inside_radius.from_i64(2);
    let maximum_residual = inside_radius.epsilon() * tolerance * e_cm;
    let (inside_value, _) = f_x_and_df_x(inside_radius);
    if !is_finite(&inside_value) || inside_value >= -maximum_residual.clone() {
        return Err(SafeguardedNewtonError::InvalidInside {
            radius: inside_radius.clone(),
            value: inside_value,
        });
    }

    let mut lower_bound = inside_radius.clone();
    let mut upper_bound = outside_radius_guess.clone();
    if !is_finite(&upper_bound) || upper_bound <= lower_bound {
        let upper_value = if is_finite(&upper_bound) {
            f_x_and_df_x(&upper_bound).0
        } else {
            upper_bound.clone()
        };
        return Err(SafeguardedNewtonError::InvalidOutside {
            radius: upper_bound,
            value: upper_value,
            bracket_expansions: 0,
        });
    }

    let (mut upper_value, _) = f_x_and_df_x(&upper_bound);
    let mut bracket_expansions = 0;
    while is_finite(&upper_value)
        && upper_value <= zero
        && bracket_expansions < max_bracket_expansions
    {
        upper_bound *= &two;
        (upper_value, _) = f_x_and_df_x(&upper_bound);
        bracket_expansions += 1;
    }

    if !is_finite(&upper_bound) || !is_finite(&upper_value) || upper_value <= zero {
        return Err(SafeguardedNewtonError::InvalidOutside {
            radius: upper_bound,
            value: upper_value,
            bracket_expansions,
        });
    }

    let mut solution = upper_bound.clone();
    let (mut value, mut derivative) = f_x_and_df_x(&solution);
    for iteration in 0..=max_iterations {
        let current_result = NewtonIterationResult {
            solution: solution.clone(),
            derivative_at_solution: derivative.clone(),
            error_of_function: value.clone(),
            num_iterations_used: iteration,
        };

        if is_finite(&value) && value.abs() <= maximum_residual {
            if is_finite(&derivative) && derivative > zero {
                return Ok(current_result);
            }
            return Err(SafeguardedNewtonError::InvalidDerivative {
                result: current_result,
            });
        }

        if iteration == max_iterations {
            return Err(SafeguardedNewtonError::DidNotConverge {
                result: current_result,
                lower_bound,
                upper_bound,
            });
        }

        if is_finite(&value) {
            if value < zero {
                lower_bound = solution.clone();
            } else {
                upper_bound = solution.clone();
            }
        }

        let midpoint = (&lower_bound + &upper_bound) / &two;
        let newton_candidate = if is_finite(&value) && is_finite(&derivative) && derivative > zero {
            Some(&solution - &value / &derivative)
        } else {
            None
        };

        solution = match newton_candidate {
            Some(candidate)
                if is_finite(&candidate) && candidate > lower_bound && candidate < upper_bound =>
            {
                candidate
            }
            _ => midpoint,
        };
        (value, derivative) = f_x_and_df_x(&solution);
    }

    unreachable!("safeguarded Newton loop always returns")
}

/// root finding, returns the derivative at the root, so that we don't have to recompute it.
/// Also returns the value of the function whose root is being found and the number of iterations used for debug information
pub(crate) fn newton_iteration_and_derivative<T: FloatLike>(
    guess: &F<T>,
    f_x_and_df_x: impl Fn(&F<T>) -> (F<T>, F<T>),
    tolerance: &F<T>,
    max_iterations: usize,
    e_cm: &F<T>,
) -> NewtonIterationResult<T> {
    let mut x = guess.clone();
    let (mut val_f_x, mut val_df_x) = f_x_and_df_x(&x);

    let mut iteration = 0;

    while iteration < max_iterations && val_f_x.abs() > guess.epsilon() * tolerance * e_cm {
        x -= val_f_x / val_df_x;
        (val_f_x, val_df_x) = f_x_and_df_x(&x);
        iteration += 1;
    }

    NewtonIterationResult {
        solution: x,
        derivative_at_solution: val_df_x,
        error_of_function: val_f_x,
        num_iterations_used: iteration,
    }
}

#[allow(dead_code)]
pub(crate) fn newton_iteration_and_derivative_dual<T: FloatLike>(
    guess: &HyperDual<F<T>>,
    f_x_and_df_x: impl Fn(&HyperDual<F<T>>) -> (HyperDual<F<T>>, HyperDual<F<T>>),
    tolerance: &F<T>,
    max_iterations: usize,
    e_cm: &F<T>,
) -> NewtonIterationResultDual<T> {
    let mut x = guess.clone();
    let (mut val_f_x, mut val_df_x) = f_x_and_df_x(&x);

    let mut iteration = 0;

    while iteration < max_iterations
        && val_f_x.values[0].abs() > guess.values[0].epsilon() * tolerance * e_cm
    {
        x -= val_f_x / val_df_x;
        (val_f_x, val_df_x) = f_x_and_df_x(&x);
        iteration += 1;
    }

    NewtonIterationResultDual {
        solution: x,
        derivative_at_solution: val_df_x,
        error_of_function: val_f_x.values[0].clone(),
        num_iterations_used: iteration,
    }
}

#[derive(Serialize, Clone, Debug)]
pub(crate) struct NewtonIterationResult<T: FloatLike> {
    pub solution: F<T>,
    pub derivative_at_solution: F<T>,
    pub error_of_function: F<T>,
    pub num_iterations_used: usize,
}

#[derive(Clone, Debug)]
#[allow(dead_code)]
pub(crate) struct NewtonIterationResultDual<T: FloatLike> {
    pub solution: HyperDual<F<T>>,
    pub derivative_at_solution: HyperDual<F<T>>,
    pub error_of_function: F<T>,
    pub num_iterations_used: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn safeguarded_newton_expands_bracket_and_converges() {
        let result = safeguarded_newton_iteration_and_derivative(
            &F(0.0),
            &F(0.5),
            |x| (x * x - F(2.0), F(2.0) * x),
            &F(8.0),
            40,
            64,
            &F(1.0),
        )
        .unwrap();

        let expected_solution = result.solution.from_i64(2).sqrt();
        let solution_tolerance = result.solution.epsilon() * result.solution.from_i64(64);
        assert!((result.solution - expected_solution).abs() <= solution_tolerance);

        let residual_tolerance =
            result.error_of_function.epsilon() * result.error_of_function.from_i64(8);
        assert!(result.error_of_function.abs() <= residual_tolerance);
        assert!(result.derivative_at_solution > result.derivative_at_solution.zero());
    }

    #[test]
    fn safeguarded_newton_accepts_a_few_ulp_residual() {
        let residual = 2.2737367544323206e-13;
        let result = safeguarded_newton_iteration_and_derivative(
            &F(0.0),
            &F(1.0),
            |x| {
                let value = if x == &F(1.0) {
                    F(residual)
                } else {
                    x - F(1.0)
                };
                (value, F(1.0))
            },
            &F(8.0),
            40,
            64,
            &F(1000.0),
        )
        .unwrap();

        assert_eq!(result.solution, F(1.0));
        assert_eq!(result.error_of_function, F(residual));
    }

    #[test]
    fn safeguarded_newton_rejects_a_large_discontinuous_residual() {
        let error = safeguarded_newton_iteration_and_derivative(
            &F(0.0),
            &F(1.0),
            |x| {
                if x < &F(1.0) {
                    (F(-1.0), F(1.0))
                } else {
                    (F(741.0), F(1.0))
                }
            },
            &F(8.0),
            20,
            64,
            &F(1000.0),
        )
        .unwrap_err();

        assert!(matches!(
            error,
            SafeguardedNewtonError::DidNotConverge { .. }
        ));
    }

    #[test]
    fn safeguarded_newton_rejects_a_non_interior_center() {
        let error = safeguarded_newton_iteration_and_derivative(
            &F(0.0),
            &F(1.0),
            |x| (*x, F(1.0)),
            &F(8.0),
            40,
            64,
            &F(1.0),
        )
        .unwrap_err();

        assert!(matches!(
            error,
            SafeguardedNewtonError::InvalidInside { .. }
        ));
    }

    #[test]
    fn safeguarded_newton_rejects_an_invalid_outside_bracket() {
        let error = safeguarded_newton_iteration_and_derivative(
            &F(0.0),
            &F(1.0),
            |_| (F(-1.0), F(1.0)),
            &F(8.0),
            40,
            4,
            &F(1.0),
        )
        .unwrap_err();

        assert!(matches!(
            error,
            SafeguardedNewtonError::InvalidOutside {
                bracket_expansions: 4,
                ..
            }
        ));
    }

    #[test]
    fn safeguarded_newton_reports_invalid_outside_endpoint_value() {
        let error = safeguarded_newton_iteration_and_derivative(
            &F(1.0),
            &F(0.5),
            |x| (x - F(2.0), F(1.0)),
            &F(8.0),
            40,
            64,
            &F(1.0),
        )
        .unwrap_err();

        assert!(matches!(
            error,
            SafeguardedNewtonError::InvalidOutside {
                radius: F(0.5),
                value: F(-1.5),
                bracket_expansions: 0,
            }
        ));
    }

    #[test]
    fn safeguarded_newton_reports_non_finite_outside_radius() {
        let nan = F(0.0) / F(0.0);
        let error = safeguarded_newton_iteration_and_derivative(
            &F(0.0),
            &nan,
            |x| (x - F(1.0), F(1.0)),
            &F(8.0),
            40,
            64,
            &F(1.0),
        )
        .unwrap_err();

        match error {
            SafeguardedNewtonError::InvalidOutside {
                radius,
                value,
                bracket_expansions: 0,
            } => {
                assert!(radius.is_nan());
                assert!(value.is_nan());
            }
            unexpected => panic!("unexpected safeguarded Newton error: {unexpected:?}"),
        }
    }

    #[test]
    fn safeguarded_newton_rejects_an_invalid_root_derivative() {
        let error = safeguarded_newton_iteration_and_derivative(
            &F(0.0),
            &F(1.0),
            |x| (x - F(1.0), F(0.0)),
            &F(8.0),
            40,
            64,
            &F(1.0),
        )
        .unwrap_err();

        assert!(matches!(
            error,
            SafeguardedNewtonError::InvalidDerivative { .. }
        ));
    }
}
