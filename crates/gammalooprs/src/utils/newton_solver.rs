use serde::Serialize;
use std::{
    collections::BTreeMap,
    fmt::{Display, Formatter},
};
use symbolica::domains::{dual::HyperDual, float::Float as SymbolicaFloat};
use tracing::debug;

use crate::utils::{F, FloatLike};

const MINIMUM_PRECISION_RESIDUAL_IMPROVEMENT: i64 = 8;

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct RadialRootIdentity(String);

impl RadialRootIdentity {
    pub(crate) fn new(description: String) -> Self {
        Self(description)
    }
}

impl Display for RadialRootIdentity {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(formatter)
    }
}

#[derive(Clone, Debug)]
struct LocalRootConsistency {
    near_radius_step: Option<SymbolicaFloat>,
    near_lower_value: Option<SymbolicaFloat>,
    near_upper_value: Option<SymbolicaFloat>,
    near_secant_derivative_ratio: Option<SymbolicaFloat>,
    far_secant_derivative_ratio: Option<SymbolicaFloat>,
    secant_ratio_agreement: Option<SymbolicaFloat>,
    is_valid: bool,
}

impl LocalRootConsistency {
    fn invalid() -> Self {
        Self {
            near_radius_step: None,
            near_lower_value: None,
            near_upper_value: None,
            near_secant_derivative_ratio: None,
            far_secant_derivative_ratio: None,
            secant_ratio_agreement: None,
            is_valid: false,
        }
    }

    fn check<T: FloatLike>(
        result: &NewtonIterationResult<T>,
        inside_radius: &F<T>,
        f_x_and_df_x: &impl Fn(&F<T>) -> (F<T>, F<T>),
        e_cm: &F<T>,
    ) -> Self {
        // The two symmetric probes below cost four function evaluations in total. Refuse rescue
        // without probing when a caller attempted fewer iterations, so this diagnostic can never
        // evaluate the E-surface more often than the failed Newton solve it is validating.
        const PROBE_EVALUATIONS: usize = 4;
        if result.num_iterations_used < PROBE_EVALUATIONS {
            return Self::invalid();
        }

        let zero = result.solution.zero();
        let two = result.solution.from_i64(2);
        let four = result.solution.from_i64(4);
        let eight = result.solution.from_i64(8);
        let radial_scale = result
            .solution
            .abs()
            .max(e_cm.abs())
            .max(result.solution.epsilon());
        let newton_correction =
            result.error_of_function.abs() / result.derivative_at_solution.abs();
        let mut near_radius_step =
            (result.solution.epsilon().sqrt() * radial_scale).max(&four * newton_correction);
        let available_radius = &result.solution - inside_radius;
        if !is_finite(&available_radius) || available_radius <= zero {
            return Self::invalid();
        }
        let maximum_far_radius_step = &available_radius / &two;
        if &near_radius_step * &eight > maximum_far_radius_step {
            near_radius_step = &maximum_far_radius_step / &eight;
        }
        let far_radius_step = &near_radius_step * &eight;
        if !is_finite(&near_radius_step)
            || near_radius_step <= zero
            || !is_finite(&far_radius_step)
            || far_radius_step <= near_radius_step
        {
            return Self::invalid();
        }

        let probe = |radius_step: &F<T>| {
            let lower_radius = &result.solution - radius_step;
            let upper_radius = &result.solution + radius_step;
            let (lower_value, _) = f_x_and_df_x(&lower_radius);
            let (upper_value, _) = f_x_and_df_x(&upper_radius);
            let secant_derivative = (&upper_value - &lower_value) / (&two * radius_step);
            let secant_derivative_ratio = &secant_derivative / &result.derivative_at_solution;
            (
                lower_value,
                upper_value,
                secant_derivative,
                secant_derivative_ratio,
            )
        };
        let (
            near_lower_value,
            near_upper_value,
            near_secant_derivative,
            near_secant_derivative_ratio,
        ) = probe(&near_radius_step);
        let (far_lower_value, far_upper_value, far_secant_derivative, far_secant_derivative_ratio) =
            probe(&far_radius_step);
        let secant_ratio_agreement = &far_secant_derivative_ratio / &near_secant_derivative_ratio;
        let half = result.solution.from_i64(1) / &two;
        let two_again = two.clone();
        let three_quarters = result.solution.from_i64(3) / &four;
        let five_quarters = result.solution.from_i64(5) / &four;
        let is_valid = is_finite(&near_lower_value)
            && is_finite(&near_upper_value)
            && is_finite(&near_secant_derivative)
            && is_finite(&near_secant_derivative_ratio)
            && is_finite(&far_lower_value)
            && is_finite(&far_upper_value)
            && is_finite(&far_secant_derivative)
            && is_finite(&far_secant_derivative_ratio)
            && is_finite(&secant_ratio_agreement)
            && near_lower_value < zero
            && near_upper_value > zero
            && far_lower_value < zero
            && far_upper_value > zero
            && near_secant_derivative > zero
            && far_secant_derivative > zero
            && near_secant_derivative_ratio >= half
            && near_secant_derivative_ratio <= two_again
            && far_secant_derivative_ratio >= half
            && far_secant_derivative_ratio <= two_again
            && secant_ratio_agreement >= three_quarters
            && secant_ratio_agreement <= five_quarters;

        Self {
            near_radius_step: Some(near_radius_step.into()),
            near_lower_value: Some(near_lower_value.into()),
            near_upper_value: Some(near_upper_value.into()),
            near_secant_derivative_ratio: Some(near_secant_derivative_ratio.into()),
            far_secant_derivative_ratio: Some(far_secant_derivative_ratio.into()),
            secant_ratio_agreement: Some(secant_ratio_agreement.into()),
            is_valid,
        }
    }
}

#[derive(Clone, Debug)]
struct RadialRootObservation {
    precision_bits: u32,
    precision_epsilon: SymbolicaFloat,
    residual: SymbolicaFloat,
    roundoff_residual_scale: SymbolicaFloat,
    maximum_residual: SymbolicaFloat,
    solution: SymbolicaFloat,
    derivative: SymbolicaFloat,
    lower_bound: Option<SymbolicaFloat>,
    upper_bound: Option<SymbolicaFloat>,
    bracket_is_valid: bool,
    local_consistency: Option<LocalRootConsistency>,
    relative_newton_correction: SymbolicaFloat,
    relative_residual_limit: SymbolicaFloat,
}

impl RadialRootObservation {
    fn new<T: FloatLike>(
        result: &NewtonIterationResult<T>,
        bracket: Option<(&F<T>, &F<T>)>,
        local_consistency: Option<LocalRootConsistency>,
        tolerance: &F<T>,
        e_cm: &F<T>,
    ) -> Self {
        let precision_epsilon_t = result.solution.epsilon();
        let residual_t = result.error_of_function.abs();
        let derivative_abs_t = result.derivative_at_solution.abs();
        let e_cm_t = e_cm.abs();
        let tolerance_t = tolerance.abs();
        let radial_scale_t = result
            .solution
            .abs()
            .max(e_cm_t.clone())
            .max(precision_epsilon_t.clone());
        let roundoff_residual_scale_t = &precision_epsilon_t * &e_cm_t;
        let maximum_residual_t = &roundoff_residual_scale_t * &tolerance_t;
        let relative_newton_correction_t = &residual_t / &derivative_abs_t / &radial_scale_t;
        let relative_residual_limit_t = &precision_epsilon_t * &tolerance_t;

        // `SymbolicaFloat` retains the source precision. These observations cross the generic
        // stability lanes, so keeping them in a precision-carrying type avoids feeding an f128 or
        // arbitrary-precision residual through binary64 during the rescue decision.
        let precision_epsilon: SymbolicaFloat = precision_epsilon_t.into();
        let precision_bits = precision_epsilon.prec();
        let residual = residual_t.into();
        let solution = result.solution.clone().into();
        let derivative = result.derivative_at_solution.clone().into();
        let roundoff_residual_scale = roundoff_residual_scale_t.into();
        let maximum_residual = maximum_residual_t.into();
        let relative_newton_correction = relative_newton_correction_t.into();
        let relative_residual_limit = relative_residual_limit_t.into();

        let (lower_bound, upper_bound, bracket_is_valid) =
            bracket.map_or((None, None, true), |(lower_bound, upper_bound)| {
                let bracket_is_valid = lower_bound <= upper_bound
                    && &result.solution >= lower_bound
                    && &result.solution <= upper_bound;
                (
                    Some(lower_bound.clone().into()),
                    Some(upper_bound.clone().into()),
                    bracket_is_valid,
                )
            });

        Self {
            precision_bits,
            precision_epsilon,
            residual,
            roundoff_residual_scale,
            maximum_residual,
            solution,
            derivative,
            // A safeguarded bracket can collapse to one representable radius after both
            // endpoints have been updated from opposite signs. That is a localized root,
            // not a loss of the bracket invariant.
            lower_bound,
            upper_bound,
            bracket_is_valid,
            local_consistency,
            relative_newton_correction,
            relative_residual_limit,
        }
    }

    fn comparison_residual(&self) -> SymbolicaFloat {
        if self.residual >= self.roundoff_residual_scale {
            self.residual.clone()
        } else {
            self.roundoff_residual_scale.clone()
        }
    }

    fn values_are_finite(&self) -> bool {
        self.precision_epsilon.is_finite()
            && self.residual.is_finite()
            && self.roundoff_residual_scale.is_finite()
            && self.maximum_residual.is_finite()
            && self.solution.is_finite()
            && self.derivative.is_finite()
            && self.relative_newton_correction.is_finite()
            && self.relative_residual_limit.is_finite()
            && self
                .lower_bound
                .as_ref()
                .is_none_or(SymbolicaFloat::is_finite)
            && self
                .upper_bound
                .as_ref()
                .is_none_or(SymbolicaFloat::is_finite)
    }

    fn precision_rescue_improvement(&self, previous: &Self) -> Option<SymbolicaFloat> {
        let zero = SymbolicaFloat::with_val(self.precision_bits.max(previous.precision_bits), 0);
        if !self.values_are_finite()
            || !previous.values_are_finite()
            || self.precision_bits <= previous.precision_bits
            || self.solution <= zero
            || self.derivative <= zero
            || !self.bracket_is_valid
            || !self
                .local_consistency
                .as_ref()
                .is_some_and(|consistency| consistency.is_valid)
            || previous.solution <= zero
            || previous.derivative <= zero
            || !previous.bracket_is_valid
            || self.residual > previous.maximum_residual
            || self.relative_newton_correction > previous.relative_residual_limit
            || self.residual <= zero
        {
            return None;
        }

        // A lower-precision residual can be exactly zero by accidental cancellation. In that
        // case its precision-scaled roundoff floor is the meaningful comparison baseline.
        let improvement = previous.comparison_residual() / self.residual.clone();
        let required_improvement =
            SymbolicaFloat::with_val(improvement.prec(), MINIMUM_PRECISION_RESIDUAL_IMPROVEMENT);
        (improvement >= required_improvement).then_some(improvement)
    }
}

/// Transient diagnostics shared by the stability levels of one sample evaluation.
///
/// A higher-precision solve may recover a residual-limited root when it reaches the
/// original lower-precision accuracy target and clearly improves the same root. All
/// structural solver failures remain errors in every precision.
///
/// For lower and higher precisions `lo` and `hi`, rescue requires
/// `|f_hi| <= tau * epsilon_lo * E_cm`,
/// `|f_hi / f'_hi| / max(|r_hi|, E_cm) <= tau * epsilon_lo`, and
/// `max(|f_lo|, epsilon_lo * E_cm) / |f_hi| >= 8`, in addition to the bracket and
/// two-scale local-root consistency checks.
#[derive(Clone, Debug, Default)]
pub(crate) struct RadialRootDiagnostics {
    observations: BTreeMap<(RadialRootIdentity, usize), RadialRootObservation>,
    current_precision_bits: Option<u32>,
    current_occurrences: BTreeMap<RadialRootIdentity, usize>,
}

impl RadialRootDiagnostics {
    fn next_call_key<T: FloatLike>(
        &mut self,
        identity: &RadialRootIdentity,
        precision_source: &F<T>,
    ) -> (RadialRootIdentity, usize) {
        let precision: SymbolicaFloat = precision_source.clone().into();
        let precision_bits = precision.prec();
        if self.current_precision_bits != Some(precision_bits) {
            self.current_precision_bits = Some(precision_bits);
            self.current_occurrences.clear();
        }

        let occurrence = self
            .current_occurrences
            .entry(identity.clone())
            .or_default();
        let key = (identity.clone(), *occurrence);
        *occurrence += 1;
        key
    }

    fn record_observation(
        &mut self,
        key: (RadialRootIdentity, usize),
        observation: RadialRootObservation,
    ) {
        match self.observations.get(&key) {
            Some(previous)
                if observation.precision_bits < previous.precision_bits
                    || (observation.precision_bits == previous.precision_bits
                        && observation.residual >= previous.residual) => {}
            _ => {
                self.observations.insert(key, observation);
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn solve<T: FloatLike>(
        &mut self,
        identity: &RadialRootIdentity,
        inside_radius: &F<T>,
        outside_radius_guess: &F<T>,
        f_x_and_df_x: impl Fn(&F<T>) -> (F<T>, F<T>),
        tolerance: &F<T>,
        max_iterations: usize,
        max_bracket_expansions: usize,
        e_cm: &F<T>,
    ) -> Result<NewtonIterationResult<T>, SafeguardedNewtonError<T>> {
        let call_key = self.next_call_key(identity, inside_radius);
        let solve_result = safeguarded_newton_iteration_and_derivative(
            inside_radius,
            outside_radius_guess,
            &f_x_and_df_x,
            tolerance,
            max_iterations,
            max_bracket_expansions,
            e_cm,
        );

        let (result, lower_bound, upper_bound) = match &solve_result {
            Ok(result) => {
                let observation = RadialRootObservation::new(result, None, None, tolerance, e_cm);
                self.record_observation(call_key, observation);
                return solve_result;
            }
            Err(SafeguardedNewtonError::DidNotConverge {
                result,
                lower_bound,
                upper_bound,
            }) => (result, lower_bound, upper_bound),
            Err(error) => {
                debug!(
                    radial_root = %identity,
                    occurrence = call_key.1,
                    error = %error,
                    "radial root failed a structural safeguarded-solver check"
                );
                return solve_result;
            }
        };

        let previous = self.observations.get(&call_key).cloned();
        let current_precision: SymbolicaFloat = result.solution.clone().into();
        let current_precision_bits = current_precision.prec();
        let local_consistency = previous
            .as_ref()
            .is_some_and(|previous| current_precision_bits > previous.precision_bits)
            .then(|| LocalRootConsistency::check(result, inside_radius, &f_x_and_df_x, e_cm));
        let current = RadialRootObservation::new(
            result,
            Some((lower_bound, upper_bound)),
            local_consistency,
            tolerance,
            e_cm,
        );
        let precision_improvement = self
            .observations
            .get(&call_key)
            .and_then(|previous| current.precision_rescue_improvement(previous));

        if let Some(improvement) = precision_improvement {
            let previous = previous.expect("a precision improvement requires a prior observation");
            debug!(
                radial_root = %identity,
                occurrence = call_key.1,
                previous_epsilon = %previous.precision_epsilon,
                current_epsilon = %current.precision_epsilon,
                previous_residual = %previous.residual,
                previous_comparison_residual = %previous.comparison_residual(),
                current_residual = %current.residual,
                residual_improvement = %improvement,
                accepted_residual = %previous.maximum_residual,
                relative_newton_correction = %current.relative_newton_correction,
                local_near_radius_step = ?current
                    .local_consistency
                    .as_ref()
                    .and_then(|consistency| consistency.near_radius_step.as_ref()),
                local_near_lower_value = ?current
                    .local_consistency
                    .as_ref()
                    .and_then(|consistency| consistency.near_lower_value.as_ref()),
                local_near_upper_value = ?current
                    .local_consistency
                    .as_ref()
                    .and_then(|consistency| consistency.near_upper_value.as_ref()),
                local_near_secant_derivative_ratio = ?current
                    .local_consistency
                    .as_ref()
                    .and_then(|consistency| consistency.near_secant_derivative_ratio.as_ref()),
                local_far_secant_derivative_ratio = ?current
                    .local_consistency
                    .as_ref()
                    .and_then(|consistency| consistency.far_secant_derivative_ratio.as_ref()),
                local_secant_ratio_agreement = ?current
                    .local_consistency
                    .as_ref()
                    .and_then(|consistency| consistency.secant_ratio_agreement.as_ref()),
                "accepted a residual-limited radial root after precision escalation"
            );
            self.record_observation(call_key, current);
            return Ok(result.clone());
        }

        if let Some(previous) = previous
            .as_ref()
            .filter(|previous| current.precision_bits > previous.precision_bits)
        {
            let residual_improvement = previous.comparison_residual() / current.residual.clone();
            debug!(
                radial_root = %identity,
                occurrence = call_key.1,
                previous_epsilon = %previous.precision_epsilon,
                current_epsilon = %current.precision_epsilon,
                previous_residual = %previous.residual,
                previous_comparison_residual = %previous.comparison_residual(),
                current_residual = %current.residual,
                residual_improvement = %residual_improvement,
                accepted_residual = %previous.maximum_residual,
                relative_newton_correction = %current.relative_newton_correction,
                accepted_relative_newton_correction = %previous.relative_residual_limit,
                solution = %current.solution,
                derivative = %current.derivative,
                lower_bound = ?current.lower_bound,
                upper_bound = ?current.upper_bound,
                bracket_is_valid = current.bracket_is_valid,
                local_consistency = ?current.local_consistency,
                "radial root did not satisfy precision-rescue criteria"
            );
        } else {
            debug!(
                radial_root = %identity,
                occurrence = call_key.1,
                current_epsilon = %current.precision_epsilon,
                current_residual = %current.residual,
                "radial root has no matching lower-precision observation"
            );
        }

        self.record_observation(call_key, current);
        solve_result
    }
}

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
    use std::cell::Cell;

    use super::*;
    use crate::utils::f128;

    const ROOT_TOLERANCE: i64 = 64;

    fn cancellation_limited_root<T: FloatLike>(
        diagnostics: &mut RadialRootDiagnostics,
        identity: &RadialRootIdentity,
        large_momentum: f64,
        target: f64,
    ) -> Result<NewtonIterationResult<T>, SafeguardedNewtonError<T>> {
        let zero = F::<T>::default();
        let one = zero.one();
        let target = F::<T>::from_f64(target);
        let k3 = F::<T>::from_f64(large_momentum);
        let k2 = F::<T>::from_f64(large_momentum + 1.0);
        let derivative = &k2 - &k3;
        let tolerance = zero.from_i64(ROOT_TOLERANCE);

        diagnostics.solve(
            identity,
            &zero,
            &one,
            |radius| (radius * &k2 - radius * &k3 - &target, derivative.clone()),
            &tolerance,
            40,
            64,
            &one,
        )
    }

    fn discontinuous_root<T: FloatLike>(
        diagnostics: &mut RadialRootDiagnostics,
        identity: &RadialRootIdentity,
        positive_residual: f64,
    ) -> Result<NewtonIterationResult<T>, SafeguardedNewtonError<T>> {
        let zero = F::<T>::default();
        let one = zero.one();
        let half = &one / one.from_i64(2);
        let positive_residual = F::<T>::from_f64(positive_residual);
        let tolerance = zero.from_i64(ROOT_TOLERANCE);

        diagnostics.solve(
            identity,
            &zero,
            &one,
            |radius| {
                let value = if radius == &zero {
                    -one.clone()
                } else if radius < &half {
                    -positive_residual.clone()
                } else {
                    positive_residual.clone()
                };
                (value, one.clone())
            },
            &tolerance,
            40,
            64,
            &one,
        )
    }

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

    #[test]
    fn local_root_consistency_has_bounded_evaluation_cost() {
        let result = NewtonIterationResult {
            solution: F(0.5),
            derivative_at_solution: F(1.0),
            error_of_function: F(1.0e-20),
            num_iterations_used: 40,
        };
        let calls = Cell::new(0);
        let consistency = LocalRootConsistency::check(
            &result,
            &F(0.0),
            &|radius| {
                calls.set(calls.get() + 1);
                (radius - F(0.5), F(1.0))
            },
            &F(1.0),
        );
        assert!(consistency.is_valid);
        assert_eq!(calls.get(), 4);

        let short_solve_result = NewtonIterationResult {
            num_iterations_used: 3,
            ..result
        };
        calls.set(0);
        let consistency = LocalRootConsistency::check(
            &short_solve_result,
            &F(0.0),
            &|radius| {
                calls.set(calls.get() + 1);
                (radius - F(0.5), F(1.0))
            },
            &F(1.0),
        );
        assert!(!consistency.is_valid);
        assert_eq!(calls.get(), 0);
    }

    #[test]
    fn higher_precision_rescues_large_nearly_equal_momenta() {
        let mut diagnostics = RadialRootDiagnostics::default();

        for exponent in [32, 40, 48, 52] {
            let identity = RadialRootIdentity::new(format!("cancellation at 2^{exponent}"));
            let large_momentum = 2.0_f64.powi(exponent);
            let f64_error = cancellation_limited_root::<f64>(
                &mut diagnostics,
                &identity,
                large_momentum,
                1.0 / 3.0,
            )
            .unwrap_err();
            assert!(matches!(
                f64_error,
                SafeguardedNewtonError::DidNotConverge { .. }
            ));

            let rescued = cancellation_limited_root::<f128>(
                &mut diagnostics,
                &identity,
                large_momentum,
                1.0 / 3.0,
            )
            .unwrap();
            let active_precision_limit = rescued.solution.epsilon()
                * rescued.solution.from_i64(ROOT_TOLERANCE)
                * rescued.solution.one();
            assert!(
                rescued.error_of_function.abs() > active_precision_limit,
                "2^{exponent} case should exercise precision rescue rather than direct convergence"
            );
        }
    }

    #[test]
    fn higher_precision_can_use_a_successful_lower_precision_baseline() {
        let mut diagnostics = RadialRootDiagnostics::default();
        let identity = RadialRootIdentity::new("successful f64 baseline".to_string());
        let target = 0.3;

        let f64_result = diagnostics
            .solve(
                &identity,
                &F(0.0),
                &F(1.0),
                |radius| (radius - F(target), F(1.0)),
                &F(ROOT_TOLERANCE as f64),
                40,
                64,
                &F(1.0),
            )
            .unwrap();
        assert!(f64_result.error_of_function.abs() > F(0.0));
        assert!(
            f64_result.error_of_function.abs()
                <= f64_result.solution.epsilon()
                    * F(ROOT_TOLERANCE as f64)
                    * f64_result.solution.one()
        );

        // This is algebraically the same root as `r - target`, but separately evaluating
        // `(K + 1) r - K r` leaves a precision-limited residual for large K.
        let rescued = cancellation_limited_root::<f128>(
            &mut diagnostics,
            &identity,
            2.0_f64.powi(52),
            target,
        )
        .unwrap();
        let active_precision_limit =
            rescued.solution.epsilon() * rescued.solution.from_i64(ROOT_TOLERANCE);
        assert!(
            rescued.error_of_function.abs() > active_precision_limit,
            "the f128 solve should require the successful f64 observation as its baseline"
        );
    }

    #[test]
    fn higher_precision_uses_the_roundoff_floor_after_an_exact_lower_precision_zero() {
        let mut diagnostics = RadialRootDiagnostics::default();
        let identity = RadialRootIdentity::new("accidental exact f64 residual".to_string());

        let f64_result = diagnostics
            .solve(
                &identity,
                &F(0.0),
                &F(1.0),
                |radius| (radius - F(0.3), F(1.0)),
                &F(0.1),
                40,
                64,
                &F(1.0),
            )
            .unwrap();
        assert_eq!(f64_result.error_of_function, F(0.0));

        let rescued =
            cancellation_limited_root::<f128>(&mut diagnostics, &identity, 2.0_f64.powi(52), 0.3)
                .unwrap();
        let active_precision_limit =
            rescued.solution.epsilon() * rescued.solution.from_i64(ROOT_TOLERANCE);
        assert!(rescued.error_of_function.abs() > active_precision_limit);
    }

    #[test]
    fn repeated_root_identities_are_paired_by_precision_local_occurrence() {
        let mut diagnostics = RadialRootDiagnostics::default();
        let identity = RadialRootIdentity::new("two LMB channel rays".to_string());

        assert_eq!(
            diagnostics.next_call_key(&identity, &F::<f64>::default()),
            (identity.clone(), 0)
        );
        assert_eq!(
            diagnostics.next_call_key(&identity, &F::<f64>::default()),
            (identity.clone(), 1)
        );
        assert_eq!(
            diagnostics.next_call_key(&identity, &F::<f128>::default()),
            (identity.clone(), 0)
        );
        assert_eq!(
            diagnostics.next_call_key(&identity, &F::<f128>::default()),
            (identity, 1)
        );
    }

    #[test]
    fn precision_rescue_rejects_a_stable_discontinuous_residual() {
        let mut diagnostics = RadialRootDiagnostics::default();
        let identity = RadialRootIdentity::new("stable discontinuity".to_string());

        assert!(matches!(
            discontinuous_root::<f64>(&mut diagnostics, &identity, 1.0),
            Err(SafeguardedNewtonError::DidNotConverge { .. })
        ));
        assert!(matches!(
            discontinuous_root::<f128>(&mut diagnostics, &identity, 1.0),
            Err(SafeguardedNewtonError::DidNotConverge { .. })
        ));
    }

    #[test]
    fn precision_rescue_rejects_an_improved_but_large_residual() {
        let mut diagnostics = RadialRootDiagnostics::default();
        let identity = RadialRootIdentity::new("shrinking discontinuity".to_string());

        assert!(matches!(
            discontinuous_root::<f64>(&mut diagnostics, &identity, 1.0),
            Err(SafeguardedNewtonError::DidNotConverge { .. })
        ));
        assert!(matches!(
            discontinuous_root::<f128>(&mut diagnostics, &identity, 0.01),
            Err(SafeguardedNewtonError::DidNotConverge { .. })
        ));
    }

    #[test]
    fn precision_rescue_rejects_an_improving_discontinuous_non_root() {
        let mut diagnostics = RadialRootDiagnostics::default();
        let identity = RadialRootIdentity::new("improving discontinuity".to_string());
        let key = (identity.clone(), 0);

        assert!(matches!(
            discontinuous_root::<f64>(&mut diagnostics, &identity, 1.0e-12),
            Err(SafeguardedNewtonError::DidNotConverge { .. })
        ));
        let lower_precision = diagnostics.observations[&key].clone();

        assert!(matches!(
            discontinuous_root::<f128>(&mut diagnostics, &identity, 1.0e-14),
            Err(SafeguardedNewtonError::DidNotConverge { .. })
        ));
        let higher_precision = &diagnostics.observations[&key];

        // Residual improvement, the inherited f64 residual bound, and the Newton-correction
        // bound alone would accept this fabricated candidate even though the function has no
        // zero. The local sign/slope consistency check must be what rejects it.
        let residual_improvement =
            lower_precision.comparison_residual() / &higher_precision.residual;
        let required_improvement = SymbolicaFloat::with_val(
            residual_improvement.prec(),
            MINIMUM_PRECISION_RESIDUAL_IMPROVEMENT,
        );
        assert!(residual_improvement >= required_improvement);
        assert!(higher_precision.residual <= lower_precision.maximum_residual);
        assert!(
            higher_precision.relative_newton_correction <= lower_precision.relative_residual_limit
        );
        assert!(higher_precision.bracket_is_valid);
        assert!(
            !higher_precision
                .local_consistency
                .as_ref()
                .expect("residual-limited candidates have local diagnostics")
                .is_valid
        );
        assert!(
            higher_precision
                .precision_rescue_improvement(&lower_precision)
                .is_none()
        );
    }
}
