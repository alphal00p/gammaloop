//! Eager and first-derivative evaluators used by compiled sampling maps.
//!
//! Sampling maps are assembled from symbolic expressions during warmup.  This
//! wrapper keeps the expression compiler and its hyper-dual variant together so
//! a map and its Jacobian always use the same Symbolica expressions.  It is a
//! deliberately small layer around the process evaluator: graph-specific map
//! construction remains in the process sampler.

use color_eyre::Result;
use eyre::eyre;
use spenso::algebra::complex::Complex;
use symbolica::{
    atom::{Atom, AtomView},
    prelude::*,
};

use crate::{
    integrands::process::{GenericEvaluator, GenericEvaluatorFloat},
    processes::EvaluatorSettings,
    settings::runtime::{HFunction, HFunctionSettings, SamplingRadialProfile},
    utils::{F, FloatLike, hyperdual_utils::DualOrNot},
};
use symbolica::evaluate::OptimizationSettings;

use super::sampling_maps::SamplingEvaluationError;

/// Values and first partial derivatives of one evaluated expression.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingDualValue<T: FloatLike = f64> {
    /// The value of the expression at the supplied parameters.
    pub value: Complex<F<T>>,
    /// First derivatives in the compiled parameter-column order supplied to
    /// [`SamplingExpressionEvaluator::new`].
    pub derivatives: Vec<Complex<F<T>>>,
}

/// Real-valued output and first-derivative matrix of a sampling expression.
///
/// Sampling maps are real maps even though the shared Symbolica evaluator uses
/// its complex domain.  This type performs that boundary conversion once and
/// computes the signed determinant of a square output/input Jacobian.  The
/// map density must use [`Self::absolute_determinant`], since a coordinate
/// chart may reverse orientation.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingJacobianEvaluation<T: FloatLike = f64> {
    /// Expression values in the order supplied to the evaluator.
    pub values: Vec<T>,
    /// Rows are outputs; columns follow the requested active-parameter order
    /// (the compiled column order when no subset was supplied).
    pub derivatives: Vec<Vec<T>>,
    /// Signed determinant of `derivatives`.
    pub determinant: T,
}

impl<T: FloatLike> SamplingJacobianEvaluation<T> {
    /// Absolute determinant required by a push-forward volume density.
    pub fn absolute_determinant(&self) -> T {
        F(self.determinant.clone()).abs().0
    }
}

/// A warmup-compiled evaluator for one or more sampling-map expressions.
///
/// The eager program is always retained. If derivative columns are requested,
/// Symbolica's hyper-dual evaluator is built from the very same expression
/// tree, yielding exact first partials suitable for a forward-map Jacobian.
/// Cloning copies the already compiled programs and their mutable buffers;
/// it does not run expression optimization or compilation again.
#[derive(Clone)]
pub struct SamplingExpressionEvaluator {
    evaluator: GenericEvaluator,
    derivative_evaluator: Option<GenericEvaluator>,
    parameter_count: usize,
    output_count: usize,
    derivative_parameters: Vec<usize>,
}

impl SamplingExpressionEvaluator {
    /// Build an eager evaluator from expressions and their ordered parameters.
    ///
    /// Parameters are ordinary Symbolica atoms (usually variables).  The
    /// expression list must be non-empty; each expression contributes one
    /// output in the returned vector. Only `derivative_parameters` are dualized;
    /// an empty slice builds a scalar-only evaluator. Scalar reconstruction and
    /// active derivatives compile from the same expressions during warmup.
    pub fn new(
        expressions: impl IntoIterator<Item = Atom>,
        parameters: impl IntoIterator<Item = Atom>,
        derivative_parameters: &[usize],
    ) -> Result<Self> {
        let expressions = expressions.into_iter().collect::<Vec<_>>();
        let parameters = parameters.into_iter().collect::<Vec<_>>();
        if expressions.is_empty() {
            return Err(eyre!("sampling evaluator requires at least one expression"));
        }
        if parameters.is_empty() {
            return Err(eyre!("sampling evaluator requires at least one parameter"));
        }
        for (column, &parameter) in derivative_parameters.iter().enumerate() {
            if parameter >= parameters.len() {
                return Err(eyre!(
                    "sampling derivative parameter {parameter} is out of range"
                ));
            }
            if derivative_parameters[..column].contains(&parameter) {
                return Err(eyre!("sampling derivative parameters must be unique"));
            }
        }
        let compile = |dual_config| {
            GenericEvaluator::new_from_raw_params(
                expressions.clone(),
                &parameters,
                &FunctionMap::new(),
                Vec::new(),
                OptimizationSettings::default(),
                dual_config,
                &EvaluatorSettings::default(),
            )
        };
        let evaluator = compile(None)?;
        let derivative_evaluator = if derivative_parameters.is_empty() {
            None
        } else {
            let zero_components = (0..parameters.len())
                .flat_map(|parameter| {
                    derivative_parameters
                        .iter()
                        .enumerate()
                        .filter(move |&(_, &active)| active != parameter)
                        .map(move |(derivative, _)| (parameter, derivative + 1))
                })
                .collect();
            Some(compile(Some((
                first_derivative_shape(derivative_parameters.len()),
                zero_components,
            )))?)
        };
        Ok(Self {
            evaluator,
            derivative_evaluator,
            parameter_count: parameters.len(),
            output_count: expressions.len(),
            derivative_parameters: derivative_parameters.to_vec(),
        })
    }

    /// Compile a declared-positive proxy once for all supported precisions.
    ///
    /// Non-real and nonpositive constants are rejected before compilation.
    /// This is not a positivity certificate for arbitrary expressions: native
    /// score evaluation must still check every result, and the graph compiler
    /// must separately establish coverage of the intended singular features.
    pub(crate) fn new_positive_proxy(
        expression: Atom,
        parameters: impl IntoIterator<Item = Atom>,
    ) -> Result<Self> {
        if let AtomView::Num(number) = expression.as_view() {
            let coefficient = number.get_coeff_view().to_owned();
            if !coefficient.is_real() || coefficient.is_zero() || coefficient.is_negative() {
                return Err(eyre!(
                    "sampling Symbolica score constant must be real and strictly positive: {expression}"
                ));
            }
        }
        Self::new([expression], parameters, &[])
    }

    /// Compile the auxiliary-scale proposal once, including its native fit and
    /// first derivatives. The two sign inputs select algebraically equivalent
    /// stable tails: every evaluated exponential has a nonpositive argument.
    /// Inputs are log(t), log(R), log(beta), focused sign and broad sign;
    /// outputs are fitted log-scale, shape, CDF, survival and log(raw radius).
    pub(crate) fn new_lu_h_profile(
        profile: &SamplingRadialProfile,
        h: &HFunctionSettings,
    ) -> Result<Self> {
        if !h.sigma.is_finite() || h.sigma <= 0.0 {
            return Err(eyre!(
                "LU-h sampling requires a positive finite h_function.sigma"
            ));
        }
        if !profile.broad_fraction.is_finite()
            || !(0.0..=1.0).contains(&profile.broad_fraction)
            || [profile.scale, profile.shape]
                .into_iter()
                .flatten()
                .any(|value| !value.is_finite() || value <= 0.0)
        {
            return Err(eyre!(
                "LU-h sampling requires broad_fraction in [0,1] and positive finite scale/shape overrides"
            ));
        }
        let power = h.power.unwrap_or(0);
        let (log_scale, shape) = match h.function {
            HFunction::PolyExponential | HFunction::PolyLeftRightExponential => {
                if ![0, 1, 3, 4, 6, 7, 9, 10, 12, 13, 15, 16].contains(&power) {
                    return Err(eyre!(
                        "unsupported LU h-function power {power}; supported powers are 0,1,3,4,6,7,9,10,12,13,15,16"
                    ));
                }
                let a = if matches!(h.function, HFunction::PolyExponential) {
                    2
                } else {
                    1
                };
                let b = format!("((1-{power})/(2*{a}))");
                // asinh(b) written with positive arguments, evaluated in the
                // requested native precision rather than fitted in binary64.
                (
                    format!("log({})+log({b}+sqrt(1+({b})^2))/{a}", h.sigma),
                    format!("2*{a}*(1+({b})^2)^(1/4)"),
                )
            }
            HFunction::Exponential => (format!("log({}/2)", h.sigma), "1".to_owned()),
            HFunction::ExponentialCT => {
                return Err(eyre!(
                    "exponential_ct is a local threshold localization, not a normalized auxiliary LU h-function"
                ));
            }
        };
        let log_scale = profile
            .scale
            .map_or(log_scale, |scale| format!("log({scale})"));
        let shape = profile.shape.map_or(shape, |shape| shape.to_string());
        let parameters = [
            "lu_profile_y",
            "lu_profile_log_root",
            "lu_profile_log_beta",
            "lu_profile_focus_sign",
            "lu_profile_broad_sign",
        ];
        let x = format!("exp(lu_profile_focus_sign*({shape})*(lu_profile_y-({log_scale})))");
        let b = "exp(lu_profile_broad_sign*(lu_profile_y-lu_profile_log_root+lu_profile_log_beta))";
        let focused_cdf =
            format!("((1+lu_profile_focus_sign)*({x})+1-lu_profile_focus_sign)/(2*(1+({x})))");
        let focused_survival =
            format!("(1+lu_profile_focus_sign+(1-lu_profile_focus_sign)*({x}))/(2*(1+({x})))");
        let broad_cdf =
            format!("((1+lu_profile_broad_sign)*({b})+1-lu_profile_broad_sign)/(2*(1+({b})))");
        let broad_survival =
            format!("(1+lu_profile_broad_sign+(1-lu_profile_broad_sign)*({b}))/(2*(1+({b})))");
        let epsilon = profile.broad_fraction;
        let expressions = [
            log_scale,
            shape,
            format!("(1-{epsilon})*({focused_cdf})+{epsilon}*({broad_cdf})"),
            format!("(1-{epsilon})*({focused_survival})+{epsilon}*({broad_survival})"),
            "lu_profile_log_root-lu_profile_y".to_owned(),
        ];
        Self::new(
            expressions
                .iter()
                .map(|expression| try_parse!(expression))
                .collect::<std::result::Result<Vec<_>, _>>()
                .map_err(|error| eyre!(error))?,
            parameters
                .iter()
                .map(|parameter| try_parse!(*parameter))
                .collect::<std::result::Result<Vec<_>, _>>()
                .map_err(|error| eyre!(error))?,
            &[0, 1, 2, 3, 4],
        )
    }

    pub fn parameter_count(&self) -> usize {
        self.parameter_count
    }

    pub fn output_count(&self) -> usize {
        self.output_count
    }

    pub fn derivative_parameters(&self) -> &[usize] {
        &self.derivative_parameters
    }

    /// Evaluate all outputs eagerly at a real parameter point.
    pub fn evaluate<T: FloatLike>(&mut self, parameters: &[T]) -> Result<Vec<Complex<F<T>>>> {
        self.validate_parameters(parameters)?;
        let input = parameters
            .iter()
            .cloned()
            .map(|value| Complex::new_re(F(value)))
            .collect::<Vec<_>>();
        let values = <T as GenericEvaluatorFloat>::get_evaluator(&mut self.evaluator)(&input);
        values
            .into_iter()
            .map(|value| match value {
                DualOrNot::NonDual(value) => Ok(value),
                DualOrNot::Dual(_) => Err(eyre!("scalar sampling evaluator returned dual outputs")),
            })
            .collect()
    }

    /// Evaluate all outputs and their compiled first partial derivatives.
    pub fn evaluate_with_derivatives<T: FloatLike>(
        &mut self,
        parameters: &[T],
    ) -> Result<Vec<SamplingDualValue<T>>> {
        self.validate_parameters(parameters)?;
        let input = self.evaluator_parameters(parameters);
        let evaluator = self
            .derivative_evaluator
            .as_mut()
            .ok_or_else(|| eyre!("sampling evaluator was built without derivative support"))?;
        let values = <T as GenericEvaluatorFloat>::get_evaluator(evaluator)(&input);
        values
            .into_iter()
            .map(|value| match value {
                DualOrNot::Dual(value) => {
                    if value.values.len() != self.derivative_parameters.len() + 1 {
                        return Err(eyre!(
                            "sampling evaluator returned {} dual components, expected {}",
                            value.values.len(),
                            self.derivative_parameters.len() + 1
                        ));
                    }
                    Ok(SamplingDualValue {
                        value: value.values[0].clone(),
                        derivatives: value.values[1..].to_vec(),
                    })
                }
                DualOrNot::NonDual(_) => Err(eyre!(
                    "sampling evaluator did not return hyper-dual outputs"
                )),
            })
            .collect()
    }

    /// Evaluate a real map and its exact Symbolica dual Jacobian with respect
    /// to the ordered active parameters. `None` selects every compiled column.
    ///
    /// This is the common boundary used by map kernels: the expression tree
    /// is compiled once during warmup, while coordinates are supplied at
    /// runtime.  A non-real output is rejected rather than silently projected
    /// onto its real part, since doing so would invalidate the advertised
    /// volume Jacobian. Parameters outside the selected columns are held fixed;
    /// their values can still change the map and its active Jacobian. This uses
    /// the same compiled active dual program with statically zero identity-seed
    /// components to supply the requested columns. Prepared-only singular
    /// derivatives do not contaminate active columns. Singular intermediates
    /// that depend on active parameters must still have finite derivatives;
    /// this is not an algebraic extension through their removable singularities.
    pub fn evaluate_with_real_jacobian<T: FloatLike>(
        &mut self,
        parameters: &[T],
        active_parameters: Option<&[usize]>,
    ) -> Result<SamplingJacobianEvaluation<T>> {
        let active_parameters =
            active_parameters.map_or_else(|| self.derivative_parameters.clone(), <[usize]>::to_vec);
        if self.output_count != active_parameters.len() {
            return Err(eyre!(
                "sampling Jacobian requires a square map ({} outputs, {} active parameters)",
                self.output_count,
                active_parameters.len()
            ));
        }
        for (column, &parameter) in active_parameters.iter().enumerate() {
            if parameter >= self.parameter_count {
                return Err(eyre!(
                    "sampling Jacobian parameter index {parameter} is out of range for {} parameters",
                    self.parameter_count
                ));
            }
            if active_parameters[..column].contains(&parameter) {
                return Err(eyre!(
                    "sampling Jacobian active parameters must be unique; index {parameter} is repeated"
                ));
            }
        }
        let columns = active_parameters
            .iter()
            .map(|parameter| {
                self.derivative_parameters
                    .iter()
                    .position(|compiled| compiled == parameter)
                    .ok_or_else(|| {
                        eyre!("sampling derivative parameter {parameter} was not compiled")
                    })
            })
            .collect::<Result<Vec<_>>>()?;
        let dual_values = self.evaluate_with_derivatives(parameters)?;
        let mut values = Vec::with_capacity(dual_values.len());
        let mut derivatives = Vec::with_capacity(dual_values.len());
        for (output_index, value) in dual_values.into_iter().enumerate() {
            values.push(Self::real_component(
                value.value,
                "value",
                output_index,
                None,
            )?);
            let row = columns
                .iter()
                .zip(&active_parameters)
                .map(|(&column, &parameter_index)| {
                    Self::real_component(
                        value.derivatives[column].clone(),
                        "derivative",
                        output_index,
                        Some(parameter_index),
                    )
                })
                .collect::<Result<Vec<_>>>()?;
            derivatives.push(row);
        }
        let determinant = signed_determinant(&derivatives)?;
        Ok(SamplingJacobianEvaluation {
            values,
            derivatives,
            determinant,
        })
    }

    fn validate_parameters<T: FloatLike>(&self, parameters: &[T]) -> Result<()> {
        if parameters.len() != self.parameter_count {
            return Err(eyre!(
                "sampling evaluator received {} parameters, expected {}",
                parameters.len(),
                self.parameter_count
            ));
        }
        if parameters.iter().any(|value| !value.is_finite()) {
            return Err(eyre!("sampling evaluator parameters must be finite"));
        }
        Ok(())
    }

    /// Symbolica vectorizes a dual expression by storing one hyper-dual block
    /// for every scalar parameter. Keep this flattening local to the shared
    /// evaluator wrapper so eager values and dual Jacobians cannot disagree on
    /// the input layout.
    fn evaluator_parameters<T: FloatLike>(&self, parameters: &[T]) -> Vec<Complex<F<T>>> {
        let dual_width = self.derivative_parameters.len() + 1;
        let mut input = Vec::with_capacity(self.parameter_count * dual_width);
        for (parameter_index, value) in parameters.iter().cloned().enumerate() {
            let value = F(value);
            input.push(Complex::new_re(value.clone()));
            for &derivative_parameter in &self.derivative_parameters {
                input.push(Complex::new_re(
                    if parameter_index == derivative_parameter {
                        value.one()
                    } else {
                        value.zero()
                    },
                ));
            }
        }
        debug_assert_eq!(input.len(), self.parameter_count * dual_width);
        input
    }

    pub(crate) fn real_component<T: FloatLike>(
        value: Complex<F<T>>,
        role: &str,
        output_index: usize,
        parameter_index: Option<usize>,
    ) -> Result<T> {
        // Expressions in the map language are required to be real. Permit only
        // roundoff-level imaginary residue from complex constants, and retain the
        // diagnostic location when a genuinely complex map is supplied. Scale
        // relative to the native output, including values smaller than f64 can hold.
        let tolerance = value.re.epsilon() * value.re.from_usize(512) * value.re.abs();
        if !value.re.0.is_finite() || !value.im.0.is_finite() {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "eager sampling expression",
                detail: format!(
                    "sampling {role} {output_index}{} is non-finite: {value:?}",
                    parameter_index
                        .map(|index| format!("/{index}"))
                        .unwrap_or_default()
                ),
            }
            .into());
        }
        if value.im.abs() > tolerance {
            return Err(eyre!(
                "sampling {role} {output_index}{} has non-negligible imaginary part {}",
                parameter_index
                    .map(|index| format!("/{index}"))
                    .unwrap_or_default(),
                value.im
            ));
        }
        Ok(value.re.0)
    }
}

fn signed_determinant<T: FloatLike>(matrix: &[Vec<T>]) -> Result<T> {
    let dimension = matrix.len();
    if dimension == 0 || matrix.iter().any(|row| row.len() != dimension) {
        return Err(eyre!(
            "sampling Jacobian determinant requires a non-empty square matrix"
        ));
    }
    let mut work = matrix
        .iter()
        .map(|row| row.iter().cloned().map(F).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    if work.iter().flatten().any(|value| !value.0.is_finite()) {
        return Err(SamplingEvaluationError::Unrepresentable {
            operation: "sampling Jacobian determinant",
            detail: "matrix contains non-finite entries".into(),
        }
        .into());
    }
    let mut determinant = work[0][0].one();
    for pivot in 0..dimension {
        if (pivot..dimension).any(|row| !work[row][pivot].0.is_finite()) {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "sampling Jacobian determinant",
                detail: "elimination produced a non-finite pivot".into(),
            }
            .into());
        }
        let (pivot_row, pivot_value) = (pivot..dimension)
            .map(|row| (row, work[row][pivot].abs()))
            .max_by(|(_, left), (_, right)| left.partial_cmp(right).expect("finite pivots"))
            .expect("non-empty pivot range");
        if pivot_value == pivot_value.zero() {
            return Ok(pivot_value.0);
        }
        if pivot_row != pivot {
            work.swap(pivot_row, pivot);
            determinant = -determinant;
        }
        let diagonal = work[pivot][pivot].clone();
        determinant *= &diagonal;
        if !determinant.0.is_finite() || determinant == determinant.zero() {
            return Err(SamplingEvaluationError::Unrepresentable {
                operation: "sampling Jacobian determinant",
                detail: "product of nonzero pivots overflowed or underflowed".into(),
            }
            .into());
        }
        let (completed, pending) = work.split_at_mut(pivot + 1);
        let pivot_values = &completed[pivot];
        for row in pending {
            let factor = &row[pivot] / &diagonal;
            for (entry, pivot_entry) in row.iter_mut().zip(pivot_values).skip(pivot + 1) {
                *entry -= &factor * pivot_entry;
            }
        }
    }
    Ok(determinant.0)
}

fn first_derivative_shape(parameter_count: usize) -> Vec<Vec<usize>> {
    let mut shape = Vec::with_capacity(parameter_count + 1);
    shape.push(vec![0; parameter_count]);
    for derivative in 0..parameter_count {
        let mut order = vec![0; parameter_count];
        order[derivative] = 1;
        shape.push(order);
    }
    shape
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::initialisation::test_initialise;

    #[test]
    fn sampling_lu_h_profile_rejects_unsupported_physical_and_proposal_settings() {
        let profile = SamplingRadialProfile::default();
        for h in [
            HFunctionSettings {
                function: HFunction::ExponentialCT,
                ..Default::default()
            },
            HFunctionSettings {
                power: Some(2),
                ..Default::default()
            },
            HFunctionSettings {
                sigma: 0.0,
                ..Default::default()
            },
        ] {
            assert!(SamplingExpressionEvaluator::new_lu_h_profile(&profile, &h).is_err());
            // Broad-only sampling still validates the actual inherited h at
            // warmup, even though its map subsequently skips the root solve.
            assert!(
                SamplingExpressionEvaluator::new_lu_h_profile(
                    &SamplingRadialProfile {
                        broad_fraction: 1.0,
                        ..profile.clone()
                    },
                    &h
                )
                .is_err()
            );
        }
        for invalid in [
            SamplingRadialProfile {
                broad_fraction: -0.1,
                ..profile.clone()
            },
            SamplingRadialProfile {
                broad_fraction: 1.1,
                ..profile.clone()
            },
            SamplingRadialProfile {
                scale: Some(0.0),
                ..profile.clone()
            },
            SamplingRadialProfile {
                shape: Some(f64::INFINITY),
                ..profile.clone()
            },
        ] {
            assert!(
                SamplingExpressionEvaluator::new_lu_h_profile(
                    &invalid,
                    &HFunctionSettings::default()
                )
                .is_err()
            );
        }
    }

    #[test]
    fn eager_evaluation_returns_all_outputs() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("sampling_test::x"));
        let y = Atom::var(symbol!("sampling_test::y"));
        let mut evaluator = SamplingExpressionEvaluator::new(
            [x.clone() + y.clone(), x.clone() * y.clone()],
            [x, y],
            &[],
        )
        .unwrap();
        let values = evaluator.evaluate(&[2.0, 3.0]).unwrap();
        assert_eq!(values.len(), 2);
        assert_eq!(values[0], Complex::new_re(F(5.0)));
        assert_eq!(values[1], Complex::new_re(F(6.0)));
    }

    #[test]
    fn dual_evaluation_returns_symbolica_partials() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("sampling_test::dual_x"));
        let y = Atom::var(symbol!("sampling_test::dual_y"));
        let mut evaluator = SamplingExpressionEvaluator::new(
            [x.clone() * x.clone() + x.clone() * y.clone()],
            [x, y],
            &[0, 1],
        )
        .unwrap();
        let values = evaluator.evaluate_with_derivatives(&[2.0, 3.0]).unwrap();
        let [value] = values.as_slice() else {
            panic!("expected one output")
        };
        assert_eq!(value.value, Complex::new_re(F(10.0)));
        assert_eq!(
            value.derivatives,
            vec![Complex::new_re(F(7.0)), Complex::new_re(F(2.0))]
        );
        let values = evaluator.evaluate(&[2.0, 3.0]).unwrap();
        assert_eq!(values, vec![Complex::new_re(F(10.0))]);
    }

    #[test]
    fn real_jacobian_uses_symbolica_duals_and_signed_determinant() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("sampling_test::jac_x"));
        let y = Atom::var(symbol!("sampling_test::jac_y"));
        let mut evaluator = SamplingExpressionEvaluator::new(
            [
                x.clone() + x.clone() * y.clone(),
                y.clone() + x.clone() * x.clone(),
            ],
            [x, y],
            &[0, 1],
        )
        .unwrap();
        let jacobian = evaluator
            .evaluate_with_real_jacobian(&[2.0, 3.0], None)
            .unwrap();
        assert_eq!(jacobian.values, vec![8.0, 7.0]);
        assert_eq!(jacobian.derivatives, vec![vec![4.0, 2.0], vec![4.0, 1.0]]);
        assert_eq!(jacobian.determinant, -4.0);
        assert_eq!(jacobian.absolute_determinant(), 4.0);
    }

    #[test]
    fn real_jacobian_holds_prepared_parameters_fixed_in_native_precision() {
        test_initialise().unwrap();
        // Interleave prepared scale/shear with the two active cube coordinates.
        // The active determinant is s + s^2*h - h^2; differentiating prepared
        // parameters instead produces a different matrix and volume factor.
        let mut evaluator = SamplingExpressionEvaluator::new(
            [
                try_parse!("sampling_active::s*sampling_active::u + sampling_active::h*sampling_active::v").unwrap(),
                try_parse!("sampling_active::h*sampling_active::u + (1+sampling_active::s*sampling_active::h)*sampling_active::v").unwrap(),
            ],
            ["sampling_active::s", "sampling_active::u", "sampling_active::h", "sampling_active::v"]
                .map(|name| try_parse!(name).unwrap()),
            &[0, 1, 2, 3],
        ).unwrap();
        let parameters = [2.0, 0.25, 3.0, 0.5];
        let evaluated = evaluator
            .evaluate_with_real_jacobian(&parameters, Some(&[1, 3]))
            .unwrap();
        assert_eq!(evaluated.values, vec![2.0, 4.25]);
        assert_eq!(evaluated.derivatives, vec![vec![2.0, 3.0], vec![3.0, 7.0]]);
        assert!((evaluated.determinant - 5.0).abs() < 1.0e-13);
        let reversed = evaluator
            .evaluate_with_real_jacobian(&parameters, Some(&[3, 1]))
            .unwrap();
        assert_eq!(reversed.values, evaluated.values);
        assert_eq!(reversed.derivatives, vec![vec![3.0, 2.0], vec![7.0, 3.0]]);
        assert!((reversed.determinant + 5.0).abs() < 1.0e-13);
        let changed = evaluator
            .evaluate_with_real_jacobian(&[5.0, 0.25, 2.0, 0.5], Some(&[1, 3]))
            .unwrap();
        assert_eq!(changed.values, vec![2.25, 6.0]);
        assert!((changed.determinant - 51.0).abs() < 1.0e-12);
        for (columns, message) in [
            (&[1, 1][..], "unique"),
            (&[1, 4][..], "out of range"),
            (&[1][..], "square map"),
        ] {
            assert!(
                evaluator
                    .evaluate_with_real_jacobian(&parameters, Some(columns))
                    .unwrap_err()
                    .to_string()
                    .contains(message)
            );
        }

        fn check<T: FloatLike>(evaluator: &mut SamplingExpressionEvaluator) {
            let one = F::<T>::default().one();
            let s = one.from_usize(2) + &one / one.from_usize(10).powi(25);
            let h = one.from_usize(3);
            let u = &one / one.from_usize(4);
            let v = &one / one.from_usize(2);
            let expected = &s + s.square() * &h - h.square();
            let evaluated = evaluator
                .evaluate_with_real_jacobian(
                    &[s.0.clone(), u.0.clone(), h.0.clone(), v.0.clone()],
                    Some(&[1, 3]),
                )
                .unwrap();
            let tolerance = one.epsilon() * one.from_usize(2048);
            assert!((F(evaluated.values[0].clone()) - (&s * &u + &h * &v)).abs() < tolerance);
            assert!(
                (F(evaluated.values[1].clone()) - (&h * &u + (&one + &s * &h) * &v)).abs()
                    < tolerance
            );
            let determinant = F(evaluated.determinant);
            assert!((&determinant - expected).abs() < tolerance);
            assert!(determinant > one.from_usize(5));
        }
        check::<crate::utils::QuadFloat>(&mut evaluator);
        check::<crate::utils::ArbPrec>(&mut evaluator);
    }

    #[test]
    fn real_jacobian_keeps_singular_prepared_subexpressions_out_of_active_columns() {
        test_initialise().unwrap();
        let parameters = [parse!("sampling_fixed::u"), parse!("sampling_fixed::m")];
        let mut prepared = SamplingExpressionEvaluator::new(
            [parse!("sampling_fixed::u + sqrt(sampling_fixed::m)")],
            parameters.clone(),
            &[0, 1],
        )
        .unwrap();
        let mut mixed = SamplingExpressionEvaluator::new(
            [parse!(
                "sampling_fixed::u + sqrt(sampling_fixed::m*sampling_fixed::u)"
            )],
            parameters,
            &[0, 1],
        )
        .unwrap();

        fn check<T: FloatLike>(
            prepared: &mut SamplingExpressionEvaluator,
            mixed: &mut SamplingExpressionEvaluator,
        ) {
            let one = F::<T>::default().one();
            let u = &one / one.from_usize(2);
            let point = [u.0.clone(), one.zero().0];
            let evaluated = prepared
                .evaluate_with_real_jacobian(&point, Some(&[0]))
                .unwrap();
            assert_eq!(evaluated.values, vec![u.0.clone()]);
            assert_eq!(evaluated.derivatives, vec![vec![one.0.clone()]]);
            assert_eq!(evaluated.determinant, one.0);
            // The m derivative remains singular when it is actually requested.
            let error = prepared
                .evaluate_with_real_jacobian(&point, Some(&[1]))
                .unwrap_err();
            assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());

            // Static seed information does not solve removable singularities
            // whose intermediate argument still depends on an active input.
            assert_eq!(mixed.evaluate(&point).unwrap()[0].re, u);
            let error = mixed
                .evaluate_with_real_jacobian(&point, Some(&[0]))
                .unwrap_err();
            assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());
        }
        check::<f64>(&mut prepared, &mut mixed);
        check::<crate::utils::QuadFloat>(&mut prepared, &mut mixed);
        check::<crate::utils::ArbPrec>(&mut prepared, &mut mixed);
    }

    #[test]
    fn compiled_active_columns_share_a_finite_scalar_program() {
        test_initialise().unwrap();
        let parameters = ["u", "m", "v", "s", "w"]
            .map(|name| Atom::var(symbol!(format!("sampling_reduced::{name}"))));
        let [u, m, v, s, w] = &parameters;
        let expressions = [s * u + m.sqrt(), v + u * u, w + u * v];
        for (columns, message) in [(&[0, 0][..], "unique"), (&[5][..], "out of range")] {
            let error =
                SamplingExpressionEvaluator::new(expressions.clone(), parameters.clone(), columns)
                    .err()
                    .expect("invalid derivative columns must fail before compilation");
            assert!(error.to_string().contains(message));
        }
        let mut full = SamplingExpressionEvaluator::new(
            expressions.clone(),
            parameters.clone(),
            &[0, 1, 2, 3, 4],
        )
        .unwrap();
        let mut evaluator =
            SamplingExpressionEvaluator::new(expressions, parameters, &[0, 2, 4]).unwrap();
        assert!(evaluator.evaluator.dual_shape.is_none());
        assert_eq!(
            evaluator
                .derivative_evaluator
                .as_ref()
                .unwrap()
                .dual_shape
                .as_ref()
                .unwrap()
                .len(),
            4
        );

        fn check<T: FloatLike>(
            evaluator: &mut SamplingExpressionEvaluator,
            full: &mut SamplingExpressionEvaluator,
        ) {
            let one = F::<T>::default().one();
            let u = &one / one.from_usize(4);
            let v = &one / one.from_usize(2);
            let w = &one * one.from_usize(3) / one.from_usize(4);
            let s = one.from_usize(2) + &one / one.from_usize(10).powi(25);
            let point = [
                u.0.clone(),
                one.zero().0,
                v.0.clone(),
                s.0.clone(),
                w.0.clone(),
            ];
            // sqrt(m) at m=0 has a finite value and no active tangent. The
            // scalar path never executes its singular unused m derivative.
            let scalar = evaluator.evaluate(&point).unwrap();
            let jacobian = evaluator.evaluate_with_real_jacobian(&point, None).unwrap();
            let previous = full
                .evaluate_with_real_jacobian(&point, Some(&[0, 2, 4]))
                .unwrap();
            assert_eq!(jacobian, previous);
            for (scalar, previous) in scalar.iter().zip(&previous.values) {
                assert_eq!(&scalar.re.0, previous);
            }
            let tolerance = one.epsilon() * one.from_usize(2048);
            for ((scalar, value), expected) in
                scalar
                    .iter()
                    .zip(&jacobian.values)
                    .zip([&s * &u, &v + u.square(), &w + &u * &v])
            {
                assert_eq!(scalar.im, one.zero());
                assert!((&scalar.re - &expected).abs() < tolerance);
                assert!((F(value.clone()) - expected).abs() < tolerance);
            }
            assert!((F(jacobian.determinant) - &s).abs() < tolerance);
            let reordered = evaluator
                .evaluate_with_real_jacobian(&point, Some(&[4, 2, 0]))
                .unwrap();
            assert!((F(reordered.determinant) + &s).abs() < tolerance);
            let error = evaluator
                .evaluate_with_real_jacobian(&point, Some(&[0, 1, 4]))
                .unwrap_err();
            assert!(error.to_string().contains("was not compiled"));
        }
        check::<f64>(&mut evaluator, &mut full);
        check::<crate::utils::QuadFloat>(&mut evaluator, &mut full);
        check::<crate::utils::ArbPrec>(&mut evaluator, &mut full);

        // Independent Cartesian differences use only the scalar program and
        // vary the three compiled input columns, keeping m=0 and s fixed.
        let point = [0.25, 0.0, 0.5, 2.0, 0.75];
        let jacobian = evaluator.evaluate_with_real_jacobian(&point, None).unwrap();
        let step = 1e-5;
        for (column, parameter) in [0, 2, 4].into_iter().enumerate() {
            let mut plus = point;
            let mut minus = point;
            plus[parameter] += step;
            minus[parameter] -= step;
            let plus = evaluator.evaluate(&plus).unwrap();
            let minus = evaluator.evaluate(&minus).unwrap();
            for (row, (plus, minus)) in plus.iter().zip(minus).enumerate() {
                let difference = (plus.re.0 - minus.re.0) / (2.0 * step);
                assert!((difference - jacobian.derivatives[row][column]).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn real_jacobian_rejects_nonsquare_maps_and_complex_outputs() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("sampling_test::jac_nonsquare_x"));
        let y = Atom::var(symbol!("sampling_test::jac_nonsquare_y"));
        let mut nonsquare =
            SamplingExpressionEvaluator::new([x.clone()], [x.clone(), y], &[0, 1]).unwrap();
        let error = nonsquare
            .evaluate_with_real_jacobian(&[1.0, 2.0], None)
            .unwrap_err();
        assert!(error.to_string().contains("square map"));

        let complex_x = Atom::var(symbol!("sampling_test::jac_complex_x"));
        let complex_y = Atom::var(symbol!("sampling_test::jac_complex_y"));
        let mut complex = SamplingExpressionEvaluator::new(
            [complex_x.clone() + Atom::i(), complex_y.clone()],
            [complex_x, complex_y],
            &[0, 1],
        )
        .unwrap();
        let error = complex
            .evaluate_with_real_jacobian(&[1.0, 2.0], None)
            .unwrap_err();
        assert!(error.to_string().contains("imaginary part"));
    }

    #[test]
    fn derivative_request_is_rejected_for_eager_only_evaluator() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("sampling_test::x_only"));
        let mut evaluator = SamplingExpressionEvaluator::new([x.clone()], [x], &[]).unwrap();
        assert!(evaluator.evaluate_with_derivatives(&[1.0]).is_err());
    }

    #[test]
    fn native_eager_duals_and_determinants_preserve_values_outside_f64() {
        test_initialise().unwrap();
        let x = try_parse!("sampling_native::x").unwrap();
        let y = try_parse!("sampling_native::y").unwrap();
        let mut evaluator = SamplingExpressionEvaluator::new(
            [
                try_parse!("10^400 * sampling_native::x").unwrap(),
                try_parse!("10^-400 * sampling_native::y").unwrap(),
            ],
            [x, y],
            &[0, 1],
        )
        .unwrap();
        let error = evaluator
            .evaluate_with_real_jacobian(&[1.0, 1.0], None)
            .unwrap_err();
        assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());

        // QuadFloat extends the mantissa using two binary64 numbers; its
        // exponent range cannot represent these coefficients either.
        let quad_one = F::<crate::utils::QuadFloat>::default().one().0;
        let error = evaluator
            .evaluate_with_real_jacobian(&[quad_one, quad_one], None)
            .unwrap_err();
        assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());

        fn check<T: FloatLike>(evaluator: &mut SamplingExpressionEvaluator) {
            let one = F::<T>::default().one();
            let large = one.from_usize(10).powi(400);
            let small = &one / &large;
            let tolerance = one.epsilon() * one.from_usize(2048);
            // Reuse the same compiled program after its f64 failure.
            let evaluated = evaluator
                .evaluate_with_real_jacobian(&[one.0.clone(), one.0.clone()], None)
                .unwrap();
            assert!((F(evaluated.values[0].clone()) / &large - &one).abs() < tolerance);
            assert!((F(evaluated.values[1].clone()) / &small - &one).abs() < tolerance);
            assert!((F(evaluated.derivatives[0][0].clone()) / &large - &one).abs() < tolerance);
            assert!((F(evaluated.derivatives[1][1].clone()) / &small - &one).abs() < tolerance);
            assert!((F(evaluated.determinant) - &one).abs() < tolerance);

            // A tiny imaginary output remains non-real despite lying below
            // every f64 reporting scale and the old fixed absolute tolerance.
            let error = SamplingExpressionEvaluator::real_component(
                Complex::new(small.clone(), small),
                "native test",
                0,
                None,
            )
            .unwrap_err();
            assert!(error.to_string().contains("imaginary part"));
            assert!(error.downcast_ref::<SamplingEvaluationError>().is_none());
        }
        check::<crate::utils::ArbPrec>(&mut evaluator);
    }

    #[test]
    fn native_dual_seeds_resolve_sub_f64_changes() {
        test_initialise().unwrap();
        let x = try_parse!("sampling_native::perturbed_x").unwrap();
        let mut evaluator =
            SamplingExpressionEvaluator::new([x.clone() * x.clone()], [x], &[0]).unwrap();
        fn check<T: FloatLike>(evaluator: &mut SamplingExpressionEvaluator) {
            let one = F::<T>::default().one();
            let delta = &one / one.from_usize(10).powi(25);
            let x = &one + &delta;
            let value = evaluator
                .evaluate_with_derivatives(std::slice::from_ref(&x.0))
                .unwrap()
                .remove(0);
            let two = one.from_usize(2);
            let tolerance = one.epsilon() * one.from_usize(64);
            assert!((value.value.re - x.square()).abs() < tolerance);
            assert!((value.derivatives[0].re.clone() - &two * &x).abs() < tolerance);
            assert!(value.derivatives[0].re > two);
        }
        check::<crate::utils::QuadFloat>(&mut evaluator);
        check::<crate::utils::ArbPrec>(&mut evaluator);
    }
}
