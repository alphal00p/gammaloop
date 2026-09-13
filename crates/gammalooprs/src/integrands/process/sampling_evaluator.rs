//! Eager and first-derivative evaluators used by compiled sampling maps.
//!
//! Sampling maps are assembled from symbolic expressions during warmup.  This
//! wrapper keeps the expression compiler and its hyper-dual variant together so
//! a map and its Jacobian always use the same Symbolica program.  It is a
//! deliberately small layer around the process evaluator: graph-specific map
//! construction remains in the process sampler.

use color_eyre::Result;
use eyre::eyre;
use spenso::algebra::complex::Complex;
use symbolica::{atom::Atom, prelude::*};

use crate::{
    integrands::process::{GenericEvaluator, GenericEvaluatorFloat},
    processes::EvaluatorSettings,
    utils::{F, FloatLike, hyperdual_utils::DualOrNot},
};
use symbolica::evaluate::OptimizationSettings;

use super::sampling_maps::SamplingEvaluationError;

/// Values and first partial derivatives of one evaluated expression.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingDualValue<T: FloatLike = f64> {
    /// The value of the expression at the supplied parameters.
    pub value: Complex<F<T>>,
    /// First derivatives in the same order as the parameters passed to
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
    /// `derivatives[output][parameter]` from Symbolica's first-order duals.
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
/// The eager program is always retained.  If `with_derivatives` is requested,
/// Symbolica's hyper-dual evaluator is built from the very same expression
/// tree, yielding exact first partials suitable for a forward-map Jacobian.
/// Cloning copies the already compiled programs and their mutable buffers;
/// it does not run expression optimization or compilation again.
#[derive(Clone)]
pub struct SamplingExpressionEvaluator {
    evaluator: GenericEvaluator,
    parameter_count: usize,
    output_count: usize,
    with_derivatives: bool,
}

impl SamplingExpressionEvaluator {
    /// Build an eager evaluator from expressions and their ordered parameters.
    ///
    /// Parameters are ordinary Symbolica atoms (usually variables).  The
    /// expression list must be non-empty; each expression contributes one
    /// output in the returned vector.
    pub fn new(
        expressions: impl IntoIterator<Item = Atom>,
        parameters: impl IntoIterator<Item = Atom>,
        with_derivatives: bool,
    ) -> Result<Self> {
        let expressions = expressions.into_iter().collect::<Vec<_>>();
        let parameters = parameters.into_iter().collect::<Vec<_>>();
        if expressions.is_empty() {
            return Err(eyre!("sampling evaluator requires at least one expression"));
        }
        if parameters.is_empty() {
            return Err(eyre!("sampling evaluator requires at least one parameter"));
        }
        let dual_shape = with_derivatives.then(|| first_derivative_shape(parameters.len()));
        let evaluator = GenericEvaluator::new_from_raw_params(
            expressions.clone(),
            &parameters,
            &FunctionMap::new(),
            Vec::new(),
            OptimizationSettings::default(),
            dual_shape,
            &EvaluatorSettings::default(),
        )?;
        Ok(Self {
            evaluator,
            parameter_count: parameters.len(),
            output_count: expressions.len(),
            with_derivatives,
        })
    }

    pub fn parameter_count(&self) -> usize {
        self.parameter_count
    }

    pub fn output_count(&self) -> usize {
        self.output_count
    }

    pub fn has_derivatives(&self) -> bool {
        self.with_derivatives
    }

    /// Evaluate all outputs eagerly at a real parameter point.
    pub fn evaluate<T: FloatLike>(&mut self, parameters: &[T]) -> Result<Vec<Complex<F<T>>>> {
        self.validate_parameters(parameters)?;
        let input = self.evaluator_parameters(parameters);
        let values = <T as GenericEvaluatorFloat>::get_evaluator(&mut self.evaluator)(&input);
        values
            .into_iter()
            .map(|value| match value {
                DualOrNot::NonDual(value) => Ok(value),
                DualOrNot::Dual(value) => value
                    .values
                    .first()
                    .cloned()
                    .ok_or_else(|| eyre!("sampling evaluator returned an empty dual value")),
            })
            .collect()
    }

    /// Evaluate all outputs and their first partial derivatives.
    pub fn evaluate_with_derivatives<T: FloatLike>(
        &mut self,
        parameters: &[T],
    ) -> Result<Vec<SamplingDualValue<T>>> {
        if !self.with_derivatives {
            return Err(eyre!(
                "sampling evaluator was built without derivative support"
            ));
        }
        self.validate_parameters(parameters)?;
        let input = self.evaluator_parameters(parameters);
        let values = <T as GenericEvaluatorFloat>::get_evaluator(&mut self.evaluator)(&input);
        values
            .into_iter()
            .map(|value| match value {
                DualOrNot::Dual(value) => {
                    if value.values.len() != self.parameter_count + 1 {
                        return Err(eyre!(
                            "sampling evaluator returned {} dual components, expected {}",
                            value.values.len(),
                            self.parameter_count + 1
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

    /// Evaluate a real square map and its exact Symbolica dual Jacobian.
    ///
    /// This is the common boundary used by map kernels: the expression tree
    /// is compiled once during warmup, while coordinates are supplied at
    /// runtime.  A non-real output is rejected rather than silently projected
    /// onto its real part, since doing so would invalidate the advertised
    /// volume Jacobian.
    pub fn evaluate_with_real_jacobian<T: FloatLike>(
        &mut self,
        parameters: &[T],
    ) -> Result<SamplingJacobianEvaluation<T>> {
        if self.output_count != self.parameter_count {
            return Err(eyre!(
                "sampling Jacobian requires a square map ({} outputs, {} parameters)",
                self.output_count,
                self.parameter_count
            ));
        }
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
            let row = value
                .derivatives
                .into_iter()
                .enumerate()
                .map(|(parameter_index, derivative)| {
                    Self::real_component(
                        derivative,
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
        if !self.with_derivatives {
            return parameters
                .iter()
                .cloned()
                .map(|value| Complex::new_re(F(value)))
                .collect();
        }
        let dual_width = self.parameter_count + 1;
        let mut input = Vec::with_capacity(self.parameter_count * dual_width);
        for (parameter_index, value) in parameters.iter().cloned().enumerate() {
            let value = F(value);
            input.push(Complex::new_re(value.clone()));
            for derivative_index in 0..self.parameter_count {
                input.push(Complex::new_re(if parameter_index == derivative_index {
                    value.one()
                } else {
                    value.zero()
                }));
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
    fn eager_evaluation_returns_all_outputs() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("sampling_test::x"));
        let y = Atom::var(symbol!("sampling_test::y"));
        let mut evaluator = SamplingExpressionEvaluator::new(
            [x.clone() + y.clone(), x.clone() * y.clone()],
            [x, y],
            false,
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
            true,
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
            true,
        )
        .unwrap();
        let jacobian = evaluator.evaluate_with_real_jacobian(&[2.0, 3.0]).unwrap();
        assert_eq!(jacobian.values, vec![8.0, 7.0]);
        assert_eq!(jacobian.derivatives, vec![vec![4.0, 2.0], vec![4.0, 1.0]]);
        assert_eq!(jacobian.determinant, -4.0);
        assert_eq!(jacobian.absolute_determinant(), 4.0);
    }

    #[test]
    fn real_jacobian_rejects_nonsquare_maps_and_complex_outputs() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("sampling_test::jac_nonsquare_x"));
        let y = Atom::var(symbol!("sampling_test::jac_nonsquare_y"));
        let mut nonsquare =
            SamplingExpressionEvaluator::new([x.clone()], [x.clone(), y], true).unwrap();
        let error = nonsquare
            .evaluate_with_real_jacobian(&[1.0, 2.0])
            .unwrap_err();
        assert!(error.to_string().contains("square map"));

        let complex_x = Atom::var(symbol!("sampling_test::jac_complex_x"));
        let complex_y = Atom::var(symbol!("sampling_test::jac_complex_y"));
        let mut complex = SamplingExpressionEvaluator::new(
            [complex_x.clone() + Atom::i(), complex_y.clone()],
            [complex_x, complex_y],
            true,
        )
        .unwrap();
        let error = complex
            .evaluate_with_real_jacobian(&[1.0, 2.0])
            .unwrap_err();
        assert!(error.to_string().contains("imaginary part"));
    }

    #[test]
    fn derivative_request_is_rejected_for_eager_only_evaluator() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("sampling_test::x_only"));
        let mut evaluator = SamplingExpressionEvaluator::new([x.clone()], [x], false).unwrap();
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
            true,
        )
        .unwrap();
        let error = evaluator
            .evaluate_with_real_jacobian(&[1.0, 1.0])
            .unwrap_err();
        assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());

        // QuadFloat extends the mantissa using two binary64 numbers; its
        // exponent range cannot represent these coefficients either.
        let quad_one = F::<crate::utils::QuadFloat>::default().one().0;
        let error = evaluator
            .evaluate_with_real_jacobian(&[quad_one, quad_one])
            .unwrap_err();
        assert!(error.downcast_ref::<SamplingEvaluationError>().is_some());

        fn check<T: FloatLike>(evaluator: &mut SamplingExpressionEvaluator) {
            let one = F::<T>::default().one();
            let large = one.from_usize(10).powi(400);
            let small = &one / &large;
            let tolerance = one.epsilon() * one.from_usize(2048);
            // Reuse the same compiled program after its f64 failure.
            let evaluated = evaluator
                .evaluate_with_real_jacobian(&[one.0.clone(), one.0.clone()])
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
            SamplingExpressionEvaluator::new([x.clone() * x.clone()], [x], true).unwrap();
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
