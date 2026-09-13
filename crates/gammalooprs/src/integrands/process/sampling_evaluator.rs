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
    utils::{F, hyperdual_utils::DualOrNot},
};
use symbolica::evaluate::OptimizationSettings;

/// Values and first partial derivatives of one evaluated expression.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingDualValue {
    /// The value of the expression at the supplied parameters.
    pub value: Complex<F<f64>>,
    /// First derivatives in the same order as the parameters passed to
    /// [`SamplingExpressionEvaluator::new`].
    pub derivatives: Vec<Complex<F<f64>>>,
}

/// Real-valued output and first-derivative matrix of a sampling expression.
///
/// Sampling maps are real maps even though the shared Symbolica evaluator uses
/// its complex domain.  This type performs that boundary conversion once and
/// computes the signed determinant of a square output/input Jacobian.  The
/// map density must use [`Self::absolute_determinant`], since a coordinate
/// chart may reverse orientation.
#[derive(Clone, Debug, PartialEq)]
pub struct SamplingJacobianEvaluation {
    /// Expression values in the order supplied to the evaluator.
    pub values: Vec<f64>,
    /// `derivatives[output][parameter]` from Symbolica's first-order duals.
    pub derivatives: Vec<Vec<f64>>,
    /// Signed determinant of `derivatives`.
    pub determinant: f64,
}

impl SamplingJacobianEvaluation {
    /// Absolute determinant required by a push-forward volume density.
    pub fn absolute_determinant(&self) -> f64 {
        self.determinant.abs()
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
    pub fn evaluate(&mut self, parameters: &[f64]) -> Result<Vec<Complex<F<f64>>>> {
        self.validate_parameters(parameters)?;
        let input = self.evaluator_parameters(parameters);
        let values = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut self.evaluator)(&input);
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
    pub fn evaluate_with_derivatives(
        &mut self,
        parameters: &[f64],
    ) -> Result<Vec<SamplingDualValue>> {
        if !self.with_derivatives {
            return Err(eyre!(
                "sampling evaluator was built without derivative support"
            ));
        }
        self.validate_parameters(parameters)?;
        let input = self.evaluator_parameters(parameters);
        let values = <f64 as GenericEvaluatorFloat>::get_evaluator(&mut self.evaluator)(&input);
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
    pub fn evaluate_with_real_jacobian(
        &mut self,
        parameters: &[f64],
    ) -> Result<SamplingJacobianEvaluation> {
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
            values.push(real_component(value.value, "value", output_index, None)?);
            let row = value
                .derivatives
                .into_iter()
                .enumerate()
                .map(|(parameter_index, derivative)| {
                    real_component(
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
        if !determinant.is_finite() {
            return Err(eyre!(
                "sampling Jacobian determinant is not finite: {determinant}"
            ));
        }
        Ok(SamplingJacobianEvaluation {
            values,
            derivatives,
            determinant,
        })
    }

    fn validate_parameters(&self, parameters: &[f64]) -> Result<()> {
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
    fn evaluator_parameters(&self, parameters: &[f64]) -> Vec<Complex<F<f64>>> {
        if !self.with_derivatives {
            return parameters
                .iter()
                .copied()
                .map(|value| Complex::new_re(F(value)))
                .collect();
        }
        let dual_width = self.parameter_count + 1;
        let mut input = Vec::with_capacity(self.parameter_count * dual_width);
        for (parameter_index, value) in parameters.iter().copied().enumerate() {
            input.push(Complex::new_re(F(value)));
            for derivative_index in 0..self.parameter_count {
                input.push(Complex::new_re(F(if parameter_index == derivative_index {
                    1.0
                } else {
                    0.0
                })));
            }
        }
        debug_assert_eq!(input.len(), self.parameter_count * dual_width);
        input
    }
}

fn real_component(
    value: Complex<F<f64>>,
    role: &str,
    output_index: usize,
    parameter_index: Option<usize>,
) -> Result<f64> {
    let imaginary = value.im.0;
    // Expressions in the map language are required to be real.  Permit only
    // roundoff-level imaginary residue from complex constants, and retain the
    // diagnostic location when a genuinely complex map is supplied.
    let tolerance = 1.0e-13 * value.re.0.abs().max(1.0);
    if !value.re.0.is_finite() || !imaginary.is_finite() {
        return Err(eyre!(
            "sampling {role} {output_index}{} is non-finite: {value:?}",
            parameter_index
                .map(|index| format!("/{index}"))
                .unwrap_or_default()
        ));
    }
    if imaginary.abs() > tolerance {
        return Err(eyre!(
            "sampling {role} {output_index}{} has non-negligible imaginary part {}",
            parameter_index
                .map(|index| format!("/{index}"))
                .unwrap_or_default(),
            imaginary
        ));
    }
    Ok(value.re.0)
}

fn signed_determinant(matrix: &[Vec<f64>]) -> Result<f64> {
    let dimension = matrix.len();
    if dimension == 0 || matrix.iter().any(|row| row.len() != dimension) {
        return Err(eyre!(
            "sampling Jacobian determinant requires a non-empty square matrix"
        ));
    }
    let mut work = matrix.to_vec();
    let mut determinant = 1.0;
    for pivot in 0..dimension {
        let (pivot_row, pivot_value) = (pivot..dimension)
            .map(|row| (row, work[row][pivot].abs()))
            .max_by(|(_, left), (_, right)| left.total_cmp(right))
            .expect("non-empty pivot range");
        if !pivot_value.is_finite() {
            return Err(eyre!("sampling Jacobian contains a non-finite pivot"));
        }
        if pivot_value == 0.0 {
            return Ok(0.0);
        }
        if pivot_row != pivot {
            work.swap(pivot_row, pivot);
            determinant = -determinant;
        }
        let diagonal = work[pivot][pivot];
        determinant *= diagonal;
        for row in (pivot + 1)..dimension {
            let factor = work[row][pivot] / diagonal;
            for column in (pivot + 1)..dimension {
                work[row][column] -= factor * work[pivot][column];
            }
        }
    }
    Ok(determinant)
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
}
