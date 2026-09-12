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

/// A warmup-compiled evaluator for one or more sampling-map expressions.
///
/// The eager program is always retained.  If `with_derivatives` is requested,
/// Symbolica's hyper-dual evaluator is built from the very same expression
/// tree, yielding exact first partials suitable for a forward-map Jacobian.
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
        let input = parameters
            .iter()
            .copied()
            .map(|value| Complex::new_re(F(value)))
            .collect::<Vec<_>>();
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
        let input = parameters
            .iter()
            .copied()
            .map(|value| Complex::new_re(F(value)))
            .collect::<Vec<_>>();
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
    }

    #[test]
    fn derivative_request_is_rejected_for_eager_only_evaluator() {
        test_initialise().unwrap();
        let x = Atom::var(symbol!("sampling_test::x_only"));
        let mut evaluator = SamplingExpressionEvaluator::new([x.clone()], [x], false).unwrap();
        assert!(evaluator.evaluate_with_derivatives(&[1.0]).is_err());
    }
}
