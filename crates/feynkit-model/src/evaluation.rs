use std::{collections::BTreeMap, error::Error};

use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{ComplexValue, EntityKind, ModelError};

/// A callable UFO function supplied to an evaluator as a transport DTO.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct EvaluationFunction {
    pub name: String,
    pub arguments: Vec<String>,
    pub expression: Option<String>,
}

/// UFO form-factor metadata supplied to an evaluator as a transport DTO.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct EvaluationFormFactor {
    pub name: String,
    pub type_name: Option<String>,
    pub value: Option<String>,
}

/// A named expression whose value must be recomputed.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct ModelExpression {
    pub name: String,
    pub expression: String,
}

/// Owned input passed to a model evaluator.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct EvaluationRequest {
    /// External parameters and expressionless internal constants.
    pub known_parameters: BTreeMap<String, ComplexValue>,
    pub internal_parameters: Vec<ModelExpression>,
    pub couplings: Vec<ModelExpression>,
    /// Callable UFO functions available while parsing expressions.
    pub functions: Vec<EvaluationFunction>,
    /// Dynamic form-factor metadata referenced by Lorentz structures.
    pub form_factors: Vec<EvaluationFormFactor>,
}

/// Owned values returned by a model evaluator.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct EvaluatedValues {
    pub internal_parameters: BTreeMap<String, ComplexValue>,
    pub couplings: BTreeMap<String, ComplexValue>,
}

/// Expression-evaluation boundary for a [`crate::Model`].
///
/// Implementations may use Symbolica, another expression engine, or an
/// external service. The request and response are owned so an evaluator may
/// retain, transfer, or asynchronously process their data.
pub trait ModelEvaluator {
    type Error: Error + 'static;

    fn evaluate(&mut self, request: EvaluationRequest) -> Result<EvaluatedValues, Self::Error>;
}

impl<F, E> ModelEvaluator for F
where
    F: FnMut(EvaluationRequest) -> Result<EvaluatedValues, E>,
    E: Error + 'static,
{
    type Error = E;

    fn evaluate(&mut self, request: EvaluationRequest) -> Result<EvaluatedValues, Self::Error> {
        self(request)
    }
}

/// Errors returned while applying inputs and recomputing dependent values.
#[derive(Debug, Error)]
pub enum RecomputeError<E>
where
    E: Error + 'static,
{
    #[error(transparent)]
    Model(#[from] ModelError),
    #[error("model evaluator failed: {0}")]
    Evaluator(#[source] E),
    #[error("model evaluator did not return a value for {kind} '{name}'")]
    MissingValue { kind: EntityKind, name: String },
    #[error("model evaluator returned an unexpected value for {kind} '{name}'")]
    UnexpectedValue { kind: EntityKind, name: String },
}
