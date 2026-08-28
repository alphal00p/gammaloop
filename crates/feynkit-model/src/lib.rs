//! Standalone particle-physics model types.
//!
//! This crate owns the validated, serializable model representation shared by
//! feynkit consumers. Import DTOs retain the strings used by UFO/JSON files,
//! while validated model records own parsed Symbolica expressions and typed
//! identifiers for every cross-reference.

#![forbid(unsafe_code)]

mod card;
mod error;
mod evaluation;
mod model;

pub use card::{ComplexValue, ParameterCard};
pub use error::{EntityKind, ModelError, ModelValidationError};
pub use evaluation::{
    EvaluatedValues, EvaluationFormFactor, EvaluationFunction, EvaluationRequest, ModelEvaluator,
    ModelExpression, RecomputeError,
};
pub use model::{
    Coupling, CouplingId, LorentzStructure, LorentzStructureId, Model, ModelFingerprint,
    ModelFormFactor, ModelFormFactorId, ModelFunction, ModelFunctionId, Order, OrderId, Parameter,
    ParameterId, ParameterNature, ParameterType, Particle, ParticleId, Propagator, PropagatorId,
    VertexRule, VertexRuleId,
};
