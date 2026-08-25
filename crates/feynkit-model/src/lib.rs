//! Standalone particle-physics model types.
//!
//! This crate owns the validated, serializable model representation shared by
//! feynkit consumers. Expressions intentionally remain strings here: parsing
//! them into a particular symbolic algebra belongs to downstream crates.

#![forbid(unsafe_code)]

mod card;
mod error;
mod evaluation;
mod model;

pub use card::{ComplexValue, ParameterCard};
pub use error::{EntityKind, ModelError, ModelValidationError};
pub use evaluation::{
    EvaluatedValues, EvaluationRequest, ModelEvaluator, ModelExpression, RecomputeError,
};
pub use model::{
    Coupling, LorentzStructure, Model, ModelDefinition, ModelFormFactor, ModelFunction, Order,
    Parameter, ParameterNature, ParameterType, Particle, ParticleId, Propagator, VertexRule,
    VertexRuleId,
};
