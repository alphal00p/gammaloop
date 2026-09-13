//! Legacy standalone integrands.
//!
//! These helpers predate process-level sampling acceptance and are retained
//! only while callers migrate. New validation should load a real amplitude or
//! cross-section process and use its reference-function acceptance path so
//! that graph routing, channel selection and Jacobians are exercised together.
pub mod h_function;
