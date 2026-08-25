//! Precision-generic tools for relativistic kinematics.
//!
//! This crate owns only language- and application-independent physics types:
//! three- and four-momenta, spatial rotations, Lorentz boosts, signed momentum
//! combinations, and generalized-kt jet clustering. Event records, particle
//! selection, phase-space sampling, and observables belong to callers.

#![forbid(unsafe_code)]

mod clustering;
mod helicity;
mod momentum;
mod signature;
mod transform;

pub use clustering::{ClusteringError, ClusteringResult, Jet, JetAlgorithm, JetDefinition};
pub use helicity::{Helicity, HelicityError};
pub use linnet::num_traits::{Sign, SignError, SignOrZero};
pub use momentum::{Energy, FourMomentum, ThreeMomentum};
pub use signature::{MomentumSignature, Signature, SignatureError};
pub use transform::{Axis, Boost, BoostError, Rotation};

use numerica::domains::float::{Real, SingleFloat};

/// Scalar operations required by precision-generic kinematics.
///
/// Numerica's real-number traits supply constants from a representative value,
/// so algorithms do not need to pass through `f64`. The blanket implementation
/// includes native floats, Numerica multiprecision values, and compatible
/// downstream precision wrappers.
pub trait KinematicScalar: Real + SingleFloat + PartialOrd {}

impl<T> KinematicScalar for T where T: Real + SingleFloat + PartialOrd {}
