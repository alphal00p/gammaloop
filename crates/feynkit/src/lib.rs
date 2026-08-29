//! The standalone FeynKit particle-physics toolbox.
//!
//! This crate contains no implementation logic. It provides a stable, feature-
//! gated facade over the focused `feynkit-*` crates, including the native
//! Spenso tensor reducer. Raw UFO loading is not a default feature because it
//! requires an attached Python interpreter.

#![forbid(unsafe_code)]

#[cfg(feature = "cff")]
pub use feynkit_cff as cff;
#[cfg(feature = "generator")]
pub use feynkit_generator as generator;
#[cfg(feature = "graph")]
pub use feynkit_graph as graph;
#[cfg(feature = "kinematics")]
pub use feynkit_kinematics as kinematics;
#[cfg(feature = "model")]
pub use feynkit_model as model;
#[cfg(feature = "tensor")]
pub use feynkit_tensor as tensor;
#[cfg(feature = "ufo")]
pub use feynkit_ufo as ufo;

#[cfg(feature = "cff")]
pub use feynkit_cff::{CffGenerator, CffOptions, CffResult};
#[cfg(feature = "generator")]
pub use feynkit_generator::{
    GenerationOptions, GenerationResult, Generator, ParticleSelector, Process,
};
#[cfg(feature = "graph")]
pub use feynkit_graph::{FeynmanDiagram, LoopMomentumBasis};
#[cfg(feature = "kinematics")]
pub use feynkit_kinematics::{
    Boost, FourMomentum, Helicity, JetDefinition, Rotation, ThreeMomentum,
};
#[cfg(feature = "model")]
pub use feynkit_model::{Model, ParameterCard};
#[cfg(feature = "tensor")]
pub use feynkit_tensor::{
    ContractionOrbit, CosetType, FeynmanDiagramTensorExt, OrthogonalWeingarten, TensorReducer,
    TensorReduction, TensorReductionError, TensorReductionTerm, WeingartenError,
};
#[cfg(feature = "ufo")]
pub use feynkit_ufo::{LoadedModel, UfoLoadOptions, UfoLoader};

#[cfg(all(
    test,
    feature = "cff",
    feature = "generator",
    feature = "graph",
    feature = "kinematics",
    feature = "model",
    feature = "tensor"
))]
mod tests {
    #[test]
    fn facade_exposes_primary_types() {
        let _ = std::any::TypeId::of::<super::Model>();
        let _ = std::any::TypeId::of::<super::FeynmanDiagram>();
        let _ = std::any::TypeId::of::<super::Generator>();
        let _ = std::any::TypeId::of::<super::CffGenerator>();
        let _ = std::any::TypeId::of::<super::FourMomentum<f64>>();
        let _ = std::any::TypeId::of::<super::TensorReducer>();
    }
}
