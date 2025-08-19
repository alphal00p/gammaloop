use bincode_trait_derive::{Decode, Encode};
use global::GenerationSettings;
use serde::{Deserialize, Serialize};

use crate::{
    integrands::IntegrandSettings,
    observables::{ObservableSettings, PhaseSpaceSelectorSettings},
    GammaLoopContext,
};

#[derive(Debug, Clone, Default, Deserialize, Serialize, Encode, Decode)]
#[trait_decode(trait= GammaLoopContext)]
pub struct GlobalSettings {
    pub generation: GenerationSettings,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize, Encode, Decode)]
#[trait_decode(trait= GammaLoopContext)]
pub struct RuntimeSettings {
    // Runtime settings
    #[serde(rename = "General")]
    pub general: GeneralSettings,
    #[serde(rename = "Integrand")]
    pub hard_coded_integrand: IntegrandSettings,
    #[serde(rename = "Kinematics")]
    pub kinematics: KinematicsSettings,
    #[serde(rename = "Integrator")]
    pub integrator: IntegratorSettings,
    #[serde(rename = "Observables")]
    pub observables: Vec<ObservableSettings>,
    #[serde(rename = "Selectors")]
    pub selectors: Vec<PhaseSpaceSelectorSettings>,
    #[serde(rename = "Stability")]
    #[serde(default = "StabilitySettings::default")]
    pub stability: StabilitySettings,
    #[serde(rename = "sampling")]
    pub sampling: SamplingSettings,
    #[serde(rename = "subtraction")]
    pub subtraction: SubtractionSettings,
}

pub mod global;
pub use runtime::{
    kinematic::KinematicsSettings, GeneralSettings, IntegratorSettings, SamplingSettings,
    StabilitySettings, SubtractionSettings,
};
pub mod runtime;
