use bincode_trait_derive::{Decode, Encode};
use global::GenerationSettings;
use serde::{Deserialize, Serialize};

use crate::{
    cli::tracing::LogLevel,
    integrands::IntegrandSettings,
    observables::{ObservableSettings, PhaseSpaceSelectorSettings},
    utils::serde_utils::IsDefault,
    GammaLoopContext,
};

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Default, Deserialize, Serialize, Encode, Decode)]
#[trait_decode(trait= GammaLoopContext)]
pub struct GlobalSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub debug_level: LogLevel,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub generation: GenerationSettings,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Default, Deserialize, Serialize, Encode, Decode)]
#[trait_decode(trait= GammaLoopContext)]
pub struct RuntimeSettings {
    // Runtime settings
    #[serde(
        rename = "General",
        default,
        skip_serializing_if = "IsDefault::is_default"
    )]
    pub general: GeneralSettings,
    #[serde(
        rename = "Integrand",
        default,
        skip_serializing_if = "IsDefault::is_default"
    )]
    pub hard_coded_integrand: IntegrandSettings,
    #[serde(
        rename = "Kinematics",
        default,
        skip_serializing_if = "IsDefault::is_default"
    )]
    pub kinematics: KinematicsSettings,
    #[serde(
        rename = "Integrator",
        default,
        skip_serializing_if = "IsDefault::is_default"
    )]
    pub integrator: IntegratorSettings,
    #[serde(
        rename = "Observables",
        default,
        skip_serializing_if = "IsDefault::is_default"
    )]
    pub observables: Vec<ObservableSettings>,
    #[serde(
        rename = "Selectors",
        default,
        skip_serializing_if = "IsDefault::is_default"
    )]
    pub selectors: Vec<PhaseSpaceSelectorSettings>,
    #[serde(rename = "Stability")]
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub stability: StabilitySettings,
    #[serde(
        rename = "sampling",
        default,
        skip_serializing_if = "IsDefault::is_default"
    )]
    pub sampling: SamplingSettings,
    #[serde(
        rename = "subtraction",
        default,
        skip_serializing_if = "IsDefault::is_default"
    )]
    pub subtraction: SubtractionSettings,
}

pub mod global;
pub use runtime::{
    kinematic::KinematicsSettings, GeneralSettings, IntegratorSettings, SamplingSettings,
    StabilitySettings, SubtractionSettings,
};
pub mod runtime;
