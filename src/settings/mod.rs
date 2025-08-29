use std::{fs::File, io::Read, path::Path};

use bincode_trait_derive::{Decode, Encode};
use color_eyre::{Result, Section};
use eyre::{eyre, Context};
use global::GenerationSettings;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::{
    cli::tracing::LogLevel,
    integrands::IntegrandSettings,
    observables::{ObservableSettings, PhaseSpaceSelectorSettings},
    utils::serde_utils::IsDefault,
    GammaLoopContext,
};

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Default, Deserialize, Serialize, Encode, Decode, JsonSchema)]
#[trait_decode(trait= GammaLoopContext)]
#[serde(deny_unknown_fields)]
pub struct GlobalSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub debug_level: LogLevel,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub generation: GenerationSettings,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Default, Deserialize, Serialize, Encode, Decode, JsonSchema)]
#[trait_decode(trait= GammaLoopContext)]
#[serde(deny_unknown_fields)]
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

impl RuntimeSettings {
    pub(crate) fn from_file(file_path: impl AsRef<Path>) -> Result<Self> {
        let mut f = File::open(file_path.as_ref())
            .wrap_err_with(|| {
                format!(
                    "Could not open runtime settings file {}",
                    file_path.as_ref().display()
                )
            })
            .suggestion("Does the path exist?")?;

        if let Some(ext) = file_path.as_ref().extension() {
            if let Some(ext) = ext.to_str() {
                match ext {
                    "json" => serde_json::from_reader(f)
                        .map_err(|e| eyre!(format!("Error parsing model json: {}", e)))
                        .suggestion("Is it a correct json file"),
                    "yaml" | "yml" => serde_yaml::from_reader(f)
                        .map_err(|e| eyre!(format!("Error parsing model yaml: {}", e)))
                        .suggestion("Is it a correct yaml file"),
                    "toml" => {
                        let mut buf = vec![];
                        f.read(&mut buf)?;
                        toml::from_slice(&buf)
                            .map_err(|e| eyre!(format!("Error parsing model toml: {}", e)))
                            .suggestion("Is it a correct toml file")
                    }

                    _ => Err(eyre!(format!("Unknown model file extension: {}", ext)))
                        .suggestion("Is it a .json or .yaml file?"),
                }
            } else {
                Err(eyre!(format!(
                    "Could not determine file extension of runtime settings file {}",
                    file_path.as_ref().display()
                )))
                .suggestion("Does the path exist?")
            }
        } else {
            Err(eyre!(format!(
                "Could not determine file extension of runtime settings file {}",
                file_path.as_ref().display()
            )))
            .suggestion("Does the path exist?")
        }
    }
}

pub mod global;
pub use runtime::{
    kinematic::KinematicsSettings, GeneralSettings, IntegratorSettings, SamplingSettings,
    StabilitySettings, SubtractionSettings,
};
pub mod runtime;

pub const fn _default_true() -> bool {
    true
}
pub const fn _default_false() -> bool {
    false
}
pub const fn _default_one() -> f64 {
    1.0
}
pub const fn _default_usize_null() -> Option<usize> {
    None
}
