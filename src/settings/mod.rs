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
#[derive(Debug, Clone, Default, Deserialize, Serialize, Encode, Decode, JsonSchema, PartialEq)]
#[trait_decode(trait= GammaLoopContext)]
#[serde(default, deny_unknown_fields)]
pub struct GlobalSettings {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub debug_level: LogLevel,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub generation: GenerationSettings,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Default, Deserialize, Serialize, Encode, Decode, JsonSchema, PartialEq)]
#[trait_decode(trait= GammaLoopContext)]
#[serde(default, deny_unknown_fields)]
pub struct RuntimeSettings {
    // Runtime settings
    #[serde(rename = "General", skip_serializing_if = "IsDefault::is_default")]
    pub general: GeneralSettings,
    #[serde(rename = "Integrand", skip_serializing_if = "IsDefault::is_default")]
    pub hard_coded_integrand: Option<IntegrandSettings>,
    #[serde(rename = "Kinematics", skip_serializing_if = "IsDefault::is_default")]
    pub kinematics: KinematicsSettings,
    #[serde(rename = "Integrator", skip_serializing_if = "IsDefault::is_default")]
    pub integrator: IntegratorSettings,
    #[serde(rename = "Observables", skip_serializing_if = "IsDefault::is_default")]
    pub observables: Vec<ObservableSettings>,
    #[serde(rename = "Selectors", skip_serializing_if = "IsDefault::is_default")]
    pub selectors: Vec<PhaseSpaceSelectorSettings>,
    #[serde(rename = "Stability")]
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub stability: StabilitySettings,
    #[serde(rename = "sampling", skip_serializing_if = "IsDefault::is_default")]
    pub sampling: SamplingSettings,
    #[serde(rename = "subtraction", skip_serializing_if = "IsDefault::is_default")]
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

#[cfg(test)]
mod tests {
    use serde::{Deserialize, Serialize};

    use crate::settings::{global::GenerationSettings, GlobalSettings, RuntimeSettings};
    use std::fmt::Debug;

    fn generic_test_settings<T>()
    where
        T: Serialize + for<'de> Deserialize<'de> + Default + PartialEq + Debug,
    {
        let default = T::default();
        let serialized = serde_yaml::to_string(&default).unwrap();
        assert_eq!(serialized, "{}\n");
        let deserialized: T = serde_yaml::from_str(&serialized).unwrap();
        assert_eq!(default, deserialized);

        let deserialized_from_empty: T = serde_yaml::from_str("").unwrap();
        assert_eq!(default, deserialized_from_empty);
    }

    #[test]
    fn global_test_serialize_deserialize() {
        generic_test_settings::<GlobalSettings>();
    }

    #[test]
    fn generation_test_serialize_deserialize() {
        generic_test_settings::<GenerationSettings>();
    }

    #[test]
    fn gammaloop_compile_options_test_serialize_deserialize() {
        use crate::settings::global::GammaloopCompileOptions;
        generic_test_settings::<GammaloopCompileOptions>();
    }

    #[test]
    fn tropical_subgraph_table_settings_test_serialize_deserialize() {
        use crate::settings::global::TropicalSubgraphTableSettings;
        generic_test_settings::<TropicalSubgraphTableSettings>();
    }

    #[test]
    fn runtime_test_serialize_deserialize() {
        generic_test_settings::<RuntimeSettings>();
    }

    #[test]
    fn subtraction_settings_test_serialize_deserialize() {
        use crate::settings::runtime::SubtractionSettings;
        generic_test_settings::<SubtractionSettings>();
    }

    #[test]
    fn test_general_settings_serialize_deserialize() {
        use crate::settings::runtime::GeneralSettings;
        generic_test_settings::<GeneralSettings>();
    }

    #[test]
    fn test_integrator_settings_serialize_deserialize() {
        use crate::settings::runtime::IntegratorSettings;
        generic_test_settings::<IntegratorSettings>();
    }

    #[test]
    fn test_parameterization_settings_serialize_deserialize() {
        use crate::settings::runtime::ParameterizationSettings;
        generic_test_settings::<ParameterizationSettings>();
    }

    #[test]
    fn test_stability_settings_serialize_deserialize() {
        use crate::settings::runtime::StabilitySettings;
        generic_test_settings::<StabilitySettings>();
    }

    #[test]
    fn test_multi_channeling_settings_serialize_deserialize() {
        use crate::settings::runtime::MultiChannelingSettings;
        generic_test_settings::<MultiChannelingSettings>();
    }

    #[test]
    fn test_gammaloop_tropical_sampling_settings_serialize_deserialize() {
        use crate::settings::runtime::GammaloopTropicalSamplingSettings;
        generic_test_settings::<GammaloopTropicalSamplingSettings>();
    }

    #[test]
    fn test_discrete_graph_sampling_settings_serialize_deserialize() {
        use crate::settings::runtime::DiscreteGraphSamplingSettings;
        generic_test_settings::<DiscreteGraphSamplingSettings>();
    }

    #[test]
    fn test_local_counter_term_settings_serialize_deserialize() {
        use crate::settings::runtime::LocalCounterTermSettings;
        generic_test_settings::<LocalCounterTermSettings>();
    }

    #[test]
    fn test_uv_localisation_settings_serialize_deserialize() {
        use crate::settings::runtime::UVLocalisationSettings;
        generic_test_settings::<UVLocalisationSettings>();
    }

    #[test]
    fn test_integrated_counterterm_settings_serialize_deserialize() {
        use crate::settings::runtime::IntegratedCounterTermSettings;
        generic_test_settings::<IntegratedCounterTermSettings>();
    }

    #[test]
    fn test_overlap_settings_serialize_deserialize() {
        use crate::settings::runtime::OverlapSettings;
        generic_test_settings::<OverlapSettings>();
    }
    #[test]
    fn test_h_function_settings_serialize_deserialize() {
        use crate::settings::runtime::HFunctionSettings;
        generic_test_settings::<HFunctionSettings>();
    }
    #[test]
    fn test_kinematics_settings_serialize_deserialize() {
        use crate::settings::KinematicsSettings;
        generic_test_settings::<KinematicsSettings>();
    }
}
