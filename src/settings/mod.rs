use bincode_trait_derive::{Decode, Encode};
use global::{GenerationSettings, Parallelisation};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::{
    GammaLoopContext,
    integrands::IntegrandSettings,
    observables::{ObservableSettings, PhaseSpaceSelectorSettings},
    settings::runtime::HFunctionSettings,
    utils::serde_utils::IsDefault,
};

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Deserialize, Serialize, Encode, Decode, JsonSchema, PartialEq)]
#[trait_decode(trait= GammaLoopContext)]
#[serde(default, deny_unknown_fields)]
pub struct GlobalSettings {
    #[serde(skip_serializing_if = "is_default_logfile_directive")]
    #[serde(default = "default_logfile_directive")]
    pub logfile_directive: String,
    #[serde(skip_serializing_if = "is_default_display_directive")]
    #[serde(default = "default_display_directive")]
    pub display_directive: String,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub generation: GenerationSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub n_cores: Parallelisation,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Default, Deserialize, Serialize, Encode, Decode, JsonSchema, PartialEq)]
#[trait_decode(trait= GammaLoopContext)]
#[serde(default, deny_unknown_fields)]
pub struct RuntimeSettings {
    // Runtime settings
    #[serde(rename = "general", skip_serializing_if = "IsDefault::is_default")]
    pub general: GeneralSettings,
    #[serde(rename = "integrand", skip_serializing_if = "IsDefault::is_default")]
    pub hard_coded_integrand: Option<IntegrandSettings>,
    #[serde(rename = "kinematics", skip_serializing_if = "IsDefault::is_default")]
    pub kinematics: KinematicsSettings,
    #[serde(rename = "integrator", skip_serializing_if = "IsDefault::is_default")]
    pub integrator: IntegratorSettings,
    #[serde(rename = "observables", skip_serializing_if = "IsDefault::is_default")]
    pub observables: Vec<ObservableSettings>,
    #[serde(rename = "selectors", skip_serializing_if = "IsDefault::is_default")]
    pub selectors: Vec<PhaseSpaceSelectorSettings>,
    #[serde(rename = "stability")]
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub stability: StabilitySettings,
    #[serde(rename = "sampling", skip_serializing_if = "IsDefault::is_default")]
    pub sampling: SamplingSettings,
    #[serde(rename = "subtraction", skip_serializing_if = "IsDefault::is_default")]
    pub subtraction: SubtractionSettings,
    #[serde(rename = "h_function", skip_serializing_if = "IsDefault::is_default")]
    pub lu_h_function: HFunctionSettings,
}

fn default_logfile_directive() -> String {
    "debug".to_string()
}

fn default_display_directive() -> String {
    "info".to_string()
}

fn is_default_logfile_directive(val: &String) -> bool {
    val == &default_logfile_directive()
}

fn is_default_display_directive(val: &String) -> bool {
    val == &default_display_directive()
}

impl Default for GlobalSettings {
    fn default() -> Self {
        Self {
            logfile_directive: default_logfile_directive(),
            display_directive: default_display_directive(),
            generation: GenerationSettings::default(),
            n_cores: Parallelisation::default(),
        }
    }
}

pub mod global;
pub use runtime::{
    GeneralSettings, IntegratorSettings, SamplingSettings, StabilitySettings, SubtractionSettings,
    kinematic::KinematicsSettings,
};
pub mod runtime;

#[cfg(test)]
mod tests {
    use serde::{Deserialize, Serialize};

    use crate::{
        improve_ps::PhaseSpaceImprovementSettings,
        momentum::{Dep, ExternalMomenta, SignOrZero},
        settings::{
            GlobalSettings, RuntimeSettings, SamplingSettings,
            global::GenerationSettings,
            runtime::{
                DiscreteGraphSamplingSettings, DiscreteGraphSamplingType,
                GammaloopTropicalSamplingSettings, kinematic::Externals,
            },
        },
        utils::{F, serde_utils::SHOWDEFAULTS},
    };
    use std::fmt::Debug;

    fn generic_test_settings<T>()
    where
        T: Serialize + for<'de> Deserialize<'de> + Default + PartialEq + Debug,
    {
        {
            let default = T::default();
            let serialized = serde_yaml::to_string(&default).unwrap();
            assert_eq!(serialized, "{}\n");
            let deserialized: T = serde_yaml::from_str(&serialized).unwrap();
            assert_eq!(default, deserialized);

            let deserialized_from_empty: T = serde_yaml::from_str("").unwrap();
            assert_eq!(default, deserialized_from_empty);
        }
        {
            let default = T::default();
            let serialized = toml::to_string_pretty(&default).unwrap();
            assert_eq!(serialized, "");
            let deserialized: T = toml::from_str(&serialized).unwrap();
            assert_eq!(default, deserialized);

            let deserialized_from_empty: T = toml::from_str("").unwrap();
            assert_eq!(default, deserialized_from_empty);
        }
        {
            let default = T::default();
            let serialized = serde_json::to_string(&default).unwrap();
            assert_eq!(serialized, "{}");
            let deserialized: T = serde_json::from_str(&serialized).unwrap();
            assert_eq!(default, deserialized);

            let deserialized_from_empty: T = serde_yaml::from_str("").unwrap();
            assert_eq!(default, deserialized_from_empty);
        }
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
    fn compile_test_serialize_deserialize() {
        use crate::settings::global::GammaloopCompileOptions;
        generic_test_settings::<GammaloopCompileOptions>();
    }

    #[test]
    fn tropical_subgraph_table_test_serialize_deserialize() {
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
    fn test_stability_level_settings_serialize_deserialize() {
        use crate::settings::runtime::StabilitySettings;
        generic_test_settings::<StabilitySettings>();
    }

    #[test]
    fn test_kinematics_settings_serialize_deserialize() {
        use crate::settings::KinematicsSettings;
        generic_test_settings::<KinematicsSettings>();

        let kinematics_settings = KinematicsSettings {
            e_cm: 100.0,
            externals: Externals::Constant {
                momenta: vec![
                    ExternalMomenta::Independent([F(1.), F(2.), F(3.), F(4.)]),
                    ExternalMomenta::Dependent(Dep::Dep),
                ],
                helicities: vec![SignOrZero::Plus, SignOrZero::Minus],
                improvement_settings: PhaseSpaceImprovementSettings::default(),
                f_64_cache: None,
                f_128_cache: None,
            },
        };

        let toml = toml::to_string_pretty(&kinematics_settings).unwrap();
        let deserialized: KinematicsSettings = toml::from_str(&toml).unwrap();
        assert_eq!(kinematics_settings, deserialized);
    }

    #[test]
    fn how_does_tropical_look() {
        SHOWDEFAULTS.store(true, std::sync::atomic::Ordering::Relaxed);
        let sampling_settings = SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
            sample_orientations: false,
            sampling_type: DiscreteGraphSamplingType::TropicalSampling(
                GammaloopTropicalSamplingSettings {
                    upcast_on_failure: false,
                    matrix_stability_test: Some(1e-5),
                },
            ),
        });
        let toml = toml::to_string_pretty(&sampling_settings).unwrap();
        println!("{toml}");
        SHOWDEFAULTS.store(false, std::sync::atomic::Ordering::Relaxed);
    }
}
