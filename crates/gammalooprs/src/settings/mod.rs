use bincode_trait_derive::{Decode, Encode};
use global::{GenerationSettings, Parallelisation};
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::{
    GammaLoopContext,
    integrands::IntegrandSettings,
    observables::{ObservablesSettings, QuantitiesSettings, SelectorsSettings},
    settings::runtime::HFunctionSettings,
    utils::{
        F, FloatLike,
        serde_utils::{IsDefault, show_defaults_helper},
        tracing::LogStyle,
    },
};

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
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
    pub log_style: LogStyle,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub generation: GenerationSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub n_cores: Parallelisation,
}

#[cfg_attr(
    feature = "python_api",
    pyo3::pyclass(from_py_object, get_all, set_all)
)]
#[derive(Debug, Clone, Default, Deserialize, Serialize, Encode, Decode, JsonSchema, PartialEq)]
#[trait_decode(trait= GammaLoopContext)]
#[serde(default, deny_unknown_fields)]
pub struct RuntimeSettings {
    // Runtime settings
    #[serde(rename = "general", skip_serializing_if = "IsDefault::is_default")]
    pub general: GeneralSettings,
    #[serde(rename = "model", skip_serializing_if = "IsDefault::is_default")]
    pub model: RuntimeModelSettings,
    #[serde(rename = "integrand", skip_serializing_if = "IsDefault::is_default")]
    pub hard_coded_integrand: Option<IntegrandSettings>,
    #[serde(rename = "kinematics", skip_serializing_if = "IsDefault::is_default")]
    pub kinematics: KinematicsSettings,
    #[serde(rename = "integrator", skip_serializing_if = "IsDefault::is_default")]
    pub integrator: IntegratorSettings,
    #[serde(rename = "quantities", skip_serializing_if = "IsDefault::is_default")]
    pub quantities: QuantitiesSettings,
    #[serde(rename = "observables", skip_serializing_if = "IsDefault::is_default")]
    pub observables: ObservablesSettings,
    #[serde(rename = "selectors", skip_serializing_if = "IsDefault::is_default")]
    pub selectors: SelectorsSettings,
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

impl RuntimeSettings {
    pub(crate) fn additional_params<T: FloatLike>(&self) -> Vec<F<T>> {
        self.general
            .additional_param_values
            .iter()
            .map(|a| F(T::from_f64(*a)))
            .collect()
    }

    pub(crate) fn should_generate_events(&self) -> bool {
        self.general.generate_events
            || self.selectors.values().any(|selector| selector.active)
            || !self.observables.is_empty()
    }

    pub(crate) fn should_buffer_generated_events(&self) -> bool {
        self.general.generate_events || !self.observables.is_empty()
    }

    pub(crate) fn should_return_generated_events(&self) -> bool {
        self.general.generate_events
    }
}

fn default_logfile_directive() -> String {
    "off".to_string()
}

fn default_display_directive() -> String {
    "info".to_string()
}

fn is_default_logfile_directive(val: &String) -> bool {
    show_defaults_helper(val == &default_logfile_directive())
}

fn is_default_display_directive(val: &String) -> bool {
    show_defaults_helper(val == &default_display_directive())
}

impl Default for GlobalSettings {
    fn default() -> Self {
        Self {
            logfile_directive: default_logfile_directive(),
            display_directive: default_display_directive(),
            log_style: LogStyle::default(),
            generation: GenerationSettings::default(),
            n_cores: Parallelisation::default(),
        }
    }
}

pub mod global;
pub use runtime::{
    GeneralSettings, IntegratorSettings, ObservablesOutputSettings, RuntimeModelSettings,
    SamplingSettings, StabilitySettings, SubtractionSettings, kinematic::KinematicsSettings,
};
pub mod runtime;

#[cfg(test)]
mod tests {
    use serde::{Deserialize, Serialize};

    use crate::{
        momentum::{Dep, ExternalMomenta, Helicity},
        settings::{
            GlobalSettings, RuntimeSettings, SamplingSettings,
            global::{GammaloopCompileOptions, GenerationSettings},
            runtime::{
                DiscreteGraphSamplingSettings, DiscreteGraphSamplingType,
                GammaloopTropicalSamplingSettings,
                kinematic::{Externals, improvement::PhaseSpaceImprovementSettings},
            },
        },
        utils::{
            F,
            serde_utils::{SHOWDEFAULTS, ShowDefaultsGuard},
        },
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
    fn evaluator_settings_partial_tables_use_evaluator_defaults() {
        use crate::processes::EvaluatorSettings;

        let parsed: EvaluatorSettings =
            toml::from_str("horner_iterations = 0\ncpe_iterations = 0\n").unwrap();

        assert_eq!(
            parsed,
            EvaluatorSettings {
                horner_iterations: 0,
                cpe_iterations: Some(0),
                ..EvaluatorSettings::default()
            }
        );

        let serialized = toml::to_string(&parsed).unwrap();
        assert!(serialized.contains("horner_iterations = 0"));
        assert!(serialized.contains("cpe_iterations = 0"));
        assert!(!serialized.contains("iterative_orientation_optimization"));
        assert!(!serialized.contains("n_cores"));
        assert!(!serialized.contains("max_horner_scheme_variables"));
        assert!(!serialized.contains("max_common_pair_cache_entries"));
        assert!(!serialized.contains("max_common_pair_distance"));
    }

    #[test]
    fn compile_test_serialize_deserialize() {
        use crate::settings::global::GammaloopCompileOptions;
        generic_test_settings::<GammaloopCompileOptions>();
    }

    #[test]
    fn compile_settings_default_to_symjit_backend() {
        use crate::settings::global::{CompilationMode, GammaloopCompileOptions};

        assert_eq!(
            GammaloopCompileOptions::default().compilation_mode,
            CompilationMode::Symjit
        );
    }

    #[test]
    fn compile_settings_resolve_frozen_mode_from_compile_flag() {
        use crate::{
            processes::EvaluatorSettings,
            settings::global::{CompilationMode, FrozenCompilationMode, GammaloopCompileOptions},
        };

        let options = GammaloopCompileOptions {
            compilation_mode: CompilationMode::Assembly,
            ..Default::default()
        };

        assert_eq!(
            options.frozen_mode(&EvaluatorSettings {
                compile: false,
                ..Default::default()
            }),
            FrozenCompilationMode::Eager
        );
        assert!(matches!(
            options.frozen_mode(&EvaluatorSettings {
                compile: true,
                ..Default::default()
            }),
            FrozenCompilationMode::Assembly(_)
        ));
    }

    #[test]
    fn compile_settings_map_custom_args_into_symbolica_compile_options() {
        use crate::settings::global::{
            CompilationMode, CompilationOptimizationLevel, GammaloopCompileOptions,
        };

        let options = GammaloopCompileOptions {
            compilation_mode: CompilationMode::Cpp,
            optimization_level: CompilationOptimizationLevel::O1,
            fast_math: false,
            unsafe_math: false,
            compiler: "clang++".to_string(),
            custom: vec!["-g".to_string(), "-Winvalid".to_string()],
        };
        let symbolica = options.to_symbolica_compile_options();

        assert_eq!(symbolica.optimization_level, 1);
        assert!(!symbolica.fast_math);
        assert!(!symbolica.unsafe_math);
        assert_eq!(symbolica.compiler, "clang++");
        assert_eq!(symbolica.args, vec!["-g", "-Winvalid"]);
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
    fn runtime_model_settings_serialize_deserialize() {
        use crate::settings::runtime::RuntimeModelSettings;
        generic_test_settings::<RuntimeModelSettings>();
    }

    #[test]
    fn runtime_settings_serializes_model_overrides_under_model_block() {
        let mut settings = RuntimeSettings::default();
        settings
            .model
            .external_parameters
            .insert("mass_scalar_2".to_string(), (F(2.0), F(0.0)));

        let serialized = toml::to_string_pretty(&settings).unwrap();
        assert!(serialized.contains("[model]"));
        assert!(serialized.contains("mass_scalar_2 = ["));

        let deserialized: RuntimeSettings = toml::from_str(&serialized).unwrap();
        assert_eq!(settings, deserialized);
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
    fn runtime_event_generation_policy() {
        let mut settings = RuntimeSettings::default();
        assert!(!settings.should_generate_events());
        assert!(!settings.should_buffer_generated_events());
        assert!(!settings.should_return_generated_events());

        settings.observables.insert(
            "observable".to_string(),
            crate::observables::ObservableSettings {
                quantity: "pt".to_string(),
                selections: Vec::new(),
                entry_selection: crate::observables::EntrySelection::All,
                entry_index: 0,
                value_transform: crate::observables::ObservableValueTransform::Identity,
                phase: crate::observables::ObservablePhase::Real,
                misbinning_max_normalized_distance: None,
                histogram: crate::observables::HistogramSettings::Continuous(
                    crate::observables::ContinuousHistogramSettings {
                        x_min: 0.0,
                        x_max: 1.0,
                        n_bins: 1,
                        log_x_axis: false,
                        log_y_axis: true,
                        title: None,
                        type_description: "AL".to_string(),
                    },
                ),
            },
        );
        assert!(settings.should_generate_events());
        assert!(settings.should_buffer_generated_events());
        assert!(!settings.should_return_generated_events());

        settings.observables.clear();

        settings.selectors.insert(
            "selector".to_string(),
            crate::observables::SelectorSettings {
                quantity: "pt".to_string(),
                active: true,
                entry_selection: crate::observables::EntrySelection::All,
                entry_index: 0,
                selector: crate::observables::SelectorDefinitionSettings::CountRange(
                    crate::observables::CountRangeSelectorSettings {
                        min_count: 1,
                        max_count: None,
                    },
                ),
            },
        );
        assert!(settings.should_generate_events());
        assert!(!settings.should_buffer_generated_events());
        assert!(!settings.should_return_generated_events());

        settings.observables.insert(
            "observable".to_string(),
            crate::observables::ObservableSettings {
                quantity: "pt".to_string(),
                selections: Vec::new(),
                entry_selection: crate::observables::EntrySelection::All,
                entry_index: 0,
                value_transform: crate::observables::ObservableValueTransform::Identity,
                phase: crate::observables::ObservablePhase::Real,
                misbinning_max_normalized_distance: None,
                histogram: crate::observables::HistogramSettings::Continuous(
                    crate::observables::ContinuousHistogramSettings {
                        x_min: 0.0,
                        x_max: 1.0,
                        n_bins: 1,
                        log_x_axis: false,
                        log_y_axis: true,
                        title: None,
                        type_description: "AL".to_string(),
                    },
                ),
            },
        );
        assert!(settings.should_generate_events());
        assert!(settings.should_buffer_generated_events());
        assert!(!settings.should_return_generated_events());

        settings.general.generate_events = true;
        assert!(settings.should_generate_events());
        assert!(settings.should_buffer_generated_events());
        assert!(settings.should_return_generated_events());
    }

    #[test]
    fn test_integrator_settings_serialize_deserialize() {
        use crate::settings::runtime::IntegratorSettings;
        generic_test_settings::<IntegratorSettings>();
    }

    #[test]
    fn test_observables_output_settings_serialize_deserialize() {
        use crate::settings::runtime::ObservablesOutputSettings;
        generic_test_settings::<ObservablesOutputSettings>();
    }

    #[test]
    fn observables_output_settings_accept_single_entry_format_lists() {
        use crate::{
            observables::ObservableFileFormat, settings::runtime::ObservablesOutputSettings,
        };

        let parsed: ObservablesOutputSettings = toml::from_str("format = [\"hwu\"]").unwrap();
        assert_eq!(parsed.format, vec![ObservableFileFormat::Hwu]);

        let serialized = toml::to_string_pretty(&parsed).unwrap();
        assert_eq!(serialized, "format = [\"hwu\"]\n");
    }

    #[test]
    fn observables_output_settings_accept_multiple_formats() {
        use crate::{
            observables::ObservableFileFormat, settings::runtime::ObservablesOutputSettings,
        };

        let parsed: ObservablesOutputSettings =
            toml::from_str("format = [\"hwu\", \"json\"]").unwrap();
        assert_eq!(
            parsed.format,
            vec![ObservableFileFormat::Hwu, ObservableFileFormat::Json]
        );

        let serialized = toml::to_string_pretty(&parsed).unwrap();
        let reparsed: ObservablesOutputSettings = toml::from_str(&serialized).unwrap();
        assert_eq!(parsed, reparsed);
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
    fn test_uv_generation_settings_serialize_deserialize() {
        use crate::uv::UVgenerationSettings;
        generic_test_settings::<UVgenerationSettings>();
    }

    #[test]
    fn test_uv_generation_settings_serializes_renormalization_rules() {
        use crate::uv::{ApproximationType, CTIdentifier, UVgenerationSettings};
        use std::collections::{BTreeMap, BTreeSet};

        let settings = UVgenerationSettings {
            softct: true,
            renormalization_schemes: BTreeMap::from([
                (
                    CTIdentifier::new(BTreeSet::from([1]), Some(BTreeSet::from([1, 22]))),
                    ApproximationType::OS,
                ),
                (
                    CTIdentifier::new(BTreeSet::from([6]), None),
                    ApproximationType::Unsubtracted,
                ),
            ]),
            ..Default::default()
        };

        let toml = toml::to_string_pretty(&GenerationSettings {
            compile: GammaloopCompileOptions {
                compiler: "symjit".to_string(),
                ..Default::default()
            },
            uv: settings.clone(),
            ..Default::default()
        })
        .unwrap();
        println!("{}", toml);
        assert!(toml.contains("[[renormalization_schemes]]"));

        let deserialized_from_toml: UVgenerationSettings = toml::from_str(&toml).unwrap();
        assert_eq!(settings, deserialized_from_toml);

        let json = serde_json::to_string(&settings).unwrap();
        let deserialized_from_json: UVgenerationSettings = serde_json::from_str(&json).unwrap();
        assert_eq!(settings, deserialized_from_json);
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
                helicities: vec![Helicity::PLUS, Helicity::MINUS],
                improvement_settings: PhaseSpaceImprovementSettings::default(),
                f_64_cache: None,
                f_128_cache: None,
            },
        };

        let toml = toml::to_string_pretty(&kinematics_settings).unwrap();
        let deserialized: KinematicsSettings = toml::from_str(&toml).unwrap();
        assert_eq!(kinematics_settings, deserialized);

        let toml_without_ecm = format!(
            "{}\n",
            toml.lines()
                .filter(|line| !line.trim_start().starts_with("e_cm = "))
                .collect::<Vec<_>>()
                .join("\n")
        );
        let deserialized_without_ecm: KinematicsSettings =
            toml::from_str(&toml_without_ecm).unwrap();
        assert_eq!(
            deserialized_without_ecm.externals,
            kinematics_settings.externals
        );
        assert_eq!(deserialized_without_ecm.e_cm, 2.5);
    }

    #[test]
    fn sampling_settings_serializes_to_parser_shape() {
        let sampling_settings = SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
            sample_orientations: true,
            sampling_type: DiscreteGraphSamplingType::DiscreteMultiChanneling(
                crate::settings::runtime::MultiChannelingSettings::default(),
            ),
        });

        let _guard = ShowDefaultsGuard::new(true);
        let toml = toml::to_string_pretty(&sampling_settings).unwrap();
        assert!(toml.contains("graphs = \"monte_carlo\""));
        assert!(toml.contains("orientations = \"monte_carlo\""));
        assert!(toml.contains("lmb_multichanneling = true"));
        assert!(toml.contains("lmb_channels = \"monte_carlo\""));
        assert!(toml.contains("coordinate_system = \"spherical\""));
        assert!(!toml.contains("type = \"discrete_graph_sampling\""));
        assert!(!toml.contains("subtype = \"discrete_multi_channeling\""));
    }

    #[test]
    fn sampling_settings_rejects_incompatible_orientation_sampling() {
        let invalid_toml = r#"
graphs = "summed"
orientations = "monte_carlo"
lmb_multichanneling = false
lmb_channels = "summed"
coordinate_system = "spherical"
mapping = "linear"
b = 1.0
"#;

        let err = toml::from_str::<SamplingSettings>(invalid_toml).unwrap_err();
        assert!(err.to_string().contains(
            "orientations can only be set to 'monte_carlo' when graphs is 'monte_carlo'"
        ));
    }

    #[test]
    fn sampling_settings_deserializes_from_parser_shape() {
        let toml = r#"
graphs = "monte_carlo"
orientations = "summed"
lmb_multichanneling = true
lmb_channels = "summed"
coordinate_system = "momentum_space"
mapping = "log"
b = 5.0
"#;

        let settings: SamplingSettings = toml::from_str(toml).unwrap();
        assert_eq!(
            settings,
            SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
                sample_orientations: false,
                sampling_type: DiscreteGraphSamplingType::MultiChanneling(
                    crate::settings::runtime::MultiChannelingSettings {
                        alpha: 3.0,
                        parameterization_settings:
                            crate::settings::runtime::ParameterizationSettings {
                                mode: crate::settings::runtime::ParameterizationMode::MomentumSpace,
                                mapping: crate::settings::runtime::ParameterizationMapping::Log,
                                b: 5.0,
                            },
                    },
                ),
            })
        );
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
