use std::fmt::Display;

use bincode_trait_derive::{Decode, Encode};
use eyre::Result;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use typed_index_collections::TiVec;

use crate::{
    improve_ps::generate_default_momenta,
    momentum::RotationMethod,
    momentum_sample::ExternalIndex,
    settings::runtime::kinematic::Externals,
    signature::SignatureLike,
    utils::{
        serde_utils::{
            is_false, is_float, is_true, is_u64, is_usize, IsDefault, _default_rotation_axis,
            _default_stability_levels, is_default_rotation_axis, is_default_stability_levels,
        },
        F,
    },
    GammaLoopContext,
};

use super::{global::OrientationPattern, RuntimeSettings};

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Encode, Decode, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct SubtractionSettings {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub local_ct_settings: LocalCounterTermSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub integrated_ct_settings: IntegratedCounterTermSettings,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub overlap_settings: OverlapSettings,
    #[serde(skip_serializing_if = "is_false")]
    pub disable_threshold_subtraction: bool,
}

#[derive(Copy, Clone)]
pub struct LockedRuntimeSettings<'a>(&'a RuntimeSettings);
impl<'a> From<&'a RuntimeSettings> for LockedRuntimeSettings<'a> {
    fn from(value: &'a RuntimeSettings) -> Self {
        LockedRuntimeSettings(value)
    }
}

impl<'a> From<LockedRuntimeSettings<'a>> for RuntimeSettings {
    fn from(value: LockedRuntimeSettings) -> Self {
        value.0.clone()
    }
}

impl<'a> LockedRuntimeSettings<'a> {
    /// If no external momenta are set in the settings, generate a default kinematic point that satisfies the phase-space constraints
    pub(crate) fn into_with_modified_kinematics(
        self,
        external_signature: &SignatureLike<ExternalIndex>,
        external_masses: &TiVec<ExternalIndex, F<f64>>,
    ) -> Result<RuntimeSettings> {
        if external_signature.is_empty() {
            Ok(self.into())
        } else {
            match &self.0.kinematics.externals {
                Externals::Constant {
                    momenta,
                    helicities,
                    ..
                } => {
                    if momenta.is_empty() && helicities.is_empty() {
                        let new_externals = generate_default_momenta(
                            external_masses,
                            external_signature,
                            &F(self.0.kinematics.e_cm),
                        )?;

                        let mut new_settigns = self.0.clone();
                        new_settigns.kinematics.externals = new_externals;
                        Ok(new_settigns)
                    } else {
                        Ok(self.into())
                    }
                }
            }
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Deserialize, Serialize, Encode, Decode, PartialEq, JsonSchema)]
#[trait_decode(trait= GammaLoopContext)]
#[serde(default, deny_unknown_fields)]
pub struct GeneralSettings {
    #[serde(skip_serializing_if = "is_false")]
    pub use_ltd: bool,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub orientation_pat: OrientationPattern,
    #[serde(skip_serializing_if = "is_false")]
    pub load_compiled_cff: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub enable_cache: bool,
    #[serde(skip_serializing_if = "is_float::<1000>")]
    pub m_uv: f64,
    #[serde(skip_serializing_if = "is_float::<1000_000>")]
    pub mu_r_sq: f64,
}

impl Default for GeneralSettings {
    fn default() -> Self {
        Self {
            use_ltd: false,
            load_compiled_cff: false,
            enable_cache: false,
            orientation_pat: OrientationPattern::default(),
            m_uv: 1000.0,
            mu_r_sq: 1000.0 * 1000.0,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(
    Debug, Copy, Clone, PartialEq, Deserialize, Default, Serialize, Encode, Decode, JsonSchema,
)]
// #[trait_decode(trait= GammaLoopContext)]
#[serde(deny_unknown_fields)]
pub enum IntegratedPhase {
    #[serde(rename = "real")]
    #[default]
    Real,
    #[serde(rename = "imag")]
    Imag,
    #[serde(rename = "both")]
    Both,
}

pub mod kinematic;

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Deserialize, Serialize, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct IntegratorSettings {
    #[serde(skip_serializing_if = "is_usize::<64>")]
    pub n_bins: usize,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub bin_number_evolution: Option<Vec<usize>>,
    #[serde(skip_serializing_if = "is_usize::<1000>")]
    pub min_samples_for_update: usize,
    #[serde(skip_serializing_if = "is_usize::<100000>")]
    pub n_start: usize,
    #[serde(skip_serializing_if = "is_usize::<10000>")]
    pub n_increase: usize,
    #[serde(skip_serializing_if = "is_usize::<10000000000>")]
    pub n_max: usize,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub integrated_phase: IntegratedPhase,
    #[serde(skip_serializing_if = "is_float::<1>")]
    pub discrete_dim_learning_rate: f64,
    #[serde(skip_serializing_if = "is_float::<1>")]
    pub continuous_dim_learning_rate: f64,
    #[serde(skip_serializing_if = "is_false")]
    pub train_on_avg: bool,
    #[serde(skip_serializing_if = "is_true")]
    pub show_max_wgt_info: bool,
    #[serde(skip_serializing_if = "is_float::<1000>")]
    pub max_prob_ratio: f64,
    #[serde(skip_serializing_if = "is_u64::<69>")]
    pub seed: u64,
}

impl Default for IntegratorSettings {
    fn default() -> Self {
        Self {
            n_bins: 64,
            bin_number_evolution: None,
            min_samples_for_update: 1000,
            n_start: 100000,
            n_increase: 10000,
            n_max: 10000000000,
            integrated_phase: IntegratedPhase::default(),
            discrete_dim_learning_rate: 1.0,
            continuous_dim_learning_rate: 1.0,
            train_on_avg: false,
            show_max_wgt_info: true,
            max_prob_ratio: 1000.0,
            seed: 69,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Deserialize, Serialize, PartialEq, Encode, Decode, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct ParameterizationSettings {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub mode: ParameterizationMode,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub mapping: ParameterizationMapping,
    #[serde(skip_serializing_if = "is_float::<1>")]
    pub b: f64,
}

impl Default for ParameterizationSettings {
    fn default() -> Self {
        Self {
            b: 1.0,
            mode: ParameterizationMode::default(),
            mapping: ParameterizationMapping::default(),
        }
    }
}

#[derive(Serialize, Deserialize, Default, Debug)]
#[serde(default, deny_unknown_fields)]
pub struct IntegrationResult {
    pub neval: usize,
    pub real_zero: usize,
    pub im_zero: usize,
    pub result: Complex<F<f64>>,
    pub error: Complex<F<f64>>,
    pub real_chisq: F<f64>,
    pub im_chisq: F<f64>,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Serialize, Deserialize, Debug, Clone, Encode, Decode, PartialEq, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct StabilitySettings {
    #[serde(skip_serializing_if = "is_default_rotation_axis")]
    pub rotation_axis: Vec<RotationSetting>,
    #[serde(skip_serializing_if = "is_false")]
    pub rotate_numerator: bool,
    #[serde(skip_serializing_if = "is_default_stability_levels")]
    pub levels: Vec<StabilityLevelSetting>,
    #[serde(skip_serializing_if = "is_false")]
    pub check_on_norm: bool,
}

impl Default for StabilitySettings {
    fn default() -> Self {
        Self {
            rotation_axis: _default_rotation_axis(),
            levels: _default_stability_levels(),
            rotate_numerator: false,
            check_on_norm: false,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Serialize, Deserialize, Debug, Clone, Copy, Encode, Decode, PartialEq, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct StabilityLevelSetting {
    pub precision: Precision,
    pub required_precision_for_re: f64,
    pub required_precision_for_im: f64,
    pub escalate_for_large_weight_threshold: f64,
}

impl StabilityLevelSetting {
    pub fn default_double() -> Self {
        Self {
            precision: Precision::Double,
            required_precision_for_re: 1e-10,
            required_precision_for_im: 1e-10,
            escalate_for_large_weight_threshold: 0.9,
        }
    }

    pub fn default_quad() -> Self {
        Self {
            precision: Precision::Quad,
            required_precision_for_re: 1e-5,
            required_precision_for_im: 1e-5,
            escalate_for_large_weight_threshold: -1.0,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq, Encode, Decode, JsonSchema)]
#[serde(tag = "type")]
#[serde(deny_unknown_fields)]
pub enum RotationSetting {
    #[serde(rename = "x")]
    Pi2X {},
    #[serde(rename = "y")]
    Pi2Y {},
    #[serde(rename = "z")]
    Pi2Z {},
    #[serde(rename = "none")]
    None {},
    #[serde(rename = "euler_angles")]
    EulerAngles { alpha: f64, beta: f64, gamma: f64 },
}

impl Default for RotationSetting {
    fn default() -> Self {
        Self::Pi2Z {}
    }
}

impl RotationSetting {
    pub(crate) fn rotation_method(&self) -> RotationMethod {
        match self {
            Self::Pi2X {} => RotationMethod::Pi2X,
            Self::Pi2Y {} => RotationMethod::Pi2Y,
            Self::Pi2Z {} => RotationMethod::Pi2Z,
            Self::None {} => RotationMethod::Identity,
            Self::EulerAngles { alpha, beta, gamma } => {
                RotationMethod::EulerAngles(*alpha, *beta, *gamma)
            }
        }
    }

    pub(crate) fn _as_str(&self) -> String {
        match self {
            Self::Pi2X {} => "x".to_owned(),
            Self::Pi2Y {} => "y".to_owned(),
            Self::Pi2Z {} => "z".to_owned(),
            Self::None {} => "none".to_owned(),
            Self::EulerAngles { alpha, beta, gamma } => {
                format!("euler {} {} {}", alpha, beta, gamma)
            }
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(
    Serialize, Deserialize, Debug, Clone, PartialEq, Copy, Hash, Eq, Encode, Decode, JsonSchema,
)]
#[serde(deny_unknown_fields)]
pub enum Precision {
    Double,
    Quad,
    Arb,
}

impl Display for Precision {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Double => write!(f, "f64"),
            Self::Quad => write!(f, "f128"),
            Self::Arb => write!(f, "arb"),
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Encode, Decode, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum ParameterizationMode {
    #[serde(rename = "cartesian")]
    Cartesian,
    #[serde(rename = "spherical")]
    #[default]
    Spherical,
    #[serde(rename = "hyperspherical")]
    HyperSpherical,
    #[serde(rename = "hyperspherical_flat")]
    HyperSphericalFlat,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct OutputMetadata {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub model_name: String,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub output_type: String,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub contents: Vec<String>,
}

impl Display for ParameterizationMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParameterizationMode::Cartesian => write!(f, "cartesian"),
            ParameterizationMode::Spherical => write!(f, "spherical"),
            ParameterizationMode::HyperSpherical => write!(f, "hyperspherical"),
            ParameterizationMode::HyperSphericalFlat => write!(f, "flat hyperspherical"),
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(Debug, Clone, Deserialize, PartialEq, Default, Serialize, Encode, Decode, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum ParameterizationMapping {
    #[serde(rename = "log")]
    Log,
    #[serde(rename = "linear")]
    #[default]
    Linear,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
#[serde(tag = "type")]
#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[serde(deny_unknown_fields)]
pub enum SamplingSettings {
    #[serde(rename = "default")]
    Default(ParameterizationSettings),
    #[serde(rename = "multi_channeling")]
    MultiChanneling(MultiChannelingSettings),
    #[serde(rename = "discrete_graph_sampling")]
    DiscreteGraphs(DiscreteGraphSamplingSettings),
}

impl Default for SamplingSettings {
    fn default() -> Self {
        Self::Default(ParameterizationSettings::default())
    }
}

impl SamplingSettings {
    pub(crate) fn get_parameterization_settings(&self) -> Option<ParameterizationSettings> {
        match self {
            SamplingSettings::Default(settings) => Some(settings.clone()),
            SamplingSettings::MultiChanneling(settings) => {
                Some(settings.parameterization_settings.clone())
            }
            SamplingSettings::DiscreteGraphs(settings) => match &settings.sampling_type {
                DiscreteGraphSamplingType::Default(settings) => Some(settings.clone()),
                DiscreteGraphSamplingType::MultiChanneling(settings) => {
                    Some(settings.parameterization_settings.clone())
                }
                DiscreteGraphSamplingType::DiscreteMultiChanneling(settings) => {
                    Some(settings.parameterization_settings.clone())
                }
                DiscreteGraphSamplingType::TropicalSampling(_) => None,
            },
        }
    }

    pub(crate) fn sample_orientations(&self) -> bool {
        match self {
            SamplingSettings::Default(_) => false,
            SamplingSettings::MultiChanneling(_) => false,
            SamplingSettings::DiscreteGraphs(settings) => settings.sample_orientations,
        }
    }

    pub(crate) fn discrete_depth(&self) -> usize {
        match self {
            SamplingSettings::Default(_) => 0,
            SamplingSettings::MultiChanneling(_) => 0,
            SamplingSettings::DiscreteGraphs(settings) => {
                let depth_from_orientations = settings.sample_orientations as usize;

                match &settings.sampling_type {
                    DiscreteGraphSamplingType::Default(_) => 1 + depth_from_orientations,
                    DiscreteGraphSamplingType::MultiChanneling(_) => 1 + depth_from_orientations,
                    DiscreteGraphSamplingType::DiscreteMultiChanneling(_) => {
                        2 + depth_from_orientations
                    }
                    DiscreteGraphSamplingType::TropicalSampling(_) => 1 + depth_from_orientations,
                }
            }
        }
    }

    pub(crate) fn describe_settings(&self) -> String {
        match self {
            SamplingSettings::Default(settings) => {
                format!("{} coordinates", settings.mode)
            }
            SamplingSettings::MultiChanneling(settings) => {
                format!(
                    "lmb multichanneling in {} coordinates",
                    settings.parameterization_settings.mode
                )
            }
            SamplingSettings::DiscreteGraphs(settings) => {
                let discrete_graph_string = "Monte Carlo over graphs";
                let orientation_sampling_string = if settings.sample_orientations {
                    "and Monte Carlo over orientations"
                } else {
                    ""
                };

                match &settings.sampling_type {
                    DiscreteGraphSamplingType::Default(settings) => {
                        format!(
                            "{} {} in {} coordinates",
                            discrete_graph_string, orientation_sampling_string, settings.mode
                        )
                    }
                    DiscreteGraphSamplingType::MultiChanneling(settings) => {
                        format!(
                            "{}, lmb multichanneling in {} coordinates {}",
                            discrete_graph_string,
                            settings.parameterization_settings.mode,
                            orientation_sampling_string,
                        )
                    }
                    DiscreteGraphSamplingType::DiscreteMultiChanneling(settings) => {
                        format!(
                            "{}, {} and monte carlo over lmbs in {} coordinates",
                            discrete_graph_string,
                            orientation_sampling_string,
                            settings.parameterization_settings.mode
                        )
                    }
                    DiscreteGraphSamplingType::TropicalSampling(_) => {
                        format!(
                            "{} {} using 🌴🥥 tropical sampling 🥥🌴",
                            discrete_graph_string, orientation_sampling_string,
                        )
                    }
                }
            }
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
#[serde(deny_unknown_fields, default)]
pub struct MultiChannelingSettings {
    #[serde(skip_serializing_if = "is_float::<3>")]
    pub alpha: f64,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub parameterization_settings: ParameterizationSettings,
}

impl Default for MultiChannelingSettings {
    fn default() -> Self {
        Self {
            alpha: 3.0,
            parameterization_settings: ParameterizationSettings::default(),
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
#[serde(deny_unknown_fields, default)]
pub struct GammaloopTropicalSamplingSettings {
    #[serde(skip_serializing_if = "is_true")]
    pub upcast_on_failure: bool,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub matrix_stability_test: Option<f64>,
}

impl Default for GammaloopTropicalSamplingSettings {
    fn default() -> Self {
        Self {
            upcast_on_failure: true,
            matrix_stability_test: None,
        }
    }
}

impl GammaloopTropicalSamplingSettings {
    pub fn into_tropical_sampling_settings(&self) -> momtrop::TropicalSamplingSettings {
        momtrop::TropicalSamplingSettings {
            matrix_stability_test: self.matrix_stability_test,
            print_debug_info: false,
            return_metadata: false,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
#[serde(tag = "subtype", deny_unknown_fields)]
pub enum DiscreteGraphSamplingType {
    #[serde(rename = "default")]
    Default(ParameterizationSettings),
    #[serde(rename = "multi_channeling")]
    MultiChanneling(MultiChannelingSettings),
    #[serde(rename = "discrete_multi_channeling")]
    DiscreteMultiChanneling(MultiChannelingSettings),
    #[serde(rename = "tropical")]
    TropicalSampling(GammaloopTropicalSamplingSettings),
}

impl Default for DiscreteGraphSamplingType {
    fn default() -> Self {
        DiscreteGraphSamplingType::Default(ParameterizationSettings::default())
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema, Default)]
#[serde(deny_unknown_fields, default)]
pub struct DiscreteGraphSamplingSettings {
    #[serde(skip_serializing_if = "is_false")]
    pub sample_orientations: bool,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub sampling_type: DiscreteGraphSamplingType,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
#[serde(deny_unknown_fields, default)]
pub struct LocalCounterTermSettings {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub uv_localisation: UVLocalisationSettings,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct UVLocalisationSettings {
    #[serde(skip_serializing_if = "is_float::<10>")]
    pub sliver_width: f64,
    #[serde(skip_serializing_if = "is_false")]
    pub dynamic_width: bool,
    #[serde(skip_serializing_if = "is_float::<1>")]
    pub gaussian_width: f64,
}

impl Default for UVLocalisationSettings {
    fn default() -> Self {
        Self {
            sliver_width: 10.0,
            dynamic_width: false,
            gaussian_width: 1.0,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
#[serde(deny_unknown_fields, default)]
pub struct IntegratedCounterTermSettings {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub range: IntegratedCounterTermRange,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
#[serde(tag = "type")]
#[serde(deny_unknown_fields)]
pub enum IntegratedCounterTermRange {
    #[serde(rename = "infinite")]
    Infinite {
        h_function_settings: HFunctionSettings,
    },
    #[serde(rename = "compact")]
    Compact {},
}

impl Default for IntegratedCounterTermRange {
    fn default() -> Self {
        Self::Infinite {
            h_function_settings: HFunctionSettings::default(),
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(get_all, set_all))]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
#[serde(deny_unknown_fields)]
#[serde(default)]
pub struct OverlapSettings {
    //v
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub force_global_center: Option<Vec<[f64; 3]>>,
    #[serde(skip_serializing_if = "is_true")]
    pub check_global_center: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub try_origin: bool,
    #[serde(skip_serializing_if = "is_false")]
    pub try_origin_all_lmbs: bool,
}

impl Default for OverlapSettings {
    fn default() -> Self {
        Self {
            force_global_center: None,
            check_global_center: true,
            try_origin: false,
            try_origin_all_lmbs: false,
        }
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
#[serde(deny_unknown_fields)]
// #[trait_decode(trait= GammaLoopContext)]
pub enum HFunction {
    #[default]
    #[serde(rename = "poly_exponential")]
    PolyExponential,
    #[serde(rename = "exponential")]
    Exponential,
    #[serde(rename = "poly_left_right_exponential")]
    PolyLeftRightExponential,
    #[serde(rename = "exponential_ct")]
    ExponentialCT,
}

#[cfg_attr(feature = "python_api", pyo3::pyclass)]
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode, JsonSchema)]
#[serde(deny_unknown_fields)]
#[serde(default)]
pub struct HFunctionSettings {
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub function: HFunction,
    #[serde(skip_serializing_if = "is_float::<1>")]
    pub sigma: f64,
    #[serde(skip_serializing_if = "is_true")]
    pub enabled_dampening: bool,
    #[serde(skip_serializing_if = "IsDefault::is_default")]
    pub power: Option<usize>,
}

impl Default for HFunctionSettings {
    fn default() -> Self {
        Self {
            sigma: 1.0,
            function: HFunction::default(),
            enabled_dampening: true,
            power: None,
        }
    }
}
