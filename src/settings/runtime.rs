use std::{fmt::Display, fs::File, path::Path};

use bincode_trait_derive::{Decode, Encode};
use color_eyre::{Result, Section};
use eyre::Context;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;

use crate::{
    momentum::{RotationMethod, ThreeMomentum},
    utils::{serde_utils::IsDefault, FloatLike, F},
    GammaLoopContext, HFunctionSettings,
};

use super::{global::OrientationPattern, RuntimeSettings};

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Encode, Decode)]
pub struct SubtractionSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub local_ct_settings: LocalCounterTermSettings,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub integrated_ct_settings: IntegratedCounterTermSettings,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub overlap_settings: OverlapSettings,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
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

impl RuntimeSettings {
    pub(crate) fn from_file(filename: impl AsRef<Path>) -> Result<RuntimeSettings> {
        let filename = filename.as_ref();
        let f = File::open(filename)
            .wrap_err_with(|| format!("Could not open settings file {}", filename.display()))
            .suggestion("Does the path exist?")?;
        serde_yaml::from_reader(f)
            .wrap_err("Could not parse settings file")
            .suggestion("Is it a correct yaml file")
    }
}

#[derive(Debug, Clone, Deserialize, Serialize, Encode, Decode, PartialEq)]
#[trait_decode(trait= GammaLoopContext)]
pub struct GeneralSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub debug: usize,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub use_ltd: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub orientation_pat: OrientationPattern,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub load_compiled_cff: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub load_compiled_numerator: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub joint_numerator_eval: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub amplitude_prefactor: Option<Complex<F<f64>>>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub load_compiled_separate_orientations: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub force_orientations: Option<Vec<usize>>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub m_uv: F<f64>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub mu_r_sq: F<f64>,
}
#[allow(clippy::derivable_impls)] // we might not want the standard defaults in the future
impl Default for GeneralSettings {
    fn default() -> Self {
        Self {
            debug: 0,
            use_ltd: false,
            load_compiled_numerator: true,
            joint_numerator_eval: true,
            load_compiled_cff: false,
            load_compiled_separate_orientations: false,
            amplitude_prefactor: Some(Complex::new(F(0.0), F(1.0))),
            force_orientations: None,
            orientation_pat: OrientationPattern::default(),
            m_uv: F(1000.0),
            mu_r_sq: F(1000.0 * 1000.0),
        }
    }
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize, Default, Serialize, Encode, Decode)]
// #[trait_decode(trait= GammaLoopContext)]
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

#[derive(Debug, Clone, Deserialize, Serialize, Encode, Decode, PartialEq)]
// #[trait_decode(trait= GammaLoopContext)]
pub struct IntegratorSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub n_bins: usize,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub bin_number_evolution: Option<Vec<usize>>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub min_samples_for_update: usize,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub n_start: usize,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub n_increase: usize,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub n_max: usize,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub integrated_phase: IntegratedPhase,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub discrete_dim_learning_rate: F<f64>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub continuous_dim_learning_rate: F<f64>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub train_on_avg: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub show_max_wgt_info: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub max_prob_ratio: F<f64>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
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
            integrated_phase: IntegratedPhase::Real,
            discrete_dim_learning_rate: F(1.5),
            continuous_dim_learning_rate: F(1.5),
            train_on_avg: false,
            show_max_wgt_info: true,
            max_prob_ratio: F(0.01),
            seed: 69,
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize, PartialEq, Encode, Decode)]
pub struct ParameterizationSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub mode: ParameterizationMode,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub mapping: ParameterizationMapping,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub b: f64,

    #[serde(
        default = "_default_input_rescaling",
        skip_serializing_if = "IsDefault::is_default"
    )]
    pub input_rescaling: Vec<Vec<(f64, f64)>>,
    #[serde(
        default = "_default_shifts",
        skip_serializing_if = "IsDefault::is_default"
    )]
    pub shifts: Vec<(f64, f64, f64, f64)>,
}

fn _default_input_rescaling() -> Vec<Vec<(f64, f64)>> {
    vec![vec![(0.0, 1.0); 3]; 15]
}
fn _default_shifts() -> Vec<(f64, f64, f64, f64)> {
    vec![(1.0, 0.0, 0.0, 0.0); 15]
}
impl Default for ParameterizationSettings {
    fn default() -> Self {
        Self {
            b: 1.0,
            mode: ParameterizationMode::Spherical,
            mapping: ParameterizationMapping::Linear,
            input_rescaling: vec![vec![(0.0, 1.0); 3]; 15],
            shifts: vec![(1.0, 0.0, 0.0, 0.0); 15],
        }
    }
}

#[derive(Serialize, Deserialize)]
pub struct IntegrationResult {
    pub neval: i64,
    pub fail: i32,
    pub result: Vec<F<f64>>,
    pub error: Vec<F<f64>>,
    pub prob: Vec<F<f64>>,
}

#[derive(Serialize, Deserialize, Debug, Clone, Encode, Decode, PartialEq)]
pub struct StabilitySettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub rotation_axis: Vec<RotationSetting>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub rotate_numerator: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub levels: Vec<StabilityLevelSetting>,
}

impl Default for StabilitySettings {
    fn default() -> Self {
        Self {
            rotation_axis: vec![RotationSetting::default()],
            levels: vec![
                StabilityLevelSetting::default_double(),
                StabilityLevelSetting::default_quad(),
            ],
            rotate_numerator: false,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy, Encode, Decode, PartialEq)]
pub struct StabilityLevelSetting {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub precision: Precision,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub required_precision_for_re: F<f64>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub required_precision_for_im: F<f64>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub escalate_for_large_weight_threshold: F<f64>,
}

impl StabilityLevelSetting {
    fn default_double() -> Self {
        Self {
            precision: Precision::Double,
            required_precision_for_re: F(1e-15),
            required_precision_for_im: F(1e-15),
            escalate_for_large_weight_threshold: F(0.9),
        }
    }

    fn default_quad() -> Self {
        Self {
            precision: Precision::Quad,
            required_precision_for_re: F(1e-5),
            required_precision_for_im: F(1e-5),
            escalate_for_large_weight_threshold: F(-1.0),
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Default, Copy, PartialEq, Encode, Decode)]
#[serde(tag = "type")]
pub enum RotationSetting {
    #[serde(rename = "x")]
    #[default]
    Pi2X,
    #[serde(rename = "y")]
    Pi2Y,
    #[serde(rename = "z")]
    Pi2Z,
    #[serde(rename = "none")]
    None,
    #[serde(rename = "euler_angles")]
    EulerAngles { alpha: f64, beta: f64, gamma: f64 },
}

impl RotationSetting {
    pub(crate) fn rotation_method(&self) -> RotationMethod {
        match self {
            Self::Pi2X => RotationMethod::Pi2X,
            Self::Pi2Y => RotationMethod::Pi2Y,
            Self::Pi2Z => RotationMethod::Pi2Z,
            Self::None => RotationMethod::Identity,
            Self::EulerAngles { alpha, beta, gamma } => {
                RotationMethod::EulerAngles(*alpha, *beta, *gamma)
            }
        }
    }

    #[allow(clippy::type_complexity)]
    pub(crate) fn rotation_function<'a, T: FloatLike + 'a>(
        &'a self,
    ) -> Box<dyn Fn(&'a ThreeMomentum<F<T>>) -> ThreeMomentum<F<T>> + 'a> {
        match self {
            Self::Pi2X => Box::new(ThreeMomentum::perform_pi2_rotation_x),
            Self::Pi2Y => Box::new(ThreeMomentum::perform_pi2_rotation_y),
            Self::Pi2Z => Box::new(ThreeMomentum::perform_pi2_rotation_z),
            Self::None => Box::new(|vector: &ThreeMomentum<F<T>>| vector.clone()),
            Self::EulerAngles { alpha, beta, gamma } => Box::new(|vector: &ThreeMomentum<F<T>>| {
                let mut cloned_vector = vector.clone();
                let alpha_t = F::<T>::from_f64(*alpha);
                let beta_t = F::<T>::from_f64(*beta);
                let gamma_t = F::<T>::from_f64(*gamma);
                cloned_vector.rotate_mut(&alpha_t, &beta_t, &gamma_t);
                cloned_vector
            }),
        }
    }

    pub(crate) fn as_str(&self) -> String {
        match self {
            Self::Pi2X => "x".to_owned(),
            Self::Pi2Y => "y".to_owned(),
            Self::Pi2Z => "z".to_owned(),
            Self::None => "none".to_owned(),
            Self::EulerAngles { alpha, beta, gamma } => {
                format!("euler {} {} {}", alpha, beta, gamma)
            }
        }
    }
}

#[derive(
    Serialize, Deserialize, Debug, Clone, Default, PartialEq, Copy, Hash, Eq, Encode, Decode,
)]
pub enum Precision {
    #[default]
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Encode, Decode)]
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OutputMetadata {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub model_name: String,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub output_type: String,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
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

#[derive(Debug, Clone, Deserialize, PartialEq, Default, Serialize, Encode, Decode)]
pub enum ParameterizationMapping {
    #[serde(rename = "log")]
    #[default]
    Log,
    #[serde(rename = "linear")]
    Linear,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode)]
#[serde(tag = "type")]
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode)]
pub struct MultiChannelingSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub alpha: f64,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode)]
pub struct GammaloopTropicalSamplingSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub upcast_on_failure: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub matrix_stability_test: Option<f64>,
}

impl Default for GammaloopTropicalSamplingSettings {
    fn default() -> Self {
        Self {
            upcast_on_failure: true,
            matrix_stability_test: Some(1.0e-5),
        }
    }
}

impl GammaloopTropicalSamplingSettings {
    pub(crate) fn into_tropical_sampling_settings(
        &self,
        debug: usize,
    ) -> momtrop::TropicalSamplingSettings {
        if self.upcast_on_failure {
            unimplemented!("upcast_on_failure removed from momtrop, implement automatic upcast of parameterization in gammaloop and then remove this crash")
        }

        momtrop::TropicalSamplingSettings {
            matrix_stability_test: self.matrix_stability_test,
            print_debug_info: debug > 0,
            return_metadata: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode)]
#[serde(tag = "subtype")]
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode)]
pub struct DiscreteGraphSamplingSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub sample_orientations: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub sampling_type: DiscreteGraphSamplingType,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Encode, Decode)]
pub struct LocalCounterTermSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub dampen_integrable_singularity: IntegrableSingularityDampener,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub uv_localisation: UVLocalisationSettings,
}

#[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq, Encode, Decode)]
#[serde(tag = "type")]
pub enum IntegrableSingularityDampener {
    #[serde(rename = "none")]
    None,
    #[default]
    #[serde(rename = "exponential")]
    Exponential,
    #[serde(rename = "powerlike")]
    Powerlike { power: f64 },
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode)]
pub struct UVLocalisationSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub sliver_width: f64,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub dynamic_width: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
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

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Encode, Decode)]
pub struct IntegratedCounterTermSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub range: IntegratedCounterTermRange,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode)]
#[serde(tag = "type")]
pub enum IntegratedCounterTermRange {
    #[serde(rename = "infinite")]
    Infinite {
        h_function_settings: HFunctionSettings,
    },
    #[serde(rename = "compact")]
    Compact,
}

impl Default for IntegratedCounterTermRange {
    fn default() -> Self {
        Self::Infinite {
            h_function_settings: HFunctionSettings::default(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Encode, Decode)]
pub struct OverlapSettings {
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub force_global_center: Option<Vec<[f64; 3]>>,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub check_global_center: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
    pub try_origin: bool,
    #[serde(default, skip_serializing_if = "IsDefault::is_default")]
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
