#![cfg_attr(feature = "fail-on-warnings", deny(warnings))]
// #![deny(clippy::all)]
// #![warn(clippy::pedantic)]
#![warn(clippy::all)]
// #![warn(clippy::restriction)]
// #![warn(clippy::nursery)]
// #![warn(clippy::cargo)]
// #![feature(min_specialization)]
pub mod api;
pub mod cff;
pub mod cli_functions;
pub mod cross_section;
pub mod evaluation_result;
pub mod gammaloop_integrand;
pub mod graph;
pub mod h_function_test;
pub mod inspect;
pub mod integrands;
pub mod integrate;
pub mod linalg;
pub mod ltd;
pub mod model;
pub mod numerator;
pub mod observables;
pub mod tensor;
pub mod tests;
pub mod tests_from_pytest;
pub mod tropical;
pub mod utils;

use color_eyre::{Help, Report};
#[allow(unused)]
use colored::Colorize;
use eyre::WrapErr;

use integrands::*;
use lorentz_vector::LorentzVector;
use num::Complex;
use observables::ObservableSettings;
use observables::PhaseSpaceSelectorSettings;
use std::fs::File;
use std::sync::atomic::AtomicBool;
use utils::FloatLike;

use serde::{Deserialize, Serialize};

pub static INTERRUPTED: AtomicBool = AtomicBool::new(false);

pub const MAX_CORES: usize = 1000;

#[cfg(not(feature = "higher_loops"))]
pub const MAX_LOOP: usize = 3;
#[cfg(feature = "higher_loops")]
pub const MAX_LOOP: usize = 6;

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
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

pub fn set_interrupt_handler() {
    INTERRUPTED.store(false, std::sync::atomic::Ordering::Relaxed);
    let _ = ctrlc::set_handler(|| {
        INTERRUPTED.store(true, std::sync::atomic::Ordering::Relaxed);
    });
}

#[inline]
pub fn is_interrupted() -> bool {
    INTERRUPTED.load(std::sync::atomic::Ordering::Relaxed)
}

#[inline]
pub fn set_interrupted(flag: bool) {
    INTERRUPTED.store(flag, std::sync::atomic::Ordering::Relaxed);
}

const fn _default_true() -> bool {
    true
}
const fn _default_false() -> bool {
    false
}
const fn _default_usize_null() -> Option<usize> {
    None
}
fn _default_input_rescaling() -> Vec<Vec<(f64, f64)>> {
    vec![vec![(0.0, 1.0); 3]; 15]
}
fn _default_shifts() -> Vec<(f64, f64, f64, f64)> {
    vec![(1.0, 0.0, 0.0, 0.0); 15]
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct HFunctionSettings {
    pub function: HFunction,
    pub sigma: f64,
    #[serde(default = "_default_true")]
    pub enabled_dampening: bool,
    #[serde(default = "_default_usize_null")]
    pub power: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
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

#[derive(Debug, Clone, Deserialize, PartialEq, Default, Serialize)]
pub enum ParameterizationMapping {
    #[serde(rename = "log")]
    #[default]
    Log,
    #[serde(rename = "linear")]
    Linear,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct GeneralSettings {
    pub debug: usize,
    pub use_ltd: bool,
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize, Default, Serialize)]
pub enum IntegratedPhase {
    #[serde(rename = "real")]
    #[default]
    Real,
    #[serde(rename = "imag")]
    Imag,
    #[serde(rename = "both")]
    Both,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct KinematicsSettings {
    pub e_cm: f64,
    #[serde(default = "Externals::default")]
    pub externals: Externals,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct IntegratorSettings {
    pub n_bins: usize,
    pub bin_number_evolution: Option<Vec<usize>>,
    pub min_samples_for_update: usize,
    pub n_start: usize,
    pub n_increase: usize,
    pub n_max: usize,
    pub integrated_phase: IntegratedPhase,
    pub learning_rate: f64,
    pub train_on_avg: bool,
    pub show_max_wgt_info: bool,
    pub max_prob_ratio: f64,
    pub seed: u64,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ParameterizationSettings {
    pub mode: ParameterizationMode,
    pub mapping: ParameterizationMapping,
    pub b: f64,
    #[serde(default = "_default_input_rescaling")]
    pub input_rescaling: Vec<Vec<(f64, f64)>>,
    #[serde(default = "_default_shifts")]
    pub shifts: Vec<(f64, f64, f64, f64)>,
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

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
pub struct Settings {
    #[serde(rename = "General")]
    pub general: GeneralSettings,
    #[serde(rename = "Integrand")]
    pub hard_coded_integrand: IntegrandSettings,
    #[serde(rename = "Kinematics")]
    pub kinematics: KinematicsSettings,
    #[serde(rename = "Parameterization")]
    pub parameterization: ParameterizationSettings,
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
}

impl Settings {
    pub fn from_file(filename: &str) -> Result<Settings, Report> {
        let f = File::open(filename)
            .wrap_err_with(|| format!("Could not open settings file {}", filename))
            .suggestion("Does the path exist?")?;
        serde_yaml::from_reader(f)
            .wrap_err("Could not parse settings file")
            .suggestion("Is it a correct yaml file")
    }
}

#[derive(Serialize, Deserialize)]
pub struct IntegrationResult {
    pub neval: i64,
    pub fail: i32,
    pub result: Vec<f64>,
    pub error: Vec<f64>,
    pub prob: Vec<f64>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct StabilitySettings {
    rotation_axis: RotationMethod,
    levels: Vec<StabilityLevelSetting>,
}

impl Default for StabilitySettings {
    fn default() -> Self {
        Self {
            rotation_axis: RotationMethod::default(),
            levels: vec![
                StabilityLevelSetting::default_double(),
                StabilityLevelSetting::default_quad(),
            ],
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct StabilityLevelSetting {
    precision: Precision,
    required_precision_for_re: f64,
    required_precision_for_im: f64,
    escalate_for_large_weight_threshold: f64,
}

impl StabilityLevelSetting {
    fn default_double() -> Self {
        Self {
            precision: Precision::Double,
            required_precision_for_re: 1e-15,
            required_precision_for_im: 1e-15,
            escalate_for_large_weight_threshold: 0.9,
        }
    }

    fn default_quad() -> Self {
        Self {
            precision: Precision::Quad,
            required_precision_for_re: 1e-15,
            required_precision_for_im: 1e-15,
            escalate_for_large_weight_threshold: -1.0,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Default, Copy)]
pub enum RotationMethod {
    #[serde(rename = "x")]
    #[default]
    Pi2X,
    #[serde(rename = "y")]
    Pi2Y,
    #[serde(rename = "z")]
    Pi2Z,
}

impl RotationMethod {
    fn rotation_function<T: FloatLike>(&self) -> impl Fn(&LorentzVector<T>) -> LorentzVector<T> {
        match self {
            RotationMethod::Pi2X => utils::perform_pi2_rotation_x,
            RotationMethod::Pi2Y => utils::perform_pi2_rotation_y,
            RotationMethod::Pi2Z => utils::perform_pi2_rotation_z,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, Default, PartialEq, Copy)]
pub enum Precision {
    Single, // this might be useful for eventual deployment on gpu
    #[default]
    Double,
    Quad,
    Arb(usize),
}

#[derive(Serialize, Deserialize, Debug, Clone)]
#[serde(tag = "type", content = "momenta")]
pub enum Externals {
    #[serde(rename = "constant")]
    Constant(Vec<[f64; 4]>),
    // add different type of pdfs here when needed
}

impl Externals {
    #[allow(unused_variables)]
    #[inline]
    pub fn get_externals(&self, x_space_point: &[f64]) -> (Vec<LorentzVector<f64>>, f64) {
        match self {
            Externals::Constant(externals) => (
                externals
                    .iter()
                    .map(|[e0, e1, e2, e3]| LorentzVector::from_args(*e0, *e1, *e2, *e3))
                    .collect(),
                1.0,
            ),
        }
    }
}

impl Default for Externals {
    fn default() -> Self {
        Externals::Constant(vec![[0.0; 4]; 15])
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(tag = "type")]
pub enum SamplingSettings {
    #[default]
    #[serde(rename = "default")]
    Default,
    #[serde(rename = "multi_channeling")]
    MultiChanneling(MultiChannelingSettings),
    #[serde(rename = "discrete_graph_sampling")]
    DiscreteGraphs(DiscreteGraphSamplingSettings),
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
pub struct MultiChannelingSettings {
    pub alpha: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(tag = "subtype")]
pub enum DiscreteGraphSamplingSettings {
    #[default]
    #[serde(rename = "default")]
    Default,
    #[serde(rename = "multi_channeling")]
    MultiChanneling(MultiChannelingSettings),
    #[serde(rename = "discrete_multi_channeling")]
    DiscreteMultiChanneling(MultiChannelingSettings),
    #[serde(rename = "tropical")]
    TropicalSampling,
}
