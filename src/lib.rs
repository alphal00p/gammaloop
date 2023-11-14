#![cfg_attr(feature = "fail-on-warnings", deny(warnings))]

pub mod api;
pub mod cff;
pub mod cli_functions;
pub mod cross_section;
pub mod graph;
pub mod tensor;
pub mod h_function_test;
pub mod inspect;
pub mod integrands;
pub mod integrate;
pub mod ltd;
pub mod model;
pub mod observables;
pub mod tests;
pub mod tests_from_pytest;
pub mod utils;

use color_eyre::{Help, Report};
#[allow(unused)]
use colored::Colorize;
use eyre::WrapErr;

use integrands::*;
use num::Complex;
use observables::ObservableSettings;
use observables::PhaseSpaceSelectorSettings;
use std::fs::File;

use serde::{Deserialize, Serialize};

pub const MAX_CORES: usize = 1000;

#[cfg(not(feature = "higher_loops"))]
pub const MAX_LOOP: usize = 3;
#[cfg(feature = "higher_loops")]
pub const MAX_LOOP: usize = 6;

#[derive(Debug, Clone, Default, Deserialize, PartialEq)]
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

#[derive(Debug, Clone, Default, Deserialize)]
pub struct HFunctionSettings {
    pub function: HFunction,
    pub sigma: f64,
    #[serde(default = "_default_true")]
    pub enabled_dampening: bool,
    #[serde(default = "_default_usize_null")]
    pub power: Option<usize>,
}

#[derive(Debug, Clone, Deserialize, PartialEq, Default)]
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

#[derive(Debug, Clone, Deserialize, PartialEq, Default)]
pub enum ParameterizationMapping {
    #[serde(rename = "log")]
    #[default]
    Log,
    #[serde(rename = "linear")]
    Linear,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct GeneralSettings {
    pub debug: usize,
}

#[derive(Debug, Copy, Clone, PartialEq, Deserialize, Default)]
pub enum IntegratedPhase {
    #[serde(rename = "real")]
    #[default]
    Real,
    #[serde(rename = "imag")]
    Imag,
    #[serde(rename = "both")]
    Both,
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct KinematicsSettings {
    pub e_cm: f64,
}

#[derive(Debug, Clone, Default, Deserialize)]
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
}

#[derive(Debug, Clone, Default, Deserialize)]
pub struct ParameterizationSettings {
    pub mode: ParameterizationMode,
    pub mapping: ParameterizationMapping,
    pub b: f64,
    #[serde(default = "_default_input_rescaling")]
    pub input_rescaling: Vec<Vec<(f64, f64)>>,
    #[serde(default = "_default_shifts")]
    pub shifts: Vec<(f64, f64, f64, f64)>,
}

#[derive(Debug, Clone, Default, Deserialize)]
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
