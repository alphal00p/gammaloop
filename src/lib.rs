#![cfg_attr(feature = "fail-on-warnings", deny(warnings))]
// #![deny(clippy::all)]
// #![warn(clippy::pedantic)]
#![warn(clippy::all)]
// #![warn(clippy::restriction)]
// #![warn(clippy::nursery)]
// #![warn(clippy::cargo)]
// #![feature(min_specialization)]
//
#[cfg(feature = "python_api")]
pub mod api;
pub mod cff;
pub mod cli_functions;
pub mod cross_section;
pub mod debug_info;
pub mod evaluation_result;
pub mod feyngen;
pub mod gammaloop_integrand;
pub mod graph;
pub mod h_function_test;
pub mod inspect;
pub mod integrands;
pub mod integrate;
pub mod ltd;
pub mod model;
pub mod momentum;
pub mod numerator;
pub mod observables;
pub mod subtraction;

pub mod tests;
pub mod tests_from_pytest;
pub mod utils;

use crate::utils::f128;
use color_eyre::{Help, Report, Result};
#[allow(unused)]
use colored::Colorize;
use cross_section::Amplitude;
use eyre::WrapErr;
use integrands::*;
use log::debug;
use model::Particle;
use momentum::Dep;
use momentum::ExternalMomenta;
use momentum::FourMomentum;
use momentum::Helicity;
use momentum::Polarization;
use momentum::Rotatable;
use momentum::RotationMethod;
use momentum::SignOrZero;
use momentum::Signature;
use momentum::ThreeMomentum;
use numerator::NumeratorSettings;

use observables::ObservableSettings;

use observables::PhaseSpaceSelectorSettings;

use spenso::complex::Complex;
use std::fmt::Display;
use std::fs::File;
use std::sync::atomic::AtomicBool;
use std::sync::Arc;
use symbolica::evaluate::CompileOptions;
use symbolica::evaluate::InlineASM;
use utils::FloatLike;
use utils::F;

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

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
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

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct GeneralSettings {
    pub debug: usize,
    pub use_ltd: bool,
    pub load_compiled_cff: bool,
    pub load_compiled_numerator: bool,
    pub joint_numerator_eval: bool,
    pub amplitude_prefactor: Option<Complex<F<f64>>>,
    pub load_compiled_separate_orientations: bool,
    pub force_orientations: Option<Vec<usize>>,
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
        }
    }
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

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct KinematicsSettings {
    pub e_cm: F<f64>,
    pub externals: Externals,
}

impl Default for KinematicsSettings {
    fn default() -> Self {
        Self {
            e_cm: F(64.),
            externals: Externals::default(),
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct IntegratorSettings {
    pub n_bins: usize,
    pub bin_number_evolution: Option<Vec<usize>>,
    pub min_samples_for_update: usize,
    pub n_start: usize,
    pub n_increase: usize,
    pub n_max: usize,
    pub integrated_phase: IntegratedPhase,
    pub discrete_dim_learning_rate: F<f64>,
    pub continuous_dim_learning_rate: F<f64>,
    pub train_on_avg: bool,
    pub show_max_wgt_info: bool,
    pub max_prob_ratio: F<f64>,
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
    // Runtime settings
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
    #[serde(rename = "subtraction")]
    pub subtraction: SubtractionSettings,
}

impl Settings {
    pub fn sync_with_amplitude(&mut self, amplitude: &Amplitude) -> Result<()> {
        let external_signature = amplitude.external_signature();
        let external_particle_spin = amplitude.external_particle_spin_and_masslessness();

        self.kinematics
            .externals
            .set_dependent_at_end(&external_signature)?;

        self.kinematics
            .externals
            .validate_helicities(&external_particle_spin)?;

        Ok(())
    }

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
    pub result: Vec<F<f64>>,
    pub error: Vec<F<f64>>,
    pub prob: Vec<F<f64>>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct StabilitySettings {
    rotation_axis: Vec<RotationSetting>,
    rotate_numerator: bool,
    levels: Vec<StabilityLevelSetting>,
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

#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct StabilityLevelSetting {
    precision: Precision,
    required_precision_for_re: F<f64>,
    required_precision_for_im: F<f64>,
    escalate_for_large_weight_threshold: F<f64>,
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

#[derive(Serialize, Deserialize, Debug, Clone, Default, Copy, PartialEq)]
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
    pub fn rotation_method(&self) -> RotationMethod {
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
    pub fn rotation_function<'a, T: FloatLike + 'a>(
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

    fn as_str(&self) -> String {
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

#[derive(Serialize, Deserialize, Debug, Clone, Default, PartialEq, Copy, Hash, Eq)]
pub enum Precision {
    Single, // this might be useful for eventual deployment on gpu
    #[default]
    Double,
    Quad,
    Arb(usize),
}

impl Precision {
    fn as_string(self) -> String {
        match self {
            Self::Single => "f32".to_owned(),
            Self::Double => "f64".to_owned(),
            Self::Quad => "f128".to_owned(),
            Self::Arb(prec) => format!("arb_prec_{}", prec),
        }
    }
}

impl Display for Precision {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.as_string())
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
#[serde(tag = "type", content = "data")]
pub enum Externals {
    #[serde(rename = "constant")]
    Constant {
        momenta: Vec<ExternalMomenta<F<f64>>>,
        helicities: Vec<Helicity>,
    },
    // add different type of pdfs here when needed
}

impl Rotatable for Externals {
    fn rotate(&self, rotation: &momentum::Rotation) -> Self {
        match self {
            Externals::Constant {
                momenta,
                helicities,
            } => {
                let momenta = momenta.iter().map(|m| m.rotate(rotation)).collect();
                Externals::Constant {
                    momenta,
                    helicities: helicities.clone(),
                }
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum Polarizations {
    Constant {
        polarizations: Vec<Polarization<Complex<F<f64>>>>,
    },
    None,
}

impl Rotatable for Polarizations {
    fn rotate(&self, rotation: &momentum::Rotation) -> Self {
        match self {
            Polarizations::Constant { polarizations } => {
                let polarizations = polarizations.iter().map(|p| p.rotate(rotation)).collect();
                Polarizations::Constant { polarizations }
            }
            Polarizations::None => Polarizations::None,
        }
    }
}

use thiserror::Error;

#[derive(Error, Debug)]
pub enum ExternalsValidationError {
    #[error("There should be exactly one dependent external momentum")]
    WrongNumberOfDependentMomenta,
    #[error("Found {0} momenta, expected {1}")]
    WrongNumberOfMomentaExpected(usize, usize),
    #[error("Found {0} helicities, expected {1}")]
    WrongNumberOfHelicities(usize, usize),
    #[error("Massless vector cannot have zero helicity: pos {0}")]
    MasslessVectorZeroHelicity(usize),
    #[error("Spinors cannot have zero helicity at pos {0}")]
    SpinorZeroHelicity(usize),
    #[error("Scalars cannot have non-zero helicity at pos {0}")]
    ScalarNonZeroHelicity(usize),
    #[error("{0} is an Unsuported external spin for pos {0}")]
    UnsupportedSpin(isize, usize),
}

impl Externals {
    pub fn validate_helicities(
        &self,
        spins: &[(isize, bool)],
    ) -> Result<(), ExternalsValidationError> {
        match self {
            Externals::Constant { helicities, .. } => {
                if helicities.len() == spins.len() {
                    for (i, (h, (s, is_massless))) in
                        helicities.iter().zip(spins.iter()).enumerate()
                    {
                        match *s {
                            1 => {
                                if !h.is_zero() {
                                    return Err(ExternalsValidationError::ScalarNonZeroHelicity(i));
                                }
                            }
                            2 => {
                                if h.is_zero() {
                                    return Err(ExternalsValidationError::SpinorZeroHelicity(i));
                                }
                            }
                            3 => {
                                if h.is_zero() && *is_massless {
                                    return Err(
                                        ExternalsValidationError::MasslessVectorZeroHelicity(i),
                                    );
                                }
                            }
                            s => return Err(ExternalsValidationError::UnsupportedSpin(s, i)),
                        }
                    }
                    Ok(())
                } else {
                    Err(ExternalsValidationError::WrongNumberOfHelicities(
                        helicities.len(),
                        spins.len(),
                    ))
                }
            }
        }
    }

    pub fn generate_polarizations(
        &self,
        external_particles: &[Arc<Particle>],
        external_signature: &Signature,
    ) -> Polarizations {
        let mut polarizations = vec![];

        let dep_ext = self.get_dependent_externals(external_signature);
        let helicities = self.get_helicities();
        for (((ext_mom, hel), p), s) in dep_ext
            .iter()
            .zip(helicities.iter())
            .zip(external_particles.iter())
            .zip(external_signature.iter())
        {
            match s {
                SignOrZero::Minus => {
                    polarizations.push(p.incoming_polarization(ext_mom, *hel));
                }
                SignOrZero::Plus => {
                    polarizations.push(p.outgoing_polarization(ext_mom, *hel));
                }
                _ => {
                    panic!("Edge type should not be virtual")
                }
            }
        }

        Polarizations::Constant { polarizations }
    }

    pub fn set_dependent_at_end(
        &mut self,
        signature: &Signature,
    ) -> Result<(), ExternalsValidationError> {
        match self {
            Externals::Constant { momenta, .. } => {
                let mut sum: FourMomentum<F<f128>> = FourMomentum::from([F(0.0); 4]).higher();
                let mut pos_dep = 0;
                let mut n_dep = 0;

                let mut dependent_sign = SignOrZero::Plus;

                for ((i, m), s) in momenta.iter().enumerate().zip(signature.iter()) {
                    if let Ok(a) = FourMomentum::try_from(*m) {
                        sum -= *s * a.higher();
                    } else {
                        pos_dep = i;
                        n_dep += 1;
                        dependent_sign = *s;
                    }
                }
                if n_dep == 1 {
                    momenta[pos_dep] = (dependent_sign * sum.lower()).into();
                    let len = momenta.len();
                    momenta[len - 1] = ExternalMomenta::Dependent(Dep::Dep);
                } else if n_dep == 0 {
                    debug!("No dependent momentum found, adding the sum at the end");
                    momenta.push(ExternalMomenta::Dependent(Dep::Dep));
                } else {
                    return Err(ExternalsValidationError::WrongNumberOfDependentMomenta);
                }

                let len = momenta.len();
                if len == signature.len() {
                    Ok(())
                } else {
                    Err(ExternalsValidationError::WrongNumberOfMomentaExpected(
                        len - 1,
                        signature.len() - 1,
                    ))
                }
            }
        }
    }

    pub fn get_dependent_externals<T: FloatLike>(
        &self,
        external_signature: &Signature,
    ) -> Vec<FourMomentum<F<T>>>
// where
    //     T::Higher: PrecisionUpgradable<Lower = T> + FloatLike,
    {
        match self {
            Externals::Constant { momenta, .. } => {
                let mut sum: FourMomentum<F<T>> = FourMomentum::from([
                    F::<T>::from_f64(0.0),
                    F::from_f64(0.0),
                    F::from_f64(0.0),
                    F::from_f64(0.0),
                ]);
                // .higher();
                let mut pos_dep = 0;

                let mut dependent_sign = SignOrZero::Plus;

                let mut dependent_momenta = vec![];

                for ((i, m), s) in momenta.iter().enumerate().zip(external_signature.iter()) {
                    if let Ok(a) = FourMomentum::try_from(*m) {
                        // println!("external{i}: {}", a);
                        let a = FourMomentum::<F<T>>::from_ff64(&a);
                        sum -= *s * a.clone(); //.higher();
                        dependent_momenta.push(a);
                    } else {
                        pos_dep = i;
                        dependent_sign = *s;
                        dependent_momenta.push(sum.clone()); //.lower());
                    }
                }

                dependent_momenta[pos_dep] = dependent_sign * sum; //.lower();
                dependent_momenta
            }
        }
    }

    #[allow(unused_variables)]
    #[inline]
    pub fn get_indep_externals(&self) -> Vec<FourMomentum<F<f64>>> {
        match self {
            Externals::Constant {
                momenta,
                helicities,
            } => {
                let momenta: Vec<FourMomentum<_>> = momenta
                    .iter()
                    .flat_map(|e| FourMomentum::try_from(*e))
                    .collect();
                momenta
            }
        }
    }

    pub fn get_helicities(&self) -> &[Helicity] {
        match self {
            Externals::Constant { helicities, .. } => helicities,
        }
    }

    pub fn pdf(&self, _x_space_point: &[F<f64>]) -> F<f64> {
        match self {
            Externals::Constant { .. } => F(1.0),
        }
    }
}

#[test]
fn external_inv() {
    let mut ext = Externals::Constant {
        momenta: vec![[F(1.), F(2.), F(3.), F(4.)].into(); 3],
        helicities: vec![Helicity::Plus; 4],
    };

    let signs: Signature = [1i8, 1, 1, 1].into_iter().collect();
    ext.set_dependent_at_end(&signs).unwrap();

    let momenta = vec![
        ExternalMomenta::Dependent(Dep::Dep),
        [F(1.), F(2.), F(3.), F(4.)].into(),
        [F(1.), F(2.), F(3.), F(4.)].into(),
        [F(-3.), F(-6.), F(-9.), F(-12.)].into(),
    ];
    let mut ext2 = Externals::Constant {
        momenta,
        helicities: vec![Helicity::Plus; 4],
    };

    ext2.set_dependent_at_end(&signs).unwrap();

    assert_eq!(ext, ext2);
}

impl Default for Externals {
    fn default() -> Self {
        Externals::Constant {
            momenta: vec![],
            helicities: vec![],
        }
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct MultiChannelingSettings {
    pub alpha: f64,
}

impl Default for MultiChannelingSettings {
    fn default() -> Self {
        Self { alpha: 3.0 }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GammaloopTropicalSamplingSettings {
    pub upcast_on_failure: bool,
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
    pub fn into_tropical_sampling_settings(
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
    TropicalSampling(GammaloopTropicalSamplingSettings),
}

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct LocalCounterTermSettings {
    pub dampen_integrable_singularity: IntegrableSingularityDampener,
    pub uv_localisation: UVLocalisationSettings,
}

#[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq)]
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct UVLocalisationSettings {
    pub sliver_width: f64,
    pub dynamic_width: bool,
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

#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq)]
pub struct IntegratedCounterTermSettings {
    range: IntegratedCounterTermRange,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct OverlapSettings {
    pub force_global_center: Option<Vec<[f64; 3]>>,
    pub check_global_center: bool,
    pub try_origin: bool,
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

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct SubtractionSettings {
    pub local_ct_settings: LocalCounterTermSettings,
    pub integrated_ct_settings: IntegratedCounterTermSettings,
    pub overlap_settings: OverlapSettings,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExportSettings {
    // Generation Time settings
    pub compile_cff: bool,
    pub numerator_settings: NumeratorSettings,
    pub cpe_rounds_cff: Option<usize>,
    pub compile_separate_orientations: bool,
    pub gammaloop_compile_options: GammaloopCompileOptions,
    pub tropical_subgraph_table_settings: TropicalSubgraphTableSettings,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GammaloopCompileOptions {
    pub inline_asm: bool,
    pub optimization_level: usize,
    pub fast_math: bool,
    pub unsafe_math: bool,
    pub compiler: String,
    pub custom: Vec<String>,
}

impl GammaloopCompileOptions {
    pub fn inline_asm(&self) -> InlineASM {
        if self.inline_asm {
            InlineASM::default()
        } else {
            InlineASM::None
        }
    }

    #[allow(clippy::needless_update)]
    pub fn to_symbolica_compile_options(&self) -> CompileOptions {
        CompileOptions {
            optimization_level: self.optimization_level,
            fast_math: self.fast_math,
            unsafe_math: self.unsafe_math,
            compiler: self.compiler.clone(),
            custom: self.custom.clone(),
            ..CompileOptions::default()
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TropicalSubgraphTableSettings {
    pub panic_on_fail: bool,
    pub target_omega: f64,
}
