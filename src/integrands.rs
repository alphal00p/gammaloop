use crate::evaluation_result::{EvaluationMetaData, EvaluationResult};
use crate::gammaloop_integrand::GammaLoopIntegrand;
use crate::h_function_test::{HFunctionTestIntegrand, HFunctionTestSettings};
use crate::observables::EventManager;
use crate::utils::FloatLike;
use crate::{utils, IntegratorSettings, Precision, Settings};
use enum_dispatch::enum_dispatch;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use lorentz_vector::LorentzVector;
use num::Complex;
use serde::{Deserialize, Serialize};
use std::fmt::{Display, Formatter};
use symbolica::numerical_integration::{ContinuousGrid, Grid, Sample};

#[derive(Debug, Clone, Serialize, Deserialize)]
#[allow(non_snake_case)]
#[serde(tag = "type")]
pub enum IntegrandSettings {
    #[serde(rename = "unit_surface")]
    UnitSurface(UnitSurfaceSettings),
    #[serde(rename = "unit_volume")]
    UnitVolume(UnitVolumeSettings),
    #[serde(rename = "h_function_test")]
    HFunctionTest(HFunctionTestSettings),
    #[serde(rename = "gamma_loop")]
    GammaLoop,
}

impl Display for IntegrandSettings {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            IntegrandSettings::UnitSurface(_) => write!(f, "unit_surface"),
            IntegrandSettings::UnitVolume(_) => write!(f, "unit_volume"),
            IntegrandSettings::HFunctionTest(_) => {
                write!(f, "h_function_test")
            }
            IntegrandSettings::GammaLoop => write!(f, "gamma_loop"),
        }
    }
}

impl Default for IntegrandSettings {
    fn default() -> IntegrandSettings {
        IntegrandSettings::UnitSurface(UnitSurfaceSettings { n_3d_momenta: 1 })
    }
}

#[enum_dispatch]
pub trait HasIntegrand {
    fn create_grid(&self) -> Grid<f64>;

    fn evaluate_sample(
        &self,
        sample: &Sample<f64>,
        wgt: f64,
        iter: usize,
        use_f128: bool,
        max_eval: f64,
    ) -> EvaluationResult;

    fn get_n_dim(&self) -> usize;

    fn get_integrator_settings(&self) -> IntegratorSettings {
        IntegratorSettings::default()
    }

    // In case your integrand supports observable, then overload this function to combine the observables
    fn merge_results<I: HasIntegrand>(&mut self, _other: &mut I, _iter: usize) {}

    // In case your integrand supports observable, then overload this function to write the observables to file
    fn update_results(&mut self, _iter: usize) {}

    // In case your integrand has an EventManager, overload this function to return it
    fn get_event_manager_mut(&mut self) -> &mut EventManager {
        panic!("This integrand does not have an EventManager");
    }
}

#[enum_dispatch(HasIntegrand)]
#[derive(Clone)]
pub enum Integrand {
    UnitSurface(UnitSurfaceIntegrand),
    UnitVolume(UnitVolumeIntegrand),
    HFunctionTest(HFunctionTestIntegrand),
    GammaLoopIntegrand(GammaLoopIntegrand),
}

pub fn integrand_factory(settings: &Settings) -> Integrand {
    match settings.hard_coded_integrand.clone() {
        IntegrandSettings::UnitSurface(integrand_settings) => Integrand::UnitSurface(
            UnitSurfaceIntegrand::new(settings.clone(), integrand_settings),
        ),
        IntegrandSettings::UnitVolume(integrand_settings) => Integrand::UnitVolume(
            UnitVolumeIntegrand::new(settings.clone(), integrand_settings),
        ),
        IntegrandSettings::HFunctionTest(integrand_settings) => Integrand::HFunctionTest(
            HFunctionTestIntegrand::new(settings.clone(), integrand_settings),
        ),
        IntegrandSettings::GammaLoop => {
            unimplemented!("unsupported integrand construction method, please use the exporter to generate the integrand");
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct UnitSurfaceSettings {
    pub n_3d_momenta: usize,
}

#[derive(Clone)]
pub struct UnitSurfaceIntegrand {
    pub settings: Settings,
    pub n_dim: usize,
    pub n_3d_momenta: usize,
    pub surface: f64,
}

#[allow(unused)]
impl UnitSurfaceIntegrand {
    pub fn new(
        settings: Settings,
        integrand_settings: UnitSurfaceSettings,
    ) -> UnitSurfaceIntegrand {
        let n_dim = utils::get_n_dim_for_n_loop_momenta(
            &settings,
            integrand_settings.n_3d_momenta,
            true,
            None,
        );
        let surface = utils::compute_surface_and_volume(
            integrand_settings.n_3d_momenta * 3 - 1,
            settings.kinematics.e_cm,
        )
        .0;
        UnitSurfaceIntegrand {
            settings,
            n_3d_momenta: integrand_settings.n_3d_momenta,
            n_dim,
            surface,
        }
    }

    fn evaluate_numerator<T: FloatLike>(&self, loop_momenta: &[LorentzVector<T>]) -> T {
        T::from_f64(1.0).unwrap()
    }

    fn parameterize<T: FloatLike>(&self, xs: &[T]) -> (Vec<[T; 3]>, T) {
        utils::global_parameterize(
            xs,
            Into::<T>::into(self.settings.kinematics.e_cm * self.settings.kinematics.e_cm),
            &self.settings,
            true,
        )
    }
}

#[allow(unused)]
impl HasIntegrand for UnitSurfaceIntegrand {
    fn create_grid(&self) -> Grid<f64> {
        Grid::Continuous(ContinuousGrid::new(
            self.n_dim,
            self.settings.integrator.n_bins,
            self.settings.integrator.min_samples_for_update,
            self.settings.integrator.bin_number_evolution.clone(),
            self.settings.integrator.train_on_avg,
        ))
    }

    fn get_n_dim(&self) -> usize {
        self.n_dim
    }

    fn evaluate_sample(
        &self,
        sample: &Sample<f64>,
        wgt: f64,
        iter: usize,
        use_f128: bool,
        max_eval: f64,
    ) -> EvaluationResult {
        let start_evaluate_sample = std::time::Instant::now();

        let xs = match sample {
            Sample::Continuous(_w, v) => v,
            _ => panic!("Wrong sample type"),
        };
        let mut sample_xs = vec![self.settings.kinematics.e_cm];
        sample_xs.extend(xs);

        let before_parameterization = std::time::Instant::now();
        let (moms, jac) = self.parameterize(sample_xs.as_slice());
        let mut loop_momenta = vec![];
        for m in &moms {
            loop_momenta.push(LorentzVector {
                t: ((m[0] + m[1] + m[2]) * (m[0] + m[1] + m[2])).sqrt(),
                x: m[0],
                y: m[1],
                z: m[2],
            });
        }

        let parameterization_time = before_parameterization.elapsed();

        let before_evaluation = std::time::Instant::now();
        let mut itg_wgt = self.evaluate_numerator(loop_momenta.as_slice());
        // Normalize the integral
        itg_wgt /= self.surface;
        if self.settings.general.debug > 1 {
            info!("Sampled loop momenta:");
            for (i, l) in loop_momenta.iter().enumerate() {
                info!(
                    "k{} = ( {:-23}, {:-23}, {:-23}, {:-23} )",
                    i,
                    format!("{:+.16e}", l.t),
                    format!("{:+.16e}", l.x),
                    format!("{:+.16e}", l.y),
                    format!("{:+.16e}", l.z)
                );
            }
            info!("Integrator weight : {:+.16e}", wgt);
            info!("Integrand weight  : {:+.16e}", itg_wgt);
            info!("Sampling jacobian : {:+.16e}", jac);
            info!("Final contribution: {:+.16e}", itg_wgt * jac);
        }

        let is_nan = itg_wgt.is_nan();

        let evaluation_time = before_evaluation.elapsed();

        let evaluation_metadata = EvaluationMetaData {
            total_timing: start_evaluate_sample.elapsed(),
            rep3d_evaluation_time: evaluation_time,
            parameterization_time,
            relative_instability_error: Complex::new(0., 0.),
            highest_precision: Precision::Double,
            is_nan,
        };

        EvaluationResult {
            integrand_result: Complex::new(itg_wgt, 0.) * jac,
            integrator_weight: wgt,
            event_buffer: vec![],
            evaluation_metadata,
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct UnitVolumeSettings {
    pub n_3d_momenta: usize,
}

#[derive(Clone)]
pub struct UnitVolumeIntegrand {
    pub settings: Settings,
    pub n_dim: usize,
    pub n_3d_momenta: usize,
    pub volume: f64,
}

#[allow(unused)]
impl UnitVolumeIntegrand {
    pub fn new(settings: Settings, integrand_settings: UnitVolumeSettings) -> UnitVolumeIntegrand {
        let n_dim = utils::get_n_dim_for_n_loop_momenta(
            &settings,
            integrand_settings.n_3d_momenta,
            false,
            None,
        );
        let volume = utils::compute_surface_and_volume(
            integrand_settings.n_3d_momenta * 3,
            settings.kinematics.e_cm,
        )
        .1;
        UnitVolumeIntegrand {
            settings,
            n_3d_momenta: integrand_settings.n_3d_momenta,
            n_dim,
            volume,
        }
    }

    fn evaluate_numerator<T: FloatLike>(&self, loop_momenta: &[LorentzVector<T>]) -> T {
        if loop_momenta
            .iter()
            .map(|l| l.spatial_squared())
            .sum::<T>()
            .sqrt()
            > Into::<T>::into(self.settings.kinematics.e_cm)
        {
            T::from_f64(0.0).unwrap()
        } else {
            T::from_f64(1.0).unwrap()
        }
    }

    fn parameterize<T: FloatLike>(&self, xs: &[T]) -> (Vec<[T; 3]>, T) {
        utils::global_parameterize(
            xs,
            Into::<T>::into(self.settings.kinematics.e_cm * self.settings.kinematics.e_cm),
            &self.settings,
            false,
        )
    }
}

#[allow(unused)]
impl HasIntegrand for UnitVolumeIntegrand {
    fn create_grid(&self) -> Grid<f64> {
        Grid::Continuous(ContinuousGrid::new(
            self.n_dim,
            self.settings.integrator.n_bins,
            self.settings.integrator.min_samples_for_update,
            self.settings.integrator.bin_number_evolution.clone(),
            self.settings.integrator.train_on_avg,
        ))
    }

    fn get_n_dim(&self) -> usize {
        self.n_dim
    }

    fn evaluate_sample(
        &self,
        sample: &Sample<f64>,
        wgt: f64,
        iter: usize,
        use_f128: bool,
        max_eval: f64,
    ) -> EvaluationResult {
        let start_evaluate_sample = std::time::Instant::now();

        let xs = match sample {
            Sample::Continuous(_w, v) => v,
            _ => panic!("Wrong sample type"),
        };

        let before_parameterization = std::time::Instant::now();

        let (moms, jac) = self.parameterize(xs);
        let mut loop_momenta = vec![];
        for m in &moms {
            loop_momenta.push(LorentzVector {
                t: 0.,
                x: m[0],
                y: m[1],
                z: m[2],
            });
        }

        let parameterization_time = before_parameterization.elapsed();

        let before_evaluation = std::time::Instant::now();
        let mut itg_wgt = self.evaluate_numerator(loop_momenta.as_slice());
        // Normalize the integral
        itg_wgt /= self.volume;
        if self.settings.general.debug > 1 {
            info!("Sampled loop momenta:");
            for (i, l) in loop_momenta.iter().enumerate() {
                info!(
                    "k{} = ( {:-23}, {:-23}, {:-23}, {:-23} )",
                    i,
                    format!("{:+.16e}", l.t),
                    format!("{:+.16e}", l.x),
                    format!("{:+.16e}", l.y),
                    format!("{:+.16e}", l.z)
                );
            }
            info!("Integrator weight : {:+.16e}", wgt);
            info!("Integrand weight  : {:+.16e}", itg_wgt);
            info!("Sampling jacobian : {:+.16e}", jac);
            info!("Final contribution: {:+.16e}", itg_wgt * jac);
        }

        let is_nan = itg_wgt.is_nan();

        let evaluation_time = before_evaluation.elapsed();

        let evaluation_metadata = EvaluationMetaData {
            total_timing: start_evaluate_sample.elapsed(),
            rep3d_evaluation_time: evaluation_time,
            parameterization_time,
            relative_instability_error: Complex::new(0., 0.),
            highest_precision: Precision::Double,
            is_nan,
        };

        EvaluationResult {
            integrand_result: Complex::new(itg_wgt, 0.) * jac,
            integrator_weight: wgt,
            event_buffer: vec![],
            evaluation_metadata,
        }
    }
}
