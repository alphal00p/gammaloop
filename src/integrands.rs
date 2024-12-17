use crate::evaluation_result::{EvaluationMetaData, EvaluationResult};
use crate::gammaloop_integrand::GammaLoopIntegrand;
use crate::h_function_test::{HFunctionTestIntegrand, HFunctionTestSettings};
use crate::momentum::FourMomentum;
use crate::observables::EventManager;
use crate::utils::{FloatLike, F};
use crate::{utils, IntegratorSettings, Precision, Settings};
use enum_dispatch::enum_dispatch;
#[allow(unused_imports)]
use log::{debug, error, info, trace, warn};
use serde::{Deserialize, Serialize};
use spenso::complex::Complex;
use std::fmt::{Display, Formatter};
use symbolica::domains::float::{NumericalFloatLike, Real};
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
        IntegrandSettings::UnitSurface(UnitSurfaceSettings { n_3d_momenta: 11 })
    }
}

#[enum_dispatch]
pub trait HasIntegrand {
    fn create_grid(&self) -> Grid<F<f64>>;

    fn evaluate_sample(
        &mut self,
        sample: &Sample<F<f64>>,
        wgt: F<f64>,
        iter: usize,
        use_f128: bool,
        max_eval: F<f64>,
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

#[derive(Clone)]
pub enum Integrand {
    UnitSurface(UnitSurfaceIntegrand),
    UnitVolume(UnitVolumeIntegrand),
    HFunctionTest(HFunctionTestIntegrand),
    GammaLoopIntegrand(GammaLoopIntegrand),
}

impl HasIntegrand for Integrand {
    fn create_grid(&self) -> Grid<F<f64>> {
        match self {
            Integrand::UnitSurface(integrand) => integrand.create_grid(),
            Integrand::UnitVolume(integrand) => integrand.create_grid(),
            Integrand::HFunctionTest(integrand) => integrand.create_grid(),
            Integrand::GammaLoopIntegrand(integrand) => integrand.create_grid(),
        }
    }

    fn evaluate_sample(
        &mut self,
        sample: &Sample<F<f64>>,
        wgt: F<f64>,
        iter: usize,
        use_f128: bool,
        max_eval: F<f64>,
    ) -> EvaluationResult {
        match self {
            Integrand::UnitSurface(integrand) => {
                integrand.evaluate_sample(sample, wgt, iter, use_f128, max_eval)
            }
            Integrand::UnitVolume(integrand) => {
                integrand.evaluate_sample(sample, wgt, iter, use_f128, max_eval)
            }
            Integrand::HFunctionTest(integrand) => {
                integrand.evaluate_sample(sample, wgt, iter, use_f128, max_eval)
            }
            Integrand::GammaLoopIntegrand(integrand) => {
                integrand.evaluate_sample(sample, wgt, iter, use_f128, max_eval)
            }
        }
    }

    fn get_n_dim(&self) -> usize {
        match self {
            Integrand::UnitSurface(integrand) => integrand.get_n_dim(),
            Integrand::UnitVolume(integrand) => integrand.get_n_dim(),
            Integrand::HFunctionTest(integrand) => integrand.get_n_dim(),
            Integrand::GammaLoopIntegrand(integrand) => integrand.get_n_dim(),
        }
    }

    fn get_integrator_settings(&self) -> IntegratorSettings {
        match self {
            Integrand::UnitSurface(integrand) => integrand.get_integrator_settings(),
            Integrand::UnitVolume(integrand) => integrand.get_integrator_settings(),
            Integrand::HFunctionTest(integrand) => integrand.get_integrator_settings(),
            Integrand::GammaLoopIntegrand(integrand) => integrand.get_integrator_settings(),
        }
    }

    fn merge_results<I: HasIntegrand>(&mut self, other: &mut I, iter: usize) {
        match self {
            Integrand::UnitSurface(integrand) => integrand.merge_results(other, iter),
            Integrand::UnitVolume(integrand) => integrand.merge_results(other, iter),
            Integrand::HFunctionTest(integrand) => integrand.merge_results(other, iter),
            Integrand::GammaLoopIntegrand(integrand) => integrand.merge_results(other, iter),
        }
    }
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
    pub surface: F<f64>,
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

    fn evaluate_numerator<T: FloatLike>(&self, loop_momenta: &[FourMomentum<F<T>>]) -> F<T> {
        loop_momenta[0].temporal.value.one()
    }

    fn parameterize<T: FloatLike>(&self, xs: &[F<T>]) -> (Vec<[F<T>; 3]>, F<T>) {
        let zero = xs[0].zero();
        utils::global_parameterize(
            xs,
            F::<T>::from_ff64(self.settings.kinematics.e_cm * self.settings.kinematics.e_cm),
            &self.settings,
            true,
        )
    }
}

#[allow(unused)]
impl HasIntegrand for UnitSurfaceIntegrand {
    fn create_grid(&self) -> Grid<F<f64>> {
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
        &mut self,
        sample: &Sample<F<f64>>,
        wgt: F<f64>,
        iter: usize,
        use_f128: bool,
        max_eval: F<f64>,
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
            loop_momenta.push(FourMomentum::from_args(
                ((m[0] + m[1] + m[2]) * (m[0] + m[1] + m[2])).sqrt(),
                m[0],
                m[1],
                m[2],
            ));
        }

        let parameterization_time = before_parameterization.elapsed();

        let before_evaluation = std::time::Instant::now();
        let mut itg_wgt = self.evaluate_numerator(loop_momenta.as_slice());
        // Normalize the integral
        itg_wgt /= self.surface;
        if self.settings.general.debug > 1 {
            info!("Sampled loop momenta:");
            for (i, l) in loop_momenta.iter().enumerate() {
                info!("k{} = ( {:-23})", i, format!("{:+.16e}", l),);
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
            relative_instability_error: Complex::new_zero(),
            highest_precision: Precision::Double,
            is_nan,
        };

        EvaluationResult {
            integrand_result: Complex::new(itg_wgt, F(0.)) * jac,
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
    pub volume: F<f64>,
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

    fn evaluate_numerator<T: FloatLike>(&self, loop_momenta: &[FourMomentum<F<T>>]) -> F<T> {
        let zero = loop_momenta[0].temporal.value.zero();
        if loop_momenta
            .iter()
            .map(|l| l.spatial.norm_squared())
            .reduce(|acc, e| acc + &e)
            .unwrap_or(zero.clone())
            .sqrt()
            > F::<T>::from_ff64(self.settings.kinematics.e_cm)
        {
            zero
        } else {
            zero.one()
        }
    }

    fn parameterize<T: FloatLike>(&self, xs: &[F<T>]) -> (Vec<[F<T>; 3]>, F<T>) {
        let zero = xs[0].zero();
        utils::global_parameterize(
            xs,
            F::<T>::from_ff64(self.settings.kinematics.e_cm * self.settings.kinematics.e_cm),
            &self.settings,
            false,
        )
    }
}

#[allow(unused)]
impl HasIntegrand for UnitVolumeIntegrand {
    fn create_grid(&self) -> Grid<F<f64>> {
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
        &mut self,
        sample: &Sample<F<f64>>,
        wgt: F<f64>,
        iter: usize,
        use_f128: bool,
        max_eval: F<f64>,
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
            loop_momenta.push(FourMomentum::new(F(0.).into(), (*m).into()));
        }

        let parameterization_time = before_parameterization.elapsed();

        let before_evaluation = std::time::Instant::now();
        let mut itg_wgt = self.evaluate_numerator(loop_momenta.as_slice());
        // Normalize the integral
        itg_wgt /= self.volume;
        if self.settings.general.debug > 1 {
            info!("Sampled loop momenta:");
            for (i, l) in loop_momenta.iter().enumerate() {
                info!("k{} = ( {:-23})", i, format!("{:+.16e}", l),);
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
            relative_instability_error: Complex::new_zero(),
            highest_precision: Precision::Double,
            is_nan,
        };

        EvaluationResult {
            integrand_result: Complex::new(itg_wgt, F(0.)) * jac,
            integrator_weight: wgt,
            event_buffer: vec![],
            evaluation_metadata,
        }
    }
}
