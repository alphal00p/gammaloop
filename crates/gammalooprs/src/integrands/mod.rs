pub mod builtin;
pub mod evaluation;
pub mod process;

use crate::integrands::evaluation::{
    EvaluationMetaData, EvaluationResult, RawBatchEvaluationResult, StabilityResult,
};
// use crate::integrands::process::ProcessIntegrandImpl;
use crate::integrands::builtin::h_function::{HFunctionTestIntegrand, HFunctionTestSettings};
use crate::integrands::process::ProcessIntegrand;
use crate::integrands::process::{amplitude, cross_section_integrand};
use crate::model::Model;
use crate::momentum::FourMomentum;
use crate::observables::{
    ObservableAccumulatorBundle, ObservableFileFormat, ObservableSnapshotBundle,
};
use crate::utils::{F, FloatLike};
use crate::{
    settings::{
        RuntimeSettings,
        runtime::{IntegratorSettings, Precision},
    },
    utils,
};

use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use enum_dispatch::enum_dispatch;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use spenso::algebra::complex::Complex;
use std::fmt::{Display, Formatter};
use std::time::Duration;
use symbolica::domains::float::Real;
use symbolica::numerical_integration::{ContinuousGrid, Grid, Sample};
#[allow(unused_imports)]
use tracing::{debug, error, info, instrument, trace, warn};

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
// #[trait_decode(trait= GammaLoopContext)]
#[allow(non_snake_case)]
#[serde(tag = "type")]
pub enum IntegrandSettings {
    #[serde(rename = "unit_surface")]
    UnitSurface(UnitSurfaceSettings),
    #[serde(rename = "unit_volume")]
    UnitVolume(UnitVolumeSettings),
    #[serde(rename = "h_function_test")]
    HFunctionTest(HFunctionTestSettings),
}

impl Display for IntegrandSettings {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            IntegrandSettings::UnitSurface(_) => write!(f, "unit_surface"),
            IntegrandSettings::UnitVolume(_) => write!(f, "unit_volume"),
            IntegrandSettings::HFunctionTest(_) => {
                write!(f, "h_function_test")
            }
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

    fn name(&self) -> String;

    fn evaluate_sample(
        &mut self,
        sample: &Sample<F<f64>>,
        model: &Model,
        wgt: F<f64>,
        iter: usize,
        use_arb_prec: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<EvaluationResult>;

    fn get_n_dim(&self) -> usize;

    fn get_integrator_settings(&self) -> IntegratorSettings {
        IntegratorSettings::default()
    }

    // In case your integrand supports observable, then overload this function to combine the observables
    fn merge_results<I: HasIntegrand>(&mut self, _other: &mut I, _iter: usize) {}

    // In case your integrand supports observable, then overload this function to write the observables to file
    fn update_results(&mut self, _iter: usize) {}
}

#[derive(Clone)]
pub enum Integrand {
    UnitSurface(UnitSurfaceIntegrand),
    UnitVolume(UnitVolumeIntegrand),
    HFunctionTest(HFunctionTestIntegrand),
    // ProcessIntegrandImpl(ProcessIntegrandImpl),
    ProcessIntegrand(ProcessIntegrand),
}

impl Integrand {
    pub fn evaluate_samples_raw(
        &mut self,
        samples: &[Sample<F<f64>>],
        model: &Model,
        iter: usize,
        use_arb_prec: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<RawBatchEvaluationResult> {
        match self {
            Integrand::ProcessIntegrand(integrand) => {
                integrand.evaluate_samples_raw(model, samples, iter, use_arb_prec, max_eval)
            }
            _ => {
                let mut results = Vec::with_capacity(samples.len());
                for sample in samples {
                    results.push(self.evaluate_sample(
                        sample,
                        model,
                        sample.get_weight(),
                        iter,
                        use_arb_prec,
                        max_eval,
                    )?);
                }
                Ok(RawBatchEvaluationResult {
                    statistics: evaluation::StatisticsCounter::from_evaluation_results(&results),
                    samples: results,
                })
            }
        }
    }

    pub fn process_evaluation_result(&mut self, result: &EvaluationResult) {
        if let Integrand::ProcessIntegrand(integrand) = self {
            integrand.process_evaluation_result(result);
        }
    }

    pub fn merge_runtime_results(&mut self, other: &mut Integrand) -> Result<()> {
        match (self, other) {
            (Integrand::ProcessIntegrand(lhs), Integrand::ProcessIntegrand(rhs)) => {
                lhs.merge_event_processing_runtime(rhs)
            }
            _ => Ok(()),
        }
    }

    pub fn update_runtime_results(&mut self, iter: usize) {
        if let Integrand::ProcessIntegrand(integrand) = self {
            integrand.update_event_processing_runtime(iter);
        }
    }

    pub fn observable_accumulator_bundle(&self) -> Option<ObservableAccumulatorBundle> {
        match self {
            Integrand::ProcessIntegrand(integrand) => integrand.observable_accumulator_bundle(),
            _ => None,
        }
    }

    pub fn observable_snapshot_bundle(&self) -> Option<ObservableSnapshotBundle> {
        match self {
            Integrand::ProcessIntegrand(integrand) => integrand.observable_snapshot_bundle(),
            _ => None,
        }
    }

    pub fn build_observable_snapshots_for_result(
        &self,
        result: &EvaluationResult,
    ) -> Option<ObservableSnapshotBundle> {
        match self {
            Integrand::ProcessIntegrand(integrand) => {
                integrand.build_observable_snapshots_for_result(result)
            }
            _ => None,
        }
    }

    pub fn write_observable_snapshots(
        &self,
        path: impl AsRef<std::path::Path>,
        format: ObservableFileFormat,
    ) -> Result<()> {
        match self {
            Integrand::ProcessIntegrand(integrand) => {
                integrand.write_observable_snapshots(path, format)
            }
            _ => Ok(()),
        }
    }
}

impl HasIntegrand for Integrand {
    fn name(&self) -> String {
        match self {
            Integrand::UnitSurface(_) => "UnitSurface".to_string(),
            Integrand::UnitVolume(_) => "UnitVolume".to_string(),
            Integrand::HFunctionTest(_) => "HFunctionTest".to_string(),
            // Integrand::ProcessIntegrandImpl(_) => "ProcessIntegrandImpl".to_string(),
            Integrand::ProcessIntegrand(i) => i.name(),
        }
    }

    fn create_grid(&self) -> Grid<F<f64>> {
        match self {
            Integrand::UnitSurface(integrand) => integrand.create_grid(),
            Integrand::UnitVolume(integrand) => integrand.create_grid(),
            Integrand::HFunctionTest(integrand) => integrand.create_grid(),
            // Integrand::ProcessIntegrandImpl(integrand) => integrand.create_grid(),
            Integrand::ProcessIntegrand(integrand) => integrand.create_grid(),
        }
    }

    fn evaluate_sample(
        &mut self,
        sample: &Sample<F<f64>>,
        model: &Model,
        wgt: F<f64>,
        iter: usize,
        use_arb_prec: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<EvaluationResult> {
        match self {
            Integrand::UnitSurface(integrand) => {
                integrand.evaluate_sample(sample, model, wgt, iter, use_arb_prec, max_eval)
            }
            Integrand::UnitVolume(integrand) => {
                integrand.evaluate_sample(sample, model, wgt, iter, use_arb_prec, max_eval)
            }
            Integrand::HFunctionTest(integrand) => {
                integrand.evaluate_sample(sample, model, wgt, iter, use_arb_prec, max_eval)
            }
            // Integrand::ProcessIntegrandImpl(integrand) => {
            //     integrand.evaluate_sample(sample,model, wgt, iter, use_f128, max_eval)
            // }
            Integrand::ProcessIntegrand(integrand) => {
                integrand.evaluate_sample(sample, model, wgt, iter, use_arb_prec, max_eval)
            }
        }
    }

    fn get_n_dim(&self) -> usize {
        match self {
            Integrand::UnitSurface(integrand) => integrand.get_n_dim(),
            Integrand::UnitVolume(integrand) => integrand.get_n_dim(),
            Integrand::HFunctionTest(integrand) => integrand.get_n_dim(),
            // Integrand::ProcessIntegrandImpl(integrand) => integrand.get_n_dim(),
            Integrand::ProcessIntegrand(integrand) => integrand.get_n_dim(),
        }
    }

    fn get_integrator_settings(&self) -> IntegratorSettings {
        match self {
            Integrand::UnitSurface(integrand) => integrand.get_integrator_settings(),
            Integrand::UnitVolume(integrand) => integrand.get_integrator_settings(),
            Integrand::HFunctionTest(integrand) => integrand.get_integrator_settings(),
            // Integrand::ProcessIntegrandImpl(integrand) => integrand.get_integrator_settings(),
            Integrand::ProcessIntegrand(integrand) => integrand.get_integrator_settings(),
        }
    }
}

pub(crate) fn integrand_factory(settings: &RuntimeSettings) -> Integrand {
    match settings.hard_coded_integrand.as_ref().unwrap().clone() {
        IntegrandSettings::UnitSurface(integrand_settings) => Integrand::UnitSurface(
            UnitSurfaceIntegrand::new(settings.clone(), integrand_settings),
        ),
        IntegrandSettings::UnitVolume(integrand_settings) => Integrand::UnitVolume(
            UnitVolumeIntegrand::new(settings.clone(), integrand_settings),
        ),
        IntegrandSettings::HFunctionTest(integrand_settings) => Integrand::HFunctionTest(
            HFunctionTestIntegrand::new(settings.clone(), integrand_settings),
        ),
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Default, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
// #[trait_decode(trait= GammaLoopContext)]
pub struct UnitSurfaceSettings {
    pub n_3d_momenta: usize,
}

#[derive(Clone)]
pub struct UnitSurfaceIntegrand {
    pub settings: RuntimeSettings,
    pub n_dim: usize,
    pub n_3d_momenta: usize,
    pub surface: F<f64>,
}

#[allow(unused)]
impl UnitSurfaceIntegrand {
    pub(crate) fn new(
        settings: RuntimeSettings,
        integrand_settings: UnitSurfaceSettings,
    ) -> UnitSurfaceIntegrand {
        let n_dim = integrand_settings.n_3d_momenta * 3 - 1;
        let surface = utils::compute_surface_and_volume(
            integrand_settings.n_3d_momenta * 3 - 1,
            F(settings.kinematics.e_cm),
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
            F::<T>::from_f64(self.settings.kinematics.e_cm * self.settings.kinematics.e_cm),
            &self
                .settings
                .sampling
                .get_parameterization_settings()
                .unwrap(),
        )
    }
}

#[allow(unused)]
impl HasIntegrand for UnitSurfaceIntegrand {
    fn name(&self) -> String {
        "UnitSurfaceIntegrand".to_string()
    }

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
        model: &Model,
        wgt: F<f64>,
        iter: usize,
        use_arb_prec: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<EvaluationResult> {
        let start_evaluate_sample = std::time::Instant::now();

        let xs = match sample {
            Sample::Continuous(_w, v) => v,
            _ => panic!("Wrong sample type"),
        };
        let mut sample_xs = vec![F(self.settings.kinematics.e_cm)];
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

        info!("Sampled loop momenta:");
        for (i, l) in loop_momenta.iter().enumerate() {
            info!("k{} = ( {:-23})", i, format!("{:+.16e}", l),);
        }
        info!("Integrator weight : {:+.16e}", wgt);
        info!("Integrand weight  : {:+.16e}", itg_wgt);
        info!("Sampling jacobian : {:+.16e}", jac);
        info!("Final contribution: {:+.16e}", itg_wgt * jac);

        let is_nan = itg_wgt.is_nan();

        let evaluation_time = before_evaluation.elapsed();

        let evaluation_metadata = EvaluationMetaData {
            total_timing: start_evaluate_sample.elapsed(),
            integrand_evaluation_time: evaluation_time,
            evaluator_evaluation_time: Duration::ZERO,
            average_evaluator_batch_size: None,
            parameterization_time,
            event_processing_time: Duration::ZERO,
            generated_event_count: 0,
            accepted_event_count: 0,
            relative_instability_error: Complex::new_zero(),
            is_nan,
            loop_momenta_escalation: None,
            stability_results: vec![StabilityResult {
                precision: Precision::Double,
                estimated_relative_accuracy: None,
                accepted_as_stable: !is_nan,
                total_time: start_evaluate_sample.elapsed(),
            }],
            evaluator_batch_size_sum: 0,
            evaluator_batch_size_count: 0,
        };

        Ok(EvaluationResult {
            integrand_result: Complex::new(itg_wgt, F(0.)),
            parameterization_jacobian: Some(jac),
            integrator_weight: wgt,
            event_groups: Default::default(),
            evaluation_metadata,
        })
    }
}

#[cfg_attr(feature = "python_api", pyo3::pyclass(from_py_object))]
#[derive(Debug, Clone, Default, Serialize, Deserialize, Encode, Decode, PartialEq, JsonSchema)]
// #[trait_decode(trait= GammaLoopContext)]
pub struct UnitVolumeSettings {
    pub n_3d_momenta: usize,
}

#[derive(Clone)]
pub struct UnitVolumeIntegrand {
    pub settings: RuntimeSettings,
    pub n_dim: usize,
    pub n_3d_momenta: usize,
    pub volume: F<f64>,
}

#[allow(unused)]
impl UnitVolumeIntegrand {
    pub(crate) fn new(
        settings: RuntimeSettings,
        integrand_settings: UnitVolumeSettings,
    ) -> UnitVolumeIntegrand {
        let n_dim = utils::get_n_dim_for_n_loop_momenta(
            &settings.sampling,
            integrand_settings.n_3d_momenta,
            None,
        );
        let volume = utils::compute_surface_and_volume(
            integrand_settings.n_3d_momenta * 3,
            F(settings.kinematics.e_cm),
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
            > F::<T>::from_f64(self.settings.kinematics.e_cm)
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
            F::<T>::from_f64(self.settings.kinematics.e_cm * self.settings.kinematics.e_cm),
            &self
                .settings
                .sampling
                .get_parameterization_settings()
                .unwrap(),
        )
    }
}

#[allow(unused)]
impl HasIntegrand for UnitVolumeIntegrand {
    fn name(&self) -> String {
        "UnitVolumeIntegrand".to_string()
    }
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
        model: &Model,
        wgt: F<f64>,
        iter: usize,
        use_arb_prec: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<EvaluationResult> {
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
        info!("Sampled loop momenta:");
        for (i, l) in loop_momenta.iter().enumerate() {
            info!("k{} = ( {:-23})", i, format!("{:+.16e}", l),);
        }
        info!("Integrator weight : {:+.16e}", wgt);
        info!("Integrand weight  : {:+.16e}", itg_wgt);
        info!("Sampling jacobian : {:+.16e}", jac);
        info!("Final contribution: {:+.16e}", itg_wgt * jac);

        let is_nan = itg_wgt.is_nan();

        let evaluation_time = before_evaluation.elapsed();

        let evaluation_metadata = EvaluationMetaData {
            total_timing: start_evaluate_sample.elapsed(),
            integrand_evaluation_time: evaluation_time,
            evaluator_evaluation_time: Duration::ZERO,
            average_evaluator_batch_size: None,
            parameterization_time,
            event_processing_time: Duration::ZERO,
            generated_event_count: 0,
            accepted_event_count: 0,
            relative_instability_error: Complex::new_zero(),
            is_nan,
            loop_momenta_escalation: None,
            stability_results: vec![StabilityResult {
                precision: Precision::Double,
                estimated_relative_accuracy: None,
                accepted_as_stable: !is_nan,
                total_time: start_evaluate_sample.elapsed(),
            }],
            evaluator_batch_size_sum: 0,
            evaluator_batch_size_count: 0,
        };

        Ok(EvaluationResult {
            integrand_result: Complex::new(itg_wgt, F(0.)),
            parameterization_jacobian: Some(jac),
            integrator_weight: wgt,
            event_groups: Default::default(),
            evaluation_metadata,
        })
    }
}
