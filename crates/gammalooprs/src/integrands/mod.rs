pub mod evaluation;
pub mod process;

use crate::integrands::evaluation::{EvaluationResult, RawBatchEvaluationResult};
use crate::integrands::process::{EvaluationTarget, ProcessIntegrand};
use crate::integrands::process::{amplitude, cross_section};
use crate::model::Model;
use crate::observables::{
    ObservableAccumulatorBundle, ObservableFileFormat, ObservableSnapshotBundle,
};
use crate::settings::runtime::IntegratorSettings;
use crate::utils::F;
#[cfg(test)]
use crate::{is_interrupted, settings::RuntimeSettings};

use color_eyre::Result;
use enum_dispatch::enum_dispatch;
use spenso::algebra::complex::Complex;
#[cfg(test)]
use symbolica::numerical_integration::ContinuousGrid;
use symbolica::numerical_integration::{Grid, Sample};

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

    fn merge_results<I: HasIntegrand>(&mut self, _other: &mut I, _iter: usize) {}

    fn update_results(&mut self, _iter: usize) {}
}

/// The only production integrand container. Sampling acceptance tests use
/// [`ProcessIntegrand`] with a reference overlay, so no standalone unit-volume
/// or profile integrand is kept as a second owner of parameterization logic.
#[derive(Clone)]
#[allow(clippy::large_enum_variant)]
pub enum Integrand {
    ProcessIntegrand(Box<ProcessIntegrand>),
    #[cfg(test)]
    TestProbe(TestProbeIntegrand),
}

impl Integrand {
    pub fn evaluate_samples_raw(
        &mut self,
        samples: &[Sample<F<f64>>],
        target: EvaluationTarget<'_>,
        iter: usize,
        use_arb_prec: bool,
        stop_on_interrupt: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<RawBatchEvaluationResult> {
        match self {
            Integrand::ProcessIntegrand(integrand) => integrand.evaluate_samples_raw(
                target,
                samples,
                iter,
                use_arb_prec,
                stop_on_interrupt,
                max_eval,
            ),
            #[cfg(test)]
            Integrand::TestProbe(_) => {
                let EvaluationTarget::Physical(model) = target else {
                    return Err(color_eyre::eyre::eyre!(
                        "reference overlays require a generated process"
                    ));
                };
                let mut results = Vec::with_capacity(samples.len());
                for sample in samples {
                    if stop_on_interrupt && is_interrupted() {
                        break;
                    }
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

    #[allow(clippy::too_many_arguments)]
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_samples_raw_with_estimate(
        &mut self,
        samples: &[Sample<F<f64>>],
        target: EvaluationTarget<'_>,
        iter: usize,
        use_arb_prec: bool,
        stop_on_interrupt: bool,
        max_eval: Complex<F<f64>>,
        integral_estimate: Option<(f64, f64)>,
    ) -> Result<RawBatchEvaluationResult> {
        match self {
            Integrand::ProcessIntegrand(integrand) => integrand.evaluate_samples_raw_with_estimate(
                target,
                samples,
                iter,
                use_arb_prec,
                stop_on_interrupt,
                max_eval,
                integral_estimate,
            ),
            #[cfg(test)]
            Integrand::TestProbe(_) => self.evaluate_samples_raw(
                samples,
                target,
                iter,
                use_arb_prec,
                stop_on_interrupt,
                max_eval,
            ),
        }
    }

    pub fn process_evaluation_result(&mut self, result: &EvaluationResult) {
        match self {
            Integrand::ProcessIntegrand(integrand) => integrand.process_evaluation_result(result),
            #[cfg(test)]
            Integrand::TestProbe(_) => {}
        }
    }

    pub fn merge_runtime_results(&mut self, other: &mut Integrand) -> Result<()> {
        match (self, other) {
            (Integrand::ProcessIntegrand(lhs), Integrand::ProcessIntegrand(rhs)) => {
                lhs.merge_event_processing_runtime(rhs)
            }
            #[cfg(test)]
            _ => Ok(()),
        }
    }

    pub fn update_runtime_results(&mut self, iter: usize) {
        match self {
            Integrand::ProcessIntegrand(integrand) => {
                integrand.update_event_processing_runtime(iter)
            }
            #[cfg(test)]
            Integrand::TestProbe(_) => {}
        }
    }

    pub fn observable_accumulator_bundle(&self) -> Option<ObservableAccumulatorBundle> {
        match self {
            Integrand::ProcessIntegrand(integrand) => integrand.observable_accumulator_bundle(),
            #[cfg(test)]
            Integrand::TestProbe(_) => None,
        }
    }

    pub fn has_observables(&self) -> bool {
        match self {
            Integrand::ProcessIntegrand(integrand) => integrand.has_observables(),
            #[cfg(test)]
            Integrand::TestProbe(_) => false,
        }
    }

    pub fn observable_snapshot_bundle(&self) -> Option<ObservableSnapshotBundle> {
        match self {
            Integrand::ProcessIntegrand(integrand) => integrand.observable_snapshot_bundle(),
            #[cfg(test)]
            Integrand::TestProbe(_) => None,
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
            #[cfg(test)]
            Integrand::TestProbe(_) => None,
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
            #[cfg(test)]
            Integrand::TestProbe(_) => Ok(()),
        }
    }
}

impl HasIntegrand for Integrand {
    fn name(&self) -> String {
        match self {
            Integrand::ProcessIntegrand(i) => i.name(),
            #[cfg(test)]
            Integrand::TestProbe(integrand) => integrand.name(),
        }
    }

    fn create_grid(&self) -> Grid<F<f64>> {
        match self {
            Integrand::ProcessIntegrand(integrand) => integrand.create_grid(),
            #[cfg(test)]
            Integrand::TestProbe(integrand) => integrand.create_grid(),
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
            Integrand::ProcessIntegrand(integrand) => {
                integrand.evaluate_sample(sample, model, wgt, iter, use_arb_prec, max_eval)
            }
            #[cfg(test)]
            Integrand::TestProbe(integrand) => {
                integrand.evaluate_sample(sample, model, wgt, iter, use_arb_prec, max_eval)
            }
        }
    }

    fn get_n_dim(&self) -> usize {
        match self {
            Integrand::ProcessIntegrand(integrand) => integrand.get_n_dim(),
            #[cfg(test)]
            Integrand::TestProbe(integrand) => integrand.get_n_dim(),
        }
    }

    fn get_integrator_settings(&self) -> IntegratorSettings {
        match self {
            Integrand::ProcessIntegrand(integrand) => integrand.get_integrator_settings(),
            #[cfg(test)]
            Integrand::TestProbe(integrand) => integrand.get_integrator_settings(),
        }
    }
}

/// Minimal test-only integrand used for exercising the integration state and
/// monitoring machinery. It intentionally has no standalone production owner;
/// process-level acceptance tests cover actual sampling and Jacobians.
#[cfg(test)]
#[derive(Clone)]
pub struct TestProbeIntegrand {
    settings: RuntimeSettings,
    n_dim: usize,
}

#[cfg(test)]
impl TestProbeIntegrand {
    pub(crate) fn new(settings: RuntimeSettings, n_dim: usize) -> Self {
        Self { settings, n_dim }
    }
}

#[cfg(test)]
impl HasIntegrand for TestProbeIntegrand {
    fn create_grid(&self) -> Grid<F<f64>> {
        Grid::Continuous(ContinuousGrid::new(
            self.n_dim,
            self.settings.integrator.n_bins,
            self.settings.integrator.min_samples_for_update,
            self.settings.integrator.bin_number_evolution.clone(),
            self.settings.integrator.train_on_avg,
        ))
    }

    fn name(&self) -> String {
        "TestProbeIntegrand".to_string()
    }

    fn evaluate_sample(
        &mut self,
        _sample: &Sample<F<f64>>,
        _model: &Model,
        wgt: F<f64>,
        _iter: usize,
        _use_arb_prec: bool,
        _max_eval: Complex<F<f64>>,
    ) -> Result<EvaluationResult> {
        let mut result = EvaluationResult::zero();
        result.integrand_result = Complex::new_re(F(1.0));
        result.parameterization_jacobian = Some(F(1.0));
        result.integrator_weight = wgt;
        Ok(result)
    }

    fn get_n_dim(&self) -> usize {
        self.n_dim
    }
}
