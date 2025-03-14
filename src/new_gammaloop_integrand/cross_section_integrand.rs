use std::time::Duration;

use itertools::Itertools;
use log::info;
use serde::Serialize;
use spenso::{complex::Complex, ufo::gamma};
use symbolica::{
    evaluate::ExpressionEvaluator,
    numerical_integration::{ContinuousGrid, Grid, Sample},
};
use typed_index_collections::TiVec;

use crate::{
    cff::esurface::Esurface,
    evaluation_result::{EvaluationMetaData, EvaluationResult},
    graph,
    integrands::HasIntegrand,
    momentum::ThreeMomentum,
    momentum_sample::{ExternalFourMomenta, LoopMomenta, MomentumSample},
    new_cs::{CrossSectionGraph, CutId},
    new_graph::{FeynmanGraph, Graph},
    signature::ExternalSignature,
    utils::{self, FloatLike, F},
    DependentMomentaConstructor, IntegratedCounterTermRange, Polarizations, Settings,
};

use super::gammaloop_sample::{self, parameterize, DiscreteGraphSample, GammaLoopSample};
const TOLERANCE: F<f64> = F(2.0);

#[derive(Clone)]
pub struct CrossSectionIntegrand {
    pub settings: Settings,
    pub polarizations: Vec<Polarizations>,
    pub graph_terms: Vec<CrossSectionGraphTerm>,
    pub n_incoming: usize,
}

#[derive(Clone)]
pub struct CrossSectionGraphTerm {
    pub bare_cff_evaluators: TiVec<CutId, ExpressionEvaluator<F<f64>>>,
    pub graph: Graph,
    pub cut_esurface: TiVec<CutId, Esurface>,
}

impl CrossSectionGraphTerm {
    fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        settings: &Settings,
    ) -> F<T> {
        let mut result = momentum_sample.zero();

        let center = LoopMomenta::from_iter(vec![
            ThreeMomentum::from([
                momentum_sample.zero(),
                momentum_sample.zero(),
                momentum_sample.zero()
            ]);
            self.graph.underlying.get_loop_number()
        ]);

        for (cut_id, esurface) in self.cut_esurface.iter_enumerated() {
            let (tstar_initial, _tstar_initial_negative) = esurface.get_radius_guess(
                &momentum_sample.loop_moms(),
                &momentum_sample.external_moms(),
                &self.graph.loop_momentum_basis,
            );

            let root_function = |t: &_| {
                esurface.compute_self_and_r_derivative(
                    t,
                    &momentum_sample.loop_moms(),
                    &center,
                    &momentum_sample.external_moms(),
                    &self.graph.loop_momentum_basis,
                    &self.graph.underlying.get_real_mass_vector(),
                )
            };

            let e_cm = F::from_ff64(settings.kinematics.e_cm);
            let tolerance = F::from_ff64(TOLERANCE);

            let newton_result = newton_iteration_and_derivative(
                &tstar_initial,
                root_function,
                &tolerance,
                20,
                &e_cm,
            );

            let rescaled_sample =
                momentum_sample.rescaled_loop_momenta(&newton_result.solution, None);

            let h_function_settings = match &settings.subtraction.integrated_ct_settings.range {
                IntegratedCounterTermRange::Compact => panic!(),
                IntegratedCounterTermRange::Infinite {
                    h_function_settings,
                } => h_function_settings,
            };

            let h_function = utils::h(&newton_result.solution, None, None, h_function_settings);

            let mut params = self
                .graph
                .underlying
                .get_energy_cache(
                    &rescaled_sample.loop_moms(),
                    &rescaled_sample.external_moms(),
                    &self.graph.loop_momentum_basis,
                )
                .into_iter()
                .map(|(_, x)| x.into_ff64())
                .collect_vec();

            params.push(newton_result.solution.into_ff64());
            params.push(h_function.into_ff64());
            params.push(newton_result.derivative_at_solution.into_ff64());

            let mut out = [F(0.0)];
            self.bare_cff_evaluators[cut_id].evaluate(&params, &mut out);

            let res = out[0];
            result += F::from_ff64(res);
        }
        result
    }
}

impl HasIntegrand for CrossSectionIntegrand {
    fn create_grid(&self) -> Grid<F<f64>> {
        Grid::Continuous(ContinuousGrid::new(
            3,
            self.settings.integrator.n_bins,
            self.settings.integrator.min_samples_for_update,
            self.settings.integrator.bin_number_evolution.clone(),
            self.settings.integrator.train_on_avg,
        ))
    }

    fn evaluate_sample(
        &mut self,
        sample: &Sample<F<f64>>,
        wgt: F<f64>,
        iter: usize,
        use_f128: bool,
        max_eval: F<f64>,
    ) -> EvaluationResult {
        let gammaloop_sample = parameterize(
            sample,
            &self.polarizations,
            DependentMomentaConstructor::CrossSection {
                n_incoming: self.n_incoming,
            },
            &self.settings,
            None,
        )
        .unwrap();

        let result = match &gammaloop_sample {
            GammaLoopSample::Default(sample) => {
                let mut res = F(0.0);
                for term in &mut self.graph_terms {
                    res += term.evaluate(sample, &self.settings);
                }
                res
            }
            GammaLoopSample::DiscreteGraph { graph_id, sample } => match sample {
                DiscreteGraphSample::Default(sample) => {
                    self.graph_terms[*graph_id].evaluate(sample, &self.settings)
                }
                _ => todo!(),
            },
            _ => todo!(),
        };

        let mut res = result * gammaloop_sample.get_default_sample().jacobian();

        let is_nan = if res.is_nan() {
            res = F(0.0);
            res.is_nan()
        } else {
            res.is_nan()
        };

        EvaluationResult {
            integrand_result: Complex::new(res, F(0.0)),
            integrator_weight: wgt,
            event_buffer: vec![],
            evaluation_metadata: EvaluationMetaData {
                total_timing: Duration::ZERO,
                rep3d_evaluation_time: Duration::ZERO,
                parameterization_time: Duration::ZERO,
                relative_instability_error: Complex::new(F(0.0), F(0.0)),
                is_nan,
                highest_precision: crate::Precision::Double,
            },
        }
    }

    fn get_n_dim(&self) -> usize {
        3
    }
}

/// root finding, returns the derivative at the root, so that we don't have to recompute it.
/// Also returns the value of the function whose root is being found and the number of iterations used for debug information
fn newton_iteration_and_derivative<T: FloatLike>(
    guess: &F<T>,
    f_x_and_df_x: impl Fn(&F<T>) -> (F<T>, F<T>),
    tolerance: &F<T>,
    max_iterations: usize,
    e_cm: &F<T>,
) -> NewtonIterationResult<T> {
    let mut x = guess.clone();
    let (mut val_f_x, mut val_df_x) = f_x_and_df_x(&x);

    let mut iteration = 0;

    while iteration < max_iterations && val_f_x.abs() > guess.epsilon() * tolerance * e_cm {
        x -= val_f_x / val_df_x;
        (val_f_x, val_df_x) = f_x_and_df_x(&x);
        iteration += 1;
    }

    NewtonIterationResult {
        solution: x,
        derivative_at_solution: val_df_x,
        error_of_function: val_f_x,
        num_iterations_used: iteration,
    }
}

#[derive(Serialize, Clone)]
struct NewtonIterationResult<T: FloatLike> {
    solution: F<T>,
    derivative_at_solution: F<T>,
    error_of_function: F<T>,
    num_iterations_used: usize,
}
