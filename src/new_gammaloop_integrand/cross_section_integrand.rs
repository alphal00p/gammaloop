use std::{cell::RefCell, time::Duration};

use itertools::Itertools;
use log::info;
use serde::Serialize;
use spenso::{complex::Complex, ufo::gamma};
use statrs::distribution::Gamma;
use symbolica::{
    domains::float::NumericalFloatLike,
    evaluate::ExpressionEvaluator,
    numerical_integration::{ContinuousGrid, DiscreteGrid, Grid, Sample},
};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cut_expression::{OrientationData, OrientationID},
        esurface::Esurface,
    },
    evaluation_result::{EvaluationMetaData, EvaluationResult},
    graph,
    integrands::HasIntegrand,
    momentum::{Rotation, ThreeMomentum},
    momentum_sample::{ExternalFourMomenta, LoopMomenta, MomentumSample},
    new_cs::{CrossSectionGraph, CutId},
    new_graph::{ExternalConnection, FeynmanGraph, Graph, LmbIndex, LoopMomentumBasis},
    signature::ExternalSignature,
    utils::{self, f128, FloatLike, F},
    DependentMomentaConstructor, DiscreteGraphSamplingSettings, DiscreteGraphSamplingType,
    IntegratedCounterTermRange, IntegratorSettings, MultiChannelingSettings, Polarizations,
    Precision, SamplingSettings, Settings,
};

use super::{
    create_stability_iterator,
    gammaloop_sample::{self, parameterize, DiscreteGraphSample, GammaLoopSample},
    stability_check, GenericEvaluator, GenericEvaluatorFloat, LmbMultiChannelingSetup,
    StabilityLevelResult,
};

const TOLERANCE: F<f64> = F(2.0);

#[derive(Clone)]
pub struct CrossSectionIntegrand {
    pub settings: Settings,
    pub polarizations: Vec<Polarizations>,
    pub rotations: Vec<Rotation>,
    pub graph_terms: Vec<CrossSectionGraphTerm>,
    pub n_incoming: usize,
    pub external_connections: Vec<ExternalConnection>,
}

#[derive(Clone)]
pub struct OrientationEvaluator {
    pub orientation_data: OrientationData,
    pub evaluators: Vec<GenericEvaluator>,
}

#[derive(Clone)]
pub struct CrossSectionGraphTerm {
    pub bare_cff_evaluators: TiVec<CutId, GenericEvaluator>,
    pub bare_cff_orientation_evaluators: TiVec<OrientationID, OrientationEvaluator>,
    pub graph: Graph,
    pub cut_esurface: TiVec<CutId, Esurface>,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
}

impl CrossSectionGraphTerm {
    fn self_get_cuts_to_evaluate(
        &self,
        orientation: Option<OrientationID>,
    ) -> Vec<(CutId, &Esurface)> {
        if let Some(orientation_id) = orientation {
            self.bare_cff_orientation_evaluators[orientation_id]
                .orientation_data
                .cuts
                .iter()
                .map(|cut_id| (*cut_id, &self.cut_esurface[*cut_id]))
                .collect()
        } else {
            self.cut_esurface.iter_enumerated().collect()
        }
    }

    fn evaluate<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        settings: &Settings,
    ) -> F<T> {
        // implementation of forced orientations, only works with sample orientation disabled
        if let Some(forced_orientations) = &settings.general.force_orientations {
            if momentum_sample.sample.orientation.is_none() {
                return forced_orientations
                    .iter()
                    .map(|orientation_usize| {
                        let mut new_sample = momentum_sample.clone();
                        new_sample.sample.orientation =
                            Some(OrientationID::from(*orientation_usize));
                        self.evaluate(&new_sample, settings)
                    })
                    .fold(momentum_sample.zero(), |sum, orientation_result| {
                        sum + orientation_result
                    });
            }
        }

        let center = LoopMomenta::from_iter(vec![
            ThreeMomentum::from([
                momentum_sample.zero(),
                momentum_sample.zero(),
                momentum_sample.zero()
            ]);
            self.graph.underlying.get_loop_number()
        ]);

        let params = self
            .self_get_cuts_to_evaluate(momentum_sample.sample.orientation)
            .into_iter()
            .map(|(cut_id, esurface)| {
                let (tstar_initial, _tstar_initial_negative) = esurface.get_radius_guess(
                    momentum_sample.loop_moms(),
                    momentum_sample.external_moms(),
                    &self.graph.loop_momentum_basis,
                );

                let root_function = |t: &_| {
                    esurface.compute_self_and_r_derivative(
                        t,
                        momentum_sample.loop_moms(),
                        &center,
                        momentum_sample.external_moms(),
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
                        rescaled_sample.loop_moms(),
                        rescaled_sample.external_moms(),
                        &self.graph.loop_momentum_basis,
                    )
                    .into_iter()
                    .map(|(_, x)| x)
                    .collect_vec();

                params.push(newton_result.solution);
                params.push(h_function);
                params.push(newton_result.derivative_at_solution);

                params
            });

        let result = match momentum_sample.sample.orientation {
            Some(orientation_id) => {
                let orientation_evaluator = &self.bare_cff_orientation_evaluators[orientation_id];
                orientation_evaluator
                    .evaluators
                    .iter()
                    .zip(params)
                    .map(|(evaluator, params)| {
                        let cut_results =
                            <T as GenericEvaluatorFloat>::get_evaluator(evaluator)(&params);
                        cut_results
                    })
                    .fold(momentum_sample.zero(), |sum, cut_result| sum + cut_result)
            }
            None => self
                .bare_cff_evaluators
                .iter()
                .zip(params)
                .map(|(evaluator, params)| {
                    let cut_results =
                        <T as GenericEvaluatorFloat>::get_evaluator(evaluator)(&params);
                    cut_results
                })
                .fold(momentum_sample.zero(), |sum, cut_result| sum + cut_result),
        };

        result
    }

    fn create_grid(
        &self,
        settings: &DiscreteGraphSamplingSettings,
        integrator_settings: &IntegratorSettings,
    ) -> Grid<F<f64>> {
        match &settings.sampling_type {
            DiscreteGraphSamplingType::Default(_)
            | DiscreteGraphSamplingType::MultiChanneling(_) => {
                let continuous_grid = self.create_default_continous_grid(integrator_settings);

                if settings.sample_orientations {
                    let continuous_grids = self
                        .bare_cff_orientation_evaluators
                        .iter()
                        .map(|_| Some(continuous_grid.clone()))
                        .collect();

                    Grid::Discrete(DiscreteGrid::new(
                        continuous_grids,
                        integrator_settings.max_prob_ratio,
                        integrator_settings.train_on_avg,
                    ))
                } else {
                    continuous_grid
                }
            }
            DiscreteGraphSamplingType::DiscreteMultiChanneling(_multichanneling_settings) => {
                let continuous_grid = self.create_default_continous_grid(integrator_settings);
                let lmb_channel_grid = Grid::Discrete(DiscreteGrid::new(
                    self.multi_channeling_setup
                        .channels
                        .iter()
                        .map(|_| Some(continuous_grid.clone()))
                        .collect_vec(),
                    integrator_settings.max_prob_ratio,
                    integrator_settings.train_on_avg,
                ));

                if settings.sample_orientations {
                    Grid::Discrete(DiscreteGrid::new(
                        self.bare_cff_orientation_evaluators
                            .iter()
                            .map(|_| Some(lmb_channel_grid.clone()))
                            .collect(),
                        integrator_settings.max_prob_ratio,
                        integrator_settings.train_on_avg,
                    ))
                } else {
                    lmb_channel_grid
                }
            }

            DiscreteGraphSamplingType::TropicalSampling(_) => todo!(),
        }
    }

    fn create_default_continous_grid(
        &self,
        integrator_settings: &IntegratorSettings,
    ) -> Grid<F<f64>> {
        Grid::Continuous(ContinuousGrid::new(
            self.graph.underlying.get_loop_number() * 3,
            integrator_settings.n_bins,
            integrator_settings.min_samples_for_update,
            integrator_settings.bin_number_evolution.clone(),
            integrator_settings.train_on_avg,
        ))
    }
}

impl HasIntegrand for CrossSectionIntegrand {
    fn create_grid(&self) -> Grid<F<f64>> {
        match &self.settings.sampling {
            SamplingSettings::Default(_) => Grid::Continuous(ContinuousGrid::new(
                self.get_n_dim(),
                self.settings.integrator.n_bins,
                self.settings.integrator.min_samples_for_update,
                self.settings.integrator.bin_number_evolution.clone(),
                self.settings.integrator.train_on_avg,
            )),
            SamplingSettings::MultiChanneling(_) => Grid::Continuous(ContinuousGrid::new(
                self.get_n_dim(),
                self.settings.integrator.n_bins,
                self.settings.integrator.min_samples_for_update,
                self.settings.integrator.bin_number_evolution.clone(),
                self.settings.integrator.train_on_avg,
            )),
            SamplingSettings::DiscreteGraphs(discrete_graph_sampling_settings) => {
                Grid::Discrete(DiscreteGrid::new(
                    self.graph_terms
                        .iter()
                        .map(|term| {
                            Some(term.create_grid(
                                discrete_graph_sampling_settings,
                                &self.settings.integrator,
                            ))
                        })
                        .collect(),
                    self.settings.integrator.max_prob_ratio,
                    self.settings.integrator.train_on_avg,
                ))
            }
        }
    }

    fn evaluate_sample(
        &mut self,
        sample: &Sample<F<f64>>,
        wgt: F<f64>,
        iter: usize,
        use_f128: bool,
        max_eval: Complex<F<f64>>,
    ) -> EvaluationResult {
        let start_eval = std::time::Instant::now();
        let stability_iterator = create_stability_iterator(&self.settings.stability, use_f128);
        let dependent_momenta_constructor = DependentMomentaConstructor::CrossSection {
            external_connections: &self.external_connections,
        };

        let mut results_of_stability_levels = Vec::with_capacity(stability_iterator.len());

        for stability_level in stability_iterator.into_iter() {
            let before_parameterization = std::time::Instant::now();
            let ((results, ltd_evaluation_time), parameterization_time) =
                match stability_level.precision {
                    Precision::Double => {
                        if let Ok(gammaloop_sample) = parameterize::<f64>(
                            sample,
                            &self.polarizations,
                            dependent_momenta_constructor,
                            &self.settings,
                            None,
                        ) {
                            let parameterization_time = before_parameterization.elapsed();
                            (
                                self.evaluate_all_rotations(&gammaloop_sample),
                                parameterization_time,
                            )
                        } else {
                            continue;
                        }
                    }
                    Precision::Quad => {
                        if let Ok(gammaloop_sample) = parameterize::<f128>(
                            sample,
                            &self.polarizations,
                            dependent_momenta_constructor,
                            &self.settings,
                            None,
                        ) {
                            let parameterization_time = before_parameterization.elapsed();
                            (
                                self.evaluate_all_rotations(&gammaloop_sample),
                                parameterization_time,
                            )
                        } else {
                            continue;
                        }
                    }
                    Precision::Arb => {
                        todo!()
                    }
                };

            let (average_result, is_stable) =
                stability_check(&self.settings, &results, &stability_level, max_eval, wgt);

            let result_of_level = StabilityLevelResult {
                result: average_result,
                stability_level_used: stability_level.precision,
                parameterization_time,
                ltd_evaluation_time,
                is_stable,
            };

            results_of_stability_levels.push(result_of_level);

            if is_stable {
                break;
            }
        }

        if let Some(stability_level_result) = results_of_stability_levels.last() {
            let re_is_nan = stability_level_result.result.re.is_nan();
            let im_is_nan = stability_level_result.result.im.is_nan();
            let is_nan = re_is_nan || im_is_nan;

            let meta_data = EvaluationMetaData {
                total_timing: start_eval.elapsed(),
                rep3d_evaluation_time: stability_level_result.ltd_evaluation_time,
                highest_precision: stability_level_result.stability_level_used,
                parameterization_time: stability_level_result.parameterization_time,
                relative_instability_error: Complex::new_zero(),
                is_nan,
            };

            let nanless_result = if re_is_nan && !im_is_nan {
                Complex::new(F(0.0), stability_level_result.result.im)
            } else if im_is_nan && !re_is_nan {
                Complex::new(stability_level_result.result.re, F(0.0))
            } else if im_is_nan && re_is_nan {
                Complex::new(F(0.0), F(0.0))
            } else {
                stability_level_result.result
            };

            EvaluationResult {
                integrand_result: nanless_result,
                integrator_weight: wgt,
                event_buffer: vec![],
                evaluation_metadata: meta_data,
            }
        } else {
            // this happens if parameterization fails at all levels
            println!("Parameterization failed at all levels");
            EvaluationResult {
                integrand_result: Complex::new(F(0.0), F(0.0)),
                integrator_weight: wgt,
                event_buffer: vec![],
                evaluation_metadata: EvaluationMetaData {
                    total_timing: Duration::ZERO,
                    rep3d_evaluation_time: Duration::ZERO,
                    parameterization_time: Duration::ZERO,
                    relative_instability_error: Complex::new(F(0.0), F(0.0)),
                    is_nan: true,
                    highest_precision: crate::Precision::Double,
                },
            }
        }
    }

    fn get_n_dim(&self) -> usize {
        assert!(self
            .graph_terms
            .iter()
            .map(|term| term.graph.underlying.get_loop_number())
            .all_equal());

        self.graph_terms[0].graph.underlying.get_loop_number() * 3
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

impl CrossSectionIntegrand {
    fn evaluate_all_rotations<T: FloatLike>(
        &self,
        gammaloop_sample: &GammaLoopSample<T>,
    ) -> (Vec<Complex<F<f64>>>, Duration) {
        // rotate the momenta for the stability tests.
        let gammaloop_samples: Vec<_> = self
            .rotations
            .iter()
            .map(|rotation| gammaloop_sample.get_rotated_sample(rotation))
            .collect();

        let start_time = std::time::Instant::now();
        let evaluation_results = gammaloop_samples
            .iter()
            .map(|gammaloop_sample| self.evaluate_single_rotation(gammaloop_sample))
            .collect_vec();
        let duration = start_time.elapsed() / gammaloop_samples.len() as u32;

        (evaluation_results, duration)
    }

    fn evaluate_single_rotation<T: FloatLike>(
        &self,
        gammaloop_sample: &GammaLoopSample<T>,
    ) -> Complex<F<f64>> {
        let result = match &gammaloop_sample {
            GammaLoopSample::Default(sample) => self
                .graph_terms
                .iter()
                .map(|term| term.evaluate(sample, &self.settings))
                .fold(gammaloop_sample.get_default_sample().zero(), |sum, term| {
                    sum + term
                }),
            GammaLoopSample::MultiChanneling { alpha, sample } => self
                .graph_terms
                .iter()
                .map(|term| {
                    let channels_samples = term
                        .multi_channeling_setup
                        .reinterpret_loop_momenta_and_compute_prefactor_all_channels(
                            sample,
                            &term.graph,
                            &term.lmbs,
                            alpha,
                        );

                    channels_samples
                        .into_iter()
                        .map(|(reparameterized_sample, prefactor)| {
                            prefactor * term.evaluate(&reparameterized_sample, &self.settings)
                        })
                        .fold(gammaloop_sample.get_default_sample().zero(), |sum, term| {
                            sum + term
                        })
                })
                .fold(gammaloop_sample.get_default_sample().zero(), |sum, term| {
                    sum + term
                }),
            GammaLoopSample::DiscreteGraph { graph_id, sample } => {
                let graph_term = &self.graph_terms[*graph_id];
                match sample {
                    DiscreteGraphSample::Default(sample) => {
                        graph_term.evaluate(sample, &self.settings)
                    }
                    DiscreteGraphSample::DiscreteMultiChanneling {
                        alpha,
                        channel_id,
                        sample,
                    } => {
                        let (reparameterized_sample, prefactor) = self.graph_terms[*graph_id]
                            .multi_channeling_setup
                            .reinterpret_loop_momenta_and_compute_prefactor(
                                *channel_id,
                                sample,
                                &graph_term.graph,
                                &graph_term.lmbs,
                                alpha,
                            );

                        prefactor * graph_term.evaluate(&reparameterized_sample, &self.settings)
                    }
                    DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                        let channel_samples = self.graph_terms[*graph_id]
                            .multi_channeling_setup
                            .reinterpret_loop_momenta_and_compute_prefactor_all_channels(
                                sample,
                                &graph_term.graph,
                                &graph_term.lmbs,
                                alpha,
                            );

                        channel_samples
                            .into_iter()
                            .map(|(reparameterized_sample, prefactor)| {
                                prefactor
                                    * graph_term.evaluate(&reparameterized_sample, &self.settings)
                            })
                            .fold(gammaloop_sample.get_default_sample().zero(), |sum, term| {
                                sum + term
                            })
                    }
                    _ => todo!(),
                }
            }
            _ => todo!(),
        };

        let res = result * gammaloop_sample.get_default_sample().jacobian();

        Complex::new_re(res.into_ff64())
    }
}
