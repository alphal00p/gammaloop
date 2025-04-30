use std::time::Duration;

use itertools::Itertools;
use serde::Serialize;
use spenso::complex::Complex;
use symbolica::numerical_integration::{Grid, Sample};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cut_expression::{CutOrientationData, OrientationID},
        esurface::Esurface,
    },
    evaluation_result::{EvaluationMetaData, EvaluationResult},
    integrands::HasIntegrand,
    momentum::{Rotation, ThreeMomentum},
    momentum_sample::{LoopMomenta, MomentumSample},
    new_cs::CutId,
    new_graph::{ExternalConnection, FeynmanGraph, Graph, LmbIndex, LoopMomentumBasis},
    utils::{self, f128, FloatLike, F},
    DependentMomentaConstructor, IntegratedCounterTermRange, Polarizations, Precision, Settings,
};

use super::{
    create_grid, create_stability_iterator, evaluate_all_rotations, gammaloop_sample::parameterize,
    stability_check, GammaloopIntegrand, GenericEvaluator, GenericEvaluatorFloat, GraphTerm,
    LmbMultiChannelingSetup, StabilityLevelResult,
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

impl GammaloopIntegrand for CrossSectionIntegrand {
    type G = CrossSectionGraphTerm;
    fn get_rotations(&self) -> impl Iterator<Item = &Rotation> {
        self.rotations.iter()
    }

    fn get_terms(&self) -> impl Iterator<Item = &Self::G> {
        self.graph_terms.iter()
    }

    fn get_settings(&self) -> &Settings {
        &self.settings
    }

    fn get_graph(&self, graph_id: usize) -> &Self::G {
        &self.graph_terms[graph_id]
    }

    fn get_polarizations(&self) -> &[Polarizations] {
        &self.polarizations
    }

    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor {
        DependentMomentaConstructor::CrossSection {
            external_connections: &self.external_connections,
        }
    }
}

#[derive(Clone)]
pub struct OrientationEvaluator {
    pub orientation_data: CutOrientationData,
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

impl GraphTerm for CrossSectionGraphTerm {
    fn evaluate<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        settings: &Settings,
    ) -> Complex<F<T>> {
        self.evaluate(momentum_sample, settings)
    }

    fn get_multi_channeling_setup(&self) -> &LmbMultiChannelingSetup {
        &self.multi_channeling_setup
    }

    fn get_graph(&self) -> &Graph {
        &self.graph
    }

    fn get_lmbs(&self) -> &TiVec<LmbIndex, LoopMomentumBasis> {
        &self.lmbs
    }

    fn get_num_orientations(&self) -> usize {
        self.bare_cff_orientation_evaluators.len()
    }
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
    ) -> Complex<F<T>> {
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
                    .fold(
                        Complex::new_re(momentum_sample.zero()),
                        |sum, orientation_result| sum + orientation_result,
                    );
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
            .map(|(_cut_id, esurface)| {
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

        Complex::new_re(result)
    }
}

impl HasIntegrand for CrossSectionIntegrand {
    fn create_grid(&self) -> Grid<F<f64>> {
        create_grid(self)
    }

    fn evaluate_sample(
        &mut self,
        sample: &Sample<F<f64>>,
        wgt: F<f64>,
        _iter: usize,
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
                                evaluate_all_rotations(self, &gammaloop_sample),
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
                                evaluate_all_rotations(self, &gammaloop_sample),
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
