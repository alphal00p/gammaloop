//! This module contains the GammaloopIntegrand struct, which represents the integrand for physical
//! amplitudes and Local Unitarity crosssections.

use core::panic;
use std::fmt::Display;
use std::time::{Duration, Instant};

use crate::cross_section::{Amplitude, AmplitudeGraph, CrossSection, IsPolarizable, SuperGraph};
use crate::debug_info::{EvalState, DEBUG_LOGGER};
use crate::evaluation_result::{EvaluationMetaData, EvaluationResult};
use crate::graph::{EdgeType, Graph, LoopMomentumBasisSpecification, SerializableGraph};
use crate::integrands::{HasIntegrand, Integrand};
use crate::integrate::UserData;
use crate::momentum::{
    FourMomentum, Polarization, Rotatable, Rotation, RotationMethod, Signature, ThreeMomentum,
};
use crate::numerator::Evaluators;
use crate::subtraction::static_counterterm::CounterTerm;
use crate::utils::{
    self, format_for_compare_digits, get_n_dim_for_n_loop_momenta, global_parameterize, FloatLike,
    PrecisionUpgradable, F,
};
use crate::{
    DiscreteGraphSamplingSettings, Externals, IntegratedPhase, Polarizations, SamplingSettings,
    Settings,
};
use crate::{Precision, StabilityLevelSetting};
use colored::Colorize;
use itertools::Itertools;
use momtrop::vector::Vector;
use serde::{Deserialize, Serialize};
use spenso::complex::Complex;
use spenso::contraction::IsZero;
use symbolica::domains::float::{NumericalFloatLike, Real};
use symbolica::numerical_integration::{ContinuousGrid, DiscreteGrid, Grid, Sample};

/// Trait to capture the common behaviour of amplitudes and cross sections
/// Mainly used to expose the properties of the underlying graph in both amplitudes and cross sections
trait GraphIntegrand {
    /// Get the underlying graph
    fn get_graph(&self) -> &Graph;

    fn get_mut_graph(&mut self) -> &mut Graph;

    /// Get the channels used for multi channeling
    fn get_multi_channeling_channels(&self) -> &[usize];

    /// Most basic form of evaluating the 3D representation of the underlying loop integral
    fn evaluate<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        rotation_for_overlap: &Rotation,
        settings: &Settings,
    ) -> Complex<F<T>>;

    /// Evaluate in a single LMB-channel
    fn evaluate_channel<T: FloatLike>(
        &mut self,
        channel_id: usize,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<F<T>>;

    /// Evaluate a sum over LMB-channels
    fn evaluate_channel_sum<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<F<T>>;

    /// Evaluate to use when tropical sampling, raises the power of the onshell energies in front of the 3D representation according to the chosen weights.
    fn evaluate_tropical<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        rotation_for_overlap: &Rotation,
        settings: &Settings,
    ) -> Complex<F<T>>;
}

impl GraphIntegrand for AmplitudeGraph<Evaluators> {
    fn get_graph(&self) -> &Graph {
        &self.graph
    }

    fn get_mut_graph(&mut self) -> &mut Graph {
        &mut self.graph
    }

    fn get_multi_channeling_channels(&self) -> &[usize] {
        &self.multi_channeling_channels
    }

    #[inline]
    #[allow(unused_variables)]
    fn evaluate_channel<T: FloatLike>(
        &mut self,
        channel_id: usize,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<F<T>> {
        if settings.general.debug > 0 {
            DEBUG_LOGGER.write("channel_id", &channel_id);
            DEBUG_LOGGER.write(
                "graph",
                &SerializableGraph::from_graph(&self.get_graph().bare_graph),
            );
        }

        let one = sample.one();
        let zero = sample.zero();
        // reference to the list of all lmbs
        let lmb_list = self
            .get_graph()
            .derived_data
            .as_ref()
            .unwrap()
            .loop_momentum_bases
            .clone()
            .unwrap();

        // if the channel list is empty, we use all channels.
        let channels = if !self.get_multi_channeling_channels().is_empty() {
            self.get_multi_channeling_channels().to_vec()
        } else {
            (0..lmb_list.len()).collect_vec()
        };

        //map the channel_id to the corresponding lmb
        let channel = channels[channel_id];
        let lmb_specification = LoopMomentumBasisSpecification::FromList(channel);
        let lmb = &lmb_list[channel]; // lmb_specification.basis(self.get_graph());

        let channels_lmbs = channels.iter().map(|&i| &lmb_list[i]);

        let rep3d = if settings.general.use_ltd {
            self.get_mut_graph()
                .evaluate_ltd_expression_in_lmb(sample, lmb, settings)
        } else {
            self.get_mut_graph().evaluate_cff_expression_in_lmb(
                sample,
                &lmb_specification,
                settings,
            )
        };

        let onshell_energies = self.get_graph().bare_graph.compute_onshell_energies_in_lmb(
            sample.loop_moms(),
            sample.external_moms(),
            lmb,
        );

        let virtual_energies = self
            .get_graph()
            .bare_graph
            .edges
            .iter()
            .enumerate()
            .filter(|(_, e)| e.edge_type == EdgeType::Virtual)
            .map(|(i, _)| i);

        let energy_product = virtual_energies
            .map(|i| one.from_i64(2) * &onshell_energies[i])
            .fold(one.clone(), |acc, x| acc * x);

        let lmb_products = channels_lmbs.map(|basis| {
            basis
                .basis
                .iter()
                .fold(one.clone(), |acc, i| acc * &onshell_energies[*i])
        });

        let denominator = lmb_products
            .map(|x| x.powf(&-F::<T>::from_f64(alpha)) * &energy_product)
            .reduce(|acc, e| acc + &e)
            .unwrap_or(zero.clone());

        let multichanneling_numerator = Complex::new(
            lmb_list[channel]
                .basis
                .iter()
                .fold(one.clone(), |acc, i| acc * &onshell_energies[*i])
                .powf(&-F::<T>::from_f64(alpha)),
            sample.zero(),
        );

        let prefactor = if let Some(p) = settings.general.amplitude_prefactor {
            p.map(|x| F::from_ff64(x))
        } else {
            Complex::new(one, zero)
        };

        let multichanneling_prefactor = &multichanneling_numerator / &denominator;

        if settings.general.debug > 0 {
            let multichanneling_prefactor_str =
                serde_json::to_string(&multichanneling_prefactor).unwrap();
            let rep3d_str = serde_json::to_string(&rep3d).unwrap();

            DEBUG_LOGGER.write("multichanneling_prefactor", &multichanneling_prefactor);
            DEBUG_LOGGER.write("rep3d", &rep3d);
        }

        multichanneling_prefactor * rep3d * prefactor
    }

    #[inline]
    fn evaluate_channel_sum<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<F<T>> {
        let zero = sample.zero();

        // a bit annoying that this is duplicated from evaluate_channel
        let lmb_list = self
            .get_graph()
            .derived_data
            .as_ref()
            .unwrap()
            .loop_momentum_bases
            .as_ref()
            .unwrap();

        let channels = if !self.get_multi_channeling_channels().is_empty() {
            self.get_multi_channeling_channels().to_vec()
        } else {
            (0..lmb_list.len()).collect_vec()
        };

        channels
            .iter()
            .map(|&i| self.evaluate_channel(i, sample, alpha, settings))
            .reduce(|acc, e| acc + &e)
            .unwrap_or(zero.clone().into())
    }

    #[inline]
    fn evaluate<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        rotation_for_overlap: &Rotation,
        settings: &Settings,
    ) -> Complex<F<T>> {
        let zero_builder = &sample.zero();

        let rep3d = if settings.general.use_ltd {
            self.get_mut_graph()
                .evaluate_ltd_expression(sample, settings)
        } else {
            self.get_mut_graph()
                .evaluate_cff_expression(sample, settings)
        };

        let energy_product = self
            .get_graph()
            .bare_graph
            .compute_energy_product(sample.loop_moms(), sample.external_moms());

        let counter_term_eval = self.get_mut_graph().evaluate_threshold_counterterm(
            sample,
            rotation_for_overlap,
            settings,
        );

        let prefactor = if let Some(p) = settings.general.amplitude_prefactor {
            p.map(|x| F::from_ff64(x))
        } else {
            Complex::new(zero_builder.one(), zero_builder.zero())
        };

        if settings.general.debug > 0 {
            DEBUG_LOGGER.write(
                "graph",
                &SerializableGraph::from_graph(&self.get_graph().bare_graph),
            );
            DEBUG_LOGGER.write("rep3d", &(&rep3d / &energy_product));
            DEBUG_LOGGER.write("ose_product", &energy_product);
            DEBUG_LOGGER.write("counter_terms", &counter_term_eval);
        }

        prefactor * (rep3d / energy_product + counter_term_eval)
    }

    #[inline]
    fn evaluate_tropical<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        rotation_for_overlap: &Rotation,
        settings: &Settings,
    ) -> Complex<F<T>> where {
        let one = sample.one();

        let rep3d = if settings.general.use_ltd {
            self.get_mut_graph()
                .evaluate_ltd_expression(sample, settings)
        } else {
            self.get_mut_graph()
                .evaluate_cff_expression(sample, settings)
        };

        let onshell_energies = self
            .get_graph()
            .bare_graph
            .compute_onshell_energies(sample.loop_moms(), sample.external_moms());

        let tropical_subgraph_table = self.get_graph().get_tropical_subgraph_table();

        let virtual_loop_energies = self
            .get_graph()
            .bare_graph
            .get_loop_edges_iterator()
            .map(|(index, _)| onshell_energies[index].clone());

        let weight_iterator = tropical_subgraph_table.iter_edge_weights();

        let energy_product = virtual_loop_energies
            .zip(weight_iterator)
            .map(|(energy, weight)| energy.powf(&F::<T>::from_f64(2. * weight - 1.)))
            .fold(one.clone(), |acc, x| acc * x);

        let tree_like_energies = self
            .get_graph()
            .bare_graph
            .get_tree_level_edges_iterator()
            .map(|(index, _)| onshell_energies[index].clone());

        let tree_product =
            tree_like_energies.fold(one.clone(), |acc, x| acc * F::<T>::from_f64(2.) * x);

        let counterterm = self.get_mut_graph().evaluate_threshold_counterterm(
            sample,
            rotation_for_overlap,
            settings,
        ) * self
            .get_graph()
            .bare_graph
            .compute_energy_product(sample.loop_moms(), sample.external_moms());

        let final_energy_product = &energy_product / &tree_product;

        let prefactor = if let Some(p) = settings.general.amplitude_prefactor {
            p.map(|x| F::from_ff64(x))
        } else {
            Complex::new(sample.one(), sample.zero())
        };

        if settings.general.debug > 0 {
            DEBUG_LOGGER.write("rep3d", &rep3d);
            DEBUG_LOGGER.write("ose_product", &energy_product);
            DEBUG_LOGGER.write("counter_terms", &counterterm);
            DEBUG_LOGGER.write(
                "graph",
                &SerializableGraph::from_graph(&self.get_graph().bare_graph),
            );
        }

        (rep3d + counterterm) * final_energy_product * prefactor
    }
}

impl GraphIntegrand for SuperGraph {
    fn get_graph(&self) -> &Graph {
        &self.graph
    }

    fn get_mut_graph(&mut self) -> &mut Graph {
        &mut self.graph
    }

    fn get_multi_channeling_channels(&self) -> &[usize] {
        todo!()
    }

    #[allow(unused)]
    fn evaluate_channel<T: FloatLike>(
        &mut self,
        channel_id: usize,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<F<T>> {
        // sum over cuts
        todo!()
    }

    #[allow(unused)]
    fn evaluate_channel_sum<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<F<T>> {
        // sum over channels
        todo!()
    }

    #[allow(unused)]
    #[inline]
    fn evaluate<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        rotate_overlap_centers: &Rotation,
        settings: &Settings,
    ) -> Complex<F<T>> {
        // sum over channels
        todo!()
    }

    #[allow(unused)]
    #[inline]
    fn evaluate_tropical<T: FloatLike>(
        &mut self,
        sample: &DefaultSample<T>,
        rotate_overlap_centers: &Rotation,
        settings: &Settings,
    ) -> Complex<F<T>> {
        // sum over channels
        todo!()
    }
}

/// Get the number of different loop momentum bases (number of spanning trees)
fn get_lmb_count<T: GraphIntegrand>(graph_integrand: &T) -> usize {
    graph_integrand
        .get_graph()
        .derived_data
        .as_ref()
        .unwrap()
        .loop_momentum_bases
        .as_ref()
        .unwrap_or_else(|| panic!("Loop momentum bases not generated"))
        .len()
}

fn get_loop_count<T: GraphIntegrand>(graph_integrand: &T) -> usize {
    graph_integrand
        .get_graph()
        .bare_graph
        .loop_momentum_basis
        .basis
        .len()
}

/// Evaluate the sample correctly according to the sample type
#[inline]
fn evaluate<I: GraphIntegrand, T: FloatLike>(
    graph_integrands: &mut [I],
    sample: &GammaLoopSample<T>,
    rotation_for_overlap: &Rotation,
    settings: &Settings,
) -> Complex<F<T>> {
    if settings.general.debug > 0 {
        DEBUG_LOGGER.write(
            "momenta_sample",
            sample.get_default_sample().possibly_rotated_sample(),
        );
    }

    let zero = sample.zero();
    match sample {
        GammaLoopSample::Default(sample) => graph_integrands
            .iter_mut()
            .map(|g| g.evaluate(sample, rotation_for_overlap, settings))
            .reduce(|acc, e| acc + &e)
            .unwrap_or(zero.clone().into()),
        GammaLoopSample::MultiChanneling { alpha, sample } => graph_integrands
            .iter_mut()
            .map(|g| g.evaluate_channel_sum(sample, *alpha, settings))
            .reduce(|acc, e| acc + &e)
            .unwrap_or(zero.clone().into()),
        GammaLoopSample::DiscreteGraph { graph_id, sample } => {
            let graph = &mut graph_integrands[*graph_id];
            match sample {
                DiscreteGraphSample::Default(sample) => {
                    graph.evaluate(sample, rotation_for_overlap, settings)
                }
                DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                    graph.evaluate_channel_sum(sample, *alpha, settings)
                }
                DiscreteGraphSample::Tropical(sample) => {
                    graph.evaluate_tropical(sample, rotation_for_overlap, settings)
                }
                DiscreteGraphSample::DiscreteMultiChanneling {
                    alpha,
                    channel_id,
                    sample,
                } => graph.evaluate_channel(*channel_id, sample, *alpha, settings),
            }
        }
    }
}

/// Create a havana grid for a single graph
fn create_grid<T: GraphIntegrand>(graph_integrand: &T, settings: &Settings) -> Grid<F<f64>> {
    let num_loops = get_loop_count(graph_integrand);

    let n_edges = graph_integrand
        .get_graph()
        .derived_data
        .as_ref()
        .unwrap()
        .tropical_subgraph_table
        .as_ref()
        .map(|t| t.get_num_edges());

    let continious_dimension = get_n_dim_for_n_loop_momenta(settings, num_loops, false, n_edges);

    let continous_grid = Grid::Continuous(ContinuousGrid::new(
        continious_dimension,
        settings.integrator.n_bins,
        settings.integrator.min_samples_for_update,
        settings.integrator.bin_number_evolution.clone(),
        settings.integrator.train_on_avg,
    ));

    match &settings.sampling {
        SamplingSettings::DiscreteGraphs(
            DiscreteGraphSamplingSettings::DiscreteMultiChanneling(_),
        ) => {
            // if the channel list is empty, we use all channels.
            let num_channels = if graph_integrand.get_multi_channeling_channels().is_empty() {
                get_lmb_count(graph_integrand)
            } else {
                graph_integrand.get_multi_channeling_channels().len()
            };

            let cont_grids = vec![Some(continous_grid); num_channels];

            Grid::Discrete(DiscreteGrid::new(
                cont_grids,
                settings.integrator.max_prob_ratio,
                settings.integrator.train_on_avg,
            ))
        }
        _ => continous_grid,
    }
}

/// Struct that represents a list of graph contirbuting to a single amplitude or cross-section.
#[derive(Clone)]
enum GraphIntegrands {
    Amplitude(Vec<AmplitudeGraph<Evaluators>>),
    CrossSection(Vec<SuperGraph>),
}

/// GammaloopIntegrand contains a list of graphs and the settings.
#[derive(Clone)]
pub struct GammaLoopIntegrand {
    pub global_data: GlobalData,
    graph_integrands: GraphIntegrands,
}

#[derive(Clone)]
pub struct GlobalData {
    pub polarizations: Vec<Polarizations>,
    pub externals: Vec<Externals>,
    pub rotations: Vec<Rotation>,
    pub settings: Settings,
}

impl GraphIntegrands {
    /// Create a havana grid for a list of graphs
    fn create_grid(&self, settings: &Settings) -> Grid<F<f64>> {
        match settings.sampling {
            SamplingSettings::DiscreteGraphs(_) => match self {
                GraphIntegrands::Amplitude(graphs) => {
                    let grids = graphs
                        .iter()
                        .map(|g| Some(create_grid(g, settings)))
                        .collect();

                    Grid::Discrete(DiscreteGrid::new(
                        grids,
                        settings.integrator.max_prob_ratio,
                        settings.integrator.train_on_avg,
                    ))
                }
                GraphIntegrands::CrossSection(graphs) => {
                    let grids = graphs
                        .iter()
                        .map(|g| Some(create_grid(g, settings)))
                        .collect();

                    Grid::Discrete(DiscreteGrid::new(
                        grids,
                        settings.integrator.max_prob_ratio,
                        settings.integrator.train_on_avg,
                    ))
                }
            },
            _ => match self {
                GraphIntegrands::Amplitude(graphs) => create_grid(&graphs[0], settings),
                GraphIntegrands::CrossSection(graphs) => create_grid(&graphs[0], settings),
            },
        }
    }
}

impl HasIntegrand for GammaLoopIntegrand {
    fn create_grid(&self) -> Grid<F<f64>> {
        self.graph_integrands
            .create_grid(&self.global_data.settings)
    }

    #[allow(unused_variables)]
    fn evaluate_sample(
        &mut self,
        sample: &symbolica::numerical_integration::Sample<F<f64>>,
        wgt: F<f64>,
        iter: usize,
        use_f128: bool,
        max_eval: F<f64>,
    ) -> EvaluationResult {
        if self.global_data.settings.general.debug > 0 {
            DEBUG_LOGGER.write("new_evaluation", &());
            DEBUG_LOGGER.write("havana_sample", sample);
        }

        let start_evaluate_sample = Instant::now();

        // setup the evaluation of the integrand in the different stability levels
        let mut results_of_stability_levels =
            Vec::with_capacity(self.global_data.settings.stability.levels.len());

        // create an iterator containing the information for evaluation at each stability level
        let stability_iterator = self.create_stability_vec(use_f128);

        let before_parameterization = std::time::Instant::now();
        let sample_point_result = self.parameterize(sample);
        let sample_point = match sample_point_result {
            Ok(sample_point) => sample_point,
            Err(_) => {
                return EvaluationResult {
                    integrand_result: Complex::new_zero(),
                    integrator_weight: F(0.0),
                    event_buffer: vec![],
                    evaluation_metadata: EvaluationMetaData {
                        total_timing: Duration::ZERO,
                        rep3d_evaluation_time: Duration::ZERO,
                        parameterization_time: Duration::ZERO,
                        relative_instability_error: Complex::new_zero(),
                        highest_precision: Precision::Double,
                        is_nan: true,
                    },
                };
            }
        };

        if self.global_data.settings.general.debug > 0 {
            DEBUG_LOGGER.write("jacobian", &sample_point.get_default_sample().jacobian());
        }

        let parameterization_time = before_parameterization.elapsed();

        // rotate the momenta for the stability tests.
        let samples: Vec<_> = self
            .global_data
            .rotations
            .iter()
            .zip(self.global_data.externals.iter().cloned())
            .zip(self.global_data.polarizations.iter().cloned())
            .map(|((r, e), p)| (sample_point.get_rotated_sample(r, e, p)))
            .collect(); //since the first rotation is the identity, the first sample is the same as the original sample

        // 1 / (2 pi )^L
        let prefactor = F(self.compute_2pi_factor().inv());

        if self.global_data.settings.general.debug > 0 {
            DEBUG_LOGGER.write("pi_prefactor", &prefactor.0);
        };

        // iterate over the stability levels, break if the point is stable
        for stability_level in &stability_iterator {
            // evaluate the integrand at the current stability level
            let (results, duration) = self.evaluate_at_prec(&samples, stability_level.precision);
            let results_scaled = results
                .iter()
                .zip(samples.iter())
                .map(|(result, sample)| result * sample.get_default_sample().jacobian() * prefactor)
                .collect_vec();

            // check for the stability
            let (avg_result, stable) = self.stability_check(
                &results_scaled,
                stability_level,
                self.global_data.settings.integrator.integrated_phase,
                max_eval,
                wgt,
            );

            results_of_stability_levels.push((
                avg_result,
                stable,
                stability_level.precision,
                duration,
            ));

            if stable {
                break;
            }
        }

        if self.global_data.settings.general.debug > 0 {
            DEBUG_LOGGER.set_state(EvalState::General);
        }

        let (res, stable, precision, duration) = results_of_stability_levels.last().unwrap_or_else(
            || panic!("No evaluation was done, perhaps the final stability level has a non-negative escalation threshold?")
        );

        if self.global_data.settings.general.debug > 0 {
            DEBUG_LOGGER.write("final_result", res);
        }

        let mut integrand_result = *res;

        let is_nan = integrand_result.re.is_nan() || integrand_result.im.is_nan() || !stable;

        if is_nan {
            integrand_result = Complex::new_zero();
        }

        let evaluation_metadata = EvaluationMetaData {
            rep3d_evaluation_time: *duration,
            parameterization_time,
            relative_instability_error: Complex::new_zero(),
            highest_precision: *precision,
            total_timing: start_evaluate_sample.elapsed(),
            is_nan,
        };

        EvaluationResult {
            integrand_result,
            integrator_weight: wgt,
            event_buffer: vec![],
            evaluation_metadata,
        }
    }

    fn get_integrator_settings(&self) -> crate::IntegratorSettings {
        self.global_data.settings.integrator.clone()
    }

    fn get_event_manager_mut(&mut self) -> &mut crate::observables::EventManager {
        todo!()
    }

    fn get_n_dim(&self) -> usize {
        match &self.graph_integrands {
            GraphIntegrands::Amplitude(graphs) => get_loop_count(&graphs[0]) * 3,
            GraphIntegrands::CrossSection(graphs) => get_loop_count(&graphs[0]) * 3,
        }
    }

    fn merge_results<I: HasIntegrand>(&mut self, _other: &mut I, _iter: usize) {}

    fn update_results(&mut self, _iter: usize) {}
}

impl GammaLoopIntegrand {
    #[inline]
    /// Evaluate the 3D representation of the integrand at a concrete floating-point precision.
    /// This function performs the evaluation twice, once for the original sample and once for the rotated sample.
    fn evaluate_at_prec(
        &mut self,
        samples: &[GammaLoopSample<f64>],
        precision: Precision,
    ) -> (Vec<Complex<F<f64>>>, Duration) {
        // measure timing if we are below the max number if we are below the max number
        let start = std::time::Instant::now();
        // cast the momenta to the relevant precision
        let results: Vec<_> = match precision {
            Precision::Single => {
                unimplemented!("From<f64> for f32 can't be implemented")
            }
            Precision::Double => match &mut self.graph_integrands {
                GraphIntegrands::Amplitude(graph_integrands) => samples
                    .iter()
                    .zip(self.global_data.rotations.iter())
                    .map(|(sample, rotation_for_overlap)| {
                        if self.global_data.settings.general.debug > 0 {
                            DEBUG_LOGGER.set_state(EvalState::PrecRot((
                                rotation_for_overlap.setting(),
                                precision,
                            )));
                        }

                        evaluate(
                            graph_integrands,
                            sample,
                            rotation_for_overlap,
                            &self.global_data.settings,
                        )
                    })
                    .collect(),
                GraphIntegrands::CrossSection(graph_integrands) => samples
                    .iter()
                    .zip(self.global_data.rotations.iter())
                    .map(|(sample, rotate_overlap_centers)| {
                        if self.global_data.settings.general.debug > 0 {
                            DEBUG_LOGGER.set_state(EvalState::PrecRot((
                                rotate_overlap_centers.setting(),
                                precision,
                            )));
                        }
                        evaluate(
                            graph_integrands,
                            sample,
                            rotate_overlap_centers,
                            &self.global_data.settings,
                        )
                    })
                    .collect(),
            },
            Precision::Quad => match &mut self.graph_integrands {
                GraphIntegrands::Amplitude(graph_integrands) => samples
                    .iter()
                    .zip(self.global_data.rotations.iter())
                    .map(|(sample, rotate_overlap_centers)| {
                        if self.global_data.settings.general.debug > 0 {
                            DEBUG_LOGGER.set_state(EvalState::PrecRot((
                                rotate_overlap_centers.setting(),
                                precision,
                            )));
                        }

                        evaluate(
                            graph_integrands,
                            &sample.higher_precision(),
                            rotate_overlap_centers,
                            &self.global_data.settings,
                        )
                        .lower()
                    })
                    .collect(),
                GraphIntegrands::CrossSection(graph_integrands) => samples
                    .iter()
                    .zip(self.global_data.rotations.iter())
                    .map(|(sample, rotate_overlap_centers)| {
                        if self.global_data.settings.general.debug > 0 {
                            DEBUG_LOGGER.set_state(EvalState::PrecRot((
                                rotate_overlap_centers.setting(),
                                precision,
                            )));
                        }

                        evaluate(
                            graph_integrands,
                            &sample.higher_precision(),
                            rotate_overlap_centers,
                            &self.global_data.settings,
                        )
                        .lower()
                    })
                    .collect(),
            },
            Precision::Arb(_prec) => {
                unimplemented!("need better traits to use arb prec")
            }
        };

        let duration = start.elapsed() / (samples.len() as u32);

        (results, duration)
    }

    /// Used to create the use_data_generator closure for havana_integrate
    pub fn user_data_generator(&self, num_cores: usize, _settings: &Settings) -> UserData {
        UserData {
            integrand: vec![Integrand::GammaLoopIntegrand(self.clone()); num_cores],
        }
    }

    /// Create an iterator which specifies the stability levels to be used for the evaluation
    #[inline]
    fn create_stability_vec(&self, use_f128: bool) -> Vec<StabilityLevelSetting> {
        if use_f128 {
            // overwrite the stability settings if use_f128 is enabled
            vec![StabilityLevelSetting {
                precision: Precision::Quad,
                required_precision_for_re: F(1e-5),
                required_precision_for_im: F(1e-5),
                escalate_for_large_weight_threshold: F(-1.),
            }]
        } else {
            self.global_data.settings.stability.levels.clone()
        }
    }

    /// Perform map from unit hypercube to 3-momenta
    #[inline]
    fn parameterize(&self, sample_point: &Sample<F<f64>>) -> Result<GammaLoopSample<f64>, String> {
        match &self.global_data.settings.sampling {
            SamplingSettings::Default => {
                let xs = unwrap_cont_sample(sample_point);
                Ok(GammaLoopSample::Default(self.default_parametrize(xs, 0)))
            }
            SamplingSettings::MultiChanneling(multichanneling_settings) => {
                let xs = unwrap_cont_sample(sample_point);
                Ok(GammaLoopSample::MultiChanneling {
                    alpha: multichanneling_settings.alpha,
                    sample: self.default_parametrize(xs, 0),
                })
            }
            SamplingSettings::DiscreteGraphs(discrete_graph_settings) => {
                match discrete_graph_settings {
                    DiscreteGraphSamplingSettings::Default => {
                        let (graph_id, xs) = unwrap_single_discrete_sample(sample_point);
                        Ok(GammaLoopSample::DiscreteGraph {
                            graph_id,
                            sample: DiscreteGraphSample::Default(
                                self.default_parametrize(xs, graph_id),
                            ),
                        })
                    }
                    DiscreteGraphSamplingSettings::MultiChanneling(multichanneling_settings) => {
                        let (graph_id, xs) = unwrap_single_discrete_sample(sample_point);
                        Ok(GammaLoopSample::DiscreteGraph {
                            graph_id,
                            sample: DiscreteGraphSample::MultiChanneling {
                                alpha: multichanneling_settings.alpha,
                                sample: self.default_parametrize(xs, graph_id),
                            },
                        })
                    }
                    DiscreteGraphSamplingSettings::TropicalSampling(tropical_sampling_settings) => {
                        let (graph_id, xs) = unwrap_single_discrete_sample(sample_point);
                        let externals = &self.global_data.settings.kinematics.externals;

                        let graph = match &self.graph_integrands {
                            GraphIntegrands::Amplitude(graphs) => &graphs[graph_id].graph,
                            GraphIntegrands::CrossSection(_graphs) => unimplemented!(), //,
                        };

                        let sampler = graph
                            .derived_data
                            .as_ref()
                            .unwrap()
                            .tropical_subgraph_table
                            .as_ref()
                            .expect("No tropical subgraph table present, disable tropical sampling or regenerate process with table");

                        let edge_data = graph
                            .bare_graph
                            .get_loop_edges_iterator()
                            .map(|(edge_id, edge)| {
                                let mass = edge.particle.mass.value;
                                let mass_re = mass.map(|complex_mass| complex_mass.re);

                                let shift = utils::compute_shift_part(
                                    &graph.bare_graph.loop_momentum_basis.edge_signatures[edge_id]
                                        .external,
                                    &externals.get_indep_externals(),
                                )
                                .spatial;

                                let shift_momtrop =
                                    Vector::from_array([shift.px, shift.py, shift.pz]);

                                (mass_re, shift_momtrop)
                            })
                            .collect_vec();

                        let sampling_result_result = sampler.generate_sample_from_x_space_point(
                            xs,
                            edge_data,
                            &tropical_sampling_settings.into_tropical_sampling_settings(
                                self.global_data.settings.general.debug,
                            ),
                            &DEBUG_LOGGER,
                        );

                        let sampling_result = match sampling_result_result {
                            Ok(sampling_result) => sampling_result,
                            Err(_) => {
                                return Err(String::from("tropical sampling failed"));
                            }
                        };

                        let loop_moms = sampling_result
                            .loop_momenta
                            .into_iter()
                            .map(Into::<ThreeMomentum<F<f64>>>::into)
                            .collect_vec();

                        let default_sample = DefaultSample::new(
                            loop_moms,
                            externals,
                            sampling_result.jacobian * externals.pdf(xs),
                            &self.global_data.polarizations[0],
                            &graph.bare_graph.external_in_or_out_signature(),
                        );
                        Ok(GammaLoopSample::DiscreteGraph {
                            graph_id,
                            sample: DiscreteGraphSample::Tropical(default_sample),
                        })
                    }
                    DiscreteGraphSamplingSettings::DiscreteMultiChanneling(
                        multichanneling_settings,
                    ) => {
                        let (graph_id, (channel_id, xs)) =
                            unwrap_double_discrete_sample(sample_point);
                        Ok(GammaLoopSample::DiscreteGraph {
                            graph_id,
                            sample: DiscreteGraphSample::DiscreteMultiChanneling {
                                alpha: multichanneling_settings.alpha,
                                channel_id,
                                sample: self.default_parametrize(xs, graph_id),
                            },
                        })
                    }
                }
            }
        }
    }

    /// Default parametrize is basically everything except tropical sampling.
    #[inline]
    fn default_parametrize(&self, xs: &[F<f64>], graph_id: usize) -> DefaultSample<f64> {
        let externals = &self.global_data.settings.kinematics.externals;

        let (loop_moms_vec, param_jacobian) = global_parameterize(
            xs,
            self.global_data.settings.kinematics.e_cm.square(),
            &self.global_data.settings,
            false,
        );

        let loop_moms = loop_moms_vec
            .into_iter()
            .map(ThreeMomentum::from)
            .collect_vec();

        let jacobian = param_jacobian * externals.pdf(xs);

        let graph = match &self.graph_integrands {
            GraphIntegrands::Amplitude(graphs) => &graphs[graph_id].graph,
            GraphIntegrands::CrossSection(_graphs) => unimplemented!(), //,
        };

        DefaultSample::new(
            loop_moms,
            externals,
            jacobian,
            &self.global_data.polarizations[0],
            &graph.bare_graph.external_in_or_out_signature(),
        )
    }

    /// Compute the average and check the accuracy of the result
    #[inline]
    fn stability_check(
        &self,
        results: &[Complex<F<f64>>],
        stability_settings: &StabilityLevelSetting,
        integrated_phase: IntegratedPhase,
        max_eval: F<f64>,
        wgt: F<f64>,
    ) -> (Complex<F<f64>>, bool) {
        if results.len() == 1 {
            return (results[0], true);
        }

        let average = results
            .iter()
            .fold(Complex::<F<f64>>::new_zero(), |acc, x| acc + x)
            / F(results.len() as f64);

        let mut errors = results.iter().map(|res| {
            let res_arr = [res.re, res.im];
            let avg_arr = [average.re, average.im];

            let (error_re, error_im) = res_arr
                .iter()
                .zip(avg_arr)
                .map(|(res_component, average_component)| {
                    if IsZero::is_zero(res_component) && IsZero::is_zero(&average_component) {
                        F(0.)
                    } else {
                        ((res_component - average_component) / average_component).abs()
                    }
                })
                .collect_tuple()
                .unwrap();
            Complex::new(error_re, error_im)
        });

        let unstable_sample = errors.position(|error| {
            error.re > stability_settings.required_precision_for_re
                || error.im > stability_settings.required_precision_for_im
        });

        if self.global_data.settings.general.debug > 1 {
            if let Some(unstable_index) = unstable_sample {
                let unstable_point = results[unstable_index];
                let rotation_axis = format!(
                    "{:?}",
                    self.global_data.settings.stability.rotation_axis[unstable_index]
                );

                let (
                    (real_formatted, rotated_real_formatted),
                    (imag_formatted, rotated_imag_formatted),
                ) = (
                    format_for_compare_digits(average.re, unstable_point.re),
                    format_for_compare_digits(average.im, unstable_point.im),
                );

                println!("{}", "\nUnstable point detected:".red());
                println!("Rotation axis: {}", rotation_axis);
                println!("\taverage result: {} + {}i", real_formatted, imag_formatted,);
                println!(
                    "\trotated result: {} + {}i",
                    rotated_real_formatted, rotated_imag_formatted,
                );
            }
        }

        let stable = unstable_sample.is_none();

        let average_for_comparison = match integrated_phase {
            IntegratedPhase::Real => average.re,
            IntegratedPhase::Imag => average.im,
            IntegratedPhase::Both => {
                unimplemented!("max wgt test unimplemented for integrated phase both")
            }
        };

        let below_wgt_threshold = if stability_settings.escalate_for_large_weight_threshold > F(0.)
            && max_eval.is_non_zero()
        {
            average_for_comparison.abs() * wgt
                < stability_settings.escalate_for_large_weight_threshold * max_eval.abs()
        } else {
            true
        };

        (average, stable && below_wgt_threshold)
    }

    #[inline]
    fn compute_2pi_factor(&self) -> f64 {
        let loop_number = match &self.graph_integrands {
            GraphIntegrands::Amplitude(graph_integrands) => get_loop_count(&graph_integrands[0]),
            GraphIntegrands::CrossSection(graph_integrands) => get_loop_count(&graph_integrands[0]),
        };
        (2. * std::f64::consts::PI).powi(loop_number as i32 * 3)
    }

    pub fn amplitude_integrand_constructor(
        mut amplitude: Amplitude<Evaluators>,
        settings: Settings,
    ) -> Self {
        // let pols = None;

        // for amplitudes we can construct counterterms beforehand if external momenta are constant
        let global_data = match settings.kinematics.externals {
            Externals::Constant { .. } => {
                let external_moms = settings.kinematics.externals.get_indep_externals();

                for amplitude_graph in amplitude.amplitude_graphs.iter_mut() {
                    let graph = &mut amplitude_graph.graph;

                    // temporary fix, rederive esurface data
                    graph
                        .generate_esurface_data()
                        .unwrap_or_else(|_| panic!("failed to generate esurface derived data"));

                    let existing_esurfaces = graph.get_existing_esurfaces(
                        &external_moms,
                        settings.kinematics.e_cm,
                        &settings,
                    );

                    if settings.general.debug > 0 {
                        println!(
                            "#{} existing esurfaces for graph {}",
                            existing_esurfaces.len(),
                            graph.bare_graph.name
                        );
                    }

                    if !existing_esurfaces.is_empty() {
                        // if settings.general.force_orientations.is_some() {
                        //     panic!("force orientations not supported with thresholds")
                        // }

                        match &settings.sampling {
                            SamplingSettings::Default => {}
                            SamplingSettings::MultiChanneling(_) => panic!(),
                            SamplingSettings::DiscreteGraphs(discrete_graph_settings) => {
                                match discrete_graph_settings {
                                    DiscreteGraphSamplingSettings::Default => {}
                                    DiscreteGraphSamplingSettings::DiscreteMultiChanneling(_) => {
                                        panic!()
                                    }
                                    DiscreteGraphSamplingSettings::MultiChanneling(_) => panic!(),
                                    DiscreteGraphSamplingSettings::TropicalSampling(_) => {}
                                }
                            }
                        }

                        let maximal_overlap = graph.get_maximal_overlap(
                            &external_moms,
                            settings.kinematics.e_cm,
                            &settings,
                        );

                        let maximal_overlap_structure = maximal_overlap
                            .overlap_groups
                            .iter()
                            .map(|overlap_group| overlap_group.existing_esurfaces.len())
                            .collect_vec();

                        if settings.general.debug > 0 {
                            println!("maximal overlap structure: {:?}", maximal_overlap_structure);
                        }

                        let counter_term =
                            CounterTerm::construct(maximal_overlap, &existing_esurfaces, graph);

                        if settings.general.debug > 1 {
                            counter_term.print_debug_data(
                                &graph.get_cff().esurfaces,
                                &external_moms,
                                &graph.bare_graph.loop_momentum_basis,
                                &graph.bare_graph.get_real_mass_vector(),
                            );
                        }

                        if let Some(derived_data) = &mut graph.derived_data {
                            derived_data.static_counterterm = Some(counter_term);
                        }
                    }
                }

                let rotations: Vec<Rotation> = Some(Rotation::new(RotationMethod::Identity))
                    .into_iter()
                    .chain(
                        settings
                            .stability
                            .rotation_axis
                            .iter()
                            .map(|axis| Rotation::new(axis.rotation_method())),
                    )
                    .collect(); // want this to include the identity rotation (i.e the first sample)

                let orig_polarizations = amplitude.polarizations(&settings.kinematics.externals);
                let polarizations = rotations
                    .iter()
                    .map(|r| orig_polarizations.rotate(r))
                    .collect();

                let externals = rotations
                    .iter()
                    .map(|r| settings.kinematics.externals.rotate(r))
                    .collect();
                GlobalData {
                    rotations,
                    polarizations,
                    externals,
                    settings,
                }
            }
        };

        Self {
            global_data,
            graph_integrands: GraphIntegrands::Amplitude(amplitude.amplitude_graphs),
        }
    }

    pub fn cross_section_integrand_constructor(
        cross_section: CrossSection,
        settings: Settings,
    ) -> Self {
        let rotations: Vec<Rotation> = settings
            .stability
            .rotation_axis
            .iter()
            .map(|axis| Rotation::new(axis.rotation_method()))
            .collect();

        let externals = rotations
            .iter()
            .map(|r| settings.kinematics.externals.rotate(r))
            .collect();
        let global_data = GlobalData {
            rotations,
            externals,
            polarizations: vec![Polarizations::None],
            settings,
        };

        Self {
            global_data,
            graph_integrands: GraphIntegrands::CrossSection(cross_section.supergraphs),
        }
    }
}

/// Sample whose structure depends on the sampling settings, and enforces these settings.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum GammaLoopSample<T: FloatLike> {
    Default(DefaultSample<T>),
    MultiChanneling {
        alpha: f64,
        sample: DefaultSample<T>,
    },
    DiscreteGraph {
        graph_id: usize,
        sample: DiscreteGraphSample<T>,
    },
}

impl GammaLoopSample<f64> {
    /// Rotation for stability checks
    #[inline]
    fn get_rotated_sample(
        &self,
        rotation: &Rotation,
        rotated_externals: Externals,
        rotated_polarizations: Polarizations,
    ) -> Self {
        if rotation.is_identity() {
            return self.clone();
        }

        let rotated_externals = rotated_externals.get_indep_externals();
        let rotated_polarizations = match rotated_polarizations {
            Polarizations::None => vec![],
            Polarizations::Constant { polarizations } => polarizations,
        };
        match self {
            GammaLoopSample::Default(sample) => {
                GammaLoopSample::Default(sample.get_rotated_sample_cached(
                    rotation,
                    rotated_externals,
                    rotated_polarizations,
                ))
            }
            GammaLoopSample::MultiChanneling { alpha, sample } => {
                GammaLoopSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.get_rotated_sample_cached(
                        rotation,
                        rotated_externals,
                        rotated_polarizations,
                    ),
                }
            }
            GammaLoopSample::DiscreteGraph { graph_id, sample } => GammaLoopSample::DiscreteGraph {
                graph_id: *graph_id,
                sample: sample.get_rotated_sample(
                    rotation,
                    rotated_externals,
                    rotated_polarizations,
                ),
            },
        }
    }
}

impl<T: FloatLike> GammaLoopSample<T> {
    pub fn zero(&self) -> F<T> {
        match self {
            GammaLoopSample::Default(sample) => sample.zero(),
            GammaLoopSample::MultiChanneling { sample, .. } => sample.zero(),
            GammaLoopSample::DiscreteGraph { sample, .. } => sample.zero(),
        }
    }

    #[allow(dead_code)]
    pub fn one(&self) -> F<T> {
        match self {
            GammaLoopSample::Default(sample) => sample.one(),
            GammaLoopSample::MultiChanneling { sample, .. } => sample.one(),
            GammaLoopSample::DiscreteGraph { sample, .. } => sample.one(),
        }
    }

    /// Cast the sample to a different precision
    #[allow(dead_code)]
    #[inline]
    fn cast_sample<T2: FloatLike>(&self) -> GammaLoopSample<T2>
    where
        F<T2>: From<F<T>>,
    {
        match self {
            GammaLoopSample::Default(sample) => GammaLoopSample::Default(sample.cast_sample()),
            GammaLoopSample::MultiChanneling { alpha, sample } => {
                GammaLoopSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.cast_sample(),
                }
            }
            GammaLoopSample::DiscreteGraph { graph_id, sample } => GammaLoopSample::DiscreteGraph {
                graph_id: *graph_id,
                sample: sample.cast_sample(),
            },
        }
    }

    fn higher_precision(&self) -> GammaLoopSample<T::Higher>
    where
        T::Higher: FloatLike,
    {
        match self {
            GammaLoopSample::Default(sample) => GammaLoopSample::Default(sample.higher_precision()),
            GammaLoopSample::MultiChanneling { alpha, sample } => {
                GammaLoopSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.higher_precision(),
                }
            }
            GammaLoopSample::DiscreteGraph { graph_id, sample } => GammaLoopSample::DiscreteGraph {
                graph_id: *graph_id,
                sample: sample.higher_precision(),
            },
        }
    }

    #[allow(dead_code)]
    fn lower_precision(&self) -> GammaLoopSample<T::Lower>
    where
        T::Lower: FloatLike,
    {
        match self {
            GammaLoopSample::Default(sample) => GammaLoopSample::Default(sample.lower_precision()),
            GammaLoopSample::MultiChanneling { alpha, sample } => {
                GammaLoopSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.lower_precision(),
                }
            }
            GammaLoopSample::DiscreteGraph { graph_id, sample } => GammaLoopSample::DiscreteGraph {
                graph_id: *graph_id,
                sample: sample.lower_precision(),
            },
        }
    }

    /// Retrieve the default sample which is contained in all types
    #[inline]
    fn get_default_sample(&self) -> &DefaultSample<T> {
        match self {
            GammaLoopSample::Default(sample) => sample,
            GammaLoopSample::MultiChanneling { sample, .. } => sample,
            GammaLoopSample::DiscreteGraph { sample, .. } => sample.get_default_sample(),
        }
    }
}

/// Sample which contains loop momenta, external momenta and the jacobian of the parameterization.
/// External momenta are part of the sample in order to facilitate the use of non-constant externals.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DefaultSample<T: FloatLike> {
    pub sample: BareSample<T>,
    pub rotated_sample: Option<BareSample<T>>,
    pub uuid: Uuid,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BareSample<T: FloatLike> {
    pub loop_moms: Vec<ThreeMomentum<F<T>>>,
    pub external_moms: Vec<FourMomentum<F<T>>>,
    pub polarizations: Vec<Polarization<Complex<F<T>>>>,
    pub jacobian: F<f64>,
}

use uuid::Uuid;

// static SAMPLECOUNT: LazyLock<Mutex<u64>> = LazyLock::new(|| Mutex::new(0));

impl<T: FloatLike> Display for DefaultSample<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Sample")?;
        write!(f, "\n\tloop momenta: ")?;
        for (index, loop_mom) in self.sample.loop_moms.iter().enumerate() {
            write!(f, "\n\t\tloop momentum {}: {}", index, loop_mom)?;
        }
        write!(f, "\n\texternal momenta: ")?;
        for (index, external_mom) in self.sample.external_moms.iter().enumerate() {
            write!(f, "\n\t\texternal momentum {}: {}", index, external_mom)?;
        }
        write!(f, "\n\tpolarizations: ")?;
        for (index, polarization) in self.sample.polarizations.iter().enumerate() {
            write!(f, "\n\t\tpolarization {}: {}", index, polarization)?;
        }
        write!(f, "\n\tjacobian: {:+e}", self.sample.jacobian)
    }
}

impl DefaultSample<f64> {}

impl<T: FloatLike> BareSample<T> {
    pub fn new(
        loop_moms: Vec<ThreeMomentum<F<T>>>,
        external_moms: &Externals,
        jacobian: F<f64>,
        polarizations: &Polarizations,
        external_signature: &Signature,
    ) -> Self {
        let polarizations = match polarizations {
            Polarizations::None => vec![],
            Polarizations::Constant { polarizations } => polarizations
                .iter()
                .map(|p| p.map(|c| c.map(|f| F::from_ff64(f))))
                .collect(),
        };

        let external_moms = external_moms.get_dependent_externals(external_signature);
        Self {
            polarizations,
            loop_moms,
            external_moms,
            jacobian,
        }
    }

    pub fn one(&self) -> F<T> {
        if let Some(f) = self.loop_moms.first() {
            f.px.one()
        } else if let Some(f) = self.external_moms.first() {
            return f.spatial.px.one();
        } else {
            panic!("No momenta in sample")
        }
    }

    pub fn zero(&self) -> F<T> {
        if let Some(f) = self.loop_moms.first() {
            f.px.zero()
        } else if let Some(f) = self.external_moms.first() {
            return f.spatial.px.zero();
        } else {
            panic!("No momenta in sample")
        }
    }

    /// Cast the sample to a different precision
    #[inline]
    fn cast_sample<T2: FloatLike>(&self) -> BareSample<T2>
    where
        F<T2>: From<F<T>>,
    {
        BareSample {
            loop_moms: self.loop_moms.iter().map(ThreeMomentum::cast).collect_vec(),
            external_moms: self
                .external_moms
                .iter()
                .map(FourMomentum::cast)
                .collect_vec(),
            polarizations: self
                .polarizations
                .iter()
                .map(Polarization::complex_cast)
                .collect_vec(),
            jacobian: self.jacobian,
        }
    }

    pub fn higher_precision(&self) -> BareSample<T::Higher>
    where
        T::Higher: FloatLike,
    {
        BareSample {
            loop_moms: self
                .loop_moms
                .iter()
                .map(ThreeMomentum::higher)
                .collect_vec(),
            external_moms: self
                .external_moms
                .iter()
                .map(FourMomentum::higher)
                .collect_vec(),
            polarizations: self
                .polarizations
                .iter()
                .map(Polarization::higher)
                .collect_vec(),
            jacobian: self.jacobian,
        }
    }

    pub fn lower_precision(&self) -> BareSample<T::Lower>
    where
        T::Lower: FloatLike,
    {
        BareSample {
            loop_moms: self
                .loop_moms
                .iter()
                .map(ThreeMomentum::lower)
                .collect_vec(),
            external_moms: self
                .external_moms
                .iter()
                .map(FourMomentum::lower)
                .collect_vec(),
            polarizations: self
                .polarizations
                .iter()
                .map(Polarization::lower)
                .collect_vec(),
            jacobian: self.jacobian,
        }
    }

    #[inline]
    /// Rotation for stability checks
    pub fn get_rotated_sample_cached(
        &self,
        rotation: &Rotation,
        rotated_externals: Vec<FourMomentum<F<T>>>,
        rotated_polarizations: Vec<Polarization<Complex<F<T>>>>,
    ) -> Self {
        Self {
            loop_moms: self
                .loop_moms
                .iter()
                .map(|l| l.rotate(rotation))
                .collect_vec(),
            external_moms: rotated_externals,
            polarizations: rotated_polarizations,
            jacobian: self.jacobian,
        }
    }

    #[inline]
    pub fn get_rotated_sample(&self, rotation: &Rotation) -> Self {
        Self {
            loop_moms: self
                .loop_moms
                .iter()
                .map(|l| l.rotate(rotation))
                .collect_vec(),
            external_moms: self
                .external_moms
                .iter()
                .map(|l| l.rotate(rotation))
                .collect_vec(),
            polarizations: self
                .polarizations
                .iter()
                .map(|l| l.rotate(rotation))
                .collect_vec(),
            jacobian: self.jacobian,
        }
    }
}

impl<T: FloatLike> DefaultSample<T> {
    pub fn possibly_rotated_sample(&self) -> &BareSample<T> {
        if let Some(rot) = self.rotated_sample.as_ref() {
            rot
        } else {
            &self.sample
        }
    }

    pub fn numerator_sample(&self, settings: &Settings) -> (&BareSample<T>, Option<Uuid>) {
        if settings.stability.rotate_numerator {
            (self.possibly_rotated_sample(), self.uuid())
        } else {
            (&self.sample, self.uuid())
        }
    }

    pub fn uuid(&self) -> Option<Uuid> {
        if self.rotated_sample.is_some() {
            None
        } else {
            Some(self.uuid)
        }
    }

    pub fn loop_moms(&self) -> &[ThreeMomentum<F<T>>] {
        if let Some(rotated_sample) = &self.rotated_sample {
            &rotated_sample.loop_moms
        } else {
            &self.sample.loop_moms
        }
    }

    #[allow(clippy::type_complexity)]
    pub fn loop_mom_pair(&self) -> (&[ThreeMomentum<F<T>>], Option<&[ThreeMomentum<F<T>>]>) {
        (
            self.sample.loop_moms.as_slice(),
            self.rotated_sample.as_ref().map(|s| s.loop_moms.as_slice()),
        )
    }

    pub fn external_moms(&self) -> &[FourMomentum<F<T>>] {
        if let Some(rotated_sample) = &self.rotated_sample {
            &rotated_sample.external_moms
        } else {
            &self.sample.external_moms
        }
    }

    #[allow(clippy::type_complexity)]
    pub fn external_mom_pair(&self) -> (&[FourMomentum<F<T>>], Option<&[FourMomentum<F<T>>]>) {
        (
            self.sample.external_moms.as_slice(),
            self.rotated_sample
                .as_ref()
                .map(|s| s.external_moms.as_slice()),
        )
    }

    pub fn polarizations(&self) -> &[Polarization<Complex<F<T>>>] {
        if let Some(rotated_sample) = &self.rotated_sample {
            &rotated_sample.polarizations
        } else {
            &self.sample.polarizations
        }
    }

    #[allow(clippy::type_complexity)]
    pub fn polarizations_pair(
        &self,
    ) -> (
        &[Polarization<Complex<F<T>>>],
        Option<&[Polarization<Complex<F<T>>>]>,
    ) {
        (
            self.sample.polarizations.as_slice(),
            self.rotated_sample
                .as_ref()
                .map(|s| s.polarizations.as_slice()),
        )
    }

    pub fn jacobian(&self) -> F<f64> {
        if let Some(rotated_sample) = &self.rotated_sample {
            rotated_sample.jacobian
        } else {
            self.sample.jacobian
        }
    }

    pub fn new(
        loop_moms: Vec<ThreeMomentum<F<T>>>,
        external_moms: &Externals,
        jacobian: F<f64>,
        polarizations: &Polarizations,
        external_signature: &Signature,
    ) -> Self {
        Self {
            sample: BareSample::new(
                loop_moms,
                external_moms,
                jacobian,
                polarizations,
                external_signature,
            ),
            rotated_sample: None,
            uuid: Uuid::new_v4(),
        }
    }

    pub fn one(&self) -> F<T> {
        self.sample.one()
    }

    pub fn zero(&self) -> F<T> {
        self.sample.zero()
    }

    /// Cast the sample to a different precision
    #[inline]
    fn cast_sample<T2: FloatLike>(&self) -> DefaultSample<T2>
    where
        F<T2>: From<F<T>>,
    {
        DefaultSample {
            sample: self.sample.cast_sample(),
            rotated_sample: self.rotated_sample.as_ref().map(|s| s.cast_sample()),
            uuid: self.uuid,
        }
    }

    pub fn higher_precision(&self) -> DefaultSample<T::Higher>
    where
        T::Higher: FloatLike,
    {
        DefaultSample {
            sample: self.sample.higher_precision(),
            rotated_sample: self.rotated_sample.as_ref().map(|s| s.higher_precision()),
            uuid: self.uuid,
        }
    }

    pub fn lower_precision(&self) -> DefaultSample<T::Lower>
    where
        T::Lower: FloatLike,
    {
        DefaultSample {
            sample: self.sample.lower_precision(),
            rotated_sample: self.rotated_sample.as_ref().map(|s| s.lower_precision()),
            uuid: self.uuid,
        }
    }

    #[inline]
    /// Rotation for stability checks
    pub fn get_rotated_sample_cached(
        &self,
        rotation: &Rotation,
        rotated_externals: Vec<FourMomentum<F<T>>>,
        rotated_polarizations: Vec<Polarization<Complex<F<T>>>>,
    ) -> Self {
        Self {
            sample: self.sample.clone(),
            rotated_sample: Some(self.sample.get_rotated_sample_cached(
                rotation,
                rotated_externals,
                rotated_polarizations,
            )),
            uuid: self.uuid,
        }
    }

    #[inline]
    pub fn get_rotated_sample(&self, rotation: &Rotation) -> Self {
        Self {
            sample: self.sample.clone(),
            rotated_sample: Some(self.sample.get_rotated_sample(rotation)),
            uuid: self.uuid,
        }
    }
}

/// This sample is used when importance sampling over graphs is used.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum DiscreteGraphSample<T: FloatLike> {
    Default(DefaultSample<T>),
    MultiChanneling {
        alpha: f64,
        sample: DefaultSample<T>,
    },
    /// This variant is equivalent to Default, but needs to be handled differently in the evaluation.
    Tropical(DefaultSample<T>),
    DiscreteMultiChanneling {
        alpha: f64,
        channel_id: usize,
        sample: DefaultSample<T>,
    },
}

impl<T: FloatLike> DiscreteGraphSample<T> {
    pub fn zero(&self) -> F<T> {
        match self {
            DiscreteGraphSample::Default(sample) => sample.zero(),
            DiscreteGraphSample::MultiChanneling { sample, .. } => sample.zero(),
            DiscreteGraphSample::Tropical(sample) => sample.zero(),
            DiscreteGraphSample::DiscreteMultiChanneling { sample, .. } => sample.zero(),
        }
    }

    pub fn one(&self) -> F<T> {
        match self {
            DiscreteGraphSample::Default(sample) => sample.one(),
            DiscreteGraphSample::MultiChanneling { sample, .. } => sample.one(),
            DiscreteGraphSample::Tropical(sample) => sample.one(),
            DiscreteGraphSample::DiscreteMultiChanneling { sample, .. } => sample.one(),
        }
    }
    /// Rotation for stability checks
    #[inline]
    fn get_rotated_sample(
        &self,
        rotation: &Rotation,
        rotated_externals: Vec<FourMomentum<F<T>>>,
        rotated_polarizations: Vec<Polarization<Complex<F<T>>>>,
    ) -> Self {
        match self {
            DiscreteGraphSample::Default(sample) => {
                DiscreteGraphSample::Default(sample.get_rotated_sample_cached(
                    rotation,
                    rotated_externals,
                    rotated_polarizations,
                ))
            }
            DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                DiscreteGraphSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.get_rotated_sample_cached(
                        rotation,
                        rotated_externals,
                        rotated_polarizations,
                    ),
                }
            }
            DiscreteGraphSample::Tropical(sample) => {
                DiscreteGraphSample::Tropical(sample.get_rotated_sample_cached(
                    rotation,
                    rotated_externals,
                    rotated_polarizations,
                ))
            }
            DiscreteGraphSample::DiscreteMultiChanneling {
                alpha,
                channel_id,
                sample,
            } => DiscreteGraphSample::DiscreteMultiChanneling {
                alpha: *alpha,
                channel_id: *channel_id,
                sample: sample.get_rotated_sample_cached(
                    rotation,
                    rotated_externals,
                    rotated_polarizations,
                ),
            },
        }
    }

    /// Cast the sample to a different precision
    #[inline]
    fn cast_sample<T2: FloatLike>(&self) -> DiscreteGraphSample<T2>
    where
        F<T2>: From<F<T>>,
    {
        match self {
            DiscreteGraphSample::Default(sample) => {
                DiscreteGraphSample::Default(sample.cast_sample())
            }
            DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                DiscreteGraphSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.cast_sample(),
                }
            }
            DiscreteGraphSample::Tropical(sample) => {
                DiscreteGraphSample::Tropical(sample.cast_sample())
            }
            DiscreteGraphSample::DiscreteMultiChanneling {
                alpha,
                channel_id,
                sample,
            } => DiscreteGraphSample::DiscreteMultiChanneling {
                alpha: *alpha,
                channel_id: *channel_id,
                sample: sample.cast_sample(),
            },
        }
    }

    fn higher_precision(&self) -> DiscreteGraphSample<T::Higher>
    where
        T::Higher: FloatLike,
    {
        match self {
            DiscreteGraphSample::Default(sample) => {
                DiscreteGraphSample::Default(sample.higher_precision())
            }
            DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                DiscreteGraphSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.higher_precision(),
                }
            }
            DiscreteGraphSample::Tropical(sample) => {
                DiscreteGraphSample::Tropical(sample.higher_precision())
            }
            DiscreteGraphSample::DiscreteMultiChanneling {
                alpha,
                channel_id,
                sample,
            } => DiscreteGraphSample::DiscreteMultiChanneling {
                alpha: *alpha,
                channel_id: *channel_id,
                sample: sample.higher_precision(),
            },
        }
    }

    fn lower_precision(&self) -> DiscreteGraphSample<T::Lower>
    where
        T::Lower: FloatLike,
    {
        match self {
            DiscreteGraphSample::Default(sample) => {
                DiscreteGraphSample::Default(sample.lower_precision())
            }
            DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                DiscreteGraphSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.lower_precision(),
                }
            }
            DiscreteGraphSample::Tropical(sample) => {
                DiscreteGraphSample::Tropical(sample.lower_precision())
            }
            DiscreteGraphSample::DiscreteMultiChanneling {
                alpha,
                channel_id,
                sample,
            } => DiscreteGraphSample::DiscreteMultiChanneling {
                alpha: *alpha,
                channel_id: *channel_id,
                sample: sample.lower_precision(),
            },
        }
    }

    /// Retrieve the default sample which is contained in all types
    #[inline]
    fn get_default_sample(&self) -> &DefaultSample<T> {
        match self {
            DiscreteGraphSample::Default(sample) => sample,
            DiscreteGraphSample::MultiChanneling { sample, .. } => sample,
            DiscreteGraphSample::Tropical(sample) => sample,
            DiscreteGraphSample::DiscreteMultiChanneling { sample, .. } => sample,
        }
    }
}

// helper functions can maybe moved to utils
#[inline]
fn unwrap_cont_sample(sample: &Sample<F<f64>>) -> &[F<f64>] {
    if let Sample::Continuous(_, xs) = sample {
        xs
    } else {
        panic!("Invalid sample structure")
    }
}

#[inline]
fn unwrap_single_discrete_sample(sample: &Sample<F<f64>>) -> (usize, &[F<f64>]) {
    if let Sample::Discrete(_, index, Some(cont_sample)) = sample {
        (*index, unwrap_cont_sample(cont_sample))
    } else {
        panic!("Invalid sample structure")
    }
}

#[inline]
fn unwrap_double_discrete_sample(sample: &Sample<F<f64>>) -> (usize, (usize, &[F<f64>])) {
    if let Sample::Discrete(_, index, Some(discrete_sample)) = sample {
        (*index, unwrap_single_discrete_sample(discrete_sample))
    } else {
        panic!("Invalid sample structure")
    }
}
