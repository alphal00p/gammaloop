//! This module contains the GammaloopIntegrand struct, which represents the integrand for physical
//! amplitudes and Local Unitarity crosssections.

use core::panic;
use std::time::{Duration, Instant};

use crate::cross_section::{Amplitude, AmplitudeGraph, CrossSection, SuperGraph};
use crate::evaluation_result::{EvaluationMetaData, EvaluationResult};
use crate::graph::{EdgeType, Graph, LoopMomentumBasisSpecification};
use crate::integrands::{HasIntegrand, Integrand};
use crate::integrate::UserData;
use crate::momentum::{FourMomentum, ThreeMomentum};
use crate::subtraction::static_counterterm::CounterTerm;
use crate::utils::{
    self, format_for_compare_digits, get_n_dim_for_n_loop_momenta, global_parameterize, FloatLike,
    PrecisionUpgradable, F,
};
use crate::{
    DiscreteGraphSamplingSettings, Externals, IntegratedPhase, SamplingSettings, Settings,
};
use crate::{Precision, StabilityLevelSetting};
use colored::Colorize;
use itertools::Itertools;
use momtrop::vector::Vector;
use spenso::complex::Complex;
use spenso::contraction::IsZero;
use symbolica::domains::float::{NumericalFloatLike, Real};
use symbolica::numerical_integration::{ContinuousGrid, DiscreteGrid, Grid, Sample};

/// Trait to capture the common behaviour of amplitudes and cross sections
/// Mainly used to expose the properties of the underlying graph in both amplitudes and cross sections
trait GraphIntegrand {
    /// Get the underlying graph
    fn get_graph(&self) -> &Graph;

    /// Get the channels used for multi channeling
    fn get_multi_channeling_channels(&self) -> &[usize];

    /// Most basic form of evaluating the 3D representation of the underlying loop integral
    fn evaluate<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        rotate_overlap_centers: Option<usize>,
        settings: &Settings,
    ) -> Complex<F<T>>;

    /// Evaluate in a single LMB-channel
    fn evaluate_channel<T: FloatLike>(
        &self,
        channel_id: usize,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<F<T>>;

    /// Evaluate a sum over LMB-channels
    fn evaluate_channel_sum<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<F<T>>;

    /// Evaluate to use when tropical sampling, raises the power of the onshell energies in front of the 3D representation according to the chosen weights.
    fn evaluate_tropical<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        rotate_overlap_centers: Option<usize>,
        settings: &Settings,
    ) -> Complex<F<T>>;
}

impl GraphIntegrand for AmplitudeGraph {
    fn get_graph(&self) -> &Graph {
        &self.graph
    }

    fn get_multi_channeling_channels(&self) -> &[usize] {
        &self.multi_channeling_channels
    }

    #[inline]
    #[allow(unused_variables)]
    fn evaluate_channel<T: FloatLike>(
        &self,
        channel_id: usize,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<F<T>> {
        let one = sample.one();
        let zero = sample.zero();
        // reference to the list of all lmbs
        let lmb_list = self
            .get_graph()
            .derived_data
            .loop_momentum_bases
            .as_ref()
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

        let channels_lmbs = channels.iter().map(|&i| &lmb_list[i]);

        let rep3d = if settings.general.use_ltd {
            self.get_graph().evaluate_ltd_expression_in_lmb(
                &sample.loop_moms,
                &sample.external_moms,
                &lmb_specification,
            )
        } else {
            self.get_graph()
                .evaluate_cff_expression_in_lmb(sample, &lmb_specification, settings)
        };

        let onshell_energies = self.get_graph().compute_onshell_energies_in_lmb(
            &sample.loop_moms,
            &sample.external_moms,
            &lmb_specification,
        );

        let virtual_energies = self
            .get_graph()
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

        multichanneling_numerator * rep3d / denominator
    }

    #[inline]
    fn evaluate_channel_sum<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<F<T>> {
        let zero = sample.zero();

        // a bit annoying that this is duplicated from evaluate_channel
        let lmb_list = self
            .get_graph()
            .derived_data
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
        &self,
        sample: &DefaultSample<T>,
        rotate_overlap_centers: Option<usize>,
        settings: &Settings,
    ) -> Complex<F<T>> {
        let zero_builder = &sample.loop_moms[0].px;

        let rep3d = if settings.general.use_ltd {
            self.get_graph()
                .evaluate_ltd_expression(&sample.loop_moms, &sample.external_moms)
        } else {
            self.get_graph().evaluate_cff_expression(sample, settings)
        };

        let energy_product = self
            .get_graph()
            .compute_energy_product(&sample.loop_moms, &sample.external_moms);

        let counter_terms = self.get_graph().derived_data.static_counterterm.as_ref();

        let counter_term_eval = match counter_terms {
            None => Complex::new(zero_builder.zero(), zero_builder.zero()),
            Some(counter_term) => counter_term.evaluate(
                &sample.loop_moms,
                &sample.external_moms,
                self.get_graph(),
                rotate_overlap_centers,
                settings,
            ),
        };

        rep3d / energy_product - counter_term_eval
    }

    #[inline]
    fn evaluate_tropical<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        rotate_overlap_centers: Option<usize>,
        settings: &Settings,
    ) -> Complex<F<T>> where {
        let one = sample.one();

        let rep3d = if settings.general.use_ltd {
            self.get_graph()
                .evaluate_ltd_expression(&sample.loop_moms, &sample.external_moms)
        } else {
            self.get_graph().evaluate_cff_expression(sample, settings)
        };

        let onshell_energies = self
            .get_graph()
            .compute_onshell_energies(&sample.loop_moms, &sample.external_moms);

        let tropical_subgraph_table = self.get_graph().get_tropical_subgraph_table();

        let virtual_loop_energies = self
            .get_graph()
            .get_loop_edges_iterator()
            .map(|(index, _)| onshell_energies[index].clone());

        let weight_iterator = tropical_subgraph_table.iter_edge_weights();

        let energy_product = virtual_loop_energies
            .zip(weight_iterator)
            .map(|(energy, weight)| energy.powf(&F::<T>::from_f64(2. * weight - 1.)))
            .fold(one.clone(), |acc, x| acc * x); // should we put Product and Sum in FloatLike?

        let tree_like_energies = self
            .get_graph()
            .get_tree_level_edges_iterator()
            .map(|(index, _)| onshell_energies[index].clone());

        let tree_product =
            tree_like_energies.fold(one.clone(), |acc, x| acc * F::<T>::from_f64(2.) * x);

        let counterterm = match &self.get_graph().derived_data.static_counterterm {
            Some(counterterm) => {
                counterterm.evaluate(
                    &sample.loop_moms,
                    &sample.external_moms,
                    self.get_graph(),
                    rotate_overlap_centers,
                    settings,
                ) * self
                    .graph
                    .compute_energy_product(&sample.loop_moms, &sample.external_moms)
            }
            None => Complex::new(one.zero(), one.zero()),
        };

        (rep3d - counterterm) * energy_product / tree_product
    }
}

impl GraphIntegrand for SuperGraph {
    fn get_graph(&self) -> &Graph {
        &self.graph
    }

    fn get_multi_channeling_channels(&self) -> &[usize] {
        todo!()
    }

    #[allow(unused)]
    fn evaluate_channel<T: FloatLike>(
        &self,
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
        &self,
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
        &self,
        sample: &DefaultSample<T>,
        rotate_overlap_centers: Option<usize>,
        settings: &Settings,
    ) -> Complex<F<T>> {
        // sum over channels
        todo!()
    }

    #[allow(unused)]
    #[inline]
    fn evaluate_tropical<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        rotate_overlap_centers: Option<usize>,
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
        .loop_momentum_bases
        .as_ref()
        .unwrap_or_else(|| panic!("Loop momentum bases not generated"))
        .len()
}

fn get_loop_count<T: GraphIntegrand>(graph_integrand: &T) -> usize {
    graph_integrand.get_graph().loop_momentum_basis.basis.len()
}

/// Evaluate the sample correctly according to the sample type
#[inline]
fn evaluate<I: GraphIntegrand, T: FloatLike>(
    graph_integrands: &[I],
    sample: &GammaLoopSample<T>,
    rotate_overlap_centers: Option<usize>,
    settings: &Settings,
) -> Complex<F<T>> {
    let zero = sample.zero();
    match sample {
        GammaLoopSample::Default(sample) => graph_integrands
            .iter()
            .map(|g| g.evaluate(sample, rotate_overlap_centers, settings))
            .reduce(|acc, e| acc + &e)
            .unwrap_or(zero.clone().into()),
        GammaLoopSample::MultiChanneling { alpha, sample } => graph_integrands
            .iter()
            .map(|g| g.evaluate_channel_sum(sample, *alpha, settings))
            .reduce(|acc, e| acc + &e)
            .unwrap_or(zero.clone().into()),
        GammaLoopSample::DiscreteGraph { graph_id, sample } => {
            let graph = &graph_integrands[*graph_id];
            match sample {
                DiscreteGraphSample::Default(sample) => {
                    graph.evaluate(sample, rotate_overlap_centers, settings)
                }
                DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                    graph.evaluate_channel_sum(sample, *alpha, settings)
                }
                DiscreteGraphSample::Tropical(sample) => {
                    graph.evaluate_tropical(sample, rotate_overlap_centers, settings)
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
    Amplitude(Vec<AmplitudeGraph>),
    CrossSection(Vec<SuperGraph>),
}

/// GammaloopIntegrand contains a list of graphs and the settings.
#[derive(Clone)]
pub struct GammaLoopIntegrand {
    pub settings: Settings,
    graph_integrands: GraphIntegrands,
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
        self.graph_integrands.create_grid(&self.settings)
    }

    #[allow(unused_variables)]
    fn evaluate_sample(
        &self,
        sample: &symbolica::numerical_integration::Sample<F<f64>>,
        wgt: F<f64>,
        iter: usize,
        use_f128: bool,
        max_eval: F<f64>,
    ) -> EvaluationResult {
        let start_evaluate_sample = Instant::now();

        // setup the evaluation of the integrand in the different stability levels
        let mut results_of_stability_levels =
            Vec::with_capacity(self.settings.stability.levels.len());

        // create an iterator containing the information for evaluation at each stability level
        let stability_iterator = self.create_stability_iterator(use_f128);

        let before_parameterization = std::time::Instant::now();
        let sample_point = self.parameterize(sample);
        let parameterization_time = before_parameterization.elapsed();

        // rotate the momenta for the stability tests.
        let rotated_sample_points = self
            .settings
            .stability
            .rotation_axis
            .iter()
            .enumerate()
            .map(|(func_index, f)| {
                (
                    sample_point.get_rotated_sample(f.rotation_function()),
                    Some(func_index),
                )
            });

        let samples = [(sample_point.clone(), None)]
            .into_iter()
            .chain(rotated_sample_points)
            .collect_vec();

        // 1 / (2 pi )^L
        let prefactor = F(self.compute_2pi_factor().inv());

        // iterate over the stability levels, break if the point is stable
        for stability_level in stability_iterator {
            // evaluate the integrand at the current stability level
            let (results, duration) = self.evaluate_at_prec(&samples, stability_level.precision);
            let results_scaled = results
                .iter()
                .zip(samples.iter())
                .map(|(result, sample)| result * sample.0.get_default_sample().jacobian * prefactor)
                .collect_vec();

            // check for the stability
            let (avg_result, stable) = self.stability_check(
                &results_scaled,
                stability_level,
                self.settings.integrator.integrated_phase,
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

        let (res, stable, precision, duration) = results_of_stability_levels.last().unwrap_or_else(
            || panic!("No evaluation was done, perhaps the final stability level has a non-negative escalation threshold?")
        );

        if self.settings.general.debug > 0 {
            println!("{}", "
            |  DEBUG -------------------------------------------------------------------------------------  |".green());
            println!();

            println!("{}", "sample from havana: ".blue());
            println!("\tsampled x: {:?}", sample);
            println!();

            println!("{}", "parametrisation result".blue());

            for (sample, rotation) in samples.iter() {
                let rotation_string = match rotation {
                    None => String::from("None"),
                    Some(index) => format!("{:?}", self.settings.stability.rotation_axis[*index]),
                };

                println!("\trotation: {}", rotation_string);
                println!("{}", "\tloop momenta: ".yellow());

                let loop_moms = &sample.get_default_sample().loop_moms;
                for (index, loop_mom) in loop_moms.iter().enumerate() {
                    println!("\t\tloop momentum {}: {:?}", index, loop_mom);
                }

                println!("{}", "\texternal momenta: ".yellow());

                let external_moms = &sample.get_default_sample().external_moms;
                for (index, external_mom) in external_moms.iter().enumerate() {
                    println!("\t\texternal momentum {}: {:?}", index, external_mom);
                }
            }

            let jacobian = sample_point.get_default_sample().jacobian;
            println!("\t{}: {:+e}", "jacobian".yellow(), jacobian);
            println!();

            println!("{}", "evaluation result: ".blue());
            println!("{}: {:+e}", "\tcff expression: ".yellow(), res);

            println!("\t{}: {:+e}", "result".yellow(), res * prefactor);
        }

        let mut integrand_result = *res;

        let is_nan = integrand_result.re.is_nan() || integrand_result.im.is_nan();

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
        self.settings.integrator.clone()
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
        &self,
        samples: &[(GammaLoopSample<f64>, Option<usize>)],
        precision: Precision,
    ) -> (Vec<Complex<F<f64>>>, Duration) {
        // measure timing if we are below the max number if we are below the max number
        let start = std::time::Instant::now();
        // cast the momenta to the relevant precision
        let results: Vec<_> = match precision {
            Precision::Single => {
                unimplemented!("From<f64> for f32 can't be implemented")
            }
            Precision::Double => match &self.graph_integrands {
                GraphIntegrands::Amplitude(graph_integrands) => samples
                    .iter()
                    .map(|(sample, rotate_overlap_centers)| {
                        evaluate(
                            graph_integrands,
                            sample,
                            *rotate_overlap_centers,
                            &self.settings,
                        )
                    })
                    .collect(),
                GraphIntegrands::CrossSection(graph_integrands) => samples
                    .iter()
                    .map(|(sample, rotate_overlap_centers)| {
                        evaluate(
                            graph_integrands,
                            sample,
                            *rotate_overlap_centers,
                            &self.settings,
                        )
                    })
                    .collect(),
            },
            Precision::Quad => match &self.graph_integrands {
                GraphIntegrands::Amplitude(graph_integrands) => samples
                    .iter()
                    .map(|(sample, rotate_overlap_centers)| {
                        evaluate(
                            graph_integrands,
                            &sample.higher_precision(),
                            *rotate_overlap_centers,
                            &self.settings,
                        )
                        .lower()
                    })
                    .collect(),
                GraphIntegrands::CrossSection(graph_integrands) => samples
                    .iter()
                    .map(|(sample, rotate_overlap_centers)| {
                        evaluate(
                            graph_integrands,
                            &sample.higher_precision(),
                            *rotate_overlap_centers,
                            &self.settings,
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
    fn create_stability_iterator(
        &self,
        use_f128: bool,
    ) -> impl Iterator<Item = &StabilityLevelSetting> {
        if use_f128 {
            // overwrite the stability settings if use_f128 is enabled
            [StabilityLevelSetting {
                precision: Precision::Quad,
                required_precision_for_re: F(1e-5),
                required_precision_for_im: F(1e-5),
                escalate_for_large_weight_threshold: F(-1.),
            }]
            .iter()
        } else {
            self.settings.stability.levels.iter()
        }
    }

    /// Perform map from unit hypercube to 3-momenta
    #[inline]
    fn parameterize(&self, sample_point: &Sample<F<f64>>) -> GammaLoopSample<f64> {
        match &self.settings.sampling {
            SamplingSettings::Default => {
                let xs = unwrap_cont_sample(sample_point);
                GammaLoopSample::Default(self.default_parametrize(xs))
            }
            SamplingSettings::MultiChanneling(multichanneling_settings) => {
                let xs = unwrap_cont_sample(sample_point);
                GammaLoopSample::MultiChanneling {
                    alpha: multichanneling_settings.alpha,
                    sample: self.default_parametrize(xs),
                }
            }
            SamplingSettings::DiscreteGraphs(discrete_graph_settings) => {
                match discrete_graph_settings {
                    DiscreteGraphSamplingSettings::Default => {
                        let (graph_id, xs) = unwrap_single_discrete_sample(sample_point);
                        GammaLoopSample::DiscreteGraph {
                            graph_id,
                            sample: DiscreteGraphSample::Default(self.default_parametrize(xs)),
                        }
                    }
                    DiscreteGraphSamplingSettings::MultiChanneling(multichanneling_settings) => {
                        let (graph_id, xs) = unwrap_single_discrete_sample(sample_point);
                        GammaLoopSample::DiscreteGraph {
                            graph_id,
                            sample: DiscreteGraphSample::MultiChanneling {
                                alpha: multichanneling_settings.alpha,
                                sample: self.default_parametrize(xs),
                            },
                        }
                    }
                    DiscreteGraphSamplingSettings::TropicalSampling => {
                        let (graph_id, xs) = unwrap_single_discrete_sample(sample_point);
                        let (external_moms, pdf) =
                            self.settings.kinematics.externals.get_externals(xs);

                        let graph = match &self.graph_integrands {
                            GraphIntegrands::Amplitude(graphs) => &graphs[graph_id].graph,
                            GraphIntegrands::CrossSection(graphs) => &graphs[graph_id].graph,
                        };

                        let sampler = graph.derived_data.tropical_subgraph_table.as_ref().unwrap();
                        let xs_f64 = xs.iter().map(|x| x.0).collect_vec();

                        let print_debug_info = self.settings.general.debug > 3;

                        let edge_data = graph
                            .get_loop_edges_iterator()
                            .map(|(edge_id, edge)| {
                                let mass = edge.particle.mass.value;
                                let mass_re = mass.map(|complex_mass| complex_mass.re.0);

                                let shift = utils::compute_shift_part(
                                    &graph.loop_momentum_basis.edge_signatures[edge_id].1,
                                    &external_moms,
                                );

                                let shift_momtrop = Vector::from_array([
                                    shift.spatial.px.0,
                                    shift.spatial.py.0,
                                    shift.spatial.pz.0,
                                ]);

                                (mass_re, shift_momtrop)
                            })
                            .collect_vec();

                        let sampling_result = sampler.generate_sample_from_x_space_point(
                            &xs_f64,
                            edge_data,
                            print_debug_info,
                        );

                        let loop_moms = sampling_result
                            .loop_momenta
                            .into_iter()
                            .map(Into::<ThreeMomentum<F<f64>>>::into)
                            .collect_vec();

                        let default_sample = DefaultSample {
                            loop_moms,
                            external_moms,
                            jacobian: F(sampling_result.jacobian) * pdf,
                        };

                        GammaLoopSample::DiscreteGraph {
                            graph_id,
                            sample: DiscreteGraphSample::Tropical(default_sample),
                        }
                    }
                    DiscreteGraphSamplingSettings::DiscreteMultiChanneling(
                        multichanneling_settings,
                    ) => {
                        let (graph_id, (channel_id, xs)) =
                            unwrap_double_discrete_sample(sample_point);
                        GammaLoopSample::DiscreteGraph {
                            graph_id,
                            sample: DiscreteGraphSample::DiscreteMultiChanneling {
                                alpha: multichanneling_settings.alpha,
                                channel_id,
                                sample: self.default_parametrize(xs),
                            },
                        }
                    }
                }
            }
        }
    }

    /// Default parametrize is basically everything except tropical sampling.
    #[inline]
    fn default_parametrize(&self, xs: &[F<f64>]) -> DefaultSample<f64> {
        let (external_moms, pdf) = self.settings.kinematics.externals.get_externals(xs);
        let (loop_moms_vec, param_jacobian) = global_parameterize(
            xs,
            self.settings.kinematics.e_cm * self.settings.kinematics.e_cm,
            &self.settings,
            false,
        );

        let loop_moms = loop_moms_vec
            .into_iter()
            .map(ThreeMomentum::from)
            .collect_vec();

        let jacobian = param_jacobian * pdf;

        DefaultSample {
            loop_moms,
            external_moms,
            jacobian,
        }
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

        let (results_for_comparison, average_for_comparison, max_wgt_for_comparison) =
            match integrated_phase {
                IntegratedPhase::Real => (
                    results.iter().map(|r| r.re).collect_vec(),
                    average.re,
                    max_eval,
                ),
                IntegratedPhase::Imag => (
                    results.iter().map(|r| r.im).collect_vec(),
                    average.im,
                    max_eval,
                ),
                IntegratedPhase::Both => unimplemented!("integrated phase both not implemented"),
            };

        let mut errors = results_for_comparison.iter().map(|res| {
            if IsZero::is_zero(res) && IsZero::is_zero(&average_for_comparison) {
                F(0.)
            } else {
                ((res - average_for_comparison) / average_for_comparison).abs()
            }
        });

        let unstable_sample = match integrated_phase {
            IntegratedPhase::Real => {
                errors.position(|error| error > stability_settings.required_precision_for_re)
            }
            IntegratedPhase::Imag => {
                errors.position(|error| error > stability_settings.required_precision_for_im)
            }
            IntegratedPhase::Both => unimplemented!("integrated phase both not implemented"),
        };

        if self.settings.general.debug > 0 {
            if let Some(unstable_index) = unstable_sample {
                let unstable_point = results[unstable_index];
                let rotation_axis = format!(
                    "{:?}",
                    self.settings.stability.rotation_axis[unstable_index]
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

        let below_wgt_threshold = if stability_settings.escalate_for_large_weight_threshold > F(0.)
            && max_wgt_for_comparison.is_non_zero()
        {
            average_for_comparison.abs() * wgt
                < stability_settings.escalate_for_large_weight_threshold
                    * max_wgt_for_comparison.abs()
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

    pub fn amplitude_integrand_constructor(mut amplitude: Amplitude, settings: Settings) -> Self {
        #[allow(irrefutable_let_patterns)]
        // for amplitudes we can construct counterterms beforehand if external momenta are constant
        if let Externals::Constant(_) = settings.kinematics.externals {
            let dummy = [];
            let (external_moms, _) = settings.kinematics.externals.get_externals(&dummy);

            for amplitude_graph in amplitude.amplitude_graphs.iter_mut() {
                let graph = &mut amplitude_graph.graph;

                // temporary fix, rederive esurface data
                graph
                    .generate_esurface_data()
                    .unwrap_or_else(|_| panic!("failed to generate esurface derived data"));

                let existing_esurfaces = graph.get_existing_esurfaces(
                    &external_moms,
                    settings.kinematics.e_cm,
                    settings.general.debug,
                );

                if settings.general.debug > 0 {
                    println!(
                        "#{} existing esurfaces for graph {}",
                        existing_esurfaces.len(),
                        graph.name
                    );
                }

                if !existing_esurfaces.is_empty() {
                    let maximal_overlap = graph.get_maximal_overlap(
                        &external_moms,
                        settings.kinematics.e_cm,
                        settings.general.debug,
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
                            &graph.loop_momentum_basis,
                            &graph.get_real_mass_vector(),
                        );
                    }

                    graph.derived_data.static_counterterm = Some(counter_term);
                }
            }
        } else {
            panic!("Only constant external momenta are supported at the moment.")
        }

        Self {
            settings,
            graph_integrands: GraphIntegrands::Amplitude(amplitude.amplitude_graphs),
        }
    }

    pub fn cross_section_integrand_constructor(
        cross_section: CrossSection,
        settings: Settings,
    ) -> Self {
        Self {
            settings,
            graph_integrands: GraphIntegrands::CrossSection(cross_section.supergraphs),
        }
    }
}

/// Sample whose structure depends on the sampling settings, and enforces these settings.
#[derive(Debug, Clone)]
enum GammaLoopSample<T: FloatLike> {
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
    /// Rotation for stability checks
    #[inline]
    fn get_rotated_sample(
        &self,
        rotation_function: impl Fn(&ThreeMomentum<F<T>>) -> ThreeMomentum<F<T>>,
    ) -> Self {
        match self {
            GammaLoopSample::Default(sample) => {
                GammaLoopSample::Default(sample.get_rotated_sample(rotation_function))
            }
            GammaLoopSample::MultiChanneling { alpha, sample } => {
                GammaLoopSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.get_rotated_sample(rotation_function),
                }
            }
            GammaLoopSample::DiscreteGraph { graph_id, sample } => GammaLoopSample::DiscreteGraph {
                graph_id: *graph_id,
                sample: sample.get_rotated_sample(rotation_function),
            },
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
#[derive(Debug, Clone)]
pub struct DefaultSample<T: FloatLike> {
    pub loop_moms: Vec<ThreeMomentum<F<T>>>,
    pub external_moms: Vec<FourMomentum<F<T>>>,
    pub jacobian: F<f64>,
}

impl<T: FloatLike> DefaultSample<T> {
    pub fn one(&self) -> F<T> {
        self.loop_moms[0].px.one()
    }

    pub fn zero(&self) -> F<T> {
        self.loop_moms[0].px.zero()
    }

    #[inline]
    /// Rotation for stability checks
    fn get_rotated_sample(
        &self,
        rotation_function: impl Fn(&ThreeMomentum<F<T>>) -> ThreeMomentum<F<T>>,
    ) -> Self {
        Self {
            loop_moms: self.loop_moms.iter().map(&rotation_function).collect_vec(),
            external_moms: self
                .external_moms
                .iter()
                .cloned()
                .map(|mut p| {
                    p.spatial = rotation_function(&p.spatial);
                    p
                })
                .collect_vec(),
            jacobian: self.jacobian,
        }
    }

    /// Cast the sample to a different precision
    #[inline]
    fn cast_sample<T2: FloatLike>(&self) -> DefaultSample<T2>
    where
        F<T2>: From<F<T>>,
    {
        DefaultSample {
            loop_moms: self.loop_moms.iter().map(ThreeMomentum::cast).collect_vec(),
            external_moms: self
                .external_moms
                .iter()
                .map(FourMomentum::cast)
                .collect_vec(),
            jacobian: self.jacobian,
        }
    }

    fn higher_precision(&self) -> DefaultSample<T::Higher>
    where
        T::Higher: FloatLike,
    {
        DefaultSample {
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
            jacobian: self.jacobian,
        }
    }

    fn lower_precision(&self) -> DefaultSample<T::Lower>
    where
        T::Lower: FloatLike,
    {
        DefaultSample {
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
            jacobian: self.jacobian,
        }
    }
}

/// This sample is used when importance sampling over graphs is used.
#[derive(Debug, Clone)]
enum DiscreteGraphSample<T: FloatLike> {
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
        rotation_function: impl Fn(&ThreeMomentum<F<T>>) -> ThreeMomentum<F<T>>,
    ) -> Self {
        match self {
            DiscreteGraphSample::Default(sample) => {
                DiscreteGraphSample::Default(sample.get_rotated_sample(rotation_function))
            }
            DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                DiscreteGraphSample::MultiChanneling {
                    alpha: *alpha,
                    sample: sample.get_rotated_sample(rotation_function),
                }
            }
            DiscreteGraphSample::Tropical(sample) => {
                DiscreteGraphSample::Tropical(sample.get_rotated_sample(rotation_function))
            }
            DiscreteGraphSample::DiscreteMultiChanneling {
                alpha,
                channel_id,
                sample,
            } => DiscreteGraphSample::DiscreteMultiChanneling {
                alpha: *alpha,
                channel_id: *channel_id,
                sample: sample.get_rotated_sample(rotation_function),
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
