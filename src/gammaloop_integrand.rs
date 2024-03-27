//! This module contains the GammaloopIntegrand struct, which represents the integrand for physical
//! amplitudes and Local Unitarity crosssections.

use core::panic;
use std::time::{Duration, Instant};

use crate::cross_section::{Amplitude, AmplitudeGraph, CrossSection, SuperGraph};
use crate::evaluation_result::{EvaluationMetaData, EvaluationResult};
use crate::graph::{EdgeType, Graph, LoopMomentumBasisSpecification};
use crate::integrands::{HasIntegrand, Integrand};
use crate::integrate::UserData;
use crate::tropical::tropical_parameterization::{self};
use crate::utils::{
    cast_complex, cast_lorentz_vector, format_for_compare_digits, get_n_dim_for_n_loop_momenta,
    global_parameterize, FloatLike,
};
use crate::{DiscreteGraphSamplingSettings, IntegratedPhase, SamplingSettings, Settings};
use crate::{Precision, StabilityLevelSetting};
use colored::Colorize;
use itertools::Itertools;
use lorentz_vector::LorentzVector;
use num::traits::{Inv, Zero};
use num::Complex;
use symbolica::numerical_integration::{ContinuousGrid, DiscreteGrid, Grid, Sample};

/// Trait to capture the common behaviour of amplitudes and cross sections
/// Mainly used to expose the properties of the underlying graph in both amplitudes and cross sections
trait GraphIntegrand {
    /// Get the underlying graph
    fn get_graph(&self) -> &Graph;

    /// Get the channels used for multi channeling
    fn get_multi_channeling_channels(&self) -> &[usize];

    /// Most basic form of evaluating the 3D representation of the underlying loop integral
    fn evaluate<T: FloatLike>(&self, sample: &DefaultSample<T>, settings: &Settings) -> Complex<T>;

    /// Evaluate in a single LMB-channel
    fn evaluate_channel<T: FloatLike>(
        &self,
        channel_id: usize,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<T>;

    /// Evaluate a sum over LMB-channels
    fn evaluate_channel_sum<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<T>;

    /// Evaluate to use when tropical sampling, raises the power of the onshell energies in front of the 3D representation according to the chosen weights.
    fn evaluate_tropical<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        settings: &Settings,
    ) -> Complex<T>;
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
    ) -> Complex<T> {
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
            self.get_graph().evaluate_cff_expression_in_lmb(
                &sample.loop_moms,
                &sample.external_moms,
                &lmb_specification,
            )
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
            .map(|i| Into::<T>::into(2.) * onshell_energies[i])
            .fold(T::one(), |acc, x| acc * x);

        let lmb_products = channels_lmbs.map(|basis| {
            basis
                .basis
                .iter()
                .fold(T::one(), |acc, i| acc * onshell_energies[*i])
        });

        let denominator = lmb_products
            .map(|x| x.powf(-Into::<T>::into(alpha)) * energy_product)
            .sum::<T>();

        let multichanneling_numerator = Complex::new(
            lmb_list[channel]
                .basis
                .iter()
                .fold(T::one(), |acc, i| acc * onshell_energies[*i])
                .powf(-Into::<T>::into(alpha)),
            T::zero(),
        );

        multichanneling_numerator * rep3d / denominator
    }

    #[inline]
    fn evaluate_channel_sum<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<T> {
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
            .sum()
    }

    #[inline]
    fn evaluate<T: FloatLike>(&self, sample: &DefaultSample<T>, settings: &Settings) -> Complex<T> {
        let rep3d = if settings.general.use_ltd {
            self.get_graph()
                .evaluate_ltd_expression(&sample.loop_moms, &sample.external_moms)
        } else {
            self.get_graph()
                .evaluate_cff_expression(&sample.loop_moms, &sample.external_moms)
        };

        let energy_product = self
            .get_graph()
            .compute_energy_product(&sample.loop_moms, &sample.external_moms);

        rep3d / energy_product
    }

    #[inline]
    fn evaluate_tropical<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        settings: &Settings,
    ) -> Complex<T> {
        let rep3d = if settings.general.use_ltd {
            self.get_graph()
                .evaluate_ltd_expression(&sample.loop_moms, &sample.external_moms)
        } else {
            self.get_graph()
                .evaluate_cff_expression(&sample.loop_moms, &sample.external_moms)
        };

        let onshell_energies = self
            .get_graph()
            .compute_onshell_energies(&sample.loop_moms, &sample.external_moms);

        let tropical_subgraph_table = self
            .get_graph()
            .derived_data
            .tropical_subgraph_table
            .as_ref()
            .unwrap();

        let virtual_loop_energies = self
            .get_graph()
            .get_loop_edges_iterator()
            .map(|(index, _)| onshell_energies[index]);

        let weight_iterator = tropical_subgraph_table
            .tropical_graph
            .topology
            .iter()
            .map(|edge| edge.weight);

        let energy_product = virtual_loop_energies
            .zip(weight_iterator)
            .map(|(energy, weight)| energy.powf(Into::<T>::into(2. * weight - 1.)))
            .fold(T::one(), |acc, x| acc * x); // should we put Product and Sum in FloatLike?

        let tree_like_energies = self
            .get_graph()
            .get_tree_level_edges_iterator()
            .map(|(index, _)| onshell_energies[index]);

        let tree_product =
            tree_like_energies.fold(T::one(), |acc, x| acc * Into::<T>::into(2.) * x);

        rep3d * energy_product / tree_product
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
    ) -> Complex<T> {
        // sum over cuts
        todo!()
    }

    #[allow(unused)]
    fn evaluate_channel_sum<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        alpha: f64,
        settings: &Settings,
    ) -> Complex<T> {
        // sum over channels
        todo!()
    }

    #[allow(unused)]
    #[inline]
    fn evaluate<T: FloatLike>(&self, sample: &DefaultSample<T>, settings: &Settings) -> Complex<T> {
        // sum over channels
        todo!()
    }

    #[allow(unused)]
    #[inline]
    fn evaluate_tropical<T: FloatLike>(
        &self,
        sample: &DefaultSample<T>,
        settings: &Settings,
    ) -> Complex<T> {
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
fn evaluate<T: GraphIntegrand, F: FloatLike>(
    graph_integrands: &[T],
    sample: &GammaLoopSample<F>,
    settings: &Settings,
) -> Complex<F> {
    match sample {
        GammaLoopSample::Default(sample) => graph_integrands
            .iter()
            .map(|g| g.evaluate(sample, settings))
            .sum::<Complex<F>>(),
        GammaLoopSample::MultiChanneling { alpha, sample } => graph_integrands
            .iter()
            .map(|g| g.evaluate_channel_sum(sample, *alpha, settings))
            .sum::<Complex<F>>(),
        GammaLoopSample::DiscreteGraph { graph_id, sample } => {
            let graph = &graph_integrands[*graph_id];
            match sample {
                DiscreteGraphSample::Default(sample) => graph.evaluate(sample, settings),
                DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                    graph.evaluate_channel_sum(sample, *alpha, settings)
                }
                DiscreteGraphSample::Tropical(sample) => graph.evaluate_tropical(sample, settings),
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
fn create_grid<T: GraphIntegrand>(graph_integrand: &T, settings: &Settings) -> Grid<f64> {
    let num_loops = get_loop_count(graph_integrand);

    let n_edges = graph_integrand
        .get_graph()
        .derived_data
        .tropical_subgraph_table
        .as_ref()
        .map(|t| t.tropical_graph.topology.len());

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
    fn create_grid(&self, settings: &Settings) -> Grid<f64> {
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
    fn create_grid(&self) -> Grid<f64> {
        self.graph_integrands.create_grid(&self.settings)
    }

    #[allow(unused_variables)]
    fn evaluate_sample(
        &self,
        sample: &symbolica::numerical_integration::Sample<f64>,
        wgt: f64,
        iter: usize,
        use_f128: bool,
        max_eval: f64,
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
        let rotation_method = self.settings.stability.rotation_axis.rotation_function();
        let rotated_sample_point = sample_point.get_rotated_sample(rotation_method);

        // 1 / (2 pi )^L
        let prefactor = self.compute_2pi_factor().inv();

        // iterate over the stability levels, break if the point is stable
        for stability_level in stability_iterator {
            // evaluate the integrand at the current stability level
            let (result, rotated_result, duration) = self.evaluate_at_prec(
                &sample_point,
                &rotated_sample_point,
                stability_level.precision,
            );

            let (result_scaled, rotated_result_scaled) = (
                result * sample_point.get_default_sample().jacobian * prefactor,
                rotated_result * rotated_sample_point.get_default_sample().jacobian * prefactor,
            );

            // check for the stability
            let (avg_result, stable) = self.stability_check(
                result_scaled,
                rotated_result_scaled,
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

            if !stable {
                if self.settings.general.debug > 0 {
                    let (
                        (real_formatted, rotated_real_formatted),
                        (imag_formatted, rotated_imag_formatted),
                    ) = (
                        format_for_compare_digits(result.re, rotated_result.re),
                        format_for_compare_digits(result.im, rotated_result.im),
                    );

                    println!("{}", "\nEscalating to next stability level:".red());
                    println!("\tresult:         {} + {}i", real_formatted, imag_formatted,);
                    println!(
                        "\trotated result: {} + {}i",
                        rotated_real_formatted, rotated_imag_formatted,
                    );
                }
            } else {
                break;
            }
        }

        let (res, stable, precision, duration) = results_of_stability_levels.last().unwrap_or_else(
            || panic!("No evaluation was done, perhaps the final stability level has a non-negative escalation threshold?")
        );

        if self.settings.general.debug > 1 {
            println!("{}", "
            |  DEBUG -------------------------------------------------------------------------------------  |".green());
            println!();

            println!("{}", "sample form havana: ".blue());
            println!("\tsampled x: {:?}", sample_point);
            println!();

            println!("{}", "parametrisation result".blue());
            println!("{}", "\tloop momenta: ".yellow());

            let loop_moms = &sample_point.get_default_sample().loop_moms;
            for (index, loop_mom) in loop_moms.iter().enumerate() {
                println!("\t\tloop momentum {}: {:?}", index, loop_mom);
            }

            println!("{}", "\trotated loop momenta: ".yellow());

            let rotated_loop_moms = &rotated_sample_point.get_default_sample().loop_moms;
            for (index, loop_mom) in rotated_loop_moms.iter().enumerate() {
                println!("\t\tloop momentum {}: {:?}", index, loop_mom);
            }

            println!("{}", "\texternal momenta: ".yellow());

            let external_moms = &sample_point.get_default_sample().external_moms;
            for (index, external_mom) in external_moms.iter().enumerate() {
                println!("\t\texternal momentum {}: {:?}", index, external_mom);
            }

            println!("{}", "\trotated external momenta: ".yellow());

            let rotated_external_moms = &rotated_sample_point.get_default_sample().external_moms;
            for (index, external_mom) in rotated_external_moms.iter().enumerate() {
                println!("\t\texternal momentum {}: {:?}", index, external_mom);
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
            integrand_result = Complex::new(0., 0.);
        }

        let evaluation_metadata = EvaluationMetaData {
            rep3d_evaluation_time: *duration,
            parameterization_time,
            relative_instability_error: Complex::new(0., 0.),
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
        sample_point: &GammaLoopSample<f64>,
        rotated_sample_point: &GammaLoopSample<f64>,
        precision: Precision,
    ) -> (Complex<f64>, Complex<f64>, Duration) {
        // measure timing if we are below the max number if we are below the max number
        let start = std::time::Instant::now();
        // cast the momenta to the relevant precision
        let (result, rotated_result) = match precision {
            Precision::Single => {
                unimplemented!("From<f64> for f32 can't be implemented")
            }
            Precision::Double => match &self.graph_integrands {
                GraphIntegrands::Amplitude(graph_integrands) => {
                    let result = evaluate(graph_integrands, sample_point, &self.settings);
                    let rotated_result =
                        evaluate(graph_integrands, rotated_sample_point, &self.settings);

                    (result, rotated_result)
                }
                GraphIntegrands::CrossSection(graph_integrands) => {
                    let result = evaluate(graph_integrands, sample_point, &self.settings);
                    let rotated_result =
                        evaluate(graph_integrands, rotated_sample_point, &self.settings);

                    (result, rotated_result)
                }
            },
            Precision::Quad => {
                let sample_point_f128 = sample_point.cast_sample::<f128::f128>();
                let rotated_sample_point_f128 = rotated_sample_point.cast_sample::<f128::f128>();

                match &self.graph_integrands {
                    GraphIntegrands::Amplitude(graph_integrands) => {
                        let result = evaluate(graph_integrands, &sample_point_f128, &self.settings);

                        let rotated_result =
                            evaluate(graph_integrands, &rotated_sample_point_f128, &self.settings);

                        (cast_complex(result), cast_complex(rotated_result))
                    }
                    GraphIntegrands::CrossSection(graph_integrands) => {
                        let result = evaluate(graph_integrands, &sample_point_f128, &self.settings);

                        let rotated_result =
                            evaluate(graph_integrands, &sample_point_f128, &self.settings);

                        (cast_complex(result), cast_complex(rotated_result))
                    }
                }
            }
            Precision::Arb(_prec) => {
                unimplemented!("need better traits to use arb prec")
            }
        };

        let duration = start.elapsed() / 2;

        (result, rotated_result, duration)
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
                required_precision_for_re: 1e-5,
                required_precision_for_im: 1e-5,
                escalate_for_large_weight_threshold: -1.,
            }]
            .iter()
        } else {
            self.settings.stability.levels.iter()
        }
    }

    /// Perform map from unit hypercube to 3-momenta
    #[inline]
    fn parameterize(&self, sample_point: &Sample<f64>) -> GammaLoopSample<f64> {
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

                        let (loop_moms_f64, jacobian_f64) =
                            tropical_parameterization::generate_tropical_sample(
                                xs,
                                &external_moms,
                                graph,
                                self.settings.general.debug,
                            )
                            .unwrap();

                        let (loop_moms, jacobian) =
                            if loop_moms_f64[0].t.is_nan() || jacobian_f64.is_nan() {
                                let xs_f128 = xs.iter().map(|x| f128::f128::from(*x)).collect_vec();
                                let external_moms_f128 =
                                    external_moms.iter().map(cast_lorentz_vector).collect_vec();

                                let (loop_moms_f128, jacobian_f128) =
                                    tropical_parameterization::generate_tropical_sample(
                                        &xs_f128,
                                        &external_moms_f128,
                                        graph,
                                        self.settings.general.debug,
                                    )
                                    .unwrap();

                                let loop_moms =
                                    loop_moms_f128.iter().map(cast_lorentz_vector).collect_vec();
                                let jacobian = Into::<f64>::into(jacobian_f128);

                                (loop_moms, jacobian)
                            } else {
                                (loop_moms_f64, jacobian_f64)
                            };

                        let default_sample = DefaultSample {
                            loop_moms,
                            external_moms,
                            jacobian: jacobian * pdf,
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
    fn default_parametrize(&self, xs: &[f64]) -> DefaultSample<f64> {
        let (external_moms, pdf) = self.settings.kinematics.externals.get_externals(xs);
        let (loop_moms_vec, param_jacobian) = global_parameterize(
            xs,
            self.settings.kinematics.e_cm * self.settings.kinematics.e_cm,
            &self.settings,
            false,
        );

        let loop_moms = loop_moms_vec
            .iter()
            .map(|x| LorentzVector::from_args(0., x[0], x[1], x[2]))
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
        result: Complex<f64>,
        rotated_result: Complex<f64>,
        stability_settings: &StabilityLevelSetting,
        integrated_phase: IntegratedPhase,
        max_eval: f64,
        wgt: f64,
    ) -> (Complex<f64>, bool) {
        let average = (result + rotated_result) / 2.;

        let (
            result_for_comparison,
            rotated_result_for_comparison,
            average_for_comparison,
            max_wgt_for_comparison,
        ) = match integrated_phase {
            IntegratedPhase::Real => (result.re, rotated_result.re, average.re, max_eval),
            IntegratedPhase::Imag => (result.im, rotated_result.im, average.im, max_eval),
            IntegratedPhase::Both => unimplemented!("integrated phase both not implemented"),
        };

        let error = if result_for_comparison.is_zero() && rotated_result_for_comparison.is_zero() {
            0.
        } else {
            ((result_for_comparison - rotated_result_for_comparison) / average_for_comparison).abs()
        };

        let stable = match integrated_phase {
            IntegratedPhase::Real => error < stability_settings.required_precision_for_re,
            IntegratedPhase::Imag => error < stability_settings.required_precision_for_im,
            IntegratedPhase::Both => unimplemented!("integrated phase both not implemented"),
        };

        let below_wgt_threshold = if stability_settings.escalate_for_large_weight_threshold > 0.
            && max_wgt_for_comparison != 0.
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

    pub fn amplitude_integrand_constructor(amplitude: Amplitude, settings: Settings) -> Self {
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
    /// Rotation for stability checks
    #[inline]
    fn get_rotated_sample(
        &self,
        rotation_function: impl Fn(&LorentzVector<T>) -> LorentzVector<T>,
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
    #[inline]
    fn cast_sample<T2: FloatLike + From<T>>(&self) -> GammaLoopSample<T2> {
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
struct DefaultSample<T: FloatLike> {
    loop_moms: Vec<LorentzVector<T>>,
    external_moms: Vec<LorentzVector<T>>,
    jacobian: f64,
}

impl<T: FloatLike> DefaultSample<T> {
    #[inline]
    /// Rotation for stability checks
    fn get_rotated_sample(
        &self,
        rotation_function: impl Fn(&LorentzVector<T>) -> LorentzVector<T>,
    ) -> Self {
        Self {
            loop_moms: self.loop_moms.iter().map(&rotation_function).collect_vec(),
            external_moms: self
                .external_moms
                .iter()
                .map(&rotation_function)
                .collect_vec(),
            jacobian: self.jacobian,
        }
    }

    /// Cast the sample to a different precision
    #[inline]
    fn cast_sample<T2: FloatLike + From<T>>(&self) -> DefaultSample<T2> {
        DefaultSample {
            loop_moms: self
                .loop_moms
                .iter()
                .map(cast_lorentz_vector::<T, T2>)
                .collect_vec(),
            external_moms: self
                .external_moms
                .iter()
                .map(cast_lorentz_vector::<T, T2>)
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
    /// Rotation for stability checks
    #[inline]
    fn get_rotated_sample(
        &self,
        rotation_function: impl Fn(&LorentzVector<T>) -> LorentzVector<T>,
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
    fn cast_sample<T2: FloatLike + From<T>>(&self) -> DiscreteGraphSample<T2> {
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
fn unwrap_cont_sample(sample: &Sample<f64>) -> &[f64] {
    if let Sample::Continuous(_, xs) = sample {
        xs
    } else {
        panic!("Invalid sample structure")
    }
}

#[inline]
fn unwrap_single_discrete_sample(sample: &Sample<f64>) -> (usize, &[f64]) {
    if let Sample::Discrete(_, index, Some(cont_sample)) = sample {
        (*index, unwrap_cont_sample(cont_sample))
    } else {
        panic!("Invalid sample structure")
    }
}

#[inline]
fn unwrap_double_discrete_sample(sample: &Sample<f64>) -> (usize, (usize, &[f64])) {
    if let Sample::Discrete(_, index, Some(discrete_sample)) = sample {
        (*index, unwrap_single_discrete_sample(discrete_sample))
    } else {
        panic!("Invalid sample structure")
    }
}
