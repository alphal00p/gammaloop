use core::panic;
use std::time::Duration;

use crate::cross_section::{Amplitude, AmplitudeGraph, SuperGraph};
use crate::graph::{EdgeType, Graph, LoopMomentumBasisSpecification};
use crate::integrands::{HasIntegrand, Integrand};
use crate::integrate::UserData;
use crate::utils::{
    cast_complex, cast_lorentz_vector, format_for_compare_digits, global_parameterize, FloatLike,
};
use crate::{IntegratedPhase, Settings};
use crate::{Precision, StabilityLevelSetting};
use colored::Colorize;
use itertools::Itertools;
use lorentz_vector::LorentzVector;
use num::Complex;
use num_traits::{Inv, Zero};
use symbolica::numerical_integration::{ContinuousGrid, DiscreteGrid, Grid, Sample};

#[derive(Debug, Copy, Clone)]
pub struct Statistics {
    max_re_eval: f64,
    max_im_eval: f64,
    num_f32_evals: usize,
    num_f64_evals: usize,
    num_f128_evals: usize,
    num_arb_evals: usize,
    total_time_in_f32: Duration,
    total_time_in_f64: Duration,
    total_time_in_f128: Duration,
    total_time_in_arb: Duration,
    num_f32_unstable: usize,
    num_f64_unstable: usize,
    num_f128_unstable: usize,
    num_arb_unstable: usize,
}

impl Statistics {
    fn update(&mut self, result: Complex<f64>, duration: Duration, stable: bool, prec: Precision) {
        self.max_re_eval = self.max_re_eval.abs().max(result.re.abs());
        self.max_im_eval = self.max_im_eval.abs().max(result.im.abs());

        match prec {
            Precision::Single => {
                self.num_f32_evals += 1;
                self.total_time_in_f32 += duration;
                if !stable {
                    self.num_f32_unstable += 1;
                }
            }
            Precision::Double => {
                self.num_f64_evals += 1;
                self.total_time_in_f64 += duration;
                if !stable {
                    self.num_f64_unstable += 1;
                }
            }
            Precision::Quad => {
                self.num_f128_evals += 1;
                self.total_time_in_f128 += duration;
                if !stable {
                    self.num_f128_unstable += 1;
                }
            }
            Precision::Arb(_) => {
                self.num_arb_evals += 1;
                self.total_time_in_arb += duration;
                if !stable {
                    self.num_arb_unstable += 1;
                }
            }
        }
    }

    fn merge(&self, other: &Self) -> Self {
        Self {
            max_re_eval: self.max_re_eval.abs().max(other.max_re_eval.abs()),
            max_im_eval: self.max_im_eval.abs().max(other.max_im_eval.abs()),
            num_f32_evals: self.num_f32_evals + other.num_f32_evals,
            num_f64_evals: self.num_f64_evals + other.num_f64_evals,
            num_f128_evals: self.num_f128_evals + other.num_f128_evals,
            num_arb_evals: self.num_arb_evals + other.num_arb_evals,
            total_time_in_f32: self.total_time_in_f32 + other.total_time_in_f32,
            total_time_in_f64: self.total_time_in_f64 + other.total_time_in_f64,
            total_time_in_f128: self.total_time_in_f128 + other.total_time_in_f128,
            total_time_in_arb: self.total_time_in_arb + other.total_time_in_arb,
            num_f32_unstable: self.num_f32_unstable + other.num_f32_unstable,
            num_f64_unstable: self.num_f64_unstable + other.num_f64_unstable,
            num_f128_unstable: self.num_f128_unstable + other.num_f128_unstable,
            num_arb_unstable: self.num_arb_unstable + other.num_arb_unstable,
        }
    }

    fn new() -> Self {
        Self {
            max_re_eval: 0.,
            max_im_eval: 0.,
            num_f32_evals: 0,
            num_f64_evals: 0,
            num_f128_evals: 0,
            num_arb_evals: 0,
            total_time_in_f32: Duration::ZERO,
            total_time_in_f64: Duration::ZERO,
            total_time_in_f128: Duration::ZERO,
            total_time_in_arb: Duration::ZERO,
            num_f32_unstable: 0,
            num_f64_unstable: 0,
            num_f128_unstable: 0,
            num_arb_unstable: 0,
        }
    }
}

// trait to capture the common behaviour of amplitudes and cross sections
trait GraphIntegrand {
    fn get_graph(&self) -> &Graph;
    fn get_multi_channeling_channels(&self) -> &[usize];
    fn evaluate<T: FloatLike>(
        &self,
        channel_id: Option<usize>,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
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
    fn evaluate<T: FloatLike>(
        &self,
        channel_id: Option<usize>,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
        settings: &Settings,
    ) -> Complex<T> {
        match channel_id {
            Some(channel_id) => {
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

                let cff = self.get_graph().evaluate_cff_expression_in_lmb(
                    loop_moms,
                    external_moms,
                    &lmb_specification,
                );

                let onshell_energies = self.get_graph().compute_onshell_energies_in_lmb(
                    loop_moms,
                    external_moms,
                    lmb_specification,
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
                    .map(|x| {
                        x.powf(-Into::<T>::into(settings.multi_channeling.alpha)) * energy_product
                    })
                    .sum::<T>();

                let multichanneling_numerator = Complex::new(
                    lmb_list[channel]
                        .basis
                        .iter()
                        .fold(T::one(), |acc, i| acc * onshell_energies[*i])
                        .powf(-Into::<T>::into(settings.multi_channeling.alpha)),
                    T::zero(),
                );

                multichanneling_numerator * cff / denominator
            }
            None => {
                if settings.multi_channeling.enabled {
                    let channels_iterator = (0..get_lmb_count(self)).map(Some);

                    channels_iterator
                        .map(|channel_id| {
                            self.evaluate(channel_id, loop_moms, external_moms, settings)
                        })
                        .sum::<Complex<T>>()
                } else {
                    let energy_product = self
                        .get_graph()
                        .compute_energy_product(loop_moms, external_moms);

                    if settings.general.use_ltd {
                        let ltd = self
                            .get_graph()
                            .evaluate_ltd_expression(loop_moms, external_moms);

                        ltd / energy_product
                    } else {
                        let cff = self
                            .get_graph()
                            .evaluate_cff_expression(loop_moms, external_moms);

                        cff / energy_product
                    }
                }
            }
        }
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
    fn evaluate<T: FloatLike>(
        &self,
        channel_id: Option<usize>,
        loop_moms: &[LorentzVector<T>],
        external_moms: &[LorentzVector<T>],
        settings: &Settings,
    ) -> Complex<T> {
        // sum over cuts
        todo!()
    }
}

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

#[inline]
fn evaluate<T: GraphIntegrand, F: FloatLike>(
    graph_integrands: &[T],
    graph_id: Option<usize>,
    channel_id: Option<usize>,
    loop_moms: &[LorentzVector<F>],
    external_moms: &[LorentzVector<F>],
    settings: &Settings,
) -> Complex<F> {
    match graph_id {
        Some(graph_id) => {
            graph_integrands[graph_id].evaluate(channel_id, loop_moms, external_moms, settings)
        }
        None => graph_integrands
            .iter()
            .map(|g| g.evaluate(channel_id, loop_moms, external_moms, settings))
            .sum(),
    }
}
#[inline]
fn compute_energy_product<T: GraphIntegrand, F: FloatLike>(
    graph_integrand: &T,
    loop_moms: &[LorentzVector<F>],
    external_moms: &[LorentzVector<F>],
) -> F {
    let graph = graph_integrand.get_graph();
    graph.compute_energy_product(loop_moms, external_moms)
}

fn create_grid<T: GraphIntegrand>(graph_integrand: &T, settings: &Settings) -> Grid<f64> {
    // assume that we are using an injective parametrisation.
    let continious_dimension = get_loop_count(graph_integrand) * 3;

    let continous_grid = Grid::Continuous(ContinuousGrid::new(
        continious_dimension,
        settings.integrator.n_bins,
        settings.integrator.min_samples_for_update,
        settings.integrator.bin_number_evolution.clone(),
        settings.integrator.train_on_avg,
    ));

    if settings.multi_channeling.discrete_sampling && settings.multi_channeling.enabled {
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
    } else {
        continous_grid
    }
}

#[derive(Clone)]
enum GraphIntegrands {
    Amplitude(Vec<AmplitudeGraph>),
    CrossSection(Vec<SuperGraph>),
}

#[derive(Clone)]
pub struct GammaLoopIntegrand {
    pub settings: Settings,
    graph_integrands: GraphIntegrands,
    pub statistics: Statistics,
}

impl GraphIntegrands {
    fn create_grid(&self, settings: &Settings) -> Grid<f64> {
        if !settings.general.discrete_sample_graphs && settings.multi_channeling.discrete_sampling {
            panic!("It is impossible to use discrete sampling over channels without discrete sampling over graphs")
        }

        if settings.general.discrete_sample_graphs {
            match self {
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
            }
        } else {
            match self {
                GraphIntegrands::Amplitude(graphs) => create_grid(&graphs[0], settings),
                GraphIntegrands::CrossSection(graphs) => create_grid(&graphs[0], settings),
            }
        }
    }
}

impl HasIntegrand for GammaLoopIntegrand {
    fn create_grid(&self) -> Grid<f64> {
        self.graph_integrands.create_grid(&self.settings)
    }

    #[allow(unused_variables)]
    fn evaluate_sample(
        &mut self,
        sample: &symbolica::numerical_integration::Sample<f64>,
        wgt: f64,
        iter: usize,
        use_f128: bool,
    ) -> num::Complex<f64> {
        // setup the evaluation of the integrand in the different stability levels
        let mut results_of_stability_levels =
            Vec::with_capacity(self.settings.stability.levels.len());

        // create an iterator containing the information for evaluation at each stability level
        let stability_iterator = self.create_stability_iterator(use_f128);

        let sample_point = self.unpack_sample(sample);

        // generate external momenta, the get_exeternals function now just returns a constant value set by the user,
        // but it may be generalized to a full fledged pdf later.
        let (external_moms, pdf_weight) = self
            .settings
            .kinematics
            .externals
            .get_externals(sample_point.sample);

        let (loop_moms, jacobian) = self.parameterize(sample_point.sample);

        // rotate the momenta for the stability tests.
        let rotation_method = self.settings.stability.rotation_axis.rotation_function();
        let rotated_loop_moms = loop_moms.iter().map(&rotation_method).collect_vec();
        let rotated_external_moms = external_moms.iter().map(&rotation_method).collect_vec();

        // iterate over the stability levels, break if the point is stable
        for stability_level in stability_iterator {
            // evaluate the integrand at the current stability level
            let (result, rotated_result, duration) = self.evaluate_at_prec(
                sample_point.graph,
                sample_point.channel,
                &loop_moms,
                &external_moms,
                &rotated_loop_moms,
                &rotated_external_moms,
                stability_level.precision,
            );

            // check for the stability
            let (avg_result, stable) = self.stability_check(
                result,
                rotated_result,
                stability_level,
                self.settings.integrator.integrated_phase,
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

        let (most_reliable_result, stable, precision, duration) = results_of_stability_levels.last().unwrap_or_else(
            || panic!("No evaluation was done, perhaps the final stability level has a non-negative escalation threshold?")
        );

        let res = most_reliable_result * jacobian * pdf_weight;
        self.statistics.update(res, *duration, *stable, *precision);

        // 1 / (2 pi )^L
        let prefactor = self.compute_2pi_factor().inv();

        if self.settings.general.debug > 1 {
            println!("{}", "
            |  DEBUG -------------------------------------------------------------------------------------  |".green());
            println!();

            println!("{}", "sample form havana: ".blue());
            println!("\tsampled x: {:?}", sample_point);
            println!();

            println!("{}", "parametrisation result".blue());
            println!("{}", "\tloop momenta: ".yellow());

            for (index, loop_mom) in loop_moms.iter().enumerate() {
                println!("\t\tloop momentum {}: {:?}", index, loop_mom);
            }

            println!("{}", "\trotated loop momenta: ".yellow());

            for (index, loop_mom) in rotated_loop_moms.iter().enumerate() {
                println!("\t\tloop momentum {}: {:?}", index, loop_mom);
            }

            println!("{}", "\texternal momenta: ".yellow());

            for (index, external_mom) in external_moms.iter().enumerate() {
                println!("\t\texternal momentum {}: {:?}", index, external_mom);
            }

            println!("{}", "\trotated external momenta: ".yellow());

            for (index, external_mom) in rotated_external_moms.iter().enumerate() {
                println!("\t\texternal momentum {}: {:?}", index, external_mom);
            }

            println!("\t{}: {:+e}", "jacobian".yellow(), jacobian);
            println!();

            println!("{}", "evaluation result: ".blue());
            println!(
                "{}: {:+e}",
                "\tcff expression: ".yellow(),
                most_reliable_result
            );

            println!("\t{}: {:+e}", "result".yellow(), res * prefactor);
        }

        res * prefactor
    }

    fn get_event_manager_mut(&mut self) -> &mut crate::observables::EventManager {
        todo!()
    }

    // this maybe needs to change to Vec<usize>
    fn get_n_dim(&self) -> usize {
        match &self.graph_integrands {
            GraphIntegrands::Amplitude(graphs) => get_loop_count(&graphs[0]) * 3,
            GraphIntegrands::CrossSection(graphs) => get_loop_count(&graphs[0]) * 3,
        }
    }

    fn merge_results<I: HasIntegrand>(&mut self, _other: &mut I, _iter: usize) {}

    fn update_results(&mut self, _iter: usize) {}
}

#[allow(clippy::too_many_arguments)]
impl GammaLoopIntegrand {
    #[inline]
    pub fn evaluate_at_prec(
        &self,
        graph_id: Option<usize>,
        channel_id: Option<usize>,
        loop_moms: &[LorentzVector<f64>],
        external_moms: &[LorentzVector<f64>],
        rotated_loop_moms: &[LorentzVector<f64>],
        rotated_external_moms: &[LorentzVector<f64>],
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
                    let result = evaluate(
                        graph_integrands,
                        graph_id,
                        channel_id,
                        loop_moms,
                        external_moms,
                        &self.settings,
                    );

                    let rotated_result = evaluate(
                        graph_integrands,
                        graph_id,
                        channel_id,
                        rotated_loop_moms,
                        rotated_external_moms,
                        &self.settings,
                    );

                    (result, rotated_result)
                }
                GraphIntegrands::CrossSection(graph_integrands) => {
                    let result = evaluate(
                        graph_integrands,
                        graph_id,
                        channel_id,
                        loop_moms,
                        external_moms,
                        &self.settings,
                    );

                    let rotated_result = evaluate(
                        graph_integrands,
                        graph_id,
                        channel_id,
                        rotated_loop_moms,
                        rotated_external_moms,
                        &self.settings,
                    );

                    (result, rotated_result)
                }
            },
            Precision::Quad => {
                let loop_moms_f128 = loop_moms
                    .iter()
                    .map(cast_lorentz_vector::<f64, f128::f128>)
                    .collect_vec();
                let rotated_loop_moms_f128 = rotated_loop_moms
                    .iter()
                    .map(cast_lorentz_vector::<f64, f128::f128>)
                    .collect_vec();
                let external_moms_f128 = external_moms
                    .iter()
                    .map(cast_lorentz_vector::<f64, f128::f128>)
                    .collect_vec();
                let rotated_external_moms_f128 = rotated_external_moms
                    .iter()
                    .map(cast_lorentz_vector::<f64, f128::f128>)
                    .collect_vec();

                match &self.graph_integrands {
                    GraphIntegrands::Amplitude(graph_integrands) => {
                        let result = evaluate(
                            graph_integrands,
                            graph_id,
                            channel_id,
                            &loop_moms_f128,
                            &external_moms_f128,
                            &self.settings,
                        );

                        let rotated_result = evaluate(
                            graph_integrands,
                            graph_id,
                            channel_id,
                            &rotated_loop_moms_f128,
                            &rotated_external_moms_f128,
                            &self.settings,
                        );

                        (cast_complex(result), cast_complex(rotated_result))
                    }
                    GraphIntegrands::CrossSection(graph_integrands) => {
                        let result = evaluate(
                            graph_integrands,
                            graph_id,
                            channel_id,
                            &loop_moms_f128,
                            &external_moms_f128,
                            &self.settings,
                        );

                        let rotated_result = evaluate(
                            graph_integrands,
                            graph_id,
                            channel_id,
                            &rotated_loop_moms_f128,
                            &rotated_external_moms_f128,
                            &self.settings,
                        );

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

    pub fn user_data_generator(&self, num_cores: usize, _settings: &Settings) -> UserData {
        UserData {
            integrand: vec![Integrand::GammaLoopIntegrand(self.clone()); num_cores],
        }
    }

    #[inline]
    fn unpack_sample<'a>(&self, sample: &'a Sample<f64>) -> SamplePoint<'a> {
        match sample {
            Sample::Discrete(_, graph_index, Some(nested_sample)) => match nested_sample.as_ref() {
                Sample::Discrete(_, channel_index, Some(cont_sample)) => {
                    let xs = if let Sample::Continuous(_, cont_sample) = cont_sample.as_ref() {
                        cont_sample
                    } else {
                        panic!("invalid sample structure")
                    };

                    SamplePoint {
                        graph: Some(*graph_index),
                        channel: Some(*channel_index),
                        sample: xs,
                    }
                }
                Sample::Continuous(_, xs) => SamplePoint {
                    graph: Some(*graph_index),
                    channel: None,
                    sample: xs,
                },
                _ => panic!("Invalid sample structure"),
            },
            Sample::Continuous(_, xs) => SamplePoint {
                graph: None,
                channel: None,
                sample: xs,
            },
            _ => panic!("Invalid sample structure"),
        }
    }

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
                accepted_radius_in_x_range: [0., 1.],
            }]
            .iter()
        } else {
            self.settings.stability.levels.iter()
        }
    }

    #[inline]
    fn parameterize(&self, x_space_point: &[f64]) -> (Vec<LorentzVector<f64>>, f64) {
        // perform the conformal map, may be generalized to include tropical sampling later.
        let (loop_moms, jacobian) = global_parameterize(
            x_space_point,
            self.settings.kinematics.e_cm * self.settings.kinematics.e_cm,
            &self.settings,
            false,
        );

        (
            loop_moms
                .into_iter()
                .map(|p| LorentzVector::from_args(0.0, p[0], p[1], p[2]))
                .collect_vec(),
            jacobian,
        )
    }

    // compute thea average and check the accuracy of the result
    #[inline]
    fn stability_check(
        &self,
        result: Complex<f64>,
        rotated_result: Complex<f64>,
        stability_settings: &StabilityLevelSetting,
        integrated_phase: IntegratedPhase,
    ) -> (Complex<f64>, bool) {
        let average = (result + rotated_result) / 2.;

        let (
            result_for_comparison,
            rotated_result_for_comparison,
            average_for_comparison,
            max_wgt_for_comparison,
        ) = match integrated_phase {
            IntegratedPhase::Real => (
                result.re,
                rotated_result.re,
                average.re,
                self.statistics.max_re_eval,
            ),
            IntegratedPhase::Imag => (
                result.im,
                rotated_result.im,
                average.im,
                self.statistics.max_im_eval,
            ),
            IntegratedPhase::Both => unimplemented!("integrated phase both not implemented"),
        };

        let error = if result_for_comparison.is_zero() && rotated_result_for_comparison.is_zero() {
            0.
        } else {
            ((result_for_comparison - rotated_result_for_comparison) / average_for_comparison).abs()
        };

        let stable = error < stability_settings.required_precision_for_re;

        let below_wgt_threshold = if stability_settings.escalate_for_large_weight_threshold > 0. {
            average_for_comparison.abs()
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
            statistics: Statistics::new(),
        }
    }
}

#[derive(Debug, Clone)]
struct SamplePoint<'a> {
    graph: Option<usize>,
    channel: Option<usize>,
    sample: &'a [f64],
}
