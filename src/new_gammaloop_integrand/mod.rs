use std::cell::RefCell;

use crate::evaluation_result::{EvaluationMetaData, EvaluationResult};
use crate::integrands::{HasIntegrand, Integrand};
use crate::integrate::UserData;
use crate::momentum::Rotation;
use crate::momentum_sample::{BareMomentumSample, LoopMomenta, MomentumSample};
use crate::new_graph::{FeynmanGraph, Graph, LmbIndex, LoopMomentumBasis};
use crate::utils::{format_for_compare_digits, get_n_dim_for_n_loop_momenta, FloatLike, F};
use bincode_trait_derive::{Decode, Encode};
use colored::Colorize;
use derive_more::{From, Into};
use enum_dispatch::enum_dispatch;
use gammaloop_sample::{parameterize, DiscreteGraphSample, GammaLoopSample};
use itertools::Itertools;
use momtrop::float::MomTropFloat;
use momtrop::SampleGenerator;
use serde::{Deserialize, Serialize};
use spenso::algebra::algebraic_traits::IsZero;
use spenso::algebra::complex::Complex;
use spenso::tensors::parametric::SerializableCompiledEvaluator;
use std::time::Duration;
use symbolica::evaluate::ExpressionEvaluator;
use symbolica::numerical_integration::{ContinuousGrid, DiscreteGrid, Grid, Sample};
use typed_index_collections::TiVec;
pub mod amplitude_integrand;
pub mod cross_section_integrand;
pub mod gammaloop_sample;
use crate::observables::EventManager;
use crate::utils::f128;
use crate::{
    DependentMomentaConstructor, DiscreteGraphSamplingSettings, DiscreteGraphSamplingType,
    IntegratorSettings, Polarizations, Precision, SamplingSettings, Settings,
    StabilityLevelSetting, StabilitySettings,
};

#[derive(Clone)]
#[enum_dispatch(HasIntegrand)]
pub enum NewIntegrand {
    Amplitude(amplitude_integrand::AmplitudeIntegrand),
    CrossSection(cross_section_integrand::CrossSectionIntegrand),
}

impl NewIntegrand {
    pub fn get_settings(&self) -> &Settings {
        match self {
            NewIntegrand::Amplitude(integrand) => &integrand.settings,
            NewIntegrand::CrossSection(integrand) => &integrand.settings,
        }
    }

    /// Used to create the use_data_generator closure for havana_integrate
    pub fn user_data_generator(&self, num_cores: usize, _settings: &Settings) -> UserData {
        UserData {
            integrand: vec![Integrand::NewIntegrand(self.clone()); num_cores],
        }
    }

    pub fn get_mut_settings(&mut self) -> &mut Settings {
        match self {
            NewIntegrand::Amplitude(integrand) => &mut integrand.settings,
            NewIntegrand::CrossSection(integrand) => &mut integrand.settings,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub enum IntegrandType {
    Amplitude,
    CrossSection,
}

fn create_stability_iterator(
    settings: &StabilitySettings,
    use_f128: bool,
) -> Vec<StabilityLevelSetting> {
    if use_f128 {
        // overwrite the stability settings if use_f128 is enabled, but attempt to use user defined settings for f128
        if let Some(f128_settings_position) =
            settings.levels.iter().position(|stability_level_setting| {
                stability_level_setting.precision == Precision::Quad
            })
        {
            vec![settings.levels[f128_settings_position]]
        } else {
            vec![StabilityLevelSetting {
                precision: Precision::Quad,
                required_precision_for_re: F(1e-5),
                required_precision_for_im: F(1e-5),
                escalate_for_large_weight_threshold: F(-1.),
            }]
        }
    } else {
        settings.levels.clone()
    }
}

#[inline]
fn stability_check(
    settings: &Settings,
    results: &[Complex<F<f64>>],
    stability_settings: &StabilityLevelSetting,
    max_eval: Complex<F<f64>>,
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

    if settings.general.debug > 1 {
        if let Some(unstable_index) = unstable_sample {
            let unstable_point = results[unstable_index];
            let rotation_axis = format!("{:?}", settings.stability.rotation_axis[unstable_index]);

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
        && max_eval.is_non_zero()
    {
        average.re.abs() * wgt
            < stability_settings.escalate_for_large_weight_threshold * max_eval.re
            || average.im.abs() * wgt
                < stability_settings.escalate_for_large_weight_threshold * max_eval.im
    } else {
        true
    };

    (average, stable && below_wgt_threshold)
}

#[derive(Clone, Encode, Decode)]
pub struct GenericEvaluator {
    pub f64_compiled: Option<RefCell<SerializableCompiledEvaluator>>,
    pub f64_eager: RefCell<ExpressionEvaluator<F<f64>>>,
    pub f128: RefCell<ExpressionEvaluator<F<f128>>>,
}

pub trait GenericEvaluatorFloat<T: FloatLike = Self> {
    fn get_evaluator(generic_evaluator: &GenericEvaluator) -> impl Fn(&[F<T>]) -> F<T>;
}

impl GenericEvaluatorFloat for f64 {
    fn get_evaluator(generic_evaluator: &GenericEvaluator) -> impl Fn(&[F<f64>]) -> F<f64> {
        |params: &[F<f64>]| {
            let mut out = vec![F(0.)];
            if let Some(compiled) = &generic_evaluator.f64_compiled {
                compiled.borrow_mut().evaluate(params, &mut out);
            } else {
                generic_evaluator
                    .f64_eager
                    .borrow_mut()
                    .evaluate(params, &mut out);
            }

            out[0]
        }
    }
}

impl GenericEvaluatorFloat for f128 {
    fn get_evaluator(generic_evaluator: &GenericEvaluator) -> impl Fn(&[F<f128>]) -> F<f128> {
        |params: &[F<f128>]| {
            let mut out = vec![params[0].zero()];
            generic_evaluator
                .f128
                .borrow_mut()
                .evaluate(params, &mut out);

            out[0].clone()
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct StabilityLevelResult {
    pub result: Complex<F<f64>>,
    pub stability_level_used: Precision,
    pub parameterization_time: Duration,
    pub ltd_evaluation_time: Duration,
    pub is_stable: bool,
}

#[derive(
    Debug, Clone, Copy, From, Into, PartialEq, Eq, Ord, PartialOrd, Hash, Serialize, Deserialize,
)]
pub struct ChannelIndex(usize);

/// Helper struct for the LMB multi-channeling setup
#[derive(Debug, Clone, Serialize, Deserialize, Encode, Decode)]
pub struct LmbMultiChannelingSetup {
    pub channels: TiVec<ChannelIndex, LmbIndex>,
}

impl LmbMultiChannelingSetup {
    fn reinterpret_loop_momenta_impl<T: FloatLike>(
        &self,
        channel_index: ChannelIndex,
        momentum_sample: &BareMomentumSample<T>,
        base_lmb: &LoopMomentumBasis,
        all_bases: &TiVec<LmbIndex, LoopMomentumBasis>,
    ) -> BareMomentumSample<T> {
        let channel_lmb = &all_bases[self.channels[channel_index]];
        let new_loop_moms: LoopMomenta<F<T>> = base_lmb
            .basis
            .iter()
            .map(|&edge_index| {
                let signature_of_edge_channel_lmb = &channel_lmb.edge_signatures[edge_index];

                signature_of_edge_channel_lmb
                    .internal
                    .apply_typed(&momentum_sample.loop_moms)
                    + signature_of_edge_channel_lmb
                        .external
                        .apply(&momentum_sample.external_moms.raw)
                        .spatial
            })
            .collect();

        BareMomentumSample {
            loop_moms: new_loop_moms,
            external_moms: momentum_sample.external_moms.clone(),
            polarizations: momentum_sample.polarizations.clone(),
            jacobian: momentum_sample.jacobian.clone(),
            orientation: momentum_sample.orientation,
        }
    }

    pub fn reinterpret_loop_momenta_and_compute_prefactor_all_channels<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        graph: &Graph,
        all_bases: &TiVec<LmbIndex, LoopMomentumBasis>,
        alpha: &F<T>,
    ) -> TiVec<ChannelIndex, (MomentumSample<T>, F<T>)> {
        self.channels
            .iter_enumerated()
            .map(|(channel_index, _)| {
                self.reinterpret_loop_momenta_and_compute_prefactor(
                    channel_index,
                    momentum_sample,
                    graph,
                    all_bases,
                    alpha,
                )
            })
            .collect()
    }

    /// This function is used to do do LMB multi-channeling without fully switching to a different lmb
    /// for each channel. The momenta provided are reinterpreted as loop momenta of the lmb corresponding to the channel_index.
    /// Then we transform these loop momenta to the fixed lmb of the graph. The prefactor is immediately computed for the requested channel
    pub fn reinterpret_loop_momenta_and_compute_prefactor<T: FloatLike>(
        &self,
        channel_index: ChannelIndex,
        momentum_sample: &MomentumSample<T>,
        graph: &Graph,
        all_bases: &TiVec<LmbIndex, LoopMomentumBasis>,
        alpha: &F<T>,
    ) -> (MomentumSample<T>, F<T>) {
        let base_lmb = &graph.loop_momentum_basis;

        let sample = MomentumSample {
            sample: self.reinterpret_loop_momenta_impl(
                channel_index,
                &momentum_sample.sample,
                base_lmb,
                all_bases,
            ),
            rotated_sample: momentum_sample
                .rotated_sample
                .as_ref()
                .map(|rotated_sample| {
                    self.reinterpret_loop_momenta_impl(
                        channel_index,
                        rotated_sample,
                        base_lmb,
                        all_bases,
                    )
                }),
            uuid: momentum_sample.uuid,
        };

        let prefactor =
            self.compute_prefactor_impl(channel_index, &sample, graph, all_bases, alpha);

        (sample, prefactor)
    }

    /// Computes the prefactor for the given channel index and momentum sample.
    fn compute_prefactor_impl<T: FloatLike>(
        &self,
        channel_index: ChannelIndex,
        momentum_sample: &MomentumSample<T>,
        graph: &Graph,
        all_bases: &TiVec<LmbIndex, LoopMomentumBasis>,
        alpha: &F<T>,
    ) -> F<T> {
        let all_energies = graph.underlying.get_energy_cache(
            &momentum_sample.sample.loop_moms,
            &momentum_sample.sample.external_moms,
            &graph.loop_momentum_basis,
        );

        let mut numerator = momentum_sample.zero();

        let denominators = self
            .channels
            .iter()
            .map(|&lmb_index| {
                let channel_product = all_bases[lmb_index]
                    .basis
                    .iter()
                    .map(|&edge_index| &all_energies[edge_index])
                    .fold(momentum_sample.one(), |product, energy| product * energy)
                    .powf(&-alpha);

                if self.channels[channel_index] == lmb_index {
                    numerator = channel_product.clone();
                }

                channel_product
            })
            .fold(momentum_sample.zero(), |sum, summand| sum + summand);

        numerator / denominators
    }
}

pub trait GammaloopIntegrand {
    type G: GraphTerm;
    fn get_rotations(&self) -> impl Iterator<Item = &Rotation>;
    fn get_terms(&self) -> impl Iterator<Item = &Self::G>;
    fn get_settings(&self) -> &Settings;
    fn get_graph(&self, graph_id: usize) -> &Self::G;
    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor;
    fn get_polarizations(&self) -> &[Polarizations];
    fn get_model_parameter_cache<T: FloatLike>(&self) -> Vec<Complex<F<T>>>;
}

fn get_global_dimension_if_exists<I: GammaloopIntegrand>(integrand: &I) -> Option<usize> {
    if integrand
        .get_settings()
        .sampling
        .get_parameterization_settings()
        .is_none()
    {
        None
    } else {
        Some(
            integrand
                .get_graph(0)
                .get_graph()
                .underlying
                .get_loop_number()
                * 3,
        )
    }
}

pub trait GraphTerm {
    fn evaluate<T: FloatLike>(
        &self,
        sample: &MomentumSample<T>,
        settings: &Settings,
        model_parameter_cache: &[Complex<F<T>>],
    ) -> Complex<F<T>>;

    fn get_multi_channeling_setup(&self) -> &LmbMultiChannelingSetup;
    fn get_graph(&self) -> &Graph;
    fn get_lmbs(&self) -> &TiVec<LmbIndex, LoopMomentumBasis>;
    fn get_num_orientations(&self) -> usize;
    fn get_tropical_sampler(&self) -> &SampleGenerator<3>;
}

fn evaluate_all_rotations<T: FloatLike, I: GammaloopIntegrand>(
    integrand: &I,
    gammaloop_sample: &GammaLoopSample<T>,
) -> (Vec<Complex<F<f64>>>, Duration) {
    // rotate the momenta for the stability tests.
    let gammaloop_samples: Vec<_> = integrand
        .get_rotations()
        .map(|rotation| gammaloop_sample.get_rotated_sample(rotation))
        .collect();

    let start_time = std::time::Instant::now();
    let evaluation_results = gammaloop_samples
        .iter()
        .map(|gammaloop_sample| evaluate_single_rotation(integrand, gammaloop_sample))
        .collect_vec();
    let duration = start_time.elapsed() / gammaloop_samples.len() as u32;

    (evaluation_results, duration)
}

fn evaluate_single_rotation<T: FloatLike, I: GammaloopIntegrand>(
    integrand: &I,
    gammaloop_sample: &GammaLoopSample<T>,
) -> Complex<F<f64>> {
    let result = match &gammaloop_sample {
        GammaLoopSample::Default(sample) => integrand
            .get_terms()
            .map(|term: &I::G| {
                term.evaluate(
                    sample,
                    integrand.get_settings(),
                    &integrand.get_model_parameter_cache(),
                )
            })
            .fold(
                Complex::new_re(gammaloop_sample.get_default_sample().zero()),
                |sum, term| sum + term,
            ),
        GammaLoopSample::MultiChanneling { alpha, sample } => integrand
            .get_terms()
            .map(|term: &I::G| {
                let channels_samples = term
                    .get_multi_channeling_setup()
                    .reinterpret_loop_momenta_and_compute_prefactor_all_channels(
                        sample,
                        term.get_graph(),
                        term.get_lmbs(),
                        alpha,
                    );

                channels_samples
                    .into_iter()
                    .map(|(reparameterized_sample, prefactor)| {
                        Complex::new_re(prefactor)
                            * term.evaluate(
                                &reparameterized_sample,
                                integrand.get_settings(),
                                &integrand.get_model_parameter_cache(),
                            )
                    })
                    .fold(
                        Complex::new_re(gammaloop_sample.get_default_sample().zero()),
                        |sum, term| sum + term,
                    )
            })
            .fold(
                Complex::new_re(gammaloop_sample.get_default_sample().zero()),
                |sum, term| sum + term,
            ),
        GammaLoopSample::DiscreteGraph { graph_id, sample } => {
            let graph_term = integrand.get_graph(*graph_id);
            match sample {
                DiscreteGraphSample::Default(sample) => graph_term.evaluate(
                    sample,
                    integrand.get_settings(),
                    &integrand.get_model_parameter_cache(),
                ),
                DiscreteGraphSample::DiscreteMultiChanneling {
                    alpha,
                    channel_id,
                    sample,
                } => {
                    let (reparameterized_sample, prefactor) = graph_term
                        .get_multi_channeling_setup()
                        .reinterpret_loop_momenta_and_compute_prefactor(
                            *channel_id,
                            sample,
                            graph_term.get_graph(),
                            graph_term.get_lmbs(),
                            alpha,
                        );

                    Complex::new_re(prefactor)
                        * graph_term.evaluate(
                            &reparameterized_sample,
                            integrand.get_settings(),
                            &integrand.get_model_parameter_cache(),
                        )
                }
                DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                    let channel_samples = graph_term
                        .get_multi_channeling_setup()
                        .reinterpret_loop_momenta_and_compute_prefactor_all_channels(
                            sample,
                            graph_term.get_graph(),
                            graph_term.get_lmbs(),
                            alpha,
                        );

                    channel_samples
                        .into_iter()
                        .map(|(reparameterized_sample, prefactor)| {
                            Complex::new_re(prefactor)
                                * graph_term.evaluate(
                                    &reparameterized_sample,
                                    integrand.get_settings(),
                                    &integrand.get_model_parameter_cache(),
                                )
                        })
                        .fold(
                            Complex::new_re(gammaloop_sample.get_default_sample().zero()),
                            |sum, term| sum + term,
                        )
                }
                DiscreteGraphSample::Tropical(sample) => {
                    let energy_cache = graph_term.get_graph().underlying.get_energy_cache(
                        &sample.loop_moms(),
                        &sample.external_moms(),
                        &graph_term.get_graph().loop_momentum_basis,
                    );

                    let prefactor = graph_term
                        .get_graph()
                        .iter_loop_edges()
                        .map(|(_, edge_index, _)| edge_index)
                        .zip(graph_term.get_tropical_sampler().iter_edge_weights())
                        .fold(sample.one(), |product, (edge_id, weight)| {
                            let energy = &energy_cache[edge_id];
                            let edge_weight = weight;
                            product * energy.powf(&F::from_f64(2. * edge_weight))
                        });

                    Complex::new_re(prefactor)
                        * graph_term.evaluate(
                            sample,
                            integrand.get_settings(),
                            &integrand.get_model_parameter_cache(),
                        )
                }
            }
        }
    };

    let f_t_result = result * gammaloop_sample.get_default_sample().jacobian();
    Complex::new(f_t_result.re.into_ff64(), f_t_result.im.into_ff64())
}

fn create_grid_for_graph<G: GraphTerm>(
    graph_term: &G,
    settings: &DiscreteGraphSamplingSettings,
    integrator_settings: &IntegratorSettings,
) -> Grid<F<f64>> {
    match &settings.sampling_type {
        DiscreteGraphSamplingType::Default(_) | DiscreteGraphSamplingType::MultiChanneling(_) => {
            let continuous_grid = create_default_continous_grid(graph_term, integrator_settings);

            if settings.sample_orientations {
                let continuous_grids = (0..graph_term.get_num_orientations())
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
            let continuous_grid = create_default_continous_grid(graph_term, integrator_settings);
            let lmb_channel_grid = Grid::Discrete(DiscreteGrid::new(
                graph_term
                    .get_multi_channeling_setup()
                    .channels
                    .iter()
                    .map(|_| Some(continuous_grid.clone()))
                    .collect_vec(),
                integrator_settings.max_prob_ratio,
                integrator_settings.train_on_avg,
            ));

            if settings.sample_orientations {
                Grid::Discrete(DiscreteGrid::new(
                    (0..graph_term.get_num_orientations())
                        .map(|_| Some(lmb_channel_grid.clone()))
                        .collect(),
                    integrator_settings.max_prob_ratio,
                    integrator_settings.train_on_avg,
                ))
            } else {
                lmb_channel_grid
            }
        }

        DiscreteGraphSamplingType::TropicalSampling(_) => {
            let dimension = get_n_dim_for_n_loop_momenta(
                &SamplingSettings::DiscreteGraphs(settings.clone()),
                graph_term.get_graph().underlying.get_loop_number(),
                false,
                Some(graph_term.get_graph().iter_loop_edges().count()),
            );

            let continious_grid = Grid::Continuous(ContinuousGrid::new(
                dimension,
                integrator_settings.n_bins,
                integrator_settings.min_samples_for_update,
                integrator_settings.bin_number_evolution.clone(),
                integrator_settings.train_on_avg,
            ));

            if settings.sample_orientations {
                let continuous_grids = (0..graph_term.get_num_orientations())
                    .map(|_| Some(continious_grid.clone()))
                    .collect();

                Grid::Discrete(DiscreteGrid::new(
                    continuous_grids,
                    integrator_settings.max_prob_ratio,
                    integrator_settings.train_on_avg,
                ))
            } else {
                continious_grid
            }
        }
    }
}

fn create_default_continous_grid<G: GraphTerm>(
    graph_term: &G,
    integrator_settings: &IntegratorSettings,
) -> Grid<F<f64>> {
    Grid::Continuous(ContinuousGrid::new(
        graph_term.get_graph().underlying.get_loop_number() * 3,
        integrator_settings.n_bins,
        integrator_settings.min_samples_for_update,
        integrator_settings.bin_number_evolution.clone(),
        integrator_settings.train_on_avg,
    ))
}

fn create_grid<I: GammaloopIntegrand>(integrand: &I) -> Grid<F<f64>> {
    let settings = integrand.get_settings();
    match &settings.sampling {
        SamplingSettings::Default(_) => Grid::Continuous(ContinuousGrid::new(
            get_global_dimension_if_exists(integrand).unwrap(),
            settings.integrator.n_bins,
            settings.integrator.min_samples_for_update,
            settings.integrator.bin_number_evolution.clone(),
            settings.integrator.train_on_avg,
        )),
        SamplingSettings::MultiChanneling(_) => Grid::Continuous(ContinuousGrid::new(
            get_global_dimension_if_exists(integrand).unwrap(),
            settings.integrator.n_bins,
            settings.integrator.min_samples_for_update,
            settings.integrator.bin_number_evolution.clone(),
            settings.integrator.train_on_avg,
        )),
        SamplingSettings::DiscreteGraphs(discrete_graph_sampling_settings) => {
            Grid::Discrete(DiscreteGrid::new(
                integrand
                    .get_terms()
                    .map(|term| {
                        Some(create_grid_for_graph(
                            term,
                            discrete_graph_sampling_settings,
                            &settings.integrator,
                        ))
                    })
                    .collect(),
                settings.integrator.max_prob_ratio,
                settings.integrator.train_on_avg,
            ))
        }
    }
}

fn evaluate_sample<I: GammaloopIntegrand>(
    integrand: &I,
    sample: &Sample<F<f64>>,
    wgt: F<f64>,
    _iter: usize,
    use_f128: bool,
    max_eval: Complex<F<f64>>,
) -> EvaluationResult {
    let start_eval = std::time::Instant::now();
    let stability_iterator =
        create_stability_iterator(&integrand.get_settings().stability, use_f128);

    let mut results_of_stability_levels = Vec::with_capacity(stability_iterator.len());

    for stability_level in stability_iterator.into_iter() {
        let before_parameterization = std::time::Instant::now();
        let ((results, ltd_evaluation_time), parameterization_time) =
            match stability_level.precision {
                Precision::Double => {
                    if let Ok(gammaloop_sample) = parameterize::<f64, I>(sample, integrand) {
                        let parameterization_time = before_parameterization.elapsed();
                        (
                            evaluate_all_rotations(integrand, &gammaloop_sample),
                            parameterization_time,
                        )
                    } else {
                        continue;
                    }
                }
                Precision::Quad => {
                    if let Ok(gammaloop_sample) = parameterize::<f128, I>(sample, integrand) {
                        let parameterization_time = before_parameterization.elapsed();
                        (
                            evaluate_all_rotations(integrand, &gammaloop_sample),
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

        let (average_result, is_stable) = stability_check(
            integrand.get_settings(),
            &results,
            &stability_level,
            max_eval,
            wgt,
        );

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
        let re_is_nan = stability_level_result.result.re.is_nan()
            || stability_level_result.result.re.is_infinite();
        let im_is_nan = stability_level_result.result.im.is_nan()
            || stability_level_result.result.im.is_infinite();
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
