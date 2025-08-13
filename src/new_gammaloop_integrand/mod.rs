use std::borrow::Cow;
use std::cell::RefCell;
use std::collections::HashMap;
use std::fmt::Display;
use std::fs;
use std::ops::Deref;
use std::path::Path;

use crate::evaluation_result::{EvaluationMetaData, EvaluationResult};
use crate::integrands::{HasIntegrand, Integrand};
use crate::integrate::UserData;
use crate::model::Model;
use crate::momentum::{Helicity, PolType, Rotation};
use crate::momentum_sample::{
    self, BareMomentumSample, ExternalFourMomenta, ExternalIndex, LoopMomenta, MomentumSample,
};
use crate::new_graph::{FeynmanGraph, Graph, LmbIndex, LoopMomentumBasis};
use crate::numerator::ParsingNet;
use crate::utils::{
    format_for_compare_digits, get_n_dim_for_n_loop_momenta, FloatLike, PrecisionUpgradable,
    ToCoefficient, F, GS, TENSORLIB,
};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::owo_colors::OwoColorize;
use colored::Colorize;
use derive_more::{From, Into};
use enum_dispatch::enum_dispatch;
use eyre::Context;
use gammaloop_sample::{parameterize, DiscreteGraphSample, GammaLoopSample};
use itertools::Itertools;
use linnet::half_edge::involution::HedgePair;
use log::debug;
use momtrop::float::MomTropFloat;
use momtrop::SampleGenerator;
use serde::{Deserialize, Serialize};
use spenso::algebra::algebraic_traits::IsZero;
use spenso::algebra::complex::Complex;
use spenso::iterators::IteratableTensor;
use spenso::network::ExecutionResult;
use spenso::structure::concrete_index::ExpandedIndex;
use spenso::tensors::parametric::{AtomViewOrConcrete, SerializableCompiledEvaluator};
use std::time::Duration;
use symbolica::atom::{Atom, AtomCore, FunctionBuilder, Symbol};
use symbolica::domains::rational::Rational;
use symbolica::evaluate::{
    CompileOptions, EvaluationFn, ExpressionEvaluator, FunctionMap, InlineASM, OptimizationSettings,
};
use symbolica::id::Replacement;
use symbolica::numerical_integration::{ContinuousGrid, DiscreteGrid, Grid, Sample};
use symbolica::parse;
use tabled::settings::Style;
use typed_index_collections::TiVec;
pub mod amplitude_integrand;
pub mod cross_section_integrand;
pub mod gammaloop_sample;
use crate::observables::EventManager;
use crate::utils::f128;
use crate::{
    DependentMomentaConstructor, DiscreteGraphSamplingSettings, DiscreteGraphSamplingType,
    GammaLoopContext, IntegratorSettings, Polarizations, Precision, SamplingSettings, Settings,
    StabilityLevelSetting, StabilitySettings,
};
use color_eyre::Result;

const HARD_CODED_M_UV: F<f64> = F(1000.0);
const HARD_CODED_M_R_SQ: F<f64> = F(1000.0);

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
#[enum_dispatch(HasIntegrand)]
pub enum NewIntegrand {
    Amplitude(amplitude_integrand::AmplitudeIntegrand),
    CrossSection(cross_section_integrand::CrossSectionIntegrand),
}

impl NewIntegrand {
    pub(crate) fn save(&self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let path = path.as_ref().join("integrand");

        let r = fs::create_dir_all(&path).with_context(|| {
            format!(
                "Trying to create directory to save amplitude {}",
                path.display()
            )
        });
        if override_existing {
            r?;
        }
        match self {
            NewIntegrand::Amplitude(integrand) => integrand.save(path, override_existing),
            NewIntegrand::CrossSection(integrand) => integrand.save(path, override_existing),
        }
    }

    pub(crate) fn get_settings(&self) -> &Settings {
        match self {
            NewIntegrand::Amplitude(integrand) => &integrand.settings,
            NewIntegrand::CrossSection(integrand) => &integrand.settings,
        }
    }

    /// Used to create the use_data_generator closure for havana_integrate
    pub(crate) fn user_data_generator(&self, num_cores: usize, _settings: &Settings) -> UserData {
        UserData {
            integrand: vec![Integrand::NewIntegrand(self.clone()); num_cores],
        }
    }

    pub(crate) fn get_mut_settings(&mut self) -> &mut Settings {
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
#[trait_decode(trait = GammaLoopContext)]
pub struct GenericEvaluator {
    pub expr: Atom,
    pub rational:
        RefCell<Option<ExpressionEvaluator<symbolica::domains::float::Complex<Rational>>>>,
    pub f64_compiled: Option<RefCell<SerializableCompiledEvaluator>>,
    pub f64_eager: RefCell<ExpressionEvaluator<Complex<F<f64>>>>,
    pub f128: RefCell<ExpressionEvaluator<Complex<F<f128>>>>,
}

pub trait GenericEvaluate<T> {
    fn evaluate(&self, params: &[T], out: &mut [T]);
    fn evaluate_single(&self, params: &[T]) -> T;
}

// impl<T: Real + Default> GenericEvaluate<T> for GenericEvaluator
// where
//     T: for<'a> From<&'a symbolica::domains::float::Complex<Rational>>,
//     symbolica::domains::float::Complex<Rational>: for<'a> From<&'a T>,
// {
//     fn evaluate(&self, params: &[T], out: &mut [T]) {
//         let rational = self.rational.take().unwrap();
//         let mut t_eval = rational.map_coeff(&|t| T::from(t));
//         t_eval.evaluate(params, out);
//         self.rational.replace(Some(t_eval.map_coeff(&|t| t.into())));
//     }

//     fn evaluate_single(&self, params: &[T]) -> T {
//         let rational = self.rational.take().unwrap();
//         let mut t_eval = rational.map_coeff(&|t| T::from(t));
//         let out = t_eval.evaluate_single(params);
//         self.rational.replace(Some(t_eval.map_coeff(&|t| t.into())));
//         out
//     }
// }

impl GenericEvaluate<Complex<F<f64>>> for GenericEvaluator {
    fn evaluate(&self, params: &[Complex<F<f64>>], out: &mut [Complex<F<f64>>]) {
        if let Some(f64_compiled) = &self.f64_compiled {
            f64_compiled.borrow_mut().evaluate(params, out);
        } else {
            self.f64_eager.borrow_mut().evaluate(params, out);
        }
    }

    fn evaluate_single(&self, params: &[Complex<F<f64>>]) -> Complex<F<f64>> {
        if let Some(f64_compiled) = &self.f64_compiled {
            let mut out = [Complex::default()];
            f64_compiled.borrow_mut().evaluate(params, &mut out);
            out[0]
        } else {
            self.f64_eager.borrow_mut().evaluate_single(params)
        }
    }
}

impl GenericEvaluate<Complex<F<f128>>> for GenericEvaluator {
    fn evaluate(&self, params: &[Complex<F<f128>>], out: &mut [Complex<F<f128>>]) {
        self.f128.borrow_mut().evaluate(params, out);
    }

    fn evaluate_single(&self, params: &[Complex<F<f128>>]) -> Complex<F<f128>> {
        self.f128.borrow_mut().evaluate_single(params)
    }
}

impl GenericEvaluator {
    pub(crate) fn compile(
        &mut self,
        filename: impl AsRef<str>,
        function_name: impl AsRef<str>,
        lib_name: impl AsRef<str>,
        inline_asm: InlineASM,
    ) {
        let compile = self
            .f64_eager
            .borrow()
            .export_cpp(filename.as_ref(), function_name.as_ref(), true, inline_asm)
            .unwrap()
            .compile(lib_name.as_ref(), CompileOptions::default())
            .unwrap()
            .load()
            .unwrap();

        self.f64_compiled = Some(RefCell::new(SerializableCompiledEvaluator {
            evaluator: compile,
            library_filename: lib_name.as_ref().to_string(),
            function_name: function_name.as_ref().to_string(),
        }));
    }

    pub(crate) fn new_from_builder(
        atom: impl AtomCore,
        builder: ParamBuilder<f64>,
        optimization_settings: OptimizationSettings,
    ) -> Self {
        let params: Vec<Atom> = (&builder)
            .into_iter()
            .flat_map(|p| p.params.clone())
            .collect();

        Self::new(atom, &builder.fn_map, &params, optimization_settings)
    }

    fn new(
        atom: impl AtomCore,
        fn_map: &FunctionMap,
        params: &[Atom],
        optimization_settings: OptimizationSettings,
    ) -> Self {
        let tree = atom
            .evaluator(&fn_map, &params, optimization_settings)
            .unwrap();

        let rational = tree.clone();
        let f64_eager = tree
            .clone()
            .map_coeff(&|r| Complex::new(F::from(&r.re), F::from(&r.im)));
        let f128 = tree.map_coeff(&|r| Complex::new(F::from(&r.re), F::from(&r.im)));

        let evaluator = GenericEvaluator {
            expr: atom.as_atom_view().to_owned(),
            rational: RefCell::new(Some(rational)),
            f64_compiled: None,
            f64_eager: RefCell::new(f64_eager),
            f128: RefCell::new(f128),
        };

        evaluator
    }
}

pub trait GenericEvaluatorFloat<T: FloatLike = Self> {
    fn get_evaluator(
        generic_evaluator: &GenericEvaluator,
    ) -> impl Fn(&[Complex<F<T>>]) -> Complex<F<T>>;

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        graph: &'a Graph,
        sample: &'a MomentumSample<T>,
        helicities: &[Helicity],
    ) -> Cow<'a, Vec<Complex<F<T>>>>;
}

impl GenericEvaluatorFloat for f64 {
    #[inline(always)]
    fn get_evaluator(
        generic_evaluator: &GenericEvaluator,
    ) -> impl Fn(&[Complex<F<f64>>]) -> Complex<F<f64>> {
        #[inline(always)]
        |params: &[Complex<F<f64>>]| {
            if let Some(compiled) = &generic_evaluator.f64_compiled {
                let mut out = [Complex::default()];
                compiled.borrow_mut().evaluate(params, &mut out);
                out[0]
            } else {
                generic_evaluator
                    .f64_eager
                    .borrow_mut()
                    .evaluate_single(params)
            }
        }
    }

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        graph: &'a Graph,
        sample: &'a MomentumSample<Self>,
        helicities: &[Helicity],
    ) -> Cow<'a, Vec<Complex<F<Self>>>> {
        param_builder.update_emr_and_get_params(sample, graph, helicities)
    }

    // fn get_debug_evaluator(
    //     generic_evaluator: &GenericEvaluatorDebug,
    // ) -> impl Fn(&[Complex<F<Self>>]) -> Complex<F<Self>> {
    //     #[inline(always)]
    //     |params: &[Complex<F<f64>>]| {
    //         // generic_evaluator
    //         //     .builder
    //         //     .borrow_mut()
    //         //     .fill_in_values(Vec::from_iter(params.iter().cloned()));

    //         // let a = generic_evaluator
    //         //     .builder
    //         //     .borrow()
    //         //     .replace(&generic_evaluator.expr);

    //         // debug!("Replaced atom:{:+>}", a);
    //         // generic_evaluator
    //         //     .expr
    //         //     .evaluate(
    //         //         |c| Complex::new_re(F::<f64>::from(c)),
    //         //         const_map,
    //         //         function_map,
    //         //     )
    //         //     .unwrap()

    //         // generic_evaluator.expr.evaluate(coeff_map, const_map, function_map)

    //         if let Some(compiled) = &generic_evaluator.f64_compiled {
    //             let mut out = [Complex::default()];
    //             compiled.borrow_mut().evaluate(params, &mut out);
    //             out[0]
    //         } else {
    //             generic_evaluator
    //                 .f64_eager
    //                 .borrow_mut()
    //                 .evaluate_single(params)
    //         }
    //     }
    // }
}

impl GenericEvaluatorFloat for f128 {
    #[inline(always)]
    fn get_evaluator(
        generic_evaluator: &GenericEvaluator,
    ) -> impl Fn(&[Complex<F<f128>>]) -> Complex<F<f128>> {
        #[inline(always)]
        |params: &[Complex<F<f128>>]| generic_evaluator.f128.borrow_mut().evaluate_single(params)
    }

    fn get_parameters<'a>(
        param_builder: &'a mut ParamBuilder,
        graph: &'a Graph,
        sample: &'a MomentumSample<Self>,
        helicities: &[Helicity],
    ) -> Cow<'a, Vec<Complex<F<Self>>>> {
        param_builder.update_emr_and_get_params(sample, graph, helicities)
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
            .loop_edges
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
            jacobian: momentum_sample.jacobian.clone(),
            orientation: momentum_sample.orientation,
        }
    }

    pub(crate) fn reinterpret_loop_momenta_and_compute_prefactor_all_channels<T: FloatLike>(
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
    pub(crate) fn reinterpret_loop_momenta_and_compute_prefactor<T: FloatLike>(
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
                    .loop_edges
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

    fn get_terms_mut(&mut self) -> impl Iterator<Item = &mut Self::G>;
    fn get_settings(&self) -> &Settings;
    fn get_graph(&self, graph_id: usize) -> &Self::G;
    fn get_graph_mut(&mut self, graph_id: usize) -> &mut Self::G;
    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor;

    // fn get_builder_cache(&self) -> &ParamBuilder<f64>;
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
        &mut self,
        sample: &MomentumSample<T>,
        settings: &Settings,
    ) -> Complex<F<T>>;

    fn get_multi_channeling_setup(&self) -> &LmbMultiChannelingSetup;
    fn get_graph(&self) -> &Graph;
    fn get_lmbs(&self) -> &TiVec<LmbIndex, LoopMomentumBasis>;
    fn get_num_orientations(&self) -> usize;
    fn get_tropical_sampler(&self) -> &SampleGenerator<3>;
}

fn evaluate_all_rotations<T: FloatLike, I: GammaloopIntegrand>(
    integrand: &mut I,
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
    integrand: &mut I,
    gammaloop_sample: &GammaLoopSample<T>,
) -> Complex<F<f64>> {
    // this is the earliest point where we can set the external momenta
    // param_builder.external_energies_value(&gammaloop_sample.get_default_sample());
    // param_builder.external_spatial_value(&gammaloop_sample.get_default_sample());

    let settings = integrand.get_settings().clone();
    let result = match &gammaloop_sample {
        GammaLoopSample::Default(sample) => integrand
            .get_terms_mut()
            .map(|term: &mut I::G| term.evaluate(sample, &settings))
            .fold(
                Complex::new_re(gammaloop_sample.get_default_sample().zero()),
                |sum, term| sum + term,
            ),
        GammaLoopSample::MultiChanneling { alpha, sample } => integrand
            .get_terms_mut()
            .map(|term: &mut I::G| {
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
                            * term.evaluate(&reparameterized_sample, &settings)
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
            let graph_term = integrand.get_graph_mut(*graph_id);
            match sample {
                DiscreteGraphSample::Default(sample) => graph_term.evaluate(sample, &settings),
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
                        * graph_term.evaluate(&reparameterized_sample, &settings)
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
                                * graph_term.evaluate(&reparameterized_sample, &settings)
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

                    Complex::new_re(prefactor) * graph_term.evaluate(sample, &settings)
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
    integrand: &mut I,
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
                    // let mut param_builder = integrand.get_builder_cache().clone();

                    // param_builder.m_uv_value(Complex::new_re(HARD_CODED_M_UV));
                    // param_builder.m_uv_atom(Atom::var(GS.m_uv));
                    // param_builder.mu_r_sq_value(Complex::new_re(HARD_CODED_M_R_SQ));
                    // param_builder.mu_r_sq_atom(Atom::var(GS.mu_r_sq));
                    // param_builder.model_parameters_atom(model);
                    // param_builder.model_parameters.values =
                    //     integrand.get_model_parameter_cache::<f64>();

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
                    // let mut param_builder = ParamBuilder::<f128>::new();

                    // param_builder.m_uv_value(Complex::new_re(F::from_ff64(HARD_CODED_M_UV)));
                    // param_builder.mu_r_sq_value(Complex::new_re(F::from_ff64(HARD_CODED_M_R_SQ)));
                    // // param_builder.model_parameters.values =
                    //     integrand.get_model_parameter_cache::<f128>();

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
        } else {
            if integrand.get_settings().general.debug > 0 {
                println!("unstable at level: {}", stability_level.precision);
                if let Ok(gammaloop_sample) = parameterize::<f64, I>(sample, integrand) {
                    let rotated_samples: Vec<_> = integrand
                        .get_rotations()
                        .map(|rotation| gammaloop_sample.get_rotated_sample(rotation))
                        .collect();

                    for (sample, result) in rotated_samples.iter().zip(&results) {
                        let default_sample = sample.get_default_sample();
                        println!(
                            "loop_moms: {}, external_moms: {}",
                            format!("{}", default_sample.loop_moms()).blue(),
                            format!("{:?}", default_sample.external_moms()).blue()
                        );

                        println!(
                            "result of current level: {}",
                            format!("{:16e}", result).blue()
                        );
                    }
                } else {
                    println!("parameterization failed");
                }
            }
        }
    }

    if integrand.get_settings().general.debug > 0 {
        println!("result at each level:");
        for level_result in results_of_stability_levels.iter() {
            println!(
                "level: {}. result: {}",
                format!("{}", level_result.stability_level_used).green(),
                format!("{:16e}", level_result.result).blue()
            );
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

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ParamValuePairs<T: FloatLike> {
    values: Vec<Complex<F<T>>>,
    params: Vec<Atom>,
}

impl<T: FloatLike> ParamValuePairs<T> {
    pub fn validate(&self) {
        assert_eq!(
            self.values.len(),
            self.params.len(),
            "Number of values and parameters must match"
        );
    }

    pub fn map_values<U: FloatLike>(
        &self,
        mut map: impl FnMut(&Complex<F<T>>) -> Complex<F<U>>,
    ) -> ParamValuePairs<U> {
        let values = self.values.iter().map(&mut map).collect();
        ParamValuePairs {
            values,
            params: self.params.clone(),
        }
    }
    pub fn replacement(&self) -> Vec<Replacement> {
        let mut replacements = Vec::new();
        for (p, v) in self.params.iter().zip_eq(self.values.iter()) {
            println!("{p}");
            let replacement = Replacement::new(
                p.clone().to_pattern(),
                Atom::num(v.clone().to_coefficient()),
            );
            replacements.push(replacement);
        }
        replacements
    }

    pub fn extract_and_fill(&mut self, values: &mut Vec<Complex<F<T>>>) {
        self.values = values.split_off(self.params.len())
    }
}

impl<T: FloatLike> Default for ParamValuePairs<T> {
    fn default() -> Self {
        Self {
            values: Vec::new(),
            params: Vec::new(),
        }
    }
}

// impl<T:FloatLike> ParamValuePairs<T>{
//     pub fn map_values<U:FloatLike>(&self,map:impl FnMut(&Complex<F<T>>)->Complex<F<U>>)->ParamValuePairs<U>{

//     }
// }

impl<T: FloatLike> ParamValuePairs<T>
where
    T::Higher: FloatLike,
    T::Lower: FloatLike,
{
    fn higher(&self) -> ParamValuePairs<T::Higher> {
        ParamValuePairs {
            values: self
                .values
                .iter()
                .map(|v| v.map_ref(|v| v.higher()))
                .collect(),
            params: self.params.clone(),
        }
    }
    fn lower(&self) -> ParamValuePairs<T::Lower> {
        ParamValuePairs {
            values: self
                .values
                .iter()
                .map(|v| v.map_ref(|v| v.lower()))
                .collect(),
            params: self.params.clone(),
        }
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ParamBuilder<T: FloatLike = f64> {
    // values: Vec<Complex<F<T>>
    m_uv: ParamValuePairs<T>,
    mu_r_sq: ParamValuePairs<T>,
    pub model_parameters: ParamValuePairs<T>,
    pub external_energies: ParamValuePairs<T>,
    pub external_spatial: ParamValuePairs<T>,
    polarizations: ParamValuePairs<T>,
    emr_spatial: ParamValuePairs<T>,
    tstar: ParamValuePairs<T>,
    h_function: ParamValuePairs<T>,
    derivative_at_tstar: ParamValuePairs<T>,
    uv_damp: ParamValuePairs<T>,
    radius: ParamValuePairs<T>,
    radius_star: ParamValuePairs<T>,
    pub reps: Vec<(Atom, Atom)>,
    // pub eager_const_map: HashMap<Atom, Complex<F<T>>>,
    // pub eager_function_map: HashMap<Symbol, EvaluationFn<Atom, Complex<F<T>>>>,
    // pub eager_fn_map:
    pub fn_map: FunctionMap,
}

pub trait UpdateAndGetParams<T: FloatLike> {
    fn update_emr_and_get_params<'a>(
        &'a mut self,
        sample: &'a MomentumSample<T>,
        graph: &'a Graph,
        helicities: &[Helicity],
    ) -> Cow<'a, Vec<Complex<F<T>>>>;
}

impl UpdateAndGetParams<f64> for ParamBuilder<f64> {
    fn update_emr_and_get_params<'a>(
        &'a mut self,
        sample: &'a MomentumSample<f64>,
        graph: &'a Graph,
        helicities: &[Helicity],
    ) -> Cow<'a, Vec<Complex<F<f64>>>> {
        let emr_spatial: Vec<_> = graph
            .iter_edges()
            .flat_map(|(pair, edge_id, _)| {
                if let HedgePair::Paired { .. } = pair {
                    let emr_vec = graph.loop_momentum_basis.edge_signatures[edge_id]
                        .compute_three_momentum_from_four(
                            sample.loop_moms(),
                            sample.external_moms(),
                        );
                    // println!("{edge_id}:{emr_vec}");
                    vec![
                        Complex::new_re(emr_vec.px),
                        Complex::new_re(emr_vec.py),
                        Complex::new_re(emr_vec.pz),
                    ]
                } else {
                    vec![]
                }
            })
            .collect();

        // parse!("s").evaluator(fn_map, params, optimization_settings).unwrap().

        self.external_spatial_value(sample);
        self.external_energies_value(sample);
        self.emr_spatial.values = emr_spatial;
        self.polarizations_values(graph, sample.external_moms(), helicities);

        println!("ParamBuilder after eval f64:\n{}", self);

        Cow::Owned(self.clone().build_values()) // ideally borrows a single vec
    }
}

impl UpdateAndGetParams<f128> for ParamBuilder<f64> {
    fn update_emr_and_get_params(
        &mut self,
        sample: &MomentumSample<f128>,
        graph: &Graph,
        helicities: &[Helicity],
    ) -> Cow<Vec<Complex<F<f128>>>> {
        let emr_spatial: Vec<_> = graph
            .iter_edges()
            .flat_map(|(pair, edge_id, _)| {
                if let HedgePair::Paired { .. } = pair {
                    let emr_vec = graph.loop_momentum_basis.edge_signatures[edge_id]
                        .compute_three_momentum_from_four(
                            sample.loop_moms(),
                            sample.external_moms(),
                        );
                    vec![
                        Complex::new_re(emr_vec.px),
                        Complex::new_re(emr_vec.py),
                        Complex::new_re(emr_vec.pz),
                    ]
                } else {
                    vec![]
                }
            })
            .collect();

        let mut new_param_builder = self.higher();
        new_param_builder.emr_spatial.values = emr_spatial;
        new_param_builder.external_spatial_value(sample);
        new_param_builder.external_energies_value(sample);
        new_param_builder.polarizations_values(graph, sample.external_moms(), helicities);

        // println!("ParamBuilder before eval f128:\n{}", self);

        Cow::Owned(new_param_builder.build_values())
    }
}

impl<T: FloatLike> ParamBuilder<T>
where
    T::Higher: FloatLike,
    T::Lower: FloatLike,
{
    fn higher(&self) -> ParamBuilder<T::Higher> {
        ParamBuilder {
            fn_map: self.fn_map.clone(),
            reps: self.reps.clone(),
            m_uv: self.m_uv.higher(),
            mu_r_sq: self.mu_r_sq.higher(),
            model_parameters: self.model_parameters.higher(),
            external_energies: self.external_energies.higher(),
            external_spatial: self.external_spatial.higher(),
            polarizations: self.polarizations.higher(),
            emr_spatial: self.emr_spatial.higher(),
            tstar: self.tstar.higher(),
            h_function: self.h_function.higher(),
            derivative_at_tstar: self.derivative_at_tstar.higher(),
            uv_damp: self.uv_damp.higher(),
            radius: self.radius.higher(),
            radius_star: self.radius_star.higher(),
        }
    }
    fn lower(&self) -> ParamBuilder<T::Lower> {
        ParamBuilder {
            fn_map: self.fn_map.clone(),
            reps: self.reps.clone(),
            m_uv: self.m_uv.lower(),
            mu_r_sq: self.mu_r_sq.lower(),
            model_parameters: self.model_parameters.lower(),
            external_energies: self.external_energies.lower(),
            external_spatial: self.external_spatial.lower(),
            polarizations: self.polarizations.lower(),
            emr_spatial: self.emr_spatial.lower(),
            tstar: self.tstar.lower(),
            h_function: self.h_function.lower(),
            derivative_at_tstar: self.derivative_at_tstar.lower(),
            uv_damp: self.uv_damp.lower(),
            radius: self.radius.lower(),
            radius_star: self.radius_star.lower(),
        }
    }
}
impl<T: FloatLike> ParamBuilder<T> {
    pub fn validate(&self) {
        debug!("Validating mu_r_sq");
        self.mu_r_sq.validate();
        debug!("Validating model_parameters");
        self.model_parameters.validate();
        debug!("Validating external_energies");
        self.external_energies.validate();
        debug!("Validating external_spatial");
        self.external_spatial.validate();
        debug!("Validating polarizations");
        self.polarizations.validate();
        debug!("Validating emr_spatial");
        self.emr_spatial.validate();
        debug!("Validating tstar");
        self.tstar.validate();
        debug!("Validating h_function");
        self.h_function.validate();
        debug!("Validating derivative_at_tstar");
        self.derivative_at_tstar.validate();
        debug!("Validating uv_damp");
        self.uv_damp.validate();
        debug!("Validating radius");
        self.radius.validate();
        debug!("Validating radius_star");
        self.radius_star.validate();
    }

    pub fn fill_in_values(&mut self, mut values: Vec<Complex<F<T>>>) {
        self.m_uv.extract_and_fill(&mut values);
        self.mu_r_sq.extract_and_fill(&mut values);
        self.model_parameters.extract_and_fill(&mut values);
        self.external_energies.extract_and_fill(&mut values);
        self.external_spatial.extract_and_fill(&mut values);
        self.polarizations.extract_and_fill(&mut values);
        self.emr_spatial.extract_and_fill(&mut values);
        self.tstar.extract_and_fill(&mut values);
        self.h_function.extract_and_fill(&mut values);
        self.derivative_at_tstar.extract_and_fill(&mut values);
        self.uv_damp.extract_and_fill(&mut values);
        self.radius.extract_and_fill(&mut values);
        self.radius_star.extract_and_fill(&mut values);
    }

    pub fn replace_non_emr(&self, atom: impl AtomCore) -> Atom {
        let reps = self
            .reps
            .iter()
            .map(|(a, b)| Replacement::new(a.clone().to_pattern(), b.clone()))
            .collect_vec();
        let mut a = atom.replace_multiple(&reps);

        for r in [
            &self.m_uv,
            &self.mu_r_sq,
            &self.model_parameters,
            &self.external_energies,
            &self.external_spatial,
            // &self.polarizations,
            // &self.emr_spatial,
            &self.tstar,
            &self.h_function,
            &self.derivative_at_tstar,
            &self.uv_damp,
            &self.radius,
            &self.radius_star,
        ]
        .into_iter()
        {
            a = a.replace_multiple(&r.replacement());
        }
        a
    }

    pub fn add_tagged_function(
        &mut self,
        name: Symbol,
        tags: Vec<Atom>,
        rename: String,
        args: Vec<Symbol>,
        body: Atom,
    ) -> Result<(), &str> {
        self.reps.push((
            FunctionBuilder::new(name)
                .add_args(&tags)
                .add_args(&args)
                .finish(),
            body.clone(),
        ));

        self.fn_map
            .add_tagged_function(name, tags, rename, args, body)

        // body.evaluate(coeff_map, const_map, function_map)
    }

    pub fn add_function(
        &mut self,
        name: Symbol,
        rename: String,
        args: Vec<Symbol>,
        body: Atom,
    ) -> Result<(), &str> {
        self.fn_map.add_function(name, rename, args, body)
    }

    pub fn add_constant(&mut self, key: Atom, value: symbolica::domains::float::Complex<Rational>) {
        self.fn_map.add_constant(key, value)
    }

    pub(crate) fn build_values(self) -> Vec<Complex<F<T>>> {
        let mut values = Vec::with_capacity(100);
        for value in self.into_iter() {
            values.extend(value.values);
        }
        values
    }
    pub(crate) fn build_params(self) -> Vec<Atom> {
        let mut values = Vec::with_capacity(100);
        for value in self.into_iter() {
            values.extend(value.params);
        }
        values
    }

    pub(crate) fn new() -> Self {
        Self {
            fn_map: FunctionMap::default(),
            m_uv: ParamValuePairs::default(),
            mu_r_sq: ParamValuePairs::default(),
            model_parameters: ParamValuePairs::default(),
            external_energies: ParamValuePairs::default(),
            polarizations: ParamValuePairs::default(),
            external_spatial: ParamValuePairs::default(),
            emr_spatial: ParamValuePairs::default(),
            tstar: ParamValuePairs::default(),
            h_function: ParamValuePairs::default(),
            derivative_at_tstar: ParamValuePairs::default(),
            uv_damp: ParamValuePairs::default(),
            radius: ParamValuePairs::default(),
            radius_star: ParamValuePairs::default(),
            reps: Vec::new(),
        }
    }

    pub(crate) fn m_uv_atom(&mut self, m_uv: Atom) {
        self.m_uv.params = vec![m_uv];
    }

    pub(crate) fn m_uv_value(&mut self, m_uv: Complex<F<T>>) {
        self.m_uv.values = vec![m_uv];
    }

    pub(crate) fn mu_r_sq_atom(&mut self, mu_r_sq: Atom) {
        self.mu_r_sq.params = vec![mu_r_sq];
    }

    pub(crate) fn mu_r_sq_value(&mut self, mu_r_sq: Complex<F<T>>) {
        self.mu_r_sq.values = vec![mu_r_sq];
    }

    pub(crate) fn model_parameters_atom(&mut self, model: &Model) {
        self.model_parameters.params = model.generate_params();
    }

    pub(crate) fn model_parameters_value(&mut self, model: &Model) {
        self.model_parameters.values = model.generate_values();
    }

    pub(crate) fn external_energies_atom(&mut self, graph: &Graph) {
        self.external_energies.params = graph.get_external_energy_atoms();
    }

    pub(crate) fn add_external_four_mom(&mut self, ext: &ExternalFourMomenta<F<T>>) {
        self.external_energies.values = ext
            .iter()
            .map(|x| Complex::new_re(x.temporal.value.clone()))
            .collect_vec();
        self.external_spatial.values = ext
            .iter()
            .flat_map(|x| x.spatial.clone().into_iter().map(|c| Complex::new_re(c)))
            .collect_vec();
    }

    pub(crate) fn polarization_params(&mut self, graph: &Graph) {
        let mut pols = graph.global_prefactor.polarizations();
        pols.sort_by(|a, b| a.0.cmp(&b.0));
        let mut params = Vec::new();

        for (p, a) in pols {
            match ParsingNet::try_from_view(a.as_view(), TENSORLIB.deref())
                .unwrap()
                .result_tensor(TENSORLIB.deref())
                .unwrap()
            {
                ExecutionResult::One => {}
                ExecutionResult::Zero => {}
                ExecutionResult::Val(a) => {
                    for (_, val) in a.iter_flat() {
                        let AtomViewOrConcrete::Atom(a) = val else {
                            panic!("SHOULD BE ATOMVIEW")
                        };

                        params.push(a.to_owned());
                    }
                }
            }
        }

        self.polarizations.params = params;

        // self.polarizations.params = graph.generate_polarization_params();
    }

    pub(crate) fn polarizations_values(
        &mut self,
        graph: &Graph,
        ext: &ExternalFourMomenta<F<T>>,
        helicities: &[Helicity],
    ) {
        let mut pols = graph.global_prefactor.polarizations();
        pols.sort_by(|a, b| a.0.cmp(&b.0));

        let mut vals = Vec::new();

        for (p, a) in pols {
            let extid = graph.loop_momentum_basis.ext_from(p.eid).unwrap();
            let hel = p.hel.unwrap_or(helicities[extid.0]);
            // println!("MOM:{}", ext[extid]);
            // println!("Pol {},{},{:?}", extid, hel, p.pol_type);
            let pol = match p.pol_type {
                PolType::Epsilon => ext[extid].pol(hel),
                PolType::EpsilonBar => ext[extid].pol(hel).bar(),
                PolType::Scalar => {
                    continue;
                }
                PolType::U => ext[extid].u(hel.try_into().unwrap()),
                PolType::V => ext[extid].v(hel.try_into().unwrap()),
                PolType::UBar => ext[extid].u(hel.try_into().unwrap()).bar(),
                PolType::VBar => ext[extid].v(hel.try_into().unwrap()).bar(),
            };

            for (_, val) in pol.tensor.iter_flat() {
                // println!("{val}");
                vals.push(val.clone());
            }
        }

        self.polarizations.values = vals;

        // self.polarizations.params = graph.generate_polarization_params();
    }

    pub(crate) fn external_energies_value(&mut self, momentum_sample: &MomentumSample<T>) {
        self.external_energies.values = momentum_sample
            .external_moms()
            .iter()
            .map(|x| Complex::new_re(x.temporal.value.clone()))
            .collect_vec();
    }

    pub(crate) fn external_spatial_atom(&mut self, graph: &Graph) {
        self.external_spatial.params = graph
            .iter_edges()
            .flat_map(|(pair, edge_id, _)| {
                if let HedgePair::Unpaired { .. } = pair {
                    vec![
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([1]))),
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([2]))),
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([3]))),
                    ]
                } else {
                    vec![]
                }
            })
            .collect();
    }

    pub(crate) fn external_spatial_value(&mut self, momentum_sample: &MomentumSample<T>) {
        self.external_spatial.values = momentum_sample
            .external_moms()
            .iter()
            .flat_map(|x| x.spatial.clone().into_iter().map(|c| Complex::new_re(c)))
            .collect_vec();
    }

    pub(crate) fn emr_spatial_atom(&mut self, graph: &Graph) {
        self.emr_spatial.params = graph
            .iter_edges()
            .flat_map(|(pair, edge_id, _)| {
                if let HedgePair::Paired { .. } = pair {
                    vec![
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([1]))),
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([2]))),
                        GS.emr_mom(edge_id, Atom::from(ExpandedIndex::from_iter([3]))),
                    ]
                } else {
                    vec![]
                }
            })
            .collect();
    }

    pub(crate) fn emr_spatial_value(&mut self, emr_spatial: Vec<Complex<F<T>>>) {
        self.emr_spatial.values = emr_spatial;
    }

    pub(crate) fn tstar_atom(&mut self, tstar: Atom) {
        self.tstar.params = vec![tstar];
    }

    pub(crate) fn tstar_value(&mut self, tstar: Complex<F<T>>) {
        self.tstar.values = vec![tstar];
    }
    pub(crate) fn h_function_atom(&mut self, h_function: Atom) {
        self.h_function.params = vec![h_function];
    }
    pub(crate) fn h_function_value(&mut self, h_function: Complex<F<T>>) {
        self.h_function.values = vec![h_function];
    }
    pub(crate) fn derivative_at_tstar_atom(&mut self, derivative_at_tstar: Atom) {
        self.derivative_at_tstar.params = vec![derivative_at_tstar];
    }
    pub(crate) fn derivative_at_tstar_value(&mut self, derivative_at_tstar: Complex<F<T>>) {
        self.derivative_at_tstar.values = vec![derivative_at_tstar];
    }

    pub(crate) fn uv_damp_atom(&mut self, uv_dampers: Vec<Atom>) {
        self.uv_damp.params = uv_dampers;
    }
    pub(crate) fn uv_damp_value(&mut self, uv_dampers: Vec<Complex<F<T>>>) {
        self.uv_damp.values = uv_dampers;
    }

    pub(crate) fn radius_atom(&mut self, radius: Atom) {
        self.radius.params = vec![radius];
    }
    pub(crate) fn radius_value(&mut self, radius: Complex<F<T>>) {
        self.radius.values = vec![radius];
    }

    pub(crate) fn radius_star_atom(&mut self, radius_star: Atom) {
        self.radius_star.params = vec![radius_star];
    }
    pub(crate) fn radius_star_value(&mut self, radius_star: Complex<F<T>>) {
        self.radius_star.values = vec![radius_star];
    }
}

impl<T: FloatLike> Display for ParamBuilder<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut table = tabled::builder::Builder::new();

        for (lhs, rhs) in &self.reps {
            table.push_record(vec![lhs.to_string(), rhs.to_string()]);
        }

        for i in self {
            for (v, p) in i.values.iter().zip(i.params.iter()) {
                table.push_record(vec![p.to_string(), v.to_string()]);
            }
        }

        table.build().with(Style::rounded()).to_string().fmt(f)
    }
}

impl<T: FloatLike> IntoIterator for ParamBuilder<T> {
    type Item = ParamValuePairs<T>;
    type IntoIter = std::array::IntoIter<Self::Item, 13>;

    fn into_iter(self) -> Self::IntoIter {
        [
            self.m_uv,
            self.mu_r_sq,
            self.model_parameters,
            self.external_energies,
            self.external_spatial,
            self.polarizations,
            self.emr_spatial,
            self.tstar,
            self.h_function,
            self.derivative_at_tstar,
            self.uv_damp,
            self.radius,
            self.radius_star,
        ]
        .into_iter()
    }
}

impl<'a, T: FloatLike> IntoIterator for &'a ParamBuilder<T> {
    type Item = &'a ParamValuePairs<T>;
    type IntoIter = std::array::IntoIter<Self::Item, 13>;

    fn into_iter(self) -> Self::IntoIter {
        [
            &self.m_uv,
            &self.mu_r_sq,
            &self.model_parameters,
            &self.external_energies,
            &self.external_spatial,
            &self.polarizations,
            &self.emr_spatial,
            &self.tstar,
            &self.h_function,
            &self.derivative_at_tstar,
            &self.uv_damp,
            &self.radius,
            &self.radius_star,
        ]
        .into_iter()
    }
}
