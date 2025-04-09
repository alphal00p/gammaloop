use std::cell::RefCell;

use crate::evaluation_result::EvaluationResult;
use crate::integrands::{HasIntegrand, Integrand};
use crate::integrate::UserData;
use crate::momentum::Energy;
use crate::momentum_sample::{BareMomentumSample, LoopMomenta, MomentumSample};
use crate::new_graph::{FeynmanGraph, Graph, LmbIndex, LoopMomentumBasis};
use crate::utils::{format_for_compare_digits, FloatLike, F};
use colored::Colorize;
use derive_more::{From, Into};
use enum_dispatch::enum_dispatch;
use itertools::Itertools;
use linnet::half_edge::involution::EdgeIndex;
use momtrop::float::MomTropFloat;
use serde::{Deserialize, Serialize};
use spenso::complex::Complex;
use spenso::contraction::IsZero;
use spenso::parametric::SerializableCompiledEvaluator;
use std::time::Duration;
use symbolica::evaluate::ExpressionEvaluator;
use symbolica::numerical_integration::{Grid, Sample};
use typed_index_collections::TiVec;
pub mod amplitude_integrand;
pub mod cross_section_integrand;
pub mod gammaloop_sample;
use crate::observables::EventManager;
use crate::utils::f128;
use crate::{
    IntegratedPhase, IntegratorSettings, Precision, Settings, StabilityLevelSetting,
    StabilitySettings,
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

    // todo provide max wgt as Complex<F<f64>>
    let average_for_comparison = match settings.integrator.integrated_phase {
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

#[derive(Clone)]
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
#[derive(Debug, Clone, Serialize, Deserialize)]
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
        }
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
