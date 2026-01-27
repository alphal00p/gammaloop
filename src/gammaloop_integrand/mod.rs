use std::fs;
use std::path::Path;

use crate::evaluation_result::{EvaluationMetaData, EvaluationResult};
use crate::graph::{FeynmanGraph, Graph, GraphGroup, GroupId, LmbIndex, LoopMomentumBasis};
use crate::integrands::{HasIntegrand, Integrand};
use crate::integrate::UserData;
use crate::model::Model;
use crate::momentum::Rotation;
use crate::momentum_sample::{BareMomentumSample, LoopMomenta, MomentumSample};
use crate::settings::GlobalSettings;
use crate::utils::{F, FloatLike, format_for_compare_digits, get_n_dim_for_n_loop_momenta};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::owo_colors::OwoColorize;
use colored::Colorize;
use derive_more::{From, Into};
use enum_dispatch::enum_dispatch;
use eyre::Context;
use gammaloop_sample::{DiscreteGraphSample, GammaLoopSample, parameterize};
use itertools::Itertools;
use linnet::half_edge::involution::EdgeVec;
use momtrop::SampleGenerator;
use momtrop::float::MomTropFloat;
use serde::{Deserialize, Serialize};
use spenso::algebra::algebraic_traits::IsZero;
use spenso::algebra::complex::Complex;
use std::time::Duration;
use symbolica::numerical_integration::{ContinuousGrid, DiscreteGrid, Grid, Sample};
use tracing::debug;
use typed_index_collections::TiVec;
pub mod amplitude_integrand;
pub mod cache_debugging;
pub mod cross_section_integrand;
pub mod gammaloop_sample;
use crate::observables::EventManager;
use crate::utils::f128;
use crate::{
    DependentMomentaConstructor, GammaLoopContext, settings::RuntimeSettings,
    settings::runtime::DiscreteGraphSamplingSettings, settings::runtime::DiscreteGraphSamplingType,
    settings::runtime::IntegratorSettings, settings::runtime::Precision,
    settings::runtime::SamplingSettings, settings::runtime::StabilityLevelSetting,
    settings::runtime::StabilitySettings,
};
use color_eyre::Result;

pub mod evaluators;
pub use evaluators::{GenericEvaluator, GenericEvaluatorFloat};

pub mod param_builder;
pub use param_builder::{ParamBuilder, ParamValuePairs, ThresholdParams, UpdateAndGetParams};

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
#[enum_dispatch(HasIntegrand)]
pub enum GLIntegrand {
    Amplitude(amplitude_integrand::AmplitudeIntegrand),
    CrossSection(cross_section_integrand::CrossSectionIntegrand),
}

impl GLIntegrand {
    pub fn warm_up(&mut self, model: &Model) -> Result<()> {
        match self {
            Self::Amplitude(a) => a.warm_up(model),
            Self::CrossSection(a) => a.warm_up(model),
        }
    }

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
            GLIntegrand::Amplitude(integrand) => integrand.save(path, override_existing),
            GLIntegrand::CrossSection(integrand) => integrand.save(path, override_existing),
        }
    }

    pub(crate) fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        settings: &GlobalSettings,
        thread_pool: &rayon::ThreadPool,
    ) -> Result<()> {
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
            GLIntegrand::Amplitude(integrand) => {
                integrand.compile(path, override_existing, settings, thread_pool)
            }
            GLIntegrand::CrossSection(integrand) => {
                integrand.compile(path, override_existing, settings, thread_pool)
            }
        }
    }

    pub fn get_settings(&self) -> &RuntimeSettings {
        match self {
            GLIntegrand::Amplitude(integrand) => &integrand.settings,
            GLIntegrand::CrossSection(integrand) => &integrand.settings,
        }
    }

    /// Used to create the use_data_generator closure for havana_integrate
    pub fn user_data_generator(&self, num_cores: usize, _settings: &RuntimeSettings) -> UserData {
        UserData {
            integrand: vec![Integrand::GLIntegrand(self.clone()); num_cores],
        }
    }

    pub fn get_mut_settings(&mut self) -> &mut RuntimeSettings {
        match self {
            GLIntegrand::Amplitude(integrand) => &mut integrand.settings,
            GLIntegrand::CrossSection(integrand) => &mut integrand.settings,
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
                required_precision_for_re: 1e-5,
                required_precision_for_im: 1e-5,
                escalate_for_large_weight_threshold: -1.,
            }]
        }
    } else {
        settings.levels.clone()
    }
}

#[inline]
fn stability_check(
    _settings: &RuntimeSettings,
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
        error.re > F(stability_settings.required_precision_for_re)
            || error.im > F(stability_settings.required_precision_for_im)
    });

    if let Some(unstable_index) = unstable_sample {
        let unstable_point = results[unstable_index];

        let ((real_formatted, rotated_real_formatted), (imag_formatted, rotated_imag_formatted)) = (
            format_for_compare_digits(average.re, unstable_point.re),
            format_for_compare_digits(average.im, unstable_point.im),
        );

        debug!("{}", "\nUnstable point detected:".red());
        debug!("\taverage result: {} + {}i", real_formatted, imag_formatted,);
        debug!(
            "\trotated result: {} + {}i",
            rotated_real_formatted, rotated_imag_formatted,
        );
    }

    let stable = unstable_sample.is_none();

    let below_wgt_threshold =
        if stability_settings.escalate_for_large_weight_threshold > 0. && max_eval.is_non_zero() {
            average.re.abs() * wgt
                < F(stability_settings.escalate_for_large_weight_threshold) * max_eval.re
                || average.im.abs() * wgt
                    < F(stability_settings.escalate_for_large_weight_threshold) * max_eval.im
        } else {
            true
        };

    (average, stable && below_wgt_threshold)
}

#[inline]
fn stability_check_on_norm(
    _settings: &RuntimeSettings,
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
        .fold(F(0.0), |acc, x| acc + x.norm_squared().sqrt())
        / F(results.len() as f64);

    let mut errors = results.iter().map(|res| {
        let res = res.norm_squared().sqrt();

        if IsZero::is_zero(&res) && IsZero::is_zero(&average) {
            F(0.)
        } else {
            ((res - average) / average).abs()
        }
    });

    let unstable_sample =
        errors.position(|error| error > F(stability_settings.required_precision_for_re));

    if let Some(unstable_index) = unstable_sample {
        let unstable_point = results[unstable_index];

        let (real_formatted, rotated_real_formatted) =
            format_for_compare_digits(average, unstable_point.re);

        debug!("{}", "\nUnstable point detected:".red());
        debug!("\tnormed average result: {}", real_formatted,);
        debug!("\tnormed rotated result: {}", rotated_real_formatted,);
    }

    let stable = unstable_sample.is_none();

    let below_wgt_threshold =
        if stability_settings.escalate_for_large_weight_threshold > 0. && max_eval.is_non_zero() {
            average.abs() * wgt
                < F(stability_settings.escalate_for_large_weight_threshold)
                    * max_eval.norm_squared().sqrt()
        } else {
            true
        };

    (results[0], stable && below_wgt_threshold)
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
#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct LmbMultiChannelingSetup {
    pub channels: TiVec<ChannelIndex, LmbIndex>,
    pub graph: Graph,
    pub all_bases: TiVec<LmbIndex, LoopMomentumBasis>,
}

impl LmbMultiChannelingSetup {
    fn reinterpret_loop_momenta_impl<T: FloatLike>(
        &self,
        channel_index: ChannelIndex,
        momentum_sample: &BareMomentumSample<T>,
        loop_mom_cache_id: usize,
    ) -> BareMomentumSample<T> {
        let channel_lmb = &self.all_bases[self.channels[channel_index]];
        let new_loop_moms: LoopMomenta<F<T>> = self
            .graph
            .loop_momentum_basis
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
            dual_loop_moms: momentum_sample.dual_loop_moms.clone().map(|dlm| todo!()),
            loop_mom_cache_id,
            loop_mom_base_cache_id: momentum_sample.loop_mom_base_cache_id,
            external_mom_cache_id: momentum_sample.external_mom_cache_id,
            external_mom_base_cache_id: momentum_sample.external_mom_base_cache_id,
            external_moms: momentum_sample.external_moms.clone(),
            jacobian: momentum_sample.jacobian.clone(),
            orientation: momentum_sample.orientation,
        }
    }

    /// Note this increments the loop_mom_cache_id of all of the returned BareMomentumSample
    #[allow(dead_code)]
    pub(crate) fn reinterpret_loop_momenta_and_compute_prefactor_all_channels<T: FloatLike>(
        &self,
        momentum_sample: &MomentumSample<T>,
        model: &Model,
        alpha: &F<T>,
        cache: bool,
    ) -> TiVec<ChannelIndex, (MomentumSample<T>, F<T>)> {
        let mut loop_mom_cache_id = momentum_sample.sample.loop_mom_cache_id;
        self.channels
            .iter_enumerated()
            .map(|(channel_index, _)| {
                if cache {
                    loop_mom_cache_id += 1;
                }
                self.reinterpret_loop_momenta_and_compute_prefactor(
                    channel_index,
                    momentum_sample,
                    loop_mom_cache_id,
                    model,
                    alpha,
                )
            })
            .collect()
    }

    /// This function is used to do do LMB multi-channeling without fully switching to a different lmb
    /// for each channel. The momenta provided are reinterpreted as loop momenta of the lmb corresponding to the channel_index.
    /// Then we transform these loop momenta to the fixed lmb of the graph. The prefactor is immediately computed for the requested channel
    ///
    /// Note this increments the loop_mom_cache_id of the returned BareMomentumSample
    pub(crate) fn reinterpret_loop_momenta_and_compute_prefactor<T: FloatLike>(
        &self,
        channel_index: ChannelIndex,
        momentum_sample: &MomentumSample<T>,
        loop_mom_cache_id: usize,
        model: &Model,
        alpha: &F<T>,
    ) -> (MomentumSample<T>, F<T>) {
        let sample = MomentumSample {
            sample: self.reinterpret_loop_momenta_impl(
                channel_index,
                &momentum_sample.sample,
                loop_mom_cache_id,
            ), // uuid: momentum_sample.uuid,
        };

        let prefactor = self.compute_prefactor_impl(channel_index, &sample, model, alpha);

        (sample, prefactor)
    }

    /// Computes the prefactor for the given channel index and momentum sample.
    fn compute_prefactor_impl<T: FloatLike>(
        &self,
        channel_index: ChannelIndex,
        momentum_sample: &MomentumSample<T>,
        model: &Model,
        alpha: &F<T>,
    ) -> F<T> {
        let all_energies = self.graph.get_energy_cache(
            model,
            &momentum_sample.sample.loop_moms,
            &momentum_sample.sample.external_moms,
            &self.graph.loop_momentum_basis,
        );

        let mut numerator = momentum_sample.zero();

        let denominators = self
            .channels
            .iter()
            .map(|&lmb_index| {
                let channel_product = self.all_bases[lmb_index]
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

    fn warm_up(&mut self, model: &Model) -> Result<()>;
    fn get_rotations(&self) -> impl Iterator<Item = &Rotation>;

    fn increment_loop_cache_id(&mut self, val: usize);
    fn loop_cache_id(&self) -> usize;

    fn increment_external_cache_id(&mut self, val: usize);
    fn external_cache_id(&self) -> usize;

    /// Signal that external momenta configuration has actually changed
    /// This increments the external cache ID to invalidate cached computations
    ///
    /// # Usage
    /// Call this method when:
    /// - You change the external momenta values in your Monte Carlo sampling
    /// - You switch to a different kinematic configuration
    /// - You want to force recomputation of cached polarizations/external quantities
    ///
    /// # Example
    /// ```rust,ignore
    /// // When you change external momenta in your sampling
    /// integrand.signal_external_momenta_changed();
    /// let new_sample = parameterize(&sample_point, &mut integrand)?;
    ///
    /// // Rotations will automatically get new cache IDs
    /// let rotated_samples = evaluate_all_rotations(&new_sample, &mut integrand, true)?;
    ///
    /// // Next iteration - revert to base configuration to reuse cache
    /// integrand.revert_to_base_external_cache_id();
    /// let next_sample = parameterize(&next_sample_point, &mut integrand)?;
    /// ```
    fn signal_external_momenta_changed(&mut self) {
        self.increment_external_cache_id(1);
    }

    /// Get the current external cache ID for reuse (doesn't increment)
    ///
    /// This returns the current cache ID without incrementing it, allowing
    /// new samples to reuse cached computations for the same external
    /// momenta configuration.
    fn get_current_external_cache_id(&self) -> usize {
        self.external_cache_id()
    }

    /// Revert to the base external cache ID for the current configuration
    ///
    /// This allows new samples to reuse the base cache ID when they represent
    /// the same underlying external momenta configuration (e.g., after rotations)
    ///
    /// # Usage
    /// Call this method when you want to create a new sample that represents
    /// the same base external momenta configuration as before, allowing reuse
    /// of cached polarizations and other external-dependent computations.
    ///
    /// # Example Cache Flow
    /// ```text
    /// 1. Initial sample:           cache_id = 0, base_cache_id = 0
    /// 2. Rotation 1:               cache_id = 1, base_cache_id = 0
    /// 3. Rotation 2:               cache_id = 2, base_cache_id = 0
    /// 4. revert_to_base():         cache_id = 0, base_cache_id = 0  // Reuse cache!
    /// 5. signal_changed():         cache_id = 3, base_cache_id = 3  // New base config
    /// 6. Rotation of new config:   cache_id = 4, base_cache_id = 3
    /// 7. revert_to_base():         cache_id = 3, base_cache_id = 3  // Reuse new base
    /// ```
    fn revert_to_base_external_cache_id(&mut self);

    /// Check if external momenta caching is beneficial
    ///
    /// Returns true if the same external cache ID has been used multiple times,
    /// indicating that caching is providing benefits.
    fn is_external_caching_beneficial(&self) -> bool {
        // Simple heuristic: if we've created samples without incrementing
        // cache ID, then caching is being used
        self.external_cache_id() < self.loop_cache_id()
    }

    /// Force a cache consistency check with detailed reporting
    fn debug_cache_state(&self, context: &str) {
        if std::env::var("GAMMALOOP_DEBUG_CACHE").is_ok() {
            let validation = self.validate_cache_consistency();
            let stats = self.get_cache_stats();

            tracing::info!("🔍 DEBUG CACHE STATE at {}", context);
            tracing::info!("   Validation: {}", validation);
            tracing::info!("   Statistics: {}", stats);

            if !validation.is_valid {
                tracing::error!("   ❌ CACHE INCONSISTENCY DETECTED!");
                panic!(
                    "Cache corruption at {}: {}",
                    context, validation.diagnostics
                );
            }

            if validation.has_rotations {
                tracing::info!(
                    "   🔄 {} rotation variants from base cache_id {}",
                    validation.current_external_cache_id - validation.base_external_cache_id,
                    validation.base_external_cache_id
                );
            }

            if stats.efficiency_ratio < 0.3 {
                tracing::warn!(
                    "   ⚠️ Very low cache efficiency: {:.1}%",
                    stats.efficiency_ratio * 100.0
                );
            } else if stats.efficiency_ratio > 0.8 {
                tracing::info!(
                    "   ✅ Excellent cache efficiency: {:.1}%",
                    stats.efficiency_ratio * 100.0
                );
            }
        }
    }

    /// Get the base external cache ID for the current configuration
    fn get_base_external_cache_id(&self) -> usize;

    /// Validate cache ID consistency and return diagnostics
    fn validate_cache_consistency(&self) -> CacheValidationResult {
        let current_id = self.external_cache_id();
        let base_id = self.get_base_external_cache_id();
        let loop_id = self.loop_cache_id();

        let is_valid = base_id <= current_id;
        let has_rotations = current_id > base_id;
        let cache_efficiency = if loop_id > 0 {
            1.0 - (current_id as f64 / loop_id as f64)
        } else {
            0.0
        };

        CacheValidationResult {
            is_valid,
            current_external_cache_id: current_id,
            base_external_cache_id: base_id,
            loop_cache_id: loop_id,
            has_rotations,
            cache_efficiency,
            diagnostics: if is_valid {
                "Cache IDs are consistent".to_string()
            } else {
                format!(
                    "ERROR: Base cache ID ({}) > Current cache ID ({})",
                    base_id, current_id
                )
            },
        }
    }

    /// Get cache usage statistics for monitoring
    fn get_cache_stats(&self) -> CacheStats {
        let validation = self.validate_cache_consistency();
        CacheStats {
            total_external_increments: validation.current_external_cache_id,
            total_loop_increments: validation.loop_cache_id,
            base_configurations: validation.base_external_cache_id + 1,
            rotational_variants: validation.current_external_cache_id
                - validation.base_external_cache_id,
            efficiency_ratio: validation.cache_efficiency,
        }
    }

    fn get_group_masters(&self) -> impl Iterator<Item = &Self::G>;

    fn get_terms_mut(&mut self) -> impl Iterator<Item = &mut Self::G>;
    fn get_settings(&self) -> &RuntimeSettings;
    fn get_master_graph(&self, group_id: GroupId) -> &Self::G;
    fn get_graph_mut(&mut self, graph_id: usize) -> &mut Self::G;
    fn get_group(&self, group_id: GroupId) -> &GraphGroup;
    fn get_group_structure(&self) -> &TiVec<GroupId, GraphGroup>;
    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor<'_>;

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
                .get_master_graph(GroupId(0))
                .get_graph()
                .get_loop_number()
                * 3,
        )
    }
}

pub trait GraphTerm {
    fn evaluate<T: FloatLike>(
        &mut self,
        sample: &MomentumSample<T>,
        model: &Model,
        settings: &RuntimeSettings,
        rotation: &Rotation,
        channel_id: Option<(ChannelIndex, F<T>)>,
    ) -> Complex<F<T>>;

    fn name(&self) -> String;

    fn warm_up(&mut self, settings: &RuntimeSettings, model: &Model) -> Result<()>;
    fn get_graph(&self) -> &Graph;
    fn get_num_channels(&self) -> usize;
    fn get_num_orientations(&self) -> usize;
    fn get_tropical_sampler(&self) -> &SampleGenerator<3>;
    fn get_mut_param_builder(&mut self) -> &mut ParamBuilder<f64>;
    fn get_real_mass_vector(&self) -> EdgeVec<Option<F<f64>>>;
}

fn evaluate_all_rotations<T: FloatLike, I: GammaloopIntegrand>(
    integrand: &mut I,
    model: &Model,
    gammaloop_sample: &GammaLoopSample<T>,
) -> (Vec<Complex<F<f64>>>, Duration) {
    let rotations = integrand.get_rotations().cloned().collect_vec();

    let cache = integrand.get_settings().general.enable_cache;

    let mut loop_mom_cache_id = integrand.loop_cache_id();
    let mut external_mom_cache_id = integrand.external_cache_id();

    // rotate the momenta for the stability tests.
    let gammaloop_samples: Vec<_> = rotations
        .iter()
        .map(|rotation| {
            if rotation.is_identity() {
                return gammaloop_sample.clone();
            }
            if cache {
                loop_mom_cache_id += 1;
                external_mom_cache_id += 1;
            }
            gammaloop_sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id)
        })
        .collect();

    let start_time = std::time::Instant::now();
    let evaluation_results = gammaloop_samples
        .iter()
        .zip(rotations.iter())
        .map(|(gammaloop_sample, rotation)| {
            debug!("Evaluating rotation: {}", rotation.method);
            evaluate_single(integrand, model, gammaloop_sample, rotation)
        })
        .collect_vec();

    if cache {
        integrand.increment_loop_cache_id(rotations.len());
        integrand.increment_external_cache_id(rotations.len());
        // After evaluating all rotations, revert to base cache ID to enable cache reuse
        // for subsequent sample points with the same base external momenta
        integrand.revert_to_base_external_cache_id();
    }

    let duration = start_time.elapsed() / gammaloop_samples.len() as u32;

    (evaluation_results, duration)
}

/// Result of cache validation checks
#[derive(Debug, Clone)]
pub struct CacheValidationResult {
    pub is_valid: bool,
    pub current_external_cache_id: usize,
    pub base_external_cache_id: usize,
    pub loop_cache_id: usize,
    pub has_rotations: bool,
    pub cache_efficiency: f64,
    pub diagnostics: String,
}

impl std::fmt::Display for CacheValidationResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Cache Validation: {} | Current: {} | Base: {} | Loop: {} | Rotations: {} | Efficiency: {:.1}%",
            if self.is_valid { "✓" } else { "✗" },
            self.current_external_cache_id,
            self.base_external_cache_id,
            self.loop_cache_id,
            if self.has_rotations { "Yes" } else { "No" },
            self.cache_efficiency * 100.0
        )
    }
}

/// Cache usage statistics
#[derive(Debug, Clone)]
pub struct CacheStats {
    pub total_external_increments: usize,
    pub total_loop_increments: usize,
    pub base_configurations: usize,
    pub rotational_variants: usize,
    pub efficiency_ratio: f64,
}

impl std::fmt::Display for CacheStats {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Cache Stats: {} base configs, {} rotations, {} loop increments, {:.1}% efficiency",
            self.base_configurations,
            self.rotational_variants,
            self.total_loop_increments,
            self.efficiency_ratio * 100.0
        )
    }
}

/// Helper macro for cache debugging
#[macro_export]
macro_rules! debug_cache {
    ($integrand:expr, $msg:expr) => {
        let validation = $integrand.validate_cache_consistency();
        tracing::debug!("{}: {}", $msg, validation);
    };
}

/// Helper macro for cache monitoring with custom conditions
#[macro_export]
macro_rules! monitor_cache {
    ($integrand:expr, $condition:expr, $msg:expr) => {
        if $condition {
            let stats = $integrand.get_cache_stats();
            tracing::info!("{}: {}", $msg, stats);
        }
    };
}

/// Helper macro for missed cache hit detection in debug mode
#[macro_export]
macro_rules! debug_cache_search {
    ($integrand:expr, $msg:expr) => {
        if std::env::var("GAMMALOOP_DEBUG_CACHE").is_ok() {
            let validation = $integrand.validate_cache_consistency();
            if !validation.is_valid {
                tracing::error!(
                    "CACHE CORRUPTION DETECTED at {}: {}",
                    $msg,
                    validation.diagnostics
                );
            } else {
                tracing::debug!("Cache search at {}: {}", $msg, validation);
            }
        }
    };
}

/// Helper macro for cache efficiency warnings
#[macro_export]
macro_rules! warn_cache_efficiency {
    ($integrand:expr, $threshold:expr, $msg:expr) => {
        let stats = $integrand.get_cache_stats();
        if stats.efficiency_ratio < $threshold {
            tracing::warn!(
                "⚠️ Low cache efficiency at {}: {:.1}% (threshold: {:.1}%)",
                $msg,
                stats.efficiency_ratio * 100.0,
                $threshold * 100.0
            );
            tracing::warn!("   Cache stats: {}", stats);
        }
    };
}

fn evaluate_single<T: FloatLike, I: GammaloopIntegrand>(
    integrand: &mut I,
    model: &Model,
    gammaloop_sample: &GammaLoopSample<T>,
    rotation: &Rotation,
) -> Complex<F<f64>> {
    let settings = integrand.get_settings().clone();
    let zero = Complex::new_re(gammaloop_sample.get_default_sample().zero());
    let loop_cache_shift = 0;
    let cache = integrand.get_settings().general.enable_cache;

    let result = match &gammaloop_sample {
        GammaLoopSample::Default(sample) => integrand
            .get_terms_mut()
            .map(|term: &mut I::G| term.evaluate(sample, model, &settings, rotation, None))
            .fold(zero.clone(), |sum, term| sum + term),
        GammaLoopSample::MultiChanneling { .. } => {
            unimplemented!(
                "deprecated due to annyoing borrow issues, just set each graph to the same group"
            );
        }
        GammaLoopSample::DiscreteGraph { group_id, sample } => {
            let group = integrand.get_group(*group_id).into_iter().collect_vec(); // collect to avoid borrowing issues

            let mut res = zero.clone();

            for graph_id in group.into_iter() {
                let graph_term_res = match sample {
                    DiscreteGraphSample::Default(sample) => integrand
                        .get_graph_mut(graph_id)
                        .evaluate(sample, model, &settings, rotation, None),
                    DiscreteGraphSample::DiscreteMultiChanneling {
                        alpha,
                        channel_id,
                        sample,
                    } => integrand.get_graph_mut(graph_id).evaluate(
                        sample,
                        model,
                        &settings,
                        rotation,
                        Some((*channel_id, alpha.clone())),
                    ),
                    DiscreteGraphSample::MultiChanneling { alpha, sample } => {
                        let num_channels = integrand.get_master_graph(*group_id).get_num_channels();
                        (0..num_channels)
                            .map(ChannelIndex::from)
                            .map(|channel_index| {
                                integrand.get_graph_mut(graph_id).evaluate(
                                    sample,
                                    model,
                                    &settings,
                                    rotation,
                                    Some((channel_index, alpha.clone())),
                                )
                            })
                            .fold(zero.clone(), |sum, term| sum + term)
                    }
                    DiscreteGraphSample::Tropical(sample) => {
                        let master_graph = integrand.get_master_graph(*group_id).get_graph();

                        let energy_cache = master_graph.get_energy_cache(
                            model,
                            sample.loop_moms(),
                            sample.external_moms(),
                            &master_graph.loop_momentum_basis,
                        );

                        let prefactor = master_graph
                            .iter_loop_edges()
                            .map(|(_, edge_index, _)| edge_index)
                            .zip(
                                integrand
                                    .get_master_graph(*group_id)
                                    .get_tropical_sampler()
                                    .iter_edge_weights(),
                            )
                            .fold(sample.one(), |product, (edge_id, weight)| {
                                let energy = &energy_cache[edge_id];
                                let edge_weight = weight;
                                product * energy.powf(&F::from_f64(2. * edge_weight))
                            });

                        Complex::new_re(prefactor)
                            * integrand
                                .get_graph_mut(graph_id)
                                .evaluate(sample, model, &settings, rotation, None)
                    }
                };

                res += graph_term_res;
            }

            res
        }
    };

    if cache {
        integrand.increment_loop_cache_id(loop_cache_shift);
    }

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
                    F(integrator_settings.max_prob_ratio),
                    integrator_settings.train_on_avg,
                ))
            } else {
                continuous_grid
            }
        }
        DiscreteGraphSamplingType::DiscreteMultiChanneling(_multichanneling_settings) => {
            let continuous_grid = create_default_continous_grid(graph_term, integrator_settings);
            let lmb_channel_grid = Grid::Discrete(DiscreteGrid::new(
                (0..graph_term.get_num_channels())
                    .map(|_| Some(continuous_grid.clone()))
                    .collect_vec(),
                F(integrator_settings.max_prob_ratio),
                integrator_settings.train_on_avg,
            ));

            if settings.sample_orientations {
                Grid::Discrete(DiscreteGrid::new(
                    (0..graph_term.get_num_orientations())
                        .map(|_| Some(lmb_channel_grid.clone()))
                        .collect(),
                    F(integrator_settings.max_prob_ratio),
                    integrator_settings.train_on_avg,
                ))
            } else {
                lmb_channel_grid
            }
        }

        DiscreteGraphSamplingType::TropicalSampling(_) => {
            let dimension = get_n_dim_for_n_loop_momenta(
                &SamplingSettings::DiscreteGraphs(settings.clone()),
                graph_term.get_graph().get_loop_number(),
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
                    F(integrator_settings.max_prob_ratio),
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
        graph_term.get_graph().get_loop_number() * 3,
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
                    .get_group_masters()
                    .map(|term| {
                        Some(create_grid_for_graph(
                            term,
                            discrete_graph_sampling_settings,
                            &settings.integrator,
                        ))
                    })
                    .collect(),
                F(settings.integrator.max_prob_ratio),
                settings.integrator.train_on_avg,
            ))
        }
    }
}

fn evaluate_sample<I: GammaloopIntegrand>(
    integrand: &mut I,
    model: &Model,
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
                        debug!("f64 parameterization succeeded");
                        debug!(
                            "jacobian: {:+16e}",
                            gammaloop_sample.get_default_sample().jacobian()
                        );
                        let parameterization_time = before_parameterization.elapsed();
                        (
                            evaluate_all_rotations(integrand, model, &gammaloop_sample),
                            parameterization_time,
                        )
                    } else {
                        continue;
                    }
                }
                Precision::Quad => {
                    if let Ok(gammaloop_sample) = parameterize::<f128, I>(sample, integrand) {
                        debug!("f128 parameterization succeeded");
                        debug!(
                            "jacobian: {:+16e}",
                            gammaloop_sample.get_default_sample().jacobian()
                        );
                        let parameterization_time = before_parameterization.elapsed();
                        (
                            evaluate_all_rotations(integrand, model, &gammaloop_sample),
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

        let (average_result, is_stable) = if integrand.get_settings().stability.check_on_norm {
            stability_check_on_norm(
                integrand.get_settings(),
                &results,
                &stability_level,
                max_eval,
                wgt,
            )
        } else {
            stability_check(
                integrand.get_settings(),
                &results,
                &stability_level,
                max_eval,
                wgt,
            )
        };

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
            debug!("unstable at level: {}", stability_level.precision);
            if let Ok(gammaloop_sample) = parameterize::<f64, I>(sample, integrand) {
                let mut loop_mom_cache_id = integrand.loop_cache_id();
                let mut external_mom_cache_id = integrand.external_cache_id();
                let mut shift = 0;

                let rotated_samples: Vec<_> = integrand
                    .get_rotations()
                    .map(|rotation| {
                        if rotation.is_identity() {
                            return gammaloop_sample.clone();
                        }
                        loop_mom_cache_id += 1;
                        shift += 1;
                        external_mom_cache_id += 1;
                        gammaloop_sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id)
                    })
                    .collect();
                integrand.increment_external_cache_id(shift);
                integrand.increment_loop_cache_id(shift);

                for (sample, result) in rotated_samples.iter().zip(&results) {
                    let default_sample = sample.get_default_sample();
                    debug!(
                        "loop_moms: {}, external_moms: {}",
                        format!("{}", default_sample.loop_moms()).blue(),
                        format!("{:?}", default_sample.external_moms()).blue()
                    );

                    debug!(
                        "result of current level: {}",
                        format!("{:16e}", result).blue()
                    );
                }
            } else {
                debug!("parameterization failed");
            }
        }
    }

    debug!("result at each level:");
    for level_result in results_of_stability_levels.iter() {
        debug!(
            "level: {}. result: {}",
            format!("{}", level_result.stability_level_used).green(),
            format!("{:16e}", level_result.result).blue()
        );
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
                highest_precision: Precision::Double,
            },
        }
    }
}
