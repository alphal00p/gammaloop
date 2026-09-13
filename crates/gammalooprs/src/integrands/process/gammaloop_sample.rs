use crate::graph::GroupId;
use crate::integrands::process::{GraphTerm, SamplingMomentumSampleContext};
use crate::momentum::sample::MomentumSample;
use crate::momentum::{Rotation, ThreeMomentum};

use crate::utils::{self, F, FloatLike, global_parameterize};
use crate::{
    DependentMomentaConstructor, settings::runtime::DiscreteGraphSamplingType,
    settings::runtime::LmbChannelWeight, settings::runtime::ParameterizationMode,
    settings::runtime::ParameterizationSettings, settings::runtime::SamplingSettings,
    settings::runtime::kinematic::KinematicsSettings,
};
use color_eyre::Result;
use eyre::eyre;
use itertools::Itertools;
use momtrop::vector::Vector;
use symbolica::numerical_integration::Sample;

use super::{ProcessIntegrandImpl, SamplingChannelId, resolve_discrete_selection_for_sampling};

// discrete dimensions, continious dimensions
fn unwrap_sample<T: FloatLike>(sample: &Sample<F<f64>>) -> (Vec<usize>, Vec<F<T>>) {
    let discrete_dimensions = Vec::new();
    unwrap_sample_impl(discrete_dimensions, sample)
}

fn unwrap_sample_impl<T: FloatLike>(
    mut discrete_dimensions: Vec<usize>,
    sample: &Sample<F<f64>>,
) -> (Vec<usize>, Vec<F<T>>) {
    match sample {
        Sample::Continuous(_, xs) => {
            let xs = xs.iter().map(|x| F::from_ff64(*x)).collect();
            (discrete_dimensions, xs)
        }
        Sample::Uniform(_, discrete, xs) => {
            let xs = xs.iter().map(|x| F::from_ff64(*x)).collect();
            (discrete.clone(), xs)
        }
        Sample::Discrete(_, index, sample) => {
            discrete_dimensions.push(*index);
            unwrap_sample_impl(
                discrete_dimensions,
                sample.as_ref().expect("invalid sample structure"),
            )
        }
    }
}

/// Sample whose structure depends on the sampling settings, and enforces these settings.
#[derive(Debug, Clone)]
pub enum GammaLoopSample<T: FloatLike> {
    Default {
        sample: MomentumSample<T>,
        use_lmb_basis: bool,
    },
    Graph {
        graph_id: usize,
        sample: MomentumSample<T>,
    },
    MultiChanneling {
        alpha: F<T>,
        channel_weight: LmbChannelWeight,
        /// Unit-cube coordinates retained for canonical summed channels.
        /// Amplitude graph terms replay these coordinates through the one
        /// canonical bridge; cross-section terms retain them for the guarded
        /// legacy route until their conditional LU/t* context is available.
        sampling_coordinates: Option<Vec<F<T>>>,
        sample: MomentumSample<T>,
    },
    DiscreteGraph {
        group_id: GroupId,
        sample: DiscreteGraphSample<T>,
    },
}

impl<T: FloatLike> GammaLoopSample<T> {
    pub(crate) fn rotate(
        &self,
        rotation: &Rotation,
        loop_mom_cache_id: usize,
        external_mom_cache_id: usize,
    ) -> Self {
        if rotation.is_identity() {
            return self.clone();
        }

        match self {
            GammaLoopSample::Default {
                sample,
                use_lmb_basis,
            } => GammaLoopSample::Default {
                sample: sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id),
                use_lmb_basis: *use_lmb_basis,
            },
            GammaLoopSample::Graph { graph_id, sample } => GammaLoopSample::Graph {
                graph_id: *graph_id,
                sample: sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id),
            },
            GammaLoopSample::MultiChanneling {
                alpha,
                channel_weight,
                sampling_coordinates,
                sample,
            } => GammaLoopSample::MultiChanneling {
                alpha: alpha.clone(),
                channel_weight: *channel_weight,
                sampling_coordinates: sampling_coordinates.clone(),
                sample: sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id),
            },
            GammaLoopSample::DiscreteGraph { group_id, sample } => GammaLoopSample::DiscreteGraph {
                group_id: *group_id,
                sample: sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id),
            },
        }
    }

    #[allow(dead_code)]
    pub(crate) fn zero(&self) -> F<T> {
        match self {
            GammaLoopSample::Default { sample, .. } => sample.zero(),
            GammaLoopSample::Graph { sample, .. } => sample.zero(),
            GammaLoopSample::MultiChanneling { sample, .. } => sample.zero(),
            GammaLoopSample::DiscreteGraph { sample, .. } => sample.zero(),
        }
    }

    #[allow(dead_code)]
    pub(crate) fn one(&self) -> F<T> {
        match self {
            GammaLoopSample::Default { sample, .. } => sample.one(),
            GammaLoopSample::Graph { sample, .. } => sample.one(),
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
            GammaLoopSample::Default {
                sample,
                use_lmb_basis,
            } => GammaLoopSample::Default {
                sample: sample.cast_sample(),
                use_lmb_basis: *use_lmb_basis,
            },
            GammaLoopSample::Graph { graph_id, sample } => GammaLoopSample::Graph {
                graph_id: *graph_id,
                sample: sample.cast_sample(),
            },
            GammaLoopSample::MultiChanneling {
                alpha,
                channel_weight,
                sampling_coordinates,
                sample,
            } => GammaLoopSample::MultiChanneling {
                alpha: alpha.clone().into(),
                channel_weight: *channel_weight,
                sampling_coordinates: sampling_coordinates
                    .as_ref()
                    .map(|coordinates| coordinates.iter().cloned().map(F::<T2>::from).collect()),
                sample: sample.cast_sample(),
            },
            GammaLoopSample::DiscreteGraph { group_id, sample } => GammaLoopSample::DiscreteGraph {
                group_id: *group_id,
                sample: sample.cast_sample(),
            },
        }
    }

    #[allow(unused)]
    fn higher_precision(&self) -> GammaLoopSample<T::Higher>
    where
        T::Lower: FloatLike + Default,
        T::Higher: FloatLike + Default,
    {
        match self {
            GammaLoopSample::Default {
                sample,
                use_lmb_basis,
            } => GammaLoopSample::Default {
                sample: sample.higher_precision(),
                use_lmb_basis: *use_lmb_basis,
            },
            GammaLoopSample::Graph { graph_id, sample } => GammaLoopSample::Graph {
                graph_id: *graph_id,
                sample: sample.higher_precision(),
            },
            GammaLoopSample::MultiChanneling {
                alpha,
                channel_weight,
                sampling_coordinates,
                sample,
            } => GammaLoopSample::MultiChanneling {
                alpha: alpha.higher(),
                channel_weight: *channel_weight,
                sampling_coordinates: sampling_coordinates.as_ref().map(|coordinates| {
                    coordinates
                        .iter()
                        .map(|coordinate| coordinate.higher())
                        .collect()
                }),
                sample: sample.higher_precision(),
            },
            GammaLoopSample::DiscreteGraph { group_id, sample } => GammaLoopSample::DiscreteGraph {
                group_id: *group_id,
                sample: sample.higher_precision(),
            },
        }
    }

    #[allow(dead_code)]
    fn lower_precision(&self) -> GammaLoopSample<T::Lower>
    where
        T::Lower: FloatLike + Default,
        T::Higher: FloatLike + Default,
    {
        match self {
            GammaLoopSample::Default {
                sample,
                use_lmb_basis,
            } => GammaLoopSample::Default {
                sample: sample.lower_precision(),
                use_lmb_basis: *use_lmb_basis,
            },
            GammaLoopSample::Graph { graph_id, sample } => GammaLoopSample::Graph {
                graph_id: *graph_id,
                sample: sample.lower_precision(),
            },
            GammaLoopSample::MultiChanneling {
                alpha,
                channel_weight,
                sampling_coordinates,
                sample,
            } => GammaLoopSample::MultiChanneling {
                alpha: alpha.lower(),
                channel_weight: *channel_weight,
                sampling_coordinates: sampling_coordinates.as_ref().map(|coordinates| {
                    coordinates
                        .iter()
                        .map(|coordinate| coordinate.lower())
                        .collect()
                }),
                sample: sample.lower_precision(),
            },
            GammaLoopSample::DiscreteGraph { group_id, sample } => GammaLoopSample::DiscreteGraph {
                group_id: *group_id,
                sample: sample.lower_precision(),
            },
        }
    }

    /// Retrieve the default sample which is contained in all types
    #[inline]
    pub(crate) fn get_default_sample(&self) -> &MomentumSample<T> {
        match self {
            GammaLoopSample::Default { sample, .. } => sample,
            GammaLoopSample::Graph { sample, .. } => sample,
            GammaLoopSample::MultiChanneling { sample, .. } => sample,
            GammaLoopSample::DiscreteGraph { sample, .. } => sample.get_default_sample(),
        }
    }
}

/// This sample is used when importance sampling over graphs is used.
#[derive(Debug, Clone)]
pub enum DiscreteGraphSample<T: FloatLike> {
    Default {
        sample: MomentumSample<T>,
        use_lmb_basis: bool,
    },
    MultiChanneling {
        alpha: F<T>,
        channel_weight: LmbChannelWeight,
        sample: MomentumSample<T>,
    },
    /// This variant is equivalent to Default, but needs to be handled differently in the evaluation.
    Tropical(MomentumSample<T>),
    /// A point generated or selected through the canonical compiled
    /// sampling-channel bridge. The bridge owns the parent-frame map and its
    /// partition; graph evaluation must therefore not reinterpret this point
    /// through an LMB a second time.
    Advanced {
        channel_id: SamplingChannelId,
        /// Original unit-cube coordinates used by the canonical channel map.
        /// They are preserved so a future physical channel can run its
        /// per-sample LU/t* preparation before applying a conditional map.
        /// Direct momentum evaluations have no such coordinates and store
        /// `None`.
        sampling_coordinates: Option<Vec<F<T>>>,
        /// Direct momentum evaluations have no parameterization Jacobian
        /// boundary, so their selected partition factor is applied after the
        /// graph term is evaluated. X-space samples fold this factor into the
        /// map Jacobian and leave this field unset.
        partition_weight: Option<F<T>>,
        sample: MomentumSample<T>,
    },
}

/// Whether the selected runtime mode still uses the legacy *summed* channel
/// estimator.  Both top-level summed sampling and discrete graph sampling can
/// request this mode; keeping the predicate central prevents the nested graph
/// form from bypassing the graph-aware-channel guard below.
fn is_summed_multichanneling(settings: &SamplingSettings) -> bool {
    match settings {
        SamplingSettings::MultiChanneling(_) => true,
        SamplingSettings::DiscreteGraphs(settings) => matches!(
            &settings.sampling_type,
            DiscreteGraphSamplingType::MultiChanneling(_)
        ),
        SamplingSettings::Default(_) => false,
    }
}

impl<T: FloatLike> DiscreteGraphSample<T> {
    #[allow(dead_code)]
    pub(crate) fn zero(&self) -> F<T> {
        match self {
            DiscreteGraphSample::Default { sample, .. } => sample.zero(),
            DiscreteGraphSample::MultiChanneling { sample, .. } => sample.zero(),
            DiscreteGraphSample::Tropical(sample) => sample.zero(),
            DiscreteGraphSample::Advanced { sample, .. } => sample.zero(),
        }
    }

    pub(crate) fn one(&self) -> F<T> {
        match self {
            DiscreteGraphSample::Default { sample, .. } => sample.one(),
            DiscreteGraphSample::MultiChanneling { sample, .. } => sample.one(),
            DiscreteGraphSample::Tropical(sample) => sample.one(),
            DiscreteGraphSample::Advanced { sample, .. } => sample.one(),
        }
    }

    /// Rotation for stability checks
    #[inline]
    fn rotate(
        &self,
        rotation: &Rotation,
        loop_mom_cache_id: usize,
        external_mom_cache_id: usize,
    ) -> Self {
        match self {
            DiscreteGraphSample::Default {
                sample,
                use_lmb_basis,
            } => DiscreteGraphSample::Default {
                sample: sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id),
                use_lmb_basis: *use_lmb_basis,
            },
            DiscreteGraphSample::MultiChanneling {
                alpha,
                channel_weight,
                sample,
            } => DiscreteGraphSample::MultiChanneling {
                alpha: alpha.clone(),
                channel_weight: *channel_weight,
                sample: sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id),
            },
            DiscreteGraphSample::Tropical(sample) => DiscreteGraphSample::Tropical(sample.rotate(
                rotation,
                loop_mom_cache_id,
                external_mom_cache_id,
            )),
            DiscreteGraphSample::Advanced {
                channel_id,
                sampling_coordinates,
                partition_weight,
                sample,
            } => DiscreteGraphSample::Advanced {
                channel_id: *channel_id,
                sampling_coordinates: sampling_coordinates.clone(),
                partition_weight: partition_weight.clone(),
                sample: sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id),
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
            DiscreteGraphSample::Default {
                sample,
                use_lmb_basis,
            } => DiscreteGraphSample::Default {
                sample: sample.cast_sample(),
                use_lmb_basis: *use_lmb_basis,
            },
            DiscreteGraphSample::MultiChanneling {
                alpha,
                channel_weight,
                sample,
            } => DiscreteGraphSample::MultiChanneling {
                alpha: Into::<F<T2>>::into(alpha.clone()),
                channel_weight: *channel_weight,
                sample: sample.cast_sample(),
            },
            DiscreteGraphSample::Tropical(sample) => {
                DiscreteGraphSample::Tropical(sample.cast_sample())
            }
            DiscreteGraphSample::Advanced {
                channel_id,
                sampling_coordinates,
                partition_weight,
                sample,
            } => DiscreteGraphSample::Advanced {
                channel_id: *channel_id,
                sampling_coordinates: sampling_coordinates
                    .as_ref()
                    .map(|coordinates| coordinates.iter().cloned().map(F::<T2>::from).collect()),
                partition_weight: partition_weight
                    .as_ref()
                    .map(|weight| F::<T2>::from(weight.clone())),
                sample: sample.cast_sample(),
            },
        }
    }

    fn higher_precision(&self) -> DiscreteGraphSample<T::Higher>
    where
        T::Higher: FloatLike + Default,
        T::Lower: FloatLike + Default,
    {
        match self {
            DiscreteGraphSample::Default {
                sample,
                use_lmb_basis,
            } => DiscreteGraphSample::Default {
                sample: sample.higher_precision(),
                use_lmb_basis: *use_lmb_basis,
            },
            DiscreteGraphSample::MultiChanneling {
                alpha,
                channel_weight,
                sample,
            } => DiscreteGraphSample::MultiChanneling {
                alpha: alpha.higher(),
                channel_weight: *channel_weight,
                sample: sample.higher_precision(),
            },
            DiscreteGraphSample::Tropical(sample) => {
                DiscreteGraphSample::Tropical(sample.higher_precision())
            }
            DiscreteGraphSample::Advanced {
                channel_id,
                sampling_coordinates,
                partition_weight,
                sample,
            } => DiscreteGraphSample::Advanced {
                channel_id: *channel_id,
                sampling_coordinates: sampling_coordinates.as_ref().map(|coordinates| {
                    coordinates
                        .iter()
                        .map(|coordinate| coordinate.higher())
                        .collect()
                }),
                partition_weight: partition_weight.as_ref().map(|weight| weight.higher()),
                sample: sample.higher_precision(),
            },
        }
    }

    fn lower_precision(&self) -> DiscreteGraphSample<T::Lower>
    where
        T::Higher: FloatLike + Default,
        T::Lower: FloatLike + Default,
    {
        match self {
            DiscreteGraphSample::Default {
                sample,
                use_lmb_basis,
            } => DiscreteGraphSample::Default {
                sample: sample.lower_precision(),
                use_lmb_basis: *use_lmb_basis,
            },
            DiscreteGraphSample::MultiChanneling {
                alpha,
                channel_weight,
                sample,
            } => DiscreteGraphSample::MultiChanneling {
                alpha: alpha.lower(),
                channel_weight: *channel_weight,
                sample: sample.lower_precision(),
            },
            DiscreteGraphSample::Tropical(sample) => {
                DiscreteGraphSample::Tropical(sample.lower_precision())
            }
            DiscreteGraphSample::Advanced {
                channel_id,
                sampling_coordinates,
                partition_weight,
                sample,
            } => DiscreteGraphSample::Advanced {
                channel_id: *channel_id,
                sampling_coordinates: sampling_coordinates.as_ref().map(|coordinates| {
                    coordinates
                        .iter()
                        .map(|coordinate| coordinate.lower())
                        .collect()
                }),
                partition_weight: partition_weight.as_ref().map(|weight| weight.lower()),
                sample: sample.lower_precision(),
            },
        }
    }

    /// Retrieve the default sample which is contained in all types
    #[inline]
    fn get_default_sample(&self) -> &MomentumSample<T> {
        match self {
            DiscreteGraphSample::Default { sample, .. } => sample,
            DiscreteGraphSample::MultiChanneling { sample, .. } => sample,
            DiscreteGraphSample::Tropical(sample) => sample,
            DiscreteGraphSample::Advanced { sample, .. } => sample,
        }
    }

    /// Unit-cube coordinates retained for a canonical advanced channel.
    /// `None` identifies a direct momentum-space sample, which has no
    /// parameterization coordinates to replay for a deferred physical map.
    #[allow(dead_code)]
    pub(crate) fn sampling_coordinates(&self) -> Option<&[F<T>]> {
        match self {
            Self::Advanced {
                sampling_coordinates,
                ..
            } => sampling_coordinates.as_deref(),
            _ => None,
        }
    }
}

#[inline]
pub(crate) fn parameterize<T: FloatLike, I: ProcessIntegrandImpl>(
    sample_point: &Sample<F<f64>>,
    integrand: &mut I,
) -> Result<GammaLoopSample<T>> {
    let (discrete_indices, xs) = unwrap_sample(sample_point);
    let settings = integrand.get_settings();
    let loop_mom_cache_id = integrand.loop_cache_id();
    let external_mom_cache_id = integrand.external_cache_id();
    let dependent_momenta_constructor = integrand.get_dependent_momenta_constructor();
    let parameterization_settings = settings.sampling.get_parameterization_settings();
    if let Some(parameterization_settings) = parameterization_settings.as_ref()
        && matches!(
            &settings.sampling,
            SamplingSettings::MultiChanneling(_) | SamplingSettings::DiscreteGraphs(_)
        )
    {
        // Resolve the canonical catalogue before decoding discrete indices.
        // This turns unsupported named/surface channels into a diagnostic
        // instead of silently treating them as an empty legacy channel axis.
        for group_id in 0..integrand.get_group_structure().len() {
            let graph = integrand.get_master_graph(GroupId(group_id));
            graph.sampling_channel_ids(parameterization_settings)?;
            if is_summed_multichanneling(&settings.sampling)
                && (!matches!(&settings.sampling, SamplingSettings::MultiChanneling(_))
                    || !graph.supports_canonical_summed_sampling())
            {
                for channel_id in graph.sampling_channel_ids(parameterization_settings)? {
                    if !graph.sampling_channel_is_lmb(channel_id, parameterization_settings)? {
                        return Err(eyre!(
                            "summed sampling multichanneling cannot use graph-aware channel {}; use discrete multi-channeling",
                            channel_id.index()
                        ));
                    }
                }
            }
        }
    }
    let (group_id, orientation_id, channel_id) = resolve_discrete_selection_for_sampling(
        &settings.sampling,
        &discrete_indices,
        integrand.get_group_structure().len(),
        |group_id| Some(integrand.get_master_graph(group_id).get_num_orientations()),
        |group_id| {
            parameterization_settings
                .as_ref()
                .map(|settings| {
                    integrand
                        .get_master_graph(group_id)
                        .sampling_channel_ids(settings)
                        .map(|channel_ids| channel_ids.len())
                })
                .transpose()
        },
    )?;

    match &settings.sampling {
        SamplingSettings::Default(parameterization_settings) => Ok(GammaLoopSample::Default {
            sample: default_parametrize(
                &xs,
                dependent_momenta_constructor,
                parameterization_settings,
                &settings.kinematics,
                None,
                loop_mom_cache_id,
                external_mom_cache_id,
            ),
            use_lmb_basis: true,
        }),
        SamplingSettings::MultiChanneling(multichanneling_settings) => {
            Ok(GammaLoopSample::MultiChanneling {
                alpha: F::from_f64(multichanneling_settings.alpha),
                channel_weight: multichanneling_settings.channel_weight,
                sampling_coordinates: Some(xs.clone()),
                sample: default_parametrize(
                    &xs,
                    dependent_momenta_constructor,
                    &multichanneling_settings.parameterization_settings,
                    &settings.kinematics,
                    None,
                    loop_mom_cache_id,
                    external_mom_cache_id,
                ),
            })
        }
        SamplingSettings::DiscreteGraphs(discrete_graph_settings) => {
            let group_id = group_id.ok_or_else(|| {
                eyre!("Internal error: missing graph-group selection for discrete graph sampling.")
            })?;

            match &discrete_graph_settings.sampling_type {
                DiscreteGraphSamplingType::Default(parameterization_settings) => {
                    Ok(GammaLoopSample::DiscreteGraph {
                        group_id,
                        sample: DiscreteGraphSample::Default {
                            sample: default_parametrize(
                                &xs,
                                dependent_momenta_constructor,
                                parameterization_settings,
                                &settings.kinematics,
                                orientation_id,
                                loop_mom_cache_id,
                                external_mom_cache_id,
                            ),
                            use_lmb_basis: true,
                        },
                    })
                }
                DiscreteGraphSamplingType::MultiChanneling(multichanneling_settings) => {
                    Ok(GammaLoopSample::DiscreteGraph {
                        group_id,
                        sample: DiscreteGraphSample::MultiChanneling {
                            alpha: F::from_f64(multichanneling_settings.alpha),
                            channel_weight: multichanneling_settings.channel_weight,
                            sample: default_parametrize(
                                &xs,
                                dependent_momenta_constructor,
                                &multichanneling_settings.parameterization_settings,
                                &settings.kinematics,
                                orientation_id,
                                loop_mom_cache_id,
                                external_mom_cache_id,
                            ),
                        },
                    })
                }
                DiscreteGraphSamplingType::TropicalSampling(tropical_sampling_settings) => {
                    let graph = integrand.get_master_graph(group_id);
                    let externals = &settings
                        .kinematics
                        .externals
                        .get_dependent_externals(dependent_momenta_constructor)?;

                    let sampler = graph.get_tropical_sampler();

                    let edge_data = graph
                        .get_graph()
                        .iter_loop_edges()
                        .map(|(_, edge_id, _edge)| {
                            let mass_re = graph.get_real_mass_vector()[edge_id].map(F::from_ff64);

                            let shift = utils::compute_shift_part(
                                &graph.get_graph().loop_momentum_basis.edge_signatures[edge_id]
                                    .external,
                                externals,
                            )
                            .spatial;

                            let shift_momtrop = Vector::from_array([shift.px, shift.py, shift.pz]);

                            (mass_re, shift_momtrop)
                        })
                        .collect_vec();

                    let sampling_result_result = sampler.generate_sample_from_x_space_point(
                        &xs,
                        edge_data,
                        &tropical_sampling_settings.into_tropical_sampling_settings(),
                        None,
                    );

                    let sampling_result = match sampling_result_result {
                        Ok(sampling_result) => sampling_result,
                        Err(_) => {
                            return Err(eyre!("tropical sampling failed"));
                        }
                    };

                    let loop_moms = sampling_result
                        .loop_momenta
                        .into_iter()
                        .map(Into::<ThreeMomentum<F<T>>>::into)
                        .collect();

                    let default_sample = MomentumSample::new(
                        loop_moms,
                        loop_mom_cache_id,
                        &settings.kinematics.externals,
                        external_mom_cache_id,
                        sampling_result.jacobian,
                        dependent_momenta_constructor,
                        orientation_id,
                    )?;
                    Ok(GammaLoopSample::DiscreteGraph {
                        group_id,
                        sample: DiscreteGraphSample::Tropical(default_sample),
                    })
                }
                DiscreteGraphSamplingType::DiscreteMultiChanneling(multichanneling_settings) => {
                    let channel_id = channel_id.ok_or_else(|| {
                        eyre!(
                            "Internal error: missing channel selection for discrete multi-channeling."
                        )
                    })?;

                    let parameterization_settings =
                        &multichanneling_settings.parameterization_settings;
                    let graph = integrand.get_master_graph(group_id);
                    let externals = settings
                        .kinematics
                        .externals
                        .get_dependent_externals::<f64>(dependent_momenta_constructor)?;
                    let external_momenta = externals
                        .iter()
                        .map(|momentum| {
                            [
                                momentum.temporal.value.0,
                                momentum.spatial.px.0,
                                momentum.spatial.py.0,
                                momentum.spatial.pz.0,
                            ]
                        })
                        .collect::<Vec<_>>();
                    let bridge = graph.compile_sampling_bridge(
                        parameterization_settings,
                        settings.kinematics.e_cm,
                        &external_momenta,
                        orientation_id,
                    )?;
                    let coordinates = xs.iter().map(|x| x.clone().into_ff64().0).collect_vec();
                    let mapped = bridge.forward(channel_id, &coordinates)?;
                    let mut sample =
                        mapped.to_momentum_sample::<T>(SamplingMomentumSampleContext {
                            loop_mom_cache_id,
                            external_moms: &settings.kinematics.externals,
                            external_mom_cache_id,
                            dependent_momenta_constructor,
                            orientation: orientation_id,
                        })?;
                    let partition_weight =
                        mapped.partition.weight(channel_id.0).ok_or_else(|| {
                            eyre!(
                                "sampling channel partition has no weight for channel {}",
                                channel_id.0
                            )
                        })?;
                    if !partition_weight.is_finite() || partition_weight <= 0.0 {
                        return Err(eyre!(
                            "sampling channel partition has invalid weight {partition_weight}"
                        ));
                    }
                    sample.sample.jacobian = sample.sample.jacobian * F::from_f64(partition_weight);
                    Ok(GammaLoopSample::DiscreteGraph {
                        group_id,
                        sample: DiscreteGraphSample::Advanced {
                            channel_id,
                            sampling_coordinates: Some(xs.clone()),
                            partition_weight: None,
                            sample,
                        },
                    })
                }
            }
        }
    }
}

/// Default parametrize is basically everything except tropical sampling.
#[inline]
fn default_parametrize<T: FloatLike>(
    xs: &[F<T>],
    dependent_momenta_constructor: DependentMomentaConstructor,
    parameterization_settings: &ParameterizationSettings,
    kinematics: &KinematicsSettings,
    orientation: Option<usize>,
    loop_mom_cache_id: usize,
    external_mom_cache_id: usize,
) -> MomentumSample<T> {
    let externals = &kinematics.externals;

    let (loop_moms_vec, param_jacobian) =
        global_parameterize(xs, F::from_f64(kinematics.e_cm), parameterization_settings);

    let loop_moms = loop_moms_vec.into_iter().map(ThreeMomentum::from).collect();

    let jacobian = param_jacobian;
    // info!("Default parametrize with ext cache id:{external_mom_cache_id}");
    let mut sample = MomentumSample::new(
        loop_moms,
        loop_mom_cache_id,
        externals,
        external_mom_cache_id,
        jacobian,
        dependent_momenta_constructor,
        orientation,
    )
    .unwrap();

    if matches!(
        parameterization_settings.mode,
        ParameterizationMode::SphericalProductCommonRadial
    ) {
        sample.sample.parameterization_branch = Some(if xs[0] < F::from_f64(0.5) { 0 } else { 1 });
    }

    sample
}

#[cfg(test)]
mod tests {
    use crate::momentum::sample::{LoopMomenta, MomentumSample};
    use crate::utils::F;
    use crate::{DependentMomentaConstructor, settings::runtime::kinematic::Externals};

    use super::is_summed_multichanneling;
    use super::{DiscreteGraphSample, SamplingChannelId, unwrap_sample};
    use crate::settings::runtime::{
        DiscreteGraphSamplingSettings, DiscreteGraphSamplingType, MultiChannelingSettings,
        ParameterizationSettings, SamplingSettings,
    };

    #[test]
    fn summed_multichanneling_guard_covers_top_level_and_graph_modes() {
        assert!(is_summed_multichanneling(
            &SamplingSettings::MultiChanneling(MultiChannelingSettings::default(),)
        ));
        assert!(is_summed_multichanneling(
            &SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
                sampling_type: DiscreteGraphSamplingType::MultiChanneling(
                    MultiChannelingSettings::default(),
                ),
                ..Default::default()
            },)
        ));
        assert!(!is_summed_multichanneling(
            &SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
                sampling_type: DiscreteGraphSamplingType::DiscreteMultiChanneling(
                    MultiChannelingSettings::default(),
                ),
                ..Default::default()
            },)
        ));
        assert!(!is_summed_multichanneling(&SamplingSettings::Default(
            ParameterizationSettings::default(),
        )));
    }

    #[test]
    fn uniform_acceptance_samples_preserve_canonical_discrete_selection_order() {
        // The saved-state acceptance harness uses `Sample::Uniform` so all
        // discrete axes remain attached to one continuous point.  Verify the
        // decoder preserves the graph/orientation/channel order consumed by
        // `resolve_discrete_selection_for_sampling`.
        let sample = symbolica::numerical_integration::Sample::Uniform(
            F(1.0),
            vec![4, 2, 7],
            vec![F(0.25), F(0.75)],
        );
        let (discrete, continuous) = unwrap_sample::<f64>(&sample);
        assert_eq!(discrete, vec![4, 2, 7]);
        assert_eq!(continuous, vec![F(0.25), F(0.75)]);
    }

    #[test]
    fn advanced_sample_coordinates_survive_precision_conversion() {
        let sample = MomentumSample::new(
            LoopMomenta::from(vec![]),
            0,
            &Externals::default(),
            0,
            F(1.0),
            DependentMomentaConstructor::CrossSection,
            None,
        )
        .expect("empty cross-section sample is valid for metadata test");
        let coordinates = vec![F(0.125), F(0.625), F(0.875)];
        let advanced = DiscreteGraphSample::Advanced {
            channel_id: SamplingChannelId::from(3),
            sampling_coordinates: Some(coordinates.clone()),
            partition_weight: None,
            sample,
        };

        assert_eq!(
            advanced.sampling_coordinates(),
            Some(coordinates.as_slice())
        );
        let cast = advanced.cast_sample::<f64>();
        assert_eq!(cast.sampling_coordinates(), Some(coordinates.as_slice()));
        let higher = advanced.higher_precision();
        assert_eq!(
            higher
                .sampling_coordinates()
                .expect("higher precision keeps coordinates")
                .iter()
                .map(|coordinate| coordinate.into_f64())
                .collect::<Vec<_>>(),
            coordinates
                .iter()
                .map(|coordinate| coordinate.into_f64())
                .collect::<Vec<_>>()
        );
        let lower = higher.lower_precision();
        assert_eq!(lower.sampling_coordinates(), Some(coordinates.as_slice()));
    }

    #[test]
    fn direct_momentum_advanced_samples_have_no_replay_coordinates() {
        // The direct-momentum route intentionally cannot replay a unit-cube
        // map; a deferred physical channel must reject that route explicitly.
        let sample = MomentumSample::new(
            LoopMomenta::from(vec![]),
            0,
            &Externals::default(),
            0,
            F(1.0),
            DependentMomentaConstructor::CrossSection,
            None,
        )
        .unwrap();
        let advanced = DiscreteGraphSample::Advanced {
            channel_id: SamplingChannelId::from(0),
            sampling_coordinates: None,
            partition_weight: None,
            sample,
        };
        assert!(advanced.sampling_coordinates().is_none());
    }
}
