use crate::graph::GroupId;
use crate::integrands::evaluation::EvaluationMetaData;
use crate::integrands::process::{GraphTerm, SamplingMomentumSampleContext};
use crate::momentum::sample::MomentumSample;
use crate::momentum::{Rotation, ThreeMomentum};

use crate::utils::{self, F, FloatLike, global_parameterize};
use crate::{
    DependentMomentaConstructor, settings::runtime::DiscreteGraphSamplingType,
    settings::runtime::ParameterizationMode, settings::runtime::ParameterizationSettings,
    settings::runtime::SamplingSettings, settings::runtime::kinematic::KinematicsSettings,
};
use color_eyre::Result;
use eyre::eyre;
use itertools::Itertools;
use momtrop::vector::Vector;
use symbolica::numerical_integration::Sample;

use super::{
    ProcessIntegrandImpl, SamplingChannelId, SamplingChannelRuntimeContexts,
    resolve_discrete_selection_for_sampling, sampling_context::PreparedLUHost,
};

// discrete dimensions, continious dimensions
pub(crate) fn unwrap_sample<T: FloatLike>(sample: &Sample<F<f64>>) -> (Vec<usize>, Vec<F<T>>) {
    let discrete_dimensions = Vec::new();
    unwrap_sample_impl(discrete_dimensions, sample)
}

fn unwrap_sample_impl<T: FloatLike>(
    mut discrete_dimensions: Vec<usize>,
    sample: &Sample<F<f64>>,
) -> (Vec<usize>, Vec<F<T>>) {
    match sample {
        Sample::Continuous(_, xs) => {
            let xs = xs
                .iter()
                .map(|x| F(T::from_f64_exact_binary(x.0)))
                .collect();
            (discrete_dimensions, xs)
        }
        Sample::Uniform(_, discrete, xs) => {
            let xs = xs
                .iter()
                .map(|x| F(T::from_f64_exact_binary(x.0)))
                .collect();
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
pub(crate) enum GammaLoopSample<T: FloatLike> {
    Default {
        sample: MomentumSample<T>,
        use_lmb_basis: bool,
    },
    Graph {
        graph_id: usize,
        sample: MomentumSample<T>,
    },
    MultiChanneling {
        /// Unit-cube coordinates retained for canonical summed channels.
        /// Every graph term replays these coordinates through the canonical
        /// bridge. Per-channel Jacobians and partitions multiply the graph
        /// result, so this outer sample carries a unit Jacobian.
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
                sampling_coordinates,
                sample,
            } => GammaLoopSample::MultiChanneling {
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
pub(crate) enum DiscreteGraphSample<T: FloatLike> {
    Default {
        sample: MomentumSample<T>,
        use_lmb_basis: bool,
    },
    MultiChanneling {
        /// Unit-cube coordinates retained for canonical per-channel routing.
        /// Both process kinds replay their maps before physical evaluation;
        /// conditional cut/side maps additionally require prepared LU/t* data.
        sampling_coordinates: Option<Vec<F<T>>>,
        sample: MomentumSample<T>,
    },
    /// This variant is equivalent to Default, but needs to be handled differently in the evaluation.
    Tropical(MomentumSample<T>),
    /// A point generated or selected through the canonical compiled
    /// sampling-channel bridge. The bridge owns the parent-frame map and its
    /// partition; graph evaluation must therefore not reinterpret this point
    /// through an LMB a second time.
    SamplingChannel {
        channel_id: SamplingChannelId,
        /// Original unit-cube coordinates used by the canonical channel map.
        /// Native retries replay the original source before preparing LU/t*
        /// data and applying the conditional map; derived data is never cast.
        /// Direct momentum evaluations have no such coordinates and store
        /// `None`.
        sampling_coordinates: Option<Vec<F<T>>>,
        /// Direct momentum evaluations have no parameterization Jacobian
        /// boundary, so their selected partition factor is applied after the
        /// graph term is evaluated. X-space samples fold this factor into the
        /// map Jacobian and leave this field unset.
        partition_weight: Option<F<T>>,
        /// Initial selected-map authority, retained in its original native frame.
        /// Stability rotations rotate the completed sample; physical adoption
        /// rotates these rays once and verifies them against that sample.
        prepared_lu_hosts: Vec<PreparedLUHost<T>>,
        sample: MomentumSample<T>,
    },
}

impl<T: FloatLike> DiscreteGraphSample<T> {
    #[allow(dead_code)]
    pub(crate) fn zero(&self) -> F<T> {
        match self {
            DiscreteGraphSample::Default { sample, .. } => sample.zero(),
            DiscreteGraphSample::MultiChanneling { sample, .. } => sample.zero(),
            DiscreteGraphSample::Tropical(sample) => sample.zero(),
            DiscreteGraphSample::SamplingChannel { sample, .. } => sample.zero(),
        }
    }

    pub(crate) fn one(&self) -> F<T> {
        match self {
            DiscreteGraphSample::Default { sample, .. } => sample.one(),
            DiscreteGraphSample::MultiChanneling { sample, .. } => sample.one(),
            DiscreteGraphSample::Tropical(sample) => sample.one(),
            DiscreteGraphSample::SamplingChannel { sample, .. } => sample.one(),
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
                sampling_coordinates,
                sample,
            } => DiscreteGraphSample::MultiChanneling {
                sampling_coordinates: sampling_coordinates.clone(),
                sample: sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id),
            },
            DiscreteGraphSample::Tropical(sample) => DiscreteGraphSample::Tropical(sample.rotate(
                rotation,
                loop_mom_cache_id,
                external_mom_cache_id,
            )),
            DiscreteGraphSample::SamplingChannel {
                channel_id,
                sampling_coordinates,
                partition_weight,
                prepared_lu_hosts,
                sample,
            } => DiscreteGraphSample::SamplingChannel {
                channel_id: *channel_id,
                sampling_coordinates: sampling_coordinates.clone(),
                partition_weight: partition_weight.clone(),
                prepared_lu_hosts: prepared_lu_hosts.clone(),
                sample: sample.rotate(rotation, loop_mom_cache_id, external_mom_cache_id),
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
            DiscreteGraphSample::SamplingChannel { sample, .. } => sample,
        }
    }

    /// Unit-cube coordinates retained for a canonical sampling channel.
    /// `None` identifies a direct momentum-space sample, which has no
    /// parameterization coordinates to replay for a deferred physical map.
    #[allow(dead_code)]
    pub(crate) fn sampling_coordinates(&self) -> Option<&[F<T>]> {
        match self {
            Self::MultiChanneling {
                sampling_coordinates,
                ..
            }
            | Self::SamplingChannel {
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
    metadata: &mut EvaluationMetaData,
) -> Result<GammaLoopSample<T>> {
    integrand.prepare_sampling_precision::<T>()?;
    let (discrete_indices, xs) = unwrap_sample(sample_point);
    let settings = integrand.get_settings();
    let loop_mom_cache_id = integrand.loop_cache_id();
    let external_mom_cache_id = integrand.external_cache_id();
    let dependent_momenta_constructor = integrand.get_dependent_momenta_constructor();
    // Validate the warmed canonical catalogue before decoding discrete indices.
    // Each complete triangular map prepares its cut/LU data internally. The
    // original-draw policy owner retains discrete choices across native replay;
    // numerical map data and Jacobians are rebuilt at the requested precision.
    let (group_id, orientation_id, channel_id) = resolve_discrete_selection_for_sampling(
        &settings.sampling,
        &discrete_indices,
        integrand.get_group_structure().len(),
        |group_id| Some(integrand.get_master_graph(group_id).get_num_orientations()),
        |group_id| {
            Ok(Some(
                integrand
                    .get_master_graph(group_id)
                    .sampling_setup()
                    .sampling_bridge::<T>()?
                    .channels()
                    .len(),
            ))
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
            let mut sample = default_parametrize(
                &xs,
                dependent_momenta_constructor,
                &multichanneling_settings.parameterization_settings,
                &settings.kinematics,
                None,
                loop_mom_cache_id,
                external_mom_cache_id,
            );
            sample.sample.jacobian = sample.one();
            Ok(GammaLoopSample::MultiChanneling {
                sampling_coordinates: Some(xs.clone()),
                sample,
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
                    let mut sample = default_parametrize(
                        &xs,
                        dependent_momenta_constructor,
                        &multichanneling_settings.parameterization_settings,
                        &settings.kinematics,
                        orientation_id,
                        loop_mom_cache_id,
                        external_mom_cache_id,
                    );
                    sample.sample.jacobian = sample.one();
                    Ok(GammaLoopSample::DiscreteGraph {
                        group_id,
                        sample: DiscreteGraphSample::MultiChanneling {
                            sampling_coordinates: Some(xs.clone()),
                            sample,
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
                DiscreteGraphSamplingType::SamplingMultiChanneling(_) => {
                    let channel_id = channel_id.ok_or_else(|| {
                        eyre!(
                            "Internal error: missing channel selection for discrete multi-channeling."
                        )
                    })?;

                    let graph = integrand.get_master_graph(group_id);
                    let bridge = graph.sampling_setup().sampling_bridge::<T>()?;
                    let coordinates = xs.iter().map(|x| x.0.clone()).collect_vec();
                    let mut contexts = SamplingChannelRuntimeContexts::for_draw(
                        bridge.channels().len(),
                        integrand.get_group(group_id).master(),
                        channel_id,
                        metadata,
                    );
                    let mapped = bridge.forward_with_runtime_contexts(
                        channel_id,
                        &coordinates,
                        &mut contexts,
                    )?;
                    // The mapped point and physical sample use the same native precision.
                    let mut sample = mapped.to_momentum_sample(SamplingMomentumSampleContext {
                        loop_mom_cache_id,
                        external_moms: &settings.kinematics.externals,
                        external_mom_cache_id,
                        dependent_momenta_constructor,
                        orientation: orientation_id,
                    })?;
                    sample.sample.jacobian = F(mapped.selected_factor()?);
                    Ok(GammaLoopSample::DiscreteGraph {
                        group_id,
                        sample: DiscreteGraphSample::SamplingChannel {
                            channel_id,
                            sampling_coordinates: Some(xs.clone()),
                            partition_weight: None,
                            prepared_lu_hosts: mapped.prepared_lu_hosts,
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

    use super::{DiscreteGraphSample, SamplingChannelId, unwrap_sample};
    use crate::settings::runtime::{
        DiscreteGraphSamplingSettings, DiscreteGraphSamplingType, MultiChannelingSettings,
        ParameterizationSettings, SamplingSettings,
    };

    #[test]
    fn sampling_channel_guard_covers_summed_and_monte_carlo_modes() {
        assert!(
            SamplingSettings::MultiChanneling(MultiChannelingSettings::default())
                .uses_sampling_channels()
        );
        for sampling_type in [
            DiscreteGraphSamplingType::MultiChanneling(MultiChannelingSettings::default()),
            DiscreteGraphSamplingType::SamplingMultiChanneling(MultiChannelingSettings::default()),
        ] {
            assert!(
                SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings {
                    sampling_type,
                    ..Default::default()
                })
                .uses_sampling_channels()
            );
        }
        assert!(
            !SamplingSettings::Default(ParameterizationSettings::default())
                .uses_sampling_channels()
        );
        assert!(
            !SamplingSettings::DiscreteGraphs(DiscreteGraphSamplingSettings::default())
                .uses_sampling_channels()
        );
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
    fn sampling_channel_coordinates_replay_at_native_precision() {
        use crate::utils::{ArbPrec, FloatLike, QuadFloat};

        // Cast the original sample coordinates to each precision. Prepared
        // rays, roots and Jacobians must instead be rebuilt from this source.
        let coordinates = vec![F(0.1), F(0.625), F(0.875)];
        let source = symbolica::numerical_integration::Sample::Uniform(
            F(1.0),
            vec![4, 2, 3],
            coordinates.clone(),
        );
        fn check<T: FloatLike>(source: &symbolica::numerical_integration::Sample<F<f64>>) {
            let (selection, replay) = unwrap_sample::<T>(source);
            assert_eq!(selection, vec![4, 2, 3]);
            let (_, original) = unwrap_sample::<f64>(source);
            for (native, original) in replay.iter().zip(original) {
                assert_eq!(native.0, T::from_f64_exact_binary(original.0));
            }
        }
        check::<f64>(&source);
        check::<QuadFloat>(&source);
        check::<ArbPrec>(&source);
    }

    #[test]
    fn direct_momentum_sampling_channels_have_no_replay_coordinates() {
        // Direct momentum inputs intentionally have no cube to replay. Their
        // channel density comes from the compiled inverse at the supplied point.
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
        let sampling_channel = DiscreteGraphSample::SamplingChannel {
            channel_id: SamplingChannelId::from(0),
            sampling_coordinates: None,
            partition_weight: None,
            prepared_lu_hosts: vec![],
            sample,
        };
        assert!(sampling_channel.sampling_coordinates().is_none());
    }

    #[test]
    fn sampling_channel_preserves_canonical_identity_and_parent() {
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
        let sampling_channel = DiscreteGraphSample::SamplingChannel {
            channel_id: SamplingChannelId::from(11),
            sampling_coordinates: Some(coordinates.clone()),
            partition_weight: None,
            prepared_lu_hosts: vec![],
            sample: sample.clone(),
        };
        assert_eq!(
            sampling_channel.sampling_coordinates(),
            Some(coordinates.as_slice())
        );
        let DiscreteGraphSample::SamplingChannel {
            channel_id,
            sample: parent,
            ..
        } = sampling_channel
        else {
            unreachable!()
        };
        assert_eq!(channel_id, SamplingChannelId::from(11));
        assert_eq!(parent.sample.jacobian, sample.sample.jacobian);
        assert_eq!(
            parent.sample.loop_mom_cache_id,
            sample.sample.loop_mom_cache_id
        );
    }

    #[test]
    fn sampling_map_rejects_non_unit_original_coordinates() {
        use crate::integrands::process::{
            SamplingMapComponent, SurfaceRadialMap, sampling_context::SamplingMapContext,
        };
        let map = SurfaceRadialMap::new(3, vec![0.0; 3], None, 1.0, 1.0).unwrap();
        for coordinates in [[0.2, 1.2, 0.3], [0.2, f64::NAN, 0.3]] {
            assert!(
                SamplingMapComponent::forward(
                    &map,
                    &coordinates,
                    &mut SamplingMapContext::detached(&[]),
                )
                .is_err()
            );
        }
    }

    #[test]
    fn summed_sampling_retains_coordinates_for_canonical_amplitude_routing() {
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
        let summed = DiscreteGraphSample::MultiChanneling {
            sampling_coordinates: Some(coordinates.clone()),
            sample,
        };
        assert_eq!(summed.sampling_coordinates(), Some(coordinates.as_slice()));
    }
}
