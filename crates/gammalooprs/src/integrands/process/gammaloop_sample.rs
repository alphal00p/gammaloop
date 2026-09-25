use crate::graph::GroupId;
use crate::integrands::evaluation::EvaluationMetaData;
use crate::integrands::process::{GraphTerm, SamplingMomentumSampleContext};
use crate::momentum::sample::MomentumSample;
use crate::momentum::{Rotation, ThreeMomentum};
use crate::processes::CutGroupId;
use crate::subtraction::lu_counterterm::LUSharedOverlaps;
use crate::utils::ArbPrec;
use std::sync::Arc;
use typed_index_collections::TiVec;

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
    resolve_discrete_selection_for_sampling, sampling::context::PreparedLUHost,
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

/// A completed draw in graph-parent coordinates, with the existing event-group
/// boundaries. Every row already contains its map/partition factor; physical
/// retries materialize these rows from one immutable canonical draw.
#[derive(Debug, Clone)]
pub(crate) struct GammaLoopSample<T: FloatLike> {
    pub(crate) groups: Vec<(Option<GroupId>, Vec<DiscreteGraphSample<T>>)>,
}

/// One graph/channel contribution. The bridge owns its parent-frame map, so
/// physical evaluation must not reinterpret this point through an LMB again.
#[derive(Debug, Clone)]
pub(crate) struct DiscreteGraphSample<T: FloatLike> {
    pub(crate) graph_id: usize,
    pub(crate) channel_id: Option<SamplingChannelId>,
    /// Initial selected-map authority remains unrotated. Physical adoption
    /// rotates its rays once and verifies them against the completed sample.
    pub(crate) prepared_lu_hosts: Vec<PreparedLUHost<T>>,
    /// Physical center decisions belong to the original point and accepted cut
    /// set. Native rows retain this unrotated Arb authority across all retries.
    pub(crate) physical_overlaps: Option<Arc<TiVec<CutGroupId, Option<LUSharedOverlaps<ArbPrec>>>>>,
    pub(crate) sample: MomentumSample<T>,
    /// Tropical compensation historically multiplies only the integrand, while
    /// map/partition factors also multiply events and reference moments.
    pub(crate) integrand_prefactor: F<T>,
}

impl<T: FloatLike> GammaLoopSample<T> {
    /// Reserve the same fraction of the requested physical accuracy for map
    /// agreement, canonical-factor materialization and host adoption.
    pub(crate) fn relative_accuracy_budget(settings: &crate::settings::RuntimeSettings) -> f64 {
        settings
            .stability
            .levels
            .iter()
            .filter(|level| T::sampling_precision() == level.precision.into())
            .map(|level| {
                level
                    .required_precision_for_re
                    .min(level.required_precision_for_im)
            })
            .reduce(f64::min)
            .map(|tolerance| 0.1 * tolerance)
            .unwrap_or_else(|| {
                if T::sampling_precision() == crate::utils::SamplingPrecision::Arb {
                    // Forced Arb uses the existing default physical requirement
                    // when no Arb level is configured; source admission reserves it too.
                    let forced = super::create_stability_iterator(&settings.stability, true);
                    0.1 * forced[0]
                        .required_precision_for_re
                        .min(forced[0].required_precision_for_im)
                } else {
                    F::<T>::default().epsilon().sqrt().into_ff64().0
                }
            })
    }

    /// Source accuracy covers every permitted physical lane, including a later
    /// rescue. Invalid requirements cannot silently disappear in a minimum.
    pub(crate) fn source_accuracy_budget(
        settings: &crate::settings::RuntimeSettings,
    ) -> Result<f64> {
        if settings.stability.levels.is_empty() {
            return Err(eyre!(
                "sampling source requires at least one stability level"
            ));
        }
        // The precise API can force Arb even when no Arb lane is configured;
        // reserve its existing default requirement as well, before any draw.
        let forced = super::create_stability_iterator(&settings.stability, true);
        let mut budget = f64::INFINITY;
        for level in settings.stability.levels.iter().chain(&forced) {
            for requirement in [
                level.required_precision_for_re,
                level.required_precision_for_im,
            ] {
                let allocated = 0.1 * requirement;
                if !requirement.is_finite() || !allocated.is_finite() || allocated <= 0.0 {
                    return Err(eyre!(
                        "sampling source requires positive finite accuracy for every stability component; {} has requirement {requirement}",
                        level.precision
                    ));
                }
                budget = budget.min(allocated);
            }
        }
        Ok(budget)
    }

    pub(crate) fn rotate(
        &self,
        rotation: &Rotation,
        loop_mom_cache_id: usize,
        external_mom_cache_id: usize,
    ) -> Self {
        if rotation.is_identity() {
            return self.clone();
        }
        let mut result = self.clone();
        for (index, row) in result
            .groups
            .iter_mut()
            .flat_map(|(_, rows)| rows)
            .enumerate()
        {
            row.sample =
                row.sample
                    .rotate(rotation, loop_mom_cache_id + index, external_mom_cache_id);
        }
        result
    }

    #[allow(dead_code)]
    pub(crate) fn zero(&self) -> F<T> {
        self.get_default_sample().zero()
    }

    #[allow(dead_code)]
    pub(crate) fn one(&self) -> F<T> {
        self.get_default_sample().one()
    }

    /// The first completed row supplies common external/orientation metadata.
    #[inline]
    pub(crate) fn get_default_sample(&self) -> &MomentumSample<T> {
        &self.groups[0].1[0].sample
    }

    pub(crate) fn row_count(&self) -> usize {
        self.groups.iter().map(|(_, rows)| rows.len()).sum()
    }

    /// Default samples start in a selected LMB only for cube input. Resolve that
    /// affine reinterpretation before any physical/reference stability body.
    pub(crate) fn from_default<I: ProcessIntegrandImpl>(
        integrand: &I,
        sample: MomentumSample<T>,
        use_lmb_basis: bool,
        group_id: Option<GroupId>,
    ) -> Result<Self> {
        let groups = if let Some(group_id) = group_id {
            vec![(
                Some(group_id),
                integrand.get_group(group_id).into_iter().collect_vec(),
            )]
        } else if integrand.groups_default_sample_events_by_graph_group() {
            integrand
                .get_group_structure()
                .iter_enumerated()
                .map(|(id, group)| (Some(id), group.into_iter().collect_vec()))
                .collect()
        } else {
            (0..integrand.graph_count())
                .map(|id| (None, vec![id]))
                .collect()
        };
        let groups = groups
            .into_iter()
            .map(|(group_id, graph_ids)| {
                let rows = graph_ids
                    .into_iter()
                    .map(|graph_id| {
                        let lmb = super::selected_lmb_basis_for_default_sampling(
                            integrand,
                            graph_id,
                            use_lmb_basis,
                        )?;
                        let sample = match lmb {
                            Some(lmb) => integrand
                                .get_graph(graph_id)
                                .sampling_setup()
                                .reinterpret_loop_momenta_for_lmb(
                                    lmb,
                                    &sample,
                                    sample.sample.loop_mom_cache_id,
                                ),
                            None => sample.clone(),
                        };
                        Ok(DiscreteGraphSample {
                            graph_id,
                            channel_id: None,
                            prepared_lu_hosts: vec![],
                            physical_overlaps: None,
                            integrand_prefactor: sample.one(),
                            sample,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok((group_id, rows))
            })
            .collect::<Result<Vec<_>>>()?;
        if groups.is_empty() {
            return Err(eyre!("Cannot prepare an integrand with no graph terms."));
        }
        Ok(Self { groups })
    }

    pub(crate) fn from_selected<I: ProcessIntegrandImpl>(
        integrand: &I,
        group_id: GroupId,
        channel_id: SamplingChannelId,
        sample: MomentumSample<T>,
        prepared_lu_hosts: Vec<PreparedLUHost<T>>,
    ) -> Result<Self> {
        let master = integrand.get_group(group_id).master();
        if prepared_lu_hosts.iter().any(|host| {
            host.source.graph_id != master
                || host.source.generating_channel != channel_id
                || host.source.target_channel != channel_id
        }) {
            return Err(eyre!(
                "selected LU host records do not belong to graph-group master {master} and generating channel {channel_id:?}"
            ));
        }
        let rows = integrand
            .get_group(group_id)
            .into_iter()
            .map(|graph_id| DiscreteGraphSample {
                graph_id,
                channel_id: Some(channel_id),
                sample: sample.clone(),
                integrand_prefactor: sample.one(),
                // Foreign group members have their own physical cuts and cannot
                // adopt the master's ray, even when cut IDs happen to coincide.
                physical_overlaps: None,
                prepared_lu_hosts: if graph_id == master {
                    prepared_lu_hosts.clone()
                } else {
                    vec![]
                },
            })
            .collect();
        Ok(Self {
            groups: vec![(Some(group_id), rows)],
        })
    }

    pub(crate) fn from_tropical<I: ProcessIntegrandImpl>(
        integrand: &I,
        group_id: GroupId,
        sample: MomentumSample<T>,
    ) -> Result<Self> {
        let graph = integrand.get_master_graph(group_id);
        let master = graph.get_graph();
        let masses = graph.get_real_mass_vector()?;
        let prefactor = master
            .iter_loop_edges()
            .zip(graph.get_tropical_sampler().iter_edge_weights())
            .fold(sample.one(), |product, ((_, edge_id, _), weight)| {
                let momentum = master.loop_momentum_basis.edge_signatures[edge_id]
                    .compute_four_momentum_from_three(sample.loop_moms(), sample.external_moms());
                let energy = momentum
                    .spatial
                    .on_shell_energy(masses[edge_id].map(F::from_ff64))
                    .value;
                product * energy.powf(&F::from_f64(2. * weight))
            });
        let mut draw = Self::from_default(integrand, sample, false, Some(group_id))?;
        for row in &mut draw.groups[0].1 {
            row.integrand_prefactor = prefactor.clone();
        }
        Ok(draw)
    }

    /// Freeze the completed numerical proposal once, while rebuilding physical
    /// external data from its original Arb owner. Retained source hosts are
    /// embedded exactly; no map, inverse or root is replayed.
    pub(crate) fn into_canonical(
        self,
        externals: &crate::settings::runtime::kinematic::Externals,
        dependent_momenta_constructor: DependentMomentaConstructor,
    ) -> Result<GammaLoopSample<ArbPrec>> {
        let groups = self.groups.into_iter().map(|(group_id, rows)| {
            let rows = rows.into_iter().map(|row| {
                let native = row.sample.sample;
                if native.dual_loop_moms.is_some() || row.physical_overlaps.is_some() {
                    return Err(eyre!("fixed source supports only scalar rows before physical overlap preparation"));
                }
                let momenta = native.loop_moms.iter().map(|momentum| {
                    Ok(ThreeMomentum::new(momentum.px.to_arb_exact()?, momentum.py.to_arb_exact()?, momentum.pz.to_arb_exact()?))
                }).collect::<Result<_>>()?;
                let mut sample = MomentumSample::<ArbPrec>::new(
                    momenta, native.loop_mom_cache_id,
                    externals,
                    native.external_mom_cache_id, native.jacobian.to_arb_exact()?,
                    dependent_momenta_constructor, native.orientation,
                )?;
                sample.sample.loop_mom_base_cache_id = native.loop_mom_base_cache_id;
                sample.sample.external_mom_base_cache_id = native.external_mom_base_cache_id;
                sample.sample.parameterization_branch = native.parameterization_branch;
                Ok(DiscreteGraphSample {
                    graph_id: row.graph_id, channel_id: row.channel_id, sample,
                    integrand_prefactor: row.integrand_prefactor.to_arb_exact()?,
                    physical_overlaps: None,
                    prepared_lu_hosts: row.prepared_lu_hosts.iter().map(PreparedLUHost::to_arb_exact).collect::<Result<_>>()?,
                })
            }).collect::<Result<_>>()?;
            Ok((group_id, rows))
        }).collect::<Result<_>>()?;
        Ok(GammaLoopSample { groups })
    }
}

impl GammaLoopSample<ArbPrec> {
    pub(crate) fn prepare_physical_overlaps<I: ProcessIntegrandImpl>(
        &mut self,
        integrand: &mut I,
        model: &crate::model::Model,
        metadata: &mut EvaluationMetaData,
    ) -> Result<()> {
        let settings = integrand.get_settings().clone();
        if settings.subtraction.disable_threshold_subtraction {
            return Ok(());
        }
        let rotation = Rotation::new(crate::momentum::RotationMethod::Identity);
        let mut runtime = integrand.take_event_processing_runtime();
        let history = std::mem::take(&mut metadata.radial_root_diagnostics);
        let started = std::time::Instant::now();
        let result = (|| {
            for row in self.groups.iter_mut().flat_map(|(_, rows)| rows) {
                row.physical_overlaps = integrand
                    .get_graph_mut(row.graph_id)
                    .prepare_physical_overlaps(
                        &row.sample,
                        super::GraphTermEvaluationContext {
                            model,
                            settings: &settings,
                            event_processing_runtime: runtime.as_mut(),
                            rotation: &rotation,
                            evaluation_metadata: metadata,
                            // An explicitly channel-dependent selector retains
                            // its user-requested metadata; geometry never uses it.
                            sampling_channel: row.channel_id,
                            graph_id: row.graph_id,
                            prepared_lu_hosts: &[],
                            canonical_sample: None,
                            sampling_accuracy_budget: Self::relative_accuracy_budget(&settings),
                        },
                    )?
                    .map(Arc::new);
            }
            Ok(())
        })();
        let elapsed = started.elapsed();
        metadata.integrand_evaluation_time += elapsed;
        metadata.canonical_physical_preparation_time += elapsed;
        // These roots describe fixed physical preparation, not native retry
        // occurrences. Selectors must not consume observable or event counts.
        metadata.radial_root_diagnostics = history;
        integrand.restore_event_processing_runtime(runtime);
        result
    }

    /// Convert directly from the retained canonical output, never from an
    /// earlier physical lane. The combined factor is not recomputed from
    /// separately rounded map and partition values.
    pub(crate) fn materialize<T: FloatLike>(&self, tolerance: f64) -> Result<GammaLoopSample<T>> {
        use super::sampling::maps::SamplingEvaluationError;
        if !tolerance.is_finite() || tolerance <= 0.0 {
            return Err(eyre!(
                "sampling materialization accuracy must be positive and finite"
            ));
        }
        let groups = self
            .groups
            .iter()
            .map(|(group_id, rows)| {
                let rows = rows
                    .iter()
                    .map(|row| {
                        let sample = row.sample.materialize::<T>()?;
                        let original = &row.sample.sample.jacobian;
                        if original <= &original.zero() {
                            return Err(SamplingEvaluationError::Unrepresentable {
                                operation: "canonical row factor",
                                detail: format!(
                                    "graph {} channel {:?}: factor must be positive",
                                    row.graph_id, row.channel_id
                                ),
                            }
                            .into());
                        }
                        // The product of these two relative conversion bounds stays
                        // within the one row budget, including tropical compensation.
                        let factor_tolerance = 0.25 * tolerance.min(1.0);
                        sample
                            .sample
                            .jacobian
                            .verify_arb_materialization(&original.0, factor_tolerance)?;
                        let integrand_prefactor = F::<T>::from_arb(&row.integrand_prefactor.0)?;
                        integrand_prefactor.verify_arb_materialization(
                            &row.integrand_prefactor.0,
                            factor_tolerance,
                        )?;
                        Ok(DiscreteGraphSample {
                            graph_id: row.graph_id,
                            channel_id: row.channel_id,
                            sample,
                            integrand_prefactor,
                            physical_overlaps: row.physical_overlaps.clone(),
                            prepared_lu_hosts: row
                                .prepared_lu_hosts
                                .iter()
                                .map(|host| host.materialize::<T>(tolerance))
                                .collect::<Result<_>>()?,
                        })
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok((*group_id, rows))
            })
            .collect::<Result<_>>()?;
        Ok(GammaLoopSample { groups })
    }
}

#[inline]
pub(crate) fn parameterize<T: FloatLike, I: ProcessIntegrandImpl>(
    sample_point: &Sample<F<f64>>,
    integrand: &mut I,
    metadata: &mut EvaluationMetaData,
) -> Result<GammaLoopSample<T>> {
    integrand.prepare_sampling_precision::<T>()?;
    let (discrete_indices, xs) = unwrap_sample::<T>(sample_point);
    let settings = integrand.get_settings();
    let loop_mom_cache_id = integrand.loop_cache_id();
    let external_mom_cache_id = integrand.external_cache_id();
    let dependent_momenta_constructor = integrand.get_dependent_momenta_constructor();
    // Validate the warmed catalogue before decoding the existing discrete axes.
    // Complete triangular maps prepare cut/LU data internally. Production calls
    // this once at fixed precision and retains every completed row.
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
    let default = |parameterization: &ParameterizationSettings| {
        default_parametrize(
            &xs,
            dependent_momenta_constructor,
            parameterization,
            &settings.kinematics,
            orientation_id,
            loop_mom_cache_id,
            external_mom_cache_id,
        )
    };
    match &settings.sampling {
        SamplingSettings::Default(parameterization) => {
            GammaLoopSample::from_default(integrand, default(parameterization), true, None)
        }
        SamplingSettings::DiscreteGraphs(discrete)
            if matches!(
                discrete.sampling_type,
                DiscreteGraphSamplingType::Default(_)
            ) =>
        {
            GammaLoopSample::from_default(
                integrand,
                default(&settings.sampling.get_parameterization_settings().unwrap()),
                true,
                group_id,
            )
        }
        SamplingSettings::DiscreteGraphs(discrete)
            if matches!(
                discrete.sampling_type,
                DiscreteGraphSamplingType::TropicalSampling(_)
            ) =>
        {
            let group_id =
                group_id.ok_or_else(|| eyre!("missing tropical graph-group selection"))?;
            let DiscreteGraphSamplingType::TropicalSampling(tropical) = &discrete.sampling_type
            else {
                unreachable!()
            };
            let graph = integrand.get_master_graph(group_id);
            let externals = settings
                .kinematics
                .externals
                .get_dependent_externals(dependent_momenta_constructor)?;
            let masses = graph.get_real_mass_vector()?;
            let edge_data = graph
                .get_graph()
                .iter_loop_edges()
                .map(|(_, edge_id, _)| {
                    let mass = masses[edge_id].map(F::from_ff64);
                    let shift = utils::compute_shift_part(
                        &graph.get_graph().loop_momentum_basis.edge_signatures[edge_id].external,
                        &externals,
                    )
                    .spatial;
                    (mass, Vector::from_array([shift.px, shift.py, shift.pz]))
                })
                .collect_vec();
            let sampled = graph
                .get_tropical_sampler()
                .generate_sample_from_x_space_point(
                    &xs,
                    edge_data,
                    &tropical.into_tropical_sampling_settings(),
                    None,
                )
                .map_err(|_| eyre!("tropical sampling failed"))?;
            let sample = MomentumSample::new(
                sampled
                    .loop_momenta
                    .into_iter()
                    .map(Into::<ThreeMomentum<F<T>>>::into)
                    .collect(),
                loop_mom_cache_id,
                &settings.kinematics.externals,
                external_mom_cache_id,
                sampled.jacobian,
                dependent_momenta_constructor,
                orientation_id,
            )?;
            GammaLoopSample::from_tropical(integrand, group_id, sample)
        }
        _ => {
            let coordinates = xs.iter().map(|x| x.0.clone()).collect_vec();
            let groups = group_id.map(|id| vec![id]).unwrap_or_else(|| {
                integrand
                    .get_group_structure()
                    .iter_enumerated()
                    .map(|(id, _)| id)
                    .collect()
            });
            let mut prepared_groups = Vec::with_capacity(groups.len());
            for group in groups {
                let mut rows = Vec::new();
                let graph_id = integrand.get_group(group).master();
                let bridge = integrand
                    .get_graph(graph_id)
                    .sampling_setup()
                    .sampling_bridge::<T>()?;
                let channels = channel_id.map(|id| vec![id]).unwrap_or_else(|| {
                    (0..bridge.channels().len())
                        .map(SamplingChannelId::from)
                        .collect()
                });
                // Explicit sums contribute J_c(x) w_c(T_c(x)) f(T_c(x)).
                // Monte Carlo supplies inverse channel probability separately;
                // there is no channel-count multiplier on these rows.
                for selected in channels {
                    let mut contexts = SamplingChannelRuntimeContexts::for_draw(
                        bridge.channels().len(),
                        graph_id,
                        selected,
                        metadata,
                    );
                    let mapped = bridge.forward_with_runtime_contexts(
                        selected,
                        &coordinates,
                        &mut contexts,
                    )?;
                    let mut sample = mapped.to_momentum_sample(SamplingMomentumSampleContext {
                        loop_mom_cache_id,
                        external_moms: &settings.kinematics.externals,
                        external_mom_cache_id,
                        dependent_momenta_constructor,
                        orientation: orientation_id,
                    })?;
                    sample.sample.jacobian = F(mapped.selected_factor()?);
                    // Both execution modes evaluate a whole graph group at
                    // the master's one physical point for this channel.
                    // Graph members must cancel before taking its absolute
                    // value; a new channel supplies a distinct point.
                    let mut selected_draw = GammaLoopSample::from_selected(
                        integrand,
                        group,
                        selected,
                        sample,
                        mapped.prepared_lu_hosts,
                    )?;
                    rows.append(&mut selected_draw.groups[0].1);
                }
                prepared_groups.push((Some(group), rows));
            }
            if prepared_groups.is_empty() || prepared_groups.iter().any(|(_, rows)| rows.is_empty())
            {
                return Err(eyre!("sampling draw has no graph/channel rows"));
            }
            Ok(GammaLoopSample {
                groups: prepared_groups,
            })
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

    use super::{DiscreteGraphSample, GammaLoopSample, SamplingChannelId, unwrap_sample};
    use crate::settings::runtime::{
        DiscreteGraphSamplingSettings, DiscreteGraphSamplingType, MultiChannelingSettings,
        ParameterizationSettings, SamplingSettings,
    };

    #[test]
    fn fixed_source_budget_reserves_all_lanes_and_rejects_invalid_requests() {
        use crate::settings::{RuntimeSettings, runtime::StabilityLevelSetting};
        let mut settings = RuntimeSettings::default();
        settings.stability.levels = vec![
            StabilityLevelSetting::default_double(),
            StabilityLevelSetting::default_quad(),
            StabilityLevelSetting::default_arb(),
        ];
        settings.stability.levels[0].required_precision_for_im = 1.0e-12;
        settings.stability.levels[2].required_precision_for_re = 1.0e-20;
        assert_eq!(
            GammaLoopSample::<crate::utils::QuadFloat>::source_accuracy_budget(&settings).unwrap(),
            0.1 * 1.0e-20
        );
        assert_eq!(
            GammaLoopSample::<crate::utils::ArbPrec>::relative_accuracy_budget(&settings),
            0.1 * 1.0e-20
        );
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
            settings.stability.levels[0].required_precision_for_im = invalid;
            assert!(
                GammaLoopSample::<crate::utils::ArbPrec>::source_accuracy_budget(&settings)
                    .is_err(),
                "{invalid}"
            );
        }
        settings.stability.levels = vec![StabilityLevelSetting {
            required_precision_for_re: 0.1,
            required_precision_for_im: 0.1,
            ..StabilityLevelSetting::default_double()
        }];
        assert_eq!(
            GammaLoopSample::<crate::utils::ArbPrec>::source_accuracy_budget(&settings).unwrap(),
            0.1 * 1.0e-5
        );
        assert_eq!(
            GammaLoopSample::<crate::utils::ArbPrec>::relative_accuracy_budget(&settings),
            0.1 * 1.0e-5
        );
        assert_eq!(
            GammaLoopSample::<crate::utils::QuadFloat>::relative_accuracy_budget(&settings),
            F::<crate::utils::QuadFloat>::default()
                .epsilon()
                .sqrt()
                .into_ff64()
                .0
        );
        settings.stability.levels.clear();
        assert!(
            GammaLoopSample::<crate::utils::ArbPrec>::source_accuracy_budget(&settings).is_err()
        );
    }

    #[test]
    fn fixed_source_rows_preserve_exact_outputs_and_original_arb_externals() {
        use crate::{
            graph::GroupId,
            momentum::{ExternalMomenta, Rotation, RotationMethod, ThreeMomentum},
            utils::{ArbPrec, FloatLike, QuadFloat, SamplingFloat},
        };
        use symbolica::domains::float::DoubleFloat;

        fn check<T: FloatLike>(low: F<T>) {
            let mut externals = Externals::default();
            let Externals::Constant { momenta, .. } = &mut externals;
            *momenta = [[2.0, 1.0, 0.0, 0.0], [2.0, -1.0, 0.0, 0.0]]
                .map(|momentum| ExternalMomenta::Independent(momentum.map(F)))
                .to_vec();
            let mut original = externals
                .get_dependent_externals::<ArbPrec>(DependentMomentaConstructor::CrossSection)
                .unwrap();
            let one = F::<ArbPrec>::default().one();
            // These separated terms exceed both Quad limbs and the 256-bit source.
            let precise = &one + one.from_i64(2).powi(-80) + one.from_i64(2).powi(-320);
            original[crate::momentum::sample::ExternalIndex(0)]
                .spatial
                .px = precise.clone();
            assert_ne!(
                F::<T>::from_arb(&precise.0)
                    .unwrap()
                    .to_arb_exact()
                    .unwrap(),
                precise
            );
            let Externals::Constant { arb_cache, .. } = &mut externals;
            arb_cache.set(original.clone());
            let mut sample = MomentumSample::new(
                vec![ThreeMomentum::new(
                    low.clone(),
                    F::<T>::default().zero(),
                    F::<T>::default().one(),
                )]
                .into(),
                17,
                &externals,
                23,
                F::<T>::default().from_usize(3),
                DependentMomentaConstructor::CrossSection,
                Some(7),
            )
            .unwrap();
            sample.sample.loop_mom_base_cache_id = 11;
            sample.sample.external_mom_base_cache_id = 19;
            sample.sample.parameterization_branch = Some(2);
            let rows = (0..3)
                .map(|index| DiscreteGraphSample {
                    graph_id: index / 2,
                    channel_id: Some(SamplingChannelId(index % 2)),
                    sample: sample.clone(),
                    integrand_prefactor: low.clone(),
                    prepared_lu_hosts: vec![],
                    physical_overlaps: None,
                })
                .collect();
            let draw = GammaLoopSample {
                groups: vec![(Some(GroupId(4)), rows)],
            };
            let canonical = draw
                .clone()
                .into_canonical(&externals, DependentMomentaConstructor::CrossSection)
                .unwrap();
            assert_eq!(canonical.groups[0].0, Some(GroupId(4)));
            for (native, row) in draw.groups[0].1.iter().zip(&canonical.groups[0].1) {
                assert_eq!(
                    (native.graph_id, native.channel_id),
                    (row.graph_id, row.channel_id)
                );
                assert_eq!(row.sample.sample.loop_mom_cache_id, 17);
                assert_eq!(row.sample.sample.loop_mom_base_cache_id, 11);
                assert_eq!(row.sample.sample.external_mom_cache_id, 23);
                assert_eq!(row.sample.sample.external_mom_base_cache_id, 19);
                assert_eq!(row.sample.sample.orientation, Some(7));
                assert_eq!(row.sample.sample.parameterization_branch, Some(2));
                assert_eq!(row.sample.external_moms(), &original);
                assert_eq!(
                    row.sample.jacobian(),
                    native.sample.jacobian().to_arb_exact().unwrap()
                );
                assert_eq!(row.sample.loop_moms().0[0].px, low.to_arb_exact().unwrap());
                assert_eq!(row.integrand_prefactor, low.to_arb_exact().unwrap());
            }
            let frozen = format!("{canonical:?}");
            let native = canonical.materialize::<T>(1.0e-13).unwrap();
            assert_eq!(native.groups[0].1[0].sample.loop_moms(), sample.loop_moms());
            let double = canonical.materialize::<f64>(1.0e-6).unwrap();
            let rotated = double.rotate(&Rotation::new(RotationMethod::Pi2Z), 101, 103);
            assert_eq!(rotated.groups[0].1[2].sample.sample.loop_mom_cache_id, 103);
            canonical.materialize::<ArbPrec>(1.0e-13).unwrap();
            assert_eq!(format!("{canonical:?}"), frozen);
            let mut invalid = draw.clone();
            invalid.groups[0].1[0].sample.sample.dual_loop_moms = Some(vec![].into());
            assert!(
                invalid
                    .into_canonical(&externals, DependentMomentaConstructor::CrossSection)
                    .is_err()
            );
            let mut invalid = draw;
            invalid.groups[0].1[0].physical_overlaps =
                Some(std::sync::Arc::new(Default::default()));
            assert!(
                invalid
                    .into_canonical(&externals, DependentMomentaConstructor::CrossSection)
                    .is_err()
            );
        }
        check(F(QuadFloat::from(DoubleFloat::from_compensated_sum(
            1.0,
            2.0_f64.powi(-400),
        ))));
        let one = F::<SamplingFloat>::default().one();
        check(&one + one.from_i64(2).powi(-180) + one.from_i64(2).powi(-240));
    }

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

        // Original source tokens retain their exact binary values. Production
        // maps consume these once at the fixed epoch precision before materialization.
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
    fn canonical_row_materializes_directly_without_losing_hidden_precision() {
        use crate::momentum::ThreeMomentum;
        use crate::utils::{ArbPrec, QuadFloat};
        let coordinate = F::<ArbPrec>::from_f64(1.0) + F::<ArbPrec>::from_f64(2.0).powi(-80);
        let sample = MomentumSample::new(
            vec![ThreeMomentum::new(
                coordinate,
                F(ArbPrec::from_f64_exact_binary(0.0)),
                F(ArbPrec::from_f64_exact_binary(0.0)),
            )]
            .into(),
            7,
            &Externals::default(),
            9,
            F(ArbPrec::from_f64_exact_binary(3.0)),
            DependentMomentaConstructor::CrossSection,
            Some(2),
        )
        .unwrap();
        let draw = GammaLoopSample {
            groups: vec![(
                None,
                vec![DiscreteGraphSample {
                    graph_id: 4,
                    channel_id: Some(SamplingChannelId(11)),
                    prepared_lu_hosts: vec![],
                    physical_overlaps: None,
                    integrand_prefactor: sample.one(),
                    sample,
                }],
            )],
        };
        let double = draw.materialize::<f64>(1.0e-6).unwrap();
        let quad = draw.materialize::<QuadFloat>(1.0e-15).unwrap();
        assert_eq!(double.get_default_sample().loop_moms().0[0].px, F(1.0));
        assert_ne!(
            quad.get_default_sample().loop_moms().0[0].px,
            F(QuadFloat::from_f64_exact_binary(1.0))
        );
        assert_eq!(quad.groups[0].1[0].channel_id, Some(SamplingChannelId(11)));
        assert_eq!(quad.get_default_sample().sample.orientation, Some(2));
        assert_eq!(quad.get_default_sample().sample.loop_mom_cache_id, 7);
        assert_eq!(
            quad.get_default_sample().jacobian(),
            F(QuadFloat::from_f64_exact_binary(3.0))
        );
    }

    #[test]
    fn canonical_row_rejects_an_unrepresentable_combined_factor() {
        use crate::utils::ArbPrec;
        let factor = F::<ArbPrec>::from_f64(1.0) + F::<ArbPrec>::from_f64(2.0).powi(-80);
        let sample = MomentumSample::new(
            LoopMomenta::from(vec![crate::momentum::ThreeMomentum::new(
                F::<ArbPrec>::from_f64(1.0),
                F::<ArbPrec>::from_f64(0.0),
                F::<ArbPrec>::from_f64(0.0),
            )]),
            0,
            &Externals::default(),
            0,
            factor,
            DependentMomentaConstructor::CrossSection,
            None,
        )
        .unwrap();
        let mut draw = GammaLoopSample {
            groups: vec![(
                None,
                vec![DiscreteGraphSample {
                    graph_id: 0,
                    channel_id: None,
                    prepared_lu_hosts: vec![],
                    physical_overlaps: None,
                    integrand_prefactor: sample.one(),
                    sample,
                }],
            )],
        };
        assert!(draw.materialize::<f64>(1.0e-20).is_ok());
        assert!(draw.materialize::<f64>(1.0e-30).is_err());
        draw.groups[0].1[0].sample.sample.jacobian *= F::<ArbPrec>::from_f64(2.0).powi(2048);
        assert!(draw.materialize::<f64>(1.0e-6).is_err());
    }

    #[test]
    fn sampling_map_rejects_non_unit_original_coordinates() {
        use crate::integrands::process::{
            SamplingMapComponent, SurfaceRadialMap, sampling::context::SamplingMapContext,
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
    fn summed_sampling_retains_distinct_completed_rows_under_rotation() {
        use crate::momentum::{Rotation, RotationMethod, ThreeMomentum};
        let sample = MomentumSample::new(
            vec![ThreeMomentum::new(F(1.0), F(2.0), F(3.0))].into(),
            0,
            &Externals::default(),
            0,
            F(2.0),
            DependentMomentaConstructor::CrossSection,
            None,
        )
        .unwrap();
        let mut second = sample.clone();
        second.sample.loop_moms.0[0].px = F(4.0);
        second.sample.jacobian = F(7.0);
        let draw = GammaLoopSample {
            groups: vec![(
                Some(crate::graph::GroupId(0)),
                [sample, second]
                    .into_iter()
                    .enumerate()
                    .map(|(id, sample)| DiscreteGraphSample {
                        graph_id: 0,
                        channel_id: Some(SamplingChannelId(id)),
                        prepared_lu_hosts: vec![],
                        physical_overlaps: None,
                        integrand_prefactor: sample.one(),
                        sample,
                    })
                    .collect(),
            )],
        };
        let rotated = draw.rotate(&Rotation::new(RotationMethod::Pi2Z), 17, 23);
        assert_eq!(rotated.row_count(), 2);
        assert_eq!(rotated.groups[0].1[0].sample.jacobian(), F(2.0));
        assert_eq!(rotated.groups[0].1[1].sample.jacobian(), F(7.0));
        assert_eq!(rotated.groups[0].1[0].sample.sample.loop_mom_cache_id, 17);
        assert_eq!(rotated.groups[0].1[1].sample.sample.loop_mom_cache_id, 18);
        assert_ne!(
            rotated.groups[0].1[0].sample.loop_moms(),
            rotated.groups[0].1[1].sample.loop_moms()
        );
    }
}
