use std::{
    collections::HashSet,
    fs::{self},
    path::Path,
};

use bincode_trait_derive::{Decode, Encode};

use color_eyre::{Help, Result};

use colored::Colorize;
use eyre::{Context, eyre};
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
    subgraph::{ModifySubSet, SubSetLike, subset::SubSet},
};
use momtrop::SampleGenerator;
use rayon::iter::{IntoParallelRefMutIterator, ParallelIterator};

use spenso::algebra::complex::Complex;
use symbolica::{
    atom::AtomCore,
    evaluate::OptimizationSettings,
    numerical_integration::{Grid, Sample},
};
use tracing::{debug, info, instrument, warn};
use typed_index_collections::{TiVec, ti_vec};

use crate::{
    DependentMomentaConstructor, F, FloatLike, GammaLoopContext, GammaLoopContextContainer,
    cff::{
        esurface::{
            EsurfaceCollection, ExistingEsurfaces, GroupEsurfaceId, RaisedEsurfaceId,
            get_representative,
        },
        expression::OrientationID,
        surface::HybridSurfaceID,
    },
    graph::{
        FeynmanGraph, Graph, GraphGroup, GraphGroupPosition, GroupId, LMBext, LmbIndex,
        LoopMomentumBasis,
    },
    integrands::{
        HasIntegrand,
        evaluation::{EvaluationResult, GraphEvaluationResult},
        process::{
            ChannelIndex, LmbChannelWeightingSettings, ParamBuilder,
            evaluators::{ActiveF64Backend, EvaluatorStack},
        },
    },
    model::Model,
    momentum::{
        Helicity, Rotation, RotationMethod, SignOrZero, ThreeMomentum,
        sample::{ExternalIndex, MomentumSample},
        signature::SignatureLike,
    },
    observables::{AdditionalWeightKey, EventProcessingRuntime, GenericEvent},
    processes::{AmplitudeGraph, GraphGenerationStats, GroupDerivedData},
    settings::{GlobalSettings, RuntimeSettings, global::FrozenCompilationMode},
    subtraction::{
        amplitude_counterterm::{
            AmplitudeCountertermAtom, AmplitudeCountertermData, AmplitudeCountertermEvaluation,
            AmplitudeLocalCountertermEvaluation, OverlapStructureWithKinematics,
        },
        overlap::{OverlapInput, SingleGraphOverlapData, find_maximal_overlap},
    },
    utils::{
        ArbPrec, W_, compute_shift_part, serde_utils::SmartSerde, symbolica_ext::LOGPRINTOPTS,
    },
};

use super::{
    GraphTerm, GraphTermEvaluationContext, LmbMultiChannelingSetup, ProcessIntegrandImpl,
    RuntimeCache, create_grid, evaluate_sample, filtered_orientation_count,
    format_lmb_channel_label, format_orientation_label, histogram_process_info_for_integrand,
    prepare_buffered_event, resolve_visible_orientation_id,
};

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeGraphTerm {
    pub original_integrand: EvaluatorStack,
    pub orientations: TiVec<OrientationID, EdgeVec<Orientation>>,
    pub orientation_filter: SubSet<OrientationID>,
    pub esurfaces: EsurfaceCollection,
    pub threshold_counterterm: AmplitudeCountertermData,
    pub multi_channeling_setup: LmbMultiChannelingSetup,
    pub lmbs: TiVec<LmbIndex, LoopMomentumBasis>,
    pub tropical_sampler: Option<SampleGenerator<3>>,
    pub graph: Graph,
    pub estimated_scale: Option<F<f64>>,
    pub param_builder: ParamBuilder,
    pub real_mass_vec: Option<EdgeVec<Option<F<f64>>>>,
    pub master_external_signature: SignatureLike<ExternalIndex>,
    pub master_external_pdgs: Vec<isize>,
}

/// Num(sigma_1,sigma_2,...)*(CFF_1 delta(edge(1),1) delta_(1,1,1,-1,1)+CFF_3 delta_(1,1,1,-1,1)+CFF_2 delta_(1,1,1,-1,1))
impl AmplitudeGraphTerm {
    pub fn kinematics_for_threshold_approach(
        &mut self,
        settings: &RuntimeSettings,
        model: &Model,
        momentum_sample: &MomentumSample<ArbPrec>,
    ) -> Result<OverlapStructureWithKinematics<ArbPrec>> {
        self.threshold_counterterm.kinematics_for_approach(
            momentum_sample,
            &self.graph,
            model,
            &self.esurfaces,
            &Rotation::new(RotationMethod::Identity),
            settings,
        )
    }
    pub fn from_amplitude_graph(
        graph: &AmplitudeGraph,
        own_group_position: GraphGroupPosition,
        esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<RaisedEsurfaceId>>>,
        _model: &Model,
        settings: &GlobalSettings,
    ) -> Result<(Self, GraphGenerationStats)> {
        let started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_graph_term_from_graph_start",
            graph = %graph.graph.name,
            "Generation timing milestone"
        );
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        let mut stats = GraphGenerationStats::default();
        let selected_generation_orientations = graph
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .filter(|orientation| settings.generation.orientation_pattern.filter(*orientation))
            .collect_vec();
        let orientations: TiVec<OrientationID, EdgeVec<Orientation>> =
            selected_generation_orientations
                .iter()
                .map(|a| a.data.orientation.clone())
                .collect();
        crate::debug_tags!(#generation, #profile, #compile, #graph, #orientation, #summary;
            stage = "amplitude_graph_term_orientations_done",
            graph = %graph.graph.name,
            orientation_count = orientations.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        if orientations.is_empty() {
            let pattern = settings
                .generation
                .orientation_pattern
                .pat
                .as_ref()
                .map(ToString::to_string)
                .unwrap_or_else(|| "<empty>".to_string());
            return Err(eyre!(
                "Generation orientation pattern {pattern} matched no orientations for graph {}",
                graph.graph.name
            ));
        }

        let selected_generation_esurfaces = selected_generation_orientations
            .iter()
            .flat_map(|orientation| {
                orientation.expression.iter_nodes().filter_map(|tree_node| {
                    if let HybridSurfaceID::Esurface(esurface_id) = tree_node.data {
                        Some(esurface_id)
                    } else {
                        None
                    }
                })
            })
            .collect::<HashSet<_>>();
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_graph_term_esurfaces_done",
            graph = %graph.graph.name,
            selected_esurface_count = selected_generation_esurfaces.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );

        crate::debug_tags!(#generation, #graph, #orientation, #compile, #dump;
            orientation_parametric_integrand = %graph.derived_data.all_mighty_integrand.printer(LOGPRINTOPTS.clone()),
            "Building evaluator for all orientations \n{}",
            graph.graph.param_builder.table()
        );

        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        let original_started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_graph_term_original_evaluator_start",
            graph = %graph.graph.name,
            orientation_count = orientations.len(),
            "Generation timing milestone"
        );
        let (original_integrand, evaluator_timings) = EvaluatorStack::new_with_timings(
            &[&graph.derived_data.all_mighty_integrand],
            &graph.graph.param_builder,
            orientations.as_slice().as_ref(),
            None,
            &settings.generation.evaluator,
        )?;
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_graph_term_original_evaluator_done",
            graph = %graph.graph.name,
            evaluator_count = original_integrand.generic_evaluator_count(),
            elapsed_ms = original_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            spenso_ms = evaluator_timings.spenso_time.as_secs_f64() * 1000.0,
            symbolica_ms = evaluator_timings.symbolica_time.as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        stats.add_evaluator_build_timings(evaluator_timings);
        stats.evaluator_count += original_integrand.generic_evaluator_count();

        let mut threshold_counterterm = AmplitudeCountertermData::new_empty(own_group_position);
        let mut threshold_evaluators =
            Vec::with_capacity(graph.derived_data.threshold_counterterms.len());
        let selected_generation_raised_esurfaces: HashSet<RaisedEsurfaceId> =
            if graph.derived_data.threshold_counterterms.is_empty() {
                HashSet::new()
            } else {
                selected_generation_esurfaces
                    .iter()
                    .map(|esurface_id| graph.derived_data.raised_esurface_ids[*esurface_id])
                    .collect()
            };
        let active_mask: TiVec<RaisedEsurfaceId, bool> = graph
            .derived_data
            .threshold_counterterms
            .iter_enumerated()
            .map(|(raised_esurface_id, _)| {
                selected_generation_raised_esurfaces.contains(&raised_esurface_id)
            })
            .collect();
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_graph_term_threshold_setup_done",
            graph = %graph.graph.name,
            threshold_count = graph.derived_data.threshold_counterterms.len(),
            active_threshold_count = active_mask.iter().filter(|active| **active).count(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        for (raised_esurface_id, ct) in graph.derived_data.threshold_counterterms.iter_enumerated()
        {
            if crate::is_interrupted() {
                return Err(eyre!("Generation interrupted by user"));
            }
            let threshold_started = std::time::Instant::now();
            crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
                stage = "amplitude_graph_term_threshold_evaluator_start",
                graph = %graph.graph.name,
                raised_esurface_id = %raised_esurface_id.0,
                active = active_mask[raised_esurface_id],
                "Generation timing milestone"
            );
            let masked_counterterm = if active_mask[raised_esurface_id] {
                ct.clone()
            } else {
                ct.zero_like()
            };
            let (evaluator, evaluator_timings) = masked_counterterm.to_evaluator_with_timings(
                &graph.graph.param_builder,
                &orientations,
                settings,
            );
            crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
                stage = "amplitude_graph_term_threshold_evaluator_done",
                graph = %graph.graph.name,
                raised_esurface_id = %raised_esurface_id.0,
                active = active_mask[raised_esurface_id],
                evaluator_count = evaluator.generic_evaluator_count(),
                elapsed_ms = threshold_started.elapsed().as_secs_f64() * 1000.0,
                total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
                spenso_ms = evaluator_timings.spenso_time.as_secs_f64() * 1000.0,
                symbolica_ms = evaluator_timings.symbolica_time.as_secs_f64() * 1000.0,
                "Generation timing milestone"
            );
            if crate::is_interrupted() {
                return Err(eyre!("Generation interrupted by user"));
            }
            stats.add_evaluator_build_timings(evaluator_timings);
            stats.evaluator_count += evaluator.generic_evaluator_count();
            threshold_evaluators.push(evaluator);
        }
        threshold_counterterm.evaluators = threshold_evaluators.into();
        threshold_counterterm.generated_mask = graph
            .derived_data
            .threshold_counterterms
            .iter()
            .map(AmplitudeCountertermAtom::is_generated)
            .collect();
        // `generated_mask` answers "does this threshold slot exist symbolically at all?",
        // while `active_mask` answers "is that generated slot compatible with the
        // selected generation-time orientation subset for this evaluator build?"
        threshold_counterterm.active_mask = active_mask;
        threshold_counterterm.raised_data = graph.derived_data.raised_data.clone();
        threshold_counterterm.helper_evaluators = graph
            .derived_data
            .raised_data
            .pass_two_evaluator
            .clone()
            .unwrap_or_default();
        stats.evaluator_count += threshold_counterterm.helper_evaluators.len();

        threshold_counterterm.local_esurface_exists = ti_vec![true; esurface_map.len()];
        threshold_counterterm.esurface_map = esurface_map;
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_graph_term_from_graph_done",
            graph = %graph.graph.name,
            evaluator_count = stats.evaluator_count,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            spenso_ms = stats.evaluator_spenso_time.as_secs_f64() * 1000.0,
            symbolica_ms = stats.evaluator_symbolica_time.as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );

        Ok((
            AmplitudeGraphTerm {
                orientation_filter: SubSet::full(orientations.len()),
                orientations,
                original_integrand,
                tropical_sampler: graph.derived_data.tropical_sampler.clone(),
                graph: graph.graph.clone(),
                multi_channeling_setup: LmbMultiChannelingSetup {
                    channels: TiVec::new(),
                    graph: graph.graph.clone(), // will be overwritten later,
                    all_bases: TiVec::new(),
                }, // to be taken from froup master
                lmbs: graph
                    .derived_data
                    .lmbs
                    .clone()
                    .expect("lmbs should have been created"),
                threshold_counterterm,
                estimated_scale: None,
                esurfaces: graph
                    .derived_data
                    .cff_expression
                    .as_ref()
                    .expect("cff_expression should have been created")
                    .surfaces
                    .esurface_cache
                    .clone(),
                param_builder: graph.graph.param_builder.clone(),
                real_mass_vec: None,
                master_external_signature: graph.graph.get_external_signature(),
                master_external_pdgs: graph
                    .graph
                    .get_external_partcles()
                    .into_iter()
                    .map(|particle| particle.pdg_code)
                    .collect(),
            },
            stats,
        ))
    }

    #[instrument(
          name = "compile",
          level = "info",
          skip(self, path, override_existing, frozen_mode),
          fields(
              graph.name = %self.graph.name,
              path = %path.as_ref().display(),
          )
      )]
    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        frozen_mode: &FrozenCompilationMode,
    ) -> Result<std::time::Duration> {
        let compile_started = std::time::Instant::now();
        let graph_path = path.as_ref().join(&self.graph.name);

        fs::create_dir_all(&graph_path).with_context(|| {
            format!(
                "Trying to create directory to save amplitude {}",
                graph_path.display()
            )
        })?;

        self.original_integrand.compile(
            "orientation_parametric_integrand",
            &graph_path,
            frozen_mode,
        )?;

        self.threshold_counterterm
            .compile(&graph_path, override_existing, frozen_mode)?;

        Ok(compile_started.elapsed())
    }

    pub(crate) fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut crate::integrands::process::GenericEvaluator) -> Result<()>,
    ) -> Result<()> {
        self.original_integrand
            .for_each_generic_evaluator_mut(&mut f)?;
        self.threshold_counterterm
            .for_each_generic_evaluator_mut(&mut f)?;
        Ok(())
    }

    pub(crate) fn generic_evaluator_count(&self) -> usize {
        self.original_integrand.generic_evaluator_count()
            + self.threshold_counterterm.generic_evaluator_count()
    }

    #[instrument(
          level = "debug",
          skip_all,
          fields(
              term.name = %self.name(),
          )
    )]
    fn generate_event<T: FloatLike>(
        &self,
        settings: &RuntimeSettings,
        orientation_id: Option<usize>,
        channel_id: Option<ChannelIndex>,
    ) -> Result<GenericEvent<T>> {
        let externals = settings
            .kinematics
            .externals
            .get_dependent_externals(DependentMomentaConstructor::Amplitude(
                &self.master_external_signature,
            ))
            .with_context(|| {
                format!(
                    "when getting master externals to build amplitude event for graph: {}",
                    self.graph.name
                )
            })?;

        if externals.len() != self.master_external_pdgs.len() {
            return Err(eyre!(
                "Amplitude graph '{}' has inconsistent master external metadata: {} momenta vs {} PDGs.",
                self.graph.name,
                externals.len(),
                self.master_external_pdgs.len()
            ));
        }

        let mut event = GenericEvent::default();
        event.cut_info.cut_id = 0;
        event.cut_info.orientation_id = orientation_id;
        event.cut_info.lmb_channel_id = channel_id.map(usize::from);
        event.cut_info.lmb_channel_edge_ids =
            channel_id.map(|channel_id| self.multi_channeling_setup.channel_edge_ids(channel_id));

        for ((sign, momentum), pdg) in self
            .master_external_signature
            .iter()
            .zip(externals.iter())
            .zip(self.master_external_pdgs.iter().copied())
        {
            match sign {
                SignOrZero::Plus => {
                    event.kinematic_configuration.0.push(momentum.clone());
                    event.cut_info.particle_pdgs.0.push(pdg);
                }
                SignOrZero::Minus => {
                    event.kinematic_configuration.1.push(momentum.clone());
                    event.cut_info.particle_pdgs.1.push(pdg);
                }
                SignOrZero::Zero => {
                    return Err(eyre!(
                        "Amplitude graph '{}' has an invalid zero-sign external momentum in its master signature.",
                        self.graph.name
                    ));
                }
            }
        }

        Ok(event)
    }

    fn evaluate_impl<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        context: &mut GraphTermEvaluationContext<'_, '_, T>,
    ) -> Result<(Complex<F<T>>, AmplitudeCountertermEvaluation<T>)> {
        let (momentum_sample, prefactor) =
            if let Some((channel_id, alpha, channel_weight)) = &context.channel_id {
                let parameterization_settings = context
                    .settings
                    .sampling
                    .get_parameterization_settings()
                    .expect("LMB multichanneling requires a parameterization.");
                let weighting_settings = LmbChannelWeightingSettings {
                    model: context.model,
                    alpha,
                    channel_weight: *channel_weight,
                    parameterization_settings: &parameterization_settings,
                    e_cm: context.settings.kinematics.e_cm,
                };

                self.multi_channeling_setup
                    .reinterpret_loop_momenta_and_compute_prefactor(
                        *channel_id,
                        momentum_sample,
                        0,
                        weighting_settings,
                    )
            } else {
                (momentum_sample.clone(), momentum_sample.one())
            };

        let hel = context.settings.kinematics.externals.get_helicities();
        let orientations =
            momentum_sample.orientations(&self.orientation_filter, &self.orientations);

        debug!("loop_moms: {}", momentum_sample.loop_moms());
        debug!("jacobian: {:16e}", momentum_sample.jacobian());
        // debug!("Og paramBuilder: \n{}", self.param_builder.table());

        let input = T::get_parameters(
            &mut self.param_builder,
            (
                context.settings.general.enable_cache,
                context.settings.general.debug_cache,
            ),
            &self.graph,
            &momentum_sample,
            hel,
            &context.settings.additional_params(),
            None,
            None,
            None,
        );
        let result = self
            .original_integrand
            .evaluate(
                input,
                orientations,
                context.settings,
                context.evaluation_metadata,
                context.record_primary_timing,
            )?
            .pop()
            .unwrap()
            .unwrap_real();
        // debug!("parambuilder 244: {}", self.param_builder);
        let counterterm_evaluation = self.threshold_counterterm.evaluate(
            &momentum_sample,
            &self.graph,
            context.model,
            &self.esurfaces,
            context.rotation,
            context.settings,
            &mut self.param_builder,
            orientations,
            context.evaluation_metadata,
            context.record_primary_timing,
        )?;
        let sum_of_cts = counterterm_evaluation.total.clone();

        crate::debug_tags!(#integration, #subtraction;
            bare_cff = format!("{result:16e}"),
            "{}: {result:16e}",
            self.graph.name
        );
        crate::debug_tags!(#integration, #subtraction;
            cts = format!("{sum_of_cts:16e}"),
            "{}",
            self.graph.name
        );
        crate::debug_tags!(#integration, #subtraction; "result: {result:16e}");
        crate::debug_tags!(#integration, #subtraction; "sum_of_cts: {sum_of_cts:16e}");

        crate::debug_tags!(#integration, #subtraction;
            value = format!("{sum_of_cts:16e}"),
            "evaluated sum of threshold counterterms"
        );

        let diff = result.clone() - sum_of_cts.clone();
        let result_is_nonfinite = result.re.is_nan()
            || result.re.is_infinite()
            || result.im.is_nan()
            || result.im.is_infinite();
        let sum_of_cts_is_nonfinite = sum_of_cts.re.is_nan()
            || sum_of_cts.re.is_infinite()
            || sum_of_cts.im.is_nan()
            || sum_of_cts.im.is_infinite();
        let diff_is_nonfinite =
            diff.re.is_nan() || diff.re.is_infinite() || diff.im.is_nan() || diff.im.is_infinite();
        crate::debug_tags!(#integration, #subtraction, #threshold, #inspect;
            stage = "amplitude_threshold_subtraction",
            graph = %self.graph.name,
            original = %format!("{:+16e}", result),
            cts = %format!("{:+16e}", sum_of_cts),
            diff = %format!("{:+16e}", diff),
            prefactor = %format!("{:+16e}", prefactor),
            original_nonfinite = result_is_nonfinite,
            cts_nonfinite = sum_of_cts_is_nonfinite,
            diff_nonfinite = diff_is_nonfinite,
            "amplitude threshold subtraction"
        );

        Ok((
            diff * prefactor.clone(),
            AmplitudeCountertermEvaluation {
                total: sum_of_cts * prefactor.clone(),
                local_counterterms: counterterm_evaluation
                    .local_counterterms
                    .into_iter()
                    .map(
                        |AmplitudeLocalCountertermEvaluation {
                             esurface_id,
                             overlap_group,
                             value,
                         }| AmplitudeLocalCountertermEvaluation {
                            esurface_id,
                            overlap_group,
                            value: value * prefactor.clone(),
                        },
                    )
                    .collect(),
            },
        ))
    }
}

impl GraphTerm for AmplitudeGraphTerm {
    #[instrument(
          skip_all,
          fields(
              term.name = %self.name(),
          ),
          err
    )]
    fn warm_up(&mut self, settings: &RuntimeSettings, model: &Model) -> Result<()> {
        self.orientation_filter = SubSet::empty(self.orientations.len());
        for (id, o) in self.orientations.iter_enumerated() {
            if settings.general.orientation_pat.filter(o) {
                self.orientation_filter.add(id);
            }
        }
        if self.orientation_filter.included_iter().next().is_none() {
            let pattern = settings
                .general
                .orientation_pat
                .pat
                .as_ref()
                .map(ToString::to_string)
                .unwrap_or_else(|| "<empty>".to_string());
            return Err(eyre!(
                "Runtime orientation pattern {pattern} matched no orientations for graph {}",
                self.graph.name
            ));
        }

        self.estimated_scale = Some(
            self.graph
                .expected_scale(F(settings.kinematics.e_cm), model),
        );

        let externals = settings
            .kinematics
            .externals
            .get_dependent_externals(DependentMomentaConstructor::Amplitude(
                &self.graph.get_external_signature(),
            ))
            .with_context(|| {
                format!("when getting externals to build amplitude graph term for integrand for graph: {}", self.graph.name)
            })?;

        self.graph
            .param_builder
            .add_external_four_mom_all_derivatives(&externals);
        let pols = self.graph.param_builder.pairs.polarizations_values(
            &self.graph,
            &externals,
            settings.kinematics.externals.get_helicities(),
        );
        self.graph.param_builder.pairs.warn_zero_polarizations(
            &self.graph,
            &externals,
            settings.kinematics.externals.get_helicities(),
        );

        for (value_index, values) in self.graph.param_builder.values.iter_mut().enumerate() {
            let multiplicative_offset = value_index + 1;
            let mut polarization_start = self
                .graph
                .param_builder
                .pairs
                .polarizations
                .value_range
                .start
                * multiplicative_offset;

            for pol_value in pols.iter() {
                values[polarization_start] = *pol_value;
                polarization_start += multiplicative_offset;
            }
        }

        self.graph
            .param_builder
            .m_uv_value(Complex::new_re(F(settings.general.m_uv)));

        self.graph
            .param_builder
            .renormalization_localization_scale_value(Complex::new_re(F(settings
                .general
                .renormalization_localization_scale)));
        self.graph
            .param_builder
            .mu_r_sq_value(Complex::new_re(F(settings.general.mu_r_sq())));
        self.graph.param_builder.update_model_values(model);

        self.param_builder = self.graph.param_builder.clone();

        let masses = self
            .graph
            .new_edgevec(|e, _, _| e.mass_value(model, &self.param_builder).map(|c| c.re));

        self.real_mass_vec = Some(masses);

        Ok(())
    }

    fn name(&self) -> String {
        self.graph.name.clone()
    }

    fn orientation_label(&self, orientation_id: usize) -> Option<String> {
        self.orientations
            .get(resolve_visible_orientation_id(
                &self.orientation_filter,
                orientation_id,
            )?)
            .map(format_orientation_label)
    }

    fn lmb_channel_label(&self, channel_id: ChannelIndex) -> Option<String> {
        Some(format_lmb_channel_label(
            &self.multi_channeling_setup.channel_edge_ids(channel_id),
        ))
    }

    fn get_graph(&self) -> &Graph {
        &self.graph
    }

    fn get_num_channels(&self) -> usize {
        self.multi_channeling_setup.channels.len()
    }

    fn evaluate<T: FloatLike>(
        &mut self,
        momentum_sample: &MomentumSample<T>,
        mut context: GraphTermEvaluationContext<'_, '_, T>,
    ) -> Result<GraphEvaluationResult<T>> {
        let event_channel_id = context
            .channel_id
            .as_ref()
            .map(|(channel_id, _, _)| *channel_id);
        let prepared_event = prepare_buffered_event(
            context.settings,
            context.rotation,
            context.event_processing_runtime.take(),
            || {
                self.generate_event(
                    context.settings,
                    momentum_sample.sample.orientation,
                    event_channel_id,
                )
            },
        )?;
        if !prepared_event.selectors_pass {
            return Ok(GraphEvaluationResult {
                integrand_result: Complex::new_re(momentum_sample.zero()),
                event_groups: crate::observables::GenericEventGroupList::default(),
                event_processing_time: prepared_event.event_processing_time,
                generated_event_count: prepared_event.generated_event_count,
                accepted_event_count: prepared_event.accepted_event_count,
            });
        }

        let (integrand_result, counterterm_evaluation) =
            self.evaluate_impl(momentum_sample, &mut context)?;

        let mut event_groups = crate::observables::GenericEventGroupList::default();
        let generated_event_count = prepared_event.generated_event_count;
        let accepted_event_count = prepared_event.accepted_event_count;
        if let Some(mut event) = prepared_event.buffered_event {
            event.weight = integrand_result.clone();

            if context.settings.general.store_additional_weights_in_event {
                let original = integrand_result.clone() + counterterm_evaluation.total;
                event
                    .additional_weights
                    .weights
                    .insert(AdditionalWeightKey::Original, original);
                for threshold_counterterm in counterterm_evaluation.local_counterterms.into_iter() {
                    event.additional_weights.weights.insert(
                        AdditionalWeightKey::AmplitudeThresholdCounterterm {
                            esurface_id: threshold_counterterm.esurface_id.0,
                            overlap_group: threshold_counterterm.overlap_group,
                        },
                        -threshold_counterterm.value,
                    );
                }
            }

            event_groups.push_singleton(event);
        }

        Ok(GraphEvaluationResult {
            integrand_result,
            event_groups,
            event_processing_time: prepared_event.event_processing_time,
            generated_event_count,
            accepted_event_count,
        })
    }

    fn get_num_orientations(&self) -> usize {
        filtered_orientation_count(&self.orientation_filter, &self.orientations)
    }

    fn get_tropical_sampler(&self) -> &SampleGenerator<3> {
        self.tropical_sampler
            .as_ref()
            .expect("Tropical sampler should be set.")
    }

    fn get_mut_param_builder(&mut self) -> &mut ParamBuilder<f64> {
        &mut self.param_builder
    }

    fn get_real_mass_vector(&self) -> EdgeVec<Option<F<f64>>> {
        self.real_mass_vec
            .as_ref()
            .expect("real mass vector should be set")
            .clone()
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeIntegrand {
    pub settings: RuntimeSettings,
    pub data: AmplitudeIntegrandData,
    pub(crate) event_processing_runtime: RuntimeCache<EventProcessingRuntime>,
    pub(crate) active_f64_backend: RuntimeCache<ActiveF64Backend>,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeIntegrandData {
    pub rotations: Option<Vec<Rotation>>,
    pub name: String,
    pub compilation: FrozenCompilationMode,
    pub loop_cache_id: usize,
    pub external_cache_id: usize,
    /// Cache ID for the base (unrotated) external momentum configuration
    pub base_external_cache_id: usize,
    pub graph_terms: Vec<AmplitudeGraphTerm>,
    pub external_signature: SignatureLike<ExternalIndex>,
    pub graph_group_structure: TiVec<GroupId, GraphGroup>,
    pub graph_to_group_id: Vec<usize>,
    pub group_derived_data: TiVec<GroupId, GroupDerivedData>,
}

pub mod export;
pub mod load;

impl AmplitudeIntegrand {
    fn threshold_esurface_specifier(
        &self,
        group_id: GroupId,
        group_esurface_id: GroupEsurfaceId,
    ) -> String {
        let details = self.data.group_derived_data[group_id].esurface_map[group_esurface_id]
            .iter_enumerated()
            .filter_map(|(graph_group_pos, raised_esurface_id)| {
                raised_esurface_id.map(|raised_esurface_id| (graph_group_pos, raised_esurface_id))
            })
            .find(|(graph_group_pos, _)| {
                let graph_id = self.data.graph_group_structure[group_id][*graph_group_pos];
                self.data.graph_terms[graph_id]
                    .threshold_counterterm
                    .local_esurface_exists[group_esurface_id]
            })
            .or_else(|| {
                get_representative(
                    &self.data.group_derived_data[group_id].esurface_map[group_esurface_id],
                )
                .ok()
            })
            .map(|(graph_group_pos, raised_esurface_id)| {
                let graph_id = self.data.graph_group_structure[group_id][graph_group_pos];
                let graph_term = &self.data.graph_terms[graph_id];
                let esurface_id = graph_term.threshold_counterterm.raised_data.raised_groups
                    [raised_esurface_id]
                    .esurface_ids[0];
                let esurface = &graph_term.esurfaces[esurface_id];
                let loop_order = esurface.energies.len().saturating_sub(1);
                let edge_ids = esurface.energies.iter().map(|edge_id| edge_id.0).join(",");

                (loop_order, format!("{}({edge_ids})", group_esurface_id.0))
            });

        match details {
            Some((1, specifier)) => specifier.green().bold().to_string(),
            Some((2, specifier)) => specifier.blue().bold().to_string(),
            Some((_, specifier)) => specifier.cyan().bold().to_string(),
            None => format!("{}(?)", group_esurface_id.0).dimmed().to_string(),
        }
    }

    fn format_overlap_structure(
        &self,
        group_id: GroupId,
        overlap: &crate::subtraction::overlap::OverlapStructure,
    ) -> String {
        overlap
            .overlap_groups
            .iter()
            .map(|overlap_group| {
                let esurfaces = overlap_group
                    .existing_esurfaces
                    .iter()
                    .map(|existing_esurface_id| {
                        let group_esurface_id = overlap.existing_esurfaces[*existing_esurface_id];
                        self.threshold_esurface_specifier(group_id, group_esurface_id)
                    })
                    .join(", ");

                format!(
                    "{}: [{}]",
                    overlap_group.existing_esurfaces.len(),
                    esurfaces
                )
            })
            .join(", ")
    }

    pub(crate) fn kinematics_for_threshold_approach(
        &mut self,
        momentum_sample: &MomentumSample<ArbPrec>,
        model: &Model,
    ) -> Result<Vec<OverlapStructureWithKinematics<ArbPrec>>> {
        self.data
            .graph_terms
            .iter_mut()
            .map(|term| {
                term.kinematics_for_threshold_approach(&self.settings, model, momentum_sample)
            })
            .try_collect()
    }

    fn warn_on_off_shell_external_states(&self, model: &Model) -> Result<()> {
        let externals = self
            .settings
            .kinematics
            .externals
            .get_dependent_externals::<f64>(DependentMomentaConstructor::Amplitude(
                &self.data.external_signature,
            ))?;
        let masses = self.data.graph_terms[0].graph.get_external_masses(model);
        let particles = self.data.graph_terms[0].graph.get_external_partcles();
        let helicities = self.settings.kinematics.externals.get_helicities();

        if externals.len() != masses.len()
            || externals.len() != particles.len()
            || externals.len() != helicities.len()
        {
            return Ok(());
        }

        let one = F::<f64>::from_f64(1.0);
        let threshold = F::<f64>::from_f64(1.0e-8)
            * (F::<f64>::from_f64(self.settings.kinematics.e_cm).square() + one);

        for (external_index, (((momentum, mass), particle), helicity)) in externals
            .iter()
            .zip(masses.iter())
            .zip(particles.iter())
            .zip(helicities.iter())
            .enumerate()
        {
            if particle.is_scalar() || !matches!(helicity, Helicity::Signed(_)) {
                continue;
            }

            let expected_mass_sq = mass.square();
            let actual_mass_sq = momentum.square();
            let mass_sq_difference = actual_mass_sq - expected_mass_sq;
            let off_shellness = if mass_sq_difference < mass_sq_difference.zero() {
                -mass_sq_difference
            } else {
                mass_sq_difference
            };

            if off_shellness > threshold {
                warn!(
                    "External state {external_index} ('{}', PDG {}, helicity {}) is off shell: p^2 = {:+e}, expected m^2 = {:+e} (|Δ| = {:+e}). Fixed-helicity spinor/vector/tensor external states may be ill-defined at this kinematic point.",
                    particle.name,
                    particle.pdg_code,
                    helicity,
                    actual_mass_sq,
                    expected_mass_sq,
                    off_shellness
                );
            }
        }

        Ok(())
    }

    pub(crate) fn frozen_compilation(&self) -> &FrozenCompilationMode {
        &self.data.compilation
    }

    pub(crate) fn active_f64_backend(&self) -> ActiveF64Backend {
        self.active_f64_backend
            .as_ref()
            .copied()
            .unwrap_or(ActiveF64Backend::Eager)
    }

    fn for_each_generic_evaluator_mut(
        &mut self,
        mut f: impl FnMut(&mut crate::integrands::process::GenericEvaluator) -> Result<()>,
    ) -> Result<()> {
        for graph_term in &mut self.data.graph_terms {
            graph_term.for_each_generic_evaluator_mut(&mut f)?;
        }
        Ok(())
    }

    fn has_complete_external_artifacts(&mut self) -> Result<bool> {
        let mut has_all = true;
        self.for_each_generic_evaluator_mut(|evaluator| {
            has_all &= evaluator.has_external_compiled_artifact();
            Ok(())
        })?;
        Ok(has_all)
    }

    pub(crate) fn prepare_runtime_backends_after_generation_with_compile_times(
        &mut self,
    ) -> Result<Vec<std::time::Duration>> {
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        match self.data.compilation {
            FrozenCompilationMode::Symjit => {
                let mut compile_times = Vec::with_capacity(self.data.graph_terms.len());
                for graph_term in &mut self.data.graph_terms {
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    let compile_started = std::time::Instant::now();
                    graph_term
                        .for_each_generic_evaluator_mut(|evaluator| evaluator.activate_symjit())?;
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    compile_times.push(compile_started.elapsed());
                }
                self.active_f64_backend.set(ActiveF64Backend::Symjit);
                Ok(compile_times)
            }
            FrozenCompilationMode::Eager
            | FrozenCompilationMode::Cpp(_)
            | FrozenCompilationMode::Assembly(_) => {
                self.for_each_generic_evaluator_mut(|evaluator| {
                    evaluator.activate_eager();
                    Ok(())
                })?;
                self.active_f64_backend.set(ActiveF64Backend::Eager);
                Ok(vec![std::time::Duration::ZERO; self.data.graph_terms.len()])
            }
        }
    }

    pub(crate) fn prepare_runtime_backends_after_generation(&mut self) -> Result<()> {
        let _ = self.prepare_runtime_backends_after_generation_with_compile_times()?;
        Ok(())
    }

    pub(crate) fn activate_runtime_backends_after_load(
        &mut self,
        allow_symjit_fallback: bool,
    ) -> Result<Option<String>> {
        match self.data.compilation.clone() {
            FrozenCompilationMode::Eager => {
                self.prepare_runtime_backends_after_generation()?;
                Ok(None)
            }
            FrozenCompilationMode::Symjit => {
                self.for_each_generic_evaluator_mut(|evaluator| evaluator.activate_symjit())?;
                self.active_f64_backend.set(ActiveF64Backend::Symjit);
                Ok(None)
            }
            FrozenCompilationMode::Cpp(_) => {
                self.activate_external_after_load(ActiveF64Backend::Cpp, allow_symjit_fallback)
            }
            FrozenCompilationMode::Assembly(_) => {
                self.activate_external_after_load(ActiveF64Backend::Assembly, allow_symjit_fallback)
            }
        }
    }

    fn activate_external_after_load(
        &mut self,
        backend: ActiveF64Backend,
        allow_symjit_fallback: bool,
    ) -> Result<Option<String>> {
        if !self.has_complete_external_artifacts()? {
            self.prepare_runtime_backends_after_generation()?;
            return Ok(None);
        }

        match self.for_each_generic_evaluator_mut(|evaluator| {
            evaluator.activate_external_from_artifact(backend)
        }) {
            Ok(()) => {
                self.active_f64_backend.set(backend);
                Ok(None)
            }
            Err(err) if allow_symjit_fallback => {
                let error_message = err.to_string();
                self.for_each_generic_evaluator_mut(|evaluator| evaluator.activate_symjit())?;
                self.active_f64_backend.set(ActiveF64Backend::Symjit);
                Ok(Some(error_message))
            }
            Err(err) => Err(err),
        }
    }

    pub(crate) fn compile(
        &mut self,
        path: impl AsRef<Path> + Sync,
        override_existing: bool,
        thread_pool: &rayon::ThreadPool,
    ) -> Result<Vec<(String, std::time::Duration)>> {
        let frozen_mode = self.data.compilation.clone();
        let compile_times = thread_pool.install(|| {
            self.data
                .graph_terms
                .par_iter_mut()
                .map(|a| {
                    a.compile(path.as_ref(), override_existing, &frozen_mode)
                        .map(|duration| (a.graph.name.clone(), duration))
                })
                .collect::<Result<Vec<_>>>()
        })?;

        self.active_f64_backend
            .set(ActiveF64Backend::from_frozen_mode(&self.data.compilation));
        Ok(compile_times)
    }

    pub(crate) fn save(&self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let binary = bincode::encode_to_vec(&self.data, bincode::config::standard())?;
        fs::write(path.as_ref().join("integrand.bin"), binary)?;

        // debug!("HE3");
        //
        self.settings
            .to_file(path.as_ref().join("settings.toml"), override_existing)
            .with_context(|| "Error saving settings.toml file for amplitude integrand")?;
        // debug!("HE");

        Ok(())
    }

    pub(crate) fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("integrand.bin"))?;
        let (data, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;

        let settings = SmartSerde::from_file(
            path.as_ref().join("settings.toml"),
            "runtime settings for amplitude integrand",
        )?;

        Ok(AmplitudeIntegrand {
            settings,
            data,
            event_processing_runtime: RuntimeCache::default(),
            active_f64_backend: RuntimeCache::default(),
        })
    }

    pub(crate) fn invalidate_event_processing_runtime(&mut self) {
        self.event_processing_runtime.invalidate();
    }

    pub(crate) fn get_existing_esurfaces(
        &mut self,
        model: &Model,
    ) -> TiVec<GroupId, ExistingEsurfaces> {
        let external_moms = self
            .settings
            .kinematics
            .externals
            .get_dependent_externals::<f64>(DependentMomentaConstructor::Amplitude(
                &self.data.external_signature,
            ))
            .expect("could not get externals");
        let e_cm = self.settings.kinematics.e_cm;

        let mut all_existing_esurfaces = TiVec::new();

        for (group_id, graph_group) in self.data.graph_group_structure.clone().iter_enumerated() {
            let group_esurface_map = self.data.group_derived_data[group_id].esurface_map.clone();

            for (_, graph_id) in graph_group.iter_enumerated() {
                self.data.graph_terms[graph_id]
                    .threshold_counterterm
                    .local_esurface_exists = ti_vec![false; group_esurface_map.len()];
            }

            let mut group_existing_esurfaces = ExistingEsurfaces::new();

            for (group_esurface_id, raised_esurface_map) in group_esurface_map.iter_enumerated() {
                let mapped_esurfaces = raised_esurface_map
                    .iter_enumerated()
                    .filter_map(|(graph_group_pos, option_raised_esurface_id)| {
                        option_raised_esurface_id
                            .map(|raised_esurface_id| (graph_group_pos, raised_esurface_id))
                    })
                    .collect_vec();

                let (representative_graph_group_pos, representative_raised_esurface_id) =
                    mapped_esurfaces
                        .first()
                        .copied()
                        .expect("no graph in group has this esurface, map corrupted");

                let mut representative_exists = false;
                let mut any_candidate_exists = false;

                for (graph_group_pos, raised_esurface_id) in mapped_esurfaces {
                    let graph_id = graph_group[graph_group_pos];
                    let candidate_exists = {
                        let graph_term = &self.data.graph_terms[graph_id];
                        let graph = &graph_term.graph;
                        let raised_group =
                            &graph_term.threshold_counterterm.raised_data.raised_groups
                                [raised_esurface_id];
                        let esurface_id = raised_group.esurface_ids[0];
                        let esurface = &graph_term.esurfaces[esurface_id];
                        let lmb = &graph.loop_momentum_basis;
                        let real_mass_vector = graph.get_real_mass_vector(model);
                        let candidate_exists =
                            esurface.exists(&external_moms, lmb, &real_mass_vector, &F(e_cm));
                        if tracing::event_enabled!(tracing::Level::DEBUG) {
                            let shift_part =
                                esurface.compute_shift_part_from_momenta(&external_moms, lmb);
                            let mass_sum: F<f64> = esurface
                                .energies
                                .iter()
                                .map(|index| &real_mass_vector[*index])
                                .fold(F::from_f64(0.0), |acc, x| acc + x);
                            let zero_vector = ThreeMomentum::new(F(0.0), F(0.0), F(0.0));
                            let shift_vector = esurface
                                .external_shift
                                .iter()
                                .map(|(index, sign)| {
                                    let external_signature = &lmb.edge_signatures[*index].external;
                                    compute_shift_part(external_signature, &external_moms).spatial
                                        * F::from_f64(*sign as f64)
                                })
                                .reduce(|acc, x| acc + x)
                                .unwrap_or_else(|| zero_vector.clone());
                            let shift_vector_sq = shift_vector.norm_squared();
                            let existence_margin = &shift_part * &shift_part
                                - &shift_vector_sq
                                - &mass_sum * &mass_sum;
                            let lmb_reps = graph.integrand_replacement(
                                &graph.full_filter(),
                                &graph.loop_momentum_basis,
                                &[W_.x___],
                            );
                            let atom = esurface.lmb_atom_simplified(graph, &lmb_reps);
                            let raw_atom = esurface.to_atom(&[]);
                            let edge_ids = esurface
                                .energies
                                .iter()
                                .map(|edge_id| edge_id.0)
                                .collect_vec();
                            let local_esurface_ids = raised_group
                                .esurface_ids
                                .iter()
                                .map(|esurface_id| esurface_id.0)
                                .collect_vec();
                            let is_representative = graph_group_pos
                                == representative_graph_group_pos
                                && raised_esurface_id == representative_raised_esurface_id;
                            let generated = graph_term
                                .threshold_counterterm
                                .generated_mask
                                .get(raised_esurface_id)
                                .copied();
                            let active = graph_term
                                .threshold_counterterm
                                .active_mask
                                .get(raised_esurface_id)
                                .copied();

                            crate::debug_tags!(#integration, #subtraction, #threshold, #inspect, #esurface;
                                stage = "amplitude_threshold_esurface_candidate",
                                group_id = group_id.0,
                                group_esurface_id = group_esurface_id.0,
                                graph = %graph.name,
                                graph_group_pos = graph_group_pos.0,
                                raised_esurface_id = raised_esurface_id.0,
                                esurface_id = esurface_id.0,
                                local_esurface_ids = ?local_esurface_ids,
                                edges = ?edge_ids,
                                representative = is_representative,
                                representative_graph_group_pos = representative_graph_group_pos.0,
                                representative_raised_esurface_id = representative_raised_esurface_id.0,
                                candidate_exists,
                                generated = ?generated,
                                active = ?active,
                                max_occurrence = raised_group.max_occurence,
                                shift_part = %format!("{:+16e}", shift_part),
                                shift_vector_sq = %format!("{:+16e}", shift_vector_sq),
                                mass_sum = %format!("{:+16e}", mass_sum),
                                existence_margin = %format!("{:+16e}", existence_margin),
                                file.atom = %atom,
                                file.raw_atom = %raw_atom,
                                "amplitude threshold esurface candidate"
                            );
                        }

                        candidate_exists
                    };

                    self.data.graph_terms[graph_id]
                        .threshold_counterterm
                        .local_esurface_exists[group_esurface_id] = candidate_exists;

                    if graph_group_pos == representative_graph_group_pos
                        && raised_esurface_id == representative_raised_esurface_id
                    {
                        representative_exists = candidate_exists;
                    }
                    any_candidate_exists |= candidate_exists;
                }

                crate::debug_tags!(#integration, #subtraction, #threshold, #inspect, #esurface;
                    stage = "amplitude_threshold_group_esurface_existence",
                    group_id = group_id.0,
                    group_esurface_id = group_esurface_id.0,
                    representative_exists,
                    any_candidate_exists,
                    "amplitude threshold group esurface existence"
                );

                if any_candidate_exists {
                    group_existing_esurfaces.push(group_esurface_id);
                }
            }

            all_existing_esurfaces.push(group_existing_esurfaces);
        }

        let groups_above_threshold = all_existing_esurfaces
            .iter()
            .filter(|existing_esurfaces| !existing_esurfaces.is_empty())
            .count();
        let existing_esurface_count: usize = all_existing_esurfaces
            .iter()
            .map(ExistingEsurfaces::len)
            .sum();
        let threshold_status = if existing_esurface_count == 0 {
            "below threshold"
        } else {
            "above threshold"
        };
        info!(
            integrand = %self.data.name,
            threshold_status,
            groups_above_threshold,
            group_count = all_existing_esurfaces.len(),
            existing_esurfaces = existing_esurface_count,
            "Input is {threshold_status}: {existing_esurface_count} existing threshold E-surfaces across {groups_above_threshold}/{} graph groups",
            all_existing_esurfaces.len()
        );

        all_existing_esurfaces
    }

    fn validate_runtime_threshold_counterterms(
        &self,
        existing_esurfaces: &TiVec<GroupId, ExistingEsurfaces>,
    ) -> Result<()> {
        for (group_id, group_existing_esurfaces) in existing_esurfaces.iter_enumerated() {
            for group_esurface_id in group_existing_esurfaces.iter().copied() {
                for (graph_group_pos, raised_esurface_id) in self.data.group_derived_data[group_id]
                    .esurface_map[group_esurface_id]
                    .iter_enumerated()
                    .filter_map(|(graph_group_pos, raised_esurface_id)| {
                        raised_esurface_id
                            .map(|raised_esurface_id| (graph_group_pos, raised_esurface_id))
                    })
                {
                    let graph_id = self.data.graph_group_structure[group_id][graph_group_pos];
                    let graph_term = &self.data.graph_terms[graph_id];
                    if !graph_term
                        .threshold_counterterm
                        .local_esurface_exists
                        .get(group_esurface_id)
                        .copied()
                        .unwrap_or(true)
                    {
                        continue;
                    }
                    let is_generated = graph_term
                        .threshold_counterterm
                        .generated_mask
                        .get(raised_esurface_id)
                        .copied()
                        .ok_or_else(|| {
                            eyre!(
                                "Threshold counterterm generation mask is inconsistent for graph '{}' and raised e-surface {}",
                                graph_term.graph.name,
                                raised_esurface_id.0
                            )
                        })?;

                    if !is_generated {
                        return Err(eyre!(
                            "Amplitude integrand '{}' was generated with specialized threshold-subtraction assumptions, but the current runtime model parameters require a trimmed threshold counterterm for graph '{}' and group e-surface {} ({})",
                            self.name(),
                            graph_term.graph.name,
                            group_esurface_id.0,
                            self.data.group_derived_data[group_id].esurface_atoms[group_esurface_id]
                        ))
                        .with_note(|| {
                            "Regenerate the integrand or restore compatible shared/per-integrand model parameters.".to_string()
                        });
                    }
                }
            }
        }

        Ok(())
    }
}

impl ProcessIntegrandImpl for AmplitudeIntegrand {
    type G = AmplitudeGraphTerm;

    fn external_cache_id(&self) -> usize {
        // info!("Getting cache id {}", self.data.external_cache_id);
        self.data.external_cache_id
    }

    fn increment_external_cache_id(&mut self, val: usize) {
        // info!(
        //     "Incrementing cache id {} by {val}",
        //     self.data.external_cache_id
        // );
        self.data.external_cache_id += val
    }

    fn signal_external_momenta_changed(&mut self) {
        self.increment_external_cache_id(1);
        // Update base cache ID when the fundamental configuration changes
        self.data.base_external_cache_id = self.data.external_cache_id;
    }

    fn get_current_external_cache_id(&self) -> usize {
        self.external_cache_id()
    }

    /// Revert to the base external cache ID for the current configuration
    fn revert_to_base_external_cache_id(&mut self) {
        // info!(
        //     "Reverting external cache id from {} to base {}",
        //     self.data.external_cache_id,
        //     self.data.base_external_cache_id
        // );
        self.data.external_cache_id = self.data.base_external_cache_id;
    }

    fn get_base_external_cache_id(&self) -> usize {
        self.data.base_external_cache_id
    }

    fn increment_loop_cache_id(&mut self, val: usize) {
        self.data.loop_cache_id += val
    }

    fn loop_cache_id(&self) -> usize {
        self.data.loop_cache_id
    }

    #[instrument(
          skip_all,
          fields(
              integrand.name = %self.name(),
          )
    )]
    fn warm_up(&mut self, model: &Model) -> Result<()> {
        self.data.rotations = Some(
            Some(Rotation::new(RotationMethod::Identity))
                .into_iter()
                .chain(
                    self.settings
                        .stability
                        .rotation_axis
                        .iter()
                        .map(|axis| Rotation::new(axis.rotation_method())),
                )
                .collect(),
        );

        for a in self.data.graph_terms.iter_mut() {
            a.warm_up(&self.settings, model)?;
        }
        let e_cm = F(self.settings.kinematics.e_cm);
        let constructor = DependentMomentaConstructor::Amplitude(&self.data.external_signature);
        let masses = self.data.graph_terms[0].graph.get_external_masses(model);

        self.settings
            .kinematics
            .externals
            .improve_and_cache(constructor, &masses, &e_cm)?;
        self.warn_on_off_shell_external_states(model)?;

        let thresholds_generated = self
            .data
            .graph_terms
            .iter()
            .all(|term| !term.threshold_counterterm.evaluators.is_empty());

        let is_tree_level = self.data.graph_terms[0].graph.get_loop_number() == 0;
        let existing_esurfaces = if thresholds_generated && !is_tree_level {
            debug!("esurface existence check");
            Some(self.get_existing_esurfaces(model))
        } else {
            None
        };

        if !thresholds_generated && !self.settings.subtraction.disable_threshold_subtraction {
            warn!(
                "Not all graphs have threshold counterterms generated, but threshold subtraction is not disabled. disable runtime threshold subtraction to remove this warning"
            );
            self.settings.subtraction.disable_threshold_subtraction = true;
        }

        if let Some(existing_esurfaces) = &existing_esurfaces {
            let existing_esurface_count: usize =
                existing_esurfaces.iter().map(ExistingEsurfaces::len).sum();
            if self.settings.subtraction.disable_threshold_subtraction
                && existing_esurface_count > 0
            {
                let groups_above_threshold = existing_esurfaces
                    .iter()
                    .filter(|existing_esurfaces| !existing_esurfaces.is_empty())
                    .count();
                warn!(
                    integrand = %self.data.name,
                    groups_above_threshold,
                    group_count = existing_esurfaces.len(),
                    existing_esurfaces = existing_esurface_count,
                    "Input is above threshold, but threshold subtraction is disabled. Turn threshold subtraction on for this input by setting subtraction.disable_threshold_subtraction=false."
                );
            }
        }

        if !self.settings.subtraction.disable_threshold_subtraction && !is_tree_level {
            let existing_esurfaces = existing_esurfaces
                .expect("threshold existence should be checked before runtime threshold setup");
            self.validate_runtime_threshold_counterterms(&existing_esurfaces)?;
            for (group_id, existing_esurfaces) in existing_esurfaces.iter_enumerated() {
                debug!(
                    "solving overlap for group {}, number of thresholds: {}",
                    group_id.0,
                    existing_esurfaces.len()
                );
                if tracing::event_enabled!(tracing::Level::DEBUG) {
                    for group_esurface_id in existing_esurfaces.iter() {
                        let Some((graph_group_pos, raised_esurface_id)) =
                            self.data.group_derived_data[group_id].esurface_map[*group_esurface_id]
                                .iter_enumerated()
                                .filter_map(|(graph_group_pos, raised_esurface_id)| {
                                    raised_esurface_id.map(|raised_esurface_id| {
                                        (graph_group_pos, raised_esurface_id)
                                    })
                                })
                                .find(|(graph_group_pos, _)| {
                                    let graph_id =
                                        self.data.graph_group_structure[group_id][*graph_group_pos];
                                    self.data.graph_terms[graph_id]
                                        .threshold_counterterm
                                        .local_esurface_exists[*group_esurface_id]
                                })
                        else {
                            continue;
                        };
                        let graph_id = self.data.graph_group_structure[group_id][graph_group_pos];
                        let graph_term = &self.data.graph_terms[graph_id];
                        let graph = &graph_term.graph;
                        let esurface_id =
                            graph_term.threshold_counterterm.raised_data.raised_groups
                                [raised_esurface_id]
                                .esurface_ids[0];
                        let lmb_reps = graph.integrand_replacement(
                            &graph.full_filter(),
                            &graph.loop_momentum_basis,
                            &[W_.x___],
                        );
                        let atom =
                            graph_term.esurfaces[esurface_id].lmb_atom_simplified(graph, &lmb_reps);
                        crate::debug_tags!(#integration, #subtraction, #threshold, #inspect, #esurface;
                            stage = "amplitude_threshold_existing_esurface",
                            group_id = group_id.0,
                            group_esurface_id = group_esurface_id.0,
                            graph = %graph.name,
                            graph_group_pos = graph_group_pos.0,
                            raised_esurface_id = raised_esurface_id.0,
                            esurface_id = esurface_id.0,
                            file.atom = %atom,
                            "amplitude threshold existing esurface"
                        );
                    }
                }

                let graph_data = self.data.graph_group_structure[group_id]
                    .into_iter()
                    .map(|graph_id| {
                        let graph = &self.data.graph_terms[graph_id];
                        SingleGraphOverlapData {
                            lmb: &graph.graph.loop_momentum_basis,
                            esurfaces: &graph.esurfaces,
                            raised_data: &graph.threshold_counterterm.raised_data,
                            edge_masses: graph.graph.get_real_mass_vector::<f64>(model),
                        }
                    })
                    .collect();
                let local_esurface_exists = self.data.graph_group_structure[group_id]
                    .into_iter()
                    .map(|graph_id| {
                        self.data.graph_terms[graph_id]
                            .threshold_counterterm
                            .local_esurface_exists
                            .clone()
                    })
                    .collect();

                let overlap_input = OverlapInput {
                    graph_data,
                    settings: &self.settings,
                    group_esurface_map: self.data.group_derived_data[group_id].esurface_map.clone(),
                    local_esurface_exists,
                };

                let external_moms = self
                    .settings
                    .kinematics
                    .externals
                    .get_dependent_externals::<f64>(DependentMomentaConstructor::Amplitude(
                        &self.data.external_signature,
                    ))
                    .expect("could not get externals");

                let overlap =
                    find_maximal_overlap(&overlap_input, existing_esurfaces, &external_moms)
                        .with_context(|| {
                            let readable_esurfaces = existing_esurfaces
                                .iter()
                                .map(|group_esurface_id| {
                                    let (graph_group_pos, raised_esurface_id) = get_representative(
                                        &self.data.group_derived_data[group_id].esurface_map
                                            [*group_esurface_id],
                                    )
                                    .unwrap();
                                    let graph_id =
                                        self.data.graph_group_structure[group_id][graph_group_pos];
                                    let graph_term = &self.data.graph_terms[graph_id];
                                    let graph = &graph_term.graph;
                                    let esurface_id =
                                        graph_term.threshold_counterterm.raised_data.raised_groups
                                            [raised_esurface_id]
                                            .esurface_ids[0];
                                    let lmb_reps = graph.integrand_replacement(
                                        &graph.full_filter(),
                                        &graph.loop_momentum_basis,
                                        &[W_.x___],
                                    );

                                    let esurface = &graph_term.esurfaces[esurface_id];
                                    let atom = esurface.lmb_atom_simplified(graph, &lmb_reps);
                                    (esurface_id, atom)
                                })
                                .collect_vec();

                            let mut msg = format!(
                                "finding overlap for group: {}, existing esurfaces:\n",
                                group_id.0
                            );

                            for readable_esurface in readable_esurfaces {
                                msg += &format!(
                                    "esurface id: {}, atom: {}\n",
                                    readable_esurface.0.0, readable_esurface.1
                                );
                            }

                            msg
                        })?;

                info!(
                    "overlap structure of group {}: [{}]",
                    group_id.0,
                    self.format_overlap_structure(group_id, &overlap)
                );

                let loop_number = self
                    .data
                    .graph_terms
                    .first()
                    .unwrap()
                    .graph
                    .get_loop_number();
                let model_params = self
                    .data
                    .graph_terms
                    .first()
                    .unwrap()
                    .param_builder
                    .pairs
                    .model_parameters
                    .params
                    .clone();

                for graph_id in self.data.graph_group_structure[group_id].into_iter() {
                    let graph_term = &mut self.data.graph_terms[graph_id];
                    let max_required_power = graph_term
                        .threshold_counterterm
                        .raised_data
                        .raised_groups
                        .iter()
                        .map(|raised_group| raised_group.max_occurence)
                        .max()
                        .unwrap_or(0) as i32
                        + 1;

                    let mut localized_overlap = overlap.localized_to_existing_surfaces(
                        &graph_term.threshold_counterterm.local_esurface_exists,
                    );
                    localized_overlap
                        .build_evaluators(
                            &self.data.group_derived_data[group_id].esurface_atoms,
                            &OptimizationSettings::default(),
                            loop_number,
                            external_moms.len(),
                            model_params.clone(),
                            max_required_power,
                        )
                        .with_context(|| {
                            format!(
                                "Failed to build graph-local multi-channeling evaluators for group {} graph {}",
                                group_id.0, graph_term.graph.name
                            )
                        })?;

                    crate::debug_tags!(#integration, #subtraction, #threshold, #inspect, #overlap;
                        stage = "amplitude_threshold_localized_overlap",
                        group_id = group_id.0,
                        graph = %graph_term.graph.name,
                        existing_esurfaces = ?localized_overlap
                            .existing_esurfaces
                            .iter()
                            .map(|group_esurface_id| group_esurface_id.0)
                            .collect_vec(),
                        overlap_group_sizes = ?localized_overlap
                            .overlap_groups
                            .iter()
                            .map(|group| group.existing_esurfaces.len())
                            .collect_vec(),
                        "amplitude threshold localized overlap"
                    );

                    graph_term.threshold_counterterm.overlap = localized_overlap;
                }
            }
        }

        self.event_processing_runtime.set(
            EventProcessingRuntime::from_settings_with_model_and_process_info(
                &self.settings,
                model,
                &histogram_process_info_for_integrand(self),
            )?,
        );

        Ok(())
    }

    fn get_rotations(&self) -> impl Iterator<Item = &Rotation> {
        self.data.rotations.as_ref().expect("forgot warmup").iter()
    }

    fn get_terms_mut(&mut self) -> impl Iterator<Item = &mut Self::G> {
        self.data.graph_terms.iter_mut()
    }

    fn graph_count(&self) -> usize {
        self.data.graph_terms.len()
    }

    fn get_master_graph(&self, group_id: GroupId) -> &Self::G {
        let group_master = self.data.graph_group_structure[group_id].master();
        &self.data.graph_terms[group_master]
    }

    fn get_group_masters(&self) -> impl Iterator<Item = &Self::G> {
        self.data
            .graph_group_structure
            .iter()
            .map(|group| &self.data.graph_terms[group.master()])
    }

    fn get_settings(&self) -> &RuntimeSettings {
        &self.settings
    }

    fn get_graph_mut(&mut self, graph_id: usize) -> &mut Self::G {
        &mut self.data.graph_terms[graph_id]
    }

    fn graph_group_id_for_graph(&self, graph_id: usize) -> Option<usize> {
        self.data.graph_to_group_id.get(graph_id).copied()
    }

    fn get_graph(&self, graph_id: usize) -> &Self::G {
        &self.data.graph_terms[graph_id]
    }

    fn get_group(&self, group_id: GroupId) -> &GraphGroup {
        &self.data.graph_group_structure[group_id]
    }

    fn get_dependent_momenta_constructor(&self) -> DependentMomentaConstructor<'_> {
        DependentMomentaConstructor::Amplitude(&self.data.external_signature)
    }

    fn get_group_structure(&self) -> &TiVec<GroupId, GraphGroup> {
        &self.data.graph_group_structure
    }

    fn take_event_processing_runtime(&mut self) -> Option<EventProcessingRuntime> {
        self.event_processing_runtime.take()
    }

    fn restore_event_processing_runtime(&mut self, runtime: Option<EventProcessingRuntime>) {
        if let Some(runtime) = runtime {
            self.event_processing_runtime.set(runtime);
        }
    }

    fn event_processing_runtime(&self) -> Option<&EventProcessingRuntime> {
        self.event_processing_runtime.as_ref()
    }

    fn event_processing_runtime_mut(&mut self) -> Option<&mut EventProcessingRuntime> {
        self.event_processing_runtime.as_mut()
    }

    fn groups_default_sample_events_by_graph_group(&self) -> bool {
        true
    }
    // fn get_builder_cache(&self) -> &ParamBuilder<f64> {
    //     &self.data.builder_cache
    // }
}

impl HasIntegrand for AmplitudeIntegrand {
    fn create_grid(&self) -> Grid<F<f64>> {
        create_grid(self)
    }

    fn name(&self) -> String {
        self.data.name.clone()
    }

    fn evaluate_sample(
        &mut self,
        sample: &Sample<F<f64>>,
        model: &Model,
        wgt: F<f64>,
        iter: usize,
        use_arb_prec: bool,
        max_eval: Complex<F<f64>>,
    ) -> Result<EvaluationResult> {
        evaluate_sample(self, model, sample, wgt, iter, use_arb_prec, max_eval)
    }

    fn get_n_dim(&self) -> usize {
        if self
            .settings
            .sampling
            .get_parameterization_settings()
            .is_some()
        {
            self.data.graph_terms[0].graph.get_loop_number() * 3
        } else {
            // let dimensions = self
            //     .data
            //     .graph_terms
            //     .iter()
            //     .map(|term| term.get_tropical_sampler().get_dimension())
            //     .sorted()
            //     .collect_vec();

            tracing::warn!(
                "get n dim called for tropical sampling, if groups are enabled this function panics, returning bs value to avoid this"
            );
            69

            //let median_dimension = dimensions[dimensions.len() / 2];
            //median_dimension
        }
    }
}
