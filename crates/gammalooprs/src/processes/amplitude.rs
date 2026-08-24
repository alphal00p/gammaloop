use std::{
    collections::{BTreeMap, HashMap},
    fmt,
    fs::{self, File},
    io::Write,
    iter,
    path::Path,
};

use ahash::AHashSet;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use momtrop::SampleGenerator;

use idenso::dirac::GammaSimplifier;
use rayon::{
    ThreadPool,
    iter::{IndexedParallelIterator, IntoParallelRefMutIterator, ParallelIterator},
};
use spenso::algebra::complex::Complex;
use tracing::{info_span, instrument};
use tracing_indicatif::{span_ext::IndicatifSpanExt, style::ProgressStyle};
use vakint::{EvaluationMethod, NumericalEvaluationResult, Vakint, vakint_symbol};

use crate::{
    GammaLoopContext, GammaLoopContextContainer,
    cff::{
        esurface::{GroupEsurfaceId, RaisedEsurfaceData, RaisedEsurfaceGroup, RaisedEsurfaceId},
        expression::{CFFExpression, OrientationID},
    },
    graph::{
        GraphGroup, GraphGroupPosition, GroupId, LMBext, LmbIndex, LoopMomentumBasis,
        cuts::{CutSet, ResidueSelector},
        threshold_counterterms::{
            ThresholdCountertermMultiplier, ThresholdCountertermSpec, ThresholdCountertermVariant,
        },
    },
    integrands::process::{
        GenericEvaluator, LmbMultiChannelingSetup,
        amplitude::{AmplitudeGraphTerm, AmplitudeIntegrand, AmplitudeIntegrandData},
        graph_to_group_id_for_group_structure,
    },
    model::ArcParticle,
    momentum::{
        sample::{ExternalIndex, SubspaceData},
        signature::SignatureLike,
    },
    processes::{
        DotExportSettings, EvaluatorSettings, GraphGenerationStats, GraphGroupSelectionPlan,
        GraphGroupSelectionSpec, NamedGraphGenerationReport,
        ResolvedThresholdCountertermAssociation, ResolvedThresholdCountertermVariant,
        ResolvedThresholdCounterterms, SingleThresholdPieces, StandaloneExportSettings,
        ThresholdCountertermOrigin, ThresholdCountertermSide, ThresholdCountertermVariantId,
        build_derivative_structure_atom, params_for_derivative_order,
    },
    settings::{
        GlobalSettings, RuntimeSettings, global::OrientationPattern, runtime::LockedRuntimeSettings,
    },
    subtraction::amplitude_counterterm::AmplitudeCountertermAtom,
    utils::{F, GS, Length, W_},
    uv::{
        RenormalizationPart, UVgenerationSettings, UltravioletGraph,
        approx::{CutStructure, OrientationProjection, integrated::to_vakint_integrand},
        settings::VakintSettings,
    },
};
use eyre::{Context, eyre};
use itertools::Itertools;
use linnet::{
    half_edge::{
        involution::{EdgeIndex, Flow, HedgePair},
        subgraph::{ModifySubSet, SuBitGraph, SubGraphLike, SubSetOps},
    },
    parser::DotGraph,
};
use spenso::shadowing::symbolica_utils::LogPrint;
use symbolica::{atom::Var, prelude::*};
use tracing::{debug, info};
use typed_index_collections::{TiVec, ti_vec};

use super::generation_progress::{self, GenerationProcessKind, GenerationProgressPhase};

use crate::{
    cff::esurface::EsurfaceID,
    graph::{FeynmanGraph, Graph},
    integrands::process::ProcessIntegrand,
    model::Model,
    settings::global::GenerationSettings,
};

use crate::graph::parse::complete_group_parsing;

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Amplitude {
    pub name: String,
    pub integrand: Option<ProcessIntegrand>,
    pub graphs: Vec<AmplitudeGraph>,
    pub graph_group_structure: TiVec<GroupId, GraphGroup>,
    pub external_particles: Vec<ArcParticle>,
    pub external_signature: SignatureLike<ExternalIndex>,
    pub group_derived_data: TiVec<GroupId, GroupDerivedData>,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct GroupDerivedData {
    pub esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<RaisedEsurfaceId>>>,
    pub esurface_atoms: TiVec<GroupEsurfaceId, Atom>,
}

impl Amplitude {
    pub fn plan_graph_group_selection(
        &self,
        spec: &GraphGroupSelectionSpec,
    ) -> Result<GraphGroupSelectionPlan> {
        spec.plan(&self.graph_group_structure, |graph_id| {
            self.graphs.get(graph_id).map(|graph| &graph.graph)
        })
    }

    pub fn validate_graph_group_selection_plan(
        &self,
        plan: &GraphGroupSelectionPlan,
    ) -> Result<()> {
        if !self.group_derived_data.is_empty()
            && self.group_derived_data.len() != self.graph_group_structure.len()
        {
            return Err(eyre!(
                "Amplitude '{}' has {} group-derived-data entries for {} graph groups.",
                self.name,
                self.group_derived_data.len(),
                self.graph_group_structure.len()
            ));
        }

        for &old_group_id in plan.retained_group_ids() {
            plan.new_group_id_for_old(old_group_id).ok_or_else(|| {
                eyre!(
                    "Selection plan is missing compact group id for old group {}.",
                    old_group_id.0
                )
            })?;
            if old_group_id.0 >= self.graph_group_structure.len() {
                return Err(eyre!(
                    "Selection plan refers to missing graph group {}.",
                    old_group_id.0
                ));
            }
            let group = &self.graph_group_structure[old_group_id];
            for old_graph_id in group {
                if old_graph_id >= self.graphs.len() {
                    return Err(eyre!(
                        "Graph group {} refers to missing graph id {}.",
                        old_group_id.0,
                        old_graph_id
                    ));
                }
            }
            let master = group.master();
            if master >= self.graphs.len() {
                return Err(eyre!(
                    "Graph group {} refers to missing master graph id {}.",
                    old_group_id.0,
                    master
                ));
            }
        }

        Ok(())
    }

    pub fn apply_graph_group_selection(&mut self, plan: &GraphGroupSelectionPlan) -> Result<()> {
        self.validate_graph_group_selection_plan(plan)?;

        let mut old_graph_to_new_group = vec![None; self.graphs.len()];
        let mut old_graph_is_master = vec![false; self.graphs.len()];
        for &old_group_id in plan.retained_group_ids() {
            let new_group_id = plan.new_group_id_for_old(old_group_id).ok_or_else(|| {
                eyre!(
                    "Selection plan is missing compact group id for old group {}.",
                    old_group_id.0
                )
            })?;
            let group = &self.graph_group_structure[old_group_id];
            for old_graph_id in group {
                old_graph_to_new_group[old_graph_id] = Some(new_group_id);
            }
            old_graph_is_master[group.master()] = true;
        }

        let mut new_graphs = self
            .graphs
            .iter()
            .cloned()
            .enumerate()
            .filter_map(|(old_graph_id, mut graph)| {
                let new_group_id = old_graph_to_new_group[old_graph_id]?;
                graph.graph.group_id = Some(new_group_id);
                graph.graph.is_group_master = old_graph_is_master[old_graph_id];
                if let Some(multi_channeling_setup) = &mut graph.derived_data.multi_channeling_setup
                {
                    multi_channeling_setup.graph.group_id = Some(new_group_id);
                    multi_channeling_setup.graph.is_group_master =
                        old_graph_is_master[old_graph_id];
                }
                Some(graph)
            })
            .collect::<Vec<_>>();

        let mut parsed_graphs = new_graphs
            .iter()
            .map(|graph| graph.graph.clone())
            .collect::<Vec<_>>();
        let new_graph_group_structure = complete_group_parsing(&mut parsed_graphs)?;
        for (graph, parsed_graph) in new_graphs.iter_mut().zip(parsed_graphs) {
            graph.graph.group_id = parsed_graph.group_id;
            graph.graph.is_group_master = parsed_graph.is_group_master;
            if let Some(multi_channeling_setup) = &mut graph.derived_data.multi_channeling_setup {
                multi_channeling_setup.graph.group_id = graph.graph.group_id;
                multi_channeling_setup.graph.is_group_master = graph.graph.is_group_master;
            }
        }

        let new_group_derived_data = if self.group_derived_data.is_empty() {
            TiVec::new()
        } else {
            plan.retained_group_ids()
                .iter()
                .map(|&old_group_id| self.group_derived_data[old_group_id].clone())
                .collect::<TiVec<GroupId, _>>()
        };

        self.graphs = new_graphs;
        self.graph_group_structure = new_graph_group_structure;
        self.group_derived_data = new_group_derived_data;
        self.integrand = None;
        Ok(())
    }

    pub fn export_standalone(
        &self,
        path: impl AsRef<Path>,
        settings: &StandaloneExportSettings,
    ) -> Result<()> {
        if let Some(integrand) = &self.integrand {
            integrand.export_standalone(path, settings)?
        } else {
            return Err(eyre!(
                "Cannot warm up amplitude {} without integrand",
                self.name
            ));
        }

        Ok(())
    }

    #[instrument(
          skip_all,
          fields(
              amplitude.name = %self.name,
          )
    )]
    pub(crate) fn warm_up(&mut self, model: &Model) -> Result<()> {
        if let Some(integrand) = &mut self.integrand {
            integrand.warm_up(model)
        } else {
            Err(eyre!(
                "Cannot warm up amplitude {} without integrand",
                self.name
            ))
        }
    }

    #[instrument(
          skip_all,
          fields(
              path = %path.as_ref().display(),
          )
    )]
    pub(crate) fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("amp.bin"))?;
        let (mut amp, _): (Self, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;

        if path.as_ref().join("integrand").exists() {
            let integrand = AmplitudeIntegrand::load(path.as_ref().join("integrand"), context)?;
            amp.integrand = Some(ProcessIntegrand::Amplitude(integrand));
        }

        Ok(amp)
    }

    #[instrument(
          skip_all,
          fields(
              amplitude.name = %self.name,
          )
    )]
    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        info!("Compiling amplitude {}", self.name);
        let p = path.as_ref().join(&self.name);

        let r = fs::create_dir_all(&p).with_context(|| {
            format!(
                "Trying to create directory to save amplitude {}",
                p.display()
            )
        });
        if override_existing {
            r?;
        }
        if let Some(integrand) = &mut self.integrand {
            let compile_times = integrand.compile(&p, override_existing, thread_pool)?;
            return Ok(compile_times
                .into_iter()
                .map(|(graph_name, duration)| NamedGraphGenerationReport {
                    integrand_name: self.name.clone(),
                    graph_name,
                    stats: GraphGenerationStats {
                        total_time: duration,
                        evaluator_compile_time: duration,
                        ..GraphGenerationStats::default()
                    },
                })
                .collect());
        };
        Ok(Vec::new())
    }

    #[instrument(
          skip_all,
          fields(
              amplitude.name = %self.name,
              path = %path.as_ref().display(),
          )
    )]
    pub fn save(&mut self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let p = path.as_ref().join(&self.name);

        let r = fs::create_dir_all(&p).with_context(|| {
            format!(
                "Trying to create directory to save amplitude {}",
                p.display()
            )
        });
        if override_existing {
            r?;
        }

        let integrand = self.integrand.take();
        if let Some(integrand) = &integrand {
            integrand.save(&p, override_existing)?;
        };

        let binary = bincode::encode_to_vec(&(*self), bincode::config::standard())?;
        if override_existing {
            fs::write(p.join("amp.bin"), binary)?;
        } else {
            let mut file = File::create_new(p.join("amp.bin"))?;
            file.write_all(&binary)?;
        }

        self.integrand = integrand;
        Ok(())
    }

    #[instrument(
        skip_all,
        fields(
             amplitude.name = %self.name,
        )
    )]
    pub fn preprocess(
        &mut self,
        model: &Model,
        settings: &GenerationSettings,
        locked_runtime_settings: &LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        // preprocess each graph individually
        let integrand_name = self.name.clone();

        let preprocess_span = if generation_progress::detailed_progress_enabled() {
            let span = info_span!("Preprocessing graphs", indicatif.pb_show = true);
            span.pb_set_style(&ProgressStyle::with_template(
                "{wide_bar} {pos}/{len} {msg}",
            )?);
            span.pb_set_length(self.graphs.len() as u64);
            span.pb_set_message("Preprocessing graphs");
            Some(span)
        } else {
            None
        };
        let preprocess_span_enter = preprocess_span.as_ref().map(|span| span.enter());

        let preprocess_reports = thread_pool.install(|| {
            let parent = preprocess_span.clone();
            self.graphs
                .par_iter_mut()
                .map(|amplitude_graph| {
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    let graph_name = amplitude_graph.graph.name.clone();
                    generation_progress::graph_started(
                        GenerationProcessKind::Amplitude,
                        &integrand_name,
                        &graph_name,
                        None,
                    );
                    let _guard = parent.as_ref().map(|span| span.enter());
                    let stats =
                        amplitude_graph.preprocess(model, settings, locked_runtime_settings);
                    if let Some(span) = &parent {
                        span.pb_inc(1);
                    }

                    let stats = stats?;
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    generation_progress::graph_finished(
                        GenerationProcessKind::Amplitude,
                        &integrand_name,
                        &graph_name,
                        &stats,
                        None,
                    );

                    Ok(NamedGraphGenerationReport {
                        integrand_name: integrand_name.clone(),
                        graph_name,
                        stats,
                    })
                })
                .collect::<Result<Vec<_>>>()
        })?;

        drop(preprocess_span_enter);
        drop(preprocess_span);

        self.generate_grouped_derived_data()?;

        Ok(preprocess_reports)
    }

    #[instrument(
        skip_all,
          fields(
              amplitude.name = %self.name,
          )
      )]
    pub fn build_integrand(
        &mut self,
        model: &Model,
        process_name: &str,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        let started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #graph, #summary;
            stage = "amplitude_build_integrand_start",
            integrand = %self.name,
            graph_count = self.graphs.len(),
            "Generation timing milestone"
        );
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        let integrand_name = self.name.clone();
        generation_progress::begin_phase(
            GenerationProgressPhase::GraphGeneration,
            GenerationProcessKind::Amplitude,
            process_name,
            &integrand_name,
            self.graphs.len(),
            None,
        );
        let mut graph_reports = Vec::new();
        let terms: Vec<_> = thread_pool.install(|| {
            self.graphs
                .par_iter_mut()
                .enumerate()
                .map(|(graph_id, graph)| {
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    let graph_started = std::time::Instant::now();
                    let group_id = graph.graph.group_id.unwrap(); // should always be set
                    let esurface_map = &self.group_derived_data[group_id].esurface_map;
                    let group_pos = self.graph_group_structure[group_id]
                        .find_position(graph_id)
                        .unwrap();

                    crate::debug_tags!(#generation, #profile, #graph, #summary;
                        stage = "amplitude_generate_term_for_graph_start",
                        integrand = %integrand_name,
                        graph = %graph.graph.name,
                        graph_id,
                        group_id = %group_id.0,
                        "Generation timing milestone"
                    );
                    generation_progress::graph_started(
                        GenerationProcessKind::Amplitude,
                        &integrand_name,
                        &graph.graph.name,
                        None,
                    );
                    let _progress_context_guard =
                        generation_progress::enter_progress_context(graph.graph.name.clone());
                    let (term, mut stats) = graph.generate_term_for_graph(
                        model,
                        group_pos,
                        esurface_map.clone(),
                        global_settings,
                    )?;
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    stats.evaluator_count = term.generic_evaluator_count();
                    stats.total_time += graph_started.elapsed();
                    crate::debug_tags!(#generation, #profile, #graph, #summary;
                        stage = "amplitude_generate_term_for_graph_done",
                        integrand = %integrand_name,
                        graph = %graph.graph.name,
                        graph_id,
                        group_id = %group_id.0,
                        elapsed_ms = graph_started.elapsed().as_secs_f64() * 1000.0,
                        "Generation timing milestone"
                    );
                    generation_progress::graph_finished(
                        GenerationProcessKind::Amplitude,
                        &integrand_name,
                        &graph.graph.name,
                        &stats,
                        None,
                    );
                    Ok((
                        term,
                        NamedGraphGenerationReport {
                            integrand_name: integrand_name.clone(),
                            graph_name: graph.graph.name.clone(),
                            stats,
                        },
                    ))
                })
                .collect::<Result<Vec<_>>>()
        })?;
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        for (_, report) in &terms {
            graph_reports.push(report.clone());
        }
        let mut terms = terms.into_iter().map(|(term, _)| term).collect::<Vec<_>>();

        for group in self.graph_group_structure.iter() {
            let master_graph_id = group.master();
            let mc_of_master = self.graphs[master_graph_id]
                .derived_data
                .multi_channeling_setup
                .as_ref()
                .unwrap();
            let master_graph = &self.graphs[master_graph_id].graph;
            let master_external_signature = master_graph.get_external_signature();
            let master_external_pdgs = master_graph
                .get_external_partcles()
                .into_iter()
                .map(|particle| particle.pdg_code)
                .collect_vec();

            for graph_id in group.into_iter() {
                terms[graph_id].multi_channeling_setup = mc_of_master.clone();
                terms[graph_id].master_external_signature = master_external_signature.clone();
                terms[graph_id].master_external_pdgs = master_external_pdgs.clone();
            }
        }

        let backend_started = std::time::Instant::now();
        let graph_count = terms.len();
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_prepare_runtime_backends_start",
            integrand = %self.name,
            graph_count,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        generation_progress::backend_started(
            GenerationProcessKind::Amplitude,
            &self.name,
            graph_count,
        );
        let mut amplitude_integrand = AmplitudeIntegrand {
            settings: runtime_default.into_with_modified_kinematics(
                &self.external_signature,
                &self.graphs[0].graph.get_external_masses(model),
            )?,
            data: AmplitudeIntegrandData {
                name: self.name.clone(),
                compilation: global_settings
                    .generation
                    .compile
                    .frozen_mode(&global_settings.generation.evaluator),
                rotations: None,
                loop_cache_id: 0,
                external_cache_id: 0,
                base_external_cache_id: 0,
                graph_terms: terms,
                external_signature: self.external_signature.clone(),
                graph_group_structure: self.graph_group_structure.clone(),
                graph_to_group_id: graph_to_group_id_for_group_structure(
                    &self.graph_group_structure,
                ),
                group_derived_data: self.group_derived_data.clone(),
            },
            event_processing_runtime: Default::default(),
            active_f64_backend: Default::default(),
        };
        let compile_times =
            amplitude_integrand.prepare_runtime_backends_after_generation_with_compile_times()?;
        crate::debug_tags!(#generation, #profile, #compile, #graph, #summary;
            stage = "amplitude_prepare_runtime_backends_done",
            integrand = %self.name,
            graph_count,
            elapsed_ms = backend_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        generation_progress::backend_finished(
            GenerationProcessKind::Amplitude,
            &self.name,
            backend_started.elapsed(),
        );
        for (report, compile_time) in graph_reports.iter_mut().zip(compile_times) {
            report.stats.evaluator_compile_time += compile_time;
            report.stats.total_time += compile_time;
        }
        self.integrand = Some(ProcessIntegrand::Amplitude(amplitude_integrand));
        crate::debug_tags!(#generation, #profile, #graph, #summary;
            stage = "amplitude_build_integrand_done",
            integrand = %self.name,
            graph_count = self.graphs.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        Ok(graph_reports)
    }

    #[instrument(
        skip_all,
          fields(
              amplitude.name = %self.name,
          )
      )]
    #[allow(dead_code)]
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::io::Error> {
        for graph in &self.graphs {
            graph.write_dot(writer, settings)?;
            writeln!(writer)?;
        }
        Ok(())
    }

    #[instrument(
        skip_all,
          fields(
              amplitude.name = %self.name,
          )
      )]
    pub fn write_dot_fmt<W: fmt::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::fmt::Error> {
        for graph in &self.graphs {
            graph.write_dot_fmt(writer, settings)?;
            writeln!(writer)?;
        }
        Ok(())
    }

    pub fn generate_grouped_derived_data(&mut self) -> Result<()> {
        // for each group we must collect all inequivalent esurfaces.

        let group_derived_data = self
            .graph_group_structure
            .iter()
            .map(|group| {
                let mut group_esurface_structure =
                    BTreeMap::<Atom, TiVec<GraphGroupPosition, Option<RaisedEsurfaceId>>>::default(
                    );

                for (graph_group_position, graph_id) in group.iter_enumerated() {
                    let amplitude_graph = &self.graphs[graph_id];
                    let lmb_reps = amplitude_graph.graph.integrand_replacement(
                        &amplitude_graph.graph.full_filter(),
                        &amplitude_graph.graph.loop_momentum_basis,
                        &[W_.x___],
                    );

                    let esurfaces = &amplitude_graph.graph.surface_cache.esurface_cache;

                    for (raised_esurface_id, raised_group) in amplitude_graph
                        .derived_data
                        .raised_data
                        .raised_groups
                        .iter_enumerated()
                    {
                        let esurface = &esurfaces[raised_group.esurface_ids[0]];
                        let esurface_atom = esurface.lmb_atom(&amplitude_graph.graph, &lmb_reps);

                        group_esurface_structure
                            .entry(esurface_atom)
                            .or_insert(ti_vec![None; group.len()])[graph_group_position] =
                            Some(raised_esurface_id);
                    }
                }

                let (surface_atoms, esurface_map) = group_esurface_structure.into_iter().unzip();

                GroupDerivedData {
                    esurface_map,
                    esurface_atoms: surface_atoms,
                }
            })
            .collect::<TiVec<GroupId, _>>();

        let _: () = self.group_derived_data = group_derived_data;
        Ok(())
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait= GammaLoopContext)]
pub struct AmplitudeGraph {
    pub graph: Graph,
    pub derived_data: AmplitudeDerivedData,
}

#[derive(Clone)]
struct ResolvedAmplitudeThresholdCountertermDraft {
    name: String,
    origin: ThresholdCountertermOrigin,
    disable: bool,
    requested_subspace: Option<Vec<EdgeIndex>>,
    requested_parent_lmb: Option<Vec<EdgeIndex>>,
    multiplier: Option<ThresholdCountertermMultiplier>,
    subspace: SubspaceData,
}

impl ResolvedAmplitudeThresholdCountertermDraft {
    fn is_compatible_with(
        &self,
        other: &Self,
        all_lmbs: &TiVec<LmbIndex, LoopMomentumBasis>,
    ) -> bool {
        self.name == other.name
            && self.disable == other.disable
            && self.multiplier == other.multiplier
            && self
                .subspace
                .has_equivalent_embedding(&other.subspace, all_lmbs)
    }
}

struct ResolvedAmplitudeThresholdAssociationDraft {
    esurface_id: EsurfaceID,
    threshold_edges: Vec<EdgeIndex>,
    origin: ThresholdCountertermOrigin,
    variant_subspaces: Vec<SubspaceData>,
}

struct ResolvedAmplitudeThresholdGroupDraft {
    raised_esurface_group: RaisedEsurfaceGroup,
    associations: Vec<ResolvedAmplitudeThresholdAssociationDraft>,
    variants: Vec<ResolvedAmplitudeThresholdCountertermDraft>,
}

struct AmplitudeThresholdCountertermBuild {
    legacy_counterterms: TiVec<RaisedEsurfaceId, AmplitudeCountertermAtom>,
    variants: TiVec<ThresholdCountertermVariantId, AmplitudeThresholdCountertermVariant>,
    resolved: ResolvedThresholdCounterterms,
    raised_esurface_ids: TiVec<EsurfaceID, RaisedEsurfaceId>,
}

pub struct AnalyticalEvaluationConfig<'a> {
    pub model: &'a Model,
    pub refresh_model_values: bool,
    pub evaluate_numerically: bool,
    pub vakint: &'a Vakint,
    pub true_settings: &'a vakint::VakintSettings,
    pub settings: &'a VakintSettings,
    pub run_time_settings: &'a RuntimeSettings,
    pub include_global_numerator: bool,
}

impl AmplitudeGraph {
    pub(crate) fn new(graph: Graph) -> Self {
        AmplitudeGraph {
            graph,
            derived_data: AmplitudeDerivedData {
                all_mighty_integrand: Atom::Zero,
                cff_expression: None,

                lmbs: None,
                tropical_sampler: None,
                multi_channeling_setup: None,
                threshold_counterterms: TiVec::new(),
                threshold_counterterm_variants: TiVec::new(),
                resolved_threshold_counterterms: None,
                raised_data: RaisedEsurfaceData {
                    raised_groups: TiVec::new(),
                    pass_two_evaluator: None,
                },
                raised_esurface_ids: TiVec::new(),
            },
        }
    }
}

impl AmplitudeGraph {
    pub fn renormalization_part(
        &mut self,
        settings: &UVgenerationSettings,
    ) -> Result<RenormalizationPart> {
        if self.derived_data.cff_expression.is_none() {
            self.generate_cff(&OrientationPattern::default())?;
        }
        let valid_orientations: Vec<_> = self
            .derived_data
            .cff_expression
            .as_ref()
            .expect("cff_expression should have been created")
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect();

        settings.orchestrator.renormalization_part(
            &mut self.graph,
            OrientationProjection::new(&valid_orientations, &OrientationPattern::default()),
            settings,
        )
    }

    #[allow(dead_code)]
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::io::Error> {
        if let Some(graph) = self.graph_with_materialized_threshold_counterterms(settings) {
            graph.dot_serialize_io(writer, settings)
        } else {
            self.graph.dot_serialize_io(writer, settings)
        }
    }

    pub(crate) fn write_dot_fmt<W: fmt::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::fmt::Error> {
        if let Some(graph) = self.graph_with_materialized_threshold_counterterms(settings) {
            graph.dot_serialize_fmt(writer, settings)
        } else {
            self.graph.dot_serialize_fmt(writer, settings)
        }
    }

    fn graph_with_materialized_threshold_counterterms(
        &self,
        settings: &DotExportSettings,
    ) -> Option<Graph> {
        if !settings.include_autogenerated_fields {
            return None;
        }
        let resolved = self.derived_data.resolved_threshold_counterterms.as_ref()?;
        let all_lmbs = self.derived_data.lmbs.as_ref()?;
        let mut graph = self.graph.clone();
        graph.threshold_counterterms = crate::graph::autogen::Autogen::explicit(
            resolved.materialized_spec(&self.graph.threshold_counterterms, all_lmbs),
        );
        Some(graph)
    }

    #[instrument(skip_all, err)]
    pub(crate) fn generate_cff(&mut self, orientation_pattern: &OrientationPattern) -> Result<()> {
        let _progress_guard = generation_progress::enter_detailed_progress_span("Generating CFF");
        let shift_rewrite = self
            .graph
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let contract_edges = self
            .graph
            .iter_edges_of(
                &self
                    .graph
                    .tree_edges
                    .subtract(&self.graph.initial_state_cut)
                    .subtract(&self.graph.external_filter::<SuBitGraph>()),
            )
            .map(|x| x.1)
            .collect_vec();

        let cff_expression =
            self.graph
                .generate_cff(&contract_edges, &shift_rewrite, orientation_pattern)?;

        self.derived_data.cff_expression = Some(cff_expression);

        Ok(())
    }

    #[instrument(skip_all, err)]
    pub(crate) fn preprocess(
        &mut self,
        model: &Model,
        settings: &GenerationSettings,
        locked_runtime_settings: &LockedRuntimeSettings,
    ) -> Result<GraphGenerationStats> {
        let _progress_guard = generation_progress::enter_detailed_progress_span("preprocessing");
        let preprocess_started = std::time::Instant::now();
        let vk = crate::utils::vakint()?;

        self.derived_data.threshold_counterterms = TiVec::new();
        self.derived_data.threshold_counterterm_variants = TiVec::new();
        self.derived_data.resolved_threshold_counterterms = None;
        self.derived_data.raised_esurface_ids = TiVec::new();
        self.derived_data.raised_data = RaisedEsurfaceData {
            raised_groups: TiVec::new(),
            pass_two_evaluator: None,
        };

        self.generate_cff(&settings.orientation_pattern)?;

        // UV orchestration can extend the graph surface cache, while raised IDs
        // belong to the CFF expression generated above.
        let raised_data = settings.threshold_subtraction.enable_thresholds.then(|| {
            self.graph.determine_raised_esurfaces_from_expression(
                self.derived_data
                    .cff_expression
                    .as_ref()
                    .expect("cff_expression should have been created"),
            )
        });

        self.build_integrands(settings, vk)?;

        if self.graph.is_group_master {
            self.build_tropical_sampler(settings)?;
        }

        self.build_lmbs();

        if self.graph.is_group_master {
            self.build_multi_channeling_channels(settings.override_lmb_heuristics);
        }

        if let Some(mut raised_data) = raised_data {
            let max_order = raised_data
                .raised_groups
                .iter()
                .map(|raised_group| raised_group.max_occurence)
                .max()
                .unwrap_or(0);
            if max_order > 1 {
                self.graph.param_builder.initialize_duals(max_order);
            }
            raised_data.pass_two_evaluator = Some(
                (1..=max_order)
                    .map(|order| {
                        threshold_counterterm_helper(
                            order as u8,
                            self.graph.get_loop_number(),
                            &settings.evaluator,
                        )
                    })
                    .collect(),
            );

            let build = self.build_threshold_counterterm_parametric_integrand(
                &raised_data,
                settings,
                vk,
                locked_runtime_settings,
                model,
            )?;
            self.derived_data.threshold_counterterms = build.legacy_counterterms;
            self.derived_data.threshold_counterterm_variants = build.variants;
            self.derived_data.resolved_threshold_counterterms = Some(build.resolved);
            self.derived_data.raised_esurface_ids = build.raised_esurface_ids;
            self.derived_data.raised_data = raised_data;
        }

        Ok(GraphGenerationStats {
            total_time: preprocess_started.elapsed(),
            ..GraphGenerationStats::default()
        })
    }

    #[instrument(skip_all)]
    fn build_multi_channeling_channels(&mut self, override_lmb_heuristics: bool) {
        let _progress_guard =
            generation_progress::enter_detailed_progress_span("Building Multi-Channeling Channels");
        let channels = self.graph.build_multi_channeling_channels(
            self.derived_data.lmbs.as_ref().unwrap(),
            override_lmb_heuristics,
        );

        self.derived_data.multi_channeling_setup = Some(channels)
    }

    //  pub fn ct_params(&self, model: &Model) -> ParamBuilder<f64> {
    //      let mut param_builder = self.param_builder_core(model);
    //      param_builder.uv_damp_atom(vec![
    //          Atom::var(GS.uv_damp_plus),
    //          Atom::var(GS.uv_damp_minus),
    //      ]);
    //      param_builder.derivative_at_tstar_atom(Atom::var(GS.deta));
    //      param_builder.radius_atom(Atom::var(GS.radius));
    //      param_builder.radius_star_atom(Atom::var(GS.radius_star));
    //      param_builder.h_function_atom(Atom::var(GS.hfunction));
    //      param_builder
    //  }

    //  pub(crate) fn param_builder_core(&self, model: &Model) -> ParamBuilder<f64> {
    //      // the float type does not matter here
    //      let mut param_builder = ParamBuilder::<f64>::new_empty();

    //      // this is wrong if we allow for vacuum graphs
    //      param_builder.external_energies_atom(&self.graph);
    //      param_builder.orientation_params(&self.graph);

    //      // param_builder.polarizations(&self.graph);
    //      // param_builder.polar

    //      // spatial components of external momenta
    //      param_builder.external_spatial_atom(&self.graph);
    //      param_builder.polarization_params(&self.graph);
    //      // spatial EMR
    //      param_builder.emr_spatial_atom(&self.graph);

    //      param_builder.model_parameters_atom(model);

    //      param_builder.m_uv_atom(Atom::var(GS.m_uv_vacuum));

    //      param_builder.mu_r_sq_atom(Atom::var(GS.mu_r_sq));

    //      self.add_function_map(&mut param_builder);

    //      param_builder
    //  }

    // pub fn fill_in_params(
    //     &self,
    //     param_builder: &mut ParamBuilder,
    //     model: &Model,
    //     settings: &Settings,
    // ) {
    //     param_builder.external_energies_value(momentum_sample);
    // }

    // fn get_eager_const_map(&self)->HashM

    pub fn to_numerical(
        numerical_result: AtomView,
        true_settings: &vakint::VakintSettings,
    ) -> Result<NumericalEvaluationResult> {
        Ok(NumericalEvaluationResult::from_atom(
            numerical_result,
            vakint_symbol!(&true_settings.epsilon_symbol),
            true_settings,
        )?)
    }

    pub fn analytical_evaluation<S: SubGraphLike<Base = SuBitGraph> + SubSetOps>(
        &self,
        component: &S,
        config: AnalyticalEvaluationConfig<'_>,
    ) -> Result<Atom> {
        let mut true_settings = config.true_settings.clone();
        true_settings.number_of_terms_in_epsilon_expansion =
            self.graph.n_loops(&self.graph.no_dummy()) as i64 + 1;
        let pysec_dec_enabled_in_vakint = true_settings.evaluation_order.0.iter().find_map(|o| {
            if let EvaluationMethod::PySecDec(opts) = o {
                Some(opts)
            } else {
                None
            }
        });

        let complex_params_vakint =
            if config.evaluate_numerically || pysec_dec_enabled_in_vakint.is_some() {
                let mut param_builder = self.graph.param_builder.clone(); //ParamBuilder::<f64>::new(&self.graph, model);
                if config.refresh_model_values {
                    param_builder.update_model_values(config.model);
                }
                param_builder.m_uv_value(Complex::new_re(F(config.run_time_settings.general.m_uv)));
                param_builder.renormalization_localization_scale_value(Complex::new_re(F(config
                    .run_time_settings
                    .general
                    .renormalization_localization_scale)));
                param_builder.mu_r_sq_value(Complex::new_re(F(config
                    .run_time_settings
                    .general
                    .mu_r_sq())));

                // println!("\nParamBuilder parameters:\n{}", param_builder);

                let mut complex_params: HashMap<String, symbolica::domains::float::Complex<f64>> =
                    HashMap::default();
                for params in param_builder.pairs.into_iter() {
                    let ps: crate::integrands::process::ParamValuePairs = params;
                    for (p_name, p_value) in
                        ps.params
                            .iter()
                            .zip(ps.value_range)
                            .map(|(a, value_index)| {
                                (
                                    a.to_canonical_string(),
                                    param_builder.values[0][value_index],
                                )
                            })
                    {
                        complex_params.insert(
                            p_name,
                            symbolica::domains::float::Complex::new(
                                p_value.re.into(),
                                p_value.im.into(),
                            ),
                        );
                    }
                }
                // Make sure to remove entries already supported by vakint, as they may not match required precision
                for atom in &[
                    Atom::Var(Var::new(Symbol::PI)),
                    Atom::Var(Var::new(vakint_symbol!("EulerGamma"))),
                    function!(Symbol::LOG, Atom::num(2)),
                ] {
                    _ = complex_params.remove(&atom.to_string());
                }

                // Make sure to properly do the upcasting to required precision in vakint settings
                config
                    .vakint
                    .params_from_complex_f64(&true_settings, &complex_params)
            } else {
                HashMap::default()
            };

        if let Some(pysec_dec_opts) = pysec_dec_enabled_in_vakint {
            true_settings.evaluation_order.adjust(
                None,
                pysec_dec_opts.relative_precision,
                &HashMap::default(),
                &complex_params_vakint,
                &HashMap::default(),
            );
        }

        let mut num = self
            .graph
            .numerator(component, &self.graph.empty_subgraph());
        if config.include_global_numerator {
            num.state.expr *= &self.graph.global_prefactor.num;
        }

        let before_gamma = num.to_d_dim(GS.dim).get_single_atom().unwrap();
        let before_gamma_plain = before_gamma.to_plain_string();
        let four_dimensional_numerator = before_gamma.simplify_gamma();
        let after_gamma_plain = four_dimensional_numerator.to_plain_string();
        crate::debug_tags!(#uv, #integrated, #vakint, #profile, #trace;
            stage = "amplitude_to_vakint_after_simplify_gamma",
            changed = before_gamma_plain != after_gamma_plain,
            before_gamma_count = %before_gamma_plain.matches("spenso::gamma").count(),
            after_gamma_count = %after_gamma_plain.matches("spenso::gamma").count(),
            before_chain_count = %before_gamma_plain.matches("spenso::chain").count(),
            after_chain_count = %after_gamma_plain.matches("spenso::chain").count(),
            log.before_gamma = before_gamma,
            log.after_gamma = four_dimensional_numerator,
            "Gamma simplification before Vakint"
        );

        let mut four_dimensional_integrand = four_dimensional_numerator
            / self
                .graph
                .denominator(component, |e| e.extra_data.vakint_edge_power.unwrap_or(1));

        // println!("Four-dimensional integrand: {}", four_dimensional_integrand);

        let component_lmb = self.graph.lmb_of(component);
        let mom_reps = self.graph.uv_wrapped_replacement(
            &self.graph.full_filter(),
            &component_lmb,
            &[W_.x___],
        );

        // println!("Reps:");
        // for r in &mom_reps {
        //     println!("{r}");
        // }

        // rewrite the inner_t as well
        four_dimensional_integrand = four_dimensional_integrand.replace_multiple(&mom_reps);

        // println!("LMB: {}", lmb);
        // let vk_mom = vakint_symbol!("k");
        // for (i, l) in lmb.loop_edges.iter().enumerate() {
        //     four_dimensional_integrand = four_dimensional_integrand
        //         .replace(function!(GS.emr_mom, usize::from(*l) as i64))
        //         .with(function!(vk_mom, i as i64 + 1))
        //         .replace(function!(
        //             GS.emr_mom,
        //             function!(vk_mom, i as i64 + 1),
        //             W_.x___
        //         ))
        //         .with(function!(vk_mom, i as i64 + 1, W_.x___));
        // }

        let mut vakint_integrand = to_vakint_integrand(
            &four_dimensional_integrand,
            &self.graph,
            &self.graph.full_filter(),
            &self.graph.empty_subgraph::<SuBitGraph>(),
            config.settings,
            false,
        )?;

        vakint_integrand.canonicalize(&true_settings, &config.vakint.topologies, false)?;
        // println!("Canonized: {}", vakint_integrand);
        vakint_integrand.tensor_reduce(config.vakint, &true_settings)?;
        // println!("Tensor Reduced {}", vakint_integrand);
        vakint_integrand.evaluate_integral(config.vakint, &true_settings)?;
        // println!("Evaluated {}", vakint_integrand);
        let analytical_evaluation: Atom = vakint_integrand.into();
        // println!(
        //     "\nVakint analytical evaluation:\n{:#}",
        //     analytical_evaluation
        // );
        if !config.evaluate_numerically {
            Ok(analytical_evaluation)
        } else {
            let (numerical_evaluation, _error) = config
                .vakint
                .numerical_evaluation(
                    &true_settings,
                    analytical_evaluation.as_view(),
                    &HashMap::default(),
                    &complex_params_vakint,
                    None,
                )
                .unwrap();

            // println!("\nVakint numerical evaluation:\n{:#}", numerical_evaluation);

            let numerical_evaluation_atom =
                numerical_evaluation.to_atom(vakint_symbol!(true_settings.epsilon_symbol.clone()));

            Ok(numerical_evaluation_atom)
        }
    }

    #[instrument(skip_all, err)]
    pub(crate) fn build_integrands(
        &mut self,
        settings: &GenerationSettings,
        vakint: &Vakint,
    ) -> Result<()> {
        let _progress_guard =
            generation_progress::enter_detailed_progress_span("Building Parametric Integrand");
        let started = std::time::Instant::now();
        crate::debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "amplitude_graph_build_integrands_start",
            graph = %self.graph.name,
            subtract_uv = settings.uv.subtract_uv,
            generate_integrated = settings.uv.generate_integrated,
            only = %settings.uv.final_integrand,
            "Generation timing milestone"
        );
        let valid_orientations: Vec<_> = self
            .derived_data
            .cff_expression
            .as_ref()
            .expect("cff_expression should have been created")
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect();
        crate::debug_tags!(#generation, #profile, #graph, #orientation, #summary;
            stage = "amplitude_graph_valid_orientations_done",
            graph = %self.graph.name,
            orientation_count = valid_orientations.len(),
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );
        let cutstructure = CutStructure::empty(&self.graph);
        let orchestration_started = std::time::Instant::now();
        let parametric_exprs = settings.uv.orchestrator.parametric_integrands(
            &mut self.graph,
            cutstructure,
            vakint,
            OrientationProjection::new(&valid_orientations, &settings.orientation_pattern),
            &settings.uv,
        )?;
        crate::debug_tags!(#generation, #profile, #uv, #graph, #summary;
            stage = "amplitude_graph_parametric_orchestration_done",
            graph = %self.graph.name,
            parametric_integrand_count = parametric_exprs.len(),
            elapsed_ms = orchestration_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );

        let normalization_started = std::time::Instant::now();
        let exprs: Vec<_> = parametric_exprs.into_iter().collect();
        crate::debug_tags!(#generation, #profile, #graph, #summary;
            stage = "amplitude_graph_cff_normalization_done",
            graph = %self.graph.name,
            expr_count = exprs.len(),
            elapsed_ms = normalization_started.elapsed().as_secs_f64() * 1000.0,
            total_elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );

        let assign_started = std::time::Instant::now();
        let integrands = exprs.into_iter().next().unwrap().integrands;
        self.derived_data.all_mighty_integrand = integrands.iter().next().unwrap().1.clone(); // should be exactly one expression
        crate::debug_tags!(#generation, #profile, #graph, #summary;
            stage = "amplitude_graph_build_integrands_done",
            graph = %self.graph.name,
            assign_elapsed_ms = assign_started.elapsed().as_secs_f64() * 1000.0,
            elapsed_ms = started.elapsed().as_secs_f64() * 1000.0,
            "Generation timing milestone"
        );

        Ok(())
    }

    fn configured_amplitude_threshold_variants(
        spec: &ThresholdCountertermSpec,
        threshold_edges: &[EdgeIndex],
    ) -> Vec<(ThresholdCountertermVariant, ThresholdCountertermOrigin)> {
        let configured = spec
            .cuts
            .iter()
            .find(|cut| cut.edges.is_empty())
            .and_then(|cut| {
                cut.thresholds
                    .iter()
                    .find(|threshold| threshold.edges == threshold_edges)
            })
            .filter(|threshold| !threshold.counterterms.is_empty());

        match configured {
            Some(threshold) => threshold
                .counterterms
                .iter()
                .cloned()
                .map(|variant| (variant, ThresholdCountertermOrigin::Explicit))
                .collect(),
            None => vec![(
                ThresholdCountertermVariant {
                    name: Some("default".to_string()),
                    center_group: None,
                    subspace: None,
                    parent_lmb: None,
                    disable: false,
                    multiplier: None,
                },
                ThresholdCountertermOrigin::Autogenerated,
            )],
        }
    }

    fn preferred_amplitude_parent_lmb(&self) -> Result<LmbIndex> {
        let all_lmbs = self
            .derived_data
            .lmbs
            .as_ref()
            .expect("amplitude threshold resolution requires loop-momentum bases");
        if let Some((lmb_index, _)) = all_lmbs
            .iter_enumerated()
            .find(|(_, lmb)| lmb.loop_edges == self.graph.loop_momentum_basis.loop_edges)
        {
            return Ok(lmb_index);
        }

        all_lmbs
            .iter_enumerated()
            .next()
            .map(|(lmb_index, _)| lmb_index)
            .ok_or_else(|| eyre!("Graph '{}' has no generated LMBs", self.graph.name))
    }

    fn resolve_amplitude_threshold_variant_subspace(
        &self,
        variant: &ThresholdCountertermVariant,
        legacy_subspace: &SubspaceData,
        context: &str,
    ) -> Result<SubspaceData> {
        let all_lmbs = self
            .derived_data
            .lmbs
            .as_ref()
            .expect("amplitude threshold resolution requires loop-momentum bases");
        let containing_subgraph = self.graph.no_dummy();
        let build_in_parent = |parent_lmb_index| match &variant.subspace {
            Some(requested_basis_edges) => SubspaceData::new_from_parent_basis_edges(
                requested_basis_edges,
                &containing_subgraph,
                parent_lmb_index,
                &self.graph,
                all_lmbs,
            ),
            None => SubspaceData::new_with_user_selected_lmb(
                containing_subgraph.clone(),
                parent_lmb_index,
                &self.graph,
                all_lmbs,
            ),
        };

        if let Some(requested_parent_lmb) = &variant.parent_lmb {
            let matching_lmbs = all_lmbs
                .iter_enumerated()
                .filter_map(|(lmb_index, lmb)| {
                    lmb.loop_edges
                        .iter()
                        .eq(requested_parent_lmb.iter())
                        .then_some(lmb_index)
                })
                .collect_vec();
            if matching_lmbs.is_empty() {
                return Err(eyre!(
                    "{context} requests parent_lmb {:?}, which is not among graph '{}' generated LMBs",
                    requested_parent_lmb,
                    self.graph.name,
                ));
            }
            let mut built = matching_lmbs
                .iter()
                .map(|&lmb_index| build_in_parent(lmb_index).map(|subspace| (lmb_index, subspace)))
                .collect::<Result<Vec<_>>>()
                .with_context(|| {
                    format!("{context} is incompatible with its requested parent_lmb")
                })?;
            let (_, selected) = built.remove(0);
            if built
                .iter()
                .all(|(_, candidate)| selected.has_equivalent_embedding(candidate, all_lmbs))
            {
                return Ok(selected);
            }
            return Err(eyre!(
                "{context} parent_lmb {:?} has genuinely different generated embeddings",
                requested_parent_lmb,
            ));
        }

        if variant.subspace.is_none() {
            return Ok(legacy_subspace.clone());
        }

        let preferred_parent = legacy_subspace.parent_lmb_index();
        if let Ok(subspace) = build_in_parent(preferred_parent) {
            return Ok(subspace);
        }

        let mut compatible = Vec::new();
        let mut rejections = Vec::new();
        for (lmb_index, lmb) in all_lmbs.iter_enumerated() {
            if lmb_index == preferred_parent {
                continue;
            }
            match build_in_parent(lmb_index) {
                Ok(subspace) => compatible.push((lmb_index, subspace)),
                Err(error) => rejections.push(format!("parent {:?}: {error:#}", lmb.loop_edges,)),
            }
        }

        match compatible.len() {
            1 => Ok(compatible.pop().unwrap().1),
            0 => Err(eyre!(
                "{context} cannot resolve subspace {:?} in any generated parent LMB. Rejections:\n{}",
                variant.subspace,
                rejections.join("\n"),
            )),
            _ => {
                let (_, selected) = compatible.remove(0);
                if compatible
                    .iter()
                    .all(|(_, candidate)| selected.has_equivalent_embedding(candidate, all_lmbs))
                {
                    Ok(selected)
                } else {
                    Err(eyre!(
                        "{context} subspace {:?} has multiple genuinely different non-preferred parent LMB embeddings {:?}; specify parent_lmb",
                        variant.subspace,
                        compatible
                            .iter()
                            .map(|(lmb_index, _)| &all_lmbs[*lmb_index].loop_edges)
                            .collect_vec(),
                    ))
                }
            }
        }
    }

    fn resolve_amplitude_threshold_counterterm_directives(
        &self,
        raised_data: &RaisedEsurfaceData,
        settings: &GenerationSettings,
        locked_runtime_settings: &LockedRuntimeSettings,
        model: &Model,
    ) -> Result<(
        ResolvedThresholdCounterterms,
        TiVec<EsurfaceID, RaisedEsurfaceId>,
    )> {
        let global_cff = self
            .derived_data
            .cff_expression
            .as_ref()
            .expect("cff_expression should have been created");
        let all_lmbs = self
            .derived_data
            .lmbs
            .as_ref()
            .expect("amplitude threshold resolution requires loop-momentum bases");
        let preferred_parent = self.preferred_amplitude_parent_lmb()?;
        let legacy_subspace = SubspaceData::new_with_user_selected_lmb(
            self.graph.no_dummy(),
            preferred_parent,
            &self.graph,
            all_lmbs,
        )
        .with_context(|| {
            format!(
                "Graph '{}' cannot construct its legacy maximal amplitude threshold subspace",
                self.graph.name,
            )
        })?;

        let mut raised_esurface_ids: TiVec<EsurfaceID, Option<RaisedEsurfaceId>> =
            ti_vec![None; global_cff.surfaces.esurface_cache.len()];
        for (raised_esurface_id, raised_group) in raised_data.raised_groups.iter_enumerated() {
            for &esurface_id in &raised_group.esurface_ids {
                raised_esurface_ids[esurface_id] = Some(raised_esurface_id);
            }
        }
        let raised_esurface_ids = raised_esurface_ids
            .into_iter()
            .enumerate()
            .map(|(esurface_id, raised_esurface_id)| {
                raised_esurface_id.ok_or_else(|| {
                    eyre!(
                        "Graph '{}' E-surface {esurface_id} is missing from raised-esurface data",
                        self.graph.name,
                    )
                })
            })
            .collect::<Result<TiVec<EsurfaceID, RaisedEsurfaceId>>>()?;

        let external_filter: SuBitGraph = self.graph.external_filter();
        let mut incoming_externals = Vec::new();
        let mut outgoing_externals = Vec::new();
        for (edge, edge_id, _) in self.graph.iter_edges_of(&external_filter) {
            match edge {
                HedgePair::Unpaired {
                    flow: Flow::Sink, ..
                } => incoming_externals.push(edge_id),
                HedgePair::Unpaired {
                    flow: Flow::Source, ..
                } => outgoing_externals.push(edge_id),
                _ => unreachable!("the external filter must contain only unpaired edges"),
            }
        }
        let masses = settings
            .threshold_subtraction
            .check_esurface_at_generation
            .then(|| self.graph.get_real_mass_vector(model));
        let selected_esurface_ids = global_cff
            .orientations
            .iter()
            .flat_map(|orientation| {
                orientation
                    .expression
                    .iter_nodes()
                    .filter_map(|node| match node.data {
                        crate::cff::surface::HybridSurfaceID::Esurface(esurface_id) => {
                            Some(esurface_id)
                        }
                        _ => None,
                    })
            })
            .collect::<AHashSet<_>>();

        let configured_thresholds = self
            .graph
            .threshold_counterterms
            .cuts
            .iter()
            .flat_map(|cut| cut.thresholds.iter().map(|threshold| &threshold.edges))
            .collect_vec();
        if !configured_thresholds.is_empty() {
            // A selected orientation need not contain every configured threshold. Validate the
            // declarations against topology-wide bond identities so those absent only from the
            // selected CFF remain dormant without accepting a valid-edge typo silently. Bonds
            // crossing an edge contracted by amplitude CFF generation cannot become E-surfaces.
            let contracted_edges = self
                .graph
                .iter_edges_of(
                    &self
                        .graph
                        .tree_edges
                        .subtract(&self.graph.initial_state_cut)
                        .subtract(&self.graph.external_filter::<SuBitGraph>()),
                )
                .map(|(_, edge_id, _)| edge_id)
                .collect::<AHashSet<_>>();
            let no_dummy = self.graph.no_dummy();
            let external_count = self
                .graph
                .iter_edges_of(&no_dummy)
                .filter(|(pair, _, _)| matches!(pair, HedgePair::Unpaired { .. }))
                .count();
            let topology_thresholds = self
                .graph
                .underlying
                .all_bonds_of(&no_dummy, &(1..))
                .into_iter()
                .filter_map(|bond| {
                    let mut edges = Vec::new();
                    let mut bond_external_count = 0;
                    for (pair, edge_id, _) in self.graph.iter_edges_of(&bond) {
                        match pair {
                            HedgePair::Split { .. } => edges.push(edge_id),
                            HedgePair::Unpaired { .. } => bond_external_count += 1,
                            HedgePair::Paired { .. } => {
                                unreachable!("a topology bond contains only boundary half-edges")
                            }
                        }
                    }
                    edges.sort_unstable();
                    (!edges.is_empty()
                        && !edges.iter().any(|edge| contracted_edges.contains(edge))
                        && bond_external_count > 0
                        && bond_external_count < external_count)
                        .then_some(edges)
                })
                .collect::<AHashSet<_>>();

            for threshold_edges in configured_thresholds {
                if !topology_thresholds.contains(threshold_edges) {
                    return Err(eyre!(
                        "Amplitude graph '{}' threshold_counterterms threshold {:?} does not match a topology-discovered E-surface",
                        self.graph.name,
                        threshold_edges,
                    ));
                }
            }
        }

        let mut group_drafts = Vec::new();
        let mut legacy_equivalent = true;
        for (raised_esurface_id, raised_group) in raised_data.raised_groups.iter_enumerated() {
            let selected_group_esurface_ids = raised_group
                .esurface_ids
                .iter()
                .copied()
                .filter(|esurface_id| selected_esurface_ids.contains(esurface_id))
                .collect_vec();
            let Some(&representative_esurface_id) = selected_group_esurface_ids.first() else {
                debug!(
                    "Leaving amplitude graph '{}' raised threshold group {} dormant because none of its E-surfaces occurs in the selected orientations",
                    self.graph.name, raised_esurface_id.0,
                );
                continue;
            };
            let representative_esurface =
                &global_cff.surfaces.esurface_cache[representative_esurface_id];
            if representative_esurface.external_shift.is_empty() {
                continue;
            }
            if let Some(masses) = &masses
                && !locked_runtime_settings.existence_check(
                    representative_esurface,
                    masses,
                    &self.graph.get_external_signature(),
                    &self.graph.loop_momentum_basis,
                    settings.threshold_subtraction.esurface_existence_threshold,
                )
            {
                continue;
            }
            if settings
                .threshold_subtraction
                .assume_positive_external_energies
                && masses.is_none()
                && !representative_esurface
                    .external_shift_is_strictly_negative_for_positive_energies(
                        &incoming_externals,
                        &outgoing_externals,
                    )
            {
                continue;
            }

            let mut associations = Vec::new();
            let mut representative_variants =
                None::<Vec<ResolvedAmplitudeThresholdCountertermDraft>>;
            for &esurface_id in &selected_group_esurface_ids {
                let mut threshold_edges = global_cff.surfaces.esurface_cache[esurface_id]
                    .energies
                    .iter()
                    .copied()
                    .collect_vec();
                threshold_edges.sort_unstable();
                let configured = Self::configured_amplitude_threshold_variants(
                    &self.graph.threshold_counterterms,
                    &threshold_edges,
                );
                let mut variants = Vec::with_capacity(configured.len());
                for (variant, origin) in configured {
                    if !variant.disable
                        && variant
                            .multiplier
                            .as_ref()
                            .is_some_and(|multiplier| multiplier.symmetrize)
                    {
                        unimplemented!(
                            "symmetrized threshold-counterterm multipliers are not implemented (amplitude graph '{}', threshold {:?}, variant '{}')",
                            self.graph.name,
                            threshold_edges,
                            variant.name.as_deref().unwrap_or("default"),
                        );
                    }
                    let name = variant
                        .name
                        .clone()
                        .unwrap_or_else(|| "default".to_string());
                    let context = format!(
                        "amplitude graph '{}' threshold {:?} variant '{}'",
                        self.graph.name, threshold_edges, name,
                    );
                    let subspace = self.resolve_amplitude_threshold_variant_subspace(
                        &variant,
                        &legacy_subspace,
                        &context,
                    )?;
                    variants.push(ResolvedAmplitudeThresholdCountertermDraft {
                        name,
                        origin,
                        disable: variant.disable,
                        requested_subspace: variant.subspace,
                        requested_parent_lmb: variant.parent_lmb,
                        multiplier: variant.multiplier,
                        subspace,
                    });
                }
                let association_origin = variants
                    .first()
                    .map(|variant| variant.origin)
                    .unwrap_or(ThresholdCountertermOrigin::Autogenerated);
                let variant_subspaces = variants
                    .iter()
                    .map(|variant| variant.subspace.clone())
                    .collect();

                if let Some(expected) = &representative_variants {
                    if variants.len() != expected.len()
                        || variants.iter().zip(expected).any(|(variant, expected)| {
                            !variant.is_compatible_with(expected, all_lmbs)
                        })
                    {
                        return Err(eyre!(
                            "Amplitude graph '{}' raised threshold {:?} resolves incompatible variants for E-surfaces {} and {}",
                            self.graph.name,
                            raised_group.esurface_ids,
                            representative_esurface_id.0,
                            esurface_id.0,
                        ));
                    }
                } else {
                    representative_variants = Some(variants);
                }
                associations.push(ResolvedAmplitudeThresholdAssociationDraft {
                    esurface_id,
                    threshold_edges,
                    origin: association_origin,
                    variant_subspaces,
                });
            }

            let variants = representative_variants.unwrap_or_default();
            if variants.len() != 1
                || variants[0].disable
                || variants[0].multiplier.is_some()
                || !variants[0]
                    .subspace
                    .has_equivalent_embedding(&legacy_subspace, all_lmbs)
            {
                legacy_equivalent = false;
            }
            group_drafts.push(ResolvedAmplitudeThresholdGroupDraft {
                raised_esurface_group: raised_group.clone(),
                associations,
                variants,
            });
        }

        let mut variants =
            TiVec::<ThresholdCountertermVariantId, ResolvedThresholdCountertermVariant>::new();
        for group in group_drafts {
            for (variant_index, variant) in group.variants.into_iter().enumerate() {
                if variant.disable {
                    continue;
                }
                let parent_lmb_index = variant.subspace.parent_lmb_index();
                variants.push(ResolvedThresholdCountertermVariant {
                    name: variant.name,
                    center_group: None,
                    cut_group_id: None,
                    associations: group
                        .associations
                        .iter()
                        .map(|association| {
                            let subspace = association.variant_subspaces[variant_index].clone();
                            ResolvedThresholdCountertermAssociation {
                                cut_id: None,
                                cut_edges: Vec::new(),
                                threshold_edges: association.threshold_edges.clone(),
                                esurface_id: association.esurface_id,
                                requires_explicit_parent_lmb: subspace.parent_lmb_index()
                                    != legacy_subspace.parent_lmb_index(),
                                subspace,
                                eligible: true,
                                origin: association.origin,
                            }
                        })
                        .collect(),
                    side: ThresholdCountertermSide::Amplitude,
                    threshold_esurface_ids: group.raised_esurface_group.esurface_ids.clone(),
                    raised_esurface_group: group.raised_esurface_group.clone(),
                    requested_subspace: variant.requested_subspace,
                    requested_parent_lmb: variant.requested_parent_lmb,
                    resolved_parent_lmb: all_lmbs[parent_lmb_index].loop_edges.clone().into(),
                    subspace_loop_count: variant.subspace.loopcount(),
                    subspace: variant.subspace,
                    multiplier: variant.multiplier,
                });
            }
        }

        Ok((
            ResolvedThresholdCounterterms {
                legacy_equivalent,
                variants,
                cross_section_cut_groups: TiVec::new(),
                center_groups: Vec::new(),
            },
            raised_esurface_ids,
        ))
    }

    #[instrument(skip_all, err)]
    fn build_threshold_counterterm_parametric_integrand(
        &mut self,
        raised_data: &RaisedEsurfaceData,
        settings: &GenerationSettings,
        vakint: &Vakint,
        locked_runtime_settings: &LockedRuntimeSettings,
        model: &Model,
    ) -> Result<AmplitudeThresholdCountertermBuild> {
        let _progress_guard =
            generation_progress::enter_detailed_progress_span("Building Threshold Counterterms");
        let valid_orientations: Vec<_> = self
            .derived_data
            .cff_expression
            .as_ref()
            .expect("cff_expression should have been created")
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect();
        let (resolved, raised_esurface_ids) = self
            .resolve_amplitude_threshold_counterterm_directives(
                raised_data,
                settings,
                locked_runtime_settings,
                model,
            )?;
        let mut cuts = Vec::with_capacity(resolved.variants.len());
        for variant in &resolved.variants {
            let mut cut_union: SuBitGraph = self.graph.empty_subgraph();
            let representative_esurface = &self
                .derived_data
                .cff_expression
                .as_ref()
                .expect("cff_expression should have been created")
                .surfaces
                .esurface_cache[variant.raised_esurface_group.esurface_ids[0]];
            for energy in &representative_esurface.energies {
                let (_, hedge_pair) = self.graph[energy];
                match hedge_pair {
                    HedgePair::Paired { source, sink } => {
                        cut_union.add(source);
                        cut_union.add(sink);
                    }
                    _ => unreachable!(),
                }
            }

            let cutset = CutSet {
                residue_selector: ResidueSelector {
                    lu_cut: None,
                    left_th_cut: Some(variant.raised_esurface_group.clone()),
                    right_th_cut: None,
                },
                union: cut_union,
                canonicalize_external_shifts: false,
            };

            cuts.push(cutset);
        }
        let cut_structure = CutStructure { cuts };
        let mut exprs = if resolved.variants.is_empty() {
            Vec::new().into_iter()
        } else {
            settings
                .uv
                .orchestrator
                .parametric_integrands(
                    &mut self.graph,
                    cut_structure,
                    vakint,
                    OrientationProjection::new(&valid_orientations, &settings.orientation_pattern),
                    &settings.uv,
                )?
                .into_iter()
        };

        let mut variants =
            TiVec::<ThresholdCountertermVariantId, AmplitudeThresholdCountertermVariant>::new();
        for (variant_id, variant) in resolved.variants.iter_enumerated() {
            let expr = exprs.next().ok_or_else(|| {
                eyre!(
                    "Threshold orchestrator returned too few amplitude counterterms for graph '{}'",
                    self.graph.name,
                )
            })?;
            let jacobian_factor =
                Atom::var(GS.radius_star_left).pow(variant.subspace_loop_count as i32 * 3 - 1);
            let expr = expr.map(|integrand| integrand * &jacobian_factor);
            let counterterm_atom = AmplitudeCountertermAtom {
                parametric: expr.integrands,
            };
            let raised_group = expr.cuts.residue_selector.left_th_cut.ok_or_else(|| {
                eyre!(
                    "Threshold orchestrator amplitude result {} for graph '{}' has no threshold residue selector",
                    variant_id.0,
                    self.graph.name,
                )
            })?;
            if raised_group != variant.raised_esurface_group {
                return Err(eyre!(
                    "Threshold orchestrator amplitude result {} for graph '{}' returned raised group {:?}, expected {:?}",
                    variant_id.0,
                    self.graph.name,
                    raised_group.esurface_ids,
                    variant.raised_esurface_group.esurface_ids,
                ));
            }
            let raised_esurface_id = raised_esurface_ids[raised_group.esurface_ids[0]];
            debug!("raised_esurface_id: {}", raised_esurface_id.0);

            for (_, integrand) in counterterm_atom.parametric.iter() {
                debug!("counterterm integrand: {}", integrand.log_print(Some(100)));
            }

            variants.push(AmplitudeThresholdCountertermVariant {
                raised_esurface_id,
                atom: counterterm_atom,
            });
        }
        if exprs.next().is_some() {
            return Err(eyre!(
                "Threshold orchestrator returned too many amplitude counterterms for graph '{}'",
                self.graph.name,
            ));
        }

        let legacy_counterterms = if resolved.legacy_equivalent {
            let mut legacy_counterterms = ti_vec![
                AmplitudeCountertermAtom::new();
                raised_data.raised_groups.len()
            ];
            for variant in &variants {
                if legacy_counterterms[variant.raised_esurface_id].is_generated() {
                    return Err(eyre!(
                        "Legacy-equivalent amplitude graph '{}' generated duplicate threshold variants for raised group {}",
                        self.graph.name,
                        variant.raised_esurface_id.0,
                    ));
                }
                legacy_counterterms[variant.raised_esurface_id] = variant.atom.clone();
            }
            legacy_counterterms
        } else {
            // The generalized runtime consumes `variants` directly and deliberately leaves the
            // homogeneous raised-surface lane empty so duplicate geometry keeps independent IDs.
            TiVec::new()
        };

        Ok(AmplitudeThresholdCountertermBuild {
            legacy_counterterms,
            variants,
            resolved,
            raised_esurface_ids,
        })
    }

    #[instrument(skip_all)]
    fn build_lmbs(&mut self) {
        let _progress_guard =
            generation_progress::enter_detailed_progress_span("Building Loop Momentum Bases");
        let lmbs = self
            .graph
            .generate_loop_momentum_bases_of(&self.graph.no_dummy());

        self.derived_data.lmbs = Some(lmbs)
    }

    #[instrument(skip_all, err)]
    fn build_tropical_sampler(&mut self, process_settings: &GenerationSettings) -> Result<()> {
        let _progress_guard =
            generation_progress::enter_detailed_progress_span("Building Tropical Sampler");
        if process_settings
            .tropical_subgraph_table
            .disable_tropical_generation
        {
            debug!("Tropical subgraph table generation is disabled.");
            return Ok(());
        }
        let num_virtual_loop_edges = self.graph.iter_loop_edges().count();

        if num_virtual_loop_edges == 0 {
            debug!("Graph has no loop edges, skipping tropical sampler generation.");
            return Ok(());
        }

        let num_loops = self.graph.loop_momentum_basis.loop_edges.len();
        let target_omega = process_settings.tropical_subgraph_table.target_omega;

        let weight = (target_omega + (3 * num_loops) as f64 / 2.) / num_virtual_loop_edges as f64;

        debug!(
            "Building tropical subgraph table with all edge weights set to: {}",
            weight
        );

        let tropical_edges = self
            .graph
            .iter_loop_edges()
            .map(|(pair, _edge_id, edge)| {
                let is_massive = edge.data.particle.is_massive();

                let vertices = match pair {
                    HedgePair::Paired { source, sink } => (
                        self.graph.underlying.node_id(source).0 as u8,
                        self.graph.underlying.node_id(sink).0 as u8,
                    ),
                    _ => unreachable!(),
                };

                momtrop::Edge {
                    is_massive,
                    weight,
                    vertices,
                }
            })
            .collect_vec();

        let mut external_vertices_pool = AHashSet::new();

        for (pair, _, _) in self.graph.iter_non_loop_edges() {
            match pair {
                HedgePair::Paired { source, sink } => {
                    let source_id = self.graph.underlying.node_id(source).0 as u8;
                    let sink_id = self.graph.underlying.node_id(sink).0 as u8;

                    external_vertices_pool.insert(source_id);
                    external_vertices_pool.insert(sink_id);
                }
                HedgePair::Unpaired { hedge, .. } => {
                    let id = self.graph.underlying.node_id(hedge).0 as u8;
                    external_vertices_pool.insert(id);
                }
                _ => unreachable!(),
            }
        }

        let mut external_vertices = vec![];

        for tropical_edge in &tropical_edges {
            if external_vertices_pool.contains(&tropical_edge.vertices.0) {
                external_vertices.push(tropical_edge.vertices.0);
            }

            if external_vertices_pool.contains(&tropical_edge.vertices.1) {
                external_vertices.push(tropical_edge.vertices.1);
            }
        }

        let tropical_graph = momtrop::Graph {
            edges: tropical_edges,
            externals: external_vertices,
        };

        let loop_part = self
            .graph
            .iter_loop_edges()
            .map(|(_, edge_id, _edge)| {
                self.graph.loop_momentum_basis.edge_signatures[edge_id]
                    .internal
                    .clone()
                    .to_momtrop_format()
            })
            .collect_vec();

        let sampler = tropical_graph
            .build_sampler(loop_part)
            .map_err(|e| eyre!(e))?;

        let _: () = self.derived_data.tropical_sampler = Some(sampler);
        Ok(())
    }

    // Expects cff_expression, esurface_data,
    #[instrument(
          name = "generate_term_for_graph",
          level = "info",
          skip(self, model, global_settings),
          fields(
              graph.name = %self.graph.name

          ),
          err
      )]
    fn generate_term_for_graph(
        &self,
        model: &Model,
        own_group_position: GraphGroupPosition,
        esurface_map: TiVec<GroupEsurfaceId, TiVec<GraphGroupPosition, Option<RaisedEsurfaceId>>>,
        global_settings: &GlobalSettings,
    ) -> Result<(AmplitudeGraphTerm, GraphGenerationStats)> {
        let _progress_guard = generation_progress::enter_detailed_progress_span(&format!(
            "Generating Evaluators for {}",
            self.graph.name
        ));
        AmplitudeGraphTerm::from_amplitude_graph(
            self,
            own_group_position,
            esurface_map,
            model,
            global_settings,
        )
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeThresholdCountertermVariant {
    pub raised_esurface_id: RaisedEsurfaceId,
    pub atom: AmplitudeCountertermAtom,
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct AmplitudeDerivedData {
    pub all_mighty_integrand: Atom,
    /// Compatibility storage used by the current homogeneous amplitude runtime.
    pub threshold_counterterms: TiVec<RaisedEsurfaceId, AmplitudeCountertermAtom>,
    /// Canonical variant-indexed symbolic storage. Duplicate geometric thresholds remain distinct.
    pub threshold_counterterm_variants:
        TiVec<ThresholdCountertermVariantId, AmplitudeThresholdCountertermVariant>,
    pub resolved_threshold_counterterms: Option<ResolvedThresholdCounterterms>,
    pub raised_data: RaisedEsurfaceData,
    pub raised_esurface_ids: TiVec<EsurfaceID, RaisedEsurfaceId>,
    pub multi_channeling_setup: Option<LmbMultiChannelingSetup>,
    pub lmbs: Option<TiVec<LmbIndex, LoopMomentumBasis>>,
    pub tropical_sampler: Option<SampleGenerator<3>>,
    pub cff_expression: Option<CFFExpression<OrientationID>>,
}

pub trait AmplitudeState:
    Clone + std::fmt::Debug + bincode::Encode + for<'a> bincode::Decode<GammaLoopContextContainer<'a>>
{
}
impl AmplitudeState for () {}

#[derive(Clone, Encode, Decode, Debug)]
pub struct Processed {}
impl AmplitudeState for Processed {}

// #[derive(Clone, Encode, Decode, Debug)]
// pub struct ReadyForTerm {}
// impl AmplitudeState for ReadyForTerm {}

impl Amplitude {
    pub fn from_dot_string<Str: AsRef<str>>(s: Str, name: String, model: &Model) -> Result<Self> {
        let graphs = Graph::from_string(s, model)?;

        let mut amp = Amplitude::new(name);
        for g in graphs {
            amp.add_graph(g)?;
        }
        Ok(amp)
    }

    pub fn from_dot_file<P>(p: P, name: String, model: &Model) -> Result<Self>
    where
        P: AsRef<Path>,
    {
        let graphs = Graph::from_file(p, model)?;

        let mut amp = Amplitude::new(name);
        for g in graphs {
            amp.add_graph(g)?;
        }
        Ok(amp)
    }

    pub fn from_graph_list(name: impl ToString, mut graphs: Vec<Graph>) -> Result<Self> {
        let mut amplitude: Amplitude = Amplitude::new(name);
        amplitude.graph_group_structure = complete_group_parsing(&mut graphs)?;

        for amplitude_graph in graphs {
            amplitude.add_graph(amplitude_graph)?;
        }
        Ok(amplitude)
    }

    pub(crate) fn new(name: impl ToString) -> Self {
        Self {
            integrand: None,
            name: name.to_string(),
            graphs: vec![],
            graph_group_structure: TiVec::new(),
            external_particles: vec![],
            external_signature: SignatureLike::from_iter(iter::empty::<i8>()),
            group_derived_data: TiVec::new(),
        }
    }

    fn add_graph(&mut self, graph: Graph) -> Result<()> {
        let new_external_particels = graph.get_external_partcles();
        let new_external_signature = graph.get_external_signature();

        if !self.graphs.is_empty() {
            if self.external_particles != new_external_particels {
                return Err(eyre!("amplitude graph has different externals")).with_context(|| {
                    format!(
                        "Found {} externals, expected {} for the graph {}",
                        new_external_particels.len(),
                        self.external_particles.len(),
                        DotGraph::from(&graph).debug_dot()
                    )
                });
            }

            if self.external_signature != new_external_signature {
                return Err(eyre!("wrong external signature"));
            }
        } else {
            self.external_particles = new_external_particels;
            self.external_signature = new_external_signature;
        }

        self.graphs.push(AmplitudeGraph::new(graph));

        //  TODO: validate that the graph is compatible
        Ok(())
    }
}

fn threshold_counterterm_helper_atoms(
    order: u8,
    loop_number: usize,
) -> (SingleThresholdPieces<Atom>, Atom) {
    let loop_3 = loop_number as i64 * 3;

    let laurent_coeff_indices = (1..=order).map(|i| -(i as i8));

    let mut laurent_coeffs = laurent_coeff_indices.map(|laurent_coeff_index| {
        build_derivative_structure_atom(order, laurent_coeff_index)
            .replace(GS.rescale_star)
            .with(GS.radius_star_left)
    });

    let i = Atom::i();

    let radius = Atom::var(GS.radius_left);
    let radius_star = Atom::var(GS.radius_star_left);
    let uv_damp_plus = Atom::var(GS.uv_damp_plus_left);
    let uv_damp_minus = Atom::var(GS.uv_damp_minus_left);
    let hfunction = Atom::var(GS.hfunction_left_th);

    let delta_r_plus = &radius - &radius_star;
    let delta_r_minus = -&radius - &radius_star;

    let jacobian_ratio = (Atom::one() / &radius).pow(loop_3 - 1);

    let local_prefactor =
        &jacobian_ratio * (uv_damp_plus / &delta_r_plus + uv_damp_minus / &delta_r_minus);

    let integrated_prefactor = -i * Atom::var(GS.pi) * &jacobian_ratio * hfunction;

    let leading_laurent_coeff = laurent_coeffs.next().unwrap();
    let mut raised_local = Atom::Zero;

    for pow in 2..=order {
        raised_local += laurent_coeffs.next().unwrap()
            * &jacobian_ratio
            * (Atom::one() / delta_r_plus.pow(pow as i64)
                + Atom::one() / delta_r_minus.pow(pow as i64));
    }

    let local = local_prefactor.clone() * &leading_laurent_coeff + &raised_local;
    let integrated = integrated_prefactor.clone() * &leading_laurent_coeff;
    // Keep the no-directive helper structurally identical to the historical single-output lane.
    let legacy = (local_prefactor + integrated_prefactor) * leading_laurent_coeff + raised_local;

    debug!(
        "Threshold counterterm helper atoms for order {} and loop number {}: local={}, integrated={}, legacy={}",
        order, loop_number, local, integrated, legacy,
    );
    (SingleThresholdPieces { local, integrated }, legacy)
}

enum ThresholdCountertermHelperOutputs {
    Legacy,
    Pieces,
    LegacyAndPieces,
}

fn build_threshold_counterterm_helper(
    order: u8,
    loop_number: usize,
    evaluator_settings: &EvaluatorSettings,
    outputs: ThresholdCountertermHelperOutputs,
) -> GenericEvaluator {
    let (pieces, legacy) = threshold_counterterm_helper_atoms(order, loop_number);
    let atoms = match outputs {
        ThresholdCountertermHelperOutputs::Legacy => vec![legacy],
        ThresholdCountertermHelperOutputs::Pieces => vec![pieces.local, pieces.integrated],
        ThresholdCountertermHelperOutputs::LegacyAndPieces => {
            vec![legacy, pieces.local, pieces.integrated]
        }
    };
    let mut fn_map = FunctionMap::default();
    fn_map
        .add_aliases([(
            GS.pi.into(),
            Atom::num(Rational::try_from(std::f64::consts::PI).unwrap()),
        )])
        .unwrap();

    let mut params = params_for_derivative_order(order)
        .into_iter()
        .map(|param| param.replace(GS.rescale_star).with(GS.radius_star_left))
        .collect_vec();

    let radius = Atom::var(GS.radius_left);
    let radius_star = Atom::var(GS.radius_star_left);
    let uv_damp_plus = Atom::var(GS.uv_damp_plus_left);
    let uv_damp_minus = Atom::var(GS.uv_damp_minus_left);
    let hfunction = Atom::var(GS.hfunction_left_th);

    params.push(radius);
    params.push(radius_star);
    params.push(uv_damp_plus);
    params.push(uv_damp_minus);
    params.push(hfunction);

    GenericEvaluator::new_from_raw_params(
        atoms,
        &params,
        &fn_map,
        vec![],
        evaluator_settings.optimization_settings(),
        None,
        evaluator_settings,
    )
    .unwrap()
    .into_eager_only()
}

pub(crate) fn threshold_counterterm_helper(
    order: u8,
    loop_number: usize,
    evaluator_settings: &EvaluatorSettings,
) -> GenericEvaluator {
    build_threshold_counterterm_helper(
        order,
        loop_number,
        evaluator_settings,
        ThresholdCountertermHelperOutputs::Legacy,
    )
}

pub(crate) fn threshold_counterterm_pieces_helper(
    order: u8,
    loop_number: usize,
    evaluator_settings: &EvaluatorSettings,
) -> GenericEvaluator {
    build_threshold_counterterm_helper(
        order,
        loop_number,
        evaluator_settings,
        ThresholdCountertermHelperOutputs::Pieces,
    )
}

pub(crate) fn threshold_counterterm_recording_helper(
    order: u8,
    loop_number: usize,
    evaluator_settings: &EvaluatorSettings,
) -> GenericEvaluator {
    build_threshold_counterterm_helper(
        order,
        loop_number,
        evaluator_settings,
        ThresholdCountertermHelperOutputs::LegacyAndPieces,
    )
}

#[cfg(test)]
pub mod test {

    use std::{
        fs,
        io::Cursor,
        time::{SystemTime, UNIX_EPOCH},
    };

    use crate::{
        DependentMomentaConstructor, GammaLoopContextContainer,
        cff::{esurface::ExistingEsurfaceId, expression::OrientationID},
        dot,
        graph::{
            FeynmanGraph, Graph, GraphGroupPosition,
            autogen::Autogen,
            parse::IntoGraph,
            threshold_counterterms::{
                ThresholdCountertermCut, ThresholdCountertermMultiplier, ThresholdCountertermSpec,
                ThresholdCountertermThreshold, ThresholdCountertermVariant,
            },
        },
        initialisation::test_initialise,
        integrands::process::{
            MomentumSpaceEvaluationInput, ProcessIntegrand, amplitude::AmplitudeGraphTerm,
        },
        momentum::{
            Dep, ExternalMomenta, Helicity, ThreeMomentum,
            sample::{LoopIndex, MomentumSample},
        },
        processes::{
            AmplitudeGraph, DotExportSettings, ThresholdCountertermComponentKind,
            ThresholdCountertermOrigin, ThresholdCountertermSide,
        },
        settings::{
            GlobalSettings, RuntimeSettings,
            global::{GenerationSettings, OrientationPattern, ThresholdSubtractionSettings},
            runtime::kinematic::{
                Externals, KinematicsSettings, improvement::PhaseSpaceImprovementSettings,
            },
        },
        subtraction::amplitude_counterterm::AmplitudeCountertermComponentEvaluation,
        utils::{ArbPrec, F, load_generic_model},
    };
    use itertools::Itertools;
    use spenso::algebra::complex::Complex;
    use symbolica::state::State;
    use typed_index_collections::TiVec;

    #[test]
    fn amplitude_directive_defaults_and_explicit_variants_use_the_empty_cut() {
        let threshold_edges = vec![super::EdgeIndex::from(0), super::EdgeIndex::from(1)];
        let mut spec = ThresholdCountertermSpec {
            schema_version: 1,
            cuts: vec![ThresholdCountertermCut {
                edges: Vec::new(),
                thresholds: vec![ThresholdCountertermThreshold {
                    edges: threshold_edges.clone(),
                    counterterms: Vec::new(),
                }],
            }],
        };

        let default =
            super::AmplitudeGraph::configured_amplitude_threshold_variants(&spec, &threshold_edges);
        assert_eq!(default.len(), 1);
        assert_eq!(default[0].0.name.as_deref(), Some("default"));
        assert_eq!(default[0].1, ThresholdCountertermOrigin::Autogenerated);

        spec.cuts[0].thresholds[0].counterterms = vec![
            ThresholdCountertermVariant {
                name: Some("disabled".to_string()),
                center_group: None,
                subspace: None,
                parent_lmb: None,
                disable: true,
                multiplier: None,
            },
            ThresholdCountertermVariant {
                name: Some("duplicate".to_string()),
                center_group: None,
                subspace: None,
                parent_lmb: None,
                disable: false,
                multiplier: None,
            },
        ];
        let explicit =
            super::AmplitudeGraph::configured_amplitude_threshold_variants(&spec, &threshold_edges);
        assert_eq!(explicit.len(), 2);
        assert!(explicit[0].0.disable);
        assert!(
            explicit
                .iter()
                .all(|(_, origin)| *origin == ThresholdCountertermOrigin::Explicit)
        );
    }

    #[test]
    fn amplitude_legacy_parent_falls_back_when_the_generation_lmb_is_not_generated() {
        test_initialise().unwrap();
        let model = load_generic_model("scalars");
        let mut graph: AmplitudeGraph =
            include_str!("../../../../tests/resources/graphs/uv_tests/dotted_sunrise.dot")
                .into_graph(&model)
                .unwrap();
        graph.build_lmbs();

        let all_lmbs = graph.derived_data.lmbs.as_ref().unwrap();
        let generation_lmb = graph.graph.loop_momentum_basis.loop_edges.clone();
        assert!(
            all_lmbs.iter().all(|lmb| lmb.loop_edges != generation_lmb),
            "the dotted-sunrise topology must exercise the legacy LMB fallback",
        );
        let fallback_parent = all_lmbs.iter_enumerated().next().unwrap().0;
        assert_eq!(
            graph.preferred_amplitude_parent_lmb().unwrap(),
            fallback_parent,
        );

        let legacy_subspace = super::SubspaceData::new_with_user_selected_lmb(
            graph.graph.no_dummy(),
            fallback_parent,
            &graph.graph,
            all_lmbs,
        )
        .unwrap();
        let error = graph
            .resolve_amplitude_threshold_variant_subspace(
                &ThresholdCountertermVariant {
                    name: Some("explicit_missing_parent".to_string()),
                    center_group: None,
                    subspace: None,
                    parent_lmb: Some(generation_lmb.iter().copied().collect()),
                    disable: false,
                    multiplier: None,
                },
                &legacy_subspace,
                "dotted-sunrise explicit-parent test",
            )
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("is not among graph 'sunrise' generated LMBs")
        );
    }

    #[test]
    fn amplitude_rejects_a_valid_edge_set_that_is_not_a_topological_threshold() {
        test_initialise().unwrap();
        let mut graph: AmplitudeGraph = dot!(
            digraph unmatched_amplitude_threshold {
                edge [particle=scalar_1]
                node [num=1]
                e [style=invis]
                e -> A:0 [id=3]
                B:1 -> e [id=2]
                A -> B [id=1]
                A -> B [id=0]
            },
            "scalars"
        )
        .unwrap();
        graph.graph.threshold_counterterms = Autogen::explicit(ThresholdCountertermSpec {
            schema_version: 1,
            cuts: vec![ThresholdCountertermCut {
                edges: Vec::new(),
                thresholds: vec![ThresholdCountertermThreshold {
                    // Edge 0 is internal and valid, but the parallel edge 1 means it is not a
                    // complete connected cut boundary and therefore cannot be an E-surface.
                    edges: vec![super::EdgeIndex(0)],
                    counterterms: Vec::new(),
                }],
            }],
        });

        let settings = GenerationSettings {
            threshold_subtraction: ThresholdSubtractionSettings {
                enable_thresholds: true,
                assume_positive_external_energies: false,
                check_esurface_at_generation: false,
                ..Default::default()
            },
            ..Default::default()
        };
        graph.generate_cff(&settings.orientation_pattern).unwrap();
        let raised_data = graph.graph.determine_raised_esurfaces_from_expression(
            graph.derived_data.cff_expression.as_ref().unwrap(),
        );
        graph.build_lmbs();

        let error = graph
            .resolve_amplitude_threshold_counterterm_directives(
                &raised_data,
                &settings,
                &(&RuntimeSettings::default()).into(),
                &load_generic_model("scalars"),
            )
            .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("does not match a topology-discovered E-surface")
        );
    }

    #[test]
    fn amplitude_orientation_excluded_topological_threshold_remains_dormant() {
        test_initialise().unwrap();
        let mut graph: AmplitudeGraph = dot!(
            digraph dormant_amplitude_threshold {
                edge [particle=scalar_1]
                node [num=1]
                e [style=invis]
                e -> A:0 [id=4]
                C:1 -> e [id=3]
                A -> B [id=0]
                B -> C [id=1]
                A -> C [id=2]
            },
            "scalars"
        )
        .unwrap();
        let selected_graph = graph.graph.clone();
        graph.generate_cff(&OrientationPattern::default()).unwrap();
        let full_cff = graph.derived_data.cff_expression.as_ref().unwrap();
        let physical_thresholds = full_cff
            .surfaces
            .esurface_cache
            .iter_enumerated()
            .filter(|(_, esurface)| !esurface.external_shift.is_empty())
            .map(|(esurface_id, esurface)| {
                (
                    esurface_id,
                    esurface.energies.iter().copied().sorted().collect_vec(),
                )
            })
            .collect_vec();
        let (orientation_pattern, dormant_threshold_edges) = full_cff
            .orientations
            .iter()
            .find_map(|orientation| {
                let present = orientation
                    .expression
                    .iter_nodes()
                    .filter_map(|node| match node.data {
                        crate::cff::surface::HybridSurfaceID::Esurface(esurface_id) => Some(
                            full_cff.surfaces.esurface_cache[esurface_id]
                                .energies
                                .iter()
                                .copied()
                                .sorted()
                                .collect_vec(),
                        ),
                        _ => None,
                    })
                    .collect::<std::collections::BTreeSet<_>>();
                physical_thresholds
                    .iter()
                    .find(|(_, edges)| !present.contains(edges))
                    .map(|(_, edges)| {
                        (
                            OrientationPattern::from_orientation(orientation),
                            edges.clone(),
                        )
                    })
            })
            .expect("the triangle must have a physical threshold absent from one orientation");

        let mut graph = AmplitudeGraph::new(selected_graph);
        graph.graph.threshold_counterterms = Autogen::explicit(ThresholdCountertermSpec {
            schema_version: 1,
            cuts: vec![ThresholdCountertermCut {
                edges: Vec::new(),
                thresholds: vec![ThresholdCountertermThreshold {
                    edges: dormant_threshold_edges.clone(),
                    counterterms: vec![ThresholdCountertermVariant {
                        name: Some("dormant".to_string()),
                        center_group: None,
                        subspace: None,
                        parent_lmb: None,
                        disable: false,
                        multiplier: None,
                    }],
                }],
            }],
        });
        let settings = GenerationSettings {
            orientation_pattern,
            threshold_subtraction: ThresholdSubtractionSettings {
                enable_thresholds: true,
                assume_positive_external_energies: false,
                check_esurface_at_generation: false,
                ..Default::default()
            },
            ..Default::default()
        };
        graph.generate_cff(&settings.orientation_pattern).unwrap();
        assert_eq!(
            graph
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .orientations
                .len(),
            1,
        );
        let raised_data = graph.graph.determine_raised_esurfaces_from_expression(
            graph.derived_data.cff_expression.as_ref().unwrap(),
        );
        graph.build_lmbs();
        let (resolved, _) = graph
            .resolve_amplitude_threshold_counterterm_directives(
                &raised_data,
                &settings,
                &(&RuntimeSettings::default()).into(),
                &load_generic_model("scalars"),
            )
            .unwrap();
        assert!(resolved.variants.iter().all(|variant| {
            variant.name != "dormant"
                && variant
                    .associations
                    .iter()
                    .all(|association| association.threshold_edges != dormant_threshold_edges)
        }));
    }

    #[test]
    fn amplitude_variant_subspace_resolves_in_the_requested_parent() {
        test_initialise().unwrap();
        let mut graph: AmplitudeGraph = dot!(
            digraph amplitude_variant_subspace {
                edge [particle=scalar_1]
                node [num=1]
                e [style=invis]
                e -> A:0 [id=3]
                B:1 -> e [id=4]
                A -> B [id=0]
                A -> B [id=1]
                A -> B [id=2]
            },
            "scalars"
        )
        .unwrap();
        graph.build_lmbs();
        let preferred_parent = graph.preferred_amplitude_parent_lmb().unwrap();
        let all_lmbs = graph.derived_data.lmbs.as_ref().unwrap();
        let parent_edges = all_lmbs[preferred_parent]
            .loop_edges
            .iter()
            .copied()
            .collect::<Vec<_>>();
        let requested_subspace = vec![parent_edges[0]];
        let legacy_subspace = super::SubspaceData::new_with_user_selected_lmb(
            graph.graph.no_dummy(),
            preferred_parent,
            &graph.graph,
            all_lmbs,
        )
        .unwrap();
        let resolved = graph
            .resolve_amplitude_threshold_variant_subspace(
                &ThresholdCountertermVariant {
                    name: Some("one_loop".to_string()),
                    center_group: None,
                    subspace: Some(requested_subspace.clone()),
                    parent_lmb: Some(parent_edges),
                    disable: false,
                    multiplier: None,
                },
                &legacy_subspace,
                "amplitude subspace test",
            )
            .unwrap();

        assert_eq!(resolved.parent_lmb_index(), preferred_parent);
        assert_eq!(resolved.loopcount(), 1);
        assert_eq!(
            resolved.iter_basis_edges(all_lmbs).collect::<Vec<_>>(),
            requested_subspace,
        );
    }

    #[test]
    fn amplitude_preprocessing_keeps_duplicate_geometric_variants_independent() {
        std::thread::Builder::new()
            .name("amplitude-threshold-variant-test".to_string())
            .stack_size(32 * 1024 * 1024)
            .spawn(|| {
                test_initialise().unwrap();
                let mut graph: AmplitudeGraph = dot!(
                    digraph duplicate_threshold_variants {
                        edge [particle=scalar_1]
                        node [num=1]
                        e [style=invis]
                        e -> A:0 [id=3]
                        B:1 -> e [id=2]
                        A -> B [id=1]
                        A -> B [id=0]
                    },
                    "scalars"
                )
                .unwrap();
                graph.generate_cff(&OrientationPattern::default()).unwrap();
                let threshold_edges = graph
                    .derived_data
                    .cff_expression
                    .as_ref()
                    .unwrap()
                    .surfaces
                    .esurface_cache
                    .iter()
                    .find(|esurface| !esurface.external_shift.is_empty())
                    .expect("the scalar bubble must contain a physical threshold")
                    .energies
                    .iter()
                    .copied()
                    .sorted()
                    .collect::<Vec<_>>();
                graph.graph.threshold_counterterms = Autogen::explicit(ThresholdCountertermSpec {
                    schema_version: 1,
                    cuts: vec![ThresholdCountertermCut {
                        edges: Vec::new(),
                        thresholds: vec![ThresholdCountertermThreshold {
                            edges: threshold_edges.clone(),
                            counterterms: vec![
                                ThresholdCountertermVariant {
                                    name: Some("first".to_string()),
                                    center_group: None,
                                    subspace: None,
                                    parent_lmb: None,
                                    disable: false,
                                    multiplier: None,
                                },
                                ThresholdCountertermVariant {
                                    name: Some("disabled".to_string()),
                                    center_group: None,
                                    subspace: None,
                                    parent_lmb: None,
                                    disable: true,
                                    multiplier: None,
                                },
                                ThresholdCountertermVariant {
                                    name: Some("second".to_string()),
                                    center_group: None,
                                    subspace: None,
                                    parent_lmb: None,
                                    disable: false,
                                    multiplier: None,
                                },
                            ],
                        }],
                    }],
                });

                let model = load_generic_model("scalars");
                let mut generation_settings = GenerationSettings {
                    threshold_subtraction: ThresholdSubtractionSettings {
                        enable_thresholds: true,
                        assume_positive_external_energies: false,
                        check_esurface_at_generation: false,
                        ..Default::default()
                    },
                    ..Default::default()
                };
                graph
                    .preprocess(
                        &model,
                        &generation_settings,
                        &(&RuntimeSettings::default()).into(),
                    )
                    .unwrap();

                let resolved = graph
                    .derived_data
                    .resolved_threshold_counterterms
                    .as_ref()
                    .unwrap();
                assert!(!resolved.legacy_equivalent);
                let mut names_by_group = std::collections::BTreeMap::<Vec<usize>, Vec<&str>>::new();
                for variant in resolved.variants.iter().filter(|variant| {
                    variant
                        .associations
                        .iter()
                        .any(|association| association.threshold_edges == threshold_edges)
                }) {
                    names_by_group
                        .entry(
                            variant
                                .threshold_esurface_ids
                                .iter()
                                .map(|esurface_id| esurface_id.0)
                                .collect(),
                        )
                        .or_default()
                        .push(&variant.name);
                }
                assert!(!names_by_group.is_empty());
                assert!(
                    names_by_group
                        .values()
                        .all(|names| names == &["first", "second"])
                );
                assert!(
                    resolved
                        .variants
                        .iter()
                        .all(|variant| variant.side == ThresholdCountertermSide::Amplitude)
                );
                assert_eq!(
                    graph.derived_data.threshold_counterterm_variants.len(),
                    resolved.variants.len(),
                );
                assert!(
                    graph
                        .derived_data
                        .threshold_counterterm_variants
                        .iter()
                        .all(|variant| variant.atom.is_generated())
                );
                assert!(graph.derived_data.threshold_counterterms.is_empty());

                let (term, _) = AmplitudeGraphTerm::from_amplitude_graph(
                    &graph,
                    GraphGroupPosition(0),
                    TiVec::new(),
                    &model,
                    &GlobalSettings {
                        generation: generation_settings.clone(),
                        ..Default::default()
                    },
                )
                .unwrap();
                let registry = term
                    .threshold_counterterm
                    .metadata_registry
                    .as_ref()
                    .expect("explicit duplicate variants must create amplitude metadata");
                assert_eq!(registry.variants.len(), resolved.variants.len());
                assert_eq!(registry.components.len(), resolved.variants.len() * 2);
                assert!(registry.evaluators.is_empty());
                assert!(term.threshold_counterterm.threshold_multipliers.is_none());

                let decomposition =
                    crate::integrands::process::amplitude::amplitude_threshold_event_info(
                        registry,
                        Complex::new_re(F(100.0)),
                        vec![
                            AmplitudeCountertermComponentEvaluation {
                                variant_id: super::ThresholdCountertermVariantId(0),
                                kind: ThresholdCountertermComponentKind::Local,
                                esurface_id: super::RaisedEsurfaceId(4),
                                overlap_group: 2,
                                multiplier_value: F(3.0),
                                bare: Some(Complex::new_re(F(-2.0))),
                                weighted: Complex::new_re(F(-6.0)),
                                evaluation_skipped: false,
                            },
                            AmplitudeCountertermComponentEvaluation {
                                variant_id: super::ThresholdCountertermVariantId(0),
                                kind: ThresholdCountertermComponentKind::Integrated,
                                esurface_id: super::RaisedEsurfaceId(4),
                                overlap_group: 2,
                                multiplier_value: F(0.0),
                                bare: None,
                                weighted: Complex::new_re(F(0.0)),
                                evaluation_skipped: true,
                            },
                        ],
                    )
                    .unwrap();
                assert_eq!(decomposition.total(), Complex::new_re(F(94.0)));
                assert_eq!(decomposition.components[0].component_id, 0);
                assert_eq!(decomposition.components[1].component_id, 1);
                assert!(decomposition.components[1].evaluation_skipped);

                graph.graph.threshold_counterterms =
                    Autogen::generated(ThresholdCountertermSpec::default());
                graph
                    .preprocess(
                        &model,
                        &generation_settings,
                        &(&RuntimeSettings::default()).into(),
                    )
                    .unwrap();
                let resolved = graph
                    .derived_data
                    .resolved_threshold_counterterms
                    .as_ref()
                    .unwrap();
                assert!(resolved.legacy_equivalent);
                assert_eq!(
                    graph.derived_data.threshold_counterterms.len(),
                    graph.derived_data.raised_data.raised_groups.len(),
                );
                assert_eq!(
                    graph
                        .derived_data
                        .threshold_counterterms
                        .iter()
                        .filter(|counterterm| counterterm.is_generated())
                        .count(),
                    graph.derived_data.threshold_counterterm_variants.len(),
                );

                generation_settings.threshold_subtraction.enable_thresholds = false;
                graph
                    .preprocess(
                        &model,
                        &generation_settings,
                        &(&RuntimeSettings::default()).into(),
                    )
                    .unwrap();
                assert!(graph.derived_data.resolved_threshold_counterterms.is_none());
                assert!(graph.derived_data.threshold_counterterms.is_empty());
                assert!(graph.derived_data.threshold_counterterm_variants.is_empty());
                assert!(graph.derived_data.raised_data.raised_groups.is_empty());
            })
            .expect("amplitude threshold-variant test thread must start")
            .join()
            .expect("amplitude threshold-variant test thread must finish successfully");
    }

    #[test]
    fn generalized_amplitude_approach_keeps_duplicate_variant_subspaces_and_complements() {
        std::thread::Builder::new()
            .name("generalized-amplitude-threshold-approach-test".to_string())
            .stack_size(32 * 1024 * 1024)
            .spawn(|| {
                test_initialise().unwrap();
                let model = load_generic_model("scalars");
                let mut graph: Graph = dot!(
                    digraph generalized_threshold_approach {
                        edge [particle=scalar_1]
                        node [num=1]
                        e [style=invis]
                        e -> A:0 [id=4]
                        B:1 -> e [id=3]
                        A -> B [id=2]
                        A -> B [id=1]
                        A -> B [id=0]
                    },
                    "scalars"
                )
                .unwrap();
                let mut parent_probe = AmplitudeGraph::new(graph.clone());
                parent_probe.build_lmbs();
                let preferred_parent = parent_probe.preferred_amplitude_parent_lmb().unwrap();
                let parent_edges = parent_probe.derived_data.lmbs.as_ref().unwrap()
                    [preferred_parent]
                    .loop_edges
                    .iter()
                    .copied()
                    .collect::<Vec<_>>();
                assert_eq!(parent_edges.len(), 2);
                graph.threshold_counterterms = Autogen::explicit(ThresholdCountertermSpec {
                    schema_version: 1,
                    cuts: vec![ThresholdCountertermCut {
                        edges: Vec::new(),
                        thresholds: vec![ThresholdCountertermThreshold {
                            edges: vec![
                                super::EdgeIndex(0),
                                super::EdgeIndex(1),
                                super::EdgeIndex(2),
                            ],
                            counterterms: vec![
                                ThresholdCountertermVariant {
                                    name: Some("one_loop".to_string()),
                                    center_group: None,
                                    subspace: Some(vec![parent_edges[0]]),
                                    parent_lmb: Some(parent_edges.clone()),
                                    disable: false,
                                    multiplier: None,
                                },
                                ThresholdCountertermVariant {
                                    name: Some("two_loop".to_string()),
                                    center_group: None,
                                    subspace: Some(parent_edges.clone()),
                                    parent_lmb: Some(parent_edges),
                                    disable: false,
                                    multiplier: None,
                                },
                            ],
                        }],
                    }],
                });

                let generation_settings = GenerationSettings {
                    threshold_subtraction: ThresholdSubtractionSettings {
                        enable_thresholds: true,
                        assume_positive_external_energies: false,
                        check_esurface_at_generation: false,
                        ..Default::default()
                    },
                    ..Default::default()
                };
                let runtime = RuntimeSettings {
                    kinematics: KinematicsSettings {
                        e_cm: 6.0,
                        externals: Externals::Constant {
                            momenta: vec![
                                ExternalMomenta::Independent([F(6.0), F(0.0), F(0.0), F(0.0)]),
                                ExternalMomenta::Dependent(Dep::Dep),
                            ],
                            helicities: vec![Helicity::ZERO; 2],
                            improvement_settings: PhaseSpaceImprovementSettings::default(),
                            f_64_cache: None,
                            f_128_cache: None,
                        },
                    },
                    ..RuntimeSettings::default()
                };
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(1)
                    .build()
                    .unwrap();
                let mut amplitude = super::Amplitude::from_graph_list(
                    "generalized_threshold_approach",
                    vec![graph],
                )
                .unwrap();
                amplitude
                    .preprocess(&model, &generation_settings, &(&runtime).into(), &pool)
                    .unwrap();
                amplitude
                    .build_integrand(
                        &model,
                        "generalized_threshold_approach",
                        &GlobalSettings {
                            generation: generation_settings,
                            ..Default::default()
                        },
                        (&runtime).into(),
                        &pool,
                    )
                    .unwrap();
                let ProcessIntegrand::Amplitude(integrand) = amplitude.integrand.as_mut().unwrap()
                else {
                    panic!("approach fixture built a non-amplitude integrand")
                };
                let external_signature =
                    integrand.data.graph_terms[0].graph.get_external_signature();
                let sample = MomentumSample::new(
                    vec![
                        ThreeMomentum::new(
                            F::<ArbPrec>::from_f64(0.4),
                            F::<ArbPrec>::from_f64(0.1),
                            F::<ArbPrec>::from_f64(-0.2),
                        ),
                        ThreeMomentum::new(
                            F::<ArbPrec>::from_f64(-0.3),
                            F::<ArbPrec>::from_f64(0.2),
                            F::<ArbPrec>::from_f64(0.15),
                        ),
                    ]
                    .into_iter()
                    .collect(),
                    0,
                    &runtime.kinematics.externals,
                    0,
                    F::<ArbPrec>::from_f64(1.0),
                    DependentMomentaConstructor::Amplitude(&external_signature),
                    None,
                )
                .unwrap();

                let term = &mut integrand.data.graph_terms[0];
                let approach = term
                    .kinematics_for_threshold_approach(&runtime, &model, &sample)
                    .unwrap();
                let variant_ids = approach
                    .variant_ids
                    .as_ref()
                    .expect("generalized approach must expose stable variant IDs");
                assert_eq!(variant_ids.len(), approach.existing_esurfaces.len());
                assert!(variant_ids.len() >= 2);
                assert_eq!(
                    variant_ids
                        .iter()
                        .copied()
                        .collect::<std::collections::BTreeSet<_>>()
                        .len(),
                    variant_ids.len(),
                );

                let counterterm = &term.threshold_counterterm;
                let mut loop_counts_by_raised =
                    std::collections::BTreeMap::<_, std::collections::BTreeSet<_>>::new();
                for &variant_id in variant_ids {
                    loop_counts_by_raised
                        .entry(counterterm.variant_raised_esurfaces[variant_id])
                        .or_default()
                        .insert(counterterm.variant_subspaces[variant_id].loopcount());
                }
                assert!(loop_counts_by_raised.values().any(|counts| {
                    counts
                        == &[1, 2]
                            .into_iter()
                            .collect::<std::collections::BTreeSet<_>>()
                }));

                let mut referenced_variants = std::collections::BTreeSet::new();
                for group in &approach.overlap_groups_with_kinematics {
                    let centers = group
                        .approach_centers_at_esurface
                        .as_ref()
                        .expect("generalized overlap groups need per-variant centers");
                    assert_eq!(centers.len(), group.overlap_group.existing_esurfaces.len());
                    assert_eq!(
                        group.loop_momenta_at_esurface.len(),
                        group.overlap_group.existing_esurfaces.len()
                    );
                    for (position, &existing_esurface_id) in
                        group.overlap_group.existing_esurfaces.iter().enumerate()
                    {
                        let variant_id = variant_ids[existing_esurface_id];
                        referenced_variants.insert(variant_id);
                        let local_id = ExistingEsurfaceId::from(position);
                        let root = group.loop_momenta_at_esurface[local_id]
                            .as_ref()
                            .expect("each generalized threshold needs an r_star sample");
                        let center = centers[local_id]
                            .as_ref()
                            .expect("each generalized threshold needs a full center");
                        for momentum in root.loop_moms().iter().chain(center.iter()) {
                            for component in [&momentum.px, &momentum.py, &momentum.pz] {
                                assert!(component.into_f64().is_finite());
                            }
                        }

                        let subspace = &counterterm.variant_subspaces[variant_id];
                        let parent_lmb = subspace.get_lmb(&counterterm.lmbs);
                        let base_in_parent =
                            sample.lmb_transform(&term.graph.loop_momentum_basis, parent_lmb);
                        let root_in_parent =
                            root.lmb_transform(&term.graph.loop_momentum_basis, parent_lmb);
                        let mut center_sample = root.clone();
                        center_sample.sample.loop_moms = center.clone();
                        let center_in_parent = center_sample
                            .lmb_transform(&term.graph.loop_momentum_basis, parent_lmb);
                        for loop_index in (0..term.graph.get_loop_number()).map(LoopIndex::from) {
                            if subspace.contains_loop_index(loop_index) {
                                continue;
                            }
                            for (actual, expected) in [
                                (
                                    &root_in_parent.loop_moms()[loop_index],
                                    &base_in_parent.loop_moms()[loop_index],
                                ),
                                (
                                    &center_in_parent.loop_moms()[loop_index],
                                    &base_in_parent.loop_moms()[loop_index],
                                ),
                            ] {
                                for (actual, expected) in [
                                    (&actual.px, &expected.px),
                                    (&actual.py, &expected.py),
                                    (&actual.pz, &expected.pz),
                                ] {
                                    assert!((actual - expected).abs().into_f64() < 1.0e-20);
                                }
                            }
                        }
                    }
                }
                assert_eq!(
                    referenced_variants,
                    variant_ids.iter().copied().collect(),
                    "overlap groups must retain every semantic duplicate independently",
                );
            })
            .expect("generalized amplitude approach test thread must start")
            .join()
            .expect("generalized amplitude approach test thread must finish successfully");
    }

    #[test]
    fn amplitude_dot_export_materializes_resolved_defaults_and_round_trips() {
        std::thread::Builder::new()
            .name("amplitude-threshold-dot-export-test".to_string())
            .stack_size(32 * 1024 * 1024)
            .spawn(|| {
                test_initialise().unwrap();
                let mut graph: AmplitudeGraph = dot!(
                    digraph amplitude_threshold_dot_export {
                        // A massive external scalar above the two-massless-particle threshold
                        // guarantees that the runtime part of this regression records L/I pieces.
                        edge [particle=scalar_0]
                        node [num=1]
                        e [style=invis]
                        e -> A:0 [id=3 particle=scalar_2]
                        B:1 -> e [id=2 particle=scalar_2]
                        A -> B [id=1]
                        A -> B [id=0]
                    },
                    "scalars"
                )
                .unwrap();
                let model = load_generic_model("scalars");
                let generation_settings = GenerationSettings {
                    threshold_subtraction: ThresholdSubtractionSettings {
                        enable_thresholds: true,
                        assume_positive_external_energies: false,
                        check_esurface_at_generation: false,
                        ..Default::default()
                    },
                    ..Default::default()
                };
                graph
                    .preprocess(
                        &model,
                        &generation_settings,
                        &(&RuntimeSettings::default()).into(),
                    )
                    .unwrap();

                assert!(graph.graph.threshold_counterterms.autogenerated);
                assert!(graph.graph.threshold_counterterms.cuts.is_empty());
                let resolved = graph
                    .derived_data
                    .resolved_threshold_counterterms
                    .as_ref()
                    .unwrap();
                assert!(resolved.legacy_equivalent);
                let expected_materialized = resolved.materialized_spec(
                    &graph.graph.threshold_counterterms,
                    graph.derived_data.lmbs.as_ref().unwrap(),
                );
                assert!(!expected_materialized.cuts.is_empty());
                assert!(
                    expected_materialized
                        .cuts
                        .iter()
                        .all(|cut| cut.edges.is_empty())
                );

                let mut ordinary = String::new();
                graph
                    .write_dot_fmt(&mut ordinary, &DotExportSettings::default())
                    .unwrap();
                assert!(!ordinary.contains("threshold_counterterms"));
                let ordinary_import = Graph::from_string(&ordinary, &model).unwrap().remove(0);
                assert!(ordinary_import.threshold_counterterms.autogenerated);
                assert!(ordinary_import.threshold_counterterms.cuts.is_empty());

                let mut normalized = String::new();
                graph
                    .write_dot_fmt(
                        &mut normalized,
                        &DotExportSettings {
                            include_autogenerated_fields: true,
                            ..DotExportSettings::default()
                        },
                    )
                    .unwrap();
                let normalized_import = Graph::from_string(&normalized, &model).unwrap().remove(0);
                assert!(!normalized_import.threshold_counterterms.autogenerated);
                assert_eq!(
                    *normalized_import.threshold_counterterms,
                    expected_materialized
                );
                assert!(graph.graph.threshold_counterterms.autogenerated);
                assert!(graph.graph.threshold_counterterms.cuts.is_empty());

                let mut round_tripped = AmplitudeGraph::new(normalized_import);
                round_tripped
                    .preprocess(
                        &model,
                        &generation_settings,
                        &(&RuntimeSettings::default()).into(),
                    )
                    .unwrap();
                assert!(
                    round_tripped
                        .derived_data
                        .resolved_threshold_counterterms
                        .as_ref()
                        .unwrap()
                        .legacy_equivalent
                );
                let (term, _) = AmplitudeGraphTerm::from_amplitude_graph(
                    &round_tripped,
                    GraphGroupPosition(0),
                    TiVec::new(),
                    &model,
                    &GlobalSettings {
                        generation: generation_settings.clone(),
                        ..Default::default()
                    },
                )
                .unwrap();
                assert!(term.threshold_counterterm.legacy_equivalent);
                assert!(term.threshold_counterterm.variant_evaluators.is_empty());
                assert!(term.threshold_counterterm.variant_subspaces.is_empty());
                assert!(term.threshold_counterterm.metadata_registry.is_some());
                assert_eq!(
                    term.threshold_counterterm.helper_evaluators.len(),
                    round_tripped
                        .derived_data
                        .raised_data
                        .pass_two_evaluator
                        .as_ref()
                        .unwrap()
                        .len(),
                );

                let normalized_graph = Graph::from_string(&normalized, &model).unwrap().remove(0);
                let loop_count = normalized_graph.get_loop_number();
                // A deterministic timelike external momentum is required here: the generic
                // random 1 -> 1 sample is lightlike and leaves this massless bubble below its
                // threshold, so it cannot exercise the detailed L/I event decomposition.
                let mut runtime = RuntimeSettings {
                    kinematics: KinematicsSettings {
                        e_cm: 4.0,
                        externals: Externals::Constant {
                            momenta: vec![
                                ExternalMomenta::Independent([F(4.0), F(0.0), F(0.0), F(0.0)]),
                                ExternalMomenta::Dependent(Dep::Dep),
                            ],
                            helicities: vec![Helicity::ZERO; 2],
                            improvement_settings: PhaseSpaceImprovementSettings::default(),
                            f_64_cache: None,
                            f_128_cache: None,
                        },
                    },
                    ..RuntimeSettings::default()
                };
                runtime.general.generate_events = true;
                runtime.general.store_additional_weights_in_event = true;
                let mut amplitude = super::Amplitude::from_graph_list(
                    "amplitude_threshold_event_roundtrip",
                    vec![normalized_graph],
                )
                .unwrap();
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(1)
                    .build()
                    .unwrap();
                amplitude
                    .preprocess(&model, &generation_settings, &(&runtime).into(), &pool)
                    .unwrap();
                amplitude
                    .build_integrand(
                        &model,
                        "amplitude_threshold_event_roundtrip",
                        &GlobalSettings {
                            generation: generation_settings.clone(),
                            ..Default::default()
                        },
                        (&runtime).into(),
                        &pool,
                    )
                    .unwrap();
                let integrand = amplitude.integrand.as_mut().unwrap();
                integrand.warm_up(&model).unwrap();
                let result = integrand
                    .evaluate_momentum_configuration(
                        &model,
                        &crate::integrands::process::MomentumSpaceEvaluationInput {
                            loop_momenta: (0..loop_count)
                                .map(|index| {
                                    crate::momentum::ThreeMomentum::new(
                                        F(0.35 + index as f64 * 0.1),
                                        F(-0.2),
                                        F(0.45),
                                    )
                                })
                                .collect(),
                            integrator_weight: F(1.0),
                            graph_id: Some(0),
                            group_id: None,
                            orientation: None,
                            channel_id: None,
                        },
                        false,
                    )
                    .unwrap();
                let events = result
                    .event_groups
                    .iter()
                    .flat_map(|group| group.iter())
                    .collect::<Vec<_>>();
                assert_eq!(events.len(), 1);
                let event = events[0];
                let decomposition = event
                    .additional_weights
                    .threshold_counterterms
                    .as_ref()
                    .expect("explicit normalized amplitude must record threshold components");
                assert!(!decomposition.components.is_empty());
                assert_eq!(decomposition.total(), event.weight);
                assert!(decomposition.components.iter().all(|component| {
                    component.multiplier_values.as_slice() == [F(1.0)]
                        && component.bare.is_some()
                        && component.bare == Some(component.weighted)
                        && !component.evaluation_skipped
                }));
                assert!(event.additional_weights.weights.keys().any(|key| matches!(
                    key,
                    crate::observables::AdditionalWeightKey::AmplitudeThresholdCounterterm { .. }
                )));
            })
            .expect("amplitude threshold DOT-export test thread must start")
            .join()
            .expect("amplitude threshold DOT-export test thread must finish successfully");
    }

    #[test]
    fn grouped_amplitude_threshold_directives_remain_member_local_through_runtime_and_save_load() {
        std::thread::Builder::new()
            .name("grouped-amplitude-threshold-ownership-test".to_string())
            .stack_size(32 * 1024 * 1024)
            .spawn(|| {
                test_initialise().unwrap();
                let model = load_generic_model("scalars");
                let mut member_a: Graph = dot!(
                    digraph grouped_threshold_member_a {
                        graph [group_id=0 is_group_master=true]
                        edge [particle=scalar_0]
                        node [num=1]
                        e [style=invis]
                        e -> A:0 [id=3 particle=scalar_2]
                        B:1 -> e [id=2 particle=scalar_2]
                        A -> B [id=1]
                        A -> B [id=0 lmb_id=0]
                    },
                    "scalars"
                )
                .unwrap();
                let mut member_b: Graph = dot!(
                    digraph grouped_threshold_member_b {
                        graph [group_id=0]
                        edge [particle=scalar_0]
                        node [num=1]
                        e [style=invis]
                        e -> A:0 [id=3 particle=scalar_2]
                        B:1 -> e [id=2 particle=scalar_2]
                        A -> B [id=1 lmb_id=0]
                        A -> B [id=0]
                    },
                    "scalars"
                )
                .unwrap();
                for (graph, name, subspace_edge, expression) in [
                    (&mut member_a, "member_a_one_loop", 0, "1"),
                    (&mut member_b, "member_b_one_loop", 1, "2"),
                ] {
                    graph.threshold_counterterms = Autogen::explicit(ThresholdCountertermSpec {
                        schema_version: 1,
                        cuts: vec![ThresholdCountertermCut {
                            edges: Vec::new(),
                            thresholds: vec![ThresholdCountertermThreshold {
                                edges: vec![super::EdgeIndex(0), super::EdgeIndex(1)],
                                counterterms: vec![ThresholdCountertermVariant {
                                    name: Some(name.to_string()),
                                    center_group: None,
                                    subspace: Some(vec![super::EdgeIndex(subspace_edge)]),
                                    parent_lmb: None,
                                    disable: false,
                                    multiplier: Some(ThresholdCountertermMultiplier {
                                        expression: expression.to_string(),
                                        symmetrize: false,
                                        opaque_derivatives: true,
                                    }),
                                }],
                            }],
                        }],
                    });
                }

                let generation_settings = GenerationSettings {
                    threshold_subtraction: ThresholdSubtractionSettings {
                        enable_thresholds: true,
                        assume_positive_external_energies: false,
                        check_esurface_at_generation: false,
                        ..Default::default()
                    },
                    ..Default::default()
                };
                // Keep the two-member bubble above threshold so both member-local registries
                // are referenced by actual detailed L/I event components.
                let mut runtime = RuntimeSettings {
                    kinematics: KinematicsSettings {
                        e_cm: 4.0,
                        externals: Externals::Constant {
                            momenta: vec![
                                ExternalMomenta::Independent([F(4.0), F(0.0), F(0.0), F(0.0)]),
                                ExternalMomenta::Dependent(Dep::Dep),
                            ],
                            helicities: vec![Helicity::ZERO; 2],
                            improvement_settings: PhaseSpaceImprovementSettings::default(),
                            f_64_cache: None,
                            f_128_cache: None,
                        },
                    },
                    ..RuntimeSettings::default()
                };
                runtime.general.generate_events = true;
                runtime.general.store_additional_weights_in_event = true;
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(1)
                    .build()
                    .unwrap();
                let mut amplitude = super::Amplitude::from_graph_list(
                    "grouped_threshold_ownership",
                    vec![member_a, member_b],
                )
                .unwrap();
                assert_eq!(amplitude.graph_group_structure.len(), 1);
                assert_eq!(
                    (&amplitude.graph_group_structure[super::GroupId(0)])
                        .into_iter()
                        .count(),
                    2
                );
                amplitude
                    .preprocess(&model, &generation_settings, &(&runtime).into(), &pool)
                    .unwrap();

                let expected = [
                    ("grouped_threshold_member_a", "member_a_one_loop", 0, "1"),
                    ("grouped_threshold_member_b", "member_b_one_loop", 1, "2"),
                ];
                for (graph, (_, name, subspace_edge, expression)) in
                    amplitude.graphs.iter().zip(expected)
                {
                    let resolved = graph
                        .derived_data
                        .resolved_threshold_counterterms
                        .as_ref()
                        .unwrap();
                    assert!(!resolved.legacy_equivalent);
                    // The bubble has two raised E-surface groups for this energy-edge set. The
                    // compact declaration resolves once per group, with IDs scoped to this graph
                    // member rather than allocated across the owning amplitude group.
                    assert_eq!(resolved.variants.len(), 2);
                    assert_eq!(
                        resolved
                            .variants
                            .keys()
                            .map(|variant_id| variant_id.0)
                            .collect::<Vec<_>>(),
                        [0, 1]
                    );
                    for variant in &resolved.variants {
                        assert_eq!(variant.name, name);
                        assert_eq!(
                            variant.requested_subspace,
                            Some(vec![super::EdgeIndex(subspace_edge)])
                        );
                        assert_eq!(
                            variant
                                .subspace
                                .iter_basis_edges(graph.derived_data.lmbs.as_ref().unwrap())
                                .collect::<Vec<_>>(),
                            vec![super::EdgeIndex(subspace_edge)]
                        );
                        assert_eq!(variant.multiplier.as_ref().unwrap().expression, expression);
                    }
                }

                for include_autogenerated_fields in [false, true] {
                    let mut exported = String::new();
                    amplitude
                        .write_dot_fmt(
                            &mut exported,
                            &DotExportSettings {
                                include_autogenerated_fields,
                                ..DotExportSettings::default()
                            },
                        )
                        .unwrap();
                    let exported_graphs = Graph::from_string(&exported, &model).unwrap();
                    assert_eq!(exported_graphs.len(), 2);
                    for (graph, (graph_name, name, subspace_edge, expression)) in
                        exported_graphs.iter().zip(expected)
                    {
                        assert_eq!(graph.name, graph_name);
                        let spec = &graph.threshold_counterterms;
                        assert!(!spec.autogenerated);
                        assert_eq!(spec.cuts.len(), 1);
                        assert!(spec.cuts[0].edges.is_empty());
                        let variants = &spec.cuts[0].thresholds[0].counterterms;
                        assert_eq!(variants.len(), 1);
                        assert_eq!(variants[0].name.as_deref(), Some(name));
                        assert_eq!(
                            variants[0].subspace,
                            Some(vec![super::EdgeIndex(subspace_edge)])
                        );
                        assert_eq!(
                            variants[0].multiplier.as_ref().unwrap().expression,
                            expression
                        );
                    }
                }

                amplitude
                    .build_integrand(
                        &model,
                        "grouped_threshold_ownership",
                        &GlobalSettings {
                            generation: generation_settings.clone(),
                            ..Default::default()
                        },
                        (&runtime).into(),
                        &pool,
                    )
                    .unwrap();
                let integrand = amplitude.integrand.as_mut().unwrap();
                integrand.warm_up(&model).unwrap();

                let inspect_ownership = |integrand: &ProcessIntegrand| {
                    let ProcessIntegrand::Amplitude(integrand) = integrand else {
                        panic!("grouped amplitude fixture built a non-amplitude integrand")
                    };
                    assert_eq!(integrand.data.graph_group_structure.len(), 1);
                    assert_eq!(
                        (&integrand.data.graph_group_structure[super::GroupId(0)])
                            .into_iter()
                            .count(),
                        2
                    );
                    for (term, (graph_name, name, subspace_edge, expression)) in
                        integrand.data.graph_terms.iter().zip(expected)
                    {
                        assert_eq!(term.graph.name, graph_name);
                        let counterterm = &term.threshold_counterterm;
                        assert!(!counterterm.legacy_equivalent);
                        assert_eq!(counterterm.variant_metadata.len(), 2);
                        for (variant_id, variant) in counterterm.variant_metadata.iter_enumerated()
                        {
                            assert!(variant_id.0 < 2);
                            assert_eq!(variant.name, name);
                            assert_eq!(
                                variant.requested_subspace,
                                Some(vec![super::EdgeIndex(subspace_edge)])
                            );
                        }
                        let multipliers = counterterm.threshold_multipliers.as_ref().unwrap();
                        assert_eq!(multipliers.left_variants().len(), 2);
                        for (variant_id, reference) in
                            multipliers.left_variants().iter().enumerate()
                        {
                            assert_eq!(reference.variant_id.0, variant_id);
                            assert_eq!(reference.evaluator_id.unwrap().0, 0);
                        }
                        let registry = counterterm.metadata_registry.as_ref().unwrap();
                        assert_eq!(registry.graph_name, graph_name);
                        assert_eq!(registry.variants.len(), 2);
                        for (variant_id, variant) in registry.variants.iter().enumerate() {
                            assert_eq!(variant.variant_id, variant_id);
                            assert_eq!(variant.name, name);
                            assert_eq!(variant.resolved_subspace, [subspace_edge]);
                        }
                        assert_eq!(registry.evaluators.len(), 1);
                        assert_eq!(registry.evaluators[0].evaluator_id, 0);
                        assert_eq!(registry.evaluators[0].variant_ids, [0, 1]);
                        assert_eq!(registry.evaluators[0].expression, expression);
                        assert!(registry.components.iter().all(|component| {
                            component.variant_ids.len() == 1
                                && component.variant_ids[0] < 2
                                && component.evaluator_ids == [Some(0)]
                        }));
                    }
                };
                inspect_ownership(integrand);

                let evaluate_members = |integrand: &mut ProcessIntegrand| {
                    (0..2)
                        .map(|graph_id| {
                            let registry = match integrand {
                                ProcessIntegrand::Amplitude(amplitude_integrand) => {
                                    amplitude_integrand.data.graph_terms[graph_id]
                                        .threshold_counterterm
                                        .metadata_registry
                                        .clone()
                                        .unwrap()
                                }
                                ProcessIntegrand::CrossSection(_) => unreachable!(),
                            };
                            let result = integrand
                                .evaluate_momentum_configuration(
                                    &model,
                                    &MomentumSpaceEvaluationInput {
                                        loop_momenta: vec![crate::momentum::ThreeMomentum::new(
                                            F(0.35),
                                            F(-0.2),
                                            F(0.45),
                                        )],
                                        integrator_weight: F(1.0),
                                        graph_id: Some(graph_id),
                                        group_id: None,
                                        orientation: None,
                                        channel_id: None,
                                    },
                                    false,
                                )
                                .unwrap();
                            let events = result
                                .event_groups
                                .iter()
                                .flat_map(|group| group.iter())
                                .collect::<Vec<_>>();
                            assert_eq!(events.len(), 1);
                            let event = events[0];
                            let decomposition = event
                                .additional_weights
                                .threshold_counterterms
                                .as_ref()
                                .expect("each grouped member must record its own decomposition");
                            assert_eq!(decomposition.total(), event.weight);
                            assert!(!decomposition.components.is_empty());
                            let expected_multiplier = F(if graph_id == 0 { 1.0 } else { 2.0 });
                            for component in &decomposition.components {
                                let metadata = &registry.components[component.component_id];
                                assert_eq!(metadata.variant_ids.len(), 1);
                                assert!(metadata.variant_ids[0] < 2);
                                assert_eq!(metadata.evaluator_ids, [Some(0)]);
                                assert_eq!(
                                    component.multiplier_values.as_slice(),
                                    [expected_multiplier]
                                );
                                assert_eq!(component.effective_multiplier, expected_multiplier);
                            }
                            (
                                event.weight,
                                decomposition
                                    .components
                                    .iter()
                                    .map(|component| {
                                        (
                                            component.component_id,
                                            component.multiplier_values.clone(),
                                            component.effective_multiplier,
                                            component.bare,
                                            component.weighted,
                                            component.evaluation_skipped,
                                        )
                                    })
                                    .collect::<Vec<_>>(),
                            )
                        })
                        .collect::<Vec<_>>()
                };
                let before_save = evaluate_members(integrand);

                let unique = SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_nanos();
                let save_root = std::env::temp_dir().join(format!(
                    "gammalooprs-grouped-amplitude-threshold-{}-{unique}",
                    std::process::id(),
                ));
                amplitude.save(&save_root, true).unwrap();
                let mut state_bytes = Vec::new();
                State::export(&mut state_bytes).unwrap();
                let state_map = State::import(&mut Cursor::new(state_bytes), None).unwrap();
                let context = GammaLoopContextContainer {
                    model: &model,
                    state_map: &state_map,
                };
                let mut loaded =
                    super::Amplitude::load(save_root.join("grouped_threshold_ownership"), context)
                        .unwrap();
                let loaded_integrand = loaded.integrand.as_mut().unwrap();
                loaded_integrand.warm_up(&model).unwrap();
                inspect_ownership(loaded_integrand);
                let after_load = evaluate_members(loaded_integrand);
                assert_eq!(after_load, before_save);
                fs::remove_dir_all(&save_root).unwrap();
            })
            .expect("grouped-amplitude threshold ownership test thread must start")
            .join()
            .expect("grouped-amplitude threshold ownership test must finish successfully");
    }

    #[test]
    fn generalized_amplitude_directive_preserves_raised_threshold_order_and_components() {
        std::thread::Builder::new()
            .name("generalized-raised-amplitude-threshold-test".to_string())
            .stack_size(32 * 1024 * 1024)
            .spawn(|| {
                test_initialise().unwrap();
                let model = crate::model::Model::from_file(
                    std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                        .join("../../assets/models/json/scalars/scalars_2p_3p.json"),
                )
                .unwrap();
                let mut graph = Graph::from_string(
                    include_str!("../../../../tests/resources/graphs/dotted_bubble_amp.dot"),
                    &model,
                )
                .unwrap()
                .remove(0);
                let mut generation_settings = GenerationSettings {
                    threshold_subtraction: ThresholdSubtractionSettings {
                        enable_thresholds: true,
                        assume_positive_external_energies: false,
                        check_esurface_at_generation: false,
                        ..Default::default()
                    },
                    ..Default::default()
                };
                generation_settings
                    .tropical_subgraph_table
                    .disable_tropical_generation = true;

                let mut legacy_graph = AmplitudeGraph::new(graph.clone());
                legacy_graph
                    .preprocess(
                        &model,
                        &generation_settings,
                        &(&RuntimeSettings::default()).into(),
                    )
                    .unwrap();
                let expected_thresholds = [
                    vec![super::EdgeIndex(1), super::EdgeIndex(3)],
                    vec![super::EdgeIndex(2), super::EdgeIndex(3)],
                ];
                let legacy_resolved = legacy_graph
                    .derived_data
                    .resolved_threshold_counterterms
                    .as_ref()
                    .unwrap();
                assert_eq!(legacy_resolved.variants.len(), 2);
                for variant in &legacy_resolved.variants {
                    assert_eq!(variant.raised_esurface_group.max_occurence, 2);
                    assert_eq!(
                        variant
                            .associations
                            .iter()
                            .map(|association| association.threshold_edges.clone())
                            .collect::<Vec<_>>(),
                        expected_thresholds
                    );
                    assert_eq!(variant.subspace_loop_count, 1);
                    assert_eq!(
                        variant
                            .subspace
                            .iter_basis_edges(legacy_graph.derived_data.lmbs.as_ref().unwrap())
                            .collect::<Vec<_>>(),
                        [super::EdgeIndex(3)]
                    );
                }

                let explicit_variant = ThresholdCountertermVariant {
                    name: Some("raised_one_loop".to_string()),
                    center_group: None,
                    subspace: Some(vec![super::EdgeIndex(3)]),
                    parent_lmb: None,
                    disable: false,
                    multiplier: Some(ThresholdCountertermMultiplier {
                        expression: "2".to_string(),
                        symmetrize: false,
                        opaque_derivatives: true,
                    }),
                };
                graph.threshold_counterterms = Autogen::explicit(ThresholdCountertermSpec {
                    schema_version: 1,
                    cuts: vec![ThresholdCountertermCut {
                        edges: Vec::new(),
                        thresholds: expected_thresholds
                            .iter()
                            .map(|edges| ThresholdCountertermThreshold {
                                edges: edges.clone(),
                                counterterms: vec![explicit_variant.clone()],
                            })
                            .collect(),
                    }],
                });
                let mut runtime = RuntimeSettings {
                    kinematics: KinematicsSettings {
                        e_cm: 4.0,
                        externals: Externals::Constant {
                            momenta: vec![
                                ExternalMomenta::Independent([
                                    F(4.0),
                                    F(0.0),
                                    F(0.0),
                                    F(0.0),
                                ]),
                                ExternalMomenta::Dependent(Dep::Dep),
                            ],
                            helicities: vec![Helicity::ZERO; 2],
                            improvement_settings: PhaseSpaceImprovementSettings::default(),
                            f_64_cache: None,
                            f_128_cache: None,
                        },
                    },
                    ..RuntimeSettings::default()
                };
                runtime.general.generate_events = true;
                runtime.general.store_additional_weights_in_event = true;
                let pool = rayon::ThreadPoolBuilder::new()
                    .num_threads(1)
                    .build()
                    .unwrap();
                let mut amplitude = super::Amplitude::from_graph_list(
                    "generalized_raised_amplitude_threshold",
                    vec![graph],
                )
                .unwrap();
                amplitude
                    .preprocess(&model, &generation_settings, &(&runtime).into(), &pool)
                    .unwrap();
                let resolved = amplitude.graphs[0]
                    .derived_data
                    .resolved_threshold_counterterms
                    .as_ref()
                    .unwrap();
                assert!(!resolved.legacy_equivalent);
                assert_eq!(resolved.variants.len(), 2);
                for variant in &resolved.variants {
                    assert_eq!(variant.name, "raised_one_loop");
                    assert_eq!(variant.raised_esurface_group.max_occurence, 2);
                    assert_eq!(variant.subspace_loop_count, 1);
                    assert_eq!(variant.multiplier.as_ref().unwrap().expression, "2");
                }
                for symbolic in &amplitude.graphs[0]
                    .derived_data
                    .threshold_counterterm_variants
                {
                    let orders = symbolic
                        .atom
                        .parametric
                        .iter()
                        .filter_map(|(index, _)| index.left_threshold_order)
                        .collect::<std::collections::BTreeSet<_>>();
                    assert_eq!(orders, [1, 2].into_iter().collect());
                }

                amplitude
                    .build_integrand(
                        &model,
                        "generalized_raised_amplitude_threshold",
                        &GlobalSettings {
                            generation: generation_settings.clone(),
                            ..Default::default()
                        },
                        (&runtime).into(),
                        &pool,
                    )
                    .unwrap();
                let inspect_raised_runtime = |integrand: &ProcessIntegrand| {
                    let ProcessIntegrand::Amplitude(integrand) = integrand else {
                        panic!("raised amplitude fixture built a non-amplitude integrand")
                    };
                    let counterterm = &integrand.data.graph_terms[0].threshold_counterterm;
                    assert!(!counterterm.legacy_equivalent);
                    assert_eq!(counterterm.variant_evaluators.len(), 2);
                    for (variant_id, evaluator) in
                        counterterm.variant_evaluators.iter_enumerated()
                    {
                        assert_eq!(
                            evaluator
                                .evaluator_stacks
                                .keys()
                                .filter_map(|index| index.left_threshold_order)
                                .collect::<std::collections::BTreeSet<_>>(),
                            [1, 2].into_iter().collect()
                        );
                        let helpers = &counterterm.variant_helper_evaluators[variant_id];
                        assert_eq!(helpers.len(), 2);
                        assert!(helpers.iter().all(|helper| helper.exprs_len == 2));
                    }

                    let registry = counterterm.metadata_registry.as_ref().unwrap();
                    assert_eq!(registry.variants.len(), 2);
                    assert!(registry.variants.iter().all(|variant| {
                        variant.name == "raised_one_loop"
                            && variant.subspace_loop_count == 1
                            && variant.resolved_subspace == [3]
                    }));
                    assert_eq!(registry.evaluators.len(), 1);
                    assert_eq!(registry.evaluators[0].variant_ids, [0, 1]);
                    assert_eq!(registry.evaluators[0].expression, "2");
                    assert_eq!(registry.components.len(), 4);
                    for variant_id in 0..2 {
                        assert_eq!(
                            registry
                                .components
                                .iter()
                                .filter(|component| component.variant_ids == [variant_id])
                                .map(|component| component.kind)
                                .collect::<std::collections::BTreeSet<_>>(),
                            [
                                ThresholdCountertermComponentKind::Local,
                                ThresholdCountertermComponentKind::Integrated,
                            ]
                            .into_iter()
                            .collect()
                        );
                    }
                    registry.clone()
                };
                let registry_before_save =
                    inspect_raised_runtime(amplitude.integrand.as_ref().unwrap());

                let evaluate_raised = |integrand: &mut ProcessIntegrand| {
                    let result = integrand
                        .evaluate_momentum_configuration(
                            &model,
                            &MomentumSpaceEvaluationInput {
                                loop_momenta: vec![crate::momentum::ThreeMomentum::new(
                                    F(0.35),
                                    F(-0.2),
                                    F(0.45),
                                )],
                                integrator_weight: F(1.0),
                                graph_id: Some(0),
                                group_id: None,
                                orientation: None,
                                channel_id: None,
                            },
                            false,
                        )
                        .unwrap();
                    let events = result
                        .event_groups
                        .iter()
                        .flat_map(|group| group.iter())
                        .collect::<Vec<_>>();
                    assert_eq!(events.len(), 1);
                    let event = events[0];
                    let decomposition = event
                        .additional_weights
                        .threshold_counterterms
                        .as_ref()
                        .expect("raised generalized threshold must record completed L/I pieces");
                    assert_eq!(decomposition.total(), event.weight);
                    assert_eq!(decomposition.components.len(), 2);
                    assert_eq!(
                        decomposition
                            .components
                            .iter()
                            .map(|component| {
                                registry_before_save.components[component.component_id].kind
                            })
                            .collect::<std::collections::BTreeSet<_>>(),
                        [
                            ThresholdCountertermComponentKind::Local,
                            ThresholdCountertermComponentKind::Integrated,
                        ]
                        .into_iter()
                        .collect()
                    );
                    for component in &decomposition.components {
                        let metadata = &registry_before_save.components[component.component_id];
                        assert_eq!(metadata.variant_ids.len(), 1);
                        assert!(metadata.variant_ids[0] < 2);
                        assert_eq!(metadata.evaluator_ids, [Some(0)]);
                        assert_eq!(component.multiplier_values.as_slice(), [F(2.0)]);
                        assert_eq!(component.effective_multiplier, F(2.0));
                        assert!(component.bare.is_some());
                        assert!(!component.evaluation_skipped);
                        let crate::observables::ThresholdCountertermComponentOccurrence::Amplitude {
                            raised_esurface_id,
                            overlap_group,
                        } = &component.occurrence
                        else {
                            panic!("amplitude fixture recorded a non-amplitude occurrence")
                        };
                        assert!(event.additional_weights.weights.contains_key(
                            &crate::observables::AdditionalWeightKey::AmplitudeThresholdCountertermVariant {
                                variant_id: metadata.variant_ids[0],
                                esurface_id: *raised_esurface_id,
                                overlap_group: *overlap_group,
                            },
                        ));
                    }
                    (
                        event.weight,
                        decomposition.original,
                        decomposition
                            .components
                            .iter()
                            .map(|component| {
                                (
                                    component.component_id,
                                    component.occurrence.clone(),
                                    component.multiplier_values.clone(),
                                    component.effective_multiplier,
                                    component.bare,
                                    component.weighted,
                                    component.evaluation_skipped,
                                )
                            })
                            .collect::<Vec<_>>(),
                    )
                };
                let integrand = amplitude.integrand.as_mut().unwrap();
                integrand.warm_up(&model).unwrap();
                let before_save = evaluate_raised(integrand);

                let unique = SystemTime::now()
                    .duration_since(UNIX_EPOCH)
                    .unwrap()
                    .as_nanos();
                let save_root = std::env::temp_dir().join(format!(
                    "gammalooprs-raised-amplitude-threshold-{}-{unique}",
                    std::process::id(),
                ));
                amplitude.save(&save_root, true).unwrap();
                let mut state_bytes = Vec::new();
                State::export(&mut state_bytes).unwrap();
                let state_map = State::import(&mut Cursor::new(state_bytes), None).unwrap();
                let context = GammaLoopContextContainer {
                    model: &model,
                    state_map: &state_map,
                };
                let mut loaded = super::Amplitude::load(
                    save_root.join("generalized_raised_amplitude_threshold"),
                    context,
                )
                .unwrap();
                let loaded_integrand = loaded.integrand.as_mut().unwrap();
                assert_eq!(
                    inspect_raised_runtime(loaded_integrand),
                    registry_before_save
                );
                loaded_integrand.warm_up(&model).unwrap();
                assert_eq!(evaluate_raised(loaded_integrand), before_save);
                fs::remove_dir_all(&save_root).unwrap();
            })
            .unwrap()
            .join()
            .unwrap();
    }

    #[test]
    fn amplitude_tree() {
        test_initialise().unwrap();
        let mut graph: AmplitudeGraph = dot!(digraph qqx_aaa_tree_1 {
                num="spenso::g(spenso::dind(spenso::cof(3, hedge(1))), spenso::cof(3, hedge(2)))/3"
                ext    [style=invis]
                ext -> v1:1 [particle="d" id=1];
                ext -> v3:2 [particle="d~" id=2];
                v1:3 -> ext [particle="a" id=3];
                v2:4 -> ext [particle="a" id=4];
                v3:0 -> ext [particle="a" id=0];
                v1 -> v2 [particle="d" id=5];
                v2 -> v3 [particle="d" id=6];
    })
    .unwrap();

        let _model = load_generic_model("sm");

        graph.generate_cff(&OrientationPattern::default()).unwrap();
        // graph.build_parametric_integrand(&GenerationSettings::default());

        let param_builder = &graph.graph.param_builder;
        println!("{param_builder}");

        // GenericEvaluator::new_from_builder(
        //     [GS.orientation_delta(&EdgeVec::from_iter(vec![Orientation::Default; 7]))],
        //     param_builder,
        //     None,
        //     OptimizationSettings::default(),
        //     true,
        // )
        // .unwrap();
        // println!(" {}", a);
    }

    #[test]
    fn generation_orientation_pattern_filters_evaluator_orientations() {
        test_initialise().unwrap();
        let mut graph: AmplitudeGraph = dot!(
            digraph bub {
                edge [particle=scalar_1]
                node [num=1]
                e [style=invis]
                e -> A:0 [id=3]
                B:1 -> e [id=2]
                A -> B [id=1]
                A -> B [id=0]
            },
            "scalars"
        )
        .unwrap();

        let model = load_generic_model("scalars");
        let generation_settings = GenerationSettings {
            threshold_subtraction: ThresholdSubtractionSettings {
                enable_thresholds: false,
                ..Default::default()
            },
            ..Default::default()
        };
        let runtime_settings = RuntimeSettings::default();
        graph
            .preprocess(&model, &generation_settings, &(&runtime_settings).into())
            .unwrap();

        assert!(
            graph
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .orientations
                .len()
                > 1
        );

        let global_settings = GlobalSettings {
            generation: GenerationSettings {
                orientation_pattern: OrientationPattern::from_orientation(
                    &graph
                        .derived_data
                        .cff_expression
                        .as_ref()
                        .unwrap()
                        .orientations[OrientationID(0)],
                ),
                threshold_subtraction: ThresholdSubtractionSettings {
                    enable_thresholds: false,
                    ..Default::default()
                },
                ..Default::default()
            },
            ..Default::default()
        };

        let (term, _stats) = AmplitudeGraphTerm::from_amplitude_graph(
            &graph,
            GraphGroupPosition(0),
            TiVec::new(),
            &model,
            &global_settings,
        )
        .unwrap();

        assert_eq!(term.orientations.len(), 1);
        assert_eq!(
            term.orientations[OrientationID(0)],
            graph
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .orientations[OrientationID(0)]
            .data
            .orientation
        );
    }
}
