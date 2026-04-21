use std::{
    fmt::{Display, Formatter},
    fs::{self, File},
    io::Write,
    ops::{Index, IndexMut},
    path::{Path, PathBuf},
};

use ahash::HashMap;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use itertools::Itertools;
use rayon::{
    ThreadPool,
    iter::{IntoParallelRefMutIterator, ParallelIterator},
};
use spenso::algebra::algebraic_traits::IsZero;
use tracing::info;
use vakint::Vakint;

use crate::{
    GammaLoopContext, GammaLoopContextContainer,
    cff::{
        esurface::{RaisedEsurfaceData, RaisedEsurfaceGroup},
        expression::{CFFExpression, OrientationID},
    },
    debug_tags, define_index,
    graph::{
        GraphGroup, GroupId, LMBext, LmbIndex, LoopMomentumBasis,
        cuts::{CutSet, ResidueSelector},
        parse::complete_group_parsing,
    },
    integrands::process::{
        GenericEvaluator, LmbMultiChannelingSetup, ParamBuilder,
        cross_section::CrossSectionIntegrandData, graph_to_group_id_for_group_structure,
    },
    model::ArcParticle,
    momentum::{
        Helicity,
        sample::{ExternalIndex, SubspaceData},
    },
    processes::{
        DotExportSettings, EvaluatorSettings, GraphGenerationStats, NamedGraphGenerationReport,
    },
    settings::{GlobalSettings, global::GenerationSettings, runtime::LockedRuntimeSettings},
    utils::{GS, hyperdual_utils::shape_for_t_derivatives},
    uv::{approx::CutStructure, forest::ParametricIntegrands, wood::CutWoods},
};
use eyre::{Context, eyre};
use linnet::half_edge::{
    involution::{EdgeVec, Orientation},
    subgraph::{
        HedgeNode, Inclusion, InternalSubGraph, OrientedCut, SuBitGraph, SubGraphLike, SubSetOps,
    },
};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    evaluate::FunctionMap,
    function, parse, symbol,
};
use tracing::{debug, warn};
use typed_index_collections::{TiVec, ti_vec};

use crate::{
    cff::esurface::{Esurface, EsurfaceID},
    graph::{ExternalConnection, FeynmanGraph, Graph},
    integrands::process::{
        ProcessIntegrand,
        cross_section::{CrossSectionGraphTerm, CrossSectionIntegrand},
    },
    model::Model,
};

use crate::processes::ProcessDefinition;

#[derive(Clone, Debug, Encode, Decode)]
pub struct IteratedCtCollection<T> {
    data: Vec<T>,
    num_left_thresholds: usize,
}

impl<T> IteratedCtCollection<T> {
    pub fn map_ref<U, F>(&self, f: F) -> IteratedCtCollection<U>
    where
        F: Fn(&T) -> U,
    {
        let data = self.data.iter().map(f).collect();
        IteratedCtCollection {
            data,
            num_left_thresholds: self.num_left_thresholds,
        }
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }

    pub(crate) fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.data.iter_mut()
    }

    pub(crate) fn num_left_thresholds(&self) -> usize {
        self.num_left_thresholds
    }
}

impl<T> Index<(LeftThresholdId, RightThresholdId)> for IteratedCtCollection<T> {
    type Output = T;

    fn index(&self, index: (LeftThresholdId, RightThresholdId)) -> &Self::Output {
        let (left_id, right_id) = index;
        &self.data[left_id.0 * self.num_left_thresholds + right_id.0]
    }
}

impl<T> IndexMut<(LeftThresholdId, RightThresholdId)> for IteratedCtCollection<T> {
    fn index_mut(&mut self, index: (LeftThresholdId, RightThresholdId)) -> &mut Self::Output {
        let (left_id, right_id) = index;
        &mut self.data[left_id.0 * self.num_left_thresholds + right_id.0]
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct LUCounterTermData {
    pub left_thresholds: TiVec<LeftThresholdId, EsurfaceID>,
    pub right_thresholds: TiVec<RightThresholdId, EsurfaceID>,
    pub left_atoms: TiVec<LeftThresholdId, ParametricIntegrands>,
    pub right_atoms: TiVec<RightThresholdId, ParametricIntegrands>,
    pub iterated: IteratedCtCollection<ParametricIntegrands>,
}

define_index! {pub struct GlobalThresholdId;}
define_index! {pub struct RightThresholdId;}
define_index! {pub struct LeftThresholdId;}
define_index! {pub struct RaisedCutId;}

use derive_more::{From, Into};
#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSection {
    pub name: String,
    pub integrand: Option<ProcessIntegrand>,
    pub supergraphs: Vec<CrossSectionGraph>,
    pub external_particles: Vec<ArcParticle>,
    pub external_connections: Vec<ExternalConnection>,
    pub n_incmoming: usize,
    pub graph_group_structure: TiVec<GroupId, GraphGroup>,
}

impl CrossSection {
    fn storage_path(&self, base: &Path) -> PathBuf {
        base.join(&self.name)
    }

    pub fn export_standalone(
        &self,
        path: impl AsRef<Path>,
        settings: &crate::processes::StandaloneExportSettings,
    ) -> Result<()> {
        if let Some(integrand) = &self.integrand {
            integrand.export_standalone(path, settings)?;
        } else {
            return Err(eyre!(
                "Cannot export standalone cross section {} without integrand",
                self.name
            ));
        }

        Ok(())
    }

    #[allow(dead_code)]
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::io::Error> {
        for graph in &self.supergraphs {
            graph.write_dot(writer, settings)?;
            writeln!(writer)?;
        }
        Ok(())
    }

    pub fn write_dot_fmt<W: std::fmt::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::fmt::Error> {
        for graph in &self.supergraphs {
            graph.write_dot_fmt(writer, settings)?;
            writeln!(writer)?;
        }
        Ok(())
    }

    pub(crate) fn new(name: String) -> Self {
        Self {
            name,
            integrand: None,
            supergraphs: vec![],
            external_connections: vec![],
            external_particles: vec![],
            n_incmoming: 0,
            graph_group_structure: TiVec::new(),
        }
    }

    pub fn from_graph_list(name: String, mut graphs: Vec<Graph>, _model: &Model) -> Result<Self> {
        let mut cross_section = CrossSection::new(name);
        cross_section.graph_group_structure = complete_group_parsing(&mut graphs)?;
        // println!("group structure: {:?}", cross_section.graph_group_structure);

        for cross_section_graph in graphs {
            // cross_section_graph.param_builder =  ParamBuilder::new(&cross_section_graph, model);
            cross_section.add_supergraph(cross_section_graph)?;
        }
        Ok(cross_section)
    }

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

    fn add_supergraph(&mut self, supergraph: Graph) -> Result<()> {
        if self.external_particles.is_empty() {
            let external_particles = supergraph.get_external_partcles();
            if !external_particles.len().is_multiple_of(2) {
                return Err(eyre!(
                    "expected even number of externals for forward scattering graph"
                ));
            }
            self.external_particles = external_particles;
            self.n_incmoming = self.external_particles.len() / 2;
        } else if self.external_particles != supergraph.get_external_partcles() {
            return Err(eyre!(
                "attempt to add supergraph with differnt external particles"
            ));
        }

        let cross_section_graph = CrossSectionGraph::new(supergraph);
        self.supergraphs.push(cross_section_graph);

        // TODO: validate that the graph is compatible
        Ok(())
    }

    pub fn preprocess(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
        global_settings: &GenerationSettings,
        runtime_default: LockedRuntimeSettings,
        generation_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        let integrand_name = self.name.clone();
        generation_pool.install(|| {
            self.supergraphs
                .par_iter_mut()
                .map(|supergraph| {
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    let stats = supergraph.preprocess(
                        model,
                        process_definition,
                        global_settings,
                        runtime_default,
                    )?;
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    Ok(NamedGraphGenerationReport {
                        integrand_name: integrand_name.clone(),
                        graph_name: supergraph.graph.name.clone(),
                        stats,
                    })
                })
                .collect::<Result<Vec<_>>>()
        })
    }

    pub fn build_integrand(
        &mut self,
        model: &Model,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        generation_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        if crate::is_interrupted() {
            return Err(eyre!("Generation interrupted by user"));
        }
        let integrand_name = self.name.clone();
        let mut graph_reports = Vec::new();
        let terms = generation_pool.install(|| {
            self.supergraphs
                .par_iter_mut()
                .map(|sg| {
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    let graph_started = std::time::Instant::now();
                    let (term, mut stats) = sg.generate_term_for_graph(model, global_settings)?;
                    if crate::is_interrupted() {
                        return Err(eyre!("Generation interrupted by user"));
                    }
                    stats.total_time += graph_started.elapsed();
                    Ok((
                        term,
                        NamedGraphGenerationReport {
                            integrand_name: integrand_name.clone(),
                            graph_name: sg.graph.name.clone(),
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
            let master = group.master();
            let mc_of_master = self.supergraphs[master]
                .derived_data
                .multi_channeling_setup
                .as_ref()
                .unwrap();

            for graph_id in group.into_iter() {
                terms[graph_id].multi_channeling_setup = mc_of_master.clone();
            }
        }

        let mut cross_section_integrand = CrossSectionIntegrand {
            settings: runtime_default.into(),
            data: CrossSectionIntegrandData {
                compilation: global_settings
                    .generation
                    .compile
                    .frozen_mode(&global_settings.generation.evaluator),
                loop_cache_id: 0,
                external_cache_id: 0,
                base_external_cache_id: 0,
                rotations: None,
                name: self.name.clone(),
                external_connections: self.external_connections.clone(),
                n_incoming: self.n_incmoming,
                // polarizations,
                graph_terms: terms,
                graph_group_structure: self.graph_group_structure.clone(),
                graph_to_group_id: graph_to_group_id_for_group_structure(
                    &self.graph_group_structure,
                ),
            },
            event_processing_runtime: Default::default(),
            active_f64_backend: Default::default(),
        };
        let compile_times = cross_section_integrand
            .prepare_runtime_backends_after_generation_with_compile_times()?;
        for (report, compile_time) in graph_reports.iter_mut().zip(compile_times) {
            report.stats.evaluator_compile_time += compile_time;
            report.stats.total_time += compile_time;
        }

        self.integrand = Some(ProcessIntegrand::CrossSection(cross_section_integrand));
        Ok(graph_reports)
    }

    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        thread_pool: &ThreadPool,
    ) -> Result<Vec<NamedGraphGenerationReport>> {
        info!("Compiling cross section {}", self.name);
        let p = self.storage_path(path.as_ref());

        let result = fs::create_dir_all(&p).with_context(|| {
            format!(
                "Trying to create directory to compile cross section {}",
                p.display()
            )
        });

        if override_existing {
            result?;
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
        }
        Ok(Vec::new())
    }

    pub fn save(&mut self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let p = self.storage_path(path.as_ref());
        let r = fs::create_dir_all(&p).with_context(|| {
            format!(
                "Trying to create directory to save cross section {}",
                p.display()
            )
        });

        if override_existing {
            r?;
        }

        let integrand = self.integrand.take();
        if let Some(integrand) = &integrand {
            integrand.save(&p, override_existing)?;
        }

        let binary = bincode::encode_to_vec(&(*self), bincode::config::standard())?;
        if override_existing {
            fs::write(p.join("cs.bin"), &binary)?;
        } else {
            let mut file = File::create_new(p.join("cs.bin"))?;
            file.write_all(&binary)?;
        }

        self.integrand = integrand;
        Ok(())
    }

    pub(crate) fn load(path: impl AsRef<Path>, context: GammaLoopContextContainer) -> Result<Self> {
        let binary = fs::read(path.as_ref().join("cs.bin"))?;
        let (mut cs, _): (Self, _) =
            bincode::decode_from_slice_with_context(&binary, bincode::config::standard(), context)?;

        if path.as_ref().join("integrand").exists() {
            let integrand = CrossSectionIntegrand::load(path.as_ref().join("integrand"), context)?;
            cs.integrand = Some(ProcessIntegrand::CrossSection(integrand));
        }

        Ok(cs)
    }
}

#[derive(Clone, bincode::Encode, bincode::Decode)]
pub struct CrossSectionCut {
    pub cut: OrientedCut,
    pub left: SuBitGraph,
    pub right: SuBitGraph,
}

impl CrossSectionCut {
    pub(crate) fn is_s_channel(&self, cross_section_graph: &CrossSectionGraph) -> Result<bool> {
        let nodes_of_left_cut: Vec<_> = cross_section_graph
            .graph
            .underlying
            .iter_nodes_of(&self.left)
            .map(|(nid, _, _)| nid)
            .collect();

        let left_node = cross_section_graph
            .graph
            .underlying
            .combine_to_single_hedgenode(&nodes_of_left_cut);
        let res = left_node.includes(&cross_section_graph.source_nodes)
            && cross_section_graph.target_nodes.weakly_disjoint(&left_node);

        if !res {
            warn!("s channel check wrong");
        }

        Ok(true)
    }

    pub(crate) fn is_valid_for_process(
        &self,
        cross_section_graph: &CrossSectionGraph,
        process: &ProcessDefinition,
        model: &Model,
    ) -> Result<bool> {
        if self.is_s_channel(cross_section_graph)? {
            let cut_content_builder = self
                .cut
                .iter_edges(&cross_section_graph.graph.underlying)
                .filter_map(|(orientation, edge_data)| {
                    Some(if orientation == Orientation::Reversed {
                        edge_data.data.particle()?.get_anti_particle(model)
                    } else {
                        edge_data.data.particle()?.clone()
                    })
                })
                .collect_vec();

            let any_pdg_list_passes = process
                .final_pdgs_lists
                .iter()
                .map(|x| {
                    x.iter()
                        .map(|pdg| model.get_particle_from_pdg(*pdg as isize))
                })
                .any(|particle_content| {
                    let mut cut_content = cut_content_builder.clone();
                    debug!(
                        "cut content: {:?}",
                        cut_content.iter().map(|p| p.name.clone()).collect_vec()
                    );

                    for particle in particle_content {
                        if let Some(index) = cut_content.iter().position(|p| p == &particle) {
                            cut_content.remove(index);
                        } else {
                            debug!("wrong particles");
                            return false;
                        }
                    }

                    let (n_unresolved, unresolved_cut_content) =
                        process.unresolved_cut_content(model);

                    if cut_content.len() > n_unresolved {
                        debug!(" too many unresolved particles");
                        return false;
                    }

                    if !cut_content
                        .iter()
                        .all(|particle| unresolved_cut_content.contains(particle))
                    {
                        debug!("wrong unresolved particles");
                        return false;
                    }

                    true
                });

            if !any_pdg_list_passes {
                debug!("wrong pdg list");
                return Ok(false);
            }

            let amplitude_couplings = process.amplitude_filters.get_coupling_orders();
            let amplitude_loop_count = process.amplitude_filters.get_loop_count_range();

            if let Some((min_loop, max_loop)) = amplitude_loop_count {
                let loop_range = min_loop..=max_loop;
                let left_internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(
                    self.left.clone(),
                    &cross_section_graph.graph.underlying,
                );

                let right_internal_subgraph = InternalSubGraph::cleaned_filter_pessimist(
                    self.right.clone(),
                    &cross_section_graph.graph.underlying,
                );

                let left_loop = cross_section_graph
                    .graph
                    .underlying
                    .cyclotomatic_number(&left_internal_subgraph);

                let right_loop = cross_section_graph
                    .graph
                    .underlying
                    .cyclotomatic_number(&right_internal_subgraph);

                let total_loops = left_loop + right_loop;

                if !loop_range.contains(&total_loops) {
                    debug!("incorrect loop count");
                    return Ok(false);
                }
            }

            if amplitude_couplings.is_some() {
                todo!("waiting for update")
            }

            Ok(true)
        } else {
            debug!("cut is not s channel");
            Ok(false)
        }
    }
}

#[derive(
    Debug, Clone, Serialize, Decode, Deserialize, From, Into, Hash, PartialEq, Copy, Eq, Encode,
)]
pub struct CutId(pub usize);

impl Display for CutId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionGraph {
    pub graph: Graph,
    pub source_nodes: HedgeNode,
    pub target_nodes: HedgeNode,
    pub cuts: TiVec<CutId, CrossSectionCut>,
    pub cut_esurface: TiVec<CutId, Esurface>,
    pub cut_esurface_id_map: TiVec<CutId, EsurfaceID>,
    pub derived_data: CrossSectionDerivedData,
}

impl CrossSectionGraph {
    pub(crate) fn new(graph: Graph) -> Self {
        let (source_node, target_node) = graph.get_source_and_target();

        Self {
            graph,
            source_nodes: source_node,
            target_nodes: target_node,
            cuts: TiVec::new(),
            cut_esurface: TiVec::new(),
            cut_esurface_id_map: TiVec::new(),
            derived_data: CrossSectionDerivedData::new_empty(),
        }
    }

    pub(crate) fn apply_spin_sum(
        &mut self,
        model: &Model,
        generation_settings: &GenerationSettings,
        locked_runtime_settings: &LockedRuntimeSettings,
    ) -> Result<()> {
        for (extid, hel) in locked_runtime_settings.helicities().iter().enumerate() {
            //println!("{extid},{hel:?}");
            let eid = self.graph.loop_momentum_basis.ext_edges[ExternalIndex(extid)];

            let Some(p) = self.graph.underlying[eid].particle() else {
                continue;
            };

            match hel {
                Helicity::Summed => {
                    let Some(p) = p.polarization_sum(
                        eid,
                        false,
                        generation_settings.vector_polarization_sum_gauge,
                    )?
                    else {
                        continue;
                    };
                    self.graph.global_prefactor.projector =
                        self.graph.global_prefactor.projector.replace_multiple(&[p]);
                }
                Helicity::SummedAveraged => {
                    let Some(p) = p.polarization_sum(
                        eid,
                        true,
                        generation_settings.vector_polarization_sum_gauge,
                    )?
                    else {
                        continue;
                    };
                    self.graph.global_prefactor.projector =
                        self.graph.global_prefactor.projector.replace_multiple(&[p]);
                }
                _ => {}
            }
        }

        self.graph.polarizations = self.graph.global_prefactor.polarizations();
        self.graph.param_builder = ParamBuilder::new(
            &self.graph,
            model,
            &self.graph.loop_momentum_basis,
            &self.graph.param_builder.pairs.additional_params.params,
        );
        Ok(())
    }

    pub(crate) fn preprocess(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
        settings: &GenerationSettings,
        runtime_default: LockedRuntimeSettings,
    ) -> Result<GraphGenerationStats> {
        let preprocess_started = std::time::Instant::now();
        let mut stats = GraphGenerationStats::default();
        self.apply_spin_sum(model, settings, &runtime_default)?;
        debug_tags!(#generation; "generating cuts");
        self.generate_cuts(model, process_definition, settings)?;
        debug_tags!(#generation; "generating esurfaces corresponding to cuts");
        self.generate_esurface_cuts();
        debug_tags!(#generation; "generating cff");
        stats.merge_in_place(&self.generate_cff(settings)?);
        debug_tags!(#generation; "building lmbs");
        self.build_lmbs()?;
        debug_tags!(#generation; "building multi channeling channels");

        if self.graph.is_group_master {
            self.build_multi_channeling_channels(settings.override_lmb_heuristics);
        }

        let vk = crate::utils::vakint()?;
        debug_tags!(#generation; "building parametric integrand");
        self.build_parametric_integrand(settings, vk)?;
        //self.build_parametric_integrand_raised_cuts(settings)?;

        if settings.threshold_subtraction.enable_thresholds {
            debug_tags!(#generation, #subtraction; "building threshold counterterm");
            self.build_subspace_data()?;
            self.build_threshold_counterterm(settings, vk)?;
        }

        stats.total_time += preprocess_started.elapsed();
        Ok(stats)
    }

    #[allow(dead_code)]
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::io::Error> {
        self.graph.dot_serialize_io(writer, settings)
    }

    pub(crate) fn write_dot_fmt<W: std::fmt::Write>(
        &self,
        writer: &mut W,
        settings: &DotExportSettings,
    ) -> Result<(), std::fmt::Error> {
        self.graph.dot_serialize_fmt(writer, settings)
    }

    fn generate_cff(&mut self, settings: &GenerationSettings) -> Result<GraphGenerationStats> {
        let canonize_esurface = self
            .graph
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let contract_edges = self
            .graph
            .iter_edges_of(
                &self
                    .graph
                    .tree_edges
                    .subtract(&self.graph.initial_state_cut),
            )
            .map(|x| x.1)
            .collect_vec();

        let global_cff = self
            .graph
            .generate_cff(&contract_edges, &canonize_esurface)?;

        let cut_esurface_map = self
            .cut_esurface
            .iter()
            .map(|esurface| {
                self.graph
                    .surface_cache
                    .esurface_cache
                    .iter()
                    .position(|e_sf| e_sf == esurface)
                    .unwrap_or_else(|| {
                        println!("esurfaces corruped");

                        println!("cut esurfaces: {:?}", self.cut_esurface);
                        println!(
                            "graph esurfaces: {:?}",
                            self.graph.surface_cache.esurface_cache
                        );
                        println!(
                            "graph hsurfaces: {:?}",
                            self.graph.surface_cache.hsurface_cache
                        );
                        panic!()
                    })
                    .into()
            })
            .collect();

        self.cut_esurface_id_map = cut_esurface_map;

        let esurface_raised_data = self
            .graph
            .determine_raised_esurfaces_from_expression(&global_cff);

        let (raised_cut_data, raised_cut_stats) = RaisedCutData::new_from_esurface(
            &esurface_raised_data,
            &self.cut_esurface_id_map,
            &settings.evaluator,
        );

        self.derived_data.global_cff_expression = Some(global_cff);
        self.derived_data.raised_data = raised_cut_data;

        Ok(raised_cut_stats)
    }

    fn generate_cuts(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
        settings: &GenerationSettings,
    ) -> Result<()> {
        info!("generating cuts for graph: {}", self.graph.name);

        let all_st_cuts = self.graph.all_st_cuts_for_cs(
            self.source_nodes.clone(),
            self.target_nodes.clone(),
            &self.graph.get_initial_state_tree().0,
        );

        info!("num s_t cuts: {}", all_st_cuts.len());

        let mut cuts: TiVec<CutId, CrossSectionCut> = all_st_cuts
            .into_iter()
            .map(|(left, cut, right)| CrossSectionCut { cut, left, right })
            .filter(|cut| cut.cut.nedges(&self.graph) > 1)
            .filter_map(
                |cut| match cut.is_valid_for_process(self, process_definition, model) {
                    Ok(true) => Some(Ok(cut)),
                    Ok(false) => None,
                    Err(e) => Some(Err(e)),
                },
            )
            .collect::<Result<_>>()?;

        cuts.sort_by(|a, b| a.cut.cmp(&b.cut));

        if !settings.force_cuts.is_empty() {
            let force_cuts_sorted = settings
                .force_cuts
                .clone()
                .into_iter()
                .map(|cut_edges| cut_edges.into_iter().sorted().collect_vec())
                .collect_vec();

            cuts.retain(|cut| {
                let edges_in_cut = self
                    .graph
                    .iter_edges_of(&cut.cut)
                    .map(|(_, _, e)| e.data.name.value.clone())
                    .sorted()
                    .collect_vec();

                force_cuts_sorted.contains(&edges_in_cut)
            });
        }

        self.cuts = cuts;

        info!(
            "found {} cuts for graph: {}",
            self.cuts.len(),
            self.graph.name
        );

        Ok(())
    }

    fn generate_esurface_cuts(&mut self) {
        debug!("generating esurfaces for cuts");

        let esurfaces: TiVec<CutId, Esurface> = self
            .cuts
            .iter()
            .map(|cut| {
                Esurface::new_from_cut_left(
                    &self.graph.underlying,
                    cut,
                    Some(&self.graph.initial_state_cut),
                )
            })
            .collect();

        debug!("generated esurfaces {:?}", esurfaces);

        self.cut_esurface = esurfaces;
    }

    pub(crate) fn build_parametric_integrand(
        &mut self,
        settings: &GenerationSettings,
        vakint: &Vakint,
    ) -> Result<()> {
        self.derived_data.cut_paramatric_integrand = self.build_integrand(settings, vakint)?;
        Ok(())
    }

    fn build_integrand(
        &mut self,
        settings: &GenerationSettings,
        vakint: &Vakint,
    ) -> Result<TiVec<RaisedCutId, ParametricIntegrands>> {
        let max_order = self
            .derived_data
            .raised_data
            .raised_cut_groups
            .iter()
            .map(|cut_group| cut_group.related_esurface_group.max_occurence)
            .max()
            .unwrap();

        self.graph
            .param_builder
            .initialize_t_derivatives(max_order - 1);

        let cuts = self
            .derived_data
            .raised_data
            .raised_cut_groups
            .iter()
            .map(|cuts| CutSet {
                residue_selector: ResidueSelector {
                    lu_cut: Some(cuts.related_esurface_group.clone()),
                    left_th_cut: None,
                    right_th_cut: None,
                },
                union: cuts
                    .cuts
                    .iter()
                    .map(|cut_id| self.cuts[*cut_id].cut.as_subgraph())
                    .reduce(|cut_1, cut_2| cut_1.union(&cut_2))
                    .unwrap_or_else(|| self.graph.empty_subgraph()),
            })
            .collect();

        let cut_structure = CutStructure { cuts };

        let cut_woods = CutWoods::new(cut_structure, &self.graph, &settings.uv);
        let valid_orientations: Vec<_> = self
            .derived_data
            .global_cff_expression
            .as_ref()
            .expect("global_cff_expression should have been created")
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect();

        let lu_prefactor = self.lu_prefactor_helper();

        let mut cut_forests = cut_woods.unfold(&self.graph);
        cut_forests.compute(&mut self.graph, vakint, &valid_orientations, &settings.uv)?;

        let parametric_integrands =
            cut_forests.orientation_parametric_exprs(&self.graph, &settings.uv)?;

        Ok(parametric_integrands
            .into_iter()
            .map(|integrand| integrand.map(|a| a * &lu_prefactor))
            .collect())
    }

    fn lu_prefactor_helper(&self) -> Atom {
        let loop_number = self.graph.cyclotomatic_number(&self.graph.full_filter())
            - self.graph.initial_state_cut.nedges(&self.graph);

        let loop_3 = loop_number as i64 * 3;
        let factors_of_pi = (Atom::num(2) * Atom::var(GS.pi)).pow(loop_3 - 1); // multiply with 2pi from energy conservation delta

        let tstar = Atom::var(GS.rescale_star);
        let tsrat_pow = tstar.pow(loop_3);
        let hfunction = Atom::var(GS.hfunction_lu_cut);
        tsrat_pow * hfunction / factors_of_pi
    }

    fn th_prefactor_helper(
        &self,
        subspace_loop_count: usize,
        is_on_right: bool,
        include_integrated: bool,
    ) -> Atom {
        let radius = if is_on_right {
            Atom::var(GS.radius_right)
        } else {
            Atom::var(GS.radius_left)
        };

        let radius_star = if is_on_right {
            Atom::var(GS.radius_star_right)
        } else {
            Atom::var(GS.radius_star_left)
        };

        let grad_eta = if is_on_right {
            Atom::var(GS.deta_right_th)
        } else {
            Atom::var(GS.deta_left_th)
        };

        let uv_damp_plus = if is_on_right {
            Atom::var(GS.uv_damp_plus_right)
        } else {
            Atom::var(GS.uv_damp_plus_left)
        };

        let uv_damp_minus = if is_on_right {
            Atom::var(GS.uv_damp_minus_right)
        } else {
            Atom::var(GS.uv_damp_minus_left)
        };

        let h_function = if is_on_right {
            Atom::var(GS.hfunction_right_th)
        } else {
            Atom::var(GS.hfunction_left_th)
        };

        let i = if is_on_right { -Atom::i() } else { Atom::i() };
        let pi = Atom::var(GS.pi);

        let jacobian_ratio = (&radius_star / &radius).pow(subspace_loop_count as i64 * 3 - 1);

        let local_prefactor = (&uv_damp_plus / (&radius - &radius_star)
            + &uv_damp_minus / (-&radius - &radius_star))
            / &grad_eta
            * &jacobian_ratio;

        let integrated_prefactor = if include_integrated {
            &h_function * &i * &pi * &jacobian_ratio / &grad_eta
        } else {
            Atom::new()
        };

        debug!(
            "th prefactor local: {}, integrated: {}",
            local_prefactor, integrated_prefactor
        );

        local_prefactor + integrated_prefactor
    }

    fn build_lmbs(&mut self) -> Result<()> {
        let mut lmbs: TiVec<LmbIndex, LoopMomentumBasis> = vec![].into();

        let externals: SuBitGraph = self.graph.empty_subgraph();
        let full_filter = self.graph.full_filter();
        let cut_graph = full_filter.subtract(&self.graph.initial_state_cut.right);

        for s in self.graph.all_spanning_forests_of(&cut_graph) {
            let mut lmb = self.graph.lmb_impl(&full_filter, &s, externals.clone())?;
            let mut exts = vec![];

            for i in lmb.loop_edges.iter() {
                let (_, p) = &self.graph[i];

                if self.graph.initial_state_cut.intersects(p) {
                    exts.push(*i);
                }
            }

            exts.sort();

            for e in exts {
                let mut loopid = None;
                for (l, s) in lmb.edge_signatures[e].internal.iter_enumerated() {
                    if s.is_non_zero() {
                        if loopid.is_none() {
                            loopid = Some(l);
                        } else {
                            panic!("external edge has multiple loop momenta")
                        }
                    }
                }
                lmb.put_loop_to_ext(loopid.unwrap());
            }
            //lmbs.push(self.graph.lmb_impl(&full_filter, &s, externals.clone()));
            lmbs.push(lmb);
        }

        let sorted_graph_loop_edges = self
            .graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .sorted()
            .collect_vec();

        let matching_lmb_index = lmbs.iter_enumerated().find_map(|(lmb_index, lmb)| {
            let sorted_loop_edges = lmb.loop_edges.iter().copied().sorted().collect_vec();
            (sorted_loop_edges == sorted_graph_loop_edges).then_some(lmb_index)
        });

        if let Some(matching_lmb_index) = matching_lmb_index {
            {
                let matching_lmb = &mut lmbs[matching_lmb_index];
                for (target_pos, target_edge) in
                    self.graph.loop_momentum_basis.loop_edges.iter_enumerated()
                {
                    let current_pos = matching_lmb
                        .loop_edges
                        .iter()
                        .find_position(|edge| *edge == target_edge)
                        .map(|(pos, _)| pos)
                        .unwrap();

                    if current_pos != target_pos.0 {
                        matching_lmb.swap_loops(current_pos.into(), target_pos);
                    }
                }
            }

            let matching_lmb_pos: usize = matching_lmb_index.into();
            if matching_lmb_pos != 0 {
                lmbs.raw[..=matching_lmb_pos].rotate_right(1);
            }
        } else {
            warn!(
                "Could not match current graph LMB against generated LMBs for graph {}",
                self.graph.name
            );
        }

        let sorted_graph_loop_edges = self
            .graph
            .loop_momentum_basis
            .loop_edges
            .iter()
            .copied()
            .sorted()
            .collect_vec();

        let matching_lmb_index = lmbs.iter_enumerated().find_map(|(lmb_index, lmb)| {
            let sorted_loop_edges = lmb.loop_edges.iter().copied().sorted().collect_vec();
            (sorted_loop_edges == sorted_graph_loop_edges).then_some(lmb_index)
        });

        if let Some(matching_lmb_index) = matching_lmb_index {
            {
                let matching_lmb = &mut lmbs[matching_lmb_index];
                for (target_pos, target_edge) in
                    self.graph.loop_momentum_basis.loop_edges.iter_enumerated()
                {
                    let current_pos = matching_lmb
                        .loop_edges
                        .iter()
                        .find_position(|edge| *edge == target_edge)
                        .map(|(pos, _)| pos)
                        .unwrap();

                    if current_pos != target_pos.0 {
                        matching_lmb.swap_loops(current_pos.into(), target_pos);
                    }
                }
            }

            let matching_lmb_pos: usize = matching_lmb_index.into();
            if matching_lmb_pos != 0 {
                lmbs.raw[..=matching_lmb_pos].rotate_right(1);
            }
        } else {
            warn!(
                "Could not match current graph LMB against generated LMBs for graph {}",
                self.graph.name
            );
        }

        self.derived_data.lmbs = Some(lmbs);
        Ok(())
    }

    fn build_multi_channeling_channels(&mut self, override_lmb_heuristics: bool) {
        let channels = self.graph.build_multi_channeling_channels(
            self.derived_data.lmbs.as_ref().unwrap(),
            override_lmb_heuristics,
        );

        self.derived_data.multi_channeling_setup = Some(channels)
    }

    fn build_threshold_counterterm(
        &mut self,
        settings: &GenerationSettings,
        vakint: &Vakint,
    ) -> Result<()> {
        // threshold enumeration as st cuts
        let all_possible_thresholds: TiVec<GlobalThresholdId, _> = {
            let mut unsorted = self.graph.all_st_cuts_for_cs(
                self.source_nodes.clone(),
                self.target_nodes.clone(),
                &self.graph.get_initial_state_tree().0,
            );
            unsorted.retain(|(_left, cut, _right)| cut.nedges(&self.graph) > 1);
            if settings.threshold_subtraction.skip_thresholds_that_are_cuts {
                unsorted.retain(|(_left, cut, _right)| {
                    !self.cuts.iter().any(|cs_cut| &cs_cut.cut == cut)
                });
            }

            unsorted.sort_by(|a, b| a.1.cmp(&b.1));
            unsorted.into()
        };

        let mut left_cut_threshold_data: TiVec<CutId, TiVec<LeftThresholdId, EsurfaceID>> =
            ti_vec![TiVec::new(); self.cuts.len()];

        let mut right_cut_threshold_data: TiVec<CutId, TiVec<RightThresholdId, EsurfaceID>> =
            ti_vec![TiVec::new(); self.cuts.len()];

        for (cut_id, cut) in self.cuts.iter_enumerated() {
            for (_threshold_id, (left_threshold_diagram, threshold_cut, right_threshold_diagram)) in
                all_possible_thresholds.iter_enumerated()
            {
                if &cut.cut == threshold_cut
                    && settings.threshold_subtraction.skip_thresholds_that_are_cuts
                {
                    continue;
                }

                // if the subgraph on the left of the threshold cut is a subgraph of the left amplitude, then the threshold is on the left of the cut
                if cut.left.includes(left_threshold_diagram) {
                    let sandwich = cut.left.intersection(right_threshold_diagram);

                    // now we must check that the threshold cuts a loop and that the sandwich is connected
                    if self.graph.underlying.cyclotomatic_number(&cut.left)
                        > self
                            .graph
                            .underlying
                            .cyclotomatic_number(left_threshold_diagram)
                        && self.graph.underlying.is_connected(&sandwich)
                        && !self
                            .graph
                            .is_always_pinch(&sandwich, &cut.cut, threshold_cut)
                    {
                        let cross_section_cut_for_threshold = CrossSectionCut {
                            cut: threshold_cut.clone(),
                            left: left_threshold_diagram.clone(),
                            right: right_threshold_diagram.clone(),
                        };

                        let threshold_esurface = Esurface::new_from_cut_left(
                            &self.graph.underlying,
                            &cross_section_cut_for_threshold,
                            Some(&self.graph.initial_state_cut),
                        );

                        let threshold_id = self
                            .graph
                            .surface_cache
                            .esurface_cache
                            .position(|esurface| esurface == &threshold_esurface)
                            .unwrap();

                        left_cut_threshold_data[cut_id].push(threshold_id);
                    }
                } else if cut.right.includes(right_threshold_diagram) {
                    let sandwich = cut.right.intersection(left_threshold_diagram);
                    if self.graph.underlying.cyclotomatic_number(&cut.right)
                        > self
                            .graph
                            .underlying
                            .cyclotomatic_number(right_threshold_diagram)
                        && self.graph.underlying.is_connected(&sandwich)
                        && !self
                            .graph
                            .is_always_pinch(&sandwich, &cut.cut, threshold_cut)
                    {
                        let cross_section_cut_for_threshold = CrossSectionCut {
                            cut: threshold_cut.clone(),
                            left: left_threshold_diagram.clone(),
                            right: right_threshold_diagram.clone(),
                        };

                        let threshold_esurface = Esurface::new_from_cut_left(
                            &self.graph.underlying,
                            &cross_section_cut_for_threshold,
                            Some(&self.graph.initial_state_cut),
                        );

                        // threshold_esurface.subspace_graph =
                        //     InternalSubGraph::cleaned_filter_pessimist(
                        //         cut.right.clone(),
                        //         &self.graph,
                        //     );

                        let threshold_id = self
                            .graph
                            .surface_cache
                            .esurface_cache
                            .position(|esurface| esurface == &threshold_esurface)
                            .unwrap();

                        right_cut_threshold_data[cut_id].push(threshold_id);
                    }
                }
            }
        }

        let mut left_raised_cut_threshold_data: TiVec<
            RaisedCutId,
            TiVec<LeftThresholdId, EsurfaceID>,
        > = TiVec::new();

        let mut right_raised_cut_threshold_data: TiVec<
            RaisedCutId,
            TiVec<RightThresholdId, EsurfaceID>,
        > = TiVec::new();

        for raised_cut_group in self.derived_data.raised_data.raised_cut_groups.iter() {
            let left_thresholds: Vec<EsurfaceID> = raised_cut_group
                .cuts
                .iter()
                .flat_map(|cut_id| left_cut_threshold_data[*cut_id].clone())
                .sorted()
                .dedup()
                .collect();

            let mut right_thresholds: Vec<EsurfaceID> = raised_cut_group
                .cuts
                .iter()
                .flat_map(|cut_id| right_cut_threshold_data[*cut_id].clone())
                .sorted()
                .dedup()
                .collect();

            right_thresholds.retain(|esurface_id| !left_thresholds.contains(esurface_id));

            left_raised_cut_threshold_data.push(left_thresholds.into());
            right_raised_cut_threshold_data.push(right_thresholds.into());
        }

        let mut cut_structure = vec![];

        for (raised_cut_id, raised_cut_group) in self
            .derived_data
            .raised_data
            .raised_cut_groups
            .iter_enumerated()
        {
            let left_thresholds = &left_raised_cut_threshold_data[raised_cut_id];
            let right_thresholds = &right_raised_cut_threshold_data[raised_cut_id];

            let cutcosky_cut_untion = raised_cut_group
                .cuts
                .iter()
                .map(|cut_id| self.cuts[*cut_id].cut.as_subgraph())
                .reduce(|a, b| a.union(&b))
                .unwrap_or(self.graph.empty_subgraph());

            for esurface_id in left_thresholds {
                let raised_esurface_group = RaisedEsurfaceGroup {
                    esurface_ids: vec![*esurface_id],
                    max_occurence: 1,
                };

                let esurface_cut_union = self.graph.surface_cache.esurface_cache[*esurface_id]
                    .energies
                    .iter()
                    .map(|edge_id| self.graph.get_edge_subgraph(*edge_id))
                    .fold(cutcosky_cut_untion.clone(), |acc, subgraph| {
                        acc.union(&subgraph)
                    });

                cut_structure.push(CutSet {
                    residue_selector: ResidueSelector {
                        lu_cut: Some(raised_cut_group.related_esurface_group.clone()),
                        left_th_cut: Some(raised_esurface_group.clone()),
                        right_th_cut: None,
                    },
                    union: esurface_cut_union,
                });
            }

            for esurface_id in right_thresholds {
                let raised_esurface_group = RaisedEsurfaceGroup {
                    esurface_ids: vec![*esurface_id],
                    max_occurence: 1,
                };

                let esurface_cut_union = self.graph.surface_cache.esurface_cache[*esurface_id]
                    .energies
                    .iter()
                    .map(|edge_id| self.graph.get_edge_subgraph(*edge_id))
                    .fold(cutcosky_cut_untion.clone(), |acc, subgraph| {
                        acc.union(&subgraph)
                    });

                cut_structure.push(CutSet {
                    residue_selector: ResidueSelector {
                        lu_cut: Some(raised_cut_group.related_esurface_group.clone()),
                        left_th_cut: None,
                        right_th_cut: Some(raised_esurface_group.clone()),
                    },
                    union: esurface_cut_union,
                });
            }

            for (left_esurface_id, right_esurface_id) in left_thresholds
                .iter()
                .cartesian_product(right_thresholds.iter())
            {
                let left_raised_esurface_group = RaisedEsurfaceGroup {
                    esurface_ids: vec![*left_esurface_id],
                    max_occurence: 1,
                };

                let right_raised_esurface_group = RaisedEsurfaceGroup {
                    esurface_ids: vec![*right_esurface_id],
                    max_occurence: 1,
                };

                let esurface_cut_union = self.graph.surface_cache.esurface_cache[*left_esurface_id]
                    .energies
                    .iter()
                    .map(|edge_id| self.graph.get_edge_subgraph(*edge_id))
                    .fold(cutcosky_cut_untion.clone(), |acc, subgraph| {
                        acc.union(&subgraph)
                    })
                    .union(
                        &self.graph.surface_cache.esurface_cache[*right_esurface_id]
                            .energies
                            .iter()
                            .map(|edge_id| self.graph.get_edge_subgraph(*edge_id))
                            .fold(
                                self.graph.empty_subgraph::<SuBitGraph>(),
                                |acc, subgraph| acc.union(&subgraph),
                            ),
                    );

                cut_structure.push(CutSet {
                    residue_selector: ResidueSelector {
                        lu_cut: Some(raised_cut_group.related_esurface_group.clone()),
                        left_th_cut: Some(left_raised_esurface_group.clone()),
                        right_th_cut: Some(right_raised_esurface_group.clone()),
                    },
                    union: esurface_cut_union,
                });
            }
        }

        let cut_structure = CutStructure {
            cuts: cut_structure,
        };

        let cut_woods = CutWoods::new(cut_structure, &self.graph, &settings.uv);
        let mut cut_forests = cut_woods.unfold(&self.graph);
        let valid_orientations: Vec<_> = self
            .derived_data
            .global_cff_expression
            .as_ref()
            .expect("global_cff_expression should have been created")
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect();

        cut_forests.compute(&mut self.graph, vakint, &valid_orientations, &settings.uv)?;

        let mut threshold_counterterms = cut_forests
            .orientation_parametric_exprs(&self.graph, &settings.uv)?
            .into_iter();

        let lu_prefactor = self.lu_prefactor_helper();

        let mut result = TiVec::<RaisedCutId, LUCounterTermData>::new();
        for (raised_cut_id, _raised_cut_group) in self
            .derived_data
            .raised_data
            .raised_cut_groups
            .iter_enumerated()
        {
            let (left_subspace, right_subspace) = &self.derived_data.subspace_data[raised_cut_id];

            let th_prefactor_left = self.th_prefactor_helper(
                left_subspace.loopcount(),
                false,
                !settings.threshold_subtraction.disable_integrated_ct,
            );

            let th_prefactor_right = self.th_prefactor_helper(
                right_subspace.loopcount(),
                true,
                !settings.threshold_subtraction.disable_integrated_ct,
            );

            let iterated_prefactor = &th_prefactor_left * &th_prefactor_right;

            let mut left_atoms = TiVec::<LeftThresholdId, _>::new();
            let mut right_atoms = TiVec::<RightThresholdId, _>::new();
            let mut iterated_atoms = vec![];

            for _ in 0..left_raised_cut_threshold_data[raised_cut_id].len() {
                left_atoms.push(
                    threshold_counterterms
                        .next()
                        .unwrap()
                        .map(|x| x * &th_prefactor_left * &lu_prefactor),
                );
            }

            for _ in 0..right_raised_cut_threshold_data[raised_cut_id].len() {
                right_atoms.push(
                    threshold_counterterms
                        .next()
                        .unwrap()
                        .map(|x| x * &th_prefactor_right * &lu_prefactor),
                );
            }

            for _ in 0..(left_raised_cut_threshold_data[raised_cut_id].len()
                * right_raised_cut_threshold_data[raised_cut_id].len())
            {
                iterated_atoms.push(
                    threshold_counterterms
                        .next()
                        .unwrap()
                        .map(|x| x * &iterated_prefactor * &lu_prefactor),
                );
            }

            let iterated_collection = IteratedCtCollection {
                data: iterated_atoms,
                num_left_thresholds: left_atoms.len(),
            };

            let counterterm_data = LUCounterTermData {
                left_thresholds: left_raised_cut_threshold_data[raised_cut_id].clone(),
                right_thresholds: right_raised_cut_threshold_data[raised_cut_id].clone(),
                left_atoms,
                right_atoms,
                iterated: iterated_collection,
            };
            result.push(counterterm_data);
        }

        self.derived_data.threshold_counterterms = result;

        Ok(())
    }

    fn build_subspace_data(&mut self) -> Result<()> {
        let all_lmbs = self.derived_data.lmbs.as_ref().unwrap();

        let subspace_data = self
            .derived_data
            .raised_data
            .raised_cut_groups
            .iter()
            .map(|cut_group| {
                let valid_subspace_lmbs = all_lmbs
                    .iter_enumerated()
                    .filter_map(|(index, lmb)| {
                        let mut edges_in_cut = self
                            .graph
                            .underlying
                            .iter_edges_of(&self.cuts[cut_group.cuts[0]].cut)
                            .map(|(_, e, _)| e)
                            .collect_vec();

                        edges_in_cut.retain(|e| !lmb.loop_edges.contains(e));
                        if edges_in_cut.len() == 1 {
                            Some(index)
                        } else {
                            None
                        }
                    })
                    .collect_vec();

                let left_subgraphs = cut_group
                    .cuts
                    .iter()
                    .map(|cut_id| self.cuts[*cut_id].left.clone())
                    .sorted_by(|subgraph_a, subgraph_b| {
                        let n_edges_a = subgraph_a.nedges(&self.graph);
                        let n_edges_b = subgraph_b.nedges(&self.graph);
                        n_edges_a.cmp(&n_edges_b)
                    })
                    .collect_vec();

                let right_subgraphs = cut_group
                    .cuts
                    .iter()
                    .map(|cut_id| self.cuts[*cut_id].right.clone())
                    .sorted_by(|subgraph_a, subgraph_b| {
                        let n_edges_a = subgraph_a.nedges(&self.graph);
                        let n_edges_b = subgraph_b.nedges(&self.graph);
                        n_edges_a.cmp(&n_edges_b)
                    })
                    .collect_vec();

                let smallest_left_subgraph = left_subgraphs.first().unwrap().clone();
                let smallest_right_subgraph = right_subgraphs.first().unwrap().clone();

                let mut possible_subspaces = valid_subspace_lmbs
                    .iter()
                    .map(|lmb_index| {
                        (
                            SubspaceData::new_with_user_selected_lmb(
                                smallest_left_subgraph.clone(),
                                *lmb_index,
                                &self.graph,
                                all_lmbs,
                            )
                            .unwrap(),
                            SubspaceData::new_with_user_selected_lmb(
                                smallest_right_subgraph.clone(),
                                *lmb_index,
                                &self.graph,
                                all_lmbs,
                            )
                            .unwrap(),
                        )
                    })
                    .collect_vec();

                possible_subspaces.sort_by_key(|(left, right)| {
                    (
                        left.iter_basis_edges(all_lmbs).collect_vec(),
                        right.iter_basis_edges(all_lmbs).collect_vec(),
                    )
                });

                Ok(possible_subspaces.first().unwrap().clone())
            })
            .collect::<Result<_>>()?;

        let _: () = self.derived_data.subspace_data = subspace_data;
        Ok(())
    }

    fn generate_term_for_graph(
        &self,
        _model: &Model,
        settings: &GlobalSettings,
    ) -> Result<(CrossSectionGraphTerm, GraphGenerationStats)> {
        CrossSectionGraphTerm::from_cross_section_graph(self, settings)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionDerivedData {
    pub orientations: Option<TiVec<OrientationID, EdgeVec<Orientation>>>,
    pub cut_paramatric_integrand: TiVec<RaisedCutId, ParametricIntegrands>,
    pub global_cff_expression: Option<CFFExpression<OrientationID>>,
    pub lmbs: Option<TiVec<LmbIndex, LoopMomentumBasis>>,
    pub multi_channeling_setup: Option<LmbMultiChannelingSetup>,
    pub threshold_counterterms: TiVec<RaisedCutId, LUCounterTermData>,
    pub subspace_data: TiVec<RaisedCutId, (SubspaceData, SubspaceData)>,
    pub raised_data: RaisedCutData,
}

#[derive(Clone, Encode, Decode, Debug)]
#[trait_decode(trait = GammaLoopContext)]
pub struct RaisedCutData {
    pub raised_cut_groups: TiVec<RaisedCutId, RaisedCutGroup>,
    pub dual_shapes: Vec<Vec<Vec<usize>>>,
    pub pass_two_evaluators: Vec<GenericEvaluator>,
}

#[derive(Clone, Encode, Decode, Debug)]
#[trait_decode(trait = GammaLoopContext)]
pub struct RaisedCutGroup {
    pub cuts: Vec<CutId>,
    pub related_esurface_group: RaisedEsurfaceGroup,
}

impl Default for RaisedCutData {
    fn default() -> Self {
        Self::new()
    }
}

impl RaisedCutData {
    pub fn new() -> Self {
        RaisedCutData {
            raised_cut_groups: TiVec::new(),
            dual_shapes: vec![],
            pass_two_evaluators: vec![],
        }
    }

    pub fn new_from_esurface(
        raised_esurface_data: &RaisedEsurfaceData,
        cut_esurface_map: &TiVec<CutId, EsurfaceID>,
        evaluator_settings: &EvaluatorSettings,
    ) -> (Self, GraphGenerationStats) {
        let mut stats = GraphGenerationStats::default();
        let reversed_map = cut_esurface_map
            .iter_enumerated()
            .map(|(cut_id, &esurface_id)| (esurface_id, cut_id))
            .collect::<HashMap<EsurfaceID, CutId>>();

        let mut groups = TiVec::new();

        for (_raised_esurface_id, raised_esurface_group) in
            raised_esurface_data.raised_groups.iter_enumerated()
        {
            if cut_esurface_map.contains(&raised_esurface_group.esurface_ids[0]) {
                let cuts = raised_esurface_group
                    .esurface_ids
                    .iter()
                    .map(|esurface_id| reversed_map[esurface_id])
                    .collect::<Vec<_>>();

                let raised_cut_group = RaisedCutGroup {
                    cuts,
                    related_esurface_group: raised_esurface_group.clone(),
                };

                groups.push(raised_cut_group);
            } else {
                continue;
            }
        }

        let global_max_occurence = groups
            .iter()
            .map(|group| group.related_esurface_group.max_occurence)
            .max()
            .unwrap_or_else(|| {
                println!("corrupted groups");
                panic!();
            });

        let dual_shapes = (1..global_max_occurence)
            .map(shape_for_t_derivatives)
            .collect();

        let pass_two_evaluators = (1..=global_max_occurence)
            .map(|i| {
                let evaluator_started = std::time::Instant::now();
                let evaluator = build_derivative_structure(i as u8, -1, evaluator_settings);
                stats.evaluator_symbolica_time += evaluator_started.elapsed();
                stats.evaluator_count += 1;
                evaluator
            })
            .collect();

        (
            Self {
                raised_cut_groups: groups,
                dual_shapes,
                pass_two_evaluators,
            },
            stats,
        )
    }
}

impl CrossSectionDerivedData {
    fn new_empty() -> Self {
        Self {
            orientations: None,
            global_cff_expression: None,
            cut_paramatric_integrand: TiVec::new(),
            lmbs: None,
            multi_channeling_setup: None,
            threshold_counterterms: TiVec::new(),
            subspace_data: TiVec::new(),
            raised_data: RaisedCutData::new(),
        }
    }
}

pub(crate) fn build_derivative_structure_atom(
    singularity_order: u8,
    laurent_coefficient: i8,
) -> Atom {
    assert!(
        laurent_coefficient <= -1,
        "only laurent coefficients up to -1 are supported"
    );

    assert!(
        singularity_order >= 1,
        "eta order must be at least 1, got {singularity_order}"
    );
    assert!(
        singularity_order >= -laurent_coefficient as u8,
        "eta order must be at least the negative of the laurent coefficient, got {singularity_order} for laurent coefficient {laurent_coefficient}"
    );

    let order = singularity_order as i32;
    let laurent_coefficient = laurent_coefficient as i32;
    let f = symbol!("f");

    let expansion = parse!("η(t)")
        .series(
            GS.rescale,
            Atom::var(GS.rescale_star),
            (order, 1).into(),
            true,
        )
        .unwrap()
        .to_atom()
        .replace(function!(symbol!("η"), GS.rescale_star))
        .level_range((0, Some(0)))
        .with(0);

    let mut expression_to_derive = function!(f, GS.rescale)
        * expansion.pow(-order)
        * (GS.rescale - GS.rescale_star).pow(order);

    for _ in 1..=(order + laurent_coefficient) {
        expression_to_derive = expression_to_derive.derivative(GS.rescale);
    }

    expression_to_derive = expression_to_derive
        .replace(GS.rescale - GS.rescale_star)
        .with(parse!("delta_t"));

    let polynomial_in_delta_t = expression_to_derive
        .series(symbol!("delta_t"), Atom::num(0), (0, 1).into(), true)
        .unwrap();

    let factorial_prefactor = (2..=(order + laurent_coefficient as i32)).product::<i32>();

    let mut expression_to_derive = polynomial_in_delta_t.to_atom() / Atom::num(factorial_prefactor);

    expression_to_derive = expression_to_derive
        .replace(GS.rescale)
        .with(GS.rescale_star);

    expression_to_derive
}

pub(crate) fn build_derivative_structure(
    singularity_order: u8,
    _laurent_coefficient: i8,
    evaluator_settings: &EvaluatorSettings,
) -> GenericEvaluator {
    let expression_to_derive =
        build_derivative_structure_atom(singularity_order, _laurent_coefficient);

    let params = params_for_derivative_order(singularity_order);

    GenericEvaluator::new_from_raw_params(
        [expression_to_derive],
        &params,
        &FunctionMap::default(),
        vec![],
        evaluator_settings.optimization_settings(),
        None,
        evaluator_settings,
    )
    .unwrap()
    .into_eager_only()
}

fn params_for_derivative_order(derivative_order: u8) -> Vec<Atom> {
    let f = symbol!("f");
    let eta = symbol!("η");

    let f_0 = function!(f, GS.rescale_star);
    let eta_1 = function!(eta, GS.rescale_star).derivative(GS.rescale_star);

    let mut f_parameters = vec![f_0.clone()];
    let mut eta_params = vec![eta_1.clone()];

    for _ in 2..=derivative_order {
        let next_f = f_parameters.last().unwrap().derivative(GS.rescale_star);
        f_parameters.push(next_f);

        let next_eta = eta_params.last().unwrap().derivative(GS.rescale_star);
        eta_params.push(next_eta);
    }

    let mut result = vec![];
    result.extend(f_parameters);
    result.extend(eta_params);
    result
}

#[cfg(test)]
mod tests {
    use std::{
        fs,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    fn fresh_temp_dir(name: &str) -> PathBuf {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let path = std::env::temp_dir().join(format!(
            "gammalooprs-{name}-{}-{unique}",
            std::process::id()
        ));
        fs::create_dir_all(&path).unwrap();
        path
    }

    #[test]
    fn cross_section_storage_path_stays_inside_process_folder() {
        let temp = fresh_temp_dir("cross-section-storage-path");
        let cross_section = super::CrossSection::new("NLO".to_string());

        assert_eq!(cross_section.storage_path(&temp), temp.join("NLO"));
        fs::remove_dir_all(temp).unwrap();
    }
}
