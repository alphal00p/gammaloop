use std::{
    collections::{BTreeMap, BTreeSet},
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
use idenso::{color::ColorSimplifier, metric::MetricSimplifier};
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
        expression::{
            OrientationID, ThreeDExpression, normalize_cut_edge_support_with_raised_edge_groups,
            normalize_three_d_expression_cut_support_with_raised_edge_groups,
            remove_ltd_global_contact_completions_from_local_residue,
            select_lu_cut_residue_for_representation,
        },
        surface::GammaLoopLinearEnergyExpr,
    },
    define_index,
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
    settings::{
        GlobalSettings,
        global::{GenerationSettings, ThreeDRepresentation},
        runtime::LockedRuntimeSettings,
    },
    utils::{
        GS, W_, external_energy_atom_from_index, hyperdual_utils::shape_for_t_derivatives,
        symbolica_ext::IsNeg,
    },
    uv::{
        approx::CutStructure,
        forest::{ParametricIntegrands, cff_explicit_sum_needs_outer_orientation_projection},
        wood::CutWoods,
    },
};
use eyre::{Context, eyre};
use linnet::half_edge::{
    involution::{EdgeIndex, EdgeVec, Orientation},
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
use three_dimensional_reps::{RepresentationMode, surface::LinearEnergyExpr};
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

fn numerator_fingerprint(numerator: &str) -> String {
    use std::hash::{DefaultHasher, Hash, Hasher};
    let mut hasher = DefaultHasher::new();
    numerator.hash(&mut hasher);
    format!("{:016x}", hasher.finish())
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
        global_settings: &GlobalSettings,
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
                explicit_orientation_sum_only: global_settings
                    .generation
                    .explicit_orientation_sum_only,
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
    pub(crate) fn positive_energy_cut_orientation_sign(&self, graph: &Graph) -> Result<i64> {
        let cut_edges = self.cut.iter_edges(graph).collect_vec();
        if cut_edges.is_empty() {
            return Err(eyre!(
                "Cannot determine Cutkosky positive-energy orientation for an empty cut in graph {}",
                graph.name
            ));
        }

        let mut cut_orientation_sign = 1;
        for (orientation, edge_data) in &cut_edges {
            match *orientation {
                Orientation::Default => {}
                Orientation::Reversed => cut_orientation_sign *= -1,
                Orientation::Undirected => {
                    return Err(eyre!(
                        "Cannot determine Cutkosky positive-energy orientation for undirected edge {} in graph {}",
                        edge_data.data.name,
                        graph.name
                    ));
                }
            }
        }
        Ok(cut_orientation_sign)
    }

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

fn loop_signature_row(lmb: &LoopMomentumBasis, edge_id: EdgeIndex, row_sign: i64) -> Vec<i64> {
    let signature = lmb.edge_signatures[edge_id].internal.to_momtrop_format();
    signature
        .iter()
        .map(|entry| row_sign * *entry as i64)
        .collect_vec()
}

fn integer_determinant_sign(rows: &[Vec<i64>]) -> Option<i64> {
    let dimension = rows.len();
    if !rows.iter().all(|row| row.len() == dimension) {
        return None;
    }
    if dimension == 0 {
        return Some(1);
    }

    let mut matrix = rows
        .iter()
        .map(|row| row.iter().map(|entry| i128::from(*entry)).collect_vec())
        .collect_vec();
    let mut determinant_sign = 1_i128;
    let mut previous_pivot = 1_i128;

    for column in 0..dimension - 1 {
        let pivot_row = (column..dimension).find(|row| matrix[*row][column] != 0)?;
        if pivot_row != column {
            matrix.swap(column, pivot_row);
            determinant_sign = -determinant_sign;
        }

        let pivot = matrix[column][column];
        for row in column + 1..dimension {
            for col in column + 1..dimension {
                matrix[row][col] = (matrix[row][col] * pivot
                    - matrix[row][column] * matrix[column][col])
                    / previous_pivot;
            }
            matrix[row][column] = 0;
        }
        previous_pivot = pivot;
    }

    let determinant = determinant_sign * matrix[dimension - 1][dimension - 1];
    Some(determinant.signum() as i64)
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
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
    ) -> Result<GraphGenerationStats> {
        global_settings.ensure_step_iii_pending_options_are_supported()?;
        let settings = &global_settings.generation;
        if global_settings.three_d_representation == ThreeDRepresentation::Ltd
            && settings.uv.subtract_uv
            && !settings.uv.local_uv_cts_from_expanded_4d_integrands
        {
            return Err(eyre!(
                "`global.3d_representation = LTD` with local UV counterterms from 3D expansions is not supported for cross-section supergraphs; set `global.generation.uv.local_uv_cts_from_expanded_4d_integrands = true` to use the representation-neutral 4D-expanded local UV construction"
            ));
        }
        let preprocess_started = std::time::Instant::now();
        let mut stats = GraphGenerationStats::default();
        self.apply_spin_sum(model, settings, &runtime_default)?;
        debug!("generating cuts");
        self.generate_cuts(model, process_definition, settings)?;
        debug!("generating esurfaces corresponding to cuts");
        self.generate_esurface_cuts();
        debug!("generating 3D expression");
        stats.merge_in_place(&self.build_3d_expression(global_settings)?);
        debug!("building lmbs");
        self.build_lmbs()?;
        debug!("building multi channeling channels");

        if self.graph.is_group_master {
            self.build_multi_channeling_channels(settings.override_lmb_heuristics);
        }

        let vk = crate::utils::vakint()?;
        debug!("building parametric integrand");
        self.build_parametric_integrand(global_settings, vk)?;
        //self.build_parametric_integrand_raised_cuts(settings)?;

        if settings.threshold_subtraction.enable_thresholds {
            debug!("building threshold counterterm");
            self.build_subspace_data()?;
            self.build_threshold_counterterm(settings, global_settings.three_d_representation, vk)?;
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

    fn build_3d_expression(&mut self, settings: &GlobalSettings) -> Result<GraphGenerationStats> {
        let canonize_esurface = self
            .graph
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let representation = settings.three_d_representation;
        let options = self
            .graph
            .production_3d_expression_options(representation, &settings.generation)?;
        let mut global_expression =
            self.graph
                .generate_3d_expression_for_integrand(&[], &canonize_esurface, &options)?;

        let raised_edge_groups = self.graph.get_raised_edge_groups();
        normalize_three_d_expression_cut_support_with_raised_edge_groups(
            &mut global_expression,
            &raised_edge_groups,
        );
        let generated_esurfaces = &global_expression.surfaces.esurface_cache;
        let normalized_generated_esurfaces =
            self.normalized_generated_esurfaces(generated_esurfaces, &raised_edge_groups);
        let cut_esurface_map = self
            .cut_esurface
            .iter()
            .map(|esurface| {
                self.find_generated_esurface_id(
                    esurface,
                    generated_esurfaces,
                    &normalized_generated_esurfaces,
                    &raised_edge_groups,
                    "Cutkosky cut",
                )
            })
            .collect::<Result<_>>()?;

        self.cut_esurface_id_map = cut_esurface_map;

        let esurface_raised_data = self
            .graph
            .determine_raised_esurfaces_from_expression(&global_expression);

        let (raised_cut_data, raised_cut_stats) = RaisedCutData::new_from_esurface(
            &esurface_raised_data,
            &self.cut_esurface_id_map,
            &settings.generation.evaluator,
        );

        let cff_projection = if representation == ThreeDRepresentation::Ltd {
            let mut cff_options = self.graph.production_3d_expression_options(
                ThreeDRepresentation::Cff,
                &settings.generation,
            )?;
            cff_options.energy_degree_bounds.clear();
            let mut cff_expression = self.graph.generate_3d_expression_for_integrand(
                &[],
                &canonize_esurface,
                &cff_options,
            )?;
            normalize_three_d_expression_cut_support_with_raised_edge_groups(
                &mut cff_expression,
                &raised_edge_groups,
            );

            let cff_generated_esurfaces = &cff_expression.surfaces.esurface_cache;
            let cff_normalized_generated_esurfaces =
                self.normalized_generated_esurfaces(cff_generated_esurfaces, &raised_edge_groups);
            let cff_cut_esurface_map = self
                .cut_esurface
                .iter()
                .map(|esurface| {
                    self.find_generated_esurface_id(
                        esurface,
                        cff_generated_esurfaces,
                        &cff_normalized_generated_esurfaces,
                        &raised_edge_groups,
                        "CFF projection Cutkosky cut",
                    )
                })
                .collect::<Result<_>>()?;

            let cff_esurface_raised_data = self
                .graph
                .determine_raised_esurfaces_from_expression(&cff_expression);
            let (cff_raised_data, _) = RaisedCutData::new_from_esurface(
                &cff_esurface_raised_data,
                &cff_cut_esurface_map,
                &settings.generation.evaluator,
            );

            Some((cff_expression, cff_cut_esurface_map, cff_raised_data))
        } else {
            None
        };

        self.derived_data.global_three_d_expression = Some(global_expression);
        self.derived_data.raised_data = raised_cut_data;
        self.derived_data.cff_projection_three_d_expression = None;
        self.derived_data.cff_projection_cut_esurface_id_map = TiVec::new();
        self.derived_data.cff_projection_raised_data = RaisedCutData::new();
        if let Some((cff_expression, cff_cut_esurface_map, cff_raised_data)) = cff_projection {
            self.derived_data.cff_projection_three_d_expression = Some(cff_expression);
            self.derived_data.cff_projection_cut_esurface_id_map = cff_cut_esurface_map;
            self.derived_data.cff_projection_raised_data = cff_raised_data;
        }

        Ok(raised_cut_stats)
    }

    fn normalized_generated_esurfaces(
        &self,
        generated_esurfaces: &TiVec<EsurfaceID, Esurface>,
        raised_edge_groups: &[Vec<linnet::half_edge::involution::EdgeIndex>],
    ) -> TiVec<EsurfaceID, Esurface> {
        generated_esurfaces
            .iter()
            .map(|esurface| {
                Graph::normalize_esurface_with_raised_edge_groups(esurface, raised_edge_groups)
            })
            .collect()
    }

    fn find_generated_esurface_id(
        &self,
        esurface: &Esurface,
        generated_esurfaces: &TiVec<EsurfaceID, Esurface>,
        normalized_generated_esurfaces: &TiVec<EsurfaceID, Esurface>,
        raised_edge_groups: &[Vec<linnet::half_edge::involution::EdgeIndex>],
        context: impl Display,
    ) -> Result<EsurfaceID> {
        let normalized_esurface =
            Graph::normalize_esurface_with_raised_edge_groups(esurface, raised_edge_groups);
        if let Some(position) = generated_esurfaces
            .iter()
            .position(|generated| generated == esurface)
        {
            return Ok(position.into());
        }
        generated_esurfaces
            .iter()
            .zip(normalized_generated_esurfaces.iter())
            .position(|(_, normalized_generated)| normalized_generated == &normalized_esurface)
            .ok_or_else(|| {
                eyre!(
                    "{context} E-surface {esurface:?} for graph {} was not present in the generated 3D expression surface cache.\n\
                     Normalized requested E-surface: {:?}\n\
                     Cut E-surfaces: {:?}\nGenerated E-surfaces: {:?}\nGenerated H-surfaces: {:?}",
                    self.graph.name,
                    normalized_esurface,
                    self.cut_esurface,
                    generated_esurfaces,
                    self.graph.surface_cache.hsurface_cache,
                )
            })
            .map(Into::into)
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
        global_settings: &GlobalSettings,
        vakint: &Vakint,
    ) -> Result<()> {
        self.derived_data.cut_paramatric_integrand =
            self.build_integrand(global_settings, vakint)?;
        Ok(())
    }

    fn build_integrand(
        &mut self,
        global_settings: &GlobalSettings,
        vakint: &Vakint,
    ) -> Result<TiVec<RaisedCutId, ParametricIntegrands>> {
        let settings = &global_settings.generation;
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

        let cuts = self.cutsets_for_raised_cuts()?;

        if !settings.uv.subtract_uv {
            return self.build_direct_3d_original_integrand(
                &cuts,
                settings,
                global_settings.three_d_representation,
            );
        }

        let representation = global_settings.three_d_representation;
        let root_expression = self
            .derived_data
            .global_three_d_expression
            .as_ref()
            .expect("global 3D expression should have been created");
        let raised_cut_groups = self.derived_data.raised_data.raised_cut_groups.clone();

        let cuts = self.cutsets_for_raised_cut_groups(&raised_cut_groups)?;
        let cut_structure = CutStructure { cuts };

        let cut_woods = CutWoods::new(cut_structure, &self.graph, &settings.uv);
        let valid_orientations: Vec<_> = root_expression
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect();

        let lu_prefactor = self.lu_prefactor_helper();
        let representation_prefactor = self.three_d_representation_lu_prefactor(representation);
        let mut cut_forests = cut_woods.unfold(&self.graph);
        cut_forests.compute(
            &mut self.graph,
            vakint,
            &valid_orientations,
            settings,
            Some(root_expression),
            representation,
            settings.explicit_orientation_sum_only,
        )?;

        let parametric_integrands =
            cut_forests.orientation_parametric_exprs(&self.graph, &settings.uv)?;

        self.dump_parametric_integrands_if_requested(
            &parametric_integrands,
            representation,
            settings.uv.local_uv_cts_from_expanded_4d_integrands,
        );

        Ok(parametric_integrands
            .into_iter()
            .map(|integrand| {
                if cff_explicit_sum_needs_outer_orientation_projection(
                    &self.graph,
                    &integrand.cuts,
                    settings,
                    representation,
                ) && !integrand.explicitly_summed_orientations
                {
                    integrand.sum_orientations_explicitly(&valid_orientations)
                } else {
                    integrand
                }
            })
            .zip(self.derived_data.raised_data.raised_cut_groups.iter())
            .map(|(integrand, _raised_cut_group)| {
                integrand.map(|a| a * &lu_prefactor * &representation_prefactor)
            })
            .collect())
    }

    fn cutsets_for_raised_cuts(&self) -> Result<Vec<CutSet>> {
        self.cutsets_for_raised_cut_groups(&self.derived_data.raised_data.raised_cut_groups)
    }

    fn cutsets_for_raised_cut_groups(
        &self,
        raised_cut_groups: &TiVec<RaisedCutId, RaisedCutGroup>,
    ) -> Result<Vec<CutSet>> {
        raised_cut_groups
            .iter()
            .map(|cuts| {
                Ok(CutSet {
                    residue_selector: self
                        .residue_selector_for_raised_cut_group(cuts, None, None)?,
                    union: cuts
                        .cuts
                        .iter()
                        .map(|cut_id| self.cuts[*cut_id].cut.as_subgraph())
                        .reduce(|cut_1, cut_2| cut_1.union(&cut_2))
                        .unwrap_or_else(|| self.graph.empty_subgraph()),
                })
            })
            .collect()
    }

    fn residue_selector_for_raised_cut_group(
        &self,
        raised_cut_group: &RaisedCutGroup,
        left_th_cut: Option<RaisedEsurfaceGroup>,
        right_th_cut: Option<RaisedEsurfaceGroup>,
    ) -> Result<ResidueSelector> {
        let (lu_cut_edge_sets, cut_orientation_signs, simple_ltd_lu_cut_signs) =
            self.lu_cut_edge_sets_with_cutkosky_signs(raised_cut_group.cuts.iter().copied())?;
        let is_simple_lu_cut = raised_cut_group.related_esurface_group.max_occurence == 1;
        let simple_ltd_lu_cut_surface_family_sign = if is_simple_lu_cut {
            self.simple_ltd_lu_cut_surface_family_sign(raised_cut_group, &lu_cut_edge_sets)?
        } else {
            1
        };
        let ltd_lu_cut_esurface_signs = if is_simple_lu_cut {
            raised_cut_group
                .cut_esurface_ids
                .iter()
                .copied()
                .zip(simple_ltd_lu_cut_signs)
                .collect()
        } else {
            // Repeated/confluent LU residues still need the selected local
            // variables oriented relative to the positive-energy Cutkosky
            // direction. The additional repeated-pole derivative parity cannot
            // be represented as a sign of one selected denominator variable,
            // so it is applied once as a residue-basis prefactor below.
            raised_cut_group
                .cut_esurface_ids
                .iter()
                .copied()
                .zip(cut_orientation_signs.iter().copied())
                .collect()
        };
        let ltd_lu_cut_residue_prefactor_sign = if is_simple_lu_cut {
            simple_ltd_lu_cut_surface_family_sign
        } else {
            self.repeated_ltd_lu_cut_residue_prefactor_sign(
                &raised_cut_group.related_esurface_group,
                &cut_orientation_signs,
            )
        };
        Ok(ResidueSelector {
            lu_cut: Some(raised_cut_group.related_esurface_group.clone()),
            lu_cut_edge_sets,
            ltd_lu_cut_esurface_signs,
            ltd_lu_cut_residue_prefactor_sign,
            left_th_cut,
            right_th_cut,
        })
    }

    fn simple_ltd_lu_cut_surface_family_sign(
        &self,
        raised_cut_group: &RaisedCutGroup,
        lu_cut_edge_sets: &[Vec<EdgeIndex>],
    ) -> Result<i64> {
        let Some(cff_expression) = self.derived_data.cff_projection_three_d_expression.as_ref()
        else {
            return Ok(1);
        };

        let signs = raised_cut_group
            .cuts
            .iter()
            .copied()
            .zip(lu_cut_edge_sets.iter())
            .map(|(cut_id, cut_edge_set)| {
                let cff_cut_esurface_id = self
                    .derived_data
                    .cff_projection_cut_esurface_id_map
                    .get(cut_id)
                    .copied()
                    .ok_or_else(|| {
                        eyre!(
                            "Cannot determine CFF LU surface-family sign for cut {} in graph {}: no CFF projection cut E-surface was recorded",
                            cut_id.0,
                            self.graph.name
                        )
                    })?;
                let cff_raised_cut_group = self
                    .derived_data
                    .cff_projection_raised_data
                    .raised_cut_groups
                    .iter()
                    .find(|group| {
                        group.cuts.contains(&cut_id)
                            && group
                                .related_esurface_group
                                .esurface_ids
                                .contains(&cff_cut_esurface_id)
                    })
                    .ok_or_else(|| {
                        eyre!(
                            "Cannot determine CFF LU surface-family sign for cut {} in graph {}: no CFF projection raised-cut group contains E-surface {:?}",
                            cut_id.0,
                            self.graph.name,
                            cff_cut_esurface_id
                        )
                    })?;

                let residues = cff_expression.clone().select_esurface_residue_with_cut_edges(
                    &cff_raised_cut_group.related_esurface_group,
                    std::slice::from_ref(cut_edge_set),
                );
                let mut signs = BTreeSet::new();
                let mut variant_count = 0usize;
                for residue in residues {
                    for orientation in residue.orientations {
                        for variant in orientation.variants {
                            variant_count += 1;
                            if variant.prefactor == Atom::Zero {
                                continue;
                            }
                            signs.insert(if variant.prefactor.is_negative() {
                                -1
                            } else {
                                1
                            });
                        }
                    }
                }

                if signs.len() != 1 {
                    return Err(eyre!(
                        "Cannot determine a unique CFF LU surface-family sign for cut {} in graph {}: found signs {:?} across {} selected CFF variants",
                        cut_id.0,
                        self.graph.name,
                        signs,
                        variant_count
                    ));
                }

                Ok(*signs.iter().next().expect("one sign was checked above"))
            })
            .collect::<Result<BTreeSet<_>>>()?;

        if signs.len() != 1 {
            return Err(eyre!(
                "Cannot determine a unique CFF LU surface-family sign for graph {} and cuts {:?}: found signs {:?}",
                self.graph.name,
                raised_cut_group.cuts,
                signs
            ));
        }

        Ok(*signs.iter().next().expect("one sign was checked above"))
    }

    fn lu_cut_edge_sets_with_cutkosky_signs(
        &self,
        cut_ids: impl IntoIterator<Item = CutId>,
    ) -> Result<(
        Vec<Vec<linnet::half_edge::involution::EdgeIndex>>,
        Vec<i64>,
        Vec<i64>,
    )> {
        let raised_edge_groups = self.graph.get_raised_edge_groups();
        let edge_sets_and_signs = cut_ids
            .into_iter()
            .map(|cut_id| {
                let cut_edges = self.cuts[cut_id].cut.iter_edges(&self.graph).collect_vec();
                if cut_edges.is_empty() {
                    return Err(eyre!(
                        "Cannot build Cutkosky residue selector for an empty cut in graph {}",
                        self.graph.name
                    ));
                }
                for (orientation, edge_data) in &cut_edges {
                    if *orientation == Orientation::Undirected {
                        return Err(eyre!(
                            "Cannot build Cutkosky residue selector for undirected edge {} in graph {}",
                            edge_data.data.name,
                            self.graph.name
                        ));
                    }
                }

                let cut_edges = cut_edges
                    .iter()
                    .map(|(_, edge_data)| {
                        self.graph
                            .edge_name_to_index(&edge_data.data.name)
                            .expect("Cut edge should belong to the graph")
                    })
                    .sorted()
                    .collect_vec();
                let cut_orientation_sign =
                    self.cuts[cut_id].positive_energy_cut_orientation_sign(&self.graph)?;
                let (coordinate_basis_sign, apply_cut_orientation_sign) =
                    self.simple_ltd_lu_cut_selected_denominator_sign(cut_id)?;
                let selected_denominator_sign = if apply_cut_orientation_sign {
                    cut_orientation_sign * coordinate_basis_sign
                } else {
                    coordinate_basis_sign
                };
                Ok((
                    normalize_cut_edge_support_with_raised_edge_groups(
                        &cut_edges,
                        &raised_edge_groups,
                    ),
                    cut_orientation_sign,
                    selected_denominator_sign,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        let mut edge_sets = Vec::with_capacity(edge_sets_and_signs.len());
        let mut cut_orientation_signs = Vec::with_capacity(edge_sets_and_signs.len());
        let mut simple_ltd_lu_cut_signs = Vec::with_capacity(edge_sets_and_signs.len());
        for (edge_set, cut_orientation_sign, simple_ltd_lu_cut_sign) in edge_sets_and_signs {
            edge_sets.push(edge_set);
            cut_orientation_signs.push(cut_orientation_sign);
            simple_ltd_lu_cut_signs.push(simple_ltd_lu_cut_sign);
        }
        Ok((edge_sets, cut_orientation_signs, simple_ltd_lu_cut_signs))
    }

    fn simple_ltd_lu_cut_selected_denominator_sign(&self, cut_id: CutId) -> Result<(i64, bool)> {
        self.simple_ltd_lu_cut_coordinate_basis_sign(cut_id)
    }

    fn simple_ltd_lu_cut_coordinate_basis_sign(&self, cut_id: CutId) -> Result<(i64, bool)> {
        let all_lmbs = self.derived_data.lmbs.as_ref().ok_or_else(|| {
            eyre!(
                "Cannot determine LTD LU cut coordinate-basis sign for graph {} before LMBs are built",
                self.graph.name
            )
        })?;
        let cut = &self.cuts[cut_id];
        let mut cut_edges = self
            .graph
            .underlying
            .iter_edges_of(&cut.cut)
            .map(|(_, edge_id, _)| edge_id)
            .collect_vec();
        cut_edges.sort_unstable();
        let cut_edge_orientations = cut
            .cut
            .iter_edges(&self.graph)
            .map(|(orientation, edge_data)| {
                let edge_id = self
                    .graph
                    .edge_name_to_index(&edge_data.data.name)
                    .expect("Cut edge should belong to the graph");
                let sign = match orientation {
                    Orientation::Default => Ok(1),
                    Orientation::Reversed => Ok(-1),
                    Orientation::Undirected => Err(eyre!(
                        "Cannot determine LTD LU cut coordinate-basis sign for undirected edge {} in graph {}",
                        edge_data.data.name,
                        self.graph.name
                    )),
                }?;
                Ok((edge_id, sign))
            })
            .collect::<Result<BTreeMap<_, _>>>()?;

        let mut candidates = all_lmbs
            .iter_enumerated()
            .filter_map(|(lmb_index, lmb)| {
                let mut omitted_cut_edges = cut_edges.clone();
                omitted_cut_edges.retain(|edge_id| !lmb.loop_edges.contains(edge_id));
                if omitted_cut_edges.len() != 1 {
                    return None;
                }
                let left = SubspaceData::new_with_user_selected_lmb(
                    cut.left.clone(),
                    lmb_index,
                    &self.graph,
                    all_lmbs,
                )
                .ok()?;
                let right = SubspaceData::new_with_user_selected_lmb(
                    cut.right.clone(),
                    lmb_index,
                    &self.graph,
                    all_lmbs,
                )
                .ok()?;
                Some((lmb_index, omitted_cut_edges[0], left, right))
            })
            .collect_vec();
        candidates.sort_by_key(|(_, _, left, right)| {
            (
                left.iter_basis_edges(all_lmbs).collect_vec(),
                right.iter_basis_edges(all_lmbs).collect_vec(),
            )
        });
        let (lmb_index, omitted_cut_edge, left, right) = candidates.first().ok_or_else(|| {
            eyre!(
                "Cannot determine LTD LU cut coordinate-basis sign for cut {} in graph {}: no LMB resolves exactly one cut edge externally",
                cut_id.0,
                self.graph.name
            )
        })?;
        let lmb = &all_lmbs[*lmb_index];
        // If the selected simple LU residue has no left-branch loop coordinate
        // in its local basis, the coordinate determinant is not anchored in the
        // left-to-right positive-energy Cutkosky convention. The Cutkosky
        // edge-flow product then converts the selected LTD surface variable to
        // that common convention. Once a left-branch loop coordinate is present,
        // the determinant already carries this orientation information.
        let apply_cut_orientation_sign = left.loopcount() == 0;

        let mut rows = Vec::new();
        for edge_id in left.iter_basis_edges(all_lmbs) {
            rows.push(loop_signature_row(lmb, edge_id, 1));
        }
        for edge_id in cut_edges
            .iter()
            .copied()
            .filter(|edge_id| edge_id != omitted_cut_edge)
        {
            let row_sign = *cut_edge_orientations
                .get(&edge_id)
                .expect("Cut edge orientation was collected above");
            rows.push(loop_signature_row(lmb, edge_id, row_sign));
        }
        for edge_id in right.iter_basis_edges(all_lmbs) {
            rows.push(loop_signature_row(lmb, edge_id, 1));
        }

        let determinant_sign = integer_determinant_sign(&rows).ok_or_else(|| {
            eyre!(
                "Cannot determine LTD LU cut coordinate-basis sign for cut {} in graph {}: singular resolved LU basis {:?}",
                cut_id.0,
                self.graph.name,
                rows
            )
        })?;
        // When the cut leaves no loop coordinate on either branch, the LU
        // residue is resolved directly against the external energy-conservation
        // variable. The determinant above orients the remaining cut-edge row;
        // the selected E-surface variable is the complementary tree-tree
        // direction and carries the opposite parity.
        let determinant_sign = if left.loopcount() == 0 && right.loopcount() == 0 {
            -determinant_sign
        } else {
            determinant_sign
        };
        Ok((determinant_sign, apply_cut_orientation_sign))
    }

    fn repeated_ltd_lu_cut_residue_prefactor_sign(
        &self,
        raised_esurface_group: &RaisedEsurfaceGroup,
        cut_orientation_signs: &[i64],
    ) -> i64 {
        let derivative_parity = if (raised_esurface_group.max_occurence - 1).is_multiple_of(2) {
            1
        } else {
            -1
        };

        // A confluent LTD denominator can represent several Cutkosky
        // alternatives by the same generated E-surface. The selector can store
        // only one sign for that selected variable, so the product of the
        // positive-energy orientations of the collapsed alternatives is a
        // residue-basis prefactor.
        derivative_parity * cut_orientation_signs.iter().product::<i64>()
    }

    fn finalize_parametric_integrand_atom(&self, atom: Atom) -> Result<Atom> {
        Ok(atom
            .replace(GS.dim)
            .with(4)
            .simplify_color()
            .replace(GS.den(W_.a_, W_.b_, W_.c_, W_.d_))
            .with(W_.d_)
            .expand_dots()?
            .simplify_metrics()
            .collect_factors())
    }

    fn ltd_lu_residue_atoms_from_explicit_sum(
        &self,
        expression: &ThreeDExpression<OrientationID>,
        numerator: &Atom,
        raised_cut_group: &RaisedCutGroup,
        cutset: &CutSet,
        settings: &GenerationSettings,
    ) -> Result<Vec<Atom>> {
        let mut expression = expression.clone();
        remove_ltd_global_contact_completions_from_local_residue(&mut expression);
        let atom = self
            .graph
            .three_d_expression_parametric_atom_with_numerator_gs(
                &expression,
                numerator,
                RepresentationMode::Ltd,
                true,
                &settings.orientation_pattern,
            );

        let residue_parameter = symbol!("ltd_lu_residue_parameter");
        let residue_parameter_atom = Atom::var(residue_parameter);
        let max_occurrence = raised_cut_group.related_esurface_group.max_occurence;
        let mut coefficients = vec![Atom::Zero; max_occurrence];

        let selected_surfaces = if max_occurrence == 1 {
            cutset.residue_selector.ltd_lu_cut_esurface_signs.clone()
        } else {
            raised_cut_group
                .related_esurface_group
                .esurface_ids
                .iter()
                .copied()
                .map(|esurface_id| (esurface_id, 1))
                .collect_vec()
        };

        for (esurface_id, selected_sign) in &selected_surfaces {
            let esurface = expression.surfaces.esurface_cache[*esurface_id].clone();
            let Some((localized_external_edge, localized_external_sign)) =
                esurface.external_shift.first().copied()
            else {
                return Err(eyre!(
                    "cannot take an LTD LU residue for graph {} on E-surface {esurface_id:?}: no external-energy shift is available to parametrize the selected surface",
                    self.graph.name
                ));
            };
            if localized_external_sign.unsigned_abs() != 1 {
                return Err(eyre!(
                    "cannot take an LTD LU residue for graph {} on E-surface {esurface_id:?}: expected a unit external-energy coefficient, found {localized_external_sign}",
                    self.graph.name
                ));
            }

            let surface_without_localized_external = esurface
                .energies
                .iter()
                .fold(LinearEnergyExpr::zero(), |acc, edge_id| {
                    acc + LinearEnergyExpr::ose(*edge_id, 1)
                });
            let surface_without_localized_external = esurface
                .external_shift
                .iter()
                .fold(
                    surface_without_localized_external,
                    |acc, (external_edge, external_sign)| {
                        if *external_edge == localized_external_edge {
                            acc
                        } else {
                            acc + LinearEnergyExpr::external(*external_edge, *external_sign)
                        }
                    },
                )
                .to_atom_gs(&[]);
            let localized_external_replacement = (Atom::num(*selected_sign)
                * residue_parameter_atom.clone()
                - surface_without_localized_external)
                / Atom::num(localized_external_sign);
            let localized_atom = atom
                .clone()
                .replace(external_energy_atom_from_index(localized_external_edge))
                .with(localized_external_replacement)
                .expand();
            let series = localized_atom.series(residue_parameter, Atom::Zero, 0)
                .map_err(|error| eyre!(
                    "failed to extract LTD LU residue for graph {} on E-surface {esurface_id:?}: {error}",
                    self.graph.name
                ))?;
            for occurrence in 1..=max_occurrence {
                let occurrence_i64 = i64::try_from(occurrence).map_err(|_| {
                    eyre!(
                        "failed to extract LTD LU residue for graph {} on E-surface {esurface_id:?}: occurrence does not fit in i64",
                        self.graph.name
                    )
                })?;
                let selected_variable_parity = if occurrence.is_multiple_of(2) {
                    1
                } else {
                    *selected_sign
                };
                coefficients[occurrence - 1] += series.coefficient((-occurrence_i64, 1).into())
                    * Atom::num(selected_variable_parity);
            }
        }

        Ok(coefficients)
    }

    fn build_direct_3d_original_integrand(
        &mut self,
        cuts: &[CutSet],
        settings: &GenerationSettings,
        representation: ThreeDRepresentation,
    ) -> Result<TiVec<RaisedCutId, ParametricIntegrands>> {
        let expression = self
            .derived_data
            .global_three_d_expression
            .as_ref()
            .expect("global 3D expression should have been created");
        let numerator = self
            .graph
            .production_numerator_atom_for_full_3d_expression();
        let lu_prefactor = self.lu_prefactor_helper();

        let parametric_integrands = self
            .derived_data
            .raised_data
            .raised_cut_groups
            .iter()
            .zip(cuts.iter())
            .map(|(raised_cut_group, cutset)| {
                if representation == ThreeDRepresentation::Ltd {
                    let integrands = self
                        .ltd_lu_residue_atoms_from_explicit_sum(
                            expression,
                            &numerator,
                            raised_cut_group,
                            cutset,
                            settings,
                        )?
                        .into_iter()
                        .map(|atom| {
                            let ltd_lu_cut_residue_prefactor =
                                Atom::num(cutset.residue_selector.ltd_residue_prefactor_sign());
                            // The direct no-UV LTD path extracts the actual
                            // local-series coefficient in the selected
                            // Cutkosky variable. That coefficient already
                            // contains the local residue Jacobian convention;
                            // the global LTD/CFF representation bridge is only
                            // needed by paths that project residues
                            // combinatorially on the oriented 3D structure.
                            self.finalize_parametric_integrand_atom(
                                atom * &lu_prefactor * ltd_lu_cut_residue_prefactor,
                            )
                        })
                        .collect::<Result<Vec<_>>>()?;
                    return Ok(ParametricIntegrands {
                        integrands,
                        cuts: cutset.clone(),
                        explicitly_summed_orientations: true,
                    });
                }
                let residues = select_lu_cut_residue_for_representation(
                    expression.clone(),
                    &raised_cut_group.related_esurface_group,
                    &cutset.residue_selector.lu_cut_edge_sets,
                    &cutset.residue_selector.ltd_lu_cut_esurface_signs,
                    representation,
                );
                self.dump_selected_residues_if_requested(
                    &residues,
                    &raised_cut_group.related_esurface_group,
                    &cutset.residue_selector.lu_cut_edge_sets,
                    representation,
                    settings.uv.local_uv_cts_from_expanded_4d_integrands,
                    &raised_cut_group.cuts,
                );
                let integrands = residues
                    .into_iter()
                    .map(|residue| {
                        let atom = self
                            .graph
                            .three_d_expression_parametric_atom_with_numerator_gs(
                                &residue,
                                &numerator,
                                RepresentationMode::Cff,
                                true,
                                &settings.orientation_pattern,
                            );
                        self.finalize_parametric_integrand_atom(atom * &lu_prefactor)
                    })
                    .collect::<Result<Vec<_>>>()?;
                Ok(ParametricIntegrands {
                    integrands,
                    cuts: cutset.clone(),
                    explicitly_summed_orientations: true,
                })
            })
            .collect::<Result<TiVec<_, _>>>()?;
        self.dump_parametric_integrands_if_requested(
            &parametric_integrands,
            representation,
            settings.uv.local_uv_cts_from_expanded_4d_integrands,
        );
        Ok(parametric_integrands)
    }

    fn dump_selected_residues_if_requested(
        &self,
        residues: &[ThreeDExpression<OrientationID>],
        raised_esurface_group: &RaisedEsurfaceGroup,
        cut_edge_sets: &[Vec<EdgeIndex>],
        representation: ThreeDRepresentation,
        local_uv_from_expanded_4d: bool,
        cut_ids: &[CutId],
    ) {
        let Ok(dump_dir) = std::env::var("GAMMALOOP_DUMP_PARAMETRIC_INTEGRANDS") else {
            return;
        };
        let representation_label = match representation {
            ThreeDRepresentation::Cff => "cff",
            ThreeDRepresentation::Ltd => "ltd",
        };
        let mode = if local_uv_from_expanded_4d {
            "local4d"
        } else {
            "local3d"
        };
        let numerator = self
            .graph
            .production_numerator_atom_for_full_3d_expression()
            .to_canonical_string();
        let numerator_fingerprint = numerator_fingerprint(&numerator);
        let cut_id = cut_ids
            .first()
            .map(|cut_id| usize::from(*cut_id))
            .unwrap_or(0);
        let path = std::path::Path::new(&dump_dir).join(format!(
            "{}_{}_{}_{}_cut{}_selected_residues.txt",
            self.graph.name, representation_label, mode, numerator_fingerprint, cut_id
        ));
        let mut dump = String::new();
        use std::fmt::Write;
        writeln!(dump, "numerator {numerator}").ok();
        writeln!(dump, "raised_esurface_group {:?}", raised_esurface_group).ok();
        for cut_id in cut_ids {
            let cut = &self.cuts[*cut_id];
            let cut_edges = cut
                .cut
                .iter_edges(&self.graph)
                .map(|(orientation, edge_data)| {
                    let edge_index = self
                        .graph
                        .edge_name_to_index(&edge_data.data.name)
                        .expect("Cut edge should belong to the graph");
                    (edge_index, orientation)
                })
                .collect_vec();
            writeln!(dump, "cut {cut_id:?} oriented_edges {:?}", cut_edges).ok();
        }
        writeln!(dump, "cut_edge_sets {:?}", cut_edge_sets).ok();
        writeln!(dump, "residue_count {}", residues.len()).ok();
        for (residue_id, residue) in residues.iter().enumerate() {
            writeln!(dump, "residue {residue_id}").ok();
            writeln!(dump, "  esurfaces {:?}", residue.surfaces.esurface_cache).ok();
            writeln!(dump, "  hsurfaces {:?}", residue.surfaces.hsurface_cache).ok();
            writeln!(dump, "  orientations {}", residue.orientations.len()).ok();
            for (orientation_id, orientation) in residue.orientations.iter_enumerated() {
                writeln!(
                    dump,
                    "  orientation {} label {:?} variants {} unfolded {}",
                    usize::from(orientation_id),
                    orientation.data.label,
                    orientation.variants.len(),
                    orientation.num_unfolded_terms()
                )
                .ok();
                writeln!(
                    dump,
                    "    loop_energy_map {:?}",
                    orientation.loop_energy_map
                )
                .ok();
                writeln!(
                    dump,
                    "    edge_energy_map {:?}",
                    orientation.edge_energy_map
                )
                .ok();
                for (variant_id, variant) in orientation.variants.iter().enumerate() {
                    writeln!(
                        dump,
                        "    variant {variant_id} origin {:?} prefactor {} half_edges {:?} denominator_edges {:?} denominator_surface_signs {:?} denominator_edge_support_signs {:?} numerator_surfaces {:?} denominator {:?}",
                        variant.origin,
                        variant.prefactor.to_canonical_string(),
                        variant.half_edges,
                        variant.denominator_edges,
                        variant.denominator_surface_signs,
                        variant.denominator_edge_support_signs,
                        variant.numerator_surfaces,
                        variant.denominator
                    )
                    .ok();
                }
            }
        }
        std::fs::create_dir_all(&dump_dir).ok();
        std::fs::write(path, dump).ok();
    }

    fn dump_parametric_integrands_if_requested<'a>(
        &self,
        parametric_integrands: impl IntoIterator<Item = &'a ParametricIntegrands>,
        representation: ThreeDRepresentation,
        local_uv_from_expanded_4d: bool,
    ) {
        let Ok(dump_dir) = std::env::var("GAMMALOOP_DUMP_PARAMETRIC_INTEGRANDS") else {
            return;
        };
        let representation_label = match representation {
            ThreeDRepresentation::Cff => "cff",
            ThreeDRepresentation::Ltd => "ltd",
        };
        let mode = if local_uv_from_expanded_4d {
            "local4d"
        } else {
            "local3d"
        };
        let numerator = self
            .graph
            .production_numerator_atom_for_full_3d_expression()
            .to_canonical_string();
        let numerator_fingerprint = numerator_fingerprint(&numerator);
        let path = std::path::Path::new(&dump_dir).join(format!(
            "{}_{}_{}_{}_parametric_integrands.txt",
            self.graph.name, representation_label, mode, numerator_fingerprint
        ));
        let mut dump = String::new();
        use std::fmt::Write;
        writeln!(dump, "numerator {numerator}").ok();
        for (cut_id, integrands) in parametric_integrands.into_iter().enumerate() {
            writeln!(dump, "cut {cut_id}").ok();
            writeln!(dump, "  cut_set {:?}", integrands.cuts).ok();
            for (derivative, atom) in integrands.integrands.iter().enumerate() {
                writeln!(
                    dump,
                    "derivative {derivative}: {}",
                    atom.to_canonical_string()
                )
                .ok();
            }
        }
        std::fs::create_dir_all(&dump_dir).ok();
        std::fs::write(path, dump).ok();
    }

    fn dump_threshold_counterterms_if_requested<'a>(
        &self,
        threshold_counterterms: impl IntoIterator<Item = &'a LUCounterTermData>,
        representation: ThreeDRepresentation,
        local_uv_from_expanded_4d: bool,
    ) {
        let Ok(dump_dir) = std::env::var("GAMMALOOP_DUMP_THRESHOLD_COUNTERTERMS") else {
            return;
        };
        let representation_label = match representation {
            ThreeDRepresentation::Cff => "cff",
            ThreeDRepresentation::Ltd => "ltd",
        };
        let mode = if local_uv_from_expanded_4d {
            "local4d"
        } else {
            "local3d"
        };
        let numerator = self
            .graph
            .production_numerator_atom_for_full_3d_expression()
            .to_canonical_string();
        let numerator_fingerprint = numerator_fingerprint(&numerator);
        let path = std::path::Path::new(&dump_dir).join(format!(
            "{}_{}_{}_{}_threshold_counterterms.txt",
            self.graph.name, representation_label, mode, numerator_fingerprint
        ));
        let mut dump = String::new();
        use std::fmt::Write;
        writeln!(dump, "numerator {numerator}").ok();
        for (cut_id, counterterm_data) in threshold_counterterms.into_iter().enumerate() {
            writeln!(dump, "cut {cut_id}").ok();
            for (left_id, atoms) in counterterm_data.left_atoms.iter_enumerated() {
                writeln!(
                    dump,
                    "left {left_id:?} esurface {:?}",
                    counterterm_data.left_thresholds[left_id]
                )
                .ok();
                for (derivative, atom) in atoms.integrands.iter().enumerate() {
                    writeln!(
                        dump,
                        "derivative {derivative}: {}",
                        atom.to_canonical_string()
                    )
                    .ok();
                }
            }
            for (right_id, atoms) in counterterm_data.right_atoms.iter_enumerated() {
                writeln!(
                    dump,
                    "right {right_id:?} esurface {:?}",
                    counterterm_data.right_thresholds[right_id]
                )
                .ok();
                for (derivative, atom) in atoms.integrands.iter().enumerate() {
                    writeln!(
                        dump,
                        "derivative {derivative}: {}",
                        atom.to_canonical_string()
                    )
                    .ok();
                }
            }
            for (index, atoms) in counterterm_data.iterated.iter().enumerate() {
                writeln!(dump, "iterated {index}").ok();
                for (derivative, atom) in atoms.integrands.iter().enumerate() {
                    writeln!(
                        dump,
                        "derivative {derivative}: {}",
                        atom.to_canonical_string()
                    )
                    .ok();
                }
            }
        }
        std::fs::create_dir_all(&dump_dir).ok();
        std::fs::write(path, dump).ok();
    }

    fn forward_scattering_loop_count(&self) -> usize {
        let loop_number = self.graph.cyclotomatic_number(&self.graph.full_filter())
            - self.graph.initial_state_cut.nedges(&self.graph);

        loop_number
    }

    fn lu_prefactor_helper(&self) -> Atom {
        let loop_number = self.forward_scattering_loop_count();

        let loop_3 = loop_number as i64 * 3;
        let factors_of_pi = (Atom::num(2) * Atom::var(GS.pi)).pow(loop_3 - 1); // multiply with 2pi from energy conservation delta

        let tstar = Atom::var(GS.rescale_star);
        let tsrat_pow = tstar.pow(loop_3);
        let hfunction = Atom::var(GS.hfunction_lu_cut);
        tsrat_pow * hfunction / factors_of_pi
    }

    fn three_d_representation_lu_prefactor(&self, representation: ThreeDRepresentation) -> Atom {
        match representation {
            ThreeDRepresentation::Cff => Atom::num(1),
            ThreeDRepresentation::Ltd => self.ltd_lu_measure_parity_factor(),
        }
    }

    fn ltd_lu_measure_parity_factor(&self) -> Atom {
        // Cross-section LU residues are normalized in the CFF loop-measure
        // convention. LTD performs the same energy integrations in the dual
        // measure convention, leaving the representation bridge (-1)^(L-1)
        // for an L-loop forward-scattering graph. This is global to the
        // representation, independent of the selected Cutkosky cut.
        let parity = if self
            .forward_scattering_loop_count()
            .saturating_sub(1)
            .is_multiple_of(2)
        {
            1
        } else {
            -1
        };
        Atom::num(i64::from(parity))
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
        representation: ThreeDRepresentation,
        vakint: &Vakint,
    ) -> Result<()> {
        let root_expression = self
            .derived_data
            .global_three_d_expression
            .as_ref()
            .expect("threshold root 3D expression should have been created");

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

        let raised_edge_groups = self.graph.get_raised_edge_groups();
        let generated_esurfaces = &root_expression.surfaces.esurface_cache;
        let normalized_generated_esurfaces =
            self.normalized_generated_esurfaces(generated_esurfaces, &raised_edge_groups);

        for (cut_id, cut) in self.cuts.iter_enumerated() {
            for (threshold_id, (left_threshold_diagram, threshold_cut, right_threshold_diagram)) in
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

                        let threshold_id = self.find_generated_esurface_id(
                            &threshold_esurface,
                            generated_esurfaces,
                            &normalized_generated_esurfaces,
                            &raised_edge_groups,
                            format!("left threshold {threshold_id:?} for Cutkosky cut {cut_id:?}"),
                        )?;
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

                        let threshold_id = self.find_generated_esurface_id(
                            &threshold_esurface,
                            generated_esurfaces,
                            &normalized_generated_esurfaces,
                            &raised_edge_groups,
                            format!("right threshold {threshold_id:?} for Cutkosky cut {cut_id:?}"),
                        )?;
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

                let esurface_cut_union = generated_esurfaces[*esurface_id]
                    .energies
                    .iter()
                    .map(|edge_id| self.graph.get_edge_subgraph(*edge_id))
                    .fold(cutcosky_cut_untion.clone(), |acc, subgraph| {
                        acc.union(&subgraph)
                    });

                cut_structure.push(CutSet {
                    residue_selector: self.residue_selector_for_raised_cut_group(
                        raised_cut_group,
                        Some(raised_esurface_group.clone()),
                        None,
                    )?,
                    union: esurface_cut_union,
                });
            }

            for esurface_id in right_thresholds {
                let raised_esurface_group = RaisedEsurfaceGroup {
                    esurface_ids: vec![*esurface_id],
                    max_occurence: 1,
                };

                let esurface_cut_union = generated_esurfaces[*esurface_id]
                    .energies
                    .iter()
                    .map(|edge_id| self.graph.get_edge_subgraph(*edge_id))
                    .fold(cutcosky_cut_untion.clone(), |acc, subgraph| {
                        acc.union(&subgraph)
                    });

                cut_structure.push(CutSet {
                    residue_selector: self.residue_selector_for_raised_cut_group(
                        raised_cut_group,
                        None,
                        Some(raised_esurface_group.clone()),
                    )?,
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

                let esurface_cut_union = generated_esurfaces[*left_esurface_id]
                    .energies
                    .iter()
                    .map(|edge_id| self.graph.get_edge_subgraph(*edge_id))
                    .fold(cutcosky_cut_untion.clone(), |acc, subgraph| {
                        acc.union(&subgraph)
                    })
                    .union(
                        &generated_esurfaces[*right_esurface_id]
                            .energies
                            .iter()
                            .map(|edge_id| self.graph.get_edge_subgraph(*edge_id))
                            .fold(
                                self.graph.empty_subgraph::<SuBitGraph>(),
                                |acc, subgraph| acc.union(&subgraph),
                            ),
                    );

                cut_structure.push(CutSet {
                    residue_selector: self.residue_selector_for_raised_cut_group(
                        raised_cut_group,
                        Some(left_raised_esurface_group.clone()),
                        Some(right_raised_esurface_group.clone()),
                    )?,
                    union: esurface_cut_union,
                });
            }
        }

        let cut_structure = CutStructure {
            cuts: cut_structure,
        };

        let cut_woods = CutWoods::new(cut_structure, &self.graph, &settings.uv);
        let mut cut_forests = cut_woods.unfold(&self.graph);
        let valid_orientations: Vec<_> = root_expression
            .orientations
            .iter()
            .map(|orientation| orientation.data.orientation.clone())
            .collect();

        cut_forests.compute(
            &mut self.graph,
            vakint,
            &valid_orientations,
            settings,
            Some(root_expression),
            representation,
            settings.explicit_orientation_sum_only,
        )?;

        let mut threshold_counterterms = cut_forests
            .orientation_parametric_exprs(&self.graph, &settings.uv)?
            .into_iter();

        let lu_prefactor = self.lu_prefactor_helper();
        let representation_prefactor = self.three_d_representation_lu_prefactor(representation);

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
                let atoms = threshold_counterterms.next().unwrap();
                let atoms = if cff_explicit_sum_needs_outer_orientation_projection(
                    &self.graph,
                    &atoms.cuts,
                    settings,
                    representation,
                ) && !atoms.explicitly_summed_orientations
                {
                    atoms.sum_orientations_explicitly(&valid_orientations)
                } else {
                    atoms
                };
                left_atoms.push(
                    atoms.map(|x| {
                        x * &th_prefactor_left * &lu_prefactor * &representation_prefactor
                    }),
                );
            }

            for _ in 0..right_raised_cut_threshold_data[raised_cut_id].len() {
                let atoms = threshold_counterterms.next().unwrap();
                let atoms = if cff_explicit_sum_needs_outer_orientation_projection(
                    &self.graph,
                    &atoms.cuts,
                    settings,
                    representation,
                ) && !atoms.explicitly_summed_orientations
                {
                    atoms.sum_orientations_explicitly(&valid_orientations)
                } else {
                    atoms
                };
                right_atoms.push(
                    atoms.map(|x| {
                        x * &th_prefactor_right * &lu_prefactor * &representation_prefactor
                    }),
                );
            }

            for _ in 0..(left_raised_cut_threshold_data[raised_cut_id].len()
                * right_raised_cut_threshold_data[raised_cut_id].len())
            {
                let atoms = threshold_counterterms.next().unwrap();
                let atoms = if cff_explicit_sum_needs_outer_orientation_projection(
                    &self.graph,
                    &atoms.cuts,
                    settings,
                    representation,
                ) && !atoms.explicitly_summed_orientations
                {
                    atoms.sum_orientations_explicitly(&valid_orientations)
                } else {
                    atoms
                };
                iterated_atoms.push(
                    atoms.map(|x| {
                        x * &iterated_prefactor * &lu_prefactor * &representation_prefactor
                    }),
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

        self.dump_threshold_counterterms_if_requested(
            &result,
            representation,
            settings.uv.local_uv_cts_from_expanded_4d_integrands,
        );
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
    pub global_three_d_expression: Option<ThreeDExpression<OrientationID>>,
    pub cff_projection_three_d_expression: Option<ThreeDExpression<OrientationID>>,
    pub cff_projection_cut_esurface_id_map: TiVec<CutId, EsurfaceID>,
    pub cff_projection_raised_data: RaisedCutData,
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
    pub cut_esurface_ids: Vec<EsurfaceID>,
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
        let mut reversed_map = HashMap::<EsurfaceID, Vec<CutId>>::default();
        for (cut_id, &esurface_id) in cut_esurface_map.iter_enumerated() {
            reversed_map.entry(esurface_id).or_default().push(cut_id);
        }

        let mut groups = TiVec::new();

        for (_raised_esurface_id, raised_esurface_group) in
            raised_esurface_data.raised_groups.iter_enumerated()
        {
            let cuts = raised_esurface_group
                .esurface_ids
                .iter()
                .filter_map(|esurface_id| reversed_map.get(esurface_id))
                .flatten()
                .copied()
                .unique()
                .collect::<Vec<_>>();

            if !cuts.is_empty() {
                let cuts = cuts
                    .into_iter()
                    .sorted_by_key(|cut_id| cut_id.0)
                    .collect_vec();
                let cut_esurface_ids = cuts
                    .iter()
                    .map(|cut_id| cut_esurface_map[*cut_id])
                    .collect_vec();
                let representative_esurface_id = cut_esurface_map[cuts[0]];
                let mut related_esurface_group = raised_esurface_group.clone();
                if let Some(position) = related_esurface_group
                    .esurface_ids
                    .iter()
                    .position(|esurface_id| *esurface_id == representative_esurface_id)
                {
                    related_esurface_group.esurface_ids.swap(0, position);
                }
                let raised_cut_group = RaisedCutGroup {
                    cuts,
                    cut_esurface_ids,
                    related_esurface_group,
                };

                groups.push(raised_cut_group);
            } else {
                continue;
            }
        }
        groups.sort_by_key(|group: &RaisedCutGroup| {
            (
                group.cuts.iter().map(|cut_id| cut_id.0).min(),
                group
                    .related_esurface_group
                    .esurface_ids
                    .iter()
                    .copied()
                    .min(),
            )
        });

        let global_max_occurence = groups
            .iter()
            .map(|group| group.related_esurface_group.max_occurence)
            .max()
            .expect("at least one raised cut group should be associated with cut E-surfaces");

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
            global_three_d_expression: None,
            cff_projection_three_d_expression: None,
            cff_projection_cut_esurface_id_map: TiVec::new(),
            cff_projection_raised_data: RaisedCutData::new(),
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
        .series(GS.rescale, Atom::var(GS.rescale_star), (order, 1))
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
        .series(symbol!("delta_t"), Atom::num(0), (0, 1))
        .unwrap();

    let factorial_prefactor = (2..=(order + laurent_coefficient)).product::<i32>();

    let mut expression_to_derive = polynomial_in_delta_t.to_atom() / Atom::num(factorial_prefactor);

    expression_to_derive = expression_to_derive
        .replace(GS.rescale)
        .with(GS.rescale_star);

    expression_to_derive
}

pub(crate) fn build_derivative_structure(
    singularity_order: u8,
    laurent_coefficient: i8,
    evaluator_settings: &EvaluatorSettings,
) -> GenericEvaluator {
    let expression_to_derive =
        build_derivative_structure_atom(singularity_order, laurent_coefficient);

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

pub(crate) fn params_for_derivative_order(derivative_order: u8) -> Vec<Atom> {
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
