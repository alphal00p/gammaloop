use std::{
    fmt::{Display, Formatter},
    fs::{self, File},
    io::Write,
    ops::{Index, IndexMut},
    path::Path,
};

use ahash::HashMap;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;
use idenso::color::ColorSimplifier;
use itertools::Itertools;
use rayon::{
    ThreadPool,
    iter::{IntoParallelRefMutIterator, ParallelIterator},
};
use spenso::{algebra::algebraic_traits::IsZero, structure::concrete_index::ExpandedIndex};
use tracing::info;
use vakint::Vakint;

use crate::{
    GammaLoopContext, GammaLoopContextContainer,
    cff::{
        esurface::{RaisedEsurfaceData, RaisedEsurfaceGroup},
        expression::{CFFExpression, OrientationID},
    },
    define_index, disable,
    graph::{
        GraphGroup, GroupId, LMBext, LmbIndex, LoopMomentumBasis, cuts::CutSet,
        parse::complete_group_parsing,
    },
    integrands::process::{
        GenericEvaluator, LmbMultiChannelingSetup,
        cross_section_integrand::CrossSectionIntegrandData,
    },
    model::ArcParticle,
    momentum::sample::SubspaceData,
    numerator::symbolica_ext::AtomCoreExt,
    processes::{DotExportSettings, EvaluatorSettings},
    settings::{GlobalSettings, global::GenerationSettings, runtime::LockedRuntimeSettings},
    utils::{GS, hyperdual_utils::shape_for_t_derivatives},
    uv::{approx::CutStructure, forest::ParametricIntegrands, uv_graph::UVE, wood::CutWoods},
};
use eyre::{Context, eyre};
use linnet::half_edge::{
    involution::{EdgeIndex, Orientation},
    subgraph::{
        HedgeNode, Inclusion, InternalSubGraph, OrientedCut, SuBitGraph, SubGraphLike, SubSetOps,
    },
};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    evaluate::{FunctionMap, OptimizationSettings},
    function,
    id::Replacement,
    parse, symbol,
};
use tracing::{debug, warn};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cut_expression::{CFFCutsExpression, CutOrientationData},
        esurface::{Esurface, EsurfaceID},
    },
    graph::{ExternalConnection, FeynmanGraph, Graph},
    integrands::process::{
        ProcessIntegrand,
        cross_section_integrand::{CrossSectionGraphTerm, CrossSectionIntegrand},
    },
    model::Model,
    numerator::ParsingNet,
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

#[derive(Clone, Debug)]
struct CsAmplitudeCTDiagram {
    left_subgraph: SuBitGraph,
    threshold_cut: OrientedCut,
    right_subgraph: SuBitGraph,
    reversed_dangling_edges: Vec<EdgeIndex>,
    network: Option<ParsingNet>,
}

#[derive(Clone, Debug, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct LUCounterTermData {
    pub left_thresholds: TiVec<LeftThresholdId, Esurface>,
    pub right_thresholds: TiVec<RightThresholdId, Esurface>,
    pub left_atoms: TiVec<LeftThresholdId, Atom>,
    pub right_atoms: TiVec<RightThresholdId, Atom>,
    pub iterated: IteratedCtCollection<Atom>,
}

impl CsAmplitudeCTDiagram {
    fn get_tensor_network_cached(
        &mut self,
        graph: &Graph,
        lu_cut: &OrientedCut,
        vakint: (&Vakint, &vakint::VakintSettings),
        add_lu_cut_feynman_rules: bool,
        settings: &GenerationSettings,
        conjugate: bool,
    ) -> ParsingNet {
        disable! {
        if let Some(network) = &self.network {
            network.clone()
        } else {
            let all_cut_edges = graph
                .iter_edges_of(&self.threshold_cut)
                .chain(graph.iter_edges_of(lu_cut))
                .map(|(_, e, _)| e)
                .collect_vec();

            let left_orientations = get_orientations_from_subgraph(
                graph,
                &self.left_subgraph,
                &self.reversed_dangling_edges,
            )
            .into_iter()
            .map(|cff_graph| cff_graph.global_orientation)
            .filter(|or| settings.orientation_pattern.alt_filter(or))
            .collect::<TiVec<OrientationID, _>>();

            let right_orientations = get_orientations_from_subgraph(
                graph,
                &self.right_subgraph,
                &self.reversed_dangling_edges,
            )
            .into_iter()
            .map(|cff_graph| cff_graph.global_orientation)
            .filter(|or| settings.orientation_pattern.alt_filter(or))
            .collect::<TiVec<OrientationID, _>>();

            let left_wood = graph.wood(&self.left_subgraph);
            let right_wood = graph.wood(&self.right_subgraph);

            let mut left_forest = left_wood.unfold(graph, &graph.loop_momentum_basis);
            let mut right_forest = right_wood.unfold(graph, &graph.loop_momentum_basis);

            let post = PostProcessingSetup {
                constraint_data: None,
                rewrite_esurfaces: None,
            };

            left_forest
                .compute(
                    graph,
                    &graph.tree_edges,
                    &self.left_subgraph,
                    vakint,
                    &left_orientations,
                    &None,
                    &all_cut_edges,
                    &graph.get_edges_in_initial_state_cut(),
                    post.clone(),
                    &settings.uv,
                    conjugate,
                )
                .unwrap();

            right_forest
                .compute(
                    graph,
                    &graph.tree_edges,
                    &self.right_subgraph,
                    vakint,
                    &right_orientations,
                    &None,
                    &all_cut_edges,
                    &graph.get_edges_in_initial_state_cut(),
                    post.clone(),
                    &settings.uv,
                    conjugate,
                )
                .unwrap();

            let left_expr = left_forest.orientation_parametric_expr(
                Some(
                    &self
                        .threshold_cut
                        .right
                        .union(&self.threshold_cut.left)
                        .subtract(&lu_cut.right)
                        .subtract(&lu_cut.left),
                ),
                graph,
                settings.uv.add_sigma,
            );

            let cut_edges_for_right = if add_lu_cut_feynman_rules {
                Some(lu_cut.left.union(&lu_cut.right))
            } else {
                None
            };

            let right_expr = right_forest.orientation_parametric_expr(
                cut_edges_for_right.as_ref(),
                graph,
                settings.uv.add_sigma,
            );

            self.network = Some(left_expr * right_expr);
            self.network.clone().unwrap()
        }
        }
        todo!()
    }
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
        generation_settings: &GenerationSettings,
        generation_pool: &ThreadPool,
    ) -> Result<()> {
        generation_pool.install(|| {
            self.supergraphs.par_iter_mut().try_for_each(|supergraph| {
                supergraph.preprocess(model, process_definition, generation_settings)
            })
        })
    }

    pub fn build_integrand(
        &mut self,
        model: &Model,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        generation_pool: &ThreadPool,
    ) -> Result<()> {
        let mut terms = generation_pool.install(|| {
            self.supergraphs
                .par_iter_mut()
                .map(|sg| sg.generate_term_for_graph(model, global_settings))
                .collect::<Result<Vec<_>>>()
        })?;

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

        let cross_section_integrand = CrossSectionIntegrand {
            settings: runtime_default.into(),
            data: CrossSectionIntegrandData {
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
            },
        };

        self.integrand = Some(ProcessIntegrand::CrossSection(cross_section_integrand));
        Ok(())
    }

    pub fn compile(
        &mut self,
        path: impl AsRef<Path>,
        override_existing: bool,
        settings: &GlobalSettings,
        thread_pool: &ThreadPool,
    ) -> Result<()> {
        info!("Compiling cross section {}", self.name);
        let p = path.as_ref().join(format!("cs_{}", self.name));

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
            integrand.compile(&p, override_existing, settings, thread_pool)?;
        }
        Ok(())
    }

    pub fn save(&mut self, path: impl AsRef<Path>, override_existing: bool) -> Result<()> {
        let p = path.as_ref().join(&self.name);
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

    pub(crate) fn preprocess(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
        settings: &GenerationSettings,
    ) -> Result<()> {
        debug!("generating cuts");
        self.generate_cuts(model, process_definition, settings)?;
        debug!("generating esurfaces corresponding to cuts");
        self.generate_esurface_cuts();
        debug!("generating cff");
        self.generate_cff(settings)?;
        debug!("building lmbs");
        self.build_lmbs();
        debug!("building multi channeling channels");

        if self.graph.is_group_master {
            self.build_multi_channeling_channels();
        }

        let vk = crate::utils::vakint()?;
        debug!("building parametric integrand");
        self.build_parametric_integrand(settings, vk)?;
        //self.build_parametric_integrand_raised_cuts(settings)?;

        if settings.threshold_subtraction.enable_thresholds {
            debug!("building threshold counterterm");
            self.build_threshold_counterterm(settings, vk)?;
            self.build_subspace_data()?;
        }

        Ok(())
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

    fn generate_cff(&mut self, _settings: &GenerationSettings) -> Result<()> {
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

        let raised_cut_data =
            RaisedCutData::new_from_esurface(&esurface_raised_data, &self.cut_esurface_id_map);

        self.derived_data.global_cff_expression = Some(global_cff);
        self.derived_data.raised_data = raised_cut_data;

        Ok(())
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
                    .map(|(_, _, e)| e.data.name.clone())
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

    fn get_initial_state_tree_data(&self) -> (ParsingNet, Vec<Replacement>) {
        let (tree_structure, props) = self.graph.get_initial_state_tree();
        let mut prop_atoms = Atom::num(1);
        let mut replacements = vec![];

        let external_energy_atoms = self
            .graph
            .loop_momentum_basis
            .ext_edges
            .iter()
            .map(|e_id| GS.emr_mom(*e_id, Atom::from(ExpandedIndex::from_iter([0]))))
            .collect_vec();

        for edge_id in props.iter() {
            let emr = (0..4)
                .map(|mu| GS.emr_mom(*edge_id, Atom::from(ExpandedIndex::from_iter([mu]))))
                .collect_vec();

            let mass = self.graph[*edge_id].mass_atom();
            let prop_atom = Atom::num(1)
                / (&emr[0] * &emr[0]
                    - &emr[1] * &emr[1]
                    - &emr[2] * &emr[2]
                    - &emr[3] * &emr[3]
                    - &mass * &mass);

            let replaced_mom = self.graph.loop_momentum_basis.edge_signatures[*edge_id]
                .external
                .apply(&external_energy_atoms);

            let replacement_1 = Replacement::new(
                GS.emr_mom(*edge_id, Atom::from(ExpandedIndex::from_iter([0])))
                    .to_pattern(),
                replaced_mom.clone(),
            );

            let replacement_2 =
                Replacement::new(GS.ose(*edge_id).to_pattern(), replaced_mom.clone());

            replacements.push(replacement_1);
            replacements.push(replacement_2);
            prop_atoms *= prop_atom;
        }

        prop_atoms = prop_atoms.replace_multiple(&replacements);

        let initial_state_tree_expr = ColorSimplifier::wrap_color(
            &(self
                .graph
                .iter_edges_of(&tree_structure)
                .fold(Atom::num(1), |acc, (_, _, edge)| acc * &edge.data.num)
                * self
                    .graph
                    .iter_nodes_of(&tree_structure)
                    .fold(Atom::num(1), |acc, (_, _, vertex)| acc * vertex.get_num())
                * &prop_atoms),
            GS.color_wrap,
        )
        .parse_into_net()
        .unwrap();

        (initial_state_tree_expr, replacements)
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
                esurfaces: Some(cuts.related_esurface_group.clone()),
                union: cuts
                    .cuts
                    .iter()
                    .map(|cut_id| self.cuts[*cut_id].cut.as_subgraph())
                    .reduce(|cut_1, cut_2| cut_1.union(&cut_2))
                    .unwrap_or_else(|| self.graph.empty_subgraph()),
            })
            .collect();

        let cut_structure = CutStructure { cuts };

        let cut_woods = CutWoods::new(cut_structure, &self.graph, &settings.uv.vakint);

        let lu_prefactor = self.lu_prefactor_helper_new();

        let mut cut_forests = cut_woods.unfold(&self.graph);
        cut_forests.compute(&mut self.graph, vakint, &settings.uv)?;
        Ok(cut_forests
            .orientation_parametric_exprs(&self.graph, settings.uv.add_sigma)?
            .into_iter()
            .map(|integrand| integrand.map(|a| a * &lu_prefactor))
            .collect())
    }

    fn lu_prefactor_helper(&self) -> Atom {
        let loop_number = self.graph.cyclotomatic_number(&self.graph.full_filter())
            - self.graph.initial_state_cut.nedges(&self.graph);

        let loop_3 = loop_number as i64 * 3;
        let grad_eta = Atom::var(GS.deta_lu_cut);
        let factors_of_pi = (Atom::num(2) * Atom::var(GS.pi)).pow(loop_3 - 1); // multiply with 2pi from energy conservation delta

        let tstar = Atom::var(GS.rescale_star);
        let tsrat_pow = tstar.pow(loop_3);
        let hfunction = Atom::var(GS.hfunction_lu_cut);
        tsrat_pow * hfunction / factors_of_pi / grad_eta
    }

    fn lu_prefactor_helper_new(&self) -> Atom {
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

    fn build_lmbs(&mut self) {
        let mut lmbs: TiVec<LmbIndex, LoopMomentumBasis> = vec![].into();

        let externals: SuBitGraph = self.graph.empty_subgraph();
        let full_filter = self.graph.full_filter();
        let cut_graph = full_filter.subtract(&self.graph.initial_state_cut.right);

        for s in self.graph.all_spanning_forests_of(&cut_graph) {
            let mut lmb = self.graph.lmb_impl(&full_filter, &s, externals.clone());
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

        self.derived_data.lmbs = Some(lmbs)
    }

    fn build_multi_channeling_channels(&mut self) {
        let channels = self
            .graph
            .build_multi_channeling_channels(self.derived_data.lmbs.as_ref().unwrap());

        self.derived_data.multi_channeling_setup = Some(channels)
    }

    fn build_threshold_counterterm(
        &mut self,
        settings: &GenerationSettings,
        vakint: &Vakint,
    ) -> Result<()> {
        disable!(
        // thershold enumeration as st cuts
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

        let (initial_state_tree, replacements) = self.get_initial_state_tree_data();

        let global_num = self.graph.global_network();

        let mut counterterms = TiVec::<CutId, LUCounterTermData>::new();

        for (cut_id, cut) in self.cuts.iter_enumerated() {
            let mut thresholds_on_the_left = TiVec::<LeftThresholdId, CsAmplitudeCTDiagram>::new();
            let mut thresholds_on_the_right =
                TiVec::<RightThresholdId, CsAmplitudeCTDiagram>::new();

            let mut threshold_esurfaces_on_the_left = TiVec::<LeftThresholdId, Esurface>::new();
            let mut threshold_esurfaces_on_the_right = TiVec::<RightThresholdId, Esurface>::new();

            let reversed_edges_in_xs_cut = cut
                .cut
                .iter_edges(&self.graph.underlying)
                .filter_map(|(orientation, edge)| match orientation {
                    Orientation::Reversed => Some(
                        self.graph
                            .edge_name_to_index(edge.data.name.as_str())
                            .unwrap(),
                    ),
                    _ => None,
                })
                .collect_vec();

            for (_threshold_id, (left_threshold_diagram, threshold_cut, right_threshold_diagram)) in
                all_possible_thresholds.iter_enumerated()
            {
                if &cut.cut == threshold_cut {
                    continue;
                }

                let mut reversed_dangling_edges = reversed_edges_in_xs_cut.clone();

                threshold_cut
                    .iter_edges(&self.graph.underlying)
                    .for_each(|(orientation, edge)| {
                        if orientation == Orientation::Reversed {
                            let edge_index = self
                                .graph
                                .edge_name_to_index(edge.data.name.as_str())
                                .unwrap();

                            if !reversed_dangling_edges.contains(&edge_index) {
                                reversed_dangling_edges.push(edge_index);
                            }
                        }
                    });

                // if the subgraph on the left of the threshold cut is a subgraph of the left amplitude, then the threshold is on the left of the cut
                if cut.left.includes(left_threshold_diagram) {
                    // now we must check that the threshold cuts a loop
                    if self.graph.underlying.cyclotomatic_number(&cut.left)
                        > self
                            .graph
                            .underlying
                            .cyclotomatic_number(left_threshold_diagram)
                    {
                        let right_subgraph = right_threshold_diagram.subtract(&cut.right);

                        let ct_diagram = CsAmplitudeCTDiagram {
                            left_subgraph: left_threshold_diagram.clone(),
                            threshold_cut: threshold_cut.clone(),
                            right_subgraph,
                            reversed_dangling_edges,
                            network: None,
                        };

                        thresholds_on_the_left.push(ct_diagram);

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

                        //threshold_esurface.subspace_graph =
                        //    InternalSubGraph::cleaned_filter_pessimist(
                        //        cut.left.clone(),
                        //        &self.graph,
                        //    );

                        threshold_esurfaces_on_the_left.push(threshold_esurface);
                    }
                } else if cut.right.includes(right_threshold_diagram)
                    && self.graph.underlying.cyclotomatic_number(&cut.right)
                        > self
                            .graph
                            .underlying
                            .cyclotomatic_number(right_threshold_diagram)
                {
                    let left_subgraph = left_threshold_diagram.subtract(&cut.left);

                    let ct_diagram = CsAmplitudeCTDiagram {
                        left_subgraph,
                        threshold_cut: threshold_cut.clone(),
                        right_subgraph: right_threshold_diagram.clone(),
                        reversed_dangling_edges,
                        network: None,
                    };

                    thresholds_on_the_right.push(ct_diagram);

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

                    threshold_esurfaces_on_the_right.push(threshold_esurface);
                }
            }

            let left_counterterms = thresholds_on_the_left
                .iter_mut()
                .map(|ct_diagram| {
                    let left_ct = ct_diagram.get_tensor_network_cached(
                        &self.graph,
                        &cut.cut,
                        vakint,
                        true,
                        settings,
                        false,
                    );

                    let mut product = left_ct
                        * self.derived_data.tensor_network_cache[cut_id].1.clone()
                        * global_num.clone()
                        * initial_state_tree.clone();

                    product
                        .execute::<Sequential, SmallestDegree, _, _, _>(
                            TENSORLIB.read().unwrap().deref(),
                            FUN_LIB.deref(),
                        )
                        .unwrap();

                    let left_scalar: Atom = product
                        .result_scalar()
                        .with_context(|| "in building threshold counterterm left")?
                        .into();

                    let mut left_integrand = left_scalar
                        .unwrap_function(GS.color_wrap)
                        .simplify_color()
                        .replace(function!(GS.energy, W_.x_))
                        .with(function!(GS.ose, W_.x_))
                        .replace(function!(GS.theta, W_.x_).pow(Atom::var(W_.n_)))
                        .with(function!(GS.theta, W_.x_))
                        .expand_dots()
                        .unwrap();

                    for (_, edge_index, _) in self
                        .graph
                        .underlying
                        .iter_edges_of(&self.graph.initial_state_cut)
                    {
                        left_integrand = left_integrand.replace(GS.ose(edge_index)).with(
                            GS.emr_mom(edge_index, Atom::from(ExpandedIndex::from_iter([0]))),
                        );
                    }

                    let lu_prefactor = self.lu_prefactor_helper();

                    let left_loop_count = self.graph.underlying.cyclotomatic_number(&cut.left);

                    let th_prefactor = self.th_prefactor_helper(
                        left_loop_count,
                        false,
                        !settings.threshold_subtraction.disable_integrated_ct,
                    );
                    left_integrand = left_integrand.replace_multiple(&replacements);

                    Ok(left_integrand * lu_prefactor * th_prefactor)
                })
                .collect::<Result<TiVec<LeftThresholdId, Atom>>>()?;

            let right_counterterms = thresholds_on_the_right
                .iter_mut()
                .map(|ct_diagram| {
                    let right_ct = ct_diagram.get_tensor_network_cached(
                        &self.graph,
                        &cut.cut,
                        vakint,
                        false,
                        settings,
                        true,
                    );

                    let mut product = self.derived_data.tensor_network_cache[cut_id].0.clone()
                        * right_ct
                        * global_num.clone()
                        * initial_state_tree.clone();

                    product
                        .execute::<Sequential, SmallestDegree, _, _, _>(
                            TENSORLIB.read().unwrap().deref(),
                            FUN_LIB.deref(),
                        )
                        .unwrap();

                    let right_scalar: Atom = product
                        .result_scalar()
                        .with_context(|| "in building threshold counterterm right")?
                        .into();

                    let mut right_integrand = right_scalar
                        .unwrap_function(GS.color_wrap)
                        .simplify_color()
                        .replace(function!(GS.energy, W_.x_))
                        .with(function!(GS.ose, W_.x_))
                        .replace(function!(GS.theta, W_.x_).pow(Atom::var(W_.n_)))
                        .with(function!(GS.theta, W_.x_))
                        .expand_dots()
                        .unwrap();

                    for (_, edge_index, _) in self
                        .graph
                        .underlying
                        .iter_edges_of(&self.graph.initial_state_cut)
                    {
                        right_integrand = right_integrand.replace(GS.ose(edge_index)).with(
                            GS.emr_mom(edge_index, Atom::from(ExpandedIndex::from_iter([0]))),
                        );
                    }

                    let lu_prefactor = self.lu_prefactor_helper();

                    let right_loop_count = self.graph.underlying.cyclotomatic_number(&cut.right);

                    let th_prefactor = self.th_prefactor_helper(
                        right_loop_count,
                        true,
                        !settings.threshold_subtraction.disable_integrated_ct,
                    );

                    right_integrand = right_integrand.replace_multiple(&replacements);

                    Ok(right_integrand * lu_prefactor * th_prefactor)
                })
                .collect::<Result<TiVec<RightThresholdId, Atom>>>()?;

            let iterated_counterterms = thresholds_on_the_left
                .iter()
                .cartesian_product(&thresholds_on_the_right)
                .map(|(left_ct_diagram, right_ct_diagram)| {
                    let mut product = left_ct_diagram.network.clone().unwrap()
                        * right_ct_diagram.network.clone().unwrap()
                        * global_num.clone()
                        * initial_state_tree.clone();

                    product
                        .execute::<Sequential, SmallestDegree, _, _, _>(
                            TENSORLIB.read().unwrap().deref(),
                            FUN_LIB.deref(),
                        )
                        .unwrap();

                    let iterated_scalar: Atom = product
                        .result_scalar()
                        .with_context(|| "in building threshold counterterm iterated")?
                        .into();

                    let mut iterated_integrand = iterated_scalar
                        .unwrap_function(GS.color_wrap)
                        .simplify_color()
                        .replace(function!(GS.energy, W_.x_))
                        .with(function!(GS.ose, W_.x_))
                        .replace(function!(GS.theta, W_.x_).pow(Atom::var(W_.n_)))
                        .with(function!(GS.theta, W_.x_))
                        .expand_dots()
                        .unwrap();

                    for (_, edge_index, _) in self
                        .graph
                        .underlying
                        .iter_edges_of(&self.graph.initial_state_cut)
                    {
                        iterated_integrand = iterated_integrand.replace(GS.ose(edge_index)).with(
                            GS.emr_mom(edge_index, Atom::from(ExpandedIndex::from_iter([0]))),
                        );
                    }

                    let lu_prefactor = self.lu_prefactor_helper();
                    let left_prefactor = self.th_prefactor_helper(
                        self.graph.underlying.cyclotomatic_number(&cut.left),
                        false,
                        !settings.threshold_subtraction.disable_integrated_ct,
                    );
                    let right_prefactor = self.th_prefactor_helper(
                        self.graph.underlying.cyclotomatic_number(&cut.right),
                        true,
                        !settings.threshold_subtraction.disable_integrated_ct,
                    );

                    iterated_integrand = iterated_integrand.replace_multiple(&replacements);

                    Ok(iterated_integrand * lu_prefactor * left_prefactor * right_prefactor)
                })
                .collect::<Result<Vec<_>>>()?;

            let iterated_counterterm = IteratedCtCollection {
                data: iterated_counterterms,
                num_left_thresholds: left_counterterms.len(),
            };

            let lu_counterterm_atom = LUCounterTermData {
                left_thresholds: threshold_esurfaces_on_the_left,
                right_thresholds: threshold_esurfaces_on_the_right,
                left_atoms: left_counterterms,
                right_atoms: right_counterterms,
                iterated: iterated_counterterm,
            };

            counterterms.push(lu_counterterm_atom);
        }

        self.derived_data.threshold_counterterms = counterterms;);
        Ok(())
    }

    fn build_subspace_data(&mut self) -> Result<()> {
        let all_lmbs = self.derived_data.lmbs.as_ref().unwrap();

        let subspace_data = self
            .cuts
            .iter()
            .map(|cut| {
                let (subspace_lmb_index, lmb) = all_lmbs
                    .iter_enumerated()
                    .find(|(_index, lmb)| {
                        let mut edges_in_cut = self
                            .graph
                            .underlying
                            .iter_edges_of(&cut.cut)
                            .map(|(_, e, _)| e)
                            .collect_vec();

                        edges_in_cut.retain(|e| !lmb.loop_edges.contains(e));
                        edges_in_cut.len() == 1
                    })
                    .unwrap();

                println!("externals: {:?}", lmb.ext_edges);

                let left_subspace = SubspaceData::new_with_user_selected_lmb(
                    cut.left.clone(),
                    subspace_lmb_index,
                    &self.graph,
                    all_lmbs,
                )?;
                let right_subspace = SubspaceData::new_with_user_selected_lmb(
                    cut.right.clone(),
                    subspace_lmb_index,
                    &self.graph,
                    all_lmbs,
                )?;

                Ok((left_subspace, right_subspace))
            })
            .collect::<Result<_>>()?;

        let _: () = self.derived_data.subspace_data = subspace_data;
        Ok(())
    }

    fn generate_term_for_graph(
        &self,
        _model: &Model,
        settings: &GlobalSettings,
    ) -> Result<CrossSectionGraphTerm> {
        CrossSectionGraphTerm::from_cross_section_graph(self, settings)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionDerivedData {
    pub orientations: Option<TiVec<OrientationID, CutOrientationData>>,
    pub cut_paramatric_integrand: TiVec<RaisedCutId, ParametricIntegrands>,
    pub global_cff_expression: Option<CFFExpression<OrientationID>>,
    pub cff_expression: Option<CFFCutsExpression>,
    pub lmbs: Option<TiVec<LmbIndex, LoopMomentumBasis>>,
    pub multi_channeling_setup: Option<LmbMultiChannelingSetup>,
    pub tensor_network_cache: TiVec<CutId, (ParsingNet, ParsingNet)>,
    pub threshold_counterterms: TiVec<CutId, LUCounterTermData>,
    pub subspace_data: TiVec<CutId, (SubspaceData, SubspaceData)>,
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
    ) -> Self {
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
            .map(|i| shape_for_t_derivatives(i))
            .collect();

        let pass_two_evaluators = (1..=global_max_occurence)
            .map(|i| build_derivative_structure(i as u8))
            .collect();

        Self {
            raised_cut_groups: groups,
            dual_shapes,
            pass_two_evaluators,
        }
    }
}

impl CrossSectionDerivedData {
    fn new_empty() -> Self {
        Self {
            orientations: None,
            global_cff_expression: None,
            cff_expression: None,
            cut_paramatric_integrand: TiVec::new(),
            lmbs: None,
            multi_channeling_setup: None,
            tensor_network_cache: TiVec::new(),
            threshold_counterterms: TiVec::new(),
            subspace_data: TiVec::new(),
            raised_data: RaisedCutData::new(),
        }
    }
}

fn build_derivative_structure(order: u8) -> GenericEvaluator {
    let order = order as i32;
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

    for _ in 1..order {
        expression_to_derive = expression_to_derive.derivative(GS.rescale);
    }

    expression_to_derive = expression_to_derive
        .replace(GS.rescale - GS.rescale_star)
        .with(parse!("delta_t"));

    let polynomial_in_delta_t = expression_to_derive
        .series(symbol!("delta_t"), Atom::num(0), (0, 1).into(), true)
        .unwrap();

    let factorial_prefactor = (2..order).product::<i32>();

    let mut expression_to_derive = polynomial_in_delta_t.to_atom() / Atom::num(factorial_prefactor);

    expression_to_derive = expression_to_derive
        .replace(GS.rescale)
        .with(GS.rescale_star);

    let params = params_for_derivative_order(order as u8);

    GenericEvaluator::new_from_raw_params(
        [expression_to_derive],
        &params,
        &FunctionMap::default(),
        vec![],
        OptimizationSettings::default(),
        None,
        &EvaluatorSettings::default(),
    )
    .unwrap()
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
