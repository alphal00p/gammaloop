use std::{
    fmt::{Display, Formatter},
    fs::{self, File},
    io::Write,
    ops::{Deref, Index, IndexMut},
    path::Path,
};

// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use color_eyre::Result;

use idenso::color::ColorSimplifier;
use rayon::{
    ThreadPool,
    iter::{IntoParallelRefMutIterator, ParallelIterator},
};
use spenso::{
    network::{Sequential, SmallestDegree},
    structure::concrete_index::ExpandedIndex,
};
use tracing::info;
use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

use crate::{
    GammaLoopContext, GammaLoopContextContainer,
    cff::{
        cut_expression::SuperGraphOrientationID, expression::AmplitudeOrientationID,
        generation::get_orientations_from_subgraph,
    },
    define_index,
    gammaloop_integrand::{
        LmbMultiChannelingSetup, ParamBuilder, cross_section_integrand::CrossSectionIntegrandData,
    },
    graph::{
        GraphGroup, GroupId, LMBext, LmbIndex, LoopMomentumBasis, parse::complete_group_parsing,
    },
    model::ArcParticle,
    numerator::symbolica_ext::AtomCoreExt,
    settings::{GlobalSettings, global::GenerationSettings, runtime::LockedRuntimeSettings},
    utils::{FUN_LIB, GS, TENSORLIB, W_},
    uv::UltravioletGraph,
};
use eyre::{Context, eyre};
use itertools::Itertools;
use linnet::half_edge::{
    involution::{EdgeIndex, Orientation},
    subgraph::{
        HedgeNode, Inclusion, InternalSubGraph, OrientedCut, SuBitGraph, SubGraphLike, SubSetOps,
    },
};
use log::{debug, warn};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    function,
};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cut_expression::{CFFCutsExpression, CutOrientationData},
        esurface::{Esurface, EsurfaceID},
        generation::generate_cff_with_cuts,
    },
    gammaloop_integrand::{
        GLIntegrand,
        cross_section_integrand::{CrossSectionGraphTerm, CrossSectionIntegrand},
    },
    graph::{ExternalConnection, FeynmanGraph, Graph},
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
        vakint: &Vakint,
        add_lu_cut_feynman_rules: bool,
        settings: &GenerationSettings,
        conjugate: bool,
    ) -> ParsingNet {
        if let Some(network) = &self.network {
            network.clone()
        } else {
            let all_cut_edges = graph
                .iter_edges_of(&self.threshold_cut)
                .chain(graph.iter_edges_of(lu_cut))
                .map(|(_, e, _)| e)
                .collect_vec();

            let left_orientations = get_orientations_from_subgraph(
                &graph,
                &self.left_subgraph,
                &self.reversed_dangling_edges,
            )
            .into_iter()
            .map(|cff_graph| cff_graph.global_orientation)
            .filter(|or| settings.orientation_pattern.alt_filter(or))
            .collect::<TiVec<SuperGraphOrientationID, _>>();

            let right_orientations = get_orientations_from_subgraph(
                &graph,
                &self.right_subgraph,
                &self.reversed_dangling_edges,
            )
            .into_iter()
            .map(|cff_graph| cff_graph.global_orientation)
            .filter(|or| settings.orientation_pattern.alt_filter(or))
            .collect::<TiVec<SuperGraphOrientationID, _>>();

            let left_wood = graph.wood(&self.left_subgraph);
            let right_wood = graph.wood(&self.right_subgraph);

            let mut left_forest = left_wood.unfold(graph, &graph.loop_momentum_basis);
            let mut right_forest = right_wood.unfold(graph, &graph.loop_momentum_basis);

            left_forest.compute(
                graph,
                &self.left_subgraph,
                vakint,
                &left_orientations,
                &None,
                &all_cut_edges,
                &settings.uv,
                conjugate,
            );

            right_forest.compute(
                graph,
                &self.right_subgraph,
                vakint,
                &right_orientations,
                &None,
                &all_cut_edges,
                &settings.uv,
                conjugate,
            );

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
            );

            let cut_edges_for_right = if add_lu_cut_feynman_rules {
                Some(lu_cut.left.union(&lu_cut.right))
            } else {
                None
            };

            let right_expr =
                right_forest.orientation_parametric_expr(cut_edges_for_right.as_ref(), graph);

            self.network = Some(left_expr * right_expr);
            self.network.clone().unwrap()
        }
    }
}

define_index! {pub struct GlobalThresholdId;}
define_index! {pub struct RightThresholdId;}
define_index! {pub struct LeftThresholdId;}

use derive_more::{From, Into};
#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSection {
    pub name: String,
    pub integrand: Option<GLIntegrand>,
    pub supergraphs: Vec<CrossSectionGraph>,
    pub external_particles: Vec<ArcParticle>,
    pub external_connections: Vec<ExternalConnection>,
    pub n_incmoming: usize,
    pub graph_group_structure: TiVec<GroupId, GraphGroup>,
}

impl CrossSection {
    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
    ) -> Result<(), std::io::Error> {
        for graph in &self.supergraphs {
            graph.write_dot(writer)?;
            writeln!(writer)?;
        }
        Ok(())
    }

    pub fn write_dot_fmt<W: std::fmt::Write>(&self, writer: &mut W) -> Result<(), std::fmt::Error> {
        for graph in &self.supergraphs {
            graph.write_dot_fmt(writer)?;
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

    pub fn from_graph_list(name: String, mut graphs: Vec<Graph>, model: &Model) -> Result<Self> {
        let mut cross_section = CrossSection::new(name);
        cross_section.graph_group_structure = complete_group_parsing(&mut graphs)?;
        println!("group structure: {:?}", cross_section.graph_group_structure);

        for mut cross_section_graph in graphs {
            cross_section_graph.param_builder = ParamBuilder::new(&cross_section_graph, model);
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
            if external_particles.len() % 2 != 0 {
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

        self.integrand = Some(GLIntegrand::CrossSection(cross_section_integrand));
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
        let p = path.as_ref().join(format!("cs_{}", self.name));
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
            file.write(&binary)?;
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
            cs.integrand = Some(GLIntegrand::CrossSection(integrand));
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
        self.generate_cuts(model, process_definition)?;
        debug!("generating esurfaces corresponding to cuts");
        self.generate_esurface_cuts();
        debug!("generating cff");
        self.generate_cff(settings)?;
        debug!("extending cut esurface cache");
        self.update_surface_cache();
        debug!("building lmbs");
        self.build_lmbs();
        debug!("building multi channeling channels");

        if self.graph.is_group_master {
            self.build_multi_channeling_channels();
        }
        debug!("building parametric integrand");
        self.build_parametric_integrand(&settings)?;

        if settings.enable_thresholds {
            debug!("building threshold counterterm");
            self.build_threshold_counteterm(settings)?;
        }

        Ok(())
    }

    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
    ) -> Result<(), std::io::Error> {
        self.graph.dot_serialize_io(writer)
    }

    pub(crate) fn write_dot_fmt<W: std::fmt::Write>(
        &self,
        writer: &mut W,
    ) -> Result<(), std::fmt::Error> {
        self.graph.dot_serialize_fmt(writer)
    }

    pub(crate) fn update_surface_cache(&mut self) {
        let esurface_cache = &mut self
            .derived_data
            .cff_expression
            .as_mut()
            .unwrap()
            .surfaces
            .esurface_cache;

        // if a cut was not generated during cff, we still add it to the surface cache such that it has an esurface_id
        for esurface in self.cut_esurface.iter() {
            if let Some(esurface_id) = esurface_cache.iter().position(|e| e == esurface) {
                self.cut_esurface_id_map.push(esurface_id.into());
            } else {
                self.cut_esurface_id_map.push(esurface_cache.len().into());
                esurface_cache.push(esurface.clone());
            }
        }
    }

    fn generate_cff(&mut self, settings: &GenerationSettings) -> Result<()> {
        // hardcorde 1 to n for now
        debug!("generating cff");

        let shift_rewrite = self
            .graph
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let cff_cut_expression =
            generate_cff_with_cuts(&self.graph.underlying, &shift_rewrite, &self.cuts)?;

        Ok(self.derived_data.cff_expression = Some(cff_cut_expression))
    }

    fn generate_cuts(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
    ) -> Result<()> {
        info!("generatig cuts for graph: {}", self.graph.name);

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
    ) -> Result<()> {
        self.derived_data.cut_paramatric_integrand =
            self.build_original_parametric_integrand(settings)?;
        Ok(())
    }

    fn build_original_parametric_integrand(
        &mut self,
        settings: &GenerationSettings,
    ) -> Result<TiVec<CutId, Atom>> {
        let global_num = self.graph.global_network();

        let (externals, props) = self.graph.get_initial_state_tree();

        let canonize_esurface = self
            .graph
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let mut integrands = TiVec::new();

        for ((cut_id, cut), esurface) in self.cuts.iter_enumerated().zip(self.cut_esurface.iter()) {
            let expression_for_cut = &self
                .derived_data
                .cff_expression
                .as_ref()
                .unwrap()
                .cut_expressions[cut_id];

            let left_orientations = expression_for_cut
                .left_amplitude
                .orientations
                .iter()
                .map(|expr| expr.data.orientation.clone())
                .filter(|o| settings.orientation_pattern.alt_filter(o))
                .collect::<TiVec<AmplitudeOrientationID, _>>();

            let right_orientations = expression_for_cut
                .right_amplitude
                .orientations
                .iter()
                .map(|expr| expr.data.orientation.clone())
                .filter(|o| settings.orientation_pattern.alt_filter(o))
                .collect::<TiVec<AmplitudeOrientationID, _>>();

            let vakint = self.new_vakint();

            let left_wood = self.graph.wood(&cut.left);
            let right_wood = self.graph.wood(&cut.right);

            let mut left_forest = left_wood.unfold(&self.graph, &self.graph.loop_momentum_basis);
            let mut right_forest = right_wood.unfold(&self.graph, &self.graph.loop_momentum_basis);

            left_forest.compute(
                &self.graph,
                &cut.left,
                &vakint,
                &left_orientations,
                &canonize_esurface,
                &esurface.energies,
                &settings.uv,
                false,
            );

            right_forest.compute(
                &self.graph,
                &cut.right,
                &vakint,
                &right_orientations,
                &canonize_esurface,
                &esurface.energies,
                &settings.uv,
                true,
            );

            let left_expr = left_forest.orientation_parametric_expr(
                Some(&esurface.bitvec(&self.graph.underlying)),
                &self.graph,
            );

            let right_expr = right_forest.orientation_parametric_expr(None, &self.graph);

            if settings.enable_thresholds {
                // we can recycle these tensor networks for the single threshold counterterms
                self.derived_data
                    .tensor_network_cache
                    .push((left_expr.clone(), right_expr.clone()));
            }

            let mut product = left_expr * right_expr * global_num.clone();

            product
                .execute::<Sequential, SmallestDegree, _, _, _>(
                    TENSORLIB.read().unwrap().deref(),
                    FUN_LIB.deref(),
                )
                .unwrap();

            let scalar: Atom = product
                .result_scalar()
                .with_context(|| "in building LU integrand")?
                .into();

            let mut integrand = scalar
                .unwrap_function(GS.color_wrap)
                .simplify_color()
                .replace(function!(GS.energy, W_.x_))
                .with(function!(GS.ose, W_.x_))
                .replace(function!(GS.theta, W_.x_).pow(Atom::var(W_.n_)))
                .with(function!(GS.theta, W_.x_));

            for (_, edge_index, _) in self
                .graph
                .underlying
                .iter_edges_of(&self.graph.initial_state_cut)
            {
                integrand = integrand
                    .replace(GS.ose(edge_index))
                    .with(GS.emr_mom(edge_index, Atom::from(ExpandedIndex::from_iter([0]))));
            }

            let prefactor = self.lu_prefactor_helper();
            let integrand_with_prefactor = prefactor * integrand;
            debug!("integrand for cut {}: {}", cut_id, integrand_with_prefactor);
            integrands.push(integrand_with_prefactor);
        }

        Ok(integrands)
    }

    fn lu_prefactor_helper(&self) -> Atom {
        let loop_number = self.graph.cyclotomatic_number(&self.graph.full_filter())
            - self.graph.initial_state_cut.nedges(&self.graph);

        let loop_3 = loop_number as i64 * 3;
        let grad_eta = Atom::var(GS.deta_lu_cut);
        let factors_of_pi = (Atom::num(2) * Atom::var(GS.pi)).npow(loop_3 - 1); // multiply with 2pi from energy conservation delta

        let tstar = Atom::var(GS.rescale_star);
        let tsrat_pow = tstar.npow(loop_3);
        let hfunction = Atom::var(GS.hfunction_lu_cut);
        tsrat_pow * hfunction / factors_of_pi / grad_eta
    }

    fn th_prefactor_helper(&self, subspace_loop_count: usize, is_on_right: bool) -> Atom {
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

        let i = Atom::i();
        let pi = Atom::var(GS.pi);

        let jacobian_ratio = (&radius_star / &radius).npow(subspace_loop_count as i64 * 3 - 1);

        let local_prefactor = (&uv_damp_plus / (&radius - &radius_star)
            + &uv_damp_minus / (-&radius - &radius_star))
            / &grad_eta
            * &jacobian_ratio;

        let integrated_prefactor = &h_function * &i * &pi * &jacobian_ratio / &grad_eta;

        local_prefactor + integrated_prefactor
    }

    fn new_vakint(&self) -> Vakint {
        Vakint::new(Some(VakintSettings {
            allow_unknown_integrals: false,
            evaluation_order: EvaluationOrder::alphaloop_only(),
            integral_normalization_factor: LoopNormalizationFactor::MSbar,
            run_time_decimal_precision: 32,
            number_of_terms_in_epsilon_expansion: self.graph.n_loops(&self.graph.no_dummy()) as i64
                + 1,
            // temporary_directory: Some("./form".into()),
            mu_r_sq_symbol: GS.mu_r_sq.get_name().to_string(),
            ..VakintSettings::default()
        }))
        .unwrap()
    }

    fn build_lmbs(&mut self) {
        let lmbs = self.graph.generate_loop_momentum_bases(
            &self
                .graph
                .underlying
                .full_filter()
                .subtract(&self.graph.initial_state_cut.left)
                .subtract(&self.graph.initial_state_cut.right),
        );

        self.derived_data.lmbs = Some(lmbs)
    }

    fn build_multi_channeling_channels(&mut self) {
        let channels = self
            .graph
            .build_multi_channeling_channels(self.derived_data.lmbs.as_ref().unwrap());

        self.derived_data.multi_channeling_setup = Some(channels)
    }

    fn build_threshold_counteterm(&mut self, settings: &GenerationSettings) -> Result<()> {
        // thershold enumeration as st cuts
        let all_possible_thresholds: TiVec<GlobalThresholdId, _> = {
            let mut unsorted = self.graph.all_st_cuts_for_cs(
                self.source_nodes.clone(),
                self.target_nodes.clone(),
                &self.graph.get_initial_state_tree().0,
            );

            unsorted.sort_by(|a, b| a.1.cmp(&b.1));
            unsorted.into()
        };

        let vakint = self.new_vakint();
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

                        let mut threshold_esurface = Esurface::new_from_cut_left(
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
                } else if cut.right.includes(right_threshold_diagram) {
                    if self.graph.underlying.cyclotomatic_number(&cut.right)
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

                        let mut threshold_esurface = Esurface::new_from_cut_left(
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
            }

            let left_counterterms = thresholds_on_the_left
                .iter_mut()
                .map(|ct_diagram| {
                    let left_ct = ct_diagram.get_tensor_network_cached(
                        &self.graph,
                        &cut.cut,
                        &vakint,
                        true,
                        settings,
                        false,
                    );

                    let mut product = left_ct
                        * self.derived_data.tensor_network_cache[cut_id].1.clone()
                        * global_num.clone();

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
                        .with(function!(GS.theta, W_.x_));

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

                    let th_prefactor = self.th_prefactor_helper(left_loop_count, false);

                    Ok(left_integrand * lu_prefactor * th_prefactor)
                })
                .collect::<Result<TiVec<LeftThresholdId, Atom>>>()?;

            let right_counterterms = thresholds_on_the_right
                .iter_mut()
                .map(|ct_diagram| {
                    let right_ct = ct_diagram.get_tensor_network_cached(
                        &self.graph,
                        &cut.cut,
                        &vakint,
                        false,
                        settings,
                        true,
                    );

                    let mut product = self.derived_data.tensor_network_cache[cut_id].0.clone()
                        * right_ct
                        * global_num.clone();

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
                        .with(function!(GS.theta, W_.x_));

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

                    let th_prefactor = self.th_prefactor_helper(right_loop_count, true);

                    Ok(right_integrand * lu_prefactor * th_prefactor)
                })
                .collect::<Result<TiVec<RightThresholdId, Atom>>>()?;

            let iterated_counterterms = thresholds_on_the_left
                .iter()
                .cartesian_product(&thresholds_on_the_right)
                .map(|(left_ct_diagram, right_ct_diagram)| {
                    let mut product = left_ct_diagram.network.clone().unwrap()
                        * right_ct_diagram.network.clone().unwrap()
                        * global_num.clone();

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
                        .with(function!(GS.theta, W_.x_));

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
                        self.graph
                            .underlying
                            .cyclotomatic_number(&left_ct_diagram.left_subgraph),
                        false,
                    );
                    let right_prefactor = self.th_prefactor_helper(
                        self.graph
                            .underlying
                            .cyclotomatic_number(&right_ct_diagram.right_subgraph),
                        true,
                    );

                    Ok(iterated_scalar * lu_prefactor * left_prefactor * right_prefactor)
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

        self.derived_data.threshold_counterterms = counterterms;
        Ok(())
    }

    fn generate_term_for_graph(
        &self,
        _model: &Model,
        settings: &GlobalSettings,
    ) -> Result<CrossSectionGraphTerm> {
        CrossSectionGraphTerm::from_cross_section_graph(&self, settings)
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionDerivedData {
    pub orientations: Option<TiVec<SuperGraphOrientationID, CutOrientationData>>,
    pub cut_paramatric_integrand: TiVec<CutId, Atom>,
    pub cff_expression: Option<CFFCutsExpression>,
    pub lmbs: Option<TiVec<LmbIndex, LoopMomentumBasis>>,
    pub multi_channeling_setup: Option<LmbMultiChannelingSetup>,
    pub tensor_network_cache: TiVec<CutId, (ParsingNet, ParsingNet)>,
    pub threshold_counterterms: TiVec<CutId, LUCounterTermData>,
}

impl CrossSectionDerivedData {
    fn new_empty() -> Self {
        Self {
            orientations: None,
            cff_expression: None,
            cut_paramatric_integrand: TiVec::new(),
            lmbs: None,
            multi_channeling_setup: None,
            tensor_network_cache: TiVec::new(),
            threshold_counterterms: TiVec::new(),
        }
    }
}
