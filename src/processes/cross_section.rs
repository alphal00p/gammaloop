use std::{
    fmt::{Display, Formatter},
    fs::{self, File},
    io::Write,
    marker::PhantomData,
    ops::Deref,
    path::Path,
};

use ahash::AHashSet;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use bitvec::vec::BitVec;
use color_eyre::Result;

use idenso::{color::ColorSimplifier, metric::MS};
use rayon::ThreadPool;
use spenso::network::{Sequential, SmallestDegree};
use tracing::info;
use vakint::{EvaluationOrder, LoopNormalizationFactor, Vakint, VakintSettings};

use crate::{
    cff::{cut_expression::SuperGraphOrientationID, expression::AmplitudeOrientationID},
    gammaloop_integrand::{
        cross_section_integrand::{CrossSectionIntegrandData, OrientationEvaluator},
        param_builder::ParamBuilderGraph,
        GenericEvaluator, LmbMultiChannelingSetup, ParamBuilder,
    },
    graph::{
        get_cff_inverse_energy_product_impl, parse::complete_group_parsing, GraphGroup, GroupId,
        LMBext, LmbIndex, LoopMomentumBasis,
    },
    model::ArcParticle,
    numerator::symbolica_ext::AtomCoreExt,
    processes::Amplitude,
    settings::{global::GenerationSettings, runtime::LockedRuntimeSettings, GlobalSettings},
    status_info,
    utils::{ose_atom_from_index, GS, TENSORLIB, W_},
    uv::UltravioletGraph,
    GammaLoopContext, GammaLoopContextContainer,
};
use eyre::{eyre, Context};
use itertools::Itertools;
use linnet::half_edge::{
    involution::{Flow, HedgePair, Orientation},
    subgraph::{HedgeNode, Inclusion, InternalSubGraph, OrientedCut},
};
use log::{debug, warn};
use serde::{Deserialize, Serialize};
use symbolica::{
    atom::{Atom, AtomCore},
    domains::rational::Rational,
    evaluate::{FunctionMap, OptimizationSettings},
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
        cross_section_integrand::{CrossSectionGraphTerm, CrossSectionIntegrand},
        GLIntegrand,
    },
    graph::{ExternalConnection, FeynmanGraph, Graph},
    model::Model,
};

use crate::processes::ProcessDefinition;

pub trait CrossSectionState:
    Clone + std::fmt::Debug + bincode::Encode + for<'a> bincode::Decode<GammaLoopContextContainer<'a>>
{
}
impl CrossSectionState for () {}

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

    pub fn from_graph_list(name: String, mut graphs: Vec<Graph>) -> Result<Self> {
        let mut cross_section = CrossSection::new(name);
        cross_section.graph_group_structure = complete_group_parsing(&mut graphs)?;

        for cross_section_graph in graphs {
            cross_section.add_supergraph(cross_section_graph)?;
        }
        Ok(cross_section)
    }

    pub(crate) fn warm_up(&mut self, _model: &Model) -> Result<()> {
        let _derived_data = self
            .supergraphs
            .iter()
            .map(|g| &g.derived_data)
            .collect_vec();

        todo!()
        //let derived_data_container = DerivedDataContainer::CrossSection(&derived_data);

        //self.integrand
        //    .as_mut()
        //    .map(|a| a.warm_up(derived_data_container));
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
        generatioon_pool: &ThreadPool,
    ) -> Result<()> {
        for supergraph in &mut self.supergraphs {
            supergraph.preprocess(model, process_definition, generation_settings)?;
        }
        Ok(())
    }

    pub fn build_integrand(
        &mut self,
        model: &Model,
        global_settings: &GlobalSettings,
        runtime_default: LockedRuntimeSettings,
        generatioon_pool: &ThreadPool,
    ) -> Result<()> {
        let terms = self
            .supergraphs
            .iter()
            .map(|sg| sg.generate_term_for_graph(model))
            .collect_vec();

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
}

#[derive(Clone, bincode::Encode, bincode::Decode)]
pub struct CrossSectionCut {
    pub cut: OrientedCut,
    #[bincode(with_serde)]
    pub left: BitVec,
    #[bincode(with_serde)]
    pub right: BitVec,
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
        let mut source_nodes = AHashSet::new();
        let mut target_nodes = AHashSet::new();

        for (hedge_pair, _, _) in graph.underlying.iter_edges_of(&graph.initial_state_cut) {
            match hedge_pair {
                HedgePair::Split { source, sink, .. } => {
                    // yes it is supposed to be like this, becuase a sink hedge goes into the node from which to construct the cuts
                    let source_node = graph.underlying.node_id(source);
                    let sink_node = graph.underlying.node_id(sink);

                    source_nodes.insert(source_node);
                    target_nodes.insert(sink_node);
                }
                _ => {
                    unreachable!();
                }
            }
        }

        assert_eq!(source_nodes.len(), target_nodes.len());

        let source_node_vec = source_nodes.into_iter().collect_vec();
        let target_node_vec = target_nodes.into_iter().collect_vec();

        let source_node = graph
            .underlying
            .combine_to_single_hedgenode(&source_node_vec);

        let target_node = graph
            .underlying
            .combine_to_single_hedgenode(&target_node_vec);

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
        self.generate_cff()?;
        debug!("extending cut esurface cache");
        self.update_surface_cache();
        debug!("building lmbs");
        self.build_lmbs();
        debug!("building multi channeling channels");
        self.build_multi_channeling_channels();
        debug!("building parametric integrand");
        self.build_parametric_integrand(&settings)?;

        Ok(())
    }

    pub(crate) fn write_dot<W: std::io::Write>(
        &self,
        writer: &mut W,
    ) -> Result<(), std::io::Error> {
        self.graph.dot_serialize_io(writer)
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

    fn generate_cff(&mut self) -> Result<()> {
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

        let all_st_cuts = self
            .graph
            .underlying
            .all_cuts(self.source_nodes.clone(), self.target_nodes.clone())
            .into_iter()
            .map(|(l, mut c, r)| {
                // remove initial state cut edges from cut
                self.graph
                    .initial_state_cut
                    .left
                    .iter()
                    .by_vals()
                    .enumerate()
                    .chain(
                        self.graph
                            .initial_state_cut
                            .right
                            .iter()
                            .by_vals()
                            .enumerate(),
                    )
                    .for_each(|(index, is_set)| {
                        if is_set {
                            c.left.set(index, false);
                            c.right.set(index, false);
                        }
                    });

                (l, c, r)
            });

        info!("num s_t cuts: {}", all_st_cuts.len());

        let mut cuts: TiVec<CutId, CrossSectionCut> = all_st_cuts
            .into_iter()
            .map(|(left, cut, right)| CrossSectionCut { cut, left, right })
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
        &self,
        settings: &GenerationSettings,
    ) -> Result<TiVec<CutId, Atom>> {
        let global_num = self.graph.global_network();

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
                .filter(|o| settings.orientation_pattern.filter(o))
                .collect::<TiVec<AmplitudeOrientationID, _>>();

            let right_orientations = expression_for_cut
                .right_amplitude
                .orientations
                .iter()
                .map(|expr| expr.data.orientation.clone())
                .filter(|o| settings.orientation_pattern.filter(o))
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
            );

            right_forest.compute(
                &self.graph,
                &cut.right,
                &vakint,
                &right_orientations,
                &canonize_esurface,
                &esurface.energies,
                &settings.uv,
            );

            let left_expr = left_forest.orientation_parametric_expr(
                Some(&esurface.bitvec(&self.graph.underlying)),
                &self.graph,
            );

            let right_expr = right_forest.orientation_parametric_expr(None, &self.graph);

            let mut product = left_expr * right_expr * global_num.clone();

            product
                .execute::<Sequential, SmallestDegree, _, _>(TENSORLIB.read().unwrap().deref())
                .unwrap();

            let scalar: Atom = product
                .result_scalar()
                .with_context(|| "in building LU integrand")?
                .into();

            let integrand = scalar
                .unwrap_function(GS.color_wrap)
                .simplify_color()
                .replace(function!(GS.energy, W_.x_))
                .with(function!(GS.ose, W_.x_));

            let loop_number = self.graph.cyclotomatic_number(&cut.left)
                + self.graph.n_loops(&cut.right)
                + (esurface.energies.len() - 1);

            let loop_3 = loop_number as i64 * 3;
            let grad_eta = Atom::var(GS.deta);
            let factors_of_pi = (Atom::num(2) * Atom::var(GS.pi)).npow(loop_3);
            let tstar = Atom::var(GS.rescale_star);
            let tsrat_pow = tstar.npow(loop_3);
            let hfunction = Atom::var(GS.hfunction);

            let prefactor = tsrat_pow * hfunction / factors_of_pi / grad_eta;
            let integrand_with_prefactor = prefactor * integrand;
            debug!("integrand for cut {}: {}", cut_id, integrand_with_prefactor);
            integrands.push(integrand_with_prefactor);
        }

        Ok(integrands)
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
        let lmbs = self
            .graph
            .underlying
            .generate_loop_momentum_bases(&self.graph.underlying.full_filter());

        self.derived_data.lmbs = Some(lmbs)
    }

    fn build_multi_channeling_channels(&mut self) {
        let channels = self
            .graph
            .build_multi_channeling_channels(self.derived_data.lmbs.as_ref().unwrap());

        self.derived_data.multi_channeling_setup = Some(channels)
    }

    fn generate_term_for_graph(&self, _model: &Model) -> CrossSectionGraphTerm {
        todo!()
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
}

impl CrossSectionDerivedData {
    fn new_empty() -> Self {
        Self {
            orientations: None,
            cff_expression: None,
            cut_paramatric_integrand: TiVec::new(),
            lmbs: None,
            multi_channeling_setup: None,
        }
    }
}
