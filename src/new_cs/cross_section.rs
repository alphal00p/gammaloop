use std::{
    fmt::{Display, Formatter},
    marker::PhantomData,
};

use ahash::AHashSet;
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use bitvec::vec::BitVec;
use color_eyre::Result;

use idenso::metric::MS;

use crate::{
    cff::{
        cut_expression::SuperGraphOrientationID,
        esurface::{generate_esurface_data, EsurfaceDerivedData},
    },
    model::ArcParticle,
    new_gammaloop_integrand::{
        cross_section_integrand::{CrossSectionIntegrandData, OrientationEvaluator},
        GenericEvaluator, LmbMultiChannelingSetup, ParamBuilder,
    },
    new_graph::{get_cff_inverse_energy_product_impl, LMBext, LmbIndex, LoopMomentumBasis},
    utils::{ose_atom_from_index, GS, W_},
    GammaLoopContext, GammaLoopContextContainer,
};
use eyre::eyre;
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
    cross_section::IsPolarizable,
    integrands::Integrand,
    model::Model,
    momentum::{Rotatable, Rotation, RotationMethod},
    new_gammaloop_integrand::{
        cross_section_integrand::{CrossSectionGraphTerm, CrossSectionIntegrand},
        NewIntegrand,
    },
    new_graph::{ExternalConnection, FeynmanGraph, Graph},
    DependentMomentaConstructor, Externals, Polarizations, Settings,
};

use super::ProcessDefinition;

pub trait CrossSectionState:
    Clone + std::fmt::Debug + bincode::Encode + for<'a> bincode::Decode<GammaLoopContextContainer<'a>>
{
}
impl CrossSectionState for () {}

use derive_more::{From, Into};
#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSection<S: CrossSectionState = ()> {
    pub name: String,
    pub integrand: Option<NewIntegrand>,
    pub supergraphs: Vec<CrossSectionGraph<S>>,
    pub external_particles: Vec<ArcParticle>,
    pub external_connections: Vec<ExternalConnection>,
    pub n_incmoming: usize,
}

impl<S: CrossSectionState> CrossSection<S> {
    pub(crate) fn new(name: String) -> Self {
        Self {
            name,
            integrand: None,
            supergraphs: vec![],
            external_connections: vec![],
            external_particles: vec![],
            n_incmoming: 0,
        }
    }

    pub(crate) fn add_supergraph(&mut self, supergraph: Graph) -> Result<()> {
        if self.external_particles.is_empty() {
            let external_particles = supergraph.underlying.get_external_partcles();
            if external_particles.len() % 2 != 0 {
                return Err(eyre!(
                    "expected even number of externals for forward scattering graph"
                ));
            }
            self.external_particles = external_particles;
            self.n_incmoming = self.external_particles.len() / 2;
        } else if self.external_particles != supergraph.underlying.get_external_partcles() {
            return Err(eyre!(
                "attempt to add supergraph with differnt external particles"
            ));
        }

        let cross_section_graph = CrossSectionGraph::new(supergraph);
        self.supergraphs.push(cross_section_graph);

        // TODO: validate that the graph is compatible
        Ok(())
    }

    pub(crate) fn preprocess(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
    ) -> Result<()> {
        for supergraph in &mut self.supergraphs {
            supergraph.preprocess(model, process_definition)?;
        }
        Ok(())
    }

    pub(crate) fn build_integrand(&mut self, settings: Settings, model: &Model) -> Result<()> {
        let terms = self
            .supergraphs
            .iter()
            .map(|sg| sg.generate_term_for_graph(&settings))
            .collect_vec();

        let rotations: Vec<Rotation> = Some(Rotation::new(RotationMethod::Identity))
            .into_iter()
            .chain(
                settings
                    .stability
                    .rotation_axis
                    .iter()
                    .map(|axis| Rotation::new(axis.rotation_method())),
            )
            .collect(); // want this to include the identity rotation (i.e the first sample)

        let orig_polarizations = self.polarizations(&settings.kinematics.externals);

        let polarizations = rotations
            .iter()
            .map(|r| orig_polarizations.rotate(r))
            .collect();

        // let model_parameter_cache = model.generate_values();

        let cross_section_integrand = CrossSectionIntegrand {
            settings,
            data: CrossSectionIntegrandData {
                rotations,
                name: self.name.clone(),
                external_connections: self.external_connections.clone(),
                n_incoming: self.n_incmoming,
                polarizations,
                graph_terms: terms,
            },
        };

        self.integrand = Some(NewIntegrand::CrossSection(cross_section_integrand));
        Ok(())
    }
}

impl<S: CrossSectionState> IsPolarizable for CrossSection<S> {
    fn polarizations(&self, externals: &Externals) -> Polarizations {
        externals.generate_polarizations(
            &self.external_particles,
            DependentMomentaConstructor::CrossSection {
                external_connections: &self.external_connections,
            },
        )
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
    pub(crate) fn is_s_channel<S: CrossSectionState>(
        &self,
        cross_section_graph: &CrossSectionGraph<S>,
    ) -> Result<bool> {
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

    pub(crate) fn is_valid_for_process<S: CrossSectionState>(
        &self,
        cross_section_graph: &CrossSectionGraph<S>,
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

                    if cut_content.len() > process.n_unresolved {
                        debug!(" too many unresolved particles");
                        return false;
                    }

                    if !cut_content
                        .iter()
                        .all(|particle| process.unresolved_cut_content.contains(particle))
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
pub struct CrossSectionGraph<S: CrossSectionState = ()> {
    pub graph: Graph,
    pub source_nodes: HedgeNode,
    pub target_nodes: HedgeNode,
    pub cuts: TiVec<CutId, CrossSectionCut>,
    pub cut_esurface: TiVec<CutId, Esurface>,
    pub cut_esurface_id_map: TiVec<CutId, EsurfaceID>,
    pub derived_data: CrossSectionDerivedData<S>,
}

impl<S: CrossSectionState> CrossSectionGraph<S> {
    pub(crate) fn new(graph: Graph) -> Self {
        let mut source_nodes = AHashSet::new();
        let mut target_nodes = AHashSet::new();

        for (hedge_pair, _, _) in graph.underlying.iter_edges() {
            match hedge_pair {
                HedgePair::Unpaired { hedge, flow } => {
                    let node_id = graph.underlying.node_id(hedge);
                    match flow {
                        Flow::Source => {
                            target_nodes.insert(node_id);
                        }
                        Flow::Sink => {
                            source_nodes.insert(node_id);
                        }
                    }
                }
                _ => continue,
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
            derived_data: CrossSectionDerivedData::<S>::new_empty(),
        }
    }

    pub(crate) fn preprocess(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
    ) -> Result<()> {
        self.generate_cuts(model, process_definition)?;
        self.generate_esurface_cuts();
        self.generate_cff()?;
        self.update_surface_cache();

        //self.build_cut_evaluators(model, None);
        //self.build_orientation_evaluators(model);
        self.build_lmbs();
        self.build_esurface_derived_data()?;
        Ok(self.build_multi_channeling_channels())
    }

    pub(crate) fn build_esurface_derived_data(&mut self) -> Result<()> {
        let lmbs = self.derived_data.lmbs.as_ref().unwrap();
        let esurfaces = &self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .esurface_cache;

        let esurface_data = generate_esurface_data(&self.graph, lmbs, esurfaces)?;
        Ok(self.derived_data.esurface_data = Some(esurface_data))
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
            .underlying
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
        debug!("generatig cuts for graph: {}", self.graph.name);

        let all_st_cuts = self
            .graph
            .underlying
            .all_cuts(self.source_nodes.clone(), self.target_nodes.clone());

        debug!("num s_t cuts: {}", all_st_cuts.len());

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

        debug!(
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
            .map(|cut| Esurface::new_from_cut_left(&self.graph.underlying, cut))
            .collect();

        self.cut_esurface = esurfaces;
    }

    fn build_left_right_amplitudes(&self, cut: CutId) -> (Atom, Atom) {
        let (left_amplitude, right_amplitude) = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .to_atom_for_cut(cut);

        let cut_edges = self
            .graph
            .underlying
            .iter_edges_of(&self.cuts[cut].cut)
            .map(|(_, id, _)| id)
            .collect_vec();

        let left_amplitude_energy_sub = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .substitute_energies(&left_amplitude, &cut_edges);

        let right_amplitude_energy_sub = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .substitute_energies(&right_amplitude, &cut_edges);

        let left_ose_product =
            get_cff_inverse_energy_product_impl(&self.graph.underlying, &self.cuts[cut].left, &[]);
        let right_ose_product =
            get_cff_inverse_energy_product_impl(&self.graph.underlying, &self.cuts[cut].right, &[]);

        (
            left_amplitude_energy_sub * left_ose_product,
            right_amplitude_energy_sub * right_ose_product,
        )
    }

    pub(crate) fn build_left_right_amplitudes_for_orientation(
        &self,
        cut: CutId,
        orientation: SuperGraphOrientationID,
    ) -> (Atom, Atom) {
        let (left_amplitude, right_amplitude) = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .get_atom_for_orientation_and_cut(orientation, cut);

        let cut_edges = self
            .graph
            .underlying
            .iter_edges_of(&self.cuts[cut].cut)
            .map(|(_, id, _)| id)
            .collect_vec();

        let left_amplitude_energy_sub = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .substitute_energies(&left_amplitude, &cut_edges);

        let right_amplitude_energy_sub = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .substitute_energies(&right_amplitude, &cut_edges);

        let left_ose_product =
            get_cff_inverse_energy_product_impl(&self.graph.underlying, &self.cuts[cut].left, &[]);
        let right_ose_product =
            get_cff_inverse_energy_product_impl(&self.graph.underlying, &self.cuts[cut].right, &[]);

        (
            left_amplitude_energy_sub * left_ose_product,
            right_amplitude_energy_sub * right_ose_product,
        )
    }

    fn build_atom_for_cut(&self, cut_id: CutId) -> Atom {
        let (left_amplitude, right_amplitude) = self.build_left_right_amplitudes(cut_id);

        let product = left_amplitude * right_amplitude;

        let ose_atom = self.add_additional_factors_to_cff_atom(&product, cut_id);
        debug!("ose atom for cut: {}: {}", cut_id, ose_atom);

        let replacements = self.graph.underlying.get_ose_replacements();
        let replaced_atom = ose_atom.replace_multiple(&replacements);
        let replace_dots = replaced_atom
            .replace(function!(
                MS.dot,
                function!(GS.emr_vec, W_.x_),
                function!(GS.emr_vec, W_.y_)
            ))
            .with(
                -(function!(GS.emr_vec, W_.x_, 1) * function!(GS.emr_vec, W_.y_, 1)
                    + function!(GS.emr_vec, W_.x_, 2) * function!(GS.emr_vec, W_.y_, 2)
                    + function!(GS.emr_vec, W_.x_, 3) * function!(GS.emr_vec, W_.y_, 3)),
            );

        // debug!("replaced atom: {}", replace_dots);
        replace_dots
    }

    fn build_atom_for_orientation_and_cut(
        &self,
        cut_id: CutId,
        orientation: SuperGraphOrientationID,
    ) -> Atom {
        let (left_amplitude, right_amplitude) =
            self.build_left_right_amplitudes_for_orientation(cut_id, orientation);

        let product = left_amplitude * right_amplitude;

        let ose_atom = self.add_additional_factors_to_cff_atom(&product, cut_id);
        let replacements = self.graph.underlying.get_ose_replacements();
        let replaced_atom = ose_atom.replace_multiple(&replacements);
        let replace_dots = replaced_atom
            .replace(function!(
                MS.dot,
                function!(GS.emr_vec, W_.x_),
                function!(GS.emr_vec, W_.y_)
            ))
            .with(
                -(function!(GS.emr_vec, W_.x_, 1) * function!(GS.emr_vec, W_.y_, 1)
                    + function!(GS.emr_vec, W_.x_, 2) * function!(GS.emr_vec, W_.y_, 2)
                    + function!(GS.emr_vec, W_.x_, 3) * function!(GS.emr_vec, W_.y_, 3)),
            );

        debug!("replaced atom: {}", replace_dots);
        replace_dots
    }

    // do not use this function together with do_replacement_rules
    pub(crate) fn add_additional_factors_to_cff_atom(
        &self,
        cut_atom: &Atom,
        cut_id: CutId,
    ) -> Atom {
        let loop_3 = self.graph.underlying.get_loop_number() as i64 * 3;
        let t_star_factor = Atom::var(GS.rescale_star).npow(loop_3);

        let h_function = Atom::var(GS.hfunction);
        let grad_eta = Atom::var(GS.deta);

        let factors_of_pi = (Atom::num(2) * Atom::var(GS.pi)).npow(loop_3);

        let cut_inverse_energy_product = Atom::num(1)
            / self
                .graph
                .underlying
                .iter_edges_of(&self.cuts[cut_id].cut)
                .map(|(_, edge_id, _)| Atom::num(2) * ose_atom_from_index(edge_id))
                .fold(Atom::num(1), |product, factor| product * factor);

        let result = cut_atom * cut_inverse_energy_product * t_star_factor * h_function
            / grad_eta
            / factors_of_pi;

        //debug!("result: {}", result);
        result
    }

    fn get_builder(&self, model: &Model) -> ParamBuilder<f64> {
        let mut params = vec![];

        // all external energies
        params.extend(self.graph.underlying.get_external_energy_atoms());

        // spatial components of external momenta
        for (pair, edge_id, _) in self.graph.underlying.iter_edges() {
            match pair {
                HedgePair::Unpaired { .. } => {
                    let i64_id = Into::<usize>::into(edge_id) as i64;
                    let external_spatial = [
                        function!(GS.external_mom, i64_id, 1),
                        function!(GS.external_mom, i64_id, 2),
                        function!(GS.external_mom, i64_id, 3),
                    ];
                    params.extend(external_spatial);
                }
                _ => {}
            }
        }

        // spatial EMR
        for (pair, edge_id, _) in self.graph.underlying.iter_edges() {
            match pair {
                HedgePair::Paired { .. } => {
                    let i64_id = Into::<usize>::into(edge_id) as i64;
                    let emr_components = [
                        function!(GS.emr_vec, i64_id, 1),
                        function!(GS.emr_vec, i64_id, 2),
                        function!(GS.emr_vec, i64_id, 3),
                    ];
                    params.extend(emr_components)
                }
                _ => {}
            }
        }

        // add model parameters
        params.extend(model.generate_params());
        // add additional parameters
        params.push(Atom::var(GS.m_uv));
        params.push(Atom::var(GS.mu_r_sq));
        params.push(Atom::var(GS.rescale_star));
        params.push(Atom::var(GS.hfunction));
        params.push(Atom::var(GS.deta));
        ParamBuilder::new()
    }

    fn get_function_map(&self) -> FunctionMap {
        let mut fn_map = FunctionMap::new();
        let pi_rational = Rational::from(std::f64::consts::PI);
        fn_map.add_constant(Atom::PI.into(), pi_rational.into());
        fn_map
    }

    pub(crate) fn build_cut_evaluators(
        &mut self,
        model: &Model,
        overwrite_atoms_for_test: Option<TiVec<CutId, Atom>>,
    ) {
        let evaluators = self
            .cuts
            .iter_enumerated()
            .map(|(cut_id, _)| {
                let atom = if let Some(atoms) = &overwrite_atoms_for_test {
                    atoms[cut_id].clone()
                } else {
                    self.build_atom_for_cut(cut_id)
                };

                let params = self.get_builder(model);

                let mut eval = GenericEvaluator::new_from_builder(
                    atom,
                    params,
                    OptimizationSettings::default(),
                );
                let filename = format!("{}_cut_{}.cpp", self.graph.name, cut_id);
                let function_name = format!("{}_cut_{}", self.graph.name, cut_id);
                let lib_name = format!("{}_cut_{}.so", self.graph.name, cut_id);

                eval.compile(
                    filename,
                    function_name,
                    lib_name,
                    symbolica::evaluate::InlineASM::None,
                );
                eval
            })
            .collect();

        self.derived_data.bare_cff_evaluators = Some(evaluators)
    }

    fn build_orientation_evaluators(&mut self, model: &Model) {
        let orientation_data = &self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientation_data;

        let substituted_energies = orientation_data
            .iter_enumerated()
            .map(|(orientation_id, data)| {
                data.cuts
                    .iter()
                    .map(|cut_id| self.build_atom_for_orientation_and_cut(*cut_id, orientation_id))
                    .collect_vec()
            })
            .collect::<TiVec<SuperGraphOrientationID, _>>();

        let orientation_evaluators = substituted_energies
            .iter()
            .zip(orientation_data)
            .map(|(cut_atoms, orientation_data)| {
                let cut_evaluators = cut_atoms
                    .iter()
                    .map(|cut_atom| {
                        let params = self.get_builder(model);

                        GenericEvaluator::new_from_builder(
                            cut_atom,
                            params,
                            OptimizationSettings::default(),
                        )
                    })
                    .collect_vec();

                OrientationEvaluator {
                    orientation_data: orientation_data.clone(),
                    evaluators: cut_evaluators,
                }
            })
            .collect();

        self.derived_data.bare_cff_orientation_evaluators = Some(orientation_evaluators);
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

    fn generate_term_for_graph(&self, settings: &Settings) -> CrossSectionGraphTerm {
        let estimated_scale = self
            .graph
            .underlying
            .expected_scale(settings.kinematics.e_cm);

        CrossSectionGraphTerm {
            multi_channeling_setup: self.derived_data.multi_channeling_setup.clone().unwrap(),
            bare_cff_evaluators: self.derived_data.bare_cff_evaluators.clone().unwrap(),
            bare_cff_orientation_evaluators: self
                .derived_data
                .bare_cff_orientation_evaluators
                .clone()
                .unwrap_or_else(|| vec![].into()),
            graph: self.graph.clone(),
            cuts: self.cuts.clone(),
            cut_esurface: self.cut_esurface.clone(),
            lmbs: self.derived_data.lmbs.clone().unwrap(),
            estimated_scale,
            param_builder: ParamBuilder::new(),
        }
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct CrossSectionDerivedData<S: CrossSectionState = ()> {
    pub orientations: Option<TiVec<SuperGraphOrientationID, CutOrientationData>>,
    pub bare_cff_evaluators: Option<TiVec<CutId, GenericEvaluator>>,
    pub bare_cff_orientation_evaluators:
        Option<TiVec<SuperGraphOrientationID, OrientationEvaluator>>,
    pub cff_expression: Option<CFFCutsExpression>,
    pub lmbs: Option<TiVec<LmbIndex, LoopMomentumBasis>>,
    pub multi_channeling_setup: Option<LmbMultiChannelingSetup>,
    pub esurface_data: Option<EsurfaceDerivedData>,
    pub _temp_numerator: Option<PhantomData<S>>,
}

impl<S: CrossSectionState> CrossSectionDerivedData<S> {
    fn new_empty() -> Self {
        Self {
            orientations: None,
            cff_expression: None,
            _temp_numerator: None,
            bare_cff_evaluators: None,
            bare_cff_orientation_evaluators: None,
            lmbs: None,
            multi_channeling_setup: None,
            esurface_data: None,
        }
    }
}
