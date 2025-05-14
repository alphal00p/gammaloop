use std::{
    cell::RefCell,
    collections::HashSet,
    fmt::{format, Display, Formatter},
    iter,
    marker::PhantomData,
    path::{Path, PathBuf},
    sync::Arc,
};

use ahash::{AHashSet, HashMap};
// use bincode::{Decode, Encode};
use bincode_trait_derive::{Decode, Encode};
use bitvec::vec::BitVec;
use color_eyre::Result;
use momtrop::SampleGenerator;
use smartstring::{LazyCompact, SmartString};
use spenso::contraction::IsZero;

use crate::{
    cff::{expression::CFFExpression, generation::generate_cff_expression},
    model::ArcParticle,
    momentum_sample::ExternalIndex,
    new_gammaloop_integrand::{
        amplitude_integrand::{AmplitudeGraphTerm, AmplitudeIntegrand},
        cross_section_integrand::OrientationEvaluator,
        GenericEvaluator, LmbMultiChannelingSetup,
    },
    new_graph::{LmbIndex, LoopMomentumBasis},
    signature::SignatureLike,
    utils::f128,
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
    evaluate::FunctionMap,
    parse,
};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cut_expression::{CFFCutExpression, CutOrientationData, OrientationID},
        esurface::{Esurface, EsurfaceID},
        generation::generate_cff_with_cuts,
    },
    cross_section::IsPolarizable,
    feyngen::{FeynGenFilters, GenerationType},
    graph::BareGraph,
    integrands::Integrand,
    model::Model,
    momentum::{Rotatable, Rotation, RotationMethod},
    new_gammaloop_integrand::{
        cross_section_integrand::{CrossSectionGraphTerm, CrossSectionIntegrand},
        NewIntegrand,
    },
    new_graph::{ExternalConnection, FeynmanGraph, Graph},
    numerator::{NumeratorState, PythonState},
    utils::F,
    DependentMomentaConstructor, Externals, Polarizations, ProcessSettings, Settings,
};

use derive_more::{From, Into};

#[derive(Debug, Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct ProcessDefinition {
    pub initial_pdgs: Vec<i64>, // Do we want a pub type Pdg = i64;?
    pub final_pdgs_lists: Vec<Vec<i64>>,
    pub n_unresolved: usize, // we need al this information to know what cuts are considered at runtime
    pub unresolved_cut_content: HashSet<ArcParticle>,
    pub amplitude_filters: FeynGenFilters,
    pub cross_section_filters: FeynGenFilters,
}

impl ProcessDefinition {
    pub fn new_empty() -> Self {
        Self {
            initial_pdgs: vec![],
            final_pdgs_lists: vec![],
            n_unresolved: 0,
            unresolved_cut_content: HashSet::new(),
            amplitude_filters: FeynGenFilters(vec![]),
            cross_section_filters: FeynGenFilters(vec![]),
        }
    }
}

#[derive(Clone)]
pub struct Process<S: NumeratorState = PythonState> {
    pub definition: ProcessDefinition,
    pub collection: ProcessCollection<S>,
}

impl<S: NumeratorState> Process<S> {
    pub fn preprocess(&mut self, model: &Model, settings: &ProcessSettings) -> Result<()> {
        self.collection
            .preprocess(model, &self.definition, settings)?;
        Ok(())
    }
}

impl Process {
    pub fn from_bare_graph_list(
        name: String,
        bare_graphs: Vec<BareGraph>,
        generation_type: GenerationType,
        definition: ProcessDefinition,
        sub_classes: Option<Vec<Vec<String>>>,
    ) -> Result<Self> {
        let graphs = bare_graphs.into_iter().map(Graph::from).collect_vec();

        match generation_type {
            GenerationType::Amplitude => {
                let mut collection: ProcessCollection<PythonState> =
                    ProcessCollection::new_amplitude();

                if let Some(_sub_classes) = sub_classes {
                    todo!("implement seperation of processes into user defined sub classes");
                } else {
                    let mut amplitude: Amplitude<PythonState> = Amplitude::new(name);

                    for amplitude_graph in graphs {
                        amplitude.add_graph(amplitude_graph)?;
                    }

                    collection.add_amplitude(amplitude);
                    Ok(Self {
                        definition,
                        collection,
                    })
                }
            }
            GenerationType::CrossSection => {
                let mut collection: ProcessCollection<PythonState> =
                    ProcessCollection::new_cross_section();

                if let Some(_sub_classes) = sub_classes {
                    todo!("implement seperation of processes into user defined sub classes");
                } else {
                    let mut cross_section: CrossSection<PythonState> = CrossSection::new(name);

                    for cross_section_graph in graphs {
                        cross_section.add_supergraph(cross_section_graph)?;
                    }

                    collection.add_cross_section(cross_section);
                    Ok(Self {
                        definition,
                        collection,
                    })
                }
            }
        }
    }

    fn generate_integrands(&self, settings: Settings) -> HashMap<String, Integrand> {
        self.collection.generate_integrands(settings)
    }
}

#[derive(Clone)]
pub struct ProcessList<S: NumeratorState = PythonState> {
    pub processes: Vec<Process<S>>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ExportSettings {
    pub root_folder: PathBuf,
}

impl Default for ProcessList {
    fn default() -> Self {
        Self::new()
    }
}

// Python Api
impl ProcessList {
    /// Generates a new empty process list
    pub fn new() -> Self {
        ProcessList { processes: vec![] }
    }

    pub fn add_process(&mut self, process: Process) {
        self.processes.push(process);
    }

    /// imports a process list from a folder
    pub fn import(_settings: ExportSettings) -> Result<Self> {
        Ok(Self::new())
    }

    ///preprocesses the process list according to the settings
    pub fn preprocess(&mut self, model: &Model, settings: ProcessSettings) -> Result<()> {
        for process in self.processes.iter_mut() {
            process.preprocess(model, &settings)?;
        }

        Ok(())
    }

    pub fn generate_integrands(&self, settings: Settings) -> HashMap<String, Integrand> {
        let mut result = HashMap::default();

        for process in self.processes.iter() {
            let integrands = process.generate_integrands(settings.clone());
            result.extend(integrands);
        }

        result
    }

    /// exports a process list to a folder
    pub fn export_amplitudes(
        &self,
        amplitude_names: &[String],
        settings: &ExportSettings,
    ) -> Result<usize> {
        let mut n_exported = 0;
        for process in self.processes.iter() {
            let process_definition_data =
                bincode::encode_to_vec(&process.definition, bincode::config::standard())?;

            let path = settings
                .root_folder
                .join("sources")
                .join("process_definition.bin");

            std::fs::write(path, &process_definition_data)?;

            n_exported += process
                .collection
                .export_amplitudes(amplitude_names, settings)?;
        }
        Ok(n_exported)
    }

    pub fn export_cross_sections(&self, settings: &ExportSettings) -> Result<()> {
        for process in self.processes.iter() {
            process.collection.export_cross_sections(settings)?;
        }
        Ok(())
    }
}

#[derive(Clone)]
pub enum ProcessCollection<S: NumeratorState = PythonState> {
    Amplitudes(Vec<Amplitude<S>>),
    CrossSections(Vec<CrossSection<S>>),
}

impl<S: NumeratorState> ProcessCollection<S> {
    fn new_amplitude() -> Self {
        Self::Amplitudes(vec![])
    }

    fn new_cross_section() -> Self {
        Self::CrossSections(vec![])
    }

    fn add_amplitude(&mut self, amplitude: Amplitude<S>) {
        match self {
            Self::Amplitudes(amplitudes) => amplitudes.push(amplitude),
            _ => panic!("Cannot add amplitude to a cross section collection"),
        }
    }

    fn add_cross_section(&mut self, cross_section: CrossSection<S>) {
        match self {
            Self::CrossSections(cross_sections) => cross_sections.push(cross_section),
            _ => panic!("Cannot add cross section to an amplitude collection"),
        }
    }

    fn preprocess(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
        settings: &ProcessSettings,
    ) -> Result<()> {
        match self {
            Self::Amplitudes(amplitudes) => {
                for amplitude in amplitudes {
                    amplitude.preprocess(model, settings)?;
                }
            }
            Self::CrossSections(cross_sections) => {
                for cross_section in cross_sections {
                    cross_section.preprocess(model, process_definition)?;
                }
            }
        }
        Ok(())
    }

    fn generate_integrands(&self, settings: Settings) -> HashMap<String, Integrand> {
        let mut result = HashMap::default();
        match self {
            Self::Amplitudes(amplitudes) => {
                let name = "default".to_owned();
                for amplitude in amplitudes {
                    let integrand = amplitude.generate_integrand(settings.clone());
                    result.insert(name.clone(), integrand);
                }
            }
            Self::CrossSections(cross_sections) => {
                let name = "default".to_owned();
                for cross_section in cross_sections {
                    let integrand = cross_section.generate_integrand(settings.clone());
                    result.insert(name.clone(), integrand);
                }
            }
        }
        result
    }

    fn export_amplitudes(
        &self,
        amplitude_names: &[String],
        settings: &ExportSettings,
    ) -> Result<usize> {
        let mut n_exported = 0;
        match self {
            Self::Amplitudes(amplitudes) => {
                for amplitude in amplitudes {
                    if amplitude_names.contains(&amplitude.name.to_string()) {
                        amplitude.export(settings)?;
                        n_exported += 1;
                    }
                }
            }
            _ => {}
        }
        Ok(n_exported)
    }

    fn export_cross_sections(&self, settings: &ExportSettings) -> Result<()> {
        match self {
            Self::CrossSections(cross_sections) => {
                for cross_section in cross_sections {
                    cross_section.export(settings)?;
                }
            }
            _ => {}
        }
        Ok(())
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait = GammaLoopContext)]
pub struct Amplitude<S: NumeratorState = PythonState> {
    name: String,
    graphs: Vec<AmplitudeGraph<S>>,
    external_particles: Vec<ArcParticle>,
    external_signature: SignatureLike<ExternalIndex>,
}

impl<S: NumeratorState> Amplitude<S> {
    pub fn preprocess(&mut self, model: &Model, settings: &ProcessSettings) -> Result<()> {
        for amplitude_graph in self.graphs.iter_mut() {
            amplitude_graph.preprocess(model, settings)?;
        }
        Ok(())
    }

    pub fn generate_integrand(&self, settings: Settings) -> Integrand {
        let terms = self
            .graphs
            .iter()
            .map(|graph| graph.generate_term_for_graph(&settings))
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
            .collect();

        let orig_polarizations = self.polarizations(&settings.kinematics.externals);

        let polarizations = rotations
            .iter()
            .map(|r| orig_polarizations.rotate(r))
            .collect();

        let amplitude_integrand = AmplitudeIntegrand {
            settings,
            rotations,
            polarizations,
            graph_terms: terms,
            external_signature: self.external_signature.clone(),
        };

        Integrand::NewIntegrand(NewIntegrand::Amplitude(amplitude_integrand))
    }

    pub fn export(&self, settings: &ExportSettings) -> Result<()> {
        let path = Path::new(&settings.root_folder)
            .join("sources")
            .join("amplitudes")
            .join(self.name.as_str());

        let data = bincode::encode_to_vec(self, bincode::config::standard())?;
        std::fs::write(
            path.clone().join(&format!("amplitude_{}.bin", self.name)),
            &data,
        )?;
        Ok(())
    }

    fn load_from_file(
        name: &str,
        root_folder: &str,
        context: GammaLoopContextContainer,
    ) -> Result<Self> {
        let path = Path::new(root_folder)
            .join("sources")
            .join("amplitudes")
            .join(name);

        let data = std::fs::read(path.join(format!("amplitude_{}.bin", name)))?;
        let (amplitude, _): (Self, _) =
            bincode::decode_from_slice_with_context(&data, bincode::config::standard(), context)?;
        Ok(amplitude)
    }
}

impl<S: NumeratorState> IsPolarizable for Amplitude<S> {
    fn polarizations(&self, externals: &Externals) -> Polarizations {
        externals.generate_polarizations(
            &self.external_particles,
            DependentMomentaConstructor::Amplitude(&self.external_signature),
        )
    }
}

#[derive(Clone, Encode, Decode)]
#[trait_decode(trait= GammaLoopContext)]
pub struct AmplitudeGraph<S: NumeratorState> {
    graph: Graph,
    derived_data: AmplitudeDerivedData<S>,
}

impl<S: NumeratorState> AmplitudeGraph<S> {
    fn new(graph: Graph) -> Self {
        AmplitudeGraph {
            graph,
            derived_data: AmplitudeDerivedData {
                cff_expression: None,
                bare_cff_evaluator: None,
                bare_cff_orientation_evaluatos: None,
                _temp_numerator: None,
                lmbs: None,
                tropical_sampler: None,
                multi_channeling_setup: None,
            },
        }
    }

    fn generate_cff(&mut self) -> Result<()> {
        let shift_rewrite = self
            .graph
            .underlying
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let cff_expression = generate_cff_expression(&self.graph.underlying, &shift_rewrite)?;
        self.derived_data.cff_expression = Some(cff_expression);

        Ok(())
    }

    fn preprocess(&mut self, _model: &Model, settings: &ProcessSettings) -> Result<()> {
        self.graph
            .loop_momentum_basis
            .set_edge_signatures(&self.graph.underlying)?;

        self.generate_cff()?;
        self.build_evaluator();
        self.build_evaluator_for_orientations()?;
        self.build_tropical_sampler(settings)?;
        self.build_loop_momentum_bases();
        self.build_multi_channeling_channels();

        Ok(())
    }

    fn build_multi_channeling_channels(&mut self) {
        let channels = self
            .graph
            .build_multi_channeling_channels(self.derived_data.lmbs.as_ref().unwrap());

        self.derived_data.multi_channeling_setup = Some(channels)
    }

    fn get_params(&self) -> Vec<Atom> {
        self.graph.underlying.get_energy_atoms()
    }

    fn get_function_map(&self) -> FunctionMap {
        let mut fn_map = FunctionMap::new();
        let pi_rational = Rational::from(std::f64::consts::PI);
        fn_map.add_constant(parse!("π").unwrap(), pi_rational);
        fn_map
    }

    fn add_additional_factors_to_cff_atom(&self, cff_atom: &Atom) -> Atom {
        let inverse_energy_product = self.graph.underlying.get_cff_inverse_energy_product();
        let factors_of_pi = parse!(&format!(
            "(2*π)^{}",
            3 * self.graph.underlying.get_loop_number()
        ))
        .unwrap();

        let result = cff_atom * inverse_energy_product * &self.graph.multiplicity / factors_of_pi;
        debug!("result: {}", result);

        result
    }

    fn build_evaluator(&mut self) {
        let atom_unsubstituted = self.derived_data.cff_expression.as_ref().unwrap().to_atom();
        let atom = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .substitute_energies(&atom_unsubstituted);

        let atom = self.add_additional_factors_to_cff_atom(&atom);

        let params = self.get_params();
        let function_map = self.get_function_map();

        let mut tree = atom.to_evaluation_tree(&function_map, &params).unwrap();
        tree.horner_scheme();
        tree.common_subexpression_elimination();

        let tree_double = tree.map_coeff::<F<f64>, _>(&|r| r.into());
        let tree_quad = tree.map_coeff::<F<f128>, _>(&|r| r.into());

        let evaluator = GenericEvaluator {
            f64_compiled: None,
            f64_eager: RefCell::new(tree_double.linearize(Some(1))),
            f128: RefCell::new(tree_quad.linearize(Some(1))),
        };

        self.derived_data.bare_cff_evaluator = Some(evaluator)
    }

    fn build_loop_momentum_bases(&mut self) {
        let lmbs = self
            .graph
            .loop_momentum_basis
            .generate_loop_momentum_bases(&self.graph.underlying);

        self.derived_data.lmbs = Some(lmbs)
    }

    fn build_tropical_sampler(&mut self, process_settings: &ProcessSettings) -> Result<()> {
        let num_virtual_loop_edges = self.graph.iter_loop_edges().count();
        let num_loops = self.graph.loop_momentum_basis.basis.len();
        let target_omega = process_settings
            .tropical_subgraph_table_settings
            .target_omega;

        let weight = (target_omega + (3 * num_loops) as f64 / 2.) / num_virtual_loop_edges as f64;

        debug!(
            "Building tropical subgraph table with all edge weights set to: {}",
            weight
        );

        let tropical_edges = self
            .graph
            .iter_loop_edges()
            .map(|(pair, _edge_id, edge)| {
                let is_massive = match edge.data.particle.0.mass.value {
                    Some(complex_mass) => complex_mass.is_non_zero(),
                    None => false,
                };

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

        Ok(self.derived_data.tropical_sampler = Some(sampler))
    }

    fn build_evaluator_for_orientations(&mut self) -> Result<()> {
        let params = self.get_params();
        let function_map = self.get_function_map();

        let evaluators = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .get_orientation_atoms()
            .iter()
            .map(|orientation_atom_unsubstituted| {
                let atom_no_prefactor = self
                    .derived_data
                    .cff_expression
                    .as_ref()
                    .unwrap()
                    .surfaces
                    .substitute_energies(orientation_atom_unsubstituted);

                let atom = self.add_additional_factors_to_cff_atom(&atom_no_prefactor);

                let mut tree = atom.to_evaluation_tree(&function_map, &params).unwrap();
                tree.horner_scheme();
                tree.common_subexpression_elimination();

                let tree_double = tree.map_coeff::<F<f64>, _>(&|r| r.into());
                let tree_quad = tree.map_coeff::<F<f128>, _>(&|r| r.into());

                GenericEvaluator {
                    f64_compiled: None,
                    f64_eager: RefCell::new(tree_double.linearize(Some(1))),
                    f128: RefCell::new(tree_quad.linearize(Some(1))),
                }
            })
            .collect();

        Ok(self.derived_data.bare_cff_orientation_evaluatos = Some(evaluators))
    }

    fn generate_term_for_graph(&self, _settings: &Settings) -> AmplitudeGraphTerm {
        AmplitudeGraphTerm {
            bare_cff_evaluator: self.derived_data.bare_cff_evaluator.clone().unwrap(),
            bare_cff_orientation_evaluators: self
                .derived_data
                .bare_cff_orientation_evaluatos
                .clone()
                .unwrap(),
            tropical_sampler: self.derived_data.tropical_sampler.clone().unwrap(),
            graph: self.graph.clone(),
            multi_channeling_setup: self.derived_data.multi_channeling_setup.clone().unwrap(),
            lmbs: self.derived_data.lmbs.clone().unwrap(),
        }
    }
}

#[derive(Clone, Encode, Decode)]
pub struct AmplitudeDerivedData<S: NumeratorState> {
    cff_expression: Option<CFFExpression>,
    bare_cff_evaluator: Option<GenericEvaluator>,
    bare_cff_orientation_evaluatos: Option<TiVec<OrientationID, GenericEvaluator>>,
    _temp_numerator: Option<PhantomData<S>>,
    lmbs: Option<TiVec<LmbIndex, LoopMomentumBasis>>,
    tropical_sampler: Option<SampleGenerator<3>>,
    multi_channeling_setup: Option<LmbMultiChannelingSetup>,
}

impl<S: NumeratorState> Amplitude<S> {
    pub fn new(name: String) -> Self {
        Self {
            name,
            graphs: vec![],
            external_particles: vec![],
            external_signature: SignatureLike::from_iter(iter::empty::<i8>()),
        }
    }

    pub fn add_graph(&mut self, graph: Graph) -> Result<()> {
        let new_external_particels = graph.underlying.get_external_partcles();
        let new_external_signature = graph.underlying.get_external_signature();

        if !self.graphs.is_empty() {
            if self.external_particles != new_external_particels {
                return Err(eyre!("amplitude graph has different number of externals"));
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

#[derive(Clone)]
pub struct CrossSection<S: NumeratorState> {
    name: String,
    supergraphs: Vec<CrossSectionGraph<S>>,
    external_particles: Vec<ArcParticle>,
    external_connections: Vec<ExternalConnection>,
    n_incmoming: usize,
}

impl<S: NumeratorState> CrossSection<S> {
    pub fn new(name: String) -> Self {
        Self {
            name,
            supergraphs: vec![],
            external_connections: vec![],
            external_particles: vec![],
            n_incmoming: 0,
        }
    }

    pub fn add_supergraph(&mut self, supergraph: Graph) -> Result<()> {
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

        if self.external_connections.is_empty() {
            self.external_connections = supergraph
                .external_connections
                .as_ref()
                .expect("the definition of a cross section requires external connections")
                .clone();
        } else if &self.external_connections
            != supergraph
                .external_connections
                .as_ref()
                .expect("the definition of a cross section required external connections")
        {
            return Err(eyre!(
                "attempt to add supergraph with different external connections"
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
    ) -> Result<()> {
        for supergraph in &mut self.supergraphs {
            supergraph.preprocess(model, process_definition)?;
        }
        Ok(())
    }

    pub fn generate_integrand(&self, settings: Settings) -> Integrand {
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

        let cross_section_integrand = CrossSectionIntegrand {
            rotations,
            external_connections: self.external_connections.clone(),
            n_incoming: self.n_incmoming,
            polarizations,
            settings,
            graph_terms: terms,
        };

        Integrand::NewIntegrand(NewIntegrand::CrossSection(cross_section_integrand))
    }

    fn export(&self, settings: &ExportSettings) -> Result<()> {
        let path = Path::new(&settings.root_folder)
            .join("sources")
            .join("cross_sections")
            .join(self.name.as_str());

        for supergraph in self.supergraphs.iter() {
            let file_name = path
                .clone()
                .join(format!("cross_section_graph_{}.bin", supergraph.graph.name));

            todo!("implement export of cross section graph")

            //            let data = bincode::encode_to_vec(supergraph, bincode::config::standard())?;
            //      std::fs::write(file_name, &data)?;
        }

        Ok(())
    }
}

impl<S: NumeratorState> IsPolarizable for CrossSection<S> {
    fn polarizations(&self, externals: &Externals) -> Polarizations {
        externals.generate_polarizations(
            &self.external_particles,
            DependentMomentaConstructor::CrossSection {
                external_connections: &self.external_connections,
            },
        )
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, From, Into, Hash, PartialEq, Copy, Eq, Encode)]
pub struct CutId(usize);

impl Display for CutId {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

#[derive(Clone)]
pub struct CrossSectionGraph<S: NumeratorState = PythonState> {
    pub graph: Graph,
    pub source_nodes: HedgeNode,
    pub target_nodes: HedgeNode,
    pub cuts: TiVec<CutId, CrossSectionCut>,
    pub cut_esurface: TiVec<CutId, Esurface>,
    pub cut_esurface_id_map: TiVec<CutId, EsurfaceID>,
    pub derived_data: CrossSectionDerivedData<S>,
}

impl<S: NumeratorState> CrossSectionGraph<S> {
    fn new(graph: Graph) -> Self {
        let mut source_nodes = AHashSet::new();
        let mut target_nodes = AHashSet::new();

        for (hedge_pair, _, _) in graph.underlying.iter_all_edges() {
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

    pub fn preprocess(
        &mut self,
        model: &Model,
        process_definition: &ProcessDefinition,
    ) -> Result<()> {
        // not sure why I need to call this again, but it seems like I do?
        self.graph
            .loop_momentum_basis
            .set_edge_signatures(&self.graph.underlying)?;
        self.generate_cuts(model, process_definition)?;
        self.generate_esurface_cuts();
        self.generate_cff();
        self.update_surface_cache();

        for (cut_id, esurface) in self.cut_esurface.iter_enumerated() {
            debug!(
                "cut_id: {:?} \n, esurface: {:#?} \n, expression: {}",
                cut_id,
                esurface,
                self.derived_data
                    .cff_expression
                    .as_ref()
                    .unwrap()
                    .to_atom_for_cut(cut_id)
            );
        }

        self.build_cut_evaluators();
        self.build_orientation_evaluators();
        self.build_lmbs();
        Ok(self.build_multi_channeling_channels())
    }

    pub fn update_surface_cache(&mut self) {
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

    fn generate_cff(&mut self) {
        // hardcorde 1 to n for now
        debug!("generating cff");

        let shift_rewrite = self
            .graph
            .underlying
            .get_esurface_canonization(&self.graph.loop_momentum_basis);

        let cff_cut_expression =
            generate_cff_with_cuts(&self.graph.underlying, &shift_rewrite, &self.cuts);

        self.derived_data.cff_expression = Some(cff_cut_expression)
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

        let cuts: TiVec<CutId, CrossSectionCut> = all_st_cuts
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

    fn build_atom_for_cut(&self, cut_id: CutId) -> Atom {
        let cut_atom = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .to_atom_for_cut(cut_id);

        let cut_atom_energy_sub = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .surfaces
            .substitute_energies(&cut_atom);

        self.add_additional_factors_to_cff_atom(&cut_atom_energy_sub)
    }

    fn add_additional_factors_to_cff_atom(&self, cut_atom: &Atom) -> Atom {
        let inverse_energy_product = self.graph.underlying.get_cff_inverse_energy_product();

        let t_star_factor = parse!(&format!(
            "tstar^({})",
            3 * self.graph.underlying.get_loop_number()
        ))
        .unwrap();

        let h_function = parse!("h").unwrap();
        let grad_eta = parse!("∇η").unwrap();

        let factors_of_pi = parse!(&format!(
            "(2*π)^{}",
            3 * self.graph.underlying.get_loop_number()
        ))
        .unwrap();

        let result = cut_atom
            * inverse_energy_product
            * t_star_factor
            * h_function
            * &self.graph.multiplicity
            / grad_eta
            / factors_of_pi;

        debug!("result: {}", result);
        result
    }

    fn get_params(&self) -> Vec<Atom> {
        let mut params = self.graph.underlying.get_energy_atoms();

        params.push(parse!("tstar").unwrap());
        params.push(parse!("h").unwrap());
        params.push(parse!("∇η").unwrap());

        params
    }

    fn get_function_map(&self) -> FunctionMap {
        let mut fn_map = FunctionMap::new();
        let pi_rational = Rational::from(std::f64::consts::PI);
        fn_map.add_constant(parse!("π").unwrap(), pi_rational);
        fn_map
    }

    fn build_cut_evaluators(&mut self) {
        let evaluators = self
            .cuts
            .iter_enumerated()
            .map(|(cut_id, _)| {
                let atom = self.build_atom_for_cut(cut_id);

                let params = self.get_params();
                let mut tree = atom
                    .to_evaluation_tree(&self.get_function_map(), &params)
                    .unwrap();

                tree.horner_scheme();
                tree.common_subexpression_elimination();

                let tree_double = tree.map_coeff::<F<f64>, _>(&|r| r.into());
                let tree_quad = tree.map_coeff::<F<f128>, _>(&|r| r.into());

                GenericEvaluator {
                    f64_compiled: None,
                    f64_eager: RefCell::new(tree_double.linearize(Some(1))),
                    f128: RefCell::new(tree_quad.linearize(Some(1))),
                }
            })
            .collect();

        self.derived_data.bare_cff_evaluators = Some(evaluators)
    }

    fn build_orientation_evaluators(&mut self) {
        let orientation_atoms = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .get_orientation_atoms();
        let substituted_energies = orientation_atoms
            .iter()
            .map(|cut_atoms| {
                cut_atoms
                    .iter()
                    .map(|atom| {
                        let substituded_atom = self
                            .derived_data
                            .cff_expression
                            .as_ref()
                            .unwrap()
                            .surfaces
                            .substitute_energies(atom);
                        substituded_atom
                    })
                    .collect_vec()
            })
            .collect::<TiVec<OrientationID, _>>();

        let orientation_datas = self
            .derived_data
            .cff_expression
            .as_ref()
            .unwrap()
            .orientations
            .iter()
            .map(|orienatation| orienatation.data.clone());

        let orientation_evaluators = substituted_energies
            .iter()
            .zip(orientation_datas)
            .map(|(cut_atoms, orientation_data)| {
                let cut_evaluators = cut_atoms
                    .iter()
                    .map(|cut_atom| {
                        let cut_atom = self.add_additional_factors_to_cff_atom(cut_atom);
                        let params = self.get_params();
                        let mut tree = cut_atom
                            .to_evaluation_tree(&self.get_function_map(), &params)
                            .unwrap();

                        tree.horner_scheme();
                        tree.common_subexpression_elimination();

                        let tree_double = tree.map_coeff::<F<f64>, _>(&|r| r.into());
                        let tree_quad = tree.map_coeff::<F<f128>, _>(&|r| r.into());

                        GenericEvaluator {
                            f64_compiled: None,
                            f64_eager: RefCell::new(tree_double.linearize(Some(1))),
                            f128: RefCell::new(tree_quad.linearize(Some(1))),
                        }
                    })
                    .collect_vec();

                OrientationEvaluator {
                    orientation_data,
                    evaluators: cut_evaluators,
                }
            })
            .collect();

        self.derived_data.bare_cff_orientation_evaluators = Some(orientation_evaluators);
    }

    fn build_lmbs(&mut self) {
        let lmbs = self
            .graph
            .loop_momentum_basis
            .generate_loop_momentum_bases(&self.graph.underlying);

        self.derived_data.lmbs = Some(lmbs)
    }

    fn build_multi_channeling_channels(&mut self) {
        let channels = self
            .graph
            .build_multi_channeling_channels(self.derived_data.lmbs.as_ref().unwrap());

        self.derived_data.multi_channeling_setup = Some(channels)
    }

    fn generate_term_for_graph(&self, _settings: &Settings) -> CrossSectionGraphTerm {
        CrossSectionGraphTerm {
            multi_channeling_setup: self.derived_data.multi_channeling_setup.clone().unwrap(),
            bare_cff_evaluators: self.derived_data.bare_cff_evaluators.clone().unwrap(),
            bare_cff_orientation_evaluators: self
                .derived_data
                .bare_cff_orientation_evaluators
                .clone()
                .unwrap(),
            graph: self.graph.clone(),
            cut_esurface: self.cut_esurface.clone(),
            lmbs: self.derived_data.lmbs.clone().unwrap(),
        }
    }
}

#[derive(Clone, Encode)]
pub struct CrossSectionDerivedData<S: NumeratorState> {
    orientations: Option<TiVec<OrientationID, CutOrientationData>>,
    bare_cff_evaluators: Option<TiVec<CutId, GenericEvaluator>>,
    bare_cff_orientation_evaluators: Option<TiVec<OrientationID, OrientationEvaluator>>,
    cff_expression: Option<CFFCutExpression>,
    lmbs: Option<TiVec<LmbIndex, LoopMomentumBasis>>,
    multi_channeling_setup: Option<LmbMultiChannelingSetup>,
    _temp_numerator: Option<PhantomData<S>>,
}

impl<S: NumeratorState> CrossSectionDerivedData<S> {
    fn new_empty() -> Self {
        Self {
            orientations: None,
            cff_expression: None,
            _temp_numerator: None,
            bare_cff_evaluators: None,
            bare_cff_orientation_evaluators: None,
            lmbs: None,
            multi_channeling_setup: None,
        }
    }
}

#[derive(Clone)]
pub struct CrossSectionCut {
    pub cut: OrientedCut,
    pub left: BitVec,
    pub right: BitVec,
}

impl CrossSectionCut {
    pub fn is_s_channel<S: NumeratorState>(
        &self,
        cross_section_graph: &CrossSectionGraph<S>,
    ) -> Result<bool> {
        let nodes_of_left_cut: Vec<_> = cross_section_graph
            .graph
            .underlying
            .iter_node_data(&self.left)
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

    pub fn is_valid_for_process<S: NumeratorState>(
        &self,
        cross_section_graph: &CrossSectionGraph<S>,
        process: &ProcessDefinition,
        model: &Model,
    ) -> Result<bool> {
        if self.is_s_channel(cross_section_graph)? {
            let mut cut_content = self
                .cut
                .iter_edges(&cross_section_graph.graph.underlying)
                .map(|(orientation, edge_data)| {
                    if orientation == Orientation::Reversed {
                        edge_data.data.particle.0.get_anti_particle(model)
                    } else {
                        edge_data.data.particle.clone()
                    }
                })
                .collect_vec();

            let particle_content = process
                .final_pdgs_lists
                .first()
                .unwrap()
                .iter()
                .map(|pdg| model.get_particle_from_pdg(*pdg as isize));

            debug!(
                "cut content: {:?}",
                cut_content.iter().map(|p| p.0.pdg_code).collect_vec()
            );

            for particle in particle_content {
                if let Some(index) = cut_content.iter().position(|p| p == &particle) {
                    cut_content.remove(index);
                } else {
                    debug!("wrong particles");
                    return Ok(false);
                }
            }

            if cut_content.len() > process.n_unresolved {
                debug!(" too many unresolved particles");
                return Ok(false);
            }

            if !cut_content
                .iter()
                .all(|particle| process.unresolved_cut_content.contains(particle))
            {
                debug!("wrong unresolved particles");
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
