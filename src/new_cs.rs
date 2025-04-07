use std::{
    cell::{Ref, RefCell},
    collections::BTreeMap,
    fmt::{Display, Formatter},
    marker::PhantomData,
    path::PathBuf,
    sync::Arc,
};

use ahash::{AHashSet, HashMap};
use bincode::de;
use bitvec::vec::BitVec;
use color_eyre::Result;

use crate::{new_gammaloop_integrand::GenericEvaluator, utils::f128};
use eyre::eyre;
use itertools::Itertools;
use linnet::half_edge::{
    builder::HedgeGraphBuilder,
    involution::{Flow, HedgePair, Orientation},
    subgraph::{HedgeNode, Inclusion, InternalSubGraph, OrientedCut},
    HedgeGraph, NodeIndex,
};
use log::{debug, warn};
use serde::{Deserialize, Serialize};
use spenso::upgrading_arithmetic::GreaterThan;
use symbolica::{
    atom::{Atom, AtomCore, Symbol},
    coefficient::ConvertToRing,
    domains::{float::NumericalFloatLike, rational::Rational},
    evaluate::{ExpressionEvaluator, FunctionMap, OptimizationSettings},
    parse, symbol,
};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cut_expression::{
            CFFCutExpression, CutOrientationExpression, OrientationData, OrientationID,
        },
        esurface::{self, Esurface, EsurfaceID},
        generation::{generate_cff_with_cuts, ShiftRewrite},
    },
    cross_section::{self, IsPolarizable},
    disable,
    feyngen::{
        diagram_generator::FeynGen, FeynGenFilter, FeynGenFilters, FeynGenOptions, GenerationType,
        NumeratorAwareGraphGroupingOption, SelfEnergyFilterOptions, SnailFilterOptions,
        TadpolesFilterOptions,
    },
    graph::{self, BareGraph, BareVertex, EdgeType},
    integrands::Integrand,
    model::{self, Model, Particle},
    momentum::{Rotatable, Rotation, RotationMethod},
    momentum_sample::ExternalIndex,
    new_gammaloop_integrand::{
        cross_section_integrand::{self, CrossSectionGraphTerm, CrossSectionIntegrand},
        NewIntegrand,
    },
    new_graph::{Edge, ExternalConnection, FeynmanGraph, Graph, Vertex},
    numerator::{GlobalPrefactor, NumeratorState, PythonState},
    utils::{F, GS},
    DependentMomentaConstructor, Externals, Polarizations, ProcessSettings, Settings,
};

use derive_more::{From, Into};

#[derive(Debug, Clone)]
pub struct ProcessDefinition {
    pub initial_pdgs: Vec<i64>, // Do we want a pub type Pdg = i64;?
    pub final_pdgs_lists: Vec<Vec<i64>>,
    pub n_unresolved: usize, // we need al this information to know what cuts are considered at runtime
    pub unresolved_cut_content: AHashSet<Arc<Particle>>,
    pub amplitude_filters: FeynGenFilters,
    pub cross_section_filters: FeynGenFilters,
}

#[derive(Clone)]
pub struct Process<S: NumeratorState = PythonState> {
    pub definition: ProcessDefinition,
    pub collection: ProcessCollection<S>,
}

impl<S: NumeratorState> Process<S> {
    pub fn preprocess(&mut self, model: &Model) -> Result<()> {
        self.collection.preprocess(model, &self.definition)?;
        Ok(())
    }
}

pub struct GenerationOptions {
    pub feyngen_options: FeynGenOptions,
    pub numerator_aware_isomorphism_grouping: NumeratorAwareGraphGroupingOption,
    pub filter_self_loop: bool,
    pub graph_prefix: String,
    pub selected_graphs: Option<Vec<String>>,
    pub vetoed_graphs: Option<Vec<String>>,
    pub loop_momentum_bases: Option<HashMap<String, Vec<String>>>,
    pub global_prefactor: GlobalPrefactor,
    pub num_threads: Option<usize>,
}

// should produce e+e- -> d d~ g
#[cfg(test)]
fn test_process_cross_section() -> GenerationOptions {
    let feyngen_options = FeynGenOptions {
        amplitude_filters: FeynGenFilters(vec![]),
        cross_section_filters: FeynGenFilters(vec![
            FeynGenFilter::SelfEnergyFilter(SelfEnergyFilterOptions::default()),
            FeynGenFilter::TadpolesFilter(TadpolesFilterOptions::default()),
            FeynGenFilter::ZeroSnailsFilter(SnailFilterOptions::default()),
        ]),
        initial_pdgs: vec![11, -11],
        final_pdgs_lists: vec![vec![1, -1, 21]],
        generation_type: GenerationType::CrossSection,
        symmetrize_final_states: false,
        symmetrize_initial_states: false,
        symmetrize_left_right_states: false,
        loop_count_range: (2, 2),
        allow_symmetrization_of_external_fermions_in_amplitudes: false,
        max_multiplicity_for_fast_cut_filter: 6,
    };

    GenerationOptions {
        feyngen_options,
        numerator_aware_isomorphism_grouping: NumeratorAwareGraphGroupingOption::OnlyDetectZeroes,
        filter_self_loop: true,
        graph_prefix: "GL".to_string(),
        selected_graphs: None,
        vetoed_graphs: None,
        loop_momentum_bases: None,
        global_prefactor: GlobalPrefactor::default(),
        num_threads: None,
    }
}

#[cfg(test)]
fn test_process_amplitude() -> GenerationOptions {
    use crate::numerator::GlobalPrefactor;

    let feyngen_options = FeynGenOptions {
        amplitude_filters: FeynGenFilters(vec![
            FeynGenFilter::SelfEnergyFilter(SelfEnergyFilterOptions::default()),
            FeynGenFilter::TadpolesFilter(TadpolesFilterOptions::default()),
            FeynGenFilter::ZeroSnailsFilter(SnailFilterOptions::default()),
        ]),
        cross_section_filters: FeynGenFilters(vec![]),
        generation_type: GenerationType::Amplitude,
        initial_pdgs: vec![22, 22],
        final_pdgs_lists: vec![vec![22, 22]],
        symmetrize_final_states: true,
        symmetrize_initial_states: true,
        symmetrize_left_right_states: false,
        allow_symmetrization_of_external_fermions_in_amplitudes: false,
        loop_count_range: (1, 1),
        max_multiplicity_for_fast_cut_filter: 6,
    };

    GenerationOptions {
        feyngen_options,
        numerator_aware_isomorphism_grouping: NumeratorAwareGraphGroupingOption::OnlyDetectZeroes,
        filter_self_loop: true,
        graph_prefix: "GL".to_string(),
        selected_graphs: None,
        vetoed_graphs: None,
        loop_momentum_bases: None,
        global_prefactor: GlobalPrefactor::default(),
        num_threads: None,
    }
}

impl Process {
    pub fn generate(options: GenerationOptions, model: &Model) -> Result<Self> {
        disable! {
        let definition = ProcessDefinition {
            initial_pdgs: options.feyngen_options.initial_pdgs.clone(),
            final_pdgs: options.feyngen_options.final_pdgs.clone(),
        };

        let diagram_generator = FeynGen::new(options.feyngen_options.clone());

        let bare_collection = diagram_generator.generate(
            model,
            &options.numerator_aware_isomorphism_grouping,
            options.filter_self_loop,
            options.graph_prefix,
            options.selected_graphs,
            options.vetoed_graphs,
            options.loop_momentum_bases,
            options.global_prefactor,
            options.num_threads,
        )?;

        let hedge_collection = bare_collection
            .into_iter()
            .map(|bare_graph| {
                let mut hedge_graph_builder = HedgeGraphBuilder::<Edge, Vertex>::new();
                for vertex in bare_graph.vertices {
                    // do not add the vertices attached to externals
                    if vertex.edges.len() != 1 {
                        hedge_graph_builder.add_node(vertex.clone().into());
                    }
                }

                for edge in bare_graph.edges {
                    match edge.edge_type {
                        EdgeType::Virtual => {
                            hedge_graph_builder.add_edge(
                                NodeIndex(edge.vertices[0]),
                                NodeIndex(edge.vertices[1]),
                                edge.into(),
                                true,
                            );
                        }
                        EdgeType::Incoming => {
                            hedge_graph_builder.add_external_edge(
                                NodeIndex(edge.vertices[1]),
                                edge.into(),
                                // I am not sure what to put for these two
                                true,
                                Flow::Source,
                            );
                        }
                        EdgeType::Outgoing => {
                            hedge_graph_builder.add_external_edge(
                                NodeIndex(edge.vertices[0]),
                                edge.into(),
                                // I am not sure what to put for these two
                                true,
                                Flow::Sink,
                            );
                        }
                    }
                }

                Graph::new(bare_graph.overall_factor, hedge_graph_builder.build())
            })
            .collect::<Result<Vec<Graph>, _>>()?;

        }
        todo!();
    }

    pub fn from_bare_graph_list(
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
                    let mut amplitude: Amplitude<PythonState> = Amplitude::new();

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
                    let mut cross_section: CrossSection<PythonState> = CrossSection::new();

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

struct ExportSettings {
    root_folder: PathBuf,
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

    /// given the process definition generates a new process and adds it to the list
    pub fn generate(&mut self, options: GenerationOptions, model: &Model) -> Result<()> {
        self.processes.push(Process::generate(options, model)?);
        Ok(())
    }

    /// imports a process list from a folder
    pub fn import(settings: ExportSettings) -> Result<Self> {
        Ok(Self::new())
    }

    ///preprocesses the process list according to the settings
    pub fn preprocess(&mut self, model: &Model, settings: ProcessSettings) -> Result<()> {
        for process in self.processes.iter_mut() {
            process.preprocess(&model)?;
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
    pub fn export(&self, settings: ExportSettings) -> Result<()> {
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

    fn preprocess(&mut self, model: &Model, process_definition: &ProcessDefinition) -> Result<()> {
        match self {
            Self::Amplitudes(amplitudes) => {
                todo!()
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
            Self::Amplitudes(_) => todo!(),
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
}

#[derive(Clone)]
pub struct Amplitude<S: NumeratorState = PythonState> {
    graphs: Vec<AmplitudeGraph<S>>,
}

#[derive(Clone)]
pub struct AmplitudeGraph<S: NumeratorState = PythonState> {
    graph: Graph,
    derived_data: AmplitudeDerivedData<S>,
}

impl<S: NumeratorState> AmplitudeGraph<S> {
    fn new(graph: Graph) -> Self {
        AmplitudeGraph {
            graph,
            derived_data: AmplitudeDerivedData {
                temp_numerator: PhantomData,
            },
        }
    }
}

#[derive(Clone)]
pub struct AmplitudeDerivedData<S: NumeratorState> {
    temp_numerator: PhantomData<S>,
}

impl<S: NumeratorState> Amplitude<S> {
    pub fn new() -> Self {
        Self { graphs: vec![] }
    }

    pub fn add_graph(&mut self, graph: Graph) -> Result<()> {
        self.graphs.push(AmplitudeGraph::new(graph));
        /// TODO: validate that the graph is compatible
        Ok(())
    }
}

#[derive(Clone)]
pub struct CrossSection<S: NumeratorState> {
    supergraphs: Vec<CrossSectionGraph<S>>,
    external_particles: Vec<Arc<Particle>>,
    external_connections: Vec<ExternalConnection>,
    n_incmoming: usize,
}

impl<S: NumeratorState> CrossSection<S> {
    pub fn new() -> Self {
        Self {
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
            external_connections: self.external_connections.clone(),
            n_incoming: self.n_incmoming,
            polarizations,
            settings,
            graph_terms: terms,
        };

        Integrand::NewIntegrand(NewIntegrand::CrossSection(cross_section_integrand))
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

#[derive(Debug, Clone, Serialize, Deserialize, From, Into, Hash, PartialEq, Copy, Eq)]
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

        let num_incoming = source_nodes.len();

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
                self.derived_data.cff_expression.to_atom_for_cut(cut_id)
            );
        }
        Ok(())
    }

    pub fn update_surface_cache(&mut self) {
        let esurface_cache = &mut self.derived_data.cff_expression.surfaces.esurface_cache;

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

        self.derived_data.cff_expression = cff_cut_expression;
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
        let cut_atom = self.derived_data.cff_expression.to_atom_for_cut(cut_id);

        let cut_atom_energy_sub = self
            .derived_data
            .cff_expression
            .surfaces
            .substitute_energies(&cut_atom);

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

        let result = cut_atom_energy_sub
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
        let mut params = self
            .graph
            .underlying
            .iter_all_edges()
            .map(|(pair, edge_id, _)| match pair {
                HedgePair::Paired { .. } => {
                    parse!(&format!("Q({}, cind(0))", Into::<usize>::into(edge_id))).unwrap()
                }
                HedgePair::Unpaired { .. } => {
                    parse!(&format!("P({}, cind(0))", Into::<usize>::into(edge_id))).unwrap()
                }
                _ => unreachable!(),
            })
            .collect_vec();

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

    fn generate_term_for_graph(&self, settings: &Settings) -> CrossSectionGraphTerm {
        let mut evaluators = TiVec::new();

        for (cut_id, _) in self.cuts.iter_enumerated() {
            let atom = self.build_atom_for_cut(cut_id);

            let params = self.get_params();
            let mut tree = atom
                .to_evaluation_tree(&self.get_function_map(), &params)
                .unwrap();

            tree.horner_scheme();
            tree.common_subexpression_elimination();

            let tree_double = tree.map_coeff::<F<f64>, _>(&|r| r.into());
            let tree_quad = tree.map_coeff::<F<f128>, _>(&|r| r.into());

            let generic_evaluator = GenericEvaluator {
                f64_compiled: None,
                f64_eager: RefCell::new(tree_double.linearize(Some(1))),
                f128: RefCell::new(tree_quad.linearize(Some(1))),
            };

            evaluators.push(generic_evaluator);
        }

        CrossSectionGraphTerm {
            bare_cff_evaluators: evaluators,
            graph: self.graph.clone(),
            cut_esurface: self.cut_esurface.clone(),
        }
    }
}

#[derive(Clone)]
pub struct CrossSectionDerivedData<S: NumeratorState = PythonState> {
    pub orientations: TiVec<OrientationID, OrientationData>,
    pub cff_expression: CFFCutExpression,
    temp_numerator: PhantomData<S>,
}

impl<S: NumeratorState> CrossSectionDerivedData<S> {
    fn new_empty() -> Self {
        Self {
            orientations: TiVec::new(),
            cff_expression: CFFCutExpression::new_empty(),
            temp_numerator: PhantomData,
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
        let nodes_of_left_cut: Option<Vec<_>> = cross_section_graph
            .graph
            .underlying
            .iter_node_data(&self.left)
            .map(|(hedge_node, _)| {
                cross_section_graph
                    .graph
                    .underlying
                    .id_from_hairs(hedge_node)
            })
            .collect();

        if let Some(nodes_left_of_cut) = nodes_of_left_cut {
            let left_node = cross_section_graph
                .graph
                .underlying
                .combine_to_single_hedgenode(&nodes_left_of_cut);
            let res = left_node.includes(&cross_section_graph.source_nodes)
                && cross_section_graph.target_nodes.weakly_disjoint(&left_node);

            if !res {
                warn!("s channel check wrong");
            }

            Ok(true)
        } else {
            Err(eyre!("Could not determine nodes in left cut"))
        }
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
                        edge_data.data.particle.get_anti_particle(model)
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
                cut_content.iter().map(|p| p.pdg_code).collect_vec()
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

#[test]
fn test_new_structure_cross_section() {
    let model = crate::tests_from_pytest::load_generic_model("sm");
    let generation_options = test_process_cross_section();

    let mut process_list = ProcessList::new();
    process_list.generate(generation_options, &model).unwrap();
}

#[test]
fn test_new_structure_ampltiude() {
    let model = crate::tests_from_pytest::load_generic_model("sm");
    let generation_options = test_process_amplitude();

    let mut process_list = ProcessList::new();
    process_list.generate(generation_options, &model).unwrap();
}
