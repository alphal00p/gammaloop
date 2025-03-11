use std::{collections::btree_map::Range, marker::PhantomData, path::PathBuf, sync::Arc};

use ahash::{AHashSet, HashMap};
use bitvec::{index, vec::BitVec};
use color_eyre::Result;

use eyre::eyre;
use hyperdual::Num;
use itertools::{partition, Itertools};
use linnet::half_edge::{
    builder::HedgeGraphBuilder,
    involution::{Flow, HedgePair, Orientation},
    subgraph::{HedgeNode, InternalSubGraph, OrientedCut},
    HedgeGraph, NodeIndex,
};
use serde::{Deserialize, Serialize};
use spenso::upgrading_arithmetic::GreaterThan;
use symbolica::{atom::Atom, parse};
use typed_index_collections::TiVec;

use crate::{
    cff::{
        cut_expression::{CutOrientationExpression, OrientationData, OrientationID},
        esurface::Esurface,
    },
    cross_section, disable,
    feyngen::{
        diagram_generator::FeynGen, FeynGenFilter, FeynGenFilters, FeynGenOptions, GenerationType,
        NumeratorAwareGraphGroupingOption, SelfEnergyFilterOptions, SnailFilterOptions,
        TadpolesFilterOptions,
    },
    graph::{self, BareGraph, BareVertex, EdgeType},
    model::{self, Model, Particle},
    new_graph::{Edge, Graph, Vertex},
    numerator::{GlobalPrefactor, NumeratorState, PythonState},
    ProcessSettings,
};

use derive_more::{From, Into};

#[derive(Debug, Clone)]
pub struct ProcessDefinition {
    pub initial_pdgs: Vec<i64>, // Do we want a pub type Pdg = i64;?
    pub final_pdgs: Vec<i64>,
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
        final_pdgs: vec![1, -1, 21],
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
        final_pdgs: vec![22, 22],
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
        let graphs = bare_graphs
            .into_iter()
            .map(|bare_graph| Graph::from(bare_graph))
            .collect_vec();

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
    pub fn preprocess(&mut self, settings: ProcessSettings) -> Result<()> {
        Ok(())
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
}

impl<S: NumeratorState> CrossSection<S> {
    pub fn new() -> Self {
        Self {
            supergraphs: vec![],
        }
    }

    pub fn add_supergraph(&mut self, supergraph: Graph) -> Result<()> {
        let mut cross_section_graph = CrossSectionGraph::new(supergraph);
        self.supergraphs.push(cross_section_graph);
        /// TODO: validate that the graph is compatible
        Ok(())
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, From, Into, Hash, PartialEq, Copy, Eq)]
pub struct CutId(usize);

#[derive(Clone)]
pub struct CrossSectionGraph<S: NumeratorState = PythonState> {
    graph: Graph,
    source_nodes: AHashSet<NodeIndex>,
    target_nodes: AHashSet<NodeIndex>,
    cuts: TiVec<CutId, CrossSectionCut>,
    cut_esurface: TiVec<CutId, Esurface>,
    derived_data: CrossSectionDerivedData<S>,
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

        Self {
            graph,
            source_nodes,
            target_nodes,
            cuts: TiVec::new(),
            cut_esurface: TiVec::new(),
            derived_data: CrossSectionDerivedData::<S>::new_empty(),
        }
    }

    fn generate_cuts(
        &mut self,
        model: &Model,
        process_definition: ProcessDefinition,
    ) -> Result<()> {
        if let (Some(&source), Some(&target)) = (
            self.source_nodes.iter().next(),
            self.target_nodes.iter().next(),
        ) {
            let cuts: TiVec<CutId, CrossSectionCut> = self
                .graph
                .underlying
                .all_cuts(source, target)
                .into_iter()
                .map(|(left, cut, right)| CrossSectionCut { cut, left, right })
                .filter_map(|cut| {
                    match cut.is_valid_for_process(self, &process_definition, model) {
                        Ok(true) => Some(Ok(cut)),
                        Ok(false) => None,
                        Err(e) => Some(Err(e)),
                    }
                })
                .collect::<Result<_>>()?;

            self.cuts = cuts;
            Ok(())
        } else {
            Err(eyre!("Could not find cuts for graph: {}", self.graph.name))
        }
    }

    fn generate_esurface_cuts(&mut self) {
        let esurfaces: TiVec<CutId, Esurface> = self
            .cuts
            .iter()
            .map(|cut| Esurface::new_from_cut_left(&self.graph.underlying, cut))
            .collect();

        self.cut_esurface = esurfaces;
    }
}

#[derive(Clone)]
pub struct CrossSectionDerivedData<S: NumeratorState = PythonState> {
    pub orientations: TiVec<OrientationID, OrientationData>,
    temp_numerator: PhantomData<S>,
}

impl<S: NumeratorState> CrossSectionDerivedData<S> {
    fn new_empty() -> Self {
        Self {
            orientations: TiVec::new(),
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
        let nodes_of_left_cut: Option<AHashSet<_>> = cross_section_graph
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
            Ok(cross_section_graph
                .source_nodes
                .is_subset(&nodes_left_of_cut)
                && cross_section_graph
                    .target_nodes
                    .intersection(&nodes_left_of_cut)
                    .count()
                    == 0)
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
                .iter_edges_relative(&cross_section_graph.graph.underlying)
                .map(|(orientation, edge_data)| {
                    if orientation == Orientation::Reversed {
                        edge_data.data.particle.get_anti_particle(model)
                    } else {
                        edge_data.data.particle.clone()
                    }
                })
                .collect_vec();

            let particle_content = process
                .final_pdgs
                .iter()
                .map(|pdg| model.get_particle_from_pdg(*pdg as isize));

            for particle in particle_content {
                if let Some(index) = cut_content.iter().position(|p| p == &particle) {
                    cut_content.remove(index);
                } else {
                    return Ok(false);
                }
            }

            if cut_content.len() > process.n_unresolved {
                return Ok(false);
            }

            if !cut_content
                .iter()
                .all(|particle| process.unresolved_cut_content.contains(particle))
            {
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
                    return Ok(false);
                }
            }

            if amplitude_couplings.is_some() {
                todo!("waiting for update")
            }

            Ok(true)
        } else {
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
