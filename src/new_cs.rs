use std::path::PathBuf;

use ahash::HashMap;
use bitvec::vec::BitVec;
use color_eyre::Result;

use itertools::Itertools;
use linnet::half_edge::{
    builder::HedgeGraphBuilder,
    involution::Flow,
    subgraph::{HedgeNode, OrientedCut},
    NodeIndex,
};
use serde::{Deserialize, Serialize};
use symbolica::{atom::Atom, parse};
use typed_index_collections::TiVec;

use crate::{
    cff::cut_expression::{OrientationData, OrientationID},
    feyngen::{
        diagram_generator::FeynGen, FeynGenFilter, FeynGenFilters, FeynGenOptions, GenerationType,
        NumeratorAwareGraphGroupingOption, SelfEnergyFilterOptions, SnailFilterOptions,
        TadpolesFilterOptions,
    },
    graph::{BareGraph, BareVertex, EdgeType},
    model::Model,
    new_graph::{Edge, Graph, Vertex},
    numerator::{GlobalPrefactor, NumeratorState, PythonState},
    ProcessSettings,
};

use derive_more::{From, Into};

#[derive(Debug, Clone)]
pub struct ProcessDefinition {
    pub initial_pdgs: Vec<i64>, // Do we want a pub type Pdg = i64;?
    pub final_pdgs: Vec<i64>,
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
                    .map(Graph::forget_type)
            })
            .collect::<Result<Vec<Graph>, _>>()?;

        todo!();
    }

    pub fn from_bare_graph_list(
        bare_graphs: Vec<BareGraph>,
        generation_type: GenerationType,
    ) -> Self {
        todo!();
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
}

#[derive(Clone)]
pub struct Amplitude<S: NumeratorState = PythonState> {
    graphs: Vec<AmplitudeGraph<S>>,
}

#[derive(Clone)]
pub struct AmplitudeGraph<S: NumeratorState = PythonState> {
    graph: Graph<S>,
    derived_data: AmplitudeDerivedData,
}

impl<S: NumeratorState> AmplitudeGraph<S> {
    fn new(graph: Graph<S>) -> Self {
        AmplitudeGraph {
            graph,
            derived_data: AmplitudeDerivedData {},
        }
    }
}

#[derive(Clone)]
pub struct AmplitudeDerivedData {}

impl<S: NumeratorState> Amplitude<S> {
    pub fn new() -> Self {
        Self { graphs: vec![] }
    }

    pub fn add_graph(&mut self, graph: Graph<S>) -> Result<()> {
        self.graphs.push(AmplitudeGraph::new(graph));
        /// TODO: validate that the graph is compatible
        Ok(())
    }
}

#[derive(Clone)]
pub struct CrossSection<S: NumeratorState = PythonState> {
    supergraphs: Vec<CrossSectionGraph<S>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, From, Into, Hash, PartialEq, Copy, Eq)]
pub struct CutId(usize);

#[derive(Clone)]
pub struct CrossSectionGraph<S: NumeratorState = PythonState> {
    graph: Graph<S>,
    cuts: TiVec<CutId, CrossSectionCut>,
    derived_data: CrossSectionDerivedData,
}

impl CrossSectionGraph {
    fn new_from_graph(graph: Graph) -> Self {
        Self {
            graph,
            cuts: TiVec::new(),
            derived_data: CrossSectionDerivedData::new_empty(),
        }
    }
}

#[derive(Clone)]
pub struct CrossSectionDerivedData {
    pub orientations: TiVec<OrientationID, OrientationData>,
}

impl CrossSectionDerivedData {
    fn new_empty() -> Self {
        Self {
            orientations: TiVec::new(),
        }
    }
}

#[derive(Clone)]
pub struct CrossSectionCut {
    pub cut: OrientedCut,
    pub left: BitVec,
    pub right: BitVec,
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
